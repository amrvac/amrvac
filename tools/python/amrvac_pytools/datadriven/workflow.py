"""Notebook-facing orchestration for the single-frame DataConstrain workflow.

This module intentionally contains orchestration only.  Magnetic-data
algorithms live in :mod:`pipeline`, AMRVAC case generation lives in
:mod:`cases`, and relaxation diagnostics live in :mod:`diagnostics`.
"""

from __future__ import print_function

import hashlib
import json
from importlib.util import find_spec
from pathlib import Path
from shlex import quote


DEFAULT_BLOCK_SIZES = (12, 14, 16, 18, 20)
CORE_PACKAGES = ("numpy", "matplotlib", "astropy")
RAW_HMI_PACKAGES = ("scipy", "sunpy")
_DEFAULT_NLFFF_METHOD = object()
NLFFF_METHODS = ("legacy_mfr", "optimization", "grad_rubin")


def check_data_driven_dependencies():
    """Validate core notebook dependencies and report optional raw-HMI extras."""

    missing_core = [name for name in CORE_PACKAGES if find_spec(name) is None]
    if missing_core:
        raise ModuleNotFoundError(
            "Missing required Python package(s): {}. Install them in this Jupyter "
            "kernel environment, then restart the kernel.".format(", ".join(missing_core))
        )
    missing_raw_hmi = [name for name in RAW_HMI_PACKAGES if find_spec(name) is None]
    if missing_raw_hmi:
        print(
            "NOTE: SHARP/CEA input is supported, but raw HMI input additionally "
            "requires: {}".format(", ".join(missing_raw_hmi))
        )
    return {"core": list(CORE_PACKAGES), "raw_hmi": list(RAW_HMI_PACKAGES)}


class StageGridConfig(object):
    """Resolution choices for one AMRVAC stage group (Python 3.6 compatible)."""

    def __init__(
        self,
        boundary_reduction_level=1,
        amrvac_refinement_level=1,
        block_sizes=DEFAULT_BLOCK_SIZES,
    ):
        self.boundary_reduction_level = int(boundary_reduction_level)
        self.amrvac_refinement_level = int(amrvac_refinement_level)
        self.block_sizes = tuple(int(value) for value in block_sizes)
        if self.boundary_reduction_level < 1:
            raise ValueError("boundary_reduction_level must be >= 1")
        if self.amrvac_refinement_level < 1:
            raise ValueError("amrvac_refinement_level must be >= 1")
        if not self.block_sizes or any(value < 1 for value in self.block_sizes):
            raise ValueError("block_sizes must contain positive integers")

    def as_dict(self):
        return {
            "boundary_reduction_level": self.boundary_reduction_level,
            "amrvac_refinement_level": self.amrvac_refinement_level,
            "block_sizes": list(self.block_sizes),
        }


class DataConstrainWorkflow(object):
    """Coordinate the single-frame Potential -> MFR -> DataConstrained chain."""

    def __init__(
        self,
        amrvac_root,
        project_dir,
        input_dir,
        relaxation_grid=None,
        evolution_grid=None,
        manifest_workflow="data_constrain",
    ):
        self.amrvac_root = Path(amrvac_root).expanduser().resolve()
        if project_dir is None:
            project_dir = self.amrvac_root / "DrivenFieldProject"
        self.project_dir = Path(project_dir).expanduser().resolve()
        self.input_dir = Path(input_dir).expanduser().resolve()
        self.relaxation_grid = relaxation_grid or StageGridConfig(2, 1)
        self.evolution_grid = evolution_grid or StageGridConfig(1, 2)
        self._manifest_workflow = str(manifest_workflow)
        if not isinstance(self.relaxation_grid, StageGridConfig):
            raise TypeError("relaxation_grid must be a StageGridConfig")
        if not isinstance(self.evolution_grid, StageGridConfig):
            raise TypeError("evolution_grid must be a StageGridConfig")
        if self.relaxation_grid.block_sizes != self.evolution_grid.block_sizes:
            raise ValueError("relaxation and evolution stages must use the same block_sizes")
        self.region = None
        self._input_state = None
        self._preview_state = None
        self._region_plot_state = None
        self._initial_field_state = None
        self._data_constrained_state = None

    @property
    def paths(self):
        root = self.project_dir
        return {
            "project": root,
            "manifest": root / "project.json",
            "prepared_magnetograms": root / "PreparedMagnetograms",
            "boundaries": root / "MagneticBoundary",
            "potential": root / "PotentialField",
            "mfr": root / "MagnetofrictionalRelaxation",
            "optimization": root / "Optimization",
            "grad_rubin": root / "GradRubin",
            "alpha_products": root / "AlphaProducts",
            "data_constrained": root / "DataConstrained",
        }

    def inspect_input(self, refresh=False):
        from .fits_io import inspect_magnetic_input, summarize_magnetic_input

        if self._input_state is not None and not refresh:
            return self._input_state
        info = inspect_magnetic_input(self.input_dir)
        summary = summarize_magnetic_input(info)
        if info["kind"] == "unknown":
            raise ValueError(info["message"])
        self._input_state = {
            "info": info,
            "summary": summary,
            "paths": self._serializable_paths(),
        }
        return self._input_state

    def preview_input(self, bmax=1000.0, show=True):
        """Load and plot the first Br-like frame without exposing plot glue."""

        from .fits_io import load_magnetic_preview

        inspected = self.inspect_input()
        preview = load_magnetic_preview(inspected["info"])
        plot = plot_magnetic_input_preview(preview, bmax=bmax, show=show)
        inspected.update({"preview": preview, "plot": plot})
        self._preview_state = inspected
        return inspected

    def plan_region(
        self,
        preview_shape=None,
        cea_patch=None,
        window=None,
        auto_trim=True,
        region_mode="auto",
        raw_hmi_cea_patch=None,
        sharp_pixel_window=None,
    ):
        """Resolve the input-appropriate region and plan both stage grids.

        ``cea_patch`` and ``window`` remain backwards-compatible aliases for
        ``raw_hmi_cea_patch`` and ``sharp_pixel_window`` respectively.
        """

        if preview_shape is None:
            preview_shape = self._require_preview()["preview"]["data"].shape
        info = self.inspect_input()["info"]
        self.region = plan_workflow_region(
            input_info=info,
            preview_shape=preview_shape,
            cea_patch=cea_patch,
            window=window,
            region_mode=region_mode,
            raw_hmi_cea_patch=raw_hmi_cea_patch,
            sharp_pixel_window=sharp_pixel_window,
            relaxation_grid=self.relaxation_grid,
            evolution_grid=self.evolution_grid,
            auto_trim=auto_trim,
        )
        self._update_manifest()
        return self.region

    def plot_region(self, preview_data=None, bmax=1000.0, show=True):
        if self.region is None:
            raise RuntimeError("plan_region must be called before plot_region")
        if preview_data is None:
            preview_data = self._require_preview()["preview"]["data"]
        info = self.inspect_input()["info"]
        self._region_plot_state = plot_data_constrain_region(
            info,
            preview_data,
            self.region,
            bmax=bmax,
            show=show,
        )
        return self._region_plot_state

    def prepare_initial_field(
        self,
        snapshot_index=0,
        cea_remap_options=None,
        geometry=None,
        preprocess=False,
        preprocessing_mode=None,
        fail_on_nonconvergence=False,
        nghost=2,
        quicklook=True,
        vmax=500.0,
        stage_potential=True,
        stage_mfr=None,
        nlfff_method=_DEFAULT_NLFFF_METHOD,
        potential_options=None,
        mfr_options=None,
        optimization_options=None,
        grad_rubin_options=None,
        gr_alpha_cleaning=False,
        gr_alpha_preset="recommended",
        gr_alpha_options=None,
        unit_length_cm=1.0e9,
        unit_magneticfield_g=100.0,
        potential_restart_file=None,
    ):
        """Prepare one frame and stage PotentialField plus one NLFFF method.

        Raw HMI sequences are deliberately remapped one selected frame at a
        time.  This method never calls the full-sequence conversion API.
        """

        if self.region is None:
            raise RuntimeError("plan_region must be called before prepare_initial_field")
        from .cases import (
            stage_grad_rubin_nlfff_case,
            stage_magnetofrictional_relaxation_case,
            stage_optimization_nlfff_case,
            stage_potential_field_case,
        )
        from .fits_io import inspect_magnetic_input
        from .pipeline import make_cea_patch, prepare_boundary_frame

        if nlfff_method is _DEFAULT_NLFFF_METHOD:
            selected_method = "legacy_mfr" if stage_mfr is not False else "potential"
        else:
            selected_method = _normalize_nlfff_method(nlfff_method)
        if selected_method != "potential" and selected_method not in NLFFF_METHODS:
            raise ValueError("NLFFF_METHOD must be one of {}".format(", ".join(NLFFF_METHODS)))

        snapshot_index = int(snapshot_index)
        if snapshot_index < 0:
            raise ValueError("snapshot_index must be >= 0")
        info = self.inspect_input()["info"]
        boundary_dir, boundary_snapshot_index, source_frame = self._prepare_selected_frame(
            info,
            snapshot_index,
            cea_remap_options or {},
            make_cea_patch,
            inspect_magnetic_input,
        )

        paths = self.paths
        relaxation_meta = prepare_boundary_frame(
            boundary_dir,
            paths["boundaries"] / "Relaxation" / "SingleFrame",
            window=self.region["window"],
            level=self.relaxation_grid.boundary_reduction_level,
            geometry=geometry,
            snapshot_index=boundary_snapshot_index,
            preprocess=preprocess,
            preprocessing_mode=preprocessing_mode,
            fail_on_nonconvergence=fail_on_nonconvergence,
            nghost=nghost,
            quicklook=quicklook,
            vmax=vmax,
        )
        evolution_meta = prepare_boundary_frame(
            boundary_dir,
            paths["boundaries"] / "Evolution" / "SingleFrame",
            window=self.region["window"],
            level=self.evolution_grid.boundary_reduction_level,
            geometry=geometry,
            snapshot_index=boundary_snapshot_index,
            preprocess=preprocess,
            preprocessing_mode=preprocessing_mode,
            fail_on_nonconvergence=fail_on_nonconvergence,
            nghost=nghost,
            quicklook=quicklook,
            vmax=vmax,
        )
        _validate_product_plan("relaxation", relaxation_meta, self.region["plans"]["relaxation"])
        _validate_product_plan("evolution", evolution_meta, self.region["plans"]["evolution"])

        potential_summary = None
        mfr_summary = None
        nlfff_summary = None
        alpha_summary = None
        relaxation_boundary = _first_boundary(relaxation_meta)
        if gr_alpha_cleaning:
            alpha_summary = _write_grad_rubin_external_alpha_from_boundary(
                relaxation_meta,
                relaxation_boundary,
                paths["alpha_products"],
                preset=gr_alpha_preset,
                alpha_options=gr_alpha_options,
                unit_length_cm=unit_length_cm,
                unit_magneticfield_g=unit_magneticfield_g,
                used_by_method=(selected_method == "grad_rubin"),
                nlfff_method=selected_method,
            )
        if stage_potential:
            options = dict(potential_options or {})
            options.setdefault("portable_paths", True)
            potential_summary = stage_potential_field_case(
                relaxation_meta,
                paths["potential"],
                amrvac_root=self.amrvac_root,
                boundary_filename=relaxation_boundary,
                refine_max_level=self.relaxation_grid.amrvac_refinement_level,
                block_nx1=self.region["plans"]["relaxation"]["x"]["block_size"],
                block_nx2=self.region["plans"]["relaxation"]["y"]["block_size"],
                block_nx3=self.region["plans"]["relaxation"]["y"]["block_size"],
                **options
            )

        if potential_restart_file is None:
            potential_restart_file = paths["potential"] / "output" / "data_driven_potential0000.dat"
        if selected_method == "legacy_mfr":
            options = dict(mfr_options or {})
            options.setdefault("portable_paths", True)
            mfr_summary = stage_magnetofrictional_relaxation_case(
                relaxation_meta,
                paths["mfr"],
                potential_restart_file=potential_restart_file,
                amrvac_root=self.amrvac_root,
                boundary_filename=relaxation_boundary,
                refine_max_level=self.relaxation_grid.amrvac_refinement_level,
                block_nx1=self.region["plans"]["relaxation"]["x"]["block_size"],
                block_nx2=self.region["plans"]["relaxation"]["y"]["block_size"],
                block_nx3=self.region["plans"]["relaxation"]["y"]["block_size"],
                **options
            )
            nlfff_summary = mfr_summary
        elif selected_method == "optimization":
            options = dict(optimization_options or {})
            options.setdefault("portable_paths", True)
            nlfff_summary = stage_optimization_nlfff_case(
                relaxation_meta,
                paths["optimization"],
                potential_restart_file=potential_restart_file,
                amrvac_root=self.amrvac_root,
                boundary_filename=relaxation_boundary,
                refine_max_level=self.relaxation_grid.amrvac_refinement_level,
                block_nx1=self.region["plans"]["relaxation"]["x"]["block_size"],
                block_nx2=self.region["plans"]["relaxation"]["y"]["block_size"],
                block_nx3=self.region["plans"]["relaxation"]["y"]["block_size"],
                **options
            )
        elif selected_method == "grad_rubin":
            options = dict(grad_rubin_options or {})
            options.setdefault("portable_paths", True)
            if gr_alpha_cleaning:
                options["gr_alpha_source"] = "external"
                options["external_alpha_filename"] = alpha_summary["external_alpha"]
            else:
                options.setdefault("gr_alpha_source", "vector_magnetogram")
            nlfff_summary = stage_grad_rubin_nlfff_case(
                relaxation_meta,
                paths["grad_rubin"],
                potential_restart_file=potential_restart_file,
                amrvac_root=self.amrvac_root,
                boundary_filename=relaxation_boundary,
                refine_max_level=self.relaxation_grid.amrvac_refinement_level,
                block_nx1=self.region["plans"]["relaxation"]["x"]["block_size"],
                block_nx2=self.region["plans"]["relaxation"]["y"]["block_size"],
                block_nx3=self.region["plans"]["relaxation"]["y"]["block_size"],
                **options
            )

        result = {
            "snapshot_index": snapshot_index,
            "source_frame": source_frame,
            "boundary_input_dir": str(boundary_dir),
            "relaxation_boundary": relaxation_meta,
            "evolution_boundary": evolution_meta,
            "potential_case": potential_summary,
            "mfr_case": mfr_summary,
            "nlfff_method": selected_method,
            "nlfff_case": nlfff_summary,
            "alpha_product": alpha_summary,
            "potential_restart_file": str(Path(potential_restart_file)),
            "summary": _initial_field_summary(
                self.project_dir,
                relaxation_meta,
                evolution_meta,
                potential_summary,
                nlfff_summary,
                selected_method,
                alpha_summary,
            ),
        }
        self._update_manifest({
            "reference_frame": {
                "snapshot_index": snapshot_index,
                "source_frame": source_frame,
                "boundary_input_dir": str(boundary_dir),
            },
            "initial_field": result["summary"],
        })
        self._initial_field_state = result
        return result

    def analyze_and_stage_data_constrained(
        self,
        evolution_boundary=None,
        restart_file=None,
        diagnostics_csv=None,
        selected_restart_number=None,
        mf_ditsave=None,
        restart_prefix="data_driven_mfr",
        show_plot=True,
        plot_lorentz_force=False,
        data_constrained_options=None,
    ):
        """Select an MFR checkpoint, plot diagnostics, and stage DataConstrained."""

        from .cases import stage_data_constrained_case
        from .diagnostics import (
            find_relaxation_restart_snapshots,
            normalize_nlfff_metrics,
            plot_relaxation_diagnostics,
            read_relaxation_diagnostics,
            select_relaxation_restart,
        )

        if self.region is None:
            raise RuntimeError("plan_region must be called before staging DataConstrained")
        if evolution_boundary is None:
            if self._initial_field_state is None:
                raise RuntimeError(
                    "prepare_initial_field must be called before staging DataConstrained"
                )
            evolution_boundary = self._initial_field_state["evolution_boundary"]
        paths = self.paths
        if mf_ditsave is None:
            mfr_case = None
            if self._initial_field_state is not None:
                mfr_case = self._initial_field_state.get("mfr_case")
            mf_ditsave = (mfr_case or {}).get("mf_ditsave", 20000)
        diagnostics = None
        restart_summary = None
        selection = None
        plot_summary = None
        metrics_summary = None
        initial = self._initial_field_state or {}
        nlfff_case = initial.get("nlfff_case")
        nlfff_method = initial.get("nlfff_method", "legacy_mfr")
        if restart_file is None:
            if nlfff_case is None:
                raise RuntimeError("prepare_initial_field did not stage an NLFFF case")
            metrics_summary = normalize_nlfff_metrics(nlfff_case)
            if nlfff_method == "legacy_mfr":
                output_dir = paths["mfr"] / "output"
                if diagnostics_csv is None:
                    diagnostics_csv = metrics_summary.get("path")
                diagnostics = read_relaxation_diagnostics(diagnostics_csv)
                restart_summary = find_relaxation_restart_snapshots(
                    output_dir=output_dir,
                    base_filename=restart_prefix,
                    diagnostics=diagnostics,
                    mf_ditsave=mf_ditsave,
                )
                selection = select_relaxation_restart(
                    diagnostics=diagnostics,
                    output_dir=output_dir,
                    base_filename=restart_prefix,
                    mf_ditsave=mf_ditsave,
                    snapshot_number=selected_restart_number,
                )
                restart_file = selection.get("restart_file")
                if restart_file is None:
                    raise RuntimeError("no magnetofrictional restart snapshot was found")
                selected = selection.get("selected_marker") or {}
                markers = [dict(marker) for marker in restart_summary["markers"]]
                for marker in markers:
                    marker["selected"] = marker.get("snapshot_number") == selected.get("snapshot_number")
                plot_columns = ["cw_sin_theta"]
                if plot_lorentz_force:
                    plot_columns.append("lorentz_force")
                plot_summary = plot_relaxation_diagnostics(
                    diagnostics,
                    output_path=output_dir / "relaxation_diagnostics_quicklook.png",
                    columns=plot_columns,
                    x_column="iteration",
                    restart_markers=markers,
                )
                if show_plot:
                    import matplotlib.pyplot as plt

                    plt.show()
            else:
                restart_file = _one_shot_nlfff_restart_file(nlfff_case)
                selection = {
                    "restart_file": restart_file,
                    "reason": "{} one-shot NLFFF output selected".format(nlfff_method),
                    "warnings": list(metrics_summary.get("warnings", [])),
                    "diagnostic_values": {},
                    "metrics": metrics_summary,
                }
        else:
            restart_file = str(Path(restart_file).expanduser().resolve())
            selection = {
                "restart_file": restart_file,
                "reason": "explicit restart file selected",
                "warnings": [],
                "diagnostic_values": {},
            }

        plan = self.region["plans"]["evolution"]
        options = dict(data_constrained_options or {})
        options.setdefault("portable_paths", True)
        case_summary = stage_data_constrained_case(
            evolution_boundary,
            paths["data_constrained"],
            restart_file=restart_file,
            amrvac_root=self.amrvac_root,
            boundary_filename=_first_boundary(evolution_boundary),
            refine_max_level=self.evolution_grid.amrvac_refinement_level,
            block_nx1=plan["x"]["block_size"],
            block_nx2=plan["y"]["block_size"],
            block_nx3=plan["y"]["block_size"],
            **options
        )
        warnings = []
        if restart_summary is not None:
            warnings.extend(restart_summary.get("warnings", []))
        warnings.extend(selection.get("warnings", []))
        result = {
            "diagnostics": diagnostics,
            "restart_summary": restart_summary,
            "selection": selection,
            "plot": plot_summary,
            "metrics": metrics_summary,
            "case": case_summary,
            "warnings": warnings,
        }
        self._update_manifest({
            "data_constrained": {
                "case_dir": case_summary.get("case_dir"),
                "restart_file": selection.get("restart_file"),
                "mhd_model": case_summary.get("mhd_model"),
                "atmosphere_model": case_summary.get("atmosphere_model"),
            }
        })
        self._data_constrained_state = result
        return result

    def case_commands(self, case_summary, nproc=4):
        return format_case_commands(case_summary, self.amrvac_root, nproc=nproc)

    def input_report(self):
        """Return a concise notebook report for the inspected magnetic input."""

        state = self._require_preview()
        summary = state["summary"]
        lines = [
            "Project directory: {}".format(self.project_dir),
            "Input kind: {}".format(summary["kind"]),
            "Frames found: {}".format(summary.get("frame_count", 1)),
            "Preview shape (ny, nx): {}".format(state["preview"]["data"].shape),
        ]
        return _report_with_warnings(lines, summary.get("warnings", []))

    def region_report(self):
        """Return grid sizes for the selected master region and both stage groups."""

        if self.region is None:
            raise RuntimeError("plan_region must be called before region_report")
        master = self.region["master"]
        lines = [
            "Master pixels: ({}, {}) -> ({}, {}); trim=({}, {})".format(
                master["original_nx"],
                master["original_ny"],
                master["adjusted_nx"],
                master["adjusted_ny"],
                master["trim_x"],
                master["trim_y"],
            )
        ]
        for label in ("relaxation", "evolution"):
            plan = self.region["plans"][label]
            lines.append(
                "{}: boundary=({}, {}), base=({}, {}), blocks=({}, {}), AMR level={}".format(
                    label.capitalize(),
                    plan["x"]["boundary_physical"],
                    plan["y"]["boundary_physical"],
                    plan["x"]["base"],
                    plan["y"]["base"],
                    plan["x"]["block_size"],
                    plan["y"]["block_size"],
                    plan["amrvac_refinement_level"],
                )
            )
        warnings = [] if self._region_plot_state is None else self._region_plot_state.get("warnings", [])
        return _report_with_warnings(lines, warnings)

    def initial_field_report(self):
        """Return a concise report for prepared boundaries and staged cases."""

        state = self._require_initial_field()
        summary = state["summary"]
        return "\n".join([
            "Reference frame: {}".format(state["source_frame"]),
            "Relaxation boundary: {} physical pixels".format(
                summary["relaxation_boundary"]["physical_grid"]
            ),
            "Evolution boundary: {} physical pixels".format(
                summary["evolution_boundary"]["physical_grid"]
            ),
            "PotentialField case: {}".format(summary["potential_case_dir"]),
            "Selected NLFFF method: {}".format(summary["nlfff_method"]),
            "Selected NLFFF case: {}".format(summary["nlfff_case_dir"]),
            "Unified NLFFF metrics: {}".format(summary["nlfff_metrics_csv"]),
            "GR external alpha: {}".format(summary.get("external_alpha") or "not used"),
        ])

    def initial_field_commands(self, nproc=4):
        """Return terminal commands for PotentialField followed by selected NLFFF."""

        state = self._require_initial_field()
        sections = []
        for title, case in (
            ("1. Potential-field extrapolation", state["potential_case"]),
            ("2. Selected NLFFF method ({})".format(state["nlfff_method"]), state["nlfff_case"]),
        ):
            if case is not None:
                sections.append(
                    "{}\n{}".format(title, "\n".join(self.case_commands(case, nproc=nproc)))
                )
        return "\n\n".join(sections)

    def normalize_initial_nlfff_metrics(self, overwrite=False):
        """Ensure/report the selected method's common ``_nlfff_metrics.csv``."""

        state = self._require_initial_field()
        case = state.get("nlfff_case")
        if case is None:
            return {
                "path": None,
                "status": "no_nlfff_case",
                "method": state.get("nlfff_method"),
            }
        from .diagnostics import normalize_nlfff_metrics

        result = normalize_nlfff_metrics(case, overwrite=overwrite)
        state["nlfff_metrics"] = result
        self._update_manifest({"initial_nlfff_metrics": result})
        return result

    def initial_nlfff_metrics_report(self):
        state = self._require_initial_field()
        metrics = state.get("nlfff_metrics")
        if metrics is None:
            metrics = self.normalize_initial_nlfff_metrics()
        lines = [
            "Selected NLFFF method: {}".format(state.get("nlfff_method")),
            "Unified metrics CSV: {}".format(metrics.get("path")),
            "Metrics status: {}".format(metrics.get("status")),
        ]
        if metrics.get("source"):
            lines.append("Metrics source: {}".format(metrics["source"]))
        if metrics.get("adapter_audit"):
            lines.append("Adapter audit: {}".format(metrics["adapter_audit"]))
        return _report_with_warnings(lines, metrics.get("warnings", []))

    def plot_initial_nlfff_metrics(self, show=False, output_path=None):
        """Create the common metrics quicklook for the selected NLFFF method.

        The public notebook path is method-neutral.  Only legacy MFR is
        allowed to add explicitly discovered checkpoint markers; Optimization
        and Grad--Rubin retain their one-shot ``0000.dat`` semantics and never
        search MFR restart files here.
        """

        state = self._require_initial_field()
        metrics = state.get("nlfff_metrics")
        if metrics is None:
            metrics = self.normalize_initial_nlfff_metrics()
        if not metrics.get("path"):
            raise RuntimeError(
                "selected NLFFF method has no common metrics CSV: {}".format(
                    metrics.get("status")
                )
            )

        from .diagnostics import (
            find_relaxation_restart_snapshots,
            plot_nlfff_metrics,
            read_nlfff_metrics,
        )

        method = str(state.get("nlfff_method", "legacy_mfr"))
        markers = []
        marker_warnings = []
        marker_source = "none"
        if method == "legacy_mfr":
            mfr_case = state.get("mfr_case") or state.get("nlfff_case") or {}
            case_dir = mfr_case.get("case_dir")
            if case_dir:
                diagnostic_data = read_nlfff_metrics(metrics["path"])
                base_filename = Path(
                    str(mfr_case.get("base_filename", "output/data_driven_mfr"))
                ).name
                restart_summary = find_relaxation_restart_snapshots(
                    output_dir=Path(case_dir).expanduser().resolve() / "output",
                    base_filename=base_filename,
                    diagnostics=diagnostic_data,
                    mf_ditsave=mfr_case.get("mf_ditsave", 20000),
                )
                markers = [dict(marker) for marker in restart_summary.get("markers", [])]
                marker_warnings.extend(restart_summary.get("warnings", []))
                marker_source = "legacy_mfr_explicit_checkpoint_search"

        plot = plot_nlfff_metrics(
            metrics,
            output_path=output_path,
            method=method,
            restart_markers=markers,
            show=show,
        )
        plot["restart_search"] = marker_source
        plot["warnings"] = list(dict.fromkeys(marker_warnings + plot.get("warnings", [])))
        state["nlfff_metrics_plot"] = plot
        manifest_plot = {
            key: plot.get(key)
            for key in (
                "output_path",
                "metrics_csv",
                "method",
                "plotted_columns",
                "skipped_columns",
                "restart_search",
                "warnings",
            )
        }
        self._update_manifest({"initial_nlfff_metrics_plot": manifest_plot})
        return plot

    def data_constrained_report(self):
        """Return the selected restart diagnostics and staged-case summary."""

        state = self._require_data_constrained()
        selection = state["selection"]
        selected_marker = selection.get("selected_marker") or {}
        lines = ["Selected restart: {}".format(selection["restart_file"])]
        if selected_marker.get("iteration") is not None:
            lines.append("Selected iteration: {}".format(selected_marker["iteration"]))
        if selection.get("reason"):
            lines.append("Reason: {}".format(selection["reason"]))
        labels = {"cw_sin_theta": "sin(theta)"}
        for key, value in selection.get("diagnostic_values", {}).items():
            lines.append("  {}: {:.6e}".format(labels.get(key, key), value))
        case = state["case"]
        lines.extend([
            "DataConstrained case: {}".format(case["case_dir"]),
            "Model: {} / {}".format(case["mhd_model"], case["atmosphere_model"]),
        ])
        return _report_with_warnings(lines, state.get("warnings", []), prepend=True)

    def data_constrained_commands(self, nproc=4):
        """Return terminal commands for the staged DataConstrained case."""

        state = self._require_data_constrained()
        return "\n".join(self.case_commands(state["case"], nproc=nproc))

    def _require_preview(self):
        if self._preview_state is None:
            raise RuntimeError("preview_input must be called first")
        return self._preview_state

    def _require_initial_field(self):
        if self._initial_field_state is None:
            raise RuntimeError("prepare_initial_field must be called first")
        return self._initial_field_state

    def _require_data_constrained(self):
        if self._data_constrained_state is None:
            raise RuntimeError("analyze_and_stage_data_constrained must be called first")
        return self._data_constrained_state

    def _prepare_selected_frame(self, info, snapshot_index, remap_options, remapper, inspector):
        kind = info["kind"]
        if kind in ("raw_hmi_vector", "raw_hmi_vector_sequence"):
            frames = info.get("raw_hmi_sequence", {}).get("complete", [])
            if snapshot_index >= len(frames):
                raise IndexError("snapshot_index {} exceeds {} raw HMI frames".format(snapshot_index, len(frames)))
            frame = frames[snapshot_index]
            output_dir = self.paths["prepared_magnetograms"] / "ReferenceFrame"
            options = dict(remap_options)
            options.setdefault("components_only", True)
            options.setdefault("quicklook", True)
            remapper(
                field=frame["field"],
                inclination=frame["inclination"],
                azimuth=frame["azimuth"],
                disambig=frame["disambig"],
                output_dir=output_dir,
                **dict(self.region["cea_patch"], **options)
            )
            converted = inspector(output_dir)
            if converted["kind"] not in ("cea_vector", "cea_vector_sequence"):
                raise RuntimeError("selected raw HMI frame did not produce Br/Bt/Bp files")
            return output_dir, 0, frame.get("key")
        if kind in ("cea_vector", "cea_vector_sequence"):
            count = len(info.get("components", {}).get("br", []))
            if snapshot_index >= count:
                raise IndexError("snapshot_index {} exceeds {} CEA frames".format(snapshot_index, count))
            frame_path = info["components"]["br"][snapshot_index]
            return Path(info["directory"]), snapshot_index, Path(frame_path).name
        raise ValueError("DataConstrain requires all three Br/Bt/Bp components")

    def _serializable_paths(self):
        return {key: str(value) for key, value in self.paths.items()}

    def _update_manifest(self, updates=None):
        """Persist portable workflow state for later notebooks and migration."""

        from .manifest import update_project_manifest

        shared = {
            "input_dir": str(self.input_dir),
            "project_dir": str(self.project_dir),
            "paths": self._serializable_paths(),
            "relaxation_grid": self.relaxation_grid.as_dict(),
            "evolution_grid": self.evolution_grid.as_dict(),
        }
        if self.region is not None:
            shared["region"] = self.region
        update_project_manifest(
            self.paths["manifest"],
            shared=shared,
            workflow=self._manifest_workflow,
            workflow_updates=updates or {},
        )
        return self.paths["manifest"]


def create_data_constrain_workflow(*args, **kwargs):
    """Lazy-constructor-friendly public entry point used by notebooks."""

    return DataConstrainWorkflow(*args, **kwargs)


def plan_data_constrain_region(
    input_info,
    preview_shape,
    cea_patch=None,
    window=None,
    relaxation_grid=None,
    evolution_grid=None,
    auto_trim=True,
):
    """Return an adjusted master region and two mutually compatible grid plans."""

    relaxation_grid = relaxation_grid or StageGridConfig(2, 1)
    evolution_grid = evolution_grid or StageGridConfig(1, 2)
    if relaxation_grid.block_sizes != evolution_grid.block_sizes:
        raise ValueError("relaxation and evolution stages must use the same block_sizes")
    kind = input_info["kind"]
    is_raw = kind in ("raw_hmi_vector", "raw_hmi_vector_sequence")
    if is_raw:
        if cea_patch is None:
            raise ValueError("raw HMI input requires cea_patch")
        from .cea import make_cea_patch_grid

        patch = dict(cea_patch)
        grid = make_cea_patch_grid(
            patch["center_lon"],
            patch["center_lat"],
            patch["width_degree"],
            patch["height_degree"],
            patch.get("resolution_degree", 0.03),
        )
        ny, nx = grid.shape
        selected_window = None
    else:
        patch = None
        selected_window = _checked_window(window, preview_shape)
        nx, ny = selected_window["nx"], selected_window["ny"]

    adjusted_nx = _combined_adjusted_dimension(nx, relaxation_grid, evolution_grid)
    adjusted_ny = _combined_adjusted_dimension(ny, relaxation_grid, evolution_grid)
    trim_x, trim_y = nx - adjusted_nx, ny - adjusted_ny
    if (trim_x or trim_y) and not auto_trim:
        raise ValueError(
            "selected region is not compatible with both AMRVAC grids; suggested size is {} x {}".format(
                adjusted_nx, adjusted_ny
            )
        )
    if is_raw:
        resolution = float(patch.get("resolution_degree", 0.03))
        patch["width_degree"] = adjusted_nx * resolution
        patch["height_degree"] = adjusted_ny * resolution
    else:
        selected_window["nx"] = adjusted_nx
        selected_window["ny"] = adjusted_ny

    plans = {
        "relaxation": _exact_stage_plan(adjusted_nx, adjusted_ny, relaxation_grid),
        "evolution": _exact_stage_plan(adjusted_nx, adjusted_ny, evolution_grid),
    }
    return {
        "input_kind": kind,
        "selection_mode": "fixed_cea_patch" if is_raw else "cea_pixel_window",
        "cea_patch": patch,
        "window": selected_window,
        "master": {
            "original_nx": int(nx),
            "original_ny": int(ny),
            "adjusted_nx": int(adjusted_nx),
            "adjusted_ny": int(adjusted_ny),
            "trim_x": int(trim_x),
            "trim_y": int(trim_y),
        },
        "plans": plans,
    }


def plan_workflow_region(
    input_info,
    preview_shape,
    cea_patch=None,
    window=None,
    relaxation_grid=None,
    evolution_grid=None,
    auto_trim=True,
    region_mode="auto",
    raw_hmi_cea_patch=None,
    sharp_pixel_window=None,
):
    """Resolve auto/manual selection and call the shared low-level planner."""

    mode = str(region_mode).strip().lower()
    if mode not in ("auto", "fixed_cea", "pixel_window"):
        raise ValueError("region_mode must be 'auto', 'fixed_cea', or 'pixel_window'")
    kind = input_info["kind"]
    is_raw = kind in ("raw_hmi_vector", "raw_hmi_vector_sequence")
    is_cea = kind in ("cea_vector", "cea_vector_sequence")
    if not (is_raw or is_cea):
        raise ValueError("region planning requires raw HMI or SHARP/CEA vector input")

    if raw_hmi_cea_patch is not None and cea_patch is not None and raw_hmi_cea_patch != cea_patch:
        raise ValueError("raw_hmi_cea_patch and legacy cea_patch disagree")
    if sharp_pixel_window is not None and window is not None and sharp_pixel_window != window:
        raise ValueError("sharp_pixel_window and legacy window disagree")
    patch = raw_hmi_cea_patch if raw_hmi_cea_patch is not None else cea_patch
    pixel_window = sharp_pixel_window if sharp_pixel_window is not None else window

    resolved_mode = ("fixed_cea" if is_raw else "pixel_window") if mode == "auto" else mode
    if is_raw and resolved_mode == "pixel_window":
        raise ValueError("raw HMI sequences only support region_mode='fixed_cea'")

    planner_info = input_info
    if resolved_mode == "fixed_cea":
        if patch is None:
            raise ValueError("region_mode='fixed_cea' requires RAW_HMI_CEA_PATCH")
        # The low-level planner identifies fixed-CEA selection by raw-HMI kind;
        # use a proxy for an explicit SHARP/CEA reprojection request.
        planner_info = dict(input_info)
        planner_info["kind"] = "raw_hmi_vector_sequence"
        selected_patch = patch
        selected_window = None
    else:
        selected_patch = None
        selected_window = pixel_window

    result = plan_data_constrain_region(
        input_info=planner_info,
        preview_shape=preview_shape,
        cea_patch=selected_patch,
        window=selected_window,
        relaxation_grid=relaxation_grid,
        evolution_grid=evolution_grid,
        auto_trim=auto_trim,
    )
    result["input_kind"] = kind
    result["requested_region_mode"] = mode
    result["resolved_region_mode"] = resolved_mode
    return result


def plot_data_constrain_region(input_info, preview_data, region, bmax=1000.0, show=True):
    """Draw the full magnetogram and selected DataConstrain master region."""

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    data = np.asarray(preview_data)
    finite = data[np.isfinite(data)]
    if bmax is None:
        bmax = float(np.nanpercentile(np.abs(finite), 99.5)) if finite.size else 1.0
    bmax = max(float(bmax), 1.0)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), dpi=120, constrained_layout=True)
    axes[0].imshow(data, origin="lower", cmap="gray", vmin=-bmax, vmax=bmax)
    is_raw_hmi = input_info["kind"] in ("raw_hmi_vector", "raw_hmi_vector_sequence")
    field_label = r"$B_{\mathrm{LOS}}$" if is_raw_hmi else r"$B_r$"
    warnings = []
    if region["selection_mode"] == "fixed_cea_patch":
        from .fits_io import load_cea_patch_preview

        try:
            patch_preview = load_cea_patch_preview(
                input_info,
                region["cea_patch"],
                preview_data=data,
            )
        except Exception as error:
            patch_preview = None
            warnings.append(
                "CEA patch preview unavailable: {}: {}".format(type(error).__name__, error)
            )
        axes[0].set_title("(a) HMI full disk")
        axes[0].set_xlabel(r"$x$ [pixel]")
        axes[0].set_ylabel(r"$y$ [pixel]")
        if patch_preview is not None:
            axes[0].plot(
                patch_preview["outline_x"],
                patch_preview["outline_y"],
                color="red",
                linewidth=1.5,
            )
            selected = patch_preview["data"]
            image = axes[1].imshow(
                selected,
                origin="lower",
                cmap="gray",
                vmin=-bmax,
                vmax=bmax,
            )
            axes[1].set_title("(b) Selected CEA domain")
            axes[1].set_xlabel(r"$x$ [pixel]")
            axes[1].set_ylabel(r"$y$ [pixel]")
        else:
            selected = None
            image = axes[0].images[0]
            axes[1].axis("off")
            axes[1].set_title("(b) CEA patch preview unavailable")
            axes[1].text(
                0.5,
                0.5,
                "The actual CEA remap can still run.\nCheck SunPy/Astropy/WCS metadata if needed.",
                ha="center",
                va="center",
                wrap=True,
                transform=axes[1].transAxes,
            )
    else:
        selected_window = region["window"]
        selected = data[
            selected_window["y0"]:selected_window["y0"] + selected_window["ny"],
            selected_window["x0"]:selected_window["x0"] + selected_window["nx"],
        ]
        axes[0].add_patch(Rectangle(
            (selected_window["x0"] - 0.5, selected_window["y0"] - 0.5),
            selected_window["nx"], selected_window["ny"],
            fill=False, edgecolor="red", linewidth=1.6,
        ))
        axes[0].set_title("(a) Full CEA map")
        axes[0].set_xlabel(r"$x$ [pixel]")
        axes[0].set_ylabel(r"$y$ [pixel]")
        image = axes[1].imshow(selected, origin="lower", cmap="gray", vmin=-bmax, vmax=bmax)
        axes[1].set_title("(b) Selected region")
        axes[1].set_xlabel(r"$x$ [pixel]")
        axes[1].set_ylabel(r"$y$ [pixel]")
    colorbar_parent = axes[1] if selected is not None else axes[0]
    divider = make_axes_locatable(colorbar_parent)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    fig.colorbar(image, cax=cax, label="{} [G]".format(field_label))
    if show:
        plt.show()
    return {
        "figure": fig,
        "axes": axes,
        "selected_shape": tuple(selected.shape) if selected is not None else None,
        "warnings": warnings,
    }


def plot_magnetic_input_preview(preview, bmax=1000.0, show=True):
    """Plot one compact Br-like input preview with symmetric color limits."""

    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    data = np.asarray(preview["data"])
    finite = data[np.isfinite(data)]
    if bmax is None:
        bmax = float(np.nanpercentile(np.abs(finite), 99.5)) if finite.size else 1.0
    bmax = max(float(bmax), 1.0)

    fig, axis = plt.subplots(figsize=(7, 5), dpi=120)
    image = axis.imshow(data, origin="lower", cmap="gray", vmin=-bmax, vmax=bmax)
    axis.set_title(preview["title"])
    axis.set_xlabel(r"$x$ [pixel]")
    axis.set_ylabel(r"$y$ [pixel]")
    divider = make_axes_locatable(axis)
    color_axis = divider.append_axes("right", size="3%", pad=0.05)
    quantity_label = preview.get("quantity_label", preview["quantity"])
    fig.colorbar(image, cax=color_axis, label="{} [G]".format(quantity_label))
    fig.tight_layout()
    if show:
        plt.show()
    return {"figure": fig, "axis": axis, "bmax": bmax}


def format_case_commands(case_summary, amrvac_root, nproc=4):
    """Return portable setup/build/run commands for one staged case."""

    nproc = int(nproc)
    if nproc < 1:
        raise ValueError("nproc must be positive")
    case_dir = Path(case_summary["case_dir"])
    return [
        "export AMRVAC_DIR={}".format(quote(str(Path(amrvac_root)))),
        "cd {}".format(quote(str(case_dir))),
        '"$AMRVAC_DIR/setup.pl" -d=3',
        "make",
        "mpirun -np {} ./amrvac -i {} {}".format(
            nproc,
            quote(Path(case_summary["base_par"]).name),
            quote(Path(case_summary["override_par"]).name),
        ),
    ]


def _report_with_warnings(lines, warnings, prepend=False):
    warning_lines = ["WARNING: {}".format(warning) for warning in warnings]
    content = list(lines)
    if prepend:
        content = warning_lines + content
    else:
        content.extend(warning_lines)
    return "\n".join(content)


def _combined_adjusted_dimension(size, relaxation_grid, evolution_grid):
    for adjusted in range(int(size), 0, -1):
        if _compatible(adjusted, relaxation_grid) and _compatible(adjusted, evolution_grid):
            return adjusted
    raise ValueError("size {} is too small for the requested AMRVAC block constraints".format(size))


def _compatible(size, config):
    reducer = 2 ** (config.boundary_reduction_level - 1)
    amr_factor = 2 ** (config.amrvac_refinement_level - 1)
    return any(size % (reducer * amr_factor * block) == 0 for block in config.block_sizes)


def _exact_stage_plan(nx, ny, config):
    from .cases import recommend_amrvac_grid

    plan = recommend_amrvac_grid(
        nx,
        ny,
        boundary_reduction_level=config.boundary_reduction_level,
        amrvac_refinement_level=config.amrvac_refinement_level,
        block_sizes=config.block_sizes,
    )
    if plan["total_trim"]:
        raise RuntimeError("internal grid planning error")
    return plan


def _checked_window(window, preview_shape):
    ny, nx = int(preview_shape[0]), int(preview_shape[1])
    if window is None:
        selected = {"x0": 0, "y0": 0, "nx": nx, "ny": ny}
    else:
        selected = {key: int(window[key]) for key in ("x0", "y0", "nx", "ny")}
    if selected["x0"] < 0 or selected["y0"] < 0 or selected["nx"] < 1 or selected["ny"] < 1:
        raise ValueError("invalid window: {}".format(selected))
    if selected["x0"] + selected["nx"] > nx or selected["y0"] + selected["ny"] > ny:
        raise ValueError("window {} exceeds preview shape {}".format(selected, preview_shape))
    return selected


def _validate_product_plan(label, product, plan):
    meta = product.get("amrvac") or {}
    actual = (meta.get("nx_physical"), meta.get("ny_physical"))
    expected = (plan["x"]["boundary_physical"], plan["y"]["boundary_physical"])
    if actual != expected:
        raise RuntimeError("{} boundary shape {} does not match grid plan {}".format(label, actual, expected))


def _first_boundary(product):
    outputs = product.get("outputs", [])
    if not outputs:
        raise ValueError("boundary product does not contain output frames")
    return outputs[0]


def _normalize_nlfff_method(value):
    value = str(value).strip().lower().replace("-", "_")
    aliases = {
        "mfr": "legacy_mfr",
        "legacy_mf": "legacy_mfr",
        "magnetofriction": "legacy_mfr",
        "magnetofrictional_relaxation": "legacy_mfr",
        "opt": "optimization",
        "weighted_optimization": "optimization",
        "gr": "grad_rubin",
        "grad_rubin": "grad_rubin",
        "grad_rubin_nlfff": "grad_rubin",
        "potential": "potential",
    }
    return aliases.get(value, value)


def _write_grad_rubin_external_alpha_from_boundary(
    boundary_metadata,
    boundary_filename,
    output_dir,
    preset,
    alpha_options,
    unit_length_cm,
    unit_magneticfield_g,
    used_by_method,
    nlfff_method,
):
    import numpy as np

    from .alpha_cleaning import AlphaCleaningConfig, clean_grad_rubin_alpha, write_alpha_cleaning_audit
    from .external_alpha import write_external_alpha
    from .writers import ensure_output_dir, read_boundary_frame, write_json

    output_dir = ensure_output_dir(output_dir)
    frame = read_boundary_frame(boundary_filename)
    meta = boundary_metadata.get("amrvac") or {}
    unit_length_cm = float(unit_length_cm)
    unit_magneticfield_g = float(unit_magneticfield_g)
    dx_code = frame["dx"] * 1.0e5 / unit_length_cm
    dy_code = frame["dy"] * 1.0e5 / unit_length_cm
    xc_code = 0.5 * (float(meta["xprobmin1"]) + float(meta["xprobmax1"]))
    yc_code = 0.5 * (float(meta["xprobmin2"]) + float(meta["xprobmax2"]))
    x = xc_code + (np.arange(frame["nx"], dtype=float) - 0.5 * frame["nx"] + 0.5) * dx_code
    y = yc_code + (np.arange(frame["ny"], dtype=float) - 0.5 * frame["ny"] + 0.5) * dy_code
    options = dict(alpha_options or {})
    config = AlphaCleaningConfig(
        preset,
        derivative_spacing_x=dx_code,
        derivative_spacing_y=dy_code,
        coordinate_unit_cm=unit_length_cm,
        unit_length_cm=unit_length_cm,
        **options
    )
    product = clean_grad_rubin_alpha(frame["bx"], frame["by"], frame["bz"], config=config)
    audit_path = write_alpha_cleaning_audit(
        output_dir / "grad_rubin_alpha_cleaning_audit.json",
        product["audit"],
    )
    audit_sha256 = _file_sha256(audit_path)
    preprocess_info = boundary_metadata.get("preprocess")
    if isinstance(preprocess_info, dict):
        preprocess_enabled = bool(preprocess_info.get("enabled", False))
    else:
        preprocess_enabled = bool(preprocess_info)
    alpha_path = write_external_alpha(
        output_dir / "grad_rubin_external_alpha_v1.dat",
        product["alpha_clean"],
        x,
        y,
        alpha_raw=product["alpha_raw"],
        weight=product["weight"],
        pil_mask=product["pil_mask"],
        valid_mask=product["valid_mask"],
        polarity_mask=product["polarity_mask"],
        unit_length_cm=unit_length_cm,
        unit_magneticfield_g=unit_magneticfield_g,
        coordinate_unit="code_length",
        coordinate_unit_cm=unit_length_cm,
        metadata={
            "created_by": "DataConstrainWorkflow.prepare_initial_field",
            "source_boundary": str(Path(boundary_filename).expanduser().resolve()),
            "source_boundary_sha256": _file_sha256(boundary_filename),
            "preprocessing_audit": (boundary_metadata.get("preprocessing") or {}).get("audit_files"),
            "preprocess": {
                "enabled": preprocess_enabled,
                "mode": (boundary_metadata.get("preprocessing") or {}).get("mode"),
            },
            "alpha_cleaning_audit": str(audit_path),
            "alpha_cleaning_audit_sha256": audit_sha256,
            "nlfff_method": nlfff_method,
            "used_by_grad_rubin": bool(used_by_method),
            "unused": not bool(used_by_method),
        },
    )
    summary = {
        "external_alpha": str(alpha_path),
        "audit": str(audit_path),
        "audit_sha256": audit_sha256,
        "file_sha256": _file_sha256(alpha_path),
        "preset": str(preset),
        "used_by_method": bool(used_by_method),
        "nlfff_method": str(nlfff_method),
        "support": product["audit"].get("support", {}),
    }
    write_json(output_dir / "grad_rubin_external_alpha_manifest.json", summary)
    return summary


def _file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).expanduser().resolve().open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _one_shot_nlfff_restart_file(case_summary):
    base = Path(str(case_summary["base_filename"]))
    case_dir = Path(case_summary["case_dir"]).expanduser().resolve()
    if not base.is_absolute():
        base = case_dir / base
    return str((base.parent / (base.name + "0000.dat")).resolve())


def _initial_field_summary(project_dir, relaxation, evolution, potential, nlfff, nlfff_method="legacy_mfr", alpha=None):
    def boundary_summary(product):
        meta = product.get("amrvac") or {}
        return {
            "output_dir": product.get("output_dir"),
            "frame": Path(_first_boundary(product)).name,
            "physical_grid": (meta.get("nx_physical"), meta.get("ny_physical")),
            "file_grid": (meta.get("nx"), meta.get("ny")),
        }

    return {
        "project_dir": str(project_dir),
        "relaxation_boundary": boundary_summary(relaxation),
        "evolution_boundary": boundary_summary(evolution),
        "potential_case_dir": potential.get("case_dir") if potential else None,
        "mfr_case_dir": nlfff.get("case_dir") if nlfff and nlfff_method == "legacy_mfr" else None,
        "nlfff_method": nlfff_method,
        "nlfff_case_dir": nlfff.get("case_dir") if nlfff else None,
        "nlfff_metrics_csv": nlfff.get("nlfff_metrics_csv") if nlfff else None,
        "method_metrics_csv": nlfff.get("method_metrics_csv") if nlfff else None,
        "external_alpha": alpha.get("external_alpha") if alpha else None,
    }


__all__ = [
    "check_data_driven_dependencies",
    "DataConstrainWorkflow",
    "StageGridConfig",
    "create_data_constrain_workflow",
    "format_case_commands",
    "plan_data_constrain_region",
    "plan_workflow_region",
    "plot_data_constrain_region",
    "plot_magnetic_input_preview",
]
