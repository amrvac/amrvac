"""Notebook-facing orchestration for B-only time-dependent workflows."""

from __future__ import print_function

import hashlib
import json
import re
from pathlib import Path

from .workflow import (
    DataConstrainWorkflow,
    StageGridConfig,
    format_case_commands,
    plan_workflow_region,
    _normalize_nlfff_method,
    _one_shot_nlfff_restart_file,
)


class DataDrivenWorkflow(object):
    """Coordinate Potential -> optional MFR -> TMF/DataDriven evolution."""

    def __init__(
        self,
        amrvac_root,
        project_dir,
        input_dir,
        relaxation_grid=None,
        evolution_grid=None,
    ):
        self._base = DataConstrainWorkflow(
            amrvac_root=amrvac_root,
            project_dir=project_dir,
            input_dir=input_dir,
            relaxation_grid=relaxation_grid,
            evolution_grid=evolution_grid,
            manifest_workflow="data_driven",
        )
        self.initial_field = "legacy_mfr"
        self.nlfff_method = "legacy_mfr"
        self.evolution_mode = "tmf"
        self.sequence_selection = None
        self._sequence_state = None
        self._restart_state = None
        self._evolution_state = None

    @property
    def amrvac_root(self):
        return self._base.amrvac_root

    @property
    def project_dir(self):
        return self._base.project_dir

    @property
    def input_dir(self):
        return self._base.input_dir

    @property
    def relaxation_grid(self):
        return self._base.relaxation_grid

    @property
    def evolution_grid(self):
        return self._base.evolution_grid

    @property
    def region(self):
        return self._base.region

    @property
    def paths(self):
        paths = dict(self._base.paths)
        paths.update({
            "tmf": self.project_dir / "TimeDependentMagnetofriction",
            "data_driven": self.project_dir / "DataDriven",
        })
        return paths

    def inspect_input(self, refresh=False):
        return self._base.inspect_input(refresh=refresh)

    def preview_input(self, *args, **kwargs):
        return self._base.preview_input(*args, **kwargs)

    def input_report(self):
        return self._base.input_report()

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
        """Resolve and plan the input-appropriate region selection.

        ``auto`` selects a fixed CEA patch for raw HMI and a pixel window for
        an existing SHARP/CEA sequence.  Explicit ``fixed_cea`` is available
        for either input kind, while explicit ``pixel_window`` remains invalid
        for raw HMI.  ``cea_patch`` and ``window`` are retained as backwards-
        compatible aliases for ``raw_hmi_cea_patch`` and
        ``sharp_pixel_window`` respectively.
        """

        if preview_shape is None:
            preview_shape = self._base._require_preview()["preview"]["data"].shape
        result = plan_workflow_region(
            input_info=self.inspect_input()["info"],
            preview_shape=preview_shape,
            cea_patch=cea_patch,
            window=window,
            auto_trim=auto_trim,
            region_mode=region_mode,
            raw_hmi_cea_patch=raw_hmi_cea_patch,
            sharp_pixel_window=sharp_pixel_window,
            relaxation_grid=self.relaxation_grid,
            evolution_grid=self.evolution_grid,
        )
        self._base.region = result
        self._base._update_manifest()
        self._update_manifest({"region_planned": True})
        return result

    def plot_region(self, *args, **kwargs):
        return self._base.plot_region(*args, **kwargs)

    def region_report(self):
        return self._base.region_report()

    def configure_sequence(
        self,
        start_snapshot_index=0,
        stop_snapshot_index=None,
        frame_stride=1,
        initial_field=None,
        nlfff_method="legacy_mfr",
        evolution_mode="tmf",
    ):
        """Validate the time range and bind its first frame to the initial field."""

        info = self.inspect_input()["info"]
        count = _input_frame_count(info)
        start = int(start_snapshot_index)
        stop = count if stop_snapshot_index is None else int(stop_snapshot_index)
        stride = int(frame_stride)
        if start < 0 or start >= count:
            raise IndexError("START_SNAPSHOT_INDEX is outside 0..{}".format(count - 1))
        if stop <= start or stop > count:
            raise IndexError("STOP_SNAPSHOT_INDEX must be in {}..{}".format(start + 1, count))
        if stride < 1:
            raise ValueError("FRAME_STRIDE must be >= 1")
        indices = list(range(start, stop, stride))
        if len(indices) < 2:
            raise ValueError("time-dependent workflows require at least two selected frames")
        if initial_field is not None:
            legacy_initial = str(initial_field).strip().lower()
            if legacy_initial not in ("potential", "mfr"):
                raise ValueError("INITIAL_FIELD must be 'potential' or 'mfr'")
            if nlfff_method == "legacy_mfr":
                nlfff_method = "legacy_mfr" if legacy_initial == "mfr" else "potential"
        nlfff_method = _normalize_nlfff_method(nlfff_method)
        if nlfff_method not in ("legacy_mfr", "optimization", "grad_rubin", "potential"):
            raise ValueError("NLFFF_METHOD must be legacy_mfr, optimization, or grad_rubin")
        evolution_mode = str(evolution_mode).strip().lower()
        if evolution_mode not in ("tmf", "data_driven"):
            raise ValueError("EVOLUTION_MODE must be 'tmf' or 'data_driven'")
        self.initial_field = nlfff_method
        self.nlfff_method = nlfff_method
        self.evolution_mode = evolution_mode
        self.sequence_selection = {
            "start": start,
            "stop": stop,
            "stride": stride,
            "indices": indices,
        }
        self._update_manifest({
            "sequence_selection": {key: value for key, value in self.sequence_selection.items() if key != "indices"},
            "nlfff_method": nlfff_method,
            "evolution_mode": evolution_mode,
        })
        return self.sequence_selection

    def sequence_configuration_report(self):
        """Return a compact summary of the selected observation interval."""

        selection = self._require_selection()
        indices = selection["indices"]
        return "\n".join([
            "Selected frames: {}".format(len(indices)),
            "Snapshot indices: {} to {} (stride {})".format(
                indices[0], indices[-1], selection["stride"]
            ),
            "First selected frame defines observation time 0 s",
            "Selected NLFFF method: {}".format(self.nlfff_method),
            "Evolution mode: {}".format(self.evolution_mode),
        ])

    def prepare_reference_frame(
        self,
        cea_remap_options=None,
        geometry=None,
        preprocess=False,
        preprocessing_mode=None,
        fail_on_nonconvergence=False,
        nghost=2,
        quicklook=True,
        vmax=500.0,
        potential_options=None,
        mfr_options=None,
        optimization_options=None,
        grad_rubin_options=None,
        gr_alpha_cleaning=False,
        gr_alpha_preset="recommended",
        gr_alpha_options=None,
        unit_length_cm=1.0e9,
        unit_magneticfield_g=100.0,
        resume=True,
        overwrite=False,
    ):
        """Prepare the first selected frame and stage the requested initial cases."""

        selection = self._require_selection()
        info = self.inspect_input()["info"]
        reference_signature_payload = _reference_signature_payload(
            info=info,
            reference_index=selection["start"],
            region=self.region,
            relaxation_grid=self.relaxation_grid,
            evolution_grid=self.evolution_grid,
            cea_remap_options=cea_remap_options or {},
            geometry=geometry,
            preprocess=preprocess,
            preprocessing_mode=preprocessing_mode,
            fail_on_nonconvergence=fail_on_nonconvergence,
            nghost=nghost,
            nlfff_method=self.nlfff_method,
            potential_options=potential_options or {},
            mfr_options=mfr_options or {},
            optimization_options=optimization_options or {},
            grad_rubin_options=grad_rubin_options or {},
            gr_alpha_cleaning=bool(gr_alpha_cleaning),
            gr_alpha_preset=gr_alpha_preset,
            gr_alpha_options=gr_alpha_options or {},
            unit_length_cm=unit_length_cm,
            unit_magneticfield_g=unit_magneticfield_g,
        )
        reference_signature = _payload_signature(reference_signature_payload)
        reference_cache = self.paths["prepared_magnetograms"] / "reference_cache.json"
        cached_reference = _read_json(reference_cache)
        if cached_reference and cached_reference.get("signature") != reference_signature and not overwrite:
            raise RuntimeError(
                "reference-frame settings changed; use a new project directory or set "
                "OVERWRITE_PREPARED_SEQUENCE=True"
            )
        if resume and not overwrite and cached_reference:
            cached_state = cached_reference.get("state")
            if cached_reference.get("signature") == reference_signature and _initial_state_exists(cached_state):
                self._base._initial_field_state = cached_state
                self._update_manifest({
                    "reference_frame": {
                        "snapshot_index": selection["start"],
                        "source_frame": cached_state["source_frame"],
                        "nlfff_method": self.nlfff_method,
                        "signature": reference_signature,
                        "reused": True,
                    }
                })
                return cached_state
        reference_index = selection["start"]
        original_input_dir = self._base.input_dir
        original_input_state = self._base._input_state
        if info["kind"] == "cea_vector_sequence" and self.region["selection_mode"] == "fixed_cea_patch":
            from .pipeline import make_fixed_cea_sequence

            fixed_dir = self.paths["prepared_magnetograms"] / "ReferenceFixedCEA"
            make_fixed_cea_sequence(
                self.input_dir,
                fixed_dir,
                self.region["cea_patch"],
                indices=[reference_index],
                progress=False,
                resume=True,
            )
            self._base.input_dir = Path(fixed_dir)
            self._base._input_state = None
            reference_index = 0
        try:
            result = self._base.prepare_initial_field(
                snapshot_index=reference_index,
                cea_remap_options=cea_remap_options,
                geometry=geometry,
                preprocess=preprocess,
                preprocessing_mode=preprocessing_mode,
                fail_on_nonconvergence=fail_on_nonconvergence,
                nghost=nghost,
                quicklook=quicklook,
                vmax=vmax,
                stage_potential=True,
                nlfff_method=self.nlfff_method,
                potential_options=potential_options,
                mfr_options=mfr_options,
                optimization_options=optimization_options,
                grad_rubin_options=grad_rubin_options,
                gr_alpha_cleaning=gr_alpha_cleaning,
                gr_alpha_preset=gr_alpha_preset,
                gr_alpha_options=gr_alpha_options,
                unit_length_cm=unit_length_cm,
                unit_magneticfield_g=unit_magneticfield_g,
            )
        finally:
            self._base.input_dir = original_input_dir
            self._base._input_state = original_input_state
        result["snapshot_index"] = selection["start"]
        result["signature"] = reference_signature
        from .writers import ensure_output_dir, write_json
        ensure_output_dir(reference_cache.parent)
        write_json(reference_cache, {
            "signature": reference_signature,
            "signature_payload": reference_signature_payload,
            "state": result,
        })
        self._update_manifest({
            "reference_frame": {
                "snapshot_index": selection["start"],
                "source_frame": result["source_frame"],
                "nlfff_method": self.nlfff_method,
                "signature": reference_signature,
                "reused": False,
            }
        })
        return result

    def prepare_magnetic_sequence(
        self,
        cea_remap_options=None,
        geometry=None,
        preprocess=False,
        preprocessing_mode=None,
        fail_on_nonconvergence=False,
        nghost=2,
        quicklook=True,
        vmax=500.0,
        resume=True,
        overwrite=False,
        progress=True,
    ):
        """Prepare the selected full sequence and write high-resolution B frames."""

        from .fits_io import inspect_magnetic_input, summarize_magnetic_input
        from .pipeline import make_cea_sequence, prepare_boundary_sequence
        from .writers import ensure_output_dir, write_json

        selection = self._require_selection()
        info = self.inspect_input()["info"]
        indices = selection["indices"]
        remap_options = dict(cea_remap_options or {})
        signature_payload = _sequence_signature_payload(
            info, indices, self.region, remap_options, self.evolution_grid,
            preprocess, nghost,
            preprocessing_mode=preprocessing_mode,
            fail_on_nonconvergence=fail_on_nonconvergence,
        )
        signature = _payload_signature(signature_payload)
        cache_root = ensure_output_dir(self.paths["prepared_magnetograms"] / "Sequence")
        cache_pointer = cache_root / "sequence_cache.json"
        cached = _read_json(cache_pointer)
        if cached and cached.get("signature") != signature and not overwrite:
            raise RuntimeError(
                "prepared-sequence settings changed; use a new project directory or set "
                "OVERWRITE_PREPARED_SEQUENCE=True"
            )

        output_dir = self.paths["boundaries"] / "Evolution" / "Sequence" / signature[:16]
        if resume and not overwrite and cached and cached.get("signature") == signature:
            metadata = _read_json(output_dir / "boundary_parameters.json")
            if _boundary_sequence_exists(metadata, expected_count=len(indices)):
                warnings = summarize_magnetic_input(info).get("warnings", [])
                state = {
                    "signature": signature,
                    "signature_payload": signature_payload,
                    "boundary_input_dir": cached.get("converted_dir", str(self.input_dir)),
                    "boundary": metadata,
                    "summary": _sequence_summary(metadata, warnings),
                    "reused": True,
                }
                self._sequence_state = state
                self._update_manifest({
                    "sequence": {
                        "signature": signature,
                        "output_dir": metadata["output_dir"],
                        "frame_count": metadata["output_frame_count"],
                        "selection": {key: value for key, value in selection.items() if key != "indices"},
                        "reused": True,
                    }
                })
                return state

        boundary_input_dir = self.input_dir
        boundary_indices = indices
        if info["kind"] == "raw_hmi_vector_sequence":
            selected_frames = [info["raw_hmi_sequence"]["complete"][index] for index in indices]
            converted_dir = cache_root / signature[:16]
            raw_input = {"raw_hmi_sequence": {"complete": selected_frames}}
            options = dict(remap_options)
            options.setdefault("components_only", True)
            options.setdefault("quicklook_first", True)
            make_cea_sequence(
                raw_input,
                converted_dir,
                progress=progress,
                resume=bool(resume and not overwrite),
                **dict(self.region["cea_patch"], **options)
            )
            converted = inspect_magnetic_input(converted_dir)
            if converted["kind"] != "cea_vector_sequence":
                raise RuntimeError("raw HMI sequence did not produce a complete Br/Bt/Bp sequence")
            boundary_input_dir = converted_dir
            boundary_indices = list(range(len(selected_frames)))
        elif info["kind"] == "cea_vector_sequence" and self.region["selection_mode"] == "fixed_cea_patch":
            from .pipeline import make_fixed_cea_sequence

            converted_dir = cache_root / signature[:16]
            make_fixed_cea_sequence(
                self.input_dir,
                converted_dir,
                self.region["cea_patch"],
                indices=indices,
                progress=progress,
                resume=bool(resume and not overwrite),
            )
            boundary_input_dir = converted_dir
            boundary_indices = list(range(len(indices)))
        elif info["kind"] != "cea_vector_sequence":
            raise ValueError("DataDriven requires a Br/Bt/Bp or raw-HMI sequence")

        metadata = prepare_boundary_sequence(
            boundary_input_dir,
            output_dir,
            window=self.region["window"],
            level=self.evolution_grid.boundary_reduction_level,
            geometry=geometry,
            preprocess=preprocess,
            preprocessing_mode=preprocessing_mode,
            fail_on_nonconvergence=fail_on_nonconvergence,
            nghost=nghost,
            quicklook=quicklook,
            vmax=vmax,
            indices=boundary_indices,
            require_timestamps=True,
            progress=progress,
        )
        from .workflow import _validate_product_plan
        _validate_product_plan("evolution", metadata, self.region["plans"]["evolution"])
        warnings = summarize_magnetic_input(info).get("warnings", [])
        summary = _sequence_summary(metadata, warnings)
        state = {
            "signature": signature,
            "signature_payload": signature_payload,
            "boundary_input_dir": str(boundary_input_dir),
            "boundary": metadata,
            "summary": summary,
            "reused": False,
        }
        write_json(cache_pointer, {
            "signature": signature,
            "converted_dir": str(boundary_input_dir),
            "selection": selection,
        })
        self._sequence_state = state
        self._update_manifest({
            "sequence": {
                "signature": signature,
                "output_dir": metadata["output_dir"],
                "frame_count": metadata["output_frame_count"],
                "selection": {key: value for key, value in selection.items() if key != "indices"},
                "reused": False,
            }
        })
        return state

    def resolve_initial_restart(
        self,
        restart_file=None,
        selected_restart_number=None,
        show_plot=True,
        plot_lorentz_force=False,
    ):
        """Resolve the selected NLFFF restart for time evolution."""

        if self._base._initial_field_state is None:
            raise RuntimeError("prepare_reference_frame must be called first")
        initial = self._base._initial_field_state
        if restart_file is not None:
            state = {"restart_file": str(Path(restart_file).expanduser().resolve()), "reason": "explicit restart"}
        elif self.nlfff_method == "potential":
            state = {
                "restart_file": str(Path(initial["potential_restart_file"]).expanduser().resolve()),
                "reason": "legacy potential-only initial field selected",
            }
        elif self.nlfff_method == "legacy_mfr":
            from .diagnostics import (
                find_relaxation_restart_snapshots,
                normalize_nlfff_metrics,
                plot_relaxation_diagnostics,
                read_relaxation_diagnostics,
                select_relaxation_restart,
            )
            mfr_case = initial.get("mfr_case")
            if mfr_case is None:
                raise RuntimeError("NLFFF_METHOD='legacy_mfr' requires a staged MFR case")
            output_dir = self.paths["mfr"] / "output"
            metrics = normalize_nlfff_metrics(mfr_case)
            diagnostics = read_relaxation_diagnostics(metrics["path"])
            ditsave = int(mfr_case.get("mf_ditsave", 20000))
            restart_summary = find_relaxation_restart_snapshots(
                output_dir=output_dir,
                base_filename="data_driven_mfr",
                diagnostics=diagnostics,
                mf_ditsave=ditsave,
            )
            selection = select_relaxation_restart(
                diagnostics=diagnostics,
                output_dir=output_dir,
                base_filename="data_driven_mfr",
                mf_ditsave=ditsave,
                snapshot_number=selected_restart_number,
            )
            selected = selection.get("selected_marker") or {}
            markers = [dict(marker) for marker in restart_summary["markers"]]
            for marker in markers:
                marker["selected"] = marker.get("snapshot_number") == selected.get("snapshot_number")
            columns = ["cw_sin_theta"]
            if plot_lorentz_force:
                columns.append("lorentz_force")
            plot_relaxation_diagnostics(
                diagnostics,
                output_path=output_dir / "relaxation_diagnostics_quicklook.png",
                columns=columns,
                x_column="iteration",
                restart_markers=markers,
            )
            if show_plot:
                import matplotlib.pyplot as plt
                plt.show()
            state = dict(selection)
            state["warnings"] = restart_summary.get("warnings", []) + selection.get("warnings", [])
            state["metrics"] = metrics
        else:
            from .diagnostics import normalize_nlfff_metrics

            nlfff_case = initial.get("nlfff_case")
            if nlfff_case is None:
                raise RuntimeError("NLFFF_METHOD='{}' requires a staged case".format(self.nlfff_method))
            metrics = normalize_nlfff_metrics(nlfff_case)
            state = {
                "restart_file": _one_shot_nlfff_restart_file(nlfff_case),
                "reason": "{} one-shot NLFFF output selected".format(self.nlfff_method),
                "warnings": list(metrics.get("warnings", [])),
                "metrics": metrics,
            }
        if not state.get("restart_file"):
            raise RuntimeError("no initial restart was resolved")
        self._restart_state = state
        self._update_manifest({"initial_restart": state})
        return state

    def stage_evolution(self, mode=None, driving_time_scale=12.0, case_options=None):
        """Stage exactly one TMF or DataDriven case from the prepared sequence."""

        if self._sequence_state is None:
            raise RuntimeError("prepare_magnetic_sequence must be called first")
        if self._restart_state is None:
            raise RuntimeError("resolve_initial_restart must be called first")
        from .cases import stage_data_driven_case, stage_time_dependent_magnetofriction_case

        mode = self.evolution_mode if mode is None else str(mode).strip().lower()
        if mode not in ("tmf", "data_driven"):
            raise ValueError("mode must be 'tmf' or 'data_driven'")
        options = dict(case_options or {})
        options.setdefault("portable_paths", True)
        plan = self.region["plans"]["evolution"]
        common = dict(
            sequence_metadata=self._sequence_state["boundary"],
            restart_file=self._restart_state["restart_file"],
            amrvac_root=self.amrvac_root,
            driving_time_scale=driving_time_scale,
            refine_max_level=self.evolution_grid.amrvac_refinement_level,
            block_nx1=plan["x"]["block_size"],
            block_nx2=plan["y"]["block_size"],
            block_nx3=plan["y"]["block_size"],
        )
        if mode == "tmf":
            case = stage_time_dependent_magnetofriction_case(
                case_dir=self.paths["tmf"], **dict(common, **options)
            )
        else:
            case = stage_data_driven_case(
                case_dir=self.paths["data_driven"], **dict(common, **options)
            )
        self.evolution_mode = mode
        self._evolution_state = {"mode": mode, "case": case}
        self._update_manifest({
            "evolution": {
                "mode": mode,
                "case_dir": case["case_dir"],
                "driving_time_scale": float(driving_time_scale),
                "restart_file": self._restart_state["restart_file"],
            }
        })
        return self._evolution_state

    def reference_report(self):
        return self._base.initial_field_report()

    def initial_field_commands(self, nproc=4):
        return self._base.initial_field_commands(nproc=nproc)

    def normalize_initial_nlfff_metrics(self, overwrite=False):
        """Expose the shared metrics normalizer to the public notebook."""

        return self._base.normalize_initial_nlfff_metrics(overwrite=overwrite)

    def initial_nlfff_metrics_report(self):
        return self._base.initial_nlfff_metrics_report()

    def plot_initial_nlfff_metrics(self, show=False, output_path=None):
        """Create the same method-neutral metrics quicklook as DataConstrain."""

        return self._base.plot_initial_nlfff_metrics(
            show=show,
            output_path=output_path,
        )

    def sequence_report(self):
        if self._sequence_state is None:
            raise RuntimeError("prepare_magnetic_sequence must be called first")
        summary = self._sequence_state["summary"]
        lines = [
            "Frames written: {}".format(summary["frame_count"]),
            "Observation time: {:.1f} to {:.1f} s".format(summary["start_time"], summary["end_time"]),
            "Cadence min/median/max: {:.1f} / {:.1f} / {:.1f} s".format(
                summary["cadence_min"], summary["cadence_median"], summary["cadence_max"]
            ),
            "Boundary directory: {}".format(self._sequence_state["boundary"]["output_dir"]),
        ]
        lines.extend("WARNING: {}".format(item) for item in summary.get("warnings", []))
        return "\n".join(lines)

    def plot_sequence_quicklook(
        self,
        show=True,
        vmax=None,
        snapshot_count=9,
        columns=3,
    ):
        """Plot evenly sampled Bz snapshots in a compact image grid."""

        if self._sequence_state is None:
            raise RuntimeError("prepare_magnetic_sequence must be called first")
        metadata = self._sequence_state["boundary"]
        outputs = [Path(path) for path in metadata.get("outputs", [])]
        selected_indices = _representative_snapshot_indices(
            len(outputs), snapshot_count=snapshot_count
        )
        if not selected_indices:
            return {"paths": [], "frame_paths": [], "figure": None}
        from .writers import read_boundary_frame
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        frames = [read_boundary_frame(outputs[index]) for index in selected_indices]
        if vmax is None:
            vmax = float(metadata.get("quicklook_vmax", 500.0))
        vmax = float(vmax)
        columns = int(columns)
        if columns < 1:
            raise ValueError("columns must be >= 1")
        rows = (len(frames) + columns - 1) // columns
        fig, axes = plt.subplots(
            rows, columns,
            figsize=(5.2 * columns, 3.6 * rows),
            squeeze=False,
        )
        flat_axes = list(axes.ravel())
        frame_metadata = metadata.get("frames", [])
        for position, (axis, frame, snapshot_index) in enumerate(zip(
            flat_axes, frames, selected_indices
        )):
            data = frame["bz"]
            image = axis.imshow(
                data, origin="lower", cmap="gray", vmin=-vmax, vmax=vmax
            )
            axis.set_box_aspect(float(data.shape[0]) / float(data.shape[1]))
            selected_metadata = (
                frame_metadata[snapshot_index]
                if snapshot_index < len(frame_metadata)
                else {}
            )
            axis.set_title(
                "Frame {}/{} — Br\n{}".format(
                    snapshot_index + 1,
                    len(outputs),
                    _observation_time_label(selected_metadata, frame),
                )
            )
            axis.set_xlabel("x pixel")
            axis.set_ylabel("y pixel")
            is_row_end = (position % columns == columns - 1) or position == len(frames) - 1
            if is_row_end:
                color_axis = inset_axes(
                    axis,
                    width="3.5%",
                    height="100%",
                    loc="lower left",
                    bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                    bbox_transform=axis.transAxes,
                    borderpad=0,
                )
                colorbar = fig.colorbar(image, cax=color_axis)
                colorbar.set_label("G")
        for axis in flat_axes[len(frames):]:
            axis.set_visible(False)
        fig.subplots_adjust(wspace=0.42, hspace=0.55)
        if show:
            plt.show()
        output_dir = Path(metadata["output_dir"])
        legacy_paths = [
            output_dir / "boundary_{}_quicklook.png".format(label)
            for label in ("first", "middle", "last")
        ]
        return {
            "paths": [str(path) for path in legacy_paths if path.exists()],
            "frame_paths": [str(outputs[index]) for index in selected_indices],
            "snapshot_indices": selected_indices,
            "snapshot_count": len(selected_indices),
            "figure": fig,
        }

    def restart_report(self):
        if self._restart_state is None:
            raise RuntimeError("resolve_initial_restart must be called first")
        lines = [
            "Initial restart: {}".format(self._restart_state["restart_file"]),
            "Reason: {}".format(self._restart_state.get("reason", "selected")),
        ]
        metrics = self._restart_state.get("metrics")
        if metrics:
            lines.append("Unified NLFFF metrics: {}".format(metrics.get("path")))
            lines.append("Metrics status: {}".format(metrics.get("status")))
        return "\n".join(lines)

    def evolution_report(self):
        if self._evolution_state is None:
            raise RuntimeError("stage_evolution must be called first")
        case = self._evolution_state["case"]
        return "Evolution mode: {}\nCase directory: {}".format(
            self._evolution_state["mode"], case["case_dir"]
        )

    def evolution_commands(self, nproc=4):
        if self._evolution_state is None:
            raise RuntimeError("stage_evolution must be called first")
        return "\n".join(format_case_commands(self._evolution_state["case"], self.amrvac_root, nproc=nproc))

    def _require_selection(self):
        if self.sequence_selection is None:
            raise RuntimeError("configure_sequence must be called first")
        if self.region is None:
            raise RuntimeError("plan_region must be called first")
        return self.sequence_selection

    def _update_manifest(self, updates):
        from .manifest import update_project_manifest

        shared = {
            "input_dir": str(self.input_dir),
            "project_dir": str(self.project_dir),
            "paths": {key: str(value) for key, value in self.paths.items()},
            "relaxation_grid": self.relaxation_grid.as_dict(),
            "evolution_grid": self.evolution_grid.as_dict(),
            "region": self.region,
        }
        return update_project_manifest(
            self.paths["manifest"], shared=shared,
            workflow="data_driven", workflow_updates=updates,
        )


def create_data_driven_workflow(*args, **kwargs):
    return DataDrivenWorkflow(*args, **kwargs)


def _input_frame_count(info):
    if info["kind"] == "raw_hmi_vector_sequence":
        return len(info["raw_hmi_sequence"]["complete"])
    if info["kind"] == "cea_vector_sequence":
        return len(info["components"]["br"])
    raise ValueError("DataDriven requires a time sequence with at least two complete frames")


def _file_identity(path):
    path = Path(path)
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": int(stat.st_mtime_ns)}


def _sequence_signature_payload(
    info,
    indices,
    region,
    remap_options,
    grid,
    preprocess,
    nghost,
    preprocessing_mode=None,
    fail_on_nonconvergence=False,
):
    if info["kind"] == "raw_hmi_vector_sequence":
        selected = []
        for index in indices:
            frame = info["raw_hmi_sequence"]["complete"][index]
            selected.append({name: _file_identity(frame[name]) for name in ("field", "inclination", "azimuth", "disambig")})
    else:
        selected = []
        for index in indices:
            selected.append({name: _file_identity(info["components"][name][index]) for name in ("br", "bt", "bp")})
    return {
        "kind": info["kind"],
        "frames": selected,
        "region": region,
        "remap_options": remap_options,
        "evolution_grid": grid.as_dict(),
        "preprocess": bool(preprocess),
        "preprocessing_mode": preprocessing_mode,
        "fail_on_nonconvergence": bool(fail_on_nonconvergence),
        "nghost": int(nghost),
    }


def _reference_signature_payload(
    info,
    reference_index,
    region,
    relaxation_grid,
    evolution_grid,
    cea_remap_options,
    geometry,
    preprocess,
    preprocessing_mode,
    fail_on_nonconvergence,
    nghost,
    nlfff_method,
    potential_options,
    mfr_options,
    optimization_options,
    grad_rubin_options,
    gr_alpha_cleaning,
    gr_alpha_preset,
    gr_alpha_options,
    unit_length_cm,
    unit_magneticfield_g,
):
    if info["kind"] == "raw_hmi_vector_sequence":
        frame = info["raw_hmi_sequence"]["complete"][reference_index]
        source = {
            name: _file_identity(frame[name])
            for name in ("field", "inclination", "azimuth", "disambig")
        }
    elif info["kind"] == "cea_vector_sequence":
        source = {
            name: _file_identity(info["components"][name][reference_index])
            for name in ("br", "bt", "bp")
        }
    else:
        raise ValueError("DataDriven requires a Br/Bt/Bp or raw-HMI sequence")
    return {
        "source": source,
        "reference_index": int(reference_index),
        "region": region,
        "relaxation_grid": relaxation_grid.as_dict(),
        "evolution_grid": evolution_grid.as_dict(),
        "cea_remap_options": cea_remap_options,
        "geometry": geometry,
        "preprocess": bool(preprocess),
        "preprocessing_mode": preprocessing_mode,
        "fail_on_nonconvergence": bool(fail_on_nonconvergence),
        "nghost": int(nghost),
        "nlfff_method": str(nlfff_method),
        "potential_options": potential_options,
        "mfr_options": mfr_options,
        "optimization_options": optimization_options,
        "grad_rubin_options": grad_rubin_options,
        "gr_alpha_cleaning": bool(gr_alpha_cleaning),
        "gr_alpha_preset": str(gr_alpha_preset),
        "gr_alpha_options": gr_alpha_options,
        "unit_length_cm": float(unit_length_cm),
        "unit_magneticfield_g": float(unit_magneticfield_g),
    }


def _payload_signature(payload):
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _read_json(path):
    path = Path(path)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (TypeError, ValueError):
        return None


def _paths_from_outputs(outputs):
    if isinstance(outputs, dict):
        return [Path(value) for value in outputs.values() if value]
    return [Path(value) for value in (outputs or []) if value]


def _initial_state_exists(state):
    if not isinstance(state, dict):
        return False
    for key in ("relaxation_boundary", "evolution_boundary"):
        metadata = state.get(key) or {}
        if not _paths_from_outputs(metadata.get("outputs")):
            return False
        if not all(path.exists() for path in _paths_from_outputs(metadata.get("outputs"))):
            return False
    for key in ("potential_case", "mfr_case", "nlfff_case"):
        case = state.get(key)
        if case is None:
            continue
        required = [case.get("base_par"), case.get("override_par")]
        if not all(value and Path(value).exists() for value in required):
            return False
    return True


def _boundary_sequence_exists(metadata, expected_count):
    if not isinstance(metadata, dict):
        return False
    if int(metadata.get("output_frame_count", 0)) != int(expected_count):
        return False
    outputs = _paths_from_outputs(metadata.get("outputs"))
    if len(outputs) != int(expected_count):
        return False
    return all(path.is_file() for path in outputs)


def _sequence_summary(metadata, warnings):
    import numpy as np

    times = np.asarray([frame["snapshot_time"] for frame in metadata["frames"]], dtype=float)
    cadence = np.diff(times)
    if len(times) < 2 or np.any(cadence <= 0.0):
        raise ValueError("output boundary times must be strictly increasing")
    return {
        "frame_count": len(times),
        "start_time": float(times[0]),
        "end_time": float(times[-1]),
        "cadence_min": float(np.min(cadence)),
        "cadence_median": float(np.median(cadence)),
        "cadence_max": float(np.max(cadence)),
        "warnings": list(warnings),
    }


def _representative_snapshot_indices(frame_count, snapshot_count=9):
    """Return uniformly spaced indices including the first and last frame."""

    frame_count = int(frame_count)
    if frame_count <= 0:
        return []
    snapshot_count = int(snapshot_count)
    if snapshot_count < 1:
        raise ValueError("snapshot_count must be >= 1")
    selected_count = min(frame_count, snapshot_count)
    if selected_count == 1:
        return [0]
    return [
        int(round(index * (frame_count - 1) / float(selected_count - 1)))
        for index in range(selected_count)
    ]


def _observation_time_label(frame_metadata, boundary_frame):
    """Format a real observation time, with relative time as a fallback."""

    candidates = [frame_metadata.get("observation_time")]
    source = frame_metadata.get("source") or {}
    candidates.extend(source.get(name) for name in ("Br", "br", "Bt", "bt"))
    pattern = re.compile(
        r"(\d{4})[.\-]?(\d{2})[.\-]?(\d{2})[_T ]"
        r"(\d{2}):?(\d{2}):?(\d{2})",
        re.IGNORECASE,
    )
    for candidate in candidates:
        if not candidate:
            continue
        text = str(candidate)
        match = pattern.search(text)
        if match:
            year, month, day, hour, minute, second = match.groups()
            upper = text.upper()
            scale = " TAI" if "TAI" in upper else (" UTC" if "UTC" in upper else "")
            return "{}-{}-{} {}:{}:{}{}".format(
                year, month, day, hour, minute, second, scale
            )
    return "t = {:.1f} s".format(boundary_frame["snapshot_time"])


__all__ = ["DataDrivenWorkflow", "create_data_driven_workflow"]
