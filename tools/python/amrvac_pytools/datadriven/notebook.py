"""Small public configuration layer for the two shipped notebooks.

The notebooks should present user choices, not the internal case names and
manifest/cache bookkeeping used by the workflow classes.  This module keeps
that translation in one place while preserving the existing workflow APIs.
"""

from __future__ import print_function

import math


PUBLIC_NLFFF_METHODS = ("legacy_mfr", "optimization", "grad_rubin")
PUBLIC_PREPROCESSING_MODES = ("none", "conservative", "recommended")
PUBLIC_GR_ALPHA_PRESETS = ("none", "conservative", "recommended")
DEFAULT_BLOCK_SIZES = (12, 14, 16, 18, 20)


def normalize_notebook_options(
    nlfff_method="legacy_mfr",
    preprocess_vector=False,
    preprocessing_mode="recommended",
    gr_alpha_cleaning=False,
    gr_alpha_preset="recommended",
    unit_length_cm=1.0e9,
    unit_magneticfield_g=100.0,
):
    """Validate public switches and return their effective workflow values.

    This function deliberately does not enable either switch based on the
    selected NLFFF method.  When a switch is off, its configured preset is
    retained for auditability but the effective operation is disabled.
    """

    method = str(nlfff_method).strip().lower()
    if method not in PUBLIC_NLFFF_METHODS:
        raise ValueError(
            "NLFFF_METHOD must be one of: {}".format(", ".join(PUBLIC_NLFFF_METHODS))
        )
    mode = str(preprocessing_mode).strip().lower()
    if mode not in PUBLIC_PREPROCESSING_MODES:
        raise ValueError(
            "PREPROCESSING_MODE must be one of: {}".format(
                ", ".join(PUBLIC_PREPROCESSING_MODES)
            )
        )
    preset = str(gr_alpha_preset).strip().lower()
    if preset not in PUBLIC_GR_ALPHA_PRESETS:
        raise ValueError(
            "GR_ALPHA_PRESET must be one of: {}".format(
                ", ".join(PUBLIC_GR_ALPHA_PRESETS)
            )
        )
    length = float(unit_length_cm)
    field = float(unit_magneticfield_g)
    if not math.isfinite(length) or length <= 0.0:
        raise ValueError("UNIT_LENGTH_CM must be finite and positive")
    if not math.isfinite(field) or field <= 0.0:
        raise ValueError("UNIT_MAGNETICFIELD_G must be finite and positive")

    preprocess = bool(preprocess_vector)
    clean_alpha = bool(gr_alpha_cleaning)
    return {
        "nlfff_method": method,
        "preprocess_vector": preprocess,
        "preprocessing_mode": mode,
        "effective_preprocessing_mode": mode if preprocess else "none",
        "gr_alpha_cleaning": clean_alpha,
        "gr_alpha_preset": preset,
        "unit_length_cm": length,
        "unit_magneticfield_g": field,
        "advisories": _notebook_advisories(method, preprocess, clean_alpha),
    }


def notebook_configuration_report(options):
    """Return the compact user-facing switch summary and advisories."""

    options = dict(options or {})
    lines = [
        "Selected NLFFF method: {}".format(options.get("nlfff_method")),
        "Shared vector preprocessing: {}".format(
            "enabled" if options.get("preprocess_vector") else "disabled"
        ),
        "Grad-Rubin alpha cleaning: {}".format(
            "enabled" if options.get("gr_alpha_cleaning") else "disabled"
        ),
    ]
    lines.extend("ADVISORY: {}".format(item) for item in options.get("advisories", ()))
    return "\n".join(lines)


def create_notebook_workflow(
    workflow_kind,
    amrvac_root,
    project_dir,
    input_dir,
    relaxation_boundary_reduction_level,
    evolution_boundary_reduction_level,
    evolution_amrvac_refinement_level,
    block_sizes=DEFAULT_BLOCK_SIZES,
):
    """Create one public workflow with validated, compatible stage grids."""

    from .workflow import StageGridConfig

    relaxation_grid = StageGridConfig(
        relaxation_boundary_reduction_level, 1, block_sizes
    )
    evolution_grid = StageGridConfig(
        evolution_boundary_reduction_level,
        evolution_amrvac_refinement_level,
        block_sizes,
    )
    kind = str(workflow_kind).strip().lower()
    if kind == "data_constrain":
        from .workflow import create_data_constrain_workflow

        return create_data_constrain_workflow(
            amrvac_root=amrvac_root,
            project_dir=project_dir,
            input_dir=input_dir,
            relaxation_grid=relaxation_grid,
            evolution_grid=evolution_grid,
        )
    if kind == "data_driven":
        from .data_driven_workflow import create_data_driven_workflow

        return create_data_driven_workflow(
            amrvac_root=amrvac_root,
            project_dir=project_dir,
            input_dir=input_dir,
            relaxation_grid=relaxation_grid,
            evolution_grid=evolution_grid,
        )
    raise ValueError("workflow_kind must be 'data_constrain' or 'data_driven'")


def default_notebook_nlfff_run_options(write_detailed_history=False):
    """Return the small, method-neutral run configuration shown to users."""

    detailed = bool(write_detailed_history)
    return {
        "legacy_mfr": {
            "max_iterations": 100000,
            "save_every": 5000,
            "cc": 0.5,
            "cy": 0.2,
            "cdivb": 0.01,
        },
        "optimization": {"max_iterations": 1000},
        "grad_rubin": {"polarity": 1, "max_iterations": 50},
        "write_detailed_history": detailed,
    }


def notebook_boundary_options(options, nghost=2, fail_on_nonconvergence=False, vmax=500.0):
    """Translate public switch values to boundary-preparation keywords."""

    options = dict(options or {})
    return {
        "nghost": int(nghost),
        "preprocess": bool(options.get("preprocess_vector", False)),
        "preprocessing_mode": options.get("effective_preprocessing_mode", "none"),
        "fail_on_nonconvergence": bool(fail_on_nonconvergence),
        "vmax": float(vmax),
    }


def notebook_single_frame_options(
    options,
    snapshot_index,
    cea_remap_options,
    potential_options,
    nlfff_run_options=None,
    gr_alpha_options=None,
    boundary_options=None,
):
    """Build ``prepare_initial_field`` keywords without exposing internals."""

    return _initial_field_options(
        options=options,
        cea_remap_options=cea_remap_options,
        potential_options=potential_options,
        nlfff_run_options=nlfff_run_options,
        gr_alpha_options=gr_alpha_options,
        boundary_options=boundary_options,
        snapshot_index=int(snapshot_index),
    )


def notebook_reference_frame_options(
    options,
    cea_remap_options,
    potential_options,
    nlfff_run_options=None,
    gr_alpha_options=None,
    boundary_options=None,
    resume=True,
    overwrite=True,
):
    """Build ``prepare_reference_frame`` keywords without exposing internals."""

    result = _initial_field_options(
        options=options,
        cea_remap_options=cea_remap_options,
        potential_options=potential_options,
        nlfff_run_options=nlfff_run_options,
        gr_alpha_options=gr_alpha_options,
        boundary_options=boundary_options,
    )
    # DataDrivenWorkflow binds the selected method in configure_sequence();
    # unlike DataConstrainWorkflow.prepare_initial_field(), its reference-frame
    # method therefore does not accept an nlfff_method keyword.
    result.pop("nlfff_method", None)
    result.update({"resume": bool(resume), "overwrite": bool(overwrite)})
    return result


def notebook_sequence_options(
    cea_remap_options,
    boundary_options,
    progress=True,
):
    """Build sequence-preparation keywords with cache policy kept internal."""

    result = {
        "cea_remap_options": dict(cea_remap_options or {}),
        "resume": True,
        "overwrite": True,
        "progress": bool(progress),
    }
    result.update(dict(boundary_options or {}))
    return result


def _initial_field_options(
    options,
    cea_remap_options,
    potential_options,
    nlfff_run_options,
    gr_alpha_options,
    boundary_options,
    snapshot_index=None,
):
    options = dict(options or {})
    run = _merge_run_options(nlfff_run_options)
    result = {
        "cea_remap_options": dict(cea_remap_options or {}),
        "potential_options": dict(potential_options or {}),
        "mfr_options": run["mfr_options"],
        "optimization_options": run["optimization_options"],
        "grad_rubin_options": run["grad_rubin_options"],
        "nlfff_method": options["nlfff_method"],
        "gr_alpha_cleaning": options["gr_alpha_cleaning"],
        "gr_alpha_preset": options["gr_alpha_preset"],
        "gr_alpha_options": dict(gr_alpha_options or {}),
        "unit_length_cm": options["unit_length_cm"],
        "unit_magneticfield_g": options["unit_magneticfield_g"],
    }
    result.update(dict(boundary_options or {}))
    if snapshot_index is not None:
        result["snapshot_index"] = int(snapshot_index)
    return result


def _merge_run_options(run_options):
    public = default_notebook_nlfff_run_options()
    if run_options:
        for method in ("legacy_mfr", "optimization", "grad_rubin"):
            public[method].update(dict(run_options.get(method, {})))
        if "write_detailed_history" in run_options:
            public["write_detailed_history"] = bool(run_options["write_detailed_history"])

    mfr = public["legacy_mfr"]
    optimization = public["optimization"]
    grad_rubin = public["grad_rubin"]
    detailed = public["write_detailed_history"]
    return {
        "mfr_options": {
            "mf_it_max": int(mfr["max_iterations"]),
            "mf_ditsave": int(mfr["save_every"]),
            "mf_cc": float(mfr["cc"]),
            "mf_cy": float(mfr["cy"]),
            "mf_cdivb": float(mfr["cdivb"]),
            "mf_write_detailed_history": detailed,
        },
        "optimization_options": {
            "nlfff_max_iterations": int(optimization["max_iterations"]),
            "nlfff_write_detailed_history": detailed,
        },
        "grad_rubin_options": {
            "gr_polarity": int(grad_rubin["polarity"]),
            "gr_max_iterations": int(grad_rubin["max_iterations"]),
            "gr_write_detailed_history": detailed,
        },
    }


def _notebook_advisories(method, preprocess, clean_alpha):
    advisories = []
    if method in ("optimization", "grad_rubin") and not preprocess:
        advisories.append(
            "shared vector preprocessing is recommended for observational "
            "Optimization/Grad-Rubin runs, but remains explicitly off"
        )
    if method == "legacy_mfr" and not preprocess:
        advisories.append(
            "shared vector preprocessing is optional for the legacy embedded-MHD "
            "MFR path and is currently off"
        )
    if method == "grad_rubin" and not clean_alpha:
        advisories.append(
            "Grad-Rubin observational alpha cleaning is recommended; "
            "alpha_source remains vector_magnetogram"
        )
    if clean_alpha and method != "grad_rubin":
        advisories.append(
            "the Grad-Rubin external-alpha product will be audited and marked "
            "unused by the selected method"
        )
    return tuple(advisories)
