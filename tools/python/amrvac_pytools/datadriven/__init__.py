"""Data-driven AMRVAC preprocessing helpers.

The package is intentionally notebook-friendly: public entry points are exposed
here, while heavier optional dependencies such as SunPy are imported lazily by
the functions that need them.
"""


def make_cea_patch(*args, **kwargs):
    from .pipeline import make_cea_patch as _make_cea_patch

    return _make_cea_patch(*args, **kwargs)


def make_cea_sequence(*args, **kwargs):
    from .pipeline import make_cea_sequence as _make_cea_sequence

    return _make_cea_sequence(*args, **kwargs)


def prepare_potential_from_br(*args, **kwargs):
    from .pipeline import prepare_potential_from_br as _prepare_potential_from_br

    return _prepare_potential_from_br(*args, **kwargs)


def prepare_nlfff_from_vector(*args, **kwargs):
    from .pipeline import prepare_nlfff_from_vector as _prepare_nlfff_from_vector

    return _prepare_nlfff_from_vector(*args, **kwargs)


def prepare_data_constrained_from_vector(*args, **kwargs):
    from .pipeline import (
        prepare_data_constrained_from_vector as _prepare_data_constrained_from_vector,
    )

    return _prepare_data_constrained_from_vector(*args, **kwargs)


def prepare_tmf_sequence(*args, **kwargs):
    from .pipeline import prepare_tmf_sequence as _prepare_tmf_sequence

    return _prepare_tmf_sequence(*args, **kwargs)


def prepare_mhd_sequence(*args, **kwargs):
    from .pipeline import prepare_mhd_sequence as _prepare_mhd_sequence

    return _prepare_mhd_sequence(*args, **kwargs)


def prepare_boundary_frame(*args, **kwargs):
    from .pipeline import prepare_boundary_frame as _prepare_boundary_frame

    return _prepare_boundary_frame(*args, **kwargs)


def prepare_boundary_sequence(*args, **kwargs):
    from .pipeline import prepare_boundary_sequence as _prepare_boundary_sequence

    return _prepare_boundary_sequence(*args, **kwargs)


def recommend_amrvac_grid(*args, **kwargs):
    from .cases import recommend_amrvac_grid as _recommend_amrvac_grid

    return _recommend_amrvac_grid(*args, **kwargs)


def stage_potential_field_case(*args, **kwargs):
    from .cases import stage_potential_field_case as _stage_potential_field_case

    return _stage_potential_field_case(*args, **kwargs)


def stage_magnetofrictional_relaxation_case(*args, **kwargs):
    from .cases import (
        stage_magnetofrictional_relaxation_case as _stage_magnetofrictional_relaxation_case,
    )

    return _stage_magnetofrictional_relaxation_case(*args, **kwargs)


def stage_data_constrained_case(*args, **kwargs):
    from .cases import stage_data_constrained_case as _stage_data_constrained_case

    return _stage_data_constrained_case(*args, **kwargs)


def stage_time_dependent_magnetofriction_case(*args, **kwargs):
    from .cases import (
        stage_time_dependent_magnetofriction_case as _stage_time_dependent_magnetofriction_case,
    )

    return _stage_time_dependent_magnetofriction_case(*args, **kwargs)


def stage_data_driven_case(*args, **kwargs):
    from .cases import stage_data_driven_case as _stage_data_driven_case

    return _stage_data_driven_case(*args, **kwargs)


def read_relaxation_diagnostics(*args, **kwargs):
    from .diagnostics import read_relaxation_diagnostics as _read_relaxation_diagnostics

    return _read_relaxation_diagnostics(*args, **kwargs)


def plot_relaxation_diagnostics(*args, **kwargs):
    from .diagnostics import plot_relaxation_diagnostics as _plot_relaxation_diagnostics

    return _plot_relaxation_diagnostics(*args, **kwargs)


def select_relaxation_restart(*args, **kwargs):
    from .diagnostics import select_relaxation_restart as _select_relaxation_restart

    return _select_relaxation_restart(*args, **kwargs)


def find_relaxation_restart_snapshots(*args, **kwargs):
    from .diagnostics import find_relaxation_restart_snapshots as _find_relaxation_restart_snapshots

    return _find_relaxation_restart_snapshots(*args, **kwargs)


def inspect_magnetic_input(*args, **kwargs):
    from .fits_io import inspect_magnetic_input as _inspect_magnetic_input

    return _inspect_magnetic_input(*args, **kwargs)


def discover_magnetic_inputs(*args, **kwargs):
    from .fits_io import discover_magnetic_inputs as _discover_magnetic_inputs

    return _discover_magnetic_inputs(*args, **kwargs)


def load_magnetic_preview(*args, **kwargs):
    from .fits_io import load_magnetic_preview as _load_magnetic_preview

    return _load_magnetic_preview(*args, **kwargs)


def load_cea_patch_preview(*args, **kwargs):
    from .fits_io import load_cea_patch_preview as _load_cea_patch_preview

    return _load_cea_patch_preview(*args, **kwargs)


def summarize_magnetic_input(*args, **kwargs):
    from .fits_io import summarize_magnetic_input as _summarize_magnetic_input

    return _summarize_magnetic_input(*args, **kwargs)


def create_data_constrain_workflow(*args, **kwargs):
    from .workflow import create_data_constrain_workflow as _create_workflow

    return _create_workflow(*args, **kwargs)


def create_data_driven_workflow(*args, **kwargs):
    from .data_driven_workflow import create_data_driven_workflow as _create_workflow

    return _create_workflow(*args, **kwargs)


def check_data_driven_dependencies(*args, **kwargs):
    from .workflow import check_data_driven_dependencies as _check_dependencies

    return _check_dependencies(*args, **kwargs)


def stage_grid_config(*args, **kwargs):
    from .workflow import StageGridConfig

    return StageGridConfig(*args, **kwargs)


def plan_data_constrain_region(*args, **kwargs):
    from .workflow import plan_data_constrain_region as _plan_region

    return _plan_region(*args, **kwargs)


def plot_data_constrain_region(*args, **kwargs):
    from .workflow import plot_data_constrain_region as _plot_region

    return _plot_region(*args, **kwargs)


def format_case_commands(*args, **kwargs):
    from .workflow import format_case_commands as _format_case_commands

    return _format_case_commands(*args, **kwargs)


__all__ = [
    "discover_magnetic_inputs",
    "inspect_magnetic_input",
    "load_magnetic_preview",
    "load_cea_patch_preview",
    "summarize_magnetic_input",
    "make_cea_patch",
    "make_cea_sequence",
    "prepare_potential_from_br",
    "prepare_nlfff_from_vector",
    "prepare_data_constrained_from_vector",
    "prepare_tmf_sequence",
    "prepare_mhd_sequence",
    "prepare_boundary_frame",
    "prepare_boundary_sequence",
    "recommend_amrvac_grid",
    "stage_potential_field_case",
    "stage_magnetofrictional_relaxation_case",
    "stage_data_constrained_case",
    "stage_time_dependent_magnetofriction_case",
    "stage_data_driven_case",
    "read_relaxation_diagnostics",
    "plot_relaxation_diagnostics",
    "select_relaxation_restart",
    "find_relaxation_restart_snapshots",
    "create_data_constrain_workflow",
    "create_data_driven_workflow",
    "check_data_driven_dependencies",
    "stage_grid_config",
    "plan_data_constrain_region",
    "plot_data_constrain_region",
    "format_case_commands",
]
