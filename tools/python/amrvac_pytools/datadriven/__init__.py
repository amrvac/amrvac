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


def clean_grad_rubin_alpha(*args, **kwargs):
    from .alpha_cleaning import clean_grad_rubin_alpha as _clean_grad_rubin_alpha

    return _clean_grad_rubin_alpha(*args, **kwargs)


def write_alpha_cleaning_audit(*args, **kwargs):
    from .alpha_cleaning import write_alpha_cleaning_audit as _write_audit

    return _write_audit(*args, **kwargs)


def alpha_cleaning_config(*args, **kwargs):
    from .alpha_cleaning import AlphaCleaningConfig as _AlphaCleaningConfig

    return _AlphaCleaningConfig(*args, **kwargs)


def AlphaCleaningConfig(*args, **kwargs):
    """Lazy public constructor retaining the class-style notebook spelling."""

    return alpha_cleaning_config(*args, **kwargs)


def preprocess_vector_magnetogram(*args, **kwargs):
    from .preprocessing import preprocess_vector_magnetogram as _preprocess_vector_magnetogram

    return _preprocess_vector_magnetogram(*args, **kwargs)


def VectorPreprocessingConfig(*args, **kwargs):
    from .preprocessing import VectorPreprocessingConfig as _VectorPreprocessingConfig

    return _VectorPreprocessingConfig(*args, **kwargs)


def recommend_vector_preprocessing(*args, **kwargs):
    from .preprocessing import recommend_vector_preprocessing as _recommend_vector_preprocessing

    return _recommend_vector_preprocessing(*args, **kwargs)


def write_external_alpha(*args, **kwargs):
    from .external_alpha import write_external_alpha as _write_external_alpha

    return _write_external_alpha(*args, **kwargs)


def read_external_alpha(*args, **kwargs):
    from .external_alpha import read_external_alpha as _read_external_alpha

    return _read_external_alpha(*args, **kwargs)


def validate_external_alpha_grid(*args, **kwargs):
    from .external_alpha import validate_external_alpha_grid as _validate_external_alpha_grid

    return _validate_external_alpha_grid(*args, **kwargs)


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


def stage_optimization_nlfff_case(*args, **kwargs):
    from .cases import stage_optimization_nlfff_case as _stage_optimization_nlfff_case

    return _stage_optimization_nlfff_case(*args, **kwargs)


def stage_grad_rubin_nlfff_case(*args, **kwargs):
    from .cases import stage_grad_rubin_nlfff_case as _stage_grad_rubin_nlfff_case

    return _stage_grad_rubin_nlfff_case(*args, **kwargs)


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


def read_nlfff_metrics(*args, **kwargs):
    from .diagnostics import read_nlfff_metrics as _read_nlfff_metrics

    return _read_nlfff_metrics(*args, **kwargs)


def plot_relaxation_diagnostics(*args, **kwargs):
    from .diagnostics import plot_relaxation_diagnostics as _plot_relaxation_diagnostics

    return _plot_relaxation_diagnostics(*args, **kwargs)


def plot_nlfff_metrics(*args, **kwargs):
    from .diagnostics import plot_nlfff_metrics as _plot_nlfff_metrics

    return _plot_nlfff_metrics(*args, **kwargs)


def select_relaxation_restart(*args, **kwargs):
    from .diagnostics import select_relaxation_restart as _select_relaxation_restart

    return _select_relaxation_restart(*args, **kwargs)


def find_relaxation_restart_snapshots(*args, **kwargs):
    from .diagnostics import find_relaxation_restart_snapshots as _find_relaxation_restart_snapshots

    return _find_relaxation_restart_snapshots(*args, **kwargs)


def nlfff_metrics_csv_path(*args, **kwargs):
    from .diagnostics import nlfff_metrics_csv_path as _nlfff_metrics_csv_path

    return _nlfff_metrics_csv_path(*args, **kwargs)


def normalize_nlfff_metrics(*args, **kwargs):
    from .diagnostics import normalize_nlfff_metrics as _normalize_nlfff_metrics

    return _normalize_nlfff_metrics(*args, **kwargs)


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


def normalize_notebook_options(*args, **kwargs):
    from .notebook import normalize_notebook_options as _normalize_options

    return _normalize_options(*args, **kwargs)


def notebook_configuration_report(*args, **kwargs):
    from .notebook import notebook_configuration_report as _configuration_report

    return _configuration_report(*args, **kwargs)


def create_notebook_workflow(*args, **kwargs):
    from .notebook import create_notebook_workflow as _create_workflow

    return _create_workflow(*args, **kwargs)


def default_notebook_nlfff_run_options(*args, **kwargs):
    from .notebook import default_notebook_nlfff_run_options as _default_options

    return _default_options(*args, **kwargs)


def notebook_boundary_options(*args, **kwargs):
    from .notebook import notebook_boundary_options as _boundary_options

    return _boundary_options(*args, **kwargs)


def notebook_single_frame_options(*args, **kwargs):
    from .notebook import notebook_single_frame_options as _frame_options

    return _frame_options(*args, **kwargs)


def notebook_reference_frame_options(*args, **kwargs):
    from .notebook import notebook_reference_frame_options as _reference_options

    return _reference_options(*args, **kwargs)


def notebook_sequence_options(*args, **kwargs):
    from .notebook import notebook_sequence_options as _sequence_options

    return _sequence_options(*args, **kwargs)


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
    "clean_grad_rubin_alpha",
    "write_alpha_cleaning_audit",
    "alpha_cleaning_config",
    "AlphaCleaningConfig",
    "preprocess_vector_magnetogram",
    "VectorPreprocessingConfig",
    "recommend_vector_preprocessing",
    "write_external_alpha",
    "read_external_alpha",
    "validate_external_alpha_grid",
    "prepare_tmf_sequence",
    "prepare_mhd_sequence",
    "prepare_boundary_frame",
    "prepare_boundary_sequence",
    "recommend_amrvac_grid",
    "stage_potential_field_case",
    "stage_magnetofrictional_relaxation_case",
    "stage_optimization_nlfff_case",
    "stage_grad_rubin_nlfff_case",
    "stage_data_constrained_case",
    "stage_time_dependent_magnetofriction_case",
    "stage_data_driven_case",
    "read_relaxation_diagnostics",
    "read_nlfff_metrics",
    "plot_relaxation_diagnostics",
    "plot_nlfff_metrics",
    "select_relaxation_restart",
    "find_relaxation_restart_snapshots",
    "nlfff_metrics_csv_path",
    "normalize_nlfff_metrics",
    "create_data_constrain_workflow",
    "create_data_driven_workflow",
    "check_data_driven_dependencies",
    "normalize_notebook_options",
    "notebook_configuration_report",
    "create_notebook_workflow",
    "default_notebook_nlfff_run_options",
    "notebook_boundary_options",
    "notebook_single_frame_options",
    "notebook_reference_frame_options",
    "notebook_sequence_options",
    "stage_grid_config",
    "plan_data_constrain_region",
    "plot_data_constrain_region",
    "format_case_commands",
]
