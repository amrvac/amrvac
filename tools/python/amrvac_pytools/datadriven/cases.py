"""AMRVAC case staging and parameter-file generation helpers."""

from __future__ import print_function

import math
import os
from pathlib import Path
import struct
import warnings


_UNIFIED_BOUNDARY_HEADER = struct.Struct("=diidd")


def recommend_amrvac_grid(
    nx,
    ny,
    boundary_reduction_level=1,
    amrvac_refinement_level=1,
    block_sizes=(12, 14, 16, 18, 20),
):
    """Recommend trim-only AMRVAC-compatible dimensions.

    ``nx`` and ``ny`` are the selected magnetic pixels before Python boundary
    reduction. The recommended adjusted sizes are the largest values not
    exceeding the inputs for which the reduced boundary size can be represented
    by an AMRVAC base grid divisible by one of ``block_sizes`` at the requested
    refinement level.
    """

    x = _recommend_amrvac_dimension(
        nx,
        boundary_reduction_level=boundary_reduction_level,
        amrvac_refinement_level=amrvac_refinement_level,
        block_sizes=block_sizes,
    )
    y = _recommend_amrvac_dimension(
        ny,
        boundary_reduction_level=boundary_reduction_level,
        amrvac_refinement_level=amrvac_refinement_level,
        block_sizes=block_sizes,
    )
    return {
        "x": x,
        "y": y,
        "boundary_reduction_level": int(boundary_reduction_level),
        "amrvac_refinement_level": int(amrvac_refinement_level),
        "block_sizes": [int(item) for item in block_sizes],
        "total_trim": x["trim"] + y["trim"],
        "compatible": x["trim"] == 0 and y["trim"] == 0,
    }


def _validate_boundary_frame_against_metadata(boundary_filename, amrvac_meta):
    """Check an explicitly selected unified boundary frame before staging.

    Unified boundary frames contain their total horizontal dimensions and
    spacing, but not an absolute horizontal center.  The latter therefore
    remains defined by the accompanying AMRVAC metadata.
    """

    path = Path(boundary_filename).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError("cannot find boundary frame: {}".format(path))

    with path.open("rb") as handle:
        raw_header = handle.read(_UNIFIED_BOUNDARY_HEADER.size)
    if len(raw_header) != _UNIFIED_BOUNDARY_HEADER.size:
        raise ValueError(
            "boundary frame {} is too short for the unified data-driven header".format(path)
        )

    snapshot_time, nx, ny, dx_km, dy_km = _UNIFIED_BOUNDARY_HEADER.unpack(raw_header)
    if nx < 1 or ny < 1 or not all(
        math.isfinite(value) for value in (snapshot_time, dx_km, dy_km)
    ):
        raise ValueError("boundary frame {} has an invalid header".format(path))
    if dx_km <= 0.0 or dy_km <= 0.0:
        raise ValueError("boundary frame {} requires positive dx and dy".format(path))

    expected_bytes = _UNIFIED_BOUNDARY_HEADER.size + nx * ny * 3 * 8
    actual_bytes = path.stat().st_size
    if actual_bytes != expected_bytes:
        raise ValueError(
            "boundary frame {} has {} bytes; its {} x {} header requires {} bytes".format(
                path, actual_bytes, nx, ny, expected_bytes
            )
        )

    required = (
        "nx",
        "ny",
        "nx_physical",
        "ny_physical",
        "dx_cm",
        "dy_cm",
        "xprobmin1",
        "xprobmax1",
        "xprobmin2",
        "xprobmax2",
    )
    missing = [key for key in required if key not in amrvac_meta]
    if missing:
        raise ValueError(
            "cannot validate explicit boundary frame because AMRVAC metadata lacks {}".format(
                ", ".join(missing)
            )
        )

    expected_nx = int(amrvac_meta["nx"])
    expected_ny = int(amrvac_meta["ny"])
    nx_physical = int(amrvac_meta["nx_physical"])
    ny_physical = int(amrvac_meta["ny_physical"])
    expected_dx_km = float(amrvac_meta["dx_cm"]) / 1.0e5
    expected_dy_km = float(amrvac_meta["dy_cm"]) / 1.0e5
    width_x_cm = (
        float(amrvac_meta["xprobmax1"]) - float(amrvac_meta["xprobmin1"])
    ) * 1.0e9
    width_y_cm = (
        float(amrvac_meta["xprobmax2"]) - float(amrvac_meta["xprobmin2"])
    ) * 1.0e9

    differences = []
    if (nx, ny) != (expected_nx, expected_ny):
        differences.append(
            "file grid is {} x {}, metadata expects {} x {}".format(
                nx, ny, expected_nx, expected_ny
            )
        )
    if not math.isclose(dx_km, expected_dx_km, rel_tol=1.0e-9, abs_tol=1.0e-12):
        differences.append(
            "file dx is {:.16g} km, metadata expects {:.16g} km".format(
                dx_km, expected_dx_km
            )
        )
    if not math.isclose(dy_km, expected_dy_km, rel_tol=1.0e-9, abs_tol=1.0e-12):
        differences.append(
            "file dy is {:.16g} km, metadata expects {:.16g} km".format(
                dy_km, expected_dy_km
            )
        )
    if not math.isclose(
        width_x_cm, nx_physical * dx_km * 1.0e5, rel_tol=1.0e-9, abs_tol=1.0e-6
    ):
        differences.append("x1 range is inconsistent with nx_physical and file dx")
    if not math.isclose(
        width_y_cm, ny_physical * dy_km * 1.0e5, rel_tol=1.0e-9, abs_tol=1.0e-6
    ):
        differences.append("x2 range is inconsistent with ny_physical and file dy")

    if differences:
        raise ValueError(
            "boundary frame {} does not match its AMRVAC metadata: {}".format(
                path, "; ".join(differences)
            )
        )


def stage_potential_field_case(
    boundary_metadata,
    case_dir,
    amrvac_root=None,
    template_dir=None,
    boundary_filename=None,
    base_par_name="amrvac.par",
    override_par_name="data_driven_boundary.par",
    domain_nx3=None,
    block_nx1=None,
    block_nx2=None,
    block_nx3=None,
    refine_max_level=1,
    allowed_block_sizes=(12, 14, 16, 18, 20),
    potential_zshift_Mm=3.0,
    potential_field_method="fft",
    fft_padding_factor=2,
    lalpha=0.0,
    fft_top_boundary="open",
    lfff_flux_treatment="strict",
    lfff_max_flux_imbalance=0.1,
    portable_paths=True,
):
    """Stage a minimal AMRVAC potential-field case for one boundary frame.

    The staged case keeps the demo ``amrvac.par`` as a reusable base file and
    writes a small override ``.par`` containing data-specific mesh and boundary
    settings. It is intended to be run as ``amrvac -i amrvac.par
    data_driven_boundary.par``.
    """

    case_dir = _prepare_case_dir(case_dir)
    template_dir = _potential_field_template_dir(amrvac_root, template_dir)

    base_par_source = template_dir / base_par_name
    if not base_par_source.exists():
        base_par_source = template_dir / "par_init" / "potential_field.par"
    copied = _copy_case_files(case_dir, [
        template_dir / "mod_usr.t",
        template_dir / "makefile",
        (base_par_source, base_par_name),
    ])

    boundary_path, boundary_was_explicit = _resolve_boundary_path(
        boundary_metadata, boundary_filename
    )
    boundary_filename = str(boundary_path)
    boundary_filename_for_par = _path_for_case_par(boundary_path, case_dir, portable_paths)

    amrvac_meta = boundary_metadata.get("amrvac")
    if not amrvac_meta:
        raise ValueError("boundary_metadata does not contain AMRVAC mesh metadata")
    if boundary_was_explicit:
        _validate_boundary_frame_against_metadata(boundary_filename, amrvac_meta)

    override_path = case_dir / override_par_name
    override_path.write_text(
        _format_potential_field_override_par(
            amrvac_meta,
            boundary_filename=boundary_filename_for_par,
            domain_nx3=domain_nx3,
            block_nx1=block_nx1,
            block_nx2=block_nx2,
            block_nx3=block_nx3,
            refine_max_level=refine_max_level,
            allowed_block_sizes=allowed_block_sizes,
            potential_zshift_Mm=potential_zshift_Mm,
            potential_field_method=potential_field_method,
            fft_padding_factor=fft_padding_factor,
            lalpha=lalpha,
            fft_top_boundary=fft_top_boundary,
            lfff_flux_treatment=lfff_flux_treatment,
            lfff_max_flux_imbalance=lfff_max_flux_imbalance,
        ),
        encoding="utf-8",
    )

    return {
        "case_dir": str(case_dir),
        "base_par": str(case_dir / base_par_name),
        "override_par": str(override_path),
        "boundary_filename": boundary_filename,
        "boundary_filename_in_par": boundary_filename_for_par,
        "potential_field_method": str(potential_field_method).strip().lower(),
        "fft_padding_factor": int(fft_padding_factor),
        "lalpha": float(lalpha),
        "fft_top_boundary": str(fft_top_boundary).strip().lower(),
        "lfff_flux_treatment": str(lfff_flux_treatment).strip().lower(),
        "lfff_max_flux_imbalance": float(lfff_max_flux_imbalance),
        "copied": copied,
        "run_command": "mpirun -np 4 ./amrvac -i {} {}".format(base_par_name, override_par_name),
    }


def stage_magnetofrictional_relaxation_case(
    boundary_metadata,
    case_dir,
    potential_restart_file,
    amrvac_root=None,
    template_dir=None,
    boundary_filename=None,
    base_par_name="amrvac.par",
    override_par_name="data_driven_mfr.par",
    domain_nx3=None,
    block_nx1=None,
    block_nx2=None,
    block_nx3=None,
    refine_max_level=1,
    allowed_block_sizes=(12, 14, 16, 18, 20),
    base_filename="output/data_driven_mfr",
    mf_it_max=100000,
    mf_ditsave=20000,
    mf_cc=0.5,
    mf_cy=0.2,
    mf_cdivb=0.01,
    mf_log_mode="auto",
    mf_log_filename="",
    portable_paths=True,
):
    """Stage an AMRVAC magnetofrictional-relaxation case.

    The case restarts from a potential-field snapshot and uses the same V1
    three-component boundary frame to hold the bottom magnetic field fixed.
    """

    case_dir = _prepare_case_dir(case_dir)
    template_dir = _magnetofrictional_template_dir(amrvac_root, template_dir)

    copied = _copy_case_files(case_dir, [
        template_dir / "mod_usr.t",
        template_dir / "makefile",
        (template_dir / base_par_name, base_par_name),
    ])

    boundary_path, _ = _resolve_boundary_path(boundary_metadata, boundary_filename)
    restart_path = Path(potential_restart_file).expanduser().resolve()
    boundary_filename = str(boundary_path)
    potential_restart_file = str(restart_path)
    boundary_filename_for_par = _path_for_case_par(boundary_path, case_dir, portable_paths)
    potential_restart_file_for_par = _path_for_case_par(restart_path, case_dir, portable_paths)

    amrvac_meta = boundary_metadata.get("amrvac")
    if not amrvac_meta:
        raise ValueError("boundary_metadata does not contain AMRVAC mesh metadata")

    override_path = case_dir / override_par_name
    override_path.write_text(
        _format_magnetofrictional_relaxation_override_par(
            amrvac_meta,
            boundary_filename=boundary_filename_for_par,
            potential_restart_file=potential_restart_file_for_par,
            domain_nx3=domain_nx3,
            block_nx1=block_nx1,
            block_nx2=block_nx2,
            block_nx3=block_nx3,
            refine_max_level=refine_max_level,
            allowed_block_sizes=allowed_block_sizes,
            base_filename=base_filename,
            mf_it_max=mf_it_max,
            mf_ditsave=mf_ditsave,
            mf_cc=mf_cc,
            mf_cy=mf_cy,
            mf_cdivb=mf_cdivb,
            mf_log_mode=mf_log_mode,
            mf_log_filename=mf_log_filename,
        ),
        encoding="utf-8",
    )

    return {
        "case_dir": str(case_dir),
        "base_par": str(case_dir / base_par_name),
        "override_par": str(override_path),
        "boundary_filename": boundary_filename,
        "boundary_filename_in_par": boundary_filename_for_par,
        "potential_restart_file": potential_restart_file,
        "potential_restart_file_in_par": potential_restart_file_for_par,
        "base_filename": str(base_filename),
        "mf_ditsave": int(mf_ditsave),
        "mf_log_mode": str(mf_log_mode).strip().lower(),
        "mf_log_filename": str(mf_log_filename),
        "copied": copied,
        "run_command": "mpirun -np 4 ./amrvac -i {} {}".format(base_par_name, override_par_name),
    }


def stage_data_constrained_case(
    boundary_metadata,
    case_dir,
    restart_file,
    amrvac_root=None,
    template_dir=None,
    boundary_filename=None,
    base_par_name="amrvac.par",
    override_par_name="data_constrained.par",
    domain_nx3=None,
    block_nx1=None,
    block_nx2=None,
    block_nx3=None,
    refine_max_level=1,
    allowed_block_sizes=(12, 14, 16, 18, 20),
    base_filename="output/data_constrained",
    portable_paths=True,
    mhd_model="zero_beta",
    atmosphere_model=None,
    atmosphere_source="hydrostatic",
    coronal_temperature_k=1.0e6,
    rho_reference_height_cm=None,
    rho_reference_numberdensity_cm3=None,
    temperature_curve="AL-C7",
    heating_amplitude_cgs=1.0e-4,
    heating_scale_height_cm=5.0e9,
    relaxed_atmosphere_file=None,
):
    """Stage a fixed-bottom, multi-physics data-constrained MHD case.

    The restart supplies only the magnetic field.  The selected atmosphere
    rebuilds all plasma variables during the initial ``firstprocess`` pass.
    """

    import shutil

    case_dir = _prepare_case_dir(case_dir)
    template_dir = _data_constrained_template_dir(amrvac_root, template_dir)

    model_options = {"zero_beta", "isothermal", "adiabatic", "thermodynamic"}
    mhd_model = str(mhd_model).strip().lower()
    if mhd_model not in model_options:
        raise ValueError("mhd_model must be one of {}".format(sorted(model_options)))

    if atmosphere_model is None:
        atmosphere_model = {
            "zero_beta": "uniform",
            "isothermal": "corona",
            "adiabatic": "corona",
            "thermodynamic": "chromosphere",
        }[mhd_model]
    atmosphere_model = str(atmosphere_model).strip().lower()
    atmosphere_source = str(atmosphere_source).strip().lower()
    if atmosphere_model not in {"uniform", "corona", "chromosphere"}:
        raise ValueError("atmosphere_model must be uniform, corona, or chromosphere")
    if atmosphere_source not in {"hydrostatic", "relaxed_table"}:
        raise ValueError("atmosphere_source must be hydrostatic or relaxed_table")
    if mhd_model == "isothermal" and (
        atmosphere_model != "corona" or atmosphere_source != "hydrostatic"
    ):
        raise ValueError("isothermal supports only corona + hydrostatic")
    if mhd_model in {"adiabatic", "thermodynamic"} and atmosphere_model == "uniform":
        raise ValueError("uniform atmosphere is supported only by zero_beta")
    if atmosphere_source == "relaxed_table" and mhd_model == "isothermal":
        raise ValueError("isothermal does not support relaxed_table")
    if atmosphere_source == "relaxed_table" and relaxed_atmosphere_file is None:
        raise ValueError("relaxed_table requires relaxed_atmosphere_file")
    if atmosphere_source != "relaxed_table" and relaxed_atmosphere_file is not None:
        raise ValueError("relaxed_atmosphere_file requires atmosphere_source='relaxed_table'")
    if atmosphere_source != "relaxed_table" and relaxed_atmosphere_file is not None:
        raise ValueError("relaxed_atmosphere_file requires atmosphere_source='relaxed_table'")

    coronal_temperature_k = float(coronal_temperature_k)
    heating_amplitude_cgs = float(heating_amplitude_cgs)
    heating_scale_height_cm = float(heating_scale_height_cm)
    if coronal_temperature_k <= 0.0:
        raise ValueError("coronal_temperature_k must be positive")
    if heating_amplitude_cgs < 0.0:
        raise ValueError("heating_amplitude_cgs must be non-negative")
    if heating_scale_height_cm <= 0.0:
        raise ValueError("heating_scale_height_cm must be positive")
    if rho_reference_height_cm is None:
        rho_reference_height_cm = 1.0e9 if atmosphere_model == "chromosphere" else 0.0
    if rho_reference_numberdensity_cm3 is None:
        rho_reference_numberdensity_cm3 = (
            5.0e8 if atmosphere_model == "chromosphere" else 1.0e9
        )
    rho_reference_height_cm = float(rho_reference_height_cm)
    rho_reference_numberdensity_cm3 = float(rho_reference_numberdensity_cm3)
    if rho_reference_numberdensity_cm3 <= 0.0:
        raise ValueError("rho_reference_numberdensity_cm3 must be positive")

    copied = _copy_case_files(case_dir, [
        template_dir / "mod_usr.t",
        template_dir / "makefile",
        template_dir / "amrvac.h",
    ])

    selected_base_par = template_dir / ("par_" + mhd_model) / base_par_name
    if not selected_base_par.exists() and mhd_model == "zero_beta":
        # Backward compatibility for custom templates made before model folders.
        selected_base_par = template_dir / base_par_name
    if not selected_base_par.exists():
        raise FileNotFoundError("cannot find model base par: {}".format(selected_base_par))
    copied.extend(_copy_case_files(case_dir, [(selected_base_par, base_par_name)]))

    boundary_path, _ = _resolve_boundary_path(boundary_metadata, boundary_filename)
    restart_path = Path(restart_file).expanduser().resolve()
    boundary_filename = str(boundary_path)
    restart_file = str(restart_path)
    boundary_filename_for_par = _path_for_case_par(boundary_path, case_dir, portable_paths)
    restart_file_for_par = _path_for_case_par(restart_path, case_dir, portable_paths)

    relaxed_atmosphere_file_for_par = ""
    relaxed_atmosphere_file_resolved = None
    if relaxed_atmosphere_file is not None:
        relaxed_source = Path(relaxed_atmosphere_file).expanduser().resolve()
        if not relaxed_source.exists():
            raise FileNotFoundError("cannot find relaxed atmosphere table: {}".format(relaxed_source))
        atmosphere_dir = case_dir / "atmosphere"
        atmosphere_dir.mkdir(exist_ok=True)
        relaxed_target = atmosphere_dir / relaxed_source.name
        if relaxed_source != relaxed_target.resolve():
            shutil.copy2(str(relaxed_source), str(relaxed_target))
            copied.append(str(relaxed_target))
        relaxed_atmosphere_file_resolved = str(relaxed_source)
        relaxed_atmosphere_file_for_par = _path_for_case_par(
            relaxed_target, case_dir, portable_paths
        )

    amrvac_meta = boundary_metadata.get("amrvac")
    if not amrvac_meta:
        raise ValueError("boundary_metadata does not contain AMRVAC mesh metadata")

    override_path = case_dir / override_par_name
    override_path.write_text(
        _format_data_constrained_override_par(
            amrvac_meta,
            boundary_filename=boundary_filename_for_par,
            restart_file=restart_file_for_par,
            domain_nx3=domain_nx3,
            block_nx1=block_nx1,
            block_nx2=block_nx2,
            block_nx3=block_nx3,
            refine_max_level=refine_max_level,
            allowed_block_sizes=allowed_block_sizes,
            base_filename=base_filename,
            mhd_model=mhd_model,
            atmosphere_model=atmosphere_model,
            atmosphere_source=atmosphere_source,
            coronal_temperature_k=coronal_temperature_k,
            rho_reference_height_cm=rho_reference_height_cm,
            rho_reference_numberdensity_cm3=rho_reference_numberdensity_cm3,
            temperature_curve=temperature_curve,
            heating_amplitude_cgs=heating_amplitude_cgs,
            heating_scale_height_cm=heating_scale_height_cm,
            relaxed_atmosphere_file=relaxed_atmosphere_file_for_par,
        ),
        encoding="utf-8",
    )

    return {
        "case_dir": str(case_dir),
        "base_par": str(case_dir / base_par_name),
        "override_par": str(override_path),
        "boundary_filename": boundary_filename,
        "boundary_filename_in_par": boundary_filename_for_par,
        "restart_file": restart_file,
        "restart_file_in_par": restart_file_for_par,
        "base_filename": str(base_filename),
        "mhd_model": mhd_model,
        "atmosphere_model": atmosphere_model,
        "atmosphere_source": atmosphere_source,
        "selected_base_par": str(selected_base_par),
        "relaxed_atmosphere_file": relaxed_atmosphere_file_resolved,
        "relaxed_atmosphere_file_in_par": relaxed_atmosphere_file_for_par,
        "copied": copied,
        "initial_run_command": "mpirun -np 4 ./amrvac -i {} {}".format(
            base_par_name, override_par_name
        ),
        "run_command": "mpirun -np 4 ./amrvac -i {} {}".format(base_par_name, override_par_name),
    }


def stage_time_dependent_magnetofriction_case(
    sequence_metadata,
    case_dir,
    restart_file,
    amrvac_root=None,
    template_dir=None,
    base_par_name="amrvac.par",
    override_par_name="time_dependent_mf.par",
    domain_nx3=None,
    block_nx1=None,
    block_nx2=None,
    block_nx3=None,
    refine_max_level=1,
    allowed_block_sizes=(12, 14, 16, 18, 20),
    driving_time_scale=12.0,
    base_filename="output/time_dependent_mf",
    portable_paths=True,
):
    """Stage a B-only time-dependent magnetofriction case."""

    return _stage_time_dependent_case(
        sequence_metadata=sequence_metadata,
        case_dir=case_dir,
        restart_file=restart_file,
        amrvac_root=amrvac_root,
        template_dir=template_dir,
        template_name="TimeDependentMagnetofriction",
        base_par_source=None,
        base_par_name=base_par_name,
        override_par_name=override_par_name,
        domain_nx3=domain_nx3,
        block_nx1=block_nx1,
        block_nx2=block_nx2,
        block_nx3=block_nx3,
        refine_max_level=refine_max_level,
        allowed_block_sizes=allowed_block_sizes,
        driving_time_scale=driving_time_scale,
        base_filename=base_filename,
        portable_paths=portable_paths,
        firstprocess=False,
        usr_lines=[],
        result_metadata={"mode": "tmf"},
    )


def stage_data_driven_case(
    sequence_metadata,
    case_dir,
    restart_file,
    amrvac_root=None,
    template_dir=None,
    base_par_name="amrvac.par",
    override_par_name="data_driven.par",
    domain_nx3=None,
    block_nx1=None,
    block_nx2=None,
    block_nx3=None,
    refine_max_level=1,
    allowed_block_sizes=(12, 14, 16, 18, 20),
    driving_time_scale=12.0,
    base_filename="output/data_driven",
    portable_paths=True,
    mhd_model="zero_beta",
    atmosphere_model=None,
    atmosphere_source="hydrostatic",
    coronal_temperature_k=1.0e6,
    rho_reference_height_cm=None,
    rho_reference_numberdensity_cm3=None,
    bottom_numberdensity_cm3=None,
    temperature_curve="AL-C7",
    heating_amplitude_cgs=1.0e-4,
    heating_scale_height_cm=5.0e9,
    relaxed_atmosphere_file=None,
):
    """Stage a B-only time-dependent MHD case with selectable physics."""

    model_options = {"zero_beta", "isothermal", "adiabatic", "thermodynamic"}
    mhd_model = str(mhd_model).strip().lower()
    if mhd_model not in model_options:
        raise ValueError("mhd_model must be one of {}".format(sorted(model_options)))
    if atmosphere_model is None:
        atmosphere_model = {
            "zero_beta": "uniform",
            "isothermal": "corona",
            "adiabatic": "corona",
            "thermodynamic": "chromosphere",
        }[mhd_model]
    atmosphere_model = str(atmosphere_model).strip().lower()
    atmosphere_source = str(atmosphere_source).strip().lower()
    if atmosphere_model == "chromosphere":
        warnings.warn(
            "Chromospheric/transition-region atmospheres can be unstable under "
            "a directly B-driven lower boundary; a coronal atmosphere is recommended.",
            RuntimeWarning,
            stacklevel=2,
        )
    if atmosphere_model not in {"uniform", "corona", "chromosphere"}:
        raise ValueError("atmosphere_model must be uniform, corona, or chromosphere")
    if atmosphere_source not in {"hydrostatic", "relaxed_table"}:
        raise ValueError("atmosphere_source must be hydrostatic or relaxed_table")
    if mhd_model == "isothermal" and (
        atmosphere_model != "corona" or atmosphere_source != "hydrostatic"
    ):
        raise ValueError("isothermal supports only corona + hydrostatic")
    if mhd_model in {"adiabatic", "thermodynamic"} and atmosphere_model == "uniform":
        raise ValueError("uniform atmosphere is supported only by zero_beta")
    if atmosphere_source == "relaxed_table" and relaxed_atmosphere_file is None:
        raise ValueError("relaxed_table requires relaxed_atmosphere_file")
    if rho_reference_height_cm is None:
        rho_reference_height_cm = 1.0e9 if atmosphere_model == "chromosphere" else 0.0
    if rho_reference_numberdensity_cm3 is None:
        rho_reference_numberdensity_cm3 = 5.0e8 if atmosphere_model == "chromosphere" else 1.0e9
    if bottom_numberdensity_cm3 is not None:
        bottom_numberdensity_cm3 = float(bottom_numberdensity_cm3)
        if not math.isfinite(bottom_numberdensity_cm3) or bottom_numberdensity_cm3 <= 0.0:
            raise ValueError("bottom_numberdensity_cm3 must be positive and finite")

    relaxed_for_par = ""
    if relaxed_atmosphere_file is not None:
        import shutil

        relaxed_source = Path(relaxed_atmosphere_file).expanduser().resolve()
        if not relaxed_source.is_file():
            raise FileNotFoundError("cannot find relaxed atmosphere table: {}".format(relaxed_source))
        atmosphere_dir = Path(case_dir).expanduser().resolve() / "atmosphere"
        atmosphere_dir.mkdir(parents=True, exist_ok=True)
        relaxed_target = atmosphere_dir / relaxed_source.name
        if relaxed_target.resolve() != relaxed_source:
            shutil.copy2(str(relaxed_source), str(relaxed_target))
        relaxed_for_par = _path_for_case_par(relaxed_target, case_dir, portable_paths)

    root = Path(amrvac_root).expanduser().resolve() if amrvac_root is not None else Path(__file__).resolve().parents[4]
    source_root = Path(template_dir).expanduser().resolve() if template_dir is not None else root / "tests/demo4/Data_Driven/DataDriven"
    selected_base = source_root / ("par_" + mhd_model) / base_par_name
    usr_lines = [
        "  physics_model='{}'".format(mhd_model),
        "  atmosphere_model='{}'".format(atmosphere_model),
        "  atmosphere_source='{}'".format(atmosphere_source),
        "  coronal_temperature_k={}".format(_fortran_d(coronal_temperature_k)),
        "  rho_reference_height_cm={}".format(_fortran_d(rho_reference_height_cm)),
        "  rho_reference_numberdensity_cm3={}".format(_fortran_d(rho_reference_numberdensity_cm3)),
        "  temperature_curve='{}'".format(str(temperature_curve).replace("'", "''")),
        "  heating_amplitude_cgs={}".format(_fortran_d(heating_amplitude_cgs)),
        "  heating_scale_height_cm={}".format(_fortran_d(heating_scale_height_cm)),
        "  relaxed_atmosphere_file='{}'".format(relaxed_for_par.replace("'", "''")),
    ]
    if bottom_numberdensity_cm3 is not None:
        usr_lines.append(
            "  bottom_numberdensity_cm3={}".format(_fortran_d(bottom_numberdensity_cm3))
        )
    return _stage_time_dependent_case(
        sequence_metadata=sequence_metadata,
        case_dir=case_dir,
        restart_file=restart_file,
        amrvac_root=amrvac_root,
        template_dir=template_dir,
        template_name="DataDriven",
        base_par_source=selected_base,
        base_par_name=base_par_name,
        override_par_name=override_par_name,
        domain_nx3=domain_nx3,
        block_nx1=block_nx1,
        block_nx2=block_nx2,
        block_nx3=block_nx3,
        refine_max_level=refine_max_level,
        allowed_block_sizes=allowed_block_sizes,
        driving_time_scale=driving_time_scale,
        base_filename=base_filename,
        portable_paths=portable_paths,
        firstprocess=True,
        usr_lines=usr_lines,
        result_metadata={
            "mode": "data_driven",
            "mhd_model": mhd_model,
            "atmosphere_model": atmosphere_model,
            "atmosphere_source": atmosphere_source,
            "bottom_numberdensity_cm3": bottom_numberdensity_cm3,
        },
    )


def _stage_time_dependent_case(
    sequence_metadata,
    case_dir,
    restart_file,
    amrvac_root,
    template_dir,
    template_name,
    base_par_source,
    base_par_name,
    override_par_name,
    domain_nx3,
    block_nx1,
    block_nx2,
    block_nx3,
    refine_max_level,
    allowed_block_sizes,
    driving_time_scale,
    base_filename,
    portable_paths,
    firstprocess,
    usr_lines,
    result_metadata,
):
    driving_time_scale = float(driving_time_scale)
    if not math.isfinite(driving_time_scale) or driving_time_scale <= 0.0:
        raise ValueError("driving_time_scale must be positive and finite")
    sequence_dir, outputs, times = _validate_boundary_sequence_metadata(sequence_metadata)
    amrvac_meta = sequence_metadata.get("amrvac")
    if not amrvac_meta:
        raise ValueError("sequence_metadata does not contain AMRVAC mesh metadata")

    case_dir = _prepare_case_dir(case_dir)
    if template_dir is None:
        root = Path(amrvac_root).expanduser().resolve() if amrvac_root is not None else Path(__file__).resolve().parents[4]
        template_dir = root / "tests/demo4/Data_Driven" / template_name
    else:
        template_dir = Path(template_dir).expanduser().resolve()
    if base_par_source is None:
        base_par_source = template_dir / base_par_name
    copied = _copy_case_files(case_dir, [
        template_dir / "mod_usr.t",
        template_dir / "makefile",
        (base_par_source, base_par_name),
    ])
    amrvac_h = template_dir / "amrvac.h"
    if amrvac_h.exists():
        copied.extend(_copy_case_files(case_dir, [amrvac_h]))

    restart_path = Path(restart_file).expanduser().resolve()
    series_for_par = _path_for_case_par(sequence_dir, case_dir, portable_paths)
    restart_for_par = _path_for_case_par(restart_path, case_dir, portable_paths)
    mesh_lines = _time_dependent_mesh_lines(
        amrvac_meta, domain_nx3, block_nx1, block_nx2, block_nx3,
        refine_max_level, allowed_block_sizes,
    )
    reset_grid = "  reset_grid=.true.\n" if int(refine_max_level) > 1 else ""
    firstprocess_line = "  firstprocess=.true.\n" if firstprocess else ""
    override_path = case_dir / override_par_name
    override_path.write_text(
        "! Generated B-only time-dependent boundary overrides.\n\n"
        "&filelist\n"
        "  base_filename='{}'\n"
        "  restart_from_file='{}'\n"
        "{}{}"
        "/\n\n"
        "&stoplist\n  reset_time=.true.\n  reset_it=.true.\n/\n\n"
        "&meshlist\n{}\n/\n\n"
        "&usr_list\n"
        "  boundary_series_dir='{}'\n"
        "  boundary_frame_count={}\n"
        "  driving_time_scale={}\n"
        "{}\n/\n".format(
            str(base_filename).replace("'", "''"),
            restart_for_par.replace("'", "''"),
            firstprocess_line,
            reset_grid,
            "\n".join(mesh_lines),
            series_for_par.replace("'", "''"),
            len(outputs),
            _fortran_d(driving_time_scale),
            "\n".join(usr_lines),
        ),
        encoding="utf-8",
    )
    result = {
        "case_dir": str(case_dir),
        "base_par": str(case_dir / base_par_name),
        "override_par": str(override_path),
        "sequence_dir": str(sequence_dir),
        "sequence_dir_in_par": series_for_par,
        "sequence_count": len(outputs),
        "snapshot_times": times,
        "restart_file": str(restart_path),
        "restart_file_in_par": restart_for_par,
        "driving_time_scale": driving_time_scale,
        "copied": copied,
        "run_command": "mpirun -np 4 ./amrvac -i {} {}".format(base_par_name, override_par_name),
    }
    result.update(result_metadata)
    return result


def _validate_boundary_sequence_metadata(sequence_metadata):
    outputs = [Path(item).expanduser().resolve() for item in sequence_metadata.get("outputs", [])]
    if len(outputs) < 2:
        raise ValueError("time-dependent evolution requires at least two boundary frames")
    expected_names = ["B_{:04d}.dat".format(index) for index in range(1, len(outputs) + 1)]
    if [path.name for path in outputs] != expected_names:
        raise ValueError("boundary frames must be a contiguous B_0001.dat sequence")
    if len({path.parent for path in outputs}) != 1:
        raise ValueError("all boundary frames must be in one directory")
    amrvac_meta = sequence_metadata.get("amrvac") or {}
    times = []
    reference = None
    for path in outputs:
        _validate_boundary_frame_against_metadata(path, amrvac_meta)
        with path.open("rb") as handle:
            header = _UNIFIED_BOUNDARY_HEADER.unpack(handle.read(_UNIFIED_BOUNDARY_HEADER.size))
        time_value, nx, ny, dx_km, dy_km = header
        current = (nx, ny, dx_km, dy_km)
        if reference is None:
            reference = current
        elif current != reference:
            raise ValueError("boundary sequence headers are inconsistent")
        if times and time_value <= times[-1]:
            raise ValueError("boundary snapshot times must be strictly increasing")
        times.append(float(time_value))
    if abs(times[0]) > 1.0e-10:
        raise ValueError("first boundary snapshot_time must be zero")
    return outputs[0].parent, outputs, times


def _time_dependent_mesh_lines(
    amrvac_meta, domain_nx3, block_nx1, block_nx2, block_nx3,
    refine_max_level, allowed_block_sizes,
):
    refine_max_level = int(refine_max_level)
    if refine_max_level < 1:
        raise ValueError("refine_max_level must be >= 1")
    amr_factor = 2 ** (refine_max_level - 1)
    base_nx1, block_nx1 = _amrvac_base_grid_and_block(
        amrvac_meta["nx_physical"], amr_factor, block_nx1, allowed_block_sizes, "x1"
    )
    base_nx2, block_nx2 = _amrvac_base_grid_and_block(
        amrvac_meta["ny_physical"], amr_factor, block_nx2, allowed_block_sizes, "x2"
    )
    if domain_nx3 is None:
        domain_nx3 = base_nx2
    if block_nx3 is None:
        block_nx3 = block_nx2
    domain_nx3, block_nx3 = int(domain_nx3), int(block_nx3)
    if domain_nx3 % block_nx3:
        raise ValueError("domain_nx3 must be divisible by block_nx3")
    xprobmin3 = float(amrvac_meta["xprobmin3"])
    xprobmax3 = float(amrvac_meta["xprobmax3"])
    if amrvac_meta.get("boundary_plane") == "first_ghost_center":
        height = float(amrvac_meta["domain_height_10mm"])
        boundary_z = float(amrvac_meta.get("boundary_z_10mm", 0.0))
        xprobmin3 = boundary_z + 0.5 * height / float(domain_nx3 * amr_factor)
        xprobmax3 = xprobmin3 + height
    values = [
        ("block_nx1", block_nx1), ("block_nx2", block_nx2), ("block_nx3", block_nx3),
        ("refine_max_level", refine_max_level), ("domain_nx1", base_nx1),
        ("domain_nx2", base_nx2), ("domain_nx3", domain_nx3),
    ]
    lines = ["  {}={}".format(key, int(value)) for key, value in values]
    for key, value in (
        ("xprobmin1", amrvac_meta["xprobmin1"]), ("xprobmax1", amrvac_meta["xprobmax1"]),
        ("xprobmin2", amrvac_meta["xprobmin2"]), ("xprobmax2", amrvac_meta["xprobmax2"]),
        ("xprobmin3", xprobmin3), ("xprobmax3", xprobmax3),
    ):
        lines.append("  {}={}".format(key, _fortran_d(value)))
    return lines


def _prepare_case_dir(case_dir):
    from .writers import ensure_output_dir

    case_dir = ensure_output_dir(case_dir)
    (case_dir / "output").mkdir(exist_ok=True)
    return case_dir


def _copy_case_files(case_dir, sources):
    import shutil

    copied = []
    for item in sources:
        if isinstance(item, tuple):
            source, target_name = item
        else:
            source, target_name = item, Path(item).name
        source = Path(source)
        if not source.exists():
            raise FileNotFoundError("cannot find template file: {}".format(source))
        target = Path(case_dir) / target_name
        shutil.copy2(str(source), str(target))
        copied.append(str(target))
    return copied


def _resolve_boundary_path(boundary_metadata, boundary_filename):
    was_explicit = boundary_filename is not None
    if boundary_filename is None:
        outputs = boundary_metadata.get("outputs", [])
        if not outputs:
            raise ValueError("boundary_metadata does not contain output boundary frames")
        if isinstance(outputs, dict):
            boundary_filename = outputs.get("boundary_frame")
            if not boundary_filename:
                raise ValueError(
                    "boundary_metadata outputs does not contain a boundary_frame"
                )
        else:
            boundary_filename = outputs[0]
    return Path(boundary_filename).expanduser().resolve(), was_explicit


def _potential_field_template_dir(amrvac_root, template_dir):
    if template_dir is not None:
        return Path(template_dir).expanduser().resolve()
    if amrvac_root is None:
        amrvac_root = Path(__file__).resolve().parents[4]
    return Path(amrvac_root).expanduser().resolve() / "tests/demo4/Data_Driven/PotentialField"


def _magnetofrictional_template_dir(amrvac_root, template_dir):
    if template_dir is not None:
        return Path(template_dir).expanduser().resolve()
    if amrvac_root is None:
        amrvac_root = Path(__file__).resolve().parents[4]
    return Path(amrvac_root).expanduser().resolve() / "tests/demo4/Data_Driven/MagnetofrictionalRelaxation"


def _data_constrained_template_dir(amrvac_root, template_dir):
    if template_dir is not None:
        return Path(template_dir).expanduser().resolve()
    if amrvac_root is None:
        amrvac_root = Path(__file__).resolve().parents[4]
    return Path(amrvac_root).expanduser().resolve() / "tests/demo4/Data_Driven/DataConstrained"


def _path_for_case_par(path, case_dir, portable_paths):
    path = Path(path).expanduser().resolve()
    if portable_paths:
        return os.path.relpath(str(path), str(Path(case_dir).expanduser().resolve()))
    return str(path)


def _format_potential_field_override_par(
    amrvac_meta,
    boundary_filename,
    domain_nx3=None,
    block_nx1=None,
    block_nx2=None,
    block_nx3=None,
    refine_max_level=1,
    allowed_block_sizes=(12, 14, 16, 18, 20),
    potential_zshift_Mm=3.0,
    potential_field_method="fft",
    fft_padding_factor=2,
    lalpha=0.0,
    fft_top_boundary="open",
    lfff_flux_treatment="strict",
    lfff_max_flux_imbalance=0.1,
):
    refine_max_level = int(refine_max_level)
    if refine_max_level != 1:
        raise ValueError("PotentialField requires refine_max_level=1")
    potential_field_method = str(potential_field_method).strip().lower()
    if potential_field_method not in ("green", "fft"):
        raise ValueError("potential_field_method must be 'green' or 'fft'")
    fft_padding_factor = int(fft_padding_factor)
    if fft_padding_factor < 1:
        raise ValueError("fft_padding_factor must be >= 1")
    lalpha = float(lalpha)
    if not math.isfinite(lalpha):
        raise ValueError("lalpha must be finite")
    fft_top_boundary = str(fft_top_boundary).strip().lower()
    if fft_top_boundary not in ("open", "closed"):
        raise ValueError("fft_top_boundary must be 'open' or 'closed'")
    lfff_flux_treatment = str(lfff_flux_treatment).strip().lower()
    if lfff_flux_treatment not in ("strict", "subtract_mean"):
        raise ValueError(
            "lfff_flux_treatment must be 'strict' or 'subtract_mean'"
        )
    lfff_max_flux_imbalance = float(lfff_max_flux_imbalance)
    if not math.isfinite(lfff_max_flux_imbalance) or not (
        0.0 <= lfff_max_flux_imbalance <= 1.0
    ):
        raise ValueError("lfff_max_flux_imbalance must be between 0 and 1")
    amr_factor = 2 ** (refine_max_level - 1)
    allowed_block_sizes = tuple(int(item) for item in allowed_block_sizes)

    base_nx1, auto_block_nx1 = _amrvac_base_grid_and_block(
        amrvac_meta["nx_physical"],
        amr_factor,
        block_nx1,
        allowed_block_sizes,
        "x1",
    )
    base_nx2, auto_block_nx2 = _amrvac_base_grid_and_block(
        amrvac_meta["ny_physical"],
        amr_factor,
        block_nx2,
        allowed_block_sizes,
        "x2",
    )
    block_nx1 = auto_block_nx1
    block_nx2 = auto_block_nx2
    if domain_nx3 is None:
        domain_nx3 = base_nx2
    domain_nx3 = int(domain_nx3)
    if block_nx3 is None:
        block_nx3 = block_nx2
    block_nx3 = int(block_nx3)
    if domain_nx3 % block_nx3:
        raise ValueError("domain_nx3={} is not divisible by block_nx3={}".format(domain_nx3, block_nx3))

    source_plane_depth = 0.0
    if amrvac_meta.get("boundary_plane") == "first_ghost_center":
        domain_height = float(amrvac_meta["domain_height_10mm"])
        boundary_z = float(amrvac_meta.get("boundary_z_10mm", 0.0))
        if domain_height <= 0.0:
            raise ValueError("domain_height_10mm must be positive")
        dz_base = domain_height / float(domain_nx3)
        dz_finest = dz_base / float(amr_factor)
        source_plane_depth = 0.5 * dz_finest
        xprobmin3 = boundary_z + source_plane_depth
        xprobmax3 = xprobmin3 + domain_height
    else:
        xprobmin3 = float(amrvac_meta["xprobmin3"])
        xprobmax3 = float(amrvac_meta["xprobmax3"])

    if potential_field_method == "fft":
        llift_code_unit = 0.0
    else:
        # Green's source plane is below xprobmin3 by llift. Include the
        # half-finest-cell distance from the physical face to the magnetogram.
        llift_code_unit = float(potential_zshift_Mm) / 10.0 + source_plane_depth

    mesh_lines = []
    for key, value in (
        ("block_nx1", block_nx1),
        ("block_nx2", block_nx2),
        ("block_nx3", block_nx3),
    ):
        if value is not None:
            mesh_lines.append("  {}={}".format(key, int(value)))
    mesh_lines.extend(
        [
            "  refine_max_level={}".format(refine_max_level),
            "  domain_nx1={}".format(base_nx1),
            "  domain_nx2={}".format(base_nx2),
            "  domain_nx3={}".format(int(domain_nx3)),
            "  xprobmin1={}".format(_fortran_d(amrvac_meta["xprobmin1"])),
            "  xprobmax1={}".format(_fortran_d(amrvac_meta["xprobmax1"])),
            "  xprobmin2={}".format(_fortran_d(amrvac_meta["xprobmin2"])),
            "  xprobmax2={}".format(_fortran_d(amrvac_meta["xprobmax2"])),
            "  xprobmin3={}".format(_fortran_d(xprobmin3)),
            "  xprobmax3={}".format(_fortran_d(xprobmax3)),
        ]
    )
    escaped_boundary = str(boundary_filename).replace("'", "''")
    return (
        "! Data-specific overrides generated by the AMRVAC data-driven notebook.\n"
        "! Run with: amrvac -i amrvac.par data_driven_boundary.par\n"
        "\n"
        "&meshlist\n"
        "{}\n"
        "/\n"
        "\n"
        "&usr_list\n"
        "  boundary_filename='{}'\n"
        "  potential_field_method='{}'\n"
        "  fft_padding_factor={}\n"
        "  lalpha={}\n"
        "  fft_top_boundary='{}'\n"
        "  lfff_flux_treatment='{}'\n"
        "  lfff_max_flux_imbalance={}\n"
        "  llift={}\n"
        "/\n"
    ).format(
        "\n".join(mesh_lines),
        escaped_boundary,
        potential_field_method,
        fft_padding_factor,
        _fortran_d(lalpha),
        fft_top_boundary,
        lfff_flux_treatment,
        _fortran_d(lfff_max_flux_imbalance),
        _fortran_d(llift_code_unit),
    )


def _format_magnetofrictional_relaxation_override_par(
    amrvac_meta,
    boundary_filename,
    potential_restart_file,
    domain_nx3=None,
    block_nx1=None,
    block_nx2=None,
    block_nx3=None,
    refine_max_level=1,
    allowed_block_sizes=(12, 14, 16, 18, 20),
    base_filename="output/data_driven_mfr",
    mf_it_max=100000,
    mf_ditsave=20000,
    mf_cc=0.5,
    mf_cy=0.2,
    mf_cdivb=0.01,
    mf_log_mode="auto",
    mf_log_filename="",
):
    refine_max_level = int(refine_max_level)
    if refine_max_level < 1:
        raise ValueError("refine_max_level must be >= 1")
    mf_log_mode = str(mf_log_mode).strip().lower()
    if mf_log_mode not in ("auto", "append", "replace"):
        raise ValueError("mf_log_mode must be 'auto', 'append', or 'replace'")
    amr_factor = 2 ** (refine_max_level - 1)
    allowed_block_sizes = tuple(int(item) for item in allowed_block_sizes)

    base_nx1, auto_block_nx1 = _amrvac_base_grid_and_block(
        amrvac_meta["nx_physical"],
        amr_factor,
        block_nx1,
        allowed_block_sizes,
        "x1",
    )
    base_nx2, auto_block_nx2 = _amrvac_base_grid_and_block(
        amrvac_meta["ny_physical"],
        amr_factor,
        block_nx2,
        allowed_block_sizes,
        "x2",
    )
    block_nx1 = auto_block_nx1
    block_nx2 = auto_block_nx2
    if domain_nx3 is None:
        domain_nx3 = base_nx2
    domain_nx3 = int(domain_nx3)
    if block_nx3 is None:
        block_nx3 = block_nx2
    block_nx3 = int(block_nx3)
    if domain_nx3 % block_nx3:
        raise ValueError("domain_nx3={} is not divisible by block_nx3={}".format(domain_nx3, block_nx3))

    if amrvac_meta.get("boundary_plane") == "first_ghost_center":
        domain_height = float(amrvac_meta["domain_height_10mm"])
        boundary_z = float(amrvac_meta.get("boundary_z_10mm", 0.0))
        if domain_height <= 0.0:
            raise ValueError("domain_height_10mm must be positive")
        dz_base = domain_height / float(domain_nx3)
        dz_finest = dz_base / float(amr_factor)
        xprobmin3 = boundary_z + 0.5 * dz_finest
        xprobmax3 = xprobmin3 + domain_height
    else:
        xprobmin3 = float(amrvac_meta["xprobmin3"])
        xprobmax3 = float(amrvac_meta["xprobmax3"])

    mesh_lines = []
    for key, value in (
        ("block_nx1", block_nx1),
        ("block_nx2", block_nx2),
        ("block_nx3", block_nx3),
    ):
        if value is not None:
            mesh_lines.append("  {}={}".format(key, int(value)))
    mesh_lines.extend(
        [
            "  refine_max_level={}".format(refine_max_level),
            "  domain_nx1={}".format(base_nx1),
            "  domain_nx2={}".format(base_nx2),
            "  domain_nx3={}".format(int(domain_nx3)),
            "  xprobmin1={}".format(_fortran_d(amrvac_meta["xprobmin1"])),
            "  xprobmax1={}".format(_fortran_d(amrvac_meta["xprobmax1"])),
            "  xprobmin2={}".format(_fortran_d(amrvac_meta["xprobmin2"])),
            "  xprobmax2={}".format(_fortran_d(amrvac_meta["xprobmax2"])),
            "  xprobmin3={}".format(_fortran_d(xprobmin3)),
            "  xprobmax3={}".format(_fortran_d(xprobmax3)),
        ]
    )

    escaped_boundary = str(boundary_filename).replace("'", "''")
    escaped_restart = str(potential_restart_file).replace("'", "''")
    escaped_base_filename = str(base_filename).replace("'", "''")
    escaped_log_filename = str(mf_log_filename).replace("'", "''")
    reset_grid_line = "  reset_grid=.true.\n" if refine_max_level > 1 else ""
    return (
        "! Data-specific overrides generated by the AMRVAC data-driven notebook.\n"
        "! Run after the potential-field case has written the restart snapshot.\n"
        "! Run with: amrvac -i amrvac.par data_driven_mfr.par\n"
        "\n"
        "&filelist\n"
        "  base_filename='{}'\n"
        "  restart_from_file='{}'\n"
        "{}"
        "/\n"
        "\n"
        "&meshlist\n"
        "{}\n"
        "/\n"
        "\n"
        "&usr_list\n"
        "  boundary_filename='{}'\n"
        "/\n"
        "\n"
        "&mf_list\n"
        "  mf_it_max={}\n"
        "  mf_ditsave={}\n"
        "  mf_cc={}\n"
        "  mf_cy={}\n"
        "  mf_cdivb={}\n"
        "  mf_log_mode='{}'\n"
        "  mf_log_filename='{}'\n"
        "/\n"
    ).format(
        escaped_base_filename,
        escaped_restart,
        reset_grid_line,
        "\n".join(mesh_lines),
        escaped_boundary,
        int(mf_it_max),
        int(mf_ditsave),
        _fortran_d(mf_cc),
        _fortran_d(mf_cy),
        _fortran_d(mf_cdivb),
        mf_log_mode,
        escaped_log_filename,
    )


def _format_data_constrained_override_par(
    amrvac_meta,
    boundary_filename,
    restart_file,
    domain_nx3=None,
    block_nx1=None,
    block_nx2=None,
    block_nx3=None,
    refine_max_level=1,
    allowed_block_sizes=(12, 14, 16, 18, 20),
    base_filename="output/data_constrained",
    mhd_model="zero_beta",
    atmosphere_model="uniform",
    atmosphere_source="hydrostatic",
    coronal_temperature_k=1.0e6,
    rho_reference_height_cm=0.0,
    rho_reference_numberdensity_cm3=1.0e9,
    temperature_curve="AL-C7",
    heating_amplitude_cgs=1.0e-4,
    heating_scale_height_cm=5.0e9,
    relaxed_atmosphere_file="",
):
    refine_max_level = int(refine_max_level)
    if refine_max_level < 1:
        raise ValueError("refine_max_level must be >= 1")
    amr_factor = 2 ** (refine_max_level - 1)
    allowed_block_sizes = tuple(int(item) for item in allowed_block_sizes)

    base_nx1, auto_block_nx1 = _amrvac_base_grid_and_block(
        amrvac_meta["nx_physical"],
        amr_factor,
        block_nx1,
        allowed_block_sizes,
        "x1",
    )
    base_nx2, auto_block_nx2 = _amrvac_base_grid_and_block(
        amrvac_meta["ny_physical"],
        amr_factor,
        block_nx2,
        allowed_block_sizes,
        "x2",
    )
    block_nx1 = auto_block_nx1
    block_nx2 = auto_block_nx2
    if domain_nx3 is None:
        domain_nx3 = base_nx2
    domain_nx3 = int(domain_nx3)
    if block_nx3 is None:
        block_nx3 = block_nx2
    block_nx3 = int(block_nx3)
    if domain_nx3 % block_nx3:
        raise ValueError("domain_nx3={} is not divisible by block_nx3={}".format(domain_nx3, block_nx3))

    if amrvac_meta.get("boundary_plane") == "first_ghost_center":
        domain_height = float(amrvac_meta["domain_height_10mm"])
        boundary_z = float(amrvac_meta.get("boundary_z_10mm", 0.0))
        if domain_height <= 0.0:
            raise ValueError("domain_height_10mm must be positive")
        dz_base = domain_height / float(domain_nx3)
        dz_finest = dz_base / float(amr_factor)
        xprobmin3 = boundary_z + 0.5 * dz_finest
        xprobmax3 = xprobmin3 + domain_height
    else:
        xprobmin3 = float(amrvac_meta["xprobmin3"])
        xprobmax3 = float(amrvac_meta["xprobmax3"])

    mesh_lines = []
    for key, value in (
        ("block_nx1", block_nx1),
        ("block_nx2", block_nx2),
        ("block_nx3", block_nx3),
    ):
        if value is not None:
            mesh_lines.append("  {}={}".format(key, int(value)))
    mesh_lines.extend(
        [
            "  refine_max_level={}".format(refine_max_level),
            "  domain_nx1={}".format(base_nx1),
            "  domain_nx2={}".format(base_nx2),
            "  domain_nx3={}".format(int(domain_nx3)),
            "  xprobmin1={}".format(_fortran_d(amrvac_meta["xprobmin1"])),
            "  xprobmax1={}".format(_fortran_d(amrvac_meta["xprobmax1"])),
            "  xprobmin2={}".format(_fortran_d(amrvac_meta["xprobmin2"])),
            "  xprobmax2={}".format(_fortran_d(amrvac_meta["xprobmax2"])),
            "  xprobmin3={}".format(_fortran_d(xprobmin3)),
            "  xprobmax3={}".format(_fortran_d(xprobmax3)),
        ]
    )

    escaped_boundary = str(boundary_filename).replace("'", "''")
    escaped_restart = str(restart_file).replace("'", "''")
    escaped_base_filename = str(base_filename).replace("'", "''")
    escaped_model = str(mhd_model).replace("'", "''")
    escaped_atmosphere = str(atmosphere_model).replace("'", "''")
    escaped_source = str(atmosphere_source).replace("'", "''")
    escaped_curve = str(temperature_curve).replace("'", "''")
    escaped_relaxed = str(relaxed_atmosphere_file or "").replace("'", "''")
    reset_grid_line = "  reset_grid=.true.\n" if refine_max_level > 1 else ""
    return (
        "! Data-specific overrides generated by the AMRVAC data-driven notebook.\n"
        "! Run after the magnetofrictional/NLFFF relaxation has produced a restart snapshot.\n"
        "! Run with: amrvac -i amrvac.par data_constrained.par\n"
        "\n"
        "&filelist\n"
        "  base_filename='{}'\n"
        "  restart_from_file='{}'\n"
        "  firstprocess=.true.\n"
        "{}"
        "/\n"
        "\n"
        "&stoplist\n"
        "  reset_time=.true.\n"
        "  reset_it=.true.\n"
        "/\n"
        "\n"
        "&meshlist\n"
        "{}\n"
        "/\n"
        "\n"
        "&usr_list\n"
        "  boundary_filename='{}'\n"
        "  physics_model='{}'\n"
        "  atmosphere_model='{}'\n"
        "  atmosphere_source='{}'\n"
        "  coronal_temperature_k={}\n"
        "  rho_reference_height_cm={}\n"
        "  rho_reference_numberdensity_cm3={}\n"
        "  temperature_curve='{}'\n"
        "  heating_amplitude_cgs={}\n"
        "  heating_scale_height_cm={}\n"
        "  relaxed_atmosphere_file='{}'\n"
        "/\n"
    ).format(
        escaped_base_filename,
        escaped_restart,
        reset_grid_line,
        "\n".join(mesh_lines),
        escaped_boundary,
        escaped_model,
        escaped_atmosphere,
        escaped_source,
        _fortran_d(coronal_temperature_k),
        _fortran_d(rho_reference_height_cm),
        _fortran_d(rho_reference_numberdensity_cm3),
        escaped_curve,
        _fortran_d(heating_amplitude_cgs),
        _fortran_d(heating_scale_height_cm),
        escaped_relaxed,
    )


def _fortran_d(value):
    return "{:.12g}d0".format(float(value))


def _recommend_amrvac_dimension(
    size,
    boundary_reduction_level,
    amrvac_refinement_level,
    block_sizes,
):
    size = int(size)
    boundary_reducer = 2 ** (int(boundary_reduction_level) - 1)
    amr_factor = 2 ** (int(amrvac_refinement_level) - 1)
    if size < 1:
        raise ValueError("size must be positive")
    if boundary_reducer < 1 or amr_factor < 1:
        raise ValueError("levels must be >= 1")
    best = None
    for block_size in (int(item) for item in block_sizes):
        unit = boundary_reducer * amr_factor * block_size
        adjusted = (size // unit) * unit
        if adjusted < unit:
            continue
        trim = size - adjusted
        candidate = {
            "original": size,
            "adjusted": adjusted,
            "trim": trim,
            "boundary_reducer": boundary_reducer,
            "amr_factor": amr_factor,
            "boundary_physical": adjusted // boundary_reducer,
            "base": adjusted // (boundary_reducer * amr_factor),
            "block_size": block_size,
        }
        if best is None or (candidate["trim"], candidate["block_size"]) < (best["trim"], best["block_size"]):
            best = candidate
    if best is None:
        raise ValueError("size {} is too small for the requested AMRVAC block constraints".format(size))
    return best


def _amrvac_base_grid_and_block(physical_size, amr_factor, block_size, allowed_block_sizes, axis_name):
    physical_size = int(physical_size)
    if physical_size % amr_factor:
        raise ValueError(
            "{} physical size {} is not divisible by AMR factor {}; adjust the selected window or reduction level".format(
                axis_name, physical_size, amr_factor
            )
        )
    base_size = physical_size // amr_factor
    if block_size is not None:
        block_size = int(block_size)
        if base_size % block_size:
            raise ValueError("{} base grid {} is not divisible by block size {}".format(axis_name, base_size, block_size))
        return base_size, block_size
    candidates = [int(item) for item in allowed_block_sizes if base_size % int(item) == 0]
    if not candidates:
        raise ValueError(
            "{} base grid {} is not divisible by any allowed block size {}; adjust the selected window in Step 3".format(
                axis_name, base_size, list(allowed_block_sizes)
            )
        )
    return base_size, candidates[0]




__all__ = [
    "recommend_amrvac_grid",
    "stage_potential_field_case",
    "stage_magnetofrictional_relaxation_case",
    "stage_data_constrained_case",
    "stage_time_dependent_magnetofriction_case",
    "stage_data_driven_case",
]
