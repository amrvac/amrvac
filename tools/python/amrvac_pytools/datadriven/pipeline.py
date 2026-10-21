"""High-level AMRVAC data-driven preprocessing API."""

from __future__ import print_function

from datetime import datetime
from pathlib import Path


SOLAR_RADIUS_CM = 6.955e10


def make_cea_patch(
    field,
    inclination,
    azimuth,
    disambig,
    output_dir,
    center_lon,
    center_lat,
    width_degree,
    height_degree,
    resolution_degree=0.03,
    sampling_mode="sharp",
    oversample_resolution_degree=0.01,
    smooth_sigma_degree=0.01,
    smooth_truncate=2.0,
    interpolation_order=3,
    disambig_method=2,
    fill_value=None,
    disk_lon=0.0,
    disk_latitude=None,
    p_angle=None,
    components_only=False,
    quicklook=True,
    vmax=1500.0,
    output_prefix=None,
):
    """Make SHARP-like CEA Br/Bt/Bp files from raw HMI vector segments."""

    try:
        import astropy.units as u
        import numpy as np
        from astropy.io import fits
        from sunpy.coordinates import sun
        from sunpy.map import Map
    except ImportError as error:
        raise ImportError("astropy, numpy, and sunpy are required for raw HMI CEA remapping") from error

    from .cea import (
        gaussian_block_reduce,
        hmi_native_vector_components,
        make_cea_patch_grid,
        native_to_heliographic_components,
    )
    from .cartesian import local_cartesian_from_heliographic
    from .hmi import hmi_disambiguate_azimuth
    from .writers import ensure_output_dir, write_json

    output_dir = ensure_output_dir(output_dir)
    fill_value = np.nan if fill_value is None else fill_value
    prefix = _clean_output_prefix(output_prefix)

    field_map = Map(field)
    inclination_map = Map(inclination)
    azimuth_map = Map(azimuth)
    disambig_map = Map(disambig)

    corrected_azimuth = hmi_disambiguate_azimuth(
        azimuth_map.data,
        disambig_map.data,
        method=disambig_method,
    )
    bxi, beta, bzeta = hmi_native_vector_components(
        field_map.data, inclination_map.data, corrected_azimuth
    )
    grid = make_cea_patch_grid(
        center_lon, center_lat, width_degree, height_degree, resolution_degree
    )
    sampled_bxi, sampled_beta, sampled_bzeta, sampling_metadata = _sample_native_components(
        field_map,
        bxi,
        beta,
        bzeta,
        grid,
        sampling_mode=sampling_mode,
        interpolation_order=interpolation_order,
        fill_value=fill_value,
        oversample_resolution_degree=oversample_resolution_degree,
        smooth_sigma_degree=smooth_sigma_degree,
        smooth_truncate=smooth_truncate,
        gaussian_block_reduce=gaussian_block_reduce,
    )

    disk_lat = disk_latitude
    if disk_lat is None:
        disk_lat = float(field_map.observer_coordinate.lat.to_value(u.deg))
    p_angle_value = p_angle
    if p_angle_value is None:
        try:
            p_angle_value = float(-sun.P(field_map.date).to_value(u.deg))
        except Exception:
            p_angle_value = float(-field_map.meta.get("crota2", field_map.meta.get("CROTA2", 0.0)))

    br, bt, bp = native_to_heliographic_components(
        sampled_bxi,
        sampled_beta,
        sampled_bzeta,
        grid.lon,
        grid.lat,
        disk_lon_degree=disk_lon,
        disk_lat_degree=disk_lat,
        p_angle_degree=p_angle_value,
    )
    bx, by, bz = local_cartesian_from_heliographic(br, bt, bp)

    component_paths = {
        "Br": _prefixed_output_path(output_dir, prefix, "Br.fits"),
        "Bt": _prefixed_output_path(output_dir, prefix, "Bt.fits"),
        "Bp": _prefixed_output_path(output_dir, prefix, "Bp.fits"),
    }
    _write_component_fits(component_paths["Br"], br, field_map, grid, "Br", disk_lon, disk_lat, p_angle_value)
    _write_component_fits(component_paths["Bt"], bt, field_map, grid, "Bt", disk_lon, disk_lat, p_angle_value)
    _write_component_fits(component_paths["Bp"], bp, field_map, grid, "Bp", disk_lon, disk_lat, p_angle_value)
    if not components_only:
        sampled_disambig = _sample_nearest_segment(disambig_map, grid, fill_value)
        _write_component_fits(
            _prefixed_output_path(output_dir, prefix, "disambig.fits"),
            sampled_disambig,
            field_map,
            grid,
            "disambig",
            disk_lon,
            disk_lat,
            p_angle_value,
        )
        np.savez_compressed(
            _prefixed_output_path(output_dir, prefix, "cea_patch.npz"),
            x=grid.x, y=grid.y, lon=grid.lon, lat=grid.lat,
            bxi=sampled_bxi, beta=sampled_beta, bzeta=sampled_bzeta,
            br=br, bt=bt, bp=bp, bx=bx, by=by, bz=bz, disambig=sampled_disambig,
        )
        np.savez_compressed(
            _prefixed_output_path(output_dir, prefix, "cartesian_boundary.npz"),
            bx=bx,
            by=by,
            bz=bz,
        )

    metadata = {
        "inputs": {
            "field": str(field),
            "inclination": str(inclination),
            "azimuth": str(azimuth),
            "disambig": str(disambig),
        },
        "grid": {
            "center_lon": center_lon,
            "center_lat": center_lat,
            "width_degree": width_degree,
            "height_degree": height_degree,
            "resolution_degree": resolution_degree,
            "shape": list(br.shape),
        },
        "geometry": {
            "disk_lon": disk_lon,
            "disk_lat": disk_lat,
            "p_angle": p_angle_value,
            "longitude_frame": "HeliographicStonyhurst",
        },
        "sampling": sampling_metadata,
        "cartesian_boundary": {"Bx": "Bp", "By": "-Bt", "Bz": "Br"},
        "outputs": {
            "Br": str(component_paths["Br"]),
            "Bt": str(component_paths["Bt"]),
            "Bp": str(component_paths["Bp"]),
        },
    }
    write_json(_prefixed_output_path(output_dir, prefix, "parameters.json"), metadata)
    if quicklook:
        from .writers import write_quicklook

        write_quicklook(_prefixed_output_path(output_dir, prefix, "quicklook.png"), bx, by, bz, vmax=vmax)
    return metadata


def make_cea_sequence(
    raw_hmi_input,
    output_dir,
    center_lon,
    center_lat,
    width_degree,
    height_degree,
    resolution_degree=0.03,
    sampling_mode="sharp",
    oversample_resolution_degree=0.01,
    smooth_sigma_degree=0.01,
    smooth_truncate=2.0,
    interpolation_order=3,
    disambig_method=2,
    fill_value=None,
    disk_lon=0.0,
    disk_latitude=None,
    p_angle=None,
    components_only=True,
    quicklook_first=True,
    vmax=1500.0,
    progress=False,
    resume=False,
):
    """Make a flat CEA-like Br/Bt/Bp sequence from raw HMI vector segments."""

    from .writers import ensure_output_dir, write_json

    output_dir = ensure_output_dir(output_dir)
    frames = _raw_hmi_frames_from_input(raw_hmi_input)
    if not frames:
        raise ValueError("no complete raw HMI frames found")

    converted = []
    for index, frame in enumerate(_progress_iterator(frames, progress, "Remapping raw HMI frames"), start=1):
        prefix = "raw_hmi_{:04d}".format(index)
        expected = {
            component: _prefixed_output_path(output_dir, prefix, "{}.fits".format(component))
            for component in ("Br", "Bt", "Bp")
        }
        if resume and all(path.exists() for path in expected.values()):
            metadata = {
                "inputs": {name: str(frame[name]) for name in ("field", "inclination", "azimuth", "disambig")},
                "outputs": {name: str(path) for name, path in expected.items()},
            }
        else:
            metadata = make_cea_patch(
                field=frame["field"],
                inclination=frame["inclination"],
                azimuth=frame["azimuth"],
                disambig=frame["disambig"],
                output_dir=output_dir,
                center_lon=center_lon,
                center_lat=center_lat,
                width_degree=width_degree,
                height_degree=height_degree,
                resolution_degree=resolution_degree,
                sampling_mode=sampling_mode,
                oversample_resolution_degree=oversample_resolution_degree,
                smooth_sigma_degree=smooth_sigma_degree,
                smooth_truncate=smooth_truncate,
                interpolation_order=interpolation_order,
                disambig_method=disambig_method,
                fill_value=fill_value,
                disk_lon=disk_lon,
                disk_latitude=disk_latitude,
                p_angle=p_angle,
                components_only=components_only,
                quicklook=bool(quicklook_first and index == 1),
                vmax=vmax,
                output_prefix=prefix,
            )
        converted.append({
            "index": index,
            "key": frame.get("key", prefix),
            "prefix": prefix,
            "inputs": metadata["inputs"],
            "outputs": metadata["outputs"],
        })

    sequence_metadata = {
        "mode": "raw_hmi_to_cea_sequence",
        "output_dir": str(output_dir),
        "frame_count": len(converted),
        "frames": converted,
        "outputs": {
            "component_pattern": "raw_hmi_####.[Br|Bt|Bp].fits",
        },
    }
    write_json(output_dir / "sequence_parameters.json", sequence_metadata)
    return sequence_metadata


def make_fixed_cea_sequence(
    input_dir,
    output_dir,
    cea_patch,
    indices=None,
    progress=False,
    resume=False,
):
    """Resample a prepared Br/Bt/Bp sequence onto one fixed CEA grid."""

    try:
        import numpy as np
        from sunpy.map import Map
    except ImportError as error:
        raise ImportError(
            "numpy and sunpy are required for fixed-CEA SHARP resampling"
        ) from error

    from .cea import (
        _sample_at_source_pixels,
        _source_pixels_from_heliographic,
        make_cea_patch_grid,
    )
    from .fits_io import discover_vector_sequence
    from .writers import ensure_output_dir, write_json

    series = discover_vector_sequence(input_dir)
    frame_count = len(series["br"])
    selected = list(range(frame_count)) if indices is None else [int(value) for value in indices]
    if not selected:
        raise ValueError("selected fixed-CEA sequence is empty")
    if selected != sorted(set(selected)):
        raise ValueError("fixed-CEA sequence indices must be unique and increasing")
    if selected[0] < 0 or selected[-1] >= frame_count:
        raise IndexError("fixed-CEA sequence index outside 0..{}".format(frame_count - 1))

    patch = dict(cea_patch)
    grid = make_cea_patch_grid(
        patch["center_lon"], patch["center_lat"],
        patch["width_degree"], patch["height_degree"],
        patch.get("resolution_degree", 0.03),
    )
    output_dir = ensure_output_dir(output_dir)
    converted = []
    iterator = _progress_iterator(selected, progress, "Resampling SHARP frames to fixed CEA")
    for output_index, source_index in enumerate(iterator, start=1):
        prefix = "fixed_cea_{:04d}".format(output_index)
        targets = {
            name: _prefixed_output_path(output_dir, prefix, "{}.fits".format(name))
            for name in ("Br", "Bt", "Bp")
        }
        if not (resume and all(path.exists() for path in targets.values())):
            maps = {
                "Br": Map(str(series["br"][source_index])),
                "Bt": Map(str(series["bt"][source_index])),
                "Bp": Map(str(series["bp"][source_index])),
            }
            for component, source_map in maps.items():
                pixel_yx = _source_pixels_from_heliographic(source_map, grid.lon, grid.lat)
                sampled = _sample_at_source_pixels(
                    source_map.data, pixel_yx, grid.shape, order=1, cval=np.nan
                )
                meta = source_map.meta
                disk_lon = float(meta.get("HGLN_OBS", 0.0))
                disk_lat = float(meta.get("HGLT_OBS", 0.0))
                p_angle = float(meta.get("P_ANGLE", meta.get("CROTA2", 0.0)))
                _write_component_fits(
                    targets[component], sampled, source_map, grid, component,
                    disk_lon, disk_lat, p_angle,
                )
        converted.append({
            "output_index": output_index,
            "source_index": source_index,
            "outputs": {name: str(path) for name, path in targets.items()},
        })
    metadata = {
        "mode": "fixed_cea_vector_sequence",
        "input_dir": str(Path(input_dir).resolve()),
        "output_dir": str(output_dir),
        "frame_count": len(converted),
        "source_indices": selected,
        "cea_patch": patch,
        "frames": converted,
    }
    write_json(output_dir / "sequence_parameters.json", metadata)
    return metadata


def prepare_potential_from_br(
    br,
    output_dir,
    window=None,
    level=1,
    geometry=None,
    nghost=2,
    quicklook=True,
    vmax=500.0,
):
    """Prepare scenario 1: single Br map for an initial potential field."""

    import numpy as np

    from .cartesian import (
        boundary_metadata,
        crop_components_with_padding,
        format_guo_amrvac_parameters,
        multigrid_reduce_components,
        multigrid_reducer,
        write_potential_boundary,
    )
    from .fits_io import load_br
    from .writers import ensure_output_dir, write_json, write_quicklook

    output_dir = ensure_output_dir(output_dir)
    bz, header = _array_or_br_fits(br, load_br)
    bx = np.zeros_like(bz)
    by = np.zeros_like(bz)
    reducer = multigrid_reducer(level)
    crop = _window_for_shape(window, bz.shape)
    bx_crop, by_crop, bz_crop, padding = crop_components_with_padding(
        bx, by, bz,
        x0=crop["x0"], y0=crop["y0"], nx=crop["nx"], ny=crop["ny"],
        pad_x=nghost * reducer, pad_y=nghost * reducer,
    )
    bx_mg, by_mg, bz_mg = multigrid_reduce_components(
        bx_crop, by_crop, bz_crop, level=level
    )
    meta = _metadata_for_boundary(
        header, geometry, crop, level, bx_mg.shape, nghost, boundary_metadata
    )
    write_potential_boundary(output_dir / "potential_boundary.dat", bz_mg, meta)
    parameters_text = format_guo_amrvac_parameters(meta)
    (output_dir / "amrvac_parameters.txt").write_text(parameters_text, encoding="utf-8")
    if quicklook:
        write_quicklook(output_dir / "potential_quicklook.png", bx_mg, by_mg, bz_mg, vmax=vmax)
    metadata = _base_metadata(
        "potential",
        crop,
        level,
        padding,
        meta,
        {"potential_boundary": str(output_dir / "potential_boundary.dat")},
    )
    write_json(output_dir / "parameters.json", metadata)
    return metadata


def prepare_nlfff_from_vector(
    br,
    bt,
    bp,
    output_dir,
    window=None,
    level=1,
    geometry=None,
    preprocess=False,
    preprocess_mu3=0.1,
    preprocess_mu4=0.1,
    preprocess_max_iter=5000,
    preprocess_tol=1.0e-4,
    nghost=2,
    quicklook=True,
    vmax=500.0,
    preprocessing_mode=None,
):
    """Prepare scenario 2: single Br/Bt/Bp map for NLFFF relaxation."""

    return _prepare_vector_static(
        "nlfff",
        br, bt, bp, output_dir, window, level, geometry,
        preprocess, preprocess_mu3, preprocess_mu4, preprocess_max_iter,
        preprocess_tol, nghost, quicklook, vmax,
        preprocessing_mode=preprocessing_mode,
    )


def prepare_data_constrained_from_vector(
    br,
    bt,
    bp,
    output_dir,
    window=None,
    level=1,
    geometry=None,
    preprocess=False,
    preprocess_mu3=0.1,
    preprocess_mu4=0.1,
    preprocess_max_iter=5000,
    preprocess_tol=1.0e-4,
    nghost=2,
    quicklook=True,
    vmax=500.0,
    preprocessing_mode=None,
):
    """Prepare scenario 3: fixed-bottom data-constrained MHD inputs."""

    return _prepare_vector_static(
        "data_constrained",
        br, bt, bp, output_dir, window, level, geometry,
        preprocess, preprocess_mu3, preprocess_mu4, preprocess_max_iter,
        preprocess_tol, nghost, quicklook, vmax,
        preprocessing_mode=preprocessing_mode,
        extra={"bottom_velocity": "zero", "time_dependence": "fixed_boundary"},
    )


def prepare_tmf_sequence(
    input_dir,
    output_dir,
    window=None,
    level=1,
    geometry=None,
    preprocess=False,
    velocity=None,
    quicklook=False,
):
    """Prepare scenario 4: B-only time-dependent magnetofriction sequence."""

    if velocity is not None:
        raise NotImplementedError("velocity input is reserved for a future DAVE adapter")
    return _prepare_vector_sequence(
        "tmf_sequence", input_dir, output_dir, window, level, geometry,
        preprocess, quicklook,
    )


def prepare_mhd_sequence(
    input_dir,
    output_dir,
    window=None,
    level=1,
    geometry=None,
    preprocess=False,
    velocity=None,
    mhd_mode="full_mhd",
    quicklook=False,
):
    """Prepare scenario 5: B-only full/isothermal MHD sequence."""

    if velocity is not None:
        raise NotImplementedError("velocity input is reserved for a future DAVE adapter")
    if mhd_mode not in ("full_mhd", "isothermal_mhd"):
        raise ValueError("mhd_mode must be 'full_mhd' or 'isothermal_mhd'")
    return _prepare_vector_sequence(
        mhd_mode, input_dir, output_dir, window, level, geometry,
        preprocess, quicklook,
    )


def prepare_boundary_frame(
    input_dir,
    output_dir,
    window=None,
    level=1,
    geometry=None,
    snapshot_index=0,
    preprocess=False,
    preprocess_mu3=0.1,
    preprocess_mu4=0.1,
    preprocess_max_iter=5000,
    preprocess_tol=1.0e-4,
    preprocessing_mode=None,
    fail_on_nonconvergence=False,
    nghost=2,
    quicklook=True,
    vmax=500.0,
):
    """Write one unified three-component magnetic-boundary frame."""

    metadata = _prepare_unified_boundary(
        input_dir=input_dir,
        output_dir=output_dir,
        window=window,
        level=level,
        geometry=geometry,
        snapshot_index=snapshot_index,
        all_frames=False,
        preprocess=preprocess,
        preprocess_mu3=preprocess_mu3,
        preprocess_mu4=preprocess_mu4,
        preprocess_max_iter=preprocess_max_iter,
        preprocess_tol=preprocess_tol,
        preprocessing_mode=preprocessing_mode,
        fail_on_nonconvergence=fail_on_nonconvergence,
        nghost=nghost,
        quicklook=quicklook,
        vmax=vmax,
    )
    return metadata


def prepare_boundary_sequence(
    input_dir,
    output_dir,
    window=None,
    level=1,
    geometry=None,
    preprocess=False,
    preprocess_mu3=0.1,
    preprocess_mu4=0.1,
    preprocess_max_iter=5000,
    preprocess_tol=1.0e-4,
    preprocessing_mode=None,
    fail_on_nonconvergence=False,
    nghost=2,
    quicklook=True,
    vmax=500.0,
    indices=None,
    require_timestamps=False,
    progress=False,
):
    """Write all frames as unified three-component magnetic-boundary files."""

    metadata = _prepare_unified_boundary(
        input_dir=input_dir,
        output_dir=output_dir,
        window=window,
        level=level,
        geometry=geometry,
        snapshot_index=0,
        all_frames=True,
        preprocess=preprocess,
        preprocess_mu3=preprocess_mu3,
        preprocess_mu4=preprocess_mu4,
        preprocess_max_iter=preprocess_max_iter,
        preprocess_tol=preprocess_tol,
        preprocessing_mode=preprocessing_mode,
        fail_on_nonconvergence=fail_on_nonconvergence,
        nghost=nghost,
        quicklook=quicklook,
        vmax=vmax,
        selected_indices=indices,
        require_timestamps=require_timestamps,
        progress=progress,
    )
    return metadata

def _prepare_vector_static(
    mode,
    br,
    bt,
    bp,
    output_dir,
    window,
    level,
    geometry,
    preprocess,
    preprocess_mu3,
    preprocess_mu4,
    preprocess_max_iter,
    preprocess_tol,
    nghost,
    quicklook,
    vmax,
    extra=None,
    preprocessing_mode=None,
):
    from .cartesian import (
        boundary_metadata,
        crop_components_with_padding,
        multigrid_reduce_components,
        multigrid_reducer,
    )
    from .fits_io import load_vector_components
    from .preprocessing import VectorPreprocessingConfig, preprocess_vector_magnetogram
    from .writers import write_boundary_frame, write_json, write_static_boundary_products

    bx, by, bz, header = _arrays_or_vector_fits(br, bt, bp, load_vector_components)
    reducer = multigrid_reducer(level)
    crop = _window_for_shape(window, bx.shape)
    bx_crop, by_crop, bz_crop, padding = crop_components_with_padding(
        bx, by, bz,
        x0=crop["x0"], y0=crop["y0"], nx=crop["nx"], ny=crop["ny"],
        pad_x=nghost * reducer, pad_y=nghost * reducer,
    )
    bx_mg, by_mg, bz_mg = multigrid_reduce_components(
        bx_crop, by_crop, bz_crop, level=level
    )
    if preprocessing_mode is None:
        preprocessing_mode = "recommended" if preprocess else "none"
    dx_cm, dy_cm = _static_spacing(header, geometry, level)
    preprocess_config = VectorPreprocessingConfig(
        mode=preprocessing_mode,
        mu3=preprocess_mu3,
        mu4=preprocess_mu4,
        max_iter=preprocess_max_iter,
        tol=preprocess_tol,
        dx=dx_cm / 1.0e5,
        dy=dy_cm / 1.0e5,
        geometry_mode="centered" if preprocessing_mode != "none" else None,
        edge_treatment="nonperiodic" if preprocessing_mode != "none" else None,
    )
    bx_mg, by_mg, bz_mg, preprocessing_audit = preprocess_vector_magnetogram(
        bx_mg, by_mg, bz_mg, config=preprocess_config
    )
    preprocess = preprocess_config.mode != "none"

    meta = _metadata_for_boundary(
        header, geometry, crop, level, bx_mg.shape, nghost, boundary_metadata
    )
    outputs = write_static_boundary_products(
        output_dir, bx_mg, by_mg, bz_mg, meta, quicklook=quicklook, vmax=vmax
    )
    preprocessing_audit_path = Path(output_dir) / "preprocessing_audit.json"
    write_json(preprocessing_audit_path, preprocessing_audit)
    outputs["preprocessing_audit"] = str(preprocessing_audit_path)
    if mode == "data_constrained":
        boundary_frame = write_boundary_frame(
            Path(output_dir) / "B_0001.dat",
            bx_mg,
            by_mg,
            bz_mg,
            snapshot_time=0.0,
            dx=meta.dx_cm / 1.0e5,
            dy=meta.dy_cm / 1.0e5,
        )
        outputs["boundary_frame"] = str(boundary_frame)
    metadata = _base_metadata(mode, crop, level, padding, meta, outputs)
    metadata["preprocess"] = {
        "enabled": bool(preprocess),
        "mu3": preprocess_mu3,
        "mu4": preprocess_mu4,
        "max_iter": preprocess_max_iter,
        "tol": preprocess_tol,
        "mode": preprocess_config.mode,
        "audit": preprocessing_audit,
    }
    if extra:
        metadata.update(extra)
    write_json(Path(output_dir) / "parameters.json", metadata)
    return metadata


def _prepare_unified_boundary(
    input_dir,
    output_dir,
    window,
    level,
    geometry,
    snapshot_index,
    all_frames,
    preprocess,
    preprocess_mu3,
    preprocess_mu4,
    preprocess_max_iter,
    preprocess_tol,
    preprocessing_mode,
    fail_on_nonconvergence,
    nghost,
    quicklook,
    vmax,
    selected_indices=None,
    require_timestamps=False,
    progress=False,
):
    from .cartesian import (
        boundary_metadata,
        crop_components_with_padding,
        multigrid_reduce_components,
        multigrid_reducer,
    )
    from .fits_io import discover_vector_sequence, load_vector_components
    from .preprocessing import VectorPreprocessingConfig, preprocess_vector_magnetogram
    from .writers import ensure_output_dir, write_boundary_outputs, write_json, write_quicklook

    output_dir = ensure_output_dir(output_dir)
    series = discover_vector_sequence(input_dir)
    frame_count = len(series["br"])
    if frame_count == 0:
        raise ValueError("no Br/Bt/Bp frames found")
    snapshot_index = int(snapshot_index)
    if snapshot_index < 0 or snapshot_index >= frame_count:
        raise IndexError("snapshot_index {} is outside 0..{}".format(snapshot_index, frame_count - 1))

    if all_frames:
        if selected_indices is None:
            selected_indices = list(range(frame_count))
        else:
            selected_indices = [int(value) for value in selected_indices]
            if not selected_indices:
                raise ValueError("selected sequence is empty")
            if len(set(selected_indices)) != len(selected_indices):
                raise ValueError("selected sequence contains duplicate frame indices")
            if any(value < 0 or value >= frame_count for value in selected_indices):
                raise IndexError("selected sequence contains an index outside 0..{}".format(frame_count - 1))
            if selected_indices != sorted(selected_indices):
                raise ValueError("selected sequence indices must be strictly increasing")
    else:
        selected_indices = [snapshot_index]
    reducer = multigrid_reducer(level)
    crop = None
    frames = []
    preprocessing_audits = []
    correction_continuity = []
    previous_correction = None
    first_time = None
    meta = None
    padding = None

    first_header = None
    if selected_indices:
        first_index = selected_indices[0]
        _, _, _, first_header = load_vector_components(
            series["br"][first_index], series["bt"][first_index], series["bp"][first_index]
        )
        first_time = _header_time_seconds(first_header)
        if all_frames and require_timestamps and first_time is None:
            raise ValueError("the first selected magnetic frame has no readable observation timestamp")

    iterator = _progress_iterator(selected_indices, progress, "Writing magnetic boundary frames")
    previous_time = None
    expected_shape = None
    expected_spacing = None
    for source_index in iterator:
        br = series["br"][source_index]
        bt = series["bt"][source_index]
        bp = series["bp"][source_index]
        bx, by, bz, header = load_vector_components(br, bt, bp)
        if crop is None:
            crop = _window_for_shape(window, bx.shape)
        bx, by, bz, frame_padding = crop_components_with_padding(
            bx, by, bz,
            x0=crop["x0"], y0=crop["y0"], nx=crop["nx"], ny=crop["ny"],
            pad_x=nghost * reducer, pad_y=nghost * reducer,
        )
        bx, by, bz = multigrid_reduce_components(bx, by, bz, level=level)
        frame_spacing = _sequence_spacing(header, geometry, level)
        if preprocessing_mode is None:
            selected_preprocessing_mode = "recommended" if preprocess else "none"
        else:
            selected_preprocessing_mode = preprocessing_mode
        before_preprocess = [item.copy() for item in (bx, by, bz)]
        preprocess_config = VectorPreprocessingConfig(
            mode=selected_preprocessing_mode,
            mu3=preprocess_mu3,
            mu4=preprocess_mu4,
            max_iter=preprocess_max_iter,
            tol=preprocess_tol,
            dx=frame_spacing[0],
            dy=frame_spacing[1],
            geometry_mode="centered" if selected_preprocessing_mode != "none" else None,
            edge_treatment="nonperiodic" if selected_preprocessing_mode != "none" else None,
            fail_on_nonconvergence=fail_on_nonconvergence,
        )
        bx, by, bz, preprocessing_audit = preprocess_vector_magnetogram(
            bx, by, bz, config=preprocess_config
        )
        preprocess = preprocess_config.mode != "none"
        correction = [new - old for old, new in zip(before_preprocess, (bx, by, bz))]
        correction_norm = float(
            sum(float((item * item).sum()) for item in correction) ** 0.5
        )
        if previous_correction is None:
            correction_delta = None
            correction_delta_relative = None
        else:
            delta = [new - old for old, new in zip(previous_correction, correction)]
            correction_delta = float(
                sum(float((item * item).sum()) for item in delta) ** 0.5
            )
            correction_delta_relative = correction_delta / max(correction_norm, 1.0e-300)
        previous_correction = correction
        if meta is None:
            meta = _metadata_for_boundary(
                header, geometry, crop, level, bx.shape, nghost, boundary_metadata
            )
            padding = frame_padding
            expected_shape = bx.shape
            expected_spacing = frame_spacing
        else:
            if bx.shape != expected_shape:
                raise ValueError(
                    "magnetic-frame grids must have one consistent shape; frame {} has {}, expected {}".format(
                        source_index, bx.shape, expected_shape
                    )
                )
            for label, value, expected in zip(("dx", "dy"), frame_spacing, expected_spacing):
                if abs(value - expected) > 1.0e-10 * max(1.0, abs(expected)):
                    raise ValueError(
                        "magnetic-frame grids must have consistent {}; frame {} has {}, expected {}".format(
                            label, source_index, value, expected
                        )
                    )
        time_seconds = _header_time_seconds(header)
        if all_frames and require_timestamps and time_seconds is None:
            raise ValueError("magnetic frame {} has no readable observation timestamp".format(source_index))
        if all_frames and time_seconds is not None:
            if previous_time is not None and time_seconds <= previous_time:
                raise ValueError("magnetic-frame timestamps must be strictly increasing")
            previous_time = time_seconds
        if first_time is not None and time_seconds is not None:
            snapshot_time = time_seconds - first_time
        elif all_frames:
            snapshot_time = float(source_index)
        else:
            snapshot_time = 0.0
        dx, dy = frame_spacing
        audit_path = output_dir / "preprocessing_audit_frame_{:04d}.json".format(
            len(preprocessing_audits) + 1
        )
        preprocessing_audit.update({
            "frame": {
                "source_index": int(source_index),
                "snapshot_time": float(snapshot_time),
                "observation_time": _header_observation_time(header),
                "dx_km": float(dx),
                "dy_km": float(dy),
            },
            "correction_norm": correction_norm,
            "correction_delta_to_previous": correction_delta,
            "correction_delta_to_previous_relative": correction_delta_relative,
        })
        write_json(audit_path, preprocessing_audit)
        preprocessing_audits.append(str(audit_path))
        correction_continuity.append({
            "source_index": int(source_index),
            "correction_norm": correction_norm,
            "delta_to_previous": correction_delta,
            "delta_to_previous_relative": correction_delta_relative,
        })
        frames.append({
            "bx": bx,
            "by": by,
            "bz": bz,
            "time": snapshot_time,
            "observation_time": _header_observation_time(header),
            "dx": dx,
            "dy": dy,
            "source_index": source_index,
            "source": {"Br": str(br), "Bt": str(bt), "Bp": str(bp)},
        })

    paths = write_boundary_outputs(output_dir, frames, prefix="B")
    if quicklook and frames:
        quicklook_frames = [("first", frames[0])]
        if len(frames) > 2:
            quicklook_frames.append(("middle", frames[len(frames) // 2]))
        if len(frames) > 1:
            quicklook_frames.append(("last", frames[-1]))
        for label, frame in quicklook_frames:
            write_quicklook(
                output_dir / "boundary_{}_quicklook.png".format(label),
                frame["bx"], frame["by"], frame["bz"], vmax=vmax,
            )

    metadata = {
        "mode": "boundary_sequence" if all_frames else "boundary_frame",
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "window": crop,
        "level": level,
        "nghost": nghost,
        "padding": padding,
        "preprocess": bool(preprocess),
        "preprocessing": {
            "mode": selected_preprocessing_mode if selected_indices else ("recommended" if preprocess else "none"),
            "audit_files": preprocessing_audits,
            "correction_continuity": correction_continuity,
        },
        "quicklook_vmax": float(vmax),
        "frame_format": "snapshot_time,nx,ny,dx,dy,Bx,By,Bz",
        "spacing_unit": "km",
        "time_unit": "s",
        "snapshot_index": None if all_frames else snapshot_index,
        "input_frame_count": frame_count,
        "output_frame_count": len(frames),
        "outputs": [str(path) for path in paths],
        "amrvac": meta.as_dict if meta is not None else None,
        "frames": [
            {
                "snapshot_time": frame["time"],
                "observation_time": frame["observation_time"],
                "dx": frame["dx"],
                "dy": frame["dy"],
                "source_index": frame["source_index"],
                "source": frame["source"],
            }
            for frame in frames
        ],
    }
    write_json(output_dir / "boundary_parameters.json", metadata)
    return metadata


def _prepare_vector_sequence(
    mode,
    input_dir,
    output_dir,
    window,
    level,
    geometry,
    preprocess,
    quicklook,
):
    from .cartesian import crop_components, multigrid_reduce_components
    from .fits_io import discover_vector_sequence, load_vector_components
    from .preprocessing import VectorPreprocessingConfig, preprocess_vector_magnetogram
    from .writers import ensure_output_dir, write_boundary_outputs, write_json, write_quicklook

    output_dir = ensure_output_dir(output_dir)
    series = discover_vector_sequence(input_dir)
    crop = None
    frames = []
    preprocessing_audits = []
    correction_continuity = []
    previous_correction = None
    first_time = None
    for index, (br, bt, bp) in enumerate(zip(series["br"], series["bt"], series["bp"])):
        bx, by, bz, header = load_vector_components(br, bt, bp)
        if crop is None:
            crop = _window_for_shape(window, bx.shape)
        bx, by, bz = crop_components(
            bx, by, bz, x0=crop["x0"], y0=crop["y0"], nx=crop["nx"], ny=crop["ny"]
        )
        bx, by, bz = multigrid_reduce_components(bx, by, bz, level=level)
        time_seconds = _header_time_seconds(header)
        if time_seconds is None:
            time_seconds = float(index)
        if first_time is None:
            first_time = time_seconds
        dx, dy = _sequence_spacing(header, geometry, level)
        before_preprocess = [item.copy() for item in (bx, by, bz)]
        preprocess_config = VectorPreprocessingConfig(
            mode="recommended" if preprocess else "none",
            dx=dx,
            dy=dy,
            geometry_mode="centered" if preprocess else None,
            edge_treatment="nonperiodic" if preprocess else None,
        )
        bx, by, bz, preprocessing_audit = preprocess_vector_magnetogram(
            bx, by, bz, config=preprocess_config
        )
        correction = [new - old for old, new in zip(before_preprocess, (bx, by, bz))]
        correction_norm = float(
            sum(float((item * item).sum()) for item in correction) ** 0.5
        )
        if previous_correction is None:
            correction_delta = None
            correction_delta_relative = None
        else:
            delta = [new - old for old, new in zip(previous_correction, correction)]
            correction_delta = float(
                sum(float((item * item).sum()) for item in delta) ** 0.5
            )
            correction_delta_relative = correction_delta / max(correction_norm, 1.0e-300)
        previous_correction = correction
        audit_path = output_dir / "preprocessing_audit_frame_{:04d}.json".format(index + 1)
        preprocessing_audit.update({
            "frame": {
                "source_index": int(index),
                "snapshot_time": float(time_seconds - first_time),
                "observation_time": _header_observation_time(header),
                "dx_km": float(dx),
                "dy_km": float(dy),
            },
            "correction_norm": correction_norm,
            "correction_delta_to_previous": correction_delta,
            "correction_delta_to_previous_relative": correction_delta_relative,
        })
        write_json(audit_path, preprocessing_audit)
        preprocessing_audits.append(str(audit_path))
        correction_continuity.append({
            "source_index": int(index),
            "correction_norm": correction_norm,
            "delta_to_previous": correction_delta,
            "delta_to_previous_relative": correction_delta_relative,
        })
        frames.append({
            "bx": bx,
            "by": by,
            "bz": bz,
            "time": time_seconds - first_time,
            "dx": dx,
            "dy": dy,
            "source": {"Br": str(br), "Bt": str(bt), "Bp": str(bp)},
        })
    paths = write_boundary_outputs(output_dir, frames, prefix="B")
    if quicklook and frames:
        write_quicklook(output_dir / "sequence_first_quicklook.png", frames[0]["bx"], frames[0]["by"], frames[0]["bz"])
    metadata = {
        "mode": mode,
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "window": crop,
        "level": level,
        "preprocess": bool(preprocess),
        "preprocessing": {
            "mode": "recommended" if preprocess else "none",
            "audit_files": preprocessing_audits,
            "correction_continuity": correction_continuity,
        },
        "sequence_format": "snapshot_time,nx,ny,dx,dy,Bx,By,Bz",
        "time_scaling": "not_written_python_v1",
        "frame_count": len(frames),
        "outputs": [str(path) for path in paths],
        "frames": [
            {"time": frame["time"], "dx": frame["dx"], "dy": frame["dy"], "source": frame["source"]}
            for frame in frames
        ],
    }
    write_json(output_dir / "sequence_parameters.json", metadata)
    return metadata


def _raw_hmi_frames_from_input(raw_hmi_input):
    if isinstance(raw_hmi_input, (str, Path)):
        from .fits_io import inspect_magnetic_input

        info = inspect_magnetic_input(raw_hmi_input)
    else:
        info = raw_hmi_input

    if isinstance(info, dict) and "raw_hmi_sequence" in info:
        frames = info["raw_hmi_sequence"].get("complete", [])
    elif isinstance(info, dict) and "raw_hmi" in info:
        frames = [info["raw_hmi"]]
    else:
        frames = info

    normalized = []
    for index, frame in enumerate(frames, start=1):
        if frame is None:
            continue
        missing = [segment for segment in ("field", "inclination", "azimuth", "disambig") if not frame.get(segment)]
        if missing:
            raise ValueError("raw HMI frame {} is missing {}".format(index, ", ".join(missing)))
        normalized.append({
            "key": frame.get("key", "frame_{:04d}".format(index)),
            "field": Path(frame["field"]),
            "inclination": Path(frame["inclination"]),
            "azimuth": Path(frame["azimuth"]),
            "disambig": Path(frame["disambig"]),
        })
    return normalized


def _progress_iterator(items, enabled, desc):
    if not enabled:
        return items
    try:
        from tqdm.auto import tqdm
    except ImportError:
        total = len(items)

        def _printing_iterator():
            for index, item in enumerate(items, start=1):
                print("{} {}/{}".format(desc, index, total))
                yield item

        return _printing_iterator()
    return tqdm(items, desc=desc, unit="frame")


def _normalize_window(window):
    if isinstance(window, dict):
        return {
            "x0": int(window["x0"]),
            "y0": int(window["y0"]),
            "nx": int(window["nx"]),
            "ny": int(window["ny"]),
        }
    if len(window) != 4:
        raise ValueError("window must be a dict or (x0, y0, nx, ny)")
    return {"x0": int(window[0]), "y0": int(window[1]), "nx": int(window[2]), "ny": int(window[3])}


def _clean_output_prefix(prefix):
    if prefix is None:
        return None
    text = str(prefix).strip()
    if not text:
        return None
    return "".join(char if char.isalnum() or char in ("-", "_") else "_" for char in text)


def _prefixed_output_path(output_dir, prefix, filename):
    if prefix is None:
        return output_dir / filename
    return output_dir / "{}.{}".format(prefix, filename)


def _window_for_shape(window, shape):
    if window is None:
        return {"x0": 0, "y0": 0, "nx": int(shape[1]), "ny": int(shape[0])}
    return _normalize_window(window)


def _array_or_br_fits(value, loader):
    try:
        import numpy as np
    except ImportError as error:
        raise ImportError("numpy is required for data-driven preprocessing") from error
    if isinstance(value, (str, Path)):
        return loader(value)
    return np.asarray(value, dtype=float), None


def _arrays_or_vector_fits(br, bt, bp, loader):
    try:
        import numpy as np
    except ImportError as error:
        raise ImportError("numpy is required for data-driven preprocessing") from error
    if all(isinstance(value, (str, Path)) for value in (br, bt, bp)):
        return loader(br, bt, bp)
    bx = np.asarray(bp, dtype=float)
    by = -np.asarray(bt, dtype=float)
    bz = np.asarray(br, dtype=float)
    return bx, by, bz, None


def _metadata_for_boundary(header, geometry, crop, level, shape, nghost, boundary_metadata):
    dx_cm, dy_cm = _static_spacing(header, geometry, level)
    xc_cm, yc_cm = _window_center(header, geometry, crop)
    geometry = geometry or {}
    boundary_z = geometry.get("boundary_z_10mm", geometry.get("z_min_10mm", 0.0))
    return boundary_metadata(
        nx=shape[1],
        ny=shape[0],
        xc_cm=xc_cm,
        yc_cm=yc_cm,
        dx_cm=dx_cm,
        dy_cm=dy_cm,
        nghost=nghost,
        z_min_10mm=boundary_z,
    )


def _static_spacing(header, geometry, level):
    import numpy as np

    geometry = geometry or {}
    reducer = 2 ** (int(level) - 1)
    if "dx_cm" in geometry:
        dx_cm = float(geometry["dx_cm"])
    elif header is not None:
        dx_cm = SOLAR_RADIUS_CM * np.radians(abs(float(header.get("CDELT1")))) * reducer
    else:
        raise ValueError("geometry['dx_cm'] is required for array inputs")
    if "dy_cm" in geometry:
        dy_cm = float(geometry["dy_cm"])
    elif header is not None:
        dy_cm = SOLAR_RADIUS_CM * np.radians(abs(float(header.get("CDELT2", header.get("CDELT1"))))) * reducer
    else:
        dy_cm = dx_cm
    return dx_cm, dy_cm


def _sequence_spacing(header, geometry, level):
    dx_cm, dy_cm = _static_spacing(header, geometry, level)
    return dx_cm / 1.0e5, dy_cm / 1.0e5


def _window_center(header, geometry, crop):
    import numpy as np

    geometry = geometry or {}
    if "xc_cm" in geometry and "yc_cm" in geometry:
        return float(geometry["xc_cm"]), float(geometry["yc_cm"])
    if header is None:
        return float(geometry.get("xc_cm", 0.0)), float(geometry.get("yc_cm", 0.0))
    x_center = crop["x0"] + 0.5 * (crop["nx"] - 1)
    y_center = crop["y0"] + 0.5 * (crop["ny"] - 1)
    crpix1 = float(header.get("CRPIX1", (int(header["NAXIS1"]) + 1.0) / 2.0))
    crpix2 = float(header.get("CRPIX2", (int(header["NAXIS2"]) + 1.0) / 2.0))
    cdelt1 = float(header.get("CDELT1"))
    cdelt2 = float(header.get("CDELT2", header.get("CDELT1")))
    xc_cm = SOLAR_RADIUS_CM * np.radians((x_center + 1.0 - crpix1) * cdelt1)
    yc_cm = SOLAR_RADIUS_CM * np.radians((y_center + 1.0 - crpix2) * cdelt2)
    return xc_cm, yc_cm


def _base_metadata(mode, crop, level, padding, meta, outputs):
    return {
        "mode": mode,
        "window": crop,
        "level": level,
        "padding": padding,
        "amrvac": meta.as_dict,
        "outputs": outputs,
    }


def _header_time_seconds(header):
    if header is None:
        return None
    text = header.get("T_OBS") or header.get("DATE-OBS") or header.get("T_REC")
    if not text:
        return None
    candidates = [str(text)]
    if len(str(text)) >= 19:
        candidates.append(str(text)[:19])
    for candidate in candidates:
        for fmt in ("%Y.%m.%d_%H:%M:%S", "%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S"):
            try:
                return datetime.strptime(candidate, fmt).timestamp()
            except ValueError:
                continue
    return None


def _header_observation_time(header):
    """Return the original observation-time text for plotting and manifests."""

    if header is None:
        return None
    value = header.get("T_OBS") or header.get("DATE-OBS") or header.get("T_REC")
    return None if value in (None, "") else str(value)


def _sample_native_components(
    source_map,
    bxi,
    beta,
    bzeta,
    grid,
    sampling_mode,
    interpolation_order,
    fill_value,
    oversample_resolution_degree,
    smooth_sigma_degree,
    smooth_truncate,
    gaussian_block_reduce,
):
    from .cea import _sample_at_source_pixels, _source_pixels_from_heliographic

    if sampling_mode == "direct":
        pixel_yx = _source_pixels_from_heliographic(source_map, grid.lon, grid.lat)
        return (
            _sample_at_source_pixels(bxi, pixel_yx, grid.shape, interpolation_order, fill_value),
            _sample_at_source_pixels(beta, pixel_yx, grid.shape, interpolation_order, fill_value),
            _sample_at_source_pixels(bzeta, pixel_yx, grid.shape, interpolation_order, fill_value),
            {"mode": "direct", "oversample_factor": 1, "interpolation_order": interpolation_order},
        )
    if sampling_mode != "sharp":
        raise ValueError("sampling_mode must be 'direct' or 'sharp'")
    from .cea import make_cea_patch_grid

    factor = int(round(grid.resolution_degree / oversample_resolution_degree))
    if factor < 1:
        raise ValueError("oversample_resolution_degree must be no larger than resolution_degree")
    high_resolution = grid.resolution_degree / factor
    high_grid = make_cea_patch_grid(
        grid.center_lon,
        grid.center_lat,
        grid.shape[1] * grid.resolution_degree,
        grid.shape[0] * grid.resolution_degree,
        high_resolution,
    )
    pixel_yx = _source_pixels_from_heliographic(source_map, high_grid.lon, high_grid.lat)
    high_bxi = _sample_at_source_pixels(bxi, pixel_yx, high_grid.shape, interpolation_order, fill_value)
    high_beta = _sample_at_source_pixels(beta, pixel_yx, high_grid.shape, interpolation_order, fill_value)
    high_bzeta = _sample_at_source_pixels(bzeta, pixel_yx, high_grid.shape, interpolation_order, fill_value)
    sigma_pixel = smooth_sigma_degree / high_resolution
    return (
        gaussian_block_reduce(high_bxi, factor, factor, sigma_y=sigma_pixel, sigma_x=sigma_pixel, truncate=smooth_truncate, cval=fill_value),
        gaussian_block_reduce(high_beta, factor, factor, sigma_y=sigma_pixel, sigma_x=sigma_pixel, truncate=smooth_truncate, cval=fill_value),
        gaussian_block_reduce(high_bzeta, factor, factor, sigma_y=sigma_pixel, sigma_x=sigma_pixel, truncate=smooth_truncate, cval=fill_value),
        {
            "mode": "sharp",
            "oversample_factor": factor,
            "interpolation_order": interpolation_order,
            "smooth_sigma_degree": smooth_sigma_degree,
            "smooth_truncate": smooth_truncate,
        },
    )


def _sample_nearest_segment(segment_map, grid, fill_value):
    from .cea import _sample_at_source_pixels, _source_pixels_from_heliographic

    pixel_yx = _source_pixels_from_heliographic(segment_map, grid.lon, grid.lat)
    return _sample_at_source_pixels(
        segment_map.data, pixel_yx, grid.shape, order=0, cval=fill_value
    )


def _write_component_fits(path, data, source_map, grid, component, disk_lon, disk_lat, p_angle):
    import numpy as np
    from astropy.io import fits

    header = fits.Header()
    header["SIMPLE"] = True
    header["BITPIX"] = -64
    header["NAXIS"] = 2
    header["NAXIS1"] = grid.shape[1]
    header["NAXIS2"] = grid.shape[0]
    header["CTYPE1"] = "HGLN-CEA"
    header["CTYPE2"] = "HGLT-CEA"
    header["CUNIT1"] = "deg"
    header["CUNIT2"] = "deg"
    header["CRPIX1"] = (grid.shape[1] + 1.0) / 2.0
    header["CRPIX2"] = (grid.shape[0] + 1.0) / 2.0
    header["CRVAL1"] = grid.center_lon
    header["CRVAL2"] = grid.center_lat
    header["CDELT1"] = grid.resolution_degree
    header["CDELT2"] = grid.resolution_degree
    header["CROTA2"] = 0.0
    header["BUNIT"] = "G"
    header["CONTENT"] = "CEA {}".format(component)
    header["COMP"] = component
    header["DATE-OBS"] = str(source_map.date.isot)
    header["T_OBS"] = source_map.date.strftime("%Y.%m.%d_%H:%M:%S")
    header["HGLN_OBS"] = float(disk_lon)
    header["HGLT_OBS"] = float(disk_lat)
    header["P_ANGLE"] = float(p_angle)
    fits.PrimaryHDU(np.asarray(data, dtype=np.float64), header=header).writeto(str(path), overwrite=True)
