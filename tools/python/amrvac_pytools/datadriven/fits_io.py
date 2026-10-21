"""FITS and component discovery helpers for data-driven preprocessing."""

from __future__ import print_function

import re
from pathlib import Path


COMPONENT_SUFFIXES = {
    "br": ("Br.fits", ".Br.fits", "_Br.fits", ".br.fits", "_br.fits"),
    "bt": ("Bt.fits", ".Bt.fits", "_Bt.fits", ".bt.fits", "_bt.fits"),
    "bp": ("Bp.fits", ".Bp.fits", "_Bp.fits", ".bp.fits", "_bp.fits"),
}


HMI_SEGMENT_NAMES = ("field", "inclination", "azimuth", "disambig")


def read_fits_data(path):
    """Read the first 2D image from a FITS file."""

    try:
        import numpy as np
        from astropy.io import fits
    except ImportError as error:
        raise ImportError("astropy and numpy are required to read FITS files") from error

    path = Path(path)
    with fits.open(str(path)) as hdul:
        for hdu in hdul:
            if hdu.data is not None and np.ndim(hdu.data) == 2:
                return np.asarray(hdu.data, dtype=float), hdu.header.copy()
    raise ValueError("no two-dimensional image data found in {}".format(path))


def resolve_component_paths(br=None, bt=None, bp=None, cea_dir=None):
    """Resolve one Br/Bt/Bp triplet from explicit paths or a CEA directory."""

    if cea_dir is not None:
        cea_dir = Path(cea_dir)
        br = br or cea_dir / "Br.fits"
        bt = bt or cea_dir / "Bt.fits"
        bp = bp or cea_dir / "Bp.fits"
    paths = {"br": br, "bt": bt, "bp": bp}
    missing = [name for name, path in paths.items() if path is None or not Path(path).exists()]
    if missing:
        raise FileNotFoundError("missing component FITS: {}".format(", ".join(missing)))
    return Path(br), Path(bt), Path(bp)


def discover_vector_sequence(directory):
    """Discover sorted Br/Bt/Bp files in a sequence directory."""

    directory = Path(directory)
    files = [path for path in directory.iterdir() if path.is_file()]
    grouped = _group_component_frames(_find_component_files(files))
    if grouped["duplicates"]:
        duplicate = grouped["duplicates"][0]
        raise ValueError(
            "duplicate {} files for magnetic frame {}".format(
                duplicate["component"], duplicate["key"]
            )
        )
    if not grouped["complete"]:
        raise ValueError("no Br/Bt/Bp sequence found in {}".format(directory))
    return {
        component: [frame[component] for frame in grouped["complete"]]
        for component in ("br", "bt", "bp")
    }


def inspect_magnetic_input(directory):
    """Inspect one directory and classify the magnetic input products.

    The notebook uses this as the first user-facing gate. It intentionally
    relies on filename conventions only; FITS headers are read later by the
    selected workflow.
    """

    directory = Path(directory)
    if not directory.exists():
        raise FileNotFoundError("input directory does not exist: {}".format(directory))
    if not directory.is_dir():
        raise NotADirectoryError("input path is not a directory: {}".format(directory))

    files = sorted([path for path in directory.iterdir() if path.is_file()])
    component_files = _find_component_files(files)
    component_sequence = _group_component_frames(component_files)
    raw_hmi_sequence = _find_raw_hmi_sequence(files)
    raw_hmi_files = _first_raw_hmi_frame(raw_hmi_sequence)

    has_vector = bool(component_sequence["complete"])
    has_raw_hmi = all(raw_hmi_files[name] is not None for name in HMI_SEGMENT_NAMES)

    if has_vector:
        counts = {name: len(component_sequence["complete"]) for name in ("br", "bt", "bp")}
        if component_sequence["duplicates"]:
            kind = "unknown"
            duplicate = component_sequence["duplicates"][0]
            message = "Duplicate {} files were found for frame {}.".format(
                duplicate["component"], duplicate["key"]
            )
        elif counts["br"] == 1:
            kind = "cea_vector"
            message = "Found one CEA/SHARP-like Br/Bt/Bp triplet."
        else:
            kind = "cea_vector_sequence"
            message = "Found a CEA/SHARP-like Br/Bt/Bp sequence with {} frames.".format(counts["br"])
        return {
            "kind": kind,
            "directory": str(directory),
            "message": message,
            "components": {
                key: [str(frame[key]) for frame in component_sequence["complete"]]
                for key in ("br", "bt", "bp")
            },
            "component_sequence": _serializable_component_sequence(component_sequence),
            "raw_hmi": {key: str(value) if value is not None else None for key, value in raw_hmi_files.items()},
            "raw_hmi_sequence": _serializable_raw_hmi_sequence(raw_hmi_sequence),
        }

    if component_files["br"] and not component_files["bt"] and not component_files["bp"]:
        kind = "br_only" if len(component_files["br"]) == 1 else "br_sequence"
        return {
            "kind": kind,
            "directory": str(directory),
            "message": "Found Br-only input. This is suitable for potential or Bn-only workflows.",
            "components": {key: [str(path) for path in value] for key, value in component_files.items()},
            "component_sequence": _serializable_component_sequence(component_sequence),
            "raw_hmi": {key: str(value) if value is not None else None for key, value in raw_hmi_files.items()},
            "raw_hmi_sequence": _serializable_raw_hmi_sequence(raw_hmi_sequence),
        }

    if has_raw_hmi:
        complete_count = len(raw_hmi_sequence["complete"])
        if raw_hmi_sequence["duplicates"]:
            duplicate = raw_hmi_sequence["duplicates"][0]
            kind = "unknown"
            message = "Duplicate {} files were found for raw HMI frame {}.".format(
                duplicate["segment"], duplicate["key"]
            )
        else:
            kind = "raw_hmi_vector" if complete_count == 1 else "raw_hmi_vector_sequence"
            message = "Found raw HMI vector segments: field, inclination, azimuth, disambig."
            if complete_count > 1:
                message = "Found a raw HMI vector segment sequence with {} complete frames.".format(complete_count)
        return {
            "kind": kind,
            "directory": str(directory),
            "message": message,
            "components": {key: [str(path) for path in value] for key, value in component_files.items()},
            "component_sequence": _serializable_component_sequence(component_sequence),
            "raw_hmi": {key: str(value) if value is not None else None for key, value in raw_hmi_files.items()},
            "raw_hmi_sequence": _serializable_raw_hmi_sequence(raw_hmi_sequence),
        }

    return {
        "kind": "unknown",
        "directory": str(directory),
        "message": (
            "Could not identify magnetic input. Expected either Br/Bt/Bp FITS files "
            "or raw HMI field/inclination/azimuth/disambig FITS files."
        ),
        "components": {key: [str(path) for path in value] for key, value in component_files.items()},
        "component_sequence": _serializable_component_sequence(component_sequence),
        "raw_hmi": {key: str(value) if value is not None else None for key, value in raw_hmi_files.items()},
        "raw_hmi_sequence": _serializable_raw_hmi_sequence(raw_hmi_sequence),
    }


def discover_magnetic_inputs(directory, include_subdirectories=True):
    """Return recognized magnetic-input candidates under one directory."""

    directory = Path(directory)
    candidates = []
    roots = [directory]
    if include_subdirectories and directory.is_dir():
        roots.extend(sorted([path for path in directory.iterdir() if path.is_dir()]))
    for root in roots:
        try:
            info = inspect_magnetic_input(root)
        except (FileNotFoundError, NotADirectoryError):
            continue
        if info["kind"] != "unknown":
            info["label"] = "{} ({})".format(root.name, info["kind"])
            candidates.append(info)
    return candidates


def summarize_magnetic_input(info, max_warnings=8):
    """Return a compact, notebook-friendly summary for inspected input."""

    kind = info["kind"]
    summary = {
        "kind": kind,
        "directory": info["directory"],
        "warnings": [],
    }
    if kind in ("cea_vector", "cea_vector_sequence", "br_only", "br_sequence", "unknown"):
        components = info["components"]
        counts = {name: len(components.get(name, [])) for name in ("br", "bt", "bp")}
        summary["component_counts"] = counts
        summary["frame_count"] = min([count for count in counts.values() if count > 0] or [0])
        br_files = components.get("br", [])
        if br_files:
            summary["first_br"] = Path(br_files[0]).name
            summary["last_br"] = Path(br_files[-1]).name
        component_sequence = info.get("component_sequence", {})
        incomplete = component_sequence.get("incomplete", [])
        if incomplete:
            summary["warnings"].append(
                "{} incomplete Br/Bt/Bp frame(s) detected.".format(len(incomplete))
            )
            for frame in incomplete[:max_warnings]:
                summary["warnings"].append(
                    "frame {} is missing {}".format(frame["key"], ", ".join(frame["missing"]))
                )
            if len(incomplete) > max_warnings:
                summary["warnings"].append(
                    "... {} more incomplete frame(s) omitted".format(len(incomplete) - max_warnings)
                )
        raw_incomplete = info.get("raw_hmi_sequence", {}).get("incomplete", [])
        if raw_incomplete:
            summary["warnings"].append(
                "{} incomplete raw HMI frame(s) detected.".format(len(raw_incomplete))
            )
            for frame in raw_incomplete[:max_warnings]:
                summary["warnings"].append(
                    "frame {} is missing {}".format(frame["key"], ", ".join(frame["missing"]))
                )
            if len(raw_incomplete) > max_warnings:
                summary["warnings"].append(
                    "... {} more incomplete raw HMI frame(s) omitted".format(
                        len(raw_incomplete) - max_warnings
                    )
                )
        return summary
    if kind in ("raw_hmi_vector", "raw_hmi_vector_sequence"):
        raw_sequence = info.get("raw_hmi_sequence", {})
        complete = raw_sequence.get("complete", [])
        incomplete = raw_sequence.get("incomplete", [])
        first_frame = complete[0] if complete else {}
        present = [name for name in HMI_SEGMENT_NAMES if first_frame.get(name) is not None]
        summary["segment_count"] = len(present)
        summary["segments"] = list(HMI_SEGMENT_NAMES)
        summary["frame_count"] = len(complete)
        if complete:
            summary["first_frame"] = complete[0].get("key")
            summary["last_frame"] = complete[-1].get("key")
        if incomplete:
            summary["warnings"].append(
                "{} incomplete raw HMI frame(s) detected.".format(len(incomplete))
            )
            for frame in incomplete[:max_warnings]:
                summary["warnings"].append(
                    "frame {} is missing {}".format(frame["key"], ", ".join(frame["missing"]))
                )
            if len(incomplete) > max_warnings:
                summary["warnings"].append(
                    "... {} more incomplete frame(s) omitted".format(len(incomplete) - max_warnings)
                )
        return summary
    return summary


def load_magnetic_preview(info):
    """Load a first-frame Br-like image for notebook preview.

    For CEA/SHARP-like inputs this is the first Br FITS image. For raw HMI
    vector segments this is the native line-of-sight component proxy
    `field*cos(inclination)`, before CEA remapping.
    """

    try:
        import numpy as np
    except ImportError as error:
        raise ImportError("numpy is required to preview magnetic inputs") from error

    kind = info["kind"]
    if kind in ("cea_vector", "cea_vector_sequence", "br_only", "br_sequence"):
        br_paths = info["components"]["br"]
        if not br_paths:
            raise ValueError("no Br file is available for preview")
        data, header = read_fits_data(br_paths[0])
        time_label = _time_label_from_name(Path(br_paths[0]).name)
        title_parts = [r"CEA magnetogram, $B_r$"]
        if time_label is not None:
            title_parts.append(time_label)
        return {
            "data": data,
            "header": header,
            "title": ", ".join(title_parts),
            "source": br_paths[0],
            "quantity": "Br",
            "quantity_label": r"$B_r$",
        }
    if kind in ("raw_hmi_vector", "raw_hmi_vector_sequence"):
        field, header = read_fits_data(info["raw_hmi"]["field"])
        inclination, _ = read_fits_data(info["raw_hmi"]["inclination"])
        proxy = field * np.cos(np.radians(inclination))
        time_label = _time_label_from_name(Path(info["raw_hmi"]["field"]).name)
        title = r"HMI magnetogram, $B_{\mathrm{LOS}}$"
        if time_label is not None:
            title = "{}, {}".format(title, time_label)
        return {
            "data": proxy,
            "header": header,
            "title": title,
            "source": info["raw_hmi"]["field"],
            "quantity": "native_bzeta_proxy",
            "quantity_label": r"$B_{\mathrm{LOS}}$",
        }
    raise ValueError("cannot preview unknown input kind: {}".format(kind))


def load_cea_patch_preview(info, cea_patch, preview_data=None):
    """Project a fixed CEA patch onto the first input frame for preview.

    This is primarily used by the notebook for raw HMI inputs: the user selects
    a heliographic CEA patch, while the first preview image is still on the
    native HMI image plane. The returned outline is in source-image pixel
    coordinates, and the returned patch image is an approximate preview sampled
    from the first-frame Br-like preview quantity.
    """

    try:
        import astropy.units as u
        import numpy as np
        from astropy.coordinates import SkyCoord
        from scipy.ndimage import map_coordinates
        from sunpy.coordinates import frames
        from sunpy.map import Map
    except ImportError as error:
        raise ImportError(
            "astropy, scipy, sunpy, and numpy are required for CEA patch preview"
        ) from error

    from .cea import (
        _sample_at_source_pixels,
        _source_pixels_from_heliographic,
        make_cea_patch_grid,
    )

    if preview_data is None:
        preview = load_magnetic_preview(info)
        data = preview["data"]
        header = preview["header"]
    else:
        preview = load_magnetic_preview(info)
        data = np.asarray(preview_data, dtype=float)
        header = preview["header"]

    source_map = Map(data, header)
    grid = make_cea_patch_grid(
        cea_patch["center_lon"],
        cea_patch["center_lat"],
        cea_patch["width_degree"],
        cea_patch["height_degree"],
        cea_patch.get("resolution_degree", 0.03),
    )

    patch_yx = _source_pixels_from_heliographic(source_map, grid.lon, grid.lat)
    sampled = _sample_at_source_pixels(
        data, patch_yx, grid.shape, order=1, cval=np.nan, map_coordinates=map_coordinates
    )

    lon_outline = np.concatenate([
        grid.lon[0, :],
        grid.lon[:, -1],
        grid.lon[-1, ::-1],
        grid.lon[::-1, 0],
    ])
    lat_outline = np.concatenate([
        grid.lat[0, :],
        grid.lat[:, -1],
        grid.lat[-1, ::-1],
        grid.lat[::-1, 0],
    ])
    outline_yx = _source_pixels_from_heliographic(
        source_map,
        lon_outline,
        lat_outline,
        SkyCoord=SkyCoord,
        u=u,
        frames=frames,
    )

    return {
        "data": sampled,
        "outline_x": outline_yx[1].reshape(-1),
        "outline_y": outline_yx[0].reshape(-1),
        "grid_shape": grid.shape,
        "title": "Fixed CEA patch preview",
    }


def load_vector_components(br, bt, bp):
    """Load Br/Bt/Bp FITS and return local Cartesian Bx/By/Bz arrays."""

    from .cartesian import local_cartesian_from_heliographic

    br_data, header = read_fits_data(br)
    bt_data, _ = read_fits_data(bt)
    bp_data, _ = read_fits_data(bp)
    bx, by, bz = local_cartesian_from_heliographic(br_data, bt_data, bp_data)
    return bx, by, bz, header


def load_br(br):
    """Load a Br FITS image."""

    return read_fits_data(br)


def _find_component_files(files):
    series = {}
    for component, suffixes in COMPONENT_SUFFIXES.items():
        series[component] = sorted(
            [path for path in files if any(path.name.endswith(suffix) for suffix in suffixes)],
            key=_sequence_sort_key,
        )
    return series


def _find_raw_hmi_files(files):
    return _first_raw_hmi_frame(_find_raw_hmi_sequence(files))


def _find_raw_hmi_sequence(files):
    frames = {}
    duplicates = []
    for path in sorted(files, key=_sequence_sort_key):
        name = path.name.lower()
        if not (name.endswith(".fits") or name.endswith(".fit") or name.endswith(".fts")):
            continue
        for segment in HMI_SEGMENT_NAMES:
            if _looks_like_hmi_segment(name, segment):
                key = _raw_hmi_frame_key(path, segment)
                frame = frames.setdefault(key, {"key": key})
                if segment in frame:
                    duplicates.append({
                        "key": key,
                        "segment": segment,
                        "paths": [frame[segment], path],
                    })
                else:
                    frame[segment] = path
                break

    complete = []
    incomplete = []
    for key in sorted(frames, key=_frame_key_sort_key):
        frame = frames[key]
        missing = [segment for segment in HMI_SEGMENT_NAMES if segment not in frame]
        if missing:
            incomplete.append({"key": key, "missing": missing, "present": dict(frame)})
        else:
            complete.append(frame)
    return {"complete": complete, "incomplete": incomplete, "duplicates": duplicates}


def _first_raw_hmi_frame(raw_hmi_sequence):
    found = {name: None for name in HMI_SEGMENT_NAMES}
    if not raw_hmi_sequence["complete"]:
        return found
    first = raw_hmi_sequence["complete"][0]
    for segment in HMI_SEGMENT_NAMES:
        found[segment] = first[segment]
    return found


def _serializable_raw_hmi_sequence(raw_hmi_sequence):
    complete = []
    for frame in raw_hmi_sequence["complete"]:
        serial = {"key": frame["key"]}
        for segment in HMI_SEGMENT_NAMES:
            serial[segment] = str(frame[segment])
        complete.append(serial)
    incomplete = []
    for frame in raw_hmi_sequence["incomplete"]:
        present = {
            segment: str(path)
            for segment, path in frame["present"].items()
            if segment in HMI_SEGMENT_NAMES
        }
        incomplete.append({
            "key": frame["key"],
            "missing": list(frame["missing"]),
            "present": present,
        })
    duplicates = [{
        "key": frame["key"],
        "segment": frame["segment"],
        "paths": [str(path) for path in frame["paths"]],
    } for frame in raw_hmi_sequence.get("duplicates", [])]
    return {"complete": complete, "incomplete": incomplete, "duplicates": duplicates}


def _group_component_frames(component_files):
    frames = {}
    duplicates = []
    for component in ("br", "bt", "bp"):
        for path in component_files.get(component, []):
            key = _component_frame_key(path, component)
            frame = frames.setdefault(key, {"key": key})
            if component in frame:
                duplicates.append({
                    "key": key,
                    "component": component,
                    "paths": [frame[component], path],
                })
            else:
                frame[component] = path
    complete = []
    incomplete = []
    for key in sorted(frames, key=_frame_key_sort_key):
        frame = frames[key]
        missing = [component for component in ("br", "bt", "bp") if component not in frame]
        if missing:
            incomplete.append({"key": key, "missing": missing, "present": dict(frame)})
        else:
            complete.append(frame)
    return {"complete": complete, "incomplete": incomplete, "duplicates": duplicates}


def _serializable_component_sequence(component_sequence):
    complete = [{
        "key": frame["key"],
        "br": str(frame["br"]),
        "bt": str(frame["bt"]),
        "bp": str(frame["bp"]),
    } for frame in component_sequence["complete"]]
    incomplete = [{
        "key": frame["key"],
        "missing": list(frame["missing"]),
        "present": {
            component: str(path)
            for component, path in frame["present"].items()
            if component in ("br", "bt", "bp")
        },
    } for frame in component_sequence["incomplete"]]
    duplicates = [{
        "key": frame["key"],
        "component": frame["component"],
        "paths": [str(path) for path in frame["paths"]],
    } for frame in component_sequence["duplicates"]]
    return {"complete": complete, "incomplete": incomplete, "duplicates": duplicates}


def _looks_like_hmi_segment(filename, segment):
    tokens = re.split(r"[^a-z0-9]+", filename.lower())
    if segment in tokens:
        return True
    return ".{}.".format(segment) in filename.lower()


def _raw_hmi_frame_key(path, segment):
    name = Path(path).name
    base = re.sub(r"\.(fits|fit|fts)(\.gz)?$", "", name, flags=re.IGNORECASE)
    key = re.sub(r"(^|[._-]){}($|[._-])".format(segment), _remove_segment_token, base, count=1, flags=re.IGNORECASE)
    key = re.sub(r"[._-]+$", "", key)
    key = re.sub(r"^[._-]+", "", key)
    key = re.sub(r"[._-]{2,}", ".", key)
    return key or base


def _remove_segment_token(match):
    prefix = match.group(1)
    suffix = match.group(2)
    if prefix and suffix:
        return prefix
    return ""


def _component_frame_keys(components):
    keys = {}
    for component, paths in components.items():
        keys[component] = {
            _component_frame_key(Path(path), component)
            for path in paths
        }
    return keys


def _component_frame_key(path, component):
    name = Path(path).name
    suffixes = sorted(COMPONENT_SUFFIXES[component], key=len, reverse=True)
    lower_name = name.lower()
    for suffix in suffixes:
        if lower_name.endswith(suffix.lower()):
            return name[:-len(suffix)].rstrip("._-")
    return name


def _frame_key_sort_key(key):
    numbers = [int(match) for match in re.findall(r"\d+", key)]
    return (numbers, key)


def _time_label_from_name(name):
    match = re.search(r"(\d{8})_(\d{6})_TAI", name)
    if not match:
        return None
    date, time = match.groups()
    return "{}-{}-{} {}:{}:{} TAI".format(
        date[:4], date[4:6], date[6:8], time[:2], time[2:4], time[4:6]
    )


def _sequence_sort_key(path):
    text = Path(path).name
    numbers = [int(match) for match in re.findall(r"\d+", text)]
    return (numbers, text)
