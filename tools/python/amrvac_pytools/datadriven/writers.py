"""Writers for AMRVAC data-driven preprocessing products."""

from __future__ import print_function

import json
from pathlib import Path


def ensure_output_dir(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path, data):
    path = Path(path)
    path.write_text(json.dumps(_json_safe(data), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def write_quicklook(path, bx, by, bz, vmax=500.0):
    """Write a simple Bx/By/Bz quicklook plot."""

    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable
    except ImportError as error:
        raise ImportError("matplotlib is required for quicklook plots") from error

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=160, constrained_layout=True)
    for ax, data, title in zip(axes, (bx, by, bz), ("Bx = Bp", "By = -Bt", "Bz = Br")):
        im = ax.imshow(data, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(title)
        ax.set_xlabel("x pixel")
        ax.set_ylabel("y pixel")
        divider = make_axes_locatable(ax)
        color_axis = divider.append_axes("right", size="4%", pad=0.06)
        fig.colorbar(im, cax=color_axis)
    fig.savefig(str(path), bbox_inches="tight")
    plt.close(fig)
    return path


def read_boundary_frame(path):
    """Read one unified magnetic-boundary frame written by this module."""

    import numpy as np

    path = Path(path)
    with path.open("rb") as handle:
        snapshot_time = np.fromfile(handle, dtype=np.float64, count=1)
        nx = np.fromfile(handle, dtype=np.int32, count=1)
        ny = np.fromfile(handle, dtype=np.int32, count=1)
        dx = np.fromfile(handle, dtype=np.float64, count=1)
        dy = np.fromfile(handle, dtype=np.float64, count=1)
        values = np.fromfile(handle, dtype=np.float64)
    if not all(array.size == 1 for array in (snapshot_time, nx, ny, dx, dy)):
        raise ValueError("incomplete magnetic-boundary header: {}".format(path))
    nx_value = int(nx[0])
    ny_value = int(ny[0])
    expected = nx_value * ny_value * 3
    if values.size != expected:
        raise ValueError(
            "magnetic-boundary payload has {} values, expected {}: {}".format(
                values.size, expected, path
            )
        )
    cube = values.reshape((nx_value, ny_value, 3), order="F")
    return {
        "path": str(path),
        "snapshot_time": float(snapshot_time[0]),
        "nx": nx_value,
        "ny": ny_value,
        "dx": float(dx[0]),
        "dy": float(dy[0]),
        "bx": cube[:, :, 0].T,
        "by": cube[:, :, 1].T,
        "bz": cube[:, :, 2].T,
    }


def write_static_boundary_products(output_dir, bx, by, bz, metadata, quicklook=True, vmax=500.0):
    """Write potential/NLFFF Guo-AMRVAC boundary products."""

    from .cartesian import (
        format_guo_amrvac_parameters,
        write_guo_allboundaries,
        write_guo_grid_ini,
        write_guo_nd_ini,
        write_nlfff_boundary,
        write_potential_boundary,
    )

    output_dir = ensure_output_dir(output_dir)
    write_potential_boundary(output_dir / "potential_boundary.dat", bz, metadata)
    write_nlfff_boundary(output_dir / "nlfff_boundary.dat", bx, by, bz, metadata)
    write_guo_allboundaries(output_dir / "allboundaries.dat", bx, by, bz)
    write_guo_grid_ini(output_dir / "grid.ini", nx=metadata.nx, ny=metadata.ny)
    write_guo_nd_ini(output_dir / "nd.ini", nxmax=metadata.nx, nymax=metadata.ny)
    parameters_text = format_guo_amrvac_parameters(metadata)
    (output_dir / "amrvac_parameters.txt").write_text(parameters_text, encoding="utf-8")
    if quicklook:
        write_quicklook(output_dir / "boundary_quicklook.png", bx, by, bz, vmax=vmax)
    return {
        "potential_boundary": str(output_dir / "potential_boundary.dat"),
        "nlfff_boundary": str(output_dir / "nlfff_boundary.dat"),
        "allboundaries": str(output_dir / "allboundaries.dat"),
        "grid_ini": str(output_dir / "grid.ini"),
        "nd_ini": str(output_dir / "nd.ini"),
        "amrvac_parameters": str(output_dir / "amrvac_parameters.txt"),
        "quicklook": str(output_dir / "boundary_quicklook.png") if quicklook else None,
    }


def write_boundary_frame(path, bx, by, bz, snapshot_time, dx, dy):
    """Write one unified magnetic-boundary frame.

    Binary layout:
    snapshot_time, nx, ny, dx, dy, then Bx/By/Bz as an (nx, ny, 3)
    Fortran-ordered stream. The file is intentionally only magnetic data and
    geometry; any driving-time scaling belongs to AMRVAC-side user logic.
    """

    import numpy as np

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    bx = np.asarray(bx, dtype=np.float64)
    by = np.asarray(by, dtype=np.float64)
    bz = np.asarray(bz, dtype=np.float64)
    if bx.shape != by.shape or bx.shape != bz.shape:
        raise ValueError("bx, by, and bz must have identical shapes")
    ny, nx = bx.shape
    data = np.stack((bx, by, bz), axis=-1).transpose(1, 0, 2)
    with path.open("wb") as handle:
        np.asarray(snapshot_time, dtype=np.float64).tofile(handle)
        np.asarray(nx, dtype=np.int32).tofile(handle)
        np.asarray(ny, dtype=np.int32).tofile(handle)
        np.asarray(dx, dtype=np.float64).tofile(handle)
        np.asarray(dy, dtype=np.float64).tofile(handle)
        np.asfortranarray(data, dtype=np.float64).ravel(order="F").tofile(handle)
    return path


def write_boundary_outputs(output_dir, frames, prefix="B"):
    """Write one or more unified magnetic-boundary frames."""

    output_dir = ensure_output_dir(output_dir)
    paths = []
    for index, frame in enumerate(frames, start=1):
        path = output_dir / "{}_{:04d}.dat".format(prefix, index)
        paths.append(
            write_boundary_frame(
                path,
                frame["bx"],
                frame["by"],
                frame["bz"],
                frame["time"],
                frame["dx"],
                frame["dy"],
            )
        )
    return paths
def _json_safe(value):
    try:
        import numpy as np
    except ImportError:
        np = None
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if np is not None and isinstance(value, np.generic):
        return _json_safe(value.item())
    if np is not None and isinstance(value, float) and not np.isfinite(value):
        return None
    if hasattr(value, "as_dict"):
        return _json_safe(value.as_dict)
    return value
