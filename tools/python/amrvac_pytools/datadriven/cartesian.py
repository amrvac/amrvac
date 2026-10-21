"""Helpers for Cartesian simulation grids."""

from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import numpy as np


SOLAR_RADIUS_CM = 6.955e10


class CartesianGrid:
    """A Cartesian simulation grid with explicit axis arrays."""

    def __init__(self, x, y, z, unit=""):
        self.x = x
        self.y = y
        self.z = z
        self.unit = unit

    @property
    def shape(self) -> Tuple[int, int, int]:
        """Return the expected data shape in `(z, y, x)` order."""

        return (len(self.z), len(self.y), len(self.x))

    @property
    def spacing(self) -> Tuple[Optional[float], Optional[float], Optional[float]]:
        """Return mean `(dx, dy, dz)` spacing, or `None` for short axes."""

        return (_mean_spacing(self.x), _mean_spacing(self.y), _mean_spacing(self.z))


def _mean_spacing(axis) -> Optional[float]:
    axis = np.asarray(axis, dtype=float)
    if axis.size < 2:
        return None
    return float(np.nanmean(np.diff(axis)))


class CartesianBoundaryMetadata:
    """Metadata used by Guo/AMRVAC Cartesian bottom-boundary files."""

    def __init__(
        self,
        nx,
        ny,
        nx_physical,
        ny_physical,
        xc_cm,
        yc_cm,
        dx_cm,
        dy_cm,
        xprobmin1,
        xprobmax1,
        xprobmin2,
        xprobmax2,
        xprobmin3,
        xprobmax3,
        boundary_z_10mm=0.0,
        domain_height_10mm=None,
        boundary_plane="first_ghost_center",
    ):
        self.nx = nx
        self.ny = ny
        self.nx_physical = nx_physical
        self.ny_physical = ny_physical
        self.xc_cm = xc_cm
        self.yc_cm = yc_cm
        self.dx_cm = dx_cm
        self.dy_cm = dy_cm
        self.xprobmin1 = xprobmin1
        self.xprobmax1 = xprobmax1
        self.xprobmin2 = xprobmin2
        self.xprobmax2 = xprobmax2
        self.xprobmin3 = xprobmin3
        self.xprobmax3 = xprobmax3
        self.boundary_z_10mm = boundary_z_10mm
        self.domain_height_10mm = (
            xprobmax3 - xprobmin3 if domain_height_10mm is None else domain_height_10mm
        )
        self.boundary_plane = boundary_plane

    @property
    def as_dict(self) -> Dict[str, Union[float, int, str]]:
        return {
            "nx": self.nx,
            "ny": self.ny,
            "nx_physical": self.nx_physical,
            "ny_physical": self.ny_physical,
            "xc_cm": self.xc_cm,
            "yc_cm": self.yc_cm,
            "dx_cm": self.dx_cm,
            "dy_cm": self.dy_cm,
            "xprobmin1": self.xprobmin1,
            "xprobmax1": self.xprobmax1,
            "xprobmin2": self.xprobmin2,
            "xprobmax2": self.xprobmax2,
            "xprobmin3": self.xprobmin3,
            "xprobmax3": self.xprobmax3,
            "boundary_z_10mm": self.boundary_z_10mm,
            "domain_height_10mm": self.domain_height_10mm,
            "boundary_plane": self.boundary_plane,
        }


def local_cartesian_from_heliographic(br, bt, bp):
    """Return local Cartesian `(Bx, By, Bz)` from CEA `Br/Btheta/Bphi`.

    The convention follows the Guo Cartesian examples and the maintained
    `make_cea_patch.py` workflow: x is along CEA longitude, y is northward,
    and z is radially outward.
    """

    br = np.asarray(br, dtype=float)
    bt = np.asarray(bt, dtype=float)
    bp = np.asarray(bp, dtype=float)
    if br.shape != bt.shape or br.shape != bp.shape:
        raise ValueError("br, bt, and bp must have the same shape")
    return bp.copy(), -bt.copy(), br.copy()


def crop_components(bx, by, bz, *, x0: int, y0: int, nx: int, ny: int):
    """Crop same-shaped `(ny, nx)` magnetic components with 0-based pixels."""

    bx = np.asarray(bx, dtype=float)
    by = np.asarray(by, dtype=float)
    bz = np.asarray(bz, dtype=float)
    if bx.shape != by.shape or bx.shape != bz.shape:
        raise ValueError("bx, by, and bz must have the same shape")

    x0 = int(x0)
    y0 = int(y0)
    nx = int(nx)
    ny = int(ny)
    if x0 < 0 or y0 < 0 or nx < 1 or ny < 1:
        raise ValueError("x0, y0 must be non-negative and nx, ny must be positive")
    if x0 + nx > bx.shape[1] or y0 + ny > bx.shape[0]:
        raise ValueError(
            f"crop ({x0}:{x0 + nx}, {y0}:{y0 + ny}) exceeds input shape {bx.shape}"
        )
    region = (slice(y0, y0 + ny), slice(x0, x0 + nx))
    return bx[region].copy(), by[region].copy(), bz[region].copy()


def crop_components_with_padding(
    bx,
    by,
    bz,
    *,
    x0: int,
    y0: int,
    nx: int,
    ny: int,
    pad_x: int,
    pad_y: Optional[int] = None,
    fill_value: float = 0.0,
):
    """Crop `(x0, y0, nx, ny)` plus padding, filling out-of-frame pixels.

    The requested `x0, y0, nx, ny` describe the physical/potential region.
    `pad_x/pad_y` pixels are added on each side before cropping. This is useful
    when Guo/AMRVAC needs two ghost cells after multigrid rebinning: for
    `level=3`, pass `pad_x=pad_y=2*4`.
    """

    bx = np.asarray(bx, dtype=float)
    by = np.asarray(by, dtype=float)
    bz = np.asarray(bz, dtype=float)
    if bx.shape != by.shape or bx.shape != bz.shape:
        raise ValueError("bx, by, and bz must have the same shape")
    if bx.ndim != 2:
        raise ValueError("magnetic components must be 2-D")

    x0 = int(x0)
    y0 = int(y0)
    nx = int(nx)
    ny = int(ny)
    pad_x = int(pad_x)
    pad_y = int(pad_x if pad_y is None else pad_y)
    if x0 < 0 or y0 < 0 or nx < 1 or ny < 1:
        raise ValueError("x0, y0 must be non-negative and nx, ny must be positive")
    if pad_x < 0 or pad_y < 0:
        raise ValueError("pad_x and pad_y must be non-negative")
    if x0 + nx > bx.shape[1] or y0 + ny > bx.shape[0]:
        raise ValueError(
            f"physical window ({x0}:{x0 + nx}, {y0}:{y0 + ny}) exceeds input shape {bx.shape}"
        )

    requested_x0 = x0 - pad_x
    requested_y0 = y0 - pad_y
    requested_x1 = x0 + nx + pad_x
    requested_y1 = y0 + ny + pad_y
    out_shape = (requested_y1 - requested_y0, requested_x1 - requested_x0)

    source_x0 = max(requested_x0, 0)
    source_y0 = max(requested_y0, 0)
    source_x1 = min(requested_x1, bx.shape[1])
    source_y1 = min(requested_y1, bx.shape[0])
    dest_x0 = source_x0 - requested_x0
    dest_y0 = source_y0 - requested_y0
    dest_x1 = dest_x0 + (source_x1 - source_x0)
    dest_y1 = dest_y0 + (source_y1 - source_y0)

    outputs = []
    for component in (bx, by, bz):
        out = np.full(out_shape, fill_value, dtype=float)
        out[dest_y0:dest_y1, dest_x0:dest_x1] = component[source_y0:source_y1, source_x0:source_x1]
        outputs.append(out)

    padding = {
        "requested_x0": requested_x0,
        "requested_y0": requested_y0,
        "requested_x1": requested_x1,
        "requested_y1": requested_y1,
        "source_x0": source_x0,
        "source_y0": source_y0,
        "source_x1": source_x1,
        "source_y1": source_y1,
        "left": max(0, -requested_x0),
        "right": max(0, requested_x1 - bx.shape[1]),
        "bottom": max(0, -requested_y0),
        "top": max(0, requested_y1 - bx.shape[0]),
        "fill_value": float(fill_value),
    }
    return (*outputs, padding)


def multigrid_reduce_components(bx, by, bz, *, level: int):
    """Apply Guo's modified multigrid rebinning for one requested level.

    `level=1` keeps the original grid. `level=2` performs 2x2 block averaging,
    `level=3` performs 4x4 block averaging, and so on.
    """

    reducer = multigrid_reducer(level)
    bx = np.asarray(bx, dtype=float)
    by = np.asarray(by, dtype=float)
    bz = np.asarray(bz, dtype=float)
    if bx.shape != by.shape or bx.shape != bz.shape:
        raise ValueError("bx, by, and bz must have the same shape")
    if bx.ndim != 2:
        raise ValueError("magnetic components must be 2-D")
    if bx.shape[0] % reducer or bx.shape[1] % reducer:
        raise ValueError(
            f"input shape {bx.shape} must be divisible by reducer {reducer}; "
            "choose nx, ny as multiples of 2^(level-1)"
        )
    if reducer == 1:
        return bx.copy(), by.copy(), bz.copy()
    return (
        _block_mean_2d(bx, reducer, reducer),
        _block_mean_2d(by, reducer, reducer),
        _block_mean_2d(bz, reducer, reducer),
    )


def multigrid_reducer(level: int) -> int:
    level = int(level)
    if level < 1:
        raise ValueError("level must be >= 1")
    return 2 ** (level - 1)


def _block_mean_2d(data, factor_y: int, factor_x: int):
    data = np.asarray(data, dtype=float)
    ny = data.shape[0] // factor_y
    nx = data.shape[1] // factor_x
    blocks = data.reshape(ny, factor_y, nx, factor_x)
    return np.nanmean(blocks, axis=(1, 3))


def preprocess_cartesian_field(
    bx,
    by,
    bz,
    *,
    mu3: float = 0.001,
    mu4: float = 0.01,
    max_iter: int = 2000,
    tol: float = 1.0e-4,
    report_interval: int = 50,
    dx: float = 1.0,
    dy: float = 1.0,
    geometry_mode: str = "legacy_unit_square",
    edge_treatment: str = "legacy_periodic",
    progress=None,
):
    """Wiegelmann-style preprocessing used by Guo's `prepro_wie/prepro.pro`.

    The implementation mirrors the IDL code's normalized force, torque, data,
    and smoothing terms.  ``edge_treatment='legacy_periodic'`` preserves the
    old IDL-style periodic Laplacian.  Notebook policy uses
    ``edge_treatment='nonperiodic'`` for observational magnetograms.
    """

    bx = np.asarray(bx, dtype=float).copy()
    by = np.asarray(by, dtype=float).copy()
    bz = np.asarray(bz, dtype=float).copy()
    bx_input = bx.copy()
    by_input = by.copy()
    bz_input = bz.copy()
    if bx.shape != by.shape or bx.shape != bz.shape:
        raise ValueError("bx, by, and bz must have the same shape")
    if bx.ndim != 2:
        raise ValueError("magnetic components must be 2-D")
    if not all(np.all(np.isfinite(item)) for item in (bx, by, bz)):
        raise ValueError("preprocessing requires finite magnetic components")

    ny, nx = bz.shape
    if nx < 2 or ny < 2:
        raise ValueError("preprocessing requires at least 2 pixels in each direction")
    max_iter = int(max_iter)
    report_interval = int(report_interval)
    if max_iter < 0:
        raise ValueError("max_iter must be non-negative")
    if tol < 0.0:
        raise ValueError("tol must be non-negative")
    dx = float(dx)
    dy = float(dy)
    if dx <= 0.0 or dy <= 0.0:
        raise ValueError("dx and dy must be positive")
    edge_treatment = str(edge_treatment).strip().lower()
    if edge_treatment not in ("legacy_periodic", "periodic", "nonperiodic"):
        raise ValueError("edge_treatment must be legacy_periodic, periodic, or nonperiodic")

    x, y, dxdy = _preprocess_geometry(nx, ny, dx, dy, geometry_mode)

    bave2d = np.nanmean(np.sqrt(bx**2 + by**2 + bz**2))
    if not np.isfinite(bave2d) or bave2d == 0.0:
        raise ValueError("cannot preprocess a zero or non-finite magnetic field")

    bx /= bave2d
    by /= bave2d
    bz /= bave2d
    bxo = bx.copy()
    byo = by.copy()
    bzo = bz.copy()

    mu1 = 0.1
    mu2 = 0.1
    mu3_input = float(mu3)
    mu4_input = float(mu4)
    if mu3_input < 0.0 or mu4_input < 0.0:
        raise ValueError("mu3 and mu4 must be non-negative")
    mu3 = mu3_input / 10.0
    mu4 = mu4_input / 10.0

    emag = np.nansum(bx**2 + by**2 + bz**2)
    r = np.sqrt(x**2 + y**2)
    ihelp = np.nansum(r * (bx**2 + by**2 + bz**2))
    if emag == 0.0 or ihelp == 0.0:
        raise ValueError("cannot preprocess a zero magnetic field")

    metrics, terms = _preprocess_metrics(
        bx, by, bz, bxo, byo, bzo, x, y, dxdy, emag, ihelp,
        edge_treatment=edge_treatment,
    )
    metrics.update({
        "iteration": 0,
        "iterations": 0,
        "dL": float(np.inf),
        "converged": False,
        "stop_reason": "max_iter_zero" if max_iter == 0 else "not_started",
        "effective_mu1": float(mu1),
        "effective_mu2": float(mu2),
        "effective_mu3": float(mu3),
        "effective_mu4": float(mu4),
        "input_mu3": float(mu3_input),
        "input_mu4": float(mu4_input),
        "tol": float(tol),
        "max_iter": int(max_iter),
        "geometry_mode": str(geometry_mode),
        "edge_treatment": edge_treatment,
        "dx": float(dx),
        "dy": float(dy),
    })
    if progress is not None and report_interval > 0:
        progress(metrics)
    if max_iter == 0:
        return bx_input, by_input, bz_input, metrics

    old_l12 = metrics["L1"] + metrics["L2"]
    old_l3 = metrics["L3"]
    old_l4 = metrics["L4"]
    stop_reason = "max_iterations"

    for iteration in range(1, max_iter + 1):
        term1a, term1b, term1c, term2a, term2b, term2c, term4a, term4b, term4c = terms
        bx_next = (
            bx
            + mu1 * (-2.0 * term1a * bz + 4.0 * term1c * bx)
            - mu3 * 2.0 * (bx - bxo)
            - mu4 * term4a
            + mu2 * (4.0 * term2a * x * bx + 4.0 * term2b * y * bx - 2.0 * term2c * y * bz)
        )
        by_next = (
            by
            + mu1 * (-2.0 * term1b * bz + 4.0 * term1c * by)
            - mu3 * 2.0 * (by - byo)
            - mu4 * term4b
            + mu2 * (4.0 * term2a * x * by + 4.0 * term2b * y * by + 2.0 * term2c * x * bz)
        )
        bz_next = bz - mu3 * 2.0 * (bz - bzo) - mu4 * term4c
        bx, by, bz = bx_next, by_next, bz_next
        if not all(np.all(np.isfinite(item)) for item in (bx, by, bz)):
            raise FloatingPointError("preprocessing produced non-finite values at iteration {}".format(iteration))

        metrics, terms = _preprocess_metrics(
            bx, by, bz, bxo, byo, bzo, x, y, dxdy, emag, ihelp,
            edge_treatment=edge_treatment,
        )
        l12 = metrics["L1"] + metrics["L2"]
        l3 = metrics["L3"]
        l4 = metrics["L4"]
        dl = (
            abs(l12 - old_l12) / max(abs(l12), np.finfo(float).tiny)
            + abs(l3 - old_l3) / max(abs(l3), np.finfo(float).tiny)
            + abs(l4 - old_l4) / max(abs(l4), np.finfo(float).tiny)
        )
        metrics.update({
            "iteration": int(iteration),
            "iterations": int(iteration),
            "dL": float(dl),
            "converged": bool(dl <= tol),
            "stop_reason": "converged" if dl <= tol else "max_iterations",
            "effective_mu1": float(mu1),
            "effective_mu2": float(mu2),
            "effective_mu3": float(mu3),
            "effective_mu4": float(mu4),
            "input_mu3": float(mu3_input),
            "input_mu4": float(mu4_input),
            "tol": float(tol),
            "max_iter": int(max_iter),
            "geometry_mode": str(geometry_mode),
            "edge_treatment": edge_treatment,
            "dx": float(dx),
            "dy": float(dy),
        })
        if progress is not None and report_interval > 0 and iteration % report_interval == 0:
            progress(metrics)
        if dl <= tol:
            stop_reason = "converged"
            break
        old_l12, old_l3, old_l4 = l12, l3, l4
    metrics["stop_reason"] = stop_reason
    metrics["converged"] = bool(stop_reason == "converged")

    return bx * bave2d, by * bave2d, bz * bave2d, metrics


def _preprocess_geometry(nx, ny, dx, dy, geometry_mode):
    mode = str(geometry_mode).strip().lower()
    if mode == "legacy_unit_square":
        x1 = np.linspace(0.0, 1.0, nx)
        y1 = np.linspace(0.0, 1.0, ny)
        dxdy = (1.0 / (nx - 1)) * (1.0 / (ny - 1))
    elif mode == "centered":
        x1 = (np.arange(nx, dtype=float) - 0.5 * (nx - 1)) * float(dx)
        y1 = (np.arange(ny, dtype=float) - 0.5 * (ny - 1)) * float(dy)
        dxdy = abs(float(dx) * float(dy))
    else:
        raise ValueError("geometry_mode must be legacy_unit_square or centered")
    x, y = np.meshgrid(x1, y1)
    return x, y, float(dxdy)


def _laplace(data, edge_treatment):
    edge_treatment = str(edge_treatment).strip().lower()
    if edge_treatment in ("legacy_periodic", "periodic"):
        return _periodic_laplace(data)
    if edge_treatment != "nonperiodic":
        raise ValueError("edge_treatment must be legacy_periodic, periodic, or nonperiodic")
    padded = np.pad(data, 1, mode="edge")
    return (
        -4.0 * data
        + padded[1:-1, 2:]
        + padded[1:-1, :-2]
        + padded[2:, 1:-1]
        + padded[:-2, 1:-1]
    )


def _preprocess_metrics(
    bx, by, bz, bxo, byo, bzo, x, y, dxdy, emag, ihelp,
    edge_treatment="legacy_periodic",
):
    term1a = np.nansum(bx * bz) / emag
    term1b = np.nansum(by * bz) / emag
    term1c = (np.nansum(bz**2) - np.nansum(bx**2 + by**2)) / emag

    term2a = (np.nansum(x * bz**2) - np.nansum(x * (bx**2 + by**2))) / ihelp
    term2b = (np.nansum(y * bz**2) - np.nansum(y * (bx**2 + by**2))) / ihelp
    term2c = (np.nansum(y * bx * bz) - np.nansum(x * by * bz)) / ihelp

    lap_bx = _laplace(bx, edge_treatment)
    lap_by = _laplace(by, edge_treatment)
    lap_bz = _laplace(bz, edge_treatment)
    term4a = 2.0 * _laplace(lap_bx, edge_treatment)
    term4b = 2.0 * _laplace(lap_by, edge_treatment)
    term4c = 2.0 * _laplace(lap_bz, edge_treatment)

    eps_force = abs(term1a) + abs(term1b) + abs(term1c)
    eps_torque = abs(term2a) + abs(term2b) + abs(term2c)
    eps_smooth = (
        np.nansum(np.abs(lap_bx)) + np.nansum(np.abs(lap_by)) + np.nansum(np.abs(lap_bz))
    ) * dxdy / np.sqrt(emag)
    metrics = {
        "L1": float(term1a**2 + term1b**2 + term1c**2),
        "L2": float(term2a**2 + term2b**2 + term2c**2),
        "L3": float(
            np.nansum((bx - bxo) ** 2 + (by - byo) ** 2 + (bz - bzo) ** 2)
            / (emag * emag)
        ),
        "L4": float(dxdy * np.nansum(lap_bx**2 + lap_by**2 + lap_bz**2) / emag),
        "eps_force": float(eps_force),
        "eps_torque": float(eps_torque),
        "eps_smooth": float(eps_smooth),
    }
    terms = (term1a, term1b, term1c, term2a, term2b, term2c, term4a, term4b, term4c)
    return metrics, terms


def cartesian_preprocess_diagnostics(
    bx, by, bz, dx=1.0, dy=1.0,
    geometry_mode="legacy_unit_square", edge_treatment="legacy_periodic",
):
    """Return force, torque, and smoothing diagnostics for a Cartesian boundary."""

    bx = np.asarray(bx, dtype=float)
    by = np.asarray(by, dtype=float)
    bz = np.asarray(bz, dtype=float)
    if bx.shape != by.shape or bx.shape != bz.shape:
        raise ValueError("bx, by, and bz must have the same shape")
    if bx.ndim != 2:
        raise ValueError("magnetic components must be 2-D")
    if not all(np.all(np.isfinite(item)) for item in (bx, by, bz)):
        raise ValueError("diagnostics require finite magnetic components")

    ny, nx = bz.shape
    if nx < 2 or ny < 2:
        raise ValueError("diagnostics require at least 2 pixels in each direction")
    dx = float(dx)
    dy = float(dy)
    if dx <= 0.0 or dy <= 0.0:
        raise ValueError("dx and dy must be positive")

    bave2d = np.nanmean(np.sqrt(bx**2 + by**2 + bz**2))
    if not np.isfinite(bave2d) or bave2d == 0.0:
        raise ValueError("cannot diagnose a zero or non-finite magnetic field")
    bx = bx / bave2d
    by = by / bave2d
    bz = bz / bave2d

    x, y, dxdy = _preprocess_geometry(nx, ny, dx, dy, geometry_mode)

    emag = np.nansum(bx**2 + by**2 + bz**2)
    r = np.sqrt(x**2 + y**2)
    ihelp = np.nansum(r * (bx**2 + by**2 + bz**2))
    if emag == 0.0 or ihelp == 0.0:
        raise ValueError("cannot diagnose a zero magnetic field")

    term1a = np.nansum(bx * bz) / emag
    term1b = np.nansum(by * bz) / emag
    term1c = (np.nansum(bz**2) - np.nansum(bx**2 + by**2)) / emag
    term2a = (np.nansum(x * bz**2) - np.nansum(x * (bx**2 + by**2))) / ihelp
    term2b = (np.nansum(y * bz**2) - np.nansum(y * (bx**2 + by**2))) / ihelp
    term2c = (np.nansum(y * bx * bz) - np.nansum(x * by * bz)) / ihelp

    lap_bx = _laplace(bx, edge_treatment)
    lap_by = _laplace(by, edge_treatment)
    lap_bz = _laplace(bz, edge_treatment)
    eps_smooth = (
        np.nansum(np.abs(lap_bx)) + np.nansum(np.abs(lap_by)) + np.nansum(np.abs(lap_bz))
    ) * dxdy / np.sqrt(emag)

    return {
        "eps_force": float(abs(term1a) + abs(term1b) + abs(term1c)),
        "eps_torque": float(abs(term2a) + abs(term2b) + abs(term2c)),
        "eps_smooth": float(eps_smooth),
        "L1": float(term1a**2 + term1b**2 + term1c**2),
        "L2": float(term2a**2 + term2b**2 + term2c**2),
        "L4": float(dxdy * np.nansum(lap_bx**2 + lap_by**2 + lap_bz**2) / emag),
        "mean_field_strength": float(bave2d),
        "Bz_flux": float(np.nansum(bz) * dxdy),
        "geometry_mode": str(geometry_mode),
        "edge_treatment": str(edge_treatment),
        "dx": float(dx),
        "dy": float(dy),
    }


def _periodic_laplace(data):
    return (
        -4.0 * data
        + np.roll(data, 1, axis=1)
        + np.roll(data, -1, axis=1)
        + np.roll(data, 1, axis=0)
        + np.roll(data, -1, axis=0)
    )


def boundary_metadata(
    *,
    nx: int,
    ny: int,
    xc_cm: float,
    yc_cm: float,
    dx_cm: float,
    dy_cm: float,
    nghost: int = 2,
    z_min_10mm: float = 0.0,
) -> CartesianBoundaryMetadata:
    """Return Guo/AMRVAC Cartesian domain parameters.

    Boundary maps include `nghost` pixels on each horizontal side. The physical
    domain is therefore `(nx - 2*nghost, ny - 2*nghost)`. Vertically, the
    magnetogram is at `z_min_10mm` on the first lower ghost-cell center; the
    staging step places the physical lower face half a finest cell above it.
    """

    nx = int(nx)
    ny = int(ny)
    nghost = int(nghost)
    nx_physical = nx - 2 * nghost
    ny_physical = ny - 2 * nghost
    if nx_physical <= 0 or ny_physical <= 0:
        raise ValueError("nx and ny must be larger than twice nghost")

    x1 = float(xc_cm) - nx_physical * float(dx_cm) / 2.0
    x2 = float(xc_cm) + nx_physical * float(dx_cm) / 2.0
    y1 = float(yc_cm) - ny_physical * float(dy_cm) / 2.0
    y2 = float(yc_cm) + ny_physical * float(dy_cm) / 2.0
    return CartesianBoundaryMetadata(
        nx=nx,
        ny=ny,
        nx_physical=nx_physical,
        ny_physical=ny_physical,
        xc_cm=float(xc_cm),
        yc_cm=float(yc_cm),
        dx_cm=float(dx_cm),
        dy_cm=float(dy_cm),
        xprobmin1=x1 * 1.0e-9,
        xprobmax1=x2 * 1.0e-9,
        xprobmin2=y1 * 1.0e-9,
        xprobmax2=y2 * 1.0e-9,
        xprobmin3=float(z_min_10mm),
        xprobmax3=float(z_min_10mm)+(y2 - y1) * 1.0e-9,
        boundary_z_10mm=float(z_min_10mm),
        domain_height_10mm=(y2 - y1) * 1.0e-9,
        boundary_plane="first_ghost_center",
    )


def write_guo_allboundaries(path, bx, by, bz):
    """Write Guo `allboundaries.dat` ASCII triples with x fastest."""

    path = Path(path)
    bx = np.asarray(bx, dtype=float)
    by = np.asarray(by, dtype=float)
    bz = np.asarray(bz, dtype=float)
    if bx.shape != by.shape or bx.shape != bz.shape:
        raise ValueError("bx, by, and bz must have the same shape")
    with path.open("w", encoding="utf-8") as handle:
        for iy in range(bx.shape[0]):
            for ix in range(bx.shape[1]):
                handle.write(f"{bx[iy, ix]:.9e}\n")
                handle.write(f"{by[iy, ix]:.9e}\n")
                handle.write(f"{bz[iy, ix]:.9e}\n")
    return path


def write_guo_grid_ini(path, *, nx: int, ny: int, nz: Optional[int] = None, nd: Optional[int] = None, mu: float = 0.1):
    """Write the `grid.ini` text file produced by Guo's multigrid helper."""

    path = Path(path)
    nx = int(nx)
    ny = int(ny)
    nz = int(min(nx, ny) if nz is None else nz)
    nd = int(min(nx, ny) // 8 if nd is None else nd)
    with path.open("w", encoding="utf-8") as handle:
        for key, value in (("nx", nx), ("ny", ny), ("nz", nz), ("mu", float(mu)), ("nd", nd)):
            handle.write(f"{key}\n")
            handle.write(f"{value}\n")
    return path


def write_guo_nd_ini(
    path,
    *,
    nxmax: int,
    nymax: int,
    nzmax: Optional[int] = None,
    nd: Optional[int] = None,
):
    """Write the `nd.ini` text file produced before multigrid reduction."""

    path = Path(path)
    nxmax = int(nxmax)
    nymax = int(nymax)
    nzmax = int(min(nxmax, nymax) if nzmax is None else nzmax)
    nd = int(min(nxmax, nymax) // 8 if nd is None else nd)
    with path.open("w", encoding="utf-8") as handle:
        for key, value in (("nd", nd), ("nxmax", nxmax), ("nymax", nymax), ("nzmax", nzmax)):
            handle.write(f"{key}\n")
            handle.write(f"{value}\n")
    return path


def write_potential_boundary(path, bz, metadata: CartesianBoundaryMetadata):
    """Write `potential_boundary.dat` expected by Guo's potential example."""

    bz = np.asarray(bz, dtype=float)
    if bz.ndim != 2:
        raise ValueError("bz must be 2-D")
    if bz.shape[0] <= 4 or bz.shape[1] <= 4:
        raise ValueError("potential boundary requires at least 5 pixels in each direction")
    core = bz[2:-2, 2:-2]
    if core.shape != (metadata.ny_physical, metadata.nx_physical):
        raise ValueError(
            f"metadata physical shape {(metadata.ny_physical, metadata.nx_physical)} "
            f"does not match cropped Bz shape {core.shape}"
        )
    path = Path(path)
    with path.open("wb") as handle:
        np.asarray([metadata.nx_physical, metadata.ny_physical], dtype=np.int32).tofile(handle)
        np.asarray([metadata.xc_cm, metadata.yc_cm, metadata.dx_cm, metadata.dy_cm], dtype=np.float64).tofile(handle)
        np.asarray(core, dtype=np.float64).ravel(order="C").tofile(handle)
    return path


def write_nlfff_boundary(path, bx, by, bz, metadata: CartesianBoundaryMetadata):
    """Write `nlfff_boundary.dat` expected by Guo's magnetofriction example."""

    bx = np.asarray(bx, dtype=float)
    by = np.asarray(by, dtype=float)
    bz = np.asarray(bz, dtype=float)
    if bx.shape != by.shape or bx.shape != bz.shape:
        raise ValueError("bx, by, and bz must have the same shape")
    if bx.shape != (metadata.ny, metadata.nx):
        raise ValueError(f"metadata shape {(metadata.ny, metadata.nx)} does not match field shape {bx.shape}")
    path = Path(path)
    with path.open("wb") as handle:
        np.asarray([metadata.nx, metadata.ny], dtype=np.int32).tofile(handle)
        np.asarray([metadata.xc_cm, metadata.yc_cm, metadata.dx_cm, metadata.dy_cm], dtype=np.float64).tofile(handle)
        for component in (bx, by, bz):
            np.asarray(component, dtype=np.float64).ravel(order="C").tofile(handle)
    return path


def format_guo_amrvac_parameters(
    metadata: CartesianBoundaryMetadata,
    *,
    domain_nx3: Optional[int] = None,
    block_nx1: Optional[int] = None,
    block_nx2: Optional[int] = None,
    block_nx3: Optional[int] = None,
):
    """Return an `amrvac.par` meshlist snippet for the Guo Cartesian examples."""

    domain_nx3 = metadata.ny_physical if domain_nx3 is None else int(domain_nx3)
    lines = []
    if block_nx1 is not None:
        lines.append(f"        block_nx1={int(block_nx1)}")
    if block_nx2 is not None:
        lines.append(f"        block_nx2={int(block_nx2)}")
    if block_nx3 is not None:
        lines.append(f"        block_nx3={int(block_nx3)}")
    lines.extend(
        [
            f"        domain_nx1={metadata.nx_physical}",
            f"        domain_nx2={metadata.ny_physical}",
            f"        domain_nx3={domain_nx3}",
            f"        xprobmin1={metadata.xprobmin1:.12g}d0",
            f"        xprobmax1={metadata.xprobmax1:.12g}d0",
            f"        xprobmin2={metadata.xprobmin2:.12g}d0",
            f"        xprobmax2={metadata.xprobmax2:.12g}d0",
            f"        xprobmin3={metadata.xprobmin3:.12g}d0",
            f"        xprobmax3={metadata.xprobmax3:.12g}d0",
        ]
    )
    return "\n".join(lines) + "\n"
