"""Helpers for SHARP-like Cylindrical Equal Area magnetic patches."""

from typing import Tuple

import numpy as np


DEGREE = np.pi / 180.0


class CeaPatchGrid:
    """A local CEA patch grid centered on heliographic longitude/latitude."""

    def __init__(self, x, y, lon, lat, center_lon, center_lat, resolution_degree):
        self.x = x
        self.y = y
        self.lon = lon
        self.lat = lat
        self.center_lon = center_lon
        self.center_lat = center_lat
        self.resolution_degree = resolution_degree

    @property
    def shape(self) -> Tuple[int, int]:
        """Return the patch shape in `(lat/CEA-y, lon/CEA-x)` order."""

        return self.lon.shape


def cea_xy_centers(nx: int, ny: int, resolution_degree: float):
    """Return centered CEA `x` and `y` coordinate arrays in degrees."""

    nx = int(nx)
    ny = int(ny)
    resolution_degree = float(resolution_degree)
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be positive")
    if resolution_degree <= 0.0:
        raise ValueError("resolution_degree must be positive")

    x = (np.arange(nx, dtype=float) - 0.5 * (nx - 1)) * resolution_degree
    y = (np.arange(ny, dtype=float) - 0.5 * (ny - 1)) * resolution_degree
    return x, y


def cea_lon_lat_from_xy(x_degree, y_degree, center_lon_degree: float, center_lat_degree: float):
    """Convert local CEA coordinates to heliographic lon/lat in degrees.

    This implements the oblique CEA inverse mapping used for SHARP CEA patches
    as described by Sun (2018). `x_degree` and `y_degree` are local patch
    coordinates centered on `(center_lon_degree, center_lat_degree)`.
    """

    x = np.radians(np.asarray(x_degree, dtype=float))
    y = np.radians(np.asarray(y_degree, dtype=float))
    lon0 = float(center_lon_degree) * DEGREE
    lat0 = float(center_lat_degree) * DEGREE

    if np.any(np.abs(y) > 1.0):
        raise ValueError("CEA y coordinate is outside the valid sine-latitude range")

    root = np.sqrt(np.maximum(1.0 - y**2, 0.0))
    sin_lat = np.cos(lat0) * y + np.sin(lat0) * root * np.cos(x)
    lat = np.arcsin(np.clip(sin_lat, -1.0, 1.0))
    lon_offset = np.arctan2(
        root * np.sin(x),
        np.cos(lat0) * root * np.cos(x) - np.sin(lat0) * y,
    )
    lon = lon0 + lon_offset
    return np.degrees(lon), np.degrees(lat)


def make_cea_patch_grid(
    center_lon: float,
    center_lat: float,
    width_degree: float,
    height_degree: float,
    resolution_degree: float,
) -> CeaPatchGrid:
    """Return a centered oblique-CEA patch grid."""

    width_degree = float(width_degree)
    height_degree = float(height_degree)
    resolution_degree = float(resolution_degree)
    if width_degree <= 0.0 or height_degree <= 0.0:
        raise ValueError("width_degree and height_degree must be positive")

    nx = int(np.round(width_degree / resolution_degree))
    ny = int(np.round(height_degree / resolution_degree))
    if nx < 1 or ny < 1:
        raise ValueError("width_degree and height_degree must span at least one pixel")

    x, y = cea_xy_centers(nx, ny, resolution_degree)
    x2, y2 = np.meshgrid(x, y)
    lon, lat = cea_lon_lat_from_xy(x2, y2, center_lon, center_lat)
    return CeaPatchGrid(
        x=x,
        y=y,
        lon=lon,
        lat=lat,
        center_lon=float(center_lon),
        center_lat=float(center_lat),
        resolution_degree=resolution_degree,
    )


def gaussian_block_reduce(
    data,
    factor_y: int,
    factor_x: int,
    *,
    sigma_y: float,
    sigma_x: float,
    truncate: float = 2.0,
    cval: float = np.nan,
):
    """Smooth and reduce a 2-D array by integer block factors.

    This is a lightweight approximation of the SHARP CEA oversample-and-resize
    step: finite values and their validity mask are smoothed separately, then
    each integer block is averaged with the smoothed validity weights.
    """

    try:
        from scipy.ndimage import gaussian_filter
    except ImportError as error:
        raise ImportError("scipy is required for SHARP-like Gaussian reduction") from error

    data = np.asarray(data, dtype=float)
    if data.ndim != 2:
        raise ValueError("data must be a 2-D array")
    factor_y = int(factor_y)
    factor_x = int(factor_x)
    if factor_y < 1 or factor_x < 1:
        raise ValueError("reduction factors must be positive")
    if data.shape[0] % factor_y != 0 or data.shape[1] % factor_x != 0:
        raise ValueError("data shape must be divisible by the reduction factors")

    valid = np.isfinite(data).astype(float)
    values = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
    sigma = (float(sigma_y), float(sigma_x))
    smoothed_values = gaussian_filter(values * valid, sigma=sigma, mode="constant", cval=0.0, truncate=truncate)
    smoothed_weight = gaussian_filter(valid, sigma=sigma, mode="constant", cval=0.0, truncate=truncate)

    with np.errstate(invalid="ignore", divide="ignore"):
        smoothed = smoothed_values / smoothed_weight

    ny = data.shape[0] // factor_y
    nx = data.shape[1] // factor_x
    smoothed = smoothed.reshape(ny, factor_y, nx, factor_x)
    weights = smoothed_weight.reshape(ny, factor_y, nx, factor_x)

    weighted_sum = np.nansum(smoothed * weights, axis=(1, 3))
    weight_sum = np.sum(weights, axis=(1, 3))
    with np.errstate(invalid="ignore", divide="ignore"):
        reduced = weighted_sum / weight_sum
    reduced[weight_sum <= 0.0] = cval
    return reduced


def _source_pixels_from_heliographic(smap, lon, lat, SkyCoord=None, u=None, frames=None):
    """Map heliographic lon/lat samples to source-map ``(y, x)`` pixels."""

    if SkyCoord is None:
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        from sunpy.coordinates import frames

    try:
        obstime = smap.date
    except Exception:
        obstime = None
    coords = SkyCoord(
        np.asarray(lon).ravel() * u.deg,
        np.asarray(lat).ravel() * u.deg,
        frame=frames.HeliographicStonyhurst,
        obstime=obstime,
        rsun=smap.rsun_meters,
        observer="earth",
    )
    hpc = coords.transform_to(smap.coordinate_frame)
    pixel_x, pixel_y = smap.world_to_pixel(hpc)
    shape = np.asarray(lon).shape
    return np.vstack([
        pixel_y.value.reshape(shape).ravel(),
        pixel_x.value.reshape(shape).ravel(),
    ])


def _sample_at_source_pixels(
    data,
    pixel_yx,
    shape,
    order,
    cval,
    map_coordinates=None,
):
    """Interpolate finite source pixels and preserve invalid/outside samples."""

    if map_coordinates is None:
        try:
            from scipy.ndimage import map_coordinates
        except ImportError as error:
            raise ImportError("scipy is required for CEA image sampling") from error

    coords = np.asarray(pixel_yx, dtype=float).copy()
    invalid = ~np.isfinite(coords).all(axis=0)
    coords[:, invalid] = -1.0
    data = np.asarray(data, dtype=float)
    finite = np.isfinite(data)
    clean = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
    sampled = map_coordinates(clean, coords, order=order, mode="constant", cval=0.0)
    weight = map_coordinates(
        finite.astype(float), coords, order=0, mode="constant", cval=0.0
    )
    sampled[invalid | (weight < 0.5)] = cval
    return sampled.reshape(shape)


def hmi_native_vector_components(field, inclination, azimuth):
    """Convert HMI `B, inclination, azimuth` to native `(Bxi, Beta, Bzeta)`.

    The convention follows the SHARP/JSOC vector-field azimuth basis: `+xi`
    is CCD +x, `+eta` is CCD +y, and `+zeta` is out of the image plane /
    line-of-sight.
    """

    field = np.asarray(field, dtype=float)
    inclination = np.asarray(inclination, dtype=float)
    azimuth = np.asarray(azimuth, dtype=float)
    if field.shape != inclination.shape or field.shape != azimuth.shape:
        raise ValueError("field, inclination, and azimuth must have the same shape")

    gamma = np.radians(inclination)
    psi = np.radians(azimuth)
    bxi = field * np.sin(gamma) * np.sin(psi)
    beta = -field * np.sin(gamma) * np.cos(psi)
    bzeta = field * np.cos(gamma)
    return bxi, beta, bzeta


def native_to_heliographic_components(
    bxi,
    beta,
    bzeta,
    lon_degree,
    lat_degree,
    *,
    disk_lon_degree: float = 0.0,
    disk_lat_degree: float = 0.0,
    p_angle_degree: float = 0.0,
):
    """Transform native image-plane vectors to `(Br, Btheta, Bphi)`.

    `Btheta` is positive southward, matching the SHARP CEA convention.
    Longitudes should use the same origin as `disk_lon_degree`; for
    Stonyhurst coordinates this is usually zero.
    """

    bxi = np.asarray(bxi, dtype=float)
    beta = np.asarray(beta, dtype=float)
    bzeta = np.asarray(bzeta, dtype=float)
    lon = np.radians(np.asarray(lon_degree, dtype=float))
    lat = np.radians(np.asarray(lat_degree, dtype=float))
    if bxi.shape != beta.shape or bxi.shape != bzeta.shape:
        raise ValueError("bxi, beta, and bzeta must have the same shape")
    if bxi.shape != lon.shape or bxi.shape != lat.shape:
        raise ValueError("field components and lon/lat must have the same shape")

    lon0 = float(disk_lon_degree) * DEGREE
    disk_lat = float(disk_lat_degree) * DEGREE
    p_angle = float(p_angle_degree) * DEGREE
    dlon = lon - lon0

    sin_b = np.sin(disk_lat)
    cos_b = np.cos(disk_lat)
    sin_p = np.sin(p_angle)
    cos_p = np.cos(p_angle)
    sin_lon = np.sin(dlon)
    cos_lon = np.cos(dlon)
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)

    k11 = cos_lat * (sin_b * sin_p * cos_lon + cos_p * sin_lon) - sin_lat * (cos_b * sin_p)
    k12 = -cos_lat * (sin_b * cos_p * cos_lon - sin_p * sin_lon) + sin_lat * (cos_b * cos_p)
    k13 = cos_lat * cos_b * cos_lon + sin_lat * sin_b

    k21 = sin_lat * (sin_b * sin_p * cos_lon + cos_p * sin_lon) + cos_lat * (cos_b * sin_p)
    k22 = -sin_lat * (sin_b * cos_p * cos_lon - sin_p * sin_lon) - cos_lat * (cos_b * cos_p)
    k23 = sin_lat * cos_b * cos_lon - cos_lat * sin_b

    k31 = -sin_b * sin_p * sin_lon + cos_p * cos_lon
    k32 = sin_b * cos_p * sin_lon + sin_p * cos_lon
    k33 = -cos_b * sin_lon

    br = k11 * bxi + k12 * beta + k13 * bzeta
    bt = k21 * bxi + k22 * beta + k23 * bzeta
    bp = k31 * bxi + k32 * beta + k33 * bzeta
    return br, bt, bp
