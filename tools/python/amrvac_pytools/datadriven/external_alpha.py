"""Versioned AMRVAC external-alpha boundary product.

The stream format is deliberately small and language-neutral.  It contains
the grid, units, raw/clean alpha, weights, and masks, so Fortran can validate
the product before using it while Python can round-trip and audit it.
"""

from __future__ import print_function

import hashlib
import json
import struct
from pathlib import Path

import numpy as np


MAGIC = b"AMRVAC_EXTERNAL_ALPHA_V1"
MAGIC_BYTES = 32
FORMAT_VERSION = 1
METADATA_BYTES = 4096
HEADER = struct.Struct("<32s4i6d2i")


def _hash_array(array):
    array = np.asarray(array)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(repr(tuple(array.shape)).encode("ascii"))
    digest.update(np.ascontiguousarray(array).view(np.uint8))
    return digest.hexdigest()


def _as_xy_array(array, shape, dtype=float):
    array = np.asarray(array, dtype=dtype)
    if array.shape != shape:
        raise ValueError("external-alpha array shape {} != {}".format(array.shape, shape))
    return np.asfortranarray(array.T)


def write_external_alpha(
    path,
    alpha,
    x,
    y,
    *,
    alpha_raw=None,
    weight=None,
    pil_mask=None,
    valid_mask=None,
    polarity_mask=None,
    unit_length_cm,
    unit_magneticfield_g,
    alpha_unit="1/code_length",
    coordinate_unit="code_length",
    coordinate_unit_cm=None,
    metadata=None,
):
    """Write ``alpha`` with Python ``(ny,nx)`` arrays and x/y coordinates."""

    path = Path(path)
    alpha = np.asarray(alpha, dtype=np.float64)
    if alpha.ndim != 2:
        raise ValueError("alpha must be a 2-D (ny,nx) array")
    ny, nx = alpha.shape
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.shape != (nx,) or y.shape != (ny,):
        raise ValueError("x/y coordinates do not match alpha shape")
    if not np.all(np.isfinite(alpha)) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("external-alpha coordinates and values must be finite")
    if not np.all(np.diff(x) > 0.0) or not np.all(np.diff(y) > 0.0):
        raise ValueError("external-alpha coordinates must be strictly increasing")
    alpha_raw = alpha.copy() if alpha_raw is None else np.asarray(alpha_raw, dtype=np.float64)
    weight = np.ones_like(alpha) if weight is None else np.asarray(weight, dtype=np.float64)
    pil_mask = np.zeros_like(alpha, dtype=np.int32) if pil_mask is None else np.asarray(pil_mask, dtype=np.int32)
    valid_mask = np.isfinite(alpha) if valid_mask is None else np.asarray(valid_mask, dtype=np.int32)
    polarity_mask = np.sign(alpha).astype(np.int32) if polarity_mask is None else np.asarray(polarity_mask, dtype=np.int32)
    for value, name in ((alpha_raw, "alpha_raw"), (weight, "weight"), (pil_mask, "pil_mask"),
                        (valid_mask, "valid_mask"), (polarity_mask, "polarity_mask")):
        if value.shape != alpha.shape:
            raise ValueError("{} shape does not match alpha".format(name))
    if not np.all(np.isfinite(alpha_raw)) or not np.all(np.isfinite(weight)):
        raise ValueError("external-alpha alpha_raw and weight must be finite")
    if np.any(weight < 0.0):
        raise ValueError("external-alpha weights must be non-negative")
    if np.any(~np.isin(polarity_mask, (-1, 0, 1))):
        raise ValueError("polarity_mask must contain only -1, 0, or +1")
    unit_length_cm = float(unit_length_cm)
    unit_magneticfield_g = float(unit_magneticfield_g)
    if unit_length_cm <= 0.0 or unit_magneticfield_g <= 0.0:
        raise ValueError("external-alpha units must be positive")
    if str(alpha_unit) != "1/code_length":
        raise ValueError("external-alpha alpha values must be converted to 1/code_length")
    if coordinate_unit_cm is not None:
        coordinate_unit_cm = float(coordinate_unit_cm)
        if coordinate_unit_cm <= 0.0:
            raise ValueError("coordinate_unit_cm must be positive when provided")
    dx = float(np.mean(np.diff(x))) if nx > 1 else 0.0
    dy = float(np.mean(np.diff(y))) if ny > 1 else 0.0
    xc = float(0.5 * (x[0] + x[-1]))
    yc = float(0.5 * (y[0] + y[-1]))
    metadata_out = dict(metadata or {})
    metadata_out.update({
        "schema": "amrvac.external_alpha.v1",
        "format_version": FORMAT_VERSION,
        "shape_nx_ny": [nx, ny],
        "alpha_unit": str(alpha_unit),
        "coordinate_unit": str(coordinate_unit),
        "coordinate_unit_cm": coordinate_unit_cm,
        "unit_length_cm": unit_length_cm,
        "unit_magneticfield_g": unit_magneticfield_g,
        "arrays_sha256": {
            "alpha_raw": _hash_array(alpha_raw),
            "alpha": _hash_array(alpha),
            "weight": _hash_array(weight),
            "pil_mask": _hash_array(pil_mask),
            "valid_mask": _hash_array(valid_mask),
            "polarity_mask": _hash_array(polarity_mask),
        },
    })
    blob = json.dumps(metadata_out, sort_keys=True, separators=(",", ":")).encode("utf-8")
    if len(blob) > METADATA_BYTES:
        raise ValueError("external-alpha metadata exceeds {} bytes".format(METADATA_BYTES))
    header = HEADER.pack(
        MAGIC.ljust(MAGIC_BYTES, b" "), FORMAT_VERSION, nx, ny, 0,
        unit_length_cm, unit_magneticfield_g, dx, dy, xc, yc, 1, 1,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        handle.write(header)
        handle.write(np.asarray(x, dtype="<f8").tobytes(order="C"))
        handle.write(np.asarray(y, dtype="<f8").tobytes(order="C"))
        for array in (alpha_raw, alpha, weight):
            handle.write(_as_xy_array(array, alpha.shape, np.float64).astype("<f8", copy=False).tobytes(order="F"))
        for array in (pil_mask, valid_mask, polarity_mask):
            handle.write(_as_xy_array(array, alpha.shape, np.int32).astype("<i4", copy=False).tobytes(order="F"))
        handle.write(struct.pack("<i", len(blob)))
        handle.write(blob.ljust(METADATA_BYTES, b"\0"))
    return path


def _read_xy_array(handle, nx, ny, dtype):
    count = nx * ny
    values = np.frombuffer(handle.read(count * np.dtype(dtype).itemsize), dtype=dtype, count=count)
    if values.size != count:
        raise ValueError("truncated external-alpha payload")
    return values.reshape((nx, ny), order="F").T.copy()


def read_external_alpha(path):
    """Read and validate a V1 product, returning arrays in ``(ny,nx)`` order."""

    path = Path(path)
    with path.open("rb") as handle:
        header_bytes = handle.read(HEADER.size)
        if len(header_bytes) != HEADER.size:
            raise ValueError("truncated external-alpha header")
        magic, version, nx, ny, _reserved, unit_length_cm, unit_b, dx, dy, xc, yc, alpha_code, mask_code = HEADER.unpack(header_bytes)
        if magic.rstrip(b" ") != MAGIC or version != FORMAT_VERSION:
            raise ValueError("unsupported external-alpha magic/version")
        if nx < 1 or ny < 1 or unit_length_cm <= 0.0 or unit_b <= 0.0:
            raise ValueError("invalid external-alpha header")
        x = np.frombuffer(handle.read(nx * 8), dtype="<f8", count=nx).copy()
        y = np.frombuffer(handle.read(ny * 8), dtype="<f8", count=ny).copy()
        if x.size != nx or y.size != ny:
            raise ValueError("truncated external-alpha coordinates")
        alpha_raw = _read_xy_array(handle, nx, ny, "<f8")
        alpha = _read_xy_array(handle, nx, ny, "<f8")
        weight = _read_xy_array(handle, nx, ny, "<f8")
        pil_mask = _read_xy_array(handle, nx, ny, "<i4")
        valid_mask = _read_xy_array(handle, nx, ny, "<i4").astype(bool)
        polarity_mask = _read_xy_array(handle, nx, ny, "<i4")
        length_bytes = handle.read(4)
        if len(length_bytes) != 4:
            raise ValueError("truncated external-alpha metadata length")
        metadata_length = struct.unpack("<i", length_bytes)[0]
        blob = handle.read(METADATA_BYTES)
        if metadata_length < 0 or metadata_length > METADATA_BYTES:
            raise ValueError("invalid external-alpha metadata length")
        metadata = json.loads(blob[:metadata_length].decode("utf-8")) if metadata_length else {}
    result = {
        "path": str(path.resolve()),
        "format_version": int(version),
        "nx": int(nx), "ny": int(ny),
        "unit_length_cm": float(unit_length_cm),
        "unit_magneticfield_g": float(unit_b),
        "dx": float(dx), "dy": float(dy), "xc": float(xc), "yc": float(yc),
        "alpha_unit_code": int(alpha_code), "mask_code": int(mask_code),
        "x": x, "y": y,
        "alpha_raw": alpha_raw, "alpha": alpha, "weight": weight,
        "pil_mask": pil_mask.astype(bool), "valid_mask": valid_mask,
        "polarity_mask": polarity_mask,
        "metadata": metadata,
    }
    if not all(np.all(np.isfinite(result[key])) for key in ("x", "y", "alpha_raw", "alpha", "weight")):
        raise ValueError("external-alpha product contains non-finite values")
    if np.any(result["weight"] < 0.0) or np.any(~np.isin(polarity_mask, (-1, 0, 1))):
        raise ValueError("external-alpha masks or weights are invalid")
    if result["metadata"].get("alpha_unit") != "1/code_length":
        raise ValueError("external-alpha metadata must declare alpha_unit='1/code_length'")
    result["file_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def validate_external_alpha_grid(product, x, y, unit_length_cm, unit_magneticfield_g, bz=None, tolerance=1.0e-12):
    """Strictly validate product geometry/units and optional Bz polarity."""

    if isinstance(product, (str, Path)):
        product = read_external_alpha(product)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if product["x"].shape != x.shape or product["y"].shape != y.shape:
        raise ValueError("external-alpha grid shape does not match magnetogram")
    if not np.allclose(product["x"], x, rtol=0.0, atol=tolerance) or not np.allclose(product["y"], y, rtol=0.0, atol=tolerance):
        raise ValueError("external-alpha coordinates do not match magnetogram")
    if not np.isclose(product["unit_length_cm"], unit_length_cm, rtol=0.0, atol=tolerance * max(1.0, abs(unit_length_cm))):
        raise ValueError("external-alpha unit_length does not match magnetogram")
    if not np.isclose(product["unit_magneticfield_g"], unit_magneticfield_g, rtol=0.0, atol=tolerance * max(1.0, abs(unit_magneticfield_g))):
        raise ValueError("external-alpha magnetic unit does not match magnetogram")
    if bz is not None:
        bz = np.asarray(bz, dtype=float)
        if bz.shape != product["alpha"].shape:
            raise ValueError("external-alpha Bz shape does not match product")
        support = product["valid_mask"] & (product["polarity_mask"] != 0)
        if np.any(np.sign(bz[support]) != product["polarity_mask"][support]):
            raise ValueError("external-alpha polarity mask disagrees with Bz")
    return True


__all__ = [
    "FORMAT_VERSION", "MAGIC", "read_external_alpha", "validate_external_alpha_grid",
    "write_external_alpha",
]
