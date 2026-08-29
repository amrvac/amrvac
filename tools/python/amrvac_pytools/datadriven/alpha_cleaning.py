"""Polarity-aware Grad--Rubin alpha conditioning for vector magnetograms.

This module is intentionally separate from shared force/torque preprocessing.
It is designed for observational vector magnetograms; analytic TD99/CFIT
comparisons should use ``preset='none'`` or the dedicated external-alpha
product path, never an analytic rope mask hidden inside this cleaner.
"""

from __future__ import print_function

import hashlib
import json
import warnings
from pathlib import Path

import numpy as np


ALPHA_CLEANING_PRESETS = {
    "none": {
        "bzero_sigma": 0.0,
        "bfull_sigma": 0.0,
        "bzero_fraction": 0.0,
        "bfull_fraction": 0.0,
        "current_sigma_zero": 0.0,
        "current_sigma_full": 0.0,
        "pil_distance_zero": 0.0,
        "pil_distance_full": 0.0,
        "local_mad_sigma": 0.0,
        "local_window_radius": 0,
        "smoothing_sigma_pixels": 0.0,
        "min_component_pixels": 1,
    },
    "conservative": {
        "bzero_sigma": 2.0,
        "bfull_sigma": 4.0,
        "bzero_fraction": 0.005,
        "bfull_fraction": 0.01,
        "current_sigma_zero": 1.0,
        "current_sigma_full": 3.0,
        "pil_distance_zero": 0.5,
        "pil_distance_full": 1.0,
        "local_mad_sigma": 6.0,
        "local_window_radius": 2,
        "smoothing_sigma_pixels": 0.5,
        "min_component_pixels": 2,
    },
    "recommended": {
        "bzero_sigma": 3.0,
        "bfull_sigma": 5.0,
        "bzero_fraction": 0.01,
        "bfull_fraction": 0.02,
        "current_sigma_zero": 2.0,
        "current_sigma_full": 4.0,
        "pil_distance_zero": 1.0,
        "pil_distance_full": 2.0,
        "local_mad_sigma": 5.0,
        "local_window_radius": 3,
        "smoothing_sigma_pixels": 1.0,
        "min_component_pixels": 4,
    },
}


class AlphaCleaningConfig(object):
    """Configuration for ``clean_grad_rubin_alpha``.

    Any keyword matching a preset field is an explicit override.  Noise
    sigmas and resolution are kept separate from the preset so the audit can
    state whether they came from supplied uncertainties or a quiet-region
    estimate.
    """

    def __init__(self, preset="none", **overrides):
        preset = str(preset).strip().lower()
        if preset not in ALPHA_CLEANING_PRESETS:
            raise ValueError("preset must be none, conservative, or recommended")
        values = dict(ALPHA_CLEANING_PRESETS[preset])
        for key in tuple(values):
            if key in overrides and overrides[key] is not None:
                values[key] = overrides.pop(key)
        self.preset = preset
        self.sigma_Bx = overrides.pop("sigma_Bx", None)
        self.sigma_By = overrides.pop("sigma_By", None)
        self.sigma_Bz = overrides.pop("sigma_Bz", None)
        self.resolution_fwhm_pixels = float(
            overrides.pop("resolution_fwhm_pixels", 1.0)
        )
        self.quiet_percentile = float(overrides.pop("quiet_percentile", 30.0))
        self.derivative_spacing_x = float(overrides.pop("derivative_spacing_x", 1.0))
        self.derivative_spacing_y = float(overrides.pop("derivative_spacing_y", 1.0))
        self.coordinate_unit_cm = overrides.pop("coordinate_unit_cm", None)
        self.unit_length_cm = overrides.pop("unit_length_cm", None)
        explicit_alpha_scale = overrides.pop("alpha_scale_to_code", None)
        if self.coordinate_unit_cm is not None:
            self.coordinate_unit_cm = float(self.coordinate_unit_cm)
        if self.unit_length_cm is not None:
            self.unit_length_cm = float(self.unit_length_cm)
        if explicit_alpha_scale is None:
            if self.coordinate_unit_cm is not None and self.unit_length_cm is not None:
                self.alpha_scale_to_code = self.unit_length_cm / self.coordinate_unit_cm
                self.alpha_scale_source = "unit_length_cm/coordinate_unit_cm"
            else:
                self.alpha_scale_to_code = 1.0
                self.alpha_scale_source = "identity"
        else:
            self.alpha_scale_to_code = float(explicit_alpha_scale)
            self.alpha_scale_source = "explicit_alpha_scale_to_code"
            if self.coordinate_unit_cm is not None and self.unit_length_cm is not None:
                expected = self.unit_length_cm / self.coordinate_unit_cm
                if not np.isclose(
                    self.alpha_scale_to_code,
                    expected,
                    rtol=1.0e-12,
                    atol=1.0e-15,
                ):
                    raise ValueError(
                        "alpha_scale_to_code disagrees with unit_length_cm/coordinate_unit_cm"
                    )
        for key, value in values.items():
            setattr(self, key, value)
        if overrides:
            raise TypeError("unknown alpha-cleaning options: {}".format(sorted(overrides)))
        if self.resolution_fwhm_pixels <= 0.0:
            raise ValueError("resolution_fwhm_pixels must be positive")
        if self.derivative_spacing_x <= 0.0 or self.derivative_spacing_y <= 0.0:
            raise ValueError("derivative spacings must be positive")
        if self.alpha_scale_to_code <= 0.0:
            raise ValueError("alpha_scale_to_code must be positive")
        if self.coordinate_unit_cm is not None and self.coordinate_unit_cm <= 0.0:
            raise ValueError("coordinate_unit_cm must be positive when provided")
        if self.unit_length_cm is not None and self.unit_length_cm <= 0.0:
            raise ValueError("unit_length_cm must be positive when provided")
        if not 0.0 <= self.quiet_percentile <= 100.0:
            raise ValueError("quiet_percentile must be in [0, 100]")
        if not (self.bzero_sigma >= 0.0 and self.bfull_sigma >= 0.0):
            raise ValueError("Bz sigma thresholds must be non-negative")
        if not (self.bzero_fraction >= 0.0 and self.bfull_fraction >= 0.0):
            raise ValueError("Bz fractional thresholds must be non-negative")
        if self.bfull_sigma < self.bzero_sigma or self.bfull_fraction < self.bzero_fraction:
            raise ValueError("Bz full thresholds must be >= zero thresholds")
        if self.current_sigma_zero < 0.0 or self.current_sigma_full < 0.0:
            raise ValueError("current significance thresholds must be non-negative")
        if self.current_sigma_full < self.current_sigma_zero:
            raise ValueError("current_sigma_full must be >= current_sigma_zero")
        if self.pil_distance_zero < 0.0 or self.pil_distance_full < 0.0:
            raise ValueError("PIL distance thresholds must be non-negative")
        if self.pil_distance_full < self.pil_distance_zero:
            raise ValueError("pil_distance_full must be >= pil_distance_zero")
        if self.local_mad_sigma < 0.0 or self.local_window_radius < 0:
            raise ValueError("local outlier parameters must be non-negative")
        if self.smoothing_sigma_pixels < 0.0 or self.min_component_pixels < 1:
            raise ValueError("smoothing and component-size parameters are invalid")

    def as_dict(self):
        keys = [
            "preset", "sigma_Bx", "sigma_By", "sigma_Bz",
            "resolution_fwhm_pixels", "quiet_percentile",
            "derivative_spacing_x", "derivative_spacing_y",
            "coordinate_unit_cm", "unit_length_cm", "alpha_scale_to_code",
            "alpha_scale_source",
            "bzero_sigma", "bfull_sigma", "bzero_fraction", "bfull_fraction",
            "current_sigma_zero", "current_sigma_full", "pil_distance_zero",
            "pil_distance_full", "local_mad_sigma", "local_window_radius",
            "smoothing_sigma_pixels", "min_component_pixels",
        ]
        return {key: getattr(self, key) for key in keys}


def _array_hash(array):
    array = np.asarray(array)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(repr(tuple(array.shape)).encode("ascii"))
    digest.update(np.ascontiguousarray(array).view(np.uint8))
    return digest.hexdigest()


def _shift_zero(array, dy, dx):
    array = np.asarray(array)
    result = np.zeros_like(array)
    y0 = max(0, -dy)
    y1 = min(array.shape[0], array.shape[0] - dy)
    x0 = max(0, -dx)
    x1 = min(array.shape[1], array.shape[1] - dx)
    if y0 < y1 and x0 < x1:
        result[y0 + dy:y1 + dy, x0 + dx:x1 + dx] = array[y0:y1, x0:x1]
    return result


def _dilate(mask, radius):
    mask = np.asarray(mask, dtype=bool)
    result = mask.copy()
    for _ in range(int(max(0, radius))):
        expanded = result.copy()
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                expanded |= _shift_zero(result, dy, dx)
        result = expanded
    return result


def _estimate_noise(component, bz, supplied, quiet_percentile):
    if supplied is not None:
        value = float(supplied)
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError("provided magnetic noise sigmas must be positive")
        return value, "provided_uncertainty"
    finite = np.isfinite(component) & np.isfinite(bz)
    if not np.any(finite):
        raise ValueError("cannot estimate noise from an empty component")
    threshold = np.nanpercentile(np.abs(bz[finite]), quiet_percentile)
    quiet = finite & (np.abs(bz) <= threshold)
    values = component[quiet]
    if values.size < 8:
        values = component[finite]
    median = np.median(values)
    mad = np.median(np.abs(values - median)) / 0.6744897501960817
    scale = max(float(np.nanmax(np.abs(component[finite]))), 1.0)
    if not np.isfinite(mad) or mad <= np.finfo(float).eps * scale:
        mad = np.std(values)
        source = "quiet_region_std_fallback"
    else:
        source = "quiet_region_mad"
    return max(float(mad), np.finfo(float).eps * scale), source


def _cosine_ramp(value, zero, full):
    value = np.asarray(value, dtype=float)
    zero = np.asarray(zero, dtype=float)
    full = np.asarray(full, dtype=float)
    denominator = full - zero
    output = np.ones_like(np.broadcast_to(value, np.broadcast(value, zero, full).shape), dtype=float)
    active = np.broadcast_to(denominator > 0.0, output.shape)
    ratio = np.zeros_like(output, dtype=float)
    np.divide(value - zero, denominator, out=ratio, where=active)
    ratio = np.clip(ratio, 0.0, 1.0)
    output = np.where(active, 0.5 * (1.0 - np.cos(np.pi * ratio)), output)
    return output


def _gaussian_smooth(values, sigma):
    values = np.asarray(values, dtype=float)
    if sigma <= 0.0:
        return values.copy()
    try:
        from scipy.ndimage import gaussian_filter

        finite = np.isfinite(values)
        numerator = gaussian_filter(np.where(finite, values, 0.0), sigma=sigma, mode="nearest")
        denominator = gaussian_filter(finite.astype(float), sigma=sigma, mode="nearest")
        output = values.copy()
        valid = denominator > np.finfo(float).eps
        output[valid] = numerator[valid] / denominator[valid]
        return output
    except ImportError:
        finite = np.isfinite(values)
        weights = finite.astype(float)
        return _separable_normalized_smooth(
            np.where(finite, values, 0.0),
            finite,
            weights,
            sigma,
        )


def _pil_mask(bz, bfull, resolution_fwhm_pixels):
    smoothing_sigma = max(float(resolution_fwhm_pixels) / 2.354820045, 0.0)
    seed_bz = _gaussian_smooth(bz, smoothing_sigma)
    finite = np.isfinite(seed_bz)
    reliable_positive = finite & (seed_bz >= bfull)
    reliable_negative = finite & (seed_bz <= -bfull)
    reliable_radius = max(1, int(np.ceil(float(resolution_fwhm_pixels))))
    positive_near = _dilate(reliable_positive, reliable_radius)
    negative_near = _dilate(reliable_negative, reliable_radius)
    crossing = finite & positive_near & negative_near
    positive_sign = finite & (seed_bz > 0.0)
    negative_sign = finite & (seed_bz < 0.0)
    for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        crossing |= (
            positive_sign & _shift_zero(negative_sign, dy, dx)
        ) | (
            negative_sign & _shift_zero(positive_sign, dy, dx)
        )
    reliable_polarity = np.where(
        reliable_positive, 1, np.where(reliable_negative, -1, 0)
    ).astype(np.int8)
    signed_polarity = np.where(bz > 0.0, 1, np.where(bz < 0.0, -1, 0)).astype(np.int8)
    metadata = {
        "pil_smoothing_sigma_pixels": float(smoothing_sigma),
        "reliable_radius_pixels": int(reliable_radius),
        "reliable_positive": int(np.count_nonzero(reliable_positive)),
        "reliable_negative": int(np.count_nonzero(reliable_negative)),
        "sign_crossing_pixels": int(np.count_nonzero(crossing)),
    }
    return crossing, reliable_polarity, signed_polarity, metadata


def _jz_noise_map(shape, sigma_bx, sigma_by, spacing_x, spacing_y):
    ny, nx = shape
    sigma_dby_dx = np.full(shape, np.sqrt(2.0) * sigma_by / (2.0 * spacing_x), dtype=float)
    sigma_dbx_dy = np.full(shape, np.sqrt(2.0) * sigma_bx / (2.0 * spacing_y), dtype=float)
    if nx > 2:
        edge_x = np.sqrt(26.0) * sigma_by / (2.0 * spacing_x)
        sigma_dby_dx[:, 0] = edge_x
        sigma_dby_dx[:, -1] = edge_x
    if ny > 2:
        edge_y = np.sqrt(26.0) * sigma_bx / (2.0 * spacing_y)
        sigma_dbx_dy[0, :] = edge_y
        sigma_dbx_dy[-1, :] = edge_y
    return np.sqrt(sigma_dby_dx * sigma_dby_dx + sigma_dbx_dy * sigma_dbx_dy)


def _effective_cleaning_pixels(config):
    resolution = float(config.resolution_fwhm_pixels)
    return {
        "local_window_radius_pixels": int(np.ceil(float(config.local_window_radius) * resolution)),
        "smoothing_sigma_pixels": float(config.smoothing_sigma_pixels) * resolution,
        "min_component_pixels": max(1, int(np.ceil(float(config.min_component_pixels)))),
        "pil_dilation_radius_pixels": int(np.ceil(resolution)),
    }


def _distance_to_mask(mask):
    if not np.any(mask):
        return np.full(mask.shape, np.inf, dtype=float)
    try:
        from scipy.ndimage import distance_transform_edt

        return distance_transform_edt(~mask)
    except ImportError:
        distance = np.full(mask.shape, np.inf, dtype=float)
        distance[mask] = 0.0
        for step in range(1, max(mask.shape) + 1):
            frontier = np.isinf(distance)
            if not np.any(frontier):
                break
            candidate = np.minimum.reduce([
                _shift_zero(distance, 1, 0), _shift_zero(distance, -1, 0),
                _shift_zero(distance, 0, 1), _shift_zero(distance, 0, -1),
            ]) + 1.0
            updated = frontier & (candidate < distance)
            distance[updated] = candidate[updated]
        return distance


def _masked_local_stats(values, mask, radius):
    if radius <= 0:
        return values.copy(), np.zeros_like(values)
    size = 2 * int(radius) + 1
    padded_values = np.pad(values, radius, mode="edge")
    padded_mask = np.pad(mask, radius, mode="constant", constant_values=False)
    windows = np.lib.stride_tricks.sliding_window_view(
        padded_values, (size, size)
    )
    masks = np.lib.stride_tricks.sliding_window_view(padded_mask, (size, size))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        masked = np.where(masks, windows, np.nan)
        median = np.nanmedian(masked, axis=(-2, -1))
        deviation = np.nanmedian(np.abs(masked - median[..., None, None]), axis=(-2, -1))
    median[~np.isfinite(median)] = 0.0
    deviation[~np.isfinite(deviation)] = 0.0
    return median, deviation


def _label_components(mask):
    mask = np.asarray(mask, dtype=bool)
    try:
        from scipy.ndimage import label

        structure = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=int)
        labels, count = label(mask, structure=structure)
        sizes = np.bincount(labels.ravel())
        return labels, sizes
    except ImportError:
        labels = np.zeros(mask.shape, dtype=np.int32)
        sizes = [0]
        label_id = 0
        for y, x in zip(*np.nonzero(mask)):
            if labels[y, x] != 0:
                continue
            label_id += 1
            stack = [(int(y), int(x))]
            labels[y, x] = label_id
            count = 0
            while stack:
                cy, cx = stack.pop()
                count += 1
                for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ny, nx = cy + dy, cx + dx
                    if (0 <= ny < mask.shape[0] and 0 <= nx < mask.shape[1]
                            and mask[ny, nx] and labels[ny, nx] == 0):
                        labels[ny, nx] = label_id
                        stack.append((ny, nx))
            sizes.append(count)
        return labels, np.asarray(sizes, dtype=int)


def _separable_normalized_smooth(values, mask, weights, sigma):
    if sigma <= 0.0:
        return values.copy()
    radius = max(1, int(np.ceil(3.0 * sigma)))
    coordinates = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (coordinates / sigma) ** 2)
    kernel /= kernel.sum()
    numerator = np.zeros_like(values, dtype=float)
    denominator = np.zeros_like(values, dtype=float)
    for iy, ky in enumerate(kernel):
        for ix, kx in enumerate(kernel):
            shifted_weight = _shift_zero(mask.astype(float) * weights, iy - radius, ix - radius)
            shifted_values = _shift_zero(mask.astype(float) * weights * values, iy - radius, ix - radius)
            factor = ky * kx
            denominator += factor * shifted_weight
            numerator += factor * shifted_values
    output = values.copy()
    valid = denominator > np.finfo(float).eps
    output[valid] = numerator[valid] / denominator[valid]
    return output


def _high_k_fraction(array):
    array = np.asarray(array, dtype=float)
    finite = np.where(np.isfinite(array), array, 0.0)
    spectrum = np.abs(np.fft.fft2(finite)) ** 2
    fy = np.fft.fftfreq(array.shape[0])
    fx = np.fft.fftfreq(array.shape[1])
    ky, kx = np.meshgrid(fy, fx, indexing="ij")
    high = np.sqrt(kx * kx + ky * ky) >= 0.25
    return float(spectrum[high].sum() / max(spectrum.sum(), np.finfo(float).tiny))


def _quantiles(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {str(q): None for q in (0, 1, 5, 25, 50, 75, 95, 99, 100)}
    return {
        str(q): float(np.percentile(values, q)) for q in (0, 1, 5, 25, 50, 75, 95, 99, 100)
    }


def clean_grad_rubin_alpha(bx, by, bz, config=None):
    """Return raw and polarity-aware cleaned alpha plus an audit dictionary."""

    if config is None:
        config = AlphaCleaningConfig()
    elif not isinstance(config, AlphaCleaningConfig):
        config = AlphaCleaningConfig(**dict(config))
    bx, by, bz = [np.asarray(item, dtype=float) for item in (bx, by, bz)]
    if bx.ndim != 2 or bx.shape != by.shape or bx.shape != bz.shape:
        raise ValueError("Bx, By, and Bz must be same-shaped 2-D arrays")
    if not all(np.all(np.isfinite(item)) for item in (bx, by, bz)):
        raise ValueError("alpha cleaning requires finite magnetic components")

    db_y = np.gradient(by, config.derivative_spacing_x, axis=1, edge_order=2)
    db_x = np.gradient(bx, config.derivative_spacing_y, axis=0, edge_order=2)
    jz = db_y - db_x
    alpha_raw_physical = np.zeros_like(bz)
    np.divide(jz, bz, out=alpha_raw_physical, where=bz != 0.0)
    alpha_raw = alpha_raw_physical * config.alpha_scale_to_code
    raw_support = np.isfinite(alpha_raw) & (bz != 0.0)
    polarity_sign = np.where(bz > 0.0, 1, np.where(bz < 0.0, -1, 0)).astype(np.int8)

    if config.preset == "none":
        ones = np.ones_like(bz)
        valid = raw_support.copy()
        audit = _build_audit(
            config, bx, by, bz, jz, alpha_raw, alpha_raw.copy(), ones, valid,
            np.zeros_like(valid), np.zeros_like(valid), np.zeros_like(valid),
            np.zeros_like(valid), np.zeros_like(valid),
            {"Bx": "not_used", "By": "not_used", "Bz": "not_used"},
            effective_pixels=_effective_cleaning_pixels(config),
        )
        return {
            "alpha_raw": alpha_raw,
            "alpha_clean": alpha_raw.copy(),
            "alpha_smooth": alpha_raw.copy(),
            "weight": ones,
            "pil_mask": np.zeros_like(valid),
            "pil_core_mask": np.zeros_like(valid),
            "valid_mask": valid,
            "polarity_mask": polarity_sign,
            "outlier_mask": np.zeros_like(valid),
            "weak_field_mask": np.zeros_like(valid),
            "small_component_mask": np.zeros_like(valid),
            "jz": jz,
            "audit": audit,
        }

    sigma_bx, source_bx = _estimate_noise(bx, bz, config.sigma_Bx, config.quiet_percentile)
    sigma_by, source_by = _estimate_noise(by, bz, config.sigma_By, config.quiet_percentile)
    sigma_bz, source_bz = _estimate_noise(bz, bz, config.sigma_Bz, config.quiet_percentile)
    bmax = max(float(np.max(np.abs(bz))), np.finfo(float).eps)
    bzero = max(config.bzero_sigma * sigma_bz, config.bzero_fraction * bmax)
    bfull = max(config.bfull_sigma * sigma_bz, config.bfull_fraction * bmax)
    if bfull <= bzero:
        bfull = np.nextafter(bzero, np.inf)
    sigma_jz_map = _jz_noise_map(
        bz.shape, sigma_bx, sigma_by,
        config.derivative_spacing_x, config.derivative_spacing_y,
    )
    sigma_jz_scalar = float(np.nanmedian(sigma_jz_map))
    b_weight = _cosine_ramp(np.abs(bz), bzero, bfull)
    weak_field_mask = raw_support & (np.abs(bz) < bzero)
    pil_core, reliable_polarity_mask, polarity_mask, pil_metadata = _pil_mask(
        bz, bfull, config.resolution_fwhm_pixels
    )
    effective_pixels = _effective_cleaning_pixels(config)
    pil_mask = _dilate(pil_core, effective_pixels["pil_dilation_radius_pixels"])
    distance = _distance_to_mask(pil_core)
    pil_weight = _cosine_ramp(
        distance,
        config.pil_distance_zero * config.resolution_fwhm_pixels,
        config.pil_distance_full * config.resolution_fwhm_pixels,
    )
    current_weight = _cosine_ramp(
        np.abs(jz),
        config.current_sigma_zero * sigma_jz_map,
        config.current_sigma_full * sigma_jz_map,
    )
    weight = b_weight * pil_weight * current_weight
    polarity_support = raw_support & (polarity_mask != 0) & (b_weight > 0.0)

    outlier_mask = np.zeros_like(raw_support)
    for sign in (1, -1):
        sign_mask = polarity_support & (polarity_mask == sign)
        median, mad = _masked_local_stats(
            alpha_raw, sign_mask, effective_pixels["local_window_radius_pixels"]
        )
        scale = np.maximum(1.4826 * mad, np.finfo(float).eps)
        outlier_mask |= sign_mask & (
            np.abs(alpha_raw - median) > config.local_mad_sigma * scale
        )

    component_mask = polarity_support & (current_weight > 0.0) & ~outlier_mask
    small_component_mask = np.zeros_like(raw_support)
    for sign in (1, -1):
        sign_mask = component_mask & (polarity_mask == sign)
        labels, sizes = _label_components(sign_mask)
        if sizes.size > 1:
            small_component_mask |= sign_mask & (
                sizes[labels] < effective_pixels["min_component_pixels"]
            )
    smooth_support = component_mask & ~small_component_mask
    alpha_smooth = alpha_raw.copy()
    for sign in (1, -1):
        sign_mask = smooth_support & (polarity_mask == sign)
        alpha_smooth = np.where(
            sign_mask,
            _separable_normalized_smooth(
                alpha_smooth,
                sign_mask,
                np.maximum(weight, 0.0),
                effective_pixels["smoothing_sigma_pixels"],
            ),
            alpha_smooth,
        )
    valid = smooth_support & (weight > 0.0) & np.isfinite(alpha_smooth)
    alpha_clean = np.zeros_like(alpha_raw)
    alpha_clean[valid] = weight[valid] * alpha_smooth[valid]
    audit = _build_audit(
        config, bx, by, bz, jz, alpha_raw, alpha_clean, weight, valid,
        pil_core, pil_mask, weak_field_mask, outlier_mask, small_component_mask,
        {"Bx": source_bx, "By": source_by, "Bz": source_bz},
        noise={"sigma_Bx": sigma_bx, "sigma_By": sigma_by, "sigma_Bz": sigma_bz,
               "sigma_Jz_median": sigma_jz_scalar,
               "sigma_Jz_max": float(np.nanmax(sigma_jz_map)),
               "Bzero": bzero, "Bfull": bfull},
        effective_pixels=effective_pixels,
        pil_metadata=pil_metadata,
        reliable_polarity_mask=reliable_polarity_mask,
    )
    return {
        "alpha_raw": alpha_raw,
        "alpha_clean": alpha_clean,
        "alpha_smooth": alpha_smooth,
        "weight": weight,
        "pil_mask": pil_mask,
        "pil_core_mask": pil_core,
        "valid_mask": valid,
        "polarity_mask": polarity_mask,
        "reliable_polarity_mask": reliable_polarity_mask,
        "outlier_mask": outlier_mask,
        "weak_field_mask": weak_field_mask,
        "small_component_mask": small_component_mask,
        "jz": jz,
        "audit": audit,
    }


def _build_audit(
    config, bx, by, bz, jz, alpha_raw, alpha_clean, weight, valid,
    pil_core_mask, pil_mask, weak_field_mask, outlier_mask, small_component_mask,
    noise_sources, noise=None, effective_pixels=None, pil_metadata=None,
    reliable_polarity_mask=None,
):
    polarity = np.sign(bz).astype(np.int8)
    raw = np.isfinite(alpha_raw) & (bz != 0.0)
    clean_jz = (alpha_clean / config.alpha_scale_to_code) * bz
    abs_raw_current = np.sum(np.abs(jz[raw]))
    abs_clean_current = np.sum(np.abs(clean_jz[valid]))
    reliable = (
        np.zeros_like(valid, dtype=bool)
        if reliable_polarity_mask is None
        else np.asarray(reliable_polarity_mask) != 0
    )
    support = {
        "raw_finite_nonzero_bz": int(np.count_nonzero(raw)),
        "clean_valid": int(np.count_nonzero(valid)),
        "pil_core": int(np.count_nonzero(pil_core_mask)),
        "pil_expanded": int(np.count_nonzero(pil_mask)),
        "reliable_polarity_seed": int(np.count_nonzero(reliable)),
        "weak_field": int(np.count_nonzero(weak_field_mask)),
        "outlier_removed": int(np.count_nonzero(outlier_mask)),
        "small_component_removed": int(np.count_nonzero(small_component_mask)),
        "positive_clean": int(np.count_nonzero(valid & (polarity > 0))),
        "negative_clean": int(np.count_nonzero(valid & (polarity < 0))),
    }
    retention = {
        "signed_current": {},
        "absolute_current": {},
    }
    for sign, name in ((1, "positive"), (-1, "negative")):
        raw_sign = raw & (polarity == sign)
        clean_sign = valid & (polarity == sign)
        raw_signed = float(np.sum(jz[raw_sign]))
        clean_signed = float(np.sum(clean_jz[clean_sign]))
        raw_abs = float(np.sum(np.abs(jz[raw_sign])))
        clean_abs = float(np.sum(np.abs(clean_jz[clean_sign])))
        retention["signed_current"][name] = {
            "raw": raw_signed, "clean": clean_signed,
            "ratio": clean_signed / raw_signed if raw_signed else None,
        }
        retention["absolute_current"][name] = {
            "raw": raw_abs, "clean": clean_abs,
            "ratio": clean_abs / raw_abs if raw_abs else None,
        }
    values = alpha_clean[valid]
    audit = {
        "schema": "amrvac.grad_rubin_alpha_audit.v1",
        "configuration": config.as_dict(),
        "effective_pixels": effective_pixels or {},
        "pil_detection": pil_metadata or {},
        "noise": noise or {},
        "noise_sources": noise_sources,
        "support": support,
        "alpha_quantiles": _quantiles(values),
        "alpha_robust_rms": float(np.sqrt(np.median(values * values))) if values.size else None,
        "alpha_minmax": [float(np.min(values)), float(np.max(values))] if values.size else [None, None],
        "current_retention": retention,
        "high_k_power_fraction": {
            "jz_raw": _high_k_fraction(jz),
            "alpha_raw": _high_k_fraction(alpha_raw),
            "alpha_clean": _high_k_fraction(alpha_clean),
        },
        "net_current_per_polarity": {
            "raw": {
                "positive": float(np.sum(jz[polarity > 0])),
                "negative": float(np.sum(jz[polarity < 0])),
            },
            "clean": {
                "positive": float(np.sum(clean_jz[valid & (polarity > 0)])),
                "negative": float(np.sum(clean_jz[valid & (polarity < 0)])),
            },
        },
        "finite": {
            "Bx": bool(np.all(np.isfinite(bx))),
            "By": bool(np.all(np.isfinite(by))),
            "Bz": bool(np.all(np.isfinite(bz))),
            "Jz": bool(np.all(np.isfinite(jz))),
            "alpha_raw": bool(np.all(np.isfinite(alpha_raw))),
            "alpha_clean": bool(np.all(np.isfinite(alpha_clean))),
            "weight": bool(np.all(np.isfinite(weight))),
        },
        "input_sha256": [_array_hash(item) for item in (bx, by, bz)],
        "output_sha256": [_array_hash(item) for item in (alpha_raw, alpha_clean, weight)],
        "clean_jz_sha256": _array_hash(clean_jz),
        "total_abs_current_retention": float(abs_clean_current / abs_raw_current)
        if abs_raw_current else None,
    }
    audit["audit_sha256"] = hashlib.sha256(
        json.dumps(audit, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return audit


def write_alpha_cleaning_audit(path, audit):
    """Write a deterministic alpha-cleaning audit JSON and return its path."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(audit, indent=2, sort_keys=True) + "\n"
    path.write_text(payload, encoding="utf-8")
    return path


__all__ = [
    "ALPHA_CLEANING_PRESETS", "AlphaCleaningConfig",
    "clean_grad_rubin_alpha", "write_alpha_cleaning_audit",
]
