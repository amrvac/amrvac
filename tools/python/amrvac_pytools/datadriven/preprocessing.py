"""Configurable shared preprocessing for vector magnetograms.

The existing Wiegelmann implementation remains the numerical kernel.  This
module supplies a small, auditable policy layer so Potential, MF,
Optimization, and Grad--Rubin workflows can share the same optional boundary
preprocessing without silently enabling it.
"""

from __future__ import print_function

import hashlib

import numpy as np


PREPROCESSING_RECOMMENDATIONS = {
    "mf": "optional",
    "magnetofriction": "optional",
    "optimization": "recommended",
    "grad_rubin": "recommended",
    "grad-rubin": "recommended",
    "gr": "recommended",
}


class VectorPreprocessingConfig(object):
    """Immutable-in-practice configuration for shared vector preprocessing."""

    def __init__(
        self,
        mode="none",
        mu3=0.1,
        mu4=0.1,
        max_iter=5000,
        tol=1.0e-4,
        dx=1.0,
        dy=1.0,
        geometry_mode=None,
        edge_treatment=None,
        fail_on_nonconvergence=False,
    ):
        mode = str(mode).strip().lower()
        if mode not in ("none", "wiegelmann", "recommended", "conservative"):
            raise ValueError(
                "mode must be none, wiegelmann, conservative, or recommended"
            )
        self.mode = mode
        self.mu3 = float(mu3)
        self.mu4 = float(mu4)
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.dx = float(dx)
        self.dy = float(dy)
        self.geometry_mode = geometry_mode
        self.edge_treatment = edge_treatment
        self.fail_on_nonconvergence = bool(fail_on_nonconvergence)
        if self.mu3 < 0.0 or self.mu4 < 0.0:
            raise ValueError("mu3 and mu4 must be non-negative")
        if self.max_iter < 0 or self.tol < 0.0:
            raise ValueError("max_iter and tol must be non-negative")
        if self.dx <= 0.0 or self.dy <= 0.0:
            raise ValueError("dx and dy must be positive")
        if self.geometry_mode is not None:
            self.geometry_mode = str(self.geometry_mode).strip().lower()
            if self.geometry_mode not in ("legacy_unit_square", "centered"):
                raise ValueError("geometry_mode must be legacy_unit_square or centered")
        if self.edge_treatment is not None:
            self.edge_treatment = str(self.edge_treatment).strip().lower()
            if self.edge_treatment not in ("legacy_periodic", "periodic", "nonperiodic"):
                raise ValueError("edge_treatment must be legacy_periodic, periodic, or nonperiodic")

    def as_dict(self):
        return {
            "mode": self.mode,
            "mu3": self.mu3,
            "mu4": self.mu4,
            "max_iter": self.max_iter,
            "tol": self.tol,
            "dx": self.dx,
            "dy": self.dy,
            "geometry_mode": self.geometry_mode,
            "edge_treatment": self.edge_treatment,
            "fail_on_nonconvergence": self.fail_on_nonconvergence,
        }


def recommend_vector_preprocessing(method):
    """Return a recommendation without changing any workflow setting."""

    key = str(method).strip().lower()
    return PREPROCESSING_RECOMMENDATIONS.get(key, "optional")


def _array_hash(array):
    array = np.asarray(array)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(repr(tuple(array.shape)).encode("ascii"))
    digest.update(np.ascontiguousarray(array).view(np.uint8))
    return digest.hexdigest()


def _finite_summary(*arrays):
    return {
        "all_finite": bool(all(np.all(np.isfinite(np.asarray(item))) for item in arrays)),
        "shapes": [list(np.asarray(item).shape) for item in arrays],
        "sha256": [_array_hash(item) for item in arrays],
    }


def preprocess_vector_magnetogram(
    bx,
    by,
    bz,
    config=None,
    progress=None,
):
    """Apply the shared optional preprocessing policy.

    ``mode='none'`` is a strict copy and is the default.  The actual
    Wiegelmann update is delegated to :func:`preprocess_cartesian_field` so
    existing numerical behavior is retained.
    """

    from .cartesian import cartesian_preprocess_diagnostics, preprocess_cartesian_field

    if config is None:
        config = VectorPreprocessingConfig()
    elif not isinstance(config, VectorPreprocessingConfig):
        config = VectorPreprocessingConfig(**dict(config))

    arrays = [np.asarray(item, dtype=float) for item in (bx, by, bz)]
    if any(item.ndim != 2 for item in arrays):
        raise ValueError("vector magnetogram components must be 2-D")
    if not (arrays[0].shape == arrays[1].shape == arrays[2].shape):
        raise ValueError("vector magnetogram components must have identical shapes")
    if not all(np.all(np.isfinite(item)) for item in arrays):
        raise ValueError("vector magnetogram contains non-finite values")

    effective = _effective_preprocessing_parameters(config)
    before = cartesian_preprocess_diagnostics(
        *arrays,
        dx=effective["dx"],
        dy=effective["dy"],
        geometry_mode=effective["geometry_mode"],
        edge_treatment=effective["edge_treatment"],
    )
    output = [item.copy() for item in arrays]
    numerical_mode = config.mode
    if config.mode == "recommended":
        numerical_mode = "wiegelmann"
    elif config.mode == "conservative":
        numerical_mode = "wiegelmann"
    if numerical_mode == "wiegelmann":
        output[0], output[1], output[2], solver_metrics = preprocess_cartesian_field(
            output[0],
            output[1],
            output[2],
            mu3=effective["mu3"],
            mu4=effective["mu4"],
            max_iter=config.max_iter,
            tol=config.tol,
            dx=effective["dx"],
            dy=effective["dy"],
            geometry_mode=effective["geometry_mode"],
            edge_treatment=effective["edge_treatment"],
            progress=progress,
        )
        if solver_metrics and not solver_metrics.get("converged", False):
            message = (
                "vector preprocessing ended with stop_reason={} after {} iterations"
                .format(solver_metrics.get("stop_reason"), solver_metrics.get("iterations"))
            )
            if config.fail_on_nonconvergence:
                raise RuntimeError(message)
    else:
        solver_metrics = {
            "iterations": 0,
            "converged": True,
            "stop_reason": "disabled",
            "effective_mu1": 0.0,
            "effective_mu2": 0.0,
            "effective_mu3": 0.0,
            "effective_mu4": 0.0,
            "tol": config.tol,
            "max_iter": config.max_iter,
            "geometry_mode": effective["geometry_mode"],
            "edge_treatment": effective["edge_treatment"],
        }
    after = cartesian_preprocess_diagnostics(
        *output,
        dx=effective["dx"],
        dy=effective["dy"],
        geometry_mode=effective["geometry_mode"],
        edge_treatment=effective["edge_treatment"],
    )
    component_changes = [new - old for old, new in zip(arrays, output)]
    warnings = []
    if config.mode != "none" and solver_metrics and not solver_metrics.get("converged", False):
        warnings.append(
            "preprocessing did not converge: {}".format(solver_metrics.get("stop_reason"))
        )
    audit = {
        "schema": "amrvac.vector_preprocessing_audit.v1",
        "configuration": config.as_dict(),
        "effective_parameters": effective,
        "enabled": config.mode != "none",
        "recommendations": {
            method: recommend_vector_preprocessing(method)
            for method in ("MF", "Optimization", "Grad-Rubin")
        },
        "before": before,
        "after": after,
        "solver_metrics": solver_metrics,
        "input": _finite_summary(*arrays),
        "output": _finite_summary(*output),
        "max_absolute_change": [
            float(np.max(np.abs(new - old))) for old, new in zip(arrays, output)
        ],
        "mean_absolute_change": [
            float(np.mean(np.abs(change))) for change in component_changes
        ],
        "correction_sha256": [_array_hash(change) for change in component_changes],
        "Bz_flux": {
            "before": before["Bz_flux"],
            "after": after["Bz_flux"],
            "change": after["Bz_flux"] - before["Bz_flux"],
        },
        "warnings": warnings,
    }
    return output[0], output[1], output[2], audit


def _effective_preprocessing_parameters(config):
    mode = config.mode
    mu3 = config.mu3
    mu4 = config.mu4
    if mode == "conservative":
        mu3 = min(mu3, 0.05)
        mu4 = min(mu4, 0.05)
    geometry_mode = config.geometry_mode
    edge_treatment = config.edge_treatment
    if geometry_mode is None:
        geometry_mode = "legacy_unit_square" if mode == "wiegelmann" else "centered"
    if edge_treatment is None:
        edge_treatment = "legacy_periodic" if mode == "wiegelmann" else "nonperiodic"
    return {
        "mode": mode,
        "mu3": float(mu3),
        "mu4": float(mu4),
        "max_iter": int(config.max_iter),
        "tol": float(config.tol),
        "dx": float(config.dx),
        "dy": float(config.dy),
        "geometry_mode": geometry_mode,
        "edge_treatment": edge_treatment,
        "fail_on_nonconvergence": bool(config.fail_on_nonconvergence),
    }


__all__ = [
    "PREPROCESSING_RECOMMENDATIONS",
    "VectorPreprocessingConfig",
    "preprocess_vector_magnetogram",
    "recommend_vector_preprocessing",
]
