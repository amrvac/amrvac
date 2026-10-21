"""Curvature-equidistributing adaptive node placement.

Self-contained n-D generalisation. Given a tabulated function on an n-D
tensor-product grid, returns adaptive node positions for each axis that
equidistribute the (axis-wise marginal) curvature density.

Theory
------
Equidistribution principle (Babuška-Rheinboldt 1976): place N nodes
{x_i} on [a, b] such that

    integral_{x_i}^{x_{i+1}} m(x) dx = constant   for all i

where m(x) is a monitor function — here we use

    m(x) = |f''(x)| + alpha * max|f''|

with alpha (the baseline floor) preventing zero-density gaps. The
inverse-CDF construction takes a fine reference grid x_ref with values
m_ref, builds the cumulative integral, and inverts at N evenly-spaced
quantiles.

For n-D (>1) on a tensor-product grid, each axis is treated
independently with its own marginal density:

    m_i(x_i) = integral_{x_{j!=i}} |partial^2 f / partial x_i^2| dx_{j!=i}

This is the standard pragmatic approach — full anisotropic n-D mesh
adaptation produces non-tensor-product grids that table interpolation
cannot consume.

Public API
----------
- monitor_density(values, axis, baseline): 1D marginal curvature density
- equidistribute(grid, density, N): inverse-CDF node placement on one axis
- adaptive_nodes(values, axis_grids, N_per_axis, baseline): top-level n-D wrapper

References
----------
- Babuška, I. & Rheinboldt, W.C. 1976, SIAM J. Numer. Anal. 15, 736
- Huang, W. & Russell, R.D. 2010, "Adaptive Moving Mesh Methods", Springer
- de Boor, C. 1973, in "Applied Math Symposia"

Limitations
-----------
- Tensor-product grids only. Diagonal features get smeared across the
  marginal projections of both axes.
- The kernel that consumes the resulting nodes (PCHIP, bicubic, etc.)
  must use actual node spacings, NOT pretend the grid is uniform.
- Constant-curvature inputs fall back to uniform sampling.
"""

from __future__ import annotations

import numpy as np


def monitor_density(values, axis, baseline=0.05):
    """1D marginal curvature density along `axis` of an n-D array.

    Computes |partial^2 values / partial axis^2| via two successive
    np.gradient calls, then sums the absolute value over all other
    axes to produce a 1D density of length values.shape[axis].

    A baseline floor of `baseline * max(density)` is added so smooth
    regions still receive some sampling. Without it, equidistribution
    would collapse multiple nodes onto plateau boundaries.

    Parameters
    ----------
    values : ndarray, shape (n_0, n_1, ..., n_{N-1})
        Function samples on a regular tensor-product grid. Spacing
        does not need to be uniform — np.gradient handles the actual
        positions if you pass coords (we don't, since for the monitor
        we only need a relative measure).
    axis : int in [0, values.ndim)
        Axis along which to measure curvature.
    baseline : float in [0, 1)
        Floor as a fraction of max curvature. Default 0.05.

    Returns
    -------
    ndarray of shape (values.shape[axis],), all entries >= 0.
    """
    if not (0 <= axis < values.ndim):
        raise ValueError(f"axis {axis} out of range for ndim={values.ndim}")
    d2 = np.abs(np.gradient(np.gradient(values, axis=axis), axis=axis))
    if values.ndim > 1:
        marginal = d2.sum(axis=tuple(i for i in range(values.ndim) if i != axis))
    else:
        marginal = d2
    if not np.isfinite(marginal).all():
        # Defensive: replace NaN/inf with the median of the finite portion
        finite = marginal[np.isfinite(marginal)]
        fill = np.median(finite) if finite.size else 1.0
        marginal = np.where(np.isfinite(marginal), marginal, fill)
    m_max = marginal.max()
    if m_max == 0.0:
        # Constant-curvature input: fall back to uniform sampling.
        return np.ones_like(marginal)
    return marginal + baseline * m_max


def equidistribute(grid, density, N):
    """Place N nodes on `grid` such that the integral of `density`
    between consecutive nodes is constant.

    Inverse-CDF sampling. Endpoints grid[0] and grid[-1] are preserved.

    Parameters
    ----------
    grid : ndarray, shape (M,), strictly increasing
        Reference axis on which `density` is sampled.
    density : ndarray, shape (M,), non-negative
        Sampling density per unit grid. Constant input -> uniform output.
    N : int >= 2
        Number of output nodes.

    Returns
    -------
    ndarray of shape (N,), strictly increasing, with nodes[0] = grid[0]
    and nodes[-1] = grid[-1].

    Notes
    -----
    Robust against zero-density plateaus: a tiny x-proportional ramp
    is added to the CDF so np.interp inversion does not collapse
    consecutive linspace targets onto the same boundary x-value.
    """
    grid = np.asarray(grid, dtype=float)
    density = np.asarray(density, dtype=float)
    if grid.shape != density.shape or grid.ndim != 1:
        raise ValueError("grid and density must be 1D arrays of equal length")
    if N < 2:
        raise ValueError("N must be at least 2")
    if np.any(np.diff(grid) <= 0):
        raise ValueError("grid must be strictly increasing")
    if np.any(density < 0):
        raise ValueError("density must be non-negative")

    # Trapezoidal integral of density, evaluated at each grid point.
    # cdf[0] = 0, cdf[-1] = total mass. This is the proper CDF mapping
    # grid -> fractional integral; using bare np.cumsum biases nodes
    # to the left because it sums the trailing edge of each cell only.
    span = grid[-1] - grid[0]
    cdf = np.zeros_like(grid)
    cdf[1:] = np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(grid))
    if cdf[-1] == 0.0:
        # Pathological: all-zero density. Fall back to uniform.
        return np.linspace(grid[0], grid[-1], N)
    cdf = cdf / cdf[-1]

    # Plateau-tie-break: add a tiny x-proportional ramp so cdf is
    # strictly increasing even where density is locally zero.
    eps = 1.0e-9
    if span > 0.0:
        cdf = cdf + eps * (grid - grid[0]) / span
        cdf = cdf / cdf[-1]

    targets = np.linspace(0.0, 1.0, N)
    nodes = np.interp(targets, cdf, grid)

    # Defensive monotonicity: nudge any surviving duplicates apart.
    eps_x = 1.0e-12 * max(span, 1.0)
    for k in range(1, len(nodes)):
        if nodes[k] <= nodes[k - 1]:
            nodes[k] = nodes[k - 1] + eps_x
    return nodes


def adaptive_nodes(values, axis_grids, N_per_axis, baseline=0.05):
    """Top-level n-D wrapper: produce per-axis adaptive node arrays.

    For each axis i in [0, N), compute the marginal curvature density
    of `values` along axis i (integrated over orthogonal axes), then
    equidistribute on `axis_grids[i]` to get N_per_axis[i] nodes.

    Parameters
    ----------
    values : ndarray of shape (n_0, n_1, ..., n_{N-1})
        Function samples on the reference grid. To resolve features
        well, the reference should be at least as fine as the target
        adaptive resolution (a 1.5-2x oversample is typical).
    axis_grids : sequence of N 1D arrays
        axis_grids[i] is the reference grid for axis i, length n_i,
        strictly increasing.
    N_per_axis : int or sequence of N ints
        Number of adaptive nodes per axis. Scalar applies to all axes.
    baseline : float
        Curvature floor (fraction of max). Default 0.05.

    Returns
    -------
    list of N 1D arrays, one per axis. Endpoints preserved.

    Examples
    --------
    >>> import numpy as np
    >>> # 2D table with a sharp Gaussian feature near the centre
    >>> x = np.linspace(0, 1, 200)
    >>> y = np.linspace(0, 1, 200)
    >>> X, Y = np.meshgrid(x, y, indexing='ij')
    >>> f = np.exp(-100 * ((X - 0.5)**2 + (Y - 0.5)**2))
    >>> nodes_x, nodes_y = adaptive_nodes(f, [x, y], 32)
    >>> # nodes cluster near 0.5 along each axis
    """
    ndim = values.ndim
    if len(axis_grids) != ndim:
        raise ValueError(
            f"axis_grids has {len(axis_grids)} entries, expected {ndim}"
        )
    for i, g in enumerate(axis_grids):
        g = np.asarray(g)
        if g.shape != (values.shape[i],):
            raise ValueError(
                f"axis_grids[{i}] has shape {g.shape}, expected ({values.shape[i]},)"
            )
    if np.isscalar(N_per_axis):
        N_per_axis = [int(N_per_axis)] * ndim
    if len(N_per_axis) != ndim:
        raise ValueError(
            f"N_per_axis has {len(N_per_axis)} entries, expected {ndim}"
        )

    return [
        equidistribute(
            np.asarray(axis_grids[i]),
            monitor_density(values, i, baseline=baseline),
            N_per_axis[i],
        )
        for i in range(ndim)
    ]


# ---- Self-tests ------------------------------------------------------

def _selftest():
    """Exercise the public API on synthetic 1D, 2D, 3D inputs."""
    rng = np.random.default_rng(0)

    # 1D: feature at x=0.3
    x = np.linspace(0, 1, 500)
    f1 = np.tanh(50 * (x - 0.3))
    nodes_1d = equidistribute(x, monitor_density(f1, 0), 32)
    assert nodes_1d[0] == x[0] and nodes_1d[-1] == x[-1]
    assert np.all(np.diff(nodes_1d) > 0), "1D nodes must be strictly increasing"
    # Check that nodes cluster near the feature:
    near_feature = ((nodes_1d > 0.25) & (nodes_1d < 0.35)).sum()
    far_from = ((nodes_1d > 0.7) & (nodes_1d < 0.9)).sum()
    assert near_feature > far_from, "1D should cluster nodes near the tanh kink"

    # 2D: Gaussian centred at (0.5, 0.5)
    x = np.linspace(0, 1, 200)
    y = np.linspace(0, 1, 200)
    X, Y = np.meshgrid(x, y, indexing="ij")
    f2 = np.exp(-100 * ((X - 0.5) ** 2 + (Y - 0.5) ** 2))
    nodes_x, nodes_y = adaptive_nodes(f2, [x, y], 32)
    for nodes in (nodes_x, nodes_y):
        assert np.all(np.diff(nodes) > 0)
        # Cluster check: 0.4-0.6 vs 0.0-0.2
        n_centre = ((nodes > 0.4) & (nodes < 0.6)).sum()
        n_edge = ((nodes < 0.2)).sum()
        assert n_centre > n_edge, "2D should cluster near (0.5, 0.5)"

    # 3D: feature on the diagonal
    n = 60
    a = np.linspace(0, 1, n)
    A, B, C = np.meshgrid(a, a, a, indexing="ij")
    f3 = np.exp(-30 * (A + B + C - 1.5) ** 2)  # sheet near A+B+C=1.5
    nlist = adaptive_nodes(f3, [a, a, a], [16, 16, 16])
    for nodes in nlist:
        assert np.all(np.diff(nodes) > 0)

    # Edge case: constant function -> uniform output
    f_const = np.ones((50, 50))
    g = np.linspace(0, 1, 50)
    nx, ny = adaptive_nodes(f_const, [g, g], 10)
    assert np.allclose(nx, np.linspace(0, 1, 10))
    assert np.allclose(ny, np.linspace(0, 1, 10))

    # Edge case: zero-density plateau (the HHe boundary scenario)
    g = np.linspace(0, 1, 100)
    d = np.zeros_like(g)
    d[40:60] = 1.0  # density only in the middle 20%
    nodes = equidistribute(g, d, 32)
    assert np.all(np.diff(nodes) > 0), "duplicate-node fix must hold"

    print("adaptive_grid self-test: ALL PASS")


if __name__ == "__main__":
    _selftest()
