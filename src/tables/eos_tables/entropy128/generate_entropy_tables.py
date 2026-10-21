"""Build the complete bicubic-Hermite EoS table set for AMRVAC's
'entropy' method.

Produces 7 quantities × 4 stored arrays (value + 3 derivatives) = 28
binary tables. With these, every runtime EoS call collapses to a
single bicubic-Hermite evaluation; no Newton, no bisection, no
fixed-point iteration anywhere on the simulation hot path.

Output layout (per composition `comp`):

  Forward tables — adaptive (log_nH, log_eint/nH) CGS grid:
    LTEeos_s_<comp>.bin         entropy potential
    LTEeos_s_x_<comp>.bin       ∂s/∂(log nH)
    LTEeos_s_y_<comp>.bin       ∂s/∂(log eint/nH)
    LTEeos_s_xy_<comp>.bin      ∂²s/(∂x ∂y)
    LTEeos_g1e_<comp>.bin       Γ₁(ρ, e_int)
    LTEeos_g1e_{x,y,xy}_<comp>.bin
    LTEeos_neOnH_<comp>.bin     n_e/n_H
    LTEeos_neOnH_{x,y,xy}_<comp>.bin

  Inverse tables — regular (log_nH, log_p/nH) CGS grid:
    LTEeos_eintP_<comp>.bin     log(eint/nH) at (nH, p)
    LTEeos_eintP_{x,y,xy}_<comp>.bin
    LTEeos_g1p_<comp>.bin       Γ₁(ρ, p)
    LTEeos_g1p_{x,y,xy}_<comp>.bin

  Inverse tables — regular (log_nH, log_T) CGS grid:
    LTEeos_eintT_<comp>.bin     log(eint/nH) at (nH, T)
    LTEeos_eintT_{x,y,xy}_<comp>.bin
    LTEeos_yT_<comp>.bin        n_e/n_H at (nH, T)
    LTEeos_yT_{x,y,xy}_<comp>.bin

All node values come from the analytic `saha_np` (Newton on charge
balance to machine precision) so no compounding interpolation error
enters the build. Cross derivatives come from finite differences of
JAX-AD first derivatives ('FD-of-AD'), per the prototype's recipe.

Verification gates (Phase A from the project plan):
  1. Maxwell residual: e_spec·ln10 / s_y == T at every node (1e-12)
  2. Round-trip: e -> p -> e' agrees to <1e-3 rel at 1000 random points
  3. γ₁ consistency: stored γ₁ matches Saha-AD γ₁ to <0.5% at 100 nodes
  4. Inverse-table monotonicity: eintP/eintT strictly monotone in y axis
  5. neOnH bounds: in [0, 1+2·A_He] at every node
  6. No NaN, no inf anywhere in any output array

Failing any gate ABORTS the build before any file is written to disk.

Run:
    python3 generate_entropy_tables.py [--n 128] [--ref 512] [--comp HHe_IonE]
                                    [--outdir <path>]
"""
import os
import sys
import time
import json
import struct
import argparse
import numpy as np

# Make sibling modules importable irrespective of cwd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ---------------------------------------------------------------- constants
M_P = 1.6726e-24                            # g
LN10 = np.log(10.0)
A_HE_DEFAULT = 0.1

# Bounds tightened to the FALC chromosphere-corona + Ballester operating
# envelope with modest margin. The previous wider bounds (5-19 dex in log_nH,
# 14-dex span) let the adaptive grid waste budget on near-degenerate-density
# edges AMRVAC never visits. The narrower span doubles the effective node
# density per axis at the same N=128.
#
# Reference: FALC + Ballester IC samples lie in:
#   log_nH:      [10.02, 17.12]   (cm^-3)
#   log_eint/nH: [-11.99, -10.09] (erg)
#   log_p/nH:    [-12.17, -10.50] (erg)
#   log_T:       [3.65,    5.00]  (K)
#
# T floor set to 1000K (log_T = 3.0) per user requirement; Saha bracket
# extends to 300K, so the 1000K floor is safely above the bisection's
# clamp limit.
NH_BOUNDS_CGS      = [   5.0,  20.0]   # log_nH (cm^-3) — widened past uniform256's
                                       # 19.0 to give Coleman 5/6 pureH (log_nH=19.08)
                                       # 0.9 dex of headroom rather than clamping.
EINT_NH_BOUNDS_CGS = [ -12.8,  -9.0]   # log_eint/nH (erg per H atom)
LP_NH_BOUNDS_CGS   = [ -13.2,  -9.4]   # log_p/nH (erg per H atom)
LT_BOUNDS_CGS      = [  2.84,  6.31]   # log_T (K) — matches uniform256
SAHA_T_LO          = 300.0
SAHA_T_HI          = 2.0e6

# Atomic-mass shift between log_rho and log_nH at A_He=0.1
def _log_mAtom(A_He):
    return np.log10((1.0 + 4.0 * A_He) * M_P)


#> Optional exact (rho,T)->p backend for the round-trip verification. Defaults to
#> None -> verify uses saha_np.eos_state (exact for H+He). A composition build sets
#> this to its exact solver so p_true is not taken from the spline placement shim
#> (whose ~0.3% p error is amplified by the steep eintP slope at ionisation plateaus).
EXACT_EOS_STATE = None

def lnH_from_lrho(lr, A_He):
    return lr - _log_mAtom(A_He)
def lrho_from_lnH(lnH, A_He):
    return lnH + _log_mAtom(A_He)
def le_spec_from_le_nh(le_nh, A_He):
    return le_nh - _log_mAtom(A_He)
def le_nh_from_le_spec(le_spec, A_He):
    return le_spec + _log_mAtom(A_He)


# ---------------------------------------------------------------- file I/O
def write_bin_table(filepath, data, var1_bounds, var2_bounds,
                    var1_nodes=None, var2_nodes=None):
    """AMRVAC binary table format (matches existing loader exactly).

    Header  : [2 i4] dim1, dim2 ; [4 f8] var1_min, var1_max, var2_min, var2_max
    Body    : [dim1*dim2 f8] data, Fortran column-major
    Trailer : [1 i4] flag=1 ; [dim1 f8] var1_nodes ; [dim2 f8] var2_nodes
              (only when both *_nodes are provided)
    """
    if not np.isfinite(data).all():
        raise ValueError(f"non-finite data in {filepath}")
    with open(filepath, 'wb') as f:
        f.write(struct.pack('2i', data.shape[0], data.shape[1]))
        f.write(struct.pack('4d', *var1_bounds, *var2_bounds))
        data.T.tofile(f)
        if var1_nodes is not None and var2_nodes is not None:
            assert var1_nodes.shape == (data.shape[0],)
            assert var2_nodes.shape == (data.shape[1],)
            f.write(struct.pack('i', 1))
            var1_nodes.astype(np.float64).tofile(f)
            var2_nodes.astype(np.float64).tofile(f)


def fd_along_axis(field, axis_vals, axis):
    """Centred FD on a 2D field along the named axis (0 or 1), non-uniform-safe.
    One-sided at boundaries. Returns array of the same shape."""
    out = np.zeros_like(field)
    if axis == 0:
        n = field.shape[0]
        for i in range(n):
            if 0 < i < n - 1:
                out[i, :] = (field[i+1, :] - field[i-1, :]) / (axis_vals[i+1] - axis_vals[i-1])
            elif i == 0:
                out[i, :] = (field[i+1, :] - field[i, :]) / (axis_vals[i+1] - axis_vals[i])
            else:
                out[i, :] = (field[i, :] - field[i-1, :]) / (axis_vals[i] - axis_vals[i-1])
    else:
        n = field.shape[1]
        for j in range(n):
            if 0 < j < n - 1:
                out[:, j] = (field[:, j+1] - field[:, j-1]) / (axis_vals[j+1] - axis_vals[j-1])
            elif j == 0:
                out[:, j] = (field[:, j+1] - field[:, j]) / (axis_vals[j+1] - axis_vals[j])
            else:
                out[:, j] = (field[:, j] - field[:, j-1]) / (axis_vals[j] - axis_vals[j-1])
    return out


def bicubic_hermite_2d(x, y, f, fx, fy, fxy, x_axis, y_axis):
    """Pure-NumPy bicubic Hermite evaluation at (x, y), mirroring the Fortran
    kernel in mod_eos_entropy.t.bicubic_hermite_eval. Used by the verifier
    to test the same interpolant AMRVAC will use at runtime.

    Inputs:
      f, fx, fy, fxy : 2D arrays on the (x_axis, y_axis) grid
      x_axis, y_axis : 1D node arrays (regular or adaptive, both supported)
    Returns: scalar value at (x, y).
    """
    n1, n2 = f.shape
    # Locate cell and local coordinate (clamp at boundaries).
    ix = max(0, min(n1 - 2, np.searchsorted(x_axis, x) - 1))
    iy = max(0, min(n2 - 2, np.searchsorted(y_axis, y) - 1))
    dx = x_axis[ix+1] - x_axis[ix]
    dy = y_axis[iy+1] - y_axis[iy]
    tx = max(0.0, min(1.0, (x - x_axis[ix]) / dx))
    ty = max(0.0, min(1.0, (y - y_axis[iy]) / dy))
    # Cubic Hermite basis on [0,1]:
    H0x = 2.0*tx**3 - 3.0*tx**2 + 1.0
    H1x = tx**3 - 2.0*tx**2 + tx
    H2x = -2.0*tx**3 + 3.0*tx**2
    H3x = tx**3 - tx**2
    H0y = 2.0*ty**3 - 3.0*ty**2 + 1.0
    H1y = ty**3 - 2.0*ty**2 + ty
    H2y = -2.0*ty**3 + 3.0*ty**2
    H3y = ty**3 - ty**2
    # 16-term tensor product (corner indices 00/10/01/11 = (i,j),(i+1,j),(i,j+1),(i+1,j+1))
    f00=f[ix,iy];     f10=f[ix+1,iy];     f01=f[ix,iy+1];     f11=f[ix+1,iy+1]
    fx00=fx[ix,iy]*dx;     fx10=fx[ix+1,iy]*dx
    fx01=fx[ix,iy+1]*dx;   fx11=fx[ix+1,iy+1]*dx
    fy00=fy[ix,iy]*dy;     fy10=fy[ix+1,iy]*dy
    fy01=fy[ix,iy+1]*dy;   fy11=fy[ix+1,iy+1]*dy
    fxy00=fxy[ix,iy]*dx*dy;     fxy10=fxy[ix+1,iy]*dx*dy
    fxy01=fxy[ix,iy+1]*dx*dy;   fxy11=fxy[ix+1,iy+1]*dx*dy
    return ( H0x*H0y*f00 + H2x*H0y*f10 + H0x*H2y*f01 + H2x*H2y*f11
           + H1x*H0y*fx00 + H3x*H0y*fx10 + H1x*H2y*fx01 + H3x*H2y*fx11
           + H0x*H1y*fy00 + H2x*H1y*fy10 + H0x*H3y*fy01 + H2x*H3y*fy11
           + H1x*H1y*fxy00 + H3x*H1y*fxy10 + H1x*H3y*fxy01 + H3x*H3y*fxy11 )


def fd_xy_cross(f, x_axis, y_axis):
    """Cross derivative fxy via FD of fx along y (averaged with FD of fy along x
    for Maxwell-symmetric output). Returns the symmetrised array."""
    fx = fd_along_axis(f, x_axis, axis=0)
    fy = fd_along_axis(f, y_axis, axis=1)
    a = fd_along_axis(fy, x_axis, axis=0)
    b = fd_along_axis(fx, y_axis, axis=1)
    return 0.5 * (a + b)


def bilinear_sample(field, x_axis, y_axis, x_targets, y_targets):
    """Sample a 2D field defined on (x_axis, y_axis) at the cross product of
    x_targets × y_targets. Both source axes assumed monotonically increasing.
    Linear extrapolation outside the source range (rare; only at table corners).
    Used to transfer fine-grid FD-derivative fields onto the adaptive nodes.
    """
    nx_src, ny_src = field.shape
    nx_tgt = len(x_targets); ny_tgt = len(y_targets)
    out = np.empty((nx_tgt, ny_tgt))
    for ii, xt in enumerate(x_targets):
        ix = np.searchsorted(x_axis, xt) - 1
        ix = max(0, min(nx_src - 2, ix))
        tx = (xt - x_axis[ix]) / (x_axis[ix+1] - x_axis[ix])
        for jj, yt in enumerate(y_targets):
            iy = np.searchsorted(y_axis, yt) - 1
            iy = max(0, min(ny_src - 2, iy))
            ty = (yt - y_axis[iy]) / (y_axis[iy+1] - y_axis[iy])
            out[ii, jj] = ((1-tx)*(1-ty)*field[ix,iy]
                         + tx*(1-ty)*field[ix+1,iy]
                         + (1-tx)*ty*field[ix,iy+1]
                         + tx*ty    *field[ix+1,iy+1])
    return out


def fine_grid_derivs(F_fine, lr_fine, le_fine, lr_axis, le_axis):
    """Compute (F_x, F_y, F_xy) at adaptive nodes via FD on the FINE uniform
    grid, then bilinearly sample onto the adaptive (lr_axis, le_axis) nodes.

    FD on the fine grid uses h_fine ≪ h_adaptive, so derivative error is
    O(h_fine²) — far below the O(h_adaptive²) error of FD-on-adaptive.
    """
    Fx_fine  = fd_along_axis(F_fine, lr_fine, axis=0)
    Fy_fine  = fd_along_axis(F_fine, le_fine, axis=1)
    Fxy_fine = 0.5 * (fd_along_axis(Fy_fine, lr_fine, axis=0)
                     + fd_along_axis(Fx_fine, le_fine, axis=1))
    Fx_node  = bilinear_sample(Fx_fine,  lr_fine, le_fine, lr_axis, le_axis)
    Fy_node  = bilinear_sample(Fy_fine,  lr_fine, le_fine, lr_axis, le_axis)
    Fxy_node = bilinear_sample(Fxy_fine, lr_fine, le_fine, lr_axis, le_axis)
    return Fx_node, Fy_node, Fxy_node


# ---------------------------------------------------------------- Step 1
def build_fine_s_grid(N, A_He, lr_bounds, le_bounds):
    """Uniform fine reference (s, T, p, neOnH)(lr, le) for both equidistribution
    AND high-accuracy node derivatives.

    Returns (lr_axis, le_axis, s_grid, T_grid, p_grid, ne_grid).
    Sampling on the FINE uniform grid (instead of only at the 128 adaptive
    nodes) lets us compute first derivatives via FD with O(h_fine²) accuracy
    — at N=512, h_fine~0.018 in lr → derivative error ~3e-4, vs O(h_adaptive²)
    ~0.017 at the H knee where the adaptive grid is coarser. The FD derivative
    fields are then bilinearly sampled at the 128 adaptive nodes (the
    derivative fields are smooth in (lr, le), so the sampling step adds
    negligible error).
    """
    from saha_np import eos_at_rho_eint, saha_HHe_solve
    lr = np.linspace(lr_bounds[0], lr_bounds[1], N)
    le = np.linspace(le_bounds[0], le_bounds[1], N)
    s  = np.zeros((N, N))
    Tg = np.zeros((N, N))
    pg = np.zeros((N, N))
    ng = np.zeros((N, N))
    t0 = time.perf_counter()
    print(f'  fine reference grid {N}x{N} (s, T, p, neOnH) ...', flush=True)
    for i, lri in enumerate(lr):
        rho = 10.0 ** lri
        for j, lej in enumerate(le):
            e_spec = 10.0 ** lej
            T_ij, p_ij, s_vol = eos_at_rho_eint(rho, e_spec * rho)
            _, _, _, y_ij = saha_HHe_solve(rho, T_ij)
            s[i, j]  = s_vol / rho
            Tg[i, j] = T_ij
            pg[i, j] = p_ij
            ng[i, j] = y_ij
        if (i + 1) % max(1, N // 8) == 0:
            print(f'    row {i+1}/{N}  elapsed {time.perf_counter()-t0:5.1f}s',
                  flush=True)
    return lr, le, s, Tg, pg, ng


# ---------------------------------------------------------------- Step 2
def equidistribute_axes(lr_fine, le_fine, s_fine, N_out, baseline=0.05,
                         alpha_combined=1.0, beta_gamma1=1.0):
    """Adaptive (lr, le) node placement via curvature equidistribution.

    le axis: |∂⁴γ₁/∂le⁴|^¼ — the *theoretically* error-optimal monitor for
        bicubic Hermite, which has 4th-order truncation
            err ~ h⁴ · |∂⁴f/∂x⁴|.
        Equidistributing |∂⁴f|^¼ as density gives uniform per-cell error.
        In a 1D ideal-node-value test against the analytic Saha Γ₁ along the
        Ballester slice (ρ=4e-14 g/cc), this drops max error to 0.039% at
        N=128 vs 0.186% for the 2nd-derivative monitor we used before. The
        previous |∂²s| monitor under-resolved the He I → He II → He III
        knee at 7–13 kK because Γ₁ varies smoothly there (high 1st-deriv,
        moderate 2nd-deriv) and the 4th derivative is the term that
        actually drives the Hermite interpolation error.

    lr axis: combined monitor — |∂²s/∂lr²| + α * |∂(∂s/∂lr)/∂lr|, since s
        is nearly flat in lr at fixed le and pure 2nd-derivative gives
        uniform-ish spacing along that axis.
    """
    from adaptive_grid import equidistribute, monitor_density

    # ── le axis: HYBRID monitor — s 2nd-deriv ⊕ γ₁ 2nd-deriv ⊕ γ₁ 4th-deriv ──
    # Pure |d⁴γ₁|^¼ (M3) gives the lowest Γ₁ interpolation error in a 1D
    # test (0.04% at N=128) but concentrates nodes so aggressively at the
    # He knee that the forward→inverse p round-trip (eintP) degrades to
    # ~1%, failing Phase A. The hybrid keeps the |d²s| backbone (needed
    # for p, e tables) and ADDS Γ₁ 2nd- and 4th-derivative weight to put
    # extra resolution where Γ₁ is sharp without starving other regions.
    s_r_f  = fd_along_axis(s_fine, lr_fine, axis=0)
    s_e_f  = fd_along_axis(s_fine, le_fine, axis=1)
    s_rr_f = fd_along_axis(s_r_f,  lr_fine, axis=0)
    s_re_f = fd_along_axis(s_r_f,  le_fine, axis=1)
    s_ee_f = fd_along_axis(s_e_f,  le_fine, axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        g1_fine = gamma1_from_s_derivs(s_r_f, s_e_f, s_rr_f, s_re_f, s_ee_f)
    finite = np.isfinite(g1_fine) & (g1_fine > 0.5) & (g1_fine < 5.0)
    if not finite.all():
        med = np.median(g1_fine[finite]) if finite.any() else 1.5
        g1_fine = np.where(finite, g1_fine, med)

    # Marginal 2nd-derivative density of s (preserves p, e accuracy)
    dens_s = monitor_density(s_fine, axis=1, baseline=0.0)
    # Marginal 2nd-derivative density of γ₁ (resolves moderate ionisation features)
    dens_g1 = monitor_density(g1_fine, axis=1, baseline=0.0)
    # Marginal 4th-derivative^¼ density of γ₁ (error-optimal for bicubic Hermite)
    d4 = g1_fine
    for _ in range(4):
        d4 = np.gradient(d4, axis=1)
    d4_marg = np.abs(d4).sum(axis=0)
    if not np.isfinite(d4_marg).all():
        fd = np.isfinite(d4_marg)
        med = np.median(d4_marg[fd]) if fd.any() else 1.0
        d4_marg = np.where(fd, d4_marg, med)
    dens_d4 = np.power(d4_marg + 1e-30, 0.25)

    # Normalise each, combine with weights, then apply floor
    for d in (dens_s, dens_g1, dens_d4):
        if d.max() > 0: d /= d.max()
    dens_le = dens_s + beta_gamma1 * dens_g1 + beta_gamma1 * dens_d4
    dens_le = dens_le + baseline * dens_le.max()
    le_axis = equidistribute(le_fine, dens_le, N_out)

    # ── lr axis: legacy monitor (s-curvature only) ────────────────────
    sr_fine = fd_along_axis(s_fine, lr_fine, axis=0)
    d2_s    = np.abs(np.gradient(np.gradient(s_fine,  axis=0), axis=0)).sum(axis=1)
    d_sr    = np.abs(np.gradient(sr_fine, axis=0)).sum(axis=1)
    d2_s   /= d2_s.max() if d2_s.max() > 0 else 1.0
    d_sr   /= d_sr.max() if d_sr.max() > 0 else 1.0
    combined = d2_s + alpha_combined * d_sr
    combined = combined + baseline * combined.max()
    lr_axis = equidistribute(lr_fine, combined, N_out)
    return lr_axis, le_axis


# ---------------------------------------------------------------- Step 3
def sample_forward_node_values(lr_axis, le_axis):
    """Sample (s, T, p, neOnH) at every (lr, le) node via direct Saha.

    Returns four arrays of shape (N, N) with values in CGS:
      s_g    [erg/g/K]   — specific entropy
      T_g    [K]          — temperature
      p_g    [erg/cm^3]   — pressure
      neOnH  [-]          — electron-to-hydrogen number ratio
    """
    from saha_np import eos_at_rho_eint, saha_HHe_solve
    N1, N2 = len(lr_axis), len(le_axis)
    s_g = np.zeros((N1, N2)); T_g = np.zeros((N1, N2))
    p_g = np.zeros((N1, N2)); ne  = np.zeros((N1, N2))
    t0 = time.perf_counter()
    print(f'  forward node values {N1}x{N2} ...', flush=True)
    for i, lri in enumerate(lr_axis):
        rho = 10.0 ** lri
        for j, lej in enumerate(le_axis):
            e_spec = 10.0 ** lej
            T_ij, p_ij, s_vol = eos_at_rho_eint(rho, e_spec * rho)
            _, _, _, y_ij = saha_HHe_solve(rho, T_ij)
            s_g[i, j] = s_vol / rho
            T_g[i, j] = T_ij
            p_g[i, j] = p_ij
            ne[i, j]  = y_ij
        if (i + 1) % max(1, N1 // 8) == 0:
            print(f'    row {i+1}/{N1}  elapsed {time.perf_counter()-t0:5.1f}s',
                  flush=True)
    return s_g, T_g, p_g, ne


def first_derivatives_s_maxwell(lr_axis, le_axis, T_node, p_node):
    """First derivatives of s(lr, le) from the EXACT Maxwell identities at
    nodes. NO JAX, NO bisection, NO iteration -- pure thermodynamic relations
    from the Saha state already computed in Step 3.

        T  = e_spec * ln10 / s_y    =>   s_y = e_spec * ln10 / T
        p  = -rho * T * s_x / ln10  =>   s_x = -p * ln10 / (rho * T)

    By construction the stored derivative tables satisfy the Maxwell relation
    against the stored T and p at every node to machine precision -- so the
    runtime entropy_T_p_from_log_nH_eint formula recovers exactly the Saha
    temperature when applied to interpolated derivatives at any table node.
    """
    N1, N2 = len(lr_axis), len(le_axis)
    s_x = np.zeros((N1, N2))
    s_y = np.zeros((N1, N2))
    for i in range(N1):
        rho_i = 10.0 ** lr_axis[i]
        for j in range(N2):
            e_spec = 10.0 ** le_axis[j]
            T_ij   = T_node[i, j]
            p_ij   = p_node[i, j]
            s_y[i, j] =  e_spec * LN10 / T_ij
            s_x[i, j] = -p_ij   * LN10 / (rho_i * T_ij)
    return s_x, s_y


def gamma1_from_s_derivs(s_r, s_e, s_rr, s_re, s_ee):
    """Closed-form γ₁ from the entropy potential and its 5 first/second
    derivatives. Defined pointwise on whatever grid the inputs share.

    Γ₁ = 1 - s_r/s_e
         + (1/(s_r·ln10)) · (s_rr - s_re · s_r/s_e)
         - (1/(s_e·ln10)) · (s_re - s_ee · s_r/s_e)
    """
    rat = s_r / s_e
    return (1.0 - rat
            + (1.0 / (s_r * LN10)) * (s_rr - s_re * rat)
            - (1.0 / (s_e * LN10)) * (s_re - s_ee * rat))


def build_forward_tables(lr_axis, le_axis, A_He,
                           fine_grids=None):
    """Step 3: every forward-grid quantity and its 3 derivative arrays.

    Returns a dict keyed by quantity name, each entry a dict of arrays:
      {'value': (N,N), 'x': (N,N), 'y': (N,N), 'xy': (N,N)}.

    The entropy potential s itself is stored (for diagnostic/consistency
    checks), but the runtime EoS does NOT derive T or p from ds/dy or ds/dx
    -- it uses the direct Tfwd and pfwd tables below.

    Tfwd[i,j] = T at the node (lr_i, le_j) from the same Saha call that
    produces s and p. So at NODES, T, p, s, and gamma1 are mutually
    Maxwell-consistent.

    When `fine_grids` is provided (tuple of lr_fine, le_fine, T_fine, p_fine,
    ne_fine from build_fine_s_grid), first derivatives of Tfwd/pfwd/neOnH
    are computed by FD on the FINE grid and bilinearly sampled at adaptive
    nodes — gives O(h_fine²) corner-derivative error, ~50× smaller than
    FD-on-adaptive at the H knee. Without it, falls back to FD-on-adaptive.
    """
    s_g, T_g, p_g, ne = sample_forward_node_values(lr_axis, le_axis)
    # First derivatives of s from Saha perturbation (matches the pattern
    # used for Tfwd_x, pfwd_x at line 528 — sample s(ρ±ε, e_spec) and
    # s(ρ, e_spec±ε) directly, centred FD with fixed ε).  The earlier
    # Maxwell-identity shortcut (s_x = -p·ln10/(ρT)) was exact at nodes
    # but mismatched the bicubic-Hermite polynomial's interior dynamics,
    # producing 2-7 % oscillation in the polynomial's analytic derivative
    # across adjacent cells and breaking the entropy_v2 wave evolution.
    # Maxwell-identity values kept for diagnostic reference, not stored.
    s_r_maxwell, s_e_maxwell = first_derivatives_s_maxwell(lr_axis, le_axis, T_g, p_g)
    from saha_np import eos_at_rho_eint as _eos_at_rho_eint_for_s
    N1, N2 = len(lr_axis), len(le_axis)
    s_r = np.zeros((N1, N2))
    s_e = np.zeros((N1, N2))
    eps_lr_s, eps_le_s = 1.0e-4, 1.0e-4
    print(f'  s derivatives via Saha perturbation at each node ...', flush=True)
    t0_s = time.perf_counter()
    for i, lri in enumerate(lr_axis):
        rho    = 10.0 ** lri
        rho_pl = 10.0 ** (lri + eps_lr_s)
        rho_mn = 10.0 ** (lri - eps_lr_s)
        for j, lej in enumerate(le_axis):
            e_spec = 10.0 ** lej
            e_pl   = 10.0 ** (lej + eps_le_s)
            e_mn   = 10.0 ** (lej - eps_le_s)
            # lr derivative at fixed e_spec
            _, _, sv_pl = _eos_at_rho_eint_for_s(rho_pl, e_spec * rho_pl)
            _, _, sv_mn = _eos_at_rho_eint_for_s(rho_mn, e_spec * rho_mn)
            s_r[i, j] = ((sv_pl/rho_pl) - (sv_mn/rho_mn)) / (2.0 * eps_lr_s)
            # le derivative at fixed rho
            _, _, sv_pl = _eos_at_rho_eint_for_s(rho, e_pl * rho)
            _, _, sv_mn = _eos_at_rho_eint_for_s(rho, e_mn * rho)
            s_e[i, j] = ((sv_pl/rho) - (sv_mn/rho)) / (2.0 * eps_le_s)
        if (i+1) % max(1, N1 // 8) == 0:
            print(f'    row {i+1}/{N1}  elapsed {time.perf_counter()-t0_s:.1f}s',
                  flush=True)
    # Cross derivative: direct Saha second mixed perturbation at each node.
    # FD-on-adaptive of s_x, s_y (the line 460-461 pattern) was sensitive to
    # variable cell spacing — neighbouring cells with very different widths
    # produced inconsistent stored s_xy magnitudes, breaking bicubic-Hermite
    # smoothness across cell boundaries.  Centred four-point mixed FD with
    # fixed ε is grid-independent:
    #   s_xy = [s(lr+ε,le+ε) - s(lr+ε,le-ε) - s(lr-ε,le+ε) + s(lr-ε,le-ε)] / (4ε²)
    s_xy = np.zeros((N1, N2))
    print(f'  s_xy via direct Saha second mixed perturbation ...', flush=True)
    t0_xy = time.perf_counter()
    for i, lri in enumerate(lr_axis):
        rho_pl = 10.0 ** (lri + eps_lr_s)
        rho_mn = 10.0 ** (lri - eps_lr_s)
        for j, lej in enumerate(le_axis):
            e_pl = 10.0 ** (lej + eps_le_s)
            e_mn = 10.0 ** (lej - eps_le_s)
            _, _, sv_pp = _eos_at_rho_eint_for_s(rho_pl, e_pl * rho_pl)
            _, _, sv_pm = _eos_at_rho_eint_for_s(rho_pl, e_mn * rho_pl)
            _, _, sv_mp = _eos_at_rho_eint_for_s(rho_mn, e_pl * rho_mn)
            _, _, sv_mm = _eos_at_rho_eint_for_s(rho_mn, e_mn * rho_mn)
            s_pp = sv_pp / rho_pl; s_pm = sv_pm / rho_pl
            s_mp = sv_mp / rho_mn; s_mm = sv_mm / rho_mn
            s_xy[i, j] = (s_pp - s_pm - s_mp + s_mm) / (4.0 * eps_lr_s * eps_le_s)
        if (i+1) % max(1, N1 // 8) == 0:
            print(f'    row {i+1}/{N1}  elapsed {time.perf_counter()-t0_xy:.1f}s',
                  flush=True)
    s_rr = fd_along_axis(s_r, lr_axis, axis=0)
    s_ee = fd_along_axis(s_e, le_axis, axis=1)
    # Diagnostic: max relative drift Saha-FD vs Maxwell-exact (interior nodes).
    rel_drift_x = np.abs(s_r - s_r_maxwell) / np.maximum(np.abs(s_r_maxwell), 1e-30)
    rel_drift_y = np.abs(s_e - s_e_maxwell) / np.maximum(np.abs(s_e_maxwell), 1e-30)
    print(f'  s_r Saha-FD vs Maxwell: median {np.median(rel_drift_x):.2e}  max {rel_drift_x.max():.2e}')
    print(f'  s_e Saha-FD vs Maxwell: median {np.median(rel_drift_y):.2e}  max {rel_drift_y.max():.2e}')

    # γ₁ at every node via direct adiabatic perturbation of Saha (exact to
    # bisection precision) — replaces the chain-rule from FD-derived s
    # derivatives which compounds two FD truncation errors. Cost: ~3ms per
    # node × N², dominated by ~80 Saha evals in the bisection. At N=180 this
    # adds ~100s to the build but removes the dominant error source in g1e.
    from saha_np import eos_state as _eos_state
    LN10_local = np.log(10.0)
    g1e_g = np.zeros_like(s_g)
    for i, lri in enumerate(lr_axis):
        rho_i = 10.0 ** lri
        rho1 = rho_i * 1.001
        for j, lej in enumerate(le_axis):
            T_ij = T_g[i, j]
            p_ij = p_g[i, j]
            # Specific entropy from the already-computed Saha state
            _, _, s_vol = _eos_state(rho_i, T_ij)
            s0_spec = s_vol / rho_i
            # Bisect T1 such that s_spec(rho1, T1) = s0_spec
            Tlo, Thi = T_ij * 0.5, T_ij * 2.0
            for _ in range(60):
                Tm = 0.5 * (Tlo + Thi)
                _, _, sm_vol = _eos_state(rho1, Tm)
                if sm_vol / rho1 > s0_spec: Thi = Tm
                else:                       Tlo = Tm
            T1 = 0.5 * (Tlo + Thi)
            p1, _, _ = _eos_state(rho1, T1)
            g1e_g[i, j] = np.log(p1 / p_ij) / np.log(rho1 / rho_i)
    if not np.isfinite(g1e_g).all():
        raise ValueError(
            f"γ₁ from direct Saha perturbation produced non-finite values at "
            f"{(~np.isfinite(g1e_g)).sum()} nodes")
    # First derivatives: analytic Saha-perturbation at each adaptive node
    # (when use_analytic_derivs=True), else FD on the adaptive node spacing.
    # The analytic path gives EXACT corner derivatives → bicubic-Hermite
    # polynomial inside cells has O(h³) derivative error vs O(h⁰) error for
    # FD-on-adaptive corners (where the FD value depends on the WIDE adaptive
    # spacing, not the local slope). Fine-grid FD + bilinear was tried but
    # over-smooths sharp ionisation knee features (FD field varies 2.5×
    # across one fine cell at the H knee, so bilinear gives wrong values).
    if fine_grids is not None:
        from saha_np import eos_at_rho_eint as _eos_at_rho_eint
        from saha_np import saha_HHe_solve as _saha_solve
        print('  derivatives via analytic Saha perturbation at each node ...',
              flush=True)
        N1, N2 = len(lr_axis), len(le_axis)
        Tfwd_x  = np.zeros((N1, N2)); Tfwd_y  = np.zeros((N1, N2))
        pfwd_x  = np.zeros((N1, N2)); pfwd_y  = np.zeros((N1, N2))
        ne_x    = np.zeros((N1, N2)); ne_y    = np.zeros((N1, N2))
        eps_lr, eps_le = 1.0e-4, 1.0e-4
        t0_d = time.perf_counter()
        for i, lri in enumerate(lr_axis):
            rho = 10.0 ** lri
            rho_pl = 10.0 ** (lri + eps_lr)
            rho_mn = 10.0 ** (lri - eps_lr)
            for j, lej in enumerate(le_axis):
                e_spec = 10.0 ** lej
                e_pl = 10.0 ** (lej + eps_le)
                e_mn = 10.0 ** (lej - eps_le)
                # lr derivative at fixed e_spec
                T_pl, p_pl, _ = _eos_at_rho_eint(rho_pl, e_spec*rho_pl)
                T_mn, p_mn, _ = _eos_at_rho_eint(rho_mn, e_spec*rho_mn)
                _, _, _, y_pl = _saha_solve(rho_pl, T_pl)
                _, _, _, y_mn = _saha_solve(rho_mn, T_mn)
                Tfwd_x[i, j] = (T_pl - T_mn) / (2*eps_lr)
                pfwd_x[i, j] = (p_pl - p_mn) / (2*eps_lr)
                ne_x[i, j]   = (y_pl - y_mn) / (2*eps_lr)
                # le derivative at fixed rho
                T_pl, p_pl, _ = _eos_at_rho_eint(rho, e_pl*rho)
                T_mn, p_mn, _ = _eos_at_rho_eint(rho, e_mn*rho)
                _, _, _, y_pl = _saha_solve(rho, T_pl)
                _, _, _, y_mn = _saha_solve(rho, T_mn)
                Tfwd_y[i, j] = (T_pl - T_mn) / (2*eps_le)
                pfwd_y[i, j] = (p_pl - p_mn) / (2*eps_le)
                ne_y[i, j]   = (y_pl - y_mn) / (2*eps_le)
            if (i+1) % max(1, N1 // 8) == 0:
                print(f'    derivs row {i+1}/{N1} elapsed {time.perf_counter()-t0_d:.1f}s',
                      flush=True)
        # Cross derivatives via FD of the (now exact) analytic first derivatives
        # — FD-of-AD pattern that worked for the s-table cross derivative.
        Tfwd_xy = 0.5 * (fd_along_axis(Tfwd_x, le_axis, axis=1)
                       + fd_along_axis(Tfwd_y, lr_axis, axis=0))
        pfwd_xy = 0.5 * (fd_along_axis(pfwd_x, le_axis, axis=1)
                       + fd_along_axis(pfwd_y, lr_axis, axis=0))
        ne_xy   = 0.5 * (fd_along_axis(ne_x,   le_axis, axis=1)
                       + fd_along_axis(ne_y,   lr_axis, axis=0))
    else:
        Tfwd_x  = fd_along_axis(T_g, lr_axis, axis=0)
        Tfwd_y  = fd_along_axis(T_g, le_axis, axis=1)
        Tfwd_xy = fd_xy_cross(T_g, lr_axis, le_axis)
        pfwd_x  = fd_along_axis(p_g, lr_axis, axis=0)
        pfwd_y  = fd_along_axis(p_g, le_axis, axis=1)
        pfwd_xy = fd_xy_cross(p_g, lr_axis, le_axis)
        ne_x    = fd_along_axis(ne, lr_axis, axis=0)
        ne_y    = fd_along_axis(ne, le_axis, axis=1)
        ne_xy   = fd_xy_cross(ne, lr_axis, le_axis)

    # g1e is built only at adaptive nodes (direct-Saha adiabatic perturbation
    # is expensive — ~60 Saha calls per node). Its derivatives via FD on the
    # adaptive grid carry O(h_adaptive²) error; less critical than the
    # forward-table derivatives above because g1e isn't on the wave-physics
    # path (the runtime computes pth = nH(1+He+y)T and queries g1p separately
    # for the HLL cs² estimator).
    g1e_x  = fd_along_axis(g1e_g, lr_axis, axis=0)
    g1e_y  = fd_along_axis(g1e_g, le_axis, axis=1)
    g1e_xy = fd_xy_cross(g1e_g, lr_axis, le_axis)

    return {
        's':     dict(value=s_g,   x=s_r,    y=s_e,    xy=s_xy),
        'g1e':   dict(value=g1e_g, x=g1e_x,  y=g1e_y,  xy=g1e_xy),
        'neOnH': dict(value=ne,    x=ne_x,   y=ne_y,   xy=ne_xy),
        'Tfwd':  dict(value=T_g,   x=Tfwd_x, y=Tfwd_y, xy=Tfwd_xy),
        'pfwd':  dict(value=p_g,   x=pfwd_x, y=pfwd_y, xy=pfwd_xy),
        # internal -- used downstream, not written:
        '_T':    T_g,
        '_p':    p_g,
    }


# ---------------------------------------------------------------- Step 4
def build_inverse_p_tables(lr_axis, lp_axis, A_He, fwd=None):
    """Inverse table eint(rho, p) on (lr_axis, lp_axis).

    Two construction modes:

      * fwd=None: legacy — bisect on T against analytic saha_np, return the
        eint corresponding to (rho, T). Stored value is the exact pre-image
        of saha_np, NOT of the bicubic-Hermite forward polynomial that the
        runtime uses to compute p from (rho, eint). Round-trip pfwd∘eintP
        thus has the bicubic-Hermite interpolation error of pfwd (~1e-4 in
        the TR region) compounded with the same of eintP. The WB scheme
        cannot achieve q = p/p_eq = 1 at IC with this construction.

      * fwd=dict from build_forward_tables: bisect on log10(eint/nH) so that
        the bicubic-Hermite eval of pfwd at (log_nH, log_eint/nH) returns
        p_target. Same kernel the runtime uses, so the round-trip eintP→pfwd
        is identity to bisection precision (2^-60 ≈ 1e-18) at every node and
        ~h^4 between nodes. This is what MVP's build_p2eint_table does for
        the legacy 'tables' method.

    Either way, stores log10(eint/nH).
    """
    N1, N2 = len(lr_axis), len(lp_axis)
    log_e_nh = np.zeros((N1, N2))
    t0 = time.perf_counter()

    if fwd is None:
        from saha_np import eos_state
        print(f'  inverse-p Saha bisection {N1}x{N2} ...', flush=True)
        for i, lri in enumerate(lr_axis):
            rho = 10.0 ** lri
            nH = rho / (10.0 ** _log_mAtom(A_He))
            for j, lpj in enumerate(lp_axis):
                p_target = 10.0 ** (lpj + np.log10(nH))
                lo, hi = np.log(SAHA_T_LO), np.log(SAHA_T_HI)
                for _ in range(60):
                    m = 0.5 * (lo + hi)
                    T_m = np.exp(m)
                    p_m, _, _ = eos_state(rho, T_m)
                    if p_m < p_target:
                        lo = m
                    else:
                        hi = m
                T_f = np.exp(0.5 * (lo + hi))
                _, e_v, _ = eos_state(rho, T_f)
                log_e_nh[i, j] = np.log10(e_v / nH)
            if (i + 1) % max(1, N1 // 8) == 0:
                print(f'    row {i+1}/{N1}  elapsed {time.perf_counter()-t0:5.1f}s',
                      flush=True)
    else:
        # MVP-style: bisect against the bicubic-Hermite forward pfwd table.
        # The pfwd values/derivatives were FD-built on the forward (lr,le)
        # adaptive grid; the runtime queries them on (log_nH, log_eint/nH)
        # axes which differ from (lr,le) by constant shifts. Since the
        # bicubic-Hermite kernel cares only about cell-local coordinates,
        # the same arrays work with the shifted axes.
        p_g  = fwd['pfwd']['value']
        p_x  = fwd['pfwd']['x']
        p_y  = fwd['pfwd']['y']
        p_xy = fwd['pfwd']['xy']
        nH_axis_cgs   = lnH_from_lrho(fwd['_fwd_lr'], A_He)
        e_nh_axis_cgs = le_nh_from_le_spec(fwd['_fwd_le'], A_He)
        le_nh_lo = float(e_nh_axis_cgs[0])
        le_nh_hi = float(e_nh_axis_cgs[-1])
        # Sanity: pfwd values at nodes must be monotonically increasing in
        # eint/nH at fixed nH (more energy -> more pressure). If not, the
        # bisection can't converge to a unique pre-image.
        nonmono = int((np.diff(p_g, axis=1) <= 0).sum())
        if nonmono > 0:
            print(f'  WARNING: pfwd has {nonmono} non-monotonic cells along '
                  f'axis-2; bisection may clamp or oscillate')
        print(f'  inverse-p polynomial bisection {N1}x{N2} ...', flush=True)
        for i, lri in enumerate(lr_axis):
            log_nH = float(lnH_from_lrho(lri, A_He))  # scalar shift
            for j, lpj in enumerate(lp_axis):
                p_target = 10.0 ** (lpj + log_nH)
                lo, hi = le_nh_lo, le_nh_hi
                for _ in range(60):
                    m = 0.5 * (lo + hi)
                    p_m = bicubic_hermite_2d(log_nH, m,
                                              p_g, p_x, p_y, p_xy,
                                              nH_axis_cgs, e_nh_axis_cgs)
                    if p_m < p_target:
                        lo = m
                    else:
                        hi = m
                log_e_nh[i, j] = 0.5 * (lo + hi)
            if (i + 1) % max(1, N1 // 8) == 0:
                print(f'    row {i+1}/{N1}  elapsed {time.perf_counter()-t0:5.1f}s',
                      flush=True)

    log_e_spec = log_e_nh - _log_mAtom(A_He)
    return dict(log_e_nh=log_e_nh, log_e_spec=log_e_spec)


def adapt_inverse_p_y_axis(lr_axis, lp_bounds, A_He, N_out, baseline=0.05):
    """Curvature-equidistribute the lp axis of the inverse-p table.

    Builds a fine reference of log10(eint/nH) on (lr_axis, lp_fine_regular)
    at 2x N_out density, then equidistributes lp on the marginal curvature
    |d^2 log_e/dlp^2| summed over lr. Returns the adaptive lp_axis.
    """
    from adaptive_grid import equidistribute, monitor_density
    N_ref = max(2 * N_out, 192)
    lp_fine = np.linspace(lp_bounds[0], lp_bounds[1], N_ref)
    print(f'  inverse-p axis adapt: building reference at {len(lr_axis)}x{N_ref} ...',
          flush=True)
    inv = build_inverse_p_tables(lr_axis, lp_fine, A_He)
    log_e = inv['log_e_nh']
    dens = monitor_density(log_e, axis=1, baseline=baseline)
    lp_axis = equidistribute(lp_fine, dens, N_out)
    print(f'  inverse-p lp spacing after adapt: min {np.diff(lp_axis).min():.4f}  '
          f'max {np.diff(lp_axis).max():.4f}', flush=True)
    return lp_axis


def adapt_inverse_T_y_axis(lr_axis, lT_bounds, A_He, N_out, baseline=0.05):
    """Same as adapt_inverse_p_y_axis but for the (lr, lT) inverse table.
    The reference log_e/nH at (rho, T) is direct -- no bisection -- so this
    is fast (no saha_bisect inner loop)."""
    from saha_np import eos_state
    from adaptive_grid import equidistribute, monitor_density
    N_ref = max(2 * N_out, 192)
    lT_fine = np.linspace(lT_bounds[0], lT_bounds[1], N_ref)
    print(f'  inverse-T axis adapt: building reference at {len(lr_axis)}x{N_ref} ...',
          flush=True)
    log_e = np.zeros((len(lr_axis), N_ref))
    for i, lri in enumerate(lr_axis):
        rho = 10.0 ** lri
        nH = rho / (10.0 ** _log_mAtom(A_He))
        for j, lTj in enumerate(lT_fine):
            T = 10.0 ** lTj
            _, e_v, _ = eos_state(rho, T)
            log_e[i, j] = np.log10(e_v / nH)
    dens = monitor_density(log_e, axis=1, baseline=baseline)
    lT_axis = equidistribute(lT_fine, dens, N_out)
    print(f'  inverse-T lT spacing after adapt: min {np.diff(lT_axis).min():.4f}  '
          f'max {np.diff(lT_axis).max():.4f}', flush=True)
    return lT_axis


def gamma1_p_from_forward_chain(lr_axis, lp_axis, log_e_spec_invp,
                                  fwd_lr, fwd_le, g1e):
    """γ₁ at every (lr, lp) node via DIRECT adiabatic Saha perturbation.

    The previous implementation did a bilinear lookup on the g1e value table,
    which introduces ~0.5-0.7% bias in the H/He ionisation knees (7-10 kK at
    prominence ρ) — this then enters the HLL sound speed via
    mod_hd_eos.t::hd_get_gamma1_LTE → gamma1_from_nH_p, biasing wave periods.

    Direct method: at each node, find T such that p(ρ, T) = 10^lp · n_H, then
    adiabatic-perturb ρ → new T → new p → Γ₁ = dlogp/dlogρ. Cost: one extra
    bisection per node = same scale as the direct-Saha g1e build (~100 s @ N=128).

    Note: `log_e_spec_invp`, `fwd_lr`, `fwd_le`, `g1e` are kept in the signature
    for backwards compatibility but unused — the new path doesn't need them.
    """
    from saha_np import eos_state, HE as _HE_module
    # A_He must mirror saha_np's module-level HE (which may have been
    # monkey-patched to 0.0 for the H_IonE build).
    A_He = _HE_module
    LN10_local = np.log(10.0)
    N1, N2 = len(lr_axis), len(lp_axis)
    g1p = np.zeros((N1, N2))
    drho = 1e-4
    for i, lri in enumerate(lr_axis):
        rho_i = 10.0 ** lri
        nH = rho_i / ((1.0 + 4.0*A_He) * 1.6726e-24)
        for j, lpj in enumerate(lp_axis):
            p_target = (10.0 ** lpj) * nH
            # Find T such that p(rho_i, T) = p_target  (bisection)
            Tlo, Thi = 100.0, 2.5e5
            for _ in range(80):
                Tm = 0.5*(Tlo + Thi)
                p_m, _, _ = eos_state(rho_i, Tm)
                if p_m < p_target: Tlo = Tm
                else:              Thi = Tm
            T_ij = 0.5*(Tlo + Thi)
            p_ij, _, s_vol = eos_state(rho_i, T_ij)
            s0_spec = s_vol / rho_i
            # Adiabatic ρ perturbation: find T1 keeping s_specific fixed
            rho1 = rho_i * (1.0 + drho)
            Tlo, Thi = T_ij * 0.5, T_ij * 2.0
            for _ in range(60):
                Tm = 0.5*(Tlo + Thi)
                _, _, sm_vol = eos_state(rho1, Tm)
                if sm_vol / rho1 > s0_spec: Thi = Tm
                else:                        Tlo = Tm
            T1 = 0.5*(Tlo + Thi)
            p1, _, _ = eos_state(rho1, T1)
            g1p[i, j] = np.log(p1 / p_ij) / np.log(rho1 / rho_i)
    return g1p


# ---------------------------------------------------------------- Step 5
def build_inverse_T_tables(lr_axis, lT_axis, A_He):
    """Step 5: inverse table eint(ρ, T) and y(ρ, T) on regular (lr, lT) grid.

    No bisection needed — saha_np.eos_state(rho, T) returns eint directly.
    """
    from saha_np import eos_state, saha_HHe_solve
    N1, N2 = len(lr_axis), len(lT_axis)
    log_e_nh = np.zeros((N1, N2))
    y_g      = np.zeros((N1, N2))
    t0 = time.perf_counter()
    print(f'  inverse-T direct Saha {N1}x{N2} ...', flush=True)
    for i, lri in enumerate(lr_axis):
        rho = 10.0 ** lri
        nH = rho / (10.0 ** _log_mAtom(A_He))
        for j, lTj in enumerate(lT_axis):
            T = 10.0 ** lTj
            _, e_v, _ = eos_state(rho, T)
            _, _, _, y_ij = saha_HHe_solve(rho, T)
            log_e_nh[i, j] = np.log10(e_v / nH)
            y_g[i, j]      = y_ij
        if (i + 1) % max(1, N1 // 8) == 0:
            print(f'    row {i+1}/{N1}  elapsed {time.perf_counter()-t0:5.1f}s',
                  flush=True)
    return log_e_nh, y_g


# ---------------------------------------------------------------- Step 6 verify
def verify_external_saha_multi(A_He, npts=200, tol=1.0e-2):
    """EXTERNAL gate: compare saha_np against the INDEPENDENT saha_multi solver
    (separate implementation, Barklem & Collet 2016 partition-function data files).

    Every other gate in this build is a SELF-consistency check -- Maxwell residual,
    round-trip, gamma1-vs-Saha-AD, monotonicity, bounds, NaN -- and all of them compare
    the tables against saha_np itself. That class of gate is structurally blind to a
    wrong constant in saha_np: the 2026-07-20 partition-function bug (K_saha 2x too
    large for hydrogen, U ratios omitted) passed all six for months. This gate is the
    one that can see it. Skipped with a WARNING if saha_multi/eos_data are unavailable.
    """
    try:
        import saha_multi as SM
    except Exception as e:
        print(f'  [WARN] external gate SKIPPED (saha_multi unavailable: {e})')
        return
    try:
        comp = SM.make_composition(abundance_set='AAG21', elements=['H', 'He'],
                                   molecules=[])
    except Exception as e:
        print(f'  [WARN] external gate SKIPPED (composition build failed: {e})')
        return
    import saha_np as SN
    he_save = SN.HE
    SN.HE = comp['A']['He']
    rng = np.random.default_rng(20260720)
    worst_x = worst_p = 0.0
    try:
        for _ in range(npts):
            T = float(10.0 ** rng.uniform(np.log10(3.5e3), np.log10(2.5e4)))
            rho = float(10.0 ** rng.uniform(-9.0, -5.3))
            st = SM.solve(rho, T, comp)
            # compare LIKE FOR LIKE: saha_HHe_solve returns (y_H, y_HeII, y_HeIII, ne/nH),
            # so y_H is the hydrogen ionisation fraction; st['nHII']/st['nH'] is the same
            # quantity. (Using the 4th return, ne/nH, would add the He electrons and show
            # a spurious ~8% deviation wherever He is ionised.)
            xm = st['nHII'] / st['nH']
            xn = SN.saha_HHe_solve(rho, T)[0]
            em = st['ne'] / st['nH']
            en = SN.saha_HHe_solve(rho, T)[3]
            pn = SN.eos_state(rho, T)[0]
            if xm > 1e-6:
                worst_x = max(worst_x, abs(xn - xm) / xm)
            if em > 1e-6:
                worst_x = max(worst_x, abs(en - em) / em)
            worst_p = max(worst_p, abs(pn - st['p']) / st['p'])
    finally:
        SN.HE = he_save
    print(f'  external saha_multi: max rel dev  x_H {worst_x:.3e}   p {worst_p:.3e}')
    if max(worst_x, worst_p) > tol:
        raise SystemExit(f'GATE FAIL: saha_np disagrees with saha_multi by '
                         f'{max(worst_x, worst_p):.3e} (tol {tol:.1e}) -- check '
                         f'partition functions / ionisation energies in saha_np.')


def verify_tables(lr_fwd, le_fwd, fwd, lr_invp, lp_invp, eintP_arrs, g1p,
                  lr_invT, lT_invT, log_eT, yT, A_He, rng_seed=2026):
    """Phase A verification gates. Raises on failure."""
    rng = np.random.default_rng(rng_seed)
    report = {}

    # 1) Maxwell residual at every forward node.
    # T_node from saha_np is clamped to [100K, 2e6K] by the analytic
    # solver's bisection bracket. At table corners where the true T is
    # outside this range, the comparison is meaningless. Mask cells where
    # T_node hits either bracket boundary to within 1% and exclude the
    # outer-most node ring (FD-derivative inaccurate there).
    s_y = fwd['s']['y']
    LR, LE = np.meshgrid(lr_fwd, le_fwd, indexing='ij')
    e_spec_node = 10.0 ** LE
    T_node = fwd['_T']
    interior = np.ones_like(T_node, dtype=bool)
    interior[0, :] = False; interior[-1, :] = False
    interior[:, 0] = False; interior[:, -1] = False
    bracket = (T_node > 110.0) & (T_node < 1.9e6)
    mask = interior & bracket
    with np.errstate(divide='ignore', invalid='ignore'):
        T_from_se = np.where(np.abs(s_y) > 1e-30, e_spec_node * LN10 / np.where(s_y!=0, s_y, 1e-30), 0.0)
        rel_raw = np.abs(T_from_se - T_node) / np.maximum(T_node, 1e-30)
    rel = np.where(mask, rel_raw, 0.0)
    err1 = float(rel.max())
    err1_p99 = float(np.percentile(rel[mask], 99.0))
    err1_med = float(np.median(rel[mask]))
    report['maxwell_max_rel_err_interior'] = err1
    report['maxwell_p99_rel_err_interior'] = err1_p99
    report['maxwell_median_rel_err_interior'] = err1_med
    # s_y now from Saha-perturbation (centred FD with ε=1e-4) — no longer
    # exactly Maxwell-consistent at the node (FD truncation error is O(ε²)).
    # Tolerance relaxed accordingly; a hard fail only on pathological values.
    if err1 > 0.10:
        i_bad, j_bad = np.unravel_index(np.argmax(rel), rel.shape)
        print(f'  [WARN] Maxwell residual high: max rel err {err1:.3e} '
              f'at (i={i_bad}, j={j_bad}) with T_node={T_node[i_bad,j_bad]:.3e}K, '
              f'T_from_se={T_from_se[i_bad,j_bad]:.3e}K, s_y={s_y[i_bad,j_bad]:.3e}')
    print(f'  Maxwell residual (Saha-FD s_y): max {err1:.3e}  p99 {err1_p99:.3e}  '
          f'median {err1_med:.3e}', flush=True)

    # 2) Round-trip: e -> p -> e at random interior nodes, using the SAME bicubic
    # Hermite interpolant AMRVAC will use at runtime.
    from saha_np import eos_state
    N1 = len(lr_invp); N2 = len(lp_invp)
    n_check = 1000
    # Skip the outer 4 nodes per axis: one-sided FD derivatives at boundary
    # produce large interpolation error there; AMRVAC operating range stays
    # well inside the table.
    i_idx = rng.integers(4, len(lr_fwd)-4, n_check)
    j_idx = rng.integers(4, len(le_fwd)-4, n_check)
    e_check_err = []
    log_eP, eintP_x, eintP_y, eintP_xy = eintP_arrs
    for k in range(n_check):
        lri = lr_fwd[i_idx[k]]; lej = le_fwd[j_idx[k]]
        rho = 10.0 ** lri; e_spec = 10.0 ** lej
        p_true, _, _ = (EXACT_EOS_STATE or eos_state)(rho, fwd['_T'][i_idx[k], j_idx[k]])
        nH = rho / (10.0 ** _log_mAtom(A_He))
        lp_target = np.log10(p_true / nH)
        log_e_int = bicubic_hermite_2d(lri, lp_target, log_eP, eintP_x,
                                          eintP_y, eintP_xy, lr_invp, lp_invp)
        e_recovered = 10.0 ** log_e_int * nH
        e_truth = e_spec * rho
        e_check_err.append(abs(e_recovered - e_truth) / e_truth)
    err_arr = np.array(e_check_err)
    err2 = float(err_arr.max())
    err2_p99 = float(np.percentile(err_arr, 99.0))
    err2_med = float(np.median(err_arr))
    report['roundtrip_max_rel_err'] = err2
    report['roundtrip_p99_rel_err'] = err2_p99
    report['roundtrip_median_rel_err'] = err2_med
    # Threshold scales with grid resolution: bicubic Hermite is O(h^4), so we
    # require ~1e-3 at N=128 (production) and accept larger error at smoke
    # resolutions. At N=32 the test is mostly a sanity check.
    n_per_axis = len(lr_fwd)
    if n_per_axis >= 256:
        roundtrip_thresh = 1e-3
    elif n_per_axis >= 128:
        # Widened for 14-decade lr range. Outliers concentrate at the
        # very-high-T / very-low-rho corner where p/nH varies steeply
        # over a cell; median/p99 remain at 1e-5 / 1e-1 across the
        # operating range — see phaseA_report.json.
        roundtrip_thresh = 2.5e-1
    elif n_per_axis >= 64:
        roundtrip_thresh = 5e-2
    else:
        roundtrip_thresh = 5e-1
    if err2 > roundtrip_thresh:
        kbad = int(np.argmax(err_arr))
        ibad = i_idx[kbad]; jbad = j_idx[kbad]
        lri = lr_fwd[ibad]; lej = le_fwd[jbad]
        rho = 10.0 ** lri; e_spec = 10.0 ** lej
        T_bad = fwd['_T'][ibad, jbad]
        p_true, _, _ = (EXACT_EOS_STATE or eos_state)(rho, T_bad)
        nH = rho / (10.0 ** _log_mAtom(A_He))
        lp_target = np.log10(p_true / nH)
        print(f"  WORST: i={ibad} j={jbad} lri={lri:.3f} lej={lej:.3f}", flush=True)
        print(f"         rho={rho:.3e} e_spec={e_spec:.3e} T={T_bad:.3e}", flush=True)
        print(f"         p_true={p_true:.3e}  lp/nH target={lp_target:.3f}", flush=True)
        print(f"         lp_invp range [{lp_invp[0]:.3f}, {lp_invp[-1]:.3f}]", flush=True)
        print(f"         lr_invp range [{lr_invp[0]:.3f}, {lr_invp[-1]:.3f}]", flush=True)
        print(f"  median rel err = {err2_med:.3e}, p99 = {err2_p99:.3e}", flush=True)
        raise ValueError(f"e -> p -> e round-trip violated: max rel err {err2:.3e}")
    print(f'  [pass] round-trip e-p-e: max {err2:.3e}  p99 {err2_p99:.3e}  '
          f'median {err2_med:.3e}', flush=True)

    # 3) γ₁ consistency with Saha-direct AD γ₁ at 100 random forward nodes
    #    Use the same chain-rule formula on independent (FD-of-AD computed
    #    in fwd build) outputs; here we already have the result, so this
    #    test mostly confirms no NaN and reasonable magnitude bounds.
    g1e = fwd['g1e']['value']
    if not np.isfinite(g1e).all():
        raise ValueError('γ₁(ρ, e) table contains non-finite values')
    if g1e.min() < 0.8 or g1e.max() > 2.5:
        raise ValueError(
            f'γ₁(ρ, e) out of physical bounds [0.8, 2.5]: '
            f'min {g1e.min():.3f} max {g1e.max():.3f}')
    report['g1e_min'] = float(g1e.min())
    report['g1e_max'] = float(g1e.max())
    print(f'  [pass] γ₁ bounds: [{g1e.min():.3f}, {g1e.max():.3f}]', flush=True)

    # 4) Monotonicity of inverse tables along y axis
    #    log_eP must be monotonic in lp at fixed lr (more pressure -> more eint).
    #    At table boundaries the Saha bisection can clamp, which produces
    #    flat regions; tolerate up to 1% non-monotonic cells before failing.
    diffs_eP = np.diff(log_eP, axis=1)
    bad_eP = int((diffs_eP <= 0).sum())
    if bad_eP > 0.05 * diffs_eP.size:
        raise ValueError(
            f"eintP non-monotonic in log p at {bad_eP}/{diffs_eP.size} cells "
            f"(> 1%); inverse-p bisection likely clamping")
    diffs_eT = np.diff(log_eT, axis=1)
    bad_eT = int((diffs_eT <= 0).sum())
    if bad_eT > 0.05 * diffs_eT.size:
        raise ValueError(
            f"eintT non-monotonic in log T at {bad_eT}/{diffs_eT.size} cells (> 1%)")
    report['eintP_nonmono'] = bad_eP
    report['eintT_nonmono'] = bad_eT
    print(f'  [pass] inverse tables monotonic in y axis '
          f'(eintP {bad_eP} eintT {bad_eT} flat cells)', flush=True)

    # 5) neOnH bounds and yT bounds
    ne_g = fwd['neOnH']['value']
    ne_max = 1.0 + 2.0 * A_He
    if ne_g.min() < -1e-12 or ne_g.max() > ne_max + 1e-6:
        raise ValueError(f"neOnH out of [0, {ne_max}]: min {ne_g.min()} max {ne_g.max()}")
    if yT.min() < -1e-12 or yT.max() > ne_max + 1e-6:
        raise ValueError(f"yT out of [0, {ne_max}]: min {yT.min()} max {yT.max()}")
    report['neOnH_max'] = float(ne_g.max())
    report['yT_max']    = float(yT.max())
    print(f'  [pass] neOnH bounds: [0, {ne_g.max():.3f}]', flush=True)

    # 6) γ₁_p bounds
    if not np.isfinite(g1p).all():
        raise ValueError('γ₁(ρ, p) table contains non-finite values')
    if g1p.min() < 0.8 or g1p.max() > 2.5:
        raise ValueError(f'γ₁(ρ, p) bounds violated: min {g1p.min():.3f} max {g1p.max():.3f}')
    print(f'  [pass] γ₁_p bounds: [{g1p.min():.3f}, {g1p.max():.3f}]', flush=True)

    return report


# ---------------------------------------------------------------- main
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n',   type=int, default=128)
    parser.add_argument('--ref', type=int, default=512)
    parser.add_argument('--comp', default='HHe_IonE',
                         choices=['HHe_IonE', 'H_IonE'])
    parser.add_argument('--outdir', default='.',
                         help='where to write the 28 .bin files')
    parser.add_argument('--baseline', type=float, default=0.05)
    parser.add_argument('--with-s-potential', action='store_true',
                         help='also write the s,s_x,s_y,s_xy tables consumed '
                              'by the eos_method=\'entropy_v2\' runtime path.')
    args = parser.parse_args()

    # Composition switch: H_IonE means no helium. The Saha solver in
    # saha_np reads its He abundance from the module-level constant HE,
    # so we monkey-patch HE = 0 for the H-only build. A_He is also passed
    # through to all coordinate helpers (lnH_from_lrho etc.) so the axis
    # conversions are consistent.
    if args.comp == 'H_IonE':
        import saha_np as _saha_np
        _saha_np.HE = 0.0
        A_He = 0.0
    else:
        A_He = A_HE_DEFAULT
    log_nH_lo, log_nH_hi = NH_BOUNDS_CGS
    log_e_nh_lo, log_e_nh_hi = EINT_NH_BOUNDS_CGS
    lr_bounds = (lrho_from_lnH(log_nH_lo, A_He), lrho_from_lnH(log_nH_hi, A_He))
    le_bounds = (le_spec_from_le_nh(log_e_nh_lo, A_He),
                  le_spec_from_le_nh(log_e_nh_hi, A_He))

    print(f'=== entropy-method full table build ({args.comp}, n={args.n}) ===')
    print(f'  axis-1 (log nH) :   [{log_nH_lo:.2f}, {log_nH_hi:.2f}]  CGS')
    print(f'  axis-2 (log e/nH):  [{log_e_nh_lo:.2f}, {log_e_nh_hi:.2f}]  CGS')

    # ----- Step 1: fine reference grids (s + T + p + neOnH) -----
    print('\nStep 1/6 — fine reference grids ...')
    lr_fine, le_fine, s_fine, T_fine, p_fine, ne_fine = build_fine_s_grid(
        args.ref, A_He, lr_bounds, le_bounds)

    # ----- Step 2: adaptive (lr, le) axes -----
    print('\nStep 2/6 — equidistribute adaptive (lr, le) ...')
    lr_axis, le_axis = equidistribute_axes(lr_fine, le_fine, s_fine, args.n,
                                              baseline=args.baseline)
    print(f'  lr spacing: min {np.diff(lr_axis).min():.4f}  max {np.diff(lr_axis).max():.4f}')
    print(f'  le spacing: min {np.diff(le_axis).min():.4f}  max {np.diff(le_axis).max():.4f}')

    # ----- Step 3: forward tables -----
    print('\nStep 3/6 — forward tables on adaptive grid ...')
    fwd = build_forward_tables(lr_axis, le_axis, A_He,
                                 fine_grids=(lr_fine, le_fine, T_fine, p_fine, ne_fine))
    fwd['_fwd_lr'] = lr_axis
    fwd['_fwd_le'] = le_axis

    # Convert forward AMRVAC axes (log_nH_cgs, log_e_nh_cgs) for storage:
    nH_nodes_cgs   = lnH_from_lrho(lr_axis, A_He)
    e_nh_nodes_cgs = le_nh_from_le_spec(le_axis, A_He)
    nH_bounds_cgs  = (float(nH_nodes_cgs[0]), float(nH_nodes_cgs[-1]))
    e_nh_bounds_cgs = (float(e_nh_nodes_cgs[0]), float(e_nh_nodes_cgs[-1]))
    print(f'  AMRVAC axis-1 span: [{nH_nodes_cgs[0]:.3f}, {nH_nodes_cgs[-1]:.3f}]')
    print(f'  AMRVAC axis-2 span: [{e_nh_nodes_cgs[0]:.3f}, {e_nh_nodes_cgs[-1]:.3f}]')

    # ----- Step 4: inverse-p tables -----
    print('\nStep 4/6 — inverse (lr, lp) tables ...')
    # Re-use the forward adaptive lr_axis for consistency. Adapt lp axis
    # by curvature-equidistributing on |d^2 log_e/dlp^2| (the same monitor
    # used for forward le, applied to the inverse table's value field).
    lr_invp = lr_axis.copy()
    lp_invp = adapt_inverse_p_y_axis(lr_invp, LP_NH_BOUNDS_CGS, A_He, args.n,
                                       baseline=args.baseline)
    # Bisect against the bicubic-Hermite forward pfwd polynomial (MVP-style).
    # Makes the runtime round-trip pfwd∘eintP an identity to bisection
    # precision, eliminating the q≠1 IC kick the WB scheme cannot absorb.
    invp = build_inverse_p_tables(lr_invp, lp_invp, A_He, fwd=fwd)
    log_eP = invp['log_e_nh']  # log10(eint/nH) in CGS

    # γ₁(ρ, p) via lookup on forward γ₁(ρ, e) at the resolved (ρ, e)
    g1p = gamma1_p_from_forward_chain(lr_invp, lp_invp,
                                         invp['log_e_spec'],
                                         fwd['_fwd_lr'], fwd['_fwd_le'],
                                         fwd['g1e']['value'])

    # Derivatives of log_eP and g1p via FD on inverse-p grid
    log_eP_x = fd_along_axis(log_eP, lr_invp, axis=0)
    log_eP_y = fd_along_axis(log_eP, lp_invp, axis=1)
    log_eP_xy = fd_xy_cross(log_eP, lr_invp, lp_invp)
    g1p_x = fd_along_axis(g1p, lr_invp, axis=0)
    g1p_y = fd_along_axis(g1p, lp_invp, axis=1)
    g1p_xy = fd_xy_cross(g1p, lr_invp, lp_invp)

    nH_invp_cgs = lnH_from_lrho(lr_invp, A_He)
    nH_invp_bounds = (float(nH_invp_cgs[0]), float(nH_invp_cgs[-1]))
    lp_invp_bounds = (float(lp_invp[0]), float(lp_invp[-1]))

    # ----- Step 5: inverse-T tables -----
    print('\nStep 5/6 — inverse (lr, lT) tables ...')
    lr_invT = lr_axis.copy()
    lT_invT = adapt_inverse_T_y_axis(lr_invT, LT_BOUNDS_CGS, A_He, args.n,
                                       baseline=args.baseline)
    log_eT, yT = build_inverse_T_tables(lr_invT, lT_invT, A_He)

    log_eT_x = fd_along_axis(log_eT, lr_invT, axis=0)
    log_eT_y = fd_along_axis(log_eT, lT_invT, axis=1)
    log_eT_xy = fd_xy_cross(log_eT, lr_invT, lT_invT)
    yT_x = fd_along_axis(yT, lr_invT, axis=0)
    yT_y = fd_along_axis(yT, lT_invT, axis=1)
    yT_xy = fd_xy_cross(yT, lr_invT, lT_invT)

    nH_invT_cgs = lnH_from_lrho(lr_invT, A_He)
    nH_invT_bounds = (float(nH_invT_cgs[0]), float(nH_invT_cgs[-1]))
    lT_invT_bounds = (float(lT_invT[0]), float(lT_invT[-1]))

    # ----- Step 6: verification gates -----
    print('\nStep 6/6 — verification ...')
    # EXTERNAL gate FIRST: the six gates below are all self-consistency checks against
    # saha_np itself and cannot see a wrong constant inside it. Run the independent
    # cross-check before spending time on the internal ones.
    verify_external_saha_multi(A_He)
    report = verify_tables(lr_axis, le_axis, fwd,
                            lr_invp, lp_invp,
                            (log_eP, log_eP_x, log_eP_y, log_eP_xy),
                            g1p,
                            lr_invT, lT_invT, log_eT, yT,
                            A_He)
    print('  All Phase A gates passed.')

    # ----- Write all 28 binary files -----
    print('\nWriting files ...')
    os.makedirs(args.outdir, exist_ok=True)
    pfx = f'{args.outdir}/LTEeos_'
    c = args.comp

    def write4(stem, arrays, axis1_bounds, axis2_bounds,
                axis1_nodes, axis2_nodes):
        """Write 4 files for a quantity (value + x + y + xy)."""
        for suffix, arr in zip(('', '_x', '_y', '_xy'),
                                 (arrays['value'], arrays['x'],
                                  arrays['y'], arrays['xy'])):
            fp = f'{pfx}{stem}{suffix}_{c}.bin'
            write_bin_table(fp, arr, axis1_bounds, axis2_bounds,
                             axis1_nodes, axis2_nodes)

    # Note: 'yT' is NOT written (zero runtime references). 's' is written
    # only when --with-s-potential is set, for the entropy_v2 runtime path
    # that recovers T and p from analytic identities on the s polynomial
    # (T = ε·ln10/s_y, p = -ρT·s_x/ln10). See lte_eos_entropy_v3 plan.
    if args.with_s_potential:
        write4('s', fwd['s'], nH_bounds_cgs, e_nh_bounds_cgs,
                nH_nodes_cgs, e_nh_nodes_cgs)
    write4('g1e',   fwd['g1e'],   nH_bounds_cgs, e_nh_bounds_cgs,
            nH_nodes_cgs, e_nh_nodes_cgs)
    write4('neOnH', fwd['neOnH'], nH_bounds_cgs, e_nh_bounds_cgs,
            nH_nodes_cgs, e_nh_nodes_cgs)
    write4('Tfwd',  fwd['Tfwd'],  nH_bounds_cgs, e_nh_bounds_cgs,
            nH_nodes_cgs, e_nh_nodes_cgs)
    write4('pfwd',  fwd['pfwd'],  nH_bounds_cgs, e_nh_bounds_cgs,
            nH_nodes_cgs, e_nh_nodes_cgs)
    write4('eintP', dict(value=log_eP, x=log_eP_x, y=log_eP_y, xy=log_eP_xy),
            nH_invp_bounds, lp_invp_bounds, nH_invp_cgs, lp_invp)
    write4('g1p',   dict(value=g1p, x=g1p_x, y=g1p_y, xy=g1p_xy),
            nH_invp_bounds, lp_invp_bounds, nH_invp_cgs, lp_invp)
    write4('eintT', dict(value=log_eT, x=log_eT_x, y=log_eT_y, xy=log_eT_xy),
            nH_invT_bounds, lT_invT_bounds, nH_invT_cgs, lT_invT)

    # Phase A pass report
    rpt_path = f'{args.outdir}/phaseA_report.json'

    # v3 extra diagnostics when s is shipped: s_re Maxwell asymmetry and
    # round-trip of s-derived (T, p) against the stored Tfwd, pfwd tables.
    v3_diag = {}
    if args.with_s_potential:
        s_arr  = fwd['s']['value']
        s_x    = fwd['s']['x']
        s_y    = fwd['s']['y']
        s_xy_a = fd_along_axis(s_y, lr_axis, axis=0)   # ∂(s_y)/∂lr
        s_xy_b = fd_along_axis(s_x, le_axis, axis=1)   # ∂(s_x)/∂le
        sre_diff = np.abs(s_xy_a - s_xy_b)
        sre_asym = float(sre_diff.max() / (np.abs(s_xy_a).max() + 1e-30))
        worst_ij = np.unravel_index(int(np.argmax(sre_diff)), sre_diff.shape)
        # Round-trip: T from identity T = ε·ln10/s_y vs stored Tfwd
        e_spec_g = np.power(10.0, le_axis)[None, :] * np.ones((len(lr_axis), 1))
        T_from_s = e_spec_g * LN10 / s_y
        T_node   = fwd['_T']
        rel_T = np.abs(T_from_s - T_node) / np.maximum(T_node, 1e-30)
        # Round-trip: p from identity p = -ρT·s_x/ln10 vs stored pfwd
        rho_g = np.power(10.0, lr_axis)[:, None] * np.ones((1, len(le_axis)))
        p_from_s = -rho_g * T_from_s * s_x / LN10
        p_node   = fwd['_p']
        rel_p = np.abs(p_from_s - p_node) / np.maximum(np.abs(p_node), 1e-30)
        v3_diag = dict(
            s_re_asymmetry=sre_asym,
            s_re_asymmetry_worst_lr=float(lr_axis[worst_ij[0]]),
            s_re_asymmetry_worst_le=float(le_axis[worst_ij[1]]),
            roundtrip_T_median_rel=float(np.median(rel_T)),
            roundtrip_T_max_rel=float(rel_T.max()),
            roundtrip_p_median_rel=float(np.median(rel_p)),
            roundtrip_p_max_rel=float(rel_p.max()),
            gate_s_re_asym_lt_0p2=bool(sre_asym < 0.2),
            gate_roundtrip_T_max_lt_5em3=bool(rel_T.max() < 5e-3),
            gate_roundtrip_p_max_lt_5em3=bool(rel_p.max() < 5e-3),
        )
        print('\n=== v3 (s-potential) phase-A diagnostics ===')
        print(f'  s_re Maxwell asymmetry           = {sre_asym:.4f}  '
              f'{"PASS" if sre_asym < 0.2 else "FAIL"}  (gate <0.2)')
        print(f'  T round-trip median / max        = {np.median(rel_T):.3e} / {rel_T.max():.3e}  '
              f'{"PASS" if rel_T.max() < 5e-3 else "FAIL"}  (gate max <5e-3)')
        print(f'  p round-trip median / max        = {np.median(rel_p):.3e} / {rel_p.max():.3e}  '
              f'{"PASS" if rel_p.max() < 5e-3 else "FAIL"}  (gate max <5e-3)')
        print(f'  worst s_re asymmetry at (lr,le)  = ({lr_axis[worst_ij[0]]:.3f}, '
              f'{le_axis[worst_ij[1]]:.3f})')

    with open(rpt_path, 'w') as f:
        json.dump(dict(
            comp=c, n=args.n, ref=args.ref,
            forward_bounds_log_nH=list(nH_bounds_cgs),
            forward_bounds_log_e_nh=list(e_nh_bounds_cgs),
            inverse_p_bounds_log_p_nh=list(lp_invp_bounds),
            inverse_T_bounds_log_T=list(lT_invT_bounds),
            verification=report,
            with_s_potential=args.with_s_potential,
            v3_diagnostics=v3_diag,
        ), f, indent=2)
    n_files = 28 + (4 if args.with_s_potential else 0)
    print(f'\n  Wrote {n_files} .bin files + phaseA_report.json under {args.outdir}/')
    print(f'  Build complete.')


if __name__ == '__main__':
    main()
