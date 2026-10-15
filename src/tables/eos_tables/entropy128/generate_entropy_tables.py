"""Build the bicubic-Hermite EoS tables of the 'entropy' method of mod_eos_lte_entropy

Every runtime EoS call is one bicubic-Hermite evaluation, so each of the seven tabulated
quantities is stored as four arrays: the value and its x, y and xy derivatives. That makes 28
binary files per composition:
  forward, adaptive (log nH, log eint/nH)  LTEeos_{Tfwd,pfwd,neOnH,g1e}{,_x,_y,_xy}_<comp>.bin
  inverse, (log nH, log p/nH)             LTEeos_{eintP,g1p}{,_x,_y,_xy}_<comp>.bin
  inverse, (log nH, log T)                LTEeos_eintT{,_x,_y,_xy}_<comp>.bin

Node values come from the H+He Saha solver in saha_np. First derivatives are exact, from the
implicit differentiation of the Saha charge balance. The eintP values are pre-images under the
pfwd polynomial, so their derivatives follow from the inverse-function rule on that polynomial.
Mixed derivatives are central differences of exact first derivatives, with step STEP_FD. Gamma_1
has no exact derivative here and is differenced from exact values, with the four-point stencil
of step STEP_MIXED for the mixed term. No derivative uses neighbouring table nodes

Checks 1 to 5 run in verify_tables before any file is written, and check 6 as each is written:
  1. the Maxwell relation e_spec T_x = (T p_y - p T_y)/rho holds at every forward node
  2. the e -> p -> e round trip through pfwd and eintP closes to 1e-3 at 1000 random interior
     forward nodes (N >= 128), excluding nodes whose pressure lies off the inverse-p axis
  3. both Gamma_1 tables are finite and inside [0.8, 2.5]
  4. eintP and eintT rise along y in all but 5 per cent of cells
  5. neOnH lies in [0, 1 + 2 a_He] at every node
  6. no output array holds a NaN or an inf

Run:
    python3 generate_entropy_tables.py [--n 128] [--ref 512] [--comp HHe_IonE]
                                       [--outdir <path>] [--baseline 0.05]
"""
import argparse
import json
import os
import struct
import sys
import time

import numpy as np

# The sibling modules are importable from any working directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from adaptive_grid import equidistribute, monitor_density
from saha_np import A_HE, M_P, SahaHHe

LN10 = np.log(10.0)

# Table bounds, CGS. The FALC chromosphere and corona and the Ballester test states lie inside
#   log nH [10.02, 17.12], log eint/nH [-11.99, -10.09],
#   log p/nH [-12.17, -10.50], log T [3.65, 5.00]
# and the bounds add a modest margin to that envelope. The coldest node (about 700 K) lies inside
# the temperature bracket of saha_np, which raises rather than clamps at its edges
LOG_NH_BOUNDS = (5.0, 20.0)        # up to 20 to give Coleman 5/6 pureH (log nH 19.08) headroom
LOG_E_NH_BOUNDS = (-12.8, -9.0)    # log eint/nH (erg per H atom)
LOG_P_NH_BOUNDS = (-12.97, -9.4)   # floor at p/nH of the coldest forward state, (2/3) 10^-12.8,
                                   # so no eintP pre-image is clamped
LOG_T_BOUNDS = (2.84, 6.31)        # log T (K), as in the uniform256 tables
STEP_FD = 2.0e-5     # log10 step of a central difference of an exact first derivative or value
STEP_MIXED = 3.0e-4  # log10 step of the four-point mixed difference of an exact value

COMPOSITIONS = {'HHe_IonE': A_HE, 'H_IonE': 0.0}


def log_m_atom(a_he):
    """log10 of the mass per hydrogen nucleus, which shifts log rho to log nH and log e_spec to
    log eint/nH"""
    return np.log10((1.0 + 4.0 * a_he) * M_P)


def progress(i, n, t0):
    """Report every eighth of the rows of a node loop"""
    if (i + 1) % max(1, n // 8) == 0:
        print(f'    row {i+1}/{n}  elapsed {time.perf_counter()-t0:5.1f}s', flush=True)


# ---------------------------------------------------------------- bicubic Hermite
def _basis(t):
    """Cubic Hermite basis h00, h10, h01, h11 at t in [0, 1]"""
    return np.array([2*t**3 - 3*t**2 + 1, t**3 - 2*t**2 + t, -2*t**3 + 3*t**2, t**3 - t**2])


def _basis_d1(t):
    """First derivatives of the basis in t"""
    return np.array([6*t**2 - 6*t, 3*t**2 - 4*t + 1, -6*t**2 + 6*t, 3*t**2 - 2*t])


def _basis_d2(t):
    """Second derivatives of the basis in t"""
    return np.array([12*t - 6, 6*t - 4, -12*t + 6, 6*t - 2])


def hermite_eval(x, y, table, x_axis, y_axis):
    """P, P_x, P_y, P_xy and P_yy of the bicubic Hermite polynomial of table at (x, y), per unit
    of the axes. It mirrors hermite_eval in mod_eos_lte_entropy.t, local coordinates clamped to
    the edge cell"""
    n1, n2 = table['value'].shape
    ix = max(0, min(n1 - 2, np.searchsorted(x_axis, x) - 1))
    iy = max(0, min(n2 - 2, np.searchsorted(y_axis, y) - 1))
    dx = x_axis[ix+1] - x_axis[ix]
    dy = y_axis[iy+1] - y_axis[iy]
    tx = max(0.0, min(1.0, (x - x_axis[ix]) / dx))
    ty = max(0.0, min(1.0, (y - y_axis[iy]) / dy))
    a, b = slice(ix, ix + 2), slice(iy, iy + 2)
    # Corner values and derivatives scaled to the unit cell, interleaved as the basis expects
    m = np.empty((4, 4))
    m[0::2, 0::2] = table['value'][a, b]
    m[0::2, 1::2] = table['y'][a, b]*dy
    m[1::2, 0::2] = table['x'][a, b]*dx
    m[1::2, 1::2] = table['xy'][a, b]*dx*dy
    hx, hy, dhx, dhy = _basis(tx), _basis(ty), _basis_d1(tx), _basis_d1(ty)
    return (hx @ m @ hy, dhx @ m @ hy / dx, hx @ m @ dhy / dy, dhx @ m @ dhy / (dx*dy),
            hx @ m @ _basis_d2(ty) / dy**2)


# ---------------------------------------------------------------- node placement
def fd_along_axis(field, coords, axis):
    """Derivative of a 2D field along axis 0 or 1 with non-uniform coordinates, centred inside
    and one-sided at the two ends"""
    f = np.moveaxis(field, axis, 0)
    out = np.empty_like(f)
    out[1:-1] = (f[2:] - f[:-2]) / (coords[2:] - coords[:-2])[:, None]
    out[0] = (f[1] - f[0]) / (coords[1] - coords[0])
    out[-1] = (f[-1] - f[-2]) / (coords[-1] - coords[-2])
    return np.moveaxis(out, 0, axis)


def gamma1_from_s_derivs(s_r, s_e, s_rr, s_re, s_ee):
    """Gamma_1 from the derivatives of the specific entropy s(log rho, log e_spec)

    Gamma_1 = 1 - s_r/s_e + (s_rr - s_re s_r/s_e)/(s_r ln10) - (s_re - s_ee s_r/s_e)/(s_e ln10)
    """
    rat = s_r / s_e
    return (1.0 - rat
            + (1.0 / (s_r * LN10)) * (s_rr - s_re * rat)
            - (1.0 / (s_e * LN10)) * (s_re - s_ee * rat))


def _normalised(d):
    return d / d.max() if d.max() > 0 else d


def _fill_nonfinite(a, fill):
    """Replace NaN and inf by the median of the finite entries, or by fill if there are none"""
    ok = np.isfinite(a)
    if ok.all():
        return a
    return np.where(ok, a, np.median(a[ok]) if ok.any() else fill)


def build_fine_s_grid(eos, n, lr_bounds, le_bounds):
    """Specific entropy on a uniform fine (log rho, log e_spec) grid, the reference from which
    the nodes are placed"""
    lr = np.linspace(lr_bounds[0], lr_bounds[1], n)
    le = np.linspace(le_bounds[0], le_bounds[1], n)
    s = np.zeros((n, n))
    t0 = time.perf_counter()
    print(f'  fine reference grid {n}x{n} (s) ...', flush=True)
    for i, lri in enumerate(lr):
        rho = 10.0 ** lri
        for j, lej in enumerate(le):
            s[i, j] = eos.T_from_e(rho, 10.0 ** lej * rho)[1]['s'] / rho
        progress(i, n, t0)
    return lr, le, s


def equidistribute_axes(lr_fine, le_fine, s_fine, n_out, baseline=0.05):
    """Adaptive (log rho, log e_spec) nodes, equidistributing a monitor on each axis

    Along log e_spec the monitor adds three normalised densities. The curvature of s keeps p and
    e accurate. The curvature of Gamma_1 resolves the ionisation features, and its fourth
    derivative to the power 1/4 is the error-optimal monitor of a fourth-order interpolant. The
    last alone gives the smallest Gamma_1 error but crowds the nodes into the He knee at 7-13 kK
    so much that the e -> p -> e round trip degrades to about 1 per cent

    Along log rho, s is nearly flat at fixed e_spec, so the curvature of s alone gives almost
    uniform spacing. The monitor adds the variation of ds/d log rho
    """
    # Gamma_1 on the fine grid, from the derivatives of s. The monitor needs only its shape, so
    # unphysical values at the edges are replaced by the median
    s_r = fd_along_axis(s_fine, lr_fine, axis=0)
    s_e = fd_along_axis(s_fine, le_fine, axis=1)
    s_rr = fd_along_axis(s_r, lr_fine, axis=0)
    s_re = fd_along_axis(s_r, le_fine, axis=1)
    s_ee = fd_along_axis(s_e, le_fine, axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        g1 = gamma1_from_s_derivs(s_r, s_e, s_rr, s_re, s_ee)
    physical = np.isfinite(g1) & (g1 > 0.5) & (g1 < 5.0)
    if not physical.all():
        g1 = np.where(physical, g1, np.median(g1[physical]) if physical.any() else 1.5)

    d4 = g1
    for _ in range(4):
        d4 = np.gradient(d4, axis=1)
    dens_d4 = np.power(_fill_nonfinite(np.abs(d4).sum(axis=0), 1.0) + 1e-30, 0.25)
    dens_le = (_normalised(monitor_density(s_fine, axis=1, baseline=0.0))
               + _normalised(monitor_density(g1, axis=1, baseline=0.0))
               + _normalised(dens_d4))
    le_axis = equidistribute(le_fine, dens_le + baseline * dens_le.max(), n_out)

    dens_lr = (_normalised(monitor_density(s_fine, axis=0, baseline=0.0))
               + _normalised(np.abs(np.gradient(s_r, axis=0)).sum(axis=1)))
    lr_axis = equidistribute(lr_fine, dens_lr + baseline * dens_lr.max(), n_out)
    return lr_axis, le_axis


def adapt_y_axis(lr_axis, y_bounds, n_out, reference, label, baseline=0.05):
    """Equidistribute the y axis of an inverse table on the curvature of log10(eint/nH) along y

    reference(lr_axis, y_fine) returns log10(eint/nH) on a regular y grid at twice the output
    density, from which the nodes are placed
    """
    n_ref = max(2 * n_out, 192)
    y_fine = np.linspace(y_bounds[0], y_bounds[1], n_ref)
    print(f'  inverse-{label} axis adapt: building reference at {len(lr_axis)}x{n_ref} ...',
          flush=True)
    dens = monitor_density(reference(lr_axis, y_fine), axis=1, baseline=baseline)
    y_axis = equidistribute(y_fine, dens, n_out)
    print(f'  inverse-{label} axis spacing after adapt: min {np.diff(y_axis).min():.4f}  '
          f'max {np.diff(y_axis).max():.4f}', flush=True)
    return y_axis


# ---------------------------------------------------------------- forward tables
def node_fwd(eos, lr, le):
    """T, p and ne/nH with their exact derivatives per log10 rho and log10 e_spec, and Gamma_1,
    at (lr, le). At fixed e_spec, dlnT/dlnrho = (E - E_r)/E_T, and at fixed rho,
    dlnT/dln e_spec = E/E_T"""
    rho = 10.0 ** lr
    T, th = eos.T_from_e(rho, 10.0 ** le * rho)
    dlnT_r = (th['e'] - th['E_r']) / th['E_T']
    dlnT_e = th['e'] / th['E_T']
    # Derivatives of each quantity per ln rho and ln T, combined below along the table axes
    dq = {'T': (0.0, T), 'p': (th['p']*th['chi_r'], th['p']*th['chi_T']),
          'neOnH': (th['ne_nH_r'], th['ne_nH_T'])}
    out = dict(T=T, p=th['p'], neOnH=th['ne_nH'], G1=th['G1'])
    for k, (q_r, q_T) in dq.items():
        out[k + '_x'] = LN10 * (q_r + q_T * dlnT_r)
        out[k + '_y'] = LN10 * q_T * dlnT_e
    return out


def mixed_fd(f, x, y):
    """Four-point mixed derivative of the exact value f at (x, y)"""
    h = STEP_MIXED
    return (f(x + h, y + h) - f(x + h, y - h) - f(x - h, y + h) + f(x - h, y - h)) / (4*h*h)


def build_forward_tables(eos, lr_axis, le_axis):
    """Forward tables on the adaptive (lr, le) grid, each a dict of the value, x, y and xy
    arrays. T, p and ne/nH have exact x and y, and xy from the central difference of the exact x
    along y. Gamma_1 has x and y from central differences and xy from the four-point stencil"""
    n1, n2 = len(lr_axis), len(le_axis)
    names = ('T', 'p', 'neOnH', 'G1')
    tab = {k: {d: np.zeros((n1, n2)) for d in ('value', 'x', 'y', 'xy')} for k in names}
    g1 = lambda lr, le: node_fwd(eos, lr, le)['G1']
    t0 = time.perf_counter()
    print(f'  forward nodes {n1}x{n2}, exact derivatives ...', flush=True)
    for i, lri in enumerate(lr_axis):
        for j, lej in enumerate(le_axis):
            c = node_fwd(eos, lri, lej)
            up = node_fwd(eos, lri, lej + STEP_FD)
            dn = node_fwd(eos, lri, lej - STEP_FD)
            for k in ('T', 'p', 'neOnH'):
                tab[k]['value'][i, j] = c[k]
                tab[k]['x'][i, j] = c[k + '_x']
                tab[k]['y'][i, j] = c[k + '_y']
                tab[k]['xy'][i, j] = (up[k + '_x'] - dn[k + '_x']) / (2*STEP_FD)
            tab['G1']['value'][i, j] = c['G1']
            tab['G1']['x'][i, j] = (g1(lri + STEP_FD, lej) - g1(lri - STEP_FD, lej)) / (2*STEP_FD)
            tab['G1']['y'][i, j] = (up['G1'] - dn['G1']) / (2*STEP_FD)
            tab['G1']['xy'][i, j] = mixed_fd(g1, lri, lej)
        progress(i, n1, t0)
    return {'g1e': tab['G1'], 'neOnH': tab['neOnH'], 'Tfwd': tab['T'], 'pfwd': tab['p']}


# ---------------------------------------------------------------- inverse tables
def saha_log_e_at_p(eos, lr_axis, lp_axis, log_m):
    """log10(eint/nH) of the Saha state at (log rho, log p/nH), the reference for the lp nodes"""
    log_e = np.zeros((len(lr_axis), len(lp_axis)))
    print(f'  inverse-p Saha pre-images {len(lr_axis)}x{len(lp_axis)} ...', flush=True)
    for i, lri in enumerate(lr_axis):
        rho = 10.0 ** lri
        n_h = rho / (10.0 ** log_m)
        for j, lpj in enumerate(lp_axis):
            log_e[i, j] = np.log10(eos.T_from_p(rho, 10.0 ** lpj * n_h)[1]['e'] / n_h)
    return log_e


def saha_log_e_at_T(eos, lr_axis, lT_axis, log_m):
    """log10(eint/nH) of the Saha state at (log rho, log T), the reference for the lT nodes"""
    log_e = np.zeros((len(lr_axis), len(lT_axis)))
    for i, lri in enumerate(lr_axis):
        rho = 10.0 ** lri
        n_h = rho / (10.0 ** log_m)
        for j, lTj in enumerate(lT_axis):
            log_e[i, j] = np.log10(eos.thermo(rho, 10.0 ** lTj)['e'] / n_h)
    return log_e


def pfwd_preimages(lr_axis, lp_axis, pfwd, nh_axis, e_nh_axis, log_m):
    """log10(eint/nH) on (lr, lp) as the pre-image under the pfwd polynomial, by bisection in
    log10(eint/nH), so that the runtime round trip eintP -> pfwd is the identity at every node"""
    n1, n2 = len(lr_axis), len(lp_axis)
    log_e = np.zeros((n1, n2))
    nonmono = int((np.diff(pfwd['value'], axis=1) <= 0).sum())
    if nonmono > 0:
        print(f'  warning: pfwd has {nonmono} non-monotonic cells along axis 2, so the bisection '
              f'may clamp or oscillate')
    t0 = time.perf_counter()
    print(f'  inverse-p polynomial bisection {n1}x{n2} ...', flush=True)
    for i, lri in enumerate(lr_axis):
        log_n_h = float(lri - log_m)
        for j, lpj in enumerate(lp_axis):
            p_target = 10.0 ** (lpj + log_n_h)
            lo, hi = float(e_nh_axis[0]), float(e_nh_axis[-1])
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if hermite_eval(log_n_h, mid, pfwd, nh_axis, e_nh_axis)[0] < p_target:
                    lo = mid
                else:
                    hi = mid
            log_e[i, j] = 0.5 * (lo + hi)
        progress(i, n1, t0)
    edge = int((np.abs(log_e - e_nh_axis[0]) < 1e-12).sum()
               + (np.abs(log_e - e_nh_axis[-1]) < 1e-12).sum())
    if edge:
        raise ValueError(f'{edge} eintP pre-images pinned to the pfwd e-axis edge')
    return log_e


def eintP_derivs(log_e, lr_axis, pfwd, nh_axis, e_nh_axis, log_m):
    """Derivatives of the eintP pre-images Y(x, lp) under the pfwd polynomial P(x, Y)

    With x = log10 nH and Y = log10(e/nH), the inverse-function rule on log10 P - x = lp gives
    Y_x = (P ln10 - P_x)/P_Y, Y_lp = P ln10/P_Y and
    Y_xlp = ln10 P (ln10 P_Y - P_xY - Y_x P_YY)/P_Y^2
    """
    y_x, y_lp, y_xlp = (np.zeros_like(log_e) for _ in range(3))
    knot = 0
    for i, lri in enumerate(lr_axis):
        x = float(lri - log_m)
        for j in range(log_e.shape[1]):
            y = log_e[i, j]
            knot += int(np.any(np.abs(e_nh_axis - y) < 1e-12))
            p, p_x, p_y, p_xy, p_yy = hermite_eval(x, y, pfwd, nh_axis, e_nh_axis)
            y_x[i, j] = (p*LN10 - p_x) / p_y
            y_lp[i, j] = p*LN10 / p_y
            y_xlp[i, j] = LN10 * p * (LN10*p_y - p_xy - y_x[i, j]*p_yy) / p_y**2
    print(f'  eintP inverse-function derivatives; pre-images on a pfwd y-knot: {knot}', flush=True)
    return dict(value=log_e, x=y_x, y=y_lp, xy=y_xlp)


def build_g1p(eos, lr_axis, lp_axis, log_m):
    """Gamma_1 on the inverse-p grid, exact at the Saha state of (rho, p). x and y come from
    central differences and xy from the four-point stencil of exact values"""
    m_h = 10.0 ** log_m

    def g1(lr, lp):
        rho = 10.0 ** lr
        return eos.T_from_p(rho, 10.0 ** lp * rho / m_h)[1]['G1']

    n1, n2 = len(lr_axis), len(lp_axis)
    v, vx, vy, vxy = (np.zeros((n1, n2)) for _ in range(4))
    for i, lri in enumerate(lr_axis):
        for j, lpj in enumerate(lp_axis):
            v[i, j] = g1(lri, lpj)
            vx[i, j] = (g1(lri + STEP_FD, lpj) - g1(lri - STEP_FD, lpj)) / (2*STEP_FD)
            vy[i, j] = (g1(lri, lpj + STEP_FD) - g1(lri, lpj - STEP_FD)) / (2*STEP_FD)
            vxy[i, j] = mixed_fd(g1, lri, lpj)
    return dict(value=v, x=vx, y=vy, xy=vxy)


def node_invT(eos, lr, lT):
    """log10(eint/nH) at (rho, T) and its exact derivatives per log10 rho and log10 T"""
    th = eos.thermo(10.0 ** lr, 10.0 ** lT)
    return np.log10(th['e'] / th['nH']), th['E_r'] / th['e'] - 1.0, th['E_T'] / th['e']


def build_inverse_T_tables(eos, lr_axis, lT_axis):
    """log10(eint/nH) on the (lr, lT) grid, with exact x and y, and xy from the central
    difference of the exact x along y"""
    n1, n2 = len(lr_axis), len(lT_axis)
    v, vx, vy, vxy = (np.zeros((n1, n2)) for _ in range(4))
    for i, lri in enumerate(lr_axis):
        for j, lTj in enumerate(lT_axis):
            v[i, j], vx[i, j], vy[i, j] = node_invT(eos, lri, lTj)
            vxy[i, j] = (node_invT(eos, lri, lTj + STEP_FD)[1]
                         - node_invT(eos, lri, lTj - STEP_FD)[1]) / (2*STEP_FD)
    return dict(value=v, x=vx, y=vy, xy=vxy)


# ---------------------------------------------------------------- verification
def verify_tables(eos, log_m, lr_fwd, le_fwd, fwd, lr_invp, lp_invp, eintP, g1p, log_e_T,
                  rng_seed=2026):
    """Checks 1 to 5 of the module docstring. Raises on the first failure, else returns a report"""
    rng = np.random.default_rng(rng_seed)
    report = {}
    T = fwd['Tfwd']['value']
    p = fwd['pfwd']['value']

    # 1. Maxwell relation d(1/T)/dx|e = d(-p/(rho T))/dy|rho, which is
    #    e_spec T_x = (T p_y - p T_y)/rho. It holds only if e, p and their exact derivatives
    #    come from one free energy
    rho = 10.0 ** lr_fwd[:, None]
    e_spec = 10.0 ** le_fwd[None, :]
    lhs = e_spec * fwd['Tfwd']['x']
    rhs = (T * fwd['pfwd']['y'] - p * fwd['Tfwd']['y']) / rho
    err1 = float((np.abs(lhs - rhs) / np.abs(T * fwd['pfwd']['y'] / rho)).max())
    report['maxwell_max_rel_err'] = err1
    if not err1 <= 1e-9:
        raise ValueError(f'Maxwell relation violated: max rel err {err1:.3e}')
    print(f'  [pass] Maxwell relation: max rel err {err1:.3e}', flush=True)

    # 2. e -> p -> e round trip at random interior nodes, through the runtime bicubic Hermite.
    #    The outer 4 nodes per axis are skipped, since the operating range stays well inside
    n_check = 1000
    i_idx = rng.integers(4, len(lr_fwd)-4, n_check)
    j_idx = rng.integers(4, len(le_fwd)-4, n_check)
    errs, kept, n_outside = [], [], 0
    for k in range(n_check):
        lri = lr_fwd[i_idx[k]]
        lej = le_fwd[j_idx[k]]
        rho_k = 10.0 ** lri
        n_h = rho_k / (10.0 ** log_m)
        lp_target = np.log10(eos.thermo(rho_k, T[i_idx[k], j_idx[k]])['p'] / n_h)
        # The forward table reaches higher energies than the inverse-p axis covers. Such states
        # are never inverted at runtime, which switches to FI above 200 kK
        if lp_target < lp_invp[0] or lp_target > lp_invp[-1]:
            n_outside += 1
            continue
        kept.append(k)
        log_e = hermite_eval(lri, lp_target, eintP, lr_invp, lp_invp)[0]
        errs.append(abs(10.0 ** log_e * n_h / (10.0 ** lej * rho_k) - 1.0))
    err_arr = np.array(errs)
    err2 = float(err_arr.max())
    report['roundtrip_max_rel_err'] = err2
    report['roundtrip_p99_rel_err'] = float(np.percentile(err_arr, 99.0))
    report['roundtrip_n_outside_invp'] = n_outside
    # Bicubic Hermite is fourth order: 1e-3 at production N >= 128, looser at smoke resolutions
    n = len(lr_fwd)
    thresh = 1e-3 if n >= 128 else (5e-2 if n >= 64 else 5e-1)
    if err2 > thresh:
        kbad = kept[int(np.argmax(err_arr))]
        print(f'  worst: lr {lr_fwd[i_idx[kbad]]:.3f} le {le_fwd[j_idx[kbad]]:.3f} '
              f'T {T[i_idx[kbad], j_idx[kbad]]:.3e}', flush=True)
        raise ValueError(f'e -> p -> e round-trip violated: max rel err {err2:.3e}')
    print(f'  [pass] round-trip e-p-e: max {err2:.3e}  ({n_outside} of {n_check} nodes outside '
          f'the inverse-p axis, excluded)', flush=True)

    # 3. Gamma_1 tables finite and physical
    for name, g in (('g1e', fwd['g1e']['value']), ('g1p', g1p['value'])):
        if not np.isfinite(g).all() or g.min() < 0.8 or g.max() > 2.5:
            raise ValueError(f'{name} non-finite or outside [0.8, 2.5]: '
                             f'[{g.min():.3f}, {g.max():.3f}]')
        report[name + '_min'] = float(g.min())
        report[name + '_max'] = float(g.max())
        print(f'  [pass] {name} bounds: [{g.min():.3f}, {g.max():.3f}]', flush=True)

    # 4. Inverse tables rise along y, since more pressure or temperature means more eint
    for name, v in (('eintP', eintP['value']), ('eintT', log_e_T)):
        bad = int((np.diff(v, axis=1) <= 0).sum())
        if bad > 0.05 * v[:, 1:].size:
            raise ValueError(f'{name} non-monotonic in y at {bad}/{v[:, 1:].size} cells')
        report[name + '_nonmono'] = bad
    print('  [pass] inverse tables monotonic in y', flush=True)

    # 5. Electrons per H between none and full ionisation of H and He
    ne = fwd['neOnH']['value']
    if ne.min() < -1e-12 or ne.max() > 1.0 + 2.0 * eos.a_he + 1e-6:
        raise ValueError(f'neOnH out of [0, {1.0 + 2.0*eos.a_he}]: min {ne.min()} max {ne.max()}')
    report['neOnH_max'] = float(ne.max())
    print(f'  [pass] neOnH bounds: [0, {ne.max():.3f}]', flush=True)
    return report


# ---------------------------------------------------------------- output
def write_bin_table(filepath, data, var1_bounds, var2_bounds, var1_nodes, var2_nodes):
    """One table in the binary format the AMRVAC loader reads

    Header   [2 i4] dim1, dim2, then [4 f8] var1_min, var1_max, var2_min, var2_max
    Body     [dim1*dim2 f8] data, Fortran column-major
    Trailer  [1 i4] flag = 1, then [dim1 f8] var1 nodes and [dim2 f8] var2 nodes
    """
    if not np.isfinite(data).all():
        raise ValueError(f"non-finite data in {filepath}")
    assert var1_nodes.shape == (data.shape[0],)
    assert var2_nodes.shape == (data.shape[1],)
    with open(filepath, 'wb') as f:
        f.write(struct.pack('2i', data.shape[0], data.shape[1]))
        f.write(struct.pack('4d', *var1_bounds, *var2_bounds))
        data.T.tofile(f)
        f.write(struct.pack('i', 1))
        var1_nodes.astype(np.float64).tofile(f)
        var2_nodes.astype(np.float64).tofile(f)


def write_quantity(outdir, stem, comp, table, nodes1, nodes2):
    """Write the value and the three derivative files of one tabulated quantity"""
    bounds1 = (float(nodes1[0]), float(nodes1[-1]))
    bounds2 = (float(nodes2[0]), float(nodes2[-1]))
    for suffix, key in zip(('', '_x', '_y', '_xy'), ('value', 'x', 'y', 'xy')):
        write_bin_table(f'{outdir}/LTEeos_{stem}{suffix}_{comp}.bin', table[key], bounds1, bounds2,
                        nodes1, nodes2)


# ---------------------------------------------------------------- main
def main():
    parser = argparse.ArgumentParser(description='Build the entropy-method EoS tables')
    parser.add_argument('--n', type=int, default=128, help='nodes per table axis')
    parser.add_argument('--ref', type=int, default=512, help='points per axis of the fine grid')
    parser.add_argument('--comp', default='HHe_IonE', choices=sorted(COMPOSITIONS))
    parser.add_argument('--outdir', default='.', help='where to write the 28 .bin files')
    parser.add_argument('--baseline', type=float, default=0.05,
                        help='uniform share of every node-placement monitor')
    args = parser.parse_args()

    eos = SahaHHe(COMPOSITIONS[args.comp])
    log_m = log_m_atom(eos.a_he)
    lr_bounds = (LOG_NH_BOUNDS[0] + log_m, LOG_NH_BOUNDS[1] + log_m)
    le_bounds = (LOG_E_NH_BOUNDS[0] - log_m, LOG_E_NH_BOUNDS[1] - log_m)

    print(f'=== entropy-method full table build ({args.comp}, n={args.n}) ===')
    print(f'  axis-1 (log nH) :   [{LOG_NH_BOUNDS[0]:.2f}, {LOG_NH_BOUNDS[1]:.2f}]  CGS')
    print(f'  axis-2 (log e/nH):  [{LOG_E_NH_BOUNDS[0]:.2f}, {LOG_E_NH_BOUNDS[1]:.2f}]  CGS')

    print('\nStep 1/6: fine reference grid ...')
    lr_fine, le_fine, s_fine = build_fine_s_grid(eos, args.ref, lr_bounds, le_bounds)

    print('\nStep 2/6: equidistribute adaptive (lr, le) ...')
    lr_axis, le_axis = equidistribute_axes(lr_fine, le_fine, s_fine, args.n, baseline=args.baseline)
    print(f'  lr spacing: min {np.diff(lr_axis).min():.4f}  max {np.diff(lr_axis).max():.4f}')
    print(f'  le spacing: min {np.diff(le_axis).min():.4f}  max {np.diff(le_axis).max():.4f}')

    # The tables are built per log rho and log e_spec. The loader reads them per log nH and
    # log eint/nH, which are constant shifts of those axes
    print('\nStep 3/6: forward tables on adaptive grid ...')
    fwd = build_forward_tables(eos, lr_axis, le_axis)
    nh_axis = lr_axis - log_m
    e_nh_axis = le_axis + log_m
    print(f'  AMRVAC axis-1 span: [{nh_axis[0]:.3f}, {nh_axis[-1]:.3f}]')
    print(f'  AMRVAC axis-2 span: [{e_nh_axis[0]:.3f}, {e_nh_axis[-1]:.3f}]')

    # The inverse tables share the forward log rho axis
    print('\nStep 4/6: inverse (lr, lp) tables ...')
    lp_axis = adapt_y_axis(lr_axis, LOG_P_NH_BOUNDS, args.n,
                           lambda lr, lp: saha_log_e_at_p(eos, lr, lp, log_m), 'p',
                           baseline=args.baseline)
    log_e_p = pfwd_preimages(lr_axis, lp_axis, fwd['pfwd'], nh_axis, e_nh_axis, log_m)
    eintP = eintP_derivs(log_e_p, lr_axis, fwd['pfwd'], nh_axis, e_nh_axis, log_m)
    g1p = build_g1p(eos, lr_axis, lp_axis, log_m)

    print('\nStep 5/6: inverse (lr, lT) tables ...')
    lT_axis = adapt_y_axis(lr_axis, LOG_T_BOUNDS, args.n,
                           lambda lr, lT: saha_log_e_at_T(eos, lr, lT, log_m), 'T',
                           baseline=args.baseline)
    eintT = build_inverse_T_tables(eos, lr_axis, lT_axis)

    print('\nStep 6/6: verification ...')
    report = verify_tables(eos, log_m, lr_axis, le_axis, fwd, lr_axis, lp_axis, eintP, g1p,
                           eintT['value'])
    print('  All checks passed.')

    print('\nWriting files ...')
    os.makedirs(args.outdir, exist_ok=True)
    for stem in ('g1e', 'neOnH', 'Tfwd', 'pfwd'):
        write_quantity(args.outdir, stem, args.comp, fwd[stem], nh_axis, e_nh_axis)
    write_quantity(args.outdir, 'eintP', args.comp, eintP, nh_axis, lp_axis)
    write_quantity(args.outdir, 'g1p', args.comp, g1p, nh_axis, lp_axis)
    write_quantity(args.outdir, 'eintT', args.comp, eintT, nh_axis, lT_axis)

    with open(f'{args.outdir}/build_report.json', 'w') as f:
        json.dump(dict(
            comp=args.comp, n=args.n, ref=args.ref,
            forward_bounds_log_nH=[float(nh_axis[0]), float(nh_axis[-1])],
            forward_bounds_log_e_nh=[float(e_nh_axis[0]), float(e_nh_axis[-1])],
            inverse_p_bounds_log_p_nh=[float(lp_axis[0]), float(lp_axis[-1])],
            inverse_T_bounds_log_T=[float(lT_axis[0]), float(lT_axis[-1])],
            verification=report,
        ), f, indent=2)
    print(f'\n  Wrote 28 .bin files + build_report.json under {args.outdir}/')
    print('  Build complete.')


if __name__ == '__main__':
    main()
