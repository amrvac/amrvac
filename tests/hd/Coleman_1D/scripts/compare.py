"""Coleman_1D comparison harness: AMRVAC dat -> cached Riemann + Snow check.

Two independent references, both invoked at autotest time:

1. **Exact Riemann reference** (cached): the precomputed (x, rho, v, eint, p)
   profile at t_final for each (test, variant). Loaded from
   ../reference_riemann/<name>.npz. Regenerable via generate_references.py
   when ICs or table set changes, but the autotest itself does not run the
   Riemann solver.

2. **Snow two-gamma jump check** (live): the autotest computes the
   compression ratio across each shock side from the analytical two-gamma
   Rankine-Hugoniot quadratic and compares against both the exact reference
   and the AMRVAC sim. This is fast (~10 ms per test) so it runs inline.

The combination catches different failure modes:
  - exact Riemann disagreement -> AMRVAC cons<->prim or time-stepping bug
  - Snow disagreement only     -> AMRVAC's runtime Gamma_1 query is off
  - both agreeing -> code is right
"""

from __future__ import annotations

import os
import struct
from pathlib import Path
import numpy as np

from snow_jump import SnowJumpSolver
from eos_reader import make_eos, MP_CGS, KB_CGS

SCRIPT_DIR = Path(__file__).resolve().parent
REF_DIR    = SCRIPT_DIR.parent / 'reference_riemann'


# -- AMRVAC dat reader (1D HD with extra g1 variable) -------------------

def read_dat_1d(path: str | Path) -> dict:
    """Minimal 1D AMRVAC v5 dat reader for HD output with saveprim=.true.

    Returns dict with keys present per the file's `var_names` plus
    {x, time, domain_nx, nw, var_names}. The Coleman parfiles save in
    primitive form, so callers should access ['rho'] for density and
    ['m1'] for the x-momentum (despite being saved as velocity -- the
    saveprim convention exposes velocity under 'm1').
    """
    path = Path(path)
    with open(path, 'rb') as f:
        def ri4(): return struct.unpack('<i', f.read(4))[0]
        def ri8(): return struct.unpack('<q', f.read(8))[0]
        def rd():  return struct.unpack('<d', f.read(8))[0]

        version       = ri4()
        offset_tree   = ri4()
        offset_blocks = ri4()
        nw     = ri4()
        ndir   = ri4()
        ndim   = ri4()
        levmax = ri4()
        nleafs   = ri4()
        nparents = ri4()
        it       = ri4()
        time     = rd()
        xmin = [rd() for _ in range(ndim)]
        xmax = [rd() for _ in range(ndim)]
        domain_nx = [ri4() for _ in range(ndim)]
        block_nx  = [ri4() for _ in range(ndim)]
        periodic  = [ri4() != 0 for _ in range(ndim)]
        geometry  = f.read(16).decode('ascii').strip().rstrip('\x00')
        stagger   = ri4() != 0
        var_names = [f.read(16).decode('ascii').strip().rstrip('\x00') for _ in range(nw)]
        phys_type = f.read(16).decode('ascii').strip().rstrip('\x00')
        n_par = ri4()
        par_values = [rd() for _ in range(n_par)]
        par_names  = [f.read(16).decode('ascii').strip().rstrip('\x00') for _ in range(n_par)]
        snapshotnext = ri4(); slicenext = ri4(); collapsenext = ri4()

        f.seek(offset_tree)
        ntotal = nleafs + nparents
        leaf_flags = [struct.unpack('<i', f.read(4))[0] != 0 for _ in range(ntotal)]
        block_lvl  = np.array([ri4() for _ in range(nleafs)])
        block_ig   = np.array([ri4() for _ in range(nleafs * ndim)]).reshape(nleafs, ndim)
        block_off  = np.array([ri8() for _ in range(nleafs)])

        nx = block_nx[0]
        blocks = []
        for ib in range(nleafs):
            f.seek(block_off[ib])
            n_ghost = [ri4() for _ in range(2 * ndim)]
            ng_lo, ng_hi = n_ghost[0], n_ghost[1]
            nx_g = nx + ng_lo + ng_hi
            w_blk = np.frombuffer(f.read(nw * nx_g * 8), dtype='<f8').reshape(nw, nx_g)
            w_blk = w_blk[:, ng_lo:nx_g - ng_hi] if ng_hi > 0 else w_blk[:, ng_lo:]
            lvl = block_lvl[ib]; ig = block_ig[ib, 0]
            dx = (xmax[0] - xmin[0]) / (domain_nx[0] * 2 ** (lvl - 1))
            x0 = xmin[0] + (ig - 1) * nx * dx
            x_block = x0 + (np.arange(nx) + 0.5) * dx
            blocks.append((x_block, w_blk))

        blocks.sort(key=lambda b: b[0][0])
        x_all = np.concatenate([b[0] for b in blocks])
        w_all = np.concatenate([b[1] for b in blocks], axis=1)
        idx = np.argsort(x_all); x_all = x_all[idx]; w_all = w_all[:, idx]

        result = {'x': x_all, 'time': time, 'var_names': var_names, 'nw': nw,
                  'domain_nx': domain_nx[0]}
        for iv, name in enumerate(var_names):
            result[name] = w_all[iv]
        return result


# -- Reference comparison ----------------------------------------------

def compare_against_riemann(dat: dict, reference_path: Path,
                             ic: dict | None = None,
                             mask_eps: float = 2e-2) -> dict:
    """Compare AMRVAC final-state arrays against a cached exact Riemann
    reference. Returns a metrics dict; masks domain edges and near-vacuum
    cells (<=10% of max rho) so the percentile metrics reflect bulk-flow
    accuracy rather than shock-smear or floor-clipped cells.

    mask_eps default 2e-2 (20-cell buffer at 1024 cells). The collision
    tests (Coleman 3, 4) at t_final = 3.776 have both outgoing shocks
    leaving the domain well before the snapshot; the continuous BCs then
    inject the initial L/R states back in, polluting the outer ~15-20
    cells. The Riemann reference assumes the wave structure extends to
    infinity, so it disagrees with the sim in that polluted zone. mask_eps
    = 2e-2 excludes this region; 5e-3 (the value commonly used for Sod
    setups where the waves never reach the boundary) is too tight.

    The velocity error is normalised by the LARGER of max|v_ref| and
    max(|v_L|,|v_R|) from the ICs (when supplied). For colliding-flow
    problems (tests 3 and 4) the post-collision state has v ~= 0 across
    the domain -- without the IC-scale fallback, sub-percent absolute
    deviations would inflate to 70 % relative error."""
    ref = np.load(reference_path)
    x_ref   = ref['x']
    rho_ref = ref['rho']
    v_ref   = ref['v']
    p_ref   = ref['p']
    eint_ref = ref['eint']

    # The sim is on a uniform grid; resample the reference onto the sim grid
    # if their resolutions differ.
    x_sim = dat['x']
    if x_ref.shape != x_sim.shape or not np.allclose(x_ref, x_sim):
        rho_ref = np.interp(x_sim, x_ref, rho_ref)
        v_ref   = np.interp(x_sim, x_ref, v_ref)
        p_ref   = np.interp(x_sim, x_ref, p_ref)
        eint_ref = np.interp(x_sim, x_ref, eint_ref)

    rho_sim = dat['rho']
    v_sim   = dat['m1'] / np.maximum(rho_sim, 1e-30)
    eint_sim = dat['e'] - 0.5 * rho_sim * v_sim ** 2

    L = x_sim[-1] - x_sim[0]
    interior = (x_sim > x_sim[0] + mask_eps * L) & (x_sim < x_sim[-1] - mask_eps * L)
    not_vacuum = np.abs(rho_ref) > 0.1 * np.abs(rho_ref).max()
    mask = interior & not_vacuum

    def rel(a, b):
        return np.abs(a - b) / np.maximum(np.abs(b), 1e-12)

    rrho = rel(rho_sim, rho_ref)[mask]
    ic_v_scale = (max(abs(ic.get('v_L', 0.0)), abs(ic.get('v_R', 0.0)))
                  if ic else 0.0)
    v_scale = max(np.abs(v_ref[mask]).max(), ic_v_scale, 1e-12)
    abs_v = (np.abs(v_sim - v_ref) / v_scale)[mask]
    reint = rel(eint_sim, eint_ref)[mask]

    return {
        'time':         float(dat['time']),
        'n_cells':      int(dat['domain_nx']),
        'n_eval':       int(mask.sum()),
        # p95 is the headline metric; it captures bulk-flow accuracy
        # while excluding the 5 % of cells dominated by FV shock and
        # contact-discontinuity smearing. p99 keeps for diagnostics --
        # on a 1024-cell run it sits inside the 2-3 cells of shock smear
        # and inflates whenever the wave structure is sharp (LTE pureH
        # collision compressions of 10x routinely give p99 = 25 % from
        # just those cells, even when the bulk wave structure is right
        # to under 1 %).
        'p95_rel_rho':  float(np.percentile(rrho, 95)),
        'p99_rel_rho':  float(np.percentile(rrho, 99)),
        'max_rel_rho':  float(rrho.max()),
        'mean_rel_rho': float(rrho.mean()),
        'p95_abs_v_norm': float(np.percentile(abs_v, 95)),
        'p99_abs_v_norm': float(np.percentile(abs_v, 99)),
        'max_abs_v_norm': float(abs_v.max()),
        'p95_rel_eint': float(np.percentile(reint, 95)),
        'p99_rel_eint': float(np.percentile(reint, 99)),
        'max_rel_eint': float(reint.max()),
    }


def snow_check(dat: dict, ic: dict, eos, shock_sides: tuple[str, ...],
               reference_path: Path | None = None) -> dict:
    """Snow two-gamma jump-condition check for a single test variant.

    For each shock side, computes the compression ratio four ways:
      - r_exact_ic    : standard RH with gamma = gamma(upstream) (single-gamma baseline)
      - r_snow_g1     : two-gamma quadratic with Gamma_1(pre) and Gamma_1(post)
      - r_snow_geff   : two-gamma quadratic with gamma_eff = 1 + p/eint
                        (captures ionisation latent heat; this is the
                        flavour that overlaps with the LTE sim post-shock
                        density, whereas Gamma_1 underpredicts it)
      - r_amrvac      : measured from the AMRVAC density profile

    Pre-shock state is the IC (rho_L, v_L, T_L) or (rho_R, v_R, T_R);
    pressure and Gamma_1 there come from direct EoS lookups so that for
    LTE p_pre != (gamma_1 - 1) eint_pre (else geff_pre == g1_pre, a
    trivial identity that defeats the whole comparison).

    Post-shock state is sampled from the cached exact Riemann reference
    (reference_path) at a cell just downstream of the shock. Without
    reference_path the post-shock state is estimated by a single-gamma
    adiabatic RH step, which is a poor approximation for LTE shocks
    through an ionisation zone -- the geff result will then be very
    close to the g1 result rather than capturing the real difference.

    For FI variants this reduces to standard RH (g1 = geff = 5/3), so
    the two markers coincide -- expected.
    """
    snow = SnowJumpSolver()
    F = eos.nH2rhoFactor
    rho_sim = dat['rho']
    v_sim   = dat['m1'] / np.maximum(rho_sim, 1e-30)
    x = dat['x']
    x_mid = 0.5 * (x[0] + x[-1])

    ref = np.load(reference_path) if reference_path is not None else None
    if ref is not None:
        rho_ref  = np.interp(x, ref['x'], ref['rho'])
        eint_ref = np.interp(x, ref['x'], ref['eint'])

    out: dict = {}
    for side in shock_sides:
        if side == 'L':
            rho_pre, v_pre, T_pre = ic['rho_L'], ic['v_L'], ic['T_L']
        else:
            rho_pre, v_pre, T_pre = ic['rho_R'], ic['v_R'], ic['T_R']

        nH_pre = rho_pre / F
        log_nH_pre   = float(np.log10(nH_pre))
        eint_nH_pre  = eos.eint_nH_from_T(log_nH_pre, T_pre)
        log_eint_pre = float(np.log10(eint_nH_pre))
        eint_pre = nH_pre * eint_nH_pre

        # Pre-shock pressure: direct EoS lookup (NOT (g1-1)*eint -- that
        # makes geff trivially equal to g1).
        p_pre  = float(eos.pressure(log_nH_pre, log_eint_pre))
        g1_pre = float(eos.gamma1(log_nH_pre, log_eint_pre))
        geff_pre = 1.0 + p_pre / eint_pre

        # Sound speeds in g1 and geff flavours
        cs_g1_pre   = float(np.sqrt(g1_pre   * p_pre / rho_pre))
        cs_geff_pre = float(np.sqrt(geff_pre * p_pre / rho_pre))

        # Locate the shock as the largest |grad v| on the appropriate
        # half; sample upstream/downstream 10 cells away. The |grad rho|
        # heuristic is unreliable because in the FI Sod-like setups the
        # CONTACT discontinuity has a larger density jump than the shock
        # (e.g. coleman 2 FI: contact rho ratio 4.7, shock rho ratio 3.3
        # for gamma=5/3). Across the contact v is continuous and rho
        # jumps; across the shock both jump. |grad v| is exactly zero at
        # the contact and peaks at the shock for any EoS.
        half = x < x_mid if side == 'L' else x > x_mid
        if not half.any():
            continue
        idx_in_half = int(np.argmax(np.abs(np.gradient(v_sim[half]))))
        global_idx = int(np.where(half)[0][idx_in_half])
        n = len(rho_sim)
        offset = -10 if side == 'L' else 10
        i_pre  = max(0, min(n - 1, global_idx + offset))
        i_post = max(0, min(n - 1, global_idx - offset // 2))
        rho_post_sim = float(rho_sim[i_post])
        rho_pre_sim  = float(rho_sim[i_pre])
        r_amrvac = rho_post_sim / max(rho_pre_sim, 1e-30)

        # Mach number in shock frame from the sim shock speed.
        v_pre_sim  = float(v_sim[i_pre])
        v_post_sim = float(v_sim[i_post])
        drho = rho_post_sim - rho_pre_sim
        if abs(drho) < 1e-12:
            continue
        s_shock = (rho_post_sim * v_post_sim - rho_pre_sim * v_pre_sim) / drho
        v_shock_frame = abs(v_pre - s_shock)
        M1_g1   = v_shock_frame / max(cs_g1_pre,   1e-12)
        M1_geff = v_shock_frame / max(cs_geff_pre, 1e-12)

        # Post-shock state: prefer the cached exact reference (correctly
        # accounts for the LTE star state); fall back to a single-gamma
        # adiabatic estimate if no reference was passed in.
        if ref is not None:
            shock_x = float(x[global_idx])
            j_post = (np.where(x < shock_x)[0][-1] if side == 'R'
                       else np.where(x > shock_x)[0][0])
            rho_post = float(rho_ref[j_post])
            eint_post = float(eint_ref[j_post])
        else:
            r_std = snow.standard_rh(M1_g1, g1_pre)
            rho_post = rho_pre * r_std
            eint_post = (p_pre * (2 * g1_pre * M1_g1 ** 2 - (g1_pre - 1))
                          / (g1_pre + 1) / (g1_pre - 1))
        nH_post = rho_post / F
        log_nH_post  = float(np.log10(nH_post))
        log_eint_post = float(np.log10(max(eint_post / nH_post, 1e-30)))
        p_post   = float(eos.pressure(log_nH_post, log_eint_post))
        g1_post  = float(eos.gamma1(log_nH_post, log_eint_post))
        geff_post = 1.0 + p_post / eint_post

        r_std       = snow.standard_rh(M1_g1, g1_pre)
        r_snow_g1   = float(snow.solve(M1_g1,   g1_pre,   g1_post))
        r_snow_geff = float(snow.solve(M1_geff, geff_pre, geff_post))

        out[side] = {
            'M1_g1':   float(M1_g1),
            'M1_geff': float(M1_geff),
            'rho_pre': float(rho_pre),
            'r_exact_ic':  float(r_std),
            'r_snow_g1':   r_snow_g1,
            'r_snow_geff': r_snow_geff,
            'r_amrvac':    r_amrvac,
            'g1_pre':   g1_pre,   'g1_post':   g1_post,
            'geff_pre': geff_pre, 'geff_post': geff_post,
            'err_amrvac_vs_std':       abs(r_amrvac - r_std)       / r_std,
            'err_amrvac_vs_snow_g1':   abs(r_amrvac - r_snow_g1)   / r_snow_g1,
            'err_amrvac_vs_snow_geff': abs(r_amrvac - r_snow_geff) / r_snow_geff,
            # Back-compat alias so existing runner code keeps working;
            # gates on the geff flavour (the one that matches LTE physics).
            'err_amrvac_vs_snow': abs(r_amrvac - r_snow_geff) / r_snow_geff,
        }
    return out


def gate(cmp: dict, max_rho: float = 0.08, max_v: float = 0.02,
         max_eint: float = 0.05) -> tuple[bool, list[str]]:
    """Pass/fail on the p95 norm of |sim-ref|/|ref| over interior cells.

    Why p95, not p99: the FV scheme inevitably smears shocks and contact
    discontinuities across 2-3 cells, producing 10-30 % cell-wise errors
    in those few cells. On 1024 cells, p99 = top 10 cells, which sits
    inside that smearing zone -- so p99 stays at 25 % even when the bulk
    wave structure is right to <1 % (LTE pureH collision compressions of
    ~10x consistently show this). p95 = top 50 cells, which lands just
    outside the smearing zone and captures the bulk-flow accuracy that
    is actually diagnostic.

    Empirical p95 on a correctly-behaving 1024-cell build (post t_final
    cleanup for tests 3-4 so shocks stay in domain):

      - Coleman 1/2/5/6: rho/eint p95 < 1 % across all EoS variants.
      - Coleman 3/4 (collisions): rho p95 up to 7 %, eint p95 up to 3 %.

    Default gates leave ~1 % headroom on the strongest cases.

    Velocity is gated on absolute error normalised by max(max|v_ref|,
    max(|v_L|, |v_R|)) so collision-test post-shock v ~= 0 doesn't
    inflate the relative error.

    p99 and the max are written to coleman_metrics.log for diagnostics
    but do NOT gate. The Snow two-gamma jump check is also diagnostic by
    default; pass --gate-snow X to enforce it."""
    fails = []
    if cmp['p95_rel_rho'] > max_rho:
        fails.append(f"rho p95 {cmp['p95_rel_rho']:.2%} > {max_rho:.2%}")
    if cmp['p95_abs_v_norm'] > max_v:
        fails.append(f"v p95 (abs/v_max) {cmp['p95_abs_v_norm']:.2%} > {max_v:.2%}")
    if cmp['p95_rel_eint'] > max_eint:
        fails.append(f"eint p95 {cmp['p95_rel_eint']:.2%} > {max_eint:.2%}")
    return (not fails, fails)
