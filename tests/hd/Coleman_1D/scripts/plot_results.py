#!/usr/bin/env python3
"""Per-test profile plots for the Coleman_1D autotest.

For every (test, variant) that has both a sim (.dat) and a cached
reference (.npz), produces a 2x3 figure overlaying sim and reference for
rho, v, eint, T, Gamma_1, and rel-rho error, plus Snow markers where
applicable. Output to ../plots/coleman<N>_<variant>.png.

Modelled on Projects/lte_eos_type2/regression/run_regression.py:plot_one.
ASCII-only labels per the project rule (use spelled-out names instead of
glyphs in titles/annotations; LaTeX math in axis labels is fine, the file
content stays ASCII because backslash sequences are not multi-byte).

Usage from the test directory:
    python3 scripts/plot_results.py [--tests N ...] [--variants V ...]
    make plot   # convenience target
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
TEST_DIR   = SCRIPT_DIR.parent
OUT_DIR    = TEST_DIR
REF_DIR    = TEST_DIR / 'reference_riemann'
PLOTS_DIR  = TEST_DIR / 'plots'

sys.path.insert(0, str(SCRIPT_DIR))
from coleman_ics import TESTS                              # noqa: E402
from compare import read_dat_1d, snow_check                # noqa: E402
from eos_reader import make_eos                            # noqa: E402

VARIANTS = ('FI_HHe', 'FI_pureH', 'LTE_HHe', 'LTE_pureH')


def find_amrvac_dir(tables_subdir: str) -> Path | None:
    import os
    env = os.environ.get('AMRVAC_DIR')
    candidates = []
    if env:
        candidates.append(Path(env).resolve())
    candidates.append(SCRIPT_DIR.parent.parent.parent.parent.resolve())
    for d in candidates:
        if ((d / 'src' / 'eos').is_dir()
                and (d / 'src' / 'tables' / 'eos_tables' / tables_subdir).is_dir()):
            return d
    return None


def find_final_dat(name: str) -> Path | None:
    candidates = sorted(OUT_DIR.glob(f"{name}[0-9]*.dat"))
    return candidates[-1] if candidates else None


def plot_one(test_num: int, variant: str, eos, ic: dict) -> Path | None:
    name = f'coleman{test_num}_{variant}'
    dat_path = find_final_dat(name)
    ref_path = REF_DIR / f'{name}.npz'
    if dat_path is None or not ref_path.is_file():
        print(f'  [skip] {name}: missing dat or reference')
        return None

    dat = read_dat_1d(str(dat_path))
    ref = np.load(ref_path)
    x = dat['x']
    rho_sim = dat['rho']
    v_sim   = dat['m1'] / np.maximum(rho_sim, 1e-30)
    eint_sim = dat['e'] - 0.5 * rho_sim * v_sim ** 2

    x_ref = ref['x']
    rho_ref  = np.interp(x, x_ref, ref['rho'])
    v_ref    = np.interp(x, x_ref, ref['v'])
    p_ref    = np.interp(x, x_ref, ref['p'])
    eint_ref = np.interp(x, x_ref, ref['eint'])

    F = eos.nH2rhoFactor
    log_nH_ref = np.log10(np.maximum(rho_ref / F, 1e-30))
    log_eint_ref = np.log10(np.maximum(eint_ref / np.maximum(rho_ref / F, 1e-30), 1e-30))
    log_nH_sim = np.log10(np.maximum(rho_sim / F, 1e-30))
    log_eint_sim = np.log10(np.maximum(eint_sim / np.maximum(rho_sim / F, 1e-30), 1e-30))

    T_ref  = np.array([eos.temperature(log_nH_ref[i],  log_eint_ref[i]) for i in range(len(x))])
    g1_ref = np.array([eos.gamma1(log_nH_ref[i], log_eint_ref[i]) for i in range(len(x))])
    T_sim  = dat.get('Te', None)
    g1_sim = dat.get('g1', None)
    if T_sim is None:
        T_sim = np.array([eos.temperature(log_nH_sim[i], log_eint_sim[i]) for i in range(len(x))])

    snow_for_this = None
    try:
        snow_for_this = snow_check(dat, ic, eos, ic['shock_sides'],
                                    reference_path=ref_path)
    except Exception as exc:
        snow_for_this = {'error': str(exc)}

    fig, ax = plt.subplots(2, 3, figsize=(14, 7.5), sharex=True)
    panels = [
        (ax[0, 0], rho_ref,    rho_sim,    r'$\rho$',           True),
        (ax[0, 1], v_ref,      v_sim,      r'$v$',              False),
        (ax[0, 2], eint_ref,   eint_sim,   r'$e_{\rm int}$',    True),
        (ax[1, 0], T_ref,      T_sim,      r'$T$ [code]',       True),
    ]
    for axi, ye, ys, ylabel, logy in panels:
        axi.plot(x, ye, 'k-',  lw=1.4, label='exact Riemann')
        if ys is not None:
            axi.plot(x, ys, 'r--', lw=1.0, label='AMRVAC')
        if logy:
            axi.set_yscale('log')
        axi.set_ylabel(ylabel)
        axi.grid(alpha=0.3, which='both')

    ax_g = ax[1, 1]
    ax_g.plot(x, g1_ref, 'k-',  lw=1.4, label=r'$\Gamma_1$ (exact)')
    if g1_sim is not None:
        ax_g.plot(x, g1_sim, 'r--', lw=1.0, label=r'$\Gamma_1$ (AMRVAC)')
    ax_g.set_ylabel(r'$\Gamma_1$')
    ax_g.grid(alpha=0.3, which='both')
    ax_g.legend(fontsize=7, loc='best')

    ax_err = ax[1, 2]
    rel_err = np.abs(rho_sim - rho_ref) / np.maximum(np.abs(rho_ref), 1e-30)
    ax_err.semilogy(x, np.maximum(rel_err, 1e-16), 'r-', lw=1.0)
    for thr, lbl in [(1e-1, '10%'), (1e-2, '1%'), (1e-3, '0.1%')]:
        ax_err.axhline(thr, color='k', ls=':', lw=0.6, alpha=0.5)
        ax_err.text(x[-1], thr, lbl, fontsize=7, va='bottom', ha='right')
    ax_err.set_ylabel(r'$|\rho_{\rm sim}-\rho_{\rm exact}|/|\rho_{\rm exact}|$')
    ax_err.grid(alpha=0.3, which='both')

    ax[0, 0].legend(fontsize=8, loc='best')

    if snow_for_this and 'error' not in snow_for_this:
        for side, d in snow_for_this.items():
            if not isinstance(d, dict) or 'error' in d:
                continue
            # |grad v_ref| locates the shock cleanly for any EoS: v is
            # continuous across the contact but jumps across the shock.
            # |grad rho_ref| picks the wrong feature whenever the contact
            # jump exceeds the shock jump (every FI Sod-like setup).
            half = (x < 0.5) if side == 'L' else (x > 0.5)
            if half.any():
                idx_local = int(np.argmax(np.abs(np.gradient(v_ref[half]))))
                shock_x = x[np.where(half)[0][idx_local]]
            else:
                shock_x = 0.5
            # Snow post-shock predictions: rho_pre comes from the IC
            # (returned by snow_check), NOT from sampling the profile
            # (which by t_final is in the star state, not the pre-shock
            # state). The Gamma_1 closure ignores ionisation latent heat
            # and underpredicts compression; the gamma_eff closure
            # captures it via the lower p/eint ratio and typically
            # overlaps with the sim post-shock value.
            rho_pre = d['rho_pre']
            rho_post_g1   = rho_pre * d['r_snow_g1']
            rho_post_geff = rho_pre * d['r_snow_geff']
            ax[0, 0].axvline(shock_x, color='grey', ls=':', lw=0.6, alpha=0.5)
            ax[0, 0].plot([shock_x], [rho_post_g1], 's', color='C2', ms=8,
                          label=r'$\Gamma_1$ closure (adiabatic)'
                                if side == ic['shock_sides'][0] else None)
            ax[0, 0].plot([shock_x], [rho_post_geff], '*', color='C4', ms=11,
                          label=r'$\gamma_{\rm eff}$ closure (Snow)'
                                if side == ic['shock_sides'][0] else None)
            ax[1, 1].axvline(shock_x, color='grey', ls=':', lw=0.6, alpha=0.5)
            ax[1, 1].plot([shock_x], [d['g1_post']],   's', color='C2', ms=8)
            ax[1, 1].plot([shock_x], [d['geff_post']], '*', color='C4', ms=11)
        handles, labels = ax[0, 0].get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax[0, 0].legend(by_label.values(), by_label.keys(), fontsize=8, loc='best')

    for axi in ax[1, :]:
        axi.set_xlabel('x')

    fig.suptitle(f'{name}  t={dat["time"]:.4f}  '
                 f'(eos_type={ic.get("description", "")})', fontsize=11)
    fig.tight_layout()
    PLOTS_DIR.mkdir(exist_ok=True)
    out = PLOTS_DIR / f'{name}.png'
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out


def variant_props(variant: str) -> dict:
    eos_type, comp = variant.split('_')
    He = 0.1 if comp == 'HHe' else 0.0
    return {'eos_type': eos_type, 'He_abundance': He, 'composition': comp}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--tests',    nargs='+', type=int, default=sorted(TESTS))
    p.add_argument('--variants', nargs='+', default=list(VARIANTS), choices=list(VARIANTS))
    p.add_argument('--tables-subdir', default='entropy128')
    args = p.parse_args()

    amrvac_dir = find_amrvac_dir(args.tables_subdir)
    if amrvac_dir is None:
        print('ERROR: cannot locate AMRVAC tree; export AMRVAC_DIR.')
        return 1
    tables_dir = str(amrvac_dir / 'src' / 'tables' / 'eos_tables' / args.tables_subdir)

    n_written = 0
    for n in args.tests:
        if n not in TESTS:
            continue
        ic = TESTS[n]
        for variant in args.variants:
            v = variant_props(variant)
            try:
                eos = make_eos(v['eos_type'], tables_dir, v['He_abundance'])
            except Exception as exc:
                print(f'  [skip] coleman{n}_{variant}: EoS load failed: {exc}')
                continue
            out = plot_one(n, variant, eos, ic)
            if out is not None:
                print(f'  -> {out.relative_to(TEST_DIR)}')
                n_written += 1
    print(f'\n[summary] wrote {n_written} plots to {PLOTS_DIR.relative_to(TEST_DIR)}/')
    return 0


if __name__ == '__main__':
    sys.exit(main())
