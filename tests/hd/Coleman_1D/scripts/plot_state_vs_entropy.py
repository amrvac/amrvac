#!/usr/bin/env python3
"""Overlay the 'state' (uniform256 direct-state tables) and 'entropy'
(entropy128 bicubic-Hermite) LTE solutions for every Coleman_1D case.

Both are LTE Saha-equilibrium EoS; the conserved dynamics should agree, with
the ionisation/temperature diagnostics differing at the table/interpolation
accuracy level. One 2x3 figure per (test, composition): rho, velocity,
pressure, Te, ne/nH, Gamma_1, plus a relative-rho panel. ASCII labels, no bold.

Run from the test dir:  python3 scripts/plot_state_vs_entropy.py
Output: ../plots/coleman<N>_<comp>_state_vs_entropy.png
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
TEST_DIR = SCRIPT_DIR.parent
PLOTS = TEST_DIR / "plots"
PLOTS.mkdir(exist_ok=True)
sys.path.insert(0, str(SCRIPT_DIR))
from compare import read_dat_1d  # noqa: E402

# The .dat stores CONSERVED variables (saveprim only affects the .vtu), so
# m1 is momentum (rho*v) and e is total energy. Derive velocity = m1/rho and
# internal energy eint = e - 0.5*m1^2/rho. Te/Ne/g1 are stored directly.
def derived(d):
    rho = d["rho"]
    out = {"rho": rho,
           "v": d["m1"] / rho,
           "eint": d["e"] - 0.5 * d["m1"] ** 2 / rho}
    for k in ("Te", "Ne", "g1"):
        if k in d:
            out[k] = d[k]
    return out

# (derived key, axis label)
PANELS = [("rho", "rho"), ("v", "velocity"), ("eint", "eint"),
          ("Te", "Te"), ("Ne", "Ne"), ("g1", "Gamma_1")]

for N in range(1, 7):
    for comp in ("HHe", "pureH"):
        ent_f = TEST_DIR / f"coleman{N}_LTE_{comp}0001.dat"
        st_f = TEST_DIR / f"coleman{N}_state_{comp}0001.dat"
        if not (ent_f.exists() and st_f.exists()):
            print(f"skip coleman{N}_{comp}: missing dat")
            continue
        x = read_dat_1d(ent_f)["x"]
        de = derived(read_dat_1d(ent_f))
        ds = derived(read_dat_1d(st_f))
        fig, axs = plt.subplots(2, 3, figsize=(13.5, 7.2))
        for ax, (key, lab) in zip(axs.flat, PANELS):
            ye, ys = de.get(key), ds.get(key)
            if ye is None or ys is None:
                ax.set_visible(False)
                continue
            ax.plot(x, ye, "-", color="C0", lw=1.3, label="entropy")
            ax.plot(x, ys, "--", color="C3", lw=1.3, label="state")
            ax.set_xlabel("x")
            ax.set_ylabel(lab)
            ax.legend(fontsize=8, frameon=False)
        fig.suptitle(f"Coleman {N} ({comp}): entropy vs state   t = "
                     f"{read_dat_1d(ent_f)['time']:.4f}")
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        out = PLOTS / f"coleman{N}_{comp}_state_vs_entropy.png"
        fig.savefig(out, dpi=110)
        plt.close(fig)
        # L1 profile norm (front-offset robust; max-pointwise is dominated by
        # Riemann-front offsets and Ne->0 regions).
        def L1(k):
            a, b = de.get(k), ds.get(k)
            if a is None:
                return float("nan")
            return np.sum(np.abs(b - a)) / max(np.sum(np.abs(a)), 1e-30)
        print(f"coleman{N}_{comp}: wrote {out.name}  "
              f"L1 rho={L1('rho'):.2e} v={L1('v'):.2e} eint={L1('eint'):.2e} "
              f"Te={L1('Te'):.2e} Ne={L1('Ne'):.2e} g1={L1('g1'):.2e}")
