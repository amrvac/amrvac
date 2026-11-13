#!/usr/bin/env python3
"""Tabulate the rho convergence study written by run_convergence.sh.

Reads every conv_<setup>_<limiter>_<stepper>_<flux>_N<N>.log in the given
directory (default ../output/convergence), and for each scheme prints
N, L1, L2, Linf with the observed order between successive resolutions,
log2(E_coarse/E_fine) for doubled N, and the relative drift of mean(rho)
between the first and last log line. Prints a summary of the
observed orders at the finest pair, checks hll against tvdlf, and writes
convergence.csv and a four-panel summary figure convergence.png (needs
matplotlib: uv sync --extra analysis) into the same directory.
"""
import csv
import math
import pathlib
import re
import sys
from collections import defaultdict

outdir = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "../output/convergence")
pat = re.compile(r"conv_(\w+?)_(\w+?)_(\w+?)_(\w+?)_N(\d+)\.log$")

runs = defaultdict(dict)
for f in sorted(outdir.glob("conv_*_N*.log")):
    m = pat.search(f.name)
    if not m:
        continue
    rows = [list(map(float, l.split())) for l in f.read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    if len(rows) < 2:
        print(f"warning: {f.name} has {len(rows)} line(s), run incomplete?")
        continue
    t, l1, l2, linf, mass = rows[-1]
    drift = abs(mass - rows[0][4]) / abs(rows[0][4])
    runs[m.groups()[:4]][int(m.group(5))] = (t, l1, l2, linf, drift)

if not runs:
    sys.exit(f"no conv_*.log files in {outdir}")

def order(e_coarse, e_fine, n_coarse, n_fine):
    if e_coarse <= 0 or e_fine <= 0:
        return float("nan")
    return math.log(e_coarse / e_fine) / math.log(n_fine / n_coarse)

table = []
for key in sorted(runs):
    setup, lim, stepper, flux = key
    res = runs[key]
    ns = sorted(res)
    print(f"\n{setup:4s}  limiter={lim:8s} stepper={stepper:9s} flux={flux}")
    print(f"{'N':>5} {'L1':>11} {'ord':>5} {'L2':>11} {'ord':>5} "
          f"{'Linf':>11} {'ord':>5} {'mass drift':>11}")
    prev = None
    for n in ns:
        t, l1, l2, linf, drift = res[n]
        if prev is None:
            o1 = o2 = oi = float("nan")
        else:
            pn, (pt, pl1, pl2, plinf, _) = prev
            o1 = order(pl1, l1, pn, n)
            o2 = order(pl2, l2, pn, n)
            oi = order(plinf, linf, pn, n)
        f = lambda o: "    -" if math.isnan(o) else f"{o:5.2f}"
        print(f"{n:5d} {l1:11.3e} {f(o1)} {l2:11.3e} {f(o2)} "
              f"{linf:11.3e} {f(oi)} {drift:11.1e}")
        table.append([setup, lim, stepper, flux, n, t, l1, o1, l2, o2,
                      linf, oi, drift])
        prev = (n, res[n])

with open(outdir / "convergence.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["setup", "limiter", "stepper", "flux", "N", "time", "L1",
                "order_L1", "L2", "order_L2", "Linf", "order_Linf",
                "mass_drift"])
    w.writerows(table)
print(f"\nwrote {outdir / 'convergence.csv'}")

# --------------------------------------------------------------------------
# Findings: hll vs tvdlf, and the finest-pair orders per scheme
# --------------------------------------------------------------------------
def finest_orders(res):
    ns = sorted(res)
    if len(ns) < 2:
        return None
    a, b = ns[-2], ns[-1]
    return tuple(order(res[a][i], res[b][i], a, b) for i in (1, 2, 3)) + (a, b)

print("\n== summary: orders between the two finest resolutions")
print(f"{'setup':5s} {'limiter':8s} {'stepper':9s} {'flux':5s} {'N':>9s} "
      f"{'L1':>5s} {'L2':>5s} {'Linf':>5s}")
for key in sorted(runs):
    fo = finest_orders(runs[key])
    if fo:
        print(f"{key[0]:5s} {key[1]:8s} {key[2]:9s} {key[3]:5s} "
              f"{fo[3]:4d}-{fo[4]:<4d} {fo[0]:5.2f} {fo[1]:5.2f} {fo[2]:5.2f}")

# hll and tvdlf both reduce to exact upwinding for linear advection, so they
# must agree to round-off. onestep (Forward Euler with a second-order
# reconstruction) is not linearly stable and is kept bounded only by the
# limiter, so it amplifies round-off differences: report it separately.
hll_diff, hll_diff_fe = [], []
for key in runs:
    if key[3] == "hll" and (key[0], key[1], key[2], "tvdlf") in runs:
        other = runs[(key[0], key[1], key[2], "tvdlf")]
        dest = hll_diff_fe if key[2] == "onestep" else hll_diff
        for n, v in runs[key].items():
            if n in other and other[n][1] > 0:
                dest.append(abs(v[1] - other[n][1]) / other[n][1])
if hll_diff:
    print(f"\nhll vs tvdlf: max relative difference in L1 over "
          f"{len(hll_diff)} twostep/threestep runs = {max(hll_diff):.1e}")
if hll_diff_fe:
    print(f"hll vs tvdlf, onestep only (round-off amplified by Forward "
          f"Euler): max {max(hll_diff_fe):.1e} over {len(hll_diff_fe)} runs")
drifts = [v[4] for r in runs.values() for v in r.values()]
print(f"mass drift: max over all {len(drifts)} runs = {max(drifts):.1e}")

# --------------------------------------------------------------------------
# Summary figure
# --------------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("\nmatplotlib not available, no plot: install it with\n"
          "  cd $AGILE_DIR && uv sync --extra analysis")
    sys.exit(0)

LIMS = ["minmod", "vanleer", "mcbeta", "koren"]
STEPS = ["onestep", "twostep", "threestep"]
COL = {"minmod": "#1b9e77", "vanleer": "#d95f02", "mcbeta": "#7570b3",
       "koren": "#e7298a"}
STEPLAB = {"onestep": "onestep (Forward Euler)",
           "twostep": "twostep (predictor-corrector)",
           "threestep": "threestep (ssprk3)"}
STEPSTY = {"onestep": ":", "twostep": "--", "threestep": "-"}


def pick(setup, lim, stepper, flux="tvdlf"):
    return runs.get((setup, lim, stepper, flux))


def slopes(ax, ns, eref, orders=(1, 2, 3)):
    """Reference lines of given orders, anchored at the coarsest point."""
    n0, n1 = min(ns), max(ns)
    for p, ls in zip(orders, (":", "--", "-.")):
        ax.loglog([n0, n1], [eref, eref * (n0 / n1) ** p], color="0.55",
                  ls=ls, lw=0.9)
        ax.text(n1 * 1.03, eref * (n0 / n1) ** p, f"$N^{{-{p}}}$",
                fontsize=8, color="0.4", va="center")


def curve(ax, res, label, color, ls="-", marker="o", idx=1):
    """Plot one error curve; the order at the finest pair goes in the legend."""
    ns = sorted(res)
    es = [res[n][idx] for n in ns]
    if len(ns) > 1:
        label += f"  ({order(es[-2], es[-1], ns[-2], ns[-1]):.2f})"
    ax.loglog(ns, es, ls=ls, marker=marker, color=color, label=label,
              ms=4, lw=1.4)
    return ns, es


def finish(ax, title, ns_all, e_all):
    if not ns_all:
        ax.text(0.5, 0.5, "no runs", ha="center", va="center",
                transform=ax.transAxes, color="0.5")
        ax.set_title(title, fontsize=10)
        return
    slopes(ax, ns_all, max(e_all) * 1.6)
    ax.set_xscale("log", base=2)
    ax.set_xticks(sorted(set(ns_all)))
    ax.set_xticklabels([str(n) for n in sorted(set(ns_all))])
    ax.minorticks_off()
    ax.set_xlim(min(ns_all) / 1.25, max(ns_all) * 1.45)
    ax.set_xlabel("N (cells per direction)")
    ax.set_title(title, fontsize=10)
    ax.grid(True, which="major", lw=0.3, color="0.85")
    ax.legend(fontsize=7, loc="lower left", framealpha=0.9,
              title="(order at finest pair)", title_fontsize=7)


fig, axes = plt.subplots(2, 2, figsize=(12, 9.5))

# (a) limiters on the quasi-1D wave
ax = axes[0, 0]
ns_all, e_all = [], []
for lim in LIMS:
    r = pick("x", lim, "threestep")
    if r:
        n, e = curve(ax, r, lim, COL[lim])
        ns_all += n; e_all += e
ax.set_ylabel("L1 error")
finish(ax, "(a) limiters, x wave (quasi-1D), ssprk3, tvdlf", ns_all, e_all)

# (b) time steppers on the quasi-1D wave, for vanleer and koren
ax = axes[0, 1]
ns_all, e_all = [], []
for lim in ("vanleer", "koren"):
    for st in STEPS:
        r = pick("x", lim, st)
        if r:
            n, e = curve(ax, r, f"{lim}, {STEPLAB[st]}", COL[lim],
                         ls=STEPSTY[st])
            ns_all += n; e_all += e
ax.set_ylabel("L1 error")
finish(ax, "(b) time steppers, x wave, tvdlf", ns_all, e_all)

# (c) the fully 3D diagonal wave: L1 solid, Linf dashed
ax = axes[1, 0]
ns_all, e_all = [], []
for lim in LIMS:
    r = pick("diag", lim, "threestep")
    if r:
        n, e = curve(ax, r, f"{lim}  L1", COL[lim])
        n2, e2 = curve(ax, r, f"{lim}  Linf", COL[lim], ls="--",
                       marker="s", idx=3)
        ns_all += n + n2; e_all += e + e2
ax.set_ylabel("error")
finish(ax, "(c) diag wave (3D), ssprk3, tvdlf: L1 solid, Linf dashed",
       ns_all, e_all)

# (d) observed orders at the finest pair, every tvdlf scheme
ax = axes[1, 1]
rows = []
for setup in ("x", "diag"):
    for st in STEPS:
        for lim in LIMS:
            r = pick(setup, lim, st)
            fo = finest_orders(r) if r else None
            if fo:
                rows.append((f"{setup}  {lim}  {st}", fo))
if rows:
    import numpy as np
    y = np.arange(len(rows))
    h = 0.27
    for i, (nm, c) in enumerate((("L1", "#4c72b0"), ("L2", "#55a868"),
                                 ("Linf", "#c44e52"))):
        ax.barh(y + (i - 1) * h, [fo[i] for _, fo in rows], height=h,
                color=c, label=nm)
    for p in (1, 2, 3):
        ax.axvline(p, color="0.5", lw=0.8, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels([nm for nm, _ in rows], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 3.3)
    ax.set_xlabel("observed order (two finest resolutions)")
    ax.legend(fontsize=7, loc="lower right")
else:
    ax.text(0.5, 0.5, "no runs", ha="center", va="center",
            transform=ax.transAxes, color="0.5")
ax.set_title("(d) observed convergence orders, tvdlf", fontsize=10)

notes = []
if hll_diff:
    notes.append(f"hll = tvdlf to {max(hll_diff):.0e} (twostep/threestep)")
notes.append(f"max mass drift {max(drifts):.0e}")
fig.suptitle("phys = 'rho': convergence of a smooth periodic wave after one "
             "period   (" + ";  ".join(notes) + ")", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(outdir / "convergence.png", dpi=130)
print(f"\nwrote {outdir / 'convergence.png'}")
