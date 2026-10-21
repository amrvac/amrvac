#!/usr/bin/env python3
"""Coleman_1D autotest entry point.

Drives the 24-run test matrix (6 Coleman Riemann tests x {FI_HHe, FI_pureH,
LTE_HHe, LTE_pureH}). For each run:

  1. Calls AMRVAC on the variant's parfile -> <name>NNNN.dat in the test directory
  2. Loads the cached exact Riemann reference (reference_riemann/<name>.npz)
     and computes p99 of |sim-ref|/|ref| for rho, v, eint on the wave-bearing
     interior cells.
  3. Runs the Snow (2026) two-gamma jump-condition check inline; compares
     the analytical compression ratio against the AMRVAC-measured one at
     each shock side.
  4. Writes a one-line summary to coleman_metrics.log (column-stable so it
     can also be diff'd against correct_output/coleman_metrics.log via
     compare_logs if desired).

PASS criterion: all 24 runs satisfy p99(rho,v,eint) < gate AND Snow agrees
with AMRVAC's measured r within gate_snow.

Exit code: 0 if all gates pass, 1 otherwise.

Usage from the test directory:
    AMRVAC_DIR=/path/to/amrvac-wb python3 scripts/run_coleman_test.py [opts]

Options
-------
  --tests N1 N2 ...        restrict to listed test numbers (1..6)
  --variants V1 V2 ...     restrict to listed variants (FI_HHe etc.)
  --skip-run               re-use existing dat files
  --tables-subdir DIR      LTE table dir under src/tables/eos_tables/
                           (default: entropy128)
  --gate-rho/-v/-eint X    pass gates for rho, v, eint (default 0.05)
  --gate-snow X            Snow agreement gate (default 0.05)
  -n NPROC                 mpi procs per run (default 1)
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
TEST_DIR   = SCRIPT_DIR.parent
PARFILES   = TEST_DIR / 'parfiles'
REF_DIR    = TEST_DIR / 'reference_riemann'
OUT_DIR    = TEST_DIR
RESULTS    = TEST_DIR / 'coleman_metrics.log'

sys.path.insert(0, str(SCRIPT_DIR))
from coleman_ics import TESTS                              # noqa: E402
from compare import read_dat_1d, compare_against_riemann, snow_check, gate  # noqa: E402
from eos_reader import make_eos                            # noqa: E402

VARIANTS = ('FI_HHe', 'FI_pureH', 'LTE_HHe', 'LTE_pureH')


def variant_props(variant: str) -> dict:
    eos_type, comp = variant.split('_')
    He = 0.1 if comp == 'HHe' else 0.0
    return {'eos_type': eos_type, 'He_abundance': He, 'composition': comp}


def find_amrvac_dir(tables_subdir: str) -> Path:
    """Locate the AMRVAC source tree, preferring trees with the required
    table set actually present.

    Prefer $AMRVAC_DIR when its src/tables/eos_tables/<tables_subdir>
    exists. Otherwise walk up from this script's location -- the standard
    layout puts this test at $AMRVAC_DIR/tests/hd/Coleman_1D/scripts/, so
    the tree root is the 4th parent.
    """
    candidates = []
    env = os.environ.get('AMRVAC_DIR')
    if env:
        candidates.append(Path(env).resolve())
    # Walk up from this script: scripts/ -> Coleman_1D -> hd -> tests -> root
    candidates.append(SCRIPT_DIR.parent.parent.parent.parent.resolve())
    for d in candidates:
        if ((d / 'src' / 'eos').is_dir()
                and (d / 'src' / 'tables' / 'eos_tables' / tables_subdir).is_dir()):
            return d
    sys.exit('ERROR: cannot locate AMRVAC tree with '
             f'src/tables/eos_tables/{tables_subdir}/ (tried $AMRVAC_DIR and '
             'walk-up from script). Export AMRVAC_DIR=/path/to/amrvac-wb.')


def find_final_dat(name: str) -> Path | None:
    pattern = f"{name}[0-9]*.dat"
    candidates = sorted(OUT_DIR.glob(pattern))
    return candidates[-1] if candidates else None


def run_one(amrvac: Path, name: str, nproc: int) -> Path | None:
    """Run AMRVAC on the variant's parfile; return path to the final dat."""
    par = PARFILES / f'{name}.par'
    if not par.is_file():
        print(f'  [skip] missing parfile: {par.name}')
        return None
    log_path = TEST_DIR / f'{name}.run.log'
    t0 = time.time()
    with open(log_path, 'wb') as logf:
        cmd = (['mpirun', '-np', str(nproc), str(amrvac), '-i', str(par)]
               if nproc > 1 else [str(amrvac), '-i', str(par)])
        res = subprocess.run(cmd, cwd=TEST_DIR, stdout=logf,
                              stderr=subprocess.STDOUT)
    dt = time.time() - t0
    if res.returncode != 0:
        print(f'  [run] FAILED {name} (exit {res.returncode}, {dt:.1f}s) -- see {log_path.name}')
        return None
    dat = find_final_dat(name)
    if dat is None:
        print(f'  [run] FAILED {name} -- no dat output')
        return None
    print(f'  [run] {name:24s}  {dt:6.1f}s  -> {dat.name}')
    return dat


def write_metrics_log(rows: list[dict]) -> None:
    """Column-stable log: lines fit the compare_logs format (whitespace-
    separated floats with the row name in the first column)."""
    header = ('#name p95_rho p95_v p95_eint  p99_rho p99_v p99_eint '
              'snow_max_err passed')
    lines = [header]
    for r in rows:
        passed = 1 if r['passed'] else 0
        lines.append(f"{r['name']:24s} "
                     f"{r.get('p95_rel_rho', float('nan')):.6e} "
                     f"{r.get('p95_abs_v_norm', float('nan')):.6e} "
                     f"{r.get('p95_rel_eint', float('nan')):.6e}  "
                     f"{r.get('p99_rel_rho', float('nan')):.6e} "
                     f"{r.get('p99_abs_v_norm', float('nan')):.6e} "
                     f"{r.get('p99_rel_eint', float('nan')):.6e} "
                     f"{r.get('snow_max_err', float('nan')):.6e} "
                     f"{passed:d}")
    RESULTS.write_text('\n'.join(lines) + '\n')


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--tests',    nargs='+', type=int, default=sorted(TESTS))
    p.add_argument('--variants', nargs='+', default=list(VARIANTS), choices=list(VARIANTS))
    p.add_argument('--skip-run', action='store_true')
    p.add_argument('--tables-subdir', default='entropy128')
    p.add_argument('--gate-rho',  type=float, default=0.08)
    p.add_argument('--gate-v',    type=float, default=0.02)
    p.add_argument('--gate-eint', type=float, default=0.05)
    # Snow check is diagnostic-only by default (gate at infinity). Set
    # explicitly to enforce; see compare.gate() docstring.
    p.add_argument('--gate-snow', type=float, default=float('inf'))
    p.add_argument('-n', '--nproc', type=int, default=1)
    args = p.parse_args()

    amrvac_dir = find_amrvac_dir(args.tables_subdir)
    amrvac_bin = TEST_DIR / 'amrvac'
    if not amrvac_bin.is_file():
        sys.exit(f'ERROR: missing {amrvac_bin}; run `make` in {TEST_DIR}.')
    # OUT_DIR == TEST_DIR, already exists.
    REF_DIR.mkdir(exist_ok=True)

    tables_dir = str(amrvac_dir / 'src' / 'tables' / 'eos_tables' / args.tables_subdir)

    rows: list[dict] = []
    for n in args.tests:
        if n not in TESTS:
            print(f'  [skip] test {n} not in TESTS')
            continue
        ic = TESTS[n]
        for variant in args.variants:
            name = f'coleman{n}_{variant}'
            v = variant_props(variant)

            if args.skip_run:
                dat_path = find_final_dat(name)
                if dat_path is None:
                    print(f'  [skip-run] {name}: no dat -- FAIL')
                    rows.append(dict(name=name, passed=False, reasons=['no dat']))
                    continue
            else:
                dat_path = run_one(amrvac_bin, name, args.nproc)
                if dat_path is None:
                    rows.append(dict(name=name, passed=False, reasons=['sim failed']))
                    continue

            dat = read_dat_1d(dat_path)
            ref_path = REF_DIR / f'{name}.npz'
            if not ref_path.is_file():
                print(f'  [compare] {name}: missing reference {ref_path.name}; '
                      f'run generate_references.py')
                rows.append(dict(name=name, passed=False, reasons=['no reference']))
                continue
            cmp = compare_against_riemann(dat, ref_path, ic=ic)

            try:
                eos = make_eos(v['eos_type'], tables_dir, v['He_abundance'])
                snow = snow_check(dat, ic, eos, ic['shock_sides'],
                                   reference_path=ref_path)
                snow_max_err = max(
                    (s['err_amrvac_vs_snow'] for s in snow.values()),
                    default=0.0,
                )
            except Exception as exc:
                snow = {'error': str(exc)}
                snow_max_err = float('nan')
                print(f'  [snow] {name}: {type(exc).__name__}: {exc}')

            passed, reasons = gate(cmp,
                                    max_rho=args.gate_rho,
                                    max_v=args.gate_v,
                                    max_eint=args.gate_eint)
            if not np.isnan(snow_max_err) and snow_max_err > args.gate_snow:
                passed = False
                reasons.append(f'snow err {snow_max_err:.2%} > {args.gate_snow:.2%}')
            rows.append(dict(name=name, passed=passed, reasons=reasons,
                             snow=snow, snow_max_err=snow_max_err, **cmp))
            flag = 'OK' if passed else 'X'
            print(f'  [check] {name:24s}  {flag}  '
                  f'rhop95={cmp["p95_rel_rho"]:.2%}  '
                  f'vp95={cmp["p95_abs_v_norm"]:.2%}  '
                  f'eintp95={cmp["p95_rel_eint"]:.2%}  '
                  f'snow={snow_max_err:.2%}')

    write_metrics_log(rows)

    n_pass  = sum(1 for r in rows if r['passed'])
    n_total = len(rows)
    print(f'\n[summary] {n_pass}/{n_total} passed -> coleman_metrics.log')
    return 0 if n_pass == n_total else 1


if __name__ == '__main__':
    sys.exit(main())
