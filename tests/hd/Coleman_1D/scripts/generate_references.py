#!/usr/bin/env python3
"""One-time generator for the cached Riemann reference solutions.

Populates ../reference_riemann/<test>_<variant>.npz with the exact
(rho, v, p, eint) profile at t_final for each (test, variant) pair.

For FI variants this uses the bundled IdealGasEoS (no tables required).
For LTE variants this loads the entropy128 bicubic tables and uses the
full GeneralEoSRiemann solver -- so the LTE references depend on the
shipped table set and need to be regenerated if those tables change.

Run once when:
  - the IC table (coleman_ics.py) is edited
  - the LTE entropy128 tables are rebuilt
  - the autotest grid resolution changes (currently 1024 cells uniform)

Usage:
    AMRVAC_DIR=/path/to/amrvac-wb python3 scripts/generate_references.py [opts]

Options
-------
  --tests N1 N2 ...     restrict to listed test numbers
  --variants V1 V2 ...  restrict to listed variants
  --tables-subdir DIR   LTE table subdir (default: entropy128)
  --n-cells N           sample grid resolution (default: 1024 -- matches parfiles)
  --force               overwrite existing .npz files
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
TEST_DIR   = SCRIPT_DIR.parent
REF_DIR    = TEST_DIR / 'reference_riemann'

sys.path.insert(0, str(SCRIPT_DIR))
from coleman_ics import TESTS                      # noqa: E402
from eos_reader import IdealGasEoS                  # noqa: E402
from eos_riemann import solve_coleman_test          # noqa: E402
from entropy_eos import EntropyEoS                  # noqa: E402

VARIANTS = ('FI_HHe', 'FI_pureH', 'LTE_HHe', 'LTE_pureH')


def build_eos(variant: str, tables_dir: str):
    eos_type, comp = variant.split('_')
    He = 0.1 if comp == 'HHe' else 0.0
    composition = 'HHe' if comp == 'HHe' else 'H'
    if eos_type == 'FI':
        return IdealGasEoS(He_abundance=He)
    return EntropyEoS(tables_dir, composition=composition, He_abundance=He)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--tests',    nargs='+', type=int, default=sorted(TESTS))
    p.add_argument('--variants', nargs='+', default=list(VARIANTS), choices=list(VARIANTS))
    p.add_argument('--tables-subdir', default='entropy128')
    p.add_argument('--n-cells', type=int, default=1024)
    p.add_argument('--force', action='store_true')
    args = p.parse_args()

    amrvac_dir = os.environ.get('AMRVAC_DIR')
    if not amrvac_dir:
        print('WARNING: AMRVAC_DIR unset -- LTE variants will be skipped',
              file=sys.stderr)
        amrvac_dir = ''
    tables_dir = (Path(amrvac_dir) / 'src' / 'tables' / 'eos_tables'
                  / args.tables_subdir) if amrvac_dir else None

    REF_DIR.mkdir(exist_ok=True)
    x = np.linspace(0.0, 1.0, args.n_cells, endpoint=False) + 0.5 / args.n_cells

    n_written = n_skipped = n_skipped_existing = 0
    for n in args.tests:
        if n not in TESTS:
            continue
        ic = TESTS[n]
        for variant in args.variants:
            out_path = REF_DIR / f'coleman{n}_{variant}.npz'
            if out_path.is_file() and not args.force:
                n_skipped_existing += 1
                print(f'  [skip] {out_path.name} exists (use --force to overwrite)')
                continue
            if variant.startswith('LTE') and tables_dir is None:
                print(f'  [skip] {out_path.name}: AMRVAC_DIR required for LTE')
                n_skipped += 1
                continue
            eos = build_eos(variant, str(tables_dir) if tables_dir else '')
            print(f'\n[solving] coleman{n}_{variant}  t={ic["t_final"]:.4f}')
            rho, v, p, eint = solve_coleman_test(
                eos,
                ic['rho_L'], ic['rho_R'],
                ic['v_L'],   ic['v_R'],
                ic['T_L'],   ic['T_R'],
                x, ic['t_final'],
            )
            np.savez_compressed(out_path,
                                 x=x, rho=rho, v=v, p=p, eint=eint,
                                 t=ic['t_final'])
            n_written += 1
            print(f'  -> {out_path.name}  ({rho.size} cells)')
    print(f'\n[summary] wrote {n_written}, skipped {n_skipped + n_skipped_existing}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
