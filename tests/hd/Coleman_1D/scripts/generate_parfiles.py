#!/usr/bin/env python3
"""Generate the 24 Coleman 1D test parfiles (6 tests x {FI,LTE} x {H,HHe}).

Run once when the test matrix or per-test parameters change. The output
lives in ../parfiles/ alongside this script.
"""

from __future__ import annotations
import sys
from pathlib import Path

# Single source of truth lives in coleman_ics.py (TESTS dict). Adapt the
# layout here: coleman_ics uses t_final, this script uses time_max plus a
# per-test dtsave_log spacing for the regression_test log.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from coleman_ics import TESTS as _IC_TESTS                              # noqa: E402

# dtsave_log: 50 log entries per test, rounded to a clean value. Adjusted
# per test so the column-stable regression_test log keeps similar length
# across tests with very different t_final.
_DTSAVE_LOG = {1: 0.01, 2: 0.01, 3: 0.03, 4: 0.02, 5: 0.01, 6: 0.01}

TESTS = {
    n: dict(rho_L=ic['rho_L'], rho_R=ic['rho_R'],
            v_L=ic['v_L'],     v_R=ic['v_R'],
            T_L=ic['T_L'],     T_R=ic['T_R'],
            time_max=ic['t_final'],
            dtsave_log=_DTSAVE_LOG.get(n, 0.01),
            description=ic['description'])
    for n, ic in _IC_TESTS.items()
}

VARIANTS = {
    'FI_HHe':   dict(eos_type='FI',  He_abundance=0.1, composition='HHe'),
    'FI_pureH': dict(eos_type='FI',  He_abundance=0.0, composition='H'),
    'LTE_HHe':   dict(eos_type='LTE', He_abundance=0.1, composition='HHe'),
    'LTE_pureH': dict(eos_type='LTE', He_abundance=0.0, composition='H'),
}


def eos_block(variant: dict) -> str:
    if variant['eos_type'] == 'FI':
        return (
            "&eos_list\n"
            "  eos_type      = 'FI'\n"
            "  gamma         = 1.6666666666666667d0\n"
            f"  He_abundance  = {variant['He_abundance']:.1f}d0\n"
            "/\n"
        )
    # LTE -- composition is auto-derived from He_abundance (HHe vs H).
    # eos_method = 'entropy' selects the bicubic-Hermite entropy128 path;
    # without it, AMRVAC defaults to the legacy uniform256 tables and the
    # cached Riemann reference (built against entropy128) will disagree.
    # The relative path resolves from tests/hd/Coleman_1D/ -> src/tables/...
    return (
        "&eos_list\n"
        "  eos_type      = 'LTE'\n"
        "  ionE          = .true.\n"
        "  table_check   = .false.\n"
        "  eos_method    = 'entropy'\n"
        "  table_location = '../../../src/tables/eos_tables/entropy128/'\n"
        f"  He_abundance  = {variant['He_abundance']:.1f}d0\n"
        "/\n"
    )


PAR_TEMPLATE = """!setup.pl -d=1
!> Coleman ({test_n}) -- {variant_name}
!> {description}

&filelist
  base_filename = 'coleman{test_n}_{variant_name}'
  saveprim      = .true.
  convert_type  = 'vtuBCCmpi'
  autoconvert   = .true.
  typefilelog   = 'regression_test'
/

&savelist
  itsave(1,1)  = 0
  itsave(1,2)  = 0
  dtsave_log   = {dtsave_log:.5g}d0
  dtsave_dat   = {time_max:.5f}d0
/

&stoplist
  time_max = {time_max:.5f}d0
/

&methodlist
  time_stepper    = 'threestep'
  flux_scheme     = 20*'hllc'
  limiter         = 20*'ppm'
  typeboundspeed  = 'pvrs'
/

&boundlist
  typeboundary_min1 = 3*'cont'
  typeboundary_max1 = 3*'cont'
/

&meshlist
  refine_max_level  = 1
  block_nx1         = 64
  domain_nx1        = 1024
  xprobmin1         = 0.0d0
  xprobmax1         = 1.0d0
/

&paramlist
  typecourant = 'maxsum'
  courantpar  = 0.4d0
/

&hd_list
  hd_energy          = .true.
  hd_gravity         = .false.
  eq_state_units     = .false.
/

{eos_block}
&usr_list
  rho_L = {rho_L:.7e}
  rho_R = {rho_R:.7e}
  v_L   = {v_L:.7e}
  v_R   = {v_R:.7e}
  T_L   = {T_L:.7e}
  T_R   = {T_R:.7e}
/
"""


def main():
    # Write parfiles at the test-directory root so the standard autotest
    # framework can find them (TESTS list in test.make references them
    # by bare name, e.g. coleman1_FI_HHe.par as the dependency for
    # coleman1_FI_HHe.log). The Python harness also reads parfiles from
    # the test root.
    here = Path(__file__).resolve().parent
    out_dir = here.parent
    n_written = 0
    for n, ic in TESTS.items():
        for variant_name, variant in VARIANTS.items():
            par_text = PAR_TEMPLATE.format(
                test_n=n,
                variant_name=variant_name,
                description=ic['description'],
                time_max=ic['time_max'],
                dtsave_log=ic['dtsave_log'],
                rho_L=ic['rho_L'], rho_R=ic['rho_R'],
                v_L=ic['v_L'], v_R=ic['v_R'],
                T_L=ic['T_L'], T_R=ic['T_R'],
                eos_block=eos_block(variant),
            )
            (out_dir / f'coleman{n}_{variant_name}.par').write_text(par_text)
            n_written += 1
    print(f'Wrote {n_written} parfiles to {out_dir}')


if __name__ == '__main__':
    main()
