# UAWSoM equation and implementation map

This note records the scientific and code provenance of the UAWSoM extension
before it is integrated into the MHD module.  The scientific authority is
McMurdo et al. (2026), A&A 705, A15,
[doi:10.1051/0004-6361/202555912](https://doi.org/10.1051/0004-6361/202555912).
The compatibility reference is commit `932ee5d` of
`maxmcmurdo/Solar_MHD_2026_tutorial_UAWSoM_in_MPIAMRVAC`.

| Paper | Fork implementation | MPI-AMRVAC 4.0 target |
| --- | --- | --- |
| Eq. 3, wave pressure in momentum | `uawsom_get_flux*`, `ptotal` | `mhd_get_flux*`, `mhd_uawsom_get_coefficients` |
| Eq. 4, total energy | `uawsom_get_flux*`, `w_add_source` | MHD total-energy flux and `mhd_add_source_uawsom` |
| Eq. 5, Alfven wave energy | `wAplus_`, `wAminus_` flux/source branches | optional MHD flux variables `wAplus`, `wAminus` |
| Eq. 6, kink wave energy | `wkplus_`, `wkminus_` flux/source branches | optional MHD flux variables `wkplus`, `wkminus` |
| Eq. 7, wave pressures | inline expressions in `uawsom_get_flux*` | `mhd_uawsom_wave_pressure_cell` |
| Eq. 8, averaged density | inline `zeta` and filling-factor expressions | coefficient callback and namelist defaults |
| Eq. 9, Alfven dissipation | `Gamma_plus`, `Gamma_minus` in `w_add_source` | `mhd_add_source_uawsom` |
| Eq. 10, kink correlation length | `Lperp` in `w_add_source` | coefficient callback/default coefficient evaluator |
| Eq. 34, Alfven reflection | commented tutorial source | runtime `mhd_uawsom_reflection` branch |

## Deliberate differences from the tutorial fork

The total-energy flux follows Eqs. 4--7 directly.  Because the conserved total
energy already contains the four wave energies, only the wave-pressure work
and field-aligned propagation parts are added to the ordinary MHD energy flux.
The tutorial fork adds the complete wave flux a second time and therefore
duplicates its bulk-velocity part.

The kink expansion-work source follows the positive right-hand-side sign in
paper Eq. 4.  The tutorial fork applies this term with the opposite sign; that
fork behavior is not retained.

The tutorial fork hard-codes the second coordinate, solar radius, filling
factor, thread radius, reference field, and correlation lengths.  The merged version
uses Cartesian dimension-independent operators and runtime parameters.  A
user callback can provide local density contrast, thread radius, and Alfven
correlation length.  The supplied solar-atmosphere example implements the
paper profiles through that callback.

The solar-atmosphere test follows the paper configuration: `zeta0=5`,
`B0=20 G`, and `R0=1 Mm`.  It returns Eq. 30's
`zeta=1+(zeta0-1) exp(-z/(5 R_sun))`, Eq. 29's magnetic expansion of the
thread radius, and the temperature-dependent closure stated after Eq. 9,
`Lperp_AW=1e8 sqrt(T_MK/B_G) m`.  The callback receives whether its state is
primitive or conserved so that both solver paths evaluate the same physical
temperature.  The fork's base `zeta=6`, `R0=0.1 Mm`, and temperature-independent
Alfvén correlation length are not retained.

Alfven reflection is an optional, conservative exchange between the two
propagation populations using the local Eq. 34 coefficient.  Kink reflection
remains zero, matching the stated limitation of McMurdo et al. (2026).

Before upstream merge, this table and the two explicitly recorded tutorial
differences (energy-flux double counting and kink-work sign) are the author
confirmation checklist.
