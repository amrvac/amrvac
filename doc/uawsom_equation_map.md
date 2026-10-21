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
| Eq. 34, Alfven reflection | commented tutorial source | `one_dimensional_gradient` compatibility branch or Cartesian `cartesian_gradient_vorticity` branch |
| Cartesian gradient/vorticity closure | commented tutorial source around `Rimb`/`Rlim` | total-field `b dot grad(ln v_A)` plus `b dot curl(v)` |
| Kink reflection | commented `Rimbk`/`Rlimk` source | kink-speed gradient only; no field-aligned vorticity term |

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

## Reflection and propagation conventions

MPI-AMRVAC labels `wAplus` and `wkplus` as the populations propagating against
the local magnetic field, while `wAminus` and `wkminus` propagate along it.  In
a B0-split run every reflection speed uses the total field
`B = B_perturbation + B0`, not the perturbation alone.  Define

`b = B/|B|`, `v_A = |B|/sqrt(rho)`, and
`v_k = |B|/sqrt(rho_e (zeta+1)/2)`, with
`rho_e = rho/(1 + f zeta - f)`.  The Cartesian gradient/vorticity closure uses

`S_A = v_A b dot grad(ln v_A)`,
`Omega_A = b dot curl(v)`, and
`S_k = v_k b dot grad(ln v_k)`.

The Alfvén reflection limiter is bounded by the larger nonlinear Alfvén
damping rate, `R_imb,A = sqrt(S_A^2 + Omega_A^2)` and
`R_lim,A = min(R_imb,A, max(Gamma_plus, Gamma_minus))`.  Kink reflection uses
`R_lim,k = min(abs(S_k), max(Gamma_kplus, Gamma_kminus))`; it deliberately
omits `Omega_A`.  Neither multidimensional limiter contains
`mhd_uawsom_sigma`.
The 4:1 population-imbalance factor is zero between the two 4:1 thresholds,
and approaches a bounded signed value for stronger imbalance.  For
`W_plus >= 4 W_minus`, `F = 1 - 2 sqrt(W_minus/W_plus) >= 0`; therefore a
positive exchange `T = R F sqrt(W_plus W_minus)` removes energy from the plus
population and adds the same amount to the minus population.  For
`W_minus >= 4 W_plus`, `F <= 0`, so the updates reverse and weaken the minus
population.  Between the thresholds `F=0`.  Thus the helper always weakens
the dominant population, enhances the minor population, and leaves balanced
states unchanged.  The donor cap `abs(T) <= W_donor/dt` uses the actual
post-damping donor state; the equal-and-opposite updates conserve the exchanged
wave energy without silently repairing a negative state produced by the
ordinary sources.

`mhd_uawsom_reflection_mode='one_dimensional_gradient'` retains the original
one-dimensional gradient source and its sign convention; only this path uses
the positive `mhd_uawsom_sigma` multiplier.  The default remains
`one_dimensional_gradient` for backward compatibility.
`cartesian_gradient_vorticity` is the multidimensional Cartesian
gradient/vorticity closure and runs whenever Alfvén reflection is enabled,
including with `mhd_uawsom_sigma=0`.  `mhd_uawsom_kink_reflection=.true.`
enables the kink-speed-gradient exchange independently of both the Alfvén
switch and sigma.

The kink expansion-work source is the positive right-hand-side term in paper
Eq. 4: `+(zeta-1)/(zeta+1) p_k div(v)` in the gas-energy equation.  Its sign is
independent of the reflection exchange.

## Supported scope

The implementation supports uniform Cartesian 1D, 2D/2.5D, and 3D
total-energy MHD with the fixed-ionization EOS, including B0 splitting.  The
current task intentionally does not implement cylindrical, polar, or
spherical coordinates; semirelativistic, internal-energy, or
hydrodynamic-energy formulations; FLD; nonuniform grids; or angular-momentum
fix, `source_geom`/`source_geom_split`, or `angmomfix` combinations.

Before upstream merge, this table and the two explicitly recorded tutorial
differences (energy-flux double counting and kink-work sign) are the author
confirmation checklist.
