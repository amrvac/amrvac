# B-only data-driven MHD

This template restarts from a PotentialField or MagnetofrictionalRelaxation
snapshot and linearly interpolates the three-component magnetic sequence at
the lower boundary.  It never reads `V_XXXX.dat`; velocity treatment belongs
to the AMRVAC user boundary implementation. The current B-only implementation
extrapolates the plasma velocity from the interior, fixes the boundary density
and pressure to their initial atmosphere, and imposes the
time-interpolated B in the innermost ghost layer. `driving_time_scale` defaults
to 12. Only the two frames bracketing the current observation time are kept in
memory.

`bottom_numberdensity_cm3` optionally renormalizes the complete initial
atmosphere so that its number density at the physical lower boundary has the
requested value. This preserves the atmosphere temperature profile and keeps
the initial interior and fixed ghost-cell boundary continuous. A value of `-1`
uses the atmosphere model's original normalization.

The four base parameter sets under `par_zero_beta`, `par_isothermal`,
`par_adiabatic`, and `par_thermodynamic` mirror the DataConstrained models.
The Python workflow selects one and writes `data_driven.par`, including the
portable boundary-series path and `driving_time_scale`.

For direct B driving, a coronal atmosphere is recommended. Chromospheric and
transition-region stratification remains available, but its steep density and
temperature gradients can make the driven lower boundary substantially less
stable; both staging and AMRVAC startup emit a warning when it is selected.
