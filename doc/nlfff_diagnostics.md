# Common NLFFF iteration diagnostics

The legacy MHD magnetofriction relaxation, weighted optimization
extrapolator, and Grad--Rubin extrapolator all write

```text
<base_filename>_nlfff_metrics.csv
```

with the same schema:

```text
iteration,CW_sin_theta,epsilon_force,epsilon_div,magnetic_energy
```

The diagnostics are evaluated on all active cells using the magnetic field in
AMRVAC code units:

```text
CW_sin_theta = integral(|J x B| / |B| dV) / integral(|J| dV)
epsilon_force^2 = integral(h^2 |J x B|^2 / |B|^2 dV) / integral(B^2 dV)
epsilon_div^2   = integral(h^2 |div B|^2 dV) / integral(B^2 dV)
magnetic_energy = 0.5 integral(B^2 dV)
h = cell-volume^(1/3)
```

`CW_sin_theta`, `epsilon_force`, and `epsilon_div` should approach zero for a
well-resolved force-free, solenoidal field. Magnetic energy is not a residual;
it records which energy branch the iteration approaches.

The meaning of `iteration` follows the method's natural outer update:

- legacy MFR: magnetofriction iteration;
- optimization: attempted optimization step;
- single-polarity Grad--Rubin: current-field update;
- self-consistent Grad--Rubin: complete positive/negative pair.

The common CSV is always written. Method-specific detailed histories are
disabled by default and can be enabled independently:

```fortran
mf_write_detailed_history=.true.      ! in mf_list
nlfff_write_detailed_history=.true.   ! Optimization usr_list
gr_write_detailed_history=.true.      ! Grad--Rubin usr_list
```

The corresponding public optimization and Grad--Rubin configuration types use
the common member name `write_detailed_history`. When enabled, the detailed
files remain the authoritative source for method-specific convergence and
control quantities: MFR artificial time step and flux metric, optimization
functional and accept/reject history, and Grad--Rubin field/energy change, P/N
discrepancy, alpha change, and field-line classifications.

## Limitations

The common file uses active cells only. The optimization-specific weighted
functional additionally includes its fixed lower node and taper weights, while
the legacy MFR log uses an inner-domain mask. Their method-specific values are
therefore not identical to the common active-volume diagnostics. Magnetic
energy is in AMRVAC code units and can only be compared directly when runs use
the same units and domain.
