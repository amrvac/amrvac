# Data-Driven Data-Constrained MHD

This template starts a fixed-bottom MHD simulation from a potential-field or
magnetofrictional snapshot.  Only the three magnetic components are inherited.
The initial density, pressure/energy, and velocity are rebuilt by `mod_usr.t`;
the velocity is reset to zero.

The normalization is shared with the potential and MF cases:

- `unit_length = 1e9 cm`
- `unit_temperature = 1e6 K`
- `unit_numberdensity = 1e9 cm^-3`

## Physics models

| `mhd_model` | Variables | Atmospheres | Extra physics |
|---|---:|---|---|
| `zero_beta` | 7 | uniform, coronal, or chromospheric density prescription | no pressure or gravity |
| `isothermal` | 7 | hydrostatic isothermal corona only | spherical solar gravity |
| `adiabatic` | 8 | hydrostatic corona/chromosphere or relaxed table | energy and gravity |
| `thermodynamic` | 8 | hydrostatic corona/chromosphere or relaxed table | energy, gravity, field-aligned conduction, cooling, heating |

The model-specific parameter files are in `par_zero_beta/`,
`par_isothermal/`, `par_adiabatic/`, and `par_thermodynamic/`.  The top-level
`amrvac.par` remains a compatibility copy of the zero-beta setup.

`hydrostatic + chromosphere` uses `mod_solar_atmosphere` and defaults to the
AL-C7 temperature curve with `n(10 Mm)=5e8 cm^-3`.  The coronal default is
`T=1 MK` and `n(z=0)=1e9 cm^-3`.  A relaxed atmosphere table must contain at
least three code-unit columns:

```text
# z  rho  p  [optional extra columns]
```

Its height must increase strictly and cover the physical vertical domain.
Density and pressure are logarithmically interpolated.

The thermodynamic background heating is
`H=H0 exp(-max(z,0)/lambda)`, with defaults `H0=1e-4 erg cm^-3 s^-1` and
`lambda=50 Mm`.  Conduction and radiative-cooling controls remain ordinary
AMRVAC settings in `par_thermodynamic/amrvac.par`.

## Staging and restart behavior

Use `stage_data_constrained_case(...)` to select the model and write the
data-specific override.  For example:

```python
stage_data_constrained_case(
    metadata,
    "DrivenFieldProject/DataConstrained",
    restart_file="DrivenFieldProject/MagnetofrictionalRelaxation/output/data_driven_mfr0000.dat",
    mhd_model="thermodynamic",
    atmosphere_model="chromosphere",
)
```

The generated initial override sets `firstprocess`, `reset_time`, and
`reset_it` to true.  Run the returned `initial_run_command` once.  For an
interrupted data-constrained run, restart from its own checkpoint with a new
override that omits those three settings; otherwise the evolved plasma would
be intentionally rebuilt and the clock reset.

The lower boundary is kept at the maximum AMR level.  Its magnetic field,
density, and pressure are fixed, and all three velocity components are zero.
Standard AMRVAC `.log`, `.dat`, and VTU output is retained.  Rank zero also
writes `output/initial_atmosphere.dat` with code-unit columns
`z rho p T gravity heating`.
