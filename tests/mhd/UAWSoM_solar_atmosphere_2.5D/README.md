# UAWSoM solar atmosphere (2.5D)

This is the full tutorial case migrated from UAWSoM fork commit `932ee5d`.
It uses the optional UAWSoM equations inside `mod_mhd`, B0 splitting, solar
gravity, anisotropic thermal conduction, radiative cooling, TRAC, and kink-wave
heating.  The case-specific `usr_uawsom_coefficients` callback follows McMurdo
et al. (2026): `zeta0=5`, `B0=20 G`, `R0=1 Mm`,
`zeta=1+(zeta0-1) exp[-z/(5 R_sun)]`, and
`Lperp_AW=1e8 sqrt(T_MK/B_G) m`.  The inconsistent legacy-fork values are not
retained; no solar radius or Cartesian direction is hard-coded in generic MHD.

The full `amrvac.par` retains the tutorial grid and duration and is intended for
manual release validation.  Its default `kink_amplitude=20` matches the
tutorial.  The short regression overrides this to `1e-4`, because the coarse
16-by-16 grid cannot resolve the full-amplitude boundary gradient.  Run the
CI-sized smoke configuration with:

```console
make -f test.make
```

Visualization is portable: `python vis_tool.py [data_dir] --snapshots 0 5 10`.
The default data directory is the case directory; use `--prefix` for renamed
snapshot series.

Scientific reference: McMurdo et al. (2026), A&A 705, A15,
doi:10.1051/0004-6361/202555912.
