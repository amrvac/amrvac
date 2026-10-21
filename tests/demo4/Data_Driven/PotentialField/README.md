# Data-Driven Potential Field

This demo initializes a Cartesian potential magnetic field from a Python V1
data-driven boundary frame.

Expected boundary-frame layout:

```text
snapshot_time, nx, ny, dx, dy, Bx, By, Bz
```

`dx` and `dy` are stored in km, and magnetic-field components are stored in
Gauss. The potential-field extrapolation uses the `Bz` component.

Two extrapolation methods are available through `usr_list`:

```fortran
potential_field_method='fft'   ! default
fft_padding_factor=2
fft_top_boundary='open'         ! 'open' or 'closed'
lalpha=0.d0                     ! nonzero for a constant-alpha FFT LFFF
lfff_flux_treatment='strict'    ! or 'subtract_mean'
lfff_max_flux_imbalance=0.1d0
```

The default `potential_field_method='fft'` uses the dependency-free spectral solver.
Set `potential_field_method='green'` to use the legacy Green-function solver.
The data-driven geometry places the magnetogram at `z=0` on the first lower
ghost-cell center. The physical lower face is therefore half a finest-cell
spacing above zero. The FFT solver uses the same half-cell source-plane depth,
and requires `llift=0`. A padding factor of 1 is a strictly periodic horizontal
solution; the default factor of 2 centers the magnetogram in a zero-padded
plane to reduce periodic-image effects.

For `fft_top_boundary='open'`, every nonzero mode decays into a half-space. A
constant-alpha LFFF in this mode requires a flux-balanced magnetogram and
`abs(lalpha)` smaller than the lowest nonzero wavenumber of the padded plane.
For `fft_top_boundary='closed'`, the nonzero modes satisfy `Bz=0` at the upper
physical face. A potential field retains the padded-domain mean flux as a
uniform `k=0` field through the top. A nonzero-alpha field still requires a
flux-balanced magnetogram and stops explicitly at a finite-height resonant
mode rather than regularizing a singular solution silently.

For a nonzero-alpha LFFF, `lfff_flux_treatment='strict'` keeps the default
behavior and rejects a relative bottom-flux imbalance above `1e-8`. The
opt-in `subtract_mean` mode first measures
`abs(sum(Bz))/sum(abs(Bz))` on the cropped physical magnetogram. It subtracts
the physical-core mean only when that value is no larger than
`lfff_max_flux_imbalance`; more strongly imbalanced or nearly unipolar maps
are rejected. The correction happens before FFT zero padding and never
overwrites the input boundary file. Potential fields (`lalpha=0`) are not
modified by this option.

This PotentialField case deliberately requires Cartesian 3D, a uniform
level-one mesh (`refine_max_level=1`), and cell-centered magnetic fields
(`stagger_grid=.false.`) for both extrapolation methods. The physical magnetogram core must match the level-one
horizontal cell centers exactly. Extra boundary ghost pixels are located and
cropped by their coordinates; the solver does not assume a fixed ghost width
and does not resample the data. Unsupported FFT sizes stop with their prime
factorization and a suggested supported size.

The case uses the standard data-driven coronal normalization
`unit_length=1e9 cm`, `unit_temperature=1e6 K`, and
`unit_numberdensity=1e9 cm^-3`. AMRVAC derives the magnetic-field unit from
these values. Density is uniformly one in code units, velocity is zero, the
energy equation and gravity are disabled, and no GLM `psi` variable is stored.

The example contains a base `amrvac.par` next to `mod_usr.t`. The Python
notebook can stage a user case from this directory and write a data-specific
override file, typically `data_driven_boundary.par`, containing the generated
`boundary_filename` and mesh extents. Run the staged case with both parameter
files, for example:

```bash
mpirun -np 4 ./amrvac -i amrvac.par data_driven_boundary.par
```

The self-contained FFT backend has a small direct-DFT and mixed-radix
round-trip test. The potential/LFFF vertical mode formulas have a separate
analytic test. Run them with:

```bash
bash run_fft_core_test.sh
bash run_lfff_fft_test.sh
```
