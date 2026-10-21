# Radiation Synthesis Validation Notes

This note records the validation scope for the synthetic-emission extensions in
`mod_thermal_emission`.

## Compatibility baseline

The default `&emissionlist` behavior remains unchanged:

- `radiation_transfer='thin'`
- `ray_method='auto'`
- `dat_resolution_mode='nominal'`
- `emission_model='auto'`
- `output_tau=.false.`
- `output_absorption_fraction=.false.`
- `instrument_postprocess=.false.`
- `radio_frequency=17.d9`
- `radio_beam_fwhm=0.d0`
- `radio_beam_pixel_size=0.d0`

Existing `convert_type='SI...'` and `WI...` inputs continue to use the
historical paths. EUV images use `ray_method='auto'` by default: Cartesian
dat-resolution EUV uses the native Cartesian ray tracer, spherical EUV uses
the native spherical ray/mesh-intersection path, and unsupported combinations
fall back to `legacy`. Set `ray_method='legacy'` to force the old projector.

## Current product matrix

The current implementation is intended to keep the historical workflow
(`convert_type` plus `&emissionlist`) while adding explicit transfer and ray
choices for Cartesian post-processing.

| Product | Legacy instrument grid | Cartesian dat-resolution | Cartesian `cart` arbitrary LOS | Thick transfer | Post-process to observational grid | Spherical |
| --- | --- | --- | --- | --- | --- | --- |
| EUV AIA/IRIS/EIS image | yes, thin | yes, axis-aligned thin/thick | yes, thin/thick | yes, H/He absorption | yes, AIA-style PSF for Cartesian native EUV images | legacy thin plus native `spherical` thin/thick EUV instrument-grid and dat-resolution |
| SXR image | yes, thin | yes, axis-aligned thin | no | no | no native post-process | legacy thin instrument-grid only |
| White light | spherical legacy | no | no | no | legacy LASCO-like PSF | yes, primary supported path |
| `radio_ff` | limited through EUV image path | yes | yes, thin/thick | yes | yes, independent Gaussian radio beam | no native support |
| `pseudo_current` | limited through EUV image path | yes | yes, thin only | no | no | no native support |

Practical recommendations:

- For cold-material EUV absorption in Cartesian data, use
  `radiation_transfer='thick'`; for non-axis-aligned views, also use
  `ray_method='cart'`.
- For dat-resolution data with non-uniform output spacing, use VTU output. VTI
  is allowed only when the emitted image-plane grid is uniform.
- For radio quick-look products, the current `radio_ff` path is sufficient for
  multi-view brightness-temperature images. Use `radio_beam_fwhm` only when an
  observational beam is needed.
- For spherical data, use the legacy instrument-resolution projection for
  optically thin SXR and white light. For native EUV ray tracing, use
  `ray_method='spherical'` with `radiation_transfer='thin'` or
  `'thick'`; both instrument-resolution and dat-resolution EUV images are
  supported.

## Local demo guards

The radiation-synthesis examples documented here are intended as local demo and
developer guards. They are not wired into the top-level autotest/CI suites by
default, because the full 3D convert cases are heavier than the regular
regression matrix.

## 3D demo case

The Cartesian demo case is `tests/demo4/RadiationSynthesis_3D`. It builds a
small TDm/prominence AMR cube and keeps two thick EUV reference views:

- axis-aligned optically thick AIA 171 image with `tau` and
  `absorption_fraction`, using `convert_171_thick.par`.
- oblique native Cartesian optically thick AIA 171 image with `tau` and
  `absorption_fraction`, using `convert_171_cart_dda_thick_oblique.par`.

Run the 3D Cartesian demo manually with:

```sh
make -C tests/demo4/RadiationSynthesis_3D -f test.make
```

The generated VTI files are intended for visual inspection in ParaView rather
than as bit-for-bit regression baselines.

## Performance baseline

The current performance-sensitive paths are:

- axis-aligned thick transfer, which streams LOS layers in small batches instead
  of allocating a full `Npix * Nlos` column buffer.
- `cart`, which preselects block footprints on the image plane and routes
  per-pixel ray segments to owner MPI ranks before sorting. It emits lightweight
  profile counters for ray tests, ray hits, segment counts, MPI exchange volume,
  and sorting work in the per-conversion logs.

For quick performance comparisons, use the same TDm/prominence demo case and
record wall time for the retained reference conversions:

- `convert_171_thick.par`
- `convert_171_cart_dda_thick_oblique.par`

The most useful next profile points are ray-box intersection counts, collected
segment counts, `MPI_ALLTOALLV` volume, and same-pixel segment sorting time.

## Spherical-coordinate support

The legacy spherical path is the historical instrument-resolution projector.
It initializes a LOS/image-plane basis from `LOS_theta`, `LOS_phi`, and
`image_rotate`, sub-samples spherical cells in `r/theta/phi`, maps each
sub-sample to the image plane, and distributes the contribution through the
instrument PSF. A native ray/mesh-intersection path is also available for EUV
images with `ray_method='spherical'`.

Supported today:

- optically thin EUV instrument-resolution images on spherical grids.
- native spherical thin and thick EUV images with
  `ray_method='spherical'` and `radiation_transfer='thin'` or `'thick'`;
  both instrument-resolution and dat-resolution image grids are supported.
- spherical thick EUV diagnostic maps `tau` and `absorption_fraction`.
- optically thin SXR instrument-resolution images on spherical grids.
- LASCO-like white-light `B` and `pB` images on spherical grids.
- EUV spectra on spherical grids through the legacy instrument-resolution
  spectral path.
- occultation/opaque-inner-region handling through `R_opt_thick` for EUV/SXR
  visibility checks and the white-light occultor parameters for LASCO products.

Not supported yet:

- spherical .dat-resolution SXR images.
- `cart` on spherical data.
- `instrument_postprocess` on spherical data.
- `radio_ff` or `pseudo_current` as validated spherical products.
- `spherical` for SXR, white light, radio, pseudo-current, polar-axis
  crossing domains, or phi-wrapping domains.

The native spherical path constructs ray intersections with spherical grid
faces, sorts path points by LOS distance, and in thick EUV mode reuses the
`j/kappa` transfer layer for strict front-to-back integration. It still rejects
polar-axis crossing, phi-wrapping, and non-EUV products.

The spherical demo case is `tests/demo4/RadiationSynthesisSphericalTDm_3D`. It
builds a lightweight spherical TDm/prominence-style AMR shell and keeps two
native `spherical` thick EUV reference views:

- top-down thick AIA 171, using `convert_171_sph_native_topdown_thick.par`.
- 45-degree top-corner thick AIA 171, using
  `convert_171_sph_native_corner45_top_thick.par`.

Run the spherical demo manually with:

```sh
make -C tests/demo4/RadiationSynthesisSphericalTDm_3D -f test.make
```

The generated VTI files contain `AIA171_thick`, `tau`, and
`absorption_fraction` and are intended for ParaView inspection.

## Remaining known scope limits

- `cart` currently requires Cartesian slab dat-resolution grids.
- Dat-resolution VTI output is restricted to uniform image-plane spacing; use
  VTU for non-uniform image grids.
- `spherical` currently supports only spherical EUV images. In
  dat-resolution mode it outputs AIA intensity, plus optional `tau` and
  `absorption_fraction` in thick mode; Doppler and instrument post-processing
  are not yet defined for this path.
- The old names `cart_dda` and `sph_intersection` are accepted as aliases for
  `cart` and `spherical`.
- `pseudo_current` is optically thin only.
- `radio_ff` is a first thermal free-free implementation with a compact
  Gaunt-factor approximation, brightness-temperature output, and optional
  Gaussian beam post-processing.
