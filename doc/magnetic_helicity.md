# Finite-volume magnetic helicity

`mod_magnetic_helicity` is the MPI post-processing backend for relative
magnetic helicity in a finite Cartesian volume.  It is intentionally separate
from `mod_magnetic_topology`: the latter traces field lines with OpenMP and
must run on one MPI rank, whereas this module distributes planes with MPI and
does not call the field-line tracer.

## Quantities and physical meaning

For a magnetic field \(\mathbf B\) in a finite volume, the reference field
\(\mathbf B_p\) is the unique current-free field with the same normal
component on every physical face.  With vector potentials satisfying the
same tangential-boundary convention, the gauge-invariant relative helicity is

\[
 H_m=\int_V(\mathbf A+\mathbf A_p)\cdot(\mathbf B-\mathbf B_p)\,dV .
\]

The implementation also reports the gauge-invariant current-carrying and
volume-threading decomposition attributed to Berger (2003) and written in a
finite-volume form by Valori, Démoulin, and Pariat (2012):

\[
 H_J=\int_V(\mathbf A-\mathbf A_p)\cdot(\mathbf B-\mathbf B_p)\,dV,
 \qquad
 H_{PJ}=2\int_V\mathbf A_p\cdot(\mathbf B-\mathbf B_p)\,dV,
\]

so that `Hm = HJ + HPJ` is checked independently.  The factor 2 is part of
the code's `HPJ` definition, not an additional post-processing multiplier.
`HJ` measures the current-carrying (non-potential) part of the field; `HPJ`
is the volume-threading/mutual term between the reference field and the
current-carrying field.  Neither is the field-line twist `Tw` from the topology
backend.  The optional
`abs_HJ_over_abs_Hm` is omitted (`ratio_is_valid=false`) when the total
helicity is at the round-off scale.

The energy diagnostics are

\[
 E=\frac12\int_V|\mathbf B|^2dV,\qquad
 E_p=\frac12\int_V|\mathbf B_p|^2dV,\qquad E_{free}=E-E_p.
\]

In a normalized CGS run, helicities are written in \(\mathrm{Mx}^2\) and
energies in erg; with `SI_unit=.true.` they are written in \(\mathrm{Wb}^2\)
and joule.  Code-unit columns are always retained.

## Reference potential field contract

The reference service `mod_magnetic_reference_fv` solves

\[
 \nabla^2\phi=0,\qquad
 \partial_n\phi=\mathbf B\cdot\hat{\mathbf n}
 \quad\text{on all six faces},\qquad
 \mathbf B_p=\nabla\phi .
\]

The boundary value is the *total* normal field.  For cell-centred storage it
is obtained by the adjacent interior/ghost interpolation; for staggered
storage the face value `ws` is used.  `B0field` is added in both paths.  A
pure-Neumann problem is compatible only when the signed sum of the six face
fluxes is zero.  The code reports the relative imbalance and stops when it
exceeds `mh_max_flux_imbalance`; it never silently rebalances the magnetogram.
The additive constant in \(\phi\) is removed by the multigrid mean
constraint, which does not change \(\nabla\phi\).  The first reference solve
uses one FMG pass and subsequent residual reductions use standalone V-cycles.
The reusable reference-service default remains `1.d-8`; the helicity user
entry point defaults to `mh_mg_tolerance=1.d-4`.

## DeVore–GV gauge and plane pipeline

The vector potentials use the DeVore–GV construction.  For a selected axis
\(s\) (`mh_gauge_axis=1`, `2`, or `3`; the default is `3`, the z direction),

\[
 A_s=A_{p,s}=0,\qquad
 \mathbf A=\mathbf b+\hat{\mathbf e}_s\times
       \int_s^{s_{max}}\mathbf B\,ds',\qquad
 \mathbf A_p=\mathbf b_p+\hat{\mathbf e}_s\times
       \int_s^{s_{max}}\mathbf B_p\,ds'.
\]

The top-face transverse vector \(\mathbf b=\mathbf b_p\) is built from a
one-dimensional integral of the top-face normal field.  Along the selected
axis, a second-order recurrence matches the centred finite-volume curl, so
the interior checks `curl_A_error` and `curl_Ap_error` directly test the
discrete reconstruction.

The analysis first remaps an AMR snapshot through `level_io` to one fixed
uniform level.  It then assembles one plane at a time.  Plane owners are
distributed round-robin across MPI ranks; neighbouring owners exchange only
the accumulated integration state.  The full three-dimensional magnetic
field is not replicated on each rank.  Work arrays scale as
\(O(N_\perp^2)\) per rank in the scalar path, in addition to the distributed
reference multigrid data.  Debug VTI output is streamed separately and is
off by default.

## Conversion and namelist

Use the independent conversion type:

```fortran
&filelist
  convert=.true.
  convert_type='magnetic_helicity'
  level_io=1
/

&magnetic_helicity_list
  mh_output_file       = ''
  mh_gauge_axis        = 3
  mh_mg_tolerance      = 1.d-4
  mh_mg_max_cycles     = 50
  mh_max_flux_imbalance= 1.d-6
  mh_write_debug_vti   = .false.
  mh_write_mg_timing   = .false.
  mh_debug_vti_file    = ''
/
```

An empty `mh_output_file` uses `<base_filename>_helicity.csv`.  Set
`mh_write_debug_vti=.true.` and optionally `mh_debug_vti_file` to stream
`Bp`, `A`, `Ap`, `h_m`, `h_J`, and `h_PJ` for numerical inspection.  The
optional `mh_write_mg_timing` prints the reference-solve multigrid timer
table without changing the CSV schema.

The CSV header contains snapshot and time fields, code and physical-unit
`Hm`, `HJ`, `HPJ`, `abs_HJ_over_abs_Hm`, `E`, `Ep`, `Efree`, net-flux and
boundary-normal errors, `epsilon_div_B`, `curl_A_error`, `curl_Ap_error`,
`decomposition_error`, `mg_residual`, `mg_cycles`, and `ratio_is_valid`.

Run with MPI and keep OpenMP/BLAS at one thread when benchmarking:

```bash
OMP_NUM_THREADS=1 OMP_DYNAMIC=false \
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun --bind-to none --map-by slot -np 4 ./amrvac -convert \
  -if output/tdm_0000.dat -i amrvac_generate.par amrvac_helicity.par
```

The directly usable two-step TDm/RBSL example is
[`tests/demo4/MagneticHelicity_Cart/README.md`](../tests/demo4/MagneticHelicity_Cart/README.md).
The topology documentation remains at
[`magnetic_topology_qsl.md`](magnetic_topology_qsl.md); it documents the
single-rank OpenMP `convert_type='magnetic_topology'` path only.

## Verified references

1. M. A. Berger & G. B. Field, “The topological properties of magnetic
   helicity,” *Journal of Fluid Mechanics* **147**, 133–148 (1984),
   doi:[10.1017/S0022112084002019](https://doi.org/10.1017/S0022112084002019).
2. J. M. Finn & T. M. Antonsen, Jr., “Magnetic helicity: what is it and what is it
   good for?”, *Comments on Plasma Physics and Controlled Fusion* **9**(3),
   111–126 (1985).  No DOI was returned by the Crossref/OpenAlex records
   checked for this historical article, so none is asserted here.  The
   Pascal–Francis landing record contains the indexing variant “FININ”; the
   author-index record and Finn’s other plasma publications identify the
   author as John M. Finn.
3. M. A. Berger, “Topological quantities in magnetohydrodynamics,” in
   A. Ferriz-Mas & M. Núñez (eds.), *Advances in Nonlinear Dynamos*,
   *The Fluid Mechanics of Astrophysics and Geophysics*, vol. 9, Taylor &
   Francis, London/New York, 345–374 (2003).  This is the original 2003
   print-edition citation; the later electronic reissue is indexed by
   Crossref as doi:[10.1201/9780203493137-10](https://doi.org/10.1201/9780203493137-10)
   and is not a change to the original publication year.
4. C. R. DeVore, “Magnetic Helicity Generation by Solar Differential Rotation,”
   *The Astrophysical Journal* **539**, 944–953 (2000),
   doi:[10.1086/309274](https://doi.org/10.1086/309274).
5. G. Valori, P. Démoulin & E. Pariat, “Comparing Values of the Relative
   Magnetic Helicity in Finite Volumes,” *Solar Physics* **278**, 347–366
   (2012), doi:[10.1007/s11207-012-9951-6](https://doi.org/10.1007/s11207-012-9951-6).
6. E. Pariat, J. E. Leake, G. Valori, M. G. Linton, F. P. Zuccarello &
   K. Dalmasse, “Relative magnetic helicity as a diagnostic of solar
   eruptivity,” *Astronomy & Astrophysics* **601**, A125 (2017),
   doi:[10.1051/0004-6361/201630043](https://doi.org/10.1051/0004-6361/201630043).
7. G. Valori, E. Pariat, S. Anfinogentov, F. Chen, M. K. Georgoulis, Y. Guo,
   Y. Liu, K. Moraitis, J. K. Thalmann & S. Yang, “Magnetic Helicity
   Estimations in Models and Observations of the Solar Magnetic Field. Part I:
   Finite Volume Methods,” *Space Science Reviews* **201**, 147–200 (2016),
   doi:[10.1007/s11214-016-0299-3](https://doi.org/10.1007/s11214-016-0299-3).

The relative-helicity definition and its boundary condition are anchored in
Berger–Field (1984) and Finn–Antonsen (1985).  The `HJ`/`HPJ` split is the
Berger (2003) current-carrying/volume-threading decomposition, while Valori,
Démoulin, and Pariat (2012) provide the finite-volume derivation and the
implementation/benchmark comparison used here.  Fu/Yu (2023) is not used as a
foundational definition citation and is not needed in this module
documentation.

## Limitations

Version 1 accepts only three-dimensional, uniformly spaced, non-stretched
Cartesian analysis meshes.  Native AMR, spherical coordinates, and online
time-step diagnostics are outside this conversion path; AMR data must be
remapped with `level_io`.  Incompatible six-face net flux is a hard error.
Debug VTI is a numerical audit product, not the primary science output.
