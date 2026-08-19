# Magnetic Topology and QSL Products

This module provides field-line topology products for Cartesian and spherical
MHD data:

- field-line length;
- twist number `Tw`;
- standard `Q` / `logQ` computed by tangent-vector transport;
- perpendicular squashing factor `Qperp`.

The user-facing `q` / `logQ` product is not the former endpoint
finite-difference Q, and it is not the Qperp scalar renamed as Q. Product
geometry determines where Q is sampled. The standard-Q mapping is determined by
the field-line endpoint boundary pair. Cartesian standard `logQ` is the
Cartesian six-face boundary mapping sampled at product points. Spherical
standard `logQ` is radial-boundary mapping sampled at product points, admitting
only `rmin -> rmin`, `rmin -> rmax`, and `rmax -> rmin` endpoint pairs. Qperp
remains an independent perpendicular squashing diagnostic using endpoint
sections perpendicular to the local magnetic field.

## Product Policy

| Workflow | Public products | Notes |
| --- | --- | --- |
| `axis_plane_full_vtu` | length, optional `Tw`, optional `logQ`, optional Qperp in minimal VTI; full diagnostics in VTU | Minimal output is a singleton-dimension VTI coordinate plane. Full VTU diagnostics still use the fixed full product set. RK2 and RK45 support uniform, stretched, and AMR Cartesian grids. |
| `axis_plane_csv` | length, optional `Tw`, optional `logQ`, optional Qperp | CSV files are split by product suffix, but requested science arrays share the seedset trace. RK2 and RK45 support uniform, stretched, and AMR Cartesian grids. |
| `seed_products` | length, optional `Tw`, optional `logQ`, optional Qperp | Cartesian point sampling and spherical point sampling. Spherical standard `logQ` uses the admitted radial-boundary endpoint policy. |
| `arbitrary_plane_products` | length, optional `Tw`, optional `logQ`, optional Qperp | The arbitrary plane is sampling geometry, not the Q mapping surface. |
| `volume_vti` | length, optional `Tw`, optional `logQ`, optional Qperp | VTI only; CSV is intentionally absent. RK2 and RK45 support uniform, stretched, and AMR Cartesian simulation grids. |
| `spherical_surface_products` | length, optional `Tw`, optional `logQ`, optional Qperp | Curved/cut spherical product sampling at rmin, rconst, theta_const, phi_const, or radial_plane points. Standard `logQ` is radial-boundary endpoint mapping Q, not a theta/phi/cut-surface Q. |
| `spherical_cloud_products` | length, optional `Tw`, optional `logQ`, optional Qperp | VTU-only spherical r/theta/phi cloud sampling. Spherical standard `logQ` uses the admitted radial-boundary endpoint policy. |

Former endpoint finite-difference Q, robust stencil modes, endpoint-stencil quality diagnostics, mixed-source Q fields, and Q-specific `*_vis` aliases are not user-facing products.

Spherical surface products, `seed_products`, and `spherical_cloud_products`
support RK2 and `rk45_spherical` for Q-only, Qperp-only, and fused combined
Q+Qperp requests.

## Minimal Visualization Output

Minimal visualization output contains only requested science arrays:

- `length_total`, when `mt_compute_length=.true.` (the default)
- `twist_total`, when requested
- `logQ`, when requested
- `logQperp`, when requested

Minimal VTU/VTI output does not include diagnostic-heavy arrays such as
validity/status flags, endpoint coordinates, boundary faces, split
forward/backward lengths, `N2`, `bfactor`, or connection masks. Full-detail
VTU/VTI and CSV paths may include those diagnostics. If raw `q` is written in
full output, it is the same Scott q0 product represented by `logQ`.

## Output Format Policy

Axis-plane minimal science output is VTI because each coordinate plane is a
structured, axis-aligned sampling lattice: `xy` uses `nx x ny x 1`, `xz` uses
`nx x 1 x nz`, and `yz` uses `1 x ny x nz`. The VTI lattice is the requested
output sampling grid and remains valid when tracing through stretched or AMR
Cartesian simulation grids. Axis-plane full diagnostics remain VTU, and
axis-plane CSV remains a diagnostic/export path.

`volume_vti` remains VTI-only. `seed_products` and `arbitrary_plane_products`
use VTU for minimal visualization because they are respectively a point cloud
and an arbitrarily oriented embedded surface.

## OpenMP Post-Processing

Magnetic-topology conversion is a single-MPI-rank backend. Run it with
`convert_type='magnetic_topology'` and use OpenMP threads inside that rank:

```sh
OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores \
  mpirun -np 1 ./amrvac -convert -if snapshot.dat -i amrvac.par topology.par
```

Do not launch this conversion with multiple MPI ranks. The separate
`convert_type='magnetic_helicity'` backend is the MPI-parallel path.

The tested OpenMP route is the gfortran `ARCH=openmp` build, which enables
`-fopenmp`. The default `ARCH=default` build is serial/non-OpenMP. No Intel
`ifx`/`ifort` OpenMP architecture is added here.

Finite-volume relative magnetic helicity is a separate MPI conversion backend;
see [`magnetic_helicity.md`](magnetic_helicity.md). It is not combined with
the OpenMP topology task in one conversion.

Build the AMRVAC4 magnetic-topology QSL demo with OpenMP, for example:

```sh
cd /home/nanami/codes/amrvac4/tests/demo/AMRVAC4_MagneticTopologyQSL_3D
AMRVAC_DIR=/home/nanami/codes/amrvac4 make ARCH=openmp -j4
```

Run one process with OpenMP threads for local post-processing:

```sh
OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores \
  ./amrvac -i par_showcase/showcase_zstretch_yz_vertical_512_rk45.par
```

The AMRVAC4 demo organizes reproducible TDm/RBSL snapshot setup and QSL
post-processing examples under
`tests/demo/AMRVAC4_MagneticTopologyQSL_3D/`: snapshot-generation par files are
in `par_init/`, showcase examples in `par_showcase/`, timing/efficiency
examples in `par_benchmark/`, and focused tracing diagnostics in
`diagnostics/`.

OpenMP currently accelerates seed-level tracing for `axis_plane_full_vtu` with
both `mt_vtk_detail = 'minimal'` and `mt_vtk_detail = 'full'`. It also covers
the seed-level product loops for `arbitrary_plane_products`, `seed_products`,
and `volume_vti`. VTU/VTI/CSV writing remains serial. Start with
`OMP_NUM_THREADS=4` or `8` and benchmark locally; avoid oversubscription such
as running many MPI ranks while each rank also uses many OpenMP threads. Very
large volume products can still be expensive because snapshot loading, final
product arrays, and VTI writing are not distributed.

## Namelist Runner

Recommended Cartesian showcase settings use RK45 with local cell-fraction
step control:

```fortran
mt_trace_integrator   = 'rk45_cartesian'
mt_step_control       = 'cell_fraction'
mt_step_fraction      = 4.0d0
mt_max_steps          = -1
mt_rk45_rtol          = 1.d-4
mt_rk45_tangent_rtol  = 1.d-3
```

The practical RK2 alternative is:

```fortran
mt_trace_integrator   = 'rk2'
mt_step_control       = 'cell_fraction'
mt_step_fraction      = 2.0d0
mt_max_steps          = -1
```

Fixed-step RK2, for example `mt_dL=0.08d0` with
`mt_max_steps=1600`, is retained in the Cartesian demo only as a
historical/debug comparison. These recommendations are Cartesian product
settings; spherical products have their own demo-local smoke/showcase pars.

Typical axis-plane minimal VTI use:

```fortran
&magnetic_topology_list
  mt_mode = 'axis_plane_full_vtu'
  mt_output_file = 'qsl_axis_xy.vti'
  mt_vtk_detail = 'minimal'
  mt_plane = 'xy'
  mt_xmin = -1.d0
  mt_xmax =  1.d0
  mt_nx = 128
  mt_ymin = -1.d0
  mt_ymax =  1.d0
  mt_ny = 128
  mt_z0 = 0.d0
  mt_trace_integrator = 'rk45_cartesian'
  mt_step_control = 'cell_fraction'
  mt_step_fraction = 4.0d0
  mt_rk45_rtol = 1.d-4
  mt_rk45_tangent_rtol = 1.d-3
  mt_max_steps = -1
  mt_compute_twist = .true.
  mt_compute_q = .true.
  mt_compute_qperp = .true.
/
```

With `mt_vtk_detail='minimal'`, `axis_plane_full_vtu` writes VTI and respects
the compute flags: `length_total` is controlled by `mt_compute_length` (default
`.true.`), and `twist_total`, `logQ`, and `logQperp` are written only when
requested. With `mt_vtk_detail='full'`, the same task writes VTU full
diagnostics and keeps the fixed full product set.

Typical axis-plane CSV use:

```fortran
&magnetic_topology_list
  mt_mode = 'axis_plane_csv'
  mt_output_prefix = 'qsl_axis_csv'
  mt_plane = 'xy'
  mt_xmin = -1.d0
  mt_xmax =  1.d0
  mt_nx = 128
  mt_ymin = -1.d0
  mt_ymax =  1.d0
  mt_ny = 128
  mt_z0 = 0.d0
  mt_dL = 0.05d0
  mt_max_steps = 1000
  mt_compute_twist = .true.
  mt_compute_q = .true.
  mt_compute_qperp = .true.
/
```

The CSV runner writes length by default, plus twist, standard `logQ`, and Qperp
when enabled. Standard `logQ` is written to the `_q.csv` product file; Qperp is
written to `_qperp_method2.csv`. When length, twist, standard `logQ`, and
Qperp are all requested, the CSV runner uses the same shared Qperp-capable
tangent trace as the generic Cartesian seedset products.

## Direct Wrappers

Coordinate-plane wrappers:

- `mt_qsl_plane_vtu_xy/xz/yz`: direct/legacy axis-plane VTU writers with Scott q0 `Q` and optional Qperp. The namelist minimal visualization path uses VTI.
- `mt_topology_plane_xy/xz/yz`: shared length/twist/mapping CSV trace.
- `mt_q_plane_xy/xz/yz`: standard Cartesian `logQ` CSV.
- `mt_qperp_plane_xy/xz/yz`: Qperp-only CSV.

Other product wrappers:

- `mt_fieldline_products_seeds`: seed-set length, optional twist, optional Cartesian `logQ`, optional Qperp.
- `mt_fieldline_products_plane_arbitrary`: arbitrary orthonormal-plane length, optional twist, optional Cartesian `logQ`, optional Qperp.
- `mt_fieldline_products_volume_vti`: Cartesian volume length, optional twist, optional Cartesian `logQ`, optional Qperp.

## Interpretation

`logQ` and `logQperp` are related but distinct diagnostics. Cartesian `logQ`
is a Scott q0 quantity using endpoint boundary-face normals and
tangent-vector transport on the Cartesian six-face boundary mapping. Spherical
`logQ` uses tangent-vector transport on admitted radial-boundary endpoint
mappings: `rmin -> rmin`, `rmin -> rmax`, and `rmax -> rmin`. It does not
enable `rmax -> rmax`, theta/phi boundary Q, mixed radial-theta/radial-phi Q,
or theta-phi Q. Product geometry, including seed clouds, coordinate planes,
arbitrary planes, and volumes, is only the sampling strategy. `logQperp` is
the perpendicular squashing factor computed from endpoint sections
perpendicular to the local magnetic field. Minimal visualization output
intentionally omits connection masks and status fields. Full-detail output and
CSV diagnostics may include connection masks, endpoint statuses, termination
faces, and validity metadata.

Unsupported or intentionally absent products:

- endpoint finite-difference endpoint-FD Q;
- endpoint-stencil diagnostics;
- mixed-source Q fields;
- Q-specific visual aliases;
- endpoint finite-difference substitutes for standard `logQ`;
- spherical standard `logQ` outside the admitted radial-boundary endpoint
  policy;
- CSV output for `volume_vti`.

MPI distributed tracing and large AMR scaling performance are outside the
axis-plane QSL 1.0 scope unless validated and documented separately.
