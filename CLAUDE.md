# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

AGILE is a GPU-enabled fork of MPI-AMRVAC, a Fortran finite-volume code for
solving hyperbolic PDEs (HD, MHD, FFHD, SRHD) with adaptive mesh refinement.
Master supports 3D Cartesian, 3D spherical `(r, theta, phi)` and 3D
cylindrical `(r, z, phi)` grids, each of the latter two also in a
logarithmically stretched radial variant (`logSpherical`, `logCylindrical`);
see "Coordinate systems" below for the limits of the curvilinear support.
Domain decomposition uses MPI; GPU offload supports two selectable backends,
OpenACC and OpenMP target offload. Every offload site in the hot loops is
written once as a call to a macro from `src/mod_gpu_directives.fpp` (e.g.
`${GPU_PARALLEL_LOOP(...)}$`), which fypp expands to the matching
`!$acc ...` or `!$omp target ...` directive depending on which backend is
selected at build time (see `OPENACC=1`/`OPENMP=1` below) — don't write raw
`!$acc`/`!$omp` directives directly in hot-loop code, add a macro call
instead.

Source is written as `.fpp` files (Fortran + `fypp` preprocessor directives,
e.g. `#:if PHYS == 'hd'`) which get preprocessed into `.f90` before
compilation. Unlike upstream MPI-AMRVAC, this fork does **not** use the old
VACPP/LASY dimension-independent `.t` preprocessor (`doc/vacpp.md` and much of
`doc/*.md` describe upstream MPI-AMRVAC and are only partially applicable
here — e.g. `setup.pl` and `arch/*.arch` files referenced there don't exist in
this repo). Loop bounds are written out explicitly per-dimension (e.g.
`ixGmin1,ixGmin2,ixGmin3`) directly in the `.fpp` source.

## Environment setup

```bash
pip install uv                       # or curl -LsSf https://astral.sh/uv/install.sh | sh
export AGILE_DIR=/path/to/this/repo  # must point at the repo root
cd $AGILE_DIR && uv sync
source $AGILE_DIR/.venv/bin/activate
```

`AGILE_DIR` and an activated venv (providing `fypp`, `fortdepend`, `f90nml`)
are required for every build. Requires Python >=3.13.

## Building

Each test/run case is its own directory with a `Makefile` and an `agile.par`
parameter file. Build from inside such a directory:

```bash
cd tests/hd/KH3D
make arch=gnu              # compiles and links ./agile executable
make arch=gnu OPENACC=1    # enable GPU build with the OpenACC backend
make arch=gnu OPENMP=1     # enable GPU build with the OpenMP target-offload backend
make arch=gnu DEBUG=1      # debug flags (-g -O0 -fcheck=all, etc.)
make clean                 # remove this case's build products
```

- `arch` selects a file from `arch/*.mk` (`gnu`, `ifx`, `nvidia`, `cray`,
  `llvm`) which sets the compiler/linker and flags. Default is `gnu`
  (`mpif90`/gfortran).
- `OPENACC` and `OPENMP` are mutually exclusive GPU offload backend switches
  (`make` errors out if both are given); neither flag gives a plain CPU
  build, where every `GPU_*` directive macro expands to nothing.
- `NOGPUDIRECT=1` (only meaningful with `OPENACC` or `OPENMP`) stages MPI
  buffers through the host instead of handing device pointers to MPI. **Use
  it when all ranks share a single GPU.** Without it, a workstation with one
  GTX 1660 Ti (nvhpc 26.9, HPC-X) runs `tests/hd/KH3D` correctly on one
  rank but loses mass deterministically on two or more (`mean(rho)` 5.42 at
  2 ranks, 5.37 at 4, against 5.5), whichever MPI build or UCX transport is
  used; with `NOGPUDIRECT=1` it passes at 1, 2 and 4 ranks. The cause is not
  understood, so treat this as a workaround, not a known limit. Multi-GPU
  machines (e.g. the GPU CI) keep GPU-direct.
- Physics-relevant compile-time options (which physics module, the coordinate
  system, tracers, gravity, cooling, thermal conduction, etc.) are declared in
  `agile.par` and turned into `config.mk`/fypp defines by
  `make/config_reader.py`, validated
  against `make/config_schema.toml`. Editing `agile.par` and rerunning `make`
  regenerates `config.mk` automatically — don't hand-edit `config.mk`.
- Builds are cached under `build/<arch>-<hash>/`, keyed by a hash of the
  enabled compile-time flags; `build/latest` symlinks to the most recent one.
  `make clean-all` wipes the whole `build/` cache.
- The actual make logic lives in `make/*.mk`, included in numeric order by
  `make/include-all.mk` (00-prelude, 05-find-agile, 10-select-arch,
  15-read-config, 20-build-dir, 25-store-flags, 30-fypp, 35-dependencies,
  40-compile, 50-link, 90-clean).

## Running tests

Tests live under `tests/{hd,mhd,ffhd,srhd,srmhd}/<case>/` and each case has a
`test.make` including `tests/test_rules.make`, which builds the case fresh
and runs it (via `mpirun`/`srun`), then numerically diffs the log output
against `tests/<physics>/<case>/correct_output/` using
`tools/fortran/compare_logs` (tolerances `1.0e-5` relative / `1.0e-8`
absolute).

```bash
cd tests
make hd            # run the basic hd/* test cases
make mhd           # run the basic mhd/* test cases
make ffhd          # run the basic ffhd/* test cases
make srhd          # run the basic srhd/* test cases
make srmhd         # run the basic srmhd/* test cases
make SUITE=all hd  # basic and extended hd/* cases
make SUITE=all list # print the directories a SUITE selects
make -j4 hd         # run in parallel
make -s hd          # suppress "Entering/leaving directory" noise
make clean          # clean all test build artifacts
```

The test directories are split into two suites in `tests/Makefile`: each
physics module has a `*_BASIC` and a `*_EXTENDED` list. `SUITE` (a make
variable or environment variable, like `arch`) selects `basic` (the default),
`extended` or `all`. The split is per directory, because a directory is one
build and builds dominate the cost. The basic suite keeps one Cartesian case
per physics plus the curvilinear cases whose checks nothing else provides
(`hd/spherical` for refluxing, `hd/log_cylindrical_pole` for `log_r0`, ...).
Naming a single directory (`make hd/spherical`) runs it regardless of `SUITE`.

In CI (`.github/workflows/tests.yml` and `gpu.yml`) pushes and PRs run the
basic suite. The full suite (`all`) runs weekly on master (Sunday 02:00
UTC), on a PR carrying the `full-tests` label (adding the label triggers the
run), and from a manual run, whose `suite` input chooses between
`basic`, `extended` and `all` (`gh workflow run tests.yml -f suite=all --ref <branch>`).

To run a single case's test directly:

```bash
cd tests/hd/KH3D
make -f test.make arch=gnu
```

This is also what CI does, via `tests/test_runner.sh <target>` (used as
`arch=gnu ./test_runner.sh hd`), which fails if the build fails or if
`**FAILED**` appears in the output log.

Add a new test case by copying an existing case directory (`mod_usr.fpp`,
`agile.par`/`test.par`, `test.make`, `correct_output/`) and registering it in
`tests/Makefile` under the appropriate `*_BASIC` or `*_EXTENDED` variable (extended unless
it covers something no basic case does).

## Coordinate systems

The coordinate system is a **compile-time** choice, because the finite-volume
kernels are fypp-specialized for it and the Cartesian kernel must stay free of
curvilinear bookkeeping. Set `geometry` in `&meshlist` of `agile.par`:

```fortran
 &meshlist
   geometry = 'spherical'   ! 'cylindrical', 'logSpherical',
                            ! 'logCylindrical', or 'Cartesian' (the default)
 /
```

That becomes the `GEOM` fypp define via `make/config_schema.toml`, and it is
the only place the coordinate system needs stating. `main()` in
`src/agile.fpp` calls `set_coordinate_system_from_config()` before `usr_init`,
so a case does **not** have to call `set_coordinate_system` itself. A case may
still do so — the call is idempotent — but the name then has to agree with
`GEOM`, and `set_coordinate_system` calls `mpistop` on a mismatch rather than
silently running the wrong kernel.

Deriving it centrally matters because `set_coordinate_system` is the sole
writer of `coordinate`, `ndir`, `r_`, `phi_` and `z_`. When it was left to
each case, forgetting the call left `coordinate` at its `-1` default, and the
failure was silent rather than loud: `read_par_files` skips the conversion of
the angular domain bounds, and the `select case (coordinate)` in
`get_surface_area` and in the cell-volume fill match nothing, so `surfaceC`
and `dvolume` are never written and the run completes with garbage.
`tests/hd/spherical_uniform_flow` omits the call and
`tests/hd/spherical_blast` keeps it, so both paths stay covered.

Three things to know when writing a spherical or cylindrical case:

- **Angular coordinates in `agile.par` are given in units of 2*pi** (an
  MPI-AMRVAC convention, applied in `src/io/mod_input_output.fpp`). So for
  spherical, `xprobmin2 = 0.125d0`, `xprobmax2 = 0.375d0` means `theta` from
  `pi/4` to `3*pi/4`. For cylindrical, `set_coordinate_system` assigns
  `r_ = 1`, `z_ = 2`, `phi_ = 3` for the 3D case (`"cylindrical"` /
  `"cylindrical_3D"`), so only `xprobmin3`/`xprobmax3` get the 2*pi
  conversion; `xprobmin2`/`xprobmax2` (`z`) is a plain length like `r`.
- `ndim` is a compile-time `parameter` of 3, so spherical or cylindrical
  means the full 3D `(r, theta, phi)` or `(r, z, phi)` system with
  `ndir = 3`.
- **`ps(igrid)%x` is the volume barycentre of the cell, not the midpoint of
  its faces**, so `x +/- dx/2` is not a face. See "Cell positions are volume
  barycentres" below.

How it works: in a curvilinear build the finite-volume update in
`src/mod_finite_volume.fpp` divides face fluxes by the cell volume weighted by
the face areas (`bgeo%surfaceC`, `bgeo%dvolume`, filled per block by
`fill_geometry_device` in `src/amr/mod_amr_solution_node.fpp`), instead
of by the uniform `rnode` cell spacing. The leftover curvature terms of the
momentum equations are added by `addsource_geometry`, a fypp macro each physics
template defines (see `src/hd/mod_hd_templates.fpp`);
`src/physics/mod_physics_dummies.fpp` provides a version that `mpistop`s, so a
physics module without one fails loudly. `setdt` likewise switches from `rnode`
spacing to the physical cell sizes `bgeo%ds`.

A block's geometry lives in `bgeo` (and `bgeoc` for the one-level-coarser
representatives), a `geo_t` declared in `src/mod_physicaldata.fpp` that holds
`x`, `ds`, `dvolume`, `surfaceC` and `dx` for *all* blocks with the grid index
last, exactly as `block_grid_t` holds the solution in
`bg%w`. `initialize_vars` allocates them once for `1:max_blocks`, so GPU
kernels index them directly (`bgeo%x(ix1,ix2,ix3,1:ndim,igrid)`) rather than
chasing a per-block pointer, and `alloc_node` only has to push the one block's
slice it just filled. The matching `state` components (`ps(igrid)%x`, `%dx`,
...) are bounds-remapped pointer views into those arrays (see
`point_at_geometry` in `src/amr/mod_amr_solution_node.fpp`) — remapping rather
than plain pointer assignment, because `surfaceC` starts at index 0 and
`dvolume`/`ds` do too when `nghostcells` is odd. Keeping the views is what lets
the large body of host code that uses `ps(igrid)%x` and friends, including
every case's `mod_usr.fpp`, work unchanged.

A block's metrics are built on the device, by `fill_geometry_device` in
`src/amr/mod_amr_solution_node.fpp`, which `alloc_node` calls once per block.
There is a single path for all three coordinate systems: a uniform block is an
analytic function of nothing but its corner and spacing in `rnode`, so every
member of `geo_t` is derived on the GPU rather than computed on the CPU and
shipped across the bus on each (re)allocation of a node — which, with AMR, is
most regrids. The per-cell expressions are fypp-specialized on `GEOM`, so a
build only carries the branch it needs.

Because the device is the producer, *all* of `geo_t` is device-resident in
every build, including members no kernel reads (`dx`, and in a Cartesian build
the metrics whose values its kernels take from `rnode` instead). They have to
live where they are written. Only the positions are returned to the host by
`alloc_node` itself, since `ps(igrid)%x` is read immediately for the initial
condition, the user hooks and the output.

That residency is not free: every member is a block-sized array for all
`max_blocks` blocks a rank could ever own, allocated up front whether or not
the blocks exist. On a 12 GB card, `max_blocks = 4096` with `block_nx = 16`
puts `bgeo` plus `bgeoc` at over 4 GB before a single block is allocated, and
`initialize_vars` fails with `CUDA_ERROR_OUT_OF_MEMORY` at its `enter data`.
That is a limit of the card, not a bug — the same case runs at
`max_blocks = 2048` — but it is why a member that nothing reads does not stay
in `geo_t` on the argument that it is cheap.

Two members were dropped for exactly that reason, and neither the fine nor the
coarse representative carries them any more:

- `dsC`, the cell-face lengths. Its only reader was
  `b_from_vector_potentialA` in `src/mod_constrained_transport.fpp`, which
  needs `stagger_grid` (unsupported here) and in fact has no callers at all.
- `surface`, the cell-centre face areas. Its only reader was the
  `Stokesbased` branch of `curlvector` in `src/mod_geometry.fpp`, reachable
  only via `typecurl = 'Stokesbased'` in `&methodlist`.

The matching `state` components `s%dsC` and `s%surface` are still declared but
are never associated, so that the two routines above still compile. Both now
call `mpistop` before they would dereference them, naming what was removed —
re-enabling either means giving it back a source for its metric, not just
deleting the guard.

Two consequences of building the metrics analytically:

- **General grid stretching is gone.** A `stretch_dim` that made the cell
  spacing vary within a block could not be expressed by the analytic
  derivation, so upstream MPI-AMRVAC's stretched-grid machinery was removed
  outright: the `stretch_dim`, `stretch_uncentered`, `qstretch_baselevel` and
  `nstretchedblocks_baselevel` keys are no longer read from `&meshlist`, the
  `read_par_files` block that filled the per-level `qstretch`/`dxfirst`/`dxmid`
  tables is deleted, and the `initialize_vars` guard that used to `mpistop` on
  a stretched dimension is gone with it (nothing can request one any more). The
  `qstretch`/`dxfirst`/`stretched_dim`/`stretch_type` globals still *exist* in
  `mod_global_parameters`, pinned to the unstretched state, only because a few
  paths this fork does not exercise (`get_igslice` in `src/io/mod_slice.fpp`,
  the particle-index lookup in `src/particle/mod_particle_base.fpp`, the
  line-of-sight sub-grids in `src/physics/mod_thermal_emission.fpp`) still
  branch on them. The one stretching that *is* supported, a logarithmic
  radius, never needed any of it — it is a uniform mesh in `ln(r + log_r0)`
  and therefore still an analytic function of the block corner and a constant
  spacing (see "Logarithmically stretched radius" below). The ~1300-line host
  path that used to handle general stretching (`fill_geometry_host`), and
  `get_surface_area` in `src/mod_geometry.fpp`, are also gone.
- `Cartesian_expansion` was already unreachable here — `set_coordinate_system`
  `mpistop`s on it unless `ndim == 1`, and `ndim` is a compile-time 3 — so the
  `usr_set_surface` hook in `src/mod_usr_methods.fpp` is likewise now unused.

### Cell positions are volume barycentres

`bgeo%x` — and therefore `ps(igrid)%x`, which every case's `mod_usr.fpp` reads
— holds the **volume barycentre** of the cell, not the midpoint of its faces.
The solution is stored as a cell average, and a cell average equals the point
value at the barycentre to second order, because the linear term of the Taylor
expansion integrates to zero only about that point. So the barycentre is where
the state, the initial condition and any position-dependent source term belong.

Only the directions with a non-uniform volume weight are affected. `phi` in
both systems, and `z` in cylindrical, carry uniform weight, so their
barycentres *are* their midpoints and nothing changed for them. The two that
differ are, for faces at `r_L`, `r_R` and `theta_L`, `theta_R`:

    spherical    r     = 3/4 * (r_R^4 - r_L^4)/(r_R^3 - r_L^3)
    cylindrical  r     = 2/3 * (r_R^3 - r_L^3)/(r_R^2 - r_L^2)
    spherical    theta = theta_c + cot(theta_c) * (sin(u) - u cos(u))/sin(u),
                         u = dtheta/2

`fill_geometry_device` evaluates all three in forms from which the common
factor has been cancelled analytically, because the offset from the midpoint is
`O(h^2)` and computing it as the difference of two `O(1)` quantities would
throw away most of the mantissa on a fine grid. The offset is about 1% of a
cell width in the radial direction of the existing test grids, and grows to
`dtheta/6` — a sixth of a cell, pushed away from the axis — for the first
polar cell off `theta = 0`.

**The consequence to remember is that `x +/- dx/2` is not a cell face.** Three
places in the tree needed that and were changed to rebuild the face from
`rnode` instead:

- the face coordinates handed to the Riemann solver in
  `src/mod_finite_volume.fpp`. All three directions are built the same way
  there, from the block corner and the cell index — that one expression is
  correct in every geometry, and reduces to the old `x -/+ dr/2` wherever `x`
  is the midpoint, so it replaced the per-geometry special cases rather than
  adding to them;
- the VTU corner grid in `calc_x`, `src/io/mod_calculate_xw.fpp`;
- `fill_geometry_device` itself, which is now written face-first: the
  `RADIAL_CELL` macro produces the two faces, the midpoint `rc` and the extent
  `ds1`, and every area and volume is then written in terms of `rc` and `ds1`.
  That is what lets those expressions stay algebraically *exact* — `rc*ds1` is
  precisely `(r_R^2 - r_L^2)/2`, and `(rc^2 + ds1^2/12)*ds1` precisely
  `(r_R^3 - r_L^3)/3`, for any face pair whatsoever — so the same expressions
  serve the uniform and the stretched grids unchanged. In a uniform build
  `rc` and `ds1` are taken straight from the logical centre and spacing rather
  than rebuilt from the faces, so those builds' metrics are bit-for-bit what
  they were before, and any change in their results is attributable to the
  barycentre and the source terms alone.

The curvature source terms take no position at all any more; see the
`addsource_geometry` bullet under "Current limits" below for why that matters
at the axis. `bgeo%dx` remains "the extent of the cell in the same units as the
matching component of `x`", so its radial entry is a physical width.

**`bgeo%ds` deliberately stays on the cell midpoint**, not the barycentre. It
is a geometric extent that `setdt` turns into the CFL length, so the argument
that puts the *state* at the barycentre says nothing about it, and following
the barycentre there would quietly relax the timestep. That is not academic at
the polar axis: the barycentre sits `dtheta/6` further out, so `sin(theta_bar)`
runs about a third above `sin(theta_c)` in the first cell — exactly the cell
whose vanishing `ds(3) = r sin(theta) dphi` sets `dt` for the whole run.
Measured on `tests/hd/spherical_pole/uflow.par`: following the barycentre in
`ds` takes the run from 24 timesteps to 20, a 20% larger `dt`, for no reason
other than a mislabelled cell size. It happens not to move that case's volume
averages — they are dominated by spatial error — so this is a silent change,
which is exactly why it is worth pinning down here.

One place was deliberately *not* fixed: the face averages `half*(q(ix)+q(jx))`
in `divvector` and `curlvector` (`src/mod_geometry.fpp`) assume the face lies
midway between two cell centres, which the volume barycentre is not. Both
routines still carry their upstream `stretched_dim .and. stretch_uncentered`
branch, now permanently dead (`stretched_dim` is pinned `.false.`); the
uncentered extrapolation it does — from `x(ix)` by `dx(ix)/2` along the actual
centre separation — is also the fix for the barycentre offset. Neither routine
is on a path this fork currently exercises, so this is recorded rather than
repaired.

### Logarithmically stretched radius

`geometry = 'logSpherical'` and `geometry = 'logCylindrical'` are the same two
coordinate systems with a radial mesh that is uniform in `ln(r + r0)` rather
than in `r`. `r0` is the `log_r0` entry of `&meshlist` and defaults to zero,
which is the plain `ln(r)`: the cell width then grows in proportion to `r` and
the relative resolution `dr/r` is constant. That is the one stretching
astrophysics usually wants, and it is the only one supported.

A positive `log_r0` moves the singularity of the map from `r = 0` out to
`r = -r0`, so the mesh is *uniform*, with cell width `r0*ds`, for `r << r0` and
logarithmic for `r >> r0`. That is what lets a stretched radial grid reach
`r = 0` at all — see "Reaching the axis with `log_r0`" below.

**It needed almost none of upstream MPI-AMRVAC's stretched-grid machinery
(since removed entirely, see "General grid stretching is gone" above), and
that is the whole design.** Take the *logical* coordinate to be `xi = ln(r)`:
the mesh is then uniform in `xi`, so it is still exactly what
`fill_geometry_device` assumes — an analytic function of the block's corner and
a constant spacing. The tree, the AMR level ladder, `rnode`, block indexing,
prolongation, coarsening and the ghost exchange all keep working unchanged in
`xi`; only the step from logical faces to physical ones takes an `exp`.
Upstream's default `q = (xprobmax1/xprobmin1)^(1/domain_nx1)` stretching *is*
this grid, and its `q_l = sqrt(q_{l-1})` refinement ladder is just
`dxi -> dxi/2`, so refinement and log stretching commute for free rather than
needing per-level `qstretch`/`dxfirst` tables.

Two implementation points follow from that framing:

- **`GEOM` keeps its three values.** `geometry = 'logSpherical'` emits
  `GEOM='spherical'` *plus* a separate `LOG_RADIUS` flag, via a new `emit`
  mapping in `make/config_schema.toml` and `make/config_reader.py`. So every
  existing `#:if GEOM == 'spherical'` branch — the four physics modules'
  curvature terms, the geometry kernel, the finite-volume update, `setdt` —
  stays correct untouched, and `coordinate` stays `spherical` at runtime so
  every `select case (coordinate)` keeps matching. A log grid *is* spherical;
  only its radial spacing differs. The build hash still distinguishes the two,
  and `geometry_name` keeps the log name so a log snapshot refuses to restart
  into a uniformly spaced build.
- **`xprobmin1`/`xprobmax1` are physical radii in the par file** and are
  converted to the logical coordinate by `read_par_files`, next to where the
  angular bounds are already converted from units of `2*pi`. `xprobmin1 <= 0`
  is rejected there unless `log_r0 > 0`, which is what puts `r = 0` on the
  grid. **A case's `mod_usr.fpp` that reads `xprobmin1`/`xprobmax1` for a
  physical radius therefore has to map them back** — with `r_of_s` from
  `mod_geometry`, not by hand, since the map is only an exponential when
  `log_r0` is zero. This is the one thing about these geometries that reliably
  catches people out. The ten `log_spherical`/`log_cylindrical` directories
  predate `log_r0` and still write `exp(xprobmin1)` inline; that is correct
  only because they run at `log_r0 = 0`, and copying one of them as the
  starting point for an offset grid is exactly how the trap gets sprung.
  `tests/hd/log_cylindrical_pole` uses `r_of_s`.

  `r_of_s` and its inverse `s_of_r` (both `pure elemental`, in
  `src/mod_geometry.fpp`) are the single host-side definition of the map, and
  every host site that converts uses them: `read_par_files`, the `slicecoord`
  and `xslice` bounds checks, `get_igslice` and `calc_x`. Device code does not
  call them — the geometry kernel writes the map out inline as a fypp macro
  (an `!$acc routine` there is what `-Minline` previously miscompiled), and
  `mod_finite_volume` uses the reduced form below.

  `slicecoord` in `&savelist` is the exception that proves the rule: it is
  **not** converted, because it is a physical coordinate everywhere in the
  slice machinery. `fill_subnode` picks the cell within a block by comparing
  it against `ps(igrid)%x`, which is physical. Only `get_igslice`, which turns
  a position into a *block index*, works in the logical coordinate, so that is
  where the logarithm is taken — locally, along with the two bounds checks
  that compare against `xprobmin1`/`xprobmax1`.

  Getting this backwards is a live trap, and it bit this work once: converting
  `slicecoord` at read time fixes the block lookup and breaks the cell lookup,
  whereupon `minloc` matches nothing, `ixslice` comes back 0 and the slice is
  written with all-zero coordinates. Nothing errors. The lesson is that
  `xslice` must mean a physical coordinate throughout the slice path — the two
  halves of that path disagree about the mesh, not about the position.

Nearly everything else was already right, because it keys off `slab_uniform`
(false for any curvilinear build) rather than off uniform spacing: `setdt`
already uses `bgeo%ds`, the flux update already uses `surfaceC`/`dvolume`,
`coarsen_grid` and `prolong_2nd` already volume-weight, `fix_conserve` already
divides by `dvolume`, and the ghost-cell prolongation already does all of its
index arithmetic in the pseudo-uniform `rnode` space — which is precisely the
logical `xi` space.

Two limits worth stating:

- **With the default `log_r0 = 0` a log-radial grid can never reach `r = 0`**,
  so such a `logCylindrical` domain has no polar axis. `set_pole`'s cylindrical
  branch returns immediately in that case, which also avoids a trap: it detects
  the axis by `abs(xprobmin1) < smalldouble`, and `ln(r_min)` is exactly zero
  for the perfectly ordinary domain `r_min = 1`. With `log_r0 > 0` the logical
  coordinate is `ln(1 + r/r0)`, which is zero if and only if `r` is, so the
  same test is not merely safe but exact and the branch runs unchanged.
  `logSpherical` keeps full `theta = 0, pi` support either way; what it does
  **not** get from `log_r0` is a usable origin, see below.
- **The reconstruction stays uniform-stencil.** The MUSCL reconstruction in
  `src/mod_finite_volume.fpp` and the limiters take a scalar spacing and assume
  equal cells, so the second-order constant degrades as `q` grows. Upstream has
  the same limitation everywhere outside its `*NM` WENO family. Keep `q`
  modest — the test cases run at `q ~ 1.04` to `1.09`. A positive `log_r0`
  helps here rather than hurting: the mesh is genuinely uniform for
  `r << r0`, which is exactly the region — next to the axis — where the
  stencil assumption would otherwise be worst.

### Reaching the axis with `log_r0`

`log_r0 = r0 > 0` in `&meshlist` selects the map

    s = ln(1 + r/r0)        r = sign(s) * r0 * (exp(|s|) - 1)

instead of `s = ln(r)`. Three properties of that particular form are what make
it work, and each of them is a deliberate choice rather than an accident:

- **The pivot is the axis.** `s = 0` is `r = 0` *exactly*, so an axis-touching
  domain has `xprobmin1 = 0` after conversion, `rnode`'s corner for the block
  on the axis is exactly zero, and the logical positions of a ghost cell and
  the interior cell it mirrors are exact negatives of one another.
- **The map is odd**, which is the whole point. The pole copy in `getbc` fills
  ghost cell *j* beyond the axis from interior cell *j* of the block half a
  revolution away, so the two have to occupy *mirrored volumes* — otherwise the
  copied cell average belongs to a cell of a different size. The obvious
  alternative `r = exp(s) - r0` is **not** odd: its ghost cells shrink outward
  (`-r0(1 - e^{-jd})`) while the cells they mirror grow (`+r0(e^{jd} - 1)`),
  which misplaces the first ghost cell by about 5% of a cell width at the
  spacings the tests use. That is an O(dx) error at the axis, and it would show
  up as a `2e-4` failure of the `1e-10` `check_pole_ghosts` assertion while
  moving the volume averages in the log not at all. Combined with the pivot,
  the odd form makes the mirror **exact in exact arithmetic**: `abs()` and
  `sign()` carry the negation through, and the barycentre and volume
  expressions are even or odd in the pair of faces. gfortran at `-O3` keeps
  that bit-for-bit; ifx at `-O3` defaults to `-fp-model=fast` and moves the
  closing bits, so `check_pole_mesh`'s `assert_exact` compares to `1e-11`
  rather than exactly (see below).
- **The odd branch is confined to the ghost layer.** `s >= 0` everywhere in the
  physical domain, so the sites that only ever see the domain — the face
  positions in `src/mod_finite_volume.fpp`, and `calc_x`'s corner grid — use
  the cheaper equivalent `r = log_ra*exp(s) + log_rb`, with `log_ra`/`log_rb`
  derived once in `read_par_files` (`1`/`0` when `log_r0 = 0`, so a plain
  logarithmic build gets **bit-for-bit** `exp(s)` and its reference logs do not
  move). Only `fill_geometry_device`, which builds ghost-cell geometry, takes
  the odd branch, on a scalar condition uniform across the kernel.

`log_r0` is a **runtime** parameter, not a fypp define. It is a floating-point
value, so making it compile-time would key the build cache on a number, and
`make/config_reader.py`'s `Value` type accepts only `int`/`string` anyway. It
lives in `mod_global_parameters` alongside `log_ra`/`log_rb` with an
`!$acc declare create`, and `read_par_files` pushes all three to the device
right after converting the domain bounds. Restart safety comes for free:
the snapshot header stores the *logical* `xprobmax1` and compares it with zero
tolerance, and `ln(1 + r_max/r0)` moves with `r0`.

Two things `log_r0` does **not** buy:

- **`logSpherical` still has no usable origin.** The map makes `r = 0`
  reachable, but nothing fills the ghost cells across it: `set_pole` knows only
  about `theta`, and the origin would need a pairing
  `(r, theta, phi) -> (|r|, pi - theta, phi + pi)` — a mirror in *two* indices
  plus a block partner, materially more machinery than the pole copy's
  single-dimension walk. The zero area of the `r = 0` face does not rescue it,
  because the 5-point MUSCL stencil reaches two ghost cells when reconstructing
  the *outer* face of the first interior cell. A `symm`/`asymm` inner boundary
  is the usual approximation and is wrong for flow through the origin.
- **A free lunch on the timestep.** `dr = (r + r0)*ds`, so `r0` *is* the
  near-axis cell width scale, and a small `r0` chosen for resolution buys a CFL
  collapse on top of the `ds(3) = r*dphi -> 0` one that already dominates pole
  runs. `tests/hd/log_cylindrical_pole` runs `r0 = 0.2` on `r` from 0 to 1,
  where the cell width varies by a factor of 5.4 and `dt` is about 2.6 times
  smaller than the uniform `tests/hd/cylindrical_pole`.

### Getting the geometry back to the host

Nothing is copied back by `alloc_node`. That routine runs on every regrid,
mostly for blocks no host routine will look at, so the fetch is placed at the
sites that actually read the geometry instead. Two routines in
`src/mod_geometry.fpp` do it:

- `sync_positions_host` — just the cell-centre positions, which are what host
  code asks for most often. Takes an optional `igrid` for one block, and
  defaults to every grid in `igrids`.
- `sync_geometry_host` — the whole of `geo_t`, positions and metrics alike.
  Takes no arguments: it always covers every grid in `igrids`.

Neither is guarded by a "host is already up to date" flag, deliberately: the
geometry is per block, so any single flag would have to be cleared after
refreshing whatever set of blocks happened to be listed at the time — after
which the next regrid or `selectgrids` can put a block that was not in that
set straight back in front of host code. The bulk form is not cheaper per
block either (the arrays are strided by `igrid`, so it issues one `update` per
block regardless); pick whichever form matches the caller's shape, and only
avoid calling the bulk form from inside a per-block loop.

`sync_geometry_host` is called where `bg%w` is pulled back for output —
`saveamrfile`, next to its `!$acc update host(ps(igrid)%w)` loop — which
covers the log, the snapshot, the slices, the collapsed views and
`autoconvert`; standalone `convert` mode bypasses `saveamrfile`, so
`src/agile.fpp` syncs before `generate_plotfile` too. Individual output
readers (`get_volume_average`, `calc_grid`, ...) therefore do not sync
themselves, and nothing outside output needs the metrics on the host at all:
`fix_conserve`, the one non-output routine that divides by `dvolume`, is a
device kernel and reads `bgeo` there like the rest of the AMR machinery.

`sync_positions_host` is called wherever host code reads `ps(igrid)%x`:
`initial_condition` (per block, in `initlevelone`'s loop, in
`coarsen_grid_siblings`, and in `refine_grids` for each newly created child),
`modify_IC`, `getintbc` (when `usr_internal_bc` is
set), `process`/`process_advanced` (when the `usr_process_*_grid` hooks are
set), `usr_modify_output` in `src/agile.fpp`, `prolong_grid`/`prolong_2nd`
(when `prolongprimitive` is on), and `alloc_node` itself
just before `set_B0_grid`/`phys_set_equi_vars`, which read positions on the
host.

**Adding a host reader of `ps(igrid)%x` or of any metric means adding the
matching sync call.** A missed one reads whatever the host copy held for the
block that previously occupied that `igrid` slot — plausible-looking numbers,
no crash. Readers in code paths this fork does not currently exercise —
`mod_thermal_emission`, `mod_particle_*`, `mod_point_searching`,
`mod_interpolation`, `mod_functions_bfield`, `mod_constrained_transport` —
will need the same treatment when they are brought back into use.

Current limits of the curvilinear (spherical and cylindrical) support:

- `phys = 'hd'`, `phys = 'mhd'`, `phys = 'srhd'`, `phys = 'srmhd'` and
  `phys = 'ffhd'` all implement `addsource_geometry` for both `geometry = 'spherical'` and
  `geometry = 'cylindrical'`, each guarded by its own `#:if GEOM == '...'`
  branch in the physics module's `mod_*_templates.fpp` (a build only compiles
  the branch matching its `GEOM` define). MHD's and SRHD's curvature terms
  (`src/mhd/mod_mhd_templates.fpp`, `src/srhd/mod_srhd_templates.fpp`) were
  ported from upstream MPI-AMRVAC's `mhd_add_source_geom` (cell-centred
  GLM-MHD branch; this fork has no `stagger_grid`/constrained-transport
  support) and `srhd_add_source_geom`, with **every** geometric prefactor
  rewritten to come from the discrete `dAdV = (A_upper - A_lower)/dV` instead
  of the continuous `2/r`, `cot(theta)/r` (spherical) or `1/r` (cylindrical)
  prefactors upstream uses directly. This is not only a well-balancing trick.
  For the metrics `fill_geometry_device` builds, those discrete factors are
  *exactly* the volume averages of the continuous ones,

      spherical:    dAdV(1) = 2*<1/r>       dAdV(2) = <cot(theta)/r>
      cylindrical:  dAdV(1) =   <1/r>

  for any face pair whatsoever, not merely in the continuum limit. Two things
  follow. The pressure terms cancel the pressure part of the flux divergence
  the update actually used, which is the well-balancing; and no curvature term
  depends on where in the cell the stored position sits, which matters because
  `bgeo%x` is the volume barycentre rather than the face midpoint. Evaluating
  `cot(theta)` at the barycentre instead would be about 33% wrong in the first
  cell off the polar axis, since the exact volume average of `cot(theta)` over
  a cell is `cot(theta_midpoint)`. `x` is consequently unused by every
  `addsource_geometry` in the tree; the argument is kept for physics that may
  want it. The well-balancing is exact only where the update's own flux
  divergence stands; refluxing overwrites it at a coarse-fine interface, see
  the refluxing bullet below. SRHD's primitive `mom(:)` slot holds the spatial
  four-velocity `u^i = lfac*v^i` rather than `v^i` itself (see
  `to_primitive`/`to_conservative` and `src/srhd/mod_con2prim.fpp`'s
  `xi = tau + D + p`, `v^2 = S^2/xi^2`), so the curvature terms use the
  primitive `xi_`/`lfac_` auxiliaries directly rather than reconstructing the
  conserved momentum density. SRMHD's (`src/srmhd/mod_srmhd_templates.fpp`)
  combine the two: the momentum terms are the connection terms of the
  non-pressure stress `T^ij = W v^i v^j - (v.B)(B^i v^j + v^i B^j) -
  B^i B^j/lfac^2`, `W = rho h lfac^2 + B^2` — the tensor `get_flux` builds —
  plus `ptot = p + (B^2/lfac^2 + (v.B)^2)/2` in the same slots srhd puts `p`,
  and the induction and GLM `psi` terms are mhd's with `v = u/lfac`. **`W`
  takes the lab-frame `B^2`, not the comoving `b^2`**, which looks like a
  slip but is not: expanding the covariant `(rho h + b^2) u^i u^j - b^i b^j`
  with `b^i = B^i/lfac + lfac (v.B) v^i` cancels the `lfac^2 (v.B)^2 v^i v^j`
  pieces, so only `ptot` carries `b^2`. The primitive `rho` is the proper
  density, so `rho h lfac^2` is formed as `rho*lfac/mu` from the `mu_ =
  1/(lfac h)` auxiliary, exactly as in `get_flux`. Two checks pinned these
  terms: with `B = 0` the srmhd uniform-flow log reproduces srhd's to every
  printed digit, and a uniform field in gas at rest (an exact equilibrium)
  leaves a residual momentum that converges at second order, whereas
  flipping the sign of one magnetic stress term leaves an `O(1e-4)` residual
  that does not converge. Cylindrical's curvature terms only touch the
  `m_r`/`m_phi` (and, for MHD and SRMHD, `B_r`/`B_phi`) components — `m_z`/`B_z` pick up
  none, since a cylindrical volume element's `z`-extent does not depend on
  `r`. FFHD's `addsource_geometry` (`src/ffhd/mod_ffhd_templates.fpp`) is
  empty for both spherical and cylindrical — matching upstream's
  `ffhd_add_source_geom`, which is likewise empty for every coordinate system
  — because its conserved quantities (`rho`, the field-aligned momentum
  `m_par`, the field-aligned energy) are genuine scalars advected along a
  user-supplied field direction `b-hat`, and Gauss's theorem means a scalar's
  flux divergence needs no curvature correction in any coordinate system,
  unlike a vector quantity's components. What *is* geometry-sensitive for
  `ffhd` is unrelated to `addsource_geometry`: the optional `PDIVB` source
  term (`p * div(b-hat)`, compile-time `ffhd_pdivb=T`, default off) computes
  `div(b-hat)` in `addsource_nonlocal` as a plain finite difference over
  `dx(idir)`, which is only correct on a slab-uniform mesh — it is missing
  the curvilinear metric factors a true divergence needs in
  `geometry = 'spherical'` or `geometry = 'cylindrical'`, and this fork has
  not fixed that. Under `LOG_RADIUS` it is further out still, since the radial
  `dx(idir)` it differences over is then a ratio in `ln(r)` rather than a
  length.
- The polar axis is supported for `phys = 'hd'`, `phys = 'mhd'`,
  `phys = 'srhd'`, `phys = 'srmhd'` and `phys = 'ffhd'`; see "The polar axis" below. `ffhd`
  reaches the axis by a different route than the others: its conserved
  quantities are all scalars, so they take no sign flip and the ordinary pole
  copy in `getbc` carries them unchanged, while its frozen field is not
  exchanged through `getbc` at all — `fill_nwextra_device` re-derives it
  analytically in every ghost cell, the axis ghosts included (see "The frozen
  field" below).
- Refluxing across refinement boundaries (`src/mod_fix_conserve.fpp`) works in
  all three coordinate systems. **What a face stores in `pflux` differs by
  geometry, and the two conventions live in the same buffer**, selected by
  fypp rather than by the run-time `slab_uniform` — the two are equivalent
  (`slab_uniform` is true exactly when `coordinate == Cartesian`), but
  `bgeo%dvolume` and `bgeo%surfaceC` are not even allocated in a Cartesian
  build, so the branch has to disappear at compile time:

      Cartesian     store  qdt*f/dx        apply  * CoFiratio = 1/2**ndim
      curvilinear   store  qdt*f*A         apply  / dvolume of the *coarse* cell

  The curvilinear form is *extensive*, which is what makes the fine side work
  without a volume ratio: the four fine face areas sum to the coarse face's
  area exactly for the metrics `fill_geometry_device` builds, so the plain sum
  of four `qdt*f*A_fine` the finite-volume kernel already computes is the
  coarse face's flux, and the coarse cell divides it by its own volume. This
  is upstream MPI-AMRVAC's scheme (`mod_fix_conserve.t`, `mod_finite_volume.t`);
  the areas are the same `bgeo%surfaceC` entries the flux-divergence update
  reads a few lines above the stores.

  Two consequences worth knowing:

  - **Refluxing breaks the exact well-balancing at a coarse-fine interface.**
    The curvature source terms cancel the pressure part of the flux divergence
    *the update actually used* (see the `addsource_geometry` bullet above);
    `fix_conserve` then replaces the coarse face's contribution with the sum of
    the fine ones, and that cancellation no longer holds for that one face. The
    error is second-order and smooth, not a discrete jump, and upstream has the
    same property — but it is why a uniform-flow AMR case such as
    `tests/hd/spherical_pole/uflow_amr.par` no longer stays exactly uniform,
    and why its reference log moved when refluxing was switched on.
  - **The exchange is aggregated per destination rank.** Every chunk bound for
    one peer occupies a contiguous run of `sendbuffer` and travels in a single
    `Isend`, so the message count is the number of neighbouring ranks rather
    than the number of (block, face, child) triples — measured on four ranks in
    `tests/hd/spherical`, 25 outgoing chunks become 2 messages and 42 incoming
    ones become 3. Both ranks order a peer's run by the key `ibuf_offset` is
    indexed by, `4**3*(igrid-1) + inc1 + 4*inc2 + 16*inc3`, which the sender
    derives from the receiving block's index and the receiver from its own
    `igrid`; sorting by it makes the two layouts agree with no handshake, and
    the key also fixes each chunk's size, since `inc_d` lies in `{1,2}` exactly
    when `i_d` is zero. That is what lets the tag drop the block index
    entirely, which in turn retired an earlier workaround that striped tags
    across duplicated communicators to stay under `MPI_TAG_UB`. The one thing
    aggregation gives up is MPI's own mismatch detection, so `fix_conserve`
    checks each arrival's `MPI_GET_COUNT` against the run length it laid out.
  - **A level jump across a polar axis is not refluxed.**
    `init_comm_fix_conserve`, `recvflux`, `sendflux` and `fix_conserve` all
    `cycle` on a non-zero `neighbor_pole`, consistently, so no buffer is sized,
    sent or applied there. Upstream does the same. Level jumps elsewhere in a
    pole case are refluxed normally, and that is enough for conservation in
    the configurations tested: a closed cylindrical domain reaching the axis,
    with the refinement of `tests/hd/cylindrical_pole/uflow_amr.par`, drifts
    in mass by `3.0e-5` with the correction off and by `1.4e-14` with it on.
    That measures the configuration, not the general case — a level jump lying
    *across* an axis face would keep its mismatch.

Validated by ten test directories, one per (physics, geometry) pair, which
is as few as the compile-time parameters allow: `phys` and `geometry` are both
fypp defines, so hd cannot share a build with mhd nor spherical with
cylindrical. Within a pair the cases differ only at run time, so they share one
build and pick their initial condition through the `setup` entry of the
`&usr_list` namelist (read in `usr_init` via `params_read_user(par_files)`, the
mechanism `tests/hd/cloud_crushing` uses):

- `uflow.par`, `setup = 'uniform'` — a uniform Cartesian state, which is an
  exact steady solution. Written in curvilinear components it is a non-trivial
  function of position, so keeping it uniform exercises every curvilinear flux
  and every geometric source term at once, and is the well-balancing test.
  What "uniform" means per physics: a Cartesian velocity for hd; that plus a
  uniform Cartesian magnetic field for mhd, which exercises the induction
  equation and the GLM `psi` too; a constant field-aligned speed along a
  uniform frozen field for ffhd; and a sub-luminal Cartesian velocity for
  srhd, whose Lorentz factor stays a single global constant because it depends
  only on `|v|`, which a rotation leaves invariant; and both of the last two
  for srmhd, whose scalars (`D`, `tau`, `mu`, `lfac`) depend only on `|v|`,
  `|B|` and `v.B` and so are global constants too, while `S = W v - (v.B) B`
  and `B` are rotated cell by cell. srmhd's `blast.par` keeps the uniform
  field, so its blast is magnetised.
- `blast.par`, `setup = 'blast'` — an over-pressured sphere in gas at rest,
  deliberately off-centre in all three coordinates so every momentum component
  is exercised. It is the same setup with the velocity zeroed and a hot spot
  added, which is what lets one initial-condition routine serve both.
- `blast_amr.par` — `tests/hd/spherical` only: the same blast, Lohner-refined.
  The shell keeps moving outward in `r`, so the grid is re-made throughout the
  run rather than settling after the first regrid, which is what puts the
  curvilinear prolongation, coarsening and ghost-cell prolongation paths under
  test.
- `blast_amr_fixgrid.par` — `tests/hd/spherical` only, and **the curvilinear
  refluxing test**. Two changes from `blast_amr.par` give it teeth that no
  other log in the tree has. The domain is *closed* (reflecting in `r` and
  `theta`, `asymm` in `phi`), so the `mean(rho)` and `mean(e)` the regression
  log reports are conserved quantities and any mismatch between a coarse face
  flux and the four fine ones it faces shows up directly as mass or energy
  appearing from nowhere. And `itfixgrid = 0` freezes the mesh after the
  initial refinement, so the shell expands *through* the level jumps.
  That second point is the one that matters: with a live grid the refined
  region is buffered and every jump sits in still-uniform gas, where the
  coarse and fine fluxes agree exactly and refluxing has nothing to correct —
  `blast_amr.par`'s log does not move by a single digit whether the correction
  is applied or not, which is exactly how a broken reflux would slip through.
  Measured with the correction disabled in `src/mod_advance.fpp`, mass drifts
  by `7.3e-4` over this run; with it on, by `4.2e-15`. **A change to the
  refluxing path should be checked against this case, not against
  `blast_amr.par`.**

`agile.par` in each directory is the build reference: `make/config_reader.py`
takes the compile-time parameters from *that file alone*, so it has to declare
the union of what every par file in the directory needs. `specialboundary` is
a runtime logical as well as a fypp define, so `blast.par` leaves it out while
`specialbound_usr` stays compiled in. **Adding a par file to such a directory
is free; adding a directory costs a full rebuild in CI.**

A further **ten** directories, `tests/{hd,mhd,ffhd,srhd,srmhd}/log_spherical` and
`.../log_cylindrical`, repeat the same setups on a logarithmically stretched
radius. They are copies of the directories above, differing in three things:
`geometry` is `'logSpherical'`/`'logCylindrical'`, the radial domain is widened
to `r` from 1 to 4 so that the cell width varies by a factor of about 3.8
across it (the stretching is real, not nominal), and every use of
`xprobmin1`/`xprobmax1` for a physical radius is wrapped in `exp` — see
"Logarithmically stretched radius" above.

Each registers `check_log_grid` as `usr_process_grid`, and **that check is the
point of these directories**, not the reference log. A stretched mesh that is
subtly wrong barely moves a volume average — the same reason the pole cases
assert on ghost cells — so the mesh is compared against its closed form cell by
cell at `it = 0`: that consecutive radial positions are in exact geometric
progression `x(i+1)/x(i) = exp(dxi)`, that `ds(1)` is the physical cell width,
that each stored position is the volume barycentre of its own two faces, and
that `dvolume` is the exact integral over those faces. The reference
expressions are deliberately *different* forms from the ones the geometry
kernel evaluates — the raw quotient of quartic differences rather than the
cancellation-free rewrite, and `cos(theta_L) - cos(theta_R)` rather than
`2 sin(theta_c) sin(dtheta/2)` — so that the two agreeing means something. The
measured agreement is 2 ULP. **A new log-geometry test needs that check, not
just a reference log.**

All ten run at `log_r0 = 0`, the plain `ln(r)`. Their domains stay away from
the polar axis and from `r = 0`; the axis itself is covered by the separate
`*_pole` directories below, of which `tests/hd/log_cylindrical_pole` is the one
that exercises `log_r0 > 0`. Pole and non-pole
cases are kept as separate builds deliberately, even though their compile-time
parameters would have allowed merging them too.

### The polar axis

A spherical domain may span the full `theta` from 0 to `pi`, and a cylindrical
one may start at `r = 0`. The axis is not a boundary: a ghost cell just past
`theta = 0` at azimuth `phi` *is* the interior cell at `theta` and `phi + pi`,
so it is filled from the block half a revolution away rather than by a
boundary condition. Setting it up takes three things in `agile.par`, and
`check_pole_setup` in `src/mod_geometry.fpp` stops the run with a specific
message if any of them is missing:

```fortran
 &boundlist
   typeboundary_min2 = 5*'pole'       ! the axis face(s)
   typeboundary_min3 = 5*'periodic'   ! phi, over a full turn
   typeboundary_max3 = 5*'periodic'
 /
 &meshlist
   xprobmin2 = 0.0d0, xprobmax2 = 0.5d0   ! theta from 0 to pi
   xprobmin3 = 0.0d0, xprobmax3 = 1.0d0   ! phi over 2*pi (2*pi units)
   domain_nx3 = 16, block_nx3 = 8         ! ng3(1) = 2, even
 /
```

- **`'pole'` on the axis face.** `read_par_files` expands it into ordinary
  per-variable `symm`/`asymm` entries in `typeboundary`, which is where the
  ghost-cell exchange later reads the sign from. It is all-or-nothing per
  face. Only `typeboundary_min2`/`max2` may be `'pole'` in spherical and only
  `typeboundary_min1` in cylindrical.
- **`phi` periodic over exactly one full turn.** `set_pole` refuses to
  recognise an axis at all unless `periodB(phi_)`, and the partner block is
  found by shifting the `phi` block index by `ng(1)/2` — half a revolution
  only if `phi` really spans `2*pi`. Upstream MPI-AMRVAC has no check for
  this and silently shifts by the wrong angle on a partial-turn domain;
  `check_pole_setup` does check.
- **An even number of level-1 blocks across `phi`**, for the same `ng(1)/2`
  reason. This one `set_pole` already enforced.

The sign table follows one rule in both coordinate systems: **the component
along the mirrored direction flips, the `phi` component flips, everything else
is symmetric.** For spherical that is `m_theta`, `m_phi` (and `B_theta`,
`B_phi`), with `m_r`/`B_r` symmetric — the same table upstream uses. For
cylindrical it is `m_r`, `m_phi` (and `B_r`, `B_phi`). Note the radial entry:
**upstream marks only `m_phi`/`B_phi` antisymmetric at the cylindrical axis,
which is wrong** — the axis maps `(r, phi)` to `(r, phi + pi)`, under which
`e_r(phi+pi) = -e_r(phi)` exactly as `e_phi(phi+pi) = -e_phi(phi)`, so the
radial component is odd too. With upstream's table the first ghost cell of
`tests/hd/cylindrical_pole_uniform_flow` holds `-0.9713` where the exact value
is `+0.9713`; with this one it is exact to round-off. The indices come from
`iw_mom`/`iw_mag` rather than upstream's positional `rho-mom-[e]-B`
arithmetic, because a `HYPERTC` build registers `q_` between `e_` and `mag`
and shifts every field index.

Everything outside `getbc` was already pole-aware, inherited from upstream:
`set_pole` sets `poleB`, `find_root_neighbor` does the pi-periodic root lookup,
`find_neighbor` exports a `pole(ndim)` flag and suppresses the sibling-index
flip across it, `build_connectivity` records
`neighbor_pole(i1,i2,i3,igrid)` — 0, or the dimension the pole lies across —
and `alloc_node` clears `is_physical_boundary` on a pole face so `bc_phys` is
correctly skipped there. What had to be added is the copy itself, at the six
places in `getbc` (`src/mod_ghostcells_update.fpp`) that move data between
blocks: the same-rank and MPI variants of the same-level, restricting and
prolonging paths. Each now reads `neighbor_pole` and, when it is non-zero:

- **walks the source backwards along that dimension** (`sbase`/`sstep` in the
  code), because the send range runs outward from the axis while the ghost
  range runs inward toward it;
- **multiplies by the sign** from `typeboundary(iw, iB_pole)`, where
  `iB_pole = 2*(ipole-1) + iside`;
- **does not mirror the destination offset in that dimension**, because the
  partner faces you from the same side. In the 0..3 child-offset index this
  reads as: the restricting path swaps 1 with 2 and leaves 0 and 3 alone
  (instead of swapping 0 with 3), and the prolonging path does not swap at all.

For the MPI paths the mirror and the sign are applied on the **send** side, as
upstream does, so the unpack stays a plain strided copy. The same-level send
therefore also has to ship the *destination* offset in its info record rather
than its own, since the receiver cannot see the mirror; the coarse and fine
paths already shipped explicit `inc` indices and only needed the pole variant
of the swap. The prolongation interpolation itself needs no pole branch: by
then `bgc(1)%w` already holds mirrored, signed values.

Two consequences worth knowing:

- **The timestep collapses near the axis.** `ds(3) = r*sin(theta)*dphi` goes
  to zero there, so `setdt` gives a severe CFL restriction. No axis filtering
  or ring averaging is implemented; keep pole cases at modest resolution.
- **Volume averages cannot validate the pole copy.** A cell at the axis has
  vanishing volume and the face lying on the axis has zero area, so a wrong
  ghost value is multiplied by zero on its way into the interior and the jump
  it leaves behind is clipped by the TVD limiter. Reversing *every* sign at
  the pole moves the `mean(w)` the log reports only in the fifth digit. The
  pole test cases therefore assert on the ghost cells directly rather than on
  the log: `check_pole_ghosts` in each case's `mod_usr.fpp` compares them
  against the analytic state at `it = 0`, where the interior is still exactly
  analytic and the copy of it therefore has to be exact to round-off. **A new
  pole test needs that check, not just a reference log.**

Enabled for `phys = 'hd'`, `phys = 'mhd'`, `phys = 'srhd'`, `phys = 'srmhd'`
and `phys = 'ffhd'`. srhd needs no change to the copy itself: its primitive
`mom(:)` slot holds the spatial four-velocity `u^i = lfac*v^i` rather than
`v^i`, but `lfac` is a scalar and therefore invariant under the pole's
pi-rotation, so `u^i` transforms exactly like an ordinary vector and takes
the same `iw_mom`-driven sign table hd's momentum does. srmhd's `u^i` is the
same argument, and its `B` takes mhd's `iw_mag` flips (`read_par_files`
applies them for `physics_type` `'mhd'` or `'srmhd'`); `psi` and the
auxiliaries `mu`/`lfac`, which srmhd does exchange through `getbc`, are
scalars and stay `symm`. Dropping srmhd from the `iw_mag` flip makes
`check_pole_ghosts` fail by ~1 in both of its pole cases. ffhd needs no sign
table at all: `rho`, the field-aligned momentum `m_par` and the energy are
scalars — `m_par = rho*(v·b-hat)` is a contraction of two vectors, invariant
under the proper (det +1) rotation the pole map is — so `read_par_files`
skips the `iw_mom` lines for `physics_type == 'ffhd'` and every
getbc-exchanged variable stays `symm`. The frozen field, the one genuine
vector, is handled entirely outside `getbc`; see "The frozen field" below.

Validated by ten test directories, laid out exactly like the off-axis ones
above and for the same reason — the suite's cost is dominated by compilation,
so cases that agree on the fypp defines share one build and differ only in
their par file. They are kept separate from the off-axis directories even
though their compile-time parameters would have allowed merging: a pole case
and an ordinary curvilinear case are different test families, and mixing them
would make each directory harder to reason about than the extra build is
worth.

- `tests/hd/spherical_pole` — three runs from one build, selected by the
  `setup` entry of the `&usr_list` namelist (read in `usr_init` via
  `params_read_user(par_files)`, the same mechanism
  `tests/hd/cloud_crushing` uses): `uflow.par` is a uniform Cartesian flow on
  a single level, `uflow_amr.par` the same with half the domain in `phi`
  refined so blocks meet across the axis at different levels, and
  `blast_amr.par` a Lohner-refined blast that starts off the axis and expands
  across it. `movie.par` is a longer, visualisable variant of the blast.
- `tests/hd/log_spherical_pole` — the same directory again with
  `geometry = 'logSpherical'`, on the same physical domains so the two are
  directly comparable. This is the only place in the tree where a
  logarithmically stretched radius meets the polar axis. The two features are
  independent by construction — `LOG_RADIUS` touches only the radial
  coordinate and the pole only `theta` and `phi` — but nothing else exercises
  them together, so `check_pole_ghosts` here also calls `check_log_grid`,
  asserting the stretched mesh and the pole copy running on it in one place.
  Its `movie.par` spans `r` from 0.4 to 2.8, a factor of seven in cell width;
  `dr/r` is 0.081 at level one against the uniform run's `dr = 0.1`, so the
  stretched grid is finer everywhere inside `r = 1.2` and coarser only far out
  where the shell is weak.
- `tests/hd/cylindrical_pole` — the same three runs onto the cylindrical axis
  at `r = 0`, with the hot spot placed off the axis in `r` and carried across
  it by a background velocity pointing along `-y`.
- `tests/hd/log_cylindrical_pole` — that directory again with
  `geometry = 'logCylindrical'` and `log_r0 = 0.2`, on the same physical
  domain. It is the only place in the tree where a stretched radius reaches
  `r = 0`, and the reason it exists is `check_pole_mesh`, which asserts to
  **`1e-11`** that each ghost cell beyond the axis has the mirrored position,
  width and volume of the interior cell `getbc` fills it from — the property
  the naive `r = exp(s) - r0` would break, by several percent of a cell
  width, and the log would not notice. The mirror is exact in exact
  arithmetic and bit-exact under gfortran; the tolerance is for ifx's default
  `-fp-model=fast` (see "Reaching the axis with `log_r0`" above). Its `check_log_grid` is the offset-map rewrite of the one
  in `tests/hd/log_cylindrical`: the stored positions are no longer in
  geometric progression (the map is not self-similar any more), so the check
  is that the *cell widths* are, and the face-based barycentre and volume
  assertions carry over unchanged.
- `tests/mhd/spherical_pole` and `tests/mhd/cylindrical_pole` — the uniform
  flow with a uniform Cartesian magnetic field, which puts `B_theta`/`B_r`,
  `B_phi` and the GLM `psi` through the same treatment. The cylindrical cases
  are what pin down the radial sign discussed above.
- `tests/srhd/spherical_pole` and `tests/srhd/cylindrical_pole` — a uniform,
  sub-luminal Cartesian flow, single level, one `uflow.par` each. `mean(rho)`
  in the log is `rho0*lfac0`, not `rho0`, since srhd's conserved density is
  `D = rho*lfac`.
- `tests/srmhd/spherical_pole` and `tests/srmhd/cylindrical_pole` — the srhd
  flow plus a uniform Cartesian magnetic field, whose ghost-cell check covers
  `B`, `psi` and the `mu`/`lfac` auxiliaries as well as the fluid.
- `tests/ffhd/spherical_pole` and `tests/ffhd/cylindrical_pole` — a uniform
  field-aligned flow along a uniform Cartesian frozen field, single level, one
  `uflow.par` each. `check_pole_ghosts` here checks two things at a same-level
  axis neighbour: the fluid ghost cells against the analytic constant (the
  ordinary `getbc` pole copy, which for ffhd only has to deliver a constant),
  and the frozen-field ghost cells rebuilt in Cartesian against the uniform
  `b0` they must reduce to (`fill_nwextra_device` at the axis, which is
  where the sign flips a vector picks up across the pole actually come from).

`agile.par` in each directory is the build reference: `make/config_reader.py`
takes the compile-time parameters from *that file alone*, so it has to declare
the union of what every par file in the directory needs. `specialboundary` and
`refine_usr` are both runtime logicals as well as fypp defines, so an
individual par file can switch them off again while the routines they guard
stay compiled in. **Adding a par file to such a directory is free; adding a
directory costs a full rebuild in CI.**

`check_pole_ghosts` compares the ghost cells against the analytic state at
`it = 0`, where the interior is still exactly analytic and the pole copy of it
therefore has to be exact. How far that can be pushed across a *level jump*
depends on the setup, because there the ghost has been restricted or prolonged
on the way and comparing it against the analytic value at a point is only
meaningful where the field is smooth on the scale of a coarse cell:

- the `uniform` setups are smooth everywhere, so both branches are checked:
  round-off (1e-10) for a same-level pole neighbour and loose (1e-1) across a
  level jump, where the measured error is 3.4e-2 — far below the 0.18 a wrong
  sign produces.
- the `blast` setup has a discontinuous spot surface, and a coarse cell
  straddling it holds the 2:1 average of hot and cold gas, which differs from
  the analytic value at its centre by a large fraction of the jump however
  correct the copy is. Widening the radial domain for `movie.par` put a pole
  face on that surface and the check fired at exactly `13.5/4`, a quarter of
  the initial jump in `e`, on a `neighbor_fine` face. So that setup checks
  same-level neighbours only and covers the restricting and prolonging paths
  through its log instead — which it can afford to do because its log is the
  one that is actually sensitive to them: reversing the `theta` sign makes
  `compare_logs` fail on `mean(m_r)`, `mean(m_theta)` and `mean(m_phi)`, where
  the same break leaves the uniform-flow logs comfortably inside tolerance.

The lesson for a new pole test: point-versus-analytic only works where the
analytic state is smooth on the scale of a coarse cell. Where it is not, check
same-level neighbours and lean on a log that has teeth.

### The frozen field (`phys = 'ffhd'`)

ffhd advects its scalar conserved quantities along a static, user-supplied
unit vector `b-hat`, stored in the `w` slots `b1,b2,b3` (registered by
`var_set_extravar` in `src/ffhd/mod_ffhd_templates.fpp`, so they are the
`nwextra` variables). It is a pure function of position and never a boundary
condition.

This is handled by a generic mechanism, not by ffhd-specific code in the AMR
core. `config_schema.toml`'s `implies` gives `phys = 'ffhd'` the
`FILL_NWEXTRA_ANALYTIC` compile flag, which switches on `fill_nwextra_device`
(`src/amr/mod_amr_solution_node.fpp`) and its call sites. That routine, on the
device, evaluates a by-name user hook

```fortran
pure subroutine usr_set_nwextra(x, wextra)   ! x(1:ndim) in, wextra(:) out
  !$acc routine seq                          ! for ffhd: b-hat, already unit
```

for every cell of every block — interior, inter-block ghosts,
physical-boundary ghosts and polar-axis ghosts alike — and writes the result
verbatim into the extra slots (so a case must return `b-hat` already
normalised; `to_spherical_unit`/`to_cylindrical_unit` do). The hook is called
by name, like `usr_refine_grid` and `gravity_field`, so it is a compile-time
dependency of any `FILL_NWEXTRA_ANALYTIC` build. It runs from a single site:
`alloc_node`, right after `fill_geometry_device`, so every block that comes
into existence — a fresh root, a refined child, a coarsened parent, a
load-balanced arrival — gets its extras there and nothing overwrites them.
`initonegrid_usr` and `usr_special_bc` must **not** set `b1,b2,b3`.

Consequences:

- `nwgc = nwflux` for ffhd (not `nwflux + nwextra`): the frozen field is not
  in the ghost exchange, so `getbc` no longer carries it and `typeboundary`
  needs no rows for it. Prolongation, coarsening, the coarsen MPI buffers and
  the initial host-to-device push of the mesh are likewise bounded by `nwgc`
  rather than `nw`, so they leave the extra slots for `fill_nwextra_device` to
  own (a no-op rewrite for hd/mhd/srhd, where `nwgc == nw`). RK substeps keep
  the extras because they copy `1:nw` wholesale and the flux update never
  writes past `nwflux`.
- **The polar axis needs no special handling for the frozen field.** `bgeo%x`
  in the ghost layer beyond `theta = 0` (or `r = 0`) carries the mirrored
  coordinate — a negative `theta`, or a negative `r` — so evaluating the
  user's analytic field there reproduces on its own the sign flips a vector
  picks up across the axis (`b_r` symmetric, `b_theta`/`b_z` and `b_phi`
  antisymmetric). This is the same identity the `analytic_state` reference in
  the hd/srhd pole tests relies on.
- The old design exchanged `b` through `getbc` (`nwgc = nwflux + nwextra`) but
  filled the *physical-boundary* ghost cells only if the case happened to fill
  the whole `ixI` in `initonegrid_usr` or rewrite `b` in `usr_special_bc`. A
  case with `cont` boundaries and an interior-only IC read uninitialised
  memory there; `tests/ffhd/spherical` and `tests/ffhd/cylindrical` had
  exactly that latent bug in their `blast.par` runs, which is why their
  `correct_output/blast.log` changed when the frozen field moved to
  `usr_set_nwextra`.

## Time integration

Set `time_stepper` and `time_integrator` in `&methodlist`. `advect` in
`src/mod_advance.fpp` implements these; the default `time_integrator` of each
`time_stepper` is in bold, and the overall default is `twostep`, so
`Predictor_Corrector`:

| `time_stepper` | `time_integrator` | order | stages | time levels |
|---|---|---|---|---|
| `onestep` | **`Forward_Euler`** | 1 | 1 | 2 |
| `twostep` | **`Predictor_Corrector`** | 2 | 2 | 2 |
| `twostep` | `ssprk2` | 2 | 2 | 2 |
| `twostep` | `RK2_alfa` (`rk2_alfa` in (0,1]) | 2 | 2 | 2 |
| `threestep` | **`ssprk3`** (`ssprk_order` 3 or 2) | 3 / 2 | 3 | 3 |
| `threestep` | `RK3_BT` (`rk3_switch` 1-4: Ralston, Wray, Heun, Nyström) | 3 | 3 | 4 |
| `fourstep` | **`ssprk4`** (`ssprk_order` 3: Spiteri-Ruuth, 2: SSP(4,2)) | 3 / 2 | 4 | 3 |
| `fourstep` | `rk4` (classical) | 4 | 4 | 4 |
| `fivestep` | **`ssprk5`** (Spiteri-Ruuth (5,4)) | 4 | 5 | 4 |

That is every explicit scheme upstream MPI-AMRVAC has. `read_par_files` also
parses the IMEX family and fills its coefficients, but stops at startup for
it: there is nothing implicit to integrate here, since every
`phys_implicit_update` is the dummy in `src/physics/mod_physics.fpp`.

The SSP schemes are stable to different Courant numbers relative to forward
Euler (ssprk4: 2 at order 3, 3 at order 2), but AGILE does not rescale
`courantpar` per scheme, so a run has to raise it itself to profit.

**Stages and time levels are different numbers, and only the second costs
memory.** `nstep` counts sub-steps (flux evaluations). `nlevel` counts the
solution copies the scheme holds at once, its registers, and `initialize_vars`
allocates `bg(1:nlevel)`. Each level is `nw` times the block size times
`max_blocks` doubles, all device-resident, so it often decides how large
`max_blocks` can be on a GPU. The two counts disagree in both directions:
ssprk4 is 4 stages in 3 levels, while RK3_BT is 3 stages needing 4 levels.
`alloc_node` used to size the levels from the integrator in one place while
`bg` was sized from `nstep` in another, which left RK3_BT and the IMEX schemes
pointing past the end of `bg`.

Both numbers now come from one table, `time_levels` in
`src/mod_global_parameters.fpp`. **A new scheme adds its row there**, a `case`
branch in `advect`, and its name to the accepted list in `read_par_files`.
Everything else follows from `nlevel`: `bg` allocation, the `ps1..ps4` views
that `alloc_node` points at `bg(2..nlevel)` through `ps_level`
(`src/mod_physicaldata.fpp`), and the guard in `alloc_state` that stops on a
view beyond `size(bg)`.

Within `advect` a scheme is written as upstream writes it, with `ps`, `ps1`,
`ps2`, `ps3` becoming time levels 1, 2, 3, 4:

- `advect1(dtfactor, ..., qtC, ka, qt, kb)` adds `dtfactor*dt` times the
  update evaluated on level `ka` to level `kb`. The two must differ, because
  the flux stencil reads neighbouring cells.
- `bg_assign(bgd, a1, bg1[, a2, bg2[, a3, bg3]])` (`bgd = sum`) and
  `bg_accum(bgd, s, a1, bg1[, ...])` (`bgd = s*bgd + sum`) are the pointwise
  combinations between stages. They are generic names over device kernels
  that fypp generates for one to three sources, so each reads only the
  arrays it uses. There are two shapes so that a destination is never also
  passed as a source argument, and is never read where its coefficient would
  be zero.

The schemes are tested on the AMR case of `tests/hd/KH3D-SMR`, not on a
uniform grid, because the static refinement is what puts refluxing at every
stage under test. `tests/hd/KH3D-SMR` itself (basic suite) runs ssprk3,
ssprk4 and the default, Predictor_Corrector: its `pc.par` sets no
`time_stepper`, so the same test also checks the default.
`tests/hd/KH3D-SMR-integrators` (extended suite) is the same case with one
par file per scheme and coefficient variant: all twelve, including
`ssprk_order = 2` for ssprk3 and ssprk4, and two `rk3_switch` tableaus
(Wray, Nyström). It is a directory of its own, rather than par files added
to `KH3D-SMR`, because the suite split is per directory.

Refluxing runs after every stage except the predictor of
`Predictor_Corrector` (as upstream), and only when the mesh has more than one
level (`fix_conserve_at_step`).

## MHD units

- `unit_magneticfield = sqrt(4 pi unit_pressure)` is set in `phys_units`
  (it was left at 1 before, which scaled the Braginskii kappa_perp through
  `htc_Cchi` by the field unit). `H_ion_fr`, `He_ion_fr`, `He_ion_fr2` in
  `&mhd_list` enter b = 1 + H_ion_fr + He_abundance (...), the particles per
  hydrogen nucleus in p = b n kB T; the defaults give the old
  b = 2 + 3 He_abundance, and `He_abundance = 0, H_ion_fr = 0` gives
  Legolas' p = n k T.

## Writing a new simulation case (`mod_usr.fpp`)

Every case directory has a `mod_usr.fpp` defining `module mod_usr` with a
`usr_init()` subroutine that:
1. registers callbacks by assigning the procedure pointers declared in
   `src/mod_usr_methods.fpp` (e.g. `usr_init_one_grid`, `usr_special_bc`,
   `usr_process_grid`, `usr_refine_threshold`, ...),
2. calls `phys_activate()` last.

Hooks that run on the GPU are not procedure pointers but routines the core
calls by name from `mod_usr`, compiled in by a switch in `agile.par`:
`gravity_field` (`*_gravity`), `addsource_usr` (`*_source_usr`; the
`usr_source` pointer is a leftover and never called), `usr_refine_grid`
(`refine_usr`) and `usr_set_nwextra` (ffhd). They are marked with
`${GPU_ROUTINE_SEQ()}$`, so `mod_usr.fpp` starts with

```fortran
#:include "mod_gpu_directives.fpp"
```

by bare file name: `make/30-fypp.mk` passes `-I $(agile)/src` to fypp (and
`make/fypp-deps.py` searches the same directories for its dependency
rules), so the include works wherever the case directory is. The tests used
to write `"../../../src/mod_gpu_directives.fpp"`, which only worked at that
depth inside the repository.

The coordinate system comes from `geometry` in `agile.par` (see "Coordinate
systems" above); `usr_init()` does not need to call `set_coordinate_system`
itself. Older cases may still call it explicitly — the call is idempotent as
long as the name agrees with `geometry` — but new cases shouldn't.

The physics module itself (hd/mhd/ffhd/srhd/srmhd) is selected via `phys` in
`agile.par`'s `&methodlist`, which the config system turns into a `PHYS`
fypp define consumed by `src/physics/mod_physics.fpp` and the per-physics
`mod_*_templates.fpp` files.

## Code architecture

- `src/agile.fpp` — program entry point (`program agile`): MPI init, GPU
  device selection (`set_openacc_device`/`set_openmp_device`, whichever
  backend is active), then `main()` which reads parameters, calls
  `usr_init()`, initializes the AMR tree, and runs `timeintegration()` (the
  main time-stepping loop: `setdt` → optional user process hooks → I/O →
  `advance` → AMR regrid → loop).
- `src/mod_global_parameters.fpp` — the large module of shared global state
  (grid geometry, ghost cells, timers, I/O settings) used throughout.
- `src/mod_gpu_directives.fpp` — fypp macro definitions (`GPU_PARALLEL_LOOP`,
  `GPU_ROUTINE_SEQ`, `GPU_ENTER_DATA_COPYIN`, `GPU_HOST_DATA_USE_DEVICE`,
  etc.) that expand to OpenACC or OpenMP-target directives depending on the
  active backend; `src/mod_gpu_utils.fpp` (module `gpu_utils`) provides the
  `copy_or_update`/`copy_or_update_pointer`/`copy_or_update_alloc` helpers
  used when re-allocating AMR grid data on the device.
- `src/mod_advance.fpp`, `src/mod_finite_volume.fpp` — the finite-volume
  update step; these contain the performance-critical GPU-offloaded loops,
  written via the `GPU_*` directive macros (`GPU_PARALLEL_LOOP`,
  `GPU_ROUTINE_SEQ`, etc.) rather than raw `!$acc`/`!$omp` directives.
- `src/amr/` — the block-based octree AMR machinery: forest/tree bookkeeping
  (`mod_forest.fpp`), refinement/coarsening, load balancing, space-filling
  curve ordering, flux correction at refinement boundaries
  (`mod_amr_fct.fpp`, `mod_fix_conserve.fpp`).
- `src/physics/` — physics-agnostic scaffolding: `mod_physics.fpp` defines
  the procedure-pointer interface (`phys_get_flux`, `phys_to_conserved`,
  etc.) that each physics module fills in via `phys_activate()`; also
  radiative cooling, thermal emission, B0-splitting helpers shared across
  physics modules.
- `src/hd/`, `src/mhd/`, `src/ffhd/`, `src/srhd/`, `src/srmhd/` — one `mod_*_templates.fpp`
  each, guarded by a top-level `#:if PHYS == '...'` fypp block. These are
  **not** standalone Fortran modules: they define fypp macros (physics-specific
  flux functions, primitive/conservative conversions) that get textually
  `#:include`d — via `src/physics/mod_physics_templates.fpp` — into the real
  module files that use them (`mod_physics.fpp`, `mod_physics_vars.fpp`,
  `mod_finite_volume.fpp`, `mod_dt.fpp`, `mod_source.fpp`,
  `mod_radiative_cooling.fpp`). Only one `PHYS` value's branch is active per
  build, so e.g. `mod_physics.fpp` compiles to a different `mod_physics`
  module depending on which physics module was selected.
- `src/io/` — parameter/config file reading (`mod_config.fpp`), snapshot and
  log I/O, `.dat`/VTU conversion (`mod_convert*.fpp`), slices, collapsed
  views.
- `src/particle/` — optional Lagrangian particle module (advection, Lorentz,
  guiding-center, sampling).
- `src/limiter/`, `src/modules/` — flux limiters and small shared utilities
  (interpolation, ODE integration, RNG, timing, LU solve).
- `src/m_octree_mg_*d.fpp` — vendored multigrid solver (from the
  `octree-mg` submodule, `external_libs/octree-mg`) used for e.g. Poisson
  solves in source terms.
- `make/config_reader.py` reads `agile.par` and, together with
  `make/config_schema.toml`, produces the fypp defines / build hash inputs
  in `config.mk`. `make/fypp-deps.py` and `fortdepend` generate `.fpp`/`.f90`
  dependency rules for incremental rebuilds.

## Bug-hunting notes

- **A reduction variable must not be passed by reference into a routine
  called from the reduced loop.** `forcedrefine_grid` combined the per-block
  refine/coarsen flags with `reduction(.or.:refineflag) ...` on the vector
  loop over a block's cells and passed those same flags as `intent(inout)`
  arguments to `usr_refine_grid`. The callee then wrote to the shared
  variable, not the lane's reduction copy, and only the first cell's value
  survived: a `usr_refine_grid` rule was applied to one cell per block, the
  mesh came out lopsided, and nothing complained. CPU builds were unaffected,
  which is what made it look like a physics asymmetry. The call now takes
  per-cell `private` flags that the loop body combines into the reduction
  variables itself; the loop stays vectorised.

- **When CPU is correct and one GPU compiler is wrong and the failure is
  deterministic, it is a miscompile — bisect the kernel's directives
  mechanically, do not reason about which construct "should" be at fault.**
  Toggle one thing at a time in the offending kernel (collapse level;
  array-section assignment vs explicit scalar loop; `vector` placement; loop
  bounds; `-O` level) and watch a bucketed error metric — e.g. split the
  ghost-region error into face / edge / corner, or per variable. This is far
  faster than theorising: the polar-axis GPU bug was chased through several
  wrong hypotheses (the MPI path, a reversed array section) before mechanical
  bisection pinned it to `collapse` in one `!$acc loop` in three toggles.
  Cross-check with `arch=gnu` (gfortran) and `arch=nvidia` *without*
  `OPENACC=1` (nvfortran, host) to confirm it is the OpenACC offload and not
  the front end.

- **In a `collapse`d `!$acc loop` inside an `!$acc routine`, a `call` in the
  loop body needs every one of the body's `private` automatic arrays to have
  a compile-time size.** nvfortran's OpenACC miscompiles a collapsed `vector`
  loop — reached from a gang-level `!$acc parallel loop` through a `!$acc
  routine vector` — when the body both `call`s another `!$acc routine` *and*
  declares a `private` automatic array whose size is a run-time value: it
  gets that array's per-lane stride wrong and scatters the results into the
  wrong cells (the ghost layer comes back index-rotated; results drift rather
  than crash). Both conditions are needed. Remove the `call` (inline the
  callee) *or* give every such `private` array a `parameter` bound and the
  same `collapse(3)` is generated correctly. The reduced reproducer is on
  branch `reproducer/issue-154` (`tests/regression/issue154/repro154.f90`);
  its variant 4 keeps `collapse(3)` and the `call` and only changes
  `wpt(1:nw)` to `wpt(1:nw_phys)`.

  This bit every hd curvilinear `*_pole` case's `specialbound_usr`, whose
  `wpt(1:nw)` private buffer is sized by the run-time `nw` — logs stayed put
  but `check_pole_ghosts` failed by ~0.4 on nvfortran while gfortran and
  nvfortran-without-OpenACC were exact. **The fix now in the tree** is
  `wpt(1:nw_phys)` (the compile-time conserved-variable count from the
  physics templates, `= nw` for these cases), `collapse(3)` kept. The
  mhd/srhd pole cases call a helper too but their `private` arrays
  (`v`, `b`, `x_loc`) are already compile-time sized, and the ffhd pole
  cases have no `call` in the loop at all, so all six also keep
  `collapse(3)`. Dropping the `collapse` (`!$acc loop vector` on the
  innermost loop, outer loops sequential) also stops the corruption and was
  the first fix shipped, but the reproducer shows it still silently reorders
  a short vector loop that has a `call` and a run-time-sized `private` — it
  only looked correct because the pole cases' boundary state does not vary
  along that axis. Collapse is fine at a top-level `!$acc parallel loop`.
  This is very likely the same defect as issue #154
  (`to_primitive`/`to_conservative` called from a `collapse`d loop in
  `specialbound_usr` returning nonsense).

- **Under NVHPC the OpenMP gang/vector macros expand to `target teams loop`
  and `loop bind(parallel)`, not `teams distribute` and `parallel do`.**
  nvfortran compiles `teams distribute` with a nested `parallel do` in
  generic mode: the inner region becomes an outlined device function (the
  `NVFORTRAN-W-0155` "Passing external procedure as argument" warning), and
  every kernel takes ~255 registers per thread and ~9 KB of shared memory,
  however simple. On a GTX 1660 Ti, KH3D at 192^3 took 58.0 s against
  OpenACC's 45.3 s, with the ghost-cell kernels 4x slower; with `loop` it
  took 46.4 s. `bind(parallel)` is not optional: without it nvfortran ran its
  own dependence analysis on the finite-volume kernel, gave up on its
  private temporaries, and launched it with one thread per team, 25x slower,
  without any error. A returning `W-0155` warning means the generic-mode
  pattern is back. Other compilers keep `distribute`/`parallel do`, and
  should until benchmarked: LLVM flang (and so AMD's `amdflang`) lowers a
  nested `loop bind(parallel)` to a worksharing `do` with no enclosing
  `parallel` (`GenericLoopConversion.cpp`), which leaves one thread per
  team, and Cray CCE documents the same `bind(parallel)` -> `do` mapping
  with no `loop`-specific optimisation. MFC's `omp_macros.fpp` likewise uses
  `loop` for NVHPC only.

## API documentation (FORD)

Source-level API docs are generated with [FORD](https://forddocs.readthedocs.io/)
from the `!>`/`!<` doc comments (see `doc/documentation.md`), replacing the
old Doxygen setup:

```bash
uv sync --extra docs
source .venv/bin/activate
doc/ford/build_docs.sh   # preprocesses src/*.fpp with fypp, then runs ford
firefox doc/ford/html/index.html
```

`build_docs.sh` preprocesses every file once per physics module
(`doc/ford/docgen_{mhd,hd,ffhd,srhd}.par`), since the `#:if PHYS == '...'`
splicing described above means one pass only shows one physics module's
variant of the shared files; `merge_physics_variants.py` keeps files that
come out identical once and suffixes the rest (`mod_physics_hd`, ...).
`doc/ford/src/` and `doc/ford/html/` are generated/gitignored, not committed.
The narrative `doc/*.md` pages are part of the site (`page_dir: ..` in
`doc/ford/project.md`), as FORD pages with a `title:` header; the time
integrators are described in `doc/time_integration.md`. Page navigation on the site is
a "Contents" drawer rendered from `doc/contents.md` itself
(`doc/ford/templates/info_page.html`), so a new page must be listed there to
be reachable.

## Code style (from `doc/code_style_guide.md`)

- Fortran 90/95/2003-era style; no implicit typing, no common blocks.
- snake_case names, self-explanatory beyond trivial loop indices.
- 2-space indent for `program`/`module`/`subroutine`/`function`/`associate`;
  3-space for `type`/`interface`/`do`/`if`/`select case`/`where`/`forall`.
- Use `dp` kind and `_dp` suffixes for double precision, not
  `double precision`/bare literals.
- Array loops ordered so the leftmost index varies fastest (Fortran
  column-major order).
- All procedure arguments must declare `intent`.
