# Cartesian finite-volume magnetic helicity

This small demo reuses the validated Cartesian TDm/RBSL magnetic rope used by
the magnetic-topology examples, but runs the independent MPI finite-volume
relative-helicity backend.  It deliberately contains only the user module and
two parameter files; `setup.pl` generates build files locally and `output/`
is only a runtime directory.

## Build and generate a snapshot

From this directory, set `AMRVAC_DIR` to the repository root and build a private
three-dimensional executable:

```bash
export AMRVAC_DIR=/absolute/path/to/amrvac
cd "$AMRVAC_DIR/tests/demo4/MagneticHelicity_Cart"
"$AMRVAC_DIR/setup.pl" -d=3 -v=3
make AMRVAC_DIR="$AMRVAC_DIR"
mkdir -p output

OMP_NUM_THREADS=1 OMP_DYNAMIC=false \
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun --bind-to none --map-by slot -np 1 ./amrvac \
  -i amrvac_generate.par
```

The light default mesh is `64 x 48 x 64`.  For a u128-style run, change the
three `domain_nx*` values to `128, 96, 128` and the block sizes to `32` in
`amrvac_generate.par`.

## Run the helicity conversion

Use any MPI size; the backend is not the OpenMP field-line topology backend.
The macOS Open MPI form below avoids core-binding failures:

```bash
OMP_NUM_THREADS=1 OMP_DYNAMIC=false \
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun --bind-to none --map-by slot -np 2 ./amrvac -convert \
  -if output/tdm_0000.dat -i amrvac_generate.par amrvac_helicity.par
```

The result is `output/tdm_helicity.csv`.  It contains code and physical-unit
`Hm`, `HJ`, `HPJ`, magnetic energy `E`, potential energy `Ep`, free energy
`Efree`, the ratio `abs_HJ_over_abs_Hm`, and quality fields for flux balance,
boundary normals, curl reconstruction, decomposition closure, divergence, and
multigrid residual/cycles.  `mh_write_debug_vti` is off by default; enable it
only when a streamed `B_p`, `A`, `A_p`, and density-field VTI is needed for
numerical inspection.

The six-face Neumann reference solve requires a uniform, unstretched,
three-dimensional Cartesian mesh and a compatible net boundary flux.  AMR
snapshots must first be remapped with `level_io` to a fixed uniform level.
`convert_type='magnetic_topology'` remains a separate single-MPI-rank,
OpenMP field-line workflow.  The sampled TDm/RBSL faces have a measured
discrete relative flux imbalance of about `1.54e-5`; the demo therefore uses
`mh_max_flux_imbalance=1.d-4` and keeps the measured value in the CSV quality
columns.  Analytic or cleaned input data can use the stricter `1.d-6` guard.
