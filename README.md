[![tests](https://github.com/amrvac/AGILE-experimental/actions/workflows/tests.yml/badge.svg)](https://github.com/amrvac/AGILE-experimental/actions/workflows/tests.yml)

# AGILE

This is the development version of AGILE, a GPU enabled fork of MPI-AMRVAC. All documentation is available on [amrvac.org](http://amrvac.org/).

## Work in progress...

To install this version, follow these steps:
- install the uv package manager `pip install uv`
- make sure `$AMRVAC_DIR` points to the repository root folder
- install the required python packages `uv sync` and activate them `source .venv/bin/activate`
- go into a test `cd test/KH3D`
- to compile, load the appropriate modules, e.g. on snellius:
```
module purge
module load 2023
module load OpenMPI/4.1.5-NVHPC-24.5-CUDA-12.1.1
```
- compile with nvfortran and activated OPENACC via `make arch=nvidia OPENACC=1`

## Currently supported features on master
- Hydrodynamics module (hd)
- Source terms (gravity) and boundary conditions (standard `symm, asymm,` etc.)
- Multi-block uniform grid simulations
- Multi-GPU (MPI) simulations
