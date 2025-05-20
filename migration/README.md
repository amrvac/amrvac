Work in progress...

This folder contains the migrated make files.

Inside the `make` directory there are a set of make files that should be sourced in numerical order.

- prelude: common definitions
- find-agile: make sure that `$(amrvac)` is defined from the `AMRVAC_DIR` variable.
- select-arch: the commandline should get a `arch=` argument, or architecture defaults to `gnu`. If there was a previous build, use the same value.
- read-config: there should be a `config.mk` located from where the make command is run.
- build-dir: prepares the build directory and symlinks to the latest build
- dependencies: runs the Python script `fortdepend` which scans dependencies and turns those into a `dependencies.mk` file
- fypp: configures fypp
- compile: defines compile rules
- link: defines linking rule(s)


To install this version, follow these steps:
- install the uv package manager `pip install uv`
- make sure `$AMRVAC_DIR` points to the repository root folder
- go into the `migration` folder `cd migration`
- install the required python packages `uv sync` and activate them `source .venv/bin/activate`
- generate the fpp files for migration by running `make` once in the migration folder (generates src folder)
- go into a test `cd test/benchmark_KH3D`
- to compile, load the appropriate modules, e.g. on snellius:
```
module purge
module load 2023
module load OpenMPI/4.1.5-NVHPC-24.5-CUDA-12.1.1
```
- compile with nvfortran and activated OPENACC via `make arch=nvidia OPENACC=1`
 
