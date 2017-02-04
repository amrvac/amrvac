# Rayleigh-Taylor test problem

# Introduction

The amrvac setup in this folder demonstrates
the
[Rayleigh-Taylor instability](https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Taylor_instability),
which occurs when a lighter fluid pushes into a heavier fluid.

# How to run

## Setting up the files

    $AMRVAC_DIR/setup.pl -d=3

## Compiling

Simply issue the `make` command:

    make

When switching between 2D and 3D, perform a `make clean`.

## Running the code

    mpirun -np 4 ./amrvac -i rt_3d.par

to run with the settings for 3D.

# Changing some of the settings

Some of the settings that you could change in the `.par` files are:

name | description
---|---
`filenameout` | Base file name for output
`dtsave(1)` | Time between log output
`dtsave(2)` | Time between dat/vtu output
`mxnest` | the maximum number of refinement levels
`dixB` | number of ghost cells, depends on spatial discretization
`nxlone[1,2,3]` | the size of the coarse grid, needs to be divisible by (block size - 2 * dixB)
`tmax` | the end time of the simulation
`typeadvance` | time discretization method, e.g., twostep, threestep, ssprk43
`courantpar` | CFL number (see also typecourant)
`typefull1` | spatial discretization method, e.g., tvd, tvdlf, hlcc
`typelimiter1` | which limiter to use in the spatial discretization, e.g., cada3, koren, woodward

A more complete list of parameters can be found [par.md](par.md).

# Changing the physics and initial conditions

Have a look at the local file `mod_usr.t`. You can for example modify the
following routines:

* `initonegrid_usr`: change the initial conditions
