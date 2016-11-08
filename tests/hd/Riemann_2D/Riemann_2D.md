# Rayleigh-Taylor test problem

# Introduction

These are 2D Riemann problems as also in Liska & Wendro , SIAM J. Sci. Comput.
25, 3, pp. 995-1017, 2003 (section 4.4). You can try to modify the setups, or
compare different schemes. This is a direct generalization of the 1D Riemann
problem discussed in the course, so you can try and interpret 1D versus 2D
aspects (involving shocks, contact discontinuities, rarefaction waves). You can
also make a 3D version.

# How to run

## Setting up the files

    $AMRVAC_DIR/setup.pl -d=22 -g=20,20 -p=hd

## Compiling

Simply issue the `make` command:

    make

## Running the code

To run with e.g. 4 processors, use

    mpirun -np 4 ./amrvac -i rm_2d.par

# Changing the parameters

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

A complete list of parameters can be found [par.md](par.md).

# Changing the physics and initial conditions

Have a look at the local file `amrvacusr.t`. You can modify the following
routines:

* `initglobaldata_usr`: change the adiabatic index \f$ \gamma \f$
* `initonegrid_usr`: change the initial conditions
