# Implosion test problem

# Introduction

This is a problem on a closed domain involving shock front formation and
interaction, enriched by possible Richtmyer-Meshkov instability development.
Quantify for the latter the deposition/growth of vorticity (its component
perpendicular to the plane). Verify that conservation is perfect when using a
conservative discretization scheme (inspect/use the log file).

The 2D setup is discussed in Liska & Wendroff, SIAM J. Sci. Comput. 25 , 3, pp.
995-1017, 2003 (section 4.7). You are invited to vary the basic parameters
(pressure density contrast) or to turn the setup in true 3D variant. You should
look especially at both short term and long-term symmetry preservation of the
solution (across the diagonal in the 2D setup). The full implosion problem (not
just the upper right quadrant) is also used in Vandenbroucke & De Rijcke,
Astronomy & Computing 16 , 109 (2016), [Figure 2] to illustrate the advantages
of a moving mesh code (shadowfax) as compared to a Eulerian grid code.

# How to run

## Setting up the files

Use energy equation:

    $AMRVAC_DIR/setup.pl -d=22 -g=20,20 -p=hd -eos=default -nf=1

Using a polytropic equation of state:

    $AMRVAC_DIR/setup.pl -d=22 -g=20,20 -p=hd -eos=iso -nf=1

When switching, perform a `make clean` before recompiling.

## Compiling

Simply issue the `make` command:

    make

When switching between 2D and 3D, perform a `make clean`.

## Running the code

To run with e.g. 4 processors use

    mpirun -np 4 ./amrvac -i implosion.par

For the polytropic version, use

    mpirun -np 4 ./amrvac -i implosion_iso.par

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
