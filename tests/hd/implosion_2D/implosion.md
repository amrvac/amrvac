# Implosion test problem: 2D HD Liska-Wendroff

\test 2D implosion Liska-Wendroff

# Introduction

This is a problem on a closed domain involving shock front formation and
interaction, enriched by possible Richtmyer-Meshkov instability development.
One can quantify for the latter instability the deposition and growth of vorticity 
on the contact discontinuity (in 2D: use its component perpendicular to the plane). 
Verify that conservation is perfect when using a
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

    $AMRVAC_DIR/setup.pl -d=2

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
`base_filename` | Base file name for output
`dtsave_log` | Time between log output
`dtsave_dat` | Time between dat/vtu output
`refine_max_level` | the maximum number of refinement levels
`nghostcells` | number of ghost cells, depends on spatial discretization
`domain_nx1,domain_nx2` | the size of the coarse grid
`time_max` | the end time of the simulation
`time_stepper` | time discretization method, e.g., twostep, threestep, fourstep
`courantpar` | CFL number (see also typecourant)
`flux_scheme` | spatial discretization method, e.g., tvd, tvdlf, hlcc
`limiter` | which limiter to use in the spatial discretization, e.g., cada3, koren, woodward

A complete list of parameters can be found [par.md](par.md).

# Changing the physics and initial conditions

Have a look at the file `mod_usr.t`. You can modify the initial conditions, or try changing the adiabatic index \f$ \gamma \f$
