# Woodward-Collela 2D test problem

\test 2D HD test: shock reflection

# Introduction

This is a standard 2D shock reflection test from Woodward & Collela (JCP 54, 1984).
It serves as an example on the usage of non-trivial boundary prescriptions: the bottom boundary is
partly fixed, partly solid wall, the top boundary uses a time-dependent prescription for 
the pre- and post-shock region. It is a nice test for AMR, as more resolution will show more details in the later stages of this shock reflection problem.

# How to run

## Setting up the files

    $AMRVAC_DIR/setup.pl -d=2

## Compiling

Simply issue the `make` command:

    make

## Running the code

To run with e.g. 4 processors, use

    mpirun -np 4 ./amrvac -i wc_2d.par

# Changing the parameters

A complete list of parameters can be found [par.md](par.md).

# Changing the physics and initial conditions

Have a look at the local file `mod_usr.t`. You can modify the 
initial conditions easily, try setting up shocks of different Mach number (and change the boundary prescription accordingly). Verify how different schemes and resolutions perform on this test. The extra variable in output is a Schlieren plot of the density (gradient).
