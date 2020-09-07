# Richtmeyer Meshkov 2D test problem

\test 2D HD test: Shock interaction with a contact discontinuity

# Introduction

This is a 2D test (as studied e.g. in Delmont et al, JFM 2009, 627, 33-53) where a planar shock impinges
on an inclined contact discontinuity. The evolution in the early stages as far as the refraction pattern is concerned, can be predicted analytically, while the long term evolution shows Richtmeyer-Meshkov instability
taking place on the shocked contact interface. The AMR will be beneficial to get highly detailed finestructure.

# How to run

## Setting up the files

    $AMRVAC_DIR/setup.pl -d=2

## Compiling

Simply issue the `make` command:

    make

## Running the code

To run with e.g. 4 processors, use

    mpirun -np 4 ./amrvac -i amrvac.par

# Changing the parameters

A complete list of parameters can be found [par.md](par.md).

# Changing the physics and initial conditions

Have a look at the local file `mod_usr.t`. You can modify the 
initial conditions easily, try setting up shocks of different Mach number or change the impact angle 
and density contrast at the density discontinuity.  Verify how different schemes and resolutions perform on this test. A relevant extra variable in output is a Schlieren plot of the density (gradient).
