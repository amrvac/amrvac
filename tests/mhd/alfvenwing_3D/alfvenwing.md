# Alfven wing test problem

# Introduction

The interaction of Io's atmosphere with Jupiter's magnetospheric plasma is a representative example of a sub-Alfvénic interaction scenario. This example shows Io's plasma interaction and the development of Alfvén wings for a very simple 3D Cartesian setup. The configuration considers a neutral gas cloud of a spherical shape which represents Io's atmosphere. The sphere is placed with its center at the origin and has a radius of one unit. The plasma flows past the gas cloud and interacts via a constant ion-neutral collision frequency with the gas cloud. Source terms representing the elastic collisions between neutrals and ions are added to the momentum and energy equations.
Io's plasma interaction is studied e.g. in Bloecker at al. JGR 2018 and the output of this test should be similar to Figure 3 a,c,e in Bloecker at al. JGR 2018.
The results of the test can be compared to the analytical model of the Alfvén wing, e.g. Saur et al. JGR 1999 (see Appendix A2) or Neubauer JGR 1980.

# How to run

## Setting up the files


    $AMRVAC_DIR/setup.pl -d=3

## Compiling

Simply issue the `make` command:

    make

## Running the test

To run the test problem on e.g. 8 cores:

    mpirun -np 8 ./amrvac -i Io.par


# Changing the parameters

Some of the settings that you could change in the `Io.par` file are:

name | description
---|---
`base_filename` | Base file name for output
`xprobmin[1,2,3]` | limits of the computational domain
`xprobmax[1,2,3]` | limits of the computational domain
`domain_nx[1,2,3]` | number of grid cells per dimension
`time_max` | the end time of the simulation

It is recommended to increase the number of grid cells per dimension `domain_nx[1,2,3]`, e.g. 160 grid cells. An increased resolution will give you better results of the plasma interaction. The simulation should be performed until the Alfvén wings reach the outer boundary of the box. 

A complete list of parameters can be found [par.md](par.md).

# Changing the physics and initial conditions

Have a look at the usr_list in file `Io.par`. You can modify the initial conditions (in SI units). You can e.g. change the direction of the magnetic field to vary the propagation direction of the Alfvén waves or the magnitude of the magnetic field in order to change the Alfvén angle. Similar plasma interaction occurs at Jupiter's moon Europa. You can try the test with initial conditions for the plasma interaction with Europa (Jupiter's moon).