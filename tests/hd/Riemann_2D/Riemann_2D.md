# Riemann 2D test problems

\test 2D HD Riemann problems

# Introduction

These are 2D Riemann problems as in Lax & Liu SIAM 1998, or also in Liska & Wendroff, SIAM J. Sci. Comput.
25, 3, pp. 995-1017, 2003 (section 4.4). You can try to modify the setups, or
compare different schemes. This is a direct generalization of the 1D Riemann
problems, so you can try and interpret 1D versus 2D
aspects (involving shocks, contact discontinuities, rarefaction waves). You can
also make a 3D version.

# How to run

## Setting up the files

    $AMRVAC_DIR/setup.pl -d=2

## Compiling

Simply issue the `make` command:

    make

## Running the code

To run with e.g. 4 processors, use

    mpirun -np 4 ./amrvac -i rm_2d.par

# Changing the parameters

A complete list of parameters can be found [par.md](par.md).

# Changing the physics and initial conditions

Have a look at the local file `mod_usr.t`. You can modify the 
initial conditions easily.
