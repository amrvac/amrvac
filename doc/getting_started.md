# Getting started

[TOC]

This page shows how to run a first problem after the @ref installation.md is
completed.

# The VAC problem

Traditionally, the first test problem is the VAC advection located in
`$AMRVAC_DIR/tests/rho/vac`. Two files are present in this folder:

* `mod_usr.t`: the user code for this problem (defining e.g. initial conditions)
* `amrvac.par`: a text file with settings

# Compilation {#vac_compilation}

In the VAC test folder, run the setup script with:

    $AMRVAC_DIR/setup.pl -d=2

This will copy a makefile to the current folder, and set the problem dimension
to two. To select different compiler settings (present in `$AMRVAC_DIR/arch/`)
you could for example use:

    $AMRVAC_DIR/setup.pl -d=2 -arch=intel
    $AMRVAC_DIR/setup.pl -d=2 -arch=debug

You can also edit the `makefile` directly.

Then its time to compile the code:

    make -j 4

This will perform a parallel compilation using 4 cores. First, the AMRVAC
library is compiled in `$AMRVAC_DIR/lib/` if it is not available already. This
is the generic part of AMRVAC that does not depend on your user code. Then your
user code is also compiled, and a binary called `amrvac` is produced.

There are a couple useful commands to know about:

    make clean        # clean the local object files
    make ARCH=debug   # do a debug build (extra error checking)
    make ARCH=<name>  # use the compilation flags from arch/<name>.defs
    make allclean     # clean the local object files and the AMRVAC library

# Running the test {#vac_run}

To run the test problem on 4 cores:

    mpirun -np 4 ./amrvac -i amrvac.par

This will run amrvac with parameter-file amrvac.par. Then you have several new
files ending in `.dat`, `.vtu` and one `.log` file. The `.dat` files are used
for restarts and the `.vtu` files contain the output data to be visualized.

For more information about command line parameters

# Visualization {#visualization}

Simulation output in `.vtu` (VTK unstructured) format can directly be
visualized [Paraview](http://www.paraview.org/)
or [Visit](https://wci.llnl.gov/simulation/computer-codes/visit).

Visualization or analysis of the results can also be done **a posteriori**, by
converting `.dat` files to e.g., IDL, openDX, Tecplot, or VTK native formats.
This requires the same amrvac executable used to run the simulation.
