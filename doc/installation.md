# Installation

[TOC]

# Requirements {#requirements}

MPI-AMRVAC is in use on Mac laptops, Unix desktops, and modern supercomputing
facilities. The following software is required:

* An MPI library, e.g. OpenMPI or MPICH
* A Fortran compiler, e.g. gfortran or ifort
* `perl` for the [VACPP](vacpp.md) preprocessor (standard on Unix-based platforms)
* `git` to download and update the code

# Getting the code {#get_code}

To install in `~/codes/amrvac`:

    cd ~/codes
    git clone https://gitlab.com/mpi-amrvac/amrvac.git

# Setting AMRVAC_DIR {#install}

To use MPI-AMRVAC, the `$AMRVAC_DIR` environment variable has to point to the
installation directory. To do so using bash, you should add the following entry
to your `~/.bashrc` file (or perhaps `~/.bash_profile`):

    # Set AMRVAC_DIR
    export AMRVAC_DIR=$HOME/codes/amrvac

You can then open a new shell, or source the file with

    source ~/.bash_profile


