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
    git clone https://github.com/amrvac/amrvac.git

# Setting path {#install}

To use MPI-AMRVAC, the installation directory has to be added to the command searching path. 
To do so using bash, you should add the following entry
to your `~/.bashrc` file (or perhaps `~/.profile` in mac OS):

    PATH="$PATH:$HOME/codes/amrvac:$HOME/codes/amrvac/tools"

You can then activate the setting by source the file, or close and open a new terminal shell, 

    source ~/.bashrc
