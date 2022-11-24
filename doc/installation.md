# Installation

[TOC]

# Requirements {#requirements}

MPI-AMRVAC is in use on Mac laptops, Unix desktops, and modern supercomputing
facilities. The following software is required:

* A Fortran compiler, e.g. gfortran or ifort
* An MPI library, e.g. OpenMPI or intel MPI
* `perl` for the [VACPP](vacpp.md) preprocessor (standard on Unix-based platforms)
* `git` to download and update the code

# Getting the code {#get_code}

To install the code in folder `~/codes/amrvac`:

    cd ~/codes
    git clone https://github.com/amrvac/amrvac.git

# Setting path {#install}

To use MPI-AMRVAC, the installation directory must be exported as an environment 
parameter called AMRVAC_DIR and this variable can be added to the command searching path.
To do so using bash, you should add the following line
to your `~/.bashrc` file (or perhaps `~/.profile` in mac OS):

    export AMRVAC_DIR=$HOME/codes/amrvac
    PATH="$PATH:$HOME/codes/amrvac:$HOME/codes/amrvac/tools"

You can then activate the setting by sourcing this file (or you can close and open a new terminal shell), 

    source ~/.bashrc

# Update the code {#update_code}

It is recommended to update the code regularly to benefit from new improvements and new features 
which are continuously developed. To update the code:

    cd ~/codes/amrvac
    git pull

If you modified any source code of amrvac, which has an updated version, "git
pull" will fail to update the code. You can give up your modification by 
"git checkout your_modified_file" after making a copy of your own version (if you want to keep it), then update the code. If you
want to contribute your modification of the code, please read @ref contributing.md for more information.     
