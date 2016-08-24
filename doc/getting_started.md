# Getting started

[TOC]

# Introduction {#started_intro}

Here we will get you started with MPI-AMRVAC. Fortunately, the code has very
little dependencies: the only thing you will really need is an implementation of
MPI, e.g. [open-mpi](http://www.open-mpi.org/) or
[mpich](http://www.mpich.org/).

# Get the code {#get_code}

We will install the code to the `~/codes/amrvac` directory. The following
commands should do the trick:

    cd ~
    mkdir codes
    cd codes

You can then clone the repository using git:

    git clone https://gitlab.com/mpi-amrvac/amrvac.git
    cd amrvac

Alternatively, you can also download a zip file of the most recent version from
[Gitlab](https://gitlab.com/mpi-amrvac/amrvac/tree/master), by clicking the
"Download zip" button.

# Installation {#install}

The only thing that should be done after unpacking is to set the `$AMRVAC_DIR`
environment variable holding the path to the code. To do so using bash, you
should add the following entry to your `~/.bash_profile` (or `~/.bashrc`) file:

    # AMRVAC:
    export AMRVAC_DIR=$HOME/codes/amrvac

To use scripts more conveniently, the following line can be placed in `~/.bash_profile`:

    PATH="$AMRVAC_DIR:$AMRVAC_DIR/tools:./:$PATH"

Don't forget to source the `.bash_profile`:

    source ~/.bash_profile

# Running a test problem {#running_test}

Traditionally, the first test problem is the VAC advection located in
`$AMRVAC_DIR/tests/rho/vac`.
In this folder, you can type

    $AMRVAC_DIR/setup.pl -d=22 -g=16,16 -p=rho -u=testrho

which sets up a 2D (-d=22) advection (-p=rho) problem for an AMR-grid with
16x16 blocks (-g=16,16). This specific problem is pre-defined (-u=testrho).
You will notice some `*.t` files in your directory:

* `amrvacsettings.t`
* `mod_indices.t`
* `amrvacusr.t`
* `amrvacusrpar.t`

Normally, you only have to care about the last two. In `amrvacusr.t`, the
whole problem for example initial conditions and boudary conditions are set
up. `amrvacusrpar.t` can be used to provide additional global varibles for
your setup. You will learn more about this later, lets first build the code.

The machine specific definitions are outsourced to the directory
`$AMRVAC_DIR/arch`. There you have a number of pre-defined make rules for
various compilers. If you are using intel, you should well be served with
`default.defs`. To tell amrvac to use these definitions:
`$AMRVAC_DIR/setup.pl -arch=default`

Then its time to compile and run the code:

    make
    mpirun -np 2 ./amrvac -i testrho_vac22 < /dev/null > out &

This will run amrvac on two cores of your machine using the pre-defined
parameter-file testrho_vac22 (-i ) and route output to the file out. This
should take around 15 seconds. Then you have several new files ending in
`.dat`, `.vtu` and one file called `amrvac.log` - this is the log file updated
during computation.
The `.dat` files are used for restarts (`./amrvac -restart n`) and the `.vtu`
files contain the output data to be visualized with e.g. Paraview (we also
provide Python interfaces though).

Congratulations, you have run your first simulation with amrvac!
