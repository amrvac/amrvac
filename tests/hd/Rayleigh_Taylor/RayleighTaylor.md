# Rayleigh-Taylor test problem

# Introduction

The amrvac setup in this folder demonstrates
the
[Rayleigh-Taylor instability](https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Taylor_instability),
which occurs when a lighter fluid pushes into a heavier fluid.

# How to run

## Setting up the files

For the 2D version call the setup script like this:

    $AMRVAC_DIR/setup.pl -d=22 -g=20,20 -p=hd -eos=default

And for the 3D version:

    $AMRVAC_DIR/setup.pl -d=33 -g=20,20,20 -p=hd -eos=default

## Compiling

Simply issue the `make` command:

    make

When switching between 2D and 3D, perform a `make clean`.

## Running the code

To run with e.g. 4 processors, use

    mpirun -np 4 ./amrvac -i rt_2d.par

to run the code with the settings for 2D, and use

    mpirun -np 4 ./amrvac -i rt_3d.par

to run with the settings for 3D.

# Changing some of the settings

Some of the settings that you could change in the `.par` files are:

name | description
---|---
`filenameout` | Base file name for output
`dtsave(1)` | Time between log output
`dtsave(2)` | Time between dat/vtu output
`dtsave(5)` | Time between analysis output
`mxnest` | the maximum number of refinement levels
`dixB` | number of ghost cells, depends on spatial discretization
`nxlone[1,2,3]` | the size of the coarse grid, needs to be divisible by (block size - 2 * dixB)
`tmax` | the end time of the simulation
`typeadvance` | time discretization method, e.g., twostep, threestep, ssprk43
`courantpar` | CFL number (see also typecourant)
`typefull1` | spatial discretization method, e.g., tvd, tvdlf, hlcc
`typelimiter1` | which limiter to use in the spatial discretization, e.g., cada3, koren, woodward

A more complete list of parameters can be found [par.md](par.md).

# Changing the physics and initial conditions

Have a look at the local file `amrvacusr.t`. You can for example modify the
following routines:

* `initglobaldata_usr`: change the gravity or the adiabatic index \f$ \gamma \f$
* `initonegrid_usr`: change the initial conditions

# Looking at volume integrated quantities over time

The local file `analysis.t` shows how you can compute the kinetic energy in the
y-direction, making use of the function `get_kin_en_y` that computes \f$ \rho
v_y^2 \f$:

    pure function get_kin_en_y(w_vec, w_size) result(val)
      integer, intent(in)          :: w_size
      double precision, intent(in) :: w_vec(w_size)
      double precision             :: val

      val = w_vec(m2_)**2 / w_vec(rho_)
    end function get_kin_en_y

Here `m2_` refers to the momentum density in the second (y) dimension, and
`rho_` to the fluid density.

By changing `analysis.t` you can look at other quantities.
