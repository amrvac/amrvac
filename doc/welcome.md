# Welcome page

[TOC]

# Introduction {#introduction}

This is the documentation for the development version of MPI-AMRVAC. The code is
hosted on [Gitlab](https://gitlab.com/mpi-amrvac/amrvac), and the documentation
on [amrvac.org](http://amrvac.org/).

# Quick links {#quick_links}

* @ref installation.md
* @ref getting_started.md
* @ref contributing.md
* [List of tests](@ref test)
* @ref publications.md
* @ref acknowledgments.md
* @ref todo
* [**Recent changes**](https://gitlab.com/mpi-amrvac/amrvac/commits/physics_modules)

# MPI-AMRVAC aims {#aims}

MPI-AMRVAC is a parallel adaptive mesh refinement framework aimed at solving
(primarily hyperbolic) partial differential equations by a number of different
numerical schemes. The emphasis is on (near) conservation laws and on
shock-dominated problems in particular. A number of *physics modules* are
included; the hydrodynamics and the magnetohydrodynamics module are most
frequently used. Users can add their own physics module or modify existing ones.
The framework supports 1D to 3D simulations, in a number of different geometries
(Cartesian, cylindrical, spherical).

MPI-AMRVAC is written in Fortran 90 and uses MPI for parallelization.
The [VACPP](vacpp.md) preprocessor is used to extend Fortran with dimensional
independent notation, but users are not required to learn the VACPP syntax.

The philosophy behind MPI-AMRVAC is to use a single versatile code with options
and switches for various problems. The advantage of such a general approach is
easier maintenance, the compatibility of different parts, and the automatic
extension of new features to existing applications. MPI-AMRVAC is not a
fool-proof black-box design. A user needs to write subroutines for initial
conditions, and for source terms or special boundary conditions when needed.

# Current development {#current_develop}

MPI-ARMVAC is currently maintained and further developed by members of
the
[Centre for mathematical Plasma-Astrophysics (CmPA)](https://wis.kuleuven.be/CmPA).

In 2016-2017 a large modernization effort was started
by [Chun Xia](https://wis.kuleuven.be/CmPA/current-members/00086290)
and [Jannis Teunissen](http://teunissen.net/), who are the current maintainers.
The changes in the new version are:

* Automatic regression tests were added
* The preprocessor is only used for the problem dimension (1D, 2D, or 3D)
* The code was modernized and re-organized into Fortran modules, and there is
  now an AMRVAC library
* The focus of MPI-AMRVAC is now on non-relativistic hydro- and
  magnetohydrodynamics
* Many smaller improvements the physics modules
