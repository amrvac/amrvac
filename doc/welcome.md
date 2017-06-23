# Welcome page

[TOC]

# Introduction {#introduction}

This is the documentation for the development version of MPI-AMRVAC,
which is now hosted on [Gitlab](https://gitlab.com/mpi-amrvac/amrvac).

# Quick links {#quick_links}

* @ref installation.md
* @ref getting_started.md
* @ref changelog.md
* @ref faq.md
* @ref contributing.md
* [List of tests](@ref test)
* @ref publications.md
* @ref acknowledgments.md
* @ref todo

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

# History {#history}

This version of the MPI-AMRVAC software was initiated in the course of 2006-2007
by **Bart van der Holst** (meanwhile at University of Michigan) in a close
collaboration with [Rony Keppens](http://perswww.kuleuven.be/~u0016541) at the
[Centre for Plasma-Astrophysics (CPA)](http://wis.kuleuven.be/cpa) at
K.U.Leuven. The code has witnessed continuous improvements and additions, with
**Zakaria Meliani** as the main core developer since 2005-2006 (for the
relativistic physics modules, parallel conversion of MPI-AMRVAC I/O to a large
variety of postprocessing formats for later visualization, and for many overall
code additions/improvements), and **Oliver Porth** joining in on the core team
since 2010. Other contributors to the project are **Allard Jan van Marle** (for
optically thin radiative loss treatments), **Peter Delmont** (some of the
visualization capabilities), **Chun Xia** ([an]isotropic thermal conduction) and
the growing number of users at CmPA and abroad.

MPI-AMRVAC is an MPI-parallelized Adaptive Mesh Refinement code, with some
heritage (in the solver part) to the [Versatile Advection Code
](http://grid.engin.umich.edu/~gtoth/VAC) or VAC, initiated by Gabor Toth at
the Astronomical Institute at Utrecht in November 1994, with help from Rony
Keppens since 1996. Previous incarnations of the Adaptive Mesh Refinement
version of VAC were of restricted use only, and have been used for basic
research in AMR strategies, or for well-targeted applications. This MPI
version uses a full octree block-based approach, and allows for general
orthogonal coordinate systems. Tests have been performed on various
supercomputing facilities throughout Europe.

# Philosophy {#philosophy}

The philosophy behind MPI-AMRVAC is using a single versatile software with
options and switches for various problems, rather than developing a different
method or version for each problem separately. The advantage of such a general
approach is a reduction of overall time for software development, easier
maintenance, compatibility of different parts, automatic extension of new
features to all existing applications. The price of the general approach is
some added complexity in the source.

MPI-AMRVAC is not a fool-proof black-box design. The user is expected to
understand how the different parameters change the behaviour of the code, and
to be able to complete user written subroutines for source terms, special
boundary conditions etc. if needed. It is essential to keep the user written
parts well separated. Please put added subroutines into a separate AMRVACUSR
module in the **src/usr** subdirectory under the **src** (source) directory.
You must then think of a suitable name, of the form
**src/usr/amrvacusr.t.MINE**.
