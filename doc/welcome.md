# Welcome page

[TOC]

# Introduction {#introduction}

This is the documentation for the development version of MPI-AMRVAC. The code is
available on [Github](https://github.com/amrvac/amrvac), and the documentation
on [amrvac.org](http://amrvac.org/). If you have questions about MPI-AMRVAC,
please send them to the mailing list: <mailto:amrvacusers@ls.kuleuven.be>, which you can also
[subscribe to](https://ls.kuleuven.be/cgi-bin/wa?SUBED1=AMRVACUSERS&A=1) and
[search](https://ls.kuleuven.be/cgi-bin/wa?A0=AMRVACUSERS).

# Quick links {#quick_links}

* @ref installation.md
* @ref getting_started.md
* @ref par.md
* @ref faq.md
* @ref contributing.md
* @ref publications.md
* @ref acknowledgments.md
* [**Recent changes**](https://github.com/amrvac/amrvac/commits/master)

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

# Development {#current_develop}

MPI-AMRVAC is developed and maintained by an international team led by professor 
[Rony Keppens](https://perswww.kuleuven.be/~u0016541/) from 
[Centre for mathematical Plasma-Astrophysics (CmPA)](https://wis.kuleuven.be/CmPA), KU Leuven.

In November 2022, we released the current 3.0 version, with the aid of Beatrice Popescu Braileanu, Niels Claes, Chun Xia, Guo Yang, Wenzhi Ruan, Fabio Bacchini, Yuhao Zhou.
The movies generated from the demo simulations (found in the tests/demo folder) can be seen [here](demo-movies.md). 

Prior, in 2016-2017, a large modernization was completed by Chun Xia and 
[Jannis Teunissen](http://teunissen.net/) marked with version 2.0 with 
important changes as following:

* Automatic regression tests were added
* The preprocessor is only used for the problem dimension (1D, 2D, or 3D)
* The code was modernized and re-organized into Fortran modules, and there is
  now an AMRVAC library
* The focus of MPI-AMRVAC is now on non-relativistic hydro- and
  magnetohydrodynamics
* Many smaller improvements to the physics modules

From 2018, more developers have joined in to contribute in various aspects. 
We explicitly mention important help from Oliver Porth (main developer of MPI-AMRVAC 1.0), 
 and Hector Olivares, who are the developers of the sister GRMHD code [BHAC](http://bhac.science/), the Black Hole Accretion Code.

