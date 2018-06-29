Comments on modifications
==================

**This package was not written by me and all credit goes to the original authors
(L. Genovese et al.).**

I needed a free-space Poisson solver in 3D, and have downloaded this
GPL-licensed solver from the [Computational Physics
website](https://comphys.unibas.ch/software.htm) of the University of Basel,
through this [link](https://comphys.unibas.ch/SOFTWARE/FELORAPO/solver.tgz). The
solver is part of the [BigDFT code](http://bigdft.org).

I've made a couple of small modifications, which are perhaps useful for others:

* Simplify the makefile
* Always require MPI for building
* Build a static library, so that the solver can easily be embedded
* By default, disable printed output (so the solver stays quiet)
* Add usage instructions to PSolver example

Requirements: gfortran (or another Fortran compiler) and MPI.

The original README is copied below:

README for PSolver
==================



<pre>
	Copyright (C) 2006 BigDFT Group

	This file is part of BigDFT.

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2, or (at your option)
	any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; see the file COPYING. If not, write to
	the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
	Boston, MA 02111-1307, USA.
</pre>

This Poisson Solver is part of the BigDFT project but can also be used separately.
All the source files of this directory

ABINIT-common/defs_basis.F90 
ABINIT-common/defs_datatypes.F90 
ABINIT-common/defs_xc.F90 
ABINIT-common/drivexc.F90 
ABINIT-common/xctetr.F90 
ABINIT-common/xcwign.F90 
ABINIT-common/xcxalp.F90 
ABINIT-common/xchelu.F90 
ABINIT-common/xcpzca.F90 
ABINIT-common/xcspol.F90 
ABINIT-common/xchcth.F90 
ABINIT-common/xclb.F90 
ABINIT-common/xcpbe.F90 
ABINIT-common/invcb.F90
Poisson_Solver.f90 
PSolver.f90 
timing.f90
dcopy.f
perfdata.inc 
lazy_100.inc 
lazy_16.inc 
lazy_24.inc 
lazy_40.inc 
lazy_60.inc 
lazy_14.inc 
lazy_20.inc 
lazy_30.inc 
lazy_50.inc 
lazy_8.inc
PSolver_Main.f90 
Build_Kernel.f90 
PSolver_Base.f90 
xcenergy.f90 
3Dgradient.f90 
fft3d.f90 
scaling_function.f90
Makefile
bench.csh 
perfdata.t 

are delivered under the GNU General Public License; please see the file COPYING 
for copying conditions.

See the file INSTALL for generic compilation and installation instructions.

This Poisson Solver code is built following the method described in the papers

 Genovese, L; Deutsch, T; Goedecker, S.
 Journal of Chemical Physics, 127 (5), 054704 (2007).
 Efficient and accurate three-dimensional Poisson solver for surface
 problems.
 http://dx.doi.org/10.1063/1.2754685

 Genovese, L; Deutsch, T; Neelov, A; Goedecker, S; Beylkin, G.
 Journal of Chemical Physics, 125 (7), 074105 (2006).
 Efficient solution of Poisson's equation with free boundary conditions
 http://dx.doi.org/10.1063/1.2335442

Citing of this references is greatly appreciated if the routines are used 
for scientific work.
