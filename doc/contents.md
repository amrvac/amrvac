# Contents {#doc-contents}

[TOC]

# Introduction for new users {#intro_new_users}

* [Features](features.md) An overview of the main features of the Message
Passing Interface - Adaptive Mesh Refinement - Versatile Advection Code.
* [Getting Started](getting_started.md) How to install the Message Passing
Interface - Adaptive Mesh Refinement - Versatile Advection Code and run your
first test problem.
* [Acknowledgments](acknowledgments.md) Information on collaboration and
financial support.
* [Changes in the GIT version](gitversion.md) How to migrate from the old svn
code to the new git version.
* [FAQ](faq.md) Frequently asked questions.

# General {#general}

* [Command line](commandline.md)      Help on command-line parameters.
* [Examples](examples.md) Description of various example simulations for which
parameter files and AMRVACUSR modules have been provided.
* [Equations](equations.md) The equations and parameters in the existing
AMRVACPHYS modules. How to create a new AMRVACPHYS module.
* [User files](amrvacusr.md) How to create a new problem, specify initial
  conditions and customize AMRVAC.

# Discretization and AMR related matter {#discretization}

* [Discretization](discretization.md) The equation and its discretization, the
basic variables, the structure of the grid, boundaries, ghost cells.
* [Using polar/cylindrical/spherical coordinates](axial.md) Some information on
simulations using non-Cartesian grids, and notes on axial versus translational
symmetry.
* [Methods](methods.md) Properties of the discretization methods like TVDLF,
TVDMU, TVD, HLL, ...
* [AMR aspects](amrstructure.md) Some essential info on global parameters and
the data structures for the block-tree AMR.

# Parameters {#parameters}

* [Parameters for AMRVAC](par.md) Description of the "amrvac.par" parameter file
for MPI-AMRVAC.
* [Auxiliary variables for MPI-AMRVAC](mpiamrvac_nw.md) Description of the
intended use for _nw, nwflux, nwaux, nwextra, nwauxio_ parameters.

# Special Sources {#special_sources}

* [Thermal conduction](thermal_conduction.md) Description of solving thermal conduction. 
* [Radiative cooling](radiative_cooling.md) Description of adding radiative cooling. 

# Source Code {#source_code}

* [Source](source.md) Description of the dimension independent source language,
which is translated to F90 by the VACPP preprocessor.
* [Variable Names](varnames.md) How variable names are formed in the source
files.
* [VACPP](vacpp.md) Making and running the VACPP preprocessor itself.

# IO and postprocessing {#io}

* [File format](fileformat.md) Description of the MPI-AMRVAC data file format
***.dat**.
* [Slices](slices.md) How to output hypersurfaces (slices) for restart or
visualization.
* [Line of sight views](collapsed.md) How to output collapsed views for
visualisation and analysis (e.g. column densities)
* [Analysis routine](analysis.md) Using the run-time analysis routine
* [Converting data files for visualization](convert.md) Brief notes on how to
convert to IDL, DX (openDX), Native Tecplot (*.plt), and VTK (*.vtu) data files
(the latter for paraview visualization).

**Postprocessing tools**

* [IDL](idl.md)
* [Python](python/index.md) \todo Include python documentation here
