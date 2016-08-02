# Overview of MPI-AMRVAC features

This is a brief overview of the main features of the MPI-AMRVAC software
package.

## Structure

This software has to be configured, preprocessed, and compiled into a single
main executable program, **amrvac** which can be run on multiple processors
using MPI. MPI-AMRVAC will initialize and advance the unknowns, and perform
automated grid refinement to follow all details of the (possibly shock-
dominated) flow. The program is split into several logical parts, and heavily
uses Fortran 90 modules and dynamic allocation. The various parts are simply
sets of subroutines and functions that belong together and they are put in a
single file.

    partname		    Function

    amrvac			    Main program and subroutines.
    amrvacio/    		Input and output.
    amrvacnul/    		Dummies
    amrvacmodules/      AMRVACUSR LIBRARY functions
    amrini 	     		Initialization.
    amrgrid     		Grid setup.
    advance 			Advancing the grid tree and scheme/method selection.
    geometry			Initializing the appropriate geometry information.
    amr_coarsen_refine, amr_fix_conserve, amr_ghostcells, amr_neighbours, amr_solution_node
                         AMR specific routines
    coarsen			    grid coarsening operations
    comm_lib			MPI communication routines
    connectivity		subroutines for neighbor determinations in the tree hierarchy
    errest			    the error estimator used for AMR regridding
    forest			    the AMR octree 
    initialize_vars		basic global variable/parameter initializations
    load_balance		the load balance strategy
    refine			    refining
    set_B0     			subtracting a background potential field in MHD simulations
    setdt 			    Determining the time step limit
    convert			    Conversion of data files to different formats for visualization (postprocessing)
    tvdlf 			    Total Variation Diminishing Lax-Friedrich and Hancock predictor method, as well as HLL(C) variants
    tvd 				Total Variation Diminishing method
    cd	     			Central difference scheme
    fd	     			Finite difference scheme
    EQUATION/amrvacpar.t  	Basic equation dependent parameters.
    EQUATION/amrvacphys.t 	Basic equation dependent subroutines.
    EQUATION/roe.t 		    Basic equation dependent subroutines related to Roe scheme.
    EQUATION/hllc.t 		Basic equation dependent subroutines related to HLLC scheme.
    EQUATION/correctaux.t 	Equation dependent subroutines for error handling.
    usr/amrvacusr.t.PROBLEM Problem dependent user written subroutines.
    modules/                shared module files

The AMRVACPHYS and AMRVACUSR modules have several versions, but only the
actual module, selected by the Perl script **setup.pl** is compiled at a time.
See AMRVAC_Man/[USAGE](USAGE.html) about the configuration procedure.

Once the configuration and compilation are done, **amrvac** can advance the
solution in time, or can be used to convert previously produced data files to
specific other formats for visualization. The data itself will be saved into
_*.dat_ data files, which each contain a single snapshot of all grids and
their unknowns in time, which can be used for eventual restarts. Hence, a
simulation can be continued from any saved time step. Once a snapshot is
available in _*.dat_ format, the result can be converted for further
visualization, e.g. to _*.vtu_ format to be visualized with Paraview, see
further info on [conversion](convert.html).

## Source Language and Compilation

The `*.t` source files of AMRVAC are written in [dimension independent
notation](source.html), known as the _LASY_ syntax. A suitably modified (and
also simplified) version of the VAC Pre-Processor, [VACPP](vacpp.html)
translates the source code files to Fortran 90. The code is to be run on
**parallel machines** using MPI, so even on a single processor laptop or
desktop, we still require compilation for MPI.

## Equations

In general, MPI-AMRVAC aims to solve a system of (primarily hyperbolic)
partial differential equations written in conservation form, with optional
source terms on the right hand side:

![](figmovdir/eq.general.gif)

All the equation dependent subroutines are collected in the corresponding
subdirectory (physics module), e.g. **src/EQUATION/amrvacphys.t**. Based on
the existing physics module, one can easily adjust the code to allow for a new
(set of) equation(s). The following equations are already defined:

Module name	| Equation
---|---
rho	| Transport equation [d(rho)/dt+div(v*rho)=0, v fixed]
nonlinear |	Nonlinear scalar, e.g. Burgers
hd | Hydrodynamics [Euler equations], ideal gas _-eos=gamma (default)_, isentropic [p=const*rho^gamma] or isothermal eos [p=const*rho] _-eos=iso_
mhd	| Magnetohydrodynamics, ideal gas _-eos=gamma (default)_, isentropic [p=const*rho^gamma] or isothermal eos [p=const*rho] _-eos=iso_, also zero-beta plasma
srhd | special relativistic hydrodynamics, ideal gas EOS
srhdeos	| special relativistic hydrodynamics, Synge-type TM EOS
srmhd | special relativistic magnetohydrodynamics, ideal gas _-eos=gamma (default)_, Synge-type TM EOS _-eos=synge_ and isentropic _-eos=iso_

The resistivity for mhd can be a function of the flow variables as well as
position. Some typical source terms have been implemented as a AMRVACUSR
LIBRARY, collected in the **src/amrvacmodules/** folder:



Library | Equation(s)  Source terms
---|---
radloss | hd, mhd  optically thin radiative losses
gravity | hd, mhd  external gravity, constant magnitude and direction
pointgrav | hd, mhd  external gravity for (several) point sources
heatconduct | hd, mhd  (an)isotropic thermal conduction
viscosity | hd, mhd  viscosity
epsinf | srmhd    synchrotron losses of the cutoff electron energy
raytracing | all      solve equations on a raygrid, coupling in both directions with the fluid
fff | mhd      linear force free field extrapolation
pfss | mhd      potential field spherical surface extrapolation

These can be modified and included into the AMRVACUSR module.

See AMRVAC_Man/[equations](equations.html) for more detail.

## Grid and Boundary

**The base grid is a 1, 2, or 3 dimensional curvilinear grid, with slab or axial symmetry for the ignored direction in less than 3D.** The grid can be Cartesian, cylindrical, spherical in 3D, or a polar grid in 2D. The vector variables are represented by their respective cartesian, cylindrical, or spherical components. Numerical conservation is ensured over the grid hierarchy by a finite volume discretization.

The **slab**, **cylindrical** or **spherical** grids differ in the definition
of volumes and surfaces for the grid cells, and there are some extra terms in
the equations, like the p/r term in the radial momentum equation for
hydrodynamics. These are defined in the **addgeometry** subroutines for each
AMRVACPHYS module. For polar grids the same geometrical source terms are used.

The boundaries are represented by ghost cells around the physical mesh. Of
course, for a grid-adaptive computation, internal boundaries are handled
appropriately, and the user is expected to interfere only with the physical
domain boundaries, which represent 2, 4, or 6 regions depending on the
dimensionality (1,2 or 3). The boundary types are defined for each region and
each variable separately. The following boundary types are available:

Type | Value of the ghost cell
---|---
periodic | Copied from opposite edge of the mesh
symm | Reflected from closeby mesh cells
asymm | Reflected from closeby cells and multiplied by -1
cont | Copied from mesh cell next to the ghost cell
special | Defined by the specialbound subroutine in AMRVACUSR module

See AMRVAC_Man/[discretization](discretization.html#Grid) and
AMRVAC_Man/[axial](axial.html) for further information.

## Spatial Discretization

**In MPI-AMRVAC, most discretizations are shock capturing conservative numerical schemes**. A non exhaustive list is given by:

Name | Description
---|---
cd | Central difference
tvdlf | MUSCL type TVD-Lax-Friedrich
tvdlf1 | First order TVD-Lax-Friedrich
tvdmu | TVD-MUSCL with Hancock predictor
tvdmu1 | First order upwind scheme
tvd | One-step temporally 2nd order TVD
tvd1 | One-step temporally 1st order TVD
fd | conservative finite differences, up to fifth order with MP5 reconstruction

In multidimensional MHD calculations the divergence of the magnetic field may
become significantly different from zero. This may cause numerical instability
or inaccurate results. There are several source-term options to fix this
problem. E.g., **Powell's** non-conservative source terms, which are
proportional to **div B**, can be used to stabilize, and to improve the
accuracy for any of the methods. We have also provide several variants of
**Dedner's** GLM scheme.

See AMRVAC_Man/[methods](methods.html) for a more detailed description.

## Temporal Discretization

With the finite volume methods, second order accurate explicit time
integration can be obtained by predictor-corrector and multistep Runge-Kutta
type discretizations. In finite differences, the overall accuracy can go up to
fourth order using the appropriate temporal discretizations. The more
important ones are:

Name | Description
---|---
onestep | 1st order Euler scheme
twostep | 2nd order predictor-corrector schemes
threestep | 3rd order (TVD) Runge-Kutta
rk4 | 4th order Runge-Kutta (classical)
ssprk43 | 3rd order, strong stability preserving four step Runge-Kutta
ssprk54 | 4th order, strong stability preserving five step Runge-Kutta

The time step can be adjusted dynamically to satisfy the stability criteria.

## Plans for Future Versions

There are several planned extensions. Some of the more important and clearly
defined directions are listed below.

**Adaptive Mesh and Algorithm Refinement (AMAR)**      MPI-AMAR-VAC will combine the flexibility of MPI-AMRVAC with the possibility to couple different physics levels across the grid hierarchy.
