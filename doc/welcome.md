# Welcome page

[TOC]

# Introduction {#introduction}

This is the documentation for the development version of MPI-AMRVAC,
which is now hosted on [Gitlab](https://gitlab.com/mpi-amrvac/amrvac).

# Quick links {#quick_links}

* [Contents of the documentation](contents.md)
* [Getting started](getting_started.md)
* [List of tests](@ref test)
* [FAQ](faq.md)
* [TODO list](@ref todo)
* [Change log](changelog.md)
* [Contributing](contributing.md)
* [Acknowledgements](acknowledgments.md)

# Project {#project}

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

# Code aims and requirements {#aims}

MPI-AMRVAC aims to advance any system of (primarily hyperbolic) partial
differential equations by a number of different numerical schemes. The
emphasis is on (near) conservation laws, with shock-dominated problems as a
main research target. The actual equations are stored in separate modules, can
be added if needed, and they can be selected by a simple configuration of the
VACPP preprocessor. The dimensionality of the problem is also set through
VACPP. The numerical schemes are able to handle discontinuities and smooth
flows as well.

For running and installing the code, you need MPI (Message Passing Interface,
[MPI-2](http://www.mpi-forum.org/index.html)) combined with a Fortran 90 (or
95) compiler, and Perl for the pre-processing, which is standard on any Unix-
based platform. To make sure you have access to the latest version of the
code, you will need [subversion](http://subversion.apache.org/), a free, open-
source version control system. No other libraries are required. The code is in
use on Mac laptops, Unix desktops, to modern supercomputing facilities.

Visualization or analysis of the results can be done **a posteriori**, by
converting data files produced during runtime. Each data file contains a
single snapshot of the entire data structure (all unknowns and the full grid
tree), and are essential for e.g. restarts, and need to be saved in the
internal `*.dat` format. Each snapshot can be converted **with the same
executable used to run the simulation (but possibly on a different platform)**
to IDL, openDX, Tecplot, or VTK native formats, which must then be visualized
with the corresponding software. We (June 2010) mostly use the freely
available [Paraview](http://www.paraview.org/) application for data
visualization.

The [VACPP](vacpp.md) preprocessor translates the dimensional independent
notation (based on the Loop Annotation Syntax, or LASY) to Fortran 90. It uses
Perl, and the Perl interpreter is installed on any UNIX-like platform at most
places, and it is freely available from the net.

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

# Selected Related Publications {#publications}

**The following publications contain relevant information, partly applicable to
the current MPI-version MPI-AMRVAC (some, especially the older papers used
previous variants of the code). Especially the first paper listed here is
specific to the current version.**

  1. `Parallel, grid-adaptive approaches for relativistic hydro and magnetohydrodynamics', R. Keppens, Z. Meliani, A.J. van Marle, P. Delmont, A. Vlasis, &amp; B. van der Holst, 2011, JCP. [doi:10.1016/j.jcp.2011.01.020](http://dx.doi.org/10.1016/j.jcp.2011.01.020).

  2. `A multidimensional grid-adaptive relativistic magnetofluid code', B. van der Holst, R. Keppens &amp; Z. Meliani, 2008, Comp. Phys. Commun. ** 179**, 617-627

  3. `Hybrid block-AMR in cartesian and curvilinear coordinates: MHD applications', B. van der Holst &amp; R. Keppens, 2007, JCP **226**, 925-946

  4. `AMRVAC and relativistic hydrodynamic simulations for GRB afterglow phases', Z. Meliani, R. Keppens, F. Casse, &amp; D. Giannios, 2007, MNRAS **376**, 1189-1200

  5. `Adaptive Mesh Refinement for conservative systems: multi-dimensional efficiency evaluation', R. Keppens, M. Nool, G. Toth, J.P. Goedbloed, 2003, Comp. Phys. Comm. **153** (No. 3, 1 july issue), 317-339.

**The following publications represent applications. Some relate to earlier versions of AMRVAC.**

  6. `On radiative cooling in numerical astrophysics: the need for adaptive mesh refinement', A.J. van Marle &amp; R. Keppens, 2011, Computers &amp; Fluids **42**, 44-53

  7. `Thin shell morphology in the circumstellar medium of massive binaries', A.J. van Marle, R. Keppens, &amp; Z. Meliani, 2011, A &amp; A **527**, A3 (DOI: 10.1051/0004-6361/201015517)

  8. `Decelerating relativistic two-component jets', Z. Meliani &amp; R. Keppens, 2009, ApJ **705**, 1594-1606

  9. `Extragalactic jets with helical magnetic fields: relativistic MHD simulations', R. Keppens, Z. Meliani, B. van der Holst &amp; F. Casse, 2008, Astron. &amp; Astrophys. **486**, 663-678

  10. `Fanaroff-Riley type I jet deceleration at density discontinuities. Relativistic hydrodynamics with realistic equation of state', Z. Meliani, R. Keppens &amp; B. Giacomazzo, 2008, Astron. &amp; Astrophys. **491**, 321-337

  11. `GRB blastwaves through wind-shaped circumburst media', Z. Meliani &amp; R. Keppens, 2007, Astron. &amp; Astrophys. **467**, L41-L44

  12. `Kelvin-Helmholtz disruptions in extended magnetized jet flows', H. Baty &amp; R. Keppens, 2006, Astron. &amp; Astrophys. **447**, 9-22

  13. `Simulations of Relativistic Astrophysical Flows', J. Bergmans, R. Keppens, D.E.A. van Odyck, &amp; A. Achterberg, 2005, in ``Adaptive Mesh Refinement -- Theory and Applications'', Eds. T. Plewa, T. Linde &amp; V.G. Weirs, Lecture Notes in Computational Science and Engineering Vol. **41**, Springer, Berlin Heidelberg New York (ISBN 3-540-21147-0, Proceedings of the Chicago Workshop on Adaptive Mesh Refinement Methods, September 3-5 2003, Chicago, USA), 223-234.

  14. `Grid-adaptive computations of magnetized jets', R. Keppens, H. Baty &amp; F. Casse, 2005, Space Science Reviews ** 121**, 65-75

  15. `The 2D MHD Kelvin-Helmholtz instability: compressibility and large-scale coalescence effects', H. Baty, R. Keppens, &amp; P. Comte, 2003, Phys. of Plasmas **10**(12), 4661-4674.

  16. `AMRVAC: A multidimensional grid-adaptive magnetofluid dynamics code', M. Nool &amp; R. Keppens, 2002, Comp. Meth. Appl. Math. **2**(No.1), 92.
