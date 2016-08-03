# Setting up an amrvacusr module

This document describes how the **src/usr/amrvacusr.t.PROBLEM** and the
optional **src/usr/amrvacusrpar.t.PROBLEM** files should be written for user
defined initial and boundary conditions, input and output file formats, and
source terms. It also shows how the library source term routines
**src/amrvacmodules/*.t** can be included into the AMRVACUSR module. These
libraries (gravity, pointgrav, cooling, radloss, etc) are in principle self-
documented.

## Purpose and Use

The AMRVACUSR modules contain the problem dependent user written subroutines.
Usually a single AMRVACUSR module can be designed to contain several different
problems that all assume the same physics module. That would be realized by
the use of the parameter **iprob**, which is to be set in the corresponding
par-file. A _select case(iprob)_ construct can be used in appropriate places
then.

The setup is represented by two files in your simulation-directory,
**amrvacusr.t** and **amrvacusrpar.t** that can be copied from specific
templates in **src/usr/amrvacusr.t.PROBLEM** and
**src/usr/amrvacusrpar.t.PROBLEM** or from a test in the folder
**tests/EQUATION/PROBLEM/**. The first approach is automatized by running
**setup.pl**:

    $AMRVAC_DIR/setup.pl -u=PROBLEM

The **src/usr/amrvacusr.t.PROBLEM** file has to exist, but the
**src/usr/amrvacusrpar.t.PROBLEM** file is optional. If it does not exist, the
**amrvacusrpar.t** will be defaulted from **src/usr/amrvacusrpar.t.nul**.

_We however recommend adapting a suitable setup from the tests folder where
also the parameter file (default amrvac.par) and anything else to go with the
setup is present. _

## Creating a New Setup

In your designated simulation directory, start by copying the
**src/usr/amrvacusr.t.nul** file (or another similar file) into the new
**amrvacusr.t** file. It consists of a few include statements. The included
**amrvacnul/special*.t** files contain the default subroutines, and some or
all need to be specified for your problem. The arguments are declared and the
purpose of the subroutines is described below. Comment out the
**INCLUDE:amrvacnul/specialSUBROUTINE.t** statement(s) for the subroutine(s)
that you intend to write, and modify the comments at the beginning and at the
end of the module for clarity.

### Specialini part

Your start file should look something like this, where we already included the
**amrvacnul/specialini.t** file which is always needed:

    !=============================================================================
    ! amrvacusr.t.MYPROBLEM
    !=============================================================================
    !INCLUDE:amrvacnul/specialini.t
    INCLUDE:amrvacnul/speciallog.t
    INCLUDE:amrvacnul/specialbound.t
    INCLUDE:amrvacnul/specialsource.t
    INCLUDE:amrvacnul/usrflags.t
    !=============================================================================
    subroutine initglobaldata_usr

    include 'amrvacdef.f'
    !-----------------------------------------------------------------------------

    end subroutine initglobaldata_usr
    !=============================================================================
    subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid within ix^L

    include 'amrvacdef.f'

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    !-----------------------------------------------------------------------------

    w(ix^S,1:nw)=zero

    end subroutine initonegrid_usr
    !=============================================================================
    ! amrvacusr.t.MYPROBLEM
    !=============================================================================

Now you should edit both subroutines according to your needs: the idea is that
in _initglobaldata_usr_ you must set the global equation parameter values
(i.e. all _eqpar(*)_ entries, note again that you can code up different cases
depending on the **iprob** parameter). In the subroutine _initonegrid_usr_,
you have to make sure that at the end of this subroutine, all conservative
variable values are provided on the full grid, i.e. you need to specify
physically meaningfull _w_ entries. You have the grid info available in the
_x_ variable.

You can write the subroutine(s) either in the dimension independent notation,
described in **AMRVAC_Man/[source](source.html)**, or in Fortran 90 if the
number of dimensions is fixed for your PROBLEM.

Below some help is provided for writing new subroutines.

An example taken from the available _tests/rho/vac/amrvacusr.t_ user module is
given below

    !=============================================================================
    ! amrvacusr.t.testrho

    ! INCLUDE:amrvacnul/specialini.t
    INCLUDE:amrvacnul/speciallog.t
    INCLUDE:amrvacnul/specialbound.t
    INCLUDE:amrvacnul/specialsource.t
    INCLUDE:amrvacnul/usrflags.t
    !=============================================================================
    subroutine initglobaldata_usr

    include 'amrvacdef.f'
    !----------------------------------------------------------------------------

    {^IFONED   eqpar(v1_)=one }
    {^IFTWOD   eqpar(v1_)=one; eqpar(v2_)=one }
    {^IFTHREED eqpar(v1_)=one; eqpar(v2_)=one; eqpar(v3_)=one }

    end subroutine initglobaldata_usr
    !=============================================================================
    subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid

    include 'amrvacdef.f'

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision:: rhoflat,rhosquare,slocx^D
    double precision :: radius, xcircle^D
    !----------------------------------------------------------------------------

    rhoflat  = 0.5d0
    rhosquare= 2.0d0
    ! iprob=1 is a pure 1D Riemann problem, solvable in 1D, 2D, 3D
    if (iprob==1) then
        slocx^D=0.2d0;
        where({^D&x;(ix^S,^D)<=slocx^D|.and.})
           w(ix^S,rho_)     = rhosquare
        elsewhere
           w(ix^S,rho_)     = rhoflat
       endwhere

    ! **** many more cases in the actual file are omitted here ***

    else if (iprob==6) then
       radius = 0.2d0
       xcircle^D=zero;
       where(radius**2> ^D&(x(ix^S,^D)-xcircle^D)**2+ )
          w(ix^S,rho_)     = rhosquare
       elsewhere
          w(ix^S,rho_)     = rhoflat
       endwhere
    else
        call mpistop("iprob not available!")
    end if

    end subroutine initonegrid_usr
    !=============================================================================
    ! amrvacusr.t.testrho
    !=============================================================================

Note the use of the rho_ index name. It is clear that the **x** coordinates
are known on entry. The subroutine above works in 1, 2 or 3D.

### Specialbound part

When the predefined boundary types provided by MPI-AMRVAC are not sufficient
the **specialbound** subroutine can solve the problem. It is called for each
boundary region and variable for which the boundary type
**typeB='special'...** is selected. It is likely that you can use the
predefined boundary types for at least some of the variables and/or regions.
The **specialbound** subroutine is called _after_ the variables with
predefined boundary types have been updated in the ghost cells. The calling
interface for the specialbound subroutine is found in the
_amrvacnul.specialbound.t_ file.

An example for the use of this _specialbound_usr_ subroutine is found in the
example user file **usr/amrvacusr.t.wchd22**, which realizes the standard 2D
hydro Woodward and Collela shock reflection problem.

### Specialsource part

There are lots of possible physical source terms for the same basic equation.
Rather than writing a new physics module for each, it is simpler to define a
problem dependent source term in the AMRVACUSR module. There are already a
number of implemented source terms in the  src/amrvacusr.LIBRARY.t files,
which can be studied for example cases.

The **specialsource** subroutine is called at least once in every time step by
MPI-AMRVAC. The number of calls depends on the time integration scheme defined
by **typeadvance**, the parameter **sourcesplit** and also on the number of
dimensions if the parameter **dimsplit=T**.

In any case the subroutine should integrate **w** from time **qt** to
**qt+qdt** for the variables listed in **iw^LIM** in the region **ixO^S**. The
source terms should be evaluated for the **wCT** array, which corresponds to
the physical time **qtC**. In case of explicit time dependence, **qtC** should
be used as time. Only elements within **ixI^S** can be used from **wCT**.

The **getdt_special** subroutine can limit **dt** for numerical stability if
the source term requires that. This subroutine is called after the CFL
condition **getdt_courant** and the **getdt** subroutine of the AMRVACPHYS
module have been executed.

The **specialeta** subroutine is used for MHD to set the resistivity array
**eta** when it is not constant in time and/or space. The **current** array
must then be computed, when anomalous resistivity is described as a function
of **J**.

The **specialrefine_grid** subroutine allows to add user controlled
(de)refinement, by setting the integers _refine,coarsen_. You have all info
available to do this refining (grid level, physical values in _w_, coordinates
in _x_, time in _qt_). Similarly, the **specialvarforerrest** subroutine
allows to compute a (local) new variable, to be stored in _var_, on which you
can then base refinement as well. This is true for the lohner error estimator
only.

### Speciallog part

The _amrvacnul/speciallog.t_ file contains additional subroutines more related
to special I/O requests. The default log-file may be altered, for which you
need to code up the _printlog_special_ subroutine. For parallel execution,
this invariably means the use of MPI constructs, so you should copy in the
default version from _amrio.t_ and then study it, and modify accordingly.

The _process_grid_usr_ is a subroutine which allows to compute auxiliary
variables which happen to be non-local (like div v), and are in no way used
for flux computations. As auxiliaries, they are also not advanced. This
functionality was added to allow for separate particle treatments using
stochastic differential equations, where the particle dynamics is only relying
on local compression values etc.

The _specialvar_output_ is extremely handy to compute variables from the
actually computed conserved variables, that can then be visualized directly.
It is only used in combination with the conversion subroutines. E.g., one may
here compute current density components using the actual code discretizations
for computing a curl, and then visualize those with any of the visualization
tools applicable. You then also need to specify a label for this variable, in
_specialvarnames_output_.

## AMRVACUSR Library

Various source terms are available as library subroutines, in particular for a
uniform external gravitational field, for an external gravitational point
mass, and for optically thin radiative losses. They will always need to be
combined with user written subroutines. To include a library into the
**amrvacusr.t** file, just add a line

    INCLUDE:amrvacmodules/LIBRARY.t

and call the appropriate library routines from the subroutines
**specialsource** and **getdt_special** according to the description of the
library file. An example of that for a constant external gravity is in the
problem file **src/usr/amrvacusr.t.testhdrt**, which includes the gravity
library. It is also possible to copy the libraries into **amrvacusr.t**
directly and modify them as necessary. The parameters of the library should be
defined in the **amrvacusrpar.t** file according to the description given in
the library file. See the [equations](equations.html#HD) description as well,
below we just list radiative loss treatments.

### Radiative losses: amrvacmodules/radloss.t and amrvacmodules/cooling.t

An optically thin gas cools due to radiative losses. This involves the energy
equation only:

![](figmovdir/eq.radloss.gif)

The thermal energy loss is proportional to density squared and a complicated
function of the temperature. The two libraries differ in the details of this
function, the more general _amrvacmodules/cooling.t_ has many frequently used
cooling tables implemented, and various ways to add this local source term.

## Special Equation Parameters

The user-defined source terms or boundary conditions may contain parameters
which are often changed and have similar meaning to the equation parameters
defined in the **src/EQUATION/amrvacpar.t** files. The **amrvacusrpar.t** file
allows the user to define extra, problem dependent, equation parameters.

The indexname and the number of the special equation parameters can be defined
in the **amrvacusrpar.t** file. The values of these parameters should be set
in the _initglobaldata_usr_ subroutine.

To prepare a new **amrvacusrpar.t** file, simply copy the
**src/usr/amrvacusrpar.t.nul** file into **amrvacusrpar.t** and edit it. This
file will be included into the variable declaration part of all subroutines,
thus it can also be used for variables to be shared by subroutines in the
AMRVACUSR module.

A simple example is the following file, taken from
_src/usr/amrvacusrpar.t.testhdrt_ which just says the code that it has
equation parameters for the constant gravitational field.

    !##############################################################################
    ! include amrvacusrpar - gravity

    INTEGER,PARAMETER:: grav0_=neqpar, grav^D_=neqpar+^D, nspecialpar=^ND
    {^IFONED   CHARACTER*5 ,PARAMETER:: specialparname='grav1'}
    {^IFTWOD   CHARACTER*11,PARAMETER:: specialparname='grav1 grav2'}
    {^IFTHREED CHARACTER*17,PARAMETER:: specialparname='grav1 grav2 grav3'}

    ! end include amrvacusrpar - gravity
    !##############################################################################

