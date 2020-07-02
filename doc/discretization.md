# Discretization

[TOC]

# Introduction {#disc-intro}

This document briefly describes how the equations are discretized.

# Equations {#disc-eq}

MPI-AMRVAC aims to solve a system of partial differential equations of the
form

    dw/dt=-dF_i(w)/dx_i+S(w)=R(w)

where w is a vector of **iw=1..nw** flow variables, **F_i** are fluxes in
**idim=1..ndim** directions for each **w**, and **S** is the source including
all terms which are not described by **F_i**. The flux and source terms
together are denoted by **R**. The discretization is of a finite volume kind, 
where cell-centered quantities represent cell-averaged densities for 
conservative variables.

# Dimensional and Source Splitting {#disc-splitting}

Some methods were originally designed to handle 1D problems, they may be best
applied to multidimensional problems by splitting the equation for the
derivatives of the flux, and solving the split equations in subsequent sweeps.
All these strategies are handled in the _mod_advance.t_ module.

For a forward Euler time integration scheme the above equation can be
discretised in time as

    w^(n+1)    = [I+dt*(Fx+Fy+Fz+S)] w^n

where Fx, Fy, and Fz are short notations for the dF_i/dx_i terms for 3D.

Dimensional splitting means the following approximations:

    I+dt*(Fx+Fy+S)    --> [I+dt*(Fx+S/2)] [I+dt*(Fy+S/2)]                   (2D)

    I+dt*(Fx+Fy+Fz+S) --> [I+dt*(Fx+S/3)] [I+dt*(Fy+S/3)] [I+dt*(Fz+S/3)]   (3D)

In case **dimsplit=T**, the default splitting strategy is according to
**typedimsplit='xyyx'**, which ensures alternation in the order of the sweeps
in every second time step. This alternation can be switched off by
**typedimsplit='xy'** which may be useful for steady state calculations.
Depending on the selected scheme, we can also invoke no dimensional splitting
at all, hence **dimsplit=F** (which implies **typedimsplit='unsplit'**).

The source terms S may also be separated from the fluxes F. A reason to use
split source terms may arise, when the sources are stiff, or unstable for an
explicit evaluation. For steady state evolutions, one would typically use
unsplit sources and no dimensional splitting at all. Splitting user sources from
fluxes is achieved by **source_split_usr=T**. How sources are then added depends on
**typesourcesplit**, which can be any of 4 pre-implemented combinations.

Setting **typesourcesplit='sfs'** results in

    I+dt*(Fx+Fy+Fz+S) --> [I+dt/2*S] [I+dt*(Fx+Fy+Fz)] [I+dt/2*S]

for the dimensionally unsplit case, and

    I+dt*(Fx+Fy+S)    --> [I+dt/2*S] [I+dt*Fx] [I+dt*Fy] [I+dt/2*S]            (2D)

    I+dt*(Fx+Fy+Fz+S) --> [I+dt/2*S] [I+dt*Fx] [I+dt*Fy] [I+dt*Fz] [I+dt/2*S]  (3D)

for the dimensionally split case. To achieve second order time accuracy, the
numerical representation of the source term **S** should be second order
accurate in time. In case of special source terms written by the user, this
may be achieved by implementing in the pointed **usr_source** subroutine.

However a simple second order Runge-Kutta evaluation is already built in.
Setting **typesourcesplit='ssfss'** gives the following evaluation for **w_S =
[I+dt/2*S] w**:

    w_1 = w + dt/4 S(w)
    w_S = w + dt/2 S(w_1)

both at the beginning and at the end of the time step, otherwise the default
**w_S = w + dt/2 S(w)** is used. Other choices for **typesourcesplit** are
**'sf', 'ssf'**.

# Time Discretization {#disc-time}

The code can use a variety of [methods](methods.md) for spatial
discretization of fluxes. A method is applied to all variables on a specific
AMR level, but the method may differ from grid level to grid level. As we saw
in the previous section, the contribution of the fluxes (dimensionally split
or unsplit, with or without the source terms) need to be evaluated to at least
2nd order accuracy in time. Let **R** denote the general operator that
represents any of the Fx, (Fz+S/3), (Fx+Fy), (Fx+Fy+Fz+S), etc. terms that
arise in the equations above depending on the number of dimensions and on the
dimensional and source splitting parameters.

There are many options to evaluate **w_R=[I+dt*R] w**, it is determined by the
value of the **time_stepper** (and the **time_integrator**) parameter with the following examples:

    w_R = w + dt*R(w)                 'onestep'

    w_1 = w + dt/2*R1(w)               'twostep'  (predictor step)
    w_R = w + dt  *R(w_1)                        (corrector step)

    w_1 = w + dt/2*R(w)	              'fourstep' (classical RK4 when `time_integrator='rk4'`)
    w_2 = w + dt/2*R(w_1)
    w_3 = w + dt  *R(w_2)
    w_4 = w + dt/6*R(w_3)
    w   = w + dt/6*[R(w_1)+2*R(w_2)+2*R(w_3)+R(w_4)]

The time stepper and integrator RK4 is fourth order accurate. Not all schemes can be combined with all
options, infact the TVD scheme should use **time_integrator='onestep'**. Most methods like
TVDLF, TVD-MUSCL, HLL, HLLC can use
**time_stepper='twostep'** or **time_stepper='threestep'** or even 'fourstep' or 'fivestep'.
For consistency, the accuracy of the time integrator should have the same order as the slope
limiter's accuracy. 3rd order accurate time integrators are found within the 'threestep' options, with its default 'ssprk3' as **time_integrator**.
4th order accurate time steppers/integrators include 'fourstep' with 'rk4', or 'fivestep' with its default 'ssprk5'.
Various SSPRK schemes are found under the 'twostep', 'threestep', 'fourstep' and 'fivestep' time_steppers.

Since in the twostep method the **R1** spatial discretization in the first
_predictor_ step can be different from **R** of the second _corrector_ step,
the twostep time integration is a predictor-corrector type scheme. In the
parameter file, the arrays **typepred1** and **flux_scheme** (which need to be
specified for each AMR grid level up to **nlevelshi**), determine the method
applied in the predictor and full step.

For the multistep RK4 integration scheme, the same **flux_scheme** method is
used in each substep for **istep=1..nstep**. Full timesteps are counted by
**it_init &lt;= it &lt;= it_max**, while the physical time is **global_time &lt;= time_max**.

# Grid and Mesh {#disc-grid}

The grid is a 1, 2 or 3D grid, either Cartesian, cylindrical (to which polar
belongs), or spherical, which is selected by call `set_coordinate_system` 
in `usr_init` subroutine of `mod_usr.t`. See [coordinate system](axial.md) for 
more information.
 The coordinates of the grid points are represented by
the array **x**, usually interpreted for Cartesian coordinates as x, y, and z;
for cylindrical as r, z, and phi; for polar grid as r, phi, and z;
for spherical grid as r, theta, and phi. 
Other useful cell-related quantities calculated from **x**, like distance
between cell interfaces, volume, surface, and surface normal, are also
computed, and this happens in the _geometry.t_ module, specifically in the
_fillgeo_ subroutine. They are stored in the structures
_pw(igrid)%surfaceC(ixG^T,^D)_, _pw(igrid)%surface(ixG^T,^D)_,
 _pw(igrid)%dvolume_ and _pw(igrid)%dx_.

The spatial indices of these arrays, denoted by the dimension independent
[syntax](source.md), have varying ranges, depending on their use (as cell
centered, or cell surface quantities). Remember that 
number of physical cells (by _block_nx^D_) plus number of ghostcells 
(by _nghostcells_) gives the total grid extent:

    ixGlo^D:ixGhi^D = ixG^T

In many places, the same range, including ghost cells describing the
boundaries, is denoted as

    ixGmin^D:ixGmax^D = ixG^S

The mesh is defined as the grid without boundary layers:

    ixMlo^D,ixMhi^D = ixGlo^D+nghostcells,ixGhi^D-nghostcells

or equivalently as ixMmin^D,ixMmax^D, or in a shorter notation

    ixM^L = ixG^L^LSUBnghostcells^L

The ghost cells are updated by the subroutine `getbc` in `mod_ghostcells_update.t`. 
When a file is read or saved by MPI-AMRVAC, the ghost cells are usually not included.

You may run the [VACPP](vacpp.md) preprocessor interactively to see how the
above expressions are translated for a given number of dimensions.

# Boundary Regions {#disc-boundaries}

The boundary cells of the grid are grouped into boundary regions indexed by
**iB=1..nB**. In 1D, there are 2 regions, at left and right of the mesh; in 2D
there are 4, left, right, bottom, top (in that order); while for 3D there are
6. These boundary regions cover the **2*ndim** edges of the grid in the _left,
right, bottom, top, front, back_, in other words, _ixmin1, ixmax1, ixmin2,
ixmax2, ixmin3, ixmax3_, order. The regions overlap at the corner ghost cells.

The boundary methods are applied to each boundary region from **iB=1** to
**nB** and for each variable from **iw=1** to **nw** according to the
descriptor string **typeboundary(iw,iB)**. An exception is made for the **special**
type variables, where all variables should be provided somehow in conservative format
by users in pointed `usr_special_bc` subroutine in `mod_usr.t`.

Please note that the boundary itself is at the interface between the real and
ghost cells, therefore setting the velocity to 0 in the ghost cell will not
make the velocity 0 at the cell interface which is an interpolated value. To
have exactly zero flux through the boundary, the **'symm'** boundary type
should be used for the scalar quantities and the tangential velocity
components. For the normal vector components the anti-symmetric **'asymm'**
boundary type should be used.

## Internal boundaries {#disc-internal}

Internal boundaries can be used to overwrite the domain variables with
specified values. Internally, these are assigned before the ghost-cells and
external boundaries are applied (in subroutine `get_bc`). The user can provide
conditions on the conserved variables depending on location or time in the
pointed `usr_internal_bc` subroutine. To
activate internal boundaries, the switch

    internalboundary=.true.

has to be set in the _boundlist_ section of the par file.

