# Physics modules and equations

[TOC]

# List of physics modules {#eq_list_physics}

This document describes the equations implemented.
Information about user defined source terms are in [user module](amrvacusr.md). In
principle, the code handles anything of generic form

![](figmovdir/eq.general.gif)

The code is configured to use the specified set of equations by activating it in 
the usr_init subroutine of user module "mod_usr.t"

    subroutine usr_init()
      ...
      call EQUATION_activate()
    end subroutine usr_init


where EQUATION is one of the implemented physics modules (rho,hd,mhd),
see below.

## Transport Equation: rho {#eq_rho}

    call rho_activate()

![](figmovdir/eq.rho.gif)

The transport equation describes the transport of a scalar field, here the
density **rho** by a prescribed velocity field. This equation is used for test
purposes.

The parameters rho_v in the rho_list of amrvac.par file
 define the components of the uniform velocity field.

For a linear scalar equation the Riemann solver is trivial, thus all TVD type
methods give identical results.

## Scalar Nonlinear Equation: nonlinear {#eq_nonlinear}

    call nonlinear_activate()

This module contains various instances of a scalar nonlinear equation, including the inviscid Burgers, inviscid nonconvex equation, as well as a possibility to handle the Korteweg-de Vries equation. It allows testing of truly nonlinear (shock steepening and formation) phenomena, in 1D to multi-D, or to test how source additions are best combined with flux prescriptions and discretizations. The equation implemented in \f$N_d\f$ dimensions is

\f[
\frac{\partial \rho}{\partial t} + \nabla \cdot \mathbf{F}(\rho,\mathbf{x},t) = -\delta^2 \sum_{i=1}^{N_d} \frac{\partial^3 \rho}{\partial x_i^3}
\f]

where the RHS is activated through the _mod_kdv.t_ module. The actual flux expression can be chosen (depending on the parameter `nonlinear_flux_type`) to be one of \f$ \mathbf{F}^{\mathrm{burgers}}  =  \frac{1}{2}{\rho^2}\mathbf{v}_0 \f$ or \f$ \mathbf{F}^{\mathrm{nonconvex}}  =  {\rho^3}\mathbf{v}_0 \f$, where we introduced \f$ \mathbf{v}_0=\sum_{i=1}^{N_d} \hat{\mathbf{e}}_i \f$

## Hydrodynamics: hd {#eq_hd}

    call hd_activate()

![](figmovdir/eq.hd.gif)

The Euler equations are solved for density **rho**, the momentum density
**m=rho*v** and the total energy density **e**. The pressure is a derived
quantity which is calculated from the conservative variables.

Parameters of hydrodynamics are read in the **hd_list** of parameter file.
There is a single equation parameter, the adiabatic index **hd_gamma**
(typical value is 5/3). 

This equation module can be combined with physical sources for
(local) optically thin radiative losses by set **hd_radiative_cooling=.true.**. 
see the [radiative cooling](radiative_cooling.md) page. Schematically, it
introduces terms as

![](figmovdir/eq.radloss.gif)

The HD module can also be combined with the external gravity module 
(_src/physics/mod_gravity.t_) for uniform gravity by set **hd_gravity=.true.**

![](figmovdir/eq.gravity.gif)

and for point gravity

![](figmovdir/eq.pointgrav.gif)

Note how the gravitational constant and the non-dimensionalization is taken
into the parameters _M_point_ and its location _x_point_.

To do adiabatic hydrodynamics (i.e. solve the hydro set without energy equation added) just add in hd_list of parameter file

    hd_energy=.false.

This special case includes the equations for pressureless dust and the Shallow Water equations, and writes generally as

![](figmovdir/eq.hdadiab.gif)

The system of adiabatic hydrodynamical equations are solved for the density
**rho** and the momentum density **m=rho*v**. The pressure is a function of
density only since an isentropic initial condition is assumed. There are two
equation parameters, the adiabatic index **hd_gamma** (the isothermal
case corresponds to **hd_gamma = 1**) and the adiabatic constant
**hd_adiab** (which should be positive or zero). It is possible to set
**hd_adiab=0** and handle the case of pressureless dust.

The system of **shallow water equations** is a special case with the following
identifications: **rho=h** represents the height of the water column,
**hd_gamma2** and the adiabatic coefficient is half of the gravitational
acceleration **hd_adiab=g/2**.

There is a Roe-type Riemann solver implemented, in _hd/mod_hd_roe.t_. Several
routines specific to HLLC are in _hd/mod_hd_hllc.t_.

## Magnetohydrodynamics: mhd {#eq_mhd}

    call mhd_activate()

![](figmovdir/eq.mhd.gif)

This is the full system of the MHD equations, with the following conservative
variables: density **rho**, momentum density **m=rho*v**, total energy density
**e** and the magnetic field **B**. The magnetic field is measured in units
for which the magnetic permeability is 1.

Parameters of magnetohydrodynamics are read in the **mhd_list** of parameter file.
The source terms on the right hand side with **eta** in them are the resistive
terms.

There are three equation parameters: the polytropic index **mhd_gamma**
(which must be larger or equal to 1), and the resistivity **mhd_eta**, and
the entropy **mhd_adiab**. Ideal MHD corresponds to **mhd_eta=0**,
positive values give a uniform resistivity, while a negative value calls the
**specialeta** procedure to determine the resistivity as a
function of the coordinates, of the conservative variables, and/or of the
current density. This subroutine is to be completed by the user.

There is a Roe-type Riemann solver implemented using arithmetic averaging, in
_mhd/mod_mhd_roe.t_, while several routines specific to HLLC are in _mhd/mod_mhd_hllc.t_.

This equation module can be combined with physical sources for
(local) optically thin [radiative losses](radiative_cooling.md) by set **mhd_radiative_cooling=.true.**. 
It can also be combined with the external gravity modules by set **mhd_gravity=.true.**.

### UAWSoM wave-energy extension

Setting **mhd_uawsom=.true.** adds four conserved wave-energy densities,
**wAplus**, **wAminus**, **wkplus**, and **wkminus**, to total-energy MHD.  The
first pair represents Alfvén waves and the second pair represents kink waves.
In MPI-AMRVAC's sign convention, the `plus` variables propagate against the
magnetic field and the `minus` variables propagate along it.  Their pressures
are

    p_A = (wAplus + wAminus)/2
    p_k = (zeta + 1) (wkplus + wkminus)/4,

and the total energy includes the sum of all four wave energies.  Their fluxes
use the bulk velocity plus or minus the local Alfvén or kink speed.  Expansion,
nonlinear dissipation, and the associated thermalization are applied as source
terms.  With **mhd_uawsom_reflection=.true.**,
**mhd_uawsom_reflection_mode='one_dimensional_gradient'** retains the original
one-dimensional gradient source.  **'cartesian_gradient_vorticity'** uses the total-field directional
gradient of the Alfvén speed and the field-aligned velocity vorticity, with a
rate limiter `R_lim,A=min(R_imb,A,max(Gamma_plus,Gamma_minus))`.  This
multidimensional closure does not use `mhd_uawsom_sigma`.  The exchange is
conservative: a signed transfer removes energy from one population and adds it
to the other.  **mhd_uawsom_kink_reflection=.true.** adds the corresponding
kink-speed directional-gradient exchange with
`R_lim,k=min(abs((V_k dot grad) ln V_k),max(Gamma_kplus,Gamma_kminus))`; it
also does not use `mhd_uawsom_sigma`, and the kink source does not include a
field-aligned vorticity term.  Both exchanges have a bounded 4:1 imbalance
factor, a donor positivity cap, and zero-energy protection.

The local density contrast `zeta`, unresolved thread radius, and Alfvén
correlation length can be supplied by **usr_uawsom_coefficients**.  Otherwise
the `mhd_uawsom_*` namelist scales are used.  B0-split runs evaluate all wave
speeds and closure lengths from the total magnetic field.  Optically thin
cooling is multiplied by the transverse-structure average

    1 + f (1-f) (zeta-1)^2 / (1 + f zeta - f)^2.

The `plus` populations propagate against **B** and the `minus` populations
along **B**.  In a B0-split run, **B** is the sum of the evolved perturbation
and the static background field.  The kink expansion-work source has the
positive paper-Eq. 4 sign, `+(zeta-1)/(zeta+1) p_k div(v)`.

See `tests/mhd/UAWSoM_1D`,
`tests/mhd/UAWSoM_reflection_2.5D`,
`tests/mhd/UAWSoM_reflection_3D`, and
`tests/mhd/UAWSoM_solar_atmosphere_2.5D`, together with the detailed
provenance map in `doc/uawsom_equation_map.md`.

We also have implemented the magnetic field splitting strategy, where a static, 
background magnetic field is assumed. This modifies the equations and brings in extra
sources and flux terms.

To run isothermal magnetohydrodynamics, add in mhd_list of parameter file

    mhd_energy=.false.

![](figmovdir/eq.mhdiso.gif)

This is the system of the MHD equations without the full energy equation, and
with the following conservative variables: density **rho**, momentum density
**m=rho*v**, and the magnetic field **B**. The magnetic field is measured in
units for which the magnetic permeability is 1. The density pressure relation
is polytropic.


# Divergence B source treatments {#eq_divB_fix}

Both the classical and the special relativistic MHD module can deal with
solenoidal magnetic field corrections through source term treatments.
Traditionally, these can be written as

![](figmovdir/eq.divb.gif)

Terms proportional to **div B** are [Powell`s fix](methods.md) for
the numerical problems related to the divergence of the magnetic field. They
are used only in more than 1D. We can also
just take the term along in the induction equation, known as Janhunen`s
approach. Another option is to use the diffusive (parabolic) approach, with
the parameter _C_d_ of order unity (up to 2). Alternatively, there is the 
[Dedner`s](methods.md) generalised Lagrange multiplier (GLM) method.

# Positivity fixes {#eq_positivity_fixes}

Another, similarly corrective, action is referred to as positivity fixing.
This is merely an additional means to handle the supposedly rare instances
where due to all nonlinearities of the scheme employed, the local conservative
to primitive transformation signals a non-physical state. Our positivity fix
approach can then be activated, and one such strategy operates as follows:
identify all cells (within the same grid block) that represent physical states
surrounding a faulty cell in a rectangular zone up to **small_values_daverage** cells
away; (2) convert those cells to primitive variables; and (3) for all but the
magnetic field components, replace the faulty cell values by the average of
surrounding physical state cells. Finally, revert to conservative variables
where needed. Obviously, in this form, strict conservation may be violated.
These fix strategies are seperated off in the _mod_small_values_ modules.
They are by default inactive, and can be controlled by the parameters
**small_values_method** and other related parameters described in
[par/PROBLEM](par.md).
