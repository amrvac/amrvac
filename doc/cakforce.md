# CAK radiation force module for HD/MHD

An overview of how to activate and work with the CAK radiation force module. This module only works for spherical problems in 1-D, 2-D, and 3-D.

# Physics

Starting from the general 7-D radiation transport equation it can be demonstrated that under the assumption of isotropic emissivities the net radiation acceleration resulting from the radiation force \f$ \mathbf{f}_\mathrm{r} \f$ is
\f[
\mathbf{g}_\mathrm{r}(\mathbf{r}) = \frac{\mathbf{f}_\mathrm{r}}{\rho(\mathbf{r})} = \frac{1}{\rho(\mathbf{r}) c} \oint_{4\pi} d{\Omega} \int_{0}^{+\infty} d{\nu}\,k_\nu(\mathbf{r},{\mathbf{n}})\,I_\nu(\mathbf{r},{\mathbf{n}})\,{\mathbf{n}}.
\f]
where \f$ I_\nu(\mathbf{r},\mathbf{{n}}) \f$ is the specific intensity in direction \f$ \mathbf{n} \f$ at a certain frequency at some spatial point \f$ \mathbf{r} \f$, \f$ k_\nu(\mathbf{r},\mathbf{n}) \f$ captures energy losses from the radiation beam, and \f$ \Omega \f$ is the solid angle where the radiation beam is directed into. 

This equation is further manipulated to arrive at an expression suitable to model stellar wind outflows from OB-stars following [Castor, Abbott, & Klein (1975)](https://ui.adsabs.harvard.edu/abs/1975ApJ...195..157C/abstract) (CAK). 
It should be noted that we do not follow the original CAK parametrisation (\f$ \alpha \f$, \f$ k \f$) for the ensemble of lines, but instead the conceptually advantageous Gayley (1995) parameterisation (\f$ \alpha, \bar{Q}, Q_0 \f$). For the conditions of OB-stars the net radiative acceleration may be written as
\f[
\mathbf{g}_\mathrm{r} = g_\mathrm{cont} \mathbf{e}_r + g_\mathrm{line}  \mathbf{e}_r = g_\mathrm{e}\,\mathbf{e}_r + f_\mathrm{d} g_\mathrm{line} \mathbf{e}_r
\f]
where \f$ g_\mathrm{e} = \kappa_\mathrm{e}L_\star/(4\pi r^2 c) \f$ is the (continuum) radiation acceleration due to Thomson scattering on free electrons with opacity \f$ \kappa_\mathrm{e}=0.34 \f$ cm\f$^2\f$/g. The radiation line force expression can take several forms depending on the physics included (see below). The classical radially streaming CAK expression for the force coming from a point star obeys
\f[
g_\mathrm{line} = \frac{\bar{Q}^{1-\alpha}}{(1-\alpha)} \frac{g_\mathrm{e}}{t_r^\alpha}
\f]
with \f$ t_r = \kappa_\mathrm{e} c \rho/(dv_r/dr) \f$ the Sobolev line-optical depth. To mimic a finite-sized star a radially-dependent pre-factor may be included
\f[
f_\mathrm{d} = \frac{(1+\sigma)^{1+\alpha} - (1+\sigma \mu_\star^2)^{1+\alpha}}{(1+\alpha)\sigma(1+\sigma)^\alpha (1-\mu_\star^2)}
\f]
where \f$ \sigma = d(\ln v_r)/d(\ln r) -1 \f$ and \f$ \mu_\star^2 = 1-(R_\star/r)^2 \f$. For numerical stability, we use in practise a rewritten form of the above expression.

Furthermore, an extension of the CAK line-force expression to a formulation valid for all optical depth ranges yields
\f[
g_\mathrm{line} = \bar{Q}g_\mathrm{e} \left[ \frac{(1+Q_0 t_r)^{1-\alpha} - 1}{(1-\alpha)Q_0 t_r} \right].
\f]
Note that the original CAK expression is recovered in the optically thick limit \f$ (t_r \gg 1) \f$ and when *assuming* \f$ \bar{Q}=Q_0 \f$. As explained below, both these line-force expressions can be employed when using the module.

If the original CAK point star assumption is relaxed, the full radiation force integral presented above has to be solved over all solid angles \f$ \Omega \f$. This *vector* force expression requires then numerical quadratures over the radiation polar \f$ \Theta \f$ and azimuthal angle \f$ \Phi \f$ (these are **not** the spherical polar and azimuthal coordinate angles!) associated with radiation rays emerging from the surface. Evidently, in this case no \f$ f_\mathrm{d} \f$ factor is required and we solve for
\f[
\mathbf{g}_\mathrm{line}(\mathbf{r}) = \frac{(\kappa_\mathrm{e} \bar{Q})^{1-\alpha}}{(1-\alpha)c^{\alpha+1}} \oint_\Omega \,d\Omega\, \mathbf{n}\,I_\star(\mathbf{n}) \left[ \mathbf{n}\cdot \nabla(\mathbf{n}\cdot \mathbf{v}) /\rho(\mathbf{r}) \right]^\alpha.
\f]
With this formula all components of the force can be computed: \f$ \mathbf{g}_\mathrm{line}(\mathbf{r})  = g_{\mathrm{line},r} \mathbf{e}_r + g_{\mathrm{line},\theta} \mathbf{e}_\theta + g_{\mathrm{line},\phi} \mathbf{e}_\phi \f$. The reader unfamiliar with the physics or who wishes to understand the details better can consult the thesis of F. Driessen (KU Leuven 2022). Particularly, section 2.1 and 2.2 present an extensive description and derivation of the radiation line force after CAK. The implementation in MPI-AMRVAC follows closely the above equations.

Finally, note that contrary to the RHD physics module of MPI-AMRVAC, the radiation force prescription presented here is a **reduced** dynamical picture. Indeed, we do not solve for the detailed radiation-energy exchange that is done via the flux-limited diffusion in the RHD module. This is justified since for the typical stellar conditions where the CAK force prescription applies (OB-star atmospheres), such detailed radiation-energy exchanges between the radiation field and gas are rather minor and instead the radiation dynamically couples to the gas only without considering the radiation energy budget.

# Practical use
In order to add the CAK force to a line-driven stellar wind simulation, the user has to adjust the .par file as indicated in this section.

## Activating the module
Depending whether the setup is in (M)HD, the the force needs to be activated in the physics-dependent namelist. This means that for a HD simulation one adds in the `&hd_list` the following:
```
&hd_list
  hd_cak_force = .true.
/
```
and a similar procedure applies for the `&mhd_list` with `mhd_cak_force`.

## Customizing the setup
Additionally, specific parameters can be set in the `&cak_list` for the source term physics treatment. In the following table, the available parameters are briefly described with their possible values. If some values are not set, they are default to a 1-D CAK force prescription for a proto-typical Galactic O-supergiant.

name | type | default | description
---|---|---|---
cak_alpha | double | 0.65 | Power-law index of the CAK line-ensemble distribution function. Allowed values are in the range \f$\alpha\in [0,1[\f$.
gayley_qbar | double | 2000 | Gayley (1995) ensemble-integrated line-strength of the line-ensemble distribution function. Number should be bigger than zero.
gayley_q0 | double | 2000 | Gayley (1995) line-strength cut-off for the line-ensemble distribution function. Extension to the original CAK formulation to avoid an inifinite line force in the optically thin limit. Number should be bigger than zero.
cak_1d_force | logical | F | If true, activate the original CAK 1-D (radial) line force prescription as source term to the (M)HD equations. Hence, computes \f$g_{\mathrm{line},r}\f$.
cak_1d_opt | integer | 1 | Switch that allows to select different CAK 1-D force prescriptions. When cak_1d_opt=0 we follow the original CAK radially streaming limit (point star). When cak_1d_opt=1 the 1-D radially streaming force gets adjusted with the finite disc factor to take into account the finite extent of the star. When cak_1d_opt=2 the 1-D radially streaming force is corrected with the finite disc factor and with the line-strength cut-off.
cak_vector_force | logical | F | If true, activate the CAK line force with components computed in all coordinate directions (radial, polar, azimuthal). Hence, computes \f$g_{\mathrm{line},r}\f$, \f$g_{\mathrm{line},\theta}\f$, and \f$g_{\mathrm{line},\phi}\f$.
nthetaray | integer | 6 | If `cak_vector_force=T`, amount of radiation ray points used to describe the radiation angle \f$\Theta\f$.
nphiray | integer | 6 | If `cak_vector_force=T`, amount of radiation ray points used to describe the radiation angle \f$\Phi\f$.
fix_vector_force_1d | logical | F | If true, the CAK vector line force will only have a non-zero radial component.
cak_split |  logical | F | To add the source term in a split or unsplit (default) fashion.

## Test problems
Two test problems are provided in the HD test folder that can be configured in several ways to illustrate the above physics, see `CAKwind_spherical_1D` and `CAKwind_spherical_2.5D`. Functionalities included are:
- 1D and 2.5D spherical isothermal or adiabatic wind models.
- Rotating frame of reference for 2.5D wind.
- Radial line force or vector line force with one or more of the previous listed combinations.

Additionally, one test problem is provided in the MHD test folder, see `CAKwind_Magnetosphere_spherical_2D`. This problem treats a non-rotating, trans-magnetosonic line-driven wind interacting with a stellar dipole field that creates a so-called "Dynamical Magnetosphere". This example problem can be run in isothermal conditions or in adiabatic conditions with radiative cooling included. For the latter see also the documentation on [radiative cooling](radiative_cooling.md).

For all test problems several .par files are provided to serve as a starting point.


# Numerics and implementation

The CAK radiation force module is at mod_cak_force.t. This module has some subroutines that are called in the mod_(m)hd_phys.t module to allow to add the CAK radiation force physics in the simulation.

An important routine that the user has to call in his/her own mod_usr.t file is the set_force_norm routine. Given the stellar radius and wind temperature it will compute the (dimensionless) line-force normalisation that will be used during the force term computation. Failing to call this routine will lead to no results and code crash.

## Primary subroutines

*cak_init:*
The subroutine cak_init is called only once, during initialization of the hd or mhd physics. It reads in parameters of cak_list and makes extra variables in the w-array for output; the radiation force component(s) and the finite disc factor. The latter is only present in 1-D force prescription models. Multiple force components are only present when the vector force is included. The direction of each coordinate is given by a Roman numeral (1,2,3).

*set_force_norm:* 
This subroutine computes the force normalisation in dimensionless units. **It has to be called once within the mod_usr.t** initialisation once the stellar radius and wind temperature are available. Essentially it computes the dimensionless stellar luminosity along with a dimensionless mass-absorption free-electron opacity and dimensionless speed of light variable.

*cak_add_source:*
The subroutine cak_add_source is called at each timestep to add the CAK line force together with the radial Thomson (continuum) force. Depending on the chosen option (`cak_1d_force` or `cak_vector_force`) either the original 1-D CAK line force or the full CAK vector force is computed. The forces are all computed in dimensionless units.

In case there is adiabatic (M)HD, the forces will also be added to the energy equation. To mimic stellar heating there is a floor temperature set (being the stellar surface temperature) to the wind. This also fixes the sometimes possible occurrence of negative thermal pressures in parts of the hypersonic wind outflow.

*cak_get_dt:*
Since an additional force is introduced to the (M)HD system, we also check for possible time step constraints in every coordinate direction. When the 1-D CAK force is used this only applies to the radial direction again, when the vector CAK force is used then it depends on how many dimensions the problem uses to determine possible additional time step limitations.

## Secondary subroutines

There are several additional subroutines that have been added to facilitate the CAK force source routine.
  
*get_cak_force_radial:* 
Computes the original 1-D CAK line force under various options set by `cak_1d_opt`. The default is the computation of the finite-disc corrected force. By default the array holding the line force will have three slots for each coordinate direction, but only the radial slot will be filled. 
This force computation can used in any dimension, but will always only modify the radial momentum equation.

*get_cak_force_vector:* 
Computes the full 3-D CAK line force. This will output all three line-force components, but with `fix_vector_force_1d` it is possible to only apply a pure radial force component to the dynamical equations of the system (this will be then similar to the original 1-D CAK force with finite disc correction, except with now a proper ray quadrature performed to get the finite extent of the star).
This force can also be used in any dimension since the non-radial force components are generally non-zero. E.g. in 2.5-D models the azimuthal line force component will actually lead to a radiative torque and spindown of the star.

*get_gelectron:*
Computes the Thomson (continuum) force from free-electron scattering. This force is purely radial and it is assumed that the mass-absorption coefficient for free-electron scattering is that of a fully ionised plasma at solar abundances. 
Note when the star becomes oblate or has a latitude-dependent surface temperature the Thomson force will also get a polar component, besides the radial component. In the present implementation such latitudinal dependency is not implemented.

*get_velocity_gradient:*
This routine computes a second-order spatially accurate gradient in a certain coordinate direction for a given scalar -- in our case the scalar will be one of the velocity field components. It uses stencil information in both forward, central, and backward direction. The formulations is done for an arbitrary non-uniform grid, which is typically used for line-driven stellar wind models in order to resolve the transonic region well. 
If a uniform grid is used, the formulae will collapse to a second-order spatially accurate gradient expression for uniform grids.

*rays_init:*
Initialisation routine for the radiation rays when the `cak_vector_force` is applied. It distributes a preselected amount of rays with radiation angles (\f$\Theta\f$, \f$\Phi\f$)$ on the (projected) stellar disc. Additionally, it assigns weights to the rays in the contribution to the radiation line force. For convenience the radiation polar angle points are recast to a formulation in terms of impact parameter points. Both the ray points and weights are set by a Gauss-Legendre quadrature.
This routine is called once in the *cak_init* routine when `cak_vector_force=T`.

*gauss_legendre_quadrature:*
Routine contains a fast and efficient algorithm to compute the Gauss-Legendre quadrature for an amount of points together with their corresponding weights. Prints the chosen ray points and weights to the screen for the user.
