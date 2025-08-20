# Rotating frame for HD/MHD/RHD

An overview of how to activate the rotating frame module for Cartesian, cylindrical, and spherical problems in 1-D, 2-D, and 3-D.

# Physics

For some problems it may be desirable to let the coordinate system rotate with a certain rate and perform computations in a rotating frame of reference. For example, one can anchor then the coordinate system to the object under study to yield simpler boundary conditions. This method will, however, lead to the introduction of fictitious or "pseudo" forces in the governing equations of the system. Particularly, the Coriolis force and the centrifugal force need to be taken into account.

When active, the rotating_frame.t module adds these fictitious sources to the momentum and energy equations. At all times the rotation axis of the frame is assumed to point along the \f$ z \f$-axis of a Cartesian coordinate system. Similarly the rotating frame rotates at a **uniform** angular frequency \f$ \omega_0 \f$, meaning no fictitious Euler force can develop.

**NOTE**: The user should be aware that the fictitious forces are **not** added in an angular momentum conserving fashion. The latter means that the azimuthal component of the Coriolis force is incorporated directly into the azimuthal flux computation within the conservative equations. Such angular momentum conserving formulation may be important for some applications, see for example [Kley (1998)](https://ui.adsabs.harvard.edu/abs/1998A%26A...338L..37K/abstract). The present implementation in MPI-AMRVAC adds all fictitious force components displayed below as source terms to the conservative equations.

## Cartesian coordinates

The extra fictitious forces (per unit mass) added to the momentum equations:
\f{eqnarray*}{
 F_{x} &=& 2\omega_0 v_y + \omega_0^2 x, \\
 F_{y} &=& -2\omega_0 v_x + \omega_0^2 y,
\f}
and if doing adiabatic (R)(M)HD simulations the fictitious work due to the centrifugal force will be added
\f[
W = (x v_x + y v_y) \omega_0^2.
\f]

## Cylindrical coordinates

The extra fictitious forces (per unit mass) added to the momentum equations:
\f{eqnarray*}{
 F_{r}    &=& 2\omega_0 v_\phi + \omega_0^2 r, \\
 F_{\phi} &=& -2\omega_0 v_r,
\f}
and if doing adiabatic (R)(M)HD simulations the fictitious work due to the centrifugal force will be added
\f[
W = v_r r \omega_0^2.
\f]

## Spherical coordinates

The extra fictitious forces (per unit mass) added to the momentum equations:
\f{eqnarray*}{
 F_{r}  &=& 2\omega_0 v_\phi \sin \theta + r (\omega_0 \sin \theta)^2, \\
 F_{\theta} &=& F_{r} \cot \theta, \\
 F_{\phi} &=& -2\omega_0 \sin \theta (v_r + v_\theta \cot \theta),
\f}
and if doing adiabatic (R)(M)HD simulations the fictitious work due to the centrifugal force will be added
\f[
W = r (\omega_0 \sin \theta)^2 (v_r + v_\theta \cot \theta).
\f]

# Practical use
In order to activate the addition of fictitious forces, the user has to adjust the .par file as indicated in this section. For applications of the rotating frame physics, see the `CAKwind_spherical_2.5D` test problem in the HD tests folder or the `icarus` test problem in the MHD tests folder.

## Activating the module
Depending whether the setup is in HD, MHD, or RHD the rotating frame needs to be activated in the physics-dependent namelist. This means that for a HD simulation one adds in the `&hd_list` the following:
```
&hd_list
  hd_rotating_frame = .true.
/
```
and a similar procedure applies when employing the `&mhd_list` or `&rhd_list`. Additionally, a `&rotating_frame_list` list has to be supplied with the (dimensionless) angular frequency of the frame:
```
&rotating_frame_list
  omega_frame = DOUBLE
/
```
This value can be overwritten with a computed value in the mod_usr.t if the uniform angular frequency happens to depend on other user input parameters.

# Numerics and implementation

The rotating frame physics module is found in mod_rotating_frame.t. This module has some subroutines that are called in the mod_hd_phys.t, mod_mhd_phys.t, and mod_rhd_phys.t module to allow to initialise, compute, and add the fictitious forces to the physics governing equations.

## Primary subroutines

*rotating_frame_add_source:*
This subroutine computes the Coriolis and centrifugal contributions to the components of the momentum equation and the energy equation. Based on the geometry (cylindrical or spherical) and the problem dimension the source terms can assume different appearances.

