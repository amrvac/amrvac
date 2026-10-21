# Thermal conduction in MPI-AMRVAC

[TOC]

MPI-AMRVAC offers two approaches: parabolic conduction with super time stepping
(STS), and hyperbolic conduction with a relaxing heat flux. This page describes
both, using MHD as the main example. Set the [normalization units](amrvac_units.md)
in `usr_init()` before activating the physics module.

## Parabolic conduction {#tc-parabolic}

Enable `mhd_thermal_conduction=.true.` in `&mhd_list`, or
`hd_thermal_conduction=.true.` in `&hd_list`. The conduction equation is advanced
separately from the fluid equations using STS, including the RKL2 method.

The `&tc_list` provides these options:

- `tc_perpendicular=.true.`: add perpendicular conduction in MHD.
- `tc_saturate=.true.`: include heat-flux saturation.
- `tc_slope_limiter`: choose the transverse-gradient limiter.

The first two options default to false. If the STS update produces oscillations,
reduce `sts_dtpar` in `&sts_list` from its default 0.5. These settings apply to
the parabolic method; hyperbolic conduction uses the parameters below.

## Hyperbolic conduction: getting started {#tc-configuration}

For field-aligned conduction, use:

```fortran
&mhd_list
  mhd_energy = .true.
  mhd_hyperbolic_tc = .true.
/
```

This selects Spitzer conduction and disables `mhd_thermal_conduction` if it was
also enabled. To add perpendicular transport in 2D/3D, include these settings
in the same namelist:

```fortran
  mhd_hyperbolic_tc_use_perp = .true.
  mhd_hyperbolic_tc_perp_mode = 'fixed_reference'
  mhd_hyperbolic_tc_kappa_perp_factor = 1.d-3
```

Here \f$10^{-3}\f$ is an illustrative prescribed conductivity ratio.
Other choices are described in [Perpendicular closures](@ref tc-hyperbolic-closures).
Initialize the heat fluxes in `mod_usr.t`, for example:

```fortran
if (mhd_hyperbolic_tc) then
  w(ixO^S,qpar_) = 0.d0
  if (mhd_hyperbolic_tc_use_perp) w(ixO^S,qperp_) = 0.d0
end if
```

Zero initial flux gives a relaxation transient. The output variables are `q`
and, when enabled, `qperp`. Their boundary types default to the density boundary
type; a `special` boundary should supply heat-flux values consistent with its
thermal condition.

## Field-aligned hyperbolic conduction {#tc-hyperbolic-parallel}

The field-aligned method is described in
[Zhou et al. (2025), Section 2.2](https://doi.org/10.3847/1538-4357/ad96af).
It evolves a scalar heat flux toward the classical conductive flux:
\f[
 \boldsymbol q=q_\parallel\boldsymbol b,\qquad
 \tau\frac{\partial q_\parallel}{\partial t}+q_\parallel
 =-f_{\rm sat}\kappa_\parallel\boldsymbol b\cdot\nabla T,
\f]
where \f$\boldsymbol b\f$ is the unit vector along the total magnetic field,
including B0. The energy equation receives \f$-\nabla\cdot\boldsymbol q\f$.

By default, \f$\kappa_\parallel=\kappa_0T^{5/2}\f$, with the Spitzer prefactor
\f$\kappa_0=8\times10^{-7}\f$ in cgs or \f$8\times10^{-12}\f$ in SI,
converted internally to code units. A positive `mhd_hyperbolic_tc_kappa`
instead selects a constant conductivity in code units.

### Relaxation and saturation {#tc-hyperbolic-relaxation}

The MHD relaxation time is
\f[
 \tau=\max\left(4\Delta t,
 \frac{f_{\rm sat}\kappa_\parallel T_{\rm cond} C_{\rm CFL}^2}
 {[p/(\gamma-1)]c_{\rm f,max}^2}\right).
\f]
Here \f$C_{\rm CFL}\f$ is `courantpar`, and \f$c_{\rm f,max}\f$ is the largest
local classical fast magnetosonic speed over coordinate directions.
Normally \f$T_{\rm cond}=T\f$; with TRAC and Spitzer conduction, the conductivity
and relaxation estimate use \f$\max(T,T_{\rm cutoff})\f$, while the gradient
remains that of \f$T\f$.

Without saturation, \f$f_{\rm sat}=1\f$. Setting
`mhd_hyperbolic_tc_sat=.true.` uses
\f[
 f_{\rm sat}=\frac{1}{1+|\boldsymbol q_{\rm cl}|/q_{\rm sat}},
 \qquad q_{\rm sat}=1.5\rho(p/\rho)^{3/2},
\f]
where \f$\boldsymbol q_{\rm cl}\f$ is the classical target flux.

The finite relaxation time avoids the explicit parabolic timestep restriction,
but affects rapid thermal transients. Check accuracy against a resolved
parabolic reference; the flux solver and reconstruction also contribute
numerical diffusion.

## Perpendicular hyperbolic conduction {#tc-hyperbolic-perpendicular}

The current MHD implementation adds a new perpendicular channel to the
field-aligned method. Define
\f[
 \nabla_\perp T=\nabla T-\boldsymbol b(\boldsymbol b\cdot\nabla T),
 \qquad \boldsymbol n_\perp=\frac{\nabla_\perp T}{|\nabla_\perp T|}.
\f]
The second scalar flux evolves as
\f[
 \tau\frac{\partial q_\perp}{\partial t}+q_\perp
 =-f_{\rm sat}\kappa_\perp|\nabla_\perp T|,\qquad
 \boldsymbol q=q_\parallel\boldsymbol b+q_\perp\boldsymbol n_\perp.
\f]
Thus `qperp_` is a signed amplitude along the perpendicular temperature
gradient. The direction is set to zero where that gradient vanishes.
Both channels share \f$\tau\f$ and, when enabled, a saturation factor based
on the magnitude of their combined classical flux. This two-channel extension
is new work in MPI-AMRVAC, separate from the method in Zhou et al. (2025).

### Perpendicular closures {#tc-hyperbolic-closures}

Choose `mhd_hyperbolic_tc_perp_mode` according to the intended model:

| Mode | \f$\kappa_\perp/\kappa_\parallel\f$ | Use |
| --- | --- | --- |
| `fixed_reference` | `mhd_hyperbolic_tc_kappa_perp_factor` | Fixed reference or prescribed effective ratio. |
| `weak_field_isotropization` | \f$B_{\min}^2/(B^2+B_{\min}^2)\f$ | Empirical transition to isotropic conduction in weak fields. |
| `electron_magnetization` | \f$1/(1+\chi_e^2)\f$ | Simplified electron-magnetization model using local B, T, and electron density. |

The electron-magnetization estimate scales as
\f$\chi_e\propto BT^{3/2}/(n_e\ln\Lambda)\f$, with electron density supplied
by the EOS and a default Coulomb logarithm of 20. It describes electron
transport, without additional neutral conduction or electron-neutral collisions.

## MHD parameter summary {#tc-mhd-parameters}

All parameters belong to `&mhd_list`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `mhd_hyperbolic_tc` | `.false.` | Enable hyperbolic conduction; requires an energy equation. |
| `mhd_hyperbolic_tc_sat` | `.false.` | Enable heat-flux saturation. |
| `mhd_hyperbolic_tc_kappa` | `0` | Zero: Spitzer conductivity. Positive: constant conductivity in code units. |
| `mhd_hyperbolic_tc_use_perp` | `.false.` | Enable the perpendicular channel in 2D/3D. |
| `mhd_hyperbolic_tc_perp_mode` | `'off'` | Select one of the three closures above when enabling the perpendicular channel. |
| `mhd_hyperbolic_tc_kappa_perp_factor` | `0` | Fixed-reference ratio; zero selects an automatic normalization-based value. |
| `mhd_hyperbolic_tc_Bmin` | `0` | Weak-field transition scale; positive values are in code magnetic-field units. Zero selects `0.1/unit_magneticfield`. |
| `mhd_hyperbolic_tc_coulomb_log` | `20` | Positive Coulomb logarithm for `electron_magnetization`. |

For `fixed_reference`, the automatic ratio is
`4.d-10*unit_numberdensity**2/unit_magneticfield**2/unit_temperature**3`
in cgs, with prefactor `4.d-30` in SI. It is fixed at initialization.
For `weak_field_isotropization`, specify a positive code-unit Bmin when a
particular physical scale is intended; the default expression above applies
literally in either unit system.

To turn off perpendicular conduction, set `mhd_hyperbolic_tc_use_perp=.false.`.
A zero reference factor selects the automatic value, while enabling the channel
with mode `'off'` is rejected.

## Other physics modules {#tc-module-support}

| Module | Enable in its physics namelist | Support |
| --- | --- | --- |
| MHD | `mhd_hyperbolic_tc` | Parallel in 1D/2D/3D; perpendicular in 2D/3D. |
| FFHD | `ffhd_hyperbolic_tc` | Parallel conduction in 2D/3D along the frozen field. |
| HD | `hd_hyperbolic_thermal_conduction` | 1D conduction with a face-based heat-flux update. |

FFHD and HD use the scalar index `q_`. Their saturation switches are
`ffhd_hyperbolic_tc_sat` and `hd_htc_sat`, respectively. FFHD uses the global
hydrodynamic signal speed for relaxation and disables its parabolic conduction
when hyperbolic conduction is enabled. HD requires the user to disable its
parabolic option. Its `hd_htc_kappa_override` changes the Spitzer prefactor,
retaining the \f$T^{5/2}\f$ dependence.

## Examples {#tc-examples}

- [Oblique Fourier conduction](https://github.com/amrvac/amrvac/tree/amrvac3.3/tests/demo4/HyperbolicConduction_Fourier):
  a uniform-grid comparison of parallel, oblique, and perpendicular modes
  against an RKL2 parabolic reference, varying resolution and relaxation time.
- [Flare current-sheet conduction](https://github.com/amrvac/amrvac/tree/amrvac3.3/tests/demo4/HyperbolicConduction_CSHKP):
  a 2.5D application comparing perpendicular closures and prescribed ratios.

Each example's README provides build, run, and analysis instructions.
