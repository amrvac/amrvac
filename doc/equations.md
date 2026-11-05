title: Physics modules and equations

# List of physics modules {: #eq_list_physics }
This document describes the equations implemented. User-defined source
terms are described in [User source terms](#eq_usr_source) below. In
principle, the code handles anything of generic form

\begin{equation}
    \partial_t \mathbf{U} + \nabla \cdot \mathbf{F}(\mathbf{U}) = \mathbf{S}(\mathbf{U},\nabla \mathbf{U},\mathbf{x},t)
\end{equation}

The code is configured to use the specified set of equations by activating it via 
the code's parameter file in the methodlist

```fortran
&methodlist
   phys            = 'mhd'
/
```

where `phys` is one of the implemented physics modules (`hd`, `mhd`, `ffhd`, `srhd`, `srmhd`),
see below. The equations below are written in flat (Cartesian) space; running
on a spherical or cylindrical mesh instead adds extra geometric source terms
to the momentum (and, for `mhd` and `srmhd`, magnetic field) equations, described
separately in [Coordinate systems](coordinate_systems.html#geom_sources).
The equations and notation below follow
[Porth et al., "Astrophysics on GPUs: introducing AGILE 1.0"
(submitted to RASTI, in review; arXiv:2607.19277)](https://arxiv.org/abs/2607.19277),
Section 2.3, which is also the citable reference for AGILE 1.0.

| Module | Available solvers |
| --- | --- |
| `hd`   | `LF`, `HLL`, `HLLC` |
| `ffhd` | `LF` |
| `mhd`  | `LF`, `HLL`, `HLLC` |
| `srhd` | `LF`, `HLL` |
| `srmhd` | `LF`, `HLL` |

The scheme is chosen with `flux_scheme` in `&methodlist` (`'tvdlf'` for
`LF`, `'hll'`, `'hllc'`; `'tvdlf'` is the default). `flux_scheme` is an
array with one entry per refinement level, as in MPI-AMRVAC, but **AGILE
uses only its first entry, for every level**: different schemes on
different levels are not implemented yet. The slope limiter, `limiter`, is
applied per level. `srhd`'s `HLL` uses the relativistic Davis (1988) min/max
signal-speed estimate in `estimate_speeds_minmax`. `HLLC` is a shared frame,
`riemann_hllc_prim`, around a per-module star-region flux, `hllc_flux`: `hd`
uses Toro's (2010) variant 2 with PVRS signal speeds, and `mhd` the solver of
Li (2005) with the Davis fast-magnetosonic speeds of `HLL`, as upstream
MPI-AMRVAC does. In `mhd` the normal field and GLM's `psi` still take the
`LF` flux, exactly as under `HLL`. `srhd` has no `HLLC`. A solver not listed for a module in the table is
not supported by it: the module has no signal-speed estimate for it, only a
placeholder that calls `mpistop`.

## Compile-time switches {: #eq_compile_time }

Most of the physics options below are **compile-time** switches: they become
fypp defines through `make/config_schema.toml`, so they select which source
terms and variables are compiled in. Set them in the case's `agile.par`,
which is the file the build reads, and rebuild. Setting one only in the
parameter file of a run, without rebuilding, has no effect. These are:

| Module | Compile-time switches in its `_list` |
| --- | --- |
| `hd`   | `hd_gravity`, `hd_radiative_cooling`, `hd_source_usr`, `hd_n_tracer` |
| `mhd`  | `mhd_gravity`, `mhd_radiative_cooling`, `mhd_resistivity`, `mhd_hyperbolic_thermal_conduction`, `mhd_hyperbolic_thermal_conduction_anisotropic`, `mhd_source_usr`, `mhd_energy_only`, `mhd_n_tracer` |
| `ffhd` | `ffhd_gravity`, `ffhd_radiative_cooling`, `ffhd_hyperbolic_thermal_conduction`, `ffhd_pdivb`, `ffhd_source_usr` |
| `srhd` | `srhd_source_usr`, `srhd_n_tracer` |

Everything else in the `_list`s, for example `hd_gamma`, `mhd_eta`,
`srhd_eos` or the [density and pressure floors](#eq_floors)
`*_small_density`/`*_small_pressure`, is an ordinary run-time parameter.

## Hydrodynamics: hd {: #eq_hd }

```fortran
&methodlist
   phys            = 'hd'
/
```

The hydrodynamic module solves the partial differential equations that
govern compressible gas dynamics. The conservative variables are density,
vectorial momentum density, hydrodynamic energy density and a
user-selected set of `n_tr` tracers, collected in
**U** = [ρ, **m**≡ρ**v**, e_hd, D_tr^i], obeying

\begin{eqnarray}
 \partial_t \rho + \nabla \cdot (\rho\mathbf{v}) & =  & 0  \,,\\
  \partial_t \mathbf{m} + \nabla \cdot \left(\mathbf{v}\mathbf{m}+p\mathbf{I}\right) & = & \rho\mathbf{g} \,,\\
  \partial_t e_{\mathrm{hd}}+\nabla\cdot\left[\mathbf{v}(e_{\mathrm{hd}}+p)\right]
     & = &
     \rho\mathbf{v}\cdot\mathbf{g} -n_e  n_H\Lambda(T) \,,\\
  \partial_t D_{tr}^i+\nabla \cdot(\mathbf{v} D_{tr}^i) & = & 0 \,.
 \end{eqnarray}

Here *p* is the gas pressure, **v** the velocity vector and **I** the
identity tensor. Tracers are handled using their conservative variable
D_tr=ρϑ_tr, such that the actual tracer value ϑ_tr obeys a pure advection
equation. The pre-implemented sources on the right-hand-side allow the
user to set the external gravitational acceleration **g** (see
[Gravity](#eq_gravity) below) and select the treatment of local and
instantaneous optically thin radiative losses (see
[Radiative losses](#eq_radloss) below).

The pressure is a derived quantity, obtained from the energy closure

\begin{equation}
    e_{\mathrm{hd}}=\frac{p}{\hat{\gamma}-1}+\frac{1}{2}\rho v^2 \,,
\end{equation}

which introduces the equation parameter **hd_gamma** representing the
constant ratio of specific heats γ̂ (typical value 5/3 for a mono-atomic
gas). An ideal gas law relates *p* = 𝓡ρT, with a fixed gas constant
𝓡 = k_B/(μ̃m_p), where k_B is the Boltzmann constant, m_p the proton mass,
and the constant mean molecular weight μ̃ = (1+4He_a)/(2+3He_a) is set by
the fixed Helium abundance **He_abundance** (He_a ≡ n_He/n_H).

Parameters of hydrodynamics are read in the **hd_list** of the parameter
file, including the equation parameter **hd_gamma** and the number of
tracers **hd_n_tracer**.

This equation module can be combined with physical sources for
(local) optically thin radiative losses by setting **hd_radiative_cooling=.true.**,
see the [radiative cooling](radiative_cooling.html) page. It can also be
combined with the external gravity module by setting **hd_gravity=.true.**,
see [Gravity](#eq_gravity) below, and with [user source
terms](#eq_usr_source) by setting **hd_source_usr=.true.**. All three are
compile-time switches, set in `agile.par` (see [Compile-time
switches](#eq_compile_time)).

`hd_list` also reads **hd_energy** and **hd_adiab**, inherited from
MPI-AMRVAC, where `hd_energy = .false.` selects an isothermal gas with
p = `hd_adiab` ρ^γ̂. AGILE does not support that: the flux and conversion
kernels always evolve the energy equation, and `hd_adiab` is not used.
Leave `hd_energy` at its default, `.true.`.

### Radiative losses {: #eq_radloss }
When **hd_radiative_cooling**/**mhd_radiative_cooling**/**ffhd_radiative_cooling**
is enabled, the energy equation of the corresponding module is augmented
with the optically thin radiative loss term

\begin{equation}
    -n_e n_H \Lambda(T) \,,
\end{equation}

where Λ(T) is a tabulated radiative loss function depending only on the
local temperature T (see [Radiative cooling](radiative_cooling.html) for
the available curves and the exact integration scheme of Townsend 2009
used to advance this term). These losses scale with the squared density:
ρ is converted to electron and Hydrogen number densities under the
assumption of constant, fully ionized conditions with a fixed Helium
abundance He_a ≡ n_He/n_H, so that charge neutrality sets

\begin{equation}
    n_e=n_H(1+2\,\mathrm{He}_a) \,,\qquad
    \rho=n_H m_p(1+4\,\mathrm{He}_a) \,.
\end{equation}

### Gravity {: #eq_gravity }
The external gravitational acceleration **g** enters the momentum and
energy sources (ρ**g** and ρ**v**·**g** respectively) of the `hd`, `mhd`
and `ffhd` modules alike, activated by setting **hd_gravity**/**mhd_gravity**/**ffhd_gravity**
`= .true.` in the corresponding `_list` of `agile.par` (a compile-time
switch, see [Compile-time switches](#eq_compile_time)). AGILE does
not ship a built-in gravity profile (e.g. no separate "uniform" vs "point
mass" switches); instead, the user supplies **g** directly for the
current cell state by implementing a `gravity_field` function in
`mod_usr.fpp`, for example a constant vertical gravity as used in the
`hd/RT3D` test case:

```fortran
pure real(dp) function gravity_field(wCT, x, idim) result(field)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)       :: wCT(nw_phys)
  real(dp), intent(in)       :: x(1:ndim)
  integer, value, intent(in) :: idim

  if (idim == 1) field =  0.0_dp
  if (idim == 2) field =  0.0_dp
  if (idim == 3) field = -1.0_dp
end function gravity_field
```

`idim` selects the spatial direction (1, 2 or 3) for which the
acceleration component is requested; a spatially- or time-varying field
(e.g. a point-mass potential) can be implemented by using `x` (and any
module-level state the user maintains) inside this function.

`gravity_field` is called on the GPU, from inside the source-term kernels,
so it must be a `pure` function marked with the generic device-routine macro
`${GPU_ROUTINE_SEQ()}$`. Do not write `!$acc routine seq` directly: the macro
expands to the matching directive of whichever offload backend the build
uses (OpenACC or OpenMP) and to nothing in a CPU build. Using it needs the
macro definitions included at the top of `mod_usr.fpp`, as every test case
does:

```fortran
#:include "mod_gpu_directives.fpp"
```

Give the file name only, without a path. The build passes `$AGILE_DIR/src`
to fypp as an include directory, so this works wherever the case directory
is. Older cases wrote a path relative to the repository, such as
`"../../../src/mod_gpu_directives.fpp"`; that still works, but only while
the case stays at that depth inside the repository, so do not copy it.

### User source terms {: #eq_usr_source }

Every physics module can add a user-defined, local source term to every
cell. It is enabled with **hd_source_usr**/**mhd_source_usr**/**ffhd_source_usr**/**srhd_source_usr**
`= .true.` in `agile.par` (a compile-time switch), and the physics module
then calls a subroutine named `addsource_usr` that the case's `mod_usr.fpp`
must provide:

```fortran
subroutine addsource_usr(qdt, qt, wCT, wCTprim, wnew, x, qsourcesplit)
  ${GPU_ROUTINE_SEQ()}$
  double precision, intent(in)    :: qdt, qt                    ! step and time
  double precision, intent(in)    :: wCT(nw_phys), wCTprim(nw_phys)
  double precision, intent(in)    :: x(ndim)                    ! cell position
  double precision, intent(inout) :: wnew(nw_phys)              ! add qdt*S here
  logical, intent(in)             :: qsourcesplit

  ! e.g. a heating rate H_user(x) on the energy equation:
  wnew(iw_e) = wnew(iw_e) + qdt * ...
end subroutine addsource_usr
```

Like `gravity_field`, it is called by name from inside the source-term
kernel on the GPU, so it is a device routine and must not assume anything
outside the one cell it is given. `wCT` and `wCTprim` are the state the
stage evaluates its update on, in conserved and primitive variables; `wnew`
is the state being updated. It is always added unsplit (`qsourcesplit` is
`.false.`), and it does not limit the time step: a stiff user source needs a
small enough `courantpar` or `dtpar`. `tests/ffhd/TI3D` is a working example.
The `usr_source` procedure pointer in
`mod_usr_methods` is a leftover of MPI-AMRVAC and is never called.

## Magnetohydrodynamics: mhd {: #eq_mhd }

```fortran
&methodlist
   phys            = 'mhd'
/
```

The MHD module combines mass conservation and the tracer equation above
with the following generalizations of momentum and energy, and adds the
induction equation for the magnetic field **B**:

\begin{equation}
\partial_t \mathbf{m} + \nabla \cdot \left(\mathbf{v}\mathbf{m}+(p+{\textstyle{\frac{1}{2}}}B^2)\mathbf{I}-\mathbf{B}\mathbf{B}\right)  = \rho\mathbf{g}
\end{equation}

\begin{equation}
  \partial_t e_{\mathrm{mhd}}+\nabla\cdot\left[\mathbf{v}(e_{\mathrm{mhd}}+p+{\textstyle{\frac{1}{2}}}B^2)-\mathbf{B}(\mathbf{v}\cdot\mathbf{B})+\mathbf{q}\right]
  = \rho\mathbf{v}\cdot\mathbf{g} -n_e  n_H\Lambda(T) + \nabla \cdot \left[ \mathbf{B} \times \eta \mathbf{J}\right]\,,
\end{equation}

\begin{equation}
\partial_t\mathbf{B}+\nabla\cdot\left[\mathbf{v}\mathbf{B}-\mathbf{B}\mathbf{v} +\psi\mathbf{I}\right] = -\nabla\times\eta \mathbf{J} \,,
\end{equation}

where the total energy density closure is

\begin{equation}
    e_{\mathrm{mhd}}=\frac{p}{\hat{\gamma}-1}+\frac{1}{2}\rho v^2+\frac{1}{2}B^2 \,,
\end{equation}

**q** is the conductive heat flux, which is zero unless [hyperbolic thermal
conduction](#eq_mhd_htc) is switched on. **J** = ∇×**B** is the current density, and the magnetic field is
measured in units for which the magnetic permeability is 1. The resistive
terms, those with η, are compiled in only with **mhd_resistivity=.true.** in
`agile.par`; they are then added with the compact stencil detailed in Porth
& Xia et al. (2014), and also limit the time step. The GLM field ψ used to
control the divergence of **B** numerically is described in [Divergence B
source treatments](#eq_divB_fix) below.

Parameters of magnetohydrodynamics are read in the **mhd_list** of the
parameter file. The equation parameters are the ratio of specific heats
**mhd_gamma** and the resistivity **mhd_eta**, which only has an effect when
`mhd_resistivity` is on.

The gas law p = n k_B T counts the particles per hydrogen nucleus as
b = 1 + `H_ion_fr` + `He_abundance` (`He_ion_fr` (`He_ion_fr2` + 1) + 1),
with the ionisation fractions `H_ion_fr` = H⁺/(H⁺+H),
`He_ion_fr` = (He²⁺+He⁺)/(He²⁺+He⁺+He) and `He_ion_fr2` = He²⁺/(He²⁺+He⁺),
all 1 (fully ionised) by default; this gives b = 2 + 3 He_a as for `hd`.
`He_abundance = 0` with `H_ion_fr = 0` gives a neutral hydrogen gas,
p = n k_B T. The magnetic field unit is set to √(4π p_unit).

This equation module can be combined with physical sources for (local)
optically thin [radiative losses](radiative_cooling.html) by setting
**mhd_radiative_cooling=.true.**, with the external gravity module by
setting **mhd_gravity=.true.**, and with [user source
terms](#eq_usr_source) by setting **mhd_source_usr=.true.**; all three are
compile-time switches. **mhd_energy_only=.true.** freezes ρ, **v** and
**B** and evolves only the energy equation, for example to follow
conduction or cooling on a fixed background.

### Hyperbolic thermal conduction {: #eq_mhd_htc }

With **mhd_hyperbolic_thermal_conduction=.true.** (compile time), `mhd`
evolves the field-aligned heat flux q as an extra variable, relaxing it
towards the Spitzer value instead of solving a parabolic equation. It enters
the energy equation above as the heat flux **q** = q **b̂**, and obeys

\begin{equation}
    \partial_t q = -\frac{q + f_{\mathrm{sat}}\,\kappa_\parallel\,\mathbf{\hat{b}}\cdot\nabla T}{\tau} \,,
\end{equation}

so that it tends to −f_sat κ_∥ **b̂**·∇T, the conductive flux. The
conductivity is κ_∥ = κ_∥,0 T^(5/2) with the Spitzer coefficient
κ_∥,0 = 8×10⁻⁷ (cgs, converted to code units in `phys_units`), or the
constant **tc_kappa_par** when that is set positive. The flux saturates as
f_sat = 1/(1 + |κ_∥ **b̂**·∇T|/q_sat), with q_sat = 1.5 ρ (p/ρ)^(3/2). The
relaxation time is τ = max(4Δt, f_sat κ_∥ T C² (γ̂−1)/(p c²)), where C is
`courantpar` and c the largest fast speed of the cell.

**mhd_hyperbolic_thermal_conduction_anisotropic=.true.** (which implies the
switch above) evolves the full heat-flux vector **q** instead, which then
enters the energy flux directly, and relaxes
its parallel and perpendicular parts separately, each with its own
saturation and τ. The perpendicular conductivity is
κ_⊥ = κ_∥/(1 + χ²), with the magnetisation χ = C_χ B T^(3/2)/ρ (code
units; C_χ is set in `phys_units`), or the constant **tc_kappa_perp** when
that is set positive. Both conductivity coefficients κ_∥,0 and C_χ are
derived from the units and are not parameters.

### Divergence B source treatments {: #eq_divB_fix }

To control magnetic monopole buildup numerically, AGILE uses the
generalized Lagrangian multiplier (GLM) approach of Dedner et al. (2002),
where ψ is an added evolved variable obeying

\begin{equation}
    \partial_t \psi+c_{h}^{2} \nabla \cdot \mathbf{B} = -\frac{c_h^2}{c_p^2}\psi \,.
\end{equation}

This equation has a hyperbolic advection speed c_h and a parabolic
factor c_p. In practice, the left-hand side is dealt with using the flux
prescription c²_max**B**, while the right-hand side (the damping term) is
solved exactly over the time step Δt_n through

\begin{equation}
    \psi^{n+1}=\psi^n\exp\left(-\Delta t_n c_{\mathrm{max}}\alpha_{\mathrm{GLM}}/\min_{i=1,3}(\Delta x_i)\right)
\end{equation}

using second-order operator splitting. This implies availability of the
maximal, global physical propagation speed c_max that also sets the time
step, while

\begin{equation}
    \alpha_{\mathrm{GLM}}=\min_{i=1,3}(\Delta x_i)\,c_h/c_p^2
\end{equation}

becomes a fixed parameter, **mhd_glm_alpha** in `&mhd_list`, set to 0.5 by
default but generally taken within [0,1] (Mignone, Tzeferacos & Bodo 2010).
Δx_i is the physical extent of the cell, which on a curvilinear grid varies
from cell to cell.

## Frozen-field hydrodynamics: ffhd {: #eq_ffhd }

```fortran
&methodlist
   phys            = 'ffhd'
/
```

As an intermediate towards full 3D MHD applications, the frozen-field
hydrodynamic (FFHD) module introduced by Zhou, Li & Keppens (2024, 2025)
solves hydrodynamics along a given *frozen* 3D magnetic topology
**B**(**x**) that does not evolve in time. Introducing the local
magnetic unit vector b̂ = **B**/B, this module evolves the set
**U** = [ρ, m_∥=ρv_∥, e_hd∥, q_∥], the magnetic-field-line-projected
variant of the hydrodynamic mass, momentum and energy equations above,
optionally augmented with a field-aligned heat flux q_∥:

\begin{eqnarray}
    \partial_t \rho+\nabla \cdot  \left (  \rho v_{\parallel} \mathbf{\hat{b}}\right ) &=& 0\,,\\
    \partial_t m_{\parallel}+\nabla \cdot \left [  \left (  \rho v^2_{\parallel} +p\right ){\mathbf{\hat{b}}}\right ] &=& \rho \mathbf{g}\cdot\mathbf{\hat{b}} + p (\nabla\cdot \mathbf{\hat{b}})\,,\\
    \partial_t e_{hd\parallel} + \nabla \cdot \left[\left((e_{hd\parallel} + p) v_{\parallel}+q_\parallel\right) \mathbf{\hat{b}}\right] &=& \rho v_{\parallel} \mathbf{g}\cdot\mathbf{\hat{b}}  -n_e  n_H\Lambda(T) + H_{\mathrm{user}}(\mathbf{x})\,,\\
    \partial_t q_\parallel & = & -\frac{q_\parallel+\kappa_\parallel\mathbf{\hat{b}}\cdot \nabla T}{\tau} \,,
\end{eqnarray}

where now e_hd∥ = p/(γ̂-1) + ρv_∥²/2. Each source term is optional:

- **The heat flux q_∥** exists only with
  **ffhd_hyperbolic_thermal_conduction=.true.** (compile time). It uses a
  hyperbolic approach to purely field-aligned thermal conduction, relaxing
  q_∥ towards the conductive flux −κ_∥ **b̂**·∇T with
  κ_∥ = κ_∥,0 T^(5/2). Note the sign: q_∥ is the heat flux itself and
  enters the energy flux with a plus sign. κ_∥,0 is not a parameter: it is
  the Spitzer coefficient 8×10⁻⁷ (cgs), converted to code units in
  `phys_units` (the variable `hypertc_kappa`). The relaxation time is
  τ = max(4Δt, κ_∥ T C² (γ̂−1)/(p c²_max)), with C = `courantpar` and c_max
  the global maximum signal speed, so q_∥ approaches its target within at
  most four time steps (Rempel 2017, Warnecke & Bingert 2020, Navarro et
  al. 2022, Zhou et al. 2025). Unlike `mhd`'s
  [conduction](#eq_mhd_htc), there is no saturation.
- **The p(∇·b̂) term** is included only with **ffhd_pdivb=.true.**
  (compile time, off by default). Its ∇·b̂ is a plain fourth-order central
  difference along each direction, which is correct on a Cartesian mesh
  only: it lacks the metric factors of a divergence in spherical or
  cylindrical coordinates.
- **The gravity, cooling and H_user terms** need **ffhd_gravity**,
  **ffhd_radiative_cooling** and **ffhd_source_usr** respectively. A
  purely local user-defined heating H_user(**x**) is coded as
  `addsource_usr` in `mod_usr.fpp`, see [User source
  terms](#eq_usr_source).

### Setting the frozen field {: #eq_ffhd_bhat }

The unit vector b̂ = **B**/B is a pure function of position: it is never
evolved and is not a boundary condition. The user supplies it by
implementing the by-name hook `usr_set_nwextra` in `mod_usr.fpp` — a
`pure` subroutine, marked as a device routine with `${GPU_ROUTINE_SEQ()}$`
like `gravity_field` above, that takes one point and returns the field's
components there:

```fortran
pure subroutine usr_set_nwextra(x, bhat)
  ${GPU_ROUTINE_SEQ()}$
  double precision, intent(in)  :: x(1:ndim)   ! one cell position
  double precision, intent(out) :: bhat(1:3)   ! b-hat there, in grid components

  ! e.g. a dipole, an arcade, a potential extrapolation ...
  bhat(1:3) = ...
  bhat(1:3) = bhat(1:3) / sqrt(sum(bhat(1:3)**2))   ! MUST be a unit vector
end subroutine usr_set_nwextra
```

Key points:

* **Return a unit vector.** The result is written verbatim into the
  `b1,b2,b3` slots (the `nwextra` variables), so any non-unit field will
  simply be wrong. Normalise explicitly, as above.
* **Use the grid's coordinate components.** In a curvilinear build `bhat`
  is in `(r, θ, φ)` / `(r, z, φ)` components at `x`, not Cartesian. The
  spherical/cylindrical `ffhd` test cases carry a small `to_spherical_unit`
  / `to_cylindrical_unit` helper that rotates a constant Cartesian
  direction into the local basis and normalises in one step.
* **Do not set `b1,b2,b3` anywhere else.** In particular
  `initonegrid_usr` and `usr_special_bc` must leave them alone. AGILE
  evaluates `usr_set_nwextra` once per cell of *every* block — interior,
  inter-block ghosts, physical-boundary ghosts and polar-axis ghosts
  alike — from `alloc_node`, after every grid change, and nothing else
  writes those slots.
* **The polar axis needs no special handling.** In an axis ghost cell the
  stored position carries the mirrored coordinate (a negative `θ` or `r`),
  so evaluating the analytic field there reproduces on its own the sign
  flips a vector picks up across the axis (`b_r` symmetric,
  `b_θ`/`b_z` and `b_φ` antisymmetric).

Note that FFHD re-uses the external gravity module from the hydro case
(see [Gravity](#eq_gravity) above), but here only the effective gravity
projected along the fixed field, **g**·b̂, enters. The radiative loss
treatment (see [Radiative losses](#eq_radloss) above) is shared with `hd`
and `mhd`, activated with **ffhd_radiative_cooling=.true.**.
The equation parameter **ffhd_gamma** sets γ̂, analogous to **hd_gamma**.

## Special relativistic hydrodynamics: srhd {: #eq_srhd }

```fortran
&methodlist
   phys            = 'srhd'
/
```

The special relativistic hydrodynamics module solves the following set
of conservation laws in flat (Minkowski) spacetime

\begin{eqnarray}
    \partial_t D + \nabla\cdot (D \mathbf{v}) &=& 0 \,,\\
    \partial_t \mathbf{m} + \nabla \cdot (\mathbf{v}\mathbf{m} + p\mathbf{I}) &=& 0 \,,\\
    \partial_t \tau + \nabla \cdot\left(\mathbf{m} - D \mathbf{v}\right) &=& 0 \,,
\end{eqnarray}

where the conserved quantities are the density D, the momentum-density
**m** and the rest-mass-subtracted total energy density τ:=e-D, all
evaluated in the lab-frame. The conserved variables are

\begin{equation}
    \mathbf{U} = \left(\begin{array}{c}
        D  \\
        \mathbf{m} \\
        \tau
    \end{array}
    \right)
    = \left(\begin{array}{c}
        \Gamma \rho \\
        \Gamma^2 \rho h \mathbf{v} \\
        \Gamma^2 \rho h - p - D
    \end{array}
    \right) \,,
\end{equation}

with h being the specific enthalpy of the gas and Γ the Lorentz factor
Γ = 1/√(1-v²). The primitive variables are

\begin{equation}
    \mathbf{P} = \left(\begin{array}{c}
        \rho  \\
        \Gamma \mathbf{v} \\
        p
    \end{array}
    \right) \,,
\end{equation}

which correspond to the fluid-frame density ρ, the spatial components of
the four-velocity Γ**v** and the fluid-frame pressure p. Since primitive
variables are required to evaluate the fluxes, a non-linear inversion
**P**(**U**) is needed; AGILE uses the one introduced by Bergmans,
Keppens, van Odyck & Achterberg (2005), which expresses the energy
equation as a 1D root-finding problem in ξ = Γ²ρh. AGILE also saves the
set of auxiliary variables

\begin{equation}
    \mathbf{A} = \left(\begin{array}{c}
        \xi  \\
        \Gamma
    \end{array}
    \right) \,,
\end{equation}

useful for quick conversion between **P** and **U** and to initialize an
initial guess for the root-finder.

As with all Euler-type equations, this set needs to be closed with an
equation of state (EOS) of the form h=h(ρ,p). The two options currently
implemented (selected via the **srhd_eos** logical parameter) are the
*ideal gas* (**srhd_eos=.false.**, the default)

\begin{equation}
    h(\rho,p) = 1 + \frac{\hat{\gamma}}{\hat{\gamma}-1}\frac{p}{\rho} \,,
\end{equation}

and the *Synge gas* (**srhd_eos=.true.**) of the single-species
trans-relativistic gas (Synge 1957), defined as

\begin{equation}
    h(\rho,p) = \frac{K_3(\theta^{-1})}{K_2(\theta^{-1})} \,,
\end{equation}

with the relativistic temperature θ := p/ρ and K_n the modified Bessel
function of the second kind. As in previous work (Meliani & Sauty 2004;
Mignone & Plewa 2005; Ryu, Chattopadhyay & Choi 2006; Porth & Olivares et
al. 2017), AGILE uses an approximation to the previous expression that
avoids the expensive evaluation of the Bessel functions (Mathews 1971;
Goedbloed, Keppens & Poedts 2010). The equation parameter **srhd_gamma**
sets γ̂ for the ideal-gas EOS (and the adiabatic index used in the Synge
gas approximation).

`srhd` has no gravity or radiative cooling. It carries **srhd_n_tracer**
passive tracers and supports [user source terms](#eq_usr_source) with
**srhd_source_usr=.true.** (both compile time). Its density and pressure
floors, **srhd_small_density** and **srhd_small_pressure**, are described
under [Density and pressure floors](#eq_floors), like those of the other
modules.


## Special relativistic magnetohydrodynamics: srmhd {: #eq_srmhd }

```fortran
&methodlist
   phys            = 'srmhd'
/
```

The special relativistic MHD module extends [srhd](#eq_srhd) with a
magnetic field **B**, measured in the lab frame. It solves

\begin{eqnarray}
    \partial_t D + \nabla\cdot (D \mathbf{v}) &=& 0 \,,\\
    \partial_t \mathbf{m} + \nabla \cdot \left(\mathbf{m}\mathbf{v} - \frac{\mathbf{B}\mathbf{B}}{\Gamma^2} - (\mathbf{v}\cdot\mathbf{B})\,\mathbf{v}\mathbf{B} + p_{\mathrm{tot}}\mathbf{I}\right) &=& 0 \,,\\
    \partial_t \tau + \nabla \cdot\left(\mathbf{m} - D \mathbf{v}\right) &=& 0 \,,\\
    \partial_t\mathbf{B}+\nabla\cdot\left[\mathbf{v}\mathbf{B}-\mathbf{B}\mathbf{v} +\psi\mathbf{I}\right] &=& 0 \,,
\end{eqnarray}

with the conserved variables

\begin{equation}
    \mathbf{U} = \left(\begin{array}{c}
        D  \\
        \mathbf{m} \\
        \tau \\
        \mathbf{B}
    \end{array}
    \right)
    = \left(\begin{array}{c}
        \Gamma \rho \\
        (\Gamma^2 \rho h + B^2)\,\mathbf{v} - (\mathbf{v}\cdot\mathbf{B})\,\mathbf{B} \\
        \Gamma^2 \rho h - p - D + \frac{1}{2}B^2(1+v^2) - \frac{1}{2}(\mathbf{v}\cdot\mathbf{B})^2 \\
        \mathbf{B}
    \end{array}
    \right) \,.
\end{equation}

The total pressure adds the magnetic pressure of the *comoving* field,

\begin{equation}
    p_{\mathrm{tot}} = p + \frac{1}{2}b^2 \,, \qquad
    b^2 = \frac{B^2}{\Gamma^2} + (\mathbf{v}\cdot\mathbf{B})^2 \,.
\end{equation}

The momentum flux is the spatial part of the covariant stress
\( (\rho h + b^2)\,u^i u^j - b^i b^j + p_{\mathrm{tot}}\delta^{ij} \), with
\( b^i = B^i/\Gamma + \Gamma(\mathbf{v}\cdot\mathbf{B})\,v^i \), written in lab-frame
quantities. When it is expanded, the \( \Gamma^2(\mathbf{v}\cdot\mathbf{B})^2 v^iv^j \)
pieces cancel. The coefficient of \( v^iv^j \) is therefore
\( \Gamma^2\rho h + B^2 \), with the lab-frame \( B^2 \). The comoving
\( b^2 \) appears only in \( p_{\mathrm{tot}} \).

The primitive variables are

\begin{equation}
    \mathbf{P} = \left(\begin{array}{c}
        \rho  \\
        \Gamma \mathbf{v} \\
        p \\
        \mathbf{B}
    \end{array}
    \right) \,,
\end{equation}

As in `srhd`, the momentum slot holds the spatial four-velocity Γ**v**. The
inversion **P**(**U**) follows Kastaun, Kalinani & Ciolfi (2021), which
reduces it to a 1D root-finding problem in μ = 1/(Γh). A Brent solver finds
the root, bracketed by a second, auxiliary root-find. The
module stores the auxiliary variables

\begin{equation}
    \mathbf{A} = \left(\begin{array}{c}
        \mu  \\
        \Gamma
    \end{array}
    \right) \,,
\end{equation}

The stored μ is the initial guess for the next inversion. The auxiliaries
are refreshed wherever the conserved state is written: after every step,
after the initial condition, and in the kernels that create a block during a
regrid.

The divergence of **B** is controlled with the same GLM field ψ as in
`mhd`, including its exactly integrated damping term; see [Divergence B
source treatments](#eq_divB_fix). ψ is the additional conserved variable.

Parameters are read from the **srmhd_list** of the parameter file:

* **srmhd_gamma** sets γ̂ of the ideal-gas EOS h = 1 + γ̂/(γ̂−1)·p/ρ. That
  EOS is the only one available: **srmhd_eos=.true.** (the Synge gas)
  stops the run.
* **srmhd_glm_alpha** is α_GLM. It defaults to 0.5.
* **srmhd_small_pressure** and **srmhd_small_density** are the floors that
  `fix_prim_state` applies to a reconstructed face state. They default to
  zero and 10⁻¹⁴ respectively. Unlike in `srhd`, the inversion itself does
  not read them.

On a spherical or cylindrical mesh, `srmhd` adds the geometric source terms
given in [Coordinate systems](coordinate_systems.html#geom_spherical_srmhd).


# Positivity fixes {: #eq_positivity_fixes }

A finite-volume update can drive the density or the pressure of a cell
negative, typically in strong rarefactions, near high-Mach shocks, or where
radiative cooling removes most of the thermal energy within a step. AGILE
handles this differently from upstream MPI-AMRVAC, and more narrowly:

- **There is no general small-values machinery.** Upstream's `small_density`,
  `small_pressure`, `check_small_values`, `fix_small_values` and
  `small_values_method` are not part of AGILE, and the cell averages of the
  evolved solution are never floored by the core code.
- **Only the reconstructed face states are floored.** Every physics module
  floors the density and pressure of the face states it hands to the
  Riemann solver (see [Density and pressure floors](#eq_floors) below). A
  cell average that the update itself drives negative is not reset, and its
  negative pressure then goes into the sound speed (`sqrt(gamma*p/rho)`)
  when the next time step is computed, which yields NaN.
- **A run that goes unphysical is not stopped.** `setdt` raises a `crash`
  flag when the time step falls below `dtmin`, but the check in
  `timeintegration` that would save a snapshot and abort on it is currently
  commented out. Watch the log for NaN.

What there is: a hook on the reconstructed face states, which every module
uses for a density and pressure floor; and the temperature floor of
radiative cooling.

## The face-state hook: fix_prim_state {: #eq_fix_prim_state }

After the MUSCL reconstruction, and before the Riemann solver, every
reconstructed primitive face state passes through

```fortran
pure subroutine fix_prim_state(u)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(inout) :: u(nw_phys)   ! one face state, in primitive variables
```

It is called from `muscl_reconstruct_prim` in `src/mod_finite_volume.fpp`
for the left and right state at both faces of a cell, in every direction,
whatever the limiter. That routine is the single producer of face states, so
this one call covers all three Riemann solvers (`tvdlf`, `hll`, `hllc`). The
reconstruction limits every slot of the state independently, so a face
state can come out with a negative density or pressure even where both
neighbouring cells are admissible. The hook is the place to repair it.

Two things it is *not*:

- **It does not touch the cell averages.** It changes only the states the
  fluxes are computed from. A cell can still go negative through the update
  itself, and nothing then resets it.
- **It is not a user hook.** `fix_prim_state` is a fypp macro of the physics
  module, not a routine in `mod_usr.fpp`. `src/physics/mod_physics_dummies.fpp`
  defines it empty, which the compiler inlines away, and a physics template
  that defines its own `#:def fix_prim_state()` replaces it. All four
  modules define one.

## Density and pressure floors {: #eq_floors }

Each module's `fix_prim_state` floors the density and pressure of a face
state, at two run-time parameters of its `_list`:

| Module | Parameters |
| --- | --- |
| `hd`   | `hd_small_density`, `hd_small_pressure` |
| `mhd`  | `mhd_small_density`, `mhd_small_pressure` |
| `ffhd` | `ffhd_small_density`, `ffhd_small_pressure` |
| `srhd` | `srhd_small_density`, `srhd_small_pressure` |

for example

```fortran
&hd_list
   hd_small_density  = 1.0d-10
   hd_small_pressure = 1.0d-12
/
```

The face state is floored as ρ ← max(ρ, `*_small_density`) and
p ← max(p, `*_small_pressure`). Both parameters default to zero, and a
negative value stops the run at startup. With the defaults the `max()` only
clips values that were already negative, and leaves every admissible state
bit-for-bit unchanged, so a run that never produces a negative face state
gives exactly the same result with or without the floors. A pressure clipped
to exactly zero still gives a zero sound speed, and a density clipped to
zero a division by zero in the velocity; set both to small positive values
when floors are actually wanted. How `max()` treats a NaN is left to the
compiler, so a face state that is already NaN is not guaranteed to survive
the floor as NaN.

These floors act on the face states only, so they keep the fluxes
computable but do not stop a cell average from going negative.

### srhd {: #eq_srhd_floors }

`srhd`'s hook does one more thing, and always: the reconstruction limits the
auxiliary variables ξ and Γ like independent variables, although they are
functions of (ρ, Γ**v**, p). A reconstructed face state therefore does not
satisfy Γ = √(1+|Γ**v**|²) or ξ = Γ²ρh in general, and the signal speeds
computed from it can exceed the speed of light. After the floors,
`fix_prim_state` recomputes both from the reconstructed primitives, with the
same expressions as `to_conservative` (the Synge-gas enthalpy included).

`srhd_small_pressure` is also the lower bound of the pressure bracket in the
primitive-variable inversion (`mod_con2prim`), so in `srhd` it constrains
the pressures recovered from the cell averages too.

### Writing the hook for a new physics module {: #eq_floor_howto }

A new physics module gets the empty default until its template defines the
macro, inside its `#:if PHYS == '...'` block. At this point `u` is in
primitive variables; for `hd`, for example, the density is in `u(iw_rho)`
and the pressure in `u(iw_e)`, and the whole hook is

```fortran
#:def fix_prim_state()
  pure subroutine fix_prim_state(u)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: u(nw_phys)

    u(iw_rho) = max(u(iw_rho), hd_small_density)
    u(iw_e)   = max(u(iw_e),   hd_small_pressure)
  end subroutine fix_prim_state
#:enddef
```

with the two floors declared next to the module's other parameters
(`${GPU_DECLARE_COPYIN(...)}$`), added to its namelist, checked for
negative values, and pushed to the device after reading
(`${GPU_UPDATE_DEVICE(...)}$`). If the module carries auxiliary variables
derived from the primitives, as `srhd` does, recompute them there as well.

## Temperature floor of radiative cooling {: #eq_cooling_floor }

The one floor that acts on the evolved solution belongs to radiative
cooling. With `Tfix = .true.` in `&rc_list`, after the cooling source has
been added, a cell is given back the thermal energy it lacks to reach
ρ𝓡T_low, where `tlow` sets T_low (in code units; it defaults to the lowest
temperature of the cooling table). The test is made on the state the stage
started from, not on the cooled state, so a single stage can still cool a
cell slightly below T_low. It only applies
where `hd_radiative_cooling`, `mhd_radiative_cooling` or
`ffhd_radiative_cooling` is on, and only to the cooling term: other source
terms can still lower the temperature below T_low. See [Radiative
cooling](radiative_cooling.html).
