# Radiative transfer with short characteristics in MPI-AMRVAC

[TOC]

Optically thin cooling describes radiation that escapes freely, which fails in
the solar photosphere, in the dense chromosphere and in the cores of prominence
threads, where the gas absorbs much of what it emits. For these layers
MPI-AMRVAC solves the radiative transfer equation in gray LTE with the method of
short characteristics (SC), and adds the resulting radiative heating or cooling
to the energy equation. The module works in 2D Cartesian geometry, on a uniform
grid within the layers where it is active. This page describes the physics, how
to use the module, and how it is organised. Set the
[normalization units](amrvac_units.md) in `usr_init()`, since the opacities and
the radiation constants are converted with them.

# Physics {#rt-physics}

Along every ray direction \f$\boldsymbol\mu\f$ the specific intensity \f$I\f$
obeys
\f[
 \boldsymbol\mu\cdot\nabla I=\chi\,(B-I),
\f]
where \f$\chi=\kappa\rho\f$ is the opacity per unit length and the source function
is the Planck function, written as an energy density \f$B=aT^4\f$. Scattering is
not included, so for a given state the intensity of each ray follows from one
sweep across the grid. The rays sample the unit sphere with a fixed quadrature,
and their weighted sums give the mean intensity \f$J\f$ and the radiative flux
\f$\boldsymbol F=c\sum_k w_k\boldsymbol\mu_k I_k\f$.

The heating per unit volume is
\f[
 Q=c\,\chi\,(J-B).
\f]
Near the top of a layer this form is accurate. Deep down, where \f$J\f$ and
\f$B\f$ become nearly equal, their difference is dominated by round-off, and the
equivalent conservative form \f$Q=-\nabla\cdot\boldsymbol F\f$ is better. The
module therefore blends the two,
\f[
 Q=w\,c\,\chi\,(J-B)+(1-w)(-\nabla\cdot\boldsymbol F),
 \qquad w=e^{-k\tau},
\f]
with \f$\tau\f$ the optical depth measured from the top of the layer and \f$k\f$
set by `rt_sc_blend_k`. Where the column above a cell is optically thin, the
heating is tapered away by the factor \f$\tau^2/(\tau^2+\tau_{\min}^2)\f$, since
the optically thin cooling describes those cells.

## Opacity and the hand-over to thin cooling {#rt-opacity}

The opacity is one of three models chosen by `rt_sc_opacity`: a tabulated
Rosseland mean, an analytic fit of the H-minus and Kramers opacities with
electron scattering, or a routine supplied by the user. The two built-in models
switch the transfer off smoothly in the hot, tenuous corona. A cell keeps its
transfer when it is cooler than `rt_sc_Tcap` or denser than the pressure
`rt_sc_p_TR`, and outside both \f$\chi\f$ and \f$B\f$ are scaled down together, so
that the corona becomes both transparent and dark to the transfer.

Where the transfer is active, the radiative losses are divided between the SC
heating and the optically thin cooling, so that no loss is counted twice. The
share of the thin cooling is
\f[
 w_{\rm thin}=\tfrac12\left[1+\tanh\frac{T-T_0}{\Delta T}\right]
 \max\left[0,\,1-(p/p_{\rm TR})^2\right]
\f]
inside an RT layer and one outside, and the SC heating is multiplied by
\f$1-w_{\rm thin}\f$. The radiative cooling module applies this weight to its own
losses automatically.

# Practical use {#rt-practical}

## Activating the module {#rt-activation}

The transfer is switched on in the physics namelist, in the same way as radiative cooling,

    &hd_list
       hd_radiative_transfer = .true.
    /

or with `mhd_radiative_transfer = .true.` in `&mhd_list`. Both need the energy equation. The
settings of the transfer are read from `&rt_list`, described below, and the module then sets
itself up. It registers the output variables `Qrad`, the heating, and `Frad`, the vertical
flux, which are written to the data files with the other variables. It solves the transfer
once per time step, before the time integration, and adds the heating to the energy equation
at every stage, limited so that one stage changes the gas energy by at most `rt_sc_cfl` of
itself. When radiative cooling is also active, the heating takes only its share
\f$1-w_{\rm thin}\f$, and the refinement routines hold every block of a band on one AMR
level. Nothing needs to be added to the user file.

Two hooks remain for setups that need them. With `rt_sc_opacity = 'usr'`, the user points
`usr_rt_sc_chiB` at a routine in `usr_init()` that returns \f$\chi\f$ and \f$B\f$ on a
block, as the test `tests/hd/sc_rt_2D` does. A setup that controls its energy input from the
radiative output can call `rt_sc_emergent_tau1`, described below.

## Parameters {#rt-parameters}

The parameters are read from `&rt_list`. Heights are in cm, pressures in
dyn/cm\f$^2\f$ and temperatures in K, and they are converted to code units when
the module is set up.

name | type | default | description
---|---|---|---
rt_sc_updim | integer | last axis | The vertical axis, 1 or 2.
rt_sc_amr_level | integer | 1 | The AMR level that the bands are held on. It cannot exceed `refine_max_level`.
rt_sc_nband | integer | 1 | Number of bands, at most 8. A bar with a photosphere at each end needs two.
rt_sc_band_min, rt_sc_band_max | double precision arrays | -1e30, 1e30 | Bottom and top of each band in cm. The default values stand for the domain ends. The top is moved down to the nearest block face, and the run stops if a band is thinner than one block.
rt_sc_bc_min, rt_sc_bc_max | string arrays | 'thermal', 'open' | Incoming intensity at the bottom and top of each band. 'thermal' sets \f$I=B\f$ of the edge cell and 'open' sets \f$I=0\f$.
rt_sc_bc_hmin, rt_sc_bc_hmax | string | 'open' | Incoming intensity at the horizontal edges when the domain is not periodic.
rt_sc_opacity | string | 'Rosseland' | Opacity model, 'Rosseland', 'Kramers' or 'usr'. With 'usr' the user points `usr_rt_sc_chiB` at a routine that returns \f$\chi\f$ and \f$B\f$.
rt_sc_kappa_table | string | 'rosseland_combined.tab' | File of the Rosseland table.
rt_sc_Tcap, rt_sc_dTcap | double precision | 2e4, 5e3 | Temperature and width in K of the coronal switch of the built-in opacities.
rt_sc_p_TR, rt_sc_dp_TR | double precision | 1, 0.5 | Pressure in dyn/cm\f$^2\f$ and width in dex of the same switch.
rt_sc_quad | string | 'GL' | Angular quadrature. 'GL' takes Gauss-Legendre polar nodes times equally spaced azimuths, and 'Carlson' the level-symmetric A4 set of 24 directions, which a 2D sweep merges into 12.
rt_sc_nmu, rt_sc_ndir | integer | 4, 6 | Polar nodes (2, 4, 6 or 8) and azimuths of the GL set.
rt_sc_flux_blend | logical | T | Blend \f$c\chi(J-B)\f$ with \f$-\nabla\cdot\boldsymbol F\f$. With false the heating is \f$c\chi(J-B)\f$ everywhere.
rt_sc_blend_k | double precision | 10 | Rate \f$k\f$ of the blend weight \f$e^{-k\tau}\f$.
rt_sc_tau_min | double precision | 1e-6 | Optical depth of the thin-column taper.
rt_sc_cfl | double precision | 0.5 | Largest fraction of the gas energy that the heating may change in one stage.
rt_sc_blend_T0, rt_sc_blend_dT | double precision | 2e4, 5e3 | Centre and width in K of the temperature gate of the thin-cooling share.
rt_sc_distributed | logical | T | Use the distributed wavefront solvers. With false every rank holds the whole band; this is the original implementation, kept for comparison only and not for production runs.
rt_sc_pipeline | logical | T | Pipeline the wavefront, so that rays of different directions travel through it together.
rt_sc_seam_on | logical | T | Carry the intensity of a ray across a periodic horizontal edge, using the outflow of the previous solve.
rt_sc_valpipe | logical | F | For validation, run the ray-by-ray and the pipelined wavefront on the same state and print their difference.
rt_sc_debug | logical | F | Print the cumulative RT timing every 10 solves.

A single photospheric band from the bottom of the domain to 5 Mm, held on level 3
with the tabulated opacity and the Carlson set, uses

    &rt_list
       rt_sc_amr_level = 3
       rt_sc_band_max = 5.0d8
       rt_sc_bc_min = 'thermal'
       rt_sc_bc_max = 'open'
       rt_sc_quad = 'Carlson'
       rt_sc_opacity = 'Rosseland'
       rt_sc_kappa_table = 'rosseland_combined.tab'
    /

## Emergent flux {#rt-emergent}

A setup that controls its energy input from the radiative output, for example
through the entropy of the inflow at the bottom of a convection simulation,
needs the flux that leaves the photosphere. `rt_sc_emergent_tau1` computes the
horizontal mean of the vertical flux on the surface \f$\tau=1\f$ of the first band
and stores it in `rt_sc_Ftau1`. The flux at the top of the band itself is not
suitable, because the radiative flux keeps growing above \f$\tau=1\f$.

# Numerics and implementation {#rt-implementation}

**Transport.** Each ray is integrated with the linear short characteristics of
Olson and Kunasz. A cell takes the intensity at the upwind point on the face the
ray crosses first, interpolated linearly along that face together with
\f$\chi\f$ and \f$B\f$, and integrates across one step with the source varying
linearly in optical depth. The rays are swept cell by cell from the upwind
corner, so every cell is visited once per direction.

**Solvers.** The default solver is a pipelined wavefront (KBA). Each rank sweeps
only its own blocks of a band, one diagonal of blocks at a time, and passes the
intensity on through the ghost cells, so no rank ever holds the whole band. The
rays that share an upwind corner follow one another through the wavefront, and
the corners run together. The ray-by-ray wavefront (`rt_sc_pipeline=.false.`)
sends one ray at a time and gives the same result bit for bit. The replicated
solver (`rt_sc_distributed=.false.`) is the original implementation, kept for
comparison and not for production runs. It assembles the whole band on every
rank, and its periodic wrap is not exact: a ray swept column by column reaches
the seam before its upwind column holds that ray, so the seam takes a vacuum
inflow for it. The wavefront carries the intensity across the seam from the
previous solve instead, which is exact at steady state.

**Bands and AMR.** The transfer assumes uniform cells, so the refinement routines
hold every block that overlaps a band on `rt_sc_amr_level`, refining coarser
blocks and coarsening finer ones. Outside the bands the mesh refines freely, and
the heating there is zero.

**Modules.** The code is in `src/physics`. `mod_rt_sc_slab` is the module that
user files and the rest of the code use. It holds the coupling to the gas, the
heating source, the thin-cooling share, the AMR hold and the emergent flux.
`mod_rt_sc_state` reads and converts the parameters and holds the shared state,
`mod_rt_sc_solvers` the three solvers, `mod_rt_sc_kernel` the quadrature, the
sweep and the heating formula, and `mod_rt_sc_opacity` the opacity models.

**Limitations.** The module is 2D and gray, with LTE source functions and no
scattering. It assumes uniform Cartesian cells: other coordinate systems are
refused, stretched grids are not supported, and each band lies on a single AMR
level. In the wavefront the corner ghost of a block can hold the next ray, so the
corner value is rebuilt from the two face neighbours. In the slab test this costs
up to 8.5% of the peak heating in the corner cells of a block.

**Tests.** `tests/hd/sc_rt_2D` checks the transport of a beam across the periodic
seam, the heating blend and taper on an absorbing slab, and the Rosseland
opacity on a stratified atmosphere, with the replicated, ray-by-ray and pipelined
solvers.
