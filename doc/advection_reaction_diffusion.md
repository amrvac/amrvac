# Advection-reaction-diffusion module

[TOC]

This page describes the current (22/06/2022) functionality of the advection-reaction-diffusion (ard) module.
This module extends the reaction-diffusion module with the possibility of adding an advection term.
The source code can be found in `src/ard/mod_ard_phys.t`.

# Physics {#label-1}
Systems of \f$ s \f$ advection-reaction-diffusion equations are of the form

\f[ \frac{\partial u_i}{\partial t} + \nabla(a_i(u_1, \ldots, u_s) u_i) = \nabla(D_i \nabla u_i) + f_i(u_1, \ldots, u_s) \qquad i \in \{1,\ldots, s\} \f]

where \f$ a_i \f$ are the advection coefficients which potentially depend on the solution values themselves. \f$ D_i \f$ are the diffusion coefficients and \f$ f_i \f$ the reaction functions.
In the ard module in MPI-AMRVAC, several such systems are available with an advection term of the form:

\f[ \nabla(a_i(u_1, \ldots, u_s) u_i) = \nabla(\frac{A_i}{p} u_i^p) = A_i u_i^{p-1} \nabla u_i, \f]

for some integer power \f$ p \f$.
The systems of interest are given below.
The prefixes relate to the names given to the systems in the code, to be supplied in the `.par` file or prefixing the variables.
Some systems have a good starting resource linked to them.

- `nr` = `no_reac`: **No reaction term** or the pure advection-diffusion equation

    \f[
        \frac{\partial u}{\partial t} + A u^{p-1} \nabla u = D \nabla^2 u
    \f]

- `lg` = `logistic`: [**Fisher-KPP**](https://people.maths.ox.ac.uk/trefethen/pdectb/fisher2.pdf) or the diffusive logistic equation

    \f[
        \frac{\partial u}{\partial t} + A u^{p-1} \nabla u = D \nabla^2 u + \lambda u(1-u) 
    \f]

- `gs` = `gray-scott`: [**Gray-Scott**](http://mrob.com/pub/comp/xmorphia/index.html)

    \f{eqnarray*}{
        \frac{\partial u}{\partial t} + A_1 u^{p-1} \nabla u &=& D_1 \nabla^2 u - uv^2 + F(1-u) \\
        \frac{\partial v}{\partial t} + A_2 v^{p-1} \nabla v &=& D_2 \nabla^2 v + uv^2 - (F+k)v 
    \f}

- `br` = `brusselator`: **Brusselator**

    \f{eqnarray*}{
        \frac{\partial u}{\partial t} + A_1 u^{p-1} \nabla u &=& D_1 \nabla^2 u + u^2 v + A - (B+1) u \\
        \frac{\partial v}{\partial t} + A_2 v^{p-1} \nabla v &=& D_2 \nabla^2 v - u^2 v + Bu 
    \f}

- `sb` = `schnakenberg`: [**Schnakenberg**](https://cbeentjes.github.io/files/Ramblings/PatternFormationSchnakenberg.pdf)

    \f{eqnarray*}{
        \frac{\partial u}{\partial t} + A_1 u^{p-1} \nabla u &=& D_1 \nabla^2 u + \gamma (\alpha + u^2 v - u) \\
        \frac{\partial v}{\partial t} + A_2 v^{p-1} \nabla v &=& D_2 \nabla^2 v + \gamma (\beta - u^2 v) 
    \f}

- `ebr` = `ext-brusselator`: **Extended Brusselator** (where \f$ f(u,v) \f$ and \f$ g(u,v) \f$ are the reaction functions in the Brusselator model (see  above))

    \f{eqnarray*}{
        \frac{\partial u}{\partial t} + A_1 u^{p-1} \nabla u &=& D_1 \nabla^2 u + f(u,v) - Cu + Dw \\
        \frac{\partial v}{\partial t} + A_2 v^{p-1} \nabla v &=& D_2 \nabla^2 v + g(u,v) \\
        \frac{\partial w}{\partial t} + A_3 w^{p-1} \nabla w &=& D_3 \nabla^2 w + Cu - Dw
    \f}

- `bzfn` = `belousov_fieldnoyes`: **Oregonator** (Field-Noyes model of the Belousov-Zhabotinski reaction)

    \f{eqnarray*}{
        \frac{\partial u}{\partial t} + A_1 u^{p-1} \nabla u &=& D_1 \nabla^2 u + \frac{1}{\epsilon} (\lambda u - u w + u - u^2) \\
        \frac{\partial v}{\partial t} + A_2 v^{p-1} \nabla v &=& D_2 \nabla^2 v + u - v \\
        \frac{\partial w}{\partial t} + A_3 w^{p-1} \nabla w &=& D_3 \nabla^2 w + \frac{1}{\delta} (-\lambda w - u w + \mu v))
    \f}

- `lor` = `lorenz`: **The Lorenz system**

    \f{eqnarray*}{
        \frac{\partial u}{\partial t} + A_1 u^{p-1} \nabla u &=& D_1 \nabla^2 u + \sigma (v-u) \\
        \frac{\partial v}{\partial t} + A_2 v^{p-1} \nabla v &=& D_2 \nabla^2 v + u(r - w) - v \\
        \frac{\partial w}{\partial t} + A_3 w^{p-1} \nabla w &=& D_3 \nabla^2 w + uv - bw 
    \f}


# Numerics {#label-2}
Unlike for the reaction-diffusion module, the flux schemes and slope limiters are applicable to this advection-reaction-diffusion module.
There is a possibility to use IMEX schemes, where the diffusion part is handled implicitly and reaction explicitly.
The implicit system is then solved using the AMR-compatible multigrid solver coupled to MPI-AMRVAC.

    <!-- J. Teunissen, R. Keppens,
    A geometric multigrid library for quadtree/octree AMR grids coupled to MPI-AMRVAC,
    Computer Physics Communications,
    Volume 245,
    2019,
    106866,
    ISSN 0010-4655,
    https://doi.org/10.1016/j.cpc.2019.106866. -->


# Practical usage {#label-3}
A number of examples can be found in `tests/ard/`.
The namelist for the advection-reaction-diffusion module is as follows:

    &rd_list
        D1, D2, D3,                  ! Diffusion coefficients
        adv_pow,                     ! Power p of the unknown in the advection term
        A1, A2, A3,                  ! Advection coefficients (d-dimensional)
        sb_alpha, sb_beta, sb_kappa, ! Parameters for Schnakenberg
        gs_F, gs_k,                  ! Parameters for Gray-Scott
        br_A, br_B,                  ! Parameters for Brusselator
        br_C, br_D,                  ! Parameters for extended Brusselator
        lg_lambda,                   ! Parameter for Fisher-KPP equation
        bzfn_epsilon, bzfn_delta, bzfn_lambda, bzfn_mu, ! Parameters for Oregonator model
        lor_r, lor_sigma, lor_b,     ! Parameters for Lorenz system
        equation_name,               ! Name of the system to simulate
        dtreacpar,                   ! Timestep restriction parameter
        ard_particles, ard_source_split 
    /
