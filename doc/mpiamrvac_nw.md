# Auxiliary variables (nw...)

[TOC]

# List of auxiliary variables {#list}

This document describes the differences and the intended use of the **nwaux,
nwextra** parameters for MPIAMRVAC. These global parameters are dependent on
the physics module, and defined in the corresponding

      amrvacpar.t.physics

file, where _physics_ is currently one of _rho, hdadiab, hd, mhd, srhd,
srhdeos, srmhd, srmhdeos_.

## nw {#nw}

This global parameter sets the total number of variables, determining the
(last dimension-independent) size of the _w(ixG^T,1:nw)_ array, and is always

     nw = nwflux+nwaux+nwextra

The idea is that there are _nwflux_ conservative variables for which
corresponding fluxes are defined in the _amrvacphys.t.physics_ module, and
that we have in addition a set of _nwaux_ auxiliary variables, together with
some extra _nwextra_ variables. For these latter two types of variables, one
does not have a corresponding flux definition in _getflux, gefluxforhllc_
subroutines, and they also have no boundary conditions imposed on them.

## nwflux {#nwflux}

The first _1:nwflux_ variables should be the conservative variables of the
physics module at hand, and these are the ones that are updated by means of
fluxes across cell boundaries, and are influenced by source terms (geometric
or real physical source terms). Also, those variables change their
meaning/content in the subroutines _primitive, conserve_. These latter
subroutines are thus switching from conservative to primitive and vice versa
at appropriate places in the code, and in this switching procedure, handy use
can be made of the possible auxiliary variables. By default, only the special
relativistic physics modules _srhd, srhdeos, srmhd, srmhdeos_ have such
auxiliary variables defined for them, i.e. for these modules we have
_nwaux>0_, and auxiliary variables are then in the slots
_nwflux+1:nwflux+nwaux_. For the first _1:nwflux_ variables, the boundary
conditions are to be imposed and existing boundary types have to be provided
(like _symm, asymm, special, cont, periodic, ..._) in the par-file. For
steady-state computations, the residual is based on only the temporal
variation of these _1:nwflux_ variables.

## nwaux {#nwaux}

When for the physics module at hand, one has _nwaux>0_, the slots in the _w_-
array corresponding to _nwflux+1:nwflux+nwaux_ contain so-called auxiliary
variables. One should distinguish between two types of such auxiliary
variables, namely those that are **local** and those that are **nonlocal**.

## Local auxiliaries {#local}

Local auxiliaries can be computed directly from instantaneous local values of
the first _1:nwflux_ variables.

In the physics modules for _srhd, srhdeos, srmhd, srmhdeos_ (special
relativistic modules), these local auxiliaries are useful to store pressure or
Lorentz factor values, which in principle are computable from the set of
conservative variables at any time. As a costly non-linear iterative (Newton-
Raphson) procedure is involved, it pays off in computing time to have
corresponding auxiliary variables available when e.g. fluxes have to be
defined/added in _getflux_ and so on. These local auxiliary variables are all
to be computed in the _getaux_ subroutine in the physics module. When
_nwaux>0_, auxiliaries are computed by means of a _getaux_ call at various
places in e.g. _tvdlf.t_. Since they are to be used for facilitating flux
specifications, they also need to be known in the ghost cells. This is taken
care of by the _fix_auxiliary_ call in the boundary condition treatment: note
that they are thus computed from local (ghostcell) values from conservative
variables (themselves computed according to their boundary type), so that we
do not need to specify a boundary type for the _nwflux+1:nwflux+nwaux_
variables. In practice, since boundary types are defined for _nw_ variables
per side, just enter e.g. the string _'dummy'_ for the auxiliaries. The local
auxiliaries can in fact also be used to trigger refinement or coarsening, so
that the corresponding _flags_ integer array may use them.

## Nonlocal auxiliaries {#nonlocal}

Nonlocal auxiliaries contain gradients of the first
_1:nwflux_ variables.

We added the possibility to compute, in every timestep, and immediately prior
to saving a snapshot (and following the _setdt_ subroutine), local variables
that depend on gradients of the first variables. These can store things like
divergence of the velocity, curl of the magnetic field (current), which can be
useful for visualization later on. These nonlocal auxiliary variables should
not be computed in the _getaux_ subroutine, since the geometric info is not
available there. Similarly to the local auxiliaries, they are not advanced,
there are no boundary conditions imposed on them (in fact, these variables can
typically not be computed meaningfully in the ghost cells due to stencil
restrictions), and they need to be computed/specified in the
_process_grid_usr_ subroutine, for which the default interface is provided in
the _amrvacnul.speciallog.t_ module. These nonlocal variables should not be
used for facilitating flux computations, nor are they very useful for
triggering refinement (as they will not be truly up-to-date at all instances
in a multi-step timestep procedure, because they are not computed in the
_getaux_ subroutine). These kind of nonlocal auxiliaries may come in handy for
particle acceleration treatments that do not react back on the flow dynamics,
which involves the seperate evolution of particle distributions according to
the velocity field and local compression.

## nwextra {#nwextra}

Finally, an extra set of _nwextra_ variables can be defined in addition to the
presently available ones in physics modules like _hd, mhd_ by setting
_nwextra>0_. Like auxiliaries, they are not advanced (no fluxes need to be
defined for them in _getflux_ and so), there are no boundary conditions
imposed on them (again, enter corresponding strings _'dummy'_ in their
boundary setting in the par file), and they should not be used for basing
refinement/derefinement on. They are not computed in the ghost cells, and they
can be used to store additional physical variables like old pressure or
temperature values, the accumulated luminosity (internal energy change) during
a timestep, etc. They were added specifically for facilitating the post-
processing of simulations involving (optically thin) radiative losses.

## nwauxio {#nwauxio}

Only for post-processing purposes on saved snapshots from the code, one may
want to compute additional auxiliary variables for visualization. Hence, only
at the [convert-stage](convert.html) (from amrvac **.dat** data file to any of
the available formats in the _convert_ subroutine) do we need to enlarge the
data size from _nw_ to _nw+nwauxio_, and their calculation is to be done in
the _specialvar_output_ subroutine, for which the default interface is given
in the _amrvacnul.speciallog.t_ module. Correspondingly, the names have to be
given in the _specialvarnames_output_ as strings to be concatenated with the
_wnames/primnames_ strings. This means that normally _nwauxio=0_, but it can
be at convert stage set to a nonzero value in the _filelist_ part of the par-
file.
