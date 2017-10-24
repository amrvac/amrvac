# Auxiliary variables (nw...)

[TOC]

# List of auxiliary variables {#list}

This document describes the differences and the intended use of the **nwaux,
nwextra** parameters. These global parameters are defined in the `mod_variables.t`

## nw {#nw}

This global parameter sets the total number of variables, determining the
(last dimension-independent) size of the _w(ixG^T,1:nw)_ array, and is always

     nw = nwflux+nwaux+nwextra

The idea is that there are _nwflux_ conservative variables for which
corresponding fluxes are defined in the physics modules, and
that we have in addition a set of _nwaux_ auxiliary variables, together with
some extra _nwextra_ variables. For these latter two types of variables, one
does not have a corresponding flux definition in _getflux, gefluxforhllc_
subroutines, and they also have no boundary conditions imposed on them.

## nwflux {#nwflux}

The first _1:nwflux_ variables should be the conservative variables of the
physics module at hand, and these are the ones that are updated by means of
fluxes across cell boundaries, and are influenced by source terms (geometric
or real physical source terms). Also, those variables change their
meaning/content in the subroutines _phys_to_primitive_, _phys_to_conserved_. 
These two subroutines are thus switching from conservative to primitive and vice versa
at appropriate places in the code. For the first _1:nwflux_ variables, the boundary
conditions are to be imposed and existing boundary types have to be provided
(like _symm, asymm, special, cont, periodic, ..._) in the par-file. For
steady-state computations, the residual is based on only the temporal
variation of these _1:nwflux_ variables.

## nwaux {#nwaux}

When for the physics module at hand, one has _nwaux>0_, the slots in the _w_-
array corresponding to _nwflux+1:nwflux+nwaux_ contain so-called auxiliary
variables. They can be computed directly from instantaneous local values of
the flux _1:nwflux_ variables. 
These local auxiliary variables are all
to be computed in the _getaux_ subroutine in the physics module.
Since they are to be used for facilitating flux specifications, they also 
need to be known in the ghost cells. This is taken
care of by the _fix_auxiliary_ call in the boundary condition treatment: note
that they are thus computed from local (ghostcell) values from conservative
variables (themselves computed according to their boundary type), so that we
do not need to specify a boundary type for the _nwflux+1:nwflux+nwaux_
variables.  The local auxiliaries can in fact also be used to trigger 
refinement or coarsening. The current (Newtonian) physics modules of AMRVAC have no real use for these auxiliary variables.

## nwextra {#nwextra}

An extra set of _nwextra_ variables can be defined in addition to the
presently available ones in physics modules by using the function 
var_set_extravar (from mod_variables.t) in usr_init of mod_usr.t. Like 
auxiliaries, they are not advanced (no fluxes need to be
defined for them in _get_flux_), there are no boundary conditions
imposed on them, and they should not be used for basing
refinement/derefinement on. They are not computed in the ghost cells, and they
can be used to store additional physical variables like old pressure or
temperature values, the accumulated luminosity (internal energy change) during
a timestep, local variables that depend on gradients of the flux variables 
(divergence of the velocity, curl of the magnetic field), etc. 
they need to be computed/updated in the pointed _usr_process_grid_ subroutine.
They work to facilitate the understanding of simulation processes.

## nwauxio {#nwauxio}

Only for post-processing purposes on saved snapshots from the code, one may
want to compute additional auxiliary variables for visualization. Hence, only
at the [convert-stage](convert.md) (from .dat data file to any of the available 
formats in the _convert_ subroutine in src/amrvacio/convert.t file) do we 
need to enlarge the data size from _nw_ to _nw+nwauxio_, and their calculation is to be done in
the pointed _usr_aux_output_ subroutine, for which the default interface is given
in the mod_usr_methods.t module. Correspondingly, the names have to be
given in the pointed _usr_add_aux_names_ subroutine strings.
This means that normally _nwauxio=0_, but it can be at convert stage set to 
a nonzero value in the _filelist_ part of the par- file.
