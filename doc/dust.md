# Dust in MPI-AMRVAC

# MPI-AMRVAC dust implementation

This is an overview of the option that allows users to model dust+gas
flows in mpi-amrvac.

The source code can be found in `src/mod_dust.t`

# Physics

Dust in mpi-amrvac is treated as a pressureless fluid that coexists
with the actual gas. Dust and gas are coupled by a dragforce. It is
calculated based on the physical characteristics of gas and dust
grains according to Kwok 1975, page 584 (between equations 8 and 9).

The sticking coefficient can be set according to dedicated parameters.

It is possible to combine multiple dust species in a single
simulation, although there is no direct interaction (collisions,
fractionation, coagulation) between them.

Because single-fluid collisions are not taken into account (i.e. dusty
fluids are pressureless), a shock happening in a dusty flow will
result in unphysical sharp features as dust is brought to a stop.  For
these reasons dust+gas simulations should be treated with caution.

# Numerics

Currently the only physics module that allows for the presence of dust
is HD, for a lack of a clear formulation of the coupling between dust
and magnetic fields. Each dust species has its own continuity and
momentum equation, but lacks an energy equation because there isn't a
pressure equivalent.  During each timestep the dragforce between dust
and gas is calculated as a local value and applied to the momentum
equations of both dust and gas. A stopping time is limiting
the next timestep in order to guarantee numerical stability.

Optionally, one can make this gas/dust coupling asymmetric or remove it
completely:

- `&dust_list: dust_method = none` simply remove the coupling
- `&dust_list: backreaction` (default is `.true.`) controls the
  application of drag forces to gas this can be useful for simple
  physics tests in the limit of very small dust to gas ratios


# Practical use

In order to add dust to a simulation, the user has to do the following:

## Connecting `hd` and `dust` modules

Your parfile should contain the following

```fortran
&hd_list
  hd_dust = .true.
/
```

## Minimal `&dust_list` example

Activate the number of dust species

```fortran
&dust_list
  dust_n_species = 1      ! number of dust fluids
/
```

## Minimal `mod_usr.t` setup

The following arrays need to be initialized in the user setup.
- `dust_size(1:dust_n_species)`
- `dust_density(1:dust_n_species)`

Those arrays define, as the names and shapes suggest, the different
grain sizes and densities, in code units, for each "dusty
fluid". For instance

```fortran
dust_size(1:dust_n_species)    = 1.0d0
dust_density(1:dust_n_species) = 1.0d0
```

Furthermore, note that using dust requires `ndir+1` additional
boundary conditions per dust species in each `typeboundary_xxxx` line 
of the parfile.

Finally, you need to initialize densities and momenta for each dust
fluid. This is a minimal example of this, implemented as part
of `usr_init_one_grid`

```fortran
integer :: n
do n=1, dust_n_species
  w(ixO^S, dust_rho(n))  = 1d0
  w(ixO^S, dust_mom(:,n) = 0d0
end do
```

# Internals of the module and further extension

`dust_list` contains a boolean `dust_small_to_zero` (default is
`.false.`).  If set `.true.` this activates calls to `set_dusttozero`,
which in turn will replace densities and momenta with exact zeros
where density has become lower than a threshold value `&dust_list:
dust_min_rho`.  
