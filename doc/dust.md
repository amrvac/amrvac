# Dust in MPI-AMRVAC

# MPI-AMRVAC dust implementation

This is an overview of the option that allows users to model dust+gas
flows in mpi-amrvac.

The source code can be found in `src/mod_dust.t`

# Physics

Dust in mpi-amrvac is treated as a pressureless fluid that coexists
with the actual gas. Dust and gas are coupled by a dragforce. It is
calculated based on the phsyical characteristics of gas and dust
grains according to Kwok 1975, page 584 (between equations 8 and 9).

The sticking coefficient can be set according to dedicated parameters.

The interaction between dust and gas is considered purely kinetic,
without the influence of magnetic fields.

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
equations of both dust and gas. It is also taken into account when
calculating the size fo the next timestep in order to guarantee
numerical stability.

Optionally, one can make this gas/dust coupling asymetric or remove it
completly:

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

None of these parameters can be omitted (but more are available)

```fortran
&dust_list
  dust_n_species = 1      ! number of dust fluids
  gas_mu = 1d0            ! molecular weight
  dust_temperature = 1d0
/
```

## Minimal `mod_usr.t` setup

The following arrays need to be initialized in a `usr_set_parameters`
routine.
- `dust_size(1:dust_n_species)`
- `dust_density(1:dust_n_species)`

Those array define, as the names and shapes suggest, the different
grain sizes and densities, in code units, for each "dusty
fluid". For instance

```fortran
dust_size(1:dust_n_species)    = 1d0
dust_density(1:dust_n_species) = 1d0
```

Furthermore, note that using dust requires `ndir+2` additional
boundary conditions per dust species in each `typeboundary_xxxx` line 
of the parfile, or only `ndir+1` if you set 
`hd_list:hd_energy=.false.`.

Finally, you need to initialize densities and momenta for each dust
fluid.  This is a minimal working example of this, implemented as part
of `usr_init_one_grid`

```fortran
integer :: n
do n=1, dust_n_species
  w(ixO^S, dust_rho(n))  = 1d0
  w(ixO^S, dust_mom(:,n) = 0d0
end do
```

In addition to the above minimal requirements, please note that some conversion factors
can be defined to translate back to physical units.
(`w_convert_factor(X)` for `X` in `rho_, mom(1, :),
gas_e_, dust_rho(:), dust_mom(:)`)

Also note that `lenght_convert_factor` is used when `dust_temperature_type == ism` or
`stellar`.

That should be enough for your code to compile *and* run.


# Internals of the module and further extension

`dust_list` contains a boolean `dust_small_to_zero` (default is
`.false.`).  If set `.true.` this activates calls to `set_dusttozero`,
which in turn will replace densities and momenta with exact zeros
where density has become lower than a threshold value `&dust_list:
dust_min_rho`.  It's notably called at the end of `dust_add_source`.
This behaviour should be imitated by every other routine add source
terms to the dust components.
