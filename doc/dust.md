# Dust in MPI-AMRVAC

# MPI-AMRVAC dust implementation

This is an overview of the option that allows users to model dust+gas flows in
mpi-amrvac.

# Physics

Dust in mpi-amrvac is treated as a pressureless fluid that coexists with the
actual gas. Dust and gas are coupled by a dragforce. It is calculated
based on the phsyical characteristics of gas and dust grains according
to Kwok 1975,page 584 (between eqn 8 and 9).

The sticking coefficient can be set according to dedicated parameters.

The interaction between dust and gas is considered purely kinetic,
without the influence of magnetic fields.

It is possible to combine multiple dust species in a single
simulation, although there's no direct interaction (collisions,
fractionation, coagulation) between them.

Because particle collisions are not taken into account (e.g. dusty
fluids are pressureless), a shock happening in a dusty flow will
result in unphysical sharp features as dust is brought to a stop.
For these reasons dust+gas simulations should be treated with caution.

# Numerics

Currently the only physics module that allows for the presence of dust is HD,
for a lack of a clear formulation of the coupling between dust and
magnetic fields. Each dust species has its own continuity and momentum
equation, but lacks an energy equation because there isn't a pressure
equivalent.
During each timestep the dragforce between dust and gas is
calculated as a local value and applied to the momentum equations of both dust
and gas. It is also taken into account when calculating the size fo the next
timestep in order to guarantee numerical stability.

Optionnaly, one can make this gas/dust coupling asymetric or remove it completly:

- ``&dust_list: dust_method = none`` simply remove the coupling
- ``&dust_list: backreaction`` (default is `.true.`) controls the
  application of drag forces to gas this can be useful for simple
  physics tests in the limit of very small dust to gas ratios


# Practical use

In order to add dust to a simulation, the user has to do the following:

# THIS IS HEAVILY DEPRECATED, WILL BE REWRITTEN


  1. In definitions.h remove the line #undefine DUST
  2. recompile with $AMRVACDIR/setup.pl -p=hd -ndust=#, with # the number of dustspecies
  3. In the subroutine initglobaldata_usr define the size of each dust species (sdust(1:^NDS)) and the internal density of each dust species (rhodust(1:^NDS) as well as the scaling factors normvar(rho^DS) and normvar(v^Cd^DS). Note that both size and density of each species have to be scaled. The user may also need to set the values Lstar and/or Tdust
  4. In the subroutine initonegrid_usr set the initial values for the dust densities (w(ix^S,rho^DS)) and velocities (w(is^S,v^Cd^DS)). This may also need to be done in a special boundary, depending on the nature of the simulation.
  5. Adjust the amrvac.par file to determine how the dust is treated. This involves the following parameters, which have to be included in methodlist.
    1. dustmethod, which can be set to 'Kwok','sticking', 'linear',or 'none'. In case of 'Kwok', the dragforce will be calculated according to  Kwok 1975,page 584 (between eqn 8 and 9). In case of 'sticking', the same equations apply, but the sticking coefficient will be calculated according to Decin et al. 2006. In case of 'linear', the dragforce will simply scale with the velocity difference and in case of 'none', no drag force will exist at all.
    2. dustzero, if set to .true. will reduce the dust density to zero if it drops below the level of smallrhod in order to avoid numerical problems at extremely low densities
    3. dustspecies, which can be 'graphite' or 'silicate' affects the dust temperature. Note that it is up to the user to make sure that the dustspecies and internal density are compatible.
    4. dusttemp, is either 'constant', 'ism', or 'stellar'. Determines the dust temperature. If 'constant', it will be set to teh values of Tdust (in Kelvin). If 'ism', it will be calculated accoridng to Tielens (2005) eqn. 5.41 and 5.42. If 'stellar', it will be calculated according to Tielens (2005), eqn. 5.44 using a stellar luminosity of Lstar in solar luminosities.
    5. smallrhod is the cutoff below which the dust density will be set to zero, if dustzero == T
    6. In addition, the suer has to account for the additional variable. For example:
      1. primnames='rho v1 v2 p rhod1 rhod2 v1d1 v1d2  v2d1 v2d2 '
      2. wnames     = 'rho m1 m2 e rhod1 rhod2 m1d1 m1d2 m2d1 m2d2'
      3. typeB=   10*'special',
               10*'cont',
               'symm','symm', 'asymm','symm',2*'symm',2*'symm',2*'asymm',
               'symm','symm', 'asymm','symm',2*'symm',2*'symm',2*'asymm'
      4. flags(11)=1
flags(1)=1
