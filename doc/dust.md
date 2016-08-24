# Dust in MPI-AMRVAC

# MPI-AMRVAC dust implementation

This is an overview of the option that allows users to model dust+gas flows in
mpi-amrvac.

# Physics

Dust in mpi-amrvac is treated as a pressureless gas that coexists with the
actual gas. Dust and gas are linked through a dragforce that is calculated
based on the phsyical characteristics of the gas and the dust grains according
to Kwok 1975,page 584 (between eqn 8 and 9), but with several options for the
sticking coefficient. The interaction between dust and gas is considered
purely kinetic, without the influence of magnetic fields. It is possible to
combine multiple dust species in a single simulation. This approach dust not
take into account the interaction between different dust species, which will
simply ignore each other's presence. Neither does it take into account the
ability of grains of a single species to interpenetrate. A collision between
two dust flows of the same dust species will therefore result in a sharp
feature that is not physical as the dust is brought to a stop by the
collision. For these reasons dust+gas simulations should be treated with
caution.

# Numerics

Currently the only physics module that allows for the presence of dust is HD,
due to the lack of a clear formulation of the coupling between dust and
magnetic fields. Each dust species has its own continuity and momentum
equation, but lacks an energy equation because of the lack of a pressure
equivalent. During each timestep the dragforce between dust and gas is
calculated as a local value and applied to the momentum equations of both dust
and gas. It is also taken into account when calculating the size fo teh next
timestep in order to guarantee numerical stability.

# Practical use

In order to add dust to a simulation, the user has to do the following:

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
