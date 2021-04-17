# Thermal conduction in MPI-AMRVAC

# MPI-AMRVAC thermal conduction implementation

This is an overview of the option that allows users to add thermal condution.

# Physics

Thermal conductivity is Spitzer conductivity for fully ionized plasma

# Numerics

Thermal conduction equation is solved seperately in an operator-split fashion from
the (M)HD equations. An explicit update and super time stepping time integrator is
used to achieve second order accuracy in time. (Slope limited) symmetric scheme is used
for discretization.

# Practical use

In order to add thermal conduction to a simulation, the user has to do the following:

  1. In the input file xxx.par, add hd_thermal_conduction=.true. in the hd_list for HD 
     problems or add mhd_thermal_conduction=.true. in the mhd_list for MHD ones.
  2. In the mod_usr.t, specify normalization units for length, number density, temperature (or velocity).
     in mod_usr.t, subroutine usr_init(), add 

           unit_length=<your length unit>
           unit_numberdensity=<your number density unit>
           unit_temperature=<your temperature unit>
           unit_velocity=<your velocity unit>

     before call (m)hd_activate()
  3. You can add the name list for thermal conduction, tc_list, in par file, where parameters tc_perpendicular, 
     tc_saturate, tc_slope_limiter can be modified. tc_perpendicular=.true. will add conduction flux
     perpendicular to magnetic field, which is not considered by default. tc_saturate=.true. to consider
     saturation effect of thermal conduction, which is the default choice. 
