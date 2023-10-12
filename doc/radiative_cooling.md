# Radiative cooling in MPI-AMRVAC 

This webpage describes the optically thin radiative cooling module in MPI-AMRVAC

# Practical use

In order to use radiative cooling in a simulation, the user has to do the
following:  
  
1) In the physical name list of par file, in mhd

    &mhd_list
       mhd_radiative_cooling=.true.
    /

   or in hd


    &hd_list
       hd_radiative_cooling=.true.
    /

2) In par file, specify parameters for the module in rc_list:  

    &rc_list
       coolcurve='JCcorona'
       coolmethod='exact'
       ncool=4000
       Tfix=.false.
    /
  

# Parameter list

Here is a decription for each parameters in rc_list

    coolcurve  (character) Name of cooling curve to be used.
               Default: 'JCcorona'
               Other options: 'MB', 'MLsolar1', 'MLcosmol', 'MLwc', 'SPEX', 
               'cloudy_solar', 'cloudy_ism', 'multi'

    coolmethod (character) Numerical scheme to be used in radiative cooling 
               Default: 'exact'
               Other options: 'explicit1', 'semiimplicit', 'implicit', 'exact'

    ncool      (integer) number of points in cooling curve that will be used.
               Default: 4000. Note that this is the number of elements in 
               the final cooling table. When using the exact integration 
               method, it is advisable to set this number large. The efficiency 
               of this method make this an acceptable choice and it increases accuracy.)
    
    cfrac      (double precision) fraction of E/L that sets upper time limit
               in explicit cooling fucntions.
               Default: 0.1
    
    Tfix       (logical) If true, a preset temperature limit (tlow) is enforced 
               as minimum temperature throughout the grid. Not physical, but 
               can help overcome negative pressure issues in radiative cooling 
               instabilities. It is strongly recommended to use this.
               Default: .false.
    
    tlow       (Double precision) Used as lower temperature limit if Tfix=.true.. 
               By default set to be the lowest value on the cooling curve.

# Structure details

The radiative cooling module is at `amrvac/src/physics/mod_radiative_cooling.t`.
The module is based on two subroutines: radiative_cooling_init and 
radiative_cooling_add_source. The subroutine radiative_cooling_init is 
called only once during initialization of (hd, mhd) 
physics. It reads in parameters of rc_list and select a pre-coded cooling curve data and 
interpolates the data into a high-temperature-resolution cooling table.
The subroutine radiative_cooling_add_source is called at each timestep to add cooling
source.

radiative_cooling_init

The subroutine radiative_cooling_init uses the parameter `coolcurve` from 
the par file to select a cooling curve data (log(Lambda(T)) and log(T)). 
Then, it interpolates in the cooling curve to obtain a table of length 
`ncool`, which is a parameter. Interpolation is done by a second
order Lagrangian polynomial routine. Doing this interpolation in second order
for a large number of points has the advantage that interpolation during the
simulation can be done linearly without loss of accuracy.  
  
In order to allow scaling of the relevant variables, the module has to know
correct unit_temperature, unit_numberdensity, unit_time, and unit_pressure. 
So user has to provide three in-dependent units in [mod_usr.t](@ref amrvacusr.md).
Note that coolingcurves are defined in cgs units, so the scaling has to go from cgs
to dimensionless.  
  
At this stage, it is possible to choose from several different cooling tables,
set in the variable `cooltable`  in amrvac.par:  
  
The Delgano & McCray 1972 table.  
  
The MacDonald & Bailey 1981 table as implemented in ZEUS-3D, augmented by
the Delgano & McCray table for low temperatures  
  
3 different tables from Mellema & Lundqvist 2002.  For zero metalicity,
solar metalicity and WC star composition. The last one should be used with
extreme caution, as it cools so strongly that cooling instabilities will tend
toward infinitly small size, limited only by the grid.  

A table generated with the SPEX code, as described in Schure et al. 2009. This
is already scaled for metallicity and ionization states.

2 tables generated with the Cloudy code for ism and solar metallicity. These
tables extend to low temperatures and take into account cooling by molecular
lines.

With the exception of the SPEX table (which is rather more complicated), all
tables have been extended to 1.e9 K, by assuming pure Bremm-Strahlung at high
temperatures (Lambda ~ sqrt(T)). This has been done to allow the different
numerical cooling methods to extend to higher temperatures.  

radiative_cooling_add_source

It uses the input variable `coolmethod` to select which subroutine to call for
the radiative cooling.

The possible methods are:  
  
explicit 1: This is a simple, explicit cooling scheme that uses a call to
getdt to set an upper limit on the timestep in order to avoid numerical errors
in the radiative cooling. It is reliable, but can be very slow.  
  
explicit 2: A multistep expicit cooling scheme that divides the hydrodynamic
timestep into a number of small cooling steps, calculating an explicit cooling
value for each sub-step. It is somewhat less reliable than 'explicit 1', but
generally faster and does not mess with the hydrodynamical timestep
conditions.  
  
semi-implicit: This routine calculates the explicit cooling value, uses it to
find a new internal energy value and calculates a second cooling value from
the new temperature. Finally, the average between the two cooling values is
subtracted from the original energy value. Not as accurate as either explicit
method, but fast.  
  
implicit: This routine uses half-step refinement to find an implicit value for
the radiative cooling. Reasonably fast and accurate, but like all implicit
cooling schemes vulnerable to the fact that the implicit cooling may have
multiple solutions.  
I used half-step refinement, rather than the more elegant Newton-Raphson
technique for several reasons. Firstly, Newton-Raphson needs the derivative of
the cooling function, which technically does not exist, as the coolingcurve is
just a collection of measurements. This can be avoided by using a numerical
rather than analytical derivative, but such a solution is far from ideal.
Secondly, The typical shape of a radiative cooling curve makes the use of
Newton-Raphson problematic, since the calculation may end up in an infinite
loop, switching back and forth between extremes. Finally, the existence of a
lower boundary to the temperature (the lowest T value in the cooling curve)
presents a problem. If the Newton Raphson loop ends up outside its predefined
boundaries, serious errors may occur. There are numerical solutions to this
problem, but they tend to be less than reliable and usually time-consuming.  
  
exact: A new method to calculate radiative cooling, based on the exact
integration method devleoped by R.H.D. Townsend
([arXiv:0901.3146](http://arxiv.org/abs/0901.3146 "Abstract" )). This method
is both faster and more reliable than an implicit scheme. The exact
integration method is only used within the limits of the cooling table. For
larger temperatures the cooling is calculated explicitly. I consider this
accetable since most cooling tables have 1.e8K as a maximum value. At higher 
temperatures, Brehmstrahlung dominates the cooling. Since this is a simple
 L~sqrt(K) relationship, the explicit value is
usually correct. In any case, areas with such a high temperature tend to have
low density, rendering the radiative cooling negligeable.  
  
Secondary routines  

There are several additional subroutines that have been added to fascilitate
the radiative cooling routine:  
  
getdt_cooling: This subroutine checks the radiative cooling and internal
energy in each grid cell and ensures that the next simestep is limited to a
pre-set fraction of `cfrac`. It is only used for the explicit1 method.  
  
getvar_cooling: This subroutine can be called from getvar_output and provides
an extra output variable:  
The energy loss per unit time per volume due to the cooling.  
  
floortemperature: This subroutine can be called if `Tfix=.true.`, 
after the radiative cooling
has been calculated and enforces the lowest temperature in the cooling curve
as the minimum temperature in each gridcell.  
  
findL: A quick search program that finds the cooling rate at a given
temperature  
  
finddLdt: Similar to findL, except that it searches in a table for the
derivative of the cooling function. This subroutine is currently not in use,
but can be used for a Newton-Raphson based implicit scheme.  
  
findY: Similar to findL, but for the dimensionless cooling time used in the
exact integration method.  
  
findT: Similar to findL, but to find the correct new temperature after th
eexact integration method has been used.  
