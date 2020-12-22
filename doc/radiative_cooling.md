# Radiative cooling in MPI-AMRVAC 

The webpage describes how users can add radiative cooling to their setup in MPI-AMRVAC

# Physics

The implemented radiative cooling is optically thin and describes locally the energy loss by radiation. 
The source is

\f$  Q = -n_i n_e \Lambda(T)  \f$, 

where \f$ n_i \f$ and \f$ n_e \f$ are the ion- and electronnumberdensities. 
\f$ \Lambda(T) \f$ is the cooling funtion or curve and represents the efficiency.
This cooling curve can be a tabulated set of points obtained by detailed calculations (see e.g. Colgan et. al. 2008). 
These tables can be interpolated to high temperature resolution.
Alternatively, some cooling curves are described as piecewise power laws (see e.g. Hildner 1974). 
They are explicit functions over the whole temperature domain.



# Practical use

In order to use radiative cooling in a simulation, the user has to adjust the .par file as indicated in this section.


## Activating the module

Depending whether the setup is in hd or mhd, the user has to activate the cooling in the physics-dependent namelist:

This can be done for mhd simulations with

    &mhd_list
       mhd_radiative_cooling=.true.
    /

or similarly for hd simulations with

    &hd_list
       hd_radiative_cooling=.true.
    /

## Customizing the setup

Additionally, specific parameters can be set in the rc_list. In the following table, the available parameters are briefly described with their possible values. More details on the different cooling tables and methods are given in section **Numerics and implementation**.

name | type | default | description
---|---|---|---
coolcurve | string | 'JCcorona' | Name of cooling curve to be used. <br>The available *interpolatable tables* are 'JCcorona', 'DM', 'MB', 'MLsolar1', 'MLcosmol', 'MLwc', 'SPEX', 'SPEX_DM', 'cloudy_solar', 'cloudy_ism', 'Dere_corona', 'Dere_corona_DM', 'Dere_photo', 'Dere_photo_DM', 'Colgan', 'Colgan_DM'. <br/>The available *piecewise power laws* are 'Hildner', 'FM', 'Rosner' and 'Klimchuk'. 
coolmethod | string | 'exact' | Numerical scheme to be used in radiative cooling. <br>The available methods are 'explicit1', 'explicit2', 'semiimplicit', 'implicit' and 'exact'.
ncool | integer | 4000 | The number of point that will be used in interpolating the cooling curve. <br>Note that this is the number of elements in the final cooling table. When using the exact integration method, it is advisable to set this number large. The efficiency of this method make this an acceptable choice and it increases accuracy.
cfrac | double precision | 0.1 | The fraction of E/L that sets upper time limit in explicit cooling functions.
Tfix | logical | F |  If true, a temperature limit (tlow) is enforced as minimum temperature throughout the grid. Not physical, but can help overcome negative pressure issues in radiative cooling instabilities. It is strongly recommended to use this. <br>(It should be noted that other source terms, such as usr_source set by the user in the mod_usr.t file, can still modify the internal energy and hence the temperature.)
tlow | double precision | lowest temperature of cooling curve | Used as lower temperature limit if Tfix=.true., should be set in dimensionless units
rc_split | logical | F | If true, the radiative cooling source term will be evaluated separately from the fluxes. This can ensure fixing of the temperature, if used in combination with Tfix and a splitting scheme ending on the splitted sources terms and no user-defined source terms are split. See [Discretization](@ref discretization.md) for more information on source splitting.



## Examples of rc_list

The rc_list to use an interpolated table:

    &rc_list
       coolcurve='JCcorona'
       coolmethod='exact'
       ncool=4000
       Tfix=.true.
    /
  
The rc_list to use a piecewise power law:

    &rc_list
       coolcurve='Rosner'
       coolmethod='exact'
       Tfix=.true.
    /


# Numerics and implementation

The radiative cooling module is at `mod_radiative_cooling.t`.
The module is primarily based on two subroutines: radiative_cooling_init and 
radiative_cooling_add_source. The subroutine radiative_cooling_init is 
called only once, during initialization of the hd or mhd physics. 
It reads in parameters of rc_list and selects pre-coded cooling curve data and sets up the interpolated tables or piecewise power law. 
The subroutine radiative_cooling_add_source is called at each timestep to add the cooling
source.


## Primary subroutines

*radiative_cooling_init*

It uses the parameter 'coolcurve' from 
the .par file to determine if the cooling curve is a piecewise power law or an interpolatable table. 
In case of the latter, the corresponding pre-coded coolingcurve data (log(Lambda(T)) and log(T)) is selected. 
Then, it interpolates the cooling curve to obtain a table of length 
'ncool'. Interpolation is done by a second order Lagrangian polynomial routine. 
Doing this interpolation in second order
for a large number of points has the advantage that interpolation during the
simulation can be done linearly without loss of accuracy. 
The table containing the temporal evolution function, used in the exact method, is also constructed. 
If a piecewise power law was set as cooling curve, the pre-coded data is used to set up the power law structure. 
The constants of integration of the temporal evolution function, used in the exact method, are also determined. 

In order to allow scaling of the relevant variables, the module has to know the unit_temperature, unit_numberdensity, unit_time, and unit_pressure. 
So the user has to provide three independent units in [mod_usr.t](@ref amrvacusr.md).
Note that cooling curves are defined in cgs units, so the scaling has to go from cgs
to dimensionless.

*radiative_cooling_add_source*

It uses the parameter 'coolmethod' to select which subroutine to call for determining
the radiative cooling source. The subroutines of these numerical methods have as name cool_'coolmethod'.


## The cooling curves

The following table gives information about the available cooling curves.

<i></i>  | name | reference
---|---|---
Interpolatable tables | 'JCcorona' | Colgan & Feldman 2008 <br>Only till 1.e4K, beware for floor T treatment
<i></i> | 'DM'       | Dalgarno & McCray 1972
<i></i> | 'MB'       | MacDonald & Bailey 1981 <br>As implemented in ZEUS-3D, with values from Dalgarno & McCray 1972 for low temperatures
<i></i> | 'MLcosmol'  | Mellema & Lundqvist 2002 for zero metallicity
<i></i>	| 'MLsolar1'  | Mellema & Lundqvist 2002 for solar metallicity
<i></i>	| 'MLwc'  | Mellema & Lundqvist 2002 for WC-star metallicity <br>Use with extreme caution. It cools so strongly that cooling instabilities will tend toward infinitly small size, limited only by the grid
<i></i> | 'cloudy_ism' | Generated with the Cloudy code (Ferland et al. 1998) for ism metallicity. <br>These tables extend to low temperatures and take into account cooling by molecular lines
<i></i> | 'cloudy_solar' | Generated with the Cloudy code for solar metallicity. <br>These tables extend to low temperatures and take into account cooling by molecular lines
<i></i> | 'SPEX' | Generated with the SPEX code, as described in Schure et al. 2009. <br>This is already scaled for solar metallicity and ionization states
<i></i> | 'SPEX_DM' | Above 1.e4K: generated with the SPEX code <br>At lower temperatures: Dalgarno & McCray 1972 with a pre-set ionization fraction of 1.e-3 
<i></i> | 'Dere_corona'  | Dere et al. 2009, generated with CHIANTI for solar coronal abundances
<i></i> | 'Dere_corona_DM'  | Dere et al. 2009, generated with CHIANTI for solar coronal abundances at high temperatures. <br> At lower temperatures: Dalgarno & McCray 1972 with a pre-set ionization fraction of 1.e-3	
<i></i> | 'Dere_photo'   | Dere et al. 2009, generated with CHIANTI for solar photospheric abundances
<i></i> | 'Dere_photo_DM'   | Dere et al. 2009, generated with CHIANTI for solar photospheric abundances at high temperatures. <br> At lower temperatures: Dalgarno & McCray 1972 with a pre-set ionization fraction of 1.e-3
<i></i> | 'Colgan'   | The original table of Colgan & Feldman 2008 
<i></i> | 'Colgan_DM'   | The original table by Colgan & Feldman 2008 <br> At lower temperatures: Dalgarno & McCray 1972 with a pre-set ionization fraction of 1.e-3
Piecewise power laws  | 'Hildner'  | Hildner 1974
<i></i> | 'FM'       | Hildner 1974, but modified with constraints on the radiative timescales as in Forbes & Malherbe 1991
<i></i> | 'Rosner'       | Rosner 1978, extended by Priest 1982 and implemented as given in the PhD thesis of Van Der Linden 1991
<i></i> | 'Klimchuk'   | Klimchuk et al. 2008

With the exception of the SPEX table (which is rather more complicated), 
all interpolatable tables have been extended to 1.e9 K, 
by assuming pure Bremsstrahlung at high temperatures (\f$ \Lambda \sim \sqrt{T} \f$). 
This has been done to allow the different numerical cooling methods to extend to higher temperatures.

All coolingcurves are in cgs. The piecewise power laws have a pre-coded temperature range of 1.e3 to 1.e10 K, for practical reasons.

## The numerical methods

The following methods are implemented in the subroutines with corresponding names cool_'coolmethod'.

*explicit1:*

This is a simple, explicit cooling scheme that uses a call to 
getdt to set an upper limit on the timestep in order to avoid numerical errors
in the radiative cooling. It is reliable, but can be very slow.  
  
*explicit2:* 

A multistep expicit cooling scheme that divides the hydrodynamic
timestep into a number of small cooling steps, calculating an explicit cooling
value for each sub-step. It is somewhat less reliable than 'explicit1', but
generally faster and does not mess with the hydrodynamical timestep
conditions.  
  
*semiimplicit:* 

This routine calculates the explicit cooling value, uses it to
find a new internal energy value and calculates a second cooling value from
the new temperature. Finally, the average between the two cooling values is
subtracted from the original energy value. Not as accurate as either explicit
method, but fast.  
  
*implicit:* 

This routine uses half-step refinement to find an implicit value for
the radiative cooling. Reasonably fast and accurate, but like all implicit
cooling schemes vulnerable to the fact that the implicit cooling may have
multiple solutions.  
Half-step refinement was used, rather than the more elegant Newton-Raphson
technique for several reasons. Firstly, Newton-Raphson needs the derivative of
the cooling function, which technically does not exist, as the coolingcurve is
just a collection of measurements or a non-continuous power law. This can be avoided by using a numerical
rather than analytical derivative, but such a solution is far from ideal.
Secondly, The typical shape of a radiative cooling curve makes the use of
Newton-Raphson problematic, since the calculation may end up in an infinite
loop, switching back and forth between extremes. Finally, the existence of a
lower boundary to the temperature (the lowest T value in the cooling curve)
presents a problem. If the Newton Raphson loop ends up outside its predefined
boundaries, serious errors may occur. There are numerical solutions to this
problem, but they tend to be less than reliable and usually time-consuming.  
  
*exact:* 

A new method to calculate radiative cooling, based on the exact
integration method developed by R.H.D. Townsend
([arXiv:0901.3146](http://arxiv.org/abs/0901.3146 "Abstract" )). This method
is both faster and more reliable than an implicit scheme. 
Note that for interpolated cooling tables, the exact integration method is only used within its limits. For
larger temperatures the cooling is calculated explicitly. This is
acceptable since most cooling tables have 1.e8K as a maximum value. At higher 
temperatures, Bremsstrahlung dominates the cooling. Since this is a simple
 \f$ \Lambda \sim \sqrt{T} \f$ relationship, the explicit value is
usually correct. For the piecewise power laws, the cooling is also calculated explicitly but with its actual value instead of Bremsstrahlung. 
In any case, areas with such a high temperature tend to have
low density, rendering the radiative cooling negligible.


## Secondary subroutines

There are several additional subroutines that have been added to fascilitate
the radiative cooling routine:  
  
*getdt_cooling:* 

This subroutine checks the radiative cooling and internal
energy in each grid cell and ensures that the next simestep is limited to a
pre-set fraction of `cfrac`. It is only used for the explicit1 method.  
  
*getvar_cooling:* 

This subroutine can be called from the user-defined, pointed usr_aux_output and provides
an extra output variable:  
The energy loss per unit time per volume due to the cooling.  

*getvar_cooling_exact:*

Similar to getvar_cooling, but uses the exact method.

*floortemperature:* 

This subroutine is called if Tfix=.true.. After the radiative cooling
has been calculated and if the new temperature is lower than tlow, it enforces tlow
as the minimum temperature in each gridcell.  
  
*findl:* 

A quick search program that finds the cooling rate at a given
temperature. It uses interpolation for the interpolatable tables and explicit formulae (Townsend 2009) for the piecewise power laws.  
  
*finddldt:* 

Similar to findl, except that it searches in a table for the
derivative of the cooling function. This subroutine is currently not in use,
but can be used for a Newton-Raphson based implicit scheme. It does not work for the piecewise power laws. 
  
*findy:* 

Similar to findl, but for the dimensionless temporal evolution function (TEF) used in the
exact integration method.  
  
*findt:* 

Similar to findl, but to find the correct new temperature after the
exact integration method has been used. It is the inverse operation of findy.  

*create_y_ppl*

Creates the constants of integration in the temporal evolution function (TEF) for piecewise power laws, 
according to eq. A6 of Townsend 2009. 

*calc_l_extended*

Calculates the cooling rate for temperatures larger than the pre-coded boundary. 
It uses Bremsstrahlung for interpolated tables and the explicit formula for piecewise power laws. 
It is used in the subroutines of the numerical methods.



