title: Radiative cooling in AGILE

The webpage describes how users can add radiative cooling to their setup in AGILE

# Physics

The implemented radiative cooling is optically thin and describes locally the energy loss by radiation. 
The source is
\[  
Q = -n_i n_e \Lambda(T)  
\] 
where \( n_i \) and \( n_e \) are the ion- and electron number densities. 
\( \Lambda(T) \) is the cooling funtion or curve and represents the efficiency.
This cooling curve can be a tabulated set of points obtained by detailed calculations (see e.g. Colgan et. al. 2008). 
These tables can be interpolated to high temperature resolution.
MPI-AMRVAC also offers cooling curves described as piecewise power laws (see
e.g. Hildner 1974); AGILE does not implement those.



# Practical use

In order to use radiative cooling in a simulation, the user has to adjust the .par file as indicated in this section.


## Activating the module

Radiative cooling is available for `hd`, `ffhd` and `mhd`, and is activated
in the physics-dependent namelist, for mhd simulations with

    &mhd_list
       mhd_radiative_cooling=.true.
    /

and similarly with `hd_radiative_cooling` in `&hd_list` or
`ffhd_radiative_cooling` in `&ffhd_list`. This is a compile-time switch: it
has to be set in the case's `agile.par`, and takes effect when the case is
rebuilt (see [Compile-time switches](equations.html#eq_compile_time)).
`srhd` has no radiative cooling.

## Customizing the setup

Additionally, specific parameters can be set in the rc_list. In the following table, the available parameters are briefly described with their possible values. More details on the different cooling tables and methods are given in section **Numerics and implementation**.

name | type | default | description
---|---|---|---
coolcurve | string | 'JCcorona' | Name of cooling curve to be used. <br>The available tables are given below; any other name stops the run with "This coolingcurve is unknown". 
coolmethod | string | 'exact' | Numerical scheme to be used in radiative cooling. <br>Only 'exact' is implemented in AGILE. MPI-AMRVAC's 'explicit1', 'explicit2', 'semiimplicit' and 'implicit' are not ported; leave this at its default.
ncool | integer | 4000 | The number of point that will be used in interpolating the cooling curve. <br>Note that this is the number of elements in the final cooling table. When using the exact integration method, it is advisable to set this number large. The efficiency of this method make this an acceptable choice and it increases accuracy.
cfrac | double precision | 0.1 | Read but not used: it limited the time step of MPI-AMRVAC's explicit methods, which are not ported.
Tfix | logical | F |  If true, a temperature limit (tlow) is enforced as minimum temperature throughout the grid. Not physical, but can help overcome negative pressure issues in radiative cooling instabilities. It is strongly recommended to use this. <br>(It should be noted that other source terms, such as a [user source term](equations.html#eq_usr_source) `addsource_usr` in `mod_usr.fpp`, can still modify the internal energy and hence the temperature.)
tlow | double precision | lowest temperature of cooling curve | Used as lower temperature limit if Tfix=.true., should be set in dimensionless units
rad_cut | logical | F | If true, the cooling is suppressed below the height `rad_cut_hgt` in the last coordinate (`x(ndim)`), by a factor exp(-(x(ndim) - rad_cut_hgt)²/rad_cut_dey²). Applied on Cartesian grids only.
rad_cut_hgt | double precision | 0.5 | Height below which `rad_cut` suppresses the cooling, in code units.
rad_cut_dey | double precision | 0.15 | Width of the Gaussian suppression of `rad_cut`, in code units.

The cooling term is always added unsplit, together with the fluxes; MPI-AMRVAC's `rc_split` switch is not available.



## Examples of rc_list

An rc_list using an interpolated table:

    &rc_list
       coolcurve='Colgan_DM'
       ncool=4000
       Tfix=.true.
    /


# Numerics and implementation

The radiative cooling module is at `src/physics/mod_radiative_cooling.fpp`.
The module is primarily based on two subroutines: radiative_cooling_init and 
radiative_cooling_add_source. The subroutine radiative_cooling_init is 
called only once, during initialization of the physics module. 
It reads in parameters of rc_list, selects pre-coded cooling curve data and sets up the interpolated tables. 
The subroutine radiative_cooling_add_source is called at each timestep to add the cooling
source.


## Primary subroutines

*radiative_cooling_init*

It uses the parameter 'coolcurve' from the .par file to select the
corresponding pre-coded coolingcurve data (log(Lambda(T)) and log(T)). 
Then, it interpolates the cooling curve to obtain a table of length 
'ncool'. Interpolation is done by a second order Lagrangian polynomial routine. 
Doing this interpolation in second order
for a large number of points has the advantage that interpolation during the
simulation can be done linearly without loss of accuracy. 
The table containing the temporal evolution function, used in the exact method, is also constructed. 

In order to allow scaling of the relevant variables, the module has to know the unit_temperature, unit_numberdensity, unit_time, and unit_pressure. 
So the user has to provide three independent units in `mod_usr`.
Note that cooling curves are defined in cgs units, so the scaling has to go from cgs
to dimensionless.

*radiative_cooling_add_source*

It adds the radiative cooling source to a cell with `cool_exact`, the exact
integration method below, and then applies `floortemperature` if
`Tfix=.true.`. It is called from the physics module's local source terms,
on the GPU, once per cell and stage.


## The cooling curves

The following table gives information about the available cooling curves.

name | reference
---|---
'JCcorona' | Colgan & Feldman 2008, the default <br>Only till 1.e4K, beware for floor T treatment
'Colgan_DM' | Colgan & Feldman 2008 at high temperatures <br>At lower temperatures: Dalgarno & McCray 1972 with a pre-set ionization fraction of 1.e-3
'TRML_narrow' | A model cooling curve for the turbulent radiative mixing layer case (`tests/hd/TRML`)
'TRML_wide' | A wider variant of 'TRML_narrow', used by `tests/hd/TRML`

MPI-AMRVAC's other tables ('DM', 'MB', 'MLcosmol', 'MLsolar1', 'MLwc',
'cloudy_ism', 'cloudy_solar', 'SPEX', 'SPEX_DM', 'Dere_corona',
'Dere_corona_DM', 'Dere_photo', 'Dere_photo_DM', 'Colgan') and its piecewise
power laws ('Hildner', 'FM', 'Rosner', 'Klimchuk', 'SPEX_DM_rough',
'SPEX_DM_fine') are not ported.

Above the highest temperature of a table, the cooling rate is extended by
assuming pure Bremsstrahlung (\( \Lambda \sim \sqrt{T} \)), see
`calc_l_extended` below.

All coolingcurves are implemented in cgs units and scaled to code units
with the units of the physics module.

## The numerical method

AGILE implements one method, `cool_exact`, selected by `coolmethod='exact'`.
MPI-AMRVAC also offers 'explicit1', 'explicit2', 'semiimplicit' and
'implicit'; those are not ported.

*exact:* 

A new method to calculate radiative cooling, based on the exact
integration method developed by R.H.D. Townsend
([arXiv:0901.3146](http://arxiv.org/abs/0901.3146 "Abstract" )). This method
is both faster and more reliable than an implicit scheme. 
Note that for interpolated cooling tables, the exact integration method is only used within its limits. For
larger temperatures the cooling is calculated explicitly. This is
acceptable since most cooling tables have 1.e8K as a maximum value. At higher 
temperatures, Bremsstrahlung dominates the cooling. Since this is a simple
 \( \Lambda \sim \sqrt{T} \) relationship, the explicit value is
usually correct. 
In any case, areas with such a high temperature tend to have
low density, rendering the radiative cooling negligible.


## Secondary subroutines

There are several additional subroutines that have been added to fascilitate
the radiative cooling routine:  
  
*getvar_cooling_exact:*

Returns the energy loss per unit time per volume due to the cooling, for
output, computed with the exact method.

*floortemperature:* 

This subroutine is called if Tfix=.true., after the radiative cooling has
been added. Where the thermal pressure of the state the stage started from
(`wCT`) is below ρ𝓡T_low, it adds the missing thermal energy to the updated
state. The test is on that starting state, not on the cooled one.  
  
*findl:* 

A quick search program that finds the cooling rate at a given
temperature, by interpolation in the table.  
  
*findy:* 

Similar to findl, but for the dimensionless temporal evolution function (TEF) used in the
exact integration method.  
  
*findt:* 

Similar to findl, but to find the correct new temperature after the
exact integration method has been used. It is the inverse operation of findy.  

*calc_l_extended*

Calculates the cooling rate for temperatures larger than the pre-coded boundary. 
It uses Bremsstrahlung, Λ ∝ √T, matched to the last entry of the table.



