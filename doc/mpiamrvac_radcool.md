

Radiative cooling for AMRVAC  



This webpage describes the radiative cooling module currently implemented in
AMRVAC.  
  



Structure  



The cooling routine is based on two subroutines: coolinit and radcool, which
interact with the rest of the AMRVAC code and a separate series of
subroutines, which contain the individual cooling routines and can be called
from radcool.  
  
The subroutine coolinit is called only once (at the start of the simulation).
It reads in a cooling curve from a separate file and transforms the data in
that file into a cooling table that can be used by the code.  
  
The subroutine radcool is called at each timestep from the AMRVAC subroutine
addsource_cooling it selects which scheme you wish to use to calculate the
cooling and calls the relevant subroutine.  
  



Coolinit  



The subroutine coolinit uses the paramete 'coolcurve' from the par file to
select a cooling curve (log(Lambda(T)) and log(T)). This curve is read into
temporary arrays from a separate module. Then, using the parameter 'ncool' the
sub routine interpolates in the cooling curve to obtain a table of length
'ncool' that can be used later by the code. Interpolation is done by a second
order Lagrangian polynomial routine. Doing this interpolation in second order
for a large number of points has the advantage that interpolation during the
simulation can be done lineairly without loss of accuracy.  
  
In order to allow scaling of the relevant variables, the user has to provide
two scaling parmeters 'eqpar(Tscale_)' and 'eqpar(Lscale_)'. These variables
serve two purposes at once:  
  
They translate the cooling curve into dimensionless parameters and change:  
  
T --&gt; P/rho  
n2 Lambda[cgs] --&gt; rho2 Lambda[dimensionless]  
  
This latter part is a matter of convenience to avoid having to use the
hydrogen mass and the  
Boltzmann constant each time the cooling subroutine is called.  
Note that coolingcurves are defined in cgs, so the scaling has to go from cgs
to dimensionless.  
  
The conversion is done by:  
  
T --&gt; eqpar(Tscale_) * T  
Lambda --&gt; eqpar(Lscale_) * Lambda  
  
As a result the the scaling parameters can be calculated as:  
  
[eqpar(Tscale_)] =  kb/mh* (vs-2)  
  
[eqpar(Lscale_)] =  1/m2h * rhos ts (vs-2)  
  
with mh and kb the hydrogen mass and the Boltzmann constant in cgs.  
vs is the velocity scaling in [cm/s], rhos the density scaling factor in
[g/cm3] and ts the times scaling factor in [s]  
  
I have excluded the metalicity / ionization from this, since these are local
phenomena, which vary from place to place. We may want to implement these
parameters at a later date, but it should be done locally. Not in the table.  
  
'Coolinit' also includes the opportunity to set multiple cooling tables. This
may become useful at a later date, once we have a way of keeping track of the
metalicity.  
  
A third  variable that needs to be set is [eqpar(Mue_)], which specifies the
mean molecular weight per free electron.  
This must be done with caution. Some cooling curves, such as the SPEX curve,
are already scaled for ionization an metallicity. It is generally safe to have
eqpar(Mue_) = 1.0.  
  



  
At this stage, it is possible to choose from several different cooling tables,
set in the variable 'cooltable'  
in amrvac.par:  
  
The Delgano &amp; McCray 1972 table.  
  
The MacDonald &amp; Bailey 1981 table as implemented in ZEUS-3D, augmented by
the Delgano &amp; McCray table for low temperatures  
  
3 different tables from Mellema &amp; Lundqvist 2002.  For zero metalicity,
solar metalicity and WC star composition. The last one should be used with
extreme caution, as it cools so strongly that cooling instabilities will tend
toward infinitly small size, limited only by the grid.  
  
A table generated with the SPEX code, as described in Schure et al. 2009. This
is already scaled for metallicity and ionization states, so  
eqpar(Mue_)==1  
  
2 tables generated with the Cloudy code for ism and solar metallicity. These
tables extend to low temperatures and take into account cooling by molecular
lines. As with the SPEX curve, they should have eqpar(Mue_)==1.  
  
With the exception of the SPEX table (which is rather more complicated), all
tables have been extended to 109 K, by assuming pure Bremm-Strahlung at high
temperatures (Lambda ~ sqrt(T)). This has been done to allow the different
numerical cooling methods to extend to higher temperatures.  
  



Radcool



Radcool as a subroutine consists only of a single 'select case' statement.  
It uses the input variable 'coolmethod' to select which subroutine to call for
the radiative cooling once the simulation is running. Admittedly, I could have
set this during the pre-compilation, but I wanted to retain the possibility of
using multiple cooling schemes in the same simulation (useful in the case of
strong shocks), or switching to a different scheme after a restart.  
  
The possible methods are:  
  
explicit 1: This is a simple, explicit cooling scheme that uses a call to
getdt to set an upper limit on the timestep in order to avoid numerical errors
in the radiative cooling. It is reliable, but can be very slow.  
  
explicit 2: A multistep expicit cooling scheme that divides the hydrodynamic
timestep into a number of small cooling steps, calculating an explicit cooling
value for each 'sub-step'. It is somewhat less reliable than 'explicit 1', but
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
the cooling function, which technically doesn't exist, as the coolingcurve is
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
larger temperatures the cooling is calculated explicitly. (I consider this
accetable since most cooling tables  
have 108K as a maximum value. At higher temperatures, Brehmstrahlung dominates
the cooling. Since this is a simple L~K1/2 relationship, the explicit value is
usually correct. In any case, areas with such a high temperature tend to have
low density, rendering the radiative cooling negligeable.  
  



Secondary routines  



There are several additional subroutines that have been added to fascilitate
the radiative cooling routine:  
  
addsource_cooling: This routine calls radcool (and if so specified
floortemperature) as part of the source term calculation.  
  
getdt_cooling: This subroutine checks the radiative cooling and internal
energy in each gridcell and ensures that the next simestep is limited to a
pre-set fraction of 'e/dedt'. It is only used for the explicit1 method.  
  
getvar_cooling: This subroutine can be called from getvar_output and provides
an extra output variable:  
The energy loss per unit time per volume due to the cooling.  
  
floortemperature: This subroutine can be called after the radiative cooling
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
  
  

  



In practice  



In order to use radiative cooling in a simulation, the user has to do the
following:  
  
1)In the usr file,  
put in the line:     INCLUDE: amrvacusr.cooling.t  
put a call to         addsource_cooling         in specialsource  
put a call to         getdt_cooling                  in getdt_special  
put a call to         coolinit                            in
initglobaldata_usr  
  
2) In the usrpar file put in the line:  
                            INCLUDE: 'amrvacusrpar.cooling.t'    
  
3) Edit the usrpar file file to make sure that the specialpar and
specialparname entries include:  
  
INTEGER,PARAMETER :: Tscale_ = neqpar + 1  
INTEGER,PARAMETER :: Lscale_ = Tscale + 1  
INTEGER,PARAMETER :: Mue_ = Lscale + 1, nspecialpar = 3  
  
CHARACTER*11,PARAMETER :: specialparname = 'Tsc Lsc Mue'  



Variables



  
This is a list of the new variables I have created in order to implement the
radiative cooling scheme for  
normal hydrodynamics.  
  
N.B. It ONLY contains those variables that have global implications.  
It does NOT contain lists of local variables  
  
tlow                        (Double precision) A lower temperature limit.  
                                Normally set in coolinit to be the lowest value on   
                                the cooling curve, but can be set by the user after   
                                the call to coolinit. This should be done with caution.   
                                Allowing temperatures below the lower limit of the cooling   
                                curve can lead to numerical problems.   
  
Lcool(1:10,000)       (Double precision) declared in amrvacusrpar.cooling.t  
                                Luminosity as function of temperature that is used   
                                in the actual cooling scheme.   
                                1000 is upper limit of points. Can be used with less.  
  
tcool(1:10,000)         (Double precision) declared in amrvacusrpar.cooling.t  
                                temperature for which cooling is given that is used   
                                in the actual cooling scheme.  
                                1000 is upper limit of points. Can be used with less.  
            
dLdtcool(1:10,000)  (Double precision) declared in amrvacusrpar.cooling.t  
                                derivative of cooling to temperature, to be used in   
                                Newton-Raphson cooling scheme, if I ever implement this.  
  
Yc(1:10,000)            (Double precision) declared in amrvacusrpar.cooling.t  
                                 dimensionless cooling time parameter used for exact integration method  
  
maxiter                  (integer,parameter) declared in
amrvacusrpar.cooling.t  
                               maximum number of iterations for implicit cooling scheme   
                               set to 100, but can be changed.  
            
e_error                 (Double precision, parameter) declared in
amrvacusrpar.cooling.t  
                              maximum relative error in implicit cooling scheme.   
                              Set to 1.0D-6, but can be changed.  
  
            
ncool                      (integer) declared in amrvacdef.t  
                               number of points in cooling curve that will be used.   
                               Default: 100, maximum 10000. Note that this is the number   
                               of the final curve, not the one that is read in.   
                               Can be set in 'methodlist'   
                               (When using the exact integration method,  it is advisable to set this number large.     
                               The efficiency of this method make this an acceptable choice and it increases accuracy.)  
            
cfrac                       (double precision) declared in amrvacdef.t  
                               fraction of E/L that sets upper time limit   
                               in explicit cooling fucntions. Is the cooling equivalent   
                               of 'courantpar'   
                               Default: 0.1   
                               Can be set in 'paramlist'  
  
Tfix                         (logical) declared in amrvacdef.t  
                                if true, a preset temperature limit (tlow) is enforced as minimum   
                                temperature throughout the grid. Not physical, but can help overcome   
                                negative pressure issues in radiative cooling instabilities. It is strongly recommended   
                                to use this.  
                                Default: .false.  
                                Can be set in 'methodlist'  
  
Tref                        (Double precision) declared in
amrvacusrpar.cooling.t  
                               Reference value for the temperature, used for exact integration  
                               cooling method. Usually the highest value in the cooling curve.  
  
  
coolcurve              (character*79) declared in amrvacdef.t  
                              Name of cooling curve to be used.   
                              Default: 'DM'   
                              Other options: 'MB', 'MLsolar1', 'MLcosmol', 'MLwc', 'SPEX', 'cloudy_solar', 'cloudy_ism', 'multi'   
                              Can be set in 'methodlist'  
            
coolmethod            (character*79) declared in amrvacdef.t  
                              Numerical scheme to be used in radiative cooling   
                              Default: 'explicit2'   
                              Other options: 'explicit1', 'semiimplicit', 'implicit', 'exact'   
                              Can be set in 'methodlist'  
             
cmulti                    (integer) declared in amrvacdef.t  
                              Number of cooling curves to be used, if coolcurve=='multi'   
                              This option is not yet implemented   
                              Default: 1, maximum: 3   
                              Can be set in 'methodlist'  
  
tcmulti(1:3, 1:10,000) (Double precision) decared in  amrvacusrpar.cooling.t  
                              replacement for tcool, if cmulti &gt; 1   
                              This option is not yet implemented  
  
Lcmulti(1:3, 1:10,000) (Double precision) decared in amrvacusrpar.cooling.t  
                            replacement for Lcool, if cmulti &gt; 1    
                            This option is not yet implemented  
  
  



Issues  



  
The high metallicity cooling curve from Mellema &amp; Lundqvist causes
considerable numerical problems. The cooling becomes so strong that shells and
clumps try to reach infinite density as their temperature  drops to zero. A
physically acceptable way to deal with this has to be found.  
This is a general problem with radiative cooling, but the high metallicity
cooling curve takes it to extremes.  

  



