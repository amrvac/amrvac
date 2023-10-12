# Bondi-Hoyle-Lyttleton accretion onto a compact object

## The BHL problem

   This setup captures the hydrodynamical solution of the Bondi-Hoyle-Lyttleton (BHL) problem [5] : the deflection by a gravitational point-source of a flow which is, at infinity, planar, has a uniform mass density, pressure and speed relatively to the point-source (hereafter, the accretor). For numerical reasons, it is more convenient to model an accretor static in the mesh and a flow in motion but physically, the BHL problem refers to the equivalent configuration of accretors in fast motion in an ambient medium - such as the runaway star Zeta Ophiuci in the left panel of the Figure below or the newly formed neutron star associated to the Lighthouse Nebula (right panel), where the ambient medium in both case is the interstellar medium. 

![](figmovdir/runaway_star_lighthouse_nebula.png)

   Although the BHL problem appears as an axisymmetric extension of the simpler isotropic problem of Bondi spherical accretion (where the accretor is static with respect to the ambient medium), the seminal papers of Bondi, Hoyle and Lyttleton were published earlier, in the late-30’s / early-40’s. However, to alleviate the major difficulty introduced by the multi-dimensionality of the problem*, a dramatic assumption was initially made : the flow was assumed to be zero-temperature, ie ballistic, ie with an infinite Mach number, ie infinitely supersonic flow. Each fluid parcel could then be treated as a test-mass whose trajectory can be determined from firsthand principles. The main guideline obtained from this approach is that all the test-mass, whatever their impact parameter, converge towards a line in the wake of the accretor, the accretion line, where viscous heating must come in. The ballistic assumption, which relies on a zero-temperature flow, can no longer hold. If we assume that all the dissipation takes place along the accretion line, the normal component of the test-mass velocity is suddenly suppressed by the dissipation. The comparison between the remaining specific (ie per mass unit) kinetic energy and the gravitational potential yields two categories of streamlines : 

#### - those with an impact parameter larger than a critical one called the accretion radius, which are free and keep flowing away from the accretor, along the accretion line.
#### - those with an impact parameter smaller than the accretion radius, which are doomed to fall onto the accretor (aka “to be accreted”), still along the acretion line.

The key-quantity, the accretion radius, quantifies the amplitude of the graviational potential produced by the accretor with respect to the specific kinetic energy of the flow whose bulk motion is due to the relative speed with the accretor. It is thus a function of the mass of the accretor and of the square of the relative speed at infinity. This analysis, reproduced in sections 2.1 and 2.2 of [5], was initially carried out by Hoyle and Lyttleton in 1939 [9]. It provided the first estimation of the mass accretion rate by considering the mass flowing at infinity through a disk section of the cylinder of radius the accretion radius. The mass accretion rate obtained is proportional to the square of the mass of the accretor, the density of the flow at infinity and to the inverse of the cubed relative speed at infinity. 

   An elegant refinement, the accretion column model, was designed by Bondi and Hoyle in 1944 [2] to introduce a first glimpse of hydrodynamical (HD) effects and derive an associated lower limit on the mass accretion rate (see section 2.3 in [5]). Yet, in comparison with the Bondi spherical accretion model, the most notable absent in this problem remains the Mach number of the flow, which compares the bulk velocity to the sound speed (proportional to the temperature, assumed to be zero until now). An inspiring exercise is to show, by making the HD equations dimensionless, that the shape of the HD solutions depend only on the Mach number*2. The Bondi spherical accretion appears as the low relative speed limit of the BHL problem and its HD mass accretion rate can be analytically derived (see section 4.2 of El Mellah's PhD manuscript [7]). It drove Bondi into empirically modifying the BHL mass accretion rate formula to account for the Mach number of the flow at infinity - the so-called “interpolation formula” which however does not match the Bondi mass accretion rate for infinitely low Mach numbers.

   When numerical simulations came into the spotlight, they offered a brand new look at the BHL problem and in particular at the HD structure of the flow. It appeared that, well before reaching the accretion line, dissipation comes into play for finite Mach numbers and a bow shock is formed ahead of the accretor, with smaller opening angles of the accretion tail (in the wake of the accretor) for larger Mach numbers. Notice that in the BHL problem, the accretor is a point-source but accounting for its extension can change everything if it is of the order of the accretion radius, like in the case of accreting stars*3. In the following, we address only the case of compact objects, whose spatial extension is much smaller than the accretion radius for realistic relative speeds*4. The main science questions which can be addressed by numerical simulations are the following :

#### - what is the mass accretion rate onto the accretor and how does it evolve with the Mach number?
#### - what is the structure of the shock (opening angle, transverse profiles…) and how does it change with the Mach number?
#### - which observables can we deduce from the temperature distribution or the cooling/heating rates?
#### - what is the time-variability introduced by inhomogeneities in the flow speed or density upstream? Can it form a transient disc like structure around the accretor?

## Numerical considerations

   The current numerical setup is a two-dimensional one, mostly to afford a large spatial dynamics (ie a small inner boundary size compared to the outer boundary size). To capture the real 3D dynamics of the flow on a 2D mesh, in particular the 1/r^2 divergence term, a spherical mesh is compulsory. The invariance of the problem by rotation around the polar axis is then summoned to work in a 2D slice at constant azimuthal East-West angle and where the azimuthal speed is initially uniformly null and will remain so. It however precludes to investigate the impact of any inhomogeneity in the flow upstream if they have no reason to be axisymmetric. The need for a large dynamics originates from the need to capture both the shock and the tail, at the scale of the accretion radius, and an inner boundary larger than the compact object but small enough to not alter the structure of the shock and the mass accretion rate (see Figures 5.3 and 5.4 in El Mellah's PhD manuscript [7]). The gravitational source term in the equations must be accounted for when computing the time step. The contrast between the spatial scales (solved by the introduction of a stretched grid) translates here into an even more constraining contrast between the time scales : the number of iterations to cover the dynamical time scale of the BHL problem rises quickly with the ratio of the outer boundary radius by the inner boundary radius. Finally, we stress on the requirement for the energy equation to capture a meaningful HD solution. Indeed, albeit we make an adiabatic assumption ie we do not account for any heating or cooling term (see the section “Possible extensions” hereafter), the problem can in no case be considered as isentropic due to the presence of the shock. The polytropic assumption applies to an ideal gas if and only if its entropy is conserved ie if it does not exchange entropy with the surroundings (which it does not thanks to the adiabatic assumption) and if it undergoes infinitely slow mechanical reversible transformations (which it does not because of the shock which induce a jump in entropy). 

   In the amrvacusr.t, the gravitational source term has been introduced in the specialsource subroutine and accounted for to compute the time step in getdt_special. The scale of length is the accretion radius, the scale of speed is the speed at infinity and the scale of density is the density at infinity. The other scales can be derived by elementary (ie without supplementary geometrical factor) combinaisons of those three scales. The only shape parameter is the Mach number*5. The specialbound_usr sets the following boundary conditions, respectively at the inner and outer radii :

#### - continuity of the terms within the divergences in the mass, linear momentum and energy equations, since we aim at making a permanent state reachable. The gravitational source term can then be rewritten to be introduced within the divergence.
#### - planar flow upstream (bypassed by the bc_int subroutine, see below) and continuous condition downstream (safe as long as the flow is supersonically flowing out of the simulation space at the outer boundary).

Also, we use the internal boundary option to force the flow upstream to match the analytical ballistic solution derived in Bisnovatyi-Kogan et al 1979. Indeed, since the flow is supersonic in this region, it matches very closely the HD solution, which can be checked by not using the internal boundary option and using a planar flow upstream. Notice that in both case, since the outer boundary is at a finite distance, the bending of the flow from infinity to the outer boundary, albeit weak for an outer boundary at several accretion radii from the accretor, must be accounted for.

## Possible extensions

In order of increasing expected difficulty
#### - account for an optically thin cooling function
#### - handle the MHD case (see papers by LEE Aaron & BURLEIGH Kaylan’s team, [4] and [10]) 

## Footnotes 

*1 The BHL problem is 2D since the 3rd dimension can be dropped due to the axisymmetry of the problem around the axis of direction the relative speed and passing through the accretor.

*2 Such as the Reynolds number in HD turbulence, or the mass ratio in the Roche potential.

*3 See Ruffert94 on the impact of the spatial extension of the accretor on the mass accretion rate.

*4 See however Lora-Clavijo and Guzman 2012 and Gracia-Linares and Guzman 2015 for simulations of a relativistic flow onto a compact object, a configuration where the two scales are similar (like in the case of an accreting star) and where relativistic effects must be accounted for.

*5 Which is no longer true if one introduces a cooling function.

## Bibliography

1.['Accretion onto a rapidly moving gravitating center', Bisnovatyi-Kogan, Kazhdan, Klypin, Lutskii, Shakura, 1979, Soviet Astronomy 23](http://adsabs.harvard.edu/abs/1979SvA....23..201B)

2.['On the mechanism of accretion by stars', Bondi and Hoyle 1944, MNRAS 104](http://adsabs.harvard.edu/abs/1944MNRAS.104..273B)

3.['On spherically symmetrical accretion', Bondi 1952, MNRAS 112](http://adsabs.harvard.edu/abs/1952MNRAS.112..195B)

4.['Bondi-Hoyle accretion in a turbulent, magnetized medium', Burleigh, McKee, Cunningham, Lee, Klein, 2017, MNRAS 468](http://adsabs.harvard.edu/abs/2017MNRAS.468..717B)

5.['A review of Bondi–Hoyle–Lyttleton accretion', Edgar, 2004, New Astronomy Reviews 48](http://adsabs.harvard.edu/abs/2004NewAR..48..843E)

6.['Numerical simulations of axisymmetric hydrodynamical Bondi-Hoyle accretion onto a compact object', El Mellah & Casse, 2015, MNRAS 454](http://adsabs.harvard.edu/abs/2015MNRAS.454.2657E)

7.['Wind accretion onto compact objects', El Mellah, 2016, PhD manuscript](http://www.apc.univ-paris7.fr/~elmellah/PhD_manuscript.pdf)

8.['Accretion of Supersonic Winds onto Black Holes in 3D: Stability of the Shock Cone', Gracia-Linares, Guzmán, 2015, ApJ 812](http://adsabs.harvard.edu/abs/2015ApJ...812...23G)

9.['The effect of interstellar matter on climatic variation', Hoyle and Lyttleton 1939, Proceedings of the Cambridge Philosophical Society 35](http://adsabs.harvard.edu/abs/1939PCPS...35..405H)

10.['Bondi-Hoyle Accretion in an Isothermal Magnetized Plasma', Lee, Cunningham, McKee, Klein, 2014, ApJ 783](http://adsabs.harvard.edu/abs/2014ApJ...783...50L)

11.['Axisymmetric Bondi-Hoyle accretion on to a Schwarzschild black hole: shock cone vibrations', Lora-Clavijo and Guzmán, 2013, MNRAS 429](http://adsabs.harvard.edu/abs/2013MNRAS.429.3144L)

12.['Three-dimensional hydrodynamic Bondi-Hoyle accretion. 1: Code validation and stationary accretors', Ruffert, 1994, ApJ 427](http://adsabs.harvard.edu/abs/1994ApJ...427..342R)
