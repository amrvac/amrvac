# Test stressedIsland

\test stressedIsland

# Setup instructions

Setup this test case with

    $AMRVAC_DIR/setup.pl -d=23 -phi=0 -z=o -g=16,16 -p=mhd -eos=default -nf=1 -ndust=0 -u=nul -arch=default

# Description

In this test we investigate the particle dynamics during the evolution of an
ideal MHD instability. The test can be set up in 2.5D or 3D, but it is advisable
to start in 2.5D due to time restrictions. In 2.5D there are two magnetic
islands in the x-y plane which are perturbed by a velocity field. The islands
begin to tilt and form an environment favourable for reconnection and particle
acceleration. In 3D these islands correspond to flux ropes in the vertical
direction. However, in a 2.5D simulation we can also set a guide field
\f$B_{z0}\f$. This will have a (de)stabilising effect on the tilt instability,
depending on the magnitude of \f$B_{z0}\f$, through the tension in the flux
ropes. And hence on the acceleration of particles in and around the flux ropes.
The MHD setup and equilibrium is defined in amrvacusr.t as usual. The strenght
of the guiding field can be chosen in the amrvac.par file via setting `iprob'.

Based on the relativistic guiding center approximation, where the gyrating
motion of the particles is neglected, the particle dynamics are analysed in each
timestep of the MHD evolution. It is interesting to identity the regimes and the
regions where particles (electrons or protons) accelerate most and where the
fastest particles end up in the simulation box. The MHD simulation can also be
done separately, after which an interesting snapshot can be chosen and particles
can be evolved just in this single snapshot. The output is the ensemble of all
particles at all timesteps chosen, giving velocity, position and Lorentz factor
\f$\gamma\f$. From this an energy spectrum can be made. The second output file
gives the trajectory of the particles which are selected beforehand in
init\_particles.t, where it can also be chosen how to and where to initialise
the particles in the simulation box. By default, electrons are initialised
randomly, with a thermal velocity.

To get started, snapshots of the MHD simulation are provided in 2.5D setup with
a guide field \f$B_{z0}=0.1\f$. You can try to make plots of the magnetic field
magnitude in the \f$(x-y)\f$-plane, to see where the current sheet has formed.
In this region it is interesting to see how particles behave, so from this
snapshot a particle simulation can be started. You can integrate 100000
particle-orbits with the guiding centre approximation, which can be proven to be
accurate for the chosen MHD conditions. Notice how the particles are accelerated
into the current channels, and in the regions between the current channels. Try
to make histograms of the particle energy or equivalently, the Lorentz factor,
which can be computed from \f$u_1 = \gamma v\f$ or \f$u_3 = \gamma\f$ as \f$E =
(\gamma-1)*m_0 *c^2\f$. Remember that for low, non-relativistic energies this is
approximated by a Taylor series of which the first term is \f$\frac{1}{2} m
v^2\f$. For more specific information on how to analyse particle data, see the
descriptions of the StressedIsland tes, the coalescenceInstability test and the
particleSnapshot test.

The guiding center approximation is described in Northrop's book on the
adiabatic motion of charged particles, which is freely available:

http://babel.hathitrust.org/cgi/pt?id=mdp.39015048204716;view=1up;seq=7

The MHD setup is described in:

*Interacting tilt and kink instabilities in repelling current channels*, Keppens, R., Porth, O. \& Xia, C. 2014, ApJ, 795, 77 `doi:10.1088/0004-637X/795/1/77 `

Recent, related work on particle dynamics in a solar environment, based on the guiding center approximation can be found in:

*Particle acceleration and transport in reconnecting twisted loops in a
stratified atmosphere*, M. Gordovskyy, M., Browning, P.K., Kontar, E.P., \&
Bian, N.H. 2014, A\&A, 516, A72 `doi:10.1051/0004-6361/201321715`


