# Relativistic test particle module

The test particle module in MPI-AMRVAC solves the **special relativistic equations of motion** for test particles (i.e. particles guided by electromagnetic fields, without kinetic feedback from the particles to the electromagnetic fields). There are three options to resolve particle dynamics, listed in the order of most demanding in time to least demanding. The first is to solve the full **special relativistic equation of motion** for charged particles guided by the Lorentz force:

\f[\frac{d^2x_i}{d\tau^2} = F_{ik}\frac{dx_k}{d\tau}\f],

in which the Einstein summation convention is used and with \f$x_i = (ct, \mathbf{r})\f$ the particles position in spacetime, \f$c\f$ the speed of light in vacuum and \f$\tau\f$ the proper time. \f$F_{ik}\f$ is the electromagnetic field tensor consisting of the components of \f$\mathbf{B}\f$ and \f$\mathbf{E}\f$. The equations can be split by choosing a certain frame of reference. The first three components are similar to the equation of motion \f$d\mathbf{p}/dt = q(\mathbf{E}+\mathbf{v}\times\mathbf{B})\f$, in three-space. Where \f$\mathbf{p}=m_0 \gamma \mathbf{v}\f$ is the relativistic momentum, with \f$m_0\f$ is the rest mass of the particle and \f$\gamma = 1/\sqrt{1-v^2/c^2}\f$ the Lorentz factor, \f$q\f$ the particles charge, \f$\mathbf{E}\f$ and \f$\mathbf{B}\f$ the electromagnetic fields guiding the particle and \f$\mathbf{v}\f$ the particles velocity. The fourth component is the rate of change of energy of the particle.

Option two is solving the **relativistic guiding centre equations of motion** for charged particles (see [Northrop - Adiabatic charged-particle motion](http://onlinelibrary.wiley.com/doi/10.1029/RG001i003p00283/abstract)), consisting of three equations describing the (change in) guiding center position \f$\mathbf{R}\f$, parallel relativistic momentum \f$p_{\|} = m_0\gamma v_{\|}\f$ and relativistic magnetic moment \f$\mu_r\f$ in three-space:

\f[\frac{d\mathbf{R}}{dt} = \frac{\left(\gamma v_{\|}\right)}{\gamma}\mathbf{\hat{b}}+\frac{\mathbf{\hat{b}}}{B\left(1-\frac{E_{\perp}^{2}}{B^2}\right)} \times \Biggl\{ -\left(1-\frac{E_{\perp}^{2}}{B^2}\right)c\mathbf{E} + \Biggr. \nonumber\f]

\f[\frac{cm_0\gamma}{q}\left(v_{\|}^{2}\left(\mathbf{\hat{b}}\cdot\nabla\right)\mathbf{\hat{b}}+v_{\|}\left(\mathbf{u_E}\cdot\nabla\right)\mathbf{\hat{b}} + v_{\|}\left(\mathbf{\hat{b}}\cdot\nabla\right)\mathbf{u_E} + \left(\mathbf{u_E}\cdot \nabla\right)\mathbf{u_E}\right) + \nonumber\f]

\f[\Biggl. \frac{\mu_r c}{\gamma q}\nabla\left[B\left(1-\frac{E_{\perp}^{2}}{B^2}\right)^{1/2}\right]  + \frac{v_{\|}E_{\|}}{c}\mathbf{u_E}  \Biggr\}\f],

\f[\frac{d \left(m_0 \gamma v_{\|}\right)}{dt} =  m_0\gamma\mathbf{u_E}\cdot \left(v_{\|}^{2}\left(\mathbf{\hat{b}}\cdot\nabla\right)\mathbf{\hat{b}}+v_{\|}\left(\mathbf{u_E}\cdot\nabla\right)\mathbf{\hat{b}}\right) +\nonumber \\ 
qE_{\|} -\frac{\mu_r}{\gamma}\mathbf{\hat{b}}\cdot\nabla\left[B\left(1-\frac{E^{2}_{\perp}}{B^2}\right)^{1/2}\right]\f],

\f[\frac{d \left(m_0 \gamma^{*2} v^{*2}_{\perp}/2B^*\right)}{dt} = \frac{d \mu_{r}^{*}}{dt} = 0\f].

with \f$\mathbf{\hat{b}}\f$ the unit vector in the direction of the magnetic field, \f$E\f$ the amplitude of the electric field vector and \f$v_{\|}\f$ the component of the particle velocity vector parallel to \f$\mathbf{\hat{b}}\f$. The drift velocity, perpendicular to \f$\mathbf{B}\f$ is written as \f$\mathbf{u_E} = \frac{c\mathbf{E}\times\mathbf{\mathbf{\hat{b}}}}{B}\f$ and \f$v^{*}_{\perp}\f$ is the perpendicular velocity the particle has, where we chose the frame of reference moving at \f$\mathbf{u_E}\f$. The magnetic field in that frame is given by \f$B^* = B(1-E^{2}_{\perp}/B^2)^{1/2}\f$ up to first order. The relativistic magnetic moment \f$\mu_{r}^{*}\f$ is  also evaluated in the frame of reference moving at \f$\mathbf{u_E}\f$. The Lorentz factor is not constant, but oscillates at the gyrofrequency, \f$\gamma = \gamma^{*}(1-E^{2}_{\perp}/B^2)^{-1/2}\f$.

The third option is to solve the **advection equation** for neutral particles.

All three modules are compatible with all the geometries available in MPI-AMRVAC. The magnetic fields are interpolated from (SR)MHD data and the electric fields and derivatives are calculated in the particle module.

The particle module can be called by defining it in the **definitions.h** file. Then a choice can be made between solving the full equation of motion, applying the guiding centre approximation or just solving the advection equation. For example, for solving the guiding centre equations of motion, add the following lines to **definitions.h**:

	#define PARTICLES
	#define PARTICLES_GCA

Or:

	#define PARTICLES
	#define PARTICLES_LORENTZ

to solve the full equation of motion.

The following files from the **src** directory have to be included in the **amrvacusr.t** file:

	INCLUDE:amrvacmodules/handle_particles.t
	INCLUDE:amrvacmodules/integrate_particles.t

Or directly from your working directory. These files evolve and integrate the particles. They call the file to initiate particles: amrvacmodules/init_particles.t to initiate the particles and the particle module in which all subroutines necessary for particle communication are written: modules/mod_particles.t.

# Parameters

In the par-file there are few parameters which are relevant for the particle module. Most parameters can be chosen in the separate particle files described below. However, in the par-file one can choose whether the particles are evolved alongside the MHD simulation, or if they are separately evolved in a static MHD snapshot. To apply this option the following convert type has to be chose in the **amrvac.par** file:

	convert_type  = 'particlesmpi'

When this convert_type is not chosen, the particles will be evolved simultaneously with the MHD or HD simulation. When no particle data file is provided for the chosen snapshot(s), the code will automatically initialise the particles from amrvacmodules/init_particles.t. The timesteps chosen for output in the par-file can be used in amrvacmodules/init_particles.t to let particle data be outputted at the same timesteps as the MHD or HD simulation, or for instance at smaller timescales, like in the next example.

# Particle initialisation

In: amrvacmodules/init_particles.t in the subroutine: `init_particle_integrator()` the maximum number of iterations, the maximum time and the timesteps (both in real, physical time) can be chosen. Typical settings are:

```{.f90}
	!=============================================================================
	subroutine init_particle_integrator()
		
	use mod_particles
	use constants
	include 'amrvacdef.f'
	!-----------------------------------------------------------------------------
	
	itmax_particles   = 1000000000
	tmax_particles    = tmax * (UNIT_LENGTH/UNIT_VELOCITY)
	ditsave_particles = 8
	dtsave_ensemble   = (dtsave(2)) * UNIT_LENGTH/UNIT_VELOCITY
	dtheta            = 2.0d0 * dpi / 60.0d0
		
	end subroutine init_particle_integrator
	!=============================================================================
```

In: amrvacmodules/init_particles.t in the subroutine: `init_particles()` the total number of particles **Npart** can be defined. One can also indicate the particles which have to be followed individually (like in the example below for particle number 174969), with a higher temporal resolution, outputted in a separate file:

```{.f90}
	!=============================================================================
	subroutine init_particles()
	! initialise the particles

	use constants
	use mod_particles
	use mod_gridvars
	use Knuth_random
	include 'amrvacdef.f'

	integer, parameter                   :: Npart=200000
	!-----------------------------------------------------------------------------

	follow(174969)   = .true.
```

And the mass (for instance **CONST_me** for an electron or positron and **CONST_mp** for a proton) and charge (adding a plus sign for a positively charged particle and a minus sign for a negatively charged particle) of a particle:

```{.f90}
	particle(nparticles)%self%q      = + CONST_e
	particle(nparticles)%self%m      =   CONST_mp
```

In the same subroutine the particles can be positioned according to the wishes of the user. The particles are initialised according to a Maxwellian velocity distribution, but this can also be edited.

Standard, for the guiding centre module, the particles mass, charge, parallel momentum, magnetic moment, Lorentz factor and x, y and z position are evolved in time. For the full equation of motion, also the x, y and z component of the velocity vector are evolved instead of the parallel momentum.

Additional variables can be defined in the same subroutine as payloads. For instance the gyroradius of the particle:

```{.f90}
	  particle(nparticles)%self%payload(1) = sqrt(2.0d0*particle(nparticles)%self%m*particle(nparticles)%self%u(2)*absB)/abs(particle(nparticles)%self%q*absB)*CONST_c
```
The total number of payloads has to be defined in: modules/mod_particles.t to avoid segmentation faults.

# Particle integration

The particle integration is done in amrvacmodules/integrate_particles.t. The main integration is done with a fourth order Runge Kutta integrator with an adaptive timestep found in modules/mod_odeint.t. A minimum timestep can be chosen by setting **hmin** in the `integrate_particles()` subroutine. 

In the subroutine `integrate_particles()` the integration is done. The timestep is determined in `set_particles_dt()`. Generally, an Euler timestep is taken to determine an initial timestep, limited by a cfl condition (which can be made more restrictive by changing **cfl** and **uparcfl**) for the particles not to jump cells. The Euler timestep is chosen based on the ratio of the step in position and the velocity and on the ratio of the velocity and the acceleration of the particle. After the Euler integration the particle is brought back to its previous position and velocity and the minimum timestep is chosen out of the Euler timestep, and the two ratios mentioned before.

In the subroutine `integrate_particles_gca()` this Euler timestep is used to estimate if a particle ends in an internal ghost cell. If a particle is nearby an internal ghost cell, the timestep is restricted by choosing an **int_factor** smaller than one. Then the Runge Kutta integrator is called to solve for the position of the particle, the parallel momentum, the magnetic moment and the Lorentz factor for the guiding centre module, and the velocity of the particle for the full equation of motion and the advection module. The electromagnetic field and its derivatives are interpolated in the subroutines `derivs_gca()`, `derivs_advect()` for the guiding centre module and the advection module respectively and subroutines `get_e()` and `get_b()` for the full equation of motion in which the general interpolator is used: modules/mod_gridvars.t. 

In `integrate_particles()` the payloads, as set in amrvacmodules/init_particles.t are also updated.

# Specific particle subroutines

The most important subroutines and functions in the particle module modules/mod_particles.t are highlighted here. They can be edited accordingly, but the particle communication is done in this file and these routines rely somewhat more on parallel computing algorithms than the particle initialisation and integration. So care should be taken before editing.

## General parameters

In the modules/mod_particles.t module a few general parameters are important to set to avoid segmentation faults. The maximum number of particles and maximum number of particles per processor should always be higher than or equal to the number of particles in your simulation. The number of payloads as set in amrvacmodules/init_particles.t should be defined here as well:

```{.f90}
integer,parameter                      :: nparticleshi = 500000, nparticles_per_cpu_hi = 500000
integer, parameter                     :: npayload=1
```

## Output and visualisation

All particle data is assembled in the **particles.dat** files. For visualisation there are several options, listed here.

The particle output is generated in several subroutines in modules/mod_particles.t. There are several options, namely to output a a file called **ensemble** containing all the particles data, to output a file for individual particles and to output the particles which are destroyed, called **destroy**, based on user defined criteria or because they left the physical domain. 

For individual particles, two output types are available. If **follow(ipart)** is called in amrvacmodules/init_particles.t for particle number **ipart**, as described above, one file is written with all the data for the particles indicated, called **followed**. A second file is written per particle in which high resolution data is written for just this particle, called **particle** with the particle number indicated in the filename. 

All these files are in **csv** format, which can be opened with paraview, visit, matlab, python, excell or another program of your choice.

## Boundary conditions

The particles generally follow the boundary conditions as set for the (M)HD simulation. However there are few exceptions possible. One can choose to destroy particles which leave
the domain or which reach a certain energy. One can also choose to inject particles at another position, or to let them reflect in a certain way. This is all done in the subroutines and functions described here.

## Particle communication

The most important and useful subroutines and functions in modules/mod_particles.t are listed here. 
The subroutine `find_particle_ipe()` can be called to find the processor and the grid block the particle is in, simply by calling it and giving the particles position vector. 

The function `particle_in_domain()` can be called to give a true or false statement to indicate whether the particle is in the physical domain or not, based on the particles position vector. The physical domain does not include the ghost cells.

The function `particle_in_igrid()` can be called to give a true of false statement to indicate whether the particle is in the grid or not. The grid does include the ghost cells.

The function `particle_in_zboundaries()` can be used in 2.5D simulations, when the particle dynamics are solved in 3D. Artificial z-boundaries for the particles can be chosen (set **zboundarymin** for the lower boundary and **zboundarymax** for the upper boundary) and this routine tells if a particle has crossed the boundaries in the invariant direction, after which the user can decide to destroy the particle or inject it at another position or with another velocity.

The function `particle_in_ghostcell()` can be called to give a true or false statement to indicate whether the particle is in the (internal) ghost cells or not. It uses the same interpolation as used for the Runge Kutta integration in modules/mod_gridvars.t 

The subroutine `apply_periodb()` applies periodic boundary conditions for a particle, in accordance with any applied periodic boundary conditions in the (M)HD simulation

The subroutine `apply_2_5D_periodb()` applies periodic boundary conditions at some user-set artifical z-boundaries in the function `particle_in_zboundaries()` and then thermalises the particle according to a Maxwellian velocity distribution and repositions it at the opposite boundary.

The subroutine `apply_periodb_thermalisation()` does the same as the `apply_periodb()` subroutine, but on top of periodically reinjecting the particle, it is also thermalised according to a Maxwellian velocity distribution

The subroutine `apply_hardwallb()` applies hard wall boundaries for a particle and reflects the particles velocity.

All these subroutines can be called in the `comm_particles()` subroutine based on the users wishes. Other subroutines can be written for particular boundary conditions a user wishes to apply. The subroutine `comm_particles()` checks if the particle should be destroyed, left or moved to a new position and grid block based on user defined criteria.

The subroutine `destroy_particles()` is called to destroy particles, either when the simulation has ended, when the particles leave the physical domain and are not reinjected or when some user defined criterion is met. 

# Contact details

For more specific questions, email bart.ripperda@kuleuven.be

