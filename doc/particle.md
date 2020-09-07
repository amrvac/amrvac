# Test particle module for HD/MHD

The test particle module in MPI-AMRVAC provides the possibility to track test particles, i.e. particles evolved according to the HD/MHD fields in a simulation without any feedback from the particles to these fields. Test particles can be added on top of an (M)HD state and evolved concurrently to the fluid. Alternatively, test particles can be evolved in a static snapshot, i.e. without evolving the underlying fluid quantities. The test particle module can also be used to sample the fluid quantities at specific locations (which may differ from the computational grid points) and output the results separately from the usual <tt>.dat</tt> and <tt>.vtu</tt> files.

The particle module can work in four different modes, depending on the user's choice:
* **Advection** (compatible with HD/MHD): In this mode particles represent Lagrangian blobs of fluid that are tracked through the physical domain. The fluid properties are interpolated at the particle location, which is advected by using the local fluid speed.
* **Charged particles: Full Lorentz dynamics** (compatible with MHD): This mode is used to solve the dynamics of charged particles, whose motion is determined by the electromagnetic fields of the underlying MHD state. The equations of motion for particles are formulated according to the full Lorentz dynamics. The user can switch between the Newtonian and the special-relativistic formulation of the equations of motion. While for the Newtonian case a standard Boris integrator is employed, several numerical integrators are available for treating the relativistic case.
* **Charged particles: Guiding centre approximation** (compatible with MHD): In this case, the dynamics of charged particles is treated according to the guiding centre approximation (GCA) formalism. The user can again switch between the Newtonian and the relativistic equations of motion.
* **Sampling** (compatible with HD/MHD): Particles represent fixed points in space where the fluid quantities are interpolated. The location and output cadence at these points can be different from those chosen for the computational grid points. This mode can be used to monitor fluid quantities at specific locations during a run.

In the next sections, the details of each mode are illustrated.


# Advection (HD/MHD)
This mode is used to track fluid parcels moving through the simulation domain with the fluid speed. At each step, a simple equation of motion is solved for each particle,
\f[ \frac{d\mathbf{x}}{dt} = \mathbf{v}, \f]
where \f$\mathbf{v}\f$ is the local fluid velocity, linearly interpolated at the particle position \f$\mathbf{x}\f$. The equation of motion above is solved by means of a fourth-order Runge-Kutta integrator with adaptive time-stepping.

As the Lagrangian motion of the fluid particles is computed, it is possible to track various fluid quantities at the particle locations, e.g. the fluid density, pressure, etc. See the "Payloads" section below for further details.


# Charged particles (MHD): Full Lorentz dynamics
A charged particle of charge \f$q\f$ and mass \f$m\f$ in electromagnetic fields evolves in time according to Newton's equations of motion
\f[ \frac{d\mathbf{x}}{dt} = \mathbf{v}, \f]
\f[ \frac{d\mathbf{v}}{dt} = \frac{q}{m}\left(\mathbf{E}+\mathbf{v}\times\mathbf{B}\right), \f]
where \f$\mathbf{x}\f$ is the particle position, \f$\mathbf{v}\f$ is the velocity, and the electric and magnetic fields \f$\mathbf{E}\f$ and \f$\mathbf{B}\f$ provide the accelerating Lorentz force for the particle. The typical particle trajectory is a superposition of a parallel motion along magnetic field lines, a gyro-motion around magnetic field lines, and various "drift" mechanisms that allow the particle to cross magnetic field lines. As the magnetic field does no work on charged particles, a change in the particle kinetic energy is only associated to the presence of electric fields.

The equations above are valid for a particle travelling at speeds much smaller than the speed of light in vacuum \f$c\f$. When \f$v\rightarrow c\f$, it is necessary to adopt the **relativistic equations of motion**
\f[ \frac{d\mathbf{x}}{dt} = \frac{\mathbf{u}}{\gamma}, \f]
\f[ \frac{d\mathbf{u}}{dt} = \frac{q}{m}\left(\mathbf{E}+\frac{\mathbf{u}}{\gamma}\times\mathbf{B}\right), \f]
where \f$\mathbf{u}=\gamma\mathbf{v}\f$ is the normalised particle momentum, with \f$\gamma=\sqrt{1+u^2/c^2}=1/\sqrt{1-v^2/c^2}\f$ the Lorentz factor.

For the relativistic case, the solution of the equations for the Lorentz dynamics can be carried out numerically in several fashions. Below we list the options available in MPI-AMRVAC:
* **Boris algorithm**: the equations of motion are solved with the Boris scheme (Boris 1970), popular for its simplicity and limited computational cost.
* **Vay algorithm**: this more recent scheme (Vay 2008) is more suitable for particles travelling at highly relativistic speeds.
* **Higuera-Cary algorithm**: this scheme exhibits smaller numerical errors in the particle gyro-phase with respect to the Boris algorithm (Higuera and Cary 2010). It provides advantages for ultra-relativistic particle motion.
* **Lapenta-Markidis algorithm**: this is a fully-implicit iterative scheme, and is therefore more expensive than the three previous schemes, which carry out the solution explicitly. The higher computational cost is compensated by the fact that this scheme introduces no numerical errors in the particle energy (Lapenta and Markidis 2011). In all cases (as well as for the Newtonian case) the particle time step is obtained by ensuring that the particle gyro-period is resolved by 60 points.


# Charged particles (MHD): Guiding centre approximation
In this mode, a "reduced" set of equations is employed to resolve the dynamics of charged particles, according to the so-called "guiding centre approximation" (GCA) paradigm. This neglects the particle gyro-motion, and may present advantages in certain physical scenarios. The equations of motion in this case are
\f[ \begin{aligned}
    \frac{d\mathbf{R}}{dt} = & \mathbf{v}_{\|} + \mathbf{v}_E + \frac{\hat{\mathbf{b}}}{B}\times\frac{\mu}{q}\nabla B \\
                             & + \frac{m\hat{\mathbf{b}}}{qB} \times\left\{ v_{\|}^2\left(\hat{\mathbf{b}}\cdot\nabla\right)\hat{\mathbf{b}} + v_{\|}\left(\mathbf{v}_E\cdot\nabla\right)\hat{\mathbf{b}} + v_{\|}\left(\hat{\mathbf{b}}\cdot\nabla\right)\mathbf{v}_E + \left(\mathbf{v}_E\cdot\nabla\right)\mathbf{v}_E \right\},
\end{aligned} \f]
\f[ \frac{dv_{\|}}{dt} = \frac{q}{m}E_{\|} - \frac{\mu}{m}\hat{\mathbf{b}}\cdot\nabla B + \mathbf{v}_E\cdot\left(v_{\|}\left(\hat{\mathbf{b}}\cdot\nabla\right)\hat{\mathbf{b}}+\left(\mathbf{v}_E\cdot\nabla\right)\hat{\mathbf{b}}\right), \f]
\f[ \frac{d\mu}{dt} = 0, \f]
where \f$\mathbf{R}\f$ is the guiding centre position and \f$\mathbf{v}_{\|}=(\mathbf{v}\cdot\hat{\mathbf{b}})\hat{\mathbf{b}}\f$ is the particle velocity in the direction parallel to the magnetic field. The unit vector \f$\hat{\mathbf{b}}=\mathbf{B}/B\f$ also defines the drift velocity \f$\mathbf{v}_E=\mathbf{E}\times\hat{\mathbf{b}}/B\f$ and the parallel electric field \f$E_{\|}=\mathbf{E}\cdot\hat{\mathbf{b}}\f$. The conservation of the adiabatic invariant \f$\mu=mv_\perp^2/(2B)\f$, with \f$v_\perp\f$ the particle velocity perpendicular to \f$\mathbf{B}\f$, is an assumption of this paradigm and may not be valid in general.

For particles travelling at velocities close to \f$c\f$, the special-relativistic formulation of the equations above reads
\f[ \begin{aligned}
    \frac{d\mathbf{R}}{dt} = & \frac{\mathbf{u}_{\|}}{\gamma} + \mathbf{v}_E 
    + \frac{\kappa^2\hat{\mathbf{b}}}{B}\times\left\{\frac{\mu_r}{q\gamma}\nabla (B/\kappa) + \frac{u_{\|}}{\gamma}E_{\|}\mathbf{v}_E\right\} \\
    & + \frac{m\kappa^2\hat{\mathbf{b}}}{qB}\times\left\{\frac{u_{\|}^2}{\gamma}\left(\hat{\mathbf{b}}\cdot\nabla\right)\hat{\mathbf{b}} + u_{\|}\left(\mathbf{v}_E\cdot\nabla\right)\hat{\mathbf{b}} + u_{\|}\left(\hat{\mathbf{b}}\cdot\nabla\right)\mathbf{v}_E + \gamma\left(\mathbf{v}_E\cdot\nabla\right)\mathbf{v}_E\right\},
\end{aligned} \f]
\f[ \frac{du_{\|}}{dt} = \frac{q}{m}E_{\|} - \frac{\mu_r}{m\gamma}\hat{\mathbf{b}}\cdot\nabla \left( B/\kappa\right) + \mathbf{v}_E\cdot\left(u_{\|}\left(\hat{\mathbf{b}}\cdot\nabla\right)\hat{\mathbf{b}}+\gamma\left(\mathbf{v}_E\cdot\nabla\right)\hat{\mathbf{b}}\right), \f]
\f[ \frac{d\mu_r}{dt} = 0, \f]
where \f$\mathbf{u}_{\|}=\mathbf{v}_{\|}\gamma\f$, and \f$\kappa=1/\sqrt{1-v_E^2/c^2}\f$ is the Lorentz factor corresponding to the frame moving at a speed equal to \f$\mathbf{v}_E\f$. These equations are obtained by averaging over the particle gyro-motion, and therefore assume the conservation of the relativistic magnetic moment \f$\mu_r=m\gamma^2 v_\perp^2/(2B)\f$. Because of the averaging, it can be shown that the expansion of the guiding centre velocity leads to the definition of the Lorentz factor \f$\gamma = \kappa\sqrt{1+(u_{\|}^2+2\mu_r B/m)/c^2} = 1/\sqrt{1-(v_{\|}^2+v_E^2+2\mu_r B/m)/c^2}\f$ up to first order. This implies that the dominant drift mechanism in the guiding centre motion is represented by \f$\mathbf{v}_E\f$. For both the Newtonian and relativistic cases, the assumption of slowly time-varying electromagnetic fields has been made, in order to simplify the equations of motion. This is justified by the reasonable assumption that the particle dynamics takes place on much faster time scales than that of the typical MHD evolution.

The solution of the equations above is carried out in MPI-AMRVAC with a fourth-order Runge-Kutta numerical integrator with adaptive time stepping.


# Scattered sampling (HD/MHD)
With this mode, the user is allowed to sample the fluid quantities at an arbitrary number of points in space, which may or may not be spatially distinguished from the numerical grid points. All quantities are retrieved at such locations via linear interpolation. By default, the quantities used in a standard HD/MHD run are computed at the sampling points (density, momentum, pressure/energy/temperature, and magnetic field components in the MHD case). The user is also free to define their own list of additional quantities to be sampled (e.g. the current).


# Usage and input file parameters
All particle-related computations can be activated by setting <tt>hd_particles=.true.</tt> in the \p hd_list (if running in HD) or <tt>mhd_particles=.true.</tt> in the \p mhd_list (if running in MHD) of the <tt>.par</tt> file. If the particle calculations are switched on, additional parameters can be specified in the \p particles_list. Below is a description of these parameters and their role in the particle modules.

## Initialisation
Below, we describe the essential steps needed to correctly set up a particle simulation.

- **Choice of mode**: the user can choose to run the particle module in the advection, Lorentz, GCA, or sampling modes. This is determined by the parameter \p physics_type_particles, which can be set to \p 'advect', \p 'Lorentz', \p 'GCA', or \p 'sample'.

- **Number of particles**: the most important parameter is the number of particles whose trajectory should be integrated during the calculations. The user can choose such a number via the parameter \p num_particles (an integer).

- **Particle position and velocity**: for all choices of \p physics_type_particles, the particles must be initialised at the desired locations in space. For the \p 'Lorentz' and \p 'GCA' modes, a velocity distribution should also be specified.
  
  The user can specify their own particle initialisation setup in the \p mod_usr.t file. Then, the user must associate the pointer \p usr_create_particles with their particle initialisation subroutine (named e.g. \p generate_particles). The format for such a routine *must* be:

  \code{.f90}
  subroutine generate_particles(n_particles, x, v, q, m, follow)
    integer, intent(in)           :: n_particles
    double precision, intent(out) :: x(3, n_particles)
    double precision, intent(out) :: v(3, n_particles)
    double precision, intent(out) :: q(n_particles)
    double precision, intent(out) :: m(n_particles)
    logical, intent(out)          :: follow(n_particles)

  end subroutine generate_particles
  \endcode

  Here, \p x is the particle position, \p v is the particle velocity, \p q is the particle charge, and \p m is the particle mass. The \p follow variable tells the routine whether a certain particle should be tracked individually during the simulation (see section "Output and visualisation" below). While the position will be used in all modes, information on the velocity, charge, and mass will only be employed if \p physics_type_particles is \p 'Lorentz' or \p 'GCA'. If the user does not provide their own particle initialisation subroutine (and/or if the latter is not associated with the \p usr_create_particles pointer), the code will by default initialise all particles at randomly generated locations inside the domain, with zero velocity, mass, and charge, and with <tt>follow=.false.</tt>.

- **Payloads**: particle simulations are especially flexible in terms of the quantities that can be dynamically stored in the particle output files. On top of tracking positions and velocities, an arbitrary number of payloads can be assigned to each particle in order to monitor additional physical aspects. As an example, in the \p 'advect' mode each particle can be assigned to track the local fluid density, which will be then stored in a payload variable and added to the output. The number of payloads is chosen by the user by setting the parameters \p npayload in the \p particles_list of the <tt>.par</tt> file. The default number of payloads is 1, and payload tracking can be suppressed by setting \p npayload=0.

  The code will by default update and store a number of payloads, depending on the running mode:
  * For the \p 'advect' mode, the fluid density at the particle location will be tracked and stored in the first payload.
  * For the \p 'Lorentz' mode, up to four payloads can be updated by default: the particle Lorentz factor (\p =1 if <tt>relativistic=.false.</tt>), the particle gyroradius, the magnetic moment, and the local value of \f$ \textbf{E}\cdot\textbf{B}\f$.
  * For the \p 'GCA' mode, there are 14 default payloads: 
    * Particle gyroradius;
    * Pitch angle \f$\tan^{-1}(v_\perp/v_{\|})\f$;
    * Perpendicular velocity \f$v_\perp\f$;
    * Four parallel acceleration terms (see right-hand side of the \f$du_{\|}/dt\f$ equation above);
    * Seven drift velocity terms (in magnitude; see right-hand side of the \f$d\textbf{R}/dt\f$ equation above).
  * For the \p 'sample' mode, by default (regardless of the value of \p npayload in the <tt>.par</tt> file) there will be a number of payloads \p n=nw, where \p nw is the number of variables in the fluid simulation. Each of these payloads samples one of the *primitive* fluid quantities, and therefore in the <tt>.csv</tt> output these payloads are named according to the names given to the primitive quantities.

  If the user wishes to define a custom payload update routine, this can be done in the \p mod_usr.t file. The user-defined routine <em>must</em> be associated with the \p usr_update_payload pointer at the beginning of \p mod_usr.t. The required format for a user-defined payload update routine (e.g. named \p update_payload) is:

  \code{.f90}
  subroutine update_payload(igrid,w,wold,xgrid,xpart,upart,qpart,mpart,payload,npayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,npayload
    double precision, intent(in)  :: w(ixG^T,1:nw),wold(ixG^T,1:nw)
    double precision, intent(in)  :: xgrid(ixG^T,1:ndim),xpart(1:ndir),upart(1:ndir),qpart,mpart,particle_time
    double precision, intent(out) :: payload(npayload)
  
  end subroutine update_payload
  \endcode

- **Boundary conditions**: by default, the boundary conditions for the particles are the same as the underlying (M)HD simulation. If a particle crosses a periodic boundary, it will be re-injected at the corresponding opposite side of the simulation box. If a particle crosses an open boundary, it will be "destroyed" (i.e., removed from the simulation). Destroyed particles are stored in a dedicated output file as the simulation progresses; see the "Output and visualisation" section below for more information on particle output.


## Running
Particle calculations in MPI-AMRVAC can be carried out in two different modes, namely i) concurrently to the (M)HD evolution, or ii) in a fixed (M)HD snapshot. Additionally, the convert operations that can be performed with (M)HD output files also affect the particle outputs. Each of these options is illustrated below.

- **Running the particle module along an HD/MHD simulation**: as mentioned above, after the (M)HD setup has been coded correctly in the \p mod_usr.t file and all parameters for the (M)HD simulation are set in the <tt>.par</tt> file, the particle module can be activated simply by including the statement <tt>hd_particles=.true.</tt> or <tt>mhd_particles=.true.</tt> in the \p mhd_list of the <tt>.par</tt> file. This will tell the code to perform the particle calculations according to the parameters specified in the \p particles_list block of the input file. The (M)HD calculation will behave as usual, and the particle results will be stored in the output according to dedicated parameters (see section "Output and visualisation" below).

- **Running the particle module in a fixed fluid snapshot**: this feature allows for running the particle module only, while keeping the (M)HD background fixed in time. This is accomplished by exploiting the \p convert functionality of the code. First, the user must set the parameters <tt>convert=.true., autoconvert=.false.</tt> in the \p filelist block of the <tt>.par</tt> file. Finally, the particle module must be activated via <tt>hd_particles=.true.</tt> or <tt>mhd_particles=.true.</tt> in the \p mhd_list. The particles will be initialised according to the user-defined (or the default) routines. Alternatively, an initial particle snapshot can be used as initial condition for the particles, if provided in the same folder of the fluid snapshot (see the "output and visualisation" section below for more info on particle outputs). In this mode, the particle integration will be carried over until \p time_max is reached; the code will assume that the initial time is the same as that stored inside the fluid snapshot. This implies that, if the snapshot was saved at time (say) \p t=9 and the user wishes to integrate particles in this snapshot for a total time of \p 1, they will have to specify \p time_max=10 in the <tt>.par</tt> file.

- **Restarting a run with particles**: restarting a run that includes particles is done in the usual way, by first selecting a restart file via the \p restart_from_file parameter in the <tt>.par</tt> file for the fluid. If <tt>hd_particles=.true.</tt> or <tt>mhd_particles=.true.</tt> in the same <tt>.par</tt> file, the code will search for a particle snapshot (<tt>.dat</tt> file) in the same directory of the fluid snapshot used for the restart (and with the same base file name - see next sections for more info on particle output). The particles will then be initialised according to the information stored in the particle snapshot. Note that if the particle snapshot is not found, then the code will restart the fluid calculations from the given fluid snapshot, and it will then initialise the particle module with the user-provided particle initialisation routines (or the default routines if the user-defined ones are not provided). In essence, this allows for either restarting a previously interrupted particle run or for starting a new particle run on top of a restarted fluid run.


## Output and visualisation
Below is a description of the various outputs associated to the use of the particle module, together with a brief introduction to the visualisation of the particle results.

- **Output files and formats**: the output of particle calculations (regardless of the chosen mode) consists of two file types, which will be stored in the same folder as the (M)HD results.

  <em>Particle snapshots</em> (format <tt>.dat</tt>) will be produced only when running the particle module along with a time-evolving (M)HD calculation, and will be produced with the same cadence of the (M)HD output. These files contain all particle information in raw binary format, and do not have any direct use beside providing a checkpoint for restarts. The base file name for particle snapshots will be constructed by concatenating the \p base_filename given in the <tt>.par</tt> file with an additional \p _particles and followed by the output number (same as the fluid <tt>.dat</tt> snapshots). Note that particle snapshots will *NOT* be produced when running the particle module in a fixed (M)HD snapshot.

  <em>Particle individual and ensemble outputs</em> (format <tt>.csv</tt>) can be used to easily analyse particle results. Whether or not these files are produced is controlled by the parameters \p write_individual and \p write_ensemble in the \p particles_list of the <tt>.par</tt> file. The output cadence to both individual and ensemble files is controlled via the \p dtsave_particles parameter in the same list. The output cadence may or may not be the same of the (M)HD output; the two can be specified completely independently.

  The standard output quantities stored in the <tt>.csv</tt> files are, for each particle:
  * Particle index (unique for each particle);
  * Current parent processor number (i.e. the processor on which the particle is found at that moment);
  * Current iteration;
  * Current time;
  * Current time step \p dt;
  * Particle position <tt>(x1,x2,x3)</tt>;
  * Particle velocity <tt>(u1,u2,u3)</tt>;
  * Payloads associated to the particle, labelled \p pl01, \p pl02, ..., \p plN (where \p N is the number of payloads specified via \p npayload).

  Note that when <tt>relativistic=.true.</tt> in the \p 'Lorentz' mode, the particle velocity will be replaced by the particle normalised momentum in the output files. In the \p 'GCA' mode, the quantities stored in \p (u1,u2,u3) are not the full particle velocity components, but rather the particle parallel velocity (or normalised parallel momentum) in \p u1 and the magnetic moment (or relativistic magnetic moment) in \p u2. The Lorentz factor will be stored in \p u3 if <tt>relativistic=.true.</tt>, otherwise \p u3=1 will be set by default.

  If <tt>write_individual=.true.</tt>, the code will produce one <tt>.csv</tt> file for each particle that the user flagged with <tt>follow=.true.</tt> at initialisation (by default, when the user does not provide a custom initial particle setup, all particles are flagged with <tt>follow=.false.</tt>). The base file name for individual particle <tt>.csv</tt>'s is obtained by concatenating \p base_filename with \p _particle_XXXXXX, where \p XXXXXX is a unique integer index (in \p %06d format) associated with each particle at initialisation. All information for a single individually followed particle will be stored over time, with cadence equal to \p dtsave_particles, inside the same <tt>.csv</tt> file, such that at the end of the run that file will contain the complete history of that single particle. Therefore, individual particle <tt>.csv</tt>'s can be used to visualise the trajectory and time evolution of the quantities associated to specific particles.

  If <tt>write_ensemble=.true.</tt>, the code will produce <em>ensemble</em> <tt>.csv</tt> files which, oppositely to individual <tt>.csv</tt> files, gather information from all particles at a specific time. A single ensemble <tt>.csv</tt> file will be produced at each output time specified via \p dtsave_particles. The base file name for ensemble <tt>.csv</tt>'s is obtained by concatenating \p base_filename with \p _ensemble_XXXXXX, where \p XXXXXX is a number (format \p %06d) corresponding to the \p n-th output time, based on the cadence specified via \p dtsave_particles. Ensemble <tt>.csv</tt>'s are useful to visualise all particles together at a specific moment in time.

  As a practical example, suppose the user has chosen \p dtsave_particles such that 10 particle outputs will be produced during a simulation. Suppose further that both <tt>write_individual=.true.</tt> and <tt>write_ensemble=.true.</tt>, \p nparticles=100, and the user has flagged particles with index 36, 47, and 99 with <tt>follow=.true.</tt> at initialisation. The <tt>.csv</tt> particle files that will be found in the output folder will then be:

  * 11 ensemble <tt>.csv</tt> files (one for each \p dtsave_particles, plus the initial state) containing 100 rows each (one row for each particle);
  * 3 individual <tt>.csv</tt> files (one for each of the three followed particles) containing 11 rows each (one row for each output time plus the initial state).

  An additional <tt>.csv</tt> file, containing the \p _destroyed label, may be present in the output folder. In this file, "destroyed" particles (i.e. particles removed from the domain) are stored as the simulation progresses.


## Additional options
Additional parameters in the \p particles_list of the .par file are available for refining the particle integration process:
- \p relativistic: if <tt>.true.</tt>, the relativistic equations of motion will be solved, instead of the Newtonian ones, when \p physics_type_particles='Lorentz' or \p \p physics_type_particles='GCA'.
- \p integrator_type_particles: if \p physics_type_particles='Lorentz' and \p relativistic=<tt>.true.</tt>, the user can set this parameter to the preferred particle integrator, choosing among \p 'Boris' (for the Boris integrator), \p 'Vay' (for the Vay integrator), \p 'HC' (for the Higuera-Cary integrator), or \p 'LM' (for the Lapenta-Markidis integrator).
- \p eta_particles, \p etah_particles: when \p physics_type_particles='Lorentz' or \p physics_type_particles='GCA', the user can replace the resistivity or the hall resistivity from the MHD with a different (constant) resistivity which will be employed for calculating the electric field that enters the particle equations of motion.
- \p const_dt_particles: the user can define a constant time step for the particle integration, which will be used instead of calculating the time step from standard procedures.


# Contact details

For more specific questions, email [Fabio Bacchini](mailto:fabio.bacchini@kuleuven.be), [Bart Ripperda](mailto:bart.ripperda@gmail.com), or [Jannis Teunissen](mailto:jannis.teunissen@cwi.nl).

