!> This module contains definitions of global parameters and variables and some
!> generic functions/subroutines used in AMRVAC.
!>
module mod_global_parameters
  use mpi
  use mod_constants
  use mod_basic_types
  use mod_physicaldata
  use mod_connectivity
  use mod_variables

  implicit none
  public

  !> minimum and maximum domain boundaries for each dimension
  double precision  :: xprob^L
  !> store unstretched cell size of current level
  double precision :: dxlevel(^ND)
  !> stretch factor between cells at AMR level 1, per dimension
  double precision ::  qstretch_baselevel(^ND)

  !> physical extent of stretched border in symmetric stretching
  double precision :: xstretch^D

  !> to monitor timeintegration loop at given wall-clock time intervals
  double precision :: time_between_print

  !> accumulated wall-clock time spent on boundary conditions
  double precision :: time_bc

  !> global time step
  double precision :: dt

  !> The Courant (CFL) number used for the simulation
  double precision :: courantpar


  !> If dtpar is positive, it sets the timestep dt, otherwise courantpar is used
  !> to limit the time step based on the Courant condition.
  double precision :: dtpar

  !> For resistive MHD, the time step is also limited by the diffusion time:
  !> \f$ dt < dtdiffpar \times dx^2/eta \f$
  double precision :: dtdiffpar

  !> The global simulation time
  double precision :: global_time

  !> Start time for the simulation
  double precision :: time_init

  !> End time for the simulation
  double precision :: time_max

  !> Ending wall time (in hours) for the simulation
  double precision :: wall_time_max

  !> Stop the simulation when the time step becomes smaller than this value
  double precision :: dtmin

  !> Conversion factor for length unit
  double precision :: length_convert_factor

  !> Conversion factor for time unit
  double precision :: time_convert_factor

  !> Fix the AMR grid after this time
  double precision :: tfixgrid

  ! Physics factors
  !> Physical scaling factor for length
  double precision :: unit_length=1.d0

  !> Physical scaling factor for time
  double precision :: unit_time=1.d0

  !> Physical scaling factor for density
  double precision :: unit_density=1.d0

  !> Physical scaling factor for velocity
  double precision :: unit_velocity=1.d0

  !> Physical scaling factor for temperature
  double precision :: unit_temperature=1.d0

  !> Physical scaling factor for pressure
  double precision :: unit_pressure=1.d0

  !> Physical scaling factor for magnetic field
  double precision :: unit_magneticfield=1.d0

  !> Physical scaling factor for number density
  double precision :: unit_numberdensity=1.d0

  !> Physical scaling factor for charge
  double precision :: unit_charge=1.d0

  !> Physical scaling factor for mass
  double precision :: unit_mass=1.d0

  !> Normalised speed of light
  double precision :: c_norm=1.d0

  !> Normalised radiation constant
  double precision :: arad_norm=1.d0

  !> Physical scaling factor for Opacity
  double precision :: unit_opacity=1.d0

  !> Physical scaling factor for radiation flux
  double precision :: unit_radflux=1.d0

  !> Physical scaling factor for radiation energy density
  double precision :: unit_Erad=1.d0

  !> Physical factors useful for radiation fld
  double precision :: const_rad_a,const_kappae,const_sigmaSB

  !> error handling
  double precision :: small_temperature,small_pressure,small_density,small_e,small_r_e
  double precision :: phys_trac_mask

  double precision :: minfip=1.d0, maxfip=4.d0

  !> amplitude of background dipolar, quadrupolar, octupolar, user's field
  double precision :: Bdip=0.d0
  double precision :: Bquad=0.d0
  double precision :: Boct=0.d0
  double precision :: Busr=0.d0
  !> Stores the memory and load imbalance, used in printlog
  double precision :: Xload, Xmemory
  double precision :: tvdlfeps

  !> RK2(alfa) method parameters from Butcher tableau
  double precision :: rk_a21,rk_b1,rk_b2
  !> IMEX-222(lambda) one-parameter family of schemes
  double precision :: imex222_lambda
  double precision :: rk_beta11,rk_beta22,rk_beta33,rk_beta44,rk_c2,rk_c3,rk_c4
  double precision :: rk_alfa21,rk_alfa22,rk_alfa31,rk_alfa33,rk_alfa41,rk_alfa44
  double precision :: rk_beta54,rk_beta55,rk_alfa53,rk_alfa54,rk_alfa55,rk_c5
  double precision :: rk3_a21,rk3_a31,rk3_a32,rk3_b1,rk3_b2,rk3_b3,rk3_c2,rk3_c3
  !> IMEX_ARS3 parameter ars_gamma
  double precision :: ars_gamma
  double precision :: imex_a21,imex_a31,imex_a32,imex_b1,imex_b2,imex_ha21,imex_ha22
  double precision :: imex_b3,imex_c2,imex_c3
  !> IMEX_CB3a extra parameters
  double precision :: imex_a22, imex_a33, imex_ha32

  !> global fastest wave speed needed in fd scheme and glm method
  double precision :: cmax_global

  !> global fastest flow speed needed in glm method
  double precision :: vmax_global

  !> times for enhancing spatial resolution for EUV image/spectra
  double precision :: instrument_resolution_factor
  !> the white light emission below it (unit=Rsun) is not visible 
  double precision :: R_occultor
  !> direction of the line of sight (LOS)
  double precision :: LOS_theta,LOS_phi
  !> rotation of image
  double precision :: image_rotate
  !> where the is the origin (X=0,Y=0) of image
  double precision :: x_origin(1:3)
  !> Base file name for synthetic EUV spectrum output
  character(len=std_len) :: filename_spectrum
  !> Synthetic emission transfer mode: thin or thick
  character(len=std_len) :: radiation_transfer
  !> Synthetic emission ray traversal method
  character(len=std_len) :: ray_method
  !> Data-resolution image spacing: nominal or minimum actual cell size
  character(len=std_len) :: dat_resolution_mode
  !> Synthetic emission physical model selector
  character(len=std_len) :: emission_model
  !> Post-process dat-resolution EUV images onto the instrument pixel grid
  logical :: instrument_postprocess
  !> Observing frequency for radio free-free synthesis in Hz
  double precision :: radio_frequency
  !> Gaussian radio beam full width at half maximum in arcsec
  double precision :: radio_beam_fwhm
  !> Output pixel size for radio beam post-processing in arcsec; <=0 uses FWHM/3
  double precision :: radio_beam_pixel_size
  !> Number of image pixels processed in one ray-segment MPI batch
  integer :: radsyn_pixel_batch
  !> Maximum ray segments per pixel batch, as a factor of radsyn_pixel_batch; <=0 uses memory budget.
  !> Default value is set from the best spherical u512 benchmark behavior seen so far.
  integer :: radsyn_segment_batch_factor
  !> Approximate per-rank temporary memory budget, in MiB, for automatic ray-segment batch sizing
  double precision :: radsyn_segment_memory_mb
  !> Maximum ray segments per segmented MPI all-to-all round, as a factor of radsyn_pixel_batch
  integer :: radsyn_segment_comm_factor
  !> Print synthetic-emission ray-tracing profiling counters
  logical :: radsyn_verbose
  !> spectral window
  double precision :: spectrum_window_min,spectrum_window_max
  !> location of the slit
  double precision :: location_slit
  !> for spherical coordinate, region below it (unit=Rsun) is treated as not transparent
  double precision :: R_opt_thick

  !> domain percentage cut off shifted from each boundary when converting data
  double precision :: writespshift(^ND,2)

  double precision, allocatable :: entropycoef(:)

  !> extent of grid blocks in domain per dimension, in array over levels
  double precision, dimension(:), allocatable :: dg^D
  !> Corner coordinates
  double precision, allocatable :: rnode(:,:)
  double precision, allocatable :: rnode_sub(:,:)

  !> spatial steps for all dimensions at all levels
  double precision, allocatable :: dx(:,:)

  !> Error tolerance for refinement decision
  double precision, allocatable :: refine_threshold(:)
  !> Error tolerance ratio for derefinement decision
  double precision, allocatable :: derefine_ratio(:)

  !> Stretching factors and first cell size for each AMR level and dimension
  double precision, allocatable :: qstretch(:,:), dxfirst(:,:),  &
                                   dxfirst_1mq(:,:), dxmid(:,:)

  !> Conversion factors the primitive variables
  double precision, allocatable :: w_convert_factor(:)

  !> Weights of variables used to calculate error for mesh refinement
  double precision, allocatable :: w_refine_weight(:)

  !> refinement: lohner estimate wavefilter setting
  double precision, allocatable :: amr_wavefilter(:)

  !> Maximum number of saves that can be defined by tsave or itsave
  integer, parameter :: nsavehi=100
  !> Number of output methods
  integer, parameter :: nfile = 5

  !> Save output of type N on times tsave(:, N)
  double precision :: tsave(nsavehi,nfile)

  double precision :: tsavelast(nfile)

  !> Repeatedly save output of type N when dtsave(N) simulation time has passed
  double precision :: dtsave(nfile)

  !> Start of read out (not counting specified read outs)
  double precision :: tsavestart(nfile)

  !> Number of spatial dimensions for grid variables
  integer, parameter :: ndim=^ND

  !> The number of MPI tasks
  integer :: npe

  !> The rank of the current MPI task
  integer :: mype

  !> The MPI communicator
  integer :: icomm

  !> A global MPI error return code
  integer :: ierrmpi

  !> Per-rank load-balance timing diagnostic toggle (off by default).
  !> When .true., per-rank wall times are gathered every step and written
  !> to <basename>_rank_timing.log.
  logical :: lb_diagnose = .false.

  !> Cost-weighted automatic load balancer toggle (off by default).
  !> When .true., the SFC partitioner cuts on cumulative measured per-block
  !> cost rather than on equal block counts.
  logical :: lb_automatic = .false.

  !> Rebalance every lb_interval cycles when lb_automatic is on.
  integer :: lb_interval = 10

  !> Exponential-moving-average decay for the per-block cost.
  !> costlist <- lb_alpha*costlist + (1-lb_alpha)*measured.
  !> 0 = no memory (volatile), 1 = no update.
  double precision :: lb_alpha = 0.9d0

  !> Bound on the memory imbalance get_Morton_range_costed may create: no rank
  !> may hold more than lb_max_block_ratio * (nleafs/npe) blocks. Cost balance
  !> and memory balance are in tension, as a rank owning cheap blocks must hold
  !> more of them to match one owning expensive blocks. 1.0 is the equal-block
  !> cut; larger values allow more freedom to follow cost.
  !> get_Morton_range_active bounds the same ratio through its wa/wp weights.
  double precision :: lb_max_block_ratio = 1.25d0

  !> Per-step per-block (per-rank, indexed by igrid) cost accumulator.
  !> Reset at start of each advance call; filled inside iigrid loops by the
  !> per-block timer wrappers in mod_advance.t and mod_supertimestepping.t.
  double precision, dimension(:), allocatable :: block_cost

  !> Per-step per-block cost of the short-characteristics sweep, indexed by
  !> igrid. Separate from block_cost because rt_sc_solve() runs before advance(),
  !> which zeroes block_cost on entry. advance() seeds block_cost from this array
  !> and clears it for the next solve, so both costs share a partition weight.
  double precision, dimension(:), allocatable :: block_cost_rt

  !> Per-rank wall time spent in the SC radiative-transfer KBA sweep this step,
  !> for the lb_diagnose rank-timing log. Filled by mod_rt_sc_solvers (which runs
  !> before advance), gathered and zeroed inside advance's diagnostic block.
  double precision :: lb_rt_accum = 0.0d0

  !> Persistent global per-Morton-leaf EWMA cost, in seconds (block_cost is a
  !> wall-clock accumulation). Sized to max_blocks*npe (the upper bound on
  !> nleafs); only the first nleafs entries are meaningful. Morton numbering
  !> is invariant across load_balance migration, so costlist values are
  !> correctly preserved when blocks change rank. Refinement events
  !> insert/remove Morton indices, after which the EWMA recovers over
  !> ~5 cycles. Indexed 1..nleafs. Updated by EWMA blend inside
  !> get_Morton_range_costed.
  double precision, dimension(:), allocatable :: costlist

  !> Whether the matching costlist slot has ever carried a measurement. Morton
  !> indices created by a refinement event have not been timed; they are seeded
  !> with the mean of the measured slots, since costlist is a wall time and any
  !> fixed seed would be a units error.
  logical, dimension(:), allocatable :: costlist_seen

  !> MPI file handle for logfile
  integer :: log_fh
  !> MPI type for block including ghost cells and its size
  integer :: type_block, size_block
  !> MPI type for block coarsened by 2, and for its children blocks
  integer :: type_coarse_block, type_sub_block(2^D&)
  !> MPI type for staggered block coarsened by 2, and for its children blocks
  integer :: type_coarse_block_stg(^ND,2^D&), type_sub_block_stg(^ND,2^D&)
  !> MPI type for IO: block excluding ghost cells
  integer :: type_block_io, size_block_io
  !> MPI type for IO of staggered variables
  integer :: type_block_io_stg, size_block_io_stg
  !> MPI type for IO: cell corner (xc) or cell center (xcc) coordinates
  integer :: type_block_xc_io,type_block_xcc_io
  !> MPI type for IO: cell corner (wc) or cell center (wcc) variables
  integer :: type_block_wc_io,type_block_wcc_io

  ! geometry and domain setups
  !> the mesh range of a physical block without ghost cells
  integer :: ixM^LL

  !> Indices for cylindrical coordinates FOR TESTS, negative value when not used:
  integer :: r_ = -1
  integer :: phi_ = -1
  integer :: z_ = -1

  !> Number of spatial dimensions (components) for vector variables
  integer :: ndir=ndim

  !> starting dimension for electric field
  {^IFONED
  integer, parameter :: sdim=3
  }
  {^IFTWOD
  integer, parameter :: sdim=3
  }
  {^IFTHREED
  integer, parameter :: sdim=1
  }

  !> number of cells for each dimension in level-one mesh
  integer :: domain_nx^D

  !> number of cells for each dimension in grid block excluding ghostcells
  integer :: block_nx^D

  !> Lower index of grid block arrays (always 1)
  integer, parameter :: {ixGlo^D = 1|, }

  !> Upper index of grid block arrays
  integer :: ixGhi^D

  !> Lower index of stagger grid block arrays (always 0)
  integer, parameter :: {ixGslo^D = 0|, }

  !> Upper index of stagger grid block arrays
  integer :: ixGshi^D

  !> Number of ghost cells surrounding a grid
  integer :: nghostcells = 2

  integer, parameter :: stretch_none = 0 !< No stretching
  integer, parameter :: stretch_uni  = 1 !< Unidirectional stretching from a side
  integer, parameter :: stretch_symm = 2 !< Symmetric stretching around the center

  !> What kind of stretching is used per dimension
  integer :: stretch_type(ndim)
  !> (even) number of (symmetrically) stretched
  !> blocks at AMR level 1, per dimension
  integer ::  nstretchedblocks_baselevel(ndim)
  !> (even) number of (symmetrically) stretched blocks per level and dimension
  integer, allocatable ::  nstretchedblocks(:,:)

  !> grid hierarchy info (level and grid indices)
  integer, parameter :: nodehi=^ND+1
  integer, parameter :: plevel_=1
  integer, parameter :: pig^D_=plevel_+^D

  integer, allocatable :: node(:,:)
  integer, allocatable :: node_sub(:,:)

  !> grid location info (corner coordinates and grid spacing)
  integer, parameter :: rnodehi=3*^ND
  integer, parameter :: rpxmin0_=0
  integer, parameter :: rpxmin^D_=rpxmin0_+^D
  integer, parameter :: rpxmax0_=^ND
  integer, parameter :: rpxmax^D_=rpxmax0_+^D
  integer, parameter :: rpdx^D_=2*^ND+^D

  !> index number of the latest existing data file
  integer :: index_latest_data

  !> Save output of type N on iterations itsave(:, N)
  integer :: itsave(nsavehi,nfile)

  integer :: itsavelast(nfile)

  !> Repeatedly save output of type N when ditsave(N) time steps have passed
  integer :: ditsave(nfile)

  integer :: isavet(nfile)

  integer :: isaveit(nfile)

  !> The level at which to produce line-integrated / collapsed output
  integer :: collapseLevel

  !> Number of saved files of each type
  integer :: n_saves(1:nfile)

  !> IO: snapshot and collapsed views output numbers/labels
  integer :: snapshotnext, collapsenext

  !> IO: slice output number/label
  integer :: slicenext

  !> Constant indicating log output
  integer, parameter :: filelog_      = 1

  !> Constant indicating regular output
  integer, parameter :: fileout_      = 2

  !> Constant indicating slice output
  integer, parameter :: fileslice_    = 3

  !> Constant indicating collapsed output
  integer, parameter :: filecollapse_ = 4

  !> Constant indicating analysis output (see @ref analysis.md)
  integer, parameter :: fileanalysis_ = 5

  !> Unit for standard input
  integer, parameter :: unitstdin=5

  !> Unit for standard output
  integer, parameter :: unitterm=6

  !> Unit for error messages
  integer, parameter :: uniterr=6

  !> file handle for IO
  integer, parameter :: unitpar=9
  integer, parameter :: unitconvert=10
  integer, parameter :: unitslice=11
  integer, parameter :: unitsnapshot=12
  integer, parameter :: unitcollapse=13
  integer, parameter :: unitanalysis=14

  !> Number of auxiliary variables that are only included in output
  integer :: nwauxio

  integer :: level_io, level_io_min, level_io_max, type_endian

  !> Resume from the snapshot with this index
  integer :: snapshotini

  !> number of equilibrium set variables, besides the mag field
  integer :: number_equi_vars = 0

  integer :: phys_trac_type=1
  integer :: phys_trac_finegrid=4
  integer :: phys_trac_nzones=1
  double precision :: phys_trac_zone_splits(10)=-1.d0

  !> integer switchers for type courant
  integer, parameter :: type_maxsum=1
  integer, parameter :: type_summax=2
  integer, parameter :: type_minimum=3

  ! AMR switches
  !> The maximum number of grid blocks in a processor
  integer :: max_blocks

  !> The maximum number of levels in the grid refinement
  integer, parameter :: nlevelshi = 20

  !> Maximal number of AMR levels
  integer :: refine_max_level

  !> Fix the AMR grid after this many time steps
  integer :: itfixgrid

  !> Reconstruct the AMR grid once every ditregrid iteration(s)
  integer :: ditregrid

  !> select types of refine criterion
  integer :: refine_criterion

  !> Number of cells as buffer zone
  integer :: nbufferx^D

  integer :: levmin
  integer :: levmax
  integer :: levmax_sub

  ! Miscellaneous
  !> problem switch allowing different setups in same usr_mod.t
  integer           :: iprob

  !> Kronecker delta tensor
  integer :: kr(3,3)

  !> Levi-Civita tensor
  integer :: lvc(3,3,3)

  !> How to compute the CFL-limited time step
  integer :: type_courant=1

  !> Number of time steps taken
  integer :: it

  !> Stop the simulation after this many time steps have been taken
  integer :: it_max

  !> initial iteration count
  integer :: it_init

  !> If > 1, then in the first slowsteps-1 time steps dt is reduced
  !> by a factor \f$ 1 - (1- step/slowsteps)^2 \f$
  integer :: slowsteps

  ! Method switches

  !> Index of the sub-step in a multi-step time integrator
  integer :: istep

  !> How many sub-steps the time integrator takes
  integer :: nstep

  !> flux schemes
  integer, parameter :: fs_hll=1
  integer, parameter :: fs_hllc=2
  integer, parameter :: fs_hlld=3
  integer, parameter :: fs_hllcd=4
  integer, parameter :: fs_tvdlf=5
  integer, parameter :: fs_tvdmu=6
  integer, parameter :: fs_tvd=7
  integer, parameter :: fs_hancock=8
  integer, parameter :: fs_cd=9
  integer, parameter :: fs_cd4=10
  integer, parameter :: fs_fd=11
  integer, parameter :: fs_source=12
  integer, parameter :: fs_nul=13

  !> time stepper type
  integer :: t_stepper=0
  integer, parameter :: onestep=1
  integer, parameter :: twostep=2
  integer, parameter :: threestep=3
  integer, parameter :: fourstep=4
  integer, parameter :: fivestep=5

  !> time integrator method
  integer :: t_integrator=0
  integer, parameter :: Forward_Euler=1
  integer, parameter :: Predictor_Corrector=2
  integer, parameter :: ssprk3=3
  integer, parameter :: ssprk4=4
  integer, parameter :: ssprk5=5

  integer, parameter :: IMEX_Euler=6
  integer, parameter :: IMEX_SP=7
  integer, parameter :: RK2_alf=8
  integer, parameter :: ssprk2=9
  integer, parameter :: IMEX_Midpoint=10
  integer, parameter :: IMEX_Trapezoidal=11
  integer, parameter :: IMEX_222=12

  integer, parameter :: RK3_BT=13
  integer, parameter :: IMEX_ARS3=14
  integer, parameter :: IMEX_232=15
  integer, parameter :: IMEX_CB3a=16

  integer, parameter :: rk4=17

  !> number of grid blocks in domain per dimension, in array over levels
  integer, dimension(:), allocatable :: ng^D

  !> Which flux scheme of spatial discretization to use (per grid level)
  integer, allocatable :: flux_method(:)

  !> The spatial discretization for the predictor step when using a two
  !> step PC method
  integer, allocatable :: typepred1(:)

  !> Type of slope limiter used for reconstructing variables on cell edges
  integer, allocatable :: type_limiter(:)

  !> Type of slope limiter used for computing gradients or divergences, when
  !> typegrad or typediv are set to 'limited'
  integer, allocatable :: type_gradient_limiter(:)

  !> Coefficient alpha of the explicit diffusive flux of Mignone et al. 2005 eq. B19-B21,
  !> the second half of PPM's dissipation algorithm. 0 = off (default); the paper says
  !> "alpha is typically set to 0.1".
  !>
  !>     F_{i+1/2} -> F_{i+1/2} + k_nu (U_i - U_{i+1}) ,   k_nu = alpha*max(-D_{i+1/2},0)
  !>
  !> D is an undivided multidimensional divergence of v at the interface, so the term is
  !> active only where the flow converges. Two consequences worth knowing:
  !>   * it is identically zero for a state at rest, so it cannot disturb a hydrostatic
  !>     column -- exact discrete balance survives it untouched;
  !>   * it is proportional to the raw cell difference U_i - U_{i+1}, not to the
  !>     reconstructed face jump, so unlike everything else in the scheme it can see an
  !>     odd-even mode that a face-value reconstruction annihilates.
  double precision :: ppm_avisc = 0.0d0

  !> Residual-jump velocity reconstruction (RJV). 0 = off (default), 1 = full strength.
  !>
  !> PPM's face operator annihilates the odd-even (2 dx) mode: the face value's leading
  !> term is (q_i + q_{i+1})/2, and adjacent cells carry opposite ripple sign, while the
  !> slope depends on (q_{i+1} - q_{i-1})/2, where the two cells have the same parity. So
  !> the mode leaves no trace in w_L or w_R, the Riemann solver sees no jump, and applies
  !> exactly zero dissipation to it. It never decays once seeded.
  !>
  !> But what the operator annihilates is exactly recoverable as its own residual,
  !>
  !>     e_i = v_i - (P_{i+1/2} + P_{i-1/2})/2 ,
  !>
  !> built from the two face values PPM has already computed -- no new stencil. For a
  !> ripple A(-1)^i, e_i = A(-1)^i at full amplitude; for smooth data e_i = -h^2 v''/8.
  !> Reinject it as a jump only, leaving the face value (and hence the flux, and hence a
  !> hydrostatic column) untouched:
  !>
  !>     w_L(i+1/2) = P_{i+1/2} + beta*s_i*e_i ,  w_R(i+1/2) = P_{i+1/2} + beta*s_{i+1}*e_{i+1}
  !>
  !> with s_i = 1 only where e alternates across three cells, which no smooth field does.
  !> At beta = 1 this reproduces the piecewise-constant jump for a pure ripple, i.e. full
  !> first-order dissipation on the null component alone.
  !>
  !> Applied to the velocity components only. That is what makes it exact rather than
  !> approximate: at rest v = 0, so e = 0 identically and the jump is exactly zero, and
  !> every flux term the velocity enters carries a factor v. Measured inert on a resolved
  !> sine (l=8 and l=32), a 3-cell tanh front, a discontinuity, and an exponential column.
  double precision :: ppm_rjv = 0.0d0

  !> background magnetic field location indicator
  integer :: b0i=0

  !> Limiter used for prolongation to refined grids and ghost cells
  integer :: prolong_limiter=0

  !> bound (left/min and right.max) speed of Riemann fan
  integer :: boundspeed

  integer                       :: nxdiffusehllc
  !> SSPRK choice of methods (both threestep and fourstep, Shu-Osher 2N* implementation)
  !> also fivestep SSPRK54
  integer                       :: ssprk_order
  !> RK3 Butcher table
  integer                       :: rk3_switch
  !> IMEX_232 choice and parameters
  integer                       :: imex_switch

  !> Array indicating the type of boundary condition per variable and per
  !> physical boundary
  integer, allocatable :: typeboundary(:, :)
  !> boundary condition types
  integer, parameter :: bc_special=1
  integer, parameter :: bc_cont=2
  integer, parameter :: bc_symm=3
  integer, parameter :: bc_asymm=4
  integer, parameter :: bc_periodic=5
  integer, parameter :: bc_aperiodic=6
  integer, parameter :: bc_noinflow=7
  integer, parameter :: bc_data=8
  integer, parameter :: bc_character=9
  integer, parameter :: bc_icarus=10

  !> wavelength for output
  integer :: wavelength
  ! minimum and maximum energy of SXR (keV)
  integer :: emin_sxr,emax_sxr
  !> wave length for spectrum
  integer :: spectrum_wl
  !> direction of the slit (for dat resolution only)
  integer :: direction_slit


  !> Cartesian geometry or not
  logical :: slab

  !> uniform Cartesian geometry or not (stretched Cartesian)
  logical :: slab_uniform

  !> each cell has its own timestep or not
  logical :: local_timestep = .false.

  !> whether or not to save an output file
  logical :: save_file(nfile)

  !> If true, adjust mod_geometry routines to account for grid stretching (but
  !> the flux computation will not)
  logical :: stretch_uncentered
  !> True if a dimension is stretched
  logical :: stretched_dim(ndim)

  !> If true, restart a previous run from the latest snapshot
  logical :: resume_previous_run

  !> If true and restart_from_file is given, convert snapshots to
  !> other file formats
  logical :: convert

  !> If true, already convert to output format during the run
  logical :: autoconvert

  !> If true, convert from conservative to primitive variables in output
  logical :: saveprim

  !> do time evolving
  logical :: time_advance

  !> Force timeloop exit when final dt < dtmin
  logical :: final_dt_exit

  !> If true, reset iteration count and global_time to original values, and
  !> start writing snapshots at index 0
  logical :: reset_time

  !> If true, reset iteration count to 0
  logical :: reset_it

  !> If true, allow final dt reduction for matching time_max on output
  logical :: final_dt_reduction

  !> If true, call initonegrid_usr upon restarting
  logical :: firstprocess

  !> If true, allow a restart from a snapshot whose ndir differs from the current
  !> run (e.g. hd ndir=2 -> mhd ndir=3). Block I/O is ndim-based only, so the data
  !> loads exactly; use usr_transform_w to place the source vars into the right slots.
  logical :: allow_ndir_change

  !> If true, wall time is up, modify snapshotnext for later overwrite
  logical :: pass_wall_time

  !> If true, do H-correction to fix the carbuncle problem at grid-aligned shocks
  logical :: H_correction=.false.

  !> If true, rebuild the AMR grid upon restarting
  logical :: reset_grid
  !> True for using stagger grid
  logical :: stagger_grid=.false.
  !> Limit outstanding ghost-cell send requests to avoid MPI request pressure.
  !> In large-scale MHD tests, the original all-at-once ghost-cell exchange could
  !> stall during spherical u512 initialization. Batching avoids
  !> that failure mode, while u256 5-step tests showed less than 0.1% total
  !> runtime difference. The original path is still available by setting this to
  !> false; if batching remains stable in broader tests, the old path can be
  !> removed later.
  logical :: ghostcell_comm_batched=.true.
  integer :: ghostcell_comm_batch_size=64
  !> True for record electric field
  logical :: record_electric_field=.false.

  !> resolution of the images
  logical :: dat_resolution
  !> output optical-depth map for synthetic emission when available
  logical :: output_tau
  !> output absorption fraction for thick/thin EUV synthesis when available
  logical :: output_absorption_fraction

  !> If collapse(DIM) is true, generate output integrated over DIM
  logical :: collapse(ndim)
  !> IO switches for conversion
  logical          :: nocartesian

  !> Use particles module or not
  logical :: use_particles=.false.

  !> Use multigrid (only available in 2D and 3D)
  logical :: use_multigrid = .false.

  !> prolongate primitive variables in level-jump ghost cells
  logical                       :: prolongprimitive=.false.

  !> coarsen primitive variables in level-jump ghost cells
  logical                       :: coarsenprimitive=.false.

  !> Save a snapshot before crash a run met unphysical values
  logical :: crash=.false.

  !> check and optionally fix unphysical small values (density, gas pressure)
  logical :: check_small_values=.true.

  !> fix small values with average or replace methods
  logical :: fix_small_values=.false.

  !> split magnetic field as background B0 field
  logical :: B0field=.false.
  logical :: B0fieldAllocCoarse=.false.
  !> Use SI units (.true.) or use cgs units (.false.)
  logical :: SI_unit=.false.

  !> Use TRAC for MHD or 1D HD
  logical :: phys_trac=.false.

  !> Use escape probability for radiative cooling modification
  logical :: phys_escape_prob=.false.

  !> Whether to apply flux conservation at refinement boundaries
  logical :: fix_conserve_global = .true.
  logical :: flux_adaptive_diffusion
  double precision :: flux_adaptive_diffusion_min
  double precision :: flux_adaptive_diffusion_scale
  logical :: flux_energy_only
  logical :: flathllc,flatcd,flatsh
  !> Use split or unsplit way to add user's source terms, default: unsplit
  logical :: source_split_usr
  !> if any normal source term is added in split fasion
  logical :: any_source_split=.false.
  logical :: dimsplit
  !> whether IMEX in use or not
  logical :: use_imex_scheme

  !> need global maximal wave speed
  logical :: need_global_cmax=.false.

  ! Boundary region parameters

  !> True for dimensions with periodic boundaries
  logical :: periodB(ndim)

  !> Indicates whether there is a pole at a boundary
  logical :: poleB(2,ndim)

  !> True for dimensions with aperiodic boundaries
  logical :: aperiodB(ndim)

  !> True for save physical boundary cells in dat files
  logical :: save_physical_boundary

  !> whether copy values instead of interpolation in ghost cells of finer blocks
  logical :: ghost_copy=.false.

  !> if there is an internal boundary
  logical :: internalboundary

  !> use arcsec as length unit of images/spectra
  logical :: activate_unit_arcsec
  !> big image
  logical :: big_image

  !> True if a block has any physical boundary
  logical, allocatable :: phyboundblock(:)

  !> if true write the w variable in output
  logical, allocatable :: w_write(:)

  logical, allocatable :: writelevel(:)

  logical, allocatable :: loglimit(:), logflag(:)

  ! Parameters
  character(len=*), parameter :: undefined = 'undefined'

  !> Names of the output methods
  character(len=40), parameter  :: output_names(nfile) = &
       ['log      ', 'normal   ', 'slice    ', 'collapsed', 'analysis ']
  !> Which format to use when converting
  !>
  !> Options are: tecplot, tecplotCC, vtu, vtuCC, vtuB, vtuBCC,
  !> tecplotmpi, tecplotCCmpi, vtumpi, vtuCCmpi, vtuBmpi, vtuBCCmpi, pvtumpi, pvtuCCmpi,
  !> pvtuBmpi, pvtuBCCmpi, tecline, teclinempi, onegrid
  character(len=std_len) :: convert_type

  character(len=std_len) :: collapse_type

  !> User parameter file
  character(len=std_len)   :: usr_filename

  !> Base file name for simulation output, which will be followed by a number
  character(len=std_len) :: base_filename

  !> If not 'unavailable', resume from snapshot with this base file name
  character(len=std_len) :: restart_from_file

  !> Which type of log to write: 'normal', 'special', 'regression_test'
  character(len=std_len) :: typefilelog

  character(len=std_len) :: typeaverage
  character(len=std_len) :: typedimsplit
  character(len=std_len) :: geometry_name='default'
  character(len=std_len) :: typepoly
  !> Base file name for synthetic EUV emission output
  character(len=std_len) :: filename_euv
  !> Base file name for synthetic SXR emission output
  character(len=std_len) :: filename_sxr
  !> Base file name for synthetic white light
  character(len=std_len) :: filename_whitelight
  !> white light observation instrument
  character(len=std_len) :: whitelight_instrument

  !> Which type of TVD method to use
  character(len=std_len) :: typetvd

  character(len=std_len) :: typediv,typegrad

  !> Which par files are used as input
  character(len=std_len), allocatable :: par_files(:)

  !> Which type of entropy fix to use with Riemann-type solvers
  character(len=std_len), allocatable :: typeentropy(:)

  !> Block pointer for using one block and its previous state
  type(state), pointer :: block

  !$OMP THREADPRIVATE(block,dxlevel,b0i)

contains

  !> Cross product of two vectors
  pure subroutine cross_product(ixI^L,ixO^L,a,b,axb)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: a(ixI^S,3), b(ixI^S,3)
    double precision, intent(out) :: axb(ixI^S,3)

    axb(ixO^S,1)=a(ixO^S,2)*b(ixO^S,3)-a(ixO^S,3)*b(ixO^S,2)
    axb(ixO^S,2)=a(ixO^S,3)*b(ixO^S,1)-a(ixO^S,1)*b(ixO^S,3)
    axb(ixO^S,3)=a(ixO^S,1)*b(ixO^S,2)-a(ixO^S,2)*b(ixO^S,1)
  end subroutine cross_product

end module mod_global_parameters
