#:mute
#:include "mod_gpu_directives.fpp"
#:endmute
!> This module contains definitions of global parameters and variables and some
!> generic functions/subroutines used in AMRVAC.
!>
module mod_global_parameters
  use mod_physicaldata
  use mod_connectivity
  use mpi
  use mod_constants
  use mod_variables
  use mod_basic_types

  implicit none
  public

  ! Parameters
  character(len=*), parameter :: undefined = 'undefined'

  !> The number of MPI tasks
  integer :: npe
  ${GPU_DECLARE_CREATE('npe')}$

  !> The rank of the current MPI task
  integer :: mype
  ${GPU_DECLARE_CREATE('mype')}$

  !> The MPI communicator
  integer :: icomm

  !> A global MPI error return code
  integer :: ierrmpi

  !> MPI file handle for logfile
  integer :: log_fh
  !> MPI type for block including ghost cells and its size
  integer :: type_block, size_block
  !> MPI type for block coarsened by 2, and for its children blocks
  integer :: type_coarse_block, type_sub_block(2,2,2)
  !> MPI type for staggered block coarsened by 2, and for its children blocks
  integer :: type_coarse_block_stg(3,2,2,2), type_sub_block_stg(3,2,2,2)
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
  integer :: ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3
  ${GPU_DECLARE_CREATE('ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3')}$

  !> minimum and maximum domain boundaries for each dimension
  double precision  :: xprobmin1,xprobmin2,xprobmin3,xprobmax1,xprobmax2,&
     xprobmax3
 ${GPU_DECLARE_CREATE('xprobmin1,xprobmin2,xprobmin3,xprobmax1,xprobmax2,xprobmax3')}$

  !> Offset of the logarithmic radial map, r0 in s = ln(1 + r/r0).
  !>
  !> Only meaningful in a LOG_RADIUS build. Zero - the default - selects the
  !> plain map s = ln(r), which cannot reach r = 0. A positive value moves the
  !> singularity of the map to r = -r0, so the mesh is uniform (with width
  !> r0*ds) for r << r0 and logarithmic for r >> r0, and r = 0 becomes an
  !> ordinary face of the grid at s = 0.
  double precision :: log_r0 = 0.0d0

  !> The same map written as r = log_ra*exp(s) + log_rb, which is what device
  !> code evaluates where s >= 0 is guaranteed. log_ra = 1, log_rb = 0 for
  !> log_r0 = 0, so a plain logarithmic build gets bit-for-bit exp(s).
  double precision :: log_ra = 1.0d0, log_rb = 0.0d0
  ${GPU_DECLARE_CREATE('log_r0, log_ra, log_rb')}$

  !> Indices for cylindrical coordinates FOR TESTS, negative value when not used:
  integer :: r_ = -1
  integer :: phi_ = -1
  integer :: z_ = -1
  ${GPU_DECLARE_COPYIN('r_, phi_, z_')}$

  !> Number of spatial dimensions for grid variables
  integer, parameter :: ndim=3

  !> Number of spatial dimensions (components) for vector variables
  integer :: ndir=ndim
  ${GPU_DECLARE_COPYIN('ndir')}$

  !> starting dimension for electric field 
  integer, parameter :: sdim=1

  !> Cartesian geometry or not
  logical :: slab = .true.
  ${GPU_DECLARE_COPYIN('slab')}$

  !> uniform Cartesian geometry or not (stretched Cartesian)
  logical :: slab_uniform=.true.
  ${GPU_DECLARE_COPYIN('slab_uniform')}$
  
  !> each cell has its own timestep or not
  logical :: local_timestep = .false.
  ${GPU_DECLARE_COPYIN('local_timestep')}$

  !> number of grid blocks in domain per dimension, in array over levels
  integer, dimension(:), allocatable :: ng1,ng2,ng3
  !> extent of grid blocks in domain per dimension, in array over levels
  double precision, dimension(:), allocatable :: dg1,dg2,dg3

  !> number of cells for each dimension in level-one mesh
  integer :: domain_nx1,domain_nx2,domain_nx3

  !> number of cells for each dimension in grid block excluding ghostcells
  integer :: block_nx1,block_nx2,block_nx3
  ${GPU_DECLARE_CREATE('block_nx1,block_nx2,block_nx3')}$

  !> Lower index of grid block arrays (always 1)
  integer, parameter :: ixGlo1 = 1, ixGlo2 = 1, ixGlo3 = 1

  !> Upper index of grid block arrays
  integer :: ixGhi1,ixGhi2,ixGhi3
  ${GPU_DECLARE_CREATE('ixGhi1,ixGhi2,ixGhi3')}$

  !> Lower index of stagger grid block arrays (always 0)
  integer, parameter :: ixGslo1 = 0, ixGslo2 = 0, ixGslo3 = 0

  !> Upper index of stagger grid block arrays
  integer :: ixGshi1,ixGshi2,ixGshi3
  ${GPU_DECLARE_CREATE('ixGshi1,ixGshi2,ixGshi3')}$

  !> Index range over which the cell metrics (dvolume, ds) are known.  This is
  !> the ghosted block range ixG, plus one extra layer on either side when the
  !> number of ghost cells is odd (or the grid is staggered), because ghost
  !> cell prolongation then needs volumes and dx-intervals one cell beyond it.
  integer :: ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,ixGextmax3

  !> Number of ghost cells surrounding a grid
  integer :: nghostcells = 2
  ${GPU_DECLARE_COPYIN('nghostcells')}$

  ! --- Vestigial general-grid-stretching state -------------------------------
  ! Upstream MPI-AMRVAC's qstretch/dxfirst stretched-grid machinery has been
  ! removed: this fork builds every cell metric analytically on the device
  ! from the block corner and a constant spacing, which general stretch_dim
  ! stretching cannot express. The one supported non-uniform radius,
  ! geometry='logSpherical'/'logCylindrical', is still a uniform mesh in the
  ! logical coordinate xi = ln(r + log_r0) and needs none of this.
  !
  ! The variables below are kept only because a handful of code paths this
  ! fork does not exercise (mod_slice's get_igslice, mod_particle_base's
  ! particle-index lookup, mod_thermal_emission's line-of-sight sub-grids)
  ! still branch on them. They are pinned to the unstretched state here and
  ! there is no longer any way to change them: stretch_dim and friends are no
  ! longer read from &meshlist.
  integer, parameter :: stretch_none = 0 !< No stretching
  integer, parameter :: stretch_uni  = 1 !< Unidirectional stretching (unused)
  integer, parameter :: stretch_symm = 2 !< Symmetric stretching (unused)

  logical :: stretch_uncentered = .false.
  logical :: stretched_dim(ndim) = .false.
  integer :: stretch_type(ndim) = stretch_none
  double precision :: qstretch_baselevel(ndim) = 1.0d0
  integer :: nstretchedblocks_baselevel(ndim) = 0
  integer, allocatable :: nstretchedblocks(:,:)
  double precision :: xstretch1 = 0.0d0, xstretch2 = 0.0d0, xstretch3 = 0.0d0
  double precision, allocatable :: qstretch(:,:), dxfirst(:,:), dxfirst_1mq(:,:),&
     dxmid(:,:)
  ! --------------------------------------------------------------------------

  !> grid hierarchy info (level and grid indices)
  integer, parameter :: nodehi=3+1
  integer, parameter :: plevel_=1
  integer, parameter :: pig1_=plevel_+1,pig2_=plevel_+2,pig3_=plevel_+3

  integer, allocatable :: node(:,:)
  integer, allocatable :: node_sub(:,:)
  ${GPU_DECLARE_CREATE('node')}$

  !> grid location info (corner coordinates and grid spacing)
  integer, parameter :: rnodehi=3*3
  integer, parameter :: rpxmin0_=0
  integer, parameter :: rpxmin1_=rpxmin0_+1,rpxmin2_=rpxmin0_+2,&
     rpxmin3_=rpxmin0_+3
  integer, parameter :: rpxmax0_=3
  integer, parameter :: rpxmax1_=rpxmax0_+1,rpxmax2_=rpxmax0_+2,&
     rpxmax3_=rpxmax0_+3
  integer, parameter :: rpdx1_=2*3+1,rpdx2_=2*3+2,rpdx3_=2*3+3

  !> Corner coordinates
  double precision, allocatable :: rnode(:,:)
  ${GPU_DECLARE_CREATE('rnode')}$
  double precision, allocatable :: rnode_sub(:,:)

  double precision, allocatable :: dx(:,:)
  double precision :: dxlevel(ndim)
  ${GPU_DECLARE_CREATE('dxlevel')}$

  ! IO related quantities

  !> Maximum number of saves that can be defined by tsave or itsave
  integer, parameter :: nsavehi=100

  !> Number of output methods
  integer, parameter :: nfile = 5

  !> index number of the latest existing data file
  integer :: index_latest_data

  !> Names of the output methods
  character(len=40), parameter  :: output_names(nfile) = ['log      ',&
      'normal   ', 'slice    ', 'collapsed', 'analysis ']

  !> User parameter file
  character(len=std_len)   :: usr_filename

  !> If collapse(DIM) is true, generate output integrated over DIM
  logical :: collapse(ndim)

  !> Save output of type N on times tsave(:, N)
  double precision :: tsave(nsavehi,nfile)

  double precision :: tsavelast(nfile)

  !> Repeatedly save output of type N when dtsave(N) simulation time has passed
  double precision :: dtsave(nfile)

  !> Save output of type N on iterations itsave(:, N)
  integer :: itsave(nsavehi,nfile)

  integer :: itsavelast(nfile)

  !> Repeatedly save output of type N when ditsave(N) time steps have passed
  integer :: ditsave(nfile)

  integer :: isavet(nfile)

  integer :: isaveit(nfile)

  !> Start of read out (not counting specified read outs)
  double precision :: tsavestart(nfile)

  !> The level at which to produce line-integrated / collapsed output
  integer :: collapseLevel

  !> Number of saved files of each type
  integer :: n_saves(1:nfile)

  !> whether or not to save an output file
  logical :: save_file(nfile)

  !> to monitor timeintegration loop at given wall-clock time intervals
  double precision :: time_between_print

  !> accumulated wall-clock time spent on boundary conditions
  double precision :: time_bc

  !> IO: snapshot and collapsed views output numbers/labels
  integer :: snapshotnext, collapsenext

  !> IO: slice output numbers/labels
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

  !> IO switches for conversion
  logical          :: nocartesian
  logical, allocatable :: w_write(:)
  logical, allocatable :: writelevel(:)
  double precision :: writespshift(ndim,2)
  integer          :: level_io, level_io_min, level_io_max, type_endian

  !> Which par files are used as input
  character(len=std_len), allocatable :: par_files(:)

  !> Base file name for simulation output, which will be followed by a number
  character(len=std_len) :: base_filename

  !> If not 'unavailable', resume from snapshot with this base file name
  character(len=std_len) :: restart_from_file

  !> Which type of log to write: 'normal', 'special', 'regression_test'
  character(len=std_len) :: typefilelog

  !> Resume from the snapshot with this index
  integer :: snapshotini

  !> If true, restart a previous run from the latest snapshot
  logical :: resume_previous_run

  !> If true and restart_from_file is given, convert snapshots to
  !> other file formats
  logical                :: convert

  !> If true, already convert to output format during the run
  logical                :: autoconvert

  !> If true, convert from conservative to primitive variables in output
  logical                :: saveprim

  !> Which format to use when converting
  !>
  !> Options are: tecplot, tecplotCC, vtu, vtuCC, vtuB, vtuBCC,
  !> tecplotmpi, tecplotCCmpi, vtumpi, vtuCCmpi, vtuBmpi, vtuBCCmpi, pvtumpi, pvtuCCmpi,
  !> pvtuBmpi, pvtuBCCmpi, tecline, teclinempi, onegrid
  character(len=std_len) :: convert_type

  character(len=std_len) :: collapse_type

  !> Conversion factors the primitive variables
  double precision, allocatable :: w_convert_factor(:)

  double precision :: length_convert_factor

  !> Conversion factor for time unit
  double precision       :: time_convert_factor

  !> Stores the memory and load imbalance, used in printlog
  double precision       :: Xload, Xmemory

  !> Save a snapshot before crash a run met unphysical values
  logical :: crash=.false.
  ${GPU_DECLARE_COPYIN('crash')}$
  
  !> type of physics to build
  character(len=std_len) :: phys

  !> coordinate system to build for: 'Cartesian', 'spherical', 'cylindrical',
  !> or the logarithmic-radius variants 'logSpherical'/'logCylindrical' (which
  !> emit GEOM='spherical'/'cylindrical' plus the LOG_RADIUS flag).
  !> Read from &meshlist and turned into the GEOM fypp define by the config
  !> system; must be consistent with the set_coordinate_system call in mod_usr.
  character(len=std_len) :: geometry='Cartesian'
  
  ! Physics factors

  !> Physical scaling factor for length
  double precision :: unit_length=1.d0
  ${GPU_DECLARE_COPYIN('unit_length')}$

  !> Physical scaling factor for time
  double precision :: unit_time=1.d0
  ${GPU_DECLARE_COPYIN('unit_time')}$

  !> Physical scaling factor for density
  double precision :: unit_density=1.d0
  ${GPU_DECLARE_COPYIN('unit_density')}$

  !> Physical scaling factor for velocity
  double precision :: unit_velocity=1.d0
  ${GPU_DECLARE_COPYIN('unit_velocity')}$

  !> Physical scaling factor for temperature
  double precision :: unit_temperature=1.d0
  ${GPU_DECLARE_COPYIN('unit_temperature')}$

  !> Physical scaling factor for pressure
  double precision :: unit_pressure=1.d0
  ${GPU_DECLARE_COPYIN('unit_pressure')}$

  !> Physical scaling factor for magnetic field
  double precision :: unit_magneticfield=1.d0
  ${GPU_DECLARE_COPYIN('unit_magneticfield')}$

  !> Physical scaling factor for number density
  double precision :: unit_numberdensity=1.d0
  ${GPU_DECLARE_COPYIN('unit_numberdensity')}$

  !> Physical scaling factor for charge
  double precision :: unit_charge=1.d0
  ${GPU_DECLARE_COPYIN('unit_charge')}$

  !> Physical scaling factor for mass
  double precision :: unit_mass=1.d0
  ${GPU_DECLARE_COPYIN('unit_mass')}$

  !> Normalised speed of light
  double precision :: c_norm=1.d0
  ${GPU_DECLARE_COPYIN('c_norm')}$

  !> Physical scaling factor for Opacity
  double precision :: unit_opacity=1.d0
  ${GPU_DECLARE_COPYIN('unit_opacity')}$

  !> Physical scaling factor for radiation flux
  double precision :: unit_radflux=1.d0
  ${GPU_DECLARE_COPYIN('unit_radflux')}$

  !> amplitude of background dipolar, quadrupolar, octupolar, user's field
  double precision :: Bdip=0.d0
  double precision :: Bquad=0.d0
  double precision :: Boct=0.d0
  double precision :: Busr=0.d0
  ${GPU_DECLARE_COPYIN('Bdip, Bquad, Boct, Busr')}$

  !> split magnetic field as background B0 field
  ! TODO these should be moved in a different file  
  logical :: B0field=.false.
  logical :: B0fieldAllocCoarse=.false.
  ${GPU_DECLARE_COPYIN('B0field, B0fieldAllocCoarse')}$

  ! number of equilibrium set variables, besides the mag field
  integer :: number_equi_vars = 0

  !> Use SI units (.true.) or use cgs units (.false.)
  logical :: SI_unit=.false.
  ${GPU_DECLARE_COPYIN('SI_unit')}$

  !> Use TRAC (Johnston 2019 ApJL, 873, L22) for MHD or 1D HD
  logical :: phys_trac=.false.
  integer :: phys_trac_type=1
  double precision :: phys_trac_mask

  !> Use particles module or not
  logical :: use_particles=.false.

  !> Use multigrid (only available in 2D and 3D)
  logical :: use_multigrid = .false.

  ! AMR switches

  !> The maximum number of grid blocks in a processor
  integer :: max_blocks=4000
  ${GPU_DECLARE_COPYIN('max_blocks')}$

  !> The maximum number of levels in the grid refinement
  integer, parameter :: nlevelshi = 20

  !> Maximal number of AMR levels
  integer :: refine_max_level
  ${GPU_DECLARE_CREATE('refine_max_level')}$

  !> Specify to use user-defined refinement criterion
  logical :: refine_usr = .false.
  
  !> Weights of variables used to calculate error for mesh refinement
  double precision, allocatable :: w_refine_weight(:)
  ${GPU_DECLARE_CREATE('w_refine_weight')}$

  !> Fix the AMR grid after this time
  double precision :: tfixgrid

  !> Whether to apply flux conservation at refinement boundaries
  logical :: fix_conserve_global = .true.

  !> Fix the AMR grid after this many time steps
  integer :: itfixgrid

  !> Reconstruct the AMR grid once every ditregrid iteration(s)
  integer :: ditregrid

  !> refinement: lohner estimate wavefilter setting
  double precision, allocatable :: amr_wavefilter(:)
  ${GPU_DECLARE_CREATE('amr_wavefilter')}$

  integer                       :: refine_criterion
  logical                       :: prolongprimitive=.false.
  logical                       :: coarsenprimitive=.false.
  ${GPU_DECLARE_COPYIN('prolongprimitive, coarsenprimitive')}$

  !> Error tolerance for refinement decision
  double precision, allocatable :: refine_threshold(:)
  double precision, allocatable :: derefine_ratio(:)
  ${GPU_DECLARE_CREATE('refine_threshold, derefine_ratio')}$

  !> If true, rebuild the AMR grid upon restarting
  logical :: reset_grid
  !> True for using stagger grid
  logical :: stagger_grid=.false.
  ${GPU_DECLARE_COPYIN('stagger_grid')}$
  
  !> True for record electric field
  logical :: record_electric_field=.false.

  !> Number of cells as buffer zone
  integer :: nbufferx1,nbufferx2,nbufferx3

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

  ! Time integration aspects

  double precision :: dt
  ${GPU_DECLARE_CREATE('kr,lvc,dt')}$

  logical :: time_advance
  ${GPU_DECLARE_CREATE('time_advance')}$

  !> The Courant (CFL) number used for the simulation
  double precision :: courantpar
  ${GPU_DECLARE_CREATE('courantpar')}$

  !> How to compute the CFL-limited time step
  integer :: type_courant=1
  !> integer switchers for type courant
  integer, parameter :: type_maxsum=1
  integer, parameter :: type_summax=2
  integer, parameter :: type_minimum=3

  !> If dtpar is positive, it sets the timestep dt, otherwise courantpar is used
  !> to limit the time step based on the Courant condition.
  double precision :: dtpar

  !> For resistive MHD, the time step is also limited by the diffusion time:
  !> \f$ dt < dtdiffpar \times dx^2/eta \f$
  double precision :: dtdiffpar
  ${GPU_DECLARE_CREATE('dtdiffpar')}$

  !> The global simulation time
  double precision :: global_time
  ${GPU_DECLARE_CREATE('global_time')}$

  !> Start time for the simulation
  double precision :: time_init

  !> End time for the simulation
  double precision :: time_max

  !> Ending wall time (in hours) for the simulation
  double precision :: wall_time_max

  !> Stop the simulation when the time step becomes smaller than this value
  double precision :: dtmin

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

  !> If true, wall time is up, modify snapshotnext for later overwrite
  logical :: pass_wall_time

  !> If true, do H-correction to fix the carbuncle problem at grid-aligned shocks
  logical :: H_correction=.false.
  ${GPU_DECLARE_COPYIN('H_correction')}$

  !> Number of time steps taken
  integer :: it
  ${GPU_DECLARE_CREATE('it')}$

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

  !> How many time levels of bg the time integrator needs: its number of
  !> registers, which is not its number of sub-steps (see time_levels)
  integer :: nlevel

  !> Which flux scheme of spatial discretization to use (per grid level)
  integer, allocatable :: flux_method(:)

  !> flux schemes
  integer, parameter :: fs_hll=1
  integer, parameter :: fs_hllc=2
  integer, parameter :: fs_tvdlf=5

  !> Finite-volume layout (fv_face_flux in &methodlist): fv_per_cell, the
  !> default, has every cell compute both of its faces, so each face twice,
  !> with no face buffer; fv_shared computes each face flux once and shares
  !> it between the two cells, at the cost of a per-block face buffer on the
  !> device that can exhaust its memory
  integer, parameter :: fv_shared=1
  integer, parameter :: fv_per_cell=2
  integer :: fv_layout = fv_per_cell

  !> Adaptive LLF diffusion (Rempel et al. 2009): phi in [flux_ad_min, 1] per face per variable
  logical  :: flux_adaptive_diffusion = .false.
  ${GPU_DECLARE_CREATE('flux_adaptive_diffusion')}$
  double precision :: flux_ad_min             = 0.0d0
  ${GPU_DECLARE_CREATE('flux_ad_min')}$
  double precision :: flux_ad_scale           = 1.0d0
  ${GPU_DECLARE_CREATE('flux_ad_scale')}$

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

  !> Type of slope limiter used for reconstructing variables on cell edges
  integer, allocatable :: type_limiter(:)
  ${GPU_DECLARE_CREATE('type_limiter')}$

  !> Type of slope limiter used for computing gradients or divergences, when
  !> typegrad or typediv are set to 'limited'
  integer, allocatable :: type_gradient_limiter(:)
  ${GPU_DECLARE_CREATE('type_gradient_limiter')}$

  !> background magnetic field location indicator
  integer :: b0i=0
  ${GPU_DECLARE_COPYIN('b0i')}$

  !> Limiter used for prolongation to refined grids and ghost cells
  integer :: prolong_limiter=0

  character(len=std_len) :: typedimsplit
  character(len=std_len) :: geometry_name='default'
  character(len=std_len) :: typepoly
  ${GPU_DECLARE_COPYIN('typedimsplit, geometry_name, typepoly')}$

  logical, allocatable          :: loglimit(:), logflag(:)
  ${GPU_DECLARE_CREATE('loglimit, logflag')}$
  logical                       :: flatcd,flatsh
  ${GPU_DECLARE_CREATE('flatcd, flatsh')}$
  !> Use split or unsplit way to add user's source terms, default: unsplit
  logical                       :: source_split_usr
  ${GPU_DECLARE_CREATE('source_split_usr')}$
  !> if any normal source term is added in split fasion
  logical                       :: any_source_split=.false.
  ${GPU_DECLARE_COPYIN('any_source_split')}$
  logical                       :: dimsplit
  ${GPU_DECLARE_CREATE('dimsplit')}$

  !> RK2(alfa) method parameters from Butcher tableau
  double precision              :: rk_a21,rk_b1,rk_b2
  !> IMEX-222(lambda) one-parameter family of schemes
  double precision              :: imex222_lambda
  !> SSPRK choice of methods (both threestep and fourstep, Shu-Osher 2N* implementation)
  !> also fivestep SSPRK54
  integer                       :: ssprk_order
  double precision              :: rk_beta11,rk_beta22,rk_beta33,rk_beta44,&
     rk_c2,rk_c3,rk_c4
  double precision              :: rk_alfa21,rk_alfa22,rk_alfa31,rk_alfa33,&
     rk_alfa41,rk_alfa44
  double precision              :: rk_beta54,rk_beta55,rk_alfa53,rk_alfa54,&
     rk_alfa55,rk_c5
 ${GPU_DECLARE_CREATE('rk_beta11,rk_beta22,rk_beta33,rk_beta44,rk_c2,rk_c3,rk_c4')}$
 ${GPU_DECLARE_CREATE('rk_alfa21,rk_alfa22,rk_alfa31,rk_alfa33,rk_alfa41,rk_alfa44')}$
 ${GPU_DECLARE_CREATE('rk_beta54,rk_beta55,rk_alfa53,rk_alfa54,rk_alfa55,rk_c5')}$
  !> RK3 Butcher table
  integer                       :: rk3_switch
  double precision              :: rk3_a21,rk3_a31,rk3_a32,rk3_b1,rk3_b2,&
     rk3_b3,rk3_c2,rk3_c3
  !> IMEX_ARS3 parameter ars_gamma
  double precision              :: ars_gamma
  !> IMEX_232 choice and parameters
  integer                       :: imex_switch
  double precision              :: imex_a21,imex_a31,imex_a32,imex_b1,imex_b2,&
     imex_ha21,imex_ha22
  double precision              :: imex_b3,imex_c2,imex_c3
  !> IMEX_CB3a extra parameters
  double precision              :: imex_a22, imex_a33, imex_ha32
  !> whether IMEX in use or not
  logical                       :: use_imex_scheme
  ${GPU_DECLARE_CREATE('use_imex_scheme')}$

  character(len=std_len) :: typediv,typegrad

  !> global fastest wave speed needed in fd scheme and glm method
  double precision :: cmax_global
  ${GPU_DECLARE_CREATE('cmax_global')}$

  !> global largest a2 for schmid scheme
  double precision :: a2max_global(ndim)
  ${GPU_DECLARE_CREATE('a2max_global')}$

  !> need global maximal wave speed
  logical :: need_global_cmax=.false.
  ${GPU_DECLARE_CREATE('need_global_cmax')}$

  !> global value for schmid scheme
  logical :: need_global_a2max=.false.
  ${GPU_DECLARE_CREATE('need_global_a2max')}$

  ! Boundary region parameters

  !> True for dimensions with periodic boundaries
  logical :: periodB(ndim)
  ${GPU_DECLARE_CREATE('periodB')}$

  !> Indicates whether there is a pole at a boundary
  logical :: poleB(2,ndim)
  ${GPU_DECLARE_CREATE('poleB')}$

  !> True where the parameter file asked for typeboundary='pole' on a face.
  !> Kept separately from poleB because the two are decided independently and
  !> at different times - read_par_files expands the 'pole' string into
  !> symm/asymm entries and loses it, while set_pole derives poleB from the
  !> domain bounds much later - and check_pole_setup has to be able to tell
  !> the user when they disagree.
  logical :: poleB_requested(2,ndim) = .false.

  !> True for dimensions with aperiodic boundaries
  logical :: aperiodB(ndim)
  ${GPU_DECLARE_CREATE('aperiodB')}$

  !> True for save physical boundary cells in dat files
  logical :: save_physical_boundary
  ${GPU_DECLARE_CREATE('save_physical_boundary')}$

  !> True if a block has any physical boundary
  logical, allocatable :: phyboundblock(:)
  ${GPU_DECLARE_CREATE('phyboundblock')}$

  !> Array indicating the type of boundary condition per variable and per
  !> physical boundary
  integer, allocatable :: typeboundary(:, :)
  ${GPU_DECLARE_CREATE('typeboundary')}$
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
  !> signal to the precompiler that we need special boundaries:
  logical            :: specialboundary = .false.
  ${GPU_DECLARE_COPYIN('specialboundary')}$

  !> whether copy values instead of interpolation in ghost cells of finer blocks
  logical :: ghost_copy=.false.
  ${GPU_DECLARE_CREATE('ghost_copy')}$

  !> if there is an internal boundary
  logical :: internalboundary
  ${GPU_DECLARE_CREATE('internalboundary')}$

  !> Base file name for synthetic EUV emission output
  character(len=std_len) :: filename_euv
  !> wavelength for output
  integer :: wavelength
  !> times for enhancing spatial resolution for EUV image/spectra
  double precision :: instrument_resolution_factor
  !> use arcsec as length unit of images/spectra
  logical :: activate_unit_arcsec
  !> Base file name for synthetic SXR emission output
  character(len=std_len) :: filename_sxr
  ! minimum and maximum energy of SXR (keV)
  integer :: emin_sxr,emax_sxr
  !> Base file name for synthetic white light
  character(len=std_len) :: filename_whitelight
  !> white light observation instrument
  character(len=std_len) :: whitelight_instrument
  !> the white light emission below it (unit=Rsun) is not visible 
  double precision :: R_occultor
  !> direction of the line of sight (LOS)
  double precision :: LOS_theta,LOS_phi
  !> rotation of image
  double precision :: image_rotate
  !> where the is the origin (X=0,Y=0) of image
  double precision :: x_origin(1:3)
  !> big image
  logical :: big_image
  !> Base file name for synthetic EUV spectrum output
  character(len=std_len) :: filename_spectrum
  !> wave length for spectrum
  integer :: spectrum_wl
  !> spectral window
  double precision :: spectrum_window_min,spectrum_window_max
  !> location of the slit
  double precision :: location_slit
  !> direction of the slit (for dat resolution only)
  integer :: direction_slit
  !> for spherical coordinate, region below it (unit=Rsun) is treated as not transparent
  double precision :: R_opt_thick
  !> resolution of the images
  logical :: dat_resolution

  !> Block pointer for using one block and its previous state
  type(state), pointer :: block
  ${GPU_DECLARE_CREATE('block')}$
contains

  !> Sub-steps and time levels of a time integrator, the single place either
  !> is defined.  A stage reads one register and updates another in place, and
  !> the two must differ because the flux stencil reads neighbours, so the
  !> level count is the number of distinct registers the scheme's sequence in
  !> advect touches - fewer than its stages for the low-storage SSP schemes,
  !> more for RK3_BT and the IMEX family.  bg is allocated with nlevels.
  pure subroutine time_levels(integrator, nstages, nlevels)
    integer, intent(in)  :: integrator
    integer, intent(out) :: nstages, nlevels

    select case (integrator)
    case (Forward_Euler, IMEX_Euler, IMEX_SP)
       nstages = 1; nlevels = 2
    case (Predictor_Corrector, RK2_alf, ssprk2)
       nstages = 2; nlevels = 2
    case (IMEX_Midpoint, IMEX_Trapezoidal, IMEX_222)
       nstages = 2; nlevels = 3
    case (ssprk3)
       nstages = 3; nlevels = 3
    case (RK3_BT, IMEX_CB3a)
       nstages = 3; nlevels = 4
    case (IMEX_ARS3, IMEX_232)
       nstages = 3; nlevels = 5
    case (ssprk4)
       nstages = 4; nlevels = 3
    case (rk4)
       nstages = 4; nlevels = 4
    case (ssprk5)
       nstages = 5; nlevels = 4
    case default
       nstages = 0; nlevels = 0
    end select
  end subroutine time_levels

  !> Cross product of two vectors
  pure subroutine cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,axb)
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,3), b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       3)
    double precision, intent(out) :: axb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,3)

    axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)=a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)-a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
    axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)=a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)-a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
    axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)=a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)-a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
  end subroutine cross_product

end module mod_global_parameters
