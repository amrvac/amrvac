!> This module contains definitions of global parameters and variables
!> \todo Move the parameters to the relevant (physics) modules
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

  integer :: ixM^LL

  integer, dimension(:), allocatable :: ng^D
  double precision, dimension(:), allocatable :: dg^D

  logical :: slab

  !> @todo Move mpi related variables to e.g. mod_comm

  !> The number of MPI tasks
  integer :: npe

  !> The rank of the current MPI task
  integer :: mype

  !> The MPI communicator
  integer :: icomm

  !> A global MPI error return code
  !> @todo Make local
  integer :: ierrmpi

  integer :: log_fh
  !> MPI IO type for block including ghost cells
  integer :: type_block, type_coarse_block, type_sub_block(2^D&)
  !> MPI IO type for block excluding ghost cells
  integer :: type_block_io, size_block_io, size_block
  !> MPI IO type for transformed data excluding ghost cells
  integer :: type_block_io_tf, size_block_io_tf
  integer :: type_block_xc_io,type_block_xcc_io
  integer :: type_block_wc_io,type_block_wcc_io
  integer :: itag
  integer :: irecv, isend
  integer, dimension(:), allocatable :: recvrequest, sendrequest
  integer, dimension(:,:), allocatable :: recvstatus, sendstatus

  integer :: snapshotnext, collapseNext, icollapse

  !> split potential or linear force-free magnetic field as background B0 field
  logical :: B0field=.false.
  !> amplitude of background dipolar, quadrupolar, octupolar, user's field
  double precision :: Bdip=0.d0
  double precision :: Bquad=0.d0
  double precision :: Boct=0.d0
  double precision :: Busr=0.d0

  !> Indices for cylindrical coordinates FOR TESTS, negative value when not used:
  integer :: r_ = -1

  integer :: phi_ = -1

  integer :: z_ = -1

  !> Indices for cylindrical coordinates FOR INDEXING, always positive
  !> \todo Check whether these are still needed
  ! integer, parameter :: pphi_=^PPHI, zz_=^ZZ

  !> Number of spatial dimensions for grid variables
  integer, parameter :: ndim=^ND

  !> Number of spatial dimensions for vector variables
  integer :: ndir=ndim

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

  !> Number of output methods
  integer, parameter :: nfile         = 5

  !> Names of the output methods
  character(len=40), parameter  :: output_names(nfile) = &
       ['log      ', 'normal   ', 'slice    ', 'collapsed', 'analysis ']

  !> Unit for standard input
  integer, parameter :: unitstdin=5

  !> Unit for standard output
  integer, parameter :: unitterm=6

  !> Unit for error messages
  integer, parameter :: uniterr=6

  !> \todo Move to mod_input_output
  integer, parameter :: unitpar=9
  integer, parameter :: unitconvert=10
  integer, parameter :: unitslice=11
  integer, parameter :: unitsnapshot=12
  integer, parameter :: unitcollapse=13
  integer, parameter :: unitanalysis=14

  !> Physical scaling factor for length
  double precision :: unit_length=1.d0

  !> Physical scaling factor for time
  double precision :: unit_time=1.d0

  !> Physical scaling factor for density
  double precision :: unit_density=1.d0

  !> Physical scaling factor for velocity
  double precision :: unit_velocity=0.d0

  !> Physical scaling factor for temperature
  double precision :: unit_temperature=1.d0

  !> Physical scaling factor for pressure
  double precision :: unit_pressure=1.d0

  !> Physical scaling factor for magnetic field
  double precision :: unit_magneticfield=1.d0

  !> Physical scaling factor for number density
  double precision :: unit_numberdensity=1.d0

  !> Elapsed time for evaluate the performance
  real :: time_elapsed=0.0

  !> Weights of variables used to calculate error for mesh refinement
  double precision, allocatable :: w_refine_weight(:)

  integer           :: iprob
  !> positions of the minimum and maximum surfaces for each dimension
  double precision  :: xprob^L

  !> Kronecker delta tensor
  integer :: kr(3,3)

  !> Levi-Civita tensor
  integer :: lvc(3,3,3)

  !> The Courant (CFL) number used for the simulation
  double precision :: courantpar

  !> How to compute the CFL-limited time step.
  !>
  !> Options are 'maxsum': max(sum(c/dx)); 'summax': sum(max(c/dx)) and
  !> 'minimum: max(c/dx), where the summations loop over the grid dimensions and
  !> c is the velocity. The default 'maxsum' is the conventiontal way of
  !> computing CFL-limited time steps.
  character(len=std_len) :: typecourant

  !> If dtpar is positive, it sets the timestep dt, otherwise courantpar is used
  !> to limit the time step based on the Courant condition.
  double precision :: dtpar

  !> For resistive MHD, the time step is also limited by the diffusion time:
  !> \f$ dt < dtdiffpar \times dx^2/eta \f$
  double precision :: dtdiffpar

  !> Under construction
  !> \todo Remove time_accurate?
  logical :: time_accurate

  ! Time parameters

  !> Maximum number of saves that can be defined by tsave or itsave
  integer, parameter :: nsavehi=100

  !> The global simulation time
  double precision :: global_time

  !> Start time for the simulation
  double precision :: time_init

  !> End time for the simulation
  double precision :: time_max

  !> Stop the simulation when the time step becomes smaller than this value
  double precision :: dtmin

  !> Save output of type N on times tsave(:, N)
  double precision :: tsave(nsavehi,nfile)

  !> \todo Move tsavelast to amrvac.t
  double precision :: tsavelast(nfile)

  !> Repeatedly save output of type N when dtsave(N) simulation time has passed
  double precision :: dtsave(nfile)

  double precision :: time_between_print

  !> If true, reset iteration count and global_time to original values, and
  !> start writing snapshots at index 0
  logical :: reset_time

  !> If true, reset iteration count to 0
  logical :: reset_it

  !> If true, call initonegrid_usr upon restarting
  logical :: firstprocess

  !> If true, rebuild the AMR grid upon restarting
  logical :: reset_grid

  !> If collapse(DIM) is true, generate output integrated over DIM
  logical :: collapse(ndim)

  !> Number of time steps taken
  integer :: it

  !> Stop the simulation after this many time steps have been taken
  integer :: it_max

  !> initial iteration count
  integer :: it_init

  !> If > 1, then in the first slowsteps-1 time steps dt is reduced
  !> by a factor \f$ 1 - (1- step/slowsteps)^2 \f$
  integer :: slowsteps

  !> Save output of type N on iterations itsave(:, N)
  integer :: itsave(nsavehi,nfile)

  !> \todo remove itsavelast?
  integer :: itsavelast(nfile)

  !> Repeatedly save output of type N when ditsave(N) time steps have passed
  integer :: ditsave(nfile)

  !> \todo Move to amrvac.t
  integer :: isavet(nfile)

  !> \todo Move to amrvac.t
  integer :: isaveit(nfile)

  !> The level at which to produce line-integrated / collapsed output
  integer :: collapseLevel

  !> Number of saved files of each type
  !> \todo Move to mod_input_output
  integer :: n_saves(1:nfile)

  !> Fix the AMR grid after this time
  double precision :: tfixgrid

  !> Fix the AMR grid after this many time steps
  integer :: itfixgrid

  !> Reconstruct the AMR grid once every ditregrid iteration(s)
  integer :: ditregrid

  !> Number of auxiliary variables that are only included in the output
  integer :: nwauxio

  !> Index of the sub-step in a multi-step time integrator
  integer :: istep

  !> How many sub-steps the time integrator takes
  integer :: nstep

  ! Method switches

  !> Which time integrator to use
  character(len=std_len) :: time_integrator

  !> What should be used as a basis for the limiting in TVD methods. Options are
  !> 'original', 'previous' and 'predictor'.
  !>
  !> By default, the original value is used in 1D and for dimensional splitting,
  !> while for dimensionally unsplit multidimensional case (dimsplit=F), TVDLF
  !> and TVD-MUSCL uses the previous value from wold for limiting.
  character(len=std_len) :: typelimited

  !> How to apply dimensional splitting to the source terms, see
  !> @ref disretization.md
  character(len=std_len) :: typesourcesplit

  !> Which spatial discretization to use (per grid level)
  character(len=std_len), allocatable :: flux_scheme(:)

  !> The spatial dicretization to use for the predictor step when using a two
  !> step method
  character(len=std_len), allocatable :: typepred1(:)

  !> Type of slope limiter used for reconstructing variables on cell edges
  integer, allocatable :: type_limiter(:)

  !> Type of slope limiter used for computing gradients or divergences, when
  !> typegrad or typediv are set to 'limited'
  integer, allocatable :: type_gradient_limiter(:)

  !> \todo Remove / replace with limiter
  integer :: typelimiter

  !> \todo Remove / replace with gradient_limiter
  integer :: typegradlimiter

  !> Limiter used for prolongation to refined grids and ghost cells
  character(len=std_len) :: typeprolonglimit

  !> Which type of entropy fix to use with Riemann-type solvers
  character(len=std_len), allocatable :: typeentropy(:)

  !> Which type of TVD method to use
  character(len=std_len) :: typetvd

  !> Which type of TVDLF method to use
  character(len=std_len) :: typeboundspeed

  character(len=std_len) :: typeaverage
  character(len=std_len) :: typedimsplit
  character(len=std_len) :: typeaxial='default'
  character(len=std_len) :: typepoly

  integer                       :: refine_criterion,nxdiffusehllc,typespherical
  double precision, allocatable :: entropycoef(:)
  double precision              :: tvdlfeps
  logical, allocatable          :: loglimit(:), logflag(:)
  logical                       :: flathllc,flatcd,flatsh,flatppm
  !> Use split or unsplit way to add user's source terms, default: unsplit
  logical                       :: source_split_usr
  logical                       :: dimsplit
  logical                       :: prolongprimitive
  logical                       :: coarsenprimitive

  double precision       :: small_temperature,small_pressure,small_density
  double precision, allocatable :: amr_wavefilter(:)
  character(len=std_len) :: typediv,typegrad

  logical          :: nocartesian
  logical, allocatable :: w_write(:)
  logical, allocatable :: writelevel(:)
  double precision :: writespshift(ndim,2)
  integer          :: level_io, level_io_min, level_io_max

  ! Boundary region parameters

  !> True for dimensions with periodic boundaries
  logical :: periodB(ndim)

  !> Indicates whether there is a pole at a boundary
  logical :: poleB(2,ndim)

  !> True for dimensions with aperiodic boundaries
  logical :: aperiodB(ndim)

  !> Array indicating the type of boundary condition per variable and per
  !> physical boundary
  character(len=std_len), allocatable :: typeboundary(:, :)

  character(len=std_len) :: typeghostfill='linear',prolongation_method
  logical :: internalboundary

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

  !> If true and restart_from_file is given, convert snapshots to
  !> other file formats
  logical                :: convert

  !> If true, already convert to output format during the run
  logical                :: autoconvert

  !> If true, convert from conservative to primitive variables in output
  logical                :: saveprim

  !> Which format to use when converting
  !>
  !> Options are: idl, tecplot, tecplotCC, vtu, vtuCC, vtuB, vtuBCC, dx,
  !> tecplotmpi, tecplotCCmpi, vtumpi, vtuCCmpi, vtuBmpi, vtuBCCmpi, pvtumpi, pvtuCCmpi, 
  !> pvtuBmpi, pvtuBCCmpi, tecline, teclinempi, onegrid
  character(len=std_len) :: convert_type

  !> Data Explorer file endianness ('msb' or 'lsb' for last or most significant
  !> bit order)
  character(len=std_len) :: dxfiletype

  character(len=std_len) :: collapse_type

  !> Conversion factors the primitive variables
  double precision, allocatable :: w_convert_factor(:)

  double precision :: length_convert_factor

  !> Conversion factor for time unit
  double precision       :: time_convert_factor

  integer                :: saveigrid

  ! Stores the memory and load imbalance, to be used in printlog:
  double precision       :: Xload, Xmemory

  double precision :: time_bc

  !> Stretching factor for log stretch grid, cell size = logG * cell-center position
  double precision :: logG
  !> Stretching factor for log stretch grid, r_i = qst * r_i-1
  double precision :: qst
  !> Store stretching factors for each AMR level
  double precision, allocatable :: logGs(:), qsts(:)
  !> Switch to use stretched grid
  logical :: stretched_grid=.false.

  integer, parameter :: nodehi=^ND+1
  integer, parameter :: plevel_=1
  integer, parameter :: pig^D_=plevel_+^D

  integer, parameter :: rnodehi=3*^ND
  integer, parameter :: rpxmin0_=0
  integer, parameter :: rpxmin^D_=rpxmin0_+^D
  integer, parameter :: rpxmax0_=^ND
  integer, parameter :: rpxmax^D_=rpxmax0_+^D
  integer, parameter :: rpdx^D_=2*^ND+^D

  ! parameters for bc_phys
  integer, parameter :: ismin^D=-1+2*^D
  integer, parameter :: ismax^D=2*^D

  !> Corner coordinates
  double precision, allocatable :: rnode(:,:)
  double precision, allocatable :: rnode_sub(:,:)

  !> Error tolerance for refinement decision
  double precision, allocatable :: refine_threshold(:)
  double precision, allocatable :: derefine_ratio(:)
  double precision, allocatable :: dx(:,:)
  double precision :: dt
  double precision, allocatable :: dt_grid(:)
  double precision :: dxlevel(ndim)

  integer, allocatable :: node(:,:)
  integer, allocatable :: node_sub(:,:)

  type fluxalloc
     double precision, dimension(:^D&,:), pointer:: flux => null()
  end type fluxalloc
  !> store flux to fix conservation 
  type(fluxalloc), dimension(:,:,:), allocatable :: pflux

  !> Number of cells as buffer zone
  !> \todo is it necessary? 
  integer :: nbufferx^D
 
  !> The maximum number of grid blocks in a processor
  integer :: max_blocks

  !> The maximum number of levels in the grid refinement
  integer, parameter :: nlevelshi = 20

  !> Maximal number of AMR levels
  integer :: refine_max_level

  !> number of cells for each dimension in level-one mesh
  integer :: domain_nx^D

  !> number of cells for each dimension in grid block excluding ghostcells
  integer :: block_nx^D

  !> Lower index of grid block arrays (always 1)
  integer, parameter :: {ixGlo^D = 1|, }

  !> Upper index of grid block arrays
  integer :: ixGhi^D

  !> Number of ghost cells surrounding a grid
  integer :: nghostcells
  integer :: levmin
  integer :: levmax
  integer :: levmax_sub

  logical :: skipfinestep
  logical, allocatable :: phyboundblock(:)
  logical :: time_advance

  !> Use SI units (.true.) or use cgs units (.false.)
  logical              :: SI_unit=.false.

  !> Solve energy equation or not
  logical :: phys_energy=.true.

  !> Solve polytropic process instead of solving total energy
  logical :: solve_internal_e=.false.

  !> Use particles module or not
  logical :: use_particles=.false.

  !> Save a snapshot before crash a run met unphysical values
  logical :: crash=.false.

  ! global fastest wave speed needed in fd scheme and glm method
  double precision :: cmax_global

  !> need global maximal wave speed
  logical :: need_global_cmax=.false.

  !$OMP THREADPRIVATE(dxlevel,logG,qst)
  !$OMP THREADPRIVATE(saveigrid)
  !$OMP THREADPRIVATE(typelimiter,typegradlimiter)
end module mod_global_parameters
