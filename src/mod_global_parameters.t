!> This module contains definitions of global parameters and variables
!> \todo Move the parameters to the relevant (physics) modules
module mod_global_parameters
  use mod_physicaldata
  use mod_connectivity
  use mpi
  use mod_config

  implicit none
  public

  ! Parameters

  !> Object with configuration settings
  type(CFG_t) :: cfg

  !> The number of interleaving sending buffers for ghost cells
  integer, parameter :: npwbuf=2

  integer :: ixM^LL

  integer, dimension(:), allocatable :: ng^D
  double precision, dimension(:), allocatable :: dg^D

  logical :: slab

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
  integer :: type_subblock_io, type_subblock_x_io
  integer :: type_block_xc_io,type_block_xcc_io
  integer :: type_block_wc_io,type_block_wcc_io
  integer :: itag
  integer :: irecv, isend
  integer, dimension(:), allocatable :: recvrequest, sendrequest
  integer, dimension(:,:), allocatable :: recvstatus, sendstatus

  integer :: snapshot, snapshotnext, slice, slicenext, collapseNext, icollapse

  logical, allocatable, dimension(:^D&) :: patchfalse

  !> split potential or linear force-free magnetic field as background B0 field
  logical :: B0field
  !> amplitude of background dipolar, quadrupolar, octupolar, user's field
  double precision :: Bdip, Bquad, Boct, Busr

  !> Default length for strings
  integer, parameter :: std_len = 131

  !> Indices for cylindrical coordinates FOR TESTS, negative value when not used:
  integer, parameter :: r_=1, phi_=^PHI, z_=^Z

  !> Indices for cylindrical coordinates FOR INDEXING, always positive
  !> \todo Check whether these are still needed
  integer, parameter :: pphi_=^PPHI, zz_=^ZZ

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

  !> Maximum number of slices
  integer, parameter :: nslicemax=1000

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

  !> A very large integer
  integer, parameter :: biginteger=10000000

  ! A very small real number (but above machine precision)
  double precision, parameter :: smalldouble=1.D-12

  !> A very large real number
  double precision, parameter :: bigdouble=1.D+99

  !> \todo Remove these
  double precision, parameter :: zero    = 0.0d0
  double precision, parameter :: one     = 1.0d0
  double precision, parameter :: two     = 2.0d0
  double precision, parameter :: half    = 0.5d0
  double precision, parameter :: quarter = 0.25d0
  double precision, parameter :: third   = 1/3.0d0

  !> Pi
  double precision, parameter :: dpi=3.141592653589793238462643383279502884197169399375105d0

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

  !> Elapsed time for evaluate the performance
  real :: time_elapsed=0.0

  ! TODO temporary placed here (Jannis)
  integer           :: nwflux
  integer           :: nwaux
  integer           :: nwextra
  integer           :: nw
  integer           :: nvector

  integer                       :: nflag_
  integer, allocatable          :: flags(:)
  double precision, allocatable :: wflags(:)

  integer, dimension(:), allocatable :: iw_vector

  integer           :: iprob
  !> positions of the minimum and maximum surfaces for each dimension
  double precision  :: xprob^L

  ! TODO: remove
  integer, parameter :: neqpar = 1
  integer, parameter :: nspecialpar = 0
  double precision :: eqpar(neqpar) = -0.0d0
  character(len=*), parameter :: eqparname = "DEPRECATED"
  character(len=*), parameter :: specialparname = "DEPRECATED"

  !> For transform variables and save selected data
  !> number of w in the transformed data
  integer :: nwtf
  !> number of equation parameters in the transformed data
  integer :: neqpartf

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

  {#IFDEF MAGNETOFRICTION
  double precision :: cmf_c,cmf_y,cmf_divb
  }

  !> How to compute the residual
  character(len=std_len) :: typeresid

  !> Under construction
  !> \todo Remove time_accurate?
  logical :: time_accurate

  !> Enable additional MPI_BARRIER calls, useful when debugging on new platforms
  !> \todo Remove addmpibarrier?
  logical :: addmpibarrier

  ! Time parameters

  !> Maximum number of saves that can be defined by tsave or itsave
  integer, parameter :: nsavehi=100

  !> The global simulation time
  !> \todo change t to a longer name
  double precision :: t

  !> End time for the simulation
  double precision :: tmax

  !> Stop the simulation when the time step becomes smaller than this value
  double precision :: dtmin

  !> \todo Remove residmin?
  double precision :: residmin

  !> \todo Remove residmax?
  double precision :: residmax

  !> \todo Remove residual?
  double precision :: residual

  {#IFDEF MAGNETOFRICTION
  !> \todo What is tmf?
  double precision :: tmf
  }

  !> Save output of type N on times tsave(:, N)
  double precision :: tsave(nsavehi,nfile)

  !> \todo Move tsavelast to amrvac.t
  double precision :: tsavelast(nfile)

  !> Repeatedly save output of type N when dtsave(N) simulation time has passed
  double precision :: dtsave(nfile)

  !> Slice coordinates, see @ref slices.md
  double precision :: slicecoord(nslicemax)

  !> If true, the last time step will be reduced so that the final time is the
  !> end time of the simulation
  logical :: tmaxexact

  !> If true, do not use the tmax stored in a snapshot when restarting
  logical :: treset

  !> If true, do not use the itmax stored in a snapshot when restarting
  logical :: itreset

  !> If true, call initonegrid_usr upon restarting
  logical :: firstprocess

  !> If true, rebuild the AMR grid upon restarting
  logical :: resetgrid

  !> Call the process subroutine before the writing of a snapshot, and following
  !> the determination of the timestep constraint by means of CFL and other
  !> restrictions
  logical :: fixprocess

  !> If collapse(DIM) is true, generate output integrated over DIM
  logical :: collapse(ndim)

  !> Number of time steps taken
  integer :: it

  !> Stop the simulation after this many time steps have been taken
  integer :: itmax

  !> \todo Why do we need itmin?
  integer :: itmin

  !> If > 1, then in the first slowsteps-1 time steps dt is reduced
  !> by a factor \f$ 1 - (1- step/slowsteps)^2 \f$
  integer :: slowsteps

  {#IFDEF MAGNETOFRICTION
  integer :: itmaxmf, ditsavemf
  }

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

  !> Number of slices to output
  integer :: nslices

  !> The slice direction for each slice
  integer :: slicedir(nslicemax)

  !> The level at which to produce line-integrated / collapsed output
  integer :: collapseLevel

  !> Number of saved files of each type
  !> \todo Move to mod_input_output
  integer :: n_saves(1:nfile)

  !> Options are 1: Parallel MPI output, 0: master-slave parallel IO, -1:
  !> master-slave IO without MPI (no MPI_FILE_WRITE, MPI_FILE_OPEN etc)
  integer :: typeparIO

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
  character(len=std_len) :: typeadvance

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
  character(len=std_len), allocatable :: typefull1(:)

  !> The spatial dicretization to use for the predictor step when using a two
  !> step method
  character(len=std_len), allocatable :: typepred1(:)

  !> Type of slope limiter used for reconstructing variables on cell edges
  character(len=std_len), allocatable :: typelimiter1(:)

  !> Type of slope limiter used for computing gradients or divergences, when
  !> typegrad or typediv are set to 'limited'
  character(len=std_len), allocatable :: typegradlimiter1(:)

  !> \todo Remove / replace with typelimiter1
  character(len=std_len) :: typelimiter

  !> \todo Remove / replace with typegradlimiter1
  character(len=std_len) :: typegradlimiter

  !> Limiter used for prolongation to refined grids and ghost cells
  character(len=std_len) :: typeprolonglimit

  !> Which type of entropy fix to use with Riemann-type solvers
  character(len=std_len), allocatable :: typeentropy(:)

  !> Which type of TVD method to use
  character(len=std_len) :: typetvd

  !> Which type of TVDLF method to use
  character(len=std_len) :: typetvdlf

  character(len=std_len) :: typeaverage
  character(len=std_len) :: typedimsplit
  character(len=std_len) :: typeaxial
  character(len=std_len) :: typecoord
  character(len=std_len) :: typepoly

  integer                       :: errorestimate,nxdiffusehllc,typespherical,ncyclemax
  double precision, allocatable :: entropycoef(:)
  double precision              :: tvdlfeps, mcbeta
  logical, allocatable          :: loglimit(:), logflag(:)
  logical                       :: flathllc,flatcd,flatsh,flatppm
  logical                       :: ssplitdivb,ssplitresis,ssplituser,useprimitive,dimsplit
  logical                       :: restrictprimitive,prolongprimitive
  logical                       :: coarsenprimitive,useprimitiveRel, amrentropy

  logical                :: divbwave,compactres,BnormLF
  double precision       :: divbdiff,smallT,smallp,smallrho
  double precision, allocatable :: amr_wavefilter(:)
  character(len=std_len) :: typedivbdiff,typedivbfix,typediv,typegrad

  !> related to primitive-conservative switch in relativistic modules
  logical          :: fixsmall,strictnr,strictsmall,strictzero,strictgetaux
  double precision :: dmaxvel,tolernr,absaccnr
  integer          :: maxitnr,nflatgetaux

  logical          :: nocartesian
  logical, allocatable :: writew(:)
  logical, allocatable :: writelevel(:)
  double precision :: writespshift(ndim,2)
  integer          :: level_io, level_io_min, level_io_max

  ! cooling related parameters
  integer                :: ncool, cmulti
  character(len=std_len) :: coolcurve,coolmethod
  double precision       :: cfrac,tlow
  logical                :: Tfix

  ! local and global fastest wave speed (computed in setdt):
  double precision :: cmax_mype, cmax_global

  ! Gravity related parameters
  double precision ::  x1ptms,x2ptms,x3ptms,ptmass

  ! Boundary region parameters

  !> Number of boundaries for grid blocks
  integer, parameter :: nhiB = 2*ndim

  !> True for dimensions with periodic boundaries
  logical :: periodB(ndim)

  !> Indicates whether there is a pole at a boundary
  logical :: poleB(2,ndim)

  !> True for dimensions with aperiodic boundaries
  logical :: aperiodB(ndim)

  !> Array indicating the type of boundary condition per variable and per
  !> physical boundary
  character(len=std_len), allocatable :: typeB(:, :)

  character(len=std_len) :: typeghostfill='linear',typegridfill
  double precision ::ratebdflux
  logical :: internalboundary

  !> Which par files are used as input
  character(len=std_len), allocatable :: par_files(:)

  !> Name of input file
  !> \todo Remove this
  character(len=std_len) :: inifile

  !> Base file name for simulation output, which will be followed by a number
  character(len=std_len) :: filenameout

  !> If not 'unavailable', resume from snapshot with this base file name
  character(len=std_len) :: filenameini

  !> Log file name (without the .log extension)
  character(len=std_len) :: filenamelog

  !> The header line in the output file
  character(len=std_len) :: fileheadout

  !> Names of the conservative variables
  character(len=1024) :: wnames

  !> Names of the primitive variables
  character(len=1024) :: primnames

  !> This variable was used to store names for the log file
  !> \todo remove this variable when printlog_special has been updated in all projects
  character(len=1024) :: wnameslog

  !> Which type of log to write: 'normal', 'special', 'regression_test'
  character(len=std_len) :: typefilelog

  !> Resume from the snapshot with this index
  integer :: snapshotini

  !> If true, enable ASCII output of slices
  logical :: sliceascii

  !> If true and filenameini and snapshotini are given, convert snapshots to
  !> other file formats
  logical                :: convert

  !> If true, already convert to output format during the run
  logical                :: autoconvert

  !> If true, convert from conservative to primitive variables in output
  logical                :: saveprim

  !> If true and doing a 1D run, use a limiter to determine corner values
  logical                :: uselimiter

  logical                :: endian_swap

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

  !> Conversion factors for length (normvar(0)) and the primitive variables
  !> (normvar(1:nw))
  double precision, allocatable :: normvar(:)

  !> Conversion factor for time unit
  double precision       :: normt

  integer                :: saveigrid

  ! Stores the memory and load imbalance, to be used in printlog:
  double precision       :: Xload, Xmemory

  double precision :: time_bc
  {#IFDEF STRETCHGRID
  ! stretching factor qst for log stretch grid
  double precision :: logG, qst
  double precision, allocatable :: logGs(:), qsts(:)
  }

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
  logical :: specialtol
  double precision, allocatable :: tol(:)
  double precision, allocatable :: tolratio(:)
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
  integer :: ngridshi

  !> The maximum number of levels in the grid refinement
  integer :: nlevelshi

  !> Maximal number of AMR levels
  integer :: mxnest

  !> number of cells for each dimension in level-one mesh
  integer :: nxlone^D

  !> number of cells for each dimension in grid block excluding ghostcells
  integer :: nxblock^D

  !> Lower index of grid block arrays (always 1)
  integer, parameter :: {ixGlo^D = 1|, }

  !> Upper index of grid block arrays
  integer :: ixGhi^D

  !> Number of ghost cells surrounding a grid
  integer :: dixB
  integer :: levmin
  integer :: levmax
  integer :: levmax_sub

  logical :: skipfinestep
  logical, allocatable :: phyboundblock(:)
  logical :: time_advance{#IFDEF MAGNETOFRICTION , mf_advance}

  !$OMP THREADPRIVATE(dxlevel{#IFDEF STRETCHGRID ,logG,qst})
  !$OMP THREADPRIVATE(saveigrid)
  !$OMP THREADPRIVATE(typelimiter,typegradlimiter)
end module mod_global_parameters
