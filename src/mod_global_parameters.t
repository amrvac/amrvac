!> This module contains definitions of global parameters and variables
!> \todo Move the parameters to the relevant (physics) modules
module mod_global_parameters
  use mod_indices
  use mod_physicaldata
  use mod_connectivity
  use mpi

  implicit none
  public

  ! Parameters

  !> Default length for strings
  integer, parameter :: std_len = 131

  !> Indices for cylindrical coordinates FOR TESTS, negative value when not used:
  integer, parameter :: r_=1, phi_=^PHI, z_=^Z

  !> Indices for cylindrical coordinates FOR INDEXING, always positive
  integer, parameter :: pphi_=^PPHI, zz_=^ZZ

  include 'amrvacpar.f'

  integer, parameter :: ndim=^ND, ndir=^NC

  include 'amrvacsettings.f'

  integer, parameter :: filelog_      = 1
  integer, parameter :: fileout_      = 2
  integer, parameter :: fileslice_    = 3
  integer, parameter :: filecollapse_ = 4
  integer, parameter :: fileanalysis_ = 5
  integer, parameter :: nfile         = 5

  !> Names of the various output methods
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

  ! Units reserved for files:
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

  !> Physical scaling factor for lengths
  double precision :: UNIT_LENGTH

  !> Physical scaling factor for densities
  double precision :: UNIT_DENSITY

  !> Physical scaling factor for velocities
  double precision :: UNIT_VELOCITY

  include 'amrvacusrpar.f'

  ! For transform variables and save selected data
  integer :: nwtf
  integer :: neqpartf

  !> Kronecker delta tensor
  integer :: kr(3,3)

  !> Levi-Civita tensor
  integer :: lvc(3,3,3)

  !> Equation and method parameters
  double precision :: eqpar(neqpar+nspecialpar)

  ! Time step control parameters
  double precision :: courantpar, dtpar, dtdiffpar, dtTCpar{#IFDEF MAGNETOFRICTION ,cmf_c,cmf_y,cmf_divb}
  character(len=std_len) :: typecourant,typeresid
  logical :: addmpibarrier

  !Time parameters

  !> Maximum number of saves that can be defined by tsave or itsave
  integer, parameter :: nsavehi=100

  double precision :: t,tmax,dtmin,residmin,residmax,residual{#IFDEF MAGNETOFRICTION ,tmf}
  double precision :: tfixgrid
  double precision :: tsave(nsavehi,nfile),tsavelast(nfile),dtsave(nfile),slicecoord(nslicemax)
  logical :: tmaxexact,treset,itreset,firstprocess,resetgrid,fixprocess,changeglobals,collapse(ndim)
  integer :: it,itmax,itmin,slowsteps{#IFDEF MAGNETOFRICTION , itmaxmf, ditsavemf}
  integer :: itsave(nsavehi,nfile),itsavelast(nfile),ditsave(nfile)
  integer :: isavet(nfile),isaveit(nfile), nslices, slicedir(nslicemax), collapseLevel
  integer :: n_saves(1:nfile)
  integer :: typeparIO
  integer :: itfixgrid,ditregrid
  integer :: nwauxio
  integer :: istep, nstep

  !Method switches
  character(len=std_len) :: typeadvance
  character(len=std_len) :: typelow1(nlevelshi),typelimited,typesourcesplit
  character(len=std_len) :: typefull1(nlevelshi), typepred1(nlevelshi)
  character(len=std_len) :: typelimiter1(nlevelshi),typegradlimiter1(nlevelshi)
  character(len=std_len) :: typelimiter,typegradlimiter,typeprolonglimit
  character(len=std_len) :: typeentropy(nw),typetvd,typetvdlf,typeaverage
  character(len=std_len) :: typedimsplit,typeaxial,typecoord,typepoly
  integer                :: errorestimate,nxdiffusehllc,typespherical,ncyclemax
  double precision       :: entropycoef(nw)
  double precision       :: tvdlfeps, mcbeta, parastsnu, TCphi
  logical                :: sourceparasts,sourceimpl
  logical                :: sourceimplcycle,conduction,TCsaturate,bcphys
  logical                :: loglimit(nw),logflag(nw),flathllc,flatcd,flatsh,flatppm
  logical                :: ssplitdust,ssplitdivb,ssplitresis,ssplituser,useprimitive,dimsplit
  logical                :: restrictprimitive,prolongprimitive
  logical                :: coarsenprimitive,useprimitiveRel, amrentropy

  logical                :: divbwave,compactres,BnormLF
  double precision       :: divbdiff,smallT,smallp,smallrho,amr_wavefilter(nlevelshi)
  character(len=std_len) :: typedivbdiff,typedivbfix,typediv,typegrad

  logical          :: fixsmall,strictnr,strictsmall,strictzero,strictgetaux
  double precision :: dmaxvel,tolernr,absaccnr,tlow
  integer          :: maxitnr,nflatgetaux

  logical          :: nocartesian
  logical          :: writew(nw),writelevel(nlevelshi)
  double precision :: writespshift(ndim,2)
  integer          :: level_io, level_io_min, level_io_max

  ! cooling related parameters
  integer                :: ncool, cmulti
  character(len=std_len) :: coolcurve,coolmethod
  double precision       :: cfrac
  logical                :: Tfix

  ! dust related paramters
  logical                :: dustzero
  double precision       :: smallrhod
  character(len=std_len) :: dustmethod,dustspecies,dusttemp

  ! local and global fastest wave speed (computed in setdt):
  double precision :: cmax_mype, cmax_global

  !Gravity related parameters
  double precision ::  x1ptms,x2ptms,x3ptms,ptmass

  !Boundary region parameters

  ! Number of boundaries for grid blocks
  integer, parameter :: nhiB=2*ndim
  logical :: periodB(ndim), poleB(2,ndim), aperiodB(ndim)
  character(len=std_len) :: typeB(nw,nhiB)
  character(len=std_len) :: typeghostfill,typegridfill
  double precision ::ratebdflux
  logical :: internalboundary

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
  !> tecplotmpi, tecplotCCmpi, vtumpi, vtuCCmpi, pvtumpi, pvtuCCmpi, tecline,
  !> teclinempi, onegrid
  character(len=std_len) :: convert_type

  !> Data Explorer file endianness ('msb' or 'lsb' for last or most significant
  !> bit order)
  character(len=std_len) :: dxfiletype

  character(len=std_len) :: collapse_type

  !> Conversion factors for length (normvar(0)) and the primitive variables
  !> (normvar(1:nw))
  double precision       :: normvar(0:nw)

  !> Conversion factor for time unit
  double precision       :: normt

  integer                :: saveigrid

  ! Stores the memory and load imbalance, to be used in printlog:
  double precision       :: Xload, Xmemory

  double precision :: time_bc
  {#IFDEF STRETCHGRID
  ! stretching factor qst for log stretch grid
  double precision :: logG, qst
  double precision :: logGs(0:nlevelshi), qsts(0:nlevelshi)
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
  double precision :: rnode(rnodehi,ngridshi)
  double precision :: rnode_sub(rnodehi,ngridshi)

  !> Error tolerance for refinement decision
  double precision :: tol(nlevelshi)
  double precision :: tolratio(nlevelshi)
  double precision :: dx(ndim,nlevelshi)
  double precision :: dt
  double precision :: dtimpl
  double precision :: dt_grid(ngridshi)
  double precision :: dxlevel(ndim)

  integer :: node(nodehi,ngridshi)
  integer :: node_sub(nodehi,ngridshi)

  !> Number of cells as buffer zone
  integer :: nbufferx^D

  !> Maximal number of AMR levels
  integer :: mxnest

  !> Number of ghost cells surrounding a grid
  integer :: dixB
  integer :: levmin
  integer :: levmax
  integer :: levmax_sub

  logical :: skipfinestep
  logical :: phyboundblock(ngridshi)
  logical :: time_advance{#IFDEF MAGNETOFRICTION , mf_advance}

  !$OMP THREADPRIVATE(dxlevel{#IFDEF STRETCHGRID ,logG,qst})
  !$OMP THREADPRIVATE(saveigrid)
  !$OMP THREADPRIVATE(typelimiter,typegradlimiter)
end module mod_global_parameters
