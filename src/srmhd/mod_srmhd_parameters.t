!> Paramter module for Special Relativistic Magneto-hydrodynamics 
module mod_srmhd_parameters
  ! made by Z. MELIANI 14/02/2018
  use mod_global_parameters, only: std_len
  implicit none

    !> Whether thermal conduction is used
  logical, public                  :: srmhd_thermal_conduction = .false.

  !> Whether radiative cooling is added
  logical, public                  :: srmhd_radiative_cooling = .false.

   !> Whether sole energy equation
  logical, public                  :: srmhd_energy = .true.
  !> Whether synge eos  is added
  logical, public                  :: srmhd_eos = .false.

  !> Whether viscosity is added
  logical, public                  :: srmhd_viscosity = .false.

  !> Whether gravity is added
  logical, public                  :: srmhd_gravity = .false.

  !> Whether Hall-MHD is used
  logical, public                  :: srmhd_Hall = .false.

  !> Whether particles module is added
  logical, public                  :: srmhd_particles = .false.

  !> Whether magnetofriction is added
  logical, public                  :: srmhd_magnetofriction = .false.

  !> Whether GLM-MHD is used
  logical, public                  :: srmhd_glm = .false.

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public                  :: source_split_divb = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: srmhd_glm_alpha = 0.5d0

  !> MHD fourth order
  logical, public                  :: srmhd_4th_order = .false.

  !> Number of tracer species
  integer, public                  :: srmhd_n_tracer = 0

  !> Index of the moving frame density (in the w array)
  integer, public                  :: rho_

   !> Index of the lab frame density (in the w array)
  integer, public                  :: d_
  !> Indices of the momentum density
  integer, allocatable, public     :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public                  :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public                  :: p_

  
  !> Indices of the magnetic field
  integer, allocatable, public     :: mag(:)

  !> Indices of the GLM psi
  integer, public     :: psi_
  !> Indices of the Lorentz factor
  integer, public     :: lfac_

  !> Indices of the inertia
  integer, public     :: xi_



  !> Indices of the tracers
  integer, allocatable, public     :: tracer(:)

  !> The adiabatic index
  double precision, public                :: srmhd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: srmhd_adiab = 1.0d0

  !> The MHD resistivity
  double precision, public                :: srmhd_eta = 0.0d0

  !> The MHD hyper-resistivity
  double precision, public                :: srmhd_eta_hyper = 0.0d0

  !> TODO: what is this?
  double precision, public                :: srmhd_etah = 0.0d0


  logical,          public         :: srmhd_checkNR   = .true.
  double precision, public         :: srmhd_absaccnr  = 1.0d-8
  double precision, public         :: srmhd_tolerNR   = 1.0d-9
  double precision, public         :: srmhd_maxdspeed = 1.0d-7
  double precision, public         :: srmhd_maxspeed  = 0.9999
  integer, public                  :: srmhd_maxiterationNR=100

  !> The small_est allowed energy
  double precision                 :: small_e


  !> The small_est allowed inertia
  double precision                 :: small_xi

  ! smaller values for speed
  double precision, public                 :: small_vec2  = 0.0
  !> The number of waves
  integer :: nwwave=8

  !> Method type to clean divergence of B
  character(len=std_len), public     :: typedivbfix  = 'linde'

  !> Method type in a integer for good performance
  integer :: type_divb

  !> Coefficient of diffusive divB cleaning
  double precision :: divbdiff     = 0.8d0

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'

  !> Use a compact way to add resistivity
  logical :: compactres   = .false.

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> Helium abundance over Hydrogen
  double precision, public      :: He_abundance=0.1d0

  !> To control divB=0 fix for boundary
  logical, public     :: boundary_divbfix(2*^ND)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public     :: boundary_divbfix_skip(2*^ND)=0

  !> B0 field is force-free
  logical, public     :: B0field_forcefree=.true.
  ! DivB cleaning methods
  integer, parameter :: divb_none          = 0
  integer, parameter :: divb_glm1          = 1
  integer, parameter :: divb_glm2          = 2
  integer, parameter :: divb_powel         = 3
  integer, parameter :: divb_janhunen      = 4
  integer, parameter :: divb_linde         = 5
  integer, parameter :: divb_lindejanhunen = 6
  integer, parameter :: divb_lindepowel    = 7
  integer, parameter :: divb_lindeglm      = 8


end module mod_srmhd_parameters
