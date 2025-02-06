!> Radiation-magneto-hydrodynamics module
module mod_rmhd_phys

#include "amrvac.h"

  use mod_global_parameters, only: std_len, const_c
  use mod_thermal_conduction, only: tc_fluid
  use mod_thermal_emission, only: te_fluid
  use mod_physics
  use mod_comm_lib, only: mpistop
  use mod_functions_bfield, only: get_divb, mag

  implicit none
  private

  !> The adiabatic index
  double precision, public                :: rmhd_gamma = 5.d0/3.0d0
  !> The adiabatic constant
  double precision, public                :: rmhd_adiab = 1.0d0
  !> The MHD resistivity
  double precision, public                :: rmhd_eta = 0.0d0
  !> The MHD hyper-resistivity
  double precision, public                :: rmhd_eta_hyper = 0.0d0
  !> Hall resistivity
  double precision, public                :: rmhd_etah = 0.0d0
  !> The small_est allowed energy
  double precision, protected             :: small_e
  !> The smallest allowed radiation energy
  double precision, public, protected     :: small_r_e = 0.d0
  !> Height of the mask used in the TRAC method
  double precision, public, protected     :: rmhd_trac_mask = 0.d0
  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: rmhd_glm_alpha = 0.5d0
  !> The thermal conductivity kappa in hyperbolic thermal conduction
  double precision, public                :: hypertc_kappa
  !> Coefficient of diffusive divB cleaning
  double precision                        :: divbdiff     = 0.8d0
  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.1d0
  !> Ionization fraction of H
  !> H_ion_fr = H+/(H+ + H)
  double precision, public, protected  :: H_ion_fr=1d0
  !> Ionization fraction of He
  !> He_ion_fr = (He2+ + He+)/(He2+ + He+ + He)
  double precision, public, protected  :: He_ion_fr=1d0
  !> Ratio of number He2+ / number He+ + He2+
  !> He_ion_fr2 = He2+/(He2+ + He+)
  double precision, public, protected  :: He_ion_fr2=1d0
  ! used for eq of state when it is not defined by units,
  ! the units do not contain terms related to ionization fraction
  ! and it is p = RR * rho * T
  double precision, public, protected  :: RR=1d0
  !> gamma minus one and its inverse
  double precision :: gamma_1, inv_gamma_1
  !> inverse of squared speed of light c0 and reduced speed of light c
  double precision :: inv_squared_c0, inv_squared_c
  !> equi vars indices in the state%equi_vars array
  integer, public :: equi_rho0_ = -1
  integer, public :: equi_pe0_ = -1
  !> Number of tracer species
  integer, public, protected              :: rmhd_n_tracer = 0
  !> Index of the density (in the w array)
  integer, public, protected              :: rho_
  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)
  !> Indices of the momentum density for the form of better vectorization
  integer, public, protected              :: ^C&m^C_
  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_
  !> Index of the radiation energy
  integer, public, protected              :: r_e
  !> Indices of the momentum density for the form of better vectorization
  integer, public, protected              :: ^C&b^C_
  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_
  !> Index of the heat flux q
  integer, public, protected :: q_
  !> Indices of the GLM psi
  integer, public, protected :: psi_
  !> Indices of temperature
  integer, public, protected :: Te_
  !> Index of the cutoff temperature for the TRAC method
  integer, public, protected              :: Tcoff_
  integer, public, protected              :: Tweight_
  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)
  !> The number of waves
  integer :: nwwave=8
  !> Method type in a integer for good performance
  integer :: type_divb
  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public, protected :: boundary_divbfix_skip(2*^ND)=0
  ! DivB cleaning methods
  integer, parameter :: divb_none          = 0
  integer, parameter :: divb_multigrid     = -1
  integer, parameter :: divb_glm           = 1
  integer, parameter :: divb_powel         = 2
  integer, parameter :: divb_janhunen      = 3
  integer, parameter :: divb_linde         = 4
  integer, parameter :: divb_lindejanhunen = 5
  integer, parameter :: divb_lindepowel    = 6
  integer, parameter :: divb_lindeglm      = 7
  integer, parameter :: divb_ct            = 8
  !> Whether an energy equation is used
  logical, public, protected              :: rmhd_energy = .true.
  !> Whether thermal conduction is used
  logical, public, protected              :: rmhd_thermal_conduction = .false.
  !> Whether thermal conduction is used
  logical, public, protected              :: rmhd_hyperbolic_thermal_conduction = .false.
  !> Whether viscosity is added
  logical, public, protected              :: rmhd_viscosity = .false.
  !> Whether gravity is added
  logical, public, protected              :: rmhd_gravity = .false.
  !> Whether particles module is added
  logical, public, protected              :: rmhd_particles = .false.
  !> Whether GLM-MHD is used to control div B
  logical, public, protected              :: rmhd_glm = .false.
  !> Whether extended GLM-MHD is used with additional sources
  logical, public, protected              :: rmhd_glm_extended = .true.
  !> Whether TRAC method is used
  logical, public, protected              :: rmhd_trac = .false.
  !> Which TRAC method is used
  integer, public, protected              :: rmhd_trac_type=1
  !> Distance between two adjacent traced magnetic field lines (in finest cell size)
  integer, public, protected              :: rmhd_trac_finegrid=4
  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.
  !> Whether plasma is partially ionized
  logical, public, protected              :: rmhd_partial_ionization = .false.
  !> Whether CAK radiation line force is activated
  logical, public, protected              :: rmhd_cak_force = .false.
  !> MHD fourth order
  logical, public, protected              :: rmhd_4th_order = .false.
  !> whether split off equilibrium density
  logical, public :: has_equi_rho0 = .false.
  !> whether split off equilibrium thermal pressure
  logical, public :: has_equi_pe0 = .false.
  logical, public :: rmhd_equi_thermal = .false.
  !> whether dump full variables (when splitting is used) in a separate dat file
  logical, public, protected :: rmhd_dump_full_vars = .false.
  !> Whether divB is computed with a fourth order approximation
  integer, public, protected :: rmhd_divb_nth = 1
  !> Use a compact way to add resistivity
  logical :: compactres = .false.
  !> Add divB wave in Roe solver
  logical, public :: divbwave = .true.
  !> clean initial divB
  logical, public :: clean_initial_divb = .false.
  !> Formalism to treat radiation
  character(len=8), public :: rmhd_radiation_formalism = 'fld'
  !> In the case of no rmhd_energy, how to compute pressure
  character(len=8), public :: rmhd_pressure = 'Trad'
  !> Treat radiation fld_Rad_force
  logical, public, protected :: rmhd_radiation_force = .true.
  !> Treat radiation-gas energy interaction
  logical, public, protected :: rmhd_energy_interact = .true.
  !> Treat radiation energy diffusion
  logical, public, protected :: rmhd_radiation_diffusion = .true.
  !> Treat radiation advection
  logical, public, protected :: rmhd_radiation_advection = .true.
  !> Do a running mean over the radiation pressure when determining dt
  logical, protected :: radio_acoustic_filter = .false.
  integer, protected :: size_ra_filter = 1
  !> kb/(m_p mu)* 1/a_rad**4,
  double precision, public :: kbmpmua4
  !> Use the speed of light to calculate the timestep, usefull for debugging
  logical :: dt_c = .false.
  ! remove the below flag  and assume default value = .false.
  ! when eq state properly implemented everywhere
  ! and not anymore through units
  logical, public, protected :: eq_state_units = .true.
  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*^ND)=.true.
  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.
  !> Whether an total energy equation is used
  logical :: total_energy = .true.
  !> Whether an internal or hydrodynamic energy equation is used
  logical, public :: partial_energy = .false.
  !> Whether gravity work is included in energy equation
  logical :: gravity_energy
  !> gravity work is calculated use density times velocity or conservative momentum
  logical :: gravity_rhov = .false.
  !> Method type to clean divergence of B
  character(len=std_len), public, protected :: typedivbfix  = 'linde'
  !> Method type of constrained transport
  character(len=std_len), public, protected :: type_ct  = 'uct_contact'
  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'
  !> type of fluid for thermal conduction
  type(tc_fluid), public, allocatable     :: tc_fl
  !> type of fluid for thermal emission synthesis
  type(te_fluid), public, allocatable     :: te_fl_rmhd

  procedure(sub_convert), pointer      :: rmhd_to_primitive  => null()
  procedure(sub_convert), pointer      :: rmhd_to_conserved  => null()
  procedure(sub_small_values), pointer :: rmhd_handle_small_values => null()
  procedure(sub_get_pthermal), pointer :: rmhd_get_pthermal  => null()
  procedure(sub_get_pthermal), pointer :: rmhd_get_Rfactor   => null()
  procedure(sub_get_pthermal), pointer :: rmhd_get_temperature=> null()
  ! Public methods
  public :: rmhd_phys_init
  public :: rmhd_get_pthermal
  public :: rmhd_get_temperature
  public :: rmhd_get_v
  public :: rmhd_get_rho
  public :: rmhd_to_conserved
  public :: rmhd_to_primitive
  public :: rmhd_e_to_ei
  public :: rmhd_ei_to_e
  public :: rmhd_face_to_center
  public :: get_divb
  public :: get_current
  public :: get_normalized_divb
  public :: b_from_vector_potential
  public :: rmhd_mag_en_all
  {^NOONED
  public :: rmhd_clean_divb_multigrid
  }
  public :: rmhd_get_pradiation
  public :: rmhd_get_pthermal_plus_pradiation
  public :: rmhd_get_tgas
  public :: rmhd_get_trad
  public :: rmhd_set_mg_bounds

contains

  !> Read this module"s parameters from a file
  subroutine rmhd_read_params(files)
    use mod_global_parameters
    use mod_particles, only: particles_eta, particles_etah
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /rmhd_list/ rmhd_energy, rmhd_n_tracer, rmhd_gamma, rmhd_adiab,&
      rmhd_eta, rmhd_eta_hyper, rmhd_etah, rmhd_glm_alpha, rmhd_glm_extended,&
      rmhd_thermal_conduction, rmhd_gravity,&
      rmhd_viscosity, rmhd_4th_order, typedivbfix, source_split_divb, divbdiff,&
      typedivbdiff, type_ct, compactres, divbwave, He_abundance,&
      H_ion_fr, He_ion_fr, He_ion_fr2, eq_state_units, SI_unit, B0field ,rmhd_dump_full_vars,&
      B0field_forcefree, Bdip, Bquad, Boct, Busr, rmhd_particles, rmhd_partial_ionization,&
      particles_eta, particles_etah,has_equi_rho0, has_equi_pe0,rmhd_equi_thermal,&
      boundary_divbfix, boundary_divbfix_skip, rmhd_divb_nth, clean_initial_divb,&
      rmhd_trac, rmhd_trac_type, rmhd_trac_mask, rmhd_trac_finegrid, rmhd_cak_force,&
      rmhd_hyperbolic_thermal_conduction, &
      rmhd_pressure, rmhd_radiation_formalism,&
      rmhd_radiation_force, rmhd_energy_interact, rmhd_radiation_diffusion,&
      rmhd_radiation_advection, radio_acoustic_filter, size_ra_filter, dt_c

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, rmhd_list, end=111)
111    close(unitpar)
    end do
 
  end subroutine rmhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine rmhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer                             :: er
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    character(len=name_len)             :: names(n_par)

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = rmhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine rmhd_write_info

  subroutine rmhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init, particles_eta, particles_etah
    use mod_fld
    use mod_afld
    use mod_supertimestepping, only: sts_init, add_sts_method,&
            set_conversion_methods_to_head, set_error_handling_to_head
    use mod_cak_force, only: cak_init
    use mod_ionization_degree
    use mod_usr_methods, only: usr_Rfactor
    {^NOONED
    use mod_multigrid_coupling
    }

    integer :: itr, idir

    call rmhd_read_params(par_files)

    if(.not.eq_state_units) then
      if(rmhd_partial_ionization) then
        rmhd_partial_ionization=.false.
        if(mype==0) write(*,*) 'WARNING: set rmhd_partial_ionization=F when eq_state_units=F'
      end if
    end if

    if(rmhd_hyperbolic_thermal_conduction) then
      rmhd_thermal_conduction=.false.
      if(mype==0) write(*,*) 'WARNING: turn off parabolic TC when using hyperbolic TC'
    end if

    physics_type = "rmhd"
    phys_energy=rmhd_energy
    phys_trac=rmhd_trac
    phys_trac_type=rmhd_trac_type
    phys_partial_ionization=rmhd_partial_ionization

    phys_gamma = rmhd_gamma
    phys_trac_finegrid=rmhd_trac_finegrid

    if(rmhd_energy) then
      partial_energy=.false.
      total_energy=.true.
    else
      total_energy=.false.
    end if
    phys_total_energy=total_energy
    if(rmhd_energy) then
      gravity_energy=.true.
      if(has_equi_rho0) then
        gravity_rhov=.true.
      end if
    else
      gravity_energy=.false.
    end if

    {^IFONED
    if(rmhd_trac .and. rmhd_trac_type .gt. 2) then
      rmhd_trac_type=1
      if(mype==0) write(*,*) 'WARNING: reset rmhd_trac_type=1 for 1D simulation'
    end if
    }
    if(rmhd_trac .and. rmhd_trac_type .le. 4) then
      rmhd_trac_mask=bigdouble
      if(mype==0) write(*,*) 'WARNING: set rmhd_trac_mask==bigdouble for global TRAC method'
    end if
    phys_trac_mask=rmhd_trac_mask

    ! set default gamma for polytropic/isothermal process
    use_particles=rmhd_particles
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
       type_divb = divb_none
    {^NOONED
    case ('multigrid')
       type_divb = divb_multigrid
       use_multigrid = .true.
       mg%operator_type = mg_laplacian
       phys_global_source_after => rmhd_clean_divb_multigrid
    }
    case ('glm')
      rmhd_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_glm
    case ('powel', 'powell')
      type_divb = divb_powel
    case ('janhunen')
      type_divb = divb_janhunen
    case ('linde')
      type_divb = divb_linde
    case ('lindejanhunen')
      type_divb = divb_lindejanhunen
    case ('lindepowel')
      type_divb = divb_lindepowel
    case ('lindeglm')
      rmhd_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_lindeglm
    case ('ct')
      type_divb = divb_ct
      stagger_grid = .true.
    case default
      call mpistop('Unknown divB fix')
    end select

    allocate(start_indices(number_species),stop_indices(number_species))
    ! set the index of the first flux variable for species 1
    start_indices(1)=1
    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)
    m^C_=mom(^C);

    ! Set index of energy variable
    if (rmhd_energy) then
      nwwave = 8
      e_     = var_set_energy() ! energy density
      p_     = e_               ! gas pressure
    else
      nwwave = 7
      e_     = -1
      p_     = -1
    end if

    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)
    b^C_=mag(^C);

    if (rmhd_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    !> set radiation energy
    r_e = var_set_radiation_energy()

    if(rmhd_hyperbolic_thermal_conduction) then
      ! hyperbolic thermal conduction flux q
      q_ = var_set_q()
      need_global_cmax=.true.
    else
      q_=-1
    end if

    allocate(tracer(rmhd_n_tracer))
    ! Set starting index of tracers
    do itr = 1, rmhd_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    !if(rmhd_hyperbolic_thermal_conduction) then
    !  ! hyperbolic thermal conduction flux q
    !  q_ = var_set_auxvar('q','q')
    !  need_global_cmax=.true.
    !else
    !  q_=-1
    !end if

    !  set temperature as an auxiliary variable to get ionization degree
    if(rmhd_partial_ionization) then
      Te_ = var_set_auxvar('Te','Te')
    else
      Te_ = -1
    end if

    ! set number of variables which need update ghostcells
    nwgc=nwflux+nwaux

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux
  
    !  Number of variables need reconstruction in w
    nw_recon=nwflux

    ! set cutoff temperature when using the TRAC method, as well as an auxiliary weight
    Tweight_ = -1
    if(rmhd_trac) then
      Tcoff_ = var_set_wextra()
      iw_Tcoff=Tcoff_
      if(rmhd_trac_type .ge. 3) then
        Tweight_ = var_set_wextra()
      endif
    else
      Tcoff_ = -1
    end if

    ! set indices of equi vars and update number_equi_vars
    number_equi_vars = 0
    if(has_equi_rho0) then
      number_equi_vars = number_equi_vars + 1
      equi_rho0_ = number_equi_vars
      iw_equi_rho = equi_rho0_
    endif
    if(has_equi_pe0) then
      number_equi_vars = number_equi_vars + 1
      equi_pe0_ = number_equi_vars
      iw_equi_p = equi_pe0_
    endif
    ! determine number of stagger variables
    nws=ndim

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?
    iw_vector(2) = mag(1) - 1   ! TODO: why like this?

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nwflux))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nwflux])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    if(nwflux>mag(ndir)) then
      ! for flux of tracers, using hll flux
      flux_type(:,mag(ndir)+1:nwflux)=flux_hll
    end if

    if(ndim>1) then
      if(rmhd_glm) then
        flux_type(:,psi_)=flux_special
        do idir=1,ndir
          flux_type(idir,mag(idir))=flux_special
        end do
      else
        do idir=1,ndir
          flux_type(idir,mag(idir))=flux_tvdlf
        end do
      end if
    end if
 
    phys_get_rho             => rmhd_get_rho
    phys_get_dt              => rmhd_get_dt
    phys_get_cmax            => rmhd_get_cmax_origin
    phys_get_a2max           => rmhd_get_a2max
    phys_get_tcutoff         => rmhd_get_tcutoff
    phys_get_H_speed         => rmhd_get_H_speed
    if(has_equi_rho0) then
      phys_get_cbounds       => rmhd_get_cbounds_split_rho
    else
      phys_get_cbounds       => rmhd_get_cbounds
    end if
    if(has_equi_rho0) then
      phys_to_primitive      => rmhd_to_primitive_split_rho
      rmhd_to_primitive      => rmhd_to_primitive_split_rho
      phys_to_conserved      => rmhd_to_conserved_split_rho
      rmhd_to_conserved      => rmhd_to_conserved_split_rho
    else
      phys_to_primitive      => rmhd_to_primitive_origin
      rmhd_to_primitive      => rmhd_to_primitive_origin
      phys_to_conserved      => rmhd_to_conserved_origin
      rmhd_to_conserved      => rmhd_to_conserved_origin
    end if
    if(B0field.or.has_equi_rho0.or.has_equi_pe0) then
      phys_get_flux            => rmhd_get_flux_split
    else
      phys_get_flux            => rmhd_get_flux
    end if
    phys_get_v                 => rmhd_get_v
    if(B0field.or.has_equi_rho0) then
      phys_add_source_geom     => rmhd_add_source_geom_split
    else
      phys_add_source_geom     => rmhd_add_source_geom
    end if
    phys_add_source          => rmhd_add_source
    phys_check_params        => rmhd_check_params
    phys_write_info          => rmhd_write_info
 
    phys_handle_small_values => rmhd_handle_small_values_origin
    rmhd_handle_small_values => rmhd_handle_small_values_origin
    phys_check_w             => rmhd_check_w_origin

    phys_set_mg_bounds       => rmhd_set_mg_bounds
    phys_get_trad            => rmhd_get_trad
    phys_get_tgas            => rmhd_get_tgas
 
    phys_get_pthermal        => rmhd_get_pthermal_origin
    rmhd_get_pthermal         => rmhd_get_pthermal_origin

    if(number_equi_vars>0) then
      phys_set_equi_vars => set_equi_vars_grid
    endif

    if(type_divb==divb_glm) then
      phys_modify_wLR => rmhd_modify_wLR
    end if

    ! choose Rfactor in ideal gas law
    if(rmhd_partial_ionization) then
      rmhd_get_Rfactor=>Rfactor_from_temperature_ionization
      phys_update_temperature => rmhd_update_temperature
    else if(associated(usr_Rfactor)) then
      rmhd_get_Rfactor=>usr_Rfactor
    else
      rmhd_get_Rfactor=>Rfactor_from_constant_ionization
    end if

    if(rmhd_partial_ionization) then
      rmhd_get_temperature => rmhd_get_temperature_from_Te
    else
      if(has_equi_pe0 .and. has_equi_rho0) then
        rmhd_get_temperature => rmhd_get_temperature_from_etot_with_equi
      else
        rmhd_get_temperature => rmhd_get_temperature_from_etot
      end if
    end if

    ! if using ct stagger grid, boundary divb=0 is not done here
    if(stagger_grid) then
      phys_get_ct_velocity => rmhd_get_ct_velocity
      phys_update_faces => rmhd_update_faces
      phys_face_to_center => rmhd_face_to_center
      phys_modify_wLR => rmhd_modify_wLR
    else if(ndim>1) then
      phys_boundary_adjust => rmhd_boundary_adjust
    end if

    {^NOONED
    ! clean initial divb
    if(clean_initial_divb) phys_clean_divb => rmhd_clean_divb_multigrid
    }

    ! derive units from basic units
    call rmhd_physical_units()

    !> Initiate radiation-closure module
    select case(rmhd_radiation_formalism)
    case('fld')
      call fld_init(He_abundance, rmhd_radiation_diffusion, rmhd_energy_interact, rmhd_gamma)
    case('afld')
      call afld_init(He_abundance, rmhd_radiation_diffusion, rmhd_gamma)
    case default
      call mpistop('Radiation formalism unknown')
    end select

    if(rmhd_hyperbolic_thermal_conduction) then
      hypertc_kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
    end if
    if(.not. rmhd_energy .and. rmhd_thermal_conduction) then
      call mpistop("thermal conduction needs rmhd_energy=T")
    end if
    if(.not. rmhd_energy .and. rmhd_hyperbolic_thermal_conduction) then
      call mpistop("hyperbolic thermal conduction needs rmhd_energy=T")
    end if

    ! initialize thermal conduction module
    if (rmhd_thermal_conduction) then
      call sts_init()
      call tc_init_params(rmhd_gamma)

      allocate(tc_fl)
      call tc_get_mhd_params(tc_fl,tc_params_read_rmhd)
      call add_sts_method(rmhd_get_tc_dt_rmhd,rmhd_sts_set_source_tc_rmhd,e_,1,e_,1,.false.)
      if(has_equi_pe0 .and. has_equi_rho0) then
        tc_fl%get_temperature_from_conserved => rmhd_get_temperature_from_etot_with_equi
      else
        tc_fl%get_temperature_from_conserved => rmhd_get_temperature_from_etot
      end if
      if(has_equi_pe0 .and. has_equi_rho0) then
        tc_fl%get_temperature_from_eint => rmhd_get_temperature_from_eint_with_equi
        if(rmhd_equi_thermal) then
          tc_fl%has_equi = .true.
          tc_fl%get_temperature_equi => rmhd_get_temperature_equi
          tc_fl%get_rho_equi => rmhd_get_rho_equi
        else
          tc_fl%has_equi = .false.
        end if
      else
        tc_fl%get_temperature_from_eint => rmhd_get_temperature_from_eint
      end if
      call set_conversion_methods_to_head(rmhd_e_to_ei, rmhd_ei_to_e)
      call set_error_handling_to_head(rmhd_tc_handle_small_e)
      tc_fl%get_rho => rmhd_get_rho
      tc_fl%e_ = e_
      tc_fl%Tcoff_ = Tcoff_
    end if

    allocate(te_fl_rmhd)
    te_fl_rmhd%get_rho=> rmhd_get_rho
    te_fl_rmhd%get_pthermal=> rmhd_get_pthermal
    te_fl_rmhd%get_var_Rfactor => rmhd_get_Rfactor
{^IFTHREED
    phys_te_images => rmhd_te_images
}
    ! Initialize viscosity module
    if (rmhd_viscosity) call viscosity_init(phys_wider_stencil)

    ! Initialize gravity module
    if(rmhd_gravity) then
      call gravity_init()
    end if

    ! Initialize particles module
    if(rmhd_particles) then
      call particles_init()
      if (particles_eta  < zero) particles_eta = rmhd_eta
      if (particles_etah < zero) particles_eta = rmhd_etah
      if(mype==0) then
         write(*,*) '*****Using particles:        with rmhd_eta, rmhd_etah :', rmhd_eta, rmhd_etah
         write(*,*) '*****Using particles: particles_eta, particles_etah :', particles_eta, particles_etah
      end if
    end if

    ! initialize ionization degree table
    if(rmhd_partial_ionization) call ionization_degree_init()

    ! Initialize CAK radiation force module
    if (rmhd_cak_force) call cak_init(rmhd_gamma)
  end subroutine rmhd_phys_init

{^IFTHREED
  subroutine rmhd_te_images
    use mod_global_parameters
    use mod_thermal_emission

    select case(convert_type)
      case('EIvtiCCmpi','EIvtuCCmpi')
        call get_EUV_image(unitconvert,te_fl_rmhd)
      case('ESvtiCCmpi','ESvtuCCmpi')
        call get_EUV_spectrum(unitconvert,te_fl_rmhd)
      case('SIvtiCCmpi','SIvtuCCmpi')
        call get_SXR_image(unitconvert,te_fl_rmhd)
      case('WIvtiCCmpi','WIvtuCCmpi')
        call get_whitelight_image(unitconvert,te_fl_rmhd)
      case default
        call mpistop("Error in synthesize emission: Unknown convert_type")
      end select
  end subroutine rmhd_te_images
}

!!start th cond
  ! wrappers for STS functions in thermal_conductivity module
  ! which take as argument the tc_fluid (defined in the physics module)
  subroutine  rmhd_sts_set_source_tc_rmhd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_mhd
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_mhd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,tc_fl)
  end subroutine rmhd_sts_set_source_tc_rmhd

  function rmhd_get_tc_dt_rmhd(w,ixI^L,ixO^L,dx^D,x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_mhd
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_mhd(w,ixI^L,ixO^L,dx^D,x,tc_fl)
  end function rmhd_get_tc_dt_rmhd

  subroutine rmhd_tc_handle_small_e(w, x, ixI^L, ixO^L, step)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer, intent(in)    :: step
    character(len=140) :: error_msg

    write(error_msg,"(a,i3)") "Thermal conduction step ", step
    call rmhd_handle_small_ei(w,x,ixI^L,ixO^L,e_,error_msg)
  end subroutine rmhd_tc_handle_small_e

  ! fill in tc_fluid fields from namelist
  subroutine tc_params_read_rmhd(fl)
    use mod_global_parameters, only: unitpar,par_files
    type(tc_fluid), intent(inout) :: fl
    double precision :: tc_k_para=0d0
    double precision :: tc_k_perp=0d0
    integer                      :: n
    ! list parameters
    logical :: tc_perpendicular=.false.
    logical :: tc_saturate=.false.
    character(len=std_len)  :: tc_slope_limiter="MC"

    namelist /tc_list/ tc_perpendicular, tc_saturate, tc_slope_limiter, tc_k_para, tc_k_perp

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, tc_list, end=111)
111     close(unitpar)
    end do

    fl%tc_perpendicular = tc_perpendicular
    fl%tc_saturate = tc_saturate
    fl%tc_k_para = tc_k_para
    fl%tc_k_perp = tc_k_perp
    select case(tc_slope_limiter)
     case ('no','none')
       fl%tc_slope_limiter = 0
     case ('MC')
       ! montonized central limiter Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       fl%tc_slope_limiter = 1
     case('minmod')
       ! minmod limiter
       fl%tc_slope_limiter = 2
     case ('superbee')
       ! Roes superbee limiter (eq.3.51i)
       fl%tc_slope_limiter = 3
     case ('koren')
       ! Barry Koren Right variant
       fl%tc_slope_limiter = 4
     case default
       call mpistop("Unknown tc_slope_limiter, choose MC, minmod")
    end select
  end subroutine tc_params_read_rmhd
!!end th cond

  !> sets the equilibrium variables
  subroutine set_equi_vars_grid_faces(igrid,x,ixI^L,ixO^L)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in) :: igrid, ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: delx(ixI^S,1:ndim)
    double precision :: xC(ixI^S,1:ndim),xshift^D
    integer :: idims, ixC^L, hxO^L, ix, idims2

    if(slab_uniform)then
      ^D&delx(ixI^S,^D)=rnode(rpdx^D_,igrid)\
    else
      ! for all non-cartesian and stretched cartesian coordinates
      delx(ixI^S,1:ndim)=ps(igrid)%dx(ixI^S,1:ndim)
    endif

    do idims=1,ndim
      hxO^L=ixO^L-kr(idims,^D);
      if(stagger_grid) then
        ! ct needs all transverse cells
        ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D); ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
      else
        ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
        ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
      end if
      ! always xshift=0 or 1/2
      xshift^D=half*(one-kr(^D,idims));
      do idims2=1,ndim
        select case(idims2)
        {case(^D)
          do ix = ixC^LIM^D
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ix^D%ixC^S,^D)=x(ix^D%ixC^S,^D)+(half-xshift^D)*delx(ix^D%ixC^S,^D)
          end do\}
        end select
      end do
      call usr_set_equi_vars(ixI^L,ixC^L,xC,ps(igrid)%equi_vars(ixI^S,1:number_equi_vars,idims))
    end do
  end subroutine set_equi_vars_grid_faces

  !> sets the equilibrium variables
  subroutine set_equi_vars_grid(igrid)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in) :: igrid

    !values at the center
    call usr_set_equi_vars(ixG^LL,ixG^LL,ps(igrid)%x,ps(igrid)%equi_vars(ixG^T,1:number_equi_vars,0))

    !values at the interfaces
    call set_equi_vars_grid_faces(igrid,ps(igrid)%x,ixG^LL,ixM^LL)
  end subroutine set_equi_vars_grid

  ! w, wnew conserved, add splitted variables back to wnew
  function convert_vars_splitting(ixI^L,ixO^L, w, x, nwc) result(wnew)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L,ixO^L, nwc
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision   :: wnew(ixO^S, 1:nwc)

    if(has_equi_rho0) then
      wnew(ixO^S,rho_)=w(ixO^S,rho_)+block%equi_vars(ixO^S,equi_rho0_,0)
    else
      wnew(ixO^S,rho_)=w(ixO^S,rho_)
    endif
    wnew(ixO^S,mom(:))=w(ixO^S,mom(:))

    if (B0field) then
      ! add background magnetic field B0 to B
      wnew(ixO^S,mag(1:ndir))=w(ixO^S,mag(1:ndir))+block%B0(ixO^S,1:ndir,0)
    else
      wnew(ixO^S,mag(1:ndir))=w(ixO^S,mag(1:ndir))
    end if

    if(rmhd_energy) then
      wnew(ixO^S,e_)=w(ixO^S,e_)
      if(has_equi_pe0) then
        wnew(ixO^S,e_)=wnew(ixO^S,e_)+block%equi_vars(ixO^S,equi_pe0_,0)*inv_gamma_1
      end if
      if(B0field .and. total_energy) then
        wnew(ixO^S,e_)=wnew(ixO^S,e_)+0.5d0*sum(block%B0(ixO^S,:,0)**2,dim=ndim+1) &
            + sum(w(ixO^S,mag(:))*block%B0(ixO^S,:,0),dim=ndim+1)
      end if
    end if
  end function convert_vars_splitting

  subroutine rmhd_check_params
    use mod_global_parameters
    use mod_usr_methods
    use mod_convert, only: add_convert_method

    ! after user parameter setting
    gamma_1=rmhd_gamma-1.d0
    if (.not. rmhd_energy) then
       if (rmhd_gamma <= 0.0d0) call mpistop ("Error: rmhd_gamma <= 0")
       if (rmhd_adiab < 0.0d0) call mpistop ("Error: rmhd_adiab < 0")
       small_pressure = rmhd_adiab*small_density**rmhd_gamma
    else
       if (rmhd_gamma <= 0.0d0 .or. rmhd_gamma == 1.0d0) &
            call mpistop ("Error: rmhd_gamma <= 0 or rmhd_gamma == 1")
       inv_gamma_1=1.d0/gamma_1
       small_e = small_pressure * inv_gamma_1
    end if

    small_r_e = small_pressure/(rmhd_gamma - 1.0d0)

    if (number_equi_vars > 0 .and. .not. associated(usr_set_equi_vars)) then
      call mpistop("usr_set_equi_vars has to be implemented in the user file")
    endif
    if(convert .or. autoconvert) then
      if(convert_type .eq. 'dat_generic_mpi') then
        if(rmhd_dump_full_vars) then
          if(mype .eq. 0) print*, " add conversion method: split -> full "
          call add_convert_method(convert_vars_splitting, nw, cons_wnames, "new")
        endif
      endif
    endif

    if (use_multigrid) call rmhd_set_mg_bounds()
  end subroutine rmhd_check_params

  !> Set the boundaries for the diffusion of E
  subroutine rmhd_set_mg_bounds
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_usr_methods
    integer :: iB

    ! Set boundary conditions for the multigrid solver
    do iB = 1, 2*ndim
       select case (typeboundary(r_e, iB))
       case (bc_symm)
       ! d/dx u = 0
           mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
           mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case (bc_asymm)
           ! u = 0
           mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
           mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case (bc_cont)
           ! d/dx u = 0
           ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
           mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
           mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case (bc_periodic)
           ! Nothing to do here
       case (bc_noinflow)
           call usr_special_mg_bc(iB)
       case (bc_special)
           call usr_special_mg_bc(iB)
       case default
           call mpistop("divE_multigrid warning: unknown b.c. ")
       end select
    end do
  end subroutine rmhd_set_mg_bounds

  subroutine rmhd_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0,c_lightspeed
    double precision :: a,b

    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
      c_lightspeed=c_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi ! G^2 cm^2 dyne^-1
      c_lightspeed=const_c
    end if
    if(eq_state_units) then
      a = 1d0 + 4d0 * He_abundance
      if(rmhd_partial_ionization) then
        b = 2+.3d0
      else
        b = 1d0 + H_ion_fr + He_abundance*(He_ion_fr*(He_ion_fr2 + 1d0)+1d0)
      end if
      RR = 1d0
    else
      a = 1d0
      b = 1d0
      RR = (1d0 + H_ion_fr + He_abundance*(He_ion_fr*(He_ion_fr2 + 1d0)+1d0))/(1d0 + 4d0 * He_abundance)
    end if
    if(unit_density/=1.d0) then
      unit_numberdensity=unit_density/(a*mp)
    else
      ! unit of numberdensity is independent by default
      unit_density=a*mp*unit_numberdensity
    end if
    if(unit_velocity/=1.d0) then
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    else if(unit_pressure/=1.d0) then
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      unit_velocity=sqrt(unit_pressure/unit_density)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    else if(unit_magneticfield/=1.d0) then
      unit_pressure=unit_magneticfield**2/miu0
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      unit_velocity=sqrt(unit_pressure/unit_density)
    else if(unit_temperature/=1.d0) then
      unit_pressure=b*unit_numberdensity*kB*unit_temperature
      unit_velocity=sqrt(unit_pressure/unit_density)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    end if
    if(unit_time/=1.d0) then
      unit_length=unit_time*unit_velocity
    else
      ! unit of length is independent by default
      unit_time=unit_length/unit_velocity
    end if
    ! Additional units needed for the particles
    c_norm=c_lightspeed/unit_velocity
    unit_charge=unit_magneticfield*unit_length**2/unit_velocity/miu0
    if (.not. SI_unit) unit_charge = unit_charge*const_c
    unit_mass=unit_density*unit_length**3

    !> Units for radiative flux and opacity
    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)
  end subroutine rmhd_physical_units

  subroutine rmhd_check_w_origin(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters
    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    logical, intent(inout) :: flag(ixI^S,1:nw)
    double precision :: tmp
    integer :: ix^D

    flag=.false.
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if(has_equi_rho0) then
        tmp=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0)
      else
        tmp=w(ix^D,rho_)
      end if
      if(tmp<small_density) flag(ix^D,rho_) = .true.
      if(primitive) then
        if(has_equi_pe0) then
          if(w(ix^D,p_)+block%equi_vars(ix^D,equi_pe0_,0)<small_pressure) flag(ix^D,e_) = .true.
        else
          if(w(ix^D,p_)<small_pressure) flag(ix^D,e_) = .true.
        end if
      else
        tmp=w(ix^D,e_)-half*((^C&w(ix^D,m^C_)**2+)/tmp+(^C&w(ix^D,b^C_)**2+))
        if(has_equi_pe0) then
          if(tmp+block%equi_vars(ix^D,equi_pe0_,0)*inv_gamma_1<small_e) flag(ix^D,e_) = .true.
        else
          if(tmp<small_e) flag(ix^D,e_) = .true.
        end if
      end if
      if(w(ix^D,r_e)<small_r_e) flag(ix^D,r_e) = .true.
    {end do\}
  end subroutine rmhd_check_w_origin

  !> Transform primitive variables into conservative ones
  subroutine rmhd_to_conserved_origin(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer :: ix^D

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! Calculate total energy from pressure, kinetic and magnetic energy
      w(ix^D,e_)=w(ix^D,p_)*inv_gamma_1&
                 +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                 +(^C&w(ix^D,b^C_)**2+))
      ! Convert velocity to momentum
      ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
    {end do\}
  end subroutine rmhd_to_conserved_origin

  !> Transform primitive variables into conservative ones
  subroutine rmhd_to_conserved_split_rho(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision :: rho
    integer :: ix^D

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      rho=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i)
      ! Calculate total energy from pressure, kinetic and magnetic energy
      w(ix^D,e_)=w(ix^D,p_)*inv_gamma_1&
                 +half*((^C&w(ix^D,m^C_)**2+)*rho&
                       +(^C&w(ix^D,b^C_)**2+))
      ! Convert velocity to momentum
      ^C&w(ix^D,m^C_)=rho*w(ix^D,m^C_)\
    {end do\}
  end subroutine rmhd_to_conserved_split_rho

  !> Transform conservative variables into primitive ones
  subroutine rmhd_to_primitive_origin(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision                :: inv_rho
    integer :: ix^D

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call rmhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'rmhd_to_primitive_origin')
    end if

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      inv_rho = 1.d0/w(ix^D,rho_)
      ! Convert momentum to velocity
      ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
      ! Calculate pressure = (gamma-1) * (e-ek-eb)
      w(ix^D,p_)=gamma_1*(w(ix^D,e_)&
                -half*(w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+)&
                  +(^C&w(ix^D,b^C_)**2+)))
    {end do\}
  end subroutine rmhd_to_primitive_origin

  !> Transform conservative variables into primitive ones
  subroutine rmhd_to_primitive_split_rho(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision :: inv_rho
    integer :: ix^D

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call rmhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'rmhd_to_primitive_split_rho')
    end if

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      inv_rho=1.d0/(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
      ! Convert momentum to velocity
      ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
      ! Calculate pressure = (gamma-1) * (e-ek-eb)
      w(ix^D,p_)=gamma_1*(w(ix^D,e_)&
                  -half*((w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))*&
                  (^C&w(ix^D,m^C_)**2+)+(^C&w(ix^D,b^C_)**2+)))
    {end do\}
  end subroutine rmhd_to_primitive_split_rho

  !> Transform internal energy to total energy
  subroutine rmhd_ei_to_e(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer :: ix^D

    if(has_equi_rho0) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Calculate e = ei + ek + eb
        w(ix^D,e_)=w(ix^D,e_)&
                  +half*((^C&w(ix^D,m^C_)**2+)/&
         (w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0))&
                  +(^C&w(ix^D,b^C_)**2+))
      {end do\}
    else
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Calculate e = ei + ek + eb
        w(ix^D,e_)=w(ix^D,e_)&
                  +half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)&
                        +(^C&w(ix^D,b^C_)**2+))
      {end do\}
    end if
  end subroutine rmhd_ei_to_e

  !> Transform total energy to internal energy
  subroutine rmhd_e_to_ei(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer :: ix^D

    if(has_equi_rho0) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Calculate ei = e - ek - eb
        w(ix^D,e_)=w(ix^D,e_)&
                  -half*((^C&w(ix^D,m^C_)**2+)/&
         (w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0))&
                        +(^C&w(ix^D,b^C_)**2+))
      {end do\}
    else
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Calculate ei = e - ek - eb
        w(ix^D,e_)=w(ix^D,e_)&
                  -half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)&
                        +(^C&w(ix^D,b^C_)**2+))
      {end do\}
    end if

    if(fix_small_values) then
      call rmhd_handle_small_ei(w,x,ixI^L,ixI^L,e_,'rmhd_e_to_ei')
    end if
  end subroutine rmhd_e_to_ei

  subroutine rmhd_handle_small_values_origin(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    double precision :: rho
    integer :: idir, ix^D
    logical :: flag(ixI^S,1:nw)

    call phys_check_w(primitive, ixI^L, ixI^L, w, flag)
    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(has_equi_rho0) then
            rho=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0)
            if(flag(ix^D,rho_)) w(ix^D,rho_)=small_density-block%equi_vars(ix^D,equi_rho0_,0)
          else
            rho=w(ix^D,rho_)
            if(flag(ix^D,rho_)) w(ix^D,rho_)=small_density
          end if
         {
          if(small_values_fix_iw(m^C_)) then
            if(flag({ix^D},rho_)) w({ix^D},m^C_)=0.0d0
          end if
          \}
          if(primitive) then
            if(has_equi_pe0) then
              if(flag(ix^D,e_)) w(ix^D,p_)=small_pressure-block%equi_vars(ix^D,equi_pe0_,0)
            else
              if(flag(ix^D,e_)) w(ix^D,p_)=small_pressure
            end if
          else
            if(has_equi_pe0) then
              if(flag(ix^D,e_)) &
                w(ix^D,e_)=small_e+half*((^C&w(ix^D,m^C_)**2+)/rho+(^C&w(ix^D,b^C_)**2+))&
                -block%equi_vars(ix^D,equi_pe0_,0)*inv_gamma_1
            else
              if(flag(ix^D,e_)) &
                w(ix^D,e_)=small_e+half*((^C&w(ix^D,m^C_)**2+)/rho+(^C&w(ix^D,b^C_)**2+))
            end if
          end if
          if(flag(ix^D,r_e)) w(ix^D,r_e)=small_r_e
       {end do\}
      case ("average")
        ! do averaging of density
        call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
        if(primitive)then
          call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
        else
          ! do averaging of internal energy
          {do ix^DB=ixImin^DB,ixImax^DB\}
            if(has_equi_rho0) then
              rho=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0)
            else
              rho=w(ix^D,rho_)
            end if
            w(ix^D,e_)=w(ix^D,e_)&
                -half*((^C&w(ix^D,m^C_)**2+)/rho+(^C&w(ix^D,b^C_)**2+))
          {end do\}
          call small_values_average(ixI^L, ixO^L, w, x, flag, e_)
          ! convert back
          {do ix^DB=ixImin^DB,ixImax^DB\}
            if(has_equi_rho0) then
              rho=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0)
            else
              rho=w(ix^D,rho_)
            end if
            w(ix^D,e_)=w(ix^D,e_)&
                +half*((^C&w(ix^D,m^C_)**2+)/rho+(^C&w(ix^D,b^C_)**2+))
          {end do\}
        end if
        call small_values_average(ixI^L, ixO^L, w, x, flag, r_e)
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! do averaging of internal energy
         {do ix^DB=ixImin^DB,ixImax^DB\}
            if(has_equi_rho0) then
              rho=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0)
            else
              rho=w(ix^D,rho_)
            end if
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)/rho\
            w(ix^D,p_)=gamma_1*(w(ix^D,e_)&
                -half*((^C&w(ix^D,m^C_)**2+)/rho+(^C&w(ix^D,b^C_)**2+)))
         {end do\}
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
  end subroutine rmhd_handle_small_values_origin

  !> Calculate v vector
  subroutine rmhd_get_v(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,ndir)
    double precision :: rho(ixI^S)
    integer :: idir

    call rmhd_get_rho(w,x,ixI^L,ixO^L,rho)
    rho(ixO^S)=1.d0/rho(ixO^S)
    ! Convert momentum to velocity
    do idir = 1, ndir
       v(ixO^S, idir) = w(ixO^S, mom(idir))*rho(ixO^S)
    end do
  end subroutine rmhd_get_v

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine rmhd_get_cmax_origin(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision :: rho, inv_rho, cfast2, AvMinCs2, b2, kmax
    integer :: ix^D

    if(B0field) then
     {do ix^DB=ixOmin^DB,ixOmax^DB \}
        if(has_equi_rho0) then
          rho=(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
        else
          rho=w(ix^D,rho_)
        end if
        inv_rho=1.d0/rho
        ! sound speed**2 
        cmax(ix^D)=rmhd_gamma*w(ix^D,p_)*inv_rho
        ! store |B|^2 in v
        b2=(^C&(w(ix^D,b^C_)+block%B0(ix^D,^C,b0i))**2+)
        cfast2=b2*inv_rho+cmax(ix^D)
        AvMinCs2=cfast2**2-4.0d0*cmax(ix^D)*(w(ix^D,mag(idim))+block%B0(ix^D,idim,b0i))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        cmax(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
        cmax(ix^D)=abs(w(ix^D,mom(idim)))+cmax(ix^D)
     {end do\}
    else
     {do ix^DB=ixOmin^DB,ixOmax^DB \}
        if(has_equi_rho0) then
          rho=(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
        else
          rho=w(ix^D,rho_)
        end if
        inv_rho=1.d0/rho
        ! sound speed**2 
        cmax(ix^D)=rmhd_gamma*w(ix^D,p_)*inv_rho
        ! store |B|^2 in v
        b2=(^C&w(ix^D,b^C_)**2+)
        cfast2=b2*inv_rho+cmax(ix^D)
        AvMinCs2=cfast2**2-4.0d0*cmax(ix^D)*w(ix^D,mag(idim))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        cmax(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
        cmax(ix^D)=abs(w(ix^D,mom(idim)))+cmax(ix^D)
     {end do\}
    end if
  end subroutine rmhd_get_cmax_origin

  subroutine rmhd_get_a2max(w,x,ixI^L,ixO^L,a2max)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: a2max(ndim)
    double precision :: a2(ixI^S,ndim,nw)
    integer :: gxO^L,hxO^L,jxO^L,kxO^L,i,j

    a2=zero
    do i = 1,ndim
      !> 4th order
      hxO^L=ixO^L-kr(i,^D);
      gxO^L=hxO^L-kr(i,^D);
      jxO^L=ixO^L+kr(i,^D);
      kxO^L=jxO^L+kr(i,^D);
      a2(ixO^S,i,1:nw)=abs(-w(kxO^S,1:nw)+16.d0*w(jxO^S,1:nw)&
         -30.d0*w(ixO^S,1:nw)+16.d0*w(hxO^S,1:nw)-w(gxO^S,1:nw))
      a2max(i)=maxval(a2(ixO^S,i,1:nw))/12.d0/dxlevel(i)**2
    end do
  end subroutine rmhd_get_a2max

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine rmhd_get_tcutoff(ixI^L,ixO^L,w,x,Tco_local,Tmax_local)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: Tco_local,Tmax_local
    ! in primitive form
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, parameter :: trac_delta=0.25d0
    double precision :: tmp1(ixI^S),Te(ixI^S),lts(ixI^S)
    double precision, dimension(ixI^S,1:ndir) :: bunitvec
    double precision, dimension(ixI^S,1:ndim) :: gradT
    double precision :: Bdir(ndim)
    double precision :: ltrc,ltrp,altr(ixI^S)
    integer :: idims,jxO^L,hxO^L,ixA^D,ixB^D,ix^D
    integer :: jxP^L,hxP^L,ixP^L,ixQ^L
    logical :: lrlt(ixI^S)

    if(rmhd_partial_ionization) then
      call rmhd_get_temperature_from_Te(w,x,ixI^L,ixI^L,Te)
    else
      call rmhd_get_Rfactor(w,x,ixI^L,ixI^L,Te)
      Te(ixI^S)=w(ixI^S,p_)/(Te(ixI^S)*w(ixI^S,rho_))
    end if
    Tco_local=zero
    Tmax_local=maxval(Te(ixO^S))

    {^IFONED
    select case(rmhd_trac_type)
    case(0)
      !> test case, fixed cutoff temperature
      block%wextra(ixI^S,Tcoff_)=2.5d5/unit_temperature
    case(1)
      hxO^L=ixO^L-1;
      jxO^L=ixO^L+1;
      lts(ixO^S)=0.5d0*abs(Te(jxO^S)-Te(hxO^S))/Te(ixO^S)
      lrlt=.false.
      where(lts(ixO^S) > trac_delta)
        lrlt(ixO^S)=.true.
      end where
      if(any(lrlt(ixO^S))) then
        Tco_local=maxval(Te(ixO^S), mask=lrlt(ixO^S))
      end if
    case(2)
      !> iijima et al. 2021, LTRAC method
      ltrc=1.5d0
      ltrp=4.d0
      ixP^L=ixO^L^LADD1;
      hxO^L=ixO^L-1;
      jxO^L=ixO^L+1;
      hxP^L=ixP^L-1;
      jxP^L=ixP^L+1;
      lts(ixP^S)=0.5d0*abs(Te(jxP^S)-Te(hxP^S))/Te(ixP^S)
      lts(ixP^S)=max(one, (exp(lts(ixP^S))/ltrc)**ltrp)
      lts(ixO^S)=0.25d0*(lts(jxO^S)+two*lts(ixO^S)+lts(hxO^S))
      block%wextra(ixO^S,Tcoff_)=Te(ixO^S)*lts(ixO^S)**0.4d0
    case default
      call mpistop("rmhd_trac_type not allowed for 1D simulation")
    end select
    }
    {^NOONED
    select case(rmhd_trac_type)
    case(0)
      !> test case, fixed cutoff temperature
      block%wextra(ixI^S,Tcoff_)=2.5d5/unit_temperature
    case(1,4,6)
      ! temperature gradient at cell centers
      do idims=1,ndim
        call gradient(Te,ixI^L,ixO^L,idims,tmp1)
        gradT(ixO^S,idims)=tmp1(ixO^S)
      end do
      ! B vector
      if(B0field) then
        bunitvec(ixO^S,:)=w(ixO^S,iw_mag(:))+block%B0(ixO^S,:,0)
      else
        bunitvec(ixO^S,:)=w(ixO^S,iw_mag(:))
      end if
      if(rmhd_trac_type .gt. 1) then
        ! B direction at cell center
        Bdir=zero
        {do ixA^D=0,1\}
          ixB^D=(ixOmin^D+ixOmax^D-1)/2+ixA^D;
          Bdir(1:ndim)=Bdir(1:ndim)+bunitvec(ixB^D,1:ndim)
        {end do\}
        if(sum(Bdir(:)**2) .gt. zero) then
          Bdir(1:ndim)=Bdir(1:ndim)/dsqrt(sum(Bdir(:)**2))
        end if
        block%special_values(3:ndim+2)=Bdir(1:ndim)
      end if
      tmp1(ixO^S)=dsqrt(sum(bunitvec(ixO^S,:)**2,dim=ndim+1))
      where(tmp1(ixO^S)/=0.d0)
        tmp1(ixO^S)=1.d0/tmp1(ixO^S)
      elsewhere
        tmp1(ixO^S)=bigdouble
      end where
      ! b unit vector: magnetic field direction vector
      do idims=1,ndim
        bunitvec(ixO^S,idims)=bunitvec(ixO^S,idims)*tmp1(ixO^S)
      end do
      ! temperature length scale inversed
      lts(ixO^S)=abs(sum(gradT(ixO^S,1:ndim)*bunitvec(ixO^S,1:ndim),dim=ndim+1))/Te(ixO^S)
      ! fraction of cells size to temperature length scale
      if(slab_uniform) then
        lts(ixO^S)=minval(dxlevel)*lts(ixO^S)
      else
        lts(ixO^S)=minval(block%ds(ixO^S,:),dim=ndim+1)*lts(ixO^S)
      end if
      lrlt=.false.
      where(lts(ixO^S) > trac_delta)
        lrlt(ixO^S)=.true.
      end where
      if(any(lrlt(ixO^S))) then
        block%special_values(1)=maxval(Te(ixO^S), mask=lrlt(ixO^S))
      else
        block%special_values(1)=zero
      end if
      block%special_values(2)=Tmax_local
    case(2)
      !> iijima et al. 2021, LTRAC method
      ltrc=1.5d0
      ltrp=4.d0
      ixP^L=ixO^L^LADD2;
      ! temperature gradient at cell centers
      do idims=1,ndim
        ixQ^L=ixP^L;
        hxP^L=ixP^L;
        jxP^L=ixP^L;
        select case(idims)
          {case(^D)
           ixQmin^D=ixQmin^D+1
           ixQmax^D=ixQmax^D-1
           hxPmax^D=ixPmin^D
           jxPmin^D=ixPmax^D
          \}
        end select
        call gradient(Te,ixI^L,ixQ^L,idims,gradT(ixI^S,idims))
        call gradientF(Te,x,ixI^L,hxP^L,idims,gradT(ixI^S,idims),nghostcells,.true.)
        call gradientF(Te,x,ixI^L,jxP^L,idims,gradT(ixI^S,idims),nghostcells,.false.)
      end do
      ! B vector
     {do ix^DB=ixPmin^DB,ixPmax^DB\}
        if(B0field) then
          ^C&bunitvec(ix^D,^C)=w(ix^D,iw_mag(^C))+block%B0(ix^D,^C,0)\
        else
          ^C&bunitvec(ix^D,^C)=w(ix^D,iw_mag(^C))\
        end if
        tmp1(ix^D)=1.d0/(dsqrt(^C&bunitvec(ix^D,^C)**2+)+smalldouble)
        ! b unit vector: magnetic field direction vector
        ^D&bunitvec({ix^D},^D)=bunitvec({ix^D},^D)*tmp1({ix^D})\
        ! temperature length scale inversed
        lts(ix^D)=abs(^D&gradT({ix^D},^D)*bunitvec({ix^D},^D)+)/Te(ix^D)
        ! fraction of cells size to temperature length scale
        if(slab_uniform) then
          lts(ix^D)=min(^D&dxlevel(^D))*lts(ix^D)
        else
          lts(ix^D)=min(^D&block%ds({ix^D},^D))*lts(ix^D)
        end if
        lts(ix^D)=max(one,(exp(lts(ix^D))/ltrc)**ltrp)
     {end do\}
      ! need one ghost layer for thermal conductivity
      ixP^L=ixO^L^LADD1;
      do idims=1,ndim
        hxO^L=ixP^L-kr(idims,^D);
        jxO^L=ixP^L+kr(idims,^D);
        if(idims==1) then
          altr(ixP^S)=0.25d0*(lts(hxO^S)+two*lts(ixP^S)+lts(jxO^S))*bunitvec(ixP^S,idims)**2
        else
          altr(ixP^S)=altr(ixP^S)+0.25d0*(lts(hxO^S)+two*lts(ixP^S)+lts(jxO^S))*bunitvec(ixP^S,idims)**2
        end if
      end do
      block%wextra(ixP^S,Tcoff_)=Te(ixP^S)*altr(ixP^S)**0.4d0
    case(3,5)
      !> do nothing here
    case default
      call mpistop("unknown rmhd_trac_type")
    end select
    }
  end subroutine rmhd_get_tcutoff

  !> get H speed for H-correction to fix the carbuncle problem at grid-aligned shock front
  subroutine rmhd_get_H_speed(wprim,x,ixI^L,ixO^L,idim,Hspeed)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wprim(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: Hspeed(ixI^S,1:number_species)

    double precision :: csound(ixI^S,ndim)
    double precision, allocatable :: tmp(:^D&)
    integer :: jxC^L, ixC^L, ixA^L, id, ix^D

    Hspeed=0.d0
    ixA^L=ixO^L^LADD1;
    allocate(tmp(ixA^S))
    do id=1,ndim
      call rmhd_get_csound_prim(wprim,x,ixI^L,ixA^L,id,tmp)
      csound(ixA^S,id)=tmp(ixA^S)
    end do
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D+kr(idim,^D)-1;
    jxCmax^D=ixCmax^D+kr(idim,^D);
    jxCmin^D=ixCmin^D+kr(idim,^D);
    Hspeed(ixC^S,1)=0.5d0*abs(wprim(jxC^S,mom(idim))+csound(jxC^S,idim)-wprim(ixC^S,mom(idim))+csound(ixC^S,idim))

    do id=1,ndim
      if(id==idim) cycle
      ixAmax^D=ixCmax^D+kr(id,^D);
      ixAmin^D=ixCmin^D+kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(wprim(ixA^S,mom(id))+csound(ixA^S,id)-wprim(ixC^S,mom(id))+csound(ixC^S,id)))
      ixAmax^D=ixCmax^D-kr(id,^D);
      ixAmin^D=ixCmin^D-kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(wprim(ixC^S,mom(id))+csound(ixC^S,id)-wprim(ixA^S,mom(id))+csound(ixA^S,id)))
    end do

    do id=1,ndim
      if(id==idim) cycle
      ixAmax^D=jxCmax^D+kr(id,^D);
      ixAmin^D=jxCmin^D+kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(wprim(ixA^S,mom(id))+csound(ixA^S,id)-wprim(jxC^S,mom(id))+csound(jxC^S,id)))
      ixAmax^D=jxCmax^D-kr(id,^D);
      ixAmin^D=jxCmin^D-kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(wprim(jxC^S,mom(id))+csound(jxC^S,id)-wprim(ixA^S,mom(id))+csound(ixA^S,id)))
    end do
    deallocate(tmp)

  end subroutine rmhd_get_H_speed

  !> Estimating bounds for the minimum and maximum signal velocities without split
  subroutine rmhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,Hspeed,cmax,cmin)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
    double precision, intent(in)    :: Hspeed(ixI^S,1:number_species)

    double precision :: wmean(ixI^S,nw), csoundL(ixO^S), csoundR(ixO^S)
    double precision :: umean, dmean, tmp1, tmp2, tmp3
    integer :: ix^D

    select case (boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      call rmhd_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call rmhd_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          tmp1=sqrt(wLp(ix^D,rho_))
          tmp2=sqrt(wRp(ix^D,rho_))
          tmp3=1.d0/(tmp1+tmp2)
          umean=(wLp(ix^D,mom(idim))*tmp1+wRp(ix^D,mom(idim))*tmp2)*tmp3
          dmean=sqrt((tmp1*csoundL(ix^D)**2+tmp2*csoundR(ix^D)**2)*tmp3+&
           half*tmp1*tmp2*tmp3**2*(wRp(ix^D,mom(idim))-wLp(ix^D,mom(idim)))**2)
          cmin(ix^D,1)=umean-dmean
          cmax(ix^D,1)=umean+dmean
       {end do\}
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          tmp1=sqrt(wLp(ix^D,rho_))
          tmp2=sqrt(wRp(ix^D,rho_))
          tmp3=1.d0/(tmp1+tmp2)
          umean=(wLp(ix^D,mom(idim))*tmp1+wRp(ix^D,mom(idim))*tmp2)*tmp3
          dmean=sqrt((tmp1*csoundL(ix^D)**2+tmp2*csoundR(ix^D)**2)*tmp3+&
           half*tmp1*tmp2*tmp3**2*(wRp(ix^D,mom(idim))-wLp(ix^D,mom(idim)))**2)
          cmax(ix^D,1)=abs(umean)+dmean
       {end do\}
      end if
    case (2)
      wmean(ixO^S,1:nwflux)=0.5d0*(wLp(ixO^S,1:nwflux)+wRp(ixO^S,1:nwflux))
      call rmhd_get_csound_prim(wmean,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          cmax(ix^D,1)=max(wmean(ix^D,mom(idim))+csoundR(ix^D),zero)
          cmin(ix^D,1)=min(wmean(ix^D,mom(idim))-csoundR(ix^D),zero)
       {end do\}
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=abs(wmean(ixO^S,mom(idim)))+csoundR(ixO^S)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call rmhd_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call rmhd_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          csoundL(ix^D)=max(csoundL(ix^D),csoundR(ix^D))
          cmin(ix^D,1)=min(wLp(ix^D,mom(idim)),wRp(ix^D,mom(idim)))-csoundL(ix^D)
          cmax(ix^D,1)=max(wLp(ix^D,mom(idim)),wRp(ix^D,mom(idim)))+csoundL(ix^D)
       {end do\}
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          csoundL(ix^D)=max(csoundL(ix^D),csoundR(ix^D))
          cmax(ix^D,1)=max(wLp(ix^D,mom(idim)),wRp(ix^D,mom(idim)))+csoundL(ix^D)
       {end do\}
      end if
    end select
  end subroutine rmhd_get_cbounds

  !> Estimating bounds for the minimum and maximum signal velocities with rho split
  subroutine rmhd_get_cbounds_split_rho(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,Hspeed,cmax,cmin)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
    double precision, intent(in)    :: Hspeed(ixI^S,1:number_species)
    double precision :: wmean(ixI^S,nw), csoundL(ixO^S), csoundR(ixO^S)
    double precision :: umean, dmean, tmp1, tmp2, tmp3
    integer :: ix^D

    select case (boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      call rmhd_get_csound_prim_split(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call rmhd_get_csound_prim_split(wRp,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          tmp1=sqrt(wLp(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
          tmp2=sqrt(wRp(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
          tmp3=1.d0/(tmp1+tmp2)
          umean=(wLp(ix^D,mom(idim))*tmp1+wRp(ix^D,mom(idim))*tmp2)*tmp3
          dmean=sqrt((tmp1*csoundL(ix^D)**2+tmp2*csoundR(ix^D)**2)*tmp3+&
           half*tmp1*tmp2*tmp3**2*(wRp(ix^D,mom(idim))-wLp(ix^D,mom(idim)))**2)
          cmin(ix^D,1)=umean-dmean
          cmax(ix^D,1)=umean+dmean
       {end do\}
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          tmp1=sqrt(wLp(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
          tmp2=sqrt(wRp(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
          tmp3=1.d0/(tmp1+tmp2)
          umean=(wLp(ix^D,mom(idim))*tmp1+wRp(ix^D,mom(idim))*tmp2)*tmp3
          dmean=sqrt((tmp1*csoundL(ix^D)**2+tmp2*csoundR(ix^D)**2)*tmp3+&
           half*tmp1*tmp2*tmp3**2*(wRp(ix^D,mom(idim))-wLp(ix^D,mom(idim)))**2)
          cmax(ix^D,1)=abs(umean)+dmean
       {end do\}
      end if
    case (2)
      wmean(ixO^S,1:nwflux)=0.5d0*(wLp(ixO^S,1:nwflux)+wRp(ixO^S,1:nwflux))
      call rmhd_get_csound_prim_split(wmean,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          cmax(ix^D,1)=max(wmean(ix^D,mom(idim))+csoundR(ix^D),zero)
          cmin(ix^D,1)=min(wmean(ix^D,mom(idim))-csoundR(ix^D),zero)
       {end do\}
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=abs(wmean(ixO^S,mom(idim)))+csoundR(ixO^S)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call rmhd_get_csound_prim_split(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call rmhd_get_csound_prim_split(wRp,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          csoundL(ix^D)=max(csoundL(ix^D),csoundR(ix^D))
          cmin(ix^D,1)=min(wLp(ix^D,mom(idim)),wRp(ix^D,mom(idim)))-csoundL(ix^D)
          cmax(ix^D,1)=max(wLp(ix^D,mom(idim)),wRp(ix^D,mom(idim)))+csoundL(ix^D)
       {end do\}
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          csoundL(ix^D)=max(csoundL(ix^D),csoundR(ix^D))
          cmax(ix^D,1)=max(wLp(ix^D,mom(idim)),wRp(ix^D,mom(idim)))+csoundL(ix^D)
       {end do\}
      end if
    end select
  end subroutine rmhd_get_cbounds_split_rho

  !> prepare velocities for ct methods
  subroutine rmhd_get_ct_velocity(vcts,wLp,wRp,ixI^L,ixO^L,idim,cmax,cmin)
    use mod_global_parameters
    integer, intent(in)              :: ixI^L, ixO^L, idim
    double precision, intent(in)     :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)     :: cmax(ixI^S)
    double precision, intent(in), optional :: cmin(ixI^S)
    type(ct_velocity), intent(inout) :: vcts
    integer                          :: idimE,idimN

    ! calculate velocities related to different UCT schemes
    select case(type_ct)
    case('average')
    case('uct_contact')
      if(.not.allocated(vcts%vnorm)) allocate(vcts%vnorm(ixI^S,1:ndim))
      ! get average normal velocity at cell faces
      vcts%vnorm(ixO^S,idim)=0.5d0*(wLp(ixO^S,mom(idim))+wRp(ixO^S,mom(idim)))
    case('uct_hll')
      if(.not.allocated(vcts%vbarC)) then
        allocate(vcts%vbarC(ixI^S,1:ndir,2),vcts%vbarLC(ixI^S,1:ndir,2),vcts%vbarRC(ixI^S,1:ndir,2))
        allocate(vcts%cbarmin(ixI^S,1:ndim),vcts%cbarmax(ixI^S,1:ndim))
      end if
      ! Store magnitude of characteristics
      if(present(cmin)) then
        vcts%cbarmin(ixO^S,idim)=max(-cmin(ixO^S),zero)
        vcts%cbarmax(ixO^S,idim)=max( cmax(ixO^S),zero)
      else
        vcts%cbarmax(ixO^S,idim)=max( cmax(ixO^S),zero)
        vcts%cbarmin(ixO^S,idim)=vcts%cbarmax(ixO^S,idim)
      end if

      idimN=mod(idim,ndir)+1 ! 'Next' direction
      idimE=mod(idim+1,ndir)+1 ! Electric field direction
      ! Store velocities
      vcts%vbarLC(ixO^S,idim,1)=wLp(ixO^S,mom(idimN))
      vcts%vbarRC(ixO^S,idim,1)=wRp(ixO^S,mom(idimN))
      vcts%vbarC(ixO^S,idim,1)=(vcts%cbarmax(ixO^S,idim)*vcts%vbarLC(ixO^S,idim,1) &
           +vcts%cbarmin(ixO^S,idim)*vcts%vbarRC(ixO^S,idim,1))&
          /(vcts%cbarmax(ixO^S,idim)+vcts%cbarmin(ixO^S,idim))
      vcts%vbarLC(ixO^S,idim,2)=wLp(ixO^S,mom(idimE))
      vcts%vbarRC(ixO^S,idim,2)=wRp(ixO^S,mom(idimE))
      vcts%vbarC(ixO^S,idim,2)=(vcts%cbarmax(ixO^S,idim)*vcts%vbarLC(ixO^S,idim,2) &
           +vcts%cbarmin(ixO^S,idim)*vcts%vbarRC(ixO^S,idim,1))&
          /(vcts%cbarmax(ixO^S,idim)+vcts%cbarmin(ixO^S,idim))
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select
  end subroutine rmhd_get_ct_velocity

  !> Calculate fast magnetosonic wave speed
  subroutine rmhd_get_csound_prim(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixO^S)
    double precision :: inv_rho, cfast2, AvMinCs2, b2, kmax
    double precision :: prad_tensor(ixO^S, 1:ndim, 1:ndim)
    double precision :: prad_max(ixO^S)
    integer :: ix^D

    call rmhd_get_pradiation(w, x, ixI^L, ixO^L, prad_tensor, nghostcells-1)
    !> filter cmax
    if(radio_acoustic_filter) then
      call rmhd_radio_acoustic_filter(x, ixI^L, ixO^L, prad_max)
    endif
    ! store |B|^2 in v
    if(B0field) then
      {do ix^DB=ixOmin^DB,ixOmax^DB \}
        inv_rho=1.d0/w(ix^D,rho_)
        prad_max(ix^D) = maxval(prad_tensor(ix^D,:,:))
        if(rmhd_energy) then
          csound(ix^D)=max(rmhd_gamma,4.d0/3.d0)*(w(ix^D,p_)+prad_max(ix^D))*inv_rho
        else
          csound(ix^D)=rmhd_gamma*rmhd_adiab*w(ix^D,rho_)**gamma_1
        end if
        b2=(^C&(w(ix^D,b^C_)+block%B0(ix^D,^C,b0i))**2+)
        cfast2=b2*inv_rho+csound(ix^D)
        AvMinCs2=cfast2**2-4.0d0*csound(ix^D)*(w(ix^D,mag(idim))+&
         block%B0(ix^D,idim,b0i))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        csound(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
      {end do\}
    else
      {do ix^DB=ixOmin^DB,ixOmax^DB \}
        inv_rho=1.d0/w(ix^D,rho_)
        prad_max(ix^D)=maxval(prad_tensor(ix^D,:,:))
        if(rmhd_energy) then
          csound(ix^D)=max(rmhd_gamma,4.d0/3.d0)*(w(ix^D,p_)+prad_max(ix^D))*inv_rho
        else
          csound(ix^D)=rmhd_gamma*rmhd_adiab*w(ix^D,rho_)**gamma_1
        end if
        b2=(^C&w(ix^D,b^C_)**2+)
        cfast2=b2*inv_rho+csound(ix^D)
        AvMinCs2=cfast2**2-4.0d0*csound(ix^D)*w(ix^D,mag(idim))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        csound(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
      {end do\}
    end if
  end subroutine rmhd_get_csound_prim

  !> Calculate fast magnetosonic wave speed
  subroutine rmhd_get_csound_prim_split(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixO^S)
    double precision :: rho, inv_rho, cfast2, AvMinCs2, b2, kmax
    double precision :: prad_tensor(ixO^S, 1:ndim, 1:ndim)
    double precision :: prad_max(ixO^S)
    integer :: ix^D

    call rmhd_get_pradiation(w, x, ixI^L, ixO^L, prad_tensor, nghostcells-1)
    !> filter cmax
    if (radio_acoustic_filter) then
      call rmhd_radio_acoustic_filter(x, ixI^L, ixO^L, prad_max)
    endif

    ! store |B|^2 in v
    if(B0field) then
      {do ix^DB=ixOmin^DB,ixOmax^DB \}
        rho=(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
        inv_rho=1.d0/rho
        prad_max(ix^D) = maxval(prad_tensor(ix^D,:,:))
        if(has_equi_pe0) then
          csound(ix^D)=max(rmhd_gamma,4.d0/3.d0)*(w(ix^D,p_)+block%equi_vars(ix^D,equi_pe0_,b0i)+prad_max(ix^D))*inv_rho
        end if
        b2=(^C&(w(ix^D,b^C_)+block%B0(ix^D,^C,b0i))**2+)
        cfast2=b2*inv_rho+csound(ix^D)
        AvMinCs2=cfast2**2-4.0d0*csound(ix^D)*(w(ix^D,mag(idim))+&
         block%B0(ix^D,idim,b0i))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        csound(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
      {end do\}
    else
      {do ix^DB=ixOmin^DB,ixOmax^DB \}
        rho=(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
        inv_rho=1.d0/rho
        prad_max(ix^D) = maxval(prad_tensor(ix^D,:,:))
        if(has_equi_pe0) then
          csound(ix^D)=max(rmhd_gamma,4.d0/3.d0)*(w(ix^D,p_)+block%equi_vars(ix^D,equi_pe0_,b0i)+prad_max(ix^D))*inv_rho
        end if
        b2=(^C&w(ix^D,b^C_)**2+)
        cfast2=b2*inv_rho+csound(ix^D)
        AvMinCs2=cfast2**2-4.0d0*csound(ix^D)*w(ix^D,mag(idim))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        csound(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
      {end do\}
    end if
  end subroutine rmhd_get_csound_prim_split

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine rmhd_get_pthermal_origin(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)

    integer :: iw, ix^D

   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if(has_equi_rho0) then
        pth(ix^D)=gamma_1*(w(ix^D,e_)-half*((^C&w(ix^D,m^C_)**2+)/(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0))&
             +(^C&w(ix^D,b^C_)**2+)))+block%equi_vars(ix^D,equi_pe0_,0)
      else
        pth(ix^D)=gamma_1*(w(ix^D,e_)-half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)&
             +(^C&w(ix^D,b^C_)**2+)))
      end if
      if(fix_small_values.and.pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
   {end do\}

    if(check_small_values.and..not.fix_small_values) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(pth(ix^D)<small_pressure) then
          write(*,*) "Error: small value of gas pressure",pth(ix^D),&
               " encountered when call rmhd_get_pthermal"
          write(*,*) "Iteration: ", it, " Time: ", global_time
          write(*,*) "Location: ", x(ix^D,:)
          write(*,*) "Cell number: ", ix^D
          do iw=1,nw
            write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
          end do
          ! use erroneous arithmetic operation to crash the run
          if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
          write(*,*) "Saving status at the previous time step"
          crash=.true.
        end if
     {end do\}
    end if
  end subroutine rmhd_get_pthermal_origin

  !> copy temperature from stored Te variable
  subroutine rmhd_get_temperature_from_Te(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S) = w(ixO^S, Te_)
  end subroutine rmhd_get_temperature_from_Te

  !> Calculate temperature=p/rho when in e_ the internal energy is stored
  subroutine rmhd_get_temperature_from_eint(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    double precision :: R(ixI^S)

    call rmhd_get_Rfactor(w,x,ixI^L,ixO^L,R)
    res(ixO^S) = gamma_1 * w(ixO^S, e_)/(w(ixO^S,rho_)*R(ixO^S))
  end subroutine rmhd_get_temperature_from_eint

  !> Calculate temperature=p/rho when in e_ the total energy is stored
  subroutine rmhd_get_temperature_from_etot(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    double precision :: R(ixI^S)

    call rmhd_get_Rfactor(w,x,ixI^L,ixO^L,R)
    call rmhd_get_pthermal(w,x,ixI^L,ixO^L,res)
    res(ixO^S)=res(ixO^S)/(R(ixO^S)*w(ixO^S,rho_))
  end subroutine rmhd_get_temperature_from_etot

  subroutine rmhd_get_temperature_from_etot_with_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    double precision :: R(ixI^S)

    call rmhd_get_Rfactor(w,x,ixI^L,ixO^L,R)
    call rmhd_get_pthermal(w,x,ixI^L,ixO^L,res)
    res(ixO^S)=res(ixO^S)/(R(ixO^S)*(w(ixO^S,rho_)+block%equi_vars(ixO^S,equi_rho0_,b0i)))
  end subroutine rmhd_get_temperature_from_etot_with_equi

  subroutine rmhd_get_temperature_from_eint_with_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    double precision :: R(ixI^S)

    call rmhd_get_Rfactor(w,x,ixI^L,ixO^L,R)
    res(ixO^S) = (gamma_1 * w(ixO^S, e_) + block%equi_vars(ixO^S,equi_pe0_,b0i)) /&
                ((w(ixO^S,rho_) +block%equi_vars(ixO^S,equi_rho0_,b0i))*R(ixO^S))
  end subroutine rmhd_get_temperature_from_eint_with_equi

  subroutine rmhd_get_temperature_equi(w,x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    double precision :: R(ixI^S)

    call rmhd_get_Rfactor(w,x,ixI^L,ixO^L,R)
    res(ixO^S)= block%equi_vars(ixO^S,equi_pe0_,b0i)/(block%equi_vars(ixO^S,equi_rho0_,b0i)*R(ixO^S))
  end subroutine rmhd_get_temperature_equi

  subroutine rmhd_get_rho_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S) = block%equi_vars(ixO^S,equi_rho0_,b0i)
  end subroutine rmhd_get_rho_equi

  subroutine rmhd_get_pe_equi(w,x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S) = block%equi_vars(ixO^S,equi_pe0_,b0i)
  end subroutine rmhd_get_pe_equi

  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine rmhd_get_p_total(w,x,ixI^L,ixO^L,p)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: p(ixI^S)

    call rmhd_get_pthermal(w,x,ixI^L,ixO^L,p)
    p(ixO^S) = p(ixO^S) + 0.5d0 * sum(w(ixO^S, mag(:))**2, dim=ndim+1)
  end subroutine rmhd_get_p_total

  !> Calculate radiation pressure within ixO^L
  subroutine rmhd_get_pradiation(w, x, ixI^L, ixO^L, prad, nth)
    use mod_global_parameters
    use mod_fld
    use mod_afld
    integer, intent(in)          :: ixI^L, ixO^L, nth
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: prad(ixO^S, 1:ndim, 1:ndim)

    select case (rmhd_radiation_formalism)
    case('fld')
      call fld_get_radpress(w, x, ixI^L, ixO^L, prad, nth)
    case('afld')
      call afld_get_radpress(w, x, ixI^L, ixO^L, prad, nth)
    case default
      call mpistop('Radiation formalism unknown')
    end select
  end subroutine rmhd_get_pradiation

  !> Calculates the sum of the gas pressure and the max Prad tensor element
  subroutine rmhd_get_pthermal_plus_pradiation(w, x, ixI^L, ixO^L, pth_plus_prad)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, 1:nw)
    double precision, intent(in)  :: x(ixI^S, 1:ndim)
    double precision              :: pth(ixI^S)
    double precision              :: prad_tensor(ixO^S, 1:ndim, 1:ndim)
    double precision              :: prad_max(ixO^S)
    double precision, intent(out) :: pth_plus_prad(ixI^S)
    integer :: ix^D

    call rmhd_get_pthermal(w, x, ixI^L, ixO^L, pth)
    call rmhd_get_pradiation(w, x, ixI^L, ixO^L, prad_tensor, nghostcells)
    {do ix^D = ixOmin^D,ixOmax^D\}
      prad_max(ix^D) = maxval(prad_tensor(ix^D,:,:))
    {enddo\}
    !> filter cmax
    if (radio_acoustic_filter) then
      call rmhd_radio_acoustic_filter(x, ixI^L, ixO^L, prad_max)
    endif
    pth_plus_prad(ixO^S) = pth(ixO^S) + prad_max(ixO^S)
  end subroutine rmhd_get_pthermal_plus_pradiation

  !> Filter peaks in cmax due to radiation energy density, used for debugging
  subroutine rmhd_radio_acoustic_filter(x, ixI^L, ixO^L, prad_max)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: prad_max(ixO^S)
    double precision :: tmp_prad(ixI^S)
    integer :: ix^D, filter, idim

    if (size_ra_filter .lt. 1) call mpistop("ra filter of size < 1 makes no sense")
    if (size_ra_filter .gt. nghostcells) call mpistop("ra filter of size < nghostcells makes no sense")

    tmp_prad(ixI^S) = zero
    tmp_prad(ixO^S) = prad_max(ixO^S)
    do filter = 1,size_ra_filter
      do idim = 1,ndim
        ! {do ix^D = ixOmin^D+filter,ixOmax^D-filter\}
        {do ix^D = ixOmin^D,ixOmax^D\}
            prad_max(ix^D) = min(tmp_prad(ix^D),tmp_prad(ix^D+filter*kr(idim,^D)))
            prad_max(ix^D) = min(tmp_prad(ix^D),tmp_prad(ix^D-filter*kr(idim,^D)))
        {enddo\}
      enddo
    enddo
  end subroutine rmhd_radio_acoustic_filter

  !> Calculates gas temperature
  subroutine rmhd_get_tgas(w, x, ixI^L, ixO^L, tgas)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: pth(ixI^S)
    double precision, intent(out):: tgas(ixI^S)

    call rmhd_get_pthermal(w, x, ixI^L, ixO^L, pth)
        tgas(ixI^S) = pth(ixI^S)/w(ixI^S,rho_)
  end subroutine rmhd_get_tgas

  !> Calculates radiation temperature
  subroutine rmhd_get_trad(w, x, ixI^L, ixO^L, trad)
    use mod_global_parameters
    use mod_constants

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: trad(ixI^S)

    trad(ixI^S) = (w(ixI^S,r_e)*unit_pressure&
    /const_rad_a)**(1.d0/4.d0)/unit_temperature
  end subroutine rmhd_get_trad

  !> Calculate fluxes within ixO^L without any splitting
  subroutine rmhd_get_flux(wC,w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)
    double precision             :: vHall(ixI^S,1:ndir)
    double precision             :: ptotal
    integer                      :: iw, ix^D

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       ! Get flux of density
       f(ix^D,rho_)=w(ix^D,mom(idim))*w(ix^D,rho_)
       ! f_i[m_k]=v_i*m_k-b_k*b_i
       ^C&f(ix^D,m^C_)=wC(ix^D,mom(idim))*w(ix^D,m^C_)-w(ix^D,mag(idim))*w(ix^D,b^C_)\
       ptotal=w(ix^D,p_)+half*(^C&w(ix^D,b^C_)**2+)
       ! normal one includes total pressure
       f(ix^D,mom(idim))=f(ix^D,mom(idim))+ptotal
       ! Get flux of total energy
       ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
       f(ix^D,e_)=w(ix^D,mom(idim))*(wC(ix^D,e_)+ptotal)&
          -w(ix^D,mag(idim))*(^C&w(ix^D,b^C_)*w(ix^D,m^C_)+)
       ! f_i[b_k]=v_i*b_k-v_k*b_i
       ^C&f(ix^D,b^C_)=w(ix^D,mom(idim))*w(ix^D,b^C_)-w(ix^D,mag(idim))*w(ix^D,m^C_)\
    {end do\}
    if(rmhd_glm) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,mag(idim))=w(ix^D,psi_)
        !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
        f(ix^D,psi_)=cmax_global**2*w(ix^D,mag(idim))
     {end do\}
    end if
    if(rmhd_radiation_advection) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,r_e)=w(ix^D,mom(idim))*wC(ix^D,r_e)
      {end do\}
    else
      f(ixO^S,r_e)=zero
    endif
    ! Get flux of tracer
    do iw=1,rmhd_n_tracer
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,tracer(iw))=w(ix^D,mom(idim))*w(ix^D,tracer(iw))
     {end do\}
    end do
    if(rmhd_hyperbolic_thermal_conduction) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,e_)=f(ix^D,e_)+w(ix^D,q_)*w(ix^D,mag(idim))/(dsqrt(^D&w({ix^D},b^D_)**2+)+smalldouble)
        f(ix^D,q_)=zero
     {end do\}
    end if
  end subroutine rmhd_get_flux

  !> Calculate fluxes within ixO^L with possible splitting
  subroutine rmhd_get_flux_split(wC,w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)
    double precision             :: vHall(ixI^S,1:ndir)
    double precision             :: ptotal, Btotal(ixO^S,1:ndir)
    integer                      :: iw, ix^D

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! Get flux of density
      if(has_equi_rho0) then
        f(ix^D,rho_)=w(ix^D,mom(idim))*(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
      else
        f(ix^D,rho_)=w(ix^D,mom(idim))*w(ix^D,rho_)
      endif
      if(rmhd_energy) then
        ptotal=w(ix^D,p_)+half*(^C&w(ix^D,b^C_)**2+)
      else
        ptotal=rmhd_adiab*w(ix^D,rho_)**rmhd_gamma+half*(^C&w(ix^D,b^C_)**2+)
        if(has_equi_pe0) then
          ptotal=ptotal-block%equi_vars(ix^D,equi_pe0_,b0i)
        end if
      end if
      if(B0field) then
        ^C&btotal(ix^D,^C)=w(ix^D,b^C_)+block%B0(ix^D,^C,idim)\
        ptotal=ptotal+(^C&w(ix^D,b^C_)*block%B0(ix^D,^C,idim)+)
        ! Get flux of momentum and magnetic field
        ! f_i[m_k]=v_i*m_k-b_k*b_i
        ^C&f(ix^D,m^C_)=wC(ix^D,mom(idim))*w(ix^D,m^C_)-&
                          Btotal(ix^D,idim)*w(ix^D,b^C_)-w(ix^D,mag(idim))*block%B0(ix^D,^C,idim)\
        f(ix^D,mom(idim))=f(ix^D,mom(idim))+ptotal
      else
        ^C&btotal(ix^D,^C)=w(ix^D,b^C_)\
        ! Get flux of momentum and magnetic field
        ! f_i[m_k]=v_i*m_k-b_k*b_i
        ^C&f(ix^D,m^C_)=wC(ix^D,mom(idim))*w(ix^D,m^C_)-w(ix^D,mag(idim))*w(ix^D,b^C_)\
        f(ix^D,mom(idim))=f(ix^D,mom(idim))+ptotal
      end if
      ! f_i[b_k]=v_i*b_k-v_k*b_i
      ^C&f(ix^D,b^C_)=w(ix^D,mom(idim))*Btotal(ix^D,^C)-Btotal(ix^D,idim)*w(ix^D,m^C_)\
      ! Get flux of energy
      ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
      if(rmhd_energy) then
        f(ix^D,e_)=w(ix^D,mom(idim))*(wC(ix^D,e_)+ptotal)&
           -Btotal(ix^D,idim)*(^C&w(ix^D,b^C_)*w(ix^D,m^C_)+)
      end if
    {end do\}
    if(rmhd_glm) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,mag(idim))=w(ix^D,psi_)
        !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
        f(ix^D,psi_) = cmax_global**2*w(ix^D,mag(idim))
      {end do\}
    end if
    if(rmhd_radiation_advection) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,r_e)=w(ix^D,mom(idim))*wC(ix^D,r_e)
      {end do\}
    else
      f(ixO^S,r_e)=zero
    endif
    ! Get flux of tracer
    do iw=1,rmhd_n_tracer
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,tracer(iw))=w(ix^D,mom(idim))*w(ix^D,tracer(iw))
     {end do\}
    end do
    if(rmhd_hyperbolic_thermal_conduction) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,e_)=f(ix^D,e_)+w(ix^D,q_)*Btotal(ix^D,idim)/(dsqrt(^C&btotal(ix^D,^C)**2+)+smalldouble)
        f(ix^D,q_)=zero
     {end do\}
    end if
  end subroutine rmhd_get_flux_split

  !> use cell-center flux to get cell-face flux
  !> and get the source term as the divergence of the flux
  subroutine get_flux_on_cell_face(ixI^L,ixO^L,ff,src)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, dimension(:^D&,:), intent(inout) :: ff
    double precision, intent(out) :: src(ixI^S)

    double precision :: ffc(ixI^S,1:ndim)
    double precision :: dxinv(ndim)
    integer :: idims, ix^D, ixA^L, ixB^L, ixC^L

    ixA^L=ixO^L^LADD1;
    dxinv=1.d0/dxlevel
    ! cell corner flux in ffc
    ffc=0.d0
    ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;
    {do ix^DB=0,1\}
      ixBmin^D=ixCmin^D+ix^D;
      ixBmax^D=ixCmax^D+ix^D;
      ffc(ixC^S,1:ndim)=ffc(ixC^S,1:ndim)+ff(ixB^S,1:ndim)
    {end do\}
    ffc(ixC^S,1:ndim)=0.5d0**ndim*ffc(ixC^S,1:ndim)
    ! flux at cell face
    ff(ixI^S,1:ndim)=0.d0
    do idims=1,ndim
      ixB^L=ixO^L-kr(idims,^D);
      ixCmax^D=ixOmax^D; ixCmin^D=ixBmin^D;
      {do ix^DB=0,1 \}
         if({ ix^D==0 .and. ^D==idims | .or.}) then
           ixBmin^D=ixCmin^D-ix^D;
           ixBmax^D=ixCmax^D-ix^D;
           ff(ixC^S,idims)=ff(ixC^S,idims)+ffc(ixB^S,idims)
         end if
      {end do\}
      ff(ixC^S,idims)=ff(ixC^S,idims)*0.5d0**(ndim-1)
    end do
    src=0.d0
    if(slab_uniform) then
      do idims=1,ndim
        ff(ixA^S,idims)=dxinv(idims)*ff(ixA^S,idims)
        ixB^L=ixO^L-kr(idims,^D);
        src(ixO^S)=src(ixO^S)+ff(ixO^S,idims)-ff(ixB^S,idims)
      end do
    else
      do idims=1,ndim
        ff(ixA^S,idims)=ff(ixA^S,idims)*block%surfaceC(ixA^S,idims)
        ixB^L=ixO^L-kr(idims,^D);
        src(ixO^S)=src(ixO^S)+ff(ixO^S,idims)-ff(ixB^S,idims)
      end do
      src(ixO^S)=src(ixO^S)/block%dvolume(ixO^S)
    end if
  end subroutine get_flux_on_cell_face

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine rmhd_add_source(qdt,dtfactor,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source
    use mod_cak_force, only: cak_add_source

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,dtfactor
    double precision, intent(in)    :: wCT(ixI^S,1:nw),wCTprim(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active

    ! TODO local_timestep support is only added for splitting
    ! but not for other nonideal terms such gravity, RC, viscosity,..
    ! it will also only work for divbfix  'linde', which does not require
    ! modification as it does not use dt in the update
    if (.not. qsourcesplit) then
      if(has_equi_pe0) then
        active = .true.
        call add_pe0_divv(qdt,dtfactor,ixI^L,ixO^L,wCTprim,w,x)
      end if
      if(rmhd_hyperbolic_thermal_conduction) then
        call add_hypertc_source(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
      end if
      ! Source for B0 splitting
      if (B0field) then
        active = .true.
        call add_source_B0split(qdt,dtfactor,ixI^L,ixO^L,wCTprim,w,x)
      end if
      ! Sources for resistivity in eqs. for e, B1, B2 and B3
      if (abs(rmhd_eta)>smalldouble)then
        active = .true.
        call add_source_res2(qdt,ixI^L,ixO^L,wCT,w,x)
      end if
      if (rmhd_eta_hyper>0.d0)then
        active = .true.
        call add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
      end if
    end if
    {^NOONED
    if(source_split_divb .eqv. qsourcesplit) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_ct)
        continue ! Do nothing
      case (divb_linde)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_glm)
        active = .true.
        call add_source_glm(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(qdt,ixI^L,ixO^L,wCTprim,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(qdt,ixI^L,ixO^L,wCTprim,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
        call add_source_janhunen(qdt,ixI^L,ixO^L,wCTprim,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
        call add_source_powel(qdt,ixI^L,ixO^L,wCTprim,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
        call add_source_glm(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_multigrid)
        continue ! Do nothing
      case (divb_none)
        ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    end if
    }
    if(rmhd_viscosity) then
      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,rmhd_energy,qsourcesplit,active)
    end if
    if(rmhd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,&
           w,x,gravity_energy,gravity_rhov,qsourcesplit,active)
    end if
    if (rmhd_cak_force) then
      call cak_add_source(qdt,ixI^L,ixO^L,wCT,w,x,rmhd_energy,qsourcesplit,active)
    end if
    !> This is where the radiation force and heating/cooling are added
    call rmhd_add_radiation_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    ! update temperature from new pressure, density, and old ionization degree
    if(rmhd_partial_ionization) then
      if(.not.qsourcesplit) then
        active = .true.
        call rmhd_update_temperature(ixI^L,ixO^L,wCT,w,x)
      end if
    end if
  end subroutine rmhd_add_source

  subroutine rmhd_add_radiation_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods
    use mod_fld
    use mod_afld
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: qsourcesplit
    logical, intent(inout) :: active
    double precision :: cmax(ixI^S)

    select case(rmhd_radiation_formalism)
    case('fld')
      if(fld_diff_scheme .eq. 'mg') call fld_get_diffcoef_central(w, wCT, x, ixI^L, ixO^L)
      !> radiation force
      if(rmhd_radiation_force) call get_fld_rad_force(qdt,ixI^L,ixO^L,wCT,w,x,rmhd_energy,qsourcesplit,active)
      call rmhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'fld_e_interact')
    case('afld')
      if(fld_diff_scheme .eq. 'mg') call afld_get_diffcoef_central(w, wCT, x, ixI^L, ixO^L)
      !> radiation force
      if(rmhd_radiation_force) call get_afld_rad_force(qdt,ixI^L,ixO^L,wCT,w,x,rmhd_energy,qsourcesplit,active)
      call rmhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'fld_e_interact')
      !> photon tiring, heating and cooling
      if(rmhd_energy) then
        if (rmhd_energy_interact) call get_afld_energy_interact(qdt,ixI^L,ixO^L,wCT,w,x,rmhd_energy,qsourcesplit,active)
      endif
    case default
      call mpistop('Radiation formalism unknown')
    end select
  end subroutine rmhd_add_radiation_source

  subroutine add_pe0_divv(qdt,dtfactor,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,dtfactor
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divv(ixI^S)

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(wCT(ixI^S,mom(1:ndir)),ixI^L,ixO^L,divv,3)
      else
        call divvector(wCT(ixI^S,mom(1:ndir)),ixI^L,ixO^L,divv,2)
      end if
    else
     call divvector(wCT(ixI^S,mom(1:ndir)),ixI^L,ixO^L,divv)
    end if
    if(local_timestep) then
      w(ixO^S,e_)=w(ixO^S,e_)-dtfactor*block%dt(ixO^S)*block%equi_vars(ixO^S,equi_pe0_,0)*divv(ixO^S)
    else
      w(ixO^S,e_)=w(ixO^S,e_)-qdt*block%equi_vars(ixO^S,equi_pe0_,0)*divv(ixO^S)
    end if
  end subroutine add_pe0_divv

  subroutine get_tau(ixI^L,ixO^L,w,Te,tau,sigT5)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:nw), intent(in) :: w
    double precision, dimension(ixI^S), intent(in) :: Te
    double precision, dimension(ixI^S), intent(out) :: tau,sigT5
    double precision :: dxmin,taumin
    double precision, dimension(ixI^S) :: sigT7,eint
    integer :: ix^D

    taumin=4.d0
    !> w supposed to be wCTprim here
    if(rmhd_trac) then
      where(Te(ixO^S) .lt. block%wextra(ixO^S,Tcoff_))
        sigT5(ixO^S)=hypertc_kappa*sqrt(block%wextra(ixO^S,Tcoff_)**5)
        sigT7(ixO^S)=sigT5(ixO^S)*block%wextra(ixO^S,Tcoff_)
      else where
        sigT5(ixO^S)=hypertc_kappa*sqrt(Te(ixO^S)**5)
        sigT7(ixO^S)=sigT5(ixO^S)*Te(ixO^S)
      end where
    else
      sigT5(ixO^S)=hypertc_kappa*sqrt(Te(ixO^S)**5)
      sigT7(ixO^S)=sigT5(ixO^S)*Te(ixO^S)
    end if
    eint(ixO^S)=w(ixO^S,p_)/(rmhd_gamma-one)
    tau(ixO^S)=max(taumin*dt,sigT7(ixO^S)/eint(ixO^S)/cmax_global**2)
  end subroutine get_tau

  subroutine add_hypertc_source(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qdt
    double precision, dimension(ixI^S,1:ndim), intent(in) :: x
    double precision, dimension(ixI^S,1:nw), intent(in) :: wCT,wCTprim
    double precision, dimension(ixI^S,1:nw), intent(inout) :: w
    double precision :: invdx
    double precision, dimension(ixI^S) :: Te,tau,sigT,htc_qsrc,Tface,R
    double precision, dimension(ixI^S) :: htc_esrc,Bsum,Bunit
    double precision, dimension(ixI^S,1:ndim) :: Btot
    integer :: idims
    integer :: hxC^L,hxO^L,ixC^L,jxC^L,jxO^L,kxC^L

    call rmhd_get_Rfactor(wCTprim,x,ixI^L,ixI^L,R)
    !Te(ixI^S)=wCTprim(ixI^S,p_)/wCT(ixI^S,rho_)
    Te(ixI^S)=wCTprim(ixI^S,p_)/(R(ixI^S)*w(ixI^S,rho_))
    call get_tau(ixI^L,ixO^L,wCTprim,Te,tau,sigT)
    htc_qsrc=zero
    do idims=1,ndim
      if(B0field) then
        Btot(ixO^S,idims)=wCT(ixO^S,mag(idims))+block%B0(ixO^S,idims,0)
      else
        Btot(ixO^S,idims)=wCT(ixO^S,mag(idims))
      endif
    enddo
    Bsum(ixO^S)=sqrt(sum(Btot(ixO^S,:)**2,dim=ndim+1))+smalldouble
    do idims=1,ndim
      invdx=1.d0/dxlevel(idims)
      ixC^L=ixO^L;
      ixCmin^D=ixOmin^D-kr(idims,^D);ixCmax^D=ixOmax^D;
      jxC^L=ixC^L+kr(idims,^D);
      kxC^L=jxC^L+kr(idims,^D);
      hxC^L=ixC^L-kr(idims,^D);
      hxO^L=ixO^L-kr(idims,^D);
      Tface(ixC^S)=(7.d0*(Te(ixC^S)+Te(jxC^S))-(Te(hxC^S)+Te(kxC^S)))/12.d0
      Bunit(ixO^S)=Btot(ixO^S,idims)/Bsum(ixO^S)
      htc_qsrc(ixO^S)=htc_qsrc(ixO^S)+sigT(ixO^S)*Bunit(ixO^S)*(Tface(ixO^S)-Tface(hxO^S))*invdx
    end do
    htc_qsrc(ixO^S)=(htc_qsrc(ixO^S)+wCT(ixO^S,q_))/tau(ixO^S)
    w(ixO^S,q_)=w(ixO^S,q_)-qdt*htc_qsrc(ixO^S)
  end subroutine add_hypertc_source

  !> Compute the Lorentz force (JxB)
  subroutine get_Lorentz_force(ixI^L,ixO^L,w,JxB)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: JxB(ixI^S,3)
    double precision                :: a(ixI^S,3), b(ixI^S,3)
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision                :: current(ixI^S,7-2*ndir:3)
    integer                         :: idir, idirmin

    b=0.0d0
    if(B0field) then
      do idir = 1, ndir
        b(ixO^S, idir) = w(ixO^S,mag(idir))+block%B0(ixO^S,idir,0)
      end do
    else
      do idir = 1, ndir
        b(ixO^S, idir) = w(ixO^S,mag(idir))
      end do
    end if
    ! store J current in a
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    a=0.0d0
    do idir=7-2*ndir,3
      a(ixO^S,idir)=current(ixO^S,idir)
    end do
    call cross_product(ixI^L,ixO^L,a,b,JxB)
  end subroutine get_Lorentz_force

  !> Compute 1/(1+v_A^2/c^2) for semirelativistic MHD, where v_A is the Alfven
  !> velocity
  subroutine rmhd_gamma2_alfven(ixI^L, ixO^L, w, gamma_A2)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision, intent(out) :: gamma_A2(ixO^S)
    double precision              :: rho(ixI^S)

    ! rmhd_get_rho cannot be used as x is not a param
    if(has_equi_rho0) then
      rho(ixO^S) = w(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,b0i)
    else
      rho(ixO^S) = w(ixO^S,rho_)
    endif
    ! Compute the inverse of 1 + B^2/(rho * c^2)
    gamma_A2(ixO^S) = 1.0d0/(1.0d0+rmhd_mag_en_all(w, ixI^L, ixO^L)/rho(ixO^S)*inv_squared_c)
  end subroutine rmhd_gamma2_alfven

  !> Compute 1/sqrt(1+v_A^2/c^2) for semirelativisitic MHD, where v_A is the
  !> Alfven velocity
  function rmhd_gamma_alfven(w, ixI^L, ixO^L) result(gamma_A)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: gamma_A(ixO^S)

    call rmhd_gamma2_alfven(ixI^L, ixO^L, w, gamma_A)
    gamma_A = sqrt(gamma_A)
  end function rmhd_gamma_alfven

  subroutine rmhd_get_rho(w,x,ixI^L,ixO^L,rho)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: rho(ixI^S)

    if(has_equi_rho0) then
      rho(ixO^S) = w(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,b0i)
    else
      rho(ixO^S) = w(ixO^S,rho_)
    endif
  end subroutine rmhd_get_rho

  !> handle small or negative internal energy
  subroutine rmhd_handle_small_ei(w, x, ixI^L, ixO^L, ie, subname)
    use mod_global_parameters
    use mod_small_values
    integer, intent(in)             :: ixI^L,ixO^L, ie
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    double precision :: rho(ixI^S)
    integer :: idir
    logical :: flag(ixI^S,1:nw)

    flag=.false.
    if(has_equi_pe0) then
      where(w(ixO^S,ie)+block%equi_vars(ixO^S,equi_pe0_,0)*inv_gamma_1<small_e)&
             flag(ixO^S,ie)=.true.
    else
      where(w(ixO^S,ie)<small_e) flag(ixO^S,ie)=.true.
    endif
    if(any(flag(ixO^S,ie))) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_pe0) then
          where(flag(ixO^S,ie)) w(ixO^S,ie)=small_e - &
                  block%equi_vars(ixO^S,equi_pe0_,0)*inv_gamma_1
        else
          where(flag(ixO^S,ie)) w(ixO^S,ie)=small_e
        endif
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag, ie)
      case default
        ! small values error shows primitive variables
        w(ixO^S,e_)=w(ixO^S,e_)*gamma_1
        call rmhd_get_rho(w,x,ixI^L,ixO^L,rho)
        do idir = 1, ndir
           w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))/rho(ixO^S)
        end do
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
  end subroutine rmhd_handle_small_ei

  subroutine rmhd_update_temperature(ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_ionization_degree

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: iz_H(ixO^S),iz_He(ixO^S), pth(ixI^S)

    call ionization_degree_from_temperature(ixI^L,ixO^L,wCT(ixI^S,Te_),iz_H,iz_He)

    call rmhd_get_pthermal(w,x,ixI^L,ixO^L,pth)

    w(ixO^S,Te_)=(2.d0+3.d0*He_abundance)*pth(ixO^S)/(w(ixO^S,rho_)*(1.d0+iz_H(ixO^S)+&
     He_abundance*(iz_He(ixO^S)*(iz_He(ixO^S)+1.d0)+1.d0)))
  end subroutine rmhd_update_temperature

  !> Source terms after split off time-independent magnetic field
  subroutine add_source_B0split(qdt,dtfactor,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, dtfactor,wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: a(ixI^S,3), b(ixI^S,3), axb(ixI^S,3)
    integer :: idir

    a=0.d0
    b=0.d0
    ! for force-free field J0xB0 =0
    if(.not.B0field_forcefree) then
      ! store B0 magnetic field in b
      b(ixO^S,1:ndir)=block%B0(ixO^S,1:ndir,0)
      ! store J0 current in a
      do idir=7-2*ndir,3
        a(ixO^S,idir)=block%J0(ixO^S,idir)
      end do
      call cross_product(ixI^L,ixO^L,a,b,axb)
      if(local_timestep) then
        do idir=1,3
          axb(ixO^S,idir)=axb(ixO^S,idir)*block%dt(ixO^S)*dtfactor
        enddo
      else
        axb(ixO^S,:)=axb(ixO^S,:)*qdt
      endif
      ! add J0xB0 source term in momentum equations
      w(ixO^S,mom(1:ndir))=w(ixO^S,mom(1:ndir))+axb(ixO^S,1:ndir)
    end if
    if(total_energy) then
      a=0.d0
      ! for free-free field -(vxB0) dot J0 =0
      b(ixO^S,:)=wCT(ixO^S,mag(:))
      ! store full magnetic field B0+B1 in b
      if(.not.B0field_forcefree) b(ixO^S,:)=b(ixO^S,:)+block%B0(ixO^S,:,0)
      ! store velocity in a
      a(ixI^S,1:ndir)=wCT(ixI^S,mom(1:ndir))
      ! -E = a x b
      call cross_product(ixI^L,ixO^L,a,b,axb)
      if(local_timestep) then
        do idir=1,3
          axb(ixO^S,idir)=axb(ixO^S,idir)*block%dt(ixO^S)*dtfactor
        enddo
      else
        axb(ixO^S,:)=axb(ixO^S,:)*qdt
      endif
      ! add -(vxB) dot J0 source term in energy equation
      do idir=7-2*ndir,3
        w(ixO^S,e_)=w(ixO^S,e_)-axb(ixO^S,idir)*block%J0(ixO^S,idir)
      end do
    end if
    if (fix_small_values) call rmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_B0')
  end subroutine add_source_B0split

  !> Add resistive source to w within ixO Uses 3 point stencil (1 neighbour) in
  !> each direction, non-conservative. If the fourthorder precompiler flag is
  !> set, uses fourth order central difference for the laplacian. Then the
  !> stencil is 5 (2 neighbours).
  subroutine add_source_res1(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ixA^L,idir,jdir,kdir,idirmin,idim,jxO^L,hxO^L,ix
    integer :: lxO^L, kxO^L
    double precision :: tmp(ixI^S),tmp2(ixI^S)
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S)
    double precision :: gradeta(ixI^S,1:ndim), Bf(ixI^S,1:ndir)

    ! Calculating resistive sources involve one extra layer
    if (rmhd_4th_order) then
      ixA^L=ixO^L^LADD2;
    else
      ixA^L=ixO^L^LADD1;
    end if
    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_res1: Non-conforming input limits")
    ! Calculate current density and idirmin
    call get_current(wCT,ixI^L,ixO^L,idirmin,current)
    if (rmhd_eta>zero)then
       eta(ixA^S)=rmhd_eta
       gradeta(ixO^S,1:ndim)=zero
    else
       call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,current,eta)
       ! assumes that eta is not function of current?
       do idim=1,ndim
          call gradient(eta,ixI^L,ixO^L,idim,tmp)
          gradeta(ixO^S,idim)=tmp(ixO^S)
       end do
    end if
    if(B0field) then
      Bf(ixI^S,1:ndir)=wCT(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,0)
    else
      Bf(ixI^S,1:ndir)=wCT(ixI^S,mag(1:ndir))
    end if
    do idir=1,ndir
       ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
       if (rmhd_4th_order) then
         tmp(ixO^S)=zero
         tmp2(ixI^S)=Bf(ixI^S,idir)
         do idim=1,ndim
            lxO^L=ixO^L+2*kr(idim,^D);
            jxO^L=ixO^L+kr(idim,^D);
            hxO^L=ixO^L-kr(idim,^D);
            kxO^L=ixO^L-2*kr(idim,^D);
            tmp(ixO^S)=tmp(ixO^S)+&
                 (-tmp2(lxO^S)+16.0d0*tmp2(jxO^S)-30.0d0*tmp2(ixO^S)+16.0d0*tmp2(hxO^S)-tmp2(kxO^S)) &
                 /(12.0d0 * dxlevel(idim)**2)
         end do
       else
         tmp(ixO^S)=zero
         tmp2(ixI^S)=Bf(ixI^S,idir)
         do idim=1,ndim
            jxO^L=ixO^L+kr(idim,^D);
            hxO^L=ixO^L-kr(idim,^D);
            tmp(ixO^S)=tmp(ixO^S)+&
                 (tmp2(jxO^S)-2.0d0*tmp2(ixO^S)+tmp2(hxO^S))/dxlevel(idim)**2
         end do
       end if
       ! Multiply by eta
       tmp(ixO^S)=tmp(ixO^S)*eta(ixO^S)
       ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
       if (rmhd_eta<zero)then
          do jdir=1,ndim; do kdir=idirmin,3
             if (lvc(idir,jdir,kdir)/=0)then
                if (lvc(idir,jdir,kdir)==1)then
                   tmp(ixO^S)=tmp(ixO^S)-gradeta(ixO^S,jdir)*current(ixO^S,kdir)
                else
                   tmp(ixO^S)=tmp(ixO^S)+gradeta(ixO^S,jdir)*current(ixO^S,kdir)
                end if
             end if
          end do; end do
       end if
       ! Add sources related to eta*laplB-grad(eta) x J to B and e
       w(ixO^S,mag(idir))=w(ixO^S,mag(idir))+qdt*tmp(ixO^S)
       if(total_energy) then
          w(ixO^S,e_)=w(ixO^S,e_)+qdt*tmp(ixO^S)*Bf(ixO^S,idir)
       end if
    end do ! idir
    if(rmhd_energy) then
      ! de/dt+=eta*J**2
      w(ixO^S,e_)=w(ixO^S,e_)+qdt*eta(ixO^S)*sum(current(ixO^S,:)**2,dim=ndim+1)
    end if
    if (fix_small_values) call rmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res1')
  end subroutine add_source_res1

  !> Add resistive source to w within ixO
  !> Uses 5 point stencil (2 neighbours) in each direction, conservative
  subroutine add_source_res2(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S),curlj(ixI^S,1:3)
    double precision :: tmpvec(ixI^S,1:3),tmp(ixO^S)
    integer :: ixA^L,idir,idirmin,idirmin1

    ixA^L=ixO^L^LADD2;
    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_res2: Non-conforming input limits")
    ixA^L=ixO^L^LADD1;
    ! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
    ! Determine exact value of idirmin while doing the loop.
    call get_current(wCT,ixI^L,ixA^L,idirmin,current)
    tmpvec=zero
    if(rmhd_eta>zero)then
      do idir=idirmin,3
        tmpvec(ixA^S,idir)=current(ixA^S,idir)*rmhd_eta
      end do
    else
      call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,current,eta)
      do idir=idirmin,3
        tmpvec(ixA^S,idir)=current(ixA^S,idir)*eta(ixA^S)
      end do
    end if
    ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
    call curlvector(tmpvec,ixI^L,ixO^L,curlj,idirmin1,1,3)
    if(stagger_grid) then
      if(ndim==2.and.ndir==3) then
        ! if 2.5D
        w(ixO^S,mag(ndir)) = w(ixO^S,mag(ndir))-qdt*curlj(ixO^S,ndir)
      end if
    else
      w(ixO^S,mag(1:ndir)) = w(ixO^S,mag(1:ndir))-qdt*curlj(ixO^S,1:ndir)
    end if
    if(rmhd_energy) then
      if(rmhd_eta>zero)then
        tmp(ixO^S)=qdt*rmhd_eta*sum(current(ixO^S,:)**2,dim=ndim+1)
      else
        tmp(ixO^S)=qdt*eta(ixO^S)*sum(current(ixO^S,:)**2,dim=ndim+1)
      end if
      if(total_energy) then
        ! de/dt= +div(B x Jeta) = eta J^2 - B dot curl(eta J)
        ! de1/dt= eta J^2 - B1 dot curl(eta J)
        w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)-&
        qdt*sum(wCT(ixO^S,mag(1:ndir))*curlj(ixO^S,1:ndir),dim=ndim+1)
      else
        ! add eta*J**2 source term in the internal energy equation
        w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)
      end if
    end if
    if (fix_small_values) call rmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res2')
  end subroutine add_source_res2

  !> Add Hyper-resistive source to w within ixO
  !> Uses 9 point stencil (4 neighbours) in each direction.
  subroutine add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    !.. local ..
    double precision                :: current(ixI^S,7-2*ndir:3)
    double precision                :: tmpvec(ixI^S,1:3),tmpvec2(ixI^S,1:3),tmp(ixI^S),ehyper(ixI^S,1:3)
    integer                         :: ixA^L,idir,jdir,kdir,idirmin,idirmin1

    ixA^L=ixO^L^LADD3;
    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_hyperres: Non-conforming input limits")
    call get_current(wCT,ixI^L,ixA^L,idirmin,current)
    tmpvec(ixA^S,1:ndir)=zero
    do jdir=idirmin,3
       tmpvec(ixA^S,jdir)=current(ixA^S,jdir)
    end do
    ixA^L=ixO^L^LADD2;
    call curlvector(tmpvec,ixI^L,ixA^L,tmpvec2,idirmin1,1,3)
    ixA^L=ixO^L^LADD1;
    tmpvec(ixA^S,1:ndir)=zero
    call curlvector(tmpvec2,ixI^L,ixA^L,tmpvec,idirmin1,1,3)
    ehyper(ixA^S,1:ndir) = - tmpvec(ixA^S,1:ndir)*rmhd_eta_hyper
    ixA^L=ixO^L;
    tmpvec2(ixA^S,1:ndir)=zero
    call curlvector(ehyper,ixI^L,ixA^L,tmpvec2,idirmin1,1,3)
    do idir=1,ndir
      w(ixO^S,mag(idir)) = w(ixO^S,mag(idir))-tmpvec2(ixO^S,idir)*qdt
    end do
    if(total_energy) then
      ! de/dt= +div(B x Ehyper)
      ixA^L=ixO^L^LADD1;
      tmpvec2(ixA^S,1:ndir)=zero
      do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
        tmpvec2(ixA^S,idir) = tmpvec(ixA^S,idir)&
        + lvc(idir,jdir,kdir)*wCT(ixA^S,mag(jdir))*ehyper(ixA^S,kdir)
      end do; end do; end do
      tmp(ixO^S)=zero
      call divvector(tmpvec2,ixI^L,ixO^L,tmp)
      w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)*qdt
    end if
    if (fix_small_values)  call rmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_hyperres')
  end subroutine add_source_hyperres

  subroutine add_source_glm(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-MHD scheme or GLM-MHD scheme
    use mod_global_parameters
    use mod_geometry
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision:: divb(ixI^S), gradPsi(ixI^S), Ba(ixO^S,1:ndir)
    integer          :: idir

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (rmhd_glm_alpha < zero) then
      w(ixO^S,psi_) = abs(rmhd_glm_alpha)*wCT(ixO^S,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab_uniform) then
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*rmhd_glm_alpha/minval(dxlevel(:)))*w(ixO^S,psi_)
      else
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*rmhd_glm_alpha/minval(block%ds(ixO^S,:),dim=ndim+1))*w(ixO^S,psi_)
      end if
    end if
    if(rmhd_glm_extended) then
      if(B0field) then
        Ba(ixO^S,1:ndir)=wCT(ixO^S,mag(1:ndir))+block%B0(ixO^S,1:ndir,0)
      else
        Ba(ixO^S,1:ndir)=wCT(ixO^S,mag(1:ndir))
      end if
      ! gradient of Psi
      if(total_energy) then
        do idir=1,ndim
          select case(typegrad)
          case("central")
            call gradient(wCT(ixI^S,psi_),ixI^L,ixO^L,idir,gradPsi)
          case("limited")
            call gradientL(wCT(ixI^S,psi_),ixI^L,ixO^L,idir,gradPsi)
          end select
          ! e  = e  -qdt (b . grad(Psi))
          w(ixO^S,e_) = w(ixO^S,e_)-qdt*Ba(ixO^S,idir)*gradPsi(ixO^S)
        end do
      end if
      ! We calculate now div B
      call get_divb(wCT,ixI^L,ixO^L,divb,rmhd_divb_nth)
      ! m = m - qdt b div b
      do idir=1,ndir
        w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-qdt*Ba(ixO^S,idir)*divb(ixO^S)
      end do
    end if
    if (fix_small_values) call rmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_glm')
  end subroutine add_source_glm

  !> Add divB related sources to w within ixO corresponding to Powel
  subroutine add_source_powel(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divb(ixI^S), Ba(1:ndir)
    integer                         :: idir, ix^D

    ! calculate div B
    call get_divb(wCT,ixI^L,ixO^L,divb,rmhd_divb_nth)
    if(B0field) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! b = b - qdt v * div b
        ^C&w(ix^D,b^C_)=w(ix^D,b^C_)-qdt*wCT(ix^D,m^C_)*divb(ix^D)\
        ! m = m - qdt b div b
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)-qdt*(wCT(ix^D,b^C_)+block%B0(ix^D,^C,0))*divb(ix^D)\
        if (total_energy) then
          ! e = e - qdt (v . b) * div b
          w(ix^D,e_)=w(ix^D,e_)-qdt*(^C&wCT(ix^D,m^C_)*(wCT(ix^D,b^C_)+block%B0(ix^D,^C,0))+)*divb(ix^D)
        end if
     {end do\}
    else
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! b = b - qdt v * div b
        ^C&w(ix^D,b^C_)=w(ix^D,b^C_)-qdt*wCT(ix^D,m^C_)*divb(ix^D)\
        ! m = m - qdt b div b
        ^C&w(ix^D,m^C_)=w(ix^D,m^C_)-qdt*wCT(ix^D,b^C_)*divb(ix^D)\
        if (total_energy) then
          ! e = e - qdt (v . b) * div b
          w(ix^D,e_)=w(ix^D,e_)-qdt*(^C&wCT(ix^D,m^C_)*wCT(ix^D,b^C_)+)*divb(ix^D)
        end if
     {end do\}
    end if
    if (fix_small_values) call rmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_powel')
  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divb(ixI^S)
    integer                         :: idir, ix^D

    ! calculate div B
    call get_divb(wCT,ixI^L,ixO^L,divb,rmhd_divb_nth)
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! b = b - qdt v * div b
      ^C&w(ix^D,b^C_)=w(ix^D,b^C_)-qdt*wCT(ix^D,m^C_)*divb(ix^D)\
    {end do\}
    if (fix_small_values) call rmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_janhunen')
  end subroutine add_source_janhunen

  subroutine add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add Linde's divB related sources to wnew within ixO
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: divb(ixI^S),graddivb(ixI^S)
    integer :: idim, idir, ixp^L, i^D, iside
    logical, dimension(-1:1^D&) :: leveljump

    ! Calculate div B
    ixp^L=ixO^L^LADD1;
    call get_divb(wCT,ixI^L,ixp^L,divb,rmhd_divb_nth)
    ! for AMR stability, retreat one cell layer from the boarders of level jump
    {do i^DB=-1,1\}
      if(i^D==0|.and.) cycle
      if(neighbor_type(i^D,block%igrid)==2 .or. neighbor_type(i^D,block%igrid)==4) then
        leveljump(i^D)=.true.
      else
        leveljump(i^D)=.false.
      end if
    {end do\}
    ixp^L=ixO^L;
    do idim=1,ndim
      select case(idim)
       {case(^D)
          do iside=1,2
            i^DD=kr(^DD,^D)*(2*iside-3);
            if (leveljump(i^DD)) then
              if (iside==1) then
                ixpmin^D=ixOmin^D-i^D
              else
                ixpmax^D=ixOmax^D-i^D
              end if
            end if
          end do
       \}
      end select
    end do
    ! Add Linde's diffusive terms
    do idim=1,ndim
       ! Calculate grad_idim(divb)
       select case(typegrad)
       case("central")
         call gradient(divb,ixI^L,ixp^L,idim,graddivb)
       case("limited")
         call gradientL(divb,ixI^L,ixp^L,idim,graddivb)
       end select
       ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab_uniform) then
          graddivb(ixp^S)=graddivb(ixp^S)*divbdiff/(^D&1.0d0/dxlevel(^D)**2+)
       else
          graddivb(ixp^S)=graddivb(ixp^S)*divbdiff &
                          /(^D&1.0d0/block%ds(ixp^S,^D)**2+)
       end if
       w(ixp^S,mag(idim))=w(ixp^S,mag(idim))+graddivb(ixp^S)

       if (typedivbdiff=='all' .and. total_energy) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixp^S,e_)=w(ixp^S,e_)+wCT(ixp^S,mag(idim))*graddivb(ixp^S)
       end if
    end do
    if (fix_small_values) call rmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_linde')
  end subroutine add_source_linde

  !> get dimensionless div B = |divB| * volume / area / |B|
  subroutine get_normalized_divb(w,ixI^L,ixO^L,divb)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision                   :: divb(ixI^S), dsurface(ixI^S)
    double precision :: invB(ixO^S)
    integer :: ixA^L,idims

    call get_divb(w,ixI^L,ixO^L,divb)
    invB(ixO^S)=sqrt(rmhd_mag_en_all(w,ixI^L,ixO^L))
    where(invB(ixO^S)/=0.d0)
      invB(ixO^S)=1.d0/invB(ixO^S)
    end where
    if(slab_uniform) then
      divb(ixO^S)=0.5d0*abs(divb(ixO^S))*invB(ixO^S)/sum(1.d0/dxlevel(:))
    else
      ixAmin^D=ixOmin^D-1;
      ixAmax^D=ixOmax^D-1;
      dsurface(ixO^S)= sum(block%surfaceC(ixO^S,:),dim=ndim+1)
      do idims=1,ndim
        ixA^L=ixO^L-kr(idims,^D);
        dsurface(ixO^S)=dsurface(ixO^S)+block%surfaceC(ixA^S,idims)
      end do
      divb(ixO^S)=abs(divb(ixO^S))*invB(ixO^S)*&
      block%dvolume(ixO^S)/dsurface(ixO^S)
    end if
  end subroutine get_normalized_divb

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixI^L,ixO^L,idirmin,current)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)  :: ixO^L, ixI^L
    double precision, intent(in) :: w(ixI^S,1:nw)
    integer, intent(out) :: idirmin
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3)
    integer :: idir, idirmin0

    idirmin0 = 7-2*ndir
    call curlvector(w(ixI^S,mag(1:ndir)),ixI^L,ixO^L,current,idirmin,idirmin0,ndir)
    if(B0field) current(ixO^S,idirmin0:3)=current(ixO^S,idirmin0:3)+&
        block%J0(ixO^S,idirmin0:3)
  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine rmhd_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt
    use mod_cak_force, only: cak_get_dt
    use mod_fld, only: fld_radforce_get_dt
    use mod_afld, only: afld_radforce_get_dt
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx^D
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision              :: dxarr(ndim)
    double precision              :: current(ixI^S,7-2*ndir:3),eta(ixI^S)
    integer                       :: idirmin,idim

    dtnew = bigdouble

    if (.not. dt_c) then
      ^D&dxarr(^D)=dx^D;
      if (rmhd_eta>zero)then
         dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/rmhd_eta
      else if (rmhd_eta<zero)then
         call get_current(w,ixI^L,ixO^L,idirmin,current)
         call usr_special_resistivity(w,ixI^L,ixO^L,idirmin,x,current,eta)
         dtnew=bigdouble
         do idim=1,ndim
           if(slab_uniform) then
             dtnew=min(dtnew,&
                  dtdiffpar/(smalldouble+maxval(eta(ixO^S)/dxarr(idim)**2)))
           else
             dtnew=min(dtnew,&
                  dtdiffpar/(smalldouble+maxval(eta(ixO^S)/block%ds(ixO^S,idim)**2)))
           end if
         end do
      end if
      if(rmhd_eta_hyper>zero) then
        if(slab_uniform) then
          dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/rmhd_eta_hyper,dtnew)
        else
          dtnew=min(dtdiffpar*minval(block%ds(ixO^S,1:ndim))**4/rmhd_eta_hyper,dtnew)
        end if
      end if
      if(rmhd_radiation_force) then
        select case(rmhd_radiation_formalism)
        case('fld')
          call fld_radforce_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
        case('afld')
          call afld_radforce_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
        case default
          call mpistop('Radiation formalism unknown')
        end select
      endif
      if(rmhd_viscosity) then
        call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
      end if
      if(rmhd_gravity) then
        call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
      end if
      if (rmhd_cak_force) then
        call cak_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
      end if
    else
      {^IFONED dtnew = dx1*unit_velocity/const_c}
      {^NOONED dtnew = min(dx^D*unit_velocity/const_c)}
    endif
  end subroutine rmhd_get_dt

  ! Add geometrical source terms to w
  subroutine rmhd_add_source_geom(qdt,dtfactor,ixI^L,ixO^L,wCT,wprim,w,x)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor,x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw),wprim(ixI^S,1:nw),w(ixI^S,1:nw)
    double precision :: tmp,tmp1,invr,cot
    integer          :: ix^D
    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_
    select case (coordinate)
    case (cylindrical)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! include dt in invr, invr is always used with qdt
        if(local_timestep) then
          invr=block%dt(ix^D) * dtfactor/x(ix^D,1)
        else
          invr=qdt/x(ix^D,1)
        end if
        if(rmhd_energy) then
          tmp=wprim(ix^D,p_)+half*(^C&wprim(ix^D,b^C_)**2+)
        else
          tmp=rmhd_adiab*wprim(ix^D,rho_)**rmhd_gamma+half*(^C&wprim(ix^D,b^C_)**2+)
        end if
        if(phi_>0) then
          w(ix^D,mr_)=w(ix^D,mr_)+invr*(tmp-&
                    wprim(ix^D,bphi_)**2+wprim(ix^D,mphi_)*wCT(ix^D,mphi_))
          w(ix^D,mphi_)=w(ix^D,mphi_)+invr*(&
                   -wCT(ix^D,mphi_)*wprim(ix^D,mr_) &
                   +wprim(ix^D,bphi_)*wprim(ix^D,br_))
          if(.not.stagger_grid) then
            w(ix^D,bphi_)=w(ix^D,bphi_)+invr*&
                     (wprim(ix^D,bphi_)*wprim(ix^D,mr_) &
                     -wprim(ix^D,br_)*wprim(ix^D,mphi_))
          end if
        else
          w(ix^D,mr_)=w(ix^D,mr_)+invr*tmp
        end if
        if(rmhd_glm) w(ix^D,br_)=w(ix^D,br_)+wprim(ix^D,psi_)*invr
     {end do\}
    case (spherical)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! include dt in invr, invr is always used with qdt
        if(local_timestep) then
          invr=block%dt(ix^D) * dtfactor/x(ix^D,1)
        else
          invr=qdt/x(ix^D,1)
        end if
        if(rmhd_energy) then
          tmp1=wprim(ix^D,p_)+half*(^C&wprim(ix^D,b^C_)**2+)
        else
          tmp1=rmhd_adiab*wprim(ix^D,rho_)**rmhd_gamma+half*(^C&wprim(ix^D,b^C_)**2+)
        end if
        ! m1
        {^IFONEC
        w(ix^D,mom(1))=w(ix^D,mom(1))+two*tmp1*invr
        }
        {^NOONEC
        w(ix^D,mom(1))=w(ix^D,mom(1))+invr*&
         (two*tmp1+(^CE&wprim(ix^D,m^CE_)*wCT(ix^D,m^CE_)-wprim(ix^D,b^CE_)**2+))
        }
        ! b1
        if(rmhd_glm) then
          w(ix^D,mag(1))=w(ix^D,mag(1))+invr*2.0d0*wprim(ix^D,psi_)
        end if
        {^IFONED
        cot=0.d0
        }
        {^NOONED
        cot=1.d0/tan(x(ix^D,2))
        }
        {^IFTWOC
        ! m2
        w(ix^D,mom(2))=w(ix^D,mom(2))+invr*(tmp1*cot-wprim(ix^D,m1_)*wCT(ix^D,m2_)&
          +wprim(ix^D,b1_)*wprim(ix^D,b2_))
        ! b2
        if(.not.stagger_grid) then
          tmp=wprim(ix^D,m1_)*wprim(ix^D,b2_)-wprim(ix^D,m2_)*wprim(ix^D,b1_)
          if(rmhd_glm) then
            tmp=tmp+wprim(ix^D,psi_)*cot
          end if
          w(ix^D,mag(2))=w(ix^D,mag(2))+tmp*invr
        end if
        }
        {^IFTHREEC
        ! m2
        w(ix^D,mom(2))=w(ix^D,mom(2))+invr*(tmp1*cot-wprim(ix^D,m1_)*wCT(ix^D,m2_)&
          +wprim(ix^D,b1_)*wprim(ix^D,b2_)&
          +(wprim(ix^D,m3_)*wCT(ix^D,m3_)-wprim(ix^D,b3_)**2)*cot)
        ! b2
        if(.not.stagger_grid) then
          tmp=wprim(ix^D,m1_)*wprim(ix^D,b2_)-wprim(ix^D,m2_)*wprim(ix^D,b1_)
          if(rmhd_glm) then
            tmp=tmp+wprim(ix^D,psi_)*cot
          end if
          w(ix^D,mag(2))=w(ix^D,mag(2))+tmp*invr
        end if
        ! m3
        w(ix^D,mom(3))=w(ix^D,mom(3))-invr*&
             (wprim(ix^D,m3_)*wCT(ix^D,m1_) &
             -wprim(ix^D,b3_)*wprim(ix^D,b1_) &
            +(wprim(ix^D,m2_)*wCT(ix^D,m3_) &
             -wprim(ix^D,b2_)*wprim(ix^D,b3_))*cot)
        ! b3
        if(.not.stagger_grid) then
          w(ix^D,mag(3))=w(ix^D,mag(3))+invr*&
             (wprim(ix^D,m1_)*wprim(ix^D,b3_) &
             -wprim(ix^D,m3_)*wprim(ix^D,b1_) &
            -(wprim(ix^D,m3_)*wprim(ix^D,b2_) &
             -wprim(ix^D,m2_)*wprim(ix^D,b3_))*cot)
        end if
        }
     {end do\}
    end select
  end subroutine rmhd_add_source_geom

  ! Add geometrical source terms to w
  subroutine rmhd_add_source_geom_split(qdt,dtfactor, ixI^L,ixO^L,wCT,wprim,w,x)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), wprim(ixI^S,1:nw),w(ixI^S,1:nw)
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),invrho(ixO^S),invr(ixO^S)
    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_
    if(has_equi_rho0) then
      invrho(ixO^S) = 1d0/(wCT(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,b0i))
    else
      invrho(ixO^S) = 1d0/wCT(ixO^S,rho_)
    end if
    ! include dt in invr, invr is always used with qdt
    if(local_timestep) then
      invr(ixO^S) = block%dt(ixO^S) * dtfactor/x(ixO^S,1)
    else
      invr(ixO^S) = qdt/x(ixO^S,1)
    end if

    select case (coordinate)
    case (cylindrical)
      call rmhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp)
      if(phi_>0) then
        w(ixO^S,mr_)=w(ixO^S,mr_)+invr(ixO^S)*(tmp(ixO^S)-&
                  wCT(ixO^S,bphi_)**2+wCT(ixO^S,mphi_)**2*invrho(ixO^S))
        w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt*invr(ixO^S)*(&
                 -wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)*invrho(ixO^S) &
                 +wCT(ixO^S,bphi_)*wCT(ixO^S,br_))
        if(.not.stagger_grid) then
          w(ixO^S,bphi_)=w(ixO^S,bphi_)+invr(ixO^S)*&
                   (wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
                   -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
                   *invrho(ixO^S)
        end if
      else
        w(ixO^S,mr_)=w(ixO^S,mr_)+invr(ixO^S)*tmp(ixO^S)
      end if
      if(rmhd_glm) w(ixO^S,br_)=w(ixO^S,br_)+wCT(ixO^S,psi_)*invr(ixO^S)
    case (spherical)
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       call rmhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp1)
       tmp(ixO^S)=tmp1(ixO^S)
       if(B0field) then
         tmp2(ixO^S)=sum(block%B0(ixO^S,:,0)*wCT(ixO^S,mag(:)),dim=ndim+1)
         tmp(ixO^S)=tmp(ixO^S)+tmp2(ixO^S)
       end if
       ! m1
       tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
                  *(block%surfaceC(ixO^S,1)-block%surfaceC(h1x^S,1))/block%dvolume(ixO^S)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(idir))**2*invrho(ixO^S)-wCT(ixO^S,mag(idir))**2
           if(B0field) tmp(ixO^S)=tmp(ixO^S)-2.0d0*block%B0(ixO^S,idir,0)*wCT(ixO^S,mag(idir))
         end do
       end if
       w(ixO^S,mom(1))=w(ixO^S,mom(1))+tmp(ixO^S)*invr(ixO^S)
       ! b1
       if(rmhd_glm) then
         w(ixO^S,mag(1))=w(ixO^S,mag(1))+invr(ixO^S)*2.0d0*wCT(ixO^S,psi_)
       end if
       {^NOONED
       ! m2
       tmp(ixO^S)=tmp1(ixO^S)
       if(B0field) then
         tmp(ixO^S)=tmp(ixO^S)+tmp2(ixO^S)
       end if
       if(local_timestep) then
         tmp1(ixO^S) = block%dt(ixO^S) * tmp(ixO^S)
       else
         tmp1(ixO^S) = qdt * tmp(ixO^S)
       endif  
       ! This will make hydrostatic p=const an exact solution
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+tmp1(ixO^S) &
            *(block%surfaceC(ixO^S,2)-block%surfaceC(h2x^S,2)) &
            /block%dvolume(ixO^S)
       tmp(ixO^S)=-(wCT(ixO^S,mom(1))*wCT(ixO^S,mom(2))*invrho(ixO^S) &
            -wCT(ixO^S,mag(1))*wCT(ixO^S,mag(2)))
       if (B0field) then
          tmp(ixO^S)=tmp(ixO^S)+block%B0(ixO^S,1,0)*wCT(ixO^S,mag(2)) &
               +wCT(ixO^S,mag(1))*block%B0(ixO^S,2,0)
       end if
       if(ndir==3) then
         tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(3))**2*invrho(ixO^S) &
              -wCT(ixO^S,mag(3))**2)*dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)-2.0d0*block%B0(ixO^S,3,0)*wCT(ixO^S,mag(3))&
                 *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
         end if
       end if
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+tmp(ixO^S)*invr(ixO^S)
       ! b2
       if(.not.stagger_grid) then
         tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(2)) &
              -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(1)))*invrho(ixO^S)
         if(B0field) then
           tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(1))*block%B0(ixO^S,2,0) &
                -wCT(ixO^S,mom(2))*block%B0(ixO^S,1,0))*invrho(ixO^S)
         end if
         if(rmhd_glm) then
           tmp(ixO^S)=tmp(ixO^S) &
                + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
         end if
         w(ixO^S,mag(2))=w(ixO^S,mag(2))+tmp(ixO^S)*invr(ixO^S)
       end if
       }
       if(ndir==3) then
         ! m3
         tmp(ixO^S)=-(wCT(ixO^S,mom(3))*wCT(ixO^S,mom(1))*invrho(ixO^S) &
              -wCT(ixO^S,mag(3))*wCT(ixO^S,mag(1))) {^NOONED &
              -(wCT(ixO^S,mom(2))*wCT(ixO^S,mom(3))*invrho(ixO^S) &
              -wCT(ixO^S,mag(2))*wCT(ixO^S,mag(3))) &
              *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)+block%B0(ixO^S,1,0)*wCT(ixO^S,mag(3)) &
                 +wCT(ixO^S,mag(1))*block%B0(ixO^S,3,0) {^NOONED &
                 +(block%B0(ixO^S,2,0)*wCT(ixO^S,mag(3)) &
                 +wCT(ixO^S,mag(2))*block%B0(ixO^S,3,0)) &
                 *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
         end if
         w(ixO^S,mom(3))=w(ixO^S,mom(3))+tmp(ixO^S)*invr(ixO^S)
         ! b3
         if(.not.stagger_grid) then
           tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(3)) &
                -wCT(ixO^S,mom(3))*wCT(ixO^S,mag(1)))*invrho(ixO^S) {^NOONED &
                -(wCT(ixO^S,mom(3))*wCT(ixO^S,mag(2)) &
                -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(3)))*dcos(x(ixO^S,2)) &
                *invrho(ixO^S)/dsin(x(ixO^S,2)) }
           if (B0field) then
              tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(1))*block%B0(ixO^S,3,0) &
                   -wCT(ixO^S,mom(3))*block%B0(ixO^S,1,0))*invrho(ixO^S){^NOONED &
                   -(wCT(ixO^S,mom(3))*block%B0(ixO^S,2,0) &
                   -wCT(ixO^S,mom(2))*block%B0(ixO^S,3,0))*dcos(x(ixO^S,2)) &
                   *invrho(ixO^S)/dsin(x(ixO^S,2)) }
           end if
           w(ixO^S,mag(3))=w(ixO^S,mag(3))+tmp(ixO^S)*invr(ixO^S)
         end if
       end if
    end select
  end subroutine rmhd_add_source_geom_split

  !> Compute 2 times total magnetic energy
  function rmhd_mag_en_all(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    if (B0field) then
      mge = sum((w(ixO^S, mag(:))+block%B0(ixO^S,:,b0i))**2, dim=ndim+1)
    else
      mge = sum(w(ixO^S, mag(:))**2, dim=ndim+1)
    end if
  end function rmhd_mag_en_all

  subroutine rmhd_modify_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    double precision, intent(inout) :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw)
    type(state)                     :: s
    double precision                :: dB(ixO^S), dPsi(ixO^S)
    integer :: ix^D

    if(stagger_grid) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        wLC(ix^D,mag(idir))=s%ws(ix^D,idir)
        wRC(ix^D,mag(idir))=s%ws(ix^D,idir)
        wLp(ix^D,mag(idir))=s%ws(ix^D,idir)
        wRp(ix^D,mag(idir))=s%ws(ix^D,idir)
     {end do\}
    else
      ! Solve the Riemann problem for the linear 2x2 system for normal
      ! B-field and GLM_Psi according to Dedner 2002:
      ! This implements eq. (42) in Dedner et al. 2002 JcP 175
      ! Gives the Riemann solution on the interface
      ! for the normal B component and Psi in the GLM-MHD system.
      ! 23/04/2013 Oliver Porth
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        dB(ix^D)=wRp(ix^D,mag(idir))-wLp(ix^D,mag(idir))
        dPsi(ix^D)=wRp(ix^D,psi_)-wLp(ix^D,psi_)
        wLp(ix^D,mag(idir))=half*(wRp(ix^D,mag(idir))+wLp(ix^D,mag(idir))-dPsi(ix^D)/cmax_global)
        wLp(ix^D,psi_)=half*(wRp(ix^D,psi_)+wLp(ix^D,psi_)-dB(ix^D)*cmax_global)
        wRp(ix^D,mag(idir))=wLp(ix^D,mag(idir))
        wRp(ix^D,psi_)=wLp(ix^D,psi_)
        if(total_energy) then
          wRC(ix^D,e_)=wRC(ix^D,e_)-half*wRC(ix^D,mag(idir))**2
          wLC(ix^D,e_)=wLC(ix^D,e_)-half*wLC(ix^D,mag(idir))**2
        end if
        wRC(ix^D,mag(idir))=wLp(ix^D,mag(idir))
        wRC(ix^D,psi_)=wLp(ix^D,psi_)
        wLC(ix^D,mag(idir))=wLp(ix^D,mag(idir))
        wLC(ix^D,psi_)=wLp(ix^D,psi_)
        ! modify total energy according to the change of magnetic field
        if(total_energy) then
          wRC(ix^D,e_)=wRC(ix^D,e_)+half*wRC(ix^D,mag(idir))**2
          wLC(ix^D,e_)=wLC(ix^D,e_)+half*wLC(ix^D,mag(idir))**2
        end if
     {end do\}
    end if
    if(associated(usr_set_wLR)) call usr_set_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)
  end subroutine rmhd_modify_wLR

  subroutine rmhd_boundary_adjust(igrid,psb)
    use mod_global_parameters
    integer, intent(in) :: igrid
    type(state), target :: psb(max_blocks)
    integer :: iB, idims, iside, ixO^L, i^D

    block=>ps(igrid)
    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
    do idims=1,ndim
       ! to avoid using as yet unknown corner info in more than 1D, we
       ! fill only interior mesh ranges of the ghost cell ranges at first,
       ! and progressively enlarge the ranges to include corners later
       do iside=1,2
          i^D=kr(^D,idims)*(2*iside-3);
          if (neighbor_type(i^D,igrid)/=1) cycle
          iB=(idims-1)*2+iside
          if(.not.boundary_divbfix(iB)) cycle
          if(any(typeboundary(:,iB)==bc_special)) then
            ! MF nonlinear force-free B field extrapolation and data driven
            ! require normal B of the first ghost cell layer to be untouched by
            ! fixdivB=0 process, set boundary_divbfix_skip(iB)=1 in par file
            select case (idims)
            {case (^D)
               if (iside==2) then
                  ! maximal boundary
                  ixOmin^DD=ixGhi^D+1-nghostcells+boundary_divbfix_skip(2*^D)^D%ixOmin^DD=ixGlo^DD;
                  ixOmax^DD=ixGhi^DD;
               else
                  ! minimal boundary
                  ixOmin^DD=ixGlo^DD;
                  ixOmax^DD=ixGlo^D-1+nghostcells-boundary_divbfix_skip(2*^D-1)^D%ixOmax^DD=ixGhi^DD;
               end if \}
            end select
            call fixdivB_boundary(ixG^LL,ixO^L,psb(igrid)%w,psb(igrid)%x,iB)
          end if
       end do
    end do
  end subroutine rmhd_boundary_adjust

  subroutine fixdivB_boundary(ixG^L,ixO^L,w,x,iB)
    use mod_global_parameters
    integer, intent(in) :: ixG^L,ixO^L,iB
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision :: dx1x2,dx1x3,dx2x1,dx2x3,dx3x1,dx3x2
    integer :: ix^D,ixF^L

    select case(iB)
     case(1)
       ! 2nd order CD for divB=0 to set normal B component better
       if(total_energy) call rmhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=w(ix1+1,ixFmin2:ixFmax2,mag(1)) &
            +dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-&
                    w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=( (w(ix1+1,ixFmin2:ixFmax2,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC(ix1,ixFmin2:ixFmax2,1)&
           +(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,ixFmin2:ixFmax2,mag(2)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,2)&
           -(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))*&
             block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,2) )&
            /block%surfaceC(ix1-1,ixFmin2:ixFmax2,1)-w(ix1,ixFmin2:ixFmax2,mag(1))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         dx1x3=dxlevel(1)/dxlevel(3)
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
                     w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)) &
             +dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))-&
                     w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2))) &
             +dx1x3*(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
          ( (w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1)&
           +(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,2)&
           -(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,2)&
           +(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,3)&
           -(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,3) )&
            /block%surfaceC(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1)-&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))
         end do
       end if
       }
       if(total_energy) call rmhd_to_conserved(ixG^L,ixO^L,w,x)
     case(2)
       if(total_energy) call rmhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=w(ix1-1,ixFmin2:ixFmax2,mag(1)) &
            -dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-&
                    w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=( (w(ix1-1,ixFmin2:ixFmax2,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC(ix1-1,ixFmin2:ixFmax2,1)&
           -(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,ixFmin2:ixFmax2,mag(2)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,2)&
           +(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))*&
             block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,2) )&
            /block%surfaceC(ix1,ixFmin2:ixFmax2,1)-w(ix1,ixFmin2:ixFmax2,mag(1))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         dx1x3=dxlevel(1)/dxlevel(3)
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
                     w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)) &
             -dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))-&
                     w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2))) &
             -dx1x3*(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
          ( (w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1)&
           -(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,2)&
           +(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,2)&
           -(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,3)&
           +(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,3) )&
            /block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1)-&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))
         end do
       end if
       }
       if(total_energy) call rmhd_to_conserved(ixG^L,ixO^L,w,x)
     case(3)
       if(total_energy) call rmhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=w(ixFmin1:ixFmax1,ix2+1,mag(2)) &
            +dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))-&
                    w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         enddo
       else
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=( (w(ixFmin1:ixFmax1,ix2+1,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2,2)&
           +(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,1)&
           -(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,1) )&
            /block%surfaceC(ixFmin1:ixFmax1,ix2-1,2)-w(ixFmin1:ixFmax1,ix2,mag(2))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         dx2x3=dxlevel(2)/dxlevel(3)
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,mag(2))=w(ixFmin1:ixFmax1,&
             ix2+1,ixFmin3:ixFmax3,mag(2)) &
             +dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1))) &
             +dx2x3*(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,mag(2))=&
          ( (w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,2)&
           +(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,1)&
           -(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,1)&
           +(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,3)&
           -(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,3) )&
            /block%surfaceC(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,2)-&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2))
         end do
       end if
       }
       if(total_energy) call rmhd_to_conserved(ixG^L,ixO^L,w,x)
     case(4)
       if(total_energy) call rmhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=w(ixFmin1:ixFmax1,ix2-1,mag(2)) &
            -dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))-&
                    w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         end do
       else
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=( (w(ixFmin1:ixFmax1,ix2-1,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2-1,2)&
           -(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,1)&
           +(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,1) )&
            /block%surfaceC(ixFmin1:ixFmax1,ix2,2)-w(ixFmin1:ixFmax1,ix2,mag(2))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         dx2x3=dxlevel(2)/dxlevel(3)
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,mag(2))=w(ixFmin1:ixFmax1,&
             ix2-1,ixFmin3:ixFmax3,mag(2)) &
             -dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1))) &
             -dx2x3*(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,mag(2))=&
          ( (w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,2)&
           -(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,1)&
           +(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,1)&
           -(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,3)&
           +(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,3) )&
            /block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,2)-&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2))
         end do
       end if
       }
       if(total_energy) call rmhd_to_conserved(ixG^L,ixO^L,w,x)
     {^IFTHREED
     case(5)
       if(total_energy) call rmhd_to_primitive(ixG^L,ixO^L,w,x)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3+1
       if(slab_uniform) then
         dx3x1=dxlevel(3)/dxlevel(1)
         dx3x2=dxlevel(3)/dxlevel(2)
         do ix3=ixFmax3,ixFmin3,-1
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,mag(3))=w(ixFmin1:ixFmax1,&
             ixFmin2:ixFmax2,ix3+1,mag(3)) &
             +dx3x1*(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1))) &
             +dx3x2*(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))-&
                     w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))
         end do
       else
         do ix3=ixFmax3,ixFmin3,-1
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,mag(3))=&
          ( (w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,mag(3))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,3)&
           +(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,1)&
           -(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,1)&
           +(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,2)&
           -(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,2) )&
            /block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,3)-&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3))
         end do
       end if
       if(total_energy) call rmhd_to_conserved(ixG^L,ixO^L,w,x)
     case(6)
       if(total_energy) call rmhd_to_primitive(ixG^L,ixO^L,w,x)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3-1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx3x1=dxlevel(3)/dxlevel(1)
         dx3x2=dxlevel(3)/dxlevel(2)
         do ix3=ixFmin3,ixFmax3
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,mag(3))=w(ixFmin1:ixFmax1,&
             ixFmin2:ixFmax2,ix3-1,mag(3)) &
             -dx3x1*(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1))) &
             -dx3x2*(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))-&
                     w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))
         end do
       else
         do ix3=ixFmin3,ixFmax3
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,mag(3))=&
          ( (w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,mag(3))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,3)&
           -(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,1)&
           +(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,1)&
           -(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,2)&
           +(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))*&
             block%surfaceC(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,2) )&
            /block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,3)-&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3))
         end do
       end if
       if(total_energy) call rmhd_to_conserved(ixG^L,ixO^L,w,x)
     }
     case default
       call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine fixdivB_boundary

  {^NOONED
  subroutine rmhd_clean_divb_multigrid(qdt, qt, active)
    use mod_forest
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_geometry
    double precision, intent(in) :: qdt    !< Current time step
    double precision, intent(in) :: qt     !< Current time
    logical, intent(inout)       :: active !< Output if the source is active
    integer                      :: id
    integer, parameter           :: max_its      = 50
    double precision             :: residual_it(max_its), max_divb
    double precision             :: tmp(ixG^T), grad(ixG^T, ndim)
    double precision             :: res
    double precision, parameter  :: max_residual = 1d-3
    double precision, parameter  :: residual_reduction = 1d-10
    integer                      :: iigrid, igrid
    integer                      :: n, nc, lvl, ix^L, ixC^L, idim
    type(tree_node), pointer     :: pnode

    mg%operator_type = mg_laplacian
    ! Set boundary conditions
    do n = 1, 2*ndim
       idim = (n+1)/2
       select case (typeboundary(mag(idim), n))
       case (bc_symm)
          ! d/dx B = 0, take phi = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_asymm)
          ! B = 0, so grad(phi) = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_cont)
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_special)
          ! Assume Dirichlet boundary conditions, derivative zero
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_periodic)
          ! Nothing to do here
       case default
          write(*,*) "rmhd_clean_divb_multigrid warning: unknown boundary type"
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       end select
    end do

    ix^L=ixM^LL^LADD1;
    max_divb = 0.0d0
    ! Store divergence of B as right-hand side
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       call get_divb(ps(igrid)%w(ixG^T, 1:nw), ixG^LL, ixM^LL, tmp, &
            rmhd_divb_nth)
       mg%boxes(id)%cc({1:nc}, mg_irhs) = tmp(ixM^T)
       max_divb = max(max_divb, maxval(abs(tmp(ixM^T))))
    end do

    ! Solve laplacian(phi) = divB
    if(stagger_grid) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, max_divb, 1, MPI_DOUBLE_PRECISION, &
           MPI_MAX, icomm, ierrmpi)

      if (mype == 0) print *, "Performing multigrid divB cleaning"
      if (mype == 0) print *, "iteration vs residual"
      ! Solve laplacian(phi) = divB
      do n = 1, max_its
         call mg_fas_fmg(mg, n>1, max_res=residual_it(n))
         if (mype == 0) write(*, "(I4,E11.3)") n, residual_it(n)
         if (residual_it(n) < residual_reduction * max_divb) exit
      end do
      if (mype == 0 .and. n > max_its) then
         print *, "divb_multigrid warning: not fully converged"
         print *, "current amplitude of divb: ", residual_it(max_its)
         print *, "multigrid smallest grid: ", &
              mg%domain_size_lvl(:, mg%lowest_lvl)
         print *, "note: smallest grid ideally has <= 8 cells"
         print *, "multigrid dx/dy/dz ratio: ", mg%dr(:, 1)/mg%dr(1, 1)
         print *, "note: dx/dy/dz should be similar"
      end if
    else
      do n = 1, max_its
         call mg_fas_vcycle(mg, max_res=res)
         if (res < max_residual) exit
      end do
      if (res > max_residual) call mpistop("divb_multigrid: no convergence")
    end if

    ! Correct the magnetic field
    do iigrid = 1, igridstail
       igrid = igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       ! Compute the gradient of phi
       tmp(ix^S) = mg%boxes(id)%cc({:,}, mg_iphi)
       if(stagger_grid) then
         do idim =1, ndim
           ixCmin^D=ixMlo^D-kr(idim,^D);
           ixCmax^D=ixMhi^D;
           call gradientF(tmp,ps(igrid)%x,ixG^LL,ixC^L,idim,grad(ixG^T,idim))
           ! Apply the correction B* = B - gradient(phi)
           ps(igrid)%ws(ixC^S,idim)=ps(igrid)%ws(ixC^S,idim)-grad(ixC^S,idim)
         end do
         ! store cell-center magnetic energy
         tmp(ixM^T) = sum(ps(igrid)%w(ixM^T, mag(1:ndim))**2, dim=ndim+1)
         ! change cell-center magnetic field
         call rmhd_face_to_center(ixM^LL,ps(igrid))
       else
         do idim = 1, ndim
            call gradient(tmp,ixG^LL,ixM^LL,idim,grad(ixG^T, idim))
         end do
         ! store cell-center magnetic energy
         tmp(ixM^T) = sum(ps(igrid)%w(ixM^T, mag(1:ndim))**2, dim=ndim+1)
         ! Apply the correction B* = B - gradient(phi)
         ps(igrid)%w(ixM^T, mag(1:ndim)) = &
              ps(igrid)%w(ixM^T, mag(1:ndim)) - grad(ixM^T, :)
       end if
       if(total_energy) then
         ! Determine magnetic energy difference
         tmp(ixM^T) = 0.5_dp * (sum(ps(igrid)%w(ixM^T, &
              mag(1:ndim))**2, dim=ndim+1) - tmp(ixM^T))
         ! Keep thermal pressure the same
         ps(igrid)%w(ixM^T, e_) = ps(igrid)%w(ixM^T, e_) + tmp(ixM^T)
       end if
    end do
    active = .true.
  end subroutine rmhd_clean_divb_multigrid
  }

  subroutine rmhd_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s,vcts)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wprim(ixI^S,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,sdim:3)

    select case(type_ct)
    case('average')
      call update_faces_average(ixI^L,ixO^L,qt,qdt,fC,fE,sCT,s)
    case('uct_contact')
      call update_faces_contact(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s,vcts)
    case('uct_hll')
      call update_faces_hll(ixI^L,ixO^L,qt,qdt,fE,sCT,s,vcts)
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select
  end subroutine rmhd_update_faces

  !> get electric field though averaging neighors to update faces in CT
  subroutine update_faces_average(ixI^L,ixO^L,qt,qdt,fC,fE,sCT,s)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,sdim:3)
    double precision                   :: circ(ixI^S,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixI^S,sdim:3) :: E_resi
    integer                            :: ix^D,ixC^L,ixA^L,i1kr^D,i2kr^D
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(bfaces=>s%ws,x=>s%x)
    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ! if there is resistivity, get eta J
    if(rmhd_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)
    do idim1=1,ndim
      iwdim1 = mag(idim1)
      i1kr^D=kr(idim1,^D);
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        i2kr^D=kr(idim2,^D);
        do idir=sdim,3! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ixCmax^D=ixOmax^D;
            ixCmin^D=ixOmin^D+kr(idir,^D)-1;
            ! average cell-face electric field to cell edges
           {do ix^DB=ixCmin^DB,ixCmax^DB\}
              fE(ix^D,idir)=quarter*&
                (fC(ix^D,iwdim1,idim2)+fC({ix^D+i1kr^D},iwdim1,idim2)&
                -fC(ix^D,iwdim2,idim1)-fC({ix^D+i2kr^D},iwdim2,idim1))
              ! add resistive electric field at cell edges E=-vxB+eta J
              if(rmhd_eta/=zero) fE(ix^D,idir)=fE(ix^D,idir)+E_resi(ix^D,idir)
              ! times time step and edge length
              fE(ix^D,idir)=fE(ix^D,idir)*qdt*s%dsC(ix^D,idir)
           {end do\}
          end if
        end do
      end do
    end do
    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) &
      call usr_set_electric_field(ixI^L,ixO^L,qt,qdt,fE,sCT)
    circ(ixI^S,1:ndim)=zero
    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        ixA^L=ixC^L-kr(idim2,^D);
        do idir=sdim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)==1) then
            ! Add line integrals in direction idir
            circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                             +(fE(ixC^S,idir)&
                              -fE(ixA^S,idir))
          else if(lvc(idim1,idim2,idir)==-1) then
            ! Add line integrals in direction idir
            circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                             -(fE(ixC^S,idir)&
                              -fE(ixA^S,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do
    end associate
  end subroutine update_faces_average

  !> update faces using UCT contact mode by Gardiner and Stone 2005 JCP 205, 509
  subroutine update_faces_contact(ixI^L,ixO^L,qt,qdt,wp,fC,fE,sCT,s,vcts)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wp(ixI^S,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,sdim:3)
    double precision                   :: circ(ixI^S,1:ndim)
    ! electric field at cell centers
    double precision                   :: ECC(ixI^S,sdim:3)
    double precision                   :: Ein(ixI^S,sdim:3)
    ! gradient of E at left and right side of a cell face
    double precision                   :: EL(ixI^S),ER(ixI^S)
    ! gradient of E at left and right side of a cell corner
    double precision                   :: ELC,ERC
    ! non-ideal electric field on cell edges
    double precision, dimension(ixI^S,sdim:3) :: E_resi
    ! current on cell edges
    double precision :: jce(ixI^S,sdim:3)
    ! location at cell faces
    double precision :: xs(ixGs^T,1:ndim)
    double precision :: gradi(ixGs^T)
    integer :: ixC^L,ixA^L
    integer :: idim1,idim2,idir,iwdim1,iwdim2,ix^D,i1kr^D,i2kr^D

    associate(bfaces=>s%ws,x=>s%x,w=>s%w,vnorm=>vcts%vnorm,wCTs=>sCT%ws)
    ! if there is resistivity, get eta J
    if(rmhd_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)
    if(B0field) then
     {do ix^DB=ixImin^DB,ixImax^DB\}
        ! Calculate electric field at cell centers
       {^IFTHREED
        ECC(ix^D,1)=(wp(ix^D,b2_)+block%B0(ix^D,2,0))*wp(ix^D,m3_)-(wp(ix^D,b3_)+block%B0(ix^D,3,0))*wp(ix^D,m2_)
        ECC(ix^D,2)=(wp(ix^D,b3_)+block%B0(ix^D,3,0))*wp(ix^D,m1_)-(wp(ix^D,b1_)+block%B0(ix^D,1,0))*wp(ix^D,m3_)
        ECC(ix^D,3)=(wp(ix^D,b1_)+block%B0(ix^D,1,0))*wp(ix^D,m2_)-(wp(ix^D,b2_)+block%B0(ix^D,2,0))*wp(ix^D,m1_)
       }
       {^IFTWOD
        ECC(ix^D,3)=wp(ix^D,b1_)*wp(ix^D,m2_)-wp(ix^D,b2_)*wp(ix^D,m1_)
       }
       {^IFONED
        ECC(ix^D,3)=0.d0
       }
     {end do\}
    else
     {do ix^DB=ixImin^DB,ixImax^DB\}
        ! Calculate electric field at cell centers
       {^IFTHREED
        ECC(ix^D,1)=wp(ix^D,b2_)*wp(ix^D,m3_)-wp(ix^D,b3_)*wp(ix^D,m2_)
        ECC(ix^D,2)=wp(ix^D,b3_)*wp(ix^D,m1_)-wp(ix^D,b1_)*wp(ix^D,m3_)
        ECC(ix^D,3)=wp(ix^D,b1_)*wp(ix^D,m2_)-wp(ix^D,b2_)*wp(ix^D,m1_)
       }
       {^IFTWOD
        ECC(ix^D,3)=wp(ix^D,b1_)*wp(ix^D,m2_)-wp(ix^D,b2_)*wp(ix^D,m1_)
       }
       {^IFONED
        ECC(ix^D,3)=0.d0
       }
     {end do\}
    end if

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ! evaluate electric field along cell edges according to equation (41)
    do idim1=1,ndim
      iwdim1 = mag(idim1)
      i1kr^D=kr(idim1,^D);
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        i2kr^D=kr(idim2,^D);
        do idir=sdim,3 ! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ixCmax^D=ixOmax^D;
            ixCmin^D=ixOmin^D+kr(idir,^D)-1;
            ! Assemble indices
            ! average cell-face electric field to cell edges
           {do ix^DB=ixCmin^DB,ixCmax^DB\}
              fE(ix^D,idir)=quarter*&
              (fC(ix^D,iwdim1,idim2)+fC({ix^D+i1kr^D},iwdim1,idim2)&
              -fC(ix^D,iwdim2,idim1)-fC({ix^D+i2kr^D},iwdim2,idim1))
              if(partial_energy) Ein(ix^D,idir)=fE(ix^D,idir)
           {end do\}
            ! add slope in idim2 direction from equation (50)
            ixAmin^D=ixCmin^D;
            ixAmax^D=ixCmax^D+i1kr^D;
           {do ix^DB=ixAmin^DB,ixAmax^DB\}
              EL(ix^D)=fC(ix^D,iwdim1,idim2)-ECC(ix^D,idir)
              ER(ix^D)=fC(ix^D,iwdim1,idim2)-ECC({ix^D+i2kr^D},idir)
           {end do\}
           {!dir$ ivdep
            do ix^DB=ixCmin^DB,ixCmax^DB\}
              if(vnorm(ix^D,idim1)>0.d0) then
                ELC=EL(ix^D)
              else if(vnorm(ix^D,idim1)<0.d0) then
                ELC=EL({ix^D+i1kr^D})
              else
                ELC=0.5d0*(EL(ix^D)+EL({ix^D+i1kr^D}))
              end if
              if(vnorm({ix^D+i2kr^D},idim1)>0.d0) then
                ERC=ER(ix^D)
              else if(vnorm({ix^D+i2kr^D},idim1)<0.d0) then
                ERC=ER({ix^D+i1kr^D})
              else
                ERC=0.5d0*(ER(ix^D)+ER({ix^D+i1kr^D}))
              end if
              fE(ix^D,idir)=fE(ix^D,idir)+0.25d0*(ELC+ERC)
           {end do\}
            ! add slope in idim1 direction from equation (50)
            ixAmin^D=ixCmin^D;
            ixAmax^D=ixCmax^D+i2kr^D;
           {do ix^DB=ixAmin^DB,ixAmax^DB\}
              EL(ix^D)=-fC(ix^D,iwdim2,idim1)-ECC(ix^D,idir)
              ER(ix^D)=-fC(ix^D,iwdim2,idim1)-ECC({ix^D+i1kr^D},idir)
           {end do\}
           {!dir$ ivdep
            do ix^DB=ixCmin^DB,ixCmax^DB\}
              if(vnorm(ix^D,idim2)>0.d0) then
                ELC=EL(ix^D)
              else if(vnorm(ix^D,idim2)<0.d0) then
                ELC=EL({ix^D+i2kr^D})
              else
                ELC=0.5d0*(EL(ix^D)+EL({ix^D+i2kr^D}))
              end if
              if(vnorm({ix^D+i1kr^D},idim2)>0.d0) then
                ERC=ER(ix^D)
              else if(vnorm({ix^D+i1kr^D},idim2)<0.d0) then
                ERC=ER({ix^D+i2kr^D})
              else
                ERC=0.5d0*(ER(ix^D)+ER({ix^D+i2kr^D}))
              end if
              fE(ix^D,idir)=fE(ix^D,idir)+0.25d0*(ELC+ERC)
              ! difference between average and upwind interpolated E
              if(partial_energy) Ein(ix^D,idir)=fE(ix^D,idir)-Ein(ix^D,idir)
              ! add resistive electric field at cell edges E=-vxB+eta J
              if(rmhd_eta/=zero) fE(ix^D,idir)=fE(ix^D,idir)+E_resi(ix^D,idir)
              ! times time step and edge length
              fE(ix^D,idir)=fE(ix^D,idir)*qdt*s%dsC(ix^D,idir)
           {end do\}
          end if
        end do
      end do
    end do
    if(partial_energy) then
      ! add upwind diffused magnetic energy back to energy
      ! calculate current density at cell edges
      jce=0.d0
      do idim1=1,ndim
        do idim2=1,ndim
          do idir=sdim,3
            if (lvc(idim1,idim2,idir)==0) cycle
            ixCmax^D=ixOmax^D;
            ixCmin^D=ixOmin^D+kr(idir,^D)-1;
            ixAmax^D=ixCmax^D-kr(idir,^D)+1;
            ixAmin^D=ixCmin^D;
            ! current at transverse faces
            xs(ixA^S,:)=x(ixA^S,:)
            xs(ixA^S,idim2)=x(ixA^S,idim2)+half*s%dx(ixA^S,idim2)
            call gradientF(wCTs(ixGs^T,idim2),xs,ixGs^LL,ixC^L,idim1,gradi)
            if (lvc(idim1,idim2,idir)==1) then
              jce(ixC^S,idir)=jce(ixC^S,idir)+gradi(ixC^S)
            else
              jce(ixC^S,idir)=jce(ixC^S,idir)-gradi(ixC^S)
            end if
          end do
        end do
      end do
      do idir=sdim,3
        ixCmax^D=ixOmax^D;
        ixCmin^D=ixOmin^D+kr(idir,^D)-1;
        ! E dot J on cell edges
        Ein(ixC^S,idir)=Ein(ixC^S,idir)*jce(ixC^S,idir)
        ! average from cell edge to cell center
       {^IFTHREED
        if(idir==1) then
         {do ix^DB=ixOmin^DB,ixOmax^DB\}
            jce(ix^D,idir)=0.25d0*(Ein(ix^D,idir)+Ein(ix1,ix2-1,ix3,idir)+Ein(ix1,ix2,ix3-1,idir)&
                          +Ein(ix1,ix2-1,ix3-1,idir))
            if(jce(ix^D,idir)<0.d0) jce(ix^D,idir)=0.d0
            w(ix^D,e_)=w(ix^D,e_)+qdt*jce(ix^D,idir)
         {end do\}
        else if(idir==2) then
         {do ix^DB=ixOmin^DB,ixOmax^DB\}
            jce(ix^D,idir)=0.25d0*(Ein(ix^D,idir)+Ein(ix1-1,ix2,ix3,idir)+Ein(ix1,ix2,ix3-1,idir)&
                          +Ein(ix1-1,ix2,ix3-1,idir))
            if(jce(ix^D,idir)<0.d0) jce(ix^D,idir)=0.d0
            w(ix^D,e_)=w(ix^D,e_)+qdt*jce(ix^D,idir)
         {end do\}
        else
         {do ix^DB=ixOmin^DB,ixOmax^DB\}
            jce(ix^D,idir)=0.25d0*(Ein(ix^D,idir)+Ein(ix1-1,ix2,ix3,idir)+Ein(ix1,ix2-1,ix3,idir)&
                          +Ein(ix1-1,ix2-1,ix3,idir))
            if(jce(ix^D,idir)<0.d0) jce(ix^D,idir)=0.d0
            w(ix^D,e_)=w(ix^D,e_)+qdt*jce(ix^D,idir)
         {end do\}
        end if
       }
       {^IFTWOD
        !idir=3
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          jce(ix^D,idir)=0.25d0*(Ein(ix^D,idir)+Ein(ix1-1,ix2,idir)+Ein(ix1,ix2-1,idir)&
                        +Ein(ix1-1,ix2-1,idir))
          if(jce(ix^D,idir)<0.d0) jce(ix^D,idir)=0.d0
          w(ix^D,e_)=w(ix^D,e_)+qdt*jce(ix^D,idir)
       {end do\}
       }
        ! save additional numerical resistive heating to an extra variable
        if(nwextra>0) then
          block%w(ixO^S,nw)=block%w(ixO^S,nw)+jce(ixO^S,idir)
        end if
      end do
    end if
    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) &
      call usr_set_electric_field(ixI^L,ixO^L,qt,qdt,fE,sCT)
    circ(ixI^S,1:ndim)=zero
    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        ixA^L=ixC^L-kr(idim2,^D);
        do idir=sdim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)==1) then
            ! Add line integrals in direction idir
            circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                             +(fE(ixC^S,idir)&
                              -fE(ixA^S,idir))
          else if(lvc(idim1,idim2,idir)==-1) then
            ! Add line integrals in direction idir
            circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                             -(fE(ixC^S,idir)&
                              -fE(ixA^S,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      where(s%surfaceC(ixC^S,idim1) > smalldouble)
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do
    end associate
  end subroutine update_faces_contact

  !> update faces
  subroutine update_faces_hll(ixI^L,ixO^L,qt,qdt,fE,sCT,s,vcts)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_usr_methods
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    double precision, intent(inout)    :: fE(ixI^S,sdim:3)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision                   :: vtilL(ixI^S,2)
    double precision                   :: vtilR(ixI^S,2)
    double precision                   :: bfacetot(ixI^S,ndim)
    double precision                   :: btilL(ixI^S,ndim)
    double precision                   :: btilR(ixI^S,ndim)
    double precision                   :: cp(ixI^S,2)
    double precision                   :: cm(ixI^S,2)
    double precision                   :: circ(ixI^S,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixI^S,sdim:3) :: E_resi
    integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,ixCm^L
    integer                            :: idim1,idim2,idir

    associate(bfaces=>s%ws,bfacesCT=>sCT%ws,x=>s%x,vbarC=>vcts%vbarC,cbarmin=>vcts%cbarmin,&
      cbarmax=>vcts%cbarmax)
    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! Loop over components of electric field

    ! idir: electric field component we need to calculate
    ! idim1: directions in which we already performed the reconstruction
    ! idim2: directions in which we perform the reconstruction

    ! if there is resistivity, get eta J
    if(rmhd_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)

    do idir=sdim,3
      ! Indices
      ! idir: electric field component
      ! idim1: one surface
      ! idim2: the other surface
      ! cyclic permutation: idim1,idim2,idir=1,2,3
      ! Velocity components on the surface
      ! follow cyclic premutations:
      ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-1+kr(idir,^D);
      ! Set indices and directions
      idim1=mod(idir,3)+1
      idim2=mod(idir+1,3)+1
      jxC^L=ixC^L+kr(idim1,^D);
      ixCp^L=ixC^L+kr(idim2,^D);
      ! Reconstruct transverse transport velocities
      call reconstruct(ixI^L,ixC^L,idim2,vbarC(ixI^S,idim1,1),&
               vtilL(ixI^S,2),vtilR(ixI^S,2))
      call reconstruct(ixI^L,ixC^L,idim1,vbarC(ixI^S,idim2,2),&
               vtilL(ixI^S,1),vtilR(ixI^S,1))
      ! Reconstruct magnetic fields
      ! Eventhough the arrays are larger, reconstruct works with
      ! the limits ixG.
      if(B0field) then
        bfacetot(ixI^S,idim1)=bfacesCT(ixI^S,idim1)+block%B0(ixI^S,idim1,idim1)
        bfacetot(ixI^S,idim2)=bfacesCT(ixI^S,idim2)+block%B0(ixI^S,idim2,idim2)
      else
        bfacetot(ixI^S,idim1)=bfacesCT(ixI^S,idim1)
        bfacetot(ixI^S,idim2)=bfacesCT(ixI^S,idim2)
      end if
      call reconstruct(ixI^L,ixC^L,idim2,bfacetot(ixI^S,idim1),&
               btilL(ixI^S,idim1),btilR(ixI^S,idim1))
      call reconstruct(ixI^L,ixC^L,idim1,bfacetot(ixI^S,idim2),&
               btilL(ixI^S,idim2),btilR(ixI^S,idim2))
      ! Take the maximum characteristic
      cm(ixC^S,1)=max(cbarmin(ixCp^S,idim1),cbarmin(ixC^S,idim1))
      cp(ixC^S,1)=max(cbarmax(ixCp^S,idim1),cbarmax(ixC^S,idim1))
      cm(ixC^S,2)=max(cbarmin(jxC^S,idim2),cbarmin(ixC^S,idim2))
      cp(ixC^S,2)=max(cbarmax(jxC^S,idim2),cbarmax(ixC^S,idim2))
      ! Calculate eletric field
      fE(ixC^S,idir)=-(cp(ixC^S,1)*vtilL(ixC^S,1)*btilL(ixC^S,idim2) &
                     + cm(ixC^S,1)*vtilR(ixC^S,1)*btilR(ixC^S,idim2) &
                     - cp(ixC^S,1)*cm(ixC^S,1)*(btilR(ixC^S,idim2)-btilL(ixC^S,idim2)))&
                     /(cp(ixC^S,1)+cm(ixC^S,1)) &
                     +(cp(ixC^S,2)*vtilL(ixC^S,2)*btilL(ixC^S,idim1) &
                     + cm(ixC^S,2)*vtilR(ixC^S,2)*btilR(ixC^S,idim1) &
                     - cp(ixC^S,2)*cm(ixC^S,2)*(btilR(ixC^S,idim1)-btilL(ixC^S,idim1)))&
                     /(cp(ixC^S,2)+cm(ixC^S,2))
      ! add resistive electric field at cell edges E=-vxB+eta J
      if(rmhd_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+E_resi(ixC^S,idir)
      fE(ixC^S,idir)=qdt*s%dsC(ixC^S,idir)*fE(ixC^S,idir)
      if (.not.slab) then
        where(abs(x(ixC^S,r_)+half*dxlevel(r_)).lt.1.0d-9)
          fE(ixC^S,idir)=zero
        end where
      end if
    end do
    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) &
      call usr_set_electric_field(ixI^L,ixO^L,qt,qdt,fE,sCT)
    circ(ixI^S,1:ndim)=zero
    ! Calculate circulation on each face: interal(fE dot dl)
    do idim1=1,ndim ! Coordinate perpendicular to face
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        do idir=sdim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)/=0) then
            hxC^L=ixC^L-kr(idim2,^D);
            ! Add line integrals in direction idir
            circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                             +lvc(idim1,idim2,idir)&
                             *(fE(ixC^S,idir)&
                              -fE(hxC^S,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do
    end associate
  end subroutine update_faces_hll

  !> calculate eta J at cell edges
  subroutine get_resistive_electric_field(ixI^L,ixO^L,sCT,s,jce)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    integer, intent(in)                :: ixI^L, ixO^L
    type(state), intent(in)            :: sCT, s
    ! current on cell edges
    double precision :: jce(ixI^S,sdim:3)
    ! current on cell centers
    double precision :: jcc(ixI^S,7-2*ndir:3)
    ! location at cell faces
    double precision :: xs(ixGs^T,1:ndim)
    ! resistivity
    double precision :: eta(ixI^S)
    double precision :: gradi(ixGs^T)
    integer :: ix^D,ixC^L,ixA^L,ixB^L,idir,idirmin,idim1,idim2

    associate(x=>s%x,dx=>s%dx,w=>s%w,wCT=>sCT%w,wCTs=>sCT%ws)
    ! calculate current density at cell edges
    jce=0.d0
    do idim1=1,ndim
      do idim2=1,ndim
        do idir=sdim,3
          if (lvc(idim1,idim2,idir)==0) cycle
          ixCmax^D=ixOmax^D;
          ixCmin^D=ixOmin^D+kr(idir,^D)-1;
          ixBmax^D=ixCmax^D-kr(idir,^D)+1;
          ixBmin^D=ixCmin^D;
          ! current at transverse faces
          xs(ixB^S,:)=x(ixB^S,:)
          xs(ixB^S,idim2)=x(ixB^S,idim2)+half*dx(ixB^S,idim2)
          call gradientF(wCTs(ixGs^T,idim2),xs,ixGs^LL,ixC^L,idim1,gradi,2)
          if (lvc(idim1,idim2,idir)==1) then
            jce(ixC^S,idir)=jce(ixC^S,idir)+gradi(ixC^S)
          else
            jce(ixC^S,idir)=jce(ixC^S,idir)-gradi(ixC^S)
          end if
        end do
      end do
    end do
    ! get resistivity
    if(rmhd_eta>zero)then
      jce(ixI^S,:)=jce(ixI^S,:)*rmhd_eta
    else
      ixA^L=ixO^L^LADD1;
      call get_current(wCT,ixI^L,ixA^L,idirmin,jcc)
      call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,jcc,eta)
      ! calcuate eta on cell edges
      do idir=sdim,3
        ixCmax^D=ixOmax^D;
        ixCmin^D=ixOmin^D+kr(idir,^D)-1;
        jcc(ixC^S,idir)=0.d0
       {do ix^DB=0,1\}
          if({ ix^D==1 .and. ^D==idir | .or.}) cycle
          ixAmin^D=ixCmin^D+ix^D;
          ixAmax^D=ixCmax^D+ix^D;
          jcc(ixC^S,idir)=jcc(ixC^S,idir)+eta(ixA^S)
       {end do\}
        jcc(ixC^S,idir)=jcc(ixC^S,idir)*0.25d0
        jce(ixC^S,idir)=jce(ixC^S,idir)*jcc(ixC^S,idir)
      end do
    end if
    end associate
  end subroutine get_resistive_electric_field

  !> calculate cell-center values from face-center values
  subroutine rmhd_face_to_center(ixO^L,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in) :: ixO^L
    type(state) :: s
    integer :: ix^D

    ! calculate cell-center values from face-center values in 2nd order
    ! because the staggered arrays have an additional place to the left.
    ! Interpolate to cell barycentre using arithmetic average
    ! This might be done better later, to make the method less diffusive.
    {!dir$ ivdep
    do ix^DB=ixOmin^DB,ixOmax^DB\}
      {^IFTHREED
      s%w(ix^D,b1_)=half/s%surface(ix^D,1)*(s%ws(ix^D,1)*s%surfaceC(ix^D,1)&
        +s%ws(ix1-1,ix2,ix3,1)*s%surfaceC(ix1-1,ix2,ix3,1))
      s%w(ix^D,b2_)=half/s%surface(ix^D,2)*(s%ws(ix^D,2)*s%surfaceC(ix^D,2)&
        +s%ws(ix1,ix2-1,ix3,2)*s%surfaceC(ix1,ix2-1,ix3,2))
      s%w(ix^D,b3_)=half/s%surface(ix^D,3)*(s%ws(ix^D,3)*s%surfaceC(ix^D,3)&
        +s%ws(ix1,ix2,ix3-1,3)*s%surfaceC(ix1,ix2,ix3-1,3))
      }
      {^IFTWOD
      s%w(ix^D,b1_)=half/s%surface(ix^D,1)*(s%ws(ix^D,1)*s%surfaceC(ix^D,1)&
        +s%ws(ix1-1,ix2,1)*s%surfaceC(ix1-1,ix2,1))
      s%w(ix^D,b2_)=half/s%surface(ix^D,2)*(s%ws(ix^D,2)*s%surfaceC(ix^D,2)&
        +s%ws(ix1,ix2-1,2)*s%surfaceC(ix1,ix2-1,2))
      }
   {end do\}
    ! calculate cell-center values from face-center values in 4th order
    !do idim=1,ndim
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);

    !  ! Interpolate to cell barycentre using fourth order central formula
    !  w(ixO^S,mag(idim))=(0.0625d0/s%surface(ixO^S,idim))*&
    !         ( -ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !     +9.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !     +9.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !           -ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) )
    !end do

    ! calculate cell-center values from face-center values in 6th order
    !do idim=1,ndim
    !  fxO^L=ixO^L-3*kr(idim,^D);
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);
    !  kxO^L=ixO^L+2*kr(idim,^D);

    !  ! Interpolate to cell barycentre using sixth order central formula
    !  w(ixO^S,mag(idim))=(0.00390625d0/s%surface(ixO^S,idim))* &
    !     (  +3.0d0*ws(fxO^S,idim)*s%surfaceC(fxO^S,idim) &
    !       -25.0d0*ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !      +150.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !      +150.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !       -25.0d0*ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) &
    !        +3.0d0*ws(kxO^S,idim)*s%surfaceC(kxO^S,idim) )
    !end do
  end subroutine rmhd_face_to_center

  !> calculate magnetic field from vector potential
  subroutine b_from_vector_potential(ixIs^L, ixI^L, ixO^L, ws, x)
    use mod_global_parameters
    use mod_constrained_transport
    integer, intent(in)                :: ixIs^L, ixI^L, ixO^L
    double precision, intent(inout)    :: ws(ixIs^S,1:nws)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: Adummy(ixIs^S,1:3)

    call b_from_vector_potentialA(ixIs^L, ixI^L, ixO^L, ws, x, Adummy)
  end subroutine b_from_vector_potential

  subroutine Rfactor_from_temperature_ionization(w,x,ixI^L,ixO^L,Rfactor)
    use mod_global_parameters
    use mod_ionization_degree
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: Rfactor(ixI^S)
    double precision :: iz_H(ixO^S),iz_He(ixO^S)

    call ionization_degree_from_temperature(ixI^L,ixO^L,w(ixI^S,Te_),iz_H,iz_He)
    ! assume the first and second ionization of Helium have the same degree
    Rfactor(ixO^S)=(1.d0+iz_H(ixO^S)+0.1d0*(1.d0+iz_He(ixO^S)*(1.d0+iz_He(ixO^S))))/(2.d0+3.d0*He_abundance)
  end subroutine Rfactor_from_temperature_ionization

  subroutine Rfactor_from_constant_ionization(w,x,ixI^L,ixO^L,Rfactor)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: Rfactor(ixI^S)

    Rfactor(ixO^S)=RR
  end subroutine Rfactor_from_constant_ionization
end module mod_rmhd_phys
