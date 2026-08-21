!> Magneto-hydrodynamics module
module mod_mhd_phys

#include "amrvac.h"

  use mod_global_parameters, only: std_len, const_c
  use mod_thermal_conduction, only: tc_fluid
  use mod_radiative_cooling, only: rc_fluid
  use mod_thermal_emission, only: te_fluid
  use mod_fld, only: fld_fluid
  use mod_physics
  use mod_eos
  use mod_comm_lib, only: mpistop
  use mod_functions_bfield, only: get_divb, mag

  implicit none
  private

  !> The adiabatic constant
  double precision, public                :: mhd_adiab = 1.0d0
  !> The MHD resistivity
  double precision, public                :: mhd_eta = 0.0d0
  !> The MHD hyper-resistivity
  double precision, public                :: mhd_eta_hyper = 0.0d0
  !> Hall resistivity
  double precision, public                :: mhd_etah = 0.0d0
  !> The MHD ambipolar coefficient
  double precision, public                :: mhd_eta_ambi = 0.0d0
  !> Height of the mask used in the TRAC method
  double precision, public, protected     :: mhd_trac_mask = 0.d0
  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: mhd_glm_alpha = 0.5d0
  !> Reduced speed of light for semirelativistic MHD: 2% of light speed
  double precision, public, protected     :: mhd_reduced_c = 0.02d0*const_c
  !> The thermal conductivity kappa in hyperbolic thermal conduction
  double precision, public                :: mhd_hyperbolic_tc_kappa = 0.d0
  logical, public                         :: mhd_hyperbolic_tc_constant = .false.
  !> Coefficient of diffusive divB cleaning
  double precision, public                :: divbdiff     = 0.8d0
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
  !> inverse of squared speed of light c0 and reduced speed of light c
  double precision :: inv_squared_c0=0.d0, inv_squared_c=0.d0
  !> equi vars indices in the state%equi_vars array
  integer, public :: equi_rho0_ = -1
  integer, public :: equi_pe0_ = -1
  !> Number of tracer species
  integer, public, protected              :: mhd_n_tracer = 0
  !> Index of the density (in the w array)
  integer, public, protected              :: rho_
  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)
  !> Indices of the momentum density for the form of better vectorization
  integer, public, protected              :: ^C&m^C_
  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_
  !> Indices of the magnetic field for the form of better vectorization
  integer, public, protected              :: ^C&b^C_
  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_
  !> Index of the electron number density for LTE module
  integer, public, protected :: Ne_
  !> Index of the field-aligned heat flux q_parallel
  integer, public, protected :: qpar_
  !> Index of the perpendicular heat flux q_perp
  integer, public, protected :: qperp_
  !> Indices of the GLM psi
  integer, public, protected :: psi_
  !> Index of the radiation energy
  integer, public, protected :: r_e
  !> Indices of temperature
  integer, public, protected :: Te_
  !> Index of the FIP passive scalar rho*fip in conserved form, fip in primitive form
  integer, public, protected :: fip_ = -1
  !> Whether FIP passive scalar is enabled
  logical, public, protected :: mhd_fip = .false.
  !> Enable the Uniturbulence and Alfven Wave Solar Model extension
  logical, public, protected :: mhd_uawsom = .false.
  !> Enable conservative Alfvén reflection (one-dimensional or Cartesian gradient-vorticity)
  logical, public, protected :: mhd_uawsom_reflection = .false.
  !> Algebra used by the Alfven-wave reflection source.  one_dimensional_gradient
  !> retains the original one-dimensional source; cartesian_gradient_vorticity
  !> is the multidimensional gradient/vorticity closure.
  character(len=std_len), public, protected :: mhd_uawsom_reflection_mode = 'one_dimensional_gradient'
  !> Enable conservative kink-wave reflection from the kink-speed gradient
  logical, public, protected :: mhd_uawsom_kink_reflection = .false.
  !> Dimensionless multiplier in the one_dimensional_gradient source only
  double precision, public, protected :: mhd_uawsom_sigma = 0.d0
  !> Cartesian direction used for the prescribed density contrast and reflection gradient
  integer, public, protected :: mhd_uawsom_height_dim = -1
  !> Base transverse density contrast and filling factor
  double precision, public, protected :: mhd_uawsom_zeta0 = 5.d0
  double precision, public, protected :: mhd_uawsom_filling_factor = 0.1d0
  !> Physical input scales (lengths and magnetic field); converted to code units at init
  double precision, public, protected :: mhd_uawsom_zeta_scale = -1.d0
  double precision, public, protected :: mhd_uawsom_thread_radius0 = -1.d0
  double precision, public, protected :: mhd_uawsom_Bref = -1.d0
  double precision, public, protected :: mhd_uawsom_alfven_corr_length0 = -1.d0
  !> Conserved wave-energy indices.  The plus variables propagate against B.
  integer, public, protected :: wAplus_=-1, wAminus_=-1, wkplus_=-1, wkminus_=-1
  !> Index of the cutoff temperature for the TRAC method
  integer, public, protected              :: Tcoff_
  integer, public, protected              :: Tweight_
  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)
  !> The number of waves
  integer :: nwwave=8
  !> Method type of divb in a integer for good performance
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
  logical, public, protected              :: mhd_energy = .true.
  !> Whether thermal conduction is used
  logical, public, protected              :: mhd_thermal_conduction = .false.
  !> Whether radiative cooling is added
  logical, public, protected              :: mhd_radiative_cooling = .false.
  !> Whether thermal conduction is used
  logical, public, protected              :: mhd_hyperbolic_tc = .false.
  !> Whether saturation is considered for hyperbolic TC. When the perpendicular
  !> channel is active, the limiter uses sqrt(q_parallel^2+q_perp^2).
  logical, public, protected              :: mhd_hyperbolic_tc_sat = .false.
  !> Whether the perpendicular hyperbolic-TC channel is enabled
  logical, public, protected              :: mhd_hyperbolic_tc_use_perp = .false.
  !> Perpendicular hyperbolic-TC closure mode:
  !> 'off' = disabled
  !> 'fixed_reference' = fixed classical ratio at the reference normalisation
  !> 'weak_field_isotropization' = empirical transition using a prescribed Bmin
  !> 'electron_magnetization' = simplified ratio 1/(1+chi_e^2)
  character(len=std_len), public, protected :: mhd_hyperbolic_tc_perp_mode = 'off'
  !> Relative perpendicular hyperbolic-TC coefficient in fixed/strong-field limit:
  !> kappa_perp0 = mhd_hyperbolic_tc_kappa_perp_factor * kappa_parallel
  double precision, public, protected     :: mhd_hyperbolic_tc_kappa_perp_factor = 0.d0
  !> Field-strength transition scale for perpendicular closure
  double precision, public, protected     :: mhd_hyperbolic_tc_Bmin = 0.d0
  !> Constant Coulomb logarithm used by the simplified electron-magnetization
  !> closure. It is a namelist parameter so changing it has no per-cell cost.
  double precision, public, protected     :: mhd_hyperbolic_tc_coulomb_log = 20.d0
  !> Whether viscosity is added
  logical, public, protected              :: mhd_viscosity = .false.
  !> Whether gravity is added
  logical, public, protected              :: mhd_gravity = .false.
  !> Whether rotating frame is activated
  logical, public, protected              :: mhd_rotating_frame = .false.
  !> Whether Hall-MHD is used
  logical, public, protected              :: mhd_hall = .false.
  !> Whether Ambipolar term is used
  logical, public, protected              :: mhd_ambipolar = .false.
  !> Whether Ambipolar term is implemented using supertimestepping
  logical, public, protected              :: mhd_ambipolar_sts = .false.
  !> Whether Ambipolar term is implemented explicitly
  logical, public, protected              :: mhd_ambipolar_exp = .false.
  !> Whether particles module is added
  logical, public, protected              :: mhd_particles = .false.
  !> Whether magnetofriction is added
  logical, public, protected              :: mhd_magnetofriction = .false.
  !> Whether GLM-MHD is used to control div B
  logical, public, protected              :: mhd_glm = .false.
  !> Whether extended GLM-MHD is used with additional sources
  logical, public, protected              :: mhd_glm_extended = .true.
  !> Whether TRAC method is used
  logical, public, protected              :: mhd_trac = .false.
  !> Which TRAC method is used
  integer, public, protected              :: mhd_trac_type=1
  !> Distance between two adjacent traced magnetic field lines (in finest cell size)
  integer, public, protected              :: mhd_trac_finegrid=4
  !> TRAC-7 (Johnston 2021, A&A 654 A2): target number of cells resolving the TR
  double precision, public, protected     :: mhd_trac_delta=0.5d0
  !> Whether internal energy is solved instead of total energy
  logical, public, protected              :: mhd_internal_e = .false.
  !> Whether hydrodynamic energy is solved instead of total energy
  logical, public, protected              :: mhd_hydrodynamic_e = .false.
  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.
  !> Whether semirelativistic MHD equations (Gombosi 2002 JCP) are solved
  logical, public, protected              :: mhd_semirelativistic = .false.
  !> Whether plasma is partially ionized
  !> Whether CAK radiation line force is activated
  logical, public, protected              :: mhd_cak_force = .false.
  !> Whether radiation-gas interaction is handled using flux limited diffusion
  logical, public, protected              :: mhd_radiation_fld = .false.
  !> Radiation fluid object (gas-EoS callbacks for FLD), wired in mhd_link_eos
  type(fld_fluid), allocatable, public    :: fld_fl
  !> whether split off equilibrium density and pressure
  logical, public :: has_equi_rho_and_p = .false.
  logical, public :: mhd_equi_thermal = .false.
  !> whether dump full variables (when splitting is used) in a separate dat file
  logical, public, protected              :: mhd_dump_full_vars = .false.
  !> Whether divB is computed with a fourth order approximation
  integer, public, protected :: mhd_divb_nth = 1
  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.
  !> clean initial divB
  logical, public :: clean_initial_divb     = .false.
  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*^ND)=.true.
  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.
  !> Whether an total energy equation is used
  logical :: total_energy = .true.
  !> Whether numerical resistive heating is included when solving partial energy equation
  logical, public :: numerical_resistive_heating = .false.
  !> Whether gravity work is included in energy equation
  logical :: gravity_energy
  !> Method type to clean divergence of B
  character(len=std_len), public, protected :: typedivbfix  = 'linde'
  !> Method type of constrained transport
  character(len=std_len), public, protected :: type_ct  = 'uct_contact'
  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'
  !> type of fluid for thermal conduction
  type(tc_fluid), public, allocatable     :: tc_fl
  !> type of fluid for thermal emission synthesis
  type(te_fluid), public, allocatable     :: te_fl_mhd
  !> type of fluid for radiative cooling
  type(rc_fluid), public, allocatable     :: rc_fl

  !define the subroutine interface for the ambipolar mask
  abstract interface

    subroutine mask_subroutine(ixI^L,ixO^L,w,x,res)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in) :: w(ixI^S,1:nw)
      double precision, intent(inout) :: res(ixI^S)
    end subroutine mask_subroutine

  end interface

  procedure(mask_subroutine), pointer  :: usr_mask_ambipolar => null()
  procedure(sub_small_values), pointer :: mhd_handle_small_values => null()
  ! Public methods
  public :: usr_mask_ambipolar
  public :: mhd_phys_init
  public :: mhd_get_v
  public :: mhd_get_rho
  public :: mhd_e_to_ei
  public :: mhd_ei_to_e
  public :: mhd_face_to_center
  public :: get_divb
  public :: get_current
  !> needed  public if we want to use the ambipolar coefficient in the user file
  public :: multiplyAmbiCoef
  public :: get_normalized_divb
  public :: b_from_vector_potential
  public :: mhd_mag_en_all
  public :: mhd_uawsom_wave_energy_cell
  public :: mhd_uawsom_rho2_factor
  public :: mhd_uawsom_wave_pressure_cell
  public :: mhd_uawsom_rho2_factor_cell
  {^NOONED
  public :: mhd_clean_divb_multigrid
  }
  ! Begin: following relevant for radiative MHD using FLD
  ! first four are local and only of interest for mod_usr applications
  ! where they can be used in diagnostics
  ! NOTE those with _prim expect primitives on entry
  public :: mhd_get_pthermal_plus_pradiation
  public :: mhd_get_csrad2
  public :: mhd_get_trad
  public :: mhd_get_pradiation_from_prim
  !    as pointer phys_get_csrad2
  public :: mhd_get_csrad2_prim
  ! End: following relevant for radiative MHD using FLD
  ! Removed orphan public declarations: mhd_get_Rfactor, mhd_get_temperature_from_prim,
  ! mhd_get_temperature_from_etot. These functions live in mod_mhd_eos.t and are
  ! reached via the eos% / phys_get_Rfactor / phys_get_tgas / phys_get_temperature
  ! procedure pointers bound by mod_mhd_eos:bind_eos_to_source.

contains

  !> Read this module"s parameters from a file
  subroutine mhd_read_params(files)
    use mod_global_parameters
    use mod_particles, only: particles_eta, particles_etah
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mhd_list/ mhd_energy, mhd_n_tracer, mhd_adiab,&
      mhd_eta, mhd_eta_hyper, mhd_etah, mhd_eta_ambi, mhd_glm_alpha, mhd_glm_extended, mhd_magnetofriction,&
      mhd_thermal_conduction, mhd_radiative_cooling, mhd_hall, mhd_ambipolar, mhd_ambipolar_sts, mhd_gravity,&
      mhd_rotating_frame,mhd_viscosity, typedivbfix, source_split_divb, divbdiff,&
      typedivbdiff, type_ct, divbwave, &
      H_ion_fr, He_ion_fr, He_ion_fr2, SI_unit, B0field ,mhd_dump_full_vars,&
      B0field_forcefree, Bdip, Bquad, Boct, Busr, mhd_particles,&
      particles_eta, particles_etah,has_equi_rho_and_p,mhd_equi_thermal,&
      boundary_divbfix, boundary_divbfix_skip, mhd_divb_nth, mhd_semirelativistic,&
      mhd_reduced_c, clean_initial_divb, mhd_internal_e, numerical_resistive_heating,&
      mhd_hydrodynamic_e, mhd_trac, mhd_trac_type, mhd_trac_mask, mhd_trac_finegrid, &
      mhd_trac_delta, mhd_cak_force, &
      mhd_hyperbolic_tc, mhd_hyperbolic_tc_sat, mhd_hyperbolic_tc_kappa, &
      mhd_hyperbolic_tc_use_perp, mhd_hyperbolic_tc_perp_mode, &
      mhd_hyperbolic_tc_kappa_perp_factor, mhd_hyperbolic_tc_Bmin, &
      mhd_hyperbolic_tc_coulomb_log, &
      mhd_radiation_fld, mhd_fip, mhd_uawsom, mhd_uawsom_reflection, &
      mhd_uawsom_reflection_mode, mhd_uawsom_kink_reflection, &
      mhd_uawsom_sigma, mhd_uawsom_height_dim, mhd_uawsom_zeta0, &
      mhd_uawsom_filling_factor, mhd_uawsom_zeta_scale, &
      mhd_uawsom_thread_radius0, mhd_uawsom_Bref, &
      mhd_uawsom_alfven_corr_length0

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mhd_list, end=111)
111    close(unitpar)
    end do

    ! He_abundance is set in eos_list and accessed via eos%He_abundance

  end subroutine mhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine mhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh

    integer                             :: er
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    character(len=name_len)             :: names(n_par)

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = eos%gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine mhd_write_info

  subroutine mhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_magnetofriction, only: magnetofriction_init
    use mod_rotating_frame, only: rotating_frame_init
    use mod_supertimestepping, only: sts_init, add_sts_method,&
            set_conversion_methods_to_head, set_error_handling_to_head
    use mod_cak_force, only: cak_init
    use mod_eos_PI_tables
    use mod_geometry
    use mod_usr_methods, only: usr_Rfactor, usr_get_heating
    {^NOONED
    use mod_multigrid_coupling
    }
    use mod_fld

    integer :: itr, idir

    call mhd_read_params(par_files)
    if(mhd_uawsom_height_dim < 0) mhd_uawsom_height_dim=ndim

    if(mhd_internal_e) then
      if(mhd_hydrodynamic_e) then
        mhd_hydrodynamic_e=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_hydrodynamic_e=F when mhd_internal_e=T'
      end if
      if(has_equi_rho_and_p) then
        has_equi_rho_and_p=.false.
        if(mype==0) write(*,*) 'WARNING: set has_equi_rho_and_p=F when mhd_internal_e=T'
      end if
    end if

    if(mhd_hydrodynamic_e) then
      if(mhd_internal_e) then
        mhd_internal_e=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_internal_e=F when mhd_hydrodynamic_e=T'
      end if
      if(B0field) then
        B0field=.false.
        if(mype==0) write(*,*) 'WARNING: set B0field=F when mhd_hydrodynamic_e=T'
      end if
      if(has_equi_rho_and_p) then
        has_equi_rho_and_p=.false.
        if(mype==0) write(*,*) 'WARNING: set has_equi_rho_and_p=F when mhd_hydrodynamic_e=T'
      end if
    end if

    if(mhd_semirelativistic) then
      if(B0field) then
        B0field=.false.
        if(mype==0) write(*,*) 'WARNING: set B0field=F when mhd_semirelativistic=T'
      endif
      if(has_equi_rho_and_p) then
        has_equi_rho_and_p=.false.
        if(mype==0) write(*,*) 'WARNING: set has_equi_rho_and_p=F when mhd_semirelativistic=T'
      end if
      if(mhd_hydrodynamic_e) then
        mhd_hydrodynamic_e=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_hydrodynamic_e=F when mhd_semirelativistic=T'
      end if
    end if

    if(.not. mhd_energy) then
      if(mhd_internal_e) then
        mhd_internal_e=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_internal_e=F when mhd_energy=F'
      end if
      if(mhd_hydrodynamic_e) then
        mhd_hydrodynamic_e=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_hydrodynamic_e=F when mhd_energy=F'
      end if
      if(mhd_thermal_conduction) then
        mhd_thermal_conduction=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_thermal_conduction=F when mhd_energy=F'
      end if
      if(mhd_hyperbolic_tc) then
        mhd_hyperbolic_tc=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_hyperbolic_tc=F when mhd_energy=F'
      end if
      if(mhd_radiative_cooling) then
        mhd_radiative_cooling=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_radiative_cooling=F when mhd_energy=F'
      end if
      if(mhd_trac) then
        mhd_trac=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_trac=F when mhd_energy=F'
      end if
      if(B0field) then
        B0field=.false.
        if(mype==0) write(*,*) 'WARNING: set B0field=F when mhd_energy=F'
      end if
      if(has_equi_rho_and_p) then
        has_equi_rho_and_p=.false.
        if(mype==0) write(*,*) 'WARNING: set has_equi_rho_and_p=F when mhd_energy=F'
      end if
    end if

    if(mhd_thermal_conduction .and. mhd_hyperbolic_tc) then
      mhd_thermal_conduction=.false.
      if(mype==0) write(*,*) 'WARNING: set either parabolic TC or hyperbolic TC to F'
      if(mype==0) write(*,*) 'WARNING: defaulting to only mhd_hyperbolic_tc=T'
    end if
    {^IFONED
    if (mhd_hyperbolic_tc .and. mhd_hyperbolic_tc_use_perp) then
      call mpistop("mhd_hyperbolic_tc_use_perp is not supported in 1D")
    end if
    }

    physics_type = "mhd"
    phys_energy=mhd_energy
    phys_internal_e=mhd_internal_e
    phys_trac=mhd_trac
    phys_trac_type=mhd_trac_type

    phys_gamma = eos%gamma
    phys_trac_finegrid=mhd_trac_finegrid

    if(mhd_energy) then
      if(mhd_internal_e.or.mhd_hydrodynamic_e) then
        total_energy=.false.
      else
        numerical_resistive_heating=.false.
        total_energy=.true.
      end if
    else
      total_energy=.false.
    end if
    phys_total_energy=total_energy
    if(mhd_energy) then
      if(mhd_internal_e) then
        gravity_energy=.false.
      else
        gravity_energy=.true.
      end if
    else
      gravity_energy=.false.
    end if

    {^IFONED
    if(mhd_trac .and. mhd_trac_type .gt. 2) then
      mhd_trac_type=1
      if(mype==0) write(*,*) 'WARNING: reset mhd_trac_type=1 for 1D simulation'
    end if
    }
    if(mhd_trac .and. mhd_trac_type .le. 4) then
      mhd_trac_mask=bigdouble
      if(mype==0) write(*,*) 'WARNING: set mhd_trac_mask==bigdouble for global TRAC method'
    end if
    if(mhd_trac .and. mhd_trac_type==7) then
      ! per-cell local TRAC (Johnston 2021) needs the background heating rate for the
      ! steady-state balance; it is field-aligned and per-cell, so no global trac_mask is used.
      if(.not. associated(usr_get_heating)) &
        call mpistop("mhd_trac_type=7 requires usr_get_heating to be set in mod_usr.t")
    end if
    phys_trac_mask=mhd_trac_mask

    use_particles=mhd_particles
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
       type_divb = divb_none
    {^NOONED
    case ('multigrid')
       if(mhd_radiation_fld) call mpistop('To verify whether mg usage for FLD versus divB can be combined')
       type_divb = divb_multigrid
       use_multigrid = .true.
       mg%operator_type = mg_laplacian
       phys_global_source_after => mhd_clean_divb_multigrid
    }
    case ('glm')
      mhd_glm          = .true.
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
      mhd_glm          = .true.
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
    if (mhd_energy) then
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

    if (mhd_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    if(mhd_hyperbolic_tc) then
      qpar_ = var_set_fluxvar('q', 'q', need_bc=.false.)
      if(mhd_hyperbolic_tc_use_perp) then
        qperp_ = var_set_fluxvar('qperp', 'qperp', need_bc=.false.)
      else
        qperp_ = -1
      end if
      need_global_cmax=.true.
    else
      qpar_ = -1
      qperp_ = -1
    end if

    if (mhd_fip) then
      fip_ = var_set_fluxvar('rho_fip', 'fip', need_bc=.false.)
    else
      fip_ = -1
    end if

    if (mhd_uawsom) then
      wAplus_  = var_set_fluxvar('wAplus',  'wAplus')
      wAminus_ = var_set_fluxvar('wAminus', 'wAminus')
      wkplus_  = var_set_fluxvar('wkplus',  'wkplus')
      wkminus_ = var_set_fluxvar('wkminus', 'wkminus')
    else
      wAplus_=-1; wAminus_=-1; wkplus_=-1; wkminus_=-1
    end if

    allocate(tracer(mhd_n_tracer))
    ! Set starting index of tracers
    do itr = 1, mhd_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do
 
    if(mhd_radiation_fld)then
       if(mhd_cak_force)then
          if(mype==0) then
            write(*,*)'Warning: CAK force addition together with FLD radiation'
          endif
       endif
       if(mhd_radiative_cooling)then
          if(mype==0) then
            write(*,*)'Warning: Optically thin cooling together with FLD radiation'
          endif
       endif
       if(.not.mhd_energy)then
          call mpistop('using FLD implies the use of an energy equation, set mhd_energy=T')
       else
          if(mhd_semirelativistic)then
            call mpistop('using FLD not yet with semirelativistic energy formalism')
          endif
          if(mhd_hydrodynamic_e.or.mhd_internal_e)then
            call mpistop('using FLD not yet with hydrodynamic or internal energy formalism')
          endif
          if(has_equi_rho_and_p)then
            call mpistop('using FLD not yet with split off rho and p')
          endif
          ! Note: so far ok with total energy equation but allow both split or unsplit B0
          !> set added variable and equation for radiation energy
          r_e = var_set_radiation_energy()
          phys_get_csrad2          => mhd_get_csrad2_prim
          !> Radiation fluid object: its EoS callbacks are wired in mhd_link_eos
          allocate(fld_fl)
          !> Initiate radiation-closure module
          call fld_init()
          !> The implicit (MG diffusion) hooks need the fld_fl object, so they
          !> are wired here to physics-module wrappers that inject it.
          if(use_multigrid)then
             phys_implicit_update   => mhd_fld_implicit_update
             phys_evaluate_implicit => mhd_fld_evaluate_implicit
          endif
       endif
    else
      r_e=-1
    endif

    ! LTE/PI electron-density and temperature aux variables. These must be registered
    ! after every flux variable (including the FLD radiation energy r_e and any tracers):
    ! var_set_radiation_energy indexes by nwflux, so if r_e is registered after these
    ! nw-indexed aux vars its nwflux slot collides with Ne_ (iw_r_e==iw_ne) and overwrites
    ! its name. Registering them here keeps the flux block contiguous.
    if (eos%eos_type == 'LTE') then
      Ne_ = var_set_ne()
      Te_ = var_set_te()
    else if (eos%eos_type == 'PI') then
      Ne_ = -1
      Te_ = var_set_te()
    else
      Ne_ = -1
      Te_ = -1
    end if

    ! set number of variables which need update ghostcells
    ! set number of variables which need update ghostcells.
    ! The EoS-derived slots Ne/Te are derived state that must stay consistent with (rho,e)
    ! wherever the conserved state is valid, so they have to travel with it. var_set_ne /
    ! var_set_te bump nw but neither nwflux nor nwaux, so the historic nwflux+nwaux silently
    ! excluded them: they were never communicated, only derived, and every ghost cell held
    ! zero until something derived it. Extend the window to whichever of them exist (they are
    ! registered contiguously just above); FI leaves both at -1 and the window is unchanged.
    ! bc_phys is unaffected -- it iterates nwflux+nwaux independently.
    nwgc=nwflux+nwaux
    if (iw_ne > 0) nwgc = max(nwgc, iw_ne)
    if (iw_te > 0) nwgc = max(nwgc, iw_te)

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    ! set cutoff temperature when using the TRAC method, as well as an auxiliary weight
    Tweight_ = -1
    if(mhd_trac) then
      Tcoff_ = var_set_wextra()
      iw_Tcoff=Tcoff_
      if(mhd_trac_type .ge. 3) then
        Tweight_ = var_set_wextra()
      endif
    else
      Tcoff_ = -1
    end if

    ! set indices of equi vars and update number_equi_vars
    number_equi_vars = 0
    if(has_equi_rho_and_p) then
      number_equi_vars = number_equi_vars + 1
      equi_rho0_ = number_equi_vars
      iw_equi_rho = equi_rho0_
      number_equi_vars = number_equi_vars + 1
      equi_pe0_ = number_equi_vars
      iw_equi_p = equi_pe0_
    endif
    ! determine number of stagger variables
    nws=ndim

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1
    iw_vector(2) = mag(1) - 1

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
      if(mhd_glm) then
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
 
    phys_get_rho             => mhd_get_rho
    phys_get_dt              => mhd_get_dt
    if(mhd_semirelativistic) then
      if(mhd_energy) then
        phys_get_cmax            => mhd_get_cmax_semirelati
      else
        phys_get_cmax            => mhd_get_cmax_semirelati_noe
      end if
    else
      if(mhd_energy) then
        phys_get_cmax            => mhd_get_cmax_origin
      else
        phys_get_cmax            => mhd_get_cmax_origin_noe
      end if
    end if
    phys_get_tcutoff         => mhd_get_tcutoff
    phys_get_H_speed         => mhd_get_H_speed
    if(has_equi_rho_and_p) then
      phys_get_cbounds         => mhd_get_cbounds_split_rho
    else if(mhd_semirelativistic) then
      phys_get_cbounds         => mhd_get_cbounds_semirelati
    else
      phys_get_cbounds         => mhd_get_cbounds
    end if

    !> EOS module: phys_to_primitive / phys_to_conserved are bound by
    !> mod_mhd_eos:bind_eos_to_source to route through the EoS layer;
    !> mhd_to_primitive / mhd_to_conserved pointers are NOT used - every
    !> call to convert state goes through eos%to_primitive / eos%to_conserved.
    !> eos%inv_squared_c{0,} are set in mhd_physical_units (called below)
    !> after inv_squared_c{0,} are computed.
    if(mhd_hydrodynamic_e) then
      phys_get_flux            => mhd_get_flux_hde
    else if(mhd_semirelativistic) then
      if(mhd_energy) then
        phys_get_flux            => mhd_get_flux_semirelati
      else
        phys_get_flux            => mhd_get_flux_semirelati_noe
      end if
    else
      if(B0field.or.has_equi_rho_and_p) then
        phys_get_flux            => mhd_get_flux_split
      else if(mhd_energy) then
        phys_get_flux            => mhd_get_flux
      else
        phys_get_flux            => mhd_get_flux_noe
      end if
    end if
    phys_get_v                 => mhd_get_v
    if(mhd_semirelativistic) then
      phys_add_source_geom     => mhd_add_source_geom_semirelati
    else if(B0field.or.has_equi_rho_and_p) then
      phys_add_source_geom     => mhd_add_source_geom_split
    else
      phys_add_source_geom     => mhd_add_source_geom
    end if
    phys_add_source          => mhd_add_source
    phys_check_params        => mhd_check_params
    phys_write_info          => mhd_write_info
 
    if(mhd_internal_e) then
      phys_handle_small_values => mhd_handle_small_values_inte
      mhd_handle_small_values  => mhd_handle_small_values_inte
      phys_check_w             => mhd_check_w_inte
    else if(mhd_hydrodynamic_e) then
      phys_handle_small_values => mhd_handle_small_values_hde
      mhd_handle_small_values  => mhd_handle_small_values_hde
      phys_check_w             => mhd_check_w_hde
    else if(mhd_semirelativistic) then
      phys_handle_small_values => mhd_handle_small_values_semirelati
      mhd_handle_small_values  => mhd_handle_small_values_semirelati
      phys_check_w             => mhd_check_w_semirelati
    else if(has_equi_rho_and_p) then
      phys_handle_small_values => mhd_handle_small_values_split
      mhd_handle_small_values  => mhd_handle_small_values_split
      phys_check_w             => mhd_check_w_split
    else if(mhd_energy) then
      phys_handle_small_values => mhd_handle_small_values_origin
      mhd_handle_small_values  => mhd_handle_small_values_origin
      phys_check_w             => mhd_check_w_origin
    else
      phys_handle_small_values => mhd_handle_small_values_noe
      mhd_handle_small_values  => mhd_handle_small_values_noe
      phys_check_w             => mhd_check_w_noe
    end if
 
    ! phys_get_pthermal is set by mhd_link_eos

    if(number_equi_vars>0) then
      phys_set_equi_vars => set_equi_vars_grid
    endif

    if(type_divb==divb_glm) then
      phys_modify_wLR => mhd_modify_wLR
    end if

    ! Rfactor / temperature / pthermal pointers are bound by
    ! mod_mhd_eos:bind_eos_to_source (called by mhd_link_eos immediately
    ! after mhd_phys_init). No EoS machinery in mod_mhd_phys.t.

    ! if using ct stagger grid, boundary divb=0 is not done here
    if(stagger_grid) then
      select case(type_ct)
      case('average')
        transverse_ghost_cells = 1
        phys_get_ct_velocity => mhd_get_ct_velocity_average
        phys_update_faces => mhd_update_faces_average
      case('uct_contact')
        transverse_ghost_cells = 1
        phys_get_ct_velocity => mhd_get_ct_velocity_contact
        phys_update_faces => mhd_update_faces_contact
      case('uct_hll')
        transverse_ghost_cells = 2
        phys_get_ct_velocity => mhd_get_ct_velocity_hll
        phys_update_faces => mhd_update_faces_hll
      case default
        call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
      end select
      phys_face_to_center => mhd_face_to_center
      phys_modify_wLR => mhd_modify_wLR
    else if(ndim>1) then
      phys_boundary_adjust => mhd_boundary_adjust
    end if

    {^NOONED
    ! clean initial divb
    if(clean_initial_divb.and.mhd_radiation_fld) &
       call mpistop('To verify whether mg usage for FLD versus divB can be combined')
    if(clean_initial_divb) phys_clean_divb => mhd_clean_divb_multigrid
    }

    ! derive units from basic units
    call mhd_physical_units()

    if(mhd_uawsom) then
      ! Namelist inputs are physical quantities, following the normal AMRVAC
      ! convention for this module (cm/G in cgs, m/T in SI).
      if(mhd_uawsom_zeta_scale == -one) &
        mhd_uawsom_zeta_scale=merge(3.4805d9,3.4805d11,SI_unit)
      if(mhd_uawsom_thread_radius0 == -one) &
        mhd_uawsom_thread_radius0=merge(1.d5,1.d7,SI_unit)
      if(mhd_uawsom_Bref == -one) &
        mhd_uawsom_Bref=merge(1.d-3,10.d0,SI_unit)
      if(mhd_uawsom_alfven_corr_length0 == -one) &
        mhd_uawsom_alfven_corr_length0=merge(1.5d7,1.5d9,SI_unit)
      mhd_uawsom_zeta_scale = mhd_uawsom_zeta_scale/unit_length
      mhd_uawsom_thread_radius0 = mhd_uawsom_thread_radius0/unit_length
      mhd_uawsom_Bref = mhd_uawsom_Bref/unit_magneticfield
      mhd_uawsom_alfven_corr_length0 = &
           mhd_uawsom_alfven_corr_length0/unit_length
    end if

    if(mhd_hyperbolic_tc) then
      if(mhd_hyperbolic_tc_kappa==0.d0) then
        if(SI_unit) then
          mhd_hyperbolic_tc_kappa=8.d-12*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
        else
          mhd_hyperbolic_tc_kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
        end if
      else
        mhd_hyperbolic_tc_constant=.true.
      end if
    
      if(mhd_hyperbolic_tc_use_perp) then
        select case(trim(mhd_hyperbolic_tc_perp_mode))
        case('fixed_reference')
          if(mhd_hyperbolic_tc_kappa_perp_factor==0.d0) then
            if(SI_unit) then
              mhd_hyperbolic_tc_kappa_perp_factor = 4.d-30 * &
                unit_numberdensity**2 / unit_magneticfield**2 / unit_temperature**3
            else
              mhd_hyperbolic_tc_kappa_perp_factor = 4.d-10 * &
                unit_numberdensity**2 / unit_magneticfield**2 / unit_temperature**3
            end if
          end if
        case('weak_field_isotropization')
          if(mhd_hyperbolic_tc_Bmin==0.d0) then
            mhd_hyperbolic_tc_Bmin=0.1d0/unit_magneticfield
          end if
        case('electron_magnetization')
          if(mhd_hyperbolic_tc_coulomb_log<=zero) then
            call mpistop("mhd_hyperbolic_tc_coulomb_log must be positive")
          end if
        case default
          call mpistop("invalid mhd_hyperbolic_tc_perp_mode: "// &
               trim(mhd_hyperbolic_tc_perp_mode))
        end select
      end if
    end if
    if(.not. mhd_energy .and. mhd_thermal_conduction) then
      call mpistop("thermal conduction needs mhd_energy=T")
    end if
    if(.not. mhd_energy .and. mhd_hyperbolic_tc) then
      call mpistop("hyperbolic thermal conduction needs mhd_energy=T")
    end if
    if(.not. mhd_energy .and. mhd_radiative_cooling) then
      call mpistop("radiative cooling needs mhd_energy=T")
    end if

    !> Cache log10(nH) in wextra for LTE+IonE TC (density invariant during STS)
    if (eos%eos_type == 'LTE' .and. eos%ionE .and. mhd_thermal_conduction) then
        iw_log_nH = var_set_wextra()
    end if

    if(mhd_equi_thermal)then
       if((.not.has_equi_rho_and_p).or.(.not.total_energy))then
          mhd_equi_thermal=.false.
          if(mype==0) write(*,*) 'WARNING: turning mhd_equi_thermal=F as no splitting or total e in use'
       else
          if(mhd_thermal_conduction.or.mhd_radiative_cooling)then
            if(mype==0) write(*,*) 'Will subtract thermal balance in TC or RC with mhd_equi_thermal=T'
          else
            mhd_equi_thermal=.false.
            if(mype==0) write(*,*) 'WARNING: turning mhd_equi_thermal=F as no TC or RC in use'
          endif
       endif
    endif

    ! initialize thermal conduction module
    if (mhd_thermal_conduction) then
      call sts_init()
      call tc_init_params(eos%gamma)

      allocate(tc_fl)
      call tc_get_mhd_params(tc_fl,tc_params_read_mhd)
      if(ndim==1) then
        call add_sts_method(mhd_get_tc_dt_hd,mhd_sts_set_source_tc_hd,e_,1,e_,1,.false.)
      else
        call add_sts_method(mhd_get_tc_dt_mhd,mhd_sts_set_source_tc_mhd,e_,1,e_,1,.false.)
      endif
      ! TC function pointers (get_temperature_from_conserved/eint, get_rho,
      ! get_temperature_equi, get_rho_equi, subtract_equi) are bound by
      ! bind_eos_to_source in mod_mhd_eos.t to the correct EoS-aware
      ! implementations. No EoS machinery in mod_mhd_phys.t.
      if(.not.mhd_internal_e) then
        if(mhd_hydrodynamic_e) then
          call set_conversion_methods_to_head(mhd_e_to_ei_hde, mhd_ei_to_e_hde)
          phys_e_to_ei => mhd_e_to_ei_hde
          phys_ei_to_e => mhd_ei_to_e_hde
        else if(mhd_semirelativistic) then
          call set_conversion_methods_to_head(mhd_e_to_ei_semirelati, mhd_ei_to_e_semirelati)
          phys_e_to_ei => mhd_e_to_ei_semirelati
          phys_ei_to_e => mhd_ei_to_e_semirelati
        else
          if (iw_log_nH > 0) then
            call set_conversion_methods_to_head(mhd_e_to_ei_and_cache_log_nH, mhd_ei_to_e)
          else
            call set_conversion_methods_to_head(mhd_e_to_ei, mhd_ei_to_e)
          end if
          phys_e_to_ei => mhd_e_to_ei
          phys_ei_to_e => mhd_ei_to_e
        end if
      end if
      call set_error_handling_to_head(mhd_tc_handle_small_e)
      tc_fl%e_ = e_
      tc_fl%Tcoff_ = Tcoff_
    end if

    ! Energy conversion pointers needed by EOS module regardless of TC method
    if(.not.mhd_internal_e .and. .not.associated(phys_e_to_ei)) then
      if(mhd_hydrodynamic_e) then
        phys_e_to_ei => mhd_e_to_ei_hde
        phys_ei_to_e => mhd_ei_to_e_hde
      else if(mhd_semirelativistic) then
        phys_e_to_ei => mhd_e_to_ei_semirelati
        phys_ei_to_e => mhd_ei_to_e_semirelati
      else
        phys_e_to_ei => mhd_e_to_ei
        phys_ei_to_e => mhd_ei_to_e
      end if
    end if

    ! Initialize radiative cooling module
    if (mhd_radiative_cooling) then
      call radiative_cooling_init_params(eos%gamma,eos%He_abundance)
      allocate(rc_fl)
      rc_fl%fip_ = fip_
      call radiative_cooling_init(rc_fl,rc_params_read)
      rc_fl%e_ = e_
      rc_fl%Tcoff_ = Tcoff_
      rc_fl%subtract_equi = has_equi_rho_and_p
      ! rc_fl EoS pointers (get_rho, get_pthermal, get_var_Rfactor,
      ! equi accessors, subtract_equi) are bound by bind_eos_to_source
      ! in mod_mhd_eos.t. No EoS machinery in mod_mhd_phys.t.
    end if
    allocate(te_fl_mhd)
    ! te_fl_mhd EoS pointers are bound by bind_eos_to_source in mod_mhd_eos.t
{^IFTHREED
    phys_te_images => mhd_te_images
}

    ! consistency check for hyperresistivity implementation
    if (mhd_eta_hyper>0.0d0) then
      if(mype==0) then
         write(*,*) '*****Using hyperresistivity:  with mhd_eta_hyper :', mhd_eta_hyper
      endif
       if(B0field) then
        ! hyperresistivity not ok yet with splitting
        call mpistop("Must have B0field=F when using hyperresistivity")
      end if
    endif
    if (mhd_eta_hyper<0.0d0) then
        call mpistop("Must have mhd_eta_hyper positive when using hyperresistivity")
    endif

    ! Initialize viscosity module
    if (mhd_viscosity) then
       call viscosity_init(phys_wider_stencil)
    end if

    ! Initialize gravity module
    if(mhd_gravity) then
      call gravity_init()
    end if

    ! Initialize rotating frame module
    if(mhd_rotating_frame) then
      if(has_equi_rho_and_p) then
        ! mod_rotating_frame does not handle splitting of density
        call mpistop("Must have has_equi_rho_and_p=F when mhd_rotating_frame=T")
      end if
      call rotating_frame_init()
    endif


    ! initialize magnetofriction module
    if(mhd_magnetofriction) then
      call magnetofriction_init()
    end if

    if(mhd_hall) then
      if(mhd_semirelativistic) then
        ! semirelativistic does not incorporate hall terms
        call mpistop("Must have mhd_hall=F when mhd_semirelativistic=T")
      end if
      if(coordinate>1)then
        ! normal unsplit case or split cases do not have geometric sources for Hall included
        call mpistop("Must have Cartesian coordinates for Hall")
      endif
      ! For Hall, we need one more reconstructed layer since currents are computed
      ! in mhd_get_flux: assuming one additional ghost layer added in nghostcells.
      phys_wider_stencil = 1
    end if

    ! The perpendicular HTC geometry evaluates a fourth-order centred
    ! temperature gradient inside the reconstructed flux layer. It therefore
    ! needs one layer beyond the default two ghost cells.
    if(mhd_hyperbolic_tc .and. mhd_hyperbolic_tc_use_perp) then
      phys_wider_stencil=max(phys_wider_stencil,1)
    end if

    if(mhd_ambipolar) then
      if(mhd_ambipolar_sts) then
        call sts_init()
        if(mhd_internal_e.or.mhd_hydrodynamic_e) then
          call add_sts_method(get_ambipolar_dt,sts_set_source_ambipolar,mag(1),&
               ndir,mag(1),ndir,.true.)
        else
          ! any total energy or no energy at all case is handled here
          call add_sts_method(get_ambipolar_dt,sts_set_source_ambipolar,mom(ndir)+1,&
               mag(ndir)-mom(ndir),mag(1),ndir,.true.)
        end if
      else
        mhd_ambipolar_exp=.true.
        ! For flux ambipolar term, we need one more reconstructed layer since currents are computed
        ! in mhd_get_flux: assuming one additional ghost layer added in nghostcells.
        phys_wider_stencil = 1
      end if
    end if

    ! ionization-degree table init now lives in eos_finalise (eos% owns
    ! thermodynamic-backend init, parallel to LTE tables); see mod_eos_PI.

    ! Initialize CAK radiation force module
    if (mhd_cak_force) then
       if(mhd_internal_e.or.mhd_semirelativistic) then
          call mpistop("CAK implementation not available in internal or semirelativistic variants")
       endif
       if(has_equi_rho_and_p) then
          call mpistop("CAK force implementation not available for split off pressure and density")
       endif
       call cak_init(eos%gamma)
    endif

  end subroutine mhd_phys_init

{^IFTHREED
  subroutine mhd_te_images
    use mod_global_parameters
    use mod_thermal_emission

    select case(convert_type)
      case('EIvtiCCmpi','EIvtuCCmpi')
        call get_EUV_image(unitconvert,te_fl_mhd)
      case('ESvtiCCmpi','ESvtuCCmpi')
        call get_EUV_spectrum(unitconvert,te_fl_mhd)
      case('SIvtiCCmpi','SIvtuCCmpi')
        call get_SXR_image(unitconvert,te_fl_mhd)
      case('WIvtiCCmpi','WIvtuCCmpi')
        call get_whitelight_image(unitconvert,te_fl_mhd)
      case default
        call mpistop("Error in synthesize emission: Unknown convert_type")
      end select
  end subroutine mhd_te_images
}

!!start th cond
  ! wrappers for STS functions in thermal_conductivity module
  ! which take as argument the tc_fluid (defined in the physics module)
  subroutine  mhd_sts_set_source_tc_mhd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_mhd
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_mhd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,tc_fl)
  end subroutine mhd_sts_set_source_tc_mhd

  subroutine  mhd_sts_set_source_tc_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_hd
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,tc_fl)
  end subroutine mhd_sts_set_source_tc_hd

  function mhd_get_tc_dt_mhd(w,ixI^L,ixO^L,dx^D,x) result(dtnew)
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
  end function mhd_get_tc_dt_mhd

  function mhd_get_tc_dt_hd(w,ixI^L,ixO^L,dx^D,x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_hd
 
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_hd(w,ixI^L,ixO^L,dx^D,x,tc_fl)
  end function mhd_get_tc_dt_hd

  subroutine mhd_tc_handle_small_e(w, x, ixI^L, ixO^L, step)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer, intent(in)    :: step
    character(len=140) :: error_msg

    ! tc_patch_eint repairs w(:, e_) in place during the source-term call,
    ! but the Chebyshev recurrence can still produce residual negative
    ! e_int in the updated state.  Run mhd_handle_small_ei every substep
    ! (patch on or off) so STS does not propagate negative e through.
    write(error_msg,"(a,i3)") "Thermal conduction step ", step
    call mhd_handle_small_ei(w,x,ixI^L,ixO^L,e_,error_msg)
  end subroutine mhd_tc_handle_small_e

  ! fill in tc_fluid fields from namelist
  subroutine tc_params_read_mhd(fl)
    use mod_global_parameters, only: unitpar,par_files,unit_temperature
    type(tc_fluid), intent(inout) :: fl

    double precision :: tc_k_para=0d0
    double precision :: tc_k_perp=0d0
    integer                      :: n
    ! list parameters
    logical :: tc_perpendicular=.false.
    logical :: tc_saturate=.false.
    logical :: tc_patch_eint=.false.
    double precision :: trac_T_floor=0.d0
    character(len=std_len)  :: tc_slope_limiter="MC"

    namelist /tc_list/ tc_perpendicular, tc_saturate, tc_slope_limiter, tc_k_para, tc_k_perp, tc_patch_eint, trac_T_floor

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, tc_list, end=111)
111     close(unitpar)
    end do

    fl%tc_perpendicular = tc_perpendicular
    fl%tc_saturate = tc_saturate
    fl%tc_patch_eint = tc_patch_eint
    fl%tc_k_para = tc_k_para
    fl%tc_k_perp = tc_k_perp
    fl%trac_T_floor = trac_T_floor / unit_temperature
    select case(tc_slope_limiter)
     case ('no','none')
       fl%tc_slope_limiter = 0
     case ('MC')
       ! monotonized central limiter Woodward and Collela limiter (eq.3.51h)
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
     case ('vanleer')
       ! van Leer limiter
       fl%tc_slope_limiter = 5
     case default
       call mpistop("Unknown tc_slope_limiter, choose MC, minmod, superbee, koren, vanleer")
    end select
  end subroutine tc_params_read_mhd
!!end th cond

!!rad cool
  subroutine rc_params_read(fl)
    use mod_global_parameters, only: unitpar,par_files
    use mod_constants, only: bigdouble
    type(rc_fluid), intent(inout) :: fl

    !> Lower limit of temperature
    double precision :: tlow=bigdouble
    double precision :: rad_damp_height=0.5d0
    double precision :: rad_damp_scale=0.15d0
    integer :: n
    ! list parameters
    integer :: ncool = 4000
    !> Fixed temperature not lower than tlow
    logical :: Tfix=.false.
    !> Add cooling source in a split way (.true.) or un-split way (.false.)
    logical :: rc_split=.false.
    logical :: rad_damp=.false.
    !> Name of cooling curve
    character(len=std_len)  :: coolcurve='JCcorona'
    logical :: rad_newton = .false.
    double precision :: rad_newton_trad = 0.006d0
    double precision :: rad_newton_rhosurf = 1.d4
    double precision :: rad_newton_pthick = 25.d0
    !> HEAD-side cooling parameters (missing from common content after merge)
    double precision :: cfrac=0.1d0
    double precision :: rad_cut_hgt=0.5d0
    double precision :: rad_cut_dey=0.15d0
    !> Variable-c_V Townsend extension (Y_mod): quadrature and sub-intervals
    character(len=8) :: rc_Y_mod_quadrature='boole'
    integer :: rc_Y_mod_N_sub=16
    !> Smooth SC<->thin blend (mirrors hd): thin cooling *= 0.5(1+tanh((T-T0)/dT))*0.5(1+tanh((y-y0)/dy))

    namelist /rc_list/ coolcurve, ncool, cfrac, tlow, Tfix, rc_split, &
                       rad_cut_hgt, rad_cut_dey, &
                       rc_Y_mod_quadrature, rc_Y_mod_N_sub, &
                       rad_newton, rad_newton_trad, rad_newton_rhosurf, &
                       rad_newton_pthick, rad_damp, rad_damp_height, rad_damp_scale

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, rc_list, end=111)
111     close(unitpar)
    end do

    fl%ncool=ncool
    fl%coolcurve=coolcurve
    fl%tlow=tlow
    fl%Tfix=Tfix
    fl%rc_split=rc_split
    fl%cfrac=cfrac
    fl%rad_cut_hgt=rad_cut_hgt
    fl%rad_cut_dey=rad_cut_dey
    fl%Y_mod_quadrature=rc_Y_mod_quadrature
    fl%Y_mod_N_sub=rc_Y_mod_N_sub
    fl%rad_damp=rad_damp
    fl%rad_damp_height=rad_damp_height
    fl%rad_damp_scale=rad_damp_scale
    fl%rad_newton=rad_newton
    fl%rad_newton_trad=rad_newton_trad
    fl%rad_newton_rhosurf=rad_newton_rhosurf
    fl%rad_newton_pthick=rad_newton_pthick
  end subroutine rc_params_read

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
    integer, intent(in)             :: ixI^L,ixO^L, nwc
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision   :: wnew(ixO^S, 1:nwc)

    if(has_equi_rho_and_p) then
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

    if(mhd_energy) then
      wnew(ixO^S,e_)=w(ixO^S,e_)
      if(has_equi_rho_and_p) then
        wnew(ixO^S,e_)=wnew(ixO^S,e_)+block%equi_vars(ixO^S,equi_pe0_,0)*eos%inv_gamma_minus_1
      end if
      if(B0field .and. total_energy) then
        wnew(ixO^S,e_)=wnew(ixO^S,e_)+0.5d0*sum(block%B0(ixO^S,:,0)**2,dim=ndim+1) &
            + sum(w(ixO^S,mag(:))*block%B0(ixO^S,:,0),dim=ndim+1)
      end if
    end if

  end function convert_vars_splitting

  subroutine mhd_check_params
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry, only: coordinate, Cartesian
    use mod_convert, only: add_convert_method
    use mod_particles, only: particles_init, particles_eta, particles_etah
    use mod_particles, only: npayload,nusrpayload, &
                ngridvars,num_particles,physics_type_particles
    use mod_fld

    double precision :: a,b,Xfrac,Yfrac

    if(mhd_uawsom) then
      if(coordinate /= Cartesian .or. .not.slab_uniform) &
        call mpistop('mhd_uawsom currently requires a uniform Cartesian grid')
      if(.not.mhd_energy .or. .not.total_energy) &
        call mpistop('mhd_uawsom requires the standard total-energy MHD formulation')
      if(mhd_internal_e .or. mhd_hydrodynamic_e .or. mhd_semirelativistic .or. &
         has_equi_rho_and_p .or. mhd_radiation_fld) &
        call mpistop('mhd_uawsom: unsupported MHD option')
      if(trim(eos%eos_type) /= 'FI') &
        call mpistop("mhd_uawsom currently supports eos_type='FI' only")
      if(mhd_uawsom_height_dim < 1 .or. mhd_uawsom_height_dim > ndim) &
        call mpistop('mhd_uawsom_height_dim must be between 1 and ndim')
      if(mhd_uawsom_zeta0 <= one .or. mhd_uawsom_filling_factor <= zero .or. &
         mhd_uawsom_filling_factor >= one) &
        call mpistop('mhd_uawsom requires zeta0>1 and 0<filling_factor<1')
      if(mhd_uawsom_zeta_scale <= zero .or. &
         mhd_uawsom_thread_radius0 <= zero .or. &
         mhd_uawsom_alfven_corr_length0 <= zero .or. mhd_uawsom_Bref <= zero) &
        call mpistop('mhd_uawsom: scales must be positive')
      if(mhd_uawsom_sigma < zero) &
        call mpistop('mhd_uawsom_sigma must be non-negative')
      select case(trim(mhd_uawsom_reflection_mode))
      case('one_dimensional_gradient','cartesian_gradient_vorticity')
        continue
      case default
        call mpistop('invalid mhd_uawsom_reflection_mode')
      end select
    end if

    ! Initialize particles module here, so all extra and user vars are sample
    if(mhd_particles) then
      call particles_init()
      if (particles_eta  < zero) particles_eta = mhd_eta
      if (particles_etah < zero) particles_eta = mhd_etah
    end if

    ! gamma, gamma_minus_1, inv_gamma_minus_1 are set by eos_init
    if (.not. mhd_energy) then
       if (eos%gamma <= 0.0d0) call mpistop ("Error: gamma <= 0")
       if (mhd_adiab < 0.0d0) call mpistop ("Error: mhd_adiab < 0")
       small_pressure = mhd_adiab*small_density**eos%gamma
    else
       if (eos%gamma <= 0.0d0 .or. eos%gamma == 1.0d0) &
            call mpistop ("Error: gamma <= 0 or gamma == 1")
       small_e = small_pressure * eos%inv_gamma_minus_1
       small_r_e = small_pressure * eos%inv_gamma_minus_1
    end if

    if (number_equi_vars > 0 .and. .not. associated(usr_set_equi_vars)) then
      call mpistop("usr_set_equi_vars has to be implemented in the user file")
    endif

    if(has_equi_rho_and_p) then
      ! The Roe and HLLC solvers read the rho slot as a physical density, which
      ! under splitting holds only the perturbation. Neither has a split variant,
      ! so refuse the combination rather than return a quietly wrong wave speed.
      if(any(flux_method(:)==fs_tvdmu) .or. any(flux_method(:)==fs_hllc) &
         .or. any(flux_method(:)==fs_hllcd)) then
        call mpistop("Must have has_equi_rho_and_p=F with the roe/hllc/hllcd flux schemes")
      end if
      ! mhd_get_tcutoff forms Te from the perturbation pressure and density.
      if(mhd_trac) then
        call mpistop("Must have has_equi_rho_and_p=F when mhd_trac=T")
      end if
    end if
    if(convert .or. autoconvert) then
      if(convert_type .eq. 'dat_generic_mpi') then
        if(mhd_dump_full_vars) then
          if(mype .eq. 0) print*, " add conversion method: split -> full "
          call add_convert_method(convert_vars_splitting, nw, cons_wnames, "new")
        endif
      endif
    endif

    if(mhd_radiation_fld) then 
        if(.not.use_imex_scheme)then
           call mpistop('select IMEX scheme for FLD radiation use')
        endif
        if(use_multigrid)then
           call phys_set_mg_bounds()
        else
           if(.not.fld_no_mg)call mpistop('multigrid must have BCs for IMEX and FLD radiation use')
        endif
        if(mype==0)then
           write(*,*)'==FLD SETUP======================'
           write(*,*)'Using FLD with settings:'
           write(*,*)'Using FLD with settings: mhd_radiation_fld=',mhd_radiation_fld
           write(*,*)'Using FLD with settings: fld_fluxlimiter=',fld_fluxlimiter
           write(*,*)'Using FLD with settings: fld_interaction_method=',fld_interaction_method
           write(*,*)'Using FLD with settings: fld_opacity_law=',fld_opacity_law
           write(*,*)'Using FLD with settings: fld_kappa0=',fld_kappa0
           write(*,*)'Using FLD with settings: fld_opal_table=',fld_opal_table
           write(*,*)'Using FLD with settings: fld_Radforce_split=',fld_Radforce_split
           write(*,*)'Using FLD with settings: fld_bisect_tol=',fld_bisect_tol
           write(*,*)'Using FLD with settings: fld_diff_tol=',fld_diff_tol
           write(*,*)'Using FLD with settings: nth_for_diff_mg=',nth_for_diff_mg
           write(*,*)'      FLD has use_imex_scheme and use_multigrid=',use_imex_scheme,use_multigrid
           print *,'const_rad_a   =',const_rad_a
           print *,'NORMALIZED arad_norm=',arad_norm
           print *,'NORMALIZED c_norm=',c_norm
           print *,'const_kappae  =',const_kappae 
           if(trim(fld_opacity_law).eq.'const_norm')then
               print *,'NORMALIZED fld_kappa0          =',fld_kappa0
               print *,'physical value (in cgs or SI)  =',fld_kappa0*unit_opacity
           endif
           if(trim(fld_opacity_law).eq.'const')then
               print *,'physical fld_kappa (in cgs or SI) =',fld_kappa0
               print *,'NORMALIZED value                  =',fld_kappa0/unit_opacity
           endif
           write(*,*)'===FLD SETUP====================='
        endif
    endif

    if(mype==0)then
           write(*,*)'====MHD run with settings===================='
           write(*,*)'Using mod_mhd_phys with settings:'
           write(*,*)'SI_unit=',SI_unit
           write(*,*)'Dimensionality   :',ndim
           write(*,*)'vector components:',ndir
           write(*,*)'coordinate set to type,slab:',coordinate,slab
           write(*,*)'number of variables          nw=',nw
           write(*,*)'    start index         iwstart=',iwstart
           write(*,*)'number of      vector variables=',nvector
           write(*,*)'number of stagger variables nws=',nws
           write(*,*)'number of    variables with BCs=',nwgc
           write(*,*)'number of      vars with fluxes=',nwflux
           write(*,*)'number of   vars with flux + BC=',nwfluxbc
           write(*,*)'number of   auxiliary variables=',nwaux
           write(*,*)'number of extra vars without flux=',nwextra
           write(*,*)'number of extra vars   for wextra=',nw_extra
           write(*,*)'number of auxiliary I/O variables=',nwauxio
           write(*,*)'number of            mhd_n_tracer=',mhd_n_tracer
           write(*,*)'    mhd_energy=',mhd_energy,' with total_energy=',total_energy
           write(*,*)'    mhd_semirelativistic=',mhd_semirelativistic
           write(*,*)'    mhd_internal_e=',mhd_internal_e
           write(*,*)'    mhd_hydrodynamic_e=',mhd_hydrodynamic_e
           write(*,*)'    mhd_gravity=',mhd_gravity
           write(*,*)'    mhd_eta=',mhd_eta,' nonzero implies resistivity'
           write(*,*)'    mhd_viscosity=',mhd_viscosity
           write(*,*)'    mhd_radiative_cooling=',mhd_radiative_cooling
           write(*,*)'    mhd_cak_force=',mhd_cak_force
           write(*,*)'    mhd_radiation_fld=',mhd_radiation_fld
           write(*,*)'    mhd_thermal_conduction=',mhd_thermal_conduction
           write(*,*)'    mhd_hyperbolic_tc=',mhd_hyperbolic_tc
           write(*,*)'    mhd_trac=',mhd_trac
           write(*,*)'    mhd_hall=',mhd_hall
           write(*,*)'    mhd_ambipolar=',mhd_ambipolar
           write(*,*)'    mhd_eta_hyper=',mhd_eta_hyper
           write(*,*)'    mhd_rotating_frame=',mhd_rotating_frame
           write(*,*)'    mhd_particles=',mhd_particles
           if(mhd_particles) then
              write(*,*) '*****Using particles:        with mhd_eta, mhd_etah :', mhd_eta, mhd_etah
              write(*,*) '*****Using particles: particles_eta, particles_etah :', particles_eta, particles_etah
              write(*,*) '*****Using particles: npayload,ngridvars :', npayload,ngridvars
              write(*,*) '*****Using particles:        nusrpayload :', nusrpayload
              write(*,*) '*****Using particles:      num_particles :', num_particles
              write(*,*) '*****Using particles: physics_type_particles=',physics_type_particles
           end if
           write(*,*)'number of             ghostcells=',nghostcells
           write(*,*)'number due to phys_wider_stencil=',phys_wider_stencil
           write(*,*)'==========================================='
           print *,'========EOS and UNITS==========='
           print *,'SI_unit       =',SI_unit
           print *,'gamma=',eos%gamma
           print *,'He_abundance  =',eos%He_abundance
           print *,'RR            =',RR
           print *,'========EOS and UNITS==========='
           print *,'unit_time          =',unit_time
           print *,'unit_length        =',unit_length
           print *,'unit_velocity      =',unit_velocity
           print *,'unit_pressure      =',unit_pressure
           print *,'unit_numberdensity =',unit_numberdensity
           print *,'unit_density       =',unit_density
           print *,'unit_temperature   =',unit_temperature
           print *,'unit_mass          =',unit_mass
           print *,'unit_Erad          =',unit_Erad
           print *,'unit_radflux       =',unit_radflux
           print *,'unit_magneticfield =',unit_magneticfield
           if(SI_unit)then
              print *,'CHECK that p_u',unit_pressure,' equals ',unit_magneticfield**2/miu0_SI
           else
              print *,'CHECK that p_u',unit_pressure,' equals ',unit_magneticfield**2/(4.0d0*dpi)
           endif
           print *, 'CHECK that p_u ',unit_pressure,' equals ',unit_density*unit_velocity**2
           print *, 'CHECK that L_u ',unit_length,' equals ',unit_velocity*unit_time
           print *, 'CHECK that M_u',unit_mass,' equals ',unit_density*unit_length**3
           print *, 'density to numberdensity has factor   ',unit_density/unit_numberdensity
           if(SI_unit)then
                print *, '                     compare  this to ',mp_SI*(1.d0+4.d0*eos%He_abundance)
           else
                print *, '                     compare  this to ',mp_cgs*(1.d0+4.d0*eos%He_abundance)
           endif
           print *, 'pressure to n T has factor            ',unit_pressure/(unit_numberdensity*unit_temperature)
           if(SI_unit)then
                print *, '                     compare  this to ',kB_SI*(2.d0+3.d0*eos%He_abundance)
                a=unit_density/unit_numberdensity/mp_SI
                b=unit_pressure/(unit_numberdensity*unit_temperature*kB_SI)
           else
                print *, '                     compare  this to ',kB_cgs*(2.d0+3.d0*eos%He_abundance)
                a=unit_density/unit_numberdensity/mp_cgs
                b=unit_pressure/(unit_numberdensity*unit_temperature*kB_cgs)
           endif
           if(eos%eos_type /= 'LTE')then
              print *, 'mean molecular weight mu is =',a/b,' = ', (1.d0+4.d0*eos%He_abundance)/(2.d0+3.d0*eos%He_abundance)
              Xfrac=1.d0/a
              Yfrac=4.d0*eos%He_abundance/(1.d0+4.d0*eos%He_abundance)
              print *, 'mass fraction hydrogen X is =',1/a,' and this equals ', 1.d0/(1.d0+4.d0*eos%He_abundance)
              print *, 'mass fraction helium   Y is =',Yfrac
              print *, ' check that 1/mu', b/a,' is equal to 2X+3Y/4=',2.d0*Xfrac+3.d0*Yfrac/4.d0
              print *, ' ratio n_e/n_p=',1.d0+2.0d0*eos%He_abundance
           endif
           print *,'========UNITS==========='
    endif

  end subroutine mhd_check_params

  subroutine mhd_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0,c_lightspeed,Xfrac,sigma_Telectron
    double precision :: a,b
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
      const_sigmaSB=sigma_SB_SI
      c_lightspeed=c_SI
      sigma_Telectron=sigma_Te_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi ! G^2 cm^2 dyne^-1
      const_sigmaSB=sigma_SB_cgs
      c_lightspeed=const_c
      sigma_Telectron=sigma_Te_cgs
    end if
    ! Normalisation dispatch keyed solely on eos%eos_type (FI is the default, so
    ! legacy parfiles land in the FI/PI absorbed-(a,b), RR=1 branch -- the former
    ! eq_state_units=.true. result).
    if (eos%eos_type == 'LTE') then
      !> Remove the assumed FI normalisation from the units and handle in EoS
      a=1d0
      b=1d0
      eos%nH2rhoFactor = 1d0+4d0*eos%He_abundance
      RR=(2d0+3d0*eos%He_abundance) / (1d0+4d0*eos%He_abundance)
      Xfrac=1.d0/(1.d0+4.d0*eos%He_abundance)
    else
      !> FI / PI: absorbed-(a,b), RR=1 (a=b=1 with RR=1 would be wrong physics
      !> for He>0). PI shares FI's normalisation exactly; the partial b lives here.
      a=1d0+4d0*eos%He_abundance
      if(eos%eos_type=='PI') then
        b=1d0+H_ion_fr+eos%He_abundance*(He_ion_fr*(He_ion_fr2+1d0)+1d0)
      else
        b=2d0+3d0*eos%He_abundance
      end if
      RR=1d0
      Xfrac=1.d0/a
    end if
    if(unit_density/=1.d0 .or. unit_numberdensity/=1.d0) then
      if(unit_density/=1.d0) then
        unit_numberdensity=unit_density/(a*mp)
      else if(unit_numberdensity/=1.d0) then
        unit_density=a*mp*unit_numberdensity
      end if
      if(unit_temperature/=1.d0) then
        unit_pressure=b*unit_numberdensity*kB*unit_temperature
        unit_velocity=sqrt(unit_pressure/unit_density)
        unit_magneticfield=sqrt(miu0*unit_pressure)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_magneticfield/=1.d0) then
        unit_pressure=unit_magneticfield**2/miu0
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        unit_velocity=sqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_pressure/=1.d0) then
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        unit_velocity=sqrt(unit_pressure/unit_density)
        unit_magneticfield=sqrt(miu0*unit_pressure)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_velocity/=1.d0) then
        unit_pressure=unit_density*unit_velocity**2
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        unit_magneticfield=sqrt(miu0*unit_pressure)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=1.d0) then
        unit_velocity=unit_length/unit_time
        unit_pressure=unit_density*unit_velocity**2
        unit_magneticfield=sqrt(miu0*unit_pressure)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    else if(unit_temperature/=1.d0) then
      ! units of temperature and velocity are dependent
      if(unit_magneticfield/=1.d0) then
        unit_pressure=unit_magneticfield**2/miu0
        unit_numberdensity=unit_pressure/(b*unit_temperature*kB)
        unit_density=a*mp*unit_numberdensity
        unit_velocity=sqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_pressure/=1.d0) then
        unit_magneticfield=sqrt(miu0*unit_pressure)
        unit_numberdensity=unit_pressure/(b*unit_temperature*kB)
        unit_density=a*mp*unit_numberdensity
        unit_velocity=sqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      end if
    else if(unit_magneticfield/=1.d0) then
      ! units of magnetic field and pressure are dependent
      if(unit_velocity/=1.d0) then
        unit_pressure=unit_magneticfield**2/miu0
        unit_numberdensity=unit_pressure/(b*unit_temperature*kB)
        unit_density=a*mp*unit_numberdensity
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=0.d0) then
        unit_pressure=unit_magneticfield**2/miu0
        unit_velocity=unit_length/unit_time
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    else if(unit_pressure/=1.d0) then
      if(unit_velocity/=1.d0) then
        unit_magneticfield=sqrt(miu0*unit_pressure)
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=0.d0) then
        unit_magneticfield=sqrt(miu0*unit_pressure)
        unit_velocity=unit_length/unit_time
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    end if
    ! Additional units needed for the particles
    c_norm=c_lightspeed/unit_velocity
    unit_charge=unit_magneticfield*unit_length**2/unit_velocity/miu0
    if (.not. SI_unit) unit_charge = unit_charge*const_c
    unit_mass=unit_density*unit_length**3

    if(mhd_semirelativistic) then
      if(mhd_reduced_c<1.d0) then
        ! dimensionless speed
        inv_squared_c0=1.d0
        inv_squared_c=1.d0/mhd_reduced_c**2
      else
        inv_squared_c0=(unit_velocity/c_lightspeed)**2
        inv_squared_c=(unit_velocity/mhd_reduced_c)**2
      end if
      ! Propagate to the EoS container. Must happen AFTER inv_squared_c{0,}
      ! are set above; the assignment earlier in mhd_phys_init runs before
      ! mhd_physical_units and would store uninitialised values.
      eos%inv_squared_c0 = inv_squared_c0
      eos%inv_squared_c  = inv_squared_c
    end if

    !> Units for radiative flux and opacity as used in FLD
    ! this is the radiation constant in either cgs or SI units
    const_rad_a=4.d0*const_sigmaSB/c_lightspeed
    ! this is the dimensionless conversion factor for Erad to Trad
    arad_norm=const_rad_a*unit_temperature**4/unit_pressure
    ! This is the Thomson scattering opacity in the correct units
    ! note that the hydrogen mass fraction X=1/a in eq_state_units
    const_kappae=sigma_Telectron*(1.d0+Xfrac)/(2.0d0*mp)
    ! these are the units
    unit_Erad = unit_pressure
    unit_radflux = unit_velocity*unit_Erad
    unit_opacity = one/(unit_density*unit_length)

  end subroutine mhd_physical_units

  subroutine mhd_check_w_semirelati(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    logical, intent(inout) :: flag(ixI^S,1:nw)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)

    double precision :: tmp,b(1:ndir),v(1:ndir),factor
    integer :: ix^D

    flag=.false.
    where(w(ixO^S,rho_) < small_density) flag(ixO^S,rho_) = .true.

    if(mhd_energy) then
      if(primitive) then
        where(w(ixO^S,p_) < small_pressure) flag(ixO^S,e_) = .true.
      else
        if(mhd_internal_e) then
         {do ix^DB=ixOmin^DB,ixOmax^DB \}
            if(w(ix^D,e_) < small_e) flag(ix^D,e_) = .true.
         {end do\}
        else
         {do ix^DB=ixOmin^DB,ixOmax^DB \}
            ! Convert momentum to velocity
            tmp=(^C&w(ix^D,b^C_)*w(ix^D,m^C_)+)*inv_squared_c
            factor=1.0d0/(w(ix^D,rho_)*(w(ix^D,rho_)+(^C&w(ix^D,b^C_)**2+)*inv_squared_c))
            ^C&v(^C)=factor*(w(ix^D,m^C_)*w(ix^D,rho_)+w(ix^D,b^C_)*tmp)\
            ! E=Bxv
            {^IFTHREEC
            b(1)=w(ix^D,b2_)*v(3)-w(ix^D,b3_)*v(2)
            b(2)=w(ix^D,b3_)*v(1)-w(ix^D,b1_)*v(3)
            b(3)=w(ix^D,b1_)*v(2)-w(ix^D,b2_)*v(1)
            }
            {^IFTWOC
            b(1)=zero
            ! switch 3 with 2 to allow ^C from 1 to 2
            b(2)=w(ix^D,b1_)*v(2)-w(ix^D,b2_)*v(1)
            }
            {^IFONEC
            b(1)=zero
            }
            ! Calculate internal e = e-eK-eB-eE
            tmp=w(ix^D,e_)-half*((^C&v(^C)**2+)*w(ix^D,rho_)&
               +(^C&w(ix^D,b^C_)**2+)+(^C&b(^C)**2+)*inv_squared_c)
            if(tmp<small_e) flag(ix^D,e_)=.true.
         {end do\}
        end if
      end if
    end if

  end subroutine mhd_check_w_semirelati

  subroutine mhd_check_w_origin(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    logical, intent(inout) :: flag(ixI^S,1:nw)

    integer :: ix^D, igrp

    flag=.false.
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if(w(ix^D,rho_)<small_density) flag(ix^D,rho_) = .true.
      if(primitive) then
        if(w(ix^D,p_)<small_pressure) flag(ix^D,e_) = .true.
      else
        if(w(ix^D,e_)-half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)+&
          (^C&w(ix^D,b^C_)**2+))-mhd_uawsom_wave_energy_cell(w(ix^D,:))&
          <small_e) flag(ix^D,e_) = .true.
        if(mhd_uawsom) then
          if(w(ix^D,wAplus_)<zero) flag(ix^D,wAplus_)=.true.
          if(w(ix^D,wAminus_)<zero) flag(ix^D,wAminus_)=.true.
          if(w(ix^D,wkplus_)<zero) flag(ix^D,wkplus_)=.true.
          if(w(ix^D,wkminus_)<zero) flag(ix^D,wkminus_)=.true.
        end if
      end if
      if(mhd_radiation_fld)then
         if(w(ix^D,r_e)<small_r_e) flag(ix^D,r_e) = .true.
      endif
   {end do\}

  end subroutine mhd_check_w_origin

  subroutine mhd_check_w_split(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    logical, intent(inout) :: flag(ixI^S,1:nw)

    double precision :: tmp
    integer :: ix^D

    flag=.false.
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      tmp=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0)
      if(tmp<small_density) flag(ix^D,rho_) = .true.
      if(primitive) then
        if(w(ix^D,p_)+block%equi_vars(ix^D,equi_pe0_,0)<small_pressure) flag(ix^D,e_) = .true.
      else
        tmp=w(ix^D,e_)-half*((^C&w(ix^D,m^C_)**2+)/tmp+(^C&w(ix^D,b^C_)**2+))
        if(tmp+block%equi_vars(ix^D,equi_pe0_,0)*eos%inv_gamma_minus_1<small_e) flag(ix^D,e_) = .true.
      end if
   {end do\}

  end subroutine mhd_check_w_split

  subroutine mhd_check_w_noe(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    logical, intent(inout) :: flag(ixI^S,1:nw)

    integer :: ix^D

    flag=.false.
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if(w(ix^D,rho_)<small_density) flag(ix^D,rho_) = .true.
   {end do\}

  end subroutine mhd_check_w_noe

  subroutine mhd_check_w_inte(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    logical, intent(inout) :: flag(ixI^S,1:nw)

    integer :: ix^D

    flag=.false.
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if(w(ix^D,rho_)<small_density) flag(ix^D,rho_) = .true.
      if(primitive) then
        if(w(ix^D,p_) < small_pressure) flag(ix^D,e_) = .true.
      else
        if(w(ix^D,e_)<small_e) flag(ix^D,e_) = .true.
      end if
   {end do\}

  end subroutine mhd_check_w_inte

  subroutine mhd_check_w_hde(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    logical, intent(inout) :: flag(ixI^S,1:nw)

    integer :: ix^D

    flag=.false.
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if(w(ix^D,rho_)<small_density) flag(ix^D,rho_) = .true.
      if(primitive) then
        if(w(ix^D,p_)<small_pressure) flag(ix^D,e_) = .true.
      else
        if(w(ix^D,e_)-half*(^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)<small_e) flag(ix^D,e_) = .true.
      end if
   {end do\}

  end subroutine mhd_check_w_hde

  subroutine mhd_bound_fip(primitive, ixI^L, ixO^L, w)
    use mod_global_parameters
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rho_safe(ixI^S), fip_prim(ixI^S)

    if (.not. mhd_fip) return

    if (primitive) then
      w(ixO^S,fip_) = min(maxfip, max(minfip, w(ixO^S,fip_)))
    else
      if (has_equi_rho_and_p) then
        rho_safe(ixO^S) = max(w(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,b0i), small_density)
      else
        rho_safe(ixO^S) = max(w(ixO^S,rho_), small_density)
      end if
      fip_prim(ixO^S) = w(ixO^S,fip_) / rho_safe(ixO^S)
      fip_prim(ixO^S) = min(maxfip, max(minfip, fip_prim(ixO^S)))
      w(ixO^S,fip_) = rho_safe(ixO^S) * fip_prim(ixO^S)
    end if
  end subroutine mhd_bound_fip

  !> Return local UAWSoM closure coefficients in code units.  User cases may
  !> replace all three profiles through usr_uawsom_coefficients.
  subroutine mhd_uawsom_get_coefficients(w,x,ixI^L,ixO^L,primitive,zeta,radius,lperp_A)
    use mod_global_parameters
    use mod_usr_methods, only: usr_uawsom_coefficients
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    logical, intent(in) :: primitive
    double precision, intent(out) :: zeta(ixI^S), radius(ixI^S), lperp_A(ixI^S)
    double precision :: Bmag(ixI^S), Btotal(ixI^S,1:ndir), height_base
    integer :: idir

    if(associated(usr_uawsom_coefficients)) then
      call usr_uawsom_coefficients(w,x,ixI^L,ixO^L,primitive,zeta,radius,lperp_A)
    else
      do idir=1,ndir
        if(B0field) then
          Btotal(ixO^S,idir)=w(ixO^S,mag(idir))+block%B0(ixO^S,idir,b0i)
        else
          Btotal(ixO^S,idir)=w(ixO^S,mag(idir))
        end if
      end do
      Bmag(ixO^S)=dsqrt(sum(Btotal(ixO^S,1:ndir)**2,dim=ndim+1))
      select case(mhd_uawsom_height_dim)
      case(1)
        height_base=xprobmin1
      {^NOONED
      case(2)
        height_base=xprobmin2
      }
      {^IFTHREED
      case(3)
        height_base=xprobmin3
      }
      end select
      zeta(ixO^S)=one+(mhd_uawsom_zeta0-one)*dexp(&
           -max(x(ixO^S,mhd_uawsom_height_dim)-height_base,zero)/&
           mhd_uawsom_zeta_scale)
      radius(ixO^S)=mhd_uawsom_thread_radius0*dsqrt(&
           mhd_uawsom_Bref/max(Bmag(ixO^S),smalldouble))
      lperp_A(ixO^S)=mhd_uawsom_alfven_corr_length0*dsqrt(&
           mhd_uawsom_Bref/max(Bmag(ixO^S),smalldouble))
    end if
    if(any(zeta(ixO^S)<=one) .or. any(radius(ixO^S)<=zero) .or. &
       any(lperp_A(ixO^S)<=zero)) &
      call mpistop('mhd_uawsom coefficients require zeta>1 and positive lengths')
  end subroutine mhd_uawsom_get_coefficients

  subroutine mhd_uawsom_rho2_factor(ixI^L,ixO^L,w,x,factor)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: factor(ixI^S)
    double precision :: zeta(ixI^S), radius(ixI^S), lperp_A(ixI^S)
    call mhd_uawsom_get_coefficients(w,x,ixI^L,ixO^L,.false.,zeta,radius,lperp_A)
    factor(ixO^S)=mhd_uawsom_rho2_factor_cell(zeta(ixO^S))
  end subroutine mhd_uawsom_rho2_factor

  elemental pure double precision function mhd_uawsom_rho2_factor_cell(zeta)
    double precision, intent(in) :: zeta
    double precision :: f
    f=mhd_uawsom_filling_factor
    mhd_uawsom_rho2_factor_cell=1.d0+f*(1.d0-f)*(zeta-1.d0)**2/&
         (1.d0+f*zeta-f)**2
  end function mhd_uawsom_rho2_factor_cell

  pure double precision function mhd_uawsom_wave_pressure_cell(wcell,zeta)
    double precision, intent(in) :: wcell(:),zeta
    mhd_uawsom_wave_pressure_cell=0.5d0*(wcell(wAplus_)+wcell(wAminus_))+&
         0.25d0*(zeta+1.d0)*(wcell(wkplus_)+wcell(wkminus_))
  end function mhd_uawsom_wave_pressure_cell

  pure double precision function mhd_uawsom_wave_energy_cell(wcell)
    use mod_constants, only: zero
    double precision, intent(in) :: wcell(:)
    if(mhd_uawsom) then
      mhd_uawsom_wave_energy_cell=wcell(wAplus_)+wcell(wAminus_)+&
           wcell(wkplus_)+wcell(wkminus_)
    else
      mhd_uawsom_wave_energy_cell=zero
    end if
  end function mhd_uawsom_wave_energy_cell

  !> Signed population imbalance used by the Cartesian gradient-vorticity reflection model.
  !> A positive value transfers W_A^+ (or W_k^+) to the minus population;
  !> a negative value transfers the minus population to the plus population.
  !> The factor is zero for a balanced state and remains bounded in [-1,1]
  !> for an arbitrarily imbalanced state.
  elemental pure double precision function mhd_uawsom_imbalance_factor(wplus,wminus)
    double precision, intent(in) :: wplus,wminus

    if(wplus<=0.d0 .or. wminus<=0.d0) then
      mhd_uawsom_imbalance_factor=0.d0
    else if(4.d0*wminus<=wplus) then
      mhd_uawsom_imbalance_factor=1.d0-2.d0*dsqrt(wminus/wplus)
    else if(4.d0*wplus<=wminus) then
      mhd_uawsom_imbalance_factor=-(1.d0-2.d0*dsqrt(wplus/wminus))
    else
      mhd_uawsom_imbalance_factor=0.d0
    end if
  end function mhd_uawsom_imbalance_factor

  !> UAWSoM compression, nonlinear damping, and optional Alfven reflection.
  !> Total energy is intentionally unchanged by damping: because it contains
  !> the four W variables, their loss is recovered as gas internal energy.
  subroutine mhd_add_source_uawsom(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x)
    use mod_global_parameters
    use mod_geometry, only: divvector, gradient, curlvector
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), &
         wCTprim(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: v(ixI^S,1:ndir), divv(ixI^S)
    double precision :: zeta(ixI^S), radius(ixI^S), lperp_A(ixI^S)
    double precision :: lperp_k(ixI^S), rho_e(ixI^S)
    double precision :: gamma_plus(ixI^S), gamma_minus(ixI^S)
    double precision :: gamma_kplus(ixI^S), gamma_kminus(ixI^S)
    double precision :: dampAp(ixI^S), dampAm(ixI^S), dampkp(ixI^S), dampkm(ixI^S)
    double precision :: Btotal(ixI^S,1:ndir), Bunit(ixI^S,1:ndir), Bmag(ixI^S)
    double precision :: va(ixI^S), vk(ixI^S), lnva(ixI^S), lnvk(ixI^S)
    double precision :: grad_component(ixI^S), grad_A(ixI^S), grad_k(ixI^S)
    double precision :: curlv(ixI^S,1:3), vorticity(ixI^S)
    double precision :: rimb_A(ixI^S), rlim_A(ixI^S), rlim_k(ixI^S)
    double precision :: refl_rate(ixI^S), donor_reflection(ixI^S)
    double precision :: transfer_A(ixI^S), transfer_k(ixI^S), donor(ixI^S)
    double precision :: imbalance_A(ixI^S), imbalance_k(ixI^S)
    double precision :: wave_plus(ixI^S), wave_minus(ixI^S)
    double precision :: wave_kplus(ixI^S), wave_kminus(ixI^S)
    double precision :: f, kink_pressure(ixI^S)
    integer :: idir, idirmin
    logical :: have_reflection, have_alfven_reflection, have_kink_reflection

    call mhd_get_v(wCT,x,ixI^L,ixI^L,v)
    call divvector(v,ixI^L,ixO^L,divv)
    ! Reflection gradients need one-cell ghost values, so obtain the closure
    ! profiles on ixI rather than only on the source update region.
    call mhd_uawsom_get_coefficients(wCT,x,ixI^L,ixI^L,.false.,zeta,radius,lperp_A)
    f=mhd_uawsom_filling_factor
    rho_e(ixI^S)=max(wCT(ixI^S,rho_),small_density)/&
         max(one+f*zeta(ixI^S)-f,smalldouble)
    lperp_k(ixI^S)=dsqrt(10.d0)*dsqrt(f*dpi)*radius(ixI^S)*&
         (zeta(ixI^S)+one-f)**1.5d0/&
         ((zeta(ixI^S)-one)*(one-f**2.5d0))
    gamma_plus(ixO^S)=two/lperp_A(ixO^S)*dsqrt(&
         max(wCT(ixO^S,wAminus_),zero)/max(wCT(ixO^S,rho_),small_density))
    gamma_minus(ixO^S)=two/lperp_A(ixO^S)*dsqrt(&
         max(wCT(ixO^S,wAplus_),zero)/max(wCT(ixO^S,rho_),small_density))
    dampAp(ixO^S)=gamma_plus(ixO^S)*max(wCT(ixO^S,wAplus_),zero)
    dampAm(ixO^S)=gamma_minus(ixO^S)*max(wCT(ixO^S,wAminus_),zero)
    dampkp(ixO^S)=max(wCT(ixO^S,wkplus_),zero)**1.5d0/&
         (dsqrt(rho_e(ixO^S))*lperp_k(ixO^S))
    dampkm(ixO^S)=max(wCT(ixO^S,wkminus_),zero)**1.5d0/&
         (dsqrt(rho_e(ixO^S))*lperp_k(ixO^S))
    gamma_kplus(ixO^S)=dsqrt(max(wCT(ixO^S,wkplus_),zero)/&
         max(rho_e(ixO^S),small_density))/max(lperp_k(ixO^S),smalldouble)
    gamma_kminus(ixO^S)=dsqrt(max(wCT(ixO^S,wkminus_),zero)/&
         max(rho_e(ixO^S),small_density))/max(lperp_k(ixO^S),smalldouble)

    w(ixO^S,wAplus_)=w(ixO^S,wAplus_)-qdt*(&
         half*divv(ixO^S)*wCT(ixO^S,wAplus_)+dampAp(ixO^S))
    w(ixO^S,wAminus_)=w(ixO^S,wAminus_)-qdt*(&
         half*divv(ixO^S)*wCT(ixO^S,wAminus_)+dampAm(ixO^S))
    w(ixO^S,wkplus_)=w(ixO^S,wkplus_)-qdt*(&
         half*divv(ixO^S)*wCT(ixO^S,wkplus_)+dampkp(ixO^S))
    w(ixO^S,wkminus_)=w(ixO^S,wkminus_)-qdt*(&
         half*divv(ixO^S)*wCT(ixO^S,wkminus_)+dampkm(ixO^S))

    kink_pressure(ixO^S)=quarter*(zeta(ixO^S)+one)*&
         (wCT(ixO^S,wkplus_)+wCT(ixO^S,wkminus_))
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*(zeta(ixO^S)-one)/&
         (zeta(ixO^S)+one)*kink_pressure(ixO^S)*divv(ixO^S)

    have_alfven_reflection=.false.
    if(mhd_uawsom_reflection) then
      select case(trim(mhd_uawsom_reflection_mode))
      case('one_dimensional_gradient')
        have_alfven_reflection=(mhd_uawsom_sigma>zero)
      case('cartesian_gradient_vorticity')
        have_alfven_reflection=.true.
      end select
    end if
    have_kink_reflection=mhd_uawsom_kink_reflection
    have_reflection=have_alfven_reflection .or. have_kink_reflection
    if(have_reflection) then
      do idir=1,ndir
        if(B0field) then
          Btotal(ixI^S,idir)=wCT(ixI^S,mag(idir))+block%B0(ixI^S,idir,b0i)
        else
          Btotal(ixI^S,idir)=wCT(ixI^S,mag(idir))
        end if
      end do
      Bmag(ixI^S)=dsqrt(sum(Btotal(ixI^S,1:ndir)**2,dim=ndim+1))
      do idir=1,ndir
        Bunit(ixI^S,idir)=Btotal(ixI^S,idir)/max(Bmag(ixI^S),smalldouble)
      end do
      va(ixI^S)=Bmag(ixI^S)/dsqrt(max(wCT(ixI^S,rho_),small_density))
      vk(ixI^S)=Bmag(ixI^S)/dsqrt(max(rho_e(ixI^S)*&
           (zeta(ixI^S)+one)/two,smalldouble))
      lnva(ixI^S)=dlog(max(va(ixI^S),smalldouble))
      lnvk(ixI^S)=dlog(max(vk(ixI^S),smalldouble))

      ! Apply exchange after the ordinary source update.  The donor cap below
      ! uses that actual post-damping state; it does not silently repair a
      ! negative state produced by compression or damping.  The standard
      ! small-value handling remains responsible for such a state.
      if(have_alfven_reflection) then
        wave_plus(ixO^S)=max(wCT(ixO^S,wAplus_),zero)
        wave_minus(ixO^S)=max(wCT(ixO^S,wAminus_),zero)

        select case(trim(mhd_uawsom_reflection_mode))
        case('one_dimensional_gradient')
          ! Retain the original one-dimensional source and sign convention:
          ! a positive transfer removes W_A^+ and adds W_A^-.
          call gradient(va,ixI^L,ixO^L,mhd_uawsom_height_dim,grad_component)
          refl_rate(ixO^S)=mhd_uawsom_sigma*&
               (wCTprim(ixO^S,mom(mhd_uawsom_height_dim))+va(ixO^S))/&
               max(va(ixO^S),smalldouble)*grad_component(ixO^S)
          donor_reflection(ixO^S)=merge(max(wCT(ixO^S,wAplus_),zero),&
               max(wCT(ixO^S,wAminus_),zero),refl_rate(ixO^S)>=zero)
          transfer_A(ixO^S)=refl_rate(ixO^S)*donor_reflection(ixO^S)
        case('cartesian_gradient_vorticity')
          grad_A(ixO^S)=zero
          do idir=1,ndim
            call gradient(lnva,ixI^L,ixO^L,idir,grad_component)
            grad_A(ixO^S)=grad_A(ixO^S)+va(ixO^S)*Bunit(ixO^S,idir)*&
                 grad_component(ixO^S)
          end do
          call curlvector(v,ixI^L,ixO^L,curlv,idirmin,1,ndir)
          vorticity(ixO^S)=zero
          do idir=1,ndir
            vorticity(ixO^S)=vorticity(ixO^S)+Bunit(ixO^S,idir)*&
                 curlv(ixO^S,idir)
          end do
          ! Sigma belongs to the one-dimensional gradient source only; it is
          ! absent from the Cartesian gradient-vorticity R_imb/R_lim closure.
          rimb_A(ixO^S)=dsqrt(grad_A(ixO^S)**2+vorticity(ixO^S)**2)
          rlim_A(ixO^S)=min(rimb_A(ixO^S),&
               max(gamma_plus(ixO^S),gamma_minus(ixO^S)))
          imbalance_A(ixO^S)=mhd_uawsom_imbalance_factor(&
               wave_plus(ixO^S),wave_minus(ixO^S))
          transfer_A(ixO^S)=rlim_A(ixO^S)*imbalance_A(ixO^S)*&
               dsqrt(wave_plus(ixO^S)*wave_minus(ixO^S))
        end select

        if(qdt>zero) then
          donor(ixO^S)=merge(max(w(ixO^S,wAplus_),zero),&
               max(w(ixO^S,wAminus_),zero),transfer_A(ixO^S)>=zero)
          transfer_A(ixO^S)=sign(min(abs(transfer_A(ixO^S)),&
               donor(ixO^S)/qdt),transfer_A(ixO^S))
        else
          transfer_A(ixO^S)=zero
        end if
        w(ixO^S,wAplus_)=w(ixO^S,wAplus_)-qdt*transfer_A(ixO^S)
        w(ixO^S,wAminus_)=w(ixO^S,wAminus_)+qdt*transfer_A(ixO^S)
      end if

      if(have_kink_reflection) then
        wave_kplus(ixO^S)=max(wCT(ixO^S,wkplus_),zero)
        wave_kminus(ixO^S)=max(wCT(ixO^S,wkminus_),zero)
        grad_k(ixO^S)=zero
        do idir=1,ndim
          call gradient(lnvk,ixI^L,ixO^L,idir,grad_component)
          grad_k(ixO^S)=grad_k(ixO^S)+vk(ixO^S)*Bunit(ixO^S,idir)*&
               grad_component(ixO^S)
        end do
        ! Kink reflection deliberately uses only the kink-speed gradient;
        ! the field-aligned velocity-vorticity term belongs to Alfvén waves,
        ! and sigma is not part of this limiter.
        rlim_k(ixO^S)=min(abs(grad_k(ixO^S)),&
             max(gamma_kplus(ixO^S),gamma_kminus(ixO^S)))
        imbalance_k(ixO^S)=mhd_uawsom_imbalance_factor(&
             wave_kplus(ixO^S),wave_kminus(ixO^S))
        transfer_k(ixO^S)=rlim_k(ixO^S)*imbalance_k(ixO^S)*&
             dsqrt(wave_kplus(ixO^S)*wave_kminus(ixO^S))
        if(qdt>zero) then
          donor(ixO^S)=merge(max(w(ixO^S,wkplus_),zero),&
               max(w(ixO^S,wkminus_),zero),transfer_k(ixO^S)>=zero)
          transfer_k(ixO^S)=sign(min(abs(transfer_k(ixO^S)),&
               donor(ixO^S)/qdt),transfer_k(ixO^S))
        else
          transfer_k(ixO^S)=zero
        end if
        w(ixO^S,wkplus_)=w(ixO^S,wkplus_)-qdt*transfer_k(ixO^S)
        w(ixO^S,wkminus_)=w(ixO^S,wkminus_)+qdt*transfer_k(ixO^S)
      end if
    end if
  end subroutine mhd_add_source_uawsom

  !> Transform internal energy to total energy
  subroutine mhd_ei_to_e(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer :: ix^D

    if(has_equi_rho_and_p) then
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
                        +(^C&w(ix^D,b^C_)**2+))&
                  +mhd_uawsom_wave_energy_cell(w(ix^D,:))
     {end do\}
    end if
  end subroutine mhd_ei_to_e

  !> Transform internal energy to hydrodynamic energy
  subroutine mhd_ei_to_e_hde(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer :: ix^D

   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! Calculate e = ei + ek
      w(ix^D,e_)=w(ix^D,e_)&
                +half*(^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)
   {end do\}

  end subroutine mhd_ei_to_e_hde

  !> Transform internal energy to total energy and velocity to momentum
  subroutine mhd_ei_to_e_semirelati(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    w(ixO^S,p_)=w(ixO^S,e_)*eos%gamma_minus_1
    ! call eos%to_conserved(ixI^L,ixO^L,w,x)
    call eos%to_conserved(ixI^L,ixO^L,w,x)

  end subroutine mhd_ei_to_e_semirelati

  !> Transform total energy to internal energy
  subroutine mhd_e_to_ei(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer :: ix^D

    if(has_equi_rho_and_p) then
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
                        +(^C&w(ix^D,b^C_)**2+))&
                  -mhd_uawsom_wave_energy_cell(w(ix^D,:))
     {end do\}
    end if

    if(fix_small_values) then
      call mhd_handle_small_ei(w,x,ixI^L,ixI^L,e_,'mhd_e_to_ei')
    end if

  end subroutine mhd_e_to_ei

  !> Wrapper: e_to_ei + cache log10(nH) in wextra for LTE TC fast path.
  !> During STS substeps density is invariant, so log10(nH) is computed once
  !> per STS cycle (in sts_before_first_cycle hook) and reused across all substeps.
  subroutine mhd_e_to_ei_and_cache_log_nH(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    call mhd_e_to_ei(ixI^L,ixO^L,w,x)
    block%wextra(ixO^S, iw_log_nH) = dlog10(w(ixO^S, rho_) / eos%nH2rhoFactor)
  end subroutine mhd_e_to_ei_and_cache_log_nH

  !> Transform hydrodynamic energy to internal energy
  subroutine mhd_e_to_ei_hde(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer :: ix^D

   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! Calculate ei = e - ek
      w(ix^D,e_)=w(ix^D,e_)&
                -half*(^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)
   {end do\}

    if(fix_small_values) then
      call mhd_handle_small_ei(w,x,ixI^L,ixI^L,e_,'mhd_e_to_ei_hde')
    end if

  end subroutine mhd_e_to_ei_hde

  !> Transform total energy to internal energy and momentum to velocity
  subroutine mhd_e_to_ei_semirelati(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    call eos%to_primitive(ixI^L,ixO^L,w,x)
    w(ixO^S,e_)=w(ixO^S,p_)*eos%inv_gamma_minus_1

  end subroutine mhd_e_to_ei_semirelati

  subroutine mhd_handle_small_values_semirelati(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    double precision :: e(ixI^S,1:ndir), pressure(ixI^S), v(ixI^S,1:ndir)
    double precision :: tmp, factor
    integer :: ix^D
    logical :: flag(ixI^S,1:nw)

    flag=.false.
    where(w(ixO^S,rho_) < small_density) flag(ixO^S,rho_) = .true.

    if(mhd_energy) then
      if(primitive) then
        where(w(ixO^S,p_) < small_pressure) flag(ixO^S,e_) = .true.
      else
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          ! Convert momentum to velocity
          tmp=(^C&w(ix^D,b^C_)*w(ix^D,m^C_)+)*inv_squared_c
          factor=1.0d0/(w(ix^D,rho_)*(w(ix^D,rho_)+(^C&w(ix^D,b^C_)**2+)*inv_squared_c))
          ^C&v(ix^D,^C)=factor*(w(ix^D,m^C_)*w(ix^D,rho_)+w(ix^D,b^C_)*tmp)\
          ! E=Bxv
          {^IFTHREEC
          e(ix^D,1)=w(ix^D,b2_)*v(ix^D,3)-w(ix^D,b3_)*v(ix^D,2)
          e(ix^D,2)=w(ix^D,b3_)*v(ix^D,1)-w(ix^D,b1_)*v(ix^D,3)
          e(ix^D,3)=w(ix^D,b1_)*v(ix^D,2)-w(ix^D,b2_)*v(ix^D,1)
          }
          {^IFTWOC
          e(ix^D,1)=zero
          e(ix^D,2)=w(ix^D,b1_)*v(ix^D,2)-w(ix^D,b2_)*v(ix^D,1)
          }
          {^IFONEC
          e(ix^D,1)=zero
          }
          ! Calculate pressure = (gamma-1) * (e-eK-eB-eE)
          pressure(ix^D)=eos%gamma_minus_1*(w(ix^D,e_)&
                     -half*((^C&v(ix^D,^C)**2+)*w(ix^D,rho_)&
                     +(^C&w(ix^D,b^C_)**2+)+(^C&e(ix^D,^C)**2+)*inv_squared_c))
          if(pressure(ix^D) < small_pressure) flag(ix^D,p_) = .true.
       {end do\}
      end if
    end if

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(flag(ix^D,rho_)) then
            w(ix^D,rho_) = small_density
         ^C&w(ix^D,m^C_)=0.d0\
          end if
          if(mhd_energy) then
            if(primitive) then
              if(flag(ix^D,e_)) w(ix^D,p_) = small_pressure
            else
              if(flag(ix^D,e_)) then
                w(ix^D,e_)=small_pressure*eos%inv_gamma_minus_1+half*((^C&v(ix^D,^C)**2+)*w(ix^D,rho_)&
                           +(^C&w(ix^D,b^C_)**2+)+(^C&e(ix^D,^C)**2+)*inv_squared_c)
              end if
            end if
          end if
       {end do\}
      case ("average")
        ! do averaging of density
        call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
        if(mhd_energy) then
          if(primitive) then
            call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
          else
            w(ixO^S,e_)=pressure(ixO^S)
            call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
               w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1+half*((^C&v(ix^D,^C)**2+)*w(ix^D,rho_)&
                          +(^C&w(ix^D,b^C_)**2+)+(^C&e(ix^D,^C)**2+)*inv_squared_c)
            {end do\}
          end if
        end if
      case default
        if(.not.primitive) then
          ! change to primitive variables
          w(ixO^S,mom(1:ndir))=v(ixO^S,1:ndir)
          if(mhd_energy) w(ixO^S,e_)=pressure(ixO^S)
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
    if (mhd_fip) call mhd_bound_fip(primitive, ixI^L, ixO^L, w)
  end subroutine mhd_handle_small_values_semirelati

  subroutine mhd_handle_small_values_origin(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: ix^D, igrp
    logical :: flag(ixI^S,1:nw)

    call phys_check_w(primitive, ixI^L, ixO^L, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(flag(ix^D,rho_)) w(ix^D,rho_)=small_density
         {
          if(small_values_fix_iw(m^C_)) then
            if(flag({ix^D},rho_)) w({ix^D},m^C_)=0.0d0
          end if
          \}
          if(primitive) then
            if(flag(ix^D,e_)) w(ix^D,p_)=small_pressure
          else
            if(mhd_uawsom) then
              if(flag(ix^D,wAplus_)) w(ix^D,wAplus_)=zero
              if(flag(ix^D,wAminus_)) w(ix^D,wAminus_)=zero
              if(flag(ix^D,wkplus_)) w(ix^D,wkplus_)=zero
              if(flag(ix^D,wkminus_)) w(ix^D,wkminus_)=zero
            end if
            if(flag(ix^D,e_)) &
              w(ix^D,e_)=small_e+half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)+(^C&w(ix^D,b^C_)**2+))&
                   +mhd_uawsom_wave_energy_cell(w(ix^D,:))
          end if
          if(mhd_radiation_fld)then
            if(small_values_fix_iw(r_e)) then
             if(flag(ix^D,r_e)) w(ix^D,r_e)=small_r_e
            endif
          endif
       {end do\}
      case ("average")
        ! do averaging of density
        call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
        if(primitive)then
          call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
        else
          if(mhd_uawsom) then
            call small_values_average(ixI^L,ixO^L,w,x,flag,wAplus_)
            call small_values_average(ixI^L,ixO^L,w,x,flag,wAminus_)
            call small_values_average(ixI^L,ixO^L,w,x,flag,wkplus_)
            call small_values_average(ixI^L,ixO^L,w,x,flag,wkminus_)
          end if
          ! do averaging of internal energy
         {do ix^DB=ixImin^DB,ixImax^DB\}
            w(ix^D,e_)=w(ix^D,e_)&
                -half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)+(^C&w(ix^D,b^C_)**2+))&
                -mhd_uawsom_wave_energy_cell(w(ix^D,:))
         {end do\}
          call small_values_average(ixI^L, ixO^L, w, x, flag, e_)
          ! convert back
         {do ix^DB=ixImin^DB,ixImax^DB\}
            w(ix^D,e_)=w(ix^D,e_)&
                +half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)+(^C&w(ix^D,b^C_)**2+))&
                +mhd_uawsom_wave_energy_cell(w(ix^D,:))
         {end do\}
        end if
        if(mhd_radiation_fld) then
           call small_values_average(ixI^L, ixO^L, w, x, flag, r_e)
        endif
      case default
        if(.not.primitive) then
          !convert w to primitive
         {do ix^DB=ixOmin^DB,ixOmax^DB\}
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)/w(ix^D,rho_)\
            w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)&
                -half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)+(^C&w(ix^D,b^C_)**2+))&
                -mhd_uawsom_wave_energy_cell(w(ix^D,:)))
         {end do\}
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
    if (mhd_fip) call mhd_bound_fip(primitive, ixI^L, ixO^L, w)
  end subroutine mhd_handle_small_values_origin

  subroutine mhd_handle_small_values_split(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    double precision :: rho
    integer :: ix^D, igrp
    logical :: flag(ixI^S,1:nw)

    call phys_check_w(primitive, ixI^L, ixO^L, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          rho=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0)
          if(flag(ix^D,rho_)) w(ix^D,rho_)=small_density-block%equi_vars(ix^D,equi_rho0_,0)
         {
          if(small_values_fix_iw(m^C_)) then
            if(flag({ix^D},rho_)) w({ix^D},m^C_)=0.0d0
          end if
          \}
          if(primitive) then
            if(flag(ix^D,e_)) w(ix^D,p_)=small_pressure-block%equi_vars(ix^D,equi_pe0_,0)
          else
            if(flag(ix^D,e_)) &
              w(ix^D,e_)=small_e+half*((^C&w(ix^D,m^C_)**2+)/rho+(^C&w(ix^D,b^C_)**2+))&
              -block%equi_vars(ix^D,equi_pe0_,0)*eos%inv_gamma_minus_1
          end if
       {end do\}
      case ("average")
        ! do averaging of density
        call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
        if(primitive)then
          call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
        else
          ! do averaging of internal energy
         {do ix^DB=ixImin^DB,ixImax^DB\}
            rho=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0)
            w(ix^D,e_)=w(ix^D,e_)&
                -half*((^C&w(ix^D,m^C_)**2+)/rho+(^C&w(ix^D,b^C_)**2+))
         {end do\}
          call small_values_average(ixI^L, ixO^L, w, x, flag, e_)
          ! convert back
         {do ix^DB=ixImin^DB,ixImax^DB\}
            rho=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0)
            w(ix^D,e_)=w(ix^D,e_)&
                +half*((^C&w(ix^D,m^C_)**2+)/rho+(^C&w(ix^D,b^C_)**2+))
         {end do\}
        end if
      case default
        if(.not.primitive) then
          !convert w to primitive
         {do ix^DB=ixOmin^DB,ixOmax^DB\}
            rho=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0)
         ^C&w(ix^D,m^C_)=w(ix^D,m^C_)/rho\
            w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)&
                -half*((^C&w(ix^D,m^C_)**2+)*rho+(^C&w(ix^D,b^C_)**2+)))
         {end do\}
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
    if (mhd_fip) call mhd_bound_fip(primitive, ixI^L, ixO^L, w)
  end subroutine mhd_handle_small_values_split

  subroutine mhd_handle_small_values_inte(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: ix^D, igrp
    logical :: flag(ixI^S,1:nw)

    call phys_check_w(primitive, ixI^L, ixO^L, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(flag(ix^D,rho_)) then
            w(ix^D,rho_)=small_density
            ^C&w(ix^D,m^C_)=0.d0\
          end if
          if(primitive) then
            if(flag(ix^D,e_)) w(ix^D,p_)=small_pressure
          else
            if(flag(ix^D,e_)) w(ix^D,e_)=small_e
          end if
       {end do\}
      case ("average")
        ! do averaging of density
        call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
        ! do averaging of internal energy
        call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
      case default
        if(.not.primitive) then
          !convert w to primitive
         {do ix^DB=ixOmin^DB,ixOmax^DB\}
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)/w(ix^D,rho_)\
            w(ix^D,p_)=eos%gamma_minus_1*w(ix^D,e_)
         {end do\}
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
    if (mhd_fip) call mhd_bound_fip(primitive, ixI^L, ixO^L, w)
  end subroutine mhd_handle_small_values_inte

  subroutine mhd_handle_small_values_noe(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: ix^D, igrp
    logical :: flag(ixI^S,1:nw)

    call phys_check_w(primitive, ixI^L, ixO^L, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(flag(ix^D,rho_)) w(ix^D,rho_)=small_density
         {
          if(small_values_fix_iw(m^C_)) then
            if(flag({ix^D},rho_)) w({ix^D},m^C_)=0.0d0
          end if
          \}
       {end do\}
      case ("average")
        ! do averaging of density
        call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
      case default
        if(.not.primitive) then
          !convert w to primitive
         {do ix^DB=ixOmin^DB,ixOmax^DB\}
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)/w(ix^D,rho_)\
         {end do\}
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
    if (mhd_fip) call mhd_bound_fip(primitive, ixI^L, ixO^L, w)
  end subroutine mhd_handle_small_values_noe

  subroutine mhd_handle_small_values_hde(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: ix^D, igrp
    logical :: flag(ixI^S,1:nw)

    call phys_check_w(primitive, ixI^L, ixO^L, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(flag(ix^D,rho_)) then
            w(ix^D,rho_)=small_density
         ^C&w(ix^D,m^C_)=0.d0\
          end if
          if(primitive) then
            if(flag(ix^D,e_)) w(ix^D,p_)=small_pressure
          else
            if(flag(ix^D,e_)) w(ix^D,e_)=small_e+half*(^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)
          end if
       {end do\}
      case ("average")
        ! do averaging of density
        call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
        ! do averaging of energy
        call small_values_average(ixI^L, ixO^L, w, x, flag, e_)
      case default
        if(.not.primitive) then
          !convert w to primitive
         {do ix^DB=ixOmin^DB,ixOmax^DB\}
         ^C&w(ix^D,m^C_)=w(ix^D,m^C_)/w(ix^D,rho_)\
            w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)-half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_))
         {end do\}
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
    if (mhd_fip) call mhd_bound_fip(primitive, ixI^L, ixO^L, w)
  end subroutine mhd_handle_small_values_hde

  !> Calculate v vector
  subroutine mhd_get_v(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,ndir)

    double precision :: rho(ixI^S)
    integer :: idir

    call mhd_get_rho(w,x,ixI^L,ixO^L,rho)

    rho(ixO^S)=1.d0/rho(ixO^S)
    ! Convert momentum to velocity
    do idir = 1, ndir
       v(ixO^S, idir) = w(ixO^S, mom(idir))*rho(ixO^S)
    end do

  end subroutine mhd_get_v

  !> Calculate csound**2 within ixO^L
  subroutine mhd_get_csound2(w,x,ixI^L,ixO^L,cs2)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cs2(ixI^S)

    double precision :: rho, inv_rho, ploc
    integer :: ix^D

    {do ix^DB=ixOmin^DB,ixOmax^DB \}
        if(has_equi_rho_and_p) then
          rho=(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0))
          ploc=(w(ix^D,p_)+block%equi_vars(ix^D,equi_pe0_,0))
        else
          rho=w(ix^D,rho_)
          ploc=w(ix^D,p_)
        end if
        inv_rho=1.d0/rho
        ! sound speed**2 
        cs2(ix^D)=eos%gamma*ploc*inv_rho
   {end do\}
  end subroutine mhd_get_csound2

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine mhd_get_cmax_origin(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)

    double precision :: rho, inv_rho, ploc, cfast2, AvMinCs2, b2, kmax
    double precision :: cs2(ixI^S)
    double precision :: uawsom_zeta(ixI^S), uawsom_radius(ixI^S), uawsom_lperp(ixI^S)
    integer :: ix^D

    if(mhd_hall) kmax = dpi/min({dxlevel(^D)},bigdouble)*half
    if(mhd_uawsom) call mhd_uawsom_get_coefficients(w,x,ixI^L,ixO^L,.true.,&
         uawsom_zeta,uawsom_radius,uawsom_lperp)

    ! Sound speed squared via EoS dispatch (LTE+ionE -> Gamma_1 table; FI -> const gamma).
    call eos%get_csound2(w, x, ixI^L, ixO^L, cs2)

    if(B0field) then
     {do ix^DB=ixOmin^DB,ixOmax^DB \}
        if(has_equi_rho_and_p) then
          rho=(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
          ploc=(w(ix^D,p_)+block%equi_vars(ix^D,equi_pe0_,b0i))
        else
          rho=w(ix^D,rho_)
          ploc=w(ix^D,p_)
        end if
        inv_rho=1.d0/rho
        cmax(ix^D)=cs2(ix^D)
        ! store |B|^2 in v
        b2=(^C&(w(ix^D,b^C_)+block%B0(ix^D,^C,b0i))**2+)
        cfast2=b2*inv_rho+cmax(ix^D)
        AvMinCs2=cfast2**2-4.0d0*cmax(ix^D)*(w(ix^D,mag(idim))+block%B0(ix^D,idim,b0i))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        cmax(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
        if(mhd_hall) then
          ! take the Hall velocity into account: most simple estimate, high k limit:
          ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
          cmax(ix^D)=max(cmax(ix^D),mhd_etah*sqrt(b2)*inv_rho*kmax)
        end if
        if(mhd_uawsom) then
          cmax(ix^D)=max(cmax(ix^D),abs(w(ix^D,mag(idim))+&
               block%B0(ix^D,idim,b0i))*dsqrt(inv_rho))
          cmax(ix^D)=max(cmax(ix^D),abs(w(ix^D,mag(idim))+&
               block%B0(ix^D,idim,b0i))/dsqrt(rho*&
               ((uawsom_zeta(ix^D)+one)/(two*(one+mhd_uawsom_filling_factor*&
               uawsom_zeta(ix^D)-mhd_uawsom_filling_factor)))))
        end if
        cmax(ix^D)=abs(w(ix^D,mom(idim)))+cmax(ix^D)
     {end do\}
    else
     {do ix^DB=ixOmin^DB,ixOmax^DB \}
        if(has_equi_rho_and_p) then
          rho=(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
          ploc=(w(ix^D,p_)+block%equi_vars(ix^D,equi_pe0_,b0i))
        else
          rho=w(ix^D,rho_)
          ploc=w(ix^D,p_)
        end if
        inv_rho=1.d0/rho
        cmax(ix^D)=cs2(ix^D)
        ! store |B|^2 in v
        b2=(^C&w(ix^D,b^C_)**2+)
        cfast2=b2*inv_rho+cmax(ix^D)
        AvMinCs2=cfast2**2-4.0d0*cmax(ix^D)*w(ix^D,mag(idim))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        cmax(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
        if(mhd_hall) then
          ! take the Hall velocity into account: most simple estimate, high k limit:
          ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
          cmax(ix^D)=max(cmax(ix^D),mhd_etah*sqrt(b2)*inv_rho*kmax)
        end if
        if(mhd_uawsom) then
          cmax(ix^D)=max(cmax(ix^D),abs(w(ix^D,mag(idim)))*dsqrt(inv_rho))
          cmax(ix^D)=max(cmax(ix^D),abs(w(ix^D,mag(idim)))/dsqrt(rho*&
               ((uawsom_zeta(ix^D)+one)/(two*(one+mhd_uawsom_filling_factor*&
               uawsom_zeta(ix^D)-mhd_uawsom_filling_factor)))))
        end if
        cmax(ix^D)=abs(w(ix^D,mom(idim)))+cmax(ix^D)
     {end do\}
    end if

  end subroutine mhd_get_cmax_origin

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine mhd_get_cmax_origin_noe(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_adiab, usr_set_gamma

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)

    double precision :: rho, inv_rho, cfast2, AvMinCs2, b2, kmax
    double precision :: adiabs(ixI^S), gammas(ixI^S)
    integer :: ix^D

    if(mhd_hall) kmax = dpi/min({dxlevel(^D)},bigdouble)*half

    if(associated(usr_set_adiab)) then
       call usr_set_adiab(w,x,ixI^L,ixO^L,adiabs)
    else
       adiabs=mhd_adiab
    end if
    if(associated(usr_set_gamma)) then
       call usr_set_gamma(w,x,ixI^L,ixO^L,gammas)
    else
       gammas=eos%gamma
    end if
   {do ix^DB=ixOmin^DB,ixOmax^DB \}
       rho=w(ix^D,rho_)
       inv_rho=1.d0/rho
       ! sound speed**2 
       cmax(ix^D)=gammas(ix^D)*adiabs(ix^D)*rho**(gammas(ix^D)-1.d0)
       ! store |B|^2 in v
       b2=(^C&w(ix^D,b^C_)**2+)
       cfast2=b2*inv_rho+cmax(ix^D)
       AvMinCs2=cfast2**2-4.0d0*cmax(ix^D)*w(ix^D,mag(idim))**2*inv_rho
       if(AvMinCs2<zero) AvMinCs2=zero
       cmax(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
       if(mhd_hall) then
         ! take the Hall velocity into account: most simple estimate, high k limit:
         ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
         cmax(ix^D)=max(cmax(ix^D),mhd_etah*sqrt(b2)*inv_rho*kmax)
       end if
       cmax(ix^D)=abs(w(ix^D,mom(idim)))+cmax(ix^D)
   {end do\}

  end subroutine mhd_get_cmax_origin_noe

  !> Calculate cmax_idim for semirelativistic MHD
  subroutine mhd_get_cmax_semirelati(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout):: cmax(ixI^S)

    double precision :: csound, AvMinCs2, idim_Alfven_speed2
    double precision :: inv_rho, Alfven_speed2, gamma2
    integer :: ix^D

   {do ix^DB=ixOmin^DB,ixOmax^DB \}
      inv_rho=1.d0/w(ix^D,rho_)
      Alfven_speed2=(^C&w(ix^D,b^C_)**2+)*inv_rho
      gamma2=1.0d0/(1.d0+Alfven_speed2*inv_squared_c)
      cmax(ix^D)=1.d0-gamma2*w(ix^D,mom(idim))**2*inv_squared_c
      ! squared sound speed
      csound=eos%gamma*w(ix^D,p_)*inv_rho
      idim_Alfven_speed2=w(ix^D,mag(idim))**2*inv_rho
      ! Va_hat^2+a_hat^2 equation (57)
      ! equation (69)
      Alfven_speed2=Alfven_speed2*cmax(ix^D)+csound*(1.d0+idim_Alfven_speed2*inv_squared_c)
      AvMinCs2=(gamma2*Alfven_speed2)**2-4.0d0*gamma2*csound*idim_Alfven_speed2*cmax(ix^D)
      if(AvMinCs2<zero) AvMinCs2=zero
      ! equation (68) fast magnetosonic wave speed
      csound = sqrt(half*(gamma2*Alfven_speed2+sqrt(AvMinCs2)))
      cmax(ix^D)=gamma2*abs(w(ix^D,mom(idim)))+csound
   {end do\}

  end subroutine mhd_get_cmax_semirelati

  !> Calculate cmax_idim for semirelativistic MHD
  subroutine mhd_get_cmax_semirelati_noe(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_adiab, usr_set_gamma

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout):: cmax(ixI^S)

    double precision :: adiabs(ixI^S), gammas(ixI^S)
    double precision :: csound, AvMinCs2, idim_Alfven_speed2
    double precision :: inv_rho, Alfven_speed2, gamma2
    integer :: ix^D

    if(associated(usr_set_adiab)) then
       call usr_set_adiab(w,x,ixI^L,ixO^L,adiabs)
    else
       adiabs=mhd_adiab
    end if
    if(associated(usr_set_gamma)) then
       call usr_set_gamma(w,x,ixI^L,ixO^L,gammas)
    else
       gammas=eos%gamma
    end if

   {do ix^DB=ixOmin^DB,ixOmax^DB \}
      inv_rho=1.d0/w(ix^D,rho_)
      Alfven_speed2=(^C&w(ix^D,b^C_)**2+)*inv_rho
      gamma2=1.0d0/(1.d0+Alfven_speed2*inv_squared_c)
      cmax(ix^D)=1.d0-gamma2*w(ix^D,mom(idim))**2*inv_squared_c
      csound=gammas(ix^D)*adiabs(ix^D)*w(ix^D,rho_)**(gammas(ix^D)-1.d0)
      idim_Alfven_speed2=w(ix^D,mag(idim))**2*inv_rho
      ! Va_hat^2+a_hat^2 equation (57)
      ! equation (69)
      Alfven_speed2=Alfven_speed2*cmax(ix^D)+csound*(1.d0+idim_Alfven_speed2*inv_squared_c)
      AvMinCs2=(gamma2*Alfven_speed2)**2-4.0d0*gamma2*csound*idim_Alfven_speed2*cmax(ix^D)
      if(AvMinCs2<zero) AvMinCs2=zero
      ! equation (68) fast magnetosonic wave speed
      csound = sqrt(half*(gamma2*Alfven_speed2+sqrt(AvMinCs2)))
      cmax(ix^D)=gamma2*abs(w(ix^D,mom(idim)))+csound
   {end do\}

  end subroutine mhd_get_cmax_semirelati_noe

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine mhd_get_tcutoff(ixI^L,ixO^L,w,x,Tco_local,Tmax_local)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_get_heating
    use mod_radiative_cooling, only: findL, calc_l_extended
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    ! in primitive form
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(out) :: Tco_local,Tmax_local

    double precision, parameter :: trac_delta=0.25d0
    double precision :: Te(ixI^S),lts(ixI^S)
    double precision, dimension(1:ndim) :: Bdir, bunitvec
    double precision, dimension(ixI^S,1:ndim) :: gradT
    double precision :: ltrc,ltrp,altr
    integer :: idims,ix^D,jxO^L,hxO^L,ixA^D,ixB^D
    integer :: jxP^L,hxP^L,ixP^L,ixQ^L
    ! TRAC-7 (Johnston et al. 2021, A&A 654 A2) field-aligned local method
    double precision :: Q_heat(ixI^S), ne(ixI^S), nH_arr(ixI^S)
    double precision :: Bvec(1:ndir)
    double precision :: bmin_c, Bmag, gnorm, Bgtd, L_T, dl_eff, vmag, Bdotv, Binp2
    double precision :: v_n, a_coeff, L1, cooling, net_cool
    double precision :: kappa_par, disc, dx_over_delta, kappa_TRAC, kappa_eff, Tcoff_eff

    if (eos%eos_type == 'LTE' .or. eos%eos_type == 'PI') then
      Te(ixI^S) = w(ixI^S, Te_)
    else
      call eos%get_Rfactor(w,x,ixI^L,ixI^L,Te)
      Te(ixI^S)=w(ixI^S,p_)/(Te(ixI^S)*w(ixI^S,rho_))
    end if
    Tco_local=zero
    Tmax_local=maxval(Te(ixO^S))

    {^IFONED
    select case(mhd_trac_type)
    case(0)
      !> test case, fixed cutoff temperature
      block%wextra(ixI^S,Tcoff_)=2.5d5/unit_temperature
    case(1)
      do ix1=ixOmin1,ixOmax1
        lts(ix1)=0.5d0*abs(Te(ix1+1)-Te(ix1-1))/Te(ix1)
        if(lts(ix1)>trac_delta) then
          Tco_local=max(Tco_local,Te(ix1))
        end if
      end do
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
      call mpistop("mhd_trac_type not allowed for 1D simulation")
    end select
    }
    {^NOONED
    select case(mhd_trac_type)
    case(0)
      !> test case, fixed cutoff temperature
      block%wextra(ixI^S,Tcoff_)=2.5d5/unit_temperature
    case(1,4,6)
      ! temperature gradient at cell centers
      do idims=1,ndim
        call gradient(Te,ixI^L,ixO^L,idims,gradT(ixI^S,idims))
      end do
      if(mhd_trac_type .gt. 1) then
        ! B direction at block center
        Bdir=zero
        if(B0field) then
          {do ixA^D=0,1\}
            ixB^D=(ixOmin^D+ixOmax^D-1)/2+ixA^D;
            Bdir(1:ndim)=Bdir(1:ndim)+w(ixB^D,iw_mag(1:ndim))+block%B0(ixB^D,1:ndim,0)
          {end do\}
        else
          {do ixA^D=0,1\}
            ixB^D=(ixOmin^D+ixOmax^D-1)/2+ixA^D;
            Bdir(1:ndim)=Bdir(1:ndim)+w(ixB^D,iw_mag(1:ndim))
          {end do\}
        end if
        {^IFTWOD
        if(Bdir(1)/=0.d0) then
          block%special_values(3)=sign(1.d0,Bdir(1))/dsqrt(1.d0+(Bdir(2)/Bdir(1))**2)
        else
          block%special_values(3)=0.d0
        end if
        if(Bdir(2)/=0.d0) then
          block%special_values(4)=sign(1.d0,Bdir(2))/dsqrt(1.d0+(Bdir(1)/Bdir(2))**2)
        else
          block%special_values(4)=0.d0
        end if
        }
        {^IFTHREED
        if(Bdir(1)/=0.d0) then
          block%special_values(3)=sign(1.d0,Bdir(1))/dsqrt(1.d0+(Bdir(2)/Bdir(1))**2+&
           (Bdir(3)/Bdir(1))**2)
        else
          block%special_values(3)=0.d0
        end if
        if(Bdir(2)/=0.d0) then
          block%special_values(4)=sign(1.d0,Bdir(2))/dsqrt(1.d0+(Bdir(1)/Bdir(2))**2+&
           (Bdir(3)/Bdir(2))**2)
        else
          block%special_values(4)=0.d0
        end if
        if(Bdir(3)/=0.d0) then
          block%special_values(5)=sign(1.d0,Bdir(3))/dsqrt(1.d0+(Bdir(1)/Bdir(3))**2+&
           (Bdir(2)/Bdir(3))**2)
        else
          block%special_values(5)=0.d0
        end if
        }
      end if
      ! b unit vector: magnetic field direction vector
      block%special_values(1)=zero
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(B0field) then
          ^D&bdir(^D)=w({ix^D},iw_mag(^D))+block%B0({ix^D},^D,0)\
        else
          ^D&bdir(^D)=w({ix^D},iw_mag(^D))\
        end if
        {^IFTWOD
        if(Bdir(1)/=0.d0) then
          bunitvec(1)=sign(1.d0,Bdir(1))/dsqrt(1.d0+(Bdir(2)/Bdir(1))**2)
        else
          bunitvec(1)=0.d0
        end if
        if(Bdir(2)/=0.d0) then
          bunitvec(2)=sign(1.d0,Bdir(2))/dsqrt(1.d0+(Bdir(1)/Bdir(2))**2)
        else
          bunitvec(2)=0.d0
        end if
        ! temperature length scale inversed
        lts(ix^D)=min(block%ds(ix^D,1),block%ds(ix^D,2))*&
          abs(^D&gradT({ix^D},^D)*bunitvec(^D)+)/Te(ix^D)
        }
        {^IFTHREED
        if(Bdir(1)/=0.d0) then
          bunitvec(1)=sign(1.d0,Bdir(1))/dsqrt(1.d0+(Bdir(2)/Bdir(1))**2+(Bdir(3)/Bdir(1))**2)
        else
          bunitvec(1)=0.d0
        end if
        if(Bdir(2)/=0.d0) then
          bunitvec(2)=sign(1.d0,Bdir(2))/dsqrt(1.d0+(Bdir(1)/Bdir(2))**2+(Bdir(3)/Bdir(2))**2)
        else
          bunitvec(2)=0.d0
        end if
        if(Bdir(3)/=0.d0) then
          bunitvec(3)=sign(1.d0,Bdir(3))/dsqrt(1.d0+(Bdir(1)/Bdir(3))**2+(Bdir(2)/Bdir(3))**2)
        else
          bunitvec(3)=0.d0
        end if
        ! temperature length scale inversed
        lts(ix^D)=min(block%ds(ix^D,1),block%ds(ix^D,2),block%ds(ix^D,3))*&
          abs(^D&gradT({ix^D},^D)*bunitvec(^D)+)/Te(ix^D)
        }
        if(lts(ix^D)>trac_delta) then
          block%special_values(1)=max(block%special_values(1),Te(ix^D))
        end if
     {end do\}
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
        ! The outermost ghost cell has no room for a fourth-order one-sided
        ! stencil.  A first-order one-sided gradient keeps the two-ghost-layer
        ! LTRAC coverage without reading outside ixI.
        call gradientF(Te,x,ixI^L,hxP^L,idims,gradT(ixI^S,idims),1,.true.)
        call gradientF(Te,x,ixI^L,jxP^L,idims,gradT(ixI^S,idims),1,.false.)
      end do
      ! b unit vector: magnetic field direction vector
      if(B0field) then
       {do ix^DB=ixPmin^DB,ixPmax^DB\}
          ^D&bdir(^D)=w({ix^D},iw_mag(^D))+block%B0({ix^D},^D,0)\
          {^IFTWOD
          if(Bdir(1)/=0.d0) then
            bunitvec(1)=sign(1.d0,Bdir(1))/dsqrt(1.d0+(Bdir(2)/Bdir(1))**2)
          else
            bunitvec(1)=0.d0
          end if
          if(Bdir(2)/=0.d0) then
            bunitvec(2)=sign(1.d0,Bdir(2))/dsqrt(1.d0+(Bdir(1)/Bdir(2))**2)
          else
            bunitvec(2)=0.d0
          end if
          }
          {^IFTHREED
          if(Bdir(1)/=0.d0) then
            bunitvec(1)=sign(1.d0,Bdir(1))/dsqrt(1.d0+(Bdir(2)/Bdir(1))**2+(Bdir(3)/Bdir(1))**2)
          else
            bunitvec(1)=0.d0
          end if
          if(Bdir(2)/=0.d0) then
            bunitvec(2)=sign(1.d0,Bdir(2))/dsqrt(1.d0+(Bdir(1)/Bdir(2))**2+(Bdir(3)/Bdir(2))**2)
          else
            bunitvec(2)=0.d0
          end if
          if(Bdir(3)/=0.d0) then
            bunitvec(3)=sign(1.d0,Bdir(3))/dsqrt(1.d0+(Bdir(1)/Bdir(3))**2+(Bdir(2)/Bdir(3))**2)
          else
            bunitvec(3)=0.d0
          end if
          }
          ! temperature length scale inversed
          lts(ix^D)=abs(^D&gradT({ix^D},^D)*bunitvec(^D)+)/Te(ix^D)
          ! fraction of cells size to temperature length scale
          lts(ix^D)=min(^D&block%ds({ix^D},^D))*lts(ix^D)
          lts(ix^D)=max(one,(exp(lts(ix^D))/ltrc)**ltrp)
       {end do\}
      else
       {do ix^DB=ixPmin^DB,ixPmax^DB\}
         {^IFTWOD
          if(w(ix^D,iw_mag(1))/=0.d0) then
            bunitvec(1)=sign(1.d0,w(ix^D,iw_mag(1)))/dsqrt(1.d0+(w(ix^D,iw_mag(2))/w(ix^D,iw_mag(1)))**2)
          else
            bunitvec(1)=0.d0
          end if
          if(w(ix^D,iw_mag(2))/=0.d0) then
            bunitvec(2)=sign(1.d0,w(ix^D,iw_mag(2)))/dsqrt(1.d0+(w(ix^D,iw_mag(1))/w(ix^D,iw_mag(2)))**2)
          else
            bunitvec(2)=0.d0
          end if
         }
         {^IFTHREED
          if(w(ix^D,iw_mag(1))/=0.d0) then
            bunitvec(1)=sign(1.d0,w(ix^D,iw_mag(1)))/dsqrt(1.d0+(w(ix^D,iw_mag(2))/w(ix^D,iw_mag(1)))**2+&
              (w(ix^D,iw_mag(3))/w(ix^D,iw_mag(1)))**2)
          else
            bunitvec(1)=0.d0
          end if
          if(w(ix^D,iw_mag(2))/=0.d0) then
            bunitvec(2)=sign(1.d0,w(ix^D,iw_mag(2)))/dsqrt(1.d0+(w(ix^D,iw_mag(1))/w(ix^D,iw_mag(2)))**2+&
              (w(ix^D,iw_mag(3))/w(ix^D,iw_mag(2)))**2)
          else
            bunitvec(2)=0.d0
          end if
          if(w(ix^D,iw_mag(3))/=0.d0) then
            bunitvec(3)=sign(1.d0,w(ix^D,iw_mag(3)))/dsqrt(1.d0+(w(ix^D,iw_mag(1))/w(ix^D,iw_mag(3)))**2+&
              (w(ix^D,iw_mag(2))/w(ix^D,iw_mag(3)))**2)
          else
            bunitvec(3)=0.d0
          end if
         }
          ! temperature length scale inversed
          lts(ix^D)=abs(^D&gradT({ix^D},^D)*bunitvec(^D)+)/Te(ix^D)
          ! fraction of cells size to temperature length scale
          lts(ix^D)=min(^D&block%ds({ix^D},^D))*lts(ix^D)
          lts(ix^D)=max(one,(exp(lts(ix^D))/ltrc)**ltrp)
       {end do\}
      end if
  
      ! need one ghost layer for thermal conductivity
      ixP^L=ixO^L^LADD1;
     {do ix^DB=ixPmin^DB,ixPmax^DB\}
        {^IFTWOD
        altr=0.25d0*((lts(ix1-1,ix2)+two*lts(ix^D)+lts(ix1+1,ix2))*bunitvec(1)**2+&
                     (lts(ix1,ix2-1)+two*lts(ix^D)+lts(ix1,ix2+1))*bunitvec(2)**2)
        block%wextra(ix^D,Tcoff_)=Te(ix^D)*altr**0.4d0
        }
        {^IFTHREED
        altr=0.25d0*((lts(ix1-1,ix2,ix3)+two*lts(ix^D)+lts(ix1+1,ix2,ix3))*bunitvec(1)**2+&
                     (lts(ix1,ix2-1,ix3)+two*lts(ix^D)+lts(ix1,ix2+1,ix3))*bunitvec(2)**2+&
                     (lts(ix1,ix2,ix3-1)+two*lts(ix^D)+lts(ix1,ix2,ix3+1))*bunitvec(3)**2)
        block%wextra(ix^D,Tcoff_)=Te(ix^D)*altr**0.4d0
        }
     {end do\}
    case(7)
      ! Johnston et al. 2021 (A&A 654, A2) local field-aligned TRAC for MHD.
      ! Anisotropic conduction broadens the TR along B, so the length scale (Eq.18) and the
      ! enthalpy flux (Eq.16) are projected on the field, regularized by b_min=0.1 G (paper
      ! value) as B->0 or B perpendicular to grad T. Per-cell adaptive parallel conductivity from the
      ! steady-state balance (Eq.11/12) with the Eq.13 selection; general-EoS cooling n_e n_H L.
      ! (HD isotropic version used the grad-T direction; here it is the magnetic field.)
      call usr_get_heating(Q_heat, ixI^L, ixO^L, w, x)
      call eos%get_ne_nH(ixI^L, ixO^L, w, x, ne, nH_arr)
      block%wextra(ixI^S,Tcoff_) = Te(ixI^S)         ! default (incl. ghost layer): no broadening
      do idims=1,ndim
        call gradient(Te,ixI^L,ixO^L,idims,gradT(ixI^S,idims))
      end do
      bmin_c = 0.1d0/unit_magneticfield              ! 0.1 G in code units
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! full field B-magnitude (Eq.18 numerator, regularized), and B.v, v-magnitude over ndir
        Bmag=bmin_c**2; Bdotv=0.d0; vmag=0.d0
        do idims=1,ndir
          Bvec(idims)=w(ix^D,iw_mag(idims))
          if(B0field) Bvec(idims)=Bvec(idims)+block%B0(ix^D,idims,0)
          Bmag =Bmag +Bvec(idims)**2
          Bdotv=Bdotv+Bvec(idims)*w(ix^D,mom(idims))
          vmag =vmag +w(ix^D,mom(idims))**2
        end do
        Bmag=dsqrt(Bmag); vmag=dsqrt(vmag)
        ! in-plane grad-T norm, in-plane field^2, and in-plane B.grad T (grad T has no z-part in 2D)
        gnorm=0.d0; Binp2=0.d0; Bgtd=0.d0
        do idims=1,ndim
          gnorm=gnorm+gradT(ix^D,idims)**2
          Binp2=Binp2+Bvec(idims)**2
          Bgtd =Bgtd +Bvec(idims)*gradT(ix^D,idims)
        end do
        gnorm=dsqrt(gnorm)
        if(gnorm<smalldouble .or. Binp2<smalldouble) then
          block%wextra(ix^D,Tcoff_)=Te(ix^D)          ! uniform T or no in-plane field: no TR
        else
          ! Eq.18 field-aligned length scale (dabs: conduction is bidirectional along B)
          L_T=Te(ix^D)*Bmag/(dabs(Bgtd)+bmin_c*gnorm)
          ! grid spacing along the (in-plane) field direction b = Bvec/Bmag
          dl_eff=0.d0
          do idims=1,ndim
            dl_eff=dl_eff+(Bvec(idims)/block%ds(ix^D,idims))**2
          end do
          dl_eff=Bmag/dsqrt(dl_eff)
          ! Eq.16 field-aligned enthalpy (mass) flux
          v_n=(dabs(Bdotv)+bmin_c*vmag)/Bmag
          a_coeff=2.5d0*w(ix^D,p_)*v_n/Te(ix^D)
          ! optically-thin cooling n_e n_H Lambda(T), guarded to the table range (positive
          ! in-range test so a non-finite Te falls to zero rather than indexing findL with NaN)
          if(Te(ix^D)>rc_fl%tcoolmin .and. Te(ix^D)<rc_fl%tcoolmax) then
            call findL(Te(ix^D),L1,rc_fl); cooling=L1*ne(ix^D)*nH_arr(ix^D)
          else if(Te(ix^D)>=rc_fl%tcoolmax) then
            call calc_l_extended(Te(ix^D),L1,rc_fl); cooling=L1*ne(ix^D)*nH_arr(ix^D)
          else
            cooling=0.d0
          end if
          net_cool=dabs(cooling-Q_heat(ix^D))
          kappa_par=tc_fl%tc_k_para*Te(ix^D)**2.5d0
          disc=a_coeff**2+4.d0*tc_fl%tc_k_para*Te(ix^D)**1.5d0*net_cool
          dx_over_delta=dl_eff/mhd_trac_delta
          ! Eq.13 selection: under-resolved (Eq.11, keep mass flux) vs over-resolved (Eq.12, limiter)
          if(L_T<=2.d0*dx_over_delta) then
            kappa_TRAC=(a_coeff+dsqrt(disc))/(2.d0/dx_over_delta)
          else
            kappa_TRAC=dsqrt(4.d0*tc_fl%tc_k_para*Te(ix^D)**1.5d0*net_cool)/(2.d0/dx_over_delta)
          end if
          kappa_eff=max(kappa_TRAC,kappa_par)
          Tcoff_eff=(kappa_eff/tc_fl%tc_k_para)**0.4d0
          block%wextra(ix^D,Tcoff_)=max(Te(ix^D),Tcoff_eff)
        end if
      {end do\}
    case(3,5)
      !> do nothing here
    case default
      call mpistop("unknown mhd_trac_type")
    end select
    }
  end subroutine mhd_get_tcutoff

  !> get H speed for H-correction to fix the carbuncle problem at grid-aligned shock front
  subroutine mhd_get_H_speed(wprim,x,ixI^L,ixO^L,idim,Hspeed)
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
      if(has_equi_rho_and_p) then
         call mhd_get_csound_prim_split(wprim,x,ixI^L,ixA^L,id,tmp)
      else
         call mhd_get_csound_prim(wprim,x,ixI^L,ixA^L,id,tmp)
      endif
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

  end subroutine mhd_get_H_speed

  !> Estimating bounds for the minimum and maximum signal velocities without split
  subroutine mhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,Hspeed,cmax,cmin)
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
      call mhd_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)
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
      call mhd_get_csound_prim(wmean,x,ixI^L,ixO^L,idim,csoundR)
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
      call mhd_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)
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

  end subroutine mhd_get_cbounds

  !> Estimating bounds for the minimum and maximum signal velocities without split
  subroutine mhd_get_cbounds_semirelati(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,Hspeed,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
    double precision, intent(in)    :: Hspeed(ixI^S,1:number_species)

    double precision, dimension(ixO^S) :: csoundL, csoundR, gamma2L, gamma2R
    integer :: ix^D

    ! Miyoshi 2005 JCP 208, 315 equation (67)
    if(mhd_energy) then
      call mhd_get_csound_semirelati(wLp,x,ixI^L,ixO^L,idim,csoundL,gamma2L)
      call mhd_get_csound_semirelati(wRp,x,ixI^L,ixO^L,idim,csoundR,gamma2R)
    else
      call mhd_get_csound_semirelati_noe(wLp,x,ixI^L,ixO^L,idim,csoundL,gamma2L)
      call mhd_get_csound_semirelati_noe(wRp,x,ixI^L,ixO^L,idim,csoundR,gamma2R)
    end if
    if(present(cmin)) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        csoundL(ix^D)=max(csoundL(ix^D),csoundR(ix^D))
        cmin(ix^D,1)=min(gamma2L(ix^D)*wLp(ix^D,mom(idim)),gamma2R(ix^D)*wRp(ix^D,mom(idim)))-csoundL(ix^D)
        cmax(ix^D,1)=max(gamma2L(ix^D)*wLp(ix^D,mom(idim)),gamma2R(ix^D)*wRp(ix^D,mom(idim)))+csoundL(ix^D)
     {end do\}
    else
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        csoundL(ix^D)=max(csoundL(ix^D),csoundR(ix^D))
        cmax(ix^D,1)=max(gamma2L(ix^D)*wLp(ix^D,mom(idim)),gamma2R(ix^D)*wRp(ix^D,mom(idim)))+csoundL(ix^D)
     {end do\}
    end if

  end subroutine mhd_get_cbounds_semirelati

  !> Estimating bounds for the minimum and maximum signal velocities with rho split
  subroutine mhd_get_cbounds_split_rho(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,Hspeed,cmax,cmin)
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
      call mhd_get_csound_prim_split(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call mhd_get_csound_prim_split(wRp,x,ixI^L,ixO^L,idim,csoundR)
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
      call mhd_get_csound_prim_split(wmean,x,ixI^L,ixO^L,idim,csoundR)
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
      call mhd_get_csound_prim_split(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call mhd_get_csound_prim_split(wRp,x,ixI^L,ixO^L,idim,csoundR)
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

  end subroutine mhd_get_cbounds_split_rho

  !> prepare velocities for ct methods
  subroutine mhd_get_ct_velocity_average(vcts,wLp,wRp,ixI^L,ixO^L,idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: cmax(ixI^S)
    double precision, intent(in), optional :: cmin(ixI^S)
    type(ct_velocity), intent(inout):: vcts

  end subroutine mhd_get_ct_velocity_average

  subroutine mhd_get_ct_velocity_contact(vcts,wLp,wRp,ixI^L,ixO^L,idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: cmax(ixI^S)
    double precision, intent(in), optional :: cmin(ixI^S)
    type(ct_velocity), intent(inout):: vcts

    if(.not.allocated(vcts%vnorm)) allocate(vcts%vnorm(ixI^S,1:ndim))
    ! get average normal velocity at cell faces
    vcts%vnorm(ixO^S,idim)=0.5d0*(wLp(ixO^S,mom(idim))+wRp(ixO^S,mom(idim)))

  end subroutine mhd_get_ct_velocity_contact

  subroutine mhd_get_ct_velocity_hll(vcts,wLp,wRp,ixI^L,ixO^L,idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: cmax(ixI^S)
    double precision, intent(in), optional :: cmin(ixI^S)
    type(ct_velocity), intent(inout):: vcts

    integer                         :: idimE,idimN

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
         +vcts%cbarmin(ixO^S,idim)*vcts%vbarRC(ixO^S,idim,2))&
        /(vcts%cbarmax(ixO^S,idim)+vcts%cbarmin(ixO^S,idim))

  end subroutine mhd_get_ct_velocity_hll

  !> Calculate modified squared sound speed for FLD
  !> NOTE: only for diagnostic purposes, unused subroutine
  subroutine mhd_get_csrad2(w,x,ixI^L,ixO^L,csound)
    use mod_global_parameters
    
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
        
    double precision :: wprim(ixI^S, nw)

    wprim(ixI^S,1:nw)=w(ixI^S,1:nw)
    call eos%to_primitive(ixI^L,ixO^L,wprim,x)
    call mhd_get_csrad2_prim(wprim,x,ixI^L,ixO^L,csound)

  end subroutine mhd_get_csrad2


  !> Calculate modified squared fast wave speed for FLD
  !> NOTE: w is primitive on entry here!
  !> NOTE: used in FLD module as phys_get_csrad2
  subroutine mhd_get_csrad2_prim(w,x,ixI^L,ixO^L,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)

    double precision :: inv_rho, b2
    double precision :: prad_tensor(ixI^S, 1:ndim, 1:ndim)
    double precision :: prad_max(ixI^S)
    integer :: ix^D

    call mhd_get_pradiation_from_prim(w, x, ixI^L, ixO^L, prad_tensor)

    if(B0field) then
     {do ix^DB=ixOmin^DB,ixOmax^DB \}
        inv_rho=1.d0/w(ix^D,rho_)
        prad_max(ix^D) = maxval(prad_tensor(ix^D,:,:))
        b2=(^C&(w(ix^D,b^C_)+block%B0(ix^D,^C,b0i))**2+)
        csound(ix^D)=(eos%gamma*w(ix^D,p_)+b2+prad_max(ix^D))*inv_rho
     {end do\}
    else
     {do ix^DB=ixOmin^DB,ixOmax^DB \}
        inv_rho=1.d0/w(ix^D,rho_)
        prad_max(ix^D) = maxval(prad_tensor(ix^D,:,:))
        b2=(^C&w(ix^D,b^C_)**2+)
        csound(ix^D)=(eos%gamma*w(ix^D,p_)+b2+prad_max(ix^D))*inv_rho
     {end do\}
    end if

   if(minval(csound(ixO^S))<smalldouble)then
     print *,'issue with squared speed and rad pressure'
     print *,minval(csound(ixO^S))
     print *,minval(prad_max(ixO^S))
     call mpistop("negative squared speed in get_csrad2 for dt")
   endif
    
  end subroutine mhd_get_csrad2_prim

  !> Calculate fast magnetosonic wave speed
  subroutine mhd_get_csound_prim(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_adiab, usr_set_gamma

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixO^S)

    double precision :: adiabs(ixI^S), gammas(ixI^S)
    double precision :: inv_rho, cfast2, AvMinCs2, b2, kmax
    double precision :: cs2(ixI^S)
    integer :: ix^D

    if(mhd_hall) kmax = dpi/min({dxlevel(^D)},bigdouble)*half

    if(.not.mhd_energy) then
      if(associated(usr_set_adiab)) then
         call usr_set_adiab(w,x,ixI^L,ixO^L,adiabs)
      else
         adiabs=mhd_adiab
      end if
      if(associated(usr_set_gamma)) then
         call usr_set_gamma(w,x,ixI^L,ixO^L,gammas)
      else
         gammas=eos%gamma
      end if
    end if

    ! Sound speed squared via EoS dispatch (honours LTE+ionE Gamma_1 table).
    if(mhd_energy) then
      call eos%get_csound2(w, x, ixI^L, ixO^L, cs2)
    end if

    ! store |B|^2 in v
    if(B0field) then
     {do ix^DB=ixOmin^DB,ixOmax^DB \}
        inv_rho=1.d0/w(ix^D,rho_)
        if(mhd_energy) then
          csound(ix^D)=cs2(ix^D)
        else
          csound(ix^D)=gammas(ix^D)*adiabs(ix^D)*w(ix^D,rho_)**(gammas(ix^D)-1.d0)
        end if
        b2=(^C&(w(ix^D,b^C_)+block%B0(ix^D,^C,b0i))**2+)
        cfast2=b2*inv_rho+csound(ix^D)
        AvMinCs2=cfast2**2-4.0d0*csound(ix^D)*(w(ix^D,mag(idim))+&
         block%B0(ix^D,idim,b0i))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        csound(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
        if(mhd_hall) then
          csound(ix^D)=max(csound(ix^D),mhd_etah*sqrt(b2)*inv_rho*kmax)
        end if
     {end do\}
    else
     {do ix^DB=ixOmin^DB,ixOmax^DB \}
        inv_rho=1.d0/w(ix^D,rho_)
        if(mhd_energy) then
          csound(ix^D)=cs2(ix^D)
        else
          csound(ix^D)=gammas(ix^D)*adiabs(ix^D)*w(ix^D,rho_)**(gammas(ix^D)-1.d0)
        end if
        b2=(^C&w(ix^D,b^C_)**2+)
        cfast2=b2*inv_rho+csound(ix^D)
        AvMinCs2=cfast2**2-4.0d0*csound(ix^D)*w(ix^D,mag(idim))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        csound(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
        if(mhd_hall) then
          csound(ix^D)=max(csound(ix^D),mhd_etah*sqrt(b2)*inv_rho*kmax)
        end if
     {end do\}
    end if

  end subroutine mhd_get_csound_prim

  !> Calculate fast magnetosonic wave speed when rho and p are split
  !> hence has_equi_rho_and_p=T
  subroutine mhd_get_csound_prim_split(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixO^S)

    double precision :: rho, inv_rho, cfast2, AvMinCs2, b2, kmax
    integer :: ix^D

    if(mhd_hall) kmax = dpi/min({dxlevel(^D)},bigdouble)*half

    ! store |B|^2 in v
    if(B0field) then
     {do ix^DB=ixOmin^DB,ixOmax^DB \}
        rho=(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
        inv_rho=1.d0/rho
        csound(ix^D)=eos%gamma*(w(ix^D,p_)+block%equi_vars(ix^D,equi_pe0_,b0i))*inv_rho
        b2=(^C&(w(ix^D,b^C_)+block%B0(ix^D,^C,b0i))**2+)
        cfast2=b2*inv_rho+csound(ix^D)
        AvMinCs2=cfast2**2-4.0d0*csound(ix^D)*(w(ix^D,mag(idim))+&
         block%B0(ix^D,idim,b0i))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        csound(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
        if(mhd_hall) then
          csound(ix^D)=max(csound(ix^D),mhd_etah*sqrt(b2)*inv_rho*kmax)
        end if
     {end do\}
    else
     {do ix^DB=ixOmin^DB,ixOmax^DB \}
        rho=(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
        inv_rho=1.d0/rho
        csound(ix^D)=eos%gamma*(w(ix^D,p_)+block%equi_vars(ix^D,equi_pe0_,b0i))*inv_rho
        b2=(^C&w(ix^D,b^C_)**2+)
        cfast2=b2*inv_rho+csound(ix^D)
        AvMinCs2=cfast2**2-4.0d0*csound(ix^D)*w(ix^D,mag(idim))**2*inv_rho
        if(AvMinCs2<zero) AvMinCs2=zero
        csound(ix^D)=sqrt(half*(cfast2+sqrt(AvMinCs2)))
        if(mhd_hall) then
          csound(ix^D)=max(csound(ix^D),mhd_etah*sqrt(b2)*inv_rho*kmax)
        end if
     {end do\}
    end if

  end subroutine mhd_get_csound_prim_split

  !> Calculate cmax_idim for semirelativistic MHD
  subroutine mhd_get_csound_semirelati(w,x,ixI^L,ixO^L,idim,csound,gamma2)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! here w is primitive variables
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixO^S), gamma2(ixO^S)

    double precision :: AvMinCs2, inv_rho, Alfven_speed2, idim_Alfven_speed2
    integer :: ix^D

   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      inv_rho = 1.d0/w(ix^D,rho_)
      ! squared sound speed
      csound(ix^D)=eos%gamma*w(ix^D,p_)*inv_rho
      Alfven_speed2=(^C&w(ix^D,b^C_)**2+)*inv_rho
      gamma2(ix^D) = 1.0d0/(1.d0+Alfven_speed2*inv_squared_c)
      AvMinCs2=1.d0-gamma2(ix^D)*w(ix^D,mom(idim))**2*inv_squared_c
      idim_Alfven_speed2=w(ix^D,mag(idim))**2*inv_rho
      ! Va_hat^2+a_hat^2 equation (57)
      ! equation (69)
      Alfven_speed2=Alfven_speed2*AvMinCs2+csound(ix^D)*(1.d0+idim_Alfven_speed2*inv_squared_c)
      AvMinCs2=(gamma2(ix^D)*Alfven_speed2)**2-4.0d0*gamma2(ix^D)*csound(ix^D)*idim_Alfven_speed2*AvMinCs2
      if(AvMinCs2<zero) AvMinCs2=zero
      ! equation (68) fast magnetosonic speed
      csound(ix^D) = sqrt(half*(gamma2(ix^D)*Alfven_speed2+sqrt(AvMinCs2)))
   {end do\}

  end subroutine mhd_get_csound_semirelati

  !> Calculate cmax_idim for semirelativistic MHD
  subroutine mhd_get_csound_semirelati_noe(w,x,ixI^L,ixO^L,idim,csound,gamma2)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_adiab, usr_set_gamma

    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! here w is primitive variables
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixO^S), gamma2(ixO^S)

    double precision :: adiabs(ixI^S), gammas(ixI^S)
    double precision :: AvMinCs2, inv_rho, Alfven_speed2, idim_Alfven_speed2
    integer :: ix^D

    if(associated(usr_set_adiab)) then
       call usr_set_adiab(w,x,ixI^L,ixO^L,adiabs)
    else
       adiabs=mhd_adiab
    end if
    if(associated(usr_set_gamma)) then
       call usr_set_gamma(w,x,ixI^L,ixO^L,gammas)
    else
       gammas=eos%gamma
    end if
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      inv_rho = 1.d0/w(ix^D,rho_)
      ! squared sound speed
      csound(ix^D)=gammas(ix^D)*adiabs(ix^D)*w(ix^D,rho_)**(gammas(ix^D)-1.d0)
      Alfven_speed2=(^C&w(ix^D,b^C_)**2+)*inv_rho
      gamma2(ix^D) = 1.0d0/(1.d0+Alfven_speed2*inv_squared_c)
      AvMinCs2=1.d0-gamma2(ix^D)*w(ix^D,mom(idim))**2*inv_squared_c
      idim_Alfven_speed2=w(ix^D,mag(idim))**2*inv_rho
      ! Va_hat^2+a_hat^2 equation (57)
      ! equation (69)
      Alfven_speed2=Alfven_speed2*AvMinCs2+csound(ix^D)*(1.d0+idim_Alfven_speed2*inv_squared_c)
      AvMinCs2=(gamma2(ix^D)*Alfven_speed2)**2-4.0d0*gamma2(ix^D)*csound(ix^D)*idim_Alfven_speed2*AvMinCs2
      if(AvMinCs2<zero) AvMinCs2=zero
      ! equation (68) fast magnetosonic speed
      csound(ix^D) = sqrt(half*(gamma2(ix^D)*Alfven_speed2+sqrt(AvMinCs2)))
   {end do\}

  end subroutine mhd_get_csound_semirelati_noe

  ! Thermal pressure and temperature subroutines (mhd_get_pthermal_noe,
  ! _inte, _origin, _semirelati, _hde, _LTE, mhd_get_temperature_from_eint,
  ! _from_etot, _from_etot_LTE, _with_equi, _equi, mhd_get_rho_equi,
  ! mhd_get_pe_equi) are defined in mod_mhd_eos.t. eos%get_temperature_from_{etot,eint}
  ! and tc_fl / rc_fl / te_fl_mhd hooks are bound by bind_eos_to_source.


  !> Calculate radiation pressure within ixO^L
  subroutine mhd_get_pradiation_from_prim(w, x, ixI^L, ixO^L, prad)
    use mod_global_parameters
    use mod_fld
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: prad(ixI^S, 1:ndim, 1:ndim)
       
    call fld_get_radpress(w, x, ixI^L, ixO^L, prad, fld_fl)

  end subroutine mhd_get_pradiation_from_prim

  !> Calculates the sum of the gas pressure and the max Prad tensor element
  subroutine mhd_get_pthermal_plus_pradiation(w, x, ixI^L, ixO^L, pth_plus_prad)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, 1:nw)
    double precision, intent(in)  :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: pth_plus_prad(ixI^S)

    double precision              :: wprim(ixI^S, 1:nw)
    double precision              :: prad_tensor(ixI^S, 1:ndim, 1:ndim)
    double precision              :: prad_max(ixI^S)
    integer :: ix^D
        
    wprim(ixI^S,1:nw)=w(ixI^S,1:nw)
    call eos%to_primitive(ixI^L,ixO^L,wprim,x)
    call mhd_get_pradiation_from_prim(wprim, x, ixI^L, ixO^L, prad_tensor)
    {do ix^D = ixOmin^D,ixOmax^D\}
      prad_max(ix^D) = maxval(prad_tensor(ix^D,:,:))
    {enddo\}
    pth_plus_prad(ixO^S) = wprim(ixO^S,p_) + prad_max(ixO^S)
  end subroutine mhd_get_pthermal_plus_pradiation

  !> Calculates radiation temperature 
  subroutine mhd_get_trad(w, x, ixI^L, ixO^L, trad)
    use mod_global_parameters
    use mod_constants

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: trad(ixI^S)

    trad(ixI^S) = (w(ixI^S,r_e)/arad_norm)**(1.d0/4.d0)

  end subroutine mhd_get_trad

  !> Calculate fluxes within ixO^L without any splitting
  subroutine mhd_get_flux(wC,w,x,ixI^L,ixO^L,idim,f)
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
    double precision             :: R(ixI^S), Te(ixI^S), rho_loc(ixI^S)
    double precision             :: Bvec(ixI^S,1:ndir)
    double precision             :: bgradT(ixI^S), gradTperp_mag(ixI^S)
    double precision             :: nperp(ixI^S,1:ndir)
    double precision             :: uawsom_zeta(ixI^S), uawsom_radius(ixI^S)
    double precision             :: uawsom_lperp(ixI^S), uawsom_denom
    double precision             :: uawsom_pwave, uawsom_Bi
    logical                      :: use_perp_flux
    integer                      :: iw, ix^D, idir

    if(mhd_uawsom) call mhd_uawsom_get_coefficients(w,x,ixI^L,ixO^L,.true.,&
         uawsom_zeta,uawsom_radius,uawsom_lperp)

    if(mhd_internal_e) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Get flux of density
        f(ix^D,rho_)=w(ix^D,mom(idim))*w(ix^D,rho_)
        ! f_i[m_k]=v_i*m_k-b_k*b_i
        ^C&f(ix^D,m^C_)=wC(ix^D,mom(idim))*w(ix^D,m^C_)-w(ix^D,mag(idim))*w(ix^D,b^C_)\
        ! normal one includes total pressure
        f(ix^D,mom(idim))=f(ix^D,mom(idim))+w(ix^D,p_)+half*(^C&w(ix^D,b^C_)**2+)
        ! Get flux of internal energy
        f(ix^D,e_)=w(ix^D,mom(idim))*wC(ix^D,e_)
        ! f_i[b_k]=v_i*b_k-v_k*b_i
        ^C&f(ix^D,b^C_)=w(ix^D,mom(idim))*w(ix^D,b^C_)-w(ix^D,mag(idim))*w(ix^D,m^C_)\
     {end do\}
    else
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! Get flux of density
        f(ix^D,rho_)=w(ix^D,mom(idim))*w(ix^D,rho_)
        ! f_i[m_k]=v_i*m_k-b_k*b_i
        ^C&f(ix^D,m^C_)=wC(ix^D,mom(idim))*w(ix^D,m^C_)-w(ix^D,mag(idim))*w(ix^D,b^C_)\
        ptotal=w(ix^D,p_)+half*(^C&w(ix^D,b^C_)**2+)
        if(mhd_uawsom) then
          uawsom_pwave=mhd_uawsom_wave_pressure_cell(w(ix^D,:),uawsom_zeta(ix^D))
          ptotal=ptotal+uawsom_pwave
        end if
        ! normal one includes total pressure
        f(ix^D,mom(idim))=f(ix^D,mom(idim))+ptotal
        ! Get flux of total energy
        ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
        f(ix^D,e_)=w(ix^D,mom(idim))*(wC(ix^D,e_)+ptotal)&
           -w(ix^D,mag(idim))*(^C&w(ix^D,b^C_)*w(ix^D,m^C_)+)
        if(mhd_uawsom) then
          uawsom_denom=dsqrt(w(ix^D,rho_)*(uawsom_zeta(ix^D)+one)/&
               (two*(one+mhd_uawsom_filling_factor*uawsom_zeta(ix^D)-&
               mhd_uawsom_filling_factor)))
          uawsom_Bi=w(ix^D,mag(idim))
          f(ix^D,e_)=f(ix^D,e_)+uawsom_Bi*(&
               (w(ix^D,wAminus_)-w(ix^D,wAplus_))/dsqrt(w(ix^D,rho_))+&
               (w(ix^D,wkminus_)-w(ix^D,wkplus_))/uawsom_denom)
        end if
        ! f_i[b_k]=v_i*b_k-v_k*b_i
        ^C&f(ix^D,b^C_)=w(ix^D,mom(idim))*w(ix^D,b^C_)-w(ix^D,mag(idim))*w(ix^D,m^C_)\
     {end do\}
    end if
    if(mhd_hall) then
      call mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(total_energy) then
          ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
          f(ix^D,e_)=f(ix^D,e_)+vHall(ix^D,idim)*(^C&w(ix^D,b^C_)**2+)&
               -w(ix^D,mag(idim))*(^C&vHall(ix^D,^C)*w(ix^D,b^C_)+)
        end if
        ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
        ^C&f(ix^D,b^C_)=f(ix^D,b^C_)+vHall(ix^D,idim)*w(ix^D,b^C_)-vHall(ix^D,^C)*w(ix^D,mag(idim))\
     {end do\}
    end if

    if(mhd_glm) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,mag(idim))=w(ix^D,psi_)
        !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
        f(ix^D,psi_) = cmax_global**2*w(ix^D,mag(idim))
     {end do\}
    end if

    if(mhd_radiation_fld) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,r_e)=w(ix^D,mom(idim))*wC(ix^D,r_e)
      {end do\}
    endif

    if (mhd_fip) then
      f(ixO^S,fip_) = w(ixO^S,mom(idim)) * wC(ixO^S,fip_)
    end if
    if(mhd_uawsom) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
      uawsom_denom=dsqrt(w(ix^D,rho_)*(uawsom_zeta(ix^D)+one)/&
           (two*(one+mhd_uawsom_filling_factor*uawsom_zeta(ix^D)-&
           mhd_uawsom_filling_factor)))
      uawsom_Bi=w(ix^D,mag(idim))
      f(ix^D,wAplus_)=w(ix^D,wAplus_)*(w(ix^D,mom(idim))-&
           uawsom_Bi/dsqrt(w(ix^D,rho_)))
      f(ix^D,wAminus_)=w(ix^D,wAminus_)*(w(ix^D,mom(idim))+&
           uawsom_Bi/dsqrt(w(ix^D,rho_)))
      f(ix^D,wkplus_)=w(ix^D,wkplus_)*(w(ix^D,mom(idim))-uawsom_Bi/uawsom_denom)
      f(ix^D,wkminus_)=w(ix^D,wkminus_)*(w(ix^D,mom(idim))+uawsom_Bi/uawsom_denom)
     {end do\}
    end if
    ! Get flux of tracer
    do iw=1,mhd_n_tracer
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,tracer(iw))=w(ix^D,mom(idim))*w(ix^D,tracer(iw))
     {end do\}
    end do

    use_perp_flux = mhd_hyperbolic_tc .and. mhd_hyperbolic_tc_use_perp .and. &
                    trim(mhd_hyperbolic_tc_perp_mode)/='off'
    if(use_perp_flux) then
      call mhd_get_rho(w,x,ixI^L,ixI^L,rho_loc)
      call eos%get_Rfactor(w,x,ixI^L,ixI^L,R)
      if(has_equi_rho_and_p) then
        Te(ixI^S)=(w(ixI^S,p_)+block%equi_vars(ixI^S,equi_pe0_,b0i)) / &
                  (R(ixI^S)*rho_loc(ixI^S))
      else
        Te(ixI^S)=w(ixI^S,p_)/(R(ixI^S)*rho_loc(ixI^S))
      end if
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        do idir=1,ndir
          Bvec(ix^D,idir)=w(ix^D,mag(idir))
        end do
      {end do\}
      call mhd_get_hyperbolic_tc_geometry(ixI^L,ixO^L,Te,Bvec,bgradT,gradTperp_mag,nperp)
    end if

    if(mhd_hyperbolic_tc) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,e_)=f(ix^D,e_)+w(ix^D,qpar_)*w(ix^D,mag(idim))/(dsqrt(^C&w(ix^D,b^C_)**2+)+smalldouble)
        f(ix^D,qpar_)=zero
        if(use_perp_flux) then
          f(ix^D,e_)=f(ix^D,e_)+w(ix^D,qperp_)*nperp(ix^D,idim)
          f(ix^D,qperp_)=zero
        end if
     {end do\}
    end if
  end subroutine mhd_get_flux

  !> Calculate fluxes within ixO^L for case without energy equation, hence without splitting
  !> and assuming polytropic closure
  subroutine mhd_get_flux_noe(wC,w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_set_adiab, usr_set_gamma

    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)

    double precision             :: vHall(ixI^S,1:ndir)
    double precision             :: adiabs(ixI^S), gammas(ixI^S)
    integer                      :: iw, ix^D

    if(associated(usr_set_adiab)) then
       call usr_set_adiab(w,x,ixI^L,ixO^L,adiabs)
    else
       adiabs=mhd_adiab
    end if
    if(associated(usr_set_gamma)) then
       call usr_set_gamma(w,x,ixI^L,ixO^L,gammas)
    else
       gammas=eos%gamma
    end if
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! Get flux of density
      f(ix^D,rho_)=w(ix^D,mom(idim))*w(ix^D,rho_)
      ! f_i[m_k]=v_i*m_k-b_k*b_i
      ^C&f(ix^D,m^C_)=wC(ix^D,mom(idim))*w(ix^D,m^C_)-w(ix^D,mag(idim))*w(ix^D,b^C_)\
      ! normal one includes total pressure
      f(ix^D,mom(idim))=f(ix^D,mom(idim))+adiabs(ix^D)*w(ix^D,rho_)**gammas(ix^D)+half*(^C&w(ix^D,b^C_)**2+)
      ! f_i[b_k]=v_i*b_k-v_k*b_i
      ^C&f(ix^D,b^C_)=w(ix^D,mom(idim))*w(ix^D,b^C_)-w(ix^D,mag(idim))*w(ix^D,m^C_)\
   {end do\}
    if(mhd_hall) then
      call mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
        ^C&f(ix^D,b^C_)=f(ix^D,b^C_)+vHall(ix^D,idim)*w(ix^D,b^C_)-vHall(ix^D,^C)*w(ix^D,mag(idim))\
     {end do\}
    end if
    if(mhd_glm) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,mag(idim))=w(ix^D,psi_)
        !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
        f(ix^D,psi_) = cmax_global**2*w(ix^D,mag(idim))
     {end do\}
    end if
    if (mhd_fip) then
      f(ixO^S,fip_) = w(ixO^S,mom(idim)) * wC(ixO^S,fip_)
    end if
    ! Get flux of tracer
    do iw=1,mhd_n_tracer
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,tracer(iw))=w(ix^D,mom(idim))*w(ix^D,tracer(iw))
     {end do\}
    end do
  end subroutine mhd_get_flux_noe

  !> Calculate fluxes with hydrodynamic energy equation
  subroutine mhd_get_flux_hde(wC,w,x,ixI^L,ixO^L,idim,f)
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
    double precision             :: R(ixI^S), Te(ixI^S), rho_loc(ixI^S)
    double precision             :: Bvec(ixI^S,1:ndir)
    double precision             :: bgradT(ixI^S), gradTperp_mag(ixI^S)
    double precision             :: nperp(ixI^S,1:ndir)
    logical                      :: use_perp_flux
    integer                      :: iw, ix^D, idir

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! Get flux of density
      f(ix^D,rho_)=w(ix^D,mom(idim))*w(ix^D,rho_)
      ! f_i[m_k]=v_i*m_k-b_k*b_i
      ^C&f(ix^D,m^C_)=wC(ix^D,mom(idim))*w(ix^D,m^C_)-w(ix^D,mag(idim))*w(ix^D,b^C_)\
      ! normal one includes total pressure
      f(ix^D,mom(idim))=f(ix^D,mom(idim))+w(ix^D,p_)+half*(^C&w(ix^D,b^C_)**2+)
      ! Get flux of energy
      f(ix^D,e_)=w(ix^D,mom(idim))*(wC(ix^D,e_)+w(ix^D,p_))
      ! f_i[b_k]=v_i*b_k-v_k*b_i
      ^C&f(ix^D,b^C_)=w(ix^D,mom(idim))*w(ix^D,b^C_)-w(ix^D,mag(idim))*w(ix^D,m^C_)\
   {end do\}
    if(mhd_hall) then
      call mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
        ^C&f(ix^D,b^C_)=f(ix^D,b^C_)+vHall(ix^D,idim)*w(ix^D,b^C_)-vHall(ix^D,^C)*w(ix^D,mag(idim))\
     {end do\}
    end if
    if(mhd_glm) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,mag(idim))=w(ix^D,psi_)
        !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
        f(ix^D,psi_) = cmax_global**2*w(ix^D,mag(idim))
     {end do\}
    end if
    if (mhd_fip) then
      f(ixO^S,fip_) = w(ixO^S,mom(idim)) * wC(ixO^S,fip_)
    end if
    ! Get flux of tracer
    do iw=1,mhd_n_tracer
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,tracer(iw))=w(ix^D,mom(idim))*w(ix^D,tracer(iw))
     {end do\}
    end do
    use_perp_flux = mhd_hyperbolic_tc .and. mhd_hyperbolic_tc_use_perp .and. &
                    trim(mhd_hyperbolic_tc_perp_mode)/='off'
    if(use_perp_flux) then
      call mhd_get_rho(w,x,ixI^L,ixI^L,rho_loc)
      call eos%get_Rfactor(w,x,ixI^L,ixI^L,R)
      Te(ixI^S)=w(ixI^S,p_)/(R(ixI^S)*rho_loc(ixI^S))
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        do idir=1,ndir
          Bvec(ix^D,idir)=w(ix^D,mag(idir))
        end do
      {end do\}
      call mhd_get_hyperbolic_tc_geometry(ixI^L,ixO^L,Te,Bvec,bgradT,gradTperp_mag,nperp)
    end if
    if(mhd_hyperbolic_tc) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,e_)=f(ix^D,e_)+w(ix^D,qpar_)*w(ix^D,mag(idim))/(dsqrt(^C&w(ix^D,b^C_)**2+)+smalldouble)
        f(ix^D,qpar_)=zero
        if(use_perp_flux) then
          f(ix^D,e_)=f(ix^D,e_)+w(ix^D,qperp_)*nperp(ix^D,idim)
          f(ix^D,qperp_)=zero
        end if
     {end do\}
    end if
  end subroutine mhd_get_flux_hde

  !> Calculate fluxes within ixO^L with possible splitting
  !> this covers four cases: B0field=T and mhd_internal_e=T (where has_equi_rho_and_p=F)
  !>                         B0field=T and has_equi_rho_and_p=F for total_energy=T
  !>                         B0field=F and has_equi_rho_and_p=T for total_energy=T
  !>                         B0field=T and has_equi_rho_and_p=T for total_energy=T
  subroutine mhd_get_flux_split(wC,w,x,ixI^L,ixO^L,idim,f)
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
    double precision             :: R(ixI^S), Te(ixI^S), rho_loc(ixI^S)
    double precision             :: Bvec(ixI^S,1:ndir)
    double precision             :: bgradT(ixI^S), gradTperp_mag(ixI^S)
    double precision             :: nperp(ixI^S,1:ndir)
    double precision             :: uawsom_zeta(ixI^S), uawsom_radius(ixI^S)
    double precision             :: uawsom_lperp(ixI^S), uawsom_denom
    double precision             :: uawsom_pwave, uawsom_Bi
    logical                      :: use_perp_flux
    integer                      :: iw, ix^D, idir

    if(mhd_uawsom) call mhd_uawsom_get_coefficients(w,x,ixI^L,ixO^L,.true.,&
         uawsom_zeta,uawsom_radius,uawsom_lperp)

   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! Get flux of density
      if(has_equi_rho_and_p) then
        f(ix^D,rho_)=w(ix^D,mom(idim))*(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
      else
        f(ix^D,rho_)=w(ix^D,mom(idim))*w(ix^D,rho_)
      end if

      ptotal=w(ix^D,p_)+half*(^C&w(ix^D,b^C_)**2+)

      if(B0field) then
        ^C&btotal(ix^D,^C)=w(ix^D,b^C_)+block%B0(ix^D,^C,idim)\
        ptotal=ptotal+(^C&w(ix^D,b^C_)*block%B0(ix^D,^C,idim)+)
        ! Get flux of momentum and magnetic field
        ! f_i[m_k]=v_i*m_k-b_k*b_i
        ^C&f(ix^D,m^C_)=wC(ix^D,mom(idim))*w(ix^D,m^C_)-&
                          btotal(ix^D,idim)*w(ix^D,b^C_)-w(ix^D,mag(idim))*block%B0(ix^D,^C,idim)\
        f(ix^D,mom(idim))=f(ix^D,mom(idim))+ptotal
      else
        ^C&btotal(ix^D,^C)=w(ix^D,b^C_)\
        ! Get flux of momentum and magnetic field
        ! f_i[m_k]=v_i*m_k-b_k*b_i
        ^C&f(ix^D,m^C_)=wC(ix^D,mom(idim))*w(ix^D,m^C_)-w(ix^D,mag(idim))*w(ix^D,b^C_)\
        f(ix^D,mom(idim))=f(ix^D,mom(idim))+ptotal
      end if
      if(mhd_uawsom) then
        uawsom_pwave=mhd_uawsom_wave_pressure_cell(w(ix^D,:),uawsom_zeta(ix^D))
        f(ix^D,mom(idim))=f(ix^D,mom(idim))+uawsom_pwave
      end if
      ! f_i[b_k]=v_i*b_k-v_k*b_i
      ^C&f(ix^D,b^C_)=w(ix^D,mom(idim))*btotal(ix^D,^C)-btotal(ix^D,idim)*w(ix^D,m^C_)\

      ! Get flux of energy
      ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
      if(mhd_internal_e) then
        f(ix^D,e_)=w(ix^D,mom(idim))*wC(ix^D,e_)
      else
        f(ix^D,e_)=w(ix^D,mom(idim))*(wC(ix^D,e_)+ptotal)&
           -btotal(ix^D,idim)*(^C&w(ix^D,b^C_)*w(ix^D,m^C_)+)
        if(mhd_uawsom) then
          uawsom_denom=dsqrt(w(ix^D,rho_)*(uawsom_zeta(ix^D)+one)/&
               (two*(one+mhd_uawsom_filling_factor*uawsom_zeta(ix^D)-&
               mhd_uawsom_filling_factor)))
          uawsom_Bi=btotal(ix^D,idim)
          f(ix^D,e_)=f(ix^D,e_)+w(ix^D,mom(idim))*uawsom_pwave+uawsom_Bi*(&
               (w(ix^D,wAminus_)-w(ix^D,wAplus_))/dsqrt(w(ix^D,rho_))+&
               (w(ix^D,wkminus_)-w(ix^D,wkplus_))/uawsom_denom)
        end if
      end if
   {end do\}

    if(mhd_glm) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,mag(idim))=w(ix^D,psi_)
        !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
        f(ix^D,psi_) = cmax_global**2*w(ix^D,mag(idim))
     {end do\}
    end if

    if(mhd_radiation_fld) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,r_e)=w(ix^D,mom(idim))*wC(ix^D,r_e)
      {end do\}
    endif

    if(mhd_hall) then
      call mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
        ^C&f(ix^D,b^C_)=f(ix^D,b^C_)+vHall(ix^D,idim)*btotal(ix^D,^C)-btotal(ix^D,idim)*vHall(ix^D,^C)\
        if(total_energy) then
          ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
          f(ix^D,e_)=f(ix^D,e_)+vHall(ix^D,idim)*(^C&w(ix^D,b^C_)*btotal(ix^D,^C)+)&
                    -btotal(ix^D,idim)*(^C&vHall(ix^D,^C)*w(ix^D,b^C_)+)
        end if
     {end do\}
    end if
    if (mhd_fip) then
      f(ixO^S,fip_) = w(ixO^S,mom(idim)) * wC(ixO^S,fip_)
    end if
    if(mhd_uawsom) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
      uawsom_denom=dsqrt(w(ix^D,rho_)*(uawsom_zeta(ix^D)+one)/&
           (two*(one+mhd_uawsom_filling_factor*uawsom_zeta(ix^D)-&
           mhd_uawsom_filling_factor)))
      uawsom_Bi=btotal(ix^D,idim)
      f(ix^D,wAplus_)=w(ix^D,wAplus_)*(w(ix^D,mom(idim))-&
           uawsom_Bi/dsqrt(w(ix^D,rho_)))
      f(ix^D,wAminus_)=w(ix^D,wAminus_)*(w(ix^D,mom(idim))+&
           uawsom_Bi/dsqrt(w(ix^D,rho_)))
      f(ix^D,wkplus_)=w(ix^D,wkplus_)*(w(ix^D,mom(idim))-uawsom_Bi/uawsom_denom)
      f(ix^D,wkminus_)=w(ix^D,wkminus_)*(w(ix^D,mom(idim))+uawsom_Bi/uawsom_denom)
     {end do\}
    end if
    ! Get flux of tracer
    do iw=1,mhd_n_tracer
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,tracer(iw))=w(ix^D,mom(idim))*w(ix^D,tracer(iw))
     {end do\}
    end do
    use_perp_flux = mhd_hyperbolic_tc .and. mhd_hyperbolic_tc_use_perp .and. &
                    trim(mhd_hyperbolic_tc_perp_mode)/='off'
    if(use_perp_flux) then
      call mhd_get_rho(w,x,ixI^L,ixI^L,rho_loc)
      call eos%get_Rfactor(w,x,ixI^L,ixI^L,R)
      if(has_equi_rho_and_p) then
        Te(ixI^S)=(w(ixI^S,p_)+block%equi_vars(ixI^S,equi_pe0_,b0i)) / &
                  (R(ixI^S)*rho_loc(ixI^S))
      else
        Te(ixI^S)=w(ixI^S,p_)/(R(ixI^S)*rho_loc(ixI^S))
      end if
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        do idir=1,ndir
          Bvec(ix^D,idir)=Btotal(ix^D,idir)
        end do
      {end do\}
      call mhd_get_hyperbolic_tc_geometry(ixI^L,ixO^L,Te,Bvec,bgradT,gradTperp_mag,nperp)
    end if
    if(mhd_hyperbolic_tc) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,e_)=f(ix^D,e_)+w(ix^D,qpar_)*Btotal(ix^D,idim)/(dsqrt(^C&btotal(ix^D,^C)**2+)+smalldouble)
        f(ix^D,qpar_)=zero
        if(use_perp_flux) then
          f(ix^D,e_)=f(ix^D,e_)+w(ix^D,qperp_)*nperp(ix^D,idim)
          f(ix^D,qperp_)=zero
        end if
     {end do\}
    end if
  end subroutine mhd_get_flux_split

  !> Calculate semirelativistic fluxes within ixO^L without any splitting
  subroutine mhd_get_flux_semirelati(wC,w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)
    double precision             :: SA(ixO^S,1:ndir),E(ixO^S,1:ndir),e2
    double precision             :: R(ixI^S), Te(ixI^S), rho_loc(ixI^S)
    double precision             :: Bvec(ixI^S,1:ndir)
    double precision             :: bgradT(ixI^S), gradTperp_mag(ixI^S)
    double precision             :: nperp(ixI^S,1:ndir)
    logical                      :: use_perp_flux
    integer                      :: iw, ix^D, idir

   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! Get flux of density
      f(ix^D,rho_)=w(ix^D,mom(idim))*w(ix^D,rho_)
      ! E=Bxv
      {^IFTHREEC
      E(ix^D,1)=w(ix^D,b2_)*w(ix^D,m3_)-w(ix^D,b3_)*w(ix^D,m2_)
      E(ix^D,2)=w(ix^D,b3_)*w(ix^D,m1_)-w(ix^D,b1_)*w(ix^D,m3_)
      E(ix^D,3)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
      }
      {^IFTWOC
      E(ix^D,1)=zero
      ! switch 2 and 3 to add 3 when ^C is from 1 to 2
      E(ix^D,2)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
      }
      {^IFONEC
      E(ix^D,1)=zero
      }
      e2=(^C&e(ix^D,^C)**2+)
      if(mhd_internal_e) then
        ! Get flux of internal energy
        f(ix^D,e_)=w(ix^D,mom(idim))*wC(ix^D,e_)
      else
        ! S=ExB
        {^IFTHREEC
        SA(ix^D,1)=E(ix^D,2)*w(ix^D,b3_)-E(ix^D,3)*w(ix^D,b2_)
        SA(ix^D,2)=E(ix^D,3)*w(ix^D,b1_)-E(ix^D,1)*w(ix^D,b3_)
        SA(ix^D,3)=E(ix^D,1)*w(ix^D,b2_)-E(ix^D,2)*w(ix^D,b1_)
        }
        {^IFTWOC
        SA(ix^D,1)=-E(ix^D,2)*w(ix^D,b2_)
        SA(ix^D,2)=E(ix^D,2)*w(ix^D,b1_)
        ! set E2 back to 0, after e^2 is stored
        E(ix^D,2)=zero
        }
        {^IFONEC
        SA(ix^D,1)=zero
        }
        ! Get flux of total energy
        f(ix^D,e_)=w(ix^D,mom(idim))*(half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+)+&
                    eos%gamma*w(ix^D,p_)*eos%inv_gamma_minus_1)+SA(ix^D,idim)
      end if
      ! Get flux of momentum
      ^C&f(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,mom(idim))*w(ix^D,m^C_)&
       -w(ix^D,mag(idim))*w(ix^D,b^C_)-E(ix^D,idim)*E(ix^D,^C)*inv_squared_c\
      ! gas pressure + magnetic pressure + electric pressure
      f(ix^D,mom(idim))=f(ix^D,mom(idim))+w(ix^D,p_)+half*((^C&w(ix^D,b^C_)**2+)+e2*inv_squared_c)
      ! compute flux of magnetic field
      ! f_i[b_k]=v_i*b_k-v_k*b_i
      ^C&f(ix^D,b^C_)=w(ix^D,mom(idim))*w(ix^D,b^C_)-w(ix^D,mag(idim))*w(ix^D,m^C_)\
   {end do\}

    if(mhd_glm) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,mag(idim))=w(ix^D,psi_)
        !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
        f(ix^D,psi_)=cmax_global**2*w(ix^D,mag(idim))
     {end do\}
    end if
    if (mhd_fip) then
      f(ixO^S,fip_) = w(ixO^S,mom(idim)) * wC(ixO^S,fip_)
    end if
    ! Get flux of tracer
    do iw=1,mhd_n_tracer
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,tracer(iw))=w(ix^D,mom(idim))*w(ix^D,tracer(iw))
     {end do\}
    end do
    use_perp_flux = mhd_hyperbolic_tc .and. mhd_hyperbolic_tc_use_perp .and. &
                    trim(mhd_hyperbolic_tc_perp_mode)/='off'
    if(use_perp_flux) then
      call mhd_get_rho(w,x,ixI^L,ixI^L,rho_loc)
      call eos%get_Rfactor(w,x,ixI^L,ixI^L,R)
      Te(ixI^S)=w(ixI^S,p_)/(R(ixI^S)*rho_loc(ixI^S))
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        do idir=1,ndir
          Bvec(ix^D,idir)=w(ix^D,mag(idir))
        end do
      {end do\}
      call mhd_get_hyperbolic_tc_geometry(ixI^L,ixO^L,Te,Bvec,bgradT,gradTperp_mag,nperp)
    end if
    if(mhd_hyperbolic_tc) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,e_)=f(ix^D,e_)+w(ix^D,qpar_)*w(ix^D,mag(idim))/(dsqrt(^C&w(ix^D,b^C_)**2+)+smalldouble)
        f(ix^D,qpar_)=zero
        if(use_perp_flux) then
          f(ix^D,e_)=f(ix^D,e_)+w(ix^D,qperp_)*nperp(ix^D,idim)
          f(ix^D,qperp_)=zero
        end if
     {end do\}
    end if
  end subroutine mhd_get_flux_semirelati

  subroutine mhd_get_flux_semirelati_noe(wC,w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_set_adiab, usr_set_gamma

    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)

    double precision             :: adiabs(ixI^S), gammas(ixI^S)
    double precision             :: E(ixO^S,1:ndir),e2
    integer                      :: iw, ix^D

    if(associated(usr_set_adiab)) then
       call usr_set_adiab(w,x,ixI^L,ixO^L,adiabs)
    else
       adiabs=mhd_adiab
    end if
    if(associated(usr_set_gamma)) then
       call usr_set_gamma(w,x,ixI^L,ixO^L,gammas)
    else
       gammas=eos%gamma
    end if
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! Get flux of density
      f(ix^D,rho_)=w(ix^D,mom(idim))*w(ix^D,rho_)
      ! E=Bxv
      {^IFTHREEC
      E(ix^D,1)=w(ix^D,b2_)*w(ix^D,m3_)-w(ix^D,b3_)*w(ix^D,m2_)
      E(ix^D,2)=w(ix^D,b3_)*w(ix^D,m1_)-w(ix^D,b1_)*w(ix^D,m3_)
      E(ix^D,3)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
      e2=(^C&e(ix^D,^C)**2+)
      }
      {^IFTWOC
      E(ix^D,1)=zero
      ! switch 2 and 3 to add 3 when ^C is from 1 to 2
      E(ix^D,2)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
      e2=E(ix^D,2)**2
      E(ix^D,2)=zero
      }
      {^IFONEC
      E(ix^D,1)=zero
      e2=zero
      }
      ! Get flux of momentum
      ^C&f(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,mom(idim))*w(ix^D,m^C_)&
       -w(ix^D,mag(idim))*w(ix^D,b^C_)-E(ix^D,idim)*E(ix^D,^C)*inv_squared_c\
      ! gas pressure + magnetic pressure + electric pressure
      f(ix^D,mom(idim))=f(ix^D,mom(idim))+adiabs(ix^D)*w(ix^D,rho_)**gammas(ix^D)+half*((^C&w(ix^D,b^C_)**2+)+e2*inv_squared_c)
      ! compute flux of magnetic field
      ! f_i[b_k]=v_i*b_k-v_k*b_i
      ^C&f(ix^D,b^C_)=w(ix^D,mom(idim))*w(ix^D,b^C_)-w(ix^D,mag(idim))*w(ix^D,m^C_)\
   {end do\}

    if(mhd_glm) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,mag(idim))=w(ix^D,psi_)
        !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
        f(ix^D,psi_)=cmax_global**2*w(ix^D,mag(idim))
     {end do\}
    end if
    if (mhd_fip) then
      f(ixO^S,fip_) = w(ixO^S,mom(idim)) * wC(ixO^S,fip_)
    end if
    ! Get flux of tracer
    do iw=1,mhd_n_tracer
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,tracer(iw))=w(ix^D,mom(idim))*w(ix^D,tracer(iw))
     {end do\}
    end do
  end subroutine mhd_get_flux_semirelati_noe

  !> Source term J.E_ambi in internal energy
  !> For the ambipolar electric field we have E_ambi = -eta_A * JxBxB= eta_A * B^2 (J_perpB)
  !> and eta_A is mhd_ambi_coef/rho^2 or is user-defined
  !> the source term J.E_ambi = eta_A * B^2 * J_perpB^2 = eta_A * [(JxB)xB]^2/B^2
  !> note that J_perpB= - (JxB)xB/B^2
  !> multiplyAmbiCoef is actually doing multiplication with -mhd_ambi_coef/rho^2
  subroutine add_source_ambipolar_internal_energy(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: tmp(ixI^S),btot2(ixI^S)
    double precision :: jxbxb(ixI^S,1:3)

    call mhd_get_jxbxb(wCT,x,ixI^L,ixO^L,jxbxb)
    ! avoiding nulls here
    btot2(ixO^S)=mhd_mag_en_all(wCT,ixI^L,ixO^L)
    where (btot2(ixO^S)>smalldouble )
       tmp(ixO^S) = sum(jxbxb(ixO^S,1:3)**2,dim=ndim+1) / btot2(ixO^S)
    elsewhere
       tmp(ixO^S) = zero
    endwhere
    call multiplyAmbiCoef(ixI^L,ixO^L,tmp,wCT,x)
    ! multiplyAmbiCoef is actually doing multiplication with -mhd_ambi_coef/rho^2
    ! hence minus sign here
    w(ixO^S,e_)=w(ixO^S,e_)- qdt*tmp(ixO^S)

  end subroutine add_source_ambipolar_internal_energy

  !> this subroutine computes -J_perpB= (J x B) x B= B(J.B) - J B^2
  subroutine mhd_get_jxbxb(w,x,ixI^L,ixO^L,res)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: res(ixI^S,1:3)

    double precision :: btot(ixI^S,1:3)
    double precision :: current(ixI^S,7-2*ndir:3)
    double precision :: tmp(ixI^S),b2(ixI^S)
    integer          :: idir, idirmin

    res=0.d0
    ! Calculate current density and idirmin
    ! current has nonzero values only for components in the range idirmin, 3
    call get_current(w,ixI^L,ixO^L,idirmin,current)
 
    btot=0.d0
    if(B0field) then
      do idir=1,ndir
        btot(ixO^S, idir) = w(ixO^S,mag(idir)) + block%B0(ixO^S,idir,b0i)
      enddo
    else
      do idir=1,ndir
        btot(ixO^S, idir) = w(ixO^S,mag(idir)) 
      enddo
    endif

    tmp(ixO^S)= sum(current(ixO^S,idirmin:3)*btot(ixO^S,idirmin:3),dim=ndim+1) !J.B
    b2(ixO^S) = sum(btot(ixO^S,1:3)**2,dim=ndim+1) !B^2
    do idir=1,idirmin-1
      res(ixO^S,idir) = btot(ixO^S,idir) * tmp(ixO^S)
    enddo
    do idir=idirmin,3
      res(ixO^S,idir) = btot(ixO^S,idir) * tmp(ixO^S) - current(ixO^S,idir) * b2(ixO^S)
    enddo

    ! avoid possible issues at nulls
    do idir=1,3
      where (b2(ixO^S)<smalldouble )
       res(ixO^S,idir) = zero
      endwhere
    enddo
  end subroutine mhd_get_jxbxb

  !> Sets the sources for the ambipolar terms for the STS method
  !> The sources are added directly (instead of fluxes as in the explicit)
  !> at the corresponding indices
  !> store_flux_var is explicitly called for each of the fluxes one by one
  subroutine sts_set_source_ambipolar(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: ixI^L,ixO^L,igrid,nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step

    double precision, dimension(ixI^S,1:3) :: tmp,ff
    double precision :: fluxall(ixI^S,1:nflux,1:ndim)
    double precision :: fE(ixI^S,sdim:3)
    double precision :: btot(ixI^S,1:3),tmp2(ixI^S)
    integer :: i, ixA^L, ie_

    ixA^L=ixO^L^LADD1;

    fluxall=zero

    ! here we compute (JxB)xB= - B^2 J_perpB
    call mhd_get_jxbxb(w,x,ixI^L,ixA^L,tmp)

    ! set ambipolar electric field in tmp: E_ambi = -eta_A * JxBxB= eta_A * B^2 (J_perpB)
    ! and eta_A is mhd_ambi_coef/rho^2 or is user-defined
    ! multiplyAmbiCoef is actually doing multiplication with -mhd_ambi_coef/rho^2
    do i=1,3
      call multiplyAmbiCoef(ixI^L,ixA^L,tmp(ixI^S,i),w,x)
    enddo

    ! Note: internal     energy case is handled through add_source_internal_e
    ! Note: hydrodynamic energy case is handled through add_source_hydrodynamic_e
    !       both of the above use     add_source_ambipolar_internal_energy
    ! 
    ! Note: total energy case without B0field split is ok here and adds div(BxE_ambi)
    ! Note: total energy case in semirelativistic variant (hence no B0field split) is ok here
    ! Note: total energy with B0field=T here adds div(B_1xE_ambi) which needs correction in add_source_B0split
    if(mhd_energy .and. .not.(mhd_internal_e.or.mhd_hydrodynamic_e))  then
      btot(ixA^S,1:3)    = 0.d0
      ! HERE: only uses B_1 if split, otherwise this is B
      btot(ixA^S,1:ndir) = w(ixA^S,mag(1:ndir))
      ! compute ff= E_ambi x B (where B can be B_1 if B0field=T)
      call cross_product(ixI^L,ixA^L,tmp,btot,ff)
      ! compute actual cell face fluxes in ff and their divergence in tmp2
      call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixI^S,1,1:ndim)=ff(ixI^S,1:ndim)
      ! - sign as the source is actually div(BxE_ambi) and we have div(E_ambi x B) in tmp2
      wres(ixO^S,e_)=-tmp2(ixO^S)
    endif

    if(stagger_grid) then
      ! always 2D or more (2.5/3D)
      if(ndir>ndim) then
        !!!Bz
        ff(ixA^S,1) = tmp(ixA^S,2)
        ff(ixA^S,2) = -tmp(ixA^S,1)
        ff(ixA^S,3) = 0.d0
        call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
        if(fix_conserve_at_step) fluxall(ixI^S,1+ndir,1:ndim)=ff(ixI^S,1:ndim)
        wres(ixO^S,mag(ndir))=-tmp2(ixO^S)
      end if
      fE=0.d0
      call update_faces_ambipolar(ixI^L,ixO^L,w,x,tmp,fE,btot)
      ixAmax^D=ixOmax^D;
      ixAmin^D=ixOmin^D-1;
      wres(ixA^S,mag(1:ndim))=-btot(ixA^S,1:ndim)
    else
      !write curl(ele) as the divergence
      !m1={0,ele[[3]],-ele[[2]]}
      !m2={-ele[[3]],0,ele[[1]]}
      !m3={ele[[2]],-ele[[1]],0}

      {^IFONED
      !!!Bx
      ff(ixA^S,1) = 0.d0
      ff(ixA^S,2) = tmp(ixA^S,3)
      ff(ixA^S,3) = -tmp(ixA^S,2)
      call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixI^S,2,1:ndim)=ff(ixI^S,1:ndim)
      !flux divergence is a source now
      wres(ixO^S,mag(1))=-tmp2(ixO^S)
      if(ndir==2.or.ndir==3)then
        !!!By
        ff(ixA^S,1) = -tmp(ixA^S,3)
        ff(ixA^S,2) = 0.d0
        ff(ixA^S,3) = tmp(ixA^S,1)
        call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
        if(fix_conserve_at_step) fluxall(ixI^S,3,1:ndim)=ff(ixI^S,1:ndim)
        wres(ixO^S,mag(2))=-tmp2(ixO^S)
      endif
      }
      {^NOONED
      !!!Bx
      ff(ixA^S,1) = 0.d0
      ff(ixA^S,2) = tmp(ixA^S,3)
      ff(ixA^S,3) = -tmp(ixA^S,2)
      call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixI^S,2,1:ndim)=ff(ixI^S,1:ndim)
      !flux divergence is a source now
      wres(ixO^S,mag(1))=-tmp2(ixO^S)
      !!!By
      ff(ixA^S,1) = -tmp(ixA^S,3)
      ff(ixA^S,2) = 0.d0
      ff(ixA^S,3) = tmp(ixA^S,1)
      call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixI^S,3,1:ndim)=ff(ixI^S,1:ndim)
      wres(ixO^S,mag(2))=-tmp2(ixO^S)
      }

      if(ndir==3) then
        !!!Bz
        ff(ixA^S,1) = tmp(ixA^S,2)
        ff(ixA^S,2) = -tmp(ixA^S,1)
        ff(ixA^S,3) = 0.d0
        call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
        if(fix_conserve_at_step) fluxall(ixI^S,1+ndir,1:ndim)=ff(ixI^S,1:ndim)
        wres(ixO^S,mag(ndir))=-tmp2(ixO^S)
      end if

    end if

    if(fix_conserve_at_step) then
      fluxall=my_dt*fluxall
      call store_flux(igrid,fluxall,1,ndim,nflux)
      if(stagger_grid) then
        call store_edge(igrid,ixI^L,my_dt*fE,1,ndim)
      end if
    end if

  end subroutine sts_set_source_ambipolar

  !> get ambipolar electric field and the integrals around cell faces
  subroutine update_faces_ambipolar(ixI^L,ixO^L,w,x,ECC,fE,circ)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    ! amibipolar electric field at cell centers
    double precision, intent(in)       :: ECC(ixI^S,1:3)
    double precision, intent(out)      :: fE(ixI^S,sdim:3)
    double precision, intent(out)      :: circ(ixI^S,1:ndim)

    integer                            :: hxC^L,ixC^L,ixA^L
    integer                            :: idim1,idim2,idir,ix^D

    fE=zero
    ! calculate ambipolar electric field on cell edges from cell centers
    do idir=sdim,3
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D+kr(idir,^D)-1;
     {do ix^DB=0,1\}
        if({ ix^D==1 .and. ^D==idir | .or.}) cycle
        ixAmin^D=ixCmin^D+ix^D;
        ixAmax^D=ixCmax^D+ix^D;
        fE(ixC^S,idir)=fE(ixC^S,idir)+ECC(ixA^S,idir)
     {end do\}
      fE(ixC^S,idir)=fE(ixC^S,idir)*0.25d0*block%dsC(ixC^S,idir)
    end do

    ! Calculate circulation on each face to get value of line integral of
    ! electric field in the positive idir direction.
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D-1;

    circ=zero
    do idim1=1,ndim ! Coordinate perpendicular to face
      do idim2=1,ndim
        do idir=sdim,3 ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
      circ(ixC^S,idim1)=circ(ixC^S,idim1)/block%surfaceC(ixC^S,idim1)
    end do

  end subroutine update_faces_ambipolar

  !> use cell-center flux vector to get cell-face flux vector
  !> which will be used to add the source term as the divergence of the flux
  !> we return fluxes at all faces as well as the divergence of the flux
  !> Note that for ndir>ndim, we do not modify the input cell center flux
  subroutine get_flux_on_cell_face(ixI^L,ixO^L,ff,src)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:3), intent(inout) :: ff
    double precision, intent(out) :: src(ixI^S)

    double precision :: ffc(ixI^S,1:ndim)
    double precision :: dxinv(ndim)
    integer :: idims, ix^D, ixA^L, ixB^L, ixC^L

    ixA^L=ixO^L^LADD1;
    dxinv=1.d0/dxlevel
    ! cell corner flux in ffc
    ! TO BE GENERALIZED FOR NON-UNIFORM NON-CARTESIAN MESH
    if (slab_uniform)then
      ffc=0.d0
      ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;
      {do ix^DB=0,1\}
        ixBmin^D=ixCmin^D+ix^D;
        ixBmax^D=ixCmax^D+ix^D;
        ffc(ixC^S,1:ndim)=ffc(ixC^S,1:ndim)+ff(ixB^S,1:ndim)
      {end do\}
      ffc(ixC^S,1:ndim)=0.5d0**ndim*ffc(ixC^S,1:ndim)
    else
      call mpistop("to generalize using volume averaging")
    endif
    ! now get flux at cell face from corner fluxes in fcc
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

  !> Calculates the explicit dt for the ambipolar term
  !> This function is used by both explicit scheme and STS method
  function get_ambipolar_dt(w,ixI^L,ixO^L,dx^D,x)  result(dtnew)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew

    double precision              :: coef
    double precision              :: dxarr(ndim)
    double precision              :: tmp(ixI^S)

    ^D&dxarr(^D)=dx^D;
    tmp(ixO^S) = mhd_mag_en_all(w, ixI^L, ixO^L)
    call multiplyAmbiCoef(ixI^L,ixO^L,tmp,w,x)
    ! now we have -mhd_eta_ambi B^2 /rho^2 in tmp
    coef = maxval(dabs(tmp(ixO^S)))
    if(coef/=0.d0) then
      coef=1.d0/coef
    else
      coef=bigdouble
    end if
    if(slab_uniform) then
      dtnew=minval(dxarr(1:ndim))**2.0d0*coef
    else
      dtnew=minval(block%ds(ixO^S,1:ndim))**2.0d0*coef
    end if

  end function get_ambipolar_dt

  !> multiply res by the ambipolar coefficient
  !> The ambipolar coefficient is calculated as -mhd_eta_ambi/rho^2
  !> The user may mask its value in the user file
  !> by implementing usr_mask_ambipolar subroutine
  subroutine multiplyAmbiCoef(ixI^L,ixO^L,res,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: res(ixI^S)
    double precision :: tmp(ixI^S)
    double precision :: rho(ixI^S)

    call mhd_get_rho(w,x,ixI^L,ixI^L,rho)
    tmp(ixI^S)=-mhd_eta_ambi/rho(ixI^S)**2
    if (associated(usr_mask_ambipolar)) then
      call usr_mask_ambipolar(ixI^L,ixO^L,w,x,tmp)
    end if
    res(ixO^S) = tmp(ixO^S) * res(ixO^S)

  end subroutine multiplyAmbiCoef

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine mhd_add_source(qdt,dtfactor,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_timing, only: time_htc_total, time_htc0
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source
    use mod_cak_force, only: cak_add_source

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,dtfactor
    double precision, intent(in)    :: wCT(ixI^S,1:nw),wCTprim(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active

    !TODO local_timestep support is only added for splitting
    ! but not for other nonideal terms such gravity, RC, viscosity,..
    ! it will also only work for divbfix  'linde', which does not require
    ! modification as it does not use dt in the update

    if (.not. qsourcesplit) then
      if(mhd_uawsom) then
        active = .true.
        call mhd_add_source_uawsom(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x)
      end if
      if(mhd_internal_e) then
        ! Source for solving internal energy
        active = .true.
        call add_source_internal_e(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
      else
        if(has_equi_rho_and_p) then
          active = .true.
          call add_equi_terms(qdt,dtfactor,ixI^L,ixO^L,wCT,w,x,wCTprim)
        end if
      end if

      if(mhd_hyperbolic_tc) then
        active = .true.
        call add_hyperbolic_tc_source(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
      end if

      ! Source for B0 splitting
      if (B0field) then
        active = .true.
        ! this adds source to momentum of type J0 x B0 and to energy equation
        ! latter always + J0 * E (electric field being E_ideal, E_hall, E_ambi)
        ! used for total energy variants
        call add_source_B0split(qdt,dtfactor,ixI^L,ixO^L,wCT,w,x,wCTprim)
      end if

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
      if (abs(mhd_eta)>smalldouble)then
        active = .true.
        call add_source_res_exp(qdt,ixI^L,ixO^L,wCT,w,x)
      end if

      if (mhd_ambipolar_exp)then
        active = .true.
        call add_source_ambi_exp(qdt,ixI^L,ixO^L,wCT,w,x)
      end if

      if (mhd_eta_hyper>0.d0)then
        active = .true.
        call add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
      end if

      if(mhd_hydrodynamic_e) then
        ! Source for solving hydrodynamic energy
        active = .true.
        call add_source_hydrodynamic_e(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
      else if (mhd_semirelativistic) then
        ! add sources for semirelativistic MHD
        active = .true.
        call add_source_semirelativistic(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
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

    if(mhd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,&
           w,x,qsourcesplit,active, rc_fl)
    end if

    if(mhd_viscosity) then
      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,&
           w,x,mhd_energy,qsourcesplit,active)
    end if

    if(mhd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,&
           w,x,gravity_energy,qsourcesplit,active)
    end if

    if (mhd_cak_force) then
      call cak_add_source(qdt,ixI^L,ixO^L,wCT,w,x,mhd_energy,qsourcesplit,active)
    end if

    ! This is where the radiation force and heating/cooling are added
    if (mhd_radiation_fld) then
       call mhd_add_radiation_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    endif

    ! update temperature from new pressure, density, and old ionization degree
    if(eos%eos_type == 'PI') then
      if(.not.qsourcesplit) then
        active = .true.
        call eos%update_eos(ixI^L,ixO^L,w,x)
      end if
    end if

  end subroutine mhd_add_source

  subroutine mhd_add_radiation_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods
    use mod_fld

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw),wCTprim(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: qsourcesplit
    logical, intent(inout) :: active

    ! add radiation force and work done by it, changes momentum and gas energy
    ! handle photon tiring, heating and cooling exchange between gas and radiation field
    call add_fld_rad_force(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active,fld_fl)

  end subroutine mhd_add_radiation_source

  !> add some source terms to total energy related to has_equi_rho_and_p=T
  subroutine add_equi_terms(qdt,dtfactor,ixI^L,ixO^L,wCT,w,x,wCTprim)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,dtfactor
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCTprim(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision                :: divv(ixI^S)
    double precision :: a(ixI^S,3), b(ixI^S,3), axb(ixI^S,3)
    double precision :: gravity_field(ixI^S,1:ndim)
    integer :: idir

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(wCTprim(ixI^S,mom(1:ndir)),ixI^L,ixO^L,divv,3)
      else
        call divvector(wCTprim(ixI^S,mom(1:ndir)),ixI^L,ixO^L,divv,2)
      end if
    else
     call divvector(wCTprim(ixI^S,mom(1:ndir)),ixI^L,ixO^L,divv)
    end if
    divv(ixO^S)=divv(ixO^S)*eos%gamma*eos%inv_gamma_minus_1
    if(local_timestep) then
      w(ixO^S,e_)=w(ixO^S,e_)-dtfactor*block%dt(ixO^S)*block%equi_vars(ixO^S,equi_pe0_,0)*divv(ixO^S)
    else
      w(ixO^S,e_)=w(ixO^S,e_)-qdt*block%equi_vars(ixO^S,equi_pe0_,0)*divv(ixO^S)
    end if
    if(B0field)then
      if(B0field_forcefree.and.mhd_gravity)then
        ! add -v dot(rho_0 g)/(gamma-1)
        call usr_gravity(ixI^L,ixO^L,wCT,x,gravity_field)
        do idir=1,ndim
           w(ixO^S,e_)=w(ixO^S,e_)-qdt*wCTprim(ixO^S,mom(idir))*block%equi_vars(ixO^S,equi_rho0_,0)*gravity_field(ixO^S,idir)*eos%inv_gamma_minus_1
        enddo
      else
        a=0.d0
        b=0.d0
        ! store B0 magnetic field in b
        b(ixO^S,1:ndir)=block%B0(ixO^S,1:ndir,0)
        ! store J0 current in a
        do idir=7-2*ndir,3
          a(ixO^S,idir)=block%J0(ixO^S,idir)
        end do
        call cross_product(ixI^L,ixO^L,a,b,axb)
        ! add -v dot(rho_0 g + J0 x B_0)/(gamma-1)
        do idir=1,ndir
           w(ixO^S,e_)=w(ixO^S,e_)-qdt*wCTprim(ixO^S,mom(idir))*axb(ixO^S,idir)*eos%inv_gamma_minus_1
        enddo
        if(mhd_gravity)then
          ! add -v dot(rho_0 g)/(gamma-1)
          call usr_gravity(ixI^L,ixO^L,wCT,x,gravity_field)
          do idir=1,ndim
           w(ixO^S,e_)=w(ixO^S,e_)-qdt*wCTprim(ixO^S,mom(idir))*block%equi_vars(ixO^S,equi_rho0_,0)*gravity_field(ixO^S,idir)*eos%inv_gamma_minus_1
          enddo
        endif
      endif
    else
      if(mhd_gravity)then
        ! add -v dot(rho_0 g)/(gamma-1)
        call usr_gravity(ixI^L,ixO^L,wCT,x,gravity_field)
        do idir=1,ndim
         w(ixO^S,e_)=w(ixO^S,e_)-qdt*wCTprim(ixO^S,mom(idir))*block%equi_vars(ixO^S,equi_rho0_,0)*gravity_field(ixO^S,idir)*eos%inv_gamma_minus_1
        enddo
      endif
    endif
  end subroutine add_equi_terms

  subroutine mhd_get_hyperbolic_tc_geometry(ixI^L,ixO^L,Te,Bvec,bgradT,gradTperp_mag,nperp)
    use mod_global_parameters
    use mod_geometry, only: gradient
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: Te(ixI^S)
    double precision, intent(in) :: Bvec(ixI^S,1:ndir)
    double precision, intent(out) :: bgradT(ixI^S), gradTperp_mag(ixI^S)
    double precision, intent(out) :: nperp(ixI^S,1:ndir)

    double precision :: Bmag, bunitvec(ndir), gradT(ndir), gradT_perp(ndir)
    double precision :: gradT_cell(ixI^S,1:ndir)
    integer :: ix^D, idir

    gradT_cell=zero
    if(.not. slab_uniform) then
      do idir=1,ndim
        call gradient(Te,ixI^L,ixO^L,idir,gradT_cell(ixI^S,idir))
      end do
    end if

   {^IFTWOD
    do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        Bmag=zero
        do idir=1,ndir
          Bmag=Bmag+Bvec(ix^D,idir)**2
        end do
        Bmag=dsqrt(Bmag)

        if(Bmag>smalldouble) then
          do idir=1,ndir
            bunitvec(idir)=Bvec(ix^D,idir)/Bmag
          end do
        else
          do idir=1,ndir
            bunitvec(idir)=zero
          end do
        end if
        if(slab_uniform) then
          gradT(1)=((8.d0*(Te(ix1+1,ix2)-Te(ix1-1,ix2))-Te(ix1+2,ix2)+Te(ix1-2,ix2))/12.d0)/block%ds(ix^D,1)
          gradT(2)=((8.d0*(Te(ix1,ix2+1)-Te(ix1,ix2-1))-Te(ix1,ix2+2)+Te(ix1,ix2-2))/12.d0)/block%ds(ix^D,2)
          if(ndir>2) gradT(3)=zero
        else
          do idir=1,ndir
            gradT(idir)=gradT_cell(ix^D,idir)
          end do
        end if

        bgradT(ix^D)=zero
        do idir=1,ndir
          bgradT(ix^D)=bgradT(ix^D)+bunitvec(idir)*gradT(idir)
        end do

        do idir=1,ndir
          gradT_perp(idir)=gradT(idir)-bgradT(ix^D)*bunitvec(idir)
        end do

        gradTperp_mag(ix^D)=zero
        do idir=1,ndir
          gradTperp_mag(ix^D)=gradTperp_mag(ix^D)+gradT_perp(idir)**2
        end do
        gradTperp_mag(ix^D)=dsqrt(gradTperp_mag(ix^D))

        if(gradTperp_mag(ix^D)>smalldouble) then
          do idir=1,ndir
            nperp(ix^D,idir)=gradT_perp(idir)/gradTperp_mag(ix^D)
          end do
        else
          gradTperp_mag(ix^D)=zero
          do idir=1,ndir
            nperp(ix^D,idir)=zero
          end do
        end if
      end do
    end do
   }
   {^IFTHREED
    do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
        do ix1=ixOmin1,ixOmax1
          Bmag=dsqrt(Bvec(ix^D,1)**2+Bvec(ix^D,2)**2+Bvec(ix^D,3)**2)
          if(Bmag>smalldouble) then
            bunitvec(1)=Bvec(ix^D,1)/Bmag
            bunitvec(2)=Bvec(ix^D,2)/Bmag
            bunitvec(3)=Bvec(ix^D,3)/Bmag
          else
            bunitvec(1)=zero
            bunitvec(2)=zero
            bunitvec(3)=zero
          end if

          if(slab_uniform) then
            gradT(1)=((8.d0*(Te(ix1+1,ix2,ix3)-Te(ix1-1,ix2,ix3))-Te(ix1+2,ix2,ix3)+Te(ix1-2,ix2,ix3))/12.d0)/block%ds(ix^D,1)
            gradT(2)=((8.d0*(Te(ix1,ix2+1,ix3)-Te(ix1,ix2-1,ix3))-Te(ix1,ix2+2,ix3)+Te(ix1,ix2-2,ix3))/12.d0)/block%ds(ix^D,2)
            gradT(3)=((8.d0*(Te(ix1,ix2,ix3+1)-Te(ix1,ix2,ix3-1))-Te(ix1,ix2,ix3+2)+Te(ix1,ix2,ix3-2))/12.d0)/block%ds(ix^D,3)
          else
            do idir=1,ndir
              gradT(idir)=gradT_cell(ix^D,idir)
            end do
          end if

          bgradT(ix^D)=zero
          do idir=1,ndir
            bgradT(ix^D)=bgradT(ix^D)+bunitvec(idir)*gradT(idir)
          end do

          do idir=1,ndir
            gradT_perp(idir)=gradT(idir)-bgradT(ix^D)*bunitvec(idir)
          end do

          gradTperp_mag(ix^D)=dsqrt(gradT_perp(1)**2+gradT_perp(2)**2+gradT_perp(3)**2)
          if(gradTperp_mag(ix^D)>smalldouble) then
            do idir=1,ndir
              nperp(ix^D,idir)=gradT_perp(idir)/gradTperp_mag(ix^D)
            end do
          else
            gradTperp_mag(ix^D)=zero
            do idir=1,ndir
              nperp(ix^D,idir)=zero
            end do
          end if
        end do
      end do
    end do
    }
  end subroutine mhd_get_hyperbolic_tc_geometry

  subroutine add_hyperbolic_tc_source(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
    use mod_global_parameters
    use mod_geometry, only: gradient
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qdt
    double precision, dimension(ixI^S,1:ndim), intent(in) :: x
    double precision, dimension(ixI^S,1:nw), intent(in) :: wCT,wCTprim
    double precision, dimension(ixI^S,1:nw), intent(inout) :: w

    double precision, dimension(ixI^S) :: R,Te,rho_loc,pth_loc
    double precision, dimension(ixI^S) :: ne_loc,nH_dummy
    double precision, dimension(ixI^S,1:ndir) :: Bvec
    double precision, dimension(ixI^S) :: bgradT, gradTperp_mag
    double precision, dimension(ixI^S,1:ndir) :: nperp
    double precision, dimension(ixI^S) :: gradT_geom
    double precision, parameter :: xe_prefac_cgs = 4.753567596681522d6
    double precision :: kappa_T5,kappa_T5_perp,kappa_T5_perp_eff
    double precision :: kappa_T7,f_sat,kappaT5_bgradT,kappaT5_gradTperp,tau,B2,fB,gradT1
    double precision :: qclass_diss
    double precision :: Bmag_loc,Tloc,Tcond,nloc_code,Cchi,chi
    double precision :: cmax(ndim),c2,cfast2,avMinCs2(ndim),inv_rho
    logical :: use_perp_source
    integer :: ix^D,idir

    ! xe_prefac_cgs expects B [G], T [K], and ne [cm^-3]. Convert the
    ! normalisation units explicitly in SI mode, while retaining code-unit
    ! state variables in the cell loop below.
    Cchi=zero
    if(trim(mhd_hyperbolic_tc_perp_mode)=='electron_magnetization') then
      if(SI_unit) then
        Cchi = (xe_prefac_cgs/mhd_hyperbolic_tc_coulomb_log) * &
               (1.d4*unit_magneticfield)*unit_temperature**1.5d0 / &
               (1.d-6*unit_numberdensity)
      else
        Cchi = (xe_prefac_cgs/mhd_hyperbolic_tc_coulomb_log) * &
               unit_magneticfield*unit_temperature**1.5d0/unit_numberdensity
      end if
    end if
    call eos%get_Rfactor(wCT,x,ixI^L,ixI^L,R)
    {do ix^DB=ixImin^DB,ixImax^DB\}
      if(has_equi_rho_and_p) then
        rho_loc(ix^D)=wCTprim(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0)
        pth_loc(ix^D)=wCTprim(ix^D,p_)+block%equi_vars(ix^D,equi_pe0_,0)
      else
        rho_loc(ix^D)=wCTprim(ix^D,rho_)
        pth_loc(ix^D)=wCTprim(ix^D,p_)
      end if
      Te(ix^D)=pth_loc(ix^D)/(R(ix^D)*rho_loc(ix^D))
    {end do\}
    ! Electron number density in code units. Under the FI/PI normalisation,
    ! rho_code equals nH_code. FI uses the composition-dependent fully ionised
    ! electron count. PI obtains the local electron count from its R factor:
    ! R*(2+3 A_He) = (n_nuclei+n_e)/n_H. LTE stores ne explicitly.
    if(trim(mhd_hyperbolic_tc_perp_mode)=='electron_magnetization') then
      if(eos%eos_type=='LTE') then
        call eos%get_ne_nH(ixI^L,ixI^L,wCT, x,ne_loc,nH_dummy)
      else if(eos%eos_type=='PI') then
        ne_loc(ixI^S)=rho_loc(ixI^S)*max(R(ixI^S)*(2.d0+3.d0*eos%He_abundance) &
             -(1.d0+eos%He_abundance),smalldouble)
      else
        ne_loc(ixI^S)=rho_loc(ixI^S)*(1.d0+2.d0*eos%He_abundance)
      end if
    end if
    use_perp_source = mhd_hyperbolic_tc_use_perp .and. &
                      trim(mhd_hyperbolic_tc_perp_mode)/='off'
    if(B0field) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        do idir=1,ndir
          Bvec(ix^D,idir)=wCT(ix^D,mag(idir))+block%B0(ix^D,idir,0)
        end do
      {end do\}
    else
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        do idir=1,ndir
          Bvec(ix^D,idir)=wCT(ix^D,mag(idir))
        end do
      {end do\}
    end if
   {^NOONED
    call mhd_get_hyperbolic_tc_geometry(ixI^L,ixO^L,Te,Bvec,bgradT,gradTperp_mag,nperp)
   }
   {^IFONED
    gradT_geom=zero
    if(.not.slab_uniform) then
      call gradient(Te,ixI^L,ixO^L,1,gradT_geom)
    end if
    do ix1=ixOmin1,ixOmax1
      if(mhd_hyperbolic_tc_constant) then
        kappa_T5=mhd_hyperbolic_tc_kappa
        kappa_T7=kappa_T5*Te(ix1)
      else
        Tcond = Te(ix1)
        if(mhd_trac) then
          Tcond = max(Tcond, block%wextra(ix1,Tcoff_))
        end if
        kappa_T5=mhd_hyperbolic_tc_kappa*sqrt(Tcond**5)
        kappa_T7=kappa_T5*Tcond
      end if
      if(slab_uniform) then
        gradT1=((8.d0*(Te(ix1+1)-Te(ix1-1))-Te(ix1+2)+Te(ix1-2))/12.d0)/block%ds(ix1,1)
      else
        gradT1=gradT_geom(ix1)
      end if
      B2=zero
      do idir=1,ndir
        B2=B2+Bvec(ix1,idir)**2
      end do
      if(B2>smalldouble**2) then
        bgradT(ix1)=Bvec(ix1,1)*gradT1/dsqrt(B2)
      else
        bgradT(ix1)=zero
      end if
      kappaT5_bgradT=kappa_T5*bgradT(ix1)
      inv_rho=1.d0/rho_loc(ix1)
      c2=eos%gamma*pth_loc(ix1)*inv_rho
      cfast2 = B2*inv_rho + c2
      avMinCs2(1) = cfast2**2 - 4.0d0*c2*Bvec(ix1,1)**2*inv_rho
      cmax(1) = sqrt(half*(cfast2 + sqrt(dabs(avMinCs2(1)))))
      if(mhd_hyperbolic_tc_sat) then
        f_sat=one/(one+dabs(kappaT5_bgradT)/(1.5d0*rho_loc(ix^D)*(pth_loc(ix^D)/rho_loc(ix^D))**1.5d0))
        tau=max(4.d0*dt, f_sat*kappa_T7*courantpar**2/(pth_loc(ix^D)*eos%inv_gamma_minus_1*cmax(1)**2))
        w(ix^D,qpar_)=w(ix^D,qpar_)-qdt*(f_sat*kappaT5_bgradT+wCT(ix^D,qpar_))/tau
      else
        w(ix^D,qpar_)=w(ix^D,qpar_)-qdt*(kappaT5_bgradT+wCT(ix^D,qpar_))/&
         max(4.d0*dt, kappa_T7*courantpar**2/(pth_loc(ix^D)*eos%inv_gamma_minus_1*cmax(1)**2))
      end if
    end do
   }
   {^IFTWOD
    do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        if(mhd_hyperbolic_tc_constant) then
          kappa_T5=mhd_hyperbolic_tc_kappa
          kappa_T7=kappa_T5*Te(ix^D)
        else
          Tcond=Te(ix^D)
          if(mhd_trac) then
            Tcond=max(Tcond, block%wextra(ix^D,Tcoff_))
          end if
          kappa_T5=mhd_hyperbolic_tc_kappa*sqrt(Tcond**5)
          kappa_T7 = kappa_T5*Tcond
        end if
        kappaT5_bgradT=kappa_T5*bgradT(ix^D)
        B2 = zero
        do idir = 1, ndir
          B2 = B2 + Bvec(ix^D,idir)**2
        end do
        if(use_perp_source) then
          select case(trim(mhd_hyperbolic_tc_perp_mode))
          case('fixed_reference')
            kappa_T5_perp=mhd_hyperbolic_tc_kappa_perp_factor*kappa_T5
          case('weak_field_isotropization')
            if(mhd_hyperbolic_tc_Bmin>zero) then
              fB=B2/(B2+mhd_hyperbolic_tc_Bmin**2)
            else
              fB=one
            end if
            kappa_T5_perp_eff=(one-fB)*kappa_T5
            kappa_T5_perp=kappa_T5_perp_eff
          case('electron_magnetization')
            Bmag_loc = dsqrt(B2)
            Tloc = max(Te(ix^D), smalldouble)
            nloc_code = max(ne_loc(ix^D), smalldouble)
            chi = Cchi*Bmag_loc*Tloc**1.5d0/nloc_code
            kappa_T5_perp_eff = kappa_T5/(one+chi**2)
            kappa_T5_perp = kappa_T5_perp_eff
          case default
            kappa_T5_perp=zero
          end select
          kappaT5_gradTperp=kappa_T5_perp*gradTperp_mag(ix^D)
        end if
        inv_rho=1.d0/rho_loc(ix^D)
        c2=eos%gamma*pth_loc(ix^D)*inv_rho
        cfast2 = B2*inv_rho + c2
        do idir=1,ndim
          avMinCs2(idir)=cfast2**2-4.0d0*c2*Bvec(ix^D,idir)**2*inv_rho
          cmax(idir)=sqrt(half*(cfast2+sqrt(dabs(avMinCs2(idir)))))\
        end do
        if(mhd_hyperbolic_tc_sat) then
          qclass_diss=dabs(kappaT5_bgradT)
          if(use_perp_source) &
            qclass_diss=dsqrt(kappaT5_bgradT**2+kappaT5_gradTperp**2)
          f_sat=one/(one+qclass_diss/(1.5d0*rho_loc(ix^D)*(pth_loc(ix^D)/rho_loc(ix^D))**1.5d0))
          tau=max(4.d0*dt, f_sat*kappa_T7*courantpar**2/(pth_loc(ix^D)*eos%inv_gamma_minus_1*maxval(cmax(:))**2))
          w(ix^D,qpar_)=w(ix^D,qpar_)-qdt*(f_sat*kappaT5_bgradT+wCT(ix^D,qpar_))/tau
          if(use_perp_source) then
            w(ix^D,qperp_)=w(ix^D,qperp_)-qdt*(f_sat*kappaT5_gradTperp+wCT(ix^D,qperp_))/tau
          end if
        else
          tau=max(4.d0*dt, kappa_T7*courantpar**2/(pth_loc(ix^D)*eos%inv_gamma_minus_1*maxval(cmax(:))**2))
          w(ix^D,qpar_)=w(ix^D,qpar_)-qdt*(kappaT5_bgradT+wCT(ix^D,qpar_))/tau
          if(use_perp_source) then
            w(ix^D,qperp_)=w(ix^D,qperp_)-qdt*(kappaT5_gradTperp+wCT(ix^D,qperp_))/tau
          end if
        end if
      end do
    end do
   }
   {^IFTHREED
    do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
        do ix1=ixOmin1,ixOmax1
          if(mhd_hyperbolic_tc_constant) then
            kappa_T5=mhd_hyperbolic_tc_kappa
            kappa_T7=kappa_T5*Te(ix^D)
          else
            Tcond=Te(ix^D)
            if(mhd_trac) then
              Tcond=max(Tcond, block%wextra(ix^D,Tcoff_))
            end if
            kappa_T5=mhd_hyperbolic_tc_kappa*sqrt(Tcond**5)
            kappa_T7 = kappa_T5*Tcond
          end if
          kappaT5_bgradT=kappa_T5*bgradT(ix^D)
          B2 = zero
          do idir = 1, ndir
            B2 = B2 + Bvec(ix^D,idir)**2
          end do
          if(use_perp_source) then
            select case(trim(mhd_hyperbolic_tc_perp_mode))
            case('fixed_reference')
              kappa_T5_perp=mhd_hyperbolic_tc_kappa_perp_factor*kappa_T5
            case('weak_field_isotropization')
              if(mhd_hyperbolic_tc_Bmin>zero) then
                fB=B2/(B2+mhd_hyperbolic_tc_Bmin**2)
              else
                fB=one
              end if
              kappa_T5_perp_eff=(one-fB)*kappa_T5
              kappa_T5_perp=kappa_T5_perp_eff
            case('electron_magnetization')
              Bmag_loc = dsqrt(B2)
              Tloc = max(Te(ix^D), smalldouble)
              nloc_code = max(ne_loc(ix^D), smalldouble)
              chi = Cchi*Bmag_loc*Tloc**1.5d0/nloc_code
              kappa_T5_perp_eff = kappa_T5/(one+chi**2)
              kappa_T5_perp = kappa_T5_perp_eff
            case default
              kappa_T5_perp=zero
            end select
            kappaT5_gradTperp=kappa_T5_perp*gradTperp_mag(ix^D)
          end if
          inv_rho=1.d0/rho_loc(ix^D)
          c2=eos%gamma*pth_loc(ix^D)*inv_rho
          cfast2 = B2*inv_rho + c2
          do idir = 1, ndim
            avMinCs2(idir)=cfast2**2-4.0d0*c2*Bvec(ix^D,idir)**2*inv_rho
            cmax(idir)=sqrt(half*(cfast2+sqrt(dabs(avMinCs2(idir)))))\
          end do
          if(mhd_hyperbolic_tc_sat) then
            qclass_diss=dabs(kappaT5_bgradT)
            if(use_perp_source) &
              qclass_diss=dsqrt(kappaT5_bgradT**2+kappaT5_gradTperp**2)
            f_sat=one/(one+qclass_diss/(1.5d0*rho_loc(ix^D)*(pth_loc(ix^D)/rho_loc(ix^D))**1.5d0))
            tau=max(4.d0*dt, f_sat*kappa_T7*courantpar**2/(pth_loc(ix^D)*eos%inv_gamma_minus_1*maxval(cmax(:))**2))
            w(ix^D,qpar_)=w(ix^D,qpar_)-qdt*(f_sat*kappaT5_bgradT+wCT(ix^D,qpar_))/tau
            if(use_perp_source) then
              w(ix^D,qperp_)=w(ix^D,qperp_)-qdt*(f_sat*kappaT5_gradTperp+wCT(ix^D,qperp_))/tau
            end if
          else
            tau=max(4.d0*dt, kappa_T7*courantpar**2/(pth_loc(ix^D)*eos%inv_gamma_minus_1*maxval(cmax(:))**2))
            w(ix^D,qpar_)=w(ix^D,qpar_)-qdt*(kappaT5_bgradT+wCT(ix^D,qpar_))/tau
            if(use_perp_source) then
              w(ix^D,qperp_)=w(ix^D,qperp_)-qdt*(kappaT5_gradTperp+wCT(ix^D,qperp_))/tau
            end if
          end if
        end do
      end do
    end do
   }
  end subroutine add_hyperbolic_tc_source

  !> Compute the Lorentz force (JxB) Note: Unused subroutine
  !> perhaps useful for post-processing when made public
  subroutine get_Lorentz_force(ixI^L,ixO^L,w,JxB)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: JxB(ixI^S,3)
    double precision                :: a(ixI^S,3), b(ixI^S,3)
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3)
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

  subroutine mhd_get_rho(w,x,ixI^L,ixO^L,rho)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: rho(ixI^S)

    if(has_equi_rho_and_p) then
      rho(ixO^S) = w(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,b0i)
    else
      rho(ixO^S) = w(ixO^S,rho_)
    endif

  end subroutine mhd_get_rho

  !> handle small or negative internal energy
  subroutine mhd_handle_small_ei(w, x, ixI^L, ixO^L, ie, subname)
    use mod_global_parameters
    use mod_small_values
    integer, intent(in)             :: ixI^L,ixO^L, ie
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    double precision              :: rho(ixI^S)
    integer :: idir
    logical :: flag(ixI^S,1:nw)

    flag=.false.
    if(has_equi_rho_and_p) then
      where(w(ixO^S,ie)+block%equi_vars(ixO^S,equi_pe0_,0)*eos%inv_gamma_minus_1<small_e)&
             flag(ixO^S,ie)=.true.
    else
      where(w(ixO^S,ie)<small_e) flag(ixO^S,ie)=.true.
    endif
    if(any(flag(ixO^S,ie))) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_rho_and_p) then
          where(flag(ixO^S,ie)) w(ixO^S,ie)=small_e - &
                  block%equi_vars(ixO^S,equi_pe0_,0)*eos%inv_gamma_minus_1
        else
          where(flag(ixO^S,ie)) w(ixO^S,ie)=small_e
        endif
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag, ie)
      case default
        ! small values error shows primitive variables
        w(ixO^S,e_)=w(ixO^S,e_)*eos%gamma_minus_1
        call mhd_get_rho(w,x,ixI^L,ixO^L,rho)
        do idir = 1, ndir
           w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))/rho(ixO^S)
        end do
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if

  end subroutine mhd_handle_small_ei


  !> Source terms after split off time-independent magnetic field
  subroutine add_source_B0split(qdt,dtfactor,ixI^L,ixO^L,wCT,w,x,wCTprim)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, dtfactor,wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(in) :: wCTprim(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: a(ixI^S,3), b(ixI^S,3), axb(ixI^S,3)
    integer :: idir

    a=0.d0
    b=0.d0
    ! for force-free field J0xB0 =0
    if((.not.B0field_forcefree).and.(.not.has_equi_rho_and_p)) then
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
      b(ixO^S,:)=wCTprim(ixO^S,mag(:))
      ! store full magnetic field B0+B1 in b
      if((.not.B0field_forcefree).and.(.not.has_equi_rho_and_p)) b(ixO^S,:)=b(ixO^S,:)+block%B0(ixO^S,:,0)
      ! store velocity in a
      a(ixI^S,1:ndir)=wCTprim(ixI^S,mom(1:ndir))
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
      ! where it is adding -J0 dot (vxB_1) when appropriate
      do idir=7-2*ndir,3
        w(ixO^S,e_)=w(ixO^S,e_)-axb(ixO^S,idir)*block%J0(ixO^S,idir)
      end do
      if(mhd_hall) then
         ! store hall velocity in a, only partial current is needed
         call mhd_getv_Hall(wCT,x,ixI^L,ixO^L,a,.true.)
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
      endif
      if(mhd_ambipolar_sts) then
        ! in STS variant of ambipolar, we added for split B the term div(B_1xE_ambi)
        ! hence needs to add J_0 dot E_ambi 
        ! to get finally the term etaA (J_perpB)^/B^2-B_1 dot (curl Eambi)
        !reuse axb
        call mhd_get_jxbxb(wCT,x,ixI^L,ixO^L,axb)
        ! source J0 * E
        do idir=sdim,3
          !set electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
          call multiplyAmbiCoef(ixI^L,ixO^L,axb(ixI^S,idir),wCT,x)
          w(ixO^S,e_)=w(ixO^S,e_)+qdt*axb(ixO^S,idir)*block%J0(ixO^S,idir)
        enddo
      endif
    end if
    

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_B0')

  end subroutine add_source_B0split

  !> Source terms for semirelativistic MHD Gombosi 2002 JCP 177, 176
  subroutine add_source_semirelativistic(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in), optional :: wCTprim(ixI^S,1:nw)

    double precision :: E(ixI^S,1:3),curlE(ixI^S,1:3),divE(ixI^S)
    integer :: idir, idirmin, ix^D

   ! if ndir<3 the source is zero
   {^IFTHREEC
   {do ix^DB=ixImin^DB,ixImax^DB\}
      ! E=Bxv
      E(ix^D,1)=w(ix^D,b2_)*wCTprim(ix^D,m3_)-w(ix^D,b3_)*wCTprim(ix^D,m2_)
      E(ix^D,2)=w(ix^D,b3_)*wCTprim(ix^D,m1_)-w(ix^D,b1_)*wCTprim(ix^D,m3_)
      E(ix^D,3)=w(ix^D,b1_)*wCTprim(ix^D,m2_)-w(ix^D,b2_)*wCTprim(ix^D,m1_)
   {end do\}
    call divvector(E,ixI^L,ixO^L,divE)
    ! curl E
    call curlvector(E,ixI^L,ixO^L,curlE,idirmin,1,3)
    ! add source term in momentum equations (1/c0^2-1/c^2)(E divE - E x curlE)
    ! equation (26) and (27)
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      w(ix^D,m1_)=w(ix^D,m1_)+qdt*(inv_squared_c0-inv_squared_c)*&
       (E(ix^D,1)*divE(ix^D)-E(ix^D,2)*curlE(ix^D,3)+E(ix^D,3)*curlE(ix^D,2))
      w(ix^D,m2_)=w(ix^D,m2_)+qdt*(inv_squared_c0-inv_squared_c)*&
       (E(ix^D,2)*divE(ix^D)-E(ix^D,3)*curlE(ix^D,1)+E(ix^D,1)*curlE(ix^D,3))
      w(ix^D,m3_)=w(ix^D,m3_)+qdt*(inv_squared_c0-inv_squared_c)*&
       (E(ix^D,3)*divE(ix^D)-E(ix^D,1)*curlE(ix^D,2)+E(ix^D,2)*curlE(ix^D,1) )
   {end do\}
   }

  end subroutine add_source_semirelativistic

  !> Source terms for internal energy version of MHD
  subroutine add_source_internal_e(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: wCTprim(ixI^S,1:nw)

    double precision                :: divv(ixI^S), tmp
    integer :: ix^D

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(wCTprim(ixI^S,mom(:)),ixI^L,ixO^L,divv,3)
      else
        call divvector(wCTprim(ixI^S,mom(:)),ixI^L,ixO^L,divv,2)
      end if
    else
      call divvector(wCTprim(ixI^S,mom(:)),ixI^L,ixO^L,divv)
    end if
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      tmp=w(ix^D,e_)
      w(ix^D,e_)=w(ix^D,e_)-qdt*wCTprim(ix^D,p_)*divv(ix^D)
      if(w(ix^D,e_)<small_e) then
        w(ix^D,e_)=tmp
      end if
   {end do\}
    if(mhd_ambipolar_sts)then
      call add_source_ambipolar_internal_energy(qdt,ixI^L,ixO^L,wCT,w,x)
    end if

    if(fix_small_values) then
      call mhd_handle_small_ei(w,x,ixI^L,ixO^L,e_,'add_source_internal_e')
    end if
  end subroutine add_source_internal_e

  !> Source terms for hydrodynamic energy version of MHD
  subroutine add_source_hydrodynamic_e(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_gravity

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in), optional :: wCTprim(ixI^S,1:nw)

    double precision :: B(ixI^S,3), J(ixI^S,3), JxB(ixI^S,3)
    double precision :: current(ixI^S,7-2*ndir:3)
    double precision :: bu(ixO^S,1:ndir), tmp(ixO^S), b2(ixO^S)
    double precision :: gravity_field(ixI^S,1:ndir), Vaoc
    integer :: idir, idirmin, idims, ix^D

    {^NOTHREED
    B=0.0d0
    do idir = 1, ndir
      B(ixO^S, idir) = wCT(ixO^S,mag(idir))
    end do

    if(slab_uniform)then
      ! get current in fourth order accuracy in Cartesian
      call curlvector(wCT(ixI^S,mag(1:ndir)),ixI^L,ixO^L,current,idirmin,7-2*ndir,ndir,.true.)
    else
      call get_current(wCT,ixI^L,ixO^L,idirmin,current)
    endif

    J=0.0d0
    do idir=7-2*ndir,3
      J(ixO^S,idir)=current(ixO^S,idir)
    end do

    ! get Lorentz force JxB
    call cross_product(ixI^L,ixO^L,J,B,JxB)
    }
    {^IFTHREED
    if(slab_uniform)then
      ! get current in fourth order accuracy in Cartesian
      call curlvector(wCT(ixI^S,mag(1:ndir)),ixI^L,ixO^L,current,idirmin,1,ndir,.true.)
    else
      call get_current(wCT,ixI^L,ixO^L,idirmin,current)
    endif
    ! get Lorentz force JxB
    call cross_product(ixI^L,ixO^L,current,wCT(ixI^S,mag(1:ndir)),JxB)
    }

    ! mhd_semirelativistic does not combine with mhd_hydrodynamic_e
    !!if(mhd_semirelativistic) then
    !!  ! (v . nabla) v
    !!  do idir=1,ndir
    !!    do idims=1,ndim
    !!      call gradient(wCTprim(ixI^S,mom(idir)),ixI^L,ixO^L,idims,J(ixI^S,idims))
    !!    end do
    !!    B(ixO^S,idir)=sum(wCTprim(ixO^S,mom(1:ndir))*J(ixO^S,1:ndir),dim=ndim+1)
    !!  end do
    !!  ! nabla p
    !!  do idir=1,ndir
    !!    call gradient(wCTprim(ixI^S,p_),ixI^L,ixO^L,idir,J(ixI^S,idir))
    !!  end do
    !!  if(mhd_gravity) then
    !!    gravity_field=0.d0
    !!    call usr_gravity(ixI^L,ixO^L,wCT,x,gravity_field(ixI^S,1:ndim))
    !!    do idir=1,ndir
    !!      B(ixO^S,idir)=wCT(ixO^S,rho_)*(B(ixO^S,idir)-gravity_field(ixO^S,idir))+J(ixO^S,idir)-JxB(ixO^S,idir)
    !!    end do
    !!  else
    !!    do idir=1,ndir
    !!      B(ixO^S,idir)=wCT(ixO^S,rho_)*B(ixO^S,idir)+J(ixO^S,idir)-JxB(ixO^S,idir)
    !!    end do
    !!  end if
    !!  b2(ixO^S)=sum(wCT(ixO^S,mag(:))**2,dim=ndim+1)
    !!  tmp(ixO^S)=sqrt(b2(ixO^S))
    !!  where(tmp(ixO^S)>smalldouble)
    !!    tmp(ixO^S)=1.d0/tmp(ixO^S)
    !!  else where
    !!    tmp(ixO^S)=0.d0
    !!  end where
    !!  ! unit vector of magnetic field
    !!  do idir=1,ndir
    !!    bu(ixO^S,idir)=wCT(ixO^S,mag(idir))*tmp(ixO^S)
    !!  end do
    !!  !b2(ixO^S)=b2(ixO^S)/w(ixO^S,rho_)*inv_squared_c
    !!  !b2(ixO^S)=b2(ixO^S)/(1.d0+b2(ixO^S))
    !!  {do ix^DB=ixOmin^DB,ixOmax^DB\}
    !!     ! Va^2/c^2
    !!     Vaoc=b2(ix^D)/w(ix^D,rho_)*inv_squared_c
    !!     ! Va^2/c^2 / (1+Va^2/c^2)
    !!     b2(ix^D)=Vaoc/(1.d0+Vaoc)
    !!  {end do\}
    !!  ! bu . F
    !!  tmp(ixO^S)=sum(bu(ixO^S,1:ndir)*B(ixO^S,1:ndir),dim=ndim+1)
    !!  ! Rempel 2017 ApJ 834, 10 equation (54)
    !!  do idir=1,ndir
    !!    J(ixO^S,idir)=b2(ixO^S)*(B(ixO^S,idir)-bu(ixO^S,idir)*tmp(ixO^S))
    !!  end do
    !!  !! Rempel 2017 ApJ 834, 10 equation (29) add SR force at momentum equation
    !!  do idir=1,ndir
    !!    w(ixO^S,mom(idir))=w(ixO^S,mom(idir))+qdt*J(ixO^S,idir)
    !!  end do
    !!  ! Rempel 2017 ApJ 834, 10 equation (30) add work of Lorentz force and SR force
    !!  w(ixO^S,e_)=w(ixO^S,e_)+qdt*sum(wCTprim(ixO^S,mom(1:ndir))*&
    !!          (JxB(ixO^S,1:ndir)+J(ixO^S,1:ndir)),dim=ndim+1)
    !!else
      ! add work of Lorentz force
      w(ixO^S,e_)=w(ixO^S,e_)+qdt*sum(wCTprim(ixO^S,mom(1:ndir))*JxB(ixO^S,1:ndir),dim=ndim+1)
    !!end if

    if(mhd_ambipolar_sts)then
      call add_source_ambipolar_internal_energy(qdt,ixI^L,ixO^L,wCT,w,x)
    end if

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_hydrodynamic_e')

  end subroutine add_source_hydrodynamic_e

  !> Add resistive source to w within ixO Uses 3 point stencil (1 neighbour) in
  !> each direction, non-conservative. Uses the generic Laplacian
  !> with fourth order central difference (on uniform cartesian) for the laplacian. Then the
  !> stencil is 5 (2 neighbours). NOTE: Unused subroutine!
  subroutine add_source_res1(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ixA^L,idir,jdir,kdir,idirmin,idim
    double precision :: tmp(ixI^S),tmp2(ixI^S)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S)
    double precision :: gradeta(ixI^S,1:ndim), Bf(ixI^S,1:ndir)
    double precision :: lapl_vec(ixI^S,1:ndir)

    ! Calculating resistive sources involves one extra layer
    ! asking here for two, so Cartesian works with 4th order CD
    ixA^L=ixO^L^LADD2;

    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixI^L,ixO^L,idirmin,current)

    if (mhd_eta>zero)then
       eta(ixA^S)=mhd_eta
       gradeta(ixO^S,1:ndim)=zero
    else
       call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,current,eta)
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

    call laplacian_of_vector(Bf,ixI^L,ixO^L,lapl_vec)

    do idir=1,ndir
       ! Multiply by eta to store eta*Laplace B_idir
       tmp(ixO^S)=lapl_vec(ixO^S,idir)*eta(ixO^S)
       
       ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
       if (mhd_eta<zero)then
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

    if(mhd_energy) then
      ! de/dt+=eta*J**2
      w(ixO^S,e_)=w(ixO^S,e_)+qdt*eta(ixO^S)*sum(current(ixO^S,:)**2,dim=ndim+1)
    end if

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res1')

  end subroutine add_source_res1

  !> Add resistive source to w within ixO in an explicit fashion
  !> Uses 5 point stencil (2 neighbours) in each direction, conservative
  subroutine add_source_res_exp(qdt,ixI^L,ixO^L,wCT,w,x)
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
         call mpistop("Error in add_source_res_exp: Non-conforming input limits")

    ixA^L=ixO^L^LADD1;
    ! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
    ! Determine exact value of idirmin while doing the loop.
    call get_current(wCT,ixI^L,ixA^L,idirmin,current)

    tmpvec=zero
    if(mhd_eta>zero)then
      do idir=idirmin,3
        tmpvec(ixA^S,idir)=current(ixA^S,idir)*mhd_eta
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

    if(mhd_energy) then
      if(mhd_eta>zero)then
        tmp(ixO^S)=qdt*mhd_eta*sum(current(ixO^S,:)**2,dim=ndim+1)
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

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res_exp')
  end subroutine add_source_res_exp


  !> Add ambipolar source to w within ixO in an explicit fashion
  !> Uses 5 point stencil (2 neighbours) in each direction, conservative
  subroutine add_source_ambi_exp(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: current(ixI^S,1:3),curlj(ixI^S,1:3)
    double precision :: tmpvec(ixI^S,1:3),tmp(ixI^S),btot2(ixI^S)
    integer :: ixA^L,idir,idirmin1

    ixA^L=ixO^L^LADD2;

    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_ambi_exp: Non-conforming input limits")

    ixA^L=ixO^L^LADD1;
    ! Calculate -J_perpB = (JxB)xB
    call mhd_get_jxbxb(wCT,x,ixI^L,ixA^L,current)

    tmpvec=current
    do idir=1,3
      !set electric field in tmpvec : E=nuA * jxbxb, where nuA=-etaA/rho^2
      !tmpvec(ixA^S,i) = -(mhd_eta_ambi/w(ixA^S, rho_)**2) * jxbxb(ixA^S,i)
      call multiplyAmbiCoef(ixI^L,ixA^L,tmpvec(ixI^S,idir),wCT,x)
    end do

    ! dB/dt= -curl(J_perpB*etaA), thus B_i=B_i-eps_ijk d_j Jeta_k
    call curlvector(tmpvec,ixI^L,ixO^L,curlj,idirmin1,1,3)
    if(stagger_grid) then
      if(ndim==2.and.ndir==3) then
        ! if 2.5D
        w(ixO^S,mag(ndir)) = w(ixO^S,mag(ndir))-qdt*curlj(ixO^S,ndir)
      end if
    else
      w(ixO^S,mag(1:ndir)) = w(ixO^S,mag(1:ndir))-qdt*curlj(ixO^S,1:ndir)
    end if

    if(mhd_energy) then
      ! compute ambipolar heating term: nuA* J_perpB^2/ B^2
      ! avoiding nulls here
      btot2(ixA^S)=mhd_mag_en_all(wCT,ixI^L,ixA^L)
      where (btot2(ixA^S)>smalldouble )
          tmp(ixA^S) = sum(current(ixA^S,1:3)**2,dim=ndim+1) / btot2(ixA^S)
      elsewhere
          tmp(ixA^S) = zero
      endwhere
      ! multiply with nuA where nuA=-etaA/rho^2
      call multiplyAmbiCoef(ixI^L,ixA^L,tmp,wCT,x)
      ! compensate - sign and add timestep
      tmp(ixO^S)=-qdt*tmp(ixO^S)
      if(total_energy) then
        ! de/dt= +div(B x E_ambi) = eta J^2 - B dot curl(eta J)
        ! de1/dt= eta J^2 - B1 dot curl(eta J)
        w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)-&
        qdt*sum(wCT(ixO^S,mag(1:ndir))*curlj(ixO^S,1:ndir),dim=ndim+1)
      else
        ! add eta*J**2 source term in the internal or hydrodynamic energy equation
        w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)
      end if
    end if

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_ambi_exp')
  end subroutine add_source_ambi_exp

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
    ehyper(ixA^S,1:ndir) = - tmpvec(ixA^S,1:ndir)*mhd_eta_hyper

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

    if (fix_small_values)  call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_hyperres')

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
    if (mhd_glm_alpha < zero) then
      w(ixO^S,psi_) = abs(mhd_glm_alpha)*wCT(ixO^S,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab_uniform) then
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(dxlevel(:)))*w(ixO^S,psi_)
      else
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(block%ds(ixO^S,:),dim=ndim+1))*w(ixO^S,psi_)
      end if
    end if

    if(mhd_glm_extended) then
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
      call get_divb(wCT,ixI^L,ixO^L,divb, mhd_divb_nth)

      ! m = m - qdt b div b
      do idir=1,ndir
        w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-qdt*Ba(ixO^S,idir)*divb(ixO^S)
      end do
    end if

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_glm')

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
    call get_divb(wCT,ixI^L,ixO^L,divb, mhd_divb_nth)

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

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision                :: divb(ixI^S)
    integer                         :: idir, ix^D

    ! calculate div B
    call get_divb(wCT,ixI^L,ixO^L,divb, mhd_divb_nth)

   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! b = b - qdt v * div b
      ^C&w(ix^D,b^C_)=w(ix^D,b^C_)-qdt*wCT(ix^D,m^C_)*divb(ix^D)\
   {end do\}

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_janhunen')

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
    call get_divb(wCT,ixI^L,ixp^L,divb, mhd_divb_nth)

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
       call gradient(divb,ixI^L,ixp^L,idim,graddivb)

      {do i^DB=ixpmin^DB,ixpmax^DB\}
         ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
         graddivb(i^D)=graddivb(i^D)*divbdiff/(^D&1.0d0/block%ds({i^D},^D)**2+)

         w(i^D,mag(idim))=w(i^D,mag(idim))+graddivb(i^D)

         if (typedivbdiff=='all' .and. total_energy) then
           ! e += B_idim*eta*grad_idim(divb)
           w(i^D,e_)=w(i^D,e_)+wCT(i^D,mag(idim))*graddivb(i^D)
         end if
      {end do\}
    end do

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_linde')

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
    invB(ixO^S)=sqrt(mhd_mag_en_all(w,ixI^L,ixO^L))
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

  !> If resistivity is not zero, check diffusion time limit for dt and similar other effects
  subroutine mhd_get_dt(wprim,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt
    use mod_cak_force, only: cak_get_dt
    use mod_fld, only: fld_radforce_get_dt

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx^D
    double precision, intent(in)    :: wprim(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision              :: dxarr(ndim)
    double precision              :: current(ixI^S,7-2*ndir:3),eta(ixI^S)
    integer                       :: idirmin,idim

    dtnew = bigdouble

    ^D&dxarr(^D)=dx^D;
    if (mhd_eta>zero)then
      if(slab_uniform) then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/mhd_eta
      else
       dtnew=dtdiffpar*minval(block%ds(ixO^S,1:ndim))**2/mhd_eta
      end if
    else if (mhd_eta<zero)then
       call get_current(wprim,ixI^L,ixO^L,idirmin,current)
       call usr_special_resistivity(wprim,ixI^L,ixO^L,idirmin,x,current,eta)
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

    if(mhd_eta_hyper>zero) then
      if(slab_uniform) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/mhd_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixO^S,1:ndim))**4/mhd_eta_hyper,dtnew)
      end if
    end if

    if(mhd_viscosity) then
      call viscosity_get_dt(wprim,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(mhd_gravity) then
      call gravity_get_dt(wprim,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(mhd_ambipolar_exp) then
      dtnew=min(dtdiffpar*get_ambipolar_dt(wprim,ixI^L,ixO^L,dx^D,x),dtnew)
    endif

    if (mhd_cak_force) then
      call cak_get_dt(wprim,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(mhd_radiation_fld) then
      call fld_radforce_get_dt(wprim,ixI^L,ixO^L,dtnew,dx^D,x,fld_fl)
    endif

  end subroutine mhd_get_dt

  !> Wrappers for the FLD implicit (MG diffusion) hooks: phys_implicit_update /
  !> phys_evaluate_implicit have fixed interfaces with no fluid argument, so
  !> these inject the module's fld_fl object into the threaded fld routines.
  subroutine mhd_fld_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_fld, only: fld_implicit_update
    type(state), target          :: psa(max_blocks)
    type(state), target          :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    call fld_implicit_update(dtfactor,qdt,qtC,psa,psb,fld_fl)
  end subroutine mhd_fld_implicit_update

  subroutine mhd_fld_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    use mod_fld, only: fld_evaluate_implicit
    type(state), target          :: psa(max_blocks)
    double precision, intent(in) :: qtC

    call fld_evaluate_implicit(qtC,psa,fld_fl)
  end subroutine mhd_fld_evaluate_implicit

  ! Add geometrical source terms to w
  ! Geometric sources to momentum and induction
  ! for the regular case, not semi-relativistic, nor any splitting active
  ! but possibly no energy equation at all
  ! NOTE: Hall terms in induction not handled yet
  subroutine mhd_add_source_geom(qdt,dtfactor,ixI^L,ixO^L,wCT,wprim,w,x)
    use mod_global_parameters
    use mod_geometry
    use mod_rotating_frame, only: rotating_frame_add_source
    use mod_usr_methods, only: usr_set_adiab, usr_set_gamma

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor,x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw),wprim(ixI^S,1:nw),w(ixI^S,1:nw)

    double precision :: adiabs(ixI^S), gammas(ixI^S)
    double precision :: tmp,tmp1,invr,cot
    integer          :: ix^D
    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    if(.not.mhd_energy) then
      if(associated(usr_set_adiab)) then
         call usr_set_adiab(w,x,ixI^L,ixO^L,adiabs)
      else
         adiabs=mhd_adiab
      end if
      if(associated(usr_set_gamma)) then
         call usr_set_gamma(w,x,ixI^L,ixO^L,gammas)
      else
         gammas=eos%gamma
      end if
    end if

    select case (coordinate)
    case (cylindrical)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! include dt in invr, invr is always used with qdt
        if(local_timestep) then
          invr=block%dt(ix^D) * dtfactor/x(ix^D,1)
        else
          invr=qdt/x(ix^D,1)
        end if
        if(mhd_energy) then
          tmp=wprim(ix^D,p_)+half*(^C&wprim(ix^D,b^C_)**2+)
        else
          tmp=adiabs(ix^D)*wprim(ix^D,rho_)**gammas(ix^D)+half*(^C&wprim(ix^D,b^C_)**2+)
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
        if(mhd_glm) w(ix^D,br_)=w(ix^D,br_)+wprim(ix^D,psi_)*invr
     {end do\}
    case (spherical)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! include dt in invr, invr is always used with qdt
        if(local_timestep) then
          invr=block%dt(ix^D) * dtfactor/x(ix^D,1)
        else
          invr=qdt/x(ix^D,1)
        end if
        if(mhd_energy) then
          tmp1=wprim(ix^D,p_)+half*(^C&wprim(ix^D,b^C_)**2+)
        else
          tmp1=adiabs(ix^D)*wprim(ix^D,rho_)**gammas(ix^D)+half*(^C&wprim(ix^D,b^C_)**2+)
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
        if(mhd_glm) then
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
          if(mhd_glm) then
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
          if(mhd_glm) then
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

    if (mhd_rotating_frame) then
       call rotating_frame_add_source(qdt,dtfactor,ixI^L,ixO^L,wprim,w,x)
    end if

  end subroutine mhd_add_source_geom

  ! Add geometrical source terms to w
  ! Geometric sources to momentum and induction
  ! for the semi-relativistic, hence no splitting active
  ! but possibly no energy equation at all
  ! NOTE: Hall terms in induction not handled yet
  subroutine mhd_add_source_geom_semirelati(qdt,dtfactor,ixI^L,ixO^L,wCT,wprim,w,x)
    use mod_global_parameters
    use mod_geometry
    use mod_rotating_frame, only: rotating_frame_add_source
    use mod_usr_methods, only: usr_set_adiab, usr_set_gamma

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor,x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw),wprim(ixI^S,1:nw),w(ixI^S,1:nw)

    double precision :: adiabs(ixI^S), gammas(ixI^S)
    double precision :: tmp,tmp1,tmp2,invr,cot,ef(ixO^S,1:ndir)
    integer          :: ix^D
    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    if(.not.mhd_energy) then
      if(associated(usr_set_adiab)) then
         call usr_set_adiab(w,x,ixI^L,ixO^L,adiabs)
      else
         adiabs=mhd_adiab
      end if
      if(associated(usr_set_gamma)) then
         call usr_set_gamma(w,x,ixI^L,ixO^L,gammas)
      else
         gammas=eos%gamma
      end if
    end if

    select case (coordinate)
    case (cylindrical)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! include dt in invr, invr is always used with qdt
        if(local_timestep) then
          invr=block%dt(ix^D) * dtfactor/x(ix^D,1)
        else
          invr=qdt/x(ix^D,1)
        end if
        if(mhd_energy) then
          tmp=wprim(ix^D,p_)
        else
          tmp=adiabs(ix^D)*wprim(ix^D,rho_)**gammas(ix^D)
        end if
        ! E=Bxv
        {^IFTHREEC
        ef(ix^D,1)=wprim(ix^D,b2_)*wprim(ix^D,m3_)-wprim(ix^D,b3_)*wprim(ix^D,m2_)
        ef(ix^D,2)=wprim(ix^D,b3_)*wprim(ix^D,m1_)-wprim(ix^D,b1_)*wprim(ix^D,m3_)
        ef(ix^D,3)=wprim(ix^D,b1_)*wprim(ix^D,m2_)-wprim(ix^D,b2_)*wprim(ix^D,m1_)
        }
        {^IFTWOC
        ef(ix^D,1)=zero
        ! store e3 in e2 to count e3 when ^C is from 1 to 2
        ef(ix^D,2)=wprim(ix^D,b1_)*wprim(ix^D,m2_)-wprim(ix^D,b2_)*wprim(ix^D,m1_)
        }
        {^IFONEC
        ef(ix^D,1)=zero
        }
        if(phi_>0) then
          w(ix^D,mr_)=w(ix^D,mr_)+invr*(tmp+&
           half*((^C&wprim(ix^D,b^C_)**2+)+(^C&ef(ix^D,^C)**2+)*inv_squared_c) -&
                    wprim(ix^D,bphi_)**2+wprim(ix^D,rho_)*wprim(ix^D,mphi_)**2)
          w(ix^D,mphi_)=w(ix^D,mphi_)+invr*(&
                   -wprim(ix^D,rho_)*wprim(ix^D,mphi_)*wprim(ix^D,mr_) &
                   +wprim(ix^D,bphi_)*wprim(ix^D,br_)+ef(ix^D,phi_)*ef(ix^D,1)*inv_squared_c)
          if(.not.stagger_grid) then
            w(ix^D,bphi_)=w(ix^D,bphi_)+invr*&
                     (wprim(ix^D,bphi_)*wprim(ix^D,mr_) &
                     -wprim(ix^D,br_)*wprim(ix^D,mphi_))
          end if
        else
          w(ix^D,mr_)=w(ix^D,mr_)+invr*(tmp+half*((^C&wprim(ix^D,b^C_)**2+)+&
             (^C&ef(ix^D,^C)**2+)*inv_squared_c))
        end if
        if(mhd_glm) w(ix^D,br_)=w(ix^D,br_)+wprim(ix^D,psi_)*invr
     {end do\}
    case (spherical)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! include dt in invr, invr is always used with qdt
        if(local_timestep) then
          invr=block%dt(ix^D)*dtfactor/x(ix^D,1)
        else
          invr=qdt/x(ix^D,1)
        end if
        ! E=Bxv
        {^IFTHREEC
        ef(ix^D,1)=wprim(ix^D,b2_)*wprim(ix^D,m3_)-wprim(ix^D,b3_)*wprim(ix^D,m2_)
        ef(ix^D,2)=wprim(ix^D,b3_)*wprim(ix^D,m1_)-wprim(ix^D,b1_)*wprim(ix^D,m3_)
        ef(ix^D,3)=wprim(ix^D,b1_)*wprim(ix^D,m2_)-wprim(ix^D,b2_)*wprim(ix^D,m1_)
        }
        {^IFTWOC
        ! store e3 in e1 to count e3 when ^C is from 1 to 2
        ef(ix^D,1)=wprim(ix^D,b1_)*wprim(ix^D,m2_)-wprim(ix^D,b2_)*wprim(ix^D,m1_)
        ef(ix^D,2)=zero
        }
        {^IFONEC
        ef(ix^D,1)=zero
        }
        if(mhd_energy) then
          tmp1=wprim(ix^D,p_)+half*((^C&wprim(ix^D,b^C_)**2+)+(^C&ef(ix^D,^C)**2+)*inv_squared_c)
        else
          tmp1=adiabs(ix^D)*wprim(ix^D,rho_)**gammas(ix^D)+half*((^C&wprim(ix^D,b^C_)**2+)+(^C&ef(ix^D,^C)**2+)*inv_squared_c)
        end if
        ! m1
        {^IFONEC
        w(ix^D,m1_)=w(ix^D,m1_)+two*tmp1*invr
        }
        {^NOONEC
        w(ix^D,m1_)=w(ix^D,m1_)+invr*&
           (two*tmp1+(^CE&wprim(ix^D,rho_)*wprim(ix^D,m^CE_)**2-&
            wprim(ix^D,b^CE_)**2-ef(ix^D,^CE)**2*inv_squared_c+))
        }
        ! b1
        if(mhd_glm) then
          w(ix^D,b1_)=w(ix^D,b1_)+invr*2.0d0*wprim(ix^D,psi_)
        end if
        {^IFONED
        cot=0.d0
        }
        {^NOONED
        cot=1.d0/tan(x(ix^D,2))
        }
        {^IFTWOC
        ! m2
        w(ix^D,m2_)=w(ix^D,m2_)+invr*(tmp1*cot-wprim(ix^D,rho_)*wprim(ix^D,m1_)*wprim(ix^D,m2_)&
            +wprim(ix^D,b1_)*wprim(ix^D,b2_)+ef(ix^D,1)*ef(ix^D,2)*inv_squared_c)
        ! b2
        if(.not.stagger_grid) then
          tmp=wprim(ix^D,m1_)*wprim(ix^D,b2_)-wprim(ix^D,m2_)*wprim(ix^D,b1_)
          if(mhd_glm) then
            tmp=tmp+wprim(ix^D,psi_)*cot
          end if
          w(ix^D,b2_)=w(ix^D,b2_)+tmp*invr
        end if
        }

        {^IFTHREEC
        ! m2
        w(ix^D,m2_)=w(ix^D,m2_)+invr*(tmp1*cot-wprim(ix^D,rho_)*wprim(ix^D,m1_)*wprim(ix^D,m2_) &
            +wprim(ix^D,b1_)*wprim(ix^D,b2_)+ef(ix^D,1)*ef(ix^D,2)*inv_squared_c&
            +(wprim(ix^D,rho_)*wprim(ix^D,m3_)**2&
            -wprim(ix^D,b3_)**2-ef(ix^D,3)**2*inv_squared_c)*cot)
        ! b2
        if(.not.stagger_grid) then
          tmp=wprim(ix^D,m1_)*wprim(ix^D,b2_)-wprim(ix^D,m2_)*wprim(ix^D,b1_)
          if(mhd_glm) then
            tmp=tmp+wprim(ix^D,psi_)*cot
          end if
          w(ix^D,b2_)=w(ix^D,b2_)+tmp*invr
        end if
        ! m3
        w(ix^D,m3_)=w(ix^D,m3_)+invr*&
            (-wprim(ix^D,m3_)*wprim(ix^D,m1_)*wprim(ix^D,rho_) &
             +wprim(ix^D,b3_)*wprim(ix^D,b1_) &
             +ef(ix^D,3)*ef(ix^D,1)*inv_squared_c&
           +(-wprim(ix^D,m2_)*wprim(ix^D,m3_)*wprim(ix^D,rho_) &
             +wprim(ix^D,b2_)*wprim(ix^D,b3_)&
             +ef(ix^D,2)*ef(ix^D,3)*inv_squared_c)*cot)
        ! b3
        if(.not.stagger_grid) then
          w(ix^D,b3_)=w(ix^D,b3_)+invr*&
             (wprim(ix^D,m1_)*wprim(ix^D,b3_) &
             -wprim(ix^D,m3_)*wprim(ix^D,b1_) &
            -(wprim(ix^D,m3_)*wprim(ix^D,b2_) &
             -wprim(ix^D,m2_)*wprim(ix^D,b3_))*cot)
        end if
        }
     {end do\}
    end select

    if (mhd_rotating_frame) then
       call rotating_frame_add_source(qdt,dtfactor,ixI^L,ixO^L,wprim,w,x)
    end if

  end subroutine mhd_add_source_geom_semirelati

  ! Add geometrical source terms to w
  ! Geometric sources to momentum and induction
  ! for those cases where any kind of splitting (B0field or has_equi_rho_and_p) is active
  ! This implies that there is an energy equation included for sure
  ! B0field impacts terms in induction equation and geometric sources for them
  ! both flags affect the terms in momentum equation, in three variants (TF, TT, FT)
  ! NOTE: Hall terms in induction not handled yet
  subroutine mhd_add_source_geom_split(qdt,dtfactor,ixI^L,ixO^L,wCT,wprim,w,x)
    use mod_global_parameters
    use mod_geometry
    use mod_rotating_frame, only: rotating_frame_add_source
    use mod_usr_methods, only: usr_set_adiab, usr_set_gamma

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor,x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw),wprim(ixI^S,1:nw),w(ixI^S,1:nw)

    double precision :: tmp,tmp1,tmp2,invr,cot
    double precision :: adiabs(ixI^S), gammas(ixI^S)
    integer          :: ix^D
    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    if(.not.mhd_energy) then
      if(associated(usr_set_adiab)) then
        call usr_set_adiab(wprim,x,ixI^L,ixO^L,adiabs)
      else
        adiabs=mhd_adiab
      end if
      if(associated(usr_set_gamma)) then
        call usr_set_gamma(wprim,x,ixI^L,ixO^L,gammas)
      else
        gammas=eos%gamma
      end if
    end if

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
        if(mhd_energy) then
          tmp=wprim(ix^D,p_)+half*(^C&wprim(ix^D,b^C_)**2+)
        else
          tmp=adiabs(ix^D)*wprim(ix^D,rho_)**gammas(ix^D)+half*(^C&wprim(ix^D,b^C_)**2+)
        end if
        if(B0field) tmp=tmp+(^C&block%B0(ix^D,^C,0)*wprim(ix^D,b^C_)+)
        if(phi_>0) then
          w(ix^D,mr_)=w(ix^D,mr_)+invr*(tmp-&
                    wprim(ix^D,bphi_)**2+wprim(ix^D,mphi_)*wCT(ix^D,mphi_))
          if(B0field) then
            w(ix^D,mr_)=w(ix^D,mr_)+invr*(-block%B0(ix^D,phi_,0)*wprim(ix^D,bphi_)-wprim(ix^D,bphi_)*block%B0(ix^D,phi_,0))
          endif
          w(ix^D,mphi_)=w(ix^D,mphi_)+invr*(&
                   -wCT(ix^D,mphi_)*wprim(ix^D,mr_) &
                   +wprim(ix^D,bphi_)*wprim(ix^D,br_))
          if(B0field) then
            w(ix^D,mphi_)=w(ix^D,mphi_)+invr*(block%B0(ix^D,phi_,0)*wprim(ix^D,br_)+wprim(ix^D,bphi_)*block%B0(ix^D,r_,0))
          endif
          if(.not.stagger_grid) then
            w(ix^D,bphi_)=w(ix^D,bphi_)+invr*&
                     (wprim(ix^D,bphi_)*wprim(ix^D,mr_) &
                     -wprim(ix^D,br_)*wprim(ix^D,mphi_))
            if(B0field) then
              w(ix^D,bphi_)=w(ix^D,bphi_)+invr*&
                     (block%B0(ix^D,phi_,0)*wprim(ix^D,mr_) &
                     -block%B0(ix^D,r_,0)*wprim(ix^D,mphi_))
            endif
          end if
        else
          w(ix^D,mr_)=w(ix^D,mr_)+invr*tmp
        end if
        if(mhd_glm) w(ix^D,br_)=w(ix^D,br_)+wprim(ix^D,psi_)*invr
     {end do\}
    case (spherical)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! include dt in invr, invr is always used with qdt
        if(local_timestep) then
          invr=block%dt(ix^D) * dtfactor/x(ix^D,1)
        else
          invr=qdt/x(ix^D,1)
        end if
        tmp1=wprim(ix^D,p_)+half*(^C&wprim(ix^D,b^C_)**2+)
        if(B0field) tmp2=(^C&block%B0(ix^D,^C,0)*wprim(ix^D,b^C_)+)
        ! m1
        {^IFONEC
        w(ix^D,mom(1))=w(ix^D,mom(1))+two*tmp1*invr
        if(B0field) w(ix^D,mom(1))=w(ix^D,mom(1))+two*tmp2*invr
        }
        {^NOONEC
        if(B0field) then
          w(ix^D,mom(1))=w(ix^D,mom(1))+invr*&
           (two*(tmp1+tmp2)+(^CE&wprim(ix^D,m^CE_)*wCT(ix^D,m^CE_)-wprim(ix^D,b^CE_)**2+)- &
            (^CE&two*block%B0(ix^D,^CE,0)*wprim(ix^D,b^CE_)+))
        else
          w(ix^D,mom(1))=w(ix^D,mom(1))+invr*&
           (two*tmp1+(^CE&wprim(ix^D,m^CE_)*wCT(ix^D,m^CE_)-wprim(ix^D,b^CE_)**2+))
        end if
        }
        ! b1
        if(mhd_glm) then
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
        if(B0field) then
          w(ix^D,mom(2))=w(ix^D,mom(2))+invr*((tmp1+tmp2)*cot-wprim(ix^D,m1_)*wCT(ix^D,m2_)&
            +wprim(ix^D,b1_)*wprim(ix^D,b2_)+block%B0(ix^D,1,0)*wprim(ix^D,b2_)&
            +wprim(ix^D,b1_)*block%B0(ix^D,2,0))
        else
          w(ix^D,mom(2))=w(ix^D,mom(2))+invr*(tmp1*cot-wprim(ix^D,m1_)*wCT(ix^D,m2_)&
            +wprim(ix^D,b1_)*wprim(ix^D,b2_))
        end if
        ! b2
        if(.not.stagger_grid) then
          if(B0field) then
            tmp=wprim(ix^D,m1_)*wprim(ix^D,b2_)-wprim(ix^D,m2_)*wprim(ix^D,b1_)&
             +wprim(ix^D,m1_)*block%B0(ix^D,2,0)-wprim(ix^D,m2_)*block%B0(ix^D,1,0)
          else
            tmp=wprim(ix^D,m1_)*wprim(ix^D,b2_)-wprim(ix^D,m2_)*wprim(ix^D,b1_)
          end if
          if(mhd_glm) then
            tmp=tmp+wprim(ix^D,psi_)*cot
          end if
          w(ix^D,mag(2))=w(ix^D,mag(2))+tmp*invr
        end if
        }
        {^IFTHREEC
        ! m2
        if(B0field) then
          w(ix^D,mom(2))=w(ix^D,mom(2))+invr*((tmp1+tmp2)*cot-wprim(ix^D,m1_)*wCT(ix^D,m2_)&
            +wprim(ix^D,b1_)*wprim(ix^D,b2_)+block%B0(ix^D,1,0)*wprim(ix^D,b2_)&
            +wprim(ix^D,b1_)*block%B0(ix^D,2,0)&
            +(wprim(ix^D,m3_)*wCT(ix^D,m3_)-wprim(ix^D,b3_)**2-two*block%B0(ix^D,3,0)*wprim(ix^D,b3_))*cot)
        else
          w(ix^D,mom(2))=w(ix^D,mom(2))+invr*(tmp1*cot-wprim(ix^D,m1_)*wCT(ix^D,m2_)&
            +wprim(ix^D,b1_)*wprim(ix^D,b2_)&
            +(wprim(ix^D,m3_)*wCT(ix^D,m3_)-wprim(ix^D,b3_)**2)*cot)
        end if
        ! b2
        if(.not.stagger_grid) then
          if(B0field) then
            tmp=wprim(ix^D,m1_)*wprim(ix^D,b2_)-wprim(ix^D,m2_)*wprim(ix^D,b1_)&
             +wprim(ix^D,m1_)*block%B0(ix^D,2,0)-wprim(ix^D,m2_)*block%B0(ix^D,1,0)
          else
            tmp=wprim(ix^D,m1_)*wprim(ix^D,b2_)-wprim(ix^D,m2_)*wprim(ix^D,b1_)
          end if
          if(mhd_glm) then
            tmp=tmp+wprim(ix^D,psi_)*cot
          end if
          w(ix^D,mag(2))=w(ix^D,mag(2))+tmp*invr
        end if
        ! m3
        if(B0field) then
          w(ix^D,mom(3))=w(ix^D,mom(3))-invr*&
               (wprim(ix^D,m3_)*wCT(ix^D,m1_) &
               -wprim(ix^D,b3_)*wprim(ix^D,b1_) &
            +block%B0(ix^D,1,0)*wprim(ix^D,b3_) &
            +wprim(ix^D,b1_)*block%B0(ix^D,3,0) &
              +(wprim(ix^D,m2_)*wCT(ix^D,m3_) &
               -wprim(ix^D,b2_)*wprim(ix^D,b3_) &
            +block%B0(ix^D,2,0)*wprim(ix^D,b3_) &
            +wprim(ix^D,b2_)*block%B0(ix^D,3,0))*cot)
        else
          w(ix^D,mom(3))=w(ix^D,mom(3))-invr*&
               (wprim(ix^D,m3_)*wCT(ix^D,m1_) &
               -wprim(ix^D,b3_)*wprim(ix^D,b1_) &
              +(wprim(ix^D,m2_)*wCT(ix^D,m3_) &
               -wprim(ix^D,b2_)*wprim(ix^D,b3_))*cot)
        end if
        ! b3
        if(.not.stagger_grid) then
          if(B0field) then
            w(ix^D,mag(3))=w(ix^D,mag(3))+invr*&
               (wprim(ix^D,m1_)*wprim(ix^D,b3_) &
               -wprim(ix^D,m3_)*wprim(ix^D,b1_) &
            +wprim(ix^D,m1_)*block%B0(ix^D,3,0) &
            -wprim(ix^D,m3_)*block%B0(ix^D,1,0) &
              -(wprim(ix^D,m3_)*wprim(ix^D,b2_) &
               -wprim(ix^D,m2_)*wprim(ix^D,b3_) &
            +wprim(ix^D,m3_)*block%B0(ix^D,2,0) &
            -wprim(ix^D,m2_)*block%B0(ix^D,3,0))*cot)
          else
            w(ix^D,mag(3))=w(ix^D,mag(3))+invr*&
               (wprim(ix^D,m1_)*wprim(ix^D,b3_) &
               -wprim(ix^D,m3_)*wprim(ix^D,b1_) &
              -(wprim(ix^D,m3_)*wprim(ix^D,b2_) &
               -wprim(ix^D,m2_)*wprim(ix^D,b3_))*cot)
          end if
        end if
        }
     {end do\}
    end select

    if (mhd_rotating_frame) then
       call rotating_frame_add_source(qdt,dtfactor,ixI^L,ixO^L,wprim,w,x)
    end if

  end subroutine mhd_add_source_geom_split

  !> Compute 2 times total magnetic energy
  function mhd_mag_en_all(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    if (B0field) then
      mge = sum((w(ixO^S, mag(:))+block%B0(ixO^S,:,b0i))**2, dim=ndim+1)
    else
      mge = sum(w(ixO^S, mag(:))**2, dim=ndim+1)
    end if
  end function mhd_mag_en_all

  subroutine mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall,partial)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: vHall(ixI^S,1:ndir)
    logical, intent(in), optional :: partial

    double precision :: current(ixI^S,7-2*ndir:3)
    double precision :: rho(ixI^S)
    integer          :: idir, idirmin, ix^D
    logical          :: use_partial

    use_partial=.false.
    if(present(partial)) use_partial=partial
    call mhd_get_rho(w,x,ixI^L,ixO^L,rho)
    if(.not.use_partial)then
       ! Calculate current density and idirmin, including J0 when split
       call get_current(w,ixI^L,ixO^L,idirmin,current)
    else
       if(slab_uniform) then
         ! fourth order CD in cartesian
         call curlvector(w(ixI^S,mag(1:ndir)),ixI^L,ixO^L,current,idirmin,7-2*ndir,ndir,.true.)
       else
         call curlvector(w(ixI^S,mag(1:ndir)),ixI^L,ixO^L,current,idirmin,7-2*ndir,ndir)
       endif
    endif
    do idir = idirmin, ndir
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         vHall(ix^D,idir)=-mhd_etah*current(ix^D,idir)/rho(ix^D)
      {end do\}
    end do

  end subroutine mhd_getv_Hall

  subroutine mhd_modify_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)
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

  end subroutine mhd_modify_wLR

  subroutine mhd_boundary_adjust(igrid,psb)
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

  end subroutine mhd_boundary_adjust

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
      !  if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
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
      !  if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
     case(2)
      !  if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
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
      !  if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
     case(3)
      !  if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
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
      !  if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
     case(4)
      !  if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
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
      !  if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
     {^IFTHREED
     case(5)
      !  if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
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
      !  if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
     case(6)
      !  if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_primitive(ixG^L,ixO^L,w,x)
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
      !  if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
       if(total_energy) call eos%to_conserved(ixG^L,ixO^L,w,x)
     }
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine fixdivB_boundary

  {^NOONED
  subroutine mhd_clean_divb_multigrid(qdt, qt, active)
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
          write(*,*) "mhd_clean_divb_multigrid warning: unknown boundary type"
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
            mhd_divb_nth)
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
         call mhd_face_to_center(ixM^LL,ps(igrid))
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

  end subroutine mhd_clean_divb_multigrid
  }

  !> get electric field through averaging neighors to update faces in CT
  subroutine mhd_update_faces_average(ixI^L,ixO^L,qt,qdt,wp,fC,fE,sCT,s,vcts)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wp(ixI^S,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,sdim:3)

    double precision                   :: circ(ixI^S,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixI^S,sdim:3) :: E_resi, E_ambi
    integer                            :: ix^D,ixC^L,ixA^L,i1kr^D,i2kr^D
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! if there is resistivity, get eta J
    if(mhd_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,wp,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(mhd_ambipolar_exp) call get_ambipolar_electric_field(ixI^L,ixO^L,sCT%w,x,E_ambi)

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
              if(mhd_eta/=zero) fE(ix^D,idir)=fE(ix^D,idir)+E_resi(ix^D,idir)
              ! add ambipolar electric field
              if(mhd_ambipolar_exp) fE(ix^D,idir)=fE(ix^D,idir)+E_ambi(ix^D,idir)

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
     {do ix^DB=ixCmin^DB,ixCmax^DB\}
        ! Divide by the area of the face to get dB/dt
        if(s%surfaceC(ix^D,idim1) > smalldouble) then
          ! Time update cell-face magnetic field component
          bfaces(ix^D,idim1)=bfaces(ix^D,idim1)-circ(ix^D,idim1)/s%surfaceC(ix^D,idim1)
        end if
     {end do\}
    end do

    end associate

  end subroutine mhd_update_faces_average

  !> update faces using UCT contact mode by Gardiner and Stone 2005 JCP 205, 509
  subroutine mhd_update_faces_contact(ixI^L,ixO^L,qt,qdt,wp,fC,fE,sCT,s,vcts)
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
    double precision, dimension(ixI^S,sdim:3) :: E_resi, E_ambi
    ! current on cell edges
    double precision :: jce(ixI^S,sdim:3)
    ! location at cell faces
    double precision :: xs(ixGs^T,1:ndim)
    double precision :: gradi(ixGs^T)
    integer :: ixC^L,ixA^L
    integer :: idim1,idim2,idir,iwdim1,iwdim2,ix^D,i1kr^D,i2kr^D

    associate(bfaces=>s%ws,x=>s%x,w=>s%w,vnorm=>vcts%vnorm,wCTs=>sCT%ws)

    ! if there is resistivity, get eta J
    if(mhd_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,wp,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(mhd_ambipolar_exp) call get_ambipolar_electric_field(ixI^L,ixO^L,sCT%w,x,E_ambi)

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
              if(numerical_resistive_heating) Ein(ix^D,idir)=fE(ix^D,idir)
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
              if(numerical_resistive_heating) Ein(ix^D,idir)=fE(ix^D,idir)-Ein(ix^D,idir)
              ! add resistive electric field at cell edges E=-vxB+eta J
              if(mhd_eta/=zero) fE(ix^D,idir)=fE(ix^D,idir)+E_resi(ix^D,idir)
              ! add ambipolar electric field
              if(mhd_ambipolar_exp) fE(ix^D,idir)=fE(ix^D,idir)+E_ambi(ix^D,idir)

              ! times time step and edge length
              fE(ix^D,idir)=fE(ix^D,idir)*qdt*s%dsC(ix^D,idir)
           {end do\}
          end if
        end do
      end do
    end do

    if(numerical_resistive_heating) then
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
       !! if(nwextra>0) then
       !!   block%w(ixO^S,nw)=block%w(ixO^S,nw)+jce(ixO^S,idir)
       !! end if
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
     {do ix^DB=ixCmin^DB,ixCmax^DB\}
        ! Divide by the area of the face to get dB/dt
        if(s%surfaceC(ix^D,idim1) > smalldouble) then
          ! Time update cell-face magnetic field component
          bfaces(ix^D,idim1)=bfaces(ix^D,idim1)-circ(ix^D,idim1)/s%surfaceC(ix^D,idim1)
        end if
     {end do\}
    end do

    end associate

  end subroutine mhd_update_faces_contact

  !> update faces
  subroutine mhd_update_faces_hll(ixI^L,ixO^L,qt,qdt,wp,fC,fE,sCT,s,vcts)
    use mod_global_parameters
    use mod_usr_methods
    use mod_constrained_transport

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wp(ixI^S,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,sdim:3)

    double precision                   :: vtilL(ixI^S,2)
    double precision                   :: vtilR(ixI^S,2)
    double precision                   :: bfacetot(ixI^S,ndim)
    double precision                   :: btilL(ixI^S,ndim)
    double precision                   :: btilR(ixI^S,ndim)
    double precision                   :: cp(ixI^S,2)
    double precision                   :: cm(ixI^S,2)
    double precision                   :: circ(ixI^S,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixI^S,sdim:3) :: E_resi, E_ambi
    integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,ixCm^L
    integer                            :: idim1,idim2,idir,ix^D

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
    if(mhd_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,wp,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(mhd_ambipolar_exp) call get_ambipolar_electric_field(ixI^L,ixO^L,sCT%w,x,E_ambi)

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
      if(mhd_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+E_resi(ixC^S,idir)
      ! add ambipolar electric field
      if(mhd_ambipolar_exp) fE(ixC^S,idir)=fE(ixC^S,idir)+E_ambi(ixC^S,idir)

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
     {do ix^DB=ixCmin^DB,ixCmax^DB\}
        ! Divide by the area of the face to get dB/dt
        if(s%surfaceC(ix^D,idim1) > smalldouble) then
          ! Time update cell-face magnetic field component
          bfaces(ix^D,idim1)=bfaces(ix^D,idim1)-circ(ix^D,idim1)/s%surfaceC(ix^D,idim1)
        end if
     {end do\}
    end do

    end associate
  end subroutine mhd_update_faces_hll

  !> calculate eta J at cell edges
  subroutine get_resistive_electric_field(ixI^L,ixO^L,wp,sCT,s,jce)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)                :: ixI^L, ixO^L
    ! cell-center primitive variables
    double precision, intent(in)       :: wp(ixI^S,1:nw)
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
    if(mhd_eta>zero)then
      jce(ixI^S,:)=jce(ixI^S,:)*mhd_eta
    else
      ixA^L=ixO^L^LADD1;
      call get_current(wCT,ixI^L,ixA^L,idirmin,jcc)
      call usr_special_resistivity(wp,ixI^L,ixA^L,idirmin,x,jcc,eta)
      ! calculate eta on cell edges
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

  !> get ambipolar electric field on cell edges
  subroutine get_ambipolar_electric_field(ixI^L,ixO^L,w,x,fE)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: fE(ixI^S,sdim:3)

    double precision :: jxbxb(ixI^S,1:3)
    integer :: idir,ixA^L,ixC^L,ix^D

    ixA^L=ixO^L^LADD1;
    call mhd_get_jxbxb(w,x,ixI^L,ixA^L,jxbxb)
    ! calculate electric field on cell edges from cell centers
    do idir=sdim,3
      ! set ambipolar electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
      ! E_ambi(ixA^S,i) = -(mhd_eta_ambi/w(ixA^S, rho_)**2) * jxbxb(ixA^S,i)
      call multiplyAmbiCoef(ixI^L,ixA^L,jxbxb(ixI^S,idir),w,x)
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D+kr(idir,^D)-1;
      fE(ixC^S,idir)=0.d0
     {do ix^DB=0,1\}
        if({ ix^D==1 .and. ^D==idir | .or.}) cycle
        ixAmin^D=ixCmin^D+ix^D;
        ixAmax^D=ixCmax^D+ix^D;
        fE(ixC^S,idir)=fE(ixC^S,idir)+jxbxb(ixA^S,idir)
     {end do\}
      fE(ixC^S,idir)=fE(ixC^S,idir)*0.25d0
    end do

  end subroutine get_ambipolar_electric_field

  !> calculate cell-center values from face-center values
  subroutine mhd_face_to_center(ixO^L,s)
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

  end subroutine mhd_face_to_center

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

end module mod_mhd_phys
