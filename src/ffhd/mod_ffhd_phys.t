!> Frozen-field hydrodynamics module
module mod_ffhd_phys

#include "amrvac.h"

  use mod_global_parameters, only: std_len, const_c
  use mod_thermal_conduction, only: tc_fluid
  use mod_radiative_cooling, only: rc_fluid
  use mod_thermal_emission, only: te_fluid
  use mod_physics
  use mod_comm_lib, only: mpistop

  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: ffhd_energy = .true.

  !> Whether thermal conduction is used
  logical, public, protected              :: ffhd_thermal_conduction = .false.
  !> Whether hyperbolic type thermal conduction is used
  logical, public, protected              :: ffhd_hyperbolic_thermal_conduction = .false.
  !> type of fluid for thermal conduction
  type(tc_fluid), public, allocatable     :: tc_fl
  !> type of fluid for thermal emission synthesis
  type(te_fluid), public, allocatable     :: te_fl_ffhd

  !> Whether radiative cooling is added
  logical, public, protected              :: ffhd_radiative_cooling = .false.
  !> type of fluid for radiative cooling
  type(rc_fluid), public, allocatable     :: rc_fl

  !> Whether viscosity is added
  logical, public, protected              :: ffhd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: ffhd_gravity = .false.

  !> Whether TRAC method is used
  logical, public, protected              :: ffhd_trac = .false.

  !> Which TRAC method is used
  integer, public, protected              :: ffhd_trac_type=1

  !> Height of the mask used in the TRAC method
  double precision, public, protected     :: ffhd_trac_mask = 0.d0

  !> Distance between two adjacent traced magnetic field lines (in finest cell size)
  integer, public, protected              :: ffhd_trac_finegrid=4

  !> Whether plasma is partially ionized
  logical, public, protected              :: ffhd_partial_ionization = .false.

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> Indices of temperature
  integer, public, protected :: Te_

  !> Index of the cutoff temperature for the TRAC method
  integer, public, protected              :: Tcoff_
  integer, public, protected              :: Tweight_
  integer, public, protected              :: q_

  !> The adiabatic index
  double precision, public                :: ffhd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: ffhd_adiab = 1.0d0

  !> The small_est allowed energy
  double precision, protected             :: small_e

  !> The thermal conductivity kappa in hyperbolic thermal conduction
  double precision, public                :: hypertc_kappa

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
  ! remove the below flag  and assume default value = .false.
  ! when eq state properly implemented everywhere
  ! and not anymore through units
  logical, public, protected :: eq_state_units = .true.

  !> gamma minus one and its inverse
  double precision :: gamma_1, inv_gamma_1

  !define the subroutine interface for the ambipolar mask
  abstract interface

    function fun_kin_en(w, ixI^L, ixO^L, inv_rho) result(ke)
      use mod_global_parameters, only: nw, ndim,block
      integer, intent(in)           :: ixI^L, ixO^L
      double precision, intent(in)  :: w(ixI^S, nw)
      double precision              :: ke(ixO^S)
      double precision, intent(in), optional :: inv_rho(ixO^S)
    end function fun_kin_en

  end interface

  procedure(sub_convert), pointer      :: ffhd_to_primitive        => null()
  procedure(sub_convert), pointer      :: ffhd_to_conserved        => null()
  procedure(sub_small_values), pointer :: ffhd_handle_small_values => null()
  procedure(sub_get_pthermal), pointer :: ffhd_get_pthermal        => null()
  procedure(sub_get_pthermal), pointer :: ffhd_get_Rfactor         => null()
  procedure(sub_get_pthermal), pointer :: ffhd_get_temperature     => null()
  procedure(sub_get_v), pointer        :: ffhd_get_v               => null()
  procedure(fun_kin_en), pointer       :: ffhd_kin_en              => null()
  ! Public methods
  public :: ffhd_phys_init
  public :: ffhd_kin_en
  public :: ffhd_get_pthermal
  public :: ffhd_get_temperature
  public :: ffhd_get_v
  public :: ffhd_get_rho
  public :: ffhd_get_v_idim
  public :: ffhd_to_conserved
  public :: ffhd_to_primitive
  public :: ffhd_get_csound2
  public :: ffhd_e_to_ei
  public :: ffhd_ei_to_e

contains

  subroutine ffhd_read_params(files)
    use mod_global_parameters
    use mod_particles, only: particles_eta, particles_etah
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /ffhd_list/ ffhd_energy, ffhd_gamma, ffhd_adiab,&
      ffhd_thermal_conduction, ffhd_hyperbolic_thermal_conduction, ffhd_radiative_cooling, ffhd_gravity,&
      ffhd_viscosity, He_abundance, H_ion_fr, He_ion_fr, He_ion_fr2, eq_state_units, SI_unit,&
      B0field, Busr, ffhd_partial_ionization, ffhd_trac, ffhd_trac_type, ffhd_trac_mask, ffhd_trac_finegrid

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, ffhd_list, end=111)
111    close(unitpar)
    end do
  end subroutine ffhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine ffhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = ffhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine ffhd_write_info

  subroutine ffhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_supertimestepping, only: sts_init, add_sts_method,&
            set_conversion_methods_to_head, set_error_handling_to_head
    use mod_ionization_degree
    use mod_usr_methods, only: usr_Rfactor
    integer :: itr, idir

    call ffhd_read_params(par_files)

    if(.not. ffhd_energy) then
      if(ffhd_thermal_conduction) then
        ffhd_thermal_conduction=.false.
        if(mype==0) write(*,*) 'WARNING: set ffhd_thermal_conduction=F when ffhd_energy=F'
      end if
      if(ffhd_thermal_conduction) then
        ffhd_hyperbolic_thermal_conduction=.false.
        if(mype==0) write(*,*) 'WARNING: set ffhd_hyperbolic_thermal_conduction=F when ffhd_energy=F'
      end if
      if(ffhd_radiative_cooling) then
        ffhd_radiative_cooling=.false.
        if(mype==0) write(*,*) 'WARNING: set ffhd_radiative_cooling=F when ffhd_energy=F'
      end if
      if(ffhd_trac) then
        ffhd_trac=.false.
        if(mype==0) write(*,*) 'WARNING: set ffhd_trac=F when ffhd_energy=F'
      end if
      if(ffhd_partial_ionization) then
        ffhd_partial_ionization=.false.
        if(mype==0) write(*,*) 'WARNING: set ffhd_partial_ionization=F when ffhd_energy=F'
      end if
    end if
    if(.not.eq_state_units) then
      if(ffhd_partial_ionization) then
        ffhd_partial_ionization=.false.
        if(mype==0) write(*,*) 'WARNING: set ffhd_partial_ionization=F when eq_state_units=F'
      end if
    end if

    if(ffhd_hyperbolic_thermal_conduction) then
      ffhd_thermal_conduction=.false.
      if(mype==0) write(*,*) 'WARNING: turn off parabolic TC when using hyperbolic TC'
    end if

    physics_type = "ffhd"
    phys_energy=ffhd_energy
    phys_internal_e=.false.
    phys_trac=ffhd_trac
    phys_trac_type=ffhd_trac_type
    phys_partial_ionization=ffhd_partial_ionization

    phys_gamma = ffhd_gamma
    phys_total_energy=ffhd_energy
    phys_trac_finegrid=ffhd_trac_finegrid

    {^IFONED
    if(ffhd_trac .and. ffhd_trac_type .gt. 2) then
      ffhd_trac_type=1
      if(mype==0) write(*,*) 'WARNING: reset ffhd_trac_type=1 for 1D simulation'
    end if
    }
    if(ffhd_trac .and. ffhd_trac_type .le. 4) then
      ffhd_trac_mask=bigdouble
      if(mype==0) write(*,*) 'WARNING: set ffhd_trac_mask==bigdouble for global TRAC method'
    end if
    phys_trac_mask=ffhd_trac_mask

    allocate(start_indices(number_species),stop_indices(number_species))
    ! set the index of the first flux variable for species 1
    start_indices(1)=1
    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(1))
    mom(:) = var_set_momentum(1)

    ! Set index of energy variable
    if(ffhd_energy) then
      e_     = var_set_energy() ! energy density
      p_     = e_               ! gas pressure
    else
      e_     = -1
      p_     = -1
    end if

    if(ffhd_hyperbolic_thermal_conduction) then
      q_ = var_set_q()
      need_global_cs2max=.true.
    else
      q_=-1
    end if

    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    !  Number of variables need reconstruction in w
    nw_recon=nwflux

    !  set temperature as an auxiliary variable to get ionization degree
    if(ffhd_partial_ionization) then
      Te_ = var_set_auxvar('Te','Te')
    else
      Te_ = -1
    end if

    ! set cutoff temperature when using the TRAC method, as well as an auxiliary weight
    Tweight_ = -1
    if(ffhd_trac) then
      Tcoff_ = var_set_wextra()
      iw_Tcoff=Tcoff_
      if(ffhd_trac_type .ge. 3) then
        Tweight_ = var_set_wextra()
        iw_Tweight=Tweight_
      end if
    else
      Tcoff_ = -1
    end if

    nvector      = 0 ! No. vector vars

    ! Check whether custom flux types have been defined
    if(.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nwflux))
       flux_type = flux_default
    else if(any(shape(flux_type) /= [ndir, nwflux])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    phys_get_dt              => ffhd_get_dt
    phys_get_cmax            => ffhd_get_cmax_origin
    phys_get_a2max           => ffhd_get_a2max
    phys_get_cs2max          => ffhd_get_cs2max
    phys_get_tcutoff         => ffhd_get_tcutoff
    phys_get_cbounds         => ffhd_get_cbounds
    phys_to_primitive        => ffhd_to_primitive_origin
    ffhd_to_primitive        => ffhd_to_primitive_origin
    phys_to_conserved        => ffhd_to_conserved_origin
    ffhd_to_conserved        => ffhd_to_conserved_origin
    phys_get_flux            => ffhd_get_flux
    phys_get_v               => ffhd_get_v_origin
    ffhd_get_v               => ffhd_get_v_origin
    phys_get_rho             => ffhd_get_rho
    ffhd_kin_en              => ffhd_kin_en_origin
    !> need to check source geom here
    !phys_add_source_geom     => ffhd_add_source_geom
    phys_add_source          => ffhd_add_source
    phys_check_params        => ffhd_check_params
    phys_write_info          => ffhd_write_info
    phys_handle_small_values => ffhd_handle_small_values_origin
    ffhd_handle_small_values => ffhd_handle_small_values_origin
    phys_check_w             => ffhd_check_w_origin
 
    if(.not.ffhd_energy) then
      phys_get_pthermal      => ffhd_get_pthermal_iso
      ffhd_get_pthermal      => ffhd_get_pthermal_iso
    else
      phys_get_pthermal      => ffhd_get_pthermal_origin
      ffhd_get_pthermal      => ffhd_get_pthermal_origin
    end if

    ! choose Rfactor in ideal gas law
    if(ffhd_partial_ionization) then
      ffhd_get_Rfactor=>Rfactor_from_temperature_ionization
      phys_update_temperature => ffhd_update_temperature
    else if(associated(usr_Rfactor)) then
      ffhd_get_Rfactor=>usr_Rfactor
    else
      ffhd_get_Rfactor=>Rfactor_from_constant_ionization
    end if

    if(ffhd_partial_ionization) then
      ffhd_get_temperature => ffhd_get_temperature_from_Te
    else
      ffhd_get_temperature => ffhd_get_temperature_from_etot
    end if

    ! derive units from basic units
    call ffhd_physical_units()

    if(ffhd_hyperbolic_thermal_conduction) then
      hypertc_kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
    end if
    if(.not. ffhd_energy .and. ffhd_thermal_conduction) then
      call mpistop("thermal conduction needs ffhd_energy=T")
    end if
    if(.not. ffhd_energy .and. ffhd_hyperbolic_thermal_conduction) then
      call mpistop("hyperbolic thermal conduction needs ffhd_energy=T")
    end if
    if(.not. ffhd_energy .and. ffhd_radiative_cooling) then
      call mpistop("radiative cooling needs ffhd_energy=T")
    end if

    ! initialize thermal conduction module
    if(ffhd_thermal_conduction) then
      phys_req_diagonal = .true.

      call sts_init()
      call tc_init_params(ffhd_gamma)

      allocate(tc_fl)
      call tc_get_ffhd_params(tc_fl,tc_params_read_ffhd)
      call add_sts_method(ffhd_get_tc_dt_ffhd,ffhd_sts_set_source_tc_ffhd,e_,1,e_,1,.false.)
      tc_fl%get_temperature_from_conserved => ffhd_get_temperature_from_etot
      tc_fl%get_temperature_from_eint => ffhd_get_temperature_from_eint
      call set_conversion_methods_to_head(ffhd_e_to_ei, ffhd_ei_to_e)
      call set_error_handling_to_head(ffhd_tc_handle_small_e)
      tc_fl%get_rho => ffhd_get_rho
      tc_fl%e_ = e_
      tc_fl%Tcoff_ = Tcoff_
    end if

    ! Initialize radiative cooling module
    if(ffhd_radiative_cooling) then
      call radiative_cooling_init_params(ffhd_gamma,He_abundance)
      allocate(rc_fl)
      call radiative_cooling_init(rc_fl,rc_params_read)
      rc_fl%get_rho => ffhd_get_rho
      rc_fl%get_pthermal => ffhd_get_pthermal
      rc_fl%get_var_Rfactor => ffhd_get_Rfactor
      rc_fl%e_ = e_
      rc_fl%Tcoff_ = Tcoff_
      rc_fl%has_equi = .false.
    end if
    allocate(te_fl_ffhd)
    te_fl_ffhd%get_rho=> ffhd_get_rho
    te_fl_ffhd%get_pthermal=> ffhd_get_pthermal
    te_fl_ffhd%get_var_Rfactor => ffhd_get_Rfactor
{^IFTHREED
    phys_te_images => ffhd_te_images
}
    ! Initialize viscosity module
    if(ffhd_viscosity) call viscosity_init(phys_wider_stencil,phys_req_diagonal)

    ! Initialize gravity module
    if(ffhd_gravity) then
      call gravity_init()
    end if

    ! initialize ionization degree table
    if(ffhd_partial_ionization) call ionization_degree_init()
  end subroutine ffhd_phys_init

{^IFTHREED
  subroutine ffhd_te_images
    use mod_global_parameters
    use mod_thermal_emission

    select case(convert_type)
      case('EIvtiCCmpi','EIvtuCCmpi')
        call get_EUV_image(unitconvert,te_fl_ffhd)
      case('ESvtiCCmpi','ESvtuCCmpi')
        call get_EUV_spectrum(unitconvert,te_fl_ffhd)
      case('SIvtiCCmpi','SIvtuCCmpi')
        call get_SXR_image(unitconvert,te_fl_ffhd)
      case('WIvtiCCmpi','WIvtuCCmpi')
        call get_whitelight_image(unitconvert,te_fl_ffhd)
      case default
        call mpistop("Error in synthesize emission: Unknown convert_type")
      end select
  end subroutine ffhd_te_images
}

  subroutine  ffhd_sts_set_source_tc_ffhd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_ffhd
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_ffhd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,tc_fl)
  end subroutine ffhd_sts_set_source_tc_ffhd

  function ffhd_get_tc_dt_ffhd(w,ixI^L,ixO^L,dx^D,x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_ffhd
 
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_ffhd(w,ixI^L,ixO^L,dx^D,x,tc_fl)
  end function ffhd_get_tc_dt_ffhd

  subroutine ffhd_tc_handle_small_e(w, x, ixI^L, ixO^L, step)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer, intent(in)    :: step
    character(len=140) :: error_msg

    write(error_msg,"(a,i3)") "Thermal conduction step ", step
    call ffhd_handle_small_ei(w,x,ixI^L,ixO^L,e_,error_msg)
  end subroutine ffhd_tc_handle_small_e

  subroutine tc_params_read_ffhd(fl)
    use mod_global_parameters, only: unitpar,par_files
    type(tc_fluid), intent(inout) :: fl
    integer                      :: n
    ! list parameters
    logical :: tc_saturate=.false.
    double precision :: tc_k_para=0d0
    character(len=std_len)  :: tc_slope_limiter="MC"

    namelist /tc_list/ tc_saturate, tc_slope_limiter, tc_k_para

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, tc_list, end=111)
111     close(unitpar)
    end do

    fl%tc_saturate = tc_saturate
    fl%tc_k_para = tc_k_para
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
  end subroutine tc_params_read_ffhd

  subroutine rc_params_read(fl)
    use mod_global_parameters, only: unitpar,par_files
    use mod_constants, only: bigdouble
    type(rc_fluid), intent(inout) :: fl
    integer                      :: n
    integer :: ncool = 4000
    double precision :: cfrac=0.1d0
  
    !> Name of cooling curve
    character(len=std_len)  :: coolcurve='JCcorona'
  
    !> Name of cooling method
    character(len=std_len)  :: coolmethod='exact'
  
    !> Fixed temperature not lower than tlow
    logical    :: Tfix=.false.
  
    !> Lower limit of temperature
    double precision   :: tlow=bigdouble
  
    !> Add cooling source in a split way (.true.) or un-split way (.false.)
    logical    :: rc_split=.false.
    logical    :: rad_cut=.false.
    double precision :: rad_cut_hgt=0.5d0
    double precision :: rad_cut_dey=0.15d0

    namelist /rc_list/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix, rc_split, rad_cut, rad_cut_hgt, rad_cut_dey

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, rc_list, end=111)
111     close(unitpar)
    end do

    fl%ncool=ncool
    fl%coolcurve=coolcurve
    fl%coolmethod=coolmethod
    fl%tlow=tlow
    fl%Tfix=Tfix
    fl%rc_split=rc_split
    fl%cfrac=cfrac
    fl%rad_cut=rad_cut
    fl%rad_cut_hgt=rad_cut_hgt
    fl%rad_cut_dey=rad_cut_dey
  end subroutine rc_params_read

  subroutine ffhd_check_params
    use mod_global_parameters
    use mod_usr_methods
    use mod_convert, only: add_convert_method

    gamma_1=ffhd_gamma-1.d0
    if (.not. ffhd_energy) then
       if (ffhd_gamma <= 0.0d0) call mpistop ("Error: ffhd_gamma <= 0")
       if (ffhd_adiab < 0.0d0) call mpistop ("Error: ffhd_adiab < 0")
       small_pressure = ffhd_adiab*small_density**ffhd_gamma
    else
       if (ffhd_gamma <= 0.0d0 .or. ffhd_gamma == 1.0d0) &
            call mpistop ("Error: ffhd_gamma <= 0 or ffhd_gamma == 1")
       inv_gamma_1=1.d0/gamma_1
       small_e = small_pressure * inv_gamma_1
    end if

    if (number_equi_vars > 0 .and. .not. associated(usr_set_equi_vars)) then
      call mpistop("usr_set_equi_vars has to be implemented in the user file")
    end if
  end subroutine ffhd_check_params

  subroutine ffhd_physical_units()
    use mod_global_parameters
    double precision :: mp,kB
    double precision :: a,b

    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
    else
      mp=mp_cgs
      kB=kB_cgs
    end if
    if(eq_state_units) then
      a = 1d0 + 4d0 * He_abundance
      if(ffhd_partial_ionization) then
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
      unit_density=a*mp*unit_numberdensity
    end if
    if(unit_velocity/=1.d0) then
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
    else if(unit_pressure/=1.d0) then
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      unit_velocity=sqrt(unit_pressure/unit_density)
    else
      unit_pressure=b*unit_numberdensity*kB*unit_temperature
      unit_velocity=sqrt(unit_pressure/unit_density)
    end if
    if(unit_time/=1.d0) then
      unit_length=unit_time*unit_velocity
    else
      unit_time=unit_length/unit_velocity
    end if
    unit_mass=unit_density*unit_length**3
  end subroutine ffhd_physical_units

  subroutine ffhd_check_w_origin(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters
    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision :: tmp(ixI^S)
    logical, intent(inout) :: flag(ixI^S,1:nw)

    flag=.false.
    where(w(ixO^S,rho_) < small_density) flag(ixO^S,rho_) = .true.

    if(ffhd_energy) then
      if(primitive) then
        where(w(ixO^S,e_) < small_pressure) flag(ixO^S,e_) = .true.
      else
        tmp(ixO^S)=w(ixO^S,e_)-ffhd_kin_en(w,ixI^L,ixO^L)
        where(tmp(ixO^S) < small_e) flag(ixO^S,e_) = .true.
      end if
    end if
  end subroutine ffhd_check_w_origin

  subroutine ffhd_to_conserved_origin(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision :: inv_gamma2(ixO^S)
    integer                         :: idir

    if(ffhd_energy) then
      w(ixO^S,e_)=w(ixO^S,p_)*inv_gamma_1+half*w(ixO^S,mom(1))**2*w(ixO^S,rho_)
    end if
    w(ixO^S,mom(1))=w(ixO^S,rho_)*w(ixO^S,mom(1))
  end subroutine ffhd_to_conserved_origin

  subroutine ffhd_to_primitive_origin(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision                :: inv_rho(ixO^S), gamma2(ixO^S)

    if(fix_small_values) then
      !> fix small values preventing NaN numbers in the following converting
      call ffhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'ffhd_to_primitive_origin')
    end if

    w(ixO^S,mom(1)) = w(ixO^S,mom(1))/w(ixO^S,rho_)
    if(ffhd_energy) then
      w(ixO^S,p_)=gamma_1*(w(ixO^S,e_)-half*w(ixO^S,rho_)*w(ixO^S,mom(1))**2)
    end if
  end subroutine ffhd_to_primitive_origin

  subroutine ffhd_ei_to_e(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    w(ixI^S,e_)=w(ixI^S,e_)+ffhd_kin_en(w,ixI^L,ixI^L)
  end subroutine ffhd_ei_to_e

  subroutine ffhd_e_to_ei(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    w(ixI^S,e_)=w(ixI^S,e_)-ffhd_kin_en(w,ixI^L,ixI^L)
    if(fix_small_values) then
      call ffhd_handle_small_ei(w,x,ixI^L,ixI^L,e_,'ffhd_e_to_ei')
    end if
  end subroutine ffhd_e_to_ei

  subroutine ffhd_handle_small_values_origin(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    logical :: flag(ixI^S,1:nw)
    double precision :: tmp2(ixI^S)

    call phys_check_w(primitive, ixI^L, ixI^L, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S,rho_)) w(ixO^S,rho_) = small_density
        if(small_values_fix_iw(mom(1))) then
          where(flag(ixO^S,rho_)) w(ixO^S, mom(1)) = 0.0d0
        end if
        if(ffhd_energy) then
          if(primitive) then
            where(flag(ixO^S,e_)) w(ixO^S,p_) = small_pressure
          else
            where(flag(ixO^S,e_))
              w(ixO^S,e_) = small_e+ffhd_kin_en(w,ixI^L,ixO^L)
            end where
          end if
        end if
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
        if(ffhd_energy) then
          if(primitive) then
            call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
          else
            w(ixI^S,e_)=w(ixI^S,e_)-ffhd_kin_en(w,ixI^L,ixI^L)
            call small_values_average(ixI^L, ixO^L, w, x, flag, e_)
            w(ixI^S,e_)=w(ixI^S,e_)+ffhd_kin_en(w,ixI^L,ixI^L)
          end if
        end if
      case default
        if(.not.primitive) then
          if(ffhd_energy) then
            w(ixO^S,p_)=gamma_1*(w(ixO^S,e_)-ffhd_kin_en(w,ixI^L,ixO^L))
          end if
          w(ixO^S,mom(1))=w(ixO^S,mom(1))/w(ixO^S,rho_)
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
  end subroutine ffhd_handle_small_values_origin

  subroutine ffhd_get_v_origin(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,ndir)
    double precision :: rho(ixI^S)
    integer :: idir

    call ffhd_get_rho(w,x,ixI^L,ixO^L,rho)
    rho(ixO^S)=1.d0/rho(ixO^S)
    do idir=1,ndir
      v(ixO^S,ndir) = w(ixO^S,mom(1))*block%B0(ixO^S,idir,0)*rho(ixO^S)
    end do
  end subroutine ffhd_get_v_origin

  subroutine ffhd_get_v_idim(w,x,ixI^L,ixO^L,idim,v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S)
    double precision              :: rho(ixI^S)

    call ffhd_get_rho(w,x,ixI^L,ixO^L,rho)
    v(ixO^S) = (w(ixO^S, mom(1))*block%B0(ixO^S,idim,0)) / rho(ixO^S)
  end subroutine ffhd_get_v_idim

  subroutine ffhd_get_cmax_origin(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! w in primitive form
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)

    if(ffhd_energy) then
      cmax(ixO^S)=sqrt(ffhd_gamma*w(ixO^S,p_)/w(ixO^S,rho_))*abs(block%B0(ixO^S,idim,idim))
    else
      cmax(ixO^S)=sqrt(ffhd_gamma*ffhd_adiab*w(ixO^S,rho_)**gamma_1)*abs(block%B0(ixO^S,idim,idim))
    end if
    cmax(ixO^S)=abs(w(ixO^S,mom(1))*block%B0(ixO^S,idim,0))+cmax(ixO^S)

  end subroutine ffhd_get_cmax_origin

  subroutine ffhd_get_cs2max(w,x,ixI^L,ixO^L,cs2max)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cs2max
    double precision :: cs2(ixI^S)

    call ffhd_get_csound2(w,x,ixI^L,ixO^L,cs2)
    cs2max=maxval(cs2(ixO^S))
  end subroutine ffhd_get_cs2max

  subroutine ffhd_get_a2max(w,x,ixI^L,ixO^L,a2max)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
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
  end subroutine ffhd_get_a2max

  subroutine ffhd_get_tcutoff(ixI^L,ixO^L,w,x,Tco_local,Tmax_local)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(out) :: Tco_local,Tmax_local
    double precision, parameter :: trac_delta=0.25d0
    double precision :: tmp1(ixI^S),Te(ixI^S),lts(ixI^S)
    double precision, dimension(ixI^S,1:ndir) :: bunitvec
    double precision, dimension(ixI^S,1:ndim) :: gradT
    double precision :: Bdir(ndim)
    double precision :: ltrc,ltrp,altr(ixI^S)
    integer :: idims,jxO^L,hxO^L,ixA^D,ixB^D
    integer :: jxP^L,hxP^L,ixP^L,ixQ^L
    logical :: lrlt(ixI^S)

    call ffhd_get_temperature(w,x,ixI^L,ixI^L,Te)
    Tco_local=zero
    Tmax_local=maxval(Te(ixO^S))

    {^IFONED
    select case(ffhd_trac_type)
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
      call mpistop("ffhd_trac_type not allowed for 1D simulation")
    end select
    }
    {^NOONED
    select case(ffhd_trac_type)
    case(0)
      !> test case, fixed cutoff temperature
      if(slab_uniform) then
        !> assume cgs units
        block%wextra(ixI^S,Tcoff_)=max(min(3.d5/unit_temperature,6.d5/unit_temperature-3.d-4/unit_temperature*unit_length*x(ixI^S,ndim)),zero)
      else
        block%wextra(ixI^S,Tcoff_)=2.5d5/unit_temperature
      end if
    case(1,4,6)
      do idims=1,ndim
        call gradient(Te,ixI^L,ixO^L,idims,tmp1)
        gradT(ixO^S,idims)=tmp1(ixO^S)
      end do
      bunitvec(ixO^S,:)=block%B0(ixO^S,:,0)
      if(ffhd_trac_type .gt. 1) then
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
      else where
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
        call gradientx(Te,x,ixI^L,hxP^L,idims,gradT(ixI^S,idims),.false.)
        call gradientq(Te,x,ixI^L,jxP^L,idims,gradT(ixI^S,idims))
      end do
      bunitvec(ixP^S,:)=block%B0(ixP^S,:,0)
      lts(ixP^S)=abs(sum(gradT(ixP^S,1:ndim)*bunitvec(ixP^S,1:ndim),dim=ndim+1))/Te(ixP^S)
      if(slab_uniform) then
        lts(ixP^S)=minval(dxlevel)*lts(ixP^S)
      else
        lts(ixP^S)=minval(block%ds(ixP^S,:),dim=ndim+1)*lts(ixP^S)
      end if
      lts(ixP^S)=max(one, (exp(lts(ixP^S))/ltrc)**ltrp)
  
      altr=zero
      ixP^L=ixO^L^LADD1;
      do idims=1,ndim
        hxO^L=ixP^L-kr(idims,^D);
        jxO^L=ixP^L+kr(idims,^D);
        altr(ixP^S)=altr(ixP^S)+0.25d0*(lts(hxO^S)+two*lts(ixP^S)+lts(jxO^S))*bunitvec(ixP^S,idims)**2
      end do
      block%wextra(ixP^S,Tcoff_)=Te(ixP^S)*altr(ixP^S)**0.4d0
    case(3,5)
      !> do nothing here
    case default
      call mpistop("unknown ffhd_trac_type")
    end select
    }
  end subroutine ffhd_get_tcutoff

  subroutine ffhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,Hspeed,cmax,cmin)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
    double precision, intent(in)    :: Hspeed(ixI^S,1:number_species)
    double precision :: wmean(ixI^S,nw)
    double precision, dimension(ixI^S) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix^D

    select case (boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      tmp1(ixO^S)=sqrt(wLp(ixO^S,rho_))
      tmp2(ixO^S)=sqrt(wRp(ixO^S,rho_))
      tmp3(ixO^S)=1.d0/(tmp1(ixO^S)+tmp2(ixO^S))
      umean(ixO^S)=(wLp(ixO^S,mom(1))*block%B0(ixO^S,idim,idim)*tmp1(ixO^S)&
                   +wRp(ixO^S,mom(1))*block%B0(ixO^S,idim,idim)*tmp2(ixO^S))*tmp3(ixO^S)
      call ffhd_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call ffhd_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)
      dmean(ixO^S)=(tmp1(ixO^S)*csoundL(ixO^S)+tmp2(ixO^S)*csoundR(ixO^S)) * &
           tmp3(ixO^S) + 0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2 * &
           (wRp(ixO^S,mom(1))*block%B0(ixO^S,idim,idim)-wLp(ixO^S,mom(1))*block%B0(ixO^S,idim,idim))**2
      dmean(ixO^S)=dsqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,1)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S,1)=umean(ixO^S)+dmean(ixO^S)
      else
        cmax(ixO^S,1)=abs(umean(ixO^S))+dmean(ixO^S)
      end if
    case (2)
      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
      tmp1(ixO^S)=wmean(ixO^S,mom(1))*block%B0(ixO^S,idim,idim)/wmean(ixO^S,rho_)
      call ffhd_get_csound(wmean,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
        cmax(ixO^S,1)=max(tmp1(ixO^S)+csoundR(ixO^S),zero)
        cmin(ixO^S,1)=min(tmp1(ixO^S)-csoundR(ixO^S),zero)
      else
        cmax(ixO^S,1)=abs(tmp1(ixO^S))+csoundR(ixO^S)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call ffhd_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call ffhd_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)
      csoundL(ixO^S)=max(csoundL(ixO^S),csoundR(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,1)=min(wLp(ixO^S,mom(1))*block%B0(ixO^S,idim,idim),&
                          wRp(ixO^S,mom(1))*block%B0(ixO^S,idim,idim))-csoundL(ixO^S)
        cmax(ixO^S,1)=max(wLp(ixO^S,mom(1))*block%B0(ixO^S,idim,idim),&
                          wRp(ixO^S,mom(1))*block%B0(ixO^S,idim,idim))+csoundL(ixO^S)
      else
        cmax(ixO^S,1)=max(wLp(ixO^S,mom(1))*block%B0(ixO^S,idim,idim),&
                          wRp(ixO^S,mom(1))*block%B0(ixO^S,idim,idim))+csoundL(ixO^S)
      end if
    end select
  end subroutine ffhd_get_cbounds

  subroutine ffhd_get_csound(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)

    call ffhd_get_csound2(w,x,ixI^L,ixO^L,csound)
    csound(ixO^S) = dsqrt(csound(ixO^S))*abs(block%B0(ixO^S,idim,idim))
  end subroutine ffhd_get_csound

  !> Calculate fast magnetosonic wave speed
  subroutine ffhd_get_csound_prim(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)

    if(ffhd_energy) then
      csound(ixO^S)=ffhd_gamma*w(ixO^S,e_)/w(ixO^S,rho_)
    else
      csound(ixO^S)=ffhd_gamma*ffhd_adiab*w(ixO^S,rho_)**gamma_1
    end if
    csound(ixO^S) = dsqrt(csound(ixO^S))
  end subroutine ffhd_get_csound_prim

  subroutine ffhd_get_pthermal_iso(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)

    call ffhd_get_rho(w,x,ixI^L,ixO^L,pth)
    pth(ixO^S)=ffhd_adiab*pth(ixO^S)**ffhd_gamma
  end subroutine ffhd_get_pthermal_iso

  subroutine ffhd_get_pthermal_origin(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)
    integer                      :: iw, ix^D

    pth(ixO^S)=gamma_1*(w(ixO^S,e_)-ffhd_kin_en(w,ixI^L,ixO^L))
    if (fix_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(pth(ix^D)<small_pressure) then
            pth(ix^D)=small_pressure
         end if
      {end do^D&\}
    elseif(check_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(pth(ix^D)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                " encountered when call ffhd_get_pthermal"
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
      {end do^D&\}
    end if
  end subroutine ffhd_get_pthermal_origin

  subroutine ffhd_get_temperature_from_Te(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)

    res(ixO^S) = w(ixO^S, Te_)
  end subroutine ffhd_get_temperature_from_Te

  subroutine ffhd_get_temperature_from_eint(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    double precision :: R(ixI^S)

    call ffhd_get_Rfactor(w,x,ixI^L,ixO^L,R)
    res(ixO^S) = gamma_1 * w(ixO^S, e_)/(w(ixO^S,rho_)*R(ixO^S))
  end subroutine ffhd_get_temperature_from_eint

  subroutine ffhd_get_temperature_from_etot(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)

    double precision :: R(ixI^S)

    call ffhd_get_Rfactor(w,x,ixI^L,ixO^L,R)
    call ffhd_get_pthermal(w,x,ixI^L,ixO^L,res)
    res(ixO^S)=res(ixO^S)/(R(ixO^S)*w(ixO^S,rho_))
  end subroutine ffhd_get_temperature_from_etot

  subroutine ffhd_get_csound2(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision    :: rho(ixI^S)
    
    call ffhd_get_rho(w,x,ixI^L,ixO^L,rho)
    if(ffhd_energy) then
      call ffhd_get_pthermal(w,x,ixI^L,ixO^L,csound2)
      csound2(ixO^S)=ffhd_gamma*csound2(ixO^S)/rho(ixO^S)
    else
      csound2(ixO^S)=ffhd_gamma*ffhd_adiab*rho(ixO^S)**gamma_1
    end if
  end subroutine ffhd_get_csound2

  subroutine ffhd_get_flux(wC,w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)
    double precision             :: ptotal(ixO^S)
    double precision             :: tmp(ixI^S)
    integer                      :: idirmin, iw, idir, jdir, kdir
    double precision, dimension(ixI^S) :: Te,tau,sigT

    f(ixO^S,rho_)=w(ixO^S,mom(1))*w(ixO^S,rho_)*block%B0(ixO^S,idim,idim)

    if(ffhd_energy) then
      ptotal(ixO^S)=w(ixO^S,p_)
    else
      ptotal(ixO^S)=ffhd_adiab*w(ixO^S,rho_)**ffhd_gamma
    end if

    ! Get flux of momentum
    f(ixO^S,mom(1))=(wC(ixO^S,mom(1))*w(ixO^S,mom(1))+ptotal(ixO^S))*block%B0(ixO^S,idim,idim)

    ! Get flux of energy
    if(ffhd_energy) then
      f(ixO^S,e_)=w(ixO^S,mom(1))*(wC(ixO^S,e_)+ptotal(ixO^S))*block%B0(ixO^S,idim,idim)
      if(ffhd_hyperbolic_thermal_conduction) then
        f(ixO^S,e_)=f(ixO^S,e_)+w(ixO^S,q_)*block%B0(ixO^S,idim,idim)
        f(ixO^S,q_)=zero
      end if
    end if
  end subroutine ffhd_get_flux

  subroutine ffhd_add_source(qdt,dtfactor,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,dtfactor
    double precision, intent(in)    :: wCT(ixI^S,1:nw),wCTprim(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    if (.not. qsourcesplit) then
      active = .true.
      call add_punitb(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
      if(ffhd_hyperbolic_thermal_conduction) then
        call add_hypertc_source(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
      end if
    end if

    if(ffhd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,&
           w,x,qsourcesplit,active, rc_fl)
    end if

    if(ffhd_viscosity) then
      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,ffhd_energy,qsourcesplit,active)
    end if

    if(ffhd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,&
           w,x,ffhd_energy,.false.,qsourcesplit,active)
    end if

    ! update temperature from new pressure, density, and old ionization degree
    if(ffhd_partial_ionization) then
      if(.not.qsourcesplit) then
        active = .true.
        call ffhd_update_temperature(ixI^L,ixO^L,wCT,w,x)
      end if
    end if
  end subroutine ffhd_add_source

  subroutine add_punitb(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCTprim(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer                         :: idims,hxO^L
    double precision                :: divb(ixI^S)

    divb=zero
    if(slab_uniform) then
      do idims=1,ndim
        hxO^L=ixO^L-kr(idims,^D);
        divb(ixO^S)=divb(ixO^S)+(block%B0(ixO^S,idims,idims)-block%B0(hxO^S,idims,idims))/dxlevel(idims)
      end do
    else
      call divvector(block%B0(ixI^S,1:ndir,0),ixI^L,ixO^L,divb)
    end if
    w(ixO^S,mom(1))=w(ixO^S,mom(1))+qdt*wCTprim(ixO^S,p_)*divb(ixO^S)
  end subroutine add_punitb

  subroutine ffhd_get_rho(w,x,ixI^L,ixO^L,rho)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: rho(ixI^S)

    rho(ixO^S) = w(ixO^S,rho_)
  end subroutine ffhd_get_rho

  subroutine ffhd_handle_small_ei(w, x, ixI^L, ixO^L, ie, subname)
    use mod_global_parameters
    use mod_small_values
    integer, intent(in)             :: ixI^L,ixO^L, ie
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    integer :: idir
    logical :: flag(ixI^S,1:nw)
    double precision              :: rho(ixI^S)

    flag=.false.
    where(w(ixO^S,ie)<small_e) flag(ixO^S,ie)=.true.
    if(any(flag(ixO^S,ie))) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S,ie)) w(ixO^S,ie)=small_e
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag, ie)
      case default
        w(ixO^S,e_)=w(ixO^S,e_)*gamma_1
        call ffhd_get_rho(w,x,ixI^L,ixO^L,rho)
        w(ixO^S,mom(1)) = w(ixO^S,mom(1))/rho(ixO^S)
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
  end subroutine ffhd_handle_small_ei

  subroutine ffhd_update_temperature(ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_ionization_degree
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: iz_H(ixO^S),iz_He(ixO^S), pth(ixI^S)

    call ionization_degree_from_temperature(ixI^L,ixO^L,wCT(ixI^S,Te_),iz_H,iz_He)
    call ffhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,Te_)=(2.d0+3.d0*He_abundance)*pth(ixO^S)/(w(ixO^S,rho_)*(1.d0+iz_H(ixO^S)+&
      He_abundance*(iz_He(ixO^S)*(iz_He(ixO^S)+1.d0)+1.d0)))
  end subroutine ffhd_update_temperature

  subroutine ffhd_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx^D
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer                       :: idirmin,idim
    double precision              :: dxarr(ndim)
    double precision              :: current(ixI^S,7-2*ndir:3),eta(ixI^S)

    dtnew = bigdouble

    if(ffhd_radiative_cooling) then
      call cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x,rc_fl)
    end if

    if(ffhd_viscosity) then
      call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(ffhd_gravity) then
      call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if
  end subroutine ffhd_get_dt

!  subroutine ffhd_add_source_geom(qdt,dtfactor,ixI^L,ixO^L,wCT,w,x)
!    use mod_global_parameters
!    use mod_geometry
!
!    integer, intent(in)             :: ixI^L, ixO^L
!    double precision, intent(in)    :: qdt, dtfactor,x(ixI^S,1:ndim)
!    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
!
!    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
!    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),invrho(ixO^S),invr(ixO^S)
!
!    integer :: mr_,mphi_ ! Polar var. names
!    integer :: br_,bphi_
!
!    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
!    br_=mag(1); bphi_=mag(1)-1+phi_
!
!    ! 1/rho
!    invrho(ixO^S)=1.d0/wCT(ixO^S,rho_)
!    ! include dt in invr, invr is always used with qdt
!    if(local_timestep) then
!      invr(ixO^S) = block%dt(ixO^S) * dtfactor/x(ixO^S,1)
!    else
!      invr(ixO^S) = qdt/x(ixO^S,1)
!    end if  
!
!
!    select case (coordinate)
!    case (cylindrical)
!      call ffhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp)
!      if(phi_>0) then
!        w(ixO^S,mr_)=w(ixO^S,mr_)+invr(ixO^S)*(tmp(ixO^S)-&
!                  wCT(ixO^S,bphi_)**2+wCT(ixO^S,mphi_)**2*invrho(ixO^S))
!        w(ixO^S,mphi_)=w(ixO^S,mphi_)+invr(ixO^S)*(&
!                 -wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)*invrho(ixO^S) &
!                 +wCT(ixO^S,bphi_)*wCT(ixO^S,br_))
!        if(.not.stagger_grid) then
!          w(ixO^S,bphi_)=w(ixO^S,bphi_)+invr(ixO^S)*&
!                   (wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
!                   -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
!                   *invrho(ixO^S)
!        end if
!      else
!        w(ixO^S,mr_)=w(ixO^S,mr_)+invr(ixO^S)*tmp(ixO^S)
!      end if
!      if(ffhd_glm) w(ixO^S,br_)=w(ixO^S,br_)+wCT(ixO^S,psi_)*invr(ixO^S)
!    case (spherical)
!       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
!       call ffhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp1)
!       ! m1
!       tmp(ixO^S)=tmp1(ixO^S)*x(ixO^S,1) &
!                  *(block%surfaceC(ixO^S,1)-block%surfaceC(h1x^S,1))/block%dvolume(ixO^S)
!       do idir=2,ndir
!         tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(idir))**2*invrho(ixO^S)-wCT(ixO^S,mag(idir))**2
!       end do
!       w(ixO^S,mom(1))=w(ixO^S,mom(1))+tmp(ixO^S)*invr(ixO^S)
!       ! b1
!       if(ffhd_glm) then
!         w(ixO^S,mag(1))=w(ixO^S,mag(1))+invr(ixO^S)*2.0d0*wCT(ixO^S,psi_)
!       end if
!
!       {^NOONED
!       ! m2
!       ! This will make hydrostatic p=const an exact solution
!       if(local_timestep) then
!          tmp(ixO^S) = block%dt(ixO^S) * tmp1(ixO^S)
!       else
!          tmp(ixO^S) = qdt * tmp1(ixO^S)
!       end if  
!       w(ixO^S,mom(2))=w(ixO^S,mom(2))+tmp(ixO^S) &
!            *(block%surfaceC(ixO^S,2)-block%surfaceC(h2x^S,2)) &
!            /block%dvolume(ixO^S)
!       tmp(ixO^S)=-(wCT(ixO^S,mom(1))*wCT(ixO^S,mom(2))*invrho(ixO^S) &
!            -wCT(ixO^S,mag(1))*wCT(ixO^S,mag(2)))
!       if(ndir==3) then
!         tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(3))**2*invrho(ixO^S) &
!              -wCT(ixO^S,mag(3))**2)*dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
!       end if
!       w(ixO^S,mom(2))=w(ixO^S,mom(2))+tmp(ixO^S)*invr(ixO^S)
!       ! b2
!       if(.not.stagger_grid) then
!         tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(2)) &
!              -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(1)))*invrho(ixO^S)
!         if(ffhd_glm) then
!           tmp(ixO^S)=tmp(ixO^S) &
!                + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
!         end if
!         w(ixO^S,mag(2))=w(ixO^S,mag(2))+tmp(ixO^S)*invr(ixO^S)
!       end if
!       }
!
!       if(ndir==3) then
!         ! m3
!         tmp(ixO^S)=-(wCT(ixO^S,mom(3))*wCT(ixO^S,mom(1))*invrho(ixO^S) &
!              -wCT(ixO^S,mag(3))*wCT(ixO^S,mag(1))) {^NOONED &
!              -(wCT(ixO^S,mom(2))*wCT(ixO^S,mom(3))*invrho(ixO^S) &
!              -wCT(ixO^S,mag(2))*wCT(ixO^S,mag(3))) &
!              *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
!         w(ixO^S,mom(3))=w(ixO^S,mom(3))+tmp(ixO^S)*invr(ixO^S)
!         ! b3
!         if(.not.stagger_grid) then
!           tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(3)) &
!                -wCT(ixO^S,mom(3))*wCT(ixO^S,mag(1)))*invrho(ixO^S) {^NOONED &
!                -(wCT(ixO^S,mom(3))*wCT(ixO^S,mag(2)) &
!                -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(3)))*dcos(x(ixO^S,2)) &
!                /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
!           w(ixO^S,mag(3))=w(ixO^S,mag(3))+tmp(ixO^S)*invr(ixO^S)
!         end if
!       end if
!    end select
!  end subroutine ffhd_add_source_geom

  function ffhd_kin_en_origin(w, ixI^L, ixO^L, inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim,block
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)
    double precision, intent(in), optional :: inv_rho(ixO^S)

    if(present(inv_rho)) then
      ke(ixO^S)=0.5d0*w(ixO^S,mom(1))**2*inv_rho(ixO^S)
    else
      ke(ixO^S)=0.5d0*w(ixO^S,mom(1))**2/w(ixO^S,rho_)
    end if
  end function ffhd_kin_en_origin

  subroutine Rfactor_from_temperature_ionization(w,x,ixI^L,ixO^L,Rfactor)
    use mod_global_parameters
    use mod_ionization_degree
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: Rfactor(ixI^S)
    double precision :: iz_H(ixO^S),iz_He(ixO^S)

    call ionization_degree_from_temperature(ixI^L,ixO^L,w(ixI^S,Te_),iz_H,iz_He)
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

  subroutine get_tau(ixI^L,ixO^L,w,Te,tau,sigT5)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:nw), intent(in) :: w
    double precision, dimension(ixI^S), intent(in) :: Te
    double precision, dimension(ixI^S), intent(out) :: tau,sigT5
    integer :: ix^D
    double precision :: dxmin,taumin
    double precision, dimension(ixI^S) :: sigT7,eint

    taumin=4.d0
    !> w supposed to be wCTprim here
    if(ffhd_trac) then
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
    eint(ixO^S)=w(ixO^S,p_)/(ffhd_gamma-one)
    tau(ixO^S)=max(taumin*dt,sigT7(ixO^S)/eint(ixO^S)/cs2max_global)
  end subroutine get_tau

  subroutine add_hypertc_source(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qdt
    double precision, dimension(ixI^S,1:ndim), intent(in) :: x
    double precision, dimension(ixI^S,1:nw), intent(in) :: wCT,wCTprim
    double precision, dimension(ixI^S,1:nw), intent(inout) :: w
    integer :: idims
    integer :: hxC^L,hxO^L,ixC^L,jxC^L,jxO^L,kxC^L
    double precision :: invdx
    double precision, dimension(ixI^S) :: Te,tau,sigT,htc_qsrc,Tface
    double precision, dimension(ixI^S) :: htc_esrc

    Te(ixI^S)=wCTprim(ixI^S,p_)/wCT(ixI^S,rho_)
    call get_tau(ixI^L,ixO^L,wCTprim,Te,tau,sigT)
    htc_qsrc=zero
    do idims=1,ndim
      invdx=1.d0/dxlevel(idims)
      ixC^L=ixO^L;
      ixCmin^D=ixOmin^D-kr(idims,^D);ixCmax^D=ixOmax^D;
      jxC^L=ixC^L+kr(idims,^D);
      kxC^L=jxC^L+kr(idims,^D);
      hxC^L=ixC^L-kr(idims,^D);
      hxO^L=ixO^L-kr(idims,^D);
      Tface(ixC^S)=(7.d0*(Te(ixC^S)+Te(jxC^S))-(Te(hxC^S)+Te(kxC^S)))/12.d0
      htc_qsrc(ixO^S)=htc_qsrc(ixO^S)+sigT(ixO^S)*block%B0(ixO^S,idims,0)*(Tface(ixO^S)-Tface(hxO^S))*invdx
    end do
    htc_qsrc(ixO^S)=(htc_qsrc(ixO^S)+wCT(ixO^S,q_))/tau(ixO^S)
    w(ixO^S,q_)=w(ixO^S,q_)-qdt*htc_qsrc(ixO^S)
  end subroutine add_hypertc_source
end module mod_ffhd_phys
