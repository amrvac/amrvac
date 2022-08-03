!> Magneto-hydrodynamics module
module mod_mhd_phys

#include "amrvac.h"

  use mod_global_parameters, only: std_len, const_c
  use mod_thermal_conduction, only: tc_fluid
  use mod_radiative_cooling, only: rc_fluid
  use mod_thermal_emission, only: te_fluid
  use mod_physics

  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: mhd_energy = .true.

  !> Whether thermal conduction is used
  logical, public, protected              :: mhd_thermal_conduction = .false.
  !> type of fluid for thermal conduction
  type(tc_fluid), public, allocatable     :: tc_fl
  type(te_fluid), public, allocatable     :: te_fl_mhd

  !> Whether radiative cooling is added
  logical, public, protected              :: mhd_radiative_cooling = .false.
  !> type of fluid for radiative cooling
  type(rc_fluid), public, allocatable     :: rc_fl

  !> Whether viscosity is added
  logical, public, protected              :: mhd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: mhd_gravity = .false.

  !> Whether Hall-MHD is used
  logical, public, protected              :: mhd_Hall = .false.

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

  !> Height of the mask used in the TRAC method
  double precision, public, protected     :: mhd_trac_mask = 0.d0

  !> Distance between two adjacent traced magnetic field lines (in finest cell size)   
  integer, public, protected              :: mhd_trac_finegrid=4

  !> Whether auxiliary internal energy is solved
  logical, public, protected              :: mhd_solve_eaux = .false.

  !> Whether internal energy is solved instead of total energy
  logical, public, protected              :: mhd_internal_e = .false.

  !TODO this does not work with the splitting: check mhd_check_w_hde and mhd_handle_small_values_hde
  !> Whether hydrodynamic energy is solved instead of total energy
  logical, public, protected              :: mhd_hydrodynamic_e = .false.

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: mhd_glm_alpha = 0.5d0

  !TODO this does not work with the splitting: check mhd_check_w_semirelati and mhd_handle_small_values_semirelati
  !> Whether semirelativistic MHD equations (Gombosi 2002 JCP) are solved 
  logical, public, protected              :: mhd_semirelativistic = .false.

  !> Whether boris simplified semirelativistic MHD equations (Gombosi 2002 JCP) are solved 
  logical, public, protected              :: mhd_boris_simplification = .false.

  !> Reduced speed of light for semirelativistic MHD
  double precision, public, protected     :: mhd_reduced_c = const_c

  !> MHD fourth order
  logical, public, protected              :: mhd_4th_order = .false.

  !> whether split off equilibrium density
  logical, public :: has_equi_rho0 = .false.  
  !> whether split off equilibrium thermal pressure
  logical, public :: has_equi_pe0 = .false.  
  logical, public :: mhd_equi_thermal = .false.  

  !> equi vars indices in the state%equi_vars array
  integer, public :: equi_rho0_ = -1
  integer, public :: equi_pe0_ = -1

  !> whether dump full variables (when splitting is used) in a separate dat file
  logical, public, protected              :: mhd_dump_full_vars = .false.

  !> Number of tracer species
  integer, public, protected              :: mhd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> Indices of the magnetic field
  integer, allocatable, public, protected :: mag(:)

  !> Indices of the GLM psi
  integer, public, protected :: psi_

  !> Indices of auxiliary internal energy
  integer, public, protected :: eaux_
  integer, public, protected :: paux_

  !> Index of the cutoff temperature for the TRAC method
  integer, public, protected              :: Tcoff_
  integer, public, protected              :: Tweight_

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> The adiabatic index
  double precision, public                :: mhd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: mhd_adiab = 1.0d0

  !> The MHD resistivity
  double precision, public                :: mhd_eta = 0.0d0

  !> The MHD hyper-resistivity
  double precision, public                :: mhd_eta_hyper = 0.0d0

  !> TODO: what is this?
  double precision, public                :: mhd_etah = 0.0d0

  !> The MHD ambipolar coefficient
  double precision, public                :: mhd_eta_ambi = 0.0d0

  !> The small_est allowed energy
  double precision, protected             :: small_e

  !> The number of waves
  integer :: nwwave=8

  !> Method type to clean divergence of B
  character(len=std_len), public, protected :: typedivbfix  = 'linde'

  !> Method type of constrained transport
  character(len=std_len), public, protected :: type_ct  = 'uct_contact'

  !> Whether divB is computed with a fourth order approximation
  logical, public, protected :: mhd_divb_4thorder = .false.

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

  !> clean initial divB
  logical, public :: clean_initial_divb     = .false.

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

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*^ND)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public, protected :: boundary_divbfix_skip(2*^ND)=0

  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.

  !> Whether an total energy equation is used
  logical :: total_energy = .true.

  !> Whether unsplit semirelativistic MHD is solved
  logical :: unsplit_semirelativistic=.false.

  !> gamma minus one and its inverse
  double precision :: gamma_1, inv_gamma_1

  !> inverse of squared speed of light c0 and reduced speed of light c
  double precision :: inv_squared_c0, inv_squared_c

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

  !define the subroutine interface for the ambipolar mask
  abstract interface

    subroutine mask_subroutine(ixI^L,ixO^L,w,x,res)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in) :: w(ixI^S,1:nw)
      double precision, intent(inout) :: res(ixI^S)
    end subroutine mask_subroutine

    function fun_kin_en(w, ixI^L, ixO^L, inv_rho) result(ke)
      use mod_global_parameters, only: nw, ndim,block
      integer, intent(in)           :: ixI^L, ixO^L
      double precision, intent(in)  :: w(ixI^S, nw)
      double precision              :: ke(ixO^S)
      double precision, intent(in), optional :: inv_rho(ixO^S)
    end function fun_kin_en

  end interface

  procedure(mask_subroutine), pointer  :: usr_mask_ambipolar => null()
  procedure(sub_convert), pointer      :: mhd_to_primitive  => null()
  procedure(sub_convert), pointer      :: mhd_to_conserved  => null()
  procedure(sub_small_values), pointer :: mhd_handle_small_values => null()
  procedure(sub_get_pthermal), pointer :: mhd_get_pthermal  => null()
  procedure(sub_get_v), pointer        :: mhd_get_v         => null()
  procedure(fun_kin_en), pointer       :: mhd_kin_en        => null()
  ! Public methods
  public :: usr_mask_ambipolar
  public :: mhd_phys_init
  public :: mhd_kin_en
  public :: mhd_get_pthermal
  public :: mhd_get_v
  public :: mhd_get_v_idim
  public :: mhd_to_conserved
  public :: mhd_to_primitive
  public :: mhd_get_csound2
  public :: mhd_face_to_center
  public :: get_divb
  public :: get_current
  !> needed  public if we want to use the ambipolar coefficient in the user file
  public :: multiplyAmbiCoef  
  public :: get_normalized_divb
  public :: b_from_vector_potential
  public :: mhd_mag_en_all
  {^NOONED
  public :: mhd_clean_divb_multigrid
  }

contains

  !> Read this module"s parameters from a file
  subroutine mhd_read_params(files)
    use mod_global_parameters
    use mod_particles, only: particles_eta, particles_etah
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mhd_list/ mhd_energy, mhd_n_tracer, mhd_gamma, mhd_adiab,&
      mhd_eta, mhd_eta_hyper, mhd_etah, mhd_eta_ambi, mhd_glm_alpha, mhd_glm_extended, mhd_magnetofriction,&
      mhd_thermal_conduction, mhd_radiative_cooling, mhd_Hall, mhd_ambipolar, mhd_ambipolar_sts, mhd_gravity,&
      mhd_viscosity, mhd_4th_order, typedivbfix, source_split_divb, divbdiff,&
      typedivbdiff, type_ct, compactres, divbwave, He_abundance, &
      H_ion_fr, He_ion_fr, He_ion_fr2, eq_state_units, SI_unit, B0field ,mhd_dump_full_vars,&
      B0field_forcefree, Bdip, Bquad, Boct, Busr, mhd_particles,&
      particles_eta, particles_etah,has_equi_rho0, has_equi_pe0,mhd_equi_thermal,&
      boundary_divbfix, boundary_divbfix_skip, mhd_divb_4thorder, mhd_semirelativistic,&
      mhd_boris_simplification, mhd_reduced_c, clean_initial_divb, mhd_solve_eaux, mhd_internal_e, &
      mhd_hydrodynamic_e, mhd_trac, mhd_trac_type, mhd_trac_mask, mhd_trac_finegrid

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mhd_list, end=111)
111    close(unitpar)
    end do
 
  end subroutine mhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine mhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = mhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine mhd_write_info

  subroutine mhd_angmomfix(fC,x,wnew,ixI^L,ixO^L,idim)
    use mod_global_parameters
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(inout)    :: fC(ixI^S,1:nwflux,1:ndim),  wnew(ixI^S,1:nw)
    integer, intent(in)                :: ixI^L, ixO^L
    integer, intent(in)                :: idim
    integer                            :: hxO^L, kxC^L, iw
    double precision                   :: inv_volume(ixI^S)

    call mpistop("to do")

  end subroutine mhd_angmomfix

  subroutine mhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init, particles_eta, particles_etah
    use mod_magnetofriction, only: magnetofriction_init
    use mod_supertimestepping, only: sts_init, add_sts_method,&
            set_conversion_methods_to_head, set_error_handling_to_head
    {^NOONED
    use mod_multigrid_coupling
    }

    integer :: itr, idir

    call mhd_read_params(par_files)

    if(mhd_internal_e) then
      if(mhd_solve_eaux) then
        mhd_solve_eaux=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_solve_eaux=F when mhd_internal_e=T'
      end if
      if(mhd_hydrodynamic_e) then
        mhd_hydrodynamic_e=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_hydrodynamic_e=F when mhd_internal_e=T'
      end if
    end if

    if(mhd_semirelativistic) then
      unsplit_semirelativistic=.true.
      if(mhd_internal_e) then
        mhd_internal_e=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_internal_e=F when mhd_semirelativistic=T'
      end if
      if(mhd_hydrodynamic_e) then
        ! use split semirelativistic MHD
        unsplit_semirelativistic=.false.
      end if
      if(mhd_boris_simplification) then
        mhd_boris_simplification=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_boris_simplification=F when mhd_semirelativistic=T'
      end if
    end if

    if(.not. mhd_energy) then
      if(mhd_internal_e) then
        mhd_internal_e=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_internal_e=F when mhd_energy=F'
      end if
      if(mhd_solve_eaux) then
        mhd_solve_eaux=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_solve_eaux=F when mhd_energy=F'
      end if
      if(mhd_hydrodynamic_e) then
        mhd_hydrodynamic_e=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_hydrodynamic_e=F when mhd_energy=F'
      end if
      if(mhd_thermal_conduction) then
        mhd_thermal_conduction=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_thermal_conduction=F when mhd_energy=F'
      end if
      if(mhd_radiative_cooling) then
        mhd_radiative_cooling=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_radiative_cooling=F when mhd_energy=F'
      end if
      if(mhd_trac) then
        mhd_trac=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_trac=F when mhd_energy=F'
      end if
    end if

    physics_type = "mhd"
    phys_energy=mhd_energy
    phys_internal_e=mhd_internal_e
    phys_solve_eaux=mhd_solve_eaux
    phys_trac=mhd_trac
    phys_trac_type=mhd_trac_type

    phys_gamma = mhd_gamma

    if(mhd_energy.and..not.mhd_internal_e.and..not.mhd_hydrodynamic_e) then
      total_energy=.true.
    else
      total_energy=.false.
    end if
    phys_total_energy=total_energy

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
    phys_trac_mask=mhd_trac_mask

    if(mhd_solve_eaux) prolongprimitive=.true.

    ! set default gamma for polytropic/isothermal process
    use_particles=mhd_particles
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
       type_divb = divb_none
    {^NOONED
    case ('multigrid')
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

    !  set auxiliary internal energy variable
    if(mhd_energy .and. mhd_solve_eaux) then
      eaux_ = var_set_internal_energy()
      paux_ = eaux_
    else
      eaux_ = -1
      paux_ = -1
    end if

    if (mhd_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    allocate(tracer(mhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, mhd_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! set number of variables which need update ghostcells
    nwgc=nwflux

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
    if(stagger_grid) nws=ndim

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
      ! for flux of eaux and tracers, using hll flux
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

    phys_get_dt              => mhd_get_dt
    if(mhd_semirelativistic) then
      phys_get_cmax            => mhd_get_cmax_semirelati
    else
      phys_get_cmax            => mhd_get_cmax_origin
    end if
    phys_get_a2max           => mhd_get_a2max
    phys_get_tcutoff         => mhd_get_tcutoff
    phys_get_H_speed         => mhd_get_H_speed
    if(has_equi_rho0) then
      phys_get_cbounds         => mhd_get_cbounds_split_rho
    else if(mhd_semirelativistic) then
      phys_get_cbounds         => mhd_get_cbounds_semirelati
    else
      phys_get_cbounds         => mhd_get_cbounds
    end if
    if(has_equi_rho0) then
      phys_to_primitive        => mhd_to_primitive_split_rho
      mhd_to_primitive         => mhd_to_primitive_split_rho
      phys_to_conserved        => mhd_to_conserved_split_rho
      mhd_to_conserved         => mhd_to_conserved_split_rho
    else if(mhd_internal_e) then
      phys_to_primitive        => mhd_to_primitive_inte
      mhd_to_primitive         => mhd_to_primitive_inte
      phys_to_conserved        => mhd_to_conserved_inte
      mhd_to_conserved         => mhd_to_conserved_inte
    else if(unsplit_semirelativistic) then
      phys_to_primitive        => mhd_to_primitive_semirelati
      mhd_to_primitive         => mhd_to_primitive_semirelati
      phys_to_conserved        => mhd_to_conserved_semirelati
      mhd_to_conserved         => mhd_to_conserved_semirelati
    else if(mhd_hydrodynamic_e) then
      phys_to_primitive        => mhd_to_primitive_hde
      mhd_to_primitive         => mhd_to_primitive_hde
      phys_to_conserved        => mhd_to_conserved_hde
      mhd_to_conserved         => mhd_to_conserved_hde
    else
      phys_to_primitive        => mhd_to_primitive_origin
      mhd_to_primitive         => mhd_to_primitive_origin
      phys_to_conserved        => mhd_to_conserved_origin
      mhd_to_conserved         => mhd_to_conserved_origin
    end if
    if(unsplit_semirelativistic) then
      phys_get_flux            => mhd_get_flux_semirelati
    else
      if(B0field.or.has_equi_rho0.or.has_equi_pe0) then
        phys_get_flux            => mhd_get_flux_split
      else if(mhd_hydrodynamic_e) then
        phys_get_flux            => mhd_get_flux_hde
      else
        phys_get_flux            => mhd_get_flux
      end if
    end if
    if(mhd_boris_simplification) then
      phys_get_v                 => mhd_get_v_boris
      mhd_get_v                  => mhd_get_v_boris
      mhd_kin_en                 => mhd_kin_en_boris
    else
      phys_get_v                 => mhd_get_v_origin
      mhd_get_v                  => mhd_get_v_origin
      mhd_kin_en                 => mhd_kin_en_origin
    end if
    if(B0field.or.has_equi_rho0) then
      phys_add_source_geom     => mhd_add_source_geom_split
    else
      phys_add_source_geom     => mhd_add_source_geom
    end if
    phys_add_source          => mhd_add_source
    phys_check_params        => mhd_check_params
    phys_write_info          => mhd_write_info
    phys_angmomfix           => mhd_angmomfix
    if(unsplit_semirelativistic) then
      phys_handle_small_values => mhd_handle_small_values_semirelati
      mhd_handle_small_values  => mhd_handle_small_values_semirelati
      phys_check_w             => mhd_check_w_semirelati
      phys_get_pthermal        => mhd_get_pthermal_semirelati
      mhd_get_pthermal         => mhd_get_pthermal_semirelati
    else if(mhd_hydrodynamic_e) then
      phys_handle_small_values => mhd_handle_small_values_hde
      mhd_handle_small_values  => mhd_handle_small_values_hde
      phys_check_w             => mhd_check_w_hde
      phys_get_pthermal        => mhd_get_pthermal_hde
      mhd_get_pthermal         => mhd_get_pthermal_hde
    else
      phys_handle_small_values => mhd_handle_small_values_origin
      mhd_handle_small_values  => mhd_handle_small_values_origin
      phys_check_w             => mhd_check_w_origin
      phys_get_pthermal        => mhd_get_pthermal_origin
      mhd_get_pthermal         => mhd_get_pthermal_origin
    end if
    phys_energy_synchro      => mhd_energy_synchro
    if(number_equi_vars>0) then
      phys_set_equi_vars => set_equi_vars_grid
    endif

    if(type_divb==divb_glm) then
      phys_modify_wLR => mhd_modify_wLR
    end if

    ! if using ct stagger grid, boundary divb=0 is not done here
    if(stagger_grid) then
      phys_get_ct_velocity => mhd_get_ct_velocity
      phys_update_faces => mhd_update_faces
      phys_face_to_center => mhd_face_to_center
      phys_modify_wLR => mhd_modify_wLR
    else if(ndim>1) then
      phys_boundary_adjust => mhd_boundary_adjust
    end if

    {^NOONED
    ! clean initial divb
    if(clean_initial_divb) phys_clean_divb => mhd_clean_divb_multigrid
    }

    ! Whether diagonal ghost cells are required for the physics
    if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! derive units from basic units
    call mhd_physical_units()

    if(mhd_semirelativistic.or.mhd_boris_simplification) then
      inv_squared_c0=(unit_velocity/const_c)**2
      inv_squared_c=(unit_velocity/mhd_reduced_c)**2
    end if

    if(.not. mhd_energy .and. mhd_thermal_conduction) then
      call mpistop("thermal conduction needs mhd_energy=T")
    end if
    if(.not. mhd_energy .and. mhd_radiative_cooling) then
      call mpistop("radiative cooling needs mhd_energy=T")
    end if

    ! resistive MHD needs diagonal ghost cells
    if(mhd_eta/=0.d0) phys_req_diagonal = .true.

    ! initialize thermal conduction module
    if (mhd_thermal_conduction) then
      phys_req_diagonal = .true.

      call sts_init()
      call tc_init_params(mhd_gamma)

      allocate(tc_fl)
      call tc_get_mhd_params(tc_fl,tc_params_read_mhd)
      call add_sts_method(mhd_get_tc_dt_mhd,mhd_sts_set_source_tc_mhd,e_,1,e_,1,.false.)
      if(phys_internal_e) then
        if(has_equi_pe0 .and. has_equi_rho0) then
          tc_fl%get_temperature_from_conserved => mhd_get_temperature_from_eint_with_equi
        else
          tc_fl%get_temperature_from_conserved => mhd_get_temperature_from_eint
        end if
      else if(mhd_hydrodynamic_e) then
        tc_fl%get_temperature_from_conserved => mhd_get_temperature_from_hde
      else
        if(has_equi_pe0 .and. has_equi_rho0) then
          tc_fl%get_temperature_from_conserved => mhd_get_temperature_from_etot_with_equi
        else
          tc_fl%get_temperature_from_conserved => mhd_get_temperature_from_etot
        end if
      end if
      if(mhd_solve_eaux) then
        call set_conversion_methods_to_head(mhd_e_to_ei_aux, mhd_ei_to_e_aux)
      else if(unsplit_semirelativistic) then
        call set_conversion_methods_to_head(mhd_e_to_ei_semirelati, mhd_ei_to_e_semirelati)
      else if(mhd_hydrodynamic_e) then
        call set_conversion_methods_to_head(mhd_e_to_ei_hde, mhd_ei_to_e_hde)
      else
        call set_conversion_methods_to_head(mhd_e_to_ei, mhd_ei_to_e)
      end if
      if(has_equi_pe0 .and. has_equi_rho0) then
        tc_fl%get_temperature_from_eint => mhd_get_temperature_from_eint_with_equi
        if(mhd_equi_thermal) then
          tc_fl%has_equi = .true.
          tc_fl%get_temperature_equi => mhd_get_temperature_equi
          tc_fl%get_rho_equi => mhd_get_rho_equi
        else  
          tc_fl%has_equi = .false.
        endif
      else
        tc_fl%get_temperature_from_eint => mhd_get_temperature_from_eint
      endif
      call set_error_handling_to_head(mhd_tc_handle_small_e)
      tc_fl%get_rho => mhd_get_rho
      tc_fl%e_ = e_
      tc_fl%Tcoff_ = Tcoff_
    end if

    ! Initialize radiative cooling module
    if (mhd_radiative_cooling) then
      call radiative_cooling_init_params(mhd_gamma,He_abundance)
      allocate(rc_fl)
      call radiative_cooling_init(rc_fl,rc_params_read)
      rc_fl%get_rho => mhd_get_rho
      rc_fl%get_pthermal => mhd_get_pthermal
      rc_fl%e_ = e_
      rc_fl%eaux_ = eaux_
      rc_fl%Tcoff_ = Tcoff_
      if(has_equi_pe0 .and. has_equi_rho0 .and. mhd_equi_thermal) then
        rc_fl%has_equi = .true.
        rc_fl%get_rho_equi => mhd_get_rho_equi
        rc_fl%get_pthermal_equi => mhd_get_pe_equi
      else
        rc_fl%has_equi = .false.
      end if
    end if
    allocate(te_fl_mhd)
    te_fl_mhd%get_rho=> mhd_get_rho
    te_fl_mhd%get_pthermal=> mhd_get_pthermal
    te_fl_mhd%Rfactor = 1d0
{^IFTHREED
    phys_te_images => mhd_te_images
}
    ! Initialize viscosity module
    if (mhd_viscosity) call viscosity_init(phys_wider_stencil,phys_req_diagonal)

    ! Initialize gravity module
    if(mhd_gravity) then
      call gravity_init()
    end if

    ! Initialize particles module
    if(mhd_particles) then
      call particles_init()
      phys_req_diagonal = .true.
      if(mype==0) then
         write(*,*) '*****Using particles:        with mhd_eta, mhd_etah :', mhd_eta, mhd_etah
         write(*,*) '*****Using particles: particles_eta, particles_etah :', particles_eta, particles_etah
      end if
    end if

    ! initialize magnetofriction module
    if(mhd_magnetofriction) then
      phys_req_diagonal = .true.
      call magnetofriction_init()
    end if

    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in mhd_get_flux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in nghostcells.
    if(mhd_hall) then
      phys_req_diagonal = .true.
      if(mhd_4th_order) then
        phys_wider_stencil = 2
      else
        phys_wider_stencil = 1
      end if
    end if

    if(mhd_ambipolar) then
      phys_req_diagonal = .true.
      if(mhd_ambipolar_sts) then
        call sts_init()
        if(mhd_internal_e) then
          call add_sts_method(get_ambipolar_dt,sts_set_source_ambipolar,mag(1),&
               ndir,mag(1),ndir,.true.)
        else
          call add_sts_method(get_ambipolar_dt,sts_set_source_ambipolar,mom(ndir)+1,&
               mag(ndir)-mom(ndir),mag(1),ndir,.true.)
        end if
      else
        mhd_ambipolar_exp=.true.
        ! For flux ambipolar term, we need one more reconstructed layer since currents are computed
        ! in mhd_get_flux: assuming one additional ghost layer (two for FOURTHORDER) was
        ! added in nghostcells.
        if(mhd_4th_order) then
          phys_wider_stencil = 2
        else
          phys_wider_stencil = 1
        end if
      end if
    end if

  end subroutine mhd_phys_init

{^IFTHREED
  subroutine mhd_te_images()
    use mod_global_parameters
    use mod_thermal_emission

    select case(convert_type)
      case('EIvtiCCmpi','EIvtuCCmpi')
        call get_EUV_image(unitconvert,te_fl_mhd)
      case('ESvtiCCmpi','ESvtuCCmpi')
        call get_EUV_spectrum(unitconvert,te_fl_mhd)
      case('SIvtiCCmpi','SIvtuCCmpi')
        call get_SXR_image(unitconvert,te_fl_mhd)
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

  subroutine mhd_tc_handle_small_e(w, x, ixI^L, ixO^L, step)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer, intent(in)    :: step
    character(len=140) :: error_msg

    write(error_msg,"(a,i3)") "Thermal conduction step ", step
    call mhd_handle_small_ei(w,x,ixI^L,ixO^L,e_,error_msg)
  end subroutine mhd_tc_handle_small_e

  ! fill in tc_fluid fields from namelist
  subroutine tc_params_read_mhd(fl)
    use mod_global_parameters, only: unitpar,par_files
    type(tc_fluid), intent(inout) :: fl

    integer                      :: n

    ! list parameters
    logical :: tc_perpendicular=.true.
    logical :: tc_saturate=.false.
    double precision :: tc_k_para=0d0
    double precision :: tc_k_perp=0d0
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
    fl%tc_slope_limiter = tc_slope_limiter
  end subroutine tc_params_read_mhd
!!end th cond

!!rad cool
  subroutine rc_params_read(fl)
    use mod_global_parameters, only: unitpar,par_files
    use mod_constants, only: bigdouble
    type(rc_fluid), intent(inout) :: fl
    integer                      :: n
    ! list parameters
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


    namelist /rc_list/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix, rc_split

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

  end subroutine rc_params_read
!! end rad cool

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
    double precision   :: rho(ixI^S)

    call  mhd_get_rho(w,x,ixI^L,ixO^L,rho(ixI^S))
    wnew(ixO^S,rho_) = rho(ixO^S)
    wnew(ixO^S,mom(:)) =  w(ixO^S,mom(:))

    if (B0field) then
      ! add background magnetic field B0 to B
      wnew(ixO^S,mag(:))=w(ixO^S,mag(:))+block%B0(ixO^S,:,0)
    else
      wnew(ixO^S,mag(:))=w(ixO^S,mag(:))
    end if

    if(mhd_energy) then
      wnew(ixO^S,e_) = w(ixO^S,e_)
      if(has_equi_pe0) then
        wnew(ixO^S,e_) = wnew(ixO^S,e_) + block%equi_vars(ixO^S,equi_pe0_,0)* inv_gamma_1
      end if
      if(B0field .and. .not. mhd_internal_e) then
          wnew(ixO^S,e_)=wnew(ixO^S,e_)+0.5d0*sum(block%B0(ixO^S,:,0)**2,dim=ndim+1) &
              + sum(w(ixO^S,mag(:))*block%B0(ixO^S,:,0),dim=ndim+1)
      end if
    end if

  end function convert_vars_splitting

  subroutine mhd_check_params
    use mod_global_parameters
    use mod_usr_methods
    use mod_convert, only: add_convert_method

    ! after user parameter setting
    gamma_1=mhd_gamma-1.d0
    if (.not. mhd_energy) then
       if (mhd_gamma <= 0.0d0) call mpistop ("Error: mhd_gamma <= 0")
       if (mhd_adiab < 0.0d0) call mpistop ("Error: mhd_adiab < 0")
       small_pressure = mhd_adiab*small_density**mhd_gamma
    else
       if (mhd_gamma <= 0.0d0 .or. mhd_gamma == 1.0d0) &
            call mpistop ("Error: mhd_gamma <= 0 or mhd_gamma == 1")
       inv_gamma_1=1.d0/gamma_1
       small_e = small_pressure * inv_gamma_1
    end if

    if (number_equi_vars > 0 .and. .not. associated(usr_set_equi_vars)) then
      call mpistop("usr_set_equi_vars has to be implemented in the user file")
    endif
    if(convert .or. autoconvert) then
      if(convert_type .eq. 'dat_generic_mpi') then
        if(mhd_dump_full_vars) then
          if(mype .eq. 0) print*, " add conversion method: split -> full "
          call add_convert_method(convert_vars_splitting, nw, cons_wnames, "new")
        endif
      endif
    endif
  end subroutine mhd_check_params

  subroutine mhd_physical_units()
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
      b = 1d0 + H_ion_fr + He_abundance*(He_ion_fr*(He_ion_fr2 + 1d0)+1d0)
      RR = 1d0
    else
      a = 1d0
      b = 1d0
      RR = (1d0 + H_ion_fr + He_abundance*(He_ion_fr*(He_ion_fr2 + 1d0)+1d0))/(1d0 + 4d0 * He_abundance)   
    endif
    if(unit_numberdensity/=1.d0) then
      unit_density=a*mp*unit_numberdensity
    else if(unit_density/=1.d0) then
      unit_numberdensity=unit_density/(a*mp)
    end if
    if(unit_temperature/=1.d0) then
      unit_pressure=b*unit_numberdensity*kB*unit_temperature
      unit_velocity=sqrt(unit_pressure/unit_density)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    else if(unit_velocity/=1.d0) then
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
    end if
    if(unit_length/=1.d0) then 
      unit_time=unit_length/unit_velocity
    else if(unit_time/=1.d0) then
      unit_length=unit_time*unit_velocity
    end if
    ! Additional units needed for the particles
    c_norm=c_lightspeed/unit_velocity
    unit_charge=unit_magneticfield*unit_length**2/unit_velocity/miu0
    if (.not. SI_unit) unit_charge = unit_charge*const_c
    unit_mass=unit_density*unit_length**3

  end subroutine mhd_physical_units

  subroutine mhd_check_w_semirelati(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision :: tmp(ixO^S),b2(ixO^S),b(ixO^S,1:ndir),Ba(ixO^S,1:ndir)
    double precision :: v(ixO^S,1:ndir),gamma2(ixO^S),inv_rho(ixO^S)
    logical, intent(inout) :: flag(ixI^S,1:nw)

    integer :: idir, jdir, kdir

    flag=.false.
    where(w(ixO^S,rho_) < small_density) flag(ixO^S,rho_) = .true.

    if(mhd_energy) then
      if(primitive) then
        where(w(ixO^S,p_) < small_pressure) flag(ixO^S,e_) = .true.
      else
        if(B0field) then
          Ba(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))+block%B0(ixO^S,1:ndir,b0i)
        else
          Ba(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))
        end if
        inv_rho(ixO^S) = 1d0/w(ixO^S,rho_)
        b2(ixO^S)=sum(Ba(ixO^S,:)**2,dim=ndim+1)
        tmp(ixO^S)=sqrt(b2(ixO^S))
        where(tmp(ixO^S)>smalldouble)
          tmp(ixO^S)=1.d0/tmp(ixO^S)
        else where
          tmp(ixO^S)=0.d0
        end where
        do idir=1,ndir
          b(ixO^S,idir)=Ba(ixO^S,idir)*tmp(ixO^S)
        end do
        tmp(ixO^S)=sum(b(ixO^S,:)*w(ixO^S,mom(:)),dim=ndim+1)
        ! Va^2/c^2
        b2(ixO^S)=b2(ixO^S)*inv_rho(ixO^S)*inv_squared_c
        ! equation (15)
        gamma2(ixO^S)=1.d0/(1.d0+b2(ixO^S))
        ! Convert momentum to velocity
        do idir = 1, ndir
           v(ixO^S,idir) = gamma2*(w(ixO^S, mom(idir))+b2*b(ixO^S,idir)*tmp)*inv_rho
        end do
        ! E=Bxv
        b=0.d0
        do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
          if(lvc(idir,jdir,kdir)==1)then
            b(ixO^S,idir)=b(ixO^S,idir)+Ba(ixO^S,jdir)*v(ixO^S,kdir)
          else if(lvc(idir,jdir,kdir)==-1)then
            b(ixO^S,idir)=b(ixO^S,idir)-Ba(ixO^S,jdir)*v(ixO^S,kdir)
          end if
        end do; end do; end do
        ! Calculate internal e = e-eK-eB-eE
        tmp(ixO^S)=w(ixO^S,e_)&
                   -half*(sum(v(ixO^S,:)**2,dim=ndim+1)*w(ixO^S,rho_)&
                   +sum(w(ixO^S,mag(:))**2,dim=ndim+1)&
                   +sum(b(ixO^S,:)**2,dim=ndim+1)*inv_squared_c)
        where(tmp(ixO^S) < small_e) flag(ixO^S,e_) = .true.
      end if
    end if

  end subroutine mhd_check_w_semirelati

  subroutine mhd_check_w_origin(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision :: tmp(ixI^S)
    logical, intent(inout) :: flag(ixI^S,1:nw)

    flag=.false.
    if(has_equi_rho0) then
      tmp(ixO^S) = w(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,0)
    else  
      tmp(ixO^S) = w(ixO^S,rho_) 
    endif
    where(tmp(ixO^S) < small_density) flag(ixO^S,rho_) = .true.

    if(mhd_energy) then
      if(primitive) then
        tmp(ixO^S) = w(ixO^S,e_)
        if(has_equi_pe0) then
          tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe0_,0)
        endif
        where(tmp(ixO^S) < small_pressure) flag(ixO^S,e_) = .true.
      else
        if(mhd_internal_e) then
          tmp(ixO^S)=w(ixO^S,e_)
          if(has_equi_pe0) then
            tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe0_,0)*inv_gamma_1
          endif
          where(tmp(ixO^S) < small_e) flag(ixO^S,e_) = .true.
        else
          tmp(ixO^S)=w(ixO^S,e_)-&
              mhd_kin_en(w,ixI^L,ixO^L)-mhd_mag_en(w,ixI^L,ixO^L)
          if(has_equi_pe0) then
            tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe0_,0)*inv_gamma_1
          endif
          where(tmp(ixO^S) < small_e) flag(ixO^S,e_) = .true.
        end if
      end if
    end if

  end subroutine mhd_check_w_origin

  subroutine mhd_check_w_hde(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision :: tmp(ixI^S)
    logical, intent(inout) :: flag(ixI^S,1:nw)

    flag=.false.
    where(w(ixO^S,rho_) < small_density) flag(ixO^S,rho_) = .true.

    if(mhd_energy) then
      if(primitive) then
        where(w(ixO^S,e_) < small_pressure) flag(ixO^S,e_) = .true.
      else
        tmp(ixO^S)=w(ixO^S,e_)-mhd_kin_en(w,ixI^L,ixO^L)
        where(tmp(ixO^S) < small_e) flag(ixO^S,e_) = .true.
      end if
    end if

  end subroutine mhd_check_w_hde

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved_origin(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision :: inv_gamma2(ixO^S)
    integer                         :: idir

    !if (fix_small_values) then
    !  call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'mhd_to_conserved')
    !end if

    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(mhd_energy) then
      w(ixO^S,e_)=w(ixO^S,p_)*inv_gamma_1&
                 +half*sum(w(ixO^S,mom(:))**2,dim=ndim+1)*w(ixO^S,rho_)&
                 +mhd_mag_en(w, ixI^L, ixO^L)
      if(mhd_solve_eaux) w(ixO^S,eaux_)=w(ixO^S,paux_)*inv_gamma_1
    end if

    if(mhd_boris_simplification) then
      inv_gamma2=1.d0+sum(w(ixO^S,mag(:))**2,dim=ndim+1)/w(ixO^S,rho_)*inv_squared_c
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixO^S, mom(idir)) = inv_gamma2*w(ixO^S,rho_)*w(ixO^S, mom(idir))
      end do
    else
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixO^S, mom(idir)) = w(ixO^S,rho_)*w(ixO^S, mom(idir))
      end do
    end if
  end subroutine mhd_to_conserved_origin

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved_hde(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer                         :: idir

    !if (fix_small_values) then
    !  call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'mhd_to_conserved')
    !end if

    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(mhd_energy) then
      w(ixO^S,e_)=w(ixO^S,p_)*inv_gamma_1&
                 +half*sum(w(ixO^S,mom(:))**2,dim=ndim+1)*w(ixO^S,rho_)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
      w(ixO^S, mom(idir)) = w(ixO^S,rho_)*w(ixO^S, mom(idir))
    end do
  end subroutine mhd_to_conserved_hde

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved_inte(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision :: inv_gamma2(ixO^S)
    integer                         :: idir

    !if (fix_small_values) then
    !  call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'mhd_to_conserved')
    !end if

    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(mhd_energy) then
      w(ixO^S,e_)=w(ixO^S,p_)*inv_gamma_1
    end if

    if(mhd_boris_simplification) then
      inv_gamma2=1.d0+sum(w(ixO^S,mag(:))**2,dim=ndim+1)/w(ixO^S,rho_)*inv_squared_c
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixO^S, mom(idir)) = inv_gamma2*w(ixO^S,rho_)*w(ixO^S, mom(idir))
      end do
    else
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixO^S, mom(idir)) = w(ixO^S,rho_)*w(ixO^S, mom(idir))
      end do
    end if
  end subroutine mhd_to_conserved_inte

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved_split_rho(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: idir
    double precision                :: rho(ixI^S)

    !if (fix_small_values) then
    !  call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'mhd_to_conserved')
    !end if

    rho(ixO^S) = w(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,b0i)
    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(mhd_energy) then
      if(mhd_internal_e) then
        w(ixO^S,e_)=w(ixO^S,p_)*inv_gamma_1
      else
        w(ixO^S,e_)=w(ixO^S,p_)*inv_gamma_1&
                   +half*sum(w(ixO^S,mom(:))**2,dim=ndim+1)*rho(ixO^S)&
                   +mhd_mag_en(w, ixI^L, ixO^L)
        if(mhd_solve_eaux) w(ixO^S,eaux_)=w(ixO^S,paux_)*inv_gamma_1
      end if
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = rho(ixO^S) * w(ixO^S, mom(idir))
    end do
  end subroutine mhd_to_conserved_split_rho

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved_semirelati(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision :: E(ixO^S,1:ndir), B(ixO^S,1:ndir), S(ixO^S,1:ndir),tmp(ixO^S), b2(ixO^S)
    integer                         :: idir, jdir, kdir

    !if (fix_small_values) then
    !  call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'mhd_to_conserved')
    !end if

    if(B0field) then
      B(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))+block%B0(ixO^S,1:ndir,b0i)
    else
      B(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))
    end if
    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(mhd_energy) then
      E=0.d0
      do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
        if(lvc(idir,jdir,kdir)==1)then
          E(ixO^S,idir)=E(ixO^S,idir)+B(ixO^S,jdir)*w(ixO^S,mom(kdir))
        else if(lvc(idir,jdir,kdir)==-1)then
          E(ixO^S,idir)=E(ixO^S,idir)-B(ixO^S,jdir)*w(ixO^S,mom(kdir))
        end if
      end do; end do; end do
      ! equation (9)
      w(ixO^S,e_)=w(ixO^S,p_)*inv_gamma_1&
                 +half*(sum(w(ixO^S,mom(:))**2,dim=ndim+1)*w(ixO^S,rho_)&
                 +sum(w(ixO^S,mag(:))**2,dim=ndim+1)&
                 +sum(E(ixO^S,:)**2,dim=ndim+1)*inv_squared_c)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S,rho_) * w(ixO^S, mom(idir))
    end do
    !b2(ixO^S)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    !tmp(ixO^S)=sqrt(b2(ixO^S))
    !where(tmp(ixO^S)>smalldouble)
    !  tmp(ixO^S)=1.d0/tmp(ixO^S)
    !else where
    !  tmp(ixO^S)=0.d0
    !end where
    !do idir=1,ndir
    !  b(ixO^S,idir)=w(ixO^S,mag(idir))*tmp(ixO^S)
    !end do
    !tmp(ixO^S)=sum(b(ixO^S,:)*w(ixO^S,mom(:)),dim=ndim+1)
    !do idir = 1, ndir
    !   w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))+b2(ixO^S)/w(ixO^S,rho_)*inv_squared_c*(w(ixO^S, mom(idir))-b(ixO^S,idir)*tmp(ixO^S))
    !end do
    ! equation (5) Poynting vector
    S=0.d0
    do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
      if(lvc(idir,jdir,kdir)==1)then
        S(ixO^S,idir)=S(ixO^S,idir)+E(ixO^S,jdir)*B(ixO^S,kdir)
      else if(lvc(idir,jdir,kdir)==-1)then
        S(ixO^S,idir)=S(ixO^S,idir)-E(ixO^S,jdir)*B(ixO^S,kdir)
      end if
    end do; end do; end do
    ! equation (9)
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))+S(ixO^S,idir)*inv_squared_c
    end do
  end subroutine mhd_to_conserved_semirelati

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive_origin(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision                :: inv_rho(ixO^S), gamma2(ixO^S)
    integer                         :: idir

    inv_rho(ixO^S) = 1d0/w(ixO^S,rho_) 

    ! Calculate pressure = (gamma-1) * (e-ek-eb)
    if(mhd_energy) then
      w(ixO^S,p_)=gamma_1*(w(ixO^S,e_)&
                  -mhd_kin_en(w,ixI^L,ixO^L,inv_rho)&
                  -mhd_mag_en(w,ixI^L,ixO^L))
      if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,eaux_)*gamma_1
    end if

    if(mhd_boris_simplification) then
      gamma2=1.d0/(1.d0+sum(w(ixO^S,mag(:))**2,dim=ndim+1)*inv_rho*inv_squared_c)
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))*inv_rho*gamma2
      end do
    else
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))*inv_rho
      end do
    end if

    if (fix_small_values) then
      call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'mhd_to_primitive_origin')
    end if

  end subroutine mhd_to_primitive_origin

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive_hde(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision                :: inv_rho(ixO^S)
    integer                         :: idir

    inv_rho(ixO^S) = 1d0/w(ixO^S,rho_) 

    ! Calculate pressure = (gamma-1) * (e-ek)
    if(mhd_energy) then
      w(ixO^S,p_)=gamma_1*(w(ixO^S,e_)-mhd_kin_en(w,ixI^L,ixO^L,inv_rho))
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))*inv_rho
    end do

    if (fix_small_values) then
      call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'mhd_to_primitive_hde')
    end if

  end subroutine mhd_to_primitive_hde

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive_inte(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision                :: inv_rho(ixO^S), gamma2(ixO^S)
    integer                         :: idir

    inv_rho(ixO^S) = 1d0/w(ixO^S,rho_) 

    ! Calculate pressure = (gamma-1) * (e-ek-eb)
    if(mhd_energy) then
      w(ixO^S,p_)=w(ixO^S,e_)*gamma_1
    end if

    if(mhd_boris_simplification) then
      gamma2=1.d0/(1.d0+sum(w(ixO^S,mag(:))**2,dim=ndim+1)*inv_rho*inv_squared_c)
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))*inv_rho*gamma2
      end do
    else
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))*inv_rho
      end do
    end if

    if (fix_small_values) then
      call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'mhd_to_primitive_inte')
    end if

  end subroutine mhd_to_primitive_inte

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive_split_rho(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision                :: inv_rho(ixO^S)
    integer                         :: idir

    inv_rho(ixO^S) = 1d0/(w(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,b0i))

    ! Calculate pressure = (gamma-1) * (e-ek-eb)
    if(mhd_energy) then
      if(mhd_internal_e) then
        w(ixO^S,p_)=w(ixO^S,e_)*gamma_1
      else
        w(ixO^S,p_)=gamma_1*(w(ixO^S,e_)&
                    -mhd_kin_en(w,ixI^L,ixO^L,inv_rho)&
                    -mhd_mag_en(w,ixI^L,ixO^L))
        if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,eaux_)*gamma_1
      end if
    end if
    
    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))*inv_rho
    end do

    if (fix_small_values) then
      call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'mhd_to_primitive_split_rho')
    end if

  end subroutine mhd_to_primitive_split_rho

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive_semirelati(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision                :: inv_rho(ixO^S)
    double precision :: b(ixO^S,1:ndir), Ba(ixO^S,1:ndir),tmp(ixO^S), b2(ixO^S), gamma2(ixO^S)
    integer                         :: idir, jdir, kdir

    if(B0field) then
      Ba(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))+block%B0(ixO^S,1:ndir,b0i)
    else
      Ba(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))
    end if

    inv_rho(ixO^S) = 1d0/w(ixO^S,rho_)

    b2(ixO^S)=sum(Ba(ixO^S,:)**2,dim=ndim+1)
    tmp(ixO^S)=sqrt(b2(ixO^S))
    where(tmp(ixO^S)>smalldouble)
      tmp(ixO^S)=1.d0/tmp(ixO^S)
    else where
      tmp(ixO^S)=0.d0
    end where
    do idir=1,ndir
      b(ixO^S,idir)=Ba(ixO^S,idir)*tmp(ixO^S)
    end do
    tmp(ixO^S)=sum(b(ixO^S,:)*w(ixO^S,mom(:)),dim=ndim+1)

    ! Va^2/c^2
    b2(ixO^S)=b2(ixO^S)*inv_rho(ixO^S)*inv_squared_c
    ! equation (15)
    gamma2(ixO^S)=1.d0/(1.d0+b2(ixO^S))
    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = gamma2*(w(ixO^S, mom(idir))+b2*b(ixO^S,idir)*tmp)*inv_rho
    end do

    if(mhd_energy) then
      ! E=Bxv
      b=0.d0
      do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
        if(lvc(idir,jdir,kdir)==1)then
          b(ixO^S,idir)=b(ixO^S,idir)+Ba(ixO^S,jdir)*w(ixO^S,mom(kdir))
        else if(lvc(idir,jdir,kdir)==-1)then
          b(ixO^S,idir)=b(ixO^S,idir)-Ba(ixO^S,jdir)*w(ixO^S,mom(kdir))
        end if
      end do; end do; end do
      ! Calculate pressure = (gamma-1) * (e-eK-eB-eE)
      w(ixO^S,p_)=gamma_1*(w(ixO^S,e_)&
                 -half*(sum(w(ixO^S,mom(:))**2,dim=ndim+1)*w(ixO^S,rho_)&
                 +sum(w(ixO^S,mag(:))**2,dim=ndim+1)&
                 +sum(b(ixO^S,:)**2,dim=ndim+1)*inv_squared_c))
    end if

    if (fix_small_values) then
      call mhd_handle_small_values_semirelati(.true., w, x, ixI^L, ixO^L, 'mhd_to_primitive_semirelati')
    end if

  end subroutine mhd_to_primitive_semirelati

  !> Transform internal energy to total energy
  subroutine mhd_ei_to_e(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    if(mhd_internal_e) return
 
    ! Calculate total energy from internal, kinetic and magnetic energy
    w(ixO^S,e_)=w(ixO^S,e_)&
               +mhd_kin_en(w,ixI^L,ixO^L)&
               +mhd_mag_en(w,ixI^L,ixO^L)

  end subroutine mhd_ei_to_e

  !> Transform internal energy to hydrodynamic energy
  subroutine mhd_ei_to_e_hde(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    ! Calculate hydrodynamic energy from internal and kinetic
    w(ixO^S,e_)=w(ixO^S,e_)+mhd_kin_en(w,ixI^L,ixO^L)

  end subroutine mhd_ei_to_e_hde

  !> Transform internal energy to total energy and velocity to momentum
  subroutine mhd_ei_to_e_semirelati(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    w(ixO^S,p_)=w(ixO^S,e_)*gamma_1
    call mhd_to_conserved_semirelati(ixI^L,ixO^L,w,x)

  end subroutine mhd_ei_to_e_semirelati

  !> Transform total energy to internal energy
  subroutine mhd_e_to_ei(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    if(mhd_internal_e) return

    ! Calculate ei = e - ek - eb
    w(ixO^S,e_)=w(ixO^S,e_)&
                -mhd_kin_en(w,ixI^L,ixO^L)&
                -mhd_mag_en(w,ixI^L,ixO^L)

    if(fix_small_values) then
      call mhd_handle_small_ei(w,x,ixI^L,ixO^L,e_,'mhd_e_to_ei')
    end if

  end subroutine mhd_e_to_ei

  !> Transform hydrodynamic energy to internal energy
  subroutine mhd_e_to_ei_hde(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    ! Calculate ei = e - ek
    w(ixO^S,e_)=w(ixO^S,e_)-mhd_kin_en(w,ixI^L,ixO^L)

    if(fix_small_values) then
      call mhd_handle_small_ei(w,x,ixI^L,ixO^L,e_,'mhd_e_to_ei_hde')
    end if

  end subroutine mhd_e_to_ei_hde

  !> Transform total energy to internal energy and momentum to velocity
  subroutine mhd_e_to_ei_semirelati(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    call mhd_to_primitive_semirelati(ixI^L,ixO^L,w,x)
    w(ixO^S,e_)=w(ixO^S,p_)*inv_gamma_1

    if(fix_small_values) then
      call mhd_handle_small_ei(w,x,ixI^L,ixO^L,e_,'mhd_e_to_ei_semirelati')
    end if

  end subroutine mhd_e_to_ei_semirelati

  !> Update eaux and transform internal energy to total energy
  subroutine mhd_ei_to_e_aux(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
 
    w(ixO^S,eaux_)=w(ixO^S,e_)
    ! Calculate total energy from internal, kinetic and magnetic energy
    w(ixO^S,e_)=w(ixO^S,e_)&
               +mhd_kin_en(w,ixI^L,ixO^L)&
               +mhd_mag_en(w,ixI^L,ixO^L)

  end subroutine mhd_ei_to_e_aux

  !> Transform total energy to internal energy via eaux as internal energy
  subroutine mhd_e_to_ei_aux(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    w(ixO^S,e_)=w(ixO^S,eaux_)

  end subroutine mhd_e_to_ei_aux

  subroutine mhd_energy_synchro(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth1(ixI^S),pth2(ixI^S),alfa(ixI^S),beta(ixI^S)
    double precision, parameter :: beta_low=0.005d0,beta_high=0.05d0

!    double precision :: vtot(ixI^S),cs2(ixI^S),mach(ixI^S)
!    double precision, parameter :: mach_low=20.d0,mach_high=200.d0

    ! get magnetic energy
    alfa(ixO^S)=mhd_mag_en(w,ixI^L,ixO^L)
    pth1(ixO^S)=gamma_1*(w(ixO^S,e_)-mhd_kin_en(w,ixI^L,ixO^L)-alfa(ixO^S))
    pth2(ixO^S)=w(ixO^S,eaux_)*gamma_1
    ! get plasma beta
    beta(ixO^S)=min(pth1(ixO^S),pth2(ixO^S))/alfa(ixO^S)

    ! whether Mach number should be another criterion ?
!    vtot(ixO^S)=sum(w(ixO^S,mom(:))**2,dim=ndim+1)
!    call mhd_get_csound2(w,x,ixI^L,ixO^L,cs2)
!    mach(ixO^S)=sqrt(vtot(ixO^S)/cs2(ixO^S))/w(ixO^S,rho_)
    where(beta(ixO^S) .ge. beta_high)
!    where(beta(ixO^S) .ge. beta_high .and. mach(ixO^S) .le. mach_low)
      w(ixO^S,eaux_)=pth1(ixO^S)*inv_gamma_1
    else where(beta(ixO^S) .le. beta_low)
!    else where(beta(ixO^S) .le. beta_low .or. mach(ixO^S) .ge. mach_high)
      w(ixO^S,e_)=w(ixO^S,e_)-pth1(ixO^S)*inv_gamma_1+w(ixO^S,eaux_)
    else where
      alfa(ixO^S)=dlog(beta(ixO^S)/beta_low)/dlog(beta_high/beta_low)
!      alfa(ixO^S)=min(dlog(beta(ixO^S)/beta_low)/dlog(beta_high/beta_low),
!                      dlog(mach_high(ixO^S)/mach(ixO^S))/dlog(mach_high/mach_low))
      w(ixO^S,eaux_)=(pth2(ixO^S)*(one-alfa(ixO^S))&
                     +pth1(ixO^S)*alfa(ixO^S))*inv_gamma_1
      w(ixO^S,e_)=w(ixO^S,e_)-pth1(ixO^S)*inv_gamma_1+w(ixO^S,eaux_)
    end where
  end subroutine mhd_energy_synchro

  subroutine mhd_handle_small_values_semirelati(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    double precision :: pressure(ixI^S), inv_rho(ixI^S), b2(ixI^S), tmp(ixI^S), gamma2(ixI^S)
    double precision :: b(ixI^S,1:ndir), v(ixI^S,1:ndir), Ba(ixI^S,1:ndir)
    integer :: idir, jdir, kdir, ix^D
    logical :: flag(ixI^S,1:nw)

    if(small_values_method == "ignore") return

    !call mhd_check_w_semirelati(primitive, ixI^L, ixO^L, w, flag)

    flag=.false.
    where(w(ixO^S,rho_) < small_density) flag(ixO^S,rho_) = .true.

    if(mhd_energy) then
      if(primitive) then
        where(w(ixO^S,p_) < small_pressure) flag(ixO^S,e_) = .true.
      else
        if(B0field) then
          Ba(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,b0i)
        else
          Ba(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
        end if
        inv_rho(ixI^S) = 1d0/w(ixI^S,rho_)
        b2(ixI^S)=sum(Ba(ixI^S,:)**2,dim=ndim+1)
        tmp(ixI^S)=sqrt(b2(ixI^S))
        where(tmp(ixI^S)>smalldouble)
          tmp(ixI^S)=1.d0/tmp(ixI^S)
        else where
          tmp(ixI^S)=0.d0
        end where
        do idir=1,ndir
          b(ixI^S,idir)=Ba(ixI^S,idir)*tmp(ixI^S)
        end do
        tmp(ixI^S)=sum(b(ixI^S,:)*w(ixI^S,mom(:)),dim=ndim+1)
        ! Va^2/c^2
        b2(ixI^S)=b2(ixI^S)*inv_rho(ixI^S)*inv_squared_c
        ! equation (15)
        gamma2(ixI^S)=1.d0/(1.d0+b2(ixI^S))
        ! Convert momentum to velocity
        do idir = 1, ndir
           v(ixI^S,idir) = gamma2*(w(ixI^S, mom(idir))+b2*b(ixI^S,idir)*tmp(ixI^S))*inv_rho(ixI^S)
        end do
        ! E=Bxv
        b=0.d0
        do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
          if(lvc(idir,jdir,kdir)==1)then
            b(ixI^S,idir)=b(ixI^S,idir)+Ba(ixI^S,jdir)*v(ixI^S,kdir)
          else if(lvc(idir,jdir,kdir)==-1)then
            b(ixI^S,idir)=b(ixI^S,idir)-Ba(ixI^S,jdir)*v(ixI^S,kdir)
          end if
        end do; end do; end do
        ! Calculate pressure p = (gamma-1) e-eK-eB-eE
        pressure(ixI^S)=gamma_1*(w(ixI^S,e_)&
                   -half*(sum(v(ixI^S,:)**2,dim=ndim+1)*w(ixI^S,rho_)&
                   +sum(w(ixI^S,mag(:))**2,dim=ndim+1)&
                   +sum(b(ixI^S,:)**2,dim=ndim+1)*inv_squared_c))
        where(pressure(ixI^S) < small_pressure) flag(ixO^S,p_) = .true.
      end if
    end if

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S,rho_)) w(ixO^S,rho_) = small_density

        if(mhd_energy) then
          if(primitive) then
            where(flag(ixO^S,e_)) w(ixO^S,p_) = small_pressure
          else
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
              if(flag(ix^D,e_)) then
                w(ix^D,e_)=small_pressure*inv_gamma_1+half*(sum(v(ix^D,:)**2)*w(ix^D,rho_)&
                           +sum(w(ix^D,mag(:))**2)+sum(b(ix^D,:)**2)*inv_squared_c)
              end if
            {end do\}
          end if
        end if
      case ("average")
        ! do averaging of density
        call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
        if(mhd_energy) then
          if(primitive) then
            call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
          else
            w(ixI^S,e_)=pressure(ixI^S)
            call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
            w(ixI^S,e_)=w(ixI^S,p_)*inv_gamma_1&
                       +half*(sum(v(ixI^S,:)**2,dim=ndim+1)*w(ixI^S,rho_)&
                       +sum(w(ixI^S,mag(:))**2,dim=ndim+1)&
                       +sum(b(ixI^S,:)**2,dim=ndim+1)*inv_squared_c)
          end if
        end if
      case default
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
  end subroutine mhd_handle_small_values_semirelati

  subroutine mhd_handle_small_values_origin(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixI^S,1:nw)
    double precision :: tmp2(ixI^S)

    if(small_values_method == "ignore") return

    call phys_check_w(primitive, ixI^L, ixO^L, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_rho0) then
          where(flag(ixO^S,rho_)) w(ixO^S,rho_) = &
                    small_density-block%equi_vars(ixO^S,equi_rho0_,0)
        else
          where(flag(ixO^S,rho_)) w(ixO^S,rho_) = small_density
        endif
        do idir = 1, ndir
          if(small_values_fix_iw(mom(idir))) then
            where(flag(ixO^S,rho_)) w(ixO^S, mom(idir)) = 0.0d0
          end if
        end do

        if(mhd_energy) then
          if(primitive) then
           if(has_equi_pe0) then 
            tmp2(ixO^S) = small_pressure - &
              block%equi_vars(ixO^S,equi_pe0_,0)
           else
            tmp2(ixO^S) = small_pressure
           endif  
           where(flag(ixO^S,e_)) w(ixO^S,p_) = tmp2(ixO^S)
          else
            ! conserved
            if(has_equi_pe0) then 
              tmp2(ixO^S) = small_e - &
                block%equi_vars(ixO^S,equi_pe0_,0)*inv_gamma_1 
            else
              tmp2(ixO^S) = small_e
            endif  
            if(mhd_internal_e) then
              where(flag(ixO^S,e_))
                w(ixO^S,e_)=tmp2(ixO^S)
              end where
            else
              where(flag(ixO^S,e_))
                w(ixO^S,e_) = tmp2(ixO^S)+&
                   mhd_kin_en(w,ixI^L,ixO^L)+&
                   mhd_mag_en(w,ixI^L,ixO^L)
              end where
              if(mhd_solve_eaux) then
                where(flag(ixO^S,e_))
                  w(ixO^S,eaux_)=tmp2(ixO^S)
                end where
              end if
            end if
          end if
        end if
      case ("average")
        if(primitive)then
          call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
          if(mhd_energy) then
            call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
          end if
        else
          ! do averaging of density
          call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
          if(mhd_energy) then
             ! do averaging of pressure
            if(mhd_internal_e) then
              w(ixI^S,p_)=w(ixI^S,e_)*gamma_1
            else
              w(ixI^S,p_)=gamma_1*(w(ixI^S,e_)&
                          -mhd_kin_en(w,ixI^L,ixI^L)&
                          -mhd_mag_en(w,ixI^L,ixI^L))
            end if
            call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
             ! convert back
            if(mhd_internal_e) then
              w(ixI^S,p_)=w(ixI^S,e_)*gamma_1
            else
              w(ixI^S,p_)=inv_gamma_1*w(ixI^S,p_)&
                          +mhd_kin_en(w,ixI^L,ixI^L)&
                          +mhd_mag_en(w,ixI^L,ixI^L)
            end if
            ! eaux
            if(mhd_solve_eaux) then
              !w(ixI^S,paux_)=w(ixI^S,eaux_)*gamma_1
              call small_values_average(ixI^L, ixO^L, w, x, flag, paux_)
              !w(ixI^S,eaux_)=w(ixI^S,paux_)*inv_gamma_1
            end if
          end if
        endif
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek-eb)
          if(mhd_energy) then
            if(mhd_internal_e) then
              w(ixO^S,p_)=w(ixO^S,e_)*gamma_1
            else
              w(ixO^S,p_)=gamma_1*(w(ixO^S,e_)&
                          -mhd_kin_en(w,ixI^L,ixO^L)&
                          -mhd_mag_en(w,ixI^L,ixO^L))
              if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,eaux_)*gamma_1
            end if
          end if
          ! Convert momentum to velocity
          if(has_equi_rho0) then
            tmp2(ixO^S) = w(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,0)
          else  
            tmp2(ixO^S) = w(ixO^S,rho_) 
          endif
          do idir = 1, ndir
             w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))/tmp2(ixO^S)
          end do
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
  end subroutine mhd_handle_small_values_origin

  subroutine mhd_handle_small_values_hde(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixI^S,1:nw)
    double precision :: tmp2(ixI^S)

    if(small_values_method == "ignore") return

    call phys_check_w(primitive, ixI^L, ixO^L, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S,rho_)) w(ixO^S,rho_) = small_density
        do idir = 1, ndir
          if(small_values_fix_iw(mom(idir))) then
            where(flag(ixO^S,rho_)) w(ixO^S, mom(idir)) = 0.0d0
          end if
        end do

        if(mhd_energy) then
          where(flag(ixO^S,e_))
            w(ixO^S,e_) = small_e+mhd_kin_en(w,ixI^L,ixO^L)
          end where
        end if
      case ("average")
        ! do averaging of density
        call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
        if(mhd_energy) then
          if(primitive) then
            call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
          else
            ! do averaging of pressure
            w(ixI^S,p_)=gamma_1*(w(ixI^S,e_)-mhd_kin_en(w,ixI^L,ixI^L))
            call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
            ! convert back
            w(ixI^S,p_)=inv_gamma_1*w(ixI^S,p_)+mhd_kin_en(w,ixI^L,ixI^L)
          end if
        end if
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek)
          if(mhd_energy) then
            w(ixO^S,p_)=gamma_1*(w(ixO^S,e_)-mhd_kin_en(w,ixI^L,ixO^L))
          end if
          ! Convert momentum to velocity
          do idir = 1, ndir
             w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))/w(ixO^S,rho_)
          end do
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if

  end subroutine mhd_handle_small_values_hde

  !> Calculate v vector
  subroutine mhd_get_v_origin(w,x,ixI^L,ixO^L,v)
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

  end subroutine mhd_get_v_origin

  !> Calculate v vector
  subroutine mhd_get_v_boris(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,ndir)

    double precision              :: rho(ixI^S), gamma2(ixO^S)
    integer :: idir

    call mhd_get_rho(w,x,ixI^L,ixO^L,rho)

    rho(ixO^S)=1.d0/rho(ixO^S)
    gamma2=1.d0/(1.d0+sum(w(ixO^S,mag(:))**2,dim=ndim+1)*rho(ixO^S)*inv_squared_c)
    ! Convert momentum to velocity
    do idir = 1, ndir
       v(ixO^S, idir) = w(ixO^S, mom(idir))*rho(ixO^S)*gamma2
    end do

  end subroutine mhd_get_v_boris

  !> Calculate v component
  subroutine mhd_get_v_idim(w,x,ixI^L,ixO^L,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S)

    double precision              :: rho(ixI^S)

    call mhd_get_rho(w,x,ixI^L,ixO^L,rho)

    if(mhd_boris_simplification) then
      v(ixO^S) = w(ixO^S, mom(idim)) / rho(ixO^S) &
       /(1.d0+sum(w(ixO^S,mag(:))**2,dim=ndim+1)/rho(ixO^S)*inv_squared_c)
    else
      v(ixO^S) = w(ixO^S, mom(idim)) / rho(ixO^S)
    end if

  end subroutine mhd_get_v_idim

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine mhd_get_cmax_origin(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision :: vel(ixI^S)

    call mhd_get_csound(w,x,ixI^L,ixO^L,idim,cmax)
    call mhd_get_v_idim(w,x,ixI^L,ixO^L,idim,vel)

    cmax(ixO^S)=abs(vel(ixO^S))+cmax(ixO^S)

  end subroutine mhd_get_cmax_origin

  !> Calculate cmax_idim for semirelativistic MHD
  subroutine mhd_get_cmax_semirelati(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout):: cmax(ixI^S)
    double precision :: wprim(ixI^S,nw)
    double precision :: csound(ixO^S), AvMinCs2(ixO^S), idim_Alfven_speed2(ixO^S)
    double precision :: inv_rho(ixO^S), Alfven_speed2(ixO^S), gamma2(ixO^S), B(ixO^S,1:ndir)

    if(B0field) then
      B(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))+block%B0(ixO^S,1:ndir,b0i)
    else
      B(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))
    end if
    inv_rho = 1.d0/w(ixO^S,rho_) 

    Alfven_speed2=sum(B(ixO^S,:)**2,dim=ndim+1)*inv_rho
    gamma2 = 1.0d0/(1.d0+Alfven_speed2*inv_squared_c)

    wprim=w
    call mhd_to_primitive(ixI^L,ixO^L,wprim,x)
    cmax(ixO^S)=1.d0-gamma2*wprim(ixO^S,mom(idim))**2*inv_squared_c
    ! equation (69)
    Alfven_speed2=Alfven_speed2*cmax(ixO^S)

    ! squared sound speed
    if(mhd_energy) then
      csound=mhd_gamma*wprim(ixO^S,p_)*inv_rho
    else
      csound=mhd_gamma*mhd_adiab*w(ixO^S,rho_)**gamma_1
    end if

    idim_Alfven_speed2=B(ixO^S,idim)**2*inv_rho

    ! Va_hat^2+a_hat^2 equation (57)
    Alfven_speed2=Alfven_speed2+csound*(1.d0+idim_Alfven_speed2*inv_squared_c)

    AvMinCs2=(gamma2*Alfven_speed2)**2-4.0d0*gamma2*csound*idim_Alfven_speed2*cmax(ixO^S)

    where(AvMinCs2<zero)
       AvMinCs2=zero
    end where

    ! equation (68) fast magnetosonic wave speed
    csound = sqrt(half*(gamma2*Alfven_speed2+sqrt(AvMinCs2)))
    cmax(ixO^S)=gamma2*abs(wprim(ixO^S,mom(idim)))+csound

  end subroutine mhd_get_cmax_semirelati

  subroutine mhd_get_a2max(w,x,ixI^L,ixO^L,a2max)
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
  end subroutine mhd_get_a2max

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine mhd_get_tcutoff(ixI^L,ixO^L,w,x,Tco_local,Tmax_local)
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
    integer :: jxP^L,hxP^L,ixP^L
    logical :: lrlt(ixI^S)

    ! reuse lts as rhoc
    call mhd_get_rho(w,x,ixI^L,ixI^L,lts)
    if(mhd_internal_e) then
      tmp1(ixI^S)=w(ixI^S,e_)*(mhd_gamma-1.d0)
    else
      call phys_get_pthermal(w,x,ixI^L,ixI^L,tmp1)
    end if
    Te(ixI^S)=tmp1(ixI^S)/lts(ixI^S)
    Tco_local=zero
    Tmax_local=maxval(Te(ixO^S))

    {^IFONED
    select case(mhd_trac_type)
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
        call gradient(Te,ixI^L,ixO^L,idims,tmp1)
        gradT(ixO^S,idims)=tmp1(ixO^S)
      end do
      ! B vector
      if(B0field) then
        bunitvec(ixO^S,:)=w(ixO^S,iw_mag(:))+block%B0(ixO^S,:,0)
      else
        bunitvec(ixO^S,:)=w(ixO^S,iw_mag(:))
      end if
      if(mhd_trac_type .gt. 1) then
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
      ixP^L=ixO^L^LADD1;
      ! temperature gradient at cell centers
      do idims=1,ndim
        call gradient(Te,ixI^L,ixP^L,idims,tmp1)
        gradT(ixP^S,idims)=tmp1(ixP^S)
      end do
      ! B vector
      if(B0field) then
        bunitvec(ixP^S,:)=w(ixP^S,iw_mag(:))+block%B0(ixP^S,:,0)
      else
        bunitvec(ixP^S,:)=w(ixP^S,iw_mag(:))
      end if
      tmp1(ixP^S)=dsqrt(sum(bunitvec(ixP^S,:)**2,dim=ndim+1))
      where(tmp1(ixP^S)/=0.d0)
        tmp1(ixP^S)=1.d0/tmp1(ixP^S)
      elsewhere
        tmp1(ixP^S)=bigdouble
      end where
      ! b unit vector: magnetic field direction vector
      do idims=1,ndim
        bunitvec(ixP^S,idims)=bunitvec(ixP^S,idims)*tmp1(ixP^S)
      end do
      ! temperature length scale inversed
      lts(ixP^S)=abs(sum(gradT(ixP^S,1:ndim)*bunitvec(ixP^S,1:ndim),dim=ndim+1))/Te(ixP^S)
      ! fraction of cells size to temperature length scale
      if(slab_uniform) then
        lts(ixP^S)=minval(dxlevel)*lts(ixP^S)
      else
        lts(ixP^S)=minval(block%ds(ixP^S,:),dim=ndim+1)*lts(ixP^S)
      end if
      lts(ixP^S)=max(one, (exp(lts(ixP^S))/ltrc)**ltrp)
  
      altr=zero
      do idims=1,ndim 
        hxO^L=ixO^L-kr(idims,^D);
        jxO^L=ixO^L+kr(idims,^D);
        altr(ixO^S)=altr(ixO^S)+0.25d0*(lts(hxO^S)+two*lts(ixO^S)+lts(jxO^S))*bunitvec(ixO^S,idims)**2
      end do
      block%wextra(ixO^S,Tcoff_)=Te(ixO^S)*altr(ixO^S)**0.4d0
      ! need one ghost layer for thermal conductivity
      {block%wextra(ixOmin^D-1^D%ixP^S,Tcoff_)=block%wextra(ixOmin^D^D%ixP^S,Tcoff_) \}
      {block%wextra(ixOmax^D+1^D%ixP^S,Tcoff_)=block%wextra(ixOmax^D^D%ixP^S,Tcoff_) \}
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

    double precision :: csound(ixI^S,ndim),tmp(ixI^S)
    integer :: jxC^L, ixC^L, ixA^L, id, ix^D

    Hspeed=0.d0
    ixA^L=ixO^L^LADD1;
    do id=1,ndim
      call mhd_get_csound_prim(wprim,x,ixI^L,ixA^L,id,tmp)
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
      umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)+wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)
      call mhd_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)
      dmean(ixO^S)=(tmp1(ixO^S)*csoundL(ixO^S)**2+tmp2(ixO^S)*csoundR(ixO^S)**2)*tmp3(ixO^S)+&
       0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2*&
       (wRp(ixO^S,mom(idim))-wLp(ixO^S,mom(idim)))**2
      dmean(ixO^S)=sqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,1)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S,1)=umean(ixO^S)+dmean(ixO^S)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=abs(umean(ixO^S))+dmean(ixO^S)
      end if
    case (2)
      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
      tmp1(ixO^S)=wmean(ixO^S,mom(idim))/wmean(ixO^S,rho_)
      call mhd_get_csound(wmean,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
        cmax(ixO^S,1)=max(tmp1(ixO^S)+csoundR(ixO^S),zero)
        cmin(ixO^S,1)=min(tmp1(ixO^S)-csoundR(ixO^S),zero)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=abs(tmp1(ixO^S))+csoundR(ixO^S)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call mhd_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)
      csoundL(ixO^S)=max(csoundL(ixO^S),csoundR(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,1)=min(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))-csoundL(ixO^S)
        cmax(ixO^S,1)=max(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=max(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
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

    ! Miyoshi 2005 JCP 208, 315 equation (67)
    call mhd_get_csound_semirelati(wLp,x,ixI^L,ixO^L,idim,csoundL,gamma2L)
    call mhd_get_csound_semirelati(wRp,x,ixI^L,ixO^L,idim,csoundR,gamma2R)
    csoundL(ixO^S)=max(csoundL(ixO^S),csoundR(ixO^S))
    if(present(cmin)) then
      cmin(ixO^S,1)=min(gamma2L*wLp(ixO^S,mom(idim)),gamma2R*wRp(ixO^S,mom(idim)))-csoundL(ixO^S)
      cmax(ixO^S,1)=max(gamma2L*wLp(ixO^S,mom(idim)),gamma2R*wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
    else
      cmax(ixO^S,1)=max(gamma2L*wLp(ixO^S,mom(idim)),gamma2R*wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
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

    double precision :: wmean(ixI^S,nw)
    double precision, dimension(ixI^S) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix^D
    double precision :: rho(ixI^S)

    select case (boundspeed) 
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      tmp1(ixO^S)=sqrt(wLp(ixO^S,rho_)+block%equi_vars(ixO^S,equi_rho0_,b0i))
      tmp2(ixO^S)=sqrt(wRp(ixO^S,rho_)+block%equi_vars(ixO^S,equi_rho0_,b0i))
      tmp3(ixO^S)=1.d0/(tmp1(ixO^S)+tmp2(ixO^S))
      umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)+wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)
      call mhd_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)
      dmean(ixO^S)=(tmp1(ixO^S)*csoundL(ixO^S)**2+tmp2(ixO^S)*csoundR(ixO^S)**2)*tmp3(ixO^S)+&
       0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2*&
       (wRp(ixO^S,mom(idim))-wLp(ixO^S,mom(idim)))**2
      dmean(ixO^S)=sqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,1)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S,1)=umean(ixO^S)+dmean(ixO^S)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=abs(umean(ixO^S))+dmean(ixO^S)
      end if
    case (2)
      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
      tmp1(ixO^S)=wmean(ixO^S,mom(idim))/(wmean(ixO^S,rho_)+block%equi_vars(ixO^S,equi_rho0_,b0i))
      call mhd_get_csound(wmean,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
        cmax(ixO^S,1)=max(tmp1(ixO^S)+csoundR(ixO^S),zero)
        cmin(ixO^S,1)=min(tmp1(ixO^S)-csoundR(ixO^S),zero)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=abs(tmp1(ixO^S))+csoundR(ixO^S)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call mhd_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)
      csoundL(ixO^S)=max(csoundL(ixO^S),csoundR(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,1)=min(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))-csoundL(ixO^S)
        cmax(ixO^S,1)=max(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=max(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
      end if
    end select

  end subroutine mhd_get_cbounds_split_rho

  !> prepare velocities for ct methods
  subroutine mhd_get_ct_velocity(vcts,wLp,wRp,ixI^L,ixO^L,idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: cmax(ixI^S)
    double precision, intent(in), optional :: cmin(ixI^S)
    type(ct_velocity), intent(inout):: vcts

    integer                         :: idimE,idimN

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

  end subroutine mhd_get_ct_velocity

  !> Calculate fast magnetosonic wave speed
  subroutine mhd_get_csound(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
    double precision :: inv_rho(ixO^S)

    if(has_equi_rho0) then
      inv_rho(ixO^S) = 1d0/(w(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,b0i))
    else  
      inv_rho(ixO^S) = 1d0/w(ixO^S,rho_) 
    endif

    call mhd_get_csound2(w,x,ixI^L,ixO^L,csound)

    ! store |B|^2 in v
    b2(ixO^S) = mhd_mag_en_all(w,ixI^L,ixO^L)

    cfast2(ixO^S)   = b2(ixO^S) * inv_rho+csound(ixO^S)
    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
         * mhd_mag_i_all(w,ixI^L,ixO^L,idim)**2 * inv_rho

    where(AvMinCs2(ixO^S)<zero)
       AvMinCs2(ixO^S)=zero
    end where

    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

    if (.not. MHD_Hall) then
       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
       if (mhd_boris_simplification) then
          ! equation (88)
          csound(ixO^S) = mhd_gamma_alfven(w, ixI^L,ixO^L) * csound(ixO^S)
       end if
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
            mhd_etah * sqrt(b2(ixO^S))*inv_rho*kmax)
    end if

  end subroutine mhd_get_csound

  !> Calculate fast magnetosonic wave speed
  subroutine mhd_get_csound_prim(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
    double precision :: inv_rho(ixO^S)
    double precision :: tmp(ixI^S)

    call mhd_get_rho(w,x,ixI^L,ixO^L,tmp)
    inv_rho(ixO^S) = 1d0/tmp(ixO^S) 


    if(mhd_energy) then
      if(has_equi_pe0) then
        csound(ixO^S) = w(ixO^S,e_) + block%equi_vars(ixO^S,equi_pe0_,b0i)
      else
        csound(ixO^S) = w(ixO^S,e_) 
      endif
      csound(ixO^S)=mhd_gamma*csound(ixO^S)*inv_rho
    else
      csound(ixO^S)=mhd_gamma*mhd_adiab*tmp(ixO^S)**gamma_1
    end if

    ! store |B|^2 in v
    b2(ixO^S)        = mhd_mag_en_all(w,ixI^L,ixO^L)
    cfast2(ixO^S)   = b2(ixO^S) * inv_rho+csound(ixO^S)
    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
         * mhd_mag_i_all(w,ixI^L,ixO^L,idim)**2 * inv_rho

    where(AvMinCs2(ixO^S)<zero)
       AvMinCs2(ixO^S)=zero
    end where

    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

    if (.not. MHD_Hall) then
       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
       if (mhd_boris_simplification) then
          ! equation (88)
          csound(ixO^S) = mhd_gamma_alfven(w, ixI^L,ixO^L) * csound(ixO^S)
       end if
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
            mhd_etah * sqrt(b2(ixO^S))*inv_rho*kmax)
    end if

  end subroutine mhd_get_csound_prim

  !> Calculate cmax_idim for semirelativistic MHD
  subroutine mhd_get_csound_semirelati(w,x,ixI^L,ixO^L,idim,csound,gamma2)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! here w is primitive variables
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixO^S), gamma2(ixO^S)
    double precision :: AvMinCs2(ixO^S), kmax
    double precision :: inv_rho(ixO^S), Alfven_speed2(ixO^S), idim_Alfven_speed2(ixO^S),B(ixO^S,1:ndir)

    if(B0field) then
      B(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))+block%B0(ixO^S,1:ndir,b0i)
    else
      B(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))
    end if

    inv_rho = 1.d0/w(ixO^S,rho_) 

    Alfven_speed2=sum(B(ixO^S,:)**2,dim=ndim+1)*inv_rho
    gamma2 = 1.0d0/(1.d0+Alfven_speed2*inv_squared_c)

    AvMinCs2=1.d0-gamma2*w(ixO^S,mom(idim))**2*inv_squared_c
    ! equatoin (69)
    Alfven_speed2=Alfven_speed2*AvMinCs2

    ! squared sound speed
    if(mhd_energy) then
      csound(ixO^S)=mhd_gamma*w(ixO^S,p_)*inv_rho
    else
      csound(ixO^S)=mhd_gamma*mhd_adiab*w(ixO^S,rho_)**gamma_1
    end if

    idim_Alfven_speed2=B(ixO^S,idim)**2*inv_rho

    ! Va_hat^2+a_hat^2 equation (57)
    Alfven_speed2=Alfven_speed2+csound(ixO^S)*(1.d0+idim_Alfven_speed2*inv_squared_c)

    AvMinCs2=(gamma2*Alfven_speed2)**2-4.0d0*gamma2*csound(ixO^S)*idim_Alfven_speed2*AvMinCs2

    where(AvMinCs2<zero)
       AvMinCs2=zero
    end where

    ! equation (68) fast magnetosonic speed
    csound(ixO^S) = sqrt(half*(gamma2*Alfven_speed2+sqrt(AvMinCs2)))

  end subroutine mhd_get_csound_semirelati

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine mhd_get_pthermal_origin(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)
    integer                      :: iw, ix^D

    if(mhd_energy) then
      if(mhd_internal_e) then
        pth(ixO^S)=gamma_1*w(ixO^S,e_)
      else
        pth(ixO^S)=gamma_1*(w(ixO^S,e_)&
           - mhd_kin_en(w,ixI^L,ixO^L)&
           - mhd_mag_en(w,ixI^L,ixO^L))
      end if
      if(has_equi_pe0) then
        pth(ixO^S) = pth(ixO^S) + block%equi_vars(ixO^S,equi_pe0_,b0i)
      endif
    else
      call mhd_get_rho(w,x,ixI^L,ixO^L,pth)
      pth(ixO^S)=mhd_adiab*pth(ixO^S)**mhd_gamma
    end if

    if (fix_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(pth(ix^D)<small_pressure) then
            pth(ix^D)=small_pressure
         end if
      {enddo^D&\}
    end if

    if (check_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(pth(ix^D)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                " encountered when call mhd_get_pthermal"
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
      {enddo^D&\}
    end if

  end subroutine mhd_get_pthermal_origin

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine mhd_get_pthermal_semirelati(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)
    integer                      :: iw, ix^D

    double precision :: wprim(ixI^S,nw)

    if(mhd_energy) then
      wprim=w
      call mhd_to_primitive_semirelati(ixI^L,ixO^L,wprim,x)
      pth(ixO^S)=wprim(ixO^S,p_)
    else
      pth(ixO^S)=mhd_adiab*w(ixO^S,rho_)**mhd_gamma
    end if

    if (check_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(pth(ix^D)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                " encountered when call mhd_get_pthermal_semirelati"
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
      {enddo^D&\}
    end if

  end subroutine mhd_get_pthermal_semirelati

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine mhd_get_pthermal_hde(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)
    integer                      :: iw, ix^D

    if(mhd_energy) then
      pth(ixO^S)=gamma_1*(w(ixO^S,e_)-mhd_kin_en(w,ixI^L,ixO^L))
    else
      pth(ixO^S)=mhd_adiab*w(ixO^S,rho_)**mhd_gamma
    end if

    if (fix_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(pth(ix^D)<small_pressure) then
            pth(ix^D)=small_pressure
         end if
      {enddo^D&\}
    end if

    if (check_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(pth(ix^D)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                " encountered when call mhd_get_pthermal_hde"
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
      {enddo^D&\}
    end if

  end subroutine mhd_get_pthermal_hde

  !> Calculate temperature=p/rho when in e_ the internal energy is stored
  subroutine mhd_get_temperature_from_eint(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S) = gamma_1 * w(ixO^S, e_) /w(ixO^S,rho_)
  end subroutine mhd_get_temperature_from_eint

  !> Calculate temperature=p/rho when in e_ the total energy is stored
  !> this does not check the values of mhd_energy and mhd_internal_e, 
  !>  mhd_energy = .true. and mhd_internal_e = .false.
  !> also check small_values is avoided
  subroutine mhd_get_temperature_from_etot(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S)=(gamma_1*(w(ixO^S,e_)&
           - mhd_kin_en(w,ixI^L,ixO^L)&
           - mhd_mag_en(w,ixI^L,ixO^L)))/w(ixO^S,rho_)
  end subroutine mhd_get_temperature_from_etot

  !> Calculate temperature from hydrodynamic energy
  subroutine mhd_get_temperature_from_hde(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S)=gamma_1*(w(ixO^S,e_)&
           - mhd_kin_en(w,ixI^L,ixO^L))/w(ixO^S,rho_)
  end subroutine mhd_get_temperature_from_hde

  subroutine mhd_get_temperature_from_eint_with_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S) = (gamma_1 * w(ixO^S, e_) + block%equi_vars(ixO^S,equi_pe0_,b0i)) /&
                (w(ixO^S,rho_) +block%equi_vars(ixO^S,equi_rho0_,b0i))
  end subroutine mhd_get_temperature_from_eint_with_equi

  subroutine mhd_get_temperature_equi(w,x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S)= block%equi_vars(ixO^S,equi_pe0_,b0i)/block%equi_vars(ixO^S,equi_rho0_,b0i)
  end subroutine mhd_get_temperature_equi

  subroutine mhd_get_rho_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S) = block%equi_vars(ixO^S,equi_rho0_,b0i)
  end subroutine mhd_get_rho_equi

  subroutine mhd_get_pe_equi(w,x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S) = block%equi_vars(ixO^S,equi_pe0_,b0i)
  end subroutine mhd_get_pe_equi

  subroutine mhd_get_temperature_from_etot_with_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S)=(gamma_1*(w(ixO^S,e_)&
           - mhd_kin_en(w,ixI^L,ixO^L)&
           - mhd_mag_en(w,ixI^L,ixO^L)) +  block%equi_vars(ixO^S,equi_pe0_,b0i))&
            /(w(ixO^S,rho_) +block%equi_vars(ixO^S,equi_rho0_,b0i))
            
  end subroutine mhd_get_temperature_from_etot_with_equi

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine mhd_get_csound2(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision    :: rho(ixI^S)
    
    call mhd_get_rho(w,x,ixI^L,ixO^L,rho)
    if(mhd_energy) then
      call mhd_get_pthermal(w,x,ixI^L,ixO^L,csound2)
      csound2(ixO^S)=mhd_gamma*csound2(ixO^S)/rho(ixO^S)
    else
      csound2(ixO^S)=mhd_gamma*mhd_adiab*rho(ixO^S)**gamma_1
    end if
  end subroutine mhd_get_csound2

  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine mhd_get_p_total(w,x,ixI^L,ixO^L,p)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: p(ixI^S)

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,p)

    p(ixO^S) = p(ixO^S) + 0.5d0 * sum(w(ixO^S, mag(:))**2, dim=ndim+1)

  end subroutine mhd_get_p_total

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

    double precision             :: ptotal(ixO^S)
    double precision             :: tmp(ixI^S)
    double precision             :: vHall(ixI^S,1:ndir)
    integer                      :: idirmin, iw, idir, jdir, kdir
    double precision, allocatable, dimension(:^D&,:) :: Jambi, btot
    double precision, allocatable, dimension(:^D&) :: tmp2, tmp3

    ! Get flux of density
    f(ixO^S,rho_)=w(ixO^S,mom(idim))*w(ixO^S,rho_)

    if(mhd_energy) then
      ptotal=w(ixO^S,p_)+0.5d0*sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    else
      ptotal=mhd_adiab*w(ixO^S,rho_)**mhd_gamma+0.5d0*sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    end if

    if (mhd_Hall) then
      call mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
    end if

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixO^S,tracer(iw))=w(ixO^S,mom(idim))*w(ixO^S,tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixO^S,mom(idir))=wC(ixO^S,mom(idir))*w(ixO^S,mom(idim))+ptotal(ixO^S)-&
                            w(ixO^S,mag(idir))*w(ixO^S,mag(idim))
      else
        f(ixO^S,mom(idir))=wC(ixO^S,mom(idir))*w(ixO^S,mom(idim))-&
                            w(ixO^S,mag(idir))*w(ixO^S,mag(idim))
      end if
    end do

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if(mhd_energy) then
      if (mhd_internal_e) then
         f(ixO^S,e_)=w(ixO^S,mom(idim))*w(ixO^S,p_)
         if (mhd_Hall) then
            call mpistop("solve internal energy not implemented for Hall MHD")
         endif
      else
        f(ixO^S,e_)=w(ixO^S,mom(idim))*(wC(ixO^S,e_)+ptotal(ixO^S))&
           -w(ixO^S,mag(idim))*sum(w(ixO^S,mag(:))*w(ixO^S,mom(:)),dim=ndim+1)
        if(mhd_solve_eaux) f(ixO^S,eaux_)=w(ixO^S,mom(idim))*wC(ixO^S,eaux_)
        if(mhd_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
          f(ixO^S,e_) = f(ixO^S,e_) + vHall(ixO^S,idim) * &
             sum(w(ixO^S, mag(:))**2,dim=ndim+1) &
             - w(ixO^S,mag(idim)) * sum(vHall(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1)
        end if
      end if
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mhd_glm) then
           f(ixO^S,mag(idir))=w(ixO^S,psi_)
        else
           f(ixO^S,mag(idir))=zero
        end if
      else
        f(ixO^S,mag(idir))=w(ixO^S,mom(idim))*w(ixO^S,mag(idir))-w(ixO^S,mag(idim))*w(ixO^S,mom(idir))
        if (mhd_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          f(ixO^S,mag(idir)) = f(ixO^S,mag(idir)) &
               - vHall(ixO^S,idir)*w(ixO^S,mag(idim)) &
               + vHall(ixO^S,idim)*w(ixO^S,mag(idir))
        end if
      end if
    end do

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixO^S,psi_)  = cmax_global**2*w(ixO^S,mag(idim))
    end if

    ! Contributions of ambipolar term in explicit scheme
    if(mhd_ambipolar_exp.and. .not.stagger_grid) then
      ! ambipolar electric field
      ! E_ambi=-eta_ambi*JxBxB=-JaxBxB=B^2*Ja-(Ja dot B)*B
      !Ja=eta_ambi*J=J * mhd_eta_ambi/rho**2
      allocate(Jambi(ixI^S,1:3))
      call mhd_get_Jambi(w,x,ixI^L,ixO^L,Jambi)
      allocate(btot(ixO^S,1:3))
      btot(ixO^S,1:3) = w(ixO^S,mag(1:3))
      allocate(tmp2(ixO^S),tmp3(ixO^S))
      !tmp2 = Btot^2
      tmp2(ixO^S) = sum(btot(ixO^S,1:3)**2,dim=ndim+1)
      !tmp3 = J_ambi dot Btot
      tmp3(ixO^S) = sum(Jambi(ixO^S,:)*btot(ixO^S,:),dim=ndim+1)

      select case(idim)
        case(1)
          tmp(ixO^S)=w(ixO^S,mag(3)) *Jambi(ixO^S,2) - w(ixO^S,mag(2)) * Jambi(ixO^S,3)
          f(ixO^S,mag(2))= f(ixO^S,mag(2)) - tmp2(ixO^S) * Jambi(ixO^S,3) + tmp3(ixO^S) * btot(ixO^S,3)
          f(ixO^S,mag(3))= f(ixO^S,mag(3)) + tmp2(ixO^S) * Jambi(ixO^S,2) - tmp3(ixO^S) * btot(ixO^S,2)
        case(2)
          tmp(ixO^S)=w(ixO^S,mag(1)) *Jambi(ixO^S,3) - w(ixO^S,mag(3)) * Jambi(ixO^S,1)
          f(ixO^S,mag(1))= f(ixO^S,mag(1)) + tmp2(ixO^S) * Jambi(ixO^S,3) - tmp3(ixO^S) * btot(ixO^S,3)
          f(ixO^S,mag(3))= f(ixO^S,mag(3)) - tmp2(ixO^S) * Jambi(ixO^S,1) + tmp3(ixO^S) * btot(ixO^S,1)
        case(3)
          tmp(ixO^S)=w(ixO^S,mag(2)) *Jambi(ixO^S,1) - w(ixO^S,mag(1)) * Jambi(ixO^S,2)
          f(ixO^S,mag(1))= f(ixO^S,mag(1)) - tmp2(ixO^S) * Jambi(ixO^S,2) + tmp3(ixO^S) * btot(ixO^S,2)
          f(ixO^S,mag(2))= f(ixO^S,mag(2)) + tmp2(ixO^S) * Jambi(ixO^S,1) - tmp3(ixO^S) * btot(ixO^S,1)
      endselect

      if(mhd_energy .and. .not. mhd_internal_e) then
        f(ixO^S,e_) = f(ixO^S,e_) + tmp2(ixO^S) *  tmp(ixO^S)
      endif

      deallocate(Jambi,btot,tmp2,tmp3)
    endif

  end subroutine mhd_get_flux

  !> Calculate fluxes within ixO^L without any splitting
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

    double precision             :: pgas(ixO^S), ptotal(ixO^S)
    double precision             :: tmp(ixI^S)
    integer                      :: idirmin, iw, idir, jdir, kdir
    double precision, allocatable, dimension(:^D&,:) :: Jambi, btot
    double precision, allocatable, dimension(:^D&) :: tmp2, tmp3

    ! Get flux of density
    f(ixO^S,rho_)=w(ixO^S,mom(idim))*w(ixO^S,rho_)
    ! pgas is time dependent only
    if(mhd_energy) then
      pgas=w(ixO^S,p_)
    else
      pgas(ixO^S)=mhd_adiab*w(ixO^S,rho_)**mhd_gamma
    end if

    ptotal = pgas + 0.5d0*sum(w(ixO^S, mag(:))**2, dim=ndim+1)

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixO^S,tracer(iw))=w(ixO^S,mom(idim))*w(ixO^S,tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixO^S,mom(idir))=wC(ixO^S,mom(idir))*w(ixO^S,mom(idim))+ptotal(ixO^S)-&
                            w(ixO^S,mag(idir))*w(ixO^S,mag(idim))
      else
        f(ixO^S,mom(idir))=wC(ixO^S,mom(idir))*w(ixO^S,mom(idim))-&
                            w(ixO^S,mag(idir))*w(ixO^S,mag(idim))
      end if
    end do

    ! Get flux of energy
    if(mhd_energy) then
      f(ixO^S,e_)=w(ixO^S,mom(idim))*(wC(ixO^S,e_)+pgas(ixO^S))
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mhd_glm) then
           f(ixO^S,mag(idir))=w(ixO^S,psi_)
        else
           f(ixO^S,mag(idir))=zero
        end if
      else
        f(ixO^S,mag(idir))=w(ixO^S,mom(idim))*w(ixO^S,mag(idir))-w(ixO^S,mag(idim))*w(ixO^S,mom(idir))
      end if
    end do

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixO^S,psi_)  = cmax_global**2*w(ixO^S,mag(idim))
    end if

    ! Contributions of ambipolar term in explicit scheme
    if(mhd_ambipolar_exp.and. .not.stagger_grid) then
      ! ambipolar electric field
      ! E_ambi=-eta_ambi*JxBxB=-JaxBxB=B^2*Ja-(Ja dot B)*B
      !Ja=eta_ambi*J=J * mhd_eta_ambi/rho**2
      allocate(Jambi(ixI^S,1:3))
      call mhd_get_Jambi(w,x,ixI^L,ixO^L,Jambi)
      allocate(btot(ixO^S,1:3))
      btot(ixO^S,1:3) = w(ixO^S,mag(1:3))
      allocate(tmp2(ixO^S),tmp3(ixO^S))
      !tmp2 = Btot^2
      tmp2(ixO^S) = sum(btot(ixO^S,1:3)**2,dim=ndim+1)
      !tmp3 = J_ambi dot Btot
      tmp3(ixO^S) = sum(Jambi(ixO^S,:)*btot(ixO^S,:),dim=ndim+1)

      select case(idim)
        case(1)
          tmp(ixO^S)=w(ixO^S,mag(3)) *Jambi(ixO^S,2) - w(ixO^S,mag(2)) * Jambi(ixO^S,3)
          f(ixO^S,mag(2))= f(ixO^S,mag(2)) - tmp2(ixO^S) * Jambi(ixO^S,3) + tmp3(ixO^S) * btot(ixO^S,3)
          f(ixO^S,mag(3))= f(ixO^S,mag(3)) + tmp2(ixO^S) * Jambi(ixO^S,2) - tmp3(ixO^S) * btot(ixO^S,2)
        case(2)
          tmp(ixO^S)=w(ixO^S,mag(1)) *Jambi(ixO^S,3) - w(ixO^S,mag(3)) * Jambi(ixO^S,1)
          f(ixO^S,mag(1))= f(ixO^S,mag(1)) + tmp2(ixO^S) * Jambi(ixO^S,3) - tmp3(ixO^S) * btot(ixO^S,3)
          f(ixO^S,mag(3))= f(ixO^S,mag(3)) - tmp2(ixO^S) * Jambi(ixO^S,1) + tmp3(ixO^S) * btot(ixO^S,1)
        case(3)
          tmp(ixO^S)=w(ixO^S,mag(2)) *Jambi(ixO^S,1) - w(ixO^S,mag(1)) * Jambi(ixO^S,2)
          f(ixO^S,mag(1))= f(ixO^S,mag(1)) - tmp2(ixO^S) * Jambi(ixO^S,2) + tmp3(ixO^S) * btot(ixO^S,2)
          f(ixO^S,mag(2))= f(ixO^S,mag(2)) + tmp2(ixO^S) * Jambi(ixO^S,1) - tmp3(ixO^S) * btot(ixO^S,1)
      endselect

      if(mhd_energy .and. .not. mhd_internal_e) then
        f(ixO^S,e_) = f(ixO^S,e_) + tmp2(ixO^S) *  tmp(ixO^S)
      endif

      deallocate(Jambi,btot,tmp2,tmp3)
    endif

  end subroutine mhd_get_flux_hde

  !> Calculate fluxes within ixO^L with possible splitting
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

    double precision             :: pgas(ixO^S), ptotal(ixO^S), B(ixO^S,1:ndir)
    double precision             :: tmp(ixI^S)
    double precision             :: vHall(ixI^S,1:ndir)
    integer                      :: idirmin, iw, idir, jdir, kdir
    double precision, allocatable, dimension(:^D&,:) :: Jambi, btot
    double precision, allocatable, dimension(:^D&) :: tmp2, tmp3
    double precision :: tmp4(ixO^S)


    call mhd_get_rho(w,x,ixI^L,ixO^L,tmp)
    ! Get flux of density
    f(ixO^S,rho_)=w(ixO^S,mom(idim))*tmp(ixO^S)
    ! pgas is time dependent only
    if(mhd_energy) then
      pgas=w(ixO^S,p_)
    else
      pgas(ixO^S)=mhd_adiab*tmp(ixO^S)**mhd_gamma
      if(has_equi_pe0) then
        pgas(ixO^S)=pgas(ixO^S)-block%equi_vars(ixO^S,equi_pe0_,b0i)
      endif
    end if

    if (mhd_Hall) then
      call mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
    end if

    if(B0field) then
      B(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))+block%B0(ixO^S,1:ndir,idim)
      pgas=pgas+sum(w(ixO^S,mag(:))*block%B0(ixO^S,:,idim),dim=ndim+1)
    else
      B(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))
    end if

    ptotal = pgas + 0.5d0*sum(w(ixO^S, mag(:))**2, dim=ndim+1)

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixO^S,tracer(iw))=w(ixO^S,mom(idim))*w(ixO^S,tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    if(B0field) then
      do idir=1,ndir
        if(idim==idir) then
          f(ixO^S,mom(idir))=wC(ixO^S,mom(idir))*w(ixO^S,mom(idim))+ptotal(ixO^S)-&
                              w(ixO^S,mag(idir))*B(ixO^S,idim)-&
                       block%B0(ixO^S,idir,idim)*w(ixO^S,mag(idim))
        else
          f(ixO^S,mom(idir))=wC(ixO^S,mom(idir))*w(ixO^S,mom(idim))-&
                              w(ixO^S,mag(idir))*B(ixO^S,idim)-&
                       block%B0(ixO^S,idir,idim)*w(ixO^S,mag(idim))
        end if
      end do
    else
      do idir=1,ndir
        if(idim==idir) then
          f(ixO^S,mom(idir))=wC(ixO^S,mom(idir))*w(ixO^S,mom(idim))+ptotal(ixO^S)-&
                              w(ixO^S,mag(idir))*B(ixO^S,idim)
        else
          f(ixO^S,mom(idir))=wC(ixO^S,mom(idir))*w(ixO^S,mom(idim))-&
                              w(ixO^S,mag(idir))*B(ixO^S,idim)
        end if
      end do
    end if

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if(mhd_energy) then
      if (mhd_internal_e) then
         f(ixO^S,e_)=w(ixO^S,mom(idim))*w(ixO^S,p_)
         if (mhd_Hall) then
            call mpistop("solve internal energy not implemented for Hall MHD")
         endif
      else
        f(ixO^S,e_)=w(ixO^S,mom(idim))*(wC(ixO^S,e_)+ptotal(ixO^S))&
           -B(ixO^S,idim)*sum(w(ixO^S,mag(:))*w(ixO^S,mom(:)),dim=ndim+1)
        if(mhd_solve_eaux) f(ixO^S,eaux_)=w(ixO^S,mom(idim))*wC(ixO^S,eaux_)

        if (mhd_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
           if (mhd_etah>zero) then
              f(ixO^S,e_) = f(ixO^S,e_) + vHall(ixO^S,idim) * &
                 sum(w(ixO^S, mag(:))*B(ixO^S,:),dim=ndim+1) &
                 - B(ixO^S,idim) * sum(vHall(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1)
           end if
        end if
      end if
      if(has_equi_pe0) then
        f(ixO^S,e_)=  f(ixO^S,e_) &
          + w(ixO^S,mom(idim)) * block%equi_vars(ixO^S,equi_pe0_,idim) * inv_gamma_1
      end if
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mhd_glm) then
           f(ixO^S,mag(idir))=w(ixO^S,psi_)
        else
           f(ixO^S,mag(idir))=zero
        end if
      else
        f(ixO^S,mag(idir))=w(ixO^S,mom(idim))*B(ixO^S,idir)-B(ixO^S,idim)*w(ixO^S,mom(idir))

        if (mhd_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (mhd_etah>zero) then
            f(ixO^S,mag(idir)) = f(ixO^S,mag(idir)) &
                 - vHall(ixO^S,idir)*B(ixO^S,idim) &
                 + vHall(ixO^S,idim)*B(ixO^S,idir)
          end if
        end if

      end if
    end do

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixO^S,psi_)  = cmax_global**2*w(ixO^S,mag(idim))
    end if

    ! Contributions of ambipolar term in explicit scheme
    if(mhd_ambipolar_exp.and. .not.stagger_grid) then
      ! ambipolar electric field
      ! E_ambi=-eta_ambi*JxBxB=-JaxBxB=B^2*Ja-(Ja dot B)*B
      !Ja=eta_ambi*J=J * mhd_eta_ambi/rho**2
      allocate(Jambi(ixI^S,1:3))
      call mhd_get_Jambi(w,x,ixI^L,ixO^L,Jambi)
      allocate(btot(ixO^S,1:3))
      if(B0field) then
        do idir=1,3
          btot(ixO^S, idir) = w(ixO^S,mag(idir)) + block%B0(ixO^S,idir,idim)
        enddo
      else
        btot(ixO^S,1:3) = w(ixO^S,mag(1:3))
      endif
      allocate(tmp2(ixO^S),tmp3(ixO^S))
      !tmp2 = Btot^2
      tmp2(ixO^S) = sum(btot(ixO^S,1:3)**2,dim=ndim+1)
      !tmp3 = J_ambi dot Btot
      tmp3(ixO^S) = sum(Jambi(ixO^S,:)*btot(ixO^S,:),dim=ndim+1)

      select case(idim)
        case(1)
          tmp(ixO^S)=w(ixO^S,mag(3)) *Jambi(ixO^S,2) - w(ixO^S,mag(2)) * Jambi(ixO^S,3)
          if(B0field) tmp4(ixO^S) = w(ixO^S,mag(2)) * btot(ixO^S,3) - w(ixO^S,mag(3)) * btot(ixO^S,2)
          f(ixO^S,mag(2))= f(ixO^S,mag(2)) - tmp2(ixO^S) * Jambi(ixO^S,3) + tmp3(ixO^S) * btot(ixO^S,3)
          f(ixO^S,mag(3))= f(ixO^S,mag(3)) + tmp2(ixO^S) * Jambi(ixO^S,2) - tmp3(ixO^S) * btot(ixO^S,2)
        case(2)
          tmp(ixO^S)=w(ixO^S,mag(1)) *Jambi(ixO^S,3) - w(ixO^S,mag(3)) * Jambi(ixO^S,1)
          if(B0field) tmp4(ixO^S) = w(ixO^S,mag(3)) * btot(ixO^S,1) - w(ixO^S,mag(1)) * btot(ixO^S,3)
          f(ixO^S,mag(1))= f(ixO^S,mag(1)) + tmp2(ixO^S) * Jambi(ixO^S,3) - tmp3(ixO^S) * btot(ixO^S,3)
          f(ixO^S,mag(3))= f(ixO^S,mag(3)) - tmp2(ixO^S) * Jambi(ixO^S,1) + tmp3(ixO^S) * btot(ixO^S,1)
        case(3)
          tmp(ixO^S)=w(ixO^S,mag(2)) *Jambi(ixO^S,1) - w(ixO^S,mag(1)) * Jambi(ixO^S,2)
          if(B0field) tmp4(ixO^S) = w(ixO^S,mag(1)) * btot(ixO^S,2) - w(ixO^S,mag(2)) * btot(ixO^S,1)
          f(ixO^S,mag(1))= f(ixO^S,mag(1)) - tmp2(ixO^S) * Jambi(ixO^S,2) + tmp3(ixO^S) * btot(ixO^S,2)
          f(ixO^S,mag(2))= f(ixO^S,mag(2)) + tmp2(ixO^S) * Jambi(ixO^S,1) - tmp3(ixO^S) * btot(ixO^S,1)
      endselect

      if(mhd_energy .and. .not. mhd_internal_e) then
        f(ixO^S,e_) = f(ixO^S,e_) + tmp2(ixO^S) *  tmp(ixO^S)
        if(B0field) f(ixO^S,e_) = f(ixO^S,e_) +  tmp3(ixO^S) *  tmp4(ixO^S)
      endif

      deallocate(Jambi,btot,tmp2,tmp3)
    endif

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

    double precision             :: pgas(ixO^S)
    double precision             :: SA(ixO^S), E(ixO^S,1:ndir), B(ixO^S,1:ndir)
    integer                      :: idirmin, iw, idir, jdir, kdir

    ! gas thermal pressure
    if(mhd_energy) then
      pgas=w(ixO^S,p_)
    else
      pgas=mhd_adiab*w(ixO^S,rho_)**mhd_gamma
    end if

    ! Get flux of density
    f(ixO^S,rho_)=w(ixO^S,mom(idim))*w(ixO^S,rho_)

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixO^S,tracer(iw))=w(ixO^S,mom(idim))*w(ixO^S,tracer(iw))
    end do
    ! E=-uxB=Bxu
    if(B0field) then
      B(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))+block%B0(ixO^S,1:ndir,idim)
      pgas=pgas+sum(w(ixO^S,mag(:))*block%B0(ixO^S,:,idim),dim=ndim+1)
    else
      B(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))
    end if
    E=0.d0
    do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
      if(lvc(idir,jdir,kdir)==1)then
        E(ixO^S,idir)=E(ixO^S,idir)+B(ixO^S,jdir)*w(ixO^S,mom(kdir))
      else if(lvc(idir,jdir,kdir)==-1)then
        E(ixO^S,idir)=E(ixO^S,idir)-B(ixO^S,jdir)*w(ixO^S,mom(kdir))
      end if
    end do; end do; end do

    pgas=pgas+half*(sum(w(ixO^S, mag(:))**2, dim=ndim+1)+&
             sum(E(ixO^S,:)**2,dim=ndim+1)*inv_squared_c)

    ! Get flux of momentum
    if(B0field) then
      do idir=1,ndir
        if(idim==idir) then
          f(ixO^S,mom(idir))=w(ixO^S,rho_)*w(ixO^S,mom(idir))*w(ixO^S,mom(idim))+pgas&
           -w(ixO^S,mag(idir))*B(ixO^S,idim)-E(ixO^S,idir)*E(ixO^S,idim)*inv_squared_c&
           -block%B0(ixO^S,idir,idim)*w(ixO^S,mag(idim))
        else
          f(ixO^S,mom(idir))=w(ixO^S,rho_)*w(ixO^S,mom(idir))*w(ixO^S,mom(idim))&
           -w(ixO^S,mag(idir))*B(ixO^S,idim)-E(ixO^S,idir)*E(ixO^S,idim)*inv_squared_c&
           -block%B0(ixO^S,idir,idim)*w(ixO^S,mag(idim))
        end if
      end do
    else
      do idir=1,ndir
        if(idim==idir) then
          f(ixO^S,mom(idir))=w(ixO^S,rho_)*w(ixO^S,mom(idir))*w(ixO^S,mom(idim))+pgas&
           -w(ixO^S,mag(idir))*B(ixO^S,idim)-E(ixO^S,idir)*E(ixO^S,idim)*inv_squared_c
        else
          f(ixO^S,mom(idir))=w(ixO^S,rho_)*w(ixO^S,mom(idir))*w(ixO^S,mom(idim))&
           -w(ixO^S,mag(idir))*B(ixO^S,idim)-E(ixO^S,idir)*E(ixO^S,idim)*inv_squared_c
        end if
      end do
    end if

    ! Get flux of energy
    if(mhd_energy) then
      SA=0.d0
      do jdir=1,ndir; do kdir=1,ndir
        if(lvc(idim,jdir,kdir)==1)then
          SA=SA+E(ixO^S,jdir)*w(ixO^S,mag(kdir))
        else if(lvc(idim,jdir,kdir)==-1) then
          SA=SA-E(ixO^S,jdir)*w(ixO^S,mag(kdir))
        end if
      end do; end do
      f(ixO^S,e_)=w(ixO^S,mom(idim))*(half*w(ixO^S,rho_)*sum(w(ixO^S,mom(:))**2,dim=ndim+1)+&
                  mhd_gamma*pgas*inv_gamma_1)+SA
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mhd_glm) then
           f(ixO^S,mag(idir))=w(ixO^S,psi_)
        else
           f(ixO^S,mag(idir))=zero
        end if
      else
        f(ixO^S,mag(idir))=w(ixO^S,mom(idim))*B(ixO^S,idir)-B(ixO^S,idim)*w(ixO^S,mom(idir))
      end if
    end do

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixO^S,psi_)  = cmax_global**2*w(ixO^S,mag(idim))
    end if

  end subroutine mhd_get_flux_semirelati

  !> Source terms J.E in internal energy. 
  !> For the ambipolar term E = ambiCoef * JxBxB=ambiCoef * B^2(-J_perpB) 
  !=> the source term J.E = ambiCoef * B^2 * J_perpB^2 = ambiCoef * JxBxB^2/B^2
  !> ambiCoef is calculated as mhd_ambi_coef/rho^2,  see also the subroutine mhd_get_Jambi
  subroutine add_source_ambipolar_internal_energy(qdt,ixI^L,ixO^L,wCT,w,x,ie)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: tmp(ixI^S)
    double precision :: jxbxb(ixI^S,1:3)

    call mhd_get_jxbxb(wCT,x,ixI^L,ixO^L,jxbxb)
    tmp(ixO^S) = sum(jxbxb(ixO^S,1:3)**2,dim=ndim+1) / mhd_mag_en_all(wCT, ixI^L, ixO^L)
    call multiplyAmbiCoef(ixI^L,ixO^L,tmp,wCT,x)   
    w(ixO^S,ie)=w(ixO^S,ie)+qdt * tmp

  end subroutine add_source_ambipolar_internal_energy

  subroutine mhd_get_jxbxb(w,x,ixI^L,ixO^L,res)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: res(:^D&,:)

    double precision  :: btot(ixI^S,1:3)
    integer          :: idir, idirmin
    double precision :: current(ixI^S,7-2*ndir:3)
    double precision :: tmp(ixI^S),b2(ixI^S)

    res=0.d0
    ! Calculate current density and idirmin
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    !!!here we know that current has nonzero values only for components in the range idirmin, 3
 
    if(B0field) then
      do idir=1,3
        btot(ixO^S, idir) = w(ixO^S,mag(idir)) + block%B0(ixO^S,idir,b0i)
      enddo
    else
      btot(ixO^S,1:3) = w(ixO^S,mag(1:3))
    endif

    tmp(ixO^S) = sum(current(ixO^S,idirmin:3)*btot(ixO^S,idirmin:3),dim=ndim+1) !J.B
    b2(ixO^S) = sum(btot(ixO^S,1:3)**2,dim=ndim+1) !B^2
    do idir=1,idirmin-1
      res(ixO^S,idir) = btot(ixO^S,idir) * tmp(ixO^S)
    enddo
    do idir=idirmin,3
      res(ixO^S,idir) = btot(ixO^S,idir) * tmp(ixO^S) - current(ixO^S,idir) * b2(ixO^S)
    enddo
  end subroutine mhd_get_jxbxb

  !> Sets the sources for the ambipolar
  !> this is used for the STS method
  ! The sources are added directly (instead of fluxes as in the explicit) 
  !> at the corresponding indices
  !>  store_flux_var is explicitly called for each of the fluxes one by one
  subroutine sts_set_source_ambipolar(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: ixI^L, ixO^L,igrid,nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step

    double precision, dimension(ixI^S,1:3) :: tmp,ff
    double precision :: fluxall(ixI^S,1:nflux,1:ndim)
    double precision :: fE(ixI^S,7-2*ndim:3)
    double precision  :: btot(ixI^S,1:3),tmp2(ixI^S)
    integer :: i, ixA^L, ie_

    ixA^L=ixO^L^LADD1;

    fluxall=zero

    call mhd_get_jxbxb(w,x,ixI^L,ixA^L,tmp)

    !set electric field in tmp: E=nuA * jxbxb, where nuA=-etaA/rho^2
    do i=1,3
      !tmp(ixA^S,i) = -(mhd_eta_ambi/w(ixA^S, rho_)**2) * tmp(ixA^S,i)
      call multiplyAmbiCoef(ixI^L,ixA^L,tmp(ixI^S,i),w,x)   
    enddo

    if(mhd_energy .and. .not.mhd_internal_e) then
      !btot should be only mag. pert.
      btot(ixA^S,1:3)=0.d0
      !if(B0field) then
      !  do i=1,ndir
      !    btot(ixA^S, i) = w(ixA^S,mag(i)) + block%B0(ixA^S,i,0)
      !  enddo
      !else
        btot(ixA^S,1:ndir) = w(ixA^S,mag(1:ndir))
      !endif
      call cross_product(ixI^L,ixA^L,tmp,btot,ff)
      call get_flux_on_cell_face(ixI^L,ixO^L,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixI^S,1,1:ndim)=ff(ixI^S,1:ndim)
      !- sign comes from the fact that the flux divergence is a source now
      wres(ixO^S,e_)=-tmp2(ixO^S)
    endif

    if(stagger_grid) then
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
    double precision, intent(out)      :: fE(ixI^S,7-2*ndim:3)
    double precision, intent(out)      :: circ(ixI^S,1:ndim)

    integer                            :: hxC^L,ixC^L,ixA^L
    integer                            :: idim1,idim2,idir,ix^D

    fE=zero
    ! calcuate ambipolar electric field on cell edges from cell centers
    do idir=7-2*ndim,3
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
        do idir=7-2*ndim,3 ! Direction of line integral
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

  !> Calculates the explicit dt for the ambipokar term
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
    coef = maxval(abs(tmp(ixO^S)))
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
  !> by implemneting usr_mask_ambipolar subroutine
  subroutine multiplyAmbiCoef(ixI^L,ixO^L,res,w,x)   
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: res(ixI^S)
    double precision :: tmp(ixI^S)
    double precision :: rho(ixI^S)

    call mhd_get_rho(w,x,ixI^L,ixO^L,rho)
    tmp=0.d0
    tmp(ixO^S)=-mhd_eta_ambi/rho(ixO^S)**2
    if (associated(usr_mask_ambipolar)) then
      call usr_mask_ambipolar(ixI^L,ixO^L,w,x,tmp)
    end if

    res(ixO^S) = tmp(ixO^S) * res(ixO^S)
  end subroutine multiplyAmbiCoef

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine mhd_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active,wCTprim)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active
    double precision, intent(in), optional :: wCTprim(ixI^S,1:nw)

    if (.not. qsourcesplit) then
      if(mhd_internal_e) then
        ! Source for solving internal energy
        active = .true.
        call internal_energy_add_source(qdt,ixI^L,ixO^L,wCT,w,x,e_)
      else 
        if(mhd_solve_eaux) then
          ! Source for auxiliary internal energy equation
          active = .true.
          call internal_energy_add_source(qdt,ixI^L,ixO^L,wCT,w,x,eaux_)
        endif
        if(has_equi_pe0) then
          active = .true.
          call add_pe0_divv(qdt,ixI^L,ixO^L,wCT,w,x)
        endif
      endif

      ! Source for B0 splitting
      if (B0field.or.B0field) then
        active = .true.
        call add_source_B0split(qdt,ixI^L,ixO^L,wCT,w,x)
      end if

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
      if (abs(mhd_eta)>smalldouble)then
        active = .true.
        call add_source_res2(qdt,ixI^L,ixO^L,wCT,w,x)
      end if

      if (mhd_eta_hyper>0.d0)then
        active = .true.
        call add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
      end if
      ! add sources for semirelativistic MHD
      if (unsplit_semirelativistic) then
        active = .true.
        call add_source_semirelativistic(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
      end if
      ! add sources for hydrodynamic energy version of MHD
      if (mhd_hydrodynamic_e) then
        active = .true.
        call add_source_hydrodynamic_e(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
      end if
    end if

      {^NOONED
    if(.not.source_split_divb .and. .not.qsourcesplit .and. istep==nstep) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
      case (divb_glm)
        active = .true.
        call add_source_glm(dt,ixI^L,ixO^L,pso(block%igrid)%w,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(dt,ixI^L,ixO^L,pso(block%igrid)%w,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(dt,ixI^L,ixO^L,pso(block%igrid)%w,w,x)
      case (divb_linde)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pso(block%igrid)%w,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pso(block%igrid)%w,w,x)
        call add_source_janhunen(dt,ixI^L,ixO^L,pso(block%igrid)%w,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pso(block%igrid)%w,w,x)
        call add_source_powel(dt,ixI^L,ixO^L,pso(block%igrid)%w,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pso(block%igrid)%w,w,x)
        call add_source_glm(dt,ixI^L,ixO^L,pso(block%igrid)%w,w,x)
      case (divb_ct)
        continue ! Do nothing
      case (divb_multigrid)
        continue ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    else if(source_split_divb .and. qsourcesplit) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
      case (divb_glm)
        active = .true.
        call add_source_glm(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_linde)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
        call add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
        call add_source_powel(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
        call add_source_glm(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_ct)
        continue ! Do nothing
      case (divb_multigrid)
        continue ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    end if
    }

    if(mhd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,qsourcesplit,active, rc_fl)
    end if

    if(mhd_viscosity) then
      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,mhd_energy,qsourcesplit,active)
    end if

    if(mhd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,total_energy,qsourcesplit,active)
    end if

  end subroutine mhd_add_source

  subroutine add_pe0_divv(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: v(ixI^S,1:ndir)
    double precision                :: divv(ixI^S)

    call mhd_get_v(wCT,x,ixI^L,ixI^L,v)

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixI^L,ixO^L,divv,sixthorder=.true.)
      else
        call divvector(v,ixI^L,ixO^L,divv,fourthorder=.true.)
      end if
    else
     call divvector(v,ixI^L,ixO^L,divv)
    end if
    w(ixO^S,e_)=w(ixO^S,e_)-qdt*block%equi_vars(ixO^S,equi_pe0_,0)*divv(ixO^S)

  end subroutine add_pe0_divv

  !> Compute the Lorentz force (JxB)
  subroutine get_Lorentz_force(ixI^L,ixO^L,w,JxB)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: JxB(ixI^S,3)
    double precision                :: a(ixI^S,3), b(ixI^S,3)
    integer                         :: idir, idirmin
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3)

    b=0.0d0
    do idir = 1, ndir
      b(ixO^S, idir) = mhd_mag_i_all(w, ixI^L, ixO^L,idir)
    end do

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
  subroutine mhd_gamma2_alfven(ixI^L, ixO^L, w, gamma_A2)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision, intent(out) :: gamma_A2(ixO^S)
    double precision              :: rho(ixI^S) 

    ! mhd_get_rho cannot be used as x is not a param
    if(has_equi_rho0) then
      rho(ixO^S) = w(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,b0i)
    else  
      rho(ixO^S) = w(ixO^S,rho_) 
    endif
    ! Compute the inverse of 1 + B^2/(rho * c^2)
    gamma_A2 = 1.0d0/(1.0d0+mhd_mag_en_all(w, ixI^L, ixO^L)/rho(ixO^S)*inv_squared_c)
  end subroutine mhd_gamma2_alfven

  !> Compute 1/sqrt(1+v_A^2/c^2) for semirelativisitic MHD, where v_A is the
  !> Alfven velocity
  function mhd_gamma_alfven(w, ixI^L, ixO^L) result(gamma_A)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: gamma_A(ixO^S)

    call mhd_gamma2_alfven(ixI^L, ixO^L, w, gamma_A)
    gamma_A = sqrt(gamma_A)
  end function mhd_gamma_alfven

  subroutine internal_energy_add_source(qdt,ixI^L,ixO^L,wCT,w,x,ie)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision                :: v(ixI^S,1:ndir),divv(ixI^S)

    call mhd_get_v(wCT,x,ixI^L,ixI^L,v)
    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixI^L,ixO^L,divv,sixthorder=.true.)
      else
        call divvector(v,ixI^L,ixO^L,divv,fourthorder=.true.)
      end if
    else
     call divvector(v,ixI^L,ixO^L,divv)
    end if
    w(ixO^S,ie)=w(ixO^S,ie)-qdt*wCT(ixO^S,ie)*gamma_1*divv(ixO^S)
    if(mhd_ambipolar)then
       call add_source_ambipolar_internal_energy(qdt,ixI^L,ixO^L,wCT,w,x,ie)
    end if
    if(fix_small_values) then
      call mhd_handle_small_ei(w,x,ixI^L,ixO^L,ie,'internal_energy_add_source')
    end if
  end subroutine internal_energy_add_source

  subroutine mhd_get_rho(w,x,ixI^L,ixO^L,rho)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: rho(ixI^S)

    if(has_equi_rho0) then
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

    integer :: idir
    logical :: flag(ixI^S,1:nw)
    double precision              :: rho(ixI^S)

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
        call mhd_get_rho(w,x,ixI^L,ixO^L,rho)
        do idir = 1, ndir
           w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))/rho(ixO^S)
        end do
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if

  end subroutine mhd_handle_small_ei

  !> Source terms after split off time-independent magnetic field
  subroutine add_source_B0split(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
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
      axb(ixO^S,:)=axb(ixO^S,:)*qdt
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
      call mhd_get_v(wCT,x,ixI^L,ixO^L,a(ixI^S,1:ndir))
      call cross_product(ixI^L,ixO^L,a,b,axb)
      axb(ixO^S,:)=axb(ixO^S,:)*qdt
      ! add -(vxB) dot J0 source term in energy equation
      do idir=7-2*ndir,3
        w(ixO^S,e_)=w(ixO^S,e_)-axb(ixO^S,idir)*block%J0(ixO^S,idir)
      end do
      if(mhd_ambipolar) then
        !reuse axb
        call mhd_get_jxbxb(wCT,x,ixI^L,ixO^L,axb)
        ! source J0 * E
        do idir=7-2*ndim,3
          !set electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
          call multiplyAmbiCoef(ixI^L,ixO^L,axb(ixI^S,idir),wCT,x)   
          w(ixO^S,e_)=w(ixO^S,e_)+axb(ixO^S,idir)*block%J0(ixO^S,idir)
        enddo
      endif
    end if

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_B0')

  end subroutine add_source_B0split

  !> Source terms for semirelativistic MHD
  subroutine add_source_semirelativistic(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in), optional :: wCTprim(ixI^S,1:nw)

    double precision :: B(ixI^S,3), v(ixI^S,3), E(ixI^S,3), divE(ixI^S)
    integer :: idir

    ! store B0 magnetic field in b
    if(B0field) then
      B(ixI^S,1:ndir)=wCT(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,0)
    else
      B(ixI^S,1:ndir)=wCT(ixI^S,mag(1:ndir))
    end if
    v(ixI^S,1:ndir)=wCTprim(ixI^S,mom(1:ndir))

    call cross_product(ixI^L,ixI^L,B,v,E)
    ! add source term in momentum equations (1/c0^2-1/c^2)E divE
    call divvector(E,ixI^L,ixO^L,divE)
    divE(ixO^S)=qdt*(inv_squared_c0-inv_squared_c)*divE(ixO^S)
    do idir=1,ndir
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))+E(ixO^S,idir)*divE(ixO^S)
    end do

  end subroutine add_source_semirelativistic

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
    integer :: idir, idirmin, idims
    double precision :: current(ixI^S,7-2*ndir:3)
    double precision :: bu(ixO^S,1:ndir), tmp(ixO^S), b2(ixO^S)
    double precision :: gravity_field(ixI^S,1:ndir)

    B=0.0d0
    do idir = 1, ndir
      B(ixO^S, idir) = wCT(ixO^S,mag(idir))
    end do

    call get_current(wCT,ixI^L,ixO^L,idirmin,current)

    J=0.0d0
    do idir=7-2*ndir,3
      J(ixO^S,idir)=current(ixO^S,idir)
    end do

    ! get Lorentz force JxB
    call cross_product(ixI^L,ixO^L,J,B,JxB)

    if(mhd_semirelativistic) then
      ! (v . nabla) v
      do idir=1,ndir
        do idims=1,ndim
          call gradient(wCTprim(ixI^S,mom(idir)),ixI^L,ixO^L,idims,J(ixI^S,idims))
        end do
        B(ixO^S,idir)=sum(wCTprim(ixO^S,mom(1:ndir))*J(ixO^S,1:ndir),dim=ndim+1)
      end do
      ! nabla p
      do idir=1,ndir
        call gradient(wCTprim(ixI^S,p_),ixI^L,ixO^L,idir,J(ixI^S,idir))
      end do

      if(mhd_gravity) then
        if (.not. associated(usr_gravity)) then
          write(*,*) "mod_usr.t: please point usr_gravity to a subroutine"
          write(*,*) "like the phys_gravity in mod_usr_methods.t"
          call mpistop("add_source_hydrodynamic_e: usr_gravity not defined")
        else
          gravity_field=0.d0
          call usr_gravity(ixI^L,ixO^L,wCT,x,gravity_field(ixI^S,1:ndim))
        end if
        do idir=1,ndir
          B(ixO^S,idir)=wCT(ixO^S,rho_)*(B(ixO^S,idir)-gravity_field(ixO^S,idir))+J(ixO^S,idir)-JxB(ixO^S,idir)
        end do
      else
        do idir=1,ndir
          B(ixO^S,idir)=wCT(ixO^S,rho_)*B(ixO^S,idir)+J(ixO^S,idir)-JxB(ixO^S,idir)
        end do
      end if

      b2(ixO^S)=sum(wCT(ixO^S,mag(:))**2,dim=ndim+1)
      tmp(ixO^S)=sqrt(b2(ixO^S))
      where(tmp(ixO^S)>smalldouble)
        tmp(ixO^S)=1.d0/tmp(ixO^S)
      else where
        tmp(ixO^S)=0.d0
      end where
      ! unit vector of magnetic field
      do idir=1,ndir
        bu(ixO^S,idir)=wCT(ixO^S,mag(idir))*tmp(ixO^S)
      end do

      ! Va^2/c^2
      b2(ixO^S)=b2(ixO^S)/w(ixO^S,rho_)*inv_squared_c
      ! Va^2/c^2 / (1+Va^2/c^2)
      b2(ixO^S)=b2(ixO^S)/(1.d0+b2(ixO^S))
      ! bu . F
      tmp(ixO^S)=sum(bu(ixO^S,1:ndir)*B(ixO^S,1:ndir),dim=ndim+1)
      ! Rempel 2017 ApJ 834, 10 equation (54)
      do idir=1,ndir
        J(ixO^S,idir)=b2(ixO^S)*(B(ixO^S,idir)-bu(ixO^S,idir)*tmp(ixO^S))
      end do
      ! Rempel 2017 ApJ 834, 10 equation (29) add SR force at momentum equation
      do idir=1,ndir
        w(ixO^S,mom(idir))=w(ixO^S,mom(idir))+qdt*J(ixO^S,idir)
      end do
      ! Rempel 2017 ApJ 834, 10 equation (30) add work of Lorentz force and SR force
      w(ixO^S,e_)=w(ixO^S,e_)+qdt*sum(wCTprim(ixO^S,mom(1:ndir))*&
              (JxB(ixO^S,1:ndir)+J(ixO^S,1:ndir)),dim=ndim+1)
    else
      ! add work of Lorentz force
      w(ixO^S,e_)=w(ixO^S,e_)+qdt*sum(wCTprim(ixO^S,mom(1:ndir))*JxB(ixO^S,1:ndir),dim=ndim+1)
    end if

  end subroutine add_source_hydrodynamic_e

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
    if (mhd_4th_order) then
      ixA^L=ixO^L^LADD2;
    else
      ixA^L=ixO^L^LADD1;
    end if

    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixI^L,ixO^L,idirmin,current)

    if (mhd_eta>zero)then
       eta(ixA^S)=mhd_eta
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
       if (mhd_4th_order) then
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
       if (mhd_energy) then
          w(ixO^S,e_)=w(ixO^S,e_)+qdt*tmp(ixO^S)*Bf(ixO^S,idir)
       end if
    end do ! idir

    if(mhd_energy) then
      ! de/dt+=eta*J**2
      tmp(ixO^S)=qdt*eta(ixO^S)*sum(current(ixO^S,:)**2,dim=ndim+1)
      w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)
      if(mhd_solve_eaux) then
        ! add eta*J**2 source term in the internal energy equation
        w(ixO^S,eaux_)=w(ixO^S,eaux_)+tmp(ixO^S)
      end if
    end if

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res1')

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

    if (mhd_eta>zero)then
       eta(ixA^S)=mhd_eta
    else
       call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,current,eta)
    end if

    ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
    tmpvec(ixA^S,1:ndir)=zero
    do idir=idirmin,3
       tmpvec(ixA^S,idir)=current(ixA^S,idir)*eta(ixA^S)
    end do
    curlj=0.d0
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
      ! de/dt= +div(B x Jeta) = eta J^2 - B dot curl(eta J)
      ! de1/dt= eta J^2 - B1 dot curl(eta J)
      tmp(ixO^S)=eta(ixO^S)*sum(current(ixO^S,:)**2,dim=ndim+1)
      w(ixO^S,e_)=w(ixO^S,e_)+qdt*(tmp(ixO^S)-&
        sum(wCT(ixO^S,mag(1:ndir))*curlj(ixO^S,1:ndir),dim=ndim+1))
      if(mhd_solve_eaux) then
        ! add eta*J**2 source term in the internal energy equation
        w(ixO^S,eaux_)=w(ixO^S,eaux_)+tmp(ixO^S)
      end if
    end if

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res2')
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
    ehyper(ixA^S,1:ndir) = - tmpvec(ixA^S,1:ndir)*mhd_eta_hyper

    ixA^L=ixO^L;
    tmpvec2(ixA^S,1:ndir)=zero
    call curlvector(ehyper,ixI^L,ixA^L,tmpvec2,idirmin1,1,3)

    do idir=1,ndir
      w(ixO^S,mag(idir)) = w(ixO^S,mag(idir))-tmpvec2(ixO^S,idir)*qdt
    end do

    if (mhd_energy) then
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
    double precision:: divb(ixI^S)
    integer          :: idim,idir
    double precision :: gradPsi(ixI^S)


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
      ! gradient of Psi
      if(total_energy) then
        do idim=1,ndim
          select case(typegrad)
          case("central")
            call gradient(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
          case("limited")
            call gradientS(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
          end select
          ! e  = e  -qdt (b . grad(Psi))
          w(ixO^S,e_) = w(ixO^S,e_)-qdt*wCT(ixO^S,mag(idim))*gradPsi(ixO^S)
        end do
      end if

      ! We calculate now div B
      call get_divb(wCT,ixI^L,ixO^L,divb, mhd_divb_4thorder)

      ! m = m - qdt b div b
      do idir=1,ndir
        w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-qdt*mhd_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
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
    double precision                :: divb(ixI^S),v(ixI^S,1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb, mhd_divb_4thorder)

    ! calculate velocity
    call mhd_get_v(wCT,x,ixI^L,ixO^L,v)

    if (total_energy) then
      ! e = e - qdt (v . b) * div b
      w(ixO^S,e_)=w(ixO^S,e_)-&
           qdt*sum(v(ixO^S,:)*wCT(ixO^S,mag(:)),dim=ndim+1)*divb(ixO^S)
    end if

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*v(ixO^S,idir)*divb(ixO^S)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-qdt*mhd_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divb(ixI^S),vel(ixI^S,1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb, mhd_divb_4thorder)

    call mhd_get_v(wCT,x,ixI^L,ixO^L,vel)
    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*vel(ixO^S,idir)*divb(ixO^S)
    end do

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_janhunen')

  end subroutine add_source_janhunen

  subroutine add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add Linde's divB related sources to wnew within ixO
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: idim, idir, ixp^L, i^D, iside
    double precision :: divb(ixI^S),graddivb(ixI^S)
    logical, dimension(-1:1^D&) :: leveljump

    ! Calculate div B
    ixp^L=ixO^L^LADD1;
    call get_divb(wCT,ixI^L,ixp^L,divb, mhd_divb_4thorder)

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
         call gradientS(divb,ixI^L,ixp^L,idim,graddivb)
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

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_linde')

  end subroutine add_source_linde

  !> Calculate div B within ixO
  subroutine get_divb(w,ixI^L,ixO^L,divb, fourthorder)

    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: divb(ixI^S)
    logical, intent(in), optional   :: fourthorder

    double precision                   :: bvec(ixI^S,1:ndir)
    double precision                   :: divb_corner(ixI^S), sign
    double precision                   :: aux_vol(ixI^S)
    integer                            :: ixC^L, idir, ic^D, ix^L

    if(stagger_grid) then
      divb=0.d0
      do idir=1,ndim
        ixC^L=ixO^L-kr(idir,^D);
        divb(ixO^S)=divb(ixO^S)+block%ws(ixO^S,idir)*block%surfaceC(ixO^S,idir)-&
                                block%ws(ixC^S,idir)*block%surfaceC(ixC^S,idir)
      end do
      divb(ixO^S)=divb(ixO^S)/block%dvolume(ixO^S)
    else
      bvec(ixI^S,:)=w(ixI^S,mag(:))
      select case(typediv)
      case("central")
        call divvector(bvec,ixI^L,ixO^L,divb,fourthorder)
      case("limited")
        call divvectorS(bvec,ixI^L,ixO^L,divb)
      end select
    end if

  end subroutine get_divb

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
    integer :: idir, idirmin0

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),bvec(ixI^S,1:ndir)

    idirmin0 = 7-2*ndir

    bvec(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))

    call curlvector(bvec,ixI^L,ixO^L,current,idirmin,idirmin0,ndir)

    if(B0field) current(ixO^S,idirmin0:3)=current(ixO^S,idirmin0:3)+&
        block%J0(ixO^S,idirmin0:3)
  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine mhd_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
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

    ^D&dxarr(^D)=dx^D;
    if (mhd_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/mhd_eta
    else if (mhd_eta<zero)then
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

    if(mhd_eta_hyper>zero) then
      if(slab_uniform) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/mhd_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixO^S,1:ndim))**4/mhd_eta_hyper,dtnew)
      end if
    end if

    if(mhd_radiative_cooling) then
      call cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x,rc_fl)
    end if

    if(mhd_viscosity) then
      call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(mhd_gravity) then
      call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(mhd_ambipolar_exp) then
      dtnew=min(dtdiffpar*get_ambipolar_dt(w,ixI^L,ixO^L,dx^D,x),dtnew)
    endif

  end subroutine mhd_get_dt

  ! Add geometrical source terms to w
  subroutine mhd_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),invrho(ixO^S),invr(ixO^S)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    ! 1/rho
    invrho(ixO^S)=1.d0/wCT(ixO^S,rho_)
    invr(ixO^S)=1.d0/x(ixO^S,1)
    select case (coordinate)
    case (cylindrical)
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in MHD")
      endif
      call mhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp)
      if(phi_>0) then
        w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*invr(ixO^S)*(tmp(ixO^S)-&
                  wCT(ixO^S,bphi_)**2+wCT(ixO^S,mphi_)**2*invrho(ixO^S))
        w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt*invr(ixO^S)*(&
                 -wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)*invrho(ixO^S) &
                 +wCT(ixO^S,bphi_)*wCT(ixO^S,br_))
        if(.not.stagger_grid) then
          w(ixO^S,bphi_)=w(ixO^S,bphi_)+qdt*invr(ixO^S)*&
                   (wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
                   -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
                   *invrho(ixO^S)
        end if
      else
        w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*invr(ixO^S)*tmp(ixO^S)
      end if
      if(mhd_glm) w(ixO^S,br_)=w(ixO^S,br_)+qdt*wCT(ixO^S,psi_)*invr(ixO^S)
    case (spherical)
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       call mhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp1)
       ! m1
       tmp(ixO^S)=tmp1(ixO^S)*x(ixO^S,1) &
                  *(block%surfaceC(ixO^S,1)-block%surfaceC(h1x^S,1))/block%dvolume(ixO^S)
       do idir=2,ndir
         tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(idir))**2*invrho(ixO^S)-wCT(ixO^S,mag(idir))**2
       end do
       w(ixO^S,mom(1))=w(ixO^S,mom(1))+qdt*tmp(ixO^S)*invr(ixO^S)
       ! b1
       if(mhd_glm) then
         w(ixO^S,mag(1))=w(ixO^S,mag(1))+qdt*invr(ixO^S)*2.0d0*wCT(ixO^S,psi_)
       end if

       {^NOONED
       ! m2
       ! This will make hydrostatic p=const an exact solution
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp1(ixO^S) &
            *(block%surfaceC(ixO^S,2)-block%surfaceC(h2x^S,2)) &
            /block%dvolume(ixO^S)
       tmp(ixO^S)=-(wCT(ixO^S,mom(1))*wCT(ixO^S,mom(2))*invrho(ixO^S) &
            -wCT(ixO^S,mag(1))*wCT(ixO^S,mag(2)))
       if(ndir==3) then
         tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(3))**2*invrho(ixO^S) &
              -wCT(ixO^S,mag(3))**2)*dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
       end if
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S)*invr(ixO^S)
       ! b2
       if(.not.stagger_grid) then
         tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(2)) &
              -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(1)))*invrho(ixO^S)
         if(mhd_glm) then
           tmp(ixO^S)=tmp(ixO^S) &
                + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
         end if
         w(ixO^S,mag(2))=w(ixO^S,mag(2))+qdt*tmp(ixO^S)*invr(ixO^S)
       end if
       }

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixO^S)=-(wCT(ixO^S,mom(3))*wCT(ixO^S,mom(1))*invrho(ixO^S) &
                -wCT(ixO^S,mag(3))*wCT(ixO^S,mag(1))) {^NOONED &
                -(wCT(ixO^S,mom(2))*wCT(ixO^S,mom(3))*invrho(ixO^S) &
                -wCT(ixO^S,mag(2))*wCT(ixO^S,mag(3))) &
                *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
           w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*tmp(ixO^S)*invr(ixO^S)
         else
           call mpistop("angmomfix not implemented yet in MHD")
         end if
         ! b3
         if(.not.stagger_grid) then
           tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(3)) &
                -wCT(ixO^S,mom(3))*wCT(ixO^S,mag(1)))*invrho(ixO^S) {^NOONED &
                -(wCT(ixO^S,mom(3))*wCT(ixO^S,mag(2)) &
                -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(3)))*dcos(x(ixO^S,2)) &
                /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
           w(ixO^S,mag(3))=w(ixO^S,mag(3))+qdt*tmp(ixO^S)*invr(ixO^S)
         end if
       end if
    end select
  end subroutine mhd_add_source_geom

  ! Add geometrical source terms to w
  subroutine mhd_add_source_geom_split(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),invrho(ixO^S),invr(ixO^S)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    if(has_equi_rho0) then
      invrho(ixO^S) = 1d0/(wCT(ixO^S,rho_) + block%equi_vars(ixO^S,equi_rho0_,b0i))
    else
      invrho(ixO^S) = 1d0/wCT(ixO^S,rho_)
    end if
    invr(ixO^S)=1d0/x(ixO^S,1)

    select case (coordinate)
    case (cylindrical)
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in MHD")
      endif
      call mhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp)
      if(phi_>0) then
        w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*invr(ixO^S)*(tmp(ixO^S)-&
                  wCT(ixO^S,bphi_)**2+wCT(ixO^S,mphi_)**2*invrho(ixO^S))
        w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt*invr(ixO^S)*(&
                 -wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)*invrho(ixO^S) &
                 +wCT(ixO^S,bphi_)*wCT(ixO^S,br_))
        if(.not.stagger_grid) then
          w(ixO^S,bphi_)=w(ixO^S,bphi_)+qdt*invr(ixO^S)*&
                   (wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
                   -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
                   *invrho(ixO^S)
        end if
      else
        w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*invr(ixO^S)*tmp(ixO^S)
      end if
      if(mhd_glm) w(ixO^S,br_)=w(ixO^S,br_)+qdt*wCT(ixO^S,psi_)*invr(ixO^S)
    case (spherical)
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       call mhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp1)
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
       w(ixO^S,mom(1))=w(ixO^S,mom(1))+qdt*tmp(ixO^S)*invr(ixO^S)
       ! b1
       if(mhd_glm) then
         w(ixO^S,mag(1))=w(ixO^S,mag(1))+qdt*invr(ixO^S)*2.0d0*wCT(ixO^S,psi_)
       end if

       {^NOONED
       ! m2
       tmp(ixO^S)=tmp1(ixO^S)
       if(B0field) then
         tmp(ixO^S)=tmp(ixO^S)+tmp2(ixO^S)
       end if
       ! This will make hydrostatic p=const an exact solution
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S) &
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
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S)*invr(ixO^S)
       ! b2
       if(.not.stagger_grid) then
         tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(2)) &
              -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(1)))*invrho(ixO^S)
         if(B0field) then
           tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(1))*block%B0(ixO^S,2,0) &
                -wCT(ixO^S,mom(2))*block%B0(ixO^S,1,0))*invrho(ixO^S)
         end if
         if(mhd_glm) then
           tmp(ixO^S)=tmp(ixO^S) &
                + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
         end if
         w(ixO^S,mag(2))=w(ixO^S,mag(2))+qdt*tmp(ixO^S)*invr(ixO^S)
       end if
       }

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
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
           w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*tmp(ixO^S)*invr(ixO^S)
         else
           call mpistop("angmomfix not implemented yet in MHD")
         end if
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
           w(ixO^S,mag(3))=w(ixO^S,mag(3))+qdt*tmp(ixO^S)*invr(ixO^S)
         end if
       end if
    end select
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

  !> Compute full magnetic field by direction
  function mhd_mag_i_all(w, ixI^L, ixO^L,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idir
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mgf(ixO^S)

    if (B0field) then
      mgf = w(ixO^S, mag(idir))+block%B0(ixO^S,idir,b0i)
    else
      mgf = w(ixO^S, mag(idir))
    end if
  end function mhd_mag_i_all

  !> Compute evolving magnetic energy
  function mhd_mag_en(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    mge = 0.5d0 * sum(w(ixO^S, mag(:))**2, dim=ndim+1)
  end function mhd_mag_en

  !> compute kinetic energy
  function mhd_kin_en_origin(w, ixI^L, ixO^L, inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim,block
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)
    double precision, intent(in), optional :: inv_rho(ixO^S)

    if (present(inv_rho)) then
      ke = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) * inv_rho
    else
      if(has_equi_rho0) then
        ke(ixO^S) = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) / (w(ixO^S, rho_) + block%equi_vars(ixO^S,equi_rho0_,0))
      else
        ke(ixO^S) = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) / w(ixO^S, rho_)
      endif
    end if
  end function mhd_kin_en_origin

  !> compute kinetic energy
  function mhd_kin_en_boris(w, ixI^L, ixO^L, inv_rho) result(ke)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)
    double precision, intent(in), optional :: inv_rho(ixO^S)

    if (present(inv_rho)) then
      ke=1.d0/(1.d0+sum(w(ixO^S,mag(:))**2,dim=ndim+1)*inv_rho*inv_squared_c)
      ke=0.5d0*sum((w(ixO^S, mom(:)))**2,dim=ndim+1)*ke**2*inv_rho
    else
      ke=1.d0/(1.d0+sum(w(ixO^S,mag(:))**2,dim=ndim+1)/w(ixO^S,rho_)*inv_squared_c)
      ke=0.5d0*sum(w(ixO^S, mom(:))**2,dim=ndim+1)*ke**2/w(ixO^S, rho_)
    end if
  end function mhd_kin_en_boris

  subroutine mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: vHall(ixI^S,1:3)

    integer          :: idir, idirmin
    double precision :: current(ixI^S,7-2*ndir:3)
    double precision :: rho(ixI^S)

    call mhd_get_rho(w,x,ixI^L,ixO^L,rho)
    ! Calculate current density and idirmin
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    vHall(ixO^S,1:3) = zero
    vHall(ixO^S,idirmin:3) = - mhd_etah*current(ixO^S,idirmin:3)
    do idir = idirmin, 3
       vHall(ixO^S,idir) = vHall(ixO^S,idir)/rho(ixO^S)
    end do

  end subroutine mhd_getv_Hall

  subroutine mhd_get_Jambi(w,x,ixI^L,ixO^L,res)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, allocatable, intent(inout) :: res(:^D&,:)


    integer          :: idir, idirmin
    double precision :: current(ixI^S,7-2*ndir:3)

    res = 0d0

    ! Calculate current density and idirmin
    call get_current(w,ixI^L,ixO^L,idirmin,current)
 
    res(ixO^S,idirmin:3)=-current(ixO^S,idirmin:3)
    do idir = idirmin, 3
      call multiplyAmbiCoef(ixI^L,ixO^L,res(ixI^S,idir),w,x)   
    enddo

  end subroutine mhd_get_Jambi

    ! COMMENTED  because we have that in cmax now:
!  subroutine mhd_getdt_Hall(w,x,ixI^L,ixO^L,dx^D,dthall)
!    use mod_global_parameters
!
!    integer, intent(in) :: ixI^L, ixO^L
!    double precision, intent(in)    :: dx^D
!    double precision, intent(in)    :: w(ixI^S,1:nw)
!    double precision, intent(in)    :: x(ixI^S,1:ndim)
!    double precision, intent(out)   :: dthall
!    !.. local ..
!    double precision :: dxarr(ndim)
!    double precision :: bmag(ixI^S)
!
!    dthall=bigdouble
!
!
!    ^D&dxarr(^D)=dx^D;
!
!    if (.not. B0field) then
!       bmag(ixO^S)=sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1))
!    else
!       bmag(ixO^S)=sqrt(sum((w(ixO^S,mag(:)) + block%B0(ixO^S,1:ndir,b0i))**2))
!    end if
!
!    if(slab_uniform) then
!      dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(mhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
!    else
!      dthall=dtdiffpar*minval(block%ds(ixO^S,1:ndim))**2.0d0/(mhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
!    end if
!
!  end subroutine mhd_getdt_Hall

  subroutine mhd_modify_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    double precision, intent(inout) :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw)
    type(state)                     :: s
    double precision                :: dB(ixI^S), dPsi(ixI^S)

    if(stagger_grid) then
      wLC(ixO^S,mag(idir))=s%ws(ixO^S,idir)
      wRC(ixO^S,mag(idir))=s%ws(ixO^S,idir)
      wLp(ixO^S,mag(idir))=s%ws(ixO^S,idir)
      wRp(ixO^S,mag(idir))=s%ws(ixO^S,idir)
    else
      ! Solve the Riemann problem for the linear 2x2 system for normal
      ! B-field and GLM_Psi according to Dedner 2002:
      ! This implements eq. (42) in Dedner et al. 2002 JcP 175
      ! Gives the Riemann solution on the interface
      ! for the normal B component and Psi in the GLM-MHD system.
      ! 23/04/2013 Oliver Porth
      dB(ixO^S)   = wRp(ixO^S,mag(idir)) - wLp(ixO^S,mag(idir))
      dPsi(ixO^S) = wRp(ixO^S,psi_) - wLp(ixO^S,psi_)

      wLp(ixO^S,mag(idir))   = 0.5d0 * (wRp(ixO^S,mag(idir)) + wLp(ixO^S,mag(idir))) &
           - 0.5d0/cmax_global * dPsi(ixO^S)
      wLp(ixO^S,psi_)       = 0.5d0 * (wRp(ixO^S,psi_) + wLp(ixO^S,psi_)) &
           - 0.5d0*cmax_global * dB(ixO^S)

      wRp(ixO^S,mag(idir)) = wLp(ixO^S,mag(idir))
      wRp(ixO^S,psi_) = wLp(ixO^S,psi_)

      if(total_energy) then
        wRC(ixO^S,e_)=wRC(ixO^S,e_)-half*wRC(ixO^S,mag(idir))**2
        wLC(ixO^S,e_)=wLC(ixO^S,e_)-half*wLC(ixO^S,mag(idir))**2
      end if
      wRC(ixO^S,mag(idir)) = wLp(ixO^S,mag(idir))
      wRC(ixO^S,psi_) = wLp(ixO^S,psi_)
      wLC(ixO^S,mag(idir)) = wLp(ixO^S,mag(idir))
      wLC(ixO^S,psi_) = wLp(ixO^S,psi_)
      ! modify total energy according to the change of magnetic field
      if(total_energy) then
        wRC(ixO^S,e_)=wRC(ixO^S,e_)+half*wRC(ixO^S,mag(idir))**2
        wLC(ixO^S,e_)=wLC(ixO^S,e_)+half*wLC(ixO^S,mag(idir))**2
      end if
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
       if(total_energy) call mhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(total_energy) call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(2)
       if(total_energy) call mhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(total_energy) call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(3)
       if(total_energy) call mhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(total_energy) call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(4)
       if(total_energy) call mhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(total_energy) call mhd_to_conserved(ixG^L,ixO^L,w,x)
     {^IFTHREED
     case(5)
       if(total_energy) call mhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(total_energy) call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(6)
       if(total_energy) call mhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(total_energy) call mhd_to_conserved(ixG^L,ixO^L,w,x)
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
    integer                      :: iigrid, igrid, id
    integer                      :: n, nc, lvl, ix^L, ixC^L, idim
    type(tree_node), pointer     :: pnode
    double precision             :: tmp(ixG^T), grad(ixG^T, ndim)
    double precision             :: res
    double precision, parameter  :: max_residual = 1d-3
    double precision, parameter  :: residual_reduction = 1d-10
    integer, parameter           :: max_its      = 50
    double precision             :: residual_it(max_its), max_divb

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
            mhd_divb_4thorder)
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
           call gradientx(tmp,ps(igrid)%x,ixG^LL,ixC^L,idim,grad(ixG^T,idim),.false.)
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

  subroutine mhd_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s,vcts)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wprim(ixI^S,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

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

  end subroutine mhd_update_faces

  !> get electric field though averaging neighors to update faces in CT
  subroutine update_faces_average(ixI^L,ixO^L,qt,qdt,fC,fE,sCT,s)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    integer                            :: hxC^L,ixC^L,jxC^L,ixCm^L
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2
    double precision                   :: circ(ixI^S,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixI^S,7-2*ndim:3) :: E_resi, E_ambi

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D-1;

    ! if there is resistivity, get eta J
    if(mhd_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(mhd_ambipolar_exp) call get_ambipolar_electric_field(ixI^L,ixO^L,sCT%w,x,E_ambi)

    fE=zero

    do idim1=1,ndim 
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=7-2*ndim,3! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ! Assemble indices
            jxC^L=ixC^L+kr(idim1,^D);
            hxC^L=ixC^L+kr(idim2,^D);
            ! Interpolate to edges
            fE(ixC^S,idir)=quarter*(fC(ixC^S,iwdim1,idim2)+fC(jxC^S,iwdim1,idim2)&
                                   -fC(ixC^S,iwdim2,idim1)-fC(hxC^S,iwdim2,idim1))

            ! add resistive electric field at cell edges E=-vxB+eta J
            if(mhd_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+E_resi(ixC^S,idir)
            ! add ambipolar electric field
            if(mhd_ambipolar_exp) fE(ixC^S,idir)=fE(ixC^S,idir)+E_ambi(ixC^S,idir)

            fE(ixC^S,idir)=qdt*s%dsC(ixC^S,idir)*fE(ixC^S,idir)

            if (.not.slab) then
              where(abs(x(ixC^S,r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixC^S,idir)=zero
              end where
            end if
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
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

    end associate

  end subroutine update_faces_average

  !> update faces using UCT contact mode by Gardiner and Stone 2005 JCP 205, 509
  subroutine update_faces_contact(ixI^L,ixO^L,qt,qdt,wp,fC,fE,sCT,s,vcts)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wp(ixI^S,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    double precision                   :: circ(ixI^S,1:ndim)
    ! electric field at cell centers
    double precision                   :: ECC(ixI^S,7-2*ndim:3)
    ! gradient of E at left and right side of a cell face
    double precision                   :: EL(ixI^S),ER(ixI^S)
    ! gradient of E at left and right side of a cell corner
    double precision                   :: ELC(ixI^S),ERC(ixI^S)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixI^S,7-2*ndim:3) :: E_resi, E_ambi
    ! total magnetic field at cell centers
    double precision                   :: Btot(ixI^S,1:ndim)
    integer                            :: hxC^L,ixC^L,jxC^L,ixA^L,ixB^L
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(bfaces=>s%ws,x=>s%x,w=>s%w,vnorm=>vcts%vnorm)

    if(B0field) then
      Btot(ixI^S,1:ndim)=wp(ixI^S,mag(1:ndim))+block%B0(ixI^S,1:ndim,0)
    else
      Btot(ixI^S,1:ndim)=wp(ixI^S,mag(1:ndim))
    end if
    ECC=0.d0
    ! Calculate electric field at cell centers
    do idim1=1,ndim; do idim2=1,ndim; do idir=7-2*ndim,3
      if(lvc(idim1,idim2,idir)==1)then
         ECC(ixI^S,idir)=ECC(ixI^S,idir)+Btot(ixI^S,idim1)*wp(ixI^S,mom(idim2))
      else if(lvc(idim1,idim2,idir)==-1) then
         ECC(ixI^S,idir)=ECC(ixI^S,idir)-Btot(ixI^S,idim1)*wp(ixI^S,mom(idim2))
      endif
    enddo; enddo; enddo

    ! if there is resistivity, get eta J
    if(mhd_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(mhd_ambipolar_exp) call get_ambipolar_electric_field(ixI^L,ixO^L,sCT%w,x,E_ambi)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    fE=zero
    ! evaluate electric field along cell edges according to equation (41)
    do idim1=1,ndim 
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ixCmax^D=ixOmax^D;
            ixCmin^D=ixOmin^D+kr(idir,^D)-1;
            ! Assemble indices
            jxC^L=ixC^L+kr(idim1,^D);
            hxC^L=ixC^L+kr(idim2,^D);
            ! average cell-face electric field to cell edges
            fE(ixC^S,idir)=quarter*&
            (fC(ixC^S,iwdim1,idim2)+fC(jxC^S,iwdim1,idim2)&
            -fC(ixC^S,iwdim2,idim1)-fC(hxC^S,iwdim2,idim1))

            ! add slope in idim2 direction from equation (50)
            ixAmin^D=ixCmin^D;
            ixAmax^D=ixCmax^D+kr(idim1,^D);
            EL(ixA^S)=fC(ixA^S,iwdim1,idim2)-ECC(ixA^S,idir)
            hxC^L=ixA^L+kr(idim2,^D);
            ER(ixA^S)=fC(ixA^S,iwdim1,idim2)-ECC(hxC^S,idir)
            where(vnorm(ixC^S,idim1)>0.d0)
              ELC(ixC^S)=EL(ixC^S)
            else where(vnorm(ixC^S,idim1)<0.d0)
              ELC(ixC^S)=EL(jxC^S)
            else where
              ELC(ixC^S)=0.5d0*(EL(ixC^S)+EL(jxC^S))
            end where
            hxC^L=ixC^L+kr(idim2,^D);
            where(vnorm(hxC^S,idim1)>0.d0)
              ERC(ixC^S)=ER(ixC^S)
            else where(vnorm(hxC^S,idim1)<0.d0)
              ERC(ixC^S)=ER(jxC^S)
            else where
              ERC(ixC^S)=0.5d0*(ER(ixC^S)+ER(jxC^S))
            end where
            fE(ixC^S,idir)=fE(ixC^S,idir)+0.25d0*(ELC(ixC^S)+ERC(ixC^S))

            ! add slope in idim1 direction from equation (50)
            jxC^L=ixC^L+kr(idim2,^D);
            ixAmin^D=ixCmin^D;
            ixAmax^D=ixCmax^D+kr(idim2,^D);
            EL(ixA^S)=-fC(ixA^S,iwdim2,idim1)-ECC(ixA^S,idir)
            hxC^L=ixA^L+kr(idim1,^D);
            ER(ixA^S)=-fC(ixA^S,iwdim2,idim1)-ECC(hxC^S,idir)
            where(vnorm(ixC^S,idim2)>0.d0)
              ELC(ixC^S)=EL(ixC^S)
            else where(vnorm(ixC^S,idim2)<0.d0)
              ELC(ixC^S)=EL(jxC^S)
            else where
              ELC(ixC^S)=0.5d0*(EL(ixC^S)+EL(jxC^S))
            end where
            hxC^L=ixC^L+kr(idim1,^D);
            where(vnorm(hxC^S,idim2)>0.d0)
              ERC(ixC^S)=ER(ixC^S)
            else where(vnorm(hxC^S,idim2)<0.d0)
              ERC(ixC^S)=ER(jxC^S)
            else where
              ERC(ixC^S)=0.5d0*(ER(ixC^S)+ER(jxC^S))
            end where
            fE(ixC^S,idir)=fE(ixC^S,idir)+0.25d0*(ELC(ixC^S)+ERC(ixC^S))

            ! add resistive electric field at cell edges E=-vxB+eta J
            if(mhd_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+E_resi(ixC^S,idir)
            ! add ambipolar electric field
            if(mhd_ambipolar_exp) fE(ixC^S,idir)=fE(ixC^S,idir)+E_ambi(ixC^S,idir)

            ! times time step and edge length 
            fE(ixC^S,idir)=fE(ixC^S,idir)*qdt*s%dsC(ixC^S,idir)
            if (.not.slab) then
              where(abs(x(ixC^S,r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixC^S,idir)=zero
              end where
            end if
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
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
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
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)
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
    double precision, dimension(ixI^S,7-2*ndim:3) :: E_resi, E_ambi
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
    if(mhd_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(mhd_ambipolar_exp) call get_ambipolar_electric_field(ixI^L,ixO^L,sCT%w,x,E_ambi)

    fE=zero

    do idir=7-2*ndim,3
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
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update
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
    double precision :: jce(ixI^S,7-2*ndim:3)

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
        do idir=7-2*ndim,3
          if (lvc(idim1,idim2,idir)==0) cycle
          ixCmax^D=ixOmax^D;
          ixCmin^D=ixOmin^D+kr(idir,^D)-1;
          ixBmax^D=ixCmax^D-kr(idir,^D)+1;
          ixBmin^D=ixCmin^D;
          ! current at transverse faces
          xs(ixB^S,:)=x(ixB^S,:)
          xs(ixB^S,idim2)=x(ixB^S,idim2)+half*dx(ixB^S,idim2)
          call gradientx(wCTs(ixGs^T,idim2),xs,ixGs^LL,ixC^L,idim1,gradi,.true.)
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
      call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,jcc,eta)
      ! calcuate eta on cell edges
      do idir=7-2*ndim,3
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
      enddo
    end if

    end associate
  end subroutine get_resistive_electric_field

  !> get ambipolar electric field on cell edges
  subroutine get_ambipolar_electric_field(ixI^L,ixO^L,w,x,fE)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: fE(ixI^S,7-2*ndim:3)

    double precision :: jxbxb(ixI^S,1:3)
    integer :: idir,ixA^L,ixC^L,ix^D

    ixA^L=ixO^L^LADD1;
    call mhd_get_jxbxb(w,x,ixI^L,ixA^L,jxbxb)
    ! calcuate electric field on cell edges from cell centers
    do idir=7-2*ndim,3
      !set electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
      !jxbxb(ixA^S,i) = -(mhd_eta_ambi/w(ixA^S, rho_)**2) * jxbxb(ixA^S,i)
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
    integer, intent(in)                :: ixO^L
    type(state)                        :: s

    integer                            :: fxO^L, gxO^L, hxO^L, jxO^L, kxO^L, idim

    associate(w=>s%w, ws=>s%ws)

    ! calculate cell-center values from face-center values in 2nd order
    do idim=1,ndim
      ! Displace index to the left
      ! Even if ixI^L is the full size of the w arrays, this is ok
      ! because the staggered arrays have an additional place to the left.
      hxO^L=ixO^L-kr(idim,^D);
      ! Interpolate to cell barycentre using arithmetic average
      ! This might be done better later, to make the method less diffusive.
      w(ixO^S,mag(idim))=half/s%surface(ixO^S,idim)*&
        (ws(ixO^S,idim)*s%surfaceC(ixO^S,idim)&
        +ws(hxO^S,idim)*s%surfaceC(hxO^S,idim))
    end do

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

    end associate

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
