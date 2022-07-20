!> Magneto-hydrodynamics module
module mod_twofl_phys

#include "amrvac.h"

  use mod_physics
  use mod_global_parameters, only: std_len
  use mod_thermal_conduction, only: tc_fluid
  use mod_radiative_cooling, only: rc_fluid
  use mod_thermal_emission, only: te_fluid
  implicit none
  private
  !! E_c = E_kin + E_mag + E_int
  !! E_n = E_kin + E_int
  integer, public, parameter              :: EQ_ENERGY_TOT=2
  !! E_c = E_int
  !! E_n = E_int
  integer, public, parameter              :: EQ_ENERGY_INT=1
  !! E_n, E_c are calculated from density as c_adiab rho^gamma
  !! No energy equation => no variable assigned for it
  integer, public, parameter              :: EQ_ENERGY_NONE=0
  !! E_c = E_kin + E_int
  !! E_n = E_kin + E_int
  integer, public, parameter              :: EQ_ENERGY_KI=3
  !! additional variable for the charges energy at index eaux_
  !! E_c (index e_) = E_kin + E_mag + E_int, E_c (index eaux_) = E_int
  !! E_n (index e_) = E_kin + E_int
  integer, public, parameter              :: EQ_ENERGY_TOT2=4

  integer, public, protected              :: twofl_eq_energy = EQ_ENERGY_TOT

  !> Whether hyperdiffusivity is used
  logical, public, protected              :: twofl_hyperdiffusivity = .false.
  logical, public, protected              :: twofl_dump_hyperdiffusivity_coef = .false.
  double precision, public, protected, allocatable :: c_shk(:)
  double precision, public, protected, allocatable :: c_hyp(:)

  !> Whether thermal conduction is used
  logical, public, protected              :: twofl_thermal_conduction_c = .false.
  !> type of TC used: 1: adapted module (mhd implementation), 2: adapted module (hd implementation)
  integer, parameter, private             :: MHD_TC =1
  integer, parameter, private             :: HD_TC =2
  integer, protected                      :: use_twofl_tc_c = MHD_TC

  !> Whether radiative cooling is added
  logical, public, protected              :: twofl_radiative_cooling_c = .false.
  type(rc_fluid), allocatable :: rc_fl_c

  !> Whether viscosity is added
  logical, public, protected              :: twofl_viscosity = .false.

  !> Whether gravity is added: common flag for charges and neutrals
  logical, public, protected              :: twofl_gravity = .false.

  !> whether dump full variables (when splitting is used) in a separate dat file
  logical, public, protected              :: twofl_dump_full_vars = .false.

  !> Whether Hall-MHD is used
  logical, public, protected              :: twofl_Hall = .false.

  type(tc_fluid), public, allocatable :: tc_fl_c
  type(te_fluid), public, allocatable :: te_fl_c

  type(tc_fluid), allocatable :: tc_fl_n
  logical, public, protected              :: twofl_thermal_conduction_n = .false.
  logical, public, protected              :: twofl_radiative_cooling_n = .false.
  type(rc_fluid), allocatable :: rc_fl_n

  !> Whether TRAC method is used
  logical, public, protected              :: twofl_trac = .false.

  !> Whether GLM-MHD is used
  logical, public, protected              :: twofl_glm = .false.

  !> Which TRAC method is used            
  integer, public, protected              :: twofl_trac_type=1

  !> Height of the mask used in the TRAC method
  double precision, public, protected              :: twofl_trac_mask = 0.d0

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: twofl_glm_alpha = 0.5d0

  !> MHD fourth order
  logical, public, protected              :: twofl_4th_order = .false.

  !> Index of the density (in the w array)
  integer, public             :: rho_c_

  !> Indices of the momentum density
  integer, allocatable, public :: mom_c(:)

  !> Index of the energy density (-1 if not present)
  integer, public             :: e_c_=-1

  !> Index of the cutoff temperature for the TRAC method
  integer, public              :: Tcoff_c_
  integer, public              :: Tweight_c_

  !> Indices of the GLM psi
  integer, public, protected :: psi_

  !> Indices of auxiliary internal energy
  integer, public :: eaux_c_

  !> Indices of the magnetic field
  integer, allocatable, public :: mag(:)

  !> equi vars flags
  logical, public :: has_equi_rho_c0 = .false.  
  logical, public :: has_equi_pe_c0 = .false.  

  !> equi vars indices in the state%equi_vars array
  integer, public :: equi_rho_c0_ = -1
  integer, public :: equi_pe_c0_ = -1
  logical, public                         :: twofl_equi_thermal_c = .false.

  !neutrals:

  integer, public              :: rho_n_
  integer, allocatable, public :: mom_n(:)
  integer, public              :: e_n_
  integer, public              :: Tcoff_n_
  integer, public              :: Tweight_n_
  logical, public :: has_equi_rho_n0 = .false. 
  logical, public :: has_equi_pe_n0 = .false.  
  integer, public :: equi_rho_n0_ = -1
  integer, public :: equi_pe_n0_ = -1

  ! related to collisions:
  !> collisional alpha
  double precision, public                :: twofl_alpha_coll = 0d0
  logical, public                         :: twofl_alpha_coll_constant = .true.
  !> whether include thermal exchange collisional terms
  logical, public                         :: twofl_coll_inc_te = .true.
  !> whether include ionization/recombination inelastic collisional terms
  logical, public                         :: twofl_coll_inc_ionrec = .false.
  logical, public                         :: twofl_equi_thermal = .true.
  logical, public                         :: twofl_equi_ionrec = .false.
  logical, public                         :: twofl_equi_thermal_n = .false.
  double precision, public                :: dtcollpar = -1d0 !negative value does not impose restriction on the timestep
  !> whether dump collisional terms in a separte dat file
  logical, public, protected              :: twofl_dump_coll_terms = .false.

  ! TODO Helium abundance not used, radiative cooling init uses it
  ! not in parameters list anymore 
  double precision, public, protected  :: He_abundance = 0d0
  ! two fluid is only H plasma
  double precision, public, protected  :: Rc = 2d0
  double precision, public, protected  :: Rn = 1d0

  !> The adiabatic index
  double precision, public                :: twofl_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: twofl_adiab = 1.0d0

  !> The MHD resistivity
  double precision, public                :: twofl_eta = 0.0d0

  !> The MHD hyper-resistivity
  double precision, public                :: twofl_eta_hyper = 0.0d0

  !> The MHD Hall coefficient
  double precision, public                :: twofl_etah = 0.0d0

  !> The small_est allowed energy
  double precision, protected             :: small_e

  !> Method type to clean divergence of B
  character(len=std_len), public, protected :: typedivbfix  = 'linde'

  !> Method type of constrained transport
  character(len=std_len), public, protected :: type_ct  = 'uct_contact'

  !> Whether divB is computed with a fourth order approximation
  logical, public, protected :: twofl_divb_4thorder = .false.

  !> Method type in a integer for good performance
  integer :: type_divb

  !> Coefficient of diffusive divB cleaning
  double precision :: divbdiff     = 0.8d0

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'

  !> clean initial divB
  logical, public :: clean_initial_divb     = .false.

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*^ND)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public, protected :: boundary_divbfix_skip(2*^ND)=0

  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.

  logical :: twofl_cbounds_species = .true.

  !> added from modules: gravity
  !> source split or not
  logical :: grav_split= .false.

  !> gamma minus one and its inverse
  double precision :: gamma_1, inv_gamma_1

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

  ! Public methods
  public :: twofl_phys_init
  public :: twofl_to_conserved
  public :: twofl_to_primitive
  public :: get_divb
  public :: get_rhoc_tot
  public :: twofl_get_v_c_idim
  ! TODO needed for the roe, see if can be used for n
  public :: twofl_get_csound2_c_from_conserved
  public :: get_rhon_tot
  public :: get_alpha_coll_plasma
  public :: get_gamma_ion_rec
  public :: twofl_get_v_n_idim
  public :: get_current
  public :: twofl_get_pthermal_c
  public :: twofl_face_to_center
  public :: get_normalized_divb
  public :: b_from_vector_potential
  {^NOONED
  public :: twofl_clean_divb_multigrid
  }

  abstract interface

    subroutine implicit_mult_factor_subroutine(ixI^L, ixO^L, step_dt, JJ, res)
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in) :: step_dt
    double precision, intent(in) :: JJ(ixI^S)
    double precision, intent(out) :: res(ixI^S)

  end subroutine implicit_mult_factor_subroutine

  end interface

   procedure (implicit_mult_factor_subroutine), pointer :: calc_mult_factor => null()
   integer, protected ::  twofl_implicit_calc_mult_method = 1

contains

  !> Read this module"s parameters from a file
  subroutine twofl_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /twofl_list/ twofl_eq_energy, twofl_gamma, twofl_adiab,&
      twofl_eta, twofl_eta_hyper, twofl_etah, twofl_glm_alpha,& 
      twofl_thermal_conduction_c, use_twofl_tc_c, twofl_radiative_cooling_c, twofl_Hall, twofl_gravity,&
      twofl_viscosity, twofl_4th_order, typedivbfix, source_split_divb, divbdiff,&
      typedivbdiff, type_ct, divbwave, SI_unit, B0field,&
      B0field_forcefree, Bdip, Bquad, Boct, Busr,twofl_equi_thermal_c,&
      twofl_dump_full_vars, has_equi_rho_c0, has_equi_pe_c0, twofl_hyperdiffusivity,twofl_dump_hyperdiffusivity_coef,&
      has_equi_pe_n0, has_equi_rho_n0, twofl_thermal_conduction_n, twofl_radiative_cooling_n,  &
      twofl_alpha_coll,twofl_alpha_coll_constant,&
      twofl_coll_inc_te, twofl_coll_inc_ionrec,twofl_equi_ionrec,twofl_equi_thermal,&
      twofl_equi_thermal_n,dtcollpar,&
      twofl_dump_coll_terms,twofl_implicit_calc_mult_method,&
      boundary_divbfix, boundary_divbfix_skip, twofl_divb_4thorder, &
      clean_initial_divb,  &
      twofl_trac, twofl_trac_type, twofl_trac_mask,twofl_cbounds_species 

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, twofl_list, end=111)
111    close(unitpar)
    end do
 
  end subroutine twofl_read_params

  subroutine twofl_init_hyper(files)
    use mod_global_parameters
    use mod_hyperdiffusivity, only: hyperdiffusivity_init
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /hyperdiffusivity_list/ c_shk, c_hyp

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hyperdiffusivity_list, end=113)
113    close(unitpar)
    end do
    
    call hyperdiffusivity_init()

    !!DEBUG
    if(mype .eq. 0) then
      print*, "Using Hyperdiffusivity"
      print*, "C_SHK ", c_shk(:)
      print*, "C_HYP ", c_hyp(:)
    endif

  end subroutine twofl_init_hyper

  !> Write this module's parameters to a snapsoht
  subroutine twofl_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = twofl_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine twofl_write_info

  subroutine twofl_angmomfix(fC,x,wnew,ixI^L,ixO^L,idim)
    use mod_global_parameters
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(inout)    :: fC(ixI^S,1:nwflux,1:ndim),  wnew(ixI^S,1:nw)
    integer, intent(in)                :: ixI^L, ixO^L
    integer, intent(in)                :: idim
    integer                            :: hxO^L, kxC^L, iw
    double precision                   :: inv_volume(ixI^S)

    call mpistop("to do")

  end subroutine twofl_angmomfix

  subroutine twofl_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    !use mod_gravity, only: gravity_init
    use mod_supertimestepping, only: sts_init, add_sts_method,&
            set_conversion_methods_to_head, set_error_handling_to_head
    {^NOONED
    use mod_multigrid_coupling
    }
    integer :: itr, idir

    call twofl_read_params(par_files)
    physics_type = "twofl"
    if (twofl_cbounds_species) then
      number_species = 2
    endif
    phys_energy=.true.
  !> Solve total energy equation or not
  ! for the two fluid the true value means 
  ! E_charges = E_mag + E_kin_charges + E_int_charges
  ! E_neutrals =  E_kin_neutrals + E_int_neutrals
    phys_total_energy=.false.

  !> Solve internal enery instead of total energy
  ! for the two fluid the true vale means 
  ! E_charges = E_int_charges
  ! E_neutrals = E_int_neutrals
    phys_internal_e=.false.

  ! For the two fluid phys_energy=.true. and phys_internal_e=.false. and phys_total_energy = .false. means
  ! E_charges = E_kin_charges + E_int_charges
  ! E_neutrals =  E_kin_neutrals + E_int_neutrals
    phys_gamma = twofl_gamma

  !> Solve internal energy and total energy equations
  ! this implies two equations of energy solved
    phys_solve_eaux=.false.

    if(twofl_eq_energy == EQ_ENERGY_INT) then
      phys_internal_e = .true.
    elseif(twofl_eq_energy == EQ_ENERGY_TOT .or. twofl_eq_energy == EQ_ENERGY_TOT2) then
      phys_total_energy = .true.
      if(twofl_eq_energy == EQ_ENERGY_TOT2) then
        phys_solve_eaux = .true.
      endif
    elseif(twofl_eq_energy == EQ_ENERGY_NONE) then
      phys_energy = .false.
    endif

    phys_trac=twofl_trac
    phys_trac_type=twofl_trac_type

    if(.not. phys_energy) then
      if(twofl_thermal_conduction_n) then
        twofl_thermal_conduction_n=.false.
        if(mype==0) write(*,*) 'WARNING: set twofl_thermal_conduction_n=F when twofl_energy=F'
      end if
      if(twofl_radiative_cooling_n) then
        twofl_radiative_cooling_n=.false.
        if(mype==0) write(*,*) 'WARNING: set twofl_radiative_cooling_n=F when twofl_energy=F'
      end if
      if(twofl_thermal_conduction_c) then
        twofl_thermal_conduction_c=.false.
        if(mype==0) write(*,*) 'WARNING: set twofl_thermal_conduction_c=F when twofl_energy=F'
      end if
      if(twofl_radiative_cooling_c) then
        twofl_radiative_cooling_c=.false.
        if(mype==0) write(*,*) 'WARNING: set twofl_radiative_cooling_c=F when twofl_energy=F'
      end if
      if(twofl_trac) then
        twofl_trac=.false.
        if(mype==0) write(*,*) 'WARNING: set twofl_trac=F when twofl_energy=F'
      end if
    end if
    {^IFONED
      if(twofl_trac .and. twofl_trac_type .gt. 1) then
        twofl_trac_type=1
        if(mype==0) write(*,*) 'WARNING: set twofl_trac_type=1 for 1D simulation'
      end if
    }
    if(twofl_trac .and. twofl_trac_type .le. 3) then
      twofl_trac_mask=bigdouble
      if(mype==0) write(*,*) 'WARNING: set twofl_trac_mask==bigdouble for global TRAC method'
    end if
    phys_trac_mask=twofl_trac_mask

    if(phys_solve_eaux) prolongprimitive=.true.

    ! set default gamma for polytropic/isothermal process
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
       type_divb = divb_none
    {^NOONED
    case ('multigrid')
       type_divb = divb_multigrid
       use_multigrid = .true.
       mg%operator_type = mg_laplacian
       phys_global_source_after => twofl_clean_divb_multigrid
    }
    case ('glm')
      twofl_glm          = .true.
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
      twofl_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_lindeglm
    case ('ct')
      type_divb = divb_ct
      stagger_grid = .true.
    case default
      call mpistop('Unknown divB fix')
    end select

    allocate(start_indices(number_species))
    allocate(stop_indices(number_species))
    start_indices(1)=1
    !allocate charges first and the same order as in mhd module
    rho_c_ = var_set_fluxvar("rho_c", "rho_c")
    !set variables from mod_variables to point to charges vars
    iw_rho = rho_c_

    allocate(mom_c(ndir))
    do idir=1,ndir
      mom_c(idir) = var_set_fluxvar("m_c","v_c",idir)
    enddo

    allocate(iw_mom(ndir))
    iw_mom(1:ndir) = mom_c(1:ndir)

    ! Set index of energy variable
    if (phys_energy) then
      e_c_ = var_set_fluxvar("e_c", "p_c")
      iw_e = e_c_
    else
      e_c_ = -1
    end if

  ! ambipolar sts assumes mag and energy charges are continuous
    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)

    if (twofl_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    !  set auxiliary internal energy variable
    if(phys_energy .and. phys_solve_eaux) then
      eaux_c_ = var_set_fluxvar("eaux_c", "paux_c",need_bc=.false.)
      iw_eaux = eaux_c_
    else
      eaux_c_ = -1
    end if

    ! set cutoff temperature when using the TRAC method, as well as an auxiliary weight
    Tweight_c_ = -1
    if(twofl_trac) then
      Tcoff_c_ = var_set_wextra()
      iw_tcoff = Tcoff_c_
      if(twofl_trac_type > 2) then
        Tweight_c_ = var_set_wextra()
      endif
    else
      Tcoff_c_ = -1
    end if

  !now allocate neutrals

    ! TODO so far number_species is only used to treat them differently
    ! in the solvers (different cbounds)
    if (twofl_cbounds_species) then
      stop_indices(1)=nwflux
      start_indices(2)=nwflux+1
    endif

    ! Determine flux variables
    rho_n_ = var_set_fluxvar("rho_n", "rho_n")
    allocate(mom_n(ndir))
    do idir=1,ndir
      mom_n(idir) = var_set_fluxvar("m_n","v_n",idir)
    enddo
    if (phys_energy) then
      e_n_ = var_set_fluxvar("e_n", "p_n")
    else
      e_n_     = -1
    end if

    Tweight_n_ = -1
    if(twofl_trac) then
      Tcoff_n_ = var_set_wextra()
      if(twofl_trac_type > 2) then
        Tweight_n_ = var_set_wextra()
      endif
    else
      Tcoff_n_ = -1
    end if

    stop_indices(number_species)=nwflux

    ! set indices of equi vars and update number_equi_vars
    number_equi_vars = 0
    if(has_equi_rho_n0) then
      number_equi_vars = number_equi_vars + 1
      equi_rho_n0_ = number_equi_vars
    endif  
    if(has_equi_pe_n0) then
      number_equi_vars = number_equi_vars + 1
      equi_pe_n0_ = number_equi_vars
    endif  
    if(has_equi_rho_c0) then
      number_equi_vars = number_equi_vars + 1
      equi_rho_c0_ = number_equi_vars
      iw_equi_rho = equi_rho_c0_
    endif  
    if(has_equi_pe_c0) then
      number_equi_vars = number_equi_vars + 1
      equi_pe_c0_ = number_equi_vars
      iw_equi_p = equi_pe_c0_
    endif  

    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! determine number of stagger variables
    if(stagger_grid) nws=ndim

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    if(ndim>1) then
      if(twofl_glm) then
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

    phys_get_dt              => twofl_get_dt
    phys_get_cmax            => twofl_get_cmax
    phys_get_a2max           => twofl_get_a2max
    !phys_get_tcutoff         => twofl_get_tcutoff_c
    if(twofl_cbounds_species) then
      if (mype .eq. 0) print*, "Using different cbounds for each species nspecies = ", number_species
      phys_get_cbounds         => twofl_get_cbounds_species
      phys_get_H_speed         => twofl_get_H_speed_species
    else
      if (mype .eq. 0) print*, "Using same cbounds for all species"
      phys_get_cbounds         => twofl_get_cbounds_one
      phys_get_H_speed         => twofl_get_H_speed_one
    endif
    phys_get_flux            => twofl_get_flux
    phys_add_source_geom     => twofl_add_source_geom
    phys_add_source          => twofl_add_source
    phys_to_conserved        => twofl_to_conserved
    phys_to_primitive        => twofl_to_primitive
    phys_check_params        => twofl_check_params
    phys_check_w             => twofl_check_w
    phys_write_info          => twofl_write_info
    phys_angmomfix           => twofl_angmomfix
    phys_handle_small_values => twofl_handle_small_values
    phys_energy_synchro      => twofl_energy_synchro
    !set equilibrium variables for the new grid
    if(number_equi_vars>0) then
      phys_set_equi_vars => set_equi_vars_grid
    endif
    ! convert_type is not known here, so associate the corresp. subroutine in check_params
    if(type_divb==divb_glm) then
      phys_modify_wLR => twofl_modify_wLR
    end if

    ! if using ct stagger grid, boundary divb=0 is not done here
    if(stagger_grid) then
      phys_get_ct_velocity => twofl_get_ct_velocity
      phys_update_faces => twofl_update_faces
      phys_face_to_center => twofl_face_to_center
      phys_modify_wLR => twofl_modify_wLR
    else if(ndim>1) then
      phys_boundary_adjust => twofl_boundary_adjust
    end if

    {^NOONED
    ! clean initial divb
    if(clean_initial_divb) phys_clean_divb => twofl_clean_divb_multigrid
    }

    ! Whether diagonal ghost cells are required for the physics
    if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! derive units from basic units
    call twofl_physical_units()

    if(.not. phys_energy .and. (twofl_thermal_conduction_c& 
        .or. twofl_thermal_conduction_n)) then
      call mpistop("thermal conduction needs twofl_energy=T")
    end if

    ! initialize thermal conduction module
    if (twofl_thermal_conduction_c &
        .or. twofl_thermal_conduction_n) then
      phys_req_diagonal = .true.
      call sts_init()
      call tc_init_params(twofl_gamma)
    endif
    if (twofl_thermal_conduction_c) then
      allocate(tc_fl_c)
      if(has_equi_pe_c0 .and. has_equi_rho_c0) then
        tc_fl_c%get_temperature_from_eint => twofl_get_temperature_from_eint_c_with_equi
        if(phys_internal_e) then
          tc_fl_c%get_temperature_from_conserved => twofl_get_temperature_from_eint_c_with_equi
        else
            if(twofl_eq_energy == EQ_ENERGY_KI) then
              tc_fl_c%get_temperature_from_conserved => twofl_get_temperature_from_eki_c_with_equi
            else
              tc_fl_c%get_temperature_from_conserved => twofl_get_temperature_from_etot_c_with_equi
            endif
        endif
        if(twofl_equi_thermal_c) then
          tc_fl_c%has_equi = .true.
          tc_fl_c%get_temperature_equi => twofl_get_temperature_c_equi
          tc_fl_c%get_rho_equi => twofl_get_rho_c_equi
        else  
          tc_fl_c%has_equi = .false.
        endif
      else
        if(phys_internal_e) then
          tc_fl_c%get_temperature_from_conserved => twofl_get_temperature_from_eint_c
         else
            if(twofl_eq_energy == EQ_ENERGY_KI) then
              tc_fl_c%get_temperature_from_conserved => twofl_get_temperature_from_eki_c
            else
              tc_fl_c%get_temperature_from_conserved => twofl_get_temperature_from_etot_c
          endif  
         endif  
        tc_fl_c%get_temperature_from_eint => twofl_get_temperature_from_eint_c
      endif
      if(use_twofl_tc_c .eq. MHD_TC) then
        call tc_get_mhd_params(tc_fl_c,tc_c_params_read_mhd)
        call add_sts_method(twofl_get_tc_dt_mhd_c,twofl_sts_set_source_tc_c_mhd,e_c_,1,e_c_,1,.false.)
      else if(use_twofl_tc_c .eq. HD_TC) then
        call tc_get_hd_params(tc_fl_c,tc_c_params_read_hd)
        call add_sts_method(twofl_get_tc_dt_hd_c,twofl_sts_set_source_tc_c_hd,e_c_,1,e_c_,1,.false.)
      endif
      if(.not. phys_internal_e) then
        call set_conversion_methods_to_head(twofl_e_to_ei_c, twofl_ei_to_e_c)
      endif
      call set_error_handling_to_head(twofl_tc_handle_small_e_c)
      tc_fl_c%get_rho => get_rhoc_tot
      tc_fl_c%e_ = e_c_
      tc_fl_c%Tcoff_ = Tcoff_c_
    end if
    if (twofl_thermal_conduction_n) then
      allocate(tc_fl_n)
      call tc_get_hd_params(tc_fl_n,tc_n_params_read_hd)
      if(has_equi_pe_n0 .and. has_equi_rho_n0) then
        tc_fl_n%get_temperature_from_eint => twofl_get_temperature_from_eint_n_with_equi
        if(twofl_equi_thermal_n) then
          tc_fl_n%has_equi = .true.
          tc_fl_n%get_temperature_equi => twofl_get_temperature_n_equi
          tc_fl_n%get_rho_equi => twofl_get_rho_n_equi
        else  
          tc_fl_n%has_equi = .false.
        endif
      else
        tc_fl_n%get_temperature_from_eint => twofl_get_temperature_from_eint_n
      endif
      if(phys_internal_e) then
        if(has_equi_pe_n0 .and. has_equi_rho_n0) then
          tc_fl_n%get_temperature_from_conserved => twofl_get_temperature_from_eint_n_with_equi
        else
          tc_fl_n%get_temperature_from_conserved => twofl_get_temperature_from_eint_n
        endif 
        call add_sts_method(twofl_get_tc_dt_hd_n,twofl_sts_set_source_tc_n_hd,e_n_,1,e_n_,1,.false.)
      else
        if(has_equi_pe_n0 .and. has_equi_rho_n0) then
          tc_fl_n%get_temperature_from_conserved => twofl_get_temperature_from_etot_n_with_equi
        else
          tc_fl_n%get_temperature_from_conserved => twofl_get_temperature_from_etot_n
        endif 
        call add_sts_method(twofl_get_tc_dt_hd_n,twofl_sts_set_source_tc_n_hd,e_n_,1,e_n_,1,.false.)
        call set_conversion_methods_to_head(twofl_e_to_ei_n, twofl_ei_to_e_n)
      endif
      call set_error_handling_to_head(twofl_tc_handle_small_e_n)
      tc_fl_n%get_rho => get_rhon_tot
      tc_fl_n%e_ = e_n_
      tc_fl_n%Tcoff_ = Tcoff_n_
    end if


    if(.not. phys_energy .and. (twofl_radiative_cooling_c& 
        .or. twofl_radiative_cooling_n)) then
      call mpistop("radiative cooling needs twofl_energy=T")
    end if

    ! initialize thermal conduction module
    if (twofl_radiative_cooling_c &
        .or. twofl_radiative_cooling_n) then
    ! Initialize radiative cooling module
      call radiative_cooling_init_params(twofl_gamma,He_abundance)
      if(twofl_radiative_cooling_c) then
        allocate(rc_fl_c)
        call radiative_cooling_init(rc_fl_c,rc_params_read_c)
        rc_fl_c%get_rho => get_rhoc_tot
        rc_fl_c%get_pthermal => twofl_get_pthermal_c
        rc_fl_c%Rfactor = Rc
        rc_fl_c%e_ = e_c_
        rc_fl_c%eaux_ = eaux_c_
        rc_fl_c%Tcoff_ = Tcoff_c_
        if(has_equi_pe_c0 .and. has_equi_rho_c0 .and. twofl_equi_thermal_c) then
          rc_fl_c%has_equi = .true.
          rc_fl_c%get_rho_equi => twofl_get_rho_c_equi
          rc_fl_c%get_pthermal_equi => twofl_get_pe_c_equi
        else
          rc_fl_c%has_equi = .false.
        end if
      end if
    end if
    allocate(te_fl_c)
    te_fl_c%get_rho=> get_rhoc_tot
    te_fl_c%get_pthermal=> twofl_get_pthermal_c
    te_fl_c%Rfactor = Rc
{^IFTHREED
    phys_te_images => twofl_te_images
}

    ! Initialize viscosity module
    !!TODO
    !if (twofl_viscosity) call viscosity_init(phys_wider_stencil,phys_req_diagonal)

    ! Initialize gravity module
    if(twofl_gravity) then
    !  call gravity_init()
       call grav_params_read(par_files)
    end if

    ! Initialize particles module
    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in getflux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in nghostcells.
    if (twofl_hall) then
       phys_req_diagonal = .true.
       if (twofl_4th_order) then
          phys_wider_stencil = 2
       else
          phys_wider_stencil = 1
       end if
    end if

    if(twofl_hyperdiffusivity) then
      allocate(c_shk(1:nwflux))
      allocate(c_hyp(1:nwflux))
      call twofl_init_hyper(par_files)
    end if

  end subroutine twofl_phys_init

{^IFTHREED
  subroutine twofl_te_images()
    use mod_global_parameters
    use mod_thermal_emission

    select case(convert_type)
      case('EIvtiCCmpi','EIvtuCCmpi')
        call get_EUV_image(unitconvert,te_fl_c)
      case('ESvtiCCmpi','ESvtuCCmpi')
        call get_EUV_spectrum(unitconvert,te_fl_c)
      case('SIvtiCCmpi','SIvtuCCmpi')
        call get_SXR_image(unitconvert,te_fl_c)
      case default
        call mpistop("Error in synthesize emission: Unknown convert_type")
      end select
  end subroutine twofl_te_images
}  

  ! wrappers for STS functions in thermal_conductivity module
  ! which take as argument the tc_fluid (defined in the physics module)
  subroutine  twofl_sts_set_source_tc_c_mhd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_mhd
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_mhd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,tc_fl_c)
  end subroutine twofl_sts_set_source_tc_c_mhd

  subroutine  twofl_sts_set_source_tc_c_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_hd
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,tc_fl_c)
  end subroutine twofl_sts_set_source_tc_c_hd

  function twofl_get_tc_dt_mhd_c(w,ixI^L,ixO^L,dx^D,x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_mhd
 
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_mhd(w,ixI^L,ixO^L,dx^D,x,tc_fl_c) 
  end function twofl_get_tc_dt_mhd_c

  function twofl_get_tc_dt_hd_c(w,ixI^L,ixO^L,dx^D,x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_hd
 
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_hd(w,ixI^L,ixO^L,dx^D,x,tc_fl_c) 
  end function twofl_get_tc_dt_hd_c

  subroutine twofl_tc_handle_small_e_c(w, x, ixI^L, ixO^L, step)
    use mod_global_parameters
    use mod_small_values

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer, intent(in)    :: step

    character(len=140) :: error_msg

    write(error_msg,"(a,i3)") "Charges thermal conduction step ", step
    call twofl_handle_small_ei_c(w,x,ixI^L,ixO^L,e_c_,error_msg)
  end subroutine twofl_tc_handle_small_e_c

  subroutine  twofl_sts_set_source_tc_n_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_hd
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,tc_fl_n)
  end subroutine twofl_sts_set_source_tc_n_hd

  subroutine twofl_tc_handle_small_e_n(w, x, ixI^L, ixO^L, step)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer, intent(in)    :: step

    character(len=140) :: error_msg

    write(error_msg,"(a,i3)") "Neutral thermal conduction step ", step
    call twofl_handle_small_ei_n(w,x,ixI^L,ixO^L,e_n_,error_msg)
  end subroutine twofl_tc_handle_small_e_n

  function twofl_get_tc_dt_hd_n(w,ixI^L,ixO^L,dx^D,x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_hd
 
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_hd(w,ixI^L,ixO^L,dx^D,x,tc_fl_n) 
  end function twofl_get_tc_dt_hd_n

  subroutine tc_n_params_read_hd(fl)
    use mod_global_parameters, only: unitpar,par_files
    use mod_global_parameters, only: unitpar
    type(tc_fluid), intent(inout) :: fl
    integer                      :: n
    logical :: tc_saturate=.false.
    double precision :: tc_k_para=0d0

    namelist /tc_n_list/ tc_saturate, tc_k_para

    do n = 1, size(par_files)
       open(unitpar, file=trim(par_files(n)), status="old")
       read(unitpar, tc_n_list, end=111)
111      close(unitpar)
    end do
    fl%tc_saturate = tc_saturate
    fl%tc_k_para = tc_k_para

  end subroutine tc_n_params_read_hd

  subroutine rc_params_read_n(fl)
    use mod_global_parameters, only: unitpar,par_files
    use mod_constants, only: bigdouble
    type(rc_fluid), intent(inout) :: fl
    integer                      :: n
    ! list parameters
    integer :: ncool = 4000
    double precision :: cfrac=0.1d0
  
    !> Name of cooling curve
    character(len=std_len)  :: coolcurve='JCorona'
  
    !> Name of cooling method
    character(len=std_len)  :: coolmethod='exact'
  
    !> Fixed temperature not lower than tlow
    logical    :: Tfix=.false.
  
    !> Lower limit of temperature
    double precision   :: tlow=bigdouble
  
    !> Add cooling source in a split way (.true.) or un-split way (.false.)
    logical    :: rc_split=.false.

    namelist /rc_list_n/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix, rc_split

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, rc_list_n, end=111)
111     close(unitpar)
    end do

    fl%ncool=ncool
    fl%coolcurve=coolcurve
    fl%coolmethod=coolmethod
    fl%tlow=tlow
    fl%Tfix=Tfix
    fl%rc_split=rc_split
    fl%cfrac=cfrac
  end subroutine rc_params_read_n

  !end wrappers

  ! fill in tc_fluid fields from namelist
  subroutine tc_c_params_read_mhd(fl)
    use mod_global_parameters, only: unitpar,par_files
    type(tc_fluid), intent(inout) :: fl

    integer                      :: n

    ! list parameters
    logical :: tc_perpendicular=.true.
    logical :: tc_saturate=.false.
    double precision :: tc_k_para=0d0
    double precision :: tc_k_perp=0d0
    character(len=std_len)  :: tc_slope_limiter="MC"

    namelist /tc_c_list/ tc_perpendicular, tc_saturate, tc_slope_limiter, tc_k_para, tc_k_perp
    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, tc_c_list, end=111)
111     close(unitpar)
    end do

    fl%tc_perpendicular = tc_perpendicular
    fl%tc_saturate = tc_saturate
    fl%tc_k_para = tc_k_para
    fl%tc_k_perp = tc_k_perp
    fl%tc_slope_limiter = tc_slope_limiter
  end subroutine tc_c_params_read_mhd

  subroutine tc_c_params_read_hd(fl)
    use mod_global_parameters, only: unitpar,par_files
    use mod_global_parameters, only: unitpar
    type(tc_fluid), intent(inout) :: fl
    integer                      :: n
    logical :: tc_saturate=.false.
    double precision :: tc_k_para=0d0

    namelist /tc_c_list/ tc_saturate, tc_k_para

    do n = 1, size(par_files)
       open(unitpar, file=trim(par_files(n)), status="old")
       read(unitpar, tc_c_list, end=111)
111      close(unitpar)
    end do
    fl%tc_saturate = tc_saturate
    fl%tc_k_para = tc_k_para

  end subroutine tc_c_params_read_hd

!! end th cond

!!rad cool
  subroutine rc_params_read_c(fl)
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


    namelist /rc_list_c/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix, rc_split

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, rc_list_c, end=111)
111     close(unitpar)
    end do

    fl%ncool=ncool
    fl%coolcurve=coolcurve
    fl%coolmethod=coolmethod
    fl%tlow=tlow
    fl%Tfix=Tfix
    fl%rc_split=rc_split
    fl%cfrac=cfrac
  end subroutine rc_params_read_c

!! end rad cool

  !> sets the equilibrium variables
  subroutine set_equi_vars_grid_faces(igrid,x,ixI^L,ixO^L)
    use mod_global_parameters
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

  ! w, wnew conserved
  function convert_vars_splitting(ixI^L,ixO^L, w, x, nwc) result(wnew)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L, nwc
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim) 
    double precision   :: wnew(ixO^S, 1:nwc)
    double precision   :: rho(ixI^S)

    call  get_rhon_tot(w,x,ixI^L,ixO^L,rho(ixI^S))
    wnew(ixO^S,rho_n_) = rho(ixO^S)
    wnew(ixO^S,mom_n(:)) =  w(ixO^S,mom_n(:))
    call  get_rhoc_tot(w,x,ixI^L,ixO^L,rho(ixI^S))
    wnew(ixO^S,rho_c_) = rho(ixO^S)
    wnew(ixO^S,mom_c(:)) =  w(ixO^S,mom_c(:))

    if (B0field) then
      ! add background magnetic field B0 to B
      wnew(ixO^S,mag(:))=w(ixO^S,mag(:))+block%B0(ixO^S,:,0)
    else 
      wnew(ixO^S,mag(:))=w(ixO^S,mag(:))
    end if

    if(phys_energy) then
      wnew(ixO^S,e_n_) = w(ixO^S,e_n_)
      if(has_equi_pe_n0) then
        wnew(ixO^S,e_n_) = wnew(ixO^S,e_n_) + block%equi_vars(ixO^S,equi_pe_n0_,0)* inv_gamma_1  
      endif
      wnew(ixO^S,e_c_) = w(ixO^S,e_c_)
      if(has_equi_pe_c0) then
        wnew(ixO^S,e_c_) = wnew(ixO^S,e_c_) + block%equi_vars(ixO^S,equi_pe_c0_,0)* inv_gamma_1  
      endif
      if(B0field .and. phys_total_energy) then
          wnew(ixO^S,e_c_)=wnew(ixO^S,e_c_)+0.5d0*sum(block%B0(ixO^S,:,0)**2,dim=ndim+1) &
              + sum(w(ixO^S,mag(:))*block%B0(ixO^S,:,0),dim=ndim+1)
      endif
    endif

  end function convert_vars_splitting

  !> copied from mod_gravity
  subroutine grav_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /grav_list/ grav_split

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grav_list, end=111)
111    close(unitpar)
    end do

  end subroutine grav_params_read

  subroutine associate_dump_hyper()
    use mod_global_parameters
    use mod_convert, only: add_convert_method
    integer :: ii
    do ii = 1,ndim 
      if(ii==1) then
        call add_convert_method(dump_hyperdiffusivity_coef_x, nw, cons_wnames(1:nw), "hyper_x")
      elseif(ii==2) then
        call add_convert_method(dump_hyperdiffusivity_coef_y, nw, cons_wnames(1:nw), "hyper_y")
      else
        call add_convert_method(dump_hyperdiffusivity_coef_z, nw, cons_wnames(1:nw), "hyper_z")
      endif
    enddo
  end subroutine associate_dump_hyper

  subroutine twofl_check_params
    use mod_global_parameters
    use mod_usr_methods
    use mod_convert, only: add_convert_method

    ! after user parameter setting
    gamma_1=twofl_gamma-1.d0
    if (.not. phys_energy) then
       if (twofl_gamma <= 0.0d0) call mpistop ("Error: twofl_gamma <= 0")
       if (twofl_adiab < 0.0d0) call mpistop ("Error: twofl_adiab < 0")
       small_pressure = twofl_adiab*small_density**twofl_gamma
    else
       if (twofl_gamma <= 0.0d0 .or. twofl_gamma == 1.0d0) &
            call mpistop ("Error: twofl_gamma <= 0 or twofl_gamma == 1")
       inv_gamma_1=1.d0/gamma_1
       small_e = small_pressure * inv_gamma_1
    end if

    ! this has to be done here as use_imex_scheme is not set in init subroutine, 
    ! but here it is
    if(use_imex_scheme) then
      if(has_collisions()) then
        ! implicit collisional terms update
        phys_implicit_update => twofl_implicit_coll_terms_update
        phys_evaluate_implicit => twofl_evaluate_implicit
        if(mype .eq. 1) then
            print*, "IMPLICIT UPDATE with calc_mult_factor", twofl_implicit_calc_mult_method
        endif
        if(twofl_implicit_calc_mult_method == 1) then
          calc_mult_factor => calc_mult_factor1
        else
          calc_mult_factor => calc_mult_factor2
        endif
      endif
    else
      ! check dtcoll par for explicit implementation of the coll. terms
      if(dtcollpar .le. 0d0 .or. dtcollpar .ge. 1d0) then
        if (mype .eq. 0) print*, "Explicit update of coll terms requires 0<dtcollpar<1, dtcollpar set to 0.8."
        dtcollpar = 0.8
      endif 
        
    endif
!    if(H_ion_fr == 0d0 .and. He_ion_fr == 0d0) then
!      call mpistop("H_ion_fr or He_ion_fr must be > 0 or use hd module")
!    endif
!    if(H_ion_fr == 1d0 .and. He_ion_fr == 1d0) then
!      call mpistop("H_ion_fr or He_ion_fr must be < 1 or use mhd module")
!    endif
    if (number_equi_vars > 0 .and. .not. associated(usr_set_equi_vars)) then
      call mpistop("usr_set_equi_vars has to be implemented in the user file")
    endif
    if(convert .or. autoconvert) then
      if(convert_type .eq. 'dat_generic_mpi') then
        if(twofl_dump_full_vars) then
          if(mype .eq. 0) print*, " add conversion method: split -> full "
          call add_convert_method(convert_vars_splitting, nw, cons_wnames, "new")
        endif
        if(twofl_dump_coll_terms) then
          if(mype .eq. 0) print*, " add conversion method: dump coll terms "
          call add_convert_method(dump_coll_terms, 3, (/"alpha    ", "gamma_rec", "gamma_ion"/), "_coll")
        endif
        if(twofl_hyperdiffusivity .and. twofl_dump_hyperdiffusivity_coef) then
          if(mype .eq. 0) print*, " add conversion method: dump hyperdiffusivity coeff. "
          call associate_dump_hyper()
        endif
      endif
    endif
  end subroutine twofl_check_params

  subroutine twofl_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0,c_lightspeed
    !double precision :: a,b,c,d
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
      miu0=4.d0*dpi
      c_lightspeed=const_c
    end if


    a=1d0  
    b=1d0
    Rc=2d0
    Rn=1d0  

    !now the unit choice:
    !unit 1 from number density or density -> mH
    !unit 2 from 

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
  end subroutine twofl_physical_units

  subroutine twofl_check_w(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision :: tmp(ixI^S)
    logical, intent(inout) :: flag(ixI^S,1:nw)

    flag=.false.
        
    if(has_equi_rho_n0) then
      tmp(ixO^S) = w(ixO^S,rho_n_) + block%equi_vars(ixO^S,equi_rho_n0_,0)
    else  
      tmp(ixO^S) = w(ixO^S,rho_n_) 
    endif
    where(tmp(ixO^S) < small_density) flag(ixO^S,rho_n_) = .true.
    if(has_equi_rho_c0) then
      tmp(ixO^S) = w(ixO^S,rho_c_) + block%equi_vars(ixO^S,equi_rho_c0_,0)
    else  
      tmp(ixO^S) = w(ixO^S,rho_c_) 
    endif
    where(tmp(ixO^S) < small_density) flag(ixO^S,rho_c_) = .true.
    if(phys_energy) then
      if(primitive) then
        tmp(ixO^S) = w(ixO^S,e_n_)
        if(has_equi_pe_n0) then
          tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe_n0_,0)
        endif
        where(tmp(ixO^S) < small_pressure) flag(ixO^S,e_n_) = .true.
        tmp(ixO^S) = w(ixO^S,e_c_)
        if(has_equi_pe_c0) then
          tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe_c0_,0)
        endif
        where(tmp(ixO^S) < small_pressure) flag(ixO^S,e_c_) = .true.
        ! TODO , also in mhd?
        !if(twofl_eq_energy == EQ_ENERGY_TOT2) then 
        !  where(w(ixO^S,eaux_c_) < small_pressure) flag(ixO^S,eaux_c_) = .true.
        !endif
      else
        if(phys_internal_e) then
          tmp(ixO^S)=w(ixO^S,e_n_)
          if(has_equi_pe_n0) then
            tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe_n0_,0)*inv_gamma_1
          endif
          where(tmp(ixO^S) < small_e) flag(ixO^S,e_n_) = .true.
          tmp(ixO^S)=w(ixO^S,e_c_)
          if(has_equi_pe_c0) then
            tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe_c0_,0)*inv_gamma_1
          endif
          where(tmp(ixO^S) < small_e) flag(ixO^S,e_c_) = .true.
        else
          !neutrals
          tmp(ixO^S)=w(ixO^S,e_n_)-&
                twofl_kin_en_n(w,ixI^L,ixO^L)
          if(has_equi_pe_n0) then
            tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe_n0_,0)*inv_gamma_1
          endif
          where(tmp(ixO^S) < small_e) flag(ixO^S,e_n_) = .true.
          if(phys_total_energy) then
            tmp(ixO^S)=w(ixO^S,e_c_)-&
                twofl_kin_en_c(w,ixI^L,ixO^L)-twofl_mag_en(w,ixI^L,ixO^L)
          else
            tmp(ixO^S)=w(ixO^S,e_c_)-&
                twofl_kin_en_c(w,ixI^L,ixO^L)
          end if
          if(has_equi_pe_c0) then
            tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe_c0_,0)*inv_gamma_1
          endif
          where(tmp(ixO^S) < small_e) flag(ixO^S,e_c_) = .true.
          if(twofl_eq_energy == EQ_ENERGY_TOT2) then 
            tmp(ixO^S)=w(ixO^S,eaux_c_)
            if(has_equi_pe_c0) then
              tmp(ixO^S) = tmp(ixO^S)+block%equi_vars(ixO^S,equi_pe_c0_,0)*inv_gamma_1
            endif
            where(tmp(ixO^S) < small_e) flag(ixO^S,e_c_) = .true.
          endif
        end if
      endif
    end if

  end subroutine twofl_check_w

  !> Transform primitive variables into conservative ones
  subroutine twofl_to_conserved(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: idir
    double precision                :: rhoc(ixI^S)
    double precision                :: rhon(ixI^S)

    !if (fix_small_values) then
    !  call twofl_handle_small_values(.true., w, x, ixI^L, ixO^L, 'twofl_to_conserved')
    !end if

    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)

    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(phys_energy) then
      if(phys_internal_e) then
        w(ixO^S,e_n_)=w(ixO^S,e_n_)*inv_gamma_1
        w(ixO^S,e_c_)=w(ixO^S,e_c_)*inv_gamma_1
      else
        w(ixO^S,e_n_)=w(ixO^S,e_n_)*inv_gamma_1&
                   +half*sum(w(ixO^S,mom_n(:))**2,dim=ndim+1)*rhon(ixO^S)
        if(phys_total_energy) then
          w(ixO^S,e_c_)=w(ixO^S,e_c_)*inv_gamma_1&
                   +half*sum(w(ixO^S,mom_c(:))**2,dim=ndim+1)*rhoc(ixO^S)&
                   +twofl_mag_en(w, ixI^L, ixO^L)
          if(twofl_eq_energy == EQ_ENERGY_TOT2) then
            w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)*inv_gamma_1
          endif
        else
          ! kinetic energy + internal energy is evolved   
          w(ixO^S,e_c_)=w(ixO^S,e_c_)*inv_gamma_1&
                   +half*sum(w(ixO^S,mom_c(:))**2,dim=ndim+1)*rhoc(ixO^S)
        endif
      end if
      !print*, "TOCONS ec ", w(1:10,e_c_)
      !print*, "TOCONS en ", w(1:10,e_n_)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixO^S, mom_n(idir)) = rhon(ixO^S) * w(ixO^S, mom_n(idir))
       w(ixO^S, mom_c(idir)) = rhoc(ixO^S) * w(ixO^S, mom_c(idir))
    end do
  end subroutine twofl_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine twofl_to_primitive(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: idir
    double precision                :: rhoc(ixI^S)
    double precision                :: rhon(ixI^S)

    if (fix_small_values) then
      call twofl_handle_small_values(.false., w, x, ixI^L, ixO^L, 'twofl_to_primitive')
    end if

    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)

    if(phys_energy) then
      if(phys_internal_e) then
        w(ixO^S,e_n_)=w(ixO^S,e_n_)*gamma_1
        w(ixO^S,e_c_)=w(ixO^S,e_c_)*gamma_1
      else
        ! neutrals evolved energy = ke + e_int 
        w(ixO^S,e_n_)=gamma_1*(w(ixO^S,e_n_)&
                    -twofl_kin_en_n(w,ixI^L,ixO^L))
        ! charges
        if(phys_total_energy) then
         ! evolved energy = ke + e_int + e_mag 
          w(ixO^S,e_c_)=gamma_1*(w(ixO^S,e_c_)&
                    -twofl_kin_en_c(w,ixI^L,ixO^L)&
                    -twofl_mag_en(w,ixI^L,ixO^L))
          if(twofl_eq_energy == EQ_ENERGY_TOT2) then
            w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)*gamma_1
          endif
        else
         ! evolved energy = ke + e_int 
          w(ixO^S,e_c_)=gamma_1*(w(ixO^S,e_c_)&
                    -twofl_kin_en_c(w,ixI^L,ixO^L))
        end if
      end if
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom_c(idir)) = w(ixO^S, mom_c(idir))/rhoc(ixO^S)
       w(ixO^S, mom_n(idir)) = w(ixO^S, mom_n(idir))/rhon(ixO^S)
    end do

  end subroutine twofl_to_primitive

!!USED IN TC
  !> Transform internal energy to total energy
  subroutine twofl_ei_to_e_c(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
 
    ! Calculate total energy from internal, kinetic and magnetic energy
    if(phys_solve_eaux) w(ixI^S,eaux_c_)=w(ixI^S,e_c_)
    if(twofl_eq_energy == EQ_ENERGY_KI) then
      w(ixO^S,e_c_)=w(ixO^S,e_c_)&
                +twofl_kin_en_c(w,ixI^L,ixO^L)
    else
      w(ixO^S,e_c_)=w(ixO^S,e_c_)&
               +twofl_kin_en_c(w,ixI^L,ixO^L)&
               +twofl_mag_en(w,ixI^L,ixO^L)
    endif
  end subroutine twofl_ei_to_e_c

  !> Transform total energy to internal energy
  subroutine twofl_e_to_ei_c(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    if(twofl_eq_energy == EQ_ENERGY_KI) then
      w(ixO^S,e_c_)=w(ixO^S,e_c_)&
                -twofl_kin_en_c(w,ixI^L,ixO^L)
    else
    ! Calculate ei = e - ek - eb
      w(ixO^S,e_c_)=w(ixO^S,e_c_)&
                -twofl_kin_en_c(w,ixI^L,ixO^L)&
                -twofl_mag_en(w,ixI^L,ixO^L)
    endif
  end subroutine twofl_e_to_ei_c

   !Neutrals 
  subroutine twofl_ei_to_e_n(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
 
    ! Calculate total energy from internal and kinetic  energy

    w(ixO^S,e_n_)=w(ixO^S,e_n_)+twofl_kin_en_n(w,ixI^L,ixO^L)

  end subroutine twofl_ei_to_e_n

  !> Transform total energy to internal energy
  subroutine twofl_e_to_ei_n(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    ! Calculate ei = e - ek 
    w(ixO^S,e_n_)=w(ixO^S,e_n_)-twofl_kin_en_n(w,ixI^L,ixO^L)

    call twofl_handle_small_ei_n(w,x,ixI^L,ixO^L,e_n_,"e_to_ei_n")
  end subroutine twofl_e_to_ei_n

  subroutine twofl_energy_synchro(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth1(ixI^S),pth2(ixI^S),alfa(ixI^S),beta(ixI^S)
    double precision, parameter :: beta_low=0.005d0,beta_high=0.05d0

!    double precision :: vtot(ixI^S),cs2(ixI^S),mach(ixI^S)
!    double precision, parameter :: mach_low=20.d0,mach_high=200.d0

    ! get magnetic energy
    alfa(ixO^S)=twofl_mag_en(w,ixI^L,ixO^L)
    pth1(ixO^S)=gamma_1*(w(ixO^S,e_c_)-twofl_kin_en_c(w,ixI^L,ixO^L)-alfa(ixO^S))
    pth2(ixO^S)=w(ixO^S,eaux_c_)*gamma_1
    ! get plasma beta
    beta(ixO^S)=min(pth1(ixO^S),pth2(ixO^S))/alfa(ixO^S)

    ! whether Mach number should be another criterion ?
!    vtot(ixO^S)=sum(w(ixO^S,mom(:))**2,dim=ndim+1)
!    call twofl_get_csound2(w,x,ixI^L,ixO^L,cs2)
!    mach(ixO^S)=sqrt(vtot(ixO^S)/cs2(ixO^S))/w(ixO^S,rho_)
    where(beta(ixO^S) .ge. beta_high)
!    where(beta(ixO^S) .ge. beta_high .and. mach(ixO^S) .le. mach_low)
      w(ixO^S,eaux_c_)=pth1(ixO^S)*inv_gamma_1
    else where(beta(ixO^S) .le. beta_low)
!    else where(beta(ixO^S) .le. beta_low .or. mach(ixO^S) .ge. mach_high)
      w(ixO^S,e_c_)=w(ixO^S,e_c_)-pth1(ixO^S)*inv_gamma_1+w(ixO^S,eaux_c_)
    else where
      alfa(ixO^S)=dlog(beta(ixO^S)/beta_low)/dlog(beta_high/beta_low)
!      alfa(ixO^S)=min(dlog(beta(ixO^S)/beta_low)/dlog(beta_high/beta_low),
!                      dlog(mach_high(ixO^S)/mach(ixO^S))/dlog(mach_high/mach_low))
      w(ixO^S,eaux_c_)=(pth2(ixO^S)*(one-alfa(ixO^S))&
                     +pth1(ixO^S)*alfa(ixO^S))*inv_gamma_1
      w(ixO^S,e_c_)=w(ixO^S,e_c_)-pth1(ixO^S)*inv_gamma_1+w(ixO^S,eaux_c_)
    end where
  end subroutine twofl_energy_synchro

  subroutine twofl_handle_small_values(primitive, w, x, ixI^L, ixO^L, subname)
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
    double precision :: tmp1(ixI^S)

    if(small_values_method == "ignore") return

    call twofl_check_w(primitive, ixI^L, ixO^L, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_rho_c0) then
          where(flag(ixO^S,rho_c_)) w(ixO^S,rho_c_) = &
                    small_density-block%equi_vars(ixO^S,equi_rho_c0_,0)
        else
          where(flag(ixO^S,rho_c_)) w(ixO^S,rho_c_) = small_density
        endif
        if(has_equi_rho_n0) then
          where(flag(ixO^S,rho_n_)) w(ixO^S,rho_n_) = &
                  small_density-block%equi_vars(ixO^S,equi_rho_n0_,0)
        else
          where(flag(ixO^S,rho_n_)) w(ixO^S,rho_n_) = small_density
        endif
        do idir = 1, ndir
          if(small_values_fix_iw(mom_n(idir))) then
            where(flag(ixO^S,rho_n_)) w(ixO^S, mom_n(idir)) = 0.0d0
          end if
          if(small_values_fix_iw(mom_c(idir))) then
            where(flag(ixO^S,rho_c_)) w(ixO^S, mom_c(idir)) = 0.0d0
          end if
        end do

        if(phys_energy) then
          if(primitive) then
           if(has_equi_pe_n0) then 
            tmp1(ixO^S) = small_pressure - &
              block%equi_vars(ixO^S,equi_pe_n0_,0)
           else
            tmp1(ixO^S) = small_pressure
           endif  
           if(has_equi_pe_c0) then 
            tmp2(ixO^S) = small_e - &
              block%equi_vars(ixO^S,equi_pe_c0_,0)
           else
            tmp2(ixO^S) = small_pressure
           endif
          else
            ! conserved  
            if(has_equi_pe_n0) then 
              tmp1(ixO^S) = small_e - &
              block%equi_vars(ixO^S,equi_pe_n0_,0)*inv_gamma_1 
            else
              tmp1(ixO^S) = small_e
            endif  
            if(has_equi_pe_c0) then 
              tmp2(ixO^S) = small_e - &
                block%equi_vars(ixO^S,equi_pe_c0_,0)*inv_gamma_1 
            else
              tmp2(ixO^S) = small_e
            endif  
            if(phys_internal_e) then
              where(flag(ixO^S,e_n_))
                w(ixO^S,e_n_)=tmp1(ixO^S)
              end where
              where(flag(ixO^S,e_c_))
                w(ixO^S,e_c_)=tmp2(ixO^S)
              end where
            else
              where(flag(ixO^S,e_n_))
                w(ixO^S,e_n_) = tmp1(ixO^S)+&
                 twofl_kin_en_n(w,ixI^L,ixO^L)
              end where
              if(phys_total_energy) then
                where(flag(ixO^S,e_c_))
                  w(ixO^S,e_c_) = tmp2(ixO^S)+&
                   twofl_kin_en_c(w,ixI^L,ixO^L)+&
                   twofl_mag_en(w,ixI^L,ixO^L)
                end where
              else
                where(flag(ixO^S,e_c_))
                  w(ixO^S,e_c_) = tmp2(ixO^S)+&
                   twofl_kin_en_c(w,ixI^L,ixO^L)
                end where
              endif
              if(phys_solve_eaux) then
                where(flag(ixO^S,e_c_))
                  w(ixO^S,eaux_c_)=tmp2(ixO^S)
                end where
              end if
            end if
          end if
        end if
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag)
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek-eb)
          if(phys_energy) then
            if(phys_internal_e) then
              w(ixO^S,e_c_)=w(ixO^S,e_c_)*gamma_1
              w(ixO^S,e_n_)=w(ixO^S,e_n_)*gamma_1
            else
             w(ixO^S,e_n_)=gamma_1*(w(ixO^S,e_n_)&
                         -twofl_kin_en_n(w,ixI^L,ixO^L))
              if(phys_total_energy) then
                w(ixO^S,e_c_)=gamma_1*(w(ixO^S,e_c_)&
                          -twofl_kin_en_c(w,ixI^L,ixO^L)&
                          -twofl_mag_en(w,ixI^L,ixO^L))
               else
                 w(ixO^S,e_c_)=gamma_1*(w(ixO^S,e_c_)&
                          -twofl_kin_en_c(w,ixI^L,ixO^L))

               endif  
              if(phys_solve_eaux) w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)*gamma_1
            end if
          end if
          ! Convert momentum to velocity
          if(has_equi_rho_n0) then
            tmp1(ixO^S) = w(ixO^S,rho_n_) + block%equi_vars(ixO^S,equi_rho_n0_,0)
          else  
            tmp1(ixO^S) = w(ixO^S,rho_n_) 
          endif
    
          if(has_equi_rho_c0) then
            tmp2(ixO^S) = w(ixO^S,rho_c_) + block%equi_vars(ixO^S,equi_rho_c0_,0)
          else  
            tmp2(ixO^S) = w(ixO^S,rho_c_) 
          endif
          do idir = 1, ndir
             w(ixO^S, mom_n(idir)) = w(ixO^S, mom_n(idir))/tmp1(ixO^S)
             w(ixO^S, mom_c(idir)) = w(ixO^S, mom_c(idir))/tmp2(ixO^S)
          end do
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
  end subroutine twofl_handle_small_values

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine twofl_get_cmax(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision                :: vc(ixI^S)
    double precision                :: cmax2(ixI^S)
    double precision                :: vn(ixI^S)

    call twofl_get_csound_c_idim(w,x,ixI^L,ixO^L,idim,cmax)
    call twofl_get_v_c_idim(w,x,ixI^L,ixO^L,idim,vc)
    call twofl_get_v_n_idim(w,x,ixI^L,ixO^L,idim,vn)
    call twofl_get_csound_n(w,x,ixI^L,ixO^L,cmax2)
    cmax(ixO^S)=max(abs(vn(ixO^S))+cmax2(ixO^S),& 
            abs(vc(ixO^S))+cmax(ixO^S))
        
  end subroutine twofl_get_cmax

  subroutine twofl_get_a2max(w,x,ixI^L,ixO^L,a2max)
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
  end subroutine twofl_get_a2max

  ! COPIED from hd/moh_hd_phys
  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine twofl_get_tcutoff_n(ixI^L,ixO^L,w,x,tco_local,Tmax_local)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim),w(ixI^S,1:nw)
    double precision, intent(out) :: tco_local, Tmax_local

    double precision, parameter :: delta=0.25d0
    double precision :: tmp1(ixI^S),Te(ixI^S),lts(ixI^S)
    integer :: jxO^L,hxO^L
    logical :: lrlt(ixI^S)

    {^IFONED
    ! reuse lts as rhon
    call get_rhon_tot(w,x,ixI^L,ixI^L,lts)
    tmp1(ixI^S)=w(ixI^S,e_n_)-0.5d0*sum(w(ixI^S,mom_n(:))**2,dim=ndim+1)/lts(ixI^S)
    Te(ixI^S)=tmp1(ixI^S)/lts(ixI^S)*(twofl_gamma-1.d0)

    Tmax_local=maxval(Te(ixO^S))

    hxO^L=ixO^L-1;
    jxO^L=ixO^L+1;
    lts(ixO^S)=0.5d0*abs(Te(jxO^S)-Te(hxO^S))/Te(ixO^S)
    lrlt=.false.
    where(lts(ixO^S) > delta)
      lrlt(ixO^S)=.true.
    end where
    tco_local=zero
    if(any(lrlt(ixO^S))) then
      tco_local=maxval(Te(ixO^S), mask=lrlt(ixO^S))
    end if
    }
  end subroutine twofl_get_tcutoff_n

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine twofl_get_tcutoff_c(ixI^L,ixO^L,w,x,Tco_local,Tmax_local)
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
    double precision :: ltr(ixI^S),ltrc,ltrp,altr(ixI^S)
    integer :: idims,jxO^L,hxO^L,ixA^D,ixB^D
    integer :: jxP^L,hxP^L,ixP^L
    logical :: lrlt(ixI^S)

    ! reuse lts as rhoc
    call get_rhoc_tot(w,x,ixI^L,ixI^L,lts)
    if(phys_internal_e) then
      tmp1(ixI^S)=w(ixI^S,e_c_)
    else
      tmp1(ixI^S)=w(ixI^S,e_c_)-0.5d0*(sum(w(ixI^S,mom_c(:))**2,dim=ndim+1)/&
                       lts(ixI^S)+sum(w(ixI^S,mag(:))**2,dim=ndim+1))
    end if
    Te(ixI^S)=tmp1(ixI^S)/lts(ixI^S)*(twofl_gamma-1.d0)
    Tmax_local=maxval(Te(ixO^S))

    {^IFONED
    select case(twofl_trac_type)
    case(0)
      !> test case, fixed cutoff temperature
      w(ixI^S,Tcoff_c_)=2.5d5/unit_temperature
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
      ltrp=2.5d0
      ixP^L=ixO^L^LADD1;
      hxO^L=ixO^L-1;
      jxO^L=ixO^L+1;
      hxP^L=ixP^L-1;
      jxP^L=ixP^L+1;
      lts(ixP^S)=0.5d0*abs(Te(jxP^S)-Te(hxP^S))/Te(ixP^S)
      ltr(ixP^S)=max(one, (exp(lts(ixP^S))/ltrc)**ltrp)
      w(ixO^S,Tcoff_c_)=Te(ixO^S)*&
        (0.25*(ltr(jxO^S)+two*ltr(ixO^S)+ltr(hxO^S)))**0.4d0
    case default
      call mpistop("twofl_trac_type not allowed for 1D simulation")
    end select
    }
    {^NOONED
    select case(twofl_trac_type)
    case(0)
      !> test case, fixed cutoff temperature
      w(ixI^S,Tcoff_c_)=2.5d5/unit_temperature
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
      if(twofl_trac_type .gt. 1) then
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
      ltr(ixP^S)=max(one, (exp(lts(ixP^S))/ltrc)**ltrp)
  
      altr(ixI^S)=zero  
      do idims=1,ndim 
        hxO^L=ixO^L-kr(idims,^D);
        jxO^L=ixO^L+kr(idims,^D);
        altr(ixO^S)=altr(ixO^S) &
          +0.25*(ltr(hxO^S)+two*ltr(ixO^S)+ltr(jxO^S))*bunitvec(ixO^S,idims)**2
        w(ixO^S,Tcoff_c_)=Te(ixO^S)*altr(ixO^S)**(0.4*ltrp)
      end do
    case(3,5)
      !> do nothing here
    case default
      call mpistop("unknown twofl_trac_type")
    end select 
    }
  end subroutine twofl_get_tcutoff_c

  !> get H speed for H-correction to fix the carbuncle problem at grid-aligned shock front
  subroutine twofl_get_H_speed_one(wprim,x,ixI^L,ixO^L,idim,Hspeed) 
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
      call twofl_get_csound_prim(wprim,x,ixI^L,ixA^L,id,tmp)
      csound(ixA^S,id)=tmp(ixA^S)
    end do
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D+kr(idim,^D)-1;
    jxCmax^D=ixCmax^D+kr(idim,^D);
    jxCmin^D=ixCmin^D+kr(idim,^D);
    Hspeed(ixC^S,1)=0.5d0*abs(&
      0.5d0 * (wprim(jxC^S,mom_c(idim))+ wprim(jxC^S,mom_n(idim))) &
      +csound(jxC^S,idim)- &
      0.5d0 * (wprim(ixC^S,mom_c(idim)) + wprim(ixC^S,mom_n(idim)))&
      +csound(ixC^S,idim))

    do id=1,ndim
      if(id==idim) cycle
      ixAmax^D=ixCmax^D+kr(id,^D);
      ixAmin^D=ixCmin^D+kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(&
        0.5d0 * (wprim(ixA^S,mom_c(id)) + wprim(ixA^S,mom_n(id)))&
        +csound(ixA^S,id)-&
        0.5d0 * (wprim(ixC^S,mom_c(id)) + wprim(ixC^S,mom_n(id)))&
        +csound(ixC^S,id)))


      ixAmax^D=ixCmax^D-kr(id,^D);
      ixAmin^D=ixCmin^D-kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(&
        0.5d0 * (wprim(ixC^S,mom_c(id)) + wprim(ixC^S,mom_n(id)))&
        +csound(ixC^S,id)-&
        0.5d0 * (wprim(ixA^S,mom_c(id)) + wprim(ixA^S,mom_n(id)))&
        +csound(ixA^S,id)))

    end do

    do id=1,ndim
      if(id==idim) cycle
      ixAmax^D=jxCmax^D+kr(id,^D);
      ixAmin^D=jxCmin^D+kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(&
        0.5d0 * (wprim(ixA^S,mom_c(id)) + wprim(ixA^S,mom_n(id)))&
        +csound(ixA^S,id)-&
        0.5d0 * (wprim(jxC^S,mom_c(id)) + wprim(jxC^S,mom_n(id)))&
        +csound(jxC^S,id)))
      ixAmax^D=jxCmax^D-kr(id,^D);
      ixAmin^D=jxCmin^D-kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(&
        0.5d0 * (wprim(jxC^S,mom_c(id)) + wprim(jxC^S,mom_n(id)))&
        +csound(jxC^S,id)-&
        0.5d0 * (wprim(ixA^S,mom_c(id)) + wprim(ixA^S,mom_n(id)))&
        +csound(ixA^S,id)))
    end do

  end subroutine twofl_get_H_speed_one

  !> get H speed for H-correction to fix the carbuncle problem at grid-aligned shock front
  subroutine twofl_get_H_speed_species(wprim,x,ixI^L,ixO^L,idim,Hspeed) 
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wprim(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: Hspeed(ixI^S,1:number_species)

    double precision :: csound(ixI^S,ndim),tmp(ixI^S)
    integer :: jxC^L, ixC^L, ixA^L, id, ix^D

    Hspeed=0.d0
    ! charges
    ixA^L=ixO^L^LADD1;
    do id=1,ndim
      call twofl_get_csound_prim_c(wprim,x,ixI^L,ixA^L,id,tmp)
      csound(ixA^S,id)=tmp(ixA^S)
    end do
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D+kr(idim,^D)-1;
    jxCmax^D=ixCmax^D+kr(idim,^D);
    jxCmin^D=ixCmin^D+kr(idim,^D);
    Hspeed(ixC^S,1)=0.5d0*abs(wprim(jxC^S,mom_c(idim))+csound(jxC^S,idim)-wprim(ixC^S,mom_c(idim))+csound(ixC^S,idim))

    do id=1,ndim
      if(id==idim) cycle
      ixAmax^D=ixCmax^D+kr(id,^D);
      ixAmin^D=ixCmin^D+kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(wprim(ixA^S,mom_c(id))+csound(ixA^S,id)-wprim(ixC^S,mom_c(id))+csound(ixC^S,id)))
      ixAmax^D=ixCmax^D-kr(id,^D);
      ixAmin^D=ixCmin^D-kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(wprim(ixC^S,mom_c(id))+csound(ixC^S,id)-wprim(ixA^S,mom_c(id))+csound(ixA^S,id)))
    end do

    do id=1,ndim
      if(id==idim) cycle
      ixAmax^D=jxCmax^D+kr(id,^D);
      ixAmin^D=jxCmin^D+kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(wprim(ixA^S,mom_c(id))+csound(ixA^S,id)-wprim(jxC^S,mom_c(id))+csound(jxC^S,id)))
      ixAmax^D=jxCmax^D-kr(id,^D);
      ixAmin^D=jxCmin^D-kr(id,^D);
      Hspeed(ixC^S,1)=max(Hspeed(ixC^S,1),0.5d0*abs(wprim(jxC^S,mom_c(id))+csound(jxC^S,id)-wprim(ixA^S,mom_c(id))+csound(ixA^S,id)))
    end do

    ! neutrals
    ixA^L=ixO^L^LADD1;
    do id=1,ndim
      call twofl_get_csound_prim_n(wprim,x,ixI^L,ixA^L,id,tmp)
      csound(ixA^S,id)=tmp(ixA^S)
    end do
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D+kr(idim,^D)-1;
    jxCmax^D=ixCmax^D+kr(idim,^D);
    jxCmin^D=ixCmin^D+kr(idim,^D);
    Hspeed(ixC^S,2)=0.5d0*abs(wprim(jxC^S,mom_n(idim))+csound(jxC^S,idim)-wprim(ixC^S,mom_n(idim))+csound(ixC^S,idim))

    do id=1,ndim
      if(id==idim) cycle
      ixAmax^D=ixCmax^D+kr(id,^D);
      ixAmin^D=ixCmin^D+kr(id,^D);
      Hspeed(ixC^S,2)=max(Hspeed(ixC^S,2),0.5d0*abs(wprim(ixA^S,mom_n(id))+csound(ixA^S,id)-wprim(ixC^S,mom_n(id))+csound(ixC^S,id)))
      ixAmax^D=ixCmax^D-kr(id,^D);
      ixAmin^D=ixCmin^D-kr(id,^D);
      Hspeed(ixC^S,2)=max(Hspeed(ixC^S,2),0.5d0*abs(wprim(ixC^S,mom_n(id))+csound(ixC^S,id)-wprim(ixA^S,mom_n(id))+csound(ixA^S,id)))
    end do

    do id=1,ndim
      if(id==idim) cycle
      ixAmax^D=jxCmax^D+kr(id,^D);
      ixAmin^D=jxCmin^D+kr(id,^D);
      Hspeed(ixC^S,2)=max(Hspeed(ixC^S,2),0.5d0*abs(wprim(ixA^S,mom_n(id))+csound(ixA^S,id)-wprim(jxC^S,mom_n(id))+csound(jxC^S,id)))
      ixAmax^D=jxCmax^D-kr(id,^D);
      ixAmin^D=jxCmin^D-kr(id,^D);
      Hspeed(ixC^S,2)=max(Hspeed(ixC^S,2),0.5d0*abs(wprim(jxC^S,mom_n(id))+csound(jxC^S,id)-wprim(ixA^S,mom_n(id))+csound(ixA^S,id)))
    end do

  end subroutine twofl_get_H_speed_species

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine twofl_get_cbounds_one(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,Hspeed,cmax,cmin)
    use mod_global_parameters
    use mod_constrained_transport

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S,number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,number_species)
    double precision, intent(in)    :: Hspeed(ixI^S,1:number_species)

    double precision :: wmean(ixI^S,nw)
    double precision :: rhon(ixI^S)
    double precision :: rhoc(ixI^S)
    double precision, dimension(ixI^S) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix^D

    select case (boundspeed) 
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      call get_rhoc_tot(wLP,x,ixI^L,ixO^L,rhoc)
      call get_rhon_tot(wLP,x,ixI^L,ixO^L,rhon)
      tmp1(ixO^S)=sqrt(abs(rhoc(ixO^S)  +rhon(ixO^S)))

      call get_rhoc_tot(wRP,x,ixI^L,ixO^L,rhoc)
      call get_rhon_tot(wRP,x,ixI^L,ixO^L,rhon)
      tmp2(ixO^S)=sqrt(abs(rhoc(ixO^S) +rhon(ixO^S)))

      tmp3(ixO^S)=1.d0/(tmp1(ixO^S)+tmp2(ixO^S))
      umean(ixO^S)=(0.5*(wLp(ixO^S,mom_n(idim))+wLp(ixO^S,mom_c(idim)))*tmp1(ixO^S) + &
                    0.5*(wRp(ixO^S,mom_n(idim))+wRp(ixO^S,mom_c(idim)))*tmp2(ixO^S))*tmp3(ixO^S)
      call twofl_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call twofl_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)

      dmean(ixO^S)=(tmp1(ixO^S)*csoundL(ixO^S)**2+tmp2(ixO^S)*csoundR(ixO^S)**2)*tmp3(ixO^S)+&
       0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2*(&
       0.5*(wRp(ixO^S,mom_n(idim))+wRp(ixO^S,mom_c(idim)))- & 
       0.5*(wLp(ixO^S,mom_n(idim))+wLp(ixO^S,mom_c(idim))))**2
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
    ! typeboundspeed=='cmaxmean'
      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
      call get_rhon_tot(wmean,x,ixI^L,ixO^L,rhon)
      tmp2(ixO^S)=wmean(ixO^S,mom_n(idim))/rhon(ixO^S)
      call get_rhoc_tot(wmean,x,ixI^L,ixO^L,rhoc)
      tmp1(ixO^S)=wmean(ixO^S,mom_c(idim))/rhoc(ixO^S)
      call twofl_get_csound(wmean,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
        cmax(ixO^S,1)=max(max(abs(tmp2(ixO^S)), abs(tmp1(ixO^S)) ) +csoundR(ixO^S),zero)
        cmin(ixO^S,1)=min(min(abs(tmp2(ixO^S)),  abs(tmp1(ixO^S)) ) -csoundR(ixO^S),zero)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)= max(abs(tmp2(ixO^S)),abs(tmp1(ixO^S)))+csoundR(ixO^S)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call twofl_get_csound(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call twofl_get_csound(wRp,x,ixI^L,ixO^L,idim,csoundR)
      csoundL(ixO^S)=max(csoundL(ixO^S),csoundR(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,1)=min(0.5*(wLp(ixO^S,mom_c(idim))+ wRp(ixO^S,mom_n(idim))),&
            0.5*(wRp(ixO^S,mom_c(idim))+ wRp(ixO^S,mom_n(idim))))-csoundL(ixO^S)
        cmax(ixO^S,1)=max(0.5*(wLp(ixO^S,mom_c(idim))+ wRp(ixO^S,mom_n(idim))),&
            0.5*(wRp(ixO^S,mom_c(idim))+ wRp(ixO^S,mom_n(idim))))+csoundL(ixO^S)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=max(0.5*(wLp(ixO^S,mom_c(idim))+ wRp(ixO^S,mom_n(idim))),&
          0.5*(wRp(ixO^S,mom_c(idim))+ wRp(ixO^S,mom_n(idim))))+csoundL(ixO^S)
      end if
    end select

  end subroutine twofl_get_cbounds_one

  !> Calculate fast magnetosonic wave speed
  subroutine twofl_get_csound_prim_c(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
    double precision :: inv_rho(ixO^S)
    double precision :: rhoc(ixI^S)

    integer :: ix1,ix2


    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
    inv_rho(ixO^S)=1.d0/rhoc(ixO^S)

    if(phys_energy) then
      call twofl_get_pthermal_c_primitive(w,x,ixI^L,ixO^L,csound)
      csound(ixO^S)=twofl_gamma*csound(ixO^S)/rhoc(ixO^S)
    else
      call twofl_get_csound2_adiab_c(w,x,ixI^L,ixO^L,csound)
    endif

    ! store |B|^2 in v
    b2(ixO^S)        = twofl_mag_en_all(w,ixI^L,ixO^L)
    cfast2(ixO^S)   = b2(ixO^S) * inv_rho(ixO^S)+csound(ixO^S)
    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
         * twofl_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
         * inv_rho(ixO^S)

    where(AvMinCs2(ixO^S)<zero)
       AvMinCs2(ixO^S)=zero
    end where

    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

    if (.not. twofl_Hall) then
       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
            twofl_etah * sqrt(b2(ixO^S))*inv_rho(ixO^S)*kmax)
    end if

  end subroutine twofl_get_csound_prim_c

  !> Calculate fast magnetosonic wave speed
  subroutine twofl_get_csound_prim_n(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: rhon(ixI^S)

    if(phys_energy) then
      call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
      call twofl_get_pthermal_n_primitive(w,x,ixI^L,ixO^L,csound)
      csound(ixO^S)=twofl_gamma*csound(ixO^S)/rhon(ixO^S)
    else
      call twofl_get_csound2_adiab_n(w,x,ixI^L,ixO^L,csound)
    endif
    csound(ixO^S) = sqrt(csound(ixO^S))

  end subroutine twofl_get_csound_prim_n

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine twofl_get_cbounds_species(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,Hspeed,cmax,cmin)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_variables

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
    double precision, intent(in)    :: Hspeed(ixI^S,1:number_species)

    double precision :: wmean(ixI^S,nw)
    double precision :: rho(ixI^S)
    double precision, dimension(ixI^S) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix^D

    select case (boundspeed) 
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      ! charges
      call get_rhoc_tot(wLP,x,ixI^L,ixO^L,rho)
      tmp1(ixO^S)=sqrt(abs(rho(ixO^S)))

      call get_rhoc_tot(wRP,x,ixI^L,ixO^L,rho)
      tmp2(ixO^S)=sqrt(abs(rho(ixO^S)))

      tmp3(ixO^S)=1.d0/(tmp1(ixO^S)+tmp2(ixO^S))
      umean(ixO^S)=(wLp(ixO^S,mom_c(idim))*tmp1(ixO^S)+wRp(ixO^S,mom_c(idim))*tmp2(ixO^S))*tmp3(ixO^S)
      call twofl_get_csound_prim_c(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call twofl_get_csound_prim_c(wRp,x,ixI^L,ixO^L,idim,csoundR)


      dmean(ixO^S)=(tmp1(ixO^S)*csoundL(ixO^S)**2+tmp2(ixO^S)*csoundR(ixO^S)**2)*tmp3(ixO^S)+&
       0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2*&
       (wRp(ixO^S,mom_c(idim)) - wLp(ixO^S,mom_c(idim)))**2
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
    
      ! neutrals

      call get_rhon_tot(wLP,x,ixI^L,ixO^L,rho)
      tmp1(ixO^S)=sqrt(abs(rho(ixO^S)))

      call get_rhon_tot(wRP,x,ixI^L,ixO^L,rho)
      tmp2(ixO^S)=sqrt(abs(rho(ixO^S)))

      tmp3(ixO^S)=1.d0/(tmp1(ixO^S)+tmp2(ixO^S))
      umean(ixO^S)=(wLp(ixO^S,mom_n(idim))*tmp1(ixO^S)+wRp(ixO^S,mom_n(idim))*tmp2(ixO^S))*tmp3(ixO^S)
      call twofl_get_csound_prim_n(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call twofl_get_csound_prim_n(wRp,x,ixI^L,ixO^L,idim,csoundR)


      dmean(ixO^S)=(tmp1(ixO^S)*csoundL(ixO^S)**2+tmp2(ixO^S)*csoundR(ixO^S)**2)*tmp3(ixO^S)+&
       0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2*&
       (wRp(ixO^S,mom_n(idim)) - wLp(ixO^S,mom_n(idim)))**2
      dmean(ixO^S)=sqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,2)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S,2)=umean(ixO^S)+dmean(ixO^S)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,2)=sign(one,cmin(ix^D,2))*max(abs(cmin(ix^D,2)),Hspeed(ix^D,2))
            cmax(ix^D,2)=sign(one,cmax(ix^D,2))*max(abs(cmax(ix^D,2)),Hspeed(ix^D,2))
          {end do\}
        end if
      else
        cmax(ixO^S,2)=abs(umean(ixO^S))+dmean(ixO^S)
      end if

    case (2)
    ! typeboundspeed=='cmaxmean'
      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
     ! charges 

      call get_rhoc_tot(wmean,x,ixI^L,ixO^L,rho)
      tmp1(ixO^S)=wmean(ixO^S,mom_c(idim))/rho(ixO^S)
      call twofl_get_csound_c_idim(wmean,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
        cmax(ixO^S,1)=max(abs(tmp1(ixO^S))+csoundR(ixO^S),zero)
        cmin(ixO^S,1)=min(abs(tmp1(ixO^S))-csoundR(ixO^S),zero)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=abs(tmp1(ixO^S))+csoundR(ixO^S)
      end if
      !neutrals
      
      call get_rhon_tot(wmean,x,ixI^L,ixO^L,rho)
      tmp1(ixO^S)=wmean(ixO^S,mom_n(idim))/rho(ixO^S)
      call twofl_get_csound_n(wmean,x,ixI^L,ixO^L,csoundR)
      if(present(cmin)) then
        cmax(ixO^S,2)=max(abs(tmp1(ixO^S))+csoundR(ixO^S),zero)
        cmin(ixO^S,2)=min(abs(tmp1(ixO^S))-csoundR(ixO^S),zero)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,2)=sign(one,cmin(ix^D,2))*max(abs(cmin(ix^D,2)),Hspeed(ix^D,2))
            cmax(ix^D,2)=sign(one,cmax(ix^D,2))*max(abs(cmax(ix^D,2)),Hspeed(ix^D,2))
          {end do\}
        end if
      else
        cmax(ixO^S,2)= abs(tmp1(ixO^S))+csoundR(ixO^S)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call twofl_get_csound_c_idim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call twofl_get_csound_c_idim(wRp,x,ixI^L,ixO^L,idim,csoundR)
      csoundL(ixO^S)=max(csoundL(ixO^S),csoundR(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,1)=min(wLp(ixO^S,mom_c(idim)),wRp(ixO^S,mom_c(idim)))-csoundL(ixO^S)
        cmax(ixO^S,1)=max(wLp(ixO^S,mom_c(idim)),wRp(ixO^S,mom_c(idim)))+csoundL(ixO^S)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=max(wLp(ixO^S,mom_c(idim)),wRp(ixO^S,mom_c(idim)))+csoundL(ixO^S)
      end if
      call twofl_get_csound_n(wLp,x,ixI^L,ixO^L,csoundL)
      call twofl_get_csound_n(wRp,x,ixI^L,ixO^L,csoundR)
      csoundL(ixO^S)=max(csoundL(ixO^S),csoundR(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,2)=min(wLp(ixO^S,mom_n(idim)),wRp(ixO^S,mom_n(idim)))-csoundL(ixO^S)
        cmax(ixO^S,2)=max(wLp(ixO^S,mom_n(idim)),wRp(ixO^S,mom_n(idim)))+csoundL(ixO^S)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,2)=sign(one,cmin(ix^D,2))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,2))
            cmax(ix^D,2)=sign(one,cmax(ix^D,2))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,2))
          {end do\}
        end if
      else
        cmax(ixO^S,2)=max(wLp(ixO^S,mom_n(idim)),wRp(ixO^S,mom_n(idim)))+csoundL(ixO^S)
      end if

    end select

  end subroutine twofl_get_cbounds_species

  !> prepare velocities for ct methods
  subroutine twofl_get_ct_velocity(vcts,wLp,wRp,ixI^L,ixO^L,idim,cmax,cmin)
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
      vcts%vnorm(ixO^S,idim)=0.5d0*(wLp(ixO^S,mom_c(idim))+wRp(ixO^S,mom_c(idim)))
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
      vcts%vbarLC(ixO^S,idim,1)=wLp(ixO^S,mom_c(idimN))
      vcts%vbarRC(ixO^S,idim,1)=wRp(ixO^S,mom_c(idimN))
      vcts%vbarC(ixO^S,idim,1)=(vcts%cbarmax(ixO^S,idim)*vcts%vbarLC(ixO^S,idim,1) &
           +vcts%cbarmin(ixO^S,idim)*vcts%vbarRC(ixO^S,idim,1))&
          /(vcts%cbarmax(ixO^S,idim)+vcts%cbarmin(ixO^S,idim))

      vcts%vbarLC(ixO^S,idim,2)=wLp(ixO^S,mom_c(idimE))
      vcts%vbarRC(ixO^S,idim,2)=wRp(ixO^S,mom_c(idimE))
      vcts%vbarC(ixO^S,idim,2)=(vcts%cbarmax(ixO^S,idim)*vcts%vbarLC(ixO^S,idim,2) &
           +vcts%cbarmin(ixO^S,idim)*vcts%vbarRC(ixO^S,idim,1))&
          /(vcts%cbarmax(ixO^S,idim)+vcts%cbarmin(ixO^S,idim))
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

  end subroutine twofl_get_ct_velocity

  subroutine twofl_get_csound_c_idim(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
    double precision :: inv_rho(ixO^S)
    double precision :: tmp(ixI^S)
#if (!defined(ONE_FLUID) || ONE_FLUID==0) && (defined(A_TOT) && A_TOT == 1)
    double precision :: rhon(ixI^S)
#endif
    call get_rhoc_tot(w,x,ixI^L,ixO^L,tmp)
#if (!defined(ONE_FLUID) || ONE_FLUID==0) && (defined(A_TOT) && A_TOT == 1)
    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    inv_rho(ixO^S) = 1d0/(rhon(ixO^S)+tmp(ixO^S)) 
#else
    inv_rho(ixO^S)=1.d0/tmp(ixO^S)
#endif

    call twofl_get_csound2_c_from_conserved(w,x,ixI^L,ixO^L,csound)

    ! store |B|^2 in v
    b2(ixO^S) = twofl_mag_en_all(w,ixI^L,ixO^L)

    cfast2(ixO^S)   = b2(ixO^S) * inv_rho(ixO^S)+csound(ixO^S)
    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
         * twofl_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
         * inv_rho(ixO^S)

    where(AvMinCs2(ixO^S)<zero)
       AvMinCs2(ixO^S)=zero
    end where

    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

    if (.not. twofl_Hall) then
       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
            twofl_etah * sqrt(b2(ixO^S))*inv_rho(ixO^S)*kmax)
    end if

  end subroutine twofl_get_csound_c_idim

  !> Calculate fast magnetosonic wave speed when cbounds_species=false
  subroutine twofl_get_csound_prim(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
    double precision :: inv_rho(ixO^S)
    double precision :: rhoc(ixI^S)
#if (defined(A_TOT) && A_TOT == 1)
    double precision :: rhon(ixI^S)
#endif
    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
#if  (defined(A_TOT) && A_TOT == 1)
    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    inv_rho(ixO^S) = 1d0/(rhon(ixO^S)+rhoc(ixO^S)) 
#else
    inv_rho(ixO^S)=1.d0/rhoc(ixO^S)
#endif

    call twofl_get_csound2_primitive(w,x,ixI^L,ixO^L,csound)

    ! store |B|^2 in v
    b2(ixO^S)        = twofl_mag_en_all(w,ixI^L,ixO^L)
    cfast2(ixO^S)   = b2(ixO^S) * inv_rho(ixO^S)+csound(ixO^S)
    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
         * twofl_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
         * inv_rho(ixO^S)

    where(AvMinCs2(ixO^S)<zero)
       AvMinCs2(ixO^S)=zero
    end where

    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

    if (.not. twofl_Hall) then
       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
            twofl_etah * sqrt(b2(ixO^S))*inv_rho(ixO^S)*kmax)
    end if

    contains
    !TODO copy it inside
    subroutine twofl_get_csound2_primitive(w,x,ixI^L,ixO^L,csound2)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: w(ixI^S,nw)
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      double precision, intent(out)   :: csound2(ixI^S)
      double precision  :: pth_c(ixI^S)
      double precision  :: pth_n(ixI^S)
  
      if(phys_energy) then
        call twofl_get_pthermal_c_primitive(w,x,ixI^L,ixO^L,pth_c)
        call twofl_get_pthermal_n_primitive(w,x,ixI^L,ixO^L,pth_n)
        call twofl_get_csound2_from_pthermal(w,x,ixI^L,ixO^L,pth_c,pth_n,csound2)
      else
        call twofl_get_csound2_adiab(w,x,ixI^L,ixO^L,csound2)
      endif
    end subroutine twofl_get_csound2_primitive

  end subroutine twofl_get_csound_prim

  subroutine twofl_get_csound2(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision  :: pth_c(ixI^S)
    double precision  :: pth_n(ixI^S)

    if(phys_energy) then
      call twofl_get_pthermal_c(w,x,ixI^L,ixO^L,pth_c)
      call twofl_get_pthermal_n(w,x,ixI^L,ixO^L,pth_n)
      call twofl_get_csound2_from_pthermal(w,x,ixI^L,ixO^L,pth_c,pth_n,csound2)
    else
      call twofl_get_csound2_adiab(w,x,ixI^L,ixO^L,csound2)
    endif
  end subroutine twofl_get_csound2

  subroutine twofl_get_csound2_adiab(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision  :: rhoc(ixI^S)
    double precision  :: rhon(ixI^S)

    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    csound2(ixO^S)=twofl_gamma*twofl_adiab*&
                  max((rhoc(ixO^S)**twofl_gamma + rhon(ixO^S)**twofl_gamma)/(rhoc(ixO^S)+ rhon(ixO^S)),&
                  rhon(ixO^S)**gamma_1,rhoc(ixO^S)**gamma_1) 
  end subroutine twofl_get_csound2_adiab

  subroutine twofl_get_csound(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
    double precision :: inv_rho(ixO^S)
    double precision :: rhoc(ixI^S)
#if (defined(A_TOT) && A_TOT == 1)
    double precision :: rhon(ixI^S)
#endif
    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
#if (defined(A_TOT) && A_TOT == 1)
    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    inv_rho(ixO^S) = 1d0/(rhon(ixO^S)+rhoc(ixO^S)) 
#else
    inv_rho(ixO^S)=1.d0/rhoc(ixO^S)
#endif

    call twofl_get_csound2(w,x,ixI^L,ixO^L,csound)

    ! store |B|^2 in v
    b2(ixO^S) = twofl_mag_en_all(w,ixI^L,ixO^L)

    cfast2(ixO^S)   = b2(ixO^S) * inv_rho(ixO^S)+csound(ixO^S)
    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
         * twofl_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
         * inv_rho(ixO^S)

    where(AvMinCs2(ixO^S)<zero)
       AvMinCs2(ixO^S)=zero
    end where

    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

    if (.not. twofl_Hall) then
       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
            twofl_etah * sqrt(b2(ixO^S))*inv_rho(ixO^S)*kmax)
    end if

  end subroutine twofl_get_csound

  subroutine twofl_get_csound2_from_pthermal(w,x,ixI^L,ixO^L,pth_c,pth_n,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: pth_c(ixI^S)
    double precision, intent(in)    :: pth_n(ixI^S)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision  :: csound1(ixI^S),rhon(ixI^S),rhoc(ixI^S)

    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
#if !defined(C_TOT) || C_TOT == 0
    csound2(ixO^S)=twofl_gamma*max((pth_c(ixO^S) + pth_n(ixO^S))/(rhoc(ixO^S) + rhon(ixO^S)),&
                      pth_n(ixO^S)/rhon(ixO^S), pth_c(ixO^S)/rhoc(ixO^S))
#else
    csound2(ixO^S)=twofl_gamma*(csound2(ixO^S) + csound1(ixO^S))/(rhoc(ixO^S) + rhon(ixO^S))

#endif
  end subroutine twofl_get_csound2_from_pthermal

! end cbounds_species=false

  subroutine twofl_get_csound_n(w,x,ixI^L,ixO^L,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: pe_n1(ixI^S)
    call twofl_get_csound2_n_from_conserved(w,x,ixI^L,ixO^L,csound)
    csound(ixO^S) = sqrt(csound(ixO^S))
  end subroutine twofl_get_csound_n

  !> separate routines so that it is faster
  !> Calculate temperature=p/rho when in e_ the internal energy is stored
  subroutine twofl_get_temperature_from_eint_n(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)

    res(ixO^S) = 1d0/Rn * gamma_1 * w(ixO^S, e_n_) /w(ixO^S,rho_n_)

  end subroutine twofl_get_temperature_from_eint_n

  subroutine twofl_get_temperature_from_eint_n_with_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)

      res(ixO^S) = 1d0/Rn * (gamma_1 * w(ixO^S, e_n_) + block%equi_vars(ixO^S,equi_pe_n0_,b0i)) /&
                (w(ixO^S,rho_n_) +block%equi_vars(ixO^S,equi_rho_n0_,b0i))
  end subroutine twofl_get_temperature_from_eint_n_with_equi

!  subroutine twofl_get_temperature_n_pert_from_tot(Te, ixI^L, ixO^L, res)
!    use mod_global_parameters
!    integer, intent(in)          :: ixI^L, ixO^L
!    double precision, intent(in) :: Te(ixI^S)
!    double precision, intent(out):: res(ixI^S)
!      res(ixO^S) = Te(ixO^S) -1d0/Rn * &
!                block%equi_vars(ixO^S,equi_pe_n0_,0)/block%equi_vars(ixO^S,equi_rho_n0_,0)
!  end subroutine twofl_get_temperature_n_pert_from_tot

  subroutine twofl_get_temperature_n_equi(w,x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
      res(ixO^S) = 1d0/Rn * &
                block%equi_vars(ixO^S,equi_pe_n0_,b0i)/block%equi_vars(ixO^S,equi_rho_n0_,b0i)
  end subroutine twofl_get_temperature_n_equi

  subroutine twofl_get_rho_n_equi(w, x,ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
      res(ixO^S) = block%equi_vars(ixO^S,equi_rho_n0_,b0i)
  end subroutine twofl_get_rho_n_equi

  subroutine twofl_get_pe_n_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
      res(ixO^S) = block%equi_vars(ixO^S,equi_pe_n0_,b0i)
  end subroutine twofl_get_pe_n_equi

  !> Calculate temperature=p/rho when in e_ the total energy is stored
  !> this does not check the values of twofl_energy and twofl_internal_e, 
  !>  twofl_energy = .true. and twofl_internal_e = .false.
  !> also check small_values is avoided
  subroutine twofl_get_temperature_from_etot_n(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S)=1d0/Rn * (gamma_1*(w(ixO^S,e_n_)&
           - twofl_kin_en_n(w,ixI^L,ixO^L)))/w(ixO^S,rho_n_)
  end subroutine twofl_get_temperature_from_etot_n

  subroutine twofl_get_temperature_from_etot_n_with_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S)=1d0/Rn * (gamma_1*(w(ixO^S,e_n_)&
           - twofl_kin_en_n(w,ixI^L,ixO^L)) +  block%equi_vars(ixO^S,equi_pe_n0_,b0i))&
            /(w(ixO^S,rho_n_) +block%equi_vars(ixO^S,equi_rho_n0_,b0i))
            
  end subroutine twofl_get_temperature_from_etot_n_with_equi

  !> separate routines so that it is faster
  !> Calculate temperature=p/rho when in e_ the internal energy is stored
  subroutine twofl_get_temperature_from_eint_c(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)

    res(ixO^S) = 1d0/Rc * gamma_1 * w(ixO^S, e_c_) /w(ixO^S,rho_c_)

  end subroutine twofl_get_temperature_from_eint_c

  subroutine twofl_get_temperature_from_eint_c_with_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
      res(ixO^S) = 1d0/Rc * (gamma_1 * w(ixO^S, e_c_) + block%equi_vars(ixO^S,equi_pe_c0_,b0i)) /&
                (w(ixO^S,rho_c_) +block%equi_vars(ixO^S,equi_rho_c0_,b0i))
  end subroutine twofl_get_temperature_from_eint_c_with_equi

!  subroutine twofl_get_temperature_c_pert_from_tot(Te, ixI^L, ixO^L, res)
!    use mod_global_parameters
!    integer, intent(in)          :: ixI^L, ixO^L
!    double precision, intent(in) :: Te(ixI^S)
!    double precision, intent(out):: res(ixI^S)
!      res(ixO^S) = Te(ixO^S) -1d0/Rc * &
!                block%equi_vars(ixO^S,equi_pe_c0_,0)/block%equi_vars(ixO^S,equi_rho_c0_,0)
!  end subroutine twofl_get_temperature_c_pert_from_tot

  subroutine twofl_get_temperature_c_equi(w,x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
      res(ixO^S) = 1d0/Rc * &
                block%equi_vars(ixO^S,equi_pe_c0_,b0i)/block%equi_vars(ixO^S,equi_rho_c0_,b0i)
  end subroutine twofl_get_temperature_c_equi

  subroutine twofl_get_rho_c_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
      res(ixO^S) = block%equi_vars(ixO^S,equi_rho_c0_,b0i)
  end subroutine twofl_get_rho_c_equi

  subroutine twofl_get_pe_c_equi(w,x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
      res(ixO^S) = block%equi_vars(ixO^S,equi_pe_c0_,b0i)
  end subroutine twofl_get_pe_c_equi

  !> Calculate temperature=p/rho when in e_ the total energy is stored
  !> this does not check the values of twofl_energy and twofl_internal_e, 
  !>  twofl_energy = .true. and twofl_internal_e = .false.
  !> also check small_values is avoided
  subroutine twofl_get_temperature_from_etot_c(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S)=1d0/Rc * (gamma_1*(w(ixO^S,e_c_)&
           - twofl_kin_en_c(w,ixI^L,ixO^L)&
           - twofl_mag_en(w,ixI^L,ixO^L)))/w(ixO^S,rho_c_)
  end subroutine twofl_get_temperature_from_etot_c
  subroutine twofl_get_temperature_from_eki_c(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S)=1d0/Rc * (gamma_1*(w(ixO^S,e_c_)&
           - twofl_kin_en_c(w,ixI^L,ixO^L)))/w(ixO^S,rho_c_)
  end subroutine twofl_get_temperature_from_eki_c

  subroutine twofl_get_temperature_from_etot_c_with_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S)=1d0/Rc * (gamma_1*(w(ixO^S,e_c_)&
           - twofl_kin_en_c(w,ixI^L,ixO^L)&
           - twofl_mag_en(w,ixI^L,ixO^L)) +  block%equi_vars(ixO^S,equi_pe_c0_,b0i))&
            /(w(ixO^S,rho_c_) +block%equi_vars(ixO^S,equi_rho_c0_,b0i))
            
  end subroutine twofl_get_temperature_from_etot_c_with_equi

  subroutine twofl_get_temperature_from_eki_c_with_equi(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)
    res(ixO^S)=1d0/Rc * (gamma_1*(w(ixO^S,e_c_)&
           - twofl_kin_en_c(w,ixI^L,ixO^L)) +  block%equi_vars(ixO^S,equi_pe_c0_,b0i))&
            /(w(ixO^S,rho_c_) +block%equi_vars(ixO^S,equi_rho_c0_,b0i))
            
  end subroutine twofl_get_temperature_from_eki_c_with_equi

  subroutine twofl_get_csound2_adiab_n(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision  :: rhon(ixI^S)

    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    csound2(ixO^S)=twofl_gamma*twofl_adiab*rhon(ixO^S)**gamma_1

  end subroutine twofl_get_csound2_adiab_n

  subroutine twofl_get_csound2_n_from_conserved(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision                :: rhon(ixI^S)

    if(phys_energy) then
      call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
      call twofl_get_pthermal_n(w,x,ixI^L,ixO^L,csound2)
      csound2(ixO^S)=twofl_gamma*csound2(ixO^S)/rhon(ixO^S)
    else
      call twofl_get_csound2_adiab_n(w,x,ixI^L,ixO^L,csound2)
    endif
  end subroutine twofl_get_csound2_n_from_conserved

  !! TO DELETE
  subroutine twofl_get_csound2_n_from_primitive(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision                :: rhon(ixI^S)

    if(phys_energy) then
      call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
      call twofl_get_pthermal_n_primitive(w,x,ixI^L,ixO^L,csound2)
      csound2(ixO^S)=twofl_gamma*csound2(ixO^S)/rhon(ixO^S)
    else
      call twofl_get_csound2_adiab_n(w,x,ixI^L,ixO^L,csound2)
    endif
  end subroutine twofl_get_csound2_n_from_primitive

  subroutine twofl_get_csound2_adiab_c(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision  :: rhoc(ixI^S)

    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
    csound2(ixO^S)=twofl_gamma*twofl_adiab* rhoc(ixO^S)**gamma_1

  end subroutine twofl_get_csound2_adiab_c

  subroutine twofl_get_csound2_c_from_conserved(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision                :: rhoc(ixI^S)

    if(phys_energy) then
      call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
      call twofl_get_pthermal_c(w,x,ixI^L,ixO^L,csound2)
      csound2(ixO^S)=twofl_gamma*csound2(ixO^S)/rhoc(ixO^S)
    else
      call twofl_get_csound2_adiab_c(w,x,ixI^L,ixO^L,csound2)
    endif
  end subroutine twofl_get_csound2_c_from_conserved

  !> Calculate fluxes within ixO^L.
  subroutine twofl_get_flux(wC,w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)

    double precision             :: pgas(ixO^S), ptotal(ixO^S),tmp(ixI^S)
    double precision, allocatable:: vHall(:^D&,:)
    integer                      :: idirmin, iw, idir, jdir, kdir

    ! value at the interfaces, idim =  block%iw0 --> b0i 
    ! reuse tmp, used afterwards
    ! value at the interface so we can't put momentum
    call get_rhoc_tot(w,x,ixI^L,ixO^L,tmp)
    ! Get flux of density
    f(ixO^S,rho_c_)=w(ixO^S,mom_c(idim))*tmp(ixO^S)
    ! pgas is time dependent only
    if(phys_energy) then
      pgas(ixO^S)=w(ixO^S,e_c_)
    else
      pgas(ixO^S)=twofl_adiab*tmp(ixO^S)**twofl_gamma
      if(has_equi_pe_c0) then
        pgas(ixO^S)=pgas(ixO^S)-block%equi_vars(ixO^S,equi_pe_c0_,b0i)
      endif
    end if

    if (twofl_Hall) then
      allocate(vHall(ixI^S,1:ndir))
      call twofl_getv_Hall(w,x,ixI^L,ixO^L,vHall)
    end if

    if(B0field) tmp(ixO^S)=sum(block%B0(ixO^S,:,idim)*w(ixO^S,mag(:)),dim=ndim+1)

    ptotal(ixO^S) = pgas(ixO^S) + 0.5d0*sum(w(ixO^S, mag(:))**2, dim=ndim+1)

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixO^S,mom_c(idir))=ptotal(ixO^S)-w(ixO^S,mag(idim))*w(ixO^S,mag(idir))
        if(B0field) f(ixO^S,mom_c(idir))=f(ixO^S,mom_c(idir))+tmp(ixO^S)
      else
        f(ixO^S,mom_c(idir))= -w(ixO^S,mag(idir))*w(ixO^S,mag(idim))
      end if
      if (B0field) then
        f(ixO^S,mom_c(idir))=f(ixO^S,mom_c(idir))&
             -w(ixO^S,mag(idir))*block%B0(ixO^S,idim,idim)&
             -w(ixO^S,mag(idim))*block%B0(ixO^S,idir,idim)
      end if
      f(ixO^S,mom_c(idir))=f(ixO^S,mom_c(idir))+w(ixO^S,mom_c(idim))*wC(ixO^S,mom_c(idir))
    end do

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if(phys_energy) then
      if (phys_internal_e) then
         f(ixO^S,e_c_)=w(ixO^S,mom_c(idim))*w(ixO^S,e_c_)
         if (twofl_Hall) then
            call mpistop("solve internal energy not implemented for Hall MHD")
         endif
      else if(twofl_eq_energy == EQ_ENERGY_KI) then

        f(ixO^S,e_c_)=w(ixO^S,mom_c(idim))*(wC(ixO^S,e_c_)+pgas(ixO^S))
      else
        f(ixO^S,e_c_)=w(ixO^S,mom_c(idim))*(wC(ixO^S,e_c_)+ptotal(ixO^S))&
           -w(ixO^S,mag(idim))*sum(w(ixO^S,mag(:))*w(ixO^S,mom_c(:)),dim=ndim+1)
        !if(phys_solve_eaux) f(ixO^S,eaux_)=w(ixO^S,mom(idim))*wC(ixO^S,eaux_)

        if (B0field) then
           f(ixO^S,e_c_) = f(ixO^S,e_c_) &
              + w(ixO^S,mom_c(idim)) * tmp(ixO^S) &
              - sum(w(ixO^S,mom_c(:))*w(ixO^S,mag(:)),dim=ndim+1) * block%B0(ixO^S,idim,idim)
        end if

        if (twofl_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
           if (twofl_etah>zero) then
              f(ixO^S,e_c_) = f(ixO^S,e_c_) + vHall(ixO^S,idim) * &
                 sum(w(ixO^S, mag(:))**2,dim=ndim+1) &
                 - w(ixO^S,mag(idim)) * sum(vHall(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1)
              if (B0field) then
                 f(ixO^S,e_c_) = f(ixO^S,e_c_) &
                    + vHall(ixO^S,idim) * tmp(ixO^S) &
                    - sum(vHall(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1) * block%B0(ixO^S,idim,idim)
              end if
           end if
        end if
      end if !total_energy
      ! add flux of equilibrium internal energy corresponding to pe_c0
      if(has_equi_pe_c0) then
#if !defined(E_RM_W0) || E_RM_W0 == 1
        f(ixO^S,e_c_)=  f(ixO^S,e_c_) &
          + w(ixO^S,mom_c(idim)) * block%equi_vars(ixO^S,equi_pe_c0_,idim) * inv_gamma_1
#else
        if(phys_internal_e) then
          f(ixO^S,e_c_)=  f(ixO^S,e_c_) &
          + w(ixO^S,mom_c(idim)) * block%equi_vars(ixO^S,equi_pe_c0_,idim) * inv_gamma_1
        else
          f(ixO^S,e_c_)=  f(ixO^S,e_c_) &
          + w(ixO^S,mom_c(idim)) * block%equi_vars(ixO^S,equi_pe_c0_,idim) * twofl_gamma * inv_gamma_1
        endif
#endif
      end if
    end if !phys_energy

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (twofl_glm) then
           f(ixO^S,mag(idir))=w(ixO^S,psi_)
        else
           f(ixO^S,mag(idir))=zero
        end if
      else
        f(ixO^S,mag(idir))=w(ixO^S,mom_c(idim))*w(ixO^S,mag(idir))-w(ixO^S,mag(idim))*w(ixO^S,mom_c(idir))

        if (B0field) then
          f(ixO^S,mag(idir))=f(ixO^S,mag(idir))&
                +w(ixO^S,mom_c(idim))*block%B0(ixO^S,idir,idim)&
                -w(ixO^S,mom_c(idir))*block%B0(ixO^S,idim,idim)
        end if

        if (twofl_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (twofl_etah>zero) then
            if (B0field) then
              f(ixO^S,mag(idir)) = f(ixO^S,mag(idir)) &
                   - vHall(ixO^S,idir)*(w(ixO^S,mag(idim))+block%B0(ixO^S,idim,idim)) &
                   + vHall(ixO^S,idim)*(w(ixO^S,mag(idir))+block%B0(ixO^S,idir,idim))
            else
              f(ixO^S,mag(idir)) = f(ixO^S,mag(idir)) &
                   - vHall(ixO^S,idir)*w(ixO^S,mag(idim)) &
                   + vHall(ixO^S,idim)*w(ixO^S,mag(idir))
            end if
          end if
        end if

      end if
    end do

    if (twofl_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixO^S,psi_)  = cmax_global**2*w(ixO^S,mag(idim))
    end if

    if (twofl_Hall) then
      deallocate(vHall)
    end if

    !!neutrals
    call get_rhon_tot(w,x,ixI^L,ixO^L,tmp)
    f(ixO^S,rho_n_)=w(ixO^S,mom_n(idim))*tmp(ixO^S)
    if(phys_energy) then
      pgas(ixO^S) = w(ixO^S, e_n_)
    else
      pgas(ixO^S)=twofl_adiab*tmp(ixO^S)**twofl_gamma
      if(has_equi_pe_n0) then
        pgas(ixO^S)=pgas(ixO^S)-block%equi_vars(ixO^S,equi_pe_n0_,b0i)
      endif
    endif
    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
        !if(idim==idir) then
        !  f(ixO^S,mom_c(idir)) = pgas(ixO^S)
        !else
        !  f(ixO^S,mom_c(idir)) = 0.0d0
        !end if
        !f(ixO^S,mom_c(idir))=f(ixO^S,mom_c(idir))+w(ixO^S,mom_c(idim))*wC(ixO^S,mom_c(idir))
       f(ixO^S, mom_n(idir)) = w(ixO^S,mom_n(idim)) * wC(ixO^S, mom_n(idir))
    end do

    f(ixO^S, mom_n(idim)) = f(ixO^S, mom_n(idim)) + pgas(ixO^S)

    if(phys_energy) then
      !reuse pgas for storing a in the term: div (u_n * a) and make multiplication at the end
      pgas(ixO^S) = wC(ixO^S,e_n_)
      if(.not. phys_internal_e) then
        ! add pressure perturbation
        pgas(ixO^S) = pgas(ixO^S) + w(ixO^S,e_n_)
      endif
      ! add flux of equilibrium internal energy corresponding to pe_n0
      if(has_equi_pe_n0) then
#if !defined(E_RM_W0) || E_RM_W0 == 1
        pgas(ixO^S) = pgas(ixO^S) + block%equi_vars(ixO^S,equi_pe_n0_,idim) * inv_gamma_1
#else
        pgas(ixO^S) = pgas(ixO^S) + block%equi_vars(ixO^S,equi_pe_n0_,idim) * twofl_gamma * inv_gamma_1
#endif
      endif
      ! add u_n * a in the flux
      f(ixO^S, e_n_) = w(ixO^S,mom_n(idim)) * pgas(ixO^S)

      ! Viscosity fluxes - viscInDiv
      !if (hd_viscosity) then
      !  call visc_get_flux_prim(w, x, ixI^L, ixO^L, idim, f, phys_energy)
      !endif
    end if

  end subroutine twofl_get_flux

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine twofl_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active,wCTprim)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    !use mod_gravity, only: gravity_add_source

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active
    double precision, intent(in), optional :: wCTprim(ixI^S,1:nw)

    if (.not. qsourcesplit) then
      ! Source for solving internal energy
      if(phys_internal_e) then
        active = .true.
        call internal_energy_add_source_n(qdt,ixI^L,ixO^L,wCT,w,x)
        call internal_energy_add_source_c(qdt,ixI^L,ixO^L,wCT,w,x,e_c_)
      else 
        if(phys_solve_eaux) then
          call internal_energy_add_source_c(qdt,ixI^L,ixO^L,wCT,w,x,eaux_c_)
        endif
#if !defined(E_RM_W0) || E_RM_W0==1
        ! add -p0 div v source terms when equi are present
        if(has_equi_pe_n0) then
          active = .true.
          call add_pe_n0_divv(qdt,ixI^L,ixO^L,wCT,w,x)
        endif
        if(has_equi_pe_c0) then
          active = .true.
          call add_pe_c0_divv(qdt,ixI^L,ixO^L,wCT,w,x)
        endif
#endif
        if(twofl_eq_energy == EQ_ENERGY_KI) then
          active = .true.
          call add_source_lorentz_work(qdt,ixI^L,ixO^L,w,wCT,x)
        endif
      endif

      ! Source for B0 splitting
      if (B0field) then
        active = .true.
        call add_source_B0split(qdt,ixI^L,ixO^L,wCT,w,x)
      end if

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
      if (abs(twofl_eta)>smalldouble)then
        active = .true.
        call add_source_res2(qdt,ixI^L,ixO^L,wCT,w,x)
      end if

      if (twofl_eta_hyper>0.d0)then
        active = .true.
        call add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
      end if
      !it is not added in a split manner
      if(.not. use_imex_scheme .and. has_collisions()) then
        active = .true.
        call  twofl_explicit_coll_terms_update(qdt,ixI^L,ixO^L,w,wCT,x)
      endif
    
      if(twofl_hyperdiffusivity) then
        active = .true.
        call  add_source_hyperdiffusive(qdt,ixI^L,ixO^L,w,wCT,x)
      endif  

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

    if(twofl_radiative_cooling_c) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,qsourcesplit,active,rc_fl_c)
    end if
    if(twofl_radiative_cooling_n) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,qsourcesplit,active,rc_fl_n)
    end if
!
!    if(twofl_viscosity) then
!      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,&
!           w,x,phys_energy,qsourcesplit,active)
!    end if
!
    if(twofl_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,twofl_eq_energy .eq. EQ_ENERGY_KI .or. phys_total_energy,qsourcesplit,active)
    end if

  end subroutine twofl_add_source

  subroutine add_pe_n0_divv(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: v(ixI^S,1:ndir)

    call twofl_get_v_n(wCT,x,ixI^L,ixI^L,v)
    call add_geom_PdivV(qdt,ixI^L,ixO^L,v,-block%equi_vars(ixI^S,equi_pe_n0_,0),w,x,e_n_)

  end subroutine add_pe_n0_divv

  subroutine add_pe_c0_divv(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: v(ixI^S,1:ndir)

    call twofl_get_v_c(wCT,x,ixI^L,ixI^L,v)
    call add_geom_PdivV(qdt,ixI^L,ixO^L,v,-block%equi_vars(ixI^S,equi_pe_c0_,0),w,x,e_c_)

  end subroutine add_pe_c0_divv

  subroutine add_geom_PdivV(qdt,ixI^L,ixO^L,v,p,w,x,ind)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L,ind
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: p(ixI^S), v(ixI^S,1:ndir), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divv(ixI^S)

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixI^L,ixO^L,divv,sixthorder=.true.)
      else
        call divvector(v,ixI^L,ixO^L,divv,fourthorder=.true.)
      end if
    else
     call divvector(v,ixI^L,ixO^L,divv)
    end if
    w(ixO^S,ind)=w(ixO^S,ind)+qdt*p(ixO^S)*divv(ixO^S)
  end subroutine add_geom_PdivV

  !> Compute the Lorentz force (JxB)
  subroutine get_lorentz(ixI^L,ixO^L,w,JxB)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: JxB(ixI^S,3)
    double precision                :: a(ixI^S,3), b(ixI^S,3), tmp(ixI^S,3)
    integer                         :: idir, idirmin
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3)

    b=0.0d0
    do idir = 1, ndir
      b(ixO^S, idir) = twofl_mag_i_all(w, ixI^L, ixO^L,idir)
    end do

    ! store J current in a
    call get_current(w,ixI^L,ixO^L,idirmin,current)

    a=0.0d0
    do idir=7-2*ndir,3
      a(ixO^S,idir)=current(ixO^S,idir)
    end do

    call cross_product(ixI^L,ixO^L,a,b,JxB)
  end subroutine get_lorentz

  subroutine  add_source_lorentz_work(qdt,ixI^L,ixO^L,w,wCT,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: a(ixI^S,3), b(ixI^S,1:ndir)
    
    call get_lorentz(ixI^L, ixO^L,wCT,a)
    call twofl_get_v_c(wCT,x,ixI^L,ixO^L,b)
    w(ixO^S,e_c_)=w(ixO^S,e_c_)+qdt*sum(a(ixO^S,1:ndir)*b(ixO^S,1:ndir),dim=ndim+1)

  end subroutine  add_source_lorentz_work

  !> Calculate v_n vector
  subroutine twofl_get_v_n(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,ndir)
    double precision              :: rhon(ixI^S)
    integer :: idir

    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)

    do idir=1,ndir
      v(ixO^S,idir) = w(ixO^S, mom_n(idir)) / rhon(ixO^S)
    end do

  end subroutine twofl_get_v_n

  subroutine get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: rhon(ixI^S)
    if(has_equi_rho_n0) then
      rhon(ixO^S) = w(ixO^S,rho_n_) + block%equi_vars(ixO^S,equi_rho_n0_,b0i)
    else  
      rhon(ixO^S) = w(ixO^S,rho_n_) 
    endif

  end subroutine get_rhon_tot

  subroutine twofl_get_pthermal_n(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(out) :: pth(ixI^S)

    integer :: ix^D, iw

    if(phys_energy) then
      if(phys_internal_e) then
        pth(ixO^S)=gamma_1*w(ixO^S,e_n_)
      else
        pth(ixO^S)=gamma_1*(w(ixO^S,e_n_)&
           - twofl_kin_en_n(w,ixI^L,ixO^L))
      end if
      if(has_equi_pe_n0) then
        pth(ixO^S) = pth(ixO^S) + block%equi_vars(ixO^S,equi_pe_n0_,b0i)
      endif
    else
      call get_rhon_tot(w,x,ixI^L,ixO^L,pth)
      pth(ixO^S)=twofl_adiab*pth(ixO^S)**twofl_gamma
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
                " encountered when call twofl_get_pthermal_n"
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

  end subroutine twofl_get_pthermal_n

  subroutine twofl_get_pthermal_n_primitive(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(out) :: pth(ixI^S)

    if(phys_energy) then
      if(has_equi_pe_n0) then
        pth(ixO^S) = w(ixO^S,e_n_) + block%equi_vars(ixO^S,equi_pe_n0_,b0i)
      else
        pth(ixO^S) = w(ixO^S,e_n_) 
      endif
    else
      call get_rhon_tot(w,x,ixI^L,ixO^L,pth)
      pth(ixO^S)=twofl_adiab*pth(ixO^S)**twofl_gamma
    end if
  end subroutine twofl_get_pthermal_n_primitive

  !> Calculate v component
  subroutine twofl_get_v_n_idim(w,x,ixI^L,ixO^L,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S)
    double precision              :: rhon(ixI^S)

    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    v(ixO^S) = w(ixO^S, mom_n(idim)) / rhon(ixO^S)

  end subroutine twofl_get_v_n_idim

  subroutine internal_energy_add_source_n(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: pth(ixI^S),v(ixI^S,1:ndir),divv(ixI^S)

    call twofl_get_pthermal_n(wCT,x,ixI^L,ixO^L,pth)
    call twofl_get_v_n(wCT,x,ixI^L,ixI^L,v)
    call add_geom_PdivV(qdt,ixI^L,ixO^L,v,-pth,w,x,e_n_)

    if(fix_small_values .and. .not. has_equi_pe_n0) then
      call twofl_handle_small_ei_n(w,x,ixI^L,ixO^L,e_n_,'internal_energy_add_source')
    end if
  end subroutine internal_energy_add_source_n

  !> Calculate v_c vector
  subroutine twofl_get_v_c(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,ndir)
    double precision              :: rhoc(ixI^S)
    integer :: idir

    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
    do idir=1,ndir
      v(ixO^S,idir) = w(ixO^S, mom_c(idir)) / rhoc(ixO^S)
    end do

  end subroutine twofl_get_v_c

  subroutine get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: rhoc(ixI^S)
    if(has_equi_rho_c0) then
      rhoc(ixO^S) = w(ixO^S,rho_c_) + block%equi_vars(ixO^S,equi_rho_c0_,b0i)
    else  
      rhoc(ixO^S) = w(ixO^S,rho_c_) 
    endif

  end subroutine get_rhoc_tot

  subroutine twofl_get_pthermal_c(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(out) :: pth(ixI^S)
    integer :: ix^D, iw

    if(phys_energy) then
      if(phys_internal_e) then
        pth(ixO^S)=gamma_1*w(ixO^S,e_c_)
      elseif(phys_total_energy) then
        pth(ixO^S)=gamma_1*(w(ixO^S,e_c_)&
           - twofl_kin_en_c(w,ixI^L,ixO^L)&
           - twofl_mag_en(w,ixI^L,ixO^L))
      else
        pth(ixO^S)=gamma_1*(w(ixO^S,e_c_)&
           - twofl_kin_en_c(w,ixI^L,ixO^L))
      end if
      if(has_equi_pe_c0) then
        pth(ixO^S) = pth(ixO^S) + block%equi_vars(ixO^S,equi_pe_c0_,b0i)
      endif
    else
      call get_rhoc_tot(w,x,ixI^L,ixO^L,pth)
      pth(ixO^S)=twofl_adiab*pth(ixO^S)**twofl_gamma
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
                " encountered when call twofl_get_pe_c1"
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

  end subroutine twofl_get_pthermal_c

  subroutine twofl_get_pthermal_c_primitive(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(out) :: pth(ixI^S)

    if(phys_energy) then
      if(has_equi_pe_c0) then
        pth(ixO^S) = w(ixO^S,e_c_) + block%equi_vars(ixO^S,equi_pe_c0_,b0i)
      else
        pth(ixO^S) = w(ixO^S,e_c_) 
      endif
    else
      call get_rhoc_tot(w,x,ixI^L,ixO^L,pth)
      pth(ixO^S)=twofl_adiab*pth(ixO^S)**twofl_gamma
    end if
  end subroutine twofl_get_pthermal_c_primitive

  !> Calculate v_c component
  subroutine twofl_get_v_c_idim(w,x,ixI^L,ixO^L,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S)
    double precision              :: rhoc(ixI^S)

    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
    v(ixO^S) = w(ixO^S, mom_c(idim)) / rhoc(ixO^S)

  end subroutine twofl_get_v_c_idim

  subroutine internal_energy_add_source_c(qdt,ixI^L,ixO^L,wCT,w,x,ie)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: pth(ixI^S),v(ixI^S,1:ndir),divv(ixI^S)

    call twofl_get_pthermal_c(wCT,x,ixI^L,ixO^L,pth)
    call twofl_get_v_c(wCT,x,ixI^L,ixI^L,v)
    call add_geom_PdivV(qdt,ixI^L,ixO^L,v,-pth,w,x,ie)
    if(fix_small_values .and. .not. has_equi_pe_c0) then
      call twofl_handle_small_ei_c(w,x,ixI^L,ixO^L,ie,'internal_energy_add_source')
    end if
  end subroutine internal_energy_add_source_c

  !> handle small or negative internal energy
  subroutine twofl_handle_small_ei_c(w, x, ixI^L, ixO^L, ie, subname)
    use mod_global_parameters
    use mod_small_values
    integer, intent(in)             :: ixI^L,ixO^L, ie
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixI^S,1:nw)
    double precision              :: rhoc(ixI^S)
    double precision              :: rhon(ixI^S)

    flag=.false.
    if(has_equi_pe_c0) then
      where(w(ixO^S,ie)+block%equi_vars(ixO^S,equi_pe_c0_,0)*inv_gamma_1<small_e)&
             flag(ixO^S,ie)=.true.
    else
      where(w(ixO^S,ie)<small_e) flag(ixO^S,ie)=.true.
    endif  
    if(any(flag(ixO^S,ie))) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_pe_c0) then
          where(flag(ixO^S,ie)) w(ixO^S,ie)=small_e - &
                  block%equi_vars(ixO^S,equi_pe_c0_,0)*inv_gamma_1
        else
          where(flag(ixO^S,ie)) w(ixO^S,ie)=small_e
        endif
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag, ie)
      case default
        ! small values error shows primitive variables
        ! to_primitive subroutine cannot be used as this error handling
        ! is also used in TC where e_to_ei is explicitly called
        w(ixO^S,e_n_)=w(ixO^S,e_n_)*gamma_1
        call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
        w(ixO^S,e_c_)=w(ixO^S,e_c_)*gamma_1
        call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
        do idir = 1, ndir
           w(ixO^S, mom_n(idir)) = w(ixO^S, mom_n(idir))/rhon(ixO^S)
           w(ixO^S, mom_c(idir)) = w(ixO^S, mom_c(idir))/rhoc(ixO^S)
        end do
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if

  end subroutine twofl_handle_small_ei_c

  !> handle small or negative internal energy
  subroutine twofl_handle_small_ei_n(w, x, ixI^L, ixO^L, ie, subname)
    use mod_global_parameters
    use mod_small_values
    integer, intent(in)             :: ixI^L,ixO^L, ie
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixI^S,1:nw)
    double precision              :: rhoc(ixI^S)
    double precision              :: rhon(ixI^S)

    flag=.false.
    if(has_equi_pe_n0) then
      where(w(ixO^S,ie)+block%equi_vars(ixO^S,equi_pe_n0_,0)*inv_gamma_1<small_e)&
                         flag(ixO^S,ie)=.true.
    else
      where(w(ixO^S,ie)<small_e) flag(ixO^S,ie)=.true.
    endif
    if(any(flag(ixO^S,ie))) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_pe_n0) then
          where(flag(ixO^S,ie)) w(ixO^S,ie)=small_e - &
                block%equi_vars(ixO^S,equi_pe_n0_,0)*inv_gamma_1
        else
          where(flag(ixO^S,ie)) w(ixO^S,ie)=small_e
        endif
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag, ie)
      case default
        ! small values error shows primitive variables
        w(ixO^S,e_n_)=w(ixO^S,e_n_)*gamma_1
        call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
        w(ixO^S,e_c_)=w(ixO^S,e_c_)*gamma_1
        call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
        do idir = 1, ndir
           w(ixO^S, mom_n(idir)) = w(ixO^S, mom_n(idir))/rhon(ixO^S)
           w(ixO^S, mom_c(idir)) = w(ixO^S, mom_c(idir))/rhoc(ixO^S)
        end do
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if

  end subroutine twofl_handle_small_ei_n

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
      w(ixO^S,mom_c(1:ndir))=w(ixO^S,mom_c(1:ndir))+axb(ixO^S,1:ndir)
    end if

    if(phys_total_energy) then
      a=0.d0
      ! for free-free field -(vxB0) dot J0 =0
      b(ixO^S,:)=wCT(ixO^S,mag(:))
      ! store full magnetic field B0+B1 in b
      if(.not.B0field_forcefree) b(ixO^S,:)=b(ixO^S,:)+block%B0(ixO^S,:,0)
      ! store velocity in a
      do idir=1,ndir
        call twofl_get_v_c_idim(wCT,x,ixI^L,ixO^L,idir,a(ixI^S,idir))
      end do
      call cross_product(ixI^L,ixO^L,a,b,axb)
      axb(ixO^S,:)=axb(ixO^S,:)*qdt
      ! add -(vxB) dot J0 source term in energy equation
      do idir=7-2*ndir,3
        w(ixO^S,e_c_)=w(ixO^S,e_c_)-axb(ixO^S,idir)*block%J0(ixO^S,idir)
      end do
    end if

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_B0')

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
    if (twofl_4th_order) then
      ixA^L=ixO^L^LADD2;
    else
      ixA^L=ixO^L^LADD1;
    end if

    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixI^L,ixO^L,idirmin,current)

    if (twofl_eta>zero)then
       eta(ixA^S)=twofl_eta
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
       if (twofl_4th_order) then
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
       if (twofl_eta<zero)then
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
       if (phys_energy) then
          w(ixO^S,e_c_)=w(ixO^S,e_c_)+qdt*tmp(ixO^S)*Bf(ixO^S,idir)
          if(phys_solve_eaux) then
            w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)+qdt*tmp(ixO^S)*Bf(ixO^S,idir)
          end if
       end if
    end do ! idir

    if (phys_energy) then
       ! de/dt+=eta*J**2
      tmp(ixO^S)=qdt*eta(ixO^S)*sum(current(ixO^S,:)**2,dim=ndim+1)
      w(ixO^S,e_c_)=w(ixO^S,e_c_)+tmp(ixO^S)
      if(phys_solve_eaux) then
        ! add eta*J**2 source term in the internal energy equation
        w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)+tmp(ixO^S)
      end if
    end if

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res1')

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

    if (twofl_eta>zero)then
       eta(ixA^S)=twofl_eta
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
    if(stagger_grid.and.ndim==2.and.ndir==3) then
      ! if 2.5D
      w(ixO^S,mag(ndir)) = w(ixO^S,mag(ndir))-qdt*curlj(ixO^S,ndir)
    else
      w(ixO^S,mag(1:ndir)) = w(ixO^S,mag(1:ndir))-qdt*curlj(ixO^S,1:ndir)
    end if

    if(phys_energy) then
      ! de/dt= +div(B x Jeta) = eta J^2 - B dot curl(eta J)
      ! de1/dt= eta J^2 - B1 dot curl(eta J)
      tmp(ixO^S)=eta(ixO^S)*sum(current(ixO^S,:)**2,dim=ndim+1)
      w(ixO^S,e_c_)=w(ixO^S,e_c_)+qdt*(tmp(ixO^S)-&
        sum(wCT(ixO^S,mag(1:ndir))*curlj(ixO^S,1:ndir),dim=ndim+1))
      if(phys_solve_eaux) then
        ! add eta*J**2 source term in the internal energy equation
        w(ixO^S,eaux_c_)=w(ixO^S,eaux_c_)+tmp(ixO^S)
      end if
    end if

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res2')
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
    ehyper(ixA^S,1:ndir) = - tmpvec(ixA^S,1:ndir)*twofl_eta_hyper

    ixA^L=ixO^L;
    tmpvec2(ixA^S,1:ndir)=zero
    call curlvector(ehyper,ixI^L,ixA^L,tmpvec2,idirmin1,1,3)

    do idir=1,ndir
      w(ixO^S,mag(idir)) = w(ixO^S,mag(idir))-tmpvec2(ixO^S,idir)*qdt
    end do

    if (phys_energy) then
      ! de/dt= +div(B x Ehyper)
      ixA^L=ixO^L^LADD1;
      tmpvec2(ixA^S,1:ndir)=zero
      do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
        tmpvec2(ixA^S,idir) = tmpvec(ixA^S,idir)&
        + lvc(idir,jdir,kdir)*wCT(ixA^S,mag(jdir))*ehyper(ixA^S,kdir)
      end do; end do; end do
      tmp(ixO^S)=zero
      call divvector(tmpvec2,ixI^L,ixO^L,tmp)
      w(ixO^S,e_c_)=w(ixO^S,e_c_)+tmp(ixO^S)*qdt
    end if

    if (fix_small_values)  call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_hyperres')

  end subroutine add_source_hyperres

  subroutine add_source_glm(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-MHD scheme
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision:: divb(ixI^S)
    integer          :: idim,idir
    double precision :: gradPsi(ixI^S)

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb, twofl_divb_4thorder)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (twofl_glm_alpha < zero) then
      w(ixO^S,psi_) = abs(twofl_glm_alpha)*wCT(ixO^S,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab_uniform) then
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*twofl_glm_alpha/minval(dxlevel(:)))*w(ixO^S,psi_)
      else
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*twofl_glm_alpha/minval(block%ds(ixO^S,:),dim=ndim+1))*w(ixO^S,psi_)
      end if
    end if

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
          call gradient(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
       case("limited")
          call gradientS(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
       end select
       if (phys_total_energy) then
       ! e  = e  -qdt (b . grad(Psi))
         w(ixO^S,e_c_) = w(ixO^S,e_c_)-qdt*wCT(ixO^S,mag(idim))*gradPsi(ixO^S)
       end if
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixO^S,mom_c(idir))=w(ixO^S,mom_c(idir))-qdt*twofl_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_glm')

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
    call get_divb(wCT,ixI^L,ixO^L,divb, twofl_divb_4thorder)

    ! calculate velocity
    call twofl_get_v_c(wCT,x,ixI^L,ixO^L,v)

    if (phys_total_energy) then
      ! e = e - qdt (v . b) * div b
      w(ixO^S,e_c_)=w(ixO^S,e_c_)-&
           qdt*sum(v(ixO^S,:)*wCT(ixO^S,mag(:)),dim=ndim+1)*divb(ixO^S)
    end if

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*v(ixO^S,idir)*divb(ixO^S)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixO^S,mom_c(idir))=w(ixO^S,mom_c(idir))-qdt*twofl_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divb(ixI^S),vel(ixI^S)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb, twofl_divb_4thorder)

    ! b = b - qdt v * div b
    do idir=1,ndir
      call twofl_get_v_c_idim(wCT,x,ixI^L,ixO^L,idir,vel)
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*vel(ixO^S)*divb(ixO^S)
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_janhunen')

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
    call get_divb(wCT,ixI^L,ixp^L,divb, twofl_divb_4thorder)

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

       if (typedivbdiff=='all' .and. phys_total_energy) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixp^S,e_c_)=w(ixp^S,e_c_)+wCT(ixp^S,mag(idim))*graddivb(ixp^S)
       end if
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_linde')

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
    invB(ixO^S)=sqrt(twofl_mag_en_all(w,ixI^L,ixO^L))
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

  ! copied from gravity
  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine gravity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
       energy,qsourcesplit,active)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active
    double precision       :: vel(ixI^S)
    integer                         :: idim

    double precision :: gravity_field(ixI^S,ndim)

    if(qsourcesplit .eqv. grav_split) then
      active = .true.

      if (.not. associated(usr_gravity)) then
        write(*,*) "mod_usr.t: please point usr_gravity to a subroutine"
        write(*,*) "like the phys_gravity in mod_usr_methods.t"
        call mpistop("gravity_add_source: usr_gravity not defined")
      else
        call usr_gravity(ixI^L,ixO^L,wCT,x,gravity_field)
      end if
  
      do idim = 1, ndim
        w(ixO^S,mom_n(idim)) = w(ixO^S,mom_n(idim)) &
              + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,rho_n_)
        w(ixO^S,mom_c(idim)) = w(ixO^S,mom_c(idim)) &
              + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,rho_c_)
        if(energy) then
#if !defined(E_RM_W0) || E_RM_W0 == 1
          call twofl_get_v_n_idim(wCT,x,ixI^L,ixO^L,idim,vel)
          w(ixO^S,e_n_)=w(ixO^S,e_n_) &
              + qdt * gravity_field(ixO^S,idim) * vel(ixO^S) * wCT(ixO^S,rho_n_)
          call twofl_get_v_c_idim(wCT,x,ixI^L,ixO^L,idim,vel)
          w(ixO^S,e_c_)=w(ixO^S,e_c_) &
              + qdt * gravity_field(ixO^S,idim) * vel(ixO^S) * wCT(ixO^S,rho_c_)
#else
          w(ixO^S,e_n_)=w(ixO^S,e_n_) &
              + qdt * gravity_field(ixO^S,idim) *  wCT(ixO^S,mom_n(idim))
          w(ixO^S,e_c_)=w(ixO^S,e_c_) &
              + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,mom_c(idim))
#endif


        end if
      end do
    end if

  end subroutine gravity_add_source

  subroutine gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision                :: dxinv(1:ndim), max_grav
    integer                         :: idim

    double precision :: gravity_field(ixI^S,ndim)

    ^D&dxinv(^D)=one/dx^D;

    if(.not. associated(usr_gravity)) then
      write(*,*) "mod_usr.t: please point usr_gravity to a subroutine"
      write(*,*) "like the phys_gravity in mod_usr_methods.t"
      call mpistop("gravity_get_dt: usr_gravity not defined")
    else
      call usr_gravity(ixI^L,ixO^L,w,x,gravity_field)
    end if

    do idim = 1, ndim
      max_grav = maxval(abs(gravity_field(ixO^S,idim)))
      max_grav = max(max_grav, epsilon(1.0d0))
      dtnew = min(dtnew, 1.0d0 / sqrt(max_grav * dxinv(idim)))
    end do

  end subroutine gravity_get_dt


  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine twofl_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt
    !use mod_viscosity, only: viscosity_get_dt
    !use mod_gravity, only: gravity_get_dt

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
    if (twofl_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/twofl_eta
    else if (twofl_eta<zero)then
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

    if(twofl_eta_hyper>zero) then
      if(slab_uniform) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/twofl_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixO^S,1:ndim))**4/twofl_eta_hyper,dtnew)
      end if
    end if

    ! the timestep related to coll terms: 1/(rho_n rho_c alpha)
    if(dtcollpar>0d0 .and. has_collisions()) then
        call coll_get_dt(w,x,ixI^L,ixO^L,dtnew)
    endif

    if(twofl_radiative_cooling_c) then
      call cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x,rc_fl_c)
    end if
    if(twofl_radiative_cooling_n) then
      call cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x,rc_fl_n)
    end if
!
!    if(twofl_viscosity) then
!      call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
!    end if
!
    if(twofl_gravity) then
      call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if
    if(twofl_hyperdiffusivity) then
      call hyperdiffusivity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if


  end subroutine twofl_get_dt

  pure function has_collisions() result(res)
    logical :: res
    res = .not. twofl_alpha_coll_constant .or. twofl_alpha_coll >0d0
  end function has_collisions

  subroutine coll_get_dt(w,x,ixI^L,ixO^L,dtnew)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: dtnew

    double precision :: rhon(ixI^S), rhoc(ixI^S), alpha(ixI^S)
    double precision, allocatable :: gamma_rec(:^D&), gamma_ion(:^D&)
    double precision :: max_coll_rate

    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)

    call get_alpha_coll(ixI^L, ixO^L, w, x, alpha)
    max_coll_rate = maxval(alpha(ixO^S) * max(rhon(ixO^S), rhoc(ixO^S)))

    if(twofl_coll_inc_ionrec) then
       allocate(gamma_ion(ixI^S), gamma_rec(ixI^S)) 
       call get_gamma_ion_rec(ixI^L, ixO^L, w, x, gamma_rec, gamma_ion)
       max_coll_rate=max(max_coll_rate, maxval(gamma_ion(ixO^S)), maxval(gamma_rec(ixO^S))) 
       deallocate(gamma_ion, gamma_rec) 
    endif
    dtnew = min(dtcollpar/max_coll_rate, dtnew)

  end subroutine coll_get_dt

  ! Add geometrical source terms to w
  subroutine twofl_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),rho(ixI^S)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    ! charges

    mr_=mom_c(1); mphi_=mom_c(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_
    call get_rhoc_tot(wCT,x,ixI^L,ixO^L,rho)

    select case (coordinate)
    case (cylindrical)
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in MHD")
      endif
      call twofl_get_p_c_total(wCT,x,ixI^L,ixO^L,tmp)

      if(phi_>0) then
        w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*(tmp(ixO^S)-&
                  wCT(ixO^S,bphi_)**2+wCT(ixO^S,mphi_)**2/rho(ixO^S))
        w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt/x(ixO^S,1)*(&
                 -wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/rho(ixO^S) &
                 +wCT(ixO^S,bphi_)*wCT(ixO^S,br_))
        if(.not.stagger_grid) then
          w(ixO^S,bphi_)=w(ixO^S,bphi_)+qdt/x(ixO^S,1)*&
                   (wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
                   -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
                   /rho(ixO^S)
        end if
      else
        w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*tmp(ixO^S)
      end if
      if(twofl_glm) w(ixO^S,br_)=w(ixO^S,br_)+qdt*wCT(ixO^S,psi_)/x(ixO^S,1)
    case (spherical)
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       call twofl_get_p_c_total(wCT,x,ixI^L,ixO^L,tmp1)
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
           tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom_c(idir))**2/rho(ixO^S)-wCT(ixO^S,mag(idir))**2
           if(B0field) tmp(ixO^S)=tmp(ixO^S)-2.0d0*block%B0(ixO^S,idir,0)*wCT(ixO^S,mag(idir))
         end do
       end if
       w(ixO^S,mom_c(1))=w(ixO^S,mom_c(1))+qdt*tmp(ixO^S)/x(ixO^S,1)
       ! b1
       if(twofl_glm) then
         w(ixO^S,mag(1))=w(ixO^S,mag(1))+qdt/x(ixO^S,1)*2.0d0*wCT(ixO^S,psi_)
       end if

       {^NOONED
       ! m2
       tmp(ixO^S)=tmp1(ixO^S)
       if(B0field) then
         tmp(ixO^S)=tmp(ixO^S)+tmp2(ixO^S)
       end if
       ! This will make hydrostatic p=const an exact solution
       w(ixO^S,mom_c(2))=w(ixO^S,mom_c(2))+qdt*tmp(ixO^S) &
            *(block%surfaceC(ixO^S,2)-block%surfaceC(h2x^S,2)) &
            /block%dvolume(ixO^S)
       tmp(ixO^S)=-(wCT(ixO^S,mom_c(1))*wCT(ixO^S,mom_c(2))/rho(ixO^S) &
            -wCT(ixO^S,mag(1))*wCT(ixO^S,mag(2)))
       if (B0field) then
          tmp(ixO^S)=tmp(ixO^S)+block%B0(ixO^S,1,0)*wCT(ixO^S,mag(2)) &
               +wCT(ixO^S,mag(1))*block%B0(ixO^S,2,0)
       end if
       if(ndir==3) then
         tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom_c(3))**2/rho(ixO^S) &
              -wCT(ixO^S,mag(3))**2)*dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)-2.0d0*block%B0(ixO^S,3,0)*wCT(ixO^S,mag(3))&
                 *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
         end if
       end if
       w(ixO^S,mom_c(2))=w(ixO^S,mom_c(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
       ! b2
       if(.not.stagger_grid) then
         tmp(ixO^S)=(wCT(ixO^S,mom_c(1))*wCT(ixO^S,mag(2)) &
              -wCT(ixO^S,mom_c(2))*wCT(ixO^S,mag(1)))/rho(ixO^S)
         if(B0field) then
           tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom_c(1))*block%B0(ixO^S,2,0) &
                -wCT(ixO^S,mom_c(2))*block%B0(ixO^S,1,0))/rho(ixO^S)
         end if
         if(twofl_glm) then
           tmp(ixO^S)=tmp(ixO^S) &
                + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
         end if
         w(ixO^S,mag(2))=w(ixO^S,mag(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
       end if
       }

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixO^S)=-(wCT(ixO^S,mom_c(3))*wCT(ixO^S,mom_c(1))/rho(ixO^S) &
                -wCT(ixO^S,mag(3))*wCT(ixO^S,mag(1))) {^NOONED &
                -(wCT(ixO^S,mom_c(2))*wCT(ixO^S,mom_c(3))/rho(ixO^S) &
                -wCT(ixO^S,mag(2))*wCT(ixO^S,mag(3))) &
                *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
           if (B0field) then
              tmp(ixO^S)=tmp(ixO^S)+block%B0(ixO^S,1,0)*wCT(ixO^S,mag(3)) &
                   +wCT(ixO^S,mag(1))*block%B0(ixO^S,3,0) {^NOONED &
                   +(block%B0(ixO^S,2,0)*wCT(ixO^S,mag(3)) &
                   +wCT(ixO^S,mag(2))*block%B0(ixO^S,3,0)) &
                   *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
           end if
           w(ixO^S,mom_c(3))=w(ixO^S,mom_c(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
         else
           call mpistop("angmomfix not implemented yet in MHD")
         end if
         ! b3
         if(.not.stagger_grid) then
           tmp(ixO^S)=(wCT(ixO^S,mom_c(1))*wCT(ixO^S,mag(3)) &
                -wCT(ixO^S,mom_c(3))*wCT(ixO^S,mag(1)))/rho(ixO^S) {^NOONED &
                -(wCT(ixO^S,mom_c(3))*wCT(ixO^S,mag(2)) &
                -wCT(ixO^S,mom_c(2))*wCT(ixO^S,mag(3)))*dcos(x(ixO^S,2)) &
                /(rho(ixO^S)*dsin(x(ixO^S,2))) }
           if (B0field) then
              tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom_c(1))*block%B0(ixO^S,3,0) &
                   -wCT(ixO^S,mom_c(3))*block%B0(ixO^S,1,0))/rho(ixO^S){^NOONED &
                   -(wCT(ixO^S,mom_c(3))*block%B0(ixO^S,2,0) &
                   -wCT(ixO^S,mom_c(2))*block%B0(ixO^S,3,0))*dcos(x(ixO^S,2)) &
                   /(rho(ixO^S)*dsin(x(ixO^S,2))) }
           end if
           w(ixO^S,mag(3))=w(ixO^S,mag(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
         end if
       end if
    end select

    ! neutrals 
    !TODO no dust: see and implement them from hd/mod_hd_phys !
    !uncomment cartesian expansion
    call get_rhon_tot(wCT,x,ixI^L,ixO^L,rho)
    call twofl_get_pthermal_n(wCT, x, ixI^L, ixO^L, tmp1)

    select case (coordinate)
!    case(Cartesian_expansion)
!      !the user provides the functions of exp_factor and del_exp_factor
!      if(associated(usr_set_surface)) call usr_set_surface(ixI^L,x,block%dx,exp_factor,del_exp_factor,exp_factor_primitive)
!      tmp(ixO^S) = tmp1(ixO^S)*del_exp_factor(ixO^S)/exp_factor(ixO^S)
!      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt*tmp(ixO^S)

    case (cylindrical)
      mr_   = mom_n(r_)
      if (phi_ > 0) then
         where (rho(ixO^S) > 0d0)
            tmp(ixO^S) = tmp1(ixO^S) + wCT(ixO^S, mphi_)**2 / rho(ixO^S)
            w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * tmp(ixO^S) / x(ixO^S, r_)
         end where
         ! s[mphi]=(-mphi*mr/rho)/radius
         if(.not. angmomfix) then
            where (rho(ixO^S) > 0d0)
               tmp(ixO^S) = -wCT(ixO^S, mphi_) * wCT(ixO^S, mr_) / rho(ixO^S)
               w(ixO^S, mphi_) = w(ixO^S, mphi_) + qdt * tmp(ixO^S) / x(ixO^S, r_)
            end where
         end if
      else
         ! s[mr]=2pthermal/radius
         w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * tmp1(ixO^S) / x(ixO^S, r_)
      end if
    case (spherical)
       if(phi_>0) mphi_ = mom_n(phi_)
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       ! s[mr]=((mtheta**2+mphi**2)/rho+2*p)/r
       tmp(ixO^S) = tmp1(ixO^S) * x(ixO^S, 1) &
            *(block%surfaceC(ixO^S, 1) - block%surfaceC(h1x^S, 1)) &
            /block%dvolume(ixO^S)
       if (ndir > 1) then
         do idir = 2, ndir
           tmp(ixO^S) = tmp(ixO^S) + wCT(ixO^S, mom_n(idir))**2 / rho(ixO^S)
         end do
       end if
       w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * tmp(ixO^S) / x(ixO^S, 1)

       {^NOONED
       ! s[mtheta]=-(mr*mtheta/rho)/r+cot(theta)*(mphi**2/rho+p)/r
       tmp(ixO^S) = tmp1(ixO^S) * x(ixO^S, 1) &
            * (block%surfaceC(ixO^S, 2) - block%surfaceC(h2x^S, 2)) &
            / block%dvolume(ixO^S)
       if (ndir == 3) then
          tmp(ixO^S) = tmp(ixO^S) + (wCT(ixO^S, mom_n(3))**2 / rho(ixO^S)) / tan(x(ixO^S, 2))
       end if
       if (.not. angmomfix) then
          tmp(ixO^S) = tmp(ixO^S) - (wCT(ixO^S, mom_n(2)) * wCT(ixO^S, mr_)) / rho(ixO^S)
       end if
       w(ixO^S, mom_n(2)) = w(ixO^S, mom_n(2)) + qdt * tmp(ixO^S) / x(ixO^S, 1)

       if (ndir == 3) then
         ! s[mphi]=-(mphi*mr/rho)/r-cot(theta)*(mtheta*mphi/rho)/r
         if (.not. angmomfix) then
           tmp(ixO^S) = -(wCT(ixO^S, mom_n(3)) * wCT(ixO^S, mr_)) / rho(ixO^S)&
                      - (wCT(ixO^S, mom_n(2)) * wCT(ixO^S, mom_n(3))) / rho(ixO^S) / tan(x(ixO^S, 2))
           w(ixO^S, mom_n(3)) = w(ixO^S, mom_n(3)) + qdt * tmp(ixO^S) / x(ixO^S, 1)
         end if
       end if
       }
    end select

!    if (hd_viscosity) call visc_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
!
!    if (hd_rotating_frame) then
!       if (hd_dust) then
!          call mpistop("Rotating frame not implemented yet with dust")
!       else
!          call rotating_frame_add_source(qdt,ixI^L,ixO^L,wCT,w,x)
!       end if
!    end if
!

    contains
      subroutine twofl_get_p_c_total(w,x,ixI^L,ixO^L,p)
        use mod_global_parameters
    
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: w(ixI^S,nw)
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(out)   :: p(ixI^S)
    
        call twofl_get_pthermal_c(w,x,ixI^L,ixO^L,p)
    
        p(ixO^S) = p(ixO^S) + 0.5d0 * sum(w(ixO^S, mag(:))**2, dim=ndim+1)
    
      end subroutine twofl_get_p_c_total

  end subroutine twofl_add_source_geom

  subroutine twofl_get_temp_c_pert_from_etot(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)

    ! store pe1 in res
    res(ixO^S)=(gamma_1*(w(ixO^S,e_c_)&
           - twofl_kin_en_c(w,ixI^L,ixO^L)&
           - twofl_mag_en(w,ixI^L,ixO^L)))
    if(has_equi_pe_c0) then
      res(ixO^S) = res(ixO^S) + block%equi_vars(ixO^S,equi_pe_c0_,b0i)
      if(has_equi_rho_c0) then
        res(ixO^S) = res(ixO^S)/(Rc * (w(ixO^S,rho_c_)+ block%equi_vars(ixO^S,equi_rho_c0_,b0i))) - &
                      block%equi_vars(ixO^S,equi_pe_c0_,b0i)/(Rc * block%equi_vars(ixO^S,equi_rho_c0_,b0i))
      else
        ! infinite equi temperature with p0 and 0 density
        res(ixO^S) = 0d0
      endif
    else
      res(ixO^S) = res(ixO^S)/(Rc * w(ixO^S,rho_c_))
    endif

  end subroutine twofl_get_temp_c_pert_from_etot

  !> Compute 2 times total magnetic energy
  function twofl_mag_en_all(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    if (B0field) then
      mge(ixO^S) = sum((w(ixO^S, mag(:))+block%B0(ixO^S,:,b0i))**2, dim=ndim+1)
    else
      mge(ixO^S) = sum(w(ixO^S, mag(:))**2, dim=ndim+1)
    end if
  end function twofl_mag_en_all

  !> Compute full magnetic field by direction
  function twofl_mag_i_all(w, ixI^L, ixO^L,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idir
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mgf(ixO^S)

    if (B0field) then
      mgf(ixO^S) = w(ixO^S, mag(idir))+block%B0(ixO^S,idir,b0i)
    else
      mgf(ixO^S) = w(ixO^S, mag(idir))
    end if
  end function twofl_mag_i_all

  !> Compute evolving magnetic energy
  function twofl_mag_en(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    mge(ixO^S) = 0.5d0 * sum(w(ixO^S, mag(:))**2, dim=ndim+1)
  end function twofl_mag_en

  !> compute kinetic energy of neutrals
  function twofl_kin_en_n(w, ixI^L, ixO^L) result(ke)
    use mod_global_parameters, only: nw, ndim,block
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)

    if(has_equi_rho_n0) then
      ke(ixO^S) = 0.5d0 * sum(w(ixO^S, mom_n(:))**2, dim=ndim+1) / (w(ixO^S, rho_n_) + block%equi_vars(ixO^S,equi_rho_n0_,0))
    else
      ke(ixO^S) = 0.5d0 * sum(w(ixO^S, mom_n(:))**2, dim=ndim+1) / w(ixO^S, rho_n_)
    endif

  end function twofl_kin_en_n

  subroutine twofl_get_temp_n_pert_from_etot(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)

    ! store pe1 in res
    res(ixO^S)=(gamma_1*(w(ixO^S,e_c_)- twofl_kin_en_c(w,ixI^L,ixO^L)))
    if(has_equi_pe_n0) then
      res(ixO^S) = res(ixO^S) + block%equi_vars(ixO^S,equi_pe_n0_,b0i)
      if(has_equi_rho_n0) then
        res(ixO^S) = res(ixO^S)/(Rn * (w(ixO^S,rho_n_)+ block%equi_vars(ixO^S,equi_rho_n0_,b0i))) - &
                      block%equi_vars(ixO^S,equi_pe_n0_,b0i)/(Rn * block%equi_vars(ixO^S,equi_rho_n0_,b0i))
      else
        ! infinite equi temperature with p0 and 0 density
        res(ixO^S) = 0d0
      endif
    else
      res(ixO^S) = res(ixO^S)/(Rn * w(ixO^S,rho_n_))
    endif

  end subroutine twofl_get_temp_n_pert_from_etot

  !> compute kinetic energy of charges
  !> w are conserved variables
  function twofl_kin_en_c(w, ixI^L, ixO^L) result(ke)
    use mod_global_parameters, only: nw, ndim,block
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)

    if(has_equi_rho_c0) then
      ke(ixO^S) = 0.5d0 * sum(w(ixO^S, mom_c(:))**2, dim=ndim+1) / (w(ixO^S, rho_c_) + block%equi_vars(ixO^S,equi_rho_c0_,0))
    else
      ke(ixO^S) = 0.5d0 * sum(w(ixO^S, mom_c(:))**2, dim=ndim+1) / w(ixO^S, rho_c_)
    endif
  end function twofl_kin_en_c

  subroutine twofl_getv_Hall(w,x,ixI^L,ixO^L,vHall)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: vHall(ixI^S,1:3)

    integer          :: idir, idirmin
    double precision :: current(ixI^S,7-2*ndir:3)
    double precision :: rho(ixI^S)

    call get_rhoc_tot(w,x,ixI^L,ixO^L,rho)
    ! Calculate current density and idirmin
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    vHall(ixO^S,1:3) = zero
    vHall(ixO^S,idirmin:3) = - twofl_etah*current(ixO^S,idirmin:3)
    do idir = idirmin, 3
       vHall(ixO^S,idir) = vHall(ixO^S,idir)/rho(ixO^S)
    end do

  end subroutine twofl_getv_Hall

! the following not used
!  subroutine twofl_getdt_Hall(w,x,ixI^L,ixO^L,dx^D,dthall)
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
!    ! because we have that in cmax now:
!    return
!
!    ^D&dxarr(^D)=dx^D;
!
!    if (.not. B0field) then
!       bmag(ixO^S)=sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1))
!       bmag(ixO^S)=sqrt(sum((w(ixO^S,mag(:)) + block%B0(ixO^S,1:ndir,b0i))**2))
!    end if
!
!    if(slab_uniform) then
!      dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(twofl_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_c_)))
!    else
!      dthall=dtdiffpar*minval(block%ds(ixO^S,1:ndim))**2.0d0/(twofl_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_c_)))
!    end if
!
!  end subroutine twofl_getdt_Hall

  subroutine twofl_modify_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)
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

      if(phys_total_energy) then
        wRC(ixO^S,e_c_)=wRC(ixO^S,e_c_)-half*wRC(ixO^S,mag(idir))**2
        wLC(ixO^S,e_c_)=wLC(ixO^S,e_c_)-half*wLC(ixO^S,mag(idir))**2
      end if
      wRC(ixO^S,mag(idir)) = wLp(ixO^S,mag(idir))
      wRC(ixO^S,psi_) = wLp(ixO^S,psi_)
      wLC(ixO^S,mag(idir)) = wLp(ixO^S,mag(idir))
      wLC(ixO^S,psi_) = wLp(ixO^S,psi_)
      ! modify total energy according to the change of magnetic field
      if(phys_total_energy) then
        wRC(ixO^S,e_c_)=wRC(ixO^S,e_c_)+half*wRC(ixO^S,mag(idir))**2
        wLC(ixO^S,e_c_)=wLC(ixO^S,e_c_)+half*wLC(ixO^S,mag(idir))**2
      end if
    end if

    if(associated(usr_set_wLR)) call usr_set_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)

  end subroutine twofl_modify_wLR

  subroutine twofl_boundary_adjust(igrid,psb)
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

  end subroutine twofl_boundary_adjust

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
     case(2)
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
     case(3)
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
     case(4)
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
     {^IFTHREED
     case(5)
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
     case(6)
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
     }
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine fixdivB_boundary

  {^NOONED
  subroutine twofl_clean_divb_multigrid(qdt, qt, active)
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
          print *, "divb_multigrid warning: unknown b.c.: ", &
               typeboundary(mag(idim), n)
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
            twofl_divb_4thorder)
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
         call twofl_face_to_center(ixM^LL,ps(igrid))
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

       if(phys_total_energy) then
         ! Determine magnetic energy difference
         tmp(ixM^T) = 0.5_dp * (sum(ps(igrid)%w(ixM^T, &
              mag(1:ndim))**2, dim=ndim+1) - tmp(ixM^T))
         ! Keep thermal pressure the same
         ps(igrid)%w(ixM^T, e_c_) = ps(igrid)%w(ixM^T, e_c_) + tmp(ixM^T)
       end if
    end do

    active = .true.

  end subroutine twofl_clean_divb_multigrid
  }

  subroutine twofl_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s,vcts)
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

  end subroutine twofl_update_faces

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
    double precision, dimension(ixI^S,7-2*ndim:3) :: E_resi

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D-1;

    ! if there is resistivity, get eta J
    if(twofl_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)

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
            if(twofl_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+E_resi(ixC^S,idir)
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
         ECC(ixI^S,idir)=ECC(ixI^S,idir)+Btot(ixI^S,idim1)*wp(ixI^S,mom_c(idim2))
      else if(lvc(idim1,idim2,idir)==-1) then
         ECC(ixI^S,idir)=ECC(ixI^S,idir)-Btot(ixI^S,idim1)*wp(ixI^S,mom_c(idim2))
      endif
    enddo; enddo; enddo

    ! if there is resistivity, get eta J
    if(twofl_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)
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

            ! add current component of electric field at cell edges E=-vxB+eta J
            if(twofl_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+E_resi(ixC^S,idir)
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
    double precision                   :: btilL(s%ixGs^S,ndim)
    double precision                   :: btilR(s%ixGs^S,ndim)
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
    if(twofl_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)
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

      ! add current component of electric field at cell edges E=-vxB+eta J
      if(twofl_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+E_resi(ixC^S,idir)
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
    if(twofl_eta>zero)then
      jce(ixI^S,:)=jce(ixI^S,:)*twofl_eta
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

  !> calculate cell-center values from face-center values
  subroutine twofl_face_to_center(ixO^L,s)
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

  end subroutine twofl_face_to_center

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

  subroutine  hyperdiffusivity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_hyperdiffusivity
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: dx^D
    double precision, intent(inout) :: dtnew

    double precision :: nu(ixI^S),tmp(ixI^S),rho(ixI^S),temp(ixI^S)
    double precision :: divv(ixI^S,1:ndim)
    double precision :: vel(ixI^S,1:ndir)
    double precision :: csound(ixI^S),csound_dim(ixI^S,1:ndim)
    double precision              :: dxarr(ndim)
    double precision :: maxCoef
    integer :: ixOO^L, hxb^L, hx^L, ii, jj


    ^D&dxarr(^D)=dx^D;
    maxCoef = smalldouble

    ! charges
    call twofl_get_v_c(w,x,ixI^L,ixI^L,vel)
    call get_rhoc_tot(w,x,ixI^L,ixI^L,rho)
    call twofl_get_csound2_c_from_conserved(w,x,ixI^L,ixI^L,csound)
    csound(ixI^S) = sqrt(csound(ixI^S)) + sqrt(twofl_mag_en_all(w,ixI^L,ixI^L) /rho(ixI^S))
    csound(ixI^S) = csound(ixI^S) + sqrt(sum(vel(ixI^S,1:ndir)**2 ,dim=ndim+1))
    do ii=1,ndim
      call div_vel_coeff(ixI^L, ixOO^L, vel, ii, divv(ixI^S,ii))
      hxmin^D=ixImin^D+1;
      hxmax^D=ixImax^D-1;
      hxb^L=hx^L-kr(ii,^D);
      csound_dim(hx^S,ii) = (csound(hxb^S)+csound(hx^S))/2d0
    enddo
    call twofl_get_temp_c_pert_from_etot(w, x, ixI^L, ixI^L, temp)
    do ii=1,ndim
      !TODO the following is copied
      !rho_c
      call hyp_coeff(ixI^L, ixOO^L, w(ixI^S,rho_c_), ii, tmp(ixI^S))
      nu(ixO^S) = c_hyp(rho_c_) * csound_dim(ixO^S,ii) * dxlevel(ii) *  tmp(ixO^S) + &
                   c_shk(rho_c_) * (dxlevel(ii)**2) *divv(ixO^S,ii)
      maxCoef = max(maxCoef,maxval(nu(ixO^S)))

      !TH c  
      call hyp_coeff(ixI^L, ixOO^L, temp(ixI^S), ii, tmp(ixI^S))
      nu(ixO^S) = c_hyp(e_c_) * csound_dim(ixO^S,ii) * dxlevel(ii) *  tmp(ixO^S) + &
                   c_shk(e_c_) * (dxlevel(ii)**2) *divv(ixO^S,ii)
      nu(ixO^S) = nu(ixO^S) * rho(ixO^S) * Rc/(twofl_gamma-1d0)
      maxCoef = max(maxCoef,maxval(nu(ixO^S)))

      !visc c
      do jj=1,ndir
        call hyp_coeff(ixI^L, ixOO^L, vel(ixI^S,jj), ii, tmp(ixI^S))
        nu(ixO^S) = c_hyp(mom_c(jj)) * csound_dim(ixO^S,ii) * dxlevel(ii) *  tmp(ixO^S) + &
                   c_shk(mom_c(jj)) * (dxlevel(ii)**2) *divv(ixO^S,ii)
        nu(ixO^S) = nu(ixO^S) * rho(ixO^S) 
        maxCoef = max(maxCoef,maxval(nu(ixO^S)))
      enddo

      ! Ohmic
      do jj=1,ndir
        if(ii .ne. jj) then
          call hyp_coeff(ixI^L, ixOO^L, w(ixI^S,mag(jj)), ii, tmp(ixI^S))
          nu(ixO^S) = c_hyp(mag(jj)) * csound_dim(ixO^S,ii) * dxlevel(ii) *  tmp(ixO^S) + &
                     c_shk(mag(jj)) * (dxlevel(ii)**2) *divv(ixO^S,ii)
          maxCoef = max(maxCoef,maxval(nu(ixO^S)))
        endif
      enddo

    enddo 
  
      !TODO the following is copied, as charges, and as in add_source!
    ! neutrals
    call twofl_get_v_n(w,x,ixI^L,ixI^L,vel)
    call  twofl_get_csound_n(w,x,ixI^L,ixI^L,csound)
    csound(ixI^S) = csound(ixI^S) + sqrt(sum(vel(ixI^S,1:ndir)**2 ,dim=ndim+1))
    do ii=1,ndim
      call div_vel_coeff(ixI^L, ixOO^L, vel, ii, divv(ixI^S,ii))
      hxmin^D=ixImin^D+1;
      hxmax^D=ixImax^D-1;
      hxb^L=hx^L-kr(ii,^D);
      csound_dim(hx^S,ii) = (csound(hxb^S)+csound(hx^S))/2d0
    enddo
    call get_rhon_tot(w,x,ixI^L,ixO^L,rho)
    call twofl_get_temp_n_pert_from_etot(w, x, ixI^L, ixI^L, temp)
    do ii=1,ndim
      !rho_n
      call hyp_coeff(ixI^L, ixOO^L, w(ixI^S,rho_n_), ii, tmp(ixI^S))
      nu(ixO^S) = c_hyp(rho_n_) * csound_dim(ixO^S,ii) * dxlevel(ii) *  tmp(ixO^S) + &
                   c_shk(rho_n_) * (dxlevel(ii)**2) *divv(ixOO^S,ii)
      maxCoef = max(maxCoef,maxval(nu(ixO^S)))

      !TH n  
      call hyp_coeff(ixI^L, ixOO^L, temp(ixI^S), ii, tmp(ixI^S))
      nu(ixO^S) = c_hyp(e_n_) * csound_dim(ixO^S,ii) * dxlevel(ii) *  tmp(ixO^S) + &
                   c_shk(e_n_) * (dxlevel(ii)**2) *divv(ixO^S,ii)
      nu(ixO^S) = nu(ixO^S) * rho(ixO^S) * Rn/(twofl_gamma-1d0)
      maxCoef = max(maxCoef,maxval(nu(ixO^S)))

      !visc n
      do jj=1,ndir
        call hyp_coeff(ixI^L, ixOO^L, vel(ixI^S,jj), ii, tmp(ixI^S))
        nu(ixO^S) = c_hyp(mom_n(jj)) * csound_dim(ixO^S,ii) * dxlevel(ii) *  tmp(ixO^S) + &
                   c_shk(mom_n(jj)) * (dxlevel(ii)**2) *divv(ixO^S,ii)
        nu(ixO^S) = nu(ixO^S) * rho(ixO^S) 
        maxCoef = max(maxCoef,maxval(nu(ixO^S)))
      enddo
    enddo 

    dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**2/maxCoef,dtnew)
  end subroutine hyperdiffusivity_get_dt

  subroutine  add_source_hyperdiffusive(qdt,ixI^L,ixO^L,w,wCT,x)
    use mod_global_parameters
    use mod_hyperdiffusivity

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: wCT(ixI^S,1:nw)

    double precision :: divv(ixI^S,1:ndim)
    double precision :: vel(ixI^S,1:ndir)
    double precision :: csound(ixI^S),csound_dim(ixI^S,1:ndim)
    integer :: ii,ixOO^L,hxb^L,hx^L
    double precision :: rho(ixI^S)

    call twofl_get_v_c(wCT,x,ixI^L,ixI^L,vel)
    call get_rhoc_tot(wCT,x,ixI^L,ixI^L,rho)
    call twofl_get_csound2_c_from_conserved(wCT,x,ixI^L,ixI^L,csound)
    csound(ixI^S) = sqrt(csound(ixI^S)) + sqrt(twofl_mag_en_all(wCT,ixI^L,ixI^L) /rho(ixI^S))
    csound(ixI^S) = csound(ixI^S) + sqrt(sum(vel(ixI^S,1:ndir)**2 ,dim=ndim+1))
    do ii=1,ndim
      call div_vel_coeff(ixI^L, ixOO^L, vel, ii, divv(ixI^S,ii))
      hxmin^D=ixImin^D+1;
      hxmax^D=ixImax^D-1;
      hxb^L=hx^L-kr(ii,^D);
      csound_dim(hx^S,ii) = (csound(hxb^S)+csound(hx^S))/2d0
    enddo
    call add_density_hyper_Source(rho_c_)
    call add_viscosity_hyper_Source(rho,mom_c(1), e_c_)
    call add_th_cond_c_hyper_Source(rho)
    call add_ohmic_hyper_Source()

    call twofl_get_v_n(wCT,x,ixI^L,ixI^L,vel)
    call  twofl_get_csound_n(wCT,x,ixI^L,ixI^L,csound)
    csound(ixI^S) = csound(ixI^S) + sqrt(sum(vel(ixI^S,1:ndir)**2 ,dim=ndim+1))
    do ii=1,ndim
      call div_vel_coeff(ixI^L, ixOO^L, vel, ii, divv(ixI^S,ii))
      hxmin^D=ixImin^D+1;
      hxmax^D=ixImax^D-1;
      hxb^L=hx^L-kr(ii,^D);
      csound_dim(hx^S,ii) = (csound(hxb^S)+csound(hx^S))/2d0
    enddo
    call add_density_hyper_Source(rho_n_)
    call get_rhon_tot(wCT,x,ixI^L,ixI^L,rho)
    call add_viscosity_hyper_Source(rho,mom_n(1), e_n_)
    call add_th_cond_n_hyper_Source(rho)

    contains

    subroutine add_density_hyper_Source(index_rho)
      integer, intent(in) :: index_rho

      double precision :: nu(ixI^S), tmp(ixI^S)

      do ii=1,ndim
        call hyp_coeff(ixI^L, ixOO^L, wCT(ixI^S,index_rho), ii, tmp(ixI^S))
        nu(ixOO^S) = c_hyp(index_rho) * csound_dim(ixOO^S,ii) * dxlevel(ii) *  tmp(ixOO^S) + &
                     c_shk(index_rho) * (dxlevel(ii)**2) *divv(ixOO^S,ii)
        !print*, "IXOO HYP ", ixOO^L, " IDIMM ", ii
        call second_same_deriv(ixI^L, ixOO^L, nu(ixI^S), wCT(ixI^S,index_rho), ii, tmp)

        w(ixO^S,index_rho) = w(ixO^S,index_rho) + qdt * tmp(ixO^S)
        !print*, "RHO ", index_rho, maxval(abs(tmp(ixO^S)))
      enddo
    end subroutine add_density_hyper_Source   

    subroutine add_th_cond_c_hyper_Source(var2)
      double precision, intent(in) :: var2(ixI^S)
      double precision :: nu(ixI^S), tmp(ixI^S), var(ixI^S)
      call twofl_get_temp_c_pert_from_etot(wCT, x, ixI^L, ixI^L, var)
      do ii=1,ndim
        call hyp_coeff(ixI^L, ixOO^L, var(ixI^S), ii, tmp(ixI^S))
        nu(ixOO^S) = c_hyp(e_c_) * csound_dim(ixOO^S,ii) * dxlevel(ii) *  tmp(ixOO^S) + &
                     c_shk(e_c_) * (dxlevel(ii)**2) *divv(ixOO^S,ii)
        call second_same_deriv2(ixI^L, ixOO^L, nu(ixI^S), var2(ixI^S) ,var(ixI^S), ii, tmp)
        w(ixO^S,e_c_) = w(ixO^S,e_c_) + qdt * tmp(ixO^S) * Rc/(twofl_gamma-1d0)
      !print*, "TH C ", maxval(abs(tmp(ixO^S)))
      enddo
    end subroutine add_th_cond_c_hyper_Source   

    subroutine add_th_cond_n_hyper_Source(var2)
      double precision, intent(in) :: var2(ixI^S)
      double precision :: nu(ixI^S), tmp(ixI^S), var(ixI^S)
      call twofl_get_temp_n_pert_from_etot(wCT, x, ixI^L, ixI^L, var)
      do ii=1,ndim
        call hyp_coeff(ixI^L, ixOO^L, var(ixI^S), ii, tmp(ixI^S))
        nu(ixOO^S) = c_hyp(e_n_) * csound_dim(ixOO^S,ii) * dxlevel(ii) *  tmp(ixOO^S) + &
                     c_shk(e_n_) * (dxlevel(ii)**2) *divv(ixOO^S,ii)
        call second_same_deriv2(ixI^L, ixOO^L, nu(ixI^S), var2(ixI^S) ,var(ixI^S), ii, tmp)
        w(ixO^S,e_n_) = w(ixO^S,e_n_) + qdt * tmp(ixO^S) * Rn/(twofl_gamma-1d0)
      !print*, "TH N ", maxval(abs(tmp(ixO^S)))
      enddo
    end subroutine add_th_cond_n_hyper_Source   

    subroutine add_viscosity_hyper_Source(rho,index_mom1, index_e)
      double precision, intent(in) :: rho(ixI^S)
      integer, intent(in) :: index_mom1, index_e

      double precision :: nu(ixI^S,1:ndir,1:ndim), tmp(ixI^S),tmp2(ixI^S)
      integer :: jj

      do jj=1,ndir
        do ii=1,ndim
          call hyp_coeff(ixI^L, ixOO^L, vel(ixI^S,jj), ii, tmp(ixI^S))
          nu(ixOO^S,jj,ii) = c_hyp(index_mom1-1+jj) * csound_dim(ixOO^S,ii) * dxlevel(ii) *  tmp(ixOO^S) + &
                     c_shk(index_mom1-1+jj) * (dxlevel(ii)**2) *divv(ixOO^S,ii)
        enddo
      enddo
      
      do jj=1,ndir
        do ii=1,ndim
          call second_same_deriv2(ixI^L, ixOO^L, nu(ixI^S,jj,ii), rho(ixI^S), vel(ixI^S,jj), ii, tmp)
          call second_same_deriv2(ixI^L, ixOO^L, nu(ixI^S,jj,ii), wCT(ixI^S,index_mom1-1+jj), vel(ixI^S,jj), ii, tmp2)
          if(ii .eq. jj) then
            w(ixO^S,index_mom1-1+jj) = w(ixO^S,index_mom1-1+jj) + qdt * tmp(ixO^S)
            w(ixO^S,index_e) = w(ixO^S,index_e) + qdt * tmp2(ixO^S)

          else
            w(ixO^S,index_mom1-1+jj) = w(ixO^S,index_mom1-1+jj) + 0.5*qdt * tmp(ixO^S)
            w(ixO^S,index_e) = w(ixO^S,index_e) + 0.5*qdt * tmp2(ixO^S)
            call second_cross_deriv2(ixI^L, ixOO^L, nu(ixI^S,ii,jj), rho(ixI^S), vel(ixI^S,ii), jj, ii, tmp)
            w(ixO^S,index_mom1-1+jj) = w(ixO^S,index_mom1-1+jj) + 0.5*qdt * tmp(ixO^S)
            call second_cross_deriv2(ixI^L, ixOO^L, nu(ixI^S,jj,ii), wCT(ixI^S,index_mom1-1+jj), vel(ixI^S,jj), ii, jj, tmp2)
            w(ixO^S,index_e) = w(ixO^S,index_e) + 0.5*qdt * tmp2(ixO^S)
          endif

        enddo
      enddo

    end subroutine add_viscosity_hyper_Source   

    subroutine add_ohmic_hyper_Source()
      double precision :: nu(ixI^S,1:ndir,1:ndim), tmp(ixI^S)
      integer :: jj

      do jj=1,ndir
        do ii=1,ndim
          if(ii .ne. jj) then
            call hyp_coeff(ixI^L, ixOO^L, wCT(ixI^S,mag(jj)), ii, tmp(ixI^S))
            nu(ixOO^S,jj,ii) = c_hyp(mag(jj)) * csound_dim(ixOO^S,ii) * dxlevel(ii) *  tmp(ixOO^S) + &
                       c_shk(mag(jj)) * (dxlevel(ii)**2) *divv(ixOO^S,ii)
          endif
        enddo
      enddo
      
      do jj=1,ndir
        do ii=1,ndim
          if(ii .ne. jj) then
            !mag field
            call second_same_deriv(ixI^L, ixOO^L, nu(ixI^S,jj,ii), wCT(ixI^S,mag(jj)), ii, tmp)
            w(ixO^S,mag(jj)) = w(ixO^S,mag(jj)) + qdt * tmp(ixO^S)
            call second_cross_deriv(ixI^L, ixOO^L, nu(ixI^S,ii,jj),  wCT(ixI^S,mag(ii)), jj, ii, tmp)
            w(ixO^S,mag(jj)) = w(ixO^S,mag(jj)) + qdt * tmp(ixO^S)
            !in the total energy
            call second_same_deriv(ixI^L, ixOO^L, nu(ixI^S,jj,ii), wCT(ixI^S,mag(jj)), ii, tmp)
            w(ixO^S,e_c_) = w(ixO^S,e_c_) + qdt * tmp(ixO^S)
            call second_cross_deriv2(ixI^L, ixOO^L, nu(ixI^S,ii,jj), wCT(ixI^S,mag(jj)), wCT(ixI^S,mag(ii)), jj, ii, tmp)
            w(ixO^S,e_c_) = w(ixO^S,e_c_) + qdt * tmp(ixO^S)
          endif

        enddo
      enddo

    end subroutine add_ohmic_hyper_Source   

  end subroutine  add_source_hyperdiffusive

  function dump_hyperdiffusivity_coef_x(ixI^L,ixO^L, w, x, nwc) result(wnew)
    use mod_global_parameters
    use mod_hyperdiffusivity
    integer, intent(in)             :: ixI^L, ixO^L, nwc
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim) 
    double precision   :: wnew(ixO^S, 1:nwc)

    if(nw .ne. nwc) call mpistop("nw != nwc")
    wnew(ixO^S,1:nw) =  dump_hyperdiffusivity_coef_dim(ixI^L,ixO^L, w, x, 1) 

  end function dump_hyperdiffusivity_coef_x

  function dump_hyperdiffusivity_coef_y(ixI^L,ixO^L, w, x, nwc) result(wnew)
    use mod_global_parameters
    use mod_hyperdiffusivity
    integer, intent(in)             :: ixI^L, ixO^L, nwc
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim) 
    double precision   :: wnew(ixO^S, 1:nwc)

    if(nw .ne. nwc) call mpistop("nw != nwc")
    wnew(ixO^S,1:nw) =  dump_hyperdiffusivity_coef_dim(ixI^L,ixO^L, w, x, 2) 

  end function dump_hyperdiffusivity_coef_y

  function dump_hyperdiffusivity_coef_z(ixI^L,ixO^L, w, x, nwc) result(wnew)
    use mod_global_parameters
    use mod_hyperdiffusivity
    integer, intent(in)             :: ixI^L, ixO^L, nwc
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim) 
    double precision   :: wnew(ixO^S, 1:nwc)

    if(nw .ne. nwc) call mpistop("nw != nwc")
    wnew(ixO^S,1:nw) =  dump_hyperdiffusivity_coef_dim(ixI^L,ixO^L, w, x, 3) 

  end function dump_hyperdiffusivity_coef_z

  function dump_hyperdiffusivity_coef_dim(ixI^L,ixOP^L, w, x, ii) result(wnew)
    use mod_global_parameters
    use mod_hyperdiffusivity
    integer, intent(in)             :: ixI^L, ixOP^L, ii
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim) 
    double precision   :: wnew(ixOP^S, 1:nw)

    double precision :: nu(ixI^S),tmp(ixI^S),rho(ixI^S),temp(ixI^S)
    double precision :: divv(ixI^S)
    double precision :: vel(ixI^S,1:ndir)
    double precision :: csound(ixI^S),csound_dim(ixI^S)
    double precision              :: dxarr(ndim)
    integer :: ixOO^L, hxb^L, hx^L,  jj, ixO^L

    ! this is done because of save_physical_boundary = true
    ixOmin^D=max(ixOPmin^D,ixImin^D+3);
    ixOmax^D=min(ixOPmax^D,ixImax^D-3);

    wnew(ixOP^S,1:nw) = 0d0

    ! charges
    call twofl_get_temp_c_pert_from_etot(w, x, ixI^L, ixI^L, temp)
    call twofl_get_v_c(w,x,ixI^L,ixI^L,vel)
    call get_rhoc_tot(w,x,ixI^L,ixI^L,rho)
    call twofl_get_csound2_c_from_conserved(w,x,ixI^L,ixI^L,csound) 
    csound(ixI^S) = sqrt(csound(ixI^S)) + sqrt(twofl_mag_en_all(w,ixI^L,ixI^L) /rho(ixI^S))
    csound(ixI^S) = csound(ixI^S) + sqrt(sum(vel(ixI^S,1:ndir)**2 ,dim=ndim+1))
    !for dim
    call div_vel_coeff(ixI^L, ixOO^L, vel, ii, divv(ixI^S))
    hxmin^D=ixImin^D+1;
    hxmax^D=ixImax^D-1;
    hxb^L=hx^L-kr(ii,^D);
    csound_dim(hx^S) = (csound(hxb^S)+csound(hx^S))/2d0

    !TODO the following is copied
    !rho_c
    call hyp_coeff(ixI^L, ixOO^L, w(ixI^S,rho_c_), ii, tmp(ixI^S))
    nu(ixO^S) = c_hyp(rho_c_) * csound_dim(ixO^S) * dxlevel(ii) *  tmp(ixO^S) + &
                 c_shk(rho_c_) * (dxlevel(ii)**2) *divv(ixO^S)

    wnew(ixO^S,rho_c_) = nu(ixO^S)

    !TH c  
    call hyp_coeff(ixI^L, ixOO^L, temp(ixI^S), ii, tmp(ixI^S))
    nu(ixO^S) = c_hyp(e_c_) * csound_dim(ixO^S) * dxlevel(ii) *  tmp(ixO^S) + &
                 c_shk(e_c_) * (dxlevel(ii)**2) *divv(ixO^S)
    nu(ixO^S) = nu(ixO^S) * rho(ixO^S) * Rc/(twofl_gamma-1d0)
    wnew(ixO^S,e_c_) = nu(ixO^S)

    !visc c
    do jj=1,ndir
      call hyp_coeff(ixI^L, ixOO^L, vel(ixI^S,jj), ii, tmp(ixI^S))
      nu(ixO^S) = c_hyp(mom_c(jj)) * csound_dim(ixO^S) * dxlevel(ii) *  tmp(ixO^S) + &
                 c_shk(mom_c(jj)) * (dxlevel(ii)**2) *divv(ixO^S)
      nu(ixO^S) = nu(ixO^S) * rho(ixO^S) 
      wnew(ixO^S,mom_c(jj)) = nu(ixO^S)
    enddo

    ! Ohmic
    do jj=1,ndir
      if(ii .ne. jj) then
        call hyp_coeff(ixI^L, ixOO^L, w(ixI^S,mag(jj)), ii, tmp(ixI^S))
        nu(ixO^S) = c_hyp(mag(jj)) * csound_dim(ixO^S) * dxlevel(ii) *  tmp(ixO^S) + &
                   c_shk(mag(jj)) * (dxlevel(ii)**2) *divv(ixO^S)
        wnew(ixO^S,mag(jj)) = nu(ixO^S)
      endif
    enddo

    !end for dim
  
    ! neutrals
    call get_rhon_tot(w,x,ixI^L,ixO^L,rho)
    call twofl_get_temp_n_pert_from_etot(w, x, ixI^L, ixI^L, temp)
    call twofl_get_v_n(w,x,ixI^L,ixI^L,vel)
    call  twofl_get_csound_n(w,x,ixI^L,ixI^L,csound)
    csound(ixI^S) = csound(ixI^S) + sqrt(sum(vel(ixI^S,1:ndir)**2 ,dim=ndim+1))
    !for dim
    call div_vel_coeff(ixI^L, ixOO^L, vel, ii, divv(ixI^S))
    hxb^L=ixOO^L-kr(ii,^D);
    csound_dim(ixOO^S) = (csound(hxb^S)+csound(ixOO^S))/2d0
    !rho_n
    call hyp_coeff(ixI^L, ixOO^L, w(ixI^S,rho_n_), ii, tmp(ixI^S))
    nu(ixO^S) = c_hyp(rho_n_) * csound_dim(ixO^S) * dxlevel(ii) *  tmp(ixO^S) + &
                 c_shk(rho_n_) * (dxlevel(ii)**2) *divv(ixOO^S)
    wnew(ixO^S,rho_n_) = nu(ixO^S)

    !TH n  
    call hyp_coeff(ixI^L, ixOO^L, temp(ixI^S), ii, tmp(ixI^S))
    nu(ixO^S) = c_hyp(e_n_) * csound_dim(ixO^S) * dxlevel(ii) *  tmp(ixO^S) + &
                 c_shk(e_n_) * (dxlevel(ii)**2) *divv(ixO^S)
    nu(ixO^S) = nu(ixO^S) * rho(ixO^S) * Rn/(twofl_gamma-1d0)
    wnew(ixO^S,e_n_) = nu(ixO^S)

    !visc n
    do jj=1,ndir
      call hyp_coeff(ixI^L, ixOO^L, vel(ixI^S,jj), ii, tmp(ixI^S))
      nu(ixO^S) = c_hyp(mom_n(jj)) * csound_dim(ixO^S) * dxlevel(ii) *  tmp(ixO^S) + &
                 c_shk(mom_n(jj)) * (dxlevel(ii)**2) *divv(ixO^S)
      nu(ixO^S) = nu(ixO^S) * rho(ixO^S) 
      wnew(ixO^S,mom_n(jj)) = nu(ixO^S)
    enddo
    !end for dim

  end function dump_hyperdiffusivity_coef_dim

  function dump_coll_terms(ixI^L,ixO^L, w, x, nwc) result(wnew)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L, nwc
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim) 
    double precision   :: wnew(ixO^S, 1:nwc)
    double precision   :: tmp(ixI^S),tmp2(ixI^S)
 
    call get_alpha_coll(ixI^L, ixO^L, w, x, tmp(ixI^S))
    wnew(ixO^S,1)= tmp(ixO^S) 
    call get_gamma_ion_rec(ixI^L, ixO^L, w, x, tmp(ixI^S), tmp2(ixI^S))
    wnew(ixO^S,2)= tmp(ixO^S) 
    wnew(ixO^S,3)= tmp2(ixO^S) 

  end function dump_coll_terms

  subroutine get_gamma_ion_rec(ixI^L, ixO^L, w, x, gamma_rec, gamma_ion)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)    ::  x(ixI^S,1:ndim)
    double precision, intent(out)       :: gamma_rec(ixI^S),gamma_ion(ixI^S)
    ! calculations are done in S.I. units
    double precision, parameter :: A = 2.91e-14, & !m3/s
                                K = 0.39, &
                                XX = 0.232, &
                                Eion = 13.6 ! eV      
    double precision, parameter :: ECHARGE=1.6022d-19 !C
    double precision        :: rho(ixI^S), tmp(ixI^S)

    call twofl_get_pthermal_c(w,x,ixI^L,ixO^L,tmp)
    call get_rhoc_tot(w,x,ixI^L,ixO^L,rho)
    tmp(ixO^S) = tmp(ixO^S)/(Rc * rho(ixO^S))

    !transform to SI units
    tmp(ixO^S) = tmp(ixO^S) * unit_temperature * kB_SI/ECHARGE !* BK/ECHARGE means  K to eV 
    !number electrons rho_c = n_e * MH, in normalized units MH=1 and n = rho
    rho(ixO^S) = rho(ixO^S) * unit_numberdensity  
    if(.not. SI_unit) then
      !1/cm^3 = 1e9/m^3
      rho(ixO^S) = rho(ixO^S) * 1d9  
    endif
    gamma_rec(ixO^S) = rho(ixO^S) /sqrt(tmp(ixO^S)) * 2.6e-19
    gamma_ion(ixO^S) = ((rho(ixO^S) * A) /(XX + Eion/tmp(ixO^S))) * ((Eion/tmp(ixO^S))**K) * exp(-Eion/tmp(ixO^S))
    ! see Voronov table: valid for temp min = 1eV(approx 11605 K), Temp max = 20KeV
    !to normalized
    gamma_rec(ixO^S) = gamma_rec(ixO^S) * unit_time  
    gamma_ion(ixO^S) = gamma_ion(ixO^S) * unit_time  

  end subroutine get_gamma_ion_rec

  subroutine get_alpha_coll(ixI^L, ixO^L, w, x, alpha)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)    ::  x(ixI^S,1:ndim)
    double precision, intent(out)       :: alpha(ixI^S)
    if(twofl_alpha_coll_constant) then
      alpha(ixO^S) = twofl_alpha_coll
    else
      call get_alpha_coll_plasma(ixI^L, ixO^L, w, x, alpha)
    endif
  end subroutine get_alpha_coll

  subroutine get_alpha_coll_plasma(ixI^L, ixO^L, w, x, alpha)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)    ::  x(ixI^S,1:ndim)
    double precision, intent(out)       :: alpha(ixI^S)
    double precision        :: pe(ixI^S),rho(ixI^S), tmp(ixI^S), tmp2(ixI^S)

    double precision :: sigma_in = 1e-19 ! m^2
    ! make calculation in SI physical units

    call twofl_get_pthermal_c(w,x,ixI^L,ixO^L,pe)
    call get_rhoc_tot(w,x,ixI^L,ixO^L,rho)
    tmp(ixO^S) = pe(ixO^S)/(Rc * rho(ixO^S))
    call twofl_get_pthermal_n(w,x,ixI^L,ixO^L,pe)
    call get_rhon_tot(w,x,ixI^L,ixO^L,rho)
    tmp2(ixO^S) = pe(ixO^S)/(Rn * rho(ixO^S))
    alpha(ixO^S) = (2d0/(mp_SI**(3d0/2) * sqrt(dpi))*sqrt(0.5*(tmp(ixO^S)+tmp2(ixO^S))*unit_temperature*kB_SI) * sigma_in)*unit_time * unit_density
    if(.not. SI_unit) then
      alpha(ixO^S) = alpha(ixO^S) * 1d3 ! this comes from unit_density: g/cm^3 = 1e-3 kg/m^3
    endif

  end subroutine get_alpha_coll_plasma

  subroutine calc_mult_factor1(ixI^L, ixO^L, step_dt, JJ, res)
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in) :: step_dt
    double precision, intent(in) :: JJ(ixI^S)
    double precision, intent(out) :: res(ixI^S)

    res(ixO^S) = step_dt/(1d0 + step_dt * JJ(ixO^S))

  end subroutine calc_mult_factor1

  subroutine calc_mult_factor2(ixI^L, ixO^L, step_dt, JJ, res)
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in) :: step_dt
    double precision, intent(in) :: JJ(ixI^S)
    double precision, intent(out) :: res(ixI^S)

    res(ixO^S) = (1d0 - exp(-step_dt * JJ(ixO^S)))/JJ(ixO^S)

  end subroutine calc_mult_factor2

  subroutine advance_implicit_grid(ixI^L, ixO^L, w, wout, x, dtfactor,qdt)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in) :: qdt
    double precision, intent(in) :: dtfactor
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)    ::  x(ixI^S,1:ndim)
    double precision, intent(out)       :: wout(ixI^S,1:nw)

    integer :: idir
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),tmp3(ixI^S),tmp4(ixI^S),tmp5(ixI^S)
    double precision :: v_c(ixI^S,ndir), v_n(ixI^S,ndir)
    double precision :: rhon(ixI^S), rhoc(ixI^S), alpha(ixI^S)
    double precision, allocatable :: gamma_rec(:^D&), gamma_ion(:^D&)

    ! copy vars at the indices which are not updated here: mag. field
    wout(ixO^S,mag(:)) = w(ixO^S,mag(:))

    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
    !update density
    if(twofl_coll_inc_ionrec) then
       allocate(gamma_ion(ixI^S), gamma_rec(ixI^S)) 
       call get_gamma_ion_rec(ixI^L, ixO^L, w, x, gamma_rec, gamma_ion)
       tmp2(ixO^S) =  gamma_rec(ixO^S) +  gamma_ion(ixO^S)
       call calc_mult_factor(ixI^L, ixO^L, dtfactor * qdt, tmp2, tmp3) 

      if(.not. twofl_equi_ionrec) then
       tmp(ixO^S) = (-gamma_ion(ixO^S) * rhon(ixO^S) + &
                                        gamma_rec(ixO^S) * rhoc(ixO^S))
      else
       ! equilibrium density does not evolve through ion/rec 
       tmp(ixO^S) = (-gamma_ion(ixO^S) * w(ixO^S,rho_n_) + &
                                        gamma_rec(ixO^S) * w(ixO^S,rho_c_))
      endif
      wout(ixO^S,rho_n_) = w(ixO^S,rho_n_) + tmp(ixO^S) * tmp3(ixO^S)
      wout(ixO^S,rho_c_) = w(ixO^S,rho_c_) - tmp(ixO^S) * tmp3(ixO^S)
    else
      wout(ixO^S,rho_n_) = w(ixO^S,rho_n_)
      wout(ixO^S,rho_c_) = w(ixO^S,rho_c_)
    endif

    call get_alpha_coll(ixI^L, ixO^L, w, x, alpha)
    
    !-J11 + J12    for momentum and kinetic energy
    tmp2(ixO^S) =  alpha(ixO^S) * (rhon(ixO^S) +  rhoc(ixO^S))     
    if(twofl_coll_inc_ionrec) then
      tmp2(ixO^S) = tmp2(ixO^S) + gamma_ion(ixO^S) + gamma_rec(ixO^S)
    endif
    call calc_mult_factor(ixI^L, ixO^L, dtfactor * qdt, tmp2, tmp3) 

    ! momentum update
    do idir=1,ndir

      tmp(ixO^S) = alpha(ixO^S)* (-rhoc(ixO^S) * w(ixO^S,mom_n(idir)) + rhon(ixO^S) * w(ixO^S,mom_c(idir)))
      if(twofl_coll_inc_ionrec) then
        tmp(ixO^S) = tmp(ixO^S) - gamma_ion(ixO^S) * w(ixO^S,mom_n(idir)) + gamma_rec(ixO^S) * w(ixO^S,mom_c(idir))
      endif

      wout(ixO^S,mom_n(idir)) = w(ixO^S,mom_n(idir)) + tmp(ixO^S) * tmp3(ixO^S)
      wout(ixO^S,mom_c(idir)) = w(ixO^S,mom_c(idir)) - tmp(ixO^S) * tmp3(ixO^S)
    enddo

    ! energy update
    
    ! kinetic energy update  
    if(.not. phys_internal_e) then
      ! E_tot includes kinetic energy
      tmp1(ixO^S) = twofl_kin_en_n(w,ixI^L,ixO^L) 
      tmp2(ixO^S) = twofl_kin_en_c(w,ixI^L,ixO^L) 
      tmp4(ixO^S) = w(ixO^S,e_n_) - tmp1(ixO^S) !E_tot - E_kin
      tmp5(ixO^S) = w(ixO^S,e_c_) - tmp2(ixO^S)
      if(phys_total_energy) then
        tmp5(ixO^S) = tmp5(ixO^S) - twofl_mag_en(w,ixI^L,ixO^L)
      endif

      !!implicit update
      tmp(ixO^S) = alpha(ixO^S)*(-rhoc(ixO^S) * tmp1(ixO^S) + rhon(ixO^S) * tmp2(ixO^S))
      if(twofl_coll_inc_ionrec) then
        tmp(ixO^S) = tmp(ixO^S) - gamma_ion(ixO^S) * tmp1(ixO^S) + gamma_rec(ixO^S) * tmp2(ixO^S)
      endif

      wout(ixO^S,e_n_) = w(ixO^S,e_n_) + tmp(ixO^S) * tmp3(ixO^S)
      wout(ixO^S,e_c_) = w(ixO^S,e_c_) - tmp(ixO^S) * tmp3(ixO^S)

     else 
      tmp4(ixO^S) = w(ixO^S,e_n_) 
      tmp5(ixO^S) = w(ixO^S,e_c_) 
      ! calculate velocities, using the already updated variables
      call twofl_get_v_n(wout,x,ixI^L,ixO^L,v_n)
      call twofl_get_v_c(wout,x,ixI^L,ixO^L,v_c)
      tmp1(ixO^S) = alpha(ixO^S) * rhoc(ixO^S) * rhon(ixO^S)
      tmp2(ixO^S) = tmp1(ixO^S) 
      if(twofl_coll_inc_ionrec) then
        tmp1(ixO^S) = tmp1(ixO^S) + rhoc(ixO^S) * gamma_rec(ixO^S)
        tmp2(ixO^S) = tmp2(ixO^S) + rhon(ixO^S) * gamma_ion(ixO^S)
      endif 

      tmp(ixO^S) = 0.5d0 * sum((v_c(ixO^S,1:ndir) - v_n(ixO^S,1:ndir))**2, dim=ndim+1) &
                   * dtfactor * qdt 
      wout(ixO^S,e_n_) = w(ixO^S,e_n_) + tmp(ixO^S)*tmp1(ixO^S)
      wout(ixO^S,e_c_) = w(ixO^S,e_c_) + tmp(ixO^S)*tmp2(ixO^S)
     endif

    !update internal energy
    if(twofl_coll_inc_te) then
     if(.not. twofl_equi_thermal) then   
        if(has_equi_pe_n0) then
          tmp4(ixO^S) = tmp4(ixO^S) + block%equi_vars(ixO^S,equi_pe_n0_,0)*inv_gamma_1  
        endif
        if(has_equi_pe_c0) then
          tmp5(ixO^S) = tmp5(ixO^S) + block%equi_vars(ixO^S,equi_pe_c0_,0)*inv_gamma_1 
        endif
      endif

      tmp(ixO^S) = alpha(ixO^S) *(-rhoc(ixO^S)/Rn * tmp4(ixO^S) + rhon(ixO^S)/Rc * tmp5(ixO^S))
      tmp2(ixO^S) =  alpha(ixO^S) * (rhon(ixO^S)/Rc +  rhoc(ixO^S)/Rn)     
      if(twofl_coll_inc_ionrec) then
        tmp2(ixO^S) =  tmp2(ixO^S) + gamma_rec(ixO^S)/Rc + gamma_ion(ixO^S)/Rn 
        tmp(ixO^S) = tmp(ixO^S) - gamma_ion(ixO^S)/Rn * tmp4(ixO^S) + gamma_rec(ixO^S)/Rc * tmp5(ixO^S)
      endif

      call calc_mult_factor(ixI^L, ixO^L, dtfactor * qdt, tmp2, tmp3) 

      wout(ixO^S,e_n_) = wout(ixO^S,e_n_)+tmp(ixO^S)*tmp3(ixO^S)
      wout(ixO^S,e_c_) = wout(ixO^S,e_c_)-tmp(ixO^S)*tmp3(ixO^S)
    endif
    if(twofl_coll_inc_ionrec) then
       deallocate(gamma_ion, gamma_rec) 
    endif
  end subroutine advance_implicit_grid

  !> Implicit solve of psb=psa+dtfactor*dt*F_im(psb)
  subroutine twofl_implicit_coll_terms_update(dtfactor,qdt,qtC,psb,psa)
    use mod_global_parameters
    use mod_ghostcells_update

    type(state), target :: psa(max_blocks)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    integer :: iigrid, igrid
    !print*, "IMPL call ", it

    call getbc(global_time,0.d0,psa,1,nw)
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      block=>psa(igrid)
      call advance_implicit_grid(ixG^LL, ixG^LL, psa(igrid)%w, psb(igrid)%w, psa(igrid)%x, dtfactor,qdt)
    end do
    !$OMP END PARALLEL DO

  end subroutine twofl_implicit_coll_terms_update 

  !> inplace update of psa==>F_im(psa)
  subroutine twofl_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC

    integer :: iigrid, igrid, level

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      block=>psa(igrid)
       call coll_terms(ixG^LL,ixM^LL,psa(igrid)%w,psa(igrid)%x)
    end do
    !$OMP END PARALLEL DO

  end subroutine twofl_evaluate_implicit

  subroutine coll_terms(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    integer :: idir
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),tmp3(ixI^S),tmp4(ixI^S),tmp5(ixI^S)
    !double precision :: v_c(ixI^S,ndir), v_n(ixI^S,ndir)
    double precision, allocatable :: v_c(:^D&,:), v_n(:^D&,:)
    double precision :: rhon(ixI^S), rhoc(ixI^S), alpha(ixI^S)
    double precision, allocatable :: gamma_rec(:^D&), gamma_ion(:^D&)


    ! get velocity before overwrite density
    call get_rhon_tot(w,x,ixI^L,ixO^L,rhon)
    call get_rhoc_tot(w,x,ixI^L,ixO^L,rhoc)
    if(phys_internal_e) then
      ! get velocity before overwrite momentum
       allocate(v_n(ixI^S,ndir), v_c(ixI^S,ndir)) 
      call twofl_get_v_n(w,x,ixI^L,ixO^L,v_n)
      call twofl_get_v_c(w,x,ixI^L,ixO^L,v_c)
    else
      ! get ke before overwrite density and momentum
      tmp1(ixO^S) = twofl_kin_en_n(w,ixI^L,ixO^L) 
      tmp2(ixO^S) = twofl_kin_en_c(w,ixI^L,ixO^L) 
    endif

    !update density
    if(twofl_coll_inc_ionrec) then
       allocate(gamma_ion(ixI^S), gamma_rec(ixI^S)) 
       call get_gamma_ion_rec(ixI^L, ixO^L, w, x, gamma_rec, gamma_ion)

       if(.not. twofl_equi_ionrec) then
        tmp(ixO^S) = -gamma_ion(ixO^S) * rhon(ixO^S) + &
                                        gamma_rec(ixO^S) * rhoc(ixO^S)
       else
       ! equilibrium density does not evolve through ion/rec 
        tmp(ixO^S) = -gamma_ion(ixO^S) * w(ixO^S,rho_n_) + &
                                        gamma_rec(ixO^S) * w(ixO^S,rho_c_)
       endif
       w(ixO^S,rho_n_) = tmp(ixO^S) 
       w(ixO^S,rho_c_) = -tmp(ixO^S) 
    else
       w(ixO^S,rho_n_) = 0d0 
       w(ixO^S,rho_c_) = 0d0
  
    endif

    call get_alpha_coll(ixI^L, ixO^L, w, x, alpha)

    ! momentum update
    do idir=1,ndir

      tmp(ixO^S) = alpha(ixO^S)* (-rhoc(ixO^S) * w(ixO^S,mom_n(idir)) + rhon(ixO^S) * w(ixO^S,mom_c(idir)))
      if(twofl_coll_inc_ionrec) then
        tmp(ixO^S) = tmp(ixO^S) - gamma_ion(ixO^S) * w(ixO^S,mom_n(idir)) + gamma_rec(ixO^S) * w(ixO^S,mom_c(idir))
      endif

      w(ixO^S,mom_n(idir)) = tmp(ixO^S) 
      w(ixO^S,mom_c(idir)) = -tmp(ixO^S) 
    enddo

    ! energy update
    
    ! kinetic energy update  
    if(.not. phys_internal_e) then
      ! E_tot includes kinetic energy
      tmp4(ixO^S) = w(ixO^S,e_n_) - tmp1(ixO^S) !E_tot - E_kin
      tmp5(ixO^S) = w(ixO^S,e_c_) - tmp2(ixO^S)
      if(phys_total_energy) then
        tmp5(ixO^S) = tmp5(ixO^S) - twofl_mag_en(w,ixI^L,ixO^L)
      endif
      ! tmp4 = eint_n, tmp5 = eint_c
      ! tmp1 = ke_n, tmp2 = ke_c
      tmp(ixO^S) = alpha(ixO^S)*(-rhoc(ixO^S) * tmp1(ixO^S) + rhon(ixO^S) * tmp2(ixO^S))
      if(twofl_coll_inc_ionrec) then
        tmp(ixO^S) = tmp(ixO^S) - gamma_ion(ixO^S) * tmp1(ixO^S) + gamma_rec(ixO^S) * tmp2(ixO^S)
      endif

      w(ixO^S,e_n_) = tmp(ixO^S) 
      w(ixO^S,e_c_) = -tmp(ixO^S) 

     else 
      tmp4(ixO^S) = w(ixO^S,e_n_) 
      tmp5(ixO^S) = w(ixO^S,e_c_) 
      tmp1(ixO^S) = alpha(ixO^S) * rhoc(ixO^S) * rhon(ixO^S)
      tmp2(ixO^S) = tmp1(ixO^S) 
      if(twofl_coll_inc_ionrec) then
        tmp1(ixO^S) = tmp1(ixO^S) + rhoc(ixO^S) * gamma_rec(ixO^S)
        tmp2(ixO^S) = tmp2(ixO^S) + rhon(ixO^S) * gamma_ion(ixO^S)
      endif 

      tmp(ixO^S) = 0.5d0 * sum((v_c(ixO^S,1:ndir) - v_n(ixO^S,1:ndir))**2, dim=ndim+1) 
      w(ixO^S,e_n_) = tmp(ixO^S)*tmp1(ixO^S)
      w(ixO^S,e_c_) = tmp(ixO^S)*tmp2(ixO^S)
     endif

    !update internal energy
    if(twofl_coll_inc_te) then
     if(.not. twofl_equi_thermal) then   
        if(has_equi_pe_n0) then
          tmp4(ixO^S) = tmp4(ixO^S) + block%equi_vars(ixO^S,equi_pe_n0_,0)*inv_gamma_1  
        endif
        if(has_equi_pe_c0) then
          tmp5(ixO^S) = tmp5(ixO^S) + block%equi_vars(ixO^S,equi_pe_c0_,0)*inv_gamma_1 
        endif
      endif

      tmp(ixO^S) = alpha(ixO^S) *(-rhoc(ixO^S)/Rn * tmp4(ixO^S) + rhon(ixO^S)/Rc * tmp5(ixO^S))
      if(twofl_coll_inc_ionrec) then
        tmp(ixO^S) = tmp(ixO^S) - gamma_ion(ixO^S)/Rn * tmp4(ixO^S) + gamma_rec(ixO^S)/Rc * tmp5(ixO^S)
      endif

      w(ixO^S,e_n_) = w(ixO^S,e_n_)+tmp(ixO^S)
      w(ixO^S,e_c_) = w(ixO^S,e_c_)-tmp(ixO^S)
    endif
    if(twofl_coll_inc_ionrec) then
       deallocate(gamma_ion, gamma_rec) 
    endif
    if(phys_internal_e) then
       deallocate(v_n, v_c) 
    endif
    !set contribution to mag field
    w(ixO^S,mag(1:ndir)) = 0d0 

  end subroutine coll_terms

  subroutine twofl_explicit_coll_terms_update(qdt,ixI^L,ixO^L,w,wCT,x)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: wCT(ixI^S,1:nw)

    integer :: idir
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),tmp3(ixI^S),tmp4(ixI^S),tmp5(ixI^S)
    double precision :: v_c(ixI^S,ndir), v_n(ixI^S,ndir)
    double precision :: rhon(ixI^S), rhoc(ixI^S), alpha(ixI^S)
    double precision, allocatable :: gamma_rec(:^D&), gamma_ion(:^D&)

    call get_rhon_tot(wCT,x,ixI^L,ixO^L,rhon)
    call get_rhoc_tot(wCT,x,ixI^L,ixO^L,rhoc)
    !update density
    if(twofl_coll_inc_ionrec) then
       allocate(gamma_ion(ixI^S), gamma_rec(ixI^S)) 
       call get_gamma_ion_rec(ixI^L, ixO^L, wCT, x, gamma_rec, gamma_ion)

      if(.not. twofl_equi_ionrec) then
        tmp(ixO^S) = qdt *(-gamma_ion(ixO^S) * rhon(ixO^S) + &
                                        gamma_rec(ixO^S) * rhoc(ixO^S))
      else
       tmp(ixO^S) = qdt * (-gamma_ion(ixO^S) * wCT(ixO^S,rho_n_) + &
                                        gamma_rec(ixO^S) * wCT(ixO^S,rho_c_))
      endif  
      w(ixO^S,rho_n_) = w(ixO^S,rho_n_) + tmp(ixO^S) 
      w(ixO^S,rho_c_) = w(ixO^S,rho_c_) - tmp(ixO^S) 
    endif

    call get_alpha_coll(ixI^L, ixO^L, wCT, x, alpha)

    ! momentum update
    do idir=1,ndir

      tmp(ixO^S) = alpha(ixO^S)* (-rhoc(ixO^S) * wCT(ixO^S,mom_n(idir)) + rhon(ixO^S) * wCT(ixO^S,mom_c(idir)))
      if(twofl_coll_inc_ionrec) then
        tmp(ixO^S) = tmp(ixO^S) - gamma_ion(ixO^S) * wCT(ixO^S,mom_n(idir)) + gamma_rec(ixO^S) * wCT(ixO^S,mom_c(idir))
      endif
      tmp(ixO^S) =tmp(ixO^S) * qdt

      w(ixO^S,mom_n(idir)) = w(ixO^S,mom_n(idir)) + tmp(ixO^S) 
      w(ixO^S,mom_c(idir)) = w(ixO^S,mom_c(idir)) - tmp(ixO^S) 
    enddo

    ! energy update
    
    ! kinetic energy update  
    if(.not. phys_internal_e) then
      ! E_tot includes kinetic energy
      tmp1(ixO^S) = twofl_kin_en_n(wCT,ixI^L,ixO^L) 
      tmp2(ixO^S) = twofl_kin_en_c(wCT,ixI^L,ixO^L) 
      tmp4(ixO^S) = wCT(ixO^S,e_n_) - tmp1(ixO^S) !E_tot - E_kin
      tmp5(ixO^S) = wCT(ixO^S,e_c_) - tmp2(ixO^S)
      if(phys_total_energy) then
        tmp5(ixO^S) = tmp5(ixO^S) - twofl_mag_en(wCT,ixI^L,ixO^L)
      endif

      tmp(ixO^S) = alpha(ixO^S)*(-rhoc(ixO^S) * tmp1(ixO^S) + rhon(ixO^S) * tmp2(ixO^S))
      if(twofl_coll_inc_ionrec) then
        tmp(ixO^S) = tmp(ixO^S) - gamma_ion(ixO^S) * tmp1(ixO^S) + gamma_rec(ixO^S) * tmp2(ixO^S)
      endif
      tmp(ixO^S) =tmp(ixO^S) * qdt

      w(ixO^S,e_n_) = w(ixO^S,e_n_) + tmp(ixO^S) 
      w(ixO^S,e_c_) = w(ixO^S,e_c_) - tmp(ixO^S) 

    else 
      tmp4(ixO^S) = w(ixO^S,e_n_) 
      tmp5(ixO^S) = w(ixO^S,e_c_) 
      call twofl_get_v_n(wCT,x,ixI^L,ixO^L,v_n)
      call twofl_get_v_c(wCT,x,ixI^L,ixO^L,v_c)
      tmp1(ixO^S) = alpha(ixO^S) * rhoc(ixO^S) * rhon(ixO^S)
      tmp2(ixO^S) = tmp1(ixO^S) 
      if(twofl_coll_inc_ionrec) then
        tmp1(ixO^S) = tmp1(ixO^S) + rhoc(ixO^S) * gamma_rec(ixO^S)
        tmp2(ixO^S) = tmp2(ixO^S) + rhon(ixO^S) * gamma_ion(ixO^S)
      endif 

      tmp(ixO^S) = 0.5d0 * sum((v_c(ixO^S,1:ndir) - v_n(ixO^S,1:ndir))**2, dim=ndim+1) * qdt 
      w(ixO^S,e_n_) = w(ixO^S,e_n_) + tmp(ixO^S)*tmp1(ixO^S)
      w(ixO^S,e_c_) = w(ixO^S,e_c_) + tmp(ixO^S)*tmp2(ixO^S)
     endif

    !update internal energy
    if(twofl_coll_inc_te) then
     if(.not. twofl_equi_thermal) then   
        if(has_equi_pe_n0) then
          tmp4(ixO^S) = tmp4(ixO^S) + block%equi_vars(ixO^S,equi_pe_n0_,0)*inv_gamma_1  
        endif
        if(has_equi_pe_c0) then
          tmp5(ixO^S) = tmp5(ixO^S) + block%equi_vars(ixO^S,equi_pe_c0_,0)*inv_gamma_1 
        endif
      endif

      tmp(ixO^S) = alpha(ixO^S) *(-rhoc(ixO^S)/Rn * tmp4(ixO^S) + rhon(ixO^S)/Rc * tmp5(ixO^S))
      if(twofl_coll_inc_ionrec) then
        tmp(ixO^S) = tmp(ixO^S) - gamma_ion(ixO^S)/Rn * tmp4(ixO^S) + gamma_rec(ixO^S)/Rc * tmp5(ixO^S)
      endif

      tmp(ixO^S) =tmp(ixO^S) * qdt

      w(ixO^S,e_n_) = w(ixO^S,e_n_)+tmp(ixO^S)
      w(ixO^S,e_c_) = w(ixO^S,e_c_)-tmp(ixO^S)
    endif
    if(twofl_coll_inc_ionrec) then
       deallocate(gamma_ion, gamma_rec) 
    endif
  end subroutine twofl_explicit_coll_terms_update

end module mod_twofl_phys
