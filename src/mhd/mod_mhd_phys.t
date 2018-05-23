!> Magneto-hydrodynamics module
module mod_mhd_phys
  use mod_global_parameters, only: std_len
  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: mhd_energy = .true.

  !> Whether thermal conduction is used
  logical, public, protected              :: mhd_thermal_conduction = .false.

  !> Whether radiative cooling is added
  logical, public, protected              :: mhd_radiative_cooling = .false.

  !> Whether viscosity is added
  logical, public, protected              :: mhd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: mhd_gravity = .false.

  !> Whether Hall-MHD is used
  logical, public, protected              :: mhd_Hall = .false.

  !> Whether particles module is added
  logical, public, protected              :: mhd_particles = .false.

  !> Whether magnetofriction is added
  logical, public, protected              :: mhd_magnetofriction = .false.

  !> Whether GLM-MHD is used
  logical, public, protected              :: mhd_glm = .false.

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: mhd_glm_alpha = 0.5d0

  !> MHD fourth order
  logical, public, protected              :: mhd_4th_order = .false.

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

  !> The small_est allowed energy
  double precision, protected             :: small_e

  !> The number of waves
  integer :: nwwave=8

  !> Method type to clean divergence of B
  character(len=std_len) :: typedivbfix  = 'linde'

  !> Whether divB is computed with a fourth order approximation
  logical, public :: mhd_divb_4thorder = .false.

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
  double precision, public, protected  :: He_abundance=0.1d0

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*^ND)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public, protected :: boundary_divbfix_skip(2*^ND)=0

  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.

  !> gamma minus one and its inverse
  double precision :: gamma_1, inv_gamma_1

  ! DivB cleaning methods
  integer, parameter :: divb_none          = 0
  integer, parameter :: divb_multigrid     = -1
  integer, parameter :: divb_glm1          = 1
  integer, parameter :: divb_glm2          = 2
  integer, parameter :: divb_glm3          = 3
  integer, parameter :: divb_powel         = 4
  integer, parameter :: divb_janhunen      = 5
  integer, parameter :: divb_linde         = 6
  integer, parameter :: divb_lindejanhunen = 7
  integer, parameter :: divb_lindepowel    = 8
  integer, parameter :: divb_lindeglm      = 9


  ! Public methods
  public :: mhd_phys_init
  public :: mhd_kin_en
  public :: mhd_get_pthermal
  public :: mhd_get_v
  public :: mhd_to_conserved
  public :: mhd_to_primitive
  public :: mhd_get_csound2
  public :: get_divb
  public :: get_current
  public :: get_normalized_divb

contains

  !> Read this module"s parameters from a file
  subroutine mhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mhd_list/ mhd_energy, mhd_n_tracer, mhd_gamma, mhd_adiab,&
      mhd_eta, mhd_eta_hyper, mhd_etah, mhd_glm_alpha, mhd_magnetofriction,&
      mhd_thermal_conduction, mhd_radiative_cooling, mhd_Hall, mhd_gravity,&
      mhd_viscosity, mhd_4th_order, typedivbfix, source_split_divb, divbdiff,&
      typedivbdiff, compactres, divbwave, He_abundance, SI_unit, B0field,&
      B0field_forcefree, Bdip, Bquad, Boct, Busr, mhd_particles,&
      boundary_divbfix, boundary_divbfix_skip

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
    ! ! shifted indexes
    ! hxO^L=ixO^L-kr(idim,^D);
    ! ! all the indexes
    ! kxCmin^D=hxOmin^D;
    ! kxCmax^D=ixOmax^D;
    !
    ! inv_volume = 1.0d0/block%dvolume(ixO^S)
    !
    ! select case(typeaxial)
    ! case ("cylindrical")
    !   do iw=1,nwflux
    !     if (idim==r_ .and. iw==iw_mom(phi_)) then
    !       fC(kxC^S,iw,idim)= fC(kxC^S,iw,idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !            (inv_volume/x(ixO^S,r_))
    !     else
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !             inv_volume
    !     endif
    !   enddo
    ! case ("spherical")
    !   do iw=1,nwflux
    !     if     (idim==r_ .and. (iw==iw_mom(2) .or. iw==iw_mom(phi_))) then
    !       fC(kxC^S,iw,idim)= fC(kxC^S,iw,idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !            (inv_volume/x(ixO^S,r_))
    !     elseif (idim==2  .and. iw==iw_mom(phi_)) then
    !       fC(kxC^S,iw,idim)=fC(kxC^S,iw,idim)*dsin(x(kxC^S,2)+half*block%dx(kxC^S,2)) ! (x(4,3,1)-x(3,3,1)))
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !            (inv_volume/dsin(x(ixO^S,2)))
    !     else
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !             inv_volume
    !     endif
    !   enddo
    !
    !   ! if (idim==r_) then
    !   !   fC(kxC^S,iw_mom(phi_),idim)= fC(kxC^S,iw_mom(phi_),idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !   !   fC(kxC^S,iw_mom(phi_),idim)= fC(kxC^S,iw_mom(phi_),idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !   !   wnew(ixO^S,iw_mom(phi_))=wnew(ixO^S,iw_mom(phi_)) + (fC(ixO^S,iw_mom(phi_),idim)-fC(hxO^S,iw_mom(phi_),idim)) * &
    !   !        (inv_volume/x(ixO^S,r_))
    !   !
    !   ! elseif (idim==2) then
    !   !   fC(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,iw,idim)=fC(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,iw,idim)*dsin(x(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,2)+half*block%dx(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,2)) ! (x(4,3,1)-x(3,3,1)))
    !   !   fC(ixO^S,iw,idim)=fC(ixO^S,iw,idim)*dsin(x(ixO^S,2)+half*block%dx(ixO^S,2)) ! (x(4,3,1)-x(3,3,1)))
    !   !   wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !   !        (inv_volume/dsin(x(ixO^S,2)))
    !   ! endif
    !
    ! end select

  end subroutine mhd_angmomfix

  subroutine mhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_magnetofriction, only: magnetofriction_init
    use mod_physics
    {^NOONED
    use mod_multigrid_coupling
    }

    integer :: itr, idir

    call mhd_read_params(par_files)

    physics_type = "mhd"
    phys_energy=mhd_energy
    ! set default gamma for polytropic/isothermal process
    if(.not.mhd_energy) mhd_gamma=1.d0
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
       mhd_divb_4thorder = .true.
       phys_global_source => mhd_clean_divb_multigrid
    }
    case ('glm1')
       mhd_glm          = .true.
       need_global_cmax = .true.
       type_divb        = divb_glm1
    case ('glm2')
       mhd_glm          = .true.
       need_global_cmax = .true.
       type_divb        = divb_glm2
    case ('glm3')
       mhd_glm          = .true.
       need_global_cmax = .true.
       type_divb        = divb_glm3
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
    case default
      call mpistop('Unknown divB fix')
    end select

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

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?
    iw_vector(2) = mag(1) - 1   ! TODO: why like this?

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if
    do idir=1,ndir
      if(ndim>1) flux_type(idir,mag(idir))=flux_tvdlf
    end do
    if(mhd_glm .and. ndim>1) flux_type(:,psi_)=flux_tvdlf

    phys_get_dt              => mhd_get_dt
    phys_get_cmax            => mhd_get_cmax
    phys_get_cbounds         => mhd_get_cbounds
    phys_get_flux            => mhd_get_flux
    phys_add_source_geom     => mhd_add_source_geom
    phys_add_source          => mhd_add_source
    phys_to_conserved        => mhd_to_conserved
    phys_to_primitive        => mhd_to_primitive
    phys_check_params        => mhd_check_params
    phys_check_w             => mhd_check_w
    phys_get_pthermal        => mhd_get_pthermal
    phys_boundary_adjust     => mhd_boundary_adjust
    phys_write_info          => mhd_write_info
    phys_angmomfix           => mhd_angmomfix
    phys_handle_small_values => mhd_handle_small_values

    ! Whether diagonal ghost cells are required for the physics
    if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! derive units from basic units
    call mhd_physical_units()

    if(.not. mhd_energy .and. mhd_thermal_conduction) then
      call mpistop("thermal conduction needs mhd_energy=T")
    end if
    if(.not. mhd_energy .and. mhd_radiative_cooling) then
      call mpistop("radiative cooling needs mhd_energy=T")
    end if

    ! initialize thermal conduction module
    if (mhd_thermal_conduction) then
      phys_req_diagonal = .true.
      call thermal_conduction_init(mhd_gamma)
    end if

    ! Initialize radiative cooling module
    if (mhd_radiative_cooling) then
      call radiative_cooling_init(mhd_gamma,He_abundance)
    end if

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
    end if

    ! initialize magnetofriction module
    if (mhd_magnetofriction) then
      phys_req_diagonal = .true.
      call magnetofriction_init()
    end if

    if (mhd_glm) then
       ! Solve the Riemann problem for the linear 2x2 system for normal
       ! B-field and GLM_Psi according to Dedner 2002:
       phys_modify_wLR => glmSolve
    end if

    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in getflux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in nghostcells.
    if (mhd_hall) then
       if (mhd_4th_order) then
          phys_wider_stencil = 2
       else
          phys_wider_stencil = 1
       end if
    end if

  end subroutine mhd_phys_init

  subroutine mhd_check_params
    use mod_global_parameters

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

  end subroutine mhd_check_params

  subroutine mhd_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi
    end if
    if(unit_velocity==0) then
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*kB*unit_temperature
      unit_velocity=sqrt(unit_pressure/unit_density)
      unit_magneticfield=sqrt(miu0*unit_pressure)
      unit_time=unit_length/unit_velocity
    else
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*kB)
      unit_magneticfield=sqrt(miu0*unit_pressure)
      unit_time=unit_length/unit_velocity
    end if

  end subroutine mhd_physical_units

  subroutine mhd_check_w(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    integer, intent(inout) :: flag(ixI^S)
    double precision :: tmp(ixI^S)

    flag(ixO^S)=0
    where(w(ixO^S, rho_) < small_density) flag(ixO^S) = rho_

    if (mhd_energy) then
       if (block%e_is_internal) then
          where(w(ixO^S, e_) < small_pressure*inv_gamma_1) flag(ixO^S) = e_
       else
         if (primitive)then
           where(w(ixO^S, e_) < small_pressure) flag(ixO^S) = e_
         else
        ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
           tmp(ixO^S)=gamma_1*(w(ixO^S,e_) - &
              mhd_kin_en(w,ixI^L,ixO^L)-mhd_mag_en(w,ixI^L,ixO^L))
           where(tmp(ixO^S) < small_pressure) flag(ixO^S) = e_
         end if
       end if
    end if
  end subroutine mhd_check_w

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: idir, itr

    if (mhd_energy) then
       ! Calculate total energy from pressure, kinetic and magnetic energy
       w(ixO^S,e_)=w(ixO^S,p_)*inv_gamma_1
      if(.not.block%e_is_internal) w(ixO^S,e_)=w(ixO^S,e_) + &
        0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) * w(ixO^S, rho_) + &
        mhd_mag_en(w, ixI^L, ixO^L)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, rho_) * w(ixO^S, mom(idir))
    end do

    if (check_small_values) call mhd_handle_small_values(.false., w, x, ixI^L, ixO^L,'mhd_to_conserved')
  end subroutine mhd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision                :: inv_rho(ixO^S)
    integer                         :: itr, idir

    inv_rho = 1.0d0 / w(ixO^S, rho_)

    if (mhd_energy) then
      if(block%e_is_internal) then
        w(ixO^S, p_) = gamma_1*w(ixO^S, e_)
      else
       ! Calculate pressure = (gamma-1) * (e-0.5*(2ek+2eb))
        w(ixO^S, p_) = gamma_1*(w(ixO^S, e_) &
              - mhd_kin_en(w, ixI^L, ixO^L, inv_rho) &
              - mhd_mag_en(w, ixI^L, ixO^L))
      end if
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))*inv_rho
    end do

    if (check_small_values) call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L,'mhd_to_primitive')
  end subroutine mhd_to_primitive

  subroutine mhd_handle_small_values(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    double precision :: smallone
    integer :: idir, flag(ixI^S)

    if (small_values_method == "ignore") return

    call mhd_check_w(primitive, ixI^L, ixO^L, w, flag)

    if (any(flag(ixO^S) /= 0)) then
       select case (small_values_method)
       case ("replace")
          where(flag(ixO^S) /= 0) w(ixO^S,rho_) = small_density

          do idir = 1, ndir
             where(flag(ixO^S) /= 0) w(ixO^S, mom(idir)) = 0.0d0
          end do

          if (mhd_energy) then
             if(primitive) then
               smallone = small_pressure
             else
               smallone = small_e
             end if
             where(flag(ixO^S) /= 0) w(ixO^S,e_) = smallone
          end if
       case ("average")
          call small_values_average(ixI^L, ixO^L, w, x, flag)
       case default
          call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
       end select
    end if
  end subroutine mhd_handle_small_values

  !> Convert energy to entropy
  subroutine e_to_rhos(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision,intent(inout)  :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    if (mhd_energy) then
      if(.not.block%e_is_internal) &
        w(ixO^S, e_) = w(ixO^S, e_) - mhd_kin_en(w, ixI^L, ixO^L) &
              - mhd_mag_en(w, ixI^L, ixO^L)
      w(ixO^S, e_) = gamma_1* w(ixO^S, rho_)**(1.0d0 - mhd_gamma) * &
            w(ixO^S, e_)
    else
      call mpistop("e_to_rhos can not be used without energy equation!")
    end if
  end subroutine e_to_rhos

  !> Convert entropy to energy
  subroutine rhos_to_e(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    if (mhd_energy) then
       w(ixO^S, e_) = w(ixO^S, rho_)**gamma_1 * w(ixO^S, e_) &
            * inv_gamma_1
       if(.not.block%e_is_internal) &
         w(ixO^S, e_) =w(ixO^S, e_) + mhd_kin_en(w, ixI^L, ixO^L) + &
            mhd_mag_en(w, ixI^L, ixO^L)
    else
       call mpistop("rhos_to_e can not be used without energy equation!")
    end if
  end subroutine rhos_to_e

  !> Calculate v vector
  subroutine mhd_get_v(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,ndir)

    integer :: idir

    do idir=1,ndir
      v(ixO^S,idir) = w(ixO^S, mom(idir)) / w(ixO^S,rho_)
    end do

  end subroutine mhd_get_v

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine mhd_get_cmax(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)

    call mhd_get_csound(w,x,ixI^L,ixO^L,idim,cmax)

    cmax(ixO^S)=abs(w(ixO^S,mom(idim))/w(ixO^S,rho_))+cmax(ixO^S)

  end subroutine mhd_get_cmax

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine mhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision :: wmean(ixI^S,nw)
    double precision, dimension(ixI^S) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3

    if (typeboundspeed/='cmaxmean') then
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.

      tmp1(ixO^S)=sqrt(wLp(ixO^S,rho_))
      tmp2(ixO^S)=sqrt(wRp(ixO^S,rho_))
      tmp3(ixO^S)=1.d0/(sqrt(wLp(ixO^S,rho_))+sqrt(wRp(ixO^S,rho_)))
      umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)+wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)
      call mhd_get_csound_prim(wLp,x,ixI^L,ixO^L,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixI^L,ixO^L,idim,csoundR)
      dmean(ixO^S)=(tmp1(ixO^S)*csoundL(ixO^S)**2+tmp2(ixO^S)*csoundR(ixO^S)**2)*tmp3(ixO^S)+&
       0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2*&
       (wRp(ixO^S,mom(idim))-wLp(ixO^S,mom(idim)))**2
      dmean(ixO^S)=sqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S)=umean(ixO^S)+dmean(ixO^S)
      else
        cmax(ixO^S)=abs(umean(ixO^S))+dmean(ixO^S)
      end if
    else
      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
      tmp1(ixO^S)=wmean(ixO^S,mom(idim))/wmean(ixO^S,rho_)
      call mhd_get_csound(wmean,x,ixI^L,ixO^L,idim,csoundR)
      if(present(cmin)) then
        cmax(ixO^S)=max(tmp1(ixO^S)+csoundR(ixO^S),zero)
        cmin(ixO^S)=min(tmp1(ixO^S)-csoundR(ixO^S),zero)
      else
        cmax(ixO^S)=abs(tmp1(ixO^S))+csoundR(ixO^S)
      end if
    end if

  end subroutine mhd_get_cbounds

  !> Calculate fast magnetosonic wave speed
  subroutine mhd_get_csound(w,x,ixI^L,ixO^L,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
    double precision :: inv_rho(ixO^S)

    inv_rho=1.d0/w(ixO^S,rho_)

    call mhd_get_csound2(w,x,ixI^L,ixO^L,csound)
    ! store |B|^2 in v
    b2(ixO^S)        = mhd_mag_en_all(w,ixI^L,ixO^L)
    cfast2(ixO^S)   = b2(ixO^S) * inv_rho+csound(ixO^S)
    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
         * mhd_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
         * inv_rho

    where(AvMinCs2(ixO^S)<zero)
       AvMinCs2(ixO^S)=zero
    end where

    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

    if (.not. MHD_Hall) then
       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
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

    inv_rho=1.d0/w(ixO^S,rho_)

    if(mhd_energy) then
      csound(ixO^S)=mhd_gamma*w(ixO^S,p_)*inv_rho
    else
      csound(ixO^S)=mhd_gamma*mhd_adiab*w(ixO^S,rho_)**gamma_1
    end if
    ! store |B|^2 in v
    b2(ixO^S)        = mhd_mag_en_all(w,ixI^L,ixO^L)
    cfast2(ixO^S)   = b2(ixO^S) * inv_rho+csound(ixO^S)
    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
         * mhd_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
         * inv_rho

    where(AvMinCs2(ixO^S)<zero)
       AvMinCs2(ixO^S)=zero
    end where

    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

    if (.not. MHD_Hall) then
       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
            mhd_etah * sqrt(b2(ixO^S))*inv_rho*kmax)
    end if

  end subroutine mhd_get_csound_prim

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)

    if(mhd_energy) then
      if(block%e_is_internal) then
        pth(ixO^S)=gamma_1*w(ixO^S,e_)
      else
        pth(ixO^S)=gamma_1*(w(ixO^S,e_)&
           - mhd_kin_en(w,ixI^L,ixO^L)&
           - mhd_mag_en(w,ixI^L,ixO^L))
      end if
    else
      pth(ixO^S)=mhd_adiab*w(ixO^S,rho_)**mhd_gamma
    end if
  end subroutine mhd_get_pthermal

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine mhd_get_csound2(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)

    if(mhd_energy) then
      call mhd_get_pthermal(w,x,ixI^L,ixO^L,csound2)
      csound2(ixO^S)=mhd_gamma*csound2(ixO^S)/w(ixO^S,rho_)
    else
      csound2(ixO^S)=mhd_gamma*mhd_adiab*w(ixO^S,rho_)**gamma_1
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

  !> Calculate fluxes within ixO^L.
  subroutine mhd_get_flux(wC,w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)

    double precision             :: ptotal(ixO^S),tmp(ixI^S)
    double precision, allocatable:: vHall(:^D&,:)
    integer                      :: idirmin, iw, idir

    if (mhd_Hall) then
      allocate(vHall(ixI^S,1:ndir))
      call mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
    end if

    if(B0field) tmp(ixO^S)=sum(block%B0(ixO^S,:,idim)*w(ixO^S,mag(:)),dim=ndim+1)

    if(mhd_energy) then
      ptotal=w(ixO^S,p_)+0.5d0*sum(w(ixO^S, mag(:))**2, dim=ndim+1)
    else
      ptotal(ixO^S)=mhd_adiab*w(ixO^S,rho_)**mhd_gamma+0.5d0*sum(w(ixO^S, mag(:))**2, dim=ndim+1)
    end if

    ! Get flux of density
    f(ixO^S,rho_)=w(ixO^S,mom(idim))*w(ixO^S,rho_)

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixO^S,tracer(iw))=w(ixO^S,mom(idim))*w(ixO^S,tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixO^S,mom(idir))=ptotal(ixO^S)-w(ixO^S,mag(idim))*w(ixO^S,mag(idir))
        if(B0field) f(ixO^S,mom(idir))=f(ixO^S,mom(idir))+tmp(ixO^S)
      else
        f(ixO^S,mom(idir))= -w(ixO^S,mag(idir))*w(ixO^S,mag(idim))
      end if
      if (B0field) then
        f(ixO^S,mom(idir))=f(ixO^S,mom(idir))&
             -w(ixO^S,mag(idir))*block%B0(ixO^S,idim,idim)&
             -w(ixO^S,mag(idim))*block%B0(ixO^S,idir,idim)
      end if
      f(ixO^S,mom(idir))=f(ixO^S,mom(idir))+w(ixO^S,mom(idim))*wC(ixO^S,mom(idir))
    end do

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if (mhd_energy) then
       if (block%e_is_internal) then
          f(ixO^S,e_)=w(ixO^S,mom(idim))*w(ixO^S,p_)
          if (mhd_Hall) then
             call mpistop("solve pthermal not designed for Hall MHD")
          endif
       else
          f(ixO^S,e_)=w(ixO^S,mom(idim))*(wC(ixO^S,e_) + ptotal(ixO^S))- &
             w(ixO^S,mag(idim))*sum(w(ixO^S,mag(:))*w(ixO^S,mom(:)),dim=ndim+1)

          if (B0field) then
             f(ixO^S,e_) = f(ixO^S,e_) &
                + w(ixO^S,mom(idim)) * tmp(ixO^S) &
                - sum(w(ixO^S,mom(:))*w(ixO^S,mag(:)),dim=ndim+1) * block%B0(ixO^S,idim,idim)
          end if

          if (mhd_Hall) then
          ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
             if (mhd_etah>zero) then
                f(ixO^S,e_) = f(ixO^S,e_) + vHall(ixO^S,idim) * &
                   sum(w(ixO^S, mag(:))**2,dim=ndim+1) &
                   - w(ixO^S,mag(idim)) * sum(vHall(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1)
                if (B0field) then
                   f(ixO^S,e_) = f(ixO^S,e_) &
                      + vHall(ixO^S,idim) * tmp(ixO^S) &
                      - sum(vHall(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1) * block%B0(ixO^S,idim,idim)
                end if
             end if
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

        if (B0field) then
          f(ixO^S,mag(idir))=f(ixO^S,mag(idir))&
                +w(ixO^S,mom(idim))*block%B0(ixO^S,idir,idim)&
                -w(ixO^S,mom(idir))*block%B0(ixO^S,idim,idim)
        end if

        if (mhd_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (mhd_etah>zero) then
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

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixO^S,psi_)  = cmax_global**2*w(ixO^S,mag(idim))
    end if

  end subroutine mhd_get_flux

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine mhd_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
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

    if (.not. qsourcesplit) then
      ! Source for solving internal energy
      if (mhd_energy .and. block%e_is_internal) then
        active = .true.
        call internal_energy_add_source(qdt,ixI^L,ixO^L,wCT,w,x)
      endif

      ! Source for B0 splitting
      if (B0field) then
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
    end if

      {^NOONED
    if(.not.source_split_divb .and. .not.qsourcesplit .and. istep==nstep) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
      case (divb_glm1)
        active = .true.
        call add_source_glm1(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
      case (divb_glm2)
        active = .true.
        call add_source_glm2(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
      case (divb_glm3)
        active = .true.
        call add_source_glm3(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
      case (divb_linde)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
        call add_source_janhunen(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
        call add_source_powel(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
        call add_source_glm3(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
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
      case (divb_glm1)
        active = .true.
        call add_source_glm1(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_glm2)
        active = .true.
        call add_source_glm2(qdt,ixI^L,ixO^L,wCT,w,x)
      case (divb_glm3)
        active = .true.
        call add_source_glm3(qdt,ixI^L,ixO^L,wCT,w,x)
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
        call add_source_glm3(qdt,ixI^L,ixO^L,wCT,w,x)
     case (divb_multigrid)
        continue ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    end if
    }

    if(mhd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,qsourcesplit,active)
    end if

    if(mhd_viscosity) then
      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,mhd_energy,qsourcesplit,active)
    end if

    if(mhd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,mhd_energy,qsourcesplit,active)
    end if

  end subroutine mhd_add_source

  subroutine internal_energy_add_source(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: pth(ixI^S),v(ixI^S,1:ndir),divv(ixI^S)

    call mhd_get_v(wCT,x,ixI^L,ixI^L,v)
    call divvector(v,ixI^L,ixO^L,divv)
    call mhd_get_pthermal(wCT,x,ixI^L,ixO^L,pth)
    w(ixO^S,e_)=w(ixO^S,e_)-qdt*pth(ixO^S)*divv(ixO^S)

  end subroutine internal_energy_add_source

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

    if(mhd_energy) then
      if(.not.block%e_is_internal) then
        a=0.d0
        ! for free-free field -(vxB0) dot J0 =0
        b(ixO^S,:)=wCT(ixO^S,mag(:))
        ! store full magnetic field B0+B1 in b
        if(.not.B0field_forcefree) b(ixO^S,:)=b(ixO^S,:)+block%B0(ixO^S,:,0)
        ! store velocity in a
        do idir=1,ndir
          a(ixO^S,idir)=wCT(ixO^S,mom(idir))/wCT(ixO^S,rho_)
        end do
        call cross_product(ixI^L,ixO^L,a,b,axb)
        axb(ixO^S,:)=axb(ixO^S,:)*qdt
        ! add -(vxB) dot J0 source term in energy equation
        do idir=7-2*ndir,3
          w(ixO^S,e_)=w(ixO^S,e_)-axb(ixO^S,idir)*block%J0(ixO^S,idir)
        end do
      end if
    end if

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_B0')

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

    if (mhd_energy) then
       ! de/dt+=eta*J**2
       tmp(ixO^S)=zero
       do idir=idirmin,3
          tmp(ixO^S)=tmp(ixO^S)+current(ixO^S,idir)**2
       end do
       w(ixO^S,e_)=w(ixO^S,e_)+qdt*eta(ixO^S)*tmp(ixO^S)
    end if

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res1')

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
    integer :: ixA^L,idir,jdir,kdir,idirmin,iw,idim,idirmin1

    double precision :: tmp(ixI^S),tmp2(ixI^S)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S),curlj(ixI^S,1:3)
    double precision :: tmpvec(ixI^S,1:3)

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
    do jdir=idirmin,3
       tmpvec(ixA^S,jdir)=current(ixA^S,jdir)*eta(ixA^S)
    end do
    call curlvector(tmpvec,ixI^L,ixO^L,curlj,idirmin1,1,3)
    do idir=1,ndir
      w(ixO^S,mag(idir)) = w(ixO^S,mag(idir))-qdt*curlj(ixO^S,idir)
    end do

    if(mhd_energy) then
      ! de/dt= +div(B x Jeta) = eta J^2 - B dot curl(eta J)
      ! de1/dt= eta J^2 - B1 dot curl(eta J)
      w(ixO^S,e_)=w(ixO^S,e_)+qdt*(sum(current(ixO^S,:)**2,dim=ndim+1)*eta(ixO^S)-&
        sum(wCT(ixO^S,mag(1:ndir))*curlj(ixO^S,1:ndir),dim=ndim+1))
    end if

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res2')

  end subroutine add_source_res2

  !> Add Hyper-resistive source to w within ixO
  !> Uses 9 point stencil (4 neighbours) in each direction.
  subroutine add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry
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

    if (check_small_values)  call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_hyperres')

  end subroutine add_source_hyperres

  subroutine add_source_glm1(qdt,ixI^L,ixO^L,wCT,w,x)
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
    call get_divb(wCT,ixI^L,ixO^L,divb, mhd_divb_4thorder)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (mhd_glm_alpha < zero) then
      w(ixO^S,psi_) = abs(mhd_glm_alpha)*wCT(ixO^S,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab) then
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(dxlevel(:)))*w(ixO^S,psi_)
      else
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(block%ds(ixO^S,:),dim=ndim+1))*w(ixO^S,psi_)
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
       if (mhd_energy .and. .not.block%e_is_internal) then
       ! e  = e  -qdt (b . grad(Psi))
         w(ixO^S,e_) = w(ixO^S,e_)-qdt*wCT(ixO^S,mag(idim))*gradPsi(ixO^S)
       end if
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-qdt*mhd_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_glm1')

  end subroutine add_source_glm1

  subroutine add_source_glm2(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 38_
    ! giving the non conservative EGLM-MHD scheme.
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision:: divb(ixI^S),v(ixI^S,1:ndir)
    integer          :: idim,idir
    double precision :: gradPsi(ixI^S)

    ! calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb, mhd_divb_4thorder)

    ! calculate velocity
    call mhd_get_v(wCT,x,ixI^L,ixO^L,v)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (mhd_glm_alpha < zero) then
      w(ixO^S,psi_) = abs(mhd_glm_alpha)*wCT(ixO^S,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab) then
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(dxlevel(:)))*w(ixO^S,psi_)
      else
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(block%ds(ixO^S,:),dim=ndim+1))*w(ixO^S,psi_)
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

      ! Psi=Psi - qdt (v . grad(Psi))
      w(ixO^S,psi_) = w(ixO^S,psi_)-qdt*v(ixO^S,idim)*gradPsi(ixO^S)

      if (mhd_energy .and. .not.block%e_is_internal) then
      ! e  = e  - qdt (b . grad(Psi))
        w(ixO^S,e_) = w(ixO^S,e_)&
             -qdt*wCT(ixO^S,mag(idim))*gradPsi(ixO^S)
      end if
    end do

    if (mhd_energy .and. .not.block%e_is_internal) then
    ! e = e - qdt (v . b) * div b
       w(ixO^S,e_)=w(ixO^S,e_) - qdt * divb(ixO^S) * &
            sum(v(ixO^S,:)*wCT(ixO^S,mag(:)),dim=ndim+1)
    end if
    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*v(ixO^S,idir)*divb(ixO^S)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-qdt*mhd_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_glm2')

  end subroutine add_source_glm2

  subroutine add_source_glm3(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation (1a), (1b), (4), (1d), 19
    ! conservative hyperbolic mixed GLM-MHD with no additional source terms.
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (mhd_glm_alpha < zero) then
      w(ixO^S,psi_) = abs(mhd_glm_alpha)*w(ixO^S,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab) then
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(dxlevel(:)))*w(ixO^S,psi_)
      else
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(block%ds(ixO^S,:),dim=ndim+1))*w(ixO^S,psi_)
      end if
    end if

  end subroutine add_source_glm3

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

    if (mhd_energy .and. .not.block%e_is_internal) then
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

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divb(ixI^S)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb, mhd_divb_4thorder)

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*wCT(ixO^S,mom(idir))/wCT(ixO^S,rho_)*divb(ixO^S)
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_janhunen')

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
      if(neighbor_type(i^D,saveigrid)==2 .or. neighbor_type(i^D,saveigrid)==4) then
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
       if (slab) then
          graddivb(ixp^S)=graddivb(ixp^S)*divbdiff/(^D&1.0d0/dxlevel(^D)**2+)
       else
          graddivb(ixp^S)=graddivb(ixp^S)*divbdiff &
                          /(^D&1.0d0/block%ds(ixp^S,^D)**2+)
       end if

       w(ixp^S,mag(idim))=w(ixp^S,mag(idim))+graddivb(ixp^S)

       if (mhd_energy .and. typedivbdiff=='all' .and. .not.block%e_is_internal) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixp^S,e_)=w(ixp^S,e_)+wCT(ixp^S,mag(idim))*graddivb(ixp^S)
       end if
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_linde')

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

    bvec(ixI^S,:)=w(ixI^S,mag(:))

    select case(typediv)
    case("central")
      call divvector(bvec,ixI^L,ixO^L,divb,fourthorder)
    case("limited")
      call divvectorS(bvec,ixI^L,ixO^L,divb)
    end select
  end subroutine get_divb

  !> get dimensionless div B = |divB| * volume / area / |B|
  subroutine get_normalized_divb(w,ixI^L,ixO^L,divb)

    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision                   :: divb(ixI^S), dsurface(ixI^S)

    integer :: ixA^L,idims

    call get_divb(w,ixI^L,ixO^L,divb)
    if(slab) then
      divb(ixO^S)=0.5d0*abs(divb(ixO^S))/sqrt(mhd_mag_en_all(w,ixI^L,ixO^L))/sum(1.d0/dxlevel(:))
    else
      ixAmin^D=ixOmin^D-1;
      ixAmax^D=ixOmax^D-1;
      dsurface(ixO^S)= sum(block%surfaceC(ixO^S,:),dim=ndim+1)
      do idims=1,ndim
        ixA^L=ixO^L-kr(idims,^D);
        dsurface(ixO^S)=dsurface(ixO^S)+block%surfaceC(ixA^S,idims)
      end do
      divb(ixO^S)=abs(divb(ixO^S))/sqrt(mhd_mag_en_all(w,ixI^L,ixO^L))*&
      block%dvolume(ixO^S)/dsurface(ixO^S)
    end if

  end subroutine get_normalized_divb

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixI^L,ixO^L,idirmin,current)
    use mod_global_parameters
    use mod_geometry

    integer :: idirmin0
    integer :: ixO^L, idirmin, ixI^L
    double precision :: w(ixI^S,1:nw)
    integer :: idir

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
         if(slab) then
           dtnew=min(dtnew,&
                dtdiffpar/(smalldouble+maxval(eta(ixO^S)/dxarr(idim)**2)))
         else
           dtnew=min(dtnew,&
                dtdiffpar/(smalldouble+maxval(eta(ixO^S)/block%ds(ixO^S,idim)**2)))
         end if
       end do
    end if

    if(mhd_eta_hyper>zero) then
      if(slab) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/mhd_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixO^S,1:ndim))**4/mhd_eta_hyper,dtnew)
      end if
    end if

    if(mhd_radiative_cooling) then
      call cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(mhd_viscosity) then
      call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(mhd_gravity) then
      call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

  end subroutine mhd_get_dt

  ! Add geometrical source terms to w
  subroutine mhd_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    select case (typeaxial)
    case ('cylindrical')
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in MHD")
      endif
       call mhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp)
       if(phi_>0) then
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*(tmp(ixO^S)-&
                   wCT(ixO^S,bphi_)**2+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_))
         w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt/x(ixO^S,1)*(&
                  -wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_) &
                  +wCT(ixO^S,bphi_)*wCT(ixO^S,br_))
         w(ixO^S,bphi_)=w(ixO^S,bphi_)+qdt/x(ixO^S,1)*&
                  (wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
                  -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
                  /wCT(ixO^S,rho_)
       else
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*tmp(ixO^S)
       end if
       if(mhd_glm) w(ixO^S,br_)=w(ixO^S,br_)+qdt*wCT(ixO^S,psi_)/x(ixO^S,1)
    case ('spherical')
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
           tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(idir))**2/wCT(ixO^S,rho_)-wCT(ixO^S,mag(idir))**2
           if(B0field) tmp(ixO^S)=tmp(ixO^S)-2.0d0*block%B0(ixO^S,idir,0)*wCT(ixO^S,mag(idir))
         end do
       end if
       w(ixO^S,mom(1))=w(ixO^S,mom(1))+qdt*tmp(ixO^S)/x(ixO^S,1)
       ! b1
       if(mhd_glm) then
         w(ixO^S,mag(1))=w(ixO^S,mag(1))+qdt/x(ixO^S,1)*2.0d0*wCT(ixO^S,psi_)
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
       tmp(ixO^S)=-(wCT(ixO^S,mom(1))*wCT(ixO^S,mom(2))/wCT(ixO^S,rho_) &
            -wCT(ixO^S,mag(1))*wCT(ixO^S,mag(2)))
       if (B0field) then
          tmp(ixO^S)=tmp(ixO^S)+block%B0(ixO^S,1,0)*wCT(ixO^S,mag(2)) &
               +wCT(ixO^S,mag(1))*block%B0(ixO^S,2,0)
       end if
       if(ndir==3) then
         tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(3))**2/wCT(ixO^S,rho_) &
              -wCT(ixO^S,mag(3))**2)*dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)-2.0d0*block%B0(ixO^S,3,0)*wCT(ixO^S,mag(3))&
                 *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
         end if
       end if
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
       ! b2
       tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(2)) &
            -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(1)))/wCT(ixO^S,rho_)
       if(B0field) then
         tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(1))*block%B0(ixO^S,2,0) &
              -wCT(ixO^S,mom(2))*block%B0(ixO^S,1,0))/wCT(ixO^S,rho_)
       end if
       if(mhd_glm) then
         tmp(ixO^S)=tmp(ixO^S) &
              + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
       end if
       w(ixO^S,mag(2))=w(ixO^S,mag(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
       }

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixO^S)=-(wCT(ixO^S,mom(3))*wCT(ixO^S,mom(1))/wCT(ixO^S,rho_) &
                -wCT(ixO^S,mag(3))*wCT(ixO^S,mag(1))) {^NOONED &
                -(wCT(ixO^S,mom(2))*wCT(ixO^S,mom(3))/wCT(ixO^S,rho_) &
                -wCT(ixO^S,mag(2))*wCT(ixO^S,mag(3))) &
                *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
           if (B0field) then
              tmp(ixO^S)=tmp(ixO^S)+block%B0(ixO^S,1,0)*wCT(ixO^S,mag(3)) &
                   +wCT(ixO^S,mag(1))*block%B0(ixO^S,3,0) {^NOONED &
                   +(block%B0(ixO^S,2,0)*wCT(ixO^S,mag(3)) &
                   +wCT(ixO^S,mag(2))*block%B0(ixO^S,3,0)) &
                   *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
           end if
           w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
         else
           call mpistop("angmomfix not implemented yet in MHD")
         end if
         ! b3
         tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(3)) &
              -wCT(ixO^S,mom(3))*wCT(ixO^S,mag(1)))/wCT(ixO^S,rho_) {^NOONED &
              -(wCT(ixO^S,mom(3))*wCT(ixO^S,mag(2)) &
              -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(3)))*dcos(x(ixO^S,2)) &
              /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(1))*block%B0(ixO^S,3,0) &
                 -wCT(ixO^S,mom(3))*block%B0(ixO^S,1,0))/wCT(ixO^S,rho_){^NOONED &
                 -(wCT(ixO^S,mom(3))*block%B0(ixO^S,2,0) &
                 -wCT(ixO^S,mom(2))*block%B0(ixO^S,3,0))*dcos(x(ixO^S,2)) &
                 /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
         end if
         w(ixO^S,mag(3))=w(ixO^S,mag(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
       end if
    end select
  end subroutine mhd_add_source_geom

  !> Compute 2 times total magnetic energy
  function mhd_mag_en_all(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    if (B0field) then
      mge = sum((w(ixO^S, mag(:))+block%B0(ixO^S,:,block%iw0))**2, dim=ndim+1)
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
      mgf = w(ixO^S, mag(idir))+block%B0(ixO^S,idir,block%iw0)
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
  function mhd_kin_en(w, ixI^L, ixO^L, inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)
    double precision, intent(in), optional :: inv_rho(ixO^S)

    if (present(inv_rho)) then
       ke = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) * inv_rho
    else
       ke = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) / w(ixO^S, rho_)
    end if
  end function mhd_kin_en

  subroutine mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: vHall(ixI^S,1:3)

    integer          :: idir, idirmin
    double precision :: current(ixI^S,7-2*ndir:3)

    ! Calculate current density and idirmin
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    vHall(ixO^S,1:3) = zero
    vHall(ixO^S,idirmin:3) = - mhd_etah*current(ixO^S,idirmin:3)
    do idir = idirmin, 3
       vHall(ixO^S,idir) = vHall(ixO^S,idir)/w(ixO^S,rho_)
    end do

  end subroutine mhd_getv_Hall

  subroutine mhd_getdt_Hall(w,x,ixI^L,ixO^L,dx^D,dthall)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: dthall
    !.. local ..
    double precision :: dxarr(ndim)
    double precision :: bmag(ixI^S)

    dthall=bigdouble

    ! because we have that in cmax now:
    return

    ^D&dxarr(^D)=dx^D;

    if (.not. B0field) then
       bmag(ixO^S)=sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1))
       bmag(ixO^S)=sqrt(sum((w(ixO^S,mag(:)) + block%B0(ixO^S,1:ndir,block%iw0))**2))
    end if

    if(slab) then
      dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(mhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
    else
      dthall=dtdiffpar*minval(block%ds(ixO^S,1:ndim))**2.0d0/(mhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
    end if

  end subroutine mhd_getdt_Hall

  !> This implements eq. (42) in Dedner et al. 2002 JcP 175
  !> Gives the Riemann solution on the interface
  !> for the normal B component and Psi in the GLM-MHD system.
  !> 23/04/2013 Oliver Porth
  subroutine glmSolve(wLC,wRC,ixI^L,ixO^L,idir)
    use mod_global_parameters
    double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision                :: dB(ixI^S), dPsi(ixI^S)

    dB(ixO^S)   = wRC(ixO^S,mag(idir)) - wLC(ixO^S,mag(idir))
    dPsi(ixO^S) = wRC(ixO^S,psi_) - wLC(ixO^S,psi_)

    wLC(ixO^S,mag(idir))   = 0.5d0 * (wRC(ixO^S,mag(idir)) + wLC(ixO^S,mag(idir))) &
         - 0.5d0/cmax_global * dPsi(ixO^S)
    wLC(ixO^S,psi_)       = 0.5d0 * (wRC(ixO^S,psi_) + wLC(ixO^S,psi_)) &
         - 0.5d0*cmax_global * dB(ixO^S)

    wRC(ixO^S,mag(idir)) = wLC(ixO^S,mag(idir))
    wRC(ixO^S,psi_) = wLC(ixO^S,psi_)

  end subroutine glmSolve

  subroutine mhd_boundary_adjust
    use mod_global_parameters
    integer :: iB, idim, iside, iigrid, igrid
    integer :: ixG^L, ixO^L, i^D

    ixG^L=ixG^LL;
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        saveigrid=igrid
        block=>pw(igrid)
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        do idim=1,ndim
           ! to avoid using as yet unknown corner info in more than 1D, we
           ! fill only interior mesh ranges of the ghost cell ranges at first,
           ! and progressively enlarge the ranges to include corners later
           do iside=1,2
              i^D=kr(^D,idim)*(2*iside-3);
              if (neighbor_type(i^D,igrid)/=1) cycle
              iB=(idim-1)*2+iside
              if(.not.boundary_divbfix(iB)) cycle
              if(any(typeboundary(:,iB)=="special")) then
                ! MF nonlinear force-free B field extrapolation and data driven
                ! require normal B of the first ghost cell layer to be untouched by
                ! fixdivB=0 process, set boundary_divbfix_skip(iB)=1 in par file
                select case (idim)
                {case (^D)
                   if (iside==2) then
                      ! maximal boundary
                      ixOmin^DD=ixGmax^D+1-nghostcells+boundary_divbfix_skip(2*^D)^D%ixOmin^DD=ixGmin^DD;
                      ixOmax^DD=ixGmax^DD;
                   else
                      ! minimal boundary
                      ixOmin^DD=ixGmin^DD;
                      ixOmax^DD=ixGmin^D-1+nghostcells-boundary_divbfix_skip(2*^D-1)^D%ixOmax^DD=ixGmax^DD;
                   end if \}
                end select
                call fixdivB_boundary(ixG^L,ixO^L,pw(igrid)%wb,pw(igrid)%x,iB)
              end if
           end do
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
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab) then
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
       if(slab) then
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
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(2)
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab) then
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
       if(slab) then
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
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(3)
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab) then
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
       if(slab) then
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
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(4)
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab) then
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
       if(slab) then
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
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_conserved(ixG^L,ixO^L,w,x)
     {^IFTHREED
     case(5)
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_primitive(ixG^L,ixO^L,w,x)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3+1
       if(slab) then
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
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(6)
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_primitive(ixG^L,ixO^L,w,x)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3-1
       ixFmax3=ixOmax3-1
       if(slab) then
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
       if(mhd_energy.and..not.block%e_is_internal) call mhd_to_conserved(ixG^L,ixO^L,w,x)
     }
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine fixdivB_boundary

  {^NOONED
  subroutine mhd_clean_divb_multigrid(qdt, qt, active)
    use mod_forest
    use mod_multigrid_coupling
    use mod_global_parameters
    use mod_geometry
    double precision, intent(in) :: qdt    !< Current time step
    double precision, intent(in) :: qt     !< Current time
    logical, intent(inout)       :: active !< Output if the source is active
    integer                      :: iigrid, igrid, id
    integer                      :: n, nc, lvl, ix^L, idim
    type(tree_node), pointer     :: pnode
    double precision             :: tmp(ixG^T), grad(ixG^T, ndim)
    double precision             :: res
    double precision, parameter  :: max_residual = 1d-3
    integer, parameter           :: max_its      = 10

    mg%operator_type = mg_laplacian

    ! Set boundary conditions
    do n = 1, 2*ndim
       idim = (n+1)/2
       select case (typeboundary(mag(idim), n))
       case ('symm')
          ! d/dx B = 0, take phi = 0
          mg%bc(n, mg_iphi)%bc_type = bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('asymm')
          ! B = 0, so grad(phi) = 0
          mg%bc(n, mg_iphi)%bc_type = bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('cont')
          mg%bc(n, mg_iphi)%bc_type = bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('special')
          ! Assume Dirichlet boundary conditions, derivative zero
          mg%bc(n, mg_iphi)%bc_type = bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('periodic')
          ! Nothing to do here
       case default
          print *, "divb_multigrid warning: unknown b.c.: ", &
               trim(typeboundary(mag(idim), n))
          mg%bc(n, mg_iphi)%bc_type = bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       end select
    end do

    ix^L=ixM^LL^LADD1;

    ! Store divergence of B as right-hand side
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Geometry subroutines expect this to be set
       block => pw(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       call get_divb(pw(igrid)%w(ixG^T, 1:nw), ixG^LL, ixM^LL, tmp, &
            mhd_divb_4thorder)
       mg%boxes(id)%cc({1:nc}, mg_irhs) = tmp(ixM^T)
    end do

    ! Solve laplacian(phi) = divB
    do n = 1, max_its
       call mg_fas_vcycle(mg, max_res=res)
       if (res < max_residual) exit
    end do

    if (res > max_residual) call mpistop("divb_multigrid: no convergence")

    ! Correct the magnetic field
    do iigrid = 1, igridstail
       igrid = igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id

       ! Geometry subroutines expect this to be set
       block => pw(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       ! Compute the gradient of phi
       tmp(ix^S) = mg%boxes(id)%cc({:,}, mg_iphi)
       do idim = 1, ndim
          call gradient(tmp,ixG^LL,ixM^LL,idim,grad(ixG^T, idim))
       end do

       ! Apply the correction B* = B - gradient(phi)
       tmp(ixM^T) = sum(pw(igrid)%w(ixM^T, mag(1:ndim))**2, dim=ndim+1)
       pw(igrid)%w(ixM^T, mag(1:ndim)) = &
            pw(igrid)%w(ixM^T, mag(1:ndim)) - grad(ixM^T, :)

       ! Determine magnetic energy difference
       tmp(ixM^T) = 0.5_dp * (sum(pw(igrid)%w(ixM^T, &
            mag(1:ndim))**2, dim=ndim+1) - tmp(ixM^T))

       if (mhd_energy) then
          ! Keep thermal pressure the same
          pw(igrid)%w(ixM^T, e_) = pw(igrid)%w(ixM^T, e_) + tmp(ixM^T)

          ! Possible alternative: keep total pressure the same
          ! pw(igrid)%w(ixM^T, e_)     = pw(igrid)%w(ixM^T, e_) + &
          !      (mhd_gamma-2) * inv_gamma_1 * tmp(ixM^T)
       end if
    end do

    active = .true.

  end subroutine mhd_clean_divb_multigrid
  }

end module mod_mhd_phys
