!> Special Relativistic Magneto-hydrodynamics module
module mod_srmhd_phys
  use mod_global_parameters, only: std_len
  use mod_srmhd_parameters
  use mod_srmhd_eos
  
  implicit none
  private

!  !> Whether an energy equation is used
!  logical, public, protected              :: srmhd_energy = .true.

!  !> Whether thermal conduction is used
!  logical, public, protected              :: srmhd_thermal_conduction = .false.

!  !> Whether radiative cooling is added
!  logical, public, protected              :: srmhd_radiative_cooling = .false.


!  !> Whether synge eos  is added
!  logical, public, protected              :: srmhd_eos = .false.

!  !> Whether viscosity is added
!  logical, public, protected              :: srmhd_viscosity = .false.

!  !> Whether gravity is added
!  logical, public, protected              :: srmhd_gravity = .false.

!  !> Whether Hall-MHD is used
!  logical, public, protected              :: srmhd_Hall = .false.

!  !> Whether particles module is added
!  logical, public, protected              :: srmhd_particles = .false.

!  !> Whether magnetofriction is added
!  logical, public, protected              :: srmhd_magnetofriction = .false.

!  !> Whether GLM-MHD is used
!  logical, public, protected              :: srmhd_glm = .false.

!  !> Whether divB cleaning sources are added splitting from fluid solver
!  logical, public, protected              :: source_split_divb = .false.

!  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
!  !> taking values within [0, 1]
!  double precision, public                :: srmhd_glm_alpha = 0.5d0

!  !> MHD fourth order
!  logical, public, protected              :: srmhd_4th_order = .false.

!  !> Number of tracer species
!  integer, public, protected              :: srmhd_n_tracer = 0

!  !> Index of the density (in the w array)
!  integer, public, protected              :: rho_
!
!  !> Indices of the momentum density
!  integer, allocatable, public, protected :: mom(:)
!
!  !> Index of the energy density (-1 if not present)
!  integer, public, protected              :: e_

!  !> Index of the gas pressure (-1 if not present) should equal e_
!  integer, public, protected              :: p_

!  !> Indices of the magnetic field
!  integer, allocatable, public, protected :: mag(:)

!  !> Indices of the GLM psi
!  integer, public, protected :: psi_

!  !> Indices of the tracers
!  integer, allocatable, public, protected :: tracer(:)
!
!  !> The adiabatic index
!  double precision, public                :: srmhd_gamma = 5.d0/3.0d0
!
!  !> The adiabatic constant
!  double precision, public                :: srmhd_adiab = 1.0d0
!
!  !> The MHD resistivity
!  double precision, public                :: srmhd_eta = 0.0d0

!  !> The MHD hyper-resistivity
!  double precision, public                :: srmhd_eta_hyper = 0.0d0

!  !> TODO: what is this?
!  double precision, public                :: srmhd_etah = 0.0d0


!  logical,          public, protected     :: srmhd_checkNR = .true.
!  double precision, public, protected     :: srmhd_absaccnr=1.0d-8
!  double precision, public, protected     :: srmhd_tolerNR =1.0d-9

!  !> The small_est allowed energy
!  double precision, protected             :: small_e


!  !> The small_est allowed inertia
!  double precision, protected             :: small_xi
  
!  ! smaller values for speed
!  double precision, public, protected             :: small_vec2  = 0.0
!  !> The number of waves
!  integer :: nwwave=8

!  !> Method type to clean divergence of B
!  character(len=std_len), public, protected :: typedivbfix  = 'linde'
!
!  !> Method type in a integer for good performance
!  integer :: type_divb
!
!  !> Coefficient of diffusive divB cleaning
!  double precision :: divbdiff     = 0.8d0
!
!  !> Update all equations due to divB cleaning
!  character(len=std_len) ::    typedivbdiff = 'all'

!  !> Use a compact way to add resistivity
!  logical :: compactres   = .false.

!  !> Add divB wave in Roe solver
!  logical, public :: divbwave     = .true.

!  !> Helium abundance over Hydrogen
!  double precision, public, protected  :: He_abundance=0.1d0

!  !> To control divB=0 fix for boundary
!  logical, public, protected :: boundary_divbfix(2*^ND)=.true.

!  !> To skip * layer of ghost cells during divB=0 fix for boundary
!  integer, public, protected :: boundary_divbfix_skip(2*^ND)=0

!  !> B0 field is force-free
!  logical, public, protected :: B0field_forcefree=.true.


!  ! DivB cleaning methods
!  integer, parameter :: divb_none          = 0
!  integer, parameter :: divb_glm1          = 1
!  integer, parameter :: divb_glm2          = 2
!  integer, parameter :: divb_powel         = 3
!  integer, parameter :: divb_janhunen      = 4
!  integer, parameter :: divb_linde         = 5
!  integer, parameter :: divb_lindejanhunen = 6
!  integer, parameter :: divb_lindepowel    = 7
!  integer, parameter :: divb_lindeglm      = 8

  ! Public methods
  public :: srmhd_phys_init
  public :: srmhd_kin_en_primitive
  public :: srmhd_get_pthermal
  public :: srmhd_get_p_mag
  public :: srmhd_get_p_total
  public :: srmhd_get_v
  public :: srmhd_get_v_idim
  public :: srmhd_to_conserved
  public :: srmhd_to_primitive
  public :: srmhd_get_csound2
  public :: srmhd_get_csound_prim
  public :: get_divb
  public :: get_current
  public :: get_normalized_divb
  public :: srmhd_get_4u_from_3v
contains

  !> Read this module"s parameters from a file
  subroutine srmhd_read_params(files)
    ! made by Z. Meliani 20/02/2018 
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srmhd_list/ srmhd_energy, srmhd_eos,srmhd_n_tracer, srmhd_gamma, &
                          srmhd_adiab, srmhd_eta, srmhd_eta_hyper, &
                          srmhd_etah, srmhd_glm_alpha, srmhd_magnetofriction,&
                          srmhd_thermal_conduction, srmhd_radiative_cooling, &
                          srmhd_Hall, srmhd_gravity, srmhd_viscosity, &
                          srmhd_4th_order, typedivbfix, source_split_divb, &
                          divbdiff, typedivbdiff, compactres, divbwave, srmhd_glm,&
                          He_abundance, SI_unit, B0field,&
                          B0field_forcefree, Bdip, Bquad, Boct, Busr, &
                          srmhd_particles, boundary_divbfix, &
                          boundary_divbfix_skip, &
                          srmhd_maxiterationNR, srmhd_absaccNR,srmhd_tolerNr,&
                          srmhd_checkNR,srmhd_maxdspeed,small_vec2
    write(*,*)'Reading srmhd_list'
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, srmhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine srmhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine srmhd_write_info(fh)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = srmhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine srmhd_write_info

  subroutine srmhd_angmomfix(fC,x,wnew,ixI^L,ixO^L,idim)
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

  end subroutine srmhd_angmomfix

  subroutine srmhd_phys_init()
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_magnetofriction, only: magnetofriction_init
    use mod_physics

    integer :: itr, idir
    unit_velocity=const_c
    call srmhd_read_params(par_files)

    physics_type = "srmhd"
    phys_energy=srmhd_energy
    ! set default gamma for polytropic/isothermal process
    if(.not.srmhd_energy) srmhd_gamma=1.d0
    use_particles=srmhd_particles
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
      type_divb = divb_none
    case ('glm1')
      srmhd_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_glm1
    case ('glm2')
      srmhd_glm          = .true.
      need_global_cmax = .true.
      need_global_vmax = .true.
      type_divb        = divb_glm2
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
      srmhd_glm          = .true.
      need_global_cmax = .true.
      need_global_vmax = .true.
      type_divb        = divb_lindeglm
    case default
      call mpistop('Unknown divB fix')
    end select

    ! Determine flux variables
    rho_ = var_set_rho()
    d_=rho_
    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (srmhd_energy) then
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
    if (srmhd_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    allocate(tracer(srmhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, srmhd_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! Set index for auxiliary variables
    xi_  = var_set_auxvar('xi')
    lfac_=var_set_auxvar('lfac_')
    !nwfluxbc=nwfluxbc+2

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
    if(srmhd_glm .and. ndim>1) flux_type(:,psi_)=flux_tvdlf

    srmhd_maxspeed=1.0-srmhd_maxdspeed


    phys_get_dt              => srmhd_get_dt
    phys_get_cmax            => srmhd_get_cmax
    phys_get_cbounds         => srmhd_get_cbounds
    phys_get_flux            => srmhd_get_flux
    phys_get_v_idim          => srmhd_get_v_idim
    phys_add_source_geom     => srmhd_add_source_geom
    phys_add_source          => srmhd_add_source
    phys_to_conserved        => srmhd_to_conserved
    phys_to_primitive        => srmhd_to_primitive
    phys_get_aux             => srmhd_get_auxiliary
    phys_check_params        => srmhd_check_params
    phys_check_w             => srmhd_check_w
    phys_get_pthermal        => srmhd_get_pthermal
    phys_boundary_adjust     => srmhd_boundary_adjust
    phys_write_info          => srmhd_write_info
    phys_angmomfix           => srmhd_angmomfix
    phys_handle_small_values => srmhd_handle_small_values

    ! Whether diagonal ghost cells are required for the physics
    if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! derive units from basic units
    call srmhd_physical_units()

    if(.not. srmhd_energy .and. srmhd_thermal_conduction) then
      call mpistop("thermal conduction needs srmhd_energy=T")
    end if
    if(.not. srmhd_energy .and. srmhd_radiative_cooling) then
      call mpistop("radiative cooling needs srmhd_energy=T")
    end if

    ! initialize thermal conduction module
    !if (srmhd_thermal_conduction) then
    !  phys_req_diagonal = .true.
    !  call thermal_conduction_init(srmhd_gamma)
    !end if

    ! Initialize radiative cooling module
    !if (srmhd_radiative_cooling) then
    !  call radiative_cooling_init(srmhd_gamma,He_abundance)
    !end if

    ! Initialize viscosity module
    !if (srmhd_viscosity) call viscosity_init(phys_wider_stencil,phys_req_diagonal)

    ! Initialize gravity module
    if(srmhd_gravity) then
      call gravity_init()
    end if

    ! Initialize particles module
    if(srmhd_particles) then
      call particles_init()
      phys_req_diagonal = .true.
    end if

    ! initialize magnetofriction module
    !if(srmhd_magnetofriction) then
    !  phys_req_diagonal = .true.
    !  call magnetofriction_init()
    !end if

    if(type_divb==divb_glm1) then
      ! Solve the Riemann problem for the linear 2x2 system for normal
      ! B-field and GLM_Psi according to Dedner 2002:
      phys_modify_wLR => glmSolve
    end if

    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in getflux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in nghostcells.
    !if (srmhd_hall) then
    !   if (srmhd_4th_order) then
    !      phys_wider_stencil = 2
    !   else
    !      phys_wider_stencil = 1
    !   end if
    !end if

  end subroutine srmhd_phys_init

  subroutine srmhd_check_params
    use mod_global_parameters

    ! after user parameter setting
    gamma_1=srmhd_gamma-1.0_dp
    inv_gamma_1=1.d0/gamma_1
    gamma_to_gamma_1=srmhd_gamma/gamma_1
    if (.not. srmhd_energy) then
       if (srmhd_gamma <= 0.0d0) call mpistop ("Error: srmhd_gamma <= 0")
       if (srmhd_adiab < 0.0d0) call mpistop ("Error: srmhd_adiab < 0")
       small_pressure = srmhd_adiab*small_density**srmhd_gamma
    else
       if (srmhd_gamma <= 0.0d0 .or. srmhd_gamma == 1.0d0) &
            call mpistop ("Error: srmhd_gamma <= 0 or srmhd_gamma == 1")
       small_e = small_pressure * inv_gamma_1
    end if
    small_xi=small_density+gamma_to_gamma_1*small_pressure

  end subroutine srmhd_check_params

  subroutine srmhd_physical_units()
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
      call mpistop ("Error: in srmhd the unit_velocity=c")
    else
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*kB)
      unit_magneticfield=sqrt(miu0*unit_pressure)
      unit_time=unit_length/unit_velocity
    end if

  end subroutine srmhd_physical_units

  subroutine srmhd_check_w(primitive,ixI^L,ixO^L,w,flag)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters

    logical, intent(in)            :: primitive
    integer, intent(in)            :: ixI^L, ixO^L
    double precision, intent(in)   :: w(ixI^S,nw)
    integer, intent(inout)         :: flag(ixI^S)
    double precision :: tmp(ixI^S)

    flag(ixO^S)=0
    where(w(ixO^S, rho_)/w(ixO^S,lfac_) < small_density) flag(ixO^S) = rho_

    cond_energy : if (srmhd_energy) then
       cond_einternal : if (block%e_is_internal) then
          where(w(ixO^S, e_) < small_pressure*inv_gamma_1) flag(ixO^S) = e_
       else cond_einternal
         is_prim : if (primitive)then
           where(w(ixO^S, p_) < small_pressure) flag(ixO^S) = p_ ! p_=e_
         else is_prim
           ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
           ! in srmhd we check if xi-p-d is positive
             ! do be done
             !where(tmp(ixO^S) < small_pressure) flag(ixO^S) = e_
         end if is_prim
       end if cond_einternal
    end if cond_energy
  end subroutine srmhd_check_w

  !> Transform primitive variables into conservative ones
  subroutine srmhd_to_conserved(ixI^L,ixO^L,w,x)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(inout)    :: w(ixI^S, nw)
    double precision, intent(in)       :: x(ixI^S, 1:ndim)
    integer                            :: idir, itr
    double precision, dimension(ixO^S) :: sqrU,B2,VdotB,rhoh,sqrV
    integer, dimension(ixO^S)          :: flag_error    
    character(len=30)                  :: subname_loc
 
   
    subname_loc= 'srmhd_to_conserved' 
    flag_error(ixO^S)=0
    where(w(ixO^S,rho_)<small_density)flag_error(ixO^S)=rho_


    sqrU(ixO^S)    = sum(w(ixO^S, mom(:))**2, dim=ndim+1)
    w(ixO^S,lfac_) = dsqrt(1.0_dp+sqrU(ixO^S))
    sqrV=sqrU/w(ixO^S,lfac_)
    ! fill the auxiliary variable xi and density D
    call srmhd_get_enthalpy(ixO^L,w(ixO^S,rho_),w(ixO^S,p_),rhoh)
    ! with enthalpy w: xi= lfac^2 rhoh
    w(ixO^S,xi_) = w(ixO^S,lfac_)**2.0D0*rhoh(ixO^S)
    ! density: d = lfac * rho
    w(ixO^S,rho_)=w(ixO^S,lfac_)*w(ixO^S,rho_) 
    
    call srmhd_get_B2andVdotB(ixI^L,ixO^L,w,.false.,B2=B2,VdotB=VdotB)
    ! Convert velocity to momentum
    ! s= (xi + B^2) * v - (v.B) * B
    ! re-use rhoh as rhoh=(xi+B^2)/lfac
    rhoh(ixO^S)= (w(ixO^S,xi_)+B2(ixO^S))/w(ixO^S,lfac_)
    Loop_idirmom : do idir = 1, ndir
       w(ixO^S, mom(idir)) = rhoh(ixO^S) * w(ixO^S, mom(idir))&
                             -VdotB(ixO^S)*w(ixO^S,mag(idir))
    end do Loop_idirmom

    cond_energy : if (srmhd_energy) then
       ! re-use sqrU=v^2 B^2 - (v.B)^2
       sqrU(ixO^S) = B2(ixO^S)*sqrV(ixO^S)-VdotB(ixO^S)**2.0_dp
       ! sqrU should positive
       where(sqrU(ixO^S)<0.0_dp) flag_error(ixO^S)=e_
       ! re-use sqrV=xi - p -D
       sqrV(ixO^S) = w(ixO^S,xi_) - w(ixO^S,p_) - w(ixO^S,d_)
       where(sqrV(ixO^S)<0.0_dp) flag_error(ixO^S)=e_
       ! E = xi - p +(B^2+v^2 B^2 - (v.B)^2)/2- D
       w(ixO^S,e_)=sqrV(ixO^S) +0.5_dp*(B2(ixO^S) + sqrU(ixO^S))
      !if(.not.block%e_is_internal) w(ixO^S,e_)=w(ixO^S,e_) + x
      if(type_divb==divb_glm2) w(ixO^S,e_)=w(ixO^S,e_) &
                      + 0.5d0*w(ixO^S,psi_)**2
    end if cond_energy
    if (check_small_values) call srmhd_handle_small_values(.false., &
                                  w, x, ixI^L, ixO^L,trim(subname_loc),&
                                  flag_error=flag_error)
  end subroutine srmhd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine srmhd_to_primitive(ixI^L,ixO^L,w,x)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision                :: inv_rho(ixI^S),B2(ixO^S),VdotB(ixO^S)
    integer                         :: itr, idir
    character(len=30)               :: subname_loc
 
    subname_loc='srmhd_to_primitive'

    ! get auxiliary variables
    call srmhd_get_auxiliary(.true.,w,x,ixI^L,ixO^L,subname_loc)
    w(ixO^S,rho_) = w(ixO^S,d_)/w(ixO^S,lfac_)
    if (srmhd_energy) then
      ! Calculate pressure = (gamma-1) * (e-ek-eb)
      if(.not.block%e_is_internal) then
        call srmhd_get_pressure_fromprimitive(ixI^L,ixI^L,w,inv_rho)
        w(ixO^S,p_)   = inv_rho(ixO^S)
      else
       w(ixO^S,p_)   = gamma_1*w(ixO^S, e_)
      end if
    end if
    ! re-use inv_rho to store inverse of the density
    inv_rho = 1.0d0 / w(ixO^S, rho_)

    call srmhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2,VdotB)
    Loop_idir : do idir=1,ndir
     w(ixO^S,mom(idir)) = w(ixO^S,lfac_)*(w(ixO^S,mom(idir))&
                      +VdotB*w(ixO^S,mag(idir)))&
                      /(w(ixO^S,xi_)+B2)
    end do Loop_idir

    if (check_small_values) call srmhd_handle_small_values(.true., w, x&
                                 , ixI^L, ixO^L,trim(subname_loc))
  end subroutine srmhd_to_primitive

  subroutine srmhd_handle_small_values(primitive, w, x, ixI^L, ixO^L, &
                                       subname,flag_error)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    integer, optional, intent(in)   :: flag_error(ixO^S)               

    double precision :: smallone
    integer :: idir, flag(ixI^S)

    if (small_values_method == "ignore") return
    if(present(flag_error)) then
     flag(ixO^S) = flag_error(ixO^S)
    else
     call srmhd_check_w(primitive, ixI^L, ixO^L, w, flag)
    end if
    if (any(flag(ixO^S) /= 0)) then
       select case (small_values_method)
       case ("replace")
          where(flag(ixO^S) /= 0) w(ixO^S,rho_) = small_density

          do idir = 1, ndir
             where(flag(ixO^S) /= 0) w(ixO^S, mom(idir)) = 0.0d0
          end do

          if (srmhd_energy) then
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
  end subroutine srmhd_handle_small_values




 !> Calculate thermal pressure for enthalpy and density
  subroutine srmhd_get_pressure_fromprimitive(ixI^L,ixO^L,w,pth)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    implicit none
    integer, intent(in)            :: ixI^L, ixO^L
    double precision, intent(in)   :: w(ixI^S,nw)
    double precision, intent(out)  :: pth(ixI^S)

    double precision               :: rhoh(ixO^S)

    rhoh  = w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srmhd_get_pressure_primitive_eos(ixI^L,ixO^L,w(ixO^S,rho_),rhoh,pth) 
    !pth(ixO^S) =  w(ixO^S,p_)
  end subroutine srmhd_get_pressure_fromprimitive


  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine srmhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    implicit none

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)

    double precision             :: rho(ixO^S),rhoh(ixO^S)

    rho        = w(ixO^S,d_)/w(ixO^S,lfac_)
    rhoh =w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srmhd_get_pthermal_eos(ixI^L,ixO^L,x,rho,rhoh,w(ixO^S,e_),pth)
  end subroutine srmhd_get_pthermal


 !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamm*p/rho
  subroutine srmhd_get_csound2_prim(w,x,ixI^L,ixO^L,csound2)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters
    implicit none
    integer, intent(in)               :: ixI^L, ixO^L
    double precision, intent(in)      :: w(ixI^S,nw)
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    double precision, intent(out)     :: csound2(ixO^S)

    double precision,dimension(ixO^S) :: rhoh

    
    rhoh=w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    
    call srmhd_get_csound2_prim_eos(ixI^L,ixO^L,x,w(ixO^S,rho_),&
                                    rhoh,w(ixO^S,p_),csound2)
  end subroutine srmhd_get_csound2_prim

  !> Convert energy to entropy
  subroutine e_to_rhos(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision,intent(inout)  :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    if (srmhd_energy) then
      if(.not.block%e_is_internal) &
        w(ixO^S, e_) = w(ixO^S, e_) - srmhd_kin_en_primitive(w, ixI^L, ixO^L) &
              - srmhd_mag_en_primitive(w, ixI^L, ixO^L)
      w(ixO^S, e_) = gamma_1* w(ixO^S, rho_)**(1.0d0 - srmhd_gamma) * &
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

    if (srmhd_energy) then
       w(ixO^S, e_) = w(ixO^S, rho_)**gamma_1 * w(ixO^S, e_) &
            * inv_gamma_1
       if(.not.block%e_is_internal) &
         w(ixO^S, e_) =w(ixO^S, e_) + srmhd_kin_en_primitive(w, ixI^L, ixO^L) + &
            srmhd_mag_en_primitive(w, ixI^L, ixO^L)
    else
       call mpistop("rhos_to_e can not be used without energy equation!")
    end if
  end subroutine rhos_to_e

  !> Calculate v vector
  subroutine srmhd_get_v(w,x,ixI^L,ixO^L,v,B2,VdotB)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,1:ndir)

    double precision, optional, intent(in)  :: VdotB(ixO^S),B2(ixO^S)
    double precision                        :: sub_VdotB(ixO^S)&
                                              ,sub_B2(ixO^S)
    integer :: idir



    is_B2_in : if(present(B2)) then
     Loop_idir_b_v1 : do idir=1,ndir
      v(ixO^S,idir) = (w(ixO^S, mom(idir)) + VdotB*w(ixO^S,mag(idir)))&
                      /(w(ixO^S,xi_)+B2)
     end do Loop_idir_b_v1
    else  is_B2_in
     call srmhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     Loop_idir_b_v2: do idir=1,ndir
      v(ixO^S,idir) = (w(ixO^S, mom(idir))+sub_VdotB*w(ixO^S,mag(idir)))&
                      /(w(ixO^S,xi_)+sub_B2)
     end do Loop_idir_b_v2
    end if is_B2_in
  end subroutine srmhd_get_v



  !> Calculate v component
  subroutine srmhd_get_v_idim(w,x,ixI^L,ixO^L,idim,v_idim)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters

    integer, intent(in)                    :: ixI^L, ixO^L, idim
    double precision, intent(in)           :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out)          :: v_idim(ixI^S)

    double precision                       :: sub_VdotB(ixO^S),sub_B2(ixO^S)

     call srmhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     v_idim(ixO^S) = (w(ixO^S, mom(idim)) + sub_VdotB*w(ixO^S,mag(idim)))&
                      /(w(ixO^S,xi_)+sub_B2)
  end subroutine srmhd_get_v_idim

  !> Calculate v component
  subroutine srmhd_get_v_idim_loc(w,x,ixI^L,ixO^L,idim,v_idim,B2,VdotB)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters

    integer, intent(in)                    :: ixI^L, ixO^L, idim
    double precision, intent(in)           :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out)          :: v_idim(ixI^S)
    double precision, optional, intent(in) :: VdotB(ixO^S),B2(ixO^S)

    double precision                       :: sub_VdotB(ixO^S),sub_B2(ixO^S)

    if(present(B2))then
     v_idim(ixO^S) = (w(ixO^S, mom(idim)) + VdotB*w(ixO^S,mag(idim)))&
                      /(w(ixO^S,xi_)+B2)
    else
     call srmhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     v_idim(ixO^S) = (w(ixO^S, mom(idim)) + sub_VdotB*w(ixO^S,mag(idim)))&
                      /(w(ixO^S,xi_)+sub_B2)
    end if
  end subroutine srmhd_get_v_idim_loc

 !> Calculate v^2
  subroutine srmhd_get_v2(w,x,ixI^L,ixO^L,v2,B2,VdotB)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v2(ixO^S)

    double precision, optional, intent(in)  :: VdotB(ixO^S),B2(ixO^S)
    double precision                        :: sub_VdotB(ixO^S),sub_B2(ixO^S)
    double precision                        :: v_num(ixO^S,1:ndir)
    integer :: idir

    if(present(B2)) then
     Loop_idir_B2: do idir=1,ndir
       v_num(ixO^S,idir)=w(ixO^S, mom(idir)) + VdotB*w(ixO^S,mag(idir))
     end do  Loop_idir_B2
     v2(ixO^S) = sum(v_num(ixO^S, 1:ndir)**2.0 ,dim=ndim+1)&
                      /(w(ixO^S,xi_)+B2)**2.0
    else 
     call srmhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     Loop_idir_noB2: do idir=1,ndir
       v_num(ixO^S,idir)=w(ixO^S, mom(idir)) + sub_VdotB*w(ixO^S,mag(idir))
     end do  Loop_idir_noB2
     
     v2(ixO^S) = sum(v_num(ixO^S, :)**2.0 ,dim=ndim+1)&
                      /(w(ixO^S,xi_)+sub_B2)**2.0
    end if
  end subroutine srmhd_get_v2

  subroutine srmhd_get_cmax_gammie(ixI^L,ixO^L,idim,x,w,rhoh,B2,&
                                   VdotB,calfven,vidim,csound2,cmax&
                                   ,patch_gammie,from_cbound,cmin)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    implicit none
    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(in), dimension(ixO^S) :: rhoh,B2,VdotB&
                                                     ,calfven, csound2
    double precision, intent(inout)          :: cmax(ixI^S)
    double precision, optional, intent(inout):: cmin(ixI^S)
    logical, optional, intent(in)            :: from_cbound
    logical, optional, intent(in)            :: patch_gammie(ixO^S)

    double precision, dimension(ixO^S):: A,B
    double precision, dimension(ixO^S):: vidim,v2
   is_patch : if(present(patch_gammie))then
    where(patch_gammie(ixO^S))
     A(ixO^S) = csound2(ixO^S)+calfven(ixO^S)-csound2(ixO^S)*calfven(ixO^S)
     B(ixO^S)=vidim(ixO^S)**2.0d0
    end where
    if(.not.present(from_cbound))&
           where(patch_gammie(ixO^S))vidim(ixO^S)=dabs(vidim(ixO^S))
    cond_onedir_patch: if(ndir==1)then
     where(patch_gammie(ixO^S))cmax(ixO^S)=(vidim(ixO^S)+dsqrt(A))&
                                           /(1.0d0+dsqrt(A*B))
     if(present(cmin))&
        where(patch_gammie(ixO^S))cmin(ixO^S)=(vidim(ixO^S)-dsqrt(A))&
                                               /(1.0d0+dsqrt(A*B))
    else cond_onedir_patch

     call srmhd_get_v2(w,x,ixI^L,ixO^L,v2,B2=B2,VdotB=VdotB)
     where(patch_gammie(ixO^S))
      cmax(ixO^S)=(vidim(ixO^S)*(1.0d0-A)&
                  +dsqrt(A*(1.0d0-v2(ixO^S))*((1.0d0-v2(ixO^S)*A)&
                         -B(ixO^S)*(1.0d0-A))))&
                 /( 1.0d0-v2(ixO^S)*A)
     end where
     if(present(cmin))&
        where(patch_gammie(ixO^S))cmin(ixO^S)=(vidim(ixO^S)*(1.0d0-A)&
                  -dsqrt(A*(1.0d0-v2(ixO^S))*(1.0d0-v2(ixO^S)*A&
                         -B(ixO^S)*(1.0d0-A))))&
                 /( 1.0d0-v2(ixO^S)*A)
    end if cond_onedir_patch    
   else is_patch 
    A = csound2(ixO^S)+calfven(ixO^S)-csound2(ixO^S)*calfven(ixO^S)
    ! use B to save vidim**2
    if(.not.present(from_cbound).or..not.present(cmin))&
                                        vidim(ixO^S)=dabs(vidim(ixO^S))
    B(ixO^S)=vidim(ixO^S)**2.0d0
    cond_onedir: if(ndir==1)then
     cmax(ixO^S)=(vidim(ixO^S)+dsqrt(A))/(1.0d0+dsqrt(A*B))
     if(present(cmin))cmin(ixO^S)=(vidim(ixO^S)-dsqrt(A))/(1.0d0+dsqrt(A*B))
    else cond_onedir

     call srmhd_get_v2(w,x,ixI^L,ixO^L,v2,B2=B2,VdotB=VdotB)
     cmax(ixO^S)=(vidim(ixO^S)*(1.0d0-A)&
                  +dsqrt(A*(1.0d0-v2(ixO^S))*((1.0d0-v2(ixO^S)*A)&
                         -B(ixO^S)*(1.0d0-A))))&
                 /( 1.0d0-v2(ixO^S)*A)
     if(present(cmin))cmin(ixO^S)=(vidim(ixO^S)*(1.0d0-A)&
                  -dsqrt(A*(1.0d0-v2(ixO^S))*(1.0d0-v2(ixO^S)*A&
                         -B(ixO^S)*(1.0d0-A))))&
                 /( 1.0d0-v2(ixO^S)*A)
    end if cond_onedir

   end if is_patch 

  end subroutine srmhd_get_cmax_gammie


  !> Calculate cmax_idim using gammie method within ixO^L
  subroutine srmhd_get_cmax(w,x,ixI^L,ixO^L,idim,cmax)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    
    double precision, dimension(ixO^S):: rhoh,B2,VdotB,calfven,vidim
    double precision, dimension(ixO^S):: v2,csound2

    rhoh=w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srmhd_get_csound2(w,x,ixI^L,ixO^L,rhoh,csound2)
    call srmhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2,VdotB)
    call srmhd_get_calfven2(ixI^L,ixO^L,x,w,.true.,calfven,rhoh=rhoh,B2=B2)
    call srmhd_get_v_idim_loc(w,x,ixI^L,ixO^L,idim,vidim,B2=B2,VdotB=VdotB)

    call srmhd_get_cmax_gammie(ixI^L,ixO^L,idim,x,w,rhoh=rhoh,B2=B2,&
                               VdotB=VdotB,calfven=calfven,vidim=vidim,&
                               csound2=csound2,cmax=cmax)


!    A = cmax(ixO^S)+B(ixO^S)-cmax(ixO^S)*B(ixO^S)
!    call srmhd_get_v_idim(w,x,ixI^L,ixO^L,idim,vidim,B2=B2,VdotB=C)
!    ! use B to save vidim**2
!    B(ixO^S)=vidim(ixO^S)**2.0d0
!    cond_onedir: if(ndir==1)then
!     cmax(ixO^S)=(dabs(vidim(ixO^S))+dsqrt(A))/(1.0d0+dsqrt(A+B))
!    else cond_onedir
!   
!     call srmhd_get_v2(w,x,ixI^L,ixO^L,v2,B2=B2,VdotB=C)
!     cmax(ixO^S)=(dabs(vidim(ixO^S))*(1.0d0-A)&
!                  +dsqrt(A*(1.0d0-v2(ixO^S))*((1.0d0-v2(ixO^S)*A)&
!                         -B(ixO^S)*(1.0d0-A))))&
!                 /( 1.0d0-v2(ixO^S)*A)
!    end if cond_onedir
    
  end subroutine srmhd_get_cmax

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine srmhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,cmax,cmin)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)              :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)              :: x(ixI^S,1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision                   :: wmean(ixI^S,nw)
    double precision, dimension(ixO^S) :: umean, dmean, csound2Lp, &
                                          csound2Rp, tmp1,tmp2,tmp3, &
                                          B2Lp,B2Rp,VdotBLp,VdotBRp, &
                                          rhohLp,rhohRp, &
                                          cmaxL,cmaxR,csound2,rhoh,&
                                          B2,VdotB, vidim,vidimLp,vidimRp,&
                                          calfvenLp,calfvenRp,calfven
    character(len=30)                  :: subname_loc

    subname_loc='srmhd_get_cbounds'
    if (typeboundspeed/='cmaxmean') then
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      rhohLp=wLp(ixO^S,xi_)/wLp(ixO^S,lfac_)**2.0d0
      rhohRp=wRp(ixO^S,xi_)/wRp(ixO^S,lfac_)**2.0d0

      call srmhd_get_csound2(wLp,x,ixI^L,ixO^L,rhohLp,csound2Lp)
      call srmhd_get_csound2(wRp,x,ixI^L,ixO^L,rhohRp,csound2Rp)
      call srmhd_get_B2andVdotB(ixI^L,ixO^L,wLp,.false.,B2=B2Lp,VdotB=VdotBLp)
      call srmhd_get_B2andVdotB(ixI^L,ixO^L,wRp,.false.,B2=B2Rp,VdotB=VdotBRp)

      call srmhd_get_calfven2(ixI^L,ixO^L,x,wLp,.false.,calfvenLp,&
                              rhoh=rhohLp,B2=B2Lp)
      call srmhd_get_calfven2(ixI^L,ixO^L,x,wRp,.false.,calfvenRp,&
                              rhoh=rhohRp,B2=B2Rp)

      call srmhd_get_v_idim_loc(wLp,x,ixI^L,ixO^L,idim,vidimLp,B2=B2Lp,&
                                VdotB=VdotBLp)
      call srmhd_get_v_idim_loc(wRp,x,ixI^L,ixO^L,idim,vidimRp,B2=B2Rp,&
                                VdotB=VdotBRp)

      tmp1(ixO^S)=sqrt(wLp(ixO^S,xi_)+B2Lp(ixO^S))
      tmp2(ixO^S)=sqrt(wRp(ixO^S,xi_)+B2Rp(ixO^S))
      tmp3(ixO^S)=1.d0/(sqrt(wLp(ixO^S,xi_)+B2Lp(ixO^S))&
                  +sqrt(wRp(ixO^S,rho_)+B2Rp(ixO^S)))

      umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)&
                    +wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)

     call srmhd_get_cmax_gammie(ixI^L,ixO^L,idim,x,wLp,rhoh=rhohLp,B2=B2Lp,&
                                VdotB=VdotBLp,calfven=calfvenLp,vidim=vidimLp,&
                                csound2=csound2Lp,cmax=cmaxL)

     call srmhd_get_cmax_gammie(ixI^L,ixO^L,idim,x,wRp,rhoh=rhohRp,B2=B2Rp,&
                                VdotB=VdotBRp,calfven=calfvenRp,vidim=vidimRp,&
                                csound2=csound2Rp,cmax=cmaxR)


      dmean(ixO^S)=(tmp1(ixO^S)*cmaxL(ixO^S)**2&
                    +tmp2(ixO^S)*cmaxR(ixO^S)**2)*tmp3(ixO^S)+&
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
      ! get auxiliary variables
      call srmhd_get_auxiliary(.true.,wmean,x,ixI^L,ixO^L,subname_loc)
      rhoh=wmean(ixO^S,xi_)/wmean(ixO^S,lfac_)**2.0d0

      call srmhd_get_csound2(wmean,x,ixI^L,ixO^L,rhoh,csound2)
      call srmhd_get_B2andVdotB(ixI^L,ixO^L,wmean,.true.,B2=B2,VdotB=VdotB)
      call srmhd_get_calfven2(ixI^L,ixO^L,x,wmean,.true.,&
                               calfven,rhoh=rhoh,B2=B2) 
      call srmhd_get_v_idim_loc(wmean,x,ixI^L,ixO^L,idim,vidim,B2=B2,&
                                VdotB=VdotB)
      call srmhd_get_cmax_gammie(ixI^L,ixO^L,idim,x,wmean,rhoh=rhoh,B2=B2,&
                                 VdotB=VdotB,calfven=calfven,&
                                 vidim=vidim,csound2=csound2,cmax=cmax&
                                 ,from_cbound=.true.,cmin=cmin)
      if(present(cmin)) then
        cmax(ixO^S)=min(max(cmax(ixO^S),0.0D0),1.0d0)
        cmin(ixO^S)=max(min(cmin(ixO^S),0.0D0),-1.0d0)
      else
        cmax(ixO^S)=min(cmax(ixO^S),1.0d0)
      end if
    end if
  end subroutine srmhd_get_cbounds
    ! made by Z. Meliani 10/02/2018 
  !> Calculate the Aflven speed
  subroutine srmhd_get_calfven2(ixI^L,ixO^L,x,w,conserve,calfven,rhoh,B2)
   use mod_global_parameters
   implicit none
   integer, intent(in)           :: ixI^L, ixO^L
   double precision, intent(in)  :: w(ixI^S, 1/nw), x(ixI^S,1:ndim)
   logical,          intent(in)  :: conserve
   double precision, intent(out) :: calfven(ixO^S)
   double precision, optional    :: rhoh(ixO^S),B2(ixO^S)

   double precision              :: sub_B2(ixO^S),sub_rhoh(ixO^S)

   if(present(B2))then
    calfven=B2/(B2+rhoh)
   else
    sub_rhoh=w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srmhd_get_B2andVdotB(ixI^L,ixO^L,w,conserve,B2=sub_B2)
    calfven=sub_B2/(sub_B2+sub_rhoh)
   end if
  end subroutine srmhd_get_calfven2


  !> Calculate fast magnetosonic wave speed
!  subroutine srmhd_get_cfast(w,x,ixI^L,ixO^L,idim,csound)
!    use mod_global_parameters

!    integer, intent(in)          :: ixI^L, ixO^L, idim
!    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
!    double precision, intent(out):: csound(ixI^S)
!    double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
!    double precision :: inv_rho(ixO^S)

!    inv_rho=1.d0/w(ixO^S,rho_)

!    call srmhd_get_csound2(w,x,ixI^L,ixO^L,csound)
!    ! store |B|^2 in v
!    b2(ixO^S)        = srmhd_mag_en_all(w,ixI^L,ixO^L)
!    cfast2(ixO^S)   = b2(ixO^S) * inv_rho+csound(ixO^S)
!    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
!         * srmhd_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
!         * inv_rho

!    where(AvMinCs2(ixO^S)<zero)
!       AvMinCs2(ixO^S)=zero
!    end where

!    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))
!
!    if (.not. MHD_Hall) then
!       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
!    else
!       ! take the Hall velocity into account:
!       ! most simple estimate, high k limit:
!       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
!       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
!       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
!            srmhd_etah * sqrt(b2(ixO^S))*inv_rho*kmax)
!    end if
!
!  end subroutine srmhd_get_csound

  !> Calculate fast magnetosonic wave speed
  subroutine srmhd_get_csound_prim(w,x,ixI^L,ixO^L,idim,csound)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)
    double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
    double precision :: inv_rho(ixO^S)

    inv_rho=1.d0/w(ixO^S,rho_)

    if(srmhd_energy) then
      csound(ixO^S)=dsqrt(srmhd_gamma*w(ixO^S,p_)*inv_rho)
    else
      csound(ixO^S)=dsqrt(srmhd_gamma*srmhd_adiab*w(ixO^S,rho_)**gamma_1)
    end if
!    ! store |B|^2 in v
!    b2(ixO^S)        = srmhd_mag_en_all(w,ixI^L,ixO^L)
!    cfast2(ixO^S)   = b2(ixO^S) * inv_rho+csound(ixO^S)
!    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
!         * srmhd_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
!         * inv_rho

!    where(AvMinCs2(ixO^S)<zero)
!       AvMinCs2(ixO^S)=zero
!    end where

!    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))
!
!    if (.not. MHD_Hall) then
!       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
!    else
!       ! take the Hall velocity into account:
!       ! most simple estimate, high k limit:
!       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
!       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
!       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
!            srmhd_etah * sqrt(b2(ixO^S))*inv_rho*kmax)
!    end if
!
  end subroutine srmhd_get_csound_prim


 !> Calculate magnetic pressure within ixO^L including magnetic pressure
  subroutine srmhd_get_p_mag(w,ixI^L,ixO^L,pmag)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(out)   :: pmag(ixI^S)

    pmag(ixO^S) = 0.5d0 * (sum(w(ixO^S, mag(:))**2, dim=ndim+1)&
            +sum(w(ixO^S, mag(:))*w(ixO^S, mom(:)), dim=ndim+1)&
               /w(ixO^S,lfac_))&
               / w(ixO^S,lfac_)

  end subroutine srmhd_get_p_mag

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine srmhd_get_csound2(w,x,ixI^L,ixO^L,rhoh,csound2)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw),rhoh(ixO^S)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixO^S)
  
    double precision                :: rho(ixO^S)

    if(srmhd_energy) then
      rho=w(ixO^S,d_)/w(ixO^S,lfac_)
      call srmhd_get_csound2_eos(ixI^L,ixO^L,x,rho,rhoh,csound2)
    else
      csound2(ixO^S)=srmhd_gamma*srmhd_adiab*rho(ixO^S)**gamma_1
    end if
  end subroutine srmhd_get_csound2



  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine srmhd_get_p_total(w,x,ixI^L,ixO^L,p,B2,VdotB)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    implicit none

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: p(ixI^S)
    double precision, optional, intent(in) :: B2(ixO^S),VdotB(ixO^S)

    double precision                :: pmag(ixI^S)
    call srmhd_get_pthermal(w,x,ixI^L,ixO^L,p)
    call srmhd_get_pmag(w,x,ixI^L,ixO^L,pmag,B2=B2,VdotB=VdotB)
    p(ixO^S) = p(ixO^S) + pmag(ixO^S)
  end subroutine srmhd_get_p_total

  subroutine srmhd_get_pmag(w,x,ixI^L,ixO^L,pmag,B2,VdotB)
    ! made by Z. Meliani 10/02/2018 
   use mod_global_parameters
   implicit none
    integer, intent(in)                    :: ixI^L, ixO^L
    double precision, intent(in)           :: w(ixI^S,nw)
    double precision, intent(in)           :: x(ixI^S,1:ndim)
    double precision, intent(out)          :: pmag(ixI^S)
    double precision, optional, intent(in) :: B2(ixO^S),VdotB(ixO^S)

    double precision, dimension(ixO^S)     :: sub_B2,sub_VdotB

    is_presentB2 : if(present(B2)) then
     pmag(ixO^S) = 0.5d0*(VdotB(ixO^S)**2.0d0  &
                       + B2(ixO^S)/w(ixO^S,lfac_)**2.0d0)
    else  is_presentB2
     call srmhd_get_B2andVdotB(ixI^L,ixO^L,w,.true.,B2=sub_B2,VdotB=sub_VdotB)
     pmag(ixO^S) = 0.5d0*(sub_VdotB(ixO^S)**2.0d0  &
                       + sub_B2(ixO^S)/w(ixO^S,lfac_)**2.0d0)
    end if is_presentB2
  end subroutine srmhd_get_pmag

  !> Calculate B2 and/ot VdotB in ixO^S from conserve or primitive variables
  subroutine srmhd_get_B2andVdotB(ixI^L,ixO^L,w,wconserve,B2,VdotB)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters
    implicit none
    integer, intent(in)                      :: ixI^L, ixO^L
    ! conservative or primitive w
    double precision, intent(in)             :: w(ixI^S,nw) 
    logical, intent(in)                      :: wconserve  
    double precision,optional, intent(inout) :: B2(ixO^S),VdotB(ixO^S)

    is_VdotB_present : if(present(VdotB))then
     is_conserve : if(wconserve)then
       VdotB = sum(w(ixO^S,mag(:))*w(ixO^S,mom(:)),dim=ndim+1)&
               /w(ixO^S,xi_)
     else is_conserve
       VdotB = sum(w(ixO^S,mag(:))*w(ixO^S,mom(:)),dim=ndim+1)&
               /w(ixO^S,lfac_)
     end if is_conserve
    end if is_VdotB_present
    if(present(B2))B2    = sum(w(ixO^S,mag(:))**2.0,dim=ndim+1)  

  end subroutine srmhd_get_B2andVdotB

  
  !> Calculate fluxes within ixO^L.
  subroutine srmhd_get_flux(wC,wP,x,ixI^L,ixO^L,idim,f)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    implicit none
    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: wP(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)

    double precision             :: ptotal(ixI^S)
    double precision             :: B2(ixO^S),VdotB(ixO^S)
    double precision             :: v(ixI^S,1:ndir)
    integer                      :: idirmin, iw, idir


    call srmhd_get_B2andVdotB(ixI^L,ixO^L,wp,.false.,B2=B2,VdotB=VdotB)
    call srmhd_get_p_total(wC,x,ixI^L,ixO^L,ptotal,B2=B2,VdotB=VdotB)
    call srmhd_get_v(wC,x,ixI^L,ixO^L,v,B2=B2,VdotB=VdotB)
    ! Get flux of density
    f(ixO^S,rho_)=wP(ixO^S,mom(idim))*wP(ixO^S,rho_)

    ! Get flux of tracer
    do iw=1,srmhd_n_tracer
      f(ixO^S,tracer(iw))=wP(ixO^S,mom(idim))*wP(ixO^S,tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-B_k*b_i [+ptotal if i==k]
    Loop_idir_mom : do idir=1,ndir
      f(ixO^S,mom(idir))= v(ixO^S,idir)*wC(ixO^S,mom(idim))&
                         -wP(ixO^S,mag(idir))*(VdotB*v(ixO^S,idim)&
                         +wP(ixO^S,mag(idim))/wP(ixO^S,lfac_)**2.0d0)

    end do Loop_idir_mom
    f(ixO^S,mom(idim))=ptotal(ixO^S)+f(ixO^S,mom(idim))

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    is_notiso : if (srmhd_energy) then
       is_internal : if (block%e_is_internal) then
          f(ixO^S,e_)=v(ixO^S,idim)*wP(ixO^S,p_)
       else is_internal
          f(ixO^S,e_)=v(ixO^S,idim)*(wC(ixO^S,e_) + ptotal(ixO^S))- &
             wP(ixO^S,mag(idim))*VdotB(ixO^S)
          if(type_divb==divb_glm2) then
            f(ixO^S,e_) = f(ixO^S,e_) + vmax_global*wP(ixO^S,psi_)&
                                         *wP(ixO^S,mag(idim))
          end if

       end if is_internal
    end if is_notiso

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    Loop_idir_Bflux : do idir=1,ndir
      is_idim_Bflux : if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        is_glm_Bflux : if (srmhd_glm) then
           if(type_divb==divb_glm1) then
             f(ixO^S,mag(idir))=wP(ixO^S,psi_)
           else
             f(ixO^S,mag(idir))=vmax_global*wP(ixO^S,psi_)
           end if
        else is_glm_Bflux
           f(ixO^S,mag(idir))=zero
        end if is_glm_Bflux
      else is_idim_Bflux
        f(ixO^S,mag(idir))=v(ixO^S,idim)*wP(ixO^S,mag(idir))&
                           -wP(ixO^S,mag(idim))*v(ixO^S,idir)

      end if is_idim_Bflux
    end do Loop_idir_Bflux

    if (srmhd_glm) then
      if(type_divb==divb_glm1) then
        !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
        f(ixO^S,psi_)  = cmax_global**2*wP(ixO^S,mag(idim))
      else
        !f_i[psi]=Ch*b_{i} Eq. 3.16e Derigs et al 2018 JCP, 364, 420 
        f(ixO^S,psi_)  = vmax_global*wP(ixO^S,mag(idim))
      end if
    end if

  end subroutine srmhd_get_flux

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine srmhd_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    ! made by Z. Meliani 10/02/2018 
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
      if (srmhd_energy .and. block%e_is_internal) then
        active = .true.
        call internal_energy_add_source(qdt,ixI^L,ixO^L,wCT,w,x)
      endif

      ! Source for B0 splitting
!      if (B0field) then
!        active = .true.
!        call add_source_B0split(qdt,ixI^L,ixO^L,wCT,w,x)
!      end if

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
!      if (abs(srmhd_eta)>smalldouble)then
!        active = .true.
!        call add_source_res2(qdt,ixI^L,ixO^L,wCT,w,x)
!      end if

!      if (srmhd_eta_hyper>0.d0)then
!        active = .true.
!        call add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
!      end if
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
        call add_source_glm2(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
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
        call add_source_glm2(dt,ixI^L,ixO^L,pw(saveigrid)%wold,w,x)
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
        call add_source_glm2(qdt,ixI^L,ixO^L,wCT,w,x)
      case default
        call mpistop('Unknown divB fix')
      end select
    end if
    }

!    if(srmhd_radiative_cooling) then
!      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,&
!           w,x,qsourcesplit,active)
!    end if

!    if(srmhd_viscosity) then
!      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,&
!           w,x,srmhd_energy,qsourcesplit,active)
!    end if

    if(srmhd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,&
           w,x,srmhd_energy,qsourcesplit,active)
    end if
    call srmhd_get_auxiliary(.true.,w,x,ixI^L,ixO^L,'srmhd_add_source')

  end subroutine srmhd_add_source

  subroutine internal_energy_add_source(qdt,ixI^L,ixO^L,wCT,w,x)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: pth(ixI^S),v(ixI^S,1:ndir),divv(ixI^S)

    call srmhd_get_v(wCT,x,ixI^L,ixI^L,v)
    call divvector(v,ixI^L,ixO^L,divv)
    call srmhd_get_pthermal(wCT,x,ixI^L,ixO^L,pth)
    w(ixO^S,e_)=w(ixO^S,e_)-qdt*pth(ixO^S)*divv(ixO^S)

  end subroutine internal_energy_add_source



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
    call get_divb(wCT,ixI^L,ixO^L,divb)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (srmhd_glm_alpha < zero) then
      w(ixO^S,psi_) = abs(srmhd_glm_alpha)*wCT(ixO^S,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab) then
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*srmhd_glm_alpha/minval(dxlevel(:)))*w(ixO^S,psi_)
      else
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*srmhd_glm_alpha/minval(block%ds(ixO^S,:),dim=ndim+1))*w(ixO^S,psi_)
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
       if (srmhd_energy .and. .not.block%e_is_internal) then
       ! e  = e  -qdt (b . grad(Psi))
         w(ixO^S,e_) = w(ixO^S,e_)-qdt*wCT(ixO^S,mag(idim))*gradPsi(ixO^S)
       end if
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))&
                        -qdt*srmhd_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    if (check_small_values) call srmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_glm1')

  end subroutine add_source_glm1

  subroutine add_source_glm2(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Eq. 3.17 Derigs et al 2018 JCP, 364, 420
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision:: divb(ixI^S),v(ixI^S,1:ndir)
    integer          :: idim,idir
    double precision :: gradPsi(ixI^S,ndim), Bf(ixI^S,1:ndir)

    ! calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb)

    ! calculate velocity
    call srmhd_get_v(wCT,x,ixI^L,ixO^L,v)

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
         call gradient(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi(ixI^S,idim))
       case("limited")
         call gradientS(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi(ixI^S,idim))
       end select
    end do

    if(B0field) then
      Bf(ixI^S,1:ndir)=wCT(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,0)
    else
      Bf(ixI^S,1:ndir)=wCT(ixI^S,mag(1:ndir))
    end if

    if (srmhd_energy .and. .not.block%e_is_internal) then
       ! e = e - qdt ( (v . b) * div b + (grad psi . v) * psi)
       w(ixO^S,e_)=w(ixO^S,e_) - qdt * (divb(ixO^S) * &
            sum(v(ixO^S,:)*Bf(ixO^S,:),dim=ndim+1) + wCT(ixO^S,psi_)*&
            sum(v(ixO^S,1:ndim)*gradPsi(ixO^S,1:ndim),dim=ndim+1))
    end if

    ! b_i = b_i - qdt * v_i * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*v(ixO^S,idir)*divb(ixO^S)
    end do

    ! m_i = m_i - qdt * b_i * div b
    do idir=1,ndir
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-qdt*Bf(ixO^S,idir)*divb(ixO^S)
    end do

    ! psi = psi - qdt * (v . grad(psi) + alpha * psi)
    w(ixO^S,psi_) = w(ixO^S,psi_)-qdt*(sum(v(ixO^S,1:ndim)*gradPsi(ixO^S,1:ndim),dim=ndim+1)+&
      vmax_global/0.18d0*wCT(ixO^S,psi_))

    if (check_small_values) call srmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_glm')

  end subroutine add_source_glm2

  !> Add divB related sources to w within ixO corresponding to Powel
  subroutine add_source_powel(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divb(ixI^S),v(ixI^S,1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb)

    ! calculate velocity
    call srmhd_get_v(wCT,x,ixI^L,ixO^L,v)

    if (srmhd_energy .and. .not.block%e_is_internal) then
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
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-qdt*srmhd_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    if (check_small_values) call srmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_powel')

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
    call get_divb(wCT,ixI^L,ixO^L,divb)

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*wCT(ixO^S,mom(idir))/wCT(ixO^S,rho_)*divb(ixO^S)
    end do

    if (check_small_values) call srmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_janhunen')

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
    call get_divb(wCT,ixI^L,ixp^L,divb)

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

       if (srmhd_energy .and. typedivbdiff=='all' .and. .not.block%e_is_internal) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixp^S,e_)=w(ixp^S,e_)+wCT(ixp^S,mag(idim))*graddivb(ixp^S)
       end if
    end do

    if (check_small_values) call srmhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_linde')

  end subroutine add_source_linde

  !> Calculate div B within ixO
  subroutine get_divb(w,ixI^L,ixO^L,divb)

    use mod_global_parameters
    use mod_geometry

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision                   :: divb(ixI^S)

    double precision                   :: bvec(ixI^S,1:ndir)

    bvec(ixI^S,:)=w(ixI^S,mag(:))

    select case(typediv)
    case("central")
      call divvector(bvec,ixI^L,ixO^L,divb)
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
      divb(ixO^S)=0.5d0*abs(divb(ixO^S))/sqrt(srmhd_mag_en_all(w,ixI^L,ixO^L))/sum(1.d0/dxlevel(:))
    else
      ixAmin^D=ixOmin^D-1;
      ixAmax^D=ixOmax^D-1;
      dsurface(ixO^S)= sum(block%surfaceC(ixO^S,:),dim=ndim+1)
      do idims=1,ndim
        ixA^L=ixO^L-kr(idims,^D);
        dsurface(ixO^S)=dsurface(ixO^S)+block%surfaceC(ixA^S,idims)
      end do
      divb(ixO^S)=abs(divb(ixO^S))/sqrt(srmhd_mag_en_all(w,ixI^L,ixO^L))*&
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
  subroutine srmhd_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
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


    if(srmhd_gravity) then
      call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

  end subroutine srmhd_get_dt

  ! Add geometrical source terms to w
  subroutine srmhd_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision :: tmp(ixI^S),ptot(ixI^S),tmp2(ixI^S)
    double precision :: B2(ixO^S),VdotB(ixO^S),VdotB2(ixO^S),xiplusB2(ixO^S)
    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_
    character(len=30)  :: subname_loc

    if(typeaxial /= 'slab')then
     subname_loc='srmhd_add_geom'
     ! get auxiliary variables
     call srmhd_get_auxiliary(.true.,wCT,x,ixI^L,ixO^L,subname_loc)
     call srmhd_get_B2andVdotB(ixI^L,ixO^L,wCT,.true.,B2=B2,VdotB=VdotB)
     call srmhd_get_p_total(wCT,x,ixI^L,ixO^L,ptot)
     VdotB2=VdotB**2.0
     xiplusB2(ixO^S) = wCT(ixO^S,xi_)**2.0+B2(ixO^S)
    
     mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
     br_=mag(1); bphi_=mag(1)-1+phi_
    end if

    select case (typeaxial)
    case ('cylindrical')
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in MHD")
      endif
       if(phi_>0) then
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*(ptot(ixO^S)+&
                   (wCT(ixO^S,mphi_)**2-VdotB2(ixO^S)*wCT(ixO^S,bphi_)**2)&
                   /xiplusB2(ixO^S)-wCT(ixO^S,bphi_)**2/wCT(ixO^S,lfac_)**2.0)

         w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt/x(ixO^S,1)*(&
                  -(wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)&
                  -VdotB2(ixO^S)*wCT(ixO^S,br_)*wCT(ixO^S,bphi_))&
                  /xiplusB2(ixO^S) &
                  +wCT(ixO^S,bphi_)*wCT(ixO^S,br_)/wCT(ixO^S,lfac_)**2)

         w(ixO^S,bphi_)=w(ixO^S,bphi_)+qdt/x(ixO^S,1)*&
                  (wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
                  -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
                  /xiplusB2(ixO^S)
       else
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*ptot(ixO^S)
       end if
       if(srmhd_glm) w(ixO^S,br_)=w(ixO^S,br_)+qdt*wCT(ixO^S,psi_)/x(ixO^S,1)
    case ('spherical')
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       
       ! m1
       tmp(ixO^S)=ptot(ixO^S)*x(ixO^S,1) &
         *(block%surfaceC(ixO^S,1)-block%surfaceC(h1x^S,1))/block%dvolume(ixO^S)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(idir))**2.0&
                         -(VdotB2(ixO^S)*B2(ixO^S)))&
                        /(xiplusB2(ixO^S))&
                      -wCT(ixO^S,mag(idir))**2/wCT(ixO^S,lfac_)**2.0
         end do
       end if
       w(ixO^S,mom(1))=w(ixO^S,mom(1))+qdt*tmp(ixO^S)/x(ixO^S,1)
       ! b1
       if(srmhd_glm) then
         w(ixO^S,mag(1))=w(ixO^S,mag(1))+qdt/x(ixO^S,1)*2.0d0*wCT(ixO^S,psi_)
       end if

       {^NOONED
       ! m2
       !tmp(ixO^S)=tmp1(ixO^S)
       ! This will make hydrostatic p=const an exact solution
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*ptot(ixO^S) &
            *(block%surfaceC(ixO^S,2)-block%surfaceC(h2x^S,2)) &
            /block%dvolume(ixO^S)

       tmp(ixO^S)=-((wCT(ixO^S,mom(1))*wCT(ixO^S,mom(2))&
                   -VdotB2(ixO^S)*wCT(ixO^S,mag(1))*wCT(ixO^S,mag(2)))&
                    /xiplusB2(ixO^S) &
            -wCT(ixO^S,mag(1))*wCT(ixO^S,mag(2))/wCT(ixO^S,lfac_)**2.0)

       if(ndir==3) then
         tmp(ixO^S)=tmp(ixO^S)+((wCT(ixO^S,mom(3))**2&
              -VdotB2(ixO^S)*wCT(ixO^S,mag(3))**2.0)/xiplusB2(ixO^S) &
              -(wCT(ixO^S,mag(3))/wCT(ixO^S,lfac_))**2.0)&
              *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
       end if
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
       ! b2
       tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(2)) &
            -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(1)))/xiplusB2(ixO^S)

       if(srmhd_glm) then
         tmp(ixO^S)=tmp(ixO^S) &
              + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
       end if
       w(ixO^S,mag(2))=w(ixO^S,mag(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
       }

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixO^S)=-((wCT(ixO^S,mom(3))*wCT(ixO^S,mom(1))&
                       -VdotB2(ixO^S)*wCT(ixO^S,mag(3))*wCT(ixO^S,mag(1)))&
                       /xiplusB2(ixO^S) &
                -wCT(ixO^S,mag(3))*wCT(ixO^S,mag(1))&
                /wCT(ixO^S,lfac_)**2.0d0) {^NOONED &
                -((wCT(ixO^S,mom(2))*wCT(ixO^S,mom(3))&
                -VdotB2(ixO^S)*wCT(ixO^S,mag(2))*wCT(ixO^S,mag(3)))&
                /xiplusB2(ixO^S) &
                -wCT(ixO^S,mag(2))*wCT(ixO^S,mag(3))/wCT(ixO^S,lfac_)**2.0d0) &
                *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
           w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
         else
           call mpistop("angmomfix not implemented yet in MHD")
         end if
         ! b3
         tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(3)) &
              -wCT(ixO^S,mom(3))*wCT(ixO^S,mag(1)))/xiplusB2(ixO^S) {^NOONED &
              -(wCT(ixO^S,mom(3))*wCT(ixO^S,mag(2)) &
              -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(3)))*dcos(x(ixO^S,2)) &
              /(xiplusB2(ixO^S)*dsin(x(ixO^S,2))) }
         w(ixO^S,mag(3))=w(ixO^S,mag(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
       end if
    end select
  end subroutine srmhd_add_source_geom

  !> Compute 2 times total magnetic energy
  function srmhd_mag_en_all(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    if (B0field) then
      mge = sum((w(ixO^S, mag(:))+block%B0(ixO^S,:,block%iw0))**2, dim=ndim+1)
    else
      mge = sum(w(ixO^S, mag(:))**2, dim=ndim+1)
    end if
  end function srmhd_mag_en_all

  !> Compute full magnetic field by direction
  function srmhd_mag_i_all(w, ixI^L, ixO^L,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idir
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mgf(ixO^S)

    if (B0field) then
      mgf = w(ixO^S, mag(idir))+block%B0(ixO^S,idir,block%iw0)
    else
      mgf = w(ixO^S, mag(idir))
    end if
  end function srmhd_mag_i_all

  !> Compute evolving magnetic energy from primitive
  function srmhd_mag_en_primitive(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    mge = sum(w(ixO^S, mag(:))**2_dp, dim=ndim+1)
    mge = 0.5_dp*(mge(ixO^S)+ mge(ixO^S)*sum(w(ixO^S, mom(:))**2_dp&
           , dim=ndim+1)/w(ixO^S,lfac_)&
          +sum(w(ixO^S, mag(:))*w(ixO^S, mom(:))))/w(ixO^S,lfac_)
  end function srmhd_mag_en_primitive

  !> compute kinetic energy from primitive
  function srmhd_kin_en_primitive(w, ixI^L, ixO^L, inv_rho) result(ke)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)
    double precision, intent(in), optional :: inv_rho(ixO^S)

    if (present(inv_rho)) then
       ke = (w(ixO^S,xi_)-1.0d0) * inv_rho
    else
       ke = (w(ixO^S,xi_)-1.0d0) / w(ixO^S, rho_)
    end if
  end function srmhd_kin_en_primitive

!  !>get Hall speed
!  subroutine srmhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
!    use mod_global_parameters

!    integer, intent(in)             :: ixI^L, ixO^L
!    double precision, intent(in)    :: w(ixI^S,nw)
!    double precision, intent(in)    :: x(ixI^S,1:ndim)
!    double precision, intent(inout) :: vHall(ixI^S,1:3)

!    integer          :: idir, idirmin
!    double precision :: current(ixI^S,7-2*ndir:3)

!    ! Calculate current density and idirmin
!    call get_current(w,ixI^L,ixO^L,idirmin,current)
!    vHall(ixO^S,1:3) = zero
!    vHall(ixO^S,idirmin:3) = - srmhd_etah*current(ixO^S,idirmin:3)
!    do idir = idirmin, 3
!       vHall(ixO^S,idir) = vHall(ixO^S,idir)/w(ixO^S,rho_)
!    end do

!  end subroutine srmhd_getv_Hall

!  subroutine srmhd_getdt_Hall(w,x,ixI^L,ixO^L,dx^D,dthall)
!    use mod_global_parameters

!    integer, intent(in) :: ixI^L, ixO^L
!    double precision, intent(in)    :: dx^D
!    double precision, intent(in)    :: w(ixI^S,1:nw)
!    double precision, intent(in)    :: x(ixI^S,1:ndim)
!    double precision, intent(out)   :: dthall
!    !.. local ..
!    double precision :: dxarr(ndim)
!    double precision :: bmag(ixI^S)

!    dthall=bigdouble

!    ! because we have that in cmax now:
!    return

!    ^D&dxarr(^D)=dx^D;

!    if (.not. B0field) then
!       bmag(ixO^S)=sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1))
!       bmag(ixO^S)=sqrt(sum((w(ixO^S,mag(:)) + block%B0(ixO^S,1:ndir,block%iw0))**2))
!    end if

!    if(slab) then
!      dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(srmhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
!    else
!      dthall=dtdiffpar*minval(block%ds(ixO^S,1:ndim))**2.0d0/(srmhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
!    end if

!  end subroutine srmhd_getdt_Hall

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

  subroutine srmhd_boundary_adjust
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

  end subroutine srmhd_boundary_adjust

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
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_conserved(ixG^L,ixO^L,w,x)
     case(2)
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_conserved(ixG^L,ixO^L,w,x)
     case(3)
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_conserved(ixG^L,ixO^L,w,x)
     case(4)
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_conserved(ixG^L,ixO^L,w,x)
     {^IFTHREED
     case(5)
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_conserved(ixG^L,ixO^L,w,x)
     case(6)
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_primitive(ixG^L,ixO^L,w,x)
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
       if(srmhd_energy.and..not.block%e_is_internal) call srmhd_to_conserved(ixG^L,ixO^L,w,x)
     }
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

   end subroutine fixdivB_boundary
   
   !> subtroutine srmhd_get_auxiliary calcule using srmhd_con2prim to calculate the enthalpy and the lorentz factor
   subroutine srmhd_get_auxiliary(clipping,w,x,ixI^L,ixO^L,subname)
    ! made by Z. Meliani 10/02/2018 
    use mod_global_parameters
    use mod_srmhd_con2prim
    implicit none

    logical, intent(in)             :: clipping
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    character(len=*), intent(in)    :: subname
  
    integer                             :: ix^D,ierror
    integer                             :: flag_error(ixO^S)
    character(len=len(subname)+30):: subname_loc
  ! print*,', is yooou ',subname,saveigrid
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
    call srmhd_con2prim(w(ix^D,d_),w(ix^D,mom(1):mom(ndir)),w(ix^D,e_),&
             w(ix^D,mag(1):mag(ndir)),w(ix^D,lfac_),w(ix^D,xi_),ierror)
    if(check_small_values)then
     if(ierror/=0) then
       flag_error(ix^D) = ierror
     else
       flag_error(ix^D) = 0
     end if 
     
    end if
   {enddo^D&\} 
   is_check_small : if(check_small_values)then
    is_flag_on : if(any(flag_error(ixO^S)/=0))then
     subname_loc='srmhd_get_auxiliary from -> '//trim(subname)
     call srmhd_handle_small_values(.false., &
                                    w, x, &
                                    ixI^L, ixO^L,subname_loc,&
                                    flag_error=flag_error)
    end if is_flag_on
   end if is_check_small
   end subroutine srmhd_get_auxiliary

   subroutine srmhd_get_4u_from_3v(ixI^L,ixO^L,vtou)
    ! made by Z. Meliani 13/02/2018 
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: vtou(ixI^S,1:ndir)

    double precision                :: lfac(ixO^S)
    integer                         :: idir


    lfac= 1.0d0/dsqrt(1.0d0-sum(vtou(ixO^S,1:ndir)**2.0,dim=ndim+1))
    Loop_idir: do idir=1,ndir
     vtou(ixO^S,idir)=lfac(ixO^S)*vtou(ixO^S,idir)
    end do Loop_idir  
   end subroutine srmhd_get_4u_from_3v
end module mod_srmhd_phys
