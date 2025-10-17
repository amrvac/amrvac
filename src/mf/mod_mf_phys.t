!> Magnetofriction module
module mod_mf_phys
  use mod_global_parameters, only: std_len
  use mod_functions_bfield, only: get_divb,mag
  use mod_comm_lib, only: mpistop


  implicit none
  private

  !> viscosity coefficient in s cm^-2 for solar corona (Cheung 2012 ApJ)
  double precision, public                :: mf_nu = 1.d-15

  !> maximal limit of magnetofrictional velocity in cm s^-1 (Pomoell 2019 A&A)
  double precision, public                :: mf_vmax = 3.d6

  !> decay scale of frictional velocity near boundaries
  double precision, public                :: mf_decay_scale(2*^ND)=0.d0

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: mf_glm_alpha = 0.5d0

  !> The resistivity
  double precision, public                :: mf_eta = 0.0d0

  !> The hyper-resistivity
  double precision, public                :: mf_eta_hyper = 0.0d0

  !> Coefficient of diffusive divB cleaning
  double precision :: divbdiff     = 0.8d0

  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.1d0

  !> Method type in a integer for good performance
  integer :: type_divb

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Indices of the momentum density for the form of better vectorization
  integer, public, protected              :: ^C&m^C_

  !> Indices of the magnetic field for the form of better vectorization
  integer, public, protected              :: ^C&b^C_

  !> Indices of the GLM psi
  integer, public, protected :: psi_

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

  !> Whether particles module is added
  logical, public, protected              :: mf_particles = .false.

  !> Whether GLM-MHD is used
  logical, public, protected              :: mf_glm = .false.

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.

  !> MHD fourth order
  logical, public, protected              :: mf_4th_order = .false.

  !> set to true if need to record electric field on cell edges
  logical, public, protected              :: mf_record_electric_field = .false.

  !> Whether divB is computed with a fourth order approximation
  integer, public, protected :: mf_divb_nth = 1

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*^ND)=.true.

  !> Use a compact way to add resistivity
  logical :: compactres   = .false.

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> clean divb in the initial condition
  logical, public, protected :: clean_initial_divb=.false.

  !> Method type to clean divergence of B
  character(len=std_len), public, protected :: typedivbfix  = 'ct'

  !> Method type of constrained transport
  character(len=std_len), public, protected :: type_ct  = 'average'

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'


  ! Public methods
  public :: mf_phys_init
  public :: mf_to_conserved
  public :: mf_to_primitive
  public :: mf_face_to_center
  public :: get_divb
  public :: get_current
  public :: get_normalized_divb
  public :: b_from_vector_potential
  public :: record_force_free_metrics
  {^NOONED
  public :: mf_clean_divb_multigrid
  }

contains

  !> Read this module's parameters from a file
  subroutine mf_read_params(files)
    use mod_global_parameters
    use mod_particles, only: particles_eta
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mf_list/ mf_nu, mf_vmax, mf_decay_scale, &
      mf_eta, mf_eta_hyper, mf_glm_alpha, mf_particles,&
      particles_eta, mf_record_electric_field,&
      mf_4th_order, typedivbfix, source_split_divb, divbdiff,&
      typedivbdiff, type_ct, compactres, divbwave, He_abundance, SI_unit, &
      Bdip, Bquad, Boct, Busr, clean_initial_divb, &
      boundary_divbfix, boundary_divbfix_skip, mf_divb_nth

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mf_list, end=111)
111    close(unitpar)
    end do
 
  end subroutine mf_read_params

  !> Write this module's parameters to a snapsoht
  subroutine mf_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "nu"
    values(1) = mf_nu
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine mf_write_info

  subroutine mf_phys_init()
    use mod_global_parameters
    use mod_physics
    use mod_particles, only: particles_init, particles_eta, particles_etah
    {^NOONED
    use mod_multigrid_coupling
    }

    integer :: itr, idir

    call mf_read_params(par_files)

    physics_type = "mf"
    phys_energy = .false.
    use_particles=mf_particles

    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
       type_divb = divb_none
    {^NOONED
    case ('multigrid')
       type_divb = divb_multigrid
       use_multigrid = .true.
       mg%operator_type = mg_laplacian
       phys_global_source_after => mf_clean_divb_multigrid
    }
    case ('glm')
      mf_glm          = .true.
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
      mf_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_lindeglm
    case ('ct')
      type_divb = divb_ct
      stagger_grid = .true.
    case default
      call mpistop('Unknown divB fix')
    end select

    ! set velocity field as flux variables
    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)
    m^C_=mom(^C);

    ! set magnetic field as flux variables
    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)
    b^C_=mag(^C);

    ! start with magnetic field and skip velocity when update ghostcells
    iwstart=mag(1)
    allocate(start_indices(number_species),stop_indices(number_species))
    ! set the index of the first flux variable for species 1
    start_indices(1)=mag(1)

    if (mf_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    ! set number of variables which need update ghostcells
    nwgc=nwflux-ndir

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    ! determine number of stagger variables
    nws=ndim

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
    if(mf_glm .and. ndim>1) flux_type(:,psi_)=flux_tvdlf

    phys_get_dt              => mf_get_dt
    phys_get_cmax            => mf_get_cmax
    phys_get_cbounds         => mf_get_cbounds
    phys_get_flux            => mf_get_flux
    phys_add_source_geom     => mf_add_source_geom
    phys_add_source          => mf_add_source
    phys_to_conserved        => mf_to_conserved
    phys_to_primitive        => mf_to_primitive
    phys_check_params        => mf_check_params
    phys_write_info          => mf_write_info
    phys_special_advance     => mf_velocity_update

    if(type_divb==divb_glm) then
      phys_modify_wLR => mf_modify_wLR
    end if

    ! pass to global variable to record electric field
    record_electric_field=mf_record_electric_field

    ! Initialize particles module
    if(mf_particles) then
      call particles_init()
      ! never allow Hall effects in particles when doing magnetofrictional
      particles_etah=0.0d0
      if(mype==0) then
         write(*,*) '*****Using particles:     with mf_eta, mf_eta_hyper :', mf_eta, mf_eta_hyper
         write(*,*) '*****Using particles: particles_eta, particles_etah :', particles_eta, particles_etah
      end if
    end if

    ! if using ct stagger grid, boundary divb=0 is not done here
    if(stagger_grid) then
      phys_face_to_center => mf_face_to_center
      select case(type_ct)
      case('average')
        transverse_ghost_cells = 1
        phys_get_ct_velocity => mf_get_ct_velocity_average
        phys_update_faces => mf_update_faces_average
      case('uct_contact')
        transverse_ghost_cells = 1
        phys_get_ct_velocity => mf_get_ct_velocity_contact
        phys_update_faces => mf_update_faces_contact
      case('uct_hll')
        transverse_ghost_cells = 2
        phys_get_ct_velocity => mf_get_ct_velocity_hll
        phys_update_faces => mf_update_faces_hll
      case default
        call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
      end select
    else if(ndim>1) then
      phys_boundary_adjust => mf_boundary_adjust
    end if

    {^NOONED
    ! clean initial divb
    if(clean_initial_divb) phys_clean_divb => mf_clean_divb_multigrid
    }

    ! derive units from basic units
    call mf_physical_units()

  end subroutine mf_phys_init

  subroutine mf_check_params
    use mod_global_parameters

  end subroutine mf_check_params

  subroutine mf_physical_units()
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
      miu0=4.d0*dpi
      c_lightspeed=const_c
    end if
    a=1.d0+4.d0*He_abundance
    b=2.d0+3.d0*He_abundance
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

    ! get dimensionless mf nu
    mf_nu=mf_nu/unit_time*unit_length**2

    ! get dimensionless maximal mf velocity limit
    mf_vmax=mf_vmax/unit_velocity
 
  end subroutine mf_physical_units

  !> Transform primitive variables into conservative ones
  subroutine mf_to_conserved(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    ! nothing to do for mf
  end subroutine mf_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine mf_to_primitive(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    ! nothing to do for mf
  end subroutine mf_to_primitive

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine mf_get_cmax(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)

    cmax(ixO^S)=abs(w(ixO^S,mom(idim)))+one

  end subroutine mf_get_cmax

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine mf_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,Hspeed,cmax,cmin)
    use mod_global_parameters
    use mod_variables

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
    double precision, intent(in)    :: Hspeed(ixI^S,1:number_species)

    double precision, dimension(ixI^S) :: tmp1

    tmp1(ixO^S)=0.5d0*(wLC(ixO^S,mom(idim))+wRC(ixO^S,mom(idim)))
    if(present(cmin)) then
      cmax(ixO^S,1)=max(tmp1(ixO^S)+one,zero)
      cmin(ixO^S,1)=min(tmp1(ixO^S)-one,zero)
    else
      cmax(ixO^S,1)=abs(tmp1(ixO^S))+one
    end if

  end subroutine mf_get_cbounds

  !> prepare velocities for ct methods
  subroutine mf_get_ct_velocity_average(vcts,wLp,wRp,ixI^L,ixO^L,idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: cmax(ixI^S)
    double precision, intent(in), optional :: cmin(ixI^S)
    type(ct_velocity), intent(inout):: vcts

  end subroutine mf_get_ct_velocity_average

  subroutine mf_get_ct_velocity_contact(vcts,wLp,wRp,ixI^L,ixO^L,idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: cmax(ixI^S)
    double precision, intent(in), optional :: cmin(ixI^S)
    type(ct_velocity), intent(inout):: vcts

    ! calculate velocities related to different UCT schemes
    if(.not.allocated(vcts%vnorm)) allocate(vcts%vnorm(ixI^S,1:ndim))
    ! get average normal velocity at cell faces
    vcts%vnorm(ixO^S,idim)=0.5d0*(wLp(ixO^S,mom(idim))+wRp(ixO^S,mom(idim)))

  end subroutine mf_get_ct_velocity_contact

  subroutine mf_get_ct_velocity_hll(vcts,wLp,wRp,ixI^L,ixO^L,idim,cmax,cmin)
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
         +vcts%cbarmin(ixO^S,idim)*vcts%vbarRC(ixO^S,idim,1))&
        /(vcts%cbarmax(ixO^S,idim)+vcts%cbarmin(ixO^S,idim))

  end subroutine mf_get_ct_velocity_hll

  !> Calculate fluxes within ixO^L.
  subroutine mf_get_flux(wC,w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)

    integer :: ix^D

   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! compute flux of magnetic field
      ! f_i[b_k]=v_i*b_k-v_k*b_i
      ^C&f(ix^D,b^C_)=w(ix^D,mom(idim))*w(ix^D,b^C_)-w(ix^D,mag(idim))*w(ix^D,m^C_)\
   {end do\}
    if(mf_glm) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,mag(idim))=w(ix^D,psi_)
        !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
        f(ix^D,psi_)=cmax_global**2*w(ix^D,mag(idim))
     {end do\}
    end if

  end subroutine mf_get_flux

  !> Add global source terms to update frictional velocity on complete domain
  subroutine mf_velocity_update(qt,psa)
    use mod_global_parameters
    use mod_ghostcells_update
    double precision, intent(in) :: qt     !< Current time
    type(state), target :: psa(max_blocks) !< Compute based on this state

    integer :: iigrid,igrid
    logical :: stagger_flag
    logical :: firstmf=.true.

    if(firstmf) then
      ! point bc mpi datatype to partial type for velocity field
      type_send_srl=>type_send_srl_p1
      type_recv_srl=>type_recv_srl_p1
      type_send_r=>type_send_r_p1
      type_recv_r=>type_recv_r_p1
      type_send_p=>type_send_p_p1
      type_recv_p=>type_recv_p_p1
      call create_bc_mpi_datatype(mom(1),ndir)
      firstmf=.false.
    end if

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      block=>psa(igrid)
      call frictional_velocity(psa(igrid)%w,psa(igrid)%x,ixG^LL,ixM^LL)
    end do
    !$OMP END PARALLEL DO

    ! only update mf velocity in ghost cells
    stagger_flag=stagger_grid
    type_send_srl=>type_send_srl_p1
    type_recv_srl=>type_recv_srl_p1
    type_send_r=>type_send_r_p1
    type_recv_r=>type_recv_r_p1
    type_send_p=>type_send_p_p1
    type_recv_p=>type_recv_p_p1
    stagger_grid=.false.
    bcphys=.false.
    call getbc(qt,0.d0,psa,mom(1),ndir)
    bcphys=.true.
    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f
    stagger_grid=stagger_flag

  end subroutine mf_velocity_update

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine mf_add_source(qdt,dtfactor,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor
    double precision, intent(in)    :: wCT(ixI^S,1:nw),wCTprim(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active

    if (.not. qsourcesplit) then

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
      if (abs(mf_eta)>smalldouble)then
        active = .true.
        call add_source_res2(qdt,ixI^L,ixO^L,wCT,w,x)
      end if

      if (mf_eta_hyper>0.d0)then
        active = .true.
        call add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
      end if

    end if

    {^NOONED
    if(source_split_divb .eqv. qsourcesplit) then
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

  end subroutine mf_add_source

  subroutine frictional_velocity(w,x,ixI^L,ixO^L)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: xmin(ndim),xmax(ndim)
    double precision :: decay(ixO^S)
    double precision :: current(ixI^S,7-2*ndir:3),tmp(ixO^S)
    integer :: ix^D, idirmin,idir,jdir,kdir

    call get_current(w,ixI^L,ixO^L,idirmin,current)
    ! extrapolate current for the outmost layer, 
    ! because extrapolation of B at boundaries introduces artificial current
{
    if(block%is_physical_boundary(2*^D-1)) then
      if(slab) then
        ! exclude boundary 5 in 3D Cartesian
        if(2*^D-1==5.and.ndim==3) then
        else
          current(ixOmin^D^D%ixO^S,:)=current(ixOmin^D+1^D%ixO^S,:)
        end if
      else
        ! exclude boundary 1 spherical/cylindrical
        if(2*^D-1>1) then
          current(ixOmin^D^D%ixO^S,:)=current(ixOmin^D+1^D%ixO^S,:)
        end if
      end if
    end if
    if(block%is_physical_boundary(2*^D)) then
      current(ixOmax^D^D%ixO^S,:)=current(ixOmax^D-1^D%ixO^S,:)
    end if
\}
    w(ixO^S,mom(:))=0.d0
    ! calculate Lorentz force
    do idir=1,ndir; do jdir=idirmin,3; do kdir=1,ndir
      if(lvc(idir,jdir,kdir)/=0) then
        if(lvc(idir,jdir,kdir)==1) then
          w(ixO^S,mom(idir))=w(ixO^S,mom(idir))+current(ixO^S,jdir)*w(ixO^S,mag(kdir))
        else
          w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-current(ixO^S,jdir)*w(ixO^S,mag(kdir))
        end if
      end if
    end do; end do; end do

    tmp(ixO^S)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    ! frictional coefficient
    where(tmp(ixO^S)/=0.d0)
      tmp(ixO^S)=1.d0/(tmp(ixO^S)*mf_nu)
    endwhere

    do idir=1,ndir
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))*tmp(ixO^S)
    end do

    ! decay frictional velocity near selected boundaries
    ^D&xmin(^D)=xprobmin^D\
    ^D&xmax(^D)=xprobmax^D\
    decay(ixO^S)=1.d0
    do idir=1,ndim
      if(mf_decay_scale(2*idir-1)>0.d0) decay(ixO^S)=decay(ixO^S)*&
         (1.d0-exp(-(x(ixO^S,idir)-xmin(idir))/mf_decay_scale(2*idir-1)))
      if(mf_decay_scale(2*idir)>0.d0) decay(ixO^S)=decay(ixO^S)*&
         (1.d0-exp((x(ixO^S,idir)-xmax(idir))/mf_decay_scale(2*idir)))
    end do

    ! saturate mf velocity at mf_vmax
    tmp(ixO^S)=sqrt(sum(w(ixO^S,mom(:))**2,dim=ndim+1))/mf_vmax+1.d-12
    tmp(ixO^S)=dtanh(tmp(ixO^S))/tmp(ixO^S)
    do idir=1,ndir
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))*tmp(ixO^S)*decay(ixO^S)
    end do

  end subroutine frictional_velocity

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

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S)
    double precision :: gradeta(ixI^S,1:ndim), Bf(ixI^S,1:ndir)
    double precision :: tmp(ixI^S),tmp2(ixI^S)
    integer :: ixA^L,idir,jdir,kdir,idirmin,idim,jxO^L,hxO^L,ix
    integer :: lxO^L, kxO^L

    ! Calculating resistive sources involve one extra layer
    if (mf_4th_order) then
      ixA^L=ixO^L^LADD2;
    else
      ixA^L=ixO^L^LADD1;
    end if

    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixI^L,ixO^L,idirmin,current)

    if (mf_eta>zero)then
       eta(ixA^S)=mf_eta
       gradeta(ixO^S,1:ndim)=zero
    else
       call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,current,eta)
       ! assumes that eta is not function of current?
       do idim=1,ndim
          call gradient(eta,ixI^L,ixO^L,idim,tmp)
          gradeta(ixO^S,idim)=tmp(ixO^S)
       end do
    end if

    Bf(ixI^S,1:ndir)=wCT(ixI^S,mag(1:ndir))

    do idir=1,ndir
       ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
       if (mf_4th_order) then
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
       if (mf_eta<zero)then
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

    end do ! idir

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
    double precision :: tmpvec(ixI^S,1:3)
    integer :: ixA^L,idir,idirmin,idirmin1

    ! resistive term already added as an electric field component
    if(stagger_grid .and. ndim==ndir) return

    ixA^L=ixO^L^LADD2;

    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_res2: Non-conforming input limits")

    ixA^L=ixO^L^LADD1;
    ! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
    ! Determine exact value of idirmin while doing the loop.
    call get_current(wCT,ixI^L,ixA^L,idirmin,current)

    if (mf_eta>zero)then
       eta(ixA^S)=mf_eta
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
    ehyper(ixA^S,1:ndir) = - tmpvec(ixA^S,1:ndir)*mf_eta_hyper

    ixA^L=ixO^L;
    tmpvec2(ixA^S,1:ndir)=zero
    call curlvector(ehyper,ixI^L,ixA^L,tmpvec2,idirmin1,1,3)

    do idir=1,ndir
      w(ixO^S,mag(idir)) = w(ixO^S,mag(idir))-tmpvec2(ixO^S,idir)*qdt
    end do

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
    double precision :: gradPsi(ixI^S)
    integer          :: idim,idir

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb,mf_divb_nth)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (mf_glm_alpha < zero) then
      w(ixO^S,psi_) = abs(mf_glm_alpha)*wCT(ixO^S,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab_uniform) then
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*mf_glm_alpha/minval(dxlevel(:)))*w(ixO^S,psi_)
      else
        w(ixO^S,psi_) = dexp(-qdt*cmax_global*mf_glm_alpha/minval(block%ds(ixO^S,:),dim=ndim+1))*w(ixO^S,psi_)
      end if
    end if

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
          call gradient(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
       case("limited")
          call gradientL(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
       end select
    end do
  end subroutine add_source_glm

  !> Add divB related sources to w within ixO corresponding to Powel
  subroutine add_source_powel(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision                :: divb(ixI^S)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb,mf_divb_nth)

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*wCT(ixO^S,mom(idir))*divb(ixO^S)
    end do

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
    call get_divb(wCT,ixI^L,ixO^L,divb,mf_divb_nth)

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*wCT(ixO^S,mom(idir))*divb(ixO^S)
    end do

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
    call get_divb(wCT,ixI^L,ixp^L,divb,mf_divb_nth)

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
    end do

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
    invB(ixO^S)=sqrt(sum(w(ixO^S, mag(:))**2, dim=ndim+1))
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
  subroutine get_current(w,ixI^L,ixO^L,idirmin,current,fourthorder)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)  :: ixO^L, ixI^L
    double precision, intent(in) :: w(ixI^S,1:nw)
    integer, intent(out) :: idirmin
    logical, intent(in), optional :: fourthorder

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3)
    integer :: idir, idirmin0

    idirmin0 = 7-2*ndir

    if(present(fourthorder)) then
      call curlvector(w(ixI^S,mag(1:ndir)),ixI^L,ixO^L,current,idirmin,idirmin0,ndir,fourthorder)
    else
      call curlvector(w(ixI^S,mag(1:ndir)),ixI^L,ixO^L,current,idirmin,idirmin0,ndir)
    end if

  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine mf_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx^D
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision              :: dxarr(ndim)
    double precision              :: current(ixI^S,7-2*ndir:3),eta(ixI^S)
    integer                       :: idirmin,idim

    dtnew = bigdouble

    ^D&dxarr(^D)=dx^D;
    if (mf_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/mf_eta
    else if (mf_eta<zero)then
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

    if(mf_eta_hyper>zero) then
      if(slab_uniform) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/mf_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixO^S,1:ndim))**4/mf_eta_hyper,dtnew)
      end if
    end if
  end subroutine mf_get_dt

  ! Add geometrical source terms to w
  subroutine mf_add_source_geom(qdt,dtfactor, ixI^L,ixO^L,wCT,wprim,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), wprim(ixI^S,1:nw), w(ixI^S,1:nw)

    double precision :: tmp,invr,cs2
    integer          :: iw,idir,ix^D

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    select case (coordinate)
    case (cylindrical)
      if(phi_>0) then
        if(.not.stagger_grid) then
         {do ix^DB=ixOmin^DB,ixOmax^DB\}
            w(ix^D,bphi_)=w(ix^D,bphi_)+qdt/x(ix^D,1)*&
                     (wCT(ix^D,bphi_)*wCT(ix^D,mr_) &
                     -wCT(ix^D,br_)*wCT(ix^D,mphi_))
         {end do\}
        end if
      end if
      if(mf_glm) w(ixO^S,br_)=w(ixO^S,br_)+qdt*wCT(ixO^S,psi_)/x(ixO^S,1)
    case (spherical)
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          invr=qdt/x(ix^D,1)
         ! b1
         if(mf_glm) then
           w(ix^D,mag(1))=w(ix^D,mag(1))+invr*2.0d0*wCT(ix^D,psi_)
         end if

         {^NOONED
         ! b2
         if(.not.stagger_grid) then
           tmp=wCT(ix^D,mom(1))*wCT(ix^D,mag(2))-wCT(ix^D,mom(2))*wCT(ix^D,mag(1))
           cs2=1.d0/tan(x(ix^D,2))
           if(mf_glm) then
             tmp=tmp+cs2*wCT(ix^D,psi_)
           end if
           w(ix^D,mag(2))=w(ix^D,mag(2))+tmp*invr
         end if
         }

         if(ndir==3) then
          {^IFONED
          ! b3
          if(.not.stagger_grid) then
            w(ix^D,mag(3))=w(ix^D,mag(3))+invr*&
               (wCT(ix^D,mom(1))*wCT(ix^D,mag(3)) &
               -wCT(ix^D,mom(3))*wCT(ix^D,mag(1)))
          end if
          }
          {^NOONED
          ! b3
          if(.not.stagger_grid) then
            w(ix^D,mag(3))=w(ix^D,mag(3))+invr*&
               (wCT(ix^D,mom(1))*wCT(ix^D,mag(3)) &
               -wCT(ix^D,mom(3))*wCT(ix^D,mag(1)) &
              -(wCT(ix^D,mom(3))*wCT(ix^D,mag(2)) &
               -wCT(ix^D,mom(2))*wCT(ix^D,mag(3)))*cs2)
          end if
          }
         end if
       {end do\}
    end select

  end subroutine mf_add_source_geom

  subroutine mf_modify_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    double precision, intent(inout) :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw)
    type(state)                     :: s
    double precision                :: dB(ixI^S), dPsi(ixI^S)

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

  end subroutine mf_modify_wLR

  subroutine mf_boundary_adjust(igrid,psb)
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

  end subroutine mf_boundary_adjust

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
  subroutine mf_clean_divb_multigrid(qdt, qt, active)
    use mod_forest
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_geometry

    double precision, intent(in) :: qdt    !< Current time step
    double precision, intent(in) :: qt     !< Current time
    logical, intent(inout)       :: active !< Output if the source is active

    double precision             :: tmp(ixG^T), grad(ixG^T, ndim)
    double precision             :: res
    double precision, parameter  :: max_residual = 1d-3
    double precision, parameter  :: residual_reduction = 1d-10
    integer                      :: id
    integer, parameter           :: max_its      = 50
    double precision             :: residual_it(max_its), max_divb
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
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_periodic)
          ! Nothing to do here
       case default
          write(*,*) "mf_clean_divb_multigrid warning: unknown boundary type"
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
            mf_divb_nth)
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
           ps(igrid)%ws(ixC^S,idim)=ps(igrid)%ws(ixC^S,idim)-grad(ixC^S,idim)
         end do
         call mf_face_to_center(ixM^LL,ps(igrid))
       else
         do idim = 1, ndim
            call gradient(tmp,ixG^LL,ixM^LL,idim,grad(ixG^T, idim))
         end do
         ! Apply the correction B* = B - gradient(phi)
         ps(igrid)%w(ixM^T, mag(1:ndim)) = &
              ps(igrid)%w(ixM^T, mag(1:ndim)) - grad(ixM^T, :)
       end if

    end do

    active = .true.

  end subroutine mf_clean_divb_multigrid
  }

  !> get electric field though averaging neighors to update faces in CT
  subroutine mf_update_faces_average(ixI^L,ixO^L,qt,qdt,wp,fC,fE,sCT,s,vcts)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wp(ixI^S,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,sdim:3)

    double precision                   :: circ(ixI^S,1:ndim)
    double precision                   :: E_resi(ixI^S,sdim:3)
    integer                            :: ix^D,ixC^L,ixA^L,i1kr^D,i2kr^D
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! if there is resistivity, get eta J
    if(mf_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,E_resi)

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
              if(mf_eta/=zero) fE(ix^D,idir)=fE(ix^D,idir)+E_resi(ix^D,idir)
              if(record_electric_field) s%we(ix^D,idir)=fE(ix^D,idir)
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
      circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      ! Time update
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

    end associate

  end subroutine mf_update_faces_average

  !> update faces using UCT contact mode by Gardiner and Stone 2005 JCP 205, 509
  subroutine mf_update_faces_contact(ixI^L,ixO^L,qt,qdt,wp,fC,fE,sCT,s,vcts)
    use mod_global_parameters
    use mod_usr_methods

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
    ! gradient of E at left and right side of a cell face
    double precision                   :: EL(ixI^S),ER(ixI^S)
    ! gradient of E at left and right side of a cell corner
    double precision                   :: ELC(ixI^S),ERC(ixI^S)
    ! current on cell edges
    double precision                   :: jce(ixI^S,sdim:3)
    integer                            :: hxC^L,ixC^L,jxC^L,ixA^L,ixB^L
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(bfaces=>s%ws,x=>s%x,w=>s%w,vnorm=>vcts%vnorm)

    ECC=0.d0
    ! Calculate electric field at cell centers
    do idim1=1,ndim; do idim2=1,ndim; do idir=sdim,3
      if(lvc(idim1,idim2,idir)==1)then
         ECC(ixI^S,idir)=ECC(ixI^S,idir)+wp(ixI^S,mag(idim1))*wp(ixI^S,mom(idim2))
      else if(lvc(idim1,idim2,idir)==-1) then
         ECC(ixI^S,idir)=ECC(ixI^S,idir)-wp(ixI^S,mag(idim1))*wp(ixI^S,mom(idim2))
      endif
    enddo; enddo; enddo

    ! if there is resistivity, get eta J
    if(mf_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,jce)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ! evaluate electric field along cell edges according to equation (41)
    do idim1=1,ndim 
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=sdim,3 ! Direction of line integral
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
            if(mf_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+jce(ixC^S,idir)

            if(record_electric_field) s%we(ixC^S,idir)=fE(ixC^S,idir)
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

  end subroutine mf_update_faces_contact

  !> update faces
  subroutine mf_update_faces_hll(ixI^L,ixO^L,qt,qdt,wp,fC,fE,sCT,s,vcts)
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
    double precision                   :: btilL(ixI^S,ndim)
    double precision                   :: btilR(ixI^S,ndim)
    double precision                   :: cp(ixI^S,2)
    double precision                   :: cm(ixI^S,2)
    double precision                   :: circ(ixI^S,1:ndim)
    ! current on cell edges
    double precision                   :: jce(ixI^S,sdim:3)
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
    if(mf_eta/=zero) call get_resistive_electric_field(ixI^L,ixO^L,sCT,s,jce)

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
      call reconstruct(ixI^L,ixC^L,idim2,bfacesCT(ixI^S,idim1),&
               btilL(ixI^S,idim1),btilR(ixI^S,idim1))

      call reconstruct(ixI^L,ixC^L,idim1,bfacesCT(ixI^S,idim2),&
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
      if(mf_eta/=zero) fE(ixC^S,idir)=fE(ixC^S,idir)+jce(ixC^S,idir)

      if(record_electric_field) s%we(ixC^S,idir)=fE(ixC^S,idir)
      ! times time step and edge length 
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
      ! Time update
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

    end associate
  end subroutine mf_update_faces_hll

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
          ixCmax^D=ixOmax^D+1;
          ixCmin^D=ixOmin^D+kr(idir,^D)-2;
          ixBmax^D=ixCmax^D-kr(idir,^D)+1;
          ixBmin^D=ixCmin^D;
          ! current at transverse faces
          xs(ixB^S,:)=x(ixB^S,:)
          xs(ixB^S,idim2)=x(ixB^S,idim2)+half*dx(ixB^S,idim2)
          call gradientF(wCTs(ixGs^T,idim2),xs,ixGs^LL,ixC^L,idim1,gradi)
          if (lvc(idim1,idim2,idir)==1) then
            jce(ixC^S,idir)=jce(ixC^S,idir)+gradi(ixC^S)
          else
            jce(ixC^S,idir)=jce(ixC^S,idir)-gradi(ixC^S)
          end if
        end do
      end do
    end do
    ! get resistivity
    if(mf_eta>zero)then
      jce(ixI^S,:)=jce(ixI^S,:)*mf_eta
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
      enddo
    end if

    end associate
  end subroutine get_resistive_electric_field

  !> calculate cell-center values from face-center values
  subroutine mf_face_to_center(ixO^L,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixO^L
    type(state)                        :: s

    integer :: ix^D

    ! calculate cell-center values from face-center values in 2nd order
   {!dir$ ivdep
    do ix^DB=ixOmin^DB,ixOmax^DB\}
      {^IFTHREED
      s%w(ix^D,mag(1))=half/s%surface(ix^D,1)*(s%ws(ix^D,1)*s%surfaceC(ix^D,1)&
        +s%ws(ix1-1,ix2,ix3,1)*s%surfaceC(ix1-1,ix2,ix3,1))
      s%w(ix^D,mag(2))=half/s%surface(ix^D,2)*(s%ws(ix^D,2)*s%surfaceC(ix^D,2)&
        +s%ws(ix1,ix2-1,ix3,2)*s%surfaceC(ix1,ix2-1,ix3,2))
      s%w(ix^D,mag(3))=half/s%surface(ix^D,3)*(s%ws(ix^D,3)*s%surfaceC(ix^D,3)&
        +s%ws(ix1,ix2,ix3-1,3)*s%surfaceC(ix1,ix2,ix3-1,3))
      }
      {^IFTWOD
      s%w(ix^D,mag(1))=half/s%surface(ix^D,1)*(s%ws(ix^D,1)*s%surfaceC(ix^D,1)&
        +s%ws(ix1-1,ix2,1)*s%surfaceC(ix1-1,ix2,1))
      s%w(ix^D,mag(2))=half/s%surface(ix^D,2)*(s%ws(ix^D,2)*s%surfaceC(ix^D,2)&
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

  end subroutine mf_face_to_center

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

  subroutine record_force_free_metrics()
    use mod_global_parameters

    double precision :: sum_jbb_ipe, sum_j_ipe, sum_l_ipe, f_i_ipe, volume_pe
    double precision :: sum_jbb, sum_j, sum_l, f_i, volume, cw_sin_theta
    integer :: iigrid, igrid, ix^D
    integer :: amode, istatus(MPI_STATUS_SIZE)
    integer, save :: fhmf
    logical :: patchwi(ixG^T)
    logical, save :: logmfopened=.false.
    character(len=800) :: filename,filehead
    character(len=800) :: line,datastr

    sum_jbb_ipe = 0.d0
    sum_j_ipe = 0.d0
    sum_l_ipe = 0.d0
    f_i_ipe = 0.d0
    volume_pe=0.d0
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      call mask_inner(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,patchwi,volume_pe)
      sum_jbb_ipe = sum_jbb_ipe+integral_grid_mf(ixG^LL,ixM^LL,ps(igrid)%w,&
        ps(igrid)%x,1,patchwi)
      sum_j_ipe   = sum_j_ipe+integral_grid_mf(ixG^LL,ixM^LL,ps(igrid)%w,&
        ps(igrid)%x,2,patchwi)
      f_i_ipe=f_i_ipe+integral_grid_mf(ixG^LL,ixM^LL,ps(igrid)%w,&
        ps(igrid)%x,3,patchwi)
      sum_l_ipe   = sum_l_ipe+integral_grid_mf(ixG^LL,ixM^LL,ps(igrid)%w,&
        ps(igrid)%x,4,patchwi)
    end do
    call MPI_REDUCE(sum_jbb_ipe,sum_jbb,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,0,icomm,ierrmpi)
    call MPI_REDUCE(sum_j_ipe,sum_j,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       icomm,ierrmpi)
    call MPI_REDUCE(f_i_ipe,f_i,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       icomm,ierrmpi)
    call MPI_REDUCE(sum_l_ipe,sum_l,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       icomm,ierrmpi)
    call MPI_REDUCE(volume_pe,volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       icomm,ierrmpi)

    if(mype==0) then
      ! current- and volume-weighted average of the sine of the angle
      ! between the magnetic field B and the current density J
      cw_sin_theta = sum_jbb/sum_j
      ! volume-weighted average of the absolute value of the fractional
      ! magnetic flux change
      f_i = f_i/volume
      sum_j=sum_j/volume
      sum_l=sum_l/volume
      if(.not.logmfopened) then
        ! generate filename
        write(filename,"(a,a)") TRIM(base_filename), "_mf_metrics.csv"

        ! Delete the log when not doing a restart run
        if (restart_from_file == undefined) then
           open(unit=20,file=trim(filename),status='replace')
           close(20, status='delete')
        end if

        amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        amode=ior(amode,MPI_MODE_APPEND)
        call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,fhmf,ierrmpi)
        logmfopened=.true.
        filehead="    it, time, CW_sin_theta, current, Lorenz_force, f_i"
        if (it == 0) then
          call MPI_FILE_WRITE(fhmf,filehead,len_trim(filehead), &
                              MPI_CHARACTER,istatus,ierrmpi)
          call MPI_FILE_WRITE(fhmf,achar(10),1,MPI_CHARACTER,istatus,ierrmpi)
        end if
      end if
      line=''
      write(datastr,'(i8,a)') it,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6,a)') global_time,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6,a)') cw_sin_theta,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6,a)') sum_j,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6,a)') sum_l,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6)') f_i
      line=trim(line)//trim(datastr)//new_line('A')
      call MPI_FILE_WRITE(fhmf,line,len_trim(line),MPI_CHARACTER,istatus,ierrmpi)
      if(.not.time_advance) then
        call MPI_FILE_CLOSE(fhmf,ierrmpi)
        logmfopened=.false.
      end if
    end if

  end subroutine record_force_free_metrics

  subroutine mask_inner(ixI^L,ixO^L,w,x,patchwi,volume_pe)
    use mod_global_parameters

    integer, intent(in)         :: ixI^L,ixO^L
    double precision, intent(in):: w(ixI^S,nw),x(ixI^S,1:ndim)
    double precision, intent(inout) :: volume_pe
    logical, intent(out)        :: patchwi(ixI^S)

    double precision            :: xO^L
    integer                     :: ix^D

    {xOmin^D = xprobmin^D + 0.05d0*(xprobmax^D-xprobmin^D)\}
    {xOmax^D = xprobmax^D - 0.05d0*(xprobmax^D-xprobmin^D)\}
    if(slab) then
      xOmin^ND = xprobmin^ND
    else
      xOmin1 = xprobmin1
    end if

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if({ x(ix^DD,^D) > xOmin^D .and. x(ix^DD,^D) < xOmax^D | .and. }) then
          patchwi(ix^D)=.true.
          volume_pe=volume_pe+block%dvolume(ix^D)
        else
          patchwi(ix^D)=.false.
        endif
    {end do\}

  end subroutine mask_inner

  function integral_grid_mf(ixI^L,ixO^L,w,x,iw,patchwi)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L,iw
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    logical, intent(in) :: patchwi(ixI^S)

    double precision, dimension(ixI^S,7-2*ndir:3) :: current
    double precision, dimension(ixI^S,1:ndir) :: bvec,qvec
    double precision :: integral_grid_mf,tmp(ixI^S),bm2
    integer :: ix^D,idirmin0,idirmin,idir,jdir,kdir

    integral_grid_mf=0.d0
    select case(iw)
     case(1)
      ! Sum(|JxB|/|B|*dvolume)
      bvec(ixO^S,:)=w(ixO^S,mag(:))
      call get_current(w,ixI^L,ixO^L,idirmin,current)
      ! calculate Lorentz force
      qvec(ixO^S,1:ndir)=zero
      do idir=1,ndir; do jdir=idirmin,3; do kdir=1,ndir
         if(lvc(idir,jdir,kdir)/=0)then
           if(lvc(idir,jdir,kdir)==1)then
             qvec(ixO^S,idir)=qvec(ixO^S,idir)+current(ixO^S,jdir)*bvec(ixO^S,kdir)
           else
             qvec(ixO^S,idir)=qvec(ixO^S,idir)-current(ixO^S,jdir)*bvec(ixO^S,kdir)
           end if
         end if
      end do; end do; end do

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) then
           bm2=sum(bvec(ix^D,:)**2)
           if(bm2/=0.d0) bm2=1.d0/bm2
           integral_grid_mf=integral_grid_mf+sqrt(sum(qvec(ix^D,:)**2)*&
                           bm2)*block%dvolume(ix^D)
         end if
      {end do\}
     case(2)
      ! Sum(|J|*dvolume)
      call get_current(w,ixI^L,ixO^L,idirmin,current)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) integral_grid_mf=integral_grid_mf+sqrt(sum(current(ix^D,:)**2))*&
                           block%dvolume(ix^D)
      {end do\}
     case(3)
      ! f_i solenoidal property of B: (dvolume |div B|)/(dsurface |B|)
      ! Sum(f_i*dvolume)
      call get_normalized_divb(w,ixI^L,ixO^L,tmp)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) integral_grid_mf=integral_grid_mf+tmp(ix^D)*block%dvolume(ix^D)
      {end do\}
     case(4)
      ! Sum(|JxB|*dvolume)
      bvec(ixO^S,:)=w(ixO^S,mag(:))
      call get_current(w,ixI^L,ixO^L,idirmin,current)
      ! calculate Lorentz force
      qvec(ixO^S,1:ndir)=zero
      do idir=1,ndir; do jdir=idirmin,3; do kdir=1,ndir
         if(lvc(idir,jdir,kdir)/=0)then
           if(lvc(idir,jdir,kdir)==1)then
             qvec(ixO^S,idir)=qvec(ixO^S,idir)+current(ixO^S,jdir)*bvec(ixO^S,kdir)
           else
             qvec(ixO^S,idir)=qvec(ixO^S,idir)-current(ixO^S,jdir)*bvec(ixO^S,kdir)
           end if
         end if
      end do; end do; end do

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) integral_grid_mf=integral_grid_mf+&
                            sqrt(sum(qvec(ix^D,:)**2))*block%dvolume(ix^D)
      {end do\}
    end select
    return
  end function integral_grid_mf

end module mod_mf_phys
