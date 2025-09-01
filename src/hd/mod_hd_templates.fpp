#:if PHYS == 'hd'
  
#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter, public              :: nw_phys=2+ndim
  
  !> Whether an energy equation is used
  logical, public                         :: hd_energy = .true.
  !$acc declare copyin(hd_energy)

  !> Index of the density (in the w array)
  integer, public                         :: rho_
  !$acc declare create(rho_)

  !> Indices of the momentum density
  integer, allocatable, public            :: mom(:)
  !$acc declare create(mom)

  !> Index of the energy density (-1 if not present)
  integer, public                         :: e_
  !$acc declare create(e_)

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public                         :: p_
  !$acc declare create(p_)

  !> The adiabatic index
  double precision, public                :: hd_gamma = 5.d0/3.0d0
  !$acc declare copyin(hd_gamma)

  !> The adiabatic constant
  double precision, public                :: hd_adiab = 1.0d0
  !$acc declare copyin(hd_adiab)

  !> The helium abundance
  double precision, public                :: He_abundance=0.1d0
  !$acc declare copyin(He_abundance)

  !> Whether plasma is partially ionized
  logical, public                         :: hd_partial_ionization = .false.
  !$acc declare copyin(hd_partial_ionization)
  
  !> Whether to use gravity
  logical, public                         :: hd_gravity = .false.
  !$acc declare copyin(hd_gravity)

  !> Allows overruling default corner filling (for debug mode, since otherwise corner primitives fail)
  logical, public                         :: hd_force_diagonal = .false.
  !$acc declare copyin(hd_force_diagonal)

  !> Whether particles module is added
  logical, public                         :: hd_particles = .false.
  !$acc declare copyin(hd_particles)

#:enddef

#:def read_params()
    !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /hd_list/ hd_energy, hd_gamma, hd_adiab, hd_partial_ionization,&
        hd_force_diagonal, hd_particles, hd_gravity

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hd_list, end=111)
111    close(unitpar)
    end do

#ifdef _OPENACC
 !$acc update device(hd_energy, hd_gamma, hd_adiab, hd_partial_ionization, hd_force_diagonal, hd_particles, hd_gravity)
#endif

  end subroutine read_params
#:enddef

#:def phys_activate() 
  subroutine phys_activate()
    call phys_init()
  end subroutine phys_activate
#:enddef
  #:def phys_units()
  subroutine phys_units()
    use mod_global_parameters
    double precision :: mp, kB
    double precision :: a,b

    !> here no SI_UNIT used by default, to be implemented
    mp = mp_cgs
    kB = kB_cgs
    !> eq_state_units by default, to be implemented
    a = 1.d0+4.d0*He_abundance
    b = 2.d0+3.d0*He_abundance

    if(unit_density/=1.d0 .or. unit_numberdensity/=1.d0) then
      if(unit_density/=1.d0) then
        unit_numberdensity=unit_density/(a*mp)
      else if(unit_numberdensity/=1.d0) then
        unit_density=a*mp*unit_numberdensity
      end if
      if(unit_temperature/=1.d0) then
        unit_pressure=b*unit_numberdensity*kB*unit_temperature
        unit_velocity=dsqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_pressure/=1.d0) then
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        unit_velocity=dsqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_velocity/=1.d0) then
        unit_pressure=unit_density*unit_velocity**2
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=1.d0) then
        unit_velocity=unit_length/unit_time
        unit_pressure=unit_density*unit_velocity**2
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    else if(unit_temperature/=1.d0) then
      ! units of temperature and velocity are dependent
      if(unit_pressure/=1.d0) then
        unit_numberdensity=unit_pressure/(b*unit_temperature*kB)
        unit_density=a*mp*unit_numberdensity
        unit_velocity=dsqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      end if
    else if(unit_pressure/=1.d0) then
      if(unit_velocity/=1.d0) then
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=0.d0) then
        unit_velocity=unit_length/unit_time
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    end if
    unit_mass=unit_density*unit_length**3

    !$acc update device(unit_density, unit_numberdensity, unit_temperature, unit_pressure, unit_velocity, unit_length, unit_time, unit_mass)
  end subroutine phys_units
#:enddef

#:def phys_init()
    !> Initialize the module
  subroutine phys_init()
    use mod_global_parameters
!    use mod_particles, only: particles_init

    call phys_units()
    call read_params(par_files)

    phys_energy  = hd_energy
    phys_total_energy  = hd_energy
    phys_internal_e = .false.
    phys_gamma = hd_gamma
    phys_partial_ionization=hd_partial_ionization
 !$acc update device(physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization)

    use_particles = hd_particles

    ! Determine flux variables
    rho_ = var_set_rho()
    !$acc update device(rho_)

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)
    !$acc update device(mom)

    ! Set index of energy variable
    if (hd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if
    !$acc update device(e_,p_)

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    if (hd_force_diagonal) then
       ! ensure corners are filled, otherwise divide by zero when getting primitives
       !  --> only for debug purposes
       phys_req_diagonal = .true.
    endif

    ! set number of variables which need update ghostcells
    nwgc=nwflux
    !$acc update device(nwgc)

! use cycle, needs to be dealt with:    
!    ! Initialize particles module
!    if (hd_particles) then
!       call particles_init()
!       phys_req_diagonal = .true.
!    end if

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1
    !$acc update device(nvector, iw_vector)
    !$acc update device(phys_req_diagonal)

  end subroutine phys_init
#:enddef

#:def phys_get_dt()
  subroutine phys_get_dt(w, x, dx, dtnew)
  !$acc routine seq
#:if defined('GRAVITY')
  use mod_usr, only: gravity_field
#:endif    
    real(dp), intent(in)   :: w(nw_phys), x(1:ndim), dx(1:ndim)
    real(dp), intent(out)  :: dtnew
    ! .. local ..
    integer                :: idim
    real(dp)               :: field

    dtnew = huge(1.0d0)
    
#:if defined('GRAVITY')
    do idim = 1, ndim
       field = gravity_field(w, x, idim)
       field = max( abs(field), epsilon(1.0d0) )
       dtnew = min( dtnew, 1_dp / sqrt( field/dx(idim) ) )
    end do
#:endif    
    
  end subroutine phys_get_dt
#:enddef  

#:def addsource_local()
subroutine addsource_local(qdt, dtfactor, qtC, wCT, wCTprim, qt, wnew, x,&
    qsourcesplit)
  !$acc routine seq
#:if defined('GRAVITY')
  use mod_usr, only: gravity_field
#:endif    
  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCT(nw_phys), wCTprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit
  ! .. local ..
  integer                  :: idim
  real(dp)                 :: field

#:if defined('GRAVITY')
  do idim = 1, ndim
     field = gravity_field(wCT, x, idim)
     wnew(iw_mom(idim)) = wnew(iw_mom(idim)) + qdt * field * wCT(iw_rho)
     wnew(iw_e)         = wnew(iw_e) + qdt * field * wCT(iw_mom(idim))
  end do
#:endif  

end subroutine addsource_local
#:enddef

#:def to_primitive()
pure subroutine to_primitive(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(nw_phys)

  
       u(iw_mom(1)) = u(iw_mom(1))/u(iw_rho)
  
  
       u(iw_mom(2)) = u(iw_mom(2))/u(iw_rho)
  
  
       u(iw_mom(3)) = u(iw_mom(3))/u(iw_rho)
  

  u(iw_e) = (hd_gamma-1.0_dp) * (u(iw_e) - 0.5_dp * u(iw_rho) * &
     sum(u(iw_mom(1:ndim))**2) )

end subroutine to_primitive
#:enddef

#:def to_conservative()  
pure subroutine to_conservative(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(nw_phys)
  real(dp)                :: inv_gamma_m1

  inv_gamma_m1 = 1.0d0/(hd_gamma - 1.0_dp)

  ! Compute energy from pressure and kinetic energy
  u(iw_e) = u(iw_e) * inv_gamma_m1 + 0.5_dp * u(iw_rho) * &
     sum(u(iw_mom(1:ndim))**2)

  ! Compute momentum from density and velocity components
  
       u(iw_mom(1)) = u(iw_rho) * u(iw_mom(1))
  
  
       u(iw_mom(2)) = u(iw_rho) * u(iw_mom(2))
  
  
       u(iw_mom(3)) = u(iw_rho) * u(iw_mom(3))
  

end subroutine to_conservative
#:enddef

#:def get_flux()
subroutine get_flux(u, xC, flux_dim, flux)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: flux(nw_phys)
  real(dp)              :: inv_gamma_m1

  inv_gamma_m1 = 1.0d0/(hd_gamma - 1.0_dp)

  ! Density flux
  flux(iw_rho) = u(iw_rho) * u(iw_mom(flux_dim))

  ! Momentum flux with pressure term
  
       flux(iw_mom(1)) = u(iw_rho) * u(iw_mom(1)) * u(iw_mom(flux_dim))
  
  
       flux(iw_mom(2)) = u(iw_rho) * u(iw_mom(2)) * u(iw_mom(flux_dim))
  
  
       flux(iw_mom(3)) = u(iw_rho) * u(iw_mom(3)) * u(iw_mom(flux_dim))
  
  flux(iw_mom(flux_dim)) = flux(iw_mom(flux_dim)) + u(iw_e)

  ! Energy flux
  flux(iw_e) = u(iw_mom(flux_dim)) * (u(iw_e) * inv_gamma_m1 + 0.5_dp * &
     u(iw_rho) * sum(u(iw_mom(1:ndim))**2) + u(iw_e))

end subroutine get_flux
#:enddef

#:def get_cmax()  
pure real(dp) function get_cmax(u, x, flux_dim) result(wC)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim

  wC = sqrt(hd_gamma * u(iw_e) / u(iw_rho)) + abs(u(iw_mom(flux_dim)))

end function get_cmax
#:enddef  

#:def get_gradientT()
pure real(dp) function get_gradientT(u, x, grad_dim) result(gradT)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys, 5)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: grad_dim
  real(dp) :: Te(5),Tface(2)

  Te(1:5)=u(iw_e,1:5)/u(iw_rho,1:5)
  Tface(1)=(7.0d0*(Te(2)+Te(3))-(Te(1)+Te(4)))/12.0d0
  Tface(2)=(7.0d0*(Te(3)+Te(4))-(Te(2)+Te(5)))/12.0d0
  gradT=Tface(2)-Tface(1)

end function get_gradientT

#:enddef
#:endif
