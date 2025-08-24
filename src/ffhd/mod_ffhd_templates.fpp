#:if PHYS == 'ffhd'
  
#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  !> Only consider the hyperbolic thermal conduction situation
  integer, parameter, public              :: nw_phys=4
  
  !> Whether an energy equation is used
  logical, public                         :: ffhd_energy = .true.

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

    !> Index of the hyperbolic flux variable
  integer, public                         :: q_
  !$acc declare create(q_)

  !> The adiabatic index
  double precision, public                :: ffhd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: ffhd_adiab = 1.0d0

  !> The helium abundance
  double precision, public                :: He_abundance=0.1d0

  !> The thermal conductivity kappa in hyperbolic thermal conduction
  double precision, public                :: hypertc_kappa
  !$acc declare copyin(hypertc_kappa)

  !> Whether plasma is partially ionized
  logical, public                         :: ffhd_partial_ionization = .false.

#:enddef

#:def read_params()
    !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /ffhd_list/ ffhd_energy, ffhd_gamma, ffhd_partial_ionization

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, ffhd_list, end=111)
111    close(unitpar)
    end do

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
    #:if defined('COOLING')
    use mod_radiative_cooling, only: radiative_cooling_init_params, radiative_cooling_init
    #:endif

    call phys_units()
    call read_params(par_files)

    phys_energy  = ffhd_energy
    phys_total_energy  = ffhd_energy
    phys_internal_e = .false.
    phys_gamma = ffhd_gamma
    phys_partial_ionization=ffhd_partial_ionization
 !$acc update device(physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization)

    ! Determine flux variables
    rho_ = var_set_rho()
    !$acc update device(rho_)

    allocate(mom(1))
    mom(:) = var_set_momentum(1)
    !$acc update device(mom)

    ! Set index of energy variable
    if (ffhd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if
    !$acc update device(e_,p_)

    ! Set index for heat flux
    q_ = var_set_q()
    !$acc update device(q_)

    ! set number of variables which need update ghostcells
    nwgc=nwflux
    !$acc update device(nwgc)

    #:if defined('COOLING')
    allocate(rc_fl)
    call radiative_cooling_init_params(phys_gamma,He_abundance)
    call radiative_cooling_init(rc_fl)
    !$acc update device(rc_fl)
    #:endif

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
#:if defined('BFIELD')
  use mod_usr, only: bfield
#:endif    
#:if defined('COOLING')
  use mod_radiative_cooling, only: radiative_cooling_add_source
#:endif

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCT(nw_phys), wCTprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit
  ! .. local ..
  integer                  :: idim
  real(dp)                 :: field, mag

#:if defined('GRAVITY')
  do idim = 1, ndim
     field = gravity_field(wCT, x, idim)
     mag = bfield(x, idim)
     wnew(iw_mom(1)) = wnew(iw_mom(1)) + qdt * field * wCT(iw_rho) * mag
     wnew(iw_e)         = wnew(iw_e) + qdt * field * wCT(iw_mom(1)) * mag
  end do
#:endif  
#:if defined('BFIELD')
!> p*divb to be added here
#:endif  

#:if defined('COOLING')
  call radiative_cooling_add_source(qdt,wCT,wCTprim,wnew,x,rc_fl)
#:endif

end subroutine addsource_local
#:enddef

#:def to_primitive()
pure subroutine to_primitive(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(nw_phys)

  u(iw_mom(1)) = u(iw_mom(1))/u(iw_rho)

  u(iw_e) = (phys_gamma-1.0_dp) * (u(iw_e) - 0.5_dp * u(iw_rho) * &
     u(iw_mom(1))**2 )

end subroutine to_primitive
#:enddef

#:def to_conservative()  
pure subroutine to_conservative(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(nw_phys)
  real(dp)                :: inv_gamma_m1

  inv_gamma_m1 = 1.0d0/(phys_gamma - 1.0_dp)

  ! Compute energy from pressure and kinetic energy
  u(iw_e) = u(iw_e) * inv_gamma_m1 + 0.5_dp * u(iw_rho) * &
     u(iw_mom(1))**2

  ! Compute momentum from density and velocity components
  u(iw_mom(1)) = u(iw_rho) * u(iw_mom(1))
  
end subroutine to_conservative
#:enddef

#:def get_flux()
subroutine get_flux(u, xC, flux_dim, flux)
  !$acc routine seq
#:if defined('BFIELD')
  use mod_usr, only: bfield
#:endif    
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: xC(1:ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: flux(nw_phys)
  real(dp)              :: inv_gamma_m1
  real(dp)              :: mag

  inv_gamma_m1 = 1.0d0/(phys_gamma - 1.0_dp)

  mag = 1.0d0
#:if defined('BFIELD')
  mag = bfield(xC, flux_dim)
#:endif    

  ! Density flux
  flux(iw_rho) = u(iw_rho) * u(iw_mom(1)) * mag

  ! Momentum flux with pressure term
  flux(iw_mom(1)) = (u(iw_rho)*u(iw_mom(1))**2 + u(iw_e)) * mag
  
  ! Energy flux with hyperbolic conduction included
  flux(iw_e) = ( u(iw_mom(1))*(u(iw_e)*inv_gamma_m1 + &
               0.5_dp*u(iw_rho)*u(iw_mom(1))**2 + u(iw_e)) + &
               u(iw_q)) * mag

  ! heat flux, to be added
  flux(iw_q) = 0.0d0

end subroutine get_flux
#:enddef

#:def get_cmax()  
pure real(dp) function get_cmax(u, x, flux_dim) result(wC)
  !$acc routine seq
#:if defined('BFIELD')
  use mod_usr, only: bfield
#:endif    
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim
  real(dp)              :: mag

  mag = 1.0d0
#:if defined('BFIELD')
  mag = bfield(x, flux_dim)
#:endif    
  
  wC = dsqrt(phys_gamma*u(iw_e)/u(iw_rho)) + abs(u(iw_mom(1))*mag)

end function get_cmax
#:enddef  

#:def get_rho()
pure double precision function get_rho(w, x) result(rho)
  !$acc routine seq
  double precision, intent(in)  :: w(nwflux)
  double precision, intent(in)  :: x(1:ndim)

  rho = w(iw_rho)
end function get_rho
#:enddef

#:def get_pthermal()
pure double precision function get_pthermal(w, x) result(pth)
  !$acc routine seq
  double precision, intent(in)  :: w(nwflux)
  double precision, intent(in)  :: x(1:ndim)

  pth = (phys_gamma-1.0_dp)*(w(iw_e)-0.5_dp*w(iw_mom(1))**2/w(iw_rho))
end function get_pthermal
#:enddef

#:def get_Rfactor()
pure double precision function get_Rfactor() result(Rfactor)
  !$acc routine seq
  Rfactor = 1.0d0
end function get_Rfactor
#:enddef

#:endif
