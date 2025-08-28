#:if PHYS == 'ffhd'
  
#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  !> Only consider ffhd with hyperbolic thermal conduction hence 4 variables
  integer, parameter, public              :: nw_phys=4
  
  !> Whether an energy equation is used
  logical, public                         :: ffhd_energy = .true.
  !$acc declare copyin(ffhd_energy)

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
  !$acc declare copyin(ffhd_gamma)

  !> The adiabatic constant
  double precision, public                :: ffhd_adiab = 1.0d0
  !$acc declare copyin(ffhd_adiab)
  
  !> The helium abundance
  double precision, public                :: He_abundance=0.1d0
  !$acc declare copyin(He_abundance)

  !> The thermal conductivity kappa in hyperbolic thermal conduction
  double precision, public                :: hypertc_kappa = -1.0d0
  !$acc declare copyin(hypertc_kappa)

  !> Whether plasma is partially ionized
  logical, public                         :: ffhd_partial_ionization = .false.
  !$acc declare copyin(ffhd_partial_ionization)

  !> Allows overruling default corner filling (for debug mode, since otherwise corner primitives fail)
  logical, public                         :: ffhd_force_diagonal = .false.
  !$acc declare copyin(ffhd_force_diagonal)

  !> Whether particles module is added
  logical, public                         :: ffhd_particles = .false.
  !$acc declare copyin(ffhd_particles)

#:enddef

#:def read_params()
    !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /ffhd_list/ ffhd_energy, ffhd_gamma, ffhd_partial_ionization,&
        ffhd_force_diagonal, ffhd_particles

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, ffhd_list, end=111)
111    close(unitpar)
    end do

#ifdef _OPENACC
 !$acc update device(ffhd_energy, ffhd_gamma, ffhd_adiab, ffhd_partial_ionization, ffhd_force_diagonal, ffhd_particles)
#endif

  end subroutine read_params
#:enddef

#:def phys_activate() 
  subroutine phys_activate()
    call phys_init()
  end subroutine phys_activate
#:enddef
  
#:def phys_init()
    !> Initialize the module
  subroutine phys_init()
    use mod_global_parameters
!    use mod_particles, only: particles_init

    call phys_units()
    call read_params(par_files)

    phys_energy  = ffhd_energy
    phys_total_energy  = ffhd_energy
    phys_internal_e = .false.
    phys_gamma = ffhd_gamma
    phys_partial_ionization=ffhd_partial_ionization
    !$acc update device(physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization)

    use_particles = ffhd_particles

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
    need_global_cs2max = .true.
    !$acc update device(q_)
    !$acc update device(need_global_cs2max)

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    if (ffhd_force_diagonal) then
       ! ensure corners are filled, otherwise divide by zero when getting primitives
       !  --> only for debug purposes
       phys_req_diagonal = .true.
    endif

    ! set number of variables which need update ghostcells
    nwgc=nwflux
    !$acc update device(nwgc)

! use cycle, needs to be dealt with:    
!    ! Initialize particles module
!    if (ffhd_particles) then
!       call particles_init()
!       phys_req_diagonal = .true.
!    end if

    nvector      = 0 ! No. vector vars
    !$acc update device(nvector)
    !$acc update device(phys_req_diagonal)

!! HERE: initialize radiative cooling module

    hypertc_kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
    !$acc update device(hypertc_kappa)

  end subroutine phys_init
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
    dx, gradT, qsourcesplit)
  !$acc routine seq
#:if defined('GRAVITY')
  use mod_usr, only: gravity_field
#:endif    
#:if defined('BFIELD')
  use mod_bfield, only: magnetic_field, magnetic_field_divergence
#:endif
  use mod_global_parameters, only: dt, cs2max_global
  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCT(nw_phys), wCTprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim),dx(1:ndim)
  real(dp), intent(in)     :: gradT(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit
  ! .. local ..
  integer                  :: idim
  real(dp)                 :: field, Bfield, divBfield
  real(dp)                 :: Te,tau,htc_qrsc,sigT,invdx,taumin

#:if defined('GRAVITY')
  do idim = 1, ndim
     field  = gravity_field(wCT, x, idim)
     Bfield = magnetic_field(x, idim)
     wnew(iw_mom(1)) = wnew(iw_mom(1)) + qdt * field * wCT(iw_rho) * Bfield
     wnew(iw_e)      = wnew(iw_e) + qdt * field * wCT(iw_mom(1)) * Bfield
  end do
#:endif
  
#:if defined('BFIELD')
  divBfield       = magnetic_field_divergence(x)
  wnew(iw_mom(1)) = wnew(iw_mom(1)) + qdt * wCTprim(iw_e) * divBfield

  Te     = wCTprim(iw_e)/wCT(iw_rho)
  sigT   = hypertc_kappa*sqrt(Te**5)
  taumin = 0.4d0

  tau = taumin
  tau = max( taumin*dt, sigT*Te*(ffhd_gamma-1.0d0)/wCTprim(iw_e)/cs2max_global )

  htc_qrsc=0.0d0
  do idim = 1, ndim
    Bfield   = magnetic_field(x, idim)
    invdx    = 1.0d0/dx(idim)
    htc_qrsc = htc_qrsc+sigT*Bfield*gradT(idim)*invdx
  enddo
  htc_qrsc   = (htc_qrsc+wCT(iw_q))/tau

  wnew(iw_q) = wnew(iw_q)-qdt*htc_qrsc
#:endif  

end subroutine addsource_local
#:enddef

#:def to_primitive()
pure subroutine to_primitive(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(nw_phys)

  u(iw_mom(1)) = u(iw_mom(1))/u(iw_rho)
  
  u(iw_e) = (ffhd_gamma-1.0_dp) * (u(iw_e) - 0.5_dp * u(iw_rho) * &
     u(iw_mom(1))**2 )

end subroutine to_primitive
#:enddef

#:def to_conservative()  
pure subroutine to_conservative(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(nw_phys)
  real(dp)                :: inv_gamma_m1

  inv_gamma_m1 = 1.0_dp/(ffhd_gamma - 1.0_dp)

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
  use mod_bfield, only: magnetic_field
#:endif    
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: xC(1:ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: flux(nw_phys)
  real(dp)              :: inv_gamma_m1
  real(dp)              :: Bfield

  inv_gamma_m1 = 1.0_dp/(ffhd_gamma - 1.0_dp)

  Bfield = 1.0d0
#:if defined('BFIELD')
  Bfield = magnetic_field(xC, flux_dim)
#:endif    

  ! Density flux
  flux(iw_rho) = u(iw_rho) * u(iw_mom(1)) * Bfield

  ! Momentum flux with pressure term
  flux(iw_mom(1)) = (u(iw_rho) * u(iw_mom(1))**2 + u(iw_e)) * Bfield
  
  ! Energy flux with hyperbolic conduction included
  flux(iw_e) = ( u(iw_mom(1)) * (u(iw_e) * inv_gamma_m1 + 0.5_dp * &
     u(iw_rho) * u(iw_mom(1))**2 + u(iw_e)) &
              + u(iw_q)) * Bfield

  ! heat flux
  flux(iw_q) = 0.0d0

end subroutine get_flux
#:enddef

#:def get_cmax()  
pure real(dp) function get_cmax(u, x, flux_dim) result(wC)
  !$acc routine seq
#:if defined('BFIELD')
  use mod_bfield, only: magnetic_field
#:endif    
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim
  real(dp)              :: Bfield

  Bfield = 1.0d0
#:if defined('BFIELD')
  Bfield = magnetic_field(x, flux_dim)
#:endif    
  
  wC = sqrt(ffhd_gamma * u(iw_e) / u(iw_rho)) + abs(u(iw_mom(1))*Bfield)

end function get_cmax
#:enddef

#:def get_cs2()
!> obtain the squared sound speed
pure real(dp) function get_cs2(u) result(cs2)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)

  cs2 = ffhd_gamma * u(iw_e)/u(iw_rho)
  
end function get_cs2
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
