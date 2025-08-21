#:if PHYS == 'ffhd'
  
#:def phys_vars()

  !> Radiative cooling fluid
  #:if defined('COOLING')
    use mod_radiative_cooling, only: rc_fluid
  #:endif

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

  !> The adiabatic index
  double precision, public                :: ffhd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: ffhd_adiab = 1.0d0

  !> The thermal conductivity kappa in hyperbolic thermal conduction
  double precision, public                :: hypertc_kappa
  !$acc declare copyin(hypertc_kappa)

  !> Whether plasma is partially ionized
  logical, public                         :: ffhd_partial_ionization = .false.

  #:if defined('COOLING')
  !> Radiative cooling fluid
  type(rc_fluid), public, allocatable :: rc_fl
  !$acc declare create(rc_fl)
  #:endif

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
  
#:def phys_init()
    !> Initialize the module
  subroutine phys_init()

    #:if defined('COOLING')
    use mod_radiative_cooling, only: radiative_cooling_init_params, radiative_cooling_init
    #:endif

    use mod_global_parameters

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

    ! htc requried here

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
  real(dp)                 :: field

#:if defined('GRAVITY')
  do idim = 1, ndim
     field = gravity_field(wCT, x, idim)
     wnew(iw_mom(idim)) = wnew(iw_mom(idim)) + qdt * field * wCT(iw_rho)
     wnew(iw_e)         = wnew(iw_e) + qdt * field * wCT(iw_mom(idim))
  end do
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
subroutine get_flux(u, flux_dim, flux)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)
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
pure real(dp) function get_cmax(u, flux_dim) result(wC)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)
  integer, intent(in)   :: flux_dim

  
  wC = sqrt(hd_gamma * u(iw_e) / u(iw_rho)) + abs(u(iw_mom(flux_dim)))

end function get_cmax
#:enddef  
#:endif

#:def get_rho()
pure real(dp) function get_rho(w, x) result(rho)
  !$acc routine seq
  real(dp), intent(in)  :: w(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)

  rho = w(iw_rho)
end function get_rho
#:endif


#:def get_pthermal()
pure real(dp) function get_pthermal(w, x) result(pth)
  !$acc routine seq
  real(dp), intent(in)  :: w(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)

  @:get_kin_en()
  pth = phys_gamma*w(iw_e)-get_kin_en(w,x)
end function get_pthermal
#:endif

#:def get_Rfactor()
pure real(dp) function get_Rfactor() result(Rfactor)
  !$acc routine seq
  Rfactor = 1.0_dp
end function get_Rfactor
#:endif

#:def get_kin_en()
pure real(dp) function get_kin_en(w, x) result(kin_en)
  !$acc routine seq
  real(dp), intent(in)  :: w(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)

  kin_en = 0.5d0*w(iw_rho)*w(iw_mom(1))**2/w(iw_rho)
end function get_kin_en
#:endif