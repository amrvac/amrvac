#:if PHYS == 'mhd'
  
#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter, public              :: nw_phys=2+2*ndim+1
  
  !> Whether an energy equation is used
  logical, public                         :: mhd_energy = .true.
  !$acc declare copyin(mhd_energy)

  !> Index of the density (in the w array)
  integer, public                         :: rho_
  !$acc declare create(rho_)

  !> Indices of the momentum density
  integer, allocatable, public            :: mom(:)
  !$acc declare create(mom)

  !> Index of the energy density (-1 if not present)
  integer, public                         :: e_
  !$acc declare create(e_)

  !> Indices of the magnetic field
  integer, allocatable, public            :: mag(:)
  !$acc declare create(mag)

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public                         :: p_
  !$acc declare create(p_)

  !> Indices of the GLM psi
  integer, public :: psi_
  !$acc declare create(psi_)

  !> The adiabatic index
  double precision, public                :: mhd_gamma = 5.d0/3.0d0
  !$acc declare copyin(mhd_gamma)

  !> The adiabatic index minus 1
  double precision, public                :: mhd_gamma_m1
  !$acc declare copyin(mhd_gamma_m1)

  !> Helium abundance over Hydrogen
  double precision, public  :: He_abundance=0.1d0
  !$acc declare copyin(He_abundance)
  !> Ionization fraction of H
  !> H_ion_fr = H+/(H+ + H)
  double precision, public  :: H_ion_fr=1d0
  !$acc declare copyin(H_ion_fr)
  !> Ionization fraction of He
  !> He_ion_fr = (He2+ + He+)/(He2+ + He+ + He)
  double precision, public  :: He_ion_fr=1d0
  !$acc declare copyin(He_ion_fr)
  !> Ratio of number He2+ / number He+ + He2+
  !> He_ion_fr2 = He2+/(He2+ + He+)
  double precision, public  :: He_ion_fr2=1d0
  !$acc declare copyin(He_ion_fr2)
  ! used for eq of state when it is not defined by units,
  ! the units do not contain terms related to ionization fraction
  ! and it is p = RR * rho * T
  double precision, public  :: RR=1d0
  !$acc declare copyin(RR)

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: mhd_glm_alpha = 0.5d0
  !$acc declare copyin(mhd_glm_alpha)

  !> Whether to use gravity
  logical, public                         :: mhd_gravity = .false.
  !$acc declare copyin(mhd_gravity)

  !> Whether plasma is partially ionized
  logical, public                         :: mhd_partial_ionization = .false.
  !$acc declare copyin(mhd_partial_ionization)

  !> Whether particles module is added
  logical, public                         :: mhd_particles = .false.
  !$acc declare copyin(mhd_particles)

#:enddef

#:def read_params()
    !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mhd_list/ mhd_gamma, mhd_glm_alpha, mhd_gravity

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mhd_list, end=111)
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
    double precision :: mp, kB, miu0
    double precision :: a,b

    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi ! G^2 cm^2 dyne^-1
    end if
    a=1d0+4d0*He_abundance
    if(mhd_partial_ionization) then
      b=1d0+H_ion_fr+He_abundance*(He_ion_fr*(He_ion_fr2+1d0)+1d0)
    else
      b=2d0+3d0*He_abundance
    end if
    RR=1d0
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

    phys_energy  = mhd_energy
    phys_total_energy  = mhd_energy
    phys_internal_e = .false.
    phys_gamma = mhd_gamma
    phys_partial_ionization=mhd_partial_ionization
    need_global_cmax=.true.
    mhd_gamma_m1=mhd_gamma-1.0_dp
 !$acc update device(physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization,need_global_cmax,mhd_gamma_m1)

    use_particles = mhd_particles

    ! Determine flux variables
    rho_ = var_set_rho()
    !$acc update device(rho_)

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)
    !$acc update device(mom)

    ! Set index of energy variable
    if (mhd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if
    !$acc update device(e_,p_)

    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)
    !$acc update device(mag)

    psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    !$acc update device(psi_)

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .true.

    ! set number of variables which need update ghostcells
    nwgc=nwflux
    !$acc update device(nwgc)

! use cycle, needs to be dealt with:    
!    ! Initialize particles module
!    if (mhd_particles) then
!       call particles_init()
!       phys_req_diagonal = .true.
!    end if

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
subroutine addsource_local(qdt, dtfactor, qtC, wCT, wCTprim, qt, wnew, x, dr, &
    qsourcesplit)
  !$acc routine seq
  use mod_global_parameters, only:cmax_global
#:if defined('GRAVITY')
  use mod_usr, only: gravity_field
#:endif
  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCT(nw_phys), wCTprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim), dr(ndim)
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
  wnew(psi_)=wnew(psi_)*dexp(-qdt*cmax_global*mhd_glm_alpha/minval(dr))

end subroutine addsource_local
#:enddef

#:def addsource_nonlocal()
subroutine addsource_nonlocal(qdt, dtfactor, qtC, wCTprim, qt, wnew, x, dx, idir, &
     qsourcesplit)
  !$acc routine seq
  use mod_usr, only: bfield
  use mod_global_parameters, only: dt, cs2max_global

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCTprim(nw_phys,5)
  real(dp), intent(in)     :: x(1:ndim), dx(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  integer, intent(in)      :: idir
  logical, intent(in)      :: qsourcesplit
  ! .. local ..
  real(dp)                 :: field, mag

end subroutine addsource_nonlocal
#:enddef

#:def to_primitive()
  pure subroutine to_primitive(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(nw_phys)

    u(iw_mom(1))=u(iw_mom(1))/u(iw_rho)
    u(iw_mom(2))=u(iw_mom(2))/u(iw_rho)
    u(iw_mom(3))=u(iw_mom(3))/u(iw_rho)
    u(iw_e)=mhd_gamma_m1*(u(iw_e)-0.5_dp*&
      (u(iw_rho)*(u(iw_mom(1))**2+u(iw_mom(2))**2+u(iw_mom(3))**2)+&
       u(iw_mag(1))**2+u(iw_mag(2))**2+u(iw_mag(3))**2))

  end subroutine to_primitive
#:enddef

#:def to_conservative()  
  pure subroutine to_conservative(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(nw_phys)

    ! Compute energy from pressure and kinetic energy
    u(iw_e)=u(iw_e)/mhd_gamma_m1+0.5_dp*&
      (u(iw_rho)*(u(iw_mom(1))**2+u(iw_mom(2))**2+u(iw_mom(3))**2)+&
       u(iw_mag(1))**2+u(iw_mag(2))**2+u(iw_mag(3))**2)
    ! Compute momentum from density and velocity components
    u(iw_mom(1))=u(iw_rho)*u(iw_mom(1))
    u(iw_mom(2))=u(iw_rho)*u(iw_mom(2))
    u(iw_mom(3))=u(iw_rho)*u(iw_mom(3))

  end subroutine to_conservative
#:enddef

#:def get_flux()
  subroutine get_flux(u, xC, flux_dim, flux)
    use mod_global_parameters, only:cmax_global
    !$acc routine seq
    real(dp), intent(in)  :: u(nw_phys)
    real(dp), intent(in)  :: xC(1:ndim)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(nw_phys)
    real(dp)              :: ptotal

    ! Density flux
    flux(iw_rho)=u(iw_rho)*u(iw_mom(flux_dim))
    ! Momentum flux with pressure term
    flux(iw_mom(1))=u(iw_rho)*u(iw_mom(1))*u(iw_mom(flux_dim))-&
      u(iw_mag(flux_dim))*u(iw_mag(1))
    flux(iw_mom(2))=u(iw_rho)*u(iw_mom(2))*u(iw_mom(flux_dim))-&
      u(iw_mag(flux_dim))*u(iw_mag(2))
    flux(iw_mom(3))=u(iw_rho)*u(iw_mom(3))*u(iw_mom(flux_dim))-&
      u(iw_mag(flux_dim))*u(iw_mag(3))
    ptotal=u(iw_e)+0.5_dp*(u(iw_mag(1))**2+u(iw_mag(2))**2+u(iw_mag(3))**2)
    flux(iw_mom(flux_dim))=flux(iw_mom(flux_dim))+ptotal
    ! Energy flux
    flux(iw_e)=u(iw_mom(flux_dim))*(u(iw_e)/mhd_gamma_m1+0.5_dp*&
      u(iw_rho)*(u(iw_mom(1))**2+u(iw_mom(2))**2+u(iw_mom(3))**2)+&
      2.0_dp*ptotal-u(iw_e))-u(iw_mag(flux_dim))*&
      (u(iw_mag(1))*u(iw_mom(1))+u(iw_mag(2))*u(iw_mom(2))+u(iw_mag(3))*u(iw_mom(3)))
    ! Magnetic flux
    flux(iw_mag(1))=u(iw_mom(flux_dim))*u(iw_mag(1))-u(iw_mag(flux_dim))*u(iw_mom(1))
    flux(iw_mag(2))=u(iw_mom(flux_dim))*u(iw_mag(2))-u(iw_mag(flux_dim))*u(iw_mom(2))
    flux(iw_mag(3))=u(iw_mom(flux_dim))*u(iw_mag(3))-u(iw_mag(flux_dim))*u(iw_mom(3))
    ! GLM psi flux
    flux(iw_mag(flux_dim))=u(psi_)
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
    flux(psi_)=cmax_global**2*u(iw_mag(flux_dim))

  end subroutine get_flux
#:enddef

#:def get_cmax()  
pure real(dp) function get_cmax(u, x, flux_dim) result(wC)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim

  real(dp) :: inv_rho, b2, cfast2

  inv_rho=1.0_dp/u(iw_rho)
  wC=mhd_gamma*u(iw_e)*inv_rho
  cfast2=(u(iw_mag(1))**2+u(iw_mag(2))**2+u(iw_mag(3))**2)*inv_rho+wC
  b2=cfast2**2-4.0_dp*wC*u(iw_mag(flux_dim))**2*inv_rho
  if(b2<0.d0) b2=0.d0
  wC=sqrt(0.5_dp*(cfast2+sqrt(b2)))+abs(u(iw_mom(flux_dim)))

end function get_cmax
#:enddef  

#:def get_cs2()
!> obtain the squared sound speed
pure real(dp) function get_cs2(u) result(cs2)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)

  cs2 = phys_gamma*u(iw_e)/u(iw_rho)
end function get_cs2
#:enddef

#:endif
