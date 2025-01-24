!> User module for launching a 2-D line-driven wind from a rotating OB-star on
!> a 2.5-D grid. The radiation force follows Castor, Abbott, & Klein (1975)
!> using the the mod_cak_force module of MPI-AMRVAC. This user module is the
!> 2-D spherical analogue of the 1-D test problem.
!>
!> Two basic wind physics options can be explored depending on the setup:
!>   1. A 1-D radial force -> allows to check the wind-compressed disc (WCD)
!>      paradigm. See models by Owocki, Cranmer, & Blondin (1994), ApJ, 424.
!>   2. A 3-D vector force -> allows to check the radiative torque exerted by
!>      the azimuthal line-force component (Gayley & Owocki 2000, ApJ, 537).
!>
!> In both cases it is important that the rotation parameter W>=0.5 to see the
!> effect.
!> NOTE: at high rotational speeds the star becomes oblate, an effect which is
!>       currently not taken into account. It can be incorporated by modifying
!>       the lower BC as in Owocki et al. (1994), sect 2.2., or the PhD thesis
!>       of P. Petrenz (1999), sect 6.3.2
!>       BUT this will not work in combination with the vector force (for now).
!>       The vector force radiation ray positions are currently only specified
!>       for a spherical star.
!>
!> Coded up by Florian Driessen (2022)
module mod_usr

  ! Include a physics module
  use mod_hd
  
  ! Get access to some CAK radiation functionality
  use mod_cak_force, only: set_cak_force_norm, cak_alpha, gayley_qbar

  implicit none

  ! User input parameters
  real(8) :: mstar, rstar, twind, rhobound, beta, Wrot

  ! Extra parameters required in computation
  real(8) :: asound, Ggrav, vrot

  ! Dimensionless variables needed throughout computations
  real(8) :: dmstar, drstar, drhobound, dmdot, dvinf, dasound, dclight
  real(8) :: dGgrav, dvrot

contains

  !=============================================================================
  ! This routine should set user methods, and activate the physics module
  !=============================================================================
  subroutine usr_init()

    call set_coordinate_system("spherical_2.5D")
    call usr_params_read(par_files)

    ! Choose independent normalisation units:
    unit_length      = rstar    ! cm
    unit_temperature = twind    ! K
    unit_density     = rhobound ! g cm^-3

    
    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initial_conditions
    usr_special_bc     => special_bound
    usr_gravity        => stellar_gravity

    ! For cgs output in convert stage if wished for
    ! if (convert .and. saveprim) then
    !   w_convert_factor(rho_)   = unit_density
    !   w_convert_factor(mom(:)) = unit_velocity
    !   w_convert_factor(gcak1_) = unit_length / unit_time**2.0d0
    !   w_convert_factor(gcak2_) = unit_length / unit_time**2.0d0
    !   w_convert_factor(gcak3_) = unit_length / unit_time**2.0d0
    !   if (hd_energy) w_convert_factor(p_) = unit_pressure
    ! endif
    call HD_activate()
    call set_cak_force_norm(rstar,twind)
    
  end subroutine usr_init

  !=============================================================================
  ! Read in the usr.par file with the problem specific list
  !=============================================================================
  subroutine usr_params_read(files)

    ! Subroutine argument
    character(len=*), intent(in) :: files(:)

    ! Local variable
    integer :: n

    namelist /star_list/ mstar, rstar, twind, rhobound, beta, Wrot

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    end do

    ! Scale to cgs units
    mstar = mstar * const_MSun
    rstar = rstar * const_RSun

  end subroutine usr_params_read

  !=============================================================================
  ! Compute some quantities of interest (in CGS) before making unitless
  !=============================================================================
  subroutine initglobaldata_usr
    real(8) :: lstar, mumol, vesc, gammae, logg, logge, heff, vrotc, vinf, mdot

    ! Stellar structure
    lstar  = 4.0d0*dpi * rstar**2.0d0 * const_sigma * twind**4.0d0
    gammae = const_kappae * lstar/(4.0d0*dpi * const_G * mstar * const_c)
    logg   = log10(const_G * mstar/rstar**2.0d0)
    logge  = logg + log10(1.0d0 - gammae)
    mumol  = (1.0d0 + 4.0d0*He_abundance)/(2.0d0 + 3.0d0*He_abundance)
    asound = sqrt(twind * kB_cgs/(mumol * mp_cgs))
    heff   = asound**2.0d0 / 10.0d0**logge
    vrotc  = sqrt(const_G * mstar*(1.0d0 - gammae)/rstar)
    vrot   = vrotc * Wrot

    ! Wind quantities in CAK theory
    vesc   = sqrt(2.0d0 * const_G * mstar * (1.0d0 - gammae)/rstar)
    vinf   = vesc * sqrt(cak_alpha/(1.0d0 - cak_alpha))
    mdot   = lstar/const_c**2.0d0 * cak_alpha/(1.0d0 - cak_alpha) &
             * (gayley_qbar * gammae/(1.0d0 - gammae))**((1.0d0 - cak_alpha)/cak_alpha)
    
    if (mype == 0) then
      print*, '======================'
      print*, '   Unity quantities   '
      print*, '======================'
      print*, 'unit length        = ', unit_length
      print*, 'unit density       = ', unit_density
      print*, 'unit velocity      = ', unit_velocity
      print*, 'unit numberdensity = ', unit_numberdensity
      print*, 'unit pressure      = ', unit_pressure
      print*, 'unit temperature   = ', unit_temperature
      print*, 'unit time          = ', unit_time
      print*, '==============================================='
      print*, '   Stellar and wind parameters in CGS units    '
      print*, '==============================================='
      print*, 'L/Lsun                 = ', lstar/const_LSun
      print*, 'M/Msun                 = ', mstar/const_MSun
      print*, 'R/Rsun                 = ', rstar/const_RSun
      print*, 'Twind                  = ', twind
      print*, 'Mean molecular weight  = ', mumol
      print*, 'log(g)                 = ', logg
      print*, 'eff. log(g)            = ', logge
      print*, 'heff/Rstar             = ', heff/rstar
      print*, 'W (vrot/vrotc)         = ', Wrot
      print*, 'critical vrot          = ', vrotc
      print*, 'vrot                   = ', vrot
      print*, 'Eddington gamma        = ', gammae
      print*, 'adiabatic gamma        = ', hd_gamma
      print*, 'isothermal asound      = ', asound
      print*, 'eff. vesc              = ', vesc
      print*, 'CAK vinf               = ', vinf
      print*, 'FD vinf                = ', 3.0d0 * vesc
      print*, 'analytic Mdot CAK      = ', mdot * (const_years/const_MSun)
      print*, '... with FD correction = ', mdot/(1.0d0 + cak_alpha)**(1.0d0/cak_alpha) * (const_years/const_MSun)
    endif

    ! Make some variables needed for computations dimensionless
    dmstar    = mstar/(unit_density * unit_length**3.0d0)
    drstar    = rstar/unit_length
    dasound   = asound/unit_velocity
    drhobound = rhobound/unit_density
    dvinf     = vinf/unit_velocity
    dmdot     = mdot * unit_time/(unit_density * unit_length**3.0d0)
    dGgrav    = const_G * unit_density * unit_time**2.0d0
    dvrot     = vrot/unit_velocity
    
    if (mype == 0) then
      print*, '========================================'
      print*, '  Dimensionless computation quantities  '
      print*, '========================================'
      print*, 'Mstar    = ', dmstar
      print*, 'Rstar    = ', drstar
      print*, 'rhobound = ', drhobound
      print*, 'Mdot CAK = ', dmdot
      print*, 'asound   = ', dasound
      print*, 'vinf     = ', dvinf
      print*, 'vrot     = ', dvrot
      print*, 'Ggrav    = ', dGgrav
    endif

  end subroutine initglobaldata_usr
  
  !=============================================================================
  ! Initial conditions start from spherically symmetric 1-D CAK beta-law wind
  !=============================================================================
  subroutine initial_conditions(ixI^L,ixO^L,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    real(8) :: sfac

    ! Small offset (~vtherm/vinf) to avoid starting at terminal wind speed
    sfac = 1.0d0 - 1.0d-3**(1.0d0/beta)

    w(ixO^S,mom(1)) = dvinf * ( 1.0d0 - sfac * drstar/x(ixO^S,1) )**beta
    w(ixO^S,rho_)   = dmdot / (4.0d0*dpi * x(ixO^S,1)**2.0d0 * w(ixO^S,mom(1)))
    w(ixO^S,mom(2)) = 0.0d0
    
    if (hd_rotating_frame) then
      w(ixO^S,mom(3)) = 0.0d0
    else
      ! Angular momentum conserving
      w(ixO^S,mom(3)) = dvrot * sin(x(ixO^S,2)) * drstar**2.0d0/x(ixO^S,1)
    endif

    if (hd_energy) then
      w(ixO^S,p_) = dasound**2.0 * drhobound * (w(ixO^S,rho_)/drhobound)**hd_gamma
    endif

    call hd_to_conserved(ixI^L,ixO^L,w,x)

    ! Initialise extra vars at 0
    w(ixO^S,nw-nwextra+1:nw) = 0.0d0

  end subroutine initial_conditions

  !=============================================================================
  ! Special user boundary conditions at inner radial boundary:
  !   vr (extrapolated); rho, vtheta, vphi (fixed)
  !=============================================================================
  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: ir

    select case (iB)
    case(1)

      w(ixB^S,rho_) = drhobound

      ! Radial velocity (constant slope extrapolation)
      w(ixBmax1+1^%1ixB^S,mom(1))=w(ixBmax1+1^%1ixB^S,mom(1))/w(ixBmax1+1^%1ixB^S,rho_)
      w(ixBmax1+2^%1ixB^S,mom(1))=w(ixBmax1+2^%1ixB^S,mom(1))/w(ixBmax1+2^%1ixB^S,rho_)
      do ir = ixBmax1,ixBmin1,-1
        w(ir^%1ixB^S,mom(1)) = w(ir+1^%1ixB^S,mom(1)) &
                               - (w(ir+2^%1ixB^S,mom(1)) - w(ir+1^%1ixB^S,mom(1))) &
                               * (x(ir+1^%1ixB^S,1) - x(ir^%1ixB^S,1))&
                                 /(x(ir+2^%1ixB^S,1) - x(ir+1^%1ixB^S,1))
      enddo

      ! Polar velocity (no-slip condition to avoid equator-ward flow)
      w(ixB^S,mom(2)) = 0.0d0

      if (hd_rotating_frame) then
        w(ixB^S,mom(3)) = 0.0d0
      else
        ! Rigid body rotation of star
        w(ixB^S,mom(3)) = dvrot * sin(x(ixB^S,2))
      endif

      ! Avoid supersonic ghost cells, also avoid overloading too much
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)), dasound)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -dasound)

      if (hd_energy) then
        w(ixB^S,p_) = dasound**2.0 * drhobound * (w(ixB^S,rho_)/drhobound)**hd_gamma
      endif
      w(ixBmax1+1^%1ixB^S,mom(1))=w(ixBmax1+1^%1ixB^S,mom(1))*w(ixBmax1+1^%1ixB^S,rho_)
      w(ixBmax1+2^%1ixB^S,mom(1))=w(ixBmax1+2^%1ixB^S,mom(1))*w(ixBmax1+2^%1ixB^S,rho_)

      call hd_to_conserved(ixI^L,ixB^L,w,x)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

  !=============================================================================
  ! Compute stellar gravity
  !=============================================================================
  subroutine stellar_gravity(ixI^L,ixO^L,wCT,x,gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    real(8), intent(in)  :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    real(8), intent(out) :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:) = 0.0d0

    ! Only in radial direction
    gravity_field(ixO^S,1) = -dGgrav * dmstar/x(ixO^S,1)**2.0d0

  end subroutine stellar_gravity

end module mod_usr
