!> User module for launching a spherically symmetric line-driven wind from an
!> OB-star after the method of Castor, Abbott, & Klein (1975) using the
!> mod_cak_force module of MPI-AMRVAC.
!> The parameterisation of the line-ensemble distribution function follows the
!> line-strength formalism of Gayley (1995), including the possible cut-off.
!>
!> Runs for both isothermal and adiabatic winds. In the latter case we include
!> a floor temperature inside mod_cak_force to mimic stellar heating (and also
!> avoid negative thermal energies in the hypersonic flow regions).
!>
!> Options for wind to be specified in the cak_list of the .par file:
!>   cak_1d_opt = 0  : radially streaming CAK wind (1975)
!>   cak_1d_opt = 1  : finite disc corrected CAK wind
!>                     (Friend & Abbott 1986; Pauldrach et al. 1986)
!>   cak_1d_opt = 2  : finite disc CAK wind with line-strength cut-off
!>                     (Gayley 1995)
!>
!> Coded up by Florian Driessen (2019)
module mod_usr

  ! Include a physics module
  use mod_HD

  ! Get access to some CAK radiation functionality
  use mod_cak_force, only: set_cak_force_norm, cak_alpha, gayley_qbar, gcak1_

  implicit none

  ! User input parameters
  real(8) :: mstar_sol, rstar_sol, twind_cgs, rhobound_cgs, beta

  ! Extra parameters required in computation
  real(8) :: mstar, rstar, rhobound, mdot, vinf, asound, Ggrav

contains

  !============================================================================
  ! This routine should set user methods, and activate the physics module
  !============================================================================
  subroutine usr_init()

    {^IFONED call set_coordinate_system("spherical")}
    {^NOONED call mpistop("Modify problem yourself if you want a >1-D model")}
    call usr_params_read(par_files)

    ! Choose independent normalisation units:
    unit_length      = rstar_sol * const_RSun ! cm
    unit_temperature = twind_cgs              ! K
    unit_density     = rhobound_cgs           ! g cm^-3

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initial_conditions
    usr_special_bc     => special_bound
    usr_gravity        => stellar_gravity

    call HD_activate()

  end subroutine usr_init

  !============================================================================
  ! Read in the usr.par file with the problem specific list
  !============================================================================
  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /star_list/ mstar_sol, rstar_sol, twind_cgs, rhobound_cgs, beta

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine usr_params_read

  !============================================================================
  ! Compute some quantities of interest (in CGS) and make some unitless
  !============================================================================
  subroutine initglobaldata_usr
    real(8) :: mstar_cgs, rstar_cgs, lstar_cgs, mumol, vesc_cgs, gammae
    real(8) :: logg_cgs, logge_cgs, heff_cgs, vinf_cgs, asound_cgs, mdot_cgs
    real(8) :: pthbound, twind

    mstar_cgs = mstar_sol * const_MSun
    rstar_cgs = rstar_sol * const_RSun

    ! Stellar structure
    lstar_cgs  = 4.0d0*dpi * rstar_cgs**2.0d0 * const_sigma * twind_cgs**4.0d0
    gammae     = const_kappae * lstar_cgs / (4.0d0*dpi * const_G * mstar_cgs * const_c)
    logg_cgs   = log10(const_G * mstar_cgs/rstar_cgs**2.0d0)
    logge_cgs  = logg_cgs + log10(1.0d0 - gammae)
    mumol      = (1.0d0 + 4.0d0*He_abundance)/(2.0d0 + 3.0d0*He_abundance)
    asound_cgs = sqrt(twind_cgs * kB_cgs/(mumol * mp_cgs))
    heff_cgs   = asound_cgs**2.0d0 / 10.0d0**logge_cgs

    ! Wind quantities in CAK theory
    vesc_cgs = sqrt(2.0d0 * const_G * mstar_cgs * (1.0d0 - gammae)/rstar_cgs)
    vinf_cgs = vesc_cgs * sqrt(cak_alpha/(1.0d0 - cak_alpha))
    mdot_cgs = lstar_cgs/const_c**2.0d0 * cak_alpha/(1.0d0 - cak_alpha) &
               * (gayley_qbar * gammae/(1.0d0 - gammae))**((1.0d0 - cak_alpha)/cak_alpha)

    call set_cak_force_norm(rstar_cgs,twind_cgs)

    if (mype == 0 .and. .not.convert) then
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
      print*, 'L/Lsun                 = ', lstar_cgs/const_LSun
      print*, 'M/Msun                 = ', mstar_cgs/const_MSun
      print*, 'R/Rsun                 = ', rstar_cgs/const_RSun
      print*, 'Twind                  = ', twind_cgs
      print*, 'Mean molecular weight  = ', mumol
      print*, 'log(g)                 = ', logg_cgs
      print*, 'eff. log(g)            = ', logge_cgs
      print*, 'heff/Rstar             = ', heff_cgs/rstar_cgs
      print*, 'Eddington gamma        = ', gammae
      print*, 'isothermal asound      = ', asound_cgs
      print*, 'eff. vesc              = ', vesc_cgs
      print*, 'CAK vinf               = ', vinf_cgs
      print*, 'FD vinf                = ', 3.0d0 * vesc_cgs
      print*, 'analytic Mdot CAK      = ', mdot_cgs * (const_years/const_MSun)
      print*, '... with FD correction = ', mdot_cgs/(1.0d0 + cak_alpha)**(1.0d0/cak_alpha) * (const_years/const_MSun)
    endif

    ! Make some variables needed for computations dimensionless
    mstar    = mstar_cgs / (unit_density * unit_length**3.0d0)
    rstar    = rstar_cgs / unit_length
    twind    = twind_cgs / unit_temperature
    asound   = asound_cgs / unit_velocity
    rhobound = rhobound_cgs / unit_density
    pthbound = rhobound * twind
    vinf     = vinf_cgs / unit_velocity
    mdot     = mdot_cgs * unit_time / (unit_density * unit_length**3.0d0)
    Ggrav    = const_G * unit_density * unit_time**2.0d0

    ! Compute correct entropy
    hd_adiab = pthbound / rhobound**hd_gamma

    if (mype == 0 .and. .not.convert) then
      print*, '========================================'
      print*, '  Dimensionless computation quantities  '
      print*, '========================================'
      print*, 'Mstar    = ', mstar
      print*, 'Rstar    = ', rstar
      print*, 'Twind    = ', twind
      print*, 'rhobound = ', rhobound
      print*, 'Mdot CAK = ', mdot
      print*, 'asound   = ', asound
      print*, 'vinf     = ', vinf
      print*, 'Ggrav    = ', Ggrav
      print*, 'hd_gamma = ', hd_gamma
      print*, 'hd_adiab = ', hd_adiab
    endif

    ! For cgs output in convert stage if wished for
    ! if (convert .and. saveprim) then
    !   w_convert_factor(rho_)   = unit_density
    !   w_convert_factor(mom(:)) = unit_velocity
    !   w_convert_factor(gcak1_) = unit_length / unit_time**2.0d0
    !   if (hd_energy) w_convert_factor(p_) = unit_pressure
    !   length_convert_factor    = unit_length
    !   time_convert_factor      = unit_time
    ! endif

  end subroutine initglobaldata_usr

  !============================================================================
  ! Initial conditions start from spherical 1d beta law and mass conservation
  !============================================================================
  subroutine initial_conditions(ixI^L,ixO^L,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: sfac

    ! Small offset (~vtherm/vinf) to avoid starting at terminal wind speed
    sfac = 1.0d0 - 1.0d-3**(1.0d0/beta)

    w(ixO^S,mom(1)) = vinf * ( 1.0d0 - sfac * rstar / x(ixO^S,1) )**beta
    w(ixO^S,rho_)   = mdot / (4.0d0*dpi * x(ixO^S,1)**2.0d0 * w(ixO^S,mom(1)))

    ! Isothermal initial condition
    if (hd_energy) w(ixO^S,p_) = w(ixO^S,rho_)

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

  !============================================================================
  ! Special user boundary conditions at inner + outer radial boundary
  !============================================================================
  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: ir

    select case (iB)
    case(1)
      ! Mass density fixed at boundary mass density
      w(ixB^S,rho_) = rhobound

      ! Radial velocity field (constant slope extrapolation)
      do ir = ixBmax1,ixBmin1,-1
        if (ir == ixBmax1) then
          w(ir,mom(1))=w(ir+1,mom(1))/w(ir+1,rho_)-(w(ir+2,mom(1))/w(ir+2,rho_)-w(ir+1,mom(1))/w(ir+1,rho_)) &
                            *(x(ir+1,1)-x(ir,1))/(x(ir+2,1)-x(ir+1,1))
        else
          w(ir,mom(1))=w(ir+1,mom(1))
        end if
      end do

      ! Avoid supersonic ghost cells, also avoid overloading too much
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)), asound)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -asound)

      if (hd_energy) w(ixB^S,p_) = hd_adiab * w(ixB^S,rho_)**hd_gamma

      call hd_to_conserved(ixI^L,ixB^L,w,x)

    case(2)
      ! Constant extrapolation of all to have continuous velocity gradient
      w(ixB^S,rho_) = w(ixBmin1-1,rho_) * (x(ixBmin1-1,1) / x(ixB^S,1))**2.0d0

      w(ixBmin1,mom(1)) = w(ixBmin1-1,mom(1))/w(ixBmin1-1,rho_)+(w(ixBmin1-1,mom(1))/w(ixBmin1-1,rho_)-w(ixBmin1-2,mom(1))/w(ixBmin1-2,rho_))
      do ir = ixBmin1+1,ixBmax1
        w(ir,mom(1)) = w(ir-1,mom(1))+(w(ixBmin1-1,mom(1))/w(ixBmin1-1,rho_)-w(ixBmin1-2,mom(1))/w(ixBmin1-2,rho_))
      enddo

      if (hd_energy) w(ixB^S,p_) = hd_adiab * w(ixB^S,rho_)**hd_gamma

      call hd_to_conserved(ixI^L,ixB^L,w,x)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

  !============================================================================
  ! Compute stellar gravity
  !============================================================================
  subroutine stellar_gravity(ixI^L,ixO^L,wCT,x,gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    real(8), intent(in)  :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    real(8), intent(out) :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:) = 0.0d0

    ! Only in radial direction
    gravity_field(ixO^S,1) = -Ggrav * mstar/x(ixO^S,1)**2.0d0

  end subroutine stellar_gravity

end module mod_usr
