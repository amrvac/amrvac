!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: e_eq
  double precision :: rho0
  double precision :: t0
  double precision :: e0
  double precision :: E_r0

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_constants

    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init

!==========================================================================================

subroutine initglobaldata_usr
  use mod_global_parameters

  !Units
  unit_numberdensity = 1.d0/((1.d0+4.d0*He_abundance)*const_mp)
  unit_pressure = 1.d0
  unit_time = 1.d0

  unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
  unit_velocity = dsqrt(unit_pressure/unit_density)
  unit_length = unit_time*unit_velocity
  unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB)

  unit_radflux = unit_velocity*unit_pressure
  unit_opacity = one/(unit_density*unit_length)

  call usr_params_read(par_files)

  e_eq = (E_r0/(const_rad_a))**(1.d0/4.d0) &
        *one/(rhd_gamma-one)*const_kB*rho0 &
        /(fld_mu*const_mp)

end subroutine initglobaldata_usr


subroutine usr_params_read(files)
  use mod_global_parameters, only: unitpar
  use mod_fld
  use mod_constants
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /test_list/ rho0, E_r0, t0

  do n = 1, size(files)
     open(unitpar, file=trim(files(n)), status="old")
     read(unitpar, test_list, end=111)
     111    close(unitpar)
  end do


end subroutine usr_params_read

!==========================================================================================

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_constants
    use mod_rhd_phys, only: rhd_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, ndim)
    double precision, intent(inout) :: w(ixI^S, nw)

    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)


    ! Set initial values for w
    w(ixI^S, rho_) = rho0
    w(ixI^S, mom(:)) = zero
    w(ixI^S,r_e) = E_r0
    w(ixI^S, e_) = t0*e_eq

    print*, 'unit_time', unit_time
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, 'unit_numberdensity', unit_numberdensity
    print*, 'unit_velocity', unit_velocity
    print*, 'unit_pressure', unit_pressure
    print*, '================================================================'

  end subroutine initial_conditions

end module mod_usr

!==========================================================================================
