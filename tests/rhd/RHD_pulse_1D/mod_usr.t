!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: rho0 = 1.2d0
  double precision :: v0 = 5.d7
  double precision :: T0 = 1.d7
  double precision :: T1 = 2.d7
  double precision :: wdth = 24.d0

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    {^IFONED call set_coordinate_system("Cartesian_1D")}
    {^IFTWOD call set_coordinate_system("Cartesian_2D")}
    {^IFTHREED call set_coordinate_system("Cartesian_3D")}

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! Boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    unit_velocity = v0 !r_arr(nghostcells) ! cm
    unit_numberdensity = rho0/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = wdth

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    print*, unit_time, 's'

    rho0 = rho0/unit_density
    v0 = v0/unit_velocity
    T0 = T0/unit_temperature
    T1 = T1/unit_temperature
    wdth = wdth/unit_length

    print*, 'u_time', unit_time
    print*, 'u_length', unit_length
    print*, 'u_density', unit_density
    print*, 'u_pressure', unit_pressure

  end subroutine initglobaldata_usr


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: temp(ixI^S), pth(ixI^S)
    double precision :: kappa(ixO^S), fld_R(ixO^S), lambda(ixO^S)

    ! v0 = 0.d0

    temp(ixI^S) = T0 + (T1-T0)*dexp(-x(ixI^S,1)**2.d0/(2*wdth**2))

    w(ixI^S,rho_) = rho0*T0/temp(ixI^S) &
    + const_rad_a*fld_mu*const_mp/(3.d0*const_kB) &
    * unit_temperature**3/unit_density &
    * (T0**4.d0/temp(ixI^S) - temp(ixI^S)**3.d0)

    w(ixI^S,mom(:)) = 0.d0
    w(ixI^S,mom(1)) = w(ixI^S,rho_)*v0

    pth(ixI^S) = temp(ixI^S)*w(ixI^S,rho_)
    w(ixI^S,e_) = pth(ixI^S)/(rhd_gamma-1.d0) + half*w(ixI^S,rho_)*v0**2
    w(ixI^S,r_e) = const_rad_a*(temp(ixI^S)*unit_temperature)**4.d0/unit_pressure

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R, nghostcells)

    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_fld


    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: temp(ixB^S), pth(ixB^S)

    select case (iB)
    case(1,2)
      temp(ixB^S) = T0 + (T1-T0)*dexp(-(x(ixB^S,1)-v0*qt)**2.d0/(2*wdth**2))
      w(ixB^S,rho_) = rho0*T0/temp(ixB^S) &
      + const_rad_a*fld_mu*const_mp/(3.d0*const_kB) &
      * unit_temperature**3/unit_density &
      * (T0**4.d0/temp(ixB^S) - temp(ixB^S)**3.d0)
      w(ixB^S,mom(:)) = 0.d0
      w(ixB^S,mom(1)) = w(ixB^S,rho_)*v0
      pth(ixB^S) = temp(ixB^S)*w(ixB^S,rho_)
      w(ixB^S,e_) = pth(ixB^S)/(rhd_gamma-1.d0) + half*w(ixB^S,rho_)*v0**2
      w(ixB^S,r_e) = const_rad_a*(temp(ixB^S)*unit_temperature)**4.d0/unit_pressure

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions


  subroutine mg_boundary_conditions(iB)

    use mod_global_parameters
    use mod_multigrid_coupling

    integer, intent(in)             :: iB

    select case (iB)
    case (1)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = const_rad_a*(T0*unit_temperature)**4.d0/unit_pressure
    case (2)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = const_rad_a*(T0*unit_temperature)**4.d0/unit_pressure

      case default
        print *, "Not a standard: ", typeboundary(r_e, iB)
        error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions

end module mod_usr
