!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rmhd
  use mod_fld

  implicit none

  double precision :: rho1
  double precision :: rho2
  double precision :: v1
  double precision :: v2
  double precision :: T1
  double precision :: T2
  double precision :: B1 
  double precision :: B2 

  double precision :: vy1, vz1, By1, Bz1 
  double precision :: vy2, vz2, By2, Bz2 

  double precision :: p1,p2,eg1,eg2,Er1,Er2

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    {^IFONED call set_coordinate_system("Cartesian_1.75D")}
    {^IFTWOD call set_coordinate_system("Cartesian_2D")}
    {^IFTHREED call set_coordinate_system("Cartesian_3D")}

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Manually refine grid near shock
    usr_refine_grid => refine_shock

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! Boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_init_vector_potential=>initvecpot_usr

    ! Active the physics module
    call rmhd_activate()
  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    call params_read(par_files)

    p1 = const_kB*T1*rho1/(const_mp*fld_mu)
    p2 = const_kB*T2*rho2/(const_mp*fld_mu)

    eg1 = p1/(rmhd_gamma-1.d0) + half*rho1*(v1*v1 + vy1*vy1 + vz1*vz1) + half*(B1*B1 + By1*By1 + Bz1*Bz1)/(4.0*dpi)
    eg2 = p2/(rmhd_gamma-1.d0) + half*rho2*(v2*v2 + vy2*vy2 + vz2*vz2) + half*(B2*B2 + By2*By2 + Bz2*Bz2)/(4.0*dpi)

    Er1 = const_rad_a*T1**4.d0
    Er2 = const_rad_a*T2**4.d0

    print*, 'M_1: ', v1/dsqrt(rmhd_gamma*p1/rho1) 
    print*, 'M_2: ', v2/dsqrt(rmhd_gamma*p2/rho2) 

    print*, 'RHD-quantity: ', 'Left', ' | ', 'Right'
    print*, 'density', rho1, ' | ', rho2
    print*, 'velocity', v1,vy1,vz2, ' | ', v2,vy2,vz2 
    print*, 'momentum', rho1*v1, ' | ', rho2*v2
    print*, 'gas pressure', p1, ' | ', p2
    print*, 'gas energy', eg1, ' | ', eg2
    print*, 'radiation energy', Er1, ' | ', Er2
    print*, 'magnetic field', B1,By1,Bz1, ' | ', B2,By2,Bz2 
    print*, 'plasma beta', p1*(8.0*dpi)/(B1*B1 + By1*By1 + Bz1*Bz1), ' | ', p2*(8.0*dpi)/(B2*B2 + By2*By2 + Bz2*Bz2) 

    unit_velocity = v1 
    unit_numberdensity = rho1/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = 1.d5

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity
    unit_magneticfield = sqrt(4.0*dpi * unit_pressure); 

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    print*, 'unit_numberdensity', unit_numberdensity
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, 'unit_v', unit_velocity
    print*, 'unit_p', unit_pressure
    print*, 'unit_magneticfield', unit_magneticfield 
    print*, 'unit_opacity', unit_opacity
    print*, 'unit_radflux', unit_radflux

    rho1 = rho1/unit_density
    rho2 = rho2/unit_density

    v1 = v1/unit_velocity
    v2 = v2/unit_velocity

    T1 = T1/unit_temperature
    T2 = T2/unit_temperature

    p1 = p1/unit_pressure
    p2 = p2/unit_pressure

    eg1 = eg1/unit_pressure
    eg2 = eg2/unit_pressure

    Er1 = Er1/unit_pressure
    Er2 = Er2/unit_pressure

    B1 = B1/unit_magneticfield 
    B2 = B2/unit_magneticfield 

    vy1 = vy1/unit_velocity 
    vy2 = vy2/unit_velocity 
    vz1 = vz1/unit_velocity 
    vz2 = vz2/unit_velocity 

    By1 = By1/unit_magneticfield 
    By2 = By2/unit_magneticfield 
    Bz1 = Bz1/unit_magneticfield 
    Bz2 = Bz2/unit_magneticfield 

    print*, 'RHD-fluxes: ', 'Left', ' | ', 'Right'
    print*, 'density', rho1*v1, ' | ', rho2*v2
    print*, 'momentum', (rho1*v1*v1 + p1 + Er1/3 + half*(B1*B1 + By1*By1 + Bz1*Bz1)), ' | ',&
    (rho2*v2*v2 + p2 + Er2/3 + half*(B2*B2 + By2*By2 + Bz2*Bz2)) 
    print*, 'gas energy', p1*v1 + eg1*v1, ' | ', p2*v2 + eg2*v2
    print*, 'radiation energy', Er1*v1, ' | ', Er2*v2
    print*, 'total energy', (p1+eg1+Er1)*v1, ' | ', (p2+eg2+Er2)*v2
    print*, 'energy jump condition?', ((p1+eg1+4.0*Er1/3.0 + half*(B1*B1+By1*By1+Bz1*Bz1))*v1 - B1*(v1*B1+vy1*By1+vz1*Bz1)), ' | ',&
    ((p2+eg2+4.0*Er2/3.0 + half*(B2*B2+By2*By2+Bz2*Bz2))*v2 - B2*(v2*B2+vy2*By2+vz2*Bz2))
    print*, 'momentum_y', (rho1*v1*vy1 - B1*By1), ' | ', (rho2*v2*vy2 - B2*By2)
    print*, 'momentum_z', (rho1*v1*vz1 - B1*Bz1), ' | ', (rho2*v2*vz2 - B2*Bz2) 
    print*, 'magfield_y', (By1*v1 - B1*vy1), ' | ', (By2*v2 - B2*vy2) 
    print*, 'magfield_y', (Bz1*v1 - B1*vz1), ' | ', (Bz2*v2 - B2*vz2) 
    print*, 'plasma beta', p1*2.0/(B1*B1 + By1*By1 + Bz1*Bz1), ' | ', p2*2.0/(B2*B2 + By2*By2 + Bz2*Bz2) 

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /shock_list/ rho1, rho2, v1, v2, T1, T2, B1, B2, vy1, vz1, By1, Bz1, vy2, vz2, By2, Bz2

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, shock_list, end=113)
113    close(unitpar)
    end do
  end subroutine params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rho(ixI^S), Temp(ixI^S), press(ixI^S), vel(ixI^S)

    double precision :: kappa(ixO^S), fld_R(ixO^S), lambda(ixO^S)

    w(ixI^S,rho_) = rho1
    w(ixI^S,mom(:)) = 0.d0
    w(ixI^S,mom(1)) = rho1*v1
    w(ixI^S,mom(2)) = rho1*vy1
    w(ixI^S,mom(3)) = rho1*vz1
    w(ixI^S,e_) = eg1
    w(ixI^S,r_e) = Er1

    w(ixI^S,mag(:)) = 0.d0
    w(ixI^S,mag(1)) = B1
    w(ixI^S,mag(2)) = By1 
    w(ixI^S,mag(3)) = Bz1 

    where (x(ixI^S,1) .gt. 0.d0)
      w(ixI^S,rho_) = rho2
      w(ixI^S,mom(1)) = rho2*v2
      w(ixI^S,mom(2)) = rho2*vy2
      w(ixI^S,mom(3)) = rho2*vz2
      w(ixI^S,e_) = eg2
      w(ixI^S,r_e) = Er2
      w(ixI^S,mag(1)) = B2
      w(ixI^S,mag(2)) = By2
      w(ixI^S,mag(3)) = Bz2
    end where
  end subroutine initial_conditions


  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_fld


    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ii

    select case (iB)

    case(1)
      w(ixB^S,rho_) = rho1
      w(ixB^S,mom(:)) = 0.d0
      w(ixB^S,mom(1)) = rho1*v1
      w(ixB^S,mom(2)) = rho1*vy1 
      w(ixB^S,mom(3)) = rho1*vz1 
      w(ixB^S,e_) = eg1
      w(ixB^S,r_e) = Er1
      w(ixB^S,mag(:)) = 0.0d0
      w(ixB^S,mag(1)) = B1
      w(ixB^S,mag(2)) = By1
      w(ixB^S,mag(3)) = Bz1

    case(2)
      w(ixB^S,rho_) = rho2
      w(ixB^S,mom(:)) = 0.d0
      w(ixB^S,mom(1)) = rho2*v2
      w(ixB^S,mom(2)) = rho2*vy2
      w(ixB^S,mom(3)) = rho2*vz2
      w(ixB^S,e_) = eg2
      w(ixB^S,r_e) = Er2
      w(ixB^S,mag(:)) = 0.0d0
      w(ixB^S,mag(1)) = B2
      w(ixB^S,mag(2)) = By2 
      w(ixB^S,mag(3)) = Bz2 

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
        mg%bc(iB, mg_iphi)%bc_value = Er1
    case (2)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = Er2

      case default
        print *, "Not a standard: ", typeboundary(r_e, iB)
        error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions

  subroutine refine_shock(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    !> Refine close to base
    coarsen = -1
    refine = -1

    if (it .gt. slowsteps) then
      if (any(x(ixG^S,1) < 2.d-1 .and. x(ixG^S,1) > -2.d-1)) refine=1
    endif

    if (it .gt. slowsteps) then
      if (any(x(ixG^S,1) < 1.d-1 .and. x(ixG^S,1) > -1.d-1)) refine=1
    endif
  end subroutine refine_shock

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the corners
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    A(ixC^S) = zero
  end subroutine initvecpot_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_fld
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: lamb(ixO^S), R(ixO^S)

    call fld_get_fluxlimiter(w,x,ixI^L,ixO^L,lamb,R,1)
    w(ixO^S,nw+1)=lamb(ixO^S)
    w(ixO^S,nw+2)=R(ixO^S)
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  character(len=*) :: varnames
  varnames='Lambda R'

  end subroutine specialvarnames_output

end module mod_usr
