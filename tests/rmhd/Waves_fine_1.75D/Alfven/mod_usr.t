!> This is a template for a new user problem
module mod_usr
  use mod_rmhd 

  implicit none

  double precision :: rho0
  double precision :: eg0
  double precision :: tau_wave
  double precision :: ampl

  double precision :: T0, a0, p0, Er0
  double precision :: B0x, B0y, B0z 
  double precision :: p0_B

  double precision :: p0_by_Er0, p0_by_p0_B
  double precision :: va, vax, cslow, cfast

  double precision :: wvl, omega, wavenumber, tau_c, tau_a
  double precision :: Bo, energy_ratio, r_Bo, ca
  double precision :: A_rho, A_v, A_p, A_e, A_Er
  double precision :: A_vy, A_vz, A_By, A_Bz

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as n-D Cartesian
    {^IFONED call set_coordinate_system("Cartesian_1.75D")}
    {^IFTWOD call set_coordinate_system("Cartesian_2D")}
    {^IFTHREED call set_coordinate_system("Cartesian_3D")}

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Drive the wave using an internal boundary
    usr_internal_bc => Initialize_Wave

    usr_init_vector_potential=>initvecpot_usr

    ! Active the physics module
    call rmhd_activate() 
  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    call params_read(par_files)


    eg0 = p0/(rmhd_gamma - one) + half*(B0x*B0x + B0y*B0y + B0z*B0z)/(4.0*dpi)
    ca = dsqrt(rmhd_gamma*p0/rho0)
    a0 = dsqrt(p0/rho0)

    T0 = const_mp*fld_mu/const_kB*(p0/rho0)
    Er0 = (const_rad_a*(T0**4.d0))
    
    p0_B = (B0x*B0x + B0y*B0y + B0z*B0z)/(8.0*dpi)
    p0_by_Er0 = p0/Er0 
    p0_by_p0_B = p0/p0_B 

    va = dsqrt((B0x*B0x + B0y*B0y + B0z*B0z)/(4.0*dpi*rho0))
    vax = dsqrt((B0x*B0x)/(4.0*dpi*rho0))

    cslow = dsqrt(0.5d0*((va*va + ca*ca) - dsqrt((va*va + ca*ca)*(va*va + ca*ca) - 4.0d0*ca*ca*vax*vax)))
    cfast = dsqrt(0.5d0*((va*va + ca*ca) + dsqrt((va*va + ca*ca)*(va*va + ca*ca) - 4.0d0*ca*ca*vax*vax)))
    

    wvl = tau_wave/(rho0*fld_kappa0)
    omega = 2.d0*dpi*vax/wvl 
    wavenumber = 2.d0*dpi/wvl

    Bo = 4*rmhd_gamma*ca*eg0/(const_c*Er0) 
    r_Bo = Er0/(4*rmhd_gamma*eg0)

    !-------------------
    tau_c = const_c*fld_kappa0*rho0/omega
    tau_a = a0*fld_kappa0*rho0/omega

    if (mype .eq. 0) then
      print*, 'rho0', rho0
      print*, 'p0', p0
      print*, 'ca', ca
      print*, 'a0', a0
      print*, 'eg0', eg0
      print*, 'T0', T0
      print*, 'va', va
      print*, 'vax', vax
      print*, 'cslow', cslow
      print*, 'cfast', cfast
      print*, 'wvl', wvl
      print*, 'omega', omega
      print*, 'wavenumber', wavenumber
      print*, 'Bo', Bo
      print*, 'r_Bo', r_Bo
      print*, 'B0x', B0x
      print*, 'B0y', B0y
      print*, 'B0z', B0z
      print*, 'p0_B', p0_B
      print*, 'Er0', Er0
      print*, 'Er for radiative equilibrium', (const_rad_a*(T0**4.d0))
      print*, 'p0_by_Er0', p0_by_Er0
      print*, 'p0_by_p0_B', p0_by_p0_B
    endif
    

    ! Choose independent normalization units if using dimensionless variables.
    unit_length = wvl ! cm
    unit_velocity   = a0 ! cm/s
    unit_numberdensity = rho0/((1.d0+4.d0*He_abundance)*mp_cgs) ! cm^-3

    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kb)
    unit_time=unit_length/unit_velocity
    unit_magneticfield = sqrt(4.0*dpi * unit_pressure); 

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    rho0 = rho0/unit_density
    a0 = a0/unit_velocity
    p0 = p0/unit_pressure
    eg0 = eg0/unit_pressure
    T0 = T0/unit_temperature
    Er0 = Er0/unit_pressure

    wvl = wvl/unit_length
    omega = omega*unit_time
    wavenumber = wavenumber*unit_length

    B0x = B0x/unit_magneticfield
    B0y = B0y/unit_magneticfield
    B0z = B0z/unit_magneticfield

    p0_B = p0_B/unit_pressure
    p0_by_Er0 = p0/Er0 
    p0_by_p0_B = p0/p0_B

    ca = ca/unit_velocity
    va = va/unit_velocity
    vax = vax/unit_velocity
    cslow = cslow/unit_velocity
    cfast = cfast/unit_velocity

    if (mype .eq. 0) then
      print*, 'unit_length', unit_length
      print*, 'unit_density', unit_density
      print*, 'unit_pressure', unit_pressure
      print*, 'unit_temperature', unit_temperature
      print*, 'unit_radflux', unit_radflux
      print*, 'unit_opacity', unit_opacity
      print*, 'unit_time', unit_time
      print*, 'unit_velocity', unit_velocity
      print*, 'unit_magneticfield', unit_magneticfield
      print*, 'rho0', rho0
      print*, 'p0', p0
      print*, 'T0', T0
      print*, 'a0', a0
      print*, 'ca', ca
      print*, 'angular omega', omega/unit_time
      print*, 'wvl', wvl
      print*, 'opt tickness 1 wvl', tau_wave
      print*, 'wave number', wavenumber
      print*, 'amplitude', ampl
      print*, 'Bo', Bo
      print*, 'r_Bo', r_Bo
      print*, 'fld_mu' , fld_mu
      print*, 'B0x', B0x
      print*, 'B0y', B0y 
      print*, 'B0z', B0z 
      print*, 'p0_B', p0_B 
      print*, 'Er0', Er0
      print*, 'Er for radiative equilibrium', (const_rad_a*((T0*unit_temperature)**4.d0)/unit_pressure)
      print*, 'p0_by_Er0', p0_by_Er0
      print*, 'p0_by_p0_B', p0_by_p0_B
      print*, 'va', va
      print*, 'vax', vax
      print*, 'cslow', cslow
      print*, 'cfast', cfast
    endif

    !A_rho = ampl
    !A_v = omega/(wavenumber*rho0)*A_rho
    !A_p = (ca**2)*A_rho
    !A_Er = Er0/rho0*A_rho

    A_rho = 0.d0
    A_v = 0.d0
    A_p = 0.d0
    A_Er = 0.d0

    A_By = ampl
    A_Bz = 0.d0
    A_vy = -(B0x*wavenumber/(omega*rho0))*A_By
    A_vz = 0.d0
    A_e = A_p/(rmhd_gamma-one) + B0y*A_By + B0z*A_Bz

    if (mype .eq. 0) then
      print*, 'A_rho', A_rho
      print*, 'A_v', A_v
      print*, 'A_p', A_p
      print*, 'A_e', A_e
      print*, 'A_Er', A_Er
      print*, 'A_By', A_By
      print*, 'A_Bz', A_Bz
      print*, 'A_vy', A_vy
      print*, 'A_vz', A_vz
    endif

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wave_list/ rho0, p0, Er0, tau_wave, ampl, B0x, B0y, B0z

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wave_list, end=113)
113    close(unitpar)
    end do

  end subroutine params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: press(ixI^S), temp(ixI^S)
    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)

    ! Set initial values for w
    w(ixI^S, rho_) = rho0
    w(ixI^S, mom(:)) = 0.d0
    w(ixI^S, e_) = eg0

    press(ixI^S) = p0
    temp(ixI^S) = (const_mp*fld_mu/const_kb)*press(ixI^S)/w(ixI^S, rho_)&
    *(unit_pressure/unit_density)/unit_temperature

    w(ixI^S, r_e) = const_rad_a*(temp(ixI^S)*unit_temperature)**4.d0/unit_pressure

    !if (mype .eq. 0) then
    !  print*, 'ICs rho0', rho0
    !  print*, 'ICs eg0', eg0
    !  print*, 'ICs p0', p0
    !  print*, 'ICs T', (const_mp*fld_mu/const_kb)*p0/rho0 *(unit_pressure/unit_density)/unit_temperature
    !  print*, 'ICs r_e', const_rad_a*((const_mp*fld_mu/const_kb)*p0/rho0*(unit_pressure/unit_density))**4.d0/unit_pressure
    !endif

    if(stagger_grid) then
      call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
      call rmhd_face_to_center(ixO^L,block)
    else
      w(ixI^S,mag(:))= 0.d0
      w(ixI^S,mag(1))= B0x
      w(ixI^S,mag(2))= B0y
      w(ixI^S,mag(3))= B0z
    end if
  end subroutine initial_conditions


  subroutine Initialize_Wave(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_fld
    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision :: press(ixI^S), temp(ixI^S)

    where (x(ixI^S,1) .lt. one)
      w(ixI^S, rho_) = rho0 + A_rho*dsin(wavenumber*x(ixI^S,1)-omega*global_time)
      {^NOONED w(ixI^S, mom(2)) = zero}
      {^IFTHREED w(ixI^S, mom(3)) = zero}
      w(ixI^S, mom(1)) = w(ixI^S, rho_)*A_v*dsin(wavenumber*x(ixI^S,1)-omega*global_time)
      w(ixI^S, e_) = eg0 + A_e*dsin(wavenumber*x(ixI^S,1)-omega*global_time)

      press(ixI^S) = p0 + A_p*dsin(wavenumber*x(ixI^S,1)-omega*global_time)
      temp(ixI^S) = (const_mp*fld_mu/const_kb)*press(ixI^S)/w(ixI^S, rho_)&
      *(unit_pressure/unit_density)/unit_temperature

      w(ixI^S, r_e) = const_rad_a*(temp(ixI^S)*unit_temperature)**4.d0/unit_pressure

      w(ixI^S, mag(1)) = B0x
      w(ixI^S, mag(2)) = B0y + A_By*dsin(wavenumber*x(ixI^S,1)-omega*global_time)
      w(ixI^S, mag(3)) = B0z + A_Bz*dsin(wavenumber*x(ixI^S,1)-omega*global_time)

      w(ixI^S, mom(2)) = w(ixI^S, rho_)*A_vy*dsin(wavenumber*x(ixI^S,1)-omega*global_time)
      w(ixI^S, mom(3)) = w(ixI^S, rho_)*A_vz*dsin(wavenumber*x(ixI^S,1)-omega*global_time)

      w(ixI^S, e_) = press(ixI^S)/(rmhd_gamma - one) + half*(B0x*B0x + &
       w(ixI^S, mag(2))*w(ixI^S, mag(2)) + w(ixI^S, mag(3))*w(ixI^S, mag(3)))
      
    endwhere

  end subroutine Initialize_Wave

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the corners
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    A(ixC^S) = zero

  end subroutine initvecpot_usr

end module mod_usr
