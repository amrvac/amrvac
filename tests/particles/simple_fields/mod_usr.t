module mod_usr
  use mod_mhd
  use mod_particles

  implicit none

  double precision :: charge = 1.0d0
  double precision :: mass = 1.0d0

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_initialize

    unit_length        = 1.d0
    unit_numberdensity = 1.d0
    unit_velocity      = 1.0d0

    usr_init_one_grid => initonegrid_usr
    usr_create_particles => generate_particles

    call set_coordinate_system("Cartesian_3D")
    call mhd_activate()
    call params_read(par_files)

    call initialize_amrvac()    ! So that we have settings available

    if (physics_type_particles == 'Lorentz') then
      particles_fill_gridvars => custom_field_Lorentz
    else if (physics_type_particles == 'gca') then
      particles_fill_gridvars => custom_field_gca
    else
      call mpistop('This type of particle mover is not supported here')
    end if

  end subroutine usr_init

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /my_list/ charge, mass

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision:: rho0,p0,b0
    logical, save :: first=.true.

    rho0 = 1.0d0
    p0 = 1.0d0
    b0 = 1.0d0

    w(ixO^S,rho_)= rho0
    w(ixO^S,mom(1))= 0.0d0
    w(ixO^S,mom(2))= 0.0d0
    w(ixO^S,mom(3))= 0.0d0
    w(ixO^S,p_)=p0
    w(ixO^S,mag(1))= 0.0d0
    w(ixO^S,mag(2))= 0.0d0
    w(ixO^S,mag(3))= b0

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

    if (first .and. mype==0 )then
      write(*,*) 'Particles in 3D homogeneous B-field'
      write(*,*) 'rho - p - b:',rho0,p0,b0
      first=.false.
      print *, "Running test case", iprob
    endif

  end subroutine initonegrid_usr

  subroutine generate_particles(n_particles, x, v, q, m, follow)
    use mod_global_parameters
    use mod_particles
    integer, intent(in)           :: n_particles
    double precision, intent(out) :: x(3, n_particles)
    double precision, intent(out) :: v(3, n_particles)
    double precision, intent(out) :: q(n_particles)
    double precision, intent(out) :: m(n_particles)
    logical, intent(out)          :: follow(n_particles)
    integer                       :: n

    do n = 1, n_particles
      call get_particle(x(:, n), v(:, n), q(n), m(n), n, iprob)
    end do

    follow(:) = .true.

    ! Scale to CGS units
    x = x * 1d2 ! m to cm

    if (physics_type_particles == 'Lorentz') then
      v = v * 1d2 / const_c     ! to cm/s, then normalize by speed of light
    else if (physics_type_particles == 'gca') then
      v = v * 1d2               ! to cm/s
    else
      call mpistop('This type of particle mover is not supported here')
    end if

    q      = q * const_c * 0.1d0 ! A unit charge converted to CGS units
    m      = m * 1d3             ! kg to gram
  end subroutine generate_particles

  ! Return field at location x (in SI units: Tesla, V/m)
  subroutine get_field(x, E, B, gradB)
    use mod_global_parameters
    double precision, intent(in)  :: x(3)
    double precision, intent(out) :: E(3), B(3), gradB(3)

    gradB = 0.0d0

    select case (iprob)
    case (1)
      ! Linear acceleration
      E = [0.0d0, 0.0d0, 1.0d0]
      B = [0.0d0, 0.0d0, 1.0d0]
    case (2)
      ! Pure gyration
      E = [0.0d0, 0.0d0, 0.0d0]
      B = [0.0d0, 0.0d0, 1.0d0]
    case (3)
      ! Force-free
      E = [-1.0d0, 0.0d0, 0.0d0]
      B = [0.0d0, 0.0d0, 1.0d0]
    case (4)
      ! ExB
      E = [1.0d0, 0.0d0, 0.0d0]
      B = [0.0d0, 0.0d0, 1.0d0]
    case (5)
      ! Gradient in B
      E = [0.0d0, 0.0d0, 0.0d0]
      B = [0.0d0, 0.0d0, 1.0d0 + 1.0d-2 * x(1)]
      gradB(3) = 1.0d-2
    case default
      call mpistop("Unknown value for iprob")
    end select
  end subroutine get_field

  ! Set particle properties (SI units)
  subroutine get_particle(x, v, q, m, ipart, iprob)
    double precision, intent(out) :: x(3), v(3), q, m
    integer, intent(in)           :: ipart, iprob

    q = charge
    m = mass

    select case (iprob)
    case (1)
      ! Linear acceleration
      x = [0.0d0, 0.0d0, 0.0d0]
      v = [0.0d0, 0.0d0, 0.0d0]
    case (2)
      ! Pure gyration
      q = 2 * dpi
      x = [0.0d0, 0.0d0, 0.0d0]
      v = [1.0d0, 0.0d0, 0.0d0]
    case (3)
      ! Force-free
      x = [0.0d0, 0.0d0, 0.0d0]
      v = [0.0d0, 1.0d0, 0.0d0]
    case (4)
      ! ExB
      x = [0.0d0, 0.0d0, 0.0d0]
      v = [0.0d0, 0.0d0, 0.0d0]
    case (5)
      ! ExB
      x = [0.0d0, 0.0d0, 0.0d0]
      v = [0.0d0, 1.0d0, 0.0d0]
    case default
      call mpistop("Unknown value for iprob")
    end select
  end subroutine get_particle

  subroutine custom_field_Lorentz()
    use mod_particle_Lorentz
    use mod_global_parameters

    integer          :: igrid, iigrid, idir
    integer          :: idirmin, i^D
    double precision :: x(3), E(3), B(3), gradB(3)

    do iigrid=1,igridstail; igrid=igrids(iigrid);
       gridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0

       {do i^D = ixGlo^D, ixGhi^D\}
         x = pw(igrid)%x(i^D, :)
         call get_field(x, E, B, gradB)

         ! Convert to CGS units, 1 T -> 1e4 Gauss
         gridvars(igrid)%w(i^D,bp(:)) = B * 1.0d4

         ! Convert to CGS units
         gridvars(igrid)%w(i^D,ep(:)) = E * 1.0d6/const_c
       {end do\}

       ! The code interpolates between two states in time (even though we don't
       ! need it here)
       if (time_advance) then
         gridvars(igrid)%wold(ixG^T,bp(:)) = &
            gridvars(igrid)%w(ixG^T,bp(:))
         gridvars(igrid)%wold(ixG^T,ep(:)) = &
            gridvars(igrid)%w(ixG^T,ep(:))
       end if
    end do

  end subroutine custom_field_Lorentz

  subroutine custom_field_gca()
    use mod_particle_gca
    use mod_global_parameters

    integer          :: igrid, iigrid, idir
    integer          :: idirmin, i^D
    double precision :: x(3), E(3), B(3), gradB(3)

    do iigrid=1,igridstail; igrid=igrids(iigrid);
       gridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0

       {do i^D = ixGlo^D, ixGhi^D\}
         x = pw(igrid)%x(i^D, :)
         call get_field(x, E, B, gradB)

         ! Convert to CGS units, 1 T -> 1e4 Gauss
         gridvars(igrid)%w(i^D,bp(:)) = B * 1.0d4

         ! Convert to CGS units
         gridvars(igrid)%w(i^D,ep(:)) = E * 1.0d6/const_c

         gridvars(igrid)%w(i^D,grad_kappa_B(:)) = gradB(:) * 1.0d4
       {end do\}

       ! The code interpolates between two states in time (even though we don't
       ! need it here)
       if (time_advance) then
         gridvars(igrid)%wold(ixG^T,bp(:)) = &
            gridvars(igrid)%w(ixG^T,bp(:))
         gridvars(igrid)%wold(ixG^T,ep(:)) = &
            gridvars(igrid)%w(ixG^T,ep(:))
       end if
    end do

  end subroutine custom_field_gca

end module mod_usr
