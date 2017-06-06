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
    usr_particle_fields => custom_field

    call set_coordinate_system("Cartesian_3D")
    call mhd_activate()
    call params_read(par_files)

    call initialize_amrvac()    ! So that we have settings available

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
  subroutine get_field(x, E, B)
    use mod_global_parameters
    double precision, intent(in)  :: x(3)
    double precision, intent(out) :: E(3), B(3)

    select case (iprob)
    case (1)
      ! Linear acceleration
      E = [0.0d0, 0.0d0, 1.0d0]
      B = [0.0d0, 0.0d0, 1.0d0]
    case (2)
      ! Pure gyration
      E = [0.0d0, 0.0d0, 0.0d0]
      B = [0.0d0, 0.0d0, 8 * dpi]
    case (3)
      ! Force-free
      E = [-1.0d0, 0.0d0, 0.0d0]
      B = [0.0d0, 0.0d0, 1.0d0]
    case (4)
      ! ExB
      E = [16 * dpi, 0.0d0, 0.0d0]
      B = [0.0d0, 0.0d0, 8 * dpi]
    case (5)
      ! Gradient in B
      E = [0.0d0, 0.0d0, 0.0d0]
      B = [0.0d0, 0.0d0, 4 * dpi * (2.0d0 + 1.0d-2 * x(1))]
    case (6)
      ! Magnetic mirror (requires longer time a.t.m.)
      E = [0.0d0, 0.0d0, 0.0d0]
      B = [-x(1) * x(3), -x(2) * x(3), 1d2 + x(3)**2] * 1d-2
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
      ! Pure gyration centered around (0.5, 0.5)
      x = [0.0d0, 0.0d0, 0.0d0]
      v = [-4 * dpi, 0.0d0, 0.0d0]
    case (3)
      ! Force-free
      x = [0.0d0, 0.0d0, 0.0d0]
      v = [0.0d0, 1.0d0, 0.0d0]
    case (4)
      ! ExB
      x = [0.0d0, 0.0d0, 0.0d0]
      v = [0.0d0, 0.0d0, 0.0d0]
    case (5)
      ! Gradient in B
      x = [0.0d0, 0.0d0, 0.0d0]
      v = [-dpi, 0.0d0, 0.0d0]
    case (6)
      ! Magnetic mirror
      x = [0.0d0, 0.0d0, 0.0d0]
      v = [0.0d0, 1.0d0, 1.0d0]
    case default
      call mpistop("Unknown value for iprob")
    end select
  end subroutine get_particle

  subroutine custom_field(x, E_field, B_field)
    use mod_global_parameters
    double precision, intent(in)  :: x(ixG^T,1:ndim)
    double precision, intent(out) :: E_field(ixG^T, 1:ndir)
    double precision, intent(out) :: B_field(ixG^T, 1:ndir)

    integer          :: i^D
    double precision :: E(3), B(3), xtmp(ndim)


    {do i^D = ixGlo^D, ixGhi^D\}
    xtmp = x(i^D, :)
    call get_field(xtmp, E, B)

    ! Convert to CGS units, 1 T -> 1e4 Gauss
    B_field(i^D, :) = B * 1.0d4

    ! Convert to CGS units
    E_field(i^D, :) = E * 1.0d6/const_c
    {end do\}

  end subroutine custom_field

end module mod_usr
