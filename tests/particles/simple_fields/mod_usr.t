module mod_usr
  use mod_mhd
  use mod_particles

  implicit none

  double precision :: charge = 1.0d0
  double precision :: mass = 1.0d0

  ! Initial position (in m)
  double precision :: x0(3) = [0.0d0, 0.0d0, 0.0d0]

  ! Initial velocity (in m/s)
  double precision :: v0(3) = [0.0d0, 0.0d0, 0.0d0]

  ! Maxwellian velocity (per component vx, vy, vz)
  double precision :: maxwellian_velocity = 0.0d0

  double precision, parameter :: not_used_value = -1.0d20
  double precision :: force_E0(3) = [not_used_value, 0.0d0, 0.0d0]
  double precision :: force_B0(3) = [not_used_value, 0.0d0, 0.0d0]

  ! Use an analytic field instead of an interpolated one
  logical :: use_analytic_field = .false.

contains

  subroutine usr_init()
    use mod_initialize

    unit_length        = 1.d0
    unit_numberdensity = 1.d0
    unit_velocity      = 1.0d0

    usr_init_one_grid => initonegrid_usr
    usr_create_particles => generate_particles
    usr_particle_fields => set_custom_field

    call set_coordinate_system("Cartesian_3D")
    call mhd_activate()
    call params_read(par_files)

    call initialize_amrvac()    ! So that we have settings available

    if (use_analytic_field) then
      if (physics_type_particles /= 'Lorentz' .and. physics_type_particles /= 'Vay' .and. physics_type_particles /= 'HC') &
           call mpistop('Analytic fields only supported with Boris, HC or Vay schemes')
      usr_particle_analytic => get_analytic_field
    end if

  end subroutine usr_init

  !> Read parameters from a file
  subroutine params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /my_list/ charge, mass, x0, v0, use_analytic_field, force_E0, &
         force_B0, maxwellian_velocity

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
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
      write(*,*) 'Test particles in 3D simple B-field'
      write(*,*) 'rho - p - b:',rho0,p0,b0
      first=.false.
      print *, "Running test case", iprob
    endif

  end subroutine initonegrid_usr

  subroutine generate_particles(n_particles, x, v, q, m, follow)
    use mod_particles
    integer, intent(in)           :: n_particles
    double precision, intent(out) :: x(3, n_particles)
    double precision, intent(out) :: v(3, n_particles)
    double precision, intent(out) :: q(n_particles)
    double precision, intent(out) :: m(n_particles)
    logical, intent(out)          :: follow(n_particles)
    integer                       :: n

    do n = 1, n_particles
      call get_particle(x(:, n), v(:, n), q(n), m(n), n, n_particles, iprob)
    end do

    follow(:) = .true.

    ! Scale to CGS units
    x = x * 1d2             ! m to cm
    v = v * 1d2             ! to cm/s
    q = q * const_c * 0.1d0 ! A unit charge converted to CGS units
    m = m * 1d3             ! kg to gram
    !print*,'mass in cgs:',m,'charge in cgs:',q,'q/m ratio',q/m
  end subroutine generate_particles

  ! Return field at location x (in SI units: Tesla, V/m)
  subroutine get_field(x, E, B)
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
      B = [0.0d0, 0.0d0, 1.0d0]
    case (3)
      ! Force-free (this only works for one particle)
      E = [-v0(2), 0.0d0, 0.0d0]
      B = [0.0d0, 0.0d0, 1.0d0]
    case (4)
      ! ExB
      E = [1.0d0, 0.0d0, 0.0d0]
      B = [0.0d0, 0.0d0, 1.0d0]
    case (5)
      ! Gradient in B
      E = [0.0d0, 0.0d0, 0.0d0]
      B = [0.0d0, 0.0d0, (1.0d0 + 1.0d-2 * x(1))]
    case (6)
      ! Magnetic mirror (requires longer time a.t.m.)
      E = [0.0d0, 0.0d0, 0.0d0]
      ! x is in cm
      !B = [-x(1) * x(3), -x(2) * x(3), 1.0d4 + x(3)**2] * 1.0d-4
      B = 1.0d6 * [-x(1) * x(3), -x(2) * x(3), 1.0d18 + x(3)**2] * 1.0d-18
    case (7)
      ! Magnetic dipole (run up to t = 100)
      E = [0.0d0, 0.0d0, 0.0d0]
      ! x is in cm, this corresponds to B = 10 T at 1 m
      B = 10 * 1d6 * [3d0 * x(1) * x(3), &
           3d0 * x(2) * x(3), &
           2d0 * x(3)**2 - x(1)**2 - x(2)**2] / &
           (x(1)**2 + x(2)**2 + x(3)**2)**2.5d0
    case (8)
      ! X-null point
      E = 0.0d0
      B = [x(2), x(1), 0.0d0] * 1d-2
    case (9)
      ! electromagnetic two-body problem. x is in cm 
      ! Q=-1 (attracting central electron surrounded by positron)
      E = -1.0d0 * 1.0d14 * [x(1), x(2), 0.0d0] / &
           (x(1)**2 + x(2)**2 + x(3)**2)**(3.0d0/2.0d0)
      B = [0.0d0, 0.0d0, 0.0d0]
    case default
      call mpistop("Unknown value for iprob")
    end select

    ! The user can override these fields
    if (force_E0(1) > not_used_value) E = force_E0
    if (force_B0(1) > not_used_value) B = force_B0

  end subroutine get_field

  ! Set particle properties (in SI units)
  subroutine get_particle(x, v, q, m, ipart, n_particles, iprob)
    double precision, intent(out) :: x(3), v(3), q, m
    integer, intent(in)           :: ipart, iprob, n_particles
    double precision              :: tmp_vec(4), phi

    q = charge
    m = mass

    x = x0

    select case (iprob)
    case (4)
       ! Add Maxwellian velocity. Random numbers come in pairs of two
       tmp_vec(1:2) = rng%two_normals()
       tmp_vec(3:4) = rng%two_normals()
       v = v0 + tmp_vec(1:3) * maxwellian_velocity
       print *, ipart, v
    case (5)
       v = (v0 * ipart) / n_particles
    case (7)
       v = v0
       q = (charge * ipart) / n_particles
       if (physics_type_particles /= 'gca') then
          ! Assume B = 10 T, and v_x = 0 initially
          x(1) = x(1) + abs(v(2)) * m / (q * 10.0d0)
       end if
    case (8)
       ! Distribute over circle, velocity inwards. Avoid pi/4.
       phi = ((ipart+0.125d0) * 2 * acos(-1.0d0)) / n_particles
       x = norm2(x0) * [cos(phi), sin(phi), 0.0d0]
       v = -x * norm2(v0)/norm2(x0)
    case default
       v = v0
    end select
  end subroutine get_particle

  subroutine set_custom_field(w, x, E_field, B_field)
    double precision, intent(in)  :: w(ixG^T,nw)
    double precision, intent(in)  :: x(ixG^T,ndim)
    double precision, intent(out) :: E_field(ixG^T, ndir)
    double precision, intent(out) :: B_field(ixG^T, ndir)

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
  !print*,'E in cgs:',E(1) * 1.0d6/const_c, 'B in cgs',B(3) * 1.0d4
  !print*,'E/B ratio < 1',(E(1)*1.0d6/const_c)/(B(3)*1.0d4)
  end subroutine set_custom_field

  subroutine get_analytic_field(ix, x, tloc, vec)
    integer, intent(in)           :: ix(ndir) !< Indices in gridvars
    double precision, intent(in)  :: x(ndir)
    double precision, intent(in)  :: tloc
    double precision, intent(out) :: vec(ndir)
    double precision              :: E(3), B(3)

    call get_field(x, E, B)

    if (ix(1) == bp(1)) then
      vec = B * 1.0d4
    else if (ix(1) == ep(1)) then
      vec = E * 1.0d6/const_c
    else
      call mpistop("get_analytic_field: unknown variable index")
    end if
  end subroutine get_analytic_field

end module mod_usr
