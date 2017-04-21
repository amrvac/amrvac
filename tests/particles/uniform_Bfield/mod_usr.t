module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    unit_length=1.d0
    unit_numberdensity=1.d0
    unit_velocity=1.0d0

    usr_init_one_grid => initonegrid_usr
    usr_create_particles => generate_particles

    call set_coordinate_system("Cartesian_3D")

    call mhd_activate()

  end subroutine usr_init

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

    if(first .and. mype==0 )then
      write(*,*) 'Particles in 3D homogeneous B-field'
      write(*,*) 'rho - p - b:',rho0,p0,b0
      first=.false.
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
      ! call get_uniform_pos(x(:, n))
      ! call get_maxwellian_velocity(v(:, n), 1.0d0/const_c)
      x(:, n) = 0.5d0
      v(:, n) = [const_e/const_me, 0.0d0, 0.0d0] / const_c
      follow(n) = .true.
      q(n)      = 1.0d0
      m(n)      = 1.0d0
    end do
  end subroutine generate_particles

end module mod_usr
