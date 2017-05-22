module mod_usr
  use mod_mhd
  use mod_particles

  implicit none

  double precision :: Efield_SI(3) = [0.0d0, 0.0d0, 1.0d0]
  double precision :: Bfield_SI(3) = [0.0d0, 0.0d0, 0.0d0]
  double precision :: init_gammav_SI(3) = [1d2, 0.0d0, 0.0d0]

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    unit_length        = 1.d0
    unit_numberdensity = 1.d0
    unit_velocity      = 1.0d0

    usr_init_one_grid => initonegrid_usr
    usr_create_particles => generate_particles
    phys_fill_gridvars => custom_field

    call set_coordinate_system("Cartesian_3D")
    call mhd_activate()
    call params_read(par_files)

  end subroutine usr_init

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /my_list/ Efield_SI, Bfield_SI, init_gammav_SI

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

      ! Velocities are assumed to be normalized by the speed of light for the
      ! Lorentz mover
      v(:, n) = 1d2 * init_gammav_SI / const_c
      follow(n) = .true.

      ! A unit charge converted to CGS units
      q(n)      = const_c * 0.1d0
      m(n)      = 1d3
    end do
  end subroutine generate_particles

  subroutine custom_field()
    use mod_global_parameters

    integer                                   :: igrid, iigrid, idir
    double precision, dimension(ixG^T,1:ndir) :: beta
    double precision, dimension(ixG^T,1:nw)   :: w,wold
    integer                                   :: idirmin

    do iigrid=1,igridstail; igrid=igrids(iigrid);
       gridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0

       ! fill with magnetic field (converted to CGS units, 1 T -> 1e4 Gauss)
       gridvars(igrid)%w(ixG^T,bp(1)) = Bfield_SI(1) * 1.0d4
       gridvars(igrid)%w(ixG^T,bp(2)) = Bfield_SI(2) * 1.0d4
       gridvars(igrid)%w(ixG^T,bp(3)) = Bfield_SI(3) * 1.0d4

       gridvars(igrid)%w(ixG^T,ep(1)) = Efield_SI(1) * 1.0d6/const_c
       gridvars(igrid)%w(ixG^T,ep(2)) = Efield_SI(2) * 1.0d6/const_c
       gridvars(igrid)%w(ixG^T,ep(3)) = Efield_SI(3) * 1.0d6/const_c

       ! The code interpolates between two states in time (even though we don't
       ! need it here)
       if (time_advance) then
         gridvars(igrid)%wold(ixG^T,bp(:)) = &
            gridvars(igrid)%w(ixG^T,bp(:))
         gridvars(igrid)%wold(ixG^T,ep(:)) = &
            gridvars(igrid)%w(ixG^T,ep(:))
       end if
    end do

  end subroutine custom_field

end module mod_usr
