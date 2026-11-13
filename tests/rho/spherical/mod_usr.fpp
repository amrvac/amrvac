#:mute
#:include "mod_gpu_directives.fpp"
#:endmute

!> Scalar advection on a spherical shell away from the poles:
!> r in [1,2], theta in [pi/4, 3pi/4], phi in [0, pi/2].
!>
!> Two runs share this build, selected by &usr_list:
!>   uflow.par - amp = 0: a uniform rho in the uniform Cartesian flow rho_v,
!>               an exact steady solution. Written in spherical components the
!>               velocity is a non-trivial function of position, so this
!>               exercises every curvilinear flux. The face-normal velocity is
!>               a point value (see mod_rho_templates.fpp), so rho stays
!>               uniform to O(h^2), not to round-off.
!>   blob.par  - amp > 0: a Gaussian blob carried by rho_v.
!> The physical boundaries hold the exact solution (special), so only the
!> interior discretisation is under test.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> blob amplitude (0: uniform), width and Cartesian centre at t = 0. Plain
  !> variables rather than parameters: specialbound_usr runs on the device.
  double precision :: amp = 0.0d0, sigma = 0.2d0, x0(3) = 0.0d0
  ${GPU_DECLARE_COPYIN('amp, sigma, x0')}$

contains

  subroutine usr_init()

    ! the coordinate system comes from `geometry` in &meshlist
    call params_read_user(par_files)
    ${GPU_UPDATE_DEVICE('amp, sigma, x0')}$

    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr

    call phys_activate()

  end subroutine usr_init

  subroutine params_read_user(files)
    use mod_global_parameters, only: unitpar
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ amp, sigma, x0

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read_user

  !> Exact solution at the spherical position x and time t: rho = 1 plus the
  !> blob, translated with rho_v. Valid at any coordinates, the mirrored ones
  !> of the ghost cells beyond an axis included.
  pure double precision function rho_exact(x, t) result(r)
#ifdef _OPENMP
    ${GPU_ROUTINE_VECTOR()}$
#else
    ${GPU_ROUTINE_SEQ()}$
#endif
    double precision, intent(in) :: x(1:ndim), t
    double precision             :: xc(3)

    xc(1) = x(1)*sin(x(2))*cos(x(3))
    xc(2) = x(1)*sin(x(2))*sin(x(3))
    xc(3) = x(1)*cos(x(2))
    r = 1.0d0 + amp*exp(-sum((xc - x0 - rho_v*t)**2)/sigma**2)

  end function rho_exact

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             w(ix1,ix2,ix3,rho_) = rho_exact(x_loc, 0.0d0)
          end do
       end do
    end do

  end subroutine initonegrid_usr

  !> keep the exact solution in the physical-boundary ghost cells
  subroutine specialbound_usr(qt, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB, w, x)
    use mod_global_parameters
    ${GPU_ROUTINE_VECTOR()}$
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3
    integer, intent(in)             :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    integer, intent(in)             :: iB
    double precision, intent(in)    :: qt
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    ! collapse(3) with a call in the body is safe here: the only private
    ! array, x_loc, has a compile-time size (see "Bug-hunting notes" in
    ! CLAUDE.md and issue #154)
#ifndef _OPENMP
    ${GPU_LOOP_VECTOR("collapse(3) private(x_loc)")}$
#endif
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             w(ix1,ix2,ix3,iw_rho) = rho_exact(x_loc, qt)
          end do
       end do
    end do

  end subroutine specialbound_usr

end module mod_usr
