!> Convergence study for scalar advection (phys='rho').
!>
!> A smooth periodic wave rho = rho0 + amp*sin(2*pi*k.(x - v t)) on the unit
!> cube. Both the initial condition and the reference solution are *exact
!> cell averages* (the sine times a sinc factor per direction), so the error
!> measured is the scheme's alone. With typefilelog='special' the log hook
!> conv_print_log writes, at every log save,
!>
!>    time   L1   L2   Linf   mean(rho)
!>
!> of rho minus the exact cell average, to <base_filename>.log.
!>
!> &usr_list setup:
!>   'diag' - k = (1,1,1): the field varies along every face (fully 3D test);
!>            with rho_v = (1,1,1) one period is t = 1/3
!>   'x'    - k = (1,0,0): effectively 1D; with rho_v = (1,0,0) one period is
!>            t = 1
!> The runs are driven by run_convergence.sh; convergence.py makes the table.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  character(len=20) :: setup = 'diag'
  double precision  :: rho0 = 1.0d0, amp = 0.5d0
  double precision  :: kvec(3) = [1.0d0, 1.0d0, 1.0d0]
  double precision, parameter :: pi = 3.14159265358979323846d0

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr
    usr_print_log     => conv_print_log

    call params_read_user(par_files)
    call phys_activate()

  end subroutine usr_init

  subroutine params_read_user(files)
    use mod_global_parameters, only: unitpar
    use mod_comm_lib, only: mpistop
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ setup, rho0, amp

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

    select case (trim(setup))
    case ('diag')
       kvec = [1.0d0, 1.0d0, 1.0d0]
    case ('x')
       kvec = [1.0d0, 0.0d0, 0.0d0]
    case default
       call mpistop("usr_list: setup must be 'diag' or 'x'")
    end select

  end subroutine params_read_user

  !> Exact average of the wave over the cell centred at xc with widths dxc
  pure double precision function exact_avg(xc, dxc, t) result(r)
    double precision, intent(in) :: xc(3), dxc(3), t
    integer          :: d
    double precision :: phase, fac, a

    phase = 0.0d0
    fac   = 1.0d0
    do d = 1, 3
       phase = phase + kvec(d)*(xc(d) - rho_v(d)*t)
       if (kvec(d) /= 0.0d0) then
          a   = pi*kvec(d)*dxc(d)
          fac = fac*sin(a)/a
       end if
    end do
    r = rho0 + amp*fac*sin(2.0d0*pi*phase)

  end function exact_avg

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    integer          :: ix1, ix2, ix3
    double precision :: dxc(3)

    ! uniform Cartesian block: the cell widths follow from the centres
    dxc(1) = x(ixGmin1+1,ixGmin2,ixGmin3,1) - x(ixGmin1,ixGmin2,ixGmin3,1)
    dxc(2) = x(ixGmin1,ixGmin2+1,ixGmin3,2) - x(ixGmin1,ixGmin2,ixGmin3,2)
    dxc(3) = x(ixGmin1,ixGmin2,ixGmin3+1,3) - x(ixGmin1,ixGmin2,ixGmin3,3)

    do ix3 = ixGmin3, ixGmax3
       do ix2 = ixGmin2, ixGmax2
          do ix1 = ixGmin1, ixGmax1
             w(ix1,ix2,ix3,rho_) = exact_avg(x(ix1,ix2,ix3,1:3), dxc, 0.0d0)
          end do
       end do
    end do

  end subroutine initonegrid_usr

  !> Error norms against the exact cell averages, and the mean of rho
  subroutine conv_print_log()
    use mod_global_parameters
    integer          :: iigrid, igrid, ix1, ix2, ix3, u
    double precision :: dxc(3), dvol, e, lsum(4), gsum(4), lmax, gmax
    logical, save    :: first = .true.

    lsum = 0.0d0
    lmax = 0.0d0
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       dxc(1) = rnode(rpdx1_,igrid)
       dxc(2) = rnode(rpdx2_,igrid)
       dxc(3) = rnode(rpdx3_,igrid)
       dvol   = dxc(1)*dxc(2)*dxc(3)
       do ix3 = ixMlo3, ixMhi3
          do ix2 = ixMlo2, ixMhi2
             do ix1 = ixMlo1, ixMhi1
                e = ps(igrid)%w(ix1,ix2,ix3,rho_) - &
                    exact_avg(ps(igrid)%x(ix1,ix2,ix3,1:3), dxc, global_time)
                lsum(1) = lsum(1) + abs(e)*dvol
                lsum(2) = lsum(2) + e*e*dvol
                lsum(3) = lsum(3) + ps(igrid)%w(ix1,ix2,ix3,rho_)*dvol
                lsum(4) = lsum(4) + dvol
                lmax    = max(lmax, abs(e))
             end do
          end do
       end do
    end do

    call MPI_ALLREDUCE(lsum, gsum, 4, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(lmax, gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, icomm, ierrmpi)

    if (mype == 0) then
       open(newunit=u, file=trim(base_filename)//'.log', status='unknown', &
            position='append', action='write')
       if (first) write(u,'(a)') '# time L1 L2 Linf mean(rho)'
       write(u,'(5es24.15)') global_time, gsum(1)/gsum(4), &
            sqrt(gsum(2)/gsum(4)), gmax, gsum(3)/gsum(4)
       close(u)
    end if
    first = .false.

  end subroutine conv_print_log

end module mod_usr
