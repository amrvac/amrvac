!> Scalar advection of a ball (or a Gaussian) across a periodic cube.
!>
!> Port of MPI-AMRVAC's tests/rho/auto_3d (ball_3d.par, iprob=6). All runs in
!> this directory share one build and pick their initial condition through the
!> `setup` entry of the &usr_list namelist:
!>   test.par         - the ball, Lohner AMR (refine_max_level=2)
!>   test-fixgrid.par - the same with itfixgrid=0, so the ball crosses frozen
!>                      level jumps: on this closed (periodic) domain mean(rho)
!>                      is conserved, so it is the refluxing test
!>   hll.par          - uniform grid, HLL (exact upwinding), koren, twostep
!>   agile.par        - build reference and a visualisable run to t=1, when
!>                      the ball is back at the centre
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> 'ball' (discontinuous, upstream iprob=6) or 'gauss'
  character(len=20) :: setup = 'ball'
  double precision  :: rho_in = 2.0d0, rho_out = 0.5d0, radius = 0.2d0

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr

    call params_read_user(par_files)
    call phys_activate()

  end subroutine usr_init

  subroutine params_read_user(files)
    use mod_global_parameters, only: unitpar
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ setup, rho_in, rho_out, radius

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read_user

  !> A ball (or Gaussian) centred in the domain on a uniform background
  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    use mod_comm_lib, only: mpistop
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision :: r2(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3)

    r2 = (x(:,:,:,1) - 0.5d0*(xprobmin1+xprobmax1))**2 &
       + (x(:,:,:,2) - 0.5d0*(xprobmin2+xprobmax2))**2 &
       + (x(:,:,:,3) - 0.5d0*(xprobmin3+xprobmax3))**2

    select case (trim(setup))
    case ('ball')
       where (r2 < radius**2)
          w(:,:,:,rho_) = rho_in
       elsewhere
          w(:,:,:,rho_) = rho_out
       end where
    case ('gauss')
       w(:,:,:,rho_) = rho_out + (rho_in - rho_out)*exp(-r2/radius**2)
    case default
       call mpistop("usr_list: setup must be 'ball' or 'gauss'")
    end select

  end subroutine initonegrid_usr

end module mod_usr
