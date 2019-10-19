!> HD 2D riemann problem: can be used to check Lax&Liu, SIAM, 1998
module mod_usr
  use mod_hd

  implicit none

  ! input 4 states in primitive variables
  double precision :: rho1, p1, u1, v1
  double precision :: rho2, p2, u2, v2
  double precision :: rho3, p3, u3, v3
  double precision :: rho4, p4, u4, v4

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ rho1,p1,u1,v1,rho2,p2,u2,v2,rho3,p3,u3,v3,rho4,p4,u4,v4

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine usr_init()
    call usr_params_read(par_files)

    usr_init_one_grid => rm2d_init_one_grid

    call hd_activate()
  end subroutine usr_init

  ! Initialize one grid
  subroutine rm2d_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    {^IFONED call mpistop("This is a multi-D HD problem") }

    ! the left top  quadrant, state 2
    where (x(ix^S,1)<0.5d0 .and. x(ix^S,2)>=0.5d0)
       w(ix^S,e_)     = p2
       w(ix^S,rho_)   = rho2
       w(ix^S,mom(1)) = u2
       w(ix^S,mom(2)) = v2
    end where

    ! the right top  quadrant, state 1
    where (x(ix^S,1)>=0.5d0 .and. x(ix^S,2)>=0.5d0)
       w(ix^S,e_)     = p1
       w(ix^S,rho_)   = rho1
       w(ix^S,mom(1)) = u1
       w(ix^S,mom(2)) = v1
    end where

    ! the left bottom quadrant, state 3
    where (x(ix^S,1)<0.5d0 .and. x(ix^S,2)<0.5d0)
       w(ix^S,e_)     = p3
       w(ix^S,rho_)   = rho3
       w(ix^S,mom(1)) = u3
       w(ix^S,mom(2)) = v3
    end where

    ! the right bottom quadrant, state 4
    where (x(ix^S,1)>=0.5d0 .and. x(ix^S,2)<0.5d0)
       w(ix^S,e_)     = p4
       w(ix^S,rho_)   = rho4
       w(ix^S,mom(1)) = u4
       w(ix^S,mom(2)) = v4
    end where

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine rm2d_init_one_grid

end module mod_usr
