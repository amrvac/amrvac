module mod_usr
  use mod_rd

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => gray_scott_init

    call rd_activate()

  end subroutine usr_init

  ! initialize one grid
  subroutine gray_scott_init(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision                :: x1, x2, urand(ix^S)

    x1 = xprobmin1 + 0.5d0 * (xprobmax1 - xprobmin1)
    x2 = xprobmin2 + 0.5d0 * (xprobmax2 - xprobmin2)

    call random_number(urand)

    where (abs(x(ix^S, 1) - x1) < 0.1d0 .and. abs(x(ix^S, 2) - x2) < 0.1d0)
       w(ix^S,u_) = 1.0d-1 * (urand - 0.5d0) + 0.5d0
       w(ix^S,v_) = 0.25d0
    elsewhere
       w(ix^S,u_) = 1.0d0
       w(ix^S,v_) = 0.0d0
    endwhere

  end subroutine gray_scott_init

end module mod_usr
