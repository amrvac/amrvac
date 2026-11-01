module mod_usr
  use mod_ard

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => schnakenberg_init

    call ard_activate()

  end subroutine usr_init

  ! initialize one grid
  subroutine schnakenberg_init(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision                :: x1, x2
    double precision                :: l1, l2

    x1 = xprobmin1 + (xprobmax1 - xprobmin1)/3
    x2 = xprobmin2 + 0.5d0 * (xprobmax2 - xprobmin2)
    l1 = xprobmax1 - xprobmin1
    l2 = xprobmax2 - xprobmin2

    w(ix^S,u_) = sb_alpha + sb_beta + 1d-3 * &
         exp(-100d0 * ((x(ix^S, 1) - x1)**2 + (x(ix^S, 2) - x2)**2))
    w(ix^S,v_) = sb_beta / (sb_alpha + sb_beta)**2

  end subroutine schnakenberg_init

end module mod_usr
