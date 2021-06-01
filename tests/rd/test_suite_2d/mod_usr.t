module mod_usr

  use mod_rd
  implicit none

  integer, parameter :: dp = kind(0.0d0)

contains

  subroutine usr_init()
     integer :: i

     call set_coordinate_system('Cartesian')
     call rd_activate()

     usr_init_one_grid => gs_init
  end subroutine usr_init

  subroutine gs_init(ixG^L,ix^L,w,x)
     integer, intent(in)             :: ixG^L, ix^L
     double precision, intent(in)    :: x(ixG^S,1:ndim)
     double precision, intent(inout) :: w(ixG^S,1:nw)
     double precision :: x1, x2, l1, l2

     ! Utility variables
     x1 = xprobmin1 + 0.5d0 * (xprobmax1 - xprobmin1)
     l1 = xprobmax1 - xprobmin1
     x2 = xprobmin2 + 0.5d0 * (xprobmax2 - xprobmin2)
     l2 = xprobmax2 - xprobmin2

     ! Center interval
     w(ix^S,u_) = 1.0d0
     w(ix^S,v_) = 0.0d0
     where ( (abs(x(ix^S, 1) - x1) < 0.1d0 * l1) .and. &
             (abs(x(ix^S, 2) - x2) < 0.1d0 * l2) )
        w(ix^S,u_) = 0.5d0
        w(ix^S,v_) = 0.25d0
     endwhere
  end subroutine gs_init

end module mod_usr
