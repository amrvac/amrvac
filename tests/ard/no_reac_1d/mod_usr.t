module mod_usr
  use mod_ard

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => no_reac_init

    call ard_activate()

  end subroutine usr_init

  ! initialize one grid
  subroutine no_reac_init(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision                :: x1, x2, urand(ix^S)
    double precision                :: l1, l2, dist2(ix^S)
    logical                         :: mymask(ix^S)

    x1 = xprobmin1 + 0.5d0 * (xprobmax1 - xprobmin1)
    l1 = xprobmax1 - xprobmin1

    ! Default: steady state
    w(ix^S,u_) = 1.0d1

    call random_number(urand)

    select case (iprob)
    case (1)
       ! Center square
       w(ix^S,u_) = 0.0d0
       where (abs(x(ix^S, 1) - x1) < 0.1d0 * l1)
          w(ix^S,u_) = 1.0d0
       endwhere
    case (2)
       ! Center square with random noise
       where (abs(x(ix^S, 1) - x1) < 0.1d0 * l1)
          w(ix^S,u_) = 1.0d-1 * (urand - 0.5d0) + 0.5d0
       endwhere
    case (3)
       ! Gaussian
       dist2 = (x(ix^S, 1) - x1)**2
       w(ix^S,u_) = w(ix^S,u_) - 0.5d0 * exp(-100 * dist2/l1**2)
    case (4)
        ! One cosine peak
        w(ix^S,u_) = 0.0d0
        where (abs(x(ix^S, 1) - x1) < 2 * atan(1.0d0))
          w(ix^S,u_) = 1.0d0 * cos(x(ix^S, 1)-x1)
       endwhere
    case (5)
       ! Triangle right
       w(ix^S,u_) = 0.0d0
       where (abs(x(ix^S, 1) - x1) < 0.1d0 * l1)
          w(ix^S,u_) = (x(ix^S,1) - x1 + 0.1d0 * l1) / (0.1d0 * l1)
       endwhere
    case (6)
       ! Triangle left
       w(ix^S,u_) = 0.0d0
       where (abs(x(ix^S, 1) - x1) < 0.1d0 * l1)
          w(ix^S,u_) = 2d0 - (x(ix^S,1) - x1 + 0.1d0 * l1) / (0.1d0 * l1)
       endwhere
    case default
       call mpistop("Unknown iprob")
    end select

  end subroutine no_reac_init

end module mod_usr
