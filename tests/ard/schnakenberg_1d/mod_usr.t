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
    double precision                :: x1
    double precision                :: l1

    x1 = xprobmin1 + (xprobmax1 - xprobmin1)/2
    l1 = xprobmax1 - xprobmin1

    ! Initialize to steady state solution
    w(ix^S,u_) = sb_alpha + sb_beta
    w(ix^S,v_) = sb_beta / (sb_alpha + sb_beta)**2

    select case (iprob)
    case (1)
       ! Gaussian perturbation of u
       w(ix^S,u_) = w(ix^S,u_) + 1d-1 * &
            exp(-100d0 * (x(ix^S, 1) - x1)**2)
    case (2)
       ! Central square perturbation of u and v
       where (abs(x(ix^S, 1) - x1) < 0.1d0 * l1)
          w(ix^S,u_) = w(ix^S,u_) + 1.0d-1
          w(ix^S,v_) = w(ix^S,v_) - 1.0d-1
       endwhere
    case default
       call mpistop("Unknown iprob")
    end select

  end subroutine schnakenberg_init

end module mod_usr
