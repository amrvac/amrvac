module mod_usr
  use mod_ard
  implicit none
contains

  subroutine usr_init()
    usr_init_one_grid => brusselator_init
    call set_coordinate_system('Cartesian')
    call ard_activate()
  end subroutine usr_init

  subroutine brusselator_init(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision                :: x1, x2, l1, l2, urand(ix^S)

    x1 = xprobmin1 + 0.5d0 * (xprobmax1 - xprobmin1)
    x2 = xprobmin2 + 0.5d0 * (xprobmax2 - xprobmin2)
    l1 = xprobmax1 - xprobmin1
    l2 = xprobmax2 - xprobmin2

    select case (iprob)
    case(1)
        ! Initialize with homogeneous steady-state (A, B/A)
        ! perturbed by random noise (urand)
        call random_number(urand)
        w(ix^S, u_) = br_A      + 1.0d-3*urand
        call random_number(urand)
        w(ix^S, v_) = br_B/br_A + 1.0d-3*urand
    case(2)
        ! Initialize with homogeneous steady-state (A, B/A)
        ! perturbed in central square
        w(ix^S,u_) = br_A - 1.0d-1
        w(ix^S,v_) = br_B - 1.0d-1
        where (abs(x(ix^S, 1) - x1) < 0.2d0 * l1 .and. &
             abs(x(ix^S, 2) - x2) < 0.2d0 * l2)
            w(ix^S,u_) = br_A + 1.0d-3
            w(ix^S,v_) = br_B/br_A + 1.0d-3
        endwhere
    case(3)
        ! Initialize with homogeneous steady-state (A, B/A)
        ! perturbed in central circle
        w(ix^S,u_) = br_A - 1.0d-1
        w(ix^S,v_) = br_B - 1.0d-1
        where ((x(ix^S, 1) - x1)**2 + (x(ix^S, 2) - x2)**2 < 0.2d0 * (l1 * l2)**0.5)
            w(ix^S,u_) = br_A + 1.0d-3
            w(ix^S,v_) = br_B/br_A + 1.0d-3
        endwhere
    case default
        call mpistop("Unknown iprob")
    end select

    
  end subroutine brusselator_init

end module mod_usr
