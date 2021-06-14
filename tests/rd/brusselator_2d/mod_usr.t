module mod_usr
  use mod_rd
  implicit none
contains

  subroutine usr_init()
    usr_init_one_grid => brusselator_init
    call set_coordinate_system('Cartesian')
    call rd_activate()
  end subroutine usr_init

  subroutine brusselator_init(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision                :: urand(ix^S)

    ! Initialize with homogeneous steady-state (A, B/A)
    ! perturbed by a noise component (urand)
    call random_number(urand)
    w(ix^S, u_) = br_A      + 1.0d-3*urand
    call random_number(urand)
    w(ix^S, v_) = br_B/br_A + 1.0d-3*urand
  end subroutine brusselator_init

end module mod_usr
