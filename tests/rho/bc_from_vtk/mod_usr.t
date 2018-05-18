module mod_usr
  use mod_rho

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => initial_conditions

    call set_coordinate_system("Cartesian")
    call rho_activate()
  end subroutine usr_init

  subroutine initial_conditions(ixG^L, ix^L, w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, 1:ndim)
    double precision, intent(inout) :: w(ixG^S, 1:nw)

    w(ix^S, rho_) = 0.0d0
  end subroutine initial_conditions

end module mod_usr
