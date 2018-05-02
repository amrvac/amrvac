!> Module activating the nonlinear scalar equation module
!> allowing to solve Burgers, nonconvex, and KdV equations
module mod_nonlinear
  use mod_nonlinear_phys
  use mod_nonlinear_roe

  use mod_amrvac

  implicit none
  public

contains

  subroutine nonlinear_activate()
    call nonlinear_phys_init()
    call nonlinear_roe_init()
  end subroutine nonlinear_activate

end module mod_nonlinear
