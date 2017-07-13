!> Module containing all routines for scalar advection
module mod_rho
  use mod_rho_phys
  use mod_rho_roe

  implicit none
  public

contains

  subroutine rho_activate()
    call rho_phys_init()
    call rho_roe_init()
  end subroutine rho_activate

end module mod_rho
