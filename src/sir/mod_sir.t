!> Module containing all routines for SIR model with diffusion
module mod_sir
  use mod_sir_phys
  use mod_amrvac

  implicit none
  public

contains

  subroutine sir_activate()
    call sir_phys_init()
  end subroutine sir_activate

end module mod_sir
