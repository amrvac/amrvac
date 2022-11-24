!> Module activating the advection-reaction-diffusion equation module
module mod_ard
  use mod_ard_phys
  use mod_amrvac

  implicit none
  public

contains

  subroutine ard_activate()
    call ard_phys_init()
  end subroutine ard_activate

end module mod_ard
