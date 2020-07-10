!> Module containing all routines for reaction-diffusion
module mod_rd
  use mod_rd_phys
  use mod_amrvac

  implicit none
  public

contains

  subroutine rd_activate()
    call rd_phys_init()
  end subroutine rd_activate

end module mod_rd
