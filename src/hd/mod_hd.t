!> Module containing all hydrodynamics
module mod_hd
  use mod_hd_phys

  use mod_amrvac

  implicit none
  public

contains

  subroutine hd_activate()
    call hd_phys_init()
  end subroutine hd_activate

end module mod_hd
