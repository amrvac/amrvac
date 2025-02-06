module mod_rmhd
  use mod_rmhd_phys
  use mod_functions_bfield, only: mag
  use mod_amrvac
  implicit none
  public

contains

  subroutine rmhd_activate()
    call rmhd_phys_init()
  end subroutine rmhd_activate

end module mod_rmhd
