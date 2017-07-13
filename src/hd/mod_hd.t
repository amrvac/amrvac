!> Module containing all hydrodynamics
module mod_hd
  use mod_hd_phys
  use mod_hd_hllc
  use mod_hd_roe
  use mod_hd_ppm

  implicit none
  public

contains

  subroutine hd_activate()
    call hd_phys_init()
    call hd_hllc_init()
    call hd_roe_init()
    call hd_ppm_init()
  end subroutine hd_activate

end module mod_hd
