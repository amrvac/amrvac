!> Module containing all hydrodynamics
module mod_rhd
  use mod_rhd_phys
  use mod_rhd_hllc
  use mod_rhd_roe
  use mod_rhd_ppm

  use mod_amrvac

  implicit none
  public

contains

  subroutine rhd_activate()
    call rhd_phys_init()
    call rhd_hllc_init()
    call rhd_roe_init()
    call rhd_ppm_init()
  end subroutine rhd_activate

end module mod_rhd
