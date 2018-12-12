module mod_srmhd
  use mod_srmhd_phys
  use mod_srmhd_hllc
  use mod_srmhd_roe
  use mod_srmhd_ppm

  use mod_amrvac

  implicit none
  public

contains

  subroutine srmhd_activate()
    call srmhd_phys_init()
    call srmhd_hllc_init()
    call srmhd_roe_init()
    call srmhd_ppm_init()
  end subroutine srmhd_activate

end module mod_srmhd
