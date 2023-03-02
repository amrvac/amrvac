module mod_mhd
  use mod_mhd_phys
  use mod_mhd_hllc
  use mod_mhd_roe

  use mod_amrvac

  implicit none
  public

contains

  subroutine mhd_activate()
    call mhd_phys_init()
    call mhd_hllc_init()
    call mhd_roe_init()
  end subroutine mhd_activate

end module mod_mhd
