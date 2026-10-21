module mod_mhd
  use mod_mhd_phys
  use mod_functions_bfield, only: mag
  use mod_mhd_hllc
  use mod_mhd_roe
  use mod_mhd_eos
  use mod_eos, only: eos

  use mod_amrvac

  implicit none
  public

contains

  subroutine mhd_activate()
    call mhd_phys_init()
    call mhd_hllc_init()
    call mhd_roe_init()
    call mhd_link_eos()
  end subroutine mhd_activate

end module mod_mhd
