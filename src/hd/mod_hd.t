!> Module containing all hydrodynamics
module mod_hd
  use mod_hd_phys
  use mod_hd_hllc
  use mod_hd_roe
  use mod_hd_eos
  use mod_eos, only: eos

  use mod_amrvac

  implicit none
  public

contains

  subroutine hd_activate()
    call hd_phys_init()
    call hd_hllc_init()
    call hd_roe_init()
    call hd_link_eos()
  end subroutine hd_activate

end module mod_hd
