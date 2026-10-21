module mod_ffhd
  use mod_ffhd_phys
  use mod_ffhd_eos
  use mod_amrvac

  implicit none
  public
contains

  subroutine ffhd_activate()
    call ffhd_phys_init()
    call ffhd_link_eos()
  end subroutine ffhd_activate
end module mod_ffhd
