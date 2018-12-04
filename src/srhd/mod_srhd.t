module mod_srhd
  use mod_srhd_phys
  use mod_srhd_hllc
  use mod_srhd_roe
  use mod_srhd_eos

  implicit none
  public

contains

  subroutine srhd_activate()
    call srhd_phys_init()
    call srhd_hllc_init()
    call srhd_roe_init()
  end subroutine srhd_activate

end module mod_srhd
