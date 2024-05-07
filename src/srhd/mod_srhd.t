!> Module containing all special relativistic hydrodynamics (with EOS)
module mod_srhd
  use mod_srhd_phys

  use mod_amrvac

  implicit none
  public

contains

  subroutine srhd_activate()
    call srhd_phys_init()
  end subroutine srhd_activate

end module mod_srhd
