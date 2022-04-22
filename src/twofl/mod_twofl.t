module mod_twofl

#include "amrvac.h"

  use mod_twofl_phys
  use mod_twofl_hllc
#if defined (ONE_FLUID) && ONE_FLUID ==1
  use mod_twofl_roe
  use mod_twofl_ppm
#endif
  use mod_amrvac

  implicit none
  public

contains

  subroutine twofl_activate()
    call twofl_phys_init()
    call twofl_hllc_init()

#if defined (ONE_FLUID) && ONE_FLUID ==1
    call twofl_roe_init()
    call twofl_ppm_init()
#endif
  end subroutine twofl_activate

end module mod_twofl
