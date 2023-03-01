module mod_twofl_ppm

#include "amrvac.h"
  use mod_twofl_phys

  implicit none
  private

  public :: twofl_ppm_init

contains

  subroutine twofl_ppm_init()
    use mod_physics_ppm
    phys_ppm_flatcd => twofl_ppm_flatcd
    phys_ppm_flatsh => twofl_ppm_flatsh
  end subroutine twofl_ppm_init

  subroutine twofl_ppm_flatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L,ixL^L,ixR^L
    double precision, intent(in)    :: w(ixI^S,nw),d2w(ixI^S,1:nwflux)
    double precision, intent(inout) :: drho(ixI^S),dp(ixI^S)

    drho(ixO^S) =twofl_gamma*dabs(d2w(ixO^S,iw_rho))&
         /min(w(ixL^S,iw_rho),w(ixR^S,iw_rho))
    dp(ixO^S) = dabs(d2w(ixO^S,iw_e))/min(w(ixL^S,iw_e),w(ixR^S,iw_e))

  end subroutine twofl_ppm_flatcd

  ! based on Mignone and Miller and Collela 2002
  ! PPM flattening at shocks: we use total pressure and not thermal pressure
  subroutine twofl_ppm_flatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(inout) :: drho(ixI^S),dp(ixI^S)

    where (dabs(w(ixRR^S,iw_e)-w(ixLL^S,iw_e))>smalldouble)
       drho(ixO^S) = dabs((w(ixR^S,iw_e)-w(ixL^S,iw_e))&
            /(w(ixRR^S,iw_e)-w(ixLL^S,iw_e)))
    elsewhere
       drho(ixO^S) = zero
    end where

    ! eq. 76, page 48, Miller and Colella 2002, JCoPh, adjusted by testing
    dp(ixO^S) = abs(w(ixR^S, iw_e)-w(ixL^S, iw_e))&
         /(twofl_gamma*0.8d0*(w(ixR^S, iw_e)+w(ixL^S, iw_e)))

  end subroutine twofl_ppm_flatsh

end module mod_twofl_ppm
