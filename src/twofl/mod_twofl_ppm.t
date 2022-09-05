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

    if(twofl_eq_energy/=0) then
      drho(ixO^S) =twofl_gamma*dabs(d2w(ixO^S,rho_c_))&
           /min(w(ixL^S,rho_c_),w(ixR^S,rho_c_))
      dp(ixO^S) = dabs(d2w(ixO^S,e_c_))/min(w(ixL^S,e_c_),w(ixR^S,e_c_))
    else
      call mpistop("PPM with flatcd=.true. can not be used without energy equation!")
    end if
  end subroutine twofl_ppm_flatcd

  ! based on Mignone and Miller and Collela 2002
  ! PPM flattening at shocks: we use total pressure and not thermal pressure
  subroutine twofl_ppm_flatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(inout) :: drho(ixI^S),dp(ixI^S),dv(ixI^S)
    double precision                :: ptot(ixI^S)

    if(twofl_eq_energy/=0) then
      ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
      ptot(ixO^S)=w(ixO^S,e_c_)+half*sum(w(ixO^S,mag(:))**2,dim=ndim+1)
      where (dabs(ptot(ixRR^S)-ptot(ixLL^S))>smalldouble)
         drho(ixO^S) = dabs((ptot(ixR^S)-ptot(ixL^S))&
              /(ptot(ixRR^S)-ptot(ixLL^S)))
      elsewhere
         drho(ixO^S) = zero
      end where

      !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
      !  use "dp" to save squared sound speed, assume primitive in w
      dp(ixO^S)=(twofl_gamma*w(ixO^S,e_c_)/w(ixO^S,rho_c_))

      dp(ixO^S)  = dabs(ptot(ixR^S)-ptot(ixL^S))&
           /(w(ixO^S,rho_c_)*dp(ixO^S))
      ! recycle ptot to store v
      ptot(ixI^S)= w(ixI^S,mom_c(idims))
      call gradient(ptot,ixI^L,ixO^L,idims,dv)
    else
      call mpistop("PPM with flatsh=.true. can not be used without energy equation!")
    end if

  end subroutine twofl_ppm_flatsh

end module mod_twofl_ppm
