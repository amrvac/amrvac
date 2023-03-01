module mod_mhd_ppm
  use mod_mhd_phys

  implicit none
  private

  public :: mhd_ppm_init

contains

  subroutine mhd_ppm_init()
    use mod_physics_ppm

    phys_ppm_flatcd => mhd_ppm_flatcd
    phys_ppm_flatsh => mhd_ppm_flatsh
  end subroutine mhd_ppm_init

  subroutine mhd_ppm_flatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L,ixL^L,ixR^L
    double precision, intent(in)    :: w(ixI^S,nw),d2w(ixI^S,1:nwflux)
    double precision, intent(inout) :: drho(ixI^S),dp(ixI^S)

    drho(ixO^S) =mhd_gamma*dabs(d2w(ixO^S,rho_))&
         /min(w(ixL^S,rho_),w(ixR^S,rho_))
    dp(ixO^S) = dabs(d2w(ixO^S,p_))/min(w(ixL^S,p_),w(ixR^S,p_))

  end subroutine mhd_ppm_flatcd

  ! based on Mignone and Miller and Collela 2002
  ! PPM flattening at shocks: we use total pressure and not thermal pressure
  subroutine mhd_ppm_flatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(inout) :: drho(ixI^S),dp(ixI^S)

    ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
    where (dabs(w(ixRR^S,p_)-w(ixLL^S,p_))>smalldouble)
       drho(ixO^S) = dabs((w(ixR^S,p_)-w(ixL^S,p_))&
            /(w(ixRR^S,p_)-w(ixLL^S,p_)))
    elsewhere
       drho(ixO^S) = zero
    end where

    ! eq. 76, page 48, Miller and Colella 2002, JCoPh, adjusted by testing
    dp(ixO^S) = abs(w(ixR^S, p_)-w(ixL^S, p_))&
         /(mhd_gamma*(w(ixR^S, p_)+w(ixL^S, p_)))

  end subroutine mhd_ppm_flatsh

end module mod_mhd_ppm
