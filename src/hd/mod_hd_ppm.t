!> Hydrodynamics PPM module
module mod_hd_ppm
  use mod_hd_phys

  implicit none
  private

  public :: hd_ppm_init

contains

  subroutine hd_ppm_init()
    use mod_physics_ppm

    phys_ppm_flatcd => hd_ppm_flatcd
    phys_ppm_flatsh => hd_ppm_flatsh
  end subroutine hd_ppm_init

  subroutine hd_ppm_flatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, ixL^L, ixR^L
    double precision, intent(in)    :: w(ixI^S, nw), d2w(ixG^T, 1:nwflux)
    double precision, intent(inout) :: drho(ixG^T), dp(ixG^T)

    drho(ixO^S) = hd_gamma*abs(d2w(ixO^S, rho_))&
         /min(w(ixL^S, rho_), w(ixR^S, rho_))
    dp(ixO^S) = abs(d2w(ixO^S, e_))/min(w(ixL^S, e_), w(ixR^S, e_))
  end subroutine hd_ppm_flatcd

  subroutine hd_ppm_flatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S, nw)
    double precision, intent(inout) :: drho(ixI^S), dp(ixI^S)

    ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
    where (abs(w(ixRR^S, p_)-w(ixLL^S, p_))>smalldouble)
       drho(ixO^S) = abs((w(ixR^S, p_)-w(ixL^S, p_))&
            /(w(ixRR^S, p_)-w(ixLL^S, p_)))
    else where
       drho(ixO^S) = zero
    end where

    ! eq. 76, page 48, Miller and Colella 2002, JCoPh
    dp(ixO^S) = abs(w(ixR^S, p_)-w(ixL^S, p_))&
         /(hd_gamma*w(ixO^S, p_))

  end subroutine hd_ppm_flatsh

end module mod_hd_ppm
