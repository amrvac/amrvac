!> Hydrodynamics PPM module
module mod_rhd_ppm
  use mod_rhd_phys

  implicit none
  private

  public :: rhd_ppm_init

contains

  subroutine rhd_ppm_init()
    use mod_physics_ppm

    phys_ppm_flatcd => rhd_ppm_flatcd
    phys_ppm_flatsh => rhd_ppm_flatsh
  end subroutine rhd_ppm_init

  subroutine rhd_ppm_flatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, ixL^L, ixR^L
    double precision, intent(in)    :: w(ixI^S, nw), d2w(ixG^T, 1:nwflux)
    double precision, intent(inout) :: drho(ixG^T), dp(ixG^T)

    if(rhd_energy) then
      drho(ixO^S) = rhd_gamma*abs(d2w(ixO^S, rho_))&
           /min(w(ixL^S, rho_), w(ixR^S, rho_))
      dp(ixO^S) = abs(d2w(ixO^S, e_))/min(w(ixL^S, e_), w(ixR^S, e_))
    else
      call mpistop("PPM with flatcd=.true. can not be used with -eos = iso !")
    end if
  end subroutine rhd_ppm_flatcd

  subroutine rhd_ppm_flatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S, nw)
    double precision, intent(inout) :: drho(ixG^T), dp(ixG^T), dv(ixG^T)
    double precision                :: v(ixG^T)

    if(rhd_energy) then
      ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
      where (abs(w(ixRR^S, e_)-w(ixLL^S, e_))>smalldouble)
         drho(ixO^S) = abs((w(ixR^S, e_)-w(ixL^S, e_))&
              /(w(ixRR^S, e_)-w(ixLL^S, e_)))
      else where
         drho(ixO^S) = zero
      end where

      !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
      !  use "dp" to save squared sound speed, assuming primitives
      dp(ixO^S) =(rhd_gamma*w(ixO^S, e_)/w(ixO^S, rho_))

      dp(ixO^S) = abs(w(ixR^S, e_)-w(ixL^S, e_))&
           /(w(ixO^S, rho_)*dp(ixO^S))
      v(ixI^S)  = w(ixI^S, mom(idims))
      call gradient(v, ixI^L, ixO^L, idims, dv)
    else
      call mpistop("PPM with flatsh=.true. can not be used with -eos = iso !")
    end if
  end subroutine rhd_ppm_flatsh

end module mod_rhd_ppm
