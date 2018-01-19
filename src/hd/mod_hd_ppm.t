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

  subroutine hd_ppm_flatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dpr)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, ixL^L, ixR^L
    double precision, intent(in)    :: w(ixI^S, nw), d2w(ixG^T, 1:nwflux)
    double precision, intent(inout) :: drho(ixG^T), dpr(ixG^T)

    if(hd_energy) then
      drho(ixO^S) = hd_gamma*abs(d2w(ixO^S, rho_))&
           /min(w(ixL^S, rho_), w(ixR^S, rho_))
      dpr(ixO^S) = abs(d2w(ixO^S, e_))/min(w(ixL^S, e_), w(ixR^S, e_))
    else
      call mpistop("PPM with flatcd=.true. can not be used with -eos = iso !")
    end if
  end subroutine hd_ppm_flatcd

  subroutine hd_ppm_flatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dpr,dv)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S, nw)
    double precision, intent(inout) :: drho(ixG^T), dpr(ixG^T), dv(ixG^T)
    double precision                :: v(ixG^T)

    if(hd_energy) then
      ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
      where (abs(w(ixRR^S, e_)-w(ixLL^S, e_))>smalldouble)
         drho(ixO^S) = abs((w(ixR^S, e_)-w(ixL^S, e_))&
              /(w(ixRR^S, e_)-w(ixLL^S, e_)))
      else where
         drho(ixO^S) = zero
      end where

      !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
      !  use "dpr" to save squared sound speed, assuming primitives
      dpr(ixO^S) =(hd_gamma*w(ixO^S, e_)/w(ixO^S, rho_))

      dpr(ixO^S) = abs(w(ixR^S, e_)-w(ixL^S, e_))&
           /(w(ixO^S, rho_)*dpr(ixO^S))
      v(ixI^S)  = w(ixI^S, mom(idims))
      call gradient(v, ixI^L, ixO^L, idims, dv)
    else
      call mpistop("PPM with flatsh=.true. can not be used with -eos = iso !")
    end if
  end subroutine hd_ppm_flatsh

end module mod_hd_ppm
