module mod_srmhd_ppm
  use mod_srmhd_parameters
  use mod_srmhd_phys
  implicit none
  private

  public :: srmhd_ppm_init

contains

  subroutine srmhd_ppm_init()
    use mod_physics_ppm

    phys_ppm_flatcd => srmhd_ppm_flatcd
    phys_ppm_flatsh => srmhd_ppm_flatsh
  end subroutine srmhd_ppm_init

  subroutine srmhd_ppm_flatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dpressure)
  ! made by Z. MELIANI 14/02/2018
    use mod_global_parameters
    use mod_srmhd_phys

    integer, intent(in)             :: ixI^L,ixO^L,ixL^L,ixR^L
    double precision, intent(in)    :: w(ixI^S,nw),d2w(ixI^S,1:nwflux)
    double precision, intent(inout) :: drho(ixI^S),dpressure(ixI^S)

    double precision                ::  pmag(ixI^S),dpressuremag(ixI^S)

    if(srmhd_energy) then
      call srmhd_get_p_mag(w,ixI^L,ixI^L,pmag)
      dpressuremag(ixO^S) = pmag(ixL^S)-pmag(ixR^S)

      drho(ixO^S) =dabs(d2w(ixO^S,rho_))&
           /min(w(ixL^S,rho_),w(ixR^S,rho_))

      dpressure(ixO^S) = dabs(d2w(ixO^S,p_)+dpressuremag(ixO^S))&
                             /min(w(ixL^S,p_),w(ixR^S,p_))
    else
      call mpistop("PPM with flatcd=.true. can not be used without &
                   energy equation!")
    end if
  end subroutine srmhd_ppm_flatcd

  ! based on Mignone and Miller and Collela 2002
  ! PPM flattening at shocks: we use total pressure and not thermal pressure
  subroutine srmhd_ppm_flatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,&
                              w,drho,dpressure,dv)
  ! made by Z. MELIANI 14/02/2018
    use mod_global_parameters
    use mod_srmhd_phys
    use mod_geometry

    integer, intent(in)             :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(inout) :: drho(ixI^S),dpressure(ixI^S),dv(ixI^S)
    double precision                :: ptot(ixI^S)

    if(srmhd_energy) then
      ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
      call srmhd_get_p_total(w,pw(saveigrid)%x,ixI^L,ixI^L,ptot)
      where (dabs(ptot(ixRR^S)-ptot(ixLL^S))>smalldouble)
         drho(ixO^S) = dabs((ptot(ixR^S)-ptot(ixL^S))&
              /(ptot(ixRR^S)-ptot(ixLL^S)))
      elsewhere
         drho(ixO^S) = zero
      end where

      !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
      !  use "dpressure" to save squared sound speed, assume primitive in w
      call srmhd_get_csound_prim(w,pw(saveigrid)%x,ixI^L,ixO^L,idims,dpressure)

      dpressure(ixO^S)  = dabs(ptot(ixR^S)-ptot(ixL^S))&
           /(w(ixO^S,rho_)*dpressure(ixO^S))
      ! recycle ptot to store v
      ptot(ixI^S)= w(ixI^S,mom(idims))/w(ixI^S,lfac_)
      call gradient(ptot,ixI^L,ixO^L,idims,dv)
    else
      call mpistop("PPM with flatsh=.true. can not be used &
                    without energy equation!")
    end if

  end subroutine srmhd_ppm_flatsh

end module mod_srmhd_ppm
