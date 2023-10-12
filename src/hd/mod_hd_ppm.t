!> hydrodynamics PPM module
module mod_hd_ppm
  use mod_physics

  implicit none
  private

  ! Public methods
  public :: hd_ppm_init

contains

  !> Initialize the module
  subroutine hd_ppm_init()
    use mod_physics_ppm
    phys_ppm_flatcd => hd_ppm_flatcd
    phys_ppm_flatsh => hd_ppm_flatsh
  end subroutine hd_ppm_init

  subroutine hd_ppm_flatcd(ixI^L, ixO^L, ixL^L, ixR^L, w, d2w, drho, dp)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, ixL^L, ixR^L
    double precision, intent(in)    :: w(ixI^S, 1:nw),d2w(ixI^S,1:nw)
    double precision, intent(inout) :: drho(ixI^S), dp(ixI^S)

    drho(ixO^S) = abs(d2w(ixO^S,rho_)) &
                  / min( w(ixL^S, rho_), w(ixR^S, rho_) )

    !dp(ixO^S) = abs(d2w(ixO^S,press_)) &
    !              / min( w(ixL^S, press_), w(ixR^S, press_) )
  end subroutine hd_ppm_flatcd

  subroutine hd_ppm_flatsh(ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L, idims, w, beta, z, dv)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: beta(ixI^S), z(ixI^S), dv(ixI^S)

    ! beta here is the shock width
    ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
    !where ( abs( w(ixRR^S, press_) - w(ixLL^S, press_) )> smalldouble )
    !   beta(ixO^S) = abs( w(ixR^S, press_) - w(ixL^S, press_) ) &
    !                 /abs( w(ixRR^S, press_) - w(ixLL^S, press_) )
    !else where
    !   beta(ixO^S) = 0.0d0
    !end where

    ! Z here is the shock strength
    ! in eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
    ! where the denominator is rho * cs2, which is gamma * P in ideal gas.
    ! However, in eq. B16, page 218, Mignone and Bodo 2005, ApJS (beta1),
    ! this should be min(p(i-1), p(i+1)).
    ! Here we follow the later approach, as usually we don't have gamma 
    ! in general
    !z(ixO^S) = abs( w(ixR^S, press_) - w(ixL^S, press_) ) &
    !                 /min( w(ixR^S, press_), w(ixL^S, press_) )

    dv(ixO^S)=w(ixR^S,W_vel(idims)) - w(ixL^S,W_vel(idims))

  end subroutine hd_ppm_flatsh

end module mod_hd_ppm
