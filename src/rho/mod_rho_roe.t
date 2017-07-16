!> Module containing Roe solver for scalar advection
module mod_rho_roe
  use mod_physics_roe
  use mod_rho_phys

  implicit none
  private

  public :: rho_roe_init

contains

  subroutine rho_roe_init()
    use mod_physics_roe

    nworkroe = 1

    phys_average         => rho_average
    phys_get_eigenjump   => rho_get_eigenjump
    phys_rtimes          => rho_rtimes
  end subroutine rho_roe_init

  subroutine rho_average(wL, wR, x, ix^L, idim, wroe, workroe)
    use mod_global_parameters
    integer, intent(in)             :: ix^L, idim
    double precision, intent(in)    :: wL(ixG^T, nw), wR(ixG^T, nw)
    double precision, intent(inout) :: wroe(ixG^T, nw)
    double precision, intent(inout) :: workroe(ixG^T, nworkroe)
    double precision, intent(in)    :: x(ixG^T, 1:^ND)

    wroe(ix^S, rho_)=half*(wL(ix^S, rho_)+wR(ix^S, rho_))
  end subroutine rho_average

  subroutine rho_get_eigenjump(wL, wR, wC, x, ix^L, il, idim, smalla, a, jump, workroe)

    ! Calculate the characteristic speed a and the jump in the
    ! characteristic variable in the idim direction within ixL.
    ! For a scalar equation the characteristic and conservative variables coincide
    ! The characteristic speed is just the velocity, but it should be averaged
    ! for the cell interfaces

    use mod_global_parameters

    integer, intent(in)                          :: ix^L, il, idim
    double precision, dimension(ixG^T, nw)       :: wL, wR, wC
    double precision, dimension(ixG^T)           :: smalla, a, jump, v
    double precision, dimension(ixG^T, nworkroe) :: workroe
    double precision, intent(in)                 :: x(ixG^T, 1:^ND)
    integer                                      :: jx^L, ixC^L

    jx^L=ix^L+kr(idim,^D);
    ixCmin^D=ixmin^D; ixCmax^D=jxmax^D;

    ! No entropy fix
    smalla(ix^S)= -one
    ! The velocity is independent of w in the transport equation,
    ! but it may depend on the location
    call rho_get_v(wL, x, ixG^LL, ixC^L, idim, v)

    a(ix^S)=(v(jx^S)+v(ix^S))/2

    jump(ix^S)=wR(ix^S, rho_)-wL(ix^S, rho_)

  end subroutine rho_get_eigenjump

  subroutine rho_rtimes(q, w, ix^L, iw, il, idim, rq, workroe)

    ! Multiply q by R(il, iw), where R is the right eigenvalue matrix at wC.
    ! For a scalar equation the R matrix is unity

    use mod_global_parameters
    integer, intent(in)             :: ix^L, iw, il, idim
    double precision, intent(in)    :: w(ixG^T, nw), q(ixG^T)
    double precision, intent(inout) :: rq(ixG^T)
    double precision, intent(inout) :: workroe(ixG^T, nworkroe)

    rq(ix^S)=q(ix^S)

  end subroutine rho_rtimes

end module mod_rho_roe
