module mod_rho
  use mod_physics

  implicit none
  private

  integer, parameter, public :: rho_ = 1

  double precision :: rho_v(^ND)

  ! Public methods
  public :: rho_activate

contains

  subroutine rho_activate
    use mod_global_parameters
    physics_type = "rho"

    b0_  = -1 ! No magnetic field
    b^C_ = -1 ! No magnetic field
    e_   = -1 ! No energy (compilation of convert)

    nwflux       = 1
    nwaux        = 0
    nwextra      = 0
    nw           = nwflux + nwaux + nwextra
    nvector      = 0
    nworkroe     = 1

    nflag_ = nw+1

    phys_get_v           => rho_get_v
    phys_get_cmax        => rho_get_cmax
    phys_get_flux        => rho_get_flux
    phys_average         => rho_average
    phys_get_eigenjump   => rho_get_eigenjump
    phys_rtimes          => rho_rtimes
    phys_add_source_geom => rho_add_source_geom
  end subroutine rho_activate

  subroutine rho_read_physics_params(file_unit, io_state)
    integer, intent(in) :: file_unit
    integer, intent(out) :: io_state

    namelist /rho_list/ rho_v

    read(file_unit, rho_list, iostat=io_state)
  end subroutine rho_read_physics_params

  subroutine rho_get_v(w, x, ixI^L, ixO^L, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(out) :: v(ixG^T)

    v(ixO^S) = rho_v(idim)
  end subroutine rho_get_v

  subroutine rho_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixG^T)
    double precision, optional, intent(inout) :: cmin(ixG^T)

    call get_v(w, x, ixI^L, ixO^L, idim, cmax)

    if (present(cmin)) then
       cmin(ixO^S)=min(cmax(ixO^S), zero)
       cmax(ixO^S)=max(cmax(ixO^S), zero)
    else
       cmax(ixO^S)=abs(cmax(ixO^S))
    endif

  end subroutine rho_get_cmax

  ! There is nothing to add to the transport flux in the transport equation
  subroutine rho_get_flux(w, x, ixI^L, ixO^L, iw, idim, f, transport)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, iw, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    double precision, intent(inout) :: f(ixG^T)
    logical, intent(out)            :: transport

    f(ixO^S)  = zero
    transport = .true.
  end subroutine rho_get_flux

  subroutine rho_average(wL, wR, x, ix^L, idim, wroe, workroe)
    use mod_global_parameters
    integer, intent(in)                          :: ix^L, idim
    double precision, dimension(ixG^T, nw)       :: wL, wR, wroe
    double precision, dimension(ixG^T, nworkroe) :: workroe
    double precision, dimension(ixG^T, 1:^ND)   :: x

    wroe(ix^S, rho_)=half*(wL(ix^S, rho_)+wR(ix^S, rho_))
  end subroutine rho_average

  subroutine rho_get_eigenjump(wL, wR, wC, x, ix^L, il, idim, smalla, a, jump, workroe)

    ! Calculate the characteristic speed a and the jump in the
    ! characteristic variable in the idim direction within ixL.
    ! For a scalar equation the characteristic and conservative variables coincide
    ! The characteristic speed is just the velocity, but it should be averaged
    ! for the cell interfaces

    use mod_global_parameters

    integer                                      :: ix^L, jx^L, ixC^L, il, idim
    double precision, dimension(ixG^T, nw)       :: wL, wR, wC
    double precision, dimension(ixG^T)           :: smalla, a, jump, v
    double precision, dimension(ixG^T, nworkroe) :: workroe
    double precision, dimension(ixG^T, 1:^ND)   :: x

    jx^L=ix^L+kr(idim,^D);
    ixCmin^D=ixmin^D; ixCmax^D=jxmax^D;

    ! No entropy fix
    smalla(ix^S)= -one
    ! The velocity is independent of w in the transport equation,
    ! but it may depend on the location
    call get_v(wL, x, ixG^LL, ixC^L, idim, v)

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

  subroutine rho_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)

    ! Add geometrical source terms to w
    ! There are no geometrical source terms in the transport equation

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)

  end subroutine rho_add_source_geom

end module mod_rho
