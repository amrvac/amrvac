module mod_nonlinear
  use mod_physics

  implicit none
  private

  integer, parameter, public :: rho_ = 1

  integer :: nonlinear_flux_type = 1

  ! Public methods
  public :: nonlinear_activate

contains

  subroutine nonlinear_activate
    use mod_global_parameters

    physics_type = "nonlinear"

    nwflux       = 1
    nwaux        = 0
    nwextra      = 0
    nw           = nwflux + nwaux + nwextra
    nvector      = 0
    nworkroe     = 1

    nflag_ = nw+1

    phys_read_params     => nonlinear_read_params
    phys_get_v           => nonlinear_get_v
    phys_get_cmax        => nonlinear_get_cmax
    phys_get_flux        => nonlinear_get_flux
    phys_average         => nonlinear_average
    phys_get_eigenjump   => nonlinear_get_eigenjump
    phys_rtimes          => nonlinear_rtimes
    phys_add_source_geom => nonlinear_add_source_geom
  end subroutine nonlinear_activate

  subroutine nonlinear_read_params(file_unit, success)
    integer, intent(in) :: file_unit
    logical, intent(out) :: success

    namelist /nonlinear_list/ nonlinear_flux_type

    ! Success indicates if read succeeds
    success = .false.
    read(file_unit, nonlinear_list)
    success = .true.

    101 return
  end subroutine nonlinear_read_params

  subroutine nonlinear_get_v(w, x, ixI^L, ixO^L, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(out) :: v(ixG^T)

    select case(nonlinear_flux_type)
    case(1)
       v(ixO^S)=w(ixO^S,rho_)
    case(2)
       v(ixO^S)=3.0d0*w(ixO^S,rho_)**2
    case default
       call mpistop('Undefined fluxtype: set eqpar to 1 or 2')
    end select

  end subroutine nonlinear_get_v

  subroutine nonlinear_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixG^T)
    double precision, optional, intent(inout) :: cmin(ixG^T)

    call phys_get_v(w, x, ixI^L, ixO^L, idim, cmax)

    if (present(cmin)) then
       cmin(ixO^S)=min(cmax(ixO^S), zero)
       cmax(ixO^S)=max(cmax(ixO^S), zero)
    else
       cmax(ixO^S)=abs(cmax(ixO^S))
    endif

  end subroutine nonlinear_get_cmax

  ! There is nothing to add to the transport flux in the transport equation
  subroutine nonlinear_get_flux(w, x, ixI^L, ixO^L, iw, idim, f, transport)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, iw, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    double precision, intent(inout) :: f(ixG^T)
    logical, intent(out)            :: transport

    select case(nonlinear_flux_type)
    case(1)
       f(ixO^S)=half*w(ixO^S,rho_)**2
    case(2)
       f(ixO^S)=w(ixO^S,rho_)**3
    case default
       call mpistop('Undefined fluxtype: set eqpar to 1 or 2')
    end select

    transport = .false.

  end subroutine nonlinear_get_flux

  subroutine nonlinear_average(wL, wR, x, ix^L, idim, wroe, workroe)
    use mod_global_parameters
    integer, intent(in)                          :: ix^L, idim
    double precision, dimension(ixG^T, nw)       :: wL, wR, wroe
    double precision, dimension(ixG^T, nworkroe) :: workroe
    double precision, dimension(ixG^T, 1:^ND)    :: x

    call mpistop("Error: nonlinear_average: not implemented!")
  end subroutine nonlinear_average

  subroutine nonlinear_get_eigenjump(wL, wR, wC, x, ix^L, il, idim, smalla, a, jump, workroe)

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
    double precision, dimension(ixG^T, 1:^ND)    :: x

    call mpistop("Error: nonlinear_geteigenjump: no roe solver implemented!")

  end subroutine nonlinear_get_eigenjump

  subroutine nonlinear_rtimes(q, w, ix^L, iw, il, idim, rq, workroe)

    ! Multiply q by R(il, iw), where R is the right eigenvalue matrix at wC.
    ! For a scalar equation the R matrix is unity

    use mod_global_parameters
    integer, intent(in)             :: ix^L, iw, il, idim
    double precision, intent(in)    :: w(ixG^T, nw), q(ixG^T)
    double precision, intent(inout) :: rq(ixG^T)
    double precision, intent(inout) :: workroe(ixG^T, nworkroe)

    call mpistop("Error: nonlinear_rtimes: no roe solver implemented!")

  end subroutine nonlinear_rtimes

  subroutine nonlinear_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)

    ! Add geometrical source terms to w
    ! There are no geometrical source terms in the transport equation

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)

  end subroutine nonlinear_add_source_geom

end module mod_nonlinear
