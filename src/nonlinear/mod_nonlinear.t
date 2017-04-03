module mod_nonlinear
  use mod_physics
  use mod_nonlinear_params

  implicit none
  public

  integer, parameter :: rho_                = 1
  integer            :: nonlinear_flux_type = 1

  ! Public methods
  ! public :: nonlinear_activate

contains

  subroutine nonlinear_params_read(par_files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: par_files(:)
    integer                      :: n

    namelist /nonlinear_list/ nonlinear_flux_type

    do n = 1, size(par_files)
       open(unitpar, file=trim(par_files(n)), status='old')
       read(unitpar, nonlinear_list)
       close(unitpar)
    end do

  end subroutine nonlinear_params_read

  subroutine nonlinear_activate(par_files)
    use mod_global_parameters

    character(len=*), intent(in) :: par_files(:)

    call nonlinear_params_read(par_files)

    physics_type = "nonlinear"

    nwflux       = 1
    nwaux        = 0
    nwextra      = 0
    nw           = nwflux + nwaux + nwextra
    nvector      = 0
    nworkroe     = 1

    phys_get_v           => nonlinear_get_v
    phys_get_cmax        => nonlinear_get_cmax
    phys_get_flux        => nonlinear_get_flux
    phys_add_source_geom => nonlinear_add_source_geom
  end subroutine nonlinear_activate

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
       call mpistop('Undefined fluxtype: set nonlinear_flux_type to 1 or 2')
    end select

  end subroutine nonlinear_get_v

  subroutine nonlinear_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixG^T)
    double precision, intent(inout), optional :: cmin(ixG^T)

    call phys_get_v(w, x, ixI^L, ixO^L, idim, cmax)

    if (present(cmin)) then
       cmin(ixO^S) = min(cmax(ixO^S), zero)
       cmax(ixO^S) = max(cmax(ixO^S), zero)
    else
       cmax(ixO^S) = abs(cmax(ixO^S))
    end if
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
       call mpistop('Undefined fluxtype: set nonlinear_flux_type to 1 or 2')
    end select

    transport = .false.

  end subroutine nonlinear_get_flux

  subroutine nonlinear_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)

    ! Add geometrical source terms to w
    ! There are no geometrical source terms in the transport equation

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)

  end subroutine nonlinear_add_source_geom

end module mod_nonlinear
