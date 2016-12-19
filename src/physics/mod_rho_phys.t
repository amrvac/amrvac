module mod_rho_phys
  use mod_physics

  implicit none
  private

  integer, parameter, public          :: rho_       = 1
  double precision, protected, public :: rho_v(^ND) = 1.0d0

  ! Public methods
  public :: rho_phys_init

contains

  subroutine rho_params_read(par_files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: par_files(:)
    integer                      :: n

    namelist /rho_list/ rho_v

    do n = 1, size(par_files)
       open(unitpar, file=trim(par_files(n)), status='old')
       read(unitpar, rho_list, end=111)
111    close(unitpar)
    end do

  end subroutine rho_params_read

  subroutine rho_phys_init(par_files)
    use mod_global_parameters

    character(len=*), intent(in) :: par_files(:)

    call rho_params_read(par_files)

    physics_type = "rho"

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
    phys_add_source_geom => rho_add_source_geom

  end subroutine rho_phys_init

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
    double precision, intent(inout), optional :: cmin(ixG^T)

    call phys_get_v(w, x, ixI^L, ixO^L, idim, cmax)

    if (present(cmin)) then
       cmin(ixO^S) = min(cmax(ixO^S), zero)
       cmax(ixO^S) = max(cmax(ixO^S), zero)
    else
       cmax(ixO^S) = abs(cmax(ixO^S))
    end if
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

  subroutine rho_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)

    ! Add geometrical source terms to w
    ! There are no geometrical source terms in the transport equation

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)

  end subroutine rho_add_source_geom

end module mod_rho_phys
