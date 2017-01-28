module mod_rho_phys

  implicit none
  private

  integer, parameter, public          :: rho_       = 1
  double precision, protected, public :: rho_v(^ND) = 1.0d0

  ! Public methods
  public :: rho_phys_init
  public :: rho_get_v

contains

  subroutine rho_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /rho_list/ rho_v

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, rho_list, end=111)
111    close(unitpar)
    end do

  end subroutine rho_params_read

  subroutine rho_phys_init()
    use mod_global_parameters
    use mod_physics

    call rho_params_read(par_files)

    physics_type = "rho"

    nwflux       = 1
    nwaux        = 0
    nwextra      = 0
    nw           = nwflux + nwaux + nwextra
    nvector      = 0
    nworkroe     = 1

    nflag_ = nw+1

    phys_get_cmax        => rho_get_cmax
    phys_get_flux        => rho_get_flux
    phys_add_source_geom => rho_add_source_geom
    phys_to_conserved    => rho_to_conserved
    phys_to_primitive    => rho_to_primitive
    phys_get_dt          => rho_get_dt
  end subroutine rho_phys_init

  subroutine rho_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for rho module)
  end subroutine rho_to_conserved

  subroutine rho_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for rho module)
  end subroutine rho_to_primitive

  subroutine rho_get_v(w, x, ixI^L, ixO^L, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(out) :: v(ixI^S)

    v(ixO^S) = rho_v(idim)
  end subroutine rho_get_v

  subroutine rho_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    call rho_get_v(w, x, ixI^L, ixO^L, idim, cmax)

    if (present(cmin)) then
       cmin(ixO^S) = min(cmax(ixO^S), zero)
       cmax(ixO^S) = max(cmax(ixO^S), zero)
    else
       cmax(ixO^S) = abs(cmax(ixO^S))
    end if
  end subroutine rho_get_cmax

  subroutine rho_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble
  end subroutine rho_get_dt

  ! There is nothing to add to the transport flux in the transport equation
  subroutine rho_get_flux(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    double precision, intent(out)   :: f(ixI^S, nwflux)
    double precision                :: v(ixI^S)

    call rho_get_v(w, x, ixI^L, ixO^L, idim, v)

    f(ixO^S, rho_) = w(ixO^S, rho_) * v(ixO^S)
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
