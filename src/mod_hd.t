module mod_hd
  use mod_physics

  implicit none
  private

  logical :: hd_energy
  integer :: hd_n_dust
  integer :: hd_n_tracer

  integer, parameter :: rho_ = 1
  integer, parameter :: m0_  = rho_
  {integer, parameter :: m^C_ = m0_+^C\}

  ! Public methods
  public :: rho_activate

contains

  subroutine rho_activate
    use mod_global_parameters

    physics_type = 'hd'
    nwflux       = ^NC*(1+^NDS)+2+^NDS+^NFL
    if (.not. hd_energy) nwflux = nwflux - 1

    integer, parameter :: nwaux   = 0
    integer, parameter :: nwextra = 0
    integer, parameter :: nw      = nwflux+nwaux+nwextra

    if (hd_energy) then
       integer, parameter:: e_=m^NC_+1
       integer, parameter:: ee_=e_
       integer, parameter:: rhos_=e_
       if (hd_n_tracer > 0) then
          INTEGER, PARAMETER:: Dtr^FL_=e_+^FL
          if (hd_n_dust > 0) then
             INTEGER, PARAMETER:: rhod0_=e_+^NFL
          end if
       end if
       if (hd_tracer <= 0) then
          if (hd_n_dust > 0) then
             INTEGER, PARAMETER:: rhod0_=e_
          end if
       end if
    end if

    if (.not. hd_energy) then
       if (hd_n_tracer > 0) then
          INTEGER, PARAMETER:: Dtr^FL_=m^NC_+^FL
          if (hd_n_dust > 0) then
             INTEGER, PARAMETER:: rhod0_=m^NC_+^NFL
          end if
       end if
       if (hd_n_tracer <= 0) then
          if (hd_n_dust > 0) then
             INTEGER, PARAMETER:: rhod0_=m^NC_
          end if
       end if
    end if

    if (hd_n_dust > 0) then
       integer, parameter:: rhod^DS_=rhod0_+^DS
       integer, parameter:: m0d^DS_=rhod^DS_
       integer, parameter:: m1d^DS_=rhod^NDS_+^DS
       {^NOONEC
       integer, parameter:: m2d^DS_=m1d^NDS_+^DS }
       {^IFTHREEC
       integer, parameter:: m3d^DS_=m2d^NDS_+^DS }
    end if

    integer, parameter:: b0_=-1 ! No magnetic field
    integer, parameter:: b^C_=-1 ! No magnetic field

    if (hd_energy) then
       INTEGER, PARAMETER:: p_=e_, pp_=ee_    ! Primitive variables
    end if

    INTEGER, PARAMETER:: v0_=m0_, v^C_=m^C_

    if (hd_n_tracer > 0) then
       integer, parameter:: tr^FL_=Dtr^FL_
    end if

    INTEGER, PARAMETER:: mr_=m0_+r_, mphi_=m0_+phi_, mz_=m0_+z_  ! Polar var. names
    if (hd_n_dust > 0) then
       integer, parameter:: v1d^DS_=rhod^NDS_+^DS
       {^NOONEC
       integer, parameter:: v2d^DS_=m1d^NDS_+^DS }
       {^IFTHREEC
       integer, parameter:: v3d^DS_=m2d^NDS_+^DS }
       INTEGER, PARAMETER:: mrd^DS_=m0d^DS_+(^NDS*r_)              ! Polar var. names
       INTEGER, PARAMETER:: mphid^DS_=m0d^DS_+(^NDS*phi_)          ! Polar var. names
       INTEGER, PARAMETER:: mzd^DS_=m0d^DS_+(^NDS*z_)              ! Polar var. names
    end if

    ! Note: dust related velocity vectors not handled here
    integer, parameter :: nvector=1                             ! No. vector vars
    integer, dimension(nvector), parameter :: iw_vector=(/ m0_ /)

    if (hd_energy) then
       ! Characteristic
       INTEGER, PARAMETER:: soundRW_=1, soundLW_=2, entropW_=3, shearW0_=3      ! waves
       INTEGER, PARAMETER:: nworkroe=3
    end if

    if (hd_iso) then
       ! Characteristic
       INTEGER, PARAMETER:: soundRW_=1, soundLW_=2, shearW0_=2        ! waves
       INTEGER, PARAMETER:: nworkroe=1
    end if

    if (hd_energy) then
       INTEGER, PARAMETER:: gamma_=1, mu_=2, neqpar=2                     ! equation parameters
    end if
    if (hd_iso) then
       INTEGER, PARAMETER:: gamma_=1, adiab_=2, mu_=3, neqpar=3                     ! equation parameters
    end if
    INTEGER, PARAMETER:: nflag_=nw+1
    COMMON, INTEGER:: flags(nflag_)
    COMMON, DOUBLE PRECISION:: wflags(nflag_)

    ! xprob: problem box; iprob: problem
    COMMON, INTEGER:: iprob
    COMMON, DOUBLE PRECISION:: xprob^L

    COMMON, DOUBLE PRECISION::{#IFDEF ENERGY smalle, } minrho, minp
    {#IFDEF DUST
    DOUBLE PRECISION:: minrhod, sdust(1:^NDS), dsdust(1:^NDS), rhodust(1:^NDS), mhcgspar, kbcgspar
    DOUBLE PRECISION :: Lstar, Tdust
    COMMON /DPDUST/ minrhod, sdust, dsdust, rhodust, mhcgspar, kbcgspar, Tdust, Lstar
    }

    b0_  = -1 ! No magnetic field
    {b^C_ = -1\} ! No magnetic field
    e_   = -1 ! No energy (compilation of convert)

    nwflux       = 1
    nwaux        = 0
    nwextra      = 0
    nw           = nwflux + nwaux + nwextra
    nvector      = 0
    nworkroe     = 1

    nflag_ = nw+1

    phys_read_params     => rho_read_params
    phys_get_v           => rho_get_v
    phys_get_cmax        => rho_get_cmax
    phys_get_flux        => rho_get_flux
    phys_average         => rho_average
    phys_get_eigenjump   => rho_get_eigenjump
    phys_rtimes          => rho_rtimes
    phys_add_source_geom => rho_add_source_geom
  end subroutine rho_activate

  subroutine rho_read_params(file_unit)
    integer, intent(in) :: file_unit

    namelist /rho_list/ rho_v

    read(file_unit, rho_list)
  end subroutine rho_read_params

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

    call phys_get_v(w, x, ixI^L, ixO^L, idim, cmax)

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
    call phys_get_v(wL, x, ixG^LL, ixC^L, idim, v)

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

end module mod_hd
