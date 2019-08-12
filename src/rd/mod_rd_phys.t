!> Reaction-diffusion module (physics routines)
module mod_rd_phys

  implicit none
  private

  integer, protected, public :: u_ = 1
  integer, protected, public :: v_ = 2

  integer            :: equation_type   = 1
  integer, parameter :: eq_gray_scott   = 1
  integer, parameter :: eq_schnakenberg = 2

  !> Diffusion coefficient for first species (u)
  double precision, public, protected :: D1 = 0.05d0
  !> Diffusion coefficient for second species (v)
  double precision, public, protected :: D2 = 1.0d0

  !> Parameter for Schnakenberg model
  double precision, public, protected :: sb_alpha = 0.1305d0
  !> Parameter for Schnakenberg model
  double precision, public, protected :: sb_beta  = 0.7695d0
  !> Parameter for Schnakenberg model
  double precision, public, protected :: sb_kappa = 100.0d0

  !> Parameter for Gray-Scott model
  double precision, public, protected :: gs_F = 0.046d0
  !> Parameter for Gray-Scott model
  double precision, public, protected :: gs_k = 0.063d0

  !> Whether to handle the source term in split fashion
  logical :: rd_source_split = .false.

  ! Public methods
  public :: rd_phys_init

contains

  subroutine rd_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n
    character(len=20)            :: equation_name

    namelist /rd_list/ D1, D2, sb_alpha, sb_beta, sb_kappa, gs_F, gs_k, &
         equation_name

    equation_name = "gray-scott"

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, rd_list, end=111)
111    close(unitpar)
    end do

    select case (equation_name)
    case ("gray-scott")
       equation_type = eq_gray_scott
    case ("schnakenberg")
       equation_type = eq_schnakenberg
    case default
       call mpistop("Unknown equation_name (not gray-scott or schnakenberg)")
    end select

  end subroutine rd_params_read

  !> Write this module's parameters to a snapsoht
  subroutine rd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 0
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er
    integer                             :: idim

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)
  end subroutine rd_write_info

  subroutine rd_phys_init()
    use mod_global_parameters
    use mod_physics

    call rd_params_read(par_files)

    physics_type = "rd"
    phys_energy  = .false.
    phys_req_diagonal = .false.

    ! Use the first variable as a density
    u_ = var_set_fluxvar("u", "u")
    v_ = var_set_fluxvar("v", "v")

    phys_get_cmax     => rd_get_cmax
    phys_get_cbounds  => rd_get_cbounds
    phys_get_flux     => rd_get_flux
    phys_to_conserved => rd_to_conserved
    phys_to_primitive => rd_to_primitive
    phys_add_source   => rd_add_source
    phys_get_dt       => rd_get_dt
    phys_write_info   => rd_write_info
    phys_check_params => rd_check_params

  end subroutine rd_phys_init

  subroutine rd_check_params
    use mod_global_parameters

    if (any(flux_scheme /= "source")) then
       call mpistop("mod_rd requires flux_scheme = source")
    end if
  end subroutine rd_check_params

  subroutine rd_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for rd module)
  end subroutine rd_to_conserved

  subroutine rd_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for rd module)
  end subroutine rd_to_primitive

  subroutine rd_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixI^S)

    cmax(ixO^S) = 0.0d0
  end subroutine rd_get_cmax

  subroutine rd_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S,nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    if (present(cmin)) then
       cmin(ixO^S) = 0.0d0
       cmax(ixO^S) = 0.0d0
    else
       cmax(ixO^S) = 0.0d0
    end if

  end subroutine rd_get_cbounds

  subroutine rd_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    ! dt < dx^2 / (2 * ndim * diffusion_coeff)
    dtnew = 0.9d0 * minval([ dx^D ])**2 / (2 * ndim * max(D1, D2))

    ! TODO: set time step for reactions
  end subroutine rd_get_dt

  ! There is nothing to add to the transport flux in the transport equation
  subroutine rd_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wC(ixI^S, 1:nw)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(out)   :: f(ixI^S, nwflux)
    double precision                :: v(ixI^S)

    f(ixO^S, :) = 0.0d0
  end subroutine rd_get_flux

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine rd_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision                :: lpl_u(ixO^S), lpl_v(ixO^S)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    if (qsourcesplit .eqv. rd_source_split) then
       call rd_laplacian(ixI^L, ixO^L, wCT(ixI^S, u_), lpl_u)
       call rd_laplacian(ixI^L, ixO^L, wCT(ixI^S, v_), lpl_v)

       select case (equation_type)
       case (eq_gray_scott)
          w(ixO^S, u_) = w(ixO^S, u_) + qdt * (D1 * lpl_u - &
               wCT(ixO^S, u_) * wCT(ixO^S, v_)**2 + &
               gs_F * (1 - wCT(ixO^S, u_)))
          w(ixO^S, v_) = w(ixO^S, v_) + qdt * (D2 * lpl_v + &
               wCT(ixO^S, u_) * wCT(ixO^S, v_)**2 - &
               (gs_F + gs_k) * wCT(ixO^S, v_))
       case (eq_schnakenberg)
          w(ixO^S, u_) = w(ixO^S, u_) + qdt * (D1 * lpl_u &
               + sb_kappa * (sb_alpha - wCT(ixO^S, u_) + &
               wCT(ixO^S, u_)**2 * wCT(ixO^S, v_)))
          w(ixO^S, v_) = w(ixO^S, v_) + qdt * (D2 * lpl_v &
               + sb_kappa * (sb_beta - wCT(ixO^S, u_)**2 * wCT(ixO^S, v_)))
       case default
          call mpistop("Unknown equation type")
       end select

       active = .true.
    end if

  end subroutine rd_add_source

  !> Compute the Laplacian using a standard second order scheme. For now this
  !> method only works in slab geometries.
  subroutine rd_laplacian(ixI^L,ixO^L,var,lpl)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: var(ixI^S)
    double precision, intent(out) :: lpl(ixO^S)
    integer                       :: idir, jxO^L, hxO^L
    double precision              :: h_inv2

    if (slab) then
       lpl(ixO^S) = 0.0d0

       do idir = 1, ndim
          hxO^L=ixO^L-kr(idir,^D);
          jxO^L=ixO^L+kr(idir,^D);
          h_inv2 = 1/dxlevel(idir)**2
          lpl(ixO^S) = lpl(ixO^S) + h_inv2 * &
               (var(jxO^S) - 2 * var(ixO^S) + var(hxO^S))
       end do
    else
       call mpistop("rd_laplacian not implemented in this geometry")
    end if
  end subroutine rd_laplacian

end module mod_rd_phys
