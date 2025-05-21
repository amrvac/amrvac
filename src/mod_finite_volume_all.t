!> Module with finite volume methods for fluxes

#:include 'hd/mod_hd_templates.fpp'

module mod_finite_volume_all

  use mod_variables
  use mod_hd_phys, only: hd_gamma
  use mod_global_parameters, only: ndim
  use mod_physicaldata
  implicit none

  private

  public :: finite_volume_local

  integer, parameter :: dp = kind(0.0d0), nw_euler=2+ndim

contains

! instantiate the templated functions here:
@:addsource_local()
@:to_primitive()
@:to_conservative()
@:get_cmax()
@:get_flux()

  subroutine finite_volume_local(method, qdt, dtfactor, ixI^L, ixO^L, idims^LIM, &
       qtC, bga, qt, bgb, fC, fE)
    use mod_global_parameters

    integer, intent(in)                                                :: method
    double precision, intent(in)                                       :: qdt, dtfactor, qtC, qt
    integer, intent(in)                                                :: ixI^L, ixO^L, idims^LIM
    ! remember, old names map as: wCT => bga, wnew => bgb
    type(block_grid_t)                                    :: bga
    type(block_grid_t)                                    :: bgb
    double precision, dimension(ixI^S, 1:nwflux, 1:ndim)  :: fC ! not yet provided
    double precision, dimension(ixI^S, sdim:3)            :: fE ! not yet provided
    ! .. local ..
    integer                :: n, iigrid, ix^D
    double precision       :: uprim(nw, ixI^S)
    real(dp)               :: tmp(nw_euler,5)
    real(dp)               :: f(nw_euler, 2)
    real(dp)               :: inv_dr(ndim)
    integer                :: typelim
    real(dp)               :: xloc(ndim)
    real(dp)               :: wprim(nw_euler), wCT(nw_euler), wnew(nw_euler)
    !-----------------------------------------------------------------------------

    !$acc parallel loop private(n, uprim, inv_dr, typelim) firstprivate(ixI^L, ixO^L) present(bga%w, bgb%w)
    do iigrid = 1, igridstail_active
       n = igrids_active(iigrid)

       inv_dr  = 1/rnode(rpdx1_:rnodehi, n)
       typelim = type_limiter(node(plevel_, n))

       !$acc loop collapse(ndim) vector
       {^D& do ix^DB=ixImin^DB,ixImax^DB \}
             ! Convert to primitive
             uprim(:, ix^D) = bga%w(ix^D, :, n)
             call to_primitive(uprim(:, ix^D))
       {^D& end do \}

       !$acc loop collapse(ndim) private(f, tmp, wnew, wCT, xloc, wprim) vector
       {^D& do ix^DB=ixOmin^DB,ixOmax^DB \}
             ! Compute fluxes in all dimensions

       {^IFONED
             tmp = uprim(:, ix1-2:ix1+2)
             call muscl_flux_prim(tmp, 1, f, typelim)
             ! Update the wnew array
             bgb%w(ix1, :, n) = bgb%w(ix1, :, n) + qdt * &
                  ( (f(:, 1) - f(:, 2)) * inv_dr(1) )
       }
       {^IFTWOD
             tmp = uprim(:, ix1-2:ix1+2, ix2)
             call muscl_flux_prim(tmp, 1, f, typelim)
             bgb%w(ix1, ix2, :, n) = bgb%w(ix1, ix2, :, n) &
                  + qdt * (f(:, 1) - f(:, 2)) * inv_dr(1)

             tmp = uprim(:, ix1, ix2-2:ix2+2)
             call muscl_flux_prim(tmp, 2, f, typelim)
             bgb%w(ix1, ix2, :, n) = bgb%w(ix1, ix2, :, n) &
                  + qdt * (f(:, 1) - f(:, 2)) * inv_dr(2)
       }
       {^IFTHREED
             tmp = uprim(:, ix1-2:ix1+2, ix2, ix3)
             call muscl_flux_prim(tmp, 1, f, typelim)
             bgb%w(ix1, ix2, ix3, :, n) = bgb%w(ix1, ix2, ix3, :, n) &
                  + qdt * (f(:, 1) - f(:, 2)) * inv_dr(1)

             tmp = uprim(:, ix1, ix2-2:ix2+2, ix3)
             call muscl_flux_prim(tmp, 2, f, typelim)
             bgb%w(ix1, ix2, ix3, :, n) = bgb%w(ix1, ix2, ix3, :, n) &
                  + qdt * (f(:, 1) - f(:, 2)) * inv_dr(2)

             tmp = uprim(:, ix1, ix2, ix3-2:ix3+2)
             call muscl_flux_prim(tmp, 3, f, typelim)
             bgb%w(ix1, ix2, ix3, :, n) = bgb%w(ix1, ix2, ix3, :, n) &
                  + qdt * (f(:, 1) - f(:, 2)) * inv_dr(3)

#:if defined('GRAVITY')
             ! Add source terms:
             xloc(1:ndim) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
             wprim        = uprim(1:nw_euler, ix1, ix2, ix3)
             wCT          = bga%w(ix1, ix2, ix3, 1:nw_euler, n)
             wnew         = bgb%w(ix1, ix2, ix3, 1:nw_euler, n)
             call addsource_local(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
                  dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), & 
                  qtC, wCT, wprim, qt, wnew, xloc, .false. )
             bgb%w(ix1, ix2, ix3, :, n) = wnew(:)
#:endif             
       }
       
       {^D& end do \}
    end do

  end subroutine finite_volume_local

  subroutine muscl_flux_prim(u, flux_dim, flux, typelim)
    !$acc routine seq

    use mod_limiter, only: limiter_minmod, limiter_vanleer

    real(dp), intent(in)  :: u(nw_euler, 5)
    integer, intent(in)   :: flux_dim, typelim
    real(dp), intent(out) :: flux(nw_euler, 2)
    real(dp)              :: uL(nw_euler), uR(nw_euler), wL, wR, wmax
    real(dp)              :: flux_l(nw_euler), flux_r(nw_euler)

    ! Construct uL, uR for first cell face
    select case (typelim)
    case (limiter_minmod)
       uL = u(:, 2) + 0.5_dp * minmod(u(:, 2) - u(:, 1), u(:, 3) - u(:, 2))
       uR = u(:, 3) - 0.5_dp * minmod(u(:, 3) - u(:, 2), u(:, 4) - u(:, 3))
    case (limiter_vanleer)
       uL = u(:, 2) + 0.5_dp * vanleer(u(:, 2) - u(:, 1), u(:, 3) - u(:, 2))
       uR = u(:, 3) - 0.5_dp * vanleer(u(:, 3) - u(:, 2), u(:, 4) - u(:, 3))
    end select

    call get_flux(uL, flux_dim, flux_l)
    call get_flux(uR, flux_dim, flux_r)

    wL = get_cmax(uL, flux_dim)
    wR = get_cmax(uR, flux_dim)
    wmax = max(wL, wR)

    call to_conservative(uL)
    call to_conservative(uR)
    flux(:, 1) = 0.5_dp * ((flux_l + flux_r) - wmax * (uR - uL))

    ! Construct uL, uR for second cell face
    select case (typelim)
    case (limiter_minmod)
       uL = u(:, 3) + 0.5_dp * minmod(u(:, 3) - u(:, 2), u(:, 4) - u(:, 3))
       uR = u(:, 4) - 0.5_dp * minmod(u(:, 4) - u(:, 3), u(:, 5) - u(:, 4))
    case (limiter_vanleer)
       uL = u(:, 3) + 0.5_dp * vanleer(u(:, 3) - u(:, 2), u(:, 4) - u(:, 3))
       uR = u(:, 4) - 0.5_dp * vanleer(u(:, 4) - u(:, 3), u(:, 5) - u(:, 4))
    end select

    call get_flux(uL, flux_dim, flux_l)
    call get_flux(uR, flux_dim, flux_r)

    wL = get_cmax(uL, flux_dim)
    wR = get_cmax(uR, flux_dim)
    wmax = max(wL, wR)

    call to_conservative(uL)
    call to_conservative(uR)
    flux(:, 2) = 0.5_dp * ((flux_l + flux_r) - wmax * (uR - uL))

  end subroutine muscl_flux_prim

  elemental pure real(dp) function vanleer(a, b) result(phi)
    !$acc routine seq
    real(dp), intent(in) :: a, b
    real(dp)             :: ab

    ab = a * b
    if (ab > 0) then
       phi = 2 * ab / (a + b)
    else
       phi = 0
    end if
  end function vanleer

  elemental pure real(dp) function minmod(a, b)
    !$acc routine seq
    real(dp), intent(in) :: a, b

    if (a * b <= 0) then
       minmod = 0.0_dp
    else if (abs(a) < abs(b)) then
       minmod = a
    else
       minmod = b
    end if
  end function minmod

end module mod_finite_volume_all
