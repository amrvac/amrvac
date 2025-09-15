!> Module with finite volume methods for fluxes

#:mute
#:include "physics/mod_physics_templates.fpp"
#:endmute

module mod_finite_volume

  use mod_variables
  use mod_physics_vars
  use mod_global_parameters, only: ndim
  use mod_physicaldata
  implicit none

  private

  public :: finite_volume_local

contains

! instantiate the templated functions here for inlining:
@:addsource_local()
@:addsource_nonlocal()
@:to_primitive()
@:to_conservative()
@:get_cmax()
@:get_flux()

  subroutine finite_volume_local(method, qdt, dtfactor, ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, idimsmin,idimsmax, qtC, bga, qt, bgb, fC, fE)
    use mod_global_parameters

    integer, intent(in)                                                :: &
       method
    double precision, intent(in)                                       :: qdt,&
        dtfactor, qtC, qt
    integer, intent(in)                                                :: &
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3, idimsmin,idimsmax
    ! remember, old names map as: wCT => bga, wnew => bgb
    type(block_grid_t)                                    :: bga
    type(block_grid_t)                                    :: bgb
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nwflux, 1:ndim)  :: fC !not yet provided
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, sdim:3)            :: fE !not yet provided
    ! .. local ..
    integer                :: n, iigrid, ix1,ix2,ix3
    double precision       :: uprim(nw, ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    real(dp)               :: tmp(nw_phys,5)
    real(dp)               :: f(nw_phys, 2)
    real(dp)               :: inv_dr(ndim)
    real(dp)               :: dr(ndim)
    integer                :: typelim
    real(dp)               :: xloc(ndim)
    real(dp)               :: xlocC(ndim,2)
    real(dp)               :: wprim(nw_phys), wCT(nw_phys), wnew(nw_phys)
    !-----------------------------------------------------------------------------
    
    !$acc parallel loop private(uprim, inv_dr, dr) present(bga%w, bgb%w)
    do iigrid = 1, igridstail_active
       n = igrids_active(iigrid)

       dr  = rnode(rpdx1_:rnodehi, n)
       inv_dr  = 1/rnode(rpdx1_:rnodehi, n)
       typelim = type_limiter(node(plevel_, n))

       !$acc loop collapse(ndim) vector
       do ix3=ixImin3,ixImax3 
          do ix2=ixImin2,ixImax2 
             do ix1=ixImin1,ixImax1 
                ! Convert to primitive
                uprim(:, ix1,ix2,ix3) = bga%w(ix1,ix2,ix3, :, n)
                call to_primitive(uprim(:, ix1,ix2,ix3))
             end do
          end do
       end do

       !$acc loop collapse(ndim) private(f, wnew, tmp, xlocC, xloc #{if defined('SOURCE_LOCAL')}#, wCT, wprim #{endif}#) vector
       do ix3=ixOmin3,ixOmax3 
          do ix2=ixOmin2,ixOmax2 
             do ix1=ixOmin1,ixOmax1 
                ! Compute fluxes in all dimensions

                tmp = uprim(:, ix1-2:ix1+2, ix2, ix3)
                xlocC(1:ndim,1) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(1:ndim,2) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(1,1) = xlocC(1,1)-0.5_dp*dr(1)
                xlocC(1,2) = xlocC(1,2)+0.5_dp*dr(1)
                call muscl_flux_prim(tmp, xlocC, 1, f, typelim)
                bgb%w(ix1, ix2, ix3, :, n) = bgb%w(ix1, ix2, ix3, :,&
                     n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(1)

                tmp = uprim(:, ix1, ix2-2:ix2+2, ix3)
                xlocC(1:ndim,1) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(1:ndim,2) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(2,1) = xlocC(2,1)-0.5_dp*dr(2)
                xlocC(2,2) = xlocC(2,2)+0.5_dp*dr(2)
                call muscl_flux_prim(tmp, xlocC, 2, f, typelim)
                bgb%w(ix1, ix2, ix3, :, n) = bgb%w(ix1, ix2, ix3, :,&
                     n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(2)

                tmp = uprim(:, ix1, ix2, ix3-2:ix3+2)
                xlocC(1:ndim,1) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(1:ndim,2) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(3,1) = xlocC(3,1)-0.5_dp*dr(3)
                xlocC(3,2) = xlocC(3,2)+0.5_dp*dr(3)
                call muscl_flux_prim(tmp, xlocC, 3, f, typelim)
                bgb%w(ix1, ix2, ix3, :, n) = bgb%w(ix1, ix2, ix3, :,&
                     n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(3)

#:if defined('SOURCE_LOCAL')
                ! Add local source terms:
                xloc(1:ndim) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                wprim        = uprim(1:nw_phys, ix1, ix2, ix3)
                wCT          = bga%w(ix1, ix2, ix3, 1:nw_phys, n)
                wnew         = bgb%w(ix1, ix2, ix3, 1:nw_phys, n)
                call addsource_local(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                     dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, wCT,&
                     wprim, qt, wnew, xloc, dr, .false. )
                bgb%w(ix1, ix2, ix3, :, n) = wnew(:)
#:endif             

#:if defined('SOURCE_NONLOCAL')
                ! Add non-local (gradient) source terms:
                xloc(1:ndim) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                wnew         = bgb%w(ix1, ix2, ix3, 1:nw_phys, n)
                
                tmp = uprim(:, ix1-2:ix1+2, ix2, ix3)
                call addsource_nonlocal(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                     dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, tmp,&
                     qt, wnew, xloc, dr, 1, .false. )

                tmp = uprim(:, ix1, ix2-2:ix2+2, ix3)
                call addsource_nonlocal(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                     dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, tmp,&
                     qt, wnew, xloc, dr, 2, .false. )

                tmp = uprim(:, ix1, ix2, ix3-2:ix3+2)
                call addsource_nonlocal(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                     dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, tmp,&
                     qt, wnew, xloc, dr, 3, .false. )
                
                bgb%w(ix1, ix2, ix3, :, n) = wnew(:)           
#:endif                
             end do
          end do
       end do
    end do

  end subroutine finite_volume_local

  subroutine muscl_flux_prim(u, xlocC, flux_dim, flux, typelim)
    !$acc routine seq

    use mod_limiter, only: limiter_minmod, limiter_vanleer

    real(dp), intent(in)  :: u(nw_phys, 5)
    real(dp), intent(in)  :: xlocC(1:ndim, 2)
    integer, intent(in)   :: flux_dim, typelim
    real(dp), intent(out) :: flux(nw_phys, 2)
    real(dp)              :: uL(nw_phys), uR(nw_phys), wL, wR, wmax
    real(dp)              :: flux_l(nw_phys), flux_r(nw_phys)
    real(dp)              :: xC(ndim)

    ! Construct uL, uR for first cell face
    select case (typelim)
    case (limiter_minmod)
       uL = u(:, 2) + 0.5_dp * minmod(u(:, 2) - u(:, 1), u(:, 3) - u(:, 2))
       uR = u(:, 3) - 0.5_dp * minmod(u(:, 3) - u(:, 2), u(:, 4) - u(:, 3))
    case (limiter_vanleer)
       uL = u(:, 2) + 0.5_dp * vanleer(u(:, 2) - u(:, 1), u(:, 3) - u(:, 2))
       uR = u(:, 3) - 0.5_dp * vanleer(u(:, 3) - u(:, 2), u(:, 4) - u(:, 3))
    end select

    xC=xlocC(:,1)
    call get_flux(uL, xC, flux_dim, flux_l)
    call get_flux(uR, xC, flux_dim, flux_r)

    wL = get_cmax(uL,xC, flux_dim)
    wR = get_cmax(uR,xC, flux_dim)
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

    xC=xlocC(:,2)
    call get_flux(uL, xC, flux_dim, flux_l)
    call get_flux(uR, xC, flux_dim, flux_r)

    wL = get_cmax(uL,xC, flux_dim)
    wR = get_cmax(uR,xC, flux_dim)
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

end module mod_finite_volume
