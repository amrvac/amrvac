! - [x] Remove everything we don't use
! - [x] Remove psa, and psb from finite_volume (ps, bga and bgb should be enough)
! - [x] Test for correctness
! - [x] Replace block with ps(igrid), pass igrid to functions that need it
! - [ ] Implement local update alternative from Jannis toycode
!> Module with finite volume methods for fluxes
module mod_finite_volume_all
#include "amrvac.h"
  use mod_variables
  use mod_hd_phys, only: hd_gamma
  use mod_global_parameters, only: ndim
  use mod_physicaldata
  implicit none

  private

  public :: finite_volume_local

  integer, parameter :: dp = kind(0.0d0), nw_euler=2+ndim

contains

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
    real(dp)               :: tmp(nw_euler,5)
    real(dp)               :: f(nw_euler, 2)
    real(dp)               :: inv_dr(ndim)
    integer                :: typelim
    !-----------------------------------------------------------------------------

 !$acc parallel loop private(n, uprim, inv_dr, typelim) firstprivate(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) present(bga%w, bgb%w)
    do iigrid = 1, igridstail_active
       n = igrids_active(iigrid)

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

       !$acc loop collapse(ndim) private(f, tmp) vector
        do ix3=ixOmin3,ixOmax3 
         do ix2=ixOmin2,ixOmax2 
         do ix1=ixOmin1,ixOmax1 
             ! Compute fluxes in all dimensions

       
       
       
             tmp = uprim(:, ix1-2:ix1+2, ix2, ix3)
             call muscl_flux_euler_prim(tmp, 1, f, typelim)
             bgb%w(ix1, ix2, ix3, :, n) = bgb%w(ix1, ix2, ix3, :,&
                 n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(1)

             tmp = uprim(:, ix1, ix2-2:ix2+2, ix3)
             call muscl_flux_euler_prim(tmp, 2, f, typelim)
             bgb%w(ix1, ix2, ix3, :, n) = bgb%w(ix1, ix2, ix3, :,&
                 n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(2)

             tmp = uprim(:, ix1, ix2, ix3-2:ix3+2)
             call muscl_flux_euler_prim(tmp, 3, f, typelim)
             bgb%w(ix1, ix2, ix3, :, n) = bgb%w(ix1, ix2, ix3, :,&
                 n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(3)
      

        end do 
         end do 
         end do 
    end do

  end subroutine finite_volume_local

  pure subroutine to_primitive(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(nw_euler)

    
    u(iw_mom(1)) = u(iw_mom(1))/u(iw_rho)
    
    
    u(iw_mom(2)) = u(iw_mom(2))/u(iw_rho)
    
    
    u(iw_mom(3)) = u(iw_mom(3))/u(iw_rho)
    

    u(iw_e) = (hd_gamma-1.0_dp) * (u(iw_e) - 0.5_dp * u(iw_rho) * &
       sum(u(iw_mom(1:ndim))**2) )

  end subroutine to_primitive

  pure subroutine to_conservative(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(nw_euler)
    real(dp)                :: inv_gamma_m1

    inv_gamma_m1 = 1.0d0/(hd_gamma - 1.0_dp)

    ! Compute energy from pressure and kinetic energy
    u(iw_e) = u(iw_e) * inv_gamma_m1 + 0.5_dp * u(iw_rho) * &
       sum(u(iw_mom(1:ndim))**2)

    ! Compute momentum from density and velocity components
    
    u(iw_mom(1)) = u(iw_rho) * u(iw_mom(1))
    
    
    u(iw_mom(2)) = u(iw_rho) * u(iw_mom(2))
    
    
    u(iw_mom(3)) = u(iw_rho) * u(iw_mom(3))
    

  end subroutine to_conservative

  subroutine muscl_flux_euler_prim(u, flux_dim, flux, typelim)
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

    call euler_flux(uL, flux_dim, flux_l)
    call euler_flux(uR, flux_dim, flux_r)

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

    call euler_flux(uL, flux_dim, flux_l)
    call euler_flux(uR, flux_dim, flux_r)

    wL = get_cmax(uL, flux_dim)
    wR = get_cmax(uR, flux_dim)
    wmax = max(wL, wR)

    call to_conservative(uL)
    call to_conservative(uR)
    flux(:, 2) = 0.5_dp * ((flux_l + flux_r) - wmax * (uR - uL))

  end subroutine muscl_flux_euler_prim

  subroutine euler_flux(u, flux_dim, flux)
    !$acc routine seq
    real(dp), intent(in)  :: u(nw_euler)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(nw_euler)
    real(dp)              :: inv_gamma_m1

    inv_gamma_m1 = 1.0d0/(hd_gamma - 1.0_dp)

    ! Density flux
    flux(iw_rho) = u(iw_rho) * u(iw_mom(flux_dim))

    ! Momentum flux with pressure term
    
    flux(iw_mom(1)) = u(iw_rho) * u(iw_mom(1)) * u(iw_mom(flux_dim))
    
    
    flux(iw_mom(2)) = u(iw_rho) * u(iw_mom(2)) * u(iw_mom(flux_dim))
    
    
    flux(iw_mom(3)) = u(iw_rho) * u(iw_mom(3)) * u(iw_mom(flux_dim))
    
    flux(iw_mom(flux_dim)) = flux(iw_mom(flux_dim)) + u(iw_e)

    ! Energy flux
    flux(iw_e) = u(iw_mom(flux_dim)) * (u(iw_e) * inv_gamma_m1 + 0.5_dp * &
       u(iw_rho) * sum(u(iw_mom(1:ndim))**2) + u(iw_e))

  end subroutine euler_flux

  pure real(dp) function get_cmax(u, flux_dim) result(wC)
    !$acc routine seq
    real(dp), intent(in)  :: u(nw_euler)
    integer, intent(in)   :: flux_dim

    wC = sqrt(hd_gamma * u(iw_e) / u(iw_rho)) + abs(u(iw_mom(flux_dim)))

  end function get_cmax

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
