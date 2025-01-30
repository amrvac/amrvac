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
  use mod_physicaldata
  implicit none

  private

  public :: finite_volume_all
  public :: finite_volume_local
  public :: reconstruct_LR_gpu

  integer, parameter :: dp = kind(0.0d0), nw_euler=4

contains

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
    integer                :: n, i, j
    double precision       :: w(1:nw_euler)
    double precision       :: uprim(ixI^S,1:nw_euler)
    real(dp)               :: tmp(5, nw_euler), u(nw_euler)
    real(dp)               :: fx(nw_euler, 2), fy(nw_euler, 2)
    real(dp)               :: inv_dr(2)
    !-----------------------------------------------------------------------------

    !$acc parallel loop, private(w, uprim, inv_dr)
    do n = 1, igridstail_active

       inv_dr = 1/rnode(rpdx1_:rnodehi, n)

       !$acc loop collapse(2), private(w), vector
       do j = ixImin2, ixImax2
          do i = ixImin1, ixImax1
             ! Convert to primitive
             w = bga%w(i, j, :, n)
             call to_primitive(w)
             uprim(i, j, :) = w
          end do
       end do

       !$acc loop collapse(2), private(fx, fy, tmp), vector
       do j = ixOmin2, ixOmax2
          do i = ixOmin1, ixOmax1
             ! Compute x and y fluxes
             tmp = uprim(i-2:i+2, j, :)
             call muscl_flux_euler_prim(tmp, 1, fx)

             tmp = uprim(i, j-2:j+2, :)
             call muscl_flux_euler_prim(tmp, 2, fy)

             ! Update the wnew array
             bgb%w(i, j, :, n) = bgb%w(i, j, :, n) + qdt * &
                  ((fx(:, 1) - fx(:, 2)) * inv_dr(1) + &
                  (fy(:, 1) - fy(:, 2)) * inv_dr(2))
          end do
       end do
    end do

  end subroutine finite_volume_local

  pure subroutine to_primitive(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(nw_euler)

    u(iw_mom(1)) = u(iw_mom(1))/u(iw_rho)
    u(iw_mom(2)) = u(iw_mom(2))/u(iw_rho)

    u(iw_e) = (hd_gamma-1.0_dp) * (u(iw_e) - &
         0.5_dp * u(iw_rho)* (u(iw_mom(1))**2 + u(iw_mom(2))**2))

  end subroutine to_primitive

  pure subroutine to_conservative(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(nw_euler)
    real(dp)                :: inv_gamma_m1

    inv_gamma_m1 = 1.0d0/(hd_gamma - 1.0_dp)

    ! Compute energy from pressure and kinetic energy
    u(iw_e) = u(iw_e) * inv_gamma_m1 + &
         0.5_dp * u(iw_rho) * (u(iw_mom(1))**2 + u(iw_mom(2))**2)

    ! Compute momentum from density and velocity components
    u(iw_mom(1)) = u(iw_rho) * u(iw_mom(1))
    u(iw_mom(2)) = u(iw_rho) * u(iw_mom(2))
  end subroutine to_conservative

  subroutine muscl_flux_euler_prim(u, flux_dim, flux)
    !$acc routine seq
    real(dp), intent(in)  :: u(5, nw_euler)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(nw_euler, 2)
    real(dp)              :: uL(nw_euler), uR(nw_euler), wL, wR, wmax
    real(dp)              :: flux_l(nw_euler), flux_r(nw_euler)

    ! Construct uL, uR for first cell face
    uL = u(2, :) + 0.5_dp * vanleer(u(2, :) - u(1, :), u(3, :) - u(2, :))
    uR = u(3, :) - 0.5_dp * vanleer(u(3, :) - u(2, :), u(4, :) - u(3, :))

    call euler_flux(uL, flux_dim, flux_l)
    call euler_flux(uR, flux_dim, flux_r)

    wL = get_cmax(uL, flux_dim)
    wR = get_cmax(uR, flux_dim)
    wmax = max(wL, wR)

    call to_conservative(uL)
    call to_conservative(uR)
    flux(:, 1) = 0.5_dp * ((flux_l + flux_r) - wmax * (uR - uL))

    ! Construct uL, uR for second cell face
    uL = u(3, :) + 0.5_dp * vanleer(u(3, :) - u(2, :), u(4, :) - u(3, :))
    uR = u(4, :) - 0.5_dp * vanleer(u(4, :) - u(3, :), u(5, :) - u(4, :))

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
    flux(iw_mom(flux_dim)) = flux(iw_mom(flux_dim)) + u(iw_e)

    ! Energy flux
    flux(iw_e) = u(iw_mom(flux_dim)) * (u(iw_e) * inv_gamma_m1 + &
         0.5_dp * u(iw_rho) * (u(iw_mom(1))**2 + u(iw_mom(2))**2) + u(iw_e))

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


  !> The non-conservative Hancock predictor for TVDLF
  !>
  !> on entry:
  !> input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
  !> one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
  !> on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
  !AGILE
  !> finite volume method computing all blocks
  subroutine finite_volume_all(method, qdt, dtfactor, ixI^L, ixO^L, idims^LIM, &
       qtC, bga, qt, bgb, fC, fE)
    use mod_physics
    use mod_variables
    use mod_global_parameters
    use mod_source, only: addsource2
    use mod_usr_methods
    use mod_comm_lib, only: mpistop
    use mod_hd_phys

    integer, intent(in)                                   :: method
    double precision, intent(in)                          :: qdt, dtfactor, qtC, qt
    integer, intent(in)                                   :: ixI^L, ixO^L, idims^LIM
    type(block_grid_t)                                    :: bga
    type(block_grid_t)                                    :: bgb
    double precision, dimension(ixI^S, 1:nwflux, 1:ndim)  :: fC
    double precision, dimension(ixI^S, sdim:3)            :: fE

    ! primitive w at cell center
    double precision, dimension(ixI^S, 1:nw)       :: wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S, 1:nw)       :: wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S, 1:nw)       :: wLp, wRp
    !$acc declare create(wprim, wLC, wRC, wLp, wRp)
    double precision, dimension(ixI^S, 1:nwflux)   :: fLC, fRC
    !$acc declare create(fLC, fRC)
    double precision, dimension(ixI^S, 1:number_species)      :: cmaxC
    double precision, dimension(ixI^S, 1:number_species)      :: cminC
    double precision, dimension(ixI^S, 1:number_species)      :: Hspeed
    !$acc declare create(cmaxC, cminC, Hspeed)
    double precision, dimension(1:ndim)        :: dxinv
    integer                                    :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L, ii
    integer                                    :: ix^D
    integer                                    :: igrid, iigrid, ia, ib
    double precision, dimension(ixI^S, 1:ndim) :: x
    double precision                           :: dxs(ndim)
    !$acc declare create(dxs,dxinv)

    ia = bga%istep; ib = bgb%istep ! remember, old names map as: wCT => bga, wnew => bgb

    !$acc parallel loop private(fC, fLC, fRC, wprim, x, wRp, wLp, wLC, wRC, cmaxC, cminC, Hspeed, dxinv, dxs)
    do iigrid = 1, igridstail_active
       igrid = igrids_active(iigrid)

       x = ps(igrid)%x
       dxs = rnode(rpdx1_:rnodehi, igrid)

       ! The flux calculation contracts by one in the idims direction it is applied.
       ! The limiter contracts the same directions by one more, so expand ixO by 2.
       ix^L=ixO^L; 
       do idims = idims^LIM
          ix^L=ix^L^LADD2*kr(idims,^D); 
       end do
       if (ixI^L^LTix^L|.or.|.or.) &
            call mpistop("Error in fv : Nonconforming input limits")

       fC = 0.d0     
       fLC = 0.d0
       fRC = 0.d0

       wprim = bg(ia)%w(:^D&, :, igrid)

       call hd_to_primitive_gpu(ixI^L, ixI^L, wprim, x)

       do idims = idims^LIM
          ! use interface value of w0 at idims

          hxO^L=ixO^L-kr(idims,^D); 
          kxCmin^D=ixImin^D;kxCmax^D=ixImax^D-kr(idims,^D); 
          kxR^L=kxC^L+kr(idims,^D); 
          ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
          ixCmax^D=ixOmax^D;ixCmin^D=hxOmin^D; 

          ! wRp and wLp are defined at the same locations, and will correspond to
          ! the left and right reconstructed values at a cell face. Their indexing
          ! is similar to cell-centered values, but in direction idims they are
          ! shifted half a cell towards the 'lower' direction.

          wRp(kxC^S,1:nw)=wprim(kxR^S,1:nw)
          wLp(kxC^S,1:nw)=wprim(kxC^S,1:nw)


          ! Determine stencil size
          {ixCRmin^D = max(ixCmin^D - phys_wider_stencil,ixGlo^D) \}
          {ixCRmax^D = min(ixCmax^D + phys_wider_stencil,ixGhi^D) \}

          ! apply limited reconstruction for left and right status at cell interfaces
          call reconstruct_LR_gpu(ixI^L, ixCR^L, ixCR^L, idims, wprim, wLC, wRC, wLp, wRp, x, dxs(idims), igrid)

          ! evaluate physical fluxes according to reconstructed status
          call hd_get_flux_gpu(wLC, wLp, x, ixI^L, ixC^L, idims, fLC)
          call hd_get_flux_gpu(wRC, wRp, x, ixI^L, ixC^L, idims, fRC)

          ! get the min and max characteristic velocities
          call hd_get_cbounds_gpu(wLC, wRC, wLp, wRp, x, ixI^L, ixC^L, idims, Hspeed, cmaxC, cminC)

          ! use approximate Riemann solver to get flux at interfaces
          do ii = 1, number_species
             do iw = start_indices(ii), stop_indices(ii)
                {do ix^DB = ixCmin^DB, ixCmax^DB\}
                if (cminC(ix^D, ii) >= zero) then
                   fC(ix^D, iw, idims) = fLC(ix^D, iw)
                else if (cmaxC(ix^D, ii) <= zero) then
                   fC(ix^D, iw, idims) = fRC(ix^D, iw)
                else
                   ! Add hll dissipation to the flux
                   fC(ix^D, iw, idims) = (cmaxC(ix^D, ii)*fLC(ix^D, iw) - cminC(ix^D, ii)*fRC(ix^D, iw) &
                        + cminC(ix^D, ii)*cmaxC(ix^D, ii)*(wRC(ix^D, iw) - wLC(ix^D, iw))) &
                        /(cmaxC(ix^D, ii) - cminC(ix^D, ii))
                end if
                {end do\}
             end do
          end do


       end do ! Next idims

       dxinv = -qdt/dxs
       do idims = idims^LIM
          hxO^L=ixO^L-kr(idims,^D); 
          fC(ixI^S, 1:nwflux, idims)=dxinv(idims)*fC(ixI^S, 1:nwflux, idims)

          bg(ib)%w(ixO^S, iwstart:nwflux, igrid) = bg(ib)%w(ixO^S, iwstart:nwflux, igrid) + &
               (fC(ixO^S, iwstart:nwflux, idims) - fC(hxO^S, iwstart:nwflux, idims))

       end do ! Next idims

    end do   ! igrid

  end subroutine finite_volume_all

  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR_gpu(ixI^L, ixL^L, ixR^L, idims, w, wLC, wRC, wLp, wRp, x, dxdim, igrid)
    use mod_physics
    use mod_global_parameters
    use mod_limiter
    use mod_comm_lib, only: mpistop
    use mod_hd_phys, only: hd_to_conserved_gpu

    !$acc routine

    integer, value, intent(in) :: ixI^L, ixL^L, ixR^L, idims, igrid
    double precision, intent(in) :: dxdim
    ! cell center w in primitive form
    double precision, dimension(ixI^S, 1:nw) :: w
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S, 1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form
    double precision, dimension(ixI^S, 1:nw) :: wLp, wRp
    double precision, dimension(ixI^S, 1:ndim) :: x

    integer            :: jxR^L, ixC^L, jxC^L, iw
    double precision   :: ldw(ixI^S), rdw(ixI^S), dwC(ixI^S)
!!!!$acc declare create(ldw, rdw, dwC)
    double precision   :: a2max

    select case (type_limiter(node(plevel_, igrid)))
    case default
       jxR^L=ixR^L+kr(idims,^D); 
       ixCmax^D=jxRmax^D;ixCmin^D=ixLmin^D-kr(idims,^D); 
       jxC^L=ixC^L+kr(idims,^D); 
       do iw = 1, nwflux
          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D, iw) = dlog10(w(ixCmin^D:jxCmax^D, iw))
             wLp(ixL^S, iw) = dlog10(wLp(ixL^S, iw))
             wRp(ixR^S, iw) = dlog10(wRp(ixR^S, iw))
          end if

          dwC(ixC^S) = w(jxC^S, iw) - w(ixC^S, iw)
          if (need_global_a2max) then
             a2max = a2max_global(idims)
          else
             select case (idims)
             case (1)
                a2max = schmid_rad1
                {^IFTWOD
             case (2)
                a2max = schmid_rad2}
                {^IFTHREED
             case (2)
                a2max = schmid_rad2
             case (3)
                a2max = schmid_rad3}
             case default
                call mpistop("idims is wrong in mod_limiter")
             end select
          end if

          ! limit flux from left and/or right
          call dwlimiter2_gpu(dwC, ixI^L, ixC^L, idims, type_limiter(node(plevel_, igrid)), ldw, rdw, a2max=a2max)
          wLp(ixL^S, iw) = wLp(ixL^S, iw) + half*ldw(ixL^S)
          wRp(ixR^S, iw) = wRp(ixR^S, iw) - half*rdw(jxR^S)

          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D, iw) = 10.0d0**w(ixCmin^D:jxCmax^D, iw)
             wLp(ixL^S, iw) = 10.0d0**wLp(ixL^S, iw)
             wRp(ixR^S, iw) = 10.0d0**wRp(ixR^S, iw)
          end if
       end do
       ! if(fix_small_values) then
       !    call phys_handle_small_values(.true.,wLp,x,ixI^L,ixL^L,'reconstruct left')
       !    call phys_handle_small_values(.true.,wRp,x,ixI^L,ixR^L,'reconstruct right')
       ! end if

    end select

    wLC(ixL^S, 1:nwflux) = wLp(ixL^S, 1:nwflux)
    wRC(ixR^S, 1:nwflux) = wRp(ixR^S, 1:nwflux)

    call hd_to_conserved_gpu(ixI^L, ixL^L, wLC, x)
    call hd_to_conserved_gpu(ixI^L, ixR^L, wRC, x)

    if (nwaux > 0) then
       wLp(ixL^S, nwflux + 1:nwflux + nwaux) = wLC(ixL^S, nwflux + 1:nwflux + nwaux)
       wRp(ixR^S, nwflux + 1:nwflux + nwaux) = wRC(ixR^S, nwflux + 1:nwflux + nwaux)
    end if

  end subroutine reconstruct_LR_gpu

end module mod_finite_volume_all
