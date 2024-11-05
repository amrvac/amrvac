! - [x] Remove everything we don't use
! - [ ] Remove psa, and psb from finite_volume (ps, bga and bgb should be enough)
! - [ ] Test for correctness
! - [ ] Replace block with ps(igrid), pass igrid to functions that need it
!> Module with finite volume methods for fluxes
module mod_finite_volume_all
#include "amrvac.h"
   implicit none
   private

   public :: finite_volume_all

contains

   !> The non-conservative Hancock predictor for TVDLF
   !>
   !> on entry:
   !> input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
   !> one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
   !> on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
   !AGILE
   !> finite volume method computing all blocks
   subroutine finite_volume_all(method, qdt, dtfactor, ixI^L, ixO^L, idims^LIM, &
                                qtC, psa, bga, qt, psb, bgb, fC, fE)
      use mod_physics
      use mod_variables
      use mod_global_parameters
      !opedit: FIXME has internal procedure pointers, used for roe solvers, not implemented:
!    use mod_tvd, only:tvdlimit2
      use mod_source, only: addsource2
      use mod_usr_methods
      use mod_comm_lib, only: mpistop
      use mod_hd_phys, only: hd_get_flux_gpu, hd_to_primitive_gpu, hd_get_cbounds_gpu, hd_to_conserved_gpu

      integer, intent(in)                                   :: method
      double precision, intent(in)                          :: qdt, dtfactor, qtC, qt
      integer, intent(in)                                   :: ixI^L, ixO^L, idims^LIM
      type(state), target                                   :: psa(max_blocks)
      type(state), target                                   :: psb(max_blocks)
      type(block_grid_t)                            :: bga
      type(block_grid_t)                            :: bgb
      double precision, dimension(ixI^S, 1:nwflux, 1:ndim)    :: fC
      double precision, dimension(ixI^S, sdim:3)             :: fE

      ! primitive w at cell center
      double precision, dimension(ixI^S, 1:nw) :: wprim
      ! left and right constructed status in conservative form
      double precision, dimension(ixI^S, 1:nw) :: wLC, wRC
      ! left and right constructed status in primitive form, needed for better performance
      double precision, dimension(ixI^S, 1:nw) :: wLp, wRp
      !$acc declare create(wprim, wLC, wRC, wLp, wRp)
      double precision, dimension(ixI^S, 1:nwflux) :: fLC, fRC
      !$acc declare create(fLC, FRC)
      double precision, dimension(ixI^S, 1:number_species)      :: cmaxC
      double precision, dimension(ixI^S, 1:number_species)      :: cminC
      double precision, dimension(ixI^S)      :: Hspeed
      !$acc declare create(cmaxC, cminC, Hspeed)
      double precision, dimension(ixO^S)      :: inv_volume
      double precision, dimension(1:ndim)     :: dxinv
      integer, dimension(ixI^S)               :: patchf
      integer :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L, ii
      logical :: active
      type(ct_velocity) :: vcts
      integer :: ix^D
      double precision, dimension(ixI^S, 1:nwflux)     :: whll, Fhll, fCD
      double precision, dimension(ixI^S)              :: lambdaCD
      integer  :: rho_, p_, e_, mom(1:ndir), igrid, iigrid

      double precision, dimension(ixI^S, 1:ndim) :: x
      double precision                          :: dxs(ndim)

      !$acc parallel loop private(block)
      do iigrid = 1, igridstail_active
         igrid = igrids_active(iigrid)
         block => ps(igrid)
         !$acc enter data attach(block)

         x = ps(igrid)%x
         dxs = rnode(rpdx1_:rnodehi, igrid)

         ! TODO wCT and therefor wprim references need same treatment as sCT%w
         associate (wCT => bga%w(:^D&, :, igrid), wnew => bgb%w(:^D&, :, igrid))
            ! The flux calculation contracts by one in the idims direction it is applied.
            ! The limiter contracts the same directions by one more, so expand ixO by 2.
            ix^L = ixO^L; 
            do idims = idims^LIM
               ix^L = ix^L^LADD2*kr(idims, ^D); 
            end do
            if (ixI^L^LTix^L| .or. | .or. ) &
               call mpistop("Error in fv : Nonconforming input limits")

            ! no longer !$acc kernels present(fC, wCT)
            fC = 0.d0     ! this is updated every loop iteration (eg l293),
            ! how to then parallelize?
            fLC = 0.d0
            fRC = 0.d0

            wprim = wCT
            ! no longer !$acc end kernels

            call hd_to_primitive_gpu(ixI^L, ixI^L, wprim, x)

            do idims = idims^LIM
               ! use interface value of w0 at idims
               b0i = idims

               hxO^L = ixO^L - kr(idims, ^D); 
               kxCmin^D = ixImin^D; kxCmax^D = ixImax^D - kr(idims, ^D); 
               kxR^L = kxC^L + kr(idims, ^D); 
               if (stagger_grid) then
                  ! ct needs all transverse cells
                  ixCmax^D = ixOmax^D + nghostcells - nghostcells*kr(idims, ^D); 
                  ixCmin^D = hxOmin^D - nghostcells + nghostcells*kr(idims, ^D); 
               else
                  ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
                  ixCmax^D = ixOmax^D; ixCmin^D = hxOmin^D; 
               end if

               ! wRp and wLp are defined at the same locations, and will correspond to
               ! the left and right reconstructed values at a cell face. Their indexing
               ! is similar to cell-centered values, but in direction idims they are
               ! shifted half a cell towards the 'lower' direction.

               wRp(kxC^S, 1:nw) = wprim(kxR^S, 1:nw)
               wLp(kxC^S, 1:nw) = wprim(kxC^S, 1:nw)

               ! Determine stencil size
               ! FIXME: here `phys_wider_stencil` is an integer, not a function pointer
               {ixCRmin^D = max(ixCmin^D - phys_wider_stencil, ixGlo^D) \}
               {ixCRmax^D = min(ixCmax^D + phys_wider_stencil, ixGhi^D) \}

               ! apply limited reconstruction for left and right status at cell interfaces
               call reconstruct_LR_gpu(ixI^L, ixCR^L, ixCR^L, idims, wprim, wLC, wRC, wLp, wRp, x, dxs(idims), igrid)

               ! evaluate physical fluxes according to reconstructed status
               call hd_get_flux_gpu(wLC, wLp, x, ixI^L, ixC^L, idims, fLC)
               call hd_get_flux_gpu(wRC, wRp, x, ixI^L, ixC^L, idims, fRC)

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

            b0i = 0
            dxinv = -qdt/dxs
            do idims = idims^LIM
               hxO^L = ixO^L - kr(idims, ^D); 
               fC(ixI^S, 1:nwflux, idims) = dxinv(idims)*fC(ixI^S, 1:nwflux, idims)
               !            end if
               wnew(ixO^S, iwstart:nwflux) = wnew(ixO^S, iwstart:nwflux) + &
                                             (fC(ixO^S, iwstart:nwflux, idims) - fC(hxO^S, iwstart:nwflux, idims))

            end do ! Next idims

            call addsource2(qdt*dble(idimsmax - idimsmin + 1)/dble(ndim), &
                            dtfactor*dble(idimsmax - idimsmin + 1)/dble(ndim), &
                            ixI^L, ixO^L, 1, nw, qtC, wCT, wprim, qt, wnew, x, .false., active)
         end associate
         !$acc exit data detach(block)
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

      !$acc declare present(node)
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
         jxR^L = ixR^L + kr(idims, ^D); 
         ixCmax^D = jxRmax^D; ixCmin^D = ixLmin^D - kr(idims, ^D); 
         jxC^L = ixC^L + kr(idims, ^D); 
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
