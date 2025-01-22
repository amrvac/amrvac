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
   public :: reconstruct_LR_gpu

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
      use mod_hd_phys!, only: hd_get_flux_gpu, hd_to_primitive_gpu, hd_get_cbounds_gpu, hd_to_conserved_gpu

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
      !$acc declare create(fLC, fRC)
      double precision, dimension(ixI^S, 1:number_species)      :: cmaxC
      double precision, dimension(ixI^S, 1:number_species)      :: cminC
      double precision, dimension(ixI^S, 1:number_species)      :: Hspeed
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
      integer  :: igrid, iigrid

      double precision, dimension(ixI^S, 1:ndim) :: x
      double precision                          :: dxs(ndim)
      !$acc declare create(dxs,dxinv)
      !opedit debug:
      integer         :: idbg^D

      !$acc parallel loop private(fC, fLC, fRC, wprim, x, wRp, wLp, wLC, wRC, cmaxC, cminC, Hspeed, dxinv, dxs)
      do iigrid = 1, igridstail_active
         igrid = igrids_active(iigrid)

         x = ps(igrid)%x
         dxs = rnode(rpdx1_:rnodehi, igrid)

         associate (wCT => bga%w(:^D&, :, igrid), wnew => bgb%w(:^D&, :, igrid))

           print *, 'fvolume_all A:', igrid, maxval(wnew(:,:,1)), maxval(wCT(:,:,1))
            ! if (any(wnew(:,:,1) .ne. wnew(:,:,1)) ) then
            !    print *, 'NaN found in wnew A'
            !    do idbg2=ixImin2,ixImax2
            !       do idbg1=ixImin1,ixImax1
            !          print *, idbg1, idbg2, wnew(idbg1,idbg2,1)
            !       end do
            !    end do
            ! end if
            if (any(wnew(:,:,1) .ne. wnew(:,:,1)) ) then
            print *, 'NaN found in wnew A'
            print *, igrid, 'wnew A',ixImin1,ixImax1,ixImin2,ixImax2
               do idbg2=ixImin2,ixImax2
                  do idbg1=ixImin1,ixImax1
                     print *, igrid, idbg1, idbg2, wnew(idbg1,idbg2,1)
                  end do
               end do
            end if
           
            ! The flux calculation contracts by one in the idims direction it is applied.
            ! The limiter contracts the same directions by one more, so expand ixO by 2.
            ix^L=ixO^L; 
            do idims = idims^LIM
               ix^L=ix^L^LADD2*kr(idims,^D); 
            end do
            if (ixI^L^LTix^L|.or.|.or.) &
               call mpistop("Error in fv : Nonconforming input limits")

            ! no longer !$acc kernels present(fC, wCT)
            fC = 0.d0     ! this is updated every loop iteration (eg l293),
            ! how to then parallelize?
            fLC = 0.d0
            fRC = 0.d0

            wprim = wCT
            ! no longer !$acc end kernels

            call hd_to_primitive_gpu(ixI^L, ixI^L, wprim, x)
            print *, 'hd_gamma:', igrid, hd_gamma
            
            if (any(wprim(ixI^S,p_) .le. 0.0d0) ) then
               print *, igrid, 'negative pressure found A-2'
               do idbg2=ixImin2,ixImax2
                  do idbg1=ixImin1,ixImax1
                     if (wprim(idbg1,idbg2,p_) .lt. 0) then
                        write(*,*), 'wprim p', igrid, idbg1, idbg2, wprim(idbg1,idbg2,p_)
                     end if
                  end do
               end do
            end if

            if (any(wprim(ixI^S,1) .ne. wprim(ixI^S,1)) ) then
               print *, igrid, 'NaN found in wprim'
               do idbg2=ixImin2,ixImax2
                  do idbg1=ixImin1,ixImax1
                     if (wprim(idbg1,idbg2,1) .ne. wprim(idbg1,idbg2,1)) then
                        write(*,*), 'wp', igrid, idbg1, idbg2, wprim(idbg1,idbg2,1)
                     end if
                  end do
               end do
            end if

            do idims = idims^LIM
               ! use interface value of w0 at idims
               b0i = idims

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

               
               if (any(wRp(kxC^S,p_) .le. 0.0d0) .or. any(wLp(kxC^S,p_) .le. 0.0d0)) then
                  print *, igrid, 'negative pressure found A-1'
                  do idbg2=kxCmin2,kxCmax2
                     do idbg1=kxCmin1,kxCmax1
                        if (wRp(idbg1,idbg2,p_) .le. 0.0d0 .or. wLp(idbg1,idbg2,p_) .le. 0.0d0) then
                           write(*,*), 'pR', igrid, idbg1, idbg2, wRp(idbg1,idbg2,p_)
                           write(*,*), 'pL', igrid, idbg1, idbg2, wLp(idbg1,idbg2,p_)
                        end if
                     end do
                  end do
               end if

               
               ! Determine stencil size
               ! FIXME: here `phys_wider_stencil` is an integer, not a function pointer
               {ixCRmin^D = max(ixCmin^D - phys_wider_stencil,ixGlo^D) \}
               {ixCRmax^D = min(ixCmax^D + phys_wider_stencil,ixGhi^D) \}

               ! apply limited reconstruction for left and right status at cell interfaces
               call reconstruct_LR_gpu(ixI^L, ixCR^L, ixCR^L, idims, wprim, wLC, wRC, wLp, wRp, x, dxs(idims), igrid)
               
               if (any(wLC(ixCR^S,1) .ne. wLC(ixCR^S,1)) .or. any(wLp(ixCR^S,1) .ne. wLp(ixCR^S,1)) ) then
                  print *, igrid, 'NaN found in wLC'
                  do idbg2=ixCRmin2,ixCRmax2
                     do idbg1=ixCRmin1,ixCRmax1
                        if (wLC(idbg1,idbg2,1) .ne. wLC(idbg1,idbg2,1) .or. wLp(idbg1,idbg2,1) .ne. wLp(idbg1,idbg2,1) ) then
                           write(*,*), 'wLC', igrid, idbg1, idbg2, wLC(idbg1,idbg2,1)
                           write(*,*), 'wLp', igrid, idbg1, idbg2, wLp(idbg1,idbg2,1)
                        end if
                     end do
                  end do
               end if
               
               if (any(wRC(ixCR^S,1) .ne. wRC(ixCR^S,1)) .or. any(wRp(ixCR^S,1) .ne. wRp(ixCR^S,1))) then
                  print *, igrid, 'NaN found in wRC'
                  do idbg2=ixCRmin2,ixCRmax2
                     do idbg1=ixCRmin1,ixCRmax1
                        if (wRC(idbg1,idbg2,1) .ne. wRC(idbg1,idbg2,1) .or. wRp(idbg1,idbg2,1) .ne. wRp(idbg1,idbg2,1) ) then
                           write(*,*), 'wRC', igrid, idbg1, idbg2, wRC(idbg1,idbg2,1)
                           write(*,*), 'wRp', igrid, idbg1, idbg2, wRp(idbg1,idbg2,1)
                        end if
                     end do
                  end do
               end if

               
               if (any(wRp(ixCR^S,p_) .le. 0.0d0) .or. any(wLp(ixCR^S,p_) .le. 0.0d0)) then
                  print *, igrid, 'negative pressure found A'
                  do idbg2=ixCRmin2,ixCRmax2
                     do idbg1=ixCRmin1,ixCRmax1
                        if (wRp(idbg1,idbg2,p_) .le. 0.0d0 .or. wLp(idbg1,idbg2,p_) .le. 0.0d0) then
                           write(*,*), 'pR', igrid, idbg1, idbg2, wRp(idbg1,idbg2,p_)
                           write(*,*), 'pL', igrid, idbg1, idbg2, wLp(idbg1,idbg2,p_)
                        end if
                     end do
                  end do
               end if

               ! evaluate physical fluxes according to reconstructed status
               call hd_get_flux_gpu(wLC, wLp, x, ixI^L, ixC^L, idims, fLC)

               if (any(fLC(ixC^S,1) .ne. fLC(ixC^S,1)) ) then
                  print *, 'NaN found in fLC', idims
                  do idbg2=ixCmin2,ixCmax2
                     do idbg1=ixCmin1,ixCmax1
                        if (fLC(idbg1,idbg2,1) .ne. fLC(idbg1,idbg2,1)) then
                           write(*,*), 'FLC', igrid, idbg1, idbg2, fLC(idbg1,idbg2,1)
                        end if
                     end do
                  end do
               end if

               call hd_get_flux_gpu(wRC, wRp, x, ixI^L, ixC^L, idims, fRC)
               
               if (any(fRC(ixC^S,1) .ne. fRC(ixC^S,1)) ) then
                  print *, 'NaN found in fRC', idims
                  do idbg2=ixCmin2,ixCmax2
                     do idbg1=ixCmin1,ixCmax1
                        if (fRC(idbg1,idbg2,1) .ne. fRC(idbg1,idbg2,1)) then
                           write(*,*), 'FRC', igrid, idbg1, idbg2, fRC(idbg1,idbg2,1)
                        end if
                     end do
                  end do
               end if
               
               if (any(wRp(ixCR^S,p_) .le. 0.0d0) .or. any(wLp(ixCR^S,p_) .le. 0.0d0)) then
                  print *, igrid, 'negative pressure found B'
                  do idbg2=ixCRmin2,ixCRmax2
                     do idbg1=ixCRmin1,ixCRmax1
                        if (wRp(idbg1,idbg2,p_) .le. 0.0d0 .or. wLp(idbg1,idbg2,p_) .le. 0.0d0) then
                           write(*,*), 'pR', igrid, idbg1, idbg2, wRp(idbg1,idbg2,p_)
                           write(*,*), 'pL', igrid, idbg1, idbg2, wLp(idbg1,idbg2,p_)
                        end if
                     end do
                  end do
               end if
               
               call hd_get_cbounds_gpu(wLC, wRC, wLp, wRp, x, ixI^L, ixC^L, idims, Hspeed, cmaxC, cminC)

               if (any(cmaxC(ixC^S,1) .ne. cmaxC(ixC^S,1)) ) then
                  print *, 'NaN found in cmaxC', idims
                  do idbg2=ixCmin2,ixCmax2
                     do idbg1=ixCmin1,ixCmax1
                        if (cmaxC(idbg1,idbg2,1) .ne. cmaxC(idbg1,idbg2,1)) then
                           write(*,*), 'cmax', igrid, idbg1, idbg2, cmaxC(idbg1,idbg2,1)
                        end if
                     end do
                  end do
               end if
               
               if (any(cminC(ixC^S,1) .ne. cminC(ixC^S,1)) ) then
                  print *, 'NaN found in cminC', idims
                  do idbg2=ixCmin2,ixCmax2
                     do idbg1=ixCmin1,ixCmax1
                        if (cminC(idbg1,idbg2,1) .ne. cminC(idbg1,idbg2,1)) then
                           write(*,*), 'cmin', igrid, idbg1, idbg2, cminC(idbg1,idbg2,1)
                        end if
                     end do
                  end do
               end if
               
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

               
            if (any(fC(ixC^S,1,idims) .ne. fc(ixC^S,1,idims)) ) then
               print *, 'NaN found in fC A', idims
               do idbg2=ixCmin2,ixCmax2
                  do idbg1=ixCmin1,ixCmax1
                     if (fC(idbg1,idbg2,1,idims) .ne. fC(idbg1,idbg2,1,idims)) then
                        write(*,*), 'F', igrid, idbg1, idbg2, fC(idbg1,idbg2,1,idims)
                        write(*,*), 'A', igrid, idbg1, idbg2, cmaxC(idbg1,idbg2,1)
                        write(*,*), 'I', igrid, idbg1, idbg2, cminC(idbg1,idbg2,1)
                     end if
                  end do
               end do
            end if
               
            end do ! Next idims

            b0i = 0
            dxinv = -qdt/dxs
            do idims = idims^LIM
               hxO^L=ixO^L-kr(idims,^D); 
               fC(ixI^S, 1:nwflux, idims)=dxinv(idims)*fC(ixI^S, 1:nwflux, idims)

               
            ! if (any(fC(ixI^S,1,idims) .ne. fc(ixI^S,1,idims)) ) then
            !    print *, 'NaN found in fC B', idims
            !    do idbg2=ixImin2,ixImax2
            !       do idbg1=ixImin1,ixImax1
            !          print *, igrid, idbg1, idbg2, fC(idbg1,idbg2,1,idims)
            !       end do
            !    end do
            ! end if
               
               !            end if
               wnew(ixO^S, iwstart:nwflux)=wnew(ixO^S, iwstart:nwflux) + &
                                             (fC(ixO^S, iwstart:nwflux, idims) - fC(hxO^S, iwstart:nwflux, idims))

            end do ! Next idims

            print *, 'fvolume_all B:', igrid, maxval(wnew(:,:,1)), maxval(wCT(:,:,1))
            if (any(wnew(:,:,1) .ne. wnew(:,:,1)) ) then
               print *, 'NaN found in wnew B'
               do idbg2=ixImin2,ixImax2
                  do idbg1=ixImin1,ixImax1
                     print *, idbg1, idbg2, wnew(idbg1,idbg2,1)
                  end do
               end do
            end if

         end associate
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
