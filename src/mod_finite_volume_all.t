!> Module with finite volume methods for fluxes
module mod_finite_volume_all
#include "amrvac.h"
  implicit none
  private

  public :: finite_volume_all
  public :: hancock
  public :: reconstruct_LR

contains

  !> The non-conservative Hancock predictor for TVDLF
  !>
  !> on entry:
  !> input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
  !> one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
  !> on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
  subroutine hancock(qdt,dtfactor,ixI^L,ixO^L,idims^LIM,qtC,sCT,qt,snew,dxs,x)
    use mod_physics
    use mod_global_parameters
    use mod_source, only: addsource2
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: ixI^L, ixO^L, idims^LIM
    double precision, intent(in) :: qdt, dtfactor,qtC, qt, dxs(ndim), x(ixI^S,1:ndim)
    type(state) :: sCT, snew

    double precision, dimension(ixI^S,1:nw) :: wprim, wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixO^S)      :: inv_volume
    double precision :: fLC(ixI^S, nwflux), fRC(ixI^S, nwflux)
    double precision :: dxinv(1:ndim)
    integer :: idims, iw, ix^L, hxO^L
    logical :: active

    associate(wCT=>sCT%w,wnew=>snew%w)
      ! Expand limits in each idims direction in which fluxes are added
      ix^L=ixO^L;
      do idims= idims^LIM
         ix^L=ix^L^LADDkr(idims,^D);
      end do
      if (ixI^L^LTix^L|.or.|.or.) &
           call mpistop("Error in Hancock: Nonconforming input limits")

      wprim=wCT
      call phys_to_primitive(ixI^L,ixI^L,wprim,x)

      dxinv=-qdt/dxs
      if(.not.slab_uniform) inv_volume = 1.d0/block%dvolume(ixO^S)
      do idims= idims^LIM
         b0i=idims
         ! Calculate w_j+g_j/2 and w_j-g_j/2
         ! First copy all variables, then upwind wLC and wRC.
         ! wLC is to the left of ixO, wRC is to the right of wCT.
         hxO^L=ixO^L-kr(idims,^D);

         wRp(hxO^S,1:nwflux)=wprim(ixO^S,1:nwflux)
         wLp(ixO^S,1:nwflux)=wprim(ixO^S,1:nwflux)

         ! apply limited reconstruction for left and right status at cell interfaces
         call reconstruct_LR(ixI^L,ixO^L,hxO^L,idims,wprim,wLC,wRC,wLp,wRp,x,dxs(idims))

         ! Calculate the fLC and fRC fluxes
         call phys_get_flux(wRC,wRp,x,ixI^L,hxO^L,idims,fRC)
         call phys_get_flux(wLC,wLp,x,ixI^L,ixO^L,idims,fLC)

         ! Advect w(iw)
         if (slab_uniform) then
            if(local_timestep) then
               do iw=1,nwflux
                  wnew(ixO^S,iw)=wnew(ixO^S,iw)-block%dt(ixO^S)*dtfactor/dxs(idims)* &
                       (fLC(ixO^S, iw)-fRC(hxO^S, iw))
               end do
            else  
               do iw=1,nwflux
                  wnew(ixO^S,iw)=wnew(ixO^S,iw)+dxinv(idims)* &
                       (fLC(ixO^S, iw)-fRC(hxO^S, iw))
               end do
            endif
         else
            if(local_timestep) then
               do iw=1,nwflux
                  wnew(ixO^S,iw)=wnew(ixO^S,iw) - block%dt(ixO^S)*dtfactor * inv_volume &
                       *(block%surfaceC(ixO^S,idims)*fLC(ixO^S, iw) &
                       -block%surfaceC(hxO^S,idims)*fRC(hxO^S, iw))
               end do
            else
               do iw=1,nwflux
                  wnew(ixO^S,iw)=wnew(ixO^S,iw) - qdt * inv_volume &
                       *(block%surfaceC(ixO^S,idims)*fLC(ixO^S, iw) &
                       -block%surfaceC(hxO^S,idims)*fRC(hxO^S, iw))
               end do
            end if
         end if
      end do ! next idims
      b0i=0

      if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,dtfactor,ixI^L,ixO^L,wCT,wnew,x)

      call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
           dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim),& 
           ixI^L,ixO^L,1,nw,qtC,wCT,wprim,qt,wnew,x,.false.,active)

      ! check and optionally correct unphysical values
      if(fix_small_values) then
         call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,'exit hancock finite_volume')
      endif
    end associate
  end subroutine hancock

  !AGILE
  ! psa was ps(igrid)
  !> finite volume method computing all blocks
  subroutine finite_volume_all(method,qdt,dtfactor,ixI^L,ixO^L,idims^LIM, &
       qtC,psa,bga,qt,psb,bgb,fC,fE)
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
    double precision, dimension(ixI^S,1:nwflux,1:ndim)    :: fC
    double precision, dimension(ixI^S,sdim:3)             :: fE

    
    ! primitive w at cell center
    double precision, dimension(ixI^S,1:nw) :: wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    !$acc declare create(wprim, wLC, wRC, wLp, wRp)
    double precision, dimension(ixI^S,1:nwflux) :: fLC, fRC
    !$acc declare create(fLC, FRC)
    double precision, dimension(ixI^S,1:number_species)      :: cmaxC
    double precision, dimension(ixI^S,1:number_species)      :: cminC
    double precision, dimension(ixI^S)      :: Hspeed
    !$acc declare create(cmaxC, cminC, Hspeed)
    double precision, dimension(ixO^S)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv
    integer, dimension(ixI^S)               :: patchf
    integer :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L, ii
    logical :: active
    type(ct_velocity) :: vcts
    integer :: ix^D
    double precision, dimension(ixI^S,1:nwflux)     :: whll, Fhll, fCD
    double precision, dimension(ixI^S)              :: lambdaCD
    integer  :: rho_, p_, e_, mom(1:ndir), igrid, iigrid
                
    double precision, dimension(ixI^S,1:ndim) :: x
    double precision                          :: dxs(ndim)

    !$acc parallel loop private(block)
    do iigrid=1,igridstail_active
        igrid=igrids_active(iigrid)
        block => ps(igrid)
        !$acc enter data attach(block)

        x = ps(igrid)%x
        dxs = rnode(rpdx1_:rnodehi,igrid)

        ! TODO wCT and therefor wprim references need same treatment as sCT%w
        associate(wCT=>bga%w(:^D&,:,igrid), wnew=>bgb%w(:^D&,:,igrid))
          ! The flux calculation contracts by one in the idims direction it is applied.
          ! The limiter contracts the same directions by one more, so expand ixO by 2.
          ix^L=ixO^L;
          do idims= idims^LIM
             ix^L=ix^L^LADD2*kr(idims,^D);
          end do
          if (ixI^L^LTix^L|.or.|.or.) &
               call mpistop("Error in fv : Nonconforming input limits")

          ! no longer !$acc kernels present(fC, wCT)
          fC=0.d0     ! this is updated every loop iteration (eg l293),
                      ! how to then parallelize?
          fLC=0.d0
          fRC=0.d0
          
          wprim=wCT
          ! no longer !$acc end kernels
          
          call hd_to_primitive_gpu(ixI^L,ixI^L,wprim,x)
             
          do idims= idims^LIM
             ! use interface value of w0 at idims
             b0i=idims

             hxO^L=ixO^L-kr(idims,^D);

             kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
             kxR^L=kxC^L+kr(idims,^D);

             if(stagger_grid) then
                ! ct needs all transverse cells
                ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D);
                ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
             else
                ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
                ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
             end if

             ! wRp and wLp are defined at the same locations, and will correspond to
             ! the left and right reconstructed values at a cell face. Their indexing
             ! is similar to cell-centered values, but in direction idims they are
             ! shifted half a cell towards the 'lower' direction.
          
             wRp(kxC^S,1:nw)=wprim(kxR^S,1:nw)
             wLp(kxC^S,1:nw)=wprim(kxC^S,1:nw)
             
             ! Determine stencil size
             ! FIXME: here `phys_wider_stencil` is an integer, not a function pointer
             {ixCRmin^D = max(ixCmin^D - phys_wider_stencil,ixGlo^D)\}
             {ixCRmax^D = min(ixCmax^D + phys_wider_stencil,ixGhi^D)\}

             ! apply limited reconstruction for left and right status at cell interfaces
             call reconstruct_LR_gpu(ixI^L,ixCR^L,ixCR^L,idims,wprim,wLC,wRC,wLp,wRp,x,dxs(idims),igrid)
        
#ifndef _OPENACC         
             ! special modification of left and right status before flux evaluation
             !call phys_modify_wLR(ixI^L,ixCR^L,qt,wLC,wRC,wLp,wRp,psa,idims)
#endif
             
             ! evaluate physical fluxes according to reconstructed status
             call hd_get_flux_gpu(wLC,wLp,x,ixI^L,ixC^L,idims,fLC)
             call hd_get_flux_gpu(wRC,wRp,x,ixI^L,ixC^L,idims,fRC)

#ifndef _OPENACC      
             if(H_correction) then
                call phys_get_H_speed(wprim,x,ixI^L,ixO^L,idims,Hspeed)
             end if
#endif
             
             ! estimating bounds for the minimum and maximum signal velocities
             ! if(method==fs_tvdlf.or.method==fs_tvdmu) then
             !    call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixC^L,idims,Hspeed,cmaxC)
             !    ! index of var  velocity appears in the induction eq. 
             !    if(stagger_grid) call phys_get_ct_velocity(vcts,wLp,wRp,ixI^L,ixC^L,idims,cmaxC(ixI^S,index_v_mag))

             ! else
                call hd_get_cbounds_gpu(wLC,wRC,wLp,wRp,x,ixI^L,ixC^L,idims,Hspeed,cmaxC,cminC)
             !   if(stagger_grid) call phys_get_ct_velocity(vcts,wLp,wRp,ixI^L,ixC^L,idims,cmaxC(ixI^S,index_v_mag),cminC(ixI^S,index_v_mag))
             ! end if


           ! use approximate Riemann solver to get flux at interfaces
           select case(method)
           case(fs_hll)
             do ii=1,number_species
              do iw=start_indices(ii),stop_indices(ii)
                 {do ix^DB=ixCmin^DB,ixCmax^DB\}
                 if(cminC(ix^D,ii) >= zero) then
                    fC(ix^D,iw,idims)=fLC(ix^D,iw)
                 else if(cmaxC(ix^D,ii) <= zero) then
                    fC(ix^D,iw,idims)=fRC(ix^D,iw)
                 else
                    ! Add hll dissipation to the flux
                    fC(ix^D,iw,idims)=(cmaxC(ix^D,ii)*fLC(ix^D, iw)-cminC(ix^D,ii)*fRC(ix^D,iw)&
                         +cminC(ix^D,ii)*cmaxC(ix^D,ii)*(wRC(ix^D,iw)-wLC(ix^D,iw)))&
                         /(cmaxC(ix^D,ii)-cminC(ix^D,ii))
                 end if
                 {end do\}
              end do
    !           call get_Riemann_flux_hll(start_indices(ii),stop_indices(ii))
      !         call get_Riemann_flux_hll_gpu(start_indices(ii),stop_indices(ii))
             end do
     !      case(fs_hllc,fs_hllcd)
     !        do ii=1,number_species
     !          call get_Riemann_flux_hllc(start_indices(ii),stop_indices(ii))
     !        end do
     !      case(fs_hlld)
     !        do ii=1,number_species
     !          if(ii==index_v_mag) then
     !            call get_Riemann_flux_hlld(start_indices(ii),stop_indices(ii))
     !          else
     !            call get_Riemann_flux_hll(start_indices(ii),stop_indices(ii))
     !          endif   
     !        end do
     !      case(fs_tvdlf)
     !        do ii=1,number_species
     !          call get_Riemann_flux_tvdlf(start_indices(ii),stop_indices(ii))
     !        end do
     !      case(fs_tvdmu)
     !        call get_Riemann_flux_tvdmu()
           case default
             call mpistop('unkown Riemann flux in finite volume')
           end select
           
          end do ! Next idims
          b0i=0
          ! if(stagger_grid) call phys_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,psa(igrid),psb(igrid),vcts)
          if(slab_uniform) then
             dxinv=-qdt/dxs
             do idims= idims^LIM
                hxO^L=ixO^L-kr(idims,^D);
                ! TODO maybe put if outside loop idims: but too much code is copy pasted
                ! this is also done in hancock and fd, centdiff in mod_finite_difference
                !FIXME:
    !            if(local_timestep) then
    !               do iw=iwstart,nwflux
    !                  fC(ixI^S,iw,idims)=-block%dt(ixI^S)*dtfactor/dxs(idims)*fC(ixI^S,iw,idims)
    !               end do
    !            else
                   ! Multiply the fluxes by -dt/dx since Flux fixing expects this
                !!!!$acc kernels present(fC, wnew)
                fC(ixI^S,1:nwflux,idims)=dxinv(idims)*fC(ixI^S,1:nwflux,idims)
    !            end if
                wnew(ixO^S,iwstart:nwflux)=wnew(ixO^S,iwstart:nwflux)+&
                     (fC(ixO^S,iwstart:nwflux,idims)-fC(hxO^S,iwstart:nwflux,idims))
                !!!!$acc end kernels
                
                ! For the MUSCL scheme apply the characteristic based limiter
                !FIXME: not implemented (needs to declare create further module variables)
    !            if(method==fs_tvdmu) then
    !               call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,wLC,wRC,wnew,x,fC,dxs)
                   !           print *, 'tvdllimit2 not yet available'
    !            end if
             end do ! Next idims
             
          else
             inv_volume = 1.d0/block%dvolume(ixO^S)
             do idims= idims^LIM
                hxO^L=ixO^L-kr(idims,^D);

                if(local_timestep) then
                   do iw=iwstart,nwflux
                      fC(ixI^S,iw,idims)=-block%dt(ixI^S)*dtfactor*fC(ixI^S,iw,idims)*block%surfaceC(ixI^S,idims)
                      wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims)) * &
                           inv_volume
                   end do
                else
                   do iw=iwstart,nwflux
                      fC(ixI^S,iw,idims)=-qdt*fC(ixI^S,iw,idims)*block%surfaceC(ixI^S,idims)
                      wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims)) * &
                           inv_volume
                   end do
                end if
                ! For the MUSCL scheme apply the characteristic based limiter
                !FIXME: not implemented (needs to declare create further module variables)
!                if (method==fs_tvdmu) then
!                   call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,wLC,wRC,wnew,x,fC,dxs)
!                end if

             end do ! Next idims
          end if

          !if (.not.slab.and.idimsmin==1) &
          !     call phys_add_source_geom(qdt,dtfactor,ixI^L,ixO^L,wCT,wnew,x)
          ! if(stagger_grid) call phys_face_to_center(ixO^L,psb(igrid))


          ! check and optionally correct unphysical values
          ! if(fix_small_values) then
          !   call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,'multi-D finite_volume')
          ! end if

          call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),& 
               dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim),&
               ixI^L,ixO^L,1,nw,qtC,wCT,wprim,qt,wnew,x,.false.,active)
        end associate
        !$acc exit data detach(block)
    end do   ! igrid

  end subroutine finite_volume_all


  subroutine reconstruct_LR(ixI^L,ixL^L,ixR^L,idims,w,wLC,wRC,wLp,wRp,x,dxdim)
    use mod_physics
    use mod_global_parameters
    use mod_limiter
    use mod_comm_lib, only: mpistop

    integer, value, intent(in) :: ixI^L, ixL^L, ixR^L, idims
    double precision, intent(in) :: dxdim
    ! cell center w in primitive form
    double precision, dimension(ixI^S,1:nw) :: w
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixI^S,1:ndim) :: x

    integer            :: igrid, jxR^L, ixC^L, jxC^L, iw
    double precision   :: ldw(ixI^S), rdw(ixI^S), dwC(ixI^S)
    double precision   :: a2max
    
    select case (type_limiter(node(plevel_,igrid)))
       case (limiter_mp5)
          call MP5limiter(ixI^L,ixL^L,idims,w,wLp,wRp)
       case (limiter_weno3)
          call WENO3limiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp,1)
       case (limiter_wenoyc3)
          call WENO3limiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp,2)
       case (limiter_weno5)
          call WENO5limiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp,1)
       case (limiter_weno5nm)
          call WENO5NMlimiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp,1)
       case (limiter_wenoz5)
          call WENO5limiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp,2)
       case (limiter_wenoz5nm)
          call WENO5NMlimiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp,2)
       case (limiter_wenozp5)
          call WENO5limiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp,3)
       case (limiter_wenozp5nm)
          call WENO5NMlimiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp,3)
       case (limiter_weno5cu6)
          call WENO5CU6limiter(ixI^L,ixL^L,idims,w,wLp,wRp)
       case (limiter_teno5ad)
          call TENO5ADlimiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp)
       case (limiter_weno7)
          call WENO7limiter(ixI^L,ixL^L,idims,w,wLp,wRp,1)
       case (limiter_mpweno7)
          call WENO7limiter(ixI^L,ixL^L,idims,w,wLp,wRp,2)
       case (limiter_venk)
          call venklimiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp) 
          if(fix_small_values) then
             call phys_handle_small_values(.true.,wLp,x,ixI^L,ixL^L,'reconstruct left')
             call phys_handle_small_values(.true.,wRp,x,ixI^L,ixR^L,'reconstruct right')
          end if
       case (limiter_ppm)
          call PPMlimiter(ixI^L,ixM^LL,idims,w,w,wLp,wRp)
          if(fix_small_values) then
             call phys_handle_small_values(.true.,wLp,x,ixI^L,ixL^L,'reconstruct left')
             call phys_handle_small_values(.true.,wRp,x,ixI^L,ixR^L,'reconstruct right')
          end if
    case default
       jxR^L=ixR^L+kr(idims,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);
       do iw=1,nwflux
          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
             wLp(ixL^S,iw)=dlog10(wLp(ixL^S,iw))
             wRp(ixR^S,iw)=dlog10(wRp(ixR^S,iw))
          end if

          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)
          if(need_global_a2max) then 
             a2max=a2max_global(idims)
          else
             select case(idims)
             case(1)
                a2max=schmid_rad1
                {^IFTWOD
             case(2)
                a2max=schmid_rad2}
                {^IFTHREED
             case(2)
                a2max=schmid_rad2
             case(3)
                a2max=schmid_rad3}
             case default
                call mpistop("idims is wrong in mod_limiter")
             end select
          end if

          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixI^L,ixC^L,idims,type_limiter(block%level),ldw,rdw,a2max=a2max)
          wLp(ixL^S,iw)=wLp(ixL^S,iw)+half*ldw(ixL^S)
          wRp(ixR^S,iw)=wRp(ixR^S,iw)-half*rdw(jxR^S)

          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
             wLp(ixL^S,iw)=10.0d0**wLp(ixL^S,iw)
             wRp(ixR^S,iw)=10.0d0**wRp(ixR^S,iw)
          end if
       end do
       if(fix_small_values) then
          call phys_handle_small_values(.true.,wLp,x,ixI^L,ixL^L,'reconstruct left')
          call phys_handle_small_values(.true.,wRp,x,ixI^L,ixR^L,'reconstruct right')
       end if

    end select
    
    wLC(ixL^S,1:nwflux) = wLp(ixL^S,1:nwflux)
    wRC(ixR^S,1:nwflux) = wRp(ixR^S,1:nwflux)
    
    call phys_to_conserved(ixI^L,ixL^L,wLC,x)
    call phys_to_conserved(ixI^L,ixR^L,wRC,x)

    if(nwaux>0)then
       wLp(ixL^S,nwflux+1:nwflux+nwaux) = wLC(ixL^S,nwflux+1:nwflux+nwaux)
       wRp(ixR^S,nwflux+1:nwflux+nwaux) = wRC(ixR^S,nwflux+1:nwflux+nwaux)
    endif

  end subroutine reconstruct_LR

  
  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR_gpu(ixI^L,ixL^L,ixR^L,idims,w,wLC,wRC,wLp,wRp,x,dxdim,igrid)
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
    double precision, dimension(ixI^S,1:nw) :: w
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixI^S,1:ndim) :: x

    integer            :: jxR^L, ixC^L, jxC^L, iw
    double precision   :: ldw(ixI^S), rdw(ixI^S), dwC(ixI^S)
    !!!!$acc declare create(ldw, rdw, dwC)
    double precision   :: a2max
    
    select case (type_limiter(node(plevel_,igrid)))
    case default
       jxR^L=ixR^L+kr(idims,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);
       do iw=1,nwflux
          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
             wLp(ixL^S,iw)=dlog10(wLp(ixL^S,iw))
             wRp(ixR^S,iw)=dlog10(wRp(ixR^S,iw))
          end if

          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)
          if(need_global_a2max) then 
             a2max=a2max_global(idims)
          else
             select case(idims)
             case(1)
                a2max=schmid_rad1
                {^IFTWOD
             case(2)
                a2max=schmid_rad2}
                {^IFTHREED
             case(2)
                a2max=schmid_rad2
             case(3)
                a2max=schmid_rad3}
             case default
                call mpistop("idims is wrong in mod_limiter")
             end select
          end if

          ! limit flux from left and/or right
          call dwlimiter2_gpu(dwC,ixI^L,ixC^L,idims,type_limiter(node(plevel_,igrid)),ldw,rdw,a2max=a2max)
          wLp(ixL^S,iw)=wLp(ixL^S,iw)+half*ldw(ixL^S)
          wRp(ixR^S,iw)=wRp(ixR^S,iw)-half*rdw(jxR^S)

          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
             wLp(ixL^S,iw)=10.0d0**wLp(ixL^S,iw)
             wRp(ixR^S,iw)=10.0d0**wRp(ixR^S,iw)
          end if
       end do
       ! if(fix_small_values) then
       !    call phys_handle_small_values(.true.,wLp,x,ixI^L,ixL^L,'reconstruct left')
       !    call phys_handle_small_values(.true.,wRp,x,ixI^L,ixR^L,'reconstruct right')
       ! end if

    end select
    
    wLC(ixL^S,1:nwflux) = wLp(ixL^S,1:nwflux)
    wRC(ixR^S,1:nwflux) = wRp(ixR^S,1:nwflux)
    
    call hd_to_conserved_gpu(ixI^L,ixL^L,wLC,x)
    call hd_to_conserved_gpu(ixI^L,ixR^L,wRC,x)

    if(nwaux>0)then
       wLp(ixL^S,nwflux+1:nwflux+nwaux) = wLC(ixL^S,nwflux+1:nwflux+nwaux)
       wRp(ixR^S,nwflux+1:nwflux+nwaux) = wRC(ixR^S,nwflux+1:nwflux+nwaux)
    endif

  end subroutine reconstruct_LR_gpu
  
end module mod_finite_volume_all
