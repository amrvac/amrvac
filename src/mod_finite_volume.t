!> Module with finite volume methods for fluxes
module mod_finite_volume
#include "amrvac.h"
  implicit none
  private

  public :: finite_volume
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

    wRp=0.d0
    wLp=0.d0
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

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,dtfactor,ixI^L,ixO^L,wCT,wprim,wnew,x)

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
          dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim),& 
         ixI^L,ixO^L,1,nw,qtC,wCT,wprim,qt,wnew,x,.false.,active)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,'exit hancock finite_volume')
    endif
    end associate
  end subroutine hancock

  !> finite volume method
  subroutine finite_volume(method,qdt,dtfactor,ixI^L,ixO^L,idims^LIM, &
       qtC,sCT,qt,snew,fC,fE,dxs,x)
    use mod_physics
    use mod_variables
    use mod_global_parameters
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2
    use mod_usr_methods
    use mod_comm_lib, only: mpistop

    integer, intent(in)                                   :: method
    double precision, intent(in)                          :: qdt, dtfactor, qtC, qt, dxs(ndim)
    integer, intent(in)                                   :: ixI^L, ixO^L, idims^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in) :: x
    double precision, dimension(ixI^S,1:nwflux,1:ndim)    :: fC
    double precision, dimension(ixI^S,sdim:3)             :: fE
    type(state)                                           :: sCT, snew

    ! primitive w at cell center
    double precision, dimension(ixI^S,1:nw) :: wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixI^S,1:nwflux) :: fLC, fRC
    double precision, dimension(ixI^S,1:number_species) :: cmaxC
    double precision, dimension(ixI^S,1:number_species) :: cminC
    double precision, dimension(ixI^S) :: Hspeed
    double precision, dimension(ixO^S) :: inv_volume
    double precision, dimension(1:ndim) :: dxinv
    integer :: idims, iw, ix^D, hx^D, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L, ii
    logical :: active
    type(ct_velocity) :: vcts

    associate(wCT=>sCT%w, wnew=>snew%w)

    ! The flux calculation contracts by one in the idims direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ix^L=ixO^L;
    do idims= idims^LIM
       ix^L=ix^L^LADD2*kr(idims,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in fv : Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)
    do idims= idims^LIM
       ! use interface value of w0 at idims
       b0i=idims

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);
       ! wRp and wLp are defined at the same locations, and will correspond to
       ! the left and right reconstructed values at a cell face. Their indexing
       ! is similar to cell-centered values, but in direction idims they are
       ! shifted half a cell towards the 'lower' direction.
       do iw=1,nwflux
        {do ix^DB=ixImin^DB,ixImax^DB\}
           ! fill all cells for averaging to fix small values
           wRp(ix^D,iw)=wprim(ix^D,iw)
           wLp(ix^D,iw)=wprim(ix^D,iw)
        {end do\}
         wRp(kxC^S,iw)=wprim(kxR^S,iw)
       end do

       hxO^L=ixO^L-kr(idims,^D);
       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D);
         ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
       end if

       ! Determine stencil size
       {ixCRmin^D = max(ixCmin^D - phys_wider_stencil,ixGlo^D)\}
       {ixCRmax^D = min(ixCmax^D + phys_wider_stencil,ixGhi^D)\}

       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(ixI^L,ixCR^L,ixCR^L,idims,wprim,wLC,wRC,wLp,wRp,x,dxs(idims))

       ! special modification of left and right status before flux evaluation
       call phys_modify_wLR(ixI^L,ixCR^L,qt,wLC,wRC,wLp,wRp,sCT,idims)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(wLC,wLp,x,ixI^L,ixC^L,idims,fLC)
       call phys_get_flux(wRC,wRp,x,ixI^L,ixC^L,idims,fRC)
       if(H_correction) then
         call phys_get_H_speed(wprim,x,ixI^L,ixO^L,idims,Hspeed)
       end if
       ! estimating bounds for the minimum and maximum signal velocities
       if(method==fs_tvdlf.or.method==fs_tvdmu) then
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixC^L,idims,Hspeed,cmaxC)
         ! index of var  velocity appears in the induction eq. 
         if(stagger_grid) call phys_get_ct_velocity(vcts,wLp,wRp,ixI^L,ixC^L,idims,cmaxC(ixI^S,index_v_mag))
       else
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixC^L,idims,Hspeed,cmaxC,cminC)
         if(stagger_grid) call phys_get_ct_velocity(vcts,wLp,wRp,ixI^L,ixC^L,idims,cmaxC(ixI^S,index_v_mag),cminC(ixI^S,index_v_mag))
       end if

       ! use approximate Riemann solver to get flux at interfaces
       select case(method)
       case(fs_hll)
         do ii=1,number_species
           call get_Riemann_flux_hll(start_indices(ii),stop_indices(ii))
         end do
       case(fs_hllc,fs_hllcd)
         do ii=1,number_species
           call get_Riemann_flux_hllc(start_indices(ii),stop_indices(ii))
         end do
       case(fs_hlld)
         do ii=1,number_species
           if(ii==index_v_mag) then
             call get_Riemann_flux_hlld(start_indices(ii),stop_indices(ii))
           else
             call get_Riemann_flux_hll(start_indices(ii),stop_indices(ii))
           endif   
         end do
       case(fs_tvdlf)
         do ii=1,number_species
           call get_Riemann_flux_tvdlf(start_indices(ii),stop_indices(ii))
         end do
       case(fs_tvdmu)
         call get_Riemann_flux_tvdmu()
       case default
         call mpistop('unkown Riemann flux in finite volume')
       end select

    end do ! Next idims
    b0i=0
    if(stagger_grid) call phys_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,snew,vcts)
    if(slab_uniform) then
      if(local_timestep) then
        dxinv(1:ndim)=-dtfactor/dxs(1:ndim)
        do idims= idims^LIM
          hx^D=kr(idims,^D)\
          hxOmin^D=ixOmin^D-hx^D\
          do iw=iwstart,nwflux
           {do ix^DB=hxOmin^DB,ixOmax^DB\}
              fC(ix^D,iw,idims)=block%dt(ix^D)*dxinv(idims)*fC(ix^D,iw,idims)
           {end do\}
           {do ix^DB=ixOmin^DB,ixOmax^DB\}
              wnew(ix^D,iw)=wnew(ix^D,iw)+fC(ix^D,iw,idims)-fC(ix^D-hx^D,iw,idims)
           {end do\}
          end do
          ! For the MUSCL scheme apply the characteristic based limiter
          if(method==fs_tvdmu) &
             call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,wLC,wRC,wnew,x,fC,dxs)
        end do
      else
        dxinv(1:ndim)=-qdt/dxs(1:ndim)
        do idims= idims^LIM
          hx^D=kr(idims,^D)\
          hxOmin^D=ixOmin^D-hx^D\
          do iw=iwstart,nwflux
           {do ix^DB=hxOmin^DB,ixOmax^DB\}
              fC(ix^D,iw,idims)=dxinv(idims)*fC(ix^D,iw,idims)
           {end do\}
           {do ix^DB=ixOmin^DB,ixOmax^DB\}
              wnew(ix^D,iw)=wnew(ix^D,iw)+fC(ix^D,iw,idims)-fC(ix^D-hx^D,iw,idims)
           {end do\}
          end do
          ! For the MUSCL scheme apply the characteristic based limiter
          if(method==fs_tvdmu) &
             call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,wLC,wRC,wnew,x,fC,dxs)
        end do
      end if
    else
      inv_volume(ixO^S) = 1.d0/block%dvolume(ixO^S)
      if(local_timestep) then
        do idims= idims^LIM
          hx^D=kr(idims,^D)\
          hxOmin^D=ixOmin^D-hx^D\
          do iw=iwstart,nwflux
           {do ix^DB=hxOmin^DB,ixOmax^DB\}
              fC(ix^D,iw,idims)=-block%dt(ix^D)*dtfactor*fC(ix^D,iw,idims)*block%surfaceC(ix^D,idims)
           {end do\}
           {do ix^DB=ixOmin^DB,ixOmax^DB\}
              wnew(ix^D,iw)=wnew(ix^D,iw)+(fC(ix^D,iw,idims)-fC(ix^D-hx^D,iw,idims))*inv_volume(ix^D)
           {end do\}
          end do
          ! For the MUSCL scheme apply the characteristic based limiter
          if (method==fs_tvdmu) &
               call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,wLC,wRC,wnew,x,fC,dxs)
        end do
      else
        do idims= idims^LIM
          hx^D=kr(idims,^D)\
          hxOmin^D=ixOmin^D-hx^D\
          do iw=iwstart,nwflux
           {do ix^DB=hxOmin^DB,ixOmax^DB\}
             fC(ix^D,iw,idims)=-qdt*fC(ix^D,iw,idims)*block%surfaceC(ix^D,idims)
           {end do\}
           {do ix^DB=ixOmin^DB,ixOmax^DB\}
              wnew(ix^D,iw)=wnew(ix^D,iw)+(fC(ix^D,iw,idims)-fC(ix^D-hx^D,iw,idims))*inv_volume(ix^D)
           {end do\}
          end do
        end do
      end if
    end if

    if (.not.slab.and.idimsmin==1) &
         call phys_add_source_geom(qdt,dtfactor,ixI^L,ixO^L,wCT,wprim,wnew,x)

    if(stagger_grid) call phys_face_to_center(ixO^L,snew)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,'multi-D finite_volume')
       if(crash) then
         ! replace erroneous values with values at previous step
         wnew=wCT
         if(stagger_grid) snew%ws=sCT%ws
       end if
    end if
 
    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),& 
         dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim),&
         ixI^L,ixO^L,1,nw,qtC,wCT,wprim,qt,wnew,x,.false.,active)

  end associate
  contains

    subroutine get_Riemann_flux_tvdmu()
      do iw=iwstart,nwflux
        fC(ixC^S,iw,idims)=half*(fLC(ixC^S,iw)+fRC(ixC^S,iw))
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf(iws,iwe)
      integer, intent(in) :: iws,iwe

      integer :: ix^D,jx^D
      double precision :: fac(ixC^S),phi

      fac(ixC^S) = -0.5d0*tvdlfeps*cmaxC(ixC^S,ii)
      do iw=iws,iwe
         fC(ixC^S,iw,idims)=0.5d0*(fLC(ixC^S, iw)+fRC(ixC^S, iw))
         ! Add TVDLF dissipation to the flux
         if(flux_type(idims, iw) /= flux_no_dissipation) then
           if(flux_adaptive_diffusion) then
            {do ix^DB=ixCmin^DB,ixCmax^DB\}
               jx^D=ix^D+kr(idims,^D)\
               !> adaptive diffusion from Rempel et al. 2009, see also Rempel et al. 2014
               !> the previous version is adopt
               if(((wRC(ix^D,iw)-wLC(ix^D,iw))*(sCT%w(jx^D,iw)-sCT%w(ix^D,iw))) .gt. 1.e-18) then
                 phi=min((wRC(ix^D,iw)-wLC(ix^D,iw))**2/((sCT%w(jx^D,iw)-sCT%w(ix^D,iw))**2+1.e-18),one)
               else
                 phi=1.d0
               end if
               fC(ix^D,iw,idims)=fC(ix^D,iw,idims)+fac(ix^D)*(wRC(ix^D,iw)-wLC(ix^D,iw))*phi
            {end do\}
           else
            {do ix^DB=ixCmin^DB,ixCmax^DB\}
               fC(ix^D,iw,idims)=fC(ix^D,iw,idims)+fac(ix^D)*(wRC(ix^D,iw)-wLC(ix^D,iw))
            {end do\}
           end if
         end if
      end do

    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll(iws,iwe)
      integer, intent(in) :: iws,iwe
      integer :: ix^D
      double precision :: phi

      if(flux_adaptive_diffusion) then
        do iw=iws,iwe
          if(flux_type(idims, iw) == flux_tvdlf) then
            if(stagger_grid) then
              ! CT MHD set zero normal B flux
              fC(ixC^S,iw,idims)=0.d0
            else
              fC(ixC^S,iw,idims)=-tvdlfeps*half*max(cmaxC(ixC^S,ii),dabs(cminC(ixC^S,ii)))*&
                   (wRC(ixC^S,iw)-wLC(ixC^S,iw))
            end if
          else
           {do ix^DB=ixCmin^DB,ixCmax^DB\}
              if(cminC(ix^D,ii) >= zero) then
                fC(ix^D,iw,idims)=fLC(ix^D,iw)
              else if(cmaxC(ix^D,ii) <= zero) then
                fC(ix^D,iw,idims)=fRC(ix^D,iw)
              else
                !> reduced diffusion is from Wang et al. 2024
                phi=max(abs(cmaxC(ix^D,ii)),abs(cminC(ix^D,ii)))/(cmaxC(ix^D,ii)-cminC(ix^D,ii))
                fC(ix^D,iw,idims)=(cmaxC(ix^D,ii)*fLC(ix^D, iw)-cminC(ix^D,ii)*fRC(ix^D,iw)&
                      +phi*cminC(ix^D,ii)*cmaxC(ix^D,ii)*(wRC(ix^D,iw)-wLC(ix^D,iw)))&
                      /(cmaxC(ix^D,ii)-cminC(ix^D,ii))
              end if
           {end do\}
          end if
        end do
      else
        do iw=iws,iwe
          if(flux_type(idims, iw) == flux_tvdlf) then
            if(stagger_grid) then
              ! CT MHD set zero normal B flux
              fC(ixC^S,iw,idims)=0.d0
            else
              fC(ixC^S,iw,idims)=-tvdlfeps*half*max(cmaxC(ixC^S,ii),dabs(cminC(ixC^S,ii)))*&
                   (wRC(ixC^S,iw)-wLC(ixC^S,iw))
            end if
          else
           {do ix^DB=ixCmin^DB,ixCmax^DB\}
              if(cminC(ix^D,ii) >= zero) then
                fC(ix^D,iw,idims)=fLC(ix^D,iw)
              else if(cmaxC(ix^D,ii) <= zero) then
                fC(ix^D,iw,idims)=fRC(ix^D,iw)
              else
                fC(ix^D,iw,idims)=(cmaxC(ix^D,ii)*fLC(ix^D, iw)-cminC(ix^D,ii)*fRC(ix^D,iw)&
                      +cminC(ix^D,ii)*cmaxC(ix^D,ii)*(wRC(ix^D,iw)-wLC(ix^D,iw)))&
                      /(cmaxC(ix^D,ii)-cminC(ix^D,ii))
              end if
           {end do\}
          end if
        end do
      end if
    end subroutine get_Riemann_flux_hll

    subroutine get_Riemann_flux_hllc(iws,iwe)
      integer, intent(in) :: iws, iwe  
      double precision, dimension(ixI^S,1:nwflux)     :: whll, Fhll, fCD
      double precision, dimension(ixI^S)              :: lambdaCD

      integer, dimension(ixI^S) :: patchf
      integer  :: rho_, p_, e_, mom(1:ndir)

      rho_ = iw_rho
      if (allocated(iw_mom)) mom(:) = iw_mom(:)
      e_ = iw_e 

      if(associated(phys_hllc_init_species)) then
       call phys_hllc_init_species(ii, rho_, mom(:), e_)
      endif  

      p_ = e_

      patchf(ixC^S) =  1
      where(cminC(ixC^S,1) >= zero)
         patchf(ixC^S) = -2
      elsewhere(cmaxC(ixC^S,1) <= zero)
         patchf(ixC^S) =  2
      endwhere
      ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4
      if(method==fs_hllcd) &
           call phys_diffuse_hllcd(ixI^L,ixC^L,idims,wLC,wRC,fLC,fRC,patchf)

      !---- calculate speed lambda at CD ----!
      if(any(patchf(ixC^S)==1)) &
           call phys_get_lCD(wLC,wRC,fLC,fRC,cminC(ixI^S,ii),cmaxC(ixI^S,ii),idims,ixI^L,ixC^L, &
           whll,Fhll,lambdaCD,patchf)

      ! now patchf may be -1 or 1 due to phys_get_lCD
      if(any(abs(patchf(ixC^S))== 1))then
         !======== flux at intermediate state ========!
         call phys_get_wCD(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,&
              cminC(ixI^S,ii),cmaxC(ixI^S,ii),ixI^L,ixC^L,idims,fCD)
      endif ! Calculate the CD flux

      do iw=iws,iwe
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixC^S,iw)=-tvdlfeps*half*max(cmaxC(ixC^S,ii),abs(cminC(ixC^S,ii))) * &
                 (wRC(ixC^S,iw) - wLC(ixC^S,iw))
         else
            where(patchf(ixC^S)==-2)
               fLC(ixC^S,iw)=fLC(ixC^S,iw)
            elsewhere(abs(patchf(ixC^S))==1)
               fLC(ixC^S,iw)=fCD(ixC^S,iw)
            elsewhere(patchf(ixC^S)==2)
               fLC(ixC^S,iw)=fRC(ixC^S,iw)
            elsewhere(patchf(ixC^S)==3)
               ! fallback option, reducing to HLL flux
               fLC(ixC^S,iw)=Fhll(ixC^S,iw)
            elsewhere(patchf(ixC^S)==4)
               ! fallback option, reducing to TVDLF flux
               fLC(ixC^S,iw) = half*((fLC(ixC^S,iw)+fRC(ixC^S,iw)) &
                    -tvdlfeps * max(cmaxC(ixC^S,ii), dabs(cminC(ixC^S,ii))) * &
                    (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
            endwhere
         end if

         fC(ixC^S,iw,idims)=fLC(ixC^S,iw)

      end do ! Next iw
    end subroutine get_Riemann_flux_hllc

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315 and Guo 2016 JCP, 327, 543
    subroutine get_Riemann_flux_hlld(iws,iwe)
      integer, intent(in) :: iws, iwe
      double precision, dimension(ixI^S,1:nwflux) :: w1R,w1L,w2R,w2L
      double precision, dimension(ixI^S) :: sm,s1R,s1L,suR,suL,Bx
      double precision, dimension(ixI^S) :: pts,ptR,ptL,signBx,r1L,r1R,tmp
      ! magnetic field from the right and the left reconstruction
      double precision, dimension(ixI^S,ndir) :: BR, BL
      integer :: ip1,ip2,ip3,idir,ix^D,^C&b^C_,^C&m^C_
      integer  :: rho_, p_, e_, mom(1:ndir), mag(1:ndir)

      associate (sR=>cmaxC,sL=>cminC)

      rho_=iw_rho
      ^C&mom(^C)=iw_mom(^C)\
      m^C_=mom(^C);
      ^C&mag(^C)=iw_mag(^C)\
      b^C_=mag(^C);
      e_ = iw_e 
      p_ = e_

      ip1=idims
      ip3=3
      if(B0field) then
        BR(ixC^S,:)=wRC(ixC^S,mag(:))+block%B0(ixC^S,:,ip1)
        BL(ixC^S,:)=wLC(ixC^S,mag(:))+block%B0(ixC^S,:,ip1)
      else
        BR(ixC^S,:)=wRC(ixC^S,mag(:))
        BL(ixC^S,:)=wLC(ixC^S,mag(:))
      end if
      if(stagger_grid) then
        Bx(ixC^S)=block%ws(ixC^S,ip1)
      else
        ! HLL estimation of normal magnetic field at cell interfaces
        ! Li, Shenghai, 2005 JCP, 203, 344, equation (33)
        Bx(ixC^S)=(sR(ixC^S,ii)*BR(ixC^S,ip1)-sL(ixC^S,ii)*BL(ixC^S,ip1))/(sR(ixC^S,ii)-sL(ixC^S,ii))
      end if
     {!DEC$ VECTOR ALWAYS
      do ix^DB=ixCmin^DB,ixCmax^DB\}
        ptR(ix^D)=wRp(ix^D,p_)+0.5d0*(^C&bR(ix^D,^C)**2+)
        ptL(ix^D)=wLp(ix^D,p_)+0.5d0*(^C&bL(ix^D,^C)**2+)
        if(iw_equi_rho>0) then
          suR(ix^D)=(sR(ix^D,ii)-wRp(ix^D,mom(ip1)))*(wRC(ix^D,rho_)+block%equi_vars(ix^D,iw_equi_rho,ip1))
          suL(ix^D)=(sL(ix^D,ii)-wLp(ix^D,mom(ip1)))*(wLC(ix^D,rho_)+block%equi_vars(ix^D,iw_equi_rho,ip1))
        else
          suR(ix^D)=(sR(ix^D,ii)-wRp(ix^D,mom(ip1)))*wRC(ix^D,rho_)
          suL(ix^D)=(sL(ix^D,ii)-wLp(ix^D,mom(ip1)))*wLC(ix^D,rho_)
        end if
        ! Miyoshi equation (38) and Guo euqation (20)
        sm(ix^D)=(suR(ix^D)*wRp(ix^D,mom(ip1))-suL(ix^D)*wLp(ix^D,mom(ip1))-&
                   ptR(ix^D)+ptL(ix^D))/(suR(ix^D)-suL(ix^D))
        ! Miyoshi equation (39) and Guo euqation (28)
        w1R(ix^D,mom(ip1))=sm(ix^D)
        w1L(ix^D,mom(ip1))=sm(ix^D)
        w2R(ix^D,mom(ip1))=sm(ix^D)
        w2L(ix^D,mom(ip1))=sm(ix^D)
        ! Guo equation (22)
        w1R(ix^D,mag(ip1))=Bx(ix^D)
        w1L(ix^D,mag(ip1))=Bx(ix^D)
        if(B0field) then
          ptR(ix^D)=wRp(ix^D,p_)+0.5d0*(^C&wRC(ix^D,b^C_)**2+)
          ptL(ix^D)=wLp(ix^D,p_)+0.5d0*(^C&wLC(ix^D,b^C_)**2+)
        end if
        ! Miyoshi equation (43) and Guo equation (27)
        w1R(ix^D,rho_)=suR(ix^D)/(sR(ix^D,ii)-sm(ix^D))
        w1L(ix^D,rho_)=suL(ix^D)/(sL(ix^D,ii)-sm(ix^D))
        ip2=mod(ip1+1,ndir)
        if(ip2==0) ip2=ndir
        r1R(ix^D)=suR(ix^D)*(sR(ix^D,ii)-sm(ix^D))-Bx(ix^D)**2
        if(r1R(ix^D)/=0.d0) r1R(ix^D)=1.d0/r1R(ix^D)
        r1L(ix^D)=suL(ix^D)*(sL(ix^D,ii)-sm(ix^D))-Bx(ix^D)**2
        if(r1L(ix^D)/=0.d0) r1L(ix^D)=1.d0/r1L(ix^D)
        ! Miyoshi equation (44)
        w1R(ix^D,mom(ip2))=wRp(ix^D,mom(ip2))-Bx(ix^D)*BR(ix^D,ip2)*&
          (sm(ix^D)-wRp(ix^D,mom(ip1)))*r1R(ix^D)
        w1L(ix^D,mom(ip2))=wLp(ix^D,mom(ip2))-Bx(ix^D)*BL(ix^D,ip2)*&
          (sm(ix^D)-wLp(ix^D,mom(ip1)))*r1L(ix^D)
        ! partial solution for later usage
        w1R(ix^D,mag(ip2))=(suR(ix^D)*(sR(ix^D,ii)-wRp(ix^D,mom(ip1)))-Bx(ix^D)**2)*r1R(ix^D)
        w1L(ix^D,mag(ip2))=(suL(ix^D)*(sL(ix^D,ii)-wLp(ix^D,mom(ip1)))-Bx(ix^D)**2)*r1L(ix^D)
        {^IFTHREEC
        ip3=mod(ip1+2,ndir)
        if(ip3==0) ip3=ndir
        ! Miyoshi equation (46)
        w1R(ix^D,mom(ip3))=wRp(ix^D,mom(ip3))-Bx(ix^D)*BR(ix^D,ip3)*&
          (sm(ix^D)-wRp(ix^D,mom(ip1)))*r1R(ix^D)
        w1L(ix^D,mom(ip3))=wLp(ix^D,mom(ip3))-Bx(ix^D)*BL(ix^D,ip3)*&
          (sm(ix^D)-wLp(ix^D,mom(ip1)))*r1L(ix^D)
        ! Miyoshi equation (47)
        w1R(ix^D,mag(ip3))=BR(ix^D,ip3)*w1R(ix^D,mag(ip2))
        w1L(ix^D,mag(ip3))=BL(ix^D,ip3)*w1L(ix^D,mag(ip2))
        }
        ! Miyoshi equation (45)
        w1R(ix^D,mag(ip2))=BR(ix^D,ip2)*w1R(ix^D,mag(ip2))
        w1L(ix^D,mag(ip2))=BL(ix^D,ip2)*w1L(ix^D,mag(ip2))
        if(B0field) then
          ! Guo equation (26)
          ^C&w1R(ix^D,b^C_)=w1R(ix^D,b^C_)-block%B0(ix^D,^C,ip1)\
          ^C&w1L(ix^D,b^C_)=w1L(ix^D,b^C_)-block%B0(ix^D,^C,ip1)\
        end if
        ! equation (48)
        if(phys_energy) then
          ! Guo equation (25) equivalent to Miyoshi equation (41)
          w1R(ix^D,p_)=suR(ix^D)*(sm(ix^D)-wRp(ix^D,mom(ip1)))+ptR(ix^D)
          w1L(ix^D,p_)=w1R(ix^D,p_)
          if(B0field) then
            ! Guo equation (32)
            w1R(ix^D,p_)=w1R(ix^D,p_)+(^C&block%B0(ix^D,^C,ip1)*(wRC(ix^D,b^C_)-w1R(ix^D,b^C_))+)
            w1L(ix^D,p_)=w1L(ix^D,p_)+(^C&block%B0(ix^D,^C,ip1)*(wLC(ix^D,b^C_)-w1L(ix^D,b^C_))+)
          end if
          ! Miyoshi equation (48) and main part of Guo euqation (31)
          w1R(ix^D,e_)=((sR(ix^D,ii)-wRp(ix^D,mom(ip1)))*wRC(ix^D,e_)-ptR(ix^D)*wRp(ix^D,mom(ip1))+&
            w1R(ix^D,p_)*sm(ix^D)+Bx(ix^D)*((^C&wRp(ix^D,m^C_)*wRC(ix^D,b^C_)+)-&
            (^C&w1R(ix^D,m^C_)*w1R(ix^D,b^C_)+)))/(sR(ix^D,ii)-sm(ix^D))
          w1L(ix^D,e_)=((sL(ix^D,ii)-wLp(ix^D,mom(ip1)))*wLC(ix^D,e_)-ptL(ix^D)*wLp(ix^D,mom(ip1))+&
            w1L(ix^D,p_)*sm(ix^D)+Bx(ix^D)*((^C&wLp(ix^D,m^C_)*wLC(ix^D,b^C_)+)-&
            (^C&w1L(ix^D,m^C_)*w1L(ix^D,b^C_)+)))/(sL(ix^D,ii)-sm(ix^D))
          if(B0field) then
            ! Guo equation (31)
            w1R(ix^D,e_)=w1R(ix^D,e_)+((^C&w1R(ix^D,b^C_)*block%B0(ix^D,^C,ip1)+)*sm(ix^D)-&
                 (^C&wRC(ix^D,b^C_)*block%B0(ix^D,^C,ip1)+)*wRp(ix^D,mom(ip1)))/(sR(ix^D,ii)-sm(ix^D))
            w1L(ix^D,e_)=w1L(ix^D,e_)+((^C&w1L(ix^D,b^C_)*block%B0(ix^D,^C,ip1)+)*sm(ix^D)-&
                 (^C&wLC(ix^D,b^C_)*block%B0(ix^D,^C,ip1)+)*wLp(ix^D,mom(ip1)))/(sL(ix^D,ii)-sm(ix^D))
          end if
          if(iw_equi_p>0) then
            w1R(ix^D,e_)=w1R(ix^D,e_)+1d0/(phys_gamma-1)*block%equi_vars(ix^D,iw_equi_p,ip1)*&
               (sm(ix^D)-wRp(ix^D,mom(ip1)))/(sR(ix^D,ii)-sm(ix^D))
            w1L(ix^D,e_)=w1L(ix^D,e_)+1d0/(phys_gamma-1)*block%equi_vars(ix^D,iw_equi_p,ip1)*&
               (sm(ix^D)-wLp(ix^D,mom(ip1)))/(sL(ix^D,ii)-sm(ix^D))
          end if
        end if

        ! Miyoshi equation (49) and Guo equation (35)
        w2R(ix^D,rho_)=w1R(ix^D,rho_)
        w2L(ix^D,rho_)=w1L(ix^D,rho_)
        w2R(ix^D,mag(ip1))=w1R(ix^D,mag(ip1))
        w2L(ix^D,mag(ip1))=w1L(ix^D,mag(ip1))
        r1R(ix^D)=sqrt(w1R(ix^D,rho_))
        r1L(ix^D)=sqrt(w1L(ix^D,rho_))
        tmp(ix^D)=1.d0/(r1R(ix^D)+r1L(ix^D))
        signBx(ix^D)=sign(1.d0,Bx(ix^D))
        ! Miyoshi equation (51) and Guo equation (33)
        s1R(ix^D)=sm(ix^D)+abs(Bx(ix^D))/r1R(ix^D)
        s1L(ix^D)=sm(ix^D)-abs(Bx(ix^D))/r1L(ix^D)
        ! Miyoshi equation (59) and Guo equation (41)
        w2R(ix^D,mom(ip2))=(r1L(ix^D)*w1L(ix^D,mom(ip2))+r1R(ix^D)*w1R(ix^D,mom(ip2))+&
            (w1R(ix^D,mag(ip2))-w1L(ix^D,mag(ip2)))*signBx(ix^D))*tmp(ix^D)
        w2L(ix^D,mom(ip2))=w2R(ix^D,mom(ip2))
        ! Miyoshi equation (61) and Guo equation (43)
        w2R(ix^D,mag(ip2))=(r1L(ix^D)*w1R(ix^D,mag(ip2))+r1R(ix^D)*w1L(ix^D,mag(ip2))+&
            r1L(ix^D)*r1R(ix^D)*(w1R(ix^D,mom(ip2))-w1L(ix^D,mom(ip2)))*signBx(ix^D))*tmp(ix^D)
        w2L(ix^D,mag(ip2))=w2R(ix^D,mag(ip2))
        {^IFTHREEC
        ! Miyoshi equation (60) and Guo equation (42)
        w2R(ix^D,mom(ip3))=(r1L(ix^D)*w1L(ix^D,mom(ip3))+r1R(ix^D)*w1R(ix^D,mom(ip3))+&
            (w1R(ix^D,mag(ip3))-w1L(ix^D,mag(ip3)))*signBx(ix^D))*tmp(ix^D)
        w2L(ix^D,mom(ip3))=w2R(ix^D,mom(ip3))
        ! Miyoshi equation (62) and Guo equation (44)
        w2R(ix^D,mag(ip3))=(r1L(ix^D)*w1R(ix^D,mag(ip3))+r1R(ix^D)*w1L(ix^D,mag(ip3))+&
            r1L(ix^D)*r1R(ix^D)*(w1R(ix^D,mom(ip3))-w1L(ix^D,mom(ip3)))*signBx(ix^D))*tmp(ix^D)
        w2L(ix^D,mag(ip3))=w2R(ix^D,mag(ip3))
        }
        ! Miyoshi equation (63) and Guo equation (45)
        if(phys_energy) then
          w2R(ix^D,e_)=w1R(ix^D,e_)+r1R(ix^D)*((^C&w1R(ix^D,m^C_)*w1R(ix^D,b^C_)+)-&
            (^C&w2R(ix^D,m^C_)*w2R(ix^D,b^C_)+))*signBx(ix^D)
          w2L(ix^D,e_)=w1L(ix^D,e_)-r1L(ix^D)*((^C&w1L(ix^D,m^C_)*w1L(ix^D,b^C_)+)-&
            (^C&w2L(ix^D,m^C_)*w2L(ix^D,b^C_)+))*signBx(ix^D)
        end if

        ! convert velocity to momentum
        ^C&w1R(ix^D,m^C_)=w1R(ix^D,m^C_)*w1R(ix^D,rho_)\
        ^C&w1L(ix^D,m^C_)=w1L(ix^D,m^C_)*w1L(ix^D,rho_)\
        ^C&w2R(ix^D,m^C_)=w2R(ix^D,m^C_)*w2R(ix^D,rho_)\
        ^C&w2L(ix^D,m^C_)=w2L(ix^D,m^C_)*w2L(ix^D,rho_)\
        if(iw_equi_rho>0) then
          w1R(ix^D,rho_)=w1R(ix^D,rho_)-block%equi_vars(ix^D,iw_equi_rho,ip1)
          w1L(ix^D,rho_)=w1L(ix^D,rho_)-block%equi_vars(ix^D,iw_equi_rho,ip1)
          w2R(ix^D,rho_)=w2R(ix^D,rho_)-block%equi_vars(ix^D,iw_equi_rho,ip1)
          w2L(ix^D,rho_)=w2L(ix^D,rho_)-block%equi_vars(ix^D,iw_equi_rho,ip1)
        end if
     {end do\}

      do iw=iws,iwe
        if(flux_type(idims, iw)==flux_special) then
          ! known flux (fLC=fRC) for normal B and psi_ in GLM method
         {!dir$ ivdep
          do ix^DB=ixCmin^DB,ixCmax^DB\}
            fC(ix^D,iw,ip1)=fLC(ix^D,iw)
         {end do\}
        else if(flux_type(idims, iw)==flux_hll) then
          ! using hll flux for tracers
         {!dir$ ivdep
          do ix^DB=ixCmin^DB,ixCmax^DB\}
            fC(ix^D,iw,ip1)=(sR(ix^D,ii)*fLC(ix^D,iw)-sL(ix^D,ii)*fRC(ix^D,iw) &
              +sR(ix^D,ii)*sL(ix^D,ii)*(wRC(ix^D,iw)-wLC(ix^D,iw)))/(sR(ix^D,ii)-sL(ix^D,ii))
         {end do\}
        else
          ! construct hlld flux
          ! Miyoshi equation (66) and Guo equation (46)
         {!dir$ ivdep
          !DEC$ VECTOR ALWAYS
          do ix^DB=ixCmin^DB,ixCmax^DB\}
            if(sL(ix^D,ii)>0.d0) then
              fC(ix^D,iw,ip1)=fLC(ix^D,iw)
            else if(s1L(ix^D)>=0.d0) then
              fC(ix^D,iw,ip1)=fLC(ix^D,iw)+sL(ix^D,ii)*(w1L(ix^D,iw)-wLC(ix^D,iw))
            else if(sm(ix^D)>=0.d0) then
              fC(ix^D,iw,ip1)=fLC(ix^D,iw)+sL(ix^D,ii)*(w1L(ix^D,iw)-wLC(ix^D,iw))+&
                s1L(ix^D)*(w2L(ix^D,iw)-w1L(ix^D,iw))
            else if(s1R(ix^D)>=0.d0) then
              fC(ix^D,iw,ip1)=fRC(ix^D,iw)+sR(ix^D,ii)*(w1R(ix^D,iw)-wRC(ix^D,iw))+&
                s1R(ix^D)*(w2R(ix^D,iw)-w1R(ix^D,iw))
            else if(sR(ix^D,ii)>=0.d0) then
              fC(ix^D,iw,ip1)=fRC(ix^D,iw)+sR(ix^D,ii)*(w1R(ix^D,iw)-wRC(ix^D,iw))
            else if(sR(ix^D,ii)<0.d0) then
              fC(ix^D,iw,ip1)=fRC(ix^D,iw)
            end if
         {end do\}
        end if
      end do

      end associate
    end subroutine get_Riemann_flux_hlld

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315 and Guo 2016 JCP, 327, 543
    !> https://arxiv.org/pdf/2108.04991.pdf
    subroutine get_Riemann_flux_hlld_mag2(iws,iwe)
      implicit none
      integer, intent(in) :: iws, iwe

      double precision, dimension(ixI^S,1:nwflux) :: w1R,w1L,f1R,f1L,f2R,f2L
      double precision, dimension(ixI^S,1:nwflux) :: w2R,w2L
      double precision, dimension(ixI^S) :: sm,s1R,s1L,suR,suL,Bx
      double precision, dimension(ixI^S) :: pts,ptR,ptL,signBx,r1L,r1R,tmp
      ! velocity from the right and the left reconstruction
      double precision, dimension(ixI^S,ndir) :: vRC, vLC
      ! magnetic field from the right and the left reconstruction
      double precision, dimension(ixI^S,ndir) :: BR, BL
      integer :: ip1,ip2,ip3,idir,ix^D
      double precision :: phiPres, thetaSM, du, dv, dw
      integer :: ixV^L, ixVb^L, ixVc^L, ixVd^L, ixVe^L, ixVf^L
      integer  :: rho_, p_, e_, mom(1:ndir), mag(1:ndir)
      double precision, parameter :: aParam = 4d0

      rho_ = iw_rho
      mom(:) = iw_mom(:)
      mag(:) = iw_mag(:) 
      p_ = iw_e
      e_ = iw_e 

      associate (sR=>cmaxC,sL=>cminC)

      f1R=0.d0
      f1L=0.d0
      ip1=idims
      ip3=3

      vRC(ixC^S,:)=wRp(ixC^S,mom(:))
      vLC(ixC^S,:)=wLp(ixC^S,mom(:))

      ! reuse s1L s1R
      call get_hlld2_modif_c(wLp,x,ixI^L,ixO^L,s1L)
      call get_hlld2_modif_c(wRp,x,ixI^L,ixO^L,s1R)
      !phiPres = min(1, maxval(max(s1L(ixO^S),s1R(ixO^S))/cmaxC(ixO^S,1)))
      phiPres = min(1d0, maxval(max(s1L(ixO^S),s1R(ixO^S)))/maxval(cmaxC(ixO^S,1)))
      phiPres = phiPres*(2D0 - phiPres)  

     !we use here not reconstructed velocity: wprim?       
     ixV^L=ixO^L;
     !first dim
     ixVmin1=ixOmin1+1  
     ixVmax1=ixOmax1+1
     du = minval(wprim(ixV^S,mom(1))-wprim(ixO^S,mom(1)))
     if(du>0d0) du=0d0  
     dv = 0d0 
     dw = 0d0 

     {^NOONED
     !second dim
     !i,j-1,k 
     ixV^L=ixO^L;
     ixVmin2=ixOmin2-1  
     ixVmax2=ixOmax2-1
 
     !i,j+1,k 
     ixVb^L=ixO^L;
     ixVbmin2=ixOmin2+1  
     ixVbmax2=ixOmax2+1

     !i+1,j,k 
     ixVc^L=ixO^L;
     ixVcmin1=ixOmin1+1  
     ixVcmax1=ixOmax1+1

     !i+1,j-1,k 
     ixVd^L=ixO^L;
     ixVdmin1=ixOmin1+1  
     ixVdmax1=ixOmax1+1
     ixVdmin2=ixOmin2-1  
     ixVdmax2=ixOmax2-1

     !i+1,j+1,k 
     ixVe^L=ixO^L;
     ixVemin1=ixOmin1+1  
     ixVemax1=ixOmax1+1
     ixVemin2=ixOmin2+1  
     ixVemax2=ixOmax2+1

     dv = minval(min(wprim(ixO^S,mom(2))-wprim(ixV^S,mom(2)),&
                  wprim(ixVb^S,mom(2))-wprim(ixO^S,mom(2)),&
                  wprim(ixVc^S,mom(2))-wprim(ixVd^S,mom(2)),&
                  wprim(ixVe^S,mom(2))-wprim(ixVc^S,mom(2))&
                ))
     if(dv>0d0) dv=0d0} 

      {^IFTHREED
     !third dim
     !i,j,k-1 
     ixV^L=ixO^L;
     ixVmin3=ixOmin3-1  
     ixVmax3=ixOmax3-1
 
     !i,j,k+1 
     ixVb^L=ixO^L;
     ixVbmin3=ixOmin3+1  
     ixVbmax3=ixOmax3+1

     !i+1,j,k 
     ixVc^L=ixO^L;
     ixVcmin1=ixOmin1+1  
     ixVcmax1=ixOmax1+1

     !i+1,j,k-1 
     ixVd^L=ixO^L;
     ixVdmin1=ixOmin1+1  
     ixVdmax1=ixOmax1+1
     ixVdmin3=ixOmin3-1  
     ixVdmax3=ixOmax3-1

     !i+1,j,k+1 
     ixVe^L=ixO^L;
     ixVemin1=ixOmin1+1  
     ixVemax1=ixOmax1+1
     ixVemin3=ixOmin3+1  
     ixVemax3=ixOmax3+1
     dw = minval(min(wprim(ixO^S,mom(3))-wprim(ixV^S,mom(3)),&
                  wprim(ixVb^S,mom(3))-wprim(ixO^S,mom(3)),&
                  wprim(ixVc^S,mom(3))-wprim(ixVd^S,mom(3)),&
                  wprim(ixVe^S,mom(3))-wprim(ixVc^S,mom(3))&
                ))
     if(dw>0d0) dw=0d0}
     thetaSM = maxval(cmaxC(ixO^S,1)) 
       
     thetaSM = (min(1d0, (thetaSM-du)/(thetaSM-min(dv,dw))))**aParam

      if(B0field) then
        BR(ixC^S,:)=wRC(ixC^S,mag(:))+block%B0(ixC^S,:,ip1)
        BL(ixC^S,:)=wLC(ixC^S,mag(:))+block%B0(ixC^S,:,ip1)
      else
        BR(ixC^S,:)=wRC(ixC^S,mag(:))
        BL(ixC^S,:)=wLC(ixC^S,mag(:))
      end if
      ! HLL estimation of normal magnetic field at cell interfaces
      Bx(ixC^S)=(sR(ixC^S,index_v_mag)*BR(ixC^S,ip1)-sL(ixC^S,index_v_mag)*BL(ixC^S,ip1)-&
                 fLC(ixC^S,mag(ip1))-fRC(ixC^S,mag(ip1)))/(sR(ixC^S,index_v_mag)-sL(ixC^S,index_v_mag))
      ptR(ixC^S)=wRp(ixC^S,p_)+0.5d0*sum(BR(ixC^S,:)**2,dim=ndim+1)
      ptL(ixC^S)=wLp(ixC^S,p_)+0.5d0*sum(BL(ixC^S,:)**2,dim=ndim+1)
      suR(ixC^S)=(sR(ixC^S,index_v_mag)-vRC(ixC^S,ip1))*wRC(ixC^S,rho_)
      suL(ixC^S)=(sL(ixC^S,index_v_mag)-vLC(ixC^S,ip1))*wLC(ixC^S,rho_)
      ! Miyoshi equation (38) and Guo euqation (20)
      sm(ixC^S)=(suR(ixC^S)*vRC(ixC^S,ip1)-suL(ixC^S)*vLC(ixC^S,ip1)-&
                 thetaSM*(ptR(ixC^S)-ptL(ixC^S)) )/(suR(ixC^S)-suL(ixC^S))
      ! Miyoshi equation (39) and Guo euqation (28)
      w1R(ixC^S,mom(ip1))=sm(ixC^S)
      w1L(ixC^S,mom(ip1))=sm(ixC^S)
      w2R(ixC^S,mom(ip1))=sm(ixC^S)
      w2L(ixC^S,mom(ip1))=sm(ixC^S)
      ! Guo equation (22)
      w1R(ixC^S,mag(ip1))=Bx(ixC^S)
      w1L(ixC^S,mag(ip1))=Bx(ixC^S)
      if(B0field) then
        ptR(ixC^S)=wRp(ixC^S,p_)+0.5d0*sum(wRC(ixC^S,mag(:))**2,dim=ndim+1)
        ptL(ixC^S)=wLp(ixC^S,p_)+0.5d0*sum(wLC(ixC^S,mag(:))**2,dim=ndim+1)
      end if

      ! Miyoshi equation (43) and Guo equation (27)
      w1R(ixC^S,rho_)=suR(ixC^S)/(sR(ixC^S,index_v_mag)-sm(ixC^S))
      w1L(ixC^S,rho_)=suL(ixC^S)/(sL(ixC^S,index_v_mag)-sm(ixC^S))

      ip2=mod(ip1+1,ndir)
      if(ip2==0) ip2=ndir
      r1R(ixC^S)=suR(ixC^S)*(sR(ixC^S,index_v_mag)-sm(ixC^S))-Bx(ixC^S)**2
      where(r1R(ixC^S)/=0.d0)
        r1R(ixC^S)=1.d0/r1R(ixC^S)
      endwhere
      r1L(ixC^S)=suL(ixC^S)*(sL(ixC^S,index_v_mag)-sm(ixC^S))-Bx(ixC^S)**2
      where(r1L(ixC^S)/=0.d0)
        r1L(ixC^S)=1.d0/r1L(ixC^S)
      endwhere
      ! Miyoshi equation (44)
      w1R(ixC^S,mom(ip2))=vRC(ixC^S,ip2)-Bx(ixC^S)*BR(ixC^S,ip2)*&
        (sm(ixC^S)-vRC(ixC^S,ip1))*r1R(ixC^S)
      w1L(ixC^S,mom(ip2))=vLC(ixC^S,ip2)-Bx(ixC^S)*BL(ixC^S,ip2)*&
        (sm(ixC^S)-vLC(ixC^S,ip1))*r1L(ixC^S)
      ! partial solution for later usage
      w1R(ixC^S,mag(ip2))=(suR(ixC^S)*(sR(ixC^S,index_v_mag)-vRC(ixC^S,ip1))-Bx(ixC^S)**2)*r1R(ixC^S)
      w1L(ixC^S,mag(ip2))=(suL(ixC^S)*(sL(ixC^S,index_v_mag)-vLC(ixC^S,ip1))-Bx(ixC^S)**2)*r1L(ixC^S)
      if(ndir==3) then
        ip3=mod(ip1+2,ndir)
        if(ip3==0) ip3=ndir
        ! Miyoshi equation (46)
        w1R(ixC^S,mom(ip3))=vRC(ixC^S,ip3)-Bx(ixC^S)*BR(ixC^S,ip3)*&
          (sm(ixC^S)-vRC(ixC^S,ip1))*r1R(ixC^S)
        w1L(ixC^S,mom(ip3))=vLC(ixC^S,ip3)-Bx(ixC^S)*BL(ixC^S,ip3)*&
          (sm(ixC^S)-vLC(ixC^S,ip1))*r1L(ixC^S)
        ! Miyoshi equation (47)
        w1R(ixC^S,mag(ip3))=BR(ixC^S,ip3)*w1R(ixC^S,mag(ip2))
        w1L(ixC^S,mag(ip3))=BL(ixC^S,ip3)*w1L(ixC^S,mag(ip2))
      end if
      ! Miyoshi equation (45)
      w1R(ixC^S,mag(ip2))=BR(ixC^S,ip2)*w1R(ixC^S,mag(ip2))
      w1L(ixC^S,mag(ip2))=BL(ixC^S,ip2)*w1L(ixC^S,mag(ip2))
      if(B0field) then
        ! Guo equation (26)
        w1R(ixC^S,mag(:))=w1R(ixC^S,mag(:))-block%B0(ixC^S,:,ip1)
        w1L(ixC^S,mag(:))=w1L(ixC^S,mag(:))-block%B0(ixC^S,:,ip1)
      end if
      ! equation (48)
      if(phys_energy) then
        ! Guo equation (25) equivalent to Miyoshi equation (41)
        w1R(ixC^S,p_)=(suR(ixC^S)*ptL(ixC^S) - suL(ixC^S)*ptR(ixC^S) +&
          phiPres * suR(ixC^S)*suL(ixC^S)*(vRC(ixC^S,ip1)-vLC(ixC^S,ip1)))/&
          (suR(ixC^S)-suL(ixC^S))
        w1L(ixC^S,p_)=w1R(ixC^S,p_)
        if(B0field) then
          ! Guo equation (32)
          w1R(ixC^S,p_)=w1R(ixC^S,p_)+sum(block%B0(ixC^S,:,ip1)*(wRC(ixC^S,mag(:))-w1R(ixC^S,mag(:))),dim=ndim+1)
          w1L(ixC^S,p_)=w1L(ixC^S,p_)+sum(block%B0(ixC^S,:,ip1)*(wLC(ixC^S,mag(:))-w1L(ixC^S,mag(:))),dim=ndim+1)
        end if
        ! Miyoshi equation (48) and main part of Guo euqation (31)
        w1R(ixC^S,e_)=((sR(ixC^S,index_v_mag)-vRC(ixC^S,ip1))*wRC(ixC^S,e_)-ptR(ixC^S)*vRC(ixC^S,ip1)+&
          w1R(ixC^S,p_)*sm(ixC^S)+Bx(ixC^S)*(sum(vRC(ixC^S,:)*wRC(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w1R(ixC^S,mom(:))*w1R(ixC^S,mag(:)),dim=ndim+1)))/(sR(ixC^S,index_v_mag)-sm(ixC^S))
        w1L(ixC^S,e_)=((sL(ixC^S,index_v_mag)-vLC(ixC^S,ip1))*wLC(ixC^S,e_)-ptL(ixC^S)*vLC(ixC^S,ip1)+&
          w1L(ixC^S,p_)*sm(ixC^S)+Bx(ixC^S)*(sum(vLC(ixC^S,:)*wLC(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w1L(ixC^S,mom(:))*w1L(ixC^S,mag(:)),dim=ndim+1)))/(sL(ixC^S,index_v_mag)-sm(ixC^S))
        if(B0field) then
          ! Guo equation (31)
          w1R(ixC^S,e_)=w1R(ixC^S,e_)+(sum(w1R(ixC^S,mag(:))*block%B0(ixC^S,:,ip1),dim=ndim+1)*sm(ixC^S)-&
               sum(wRC(ixC^S,mag(:))*block%B0(ixC^S,:,ip1),dim=ndim+1)*vRC(ixC^S,ip1))/(sR(ixC^S,index_v_mag)-sm(ixC^S))
          w1L(ixC^S,e_)=w1L(ixC^S,e_)+(sum(w1L(ixC^S,mag(:))*block%B0(ixC^S,:,ip1),dim=ndim+1)*sm(ixC^S)-&
               sum(wLC(ixC^S,mag(:))*block%B0(ixC^S,:,ip1),dim=ndim+1)*vLC(ixC^S,ip1))/(sL(ixC^S,index_v_mag)-sm(ixC^S))
        end if
      end if

      ! Miyoshi equation (49) and Guo equation (35)
      w2R(ixC^S,rho_)=w1R(ixC^S,rho_)
      w2L(ixC^S,rho_)=w1L(ixC^S,rho_)
      w2R(ixC^S,mag(ip1))=w1R(ixC^S,mag(ip1))
      w2L(ixC^S,mag(ip1))=w1L(ixC^S,mag(ip1))

      r1R(ixC^S)=sqrt(w1R(ixC^S,rho_))
      r1L(ixC^S)=sqrt(w1L(ixC^S,rho_))
      tmp(ixC^S)=1.d0/(r1R(ixC^S)+r1L(ixC^S))
      signBx(ixC^S)=sign(1.d0,Bx(ixC^S))
      ! Miyoshi equation (51) and Guo equation (33)
      s1R(ixC^S)=sm(ixC^S)+abs(Bx(ixC^S))/r1R(ixC^S)
      s1L(ixC^S)=sm(ixC^S)-abs(Bx(ixC^S))/r1L(ixC^S)
      ! Miyoshi equation (59) and Guo equation (41)
      w2R(ixC^S,mom(ip2))=(r1L(ixC^S)*w1L(ixC^S,mom(ip2))+r1R(ixC^S)*w1R(ixC^S,mom(ip2))+&
          (w1R(ixC^S,mag(ip2))-w1L(ixC^S,mag(ip2)))*signBx(ixC^S))*tmp(ixC^S)
      w2L(ixC^S,mom(ip2))=w2R(ixC^S,mom(ip2))
      ! Miyoshi equation (61) and Guo equation (43)
      w2R(ixC^S,mag(ip2))=(r1L(ixC^S)*w1R(ixC^S,mag(ip2))+r1R(ixC^S)*w1L(ixC^S,mag(ip2))+&
          r1L(ixC^S)*r1R(ixC^S)*(w1R(ixC^S,mom(ip2))-w1L(ixC^S,mom(ip2)))*signBx(ixC^S))*tmp(ixC^S)
      w2L(ixC^S,mag(ip2))=w2R(ixC^S,mag(ip2))
      if(ndir==3) then
        ! Miyoshi equation (60) and Guo equation (42)
        w2R(ixC^S,mom(ip3))=(r1L(ixC^S)*w1L(ixC^S,mom(ip3))+r1R(ixC^S)*w1R(ixC^S,mom(ip3))+&
            (w1R(ixC^S,mag(ip3))-w1L(ixC^S,mag(ip3)))*signBx(ixC^S))*tmp(ixC^S)
        w2L(ixC^S,mom(ip3))=w2R(ixC^S,mom(ip3))
        ! Miyoshi equation (62) and Guo equation (44)
        w2R(ixC^S,mag(ip3))=(r1L(ixC^S)*w1R(ixC^S,mag(ip3))+r1R(ixC^S)*w1L(ixC^S,mag(ip3))+&
            r1L(ixC^S)*r1R(ixC^S)*(w1R(ixC^S,mom(ip3))-w1L(ixC^S,mom(ip3)))*signBx(ixC^S))*tmp(ixC^S)
        w2L(ixC^S,mag(ip3))=w2R(ixC^S,mag(ip3))
      end if
      ! Miyoshi equation (63) and Guo equation (45)
      if(phys_energy) then
        w2R(ixC^S,e_)=w1R(ixC^S,e_)+r1R(ixC^S)*(sum(w1R(ixC^S,mom(:))*w1R(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w2R(ixC^S,mom(:))*w2R(ixC^S,mag(:)),dim=ndim+1))*signBx(ixC^S)
        w2L(ixC^S,e_)=w1L(ixC^S,e_)-r1L(ixC^S)*(sum(w1L(ixC^S,mom(:))*w1L(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w2L(ixC^S,mom(:))*w2L(ixC^S,mag(:)),dim=ndim+1))*signBx(ixC^S)
      end if

      ! convert velocity to momentum
      do idir=1,ndir
        w1R(ixC^S,mom(idir))=w1R(ixC^S,mom(idir))*w1R(ixC^S,rho_)
        w1L(ixC^S,mom(idir))=w1L(ixC^S,mom(idir))*w1L(ixC^S,rho_)
        w2R(ixC^S,mom(idir))=w2R(ixC^S,mom(idir))*w2R(ixC^S,rho_)
        w2L(ixC^S,mom(idir))=w2L(ixC^S,mom(idir))*w2L(ixC^S,rho_)
      end do

      ! get fluxes of intermedate states
      do iw=iws,iwe
        ! CT MHD does not need normal B flux
        if(stagger_grid .and. flux_type(idims, iw) == flux_tvdlf) cycle
        if(flux_type(idims, iw) == flux_special) then
          ! known flux (fLC=fRC) for normal B and psi_ in GLM method
          f1L(ixC^S,iw)=fLC(ixC^S,iw)
          f1R(ixC^S,iw)=f1L(ixC^S,iw)
          f2L(ixC^S,iw)=f1L(ixC^S,iw)
          f2R(ixC^S,iw)=f1L(ixC^S,iw)
        else if(flux_type(idims, iw) == flux_hll) then
          ! using hll flux for tracers
          f1L(ixC^S,iw)=(sR(ixC^S,index_v_mag)*fLC(ixC^S, iw)-sL(ixC^S,index_v_mag)*fRC(ixC^S, iw) &
                    +sR(ixC^S,index_v_mag)*sL(ixC^S,index_v_mag)*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))/(sR(ixC^S,index_v_mag)-sL(ixC^S,index_v_mag))
          f1R(ixC^S,iw)=f1L(ixC^S,iw)
          f2L(ixC^S,iw)=f1L(ixC^S,iw)
          f2R(ixC^S,iw)=f1L(ixC^S,iw)
        else
          f1L(ixC^S,iw)=fLC(ixC^S,iw)+sL(ixC^S,index_v_mag)*(w1L(ixC^S,iw)-wLC(ixC^S,iw))
          f1R(ixC^S,iw)=fRC(ixC^S,iw)+sR(ixC^S,index_v_mag)*(w1R(ixC^S,iw)-wRC(ixC^S,iw))
          f2L(ixC^S,iw)=f1L(ixC^S,iw)+s1L(ixC^S)*(w2L(ixC^S,iw)-w1L(ixC^S,iw))
          f2R(ixC^S,iw)=f1R(ixC^S,iw)+s1R(ixC^S)*(w2R(ixC^S,iw)-w1R(ixC^S,iw))
        end if
      end do

      ! Miyoshi equation (66) and Guo equation (46)
     {do ix^DB=ixCmin^DB,ixCmax^DB\}
        if(sL(ix^D,index_v_mag)>0.d0) then
          fC(ix^D,iws:iwe,ip1)=fLC(ix^D,iws:iwe)
        else if(s1L(ix^D)>=0.d0) then
          fC(ix^D,iws:iwe,ip1)=f1L(ix^D,iws:iwe)
        else if(sm(ix^D)>=0.d0) then
          fC(ix^D,iws:iwe,ip1)=f2L(ix^D,iws:iwe)
        else if(s1R(ix^D)>=0.d0) then
          fC(ix^D,iws:iwe,ip1)=f2R(ix^D,iws:iwe)
        else if(sR(ix^D,index_v_mag)>=0.d0) then
          fC(ix^D,iws:iwe,ip1)=f1R(ix^D,iws:iwe)
        else if(sR(ix^D,index_v_mag)<0.d0) then
          fC(ix^D,iws:iwe,ip1)=fRC(ix^D,iws:iwe)
        end if
     {end do\}

      end associate
    end subroutine get_Riemann_flux_hlld_mag2

    !> Calculate fast magnetosonic wave speed
    subroutine get_hlld2_modif_c(w,x,ixI^L,ixO^L,csound)
      use mod_global_parameters

      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
      double precision, intent(out):: csound(ixI^S)
      double precision :: cfast2(ixI^S), AvMinCs2(ixI^S), b2(ixI^S), kmax
      double precision :: inv_rho(ixO^S), gamma_A2(ixO^S)
      integer  :: rho_, p_, e_, mom(1:ndir), mag(1:ndir)

        rho_ = iw_rho
        mom(:) = iw_mom(:)
        mag(:) = iw_mag(:) 
        p_ = iw_e
        e_ = iw_e 

      inv_rho=1.d0/w(ixO^S,rho_)

      ! store |B|^2 in v

      if (B0field) then
        b2(ixO^S) = sum((w(ixO^S, mag(:))+block%B0(ixO^S,:,b0i))**2, dim=ndim+1)
      else
        b2(ixO^S) = sum(w(ixO^S, mag(:))**2, dim=ndim+1)
      end if


      if (B0field) then
        AvMinCs2= w(ixO^S, mag(idims))+block%B0(ixO^S,idims,b0i)
      else
        AvMinCs2= w(ixO^S, mag(idims))
      end if


      csound(ixO^S) = sum(w(ixO^S, mom(:))**2, dim=ndim+1)

      cfast2(ixO^S)   = b2(ixO^S) * inv_rho+csound(ixO^S)
      AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
           * AvMinCs2(ixO^S)**2 &
           * inv_rho 

      where(AvMinCs2(ixO^S)<zero)
         AvMinCs2(ixO^S)=zero
      end where

      AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

     csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))

    end subroutine get_hlld2_modif_c

  end subroutine finite_volume

  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR(ixI^L,ixL^L,ixR^L,idims,w,wLC,wRC,wLp,wRp,x,dxdim)
    use mod_physics
    use mod_global_parameters
    use mod_limiter
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: ixI^L, ixL^L, ixR^L, idims
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
    double precision   :: a2max

    select case (type_limiter(block%level))
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

   wLC(ixL^S,1:nwflux)=wLp(ixL^S,1:nwflux)
   wRC(ixR^S,1:nwflux)=wRp(ixR^S,1:nwflux)
   call phys_to_conserved(ixI^L,ixL^L,wLC,x)
   call phys_to_conserved(ixI^L,ixR^L,wRC,x)
   if(nwaux>0)then
      wLp(ixL^S,nwflux+1:nwflux+nwaux)=wLC(ixL^S,nwflux+1:nwflux+nwaux)
      wRp(ixR^S,nwflux+1:nwflux+nwaux)=wRC(ixR^S,nwflux+1:nwflux+nwaux)
   endif

  end subroutine reconstruct_LR

end module mod_finite_volume
