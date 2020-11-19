!> Module with finite volume methods for fluxes
module mod_finite_volume
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
  subroutine hancock(qdt,ixI^L,ixO^L,idims^LIM,qtC,sCT,qt,snew,dx^D,x)
    use mod_physics
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixI^L, ixO^L, idims^LIM
    double precision, intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
    type(state) :: sCT, snew

    double precision, dimension(ixI^S,1:nw) :: wprim, wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixO^S)      :: inv_volume
    double precision :: fLC(ixI^S, nwflux), fRC(ixI^S, nwflux)
    double precision :: dxinv(1:ndim),dxdim(1:ndim)
    integer :: idims, iw, ix^L, hxO^L

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

    ^D&dxinv(^D)=-qdt/dx^D;
    ^D&dxdim(^D)=dx^D;
    do idims= idims^LIM
      block%iw0=idims
      ! Calculate w_j+g_j/2 and w_j-g_j/2
      ! First copy all variables, then upwind wLC and wRC.
      ! wLC is to the left of ixO, wRC is to the right of wCT.
      hxO^L=ixO^L-kr(idims,^D);

      wRp(hxO^S,1:nwflux)=wprim(ixO^S,1:nwflux)
      wLp(ixO^S,1:nwflux)=wprim(ixO^S,1:nwflux)

      ! apply limited reconstruction for left and right status at cell interfaces
      call reconstruct_LR(ixI^L,ixO^L,hxO^L,idims,wprim,wLC,wRC,wLp,wRp,x,dxdim(idims))

      ! Calculate the fLC and fRC fluxes
      call phys_get_flux(wRC,wRp,x,ixI^L,hxO^L,idims,fRC)
      call phys_get_flux(wLC,wLp,x,ixI^L,ixO^L,idims,fLC)

      ! Advect w(iw)
      do iw=1,nwflux
        if (slab_uniform) then
          if (associated(phys_iw_methods(iw)%inv_capacity)) then
            call phys_iw_methods(iw)%inv_capacity(ixI^L, ixO^L, wnew, inv_volume)
            wnew(ixO^S,iw)=wnew(ixO^S,iw)+dxinv(idims) * inv_volume * &
                 (fLC(ixO^S, iw)-fRC(hxO^S, iw))
           else
             wnew(ixO^S,iw)=wnew(ixO^S,iw)+dxinv(idims)* &
                  (fLC(ixO^S, iw)-fRC(hxO^S, iw))
           end if
        else
          if (associated(phys_iw_methods(iw)%inv_capacity)) then
            call phys_iw_methods(iw)%inv_capacity(ixI^L, ixO^L, wnew, inv_volume)
          else
            inv_volume = 1.0d0
          end if
          inv_volume = inv_volume/block%dvolume(ixO^S)

          wnew(ixO^S,iw)=wnew(ixO^S,iw) - qdt * inv_volume &
               *(block%surfaceC(ixO^S,idims)*fLC(ixO^S, iw) &
               -block%surfaceC(hxO^S,idims)*fRC(hxO^S, iw))
        end if
      end do
    end do ! next idims
    block%iw0=0

    do iw = 1, nwflux
      if (associated(phys_iw_methods(iw)%inv_capacity)) then
        ! Copy state before adding source terms
        wprim(ixO^S, iw) = wnew(ixO^S, iw)
      end if
    end do

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

    ! If there are capacity functions, now correct the added source terms
    do iw = 1, nwflux
      if (associated(phys_iw_methods(iw)%inv_capacity)) then
        call phys_iw_methods(iw)%inv_capacity(ixI^L, ixO^L, wnew, inv_volume)
        wnew(ixO^S, iw) = wprim(ixO^S, iw) + inv_volume * &
             (wnew(ixO^S, iw) - wprim(ixO^S, iw))
      end if
    end do

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,'finite_volume')
    endif
    end associate
  end subroutine hancock

  !> finite volume method
  subroutine finite_volume(method,qdt,ixI^L,ixO^L,idims^LIM, &
       qtC,sCT,qt,snew,sold,fC,fE,dx^D,x)
    use mod_physics
    use mod_global_parameters
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2
    use mod_usr_methods

    character(len=*), intent(in)                          :: method
    double precision, intent(in)                          :: qdt, qtC, qt, dx^D
    integer, intent(in)                                   :: ixI^L, ixO^L, idims^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
    type(state)                                           :: sCT, snew, sold
    double precision, dimension(ixI^S,1:nwflux,1:ndim)    :: fC
    double precision, dimension(ixI^S,7-2*ndim:3)         :: fE

    ! primitive w at cell center
    double precision, dimension(ixI^S,1:nw) :: wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixI^S,1:nwflux) :: fLC, fRC
    double precision, dimension(ixI^S)      :: cmaxC
    double precision, dimension(ixI^S)      :: cminC
    double precision, dimension(ixO^S)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv, dxdim
    ! cell-face location coordinates
    double precision, dimension(ixI^S,1:ndim) :: xi
    integer, dimension(ixI^S)               :: patchf
    integer :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L

    associate(wCT=>sCT%w, wnew=>snew%w, wold=>sold%w)

    fC=0.d0

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

    ^D&dxinv(^D)=-qdt/dx^D;
    ^D&dxdim(^D)=dx^D;
    do idims= idims^LIM
       ! use interface value of w0 at idims
       block%iw0=idims

       hxO^L=ixO^L-kr(idims,^D);

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);

       if(stagger_grid) then
          ! ct needs all transverse cells
          ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D); ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
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
       {ixCRmin^D = max(ixCmin^D - phys_wider_stencil,ixGlo^D)\}
       {ixCRmax^D = min(ixCmax^D + phys_wider_stencil,ixGhi^D)\}

       ! get cell-face coordinates
       xi=x
       xi(ixI^S,idims)=xi(ixI^S,idims)+0.5d0*sCT%dx(ixI^S,idims)

       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(ixI^L,ixCR^L,ixCR^L,idims,wprim,wLC,wRC,wLp,wRp,xi,dxdim(idims))

       ! special modification of left and right status before flux evaluation
       call phys_modify_wLR(ixI^L,ixCR^L,qt,wLC,wRC,wLp,wRp,sCT,idims)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(wLC,wLp,xi,ixI^L,ixC^L,idims,fLC)
       call phys_get_flux(wRC,wRp,xi,ixI^L,ixC^L,idims,fRC)

       ! estimating bounds for the minimum and maximum signal velocities
       if(method=='tvdlf'.or.method=='tvdmu') then
         call phys_get_cbounds(wLC,wRC,wLp,wRp,xi,ixI^L,ixC^L,idims,cmaxC)
       else
         call phys_get_cbounds(wLC,wRC,wLp,wRp,xi,ixI^L,ixC^L,idims,cmaxC,cminC)
       end if

       ! use approximate Riemann solver to get flux at interfaces
       select case(method)
       case('tvdmu')
         call get_Riemann_flux_tvdmu()
       case('tvdlf')
         call get_Riemann_flux_tvdlf()
       case('hll')
         call get_Riemann_flux_hll()
       case('hllc','hllcd')
         call get_Riemann_flux_hllc()
       case('hlld')
         call get_Riemann_flux_hlld()
       case default
         call mpistop('unkown Riemann flux')
       end select

       if(associated(usr_set_flux)) call usr_set_flux(ixI^L,ixC^L,qt,wLC,wRC,wLp,wRp,sCT,idims,fC)

    end do ! Next idims
    block%iw0=0

    if(stagger_grid) call phys_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,snew)

    do idims= idims^LIM
       hxO^L=ixO^L-kr(idims,^D);

       ! Multiply the fluxes by -dt/dx since Flux fixing expects this
       if (slab_uniform) then
          fC(ixI^S,1:nwflux,idims)=dxinv(idims)*fC(ixI^S,1:nwflux,idims)

          do iw = iwstart, nwflux
            if (associated(phys_iw_methods(iw)%inv_capacity)) then
              call phys_iw_methods(iw)%inv_capacity(ixI^L, ixO^L, wnew, inv_volume)
              wnew(ixO^S,iw)=wnew(ixO^S,iw) + inv_volume * &
                   (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
            else
              wnew(ixO^S,iw)=wnew(ixO^S,iw) &
                   + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
            end if
          end do
       else
          if (.not. angmomfix) then ! default case
            if (associated(phys_iw_methods(iw)%inv_capacity)) then
              call phys_iw_methods(iw)%inv_capacity(ixI^L, ixO^L, wnew, inv_volume)
            else
              inv_volume = 1.0d0
            end if
            inv_volume = inv_volume/block%dvolume(ixO^S)

            do iw=iwstart,nwflux
              fC(ixI^S,iw,idims)=-qdt*fC(ixI^S,iw,idims)*block%surfaceC(ixI^S,idims)
              wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims)) * &
                  inv_volume
            enddo
          else
            ! If angular momentum conserving way to solve the equations,
            ! some fluxes additions need to be treated specifically
            call phys_angmomfix(fC,x,wnew,ixI^L,ixO^L,idims)
          endif
       end if

       ! For the MUSCL scheme apply the characteristic based limiter
       if (method=='tvdmu') &
            call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,wLC,wRC,wnew,x,fC,dx^D)

    end do ! Next idims

    if (.not.slab.and.idimsmin==1) &
         call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)

    if(stagger_grid) call phys_face_to_center(ixO^L,snew)
 
    if(phys_solve_eaux) then
      call phys_energy_synchro(qdt,ixI^L,ixO^L,wCT,wnew,x)
    endif

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,'finite_volume')
    endif

  end associate
  contains

    subroutine get_Riemann_flux_tvdmu()

      do iw=iwstart,nwflux
         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixC^S, iw)=half*(fLC(ixC^S, iw)+fRC(ixC^S, iw))
         fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf()
      double precision :: fac(ixC^S)

      fac = -0.5d0*tvdlfeps*cmaxC(ixC^S)

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=iwstart,nwflux

         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixC^S, iw)=0.5d0*(fLC(ixC^S, iw)+fRC(ixC^S, iw))

         ! Add TVDLF dissipation to the flux
         if (flux_type(idims, iw) /= flux_no_dissipation) then
            fLC(ixC^S, iw)=fLC(ixC^S, iw) + fac*(wRC(ixC^S,iw)-wLC(ixC^S,iw))
         end if

         fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
      end do ! Next iw
    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll()

      double precision :: fac(ixC^S), div(ixC^S)

      where(cminC(ixC^S) >= zero)
        patchf(ixC^S) = -2
      elsewhere(cmaxC(ixC^S) <= zero)
        patchf(ixC^S) =  2
      elsewhere
        patchf(ixC^S) =  1
      endwhere

      fac = tvdlfeps*cminC(ixC^S)*cmaxC(ixC^S)
      div = 1/(cmaxC(ixC^S)-cminC(ixC^S))

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=iwstart,nwflux
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixC^S, iw) = half*(fLC(ixC^S, iw) + fRC(ixC^S, iw) &
                 -tvdlfeps*max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                 (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
         else
            where(patchf(ixC^S)==1)
               ! Add hll dissipation to the flux
               fLC(ixC^S, iw) = (cmaxC(ixC^S)*fLC(ixC^S, iw)-cminC(ixC^S) * fRC(ixC^S, iw) &
                    +fac*(wRC(ixC^S,iw)-wLC(ixC^S,iw))) * div
            elsewhere(patchf(ixC^S)== 2)
               fLC(ixC^S, iw)=fRC(ixC^S, iw)
            elsewhere(patchf(ixC^S)==-2)
               fLC(ixC^S, iw)=fLC(ixC^S, iw)
            endwhere
         endif

         fC(ixC^S,iw,idims)=fLC(ixC^S, iw)

      end do ! Next iw
    end subroutine get_Riemann_flux_hll

    subroutine get_Riemann_flux_hllc()
      use mod_mhd_phys
      implicit none
      double precision, dimension(ixI^S,1:nwflux)     :: whll, Fhll, fCD
      double precision, dimension(ixI^S)              :: lambdaCD

      patchf(ixC^S) =  1
      where(cminC(ixC^S) >= zero)
         patchf(ixC^S) = -2
      elsewhere(cmaxC(ixC^S) <= zero)
         patchf(ixC^S) =  2
      endwhere
      ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4
      if(method=='hllcd') &
           call phys_diffuse_hllcd(ixI^L,ixC^L,idims,wLC,wRC,fLC,fRC,patchf)

      !---- calculate speed lambda at CD ----!
      if(any(patchf(ixC^S)==1)) &
           call phys_get_lCD(wLC,wRC,fLC,fRC,cminC,cmaxC,idims,ixI^L,ixC^L, &
           whll,Fhll,lambdaCD,patchf)

      ! now patchf may be -1 or 1 due to phys_get_lCD
      if(any(abs(patchf(ixC^S))== 1))then
         !======== flux at intermediate state ========!
         call phys_get_wCD(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,&
              cminC,cmaxC,ixI^L,ixC^L,idims,fCD)
      endif ! Calculate the CD flux

      ! use hll flux for the auxiliary internal e
      if(mhd_energy.and.mhd_solve_eaux) then
        iw=eaux_
        fCD(ixC^S, iw) = (cmaxC(ixC^S)*fLC(ixC^S, iw)-cminC(ixC^S) * fRC(ixC^S, iw) &
             +cminC(ixC^S)*cmaxC(ixC^S)*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))/(cmaxC(ixC^S)-cminC(ixC^S))
      end if

      do iw=iwstart,nwflux
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixC^S,iw) = 0.5d0 * (fLC(ixC^S,iw) + fRC(ixC^S,iw) - tvdlfeps * &
                 max(cmaxC(ixC^S), abs(cminC(ixC^S))) * &
                 (wRC(ixC^S,iw) - wLC(ixC^S,iw)))
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
                    -tvdlfeps * max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                    (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
            endwhere
         end if

         fC(ixC^S,iw,idims)=fLC(ixC^S,iw)

      end do ! Next iw
    end subroutine get_Riemann_flux_hllc

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315 and Guo 2016 JCP, 327, 543
    subroutine get_Riemann_flux_hlld()
      use mod_mhd_phys
      implicit none
      double precision, dimension(ixI^S,1:nwflux) :: w1R,w1L,f1R,f1L,f2R,f2L
      double precision, dimension(ixI^S,1:nwflux) :: w2R,w2L
      double precision, dimension(ixI^S) :: sm,s1R,s1L,suR,suL,Bx
      double precision, dimension(ixI^S) :: pts,ptR,ptL,signBx,r1L,r1R,tmp
      ! velocity from the right and the left reconstruction
      double precision, dimension(ixI^S,ndir) :: vRC, vLC
      ! magnetic field from the right and the left reconstruction
      double precision, dimension(ixI^S,ndir) :: BR, BL
      integer :: ip1,ip2,ip3,idir,ix^D

      associate (sR=>cmaxC,sL=>cminC)

      f1R=0.d0
      f1L=0.d0
      ip1=idims
      ip3=3
      vRC(ixC^S,:)=wRp(ixC^S,mom(:))
      vLC(ixC^S,:)=wLp(ixC^S,mom(:))
      if(B0field) then
        BR(ixC^S,:)=wRC(ixC^S,mag(:))+block%B0(ixC^S,:,ip1)
        BL(ixC^S,:)=wLC(ixC^S,mag(:))+block%B0(ixC^S,:,ip1)
      else
        BR(ixC^S,:)=wRC(ixC^S,mag(:))
        BL(ixC^S,:)=wLC(ixC^S,mag(:))
      end if
      ! HLL estimation of normal magnetic field at cell interfaces
      Bx(ixC^S)=(sR(ixC^S)*BR(ixC^S,ip1)-sL(ixC^S)*BL(ixC^S,ip1)-&
                 fLC(ixC^S,mag(ip1))-fRC(ixC^S,mag(ip1)))/(sR(ixC^S)-sL(ixC^S))
      ptR(ixC^S)=wRp(ixC^S,p_)+0.5d0*sum(BR(ixC^S,:)**2,dim=ndim+1)
      ptL(ixC^S)=wLp(ixC^S,p_)+0.5d0*sum(BL(ixC^S,:)**2,dim=ndim+1)
      suR(ixC^S)=(sR(ixC^S)-vRC(ixC^S,ip1))*wRC(ixC^S,rho_)
      suL(ixC^S)=(sL(ixC^S)-vLC(ixC^S,ip1))*wLC(ixC^S,rho_)
      ! Miyoshi equation (38) and Guo euqation (20)
      sm(ixC^S)=(suR(ixC^S)*vRC(ixC^S,ip1)-suL(ixC^S)*vLC(ixC^S,ip1)-&
                 ptR(ixC^S)+ptL(ixC^S))/(suR(ixC^S)-suL(ixC^S))
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
      w1R(ixC^S,rho_)=suR(ixC^S)/(sR(ixC^S)-sm(ixC^S))
      w1L(ixC^S,rho_)=suL(ixC^S)/(sL(ixC^S)-sm(ixC^S))

      ip2=mod(ip1+1,ndir)
      if(ip2==0) ip2=ndir
      r1R(ixC^S)=suR(ixC^S)*(sR(ixC^S)-sm(ixC^S))-Bx(ixC^S)**2
      where(r1R(ixC^S)/=0.d0)
        r1R(ixC^S)=1.d0/r1R(ixC^S)
      endwhere
      r1L(ixC^S)=suL(ixC^S)*(sL(ixC^S)-sm(ixC^S))-Bx(ixC^S)**2
      where(r1L(ixC^S)/=0.d0)
        r1L(ixC^S)=1.d0/r1L(ixC^S)
      endwhere
      ! Miyoshi equation (44)
      w1R(ixC^S,mom(ip2))=vRC(ixC^S,ip2)-Bx(ixC^S)*BR(ixC^S,ip2)*&
        (sm(ixC^S)-vRC(ixC^S,ip1))*r1R(ixC^S)
      w1L(ixC^S,mom(ip2))=vLC(ixC^S,ip2)-Bx(ixC^S)*BL(ixC^S,ip2)*&
        (sm(ixC^S)-vLC(ixC^S,ip1))*r1L(ixC^S)
      ! partial solution for later usage
      w1R(ixC^S,mag(ip2))=(suR(ixC^S)*(sR(ixC^S)-vRC(ixC^S,ip1))-Bx(ixC^S)**2)*r1R(ixC^S)
      w1L(ixC^S,mag(ip2))=(suL(ixC^S)*(sL(ixC^S)-vLC(ixC^S,ip1))-Bx(ixC^S)**2)*r1L(ixC^S)
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
      if(mhd_energy) then
        ! Guo equation (25) equivalent to Miyoshi equation (41)
        w1R(ixC^S,p_)=suR(ixC^S)*(sm(ixC^S)-vRC(ixC^S,ip1))+ptR(ixC^S)
        w1L(ixC^S,p_)=suL(ixC^S)*(sm(ixC^S)-vLC(ixC^S,ip1))+ptL(ixC^S)
        if(B0field) then
          ! Guo equation (32)
          w1R(ixC^S,p_)=w1R(ixC^S,p_)+sum(block%B0(ixC^S,:,ip1)*(wRC(ixC^S,mag(:))-w1R(ixC^S,mag(:))),dim=ndim+1)
          w1L(ixC^S,p_)=w1L(ixC^S,p_)+sum(block%B0(ixC^S,:,ip1)*(wLC(ixC^S,mag(:))-w1L(ixC^S,mag(:))),dim=ndim+1)
        end if
        ! Miyoshi equation (48) and main part of Guo euqation (31)
        w1R(ixC^S,e_)=((sR(ixC^S)-vRC(ixC^S,ip1))*wRC(ixC^S,e_)-ptR(ixC^S)*vRC(ixC^S,ip1)+&
          w1R(ixC^S,p_)*sm(ixC^S)+Bx(ixC^S)*(sum(vRC(ixC^S,:)*wRC(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w1R(ixC^S,mom(:))*w1R(ixC^S,mag(:)),dim=ndim+1)))/(sR(ixC^S)-sm(ixC^S))
        w1L(ixC^S,e_)=((sL(ixC^S)-vLC(ixC^S,ip1))*wLC(ixC^S,e_)-ptL(ixC^S)*vLC(ixC^S,ip1)+&
          w1L(ixC^S,p_)*sm(ixC^S)+Bx(ixC^S)*(sum(vLC(ixC^S,:)*wLC(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w1L(ixC^S,mom(:))*w1L(ixC^S,mag(:)),dim=ndim+1)))/(sL(ixC^S)-sm(ixC^S))
        if(B0field) then
          ! Guo equation (31)
          w1R(ixC^S,e_)=w1R(ixC^S,e_)+(sum(w1R(ixC^S,mag(:))*block%B0(ixC^S,:,ip1),dim=ndim+1)*sm(ixC^S)-&
               sum(wRC(ixC^S,mag(:))*block%B0(ixC^S,:,ip1),dim=ndim+1)*vRC(ixC^S,ip1))/(sR(ixC^S)-sm(ixC^S))
          w1L(ixC^S,e_)=w1L(ixC^S,e_)+(sum(w1L(ixC^S,mag(:))*block%B0(ixC^S,:,ip1),dim=ndim+1)*sm(ixC^S)-&
               sum(wLC(ixC^S,mag(:))*block%B0(ixC^S,:,ip1),dim=ndim+1)*vLC(ixC^S,ip1))/(sL(ixC^S)-sm(ixC^S))
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
      if(mhd_energy) then
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
      do iw=iwstart,nwflux
        if (flux_type(idims, iw) == flux_tvdlf) then
          !! hll flux for normal B
          !f1L(ixC^S,iw)=(sR(ixC^S)*fLC(ixC^S, iw)-sL(ixC^S)*fRC(ixC^S, iw) &
          !          +sR(ixC^S)*sL(ixC^S)*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))/(sR(ixC^S)-sL(ixC^S))
          ! tvldf flux for normal B
          f1L(ixC^S,iw)= half*(fLC(ixC^S, iw) + fRC(ixC^S, iw) &
                 -tvdlfeps*max(sR(ixC^S), dabs(sL(ixC^S))) * &
                 (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
          f1R(ixC^S,iw)=f1L(ixC^S,iw)
          f2L(ixC^S,iw)=f1L(ixC^S,iw)
          f2R(ixC^S,iw)=f1L(ixC^S,iw)
        else if(flux_type(idims, iw) == flux_special) then
          ! known flux (fLC=fRC) for normal B and psi_ in GLM method
          f1L(ixC^S,iw)=fLC(ixC^S,iw)
          f1R(ixC^S,iw)=f1L(ixC^S,iw)
          f2L(ixC^S,iw)=f1L(ixC^S,iw)
          f2R(ixC^S,iw)=f1L(ixC^S,iw)
        else
          f1L(ixC^S,iw)=fLC(ixC^S,iw)+sL(ixC^S)*(w1L(ixC^S,iw)-wLC(ixC^S,iw))
          f1R(ixC^S,iw)=fRC(ixC^S,iw)+sR(ixC^S)*(w1R(ixC^S,iw)-wRC(ixC^S,iw))
          f2L(ixC^S,iw)=f1L(ixC^S,iw)+s1L(ixC^S)*(w2L(ixC^S,iw)-w1L(ixC^S,iw))
          f2R(ixC^S,iw)=f1R(ixC^S,iw)+s1R(ixC^S)*(w2R(ixC^S,iw)-w1R(ixC^S,iw))
        end if
      end do

      ! use hll flux for the auxiliary internal e
      if(mhd_energy.and.mhd_solve_eaux) then
        iw=eaux_
        f1L(ixC^S,iw)=(sR(ixC^S)*fLC(ixC^S, iw)-sL(ixC^S)*fRC(ixC^S, iw) &
                  +sR(ixC^S)*sL(ixC^S)*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))/(sR(ixC^S)-sL(ixC^S))
        f1R(ixC^S,iw)=f1L(ixC^S,iw)
        f2L(ixC^S,iw)=f1L(ixC^S,iw)
        f2R(ixC^S,iw)=f1L(ixC^S,iw)
      end if

      ! Miyoshi equation (66) and Guo equation (46)
     {do ix^DB=ixCmin^DB,ixCmax^DB\}
        if(sL(ix^D)>0.d0) then
          fC(ix^D,1:nwflux,ip1)=fLC(ix^D,1:nwflux)
        else if(s1L(ix^D)>=0.d0) then
          fC(ix^D,1:nwflux,ip1)=f1L(ix^D,1:nwflux)
        else if(sm(ix^D)>=0.d0) then
          fC(ix^D,1:nwflux,ip1)=f2L(ix^D,1:nwflux)
        else if(s1R(ix^D)>=0.d0) then
          fC(ix^D,1:nwflux,ip1)=f2R(ix^D,1:nwflux)
        else if(sR(ix^D)>=0.d0) then
          fC(ix^D,1:nwflux,ip1)=f1R(ix^D,1:nwflux)
        else if(sR(ix^D)<0.d0) then
          fC(ix^D,1:nwflux,ip1)=fRC(ix^D,1:nwflux)
        end if
     {end do\}

      end associate
    end subroutine get_Riemann_flux_hlld

  end subroutine finite_volume

  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR(ixI^L,ixL^L,ixR^L,idims,w,wLC,wRC,wLp,wRp,x,dxdim)
    use mod_physics
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: ixI^L, ixL^L, ixR^L, idims
    double precision, intent(in) :: dxdim
    double precision, dimension(ixI^S,1:nw) :: w
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixI^S,1:ndim) :: x

    integer            :: jxR^L, ixC^L, jxC^L, iw
    double precision   :: ldw(ixI^S), rdw(ixI^S), dwC(ixI^S)
    double precision   :: a2max

    select case (typelimiter)
    case (limiter_venk)
       call venklimiter(ixI^L,ixL^L,idims,dxdim,w,wLp,wRp) 
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
    case (limiter_exeno7)
       call exENO7limiter(ixI^L,ixL^L,idims,w,wLp,wRp)
    case (limiter_ppm)
       call PPMlimiter(ixI^L,ixM^LL,idims,w,w,wLp,wRp)
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
          call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter,ldw,rdw,a2max=a2max)
          wLp(ixL^S,iw)=wLp(ixL^S,iw)+half*ldw(ixL^S)
          wRp(ixR^S,iw)=wRp(ixR^S,iw)-half*rdw(jxR^S)

          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
             wLp(ixL^S,iw)=10.0d0**wLp(ixL^S,iw)
             wRp(ixR^S,iw)=10.0d0**wRp(ixR^S,iw)
          end if
       end do
    end select

    if(fix_small_values) then
      call phys_handle_small_values(.true.,wLp,x,ixI^L,ixL^L,'reconstruct left')
      call phys_handle_small_values(.true.,wRp,x,ixI^L,ixR^L,'reconstruct right')
    end if
    wLC(ixL^S,1:nw)=wLp(ixL^S,1:nw)
    wRC(ixR^S,1:nw)=wRp(ixR^S,1:nw)
    call phys_to_conserved(ixI^L,ixL^L,wLC,x)
    call phys_to_conserved(ixI^L,ixR^L,wRC,x)

  end subroutine reconstruct_LR

end module mod_finite_volume
