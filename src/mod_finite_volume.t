!> Module with finite volume methods for fluxes
module mod_finite_volume
  implicit none
  private

  public :: finite_volume 
  public :: hancock
  public :: upwindLR

contains

  !> The non-conservative Hancock predictor for TVDLF
  !>
  !> on entry:
  !> input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
  !> one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
  !> on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
  subroutine hancock(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,wnew,dx^D,x)
    use mod_physics
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixI^L, ixO^L, idim^LIM
    double precision, intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), wnew(ixI^S,1:nw)

    double precision, dimension(ixI^S,1:nw) :: wprim, wLC, wRC
    double precision :: fLC(ixI^S, nwflux), fRC(ixI^S, nwflux)
    double precision :: dxinv(1:ndim),dxdim(1:ndim)
    integer :: idim, iw, ix^L, hxO^L

    ! Expand limits in each idim direction in which fluxes are added
    ix^L=ixO^L;
    do idim= idim^LIM
       ix^L=ix^L^LADDkr(idim,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in Hancock: Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    ^D&dxinv(^D)=-qdt/dx^D;
    ^D&dxdim(^D)=dx^D;
    do idim= idim^LIM
       block%iw0=idim
       ! Calculate w_j+g_j/2 and w_j-g_j/2
       ! First copy all variables, then upwind wLC and wRC.
       ! wLC is to the left of ixO, wRC is to the right of wCT.
       hxO^L=ixO^L-kr(idim,^D);

       wRC(hxO^S,1:nwflux)=wprim(ixO^S,1:nwflux)
       wLC(ixO^S,1:nwflux)=wprim(ixO^S,1:nwflux)

       call upwindLR(ixI^L,ixO^L,hxO^L,idim,wprim,wprim,wLC,wRC,x,.false.,dxdim(idim))

       ! Calculate the fLC and fRC fluxes
       call phys_get_flux(wRC,x,ixI^L,hxO^L,idim,fRC)
       call phys_get_flux(wLC,x,ixI^L,ixO^L,idim,fLC)

       ! Advect w(iw)
       do iw=1,nwflux
          if (slab) then
             wnew(ixO^S,iw)=wnew(ixO^S,iw)+dxinv(idim)* &
                  (fLC(ixO^S, iw)-fRC(hxO^S, iw))
          else
             select case (idim)
                {case (^D)
                wnew(ixO^S,iw)=wnew(ixO^S,iw)-qdt/block%dvolume(ixO^S) &
                     *(block%surfaceC^D(ixO^S)*fLC(ixO^S, iw) &
                     -block%surfaceC^D(hxO^S)*fRC(hxO^S, iw))\}
             end select
          end if
       end do
    end do ! next idim
    block%iw0=0

    if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

  end subroutine hancock

  !> finite volume method
  subroutine finite_volume(method,qdt,ixI^L,ixO^L,idim^LIM, &
       qtC,wCT,qt,wnew,wold,fC,dx^D,x)

    use mod_physics
    use mod_global_parameters
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2

    character(len=*), intent(in)                         :: method
    double precision, intent(in)                         :: qdt, qtC, qt, dx^D
    integer, intent(in)                                  :: ixI^L, ixO^L, idim^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
    double precision, dimension(ixI^S,1:ndim)             :: xi
    double precision, dimension(ixI^S,1:nw)               :: wCT, wnew, wold
    double precision, dimension(ixI^S,1:nwflux,1:ndim)  :: fC

    double precision, dimension(ixI^S,1:nw) :: wprim, wLC, wRC, wmean
    double precision, dimension(ixI^S, nwflux) :: fLC, fRC
    double precision, dimension(ixI^S)      :: cmaxC, cmaxRC, cmaxLC
    double precision, dimension(ixI^S)      :: cminC, cminRC, cminLC
    double precision, dimension(1:ndim)     :: dxinv, dxdim
    integer, dimension(ixI^S)               :: patchf
    integer :: idim, iw, ix^L, hxO^L, ixC^L, ixCR^L, jxC^L, kxC^L, kxR^L

    if (idimmax>idimmin .and. typelimited=='original')&
         call mpistop("Error in fv: Unsplit dim. and original is limited")

    fC=0.d0

    ! The flux calculation contracts by one in the idim direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ix^L=ixO^L;
    do idim= idim^LIM
       ix^L=ix^L^LADD2*kr(idim,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in fv : Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    ^D&dxinv(^D)=-qdt/dx^D;
    ^D&dxdim(^D)=dx^D;
    do idim= idim^LIM
       ! use interface value of w0 at idim
       block%iw0=idim

       hxO^L=ixO^L-kr(idim,^D);
       ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
       ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idim,^D);
       kxR^L=kxC^L+kr(idim,^D);

       wRC(kxC^S,1:nwflux)=wprim(kxR^S,1:nwflux)
       wLC(kxC^S,1:nwflux)=wprim(kxC^S,1:nwflux)

       ! Determine stencil size
       {ixCRmin^D = ixCmin^D - phys_wider_stencil\}
       {ixCRmax^D = ixCmax^D + phys_wider_stencil\}

       ! for second order scheme: apply limiting
       select case (typelimited)
       case ('previous')
          call upwindLR(ixI^L,ixCR^L,ixCR^L,idim,wold,wprim,wLC,wRC,x,.true.,dxdim(idim))
       case ('predictor')
          call upwindLR(ixI^L,ixCR^L,ixCR^L,idim,wprim,wprim,wLC,wRC,x,.false.,dxdim(idim))
       case ('original')
          call upwindLR(ixI^L,ixCR^L,ixCR^L,idim,wnew,wprim,wLC,wRC,x,.true.,dxdim(idim))
       case default
          call mpistop("Error in reconstruction: no such base for limiter")
       end select

       ! TODO wLC,wRC take primitive values may save a lot of calculation

       ! estimating bounds for the minimum and maximum signal velocities
       if(method=='tvdlf'.or.method=='tvdmu') then
         wmean(ixC^S,1:nwflux)=0.5d0*(wLC(ixC^S,1:nwflux)+wRC(ixC^S,1:nwflux))
         call phys_get_cmax(wmean,x,ixI^L,ixC^L,idim,cmaxC)
       else
         call phys_get_cbounds(wLC,wRC,x,ixI^L,ixC^L,idim,cmaxC,cminC)
       end if

       call phys_modify_wLR(wLC, wRC, ixI^L, ixC^L, idim)

       call phys_get_flux(wLC,x,ixI^L,ixC^L,idim,fLC)
       call phys_get_flux(wRC,x,ixI^L,ixC^L,idim,fRC)

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

    end do ! Next idim
    block%iw0=0

    do idim= idim^LIM
       hxO^L=ixO^L-kr(idim,^D);
       do iw=1,nwflux

          ! Multiply the fluxes by -dt/dx since Flux fixing expects this
          if (slab) then
             fC(ixI^S,iw,idim)=dxinv(idim)*fC(ixI^S,iw,idim)
             wnew(ixO^S,iw)=wnew(ixO^S,iw) &
                  + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim))
          else
             select case (idim)
                {case (^D)
                fC(ixI^S,iw,^D)=-qdt*fC(ixI^S,iw,idim)
                wnew(ixO^S,iw)=wnew(ixO^S,iw) &
                     + (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/block%dvolume(ixO^S)\}
             end select
          end if

       end do ! Next iw
       ! For the MUSCL scheme apply the characteristic based limiter
       if (method=='tvdmu') &
            call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idim,wLC,wRC,wnew,x,fC,dx^D)

    end do ! Next idim

    if (.not.slab.and.idimmin==1) &
         call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)

    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

  contains

    subroutine get_Riemann_flux_tvdmu()

      do iw=1,nwflux
         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixC^S, iw)=half*(fLC(ixC^S, iw)+fRC(ixC^S, iw))
         if (slab) then
            fC(ixC^S,iw,idim)=fLC(ixC^S, iw)
         else
            select case (idim)
               {case (^D)
               fC(ixC^S,iw,^D)=block%surfaceC^D(ixC^S)*fLC(ixC^S, iw)\}
            end select
         end if
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf()

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=1,nwflux

         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixC^S, iw)=half*(fLC(ixC^S, iw)+fRC(ixC^S, iw))

         ! Add TVDLF dissipation to the flux
         if (flux_type(idim, iw) == flux_no_dissipation) then
            fRC(ixC^S, iw)=0.d0
         else
            ! To save memory we use fRC to store -cmax*half*(w_R-w_L)
            fRC(ixC^S, iw)=-tvdlfeps*cmaxC(ixC^S)*half*(wRC(ixC^S,iw)-wLC(ixC^S,iw))
         end if

         ! fLC contains physical+dissipative fluxes
         fLC(ixC^S, iw)=fLC(ixC^S, iw)+fRC(ixC^S, iw)

         if (slab) then
            fC(ixC^S,iw,idim)=fLC(ixC^S, iw)
         else
            select case (idim)
               {case (^D)
               fC(ixC^S,iw,^D)=block%surfaceC^D(ixC^S)*fLC(ixC^S, iw)\}
            end select
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll()

      patchf(ixC^S) =  1
      where(cminC(ixC^S) >= zero)
         patchf(ixC^S) = -2
      elsewhere(cmaxC(ixC^S) <= zero)
         patchf(ixC^S) =  2
      endwhere

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=1,nwflux
         if (flux_type(idim, iw) == flux_tvdlf) then
            fLC(ixC^S, iw) = half*((fLC(ixC^S, iw) + fRC(ixC^S, iw)) &
                 -tvdlfeps*max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                 (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
         else
            where(patchf(ixC^S)==1)
               ! Add hll dissipation to the flux
               fLC(ixC^S, iw) = (cmaxC(ixC^S)*fLC(ixC^S, iw)-cminC(ixC^S) * fRC(ixC^S, iw) &
                    +tvdlfeps*cminC(ixC^S)*cmaxC(ixC^S)*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))&
                    /(cmaxC(ixC^S)-cminC(ixC^S))
            elsewhere(patchf(ixC^S)== 2)
               fLC(ixC^S, iw)=fRC(ixC^S, iw)
            elsewhere(patchf(ixC^S)==-2)
               fLC(ixC^S, iw)=fLC(ixC^S, iw)
            endwhere
         endif

         if (slab) then
            fC(ixC^S,iw,idim)=fLC(ixC^S, iw)
         else
            select case (idim)
               {case (^D)
               fC(ixC^S,iw,^D)=block%surfaceC^D(ixC^S)*fLC(ixC^S, iw)\}
            end select
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_hll

    subroutine get_Riemann_flux_hllc()
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
           call phys_diffuse_hllcd(ixI^L,ixC^L,idim,wLC,wRC,fLC,fRC,patchf)

      !---- calculate speed lambda at CD ----!
      if(any(patchf(ixC^S)==1)) &
           call phys_get_lCD(wLC,wRC,fLC,fRC,cminC,cmaxC,idim,ixI^L,ixC^L, &
           whll,Fhll,lambdaCD,patchf)

      ! now patchf may be -1 or 1 due to phys_get_lCD
      if(any(abs(patchf(ixC^S))== 1))then
         !======== flux at intermediate state ========!
         call phys_get_wCD(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,&
              cminC,cmaxC,ixI^L,ixC^L,idim,fCD)
      endif ! Calculate the CD flux

      do iw=1,nwflux
         if (flux_type(idim, iw) == flux_tvdlf) then
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

         if (slab) then
            fC(ixC^S,iw,idim)=fLC(ixC^S,iw)
         else
            select case (idim)
               {case (^D)
               fC(ixC^S,iw,^D)=block%surfaceC^D(ixC^S)*fLC(ixC^S,iw)\}
            end select
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_hllc

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315
    subroutine get_Riemann_flux_hlld()
      use mod_mhd_phys
      implicit none
      double precision, dimension(ixI^S,1:nwflux) :: w1R,w1L,f1R,f1L
      double precision, dimension(ixI^S,1:nwflux) :: w2R,w2L,f2R,f2L
      double precision, dimension(ixI^S) :: sm,s1R,s1L,suR,suL,Bx
      double precision, dimension(ixI^S) :: pts,ptR,ptL,signBx,r1L,r1R,tmp 
      double precision, dimension(ixI^S,ndir) :: vRC, vLC
      integer :: ip1,ip2,ip3,idir,ix^D

      f1R=0.d0
      f1L=0.d0
      f2R=0.d0
      f2L=0.d0
      ip1=idim
      ip3=3
      call mhd_get_v(wRC,x,ixI^L,ixC^L,vRC)
      call mhd_get_v(wLC,x,ixI^L,ixC^L,vLC)
      Bx(ixC^S)=0.5d0*(wRC(ixC^S,mag(ip1))+wLC(ixC^S,mag(ip1)))
      suR(ixC^S)=(cmaxC(ixC^S)-vRC(ixC^S,ip1))*wRC(ixC^S,rho_)
      suL(ixC^S)=(cminC(ixC^S)-vLC(ixC^S,ip1))*wLC(ixC^S,rho_)
      call mhd_get_pthermal(wRC,x,ixI^L,ixC^L,ptR)
      call mhd_get_pthermal(wLC,x,ixI^L,ixC^L,ptL)
      ptR(ixC^S)=ptR(ixC^S)+0.5d0*sum(wRC(ixC^S,mag(:))**2,dim=ndim+1)
      ptL(ixC^S)=ptL(ixC^S)+0.5d0*sum(wLC(ixC^S,mag(:))**2,dim=ndim+1)
      ! equation (38)
      sm(ixC^S)=(suR(ixC^S)*vRC(ixC^S,ip1)-suL(ixC^S)*vLC(ixC^S,ip1)-&
                 ptR(ixC^S)+ptL(ixC^S))/(suR(ixC^S)-suL(ixC^S))
      ! equation (39)
      w1R(ixC^S,mom(ip1))=sm(ixC^S)
      w1L(ixC^S,mom(ip1))=sm(ixC^S)
      w2R(ixC^S,mom(ip1))=sm(ixC^S)
      w2L(ixC^S,mom(ip1))=sm(ixC^S)
      w1R(ixC^S,mag(ip1))=Bx(ixC^S)
      w1L(ixC^S,mag(ip1))=Bx(ixC^S)
      w2R(ixC^S,mag(ip1))=Bx(ixC^S)
      w2L(ixC^S,mag(ip1))=Bx(ixC^S)
      ! equation (41)
      pts(ixC^S)=(suR(ixC^S)*ptL(ixC^S)-suL(ixC^S)*ptR(ixC^S)+suR(ixC^S)*suL(ixC^S)*&
                 (vRC(ixC^S,ip1)-vLC(ixC^S,ip1)))/(suR(ixC^S)-suL(ixC^S))
      ! equation (43)
      w1R(ixC^S,rho_)=suR(ixC^S)/(cmaxC(ixC^S)-sm(ixC^S))
      w1L(ixC^S,rho_)=suL(ixC^S)/(cminC(ixC^S)-sm(ixC^S))
      ! equation (44) ~ (47)
      ip2=mod(ip1+1,ndir)
      if(ip2==0) ip2=ndir
      r1R(ixC^S)=suR(ixC^S)*(cmaxC(ixC^S)-sm(ixC^S))-Bx(ixC^S)**2
      where(r1R(ixC^S)/=0.d0)
        r1R(ixC^S)=1.d0/r1R(ixC^S)
      endwhere
      r1L(ixC^S)=suL(ixC^S)*(cminC(ixC^S)-sm(ixC^S))-Bx(ixC^S)**2
      where(r1L(ixC^S)/=0.d0)
        r1L(ixC^S)=1.d0/r1L(ixC^S)
      endwhere
      w1R(ixC^S,mom(ip2))=vRC(ixC^S,ip2)-Bx(ixC^S)*wRC(ixC^S,mag(ip2))*&
        (sm(ixC^S)-vRC(ixC^S,ip1))*r1R(ixC^S)
      w1L(ixC^S,mom(ip2))=vLC(ixC^S,ip2)-Bx(ixC^S)*wLC(ixC^S,mag(ip2))*&
        (sm(ixC^S)-vLC(ixC^S,ip1))*r1L(ixC^S)
      w1R(ixC^S,mag(ip2))=(suR(ixC^S)*(cmaxC(ixC^S)-vRC(ixC^S,ip1))-Bx(ixC^S)**2)*r1R(ixC^S)
      w1L(ixC^S,mag(ip2))=(suL(ixC^S)*(cminC(ixC^S)-vLC(ixC^S,ip1))-Bx(ixC^S)**2)*r1L(ixC^S)
      if(ndir==3) then
        ip3=mod(ip1+2,ndir)
        if(ip3==0) ip3=ndir
        w1R(ixC^S,mom(ip3))=vRC(ixC^S,ip3)-Bx(ixC^S)*wRC(ixC^S,mag(ip3))*&
          (sm(ixC^S)-vRC(ixC^S,ip1))*r1R(ixC^S)
        w1L(ixC^S,mom(ip3))=vLC(ixC^S,ip3)-Bx(ixC^S)*wLC(ixC^S,mag(ip3))*&
          (sm(ixC^S)-vLC(ixC^S,ip1))*r1L(ixC^S)
        w1R(ixC^S,mag(ip3))=wRC(ixC^S,mag(ip3))*w1R(ixC^S,mag(ip2))
        w1L(ixC^S,mag(ip3))=wLC(ixC^S,mag(ip3))*w1L(ixC^S,mag(ip2))
      end if
      w1R(ixC^S,mag(ip2))=wRC(ixC^S,mag(ip2))*w1R(ixC^S,mag(ip2))
      w1L(ixC^S,mag(ip2))=wLC(ixC^S,mag(ip2))*w1L(ixC^S,mag(ip2))
      ! equation (48)
      if(mhd_energy) then
        w1R(ixC^S,e_)=((cmaxC(ixC^S)-vRC(ixC^S,ip1))*wRC(ixC^S,e_)-ptR(ixC^S)*vRC(ixC^S,ip1)+&
          pts(ixC^S)*sm(ixC^S)+Bx(ixC^S)*(sum(vRC(ixC^S,:)*wRC(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w1R(ixC^S,mom(:))*w1R(ixC^S,mag(:)),dim=ndim+1)))/(cmaxC(ixC^S)-sm(ixC^S))
        w1L(ixC^S,e_)=((cminC(ixC^S)-vLC(ixC^S,ip1))*wLC(ixC^S,e_)-ptL(ixC^S)*vLC(ixC^S,ip1)+&
          pts(ixC^S)*sm(ixC^S)+Bx(ixC^S)*(sum(vLC(ixC^S,:)*wLC(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w1L(ixC^S,mom(:))*w1L(ixC^S,mag(:)),dim=ndim+1)))/(cminC(ixC^S)-sm(ixC^S))
      end if
      ! equation (49)
      w2R(ixC^S,rho_)=w1R(ixC^S,rho_)
      w2L(ixC^S,rho_)=w1L(ixC^S,rho_)
      r1R(ixC^S)=sqrt(w1R(ixC^S,rho_))
      r1L(ixC^S)=sqrt(w1L(ixC^S,rho_))
      tmp(ixC^S)=1.d0/(r1R(ixC^S)+r1L(ixC^S))
      signBx(ixC^S)=sign(1.d0,Bx(ixC^S))
      ! equation (51)
      s1R(ixC^S)=sm(ixC^S)+abs(Bx(ixC^S))/r1R(ixC^S)
      s1L(ixC^S)=sm(ixC^S)-abs(Bx(ixC^S))/r1L(ixC^S)
      ! equation (59)
      w2R(ixC^S,mom(ip2))=(r1L(ixC^S)*w1L(ixC^S,mom(ip2))+r1R(ixC^S)*w1R(ixC^S,mom(ip2))+&
          (w1R(ixC^S,mag(ip2))-w1L(ixC^S,mag(ip2)))*signBx(ixC^S))*tmp(ixC^S)
      w2L(ixC^S,mom(ip2))=w2R(ixC^S,mom(ip2))
      ! equation (61)
      w2R(ixC^S,mag(ip2))=(r1L(ixC^S)*w1R(ixC^S,mag(ip2))+r1R(ixC^S)*w1L(ixC^S,mag(ip2))+&
          r1L(ixC^S)*r1R(ixC^S)*(w1R(ixC^S,mom(ip2))-w1L(ixC^S,mom(ip2)))*signBx(ixC^S))*tmp(ixC^S)
      w2L(ixC^S,mag(ip2))=w2R(ixC^S,mag(ip2))
      if(ndir==3) then
        ! equation (60)
        w2R(ixC^S,mom(ip3))=(r1L(ixC^S)*w1L(ixC^S,mom(ip3))+r1R(ixC^S)*w1R(ixC^S,mom(ip3))+&
            (w1R(ixC^S,mag(ip3))-w1L(ixC^S,mag(ip3)))*signBx(ixC^S))*tmp(ixC^S)
        w2L(ixC^S,mom(ip3))=w2R(ixC^S,mom(ip3))
        ! equation (62)
        w2R(ixC^S,mag(ip3))=(r1L(ixC^S)*w1R(ixC^S,mag(ip3))+r1R(ixC^S)*w1L(ixC^S,mag(ip3))+&
            r1L(ixC^S)*r1R(ixC^S)*(w1R(ixC^S,mom(ip3))-w1L(ixC^S,mom(ip3)))*signBx(ixC^S))*tmp(ixC^S)
        w2L(ixC^S,mag(ip3))=w2R(ixC^S,mag(ip3))
      end if
      ! equation (63)
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
      do iw=1,nwflux
        if(iw==mag(ip1)) then
          if(flux_type(idim, iw) == flux_tvdlf) &
            fC(ixC^S,iw,ip1) = 0.5d0 * (fLC(ixC^S,iw) + fRC(ixC^S,iw) - tvdlfeps * &
                 max(cmaxC(ixC^S), abs(cminC(ixC^S))) * &
                 (wRC(ixC^S,iw) - wLC(ixC^S,iw)))
          cycle
        end if
        f1L(ixC^S,iw)=fLC(ixC^S,iw)+cminC(ixC^S)*(w1L(ixC^S,iw)-wLC(ixC^S,iw))
        f1R(ixC^S,iw)=fRC(ixC^S,iw)+cmaxC(ixC^S)*(w1R(ixC^S,iw)-wRC(ixC^S,iw))
        f2L(ixC^S,iw)=f1L(ixC^S,iw)+s1L(ixC^S)*(w2L(ixC^S,iw)-w1L(ixC^S,iw))
        f2R(ixC^S,iw)=f1R(ixC^S,iw)+s1R(ixC^S)*(w2R(ixC^S,iw)-w1R(ixC^S,iw))
        where(cminC(ixC^S)>0.d0)
          fC(ixC^S,iw,ip1)=fLC(ixC^S,iw)
        else where(s1L(ixC^S)>0.d0)
          fC(ixC^S,iw,ip1)=f1L(ixC^S,iw)
        else where(sm(ixC^S)>0.d0)
          fC(ixC^S,iw,ip1)=f2L(ixC^S,iw)
        else where(s1R(ixC^S)>0.d0)
          fC(ixC^S,iw,ip1)=f2R(ixC^S,iw)
        else where(cmaxC(ixC^S)>=0.d0)
          fC(ixC^S,iw,ip1)=f1R(ixC^S,iw)
        else where(cmaxC(ixC^S)<0.d0)
          fC(ixC^S,iw,ip1)=fRC(ixC^S,iw)
        end where
        if(.not.slab) then
          select case (ip1)
            {case (^D)
            fC(ixC^S,iw,^D)=block%surfaceC^D(ixC^S)*fC(ixC^S,iw,^D)\}
          end select
        end if
      end do

    end subroutine get_Riemann_flux_hlld

  end subroutine finite_volume

  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine upwindLR(ixI^L,ixL^L,ixR^L,idim,w,wCT,wLC,wRC,x,needprim,dxdim)
    use mod_physics
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: ixI^L, ixL^L, ixR^L, idim
    logical, intent(in) :: needprim
    double precision, intent(in) :: dxdim
    double precision, dimension(ixI^S,1:nw) :: w, wCT
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    double precision, dimension(ixI^S,1:ndim) :: x

    integer :: jxR^L, ixC^L, jxC^L, iw
    double precision :: wLtmp(ixI^S,1:nw), wRtmp(ixI^S,1:nw)
    double precision :: ldw(ixI^S), dwC(ixI^S)
    integer :: flagL(ixI^S), flagR(ixI^S)
    character(std_len) :: qtypelimiter

    ! Transform w,wL,wR to primitive variables
    if (needprim) then
       call phys_to_primitive(ixI^L,ixI^L,w,x)
    end if

    if (typelimiter == 'mp5') then
       call MP5limiter(ixI^L,ixL^L,idim,w,wLC,wRC)
    else if (typelimiter == 'ppm') then
       call PPMlimiter(ixI^L,ixM^LL,idim,w,wCT,wLC,wRC)
    else
       jxR^L=ixR^L+kr(idim,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idim,^D);
       jxC^L=ixC^L+kr(idim,^D);

       qtypelimiter=typelimiter
       do iw=1,nwflux
          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
             wLC(ixL^S,iw)=dlog10(wLC(ixL^S,iw))
             wRC(ixR^S,iw)=dlog10(wRC(ixR^S,iw))
          end if

          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

          ! limit flux from left
          if(typelimiter=='koren') qtypelimiter='korenL'
          if(typelimiter=='cada')  qtypelimiter='cadaL'
          if(typelimiter=='cada3') qtypelimiter='cada3L'
          call dwlimiter2(dwC,ixI^L,ixC^L,idim,ldw,dxdim,qtypelimiter)
          wLtmp(ixL^S,iw)=wLC(ixL^S,iw)+half*ldw(ixL^S)

          ! limit flux from right
          if(typelimiter=='koren') qtypelimiter='korenR'
          if(typelimiter=='cada')  qtypelimiter='cadaR'
          if(typelimiter=='cada3') qtypelimiter='cada3R'
          call dwlimiter2(dwC,ixI^L,ixC^L,idim,ldw,dxdim,qtypelimiter)
          wRtmp(ixR^S,iw)=wRC(ixR^S,iw)-half*ldw(jxR^S)

          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
             wLtmp(ixL^S,iw)=10.0d0**wLtmp(ixL^S,iw)
             wRtmp(ixR^S,iw)=10.0d0**wRtmp(ixR^S,iw)
          end if
       end do

       call phys_check_w(.true., ixI^L, ixL^L, wLtmp, flagL)
       call phys_check_w(.true., ixI^L, ixR^L, wRtmp, flagR)

       do iw=1,nwflux
          where (flagL(ixL^S) == 0 .and. flagR(ixR^S) == 0)
             wLC(ixL^S,iw)=wLtmp(ixL^S,iw)
             wRC(ixR^S,iw)=wRtmp(ixR^S,iw)
          end where

          ! Elsewhere, we still need to convert back when using loglimit
          if (loglimit(iw)) then
             where (flagL(ixL^S) /= 0 .or. flagR(ixR^S) /= 0)
                wLC(ixL^S,iw)=10.0d0**wLC(ixL^S,iw)
                wRC(ixR^S,iw)=10.0d0**wRC(ixR^S,iw)
             end where
          end if
       enddo
    endif

    ! Transform w,wL,wR back to conservative variables
    if(needprim)then
       call phys_to_conserved(ixI^L,ixI^L,w,x)
    endif
    call phys_to_conserved(ixI^L,ixL^L,wLC,x)
    call phys_to_conserved(ixI^L,ixR^L,wRC,x)

  end subroutine upwindLR

end module mod_finite_volume
