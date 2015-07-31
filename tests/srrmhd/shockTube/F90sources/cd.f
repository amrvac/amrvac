!#############################################################################
! module cd
! Centered difference scheme
!=============================================================================
subroutine centdiff(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,idimmin,idimmax,qtC,&
   wCT,qt,w,fC,dx1,x)

! Advance the iws flow variables from t to t+qdt within ixO^L by centered 
! differencing in space the dw/dt+dF_i(w)/dx_i=S type equation. 
! wCT contains the time centered variables at time qtC for flux and source.
! w is the old value at qt on input and the new value at qt+qdt on output.

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, idimmin,idimmax
double precision, intent(in) :: qdt, qtC, qt, dx1
double precision :: wCT(ixImin1:ixImax1,1:nw), w(ixImin1:ixImax1,1:nw)
double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
double precision :: fC(ixImin1:ixImax1,1:nwflux,1:ndim)

double precision :: v(ixGlo1:ixGhi1,ndim), f(ixGlo1:ixGhi1)
double precision :: dxinv(1:ndim)
integer :: idims, iw, ixmin1,ixmax1, hxOmin1,hxOmax1, ixCmin1,ixCmax1,&
    jxCmin1,jxCmax1
logical :: transport
!-----------------------------------------------------------------------------

! An extra layer is needed in each direction for which fluxes are added.
ixmin1=ixOmin1;ixmax1=ixOmax1;
do idims= idimmin,idimmax
   ixmin1=ixmin1-kr(idims,1);ixmax1=ixmax1+kr(idims,1);
end do
if (ixImin1>ixmin1.or.ixImax1<ixmax1) then
   call mpistop("Error in CentDiff: Non-conforming input limits")
end if

dxinv(1)=-qdt/dx1;

! Add fluxes to w
do idims= idimmin,idimmax
   ixmin1=ixOmin1-kr(idims,1);ixmax1=ixOmax1+kr(idims,1); ixCmin1=ixmin1
   ixCmax1=ixOmax1;
   jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1); 
   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);

   call getv(wCT,x,ixImin1,ixImax1,ixmin1,ixmax1,idims,f)
   v(ixmin1:ixmax1,idims)=f(ixmin1:ixmax1)

   do iw=1,nwflux
      ! Get non-transported flux
      call getflux(wCT,x,ixImin1,ixImax1,ixmin1,ixmax1,iw,idims,f,transport)
      ! Add transport flux
      if (transport) f(ixmin1:ixmax1)=f(ixmin1:ixmax1)+v(ixmin1:ixmax1,idims)&
         *wCT(ixmin1:ixmax1,iw)
      ! Center flux to interface
      fC(ixCmin1:ixCmax1,iw,idims)=half*(f(ixCmin1:ixCmax1)&
         +f(jxCmin1:jxCmax1))
      if (slab) then
         fC(ixCmin1:ixCmax1,iw,idims)=dxinv(idims)*fC(ixCmin1:ixCmax1,iw,&
            idims)
         w(ixOmin1:ixOmax1,iw)=w(ixOmin1:ixOmax1,iw)+(fC(ixOmin1:ixOmax1,iw,&
            idims)-fC(hxOmin1:hxOmax1,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,iw,1)=-qdt*mygeo%surfaceC1(ixCmin1:ixCmax1)&
               *fC(ixCmin1:ixCmax1,iw,1)
            w(ixOmin1:ixOmax1,iw)=w(ixOmin1:ixOmax1,iw)+ &
                 (fC(ixOmin1:ixOmax1,iw,1)-fC(hxOmin1:hxOmax1,iw,1))&
                    /mygeo%dvolume(ixOmin1:ixOmax1)
         end select
      end if
   end do    !next iw
end do       !next idims

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImax1,ixOmin1,&
   ixOmax1,wCT,w,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImax1,&
   ixOmin1,ixOmax1,1,nw,qtC,wCT,qt,w,x,.false.)

end subroutine centdiff
!=============================================================================
subroutine centdiff4(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,idimmin,idimmax,qtC,&
   wCT,qt,w,fC,dx1,x)

! Advance the flow variables from t to t+qdt within ixO^L by
! fourth order centered differencing in space 
! for the dw/dt+dF_i(w)/dx_i=S type equation.
! wCT contains the time centered variables at time qtC for flux and source.
! w is the old value at qt on input and the new value at qt+qdt on output.

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, idimmin,idimmax
double precision, intent(in) :: qdt, qtC, qt, dx1
double precision :: wCT(ixImin1:ixImax1,1:nw), w(ixImin1:ixImax1,1:nw)
double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
double precision, dimension(ixImin1:ixImax1,1:ndim)             ::  xi
double precision :: fC(ixImin1:ixImax1,1:nwflux,1:ndim)

double precision :: v(ixGlo1:ixGhi1,ndim), f(ixGlo1:ixGhi1)

double precision, dimension(ixGlo1:ixGhi1,1:nw) :: wLC, wRC
double precision, dimension(ixGlo1:ixGhi1)      :: vLC, vRC, phi, cmaxLC,&
    cmaxRC

double precision :: dxinv(1:ndim), dxdim
integer :: idims, iw, ixmin1,ixmax1, hxOmin1,hxOmax1, ixCmin1,ixCmax1,&
    jxCmin1,jxCmax1, hxCmin1,hxCmax1, kxCmin1,kxCmax1, kkxCmin1,kkxCmax1,&
    kkxRmin1,kkxRmax1
logical :: transport, new_cmax, patchw(ixGlo1:ixGhi1)
!-----------------------------------------------------------------------------
! two extra layers are needed in each direction for which fluxes are added.
ixmin1=ixOmin1;ixmax1=ixOmax1;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmax1=ixmax1+2*kr(idims,1);
end do

if (ixImin1>ixmin1.or.ixImax1<ixmax1) then
   call mpistop("Error in CentDiff4: Non-conforming input limits")
end if

dxinv(1)=-qdt/dx1;

if (useprimitive) then
   ! primitive limiting:
   ! this call ensures wCT is primitive with updated auxiliaries
   call primitive(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x)
endif

! Add fluxes to w
do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      end select
   end if

   ixmin1=ixOmin1-2*kr(idims,1);ixmax1=ixOmax1+2*kr(idims,1); 
   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);

   ixCmin1=hxOmin1; ixCmax1=ixOmax1;
   hxCmin1=ixCmin1-kr(idims,1);hxCmax1=ixCmax1-kr(idims,1); 
   jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1); 
   kxCmin1=ixCmin1+2*kr(idims,1);kxCmax1=ixCmax1+2*kr(idims,1); 

   kkxCmin1=ixImin1; kkxCmax1=ixImax1-kr(idims,1);
   kkxRmin1=kkxCmin1+kr(idims,1);kkxRmax1=kkxCmax1+kr(idims,1);
   wRC(kkxCmin1:kkxCmax1,1:nwflux)=wCT(kkxRmin1:kkxRmax1,1:nwflux)
   wLC(kkxCmin1:kkxCmax1,1:nwflux)=wCT(kkxCmin1:kkxCmax1,1:nwflux)

   ! Get interface positions:
   xi(kkxCmin1:kkxCmax1,1:ndim) = x(kkxCmin1:kkxCmax1,1:ndim)
   xi(kkxCmin1:kkxCmax1,idims) = half* ( x(kkxRmin1:kkxRmax1,idims)&
      +x(kkxCmin1:kkxCmax1,idims) )

   dxdim=-qdt/dxinv(idims)
   call upwindLR(ixImin1,ixImax1,ixCmin1,ixCmax1,ixCmin1,ixCmax1,idims,wCT,&
      wCT,wLC,wRC,x,.false.,dxdim)

   ! get auxiliaries for L and R states
   if (nwaux>0.and.(.not.(useprimitive))) then
         call getaux(.true.,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,'cd4_wLC')
         call getaux(.true.,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,'cd4_wRC')
   end if

   ! Calculate velocities from upwinded values
   new_cmax=.true.
   call getcmax(new_cmax,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,cmaxLC,&
      vLC,.false.)
   call getcmax(new_cmax,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,cmaxRC,&
      vLC,.false.)
   ! now take the maximum of left and right states
   vLC(ixCmin1:ixCmax1)=max(cmaxRC(ixCmin1:ixCmax1),cmaxLC(ixCmin1:ixCmax1))

   ! Calculate velocities for centered values
   call getv(wCT,xi,ixImin1,ixImax1,ixmin1,ixmax1,idims,f)
   v(ixmin1:ixmax1,idims)=f(ixmin1:ixmax1)

if (useprimitive) then
   ! primitive limiting:
   ! this call ensures wCT is primitive with updated auxiliaries
   patchw(ixImin1:ixImax1) = .false.
   call conserve(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x,patchw)
endif

   do iw=1,nwflux
      ! Get non-transported flux
      call getflux(wCT,xi,ixImin1,ixImax1,ixmin1,ixmax1,iw,idims,f,transport)
      ! Add transport flux
      if (transport) f(ixmin1:ixmax1)=f(ixmin1:ixmax1)+v(ixmin1:ixmax1,idims)&
         *wCT(ixmin1:ixmax1,iw)
      ! Center flux to interface
      ! f_i+1/2= (-f_(i+2) +7 f_(i+1) + 7 f_i - f_(i-1))/12
      fC(ixCmin1:ixCmax1,iw,idims)=(-f(kxCmin1:kxCmax1)+7.0d0&
         *(f(jxCmin1:jxCmax1)+f(ixCmin1:ixCmax1))-f(hxCmin1:hxCmax1))/12.0d0
      ! add rempel dissipative flux, only second order version for now
      fC(ixCmin1:ixCmax1,iw,idims)=fC(ixCmin1:ixCmax1,iw,idims)&
         -half*vLC(ixCmin1:ixCmax1) *(wRC(ixCmin1:ixCmax1,iw)&
         -wLC(ixCmin1:ixCmax1,iw))

      if (slab) then
         fC(ixCmin1:ixCmax1,iw,idims)=dxinv(idims)*fC(ixCmin1:ixCmax1,iw,&
            idims)
         ! result: f_(i+1/2)-f_(i-1/2) = [-f_(i+2)+8(f_(i+1)+f_(i-1))-f_(i-2)]/12
         w(ixOmin1:ixOmax1,iw)=w(ixOmin1:ixOmax1,iw)+(fC(ixOmin1:ixOmax1,iw,&
            idims)-fC(hxOmin1:hxOmax1,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,iw,1)=-qdt*mygeo%surfaceC1(ixCmin1:ixCmax1)&
               *fC(ixCmin1:ixCmax1,iw,1)
            w(ixOmin1:ixOmax1,iw)=w(ixOmin1:ixOmax1,iw)+ &
                 (fC(ixOmin1:ixOmax1,iw,1)-fC(hxOmin1:hxOmax1,iw,1))&
                    /mygeo%dvolume(ixOmin1:ixOmax1)
         end select
      end if
   end do    !next iw
end do       !next idims

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImax1,ixOmin1,&
   ixOmax1,wCT,w,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImax1,&
   ixOmin1,ixOmax1,1,nw,qtC,wCT,qt,w,x,.false.)

end subroutine centdiff4
!=============================================================================
! end module cd
!#############################################################################
