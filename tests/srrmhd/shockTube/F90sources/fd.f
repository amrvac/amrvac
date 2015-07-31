!=============================================================================
subroutine fd(method,qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,idimmin,idimmax, qtC,&
   wCT,qt,wnew,wold,fC,dx1,x)

include 'amrvacdef.f'

character(len=*), intent(in)                                     :: method
double precision, intent(in)                                     :: qdt, qtC,&
    qt, dx1
integer, intent(in)                                              :: ixImin1,&
   ixImax1, ixOmin1,ixOmax1, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,1:ndim), intent(in)            :: &
   x

double precision, dimension(ixImin1:ixImax1,1:nw), intent(inout)           :: &
   wCT, wnew, wold
double precision, dimension(ixImin1:ixImax1,1:nwflux,1:ndim),&
    intent(out)  :: fC

double precision, dimension(ixGlo1:ixGhi1,1:nwflux)                      :: &
   fCT
double precision, dimension(ixGlo1:ixGhi1,1:nw)                          :: &
   fm, fp, fmR, fpL
double precision, dimension(ixGlo1:ixGhi1)                               :: v
double precision                                                 :: &
   dxinv(1:ndim), dxdims
logical                                                          :: transport
integer                                                          :: idims, iw,&
    ixCmin1,ixCmax1, ixmin1,ixmax1, hxOmin1,hxOmax1, ixCRmin1,ixCRmax1
!-----------------------------------------------------------------------------


dxinv(1)=-qdt/dx1;
do idims= idimmin,idimmax

   select case (idims)
      case (1) 
      dxdims = dx1
   end select
   if (B0field) then
      select case (idims)
         case (1)
         myB0 => myB0_face1
      end select
   end if

   ! Get fluxes for the whole grid (mesh+dixB)
    ixCmin1 = ixOmin1 - dixB * kr(idims,1)
    ixCmax1 = ixOmax1 + dixB * kr(idims,1)

   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
   ! ix is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
   ixmax1=ixOmax1; ixmin1=hxOmin1;


   ixCRmin1=ixCmin1;ixCRmax1=ixCmax1;


   ! Calculate velocities for transport fluxes
   call getv(wCT,x,ixGlo1,ixGhi1,ixCRmin1,ixCRmax1,idims,v)
   
   do iw=1,nwflux
      call getfluxforhllc(wCT,x,ixGlo1,ixGhi1,ixCRmin1,ixCRmax1,iw,idims,fCT,&
         transport)
      if (transport) fCT(ixCRmin1:ixCRmax1,iw) = fCT(ixCRmin1:ixCRmax1,iw) &
         + v(ixCRmin1:ixCRmax1) * wCT(ixCRmin1:ixCRmax1,iw)
      ! Lax-Friedrich splitting:
      fp(ixCRmin1:ixCRmax1,iw) = half * (fCT(ixCRmin1:ixCRmax1,iw) &
         + tvdlfeps * cmax_global * wCT(ixCRmin1:ixCRmax1,iw))
      fm(ixCRmin1:ixCRmax1,iw) = half * (fCT(ixCRmin1:ixCRmax1,iw) &
         - tvdlfeps * cmax_global * wCT(ixCRmin1:ixCRmax1,iw))
   end do ! iw loop
  
   ! now do the reconstruction of fp and fm:
   call reconstructL(ixImin1,ixImax1,ixmin1,ixmax1,idims,fp,fpL,dxdims)
   call reconstructR(ixImin1,ixImax1,ixmin1,ixmax1,idims,fm,fmR,dxdims)

   do iw=1,nwflux
      if (slab) then
         fC(ixmin1:ixmax1,iw,idims) = dxinv(idims) * (fpL(ixmin1:ixmax1,iw) &
            + fmR(ixmin1:ixmax1,iw))
         wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw)+ &
            (fC(ixOmin1:ixOmax1,iw,idims)-fC(hxOmin1:hxOmax1,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixmin1:ixmax1,iw,1)=-qdt*mygeo%surfaceC1(ixmin1:ixmax1) &
               * (fpL(ixmin1:ixmax1,iw) + fmR(ixmin1:ixmax1,iw))
            wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw)+ &
              (fC(ixOmin1:ixOmax1,iw,1)-fC(hxOmin1:hxOmax1,iw,1))&
                 /mygeo%dvolume(ixOmin1:ixOmax1)
         end select
      end if
   end do ! iw loop

end do !idims loop

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImax1,ixOmin1,&
   ixOmax1,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImax1,&
   ixOmin1,ixOmax1,1,nw,qtC,wCT,qt,wnew,x,.false.)


end subroutine fd
!=============================================================================
subroutine reconstructL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC,dxdims)
include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw), dxdims

double precision, intent(out)   :: wLC(ixGlo1:ixGhi1,1:nw) 

double precision                :: ldw(ixGlo1:ixGhi1), dwC(ixGlo1:ixGhi1)
integer                         :: jxRmin1,jxRmax1, ixCmin1,ixCmax1, jxCmin1,&
   jxCmax1, kxCmin1,kxCmax1, iw
character*79                    :: savetypelimiter
!-----------------------------------------------------------------------------
select case (typelimiter)
case ('mp5')
   call MP5limiterL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC)
case default 

   kxCmin1=ixImin1; kxCmax1=ixImax1-kr(idims,1);

   wLC(kxCmin1:kxCmax1,1:nwflux) = w(kxCmin1:kxCmax1,1:nwflux)

   jxRmin1=iLmin1+kr(idims,1);jxRmax1=iLmax1+kr(idims,1);

   ixCmax1=jxRmax1; ixCmin1=iLmin1-kr(idims,1);
   jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);

   do iw=1,nwflux
      dwC(ixCmin1:ixCmax1)=w(jxCmin1:jxCmax1,iw)-w(ixCmin1:ixCmax1,iw)
      
      savetypelimiter=typelimiter
      if(savetypelimiter=='koren') typelimiter='korenL'
      if(savetypelimiter=='cada')  typelimiter='cadaL'
      if(savetypelimiter=='cada3') typelimiter='cada3L'
      call dwlimiter2(dwC,ixCmin1,ixCmax1,iw,idims,ldw,dxdims)
      typelimiter=savetypelimiter

      wLC(iLmin1:iLmax1,iw)=wLC(iLmin1:iLmax1,iw)+half*ldw(iLmin1:iLmax1)
   end do
end select

end subroutine reconstructL
!=============================================================================
subroutine reconstructR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC,dxdims)
include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw), dxdims

double precision, intent(out)   :: wRC(ixGlo1:ixGhi1,1:nw) 

double precision                :: ldw(ixGlo1:ixGhi1), dwC(ixGlo1:ixGhi1)
integer                         :: jxRmin1,jxRmax1, ixCmin1,ixCmax1, jxCmin1,&
   jxCmax1, kxCmin1,kxCmax1, kxRmin1,kxRmax1, iw
character*79                    :: savetypelimiter
!-----------------------------------------------------------------------------
select case (typelimiter)
case ('mp5')
   call MP5limiterR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC)
case default 

   kxCmin1=ixImin1; kxCmax1=ixImax1-kr(idims,1);
   kxRmin1=kxCmin1+kr(idims,1);kxRmax1=kxCmax1+kr(idims,1);

   wRC(kxCmin1:kxCmax1,1:nwflux)=w(kxRmin1:kxRmax1,1:nwflux)

   jxRmin1=iLmin1+kr(idims,1);jxRmax1=iLmax1+kr(idims,1);
   ixCmax1=jxRmax1; ixCmin1=iLmin1-kr(idims,1);
   jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);

   do iw=1,nwflux
      dwC(ixCmin1:ixCmax1)=w(jxCmin1:jxCmax1,iw)-w(ixCmin1:ixCmax1,iw)
      
      savetypelimiter=typelimiter
      if(savetypelimiter=='koren') typelimiter='korenR'
      if(savetypelimiter=='cada')  typelimiter='cadaR'
      if(savetypelimiter=='cada3') typelimiter='cada3R'
      call dwlimiter2(dwC,ixCmin1,ixCmax1,iw,idims,ldw,dxdims)
      typelimiter=savetypelimiter

      wRC(iLmin1:iLmax1,iw)=wRC(iLmin1:iLmax1,iw)-half*ldw(jxRmin1:jxRmax1)
   end do
end select

end subroutine reconstructR
!=============================================================================
