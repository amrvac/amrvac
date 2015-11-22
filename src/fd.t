!=============================================================================
subroutine fd(method,qdt,ixI^L,ixO^L,idim^LIM, &
                     qtC,wCT,qt,wnew,wold,fC,dx^D,x)

include 'amrvacdef.f'

character(len=*), intent(in)                                     :: method
double precision, intent(in)                                     :: qdt, qtC, qt, dx^D
integer, intent(in)                                              :: ixI^L, ixO^L, idim^LIM
double precision, dimension(ixI^S,1:ndim), intent(in)            :: x

double precision, dimension(ixI^S,1:nw), intent(inout)           :: wCT, wnew, wold
double precision, dimension(ixI^S,1:nwflux,1:ndim), intent(out)  :: fC

double precision, dimension(ixI^S,1:nwflux)                      :: fCT
double precision, dimension(ixI^S,1:nw)                          :: fm, fp, fmR, fpL
double precision, dimension(ixI^S)                               :: v
double precision                                                 :: dxinv(1:ndim), dxdims
logical                                                          :: transport
integer                                                          :: idims, iw, ixC^L, ix^L, hxO^L, ixCR^L
!-----------------------------------------------------------------------------


^D&dxinv(^D)=-qdt/dx^D;
do idims= idim^LIM

   select case (idims)
      {case (^D) 
      dxdims = dx^D\}
   end select
   if (B0field) then
      myB0 => myB0_cell
   end if

   ! Get fluxes for the whole grid (mesh+dixB)
   {^D& ixCmin^D = ixOmin^D - dixB * kr(idims,^D)\}
   {^D& ixCmax^D = ixOmax^D + dixB * kr(idims,^D)\}

   hxO^L=ixO^L-kr(idims,^D);
   ! ix is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
   ixmax^D=ixOmax^D; ixmin^D=hxOmin^D;

{#IFDEF HALL
   ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
   ! assuming one additional ghost layer was added in dixB.
{#IFNDEF FOURTHORDER
    {^D& ixCRmin^D = ixCmin^D + kr(idims,^D)\}
    {^D& ixCRmax^D = ixCmax^D - kr(idims,^D)\}
}
{#IFDEF FOURTHORDER
    {^D& ixCRmin^D = ixCmin^D + 2 * kr(idims,^D)\}
    {^D& ixCRmax^D = ixCmax^D - 2 * kr(idims,^D)\}
}
}{#IFNDEF HALL
   {ixCR^L=ixC^L;}
}

   ! Calculate velocities for transport fluxes
   call getv(wCT,x,ixG^LL,ixCR^L,idims,v)
   
   do iw=1,nwflux
      call getfluxforhllc(wCT,x,ixG^LL,ixCR^L,iw,idims,fCT,transport)
      if (transport) fCT(ixCR^S,iw) = fCT(ixCR^S,iw) + v(ixCR^S) * wCT(ixCR^S,iw)
      ! Lax-Friedrich splitting:
      fp(ixCR^S,iw) = half * (fCT(ixCR^S,iw) + tvdlfeps * cmax_global * wCT(ixCR^S,iw))
      fm(ixCR^S,iw) = half * (fCT(ixCR^S,iw) - tvdlfeps * cmax_global * wCT(ixCR^S,iw))
   end do ! iw loop
  
   ! now do the reconstruction of fp and fm:
   call reconstructL(ixI^L,ix^L,idims,fp,fpL,dxdims)
   call reconstructR(ixI^L,ix^L,idims,fm,fmR,dxdims)

   do iw=1,nwflux
      if (slab) then
         fC(ix^S,iw,idims) = dxinv(idims) * (fpL(ix^S,iw) + fmR(ix^S,iw))
         wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
              (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
      else
         select case (idims)
         {case (^D)
            fC(ix^S,iw,^D)=-qdt*mygeo%surfaceC^D(ix^S) * (fpL(ix^S,iw) + fmR(ix^S,iw))
            wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
              (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/mygeo%dvolume(ixO^S)\}
         end select
      end if
   end do ! iw loop

end do !idims loop

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixI^L,ixO^L,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
                                   ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)


end subroutine fd
!=============================================================================
subroutine reconstructL(ixI^L,iL^L,idims,w,wLC,dxdims)
include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, iL^L, idims
double precision, intent(in)    :: w(ixI^S,1:nw), dxdims

double precision, intent(out)   :: wLC(ixI^S,1:nw) 

double precision                :: ldw(ixI^S), dwC(ixI^S)
integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, iw
character*79                    :: savetypelimiter
!-----------------------------------------------------------------------------
select case (typelimiter)
case ('mp5')
   call MP5limiterL(ixI^L,iL^L,idims,w,wLC)
case default 

   kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);

   wLC(kxC^S,1:nwflux) = w(kxC^S,1:nwflux)

   jxR^L=iL^L+kr(idims,^D);

   ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
   jxC^L=ixC^L+kr(idims,^D);

   do iw=1,nwflux
      dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)
      
      savetypelimiter=typelimiter
      if(savetypelimiter=='koren') typelimiter='korenL'
      if(savetypelimiter=='cada')  typelimiter='cadaL'
      if(savetypelimiter=='cada3') typelimiter='cada3L'
      call dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdims)
      typelimiter=savetypelimiter

      wLC(iL^S,iw)=wLC(iL^S,iw)+half*ldw(iL^S)
   end do
end select

end subroutine reconstructL
!=============================================================================
subroutine reconstructR(ixI^L,iL^L,idims,w,wRC,dxdims)
include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, iL^L, idims
double precision, intent(in)    :: w(ixI^S,1:nw), dxdims

double precision, intent(out)   :: wRC(ixI^S,1:nw) 

double precision                :: ldw(ixI^S), dwC(ixI^S)
integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, kxR^L, iw
character*79                    :: savetypelimiter
!-----------------------------------------------------------------------------
select case (typelimiter)
case ('mp5')
   call MP5limiterR(ixI^L,iL^L,idims,w,wRC)
case default 

   kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
   kxR^L=kxC^L+kr(idims,^D);

   wRC(kxC^S,1:nwflux)=w(kxR^S,1:nwflux)

   jxR^L=iL^L+kr(idims,^D);
   ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
   jxC^L=ixC^L+kr(idims,^D);

   do iw=1,nwflux
      dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)
      
      savetypelimiter=typelimiter
      if(savetypelimiter=='koren') typelimiter='korenR'
      if(savetypelimiter=='cada')  typelimiter='cadaR'
      if(savetypelimiter=='cada3') typelimiter='cada3R'
      call dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdims)
      typelimiter=savetypelimiter

      wRC(iL^S,iw)=wRC(iL^S,iw)-half*ldw(jxR^S)
   end do
end select

end subroutine reconstructR
!=============================================================================
