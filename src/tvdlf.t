INCLUDE:ppm.t
INCLUDE:mp5.t
!=============================================================================
subroutine hancock(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,wnew,dx^D,x)

! The non-conservative Hancock predictor for TVDLF

! on entry:
! input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBdixB

! one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2

! on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2


! FCT not implemented here
use mod_physics
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, idim^LIM
double precision, intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), wnew(ixI^S,1:nw)

double precision, dimension(ixI^S,1:nw) :: wLC, wRC
double precision, dimension(ixI^S) :: fLC, fRC
double precision, dimension(ixI^S) :: vLC, vRC
double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ix^L, hxO^L, ixtest^L
logical :: transport
!-----------------------------------------------------------------------------

! Expand limits in each idims direction in which fluxes are added
ix^L=ixO^L;
do idims= idim^LIM
   ix^L=ix^L^LADDkr(idims,^D);
end do
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Error in Hancock: Nonconforming input limits")

if (useprimitive) then  
   call phys_to_primitive(ixI^L,ixI^L,wCT,x)
endif 

^D&dxinv(^D)=-qdt/dx^D;
^D&dxdim(^D)=dx^D;
do idims= idim^LIM
   if (B0field) then
      select case (idims)
      {case (^D)
         myB0 => myB0_face^D\}
      end select
   end if

   ! Calculate w_j+g_j/2 and w_j-g_j/2
   ! First copy all variables, then upwind wLC and wRC.
   ! wLC is to the left of ixO, wRC is to the right of wCT.
   hxO^L=ixO^L-kr(idims,^D);

   wRC(hxO^S,1:nwflux)=wCT(ixO^S,1:nwflux)
   wLC(ixO^S,1:nwflux)=wCT(ixO^S,1:nwflux)

   call upwindLR(ixI^L,ixO^L,hxO^L,idims,wCT,wCT,wLC,wRC,x,.false.,dxdim(idims))

   if (nwaux>0.and.(.not.(useprimitive))) then
   !!if (nwaux>0) then
      call phys_get_aux(.true.,wLC,x,ixI^L,ixO^L,'hancock_wLC')
      call phys_get_aux(.true.,wRC,x,ixI^L,hxO^L,'hancock_wRC')
   end if

   ! Calculate vLC and vRC velocities
   call phys_get_v(wRC,x,ixI^L,hxO^L,idims,vRC)
   call phys_get_v(wLC,x,ixI^L,ixO^L,idims,vLC)

   ! Advect w(iw)
   do iw=1,nwflux
      ! Calculate the fLC and fRC fluxes
      call phys_get_flux(wRC,x,ixI^L,hxO^L,iw,idims,fRC,transport)
      call phys_get_flux(wLC,x,ixI^L,ixO^L,iw,idims,fLC,transport)
      if (transport) then
         fRC(hxO^S)=fRC(hxO^S)+vRC(hxO^S)*wRC(hxO^S,iw)
         fLC(ixO^S)=fLC(ixO^S)+vLC(ixO^S)*wLC(ixO^S,iw)
      end if

      ! Advect w(iw)
      if (slab) then
         wnew(ixO^S,iw)=wnew(ixO^S,iw)+dxinv(idims)* &
                          (fLC(ixO^S)-fRC(hxO^S))
      else
         select case (idims)
         {case (^D)
            wnew(ixO^S,iw)=wnew(ixO^S,iw)-qdt/mygeo%dvolume(ixO^S) &
                  *(mygeo%surfaceC^D(ixO^S)*fLC(ixO^S) &
                   -mygeo%surfaceC^D(hxO^S)*fRC(hxO^S))\}
         end select
      end if
   end do
end do ! next idims

if (useprimitive) then  
    call phys_to_conserved(ixI^L,ixI^L,wCT,x)
else
   if(nwaux>0) call phys_get_aux(.true.,wCT,x,ixI^L,ixI^L,'hancock_wCT')
endif

if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
                                   ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

end subroutine hancock
!============================================================================
subroutine upwindLR(ixI^L,ixL^L,ixR^L,idims,w,wCT,wLC,wRC,x,needprim,dxdim)

! Determine the upwinded wLC(ixL) and wRC(ixR) from w. 
! the wCT is only used when PPM is exploited.
use mod_physics
use mod_global_parameters

integer, intent(in) :: ixI^L, ixL^L, ixR^L, idims
logical, intent(in) :: needprim
double precision, intent(in) :: dxdim
double precision, dimension(ixI^S,1:nw) :: w, wCT
double precision, dimension(ixI^S,1:nw) :: wLC, wRC
double precision, dimension(ixI^S,1:ndim) :: x

integer :: jxR^L, ixC^L, jxC^L, iw, ixtest^L
double precision :: wLtmp(ixI^S,1:nw), wRtmp(ixI^S,1:nw)
double precision :: ldw(ixI^S), dwC(ixI^S)
logical, dimension(ixI^S) :: flagL, flagR
character*79 :: savetypelimiter
!-----------------------------------------------------------------------------

! Transform w,wL,wR to primitive variables
if (needprim.and.useprimitive) then
   call phys_to_primitive(ixI^L,ixI^L,w,x)
end if

if(typelimiter/='ppm' .and. typelimiter /= 'mp5')then
 jxR^L=ixR^L+kr(idims,^D);
 ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idims,^D);
 jxC^L=ixC^L+kr(idims,^D);

 do iw=1,nwflux
   if (loglimit(iw)) then
      w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
      wLC(ixL^S,iw)=dlog10(wLC(ixL^S,iw))
      wRC(ixR^S,iw)=dlog10(wRC(ixR^S,iw))
   end if

   dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

   savetypelimiter=typelimiter
   if(savetypelimiter=='koren') typelimiter='korenL'
   if(savetypelimiter=='cada')  typelimiter='cadaL'
   if(savetypelimiter=='cada3') typelimiter='cada3L'
   call dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdim)

   wLtmp(ixL^S,iw)=wLC(ixL^S,iw)+half*ldw(ixL^S)
   if(savetypelimiter=='koren')then
     typelimiter='korenR'
     call dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdim)
   endif
   if(savetypelimiter=='cada')then
     typelimiter='cadaR'
     call dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdim)
   endif
   if(savetypelimiter=='cada3')then
     typelimiter='cada3R'
     call dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdim)
   endif
   wRtmp(ixR^S,iw)=wRC(ixR^S,iw)-half*ldw(jxR^S)
   typelimiter=savetypelimiter

   if (loglimit(iw)) then
      w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
      wLtmp(ixL^S,iw)=10.0d0**wLtmp(ixL^S,iw)
      wRtmp(ixR^S,iw)=10.0d0**wRtmp(ixR^S,iw)
   end if
 end do

 call phys_check_w(useprimitive,ixI^L,ixL^L,wLtmp,flagL)
 call phys_check_w(useprimitive,ixI^L,ixR^L,wRtmp,flagR)

 do iw=1,nwflux
   where (flagL(ixL^S).and.flagR(ixR^S))
      wLC(ixL^S,iw)=wLtmp(ixL^S,iw)
      wRC(ixR^S,iw)=wRtmp(ixR^S,iw)
   end where

   if (loglimit(iw)) then
      where (.not.(flagL(ixL^S).and.flagR(ixR^S)))
         wLC(ixL^S,iw)=10.0d0**wLC(ixL^S,iw)
         wRC(ixR^S,iw)=10.0d0**wRC(ixR^S,iw)
      end where
   end if
 enddo
else if (typelimiter .eq. 'ppm') then
 call PPMlimiter(ixI^L,ixM^LL,idims,w,wCT,wLC,wRC)
else
 call MP5limiter(ixI^L,ixL^L,idims,w,wCT,wLC,wRC)
endif

! Transform w,wL,wR back to conservative variables
if (useprimitive) then
   if(needprim)then 
      call phys_to_conserved(ixI^L,ixI^L,w,x)
   endif
   call phys_to_conserved(ixI^L,ixL^L,wLC,x)
   call phys_to_conserved(ixI^L,ixR^L,wRC,x)
end if

end subroutine upwindLR
!============================================================================
subroutine dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdim)

! Limit the centered dwC differences within ixC for iw in direction idim.
! The limiter is chosen according to typelimiter.

! Note that this subroutine is called from upwindLR (hence from methods 
! like tvdlf, hancock, hll(c) etc) or directly from tvd.t,
! but also from the gradientS and divvectorS subroutines in geometry.t
! Accordingly, the typelimiter here corresponds to one of typelimiter1
! or one of typegradlimiter1.

! note: there is no iw dependence here...

use mod_global_parameters

integer, intent(in) :: ixI^L, ixC^L, iw, idims
double precision, intent(in) :: dxdim
double precision, intent(in) :: dwC(ixI^S)
double precision, intent(out) :: ldw(ixI^S)

double precision :: tmp(ixI^S)
integer :: ixO^L, hxO^L
double precision, parameter :: qsmall=1.d-12, qsmall2=2.d-12

! cada limiter parameter values
double precision, parameter :: cadalfa=0.5d0, cadbeta=2.0d0, cadgamma=1.6d0
! full third order cada limiter
double precision :: rdelinv
double precision :: ldwA(ixI^S),ldwB(ixI^S),tmpeta(ixI^S)
double precision, parameter :: cadepsilon=1.d-14, invcadepsilon=1.d14,cadradius=0.1d0
!-----------------------------------------------------------------------------

! Contract indices in idim for output.
ixOmin^D=ixCmin^D+kr(idims,^D); ixOmax^D=ixCmax^D;
hxO^L=ixO^L-kr(idims,^D);

! Store the sign of dwC in tmp
tmp(ixO^S)=sign(one,dwC(ixO^S))
rdelinv=one/(cadradius*dxdim)**2

select case (typelimiter)
case ('minmod')
   ! Minmod limiter eq(3.51e) and (eq.3.38e) with omega=1
   ldw(ixO^S)=tmp(ixO^S)* &
      max(zero,min(dabs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S)))
case ('woodward')
   ! Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
   ldw(ixO^S)=two*tmp(ixO^S)* &
      max(zero,min(dabs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S),&
                  tmp(ixO^S)*quarter*(dwC(hxO^S)+dwC(ixO^S))))
case ('mcbeta')
   ! Woodward and Collela limiter, with factor beta
   ldw(ixO^S)=tmp(ixO^S)* &
      max(zero,min(mcbeta*dabs(dwC(ixO^S)),mcbeta*tmp(ixO^S)*dwC(hxO^S),&
                  tmp(ixO^S)*half*(dwC(hxO^S)+dwC(ixO^S))))
case ('superbee')
   ! Roes superbee limiter (eq.3.51i)
   ldw(ixO^S)=tmp(ixO^S)* &
      max(zero,min(two*dabs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S)),&
              min(dabs(dwC(ixO^S)),two*tmp(ixO^S)*dwC(hxO^S)))
case ('vanleer')
  ! van Leer limiter (eq 3.51f), but a missing delta2=1.D-12 is added
  ldw(ixO^S)=two*max(dwC(hxO^S)*dwC(ixO^S),zero) &
             /(dwC(ixO^S)+dwC(hxO^S)+qsmall)
case ('albada')
  ! Albada limiter (eq.3.51g) with delta2=1D.-12
  ldw(ixO^S)=(dwC(hxO^S)*(dwC(ixO^S)**2+qsmall)&
             +dwC(ixO^S)*(dwC(hxO^S)**2+qsmall))&
             /(dwC(ixO^S)**2+dwC(hxO^S)**2+qsmall2)
case ('korenR')
   ! Barry Koren Right variant
   ldw(ixO^S)=tmp(ixO^S)* &
      max(zero,min(two*dabs(dwC(ixO^S)),two*tmp(ixO^S)*dwC(hxO^S),&
                   (two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third))
case ('korenL')
   ! Barry Koren Left variant
   ldw(ixO^S)=tmp(ixO^S)* &
      max(zero,min(two*dabs(dwC(ixO^S)),two*tmp(ixO^S)*dwC(hxO^S),&
                   (dwC(hxO^S)*tmp(ixO^S)+two*dabs(dwC(ixO^S)))*third))
case ('cadaR')
   ! Cada Right variant
   ldw(ixO^S)=tmp(ixO^S)* &
      max(zero,min((two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
                    max(-cadalfa*dabs(dwC(ixO^S)),                     &
                        min(cadbeta*dabs(dwC(ixO^S)),                  &
                            (two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
                            cadgamma*tmp(ixO^S)*dwC(hxO^S)))))
case ('cadaL')
   ! Cada Left variant
   ldw(ixO^S)=tmp(ixO^S)* &
      max(zero,min((two*dabs(dwC(ixO^S))+tmp(ixO^S)*dwC(hxO^S))*third, &
                    max(-cadalfa*tmp(ixO^S)*dwC(hxO^S),                     &
                        min(cadbeta*tmp(ixO^S)*dwC(hxO^S),                  &
                            (two*dabs(dwC(ixO^S))+tmp(ixO^S)*dwC(hxO^S))*third, &
                            cadgamma*dabs(dwC(ixO^S))))))
case ('cada3R')
   tmpeta(ixO^S)=(dwC(ixO^S)**2+dwC(hxO^S)**2)*rdelinv
   ldwA(ixO^S)=(two*dwC(hxO^S)+dwC(ixO^S))*third
   ldwB(ixO^S)=tmp(ixO^S)* &
      max(zero,min((two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
                    max(-cadalfa*dabs(dwC(ixO^S)),                     &
                        min(cadbeta*dabs(dwC(ixO^S)),                  &
                            (two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
                            cadgamma*tmp(ixO^S)*dwC(hxO^S)))))
   where(tmpeta(ixO^S)<=one-cadepsilon)
     ldw(ixO^S)=ldwA(ixO^S)
   elsewhere(tmpeta(ixO^S)>=one+cadepsilon)
     ldw(ixO^S)=ldwB(ixO^S)
   elsewhere
     tmp(ixO^S)=(tmpeta(ixO^S)-one)*invcadepsilon
     ldw(ixO^S)=half*( (one-tmp(ixO^S))*ldwA(ixO^S) &
                      +(one+tmp(ixO^S))*ldwB(ixO^S))
   endwhere
case ('cada3L')
   tmpeta(ixO^S)=(dwC(ixO^S)**2+dwC(hxO^S)**2)*rdelinv
   ldwA(ixO^S)=(two*dwC(ixO^S)+dwC(hxO^S))*third
   ldwB(ixO^S)=tmp(ixO^S)* &
      max(zero,min((two*dabs(dwC(ixO^S))+tmp(ixO^S)*dwC(hxO^S))*third, &
                    max(-cadalfa*tmp(ixO^S)*dwC(hxO^S),                     &
                        min(cadbeta*tmp(ixO^S)*dwC(hxO^S),                  &
                            (two*dabs(dwC(ixO^S))+tmp(ixO^S)*dwC(hxO^S))*third, &
                            cadgamma*dabs(dwC(ixO^S))))))
   where(tmpeta(ixO^S)<=one-cadepsilon)
     ldw(ixO^S)=ldwA(ixO^S)
   elsewhere(tmpeta(ixO^S)>=one+cadepsilon)
     ldw(ixO^S)=ldwB(ixO^S)
   elsewhere
     tmp(ixO^S)=(tmpeta(ixO^S)-one)*invcadepsilon
     ldw(ixO^S)=half*( (one-tmp(ixO^S))*ldwA(ixO^S) &
                      +(one+tmp(ixO^S))*ldwB(ixO^S))
   endwhere
case default
   write(*,*)'Unknown limiter:',typelimiter
   call mpistop("Error in dwLimiter: No such TVD limiter")
end select

end subroutine dwlimiter2
!=============================================================================
subroutine tvdlf(method,qdt,ixI^L,ixO^L,idim^LIM, &
                     qtC,wCT,qt,wnew,wold,fC,dx^D,x)

!> method=='tvdlf'  --> 2nd order TVD-Lax-Friedrich scheme.
!> method=='tvdlf1' --> 1st order TVD-Lax-Friedrich scheme.

use mod_physics
use mod_global_parameters

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx^D
integer, intent(in)                                  :: ixI^L, ixO^L, idim^LIM
double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
double precision, dimension(ixI^S,1:ndim)             ::  xi
double precision, dimension(ixI^S,1:nw)               :: wCT, wnew, wold
double precision, dimension(ixI^S,1:nwflux,1:ndim)        :: fC

double precision, dimension(ixI^S,1:nw) :: wLC, wRC
double precision, dimension(ixI^S)      :: fLC, fRC, vLC, vRC
double precision, dimension(ixI^S)      :: cmaxC, cmaxRC, cmaxLC
double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L, ixtest^L
logical :: transport, new_cmax, CmaxMeanState
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')

if (idimmax>idimmin .and. typelimited=='original' .and. &
   method/='tvdlf1')&
   call mpistop("Error in TVDMUSCLF: Unsplit dim. and original is limited")

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ix^L=ixO^L;
do idims= idim^LIM
   ix^L=ix^L^LADD2*kr(idims,^D);
end do
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Error in TVDLF: Nonconforming input limits")


if ((method=='tvdlf').and.useprimitive) then  
   ! second order methods with primitive limiting: 
   ! this call ensures wCT is primitive with updated auxiliaries
   call phys_to_primitive(ixI^L,ixI^L,wCT,x)
endif 


^D&dxinv(^D)=-qdt/dx^D;
^D&dxdim(^D)=dx^D;
do idims= idim^LIM
   if (B0field) then
      select case (idims)
      {case (^D)
         myB0 => myB0_face^D\}
      end select
   end if

   hxO^L=ixO^L-kr(idims,^D);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
!   ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;

{#IFDEF FCT
! Flux-interpolated constrained transport needs one more layer:
   ixCmax^D=ixOmax^D+1; ixCmin^D=hxOmin^D-1;
}{#IFNDEF FCT
   ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
}

   kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
   kxR^L=kxC^L+kr(idims,^D);
{#IFDEF HALL
   ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
   ! assuming one additional ghost layer (two for FOURTHORDER) was added in dixB.
{#IFNDEF FOURTHORDER
   {ixCR^L=ixC^L^LADD1;}
}
{#IFDEF FOURTHORDER
   {ixCR^L=ixC^L^LADD2;}
}
}{#IFNDEF HALL
   {ixCR^L=ixC^L;}
}
   wRC(kxC^S,1:nwflux)=wCT(kxR^S,1:nwflux)
   wLC(kxC^S,1:nwflux)=wCT(kxC^S,1:nwflux)

   ! Get interface positions:
   xi(kxC^S,1:ndim) = x(kxC^S,1:ndim)
   xi(kxC^S,idims) = half* ( x(kxR^S,idims)+x(kxC^S,idims) )
{#IFDEF STRETCHGRID
   if(idims==1) xi(kxC^S,1)=x(kxC^S,1)*(one+half*logG)
}

   ! for tvdlf (second order scheme): apply limiting
   if (method=='tvdlf') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wold,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wCT,wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wnew,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in TVDMUSCLF: no such base for limiter")
      end select
   end if

   ! For the high order Lax-Friedrich TVDLF scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   if (CmaxMeanState) then
      ! determine mean state and store in wLC
      wLC(ixC^S,1:nwflux)= &
           half*(wLC(ixC^S,1:nwflux)+wRC(ixC^S,1:nwflux))
      ! get auxilaries for mean state
      if (nwaux>0) then
         call phys_get_aux(.true.,wLC,xi,ixI^L,ixC^L,'tvdlf_cmaxmeanstate')
      end if

      call phys_get_cmax(wLC,xi,ixI^L,ixC^L,idims,cmaxC)

      ! We regain wLC for further use
      wLC(ixC^S,1:nwflux)=two*wLC(ixC^S,1:nwflux)-wRC(ixC^S,1:nwflux)
      if (nwaux>0) then
         call phys_get_aux(.true.,wLC,xi,ixI^L,ixC^L,'tvdlf_wLC_A')
      endif
      if (nwaux>0.and.(.not.(useprimitive).or.method=='tvdlf1')) then
         call phys_get_aux(.true.,wRC,xi,ixI^L,ixC^L,'tvdlf_wRC_A')
      end if
   else
      ! get auxilaries for L and R states
      if (nwaux>0.and.(.not.(useprimitive).or.method=='tvdlf1')) then
         !!if (nwaux>0) then
         call phys_get_aux(.true.,wLC,xi,ixI^L,ixC^L,'tvdlf_wLC')
         call phys_get_aux(.true.,wRC,xi,ixI^L,ixC^L,'tvdlf_wRC')
      end if

      call phys_get_cmax(wLC,xi,ixI^L,ixC^L,idims,cmaxLC)
      call phys_get_cmax(wRC,xi,ixI^L,ixC^L,idims,cmaxRC)
      ! now take the maximum of left and right states
      cmaxC(ixC^S)=max(cmaxRC(ixC^S),cmaxLC(ixC^S))
   end if
   
   ! Calculate velocities for transport fluxes
   call phys_get_v(wLC,xi,ixI^L,ixC^L,idims,vLC)
   call phys_get_v(wRC,xi,ixI^L,ixC^L,idims,vRC)


{#IFDEF GLM
! Solve the Riemann problem for the linear 2x2 system for normal
! B-field and GLM_Psi according to Dedner 2002:
call glmSolve(wLC,wRC,ixI^L,ixC^L,idims)
}


   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
      call phys_get_flux(wLC,xi,ixI^L,ixC^L,iw,idims,fLC,transport)
      call phys_get_flux(wRC,xi,ixI^L,ixC^L,iw,idims,fRC,transport)
      if (transport) then
         fLC(ixC^S)=fLC(ixC^S)+vLC(ixC^S)*wLC(ixC^S,iw)
         fRC(ixC^S)=fRC(ixC^S)+vRC(ixC^S)*wRC(ixC^S,iw)
      end if
      ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
      fLC(ixC^S)=half*(fLC(ixC^S)+fRC(ixC^S))

      ! Add TVDLF dissipation to the flux
      if ((.not.BnormLF) .and. (iw==b0_+idims{#IFDEF GLM .or.iw==psi_}) .and. b0_>0) then
         fRC(ixC^S)=0.d0
      else
         ! To save memory we use fRC to store -cmax*half*(w_R-w_L)
         fRC(ixC^S)=-tvdlfeps*cmaxC(ixC^S)*half*(wRC(ixC^S,iw)-wLC(ixC^S,iw))
      end if
      ! fLC contains physical+dissipative fluxes
      fLC(ixC^S)=fLC(ixC^S)+fRC(ixC^S)

      if (slab) then
         fC(ixC^S,iw,idims)=fLC(ixC^S)
      else
         select case (idims)
         {case (^D)
            fC(ixC^S,iw,^D)=mygeo%surfaceC^D(ixC^S)*fLC(ixC^S)\}
         end select
      end if

   end do ! Next iw

end do ! Next idims

{#IFDEF FCT
call fct_average(ixI^L,ixO^L,fC)
}

!Now update the state:
do idims= idim^LIM
   hxO^L=ixO^L-kr(idims,^D);
   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixI^S,iw,idims)=dxinv(idims)*fC(ixI^S,iw,idims)
         wnew(ixO^S,iw)=wnew(ixO^S,iw) &
              + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
      else
         select case (idims)
         {case (^D)
            fC(ixI^S,iw,^D)=-qdt*fC(ixI^S,iw,idims)
            wnew(ixO^S,iw)=wnew(ixO^S,iw) &
              + (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/mygeo%dvolume(ixO^S)\}
         end select
      end if

   end do ! Next iw

end do ! Next idims

if ((method=='tvdlf').and.useprimitive) then  
    call phys_to_conserved(ixI^L,ixI^L,wCT,x)
else
   if(nwaux>0) call phys_get_aux(.true.,wCT,x,ixI^L,ixI^L,'tvdlf_wCT')
endif

if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),ixI^L,ixO^L,1,nw,qtC,&
                wCT,qt,wnew,x,.false.)

end subroutine tvdlf
!=============================================================================
subroutine hll(method,qdt,ixI^L,ixO^L,idim^LIM, &
                     qtC,wCT,qt,wnew,wold,fC,dx^D,x)

!> method=='hll'  --> 2nd order HLL scheme.
!> method=='hll1' --> 1st order HLL scheme.
use mod_physics
use mod_global_parameters

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx^D
integer, intent(in)                                  :: ixI^L, ixO^L, idim^LIM
double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
double precision, dimension(ixI^S,1:ndim)             :: xi
double precision, dimension(ixI^S,1:nw)               :: wCT, wnew, wold
double precision, dimension(ixI^S,1:nwflux,1:ndim)  :: fC

double precision, dimension(ixI^S,1:nw) :: wLC, wRC
double precision, dimension(ixI^S)      :: fLC, fRC, vLC, vRC
double precision, dimension(ixI^S)      :: cmaxC, cmaxRC, cmaxLC
double precision, dimension(ixI^S)      :: cminC, cminRC, cminLC
double precision, dimension(1:ndim)     :: dxinv, dxdim
integer, dimension(ixI^S)               :: patchf
integer :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, jxC^L, kxC^L, kxR^L
logical :: transport, new_cmax, CmaxMeanState, logiB
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')
logiB=(BnormLF.and.b0_>0)

if (idimmax>idimmin .and. typelimited=='original' .and. &
   method/='hll1')&
   call mpistop("Error in hll: Unsplit dim. and original is limited")


! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ix^L=ixO^L;
do idims= idim^LIM
   ix^L=ix^L^LADD2*kr(idims,^D);
end do
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Error in hll : Nonconforming input limits")

if (method=='hll'.and.useprimitive) then  
   ! second order methods with primitive limiting: 
   ! this call ensures wCT is primitive with updated auxiliaries
   call phys_to_primitive(ixI^L,ixI^L,wCT,x)
endif 

^D&dxinv(^D)=-qdt/dx^D;
^D&dxdim(^D)=dx^D;
do idims= idim^LIM
   if (B0field) then
      select case (idims)
      {case (^D)
         myB0 => myB0_face^D\}
      end select
   end if

   hxO^L=ixO^L-kr(idims,^D);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
{#IFDEF FCT
! Flux-interpolated constrained transport needs one more layer:
   ixCmax^D=ixOmax^D+1; ixCmin^D=hxOmin^D-1;
}{#IFNDEF FCT
   ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
}

   ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2 
   jxC^L=ixC^L+kr(idims,^D);

   kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D); 
   kxR^L=kxC^L+kr(idims,^D);                        
{#IFDEF HALL
   ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
   ! assuming one additional ghost layer (two for FOURTHORDER) was added in dixB.
{#IFNDEF FOURTHORDER
   {ixCR^L=ixC^L^LADD1;}
}
{#IFDEF FOURTHORDER
   {ixCR^L=ixC^L^LADD2;}
}
}{#IFNDEF HALL
   {ixCR^L=ixC^L;}
}
   wRC(kxC^S,1:nwflux)=wCT(kxR^S,1:nwflux)
   wLC(kxC^S,1:nwflux)=wCT(kxC^S,1:nwflux)

   ! Get interface positions:
   xi(kxC^S,1:ndim) = x(kxC^S,1:ndim)
   xi(kxC^S,idims) = half* ( x(kxR^S,idims)+x(kxC^S,idims) )
{#IFDEF STRETCHGRID
   if(idims==1) xi(kxC^S,1)=x(kxC^S,1)*(one+half*logG)
}

   ! for hll (second order scheme): apply limiting
   if (method=='hll') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wold,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wCT,wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wnew,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in hll: no such base for limiter")
      end select
   end if

   ! For the high order hll scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   if (CmaxMeanState) then
         ! determine mean state and store in wLC
         wLC(ixC^S,1:nwflux)= &
               half*(wLC(ixC^S,1:nwflux)+wRC(ixC^S,1:nwflux))
         ! get auxilaries for mean state
         if (nwaux>0) then
            call phys_get_aux(.true.,wLC,xi,ixI^L,ixC^L,'hll_cmaxmeanstate')
         end if

         call phys_get_cmax(wLC,xi,ixI^L,ixC^L,idims,cmaxC,cminC)

         ! We regain wLC for further use
         wLC(ixC^S,1:nwflux)=two*wLC(ixC^S,1:nwflux)-wRC(ixC^S,1:nwflux)
         if (nwaux>0) then
            call phys_get_aux(.true.,wLC,xi,ixI^L,ixC^L,'hll_wLC_B')
         endif
         if (nwaux>0.and.(.not.(useprimitive).or.method=='hll1')) then
            call phys_get_aux(.true.,wRC,xi,ixI^L,ixC^L,'hll_wRC_B')
         end if
    else
         ! get auxilaries for L and R states
         if (nwaux>0.and.(.not.(useprimitive).or.method=='hll1')) then
         !!if (nwaux>0) then
            call phys_get_aux(.true.,wLC,xi,ixI^L,ixC^L,'hll_wLC')
            call phys_get_aux(.true.,wRC,xi,ixI^L,ixC^L,'hll_wRC')
         end if

         call phys_get_cmax(wLC,xi,ixI^L,ixC^L,idims,cmaxLC,cminLC)
         call phys_get_cmax(wRC,xi,ixI^L,ixC^L,idims,cmaxRC,cminRC)
         ! now take the maximum of left and right states 
         ! S.F. Davis, SIAM J. Sci. Statist. Comput. 1988, 9, 445
         cmaxC(ixC^S)=max(cmaxRC(ixC^S),cmaxLC(ixC^S))
         cminC(ixC^S)=min(cminRC(ixC^S),cminLC(ixC^S))
   end if 

   patchf(ixC^S) =  1
   where(cminC(ixC^S) >= zero)
        patchf(ixC^S) = -2
   elsewhere(cmaxC(ixC^S) <= zero)
        patchf(ixC^S) =  2
   endwhere

   ! Calculate velocities for transport fluxes
   if(any(patchf(ixC^S)/= 2).or.(logiB)) call phys_get_v(wLC,xi,ixI^L,ixC^L,idims,vLC)
   if(any(patchf(ixC^S)/=-2).or.(logiB)) call phys_get_v(wRC,xi,ixI^L,ixC^L,idims,vRC)

{#IFDEF GLM
! Solve the Riemann problem for the linear 2x2 system for normal
! B-field and GLM_Psi according to Dedner 2002:
call glmSolve(wLC,wRC,ixI^L,ixC^L,idims)
}


   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
     if(any(patchf(ixC^S)/= 2).or.(logiB.and.(iw==b0_+idims{#IFDEF GLM .or.iw==psi_}))) then 
        call phys_get_flux(wLC,xi,ixI^L,ixC^L,iw,idims,fLC,transport)
        if (transport) fLC(ixC^S)=fLC(ixC^S)+vLC(ixC^S)*wLC(ixC^S,iw)
     end if
     if(any(patchf(ixC^S)/=-2).or.(logiB.and.(iw==b0_+idims{#IFDEF GLM .or.iw==psi_}))) then 
        call phys_get_flux(wRC,xi,ixI^L,ixC^L,iw,idims,fRC,transport)
        if (transport) fRC(ixC^S)=fRC(ixC^S)+vRC(ixC^S)*wRC(ixC^S,iw)
     end if


     if (logiB.and.(iw==b0_+idims{#IFDEF GLM .or.iw==psi_})) then
        if (BnormLF) then
           ! flat B norm using tvdlf
           fLC(ixC^S)= half*((fLC(ixC^S)+fRC(ixC^S)) &
                         -tvdlfeps*max(cmaxC(ixC^S)&
                         ,dabs(cminC(ixC^S)))*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))
         else
           fLC(ixC^S)=zero
         endif
     else
       where(patchf(ixC^S)==1)
         ! Add hll dissipation to the flux
         fLC(ixC^S)=(cmaxC(ixC^S)*fLC(ixC^S)-cminC(ixC^S)*fRC(ixC^S) &
                   +tvdlfeps*cminC(ixC^S)*cmaxC(ixC^S)*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))&
                 /(cmaxC(ixC^S)-cminC(ixC^S))
       elsewhere(patchf(ixC^S)== 2)
         fLC(ixC^S)=fRC(ixC^S)
       elsewhere(patchf(ixC^S)==-2)
         fLC(ixC^S)=fLC(ixC^S)
       endwhere
     endif

     if (slab) then
         fC(ixC^S,iw,idims)=fLC(ixC^S)
     else
         select case (idims)
         {case (^D)
            fC(ixC^S,iw,^D)=mygeo%surfaceC^D(ixC^S)*fLC(ixC^S)\}
         end select
     end if

   end do ! Next iw
end do ! Next idims


{#IFDEF FCT
call fct_average(ixI^L,ixO^L,fC)
}
do idims= idim^LIM
   hxO^L=ixO^L-kr(idims,^D);
   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixI^S,iw,idims)=dxinv(idims)*fC(ixI^S,iw,idims)
         wnew(ixO^S,iw)=wnew(ixO^S,iw) &
              + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
      else
         select case (idims)
         {case (^D)
            fC(ixI^S,iw,^D)=-qdt*fC(ixI^S,iw,idims)
            wnew(ixO^S,iw)=wnew(ixO^S,iw) &
              + (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/mygeo%dvolume(ixO^S)\}
         end select
      end if

   end do ! Next iw

end do ! Next idims

if (method=='hll'.and.useprimitive) then  
   call phys_to_conserved(ixI^L,ixI^L,wCT,x)
else
   if(nwaux>0) call phys_get_aux(.true.,wCT,x,ixI^L,ixI^L,'hll_wCT')
endif

if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
                                   ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

end subroutine hll
!=============================================================================
subroutine hllc(method,qdt,ixI^L,ixO^L,idim^LIM, &
                     qtC,wCT,qt,wnew,wold,fC,dx^D,x)

! method=='hllc'  --> 2nd order HLLC scheme.
! method=='hllc1' --> 1st order HLLC scheme.
! method=='hllcd' --> 2nd order HLLC+tvdlf scheme.
! method=='hllcd1'--> 1st order HLLC+tvdlf scheme.

use mod_physics
use mod_global_parameters

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx^D
integer, intent(in)                                  :: ixI^L, ixO^L, idim^LIM
double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
double precision, dimension(ixI^S,1:ndim)             ::  xi
double precision, dimension(ixI^S,1:nw)               :: wCT, wnew, wold
double precision, dimension(ixI^S,1:nwflux,1:ndim)  :: fC

double precision, dimension(ixI^S,1:nw)            :: wLC, wRC
double precision, dimension(ixI^S)                 :: vLC, vRC
double precision, dimension(ixI^S)                 :: cmaxC,cminC

double precision, dimension(1:ndim)                :: dxinv, dxdim

integer, dimension(ixI^S)                          :: patchf
integer :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, jxC^L, kxC^L, kxR^L
logical :: transport, new_cmax, CmaxMeanState, logiB, firstordermethod

!=== specific to HLLC and HLLCD ===!
double precision, dimension(ixI^S,1:nwflux)     :: fLC, fRC
double precision, dimension(ixI^S,1:nwflux)     :: whll, Fhll, fCD
double precision, dimension(ixI^S)              :: lambdaCD 
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')
logiB=(BnormLF.and.b0_>0)
logiB=.false.

if (idimmax>idimmin .and. typelimited=='original' .and. &
   method/='hllc1' .and. method/='hllcd1')&
   call mpistop("Error in hllc: Unsplit dim. and original is limited")



! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ix^L=ixO^L;
do idims= idim^LIM
   ix^L=ix^L^LADD2*kr(idims,^D);
end do
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Error in hllc : Nonconforming input limits")

if ((method=='hllc'.or.method=='hllcd').and.useprimitive) then
   call phys_to_primitive(ixI^L,ixI^L,wCT,x)
endif 
firstordermethod=(method=='hllc1'.or.method=='hllcd1')

^D&dxinv(^D)=-qdt/dx^D;
^D&dxdim(^D)=dx^D;
do idims= idim^LIM
   if (B0field) then
      select case (idims)
      {case (^D)
         myB0 => myB0_face^D\}
      end select
   end if

   hxO^L=ixO^L-kr(idims,^D);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
{#IFDEF FCT
! Flux-interpolated constrained transport needs one more layer:
   ixCmax^D=ixOmax^D+1; ixCmin^D=hxOmin^D-1;
}{#IFNDEF FCT
   ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
}

   ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2 
   jxC^L=ixC^L+kr(idims,^D);

   ! enlarged for ppm purposes
   kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D); 
   kxR^L=kxC^L+kr(idims,^D);                         
{#IFDEF HALL
   ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
   ! assuming one additional ghost layer (two for FOURTHORDER) was added in dixB.
{#IFNDEF FOURTHORDER
   {ixCR^L=ixC^L^LADD1;}
}
{#IFDEF FOURTHORDER
   {ixCR^L=ixC^L^LADD2;}
}
}{#IFNDEF HALL
   {ixCR^L=ixC^L;}
}
   wRC(kxC^S,1:nwflux)=wCT(kxR^S,1:nwflux)
   wLC(kxC^S,1:nwflux)=wCT(kxC^S,1:nwflux)

   ! Get interface positions:
   xi(kxC^S,1:ndim) = x(kxC^S,1:ndim)
   xi(kxC^S,idims) = half* ( x(kxR^S,idims)+x(kxC^S,idims) )
{#IFDEF STRETCHGRID
   if(idims==1) xi(kxC^S,1)=x(kxC^S,1)*(one+half*logG)
}

   ! for hllc and hllcd (second order schemes): apply limiting
   if (method=='hllc'.or.method=='hllcd') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wold,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wCT,wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wnew,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in hllc: no such base for limiter")
      end select
   end if

   ! For the high order hllc scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   if (CmaxMeanState) then
         ! determine mean state and store in wLC
         wLC(ixC^S,1:nwflux)= &
               half*(wLC(ixC^S,1:nwflux)+wRC(ixC^S,1:nwflux))
         ! get auxilaries for mean state
         if (nwaux>0) then
            call phys_get_aux(.true.,wLC,xi,ixI^L,ixC^L,'hllc_cmaxmeanstate')
         end if

         call phys_get_cmax(wLC,xi,ixI^L,ixC^L,idims,cmaxC,cminC)

         ! We regain wLC for further use
         wLC(ixC^S,1:nwflux)=two*wLC(ixC^S,1:nwflux)-wRC(ixC^S,1:nwflux)
         if (nwaux>0) then
            call phys_get_aux(.true.,wLC,xi,ixI^L,ixC^L,'hllc_wLC_B')
         end if
         if (nwaux>0.and.(.not.(useprimitive).or.firstordermethod)) then
            call phys_get_aux(.true.,wRC,xi,ixI^L,ixC^L,'hllc_wRC_B')
         end if
   else
         ! get auxilaries for L and R states
         if (nwaux>0.and.(.not.(useprimitive).or.firstordermethod)) then
         !!if (nwaux>0) then
            call phys_get_aux(.true.,wLC,xi,ixI^L,ixC^L,'hllc_wLC')
            call phys_get_aux(.true.,wRC,xi,ixI^L,ixC^L,'hllc_wRC')
         end if

         ! to save memory, use cmaxC and lambdaCD for cmacRC and cmaxLC respectively
         ! to save memory, use vLC   and vRC      for cminRC and cminLC respectively
         call phys_get_cmax(wLC,xi,ixI^L,ixC^L,idims,cmaxC,vLC)
         call phys_get_cmax(wRC,xi,ixI^L,ixC^L,idims,lambdaCD,vRC)
         ! now take the maximum of left and right states
         cmaxC(ixC^S)=max(lambdaCD(ixC^S),cmaxC(ixC^S))
         cminC(ixC^S)=min(vRC(ixC^S),vLC(ixC^S))
   end if

   patchf(ixC^S) =  1
   where(cminC(ixC^S) >= zero)
        patchf(ixC^S) = -2
   elsewhere(cmaxC(ixC^S) <= zero)
        patchf(ixC^S) =  2
   endwhere

   ! Calculate velocities for transport fluxes
   if(any(patchf(ixC^S)/= 2).or.(logiB)) call phys_get_v(wLC,xi,ixI^L,ixC^L,idims,vLC)
   if(any(patchf(ixC^S)/=-2).or.(logiB)) call phys_get_v(wRC,xi,ixI^L,ixC^L,idims,vRC)

{#IFDEF GLM
! Solve the Riemann problem for the linear 2x2 system for normal
! B-field and GLM_Psi according to Dedner 2002:
call glmSolve(wLC,wRC,ixI^L,ixC^L,idims)
}

   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
     if(any(patchf(ixC^S)/= 2).or.(logiB.and.(iw==b0_+idims{#IFDEF GLM .or.iw==psi_}))) then
        call phys_get_flux(wLC,xi,ixI^L,ixC^L,iw,idims,fLC(ixI^S,iw),transport)
        if (transport)  fLC(ixC^S,iw)=fLC(ixC^S,iw)+vLC(ixC^S)*wLC(ixC^S,iw)
     end if
     if(any(patchf(ixC^S)/=-2).or.(logiB.and.(iw==b0_+idims{#IFDEF GLM .or.iw==psi_}))) then
        call phys_get_flux(wRC,xi,ixI^L,ixC^L,iw,idims,fRC(ixI^S,iw),transport)
        if (transport)   fRC(ixC^S,iw)=fRC(ixC^S,iw)+vRC(ixC^S)*wRC(ixC^S,iw)
     end if
   end do

   ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4 
   if(method=='hllcd' .or. method=='hllcd1') &
     call phys_diffuse_hllcd(ixI^L,ixC^L,idims,wLC,wRC,fLC,fRC,patchf)

   !---- calculate speed lambda at CD ----!
   if(any(patchf(ixC^S)==1)) &
     call phys_get_lCD(wLC,wRC,fLC,fRC,cminC,cmaxC,idims,ixI^L,ixC^L, &
            whll,Fhll,lambdaCD,patchf)

   ! now patchf may be -1 or 1 due to phys_get_lCD 
   if(any(abs(patchf(ixC^S))== 1))then
      !======== flux at intermediate state ========!
      call phys_get_wCD(wLC,wRC,whll,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaCD,&
                  cminC,cmaxC,ixI^L,ixC^L,idims,fCD)
   endif ! Calculate the CD flux

   do iw=1,nwflux
     if (logiB.and.(iw==b0_+idims{#IFDEF GLM .or.iw==psi_})) then
        if (BnormLF) then
           ! flat B norm using tvdlf
           fLC(ixC^S,iw) = half*((fLC(ixC^S,iw)+fRC(ixC^S,iw)) &
                           -tvdlfeps*max(cmaxC(ixC^S)&
                           ,dabs(cminC(ixC^S)))*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))
        else
           fLC(ixC^S,iw)=zero
        end if
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
          -tvdlfeps*max(cmaxC(ixC^S),dabs(cminC(ixC^S)))*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))
       endwhere
     endif

     if (slab) then
         fC(ixC^S,iw,idims)=fLC(ixC^S,iw)
     else
         select case (idims)
         {case (^D)
            fC(ixC^S,iw,^D)=mygeo%surfaceC^D(ixC^S)*fLC(ixC^S,iw)\}
         end select
      end if

   end do ! Next iw

end do ! Next idims


{#IFDEF FCT
call fct_average(ixI^L,ixO^L,fC)
}


!Now update the state
do idims= idim^LIM
   hxO^L=ixO^L-kr(idims,^D);
   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixI^S,iw,idims)=dxinv(idims)*fC(ixI^S,iw,idims)
         wnew(ixO^S,iw)=wnew(ixO^S,iw) &
              + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
      else
         select case (idims)
         {case (^D)
            fC(ixI^S,iw,^D)=-qdt*fC(ixI^S,iw,idims)
            wnew(ixO^S,iw)=wnew(ixO^S,iw) &
              + (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/mygeo%dvolume(ixO^S)\}
         end select
      end if

   end do ! Next iw

end do ! Next idims

if ((method=='hllc'.or.method=='hllcd').and.useprimitive) then  
   call phys_to_conserved(ixI^L,ixI^L,wCT,x)
else
   if(nwaux>0) call phys_get_aux(.true.,wCT,x,ixI^L,ixI^L,'hllc_wCT')
endif

if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
                                   ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

end subroutine hllc
!=============================================================================
subroutine tvdmusclf(method,qdt,ixI^L,ixO^L,idim^LIM, &
                     qtC,wCT,qt,wnew,wold,fC,dx^D,x)

! method=='tvdmu'  --> 2nd order (may be 3rd order in 1D) TVD-MUSCL scheme.
! method=='tvdmu1' --> 1st order TVD-MUSCL scheme (upwind per charact. var.)
! FCT not implemented here.
use mod_physics
use mod_global_parameters

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx^D
integer, intent(in)                                  :: ixI^L, ixO^L, idim^LIM
double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
double precision, dimension(ixI^S,1:ndim)             ::  xi
double precision, dimension(ixI^S,1:nw)               :: wCT, wnew, wold
double precision, dimension(ixI^S,1:nwflux,1:ndim)        :: fC

double precision, dimension(ixI^S,1:nw) :: wLC, wRC
double precision, dimension(ixI^S)      :: fLC, fRC, vLC, vRC
double precision, dimension(ixI^S)      :: cmaxC, cmaxRC, cmaxLC
double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L, ixtest^L
logical :: transport, new_cmax, CmaxMeanState
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')

if (idimmax>idimmin .and. typelimited=='original' .and. &
   method/='tvdmu1')&
   call mpistop("Error in TVDMUSCLF: Unsplit dim. and original is limited")

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ix^L=ixO^L;
do idims= idim^LIM
   ix^L=ix^L^LADD2*kr(idims,^D);
end do
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Error in TVDMUSCLF: Nonconforming input limits")


if ((method=='tvdmu').and.useprimitive) then  
   ! second order methods with primitive limiting: 
   ! this call ensures wCT is primitive with updated auxiliaries
   call phys_to_primitive(ixI^L,ixI^L,wCT,x)
endif 


^D&dxinv(^D)=-qdt/dx^D;
^D&dxdim(^D)=dx^D;
do idims= idim^LIM
   if (B0field) then
      select case (idims)
      {case (^D)
         myB0 => myB0_face^D\}
      end select
   end if

   hxO^L=ixO^L-kr(idims,^D);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
   ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;

   kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
   kxR^L=kxC^L+kr(idims,^D);
{#IFDEF HALL
   ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
   ! assuming one additional ghost layer (two for FOURTHORDER) was added in dixB.
{#IFNDEF FOURTHORDER
   {ixCR^L=ixC^L^LADD1;}
}
{#IFDEF FOURTHORDER
   {ixCR^L=ixC^L^LADD2;}
}
}{#IFNDEF HALL
   {ixCR^L=ixC^L;}
}
   wRC(kxC^S,1:nwflux)=wCT(kxR^S,1:nwflux)
   wLC(kxC^S,1:nwflux)=wCT(kxC^S,1:nwflux)

   ! Get interface positions:
   xi(kxC^S,1:ndim) = x(kxC^S,1:ndim)
   xi(kxC^S,idims) = half* ( x(kxR^S,idims)+x(kxC^S,idims) )
{#IFDEF STRETCHGRID
   if(idims==1) xi(kxC^S,1)=x(kxC^S,1)*(one+half*logG)
}

   ! for tvdlf and tvdmu (second order schemes): apply limiting
   if (method=='tvdmu') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wold,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wCT,wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wnew,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in TVDMUSCLF: no such base for limiter")
      end select
   end if

   ! handle all other methods than tvdlf, namely tvdmu and tvdmu1 here
   if (nwaux>0.and.(.not.(useprimitive).or.method=='tvdmu1')) then
      !!if (nwaux>0) then
      call phys_get_aux(.true.,wLC,xi,ixI^L,ixC^L,'tvdlf_wLC_B')
      call phys_get_aux(.true.,wRC,xi,ixI^L,ixC^L,'tvdlf_wRC_B')
   end if
   

   ! Calculate velocities for transport fluxes
   call phys_get_v(wLC,xi,ixI^L,ixC^L,idims,vLC)
   call phys_get_v(wRC,xi,ixI^L,ixC^L,idims,vRC)


{#IFDEF GLM
! Solve the Riemann problem for the linear 2x2 system for normal
! B-field and GLM_Psi according to Dedner 2002:
call glmSolve(wLC,wRC,ixI^L,ixC^L,idims)
}

   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
      call phys_get_flux(wLC,xi,ixI^L,ixC^L,iw,idims,fLC,transport)
      call phys_get_flux(wRC,xi,ixI^L,ixC^L,iw,idims,fRC,transport)
      if (transport) then
         fLC(ixC^S)=fLC(ixC^S)+vLC(ixC^S)*wLC(ixC^S,iw)
         fRC(ixC^S)=fRC(ixC^S)+vRC(ixC^S)*wRC(ixC^S,iw)
      end if
      ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
      fLC(ixC^S)=half*(fLC(ixC^S)+fRC(ixC^S))

      if (slab) then
         fC(ixC^S,iw,idims)=dxinv(idims)*fLC(ixC^S)
         wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
              (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
      else
         select case (idims)
         {case (^D)
            fC(ixC^S,iw,^D)=-qdt*mygeo%surfaceC^D(ixC^S)*fLC(ixC^S)
            wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
              (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/mygeo%dvolume(ixO^S)\}
         end select
      end if

   end do ! Next iw

   ! For the MUSCL scheme apply the characteristic based limiter
   if (method=='tvdmu'.or.method=='tvdmu1') &
      call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,wLC,wRC,wnew,x,fC,dx^D)

end do ! Next idims

if ((method=='tvdmu').and.useprimitive) then  
    call phys_to_conserved(ixI^L,ixI^L,wCT,x)
else
   if(nwaux>0) call phys_get_aux(.true.,wCT,x,ixI^L,ixI^L,'tvdlf_wCT')
endif

if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),ixI^L,ixO^L,1,nw,qtC,&
                wCT,qt,wnew,x,.false.)

end subroutine tvdmusclf
!=============================================================================
{#IFDEF GLM
subroutine glmSolve(wLC,wRC,ixI^L,ixO^L,idir)
use mod_global_parameters
double precision, dimension(ixI^S,1:nw), intent(inout) :: wLC, wRC
integer, intent(in)                                  :: ixI^L, ixO^L, idir
double precision, dimension(ixI^S)                   :: dB, dPsi
!-----------------------------------------------------------------------------

! This implements eq. (42) in Dedner et al. 2002 JcP 175
! Gives the Riemann solution on the interface 
! for the normal B component and Psi in the GLM-MHD system.

! 23/04/2013 Oliver Porth

!BL(ixO^S) = wLC(ixO^S,b0_+idir)
!BR(ixO^S) = wRC(ixO^S,b0_+idir)

!PsiL(ixO^S) = wLC(ixO^S,psi_)
!PsiR(ixO^S) = wRC(ixO^S,psi_)

dB(ixO^S)   = wRC(ixO^S,b0_+idir) - wLC(ixO^S,b0_+idir)
dPsi(ixO^S) = wRC(ixO^S,psi_) - wLC(ixO^S,psi_)

wLC(ixO^S,b0_+idir)   = half * (wRC(ixO^S,b0_+idir) + wLC(ixO^S,b0_+idir)) - half/cmax_global * dPsi(ixO^S)
wLC(ixO^S,psi_)       = half * (wRC(ixO^S,psi_) + wLC(ixO^S,psi_))         - half*cmax_global * dB(ixO^S)

wRC(ixO^S,b0_+idir) = wLC(ixO^S,b0_+idir)
wRC(ixO^S,psi_) = wLC(ixO^S,psi_)

end subroutine glmSolve
!=============================================================================
}
