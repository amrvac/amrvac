!##############################################################################
! module vactvd
! Subroutines for TVD-MUSCL schemes
!=============================================================================
subroutine tvdlimit(method,qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,idimmin,&
   idimmax,w,qt,wnew,fC,dx1,x)

include 'amrvacdef.f'

character(len=*), intent(in) :: method
double precision, intent(in) :: qdt, qt, dx1
integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,nw) :: w, wnew
double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
double precision :: fC(ixImin1:ixImax1,1:nwflux,1:ndim)

integer :: idims, ixICmin1,ixICmax1, jxICmin1,jxICmax1
double precision, dimension(ixGlo1:ixGhi1,nw) :: wR, wL
!-----------------------------------------------------------------------------

do idims= idimmin,idimmax
   ixICmax1=ixOmax1+kr(idims,1); ixICmin1=ixOmin1-2*kr(idims,1);
   wL(ixICmin1:ixICmax1,1:nw)=w(ixICmin1:ixICmax1,1:nw)
   jxICmin1=ixICmin1+kr(idims,1);jxICmax1=ixICmax1+kr(idims,1);
   wR(ixICmin1:ixICmax1,1:nw)=w(jxICmin1:jxICmax1,1:nw)
   call tvdlimit2(method,qdt,ixImin1,ixImax1,ixICmin1,ixICmax1,ixOmin1,&
      ixOmax1,idims,wL,wR,wnew,x,fC,dx1)
end do

end subroutine tvdlimit
!=============================================================================
subroutine tvdlimit2(method,qdt,ixImin1,ixImax1,ixICmin1,ixICmax1,ixOmin1,&
   ixOmax1,idims,wL,wR,wnew,x,fC,dx1)

! Limit the flow variables in wnew according to typetvd. 
! wroeC is based on wL and wR.
! If method=='tvd' an extra adtdx**2*jumpC is added to phiC for 2nd order
! accuracy in time.

include 'amrvacdef.f'

character(len=*), intent(in) :: method
double precision, intent(in) :: qdt, dx1
integer, intent(in) :: ixImin1,ixImax1, ixICmin1,ixICmax1, ixOmin1,ixOmax1,&
    idims
double precision, dimension(ixGlo1:ixGhi1,nw) :: wL, wR
double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
double precision :: wnew(ixImin1:ixImax1,1:nw)
double precision :: fC(ixImin1:ixImax1,1:nwflux,1:ndim)

double precision:: workroe(ixGlo1:ixGhi1,1:nworkroe)
double precision, dimension(ixGlo1:ixGhi1,nw) :: wroeC
double precision, dimension(ixGlo1:ixGhi1) :: phiC, rphiC, jumpC, adtdxC,&
    smallaC
double precision :: dxinv(1:ndim),dxdim
integer :: hxOmin1,hxOmax1, ixCmin1,ixCmax1, jxCmin1,jxCmax1, jxICmin1,&
   jxICmax1, iw, il
!-----------------------------------------------------------------------------

hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
ixCmax1=ixOmax1; ixCmin1=hxOmin1; 

jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);
jxICmin1=ixICmin1+kr(idims,1);jxICmax1=ixICmax1+kr(idims,1);

call average(wL,wR,x,ixICmin1,ixICmax1,idims,wroeC,workroe)

dxinv(1)=qdt/dx1;

! A loop on characteristic variables to calculate the dissipative flux phiC.
do il=1,nwflux
   !Calculate the jump in the il-th characteristic variable: L(wroe)*dw
   call geteigenjump(wL,wR,wroeC,x,ixICmin1,ixICmax1,il,idims,smallaC,adtdxC,&
      jumpC,workroe)

   ! Normalize the eigenvalue "a" (and its limit "smalla" if needed):
   adtdxC(ixICmin1:ixICmax1)=adtdxC(ixICmin1:ixICmax1)*dxinv(idims)
   if (typeentropy(il)=='harten' .or. typeentropy(il)=='powell')smallaC&
      (ixICmin1:ixICmax1)=smallaC(ixICmin1:ixICmax1)*dxinv(idims)

   ! Calculate the flux limiter function phi
   dxdim=qdt/dxinv(idims)
   call getphi(method,jumpC,adtdxC,smallaC,ixICmin1,ixICmax1,ixCmin1,ixCmax1,&
      il,idims,phiC,dxdim)

   if (.not.slab) call mpistop("geometry need still be implemented in tvd")

   !Add R(iw,il)*phiC(il) to each variable iw in wnew
   do iw=1,nwflux
      call rtimes(phiC,wroeC,ixCmin1,ixCmax1,iw,il,idims,rphiC,workroe)

      rphiC(ixCmin1:ixCmax1)=rphiC(ixCmin1:ixCmax1)*half

      fC(ixCmin1:ixCmax1,iw,idims)=fC(ixCmin1:ixCmax1,iw,idims)&
         +rphiC(ixCmin1:ixCmax1)

      wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw)+rphiC&
         (ixOmin1:ixOmax1)-rphiC(hxOmin1:hxOmax1)
   end do  !iw
end do     !il

end subroutine tvdlimit2
!=============================================================================
subroutine getphi(method,jumpC,adtdxC,smallaC,ixICmin1,ixICmax1,ixCmin1,&
   ixCmax1,il,idims,phiC,dxdim)

! Calculate the dissipative flux from jumpC=L*dw and adtdx=eigenvalue*dt/dx.
! Add Lax-Wendroff type correction if method=='tvd'.
! Limit according to method and typetvd.

include 'amrvacdef.f'

character(len=*), intent(in) :: method
integer, intent(in) :: ixICmin1,ixICmax1, ixCmin1,ixCmax1, il, idims
double precision, intent(in) :: dxdim
double precision, dimension(ixGlo1:ixGhi1) :: jumpC, adtdxC, smallaC, phiC

double precision, dimension(ixGlo1:ixGhi1) :: ljumpC, tmp
integer :: jxCmin1,jxCmax1, ixmin1,ixmax1, hxmin1,hxmax1
!-----------------------------------------------------------------------------

if(method=='tvdmu'.or.method=='tvdmu1')then
   ! In the MUSCL scheme phi=|a|*jump, apply entropy fix to it
   if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
      phiC(ixCmin1:ixCmax1)=abs(adtdxC(ixCmin1:ixCmax1))*jumpC&
         (ixCmin1:ixCmax1)
   else
      where(abs(adtdxC(ixCmin1:ixCmax1))>=smallaC(ixCmin1:ixCmax1))
         phiC(ixCmin1:ixCmax1)=abs(adtdxC(ixCmin1:ixCmax1))&
            *jumpC(ixCmin1:ixCmax1)
      elsewhere
         phiC(ixCmin1:ixCmax1)=half*(smallaC(ixCmin1:ixCmax1)&
            +adtdxC(ixCmin1:ixCmax1)**2/smallaC(ixCmin1:ixCmax1))&
            *jumpC(ixCmin1:ixCmax1)
      endwhere
   endif
   ! That's all for the MUSCL scheme
   return
endif

if(method=='tvd')then
   !Entropy fix to |a|-a**2
   if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
     phiC(ixICmin1:ixICmax1)=abs(adtdxC(ixICmin1:ixICmax1))&
        -adtdxC(ixICmin1:ixICmax1)**2
   else
      where(abs(adtdxC(ixICmin1:ixICmax1))>=smallaC(ixICmin1:ixICmax1))
         phiC(ixICmin1:ixICmax1)=abs(adtdxC(ixICmin1:ixICmax1))&
            -adtdxC(ixICmin1:ixICmax1)**2
      elsewhere
         phiC(ixICmin1:ixICmax1)=half*smallaC(ixICmin1:ixICmax1)&
            +(half/smallaC(ixICmin1:ixICmax1)-one)*adtdxC(ixICmin1:ixICmax1)&
            **2
      endwhere
   endif
else
   !Entropy fix to |a|
   if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
      phiC(ixICmin1:ixICmax1)=abs(adtdxC(ixICmin1:ixICmax1))
   else
      where(abs(adtdxC(ixICmin1:ixICmax1))>=smallaC(ixICmin1:ixICmax1))
         phiC(ixICmin1:ixICmax1)=abs(adtdxC(ixICmin1:ixICmax1))
      elsewhere
         phiC(ixICmin1:ixICmax1)=half*smallaC(ixICmin1:ixICmax1)&
            +half/smallaC(ixICmin1:ixICmax1)*adtdxC(ixICmin1:ixICmax1)**2
      endwhere
   endif
endif

jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);
hxmin1=ixICmin1; hxmax1=ixICmax1-kr(idims,1);
ixmin1=hxmin1+kr(idims,1);ixmax1=hxmax1+kr(idims,1);

select case(typetvd)
case('roe')
   call dwlimiter2(jumpC,ixICmin1,ixICmax1,il,idims,ljumpC,dxdim)
   where(adtdxC(ixCmin1:ixCmax1)<=0)
      phiC(ixCmin1:ixCmax1)=phiC(ixCmin1:ixCmax1)*(jumpC(ixCmin1:ixCmax1)&
         -ljumpC(jxCmin1:jxCmax1))
   elsewhere
      phiC(ixCmin1:ixCmax1)=phiC(ixCmin1:ixCmax1)*(jumpC(ixCmin1:ixCmax1)&
         -ljumpC(ixCmin1:ixCmax1))
   end where
   !extra (a*lambda)**2*delta
   if(method=='tvd')phiC(ixCmin1:ixCmax1)=phiC(ixCmin1:ixCmax1)&
      +adtdxC(ixCmin1:ixCmax1)**2*jumpC(ixCmin1:ixCmax1)
case('sweby')
   !Sweby eqs.4.11-4.15, but no 0.5 ?!
   phiC(ixICmin1:ixICmax1)=phiC(ixICmin1:ixICmax1)*jumpC(ixICmin1:ixICmax1)
   call dwlimiter2(phiC,ixICmin1,ixICmax1,il,idims,ljumpC,dxdim)
   where(adtdxC(ixCmin1:ixCmax1)<=0)
      phiC(ixCmin1:ixCmax1)=phiC(ixCmin1:ixCmax1)-ljumpC(jxCmin1:jxCmax1)
   elsewhere
      phiC(ixCmin1:ixCmax1)=phiC(ixCmin1:ixCmax1)-ljumpC(ixCmin1:ixCmax1)
   end where
   !extra (a*lambda)**2*delta
   if(method=='tvd')phiC(ixCmin1:ixCmax1)=phiC(ixCmin1:ixCmax1)&
      +adtdxC(ixCmin1:ixCmax1)**2*jumpC(ixCmin1:ixCmax1)
case('yee')
   !eq.3.51 with correction
   call dwlimiter2(jumpC,ixICmin1,ixICmax1,il,idims,ljumpC,dxdim)

   !Use phiC as 0.5*(|nu|-nu**2) eq.3.45e for tvd otherwise 0.5*|nu|
   phiC(ixCmin1:ixCmax1)=half*phiC(ixCmin1:ixCmax1)
   !gamma*lambda eq.3.51d, use tmp to store agdtdxC
   where(abs(jumpC(ixCmin1:ixCmax1))>smalldouble)
      tmp(ixCmin1:ixCmax1)=adtdxC(ixCmin1:ixCmax1)+phiC(ixCmin1:ixCmax1)&
         *(ljumpC(jxCmin1:jxCmax1)-ljumpC(ixCmin1:ixCmax1))&
         /jumpC(ixCmin1:ixCmax1)
   elsewhere
      tmp(ixCmin1:ixCmax1)=adtdxC(ixCmin1:ixCmax1)
   end where

   !eq.3.51a
   if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
      phiC(ixCmin1:ixCmax1)=-phiC(ixCmin1:ixCmax1)*(ljumpC(jxCmin1:jxCmax1)&
         +ljumpC(ixCmin1:ixCmax1))+abs(tmp(ixCmin1:ixCmax1))&
         *jumpC(ixCmin1:ixCmax1)
   else
      where(abs(tmp(ixCmin1:ixCmax1))>=smallaC(ixCmin1:ixCmax1))
         phiC(ixCmin1:ixCmax1)=-phiC(ixCmin1:ixCmax1)*(ljumpC&
            (jxCmin1:jxCmax1)+ljumpC(ixCmin1:ixCmax1))+abs(tmp&
            (ixCmin1:ixCmax1))*jumpC(ixCmin1:ixCmax1)
      elsewhere
         phiC(ixCmin1:ixCmax1)=-phiC(ixCmin1:ixCmax1)*(ljumpC&
            (jxCmin1:jxCmax1)+ljumpC(ixCmin1:ixCmax1))+(half&
            *smallaC(ixCmin1:ixCmax1)+half/smallaC(ixCmin1:ixCmax1)&
            *tmp(ixCmin1:ixCmax1)**2)*jumpC(ixCmin1:ixCmax1)
      endwhere
   endif
case('harten')
   !See Ryu, section 2.3
   !Use phiC as 0.5*(|nu|-nu**2)*jumpC eq.3.45b,e
   phiC(ixICmin1:ixICmax1)=half*phiC(ixICmin1:ixICmax1)*jumpC&
      (ixICmin1:ixICmax1)
   call dwlimiter2(phiC,ixICmin1,ixICmax1,il,idims,ljumpC,dxdim)

   !gamma*lambda eq.3.45d, use tmp as agdtdxC
   where(abs(jumpC(ixCmin1:ixCmax1))>smalldouble)
      tmp(ixCmin1:ixCmax1)=adtdxC(ixCmin1:ixCmax1)+(ljumpC(jxCmin1:jxCmax1)&
         -ljumpC(ixCmin1:ixCmax1))/jumpC(ixCmin1:ixCmax1)
   elsewhere
      tmp(ixCmin1:ixCmax1)=adtdxC(ixCmin1:ixCmax1)
   end where
   !eq.3.45a with correction
   if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
      phiC(ixCmin1:ixCmax1)=-ljumpC(jxCmin1:jxCmax1)-ljumpC(ixCmin1:ixCmax1)&
         +jumpC(ixCmin1:ixCmax1)*abs(tmp(ixCmin1:ixCmax1))
   else
      where(abs(tmp(ixCmin1:ixCmax1))>=smallaC(ixCmin1:ixCmax1))
         phiC(ixCmin1:ixCmax1)=-ljumpC(jxCmin1:jxCmax1)-ljumpC&
            (ixCmin1:ixCmax1)+jumpC(ixCmin1:ixCmax1)*abs(tmp(ixCmin1:ixCmax1))
      elsewhere
         phiC(ixCmin1:ixCmax1)=-ljumpC(jxCmin1:jxCmax1)-ljumpC&
            (ixCmin1:ixCmax1)+jumpC(ixCmin1:ixCmax1)*(half*smallaC&
            (ixCmin1:ixCmax1)+half/smallaC(ixCmin1:ixCmax1)&
            *tmp(ixCmin1:ixCmax1)**2)
      endwhere
   endif
   !extra -(a*lambda)**2*delta
case default
   call mpistop("Error in TVDLimit: Unknown TVD type")
end select

end subroutine getphi
!=============================================================================
subroutine entropyfix(ixmin1,ixmax1,il,aL,aR,a,smalla)

! Apply entropyfix based on typeentropy(il),aL,aR, and a
! Calculate "smalla" (Harten,Powell) or modify "a" (ratio)

include 'amrvacdef.f'

integer, intent(in) :: ixmin1,ixmax1, il
double precision, dimension(ixGlo1:ixGhi1) :: aL, aR, a, smalla
!-----------------------------------------------------------------------------

select case(typeentropy(il))
case('harten')
   smalla(ixmin1:ixmax1)=max(zero,a(ixmin1:ixmax1)-aL(ixmin1:ixmax1),&
      aR(ixmin1:ixmax1)-a(ixmin1:ixmax1))
case('powell')
   smalla(ixmin1:ixmax1)=max(zero,two*(aR(ixmin1:ixmax1)-aL(ixmin1:ixmax1)))
!!case('ratio')
!!   where(aL(ix^S)<zero .and. aR(ix^S)>zero)&
!!      a(ix^S)=a(ix^S)-2*aR(ix^S)*aL(ix^S)/(aR(ix^S)-aL(ix^S))
case('yee')
   ! This has been done in geteigenjump already
case('nul')
   ! No entropyfix is applied
case default
   call mpistop("No such type of entropy fix")
end select

end subroutine entropyfix
!=============================================================================
! end module vactvd
!##############################################################################
