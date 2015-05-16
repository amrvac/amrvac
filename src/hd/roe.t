!##############################################################################
! module amrvacphys.hdroe - subroutines for Roe-type Riemann solver for HD
!=============================================================================
{#IFDEF GAMMA
subroutine average(wL,wR,x,ix^L,idim,wroe,workroe)

! Calculate the Roe average of w, assignment of variables:
! rho -> rho, m -> v, e -> h

include 'amrvacdef.f'

integer:: ix^L,idim,idir
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, dimension(ixG^T,ndim), intent(in) :: x
double precision, dimension(ixG^T,nworkroe):: workroe
!-----------------------------------------------------------------------------
call average2(wL,wR,x,ix^L,idim,wroe,workroe(ixG^T,1),workroe(ixG^T,2))

end subroutine average
!=============================================================================
subroutine average2(wL,wR,x,ix^L,idim,wroe,tmp,tmp2)

! Calculate the Roe average of w, assignment of variables:
! rho -> rho, m -> v, e -> h

include 'amrvacdef.f'

integer:: ix^L,idim,idir
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, dimension(ixG^T,ndim), intent(in) :: x
double precision, dimension(ixG^T):: tmp,tmp2
!-----------------------------------------------------------------------------
oktest=index(teststr,'average')>=1

tmp(ix^S) =sqrt(wL(ix^S,rho_))
tmp2(ix^S)=sqrt(wR(ix^S,rho_))
! The averaged density is sqrt(rhoL*rhoR)
wroe(ix^S,rho_)=tmp(ix^S)*tmp2(ix^S)

! Now the ratio sqrt(rhoL/rhoR) is put into tmp
tmp(ix^S)=tmp(ix^S)/tmp2(ix^S)
! Roe-average velocities
do idir=1,ndir
   wroe(ix^S,m0_+idir)=(wL(ix^S,m0_+idir)/wL(ix^S,rho_)*tmp(ix^S)+&
                        wR(ix^S,m0_+idir)/wR(ix^S,rho_))/(one+tmp(ix^S))
end do
! Calculate enthalpyL, then enthalpyR, then Roe-average. Use tmp2 for pressure.
call getpthermal(wL,x,ixG^LL,ix^L,tmp2)
if(oktest)write(*,*)'pL:',tmp2(ixtest^D)
wroe(ix^S,e_)=(tmp2(ix^S)+wL(ix^S,e_))/wL(ix^S,rho_)
if(oktest)write(*,*)'hL:',wroe(ixtest^D,e_)
call getpthermal(wR,x,ixG^LL,ix^L,tmp2)
if(oktest)write(*,*)'pR:',tmp2(ixtest^D)
tmp2(ix^S)=(tmp2(ix^S)+wR(ix^S,e_))/wR(ix^S,rho_)
if(oktest)write(*,*)'hR:',tmp2(ixtest^D)
wroe(ix^S,e_)=(wroe(ix^S,e_)*tmp(ix^S)+tmp2(ix^S))/(one+tmp(ix^S))
if(oktest)write(*,*)'weight,h :',tmp(ixtest^D),wroe(ixtest^D,e_)

if(oktest)write(*,*)'Roe average in direction',idim
if(oktest)write(*,*)'rho,u,h:',&
    wroe(ixtest^D,rho_),wroe(ixtest^D,m0_+idim),wroe(ixtest^D,e_)

end subroutine average2
!=============================================================================
subroutine geteigenjump(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,workroe)

! Calculate the il-th characteristic speed and the jump in the il-th 
! characteristic variable in the idim direction within ixL. 
! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

include 'amrvacdef.f'

integer:: ix^L,il,idim
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, dimension(ixG^T,ndim), intent(in) :: x
double precision, dimension(ixG^T)   :: smalla,a,jump
double precision, dimension(ixG^T,nworkroe) :: workroe
!-----------------------------------------------------------------------------
call geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump, &
   workroe(ixG^T,1),workroe(ixG^T,2),workroe(ixG^T,3))

end subroutine geteigenjump
!=============================================================================
subroutine geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,&
    csound,dpperc2,dvperc)

! Calculate the il-th characteristic speed and the jump in the il-th 
! characteristic variable in the idim direction within ixL. 
! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

include 'amrvacdef.f'

integer:: ix^L,il,idim,idir
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, dimension(ixG^T,ndim), intent(in) :: x
double precision, dimension(ixG^T)   :: smalla,a,jump,tmp,tmp2
double precision, dimension(ixG^T)   :: csound,dpperc2,dvperc
!!save dpperc2,dvperc
!!common /roe/ csound
!-----------------------------------------------------------------------------

oktest=index(teststr,'geteigenjump')>=1

if(il==1)then
   !First calculate the square of the sound speed: c**2=(gamma-1)*(h-0.5*v**2)
   csound(ix^S)=(eqpar(gamma_)-one)*(wroe(ix^S,e_)-half*(^C&wroe(ix^S,m^C_)**2+))
   if(oktest)write(*,*)'csound**2:',csound(ixtest^D)
   ! Make sure that csound**2 is positive
   csound(ix^S)=max(eqpar(gamma_)*smalldouble/wroe(ix^S,rho_),csound(ix^S))

   ! Calculate (pR-pL)/c**2
   ! To save memory we use tmp amnd tmp2 for pL and pR (getpthermal is OK)
   call getpthermal(wL,x,ixG^LL,ix^L,tmp)
   call getpthermal(wR,x,ixG^LL,ix^L,tmp2)
   dpperc2(ix^S)=(tmp2(ix^S)-tmp(ix^S))/csound(ix^S)

   !Now get the correct sound speed
   csound(ix^S)=sqrt(csound(ix^S))

   ! Calculate (vR_idim-vL_idim)/c
   dvperc(ix^S)=(wR(ix^S,m0_+idim)/wR(ix^S,rho_)-&
                 wL(ix^S,m0_+idim)/wL(ix^S,rho_))/csound(ix^S)

   if(oktest)write(*,*)'gamma,h,u,v',&
       eqpar(gamma_),wroe(ixtest^D,e_),^C&wroe(ixtest^D,m^C_)

   if(oktest)write(*,*)'csound :',csound(ixtest^D)
   if(oktest)write(*,*)'v_idim :',wroe(ixtest^D,m0_+idim)
   if(oktest)write(*,*)'dpperc2:',dpperc2(ixtest^D)
   if(oktest)write(*,*)'dvperc :',dvperc(ixtest^D)
endif

select case(il)
   case(soundRW_)
      a(ix^S)=wroe(ix^S,m0_+idim)+csound(ix^S)
      jump(ix^S)=half*(dpperc2(ix^S)+wroe(ix^S,rho_)*dvperc(ix^S))
   case(soundLW_)
      a(ix^S)=wroe(ix^S,m0_+idim)-csound(ix^S)
      jump(ix^S)=half*(dpperc2(ix^S)-wroe(ix^S,rho_)*dvperc(ix^S))
   case(entropW_)
      a(ix^S)=wroe(ix^S,m0_+idim)
      jump(ix^S)=-dpperc2(ix^S)+wR(ix^S,rho_)-wL(ix^S,rho_)
   case default
      !Determine the direction of the shear wave
      idir=il-shearW0_; if(idir>=idim)idir=idir+1
      a(ix^S)=wroe(ix^S,m0_+idim)
      jump(ix^S)=wroe(ix^S,rho_)*&
         (wR(ix^S,m0_+idir)/wR(ix^S,rho_)-wL(ix^S,m0_+idir)/wL(ix^S,rho_))
end select

! Calculate "smalla" or modify "a" based on the "typeentropy" switch
! Put left and right eigenvalues, if needed, into tmp and tmp2
! OK, since subroutines getpthermal and entropyfix do not use tmp and tmp2

select case(typeentropy(il))
case('yee')
   ! Based on Yee JCP 68,151 eq 3.23
   smalla(ix^S)=entropycoef(il)
case('harten','powell')
   ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
   select case(il)
   case(soundRW_)
      call getpthermal(wL,x,ixG^LL,ix^L,tmp)
      tmp(ix^S)=wL(ix^S,m0_+idim)/wL(ix^S,rho_)&
           + sqrt(eqpar(gamma_)*tmp(ix^S)/wL(ix^S,rho_))
      call getpthermal(wR,x,ixG^LL,ix^L,tmp2)
      tmp2(ix^S)=wR(ix^S,m0_+idim)/wR(ix^S,rho_)&
           + sqrt(eqpar(gamma_)*tmp2(ix^S)/wR(ix^S,rho_))
   case(soundLW_)
      call getpthermal(wL,x,ixG^LL,ix^L,tmp)
      tmp(ix^S)=wL(ix^S,m0_+idim)/wL(ix^S,rho_)&
           - sqrt(eqpar(gamma_)*tmp(ix^S)/wL(ix^S,rho_))
      call getpthermal(wR,x,ixG^LL,ix^L,tmp2)
      tmp2(ix^S)=wR(ix^S,m0_+idim)/wR(ix^S,rho_)&
           - sqrt(eqpar(gamma_)*tmp2(ix^S)/wR(ix^S,rho_))
   case default
      tmp(ix^S) =wL(ix^S,m0_+idim)/wL(ix^S,rho_)
      tmp2(ix^S)=wR(ix^S,m0_+idim)/wR(ix^S,rho_)
   end select
end select

call entropyfix(ix^L,il,tmp,tmp2,a,smalla)

end subroutine geteigenjump2
!=============================================================================
subroutine rtimes(q,wroe,ix^L,iw,il,idim,rq,workroe)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

include 'amrvacdef.f'

integer::          ix^L,iw,il,idim
double precision:: wroe(ixG^T,nw)
double precision, dimension(ixG^T):: q,rq
double precision, dimension(ixG^T,nworkroe):: workroe
!-----------------------------------------------------------------------------
call rtimes2(q,wroe,ix^L,iw,il,idim,rq,workroe(ixG^T,1))

end subroutine rtimes
!=============================================================================
subroutine rtimes2(q,wroe,ix^L,iw,il,idim,rq,csound)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

include 'amrvacdef.f'

integer::          ix^L,iw,il,idim,idir
double precision:: wroe(ixG^T,nw)
double precision, dimension(ixG^T):: q,rq,csound
logical:: shearwave
!!common /roe/ csound
!-----------------------------------------------------------------------------

shearwave=il>shearW0_
if(shearwave)then
   ! Direction of shearwave increases with il plus idir==idim is jumped over
   idir=il-shearW0_; if(idir>=idim)idir=idir+1
endif

select case(iw)
case(rho_)
   if(shearwave)then 
      rq(ix^S)=zero
   else
      rq(ix^S)=q(ix^S)
   endif
case(e_)
   select case(il)
   case(soundRW_)
      rq(ix^S)=q(ix^S)*(wroe(ix^S,e_)+wroe(ix^S,m0_+idim)*csound(ix^S))
   case(soundLW_)
      rq(ix^S)=q(ix^S)*(wroe(ix^S,e_)-wroe(ix^S,m0_+idim)*csound(ix^S))
   case(entropW_)
      rq(ix^S)=q(ix^S)*half*( ^C&wroe(ix^S,m^C_)**2+ )
   case default
      rq(ix^S)=q(ix^S)*wroe(ix^S,m0_+idir)
   end select
case default
   if(iw==m0_+idim)then
      select case(il)
      case(soundRW_)
         rq(ix^S)=q(ix^S)*(wroe(ix^S,m0_+idim)+csound(ix^S))
      case(soundLW_)
         rq(ix^S)=q(ix^S)*(wroe(ix^S,m0_+idim)-csound(ix^S))
      case(entropW_)
         rq(ix^S)=q(ix^S)*wroe(ix^S,m0_+idim)
      case default
         rq(ix^S)=zero
      end select
   else
      if(shearwave)then
         if(iw==m0_+idir)then
            rq(ix^S)=q(ix^S)
         else
            rq(ix^S)=zero
         endif
      else
         rq(ix^S)=q(ix^S)*wroe(ix^S,iw)
      endif
   endif
end select

end subroutine rtimes2
!=============================================================================
! end module amrvacphys.hdroe
!##############################################################################
}{#IFDEF ISO
!##############################################################################
! module amrvacphys.hdadiabroe - subroutines for Roe-type Riemann solver for HDadiab
!=============================================================================
subroutine average(wL,wR,x,ix^L,idim,wroe,workroe)

! Calculate the Roe average of w, assignment of variables:
! rho -> rho, m -> v

include 'amrvacdef.f'

integer:: ix^L,idim
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, intent(in)    :: x(ixG^T,1:ndim)
double precision, dimension(ixG^T,nworkroe):: workroe
!-----------------------------------------------------------------------------
call average2(wL,wR,x,ix^L,idim,wroe,workroe(ixG^T,1))

end subroutine average
!=============================================================================
subroutine average2(wL,wR,x,ix^L,idim,wroe,tmp)

! Calculate the Roe average of w, assignment of variables:
! rho -> rho, m -> v

include 'amrvacdef.f'

integer:: ix^L,idim,idir
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, intent(in)    :: x(ixG^T,1:ndim)
double precision, dimension(ixG^T):: tmp
!-----------------------------------------------------------------------------

select case (typeaverage)
case ('arithmetic')
  ! This is the simple arithmetic average
   wroe(ix^S,rho_)=half*(wL(ix^S,rho_)+wR(ix^S,rho_))
   {^C& wroe(ix^S,m^C_)=&
      half*(wL(ix^S,m^C_)/wL(ix^S,rho_)+wR(ix^S,m^C_)/wR(ix^S,rho_)) \}
case ('roe','default')
  ! Calculate the Roe-average
  wroe(ix^S,rho_)=sqrt(wL(ix^S,rho_)*wR(ix^S,rho_))
  ! Roe-average velocities
  tmp(ix^S)=sqrt(wL(ix^S,rho_)/wR(ix^S,rho_))
  do idir=1,ndir
    wroe(ix^S,m0_+idir)=(wL(ix^S,m0_+idir)/wL(ix^S,rho_)*tmp(ix^S)+&
                        wR(ix^S,m0_+idir)/wR(ix^S,rho_))/(one+tmp(ix^S))
  end do
end select

end subroutine average2
!=============================================================================
subroutine geteigenjump(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,workroe)

! Calculate the il-th characteristic speed and the jump in the il-th 
! characteristic variable in the idim direction within ixL. 
! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

include 'amrvacdef.f'

integer:: ix^L,il,idim
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, intent(in)    :: x(ixG^T,1:ndim)
double precision, dimension(ixG^T)   :: smalla,a,jump
double precision, dimension(ixG^T,nworkroe) :: workroe
!-----------------------------------------------------------------------------
call geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,workroe(ixG^T,1))

end subroutine geteigenjump
!=============================================================================
subroutine geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,csound)

! Calculate the il-th characteristic speed and the jump in the il-th 
! characteristic variable in the idim direction within ixL. 
! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

include 'amrvacdef.f'

integer:: ix^L,il,idim,idir
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, intent(in)    :: x(ixG^T,1:ndim)
double precision, dimension(ixG^T)   :: smalla,a,jump,tmp,tmp2
double precision, dimension(ixG^T)   :: csound
DOUBLE PRECISION,PARAMETER:: qsmall=1.D-6
!-----------------------------------------------------------------------------

select case (typeaverage)
case ('arithmetic')
   csound(ix^S)=sqrt(eqpar(adiab_)*eqpar(gamma_)*&
      wroe(ix^S,rho_)**(eqpar(gamma_)-one))
   ! This is the original simple Roe-solver
   select case(il)
   case(soundRW_)
      a(ix^S)=wroe(ix^S,m0_+idim)+csound(ix^S)
      jump(ix^S)=half*((one-wroe(ix^S,m0_+idim)/csound(ix^S))*&
               (wR(ix^S,rho_)-wL(ix^S,rho_))&
              +(wR(ix^S,m0_+idim)-wL(ix^S,m0_+idim))/csound(ix^S))
   case(soundLW_)
      a(ix^S)=wroe(ix^S,m0_+idim)-csound(ix^S)
      jump(ix^S)=half*((one+wroe(ix^S,m0_+idim)/csound(ix^S))*&
               (wR(ix^S,rho_)-wL(ix^S,rho_))&
               -(wR(ix^S,m0_+idim)-wL(ix^S,m0_+idim))/csound(ix^S))
   case default
      ! Determine direction of shear wave
      idir=il-shearW0_; if(idir>=idim)idir=idir+1
      a(ix^S)=wroe(ix^S,m0_+idim)
      jump(ix^S)=-wroe(ix^S,m0_+idir)*(wR(ix^S,rho_)-wL(ix^S,rho_))&
               +(wR(ix^S,m0_+idir)-wL(ix^S,m0_+idir))
   end select
case ('roe','default')
   where(abs(wL(ix^S,rho_)-wR(ix^S,rho_))<=qsmall*(wL(ix^S,rho_)+wR(ix^S,rho_)))
     csound(ix^S)=sqrt(eqpar(adiab_)*eqpar(gamma_)*&
         wroe(ix^S,rho_)**(eqpar(gamma_)-one))
   elsewhere
     csound(ix^S)=sqrt(eqpar(adiab_)*(wR(ix^S,rho_)**eqpar(gamma_)-&
              wL(ix^S,rho_)**eqpar(gamma_))/(wR(ix^S,rho_)-wL(ix^S,rho_)))
   end where
   ! This is the Roe solver by Glaister
   ! based on P. Glaister JCP 93, 477-480 (1991)
   select case(il)
   case(soundRW_)
      a(ix^S)=wroe(ix^S,m0_+idim)+csound(ix^S)
      jump(ix^S)=half*((wR(ix^S,rho_)-wL(ix^S,rho_))+&
              wroe(ix^S,rho_)/csound(ix^S)*(wR(ix^S,m0_+idim)/wR(ix^S,rho_)-&
              wL(ix^S,m0_+idim)/wL(ix^S,rho_)))
   case(soundLW_)
      a(ix^S)=wroe(ix^S,m0_+idim)-csound(ix^S)
      jump(ix^S)=half*((wR(ix^S,rho_)-wL(ix^S,rho_))-&
              wroe(ix^S,rho_)/csound(ix^S)*(wR(ix^S,m0_+idim)/wR(ix^S,rho_)-&
              wL(ix^S,m0_+idim)/wL(ix^S,rho_)))
   case default
      ! Determine direction of shear wave
      idir=il-shearW0_; if(idir>=idim)idir=idir+1
      a(ix^S)=wroe(ix^S,m0_+idim)
      jump(ix^S)=wroe(ix^S,rho_)*(wR(ix^S,m0_+idir)/wR(ix^S,rho_)-&
               wL(ix^S,m0_+idir)/wL(ix^S,rho_))
   end select
end select

! Calculate "smalla" or modify "a" based on the "typeentropy" switch
! Use tmp and tmp2 for the left and right eigenvalues if needed
select case(typeentropy(il))
case('yee')
   ! Based on Yee JCP 68,151 eq 3.23
   smalla(ix^S)=entropycoef(il)
case('harten','powell')
   ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
   select case(il)
   case(soundRW_)
      tmp(ix^S) =wL(ix^S,m0_+idim)/wL(ix^S,rho_)&
         + sqrt(eqpar(adiab_)*eqpar(gamma_)*wL(ix^S,rho_)**(eqpar(gamma_)-one))
      tmp2(ix^S)=wR(ix^S,m0_+idim)/wR(ix^S,rho_)&
         + sqrt(eqpar(adiab_)*eqpar(gamma_)*wR(ix^S,rho_)**(eqpar(gamma_)-one))
   case(soundLW_)
      tmp(ix^S) =wL(ix^S,m0_+idim)/wL(ix^S,rho_)&
         - sqrt(eqpar(adiab_)*eqpar(gamma_)*wL(ix^S,rho_)**(eqpar(gamma_)-one))
      tmp2(ix^S)=wR(ix^S,m0_+idim)/wR(ix^S,rho_)&
         - sqrt(eqpar(adiab_)*eqpar(gamma_)*wR(ix^S,rho_)**(eqpar(gamma_)-one))
   case default
      tmp(ix^S) =wL(ix^S,m0_+idim)/wL(ix^S,rho_)
      tmp2(ix^S)=wR(ix^S,m0_+idim)/wR(ix^S,rho_)
   end select
end select

call entropyfix(ix^L,il,tmp,tmp2,a,smalla)

end subroutine geteigenjump2
!=============================================================================
subroutine rtimes(q,wroe,ix^L,iw,il,idim,rq,workroe)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

include 'amrvacdef.f'

integer::          ix^L,iw,il,idim
double precision:: wroe(ixG^T,nw)
double precision, dimension(ixG^T):: q,rq
double precision, dimension(ixG^T,nworkroe):: workroe
!-----------------------------------------------------------------------------

call rtimes2(q,wroe,ix^L,iw,il,idim,rq,workroe(ixG^T,1))

end subroutine rtimes
!=============================================================================
subroutine rtimes2(q,wroe,ix^L,iw,il,idim,rq,csound)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

include 'amrvacdef.f'

integer::          ix^L,iw,il,idim,idir
double precision:: wroe(ixG^T,nw)
double precision, dimension(ixG^T):: q,rq,csound
!-----------------------------------------------------------------------------

if(iw==rho_)then
   select case(il)
   case(soundRW_,soundLW_)
      rq(ix^S)=q(ix^S)
   case default
      rq(ix^S)=zero
   end select
else if(iw==m0_+idim)then
   select case(il)
   case(soundRW_)
      rq(ix^S)=q(ix^S)*(wroe(ix^S,m0_+idim)+csound(ix^S))
   case(soundLW_)
      rq(ix^S)=q(ix^S)*(wroe(ix^S,m0_+idim)-csound(ix^S))
   case default
      rq(ix^S)=zero
   end select
else
   select case(il)
   case(soundRW_,soundLW_)
      rq(ix^S)=q(ix^S)*wroe(ix^S,iw)
   case default
      !Determine direction of shear wave
      idir=il-shearW0_; if(idir>=idim)idir=idir+1
      if(iw==m0_+idir) then
         rq(ix^S)=q(ix^S)
      else
         rq(ix^S)=zero
      endif
   end select
endif

end subroutine rtimes2
!=============================================================================
! end module amrvacphys.hdadiabroe
!##############################################################################
}
