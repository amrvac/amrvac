module mod_srhd_roe
! module vacphys.srhdroe - subroutines for Roe-type Riemann solver for SRHD

  integer,parameter :: soundRW_ = 1,soundLW_=2,entropW_=3,shearW0_=3 ! waves
  integer,parameter :: nworkroe = 3

contains
  
subroutine average(wL,wR,x,ix^L,idim,wroe,workroe)

! Calculate the Roe average of w, assignment of variables:
! rho -> v0, m -> v, e -> v4
! check thesis Eulderink pg.18

  use mod_global_parameters

  integer:: ix^L,idim,idir
  double precision, dimension(ixG^T,nw):: wL,wR,wroe
  double precision, dimension(ixG^T,nworkroe):: workroe
  double precision, dimension(ixG^T,1:ndim):: x
  !---------------------------------------------------------------------------
  call average2(wL,wR,x,ix^L,idim,wroe,workroe(ixG^T,1),workroe(ixG^T,2))

end subroutine average
!=============================================================================
subroutine average2(wL,wR,x,ix^L,idim,wroe,tmp,tmp2)

! Calculate the Roe average of w, assignment of variables:
! rho -> v0, m -> v, e -> v4
! check thesis Eulderink pg.18

  use mod_global_parameters
  
  integer:: ixG^L,ix^L,idim,idir
  double precision, dimension(ixG^T,nw):: wL,wR,wroe
  double precision, dimension(ixG^T,1:ndim):: x
  double precision, dimension(ixG^T):: tmp,tmp2,lfL,lfR
  !---------------------------------------------------------------------------
  
  call getaux(.true.,wL,x,ixG^LL,ix^L,'average2_wL')
  call getaux(.true.,wR,x,ixG^LL,ix^L,'average2_wR')
  
  ! Calculate K_L
  tmp(ix^S) =sqrt((wL(ix^S,d_)/wL(ix^S,lfac_))+ &
       eqpar(gamma_)*wL(ix^S,p_)/(eqpar(gamma_)-1) )
  ! Calculate K_R, K=sqrt(rho*h)
  tmp2(ix^S) =sqrt((wR(ix^S,d_)/wR(ix^S,lfac_))+ &
       eqpar(gamma_)*wR(ix^S,p_)/(eqpar(gamma_)-1) )

  !!! Lorentz factor
  !!lfL(ix^S)=1/sqrt(1-(^C&wL(ix^S,v^C_)**2+))
  !!lfR(ix^S)=1/sqrt(1-(^C&wR(ix^S,v^C_)**2+))
  
  ! V^0, see thesis Eulderink
  wroe(ix^S,d_)=(wR(ix^S,lfac_)*tmp2(ix^S)+tmp(ix^S)*wL(ix^S,lfac_))/&
       (tmp2(ix^S)+tmp(ix^S))
  
  ! V^i, see thesis Eulderink, equivalent to Roe-average velocities
  do idir=1,ndir
     wroe(ix^S,s0_+idir)=( (wR(ix^S,s0_+idir)/(wR(ix^S,lfac_)*tmp2(ix^S))) &
                          +(wL(ix^S,s0_+idir)/(wL(ix^S,lfac_)*tmp(ix^S))) &
                         )/(tmp2(ix^S)+tmp(ix^S))
  end do
  
  ! V^4, see thesis Eulderink
  wroe(ix^S,tau_)=(wR(ix^S,p_)/tmp2(ix^S)+wL(ix^S,p_)/tmp(ix^S))/&
       (tmp2(ix^S)+tmp(ix^S))
  
end subroutine average2
!=============================================================================
subroutine geteigenjump(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,workroe)

! Calculate the il-th characteristic speed and the jump in the il-th 
! characteristic variable in the idim direction within ixL. 
! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

  use mod_global_parameters
  
  integer:: ix^L,il,idim
  double precision, dimension(ixG^T,nw):: wL,wR,wroe
  double precision, dimension(ixG^T)   :: smalla,a,jump
  double precision, dimension(ixG^T,1:ndim):: x
  double precision, dimension(ixG^T,nworkroe) :: workroe
  !---------------------------------------------------------------------------
  call geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump, &
       workroe(ixG^T,1),workroe(ixG^T,2),workroe(ixG^T,3))

end subroutine geteigenjump
!=============================================================================
subroutine geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump, &
     csound,del,dv)

! Calculate the il-th characteristic speed and the jump in the il-th 
! characteristic variable in the idim direction within ixL. 
! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

  use mod_global_parameters
  
  integer:: ix^L,il,idim,idir
  double precision, dimension(ixG^T,nw):: wL,wR,wroe
  double precision, dimension(ixG^T)   :: smalla,a,jump,tmp,tmp2
  double precision, dimension(ixG^T,1:ndim):: x
  double precision, dimension(ixG^T)   :: csound,del,dv,del0,cp,e,k,y2
  !!save dpperc2,dvperc
  !!common /roe/ csound
  !---------------------------------------------------------------------------

  if(il==1)then
     !Square of sound speed: s^2=0.5*gam*v4*(1+v0^2-v^2)-0.5(gam-1)(1-v0^2+v^2)
     csound(ix^S)=half*eqpar(gamma_)*wroe(ix^S,tau_)*(one+ &
     wroe(ix^S,d_)*wroe(ix^S,d_)-(^C&wroe(ix^S,s^C_)**2+))-half* &
     (eqpar(gamma_)-one)*(one-wroe(ix^S,d_)*wroe(ix^S,d_)+(^C&wroe(ix^S,s^C_)**2+))
     
     ! Make sure that csound**2 is positive
     !csound(ix^S)=max(eqpar(gamma_)*smalldouble/wroe(ix^S,d_),csound(ix^S))
     
     ! Calculate uR-uL
     del(ix^S)=wR(ix^S,d_)-wL(ix^S,d_)
     dv(ix^S)=wR(ix^S,s0_+idim)-wL(ix^S,s0_+idim)
     del0(ix^S)=wR(ix^S,tau_)-wL(ix^S,tau_)
     
     !Now get the correct sound speed
     csound(ix^S)=sqrt(csound(ix^S))
  endif

  !Some help variables
  cp(ix^S)=one+eqpar(gamma_)*wroe(ix^S,tau_)/(eqpar(gamma_)-one)
  e(ix^S)=wroe(ix^S,d_)*wroe(ix^S,d_)-wroe(ix^S,s0_+idim)*wroe(ix^S,s0_+idim)
  k(ix^S)=wroe(ix^S,d_)*(del0(ix^S)+del(ix^S))-wroe(ix^S,s0_+idim)*dv(ix^S)
  y2(ix^S)=(one-eqpar(gamma_)*wroe(ix^S,tau_))*e(ix^S)+csound(ix^S)*csound(ix^S)

  select case(il)
  case(soundRW_)
     !lambda+=lambda2=((1-g*v4)*v0*v1+s*y)/((1-g*v4)*v0*v0+s^2)
     a(ix^S)=((one-eqpar(gamma_)*wroe(ix^S,tau_))*wroe(ix^S,d_)*wroe(ix^S,s0_+idim)+&
          csound(ix^S)*sqrt(y2(ix^S)))/((one-eqpar(gamma_)*wroe(ix^S,tau_))*&
          wroe(ix^S,d_)*wroe(ix^S,d_)+csound(ix^S)*csound(ix^S))
     !alp2=(s^2*k-s*y*(v0*dv-v1*(del0+del)+(g-1)*e*(del+cp*(-v0*(del0+del)+v1*dv)))/(-2*e*s^2)
     jump(ix^S)=(csound(ix^S)*csound(ix^S)*k(ix^S)-csound(ix^S)*sqrt(y2(ix^S))*&
          (wroe(ix^S,d_)*dv(ix^S)-wroe(ix^S,s0_+idim)*(del0(ix^S)+&
          del(ix^S)))+(eqpar(gamma_)-one)*e(ix^S)*(del(ix^S)+&
          cp(ix^S)*(-wroe(ix^S,d_)*(del0(ix^S)+del(ix^S))+&
          (^C&wroe(ix^S,s^C_)*(wR(ix^S,s^C_)-wL(ix^S,s^C_))+)  )))/&
          (-2*e(ix^S)*csound(ix^S)*csound(ix^S))
  case(soundLW_)
     !lambda-=lambda1=((1-g*v4)*v0*v1-s*y)/((1-g*v4)*v0*v0+s^2)
     a(ix^S)=((one-eqpar(gamma_)*wroe(ix^S,tau_))*wroe(ix^S,d_)*wroe(ix^S,s0_+idim)-&
          csound(ix^S)*sqrt(y2(ix^S)))/((one-eqpar(gamma_)*wroe(ix^S,tau_))*&
          wroe(ix^S,d_)*wroe(ix^S,d_)+csound(ix^S)*csound(ix^S))
     !alp1=(s^2*k+s*y*(v0*dv-v1*(del0+del)+(g-1)*e*(del+cp*(-v0*(del0+del)+v1*dv)))/(-2*e*s^2)
     jump(ix^S)=(csound(ix^S)*csound(ix^S)*k(ix^S)+csound(ix^S)*sqrt(y2(ix^S))*&
          (wroe(ix^S,d_)*dv(ix^S)-wroe(ix^S,s0_+idim)*(del0(ix^S)+&
          del(ix^S)))+(eqpar(gamma_)-one)*e(ix^S)*(del(ix^S)+&
          cp(ix^S)*(-wroe(ix^S,d_)*(del0(ix^S)+del(ix^S))+&
          (^C&wroe(ix^S,s^C_)*(wR(ix^S,s^C_)-wL(ix^S,s^C_))+) )))/&
          (-two*e(ix^S)*csound(ix^S)*csound(ix^S))
  case(entropW_)
     !lambda0=lambda3=v1/v0
     a(ix^S)=wroe(ix^S,s0_+idim)/wroe(ix^S,d_)
     !alp3=(2*s^2*k+(g-1)*e*(del+cp*(-v0*(del0+del)+v1*dv)))/(e*s^2)
     jump(ix^S)=(two*csound(ix^S)*csound(ix^S)*k(ix^S)+(eqpar(gamma_)-one)*e(ix^S)*&
          (del(ix^S)+cp(ix^S)*(-wroe(ix^S,d_)*(del0(ix^S)+del(ix^S))+&
          (^C&wroe(ix^S,s^C_)*(wR(ix^S,s^C_)-wL(ix^S,s^C_))+) )))/&
          (e(ix^S)*csound(ix^S)*csound(ix^S))
  case default
     !Determine the direction of the shear wave
     idir=il-shearW0_; if(idir>=idim)idir=idir+1
     a(ix^S)=wroe(ix^S,s0_+idim)/wroe(ix^S,d_)
     !alp4_5=del2_3-kv2_3/e
     jump(ix^S)=wR(ix^S,s0_+idir)-wL(ix^S,s0_+idir)-k(ix^S)*wroe(ix^S,s0_+idir)/e(ix^S)
  end select

  ! Calculate "smalla" or modify "a" based on the "typeentropy" switch
  ! Put left and right eigenvalues, if needed, into tmp and tmp2
  ! OK, since subroutines getpthermal and entropyfix do not use tmp and tmp2

  select case(typeentropy(il))
  case('yee')
     ! Based on Yee JCP 68,151 eq 3.23
     smalla(ix^S)=entropycoef(il)
  case('harten','powell')
     call getaux(.true.,wL,x,ixG^LL,ix^L,'geteigenjump2_wL')
     call getaux(.true.,wR,x,ixG^LL,ix^L,'geteigenjump2_wR')
     ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
     select case(il)
     case(soundRW_)
        tmp(ix^S)=wL(ix^S,s0_+idim)/wL(ix^S,d_)&
             + sqrt(eqpar(gamma_)*wL(ix^S,p_)/wL(ix^S,d_))
        tmp2(ix^S)=wR(ix^S,s0_+idim)/wR(ix^S,d_)&
             + sqrt(eqpar(gamma_)*wR(ix^S,p_)/wR(ix^S,d_))
     case(soundLW_)
        tmp(ix^S)=wL(ix^S,s0_+idim)/wL(ix^S,d_)&
             - sqrt(eqpar(gamma_)*wL(ix^S,p_)/wL(ix^S,d_))
        tmp2(ix^S)=wR(ix^S,s0_+idim)/wR(ix^S,d_)&
             - sqrt(eqpar(gamma_)*wR(ix^S,p_)/wR(ix^S,d_))
     case default
        tmp(ix^S) =wL(ix^S,s0_+idim)/wL(ix^S,d_)
        tmp2(ix^S)=wR(ix^S,s0_+idim)/wR(ix^S,d_)
     end select
  end select
  
  call entropyfix(ix^L,il,tmp,tmp2,a,smalla)
  
end subroutine geteigenjump2
!=============================================================================
subroutine rtimes(q,wroe,ix^L,iw,il,idim,rq,workroe)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

  use mod_global_parameters
  
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

  use mod_global_parameters

  integer::          ix^L,iw,il,idim,idir
  double precision:: wroe(ixG^T,nw)
  double precision, dimension(ixG^T):: q,rq,csound,cm,cp,e,y
  logical:: shearwave
  !!common /roe/ csound
!-----------------------------------------------------------------------------

  shearwave=il>shearW0_
  if(shearwave)then
     ! Direction of shearwave increases with il plus idir==idim is jumped over
     idir=il-shearW0_; if(idir>=idim)idir=idir+1
  endif

  cm(ix^S)=one-eqpar(gamma_)*wroe(ix^S,tau_)/(eqpar(gamma_)-1)
  cp(ix^S)=one+eqpar(gamma_)*wroe(ix^S,tau_)/(eqpar(gamma_)-1)
  e(ix^S)=wroe(ix^S,d_)*wroe(ix^S,d_)-&
       wroe(ix^S,s0_+idim)*wroe(ix^S,s0_+idim)
  y(ix^S)=sqrt((one-eqpar(gamma_)*wroe(ix^S,tau_))*e(ix^S)+&
       csound(ix^S)*csound(ix^S))
  
  select case(iw)
  case(d_)
     select case(il)
     case(soundRW_)
        rq(ix^S)=q(ix^S)*cm(ix^S)
     case(soundLW_)
        rq(ix^S)=q(ix^S)*cm(ix^S)
     case(entropW_)
        rq(ix^S)=q(ix^S)*(cm(ix^S)+csound(ix^S)*csound(ix^S)&
             /(eqpar(gamma_)-one))
     case default
        rq(ix^S)=-q(ix^S)*wroe(ix^S,s0_+idir)*cp(ix^S)
     end select
  case(tau_)
     select case(il)
     case(soundRW_)
        rq(ix^S)=q(ix^S)*(wroe(ix^S,d_)+csound(ix^S)*&
             wroe(ix^S,s0_+idim)/y(ix^S)-cm(ix^S))
     case(soundLW_)
        rq(ix^S)=q(ix^S)*(wroe(ix^S,d_)-csound(ix^S)*&
             wroe(ix^S,s0_+idim)/y(ix^S)-cm(ix^S))
     case(entropW_)
        rq(ix^S)=q(ix^S)*(wroe(ix^S,d_)-cm(ix^S)-&
             csound(ix^S)*csound(ix^S)/(eqpar(gamma_)-1))
     case default
        rq(ix^S)=q(ix^S)*wroe(ix^S,s0_+idir)*cp(ix^S)
     end select
  case default
     if(iw==s0_+idim)then
        select case(il)
        case(soundRW_)
           rq(ix^S)=q(ix^S)*(wroe(ix^S,s0_+idim)+&
                csound(ix^S)*wroe(ix^S,d_)/y(ix^S))
        case(soundLW_)
           rq(ix^S)=q(ix^S)*(wroe(ix^S,s0_+idim)-&
                csound(ix^S)*wroe(ix^S,d_)/y(ix^S))
        case(entropW_)
           rq(ix^S)=q(ix^S)*wroe(ix^S,s0_+idim)
        case default
           rq(ix^S)=zero
        end select
     else
        if(shearwave)then
           if(iw==s0_+idir)then
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

end module mod_srhd_roe
