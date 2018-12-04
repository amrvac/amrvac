!HLLC for SRHD  (ZM)
module mod_srhd_hllc

contains

subroutine diffuse_hllcd(ixI^L,ixO^L,idims,wLC,wRC,fLC,fRC,patchf)

! when method is hllcd or hllcd1 then:

! this subroutine is to impose enforce regions where we AVOID HLLC 
! and use TVDLF instead: this is achieved by setting patchf to 4 in
! certain regions. An additional input parameter is nxdiffusehllc
! which sets the size of the fallback region.

use mod_global_parameters

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixI^S,1:nw), intent(in)      :: wRC,wLC
double precision, dimension(ixI^S,1:nwflux),intent(in) :: fLC, fRC

integer         , dimension(ixG^T), intent(inout)        :: patchf

integer                                           :: ixOO^D,TxOO^L
integer                                           :: iw
!-----------------------------------

if(all(abs(patchf(ixO^S))==1) &
  .and.all(wLC(ixO^S,1:nwflux)-wRC(ixO^S,1:nwflux)==zero) &
  .and.typeaxial=='slab')then
  ! determine parts of domain where solution is constant spatially
  ! and then enforce the switch to TVDLF
  patchf(ixO^S) = 4
else
  ! In a user-controlled region around any point with flux sign change between
  ! left and right, ensure fallback to TVDLF
  {do ixOO^D= ixO^LIM^D\}
    {
    TxOOmin^D= max(ixOO^D - nxdiffusehllc*kr(idims,^D), ixOmin^D);
    TxOOmax^D= min(ixOO^D + nxdiffusehllc*kr(idims,^D), ixOmax^D);
    \}
    if(abs(patchf(ixOO^D))==1 .or. abs(patchf(ixOO^D))==4)then
       if(any(fRC(ixOO^D,1:nwflux)*fLC(ixOO^D,1:nwflux)<-smalldouble))then
          where(abs(patchf(TxOO^S))==1)
            patchf(TxOO^S) = 4
          endwhere
       endif
    endif
  {enddo^D&\}
endif 

end subroutine diffuse_hllcd
!=============================================================================
subroutine getlCD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixI^L,ixO^L,&
                  whll,Fhll,lambdaCD,patchf)
  
! Calculate lambda at CD and set the patchf to know the orientation
! of the riemann fan and decide on the flux choice
! We also compute here the HLL flux and w value, for fallback strategy
  
use mod_global_parameters
  
integer, intent(in)                                        :: ixI^L,ixO^L,idims
double precision, dimension(ixI^S,1:nw), intent(in)        :: wLC,wRC
double precision, dimension(ixG^T,1:nwflux), intent(in)  :: fLC,fRC
double precision, dimension(ixG^T), intent(in)             :: cmax,cmin

integer         , dimension(ixG^T), intent(inout)          :: patchf

double precision, dimension(ixG^T,1:nwflux), intent(out) :: Fhll,whll
double precision, dimension(ixG^T), intent(out)            :: lambdaCD

double precision, dimension(ixG^T)      :: Aco,Bco,Cco,Delta
logical         , dimension(ixG^T)      :: Cond_patchf
integer                                 :: iw
double precision                        :: Epsilon
!--------------------------------------------
! on entry, patch is preset to contain values from -2,1,2,4
!      -2: take left flux, no computation here
!      +2: take right flux, no computation here
!      +4: take TVDLF flux, no computation here
!       1: compute the characteristic speed for the CD
Cond_patchf(ixO^S)=(abs(patchf(ixO^S))==1)

do iw=1,nwflux
  where(Cond_patchf(ixO^S))
  !============= compute HLL flux ==============!
  Fhll(ixO^S,iw)= (cmax(ixO^S)*fLC(ixO^S,iw)-cmin(ixO^S)*fRC(ixO^S,iw) &
                   + cmin(ixO^S)*cmax(ixO^S)*(wRC(ixO^S,iw)-wLC(ixO^S,iw)))&
                  /(cmax(ixO^S)-cmin(ixO^S))
  !======== compute intermediate HLL state =======!
  whll(ixO^S,iw) = (cmax(ixO^S)*wRC(ixO^S,iw)-cmin(ixO^S)*wLC(ixO^S,iw)&
                    +fLC(ixO^S,iw)-fRC(ixO^S,iw))&
                   /(cmax(ixO^S)-cmin(ixO^S))
  endwhere
enddo


!Calculate the Characteristic speed at the contact
! part specific to SRHD and HLLC
! Eq. 18, page 4, Mignone & Bodo, MNRAS 2005
!  write this eq as:  Aco lambda^2 - Bco lambda + Cco = 0
where(Cond_patchf(ixO^S))
  Aco(ixO^S)  =(Fhll(ixO^S,tau_)+Fhll(ixO^S,d_))
  Bco(ixO^S)  =whll(ixO^S,tau_)+whll(ixO^S,d_)+Fhll(ixO^S, s0_+idims)
  Cco(ixO^S)  =whll(ixO^S,s0_+idims)
  Delta(ixO^S)=Bco(ixO^S)**2.0d0-4.0d0*Aco(ixO^S)*Cco(ixO^S)
  where(Aco(ixO^S)/=zero.and.Delta(ixO^S)>=zero)
    !Calculate the Characteristic speed at the contact
    ! only the minus sign is between [-1,1]
    lambdaCD(ixO^S) =(Bco(ixO^S)-dsqrt(Delta(ixO^S)))/(2.0d0*Aco(ixO^S))
  elsewhere(Aco(ixO^S)==zero.and.Bco(ixO^S)/=zero)
     lambdaCD(ixO^S) = Cco(ixO^S)/Bco(ixO^S) 
  elsewhere(Delta(ixO^S)<zero)
     lambdaCD(ixO^S)= zero
     ! we will fall back to HLL flux case in this degeneracy
     patchf(ixO^S)  = 3
  endwhere
endwhere

where(patchf(ixO^S) ==  3)
  Cond_patchf(ixO^S)=.false.
end where

where(Cond_patchf(ixO^S))
  ! double check whether obtained speed is in between min and max speeds given
  ! and identify in which part of the Riemann fan the time-axis is
  where(cmin(ixO^S)<zero.and.lambdaCD(ixO^S)>zero.and.lambdaCD(ixO^S)<cmax(ixO^S))
    patchf(ixO^S) = -1
  elsewhere(cmax(ixO^S)>zero.and.lambdaCD(ixO^S)<zero.and.lambdaCD(ixO^S)>cmin(ixO^S))
    patchf(ixO^S) =  1
  elsewhere(lambdaCD(ixO^S)>=cmax(ixO^S).or.lambdaCD(ixO^S) <= cmin(ixO^S))
    lambdaCD(ixO^S) = zero
    ! we will fall back to HLL flux case in this degeneracy
    patchf(ixO^S) =  3
  endwhere
endwhere

where(patchf(ixO^S) ==  3)
  Cond_patchf(ixO^S)=.false.
end where

! handle the specific case where the time axis is exactly on the CD 
if(any(lambdaCD(ixO^S)==zero.and.Cond_patchf(ixO^S)))then
  ! determine which sector (forward or backward) of the Riemann fan is smallest
  ! and select left or right flux accordingly
  where(lambdaCD(ixO^S)==zero.and.Cond_patchf(ixO^S))
    where(-cmin(ixO^S)>=cmax(ixO^S))
      patchf(ixO^S) =  1
    elsewhere
      patchf(ixO^S) = -1
    endwhere
  endwhere
endif

!-----------------------------------------------------------------------!
! eigenvalue lambda for contact is near zero: decrease noise by this trick
if(flathllc)then
  Epsilon=1.0d-6
  where(Cond_patchf(ixO^S).and. &
    dabs(lambdaCD(ixO^S))/max(cmax(ixO^S),Epsilon)< Epsilon  .and. &
    dabs(lambdaCD(ixO^S))/max(dabs(cmin(ixO^S)),Epsilon)< Epsilon)
    lambdaCD(ixO^S) =  zero
  end where
end if 
!-----------------------------------------------------------------------!

! done with computations, remaining just for emergencies...
if(any(dabs(lambdaCD(ixO^S))>one .and. Cond_patchf(ixO^S)))then
 call mpistop("problems with lambdaCD>1")
endif

! next should never happen
if(any(patchf(ixO^S)==0))then
 call mpistop("patchf=0")
endif

return 
end subroutine getlCD
!=============================================================================
subroutine getwCD(wLC,wRC,whll,vLC,vRC,fRC,fLC,Fhll,patchf,&
                  lambdaCD,cmin,cmax,ixI^L,ixO^L,idims,f)

! compute the intermediate state U*
! only needed where patchf=-1/1

! For SRHD: compute D*, S*, tau*, p* and v* 
  
use mod_global_parameters

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixI^S,1:nw), intent(in)      :: wRC,wLC
double precision, dimension(ixG^T,1:nwflux), intent(in)  :: whll, Fhll
double precision, dimension(ixG^T), intent(in)           :: vRC,vLC,lambdaCD
double precision, dimension(ixG^T), intent(in)           :: cmax,cmin
double precision, dimension(ixG^T,1:nwflux), intent(in)  :: fRC,fLC
double precision, dimension(ixG^T,1:nwflux), intent(out) :: f

integer         , dimension(ixG^T), intent(in)           :: patchf

double precision, dimension(ixG^T,1:nw)      :: wCD,wSub
double precision, dimension(ixG^T)           :: cspeed,vsub,Ratio_CD
integer                                      :: iw
!--------------------------------------------

!-------------- auxiliary Speed and array-------------!
where(patchf(ixO^S)== 1)
  cspeed(ixO^S) = cmax(ixO^S)
  vSub(ixO^S)   = vRC(ixO^S)
elsewhere(patchf(ixO^S)== -1)
  cspeed(ixO^S) = cmin(ixO^S)
  vSub(ixO^S)   = vLC(ixO^S)
endwhere

do iw=1,nw
 where(patchf(ixO^S) == 1)
   wSub(ixO^S,iw)=wRC(ixO^S,iw)
 elsewhere(patchf(ixO^S) == -1)
   wSub(ixO^S,iw)=wLC(ixO^S,iw)
 endwhere
enddo

! compute intermediate state U* from eq. 16-17 in Mignone-Bodo, MNRAS 2008
where(abs(patchf(ixO^S))==1)
   Ratio_CD(ixO^S) = (cspeed(ixO^S)-vSub(ixO^S))/(cspeed(ixO^S)-lambdaCD(ixO^S))
   !--- Density ---!
   wCD(ixO^S,d_) = wSub(ixO^S,d_)*Ratio_CD(ixO^S)
{#IFDEF TRACER
   !--- Tracer ---!
{^FL& wCD(ixO^S,Dtr^FL_)   = wSub(ixO^S,Dtr^FL_)*Ratio_CD(ixO^S) \}
}
{#IFDEF EPSINF
  !--- Particle evolution ---!
  wCD(ixO^S,epsinf_)   = wSub(ixO^S,epsinf_)*Ratio_CD(ixO^S)
  wCD(ixO^S,ne_)        = wSub(ixO^S,ne_)*Ratio_CD(ixO^S)
  wCD(ixO^S,ne0_)       = wSub(ixO^S,ne0_)*Ratio_CD(ixO^S)
}
   !--- Pressure ---!
   wCD(ixO^S,p_) = (lambdaCD(ixO^S)*cspeed(ixO^S)&
                * (wSub(ixO^S,tau_)+wSub(ixO^S,d_))+wSub(ixO^S, p_)&
                - wSub(ixO^S,s0_+idims) &
                * (cspeed(ixO^S)+lambdaCD(ixO^S)-vSub(ixO^S)))&
                  /(one-cspeed(ixO^S)*lambdaCD(ixO^S))
   wCD(ixO^S,tau_) = (wSub(ixO^S,tau_)*(cspeed(ixO^S)-vSub(ixO^S))&
                 +wCD(ixO^S,p_)*lambdaCD(ixO^S)-wSub(ixO^S,p_)*vSub(ixO^S))&
                 /(cspeed(ixO^S)-lambdaCD(ixO^S))
endwhere

do iw =s1_,s^NC_
   where(abs(patchf(ixO^S))==1)
     wCD(ixO^S,iw) = wSub(ixO^S,iw)*Ratio_CD(ixO^S)
   endwhere
   if(iw==s0_+idims)then
      where(abs(patchf(ixO^S))==1)
        wCD(ixO^S,iw)=wCD(ixO^S,iw) &
          +(wCD(ixO^S,p_)-wSub(ixO^S,p_))/(cspeed(ixO^S)-lambdaCD(ixO^S))
      endwhere
   endif
enddo

do iw=1,nwflux
   where(patchf(ixO^S) == 1)
      f(ixO^S,iw)=fRC(ixO^S,iw)+cspeed(ixO^S)*(wCD(ixO^S,iw)-wsub(ixO^S,iw))
   elsewhere(patchf(ixO^S) == -1)
      f(ixO^S,iw)=fLC(ixO^S,iw)+cspeed(ixO^S)*(wCD(ixO^S,iw)-wsub(ixO^S,iw))
   endwhere
end do

return 
end subroutine getwCD

end module mod_srhd_hllc
