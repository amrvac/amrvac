!##############################################################################
! module amrvacphys.mhdhllc.t 
!=============================================================================
subroutine diffuse_hllcd(ixI^L,ixO^L,idims,wLC,wRC,fLC,fRC,patchf)

! when method is hllcd or hllcd1 then: 

! this subroutine is to enforce regions where we AVOID HLLC
! and use TVDLF instead: this is achieved by setting patchf to 4 in
! certain regions. An additional input parameter is nxdiffusehllc
! which sets the size of the fallback region.

include 'amrvacdef.f'

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixG^T,1:nw), intent(in)      :: wRC,wLC
double precision, dimension(ixG^T,1:nwflux),intent(in) :: fLC, fRC

integer         , dimension(ixG^T), intent(inout)        :: patchf

integer                                           :: ixOO^D,TxOO^L
!-----------------------------------

! In a user-controlled region around any point with flux sign change between
! left and right, ensure fallback to TVDLF
{do ixOO^D= ixO^LIM^D\}
  {
  TxOOmin^D= max(ixOO^D - nxdiffusehllc*kr(idims,^D), ixOmin^D);
  TxOOmax^D= min(ixOO^D + nxdiffusehllc*kr(idims,^D), ixOmax^D);
  \}
  if(abs(patchf(ixOO^D)) == 1 .or. abs(patchf(ixOO^D)) == 4)Then
     if(any(fRC(ixOO^D,1:nwflux)*fLC(ixOO^D,1:nwflux)<-smalldouble))Then
       where(Abs(patchf(TxOO^S))==1)
         patchf(TxOO^S) = 4
       endwhere
     endif
  endif
{enddo^D&\}

end subroutine diffuse_hllcd
!=============================================================================
subroutine getlCD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixI^L,ixO^L, &
                  whll,Fhll,lambdaCD,patchf)

! Calculate lambda at CD and set the patchf to know the orientation
! of the riemann fan and decide on the flux choice
! We also compute here the HLL flux and w value, for fallback strategy

include 'amrvacdef.f'

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixG^T,1:nw), intent(in)      :: wLC,wRC
double precision, dimension(ixG^T,1:nwflux), intent(in):: fLC,fRC
double precision, dimension(ixG^T), intent(in)           :: cmax,cmin

integer         , dimension(ixG^T), intent(inout)        :: patchf

double precision, dimension(ixG^T,1:nwflux), intent(out) :: Fhll,whll
double precision, dimension(ixG^T), intent(out)            :: lambdaCD


logical         , dimension(ixG^T)     :: Cond_patchf
double precision                       :: Epsilon
integer                                :: iw
!--------------------------------------------
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
                    +fLC(ixO^S,iw)-fRC(ixO^S,iw))/(cmax(ixO^S)-cmin(ixO^S))
  endwhere
enddo

! deduce the characteristic speed at the CD
where(Cond_patchf(ixO^S))
  lambdaCD(ixO^S)=whll(ixO^S,m0_+idims)/whll(ixO^S,rho_)
end where


where(Cond_patchf(ixO^S))
  ! double check whether obtained speed is in between min and max speeds given
  ! and identify in which part of the Riemann fan the time-axis is
  where(cmin(ixO^S)<zero.and.lambdaCD(ixO^S)>zero&
        .and.lambdaCD(ixO^S)<cmax(ixO^S))
    patchf(ixO^S) = -1
  elsewhere(cmax(ixO^S)>zero.and.lambdaCD(ixO^S)<zero&
            .and.lambdaCD(ixO^S)>cmin(ixO^S))
    patchf(ixO^S) =  1
  elsewhere(lambdaCD(ixO^S)>=cmax(ixO^S).or.lambdaCD(ixO^S) <= cmin(ixO^S))
    lambdaCD(ixO^S) = zero
    ! we will fall back to HLL flux case in this degeneracy
    patchf(ixO^S) =  3
  endwhere
endwhere


where(patchf(ixO^S)== 3)
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

return
end subroutine getlCD
!=============================================================================
subroutine getwCD(wLC,wRC,whll,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaCD,cmin,cmax,&
                  ixI^L,ixO^L,idims,f)

! compute the intermediate state U*
! only needed where patchf=-1/1

! reference Li S., JCP, 203, 2005, 344-357
! reference T. Miyoski, Kusano JCP, 2008, 2005.

include 'amrvacdef.f'

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixG^T,1:nw), intent(in)      :: wRC,wLC
double precision, dimension(ixG^T,1:nwflux), intent(in):: whll, Fhll
double precision, dimension(ixG^T), intent(in)           :: vRC, vLC,lambdaCD
double precision, dimension(ixG^T), intent(in)           :: cmax,cmin
double precision, dimension(ixG^T,1:nwflux), intent(in):: fRC,fLC
double precision, dimension(ixG^T,1:nwflux),intent(out):: f

integer         , dimension(ixG^T), intent(in)           :: patchf

double precision, dimension(ixG^T,1:nw)        :: wCD,wSub
double precision, dimension(ixG^T,1:nwflux)    :: fSub
double precision, dimension(ixG^T)             :: vSub,cspeed,pCD,VdotBCD
integer                                        :: iw
!--------------------------------------------

!-------------- auxiliary Speed and array-------------!
where(patchf(ixO^S) == 1)
  cspeed(ixO^S) = cmax(ixO^S)
  vSub(ixO^S)   = vRC(ixO^S)
elsewhere(patchf(ixO^S) == -1)
  cspeed(ixO^S) = cmin(ixO^S)
  vSub(ixO^S)   = vLC(ixO^S)
endwhere

do iw=1,nw
  if(iw /= b0_+idims)then
    where(patchf(ixO^S) == 1)
      wSub(ixO^S,iw) =  wRC(ixO^S,iw)
    elsewhere(patchf(ixO^S) == -1)
      wSub(ixO^S,iw) =  wLC(ixO^S,iw)
    endwhere
  endif
enddo

do iw=1,nwflux
  if(iw /= b0_+idims)then
    where(patchf(ixO^S) == 1)
      fSub(ixO^S,iw) =  fRC(ixO^S,iw)
    elsewhere(patchf(ixO^S) == -1)
      fSub(ixO^S,iw) =  fLC(ixO^S,iw)
    endwhere
  endif
enddo

where(abs(patchf(ixO^S))==1)
  wCD(ixO^S,rho_) = wSub(ixO^S,rho_)&
                   *(cspeed(ixO^S)-vSub(ixO^S))/(cspeed(ixO^S)-lambdaCD(ixO^S))
{#IFDEF TRACER
{^FL& wCD(ixO^S,Dtr^FL_)   = wSub(ixO^S,Dtr^FL_)&
                   *(cspeed(ixO^S)-vSub(ixO^S))/(cspeed(ixO^S)-lambdaCD(ixO^S)) \}
}
endwhere

!==== Magnetic field ====!
do iw =b1_,b^NC_
  where(abs(patchf(ixO^S))==1)
    ! case from eq 31
    wCD(ixO^S,iw) = whll(ixO^S,iw)
  endwhere
enddo 

!------- Momentum ------!
do iw =1, ^NC
  if(iw /= idims)then
     where(abs(patchf(ixO^S))==1)
       ! eq. 21 22
      wCD(ixO^S,m0_+iw)=(cspeed(ixO^S)*wSub(ixO^S,m0_+iw)-fSub(ixO^S,m0_+iw)  &
                       -wCD(ixO^S,b0_+idims)*wCD(ixO^S,b0_+iw)) /&
                       (cspeed(ixO^S)-lambdaCD(ixO^S))
     end where
  else
     where(abs(patchf(ixO^S))==1)
        ! eq. 20
        wCD(ixO^S,m0_+iw) =  wCD(ixO^S,rho_) * lambdaCD(ixO^S)
     endwhere
  endif
enddo

where(abs(patchf(ixO^S))==1)
{#IFDEF ENERGY
   VdotBCD(ixO^S) = ({^C&whll(ixO^S,m^C_)*whll(ixO^S,b^C_)+})/whll(ixO^S,rho_)
   ! Eq 17
   pCD(ixO^S)  = wsub(ixO^S,rho_)*(cspeed(ixO^S)-vSub(ixO^S))&
                  *(lambdaCD(ixO^S)-vSub(ixO^S))&
                  +fSub(ixO^S,m0_+idims)- wsub(ixO^S,m0_+idims)*vSub(ixO^S)&
                  + wCD(ixO^S,b0_+idims)**2.0d0
   ! Eq 31
   wCD(ixO^S,e_) = (cspeed(ixO^S) * wSub(ixO^S,e_) &
                 - fSub(ixO^S,e_) +lambdaCD(ixO^S)*pCD(ixO^S)&
                 -VdotBCD(ixO^S)*wCD(ixO^S,b0_+idims))&
                 /(cspeed(ixO^S)-lambdaCD(ixO^S))
}
end where

do iw=1,nwflux
 if(iw /= b0_+idims{#IFDEF GLM .and. iw/=psi_})then
   where(abs(patchf(ixO^S))==1)
       ! f_i=fsub+lambda (wCD-wSub)
       f(ixO^S,iw)=fsub(ixO^S,iw)+cspeed(ixO^S)*(wCD(ixO^S,iw)-wsub(ixO^S,iw))
   endwhere
 else
       f(ixO^S,iw)=zero
 end if
end do

return
end subroutine getwCD
!=============================================================================
! end module amrvacphys.mhdhllc.t 
!##############################################################################
