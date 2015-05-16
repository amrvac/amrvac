!=============================================================================
!HLLC for SRMHD (ZM)
!=============================================================================
subroutine diffuse_hllcd(ixI^L,ixO^L,idims,wLC,wRC,fLC,fRC,patchf)

! when method is hllcd or hllcd1 then:

! this subroutine is to impose enforce regions where we AVOID HLLC 
! and use TVDLF instead: this is achieved by setting patchf to 4 in
! certain regions. An additional input parameter is nxdiffusehllc
! which sets the size of the fallback region.

include 'amrvacdef.f'

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixI^S,1:nw), intent(in)      :: wRC,wLC
double precision, dimension(ixI^S,1:nwflux),intent(in) :: fLC, fRC

integer         , dimension(ixG^T), intent(inout)        :: patchf

integer                                           :: ixOO^D,TxOO^L
!-----------------------------------

! In a user-controlled region around any point with flux sign change between
! left and right, ensure fallback to TVDLF
call mpistop('hllc not implemented in srmhdiso')

end subroutine diffuse_hllcd
!=============================================================================
subroutine getlCD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixI^L,ixO^L,&
                  whll,Fhll,lambdaCD,patchf)
  
! Calculate lambda at CD and set the patchf to know the orientation
! of the riemann fan and decide on the flux choice
! We also compute here the HLL flux and w value, for fallback strategy

include 'amrvacdef.f'
  
integer, intent(in)                                        :: ixI^L,ixO^L,idims
double precision, dimension(ixI^S,1:nw), intent(in)        :: wLC,wRC
double precision, dimension(ixG^T,1:nwflux), intent(in)  :: fLC,fRC
double precision, dimension(ixG^T), intent(in)             :: cmax,cmin

integer         , dimension(ixG^T), intent(inout)          :: patchf

double precision, dimension(ixG^T), intent(inout)            :: lambdaCD
double precision, dimension(ixG^T,1:nwflux), intent(inout) :: Fhll,whll

double precision, dimension(ixG^T)     :: Aco,Bco,Cco,Delta, BpdotfBp, Bp2,fBp2
logical         , dimension(ixG^T)     :: Cond_patchf, Cond_Bidimhll
double precision                       :: Epsilon
integer                                :: iw
!--------------------------------------------
! on entry, patch is preset to contain values from -2,1,2,4
!      -2: take left flux, no computation here
!      +2: take right flux, no computation here
!      +4: take TVDLF flux, no computation here
!       1: compute the characteristic speed for the CD
call mpistop('hllc not implemented in srmhdiso')

end subroutine getlCD
!=============================================================================
subroutine getwCD(wLC,wRC,whll,vLC,vRC,fRC,fLC,Fhll,patchf,&
                  lambdaCD,cmin,cmax,ixI^L,ixO^L,idims,f)
  

! compute the intermediate state U*
! only needed where patchf=-1/1

! For SRMHD compute D*, S*, tau* in wCD and v* in vCD  and p* in PCD

include 'amrvacdef.f'

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixI^S,1:nw), intent(in)      :: wLC,wRC
double precision, dimension(ixG^T,1:nwflux), intent(in)  :: fRC,fLC
double precision, dimension(ixG^T,1:nwflux), intent(in)  :: whll,Fhll
double precision, dimension(ixG^T), intent(in)           :: vRC,vLC,lambdaCD
double precision, dimension(ixG^T), intent(in)           :: cmax,cmin
double precision, dimension(ixG^T,1:nwflux), intent(inout) :: f

integer         , dimension(ixG^T), intent(inout)        :: patchf

double precision, dimension(ixG^T,s1_:s^NC_):: vCD
double precision, dimension(ixG^T,1:nw)     :: wCD,wSub
double precision, dimension(ixG^T,1:nwflux) :: fSub
double precision, dimension(ixG^T)          :: vSub,cspeed,pCD,VdotB,Ratio_CD

logical         , dimension(ixG^T)           :: Cond_Bidimhll
integer                                      :: iw
!--------------------------------------------
call mpistop('hllc not implemented in srmhdiso')

end subroutine getwCD
!============================================================================
!End HLLC for SRMHD 
!=============================================================================
