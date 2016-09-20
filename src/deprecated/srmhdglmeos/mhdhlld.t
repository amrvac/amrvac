!##############################################################################
! module amrvacphys.mhdhlld.t 
!=============================================================================
! Aim : HLLD
! Description :
!=============================================================================
subroutine diffuse_hlldd(ixI^L,ixO^L,idims,wLC,wRC,fLC,fRC,patchf)

! when method is hllcd or hllcd1 then: 

! this subroutine is to enforce regions where we AVOID HLLC
! and use TVDLF instead: this is achieved by setting patchf to 4 in
! certain regions. An additional input parameter is nxdiffusehllc
! which sets the size of the fallback region.

use mod_global_parameters

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

end subroutine diffuse_hlldd
!=============================================================================
subroutine getlD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixI^L,ixO^L, &
                  whll,Fhll,lambdaD,patchf)

! Calculate lambda at CD and set the patchf to know the orientation
! of the riemann fan and decide on the flux choice
! We also compute here the HLL flux and w value, for fallback strategy

use mod_global_parameters

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixG^T,1:nw), intent(in)      :: wLC,wRC
double precision, dimension(ixG^T,1:nwflux), intent(in):: fLC,fRC
double precision, dimension(ixG^T), intent(in)           :: cmax,cmin

integer         , dimension(ixG^T), intent(inout)        :: patchf

double precision, dimension(ixG^T,1:nwflux), intent(out) :: Fhll,whll
double precision, dimension(ixG^T,-1:1), intent(out)     :: lambdaD
!--------------------------------------------

call mpistop("hlld not yet implemeted")
end subroutine getlD
!=============================================================================
subroutine getwD(wLC,wRC,whll,x,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaD,cmin,cmax,&
                  ixI^L,ixO^L,idims,f)

! compute the intermediate state U*
! only needed where patchf=-1/1

! reference Li S., JCP, 203, 2005, 344-357
! reference T. Miyoski, Kusano JCP, 2008, 2005.

use mod_global_parameters

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixG^T,1:nw), intent(in)      :: wRC,wLC
double precision, dimension(ixG^T,1:ndim), intent(in)    :: x
double precision, dimension(ixG^T,1:nwflux), intent(in):: whll, Fhll
double precision, dimension(ixG^T), intent(in)           :: vRC, vLC
double precision, dimension(ixG^T,-1:1), intent(in)      ::lambdaD
double precision, dimension(ixG^T), intent(in)           :: cmax,cmin
double precision, dimension(ixG^T,1:nwflux), intent(in)  :: fRC,fLC
double precision, dimension(ixG^T,1:nwflux),intent(out)  :: f
integer         , dimension(ixG^T), intent(in)           :: patchf
!--------------------------------------------

call mpistop("hlld not yet implemeted")
end subroutine getwD
!=============================================================================
! end module amrvacphys.mhdhlld.t 
!##############################################################################
