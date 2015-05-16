!##############################################################################
! module amrvacnul.hllc.t - nul-routines  
!=============================================================================
subroutine diffuse_hllcd(ixI^L,ixO^L,idims,wLC,wRC,fLC,fRC,patchf)

! when method is hllcd or hllcd1 then: 

! this subroutine is to impose enforce regions where we AVOID HLLC
! and use TVDLF instead: this is achieved by setting patchf to 4 in
! certain regions. An additional input parameter is nxdiffusehllc
! which sets the size of the fallback region.

! This nul version enforces TVDLF everywhere!!!

include 'amrvacdef.f'

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixG^T,1:nw), intent(in)      :: wRC,wLC
double precision, dimension(ixG^T,1:nwflux),intent(in) :: fLC, fRC

integer         , dimension(ixG^T), intent(inout)        :: patchf
!-----------------------------------

! enforce TVDLF everywhere
patchf(ixO^S) = 4

end subroutine diffuse_hllcd
!=============================================================================
subroutine getlCD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixI^L,ixO^L, &
                  whll,Fhll,lambdaCD,patchf)

! Calculate lambda at CD and set the patchf to know the orientation
! of the riemann fan and decide on the flux choice
! We also compute here the HLL flux and w value, for fallback strategy

! In this nul version, we simply compute nothing and ensure TVDLF fallback

include 'amrvacdef.f'

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixG^T,1:nw), intent(in)      :: wLC,wRC
double precision, dimension(ixG^T,1:nwflux), intent(in):: fLC,fRC
double precision, dimension(ixG^T), intent(in)           :: cmax,cmin

integer         , dimension(ixG^T), intent(inout)        :: patchf

double precision, dimension(ixG^T,1:nwflux), intent(out) :: Fhll,whll
double precision, dimension(ixG^T), intent(out)            :: lambdaCD
!--------------------------------------------

! Next must normally be computed
Fhll(ixO^S,1:nwflux)=zero
whll(ixO^S,1:nwflux)=zero
lambdaCD(ixO^S)=zero

! this actually ensures fallback to TVDLF
patchf(ixO^S)=4

return
end subroutine getlCD
!=============================================================================
subroutine getwCD(wLC,wRC,whll,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaCD,cmin,cmax,&
                  ixI^L,ixO^L,idims,f)

! compute the intermediate state U*
! only needed where patchf=-1/1

! This nul version simply nullifies all values

include 'amrvacdef.f'

integer, intent(in)                                      :: ixI^L,ixO^L,idims
double precision, dimension(ixG^T,1:nw), intent(in)      :: wRC,wLC
double precision, dimension(ixG^T,1:nwflux), intent(in):: whll, Fhll
double precision, dimension(ixG^T), intent(in)           :: vRC, vLC,lambdaCD
double precision, dimension(ixG^T), intent(in)           :: cmax,cmin
double precision, dimension(ixG^T,1:nwflux), intent(in):: fRC,fLC

integer         , dimension(ixG^T), intent(in)           :: patchf

double precision, dimension(ixG^T,1:nwflux),intent(out) :: f
!--------------------------------------------

! Next must normally be computed
f(ixO^S,1:nwflux)  = zero

return
end subroutine getwCD
!=============================================================================
! end module amrvacnul.hllc.t
!##############################################################################
