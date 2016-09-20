!=============================================================================
subroutine getaux(clipping,w,x,ixI^L,ixO^L,subname)

! Calculate **LOCAL** auxiliary variables ixO^L from non-auxiliary entries in w
! clipping can be set to .true. to e.g. correct unphysical pressures, 
! densities, v>c,  etc.

use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,nw)
logical                         :: clipping
character(len=*)                :: subname
!-----------------------------------------------------------------------------

end subroutine getaux
!=============================================================================
subroutine checkw(checkprimitive,ixI^L,ixO^L,w,flag)
!
! Check if the conserved (primitive if checkprimitive=T) variables 
! in w represent physical values.
! Where this requirement is violated, flag is set to false
! Do *not* use the auxiliary variables  (iw>nw-nwaux-nwextra)
! Because this routine is meant to check state before calculating 
! the auxiliaries, they might not be up to date.
!
use mod_global_parameters
  
logical:: checkprimitive
integer, intent(in) :: ixI^L, ixO^L
double precision:: w(ixI^S,nw)
logical:: flag(ixG^T)
!-----------------------------------------

flag(ixO^S)=.true.

end subroutine checkw
!=============================================================================
subroutine ppmflatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)

! This nul version simply nullifies all values

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L,ixL^L,ixR^L
double precision, intent(in)  :: w(ixI^S,1:nw),d2w(ixG^T,1:nwflux)

double precision              :: drho(ixG^T),dp(ixG^T)
!-----------------------------------------------------------------------------
drho(ixO^S)=zero
dp(ixO^S)=zero

end subroutine ppmflatcd
!=============================================================================
subroutine ppmflatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)

! This nul version simply nullifies all values

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
integer, intent(in)           :: idims
double precision, intent(in)  :: w(ixI^S,nw)

double precision              :: drho(ixG^T),dp(ixG^T),dv(ixG^T)
!-----------------------------------------------------------------------------
drho(ixO^S)=zero
dp(ixO^S)=zero
dv(ixO^S)=zero

end subroutine ppmflatsh
!=============================================================================
