!##############################################################################
! module amrvacnul.roe.t - nul-routines  
!=============================================================================
subroutine average(wL,wR,x,ix^L,idims,wroe,workroe)

use mod_global_parameters

integer:: ix^L,idims,idir
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, intent(in)      :: x(ixG^T,1:ndim)
double precision, dimension(ixG^T,nworkroe):: workroe
!-----------------------------------------------------------------------------

call mpistop("Error:average: no roe solver implemented!")

end subroutine average
!=============================================================================
subroutine geteigenjump(wL,wR,wroe,x,ix^L,il,idims,smalla,a,jump,workroe)

use mod_global_parameters

integer:: ix^L,il,idims
double precision, dimension(ixG^T,nw):: wL,wR,wroe
double precision, intent(in)         :: x(ixG^T,1:ndim)
double precision, dimension(ixG^T)   :: smalla,a,jump
double precision, dimension(ixG^T,nworkroe) :: workroe
!-----------------------------------------------------------------------------

call mpistop("Error:geteigenjump: no roe solver implemented!")

end subroutine geteigenjump
!=============================================================================
subroutine rtimes(q,wroe,ix^L,iw,il,idims,rq,workroe)

use mod_global_parameters

integer::          ix^L,iw,il,idims
double precision:: wroe(ixG^T,nw)
double precision, dimension(ixG^T):: q,rq
double precision, dimension(ixG^T,nworkroe):: workroe
!-----------------------------------------------------------------------------

call mpistop("Error:rtimes: no roe solver implemented!")

end subroutine rtimes
!=============================================================================
! end module amrvacnul.roe.t
!##############################################################################
