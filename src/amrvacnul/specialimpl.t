!=============================================================================
subroutine specialsource_impl(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), wCT(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine specialsource_impl
!=============================================================================
subroutine getdt_impl(w,ixG^L,ix^L,dtnew,dx^D,x)

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
! note that depending on strictsmall etc, w values may change
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

end subroutine getdt_impl
!=============================================================================
