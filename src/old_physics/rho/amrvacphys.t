!#############################################################################
! module amrvacphys/- rho

INCLUDE:amrvacnul/addsource.t
INCLUDE:amrvacnul/hllc.t
INCLUDE:amrvacnul/getdt.t
INCLUDE:amrvacnul/getaux.t
!=============================================================================
subroutine checkglobaldata

use mod_global_parameters
!-----------------------------------------------------------------------------

! nothing to check

end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

use mod_global_parameters
!-----------------------------------------------------------------------------

! For transport nothing to do in terms of entropy fixes for Roe solver

end subroutine initglobaldata
!=============================================================================
subroutine conserve(ixI^L,ixO^L,w,x,patchw)

! Transform primitive variables into conservative ones

use mod_global_parameters

integer, intent(in)           :: ixI^L, ixO^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision              :: w(ixI^S,nw)
logical                       :: patchw(ixG^T)
!-----------------------------------------------------------------------------

! For the transport equation primitive and conservative variables coincide

end subroutine conserve
!=============================================================================
subroutine primitive(ixI^L,ixO^L,w,x)

! Transform conservative variables into primitive ones

use mod_global_parameters

integer, intent(in)           :: ixI^L, ixO^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision              :: w(ixI^S,nw)
!-----------------------------------------------------------------------------

! For the transport equation primitive and conservative variables coincide

end subroutine primitive
!=============================================================================
subroutine e_to_rhos(ixI^L,ixO^L,w,x)

use mod_global_parameters

integer, intent(in)           :: ixI^L, ixO^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision              :: w(ixI^S,nw)
!-----------------------------------------------------------------------------

! No e and rhos available

end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixI^L,ixO^L,w,x)

use mod_global_parameters

integer, intent(in)           :: ixI^L, ixO^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision              :: w(ixI^S,nw)
!-----------------------------------------------------------------------------

! No e and rhos available

end subroutine rhos_to_e
!=============================================================================
subroutine getv(w,x,ixI^L,ixO^L,idim,v)

use mod_global_parameters

integer, intent(in)           :: ixI^L, ixO^L, idim
double precision, intent(in)  ::  w(ixI^S,nw), x(ixI^S,1:ndim)
double precision, intent(out) :: v(ixG^T)
!-----------------------------------------------------------------------------

v(ixO^S) = eqpar(v0_+idim)

end subroutine getv
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixI^L,ixO^L,idim,cmax,cmin,needcmin)

use mod_global_parameters

logical :: new_cmax,needcmin
integer, intent(in) :: ixI^L, ixO^L, idim
double precision :: w(ixI^S,nw), x(ixI^S,1:ndim), cmax(ixG^T),cmin(ixG^T)
!-----------------------------------------------------------------------------
call getv(w,x,ixI^L,ixO^L,idim,cmax)

if(needcmin)then
  cmin(ixO^S)=min(cmax(ixO^S),zero)
  cmax(ixO^S)=max(cmax(ixO^S),zero)
else
  cmax(ixO^S)=abs(cmax(ixO^S))
endif 

end subroutine getcmax
!=============================================================================
subroutine getfluxforhllc(w,x,ixI^L,ixO^L,iw,idim,f,transport)

! There is nothing to add to the transport flux in the transport equation

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw, idim
double precision :: w(ixI^S,1:nw), x(ixI^S,1:ndim), f(ixG^T,1:nwflux)
logical :: transport
!-----------------------------------------------------------------------------

f(ixO^S,iw)=zero
transport=.true.

end subroutine getfluxforhllc
!=============================================================================
subroutine getflux(w,x,ixI^L,ixO^L,iw,idim,f,transport)

! There is nothing to add to the transport flux in the transport equation

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw, idim
double precision :: w(ixI^S,1:nw), x(ixI^S,1:ndim), f(ixG^T)
logical :: transport
double precision :: v(ixI^S)
!-----------------------------------------------------------------------------
!call getv(w,x,ixI^L,ixO^L,idim,v)
f(ixO^S)=zero
transport=.true.

end subroutine getflux
!=============================================================================
subroutine average(wL,wR,x,ix^L,idim,wroe,workroe)

! Calculate the algebraic mean of primitive variables, assignment:
! rho -> rho

use mod_global_parameters

integer, intent(in) :: ix^L, idim
double precision, dimension(ixG^T,nw) :: wL, wR, wroe
double precision, dimension(ixG^T,nworkroe) :: workroe
double precision, dimension(ixG^T,1:ndim):: x
!-----------------------------------------------------------------------------

wroe(ix^S,rho_)=half*(wL(ix^S,rho_)+wR(ix^S,rho_))

end subroutine average
!=============================================================================
subroutine geteigenjump(wL,wR,wC,x,ix^L,il,idim,smalla,a,jump,workroe)

! Calculate the characteristic speed a and the jump in the
! characteristic variable in the idim direction within ixL.
! For a scalar equation the characteristic and conservative variables coincide
! The characteristic speed is just the velocity, but it should be averaged
! for the cell interfaces

use mod_global_parameters

integer:: ix^L,jx^L,ixC^L,il,idim
double precision, dimension(ixG^T,nw):: wL,wR,wC
double precision, dimension(ixG^T)   :: smalla,a,jump,v
double precision, dimension(ixG^T,nworkroe)   :: workroe
double precision, dimension(ixG^T,1:ndim):: x
!-----------------------------------------------------------------------------

jx^L=ix^L+kr(idim,^D);
ixCmin^D=ixmin^D; ixCmax^D=jxmax^D;

! No entropy fix
smalla(ix^S)= -one
! The velocity is independent of w in the transport equation,
! but it may depend on the location
call getv(wL,x,ixG^LL,ixC^L,idim,v)

a(ix^S)=(v(jx^S)+v(ix^S))/2

jump(ix^S)=wR(ix^S,rho_)-wL(ix^S,rho_)

end subroutine geteigenjump
!=============================================================================
subroutine rtimes(q,w,ix^L,iw,il,idim,rq,workroe)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wC.
! For a scalar equation the R matrix is unity

use mod_global_parameters

integer::          ix^L,iw,il,idim
double precision:: w(ixG^T,nw),q(ixG^T),rq(ixG^T)
double precision, dimension(ixG^T,nworkroe)   :: workroe
!-----------------------------------------------------------------------------

rq(ix^S)=q(ix^S)

end subroutine rtimes
!=============================================================================
subroutine addgeometry(qdt,ixI^L,ixO^L,wCT,w,x)

! Add geometrical source terms to w
! There are no geometrical source terms in the transport equation

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: qdt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine addgeometry
!=============================================================================
! end module amrvacphys/- rho
!#############################################################################
