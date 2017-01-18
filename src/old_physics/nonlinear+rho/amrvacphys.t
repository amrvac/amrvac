!#############################################################################
! module amrvacphys/- nonlinear

INCLUDE:amrvacnul/addsource.t
INCLUDE:amrvacnul/hllc.t
INCLUDE:amrvacnul/getdt.t
INCLUDE:amrvacnul/getaux.t
INCLUDE:amrvacnul/roe.t
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

! For burgers/nonconvex nothing to do in terms of entropy fixes for Roe solver

end subroutine initglobaldata
!=============================================================================
subroutine conserve(ixI^L,ixO^L,w,x,patchw)

! Transform primitive variables into conservative ones

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
logical :: patchw(ixG^T)
!-----------------------------------------------------------------------------

! For the burgers/nonconvex equation primitive and conservative variables coincide

end subroutine conserve
!=============================================================================
subroutine primitive(ixI^L,ixO^L,w,x)

! Transform conservative variables into primitive ones

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

! For the burgers/nonconvex equation primitive and conservative variables coincide

end subroutine primitive
!=============================================================================
subroutine e_to_rhos(ixI^L,ixO^L,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

! No e and rhos available

end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixI^L,ixO^L,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

! No e and rhos available

end subroutine rhos_to_e
!=============================================================================
subroutine getv(w,x,ixI^L,ixO^L,idim,v)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, idim
double precision :: w(ixI^S,nw), v(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

v(ixO^S)=eqpar(v0_+idim)

end subroutine getv
!=============================================================================
subroutine getv2(w,x,ixI^L,ixO^L,idim,v)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, idim
double precision :: w(ixI^S,nw), v(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

v(ixO^S)=w(ixO^S,rho_)

end subroutine getv2
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixI^L,ixO^L,idim,cmax,cmin,needcmin)

use mod_global_parameters

logical :: new_cmax,needcmin
integer, intent(in) :: ixI^L, ixO^L, idim
double precision :: w(ixI^S,nw), cmax(ixG^T),cmin(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

call getv(w,x,ixI^L,ixO^L,idim,cmax)
if(needcmin)then
  cmin(ixO^S)=min(cmax(ixO^S),zero)
  cmax(ixO^S)=max(cmax(ixO^S),zero)
else
  cmax(ixO^S)=dabs(cmax(ixO^S))
endif 

end subroutine getcmax
!=============================================================================
subroutine getcmax_multi(new_cmax,w,x,ixI^L,ixO^L,idim,cmax,cmin,needcmin,phys)

use mod_global_parameters

logical :: new_cmax,needcmin
integer, intent(in) :: ixI^L, ixO^L, idim, phys
double precision :: w(ixI^S,nw), cmax(ixG^T),cmin(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

select case (phys)
case(1)
   call getv(w,x,ixI^L,ixO^L,idim,cmax)
case(2)
   call getv2(w,x,ixI^L,ixO^L,idim,cmax)
end select

if(needcmin)then
   cmin(ixO^S)=min(cmax(ixO^S),zero)
   cmax(ixO^S)=max(cmax(ixO^S),zero)
else
   cmax(ixO^S)=dabs(cmax(ixO^S))
endif

end subroutine getcmax_multi
!=============================================================================
subroutine getfluxforhllc(w,x,ixI^L,ixO^L,iw,idim,f,transport)

! There is nothing to add to the transport flux in the transport equation

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw, idim
double precision :: w(ixI^S,1:nw), f(ixG^T,1:nwflux)
double precision, intent(in)    :: x(ixI^S,1:ndim)
logical :: transport
!-----------------------------------------------------------------------------

select case(nint(eqpar(fluxtype_)))
 case(1)
    f(ixO^S,iw)=half*w(ixO^S,rho_)**2
 case(2)
    f(ixO^S,iw)=w(ixO^S,rho_)**3
 case default
    call mpistop('Undefined fluxtype: set eqpar to 1 or 2')
end select

transport=.false.

end subroutine getfluxforhllc
!=============================================================================
subroutine getflux(w,x,ixI^L,ixO^L,iw,idim,f,transport)

! There is nothing to add to the transport flux in the transport equation

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw, idim
double precision :: w(ixI^S,1:nw), f(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
logical :: transport
!-----------------------------------------------------------------------------

select case(nint(eqpar(fluxtype_)))
 case(1)
    f(ixO^S)=half*w(ixO^S,rho_)**2
 case(2)
    f(ixO^S)=w(ixO^S,rho_)**3
 case default
    call mpistop('Undefined fluxtype: set eqpar to 1 or 2')
end select
transport=.false.

end subroutine getflux
!=============================================================================
subroutine getflux_multi(w,x,ixI^L,ixO^L,iw,idim,f,transport,phys)

! There is nothing to add to the transport flux in the transport equation

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw, idim, phys
double precision :: w(ixI^S,1:nw), f(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
logical :: transport
!-----------------------------------------------------------------------------

select case(phys)
case(1)
   f(ixO^S) = w(ixO^S,rho_)*eqpar(v0_+idim)
case(2)
   f(ixO^S)=half*w(ixO^S,rho_)**2
end select

transport=.false.

end subroutine getflux_multi
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
! end module amrvacphys/- nonlinear
!#############################################################################
