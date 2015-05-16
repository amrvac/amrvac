!#############################################################################
! module amrvacphys/- nonlinear

INCLUDE:amrvacnul/addsource.t
INCLUDE:amrvacnul/hllc.t
INCLUDE:amrvacnul/getdt.t
INCLUDE:amrvacnul/getaux.t
INCLUDE:amrvacnul/roe.t
!=============================================================================
subroutine checkglobaldata

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

! nothing to check

end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

! For burgers/nonconvex nothing to do in terms of entropy fixes for Roe solver

end subroutine initglobaldata
!=============================================================================
subroutine conserve(ixI^L,ixO^L,w,x,patchw)

! Transform primitive variables into conservative ones

include 'amrvacdef.f'

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

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

! For the burgers/nonconvex equation primitive and conservative variables coincide

end subroutine primitive
!=============================================================================
subroutine e_to_rhos(ixI^L,ixO^L,w,x)

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

! No e and rhos available

end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixI^L,ixO^L,w,x)

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

! No e and rhos available

end subroutine rhos_to_e
!=============================================================================
subroutine getv(w,x,ixI^L,ixO^L,idim,v)

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, idim
double precision :: w(ixI^S,nw), v(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

select case(nint(eqpar(fluxtype_)))
 case(1)
   v(ixO^S)=w(ixO^S,rho_)
 case(2)
   v(ixO^S)=3.0d0*w(ixO^S,rho_)**2
 case default
    call mpistop('Undefined fluxtype: set eqpar to 1 or 2')
end select

end subroutine getv
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixI^L,ixO^L,idim,cmax,cmin,needcmin)

include 'amrvacdef.f'

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
subroutine getfluxforhllc(w,x,ixI^L,ixO^L,iw,idim,f,transport)

! There is nothing to add to the transport flux in the transport equation

include 'amrvacdef.f'

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

include 'amrvacdef.f'

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
subroutine addgeometry(qdt,ixI^L,ixO^L,wCT,w,x)

! Add geometrical source terms to w
! There are no geometrical source terms in the transport equation

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: qdt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine addgeometry
!=============================================================================
! end module amrvacphys/- nonlinear
!#############################################################################
