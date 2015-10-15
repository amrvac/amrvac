{#IFDEF HALL
!=============================================================================
subroutine getvh(w,x,ixI^L,ixO^L,vh)
include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(inout) :: vh(ixI^S,1:3)
!.. local ..
integer          :: idir, idirmin
double precision :: current(ixI^S,7-2*ndir:3)

! Calculate current density and idirmin
call getcurrent(w,ixI^L,ixO^L,idirmin,current)
vh(ixI^S,1:3) = zero
vh(ixO^S,idirmin:3) = - eqpar(etah_)*current(ixO^S,idirmin:3)
do idir = idirmin, 3
   vh(ixO^S,idir) = vh(ixO^S,idir)/w(ixO^S,rho_)
end do

end subroutine getvh
!=============================================================================
subroutine getdthall(w,x,ixI^L,ixO^L,dx^D,dthall)
include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in)    :: dx^D
double precision, intent(in)    :: w(ixI^S,1:nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(out)   :: dthall
!.. local ..
double precision :: dxarr(ndim)
double precision :: bmag(ixI^S)

dthall=bigdouble

! because we have that in cmax now:
return

^D&dxarr(^D)=dx^D;

if (.not. B0field) then
   bmag(ixO^S)=sqrt({^C& w(ixO^S,b^C_)**2 +})
else
   bmag(ixO^S)=sqrt({^C& (w(ixO^S,b^C_)+myB0%w(ixO^S,^C))**2 +})
end if

dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(eqpar(etah_)*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
end subroutine getdthall
!=============================================================================
}
