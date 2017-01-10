!#############################################################################
! module amrvacusr - testkh3d

!INCLUDE:amrvacnul/specialini.t
INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
!INCLUDE:amrvacnul/specialsource.t
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
eqpar(gamma_)=5.0d0/3.0d0

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: rhodens,rholight,kx,kz,pint,dlin,vextra,sigma
logical :: patchw(ixG^T)

logical::          first
data first/.true./
!----------------------------------------------------------------------------

rhodens=10.0d0
rholight=1.0d0
! kx=2 pi
kx=4.0d0*dpi
{^IFTHREED
! kz=8 pi
kz=32d0*atan(one)
}

select case(iprob)
case(1)
vextra=0.0d0
case(2)
vextra=10.0d0
case(3)
vextra=100.0d0
endselect

if (first) then
   if (mype==0) then
      print *,'3D HD KH'
      print *,'  --assuming y ranging from 0-1!'
      print *,'  --density ratio:',rhodens/rholight
      print *,'  --kx:',kx
      print *,'  --kz:',kz
      print *,'  --vextra:',vextra
   end if
   first=.false.
end if

! pressure at interface
pint=2.5d0

dlin=0.025d0
sigma=0.00125d0

where(x(ixG^S,2)<0.25d0.or.x(ixG^S,2)>0.75d0)
   w(ixG^S,rho_)=rholight
   !w(ixG^S,tr1_)=0.0d0
elsewhere
   w(ixG^S,rho_)=rhodens
   !w(ixG^S,tr1_)=1.0d0
endwhere

w(ixG^S,p_)=pint
where((x(ixG^S,2)<=(0.25d0-dlin)).or.(x(ixG^S,2)>=(0.75d0+dlin)))
   w(ixG^S,v1_)=vextra-0.5d0
endwhere
where((x(ixG^S,2)>=(0.25d0+dlin)).and.(x(ixG^S,2)<=(0.75d0-dlin)))
   w(ixG^S,v1_)=vextra+0.5d0
endwhere
where((x(ixG^S,2)>(0.25d0-dlin)).and.(x(ixG^S,2)<(0.25d0+dlin)))
   w(ixG^S,v1_)=(x(ixG^S,2)-(0.25d0-dlin))* &
                ((vextra+0.5d0)-(vextra-0.5d0))/(2.0d0*dlin)+(vextra-0.5d0)
endwhere
where((x(ixG^S,2)>(0.75d0-dlin)).and.(x(ixG^S,2)<(0.75d0+dlin)))
   w(ixG^S,v1_)=(x(ixG^S,2)-(0.75d0-dlin))* &
                ((vextra-0.5d0)-(vextra+0.5d0))/(2.0d0*dlin)+(vextra+0.5d0)
endwhere

w(ixG^S,v2_)=0.01d0*dsin(kx*x(ixG^S,1))* &
     (dexp(-0.5d0*(x(ixG^S,2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixG^S,2)-0.75d0)**2/sigma))
{^IFTHREED
w(ixG^S,v3_)=0.1d0*dsin(kz*x(ixG^S,3))* &
     (dexp(-0.5d0*(x(ixG^S,2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixG^S,2)-0.75d0)**2/sigma))
}

patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)


end subroutine initonegrid_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
!-----------------------------------------------------------------------------


end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------


end subroutine getdt_special
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

include 'amrvacdef.f'

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------


end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

include 'amrvacdef.f'

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')

var(ixI^S) = zero

end subroutine specialvarforerrest

!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0

!=============================================================================
! amrvacusr.t.testkh3d
!=============================================================================
