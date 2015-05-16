!#############################################################################
! module amrvacusr - testhdrt
! setup.pl -d=33 -phi=0 -z=0 -g=16,16,16 -p=hd -eos=default -nf=0 -ndust=0 -u=nul -arch=default
! play around with the subroutine initonegrid_usr which setups the initial value

INCLUDE:amrvacmodules/gravity.t
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
eqpar(grav1_)=zero
eqpar(grav2_)=-one
eqpar(grav3_)=zero
end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

integer :: ix^D

double precision:: y0,epsilon,rhodens,rholight,kx,kz,pint,dely
logical::          first
data first/.true./
!----------------------------------------------------------------------------

! the location of demarcation line`
y0=0.8d0

! density of two types
rhodens=one
rholight=0.1d0

! setup the perturbation
epsilon=0.05d0
! kx=2 pi
kx=8.0d0*atan(one)
! kz=8 pi
kz=32d0*atan(one)

! print out the info
if (first) then
   if (mype==0) then
      print *,'HD Rayleigh Taylor problem'
      print *,'  --assuming y ranging from 0-1!'
      print *,'  --interface y0-epsilon:',y0,epsilon
      print *,'  --density ratio:',rhodens/rholight
      print *,'  --kx:',kx
      print *,'  --kz:',kz
   end if
   first=.false.
end if

! initialize the density
if(kx*kz/=zero)then
   where(x(ixG^S,2)>y0+epsilon*sin(kx*x(ixG^S,1))*sin(kz*x(ixG^S,3)))
      w(ixG^S,rho_)=rhodens
   elsewhere
      w(ixG^S,rho_)=rholight
   endwhere
else
   if(kx==zero)then
      where(x(ixG^S,2)>y0+epsilon*sin(kz*x(ixG^S,3)))
         w(ixG^S,rho_)=rhodens
      elsewhere
         w(ixG^S,rho_)=rholight
      endwhere
   else
      where(x(ixG^S,2)>y0+epsilon*sin(kx*x(ixG^S,1)))
         w(ixG^S,rho_)=rhodens
      elsewhere
         w(ixG^S,rho_)=rholight
      endwhere
   endif 
endif

! set all velocity to zero
{^C&w(ixG^S,m^C_)=zero \}

! pressure at interface
pint=one
{#IFDEF ENERGY
w(ixG^S,e_)=pint-w(ixG^S,rho_)*(x(ixG^S,2)-y0)
w(ixG^S,e_)=w(ixG^S,e_)/(eqpar(gamma_)-one)
}
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

call addsource_grav(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

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

call getdt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

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
! amrvacusr.t.testhdrt
!=============================================================================
