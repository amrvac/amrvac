!=============================================================================
! INCLUDE:amrvacnul/specialini.t
INCLUDE:amrvacnul/speciallog.t
!INCLUDE:amrvacnul/specialbound.t
!INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------
eqpar(gamma_)    = 5./3.
{#IFDEF EPSINF
eqpar(epsfloor_)    = 1.0d-6
eqpar(rho0floor_)   = 1.0d-6
eqpar(rho1e_)       = 1.0d-6
}
! Jet velocity:
eqpar(vmax_)        =  0.95d0
! Uniform density and pressure:
eqpar(p0_)            =  0.1d0
eqpar(rho0_)         =  1.0d0
end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
!..local..
logical:: patchw(ixG^T)
!----------------------------------------------------------------------------
{^IFONED call mpistop("This is a multi-D SRHD problem") }
if(useprimitiveRel)then
   call mpistop('only implemented for useprimitiveRel=F')
endif


where(dabs(x(ix^S,1))<0.05d0.and.x(ix^S,2)<0.00d0)
   w(ix^S,rho_)=eqpar(rho0_)
   w(ix^S,u1_)=zero
   w(ix^S,u2_)=eqpar(vmax_)
   w(ix^S,pp_)=eqpar(p0_)
{#IFDEF TRACER
   {^FL&w(ix^S,tr^FL_)=one\}
}
{#IFDEF EPSINF
   w(ix^S,epsinf_)=one
   w(ix^S,ne_)=one
   w(ix^S,ne0_)=one
}
else where
   w(ix^S,rho_)=eqpar(rho0_)
   w(ix^S,u1_)=zero
   w(ix^S,u2_)=zero
   w(ix^S,pp_)=eqpar(p0_)
{#IFDEF EPSINF
   w(ix^S,epsinf_)=eqpar(epsfloor_)*10.
   w(ix^S,ne_)=eqpar(rho1e_)*10.
   w(ix^S,ne0_)=eqpar(rho1e_)*10.
}
end where


patchw(ix^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)
end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixB^L,iw,iB,w,x)

! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixG^L, ixB^L, iw, iB
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
integer :: ixI^L, ix2
!..local..
logical:: patchw(ixG^T)
!----------------------------------------------------------------------------
patchw(ixG^T) = .true.
!----------------------------------------------------------------------------

ixImin^DD=ixBmin^DD;
ixImax^DD=ixBmin^D-1+dixB^D%ixImax^DD=ixBmax^DD;
! Outflow:
    do ix2=ixImin2,ixImax2
   w(ixImin1:ixImax1,ix2,rho_) = w(ixImin1:ixImax1,ixImax2+1,rho_) 
   w(ixImin1:ixImax1,ix2,pp_)   = w(ixImin1:ixImax1,ixImax2+1,pp_) 
   w(ixImin1:ixImax1,ix2,u1_)  = w(ixImin1:ixImax1,ixImax2+1,u1_)
   w(ixImin1:ixImax1,ix2,u2_)  = w(ixImin1:ixImax1,ixImax2+1,u2_)
{#IFDEF EPSINF
   w(ixImin1:ixImax1,ix2,epsinf_)   = w(ixImin1:ixImax1,ixImax2+1,epsinf_) 
   w(ixImin1:ixImax1,ix2,ne_)  = w(ixImin1:ixImax1,ixImax2+1,ne_)
   w(ixImin1:ixImax1,ix2,ne0_)  = w(ixImin1:ixImax1,ixImax2+1,ne0_)
}
{#IFDEF TRACER
   {^FL&w(ixImin1:ixImax1,ix2,tr^FL_)=w(ixImin1:ixImax1,ixImax2+1,tr^FL_)\}
}
   end do
where(dabs(x(ixI^S,1))<0.05d0)
patchw(ixI^S)=.false.
   w(ixI^S,rho_)=eqpar(rho0_)
   w(ixI^S,u1_)=zero
   w(ixI^S,u2_)=eqpar(vmax_)
   w(ixI^S,pp_)=eqpar(p0_)
{#IFDEF EPSINF
   w(ixI^S,epsinf_)=one
   w(ixI^S,ne_)=one
   w(ixI^S,ne0_)=one
}
{#IFDEF TRACER
   {^FL&w(ixI^S,tr^FL_)=one\}
}
else where
! Reflective:
!   w(ixI^S,rho_) = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,rho_) 
!   w(ixI^S,e_) = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,e_) 
!   w(ixI^S,m1_) = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,m1_)
!   w(ixI^S,m2_) =-w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,m2_)
end where
  

call conserve(ixG^L,ixI^L,w,x,patchw)
end subroutine specialbound_usr
!=============================================================================
subroutine bc_int(qt,ixG^L,ixO^L,w,x)

! internal boundary, user defined
!
! This subroutine can be used to artificially overwrite ALL conservative 
! variables in a user-selected region of the mesh, and thereby act as
! an internal boundary region. It is called just before external (ghost cell)
! boundary regions will be set by the BC selection. Here, you could e.g. 
! want to introduce an extra variable (nwextra, to be distinguished from nwaux)
! which can be used to identify the internal boundary region location.
! Its effect should always be local as it acts on the mesh.
!

use mod_global_parameters

integer, intent(in) :: ixG^L,ixO^L
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

end subroutine bc_int
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

! integer :: iw
! double precision :: s(ixG^T)
!-----------------------------------------------------------------------------

! do iw= iw^LIM
!    select case(iw)
!    case(m1_)
!       ! The source is based on the time centered wCT
!       call getmyforce(wCT,ixO^L,s)
!       w(ixO^S,m1_)=w(ixO^S,m1_) + qdt*s(ixO^S)
!    case(e_)
!       call getmyheating(wCT,ixO^L,s)
!       w(ixO^S,e_) =w(ixO^S,e_)  + qdt*s(ixO^S)
!    end select
! end do

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew=bigdouble

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_global_parameters

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

!  eta(ix^S)=...

call mpistop("specialeta is not defined")

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

! you must set consistent values for integers refine/coarsen:

! refine = -1 enforce to not refine
! refine =  0 doesn't enforce anything
! refine =  1 enforce refinement

! coarsen = -1 enforce to not coarsen
! coarsen =  0 doesn't enforce anything
! coarsen =  1 enforce coarsen

use mod_global_parameters

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------

! e.g. refine for negative first coordinate x < 0 as
!
if (minval(dabs(x(ix^S,1))) < 0.1.and.minval(dabs(x(ix^S,2))) < 0.1) refine=1

end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

use mod_global_parameters

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

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0
!=============================================================================
!=============================================================================
! amrvacusr.t.liska
!=============================================================================
