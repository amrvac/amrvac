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

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
eqpar(gamma_)=1.4d0
eqpar(rhoj_)=eqpar(gamma_)
eqpar(eta_)=10.0d0
eqpar(vj_)=800.0d0


end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
!----------------------------------------------------------------------------

{^IFONED call mpistop("This is a multi-D HD problem") }

where(dabs(x(ix^S,1))<0.05d0.and.x(ix^S,2)<0.00d0)
   w(ix^S,rho_)=eqpar(rhoj_)
   w(ix^S,m1_)=zero
   w(ix^S,m2_)=eqpar(rhoj_)*eqpar(vj_)
   w(ix^S,e_)=one/(eqpar(gamma_)-one)+0.5d0*eqpar(rhoj_)*eqpar(vj_)**2.0d0
else where
   w(ix^S,rho_) = eqpar(rhoj_)/eqpar(eta_)
   w(ix^S,e_) = one/(eqpar(gamma_)-one)
   w(ix^S,m1_) = zero
   w(ix^S,m2_) = zero
end where

end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixB^L,iw,iB,w,x)

! special boundary types, user defined

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ixB^L, iw, iB
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
integer :: ixI^L, ix2
!----------------------------------------------------------------------------
ixImin^DD=ixBmin^DD;
ixImax^DD=ixBmin^D-1+dixB^D%ixImax^DD=ixBmax^DD;
! Outflow:
    do ix2=ixImin2,ixImax2
   w(ixImin1:ixImax1,ix2,rho_) = w(ixImin1:ixImax1,ixImax2+1,rho_) 
   w(ixImin1:ixImax1,ix2,e_)   = w(ixImin1:ixImax1,ixImax2+1,e_) 
   w(ixImin1:ixImax1,ix2,m1_)  = w(ixImin1:ixImax1,ixImax2+1,m1_)
   w(ixImin1:ixImax1,ix2,m2_)  = w(ixImin1:ixImax1,ixImax2+1,m2_)
   end do
where(dabs(x(ixI^S,1))<0.05d0)
   w(ixI^S,rho_)=eqpar(rhoj_)
   w(ixI^S,m1_)=zero
   w(ixI^S,m2_)=eqpar(rhoj_)*eqpar(vj_)
   w(ixI^S,e_)=one/(eqpar(gamma_)-one)+0.5d0*eqpar(rhoj_)*eqpar(vj_)**2.0d0
else where
! Reflective:
!   w(ixI^S,rho_) = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,rho_) 
!   w(ixI^S,e_) = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,e_) 
!   w(ixI^S,m1_) = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,m1_)
!   w(ixI^S,m2_) =-w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,m2_)
end where

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

include 'amrvacdef.f'

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

include 'amrvacdef.f'

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

include 'amrvacdef.f'

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

include 'amrvacdef.f'

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

include 'amrvacdef.f'

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
