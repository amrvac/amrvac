!=============================================================================
! amrvacusr.t.friedrich
!=============================================================================
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
double precision:: bval
!----------------------------------------------------------------------------
eqpar(gamma_)=5.0d0/3.0d0

select case(iprob)
 case(1)
   bval=one/0.8d0
 case(2)
   bval=one
 case(3)
   bval=one/1.2d0
 case default
   call mpistop('c/b set with iprob=1,2,3, change in parfile or in user file')
endselect

eqpar(bxini_) =bval 
eqpar(dvz_) = 0.1d0
eqpar(dp_)  = 0.1d0

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision  :: rpert, dp, dvz, pr0, rho0, Bx0
logical :: patchw(ixG^T)
logical, save :: first=.true.
!----------------------------------------------------------------------------
patchw(ixG^S)=.false.

rpert=max(0.001d0,2.0*maxval(dx(1:ndim,mxnest)))
if(mype==0.and.first)then
 print *,'2.5D MHD Friedrich diagram test with ratio sound speed to alfven speed =',one/eqpar(bxini_)
 print *,'radius for perturbed zone=',rpert
 first=.false.
endif
pr0=one/eqpar(gamma_)
rho0=one

^C&w(ix^S,v^C_)=zero;
^C&w(ix^S,b^C_)=zero;
w(ix^S,rho_) = rho0
w(ix^S,pp_) = pr0
w(ix^S,b1_) = eqpar(bxini_)


where ({^D&x(ix^S,^D)**2.0d0+} <=rpert**2.0d0)
   w(ix^S,pp_) = pr0+eqpar(dp_)
   w(ix^S,v3_) = eqpar(dvz_)
endwhere

call conserve(ixG^L,ix^L,w,x,patchw)

end subroutine initonegrid_usr
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

double precision :: rref
!--------------------------------------------------------------

rref=0.05d0

if (minval({^D&x(ix^S,^D)**2.0d0+}) <=rref**2.0d0 .and. qt<0.1d0) refine=1

end subroutine specialrefine_grid
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
subroutine printlog_special

use mod_global_parameters
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

use mod_global_parameters

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
!=============================================================================
subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
! corresponding normalization values (default value 1)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)

double precision :: wloc(ixI^S,1:nw)
double precision :: pth(ixG^T), rho(ixG^T)
!-----------------------------------------------------------------------------

wloc(ixI^S,1:nw)=w(ixI^S,1:nw)
if(saveprim)then
  pth(ixO^S)=wloc(ixO^S,p_)
else
  call getpthermal(wloc,x,ixI^L,ixO^L,pth)
endif
rho(ixO^S)=wloc(ixO^S,rho_)
w(ixO^S,nw+1)=pth(ixO^S)*rho(ixO^S)**(-eqpar(gamma_))


end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'S'
wnames=TRIM(wnames)//' '//'S'

end subroutine specialvarnames_output
!=============================================================================
{#IFDEF UCONVERT
subroutine userspecialconvert(qunitconvert)

! Allow user to use their own data-converting procedures

use mod_global_parameters
integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------

end subroutine userspecialconvert
!=============================================================================
}
{#IFDEF TRANSFORMW
subroutine transformw_usr(w,wtf,eqpar_tf,ixI^L,ixO^L)
! regenerate w and eqpar arrays to output into *tf.dat
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: wtf(ixI^S,1:nwtf)
double precision, intent(out):: eqpar_tf(neqpartf)

!-----------------------------------------------------------------------------

end subroutine transformw_usr
!=============================================================================
}
