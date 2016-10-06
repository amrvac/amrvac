!=============================================================================
! amrvacusr.t.mhdconvection
!=============================================================================
!INCLUDE:amrvacmodules/heatconduct.t
!INCLUDE:amrvacmodules/gravity.t
!INCLUDE:amrvacmodules/viscosity.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters

!-----------------------------------------------------------------------------

! equilibrium parameters
eqpar(gamma_)=5.0d0/3.0d0

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision       :: A, phi(ixG^S), ca, alpha, beta
double precision       :: v(ixG^S,1:3), b(ixG^S,1:3), vrot(ixG^S,1:3), brot(ixG^S,1:3), k(1:^ND)
logical :: patchw(ixG^T)
logical, save :: first=.true.
!----------------------------------------------------------------------------

A     = 0.1d0
ca    = 1.0d0
k(1)  = 2.0d0 * dpi

{^IFONED
alpha = 0.0d0
beta  = 0.0d0
phi(ix^S) = k(1) * x(ix^S,1)
}

{^IFTWOD 
alpha = atan(2.0d0)
beta  = 0.0d0
k(2)  = k(1) * tan(alpha)
phi(ix^S) = k(1) *x(ix^S,1) + k(2) * x(ix^S,2)
}

{^IFTHREED
alpha = atan(2.0d0)
beta  = atan(2.0d0)
k(2)  = k(1) * tan(alpha)
k(3)  = k(1) * tan(beta)
phi(ix^S) = k(1) *x(ix^S,1) + k(2) * x(ix^S,2) + k(3) * x(ix^S,3)
}

w(ix^S,rho_) = 1.0d0
w(ix^S,p_)   = 0.1d0

v(ix^S,1)  = 0.0d0
v(ix^S,2)  = A * sin(phi(ix^S))
v(ix^S,3)  = A * cos(phi(ix^S))

b(ix^S,1)  =   sqrt(w(ix^S,rho_)) * ca
b(ix^S,2)  = - sqrt(w(ix^S,rho_)) * ca * A * sin(phi(ix^S))
b(ix^S,3)  = - sqrt(w(ix^S,rho_)) * ca * A * cos(phi(ix^S))

call rotate(ixG^L,ix^L,v,vrot,alpha,beta)
call rotate(ixG^L,ix^L,b,brot,alpha,beta)

{^C& w(ix^S,v^C_) = vrot(ix^S,^C) \}
{^C& w(ix^S,b^C_) = brot(ix^S,^C) \}

patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

end subroutine initonegrid_usr
!=============================================================================
subroutine rotate(ixG^L,ix^L,v,vrot,alpha,beta)
use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: v(ixG^S,1:3)
double precision, intent(out) :: vrot(ixG^S,1:3)
double precision, intent(in) :: alpha, beta

double precision       :: gamma
!-----------------------------------------------------------------------------

gamma = atan(cos(alpha)*tan(beta))

vrot(ix^S,1) = cos(gamma) * cos(alpha) * v(ix^S,1) - sin(alpha) * v(ix^S,2) - sin(gamma) * cos(alpha) * v(ix^S,3)
vrot(ix^S,2) = cos(gamma) * sin(alpha) * v(ix^S,1) + cos(alpha) * v(ix^S,2) - sin(gamma) * sin(alpha) * v(ix^S,3)
vrot(ix^S,3) = sin(gamma)              * v(ix^S,1)                          + cos(gamma) * v(ix^S,3)

end subroutine rotate
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixO^L, iw, iB, ixG^L
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: tempb(ixG^T)
double precision :: delydelx,delydelz
integer :: ix^D,ixIM^L

logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("Special boundary is not defined")


end subroutine specialbound_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)
!----------------------------------------------------------------------------


end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
! note that depending on strictsmall etc, w values may change
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------

end subroutine getdt_special
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

end subroutine getdt_impl
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the common "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_global_parameters

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

use mod_global_parameters

integer, intent(in) :: igrid, level, ix^L, ixG^L
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

use mod_global_parameters

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')

var(ixI^S) = zero

end subroutine specialvarforerrest
!=============================================================================
subroutine printlog_special

use mod_global_parameters
!-----------------------------------------------------------------------------

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

!double precision :: tmp(ixG^T),tmp1(ixG^T),tmp2(ixG^T)
!double precision :: qvec(ixG^T,1:ndir),curlvec(ixG^T,1:ndir)
!integer :: idirmin
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

!-----------------------------------------------------------------------------

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the wnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

end subroutine specialvarnames_output
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------

end subroutine specialset_B0
!=============================================================================
subroutine bc_int(level,qt,ixG^L,ixO^L,w,x)

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

integer, intent(in) :: ixG^L,ixO^L,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

end subroutine bc_int
!=============================================================================
subroutine userspecialconvert(qunitconvert)

use mod_global_parameters
!double precision, intent(inout) :: w(ixI^S,nw+nwauxio)
!integer, intent(in)           :: ixI^L,ixO^L

integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------

end subroutine userspecialconvert
!=============================================================================
! amrvacusr.t.mhdconvection
!=============================================================================
