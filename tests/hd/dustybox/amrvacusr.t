!=============================================================================
! amrvacusr.t.KHDustM

!INCLUDE:amrvacnul/specialini.t
INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
!INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters

integer :: i
DOUBLE PRECISION :: r(0:^NDS)
logical, save :: first=.true.
!-----------------------------------------------------------------------------

eqpar(gamma_) = 5.0d0/3.0d0
eqpar(mu_)    = 1.00d0      !average particle mass

rhodust(1:^NDS)  = 3.3d0     ! density in the dust
eqpar(min_ar_)=5.0d-8
eqpar(max_ar_)=250.0d-8


eqpar(rho1_) = 1.0d-20
eqpar(vel1_) = 5.0d4
eqpar(T1_)   = 1.0d2

normvar(0)     = 1.0d18          ! normalization for distance
normvar(rho_)  = 1.0d-21         ! normalization for rho
normvar(v1_)   = 1.0d7            ! normalization for speed

normt          = normvar(0)/normvar(v1_)
normvar(p_)    = normvar(rho_)*(normvar(v1_)**2)         
{normvar(rhod^DS_)   = normvar(rho_)\}
{^DS&{^C&normvar(v^Cd^DS_) = normvar(v^C_);}\}

rhodust(1:^NDS) = rhodust(1:^NDS)/normvar(rhod1_)
eqpar(min_ar_)  = eqpar(min_ar_)/normvar(0)
eqpar(max_ar_)  = eqpar(max_ar_)/normvar(0)


! if not using "dustmethod='linear'", define rhodust(1:^NDS), 
! sdust(1:^NDS) and eqpar(mu_)
!-------------------------------

! here the dust sizes are defined. Note that this is 
! now done differently for the method used in (van Marle et al. (2011)). 

! first dust sizes in Ndust bins, with all bins having equal total mass.
! To do this, assume the particle distribution goes as r^-3.5

r(0) = eqpar(min_ar_)
do i=1,^NDS
    r(i) = (dsqrt(r(i-1)) +(dsqrt(eqpar(max_ar_))- &
        dsqrt(eqpar(min_ar_)))/^NDS)**2.0d0
    dsdust(i) = r(i)-r(i-1)
end do    
! now calculate the weigthed mean size of each bin, again assuming n goes as r^-3.5
do i=1,^NDS
    sdust(i) = (5.0d0/3.0d0)*(r(i)**(-1.5d0) - r(i-1)**(-1.5d0)) &
        /(r(i)**(-2.5d0) - r(i-1)**(-2.5d0))
end do

if(first)then
  if(mype==0)then
    do i=1,^NDS
        write(*,*) 'Dust type ',i,': grain radius r=',sdust(i)
    end do
  endif
  first=.false.
endif
!-------------------------------------


do i=1,nw
   if(loglimit(i) .and. (.not. (i==rho_ .or. i==p_))) then
   call mpistop( ' Taking the logarithm of a negative number is a REALLY STUPID IDEA ')
   endif
enddo

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: delr
logical :: patchw(ixG^T)
double precision :: vsound
!----------------------------------------------------------------------------

{^IFTWOD call mpistop("This is a 1D HDDust Riemann problem!!!")}
{^IFTHREED call mpistop("This is a 1D HDDust Riemann problem!!!")}
{^IFONED



w(ix^S,v1_)=0.0d0
w(ix^S,rho_)=eqpar(rho1_)/normvar(rho_)
w(ixG^S,p_)   = eqpar(rho1_)*(kbcgspar/(mhcgspar*eqpar(mu_)))*eqpar(T1_) / normvar(p_)
{^DS&w(ixG^S,rhod^DS_)=0.01d0*eqpar(rho1_)/(normvar(rho_)*^NDS);}
{^DS&w(ixG^S,v1d^DS_) = eqpar(vel1_)/normvar(v1_);}


}

patchw(ix^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)



end subroutine initonegrid_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

integer :: iw
!-----------------------------------------------------------------------------

! call addsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

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

integer :: idir
integer                                             :: idims,idust
double precision, dimension(1:^NDS)                 :: dtdust
double precision, dimension(ixG^T,1:^NC,1:^NDS)     :: fdrag
double precision, dimension(ixG^T) :: vt2,deltav,tstop,ptherm,vdust,vgas
!-----------------------------------------------------------------------------



end subroutine getdt_special
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

double precision :: rad(ixG^T), mid_y
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
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)

end subroutine specialset_B0
!=============================================================================
! amrvacusr.t.RiemannHDDust
!=============================================================================
