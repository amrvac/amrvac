!=============================================================================
subroutine errest
 
use mod_forest, only: refine, buffer
include 'amrvacdef.f'

integer :: igrid, iigrid, ixCoGmin1,ixCoGmax1
double precision :: factor
!-----------------------------------------------------------------------------
if (igridstail==0) return

select case (errorestimate)
case (0) 
   ! all refinement solely based on user routine specialrefine_grid
case (1) 
   ! Richardson procedure: compare coarse-integrate versus integrate-coarse
   ! done with low order, dimensionally unsplit scheme typelow1
   ! For error estimate: compare 1st order, coarse 2*dx,2*dt step  with
   ! 1st order dt step, since otherwise wCT can not be filled by interpolation

   ! Note:when there are sources: only unsplit sources are taken into account

   ! Note: the low order solutions are obtained with dimsplit=F.
   !       When overall scheme uses dimsplit=T
   !       and courantpar>0.5, the low order can be unstable. In
   !       that case, simplify this scheme, set low order step of size dt
   !       and compare coarse with available solution at t_n. Enforce this
   !       through `skipfinestep' 


   if (skipfinestep) then
      factor=one
   else
      factor=two
   end if

   ! call the integrator on a grid twice as coarse.

   ixCoGmin1=1;
   ixCoGmax1=ixGhi1/2+dixB;

   call createCoarse(ixCoGmin1,ixCoGmax1)
   call getbc(t,ixCoGmin1,ixCoGmax1,pwCoarse,pwCoCo,pgeoCoarse,pgeoCoCo,&
      .true.)
   call advectCoarse(ixCoGmin1,ixCoGmax1,factor)

   ! Now advance full grid and obtain relative error
   ! between coarse and full grid sln at time t_n+1
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call errest1_grid(igrid,pw(igrid)%w)
   end do

case (2) 
   ! simply compare w_n-1 with w_n and trigger refinement on relative
   ! differences
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call compare1_grid(igrid,pwold(igrid)%w,pw(igrid)%w)
   end do

case (3)
   ! Error estimation is based on Lohner's scheme
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call lohner_grid(igrid)
   end do

case (4)
   ! Error estimation is based on Lohner's original scheme
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call lohner_orig_grid(igrid)
   end do


case default
   call mpistop("Unknown error estimator")
end select

! enforce additional refinement on e.g. coordinate and/or time info here
if (nbufferx1/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,refine,ngridshi&
   *npe,MPI_LOGICAL,MPI_LOR, icomm,ierrmpi)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call forcedrefine_grid(igrid,pw(igrid)%w)
end do

if (nbufferx1/=0) buffer=.false.

end subroutine errest
!=============================================================================
subroutine lohner_grid(igrid)
use mod_forest, only: coarsen, refine
include 'amrvacdef.f'

integer, intent(in) :: igrid

integer :: iiflag, iflag, idims, idims2, level
integer :: ixmin1,ixmax1, hxmin1,hxmax1, jxmin1,jxmax1, h2xmin1,h2xmax1,&
    j2xmin1,j2xmax1, ix1
double precision :: epsilon, tolerance
double precision, dimension(ixMlo1:ixMhi1) :: numerator, denominator, error
double precision, dimension(ixGlo1:ixGhi1) :: tmp, tmp1, tmp2
logical, dimension(ixGlo1:ixGhi1) :: refineflag, coarsenflag
!-----------------------------------------------------------------------------
epsilon=1.0d-6
level=node(plevel_,igrid)
ixmin1=ixMlo1-1;ixmax1=ixMhi1+1;

error=zero
do iiflag=1,flags(nflag_); iflag=flags(iiflag);
   numerator=zero
   if (iflag>nw)call specialvarforerrest(ixGlo1,ixGhi1,ixGlo1,ixGhi1,iflag,&
      pw(igrid)%w,tmp1)
   do idims=1,ndim
      hxmin1=ixmin1-kr(1,idims);hxmax1=ixmax1-kr(1,idims);
      jxmin1=ixmin1+kr(1,idims);jxmax1=ixmax1+kr(1,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          tmp(ixmin1:ixmax1)=dlog10(pw(igrid)%w(jxmin1:jxmax1,iflag))&
             -dlog10(pw(igrid)%w(hxmin1:hxmax1,iflag))
        else
          tmp(ixmin1:ixmax1)=pw(igrid)%w(jxmin1:jxmax1,iflag)&
             -pw(igrid)%w(hxmin1:hxmax1,iflag)
        end if
      else
        if (logflag(iiflag)) then
          tmp(ixmin1:ixmax1)=dlog10(tmp1(jxmin1:jxmax1))-dlog10(tmp1&
             (hxmin1:hxmax1))
        else
          tmp(ixmin1:ixmax1)=tmp1(jxmin1:jxmax1)-tmp1(hxmin1:hxmax1)
        end if
      end if
      do idims2=1,ndim
         h2xmin1=ixMlo1-kr(1,idims2);h2xmax1=ixMhi1-kr(1,idims2);
         j2xmin1=ixMlo1+kr(1,idims2);j2xmax1=ixMhi1+kr(1,idims2);
         numerator=numerator+(tmp(j2xmin1:j2xmax1)-tmp(h2xmin1:h2xmax1))&
            **2.0d0
      end do
   end do
   denominator=zero
   do idims=1,ndim
      if (iflag<=nw) then
         if (logflag(iiflag)) then
          tmp=dabs(dlog10(pw(igrid)%w(ixGlo1:ixGhi1,iflag)))
         else
          tmp=dabs(pw(igrid)%w(ixGlo1:ixGhi1,iflag))
         end if
      else
         if (logflag(iiflag)) then
          tmp=dabs(dlog10(tmp1(ixGlo1:ixGhi1)))
         else
          tmp=dabs(tmp1(ixGlo1:ixGhi1))
         end if
      end if
      hxmin1=ixmin1-kr(1,idims);hxmax1=ixmax1-kr(1,idims);
      jxmin1=ixmin1+kr(1,idims);jxmax1=ixmax1+kr(1,idims);
      tmp2(ixmin1:ixmax1)=tmp(jxmin1:jxmax1)+tmp(hxmin1:hxmax1)
      hxmin1=ixMlo1-2*kr(1,idims);hxmax1=ixMhi1-2*kr(1,idims);
      jxmin1=ixMlo1+2*kr(1,idims);jxmax1=ixMhi1+2*kr(1,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          tmp(ixMlo1:ixMhi1)=dabs(dlog10(pw(igrid)%w(jxmin1:jxmax1,iflag))&
             -dlog10(pw(igrid)%w(ixMlo1:ixMhi1,iflag))) +dabs(dlog10(pw&
             (igrid)%w(ixMlo1:ixMhi1,iflag))-dlog10(pw(igrid)%w(hxmin1:hxmax1,&
             iflag)))
        else
           tmp(ixMlo1:ixMhi1)=dabs(pw(igrid)%w(jxmin1:jxmax1,iflag)&
              -pw(igrid)%w(ixMlo1:ixMhi1,iflag)) +dabs(pw(igrid)%w&
              (ixMlo1:ixMhi1,iflag)-pw(igrid)%w(hxmin1:hxmax1,iflag))
        end if
      else
        if (logflag(iiflag)) then
          tmp(ixMlo1:ixMhi1)=dabs(dlog10(tmp1(jxmin1:jxmax1))&
             -dlog10(tmp1(ixMlo1:ixMhi1))) +dabs(dlog10(tmp1(ixMlo1:ixMhi1))&
             -dlog10(tmp1(hxmin1:hxmax1)))
        else
           tmp(ixMlo1:ixMhi1)=dabs(tmp1(jxmin1:jxmax1)-tmp1(ixMlo1:ixMhi1)) &
              +dabs(tmp1(ixMlo1:ixMhi1)-tmp1(hxmin1:hxmax1))
        end if
      end if
      do idims2=1,ndim
         h2xmin1=ixMlo1-kr(1,idims2);h2xmax1=ixMhi1-kr(1,idims2);
         j2xmin1=ixMlo1+kr(1,idims2);j2xmax1=ixMhi1+kr(1,idims2);
         denominator=denominator +(tmp(ixMlo1:ixMhi1)+amr_wavefilter(level)&
            *(tmp2(j2xmin1:j2xmax1)+tmp2(h2xmin1:h2xmax1)))**2
      end do
   end do
   error=error+wflags(iiflag)*dsqrt(numerator/max(denominator,epsilon))
end do

refineflag=.false.
coarsenflag=.false.
tolerance=tol(level)
do ix1=ixMlo1,ixMhi1

   if (error(ix1) >= tolerance) then
      refineflag(ix1) = .true.
   else if (error(ix1) <= tolratio(level)*tolerance) then
      coarsenflag(ix1) = .true.
   end if
end do


if (any(refineflag(ixMlo1:ixMhi1)).and.level<mxnest) refine(igrid,mype)=.true.
if (all(coarsenflag(ixMlo1:ixMhi1)).and.level>1) coarsen(igrid,mype)=.true.

end subroutine lohner_grid
!=============================================================================
subroutine lohner_orig_grid(igrid)
use mod_forest, only: coarsen, refine
include 'amrvacdef.f'

integer, intent(in) :: igrid

integer :: iiflag, iflag, idims, level
integer :: ixmin1,ixmax1, hxmin1,hxmax1, jxmin1,jxmax1, ix1
double precision :: epsilon
double precision, dimension(ixMlo1:ixMhi1) :: numerator, denominator, error
double precision, dimension(ixGlo1:ixGhi1) :: dp, dm, dref, tmp1
logical, dimension(ixGlo1:ixGhi1) :: refineflag, coarsenflag
!-----------------------------------------------------------------------------
epsilon=1.0d-6
level=node(plevel_,igrid)
ixmin1=ixMlo1;ixmax1=ixMhi1;

error=zero
do iiflag=1,flags(nflag_); iflag=flags(iiflag);
   numerator=zero
   denominator=zero
   if (iflag>nw)call specialvarforerrest(ixGlo1,ixGhi1,ixGlo1,ixGhi1,iflag,&
      pw(igrid)%w,tmp1)
   do idims=1,ndim
      hxmin1=ixmin1-kr(1,idims);hxmax1=ixmax1-kr(1,idims);
      jxmin1=ixmin1+kr(1,idims);jxmax1=ixmax1+kr(1,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          dp(ixmin1:ixmax1)=dlog10(pw(igrid)%w(jxmin1:jxmax1,iflag))&
             -dlog10(pw(igrid)%w(ixmin1:ixmax1,iflag))
          dm(ixmin1:ixmax1)=dlog10(pw(igrid)%w(ixmin1:ixmax1,iflag))&
             -dlog10(pw(igrid)%w(hxmin1:hxmax1,iflag))
          dref(ixMlo1:ixMhi1)=dabs(dlog10(pw(igrid)%w(jxmin1:jxmax1,iflag)))&
             + 2.0d0 * dabs(dlog10(pw(igrid)%w(ixMlo1:ixMhi1,iflag))) &
             + dabs(dlog10(pw(igrid)%w(hxmin1:hxmax1,iflag)))
        else
          dp(ixmin1:ixmax1)=pw(igrid)%w(jxmin1:jxmax1,iflag)&
             -pw(igrid)%w(ixmin1:ixmax1,iflag)
          dp(ixmin1:ixmax1)=pw(igrid)%w(ixmin1:ixmax1,iflag)&
             -pw(igrid)%w(hxmin1:hxmax1,iflag)
          dref(ixMlo1:ixMhi1)=dabs(pw(igrid)%w(jxmin1:jxmax1,iflag))&
             +2.0d0*dabs(pw(igrid)%w(ixMlo1:ixMhi1,iflag)) &
             +dabs(pw(igrid)%w(hxmin1:hxmax1,iflag))
        end if
      else
        if (logflag(iiflag)) then
          dp(ixmin1:ixmax1)=dlog10(tmp1(jxmin1:jxmax1))-dlog10(tmp1&
             (ixmin1:ixmax1))
          dm(ixmin1:ixmax1)=dlog10(tmp1(ixmin1:ixmax1))-dlog10(tmp1&
             (hxmin1:hxmax1))
          dref(ixmin1:ixmax1)=dabs(dlog10(tmp1(jxmin1:jxmax1)))&
             + 2.0d0 * dabs(dlog10(tmp1(ixmin1:ixmax1))) + &
             dabs(dlog10(tmp1(hxmin1:hxmax1)))
        else
          dp(ixmin1:ixmax1)=tmp1(jxmin1:jxmax1)-tmp1(ixmin1:ixmax1)
          dm(ixmin1:ixmax1)=tmp1(ixmin1:ixmax1)-tmp1(hxmin1:hxmax1)
          dref(ixmin1:ixmax1)=dabs(tmp1(jxmin1:jxmax1))+2.0d0&
             *dabs(tmp1(ixmin1:ixmax1)) +dabs(tmp1(hxmin1:hxmax1))
        end if
      end if

      numerator(ixMlo1:ixMhi1)=numerator+(dp(ixMlo1:ixMhi1)&
         -dm(ixMlo1:ixMhi1))**2.0d0

      denominator(ixMlo1:ixMhi1)=denominator + (dabs(dp(ixMlo1:ixMhi1)) &
         + dabs(dm(ixMlo1:ixMhi1)) + amr_wavefilter(level)*dref&
         (ixMlo1:ixMhi1))**2.0d0

   end do
   error=error+wflags(iiflag)*dsqrt(numerator/max(denominator,epsilon))
end do

refineflag=.false.
coarsenflag=.false.

do ix1=ixMlo1,ixMhi1
   if (error(ix1) >= tol(level)) then
      refineflag(ix1) = .true.
   else if (error(ix1) <= tolratio(level)*tol(level)) then
      coarsenflag(ix1) = .true.
   end if
end do

if (any(refineflag(ixMlo1:ixMhi1)).and.level<mxnest) refine(igrid,mype)=.true.
if (all(coarsenflag(ixMlo1:ixMhi1)).and.level>1) coarsen(igrid,mype)=.true.

end subroutine lohner_orig_grid
!=============================================================================
subroutine compare1_grid(igrid,wold,w)
use mod_forest, only: coarsen, refine
include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: wold(ixGlo1:ixGhi1,1:nw), w(ixGlo1:ixGhi1,&
   1:nw)

integer :: ix1, iiflag, iflag, level
double precision :: epsilon
double precision :: average, error
double precision :: averages(nflag_)
logical, dimension(ixGlo1:ixGhi1) :: refineflag, coarsenflag
!-----------------------------------------------------------------------------
! identify the points to be flagged in two steps:
!  step I: compare w_n-1 with w_n solution, store flags in auxiliary
!  step II: transfer flags from auxiliary to refine and coarsen

epsilon=1.0d-6

refineflag(ixMlo1:ixMhi1) = .false.
coarsenflag(ixMlo1:ixMhi1) = .false.
level=node(plevel_,igrid)

do ix1=ixMlo1,ixMhi1 
   average=zero
   error=zero
   do iiflag=1,flags(nflag_); iflag=flags(iiflag);
      averages(iflag) = w(ix1,iflag)-wold(ix1,iflag)
      average=average+wflags(iiflag)*abs(averages(iflag))
      if (abs(wold(ix1,iflag))<smalldouble)then
         error=error+wflags(iiflag)* abs(averages(iflag))/(abs(wold(ix1,&
            iflag))+epsilon)
      else
         error=error+wflags(iiflag)* abs(averages(iflag))/(abs(wold(ix1,&
            iflag)))
      end if
   end do
   if (error >= tol(level)) then
      refineflag(ix1) = .true.
   else if (error <= tolratio(level)*tol(level)) then
      coarsenflag(ix1) = .true.
   end if
end do

if (any(refineflag(ixMlo1:ixMhi1))) then
   if (level<mxnest) refine(igrid,mype)=.true.
end if
if (time_advance) then
   if (all(coarsenflag(ixMlo1:ixMhi1)).and.level>1) coarsen(igrid,mype)=.true.
end if

end subroutine compare1_grid
!=============================================================================
subroutine createCoarse(ixCoGmin1,ixCoGmax1)
include 'amrvacdef.f'

integer, intent(in) :: ixCoGmin1,ixCoGmax1

integer :: iigrid, igrid
integer :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3
!-----------------------------------------------------------------------------
ixGmin1=ixCoGmin1;ixGmax1=ixCoGmax1;
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call createCoarse_grid(igrid,pwCoarse(igrid),pxCoarse(igrid),ixGmin1,&
      ixGmax1,pwold(igrid)%w,px(igrid)%x)
end do

end subroutine createCoarse
!=============================================================================
subroutine createCoarse_grid(igrid,pwCo,pxCo,ixCoGmin1,ixCoGmax1,wold,xold)

include 'amrvacdef.f'

integer, intent(in) :: igrid, ixCoGmin1,ixCoGmax1
double precision :: wold(ixGlo1:ixGhi1,1:nw), xold(ixGlo1:ixGhi1,1:ndim)
type(walloc) pwCo
type(xalloc) pxCo

integer :: ixCoMmin1,ixCoMmax1
!-----------------------------------------------------------------------------
dxlevel(1)=rnode(rpdx1_,igrid);

! now coarsen by 2 in every direction - conservatively
! coarse grid dimension half its size
ixCoMmin1=ixCoGmin1+dixB;ixCoMmax1=ixCoGmax1-dixB;
call coarsen_grid(wold,xold,ixGlo1,ixGhi1,ixMlo1,ixMhi1,pwCo%w,pxCo%x,&
   ixCoGmin1,ixCoGmax1,ixCoMmin1,ixCoMmax1, pgeo(igrid),pgeoCoarse(igrid),&
   coarsenprimitive,.true.)

end subroutine createCoarse_grid
!=============================================================================
subroutine advectCoarse(ixCoGmin1,ixCoGmax1,factor)
include 'amrvacdef.f'

integer :: ixCoGmin1,ixCoGmax1
double precision, intent(in) :: factor

integer :: iigrid, igrid
!-----------------------------------------------------------------------------
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call advectCoarse_grid(igrid,pwCoarse(igrid),ixCoGmin1,ixCoGmax1,factor)
end do

end subroutine advectCoarse
!=============================================================================
subroutine advectCoarse_grid(igrid,pwCo,ixCoGmin1,ixCoGmax1,factor)

include 'amrvacdef.f'

integer, intent(in) :: igrid, ixCoGmin1,ixCoGmax1
double precision, intent(in) :: factor
type(walloc) pwCo

double precision :: qdt, dx1
double precision, dimension(:,:), allocatable :: wCo1
double precision :: fC(ixCoGmin1:ixCoGmax1,1:nwflux,1:ndim)
integer :: level
!-----------------------------------------------------------------------------
dxlevel(1)=two*rnode(rpdx1_,igrid);
dx1=dxlevel(1);

! here we integrate on the coarse grid
allocate(wCo1(ixCoGmin1:ixCoGmax1,1:nw))
wCo1(ixCoGmin1:ixCoGmax1,1:nwflux)=pwCo%w(ixCoGmin1:ixCoGmax1,1:nwflux)

! 1st order scheme: do coarse time step of
!    size 2*dt starting from t_n-1 solution in pwCoarse
!    to arrive at t_n+1 (n-index from normal uncoarsened grid)
! result in pwCoarse: coarse solution at t_n+1

if (.not.slab) mygeo => pgeoCoarse(igrid)

qdt=factor*dt_grid(igrid)
level=node(plevel_,igrid)
call advect1_grid(typelow1(level),qdt,ixCoGmin1,ixCoGmax1,1,ndim,t,wCo1,t,&
    pwCo%w,wCo1,fC,dx1,pxCoarse(igrid)%x)

deallocate(wCo1)

end subroutine advectCoarse_grid
!=============================================================================
subroutine errest1_grid(igrid,w)

include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: w(ixGlo1:ixGhi1,nw)

integer :: level, ixCoGmin1,ixCoGmax1
double precision :: dx1, qdt
double precision :: fC(ixGlo1:ixGhi1,1:nwflux,1:ndim), wFi(ixGlo1:ixGhi1,1:nw)
!-----------------------------------------------------------------------------
level=node(plevel_,igrid)

dx1=dx(1,level);
dxlevel(1)=dx1;

wFi(ixGlo1:ixGhi1,1:nwflux)=w(ixGlo1:ixGhi1,1:nwflux)

if (.not.skipfinestep) then
   if (.not.slab) mygeo => pgeo(igrid)

   qdt=dt_grid(igrid)
   call advect1_grid(typelow1(level),qdt,ixGlo1,ixGhi1,1,ndim,t+qdt,w, t&
      +qdt,wFi,w,fC,dx1,px(igrid)%x)
end if

ixCoGmin1=1;
ixCoGmax1=ixGhi1/2+dixB;

call flagbadpoints(wFi,pwCoarse(igrid)%w,ixCoGmin1,ixCoGmax1,igrid,level)

end subroutine errest1_grid
!=============================================================================
subroutine flagbadpoints(w,wCo,ixCoGmin1,ixCoGmax1,igrid,level)

! compare error between coarse and fine solution in wCo, w 
! We base the comparison on the physical field selected by the index flag_ 
!
! on entry:
!  w:       normal time integration 
!  wCo:     time integration on coarsened grid (2*dx)

use mod_forest, only: coarsen, refine
include 'amrvacdef.f'

integer, intent(in)         :: igrid, ixCoGmin1,ixCoGmax1, level
double precision,intent(in) :: w(ixGlo1:ixGhi1,nw), wCo(ixCoGmin1:ixCoGmax1,&
   nw)
double precision            :: specialvar(ixGlo1:ixGhi1), specialvarCo&
   (ixCoGmin1:ixCoGmax1)

logical :: needgetaux
integer :: iCo1, iFi1, ixCoMmin1,ixCoMmax1, iiflag, iflag
double precision :: average, error
double precision :: averages(nflag_)
logical, dimension(ixGlo1:ixGhi1) :: refineflag, coarsenflag
!-----------------------------------------------------------------------------
ixCoMmin1=ixCoGmin1+dixB;ixCoMmax1=ixCoGmax1-dixB;

needgetaux=.false.
do iiflag=1,flags(nflag_); iflag=flags(iiflag);
  if (iflag>nwflux) needgetaux=.true.
end do
if (nwaux>0.and.needgetaux) then
   saveigrid=igrid
   if(.not.slab)mygeo=>pgeo(igrid)
   call getaux(.true.,w,px(igrid)%x,ixGlo1,ixGhi1,ixMlo1,ixMhi1,&
      'flagbadpoints')
   call getaux(.true.,wCo,pxCoarse(igrid)%x,ixCoGmin1,ixCoGmax1,ixCoMmin1,&
      ixCoMmax1,'flagbadpointsCo')
end if

! identify the points to be flagged in two steps (needed!):
!  step I: compare coarse with fine solution, store flags in fine auxiliary
!  step II: transfer flags from auxiliary to refine and coarsen

refineflag(ixMlo1:ixMhi1) = .false.
coarsenflag(ixMlo1:ixMhi1) = .false.


do iiflag=1,flags(nflag_); iflag=flags(iiflag);
   if (iflag>nw) then
      call specialvarforerrest(ixCoGmin1,ixCoGmax1,ixCoGmin1,ixCoGmax1,iflag,&
         wCo,specialvarCo)
      call specialvarforerrest(ixGlo1,ixGhi1,ixGlo1,ixGhi1,iflag,w,specialvar)
   end if
end do

iFi1 = ixMlo1
do iCo1 = ixCoMmin1,ixCoMmax1 
   average=zero
   error=zero
   do iiflag=1,flags(nflag_); iflag=flags(iiflag);
      if (slab) then
         if (iflag<=nw) averages(iflag)=sum(w(iFi1:iFi1+1,iflag))/two**ndim
         if (iflag>nw)  averages(iflag)=sum(specialvar(iFi1:iFi1+1))/two**ndim
      else
         if (iflag<=nw) averages(iflag)=sum(pgeo(igrid)%dvolume(iFi1:iFi1&
            +1) *w(iFi1:iFi1+1,iflag))/pgeoCoarse(igrid)%dvolume(iCo1)
         if (iflag>nw)  averages(iflag)=sum(pgeo(igrid)%dvolume(iFi1:iFi1&
            +1) *specialvar(iFi1:iFi1+1))/pgeoCoarse(igrid)%dvolume(iCo1)
      end if
      average=average+wflags(iiflag)*abs(averages(iflag))
      if (iflag<=nw) error=error+wflags(iiflag)*abs(averages(iflag)&
         -wCo(iCo1,iflag))
      if (iflag> nw) error=error+wflags(iiflag)*abs(averages(iflag)&
         -specialvarCo(iCo1))
   end do
   if (abs(average)>smalldouble) then
      error=error/average
   else
      write(unitterm,*)'Warning from flagbadpoints: zero average:',average
      write(unitterm,*)'   wCo(iCo^D,1:nw):',wCo(iCo1,1:nw),' indices:',iCo1
      write(unitterm,*)'On grid:',igrid,' at level ',level
      write(unitterm,*)'   and grid indices : ',node(pig1_,igrid)
      write(unitterm,*)'cell indices : ',iCo1
      call mpistop("")
   end if
   if (error >= tol(level)) then
      refineflag(iFi1:iFi1+1) = .true.
   else if (error <= tolratio(level)*tol(level)) then
      coarsenflag(iFi1:iFi1+1) = .true.
   end if
   iFi1 = iFi1+2
end do

iFi1 = ixMlo1
do iCo1 = ixCoMmin1,ixCoMmax1 
   if (error >= tol(level)) then
      refineflag(iFi1:iFi1+1) = .true.
   else if (error <= tolratio(level)*tol(level)) then
      coarsenflag(iFi1:iFi1+1) = .true.
   end if
   iFi1 = iFi1+2
end do

if (any(refineflag(ixMlo1:ixMhi1))) then
   if (level<mxnest) refine(igrid,mype)=.true.
end if
if (time_advance) then
   if (all(coarsenflag(ixMlo1:ixMhi1)).and.level>1) coarsen(igrid,mype)=.true.
end if

end subroutine flagbadpoints
!=============================================================================
subroutine forcedrefine_grid(igrid,w)
use mod_forest, only: coarsen, refine, buffer
include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: w(ixGlo1:ixGhi1,nw)

integer :: level
integer :: my_refine, my_coarsen
double precision :: qt
logical, dimension(ixGlo1:ixGhi1) :: refineflag
!-----------------------------------------------------------------------------
level=node(plevel_,igrid)

! initialize to 0
my_refine   = 0
my_coarsen  = 0

if (time_advance) then
   qt=t+dt
else
   if (errorestimate==1.or.errorestimate==2) then
      qt=t+dt
   else
      qt=t
   end if
end if
   
call specialrefine_grid(igrid,level,ixGlo1,ixGhi1,ixMlo1,ixMhi1,qt,w,&
   px(igrid)%x, my_refine,my_coarsen)

if (my_coarsen==1) then
   if (level>1) then
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.true.
   else
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.false.
   end if
endif

if (my_coarsen==-1)then
   coarsen(igrid,mype)=.false.
end if

if (my_refine==1) then
   if (level<mxnest) then
      refine(igrid,mype)=.true.
      coarsen(igrid,mype)=.false.
   else
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.false.
   end if
end if

if (my_refine==-1) then
  refine(igrid,mype)=.false.
end if

if (nbufferx1/=0) then
   if (refine(igrid,mype) .and. .not.buffer(igrid,mype)) then
      refineflag(ixMlo1:ixMhi1)=.true.
      call refinebuffer(igrid,refineflag)
   end if
end if

end subroutine forcedrefine_grid
!=============================================================================
subroutine forcedrefine_grid_io(igrid,w)
use mod_forest, only: coarsen, refine
include 'amrvacdef.f'

integer, intent(in)          :: igrid
double precision, intent(in) :: w(ixGlo1:ixGhi1,nw)

integer                   :: level, my_levmin, my_levmax
logical, dimension(ixGlo1:ixGhi1) :: refineflag
!-----------------------------------------------------------------------------
level=node(plevel_,igrid)

if (level_io > 0) then
   my_levmin = level_io
   my_levmax = level_io
else
   my_levmin = max(1,level_io_min)
   my_levmax = min(mxnest,level_io_max)
end if


if (level>my_levmax) then
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.true.
elseif (level<my_levmin) then
      refine(igrid,mype)=.true.
      coarsen(igrid,mype)=.false.
end if

if (level==my_levmin .or. level==my_levmax) then
  refine(igrid,mype)=.false.
  coarsen(igrid,mype)=.false.
end if


if(refine(igrid,mype).and.level>=mxnest)refine(igrid,mype)=.false.
if(coarsen(igrid,mype).and.level<=1)coarsen(igrid,mype)=.false.

end subroutine forcedrefine_grid_io
!=============================================================================
subroutine refinebuffer(igrid,refineflag)
use mod_forest, only: refine, buffer
include 'amrvacdef.f'

integer, intent(in) :: igrid
logical, dimension(ixGlo1:ixGhi1), intent(in) :: refineflag

integer :: ishiftbuf1, i1, ixmin1,ixmax1, ineighbor, ipe_neighbor, level
!-----------------------------------------------------------------------------
ishiftbuf1=ixMhi1-ixMlo1-nbufferx1+1;
do i1=-1,1
   ixmin1=max(ixMlo1,ixMlo1+i1*ishiftbuf1);
   ixmax1=min(ixMhi1,ixMhi1+i1*ishiftbuf1);
   if (ixmax1<ixmin1) cycle
   if (any(refineflag(ixmin1:ixmax1))) then
      select case (neighbor_type(i1,igrid))
      case (2)
         ineighbor=neighbor(1,i1,igrid)
         ipe_neighbor=neighbor(2,i1,igrid)
         if (.not.refine(ineighbor,ipe_neighbor)) then
            buffer(ineighbor,ipe_neighbor)=.true.
            refine(ineighbor,ipe_neighbor)=.true.
         end if
      case (3)
         level=node(plevel_,igrid)
         if (level<mxnest) then
            ineighbor=neighbor(1,i1,igrid)
            ipe_neighbor=neighbor(2,i1,igrid)
            if (.not.refine(ineighbor,ipe_neighbor)) then
               buffer(ineighbor,ipe_neighbor)=.true.
               refine(ineighbor,ipe_neighbor)=.true.
            end if
         end if
      end select
   end if
end do

end subroutine refinebuffer
!=============================================================================
