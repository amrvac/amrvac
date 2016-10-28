!=============================================================================
subroutine errest
 
use mod_forest, only: refine, buffer
use mod_global_parameters

integer :: igrid, iigrid, ixCoG^L
double precision :: factor
!-----------------------------------------------------------------------------
if (igridstail==0) return

select case (errorestimate)
case (0) 
   ! all refinement solely based on user routine specialrefine_grid
case (1) 
   ! simply compare w_n-1 with w_n and trigger refinement on relative
   ! differences
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call compare1_grid(igrid,pwold(igrid)%w,pw(igrid)%w)
   end do
!$OMP END PARALLEL DO
case (2)
   ! Error estimation is based on Lohner's original scheme
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call lohner_orig_grid(igrid)
   end do
!$OMP END PARALLEL DO

case (3)
   ! Error estimation is based on Lohner's scheme
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call lohner_grid(igrid)
   end do
!$OMP END PARALLEL DO



case default
   call mpistop("Unknown error estimator")
end select

! enforce additional refinement on e.g. coordinate and/or time info here
if (nbufferx^D/=0|.or.) &
   call MPI_ALLREDUCE(MPI_IN_PLACE,refine,ngridshi*npe,MPI_LOGICAL,MPI_LOR, &
                      icomm,ierrmpi)
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call forcedrefine_grid(igrid,pw(igrid)%w)
end do
!$OMP END PARALLEL DO

if (nbufferx^D/=0|.or.) &
buffer=.false.

end subroutine errest
!=============================================================================
subroutine lohner_grid(igrid)
use mod_forest, only: coarsen, refine
use mod_global_parameters

integer, intent(in) :: igrid

integer :: iiflag, iflag, idims, idims2, level
integer :: ix^L, hx^L, jx^L, h2x^L, j2x^L, ix^D
double precision :: epsilon, tolerance, wtol(1:nw),xtol(1:ndim)
double precision, dimension(ixM^T) :: numerator, denominator, error
double precision, dimension(ixG^T) :: tmp, tmp1, tmp2
logical, dimension(ixG^T) :: refineflag, coarsenflag
!-----------------------------------------------------------------------------
epsilon=1.0d-6
level=node(plevel_,igrid)
ix^L=ixM^LL^LADD1;

error=zero
do iiflag=1,flags(nflag_); iflag=flags(iiflag);
   numerator=zero
   if (iflag>nw)call specialvarforerrest(ixG^LL,ixG^LL,iflag,pw(igrid)%w,tmp1)
   do idims=1,ndim
      hx^L=ix^L-kr(^D,idims);
      jx^L=ix^L+kr(^D,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          tmp(ix^S)=dlog10(pw(igrid)%w(jx^S,iflag))-dlog10(pw(igrid)%w(hx^S,iflag))
        else
          tmp(ix^S)=pw(igrid)%w(jx^S,iflag)-pw(igrid)%w(hx^S,iflag)
        end if
      else
        if (logflag(iiflag)) then
          tmp(ix^S)=dlog10(tmp1(jx^S))-dlog10(tmp1(hx^S))
        else
          tmp(ix^S)=tmp1(jx^S)-tmp1(hx^S)
        end if
      end if
      do idims2=1,ndim
         h2x^L=ixM^LL-kr(^D,idims2);
         j2x^L=ixM^LL+kr(^D,idims2);
         numerator=numerator+(tmp(j2x^S)-tmp(h2x^S))**2.0d0
      end do
   end do
   denominator=zero
   do idims=1,ndim
      if (iflag<=nw) then
         if (logflag(iiflag)) then
          tmp=dabs(dlog10(pw(igrid)%w(ixG^T,iflag)))
         else
          tmp=dabs(pw(igrid)%w(ixG^T,iflag))
         end if
      else
         if (logflag(iiflag)) then
          tmp=dabs(dlog10(tmp1(ixG^T)))
         else
          tmp=dabs(tmp1(ixG^T))
         end if
      end if
      hx^L=ix^L-kr(^D,idims);
      jx^L=ix^L+kr(^D,idims);
      tmp2(ix^S)=tmp(jx^S)+tmp(hx^S)
      hx^L=ixM^LL-2*kr(^D,idims);
      jx^L=ixM^LL+2*kr(^D,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          tmp(ixM^T)=dabs(dlog10(pw(igrid)%w(jx^S,iflag))&
                     -dlog10(pw(igrid)%w(ixM^T,iflag))) &
                     +dabs(dlog10(pw(igrid)%w(ixM^T,iflag))&
                     -dlog10(pw(igrid)%w(hx^S,iflag)))
        else
           tmp(ixM^T)=dabs(pw(igrid)%w(jx^S,iflag)-pw(igrid)%w(ixM^T,iflag)) &
                      +dabs(pw(igrid)%w(ixM^T,iflag)-pw(igrid)%w(hx^S,iflag))
        end if
      else
        if (logflag(iiflag)) then
          tmp(ixM^T)=dabs(dlog10(tmp1(jx^S))-dlog10(tmp1(ixM^T))) &
                    +dabs(dlog10(tmp1(ixM^T))-dlog10(tmp1(hx^S)))
        else
           tmp(ixM^T)=dabs(tmp1(jx^S)-tmp1(ixM^T)) &
                      +dabs(tmp1(ixM^T)-tmp1(hx^S))
        end if
      end if
      do idims2=1,ndim
         h2x^L=ixM^LL-kr(^D,idims2);
         j2x^L=ixM^LL+kr(^D,idims2);
         denominator=denominator &
                    +(tmp(ixM^T)+amr_wavefilter(level)*(tmp2(j2x^S)+tmp2(h2x^S)))**2
      end do
   end do
   error=error+wflags(iiflag)*dsqrt(numerator/max(denominator,epsilon))
end do

refineflag=.false.
coarsenflag=.false.
tolerance=tol(level)
{do ix^DB=ixMlo^DB,ixMhi^DB\}
   if (specialtol) then
      wtol(1:nw)=pw(igrid)%w(ix^D,1:nw)
      xtol(1:ndim)=px(igrid)%x(ix^D,1:ndim)
      call special_tolerance(wtol,xtol,tolerance,t)
   end if
   if (error(ix^D) >= tolerance) then
      refineflag(ix^D) = .true.
   else if (error(ix^D) <= tolratio(level)*tolerance) then
      coarsenflag(ix^D) = .true.
   end if
{end do\}


if (any(refineflag(ixM^T)).and.level<mxnest) refine(igrid,mype)=.true.
if (all(coarsenflag(ixM^T)).and.level>1) coarsen(igrid,mype)=.true.

end subroutine lohner_grid
!=============================================================================
subroutine lohner_orig_grid(igrid)
use mod_forest, only: coarsen, refine
use mod_global_parameters

integer, intent(in) :: igrid

integer :: iiflag, iflag, idims, level
integer :: ix^L, hx^L, jx^L, ix^D
double precision :: epsilon, tolerance, wtol(1:nw), xtol(1:ndim)
double precision, dimension(ixM^T) :: numerator, denominator, error
double precision, dimension(ixG^T) :: dp, dm, dref, tmp1
logical, dimension(ixG^T) :: refineflag, coarsenflag
!-----------------------------------------------------------------------------
epsilon=1.0d-6
level=node(plevel_,igrid)
ix^L=ixM^LL;

error=zero
do iiflag=1,flags(nflag_); iflag=flags(iiflag);
   numerator=zero
   denominator=zero
   if (iflag>nw)call specialvarforerrest(ixG^LL,ixG^LL,iflag,pw(igrid)%w,tmp1)
   do idims=1,ndim
      hx^L=ix^L-kr(^D,idims);
      jx^L=ix^L+kr(^D,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          dp(ix^S)=dlog10(pw(igrid)%w(jx^S,iflag))-dlog10(pw(igrid)%w(ix^S,iflag))
          dm(ix^S)=dlog10(pw(igrid)%w(ix^S,iflag))-dlog10(pw(igrid)%w(hx^S,iflag))
          dref(ixM^T)=dabs(dlog10(pw(igrid)%w(jx^S,iflag)))&
                     + 2.0d0 * dabs(dlog10(pw(igrid)%w(ixM^T,iflag))) &
                     + dabs(dlog10(pw(igrid)%w(hx^S,iflag)))
        else
          dp(ix^S)=pw(igrid)%w(jx^S,iflag)-pw(igrid)%w(ix^S,iflag)
          dp(ix^S)=pw(igrid)%w(ix^S,iflag)-pw(igrid)%w(hx^S,iflag)
          dref(ixM^T)=dabs(pw(igrid)%w(jx^S,iflag))+2.0d0*dabs(pw(igrid)%w(ixM^T,iflag)) &
                      +dabs(pw(igrid)%w(hx^S,iflag))
        end if
      else
        if (logflag(iiflag)) then
          dp(ix^S)=dlog10(tmp1(jx^S))-dlog10(tmp1(ix^S))
          dm(ix^S)=dlog10(tmp1(ix^S))-dlog10(tmp1(hx^S))
          dref(ix^S)=dabs(dlog10(tmp1(jx^S)))&
                     + 2.0d0 * dabs(dlog10(tmp1(ix^S))) &
                     + dabs(dlog10(tmp1(hx^S)))
        else
          dp(ix^S)=tmp1(jx^S)-tmp1(ix^S)
          dm(ix^S)=tmp1(ix^S)-tmp1(hx^S)
          dref(ix^S)=dabs(tmp1(jx^S))+2.0d0*dabs(tmp1(ix^S)) &
                      +dabs(tmp1(hx^S))
        end if
      end if

      numerator(ixM^T)=numerator+(dp(ixM^T)-dm(ixM^T))**2.0d0

      denominator(ixM^T)=denominator &
           + (dabs(dp(ixM^T)) + dabs(dm(ixM^T)) + amr_wavefilter(level)*dref(ixM^T))**2.0d0

   end do
   error=error+wflags(iiflag)*dsqrt(numerator/max(denominator,epsilon))
end do

refineflag=.false.
coarsenflag=.false.

tolerance=tol(level)
{do ix^DB=ixMlo^DB,ixMhi^DB\}
   if (specialtol) then
      wtol(1:nw)=pw(igrid)%w(ix^D,1:nw)
      xtol(1:ndim)=px(igrid)%x(ix^D,1:ndim)
      call special_tolerance(wtol,xtol,tolerance,t)
   end if
   if (error(ix^D) >= tolerance) then
      refineflag(ix^D) = .true.
   else if (error(ix^D) <= tolratio(level)*tolerance) then
      coarsenflag(ix^D) = .true.
   end if
{end do\}

if (any(refineflag(ixM^T)).and.level<mxnest) refine(igrid,mype)=.true.
if (all(coarsenflag(ixM^T)).and.level>1) coarsen(igrid,mype)=.true.

end subroutine lohner_orig_grid
!=============================================================================
subroutine compare1_grid(igrid,wold,w)
use mod_forest, only: coarsen, refine
use mod_global_parameters

integer, intent(in) :: igrid
double precision, intent(in) :: wold(ixG^T,1:nw), w(ixG^T,1:nw)

integer :: ix^D, iiflag, iflag, level
double precision :: epsilon, tolerance, wtol(1:nw), xtol(1:ndim)
double precision :: average, error
double precision :: averages(nflag_)
logical, dimension(ixG^T) :: refineflag, coarsenflag
!-----------------------------------------------------------------------------
! identify the points to be flagged in two steps:
!  step I: compare w_n-1 with w_n solution, store flags in auxiliary
!  step II: transfer flags from auxiliary to refine and coarsen

epsilon=1.0d-6

refineflag(ixM^T) = .false.
coarsenflag(ixM^T) = .false.
level=node(plevel_,igrid)
tolerance=tol(level)
{do ix^DB=ixMlo^DB,ixMhi^DB \}
   average=zero
   error=zero
   do iiflag=1,flags(nflag_); iflag=flags(iiflag);
      averages(iflag) = w(ix^D,iflag)-wold(ix^D,iflag)
      average=average+wflags(iiflag)*abs(averages(iflag))
      if (abs(wold(ix^D,iflag))<smalldouble)then
         error=error+wflags(iiflag)* &
            abs(averages(iflag))/(abs(wold(ix^D,iflag))+epsilon)
      else
         error=error+wflags(iiflag)* &
            abs(averages(iflag))/(abs(wold(ix^D,iflag)))
      end if
   end do
   if (specialtol) then
      wtol(1:nw)=pw(igrid)%w(ix^D,1:nw)
      xtol(1:ndim)=px(igrid)%x(ix^D,1:ndim)
      call special_tolerance(wtol,xtol,tolerance,t)
   end if
   if (error >= tolerance) then
      refineflag(ix^D) = .true.
   else if (error <= tolratio(level)*tolerance) then
      coarsenflag(ix^D) = .true.
   end if
{end do\}

if (any(refineflag(ixM^T))) then
   if (level<mxnest) refine(igrid,mype)=.true.
end if
if (time_advance) then
   if (all(coarsenflag(ixM^T)).and.level>1) coarsen(igrid,mype)=.true.
end if

end subroutine compare1_grid
!=============================================================================
subroutine forcedrefine_grid(igrid,w)
use mod_forest, only: coarsen, refine, buffer
use mod_global_parameters

integer, intent(in) :: igrid
double precision, intent(in) :: w(ixG^T,nw)

integer :: level
integer :: my_refine, my_coarsen
double precision :: qt
logical, dimension(ixG^T) :: refineflag
!-----------------------------------------------------------------------------
level=node(plevel_,igrid)

! initialize to 0
my_refine   = 0
my_coarsen  = 0

if (time_advance) then
   qt=t+dt
else
   if (errorestimate==2) then
      qt=t+dt
   else
      qt=t
   end if
end if
   
call specialrefine_grid(igrid,level,ixG^LL,ixM^LL,qt,w,px(igrid)%x, &
                        my_refine,my_coarsen)

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

if (nbufferx^D/=0|.or.) then
   if (refine(igrid,mype) .and. .not.buffer(igrid,mype)) then
      refineflag(ixM^T)=.true.
      call refinebuffer(igrid,refineflag)
   end if
end if

end subroutine forcedrefine_grid
!=============================================================================
subroutine forcedrefine_grid_io(igrid,w)
use mod_forest, only: coarsen, refine
use mod_global_parameters

integer, intent(in)          :: igrid
double precision, intent(in) :: w(ixG^T,nw)

integer                   :: level, my_levmin, my_levmax
logical, dimension(ixG^T) :: refineflag
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
use mod_global_parameters

integer, intent(in) :: igrid
logical, dimension(ixG^T), intent(in) :: refineflag

integer :: ishiftbuf^D, i^D, ix^L, ineighbor, ipe_neighbor, level
!-----------------------------------------------------------------------------
ishiftbuf^D=ixMhi^D-ixMlo^D-nbufferx^D+1;
{do i^DB=-1,1\}
   ixmin^D=max(ixMlo^D,ixMlo^D+i^D*ishiftbuf^D);
   ixmax^D=min(ixMhi^D,ixMhi^D+i^D*ishiftbuf^D);
   if (ixmax^D<ixmin^D|.or.) cycle
   if (any(refineflag(ix^S))) then
      select case (neighbor_type(i^D,igrid))
      case (2)
         ineighbor=neighbor(1,i^D,igrid)
         ipe_neighbor=neighbor(2,i^D,igrid)
         if (.not.refine(ineighbor,ipe_neighbor)) then
            buffer(ineighbor,ipe_neighbor)=.true.
            refine(ineighbor,ipe_neighbor)=.true.
         end if
      case (3)
         level=node(plevel_,igrid)
         if (level<mxnest) then
            ineighbor=neighbor(1,i^D,igrid)
            ipe_neighbor=neighbor(2,i^D,igrid)
            if (.not.refine(ineighbor,ipe_neighbor)) then
               buffer(ineighbor,ipe_neighbor)=.true.
               refine(ineighbor,ipe_neighbor)=.true.
            end if
         end if
      end select
   end if
{end do\}

end subroutine refinebuffer
!=============================================================================
