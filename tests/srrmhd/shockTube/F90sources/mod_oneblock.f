!=============================================================================
! Module oneblock to hold a datastructure as output by convert_type='oneblock'. 
! Currently only ASCII
! Can be handy to read in initial conditions.
! 2015-02-18 Oliver Porth
!=============================================================================
module mod_oneblock
implicit none
save

double precision, dimension(:,:), allocatable       :: woneblock
double precision, dimension(:,:), allocatable       :: xoneblock
integer                                                    :: nc1
integer                                                    :: unit=15 !file unit to read on
!-----------------------------------------------------------------------------


contains
!=============================================================================
subroutine read_oneblock(filename)

include 'amrvacdef.f'

character(len=*), intent(in)             :: filename
! .. local ..
character(len=1024)                      :: outfilehead
integer                                  :: nctot, ix1, ixp1
double precision                         :: time
integer                                  :: idim
integer,dimension(1)                   :: sendbuff
!-----------------------------------------------------------------------------

if (mype == 0) write(*,*), 'mod_oneblock: reading ',nw,' variables on unit:',&
    unit

!----------------------------------------
! Root does the reading:
!----------------------------------------
if (mype == 0) then 
   open(unit,file=filename,status='unknown')
   ! The header information:
   read(unit,'(A)') outfilehead
   read(unit,*) nctot,nc1
   read(unit,*) time
   
   ! Allocate and read the grid and variables:
   allocate(xoneblock(nc1,1:1))
   allocate(woneblock(nc1,1:nw))
   
   do ix1=1,nc1  
   
   read(unit,*) xoneblock(ix1,1:1), woneblock(ix1,1:nw)
   
   end do
   
! Close the file
   close(unit)
end if! mype==0

!----------------------------------------
! Broadcast what mype=0 read:
!----------------------------------------
if (npe>1) then 
   sendbuff(1)=nc1;
   call MPI_BCAST(sendbuff,1,MPI_INTEGER,0,icomm,ierrmpi)
   if (mype .ne. 0) then 
      nc1=sendbuff(1);
      ! Allocate the grid and variables:
      allocate(xoneblock(nc1,1:1))
      allocate(woneblock(nc1,1:nw))
   end if
   call MPI_BCAST(xoneblock,nc1*1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
   call MPI_BCAST(woneblock,nc1*nw,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
end if! npe>1



end subroutine read_oneblock
!=============================================================================
subroutine interpolate_oneblock(x,iw,out)

double precision, dimension(1),intent(in)             :: x 
integer, intent(in)                                     :: iw
double precision, intent(out)                           :: out
! .. local ..
double precision, dimension(1)                        :: xloc
integer                                                 :: ic1, ic11, ic21
double precision                                        :: xd1


integer                                                 :: ipivot1, idir
!-----------------------------------------------------------------------------

xloc=x

!--------------------------------------------
! Hunt for the index closest to the point
! This is a bit slow but allows for stretched grids
! (still need to be orthogonal for interpolation though)
!--------------------------------------------
ipivot1=1;

ic1 = minloc(dabs(xloc(1)-xoneblock(:,1)),1,mask=.true.)




! flat interpolation would simply be:
!out = woneblock(ic^D,iw)
!return

!-------------------------------------------
! Get the left and right indices
!-------------------------------------------

if (xoneblock(ic1,1) .lt. xloc(1)) then
   ic11 = ic1
else
   ic11 = ic1 -1
end if
ic21 = ic11 + 1


!--------------------------------------------
! apply flat interpolation if outside of range, 
! change point-location to make this easy!
!--------------------------------------------

if (ic11 .lt. 1) then
   ic11 = 1
   ic21 = ic11 + 1
   xloc(1) = xoneblock(ic1,1)
end if


if (ic21 .gt. nc1) then
   ic21 = nc1
   ic11 = ic21 - 1
   xloc(1) = xoneblock(ic1,1)
end if



!-------------------------------------------
! linear, bi- and tri- linear interpolations
!-------------------------------------------

xd1 = (xloc(1)-xoneblock(ic11,1)) / (xoneblock(ic21,1) - xoneblock(ic11,1))
out = woneblock(ic11,iw) * (1.0d0 - xd1) + woneblock(ic21,iw) * xd1




end subroutine interpolate_oneblock
!=============================================================================
end module mod_oneblock
!=============================================================================
