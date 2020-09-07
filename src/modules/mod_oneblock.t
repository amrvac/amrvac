!=============================================================================
! Module oneblock to hold a datastructure as output by convert_type='oneblock'. 
! Currently only ASCII
! Can be handy to read in initial conditions.
! 2015-02-18 Oliver Porth
!=============================================================================
module mod_oneblock
implicit none
save

double precision, dimension({^D&:|,},:), allocatable       :: woneblock
double precision, dimension({^D&:|,},:), allocatable       :: xoneblock
integer                                                    :: nc^D
integer                                                    :: unit=15 !file unit to read on
!-----------------------------------------------------------------------------


contains
!=============================================================================
subroutine read_oneblock(filename)

use mod_global_parameters

character(len=*), intent(in)             :: filename
! .. local ..
character(len=1024)                      :: outfilehead
integer                                  :: nctot, ix^D, ixp^D
double precision                         :: time
integer                                  :: idim
integer,dimension(^ND)                   :: sendbuff
!-----------------------------------------------------------------------------

if (mype == 0) write(*,*) 'mod_oneblock: reading ',nw,' variables on unit:', unit

!----------------------------------------
! Root does the reading:
!----------------------------------------
if (mype == 0) then 
   open(unit,file=filename,status='unknown')
   ! The header information:
   read(unit,'(A)') outfilehead
   read(unit,*) nctot,nc^D
   read(unit,*) time
   
   ! Allocate and read the grid and variables:
   allocate(xoneblock(nc^D,1:^ND))
   allocate(woneblock(nc^D,1:nw))
   
   {do ix^DB=1,nc^DB\}  
   
   read(unit,*) xoneblock(ix^D,1:^ND), woneblock(ix^D,1:nw)
   
   {end do\}
   
! Close the file
   close(unit)
end if! mype==0

!----------------------------------------
! Broadcast what mype=0 read:
!----------------------------------------
if (npe>1) then 
   {sendbuff(^D)=nc^D;}
   call MPI_BCAST(sendbuff,^ND,MPI_INTEGER,0,icomm,ierrmpi)
   if (mype .ne. 0) then 
      {nc^D=sendbuff(^D);}
      ! Allocate the grid and variables:
      allocate(xoneblock(nc^D,1:^ND))
      allocate(woneblock(nc^D,1:nw))
   end if
   call MPI_BCAST(xoneblock,{nc^D|*}*^ND,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
   call MPI_BCAST(woneblock,{nc^D|*}*nw,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
end if! npe>1



end subroutine read_oneblock
!=============================================================================
subroutine interpolate_oneblock(x,iw,out)

double precision, dimension(^ND),intent(in)             :: x 
integer, intent(in)                                     :: iw
double precision, intent(out)                           :: out
! .. local ..
double precision, dimension(^ND)                        :: xloc
integer                                                 :: ic^D, ic1^D, ic2^D
double precision                                        :: xd^D
{^IFTWOD
double precision                                        :: c00, c10
}
{^IFTHREED
double precision                                        :: c0, c1, c00, c10, c01, c11
}
integer                                                 :: ipivot^D, idir
!-----------------------------------------------------------------------------

xloc=x

!--------------------------------------------
! Hunt for the index closest to the point
! This is a bit slow but allows for stretched grids
! (still need to be orthogonal for interpolation though)
!--------------------------------------------
ipivot^D=1;
{^IFONED
ic1 = minloc(dabs(xloc(1)-xoneblock(:,1)),1,mask=.true.)
}
{^IFTWOD
do idir = 1, ^ND
   select case (idir)
   case (1)
      ic1 = minloc(dabs(xloc(1)-xoneblock(:,ipivot2,idir)),1,mask=.true.)
   case (2)
      ic2 = minloc(dabs(xloc(2)-xoneblock(ipivot1,:,idir)),1,mask=.true.)
   case default
      call mpistop("error1 in interpolate_oneblock")
   end select
end do
}
{^IFTHREED
do idir = 1, ^ND
   select case (idir)
   case (1)
      ic1 = minloc(dabs(xloc(1)-xoneblock(:,ipivot2,ipivot3,idir)),1, &
           mask=.true.)
   case (2)
      ic2 = minloc(dabs(xloc(2)-xoneblock(ipivot1,:,ipivot3,idir)),1, &
           mask=.true.)
   case (3)
      ic3 = minloc(dabs(xloc(3)-xoneblock(ipivot1,ipivot2,:,idir)),1, &
           mask=.true.)
   case default
      call mpistop("error1 in interpolate_oneblock")
   end select
end do
}

! flat interpolation would simply be:
!out = woneblock(ic^D,iw)
!return

!-------------------------------------------
! Get the left and right indices
!-------------------------------------------
{
if (xoneblock({ic^DD},^D) .lt. xloc(^D)) then
   ic1^D = ic^D
else
   ic1^D = ic^D -1
end if
ic2^D = ic1^D + 1
\}

!--------------------------------------------
! apply flat interpolation if outside of range, 
! change point-location to make this easy!
!--------------------------------------------
{
if (ic1^D .lt. 1) then
   ic1^D = 1
   ic2^D = ic1^D + 1
   xloc(^D) = xoneblock(ic^DD,^D)
end if
\}
{
if (ic2^D .gt. nc^D) then
   ic2^D = nc^D
   ic1^D = ic2^D - 1
   xloc(^D) = xoneblock(ic^DD,^D)
end if
\}


!-------------------------------------------
! linear, bi- and tri- linear interpolations
!-------------------------------------------
{^IFONED
xd1 = (xloc(1)-xoneblock(ic11,1)) / (xoneblock(ic21,1) - xoneblock(ic11,1))
out = woneblock(ic11,iw) * (1.0d0 - xd1) + woneblock(ic21,iw) * xd1
}
{^IFTWOD
xd1 = (xloc(1)-xoneblock(ic11,ic12,1)) / (xoneblock(ic21,ic12,1) - xoneblock(ic11,ic12,1))      
xd2 = (xloc(2)-xoneblock(ic11,ic12,2)) / (xoneblock(ic11,ic22,2) - xoneblock(ic11,ic12,2))
c00 = woneblock(ic11,ic12,iw) * (1.0d0 - xd1) + woneblock(ic21,ic12,iw) * xd1
c10 = woneblock(ic11,ic22,iw) * (1.0d0 - xd1) + woneblock(ic21,ic22,iw) * xd1
out = c00 * (1.0d0 - xd2) + c10 * xd2
}
{^IFTHREED
xd1 = (xloc(1)-xoneblock(ic11,ic12,ic13,1)) / (xoneblock(ic21,ic12,ic13,1) - xoneblock(ic11,ic12,ic13,1))      
xd2 = (xloc(2)-xoneblock(ic11,ic12,ic13,2)) / (xoneblock(ic11,ic22,ic13,2) - xoneblock(ic11,ic12,ic13,2))      
xd3 = (xloc(3)-xoneblock(ic11,ic12,ic13,3)) / (xoneblock(ic11,ic12,ic23,3) - xoneblock(ic11,ic12,ic13,3))    

c00 = woneblock(ic11,ic12,ic13,iw) * (1.0d0 - xd1) + woneblock(ic21,ic12,ic13,iw) * xd1
c10 = woneblock(ic11,ic22,ic13,iw) * (1.0d0 - xd1) + woneblock(ic21,ic22,ic13,iw) * xd1
c01 = woneblock(ic11,ic12,ic23,iw) * (1.0d0 - xd1) + woneblock(ic21,ic12,ic23,iw) * xd1
c11 = woneblock(ic11,ic22,ic23,iw) * (1.0d0 - xd1) + woneblock(ic21,ic22,ic23,iw) * xd1

c0  = c00 * (1.0d0 - xd2) + c10 * xd2
c1  = c01 * (1.0d0 - xd2) + c11 * xd2

out = c0 * (1.0d0 - xd3) + c1 * xd3
}

end subroutine interpolate_oneblock
!=============================================================================
end module mod_oneblock
!=============================================================================
