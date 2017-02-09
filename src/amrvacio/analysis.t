subroutine write_analysis
! This is an example file how to use the analysis capability.  You can schedule 
! this routine using the slot 5 in itsave, dtsave and ditsave.  
! To use, just copy this file to your working directory and make your modifications.
use constants
use mod_global_parameters
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------
logical :: fileopen
integer :: iigrid, igrid, level, iw
double precision :: volume(1:nlevelshi),  voltotal
double precision :: re, ke, me, te
double precision :: dvolume(ixG^T), inBubble(ixG^T), lfac_pw(ixG^T), volumeflat(1:nlevelshi)
integer :: numlevels
integer, dimension(1:nlevelshi) :: isum_send
double precision, dimension(1:4) :: dsum_send, dsum_recv
character(len=80) :: filename
character(len=1024) :: line
logical, save :: opened=.false.
logical, save :: file_exists=.false.
integer :: amode, status(MPI_STATUS_SIZE)
double precision :: trcut

end subroutine write_analysis

