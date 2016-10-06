!=============================================================================
subroutine write_analysis
! This is an example file how to use the analysis capability.  You can schedule 
! this routine using the slot 5 in itsave, dtsave and ditsave.  
! To use, just copy this file to your working directory and make your modifications.
use constants
use mod_global_parameters
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------
logical :: fileopen
integer :: iigrid, igrid, iw
double precision, dimension(1:nw+ndim) :: send_buff, recv_buff
double precision :: w(ixG^T,1:nw), rmin,rmax,zmin,zmax, myw(1:nw), myx(1:ndim)
character(len=80) :: filename
logical, save :: opened=.false.
logical, save :: file_exists=.false.
integer :: amode, ip, jp, ipe, ipehasvalue
double precision, parameter  :: rp=2.0d0, zp=0.2d0
logical                      :: mask(ixG^T), logbuff, havepoint
integer                      :: status(MPI_STATUS_SIZE)
!-----------------------------------------------------------------------------

havepoint=.false.
do iigrid=1,igridstail; igrid=igrids(iigrid);

! Check wether we have the pivot point:
   rmin = minval(px(igrid)%x(ixM^T,1))
   rmax = maxval(px(igrid)%x(ixM^T,1))
   zmin = minval(px(igrid)%x(ixM^T,2))
   zmax = maxval(px(igrid)%x(ixM^T,2))
   
   if (.not.(rp .gt. rmin .and. rp .le. rmax &
        .and. zp .gt. zmin .and. zp .le. zmax)) then 
      cycle
   else 
      ! Ok I have the point, everyone else cycled away.  
      havepoint = .true.

      
      mask(ixM^T)=.true.
      ip = minloc(abs(rp-px(igrid)%x(ixMlo1:ixMhi1,ixMlo2,1)),1,mask(ixMlo1:ixMhi1,ixMlo2))+dixB
      jp = minloc(abs(zp-px(igrid)%x(ixMlo1,ixMlo2:ixMhi2,2)),1,mask(ixMlo1,ixMlo2:ixMhi2))+dixB
      
      w(ixG^T,1:nw) = pw(igrid)%w(ixG^T,1:nw)
      call primitive(ixG^LL,ixM^LL,w,px(igrid)%x)
      
      ! Normalize to cgs units:
      do iw = 1, nw
         w(ip,jp,iw) = normvar(iw) * w(ip,jp,iw)
      end do
      
      send_buff(1:nw)=w(ip,jp,1:nw)
      {^D&send_buff(nw+^D)=px(igrid)%x(ip,jp,^D)\}
      if (mype .ne. 0) call MPI_SEND(send_buff,nw+ndim,MPI_DOUBLE_PRECISION,0,npe+1,icomm,ierrmpi)
   end if
end do
if (mype .ne. 0) call MPI_SEND(havepoint,1,MPI_LOGICAL,0,mype,icomm,ierrmpi)

if (mype == 0) then
! ============= Communicate =============
ipehasvalue=0
do ipe=1,npe-1
   call MPI_RECV(logbuff,1,MPI_LOGICAL,ipe,ipe,&
        icomm,status,ierrmpi)
   if (logbuff) ipehasvalue=ipe
end do
if (ipehasvalue .ne. 0) then
   call MPI_RECV(recv_buff,nw+ndim,MPI_DOUBLE_PRECISION,ipehasvalue,npe+1,icomm,status,ierrmpi)
   myw(1:nw) = recv_buff(1:nw)
 {^D& myx(^D) = recv_buff(nw+^D)\}
else
   myw(1:nw) = w(ip,jp,1:nw)
 {^D& myx(^D) = px(igrid)%x(ip,jp,^D)\}
end if
! ============= Done communicating =============

if (.not.opened) then
   ! generate filename
   write(filename,"(a,a)") TRIM("analysis"),".csv"
   INQUIRE(FILE=filename, EXIST=file_exists)
   
   if (.not. file_exists) then
      open(unit=unitanalysis,file=filename,status='unknown',access='append')
      write(unitanalysis,"(a)",advance="no") trim('# t [years] rho u1 u2 u3 p b1 b2 b3 tr2 lfac xi evolve at (r,z)=(')
      write(unitanalysis,'(2(es14.6))',advance="no") myx(1)*normvar(0),myx(2)*normvar(0)
      write(unitanalysis,"(a)") trim(') [cm]')
   else
      open(unit=unitanalysis,file=filename,status='unknown',access='append')
   end if
   opened=.true.
end if
   
write(unitanalysis,'(13(es14.6))')t*UNIT_LENGTH/UNIT_VELOCITY/CONST_years, &
     (myw(iw), iw=1,nw)

end if ! mype == 0

call MPI_BARRIER(icomm,ierrmpi)

end subroutine write_analysis
!=============================================================================
