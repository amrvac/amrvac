!=============================================================================
subroutine init_comm_fix_conserve(idim^LIM)
use mod_fix_conserve
use mod_global_parameters

integer, intent(in) :: idim^LIM

integer :: iigrid, igrid, idims, iside, i^D, nxCo^D
integer :: ic^D, inc^D, ipe_neighbor
integer :: ibuf, recvsize
!-----------------------------------------------------------------------------
nsend=0
nrecv=0
recvsize=0

do idims= idim^LIM
   select case (idims)
   {case (^D)
      nrecv=nrecv+nrecv_fc(^D)
      nsend=nsend+nsend_fc(^D)
      nxCo^D=1^D%nxCo^DD=ixGhi^DD/2-dixB;
      isize(^D)={nxCo^DD*}*(nwflux)
      recvsize=recvsize+nrecv_fc(^D)*isize(^D)\}
   end select
end do

if (nrecv>0) then
   allocate(recvbuffer(recvsize),recvstatus(MPI_STATUS_SIZE,nrecv), &
            recvrequest(nrecv))

   recvrequest=MPI_REQUEST_NULL
   ibuf=1
   irecv=0

   do iigrid=1,igridstail; igrid=igrids(iigrid);
      do idims= idim^LIM
         do iside=1,2
            i^D=kr(^D,idims)*(2*iside-3);
            {^IFPHI if (neighbor_pole(i^D,igrid)/=0) cycle}
            if (neighbor_type(i^D,igrid)/=4) cycle
            {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
               inc^DB=2*i^DB+ic^DB\}
               ipe_neighbor=neighbor_child(2,inc^D,igrid)
               if (ipe_neighbor/=mype) then
                  irecv=irecv+1
                  itag=4**^ND*(igrid-1)+{inc^D*4**(^D-1)+}
                  call MPI_IRECV(recvbuffer(ibuf),isize(idims), &
                                 MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                                 icomm,recvrequest(irecv),ierrmpi)
                  ibuf=ibuf+isize(idims)
               end if
            {end do\}
         end do
      end do
   end do
end if

if (nsend>0) then
   allocate(sendstatus(MPI_STATUS_SIZE,nsend),sendrequest(nsend))
   sendrequest=MPI_REQUEST_NULL
   isend=0
end if

end subroutine init_comm_fix_conserve
!=============================================================================
subroutine allocateBflux
use mod_fix_conserve
use mod_global_parameters

integer :: iigrid, igrid, iside, i^D, nx^D, nxCo^D
!-----------------------------------------------------------------------------
nx^D=ixMhi^D-ixMlo^D+1;
nxCo^D=nx^D/2;

do iigrid=1,igridstail; igrid=igrids(iigrid);
   {do iside=1,2
      i^DD=kr(^DD,^D)*(2*iside-3);
      {^IFPHI if (neighbor_pole(i^DD,igrid)/=0) cycle}
      select case (neighbor_type(i^DD,igrid))
      case (4)
         allocate(pflux(iside,^D,igrid)%flux(1^D%1:nx^DD,1:nwflux))
      case (2)
         allocate(pflux(iside,^D,igrid)%flux(1^D%1:nxCo^DD,1:nwflux))
      end select
   end do\}
end do

end subroutine allocateBflux
!=============================================================================
subroutine deallocateBflux
use mod_fix_conserve
use mod_global_parameters

integer :: iigrid, igrid, iside
!-----------------------------------------------------------------------------
do iigrid=1,igridstail; igrid=igrids(iigrid);
   {do iside=1,2
      if (associated(pflux(iside,^D,igrid)%flux)) then
         deallocate(pflux(iside,^D,igrid)%flux)
         nullify(pflux(iside,^D,igrid)%flux)
      end if
   end do\}
end do

end subroutine deallocateBflux
!=============================================================================
subroutine fix_conserve(pwuse,idim^LIM)
use mod_fix_conserve
use mod_global_parameters

type(walloc) :: pwuse(ngridshi)
integer, intent(in) :: idim^LIM

integer :: iigrid, igrid, idims, iside, iotherside, i^D, ic^D, inc^D, ix^L
integer :: nxCo^D, iw, ix, ipe_neighbor, ineighbor, ibuf, ibufnext
double precision :: CoFiratio
!-----------------------------------------------------------------------------
if (slab) then
   ! The flux is divided by volume of fine cell. We need, however,
   ! to divide by volume of coarse cell => muliply by volume ratio
   CoFiratio=one/dble(2**ndim)
end if

if (nrecv>0) then
   call MPI_WAITALL(nrecv,recvrequest,recvstatus,ierrmpi)
   ibuf=1
end if

nxCo^D=(ixMhi^D-ixMlo^D+1)/2;

!if (.false.) then 

! for all grids: perform flux update at Coarse-Fine interfaces
do iigrid=1,igridstail; igrid=igrids(iigrid);
   do idims= idim^LIM
      select case (idims)
      {case (^D)
         do iside=1,2
            i^DD=kr(^DD,^D)*(2*iside-3);
            {^IFPHI if (neighbor_pole(i^DD,igrid)/=0) cycle}
            if (neighbor_type(i^DD,igrid)/=4) cycle

! opedit: skip over active/passive interface since flux for passive ones is 
! not computed, keep the buffer counter up to date:
            if (.not.neighbor_active(i^DD,igrid).or.&
                .not.neighbor_active(0^DD&,igrid) ) then
               {do ic^DDB=1+int((1-i^DDB)/2),2-int((1+i^DDB)/2)
               inc^DDB=2*i^DDB+ic^DDB\}
               ipe_neighbor=neighbor_child(2,inc^DD,igrid)
               if (ipe_neighbor/=mype) then
                  ibufnext=ibuf+isize(^D)
                  ibuf=ibufnext
                  end if
               {end do\}
               cycle
            end if
!

            select case (iside)
            case (1)
               ix=ixMlo^D
            case (2)
               ix=ixMhi^D
            end select

            ! remove coarse flux
            if (slab) then
               pwuse(igrid)%w(ix^D%ixM^T,1:nwflux) &
                  = pwuse(igrid)%w(ix^D%ixM^T,1:nwflux) &
                   -pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwflux)
            else
               do iw=1,nwflux
                  pwuse(igrid)%w(ix^D%ixM^T,iw)=pwuse(igrid)%w(ix^D%ixM^T,iw)&
                     -pflux(iside,^D,igrid)%flux(1^D%:^DD&,iw) &
                     /pgeo(igrid)%dvolume(ix^D%ixM^T)
               end do
            end if


            ! add fine flux
            {do ic^DDB=1+int((1-i^DDB)/2),2-int((1+i^DDB)/2)
               inc^DDB=2*i^DDB+ic^DDB\}
               ineighbor=neighbor_child(1,inc^DD,igrid)
               ipe_neighbor=neighbor_child(2,inc^DD,igrid)
               ixmin^D=ix^D%ixmin^DD=ixMlo^DD+(ic^DD-1)*nxCo^DD;
               ixmax^D=ix^D%ixmax^DD=ixmin^DD-1+nxCo^DD;
               if (ipe_neighbor==mype) then
                  iotherside=3-iside
                  if (slab) then
                     pwuse(igrid)%w(ix^S,1:nwflux) &
                       = pwuse(igrid)%w(ix^S,1:nwflux) &
                       + pflux(iotherside,^D,ineighbor)%flux(:^DD&,1:nwflux)&
                       * CoFiratio
                  else
                     do iw=1,nwflux
                        pwuse(igrid)%w(ix^S,iw)=pwuse(igrid)%w(ix^S,iw) &
                            +pflux(iotherside,^D,ineighbor)%flux(:^DD&,iw) &
                            /pgeo(igrid)%dvolume(ix^S)
                     end do
                  end if
               else
                  if (slab) then
                     ibufnext=ibuf+isize(^D)
                     pwuse(igrid)%w(ix^S,1:nwflux) &
                         = pwuse(igrid)%w(ix^S,1:nwflux)+CoFiratio &
                          *reshape(source=recvbuffer(ibuf:ibufnext-1), &
                                 shape=shape(pwuse(igrid)%w(ix^S,1:nwflux)))
                     ibuf=ibufnext
                  else
                     do iw=1,nwflux
                        ibufnext=ibuf+isize(^D)/(nwflux)
                        pwuse(igrid)%w(ix^S,iw)=pwuse(igrid)%w(ix^S,iw) &
                           +reshape(source=recvbuffer(ibuf:ibufnext-1), &
                                    shape=shape(pwuse(igrid)%w(ix^S,iw))) &
                           /pgeo(igrid)%dvolume(ix^S)
                        ibuf=ibufnext
                     end do
                  end if
               end if
            {end do\}
         end do\}
      end select
   end do
end do

!end if

if (nrecv>0) deallocate(recvbuffer,recvstatus,recvrequest)

if (nsend>0) then
   call MPI_WAITALL(nsend,sendrequest,sendstatus,ierrmpi)
   deallocate(sendstatus,sendrequest)
end if

end subroutine fix_conserve
!=============================================================================
subroutine storeflux(igrid,fC,idim^LIM)
use mod_fix_conserve
use mod_global_parameters

integer, intent(in)          :: igrid, idim^LIM
double precision, intent(in) :: fC(ixG^T,1:nwflux,1:ndim)

integer :: idims, iside, i^D, ic^D, inc^D, ix^D, ixCo^D, nxCo^D, iw
!integer :: ineighbor, ipe_neighbor
!-----------------------------------------------------------------------------
do idims = idim^LIM
   select case (idims)
   {case (^D)
      do iside=1,2
         i^DD=kr(^DD,^D)*(2*iside-3);
         {^IFPHI if (neighbor_pole(i^DD,igrid)/=0) cycle}
         select case (neighbor_type(i^DD,igrid))
         case (4)
            select case (iside)
            case (1)
               pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwflux) = &
                  -fC(dixB^D%ixM^T,1:nwflux,^D)
            case (2)
               pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwflux) = &
                  fC(ixMhi^D^D%ixM^T,1:nwflux,^D)
            end select
         case (2)
            nxCo^D=1^D%nxCo^DD=ixGhi^DD/2-dixB;
            select case (iside)
            case (1)
               do iw=1,nwflux
                  {do ixCo^DDB=1,nxCo^DDB\}
                     ix^D=dixB^D%ix^DD=ixMlo^DD+2*(ixCo^DD-1);
                     pflux(iside,^D,igrid)%flux(ixCo^DD,iw) &
                        = {^NOONEDsum}(fC(ix^D^D%ix^DD:ix^DD+1,iw,^D))
                  {end do\}
               end do
            case (2)
               do iw=1,nwflux
                  {do ixCo^DDB=1,nxCo^DDB\}
                     ix^D=ixMhi^D^D%ix^DD=ixMlo^DD+2*(ixCo^DD-1);
                     pflux(iside,^D,igrid)%flux(ixCo^DD,iw) &
                        =-{^NOONEDsum}(fC(ix^D^D%ix^DD:ix^DD+1,iw,^D))
                  {end do\}
               end do
            end select

!            ineighbor=neighbor(1,i^DD,igrid)
!            ipe_neighbor=neighbor(2,i^DD,igrid)
!            if (ipe_neighbor/=mype) then
!               ic^DD=1+modulo(node(pig^DD_,igrid)-1,2);
!               inc^D=-2*i^D+ic^D^D%inc^DD=ic^DD;
!               itag=4**^ND*(ineighbor-1)+{inc^DD*4**(^DD-1)+}
!               isend=isend+1
!               call MPI_ISEND(pflux(iside,^D,igrid)%flux,isize(^D), &
!                              MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
!                              icomm,sendrequest(isend),ierrmpi)
!            end if

         end select
      end do\}
   end select
end do

end subroutine storeflux
!=============================================================================
subroutine sendflux(igrid,idim^LIM)
use mod_fix_conserve
use mod_global_parameters

integer, intent(in)          :: igrid, idim^LIM

integer :: idims, iside, i^D, ic^D, inc^D, ix^D, ixCo^D, nxCo^D, iw
integer :: ineighbor, ipe_neighbor
!----------------------------------------------------------------------------
do idims = idim^LIM
   select case (idims)
   {case (^D)
      do iside=1,2
         i^DD=kr(^DD,^D)*(2*iside-3);
         {^IFPHI if (neighbor_pole(i^DD,igrid)/=0) cycle}
         
         if (neighbor_type(i^DD,igrid)==2) then

            ineighbor=neighbor(1,i^DD,igrid)
            ipe_neighbor=neighbor(2,i^DD,igrid)
            if (ipe_neighbor/=mype) then
               ic^DD=1+modulo(node(pig^DD_,igrid)-1,2);
               inc^D=-2*i^D+ic^D^D%inc^DD=ic^DD;
               itag=4**^ND*(ineighbor-1)+{inc^DD*4**(^DD-1)+}
               isend=isend+1

               call MPI_ISEND(pflux(iside,^D,igrid)%flux,isize(^D), &
                              MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                              icomm,sendrequest(isend),ierrmpi)
            end if
         end if
      end do\}
   end select
end do
end subroutine sendflux
!=============================================================================
