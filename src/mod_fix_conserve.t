module mod_fix_conserve
  implicit none
  private

  integer, save                        :: nrecv, nsend
  double precision, allocatable, save  :: recvbuffer(:)
  integer, save                        :: isize(^ND)
  integer, dimension(:), allocatable   :: fc_recvreq, fc_sendreq
  integer, dimension(:,:), allocatable :: fc_recvstat, fc_sendstat

  public :: init_comm_fix_conserve
  public :: allocateBflux
  public :: deallocateBflux
  public :: sendflux
  public :: recvflux
  public :: storeflux
  public :: fix_conserve

 contains

   subroutine init_comm_fix_conserve(idim^LIM)
     use mod_global_parameters

     integer, intent(in) :: idim^LIM

     integer :: iigrid, igrid, idims, iside, i^D, nxCo^D
     integer :: ic^D, inc^D, ipe_neighbor
     integer :: ibuf, recvsize

     nsend    = 0
     nrecv    = 0
     recvsize = 0

     do idims= idim^LIM
       select case (idims)
         {case (^D)
         nrecv=nrecv+nrecv_fc(^D)
         nsend=nsend+nsend_fc(^D)
         nxCo^D=1^D%nxCo^DD=ixGhi^DD/2-nghostcells;
         isize(^D)={nxCo^DD*}*(nwflux)
         recvsize=recvsize+nrecv_fc(^D)*isize(^D)\}
       end select
     end do

     ! Reallocate buffers when size differs
     if (allocated(recvbuffer)) then
       if (recvsize /= size(recvbuffer)) then
         deallocate(recvbuffer)
         allocate(recvbuffer(recvsize))
       end if
     else
       allocate(recvbuffer(recvsize))
     end if

     if (allocated(fc_recvreq)) then
       if (nrecv /= size(fc_recvreq)) then
         deallocate(fc_recvreq, fc_recvstat)
         allocate(fc_recvstat(MPI_STATUS_SIZE,nrecv), fc_recvreq(nrecv))
       end if
     else
       allocate(fc_recvstat(MPI_STATUS_SIZE,nrecv), fc_recvreq(nrecv))
     end if

     if (allocated(fc_sendreq)) then
       if (nrecv /= size(fc_sendreq)) then
         deallocate(fc_sendreq, fc_sendstat)
         allocate(fc_sendstat(MPI_STATUS_SIZE,nsend), fc_sendreq(nsend))
       end if
     else
       allocate(fc_sendstat(MPI_STATUS_SIZE,nsend), fc_sendreq(nsend))
     end if

   end subroutine init_comm_fix_conserve

   subroutine recvflux(idim^LIM)
     use mod_global_parameters

     integer, intent(in) :: idim^LIM

     integer :: iigrid, igrid, idims, iside, i^D, nxCo^D
     integer :: ic^D, inc^D, ipe_neighbor
     integer :: ibuf, recvsize

     if (nrecv>0) then
       fc_recvreq=MPI_REQUEST_NULL
       ibuf=1
       irecv=0

       do iigrid=1,igridstail; igrid=igrids(iigrid);
         do idims= idim^LIM
           do iside=1,2
             i^D=kr(^D,idims)*(2*iside-3);

             if (phi_ > 0) then
               if (neighbor_pole(i^D,igrid)/=0) cycle
             end if

             if (neighbor_type(i^D,igrid)/=4) cycle
             {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
             inc^DB=2*i^DB+ic^DB\}
             ipe_neighbor=neighbor_child(2,inc^D,igrid)
             if (ipe_neighbor/=mype) then
               irecv=irecv+1
               itag=4**^ND*(igrid-1)+{inc^D*4**(^D-1)+}
               call MPI_IRECV(recvbuffer(ibuf),isize(idims), &
                    MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                    icomm,fc_recvreq(irecv),ierrmpi)
               ibuf=ibuf+isize(idims)
             end if
             {end do\}
           end do
         end do
       end do
     end if

   end subroutine recvflux

   subroutine allocateBflux
     use mod_global_parameters

     integer :: iigrid, igrid, iside, i^D, nx^D, nxCo^D
     !-----------------------------------------------------------------------------
     nx^D=ixMhi^D-ixMlo^D+1;
     nxCo^D=nx^D/2;

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       {do iside=1,2
         i^DD=kr(^DD,^D)*(2*iside-3);

         if (phi_ > 0) then
           if (neighbor_pole(i^DD,igrid)/=0) cycle
         end if

         select case (neighbor_type(i^DD,igrid))
         case (4)
           allocate(pflux(iside,^D,igrid)%flux(1^D%1:nx^DD,1:nwflux))
         case (2)
           allocate(pflux(iside,^D,igrid)%flux(1^D%1:nxCo^DD,1:nwflux))
         end select
       end do\}
     end do

   end subroutine allocateBflux

   subroutine deallocateBflux
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

   subroutine fix_conserve(idim^LIM)
     use mod_global_parameters

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
       call MPI_WAITALL(nrecv,fc_recvreq,fc_recvstat,ierrmpi)
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

             if (phi_ > 0) then
               if (neighbor_pole(i^DD,igrid)/=0) cycle
             end if

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
               pw(igrid)%wb(ix^D%ixM^T,1:nwflux) &
                    = pw(igrid)%wb(ix^D%ixM^T,1:nwflux) &
                    -pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwflux)
             else
               do iw=1,nwflux
                 pw(igrid)%wb(ix^D%ixM^T,iw)=pw(igrid)%wb(ix^D%ixM^T,iw)&
                      -pflux(iside,^D,igrid)%flux(1^D%:^DD&,iw) &
                      /pw(igrid)%dvolume(ix^D%ixM^T)
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
                 pw(igrid)%wb(ix^S,1:nwflux) &
                      = pw(igrid)%wb(ix^S,1:nwflux) &
                      + pflux(iotherside,^D,ineighbor)%flux(:^DD&,1:nwflux)&
                      * CoFiratio
               else
                 do iw=1,nwflux
                   pw(igrid)%wb(ix^S,iw)=pw(igrid)%wb(ix^S,iw) &
                        +pflux(iotherside,^D,ineighbor)%flux(:^DD&,iw) &
                        /pw(igrid)%dvolume(ix^S)
                 end do
               end if
             else
               if (slab) then
                 ibufnext=ibuf+isize(^D)
                 pw(igrid)%wb(ix^S,1:nwflux) &
                      = pw(igrid)%wb(ix^S,1:nwflux)+CoFiratio &
                      *reshape(source=recvbuffer(ibuf:ibufnext-1), &
                      shape=shape(pw(igrid)%wb(ix^S,1:nwflux)))
                 ibuf=ibufnext
               else
                 do iw=1,nwflux
                   ibufnext=ibuf+isize(^D)/(nwflux)
                   pw(igrid)%wb(ix^S,iw)=pw(igrid)%wb(ix^S,iw) &
                        +reshape(source=recvbuffer(ibuf:ibufnext-1), &
                        shape=shape(pw(igrid)%wb(ix^S,iw))) &
                        /pw(igrid)%dvolume(ix^S)
                   ibuf=ibufnext
                 end do
               end if
             end if
             {end do\}
           end do\}
         end select
       end do
     end do

     if (nsend>0) then
       call MPI_WAITALL(nsend,fc_sendreq,fc_sendstat,ierrmpi)
     end if

   end subroutine fix_conserve

   subroutine storeflux(igrid,fC,idim^LIM)
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

           if (phi_ > 0) then
             if (neighbor_pole(i^DD,igrid)/=0) cycle
           end if

           select case (neighbor_type(i^DD,igrid))
           case (4)
             select case (iside)
             case (1)
               pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwflux) = &
                    -fC(nghostcells^D%ixM^T,1:nwflux,^D)
             case (2)
               pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwflux) = &
                    fC(ixMhi^D^D%ixM^T,1:nwflux,^D)
             end select
           case (2)
             nxCo^D=1^D%nxCo^DD=ixGhi^DD/2-nghostcells;
             select case (iside)
             case (1)
               do iw=1,nwflux
                 {do ixCo^DDB=1,nxCo^DDB\}
                 ix^D=nghostcells^D%ix^DD=ixMlo^DD+2*(ixCo^DD-1);
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
             !                              icomm,fc_sendreq(isend),ierrmpi)
             !            end if

           end select
         end do\}
       end select
     end do

   end subroutine storeflux

   subroutine sendflux(idim^LIM)
     use mod_global_parameters

     integer, intent(in) :: idim^LIM

     integer :: idims, iside, i^D, ic^D, inc^D, ix^D, ixCo^D, nxCo^D, iw
     integer :: ineighbor, ipe_neighbor, igrid, iigrid

     fc_sendreq = MPI_REQUEST_NULL
     isend      = 0

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims = idim^LIM
         select case (idims)
           {case (^D)
           do iside=1,2
             i^DD=kr(^DD,^D)*(2*iside-3);

             if (phi_ > 0) then
               if (neighbor_pole(i^DD,igrid)/=0) cycle
             end if

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
                      icomm,fc_sendreq(isend),ierrmpi)
               end if
             end if
           end do\}
         end select
       end do
     end do
   end subroutine sendflux

end module mod_fix_conserve
