!> Module for flux conservation near refinement boundaries
module mod_fix_conserve
  implicit none
  private

  type fluxalloc
     double precision, dimension(:^D&,:), pointer:: flux => null()
     double precision, dimension(:^D&,:), pointer:: edge => null()
  end type fluxalloc
  !> store flux to fix conservation
  type(fluxalloc), dimension(:,:,:), allocatable, public :: pflux

  integer, save                        :: nrecv, nsend
  double precision, allocatable, save  :: recvbuffer(:), sendbuffer(:)
  integer, dimension(:), allocatable   :: fc_recvreq, fc_sendreq
  integer, dimension(:,:), allocatable :: fc_recvstat, fc_sendstat
  integer, dimension(^ND), save        :: isize
  integer                              :: ibuf, ibuf_send
  ! ct for corner total
  integer, save                        :: nrecv_ct, nsend_ct
  ! buffer for corner coarse
  double precision, allocatable, save  :: recvbuffer_cc(:), sendbuffer_cc(:)
  integer, dimension(:), allocatable   :: cc_recvreq, cc_sendreq
  integer, dimension(:,:), allocatable :: cc_recvstat, cc_sendstat
  integer, dimension(^ND), save        :: isize_stg
  integer                              :: ibuf_cc, ibuf_cc_send
  integer                              :: itag_cc, isend_cc

  public :: init_comm_fix_conserve
  public :: allocateBflux
  public :: deallocateBflux
  public :: sendflux
  public :: recvflux
  public :: storeflux
  public :: fix_conserve

 contains

   subroutine init_comm_fix_conserve(idim^LIM,nwfluxin)
     use mod_global_parameters

     integer, intent(in) :: idim^LIM,nwfluxin

     integer :: iigrid, igrid, idims, iside, i^D, nxCo^D
     integer :: ic^D, inc^D, ipe_neighbor
     integer :: recvsize, sendsize
     integer :: recvsize_cc, sendsize_cc

     nsend    = 0
     nrecv    = 0
     recvsize = 0
     sendsize = 0
     if(stagger_grid) then
       ! Special communication for diagonal 'coarse corners'
       ! nrecv/send_cc (for 'coarse corners' is a dim=ndim-1 array which
       ! stores the faces that must be communicated in each direction.
       ! nrecv/send_ct (for 'corners total' is the total number of
       ! necessary communications. These special cases have their own
       ! send and receive buffers (send/recvbuffer_cc), their tags, etc.
       nsend_ct=0
       nrecv_ct=0
       recvsize_cc=0
       sendsize_cc=0
     end if

     do idims= idim^LIM
       select case (idims)
         {case (^D)
         nrecv=nrecv+nrecv_fc(^D)
         nsend=nsend+nsend_fc(^D)
         nxCo^D=1^D%nxCo^DD=ixGhi^DD/2-nghostcells;
         isize(^D)={nxCo^DD*}*(nwfluxin)
         recvsize=recvsize+nrecv_fc(^D)*isize(^D)
         sendsize=sendsize+nsend_fc(^D)*isize(^D)
         if(stagger_grid) then
           ! This does not consider the 'coarse corner' case
           nxCo^D=1^D%nxCo^DD=ixGhi^DD/2-nghostcells+1;
           isize_stg(^D)={nxCo^DD*}*(^ND-1)
           ! the whole size is used (cell centered and staggered)
           isize(^D)=isize(^D)+isize_stg(^D)      
           recvsize=recvsize+nrecv_fc(^D)*isize_stg(^D)
           sendsize=sendsize+nsend_fc(^D)*isize_stg(^D)
           ! Coarse corner case
           nrecv_ct=nrecv_ct+nrecv_cc(^D)
           nsend_ct=nsend_ct+nsend_cc(^D)
           recvsize_cc=recvsize_cc+nrecv_cc(^D)*isize_stg(^D)
           sendsize_cc=sendsize_cc+nsend_cc(^D)*isize_stg(^D)
         end if
         \}
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

     if (allocated(sendbuffer)) then
       if (sendsize /= size(sendbuffer)) then
         deallocate(sendbuffer)
         allocate(sendbuffer(sendsize))
       end if
     else
       allocate(sendbuffer(sendsize))
     end if

     if (allocated(fc_sendreq)) then
       if (nsend /= size(fc_sendreq)) then
         deallocate(fc_sendreq, fc_sendstat)
         allocate(fc_sendstat(MPI_STATUS_SIZE,nsend), fc_sendreq(nsend))
       end if
     else
       allocate(fc_sendstat(MPI_STATUS_SIZE,nsend), fc_sendreq(nsend))
     end if

     if(stagger_grid) then
       if (allocated(recvbuffer_cc)) then
         if (recvsize_cc /= size(recvbuffer_cc)) then
           deallocate(recvbuffer_cc)
           allocate(recvbuffer_cc(recvsize_cc))
         end if
       else
         allocate(recvbuffer_cc(recvsize_cc))
       end if

       if (allocated(cc_recvreq)) then
         if (nrecv_ct /= size(cc_recvreq)) then
           deallocate(cc_recvreq, cc_recvstat)
           allocate(cc_recvstat(MPI_STATUS_SIZE,nrecv_ct), cc_recvreq(nrecv_ct))
         end if
       else
         allocate(cc_recvstat(MPI_STATUS_SIZE,nrecv_ct), cc_recvreq(nrecv_ct))
       end if

       if (allocated(sendbuffer_cc)) then
         if (sendsize_cc /= size(sendbuffer_cc)) then
           deallocate(sendbuffer_cc)
           allocate(sendbuffer_cc(sendsize_cc))
         end if
       else
         allocate(sendbuffer_cc(sendsize_cc))
       end if

       if (allocated(cc_sendreq)) then
         if (nsend_ct /= size(cc_sendreq)) then
           deallocate(cc_sendreq, cc_sendstat)
           allocate(cc_sendstat(MPI_STATUS_SIZE,nsend_ct), cc_sendreq(nsend_ct))
         end if
       else
         allocate(cc_sendstat(MPI_STATUS_SIZE,nsend_ct), cc_sendreq(nsend_ct))
       end if
     end if

   end subroutine init_comm_fix_conserve

   subroutine recvflux(idim^LIM)
     use mod_global_parameters

     integer, intent(in) :: idim^LIM

     integer :: iigrid, igrid, idims, iside, i^D, nxCo^D
     integer :: ic^D, inc^D, ipe_neighbor
     integer :: pi^D,mi^D,ph^D,mh^D,idir

     if (nrecv>0) then
       fc_recvreq=MPI_REQUEST_NULL
       ibuf=1
       irecv=0

       do iigrid=1,igridstail; igrid=igrids(iigrid);
         do idims= idim^LIM
           do iside=1,2
             i^D=kr(^D,idims)*(2*iside-3);

             if (neighbor_pole(i^D,igrid)/=0) cycle

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

     if(stagger_grid) then
     ! receive corners
       if (nrecv_ct>0) then
         cc_recvreq=MPI_REQUEST_NULL
         ibuf_cc=1
         irecv=0

         do iigrid=1,igridstail; igrid=igrids(iigrid);
           do idims= idim^LIM
             do iside=1,2
               i^D=kr(^D,idims)*(2*iside-3);
               ! Check if there are special corners
               ! (Coarse block diagonal to a fine block)
               ! If there are, receive.
               ! Tags are calculated in the same way as for
               ! normal fluxes, but should not overlap because
               ! inc^D are different
               if (neighbor_type(i^D,igrid)==3) then
                 do idir=idims+1,ndim
                   pi^D=i^D+kr(idir,^D);
                   mi^D=i^D-kr(idir,^D);
                   ph^D=pi^D-kr(idims,^D)*(2*iside-3);
                   mh^D=mi^D-kr(idims,^D)*(2*iside-3);

                   if (neighbor_type(pi^D,igrid)==4.and.&
                       neighbor_type(ph^D,igrid)==3.and.&
                       neighbor_pole(pi^D,igrid)==0) then
                      ! Loop on children (several in 3D)
                    {do ic^DB=1+int((1-pi^DB)/2),2-int((1+pi^DB)/2)
                        inc^DB=2*pi^DB+ic^DB\}
                        ipe_neighbor=neighbor_child(2,inc^D,igrid)
                       if (mype/=ipe_neighbor) then
                         irecv=irecv+1
                         itag_cc=4**^ND*(igrid-1)+{inc^D*4**(^D-1)+}
                         call MPI_IRECV(recvbuffer_cc(ibuf_cc),isize_stg(idims),&
                                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                                        icomm,cc_recvreq(irecv),ierrmpi)
                         ibuf_cc=ibuf_cc+isize_stg(idims)
                       end if
                    {end do\}
                   end if

                   if (neighbor_type(mi^D,igrid)==4.and.&
                       neighbor_type(mh^D,igrid)==3.and.&
                       neighbor_pole(mi^D,igrid)==0) then
                      ! Loop on children (several in 3D)
                    {do ic^DB=1+int((1-mi^DB)/2),2-int((1+mi^DB)/2)
                        inc^DB=2*mi^DB+ic^DB\}
                       ipe_neighbor=neighbor_child(2,inc^D,igrid)
                       if (mype/=ipe_neighbor) then
                         irecv=irecv+1
                         itag_cc=4**^ND*(igrid-1)+{inc^D*4**(^D-1)+}
                         call MPI_IRECV(recvbuffer_cc(ibuf_cc),isize_stg(idims),&
                                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                                        icomm,cc_recvreq(irecv),ierrmpi)
                         ibuf_cc=ibuf_cc+isize_stg(idims)
                       end if
                    {end do\}
                   end if
                 end do
               end if
             end do
           end do
         end do
       end if
     end if ! end if stagger grid

   end subroutine recvflux

   subroutine allocateBflux
     use mod_global_parameters

     integer :: iigrid, igrid, iside, i^D, nx^D, nxCo^D
     integer :: idir,idim,pi^D, mi^D, ph^D, mh^D ! To detect corners

     nx^D=ixMhi^D-ixMlo^D+1;
     nxCo^D=nx^D/2;

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! For every grid,
       ! arrays for the fluxes are allocated for every face direction(^D)
       ! and every side (1=left, 2=right)
       {do iside=1,2
         i^DD=kr(^DD,^D)*(2*iside-3);

         if (neighbor_pole(i^DD,igrid)/=0) cycle

         select case (neighbor_type(i^DD,igrid))
         case(neighbor_fine)
           allocate(pflux(iside,^D,igrid)%flux(1^D%1:nx^DD,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,^D,igrid)%edge(1^D%0:nx^DD,1:ndim-1))
         case(neighbor_coarse)
           allocate(pflux(iside,^D,igrid)%flux(1^D%1:nxCo^DD,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,^D,igrid)%edge(1^D%0:nxCo^DD,1:ndim-1))
         case(neighbor_sibling)
           if(stagger_grid) then
             idim=^D
             do idir=idim+1,ndim
               pi^DD=i^DD+kr(idir,^DD);
               mi^DD=i^DD-kr(idir,^DD);
               ph^DD=pi^DD-kr(^D,^DD)*(2*iside-3);
               mh^DD=mi^DD-kr(^D,^DD)*(2*iside-3);
               if ((neighbor_type(pi^DD,igrid)==4&
                   .and.neighbor_type(ph^DD,igrid)==3)&
                   .or.(neighbor_type(mi^DD,igrid)==4&
                   .and.neighbor_type(mh^DD,igrid)==3)) then
                 allocate(pflux(iside,^D,igrid)%edge(1^D%0:nx^DD,1:ndim-1))
                 exit
               end if
             end do
           end if
         end select
       end do\}
     end do

   end subroutine allocateBflux

   subroutine deallocateBflux
     use mod_global_parameters

     integer :: iigrid, igrid, iside

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       {do iside=1,2
         if (associated(pflux(iside,^D,igrid)%flux)) then
           deallocate(pflux(iside,^D,igrid)%flux)
           nullify(pflux(iside,^D,igrid)%flux)
         end if
         if (associated(pflux(iside,^D,igrid)%edge)) then
           deallocate(pflux(iside,^D,igrid)%edge)
           nullify(pflux(iside,^D,igrid)%edge)
         end if
       end do\}
     end do

   end subroutine deallocateBflux

   subroutine fix_conserve(psb,idim^LIM,nw0,nwfluxin)
     use mod_global_parameters

     integer, intent(in) :: idim^LIM, nw0, nwfluxin
     type(state) :: psb(max_blocks)

     integer :: iigrid, igrid, idims, iside, iotherside, i^D, ic^D, inc^D, ix^L
     integer :: nxCo^D, iw, ix, ipe_neighbor, ineighbor, nbuf, ibufnext, nw1
     double precision :: CoFiratio

     nw1=nw0-1+nwfluxin
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

     ! for all grids: perform flux update at Coarse-Fine interfaces
     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims= idim^LIM
         select case (idims)
           {case (^D)
           do iside=1,2
             i^DD=kr(^DD,^D)*(2*iside-3);

             if (neighbor_pole(i^DD,igrid)/=0) cycle

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
               psb(igrid)%w(ix^D%ixM^T,nw0:nw1) &
                    = psb(igrid)%w(ix^D%ixM^T,nw0:nw1) &
                    -pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwfluxin)
             else
               do iw=nw0,nw1
                 psb(igrid)%w(ix^D%ixM^T,iw)=psb(igrid)%w(ix^D%ixM^T,iw)&
                      -pflux(iside,^D,igrid)%flux(1^D%:^DD&,iw-nw0+1) &
                      /ps(igrid)%dvolume(ix^D%ixM^T)
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
                   psb(igrid)%w(ix^S,nw0:nw1) &
                        = psb(igrid)%w(ix^S,nw0:nw1) &
                        + pflux(iotherside,^D,ineighbor)%flux(:^DD&,1:nwfluxin)&
                        * CoFiratio
                 else
                   do iw=nw0,nw1
                     psb(igrid)%w(ix^S,iw)=psb(igrid)%w(ix^S,iw) &
                          +pflux(iotherside,^D,ineighbor)%flux(:^DD&,iw-nw0+1) &
                          /ps(igrid)%dvolume(ix^S)
                   end do
                 end if
               else
                 if (slab) then
                   ibufnext=ibuf+isize(^D)
                   if(stagger_grid) ibufnext=ibufnext-isize_stg(^D)
                   psb(igrid)%w(ix^S,nw0:nw1) &
                        = psb(igrid)%w(ix^S,nw0:nw1)+CoFiratio &
                        *reshape(source=recvbuffer(ibuf:ibufnext-1), &
                        shape=shape(psb(igrid)%w(ix^S,nw0:nw1)))
                   ibuf=ibuf+isize(^D)
                 else
                   if(stagger_grid) then
                     nbuf=(isize(^D)-isize_stg(^D))/nwfluxin
                   else
                     nbuf=isize(^D)/nwfluxin
                   end if
                   do iw=nw0,nw1
                     ibufnext=ibuf+nbuf
                     psb(igrid)%w(ix^S,iw)=psb(igrid)%w(ix^S,iw) &
                          +reshape(source=recvbuffer(ibuf:ibufnext-1), &
                          shape=shape(psb(igrid)%w(ix^S,iw))) &
                          /ps(igrid)%dvolume(ix^S)
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

   subroutine storeflux(igrid,fC,idim^LIM,nwfluxin)
     use mod_global_parameters

     integer, intent(in)          :: igrid, idim^LIM, nwfluxin
     double precision, intent(in) :: fC(ixG^T,1:nwfluxin,1:ndim)

     integer :: idims, iside, i^D, ic^D, inc^D, ix^D, ixCo^D, nxCo^D, iw

     do idims = idim^LIM
       select case (idims)
         {case (^D)
         do iside=1,2
           i^DD=kr(^DD,^D)*(2*iside-3);

           if (neighbor_pole(i^DD,igrid)/=0) cycle

           select case (neighbor_type(i^DD,igrid))
           case (neighbor_fine)
             select case (iside)
             case (1)
               pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwfluxin) = &
                    -fC(nghostcells^D%ixM^T,1:nwfluxin,^D)
             case (2)
               pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwfluxin) = &
                    fC(ixMhi^D^D%ixM^T,1:nwfluxin,^D)
             end select
           case (neighbor_coarse)
             nxCo^D=1^D%nxCo^DD=ixGhi^DD/2-nghostcells;
             select case (iside)
             case (1)
               do iw=1,nwfluxin
                 {do ixCo^DDB=1,nxCo^DDB\}
                 ix^D=nghostcells^D%ix^DD=ixMlo^DD+2*(ixCo^DD-1);
                 pflux(iside,^D,igrid)%flux(ixCo^DD,iw) &
                      = {^NOONEDsum}(fC(ix^D^D%ix^DD:ix^DD+1,iw,^D))
                 {end do\}
               end do
             case (2)
               do iw=1,nwfluxin
                 {do ixCo^DDB=1,nxCo^DDB\}
                 ix^D=ixMhi^D^D%ix^DD=ixMlo^DD+2*(ixCo^DD-1);
                 pflux(iside,^D,igrid)%flux(ixCo^DD,iw) &
                      =-{^NOONEDsum}(fC(ix^D^D%ix^DD:ix^DD+1,iw,^D))
                 {end do\}
               end do
             end select
           end select
         end do\}
       end select
     end do

   end subroutine storeflux

   subroutine sendflux(idim^LIM)
     use mod_global_parameters

     integer, intent(in) :: idim^LIM

     integer :: idims, iside, i^D, ic^D, inc^D, ix^D, ixCo^D, nxCo^D, iw
     integer :: ineighbor, ipe_neighbor, igrid, iigrid, ibuf_send_next
     integer :: idir, ibuf_cc_send_next, pi^D, ph^D, mi^D, mh^D

     fc_sendreq = MPI_REQUEST_NULL
     isend      = 0
     ibuf_send  = 1
     if(stagger_grid) then
       cc_sendreq=MPI_REQUEST_NULL
       isend_cc=0
       ibuf_cc_send=1
     end if

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims = idim^LIM
         select case (idims)
        {case (^D)
           do iside=1,2
             i^DD=kr(^DD,^D)*(2*iside-3);

             if (neighbor_pole(i^DD,igrid)/=0) cycle

             if (neighbor_type(i^DD,igrid)==neighbor_coarse) then
               ! send flux to coarser neighbor
               ineighbor=neighbor(1,i^DD,igrid)
               ipe_neighbor=neighbor(2,i^DD,igrid)
               if (ipe_neighbor/=mype) then
                 ic^DD=1+modulo(node(pig^DD_,igrid)-1,2);
                 inc^D=-2*i^D+ic^D^D%inc^DD=ic^DD;
                 itag=4**^ND*(ineighbor-1)+{inc^DD*4**(^DD-1)+}
                 isend=isend+1

                 ibuf_send_next=ibuf_send+isize(^D)
                 if(stagger_grid) then
                   sendbuffer(ibuf_send:ibuf_send_next-isize_stg(^D)-1)=&
                   reshape(pflux(iside,^D,igrid)%flux,(/isize(^D)-isize_stg(^D)/))

                   sendbuffer(ibuf_send_next-isize_stg(^D):ibuf_send_next-1)=&
                   reshape(pflux(iside,^D,igrid)%edge,(/isize_stg(^D)/))
                 else
                   sendbuffer(ibuf_send:ibuf_send_next-1)=&
                   reshape(pflux(iside,^D,igrid)%flux,(/isize(^D)/))
                 end if
                 call MPI_ISEND(sendbuffer(ibuf_send),isize(^D), &
                      MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                      icomm,fc_sendreq(isend),ierrmpi)
                 ibuf_send=ibuf_send_next
               end if

               if(stagger_grid) then
                 ! If we are in a fine block surrounded by coarse blocks
                 do idir=idims+1,ndim
                   pi^DD=i^DD+kr(idir,^DD);
                   mi^DD=i^DD-kr(idir,^DD);
                   ph^DD=pi^DD-kr(idims,^DD)*(2*iside-3);
                   mh^DD=mi^DD-kr(idims,^DD)*(2*iside-3);

                   if (neighbor_type(pi^DD,igrid)==2.and.&
                       neighbor_type(ph^DD,igrid)==2.and.&
                       mype/=neighbor(2,pi^DD,igrid).and.&
                       neighbor_pole(pi^DD,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,pi^DD,igrid)
                     ipe_neighbor=neighbor(2,pi^DD,igrid)
                     ic^DD=1+modulo(node(pig^DD_,igrid)-1,2);
                     inc^DD=-2*pi^DD+ic^DD;
                     itag_cc=4**^ND*(ineighbor-1)+{inc^DD*4**(^DD-1)+}
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(^D)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-1)=&
                     reshape(pflux(iside,^D,igrid)%edge,shape=(/isize_stg(^D)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(^D),&
                                    MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                                    icomm,cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
       
                   if (neighbor_type(mi^DD,igrid)==2.and.&
                       neighbor_type(mh^DD,igrid)==2.and.&
                       mype/=neighbor(2,mi^DD,igrid).and.&
                       neighbor_pole(mi^DD,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,mi^DD,igrid)
                     ipe_neighbor=neighbor(2,mi^DD,igrid)
                     ic^DD=1+modulo(node(pig^DD_,igrid)-1,2);
                     inc^DD=-2*pi^DD+ic^DD;
                     inc^DD=-2*mi^DD+ic^DD;
                     itag_cc=4**^ND*(ineighbor-1)+{inc^DD*4**(^DD-1)+}
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(^D)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-1)=&
                     reshape(pflux(iside,^D,igrid)%edge,shape=(/isize_stg(^D)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(^D),&
                                    MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                                    icomm,cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
                 end do
               end if ! end if stagger grid

             end if
           end do\}
         end select
       end do
     end do
   end subroutine sendflux

end module mod_fix_conserve
