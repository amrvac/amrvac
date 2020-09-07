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
  integer                              :: itag, itag_cc, isend, isend_cc, irecv, irecv_cc

  public :: init_comm_fix_conserve
  public :: allocateBflux
  public :: deallocateBflux
  public :: sendflux
  public :: recvflux
  public :: store_flux, store_flux_var
  public :: store_edge
  public :: fix_conserve, fix_conserve_vars
  public :: fix_edges

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

     if (allocated(fc_sendreq)) then
       if (nsend /= size(fc_sendreq)) then
         deallocate(fc_sendreq, fc_sendstat)
         allocate(fc_sendstat(MPI_STATUS_SIZE,nsend), fc_sendreq(nsend))
       end if
     else
       allocate(fc_sendstat(MPI_STATUS_SIZE,nsend), fc_sendreq(nsend))
     end if

     if(stagger_grid) then

       if (allocated(sendbuffer)) then
         if (sendsize /= size(sendbuffer)) then
           deallocate(sendbuffer)
           allocate(sendbuffer(sendsize))
         end if
       else
         allocate(sendbuffer(sendsize))
       end if

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
         irecv_cc=0

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
                         irecv_cc=irecv_cc+1
                         itag_cc=4**^ND*(igrid-1)+{inc^D*4**(^D-1)+}
                         call MPI_IRECV(recvbuffer_cc(ibuf_cc),isize_stg(idims),&
                                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                                        icomm,cc_recvreq(irecv_cc),ierrmpi)
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
                         irecv_cc=irecv_cc+1
                         itag_cc=4**^ND*(igrid-1)+{inc^D*4**(^D-1)+}
                         call MPI_IRECV(recvbuffer_cc(ibuf_cc),isize_stg(idims),&
                                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,&
                                        icomm,cc_recvreq(irecv_cc),ierrmpi)
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

   subroutine sendflux(idim^LIM)
     use mod_global_parameters

     integer, intent(in) :: idim^LIM

     integer :: idims, iside, i^D, ic^D, inc^D, ix^D, ixCo^D, nxCo^D, iw
     integer :: ineighbor, ipe_neighbor, igrid, iigrid, ibuf_send_next
     integer :: idir, ibuf_cc_send_next, pi^D, ph^D, mi^D, mh^D

     fc_sendreq = MPI_REQUEST_NULL
     isend      = 0
     if(stagger_grid) then
       ibuf_send  = 1
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

                 if(stagger_grid) then
                   ibuf_send_next=ibuf_send+isize(^D)
                   sendbuffer(ibuf_send:ibuf_send_next-isize_stg(^D)-1)=&
                   reshape(pflux(iside,^D,igrid)%flux,(/isize(^D)-isize_stg(^D)/))

                   sendbuffer(ibuf_send_next-isize_stg(^D):ibuf_send_next-1)=&
                   reshape(pflux(iside,^D,igrid)%edge,(/isize_stg(^D)/))
                   call MPI_ISEND(sendbuffer(ibuf_send),isize(^D), &
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                        icomm,fc_sendreq(isend),ierrmpi)
                   ibuf_send=ibuf_send_next
                 else
                   call MPI_ISEND(pflux(iside,^D,igrid)%flux,isize(^D), &
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                        icomm,fc_sendreq(isend),ierrmpi)
                 end if
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
             !do idir=min(idim+1,ndim),ndim
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
     call fix_conserve1(psb,idim^LIM,nw0,1,nwfluxin)

   end subroutine fix_conserve

  !> modified fix_conserve in order to pass param for nwfluxstart = start index in pflux var
  !>compatible with store_flux1
   subroutine fix_conserve1(psb,idim^LIM,nw0,nwfluxstart,nwfluxin)
     use mod_global_parameters

     integer, intent(in) :: idim^LIM, nw0, nwfluxstart, nwfluxin
     type(state) :: psb(max_blocks)

     integer :: iigrid, igrid, idims, iside, iotherside, i^D, ic^D, inc^D, ix^L
     integer :: nxCo^D, iw, ix, ipe_neighbor, ineighbor, nbuf, ibufnext, nw1
     double precision :: CoFiratio

     nw1=nw0-1+nwfluxin
     if (slab_uniform) then
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
             if (slab_uniform) then
               psb(igrid)%w(ix^D%ixM^T,nw0:nw1) &
                    = psb(igrid)%w(ix^D%ixM^T,nw0:nw1) &
                    -pflux(iside,^D,igrid)%flux(1^D%:^DD&,nwfluxstart:nwfluxin+nwfluxstart-1)
             else
               do iw=nw0,nw1
                 psb(igrid)%w(ix^D%ixM^T,iw)=psb(igrid)%w(ix^D%ixM^T,iw)&
                      -pflux(iside,^D,igrid)%flux(1^D%:^DD&,iw-nw0+nwfluxstart) &
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
                 if (slab_uniform) then
                   psb(igrid)%w(ix^S,nw0:nw1) &
                        = psb(igrid)%w(ix^S,nw0:nw1) &
                        + pflux(iotherside,^D,ineighbor)%flux(:^DD&,nwfluxstart:nwfluxin+nwfluxstart-1)&
                        * CoFiratio
                 else
                   do iw=nw0,nw1
                     psb(igrid)%w(ix^S,iw)=psb(igrid)%w(ix^S,iw) &
                          +pflux(iotherside,^D,ineighbor)%flux(:^DD&,iw-nw0+nwfluxstart) &
                          /ps(igrid)%dvolume(ix^S)
                   end do
                 end if
               else
                 if (slab_uniform) then
                   ibufnext=ibuf+isize(^D)
                   if(stagger_grid) ibufnext=ibufnext-isize_stg(^D)
                   psb(igrid)%w(ix^S,nw0:nw1) &
                        = psb(igrid)%w(ix^S,nw0:nw1)+CoFiratio &
                        *reshape(source=recvbuffer(ibuf:ibufnext-1), &
                        shape=shape(psb(igrid)%w(ix^S,nw0:nw1)))
                   ibuf=ibuf+isize(^D)
                 else
                   ibufnext=ibuf+isize(^D)
                   if(stagger_grid) then
                     nbuf=(isize(^D)-isize_stg(^D))/nwfluxin
                   else
                     nbuf=isize(^D)/nwfluxin
                   end if
                   do iw=nw0,nw1
                     psb(igrid)%w(ix^S,iw)=psb(igrid)%w(ix^S,iw) &
                          +reshape(source=recvbuffer(ibuf:ibufnext-1), &
                          shape=shape(psb(igrid)%w(ix^S,iw))) &
                          /ps(igrid)%dvolume(ix^S)
                     ibuf=ibuf+nbuf
                   end do
                   ibuf=ibufnext
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

   end subroutine fix_conserve1

   subroutine store_flux(igrid,fC,idim^LIM,nwfluxin)
     use mod_global_parameters

     integer, intent(in)          :: igrid, idim^LIM, nwfluxin
     double precision, intent(in) :: fC(ixG^T,1:nwfluxin,1:ndim)
     call store_flux1(igrid,fC,idim^LIM,1,nwfluxin) 
   end subroutine store_flux


   !> old store_flux modified in order to pass parameter nwfstart1 = start index in pflux
   subroutine store_flux1(igrid,fC,idim^LIM,nwfstart1,nwfluxin)
     use mod_global_parameters

     integer, intent(in)          :: igrid, idim^LIM, nwfstart1, nwfluxin
     integer                      :: nwfend1
     double precision, intent(in) :: fC(ixG^T,1:nwfluxin,1:ndim)

     integer :: idims, iside, i^D, ic^D, inc^D, ix^D, ixCo^D, nxCo^D, iw
     nwfend1 = nwfstart1-1 + nwfluxin
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
               pflux(iside,^D,igrid)%flux(1^D%:^DD&,nwfstart1:nwfend1) = &
                    -fC(nghostcells^D%ixM^T,1:nwfluxin,^D)
             case (2)
               pflux(iside,^D,igrid)%flux(1^D%:^DD&,nwfstart1:nwfend1) = &
                    fC(ixMhi^D^D%ixM^T,1:nwfluxin,^D)
             end select
           case (neighbor_coarse)
             nxCo^D=1^D%nxCo^DD=ixGhi^DD/2-nghostcells;
             select case (iside)
             case (1)
               do iw=1,nwfluxin
                {do ixCo^DDB=1,nxCo^DDB\}
                   ix^D=nghostcells^D%ix^DD=ixMlo^DD+2*(ixCo^DD-1);
                   pflux(iside,^D,igrid)%flux(ixCo^DD,nwfstart1-1+iw) &
                        = {^NOONEDsum}(fC(ix^D^D%ix^DD:ix^DD+1,iw,^D))
                {end do\}
               end do
             case (2)
               do iw=1,nwfluxin
                {do ixCo^DDB=1,nxCo^DDB\}
                   ix^D=ixMhi^D^D%ix^DD=ixMlo^DD+2*(ixCo^DD-1);
                   pflux(iside,^D,igrid)%flux(ixCo^DD,nwfstart1-1+iw) &
                        =-{^NOONEDsum}(fC(ix^D^D%ix^DD:ix^DD+1,iw,^D))
                {end do\}
               end do
             end select
           end select
         end do\}
       end select
     end do

   end subroutine store_flux1

   subroutine store_edge(igrid,ixI^L,fE,idim^LIM)
     use mod_global_parameters
     
     integer, intent(in)          :: igrid, ixI^L, idim^LIM
     double precision, intent(in) :: fE(ixI^S,7-2*ndim:3)
     
     integer :: idims, idir, iside, i^D
     integer :: pi^D, mi^D, ph^D, mh^D ! To detect corners
     integer :: ixMc^L

     do idims = idim^LIM  !loop over face directions
       !! Loop over block faces
       do iside=1,2 
         i^D=kr(^D,idims)*(2*iside-3);
         if (neighbor_pole(i^D,igrid)/=0) cycle
         select case (neighbor_type(i^D,igrid))
         case (neighbor_fine)
           ! The neighbour is finer
           ! Face direction, side (left or right), restrict ==ired?, fE
           call flux_to_edge(igrid,ixI^L,idims,iside,.false.,fE)
         case(neighbor_coarse)
           ! The neighbour is coarser
           call flux_to_edge(igrid,ixI^L,idims,iside,.true.,fE)
         case(neighbor_sibling)
           ! If the neighbour is at the same level,
           ! check if there are corners
           ! If there is any corner, store the fluxes from that side
           do idir=idims+1,ndim
             pi^D=i^D+kr(idir,^D);
             mi^D=i^D-kr(idir,^D);
             ph^D=pi^D-kr(idims,^D)*(2*iside-3);
             mh^D=mi^D-kr(idims,^D)*(2*iside-3);
             if (neighbor_type(pi^D,igrid)==4&
               .and.neighbor_type(ph^D,igrid)==3) then
               call flux_to_edge(igrid,ixI^L,idims,iside,.false.,fE)
             end if
             if (neighbor_type(mi^D,igrid)==4&
               .and.neighbor_type(mh^D,igrid)==3) then
               call flux_to_edge(igrid,ixI^L,idims,iside,.false.,fE)
             end if
           end do
         end select
       end do
     end do

   end subroutine store_edge

   subroutine flux_to_edge(igrid,ixI^L,idims,iside,restrict,fE)
     use mod_global_parameters

     integer                      :: igrid,ixI^L,idims,iside
     logical                      :: restrict
     double precision, intent(in) :: fE(ixI^S,7-2*ndim:3)

     integer                      :: idir1,idir2
     integer                      :: ixE^L,ixF^L{^IFTHREED, jxF^L,}, nx^D,nxCo^D

     nx^D=ixMhi^D-ixMlo^D+1;
     nxCo^D=nx^D/2;
     ! ixE are the indices on the 'edge' array.
     ! ixF are the indices on the 'fE' array
     ! jxF are indices advanced to perform the flux restriction (sum) in 3D
     ! A line integral of the electric field on the coarse side
     ! lies over two edges on the fine side. So, in 3D we restrict by summing
     ! over two cells on the fine side.

     do idir1=1,ndim-1
      {^IFTHREED ! 3D: rotate indices among 1 and 2 to save space 
       idir2=mod(idir1+idims-1,3)+1}
      {^IFTWOD ! Assign the only field component (3) to the only index (1)
       idir2=3}

       if (restrict) then
         ! Set up indices for restriction
         ixFmin^D=ixMlo^D-1+kr(^D,idir2);
         ixFmax^D=ixMhi^D-kr(^D,idir2);
         {^IFTHREED
         jxF^L=ixF^L+kr(^D,idir2);}

         ixEmin^D=0+kr(^D,idir2);
         ixEmax^D=nxCo^D;
         select case(idims)
        {case(^D)
           ixEmin^D=1;ixEmax^D=1;
           select case(iside)
           case(1)
             ixFmax^D=ixFmin^D
             {^IFTHREEDjxFmax^D=ixFmin^D}
           case(2)
             ixFmin^D=ixFmax^D
             {^IFTHREEDjxFmin^D=ixFmax^D}
           end select
        \}
         end select

       pflux(iside,idims,igrid)%edge(ixE^S,idir1)=&
         fE(ixFmin^D:ixFmax^D:2,idir2){^IFTHREED +&
         fE(jxFmin^D:jxFmax^D:2,idir2)};

       else
         ! Set up indices for copying 
         ixFmin^D=ixMlo^D-1+kr(^D,idir2);
         ixFmax^D=ixMhi^D;
         ixEmin^D=0+kr(^D,idir2);
         ixEmax^D=nx^D;

         select case(idims)
        {case(^D)
           ixEmin^D=1;ixEmax^D=1;
           select case(iside)
           case(1)
             ixFmax^D=ixFmin^D
           case(2)
             ixFmin^D=ixFmax^D
           end select
        \}
         end select

         pflux(iside,idims,igrid)%edge(ixE^S,idir1)=fE(ixF^S,idir2)

       end if

     end do

   end subroutine flux_to_edge

   subroutine fix_edges(psuse,idim^LIM)
     use mod_global_parameters

     type(state) :: psuse(max_blocks)
     integer, intent(in) :: idim^LIM

     integer :: iigrid, igrid, idims, iside, iotherside, i^D, ic^D, inc^D, ixMc^L
     integer :: nbuf, ibufnext
     integer :: ibufnext_cc
     integer :: pi^D, mi^D, ph^D, mh^D ! To detect corners
     integer :: ixE^L(1:3), ixtE^L, ixF^L(1:ndim), ixfE^L(1:3)
     integer :: nx^D, idir, ix, ipe_neighbor, ineighbor
     logical :: pcorner(1:ndim),mcorner(1:ndim)

     if (nrecv_ct>0) then
        call MPI_WAITALL(nrecv_ct,cc_recvreq,cc_recvstat,ierrmpi)
     end if

     ! Initialise buffer counter again
     ibuf=1
     ibuf_cc=1
     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims= idim^LIM
         do iside=1,2
           i^D=kr(^D,idims)*(2*iside-3);
           if (neighbor_pole(i^D,igrid)/=0) cycle
           select case(neighbor_type(i^D,igrid))
           case(neighbor_fine)
             ! The first neighbour is finer
             if (.not.neighbor_active(i^D,igrid).or.&
                 .not.neighbor_active(0^D&,igrid) ) then
               {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
                  inc^DB=2*i^DB+ic^DB\}
                  ipe_neighbor=neighbor_child(2,inc^D,igrid)
                  !! When the neighbour is in a different process
                  if (ipe_neighbor/=mype) then
                     ibufnext=ibuf+isize(idims)
                     ibuf=ibufnext
                     end if
               {end do\}
                cycle
             end if

             ! Check if there are corners
             pcorner=.false.
             mcorner=.false.
             do idir=1,ndim
               pi^D=i^D+kr(idir,^D);
               mi^D=i^D-kr(idir,^D);
               ph^D=pi^D-kr(idims,^D)*(2*iside-3);
               mh^D=mi^D-kr(idims,^D)*(2*iside-3);
               if (neighbor_type(ph^D,igrid)==neighbor_fine) pcorner(idir)=.true.
               if (neighbor_type(mh^D,igrid)==neighbor_fine) mcorner(idir)=.true.
             end do
             ! Calculate indices range
             call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.false.,.false.,0^D&,pcorner,mcorner)
             ! Remove coarse part of circulation
             call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iside,idims,igrid)%edge,idims,iside,.false.,psuse(igrid))
             ! Add fine part of the circulation
            {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
               inc^DB=2*i^DB+ic^DB\}
               ineighbor=neighbor_child(1,inc^D,igrid)
               ipe_neighbor=neighbor_child(2,inc^D,igrid)
               iotherside=3-iside
               nx^D=(ixMhi^D-ixMlo^D+1)/2;
               call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.true.,.false.,inc^D,pcorner,mcorner)
               if (ipe_neighbor==mype) then
                 call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iotherside,idims,ineighbor)%edge,idims,iside,.true.,psuse(igrid))
               else
                 ibufnext=ibuf+isize(idims)
                 call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,&
                   reshape(source=recvbuffer(ibufnext-isize_stg(idims):ibufnext-1),&
                   shape=(/ ixtEmax^D-ixtEmin^D+1 ,^ND-1 /)),&
                   idims,iside,.true.,psuse(igrid))
                 ibuf=ibufnext
               end if
            {end do\}

           case(neighbor_sibling)
             ! The first neighbour is at the same level
             ! Check if there are corners
             do idir=idims+1,ndim
               pcorner=.false.
               mcorner=.false.
               pi^D=i^D+kr(idir,^D);
               mi^D=i^D-kr(idir,^D);
               ph^D=pi^D-kr(idims,^D)*(2*iside-3);
               mh^D=mi^D-kr(idims,^D)*(2*iside-3);
               if (neighbor_type(pi^D,igrid)==neighbor_fine&
                 .and.neighbor_type(ph^D,igrid)==neighbor_sibling&
                 .and.neighbor_pole(pi^D,igrid)==0) then
                 pcorner(idir)=.true.
                 ! Remove coarse part
                 ! Set indices
                 call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.false.,.true.,0^D&,pcorner,mcorner)
                 ! Remove
                 call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iside,idims,igrid)%edge,idims,iside,.false.,psuse(igrid))
                 ! Add fine part
                 ! Find relative position of finer grid
      {^IFTHREED do ix=1,2}
                 inc^D=kr(idims,^D)*3*(iside-1)+3*kr(idir,^D){^IFTHREED+kr(6-idir-idims,^D)*ix};
                 ineighbor=neighbor_child(1,inc^D,igrid)
                 ipe_neighbor=neighbor_child(2,inc^D,igrid)
                 iotherside=3-iside
                 ! Set indices
                 call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.true.,.true.,inc^D,pcorner,mcorner)
                 ! add
                 if (ipe_neighbor==mype) then
                   call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iotherside,idims,ineighbor)%edge,idims,iside,.true.,psuse(igrid))
                 else
                   ibufnext_cc=ibuf_cc+isize_stg(idims)
                   call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,&
                     reshape(source=recvbuffer_cc(ibuf_cc:ibufnext_cc-1),&
                     shape=(/ ixtEmax^D-ixtEmin^D+1 ,^ND-1 /)),&
                     idims,iside,.true.,psuse(igrid))
                   ibuf_cc=ibufnext_cc
                 end if
      {^IFTHREED end do}
               ! Set CoCorner to false again for next step
                 pcorner(idir)=.false.
               end if

               if (neighbor_type(mi^D,igrid)==neighbor_fine&
                 .and.neighbor_type(mh^D,igrid)==neighbor_sibling&
                 .and.neighbor_pole(mi^D,igrid)==0) then
                   mcorner(idir)=.true.
                   ! Remove coarse part
                   ! Set indices
                   call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.false.,.true.,0^D&,pcorner,mcorner)
                   ! Remove
                   call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iside,idims,igrid)%edge,idims,iside,.false.,psuse(igrid))
                   ! Add fine part
                   ! Find relative position of finer grid
        {^IFTHREED do ix=1,2}
                   inc^D=kr(idims,^D)*3*(iside-1){^IFTHREED+kr(6-idir-idims,^D)*ix};
                   ineighbor=neighbor_child(1,inc^D,igrid)
                   ipe_neighbor=neighbor_child(2,inc^D,igrid)
                   iotherside=3-iside
                   ! Set indices
                   call set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,.true.,.true.,inc^D,pcorner,mcorner)
                   ! add
                   if (ipe_neighbor==mype) then
                     call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,pflux(iotherside,idims,ineighbor)%edge,idims,iside,.true.,psuse(igrid))
                   else
                     ibufnext_cc=ibuf_cc+isize_stg(idims)
                     call add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,&
                       reshape(source=recvbuffer_cc(ibuf_cc:ibufnext_cc-1),&
                       shape=(/ ixtEmax^D-ixtEmin^D+1 ,^ND-1 /)),&
                       idims,iside,.true.,psuse(igrid))
                     ibuf_cc=ibufnext_cc
                   end if
        {^IFTHREED end do}
                 ! Set CoCorner to false again for next step
                  mcorner(idir)=.false.
               end if
             end do
           end select
         end do
       end do
     end do

     if (nsend_ct>0) call MPI_WAITALL(nsend_ct,cc_sendreq,cc_sendstat,ierrmpi)

   end subroutine fix_edges

   !> This routine sets the indexes for the correction
   !> of the circulation according to several different
   !> cases, as grids located in different cpus,
   !> presence of corners, and different relative locations
   !> of the fine grid respect to the coarse one
   subroutine set_ix_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,igrid,idims,iside,add,CoCorner,inc^D,pcorner,mcorner)
     use mod_global_parameters
     
     integer,intent(in)    :: igrid,idims,iside,inc^D
     logical,intent(in)    :: add,CoCorner
     logical,intent(inout) :: pcorner(1:ndim),mcorner(1:ndim)
     integer,intent(out)   :: ixF^L(1:ndim),ixtE^L,ixE^L(1:3),ixfE^L(1:3) ! Indices for faces and edges
     integer               :: icor^D,idim1,idir,nx^D,middle^D
     integer               :: ixtfE^L

     ! ixF -> Indices for the _F_aces, and
     ! depends on the field component
     ! ixtE -> are the _t_otal range of the 'edge' array
     ! ixE -> are the ranges of the edge array,
     ! depending on the component
     ! ixfE -> are the ranges of the fE array (3D),
     ! and also depend on the component

     ! ... General ...
     ! Assign indices for the size of the E field array

     ixtfEmin^D=ixMlo^D-1;
     ixtfEmax^D=ixMhi^D;

     if(add) then
       nx^D=(ixMhi^D-ixMlo^D+1)/2;
     else
       nx^D=ixMhi^D-ixMlo^D+1;
     end if

     do idim1=1,ndim
       ixtEmin^D=0;
       ixtEmax^D=nx^D;
       select case(idims)
       {case(^D)
         ixtEmin^D=1;ixtEmax^D=1;
         if (iside==1) ixtfEmax^D=ixtfEmin^D;
         if (iside==2) ixtfEmin^D=ixtfEmax^D;
       \}
       end select
     end do

     ! Assign indices, considering only the face
     ! (idims and iside)
     do idim1=1,ndim
       ixFmin^D(idim1)=ixMlo^D-kr(idim1,^D);
       ixFmax^D(idim1)=ixMhi^D;
       select case(idims)
       {case(^D)
          select case(iside)
          case(1)
          ixFmax^D(idim1)=ixFmin^D(idim1)
          case(2)
          ixFmin^D(idim1)=ixFmax^D(idim1)
          end select
       \}
       end select
     end do
     ! ... Relative position ...
     ! Restrict range using relative position
     if(add) then
       middle^D=(ixMhi^D+ixMlo^D)/2;
       {
       if(inc^D==1) then
         ixFmax^D(:)=middle^D
         ixtfEmax^D=middle^D
       end if
       if(inc^D==2) then
         ixFmin^D(:)=middle^D+1
         ixtfEmin^D=middle^D
       end if
       \}
     end if
     ! ... Adjust ranges of edges according to direction ...
     do idim1=1,3
       ixfEmax^D(idim1)=ixtfEmax^D;
       ixEmax^D(idim1)=ixtEmax^D;
       ixfEmin^D(idim1)=ixtfEmin^D+kr(idim1,^D);
       ixEmin^D(idim1)=ixtEmin^D+kr(idim1,^D);
     end do
     ! ... Corners ...
     ! 'Coarse' corners
     if (CoCorner) then
       do idim1=idims+1,ndim
         if (pcorner(idim1)) then
           do idir=1,3!Index arrays have size ndim
             if (idir==6-idim1-idims) then
              !!! Something here has to change
              !!! Array ixfE must have size 3, while
              !!! ixE must have size ndim
              {if (^D==idim1) then
                 ixfEmin^D(idir)=ixfEmax^D(idir)
                 if (add) then
                   ixEmax^D(idir) =ixEmin^D(idir)
                 else
                   ixEmin^D(idir) =ixEmax^D(idir)
                 end if
               end if\}
             else
               ixEmin^D(idir)=1;
               ixEmax^D(idir)=0;
               ixfEmin^D(idir)=1;
               ixfEmax^D(idir)=0;
             end if
           end do
         end if
         if (mcorner(idim1)) then
           do idir=1,3
             if (idir==6-idim1-idims) then
              {if (^D==idim1) then
                 ixfEmax^D(idir)=ixfEmin^D(idir)
                 if (add) then
                   ixEmin^D(idir) =ixEmax^D(idir)
                 else
                   ixEmax^D(idir) =ixEmin^D(idir)
                 end if
               end if\}
             else
               ixEmin^D(idir)=1;
               ixEmax^D(idir)=0;
               ixfEmin^D(idir)=1;
               ixfEmax^D(idir)=0;
             end if
           end do
         end if
       end do
     else
     ! Other kinds of corners
     ! Crop ranges to account for corners
     ! When the fine fluxes are added, we consider 
     ! whether they come from the same cpu or from
     ! a different one, in order to minimise the 
     ! amount of communication
     ! Case for different processors still not implemented!!!
      {if((idims.gt.^D).and.pcorner(^D)) then
         if((.not.add).or.(inc^D==2)) then
           !ixFmax^DD(:)=ixFmax^DD(:)-kr(^D,^DD);
           do idir=1,3
             if ((idir==idims).or.(idir==^D)) cycle
               ixfEmax^D(idir)=ixfEmax^D(idir)-1
               ixEmax^D(idir)=ixEmax^D(idir)-1
           end do
         end if
       end if\}
      {if((idims>^D).and.mcorner(^D)) then
         if((.not.add).or.(inc^D==1)) then
           !ixFmin^DD(:)=ixFmin^DD(:)+kr(^D,^DD);
           do idir=1,3
             if ((idir==idims).or.(idir==^D)) cycle
               ixfEmin^D(idir)=ixfEmin^D(idir)+1
               ixEmin^D(idir)=ixEmin^D(idir)+1
           end do
         end if
       end if\}
     end if

   end subroutine set_ix_circ

   subroutine add_sub_circ(ixF^L,ixtE^L,ixE^L,ixfE^L,edge,idims,iside,add,s)
     use mod_global_parameters

     type(state)        :: s
     integer,intent(in) :: idims,iside
     integer            :: ixF^L(1:ndim),ixtE^L,ixE^L(1:3),ixfE^L(1:3)
     double precision   :: edge(ixtE^S,1:ndim-1)
     logical,intent(in) :: add

     integer            :: idim1,idim2,idir,middle^D
     integer            :: ixfEC^L,ixEC^L
     double precision   :: fE(ixG^T,7-2*ndim:3) !!!!!!!!
     double precision   :: circ(ixG^T,1:ndim) !!!!!!!!
     integer            :: ix^L,hx^L,ixC^L,hxC^L ! Indices for edges

     ! ixF -> Indices for the faces, depends on the field component
     ! ixE -> Total range for the edges
     ! ixfE -> Edges in fE (3D) array
     ! ix,hx,ixC,hxC -> Auxiliary indices
     ! Assign quantities stored ad edges to make it as similar as 
     ! possible to the routine updatefaces.
     fE(:^D&,:)=zero
     do idim1=1,ndim-1
       {^IFTHREED ! 3D: rotate indices (see routine flux_to_edge)
       idir=mod(idim1+idims-1,3)+1}
       {^IFTWOD   ! 2D: move E back to directon 3
       idir=3}
       ixfEC^L=ixfE^L(idir);
       ixEC^L=ixE^L(idir);
       fE(ixfEC^S,idir)=edge(ixEC^S,idim1)
     end do

     ! Calculate part of circulation needed
     circ=zero
     do idim1=1,ndim
        do idim2=1,ndim
           do idir=7-2*ndim,3
             if (lvc(idim1,idim2,idir)==0) cycle
             ! Assemble indices
             ixC^L=ixF^L(idim1);
             hxC^L=ixC^L-kr(idim2,^D);
             if(idim1==idims) then
               circ(ixC^S,idim1)=circ(ixC^S,idim1)+lvc(idim1,idim2,idir)&
                                 *(fE(ixC^S,idir)-fE(hxC^S,idir))
             else
               select case(iside)
               case(2)
                 circ(ixC^S,idim1)=circ(ixC^S,idim1)+lvc(idim1,idim2,idir)*fE(ixC^S,idir)
               case(1)
                 circ(ixC^S,idim1)=circ(ixC^S,idim1)-lvc(idim1,idim2,idir)*fE(hxC^S,idir)
               end select
             end if
           end do
        end do
     end do

     ! Divide circulation by surface and add
     do idim1=1,ndim
        ixC^L=ixF^L(idim1);
        where(s%surfaceC(ixC^S,idim1)>1.0d-9*s%dvolume(ixC^S))
          circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
        elsewhere
          circ(ixC^S,idim1)=zero
        end where
        ! Add/subtract to field at face
        if (add) then
          s%ws(ixC^S,idim1)=s%ws(ixC^S,idim1)-circ(ixC^S,idim1)
        else
          s%ws(ixC^S,idim1)=s%ws(ixC^S,idim1)+circ(ixC^S,idim1)
        end if
     end do

   end subroutine add_sub_circ


  !!this conserves fluxed put by store_flux_var
  !used for the sts methods where not all the fluxes are stored and fixed.
  subroutine fix_conserve_vars(tmpPs, indexChangeStart, indexChangeN, indexChangeFixC)
    use mod_global_parameters
    type(state), target               :: tmpPs(max_blocks)
    integer, intent(in), dimension(:) :: indexChangeStart, indexChangeN
    logical, intent(in), dimension(:) :: indexChangeFixC
    integer :: i,storeIndex,j

      call recvflux(1,ndim)
      call sendflux(1,ndim)
      storeIndex = 1
      !!This is consistent with the subroutine in mod_mhs_phys which gets the indices where
      !!to store the fluxes in mod_fix_conserve, set_sts_sources_ambipolar 
      do i = 1,size(indexChangeStart)
        if (indexChangeFixC(i)) then 
          !this has to be done one by one, as they are not stored at contiguous locations
          do j = 0,indexChangeN(i)-1
            call fix_conserve1(tmpPs,1,ndim,indexChangeStart(i),storeIndex+j, 1)
          enddo
          storeIndex = storeIndex + indexChangeN(i)
        endif
      end do

  end subroutine fix_conserve_vars


  !!this stores the fluxes one at a time, even if ..
  subroutine store_flux_var(flux,indexVar,my_dt, igrid,indexChangeStart, indexChangeN, indexChangeFixC)
    use mod_global_parameters
    double precision, allocatable, intent(in), dimension(:^D&,:) :: flux
    integer, intent(in) :: indexVar,igrid
    double precision, intent(in) :: my_dt
    integer, intent(in), dimension(:) :: indexChangeStart, indexChangeN
    logical, intent(in), dimension(:) :: indexChangeFixC

    integer :: storeIndex,i
    double precision, allocatable, dimension(:^D&,:,:) :: fluxC
    logical :: found

    storeIndex = 0
    found = .false.
    do while (.not. found .and. i .le. size(indexChangeStart))
      if(indexChangeStart(i) .le. indexVar .and. indexVar .le. indexChangeStart(i) + indexChangeN(i) ) then
        storeIndex = storeIndex + indexVar - indexChangeStart(i) + 1 
        found = .true.
      else
        if (indexChangeFixC(i)) then 
            storeIndex = storeIndex + indexChangeN(i)
        endif  
        i=i+1
      endif  
    end do

    if(storeIndex>0) then
      !allocate(fluxC(ixI^S,1,1:ndim))
      allocate(fluxC(^D&size(flux,^D),1,1:ndim))
      fluxC(:^D&,1,1:ndim) = my_dt * flux(:^D&,1:ndim)
      call store_flux1(igrid,fluxC,1,ndim,storeIndex,1)
      deallocate(fluxC)
    endif
  end subroutine store_flux_var





end module mod_fix_conserve
