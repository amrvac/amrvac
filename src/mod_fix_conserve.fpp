!> Module for flux conservation near refinement boundaries
module mod_fix_conserve
  implicit none
  private

  type fluxalloc
     double precision, dimension(:,:,:,:), pointer:: flux => null()
     double precision, dimension(:,:,:,:), pointer:: edge => null()
  end type fluxalloc
  !> store flux to fix conservation
  type(fluxalloc), dimension(:,:,:), allocatable, public :: pflux

  integer, save                        :: nrecv, nsend
  double precision, allocatable, save  :: recvbuffer(:), sendbuffer(:)
  integer, dimension(:), allocatable   :: fc_recvreq, fc_sendreq
  integer, dimension(:,:), allocatable :: fc_recvstat, fc_sendstat
  integer, dimension(3), save        :: isize
  integer                              :: ibuf, ibuf_send
  ! ct for corner total
  integer, save                        :: nrecv_ct, nsend_ct
  ! buffer for corner coarse
  double precision, allocatable, save  :: recvbuffer_cc(:), sendbuffer_cc(:)
  integer, dimension(:), allocatable   :: cc_recvreq, cc_sendreq
  integer, dimension(:,:), allocatable :: cc_recvstat, cc_sendstat
  integer, dimension(3), save        :: isize_stg
  integer                              :: ibuf_cc, ibuf_cc_send
  integer                              :: itag, itag_cc, isend, isend_cc,&
      irecv, irecv_cc

  public :: init_comm_fix_conserve
  public :: allocateBflux
  public :: deallocateBflux
  public :: sendflux
  public :: recvflux
  public :: store_flux
  public :: store_edge
  public :: fix_conserve
  public :: fix_edges

 contains

   subroutine init_comm_fix_conserve(idimmin,idimmax,nwfluxin)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax,nwfluxin

     integer :: iigrid, igrid, idims, iside, i1,i2,i3, nxCo1,nxCo2,nxCo3
     integer :: ic1,ic2,ic3, inc1,inc2,inc3, ipe_neighbor
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

     do idims= idimmin,idimmax
       select case (idims)
         case (1)
         nrecv=nrecv+nrecv_fc(1)
         nsend=nsend+nsend_fc(1)
         nxCo1=1;nxCo2=ixGhi2/2-nghostcells;nxCo3=ixGhi3/2-nghostcells;
         isize(1)=nxCo1*nxCo2*nxCo3*(nwfluxin)
         recvsize=recvsize+nrecv_fc(1)*isize(1)
         sendsize=sendsize+nsend_fc(1)*isize(1)
         if(stagger_grid) then
           ! This does not consider the 'coarse corner' case
           nxCo1=1;nxCo2=ixGhi2/2-nghostcells+1;nxCo3=ixGhi3/2-nghostcells+1;
           isize_stg(1)=nxCo1*nxCo2*nxCo3*(3-1)
           ! the whole size is used (cell centered and staggered)
           isize(1)=isize(1)+isize_stg(1)      
           recvsize=recvsize+nrecv_fc(1)*isize_stg(1)
           sendsize=sendsize+nsend_fc(1)*isize_stg(1)
           ! Coarse corner case
           nrecv_ct=nrecv_ct+nrecv_cc(1)
           nsend_ct=nsend_ct+nsend_cc(1)
           recvsize_cc=recvsize_cc+nrecv_cc(1)*isize_stg(1)
           sendsize_cc=sendsize_cc+nsend_cc(1)*isize_stg(1)
         end if
         
         case (2)
         nrecv=nrecv+nrecv_fc(2)
         nsend=nsend+nsend_fc(2)
         nxCo1=ixGhi1/2-nghostcells;nxCo2=1;nxCo3=ixGhi3/2-nghostcells;
         isize(2)=nxCo1*nxCo2*nxCo3*(nwfluxin)
         recvsize=recvsize+nrecv_fc(2)*isize(2)
         sendsize=sendsize+nsend_fc(2)*isize(2)
         if(stagger_grid) then
           ! This does not consider the 'coarse corner' case
           nxCo1=ixGhi1/2-nghostcells+1;nxCo2=1;nxCo3=ixGhi3/2-nghostcells+1;
           isize_stg(2)=nxCo1*nxCo2*nxCo3*(3-1)
           ! the whole size is used (cell centered and staggered)
           isize(2)=isize(2)+isize_stg(2)      
           recvsize=recvsize+nrecv_fc(2)*isize_stg(2)
           sendsize=sendsize+nsend_fc(2)*isize_stg(2)
           ! Coarse corner case
           nrecv_ct=nrecv_ct+nrecv_cc(2)
           nsend_ct=nsend_ct+nsend_cc(2)
           recvsize_cc=recvsize_cc+nrecv_cc(2)*isize_stg(2)
           sendsize_cc=sendsize_cc+nsend_cc(2)*isize_stg(2)
         end if
         
         case (3)
         nrecv=nrecv+nrecv_fc(3)
         nsend=nsend+nsend_fc(3)
         nxCo1=ixGhi1/2-nghostcells;nxCo2=ixGhi2/2-nghostcells;nxCo3=1;
         isize(3)=nxCo1*nxCo2*nxCo3*(nwfluxin)
         recvsize=recvsize+nrecv_fc(3)*isize(3)
         sendsize=sendsize+nsend_fc(3)*isize(3)
         if(stagger_grid) then
           ! This does not consider the 'coarse corner' case
           nxCo1=ixGhi1/2-nghostcells+1;nxCo2=ixGhi2/2-nghostcells+1;nxCo3=1;
           isize_stg(3)=nxCo1*nxCo2*nxCo3*(3-1)
           ! the whole size is used (cell centered and staggered)
           isize(3)=isize(3)+isize_stg(3)      
           recvsize=recvsize+nrecv_fc(3)*isize_stg(3)
           sendsize=sendsize+nsend_fc(3)*isize_stg(3)
           ! Coarse corner case
           nrecv_ct=nrecv_ct+nrecv_cc(3)
           nsend_ct=nsend_ct+nsend_cc(3)
           recvsize_cc=recvsize_cc+nrecv_cc(3)*isize_stg(3)
           sendsize_cc=sendsize_cc+nsend_cc(3)*isize_stg(3)
         end if
         
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
           allocate(cc_recvstat(MPI_STATUS_SIZE,nrecv_ct),&
               cc_recvreq(nrecv_ct))
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
           allocate(cc_sendstat(MPI_STATUS_SIZE,nsend_ct),&
               cc_sendreq(nsend_ct))
         end if
       else
         allocate(cc_sendstat(MPI_STATUS_SIZE,nsend_ct), cc_sendreq(nsend_ct))
       end if
     end if

   end subroutine init_comm_fix_conserve

   subroutine recvflux(idimmin,idimmax)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax

     integer :: iigrid, igrid, idims, iside, i1,i2,i3, nxCo1,nxCo2,nxCo3
     integer :: ic1,ic2,ic3, inc1,inc2,inc3, ipe_neighbor
     integer :: pi1,pi2,pi3,mi1,mi2,mi3,ph1,ph2,ph3,mh1,mh2,mh3,idir

     if (nrecv>0) then
       fc_recvreq=MPI_REQUEST_NULL
       ibuf=1
       irecv=0

       do iigrid=1,igridstail; igrid=igrids(iigrid);
         do idims= idimmin,idimmax
           do iside=1,2
             i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
             i3=kr(3,idims)*(2*iside-3);

             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

             if (neighbor_type(i1,i2,i3,igrid)/=4) cycle
             do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
             do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
             do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               if (ipe_neighbor/=mype) then
                 irecv=irecv+1
                 itag=4**3*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)+inc3*4**(3-1)
                 call MPI_IRECV(recvbuffer(ibuf),isize(idims),&
                     MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                    fc_recvreq(irecv),ierrmpi)
                 ibuf=ibuf+isize(idims)
               end if
             end do
             end do
             end do
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
           do idims= idimmin,idimmax
             do iside=1,2
               i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
               i3=kr(3,idims)*(2*iside-3);
               ! Check if there are special corners
               ! (Coarse block diagonal to a fine block)
               ! If there are, receive.
               ! Tags are calculated in the same way as for
               ! normal fluxes, but should not overlap because
               ! inc^D are different
               if (neighbor_type(i1,i2,i3,igrid)==3) then
                 do idir=idims+1,ndim
                   pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
                   mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
                   ph1=pi1-kr(idims,1)*(2*iside-3)
                   ph2=pi2-kr(idims,2)*(2*iside-3)
                   ph3=pi3-kr(idims,3)*(2*iside-3);
                   mh1=mi1-kr(idims,1)*(2*iside-3)
                   mh2=mi2-kr(idims,2)*(2*iside-3)
                   mh3=mi3-kr(idims,3)*(2*iside-3);

                   if (neighbor_type(pi1,pi2,pi3,&
                      igrid)==4.and.neighbor_type(ph1,ph2,ph3,&
                      igrid)==3.and.neighbor_pole(pi1,pi2,pi3,igrid)==0) then
                      ! Loop on children (several in 3D)
                    do ic3=1+int((1-pi3)/2),2-int((1+pi3)/2)
                       inc3=2*pi3+ic3
                    do ic2=1+int((1-pi2)/2),2-int((1+pi2)/2)
                       inc2=2*pi2+ic2
                    do ic1=1+int((1-pi1)/2),2-int((1+pi1)/2)
                       inc1=2*pi1+ic1
                       ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                       if (mype/=ipe_neighbor) then
                         irecv_cc=irecv_cc+1
                         itag_cc=4**3*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                            inc3*4**(3-1)
                         call MPI_IRECV(recvbuffer_cc(ibuf_cc),&
                            isize_stg(idims),MPI_DOUBLE_PRECISION,ipe_neighbor,&
                            itag_cc,icomm,cc_recvreq(irecv_cc),ierrmpi)
                         ibuf_cc=ibuf_cc+isize_stg(idims)
                       end if
                    end do
                    end do
                    end do
                   end if

                   if (neighbor_type(mi1,mi2,mi3,&
                      igrid)==4.and.neighbor_type(mh1,mh2,mh3,&
                      igrid)==3.and.neighbor_pole(mi1,mi2,mi3,igrid)==0) then
                      ! Loop on children (several in 3D)
                    do ic3=1+int((1-mi3)/2),2-int((1+mi3)/2)
                        inc3=2*mi3+ic3
                    do ic2=1+int((1-mi2)/2),2-int((1+mi2)/2)
                        inc2=2*mi2+ic2
                    do ic1=1+int((1-mi1)/2),2-int((1+mi1)/2)
                        inc1=2*mi1+ic1
                       ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                       if (mype/=ipe_neighbor) then
                         irecv_cc=irecv_cc+1
                         itag_cc=4**3*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                            inc3*4**(3-1)
                         call MPI_IRECV(recvbuffer_cc(ibuf_cc),&
                            isize_stg(idims),MPI_DOUBLE_PRECISION,ipe_neighbor,&
                            itag_cc,icomm,cc_recvreq(irecv_cc),ierrmpi)
                         ibuf_cc=ibuf_cc+isize_stg(idims)
                       end if
                    end do
                    end do
                    end do
                   end if
                 end do
               end if
             end do
           end do
         end do
       end if
     end if ! end if stagger grid

   end subroutine recvflux

   subroutine sendflux(idimmin,idimmax)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax

     integer :: idims, iside, i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3, ix1,ix2,&
        ix3, ixCo1,ixCo2,ixCo3, nxCo1,nxCo2,nxCo3, iw
     integer :: ineighbor, ipe_neighbor, igrid, iigrid, ibuf_send_next
     integer :: idir, ibuf_cc_send_next, pi1,pi2,pi3, ph1,ph2,ph3, mi1,mi2,mi3,&
         mh1,mh2,mh3

     fc_sendreq = MPI_REQUEST_NULL
     isend      = 0
     if(stagger_grid) then
       ibuf_send  = 1
       cc_sendreq=MPI_REQUEST_NULL
       isend_cc=0
       ibuf_cc_send=1
     end if

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims = idimmin,idimmax
         select case (idims)
        case (1)
           do iside=1,2
             i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3)
             i3=kr(3,1)*(2*iside-3);

             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

             if (neighbor_type(i1,i2,i3,igrid)==neighbor_coarse) then
               ! send flux to coarser neighbor
               ineighbor=neighbor(1,i1,i2,i3,igrid)
               ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
               if (ipe_neighbor/=mype) then
                 ic1=1+modulo(node(pig1_,igrid)-1,2)
                 ic2=1+modulo(node(pig2_,igrid)-1,2)
                 ic3=1+modulo(node(pig3_,igrid)-1,2);
                 inc1=-2*i1+ic1;inc2=ic2;inc3=ic3;
                 itag=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                    inc3*4**(3-1)
                 isend=isend+1

                 if(stagger_grid) then
                   ibuf_send_next=ibuf_send+isize(1)
                   sendbuffer(ibuf_send:ibuf_send_next-isize_stg(1)-&
                      1)=reshape(pflux(iside,1,igrid)%flux,&
                      (/isize(1)-isize_stg(1)/))

                   sendbuffer(ibuf_send_next-isize_stg(1):ibuf_send_next-&
                      1)=reshape(pflux(iside,1,igrid)%edge,(/isize_stg(1)/))
                   call MPI_ISEND(sendbuffer(ibuf_send),isize(1),&
                       MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                      fc_sendreq(isend),ierrmpi)
                   ibuf_send=ibuf_send_next
                 else
                   call MPI_ISEND(pflux(iside,1,igrid)%flux,isize(1),&
                       MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                      fc_sendreq(isend),ierrmpi)
                 end if
               end if

               if(stagger_grid) then
                 ! If we are in a fine block surrounded by coarse blocks
                 do idir=idims+1,ndim
                   pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
                   mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
                   ph1=pi1-kr(idims,1)*(2*iside-3)
                   ph2=pi2-kr(idims,2)*(2*iside-3)
                   ph3=pi3-kr(idims,3)*(2*iside-3);
                   mh1=mi1-kr(idims,1)*(2*iside-3)
                   mh2=mi2-kr(idims,2)*(2*iside-3)
                   mh3=mi3-kr(idims,3)*(2*iside-3);

                   if (neighbor_type(pi1,pi2,pi3,&
                      igrid)==2.and.neighbor_type(ph1,ph2,ph3,&
                      igrid)==2.and.mype/=neighbor(2,pi1,pi2,pi3,&
                      igrid).and.neighbor_pole(pi1,pi2,pi3,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,pi1,pi2,pi3,igrid)
                     ipe_neighbor=neighbor(2,pi1,pi2,pi3,igrid)
                     ic1=1+modulo(node(pig1_,igrid)-1,2)
                     ic2=1+modulo(node(pig2_,igrid)-1,2)
                     ic3=1+modulo(node(pig3_,igrid)-1,2);
                     inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                     itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                        inc3*4**(3-1)
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(1)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
                        1)=reshape(pflux(iside,1,igrid)%edge,&
                        shape=(/isize_stg(1)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(1),&
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
                        cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
       
                   if (neighbor_type(mi1,mi2,mi3,&
                      igrid)==2.and.neighbor_type(mh1,mh2,mh3,&
                      igrid)==2.and.mype/=neighbor(2,mi1,mi2,mi3,&
                      igrid).and.neighbor_pole(mi1,mi2,mi3,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,mi1,mi2,mi3,igrid)
                     ipe_neighbor=neighbor(2,mi1,mi2,mi3,igrid)
                     ic1=1+modulo(node(pig1_,igrid)-1,2)
                     ic2=1+modulo(node(pig2_,igrid)-1,2)
                     ic3=1+modulo(node(pig3_,igrid)-1,2);
                     inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                     inc1=-2*mi1+ic1;inc2=-2*mi2+ic2;inc3=-2*mi3+ic3;
                     itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                        inc3*4**(3-1)
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(1)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
                        1)=reshape(pflux(iside,1,igrid)%edge,&
                        shape=(/isize_stg(1)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(1),&
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
                        cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
                 end do
               end if ! end if stagger grid

             end if
           end do
        case (2)
           do iside=1,2
             i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3)
             i3=kr(3,2)*(2*iside-3);

             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

             if (neighbor_type(i1,i2,i3,igrid)==neighbor_coarse) then
               ! send flux to coarser neighbor
               ineighbor=neighbor(1,i1,i2,i3,igrid)
               ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
               if (ipe_neighbor/=mype) then
                 ic1=1+modulo(node(pig1_,igrid)-1,2)
                 ic2=1+modulo(node(pig2_,igrid)-1,2)
                 ic3=1+modulo(node(pig3_,igrid)-1,2);
                 inc1=ic1;inc2=-2*i2+ic2;inc3=ic3;
                 itag=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                    inc3*4**(3-1)
                 isend=isend+1

                 if(stagger_grid) then
                   ibuf_send_next=ibuf_send+isize(2)
                   sendbuffer(ibuf_send:ibuf_send_next-isize_stg(2)-&
                      1)=reshape(pflux(iside,2,igrid)%flux,&
                      (/isize(2)-isize_stg(2)/))

                   sendbuffer(ibuf_send_next-isize_stg(2):ibuf_send_next-&
                      1)=reshape(pflux(iside,2,igrid)%edge,(/isize_stg(2)/))
                   call MPI_ISEND(sendbuffer(ibuf_send),isize(2),&
                       MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                      fc_sendreq(isend),ierrmpi)
                   ibuf_send=ibuf_send_next
                 else
                   call MPI_ISEND(pflux(iside,2,igrid)%flux,isize(2),&
                       MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                      fc_sendreq(isend),ierrmpi)
                 end if
               end if

               if(stagger_grid) then
                 ! If we are in a fine block surrounded by coarse blocks
                 do idir=idims+1,ndim
                   pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
                   mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
                   ph1=pi1-kr(idims,1)*(2*iside-3)
                   ph2=pi2-kr(idims,2)*(2*iside-3)
                   ph3=pi3-kr(idims,3)*(2*iside-3);
                   mh1=mi1-kr(idims,1)*(2*iside-3)
                   mh2=mi2-kr(idims,2)*(2*iside-3)
                   mh3=mi3-kr(idims,3)*(2*iside-3);

                   if (neighbor_type(pi1,pi2,pi3,&
                      igrid)==2.and.neighbor_type(ph1,ph2,ph3,&
                      igrid)==2.and.mype/=neighbor(2,pi1,pi2,pi3,&
                      igrid).and.neighbor_pole(pi1,pi2,pi3,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,pi1,pi2,pi3,igrid)
                     ipe_neighbor=neighbor(2,pi1,pi2,pi3,igrid)
                     ic1=1+modulo(node(pig1_,igrid)-1,2)
                     ic2=1+modulo(node(pig2_,igrid)-1,2)
                     ic3=1+modulo(node(pig3_,igrid)-1,2);
                     inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                     itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                        inc3*4**(3-1)
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(2)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
                        1)=reshape(pflux(iside,2,igrid)%edge,&
                        shape=(/isize_stg(2)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(2),&
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
                        cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
       
                   if (neighbor_type(mi1,mi2,mi3,&
                      igrid)==2.and.neighbor_type(mh1,mh2,mh3,&
                      igrid)==2.and.mype/=neighbor(2,mi1,mi2,mi3,&
                      igrid).and.neighbor_pole(mi1,mi2,mi3,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,mi1,mi2,mi3,igrid)
                     ipe_neighbor=neighbor(2,mi1,mi2,mi3,igrid)
                     ic1=1+modulo(node(pig1_,igrid)-1,2)
                     ic2=1+modulo(node(pig2_,igrid)-1,2)
                     ic3=1+modulo(node(pig3_,igrid)-1,2);
                     inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                     inc1=-2*mi1+ic1;inc2=-2*mi2+ic2;inc3=-2*mi3+ic3;
                     itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                        inc3*4**(3-1)
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(2)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
                        1)=reshape(pflux(iside,2,igrid)%edge,&
                        shape=(/isize_stg(2)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(2),&
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
                        cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
                 end do
               end if ! end if stagger grid

             end if
           end do
        case (3)
           do iside=1,2
             i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3)
             i3=kr(3,3)*(2*iside-3);

             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

             if (neighbor_type(i1,i2,i3,igrid)==neighbor_coarse) then
               ! send flux to coarser neighbor
               ineighbor=neighbor(1,i1,i2,i3,igrid)
               ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
               if (ipe_neighbor/=mype) then
                 ic1=1+modulo(node(pig1_,igrid)-1,2)
                 ic2=1+modulo(node(pig2_,igrid)-1,2)
                 ic3=1+modulo(node(pig3_,igrid)-1,2);
                 inc1=ic1;inc2=ic2;inc3=-2*i3+ic3;
                 itag=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                    inc3*4**(3-1)
                 isend=isend+1

                 if(stagger_grid) then
                   ibuf_send_next=ibuf_send+isize(3)
                   sendbuffer(ibuf_send:ibuf_send_next-isize_stg(3)-&
                      1)=reshape(pflux(iside,3,igrid)%flux,&
                      (/isize(3)-isize_stg(3)/))

                   sendbuffer(ibuf_send_next-isize_stg(3):ibuf_send_next-&
                      1)=reshape(pflux(iside,3,igrid)%edge,(/isize_stg(3)/))
                   call MPI_ISEND(sendbuffer(ibuf_send),isize(3),&
                       MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                      fc_sendreq(isend),ierrmpi)
                   ibuf_send=ibuf_send_next
                 else
                   call MPI_ISEND(pflux(iside,3,igrid)%flux,isize(3),&
                       MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                      fc_sendreq(isend),ierrmpi)
                 end if
               end if

               if(stagger_grid) then
                 ! If we are in a fine block surrounded by coarse blocks
                 do idir=idims+1,ndim
                   pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
                   mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
                   ph1=pi1-kr(idims,1)*(2*iside-3)
                   ph2=pi2-kr(idims,2)*(2*iside-3)
                   ph3=pi3-kr(idims,3)*(2*iside-3);
                   mh1=mi1-kr(idims,1)*(2*iside-3)
                   mh2=mi2-kr(idims,2)*(2*iside-3)
                   mh3=mi3-kr(idims,3)*(2*iside-3);

                   if (neighbor_type(pi1,pi2,pi3,&
                      igrid)==2.and.neighbor_type(ph1,ph2,ph3,&
                      igrid)==2.and.mype/=neighbor(2,pi1,pi2,pi3,&
                      igrid).and.neighbor_pole(pi1,pi2,pi3,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,pi1,pi2,pi3,igrid)
                     ipe_neighbor=neighbor(2,pi1,pi2,pi3,igrid)
                     ic1=1+modulo(node(pig1_,igrid)-1,2)
                     ic2=1+modulo(node(pig2_,igrid)-1,2)
                     ic3=1+modulo(node(pig3_,igrid)-1,2);
                     inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                     itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                        inc3*4**(3-1)
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(3)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
                        1)=reshape(pflux(iside,3,igrid)%edge,&
                        shape=(/isize_stg(3)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(3),&
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
                        cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
       
                   if (neighbor_type(mi1,mi2,mi3,&
                      igrid)==2.and.neighbor_type(mh1,mh2,mh3,&
                      igrid)==2.and.mype/=neighbor(2,mi1,mi2,mi3,&
                      igrid).and.neighbor_pole(mi1,mi2,mi3,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,mi1,mi2,mi3,igrid)
                     ipe_neighbor=neighbor(2,mi1,mi2,mi3,igrid)
                     ic1=1+modulo(node(pig1_,igrid)-1,2)
                     ic2=1+modulo(node(pig2_,igrid)-1,2)
                     ic3=1+modulo(node(pig3_,igrid)-1,2);
                     inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
                     inc1=-2*mi1+ic1;inc2=-2*mi2+ic2;inc3=-2*mi3+ic3;
                     itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                        inc3*4**(3-1)
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(3)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
                        1)=reshape(pflux(iside,3,igrid)%edge,&
                        shape=(/isize_stg(3)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(3),&
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
                        cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
                 end do
               end if ! end if stagger grid

             end if
           end do
         end select
       end do
     end do
   end subroutine sendflux

   subroutine allocateBflux
     use mod_global_parameters

     integer :: iigrid, igrid, iside, i1,i2,i3, nx1,nx2,nx3, nxCo1,nxCo2,nxCo3
     integer :: idir,idim,pi1,pi2,pi3, mi1,mi2,mi3, ph1,ph2,ph3, mh1,mh2,mh3 !To detect corners

     nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
     nxCo1=nx1/2;nxCo2=nx2/2;nxCo3=nx3/2;

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! For every grid,
       ! arrays for the fluxes are allocated for every face direction(^D)
       ! and every side (1=left, 2=right)
       do iside=1,2
         i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);i3=kr(3,1)*(2*iside-3);

         if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

         select case (neighbor_type(i1,i2,i3,igrid))
         case(neighbor_fine)
           allocate(pflux(iside,1,igrid)%flux(1,1:nx2,1:nx3,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,1,igrid)%edge(1,0:nx2,0:nx3,&
              1:ndim-1))
         case(neighbor_coarse)
           allocate(pflux(iside,1,igrid)%flux(1,1:nxCo2,1:nxCo3,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,1,igrid)%edge(1,0:nxCo2,&
              0:nxCo3,1:ndim-1))
         case(neighbor_sibling)
           if(stagger_grid) then
             idim=1
             do idir=idim+1,ndim
             !do idir=min(idim+1,ndim),ndim
               pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
               mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
               ph1=pi1-kr(1,1)*(2*iside-3);ph2=pi2-kr(1,2)*(2*iside-3)
               ph3=pi3-kr(1,3)*(2*iside-3);
               mh1=mi1-kr(1,1)*(2*iside-3);mh2=mi2-kr(1,2)*(2*iside-3)
               mh3=mi3-kr(1,3)*(2*iside-3);
               if ((neighbor_type(pi1,pi2,pi3,igrid)==4.and.neighbor_type(ph1,&
                  ph2,ph3,igrid)==3).or.(neighbor_type(mi1,mi2,mi3,&
                  igrid)==4.and.neighbor_type(mh1,mh2,mh3,igrid)==3)) then
                 allocate(pflux(iside,1,igrid)%edge(1,0:nx2,0:nx3,1:ndim-1))
                 exit
               end if
             end do
           end if
         end select
       end do
       do iside=1,2
         i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);i3=kr(3,2)*(2*iside-3);

         if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

         select case (neighbor_type(i1,i2,i3,igrid))
         case(neighbor_fine)
           allocate(pflux(iside,2,igrid)%flux(1:nx1,1,1:nx3,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,2,igrid)%edge(0:nx1,1,0:nx3,&
              1:ndim-1))
         case(neighbor_coarse)
           allocate(pflux(iside,2,igrid)%flux(1:nxCo1,1,1:nxCo3,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,2,igrid)%edge(0:nxCo1,1,&
              0:nxCo3,1:ndim-1))
         case(neighbor_sibling)
           if(stagger_grid) then
             idim=2
             do idir=idim+1,ndim
             !do idir=min(idim+1,ndim),ndim
               pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
               mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
               ph1=pi1-kr(2,1)*(2*iside-3);ph2=pi2-kr(2,2)*(2*iside-3)
               ph3=pi3-kr(2,3)*(2*iside-3);
               mh1=mi1-kr(2,1)*(2*iside-3);mh2=mi2-kr(2,2)*(2*iside-3)
               mh3=mi3-kr(2,3)*(2*iside-3);
               if ((neighbor_type(pi1,pi2,pi3,igrid)==4.and.neighbor_type(ph1,&
                  ph2,ph3,igrid)==3).or.(neighbor_type(mi1,mi2,mi3,&
                  igrid)==4.and.neighbor_type(mh1,mh2,mh3,igrid)==3)) then
                 allocate(pflux(iside,2,igrid)%edge(0:nx1,1,0:nx3,1:ndim-1))
                 exit
               end if
             end do
           end if
         end select
       end do
       do iside=1,2
         i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3);i3=kr(3,3)*(2*iside-3);

         if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

         select case (neighbor_type(i1,i2,i3,igrid))
         case(neighbor_fine)
           allocate(pflux(iside,3,igrid)%flux(1:nx1,1:nx2,1,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,3,igrid)%edge(0:nx1,0:nx2,1,&
              1:ndim-1))
         case(neighbor_coarse)
           allocate(pflux(iside,3,igrid)%flux(1:nxCo1,1:nxCo2,1,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,3,igrid)%edge(0:nxCo1,0:nxCo2,&
              1,1:ndim-1))
         case(neighbor_sibling)
           if(stagger_grid) then
             idim=3
             do idir=idim+1,ndim
             !do idir=min(idim+1,ndim),ndim
               pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
               mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
               ph1=pi1-kr(3,1)*(2*iside-3);ph2=pi2-kr(3,2)*(2*iside-3)
               ph3=pi3-kr(3,3)*(2*iside-3);
               mh1=mi1-kr(3,1)*(2*iside-3);mh2=mi2-kr(3,2)*(2*iside-3)
               mh3=mi3-kr(3,3)*(2*iside-3);
               if ((neighbor_type(pi1,pi2,pi3,igrid)==4.and.neighbor_type(ph1,&
                  ph2,ph3,igrid)==3).or.(neighbor_type(mi1,mi2,mi3,&
                  igrid)==4.and.neighbor_type(mh1,mh2,mh3,igrid)==3)) then
                 allocate(pflux(iside,3,igrid)%edge(0:nx1,0:nx2,1,1:ndim-1))
                 exit
               end if
             end do
           end if
         end select
       end do
     end do

   end subroutine allocateBflux

   subroutine deallocateBflux
     use mod_global_parameters

     integer :: iigrid, igrid, iside

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do iside=1,2
         if (associated(pflux(iside,1,igrid)%flux)) then
           deallocate(pflux(iside,1,igrid)%flux)
           nullify(pflux(iside,1,igrid)%flux)
         end if
         if (associated(pflux(iside,1,igrid)%edge)) then
           deallocate(pflux(iside,1,igrid)%edge)
           nullify(pflux(iside,1,igrid)%edge)
         end if
       end do
       do iside=1,2
         if (associated(pflux(iside,2,igrid)%flux)) then
           deallocate(pflux(iside,2,igrid)%flux)
           nullify(pflux(iside,2,igrid)%flux)
         end if
         if (associated(pflux(iside,2,igrid)%edge)) then
           deallocate(pflux(iside,2,igrid)%edge)
           nullify(pflux(iside,2,igrid)%edge)
         end if
       end do
       do iside=1,2
         if (associated(pflux(iside,3,igrid)%flux)) then
           deallocate(pflux(iside,3,igrid)%flux)
           nullify(pflux(iside,3,igrid)%flux)
         end if
         if (associated(pflux(iside,3,igrid)%edge)) then
           deallocate(pflux(iside,3,igrid)%edge)
           nullify(pflux(iside,3,igrid)%edge)
         end if
       end do
     end do

   end subroutine deallocateBflux

   subroutine fix_conserve(psb,idimmin,idimmax,nw0,nwfluxin)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax, nw0, nwfluxin
     type(state) :: psb(max_blocks)

     integer :: iigrid, igrid, idims, iside, iotherside, i1,i2,i3, ic1,ic2,ic3,&
         inc1,inc2,inc3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
     integer :: nxCo1,nxCo2,nxCo3, iw, ix, ipe_neighbor, ineighbor, nbuf,&
         ibufnext, nw1
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

     nxCo1=(ixMhi1-ixMlo1+1)/2;nxCo2=(ixMhi2-ixMlo2+1)/2
     nxCo3=(ixMhi3-ixMlo3+1)/2;

     ! for all grids: perform flux update at Coarse-Fine interfaces
     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims= idimmin,idimmax
         select case (idims)
           case (1)
           do iside=1,2
             i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3)
             i3=kr(3,1)*(2*iside-3);

             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

             if (neighbor_type(i1,i2,i3,igrid)/=4) cycle

 !opedit: skip over active/passive interface since flux for passive ones is 
             ! not computed, keep the buffer counter up to date:
             if (.not.neighbor_active(i1,i2,i3,&
                igrid).or..not.neighbor_active(0,0,0,igrid) ) then
               do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
           do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               if (ipe_neighbor/=mype) then
                 ibufnext=ibuf+isize(1)
                 ibuf=ibufnext
               end if
               end do
           end do
           end do
               cycle
             end if
             !

             select case (iside)
             case (1)
               ix=ixMlo1
             case (2)
               ix=ixMhi1
             end select

             ! remove coarse flux
             if (slab_uniform) then
               psb(igrid)%w(ix,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
                  nw0:nw1) = psb(igrid)%w(ix,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
                  nw0:nw1) -pflux(iside,1,igrid)%flux(1,:,:,1:nwfluxin)
             else
               do iw=nw0,nw1
                 psb(igrid)%w(ix,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
                    iw)=psb(igrid)%w(ix,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
                    iw)-pflux(iside,1,igrid)%flux(1,:,:,&
                    iw-nw0+1) /ps(igrid)%dvolume(ix,ixMlo2:ixMhi2,&
                    ixMlo3:ixMhi3)
               end do
             end if


             ! add fine flux
            do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
           do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               ixmin1=ix;ixmin2=ixMlo2+(ic2-1)*nxCo2
               ixmin3=ixMlo3+(ic3-1)*nxCo3;
               ixmax1=ix;ixmax2=ixmin2-1+nxCo2;ixmax3=ixmin3-1+nxCo3;
               if (ipe_neighbor==mype) then
                 iotherside=3-iside
                 if (slab_uniform) then
                   psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                      nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      ixmin3:ixmax3,nw0:nw1) + pflux(iotherside,1,&
                      ineighbor)%flux(:,:,:,1:nwfluxin)* CoFiratio
                 else
                   do iw=nw0,nw1
                     psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                        iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,iw) +pflux(iotherside,1,&
                        ineighbor)%flux(:,:,:,&
                        iw-nw0+1) /ps(igrid)%dvolume(ixmin1:ixmax1,&
                        ixmin2:ixmax2,ixmin3:ixmax3)
                   end do
                 end if
               else
                 if (slab_uniform) then
                   ibufnext=ibuf+isize(1)
                   if(stagger_grid) ibufnext=ibufnext-isize_stg(1)
                   psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                      nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      ixmin3:ixmax3,nw0:nw1)+CoFiratio &
                      *reshape(source=recvbuffer(ibuf:ibufnext-1),&
                       shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      ixmin3:ixmax3,nw0:nw1)))
                   ibuf=ibuf+isize(1)
                 else
                   ibufnext=ibuf+isize(1)
                   if(stagger_grid) then
                     nbuf=(isize(1)-isize_stg(1))/nwfluxin
                   else
                     nbuf=isize(1)/nwfluxin
                   end if
                   do iw=nw0,nw1
                     psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                        iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,iw) &
                        +reshape(source=recvbuffer(ibuf:ibufnext-1),&
                         shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,iw))) /ps(igrid)%dvolume(ixmin1:ixmax1,&
                        ixmin2:ixmax2,ixmin3:ixmax3)
                     ibuf=ibuf+nbuf
                   end do
                   ibuf=ibufnext
                 end if
               end if
            end do
           end do
           end do
           end do
           case (2)
           do iside=1,2
             i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3)
             i3=kr(3,2)*(2*iside-3);

             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

             if (neighbor_type(i1,i2,i3,igrid)/=4) cycle

 !opedit: skip over active/passive interface since flux for passive ones is 
             ! not computed, keep the buffer counter up to date:
             if (.not.neighbor_active(i1,i2,i3,&
                igrid).or..not.neighbor_active(0,0,0,igrid) ) then
               do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
           do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               if (ipe_neighbor/=mype) then
                 ibufnext=ibuf+isize(2)
                 ibuf=ibufnext
               end if
               end do
           end do
           end do
               cycle
             end if
             !

             select case (iside)
             case (1)
               ix=ixMlo2
             case (2)
               ix=ixMhi2
             end select

             ! remove coarse flux
             if (slab_uniform) then
               psb(igrid)%w(ixMlo1:ixMhi1,ix,ixMlo3:ixMhi3,&
                  nw0:nw1) = psb(igrid)%w(ixMlo1:ixMhi1,ix,ixMlo3:ixMhi3,&
                  nw0:nw1) -pflux(iside,2,igrid)%flux(:,1,:,1:nwfluxin)
             else
               do iw=nw0,nw1
                 psb(igrid)%w(ixMlo1:ixMhi1,ix,ixMlo3:ixMhi3,&
                    iw)=psb(igrid)%w(ixMlo1:ixMhi1,ix,ixMlo3:ixMhi3,&
                    iw)-pflux(iside,2,igrid)%flux(:,1,:,&
                    iw-nw0+1) /ps(igrid)%dvolume(ixMlo1:ixMhi1,ix,&
                    ixMlo3:ixMhi3)
               end do
             end if


             ! add fine flux
            do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
           do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               ixmin1=ixMlo1+(ic1-1)*nxCo1;ixmin2=ix
               ixmin3=ixMlo3+(ic3-1)*nxCo3;
               ixmax1=ixmin1-1+nxCo1;ixmax2=ix;ixmax3=ixmin3-1+nxCo3;
               if (ipe_neighbor==mype) then
                 iotherside=3-iside
                 if (slab_uniform) then
                   psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                      nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      ixmin3:ixmax3,nw0:nw1) + pflux(iotherside,2,&
                      ineighbor)%flux(:,:,:,1:nwfluxin)* CoFiratio
                 else
                   do iw=nw0,nw1
                     psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                        iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,iw) +pflux(iotherside,2,&
                        ineighbor)%flux(:,:,:,&
                        iw-nw0+1) /ps(igrid)%dvolume(ixmin1:ixmax1,&
                        ixmin2:ixmax2,ixmin3:ixmax3)
                   end do
                 end if
               else
                 if (slab_uniform) then
                   ibufnext=ibuf+isize(2)
                   if(stagger_grid) ibufnext=ibufnext-isize_stg(2)
                   psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                      nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      ixmin3:ixmax3,nw0:nw1)+CoFiratio &
                      *reshape(source=recvbuffer(ibuf:ibufnext-1),&
                       shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      ixmin3:ixmax3,nw0:nw1)))
                   ibuf=ibuf+isize(2)
                 else
                   ibufnext=ibuf+isize(2)
                   if(stagger_grid) then
                     nbuf=(isize(2)-isize_stg(2))/nwfluxin
                   else
                     nbuf=isize(2)/nwfluxin
                   end if
                   do iw=nw0,nw1
                     psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                        iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,iw) &
                        +reshape(source=recvbuffer(ibuf:ibufnext-1),&
                         shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,iw))) /ps(igrid)%dvolume(ixmin1:ixmax1,&
                        ixmin2:ixmax2,ixmin3:ixmax3)
                     ibuf=ibuf+nbuf
                   end do
                   ibuf=ibufnext
                 end if
               end if
            end do
           end do
           end do
           end do
           case (3)
           do iside=1,2
             i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3)
             i3=kr(3,3)*(2*iside-3);

             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

             if (neighbor_type(i1,i2,i3,igrid)/=4) cycle

 !opedit: skip over active/passive interface since flux for passive ones is 
             ! not computed, keep the buffer counter up to date:
             if (.not.neighbor_active(i1,i2,i3,&
                igrid).or..not.neighbor_active(0,0,0,igrid) ) then
               do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
           do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               if (ipe_neighbor/=mype) then
                 ibufnext=ibuf+isize(3)
                 ibuf=ibufnext
               end if
               end do
           end do
           end do
               cycle
             end if
             !

             select case (iside)
             case (1)
               ix=ixMlo3
             case (2)
               ix=ixMhi3
             end select

             ! remove coarse flux
             if (slab_uniform) then
               psb(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ix,&
                  nw0:nw1) = psb(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ix,&
                  nw0:nw1) -pflux(iside,3,igrid)%flux(:,:,1,1:nwfluxin)
             else
               do iw=nw0,nw1
                 psb(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ix,&
                    iw)=psb(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ix,&
                    iw)-pflux(iside,3,igrid)%flux(:,:,1,&
                    iw-nw0+1) /ps(igrid)%dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                    ix)
               end do
             end if


             ! add fine flux
            do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
           do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               ixmin1=ixMlo1+(ic1-1)*nxCo1;ixmin2=ixMlo2+(ic2-1)*nxCo2
               ixmin3=ix;
               ixmax1=ixmin1-1+nxCo1;ixmax2=ixmin2-1+nxCo2;ixmax3=ix;
               if (ipe_neighbor==mype) then
                 iotherside=3-iside
                 if (slab_uniform) then
                   psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                      nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      ixmin3:ixmax3,nw0:nw1) + pflux(iotherside,3,&
                      ineighbor)%flux(:,:,:,1:nwfluxin)* CoFiratio
                 else
                   do iw=nw0,nw1
                     psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                        iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,iw) +pflux(iotherside,3,&
                        ineighbor)%flux(:,:,:,&
                        iw-nw0+1) /ps(igrid)%dvolume(ixmin1:ixmax1,&
                        ixmin2:ixmax2,ixmin3:ixmax3)
                   end do
                 end if
               else
                 if (slab_uniform) then
                   ibufnext=ibuf+isize(3)
                   if(stagger_grid) ibufnext=ibufnext-isize_stg(3)
                   psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                      nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      ixmin3:ixmax3,nw0:nw1)+CoFiratio &
                      *reshape(source=recvbuffer(ibuf:ibufnext-1),&
                       shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      ixmin3:ixmax3,nw0:nw1)))
                   ibuf=ibuf+isize(3)
                 else
                   ibufnext=ibuf+isize(3)
                   if(stagger_grid) then
                     nbuf=(isize(3)-isize_stg(3))/nwfluxin
                   else
                     nbuf=isize(3)/nwfluxin
                   end if
                   do iw=nw0,nw1
                     psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                        iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,iw) &
                        +reshape(source=recvbuffer(ibuf:ibufnext-1),&
                         shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        ixmin3:ixmax3,iw))) /ps(igrid)%dvolume(ixmin1:ixmax1,&
                        ixmin2:ixmax2,ixmin3:ixmax3)
                     ibuf=ibuf+nbuf
                   end do
                   ibuf=ibufnext
                 end if
               end if
            end do
           end do
           end do
           end do
         end select
       end do
     end do

     if (nsend>0) then
       call MPI_WAITALL(nsend,fc_sendreq,fc_sendstat,ierrmpi)
     end if

   end subroutine fix_conserve

   subroutine store_flux(igrid,fC,idimmin,idimmax,nwfluxin)
     use mod_global_parameters

     integer, intent(in)          :: igrid, idimmin,idimmax, nwfluxin
     double precision, intent(in) :: fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
        ixGlo3:ixGhi3,1:nwfluxin,1:ndim)

     integer :: idims, iside, i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3, ix1,ix2,&
        ix3, ixCo1,ixCo2,ixCo3, nxCo1,nxCo2,nxCo3, iw

     do idims = idimmin,idimmax
       select case (idims)
         case (1)
         do iside=1,2
           i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3)
           i3=kr(3,1)*(2*iside-3);

           if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

           select case (neighbor_type(i1,i2,i3,igrid))
           case (neighbor_fine)
             select case (iside)
             case (1)
               pflux(iside,1,igrid)%flux(1,:,:,1:nwfluxin) = -fC(nghostcells,&
                  ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:nwfluxin,1)
             case (2)
               pflux(iside,1,igrid)%flux(1,:,:,1:nwfluxin) = fC(ixMhi1,&
                  ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:nwfluxin,1)
             end select
           case (neighbor_coarse)
             nxCo1=1;nxCo2=ixGhi2/2-nghostcells;nxCo3=ixGhi3/2-nghostcells;
             select case (iside)
             case (1)
               do iw=1,nwfluxin
                do ixCo3=1,nxCo3
         do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                   ix1=nghostcells;ix2=ixMlo2+2*(ixCo2-1)
                   ix3=ixMlo3+2*(ixCo3-1);
                   pflux(iside,1,igrid)%flux(ixCo1,ixCo2,ixCo3,&
                      iw) = sum(fC(ix1,ix2:ix2+1,ix3:ix3+1,iw,1))
                end do
         end do
         end do
               end do
             case (2)
               do iw=1,nwfluxin
                do ixCo3=1,nxCo3
         do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                   ix1=ixMhi1;ix2=ixMlo2+2*(ixCo2-1);ix3=ixMlo3+2*(ixCo3-1);
                   pflux(iside,1,igrid)%flux(ixCo1,ixCo2,ixCo3,&
                      iw) =-sum(fC(ix1,ix2:ix2+1,ix3:ix3+1,iw,1))
                end do
         end do
         end do
               end do
             end select
           end select
         end do
         case (2)
         do iside=1,2
           i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3)
           i3=kr(3,2)*(2*iside-3);

           if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

           select case (neighbor_type(i1,i2,i3,igrid))
           case (neighbor_fine)
             select case (iside)
             case (1)
               pflux(iside,2,igrid)%flux(:,1,:,1:nwfluxin) = -fC(ixMlo1:ixMhi1,&
                  nghostcells,ixMlo3:ixMhi3,1:nwfluxin,2)
             case (2)
               pflux(iside,2,igrid)%flux(:,1,:,1:nwfluxin) = fC(ixMlo1:ixMhi1,&
                  ixMhi2,ixMlo3:ixMhi3,1:nwfluxin,2)
             end select
           case (neighbor_coarse)
             nxCo1=ixGhi1/2-nghostcells;nxCo2=1;nxCo3=ixGhi3/2-nghostcells;
             select case (iside)
             case (1)
               do iw=1,nwfluxin
                do ixCo3=1,nxCo3
         do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                   ix1=ixMlo1+2*(ixCo1-1);ix2=nghostcells
                   ix3=ixMlo3+2*(ixCo3-1);
                   pflux(iside,2,igrid)%flux(ixCo1,ixCo2,ixCo3,&
                      iw) = sum(fC(ix1:ix1+1,ix2,ix3:ix3+1,iw,2))
                end do
         end do
         end do
               end do
             case (2)
               do iw=1,nwfluxin
                do ixCo3=1,nxCo3
         do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                   ix1=ixMlo1+2*(ixCo1-1);ix2=ixMhi2;ix3=ixMlo3+2*(ixCo3-1);
                   pflux(iside,2,igrid)%flux(ixCo1,ixCo2,ixCo3,&
                      iw) =-sum(fC(ix1:ix1+1,ix2,ix3:ix3+1,iw,2))
                end do
         end do
         end do
               end do
             end select
           end select
         end do
         case (3)
         do iside=1,2
           i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3)
           i3=kr(3,3)*(2*iside-3);

           if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

           select case (neighbor_type(i1,i2,i3,igrid))
           case (neighbor_fine)
             select case (iside)
             case (1)
               pflux(iside,3,igrid)%flux(:,:,1,1:nwfluxin) = -fC(ixMlo1:ixMhi1,&
                  ixMlo2:ixMhi2,nghostcells,1:nwfluxin,3)
             case (2)
               pflux(iside,3,igrid)%flux(:,:,1,1:nwfluxin) = fC(ixMlo1:ixMhi1,&
                  ixMlo2:ixMhi2,ixMhi3,1:nwfluxin,3)
             end select
           case (neighbor_coarse)
             nxCo1=ixGhi1/2-nghostcells;nxCo2=ixGhi2/2-nghostcells;nxCo3=1;
             select case (iside)
             case (1)
               do iw=1,nwfluxin
                do ixCo3=1,nxCo3
         do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                   ix1=ixMlo1+2*(ixCo1-1);ix2=ixMlo2+2*(ixCo2-1)
                   ix3=nghostcells;
                   pflux(iside,3,igrid)%flux(ixCo1,ixCo2,ixCo3,&
                      iw) = sum(fC(ix1:ix1+1,ix2:ix2+1,ix3,iw,3))
                end do
         end do
         end do
               end do
             case (2)
               do iw=1,nwfluxin
                do ixCo3=1,nxCo3
         do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                   ix1=ixMlo1+2*(ixCo1-1);ix2=ixMlo2+2*(ixCo2-1);ix3=ixMhi3;
                   pflux(iside,3,igrid)%flux(ixCo1,ixCo2,ixCo3,&
                      iw) =-sum(fC(ix1:ix1+1,ix2:ix2+1,ix3,iw,3))
                end do
         end do
         end do
               end do
             end select
           end select
         end do
       end select
     end do

   end subroutine store_flux

   subroutine store_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      fE,idimmin,idimmax)
     use mod_global_parameters
     
     integer, intent(in)          :: igrid, ixImin1,ixImin2,ixImin3,ixImax1,&
        ixImax2,ixImax3, idimmin,idimmax
     double precision, intent(in) :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
        ixImin3:ixImax3,sdim:3)
     
     integer :: idims, idir, iside, i1,i2,i3
     integer :: pi1,pi2,pi3, mi1,mi2,mi3, ph1,ph2,ph3, mh1,mh2,mh3 !To detect corners
     integer :: ixMcmin1,ixMcmin2,ixMcmin3,ixMcmax1,ixMcmax2,ixMcmax3

     do idims = idimmin,idimmax  !loop over face directions
       !! Loop over block faces
       do iside=1,2 
         i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
         i3=kr(3,idims)*(2*iside-3);
         if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
         select case (neighbor_type(i1,i2,i3,igrid))
         case (neighbor_fine)
           ! The neighbour is finer
           ! Face direction, side (left or right), restrict ==ired?, fE
           call flux_to_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
              ixImax3,idims,iside,.false.,fE)
         case(neighbor_coarse)
           ! The neighbour is coarser
           call flux_to_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
              ixImax3,idims,iside,.true.,fE)
         case(neighbor_sibling)
           ! If the neighbour is at the same level,
           ! check if there are corners
           ! If there is any corner, store the fluxes from that side
           do idir=idims+1,ndim
             pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
             mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
             ph1=pi1-kr(idims,1)*(2*iside-3);ph2=pi2-kr(idims,2)*(2*iside-3)
             ph3=pi3-kr(idims,3)*(2*iside-3);
             mh1=mi1-kr(idims,1)*(2*iside-3);mh2=mi2-kr(idims,2)*(2*iside-3)
             mh3=mi3-kr(idims,3)*(2*iside-3);
             if (neighbor_type(pi1,pi2,pi3,igrid)==4.and.neighbor_type(ph1,ph2,&
                ph3,igrid)==3) then
               call flux_to_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                  ixImax3,idims,iside,.false.,fE)
             end if
             if (neighbor_type(mi1,mi2,mi3,igrid)==4.and.neighbor_type(mh1,mh2,&
                mh3,igrid)==3) then
               call flux_to_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                  ixImax3,idims,iside,.false.,fE)
             end if
           end do
         end select
       end do
     end do

   end subroutine store_edge

   subroutine flux_to_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
      ixImax3,idims,iside,restrict,fE)
     use mod_global_parameters

     integer                      :: igrid,ixImin1,ixImin2,ixImin3,ixImax1,&
        ixImax2,ixImax3,idims,iside
     logical                      :: restrict
     double precision, intent(in) :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
        ixImin3:ixImax3,sdim:3)

     integer                      :: idir1,idir2
     integer                      :: ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,&
        ixEmax3,ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3, jxFmin1,&
        jxFmin2,jxFmin3,jxFmax1,jxFmax2,jxFmax3, nx1,nx2,nx3,nxCo1,nxCo2,nxCo3

     nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
     nxCo1=nx1/2;nxCo2=nx2/2;nxCo3=nx3/2;
     ! ixE are the indices on the 'edge' array.
     ! ixF are the indices on the 'fE' array
     ! jxF are indices advanced to perform the flux restriction (sum) in 3D
     ! A line integral of the electric field on the coarse side
     ! lies over two edges on the fine side. So, in 3D we restrict by summing
     ! over two cells on the fine side.

     do idir1=1,ndim-1
       ! 3D: rotate indices among 1 and 2 to save space 
       idir2=mod(idir1+idims-1,3)+1
      

       if (restrict) then
         ! Set up indices for restriction
         ixFmin1=ixMlo1-1+kr(1,idir2);ixFmin2=ixMlo2-1+kr(2,idir2)
         ixFmin3=ixMlo3-1+kr(3,idir2);
         ixFmax1=ixMhi1-kr(1,idir2);ixFmax2=ixMhi2-kr(2,idir2)
         ixFmax3=ixMhi3-kr(3,idir2);
         
         jxFmin1=ixFmin1+kr(1,idir2);jxFmin2=ixFmin2+kr(2,idir2)
         jxFmin3=ixFmin3+kr(3,idir2);jxFmax1=ixFmax1+kr(1,idir2)
         jxFmax2=ixFmax2+kr(2,idir2);jxFmax3=ixFmax3+kr(3,idir2);

         ixEmin1=0+kr(1,idir2);ixEmin2=0+kr(2,idir2);ixEmin3=0+kr(3,idir2);
         ixEmax1=nxCo1;ixEmax2=nxCo2;ixEmax3=nxCo3;
         select case(idims)
        case(1)
           ixEmin1=1;ixEmax1=1;
           select case(iside)
           case(1)
             ixFmax1=ixFmin1
             jxFmax1=ixFmin1
           case(2)
             ixFmin1=ixFmax1
             jxFmin1=ixFmax1
           end select
        
        case(2)
           ixEmin2=1;ixEmax2=1;
           select case(iside)
           case(1)
             ixFmax2=ixFmin2
             jxFmax2=ixFmin2
           case(2)
             ixFmin2=ixFmax2
             jxFmin2=ixFmax2
           end select
        
        case(3)
           ixEmin3=1;ixEmax3=1;
           select case(iside)
           case(1)
             ixFmax3=ixFmin3
             jxFmax3=ixFmin3
           case(2)
             ixFmin3=ixFmax3
             jxFmin3=ixFmax3
           end select
        
         end select

       pflux(iside,idims,igrid)%edge(ixEmin1:ixEmax1,ixEmin2:ixEmax2,&
          ixEmin3:ixEmax3,idir1)=fE(ixFmin1:ixFmax1:2,ixFmin2:ixFmax2:2,&
          ixFmin3:ixFmax3:2,idir2) +fE(jxFmin1:jxFmax1:2,jxFmin2:jxFmax2:2,&
          jxFmin3:jxFmax3:2,idir2);

       else
         ! Set up indices for copying 
         ixFmin1=ixMlo1-1+kr(1,idir2);ixFmin2=ixMlo2-1+kr(2,idir2)
         ixFmin3=ixMlo3-1+kr(3,idir2);
         ixFmax1=ixMhi1;ixFmax2=ixMhi2;ixFmax3=ixMhi3;
         ixEmin1=0+kr(1,idir2);ixEmin2=0+kr(2,idir2);ixEmin3=0+kr(3,idir2);
         ixEmax1=nx1;ixEmax2=nx2;ixEmax3=nx3;

         select case(idims)
        case(1)
           ixEmin1=1;ixEmax1=1;
           select case(iside)
           case(1)
             ixFmax1=ixFmin1
           case(2)
             ixFmin1=ixFmax1
           end select
        
        case(2)
           ixEmin2=1;ixEmax2=1;
           select case(iside)
           case(1)
             ixFmax2=ixFmin2
           case(2)
             ixFmin2=ixFmax2
           end select
        
        case(3)
           ixEmin3=1;ixEmax3=1;
           select case(iside)
           case(1)
             ixFmax3=ixFmin3
           case(2)
             ixFmin3=ixFmax3
           end select
        
         end select

         pflux(iside,idims,igrid)%edge(ixEmin1:ixEmax1,ixEmin2:ixEmax2,&
            ixEmin3:ixEmax3,idir1)=fE(ixFmin1:ixFmax1,ixFmin2:ixFmax2,&
            ixFmin3:ixFmax3,idir2)

       end if

     end do

   end subroutine flux_to_edge

   subroutine fix_edges(psuse,idimmin,idimmax)
     use mod_global_parameters

     type(state) :: psuse(max_blocks)
     integer, intent(in) :: idimmin,idimmax

     integer :: iigrid, igrid, idims, iside, iotherside, i1,i2,i3, ic1,ic2,ic3,&
         inc1,inc2,inc3, ixMcmin1,ixMcmin2,ixMcmin3,ixMcmax1,ixMcmax2,ixMcmax3
     integer :: nbuf, ibufnext
     integer :: ibufnext_cc
     integer :: pi1,pi2,pi3, mi1,mi2,mi3, ph1,ph2,ph3, mh1,mh2,mh3 !To detect corners
     integer :: ixEmin1(1:3),ixEmin2(1:3),ixEmin3(1:3),ixEmax1(1:3),&
        ixEmax2(1:3),ixEmax3(1:3), ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,&
        ixtEmax2,ixtEmax3, ixFmin1(1:ndim),ixFmin2(1:ndim),ixFmin3(1:ndim),&
        ixFmax1(1:ndim),ixFmax2(1:ndim),ixFmax3(1:ndim), ixfEmin1(1:3),&
        ixfEmin2(1:3),ixfEmin3(1:3),ixfEmax1(1:3),ixfEmax2(1:3),ixfEmax3(1:3)
     integer :: nx1,nx2,nx3, idir, ix, ipe_neighbor, ineighbor
     logical :: pcorner(1:ndim),mcorner(1:ndim)

     if (nrecv_ct>0) then
        call MPI_WAITALL(nrecv_ct,cc_recvreq,cc_recvstat,ierrmpi)
     end if

     ! Initialise buffer counter again
     ibuf=1
     ibuf_cc=1
     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims= idimmin,idimmax
         do iside=1,2
           i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
           i3=kr(3,idims)*(2*iside-3);
           if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
           select case(neighbor_type(i1,i2,i3,igrid))
           case(neighbor_fine)
             ! The first neighbour is finer
             if (.not.neighbor_active(i1,i2,i3,&
                igrid).or..not.neighbor_active(0,0,0,igrid) ) then
               do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
                  inc3=2*i3+ic3
               do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
                  inc2=2*i2+ic2
               do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
                  inc1=2*i1+ic1
                  ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                  !! When the neighbour is in a different process
                  if (ipe_neighbor/=mype) then
                     ibufnext=ibuf+isize(idims)
                     ibuf=ibufnext
                     end if
               end do
               end do
               end do
                cycle
             end if

             ! Check if there are corners
             pcorner=.false.
             mcorner=.false.
             do idir=1,ndim
               pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
               mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
               ph1=pi1-kr(idims,1)*(2*iside-3)
               ph2=pi2-kr(idims,2)*(2*iside-3)
               ph3=pi3-kr(idims,3)*(2*iside-3);
               mh1=mi1-kr(idims,1)*(2*iside-3)
               mh2=mi2-kr(idims,2)*(2*iside-3)
               mh3=mi3-kr(idims,3)*(2*iside-3);
               if (neighbor_type(ph1,ph2,ph3,&
                  igrid)==neighbor_fine) pcorner(idir)=.true.
               if (neighbor_type(mh1,mh2,mh3,&
                  igrid)==neighbor_fine) mcorner(idir)=.true.
             end do
             ! Calculate indices range
             call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
                ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
                ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
                ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,igrid,idims,iside,.false.,&
                .false.,0,0,0,pcorner,mcorner)
             ! Remove coarse part of circulation
             call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
                ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
                ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
                ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,pflux(iside,idims,&
                igrid)%edge,idims,iside,.false.,psuse(igrid))
             ! Add fine part of the circulation
            do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
            do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
            do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               iotherside=3-iside
               nx1=(ixMhi1-ixMlo1+1)/2;nx2=(ixMhi2-ixMlo2+1)/2
               nx3=(ixMhi3-ixMlo3+1)/2;
               call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                  ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                  ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
                  ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,igrid,&
                  idims,iside,.true.,.false.,inc1,inc2,inc3,pcorner,mcorner)
               if (ipe_neighbor==mype) then
                 call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                    ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                    ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
                    ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
                    pflux(iotherside,idims,ineighbor)%edge,idims,iside,.true.,&
                    psuse(igrid))
               else
                 ibufnext=ibuf+isize(idims)
                 call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                    ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                    ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
                    ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
                    reshape(source=recvbuffer(ibufnext-&
                    isize_stg(idims):ibufnext-1),shape=(/ ixtEmax1-ixtEmin1+1,&
                    ixtEmax2-ixtEmin2+1,ixtEmax3-ixtEmin3+1 ,3-1 /)),idims,&
                    iside,.true.,psuse(igrid))
                 ibuf=ibufnext
               end if
            end do
            end do
            end do

           case(neighbor_sibling)
             ! The first neighbour is at the same level
             ! Check if there are corners
             do idir=idims+1,ndim
               pcorner=.false.
               mcorner=.false.
               pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
               mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
               ph1=pi1-kr(idims,1)*(2*iside-3)
               ph2=pi2-kr(idims,2)*(2*iside-3)
               ph3=pi3-kr(idims,3)*(2*iside-3);
               mh1=mi1-kr(idims,1)*(2*iside-3)
               mh2=mi2-kr(idims,2)*(2*iside-3)
               mh3=mi3-kr(idims,3)*(2*iside-3);
               if (neighbor_type(pi1,pi2,pi3,&
                  igrid)==neighbor_fine.and.neighbor_type(ph1,ph2,ph3,&
                  igrid)==neighbor_sibling.and.neighbor_pole(pi1,pi2,pi3,&
                  igrid)==0) then
                 pcorner(idir)=.true.
                 ! Remove coarse part
                 ! Set indices
                 call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                    ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                    ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
                    ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
                    igrid,idims,iside,.false.,.true.,0,0,0,pcorner,mcorner)
                 ! Remove
                 call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                    ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                    ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
                    ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
                    pflux(iside,idims,igrid)%edge,idims,iside,.false.,&
                    psuse(igrid))
                 ! Add fine part
                 ! Find relative position of finer grid
       do ix=1,2
                 inc1=kr(idims,1)*3*(iside-1)+3*kr(idir,1)+kr(6-idir-idims,&
                    1)*ix
                 inc2=kr(idims,2)*3*(iside-1)+3*kr(idir,2)+kr(6-idir-idims,&
                    2)*ix
                 inc3=kr(idims,3)*3*(iside-1)+3*kr(idir,3)+kr(6-idir-idims,&
                    3)*ix;
                 ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
                 ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                 iotherside=3-iside
                 ! Set indices
                 call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                    ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                    ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
                    ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
                    igrid,idims,iside,.true.,.true.,inc1,inc2,inc3,pcorner,&
                    mcorner)
                 ! add
                 if (ipe_neighbor==mype) then
                   call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                      ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                      ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
                      ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
                      pflux(iotherside,idims,ineighbor)%edge,idims,iside,&
                      .true.,psuse(igrid))
                 else
                   ibufnext_cc=ibuf_cc+isize_stg(idims)
                   call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                      ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                      ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
                      ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
                      reshape(source=recvbuffer_cc(ibuf_cc:ibufnext_cc-1),&
                      shape=(/ ixtEmax1-ixtEmin1+1,ixtEmax2-ixtEmin2+1,&
                      ixtEmax3-ixtEmin3+1 ,3-1 /)),idims,iside,.true.,&
                      psuse(igrid))
                   ibuf_cc=ibufnext_cc
                 end if
       end do
               ! Set CoCorner to false again for next step
                 pcorner(idir)=.false.
               end if

               if (neighbor_type(mi1,mi2,mi3,&
                  igrid)==neighbor_fine.and.neighbor_type(mh1,mh2,mh3,&
                  igrid)==neighbor_sibling.and.neighbor_pole(mi1,mi2,mi3,&
                  igrid)==0) then
                   mcorner(idir)=.true.
                   ! Remove coarse part
                   ! Set indices
                   call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                      ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                      ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
                      ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
                      igrid,idims,iside,.false.,.true.,0,0,0,pcorner,mcorner)
                   ! Remove
                   call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                      ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                      ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
                      ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
                      pflux(iside,idims,igrid)%edge,idims,iside,.false.,&
                      psuse(igrid))
                   ! Add fine part
                   ! Find relative position of finer grid
         do ix=1,2
                   inc1=kr(idims,1)*3*(iside-1)+kr(6-idir-idims,1)*ix
                   inc2=kr(idims,2)*3*(iside-1)+kr(6-idir-idims,2)*ix
                   inc3=kr(idims,3)*3*(iside-1)+kr(6-idir-idims,3)*ix;
                   ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
                   ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                   iotherside=3-iside
                   ! Set indices
                   call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                      ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                      ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
                      ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
                      igrid,idims,iside,.true.,.true.,inc1,inc2,inc3,pcorner,&
                      mcorner)
                   ! add
                   if (ipe_neighbor==mype) then
                     call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                        ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                        ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,&
                        ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,&
                        ixfEmax3,pflux(iotherside,idims,ineighbor)%edge,idims,&
                        iside,.true.,psuse(igrid))
                   else
                     ibufnext_cc=ibuf_cc+isize_stg(idims)
                     call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
                        ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
                        ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,&
                        ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,&
                        ixfEmax3,reshape(source=recvbuffer_cc(&
                        ibuf_cc:ibufnext_cc-1),shape=(/ ixtEmax1-ixtEmin1+1,&
                        ixtEmax2-ixtEmin2+1,ixtEmax3-ixtEmin3+1 ,3-1 /)),idims,&
                        iside,.true.,psuse(igrid))
                     ibuf_cc=ibufnext_cc
                   end if
         end do
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
   subroutine set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
      ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,ixEmin2,&
      ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,&
      ixfEmax2,ixfEmax3,igrid,idims,iside,add,CoCorner,inc1,inc2,inc3,pcorner,&
      mcorner)
     use mod_global_parameters
     
     integer,intent(in)    :: igrid,idims,iside,inc1,inc2,inc3
     logical,intent(in)    :: add,CoCorner
     logical,intent(inout) :: pcorner(1:ndim),mcorner(1:ndim)
     integer,intent(out)   :: ixFmin1(1:ndim),ixFmin2(1:ndim),ixFmin3(1:ndim),&
        ixFmax1(1:ndim),ixFmax2(1:ndim),ixFmax3(1:ndim),ixtEmin1,ixtEmin2,&
        ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1(1:3),ixEmin2(1:3),&
        ixEmin3(1:3),ixEmax1(1:3),ixEmax2(1:3),ixEmax3(1:3),ixfEmin1(1:3),&
        ixfEmin2(1:3),ixfEmin3(1:3),ixfEmax1(1:3),ixfEmax2(1:3),ixfEmax3(1:3) !Indices for faces and edges
     integer               :: icor1,icor2,icor3,idim1,idir,nx1,nx2,nx3,middle1,&
        middle2,middle3
     integer               :: ixtfEmin1,ixtfEmin2,ixtfEmin3,ixtfEmax1,&
        ixtfEmax2,ixtfEmax3

     ! ixF -> Indices for the _F_aces, and
     ! depends on the field component
     ! ixtE -> are the _t_otal range of the 'edge' array
     ! ixE -> are the ranges of the edge array,
     ! depending on the component
     ! ixfE -> are the ranges of the fE array (3D),
     ! and also depend on the component

     ! ... General ...
     ! Assign indices for the size of the E field array

     ixtfEmin1=ixMlo1-1;ixtfEmin2=ixMlo2-1;ixtfEmin3=ixMlo3-1;
     ixtfEmax1=ixMhi1;ixtfEmax2=ixMhi2;ixtfEmax3=ixMhi3;

     if(add) then
       nx1=(ixMhi1-ixMlo1+1)/2;nx2=(ixMhi2-ixMlo2+1)/2
       nx3=(ixMhi3-ixMlo3+1)/2;
     else
       nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
     end if

     do idim1=1,ndim
       ixtEmin1=0;ixtEmin2=0;ixtEmin3=0;
       ixtEmax1=nx1;ixtEmax2=nx2;ixtEmax3=nx3;
       select case(idims)
       case(1)
         ixtEmin1=1;ixtEmax1=1;
         if (iside==1) ixtfEmax1=ixtfEmin1;
         if (iside==2) ixtfEmin1=ixtfEmax1;
       
       case(2)
         ixtEmin2=1;ixtEmax2=1;
         if (iside==1) ixtfEmax2=ixtfEmin2;
         if (iside==2) ixtfEmin2=ixtfEmax2;
       
       case(3)
         ixtEmin3=1;ixtEmax3=1;
         if (iside==1) ixtfEmax3=ixtfEmin3;
         if (iside==2) ixtfEmin3=ixtfEmax3;
       
       end select
     end do

     ! Assign indices, considering only the face
     ! (idims and iside)
     do idim1=1,ndim
       ixFmin1(idim1)=ixMlo1-kr(idim1,1);ixFmin2(idim1)=ixMlo2-kr(idim1,2)
       ixFmin3(idim1)=ixMlo3-kr(idim1,3);
       ixFmax1(idim1)=ixMhi1;ixFmax2(idim1)=ixMhi2;ixFmax3(idim1)=ixMhi3;
       select case(idims)
       case(1)
          select case(iside)
          case(1)
          ixFmax1(idim1)=ixFmin1(idim1)
          case(2)
          ixFmin1(idim1)=ixFmax1(idim1)
          end select
       
       case(2)
          select case(iside)
          case(1)
          ixFmax2(idim1)=ixFmin2(idim1)
          case(2)
          ixFmin2(idim1)=ixFmax2(idim1)
          end select
       
       case(3)
          select case(iside)
          case(1)
          ixFmax3(idim1)=ixFmin3(idim1)
          case(2)
          ixFmin3(idim1)=ixFmax3(idim1)
          end select
       
       end select
     end do
     ! ... Relative position ...
     ! Restrict range using relative position
     if(add) then
       middle1=(ixMhi1+ixMlo1)/2;middle2=(ixMhi2+ixMlo2)/2
       middle3=(ixMhi3+ixMlo3)/2;
       
       if(inc1==1) then
         ixFmax1(:)=middle1
         ixtfEmax1=middle1
       end if
       if(inc1==2) then
         ixFmin1(:)=middle1+1
         ixtfEmin1=middle1
       end if
       
       
       if(inc2==1) then
         ixFmax2(:)=middle2
         ixtfEmax2=middle2
       end if
       if(inc2==2) then
         ixFmin2(:)=middle2+1
         ixtfEmin2=middle2
       end if
       
       
       if(inc3==1) then
         ixFmax3(:)=middle3
         ixtfEmax3=middle3
       end if
       if(inc3==2) then
         ixFmin3(:)=middle3+1
         ixtfEmin3=middle3
       end if
       
     end if
     ! ... Adjust ranges of edges according to direction ...
     do idim1=1,3
       ixfEmax1(idim1)=ixtfEmax1;ixfEmax2(idim1)=ixtfEmax2
       ixfEmax3(idim1)=ixtfEmax3;
       ixEmax1(idim1)=ixtEmax1;ixEmax2(idim1)=ixtEmax2
       ixEmax3(idim1)=ixtEmax3;
       ixfEmin1(idim1)=ixtfEmin1+kr(idim1,1)
       ixfEmin2(idim1)=ixtfEmin2+kr(idim1,2)
       ixfEmin3(idim1)=ixtfEmin3+kr(idim1,3);
       ixEmin1(idim1)=ixtEmin1+kr(idim1,1)
       ixEmin2(idim1)=ixtEmin2+kr(idim1,2)
       ixEmin3(idim1)=ixtEmin3+kr(idim1,3);
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
              if (1==idim1) then
                 ixfEmin1(idir)=ixfEmax1(idir)
                 if (add) then
                   ixEmax1(idir) =ixEmin1(idir)
                 else
                   ixEmin1(idir) =ixEmax1(idir)
                 end if
               end if
              if (2==idim1) then
                 ixfEmin2(idir)=ixfEmax2(idir)
                 if (add) then
                   ixEmax2(idir) =ixEmin2(idir)
                 else
                   ixEmin2(idir) =ixEmax2(idir)
                 end if
               end if
              if (3==idim1) then
                 ixfEmin3(idir)=ixfEmax3(idir)
                 if (add) then
                   ixEmax3(idir) =ixEmin3(idir)
                 else
                   ixEmin3(idir) =ixEmax3(idir)
                 end if
               end if
             else
               ixEmin1(idir)=1;ixEmin2(idir)=1;ixEmin3(idir)=1;
               ixEmax1(idir)=0;ixEmax2(idir)=0;ixEmax3(idir)=0;
               ixfEmin1(idir)=1;ixfEmin2(idir)=1;ixfEmin3(idir)=1;
               ixfEmax1(idir)=0;ixfEmax2(idir)=0;ixfEmax3(idir)=0;
             end if
           end do
         end if
         if (mcorner(idim1)) then
           do idir=1,3
             if (idir==6-idim1-idims) then
              if (1==idim1) then
                 ixfEmax1(idir)=ixfEmin1(idir)
                 if (add) then
                   ixEmin1(idir) =ixEmax1(idir)
                 else
                   ixEmax1(idir) =ixEmin1(idir)
                 end if
               end if
              if (2==idim1) then
                 ixfEmax2(idir)=ixfEmin2(idir)
                 if (add) then
                   ixEmin2(idir) =ixEmax2(idir)
                 else
                   ixEmax2(idir) =ixEmin2(idir)
                 end if
               end if
              if (3==idim1) then
                 ixfEmax3(idir)=ixfEmin3(idir)
                 if (add) then
                   ixEmin3(idir) =ixEmax3(idir)
                 else
                   ixEmax3(idir) =ixEmin3(idir)
                 end if
               end if
             else
               ixEmin1(idir)=1;ixEmin2(idir)=1;ixEmin3(idir)=1;
               ixEmax1(idir)=0;ixEmax2(idir)=0;ixEmax3(idir)=0;
               ixfEmin1(idir)=1;ixfEmin2(idir)=1;ixfEmin3(idir)=1;
               ixfEmax1(idir)=0;ixfEmax2(idir)=0;ixfEmax3(idir)=0;
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
      if((idims.gt.1).and.pcorner(1)) then
         if((.not.add).or.(inc1==2)) then
 !ixFmax1(:)=ixFmax1(:)-kr(1,1);!ixFmax2(:)=ixFmax2(:)-kr(1,2);!ixFmax3(:)=ixFmax3(:)-kr(1,3);
           do idir=1,3
             if ((idir==idims).or.(idir==1)) cycle
               ixfEmax1(idir)=ixfEmax1(idir)-1
               ixEmax1(idir)=ixEmax1(idir)-1
           end do
         end if
       end if
      if((idims.gt.2).and.pcorner(2)) then
         if((.not.add).or.(inc2==2)) then
 !ixFmax1(:)=ixFmax1(:)-kr(2,1);!ixFmax2(:)=ixFmax2(:)-kr(2,2);!ixFmax3(:)=ixFmax3(:)-kr(2,3);
           do idir=1,3
             if ((idir==idims).or.(idir==2)) cycle
               ixfEmax2(idir)=ixfEmax2(idir)-1
               ixEmax2(idir)=ixEmax2(idir)-1
           end do
         end if
       end if
      if((idims.gt.3).and.pcorner(3)) then
         if((.not.add).or.(inc3==2)) then
 !ixFmax1(:)=ixFmax1(:)-kr(3,1);!ixFmax2(:)=ixFmax2(:)-kr(3,2);!ixFmax3(:)=ixFmax3(:)-kr(3,3);
           do idir=1,3
             if ((idir==idims).or.(idir==3)) cycle
               ixfEmax3(idir)=ixfEmax3(idir)-1
               ixEmax3(idir)=ixEmax3(idir)-1
           end do
         end if
       end if
      if((idims>1).and.mcorner(1)) then
         if((.not.add).or.(inc1==1)) then
 !ixFmin1(:)=ixFmin1(:)+kr(1,1);!ixFmin2(:)=ixFmin2(:)+kr(1,2);!ixFmin3(:)=ixFmin3(:)+kr(1,3);
           do idir=1,3
             if ((idir==idims).or.(idir==1)) cycle
               ixfEmin1(idir)=ixfEmin1(idir)+1
               ixEmin1(idir)=ixEmin1(idir)+1
           end do
         end if
       end if
      if((idims>2).and.mcorner(2)) then
         if((.not.add).or.(inc2==1)) then
 !ixFmin1(:)=ixFmin1(:)+kr(2,1);!ixFmin2(:)=ixFmin2(:)+kr(2,2);!ixFmin3(:)=ixFmin3(:)+kr(2,3);
           do idir=1,3
             if ((idir==idims).or.(idir==2)) cycle
               ixfEmin2(idir)=ixfEmin2(idir)+1
               ixEmin2(idir)=ixEmin2(idir)+1
           end do
         end if
       end if
      if((idims>3).and.mcorner(3)) then
         if((.not.add).or.(inc3==1)) then
 !ixFmin1(:)=ixFmin1(:)+kr(3,1);!ixFmin2(:)=ixFmin2(:)+kr(3,2);!ixFmin3(:)=ixFmin3(:)+kr(3,3);
           do idir=1,3
             if ((idir==idims).or.(idir==3)) cycle
               ixfEmin3(idir)=ixfEmin3(idir)+1
               ixEmin3(idir)=ixEmin3(idir)+1
           end do
         end if
       end if
     end if

   end subroutine set_ix_circ

   subroutine add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
      ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,ixEmin2,&
      ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,&
      ixfEmax2,ixfEmax3,edge,idims,iside,add,s)
     use mod_global_parameters

     type(state)        :: s
     integer,intent(in) :: idims,iside
     integer            :: ixFmin1(1:ndim),ixFmin2(1:ndim),ixFmin3(1:ndim),&
        ixFmax1(1:ndim),ixFmax2(1:ndim),ixFmax3(1:ndim),ixtEmin1,ixtEmin2,&
        ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1(1:3),ixEmin2(1:3),&
        ixEmin3(1:3),ixEmax1(1:3),ixEmax2(1:3),ixEmax3(1:3),ixfEmin1(1:3),&
        ixfEmin2(1:3),ixfEmin3(1:3),ixfEmax1(1:3),ixfEmax2(1:3),ixfEmax3(1:3)
     double precision   :: edge(ixtEmin1:ixtEmax1,ixtEmin2:ixtEmax2,&
        ixtEmin3:ixtEmax3,1:ndim-1)
     logical,intent(in) :: add

     integer            :: idim1,idim2,idir,middle1,middle2,middle3
     integer            :: ixfECmin1,ixfECmin2,ixfECmin3,ixfECmax1,ixfECmax2,&
        ixfECmax3,ixECmin1,ixECmin2,ixECmin3,ixECmax1,ixECmax2,ixECmax3
     double precision   :: fE(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
        sdim:3)
     double precision   :: circ(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
        1:ndim)
     integer            :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,hxmin1,&
        hxmin2,hxmin3,hxmax1,hxmax2,hxmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
        ixCmax2,ixCmax3,hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,hxCmax3 !Indices for edges

     ! ixF -> Indices for the faces, depends on the field component
     ! ixE -> Total range for the edges
     ! ixfE -> Edges in fE (3D) array
     ! ix,hx,ixC,hxC -> Auxiliary indices
     ! Assign quantities stored ad edges to make it as similar as 
     ! possible to the routine updatefaces.
     fE(:,:,:,:)=zero
     do idim1=1,ndim-1
        ! 3D: rotate indices (see routine flux_to_edge)
       idir=mod(idim1+idims-1,3)+1
       
       ixfECmin1=ixfEmin1(idir);ixfECmin2=ixfEmin2(idir)
       ixfECmin3=ixfEmin3(idir);ixfECmax1=ixfEmax1(idir)
       ixfECmax2=ixfEmax2(idir);ixfECmax3=ixfEmax3(idir);
       ixECmin1=ixEmin1(idir);ixECmin2=ixEmin2(idir);ixECmin3=ixEmin3(idir)
       ixECmax1=ixEmax1(idir);ixECmax2=ixEmax2(idir);ixECmax3=ixEmax3(idir);
       fE(ixfECmin1:ixfECmax1,ixfECmin2:ixfECmax2,ixfECmin3:ixfECmax3,&
          idir)=edge(ixECmin1:ixECmax1,ixECmin2:ixECmax2,ixECmin3:ixECmax3,&
          idim1)
     end do

     ! Calculate part of circulation needed
     circ=zero
     do idim1=1,ndim
        do idim2=1,ndim
           do idir=sdim,3
             if (lvc(idim1,idim2,idir)==0) cycle
             ! Assemble indices
             ixCmin1=ixFmin1(idim1);ixCmin2=ixFmin2(idim1)
             ixCmin3=ixFmin3(idim1);ixCmax1=ixFmax1(idim1)
             ixCmax2=ixFmax2(idim1);ixCmax3=ixFmax3(idim1);
             hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
             hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
             hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
             if(idim1==idims) then
               circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)-fE(hxCmin1:hxCmax1,&
                  hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir))
             else
               select case(iside)
               case(2)
                 circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                    idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3,idim1)+lvc(idim1,idim2,&
                    idir)*fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                    idir)
               case(1)
                 circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                    idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3,idim1)-lvc(idim1,idim2,&
                    idir)*fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
                    idir)
               end select
             end if
           end do
        end do
     end do

     ! Divide circulation by surface and add
     do idim1=1,ndim
        ixCmin1=ixFmin1(idim1);ixCmin2=ixFmin2(idim1);ixCmin3=ixFmin3(idim1)
        ixCmax1=ixFmax1(idim1);ixCmax2=ixFmax2(idim1);ixCmax3=ixFmax3(idim1);
        where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)>1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3))
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)
        elsewhere
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
        end where
        ! Add/subtract to field at face
        if (add) then
          s%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)=s%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)
        else
          s%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)=s%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)+circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)
        end if
     end do

   end subroutine add_sub_circ

end module mod_fix_conserve
