module mod_load_balance
  implicit none

contains
  !> reallocate blocks into processors for load balance
  subroutine load_balance
    use mod_forest
    use mod_global_parameters
    use mod_space_filling_curve

    integer :: Morton_no, recv_igrid, recv_ipe, send_igrid, send_ipe, igrid, ipe
    !> MPI recv send variables for AMR
    integer :: itag, irecv, isend
    integer, dimension(:), allocatable :: recvrequest, sendrequest
    integer, dimension(:,:), allocatable :: recvstatus, sendstatus
    !> MPI recv send variables for staggered-variable AMR
    integer :: itag_stg
    integer, dimension(:), allocatable :: recvrequest_stg, sendrequest_stg
    integer, dimension(:,:), allocatable :: recvstatus_stg, sendstatus_stg

    integer, external :: getnode

    ! Jannis: for now, not using version for passive/active blocks
    call get_Morton_range()

    if (npe==1) then
       sfc_to_igrid(:)=sfc(1,Morton_start(mype):Morton_stop(mype))
       return
    end if

    irecv=0
    isend=0
    allocate(recvstatus(MPI_STATUS_SIZE,max_blocks),recvrequest(max_blocks), &
             sendstatus(MPI_STATUS_SIZE,max_blocks),sendrequest(max_blocks))
    recvrequest=MPI_REQUEST_NULL
    sendrequest=MPI_REQUEST_NULL

    if(stagger_grid) then
      allocate(recvstatus_stg(MPI_STATUS_SIZE,max_blocks*^ND),recvrequest_stg(max_blocks*^ND), &
             sendstatus_stg(MPI_STATUS_SIZE,max_blocks*^ND),sendrequest_stg(max_blocks*^ND))
      recvrequest_stg=MPI_REQUEST_NULL
      sendrequest_stg=MPI_REQUEST_NULL
    end if

    do ipe=0,npe-1; do Morton_no=Morton_start(ipe),Morton_stop(ipe)
       recv_ipe=ipe

       send_igrid=sfc(1,Morton_no)
       send_ipe=sfc(2,Morton_no)

       if (recv_ipe/=send_ipe) then
          ! get an igrid number for the new node in recv_ipe processor
          recv_igrid=getnode(recv_ipe)
          ! update node igrid and ipe on the tree
          call change_ipe_tree_leaf(recv_igrid,recv_ipe,send_igrid,send_ipe)
          ! receive physical data of the new node in recv_ipe processor
          if (recv_ipe==mype) call lb_recv
          ! send physical data of the old node in send_ipe processor
          if (send_ipe==mype) call lb_send
       end if
       if (recv_ipe==mype) then
          if (recv_ipe==send_ipe) then
             sfc_to_igrid(Morton_no)=send_igrid
          else
             sfc_to_igrid(Morton_no)=recv_igrid
          end if
       end if
    end do; end do

    if (irecv>0) then
      call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
      if(stagger_grid) call MPI_WAITALL(irecv,recvrequest_stg,recvstatus_stg,ierrmpi)
    end if
    if (isend>0) then
      call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
      if(stagger_grid) call MPI_WAITALL(isend,sendrequest_stg,sendstatus_stg,ierrmpi)
    end if

    deallocate(recvstatus,recvrequest,sendstatus,sendrequest)
    if(stagger_grid) deallocate(recvstatus_stg,recvrequest_stg,sendstatus_stg,sendrequest_stg)

    ! post processing
    do ipe=0,npe-1; do Morton_no=Morton_start(ipe),Morton_stop(ipe)
       recv_ipe=ipe

       send_igrid=sfc(1,Morton_no)
       send_ipe=sfc(2,Morton_no)

       if (recv_ipe/=send_ipe) then
          !if (send_ipe==mype) call dealloc_node(send_igrid)
          call putnode(send_igrid,send_ipe)
       end if
    end do; end do
    {#IFDEF EVOLVINGBOUNDARY
    ! mark physical-boundary blocks on space-filling curve
    do Morton_no=Morton_start(mype),Morton_stop(mype)
       igrid=sfc_to_igrid(Morton_no)
       if (phyboundblock(igrid)) sfc_phybound(Morton_no)=1
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,sfc_phybound,nleafs,MPI_INTEGER,&
                       MPI_SUM,icomm,ierrmpi)
    }

    ! Update sfc array: igrid and ipe info in space filling curve
    call amr_Morton_order()

    contains

      subroutine lb_recv

        call alloc_node(recv_igrid)

        itag=recv_igrid
        irecv=irecv+1
        {#IFDEF EVOLVINGBOUNDARY
        if (phyboundblock(recv_igrid)) then
           call MPI_IRECV(ps(recv_igrid)%w,1,type_block,send_ipe,itag, &
                          icomm,recvrequest(irecv),ierrmpi)
        else
           call MPI_IRECV(ps(recv_igrid)%w,1,type_block_io,send_ipe,itag, &
                          icomm,recvrequest(irecv),ierrmpi)
        end if
        }{#IFNDEF EVOLVINGBOUNDARY
        call MPI_IRECV(ps(recv_igrid)%w,1,type_block_io,send_ipe,itag, &
                       icomm,recvrequest(irecv),ierrmpi)
        }
        if(stagger_grid) then
          itag=recv_igrid+max_blocks
          call MPI_IRECV(ps(recv_igrid)%ws,1,type_block_io_stg,send_ipe,itag, &
               icomm,recvrequest_stg(irecv),ierrmpi)
        end if

      end subroutine lb_recv

      subroutine lb_send

        itag=recv_igrid
        isend=isend+1
        {#IFDEF EVOLVINGBOUNDARY
        if (phyboundblock(send_igrid)) then
           call MPI_ISEND(ps(send_igrid)%w,1,type_block,recv_ipe,itag, &
                          icomm,sendrequest(isend),ierrmpi)
        else
           call MPI_ISEND(ps(send_igrid)%w,1,type_block_io,recv_ipe,itag, &
                          icomm,sendrequest(isend),ierrmpi)
        end if
        }{#IFNDEF EVOLVINGBOUNDARY
        call MPI_ISEND(ps(send_igrid)%w,1,type_block_io,recv_ipe,itag, &
                       icomm,sendrequest(isend),ierrmpi)
        }
        if(stagger_grid) then
          itag=recv_igrid+max_blocks
          call MPI_ISEND(ps(send_igrid)%ws,1,type_block_io_stg,recv_ipe,itag, &
                         icomm,sendrequest_stg(isend),ierrmpi)
        end if

      end subroutine lb_send

  end subroutine load_balance

end module mod_load_balance
