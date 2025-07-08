!> Module to coarsen and refine grids for AMR
module mod_coarsen_refine
  implicit none
  private
  !> MPI recv send variables for AMR
  integer :: itag, irecv, isend
  integer, dimension(:), allocatable :: recvrequest, sendrequest
  integer, dimension(:,:), allocatable :: recvstatus, sendstatus
  !> MPI recv send variables for staggered-variable AMR
  integer :: itag_stg
  integer, dimension(:), allocatable :: recvrequest_stg, sendrequest_stg
  integer, dimension(:,:), allocatable :: recvstatus_stg, sendstatus_stg

  ! Public subroutines
  public :: amr_coarsen_refine

contains

  !> coarsen and refine blocks to update AMR grid
  subroutine amr_coarsen_refine
    use mod_forest
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_usr_methods, only: usr_after_refine
    use mod_amr_fct
    use mod_space_filling_curve
    use mod_load_balance
    use mod_functions_connectivity, only: get_level_range,getigrids,&
       build_connectivity
    use mod_amr_solution_node, only: getnode, putnode
    use mod_functions_forest, only: coarsen_tree_leaf,refine_tree_leaf
    use mod_selectgrids, only: selectgrids
    use mod_refine, only: refine_grids


    
    use mod_multigrid_coupling
   

    integer :: iigrid, igrid, ipe, igridCo, ipeCo, level, ic1,ic2,ic3
    integer, dimension(2,2,2) :: igridFi, ipeFi
    integer :: n_coarsen, n_refine
    type(tree_node_ptr) :: tree, sibling
    logical             :: active

    call proper_nesting

    if(stagger_grid) then
      call store_faces
      call comm_faces
    end if

    n_coarsen = count(coarsen(:, :))
    n_refine = count(refine(:, :))

    ! to save memory: first coarsen then refine
    irecv=0
    isend=0
    allocate(recvstatus(MPI_STATUS_SIZE,max_blocks),recvrequest(max_blocks),&
        sendstatus(MPI_STATUS_SIZE,max_blocks),sendrequest(max_blocks))
    recvrequest=MPI_REQUEST_NULL
    sendrequest=MPI_REQUEST_NULL

    if(stagger_grid) then
      allocate(recvstatus_stg(MPI_STATUS_SIZE,max_blocks*3),&
         recvrequest_stg(max_blocks*3), sendstatus_stg(MPI_STATUS_SIZE,&
         max_blocks*3),sendrequest_stg(max_blocks*3))
      recvrequest_stg=MPI_REQUEST_NULL
      sendrequest_stg=MPI_REQUEST_NULL
    end if

    do ipe=0,npe-1
       do igrid=1,max_blocks
          if (coarsen(igrid,ipe)) then
             if (.not.associated(igrid_to_node(igrid,ipe)%node)) cycle

             tree%node => igrid_to_node(igrid,ipe)%node%parent%node
             do ic3=1,2
             do ic2=1,2
             do ic1=1,2
                sibling%node => tree%node%child(ic1,ic2,ic3)%node
                ipeFi(ic1,ic2,ic3)=sibling%node%ipe
                igridFi(ic1,ic2,ic3)=sibling%node%igrid
             end do
             end do
             end do

             ipeCo=ipeFi(1,1,1)
             igridCo=getnode(ipeCo)

             call coarsen_tree_leaf(igridCo,ipeCo,igridFi,ipeFi,active)

             call coarsen_grid_siblings(igridCo,ipeCo,igridFi,ipeFi,active)

             ! local coarsening done
             do ic3=1,2
             do ic2=1,2
             do ic1=1,2
                if (ipeFi(ic1,ic2,ic3)==ipeCo) then
                   call putnode(igridFi(ic1,ic2,ic3),ipeFi(ic1,ic2,ic3))
                   coarsen(igridFi(ic1,ic2,ic3),ipeFi(ic1,ic2,ic3))=.false.
                end if
             end do
             end do
             end do
          end if
       end do
    end do

    if (irecv>0) then
      call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
      if(stagger_grid) call MPI_WAITALL(irecv,recvrequest_stg,recvstatus_stg,&
         ierrmpi)
    end if
    if (isend>0) then
      call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
      if(stagger_grid) call MPI_WAITALL(isend,sendrequest_stg,sendstatus_stg,&
         ierrmpi)
    end if

    deallocate(recvstatus,recvrequest,sendstatus,sendrequest)
    if(stagger_grid) deallocate(recvstatus_stg,recvrequest_stg,sendstatus_stg,&
       sendrequest_stg)

    ! non-local coarsening done
    do ipe=0,npe-1
       do igrid=1,max_blocks
          if (coarsen(igrid,ipe)) then
             !if (ipe==mype) call dealloc_node(igrid) ! do not deallocate node
             ! memory preventing fragmentization of system memory as a result
             ! of frequent allocating and deallocating memory

             ! put the node (igrid number) into unused.
             call putnode(igrid,ipe)
             coarsen(igrid,ipe)=.false.
          end if
       end do
    end do

    do ipe=0,npe-1
       do igrid=1,max_blocks
          if (refine(igrid,ipe)) then

             do ic3=1,2
             do ic2=1,2
             do ic1=1,2
                igridFi(ic1,ic2,ic3)=getnode(ipe)
                ipeFi(ic1,ic2,ic3)=ipe
             end do
             end do
             end do

             call refine_tree_leaf(igridFi,ipeFi,igrid,ipe,active)

             if (ipe==mype) call refine_grids(igridFi,ipeFi,igrid,ipe,active)

             ! refinement done
             call putnode(igrid,ipe)
             refine(igrid,ipe)=.false.
          end if
       end do
    end do

    ! A crash occurs in later MPI_WAITALL when initial condition comsumes too 
    ! much time to filling new blocks with both gfortran and intel fortran compiler.
    ! This barrier cure this problem
    !TODO to find the reason
    if(.not.time_advance) call MPI_BARRIER(icomm,ierrmpi)

    if(stagger_grid) call end_comm_faces

    call get_level_range

    ! Update sfc array: igrid and ipe info in space filling curve
    call amr_Morton_order()

    call load_balance

    ! Rebuild tree connectivity
    call getigrids
    call build_connectivity

    ! Update the list of active grids
    call selectgrids
    !  grid structure now complete again.

    ! since we only filled mesh values, and advance assumes filled
    ! ghost cells, do boundary filling for the new levels
    if (time_advance) then
       call getbc(global_time+dt,0.d0,ps,iwstart,nwgc)
    else
       call getbc(global_time,0.d0,ps,iwstart,nwgc)
    end if

    
    if (use_multigrid) call mg_update_refinement(n_coarsen, n_refine)
   

    if (associated(usr_after_refine)) then
       call usr_after_refine(n_coarsen, n_refine)
    end if

  end subroutine amr_coarsen_refine

  !> For all grids on all processors, do a check on refinement flags. Make
  !> sure that neighbors will not differ more than one level of refinement.
  subroutine proper_nesting
    use mod_forest
    use mod_global_parameters
    use mod_amr_neighbors, only: find_neighbor

    logical, dimension(:,:), allocatable :: refine2
    integer :: iigrid, igrid, level, ic1,ic2,ic3, inp1,inp2,inp3, i1,i2,i3,&
        my_neighbor_type,ipe
    logical :: coarsening, pole(ndim), sendbuf(max_blocks)
    type(tree_node_ptr) :: tree, p_neighbor, my_parent, sibling, my_neighbor,&
        neighborchild

    if (nbufferx1/=0.or.nbufferx2/=0.or.nbufferx3/=0) then
       allocate(refine2(max_blocks,npe))
       call MPI_ALLREDUCE(refine,refine2,max_blocks*npe,MPI_LOGICAL,MPI_LOR,&
           icomm,ierrmpi)
       refine=refine2
    else
       sendbuf(:)=refine(:,mype)
       call MPI_ALLGATHER(sendbuf,max_blocks,MPI_LOGICAL,refine,max_blocks,&
           MPI_LOGICAL,icomm,ierrmpi)
    end if

    do level=min(levmax,refine_max_level-1),levmin+1,-1
       tree%node => level_head(level)%node
       do
          if (.not.associated(tree%node)) exit

          if (refine(tree%node%igrid,tree%node%ipe)) then
             ic1=1+modulo(tree%node%ig1-1,2);ic2=1+modulo(tree%node%ig2-1,2)
             ic3=1+modulo(tree%node%ig3-1,2);
             do inp3=ic3-2,ic3-1
             do inp2=ic2-2,ic2-1
             do inp1=ic1-2,ic1-1
                if (inp1==0.and.inp2==0.and.inp3==0) cycle
                p_neighbor%node => tree%node%parent%node
                if (inp1/=0) then
                   p_neighbor%node => p_neighbor%node%neighbor(ic1,1)%node
                   if (.not.associated(p_neighbor%node)) cycle
                end if
                if (inp2/=0) then
                   p_neighbor%node => p_neighbor%node%neighbor(ic2,2)%node
                   if (.not.associated(p_neighbor%node)) cycle
                end if
                if (inp3/=0) then
                   p_neighbor%node => p_neighbor%node%neighbor(ic3,3)%node
                   if (.not.associated(p_neighbor%node)) cycle
                end if
                if (p_neighbor%node%leaf) then
                   refine(p_neighbor%node%igrid,p_neighbor%node%ipe)=.true.
                end if
             end do
             end do
             end do
          end if

          tree%node => tree%node%next%node
       end do
    end do

    ! On each processor locally, check if grids set for coarsening are already
    ! set for refinement.

    do iigrid=1,igridstail; igrid=igrids(iigrid);
       if (refine(igrid,mype).and.coarsen(igrid,mype)) coarsen(igrid,&
          mype)=.false.
    end do

    ! For all grids on all processors, do a check on coarse refinement flags
    sendbuf(:)=coarsen(:,mype)
    call MPI_ALLGATHER(sendbuf,max_blocks,MPI_LOGICAL,coarsen,max_blocks,&
        MPI_LOGICAL,icomm,ierrmpi)

    do level=levmax,max(2,levmin),-1
       tree%node => level_head(level)%node
       do
          if (.not.associated(tree%node)) exit

          if (coarsen(tree%node%igrid,tree%node%ipe)) then
             coarsening=.true.
             my_parent%node => tree%node%parent%node

             ! are all siblings flagged for coarsen ?
    check1:  do ic3=1,2
    do ic2=1,2
    do ic1=1,2
                sibling%node => my_parent%node%child(ic1,ic2,ic3)%node
                if (sibling%node%leaf) then
                   if (coarsen(sibling%node%igrid,sibling%node%ipe)) cycle
                end if
                call unflag_coarsen_siblings
                exit check1
             end do
             end do
             end do check1

             ! Make sure that neighbors will not differ more than one level of
             ! refinement, otherwise unflag all siblings
             if (coarsening) then
    check2:     do ic3=1,2
    do ic2=1,2
    do ic1=1,2
                   sibling%node => my_parent%node%child(ic1,ic2,ic3)%node
                   do i3=ic3-2,ic3-1
                   do i2=ic2-2,ic2-1
                   do i1=ic1-2,ic1-1
                      if (i1==0.and.i2==0.and.i3==0) cycle
                      call find_neighbor(my_neighbor,my_neighbor_type, sibling,&
                         i1,i2,i3,pole)
                      select case (my_neighbor_type)
                      case (neighbor_sibling)
                         if (refine(my_neighbor%node%igrid,&
                             my_neighbor%node%ipe)) then
                            call unflag_coarsen_siblings
                            exit check2
                         else
                            cycle
                         end if
                      case (neighbor_fine)
                         neighborchild%node=>my_neighbor%node%child(1,1,&
                            1)%node
                         if (neighborchild%node%leaf) then
                            if (coarsen(neighborchild%node%igrid,&
                                neighborchild%node%ipe)) then
                               cycle
                            end if
                         end if
                         call unflag_coarsen_siblings
                         exit check2
                      end select
                   end do
                   end do
                   end do
                end do
                end do
                end do check2
             end if

          end if

          tree%node => tree%node%next%node
       end do
    end do

    contains

      subroutine unflag_coarsen_siblings

      integer :: ic1,ic2,ic3
      type(tree_node_ptr) :: sibling

      do ic3=1,2
      do ic2=1,2
      do ic1=1,2
         sibling%node => my_parent%node%child(ic1,ic2,ic3)%node
         if (sibling%node%leaf) then
            coarsen(sibling%node%igrid,sibling%node%ipe)=.false.
         end if
      end do
      end do
      end do
      coarsening=.false.

      end subroutine unflag_coarsen_siblings

  end subroutine proper_nesting

  !> coarsen sibling blocks into one block
  subroutine coarsen_grid_siblings(igrid,ipe,child_igrid,child_ipe,active)
    use mod_global_parameters
    use mod_coarsen, only: coarsen_grid
    use mod_initialize_amr, only: initial_condition
    use mod_amr_solution_node, only: alloc_node

    integer, intent(in) :: igrid, ipe
    integer, dimension(2,2,2), intent(in) :: child_igrid, child_ipe
    logical, intent(in) :: active

    integer :: igridFi, ipeFi, ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
       ixComax3, ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
        ixCoMmin1,ixCoMmin2,ixCoMmin3,ixCoMmax1,ixCoMmax2,ixCoMmax3, ic1,ic2,&
       ic3, idir

    if (ipe==mype) call alloc_node(igrid)

    ! New passive cell, coarsen from initial condition:
    if (.not. active) then
       if (ipe == mype) then
          call initial_condition(igrid)
          do ic3=1,2
          do ic2=1,2
          do ic1=1,2
          igridFi=child_igrid(ic1,ic2,ic3)
          ipeFi=child_ipe(ic1,ic2,ic3)
          !if (ipeFi==mype) then
          !   ! remove solution space of child      
          !   call dealloc_node(igridFi)
          !end if
          end do
          end do
          end do
       end if
       return
    end if

    do ic3=1,2
    do ic2=1,2
    do ic1=1,2
       igridFi=child_igrid(ic1,ic2,ic3)
       ipeFi=child_ipe(ic1,ic2,ic3)

       if (ipeFi==mype) then
          dxlevel(1)=rnode(rpdx1_,igridFi);dxlevel(2)=rnode(rpdx2_,igridFi)
          dxlevel(3)=rnode(rpdx3_,igridFi);
          if (ipe==mype) then
             ixComin1=ixMlo1+(ic1-1)*(ixMhi1-ixMlo1+1)/2
             ixComin2=ixMlo2+(ic2-1)*(ixMhi2-ixMlo2+1)/2
             ixComin3=ixMlo3+(ic3-1)*(ixMhi3-ixMlo3+1)/2;
             ixComax1=ixMhi1+(ic1-2)*(ixMhi1-ixMlo1+1)/2
             ixComax2=ixMhi2+(ic2-2)*(ixMhi2-ixMlo2+1)/2
             ixComax3=ixMhi3+(ic3-2)*(ixMhi3-ixMlo3+1)/2;

             call coarsen_grid(ps(igridFi),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
                ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid),&
                ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixComin1,ixComin2,&
                ixComin3,ixComax1,ixComax2,ixComax3)
             ! remove solution space of child
             !call dealloc_node(igridFi)
          else
             ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1;
             ixCoGmax1=ixGhi1/2+nghostcells;ixCoGmax2=ixGhi2/2+nghostcells
             ixCoGmax3=ixGhi3/2+nghostcells;
             ixCoMmin1=ixCoGmin1+nghostcells;ixCoMmin2=ixCoGmin2+nghostcells
             ixCoMmin3=ixCoGmin3+nghostcells;ixCoMmax1=ixCoGmax1-nghostcells
             ixCoMmax2=ixCoGmax2-nghostcells;ixCoMmax3=ixCoGmax3-nghostcells;
             call coarsen_grid(ps(igridFi),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
                ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,psc(igridFi),&
                 ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
                ixCoMmin1,ixCoMmin2,ixCoMmin3,ixCoMmax1,ixCoMmax2,ixCoMmax3)

             !itag=ipeFi*max_blocks+igridFi
             itag=ipeFi+igridFi
             isend=isend+1
             call MPI_ISEND(psc(igridFi)%w,1,type_coarse_block,ipe,itag, icomm,&
                sendrequest(isend),ierrmpi)
             if(stagger_grid) then
               do idir=1,ndim
                 !itag_stg=(npe+ipeFi+1)*max_blocks+igridFi*(ndir-1+idir)
                 itag_stg=(npe+ipeFi+1)+igridFi*(ndir-1+idir)
                 call MPI_ISEND(psc(igridFi)%ws,1,type_coarse_block_stg(idir,&
                    ic1,ic2,ic3),ipe,itag_stg, icomm,sendrequest_stg(isend),&
                    ierrmpi)
               end do
             end if
          end if
       else
          if (ipe==mype) then
             !itag=ipeFi*max_blocks+igridFi
             itag=ipeFi+igridFi
             irecv=irecv+1
             call MPI_IRECV(ps(igrid)%w,1,type_sub_block(ic1,ic2,ic3),ipeFi,&
                itag, icomm,recvrequest(irecv),ierrmpi)
             if(stagger_grid) then
               do idir=1,ndim
                 !itag_stg=(npe+ipeFi+1)*max_blocks+igridFi*(ndir-1+idir)
                 itag_stg=(npe+ipeFi+1)+igridFi*(ndir-1+idir)
                 call MPI_IRECV(ps(igrid)%ws,1,type_sub_block_stg(idir,ic1,ic2,&
                    ic3),ipeFi,itag_stg, icomm,recvrequest_stg(irecv),ierrmpi)
               end do
             end if
          end if
       end if
    end do
    end do
    end do

  end subroutine coarsen_grid_siblings

end module mod_coarsen_refine
