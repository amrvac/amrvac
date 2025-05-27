module mod_functions_forest

  implicit none
  private

  public :: write_forest 
  public :: read_forest 
  public :: coarsen_tree_leaf 
  public :: refine_tree_leaf 
  public :: change_ipe_tree_leaf
  public :: init_forest_root


contains


  !> build root forest
  subroutine init_forest_root
    use mod_forest
    use mod_global_parameters
    use mod_space_filling_curve
    use mod_amr_solution_node, only: getnode 
    use mod_amr_neighbors, only: find_root_neighbor
 
    integer :: ig1,ig2,ig3, level, igrid, ipe
    integer :: iside, i1,i2,i3, Morton_no, isfc
  
  
    level=1
    Morton_no=0
    ipe=0
    nparents=0
    nleafs=ng1(1)*ng2(1)*ng3(1)
    nleafs_active=nleafs
    nleafs_level(1)=ng1(1)*ng2(1)*ng3(1)
    nleafs_level(2:nlevelshi)=0
    call get_Morton_range
    ! Generate Morton-order space-filling curve to connect level 1 blocks
    call level1_Morton_order
    do isfc=1,nglev1
       ig1=sfc_iglevel1(1,isfc)
       ig2=sfc_iglevel1(2,isfc)
       ig3=sfc_iglevel1(3,isfc)
       Morton_no=Morton_no+1
       if (Morton_no>Morton_stop(ipe)) ipe=ipe+1
       igrid=getnode(ipe)
       if (ipe==mype) then
          sfc_to_igrid(Morton_no)=igrid
          igrid_to_sfc(igrid)=Morton_no
       end if
       call init_tree_leaf(tree_root(ig1,ig2,ig3),ig1,ig2,ig3,level,igrid,ipe,&
          .true.)
    end do
  
    ! update root neighbor
    do ig3=1,ng3(1)
    do ig2=1,ng2(1)
    do ig1=1,ng1(1)
       do iside=1,2
          i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3)
          i3=kr(3,1)*(2*iside-3);
          call find_root_neighbor(tree_root(ig1,ig2,ig3)%node%neighbor(iside,&
             1), tree_root(ig1,ig2,ig3),i1,i2,i3)
       end do
       do iside=1,2
          i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3)
          i3=kr(3,2)*(2*iside-3);
          call find_root_neighbor(tree_root(ig1,ig2,ig3)%node%neighbor(iside,&
             2), tree_root(ig1,ig2,ig3),i1,i2,i3)
       end do
       do iside=1,2
          i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3)
          i3=kr(3,3)*(2*iside-3);
          call find_root_neighbor(tree_root(ig1,ig2,ig3)%node%neighbor(iside,&
             3), tree_root(ig1,ig2,ig3),i1,i2,i3)
       end do
    end do
    end do
    end do
  
    ! This call is here to ensure the sfc array is initialized
    call amr_Morton_order()
  
  end subroutine init_forest_root
  
  subroutine init_tree_leaf(tree,ig1,ig2,ig3,level,igrid,ipe,active)
    use mod_forest
    implicit none
  
    type(tree_node_ptr) :: tree
    integer, intent(in) :: ig1,ig2,ig3, level, igrid, ipe
    logical, intent(in) :: active
    integer :: ic1,ic2,ic3
  
    allocate(tree%node)
  
    tree%node%ig1=ig1;tree%node%ig2=ig2;tree%node%ig3=ig3;
    tree%node%level=level
    tree%node%igrid=igrid
    tree%node%ipe=ipe
  
    tree%node%leaf=.true.
    tree%node%active=active
  
    nullify(tree%node%parent%node)
    do ic3=1,2
    do ic2=1,2
    do ic1=1,2
       nullify(tree%node%child(ic1,ic2,ic3)%node)
    end do
    end do
    end do
  
    call add_to_linked_list(level,tree)
  
    ! initialize neighbor pointers
    nullify(tree%node%neighbor(1,1)%node,tree%node%neighbor(1,2)%node,&
       tree%node%neighbor(1,3)%node,tree%node%neighbor(2,1)%node,&
       tree%node%neighbor(2,2)%node,tree%node%neighbor(2,3)%node)
  
    igrid_to_node(igrid,ipe)%node => tree%node
  
  end subroutine init_tree_leaf
  
  subroutine coarsen_tree_leaf(igrid,ipe,child_igrid,child_ipe,active)
    use mod_forest
    implicit none
  
    integer, intent(in) :: igrid, ipe
    integer, dimension(2,2,2), intent(in) :: child_igrid, child_ipe
    logical, intent(out) :: active
  
  
    integer :: level, ic1,ic2,ic3, child_level, iside, iotherside, vote
    type(tree_node_ptr) :: tree, child, child_neighbor
  
    tree%node => igrid_to_node(child_igrid(1,1,1),child_ipe(1,1,&
       1))%node%parent%node
    level=tree%node%level
  
    call add_to_linked_list(level,tree)
  
    child_level=level+1
    vote=0
  
    do ic3=1,2
    do ic2=1,2
    do ic1=1,2
       child%node => tree%node%child(ic1,ic2,ic3)%node
  
    !  vote for active:
       if(child%node%active) vote=vote+1
  
       call delete_from_linked_list(child_level,child)
  
       ! update neighbor pointers
       iside=ic1
       child_neighbor%node => child%node%neighbor(iside,1)%node
       if (associated(child_neighbor%node)) then
          if (child%node%ig1==child_neighbor%node%ig1) then ! pole
             nullify(child_neighbor%node%neighbor(iside,1)%node)
          else
             iotherside=3-iside
             nullify(child_neighbor%node%neighbor(iotherside,1)%node)
          end if
       end if
       iside=ic2
       child_neighbor%node => child%node%neighbor(iside,2)%node
       if (associated(child_neighbor%node)) then
          if (child%node%ig2==child_neighbor%node%ig2) then ! pole
             nullify(child_neighbor%node%neighbor(iside,2)%node)
          else
             iotherside=3-iside
             nullify(child_neighbor%node%neighbor(iotherside,2)%node)
          end if
       end if
       iside=ic3
       child_neighbor%node => child%node%neighbor(iside,3)%node
       if (associated(child_neighbor%node)) then
          if (child%node%ig3==child_neighbor%node%ig3) then ! pole
             nullify(child_neighbor%node%neighbor(iside,3)%node)
          else
             iotherside=3-iside
             nullify(child_neighbor%node%neighbor(iotherside,3)%node)
          end if
       end if
  
       nullify(tree%node%child(ic1,ic2,ic3)%node)
       deallocate(igrid_to_node(child_igrid(ic1,ic2,ic3),child_ipe(ic1,ic2,&
          ic3))%node)
    end do
    end do
    end do
  
    tree%node%leaf=.true.
    tree%node%igrid=igrid
    tree%node%ipe=ipe
    igrid_to_node(igrid,ipe)%node => tree%node
  
    !  Count the vote and set active/passive state:
  
    if (vote /= 2**3) then 
    !if (vote == 0) then 
       tree%node%active = .false.
       nleafs_active = nleafs_active - vote
    else
       tree%node%active = .true.
       nleafs_active = nleafs_active - vote + 1
    end if
    active = tree%node%active
  
    nleafs=nleafs-2**3+1
    nparents=nparents-1
    nleafs_level(child_level)=nleafs_level(child_level)-2**3
    nleafs_level(level)=nleafs_level(level)+1
  
  end subroutine coarsen_tree_leaf
  
  subroutine refine_tree_leaf(child_igrid,child_ipe,igrid,ipe,active)
    use mod_forest
    use mod_global_parameters
    use mod_amr_neighbors, only: find_neighbor  
 
    integer, dimension(2,2,2), intent(in) :: child_igrid, child_ipe
    integer, intent(in) :: igrid, ipe
    logical, intent(out):: active
  
    integer :: ig1,ig2,ig3, level, i1,i2,i3, ic1,ic2,ic3, child_ig1,child_ig2,&
       child_ig3, child_level, iside
    integer :: my_neighbor_type
    logical, dimension(ndim) :: pole
    type(tree_node_ptr) :: tree, child, my_neighbor
  
    tree%node => igrid_to_node(igrid,ipe)%node
    ig1=tree%node%ig1;ig2=tree%node%ig2;ig3=tree%node%ig3;
    level=tree%node%level
    active=tree%node%active
  
    tree%node%ipe=-1
    tree%node%igrid=0
    tree%node%leaf=.false.
    tree%node%active=.true.
  
    call delete_from_linked_list(level,tree)
  
    child_level=level+1
  
    do ic3=1,2
    do ic2=1,2
    do ic1=1,2
       child_ig1=2*(ig1-1)+ic1;child_ig2=2*(ig2-1)+ic2
       child_ig3=2*(ig3-1)+ic3;
       call init_tree_leaf(child,child_ig1,child_ig2,child_ig3,child_level,&
           child_igrid(ic1,ic2,ic3),child_ipe(ic1,ic2,ic3),active)
  
       igrid_to_node(child_igrid(ic1,ic2,ic3),child_ipe(ic1,ic2,&
          ic3))%node => child%node
  
       tree%node%child(ic1,ic2,ic3)%node => child%node
       child%node%parent%node => tree%node
    end do
    end do
    end do
  
    ! update neighbor pointers
    do ic3=1,2
    do ic2=1,2
    do ic1=1,2
       child%node => tree%node%child(ic1,ic2,ic3)%node
       iside=ic1
       i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);i3=kr(3,1)*(2*iside-3);
       call find_neighbor(my_neighbor,my_neighbor_type,child,i1,i2,i3,pole)
       select case (my_neighbor_type)
       case (neighbor_sibling, neighbor_fine)
          child%node%neighbor(iside,1)%node => my_neighbor%node
          if (pole(1)) then
             my_neighbor%node%neighbor(iside,1)%node => child%node
          else
             my_neighbor%node%neighbor(3-iside,1)%node => child%node
          end if
       case default
          nullify(child%node%neighbor(iside,1)%node)
       end select
       child%node%neighbor(3-ic1,1)%node=>tree%node%child(3-ic1,ic2,ic3)%node
       iside=ic2
       i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);i3=kr(3,2)*(2*iside-3);
       call find_neighbor(my_neighbor,my_neighbor_type,child,i1,i2,i3,pole)
       select case (my_neighbor_type)
       case (neighbor_sibling, neighbor_fine)
          child%node%neighbor(iside,2)%node => my_neighbor%node
          if (pole(2)) then
             my_neighbor%node%neighbor(iside,2)%node => child%node
          else
             my_neighbor%node%neighbor(3-iside,2)%node => child%node
          end if
       case default
          nullify(child%node%neighbor(iside,2)%node)
       end select
       child%node%neighbor(3-ic2,2)%node=>tree%node%child(ic1,3-ic2,ic3)%node
       iside=ic3
       i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3);i3=kr(3,3)*(2*iside-3);
       call find_neighbor(my_neighbor,my_neighbor_type,child,i1,i2,i3,pole)
       select case (my_neighbor_type)
       case (neighbor_sibling, neighbor_fine)
          child%node%neighbor(iside,3)%node => my_neighbor%node
          if (pole(3)) then
             my_neighbor%node%neighbor(iside,3)%node => child%node
          else
             my_neighbor%node%neighbor(3-iside,3)%node => child%node
          end if
       case default
          nullify(child%node%neighbor(iside,3)%node)
       end select
       child%node%neighbor(3-ic3,3)%node=>tree%node%child(ic1,ic2,3-ic3)%node
    end do
    end do
    end do
  
    nleafs=nleafs+2**3-1
    nparents=nparents+1
    nleafs_level(child_level)=nleafs_level(child_level)+2**3
    nleafs_level(level)=nleafs_level(level)-1
  
    if (active) nleafs_active = nleafs_active + 2**3-1
  
  end subroutine refine_tree_leaf
  
  subroutine change_ipe_tree_leaf(recv_igrid,recv_ipe,send_igrid,send_ipe)
    use mod_forest
    implicit none
  
    integer, intent(in) :: recv_igrid, recv_ipe, send_igrid, send_ipe
  
    type(tree_node_ptr) :: tree
  
    tree%node => igrid_to_node(send_igrid,send_ipe)%node
  
    tree%node%igrid=recv_igrid
    tree%node%ipe=recv_ipe
  
    nullify(igrid_to_node(send_igrid,send_ipe)%node)
    igrid_to_node(recv_igrid,recv_ipe)%node => tree%node
  
  end subroutine change_ipe_tree_leaf
  
  subroutine add_to_linked_list(level,tree)
    use mod_forest
    implicit none
  
    integer, intent(in) :: level
    type(tree_node_ptr) :: tree
  
    nullify(tree%node%next%node)
    if (associated(level_head(level)%node)) then
       tree%node%prev%node => level_tail(level)%node
       level_tail(level)%node%next%node => tree%node
       level_tail(level)%node => tree%node
    else
       level_head(level)%node => tree%node
       level_tail(level)%node => tree%node
       nullify(tree%node%prev%node)
    end if
  
  end subroutine add_to_linked_list
  
  subroutine delete_from_linked_list(level,tree)
    use mod_forest
    implicit none
  
    integer, intent(in) :: level
    type(tree_node_ptr) :: tree
  
    type(tree_node_ptr) :: next, prev
  
    prev%node => tree%node%prev%node
    next%node => tree%node%next%node
    if (associated(next%node).and.associated(prev%node)) then
       prev%node%next%node => next%node
       next%node%prev%node => prev%node
    else if (associated(prev%node)) then
       level_tail(level)%node => prev%node
       nullify(prev%node%next%node)
    else if (associated(next%node)) then
       level_head(level)%node => next%node
       nullify(next%node%prev%node)
    else
       nullify(level_head(level)%node)
       nullify(level_tail(level)%node)
    end if
  
  end subroutine delete_from_linked_list
  
  subroutine write_forest(file_handle)
    use mod_forest
    use mod_global_parameters
  
    integer, intent(in) :: file_handle
  
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: ig1,ig2,ig3,isfc
  
    do isfc=1,nglev1
       ig1=sfc_iglevel1(1,isfc)
       ig2=sfc_iglevel1(2,isfc)
       ig3=sfc_iglevel1(3,isfc)
       call write_node(tree_root(ig1,ig2,ig3))
    end do
  
    contains
  
      recursive subroutine write_node(tree)
        implicit none
  
        type(tree_node_ptr) :: tree
  
        integer :: ic1,ic2,ic3
  
        call MPI_FILE_WRITE(file_handle,tree%node%leaf,1,MPI_LOGICAL,status,&
           ierrmpi)
  
        if (.not.tree%node%leaf) then
           do ic3=1,2
           do ic2=1,2
           do ic1=1,2
              call write_node(tree%node%child(ic1,ic2,ic3))
           end do
           end do
           end do
        end if
  
      end subroutine write_node
  
  end subroutine write_forest
  
  subroutine read_forest(file_handle)
    use mod_forest
    use mod_global_parameters
    use mod_space_filling_curve
    use mod_amr_solution_node, only: getnode 
    use mod_functions_connectivity, only: get_level_range,build_connectivity,&
       getigrids 
    use mod_amr_neighbors, only: asign_tree_neighbor
   
    integer, intent(in) :: file_handle
  
    integer, dimension(MPI_STATUS_SIZE) :: status
    !integer :: ig^D, level, size_logical, Morton_no, igrid, ipe
    integer :: ig1,ig2,ig3, level, Morton_no, igrid, ipe, isfc
  
    Morton_no=0
    ipe=0
    level=1
    nleafs_level(1:nlevelshi) = 0
    nparents = 0
  
    call get_Morton_range
    call level1_Morton_order
    do isfc=1,nglev1
       ig1=sfc_iglevel1(1,isfc)
       ig2=sfc_iglevel1(2,isfc)
       ig3=sfc_iglevel1(3,isfc)
       allocate(tree_root(ig1,ig2,ig3)%node)
       nullify(tree_root(ig1,ig2,ig3)%node%parent%node)
       call read_node(tree_root(ig1,ig2,ig3),ig1,ig2,ig3,level)
    end do
  
    call get_level_range
  
    ! Rebuild tree connectivity
    call getigrids
    call build_connectivity
  
    ! This call is here to ensure the sfc array is initialized
    call amr_Morton_order()
  
    contains
  
      recursive subroutine read_node(tree,ig1,ig2,ig3,level)
        implicit none
  
        type(tree_node_ptr) :: tree
        integer, intent(in) :: ig1,ig2,ig3, level
  
        logical :: leaf
        integer :: ic1,ic2,ic3, child_ig1,child_ig2,child_ig3, child_level
  
        if (mype==0) then
          call MPI_FILE_READ(file_handle,leaf,1,MPI_LOGICAL, status,ierrmpi)
        end if
        if (npe>1)  call MPI_BCAST(leaf,1,MPI_LOGICAL,0,icomm,ierrmpi)
  
        tree%node%leaf=leaf
        tree%node%ig1=ig1;tree%node%ig2=ig2;tree%node%ig3=ig3;
        tree%node%level=level
        tree%node%active=.true. .and. leaf
  
        do ic3=1,2
        do ic2=1,2
        do ic1=1,2
           nullify(tree%node%child(ic1,ic2,ic3)%node)
        end do
        end do
        end do
        nullify(tree%node%neighbor(1,1)%node,tree%node%neighbor(1,2)%node,&
           tree%node%neighbor(1,3)%node,tree%node%neighbor(2,1)%node,&
           tree%node%neighbor(2,2)%node,tree%node%neighbor(2,3)%node)
        nullify(tree%node%next%node,tree%node%prev%node)
  
        call asign_tree_neighbor(tree)
  
        if (leaf) then
           call add_to_linked_list(level,tree)
           nleafs_level(level) = nleafs_level(level) + 1
  
           Morton_no=Morton_no+1
           if (Morton_no>Morton_stop(ipe)) ipe=ipe+1
           igrid=getnode(ipe)
           tree%node%igrid=igrid
           tree%node%ipe=ipe
           igrid_to_node(igrid,ipe)%node => tree%node
           if (ipe==mype) sfc_to_igrid(Morton_no)=igrid
         else
           nparents = nparents + 1
           tree%node%igrid=0
           tree%node%ipe=-1
           child_level=level+1
           do ic3=1,2
           do ic2=1,2
           do ic1=1,2
              child_ig1=2*(ig1-1)+ic1;child_ig2=2*(ig2-1)+ic2
              child_ig3=2*(ig3-1)+ic3;
              allocate(tree%node%child(ic1,ic2,ic3)%node)
              tree%node%child(ic1,ic2,ic3)%node%parent%node => tree%node
              call read_node(tree%node%child(ic1,ic2,ic3),child_ig1,child_ig2,&
                 child_ig3,child_level)
           end do
           end do
           end do
        end if
  
      end subroutine read_node
  
  end subroutine read_forest


end module mod_functions_forest
