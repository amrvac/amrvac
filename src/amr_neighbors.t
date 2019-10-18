!> find neighors of level-one root blocks
subroutine find_root_neighbor(tree_neighbor,tree,i^D)
  use mod_forest
  use mod_global_parameters
  use mod_geometry

  type(tree_node_ptr) :: tree_neighbor, tree
  integer, intent(in) :: i^D

  integer :: jg^D

  jg^D=tree%node%ig^D+i^D;

  ! find the periodic grid indices, modulo(-1,10)=9
  {if (periodB(^D)) jg^D=1+modulo(jg^D-1,ng^D(1))\}

  ! pi-periodicity at pole
  select case (coordinate)
  case (spherical) {^IFTHREED
     if (poleB(1,2).and.jg2==0) then ! northpole (theta=0)
        jg2=1; jg3=1+modulo(jg3+ng3(1)/2-1,ng3(1))
     end if
     if (poleB(2,2).and.jg2==ng2(1)+1) then ! southpole (theta=pi)
        jg2=ng2(1); jg3=1+modulo(jg3+ng3(1)/2-1,ng3(1))
     end if}
  case (cylindrical)
     if (poleB(1,1).and.jg1==0) then ! cylindrical axis
        jg1=1
        {if (^D==phi_) jg^D=1+modulo(jg^D+ng^D(1)/2-1,ng^D(1))\}
     end if
  end select

  if (jg^D>=1.and.jg^D<=ng^D(1)|.and.) then
     tree_neighbor%node => tree_root(jg^D)%node
  else
     nullify(tree_neighbor%node)
  end if

end subroutine find_root_neighbor

!> find neighors of all blocks
subroutine find_neighbor(my_neighbor,my_neighbor_type,tree,i^D,pole)
  use mod_forest
  use mod_global_parameters

  type(tree_node_ptr) :: tree, my_neighbor
  integer, intent(in) :: i^D
  integer, intent(out) :: my_neighbor_type
  logical, dimension(ndim), intent(out) :: pole

  integer :: level, ig^D, ic^D, n_ic^D, inp^D

  pole=.false.
  level=tree%node%level
  if (level==1) then
     call find_root_neighbor(my_neighbor,tree,i^D)
     if (associated(my_neighbor%node)) then
       if (phi_ > 0) then
         ig^D=tree%node%ig^D;
         {if ((poleB(2,^D).and.ig^D==ng^D(1).and.i^D==1) .or. &
              (poleB(1,^D).and.ig^D==1.and.i^D==-1)) pole(^D)=.true.\}
       end if
        if (my_neighbor%node%leaf) then
           my_neighbor_type=3
        else
           my_neighbor_type=4
        end if
     else
        my_neighbor_type=1
        return
     end if
  else
     ig^D=tree%node%ig^D;

     if (phi_ > 0) then
       {if ((poleB(2,^D).and.ig^D==ng^D(level).and.i^D==1) .or. &
            (poleB(1,^D).and.ig^D==1.and.i^D==-1)) pole(^D)=.true.\}
     end if

     ! ic^D is 1 when ig^D is odd, is 2 when ig^D is even
     ic^D=1+modulo(ig^D-1,2);
     inp^D=int((ic^D+i^D+1)/2)-1;
     my_neighbor%node => tree%node%parent%node
     {if (inp^D/=0) then
        my_neighbor%node => my_neighbor%node%neighbor(ic^D,^D)%node
        if (.not.associated(my_neighbor%node)) then
           my_neighbor_type=1
           return
        end if
     end if\}
     if (my_neighbor%node%leaf) then
        my_neighbor_type=2
     else
        {if (i^D==0 .or. pole(^D)) then
           n_ic^D=ic^D
        else
           n_ic^D=3-ic^D ! switch 1 <--> 2
        end if\}
        my_neighbor%node => my_neighbor%node%child(n_ic^D)%node
        if (associated(my_neighbor%node)) then
           if (my_neighbor%node%leaf) then
              my_neighbor_type=3
           else
              my_neighbor_type=4
           end if
        else
           my_neighbor_type=0
        end if
     end if
  end if

end subroutine find_neighbor

!> asign tree node neighor
subroutine asign_tree_neighbor(tree)
  use mod_forest
  use mod_global_parameters

  type(tree_node_ptr) :: tree

  logical, dimension(ndim) :: pole
  integer :: my_neighbor_type, i^D, iside
  type(tree_node_ptr) :: my_neighbor

  {do iside=1,2
     i^DD=kr(^DD,^D)*(2*iside-3);
     call find_neighbor(my_neighbor,my_neighbor_type,tree,i^DD,pole)
     select case (my_neighbor_type)
     case (neighbor_sibling, neighbor_fine)
        tree%node%neighbor(iside,^D)%node => my_neighbor%node
        if (associated(my_neighbor%node)) then
           if (pole(^D)) then
              my_neighbor%node%neighbor(iside,^D)%node => tree%node
           else
              my_neighbor%node%neighbor(3-iside,^D)%node => tree%node
           end if
        end if
     case default
        nullify(tree%node%neighbor(iside,^D)%node)
     end select
  end do\}

end subroutine asign_tree_neighbor
