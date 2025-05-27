module mod_amr_neighbors

  implicit none
  private

  public :: find_neighbor 
  public :: find_root_neighbor 
  public :: asign_tree_neighbor 
 
contains


  !> find neighors of level-one root blocks
  subroutine find_root_neighbor(tree_neighbor,tree,i1,i2,i3)
    use mod_forest
    use mod_global_parameters
    use mod_geometry
  
    type(tree_node_ptr) :: tree_neighbor, tree
    integer, intent(in) :: i1,i2,i3
  
    integer :: jg1,jg2,jg3
  
    jg1=tree%node%ig1+i1;jg2=tree%node%ig2+i2;jg3=tree%node%ig3+i3;
  
    ! find the periodic grid indices, modulo(-1,10)=9
    if (periodB(1)) jg1=1+modulo(jg1-1,ng1(1))
    if (periodB(2)) jg2=1+modulo(jg2-1,ng2(1))
    if (periodB(3)) jg3=1+modulo(jg3-1,ng3(1))
  
    ! pi-periodicity at pole
    select case (coordinate)
    case (spherical) 
       if (poleB(1,2).and.jg2==0) then ! northpole (theta=0)
          jg2=1; jg3=1+modulo(jg3+ng3(1)/2-1,ng3(1))
       end if
       if (poleB(2,2).and.jg2==ng2(1)+1) then ! southpole (theta=pi)
          jg2=ng2(1); jg3=1+modulo(jg3+ng3(1)/2-1,ng3(1))
       end if
    case (cylindrical)
       if (poleB(1,1).and.jg1==0) then ! cylindrical axis
          jg1=1
          if (1==phi_) jg1=1+modulo(jg1+ng1(1)/2-1,ng1(1))
          if (2==phi_) jg2=1+modulo(jg2+ng2(1)/2-1,ng2(1))
          if (3==phi_) jg3=1+modulo(jg3+ng3(1)/2-1,ng3(1))
       end if
    end select
  
    if (jg1>=1.and.jg1<=ng1(1).and.jg2>=1.and.jg2<=ng2(1).and.jg3>=&
       1.and.jg3<=ng3(1)) then
       tree_neighbor%node => tree_root(jg1,jg2,jg3)%node
    else
       nullify(tree_neighbor%node)
    end if
  
  end subroutine find_root_neighbor
  
  !> find neighors of all blocks
  subroutine find_neighbor(my_neighbor,my_neighbor_type,tree,i1,i2,i3,pole)
    use mod_forest
    use mod_global_parameters
  
    type(tree_node_ptr) :: tree, my_neighbor
    integer, intent(in) :: i1,i2,i3
    integer, intent(out) :: my_neighbor_type
    logical, dimension(ndim), intent(out) :: pole
  
    integer :: level, ig1,ig2,ig3, ic1,ic2,ic3, n_ic1,n_ic2,n_ic3, inp1,inp2,&
       inp3
  
    pole=.false.
    level=tree%node%level
    if (level==1) then
       call find_root_neighbor(my_neighbor,tree,i1,i2,i3)
       if (associated(my_neighbor%node)) then
         if (phi_ > 0) then
           ig1=tree%node%ig1;ig2=tree%node%ig2;ig3=tree%node%ig3;
           if ((poleB(2,1).and.ig1==ng1(1).and.i1==1) .or. (poleB(1,&
              1).and.ig1==1.and.i1==-1)) pole(1)=.true.
           if ((poleB(2,2).and.ig2==ng2(1).and.i2==1) .or. (poleB(1,&
              2).and.ig2==1.and.i2==-1)) pole(2)=.true.
           if ((poleB(2,3).and.ig3==ng3(1).and.i3==1) .or. (poleB(1,&
              3).and.ig3==1.and.i3==-1)) pole(3)=.true.
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
       ig1=tree%node%ig1;ig2=tree%node%ig2;ig3=tree%node%ig3;
  
       if (phi_ > 0) then
         if ((poleB(2,1).and.ig1==ng1(level).and.i1==1) .or. (poleB(1,&
            1).and.ig1==1.and.i1==-1)) pole(1)=.true.
         if ((poleB(2,2).and.ig2==ng2(level).and.i2==1) .or. (poleB(1,&
            2).and.ig2==1.and.i2==-1)) pole(2)=.true.
         if ((poleB(2,3).and.ig3==ng3(level).and.i3==1) .or. (poleB(1,&
            3).and.ig3==1.and.i3==-1)) pole(3)=.true.
       end if
  
       ! ic^D is 1 when ig^D is odd, is 2 when ig^D is even
       ic1=1+modulo(ig1-1,2);ic2=1+modulo(ig2-1,2);ic3=1+modulo(ig3-1,2);
       inp1=int((ic1+i1+1)/2)-1;inp2=int((ic2+i2+1)/2)-1
       inp3=int((ic3+i3+1)/2)-1;
       my_neighbor%node => tree%node%parent%node
       if (inp1/=0) then
          my_neighbor%node => my_neighbor%node%neighbor(ic1,1)%node
          if (.not.associated(my_neighbor%node)) then
             my_neighbor_type=1
             return
          end if
       end if
       if (inp2/=0) then
          my_neighbor%node => my_neighbor%node%neighbor(ic2,2)%node
          if (.not.associated(my_neighbor%node)) then
             my_neighbor_type=1
             return
          end if
       end if
       if (inp3/=0) then
          my_neighbor%node => my_neighbor%node%neighbor(ic3,3)%node
          if (.not.associated(my_neighbor%node)) then
             my_neighbor_type=1
             return
          end if
       end if
       if (my_neighbor%node%leaf) then
          my_neighbor_type=2
       else
          if (i1==0 .or. pole(1)) then
             n_ic1=ic1
          else
             n_ic1=3-ic1 ! switch 1 <--> 2
          end if
          if (i2==0 .or. pole(2)) then
             n_ic2=ic2
          else
             n_ic2=3-ic2 ! switch 1 <--> 2
          end if
          if (i3==0 .or. pole(3)) then
             n_ic3=ic3
          else
             n_ic3=3-ic3 ! switch 1 <--> 2
          end if
          my_neighbor%node => my_neighbor%node%child(n_ic1,n_ic2,n_ic3)%node
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
    integer :: my_neighbor_type, i1,i2,i3, iside
    type(tree_node_ptr) :: my_neighbor
  
    do iside=1,2
       i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);i3=kr(3,1)*(2*iside-3);
       call find_neighbor(my_neighbor,my_neighbor_type,tree,i1,i2,i3,pole)
       select case (my_neighbor_type)
       case (neighbor_sibling, neighbor_fine)
          tree%node%neighbor(iside,1)%node => my_neighbor%node
          if (associated(my_neighbor%node)) then
             if (pole(1)) then
                my_neighbor%node%neighbor(iside,1)%node => tree%node
             else
                my_neighbor%node%neighbor(3-iside,1)%node => tree%node
             end if
          end if
       case default
          nullify(tree%node%neighbor(iside,1)%node)
       end select
    end do
    do iside=1,2
       i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);i3=kr(3,2)*(2*iside-3);
       call find_neighbor(my_neighbor,my_neighbor_type,tree,i1,i2,i3,pole)
       select case (my_neighbor_type)
       case (neighbor_sibling, neighbor_fine)
          tree%node%neighbor(iside,2)%node => my_neighbor%node
          if (associated(my_neighbor%node)) then
             if (pole(2)) then
                my_neighbor%node%neighbor(iside,2)%node => tree%node
             else
                my_neighbor%node%neighbor(3-iside,2)%node => tree%node
             end if
          end if
       case default
          nullify(tree%node%neighbor(iside,2)%node)
       end select
    end do
    do iside=1,2
       i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3);i3=kr(3,3)*(2*iside-3);
       call find_neighbor(my_neighbor,my_neighbor_type,tree,i1,i2,i3,pole)
       select case (my_neighbor_type)
       case (neighbor_sibling, neighbor_fine)
          tree%node%neighbor(iside,3)%node => my_neighbor%node
          if (associated(my_neighbor%node)) then
             if (pole(3)) then
                my_neighbor%node%neighbor(iside,3)%node => tree%node
             else
                my_neighbor%node%neighbor(3-iside,3)%node => tree%node
             end if
          end if
       case default
          nullify(tree%node%neighbor(iside,3)%node)
       end select
    end do
  
  end subroutine asign_tree_neighbor

end module mod_amr_neighbors
