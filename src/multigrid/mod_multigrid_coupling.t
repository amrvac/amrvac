!> Module to couple the octree-mg library to AMRVAC. This file uses the VACPP
!> preprocessor, but its use is kept to a minimum.
module mod_multigrid_coupling
  use m_octree_mg

  implicit none
  public

  !> Data structure containing the multigrid tree.
  type(mg_t) :: mg

contains

  !> Setup multigrid for usage
  subroutine mg_setup_multigrid()
    use mod_global_parameters

    if (ndim == 1) &
         error stop "Multigrid not available in 1D"

    if (ndim /= mg_ndim) &
         error stop "Multigrid module was compiled for different ndim"

    if (typeaxial /= "slab") &
         error stop "Multigrid only supports slab geometry"

    if (any([ block_nx^D ] /= block_nx1)) &
         error stop "Multigrid requires all block_nx to be equal"

    call mg_comm_init(mg)
    call mg_tree_from_amrvac(mg)
  end subroutine mg_setup_multigrid

  !> Set multigrid boundary conditions according to variable iw
  subroutine mg_copy_boundary_conditions(mg, iw)
    use mod_global_parameters
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: iw
    character(len=std_len)    :: bnd_name(num_neighbors)
    integer                   :: n

    do n = 1, num_neighbors
       select case (typeboundary(iw, n))
       case ('symm')
          mg%bc(n)%bc_type = bc_neumann
          mg%bc(n)%bc_value = 0.0_dp
       case ('asymm')
          mg%bc(n)%bc_type = bc_dirichlet
          mg%bc(n)%bc_value = 0.0_dp
       case ('cont')
          mg%bc(n)%bc_type = bc_continuous
          mg%bc(n)%bc_value = 0.0_dp ! Not needed
       case ('periodic')
          ! Nothing to do here
       case default
          print *, "Not a standard: ", trim(typeboundary(iw, n))
          error stop "You have to set a user-defined boundary method"
       end select
    end do
  end subroutine mg_copy_boundary_conditions

  !> If the grid has changed, rebuild the full multigrid tree
  subroutine mg_update_refinement(n_coarsen, n_refine)
    use mod_global_parameters
    integer, intent(in) :: n_coarsen
    integer, intent(in) :: n_refine

    ! Don't build multigrid tree while doing initial refinement
    if (.not. time_advance) return

    if (.not. mg%is_allocated) then
       call mg_tree_from_amrvac(mg)
    else if (n_coarsen + n_refine > 0) then
       call mg_deallocate_storage(mg)
       call mg_tree_from_amrvac(mg)
    end if
  end subroutine mg_update_refinement

  !> Copy a variable to the multigrid tree
  subroutine mg_copy_to_tree(iw_from, iw_to)
    use mod_global_parameters
    use mod_forest
    integer, intent(in)      :: iw_from !< Variable to use as right-hand side
    integer, intent(in)      :: iw_to   !< Copy to this variable
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl
    type(tree_node), pointer :: pnode

    if (.not. mg%is_allocated) &
         error stop "multigrid tree not allocated yet"

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

{^IFTWOD
       mg%boxes(id)%cc(1:nc, 1:nc, iw_to) = &
            pw(igrid)%w(ixMlo1:ixMhi1, ixMlo2:ixMhi2, iw_from)
}
{^IFTHREED
       mg%boxes(id)%cc(1:nc, 1:nc, 1:nc, iw_to) = &
            pw(igrid)%w(ixMlo1:ixMhi1, ixMlo2:ixMhi2, ixMlo3:ixMhi3, iw_from)
}
    end do
  end subroutine mg_copy_to_tree

  !> Copy a variable from the multigrid tree
  subroutine mg_copy_from_tree(iw_from, iw_to)
    use mod_global_parameters
    use mod_forest
    integer, intent(in)      :: iw_from !< Variable to use as right-hand side
    integer, intent(in)      :: iw_to   !< Copy to this variable
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl
    type(tree_node), pointer :: pnode

    if (.not. mg%is_allocated) &
         error stop "multigrid tree not allocated yet"

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

{^IFTWOD
       pw(igrid)%w(ixMlo1:ixMhi1, ixMlo2:ixMhi2, iw_to) = &
            mg%boxes(id)%cc(1:nc, 1:nc, iw_from)
}
{^IFTHREED
       pw(igrid)%w(ixMlo1:ixMhi1, ixMlo2:ixMhi2, ixMlo3:ixMhi3, iw_to) = &
            mg%boxes(id)%cc(1:nc, 1:nc, 1:nc, iw_from)
}
    end do
  end subroutine mg_copy_from_tree

  !> Copy from multigrid tree with one layer of ghost cells. Corner ghost cells
  !> are not used/set.
  subroutine mg_copy_from_tree_gc(iw_from, iw_to)
    use mod_global_parameters
    use mod_forest
    integer, intent(in)      :: iw_from !< Variable to use as right-hand side
    integer, intent(in)      :: iw_to   !< Copy to this variable
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl
    type(tree_node), pointer :: pnode

    if (.not. mg%is_allocated) &
         error stop "multigrid tree not allocated yet"

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

{^IFTWOD
       pw(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_from)
}
{^IFTHREED
       pw(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
            ixMlo3-1:ixMhi3+1, iw_to) = &
            mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_from)
}
    end do
  end subroutine mg_copy_from_tree_gc

  !> Generate a multigrid tree that includes the amrvac tree, but also contains
  !> coarser grid levels. A number of checks has already been performed in
  !> mg_setup_multigrid, so we don't repeat these checks here.
  subroutine mg_tree_from_amrvac(mg)
    use mod_forest
    use mod_global_parameters
    type(mg_t), intent(inout)        :: mg
    integer                          :: i, n, id, ix(ndim)
    integer                          :: n_boxes_total, i_c, c_id, c_ix(ndim)
    integer                          :: min_lvl, lvl
    integer                          :: nb, nb_ix, nb_dim
    integer                          :: n_finer
    type(tree_node), pointer         :: pnode, pnode_ch
    type(tree_node_ptr), allocatable :: id_to_node(:)
    real(dp)                         :: dr_coarse

    ! Estimate number of finer blocks
    n_finer = nparents+nleafs

    call mg_build_rectangle(mg, [ domain_nx^D ], block_nx1, dx(:,1), &
         [ xprobmin^D ], periodB, n_finer)

    mg%highest_lvl = levmax
    n_boxes_total = mg%n_boxes + n_finer

    ! To link the two trees
    allocate(id_to_node(n_boxes_total))

    ! Link base level
    do i = 1, size(mg%lvls(1)%ids)
       id = mg%lvls(1)%ids(i)
       ix = mg%boxes(id)%ix

       pnode               => tree_root({ix(^D)})%node
       pnode%id            =  id
       id_to_node(id)%node => pnode
       mg%boxes(id)%rank   =  pnode%ipe
    end do

    ! Add refinement
    do lvl = 1, mg%highest_lvl
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          pnode => id_to_node(id)%node

          if (.not. pnode%leaf) then
             call mg_add_children(mg, id)

             do i_c = 1, num_children
                c_id = mg%boxes(id)%children(i_c)
                c_ix = child_dix(:, i_c) + 1
                pnode_ch => pnode%child({c_ix(^D)})%node
                id_to_node(c_id)%node => pnode_ch
                pnode_ch%id = c_id
                mg%boxes(c_id)%rank = pnode_ch%ipe
             end do
          end if
       end do

       call mg_set_leaves_parents(mg%boxes, mg%lvls(lvl))

       if (lvl < mg%highest_lvl) then
          call mg_set_next_level_ids(mg, lvl)
          call mg_set_neighbors_lvl(mg, lvl+1)
       end if
    end do

    ! Store boxes with refinement boundaries (from the coarse side)
    do lvl = 1, mg%highest_lvl
       call mg_set_refinement_boundaries(mg%boxes, mg%lvls(lvl))
    end do

    ! Assign boxes to MPI processes
    call mg_load_balance_parents(mg)

    ! Allocate storage for boxes owned by this process
    call mg_allocate_storage(mg)

  end subroutine mg_tree_from_amrvac

end module mod_multigrid_coupling
