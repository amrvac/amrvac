!> Module to couple the octree-mg library to AMRVAC. This file uses the VACPP
!> preprocessor, but its use is kept to a minimum.
module mod_multigrid_coupling
  {^IFONED
  use m_octree_mg_1d
  }
  {^IFTWOD
  use m_octree_mg_2d
  }
  {^IFTHREED
  use m_octree_mg_3d
  }


  implicit none
  public

  !> Data structure containing the multigrid tree.
  type(mg_t) :: mg

  !> If defined, this routine is called after a new multigrid tree is
  !> constructed.
  procedure(after_new_tree), pointer :: mg_after_new_tree => null()

  interface
     subroutine after_new_tree()
     end subroutine after_new_tree
  end interface

contains

  !> Setup multigrid for usage
  subroutine mg_setup_multigrid()
    use mod_global_parameters
    use mod_geometry

    if (ndim /= mg_ndim) &
         error stop "Multigrid module was compiled for different ndim"

    select case (coordinate)
    case (Cartesian)
       continue
    case (cylindrical)
       if (ndim == 3) error stop "Multigrid does not support cylindrical 3D"
       mg%geometry_type = mg_cylindrical
    case default
       error stop "Multigrid does not support your geometry"
    end select

    if (any([ block_nx^D ] /= block_nx1)) &
         error stop "Multigrid requires all block_nx to be equal"

    call mg_comm_init(mg)
    call mg_set_methods(mg)
    call mg_tree_from_amrvac(mg)
  end subroutine mg_setup_multigrid

  !> Set multigrid boundary conditions for the solution according to variable iw
  subroutine mg_copy_boundary_conditions(mg, iw)
    use mod_global_parameters
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: iw
    character(len=std_len)    :: bnd_name(mg_num_neighbors)
    integer                   :: n

    do n = 1, mg_num_neighbors
       select case (typeboundary(iw, n))
       case (bc_symm)
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_asymm)
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_cont)
          mg%bc(n, mg_iphi)%bc_type = mg_bc_continuous
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp ! Not needed
       case (bc_periodic)
          ! Nothing to do here
       case default
          write(*,*) "mg_copy_boundary_conditions: unknown boundary type"
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

  !> Copy a variable to the multigrid tree, including a layer of ghost cells
  subroutine mg_copy_to_tree(iw_from, iw_to, restrict, restrict_gc, factor, state_from)
    use mod_global_parameters
    use mod_forest
    integer, intent(in)            :: iw_from     !< Variable to use as right-hand side
    integer, intent(in)            :: iw_to       !< Copy to this variable
    logical, intent(in), optional  :: restrict    !< Restrict variable on multigrid tree
    logical, intent(in), optional  :: restrict_gc !< Fill ghost cells after restrict
    real(dp), intent(in), optional :: factor      !< out = factor * in
    !> If present, use this as the input state
    type(state), intent(in), optional, target :: state_from(max_blocks)
    integer                        :: iigrid, igrid, id
    integer                        :: nc, lvl, ilo^D, ihi^D
    type(tree_node), pointer       :: pnode
    real(dp)                       :: fac
    logical                        :: do_restrict, do_gc

    if (.not. mg%is_allocated) &
         error stop "mg_copy_to_tree: tree not allocated yet"

    fac = 1.0_dp; if (present(factor)) fac = factor
    do_restrict = .false.; if (present(restrict)) do_restrict = restrict
    do_gc = .false.; if (present(restrict_gc)) do_gc = restrict_gc

    ilo^D=ixMlo^D-1;
    ihi^D=ixMhi^D+1;

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Include one layer of ghost cells on grid leaves
       if (present(state_from)) then
          mg%boxes(id)%cc({0:nc+1}, iw_to) = &
               fac * state_from(igrid)%w({ilo^D:ihi^D}, iw_from)
       else
          mg%boxes(id)%cc({0:nc+1}, iw_to) = fac * ps(igrid)%w({ilo^D:ihi^D}, iw_from)
       end if
    end do

    if (do_restrict) then
       call mg_restrict(mg, iw_to)
       if (do_gc) call mg_fill_ghost_cells(mg, iw_to)
    end if

  end subroutine mg_copy_to_tree

  !> Copy a variable from the multigrid tree
  subroutine mg_copy_from_tree(iw_from, iw_to, state_to)
    use mod_global_parameters
    use mod_forest
    integer, intent(in)      :: iw_from !< Variable to use as right-hand side
    integer, intent(in)      :: iw_to   !< Copy to this variable
    !> If present, use this as the output state
    type(state), intent(inout), optional, target :: state_to(max_blocks)
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl
    type(tree_node), pointer :: pnode

    if (.not. mg%is_allocated) &
         error stop "mg_copy_from_tree: tree not allocated yet"

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       if (present(state_to)) then
          state_to(igrid)%w(ixM^T, iw_to) = mg%boxes(id)%cc({1:nc}, iw_from)
       else
          ps(igrid)%w(ixM^T, iw_to) = mg%boxes(id)%cc({1:nc}, iw_from)
       end if
    end do
  end subroutine mg_copy_from_tree

  !> Copy from multigrid tree with one layer of ghost cells. Corner ghost cells
  !> are not used/set.
  subroutine mg_copy_from_tree_gc(iw_from, iw_to, state_to)
    use mod_global_parameters
    use mod_forest
    integer, intent(in)      :: iw_from !< Variable to use as right-hand side
    integer, intent(in)      :: iw_to   !< Copy to this variable
    !> If present, use this as the output state
    type(state), intent(inout), optional, target :: state_to(max_blocks)
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl, ilo^D, ihi^D
    type(tree_node), pointer :: pnode

    if (.not. mg%is_allocated) &
         error stop "mg_copy_from_tree_gc: tree not allocated yet"

    ilo^D=ixMlo^D-1;
    ihi^D=ixMhi^D+1;

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       if (present(state_to)) then
          state_to(igrid)%w({ilo^D:ihi^D}, iw_to) = mg%boxes(id)%cc({0:nc+1}, iw_from)
       else
          ps(igrid)%w({ilo^D:ihi^D}, iw_to) = mg%boxes(id)%cc({0:nc+1}, iw_from)
       end if
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

             do i_c = 1, mg_num_children
                c_id = mg%boxes(id)%children(i_c)
                c_ix = mg_child_dix(:, i_c) + 1
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

    if (associated(mg_after_new_tree)) then
       call mg_after_new_tree()
    end if

  end subroutine mg_tree_from_amrvac

end module mod_multigrid_coupling

