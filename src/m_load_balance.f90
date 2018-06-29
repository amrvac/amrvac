!> Module for load balancing a tree (that already has been constructed). The
!> load balancing determines which ranks (MPI processes) allocated physical
!> storage for boxes. The tree structure itself is present on all processes.
module m_load_balance
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: mg_load_balance
  public :: mg_load_balance_parents

contains

  !> Load balance all boxes in the multigrid tree. Note that in a typical
  !> application the load balancing of the leaves is already determined, then
  !> mg_load_balance_parents can be used.
  subroutine mg_load_balance(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: i, id, lvl, single_cpu_lvl
    integer                   :: work_left, my_work, i_cpu

    ! Up to this level, all boxes have to be on a single processor because they
    ! have a different size and the communication routines don't support this
    single_cpu_lvl = max(mg%first_normal_lvl-1, mg%lowest_lvl)

    do lvl = mg%lowest_lvl, single_cpu_lvl
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          mg%boxes(id)%rank = 0
       end do
    end do

    ! Distribute the boxes equally. Due to the way the mesh is constructed, the
    ! mg%lvls(lvl)%ids array already contains a Morton-like ordering.
    do lvl = single_cpu_lvl+1, mg%highest_lvl
       work_left = size(mg%lvls(lvl)%ids)
       my_work   = 0
       i_cpu     = 0

       do i = 1, size(mg%lvls(lvl)%ids)
          if ((mg%n_cpu - i_cpu - 1) * my_work >= work_left) then
             i_cpu   = i_cpu + 1
             my_work = 0
          end if

          my_work = my_work + 1
          work_left = work_left - 1

          id = mg%lvls(lvl)%ids(i)
          mg%boxes(id)%rank = i_cpu
       end do
    end do

    do lvl = mg%lowest_lvl, mg%highest_lvl
       call update_lvl_info(mg, mg%lvls(lvl))
    end do

  end subroutine mg_load_balance

  !> Load balance the parents (non-leafs). Assign them to the rank that has most
  !> children.
  !>
  !> @todo In case of ties, assign to the rank with least boxes.
  subroutine mg_load_balance_parents(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: i, id, lvl
    integer                   :: c_ids(mg_num_children)
    integer                   :: c_ranks(mg_num_children)
    integer                   :: single_cpu_lvl, coarse_rank

    ! Up to this level, all boxes have to be on a single processor because they
    ! have a different size and the communication routines don't support this
    single_cpu_lvl = max(mg%first_normal_lvl-1, mg%lowest_lvl)

    do lvl = mg%highest_lvl-1, single_cpu_lvl+1, -1
       do i = 1, size(mg%lvls(lvl)%parents)
          id = mg%lvls(lvl)%parents(i)

          c_ids = mg%boxes(id)%children
          c_ranks = mg%boxes(c_ids)%rank
          mg%boxes(id)%rank = most_popular(c_ranks)
       end do
    end do

    ! Determine most popular CPU for coarse grids
    coarse_rank = most_popular(mg%boxes(&
         mg%lvls(single_cpu_lvl+1)%ids)%rank)

    do lvl = mg%lowest_lvl, single_cpu_lvl
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          mg%boxes(id)%rank = coarse_rank
       end do
    end do

    do lvl = mg%lowest_lvl, mg%highest_lvl
       call update_lvl_info(mg, mg%lvls(lvl))
    end do

  end subroutine mg_load_balance_parents

  pure integer function most_popular(list)
    integer, intent(in) :: list(:)
    integer             :: i, best_count, current_count

    best_count   = 0
    most_popular = -1

    do i = 1, size(list)
       current_count = count(list == list(i))

       if (current_count > best_count) then
          most_popular = list(i)
       end if
    end do

  end function most_popular

  subroutine update_lvl_info(mg, lvl)
    type(mg_t), intent(inout)     :: mg
    type(mg_lvl_t), intent(inout) :: lvl

    lvl%my_ids = pack(lvl%ids, &
         mg%boxes(lvl%ids)%rank == mg%my_rank)
    lvl%my_leaves = pack(lvl%leaves, &
         mg%boxes(lvl%leaves)%rank == mg%my_rank)
    lvl%my_parents = pack(lvl%parents, &
         mg%boxes(lvl%parents)%rank == mg%my_rank)
    lvl%my_ref_bnds = pack(lvl%ref_bnds, &
         mg%boxes(lvl%ref_bnds)%rank == mg%my_rank)
  end subroutine update_lvl_info

end module m_load_balance
