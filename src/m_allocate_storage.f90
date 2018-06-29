#include "cpp_macros.h"
module m_allocate_storage
  use m_data_structures

  implicit none
  private

  public :: mg_allocate_storage
  public :: mg_deallocate_storage

contains

  !> Deallocate all allocatable arrays
  subroutine mg_deallocate_storage(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: lvl

    if (.not. mg%is_allocated) &
         error stop "deallocate_storage: tree is not allocated"

    deallocate(mg%boxes)
    deallocate(mg%buf)

    deallocate(mg%comm_restrict%n_send)
    deallocate(mg%comm_restrict%n_recv)

    deallocate(mg%comm_prolong%n_send)
    deallocate(mg%comm_prolong%n_recv)

    deallocate(mg%comm_ghostcell%n_send)
    deallocate(mg%comm_ghostcell%n_recv)

    do lvl = mg%lowest_lvl, mg%highest_lvl
       deallocate(mg%lvls(lvl)%ids)
       deallocate(mg%lvls(lvl)%leaves)
       deallocate(mg%lvls(lvl)%parents)
       deallocate(mg%lvls(lvl)%ref_bnds)
       deallocate(mg%lvls(lvl)%my_ids)
       deallocate(mg%lvls(lvl)%my_leaves)
       deallocate(mg%lvls(lvl)%my_parents)
       deallocate(mg%lvls(lvl)%my_ref_bnds)
    end do

    mg%is_allocated       = .false.
    mg%n_boxes            = 0
    mg%phi_bc_data_stored = .false.
  end subroutine mg_deallocate_storage

  !> Allocate communication buffers and local boxes for a tree that has already
  !> been created
  subroutine mg_allocate_storage(mg)
    use m_ghost_cells, only: mg_ghost_cell_buffer_size
    use m_restrict, only: mg_restrict_buffer_size
    use m_prolong, only: mg_prolong_buffer_size
    type(mg_t), intent(inout) :: mg
    integer                   :: i, id, lvl, nc
    integer                   :: n_send(0:mg%n_cpu-1, 3)
    integer                   :: n_recv(0:mg%n_cpu-1, 3)
    integer                   :: dsize(3)
    integer                   :: n_in, n_out, n_id

    if (.not. mg%tree_created) &
         error stop "allocate_storage: tree is not yet created"

    if (mg%is_allocated) &
         error stop "allocate_storage: tree is already allocated"

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          allocate(mg%boxes(id)%cc(DTIMES(0:nc+1), &
               mg_num_vars + mg%n_extra_vars))

          ! Set all initial values to zero
          mg%boxes(id)%cc(DTIMES(:), :) = 0.0_dp
       end do
    end do

    allocate(mg%buf(0:mg%n_cpu-1))

    call mg_ghost_cell_buffer_size(mg, n_send(:, 1), &
         n_recv(:, 1), dsize(1))
    call mg_restrict_buffer_size(mg, n_send(:, 2), &
         n_recv(:, 2), dsize(2))
    call mg_prolong_buffer_size(mg, n_send(:, 3), &
         n_recv(:, 3), dsize(3))

    do i = 0, mg%n_cpu-1
       n_out = maxval(n_send(i, :) * dsize(:))
       n_in = maxval(n_recv(i, :) * dsize(:))
       n_id = maxval(n_send(i, :))
       allocate(mg%buf(i)%send(n_out))
       allocate(mg%buf(i)%recv(n_in))
       allocate(mg%buf(i)%ix(n_id))
    end do

    mg%is_allocated = .true.
  end subroutine mg_allocate_storage

end module m_allocate_storage
