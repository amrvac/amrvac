#include "cpp_macros.h"
module m_restrict
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: mg_restrict
  public :: mg_restrict_lvl
  public :: mg_restrict_buffer_size

contains

  !> Specify minimum buffer size (per process) for communication
  subroutine mg_restrict_buffer_size(mg, n_send, n_recv, dsize)
    type(mg_t), intent(inout) :: mg
    integer, intent(out)      :: n_send(0:mg%n_cpu-1)
    integer, intent(out)      :: n_recv(0:mg%n_cpu-1)
    integer, intent(out)      :: dsize
    integer :: n_out(0:mg%n_cpu-1, mg%first_normal_lvl:mg%highest_lvl)
    integer :: n_in(0:mg%n_cpu-1, mg%first_normal_lvl:mg%highest_lvl)
    integer                   :: lvl, i, id, p_id, p_rank
    integer                   :: i_c, c_id, c_rank, min_lvl

    n_out(:, :) = 0
    n_in(:, :)  = 0
    min_lvl = max(mg%lowest_lvl+1, mg%first_normal_lvl)

    do lvl = min_lvl, mg%highest_lvl
       ! Number of messages to receive (at lvl-1)
       do i = 1, size(mg%lvls(lvl-1)%my_parents)
          id = mg%lvls(lvl-1)%my_parents(i)
          do i_c = 1, mg_num_children
             c_id = mg%boxes(id)%children(i_c)

             if (c_id > mg_no_box) then
                c_rank = mg%boxes(c_id)%rank
                if (c_rank /= mg%my_rank) then
                   n_in(c_rank, lvl) = n_in(c_rank, lvl) + 1
                end if
             end if
          end do
       end do

       ! Number of messages to send
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)

          p_id = mg%boxes(id)%parent
          p_rank = mg%boxes(p_id)%rank
          if (p_rank /= mg%my_rank) then
             n_out(p_rank, lvl) = n_out(p_rank, lvl) + 1
          end if

       end do
    end do

    allocate(mg%comm_restrict%n_send(0:mg%n_cpu-1, &
         mg%first_normal_lvl:mg%highest_lvl))
    allocate(mg%comm_restrict%n_recv(0:mg%n_cpu-1, &
         mg%first_normal_lvl:mg%highest_lvl))
    mg%comm_restrict%n_send = n_out
    mg%comm_restrict%n_recv = n_in

    dsize  = (mg%box_size/2)**NDIM
    n_send = maxval(n_out, dim=2)
    n_recv = maxval(n_in, dim=2)
  end subroutine mg_restrict_buffer_size

  !> Restrict all levels
  subroutine mg_restrict(mg, iv)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: iv
    integer                   :: lvl

    do lvl = mg%highest_lvl, mg%lowest_lvl+1, -1
       call mg_restrict_lvl(mg, iv, lvl)
    end do
  end subroutine mg_restrict

  !> Restrict from lvl to lvl-1
  subroutine mg_restrict_lvl(mg, iv, lvl)
    use m_communication
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: iv
    integer, intent(in)       :: lvl
    integer                   :: i, id, dsize, nc

    if (lvl <= mg%lowest_lvl) error stop "cannot restrict lvl <= lowest_lvl"

    nc = mg%box_size_lvl(lvl)

    if (lvl >= mg%first_normal_lvl) then
       dsize = (nc/2)**NDIM

       mg%buf(:)%i_send = 0
       mg%buf(:)%i_ix   = 0

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call restrict_set_buffer(mg, id, iv)
       end do

       mg%buf(:)%i_recv = mg%comm_restrict%n_recv(:, lvl) * dsize
       call sort_and_transfer_buffers(mg, dsize)
       mg%buf(:)%i_recv = 0
    end if

    do i = 1, size(mg%lvls(lvl-1)%my_parents)
       id = mg%lvls(lvl-1)%my_parents(i)
       call restrict_onto(mg, id, nc, iv)
    end do
  end subroutine mg_restrict_lvl

  subroutine restrict_set_buffer(mg, id, iv)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer, intent(in)          :: iv
    integer                      :: IJK, n, hnc, p_id, p_rank
    real(dp) :: tmp(DTIMES(mg%box_size/2))

    hnc    = mg%box_size/2
    p_id   = mg%boxes(id)%parent
    p_rank = mg%boxes(p_id)%rank

    if (p_rank /= mg%my_rank) then
#if NDIM == 2
       do j = 1, hnc
          do i = 1, hnc
             tmp(i, j) = 0.25_dp * &
                  sum(mg%boxes(id)%cc(2*i-1:2*i, 2*j-1:2*j, iv))
          end do
       end do
#elif NDIM == 3
       do k = 1, hnc
          do j = 1, hnc
             do i = 1, hnc
                tmp(i, j, k) = 0.125_dp * sum(mg%boxes(id)%cc(2*i-1:2*i, &
                     2*j-1:2*j, 2*k-1:2*k, iv))
             end do
          end do
       end do
#endif

       ! Buffer
       n = size(tmp)
       i = mg%buf(p_rank)%i_send
       mg%buf(p_rank)%send(i+1:i+n) = pack(tmp, .true.)
       mg%buf(p_rank)%i_send = mg%buf(p_rank)%i_send + n

       ! To later sort the send buffer according to parent order
       i = mg%buf(p_rank)%i_ix
       n = mg_ix_to_ichild(mg%boxes(id)%ix)
       mg%buf(p_rank)%ix(i+1) = mg_num_children * p_id + n
       mg%buf(p_rank)%i_ix = mg%buf(p_rank)%i_ix + 1
    end if
  end subroutine restrict_set_buffer

  subroutine restrict_onto(mg, id, nc, iv)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    integer                   :: IJK, hnc, dsize, i_c, c_id
    integer                   :: c_rank, dix(NDIM)

    hnc   = nc/2
    dsize = hnc**NDIM

    do i_c = 1, mg_num_children
       c_id   = mg%boxes(id)%children(i_c)
       if (c_id == mg_no_box) cycle ! For coarsened grid
       c_rank = mg%boxes(c_id)%rank
       dix    = mg_get_child_offset(mg, c_id)

       if (c_rank == mg%my_rank) then
          do KJI_DO(1, hnc)
#if NDIM == 2
             mg%boxes(id)%cc(dix(1)+i, dix(2)+j, iv) = 0.25_dp * &
                  sum(mg%boxes(c_id)%cc(2*i-1:2*i, 2*j-1:2*j, iv))
#elif NDIM == 3
             mg%boxes(id)%cc(dix(1)+i, dix(2)+j, dix(3)+k, iv) = &
                  0.125_dp * sum(mg%boxes(c_id)%cc(2*i-1:2*i, &
                  2*j-1:2*j, 2*k-1:2*k, iv))
#endif
          end do; CLOSE_DO
       else
          i = mg%buf(c_rank)%i_recv
#if NDIM == 2
          mg%boxes(id)%cc(dix(1)+1:dix(1)+hnc, &
               dix(2)+1:dix(2)+hnc, iv) = &
               reshape(mg%buf(c_rank)%recv(i+1:i+dsize), [hnc, hnc])
#elif NDIM == 3
          mg%boxes(id)%cc(dix(1)+1:dix(1)+hnc, &
               dix(2)+1:dix(2)+hnc, dix(3)+1:dix(3)+hnc, iv) = &
               reshape(mg%buf(c_rank)%recv(i+1:i+dsize), [hnc, hnc, hnc])
#endif
          mg%buf(c_rank)%i_recv = mg%buf(c_rank)%i_recv + dsize
       end if
    end do

  end subroutine restrict_onto

end module m_restrict
