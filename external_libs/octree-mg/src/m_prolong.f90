#include "cpp_macros.h"
module m_prolong
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: mg_prolong
  public :: mg_prolong_buffer_size
  public :: mg_prolong_sparse

contains

  !> Specify minimum buffer size (per process) for communication
  subroutine mg_prolong_buffer_size(mg, n_send, n_recv, dsize)
    type(mg_t), intent(inout) :: mg
    integer, intent(out)      :: n_send(0:mg%n_cpu-1)
    integer, intent(out)      :: n_recv(0:mg%n_cpu-1)
    integer, intent(out)      :: dsize
    integer                   :: lvl, min_lvl

    if (.not. allocated(mg%comm_restrict%n_send)) then
       error stop "Call restrict_buffer_size before prolong_buffer_size"
    end if

    min_lvl = max(mg%first_normal_lvl-1, mg%lowest_lvl)
    allocate(mg%comm_prolong%n_send(0:mg%n_cpu-1, &
         min_lvl:mg%highest_lvl))
    allocate(mg%comm_prolong%n_recv(0:mg%n_cpu-1, &
         min_lvl:mg%highest_lvl))

    mg%comm_prolong%n_recv(:, mg%highest_lvl) = 0
    mg%comm_prolong%n_send(:, mg%highest_lvl) = 0

    do lvl = min_lvl, mg%highest_lvl-1
       mg%comm_prolong%n_recv(:, lvl) = &
            mg%comm_restrict%n_send(:, lvl+1)
       mg%comm_prolong%n_send(:, lvl) = &
            mg%comm_restrict%n_recv(:, lvl+1)
    end do

    ! Send fine grid points, because this is more flexible than sending coarse
    ! grid points (e.g., when multiple variables are used for interpolation)
    dsize = (mg%box_size)**NDIM
    n_send = maxval(mg%comm_prolong%n_send, dim=2)
    n_recv = maxval(mg%comm_prolong%n_recv, dim=2)
  end subroutine mg_prolong_buffer_size

  !> Prolong variable iv from lvl to variable iv_to at lvl+1
  subroutine mg_prolong(mg, lvl, iv, iv_to, method, add)
    use m_communication
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl   !< Level to prolong from
    integer, intent(in)       :: iv    !< Source variable
    integer, intent(in)       :: iv_to !< Target variable
    procedure(mg_box_prolong) :: method !< Prolongation method
    logical, intent(in)       :: add   !< If true, add to current values
    integer                   :: i, id, dsize, nc

    if (lvl == mg%highest_lvl) error stop "cannot prolong highest level"
    if (lvl < mg%lowest_lvl) error stop "cannot prolong below lowest level"

    ! Below the first normal level, all boxes are on the same CPU
    if (lvl >= mg%first_normal_lvl-1) then
       dsize            = mg%box_size**NDIM
       mg%buf(:)%i_send = 0
       mg%buf(:)%i_ix   = 0

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call prolong_set_buffer(mg, id, mg%box_size, iv, method)
       end do

       mg%buf(:)%i_recv = mg%comm_prolong%n_recv(:, lvl) * dsize
       call sort_and_transfer_buffers(mg, dsize)
       mg%buf(:)%i_recv = 0
    end if

    nc = mg%box_size_lvl(lvl+1)
    do i = 1, size(mg%lvls(lvl+1)%my_ids)
       id = mg%lvls(lvl+1)%my_ids(i)
       call prolong_onto(mg, id, nc, iv, iv_to, add, method)
    end do
  end subroutine mg_prolong

  !> In case the fine grid is on a different CPU, perform the prolongation and
  !> store the fine-grid values in the send buffer.
  subroutine prolong_set_buffer(mg, id, nc, iv, method)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    procedure(mg_box_prolong) :: method
    integer                   :: i, dix(NDIM)
    integer                   :: i_c, c_id, c_rank, dsize
    real(dp)                  :: tmp(DTIMES(nc))

    dsize = nc**NDIM

    do i_c = 1, mg_num_children
       c_id = mg%boxes(id)%children(i_c)
       if (c_id > mg_no_box) then
          c_rank = mg%boxes(c_id)%rank
          if (c_rank /= mg%my_rank) then
             dix = mg_get_child_offset(mg, c_id)
             call method(mg, id, dix, nc, iv, tmp)

             i   = mg%buf(c_rank)%i_send
             mg%buf(c_rank)%send(i+1:i+dsize) = pack(tmp, .true.)
             mg%buf(c_rank)%i_send  = mg%buf(c_rank)%i_send + dsize

             i                      = mg%buf(c_rank)%i_ix
             mg%buf(c_rank)%ix(i+1) = c_id
             mg%buf(c_rank)%i_ix    = mg%buf(c_rank)%i_ix + 1
          end if
       end if
    end do
  end subroutine prolong_set_buffer

  !> Prolong onto a child box
  subroutine prolong_onto(mg, id, nc, iv, iv_to, add, method)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv    !< Prolong from this variable
    integer, intent(in)       :: iv_to !< Prolong to this variable
    logical, intent(in)       :: add   !< If true, add to current values
    procedure(mg_box_prolong) :: method
    integer                   :: hnc, p_id, p_rank, i, dix(NDIM), dsize
    real(dp)                  :: tmp(DTIMES(nc))

    hnc    = nc/2
    p_id   = mg%boxes(id)%parent
    p_rank = mg%boxes(p_id)%rank

    if (p_rank == mg%my_rank) then
       dix    = mg_get_child_offset(mg, id)
       call method(mg, p_id, dix, nc, iv, tmp)
    else
       dsize  = nc**NDIM
       i = mg%buf(p_rank)%i_recv
       tmp = reshape(mg%buf(p_rank)%recv(i+1:i+dsize), [DTIMES(nc)])
       mg%buf(p_rank)%i_recv = mg%buf(p_rank)%i_recv + dsize
    end if

    if (add) then
       mg%boxes(id)%cc(DTIMES(1:nc), iv_to) = &
            mg%boxes(id)%cc(DTIMES(1:nc), iv_to) + tmp
    else
       mg%boxes(id)%cc(DTIMES(1:nc), iv_to) = tmp
    end if

  end subroutine prolong_onto

  !> Prolong from a parent to a child with index offset dix
  subroutine mg_prolong_sparse(mg, p_id, dix, nc, iv, fine)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: p_id             !< Id of parent
    integer, intent(in)       :: dix(NDIM)        !< Offset of child in parent grid
    integer, intent(in)       :: nc               !< Child grid size
    integer, intent(in)       :: iv               !< Prolong from this variable
    real(dp), intent(out)     :: fine(DTIMES(nc)) !< Prolonged values

    integer  :: IJK, hnc
#if NDIM == 2
    integer  :: ic, jc
    real(dp) :: f0, flx, fhx, fly, fhy
#elif NDIM == 3
    integer  :: ic, jc, kc
    real(dp) :: f0, flx, fhx, fly, fhy, flz, fhz
#endif

    hnc = nc/2

    associate (crs => mg%boxes(p_id)%cc)
#if NDIM == 2
      do j = 1, hnc
         jc = j + dix(2)
         do i = 1, hnc
            ic = i + dix(1)

            f0  = 0.5_dp * crs(ic, jc, iv)
            flx = 0.25_dp * crs(ic-1, jc, iv)
            fhx = 0.25_dp * crs(ic+1, jc, iv)
            fly = 0.25_dp * crs(ic, jc-1, iv)
            fhy = 0.25_dp * crs(ic, jc+1, iv)

            fine(2*i-1, 2*j-1) = f0 + flx + fly
            fine(2*i  , 2*j-1) = f0 + fhx + fly
            fine(2*i-1, 2*j)   = f0 + flx + fhy
            fine(2*i  , 2*j)   = f0 + fhx + fhy
         end do
      end do
#elif NDIM == 3
      do k = 1, hnc
         kc = k + dix(3)
         do j = 1, hnc
            jc = j + dix(2)
            do i = 1, hnc
               ic = i + dix(1)

               f0  = 0.25_dp * crs(ic, jc, kc, iv)
               flx = 0.25_dp * crs(ic-1, jc, kc, iv)
               fhx = 0.25_dp * crs(ic+1, jc, kc, iv)
               fly = 0.25_dp * crs(ic, jc-1, kc, iv)
               fhy = 0.25_dp * crs(ic, jc+1, kc, iv)
               flz = 0.25_dp * crs(ic, jc, kc-1, iv)
               fhz = 0.25_dp * crs(ic, jc, kc+1, iv)

               fine(2*i-1, 2*j-1, 2*k-1) = f0 + flx + fly + flz
               fine(2*i, 2*j-1, 2*k-1)   = f0 + fhx + fly + flz
               fine(2*i-1, 2*j, 2*k-1)   = f0 + flx + fhy + flz
               fine(2*i, 2*j, 2*k-1)     = f0 + fhx + fhy + flz
               fine(2*i-1, 2*j-1, 2*k)   = f0 + flx + fly + fhz
               fine(2*i, 2*j-1, 2*k)     = f0 + fhx + fly + fhz
               fine(2*i-1, 2*j, 2*k)     = f0 + flx + fhy + fhz
               fine(2*i, 2*j, 2*k)       = f0 + fhx + fhy + fhz
            end do
         end do
      end do
#endif
    end associate
  end subroutine mg_prolong_sparse

end module m_prolong
