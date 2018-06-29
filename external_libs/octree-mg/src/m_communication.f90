module m_communication
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: mg_comm_init
  public :: sort_and_transfer_buffers

contains

  !> Initialize MPI if needed, and store MPI information
  subroutine mg_comm_init(mg, comm)
    use mpi
    type(mg_t), intent(inout)     :: mg
    !> MPI communicator (default: MPI_COMM_WORLD)
    integer, intent(in), optional :: comm
    integer                       :: ierr
    logical                       :: initialized

    call mpi_initialized(initialized, ierr)
    if (.not. initialized) then
       call mpi_init(ierr)
    end if

    if (present(comm)) then
       mg%comm = comm
    else
       mg%comm = MPI_COMM_WORLD
    end if

    call mpi_comm_rank(mg%comm, mg%my_rank, ierr)
    call mpi_comm_size(mg%comm, mg%n_cpu, ierr)
  end subroutine mg_comm_init

  subroutine sort_and_transfer_buffers(mg, dsize)
    use mpi
    type(mg_t), intent(inout)    :: mg
    integer, intent(in)          :: dsize
    integer                      :: i, n_send, n_recv
    integer                      :: send_req(mg%n_cpu)
    integer                      :: recv_req(mg%n_cpu)
    integer                      :: ierr

    n_send = 0
    n_recv = 0

    do i = 0, mg%n_cpu - 1
       if (mg%buf(i)%i_send > 0) then
          n_send = n_send + 1
          call sort_sendbuf(mg%buf(i), dsize)
          call mpi_isend(mg%buf(i)%send, mg%buf(i)%i_send, MPI_DOUBLE, &
               i, 0, mg%comm, send_req(n_send), ierr)
       end if
       if (mg%buf(i)%i_recv > 0) then
          n_recv = n_recv + 1
          call mpi_irecv(mg%buf(i)%recv, mg%buf(i)%i_recv, MPI_DOUBLE, &
               i, 0, mg%comm, recv_req(n_recv), ierr)
       end if
    end do

    call mpi_waitall(n_recv, recv_req(1:n_recv), MPI_STATUSES_IGNORE, ierr)
    call mpi_waitall(n_send, send_req(1:n_send), MPI_STATUSES_IGNORE, ierr)

  end subroutine sort_and_transfer_buffers

  !> Sort send buffers according to the idbuf array
  subroutine sort_sendbuf(gc, dsize)
    use m_mrgrnk
    type(mg_buf_t), intent(inout) :: gc
    integer, intent(in)           :: dsize !< Size of send buffer elements
    integer                       :: ix_sort(gc%i_ix)
    real(dp)                      :: buf_cpy(gc%i_send)
    integer                       :: i, j, n

    call mrgrnk(gc%ix(1:gc%i_ix), ix_sort)

    buf_cpy = gc%send(1:gc%i_send)

    do n = 1, gc%i_ix
       i = (n-1) * dsize
       j = (ix_sort(n)-1) * dsize
       gc%send(i+1:i+dsize) = buf_cpy(j+1:j+dsize)
    end do
    gc%ix(1:gc%i_ix) = gc%ix(ix_sort)

  end subroutine sort_sendbuf

end module m_communication
