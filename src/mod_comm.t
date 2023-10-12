module mod_comm
  use mpi

  implicit none
  public

  !> The number of MPI tasks
  integer :: npe

  !> The rank of the current MPI task
  integer :: mype

  !> The MPI communicator
  integer :: icomm
  
  contains
  
  !> Initialize the MPI environment
  !> @todo Check for errors in return code,
  !> maybe add a err checker
  subroutine comm_start
    use mod_basic_types
  
    integer(kind=MPI_ADDRESS_KIND) :: lb
    integer(kind=MPI_ADDRESS_KIND) :: sizes
    integer                        :: ierrmpi !> @todo Check for errors in return code
  
    ! Initialize MPI
    call MPI_INIT(ierrmpi)
  
    ! Each process stores its rank, which ranges from 0 to N-1, where N is the
    ! number of processes.
    call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierrmpi)
  
    ! Store the number of processes
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierrmpi)
  
    ! Use the default communicator, which contains all the processes
    icomm = MPI_COMM_WORLD
  
    ! Get size of double/integer
    call MPI_TYPE_GET_EXTENT(MPI_REAL,lb,sizes,ierrmpi)
    if (sizes /= size_real) call mpistop("Incompatible real size")
    call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,sizes,ierrmpi)
    if (sizes /= size_double) call mpistop("Incompatible double size")
    call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,sizes,ierrmpi)
    if (sizes /= size_int) call mpistop("Incompatible integer size")
    call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,sizes,ierrmpi)
    if (sizes /= size_logical) call mpistop("Incompatible logical size")
  
  end subroutine comm_start
  
  !> Finalize (or shutdown) the MPI environment
  subroutine comm_finalize
    integer :: ierrmpi
    call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
    call MPI_FINALIZE(ierrmpi)
  end subroutine comm_finalize
  
  !> Exit GMUNU with an error message
  subroutine mpistop(message)
    character(len=*), intent(in) :: message !< The error message
    integer                      :: ierrcode, ierrmpi
    write(*, *) "ERROR for processor", mype, ":"
    write(*, *) trim(message)
    call MPI_ABORT(icomm, ierrcode, ierrmpi)
  end subroutine mpistop

end module mod_comm
