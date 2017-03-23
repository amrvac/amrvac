!> \file
!> @todo Convert this file to a module

!> Initialize the MPI environment
!> @todo Check for errors in return code
subroutine comm_start
  use mod_global_parameters

  integer(kind=MPI_ADDRESS_KIND) :: lb
  integer(kind=MPI_ADDRESS_KIND) :: size

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
  call MPI_TYPE_GET_EXTENT(MPI_REAL,lb,size,ierrmpi)
  if (size /= size_real) call mpistop("Incompatible real size")
  call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size,ierrmpi)
  if (size /= size_double) call mpistop("Incompatible double size")
  call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,size,ierrmpi)
  if (size /= size_int) call mpistop("Incompatible integer size")
  call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,size,ierrmpi)
  if (size /= size_logical) call mpistop("Incompatible logical size")

end subroutine comm_start


!> Finalize (or shutdown) the MPI environment
subroutine comm_finalize

  use mod_global_parameters
  use mod_ghostcells_update
  use mod_particles, only: particles_finish

  if(use_particles) call particles_finish
  call put_bc_comm_types
  call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
  call MPI_FINALIZE(ierrmpi)

end subroutine comm_finalize


!> Create and store the MPI types that will be used for parallel communication
subroutine init_comm_types

use mod_global_parameters

integer, dimension(ndim+1) :: sizes, subsizes, start
!integer :: i^D, ic^D, nx^D, nxCo^D, size_double
integer :: i^D, ic^D, nx^D, nxCo^D, nxG^D
!-----------------------------------------------------------------------------
nx^D=ixMhi^D-ixMlo^D+1;
nxG^D=ixGhi^D-ixGlo^D+1;
nxCo^D=nx^D/2;

^D&sizes(^D)=ixGhi^D;
sizes(ndim+1)=nw
^D&subsizes(^D)=nxG^D;
subsizes(ndim+1)=nw
^D&start(^D)=ixGlo^D-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block,ierrmpi)
call MPI_TYPE_COMMIT(type_block,ierrmpi)
size_block={nxG^D*}*nw*size_double

^D&sizes(^D)=ixGhi^D/2+nghostcells;
sizes(ndim+1)=nw
^D&subsizes(^D)=nxCo^D;
subsizes(ndim+1)=nw
^D&start(^D)=ixMlo^D-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_coarse_block,ierrmpi)
call MPI_TYPE_COMMIT(type_coarse_block,ierrmpi)

^D&sizes(^D)=ixGhi^D;
sizes(ndim+1)=nw
{do ic^DB=1,2\}
   ^D&subsizes(^D)=nxCo^D;
   subsizes(ndim+1)=nw
   ^D&start(^D)=ixMlo^D-1+(ic^D-1)*nxCo^D;
   start(ndim+1)=0
   call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                                 MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                 type_sub_block(ic^D),ierrmpi)
   call MPI_TYPE_COMMIT(type_sub_block(ic^D),ierrmpi)
{end do\}

^D&sizes(^D)=ixGhi^D;
sizes(ndim+1)=nw
^D&subsizes(^D)=nx^D;
subsizes(ndim+1)=nw
^D&start(^D)=ixMlo^D-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_io,ierrmpi)
size_block_io={nx^D*}*nw*size_double

if(nwtf>0) then
  ^D&sizes(^D)=ixGhi^D;
  sizes(ndim+1)=nwtf
  ^D&subsizes(^D)=nx^D;
  subsizes(ndim+1)=nwtf
  ^D&start(^D)=ixMlo^D-1;
  start(ndim+1)=0
  call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                                MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                                type_block_io_tf,ierrmpi)
  call MPI_TYPE_COMMIT(type_block_io_tf,ierrmpi)
  
  size_block_io_tf={nx^D*}*nwtf*size_double
endif

^D&sizes(^D)=ixMhi^D-ixMlo^D+1;
sizes(ndim+1)=^ND
^D&subsizes(^D)=sizes(^D);
subsizes(ndim+1)=^ND
^D&start(^D)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_xcc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_xcc_io,ierrmpi)

^D&sizes(^D)=ixMhi^D-ixMlo^D+2;
sizes(ndim+1)=^ND
^D&subsizes(^D)=sizes(^D);
subsizes(ndim+1)=^ND
^D&start(^D)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_xc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_xc_io,ierrmpi)

^D&sizes(^D)=ixMhi^D-ixMlo^D+1;
sizes(ndim+1)=nw+nwauxio
^D&subsizes(^D)=sizes(^D);
subsizes(ndim+1)=nw+nwauxio
^D&start(^D)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_wcc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_wcc_io,ierrmpi)

^D&sizes(^D)=ixMhi^D-ixMlo^D+2;
sizes(ndim+1)=nw+nwauxio
^D&subsizes(^D)=sizes(^D);
subsizes(ndim+1)=nw+nwauxio
^D&start(^D)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
                              type_block_wc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_wc_io,ierrmpi)

end subroutine init_comm_types


!> Exit MPI-AMRVAC with an error message
subroutine mpistop(message)
  use mod_global_parameters

  character(len=*), intent(in) :: message !< The error message
  integer                      :: ierrcode

  write(*, *) "ERROR for processor", mype, ":"
  write(*, *) trim(message)

  call MPI_ABORT(icomm, ierrcode, ierrmpi)

end subroutine mpistop

