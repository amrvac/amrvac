!=============================================================================
subroutine comm_start

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
call MPI_INIT(ierrmpi)
call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierrmpi)
call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierrmpi)

icomm=MPI_COMM_WORLD

end subroutine comm_start
!=============================================================================
subroutine comm_finalize

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
call MPI_FINALIZE(ierrmpi)

end subroutine comm_finalize
!=============================================================================
subroutine init_comm_types

include 'amrvacdef.f'

integer, dimension(ndim+1) :: sizes, subsizes, start
!integer :: i^D, ic^D, nx^D, nxCo^D, size_double
integer :: i1, ic1, nx1, nxCo1

 integer(kind=MPI_ADDRESS_KIND):: size_double, lb
!-----------------------------------------------------------------------------
nx1=ixMhi1-ixMlo1+1;
nxCo1=nx1/2;

sizes(1)=ixGhi1;
sizes(ndim+1)=nw
subsizes(1)=nx1;
subsizes(ndim+1)=nwflux
start(1)=ixMlo1-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block,ierrmpi)
call MPI_TYPE_COMMIT(type_block,ierrmpi)

sizes(1)=ixGhi1/2+dixB;
sizes(ndim+1)=nw
subsizes(1)=nxCo1;
subsizes(ndim+1)=nw
start(1)=ixMlo1-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_coarse_block,ierrmpi)
call MPI_TYPE_COMMIT(type_coarse_block,ierrmpi)

sizes(1)=ixGhi1;
sizes(ndim+1)=nw
do ic1=1,2
   subsizes(1)=nxCo1;
   subsizes(ndim+1)=nw
   start(1)=ixMlo1-1+(ic1-1)*nxCo1;
   start(ndim+1)=0
   call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_sub_block(ic1),ierrmpi)
   call MPI_TYPE_COMMIT(type_sub_block(ic1),ierrmpi)
end do

sizes(1)=ixGhi1;
sizes(ndim+1)=nw
subsizes(1)=nx1;
subsizes(ndim+1)=nw
start(1)=ixMlo1-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_io,ierrmpi)

!call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,size_double,ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
size_block_io=nx1*nw*size_double



sizes(1)=ixMhi1-ixMlo1+1;
sizes(ndim+1)=1
subsizes(1)=sizes(1);
subsizes(ndim+1)=1
start(1)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_xcc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_xcc_io,ierrmpi)

sizes(1)=ixMhi1-ixMlo1+2;
sizes(ndim+1)=1
subsizes(1)=sizes(1);
subsizes(ndim+1)=1
start(1)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_xc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_xc_io,ierrmpi)

sizes(1)=ixMhi1-ixMlo1+1;
sizes(ndim+1)=nw+nwauxio
subsizes(1)=sizes(1);
subsizes(ndim+1)=nw+nwauxio
start(1)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_wcc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_wcc_io,ierrmpi)

sizes(1)=ixMhi1-ixMlo1+2;
sizes(ndim+1)=nw+nwauxio
subsizes(1)=sizes(1);
subsizes(ndim+1)=nw+nwauxio
start(1)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_wc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_wc_io,ierrmpi)

end subroutine init_comm_types
!=============================================================================
subroutine mpistop(message)

include 'amrvacdef.f'

character(len=*), intent(in) :: message

integer :: ierrcode
!------------------------------------------------------------------------------
write(*,*) "ERROR for processor",mype,":"
write(*,*) message
call MPI_ABORT(icomm,ierrcode,ierrmpi)

end subroutine mpistop
!==============================================================================
