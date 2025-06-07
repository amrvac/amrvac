module mod_comm_lib

  implicit none
  private


  public :: comm_start
  public :: comm_finalize
  public :: init_comm_types
  public :: mpistop
  public :: mpistop_gpu

contains


  !> Initialize the MPI environment
  subroutine comm_start
    use mod_global_parameters
  
    integer(kind=MPI_ADDRESS_KIND) :: lb
    integer(kind=MPI_ADDRESS_KIND) :: sizes
  
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

#ifdef _OPENACC
    call set_openacc_device
#endif
    
  end subroutine comm_start
  
#ifdef _OPENACC
  !> Set the device to be used by OpenACC. Code based on:
  !> https://docs.nvidia.com/hpc-sdk/compilers/openacc-mpi-tutorial/
  !> OpenMPI provides an environment variable, but this is not portable
  !> Taken from Jannis Teunissen foap4
  subroutine set_openacc_device
    use mod_global_parameters
    use openacc

    ! f4%mpisize -> npe
    ! f4%mpirank -> mype
    ! f4%mpicomm -> MPI_COMM_WORLD

    interface
       ! Get a unique number to identify the host
       function gethostid() bind(C)
         import C_int
         integer (C_int) :: gethostid
       end function gethostid
    end interface

    integer :: hostids(0:npe-1), local_procs(0:npe-1)
    integer :: hostid, ierr, num_devices, my_device, rank, num_local_procs
    integer(acc_device_kind) :: dev_type

    dev_type = ACC_DEVICE_DEFAULT

    ! Get the hostids to determine how many processes are on this host
    hostid = gethostid()
    call MPI_Allgather(hostid, 1, MPI_INTEGER, hostids, 1, MPI_INTEGER, &
         MPI_COMM_WORLD, ierr)

    ! Determine the local MPI ranks and number them, starting at zero
    num_local_procs = 0
    local_procs     = 0

    do rank = 0, npe-1
       if (hostid == hostids(rank)) then
          local_procs(rank) = num_local_procs
          num_local_procs = num_local_procs+1
       endif
    enddo

    num_devices = acc_get_num_devices(dev_type)
    if (num_devices < 1) error stop "No devices available on host"

    if (num_devices < num_local_procs) then
       ! Print warning only for first local process
       if (local_procs(mype) == 0) then
          write(*, "(A,I0,A,I0,A,I0,A)") "WARNING from ", mype, &
               ": more local processes (", num_local_procs, &
               ") than GPUs (", num_devices, ")"
       endif

       my_device = mod(local_procs(mype), num_devices)
    else
       my_device = local_procs(mype)
    endif

    call acc_set_device_num(my_device, dev_type)

  end subroutine set_openacc_device
#endif
  
  !> Finalize (or shutdown) the MPI environment
  subroutine comm_finalize
    use mod_global_parameters
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
    call MPI_FINALIZE(ierrmpi)
  
  end subroutine comm_finalize
  
  !> Create and store the MPI types that will be used for parallel communication
  subroutine init_comm_types
    use mod_global_parameters
  
    integer, dimension(ndim+1) :: sizes, subsizes, start
    integer :: i1,i2,i3, ic1,ic2,ic3, nx1,nx2,nx3, nxCo1,nxCo2,nxCo3, nxG1,&
       nxG2,nxG3, idir
  
    nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
    nxG1=ixGhi1-ixGlo1+1;nxG2=ixGhi2-ixGlo2+1;nxG3=ixGhi3-ixGlo3+1;
    nxCo1=nx1/2;nxCo2=nx2/2;nxCo3=nx3/2;
  
    sizes(1)=ixGhi1;sizes(2)=ixGhi2;sizes(3)=ixGhi3;
    sizes(ndim+1)=nw
    subsizes(1)=nxG1;subsizes(2)=nxG2;subsizes(3)=nxG3;
    subsizes(ndim+1)=nw
    start(1)=ixGlo1-1;start(2)=ixGlo2-1;start(3)=ixGlo3-1;
    start(ndim+1)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_block,ierrmpi)
    call MPI_TYPE_COMMIT(type_block,ierrmpi)
    size_block=nxG1*nxG2*nxG3*nw*size_double
  
    sizes(1)=ixGhi1/2+nghostcells;sizes(2)=ixGhi2/2+nghostcells
    sizes(3)=ixGhi3/2+nghostcells;
    sizes(ndim+1)=nw
    subsizes(1)=nxCo1;subsizes(2)=nxCo2;subsizes(3)=nxCo3;
    subsizes(ndim+1)=nw
    start(1)=ixMlo1-1;start(2)=ixMlo2-1;start(3)=ixMlo3-1;
    start(ndim+1)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_coarse_block,ierrmpi)
    call MPI_TYPE_COMMIT(type_coarse_block,ierrmpi)
  
    if(stagger_grid) then
      sizes(1)=ixGhi1+1;sizes(2)=ixGhi2+1;sizes(3)=ixGhi3+1;
      sizes(ndim+1)=nws
      subsizes(1)=nx1+1;subsizes(2)=nx2+1;subsizes(3)=nx3+1;
      subsizes(ndim+1)=nws
      start(1)=ixMlo1-1;start(2)=ixMlo2-1;start(3)=ixMlo3-1;
      start(ndim+1)=0
      call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
          MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_block_io_stg,ierrmpi)
      call MPI_TYPE_COMMIT(type_block_io_stg,ierrmpi)
      size_block_io_stg=(nx1+1)*(nx2+1)*(nx3+1)*nws*size_double
  
      sizes(1)=ixGhi1/2+nghostcells+1;sizes(2)=ixGhi2/2+nghostcells+1
      sizes(3)=ixGhi3/2+nghostcells+1;
      sizes(ndim+1)=nws
     do ic3=1,2
     do ic2=1,2
     do ic1=1,2
        do idir=1,ndim
          subsizes(1)=nxCo1+kr(ic1,1)*kr(idir,1)
          subsizes(2)=nxCo2+kr(ic2,1)*kr(idir,2)
          subsizes(3)=nxCo3+kr(ic3,1)*kr(idir,3);
          subsizes(ndim+1)=1
          start(1)=ixMlo1-kr(ic1,1)*kr(idir,1)
          start(2)=ixMlo2-kr(ic2,1)*kr(idir,2)
          start(3)=ixMlo3-kr(ic3,1)*kr(idir,3);
          start(ndim+1)=idir-1
  
          call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
              type_coarse_block_stg(idir,ic1,ic2,ic3),ierrmpi)
          call MPI_TYPE_COMMIT(type_coarse_block_stg(idir,ic1,ic2,ic3),&
             ierrmpi)
        end do
     end do
     end do
     end do
  
      sizes(1)=ixGhi1+1;sizes(2)=ixGhi2+1;sizes(3)=ixGhi3+1;
      sizes(ndim+1)=nws
     do ic3=1,2
     do ic2=1,2
     do ic1=1,2
        do idir=1,3
          subsizes(1)=nxCo1+kr(ic1,1)*kr(idir,1)
          subsizes(2)=nxCo2+kr(ic2,1)*kr(idir,2)
          subsizes(3)=nxCo3+kr(ic3,1)*kr(idir,3);
          subsizes(ndim+1)=1
          start(1)=ixMlo1-kr(ic1,1)*kr(idir,1)+(ic1-1)*nxCo1
          start(2)=ixMlo2-kr(ic2,1)*kr(idir,2)+(ic2-1)*nxCo2
          start(3)=ixMlo3-kr(ic3,1)*kr(idir,3)+(ic3-1)*nxCo3;
          start(ndim+1)=idir-1
          call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_sub_block_stg(idir,&
             ic1,ic2,ic3),ierrmpi)
          call MPI_TYPE_COMMIT(type_sub_block_stg(idir,ic1,ic2,ic3),ierrmpi)
        end do
     end do
     end do
     end do
    end if
  
    sizes(1)=ixGhi1;sizes(2)=ixGhi2;sizes(3)=ixGhi3;
    sizes(ndim+1)=nw
    do ic3=1,2
    do ic2=1,2
    do ic1=1,2
       subsizes(1)=nxCo1;subsizes(2)=nxCo2;subsizes(3)=nxCo3;
       subsizes(ndim+1)=nw
       start(1)=ixMlo1-1+(ic1-1)*nxCo1;start(2)=ixMlo2-1+(ic2-1)*nxCo2
       start(3)=ixMlo3-1+(ic3-1)*nxCo3;
       start(ndim+1)=0
       call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
           MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_sub_block(ic1,ic2,ic3),&
          ierrmpi)
       call MPI_TYPE_COMMIT(type_sub_block(ic1,ic2,ic3),ierrmpi)
    end do
    end do
    end do
  
    sizes(1)=ixGhi1;sizes(2)=ixGhi2;sizes(3)=ixGhi3;
    sizes(ndim+1)=nw
    subsizes(1)=nx1;subsizes(2)=nx2;subsizes(3)=nx3;
    subsizes(ndim+1)=nw
    start(1)=ixMlo1-1;start(2)=ixMlo2-1;start(3)=ixMlo3-1;
    start(ndim+1)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_block_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_block_io,ierrmpi)
    size_block_io=nx1*nx2*nx3*nw*size_double
  
    sizes(1)=ixMhi1-ixMlo1+1;sizes(2)=ixMhi2-ixMlo2+1
    sizes(3)=ixMhi3-ixMlo3+1;
    sizes(ndim+1)=3
    subsizes(1)=sizes(1);subsizes(2)=sizes(2);subsizes(3)=sizes(3);
    subsizes(ndim+1)=3
    start(1)=0;start(2)=0;start(3)=0;
    start(ndim+1)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_block_xcc_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_block_xcc_io,ierrmpi)
  
    sizes(1)=ixMhi1-ixMlo1+2;sizes(2)=ixMhi2-ixMlo2+2
    sizes(3)=ixMhi3-ixMlo3+2;
    sizes(ndim+1)=3
    subsizes(1)=sizes(1);subsizes(2)=sizes(2);subsizes(3)=sizes(3);
    subsizes(ndim+1)=3
    start(1)=0;start(2)=0;start(3)=0;
    start(ndim+1)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_block_xc_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_block_xc_io,ierrmpi)
  
    sizes(1)=ixMhi1-ixMlo1+1;sizes(2)=ixMhi2-ixMlo2+1
    sizes(3)=ixMhi3-ixMlo3+1;
    sizes(ndim+1)=nw+nwauxio
    subsizes(1)=sizes(1);subsizes(2)=sizes(2);subsizes(3)=sizes(3);
    subsizes(ndim+1)=nw+nwauxio
    start(1)=0;start(2)=0;start(3)=0;
    start(ndim+1)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_block_wcc_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_block_wcc_io,ierrmpi)
  
    sizes(1)=ixMhi1-ixMlo1+2;sizes(2)=ixMhi2-ixMlo2+2
    sizes(3)=ixMhi3-ixMlo3+2;
    sizes(ndim+1)=nw+nwauxio
    subsizes(1)=sizes(1);subsizes(2)=sizes(2);subsizes(3)=sizes(3);
    subsizes(ndim+1)=nw+nwauxio
    start(1)=0;start(2)=0;start(3)=0;
    start(ndim+1)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_block_wc_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_block_wc_io,ierrmpi)
  
  end subroutine init_comm_types

  !> Exit MPI-AMRVAC with an error message
  ! cray does not allow char type on the GPU
  subroutine mpistop(message)
#ifndef _CRAYFTN
    !$acc routine
#endif
    use mod_global_parameters
  
    character(len=*), intent(in) :: message !< The error message
    integer                      :: ierrcode
  
    write(*, *) "ERROR for processor", mype, ":"

#ifdef _OPENACC
    write(*, *) message
    STOP
#else
    write(*, *) trim(message)
    call MPI_ABORT(icomm, ierrcode, ierrmpi)
#endif
  
  end subroutine mpistop

  subroutine mpistop_gpu()
    !$acc routine seq
    use mod_global_parameters

    integer                      :: ierrcode

    write(*, *) "ERROR for processor", mype
    STOP
  end subroutine mpistop_gpu
        
  
end module mod_comm_lib
