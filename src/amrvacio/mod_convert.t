module mod_convert

  use mpi
  use mod_variables, only: max_nw 

  implicit none
  public

  abstract interface

     function sub_convert_vars(ixI^L, ixO^L, w, x, nwc) result(wnew)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L,ixO^L, nwc
       double precision, intent(in)    :: w(ixI^S, 1:nw)
       double precision, intent(in)    :: x(ixI^S,1:ndim) 
       double precision    :: wnew(ixO^S, 1:nwc)
     end function sub_convert_vars


   end interface
  type convert_vars_method
    procedure (sub_convert_vars), pointer, nopass :: phys_convert_vars
    integer :: nwc
    character(len=40) :: file_suffix
    character(len=40) :: dataset_names(max_nw)
    type(convert_vars_method), pointer :: next
  end type convert_vars_method
  type(convert_vars_method), pointer :: head_convert_vars_methods 
contains


  subroutine init_convert()
    !initialize the head of the linked list of convert methods
    nullify(head_convert_vars_methods)
  end subroutine init_convert



   function get_names_from_string(aux_variable_names,nwc) result(names) 
    use mod_global_parameters
    character(len=*),intent(in):: aux_variable_names
    integer, intent(in) :: nwc
    character(len=name_len)   :: names(1:nwc)

    integer::  space_position,iw
    character(len=name_len)::  wname
    character(len=std_len)::  scanstring


    ! copied from subroutine getheadernames in calculate_xw
    scanstring=TRIM(adjustl(aux_variable_names))
    space_position=index(scanstring,' ')
    do iw=1,nwc
       do while (space_position==1)
         scanstring=scanstring(2:)
         space_position=index(scanstring,' ')
       enddo
       wname=scanstring(:space_position-1)
       scanstring=scanstring(space_position+1:)
       space_position=index(scanstring,' ')

       names(iw)=TRIM(adjustl(wname))
    enddo
  end function get_names_from_string  


  ! shortcut
  subroutine add_convert_method2(phys_convert_vars, nwc, aux_variable_names, file_suffix)
    use mod_global_parameters
    integer, intent(in) :: nwc
    character(len=*),intent(in):: aux_variable_names
    character(len=*), intent(in) :: file_suffix

    interface
     function phys_convert_vars(ixI^L, ixO^L, w, x, nwc) result(wnew)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, nwc
       double precision, intent(in)    :: w(ixI^S, 1:nw)
       double precision, intent(in)    :: x(ixI^S,1:ndim) 
       double precision    :: wnew(ixO^S, 1:nwc)
     end function phys_convert_vars
    end interface

    call add_convert_method(phys_convert_vars, nwc, get_names_from_string(aux_variable_names,nwc), file_suffix)

  end subroutine add_convert_method2

  subroutine add_convert_method(phys_convert_vars, nwc, dataset_names, file_suffix)
    integer, intent(in) :: nwc
    character(len=*), intent(in) :: dataset_names(:)
    character(len=*), intent(in) :: file_suffix

    interface
     function phys_convert_vars(ixI^L, ixO^L, w, x, nwc) result(wnew)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, nwc
       double precision, intent(in)    :: w(ixI^S, 1:nw)
       double precision, intent(in)    :: x(ixI^S,1:ndim) 
       double precision    :: wnew(ixO^S, 1:nwc)
     end function phys_convert_vars
    end interface

    type(convert_vars_method), pointer  :: temp

    if(nwc .gt. max_nw) then
      call mpistop("INCREASE max_nw ")
    endif  

    allocate(temp)
    temp%phys_convert_vars => phys_convert_vars
    temp%nwc = nwc
    temp%file_suffix = file_suffix
    temp%dataset_names(1:nwc) = dataset_names(1:nwc)
    temp%next => head_convert_vars_methods
    head_convert_vars_methods => temp
  end subroutine add_convert_method


  subroutine convert_all()
    type(convert_vars_method), pointer  :: temp
    temp => head_convert_vars_methods
    do while(associated(temp))
      call convert_dat_generic(temp%nwc, temp%dataset_names, temp%file_suffix, temp%phys_convert_vars)
      temp=>temp%next
    end do
  end subroutine convert_all

  ! Copied from subroutine write_snapshot in amrvacio/mod_input_output
  ! the staggered values are not saved in this subroutine!
  subroutine convert_dat_generic(nwc, dataset_names, file_suffix, convert_vars)

    use mod_forest
    use mod_global_parameters
    use mod_input_output, only: count_ix, create_output_file, snapshot_write_header1, block_shape_io

    integer                       :: file_handle, igrid, Morton_no, iwrite
    integer                       :: ipe, ix_buffer(2*ndim+1), n_values
    integer                       :: ixI^L, ixO^L, n_ghost(2*ndim)
    integer                       :: ixOs^L,n_values_stagger
    integer                       :: iorecvstatus(MPI_STATUS_SIZE)
    integer                       :: ioastatus(MPI_STATUS_SIZE)
    integer                       :: igrecvstatus(MPI_STATUS_SIZE)
    integer                       :: istatus(MPI_STATUS_SIZE)
    type(tree_node), pointer      :: pnode
    integer(kind=MPI_OFFSET_KIND) :: offset_tree_info
    integer(kind=MPI_OFFSET_KIND) :: offset_block_data
    integer(kind=MPI_OFFSET_KIND) :: offset_offsets
    double precision, allocatable :: w_buffer(:)
    double precision, allocatable:: converted_vars(:^D&,:)

    integer :: idim, itag
    integer, allocatable                       :: block_ig(:, :)
    integer, allocatable                       :: block_lvl(:)
    integer(kind=MPI_OFFSET_KIND), allocatable :: block_offset(:)

    integer, intent(in) :: nwc
    character(len=*), intent(in) :: dataset_names(:) 
    character(len=*), intent(in) :: file_suffix 
    interface

      function convert_vars(ixI^L, ixO^L,w,x,nwc) result(wres)
        use mod_global_parameters
        integer, intent(in) :: ixI^L, ixO^L, nwc
        double precision, intent(in) ::  x(ixI^S,1:ndim)
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision  :: wres(ixO^S,1:nwc)
      end function convert_vars

    end interface


    call MPI_BARRIER(icomm, ierrmpi)

    ! Allocate send/receive buffer
    n_values = count_ix(ixG^LL) * nw
    allocate(w_buffer(n_values))

    ! Allocate arrays with information about grid blocks
    allocate(block_ig(ndim, nleafs))
    allocate(block_lvl(nleafs))
    allocate(block_offset(nleafs+1))

    ! master processor
    if (mype==0) then
      call create_output_file(file_handle, snapshotnext, ".dat", file_suffix)

      ! Don't know offsets yet, we will write header again later
      offset_tree_info = -1
      offset_block_data = -1
      call snapshot_write_header1(file_handle, offset_tree_info, &
           offset_block_data, dataset_names, nwc)

      call MPI_File_get_position(file_handle, offset_tree_info, ierrmpi)

      call write_forest(file_handle)

      ! Collect information about the spatial index (ig^D) and refinement level
      ! of leaves
      do Morton_no = Morton_start(0), Morton_stop(npe-1)
         igrid = sfc(1, Morton_no)
         ipe = sfc(2, Morton_no)
         pnode => igrid_to_node(igrid, ipe)%node

         block_ig(:, Morton_no) = [ pnode%ig^D ]
         block_lvl(Morton_no) = pnode%level
         block_offset(Morton_no) = 0 ! Will be determined later
      end do

      call MPI_FILE_WRITE(file_handle, block_lvl, size(block_lvl), &
           MPI_INTEGER, istatus, ierrmpi)

      call MPI_FILE_WRITE(file_handle, block_ig, size(block_ig), &
           MPI_INTEGER, istatus, ierrmpi)

      ! Block offsets are currently unknown, but will be overwritten later
      call MPI_File_get_position(file_handle, offset_offsets, ierrmpi)
      call MPI_FILE_WRITE(file_handle, block_offset(1:nleafs), nleafs, &
           MPI_OFFSET, istatus, ierrmpi)

      call MPI_File_get_position(file_handle, offset_block_data, ierrmpi)

      ! Check whether data was written as expected
      if (offset_block_data - offset_tree_info /= &
           (nleafs + nparents) * size_logical + &
           nleafs * ((1+ndim) * size_int + 2 * size_int)) then
        if (mype == 0) then
          print *, "Warning: MPI_OFFSET type /= 8 bytes"
          print *, "This *could* cause problems when reading .dat files"
        end if
      end if

      block_offset(1) = offset_block_data
      iwrite = 0
    end if

    do Morton_no=Morton_start(mype), Morton_stop(mype)
      igrid  = sfc_to_igrid(Morton_no)
      itag   = Morton_no
      block=>ps(igrid)
      ! this might be used in convert function, 
      ! it was not used when the output is already computed vars  (write_snapshot)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

      ! start copied from block_shape_io, 
      ! because nwc is needed  as parameter
      ! TODO check if this will be used elsewehere and put it in separate subroutine
      n_ghost(:) = 0
  
      if(save_physical_boundary) then
        do idim=1,ndim
          ! Include ghost cells on lower boundary
          if(ps(igrid)%is_physical_boundary(2*idim-1)) n_ghost(idim)=nghostcells
          ! Include ghost cells on upper boundary
          if(ps(igrid)%is_physical_boundary(2*idim)) n_ghost(ndim+idim)=nghostcells
        end do
      end if
  
      {ixOmin^D = ixMlo^D - n_ghost(^D)\}
      {ixOmax^D = ixMhi^D + n_ghost(ndim+^D)\}
 
      n_values = count_ix(ixO^L) * nwc
      ! end copied from block_shape_io

      {ixImin^D = ixGlo^D\}
      {ixImax^D = ixGhi^D\}

      w_buffer(1:n_values) = pack(convert_vars(ixI^L, ixO^L,ps(igrid)%w(ixI^S, 1:nw), ps(igrid)%x(ixI^S, 1:ndim),nwc), .true.)

      ix_buffer(1) = n_values
      ix_buffer(2:) = n_ghost

      if (mype /= 0) then
        call MPI_SEND(ix_buffer, 2*ndim+1, &
             MPI_INTEGER, 0, itag, icomm, ierrmpi)
        call MPI_SEND(w_buffer, n_values, &
             MPI_DOUBLE_PRECISION, 0, itag, icomm, ierrmpi)
      else
        iwrite = iwrite+1
        call MPI_FILE_WRITE(file_handle, ix_buffer(2:), &
             2*ndim, MPI_INTEGER, istatus, ierrmpi)
        call MPI_FILE_WRITE(file_handle, w_buffer, &
             n_values, MPI_DOUBLE_PRECISION, istatus, ierrmpi)

        ! Set offset of next block
        block_offset(iwrite+1) = block_offset(iwrite) + &
             int(n_values, MPI_OFFSET_KIND) * size_double + &
             2 * ndim * size_int
      end if
    end do

    ! Write data communicated from other processors
    if (mype == 0) then
      do ipe = 1, npe-1
        do Morton_no=Morton_start(ipe), Morton_stop(ipe)
          iwrite=iwrite+1
          itag=Morton_no

          call MPI_RECV(ix_buffer, 2*ndim+1, MPI_INTEGER, ipe, itag, icomm,&
               igrecvstatus, ierrmpi)
          n_values = ix_buffer(1)

          call MPI_RECV(w_buffer, n_values, MPI_DOUBLE_PRECISION,&
               ipe, itag, icomm, iorecvstatus, ierrmpi)

          call MPI_FILE_WRITE(file_handle, ix_buffer(2:), &
             2*ndim, MPI_INTEGER, istatus, ierrmpi)
          call MPI_FILE_WRITE(file_handle, w_buffer, &
             n_values, MPI_DOUBLE_PRECISION, istatus, ierrmpi)

          ! Set offset of next block
          block_offset(iwrite+1) = block_offset(iwrite) + &
               int(n_values, MPI_OFFSET_KIND) * size_double + &
               2 * ndim * size_int
        end do
      end do

      ! Write block offsets (now we know them)
      call MPI_FILE_SEEK(file_handle, offset_offsets, MPI_SEEK_SET, ierrmpi)
      call MPI_FILE_WRITE(file_handle, block_offset(1:nleafs), nleafs, &
           MPI_OFFSET, istatus, ierrmpi)

      ! Write header again, now with correct offsets
      call MPI_FILE_SEEK(file_handle, 0_MPI_OFFSET_KIND, MPI_SEEK_SET, ierrmpi)
      call snapshot_write_header1(file_handle, offset_tree_info, &
           offset_block_data, dataset_names, nwc)

      call MPI_FILE_CLOSE(file_handle, ierrmpi)
    end if

    call MPI_BARRIER(icomm, ierrmpi)

  end subroutine convert_dat_generic

end module mod_convert
