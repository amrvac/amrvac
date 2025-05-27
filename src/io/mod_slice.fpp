!> Writes D-1 slice, can do so in various formats, depending on slice_type
module mod_slice
  use mod_basic_types
  use mod_comm_lib, only: mpistop
  implicit none

  !> Maximum number of slices
  integer, parameter :: nslicemax=1000

  !> Slice coordinates, see @ref slices.md
  double precision :: slicecoord(nslicemax)

  !> the file number of slices
  integer :: slicenext

  !> Number of slices to output
  integer :: nslices

  !> The slice direction for each slice
  integer :: slicedir(nslicemax)

  !> choose data type of slice: vtu, vtuCC, dat, or csv
  character(len=std_len) :: slice_type

  !> tag for MPI message
  integer, private :: itag

contains

  subroutine write_slice
    use mod_global_parameters
    ! Writes a D-1 slice 
    ! by Oliver Porth
    ! 22.Nov 2011
    integer :: islice

    do islice=1,nslices
       call put_slice(slicedir(islice),slicecoord(islice))
    end do

    slicenext=slicenext+1
  end subroutine write_slice

  subroutine put_slice(dir,xslice)
    use mod_forest, only: Morton_sub_start, Morton_sub_stop
    use mod_global_parameters
    ! Writes a D-1 slice 
    ! For ONED simulations, output will be appended to one csv-file per slice
    ! slices are sensitive to the saveprim switch and 
    ! can contain auxiliary io variables (nwauxio)
    ! Thus csv or vtu(CC)-files with primitive variables are obtained.  
    ! when slice_type='dat', we save a D-1 dat file for potential restarts
    ! by Oliver Porth
    ! 22.Nov 2011
    integer, intent(in) :: dir
    double precision, intent(in) :: xslice
    ! .. local ..
    integer :: Njgrid, jgrid
    integer, dimension(ndim-1) :: ixsubGlo, ixsubGhi
    integer, dimension(ndim-1) :: ixsubMlo, ixsubMhi
    integer :: size_subblock_io, nx1,nx2,nx3, slice_fh, nwexpand
    integer :: type_subblock_io, type_subblockC_io, type_subblock_x_io,&
        type_subblockC_x_io
    integer, dimension(ndim) :: sizes, subsizes, start
    double precision,dimension(0:nw+nwauxio)          :: normconv 
  
    ! Preamble: 
    nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
    slice_fh=unitslice

    if (ndim==1) then
       nwexpand = nwauxio
    else
       if(slice_type/='dat')then
          nwexpand = nwauxio
       else
          nwexpand = 0
       endif
    end if

    ! Do a last consistency check:
    select case(dir)
       case(1)
       if(xslice<xprobmin1.or.xslice>xprobmax1) call &
          mpistop("slice out of bounds")
       
       case(2)
       if(xslice<xprobmin2.or.xslice>xprobmax2) call &
          mpistop("slice out of bounds")
       
       case(3)
       if(xslice<xprobmin3.or.xslice>xprobmax3) call &
          mpistop("slice out of bounds")
       
    end select

    ! Traverse the forest and fill nodes:
    call select_slice(dir,xslice,.false.,slice_fh,normconv)

    ! Create the MPI-datatype and select indices:
    
    select case(dir)
    case (1)
       sizes(1) = ixGhi2; sizes(2) = ixGhi3;
       ixsubGlo(1) = ixGlo2; ixsubGlo(2) = ixGlo3;
       ixsubGhi(1) = ixGhi2; ixsubGhi(2) = ixGhi3;
       subsizes(1)=nx2;subsizes(2)=nx3;
       start(1)=ixMlo2-1;start(2)=ixMlo3-1;
       size_subblock_io=nx2*nx3*(nw+nwexpand)*size_double
    case (2)
       sizes(1) = ixGhi1; sizes(2) = ixGhi3;
       ixsubGlo(1) = ixGlo1; ixsubGlo(2) = ixGlo3;
       ixsubGhi(1) = ixGhi1; ixsubGhi(2) = ixGhi3;
       subsizes(1)=nx1;subsizes(2)=nx3;
       start(1)=ixMlo1-1;start(2)=ixMlo3-1;
       size_subblock_io=nx1*nx3*(nw+nwexpand)*size_double
    case (3)
       ixsubGlo(1) = ixGlo1; ixsubGlo(2) = ixGlo2;
       ixsubGhi(1) = ixGhi1; ixsubGhi(2) = ixGhi2;
       sizes(1) = ixGhi1; sizes(2) = ixGhi2;
       subsizes(1)=nx1;subsizes(2)=nx2;
       start(1)=ixMlo1-1;start(2)=ixMlo2-1;
       size_subblock_io=nx1*nx2*(nw+nwexpand)*size_double
    case default
       call mpistop("slice direction not clear in put_slice")
    end select
   
    
    

    
    ixsubMlo(2-1) = ixsubGlo(2-1)+nghostcells
    ixsubMlo(3-1) = ixsubGlo(3-1)+nghostcells;
    ixsubMhi(2-1) = ixsubGhi(2-1)-nghostcells
    ixsubMhi(3-1) = ixsubGhi(3-1)-nghostcells;
   

    sizes(ndim)=nw+nwexpand
    subsizes(ndim)=nw+nwexpand
    start(ndim)=0

    ! Types for center variables:
    call MPI_TYPE_CREATE_SUBARRAY(ndim,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
       MPI_DOUBLE_PRECISION, type_subblock_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_subblock_io,ierrmpi)

    sizes(ndim)=3
    subsizes(ndim)=3
    start(ndim)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
       MPI_DOUBLE_PRECISION, type_subblock_x_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_subblock_x_io,ierrmpi)


    ! Types for corner variables:
    subsizes(1:ndim-1) = subsizes(1:ndim-1) + 1
    start(1:ndim-1)    = start(1:ndim-1) - 1 
    sizes(ndim)=nw+nwexpand
    subsizes(ndim)=nw+nwexpand
    start(ndim)=0
    
    call MPI_TYPE_CREATE_SUBARRAY(ndim,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
       MPI_DOUBLE_PRECISION, type_subblockC_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_subblockC_io,ierrmpi)

    sizes(ndim)=3
    subsizes(ndim)=3
    start(ndim)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
       MPI_DOUBLE_PRECISION, type_subblockC_x_io,ierrmpi)
    call MPI_TYPE_COMMIT(type_subblockC_x_io,ierrmpi)

    
    ! local number of sub-grids:
    Njgrid=Morton_sub_stop(mype)-Morton_sub_start(mype)+1

    ! Now output using various schemes: 
    if (ndim==1) then 
       call put_slice_zerod
    else
       select case(slice_type)
       case ('csv')
          call put_slice_csv
       case ('dat')
          call put_slice_dat
       case ('vtu', 'vtuCC')
          call put_slice_vtu
       end select
    end if

    do jgrid=1,Njgrid
       call dealloc_subnode(jgrid)
    end do

    call MPI_TYPE_FREE(type_subblock_io,ierrmpi)
    call MPI_TYPE_FREE(type_subblock_x_io,ierrmpi)
    call MPI_TYPE_FREE(type_subblockC_io,ierrmpi)
    call MPI_TYPE_FREE(type_subblockC_x_io,ierrmpi)


  contains

    subroutine put_slice_vtu

      use mod_calculate_xw
      character(len=1024) :: filename, xlabel
      character(len=79)   :: xxlabel
      logical             :: fileopen
      character(len=name_len) :: wnamei(1:nw+nwauxio),&
         xandwnamei(1:ndim+nw+nwauxio)
      character(len=1024) :: outfilehead
      integer :: status(MPI_STATUS_SIZE), ipe

      if (mype==0) then

         inquire(slice_fh,opened=fileopen)
         if(.not.fileopen)then
            ! generate filename: 
            write(xlabel,"(D9.2)")xslice
            xxlabel=trim(xlabel)
            if(xslice>=zero)then
               write(xxlabel(1:1),"(a)") "+"
            endif
            write(filename,"(a,i1.1,a,i4.4,a)") TRIM(base_filename)//'_d',dir,&
               '_x'//trim(xxlabel)//'_n',slicenext,'.vtu'
            open(slice_fh,file=filename,status='unknown',form='formatted')
         end if
         ! get and write the header: 
         call getheadernames(wnamei,xandwnamei,outfilehead)
         ! generate xml header
         write(slice_fh,'(a)')'<?xml version="1.0"?>'
         write(slice_fh,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
         if(type_endian==1)then
            write(slice_fh,'(a)')' version="0.1" byte_order="LittleEndian">'
         else
            write(slice_fh,'(a)')' version="0.1" byte_order="BigEndian">'
         endif
         write(slice_fh,'(a)')'  <UnstructuredGrid>'
         write(slice_fh,'(a)')'<FieldData>'
         write(slice_fh,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
            'NumberOfTuples="1" format="ascii">'
         write(slice_fh,*) real(global_time*time_convert_factor)
         write(slice_fh,'(a)')'</DataArray>'
         write(slice_fh,'(a)')'</FieldData>'

         ! write to file:
         do jgrid=1, Njgrid
            call write_slice_vtk(jgrid,slice_fh,wnamei)
         end do

         ! create a recv buffer using allocate, will be deallocated at the end of the routine:
         call alloc_subnode(Njgrid+1,dir,nwauxio)

      end if

      ! Also communicate the normconv array since processor zero might not have it yet:
      if (npe>1) then
         do ipe=1,npe-1
            do jgrid=1,Morton_sub_stop(ipe)-Morton_sub_start(ipe)+1
               itag=Morton_sub_start(ipe)+jgrid-1
               itag=itag*5
               if (ipe == mype ) then 
                  call MPI_SEND(ps_sub(jgrid)%x,1,type_subblock_x_io,0,itag,&
                     icomm,ierrmpi)
                  call MPI_SEND(ps_sub(jgrid)%w,1,type_subblock_io,0,itag+1,&
                     icomm,ierrmpi)
                  call MPI_SEND(ps_sub(jgrid)%xC,1,type_subblockC_x_io,0,&
                     itag+2,icomm,ierrmpi)
                  call MPI_SEND(ps_sub(jgrid)%wC,1,type_subblockC_io,0,itag+3,&
                     icomm,ierrmpi)
                  call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION,0,&
                     itag+4,icomm,ierrmpi)
               end if
               if (mype == 0) then
                  call MPI_RECV(ps_sub(Njgrid+1)%x,1,type_subblock_x_io,ipe,&
                     itag,icomm,status,ierrmpi)
                  call MPI_RECV(ps_sub(Njgrid+1)%w,1,type_subblock_io,ipe,&
                     itag+1,icomm,status,ierrmpi)
                  call MPI_RECV(ps_sub(Njgrid+1)%xC,1,type_subblockC_x_io,ipe,&
                     itag+2,icomm,status,ierrmpi)
                  call MPI_RECV(ps_sub(Njgrid+1)%wC,1,type_subblockC_io,ipe,&
                     itag+3,icomm,status,ierrmpi)
                  call MPI_RECV(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION,ipe,&
                     itag+4,icomm,status,ierrmpi)
                  call write_slice_vtk(Njgrid+1,slice_fh,wnamei)
               end if
            end do
         end do
      endif

      if (mype==0) then

         write(slice_fh,'(a)')'</UnstructuredGrid>'
         write(slice_fh,'(a)')'</VTKFile>'
         close(slice_fh)
         call dealloc_subnode(Njgrid+1)

      end if

    end subroutine put_slice_vtu

    subroutine write_slice_vtk(jgrid,slice_fh,wnamei)

      ! this only works for 2D and 3D, 1D reduction (line to point) not allowed
      integer, intent(in)           :: jgrid, slice_fh
      character(len=name_len), intent(in) :: wnamei(1:nw+nwauxio)
      ! This remainder part only for more than 1D, but nesting with NOONED gives problems 
      
      ! .. local ..
      integer                       :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,&
         ixCmax3, ixCCmin1,ixCCmin2,ixCCmin3,ixCCmax1,ixCCmax2,ixCCmax3, nc,&
          np, iw
      integer                       :: nx1,nx2, nxC1,nxC2, icell, ix1,ix2
      double precision              :: x_VTK(1:3)
      integer                       :: VTK_type
      double precision, parameter   :: minvalue = 1.0d-99, maxvalue = 1.0d+99

      ixCCmin1 = ixsubMlo(1);ixCCmin2 = ixsubMlo(2);
      ixCCmax1 = ixsubMhi(1);ixCCmax2 = ixsubMhi(2);
      ixCmin1  = ixsubMlo(1)-1;ixCmin2  = ixsubMlo(2)-1;
      ixCmax1  = ixsubMhi(1);ixCmax2  = ixsubMhi(2);

      nx1=ixCCmax1-ixCCmin1+1;nx2=ixCCmax2-ixCCmin2+1;
      nxC1=nx1+1;nxC2=nx2+1;
      nc=nx1*nx2      ! Number of cells per subgrid
      np=nxC1*nxC2     ! Number of corner points per subgrid

      ! we write out every grid as one VTK PIECE
      write(slice_fh,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
         '" NumberOfCells="',nc,'">'

      !==============================
      ! celldata or pointdata?
      !==============================
      select case(slice_type)

      case('vtuCC') ! celldata
         write(slice_fh,'(a)')'<CellData>'
         do iw=1,nw+nwauxio
            if(iw<=nw) then 
               if(.not.w_write(iw)) cycle
            endif
            write(slice_fh,'(a,a,a)')'<DataArray type="Float64" Name="',&
               TRIM(wnamei(iw)),'" format="ascii">'
            write(slice_fh,'(200(1pe14.6))') &
               ((roundoff_minmax(ps_sub(jgrid)%w(ix1,ix2,iw)*normconv(iw),&
               minvalue,maxvalue),ix1=ixCCmin1,ixCCmax1),ix2=ixCCmin2,&
               ixCCmax2)
            write(slice_fh,'(a)')'</DataArray>'
         enddo
         write(slice_fh,'(a)')'</CellData>'


      case('vtu') ! pointdata
         write(slice_fh,'(a)')'<PointData>'
         do iw=1,nw+nwauxio
            if(iw<=nw) then 
               if(.not.w_write(iw)) cycle
            endif
            write(slice_fh,'(a,a,a)')'<DataArray type="Float64" Name="',&
               TRIM(wnamei(iw)),'" format="ascii">'
            write(slice_fh,'(200(1pe14.6))') &
               ((roundoff_minmax(ps_sub(jgrid)%wC(ix1,ix2,iw)*normconv(iw),&
               minvalue,maxvalue),ix1=ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2)
            write(slice_fh,'(a)')'</DataArray>'
         enddo
         write(slice_fh,'(a)')'</PointData>'
         
      
      end select
      !==============================
      ! Done: celldata or pointdata?
      !==============================

      !==============================
      ! output Cornerpoints
      !==============================
      write(slice_fh,'(a)')'<Points>'
      write(slice_fh,'(a)'&
         )'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
 !write cell corner coordinates in a backward dimensional loop, always 3D output
      do ix2=ixCmin2,ixCmax2 
      do ix1=ixCmin1,ixCmax1 
            x_VTK(1:3)=zero;
            x_VTK(1:ndim)=ps_sub(jgrid)%xC(ix1,ix2,1:ndim)*normconv(0);
            write(slice_fh,'(3(1pe14.6))') x_VTK
      end do 
      end do 
      write(slice_fh,'(a)')'</DataArray>'
      write(slice_fh,'(a)')'</Points>'
      !==============================
      ! Done: output Cornerpoints
      !==============================

      !==============================
      ! cell Metainformation
      !==============================
      write(slice_fh,'(a)')'<Cells>'

      ! connectivity part
      write(slice_fh,'(a)'&
         )'<DataArray type="Int32" Name="connectivity" format="ascii">'

       do ix2=1,nx2
       do ix1=1,nx1
      
      write(slice_fh,'(4(i7,1x))')(ix2-1)*nxC1+ix1-1, (ix2-1)*nxC1+ix1,&
         ix2*nxC1+ix1-1,ix2*nxC1+ix1
      end do
       end do

      write(slice_fh,'(a)')'</DataArray>'

      ! offsets data array
      write(slice_fh,'(a)'&
         )'<DataArray type="Int32" Name="offsets" format="ascii">'
      do icell=1,nc
         write(slice_fh,'(i7)') icell*(2**(3-1))
      end do
      write(slice_fh,'(a)')'</DataArray>'

      ! VTK cell type data array
      write(slice_fh,'(a)'&
         )'<DataArray type="Int32" Name="types" format="ascii">'
      ! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax
      
       VTK_type=8 
      do icell=1,nc
         write(slice_fh,'(i2)') VTK_type
      enddo
      write(slice_fh,'(a)')'</DataArray>'
      
      write(slice_fh,'(a)')'</Cells>'
      !==============================
      ! Done: cell Metainformation
      !==============================
      write(slice_fh,'(a)')'</Piece>'

      
    end subroutine write_slice_vtk

    subroutine put_slice_csv

      use mod_calculate_xw
      character(len=1024)           :: filename, xlabel
      character(len=79)             :: xxlabel
      logical                       :: fileopen
      character(len=name_len)       :: wnamei(1:nw+nwauxio),&
         xandwnamei(1:ndim+nw+nwauxio)
      character(len=1024)           :: outfilehead
      integer                       :: iw, ipe, itag
      character(len=1024)           :: line
      integer                       :: status(MPI_STATUS_SIZE)

      if (mype==0) then
         inquire(slice_fh,opened=fileopen)
         if(.not.fileopen)then
            ! generate filename: 
            write(xlabel,"(D9.2)")xslice
            xxlabel=trim(xlabel)
            if(xslice>=zero)then
               write(xxlabel(1:1),"(a)") "+"
            endif
            write(filename,"(a,i1.1,a,i4.4,a)") TRIM(base_filename)//'_d',dir,&
               '_x'//trim(xxlabel)//'_n',slicenext,'.csv'
            open(slice_fh,file=filename,status='unknown',form='formatted')
         end if
         ! get and write the header: 
         call getheadernames(wnamei,xandwnamei,outfilehead)
         line=''
         do iw=1,ndim+nw+nwauxio-1
            line = trim(line)//trim(xandwnamei(iw))//', '
         end do
         line = trim(line)//trim(xandwnamei(ndim+nw+nwauxio))
         write(slice_fh,'(a)')trim(line)
         ! create a recv buffer using allocate, will be deallocated at the end of the routine:
         call alloc_subnode(Njgrid+1,dir,nwauxio)

         ! write to file:
         do jgrid=1, Njgrid
            call put_slice_line(jgrid,slice_fh)
         end do
      end if

      ! Also communicate the normconv array since processor zero might not have it yet:
      if (npe>1) then
         do ipe=1,npe-1
            do jgrid=1,Morton_sub_stop(ipe)-Morton_sub_start(ipe)+1
               itag=Morton_sub_start(ipe)+jgrid-1
               if (ipe == mype ) then 
                  call MPI_SEND(ps_sub(jgrid)%x,1,type_subblock_x_io,0,itag,&
                     icomm,ierrmpi)
                  call MPI_SEND(ps_sub(jgrid)%w,1,type_subblock_io,0,itag,&
                     icomm,ierrmpi)
                  call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION,0,&
                     itag,icomm,ierrmpi)
               end if
               if (mype == 0) then
                  call MPI_RECV(ps_sub(Njgrid+1)%x,1,type_subblock_x_io,ipe,&
                     itag,icomm,status,ierrmpi)
                  call MPI_RECV(ps_sub(Njgrid+1)%w,1,type_subblock_io,ipe,itag,&
                     icomm,status,ierrmpi)
                  call MPI_RECV(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION,ipe,&
                     itag,icomm,status,ierrmpi)
                  call put_slice_line(Njgrid+1,slice_fh)
               end if
            end do
         end do
      endif

      if (mype==0) then
         close(slice_fh)
         call dealloc_subnode(Njgrid+1)
      end if

    end subroutine put_slice_csv

    subroutine put_slice_line(jout,file_handle)
      integer, intent(in) :: jout, file_handle
      ! .. local ..
      character(len=1024) ::line, data
      integer :: ix1,ix2,ix3,idir,iw
      double precision, parameter :: minvalue = 1.0d-99, maxvalue = 1.0d+99

      
      do ix2=ixsubMlo(2),ixsubMhi(2)
         do ix1=ixsubMlo(1),ixsubMhi(1)
            
            
               ! Format the line:
               line = ''
               do idir=1,ndim
                  
                  write(data,"(es14.6)")roundoff_minmax(ps_sub(jout)%x(ix1,ix2,&
                     idir),minvalue,maxvalue)
                 
                  
                  

                  line = trim(line)//trim(data)//', '
               end do
               do iw = 1,nw+nwauxio-1
                  
                  write(data,"(es14.6)")roundoff_minmax(ps_sub(jout)%w(ix1,ix2,&
                     iw)*normconv(iw),minvalue,maxvalue)
                 
                  
                  
                  line = trim(line)//trim(data)//', '
               end do
               
               write(data,"(es14.6)")roundoff_minmax(ps_sub(jout)%w(ix1,ix2,&
                  nw+nwauxio)*normconv(nw+nwauxio),minvalue,maxvalue)
              
               
               line = trim(line)//trim(data)
               write(file_handle,'(a)')trim(line)
               
            
         end do
      end do
      

    end subroutine put_slice_line

    subroutine put_slice_dat

      integer, dimension(max_blocks) :: iorequest
      integer, dimension(MPI_STATUS_SIZE,max_blocks) :: iostatus
      integer(kind=MPI_OFFSET_KIND) :: offset
      integer :: nsubleafs
      character(len=1024) :: filename, xlabel
      character(len=79)   :: xxlabel
      integer :: amode, status(MPI_STATUS_SIZE), iwrite

      nsubleafs=Morton_sub_stop(npe-1)
      ! generate filename
      write(xlabel,"(D9.2)")xslice
      xxlabel=trim(xlabel)
      if(xslice>=zero)then
         write(xxlabel(1:1),"(a)") "+"
      endif
      write(filename,"(a,i1.1,a,i4.4,a)") TRIM(base_filename)//'_d',dir,&
         '_x'//trim(xxlabel)//'_n',slicenext,'.dat'

      if(mype==0) then
         open(unit=slice_fh,file=filename,status='replace')
         close(unit=slice_fh)
      end if

      amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call MPI_FILE_OPEN(icomm,filename,amode,MPI_INFO_NULL,slice_fh,ierrmpi)
      iorequest=MPI_REQUEST_NULL
      iwrite=0

      do jgrid=1,Njgrid
         iwrite=iwrite+1
         offset=int(size_subblock_io,kind=MPI_OFFSET_KIND) &
            *int(Morton_sub_start(mype)+jgrid-2,kind=MPI_OFFSET_KIND)
         call MPI_FILE_IWRITE_AT(slice_fh,offset,ps_sub(jgrid)%w,1,&
            type_subblock_io, iorequest(iwrite),ierrmpi)
      end do

      if (iwrite>0) call MPI_WAITALL(iwrite,iorequest,iostatus,ierrmpi)
      call MPI_BARRIER(icomm, ierrmpi)
      call MPI_FILE_CLOSE(slice_fh,ierrmpi)

      if (mype==0) then
         amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
         call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,&
             slice_fh,ierrmpi)

         call select_slice(dir,xslice,.true.,slice_fh,normconv)

         call MPI_FILE_WRITE(slice_fh,subsizes(2-1),1,MPI_INTEGER,status,&
            ierrmpi)
         call MPI_FILE_WRITE(slice_fh,subsizes(3-1),1,MPI_INTEGER,status,&
            ierrmpi)
!         call MPI_FILE_WRITE(slice_fh,eqpar,neqpar+nspecialpar, &
!              MPI_DOUBLE_PRECISION,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,nsubleafs,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,levmax_sub,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,ndim-1,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,ndir,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,nw,1,MPI_INTEGER,status,ierrmpi)
!         call MPI_FILE_WRITE(slice_fh,neqpar+nspecialpar,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,it,1,MPI_INTEGER,status,ierrmpi)
         call MPI_FILE_WRITE(slice_fh,global_time,1,MPI_DOUBLE_PRECISION,&
            status,ierrmpi)

         call MPI_FILE_CLOSE(slice_fh,ierrmpi)
      end if

    end subroutine put_slice_dat

    subroutine put_slice_zerod

      use mod_calculate_xw
      integer::  iw
      character(len=name_len) :: wnamei(1:nw+nwauxio),&
         xandwnamei(1:ndim+nw+nwauxio)
      character(len=1024) :: outfilehead
      logical, save :: opened=.false.
      character(len=1024) ::line, data
      character(len=1024) :: filename, xlabel
      character(len=79)   :: xxlabel
      integer :: amode, iwrite, status(MPI_STATUS_SIZE)

      
    end subroutine put_slice_zerod

  end subroutine put_slice

  subroutine select_slice(dir,xslice,writeonly,file_handle,normconv)
    use mod_forest, only: tree_node_ptr, tree_root, Morton_sub_start,&
        Morton_sub_stop
    use mod_global_parameters
    integer, intent(in) :: dir
    double precision, intent(in) :: xslice
    integer, intent(in) :: file_handle
    logical, intent(in) :: writeonly
    double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv 
    ! .. local ..
    integer :: ig1,ig2,ig3, jgrid, slice_fh, ipe, mylevmax
    integer, dimension(nlevelshi) :: igslice

    jgrid = 0
    mylevmax = 0

    ! Find the global slice index for every level:
    call get_igslice(dir,xslice,igslice)

    ! Traverse forest to find grids indicating slice:
    
    select case(dir)
    case (1)
       ig1 = igslice(1)
       do ig3=1,ng3(1)
          do ig2=1,ng2(1)
             call traverse_slice(tree_root(ig1,ig2,ig3))
          end do
       end do
    case (2)
       ig2 = igslice(1)
       do ig3=1,ng3(1)
          do ig1=1,ng1(1)
             call traverse_slice(tree_root(ig1,ig2,ig3))
          end do
       end do
    case (3)
       ig3 = igslice(1)
       do ig2=1,ng2(1)
          do ig1=1,ng1(1)
             call traverse_slice(tree_root(ig1,ig2,ig3))
          end do
       end do
    case default
       call mpistop("slice direction not clear in select_slice")
    end select
   
    
    

    if (.not.writeonly) then
       ! Synchronize the levmax_sub for output (only rank 0 needs it): 
       levmax_sub = mylevmax
       call MPI_ALLREDUCE(MPI_IN_PLACE,levmax_sub,1,MPI_INTEGER,MPI_MAX,icomm,&
          ierrmpi)

       ! Communicate the subgrid indices according to new Morton sub-sfc:
       Morton_sub_start(:) = 1
       do ipe=0,npe-1
          call MPI_GATHER(jgrid,1,MPI_INTEGER,Morton_sub_stop,1,MPI_INTEGER,&
             ipe,icomm,ierrmpi)
       end do

       do ipe = 0, npe-2
          Morton_sub_start(ipe+1) = Morton_sub_stop(ipe)+Morton_sub_start(ipe+&
             1)
          Morton_sub_stop(ipe+1)  = Morton_sub_stop(ipe)+Morton_sub_stop(ipe+&
             1)
       end do
    end if

  contains

    recursive subroutine traverse_slice(tree)
      implicit none
      type(tree_node_ptr) :: tree
      integer :: ic1,ic2,ic3
      integer, dimension(MPI_STATUS_SIZE) :: status

      if (writeonly) then
         call MPI_FILE_WRITE(file_handle,tree%node%leaf,1,MPI_LOGICAL, status,&
            ierrmpi)
      end if

      if (tree%node%leaf) then
         if (tree%node%ipe == mype.and..not.writeonly) then
            mylevmax = max(mylevmax,tree%node%level)
            call fill_subnode(tree%node%igrid,tree%node%active,jgrid,dir,&
               xslice,normconv)
         end if
         return
      end if
      ! We are out for leaves now, continue for branches

      ! Get the correct child:
      select case (dir)
         case (1)
         ic1 = igslice(tree%node%level+1) - 2 * tree%node%ig1 + 2
         
         case (2)
         ic2 = igslice(tree%node%level+1) - 2 * tree%node%ig2 + 2
         
         case (3)
         ic3 = igslice(tree%node%level+1) - 2 * tree%node%ig3 + 2
         
      case default
         call mpistop("slice direction not clear in traverse_slice")
      end select

      ! Recursively descend into the correct child branch:
      
      select case(dir)
      case (1)
         do ic3=1,2
            do ic2=1,2
               call traverse_slice(tree%node%child(ic1,ic2,ic3))
            end do
         end do
      case (2)
         do ic3=1,2
            do ic1=1,2
               call traverse_slice(tree%node%child(ic1,ic2,ic3))
            end do
         end do
      case (3)
         do ic2=1,2
            do ic1=1,2
               call traverse_slice(tree%node%child(ic1,ic2,ic3))
            end do
         end do
      case default
         call mpistop("slice direction not clear in traverse_slice")
      end select
     
      
      

    end subroutine traverse_slice

  end subroutine select_slice

  subroutine fill_subnode(igrid,active,jgrid,dir,xslice,normconv)
    use mod_global_parameters
    use mod_calculate_xw
    integer, intent(in)                                       :: igrid, dir
    integer, intent(inout)                                    :: jgrid
    logical, intent(in)                                       :: active
    double precision, intent(in)                              :: xslice
    double precision,dimension(0:nw+nwauxio),intent(out)      :: normconv 
    ! .. local ..
    double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
       ixMlo3-1:ixMhi3,ndim)       :: xC_TMP, xC
    double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
       ndim)         :: xCC_TMP, xCC
    double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
       ixMlo3-1:ixMhi3,nw+nwauxio) :: wC_TMP
    double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
       nw+nwauxio)   :: wCC_TMP
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)                        :: x_save
    integer                :: ixslice, nwexpand, ixCmin1,ixCmin2,ixCmin3,&
       ixCmax1,ixCmax2,ixCmax3, ixCCmin1,ixCCmin2,ixCCmin3,ixCCmax1,ixCCmax2,&
       ixCCmax3
    logical                :: mask(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)

    mask=.false.
    ! increase grid-count:
    jgrid=jgrid+1
    ! Allocate subdim solution array:
    if (ndim==1) then
       nwexpand = nwauxio
    else
       if(slice_type/='dat')then
          nwexpand = nwauxio
       else
          nwexpand = 0
       endif
    end if
    call alloc_subnode(jgrid,dir,nwexpand)
    call fill_subnode_info(igrid,jgrid,dir)

    mask(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=.true.

    ! Now hunt for the index closest to the slice:
    
    
    
    select case (dir)
    case (1)
       ixslice = minloc(dabs(xslice-ps(igrid)%x(:,ixMlo2,ixMlo3,dir)),1,mask(:,&
          ixMlo2,ixMlo3))
    case (2)
       ixslice = minloc(dabs(xslice-ps(igrid)%x(ixMlo1,:,ixMlo3,dir)),1,&
          mask(ixMlo1,:,ixMlo3))
    case (3)
       ixslice = minloc(dabs(xslice-ps(igrid)%x(ixMlo1,ixMlo2,:,dir)),1,&
          mask(ixMlo1,ixMlo2,:))
    case default
       call mpistop("slice direction not clear in fill_subnode")
    end select
   

    call calc_x(igrid,xC,xCC)
    ! Set the coordinate to be exactly on the slice:
    xC(:,:,:,dir)  = xslice
    xCC(:,:,:,dir) = xslice
    call calc_grid(unitslice,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,&
       normconv,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,&
       ixCCmin2,ixCCmin3,ixCCmax1,ixCCmax2,ixCCmax3,.true.)
    ! CC stands for CellCenter
    ! C  stands for Corner
    
    ! Fill the subdimensional solution and position array:
    
    
    
    select case (dir)
    case (1)
       ps_sub(jgrid)%w(ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,&
          1:nw+nwexpand) = wCC_TMP(ixslice,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,&
          1:nw+nwexpand)
       ps_sub(jgrid)%x(ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,&
          1:ndim) = xCC_TMP(ixslice,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,&
          1:ndim)
       ps_sub(jgrid)%wC(ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
          1:nw+nwexpand) = wC_TMP(ixslice,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
          1:nw+nwexpand)
       ps_sub(jgrid)%xC(ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
          1:ndim) = xC_TMP(ixslice,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:ndim)
    case (2)
       ps_sub(jgrid)%w(ixCCmin1:ixCCmax1,ixCCmin3:ixCCmax3,&
          1:nw+nwexpand) = wCC_TMP(ixCCmin1:ixCCmax1,ixslice,ixCCmin3:ixCCmax3,&
          1:nw+nwexpand)
       ps_sub(jgrid)%x(ixCCmin1:ixCCmax1,ixCCmin3:ixCCmax3,&
          1:ndim) = xCC_TMP(ixCCmin1:ixCCmax1,ixslice,ixCCmin3:ixCCmax3,&
          1:ndim)
       ps_sub(jgrid)%wC(ixCmin1:ixCmax1,ixCmin3:ixCmax3,&
          1:nw+nwexpand) = wC_TMP(ixCmin1:ixCmax1,ixslice,ixCmin3:ixCmax3,&
          1:nw+nwexpand)
       ps_sub(jgrid)%xC(ixCmin1:ixCmax1,ixCmin3:ixCmax3,&
          1:ndim) = xC_TMP(ixCmin1:ixCmax1,ixslice,ixCmin3:ixCmax3,1:ndim) 
    case (3)
       ps_sub(jgrid)%w(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,&
          1:nw+nwexpand) = wCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixslice,&
          1:nw+nwexpand) 
       ps_sub(jgrid)%x(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,&
          1:ndim) = xCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixslice,&
          1:ndim)
       ps_sub(jgrid)%wC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
          1:nw+nwexpand) = wC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixslice,&
          1:nw+nwexpand) 
       ps_sub(jgrid)%xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
          1:ndim) = xC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixslice,1:ndim) 
    case default
       call mpistop("slice direction not clear in fill_subnode")
    end select
   

  end subroutine fill_subnode

  subroutine alloc_subnode(jgrid,dir,nwexpand)
    use mod_global_parameters
    integer, intent(in) :: jgrid, dir, nwexpand

    ! take care, what comes out is not necessarily a right handed system!
    
    
    
    select case (dir)
    case (1)
       allocate(ps_sub(jgrid)%w(ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw+nwexpand),&
          ps_sub(jgrid)%x(ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim),&
           ps_sub(jgrid)%wC(ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw+nwexpand),&
          ps_sub(jgrid)%xC(ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim))
    case (2)
       allocate(ps_sub(jgrid)%w(ixGlo1:ixGhi1,ixGlo3:ixGhi3,1:nw+nwexpand),&
          ps_sub(jgrid)%x(ixGlo1:ixGhi1,ixGlo3:ixGhi3,1:ndim),&
           ps_sub(jgrid)%wC(ixGlo1:ixGhi1,ixGlo3:ixGhi3,1:nw+nwexpand),&
          ps_sub(jgrid)%xC(ixGlo1:ixGhi1,ixGlo3:ixGhi3,1:ndim))
    case (3)
       allocate(ps_sub(jgrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw+nwexpand),&
          ps_sub(jgrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim),&
           ps_sub(jgrid)%wC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw+nwexpand),&
          ps_sub(jgrid)%xC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim))
    case default
       call mpistop("slice direction not clear in alloc_subnode")
    end select
   
  end subroutine alloc_subnode

  subroutine dealloc_subnode(jgrid)
    use mod_global_parameters
    integer, intent(in) :: jgrid

    if (jgrid==0) then
       call mpistop("trying to delete a non-existing grid in dealloc_subnode")
    end if

    deallocate(ps_sub(jgrid)%w,ps_sub(jgrid)%x,ps_sub(jgrid)%wC,&
       ps_sub(jgrid)%xC)

    ! reset the global node info:
    node_sub(:,jgrid)=0
    rnode_sub(:,jgrid)=zero

  end subroutine dealloc_subnode

  subroutine fill_subnode_info(igrid,jgrid,dir)
    use mod_global_parameters
    integer, intent(in) :: igrid,jgrid,dir

    node_sub(plevel_,jgrid)=node(plevel_,igrid)
    
    select case(dir)
    case (1)
       node_sub(pig1_,jgrid)=node(pig2_,igrid)
       node_sub(pig2_,jgrid)=node(pig3_,igrid)
       rnode_sub(rpdx1_,jgrid)=rnode(rpdx2_,igrid)
       rnode_sub(rpdx2_,jgrid)=rnode(rpdx3_,igrid)
       rnode_sub(rpxmin1_,jgrid)=rnode(rpxmin2_,igrid)
       rnode_sub(rpxmin2_,jgrid)=rnode(rpxmin3_,igrid)
    case (2)
       node_sub(pig1_,jgrid)=node(pig1_,igrid)
       node_sub(pig2_,jgrid)=node(pig3_,igrid)
       rnode_sub(rpdx1_,jgrid)=rnode(rpdx1_,igrid)
       rnode_sub(rpdx2_,jgrid)=rnode(rpdx3_,igrid)
       rnode_sub(rpxmin1_,jgrid)=rnode(rpxmin1_,igrid)
       rnode_sub(rpxmin2_,jgrid)=rnode(rpxmin3_,igrid)
    case (3)
       node_sub(pig1_,jgrid)=node(pig1_,igrid)
       node_sub(pig2_,jgrid)=node(pig2_,igrid)
       rnode_sub(rpdx1_,jgrid)=rnode(rpdx1_,igrid)
       rnode_sub(rpdx2_,jgrid)=rnode(rpdx2_,igrid)
       rnode_sub(rpxmin1_,jgrid)=rnode(rpxmin1_,igrid)
       rnode_sub(rpxmin2_,jgrid)=rnode(rpxmin2_,igrid)
    case default
       call mpistop("slice direction not clear in fill_subnode_info")
    end select
   
    
    

  end subroutine fill_subnode_info

  subroutine get_igslice(dir,x,igslice)
    use mod_global_parameters
    integer, intent(in) :: dir
    double precision, intent(in) :: x
    integer, dimension(nlevelshi), intent(out) :: igslice
    ! .. local ..
    integer :: level
    double precision :: distance
    !double precision :: xsgrid(ndim,nlevelshi),qs(ndim),xmgrid(ndim,3),xnew,xlgrid(ndim,3)
    !integer :: nbefore

    if (x.ne.x) call mpistop("get_igslice: your slice position is NaN!")

    select case (dir)
    case (1)
      select case (stretch_type(1))
      case (stretch_none) ! This dimension is unstretched
        if (x<=xprobmin1) then
          igslice=1
        else if (x>=xprobmax1) then
          do level=1,refine_max_level
            igslice(level)=ng1(level)
          end do
        else
          distance=x-xprobmin1
          do level = 1, refine_max_level
            igslice(level) = ceiling(distance/dg1(level))
          end do
        end if

      case (stretch_uni) ! Uniform stretching

        if (x<=xprobmin1) then
          igslice=1
        else if (x>=xprobmax1) then
          do level=1,refine_max_level
            igslice(level)=ng1(level)
          end do
        else
          distance=x-xprobmin1
          do level=1,refine_max_level
            igslice(level)= ceiling(dlog(distance/dxfirst(level,&
               1)*(qstretch(level,1)-1.d0)+1.d0)/dlog(qstretch(level,&
               1))/dble(block_nx1))
          end do
        end if

      case (stretch_symm) ! Symmetric stretching
        ! symmetric stretch about 0.5*(xprobmin+xprobmax)
        if (x<=xprobmin1) then
          igslice=1
        else if (x>=xprobmax1) then
          do level=1,refine_max_level
            igslice(level)=ng1(level)
          end do
        else if(x<xprobmin1+xstretch1) then
          distance=xprobmin1+xstretch1-x
          ! stretch to left from xprobmin+xstretch
          do level=1,refine_max_level
            igslice(level) = nstretchedblocks(level,&
               1)/2-int(dlog(distance*(qstretch(level,1)-one)/dxfirst(level,&
               1)+one)/dlog(qstretch(level,1))/dble(block_nx1))
          end do
        else if(x>xprobmax1-xstretch1) then
          distance=x-xprobmax1+xstretch1
          ! stretch to right from xprobmax-xstretch
          do level=1,refine_max_level
            igslice(level) = ceiling(dlog(distance*(qstretch(level,&
               1)-one)/dxfirst(level,1)+one)/dlog(qstretch(level,&
               1))/dble(block_nx1))+ng1(level)-nstretchedblocks(level,1)/2
          end do
        else
          ! possible non-stretched central part
          distance=x-xprobmin1-xstretch1
          do level=1,refine_max_level
            igslice(level)=nstretchedblocks(level,&
               1)/2+ceiling(distance/dg1(level))
          end do
        end if
 !xsgrid(1,1)=half*(xprobmax1-xprobmin1)!        /(half*domain_nx1-half*nstretchedblocks_baselevel(1)*block_nx1 !        +(one-qstretch_baselevel(1)**(half*nstretchedblocks_baselevel(1)*block_nx1)) !        /(one-qstretch_baselevel(1)))
 !xmgrid(1,1)=xsgrid(1,1)*(domain_nx1-nstretchedblocks_baselevel(1)*block_nx1)
        !xlgrid(1,1)=xsgrid(1,1)

        !if (x .ge. xmgrid(1,1)*half) then
        !  xnew=x-xmgrid(1,1)*half
        !  do level=1,refine_max_level
        !    qs(1)=qstretch_baselevel(1)**(one/2.d0**(level-1))
        !    if (level .gt. 1) xsgrid(1,level)=xsgrid(1,level-1)/(one+qs(1))
 !nbefore=(domain_nx1/block_nx1-nstretchedblocks_baselevel(1)/2)*2**(level-1)
 !igslice(level)=floor((dlog(one-xnew/xsgrid(1,level)*(one-qs(1)))/dlog(qs(1))-1.e-16)/block_nx1)+1+nbefore
        !    if (x>=xprobmax1) igslice(level)=domain_nx1*2**(level-1)
        !  end do
        !else if(x .ge. -xmgrid(1,1)*half) then
        !  xnew=x+xmgrid(1,1)*half
        !  do level=1,refine_max_level
        !    nbefore=nstretchedblocks_baselevel(1)/2*2**(level-1)
        !    if (level .gt. 1) xlgrid(1,level)=xlgrid(1,level-1)/2.d0
 !igslice(level)=floor(((xnew-1.e-16)/xlgrid(1,level))/block_nx1)+1+nbefore
        !  end do
        !else
        !  xnew=xmgrid(1,1)*half-x
        !  do level=1,refine_max_level
        !    qs(1)=qstretch_baselevel(1)**(one/2.d0**(level-1))
        !    if (level .gt. 1) xsgrid(1,level)=xsgrid(1,level-1)/(one+qs(1))
 !igslice(level)=nstretchedblocks_baselevel(1)/2*2**(level-1) !      -floor((dlog(one-xnew/xsgrid(1,level)*(one-qs(1)))/dlog(qs(1)))/block_nx1)
        !    if (x<=xprobmin1) igslice(level)=1
        !  end do
        !end if

      case default
        call mpistop("stretch type not supported by get_igslice")
      end select 
    
    case (2)
      select case (stretch_type(2))
      case (stretch_none) ! This dimension is unstretched
        if (x<=xprobmin2) then
          igslice=1
        else if (x>=xprobmax2) then
          do level=1,refine_max_level
            igslice(level)=ng2(level)
          end do
        else
          distance=x-xprobmin2
          do level = 1, refine_max_level
            igslice(level) = ceiling(distance/dg2(level))
          end do
        end if

      case (stretch_uni) ! Uniform stretching

        if (x<=xprobmin2) then
          igslice=1
        else if (x>=xprobmax2) then
          do level=1,refine_max_level
            igslice(level)=ng2(level)
          end do
        else
          distance=x-xprobmin2
          do level=1,refine_max_level
            igslice(level)= ceiling(dlog(distance/dxfirst(level,&
               2)*(qstretch(level,2)-1.d0)+1.d0)/dlog(qstretch(level,&
               2))/dble(block_nx2))
          end do
        end if

      case (stretch_symm) ! Symmetric stretching
        ! symmetric stretch about 0.5*(xprobmin+xprobmax)
        if (x<=xprobmin2) then
          igslice=1
        else if (x>=xprobmax2) then
          do level=1,refine_max_level
            igslice(level)=ng2(level)
          end do
        else if(x<xprobmin2+xstretch2) then
          distance=xprobmin2+xstretch2-x
          ! stretch to left from xprobmin+xstretch
          do level=1,refine_max_level
            igslice(level) = nstretchedblocks(level,&
               2)/2-int(dlog(distance*(qstretch(level,2)-one)/dxfirst(level,&
               2)+one)/dlog(qstretch(level,2))/dble(block_nx2))
          end do
        else if(x>xprobmax2-xstretch2) then
          distance=x-xprobmax2+xstretch2
          ! stretch to right from xprobmax-xstretch
          do level=1,refine_max_level
            igslice(level) = ceiling(dlog(distance*(qstretch(level,&
               2)-one)/dxfirst(level,2)+one)/dlog(qstretch(level,&
               2))/dble(block_nx2))+ng2(level)-nstretchedblocks(level,2)/2
          end do
        else
          ! possible non-stretched central part
          distance=x-xprobmin2-xstretch2
          do level=1,refine_max_level
            igslice(level)=nstretchedblocks(level,&
               2)/2+ceiling(distance/dg2(level))
          end do
        end if
 !xsgrid(2,1)=half*(xprobmax2-xprobmin2)!        /(half*domain_nx2-half*nstretchedblocks_baselevel(2)*block_nx2 !        +(one-qstretch_baselevel(2)**(half*nstretchedblocks_baselevel(2)*block_nx2)) !        /(one-qstretch_baselevel(2)))
 !xmgrid(2,1)=xsgrid(2,1)*(domain_nx2-nstretchedblocks_baselevel(2)*block_nx2)
        !xlgrid(2,1)=xsgrid(2,1)

        !if (x .ge. xmgrid(2,1)*half) then
        !  xnew=x-xmgrid(2,1)*half
        !  do level=1,refine_max_level
        !    qs(2)=qstretch_baselevel(2)**(one/2.d0**(level-1))
        !    if (level .gt. 1) xsgrid(2,level)=xsgrid(2,level-1)/(one+qs(2))
 !nbefore=(domain_nx2/block_nx2-nstretchedblocks_baselevel(2)/2)*2**(level-1)
 !igslice(level)=floor((dlog(one-xnew/xsgrid(2,level)*(one-qs(2)))/dlog(qs(2))-1.e-16)/block_nx2)+1+nbefore
        !    if (x>=xprobmax2) igslice(level)=domain_nx2*2**(level-1)
        !  end do
        !else if(x .ge. -xmgrid(2,1)*half) then
        !  xnew=x+xmgrid(2,1)*half
        !  do level=1,refine_max_level
        !    nbefore=nstretchedblocks_baselevel(2)/2*2**(level-1)
        !    if (level .gt. 1) xlgrid(2,level)=xlgrid(2,level-1)/2.d0
 !igslice(level)=floor(((xnew-1.e-16)/xlgrid(2,level))/block_nx2)+1+nbefore
        !  end do
        !else
        !  xnew=xmgrid(2,1)*half-x
        !  do level=1,refine_max_level
        !    qs(2)=qstretch_baselevel(2)**(one/2.d0**(level-1))
        !    if (level .gt. 1) xsgrid(2,level)=xsgrid(2,level-1)/(one+qs(2))
 !igslice(level)=nstretchedblocks_baselevel(2)/2*2**(level-1) !      -floor((dlog(one-xnew/xsgrid(2,level)*(one-qs(2)))/dlog(qs(2)))/block_nx2)
        !    if (x<=xprobmin2) igslice(level)=1
        !  end do
        !end if

      case default
        call mpistop("stretch type not supported by get_igslice")
      end select 
    
    case (3)
      select case (stretch_type(3))
      case (stretch_none) ! This dimension is unstretched
        if (x<=xprobmin3) then
          igslice=1
        else if (x>=xprobmax3) then
          do level=1,refine_max_level
            igslice(level)=ng3(level)
          end do
        else
          distance=x-xprobmin3
          do level = 1, refine_max_level
            igslice(level) = ceiling(distance/dg3(level))
          end do
        end if

      case (stretch_uni) ! Uniform stretching

        if (x<=xprobmin3) then
          igslice=1
        else if (x>=xprobmax3) then
          do level=1,refine_max_level
            igslice(level)=ng3(level)
          end do
        else
          distance=x-xprobmin3
          do level=1,refine_max_level
            igslice(level)= ceiling(dlog(distance/dxfirst(level,&
               3)*(qstretch(level,3)-1.d0)+1.d0)/dlog(qstretch(level,&
               3))/dble(block_nx3))
          end do
        end if

      case (stretch_symm) ! Symmetric stretching
        ! symmetric stretch about 0.5*(xprobmin+xprobmax)
        if (x<=xprobmin3) then
          igslice=1
        else if (x>=xprobmax3) then
          do level=1,refine_max_level
            igslice(level)=ng3(level)
          end do
        else if(x<xprobmin3+xstretch3) then
          distance=xprobmin3+xstretch3-x
          ! stretch to left from xprobmin+xstretch
          do level=1,refine_max_level
            igslice(level) = nstretchedblocks(level,&
               3)/2-int(dlog(distance*(qstretch(level,3)-one)/dxfirst(level,&
               3)+one)/dlog(qstretch(level,3))/dble(block_nx3))
          end do
        else if(x>xprobmax3-xstretch3) then
          distance=x-xprobmax3+xstretch3
          ! stretch to right from xprobmax-xstretch
          do level=1,refine_max_level
            igslice(level) = ceiling(dlog(distance*(qstretch(level,&
               3)-one)/dxfirst(level,3)+one)/dlog(qstretch(level,&
               3))/dble(block_nx3))+ng3(level)-nstretchedblocks(level,3)/2
          end do
        else
          ! possible non-stretched central part
          distance=x-xprobmin3-xstretch3
          do level=1,refine_max_level
            igslice(level)=nstretchedblocks(level,&
               3)/2+ceiling(distance/dg3(level))
          end do
        end if
 !xsgrid(3,1)=half*(xprobmax3-xprobmin3)!        /(half*domain_nx3-half*nstretchedblocks_baselevel(3)*block_nx3 !        +(one-qstretch_baselevel(3)**(half*nstretchedblocks_baselevel(3)*block_nx3)) !        /(one-qstretch_baselevel(3)))
 !xmgrid(3,1)=xsgrid(3,1)*(domain_nx3-nstretchedblocks_baselevel(3)*block_nx3)
        !xlgrid(3,1)=xsgrid(3,1)

        !if (x .ge. xmgrid(3,1)*half) then
        !  xnew=x-xmgrid(3,1)*half
        !  do level=1,refine_max_level
        !    qs(3)=qstretch_baselevel(3)**(one/2.d0**(level-1))
        !    if (level .gt. 1) xsgrid(3,level)=xsgrid(3,level-1)/(one+qs(3))
 !nbefore=(domain_nx3/block_nx3-nstretchedblocks_baselevel(3)/2)*2**(level-1)
 !igslice(level)=floor((dlog(one-xnew/xsgrid(3,level)*(one-qs(3)))/dlog(qs(3))-1.e-16)/block_nx3)+1+nbefore
        !    if (x>=xprobmax3) igslice(level)=domain_nx3*2**(level-1)
        !  end do
        !else if(x .ge. -xmgrid(3,1)*half) then
        !  xnew=x+xmgrid(3,1)*half
        !  do level=1,refine_max_level
        !    nbefore=nstretchedblocks_baselevel(3)/2*2**(level-1)
        !    if (level .gt. 1) xlgrid(3,level)=xlgrid(3,level-1)/2.d0
 !igslice(level)=floor(((xnew-1.e-16)/xlgrid(3,level))/block_nx3)+1+nbefore
        !  end do
        !else
        !  xnew=xmgrid(3,1)*half-x
        !  do level=1,refine_max_level
        !    qs(3)=qstretch_baselevel(3)**(one/2.d0**(level-1))
        !    if (level .gt. 1) xsgrid(3,level)=xsgrid(3,level-1)/(one+qs(3))
 !igslice(level)=nstretchedblocks_baselevel(3)/2*2**(level-1) !      -floor((dlog(one-xnew/xsgrid(3,level)*(one-qs(3)))/dlog(qs(3)))/block_nx3)
        !    if (x<=xprobmin3) igslice(level)=1
        !  end do
        !end if

      case default
        call mpistop("stretch type not supported by get_igslice")
      end select 
    
    case default
      call mpistop("slice direction not clear in get_igslice")
    end select

  end subroutine get_igslice

  double precision function roundoff_minmax(val,minval,maxval)
    implicit none
    double precision,intent(in)         :: val, minval, maxval

    roundoff_minmax = val

    if (abs(roundoff_minmax) .le. minval) then
       roundoff_minmax = 0.0d0
    end if

    if (roundoff_minmax .gt. maxval) then
       roundoff_minmax = maxval
    else if (roundoff_minmax .lt. -maxval) then
       roundoff_minmax = -maxval
    end if

    ! Replace NaN with maxval (e.g. Paraview chokes on ASCII NaN): 
    if (roundoff_minmax /= roundoff_minmax) roundoff_minmax = maxval

  end function roundoff_minmax

end module mod_slice

