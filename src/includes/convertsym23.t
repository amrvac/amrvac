{^IFTWOD
! subroutines to convert 2.5D data to 3D data
!====================================================================================
subroutine unstructuredvtk23(qunit,userconvert_type)

! output for vtu format to paraview
! not parallel, uses calc_grid to compute nwauxio variables

use mod_global_parameters

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,3) :: xC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,&
   3)   :: xCC_TMP

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,nw+nwauxio)   :: wC_TMP
double precision,dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,nw&
   +nwauxio)     :: wCC_TMP
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo1:ixGhi1,1:nw&
   +nwauxio)   :: w
double precision, dimension(0:nw+nwauxio)                   :: normconv
integer:: igrid,iigrid,level,igonlevel,icel,ixCmin1,ixCmin2,&
   ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,ixCCmin2,ixCCmin3,ixCCmax1,&
   ixCCmax2,ixCCmax3,iw
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx1,nx2,nx3,nxC1,nxC2,nxC3,nodesonlevel,elemsonlevel,nc,np,&
   VTK_type,ix1,ix2,ix3
integer :: i3grid,n3grid
double precision ::d3grid

character(len=80)::  filename
character(len=20):: userconvert_type
character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:3+nw+nwauxio)
character(len=1024) :: outfilehead
double precision :: zlength,zlengsc,zgridsc

logical ::      fileopen
!-----------------------------------------------------------------------------
if(npe>1)then
 if(mype==0) PRINT *,'unstructuredvtk23 not parallel, use vtumpi'
 call mpistop('npe>1, unstructuredvtk23')
end if

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
  ! generate filename 
  write(filename,'(a,a,i4.4,a)') TRIM(base_filename),"3D",snapshotini,".vtu"
  ! Open the file for the header part
  open(qunit,file=filename,status='unknown',form='formatted')
endif

call getheadernames(wnamei,xandwnamei,outfilehead)

! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'<UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(qunit,'(f10.2)') real(global_time*time_convert_factor)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'

! number of cells, number of corner points, per grid.
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi1-ixMlo1+1;
nxC1=nx1+1;nxC2=nx2+1;nxC3=nx3+1;
nc=nx1*nx2*nx3
np=nxC1*nxC2*nxC3

! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
zgridsc=2.d0
zlengsc=2.d0*zgridsc
zlength=zlengsc*(xprobmax1-xprobmin1)
do level=levmin,levmax
 if (writelevel(level)) then
   do iigrid=1,igridstail; igrid=igrids(iigrid);
    if (node(plevel_,igrid)/=level) cycle
    block=>ps(igrid)
    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    if ((rnode(rpxmin1_,igrid)>=xprobmin1+(xprobmax1-xprobmin1)&
       *writespshift(1,1).and.rnode(rpxmin2_,igrid)>=xprobmin2&
       +(xprobmax2-xprobmin2)*writespshift(2,1))&
       .and.(rnode(rpxmax1_,igrid)<=xprobmax1-(xprobmax1-xprobmin1)&
       *writespshift(1,2).and.rnode(rpxmax2_,igrid)<=xprobmax2-(xprobmax2&
       -xprobmin2)*writespshift(2,2))) then
      d3grid=zgridsc*(rnode(rpxmax1_,igrid)-rnode(rpxmin1_,igrid))
      n3grid=dnint(zlength/d3grid)
      ! In case primitives to be saved: use primitive subroutine
      !  extra layer around mesh only needed when storing corner values and averaging
      if(saveprim) then
       call primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
        ixGlo1,ixGlo2,ixGhi1,ixGhi2,ps(igrid)%w,ps(igrid)%x)
      endif
      ! using array w so that new output auxiliaries can be calculated by the user
      ! extend 2D data to 3D insuring variables are independent on the third coordinate
      {^IFTWOD
      do ix3=ixGlo1,ixGhi1
       w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix3,1:nw)=ps(igrid)%w(ixGlo1:ixGhi1,&
         ixGlo2:ixGhi2,1:nw)
      end do
      \}
      do i3grid=1,n3grid
       call calc_grid23(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
          ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,ixCCmin2,&
          ixCCmin3,ixCCmax1,ixCCmax2,ixCCmax3,.true.,i3grid,d3grid,w,zlength,zgridsc)
       select case(userconvert_type)
        case('vtu23')
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
             '" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<PointData>'
          do iw=1,nw
             if(iw<=nw.and.(.not.w_write(iw)))cycle

             write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                TRIM(wnamei(iw)),'" format="ascii">'
             write(qunit,'(200(1pe14.6))') (((wC_TMP(ix1,ix2,ix3,iw)&
                *normconv(iw),ix1=ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2),ix3&
                =ixCmin3,ixCmax3)
             write(qunit,'(a)')'</DataArray>'
          enddo
          if(nwauxio>0) then
           do iw=nw+1,nw+nwauxio
              write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                 TRIM(wnamei(iw)),'" format="ascii">'
              write(qunit,'(200(1pe14.6))') (((wC_TMP(ix1,ix2,ix3,iw)&
                 *normconv(iw),ix1=ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2),ix3&
                 =ixCmin3,ixCmax3)
              write(qunit,'(a)')'</DataArray>'
           enddo
          endif
          write(qunit,'(a)')'</PointData>'

          write(qunit,'(a)')'<Points>'
          write(qunit,'(a)')&
             '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          ! write cell corner coordinates in a backward dimensional loop, always
          ! 3D output
          do ix3=ixCmin3,ixCmax3
          do ix2=ixCmin2,ixCmax2
          do ix1=ixCmin1,ixCmax1
             x_VTK(1:3)=zero;
             x_VTK(1:3)=xC_TMP(ix1,ix2,ix3,1:3)*normconv(0);
             write(qunit,'(3(1pe14.6))') x_VTK
          end do
          end do
          end do
          write(qunit,'(a)')'</DataArray>'
          write(qunit,'(a)')'</Points>'
        case('vtuCC23')
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
             '" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<CellData>'
          do iw=1,nw
             if(iw<=nw.and.(.not.w_write(iw)))cycle

             write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                TRIM(wnamei(iw)),'" format="ascii">'
             write(qunit,'(200(1pe14.6))') (((wCC_TMP(ix1,ix2,ix3,iw)&
                *normconv(iw),ix1=ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2),&
                ix3=ixCCmin3,ixCCmax3)
             write(qunit,'(a)')'</DataArray>'
          enddo
          if(nwauxio>0) then
           do iw=nw+1,nw+nwauxio
              write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                 TRIM(wnamei(iw)),'" format="ascii">'
              write(qunit,'(200(1pe14.6))') (((wCC_TMP(ix1,ix2,ix3,iw)&
                 *normconv(iw),ix1=ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2),&
                 ix3=ixCCmin3,ixCCmax3)
              write(qunit,'(a)')'</DataArray>'
           enddo
          endif
          write(qunit,'(a)')'</CellData>'

          write(qunit,'(a)')'<Points>'
          write(qunit,'(a)')&
             '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          ! write cell corner coordinates in a backward dimensional loop, always
          ! 3D output
          do ix3=ixCmin3,ixCmax3
          do ix2=ixCmin2,ixCmax2
          do ix1=ixCmin1,ixCmax1
             x_VTK(1:3)=zero;
             x_VTK(1:3)=xC_TMP(ix1,ix2,ix3,1:3)*normconv(0);
             write(qunit,'(3(1pe14.6))') x_VTK
          end do
          end do
          end do
          write(qunit,'(a)')'</DataArray>'
          write(qunit,'(a)')'</Points>'
       end select

       write(qunit,'(a)')'<Cells>'

       ! connectivity part
       write(qunit,'(a)')&
          '<DataArray type="Int32" Name="connectivity" format="ascii">'
       call save_connvtk23(qunit,igrid)
       write(qunit,'(a)')'</DataArray>'

       ! offsets data array
       write(qunit,'(a)')&
          '<DataArray type="Int32" Name="offsets" format="ascii">'
       do icel=1,nc
          write(qunit,'(i7)') icel*(2**3)
       end do
       write(qunit,'(a)')'</DataArray>'

       ! VTK cell type data array
       write(qunit,'(a)')'<DataArray type="Int32" Name="types" format="ascii">'
       ! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax


        VTK_type=11
       do icel=1,nc
          write(qunit,'(i2)') VTK_type
       enddo
       write(qunit,'(a)')'</DataArray>'

       write(qunit,'(a)')'</Cells>'

       write(qunit,'(a)')'</Piece>'
      enddo
    endif
   enddo
 endif
enddo
write(qunit,'(a)')'</UnstructuredGrid>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)
end subroutine unstructuredvtk23
!====================================================================================
subroutine unstructuredvtksym23(qunit,userconvert_type)

! output for vtu format to paraview
! not parallel, uses calc_grid to compute nwauxio variables
!
! use this subroutine  when the physical domain is symmetric/asymmetric about (0,y,z) 
! plane, xprobmin1=0 and the computational domain is a half of the physical domain

use mod_global_parameters

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,3) :: xC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,&
   3)   :: xCC_TMP

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,nw+nwauxio)   :: wC_TMP
double precision,dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,nw&
   +nwauxio)     :: wCC_TMP
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo1:ixGhi1,1:nw&
   +nwauxio)   :: w
double precision, dimension(0:nw+nwauxio)                   :: normconv
integer:: igrid,iigrid,level,igonlevel,icel,ixCmin1,ixCmin2,&
   ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,ixCCmin2,ixCCmin3,ixCCmax1,&
   ixCCmax2,ixCCmax3,iw
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx1,nx2,nx3,nxC1,nxC2,nxC3,nodesonlevel,elemsonlevel,nc,np,&
   VTK_type,ix1,ix2,ix3
integer :: i3grid,n3grid
double precision ::d3grid,zlengsc,zgridsc

character(len=80)::  filename
character(len=20):: userconvert_type
character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:3+nw+nwauxio)
character(len=1024) :: outfilehead

double precision :: zlength

logical ::      fileopen
!-----------------------------------------------------------------------------
if(npe>1)then
 if(mype==0) PRINT *,'unstructuredvtksym23 not parallel, use vtumpi'
 call mpistop('npe>1, unstructuredvtksym23')
end if

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
  ! generate filename 
  write(filename,'(a,a,i4.4,a)') TRIM(base_filename),"3D",snapshotini,".vtu"
  ! Open the file for the header part
  open(qunit,file=filename,status='unknown',form='formatted')
endif

call getheadernames(wnamei,xandwnamei,outfilehead)

! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'<UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(qunit,'(f10.2)') real(global_time*time_convert_factor)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'

! number of cells, number of corner points, per grid.
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi1-ixMlo1+1;
nxC1=nx1+1;nxC2=nx2+1;nxC3=nx3+1;
nc=nx1*nx2*nx3
np=nxC1*nxC2*nxC3

! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
zgridsc=2.d0
zlengsc=2.d0*zgridsc
zlength=zlengsc*(xprobmax1-xprobmin1)
do level=levmin,levmax
 if (writelevel(level)) then
   do iigrid=1,igridstail; igrid=igrids(iigrid);
    if (node(plevel_,igrid)/=level) cycle
    block=>ps(igrid)
    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    if ((rnode(rpxmin1_,igrid)>=xprobmin1+(xprobmax1-xprobmin1)&
       *writespshift(1,1).and.rnode(rpxmin2_,igrid)>=xprobmin2&
       +(xprobmax2-xprobmin2)*writespshift(2,1))&
       .and.(rnode(rpxmax1_,igrid)<=xprobmax1-(xprobmax1-xprobmin1)&
       *writespshift(1,2).and.rnode(rpxmax2_,igrid)<=xprobmax2-(xprobmax2&
       -xprobmin2)*writespshift(2,2))) then
      d3grid=zgridsc*(rnode(rpxmax1_,igrid)-rnode(rpxmin1_,igrid))
      n3grid=dnint(zlength/d3grid)
      ! In case primitives to be saved: use primitive subroutine
      !  extra layer around mesh only needed when storing corner values and averaging
      if(saveprim) then
       call primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
        ixGlo1,ixGlo2,ixGhi1,ixGhi2,ps(igrid)%w,ps(igrid)%x)
      endif
      ! using array w so that new output auxiliaries can be calculated by the user
      ! extend 2D data to 3D insuring variables are independent on the third coordinate
      do ix3=ixGlo1,ixGhi1
       w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix3,1:nw)=ps(igrid)%w(ixGlo1:ixGhi1,&
         ixGlo2:ixGhi2,1:nw)
      end do
      do i3grid=1,n3grid
       call calc_grid23(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
          ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,ixCCmin2,&
          ixCCmin3,ixCCmax1,ixCCmax2,ixCCmax3,.true.,i3grid,d3grid,w,zlength,zgridsc)
       !! original domain ----------------------------------start 
       select case(userconvert_type)
        case('vtusym23')
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
             '" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<PointData>'
          do iw=1,nw
             if(iw<=nw.and.(.not.w_write(iw)))cycle

             write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                TRIM(wnamei(iw)),'" format="ascii">'
             write(qunit,'(200(1pe14.6))') (((wC_TMP(ix1,ix2,ix3,iw)&
                *normconv(iw),ix1=ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2),ix3&
                =ixCmin3,ixCmax3)
             write(qunit,'(a)')'</DataArray>'
          enddo
          if(nwauxio>0) then
           do iw=nw+1,nw+nwauxio
              write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                 TRIM(wnamei(iw)),'" format="ascii">'
              write(qunit,'(200(1pe14.6))') (((wC_TMP(ix1,ix2,ix3,iw)&
                 *normconv(iw),ix1=ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2),ix3&
                 =ixCmin3,ixCmax3)
              write(qunit,'(a)')'</DataArray>'
           enddo
          endif
          write(qunit,'(a)')'</PointData>'

          write(qunit,'(a)')'<Points>'
          write(qunit,'(a)')&
             '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          ! write cell corner coordinates in a backward dimensional loop, always
          ! 3D output
          do ix3=ixCmin3,ixCmax3
          do ix2=ixCmin2,ixCmax2
          do ix1=ixCmin1,ixCmax1
             x_VTK(1:3)=zero;
             x_VTK(1:3)=xC_TMP(ix1,ix2,ix3,1:3)*normconv(0);
             write(qunit,'(3(1pe14.6))') x_VTK
          end do
          end do
          end do
          write(qunit,'(a)')'</DataArray>'
          write(qunit,'(a)')'</Points>'
        case('vtuCCsym23')
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
             '" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<CellData>'
          do iw=1,nw
             if(iw<=nw.and.(.not.w_write(iw)))cycle

             write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                TRIM(wnamei(iw)),'" format="ascii">'
             write(qunit,'(200(1pe14.6))') (((wCC_TMP(ix1,ix2,ix3,iw)&
                *normconv(iw),ix1=ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2),&
                ix3=ixCCmin3,ixCCmax3)
             write(qunit,'(a)')'</DataArray>'
          enddo
          if(nwauxio>0) then
           do iw=nw+1,nw+nwauxio
              write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                 TRIM(wnamei(iw)),'" format="ascii">'
              write(qunit,'(200(1pe14.6))') (((wCC_TMP(ix1,ix2,ix3,iw)&
                 *normconv(iw),ix1=ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2),&
                 ix3=ixCCmin3,ixCCmax3)
              write(qunit,'(a)')'</DataArray>'
           enddo
          endif
          write(qunit,'(a)')'</CellData>'

          write(qunit,'(a)')'<Points>'
          write(qunit,'(a)')&
             '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          ! write cell corner coordinates in a backward dimensional loop, always
          ! 3D output
          do ix3=ixCmin3,ixCmax3
          do ix2=ixCmin2,ixCmax2
          do ix1=ixCmin1,ixCmax1
             x_VTK(1:3)=zero;
             x_VTK(1:3)=xC_TMP(ix1,ix2,ix3,1:3)*normconv(0);
             write(qunit,'(3(1pe14.6))') x_VTK
          end do
          end do
          end do
          write(qunit,'(a)')'</DataArray>'
          write(qunit,'(a)')'</Points>'
       end select

       write(qunit,'(a)')'<Cells>'

       ! connectivity part
       write(qunit,'(a)')&
          '<DataArray type="Int32" Name="connectivity" format="ascii">'
       call save_connvtk23(qunit,igrid)
       write(qunit,'(a)')'</DataArray>'

       ! offsets data array
       write(qunit,'(a)')&
          '<DataArray type="Int32" Name="offsets" format="ascii">'
       do icel=1,nc
          write(qunit,'(i7)') icel*(2**3)
       end do
       write(qunit,'(a)')'</DataArray>'

       ! VTK cell type data array
       write(qunit,'(a)')'<DataArray type="Int32" Name="types" format="ascii">'
       ! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax


        VTK_type=11
       do icel=1,nc
          write(qunit,'(i2)') VTK_type
       enddo
       write(qunit,'(a)')'</DataArray>'

       write(qunit,'(a)')'</Cells>'

       write(qunit,'(a)')'</Piece>'
       !! original domain ----------------------------------end 
       !! symmetric/asymmetric mirror ----------------------start
       select case(userconvert_type)
        case('vtusym23')
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
             '" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<PointData>'
          do iw=1,nw
            if(iw<=nw.and.(.not.w_write(iw)))cycle
            if(iw==2 .or. iw==4 .or. iw==7) then
             wC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)=&
             -wC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
            endif
             write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                TRIM(wnamei(iw)),'" format="ascii">'
             write(qunit,'(200(1pe14.6))') (((wC_TMP(ix1,ix2,ix3,iw)&
                *normconv(iw),ix1=ixCmax1,ixCmin1,-1),ix2=ixCmin2,ixCmax2),ix3&
                =ixCmin3,ixCmax3)
             write(qunit,'(a)')'</DataArray>'
          enddo
          if(nwauxio>0) then
           do iw=nw+1,nw+nwauxio
              write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                 TRIM(wnamei(iw)),'" format="ascii">'
              write(qunit,'(200(1pe14.6))') (((wC_TMP(ix1,ix2,ix3,iw)&
                 *normconv(iw),ix1=ixCmax1,ixCmin1,-1),ix2=ixCmin2,ixCmax2),ix3&
                 =ixCmin3,ixCmax3)
              write(qunit,'(a)')'</DataArray>'
           enddo
          endif
          write(qunit,'(a)')'</PointData>'

          write(qunit,'(a)')'<Points>'
          write(qunit,'(a)')&
             '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          ! write cell corner coordinates in a backward dimensional loop, always
          ! 3D output
          do ix3=ixCmin3,ixCmax3
          do ix2=ixCmin2,ixCmax2
          do ix1=ixCmax1,ixCmin1,-1
             x_VTK(1:3)=zero;
             x_VTK(1:3)=xC_TMP(ix1,ix2,ix3,1:3)*normconv(0);
             x_VTK(1)=-x_VTK(1)
             write(qunit,'(3(1pe14.6))') x_VTK
          end do
          end do
          end do
          write(qunit,'(a)')'</DataArray>'
          write(qunit,'(a)')'</Points>'
        case('vtuCCsym23')
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
             '" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<CellData>'
          do iw=1,nw
            if(iw<=nw.and.(.not.w_write(iw)))cycle
            if(iw==2 .or. iw==4 .or. iw==7) then
             wCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,iw)=&
            -wCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,iw)
            endif

             write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                TRIM(wnamei(iw)),'" format="ascii">'
             write(qunit,'(200(1pe14.6))') (((wCC_TMP(ix1,ix2,ix3,iw)&
                *normconv(iw),ix1=ixCCmax1,ixCCmin1,-1),ix2=ixCCmin2,ixCCmax2),&
                ix3=ixCCmin3,ixCCmax3)
             write(qunit,'(a)')'</DataArray>'
          enddo
          if(nwauxio>0) then
           do iw=nw+1,nw+nwauxio
              write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                 TRIM(wnamei(iw)),'" format="ascii">'
              write(qunit,'(200(1pe14.6))') (((wCC_TMP(ix1,ix2,ix3,iw)&
                 *normconv(iw),ix1=ixCCmax1,ixCCmin1,-1),ix2=ixCCmin2,ixCCmax2),&
                 ix3=ixCCmin3,ixCCmax3)
              write(qunit,'(a)')'</DataArray>'
           enddo
          endif
          write(qunit,'(a)')'</CellData>'

          write(qunit,'(a)')'<Points>'
          write(qunit,'(a)')&
             '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          ! write cell corner coordinates in a backward dimensional loop, always
          ! 3D output
          do ix3=ixCmin3,ixCmax3
          do ix2=ixCmin2,ixCmax2
          do ix1=ixCmax1,ixCmin1,-1
             x_VTK(1:3)=zero;
             x_VTK(1:3)=xC_TMP(ix1,ix2,ix3,1:3)*normconv(0);
             x_VTK(1)=-x_VTK(1)
             write(qunit,'(3(1pe14.6))') x_VTK
          end do
          end do
          end do
          write(qunit,'(a)')'</DataArray>'
          write(qunit,'(a)')'</Points>'
       end select

       write(qunit,'(a)')'<Cells>'

       ! connectivity part
       write(qunit,'(a)')&
          '<DataArray type="Int32" Name="connectivity" format="ascii">'
       call save_connvtk23(qunit,igrid)
       write(qunit,'(a)')'</DataArray>'

       ! offsets data array
       write(qunit,'(a)')&
          '<DataArray type="Int32" Name="offsets" format="ascii">'
       do icel=1,nc
          write(qunit,'(i7)') icel*(2**3)
       end do
       write(qunit,'(a)')'</DataArray>'

       ! VTK cell type data array
       write(qunit,'(a)')'<DataArray type="Int32" Name="types" format="ascii">'
       ! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax


        VTK_type=11
       do icel=1,nc
          write(qunit,'(i2)') VTK_type
       enddo
       write(qunit,'(a)')'</DataArray>'

       write(qunit,'(a)')'</Cells>'

       write(qunit,'(a)')'</Piece>'
       !! symmetric/asymmetric mirror ----------------------end
      enddo
    endif
   enddo
 endif
enddo
write(qunit,'(a)')'</UnstructuredGrid>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)
end subroutine unstructuredvtksym23
!====================================================================================
subroutine unstructuredvtkB23(qunit,userconvert_type)

! output for vtu format to paraview, binary version output
! not parallel, uses calc_grid to compute nwauxio variables

use mod_global_parameters

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,3) :: xC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,&
   3)   :: xCC_TMP

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,nw&
   +nwauxio)     :: wCC_TMP
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo1:ixGhi1,1:nw&
   +nwauxio)   :: w
integer::               igrid,iigrid,level,igonlevel,icel,ixCmin1,ixCmin2,&
   ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,ixCCmin2,ixCCmin3,ixCCmax1,&
   ixCCmax2,ixCCmax3
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx1,nx2,nx3,nxC1,nxC2,nxC3,nodesonlevel,elemsonlevel,nc,np,&
   VTK_type,ix1,ix2,ix3
double precision :: normconv(0:nw+nwauxio)
character(len=80)::  filename
character(len=20):: userconvert_type
character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:3+nw+nwauxio)
character(len=1024) :: outfilehead
integer*8 :: offset
integer :: size_int,size_double,size_length,recsep,k,iw
integer :: length,lengthcc,offset_points,offset_cells, length_coords,&
   length_conn,length_offsets
integer :: i3grid,n3grid
double precision ::d3grid,zlengsc,zgridsc
character::  buffer
character(len=6)::  bufform

double precision :: zlength

logical ::   fileopen
!-----------------------------------------------------------------------------
if(npe>1)then
 if(mype==0) PRINT *,'unstructuredvtkB23 not parallel, use vtumpi'
 call mpistop('npe>1, unstructuredvtkB23')
end if

offset=0
recsep=4
size_double=4
size_length=4
size_int=size_length

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
  ! generate filename 
  write(filename,'(a,a,i4.4,a)') TRIM(base_filename),"3D",snapshotini,".vtu"
  ! Open the file for the header part
  open(qunit,file=filename,status='unknown')
endif

call getheadernames(wnamei,xandwnamei,outfilehead)
! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'<UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(qunit,'(f10.2)') real(global_time*time_convert_factor)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'

! number of cells, number of corner points, per grid.
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi1-ixMlo1+1;
nxC1=nx1+1;nxC2=nx2+1;nxC3=nx3+1;
nc=nx1*nx2*nx3
np=nxC1*nxC2*nxC3

length=np*size_double
lengthcc=nc*size_double

length_coords=3*length
length_conn=2**3*size_int*nc
length_offsets=nc*size_int

! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
zgridsc=2.d0
zlengsc=2.d0*zgridsc
zlength=zlengsc*(xprobmax1-xprobmin1)
do level=levmin,levmax
 if (writelevel(level)) then
   do iigrid=1,igridstail; igrid=igrids(iigrid);
    if (node(plevel_,igrid)/=level) cycle
    block=>ps(igrid)
    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    if ((rnode(rpxmin1_,igrid)>=xprobmin1+(xprobmax1-xprobmin1)&
       *writespshift(1,1).and.rnode(rpxmin2_,igrid)>=xprobmin2&
       +(xprobmax2-xprobmin2)*writespshift(2,1))&
       .and.(rnode(rpxmax1_,igrid)<=xprobmax1-(xprobmax1-xprobmin1)&
       *writespshift(1,2).and.rnode(rpxmax2_,igrid)<=xprobmax2-(xprobmax2&
       -xprobmin2)*writespshift(2,2))) then
      d3grid=zgridsc*(rnode(rpxmax1_,igrid)-rnode(rpxmin1_,igrid))
      n3grid=dnint(zlength/d3grid)
      do i3grid=1,n3grid !subcycles
      select case(userconvert_type)
       case('vtuB23')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
            '" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<PointData>'
         do iw=1,nw
            if(.not.w_write(iw))cycle
            write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
               TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+length+size_length
         enddo
         if(nwauxio>0)then
          do iw=nw+1,nw+nwauxio
             write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
                TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+length+size_length
          enddo
         endif
         write(qunit,'(a)')'</PointData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
            offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_length
         write(qunit,'(a)')'</Points>'
       case('vtuBCC23')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
            '" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<CellData>'
         do iw=1,nw
            if(.not.w_write(iw))cycle
            write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
               TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+lengthcc+size_length
         enddo
         if(nwauxio>0)then
          do iw=nw+1,nw+nwauxio
             write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
                TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+lengthcc+size_length
          enddo
         endif
         write(qunit,'(a)')'</CellData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
            offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_length
         write(qunit,'(a)')'</Points>'
      end select


      write(qunit,'(a)')'<Cells>'

      ! connectivity part
      write(qunit,'(a,i16,a)')&
'<DataArray type="Int32" Name="connectivity" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_conn+size_length

      ! offsets data array
      write(qunit,'(a,i16,a)') &
         '<DataArray type="Int32" Name="offsets" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_offsets+size_length

      ! VTK cell type data array
      write(qunit,'(a,i16,a)') &
         '<DataArray type="Int32" Name="types" format="appended" offset="',&
         offset,'"/>'
      offset=offset+size_length+nc*size_int

      write(qunit,'(a)')'</Cells>'

      write(qunit,'(a)')'</Piece>'
      enddo !subcycles
    endif
   enddo
 endif
enddo

write(qunit,'(a)')'</UnstructuredGrid>'
write(qunit,'(a)')'<AppendedData encoding="raw">'

close(qunit)
open(qunit,file=filename,form='unformatted',access='stream',status='old',position='append')
buffer='_'
write(qunit) TRIM(buffer)

do level=levmin,levmax
 if (writelevel(level)) then
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (node(plevel_,igrid)/=level) cycle
      block=>ps(igrid)
      ! only output a grid when fully within clipped region selected
      ! by writespshift array
      if ((rnode(rpxmin1_,igrid)>=xprobmin1+(xprobmax1-xprobmin1)&
         *writespshift(1,1).and.rnode(rpxmin2_,igrid)>=xprobmin2&
         +(xprobmax2-xprobmin2)*writespshift(2,1))&
         .and.(rnode(rpxmax1_,igrid)<=xprobmax1-(xprobmax1-xprobmin1)&
         *writespshift(1,2).and.rnode(rpxmax2_,igrid)<=xprobmax2-(xprobmax2&
         -xprobmin2)*writespshift(2,2))) then
       d3grid=zgridsc*(rnode(rpxmax1_,igrid)-rnode(rpxmin1_,igrid))
       n3grid=dnint(zlength/d3grid)
       ! In case primitives to be saved: use primitive subroutine
       !  extra layer around mesh only needed when storing corner values and averaging
       if(saveprim) then
        call primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
         ixGlo1,ixGlo2,ixGhi1,ixGhi2,ps(igrid)%w,ps(igrid)%x)
       endif
       ! using array w so that new output auxiliaries can be calculated by the user
       ! extend 2D data to 3D insuring variables are independent on the third coordinate
       do ix3=ixGlo1,ixGhi1
        w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix3,1:nw)=ps(igrid)%w(ixGlo1:ixGhi1,&
          ixGlo2:ixGhi2,1:nw)
       end do
       do i3grid=1,n3grid !subcycles
        call calc_grid23(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
           ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,ixCCmin2,&
           ixCCmin3,ixCCmax1,ixCCmax2,ixCCmax3,.true.,i3grid,d3grid,w,zlength,zgridsc)
        do iw=1,nw
          if(.not.w_write(iw))cycle
          select case(userconvert_type)
            case('vtuB23')
              write(qunit) length
              write(qunit) (((real(wC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                 =ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2),ix3=ixCmin3,ixCmax3)
            case('vtuBCC23')
              write(qunit) lengthcc
              write(qunit) (((real(wCC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                 =ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2),ix3&
                 =ixCCmin3,ixCCmax3)
          end select
        enddo
        if(nwauxio>0)then
         do iw=nw+1,nw+nwauxio
           select case(userconvert_type)
             case('vtuB23')
               write(qunit) length
               write(qunit) (((real(wC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                  =ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2),ix3=ixCmin3,ixCmax3)
             case('vtuBCC23')
               write(qunit) lengthcc
               write(qunit) (((real(wCC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                  =ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2),ix3&
                  =ixCCmin3,ixCCmax3)
           end select
         enddo
        endif
        write(qunit) length_coords
        do ix3=ixCmin3,ixCmax3
        do ix2=ixCmin2,ixCmax2
        do ix1=ixCmin1,ixCmax1
          x_VTK(1:3)=zero;
          x_VTK(1:3)=xC_TMP(ix1,ix2,ix3,1:3)*normconv(0);
          do k=1,3
           write(qunit) reaL(x_VTK(k))
          end do
        end do
        end do
        end do

        write(qunit) length_conn
        do ix3=1,nx3
        do ix2=1,nx2
        do ix1=1,nx1



        write(qunit)&
        (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1-1, &
        (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
        (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1-1,&
        (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1,&
         ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1-1,&
         ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
         ix3*nxC2*nxC1+    ix2*nxC1+ix1-1,&
         ix3*nxC2*nxC1+    ix2*nxC1+ix1

        end do
        end do
        end do
        write(qunit) length_offsets
        do icel=1,nc
           write(qunit) icel*(2**3)
        end do




        VTK_type=11
        write(qunit) size_int*nc
        do icel=1,nc
         write(qunit) VTK_type
        enddo
       enddo !subcycles
      endif
   end do
 endif
end do

close(qunit)
open(qunit,file=filename,status='unknown',form='formatted',position='append')

write(qunit,'(a)')'</AppendedData>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)
end subroutine unstructuredvtkB23
!====================================================================================
\}
{^IFTWOD 
subroutine unstructuredvtkBsym23(qunit,userconvert_type)

! output for vtu format to paraview, binary version output
! not parallel, uses calc_grid to compute nwauxio variables
! use this subroutine  when the physical domain is symmetric/asymmetric about (0,y,z) 
! plane, xprobmin1=0 and the computational domain is a half of the physical domain

use mod_global_parameters

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,3) :: xC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,&
   3)   :: xCC_TMP

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,nw&
   +nwauxio)     :: wCC_TMP
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo1:ixGhi1,1:nw&
   +nwauxio)   :: w
integer::               igrid,iigrid,level,igonlevel,icel,ixCmin1,ixCmin2,&
   ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,ixCCmin2,ixCCmin3,ixCCmax1,&
   ixCCmax2,ixCCmax3
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx1,nx2,nx3,nxC1,nxC2,nxC3,nodesonlevel,elemsonlevel,nc,np,&
   VTK_type,ix1,ix2,ix3
double precision :: normconv(0:nw+nwauxio)
character(len=80)::  filename
character(len=20):: userconvert_type
character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:3+nw+nwauxio)
character(len=1024) :: outfilehead
integer*8 :: offset
integer :: size_int,size_double,size_length,recsep,k,iw
integer :: length,lengthcc,offset_points,offset_cells, length_coords,&
   length_conn,length_offsets
integer :: i3grid,n3grid
double precision ::d3grid,zlengsc,zgridsc
character::  buffer
character(len=6)::  bufform

double precision :: zlength

logical ::   fileopen
!-----------------------------------------------------------------------------
if(npe>1)then
 if(mype==0) PRINT *,'unstructuredvtkBsym23 not parallel, use vtumpi'
 call mpistop('npe>1, unstructuredvtkBsym23')
end if

offset=0
recsep=4
size_double=4
size_length=4
size_int=size_length

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
  ! generate filename 
  write(filename,'(a,a,i4.4,a)') TRIM(base_filename),"3D",snapshotini,".vtu"
  ! Open the file for the header part
  open(qunit,file=filename,status='unknown')
endif

call getheadernames(wnamei,xandwnamei,outfilehead)
! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'<UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(qunit,'(f10.2)') real(global_time*time_convert_factor)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'

! number of cells, number of corner points, per grid.
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi1-ixMlo1+1;
nxC1=nx1+1;nxC2=nx2+1;nxC3=nx3+1;
nc=nx1*nx2*nx3
np=nxC1*nxC2*nxC3

length=np*size_double
lengthcc=nc*size_double

length_coords=3*length
length_conn=2**3*size_int*nc
length_offsets=nc*size_int

! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
zlengsc=4.d0
zgridsc=2.d0
zlength=zlengsc*(xprobmax1-xprobmin1)
do level=levmin,levmax
 if (writelevel(level)) then
   do iigrid=1,igridstail; igrid=igrids(iigrid);
    if (node(plevel_,igrid)/=level) cycle
    block=>ps(igrid)
    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    if ((rnode(rpxmin1_,igrid)>=xprobmin1+(xprobmax1-xprobmin1)&
       *writespshift(1,1).and.rnode(rpxmin2_,igrid)>=xprobmin2&
       +(xprobmax2-xprobmin2)*writespshift(2,1))&
       .and.(rnode(rpxmax1_,igrid)<=xprobmax1-(xprobmax1-xprobmin1)&
       *writespshift(1,2).and.rnode(rpxmax2_,igrid)<=xprobmax2-(xprobmax2&
       -xprobmin2)*writespshift(2,2))) then
     d3grid=zgridsc*(rnode(rpxmax1_,igrid)-rnode(rpxmin1_,igrid))
     n3grid=dnint(zlength/d3grid)
     do i3grid=1,n3grid !subcycles
      !! original domain ----------------------------------start 
      select case(userconvert_type)
       case('vtuBsym23')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
            '" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<PointData>'
         do iw=1,nw
            if(.not.w_write(iw))cycle
            write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
               TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+length+size_length
         enddo
         if(nwauxio>0)then
          do iw=nw+1,nw+nwauxio
             write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
                TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+length+size_length
          enddo
         endif
         write(qunit,'(a)')'</PointData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
            offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_length
         write(qunit,'(a)')'</Points>'
       case('vtuBCCsym23')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
            '" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<CellData>'
         do iw=1,nw
            if(.not.w_write(iw))cycle
            write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
               TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+lengthcc+size_length
         enddo
         if(nwauxio>0)then
          do iw=nw+1,nw+nwauxio
             write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
                TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+lengthcc+size_length
          enddo
         endif
         write(qunit,'(a)')'</CellData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
            offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_length
         write(qunit,'(a)')'</Points>'
      end select


      write(qunit,'(a)')'<Cells>'

      ! connectivity part
      write(qunit,'(a,i16,a)')&
'<DataArray type="Int32" Name="connectivity" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_conn+size_length

      ! offsets data array
      write(qunit,'(a,i16,a)') &
         '<DataArray type="Int32" Name="offsets" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_offsets+size_length

      ! VTK cell type data array
      write(qunit,'(a,i16,a)') &
         '<DataArray type="Int32" Name="types" format="appended" offset="',&
         offset,'"/>'
      offset=offset+size_length+nc*size_int

      write(qunit,'(a)')'</Cells>'

      write(qunit,'(a)')'</Piece>'
      !! original domain ----------------------------------end 
      !! symetric/asymetric mirror domain -----------------start
      select case(userconvert_type)
       case('vtuBsym23')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
            '" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<PointData>'
         do iw=1,nw
            if(.not.w_write(iw))cycle
            write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
               TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+length+size_length
         enddo
         if(nwauxio>0)then
          do iw=nw+1,nw+nwauxio
             write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
                TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+length+size_length
          enddo
         endif
         write(qunit,'(a)')'</PointData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
            offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_length
         write(qunit,'(a)')'</Points>'
       case('vtuBCCsym23')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
            '" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<CellData>'
         do iw=1,nw
            if(.not.w_write(iw))cycle
            write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
               TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+lengthcc+size_length
         enddo
         if(nwauxio>0)then
          do iw=nw+1,nw+nwauxio
             write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
                TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+lengthcc+size_length
          enddo
         endif
         write(qunit,'(a)')'</CellData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
            offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_length
         write(qunit,'(a)')'</Points>'
      end select


      write(qunit,'(a)')'<Cells>'

      ! connectivity part
      write(qunit,'(a,i16,a)')&
'<DataArray type="Int32" Name="connectivity" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_conn+size_length

      ! offsets data array
      write(qunit,'(a,i16,a)') &
         '<DataArray type="Int32" Name="offsets" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_offsets+size_length

      ! VTK cell type data array
      write(qunit,'(a,i16,a)') &
         '<DataArray type="Int32" Name="types" format="appended" offset="',&
         offset,'"/>'
      offset=offset+size_length+nc*size_int

      write(qunit,'(a)')'</Cells>'

      write(qunit,'(a)')'</Piece>'
      !! symetric/asymetric mirror domain -----------------end
     enddo !subcycles
    endif
   enddo
 endif
enddo

write(qunit,'(a)')'</UnstructuredGrid>'
write(qunit,'(a)')'<AppendedData encoding="raw">'

close(qunit)
open(qunit,file=filename,form='unformatted',access='stream',status='old',position='append')
buffer='_'
write(qunit) TRIM(buffer)
do level=levmin,levmax
 if (writelevel(level)) then
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (node(plevel_,igrid)/=level) cycle
      block=>ps(igrid)
      ! only output a grid when fully within clipped region selected
      ! by writespshift array
      if ((rnode(rpxmin1_,igrid)>=xprobmin1+(xprobmax1-xprobmin1)&
         *writespshift(1,1).and.rnode(rpxmin2_,igrid)>=xprobmin2&
         +(xprobmax2-xprobmin2)*writespshift(2,1))&
         .and.(rnode(rpxmax1_,igrid)<=xprobmax1-(xprobmax1-xprobmin1)&
         *writespshift(1,2).and.rnode(rpxmax2_,igrid)<=xprobmax2-(xprobmax2&
         -xprobmin2)*writespshift(2,2))) then
       d3grid=zgridsc*(rnode(rpxmax1_,igrid)-rnode(rpxmin1_,igrid))
       n3grid=dnint(zlength/d3grid)
       ! In case primitives to be saved: use primitive subroutine
       !  extra layer around mesh only needed when storing corner values and averaging
       if(saveprim) then
        call primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
         ixGlo1,ixGlo2,ixGhi1,ixGhi2,ps(igrid)%w,ps(igrid)%x)
       endif
       ! using array w so that new output auxiliaries can be calculated by the user
       ! extend 2D data to 3D insuring variables are independent on the third coordinate
       do ix3=ixGlo1,ixGhi1
        w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix3,1:nw)=ps(igrid)%w(ixGlo1:ixGhi1,&
          ixGlo2:ixGhi2,1:nw)
       end do
       do i3grid=1,n3grid !subcycles
        call calc_grid23(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
           ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,ixCCmin2,&
           ixCCmin3,ixCCmax1,ixCCmax2,ixCCmax3,.true.,i3grid,d3grid,w,zlength,zgridsc)
        !! original domain ----------------------------------start 
        do iw=1,nw
          if(.not.w_write(iw))cycle
          select case(userconvert_type)
            case('vtuBsym23')
              write(qunit) length
              write(qunit) (((real(wC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                 =ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2),ix3=ixCmin3,ixCmax3)
            case('vtuBCCsym23')
              write(qunit) lengthcc
              write(qunit) (((real(wCC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                 =ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2),ix3&
                 =ixCCmin3,ixCCmax3)
          end select
        enddo
        if(nwauxio>0)then
         do iw=nw+1,nw+nwauxio
           select case(userconvert_type)
             case('vtuBsym23')
               write(qunit) length
               write(qunit) (((real(wC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                  =ixCmin1,ixCmax1),ix2=ixCmin2,ixCmax2),ix3=ixCmin3,ixCmax3)
             case('vtuBCCsym23')
               write(qunit) lengthcc
               write(qunit) (((real(wCC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                  =ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2),ix3&
                  =ixCCmin3,ixCCmax3)
           end select
         enddo
        endif
        write(qunit) length_coords
        do ix3=ixCmin3,ixCmax3
        do ix2=ixCmin2,ixCmax2
        do ix1=ixCmin1,ixCmax1
          x_VTK(1:3)=zero;
          x_VTK(1:3)=xC_TMP(ix1,ix2,ix3,1:3)*normconv(0);
          do k=1,3
           write(qunit) real(x_VTK(k))
          end do
        end do
        end do
        end do

        write(qunit) length_conn
        do ix3=1,nx3
        do ix2=1,nx2
        do ix1=1,nx1



        write(qunit)&
        (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1-1, &
        (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
        (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1-1,&
        (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1,&
         ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1-1,&
         ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
         ix3*nxC2*nxC1+    ix2*nxC1+ix1-1,&
         ix3*nxC2*nxC1+    ix2*nxC1+ix1

        end do
        end do
        end do
        write(qunit) length_offsets
        do icel=1,nc
           write(qunit) icel*(2**3)
        end do




        VTK_type=11
        write(qunit) size_int*nc
        do icel=1,nc
         write(qunit) VTK_type
        enddo
        !! original domain ----------------------------------end 
        !! symetric/asymetric mirror domain -----------------start
        do iw=1,nw
          if(.not.w_write(iw))cycle
          if(iw==2 .or. iw==4 .or. iw==7) then
           wC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)=&
           -wC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
           wCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,iw)=&
           -wCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,iw)
          endif
          select case(userconvert_type)
            case('vtuBsym23')
              write(qunit) length
              write(qunit) (((real(wC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                 =ixCmax1,ixCmin1,-1),ix2=ixCmin2,ixCmax2),ix3=ixCmin3,ixCmax3)
            case('vtuBCCsym23')
              write(qunit) lengthcc
              write(qunit) (((real(wCC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                 =ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2),ix3&
                 =ixCCmin3,ixCCmax3)
          end select
        enddo
        if(nwauxio>0)then
         do iw=nw+1,nw+nwauxio
           select case(userconvert_type)
             case('vtuBsym23')
               write(qunit) length
               write(qunit) (((real(wC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                  =ixCmax1,ixCmin1,-1),ix2=ixCmin2,ixCmax2),ix3=ixCmin3,ixCmax3)
             case('vtuBCCsym23')
               write(qunit) lengthcc
               write(qunit) (((real(wCC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                  =ixCCmin1,ixCCmax1),ix2=ixCCmin2,ixCCmax2),ix3&
                  =ixCCmin3,ixCCmax3)
           end select
         enddo
        endif
        write(qunit) length_coords
        do ix3=ixCmin3,ixCmax3
        do ix2=ixCmin2,ixCmax2
        do ix1=ixCmax1,ixCmin1,-1
          x_VTK(1:3)=zero;
          x_VTK(1:3)=xC_TMP(ix1,ix2,ix3,1:3)*normconv(0);
          x_VTK(1)=-x_VTK(1)
          do k=1,3
           write(qunit) real(x_VTK(k))
          end do
        end do
        end do
        end do

        write(qunit) length_conn
        do ix3=1,nx3
        do ix2=1,nx2
        do ix1=nx1,1,-1



        write(qunit)&
        (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1-1, &
        (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
        (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1-1,&
        (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1,&
         ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1-1,&
         ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
         ix3*nxC2*nxC1+    ix2*nxC1+ix1-1,&
         ix3*nxC2*nxC1+    ix2*nxC1+ix1

        end do
        end do
        end do
        write(qunit) length_offsets
        do icel=1,nc
           write(qunit) icel*(2**3)
        end do




        VTK_type=11
        write(qunit) size_int*nc
        do icel=1,nc
         write(qunit) VTK_type
        enddo
        !! symetric/asymetric mirror domain -----------------end
       enddo !subcycles
      endif
   end do
 endif
end do

close(qunit)
open(qunit,file=filename,status='unknown',form='formatted',position='append')

write(qunit,'(a)')'</AppendedData>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)
end subroutine unstructuredvtkBsym23
!=============================================================================
subroutine calc_grid23(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,ixCCmin2,ixCCmin3,&
 ixCCmax1,ixCCmax2,ixCCmax3,first,i3grid,d3grid,w,zlength,zgridsc)

! this subroutine computes both corner as well as cell-centered values
! it handles how we do the center to corner averaging, as well as 
! whether we switch to cartesian or want primitive or conservative output,
! handling the addition of B0 in B0+B1 cases, ...
!
! the normconv is passed on to specialvar_output for extending with
! possible normalization values for the nw+1:nw+nwauxio entries

use mod_global_parameters
integer, intent(in) :: qunit, igrid,i3grid
logical, intent(in) :: first

integer :: nx1,nx2,nx3, nxC1,nxC2,nxC3, ix1,ix2,ix3, ix, iw, level, idir
integer :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,&
   ixCCmin2,ixCCmin3,ixCCmax1,ixCCmax2,ixCCmax3,nxCC1,nxCC2,nxCC3
double precision :: dx1,dx2,dx3,d3grid,zlength,zgridsc

integer :: idims,jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3
double precision :: ldw(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo1:ixGhi1),&
    dwC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo1:ixGhi1)

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,3) :: xC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,&
   3)   :: xCC

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,nw+nwauxio)   :: wC
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,nw&
   +nwauxio)     :: wCC

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,3) :: xC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,&
   3)   :: xCC_TMP

double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,ixMlo1&
   -1:ixMhi1,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo1:ixMhi1,nw&
   +nwauxio)     :: wCC_TMP
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo1:ixGhi1,1:nw&
   +nwauxio)   :: w

double precision,dimension(0:nw+nwauxio)       :: normconv
integer ::iwe,iwb1,iwb2,iwb3
logical, save :: subfirst=.true.
!-----------------------------------------------------------------------------
! following only for allowing compiler to go through with debug on
{#IFDEF ENERGY
iwe=e_
}
iwb1=mag(1)
iwb2=mag(2)
iwb3=mag(2)+1

nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi1-ixMlo1+1;
level=node(plevel_,igrid)
dx1=dx(1,level);dx2=dx(2,level);dx3=zgridsc*dx(1,level);

! for normalization within the code
if(saveprim) then
  normconv(0) = length_convert_factor
  normconv(1:nw) = w_convert_factor
else
  normconv(0)=length_convert_factor
  ! assuming density
  normconv(1)=w_convert_factor(1)
  ! assuming momentum=density*velocity
  if (nw>=2) normconv(2:2+3)=w_convert_factor(1)*w_convert_factor(2:2+3)
  ! assuming energy/pressure and magnetic field
  if (nw>=2+3) normconv(2+3:nw)=w_convert_factor(2+3:nw)
end if

! coordinates of cell centers
nxCC1=nx1;nxCC2=nx2;nxCC3=nx3;
ixCCmin1=ixMlo1;ixCCmin2=ixMlo2;ixCCmin3=ixMlo1; ixCCmax1=ixMhi1
ixCCmax2=ixMhi2;ixCCmax3=ixMhi1;
do ix=ixCCmin1,ixCCmax1
    xCC(ix,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,1)=rnode(rpxmin1_,igrid)&
       +(dble(ix-ixCCmin1)+half)*dx1
end do
do ix=ixCCmin2,ixCCmax2
    xCC(ixCCmin1:ixCCmax1,ix,ixCCmin3:ixCCmax3,2)=rnode(rpxmin2_,igrid)&
       +(dble(ix-ixCCmin2)+half)*dx2
end do
do ix=ixCCmin3,ixCCmax3
    xCC(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ix,3)=-zlength/two+&
       dble(i3grid-1)*d3grid+(dble(ix-ixCCmin3)+half)*dx3
end do

! coordinates of cell corners
nxC1=nx1+1;nxC2=nx2+1;nxC3=nx3+1;
ixCmin1=ixMlo1-1;ixCmin2=ixMlo2-1;ixCmin3=ixMlo1-1; ixCmax1=ixMhi1
ixCmax2=ixMhi2;ixCmax3=ixMhi1;
do ix=ixCmin1,ixCmax1
    xC(ix,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=rnode(rpxmin1_,igrid)&
       +dble(ix-ixCmin1)*dx1
end do
do ix=ixCmin2,ixCmax2
    xC(ixCmin1:ixCmax1,ix,ixCmin3:ixCmax3,2)=rnode(rpxmin2_,igrid)&
       +dble(ix-ixCmin2)*dx2
end do
do ix=ixCmin3,ixCmax3
    xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ix,3)=-zlength/two+&
      dble(i3grid-1)*d3grid+dble(ix-ixCmin3)*dx3
end do


if (nwextra>0) then
 ! here we actually fill the ghost layers for the nwextra variables using 
 ! continuous extrapolation (as these values do not exist normally in ghost
 ! cells)
 do idims=1,3
  select case(idims)
   case(1)
     jxCmin1=ixGhi1+1-nghostcells;jxCmin2=ixGlo2;jxCmin3=ixGlo1;
     jxCmax1=ixGhi1;jxCmax2=ixGhi2;jxCmax3=ixGhi1;
     do ix1=jxCmin1,jxCmax1
         w(ix1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,nw-nwextra+1:nw) = w(jxCmin1&
            -1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,nw-nwextra+1:nw)
     end do
     jxCmin1=ixGlo1;jxCmin2=ixGlo2;jxCmin3=ixGlo1;
     jxCmax1=ixGlo1-1+nghostcells;jxCmax2=ixGhi2;jxCmax3=ixGhi1;
     do ix1=jxCmin1,jxCmax1
         w(ix1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,nw-nwextra+1:nw) = w(jxCmax1&
            +1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,nw-nwextra+1:nw)
     end do
   case(2)
     jxCmin1=ixGlo1;jxCmin2=ixGhi2+1-nghostcells;jxCmin3=ixGlo1;
     jxCmax1=ixGhi1;jxCmax2=ixGhi2;jxCmax3=ixGhi1;
     do ix2=jxCmin2,jxCmax2
         w(jxCmin1:jxCmax1,ix2,jxCmin3:jxCmax3,nw-nwextra+1:nw) &
            = w(jxCmin1:jxCmax1,jxCmin2-1,jxCmin3:jxCmax3,nw-nwextra+1:nw)
     end do
     jxCmin1=ixGlo1;jxCmin2=ixGlo2;jxCmin3=ixGlo1;
     jxCmax1=ixGhi1;jxCmax2=ixGlo2-1+nghostcells;jxCmax3=ixGhi1;
     do ix2=jxCmin2,jxCmax2
         w(jxCmin1:jxCmax1,ix2,jxCmin3:jxCmax3,nw-nwextra+1:nw) &
            = w(jxCmin1:jxCmax1,jxCmax2+1,jxCmin3:jxCmax3,nw-nwextra+1:nw)
     end do
   case(3)
     jxCmin1=ixGlo1;jxCmin2=ixGlo2;jxCmin3=ixGhi1+1-nghostcells;
     jxCmax1=ixGhi1;jxCmax2=ixGhi2;jxCmax3=ixGhi1;
     do ix3=jxCmin3,jxCmax3
         w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,ix3,nw-nwextra+1:nw) &
            = w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3-1,nw-nwextra+1:nw)
     end do
     jxCmin1=ixGlo1;jxCmin2=ixGlo2;jxCmin3=ixGlo1;
     jxCmax1=ixGhi1;jxCmax2=ixGhi2;jxCmax3=ixGlo1-1+nghostcells;
     do ix3=jxCmin3,jxCmax3
         w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,ix3,nw-nwextra+1:nw) &
            = w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmax3+1,nw-nwextra+1:nw)
     end do
  end select
 end do
end if
! next lines needed when specialvar_output uses gradients
! and later on when dwlimiter2 is used 
if(nwauxio>0)then
  ! auxiliary io variables can be computed and added by user
  ! next few lines ensure correct usage of routines like divvector etc
  dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
  if (B0field) then
    myB0_cell => pB0_cell(igrid)
    myB0      => pB0_cell(igrid)
    myB0_face1 => pB0_face1(igrid)
    myB0_face2 => pB0_face2(igrid)
  !  myB0_face3 => pB0_face3(igrid)
  end if
  ! default (no) normalization for auxiliary variables
  normconv(nw+1:nw+nwauxio)=one
  ! maybe need for restriction to ixG^LL^LSUB1 
  call specialvar_output23(ixGlo1,ixGlo2,ixGlo1,ixGhi1,ixGhi2,ixGhi1,ixGlo1&
     +1,ixGlo2+1,ixGlo1+1,ixGhi1-1,ixGhi2-1,ixGhi1-1,w,xCC,normconv)
endif

! compute the cell-center values for w first
!===========================================
! cell center values obtained from mere copy, while B0+B1 split handled here
do iw=1,nw+nwauxio
   if (B0field.and.iw>b0_.and.iw<=b0_+ndir) then
      idir=iw-b0_
      do ix3=ixCCmin3,ixCCmax3
      do ix2=ixCCmin2,ixCCmax2
      do ix1=ixCCmin1,ixCCmax1
         wCC(ix1,ix2,ix3,iw)=w(ix1,ix2,ix3,iw)+pB0_cell(igrid)%w(ix1,ix2,&
            idir)
      end do
      end do
      end do
   else
      do ix3=ixCCmin3,ixCCmax3
      do ix2=ixCCmin2,ixCCmax2
      do ix1=ixCCmin1,ixCCmax1
          wCC(ix1,ix2,ix3,iw)=w(ix1,ix2,ix3,iw)
      end do
      end do
      end do
   end if
end do
{#IFDEF ENERGY
if((.not.saveprim) .and. B0field) then
   do ix3=ixCCmin3,ixCCmax3
   do ix2=ixCCmin2,ixCCmax2
   do ix1=ixCCmin1,ixCCmax1
       wCC(ix1,ix2,ix3,iwe)=w(ix1,ix2,ix3,iwe) +half*( pB0_cell(igrid)%w(ix1,&
          ix2,iwb1-b0_)**2+pB0_cell(igrid)%w(ix1,ix2,iwb2-b0_)**2&
          +pB0_cell(igrid)%w(ix1,ix2,iwb3-b0_)**2 ) + ( w(ix1,ix2,ix3,&
          iwb1)*pB0_cell(igrid)%w(ix1,ix2,iwb1-b0_)+w(ix1,ix2,ix3,iwb2)&
          *pB0_cell(igrid)%w(ix1,ix2,iwb2-b0_)+w(ix1,ix2,ix3,iwb3)&
          *pB0_cell(igrid)%w(ix1,ix2,iwb3-b0_) )
   end do
   end do
   end do
endif
}
! compute the corner values for w now by averaging
!=================================================
if(slab_uniform)then
   ! for slab symmetry: no geometrical info required
   do iw=1,nw+nwauxio
      if (B0field.and.iw>b0_.and.iw<=b0_+ndir) then
         idir=iw-b0_
         do ix3=ixCmin3,ixCmax3
         do ix2=ixCmin2,ixCmax2
         do ix1=ixCmin1,ixCmax1
           wC(ix1,ix2,ix3,iw)=sum(w(ix1:ix1+1,ix2:ix2+1,ix3,iw) &
              +pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1&
              ,idir))/dble(2**3)+&
              sum(w(ix1:ix1+1,ix2:ix2+1,ix3+1,iw) &
              +pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1&
              ,idir))/dble(2**3)
         end do
         end do
         end do
      else
         do ix3=ixCmin3,ixCmax3
         do ix2=ixCmin2,ixCmax2
         do ix1=ixCmin1,ixCmax1
            wC(ix1,ix2,ix3,iw)=sum(w(ix1:ix1+1,ix2:ix2+1,ix3:ix3&
               +1,iw))/dble(2**3)
         end do
         end do
         end do
      end if
   end do
{#IFDEF ENERGY
   if((.not.saveprim) .and. B0field) then
      do ix3=ixCmin3,ixCmax3
      do ix2=ixCmin2,ixCmax2
      do ix1=ixCmin1,ixCmax1
         wC(ix1,ix2,ix3,iwe)=sum(w(ix1:ix1+1,ix2:ix2+1,ix3,iwe) &
            +half*( pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1&
            ,iwb1-b0_)**2+pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1&
            ,iwb2-b0_)**2+pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1&
            ,iwb3-b0_)**2 ) + ( w(ix1:ix1+1,ix2:ix2+1,ix3&
            ,iwb1)*pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1&
            ,iwb1-b0_)+w(ix1:ix1+1,ix2:ix2+1,ix3,iwb2)&
            *pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1,iwb2-b0_)&
            +w(ix1:ix1+1,ix2:ix2+1,ix3,iwb3)*pB0_cell(igrid)%w(ix1:ix1&
            +1,ix2:ix2+1,iwb3-b0_) ) ) /dble(2**3)+&
            sum(w(ix1:ix1+1,ix2:ix2+1,ix3+1,iwe) &
            +half*( pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1&
            ,iwb1-b0_)**2+pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1&
            ,iwb2-b0_)**2+pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1&
            ,iwb3-b0_)**2 ) + ( w(ix1:ix1+1,ix2:ix2+1,ix3&
            +1,iwb1)*pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1&
            ,iwb1-b0_)+w(ix1:ix1+1,ix2:ix2+1,ix3+1,iwb2)&
            *pB0_cell(igrid)%w(ix1:ix1+1,ix2:ix2+1,iwb2-b0_)&
            +w(ix1:ix1+1,ix2:ix2+1,ix3+1,iwb3)*pB0_cell(igrid)%w(ix1:ix1&
            +1,ix2:ix2+1,iwb3-b0_) ) ) /dble(2**3)
      end do
      end do
      end do
   endif
}
endif
! keep the coordinate and vector components
xC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:3)          &
   = xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:3)
wC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:nw&
   +nwauxio)    = wC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1:nw&
   +nwauxio)
xCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,&
   1:3)        = xCC(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,&
   ixCCmin3:ixCCmax3,1:3)
wCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,1:nw&
   +nwauxio)  = wCC(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,&
   1:nw+nwauxio)
end subroutine calc_grid23
!=============================================================================
subroutine save_connvtk23(qunit,igrid)

! this saves the basic line, pixel and voxel connectivity,
! as used by VTK file outputs for unstructured grid

use mod_global_parameters

integer, intent(in) :: qunit, igrid

integer :: nx1,nx2,nx3, nxC1,nxC2,nxC3, ix1,ix2,ix3
!-----------------------------------------------------------------------------
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi1-ixMlo1+1;
nxC1=nx1+1;nxC2=nx2+1;nxC3=nx3+1;

do ix3=1,nx3
do ix2=1,nx2
do ix1=1,nx1



        write(qunit,'(8(i7,1x))')&
                   (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1-1, &
                   (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
                   (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1-1,&
                   (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1,&
                       ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1-1,&
                       ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
                       ix3*nxC2*nxC1+    ix2*nxC1+ix1-1,&
                       ix3*nxC2*nxC1+    ix2*nxC1+ix1

end do
end do
end do
end subroutine save_connvtk23
!=============================================================================
subroutine specialvar_output23(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
 ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
! corresponding normalization values (default value 1)

use mod_global_parameters

integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:3)
double precision              :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,nw+nwauxio)
double precision              :: normconv(0:nw+nwauxio)

double precision :: qvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo1:ixGhi1,1:ndir),&
   curlvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo1:ixGhi1,1:ndir)
integer                       :: idirmin
!-----------------------------------------------------------------------------
! output Te
if(saveprim)then
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3,rho_)
endif
!!! store current
! qvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1)=w(ixImin1:ixImax1,&
!    ixImin2:ixImax2,ixImin3:ixImax3,mag(1))
! qvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2)=w(ixImin1:ixImax1,&
!    ixImin2:ixImax2,ixImin3:ixImax3,mag(2))
! qvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3)=w(ixImin1:ixImax1,&
!    ixImin2:ixImax2,ixImin3:ixImax3,mag(3));
!call curlvector3D(qvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
!   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,curlvec,idirmin,1,ndir)
!w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+2)=curlvec&
!   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
!w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+3)=curlvec&
!   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
!w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+4)=curlvec&
!   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3);
end subroutine specialvar_output23 \}
!====================================================================================
subroutine unstructuredvtkBsym(qunit,userconvert_type)

! output for vtu format to paraview, binary version output
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using w_convert_factor-array

use mod_global_parameters

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP

integer::               igrid,iigrid,level,igonlevel,icel,ixC^L,ixCC^L
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx^D,nxC^D,nodesonlevel,elemsonlevel,nc,np,VTK_type,ix^D
double precision :: normconv(0:nw+nwauxio)
character(len=80)::  filename
character(len=20):: userconvert_type
character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

integer*8 :: offset
integer ::  size_int,size_double,size_length,recsep,k,iw
integer ::  length,lengthcc,offset_points,offset_cells, &
           length_coords,length_conn,length_offsets
character::  buffer
character(len=6)::  bufform

logical ::   fileopen
!-----------------------------------------------------------------------------

if(npe>1)then
 if(mype==0) PRINT *,'unstructuredvtkBsym not parallel, use vtumpi'
 call mpistop('npe>1, unstructuredvtkBsym')
end if

offset=0
recsep=4
size_double=4
size_length=4
size_int=size_length

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
  ! generate filename 
  write(filename,'(a,i4.4,a)') TRIM(base_filename),snapshotini,".vtu"
  ! Open the file for the header part
  open(qunit,file=filename,status='unknown')
endif

call getheadernames(wnamei,xandwnamei,outfilehead)

! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'<UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(qunit,'(f10.2)') real(global_time*time_convert_factor)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'

! number of cells, number of corner points, per grid.
nx^D=ixMhi^D-ixMlo^D+1;
nxC^D=nx^D+1;
nc={nx^D*}
np={nxC^D*}

length=np*size_double
lengthcc=nc*size_double

length_coords=3*length
length_conn=2**^ND*size_int*nc
length_offsets=nc*size_int

! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
do level=levmin,levmax
 if (writelevel(level)) then
   do iigrid=1,igridstail; igrid=igrids(iigrid);
    if (node(plevel_,igrid)/=level) cycle
    block=>ps(igrid)
    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    if (({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
          *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
         <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})) then
      select case(userconvert_type)
       case('vtuBsym')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') &
            '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<PointData>'
         do iw=1,nw
            if(.not.w_write(iw))cycle
            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+length+size_length
         enddo
         if(nwauxio>0) then
          do iw=nw+1,nw+nwauxio
             write(qunit,'(a,a,a,i16,a)')&
                 '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                 '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+length+size_length
          enddo
         endif
         write(qunit,'(a)')'</PointData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
     '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_length
         write(qunit,'(a)')'</Points>'
       case('vtuBCCsym')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') &
            '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<CellData>'
         do iw=1,nw
            if(.not.w_write(iw))cycle
            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+lengthcc+size_length
         enddo
         if(nwauxio>0) then
          do iw=nw+1,nw+nwauxio
             write(qunit,'(a,a,a,i16,a)')&
                 '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                 '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+lengthcc+size_length
          enddo
         endif
         write(qunit,'(a)')'</CellData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
     '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_length
         write(qunit,'(a)')'</Points>'
      end select

   
      write(qunit,'(a)')'<Cells>'

      ! connectivity part
      write(qunit,'(a,i16,a)')&
        '<DataArray type="Int32" Name="connectivity" format="appended" offset="',offset,'"/>'
      offset=offset+length_conn+size_length    

      ! offsets data array
      write(qunit,'(a,i16,a)') &
        '<DataArray type="Int32" Name="offsets" format="appended" offset="',offset,'"/>'
      offset=offset+length_offsets+size_length    

      ! VTK cell type data array
      write(qunit,'(a,i16,a)') &
        '<DataArray type="Int32" Name="types" format="appended" offset="',offset,'"/>' 
      offset=offset+size_length+nc*size_int

      write(qunit,'(a)')'</Cells>'

      write(qunit,'(a)')'</Piece>'
      !! symmetry axis/plane(x=0) mirrored part --------------------------------start
      if(userconvert_type=='vtuBsym' .or. userconvert_type=='vtuBCCsym') then
       select case(userconvert_type)
        case('vtuBsym')
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') &
             '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<PointData>'
          do iw=1,nw
             if(.not.w_write(iw))cycle
             write(qunit,'(a,a,a,i16,a)')&
                 '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                 '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+length+size_length
          enddo
          if(nwauxio>0) then
           do iw=nw+1,nw+nwauxio
               write(qunit,'(a,a,a,i16,a)')&
                  '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                  '" format="appended" offset="',offset,'">'
              write(qunit,'(a)')'</DataArray>'
              offset=offset+length+size_length
           enddo
          endif
          write(qunit,'(a)')'</PointData>'

          write(qunit,'(a)')'<Points>'
          write(qunit,'(a,i16,a)') &
      '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
          ! write cell corner coordinates in a backward dimensional loop, always 3D output
          offset=offset+length_coords+size_length
          write(qunit,'(a)')'</Points>'
        case('vtuBCCsym')
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') &
             '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<CellData>'
          do iw=1,nw
             if(.not.w_write(iw))cycle
             write(qunit,'(a,a,a,i16,a)')&
                 '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                 '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+lengthcc+size_length
          enddo
          if(nwauxio>0) then
           do iw=nw+1,nw+nwauxio
              write(qunit,'(a,a,a,i16,a)')&
                  '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                  '" format="appended" offset="',offset,'">'
              write(qunit,'(a)')'</DataArray>'
              offset=offset+lengthcc+size_length
           enddo
          endif
          write(qunit,'(a)')'</CellData>'

          write(qunit,'(a)')'<Points>'
          write(qunit,'(a,i16,a)') &
      '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
          ! write cell corner coordinates in a backward dimensional loop, always 3D output
          offset=offset+length_coords+size_length
          write(qunit,'(a)')'</Points>'
       end select

   
       write(qunit,'(a)')'<Cells>'

       ! connectivity part
       write(qunit,'(a,i16,a)')&
         '<DataArray type="Int32" Name="connectivity" format="appended" offset="',offset,'"/>'
       offset=offset+length_conn+size_length    

       ! offsets data array
       write(qunit,'(a,i16,a)') &
         '<DataArray type="Int32" Name="offsets" format="appended" offset="',offset,'"/>'
       offset=offset+length_offsets+size_length    

       ! VTK cell type data array
       write(qunit,'(a,i16,a)') &
         '<DataArray type="Int32" Name="types" format="appended" offset="',offset,'"/>' 
       offset=offset+size_length+nc*size_int

       write(qunit,'(a)')'</Cells>'

       write(qunit,'(a)')'</Piece>'
      endif
      !! symmetry axis/plane x=0 mirrored part ----------------------------------end 
    endif
   enddo
 endif
enddo

write(qunit,'(a)')'</UnstructuredGrid>'
write(qunit,'(a)')'<AppendedData encoding="raw">'

close(qunit)
open(qunit,file=filename,form='unformatted',access='stream',status='old',position='append')
buffer='_'
write(qunit) TRIM(buffer)

do level=levmin,levmax
 if (writelevel(level)) then
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (node(plevel_,igrid)/=level) cycle
      block=>ps(igrid)
      ! only output a grid when fully within clipped region selected
      ! by writespshift array
      if (({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
          *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
         <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})) then
        call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                       ixC^L,ixCC^L,.true.)
        do iw=1,nw
          if(.not.w_write(iw))cycle
          select case(userconvert_type)
            case('vtuBsym')
              write(qunit) length
              write(qunit) {(|}real(wC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCmin^D,ixCmax^D)}
            case('vtuBCCsym')
              write(qunit) lengthcc
              write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
          end select 
        enddo
        if(nwauxio>0) then
         do iw=nw+1,nw+nwauxio
           select case(userconvert_type)
             case('vtuBsym')
               write(qunit) length
               write(qunit) {(|}real(wC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCmin^D,ixCmax^D)}
             case('vtuBCCsym')
               write(qunit) lengthcc
               write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
           end select 
         enddo
        endif

        write(qunit) length_coords
        {do ix^DB=ixCmin^DB,ixCmax^DB \}
          x_VTK(1:3)=zero;
          x_VTK(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0);
          do k=1,3
           write(qunit) real(x_VTK(k))
          end do
        {end do \}

        write(qunit) length_conn
        {do ix^DB=1,nx^DB\}
        {^IFONED write(qunit)ix1-1,ix1 \}
        {^IFTWOD
        write(qunit)(ix2-1)*nxC1+ix1-1, &
        (ix2-1)*nxC1+ix1,ix2*nxC1+ix1-1,ix2*nxC1+ix1
         \}
        {^IFTHREED
        write(qunit)&
        (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1-1, &
        (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
        (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1-1,&
        (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1,&
         ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1-1,&
         ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
         ix3*nxC2*nxC1+    ix2*nxC1+ix1-1,&
         ix3*nxC2*nxC1+    ix2*nxC1+ix1
         \}
        {end do\}

        write(qunit) length_offsets
        do icel=1,nc
           write(qunit) icel*(2**^ND)
        end do


       {^IFONED VTK_type=3 \}
       {^IFTWOD VTK_type=8 \}
       {^IFTHREED VTK_type=11 \}
        write(qunit) size_int*nc
        do icel=1,nc
         write(qunit) VTK_type
        enddo
        !! symmetry axis/plane x=0 mirrored part ----------------------------start
       if(userconvert_type=='vtuBsym' .or. userconvert_type=='vtuBCCsym') then
        do iw=1,nw
          if(.not.w_write(iw))cycle
          if(ndir==1) then
           if(iw==2) then
            wC_TMP(ixC^S,iw)=-wC_TMP(ixC^S,iw)
            wCC_TMP(ixCC^S,iw)=-wCC_TMP(ixCC^S,iw)
           endif
          else if(ndir==2) then
           if(iw==2 .or. iw==6) then
            wC_TMP(ixC^S,iw)=-wC_TMP(ixC^S,iw)
            wCC_TMP(ixCC^S,iw)=-wCC_TMP(ixCC^S,iw)
           endif
          else if(ndir==3) then
           if(iw==2 .or. iw==4 .or. iw==7) then
            wC_TMP(ixC^S,iw)=-wC_TMP(ixC^S,iw)
            wCC_TMP(ixCC^S,iw)=-wCC_TMP(ixCC^S,iw)
           endif
          endif
          select case(userconvert_type)
            case('vtuBsym')
              write(qunit) length
     {^IFONED write(qunit) (real(wC_TMP(ix1,iw)*normconv(iw)),ix1=ixCmax1,ixCmin1,-1)\}
     {^IFTWOD write(qunit) ((real(wC_TMP(ix1,ix2,iw)*normconv(iw)),ix1&
                =ixCmax1,ixCmin1,-1),ix2=ixCmin2,ixCmax2)\}
   {^IFTHREED write(qunit) (((real(wC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                =ixCmax1,ixCmin1,-1),ix2=ixCmin2,ixCmax2),ix3=ixCmin3,ixCmax3)\}
            case('vtuBCCsym')
              write(qunit) lengthcc
              write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
          end select 
        enddo
        if(nwauxio>0) then
         do iw=nw+1,nw+nwauxio
           select case(userconvert_type)
             case('vtuBsym')
               write(qunit) length
      {^IFONED write(qunit) (real(wC_TMP(ix1,iw)*normconv(iw)),ix1=ixCmax1,ixCmin1,-1)\}
      {^IFTWOD write(qunit) ((real(wC_TMP(ix1,ix2,iw)*normconv(iw)),ix1&
                 =ixCmax1,ixCmin1,-1),ix2=ixCmin2,ixCmax2)\}
    {^IFTHREED write(qunit) (((real(wC_TMP(ix1,ix2,ix3,iw)*normconv(iw)),ix1&
                 =ixCmax1,ixCmin1,-1),ix2=ixCmin2,ixCmax2),ix3=ixCmin3,ixCmax3)\}
             case('vtuBCCsym')
               write(qunit) lengthcc
               write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
           end select 
         enddo
        endif

        write(qunit) length_coords
        {^IFONED 
        do ix1=ixCmax1,ixCmin1,-1
          x_VTK(1:3)=zero;
          x_VTK(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0);
          x_VTK(1)=-x_VTK(1)
          do k=1,3
           write(qunit) real(x_VTK(k))
          end do
        end do
        write(qunit) length_conn
        do ix1=nx1,1,-1
         write(qunit)ix1-1,ix1
        end do
        \}
        {^IFTWOD
        do ix2=ixCmin2,ixCmax2
        do ix1=ixCmax1,ixCmin1,-1
          x_VTK(1:3)=zero;
          x_VTK(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0);
          x_VTK(1)=-x_VTK(1)
          do k=1,3
           write(qunit) real(x_VTK(k))
          end do
        end do
        end do
        write(qunit) length_conn
        do ix2=1,nx2
        do ix1=nx1,1,-1
         write(qunit)(ix2-1)*nxC1+ix1-1, &
         (ix2-1)*nxC1+ix1,ix2*nxC1+ix1-1,ix2*nxC1+ix1
        end do
        end do
        \}
        {^IFTHREED
        do ix3=ixCmin3,ixCmax3
        do ix2=ixCmin2,ixCmax2
        do ix1=ixCmax1,ixCmin1,-1
          x_VTK(1:3)=zero;
          x_VTK(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0);
          x_VTK(1)=-x_VTK(1)
          do k=1,3
           write(qunit) real(x_VTK(k))
          end do
        end do
        end do
        end do
        write(qunit) length_conn
        do ix3=1,nx3
        do ix2=1,nx2
        do ix1=nx1,1,-1
         write(qunit)&
         (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1-1, &
         (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
         (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1-1,&
         (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1,&
          ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1-1,&
          ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
          ix3*nxC2*nxC1+    ix2*nxC1+ix1-1,&
          ix3*nxC2*nxC1+    ix2*nxC1+ix1
        end do
        end do
        end do
        \}

        write(qunit) length_offsets
        do icel=1,nc
           write(qunit) icel*(2**^ND)
        end do


        {^IFONED VTK_type=3 \}
        {^IFTWOD VTK_type=8 \}
        {^IFTHREED VTK_type=11 \}
        write(qunit) size_int*nc
        do icel=1,nc
         write(qunit) VTK_type
        enddo
       endif
       !! symmetry axis/plane x=0 mirrored part------------------------------end
    endif
  end do
 endif
end do

close(qunit)
open(qunit,file=filename,status='unknown',form='formatted',position='append')

write(qunit,'(a)')'</AppendedData>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)

end subroutine unstructuredvtkBsym
!=============================================================================
subroutine oneblocksym(qunit,userconvert_type)

! this is for turning an AMR run into a single block
! the data will be all on selected level level_io

! this version should work for any dimension
! only writes w_write selected 1:nw variables, also nwauxio
! may use saveprim to switch to primitives
! this version can not work on multiple CPUs
! does not renormalize variables

! header info differs from onegrid below 

! ASCII or binary output

use mod_forest
use mod_global_parameters
integer, intent(in) :: qunit

integer             :: Morton_no,igrid,ix^D,ig^D,level
integer, pointer    :: ig_to_igrid(:^D&,:)
logical             :: fileopen
character(len=80)   :: filename
character(len=20)   :: userconvert_type
integer             :: filenr

character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

double precision :: wval1,xval1
double precision, dimension({^D&1:1},1:nw+nwauxio)   :: wval
double precision, dimension({^D&1:1},1:ndim)         :: xval
double precision:: normconv(0:nw+nwauxio)

integer           :: iw,iiw,writenw,iwrite(1:nw+nwauxio),iigrid
logical :: patchw(ixG^T)
!-----------------------------------------------------------------------------

if(levmin/=levmax.or.level_io<1)then
 if(mype==0) PRINT *,'ONEBLOCK is used only when levmin=levmax'&
       , 'or when level_io is fixed in the parfile'
 call mpistop('level_io<1, oneblocksym')
end if

if(npe>1)then
 if(mype==0) PRINT *,'ONEBLOCK as yet to be parallelized'
 call mpistop('npe>1, oneblocksym')
end if

! only variables selected by w_write will be written out
normconv(0:nw+nwauxio)=one
normconv(0) = length_convert_factor
normconv(1:nw) = w_convert_factor
writenw=count(w_write(1:nw))+nwauxio
iiw=0
do iw =1,nw
 if (.not.w_write(iw))cycle
 iiw=iiw+1
 iwrite(iiw)=iw
end do
if(nwauxio>0)then
  do iw =nw+1,nw+nwauxio
   iiw=iiw+1
   iwrite(iiw)=iw
  end do
endif

if (level_io>0) then
  allocate(ig_to_igrid(ng^D(level_io),0:npe-1))
  ig_to_igrid(:^D&,:)=-1 ! initialize
end if

do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  level=node(plevel_,igrid)
  ig^D=igrid_to_node(igrid,mype)%node%ig^D;
  ig_to_igrid(ig^D,mype)=igrid
end do

call getheadernames(wnamei,xandwnamei,outfilehead)

if (saveprim) then
 do iigrid=1,igridstail; igrid=igrids(iigrid)
  call primitive(ixG^LL,ixG^LL^LSUB1,ps(igrid)%w,ps(igrid)%x)
 end do
else
 if (nwaux>0) then
  do iigrid=1,igridstail; igrid=igrids(iigrid)
   call getaux(.true.,ps(igrid)%w,ixG^LL,ixG^LL^LSUB1,"oneblocksym")
  end do
 end if
end if


Master_cpu_open : if (mype == 0) then
 inquire(qunit,opened=fileopen)
 if (.not.fileopen) then
   ! generate filename
    filenr=snapshotini
    if (.not.convert) filenr=snapshot-1
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".blk"
   select case(userconvert_type)
    case("oneblocksym")
     open(qunit,file=filename,status='unknown')
     write(qunit,*) TRIM(outfilehead)
{^IFONED 
     write(qunit,*)( {^D&(ixMhi^D-ixMlo^D+1)*})*(Morton_stop(npe-1)-Morton_start(0)+1)*2,&
                   ng1(level_io)*(ixMhi1-ixMlo1+1)*2
}
{^IFTWOD 
     write(qunit,*)( {^D&(ixMhi^D-ixMlo^D+1)*})*(Morton_stop(npe-1)-Morton_start(0)+1)*2,&
                   ng1(level_io)*(ixMhi1-ixMlo1+1)*2,ng2(level_io)*(ixMhi2-ixMlo2+1)
}
{^IFTHREED 
     write(qunit,*)( {^D&(ixMhi^D-ixMlo^D+1)*})*(Morton_stop(npe-1)-Morton_start(0)+1)*2,&
                   ng1(level_io)*(ixMhi1-ixMlo1+1)*2,ng2(level_io)*(ixMhi2-ixMlo2+1),&
                   ng3(level_io)*(ixMhi3-ixMlo3+1)
}
     write(qunit,*) global_time*time_convert_factor
    case("oneblocksymB")
     open(qunit,file=filename,form='unformatted',status='unknown')
     write(qunit) outfilehead
{^IFONED 
     write(qunit)( {^D&(ixMhi^D-ixMlo^D+1)*})*(Morton_stop(npe-1)-Morton_start(0)+1)*2,&
                   ng1(level_io)*(ixMhi1-ixMlo1+1)*2
}
{^IFTWOD 
     write(qunit)( {^D&(ixMhi^D-ixMlo^D+1)*})*(Morton_stop(npe-1)-Morton_start(0)+1)*2,&
                   ng1(level_io)*(ixMhi1-ixMlo1+1)*2,ng2(level_io)*(ixMhi2-ixMlo2+1)
}
{^IFTHREED 
     write(qunit)( {^D&(ixMhi^D-ixMlo^D+1)*})*(Morton_stop(npe-1)-Morton_start(0)+1)*2,&
                   ng1(level_io)*(ixMhi1-ixMlo1+1)*2,ng2(level_io)*(ixMhi2-ixMlo2+1),&
                   ng3(level_io)*(ixMhi3-ixMlo3+1)
}
     write(qunit) global_time*time_convert_factor
   end select
 end if
end if Master_cpu_open

{^IFTHREED 
do ig3=1,ng3(level_io)}
   {^NOONED
   do ig2=1,ng2(level_io)}
       do ig1=1,ng1(level_io)
         igrid=ig_to_igrid(ig^D,mype)
         block=>ps
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         ! default (no) normalization for auxiliary variables
         allocate(pwio(igrid)%w(ixG^T,1:nw+nwauxio))
         pwio(igrid)%w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)
         if(nwauxio>=1)then
            call specialvar_output(ixG^LL,ixG^LL^LSUB1,pwio(igrid)%w,ps(igrid)%x,normconv)
         endif
         where(dabs(pwio(igrid)%w(ixG^T,1:nw+nwauxio))<smalldouble**2)
            pwio(igrid)%w(ixG^T,1:nw+nwauxio)=zero
         endwhere
       end do
   {^NOONED
   end do}
{^IFTHREED
end do}

{^IFTHREED
do ig3=1,ng3(level_io)
 do ix3=ixMlo3,ixMhi3}

   {^NOONED
   do ig2=1,ng2(level_io)
     do ix2=ixMlo2,ixMhi2}

!<start---x=0 axisymmetric part--------------------
       do ig1=ng1(level_io),1,-1
         do ix1=ixMhi1,ixMlo1,-1
           igrid=ig_to_igrid(ig^D,mype)
           do iw=1,nw
             if(iwrite(iw)==2 .or. iwrite(iw)==4 .or. iwrite(iw)==7) then
               pwio(igrid)%w(ix^D,iwrite(iw))=-pwio(igrid)%w(ix^D,iwrite(iw))
             endif
           enddo
           Master_write : if(mype==0) then
             select case(userconvert_type)
               case("oneblocksym")
{^IFONED
                 write(qunit,fmt="(100(e14.6))") &
                  xprobmin1-ps(igrid)%x(ix^D,1)*normconv(0),&
                  (pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw)),iw=1,writenw)
}
{^IFTWOD
                 write(qunit,fmt="(100(e14.6))") &
                  xprobmin1-ps(igrid)%x(ix^D,1)*normconv(0),&
                  ps(igrid)%x(ix^D,2)*normconv(0),&
                  (pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw)),iw=1,writenw)
}
{^IFTHREED
                 write(qunit,fmt="(100(e14.6))") &
                  xprobmin1-ps(igrid)%x(ix^D,1)*normconv(0),&
                  ps(igrid)%x(ix^D,2)*normconv(0),ps(igrid)%x(ix^D,3)*normconv(0),&
                  (pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw)),iw=1,writenw)
}
               case("oneblocksymB")
{^IFONED
                 write(qunit) real(xprobmin1-ps(igrid)%x(ix^D,1)*normconv(0)),&
                  (real(pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw))),iw=1,writenw)
}
{^IFTWOD
                 write(qunit) real(xprobmin1-ps(igrid)%x(ix^D,1)*normconv(0)),&
                  real(ps(igrid)%x(ix^D,2)*normconv(0)),&
                  (real(pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw))),iw=1,writenw)
}
{^IFTHREED
                 write(qunit) real(xprobmin1-ps(igrid)%x(ix^D,1)*normconv(0)),&
                  real(ps(igrid)%x(ix^D,2)*normconv(0)),real(ps(igrid)%x(ix^D,3)*normconv(0)),&
                  (real(pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw))),iw=1,writenw)
}
             end select
           end if Master_write
         end do
       end do
!--------x=0 axisymmetric part-----------------end>
       do ig1=1,ng1(level_io)
         do ix1=ixMlo1,ixMhi1
           igrid=ig_to_igrid(ig^D,mype)
           do iw=1,nw
             if(iwrite(iw)==2 .or. iwrite(iw)==4 .or. iwrite(iw)==7) then
               pwio(igrid)%w(ix^D,iwrite(iw))=-pwio(igrid)%w(ix^D,iwrite(iw))
             endif
           enddo
           if(mype==0) then
             select case(userconvert_type)
               case("oneblocksym")
                 write(qunit,fmt="(100(e14.6))") &
                  ps(igrid)%x(ix^D,1:ndim)*normconv(0),&
                  (pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw)),iw=1,writenw)
               case("oneblocksymB")
                 write(qunit) real(ps(igrid)%x(ix^D,1:ndim)*normconv(0)),&
                  (real(pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw))),iw=1,writenw)
             end select
           end if
         end do
       end do
    {^NOONED
     end do
   end do}
 {^IFTHREED
 end do
end do}

{^IFTHREED 
do ig3=1,ng3(level_io)}
   {^NOONED
   do ig2=1,ng2(level_io)}
       do ig1=1,ng1(level_io)
         igrid=ig_to_igrid(ig^D,mype)
         deallocate(pwio(igrid)%w)
       end do
   {^NOONED
   end do}
{^IFTHREED
end do}

close(qunit)

if (saveprim) then
 patchw(ixG^T)=.false.
 do iigrid=1,igridstail; igrid=igrids(iigrid)
  call conserve(ixG^LL,ixG^LL^LSUB1,ps(igrid)%w,ps(igrid)%x,patchw)
 end do
endif

end subroutine oneblocksym
!=============================================================================
subroutine oneblocksym23(qunit,userconvert_type)

! this is for turning an AMR run into a single block
! the data will be all on selected level level_io

! this version should work for any dimension
! only writes w_write selected 1:nw variables, also nwauxio
! may use saveprim to switch to primitives
! this version can not work on multiple CPUs
! does not renormalize variables

! header info differs from onegrid below 

! ASCII or binary output

use mod_forest
use mod_global_parameters
integer, intent(in) :: qunit

integer             :: Morton_no,igrid,ix1,ix2,ix3,ig1,ig2,ig3,level
integer, pointer    :: ig_to_igrid(:^D&,:)
logical             :: fileopen
character(len=80)   :: filename
character(len=20)   :: userconvert_type
integer             :: filenr

character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

double precision :: wval1,xval1
double precision, dimension({^D&1:1},1:nw+nwauxio)   :: wval
double precision, dimension({^D&1:1},1:ndim)         :: xval
double precision :: normconv(0:nw+nwauxio)
double precision :: dx3,lblock3

integer           :: iw,iiw,writenw,iwrite(1:nw+nwauxio),iigrid,nblock3
logical :: patchw(ixG^T)
!-----------------------------------------------------------------------------

if(levmin/=levmax.or.level_io<1)then
 if(mype==0) PRINT *,'ONEBLOCK is used only when levmin=levmax'&
       , 'or when level_io is fixed in the parfile'
 call mpistop('level_io<1, oneblocksym23')
end if

if(npe>1)then
 if(mype==0) PRINT *,'ONEBLOCK as yet to be parallelized'
 call mpistop('npe>1, oneblocksym23')
end if

if(ndim/=2) then
  if(mype==0) print *,'oneblocksym23 is only used in 2D'
  call mpistop('ndim/=2,oneblocksym23')
end if
! only variables selected by w_write will be written out
normconv(0:nw+nwauxio)=one
normconv(0) = length_convert_factor
normconv(1:nw) = w_convert_factor
writenw=count(w_write(1:nw))+nwauxio
iiw=0
do iw =1,nw
 if (.not.w_write(iw))cycle
 iiw=iiw+1
 iwrite(iiw)=iw
end do
if(nwauxio>0)then
  do iw =nw+1,nw+nwauxio
   iiw=iiw+1
   iwrite(iiw)=iw
  end do
endif

if (level_io>0) then
  allocate(ig_to_igrid(ng^D(level_io),0:npe-1))
  ig_to_igrid(:^D&,:)=-1 ! initialize
end if

do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  level=node(plevel_,igrid)
  ig^D=igrid_to_node(igrid,mype)%node%ig^D;
  ig_to_igrid(ig^D,mype)=igrid
end do

call getheadernames(wnamei,xandwnamei,outfilehead)

if (saveprim) then
 do iigrid=1,igridstail; igrid=igrids(iigrid)
  call primitive(ixG^LL,ixG^LL^LSUB1,ps(igrid)%w,ps(igrid)%x)
 end do
else
 if (nwaux>0) then
  do iigrid=1,igridstail; igrid=igrids(iigrid)
   call getaux(.true.,ps(igrid)%w,ixG^LL,ixG^LL^LSUB1,"oneblocksym")
  end do
 end if
end if

nblock3=2*ng1(level_io)
Master_cpu_open : if (mype == 0) then
 inquire(qunit,opened=fileopen)
 if (.not.fileopen) then
   ! generate filename
    filenr=snapshotini
    if (.not.convert) filenr=snapshot-1
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".blk"
   select case(userconvert_type)
    case("oneblocksym23")
     open(qunit,file=filename,status='unknown')
     write(qunit,*) TRIM(outfilehead)
     write(qunit,*)( {^D&(ixMhi^D-ixMlo^D+1)*})*(Morton_stop(npe-1)-Morton_start(0)+1)*&
       2*(ixMhi1-ixMlo1+1)*nblock3,ng1(level_io)*(ixMhi1-ixMlo1+1)*2,ng2(level_io)*&
       (ixMhi2-ixMlo2+1),nblock3*(ixMhi1-ixMlo1+1)
     write(qunit,*)global_time*time_convert_factor
    case("oneblocksym23B")
     open(qunit,file=filename,form='unformatted',status='unknown')
     write(qunit) outfilehead
     write(qunit)( {^D&(ixMhi^D-ixMlo^D+1)*})*(Morton_stop(npe-1)-Morton_start(0)+1)*&
       2*(ixMhi1-ixMlo1+1)*nblock3,ng1(level_io)*(ixMhi1-ixMlo1+1)*2,ng2(level_io)*&
       (ixMhi2-ixMlo2+1),nblock3*(ixMhi1-ixMlo1+1)
     write(qunit)global_time*time_convert_factor
   end select
 end if
end if Master_cpu_open

do ig2=1,ng2(level_io)
    do ig1=1,ng1(level_io)
      igrid=ig_to_igrid(ig^D,mype)
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      ! default (no) normalization for auxiliary variables
      allocate(pwio(igrid)%w(ixG^T,1:nw+nwauxio))
      pwio(igrid)%w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)
      if(nwauxio>=1)then
         call specialvar_output(ixG^LL,ixG^LL^LSUB1,pwio(igrid)%w,ps(igrid)%x,normconv)
      endif
      where(dabs(pwio(igrid)%w(ixG^T,1:nw+nwauxio))<smalldouble**2)
         pwio(igrid)%w(ixG^T,1:nw+nwauxio)=zero
      endwhere
    end do
end do
dx3=dxlevel(1)
lblock3=dble(ixMhi1-ixMlo1+1)*dx3
do ig3=1,nblock3
 do ix3=ixMlo1,ixMhi1
   do ig2=1,ng2(level_io)
     do ix2=ixMlo2,ixMhi2

!<start---x=0 axisymmetric part--------------------
       do ig1=ng1(level_io),1,-1
         do ix1=ixMhi1,ixMlo1,-1
           igrid=ig_to_igrid(ig^D,mype)
           do iw=1,nw
             if(iwrite(iw)==2 .or. iwrite(iw)==4 .or. iwrite(iw)==7) then
               pwio(igrid)%w(ix^D,iwrite(iw))=-pwio(igrid)%w(ix^D,iwrite(iw))
             endif
           enddo
           Master_write : if(mype==0) then
             select case(userconvert_type)
               case("oneblocksym23")
                 write(qunit,fmt="(100(e14.6))") &
                  xprobmin1-ps(igrid)%x(ix^D,1)*normconv(0),&
                  ps(igrid)%x(ix^D,2)*normconv(0),&
                  dble(ig3-1)*lblock3+(dble(ix3-nghostcells)-0.5d0)*dx3,&
                  (pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw)),iw=1,writenw)
               case("oneblocksym23B")
                 write(qunit) &
                  real(xprobmin1-ps(igrid)%x(ix^D,1)*normconv(0)),&
                  real(ps(igrid)%x(ix^D,2)*normconv(0)),&
                  real(dble(ig3-1)*lblock3+(dble(ix3-nghostcells)-0.5d0)*dx3),&
                  (real(pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw))),iw=1,writenw)
             end select
           end if Master_write
         end do
       end do
!--------x=0 axisymmetric part-----------------end>
       do ig1=1,ng1(level_io)
         do ix1=ixMlo1,ixMhi1
           igrid=ig_to_igrid(ig^D,mype)
           do iw=1,nw
             if(iwrite(iw)==2 .or. iwrite(iw)==4 .or. iwrite(iw)==7) then
               pwio(igrid)%w(ix^D,iwrite(iw))=-pwio(igrid)%w(ix^D,iwrite(iw))
             endif
           enddo
           if(mype==0) then
             select case(userconvert_type)
               case("oneblocksym23")
                 write(qunit,fmt="(100(e14.6))") &
                  ps(igrid)%x(ix^D,1:2)*normconv(0),&
                  dble(ig3-1)*lblock3+(dble(ix3-nghostcells)-0.5d0)*dx3,&
                  (pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw)),iw=1,writenw)
               case("oneblocksym23B")
                 write(qunit) real(ps(igrid)%x(ix^D,1:2)*normconv(0)),&
                  real(dble(ig3-1)*lblock3+(dble(ix3-nghostcells)-0.5d0)*dx3),&
                  (real(pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw))),iw=1,writenw)
             end select
           end if
         end do
       end do
     end do
   end do
 end do
end do

{^IFTHREED 
do ig3=1,ng3(level_io)}
   {^NOONED
   do ig2=1,ng2(level_io)}
       do ig1=1,ng1(level_io)
         igrid=ig_to_igrid(ig^D,mype)
         deallocate(pwio(igrid)%w)
       end do
   {^NOONED
   end do}
{^IFTHREED
end do}

close(qunit)

if (saveprim) then
 patchw(ixG^T)=.false.
 do iigrid=1,igridstail; igrid=igrids(iigrid)
  call conserve(ixG^LL,ixG^LL^LSUB1,ps(igrid)%w,ps(igrid)%x,patchw)
 end do
endif

end subroutine oneblocksym23
!=============================================================================
