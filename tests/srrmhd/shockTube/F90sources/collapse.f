!=============================================================================
subroutine write_collapsed
include 'amrvacdef.f'
! Writes a collapsed view of the data integrated over one grid-direction.  
! E.g. column density maps.  
! Uses flat interpolation throughout.
! by Oliver Porth
! 6.Nov 2013
integer :: idir
logical, save :: firstcollapse=.true.
!-----------------------------------------------------------------------------
if (firstcollapse) then
   icollapse=collapseNext
   firstcollapse=.false.
end if

do idir=1, ndim
   if (collapse(idir)) call put_collapse(idir)
end do

icollapse=icollapse+1
end subroutine write_collapsed
!=============================================================================
subroutine put_collapse(dir)
use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
include 'amrvacdef.f'
integer, intent(in)                               :: dir
! .. local ..
integer                                           :: jgrid, igrid, Morton_no
double precision,dimension(0:nw+nwauxio)          :: normconv 
!-----------------------------------------------------------------------------

call allocate_collapsed(dir)

jgrid=0
do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   jgrid=jgrid+1
   call alloc_subnode(jgrid,dir,nwauxio)
   call collapse_subnode(igrid,jgrid,dir,normconv)
   call integrate_subnode(igrid,jgrid,dir)
end do

! Reduce to head-node:
if (mype==0) then
   call MPI_REDUCE(MPI_IN_PLACE,collapsedData,size(collapsedData),&
      MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
else
   call MPI_REDUCE(collapsedData,collapsedData,size(collapsedData),&
      MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
end if
call MPI_BARRIER(icomm, ierrmpi)


call mpistop("sorry, 1D collapsed output routine not yet implemented (should &
   be easy)...")


! If we need the subnodes later, remove deallocation here:
do jgrid=1,Morton_stop(mype)-Morton_start(mype)+1
   call dealloc_subnode(jgrid)
end do
deallocate(collapsedData)

end subroutine put_collapse
!=============================================================================
subroutine output_collapsed_csv(dir,normconv)
include 'amrvacdef.f'
integer, intent(in)                               :: dir
double precision,dimension(0:nw+nwauxio),intent(in):: normconv 
character(len=1024) :: filename, outfilehead, line
logical             :: fileopen
integer                                           :: iw
character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
integer, dimension(ndim)                          :: myshape

!-----------------------------------------------------------------------------
if (mype==0) then

! Get coordinates:



 inquire(unitcollapse,opened=fileopen)
 if(.not.fileopen)then
      ! generate filename: 
      write(filename,"(a,i1.1,a,i1.1,a,i4.4,a)") TRIM(filenameout)//'_d',dir,'_l',collapseLevel,'_n',icollapse,'.csv'
      open(unitcollapse,file=filename,status='unknown',form='formatted')
   end if
   ! get and write the header: 
   call getheadernames(wnamei,xandwnamei,outfilehead)
   line=''
   do iw=1,ndim+nw+nwauxio-1
      if (iw .eq. dir) cycle
      line = trim(line)//trim(xandwnamei(iw))//', '
   end do
   line = trim(line)//trim(xandwnamei(ndim+nw+nwauxio))
   write(unitcollapse,'(a)')trim(line)
   myshape = shape(collapsedData)

close(unitcollapse)

end if
end subroutine output_collapsed_csv
!=============================================================================
subroutine output_collapsed_vti(dir,normconv)
include 'amrvacdef.f'
integer, intent(in)                               :: dir
double precision,dimension(0:nw+nwauxio),intent(in):: normconv 
character(len=1024) :: filename, outfilehead, line
logical             :: fileopen
integer                                           :: iw
character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
integer, dimension(ndim)                          :: myshape

double precision                                  :: origin(1:3), spacing(1:3)
integer                                           :: wholeExtent(1:6),&
    size_single, length, size_length
integer*8                                         :: offset
character                                         :: buf
!-----------------------------------------------------------------------------
if (mype==0) then
offset=0
size_single=4
size_length=4
! Get coordinates:



origin=0
spacing=zero
wholeExtent=0
myshape = shape(collapsedData)


 inquire(unitcollapse,opened=fileopen)
 if(.not.fileopen)then
      ! generate filename: 
      write(filename,"(a,i1.1,a,i1.1,a,i4.4,a)") TRIM(filenameout)//'_d',dir,'_l',collapseLevel,'_n',icollapse,'.vti'
      open(unitcollapse,file=filename,status='unknown',form='formatted')
 end if
! get the header: 
call getheadernames(wnamei,xandwnamei,outfilehead)

! generate xml header
write(unitcollapse,'(a)')'<?xml version="1.0"?>'
write(unitcollapse,'(a)',advance='no') '<VTKFile type="ImageData"'

 write(unitcollapse,'(a)')' version="0.1" byte_order="LittleEndian">'
write(unitcollapse,'(a,3(1pe14.6),a,6(i10),a,3(1pe14.6),a)')&
   '  <ImageData Origin="',origin,'" WholeExtent="',wholeExtent,'" Spacing="',&
   spacing,'">'
write(unitcollapse,'(a)')'<FieldData>'
write(unitcollapse,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
   'NumberOfTuples="1" format="ascii">'
write(unitcollapse,*) real(t*normt)
write(unitcollapse,'(a)')'</DataArray>'
write(unitcollapse,'(a)')'</FieldData>'

! we write one VTK PIECE
write(unitcollapse,'(a,6(i10),a)') '<Piece Extent="',wholeExtent,'">'
write(unitcollapse,'(a)')'<CellData>'

do iw=1,nw+nwauxio
  if(iw<=nw) then 
    if(.not.writew(iw)) cycle
  endif
  write(unitcollapse,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
     TRIM(wnamei(iw)),'" format="appended" offset="',offset,'"/>'
  offset = offset + length + size_length
enddo

write(unitcollapse,'(a)')'</CellData>'
write(unitcollapse,'(a)')'</Piece>'
write(unitcollapse,'(a)')'</ImageData>'
write(unitcollapse,'(a)')'<AppendedData encoding="raw">'
close(unitcollapse)

open(unitcollapse,file=filename,access='stream',form='unformatted',position&
   ='append')
buf='_'
write(unitcollapse) TRIM(buf)

do iw=1,nw+nwauxio
  if(iw<=nw) then 
    if(.not.writew(iw)) cycle
  endif

enddo

close(unitcollapse)
open(unitcollapse,file=filename,status='unknown',form='formatted',position&
   ='append')
write(unitcollapse,'(a)')'</AppendedData>'
write(unitcollapse,'(a)')'</VTKFile>'
close(unitcollapse)
 
end if

end subroutine output_collapsed_vti
!=============================================================================
subroutine allocate_collapsed(dir)
include 'amrvacdef.f'
integer, intent(in)                               :: dir
integer                                           :: dim1
integer                                           :: nx1
!-----------------------------------------------------------------------------
! allocate array for the collapsed data:
! number of cells per grid.
nx1=ixMhi1-ixMlo1+1;
dim1=ng1(1)*2**(collapselevel-1)*nx1



   allocate(collapsedData(nw+nwauxio))

collapsedData = zero

end subroutine allocate_collapsed
!=============================================================================
subroutine integrate_subnode(igrid,jgrid,dir)
use mod_forest, only: tree_node_ptr, igrid_to_node
include 'amrvacdef.f'
integer, intent(in)                        :: igrid, jgrid, dir
! .. local ..
type(tree_node_ptr)                        :: tree
integer                                    :: nx1
integer                                    :: ig1, level, ig1targetmin,&
    ig1targetmax
integer                                    :: ix1targetmin,ix1targetmax,&
    idim1targetmin,idim1targetmax
integer                                    :: ixMdimlo1,ixMdimhi1, ix1orig,&
    ix1

!-----------------------------------------------------------------------------
tree%node => igrid_to_node(igrid, mype)%node
 ig1 = tree%node%ig1;
level = tree%node%level
! number of cells per grid.
nx1=ixMhi1-ixMlo1+1;


if (level .gt. collapseLevel) then
   ig1targetmin = int(dble(ig1-1)/2.0d0**(level-collapseLevel))+1
   ig1targetmax = ig1targetmin
else if (level .lt. collapseLevel) then
   ig1targetmin = int(2.0d0**(collapseLevel-level))*(ig1-1)+1
   ig1targetmax = int(2.0d0**(collapseLevel-level))*ig1
else
   ig1targetmin = ig1
   ig1targetmax = ig1
end if
ix1targetmin = nx1*(ig1targetmin-1)+1
ix1targetmax = nx1*ig1targetmax





collapsedData(1:nw+nwauxio) = collapsedData(1:nw+nwauxio) + &
   pw_sub(jgrid)%w(1:nw+nwauxio)


end subroutine integrate_subnode
!=============================================================================
subroutine collapse_subnode(igrid,jgrid,dir,normconv)
include 'amrvacdef.f'
integer, intent(in) :: igrid, jgrid, dir
double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv 
! .. local ..
integer                :: ix, iw
double precision       :: dx1
double precision, dimension(ixGlo1:ixGhi1,1:nw+nwauxio)   :: w
!-----------------------------------------------------------------------------
pw_sub(jgrid)%w=zero
dx1=rnode(rpdx1_,igrid);

w(ixGlo1:ixGhi1,1:nw)=pw(igrid)%w(ixGlo1:ixGhi1,1:nw)
if(saveprim) then 
   call primitive(ixGlo1,ixGhi1,ixGlo1,ixGhi1,w,px(igrid)%x)
   normconv(0:nw)=normvar(0:nw)
else
  normconv(0:nw)=one
end if

if(nwauxio>0)then
! auxiliary io variables can be computed and added by user
typelimiter=typelimiter1(node(plevel_,igrid))
typegradlimiter=typegradlimiter1(node(plevel_,igrid))
dxlevel(1)=rnode(rpdx1_,igrid);
  if (.not.slab) mygeo => pgeo(igrid)
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
     myB0      => pB0_cell(igrid)
     myB0_face1 => pB0_face1(igrid)
  end if
  ! default (no) normalization for auxiliary variables
  normconv(nw+1:nw+nwauxio)=one
  call specialvar_output(ixGlo1,ixGhi1,ixMlo1,ixMhi1,w,px(igrid)%x,normconv)
endif



   
do ix=ixMlo1,ixMhi1
   pw_sub(jgrid)%w(1:nw+nwauxio) = pw_sub(jgrid)%w(1:nw+nwauxio) + w(ix,1:nw&
      +nwauxio) * dx1
end do
px_sub(jgrid)%x(1:ndim) = px(igrid)%x(ixMlo1,1:ndim)


end subroutine collapse_subnode
!=============================================================================
