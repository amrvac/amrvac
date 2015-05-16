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
   call MPI_REDUCE(MPI_IN_PLACE,collapsedData,size(collapsedData),MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
else
   call MPI_REDUCE(collapsedData,collapsedData,size(collapsedData),MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
end if
call MPI_BARRIER(icomm, ierrmpi)

{^NOONED
select case(collapse_type)
case('vti')
   call output_collapsed_vti(dir,normconv)
case('csv')
   call output_collapsed_csv(dir,normconv)
case('default')
   call mpistop("Unknown filetype for collapsed views")
end select
   }{^IFONED
call mpistop("sorry, 1D collapsed output routine not yet implemented (should be easy)...")
}

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
{^NOONED
integer                                           :: ix^DM
double precision                                  :: dxdim^DM, xdim^LIM^DM
}
!-----------------------------------------------------------------------------
if (mype==0) then

! Get coordinates:
{^IFTHREED
select case(dir)
case (1)
   dxdim1 = dx(2,collapseLevel)
   dxdim2 = dx(3,collapseLevel)
   xdim^LIM1=xprob^LIM2; 
   xdim^LIM2=xprob^LIM3; 
case (2)
   dxdim1 = dx(1,collapseLevel)
   dxdim2 = dx(3,collapseLevel)
   xdim^LIM1=xprob^LIM1; 
   xdim^LIM2=xprob^LIM3; 
case (3)
   dxdim1 = dx(1,collapseLevel)
   dxdim2 = dx(2,collapseLevel)
   xdim^LIM1=xprob^LIM1; 
   xdim^LIM2=xprob^LIM2; 
case default
   call mpistop("slice direction not clear in output_collapsed_csv")
end select
}
{^IFTWOD
select case(dir)
case (1)
   dxdim1 = dx(2,collapseLevel)
   xdim^LIM1=xprob^LIM2; 
case (2)
   dxdim1 = dx(1,collapseLevel)
   xdim^LIM1=xprob^LIM1; 
case default
   call mpistop("slice direction not clear in output_collapsed_csv")
end select
}

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
{^NOONED
{^DM& do ix^DMB = 1,myshape(^DMB)\}
   write(unitcollapse,'(200(1pe20.12))') &
        {^DM& dxdim^DM*dble(ix^DM)+xdimmin^DM}, &
        (collapsedData(ix^DM,iw)*normconv(iw),iw=1,nw+nwauxio)
{^DM& enddo\}
}
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
{^NOONED
integer                                           :: ix^DM
double precision                                  :: dxdim^DM, xdim^LIM^DM
}
double precision                                  :: origin(1:3), spacing(1:3)
integer                                           :: wholeExtent(1:6), size_single, length, size_length
integer*8                                         :: offset
character                                         :: buf
!-----------------------------------------------------------------------------
if (mype==0) then
offset=0
size_single=4
size_length=4
! Get coordinates:
{^IFTHREED
select case(dir)
case (1)
   dxdim1 = dx(2,1)*2.0d0**(1-collapseLevel)
   dxdim2 = dx(3,1)*2.0d0**(1-collapseLevel)
   xdim^LIM1=xprob^LIM2; 
   xdim^LIM2=xprob^LIM3; 
case (2)
   dxdim1 = dx(1,1)*2.0d0**(1-collapseLevel)
   dxdim2 = dx(3,1)*2.0d0**(1-collapseLevel)
   xdim^LIM1=xprob^LIM1; 
   xdim^LIM2=xprob^LIM3; 
case (3)
   dxdim1 = dx(1,1)*2.0d0**(1-collapseLevel)
   dxdim2 = dx(2,1)*2.0d0**(1-collapseLevel)
   xdim^LIM1=xprob^LIM1; 
   xdim^LIM2=xprob^LIM2; 
case default
   call mpistop("slice direction not clear in output_collapsed_vti")
end select
}
{^IFTWOD
select case(dir)
case (1)
   dxdim1 = dx(2,1)*2.0d0**(1-collapseLevel)
   xdim^LIM1=xprob^LIM2; 
case (2)
   dxdim1 = dx(1,1)*2.0d0**(1-collapseLevel)
   xdim^LIM1=xprob^LIM1; 
case default
   call mpistop("slice direction not clear in output_collapsed_vti")
end select
}

origin=0
spacing=zero
wholeExtent=0
myshape = shape(collapsedData)
{^NOONED
length = ^DM&myshape(^DM)*
length = length*size_single
{^DM&wholeExtent(^DM*2)=myshape(^DM); }
{^DM&spacing(^DM) = dxdim^DM; }
{^DM&origin(^D) = xdimmin^DM; }
}

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
{#IFDEF BIGENDIAN write(unitcollapse,'(a)')' version="0.1" byte_order="BigEndian">'}
{#IFNDEF BIGENDIAN write(unitcollapse,'(a)')' version="0.1" byte_order="LittleEndian">'}
write(unitcollapse,'(a,3(1pe14.6),a,6(i10),a,3(1pe14.6),a)')'  <ImageData Origin="',origin,'" WholeExtent="',wholeExtent,'" Spacing="',spacing,'">'
write(unitcollapse,'(a)')'<FieldData>'
write(unitcollapse,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(unitcollapse,*) real(t*normt)
write(unitcollapse,'(a)')'</DataArray>'
write(unitcollapse,'(a)')'</FieldData>'

! we write one VTK PIECE
write(unitcollapse,'(a,6(i10),a)') &
     '<Piece Extent="',wholeExtent,'">'
write(unitcollapse,'(a)')'<CellData>'

do iw=1,nw+nwauxio
  if(iw<=nw) then 
    if(.not.writew(iw)) cycle
  endif
  write(unitcollapse,'(a,a,a,i16,a)')&
       '<DataArray type="Float32" Name="',TRIM(wnamei(iw)),'" format="appended" offset="',offset,'"/>'
  offset = offset + length + size_length
enddo

write(unitcollapse,'(a)')'</CellData>'
write(unitcollapse,'(a)')'</Piece>'
write(unitcollapse,'(a)')'</ImageData>'
write(unitcollapse,'(a)')'<AppendedData encoding="raw">'
close(unitcollapse)

open(unitcollapse,file=filename,access='stream',form='unformatted',position='append')
buf='_'
write(unitcollapse) TRIM(buf)

do iw=1,nw+nwauxio
  if(iw<=nw) then 
    if(.not.writew(iw)) cycle
  endif
{#IFNDEF D1
   write(unitcollapse) length
   write(unitcollapse) {^DM&(|}real(collapsedData(ix^DM,iw)*normconv(iw)),{ix^DM=1,myshape(^DM))}
}
enddo

close(unitcollapse)
open(unitcollapse,file=filename,status='unknown',form='formatted',position='append')
write(unitcollapse,'(a)')'</AppendedData>'
write(unitcollapse,'(a)')'</VTKFile>'
close(unitcollapse)
 
end if

end subroutine output_collapsed_vti
!=============================================================================
subroutine allocate_collapsed(dir)
include 'amrvacdef.f'
integer, intent(in)                               :: dir
integer                                           :: dim^D
integer                                           :: nx^D
!-----------------------------------------------------------------------------
! allocate array for the collapsed data:
! number of cells per grid.
nx^D=ixMhi^D-ixMlo^D+1;
{dim^D=ng^D(1)*2**(collapselevel-1)*nx^D\}
{^IFTHREED
select case(dir)
case (1)
   allocate(collapsedData(&
        dim2,dim3,nw+nwauxio))
case (2)
   allocate(collapsedData(&
        dim1,dim3,nw+nwauxio))
case (3)
   allocate(collapsedData(&
        dim1,dim2,nw+nwauxio))
case default
   call mpistop("slice direction not clear in allocate_collapsed")
end select
}
{^IFTWOD
select case(dir)
case (1)
   allocate(collapsedData(&
        dim2,nw+nwauxio))
case (2)
   allocate(collapsedData(&
        dim1,nw+nwauxio))
case default
   call mpistop("slice direction not clear in allocate_collapsed")
end select
}
{^IFONED
   allocate(collapsedData(nw+nwauxio))
}
collapsedData = zero

end subroutine allocate_collapsed
!=============================================================================
subroutine integrate_subnode(igrid,jgrid,dir)
use mod_forest, only: tree_node_ptr, igrid_to_node
include 'amrvacdef.f'
integer, intent(in)                        :: igrid, jgrid, dir
! .. local ..
type(tree_node_ptr)                        :: tree
integer                                    :: nx^D
integer                                    :: ig^D, level, ig^Dtargetmin, ig^Dtargetmax
integer                                    :: ix^Dtarget^LIM, idim^Dtarget^LIM
integer                                    :: ixMdim^LLIM^D, ix^Dorig, ix^D
{^NOONED
integer                                    ::  igdim^DM
}
!-----------------------------------------------------------------------------
tree%node => igrid_to_node(igrid, mype)%node
{^D& ig^D = tree%node%ig^D; }
level = tree%node%level
! number of cells per grid.
nx^D=ixMhi^D-ixMlo^D+1;

{^D&
if (level .gt. collapseLevel) then
   ig^Dtargetmin = int(dble(ig^D-1)/2.0d0**(level-collapseLevel))+1
   ig^Dtargetmax = ig^Dtargetmin
else if (level .lt. collapseLevel) then
   ig^Dtargetmin = int(2.0d0**(collapseLevel-level))*(ig^D-1)+1
   ig^Dtargetmax = int(2.0d0**(collapseLevel-level))*ig^D
else
   ig^Dtargetmin = ig^D
   ig^Dtargetmax = ig^D
end if
ix^Dtargetmin = nx^D*(ig^Dtargetmin-1)+1
ix^Dtargetmax = nx^D*ig^Dtargetmax
\}

{^IFTHREED
select case(dir)
case (1)
   igdim1=(ig2-1)*nx2
   igdim2=(ig3-1)*nx3
   idim1target^LIM=ix2target^LIM; 
   idim2target^LIM=ix3target^LIM;
   ixMdim^LLIM1=ixM^LLIM2;
   ixMdim^LLIM2=ixM^LLIM3;
case (2)
   igdim1=(ig1-1)*nx1
   igdim2=(ig3-1)*nx3
   idim1target^LIM=ix1target^LIM; 
   idim2target^LIM=ix3target^LIM;
   ixMdim^LLIM1=ixM^LLIM1;
   ixMdim^LLIM2=ixM^LLIM3;
case (3)
   igdim1=(ig1-1)*nx1
   igdim2=(ig2-1)*nx2
   idim1target^LIM=ix1target^LIM; 
   idim2target^LIM=ix2target^LIM;
   ixMdim^LLIM1=ixM^LLIM1;
   ixMdim^LLIM2=ixM^LLIM2;
case default
   call mpistop("slice direction not clear in integrate_subnode")
end select

if (level .ge. collapseLevel) then
   do ix1orig = ixMdimlo1,ixMdimhi1
      do ix2orig = ixMdimlo2,ixMdimhi2
{^DM& ix^DM = int(dble(ix^DMorig-dixB+igdim^DM-1)*2.0d0**(collapseLevel-level))+1\}
         collapsedData(ix1,ix2,1:nw+nwauxio) = collapsedData(ix1,ix2,1:nw+nwauxio) &
              + pw_sub(jgrid)%w(ix1orig,ix2orig,1:nw+nwauxio) / 2.0d0**(2*(level-collapseLevel))
      end do
   end do
else
{^DM&
   do ix^DM = idim^DMtargetmin,idim^DMtargetmax\}
 {^DM& ix^DMorig = int(dble(ix^DM-idim^DMtargetmin)/2.0d0**(collapseLevel-level))+1+dixB \}
 collapsedData(ix1,ix2,1:nw+nwauxio) = collapsedData(ix1,ix2,1:nw+nwauxio) + pw_sub(jgrid)%w(ix1orig,ix2orig,1:nw+nwauxio)
 {^DM& enddo\}
end if
}
{^IFTWOD
select case(dir)
case (1)
   igdim1=(ig2-1)*nx2
   idim1target^LIM=ix2target^LIM; 
   ixMdim^LLIM1=ixM^LLIM2;
case (2)
   igdim1=(ig1-1)*nx1
   idim1target^LIM=ix1target^LIM; 
   ixMdim^LLIM1=ixM^LLIM1;
case default
   call mpistop("slice direction not clear in integrate_subnode")
end select

if (level .ge. collapseLevel) then
   do ix1orig = ixMdimlo1,ixMdimhi1
{^DM& ix^DM = int(dble(ix^DMorig-dixB+igdim^DM-1)*2.0d0**(collapseLevel-level))+1\}
         collapsedData(ix1,1:nw+nwauxio) = collapsedData(ix1,1:nw+nwauxio) &
              + pw_sub(jgrid)%w(ix1orig,1:nw+nwauxio) / 2.0d0**(level-collapseLevel)
   end do
else
{^DM&
   do ix^DM = idim^DMtargetmin,idim^DMtargetmax\}
 {^DM& ix^DMorig = int(dble(ix^DM-idim^DMtargetmin)/2.0d0**(collapseLevel-level))+1+dixB \}
 collapsedData(ix1,1:nw+nwauxio) = collapsedData(ix1,1:nw+nwauxio) + pw_sub(jgrid)%w(ix1orig,1:nw+nwauxio)
 {^DM& enddo\}
end if
}
{^IFONED
collapsedData(1:nw+nwauxio) = collapsedData(1:nw+nwauxio) + pw_sub(jgrid)%w(1:nw+nwauxio)
}

end subroutine integrate_subnode
!=============================================================================
subroutine collapse_subnode(igrid,jgrid,dir,normconv)
include 'amrvacdef.f'
integer, intent(in) :: igrid, jgrid, dir
double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv 
! .. local ..
integer                :: ix, iw
double precision       :: dx^D
double precision, dimension(ixG^T,1:nw+nwauxio)   :: w
!-----------------------------------------------------------------------------
pw_sub(jgrid)%w=zero
dx^D=rnode(rpdx^D_,igrid);

w(ixG^T,1:nw)=pw(igrid)%w(ixG^T,1:nw)
if(saveprim) then 
   call primitive(ixG^LL,ixG^LL,w,px(igrid)%x)
   normconv(0:nw)=normvar(0:nw)
else
  normconv(0:nw)=one
end if

if(nwauxio>0)then
! auxiliary io variables can be computed and added by user
typelimiter=typelimiter1(node(plevel_,igrid))
typegradlimiter=typegradlimiter1(node(plevel_,igrid))
^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  if (.not.slab) mygeo => pgeo(igrid)
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
     myB0      => pB0_cell(igrid)
     {^D&myB0_face^D => pB0_face^D(igrid)\}
  end if
  ! default (no) normalization for auxiliary variables
  normconv(nw+1:nw+nwauxio)=one
  call specialvar_output(ixG^LL,ixM^LL,w,px(igrid)%x,normconv)
endif
{^IFMHDPHYS
if (B0field) w(ixM^T,b1_:b0_+ndir) = w(ixM^T,b1_:b0_+ndir) + pB0_cell(igrid)%w(ixM^T,1:ndir)
{#IFDEF ENERGY
if((.not.saveprim) .and. B0field) then
   w(ixM^T,e_)=w(ixM^T,e_) &
           +half*( ^C&pB0_cell(igrid)%w(ixM^T,b^C_-b0_)**2+ ) &
           + ( ^C&w(ixM^T,b^C_)*pB0_cell(igrid)%w(ixM^T,b^C_-b0_)+ )
endif
}
}
{^IFTHREED
select case (dir)
case (1)
   do ix=ixMlo1,ixMhi1
      pw_sub(jgrid)%w(ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:nw+nwauxio) = &
           pw_sub(jgrid)%w(ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:nw+nwauxio) &
           + w(ix,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:nw+nwauxio) * dx1
   end do
   px_sub(jgrid)%x(ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:ndim) = &
        px(igrid)%x(ixMlo1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:ndim)
case (2)
   do ix=ixMlo2,ixMhi2
      pw_sub(jgrid)%w(ixMlo1:ixMhi1,ixMlo3:ixMhi3,1:nw+nwauxio) = &
           pw_sub(jgrid)%w(ixMlo1:ixMhi1,ixMlo3:ixMhi3,1:nw+nwauxio) &
           + w(ixMlo1:ixMhi1,ix,ixMlo3:ixMhi3,1:nw+nwauxio) * dx2
   end do
   px_sub(jgrid)%x(ixMlo1:ixMhi1,ixMlo3:ixMhi3,1:ndim) = &
        px(igrid)%x(ixMlo1:ixMhi1,ixMlo2,ixMlo3:ixMhi3,1:ndim) 
case (3)
   do ix=ixMlo3,ixMhi3
      pw_sub(jgrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,1:nw+nwauxio) = &
           pw_sub(jgrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,1:nw+nwauxio) &
           + w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ix,1:nw+nwauxio) * dx3
   end do
      px_sub(jgrid)%x(ixMlo1:ixMhi1,ixMlo2:ixMhi2,1:ndim) = &
           px(igrid)%x(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3,1:ndim) 
case default
   print*, 'subnode, dir: ', dir
   call mpistop("slice direction not clear in collapse_subnode")
end select
}
{^IFTWOD
select case (dir)
case (1)
   do ix=ixMlo1,ixMhi1
      pw_sub(jgrid)%w(ixMlo2:ixMhi2,1:nw+nwauxio) = &
           pw_sub(jgrid)%w(ixMlo2:ixMhi2,1:nw+nwauxio) &
           + w(ix,ixMlo2:ixMhi2,1:nw+nwauxio) * dx1
   end do
   px_sub(jgrid)%x(ixMlo2:ixMhi2,1:ndim) = &
        px(igrid)%x(ixMlo1,ixMlo2:ixMhi2,1:ndim)
case (2)
   do ix=ixMlo2,ixMhi2
      pw_sub(jgrid)%w(ixMlo1:ixMhi1,1:nw+nwauxio) = &
           pw_sub(jgrid)%w(ixMlo1:ixMhi1,1:nw+nwauxio) &
           + w(ixMlo1:ixMhi1,ix,1:nw+nwauxio) * dx2
   end do
   px_sub(jgrid)%x(ixMlo1:ixMhi1,1:ndim) = &
        px(igrid)%x(ixMlo1:ixMhi1,ixMlo2,1:ndim) 
case default
   call mpistop("slice direction not clear in collapse_subnode")
end select
}
{^IFONED   
do ix=ixMlo1,ixMhi1
   pw_sub(jgrid)%w(1:nw+nwauxio) = pw_sub(jgrid)%w(1:nw+nwauxio) + w(ix,1:nw+nwauxio) * dx1
end do
px_sub(jgrid)%x(1:ndim) = px(igrid)%x(ixMlo1,1:ndim)
}

end subroutine collapse_subnode
!=============================================================================
