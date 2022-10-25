!=============================================================================
subroutine generate_plotfile
use mod_usr_methods, only: usr_special_convert
use mod_global_parameters
use mod_ghostcells_update
use mod_physics, only: phys_req_diagonal,phys_te_images
use mod_convert, only: convert_all
use mod_thermal_emission
!-----------------------------------------------------------------------------

if(mype==0.and.level_io>0) write(unitterm,*)'reset tree to fixed level=',level_io
if(level_io>0 .or. level_io_min.ne.1 .or. level_io_max.ne.nlevelshi) then 
   call resettree_convert
else if(.not. phys_req_diagonal) then
   call getbc(global_time,0.d0,ps,iwstart,nwgc)
end if

select case(convert_type)
  case('tecplot','tecplotCC','tecline')
   call tecplot(unitconvert)
  case('tecplotmpi','tecplotCCmpi','teclinempi')
   call tecplot_mpi(unitconvert)
  case('vtu','vtuCC')
   call unstructuredvtk(unitconvert)
  case('vtumpi','vtuCCmpi')
   call unstructuredvtk_mpi(unitconvert)
  case('vtuB','vtuBCC','vtuBmpi','vtuBCCmpi')
   call unstructuredvtkB(unitconvert)
  case('vtuB64','vtuBCC64','vtuBmpi64','vtuBCCmpi64')
   call unstructuredvtkB64(unitconvert)
  {^IFTWOD
  case('vtuB23','vtuBCC23')
   call unstructuredvtkB23(unitconvert)
  case('vtuBsym23','vtuBCCsym23')
   call unstructuredvtkBsym23(unitconvert)
  }
  case('pvtumpi','pvtuCCmpi')
   call punstructuredvtk_mpi(unitconvert)
  case('pvtuBmpi','pvtuBCCmpi')
   call punstructuredvtkB_mpi(unitconvert)
  case('vtimpi','vtiCCmpi')
   call ImageDataVtk_mpi(unitconvert)
  case('onegrid','onegridmpi')
   call onegrid(unitconvert)
  case('oneblock','oneblockB')
   call oneblock(unitconvert)
  case('EIvtiCCmpi','ESvtiCCmpi','SIvtiCCmpi','EIvtuCCmpi','ESvtuCCmpi','SIvtuCCmpi')
    ! output synthetic euv emission
    if (ndim==3 .and. slab .and. associated(phys_te_images)) then
      call phys_te_images
    endif
  case('dat_generic_mpi')
    ! it is the responsability of the physcics module to pouplate the array of 
    ! conversion methods by calling  
    call convert_all()
  case('user','usermpi')
     if (.not. associated(usr_special_convert)) then
        call mpistop("usr_special_convert not defined")
     else
        call usr_special_convert(unitconvert)
     end if
  case default
   call mpistop("Error in generate_plotfile: Unknown convert_type")
end select



end subroutine generate_plotfile
!=============================================================================
subroutine oneblock(qunit)
! this is for turning an AMR run into a single block
! the data will be all on selected level level_io

! this version should work for any dimension
! only writes w_write selected 1:nw variables, also nwauxio
! may use saveprim to switch to primitives
! this version can not work on multiple CPUs
! does not renormalize variables

! header info differs from onegrid below 

! ASCII or binary output

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid, igrid_to_node
use mod_global_parameters
use mod_usr_methods, only: usr_aux_output
use mod_physics
use mod_calculate_xw
integer, intent(in) :: qunit

integer             :: Morton_no,igrid,ix^D,ig^D,level
integer, pointer    :: ig_to_igrid(:^D&,:)
logical             :: fileopen,writeblk(max_blocks)
character(len=80)   :: filename
integer             :: filenr,ncells^D,ncellx^D,jg^D,jig^D

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

double precision :: wval1,xval1
double precision, dimension({^D&1:1},1:nw+nwauxio)   :: wval
double precision, dimension({^D&1:1},1:ndim)         :: xval
double precision:: normconv(0:nw+nwauxio)

integer           :: iw,iiw,writenw,iwrite(1:nw+nwauxio),iigrid,idim
logical :: patchw(ixG^T)
!-----------------------------------------------------------------------------

if(level_io<1)then
 call mpistop('please specify level_io>0 for usage with oneblock')
end if

if(autoconvert)then
 call mpistop('Set autoconvert=F and convert oneblock data manually')
end if

if(npe>1)then
 if(mype==0) PRINT *,'ONEBLOCK as yet to be parallelized'
 call mpistop('npe>1, oneblock')
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

allocate(ig_to_igrid(ng^D(level_io),0:npe-1))
ig_to_igrid=-1
writeblk=.false.
do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  level=node(plevel_,igrid)
  ig^D=igrid_to_node(igrid,mype)%node%ig^D;
  ig_to_igrid(ig^D,mype)=igrid
  if(({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
        *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
       <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})) then
    writeblk(igrid)=.true.
  end if
end do

call getheadernames(wnamei,xandwnamei,outfilehead)
ncells^D=0;
ncellx^D=ixMhi^D-ixMlo^D+1\
{do ig^D=1,ng^D(level_io)\}
  igrid=ig_to_igrid(ig^D,mype)
  if(writeblk(igrid)) go to 20
{end do\}
20 continue
jg^D=ig^D;
{
jig^DD=jg^DD;
do ig^D=1,ng^D(level_io)
  jig^D=ig^D
  igrid=ig_to_igrid(jig^DD,mype)
  if(writeblk(igrid)) ncells^D=ncells^D+ncellx^D
end do
\}

do iigrid=1,igridstail; igrid=igrids(iigrid)
  if(.not.writeblk(igrid)) cycle
  block=>ps(igrid)
  if (nwauxio > 0) then
    if (.not. associated(usr_aux_output)) then
      call mpistop("usr_aux_output not defined")
    else
      call usr_aux_output(ixG^LL,ixM^LL^LADD1, &
            ps(igrid)%w,ps(igrid)%x,normconv)
    end if
  end if
end do

if (saveprim) then
  do iigrid=1,igridstail; igrid=igrids(iigrid)
    if (.not.writeblk(igrid)) cycle
    block=>ps(igrid)
    call phys_to_primitive(ixG^LL,ixG^LL^LSUB1,ps(igrid)%w,ps(igrid)%x)
    if(B0field) then
      ! add background magnetic field B0 to B
      ps(igrid)%w(ixG^T,iw_mag(:))=ps(igrid)%w(ixG^T,iw_mag(:))+ps(igrid)%B0(ixG^T,:,0)
    end if
  end do
else
  do iigrid=1,igridstail; igrid=igrids(iigrid)
    if (.not.writeblk(igrid)) cycle
    block=>ps(igrid)
    if (B0field) then
      ! add background magnetic field B0 to B
      if(phys_energy) &
        ps(igrid)%w(ixG^T,iw_e)=ps(igrid)%w(ixG^T,iw_e)+0.5d0*sum(ps(igrid)%B0(ixG^T,:,0)**2,dim=ndim+1) &
            + sum(ps(igrid)%w(ixG^T,iw_mag(:))*ps(igrid)%B0(ixG^T,:,0),dim=ndim+1)
      ps(igrid)%w(ixG^T,iw_mag(:))=ps(igrid)%w(ixG^T,iw_mag(:))+ps(igrid)%B0(ixG^T,:,0)
    end if
  end do
end if

Master_cpu_open : if (mype == 0) then
 inquire(qunit,opened=fileopen)
 if (.not.fileopen) then
   ! generate filename
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".blk"
   select case(convert_type)
    case("oneblock")
     open(qunit,file=filename,status='unknown')
     write(qunit,*) TRIM(outfilehead)
     write(qunit,*) ncells^D
     write(qunit,*) real(global_time*time_convert_factor)
    case("oneblockB")
     open(qunit,file=filename,form='unformatted',status='unknown')
     write(qunit) outfilehead
     write(qunit) ncells^D
     write(qunit) real(global_time*time_convert_factor)
   end select
 end if
end if Master_cpu_open

{^IFTHREED
do ig3=1,ng3(level_io)
 do ix3=ixMlo3,ixMhi3}

   {^NOONED
   do ig2=1,ng2(level_io)
     do ix2=ixMlo2,ixMhi2}

       do ig1=1,ng1(level_io)
         igrid=ig_to_igrid(ig^D,mype)
         if(.not.writeblk(igrid)) cycle
         do ix1=ixMlo1,ixMhi1
           Master_write : if(mype==0) then
             select case(convert_type)
               case("oneblock")
                 write(qunit,fmt="(100(e14.6))") &
                  ps(igrid)%x(ix^D,1:ndim)*normconv(0),&
                  (ps(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw)),iw=1,writenw)
               case("oneblockB")
                 write(qunit) real(ps(igrid)%x(ix^D,1:ndim)*normconv(0)),&
                  (real(ps(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw))),iw=1,writenw)
             end select
           end if Master_write
         end do
       end do
    {^NOONED
     end do
   end do}
 {^IFTHREED
 end do
end do}

close(qunit)

end subroutine oneblock
!=============================================================================
subroutine onegrid(qunit)

! this is for turning an AMR run into a single grid
! this version should work for any dimension, can be in parallel
! in 1D, should behave much like oneblock, except for header info

! only writes all 1:nw variables, no nwauxio
! may use saveprim to switch to primitives
! this version can work on multiple CPUs
! does not renormalize variables
! ASCII output

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters
use mod_physics
use mod_calculate_xw
integer, intent(in) :: qunit

integer             :: itag,Morton_no,igrid,ix^D,iw
logical             :: fileopen
character(len=80)   :: filename
integer             :: filenr

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

!.. MPI variables ..
integer           :: igrid_recv,ipe
double precision  :: w_recv(ixG^T,1:nw),x_recv(ixG^T,1:ndim)
integer, allocatable :: intstatus(:,:)

!-----------------------------------------------------------------------------

if(nwauxio>0)then
 if(mype==0) PRINT *,'ONEGRID to be used without nwauxio'
 call mpistop('nwauxio>0, onegrid')
end if

if(saveprim)then
 if(mype==0.and.nwaux>0) PRINT *,'warning: ONEGRID used with saveprim, check auxiliaries'
end if



Master_cpu_open : if (mype == 0) then
 call getheadernames(wnamei,xandwnamei,outfilehead)
 write(outfilehead,'(a)') "#"//" "//TRIM(outfilehead)
 inquire(qunit,opened=fileopen)
 if (.not.fileopen) then
   ! generate filename
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".blk"
   open(qunit,file=filename,status='unknown')
 end if
 write(qunit,"(a)")outfilehead
 write(qunit,"(i7)") ( {^D&(ixMhi^D-ixMlo^D+1)*} )*(Morton_stop(npe-1)-Morton_start(0)+1)
end if Master_cpu_open

do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  if(saveprim) call phys_to_primitive(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x)
  if (mype/=0)then
      itag=Morton_no
      call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
      call MPI_SEND(ps(igrid)%x,1,type_block_xcc_io, 0,itag,icomm,ierrmpi)
      itag=igrid
      call MPI_SEND(ps(igrid)%w,1,type_block_io, 0,itag,icomm,ierrmpi)
  else
   {do ix^DB=ixMlo^DB,ixMhi^DB\}
      do iw=1,nw
        if( dabs(ps(igrid)%w(ix^D,iw)) < 1.0d-32 ) ps(igrid)%w(ix^D,iw) = zero
      enddo
       write(qunit,fmt="(100(e14.6))") ps(igrid)%x(ix^D,1:ndim)&
                                     ,ps(igrid)%w(ix^D,1:nw)
   {end do\}
  end if
end do   

if(mype==0.and.npe>1) allocate(intstatus(MPI_STATUS_SIZE,1))

Manycpu : if (npe>1) then
 if (mype==0) then
  loop_cpu : do ipe =1, npe-1
   loop_Morton : do Morton_no=Morton_start(ipe),Morton_stop(ipe)
         itag=Morton_no
         call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
         call MPI_RECV(x_recv,1,type_block_xcc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
         itag=igrid_recv
         call MPI_RECV(w_recv,1,type_block_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
         {do ix^DB=ixMlo^DB,ixMhi^DB\}
            do iw=1,nw
              if( dabs(ps(igrid)%w(ix^D,iw)) < smalldouble ) ps(igrid)%w(ix^D,iw) = zero
            enddo
            write(qunit,fmt="(100(e14.6))") x_recv(ix^D,1:ndim)&
                                            ,w_recv(ix^D,1:nw)
         {end do\}
   end do loop_Morton
  end do loop_cpu
 end if 
end if Manycpu

if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
  if(mype==0)deallocate(intstatus)
endif

if(mype==0) close(qunit)
end subroutine onegrid 
!============================================================================
subroutine tecplot(qunit)

! output for tecplot (ASCII format)
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors

use mod_global_parameters
use mod_calculate_xw

integer, intent(in) :: qunit

integer::               igrid,iigrid,level,igonlevel,iw,idim,ix^D
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx^D,nxC^D,nodesonlevel,elemsonlevel,ixC^L,ixCC^L

integer ::              nodes, elems


double precision :: x_TEC(ndim), w_TEC(nw+nwauxio)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP
double precision, dimension(0:nw+nwauxio)                   :: normconv
logical :: fileopen,first
character(len=80) :: filename
integer  :: filenr

!!! possible length conflict
character(len=1024) :: tecplothead

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead
!-----------------------------------------------------------------------------

if(npe>1)then
 if(mype==0) PRINT *,'tecplot not parallel, use tecplotmpi'
 call mpistop('npe>1, tecplot')
end if

if(nw/=count(w_write(1:nw)))then
 if(mype==0) PRINT *,'tecplot does not use w_write=F'
 call mpistop('w_write, tecplot')
end if

if(nocartesian)then
 if(mype==0) PRINT *,'tecplot with nocartesian'
endif

inquire(qunit,opened=fileopen)
if (.not.fileopen) then
   ! generate filename    
   filenr=snapshotini
   if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".plt"
   open(qunit,file=filename,status='unknown')
end if

call getheadernames(wnamei,xandwnamei,outfilehead)

write(tecplothead,'(a)') "VARIABLES = "//TRIM(outfilehead)

write(qunit,'(a)') tecplothead(1:len_trim(tecplothead))

NumGridsOnLevel(1:nlevelshi)=0
do level=levmin,levmax
   NumGridsOnLevel(level)=0
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (node(plevel_,igrid)/=level) cycle
      NumGridsOnLevel(level)=NumGridsOnLevel(level)+1
   end do
end do

nx^D=ixMhi^D-ixMlo^D+1;
nxC^D=nx^D+1;

{^IFONED
if(convert_type=='tecline') then
   nodes=0
   elems=0
   do level=levmin,levmax
      nodes=nodes + NumGridsOnLevel(level)*{nxC^D*}
      elems=elems + NumGridsOnLevel(level)*{nx^D*}
   enddo

   write(qunit,"(a,i7,a,1pe12.5,a)") &
         'ZONE T="all levels", I=',elems, &
         ', SOLUTIONTIME=',global_time*time_convert_factor,', F=POINT' 

   igonlevel=0
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      call calc_x(igrid,xC,xCC)
      call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,ixC^L,ixCC^L,.true.)
          {do ix^DB=ixCCmin^DB,ixCCmax^DB\}
            x_TEC(1:ndim)=xCC_TMP(ix^D,1:ndim)*normconv(0)
            w_TEC(1:nw+nwauxio)=wCC_TMP(ix^D,1:nw+nwauxio)*normconv(1:nw+nwauxio)
           !write(qunit,fmt="(100(e14.6))") x_TEC, w_TEC
           write(qunit,fmt="(100(e24.16))") x_TEC, w_TEC
       {end do\}
    enddo
    close(qunit)
else
}
do level=levmin,levmax
   nodesonlevel=NumGridsOnLevel(level)*{nxC^D*}
   elemsonlevel=NumGridsOnLevel(level)*{nx^D*}
   ! for all tecplot variants coded up here, we let the TECPLOT ZONES coincide
   ! with the AMR grid LEVEL. Other options would be
   !    let each grid define a zone: inefficient for TECPLOT internal workings
   !       hence not implemented
   !    let entire octree define 1 zone: no difference in interpolation 
   !       properties across TECPLOT zones detected as yet, hence not done
   select case(convert_type)
     case('tecplot')
       ! in this option, we store the corner coordinates, as well as the corner
       ! values of all variables (obtained by averaging). This allows POINT packaging, 
       ! and thus we can save full grid info by using one call to calc_grid
       write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevel,', E=',elemsonlevel, &
            ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=POINT, ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
       do iigrid=1,igridstail; igrid=igrids(iigrid);
         if (node(plevel_,igrid)/=level) cycle
         block=>ps(igrid)
         call calc_x(igrid,xC,xCC)
         call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                        ixC^L,ixCC^L,.true.)
         {do ix^DB=ixCmin^DB,ixCmax^DB\}
            x_TEC(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0)
            w_TEC(1:nw+nwauxio)=wC_TMP(ix^D,1:nw+nwauxio)*normconv(1:nw+nwauxio)
            write(qunit,fmt="(100(e14.6))") x_TEC, w_TEC
         {end do\}
       enddo
     case('tecplotCC')
       ! in this option, we store the corner coordinates, and the cell center
       ! values of all variables. Due to this mix of corner/cell center, we must 
       ! use BLOCK packaging, and thus we have enormous overhead by using 
       ! calc_grid repeatedly to merely fill values of cell corner coordinates 
       ! and cell center values per dimension, per variable
       if(ndim+nw+nwauxio>99) call mpistop("adjust format specification in writeout")
       if(nw+nwauxio==1)then
         ! to make tecplot happy: avoid [ndim+1-ndim+1] in varlocation varset
         ! and just set [ndim+1]
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevel,', E=',elemsonlevel, &
            ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
       else
        if(ndim+nw+nwauxio<10) then
         ! difference only in length of integer format specification for ndim+nw+nwauxio
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i1,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevel,', E=',elemsonlevel, &
            ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,'-',ndim+nw+nwauxio,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
        else
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i2,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevel,', E=',elemsonlevel, &
            ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,'-',ndim+nw+nwauxio,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
        endif
       endif
       do idim=1,ndim
         first=(idim==1) 
         do iigrid=1,igridstail; igrid=igrids(iigrid);
            if (node(plevel_,igrid)/=level) cycle
            block=>ps(igrid)
            call calc_x(igrid,xC,xCC)
            call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                           ixC^L,ixCC^L,first)
            write(qunit,fmt="(100(e14.6))") xC_TMP(ixC^S,idim)*normconv(0)
         enddo
       enddo
       do iw=1,nw+nwauxio
         do iigrid=1,igridstail; igrid=igrids(iigrid);
            if (node(plevel_,igrid)/=level) cycle
            block=>ps(igrid)
            call calc_x(igrid,xC,xCC)
            call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                           ixC^L,ixCC^L,.true.)
            write(qunit,fmt="(100(e14.6))") wCC_TMP(ixCC^S,iw)*normconv(iw)
         enddo
       enddo
     case default
       call mpistop('no such tecplot type')
   end select
   igonlevel=0
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (node(plevel_,igrid)/=level) cycle
      block=>ps(igrid)
      igonlevel=igonlevel+1
      call save_conntec(qunit,igrid,igonlevel)
   end do
end do
{^IFONED endif}

close(qunit)

end subroutine tecplot
!=============================================================================
subroutine save_conntec(qunit,igrid,igonlevel)

! this saves the basic line, quad and brick connectivity,
! as used by TECPLOT file outputs for unstructured grid

use mod_global_parameters

integer, intent(in) :: qunit, igrid, igonlevel

integer :: nx^D, nxC^D, ix^D
{^IFONED   integer, external:: nodenumbertec1D \}
{^IFTWOD   integer, external:: nodenumbertec2D \}
{^IFTHREED integer, external:: nodenumbertec3D \}
!-----------------------------------------------------------------------------
nx^D=ixMhi^D-ixMlo^D+1;
nxC^D=nx^D+1;

! connectivity list
{do ix^DB=1,nx^DB\}
   {^IFTHREED
   ! basic brick connectivity
   write(qunit,'(8(i7,1x))') &
      nodenumbertec3D(ix1,  ix2-1,ix3-1,nxC1,nxC2,nxC3,igonlevel,igrid),&
      nodenumbertec3D(ix1+1,ix2-1,ix3-1,nxC1,nxC2,nxC3,igonlevel,igrid),&
      nodenumbertec3D(ix1+1,ix2  ,ix3-1,nxC1,nxC2,nxC3,igonlevel,igrid),&
      nodenumbertec3D(ix1  ,ix2  ,ix3-1,nxC1,nxC2,nxC3,igonlevel,igrid),&
      nodenumbertec3D(ix1  ,ix2-1,ix3  ,nxC1,nxC2,nxC3,igonlevel,igrid),&
      nodenumbertec3D(ix1+1,ix2-1,ix3  ,nxC1,nxC2,nxC3,igonlevel,igrid),&
      nodenumbertec3D(ix1+1,ix2  ,ix3  ,nxC1,nxC2,nxC3,igonlevel,igrid),&
      nodenumbertec3D(ix1  ,ix2  ,ix3  ,nxC1,nxC2,nxC3,igonlevel,igrid)
   }
   {^IFTWOD
   ! basic quadrilateral connectivity
   write(qunit,'(4(i7,1x))') &
      nodenumbertec2D(ix1,  ix2-1,nxC1,nxC2,igonlevel,igrid),&
      nodenumbertec2D(ix1+1,ix2-1,nxC1,nxC2,igonlevel,igrid),&
      nodenumbertec2D(ix1+1,ix2  ,nxC1,nxC2,igonlevel,igrid),&
      nodenumbertec2D(ix1  ,ix2  ,nxC1,nxC2,igonlevel,igrid)
   }
   {^IFONED
   ! basic line connectivity
   write(qunit,'(2(i7,1x))') nodenumbertec1D(ix1,nxC1,igonlevel,igrid),&
                             nodenumbertec1D(ix1+1,nxC1,igonlevel,igrid)
   }
{end do\}

end subroutine save_conntec
!=============================================================================
integer function nodenumbertec1D(i1,nx1,ig,igrid)

integer, intent(in):: i1,nx1,ig,igrid
!-----------------------------------------------------------------------------
nodenumbertec1D=i1+(ig-1)*nx1

if(nodenumbertec1D>9999999)call mpistop("too large nodenumber")
end function nodenumbertec1D
!====================================================================================
integer function nodenumbertec2D(i1,i2,nx1,nx2,ig,igrid)

integer, intent(in):: i1,i2,nx1,nx2,ig,igrid
!-----------------------------------------------------------------------------
nodenumbertec2D=i1+i2*nx1+(ig-1)*nx1*nx2

if(nodenumbertec2D>9999999)call mpistop("too large nodenumber")
end function nodenumbertec2D
!====================================================================================
integer function nodenumbertec3D(i1,i2,i3,nx1,nx2,nx3,ig,igrid)

integer, intent(in):: i1,i2,i3,nx1,nx2,nx3,ig,igrid
!-----------------------------------------------------------------------------
nodenumbertec3D=i1+i2*nx1+i3*nx1*nx2+(ig-1)*nx1*nx2*nx3

if(nodenumbertec3D>9999999)call mpistop("too large nodenumber")
end function nodenumbertec3D
!=============================================================================
subroutine unstructuredvtk(qunit)

! output for vtu format to paraview
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors

use mod_global_parameters
use mod_calculate_xw

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP
double precision, dimension(0:nw+nwauxio)                   :: normconv
integer::               igrid,iigrid,level,igonlevel,icel,ixC^L,ixCC^L,iw
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx^D,nxC^D,nodesonlevel,elemsonlevel,nc,np,VTK_type,ix^D

character(len=80)::  filename
integer          :: filenr

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical :: fileopen
!-----------------------------------------------------------------------------

if(npe>1)then
  if(mype==0) PRINT *,'unstructuredvtk not parallel, use vtumpi'
  call mpistop('npe>1, unstructuredvtk')
end if

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
  ! generate filename 
  filenr=snapshotini
  if (autoconvert) filenr=snapshotnext
  write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".vtu"
  ! Open the file for the header part
  open(qunit,file=filename,status='unknown')
end if

call getheadernames(wnamei,xandwnamei,outfilehead)

! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'<UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(qunit,*) dble(dble(global_time)*time_convert_factor)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'

! number of cells, number of corner points, per grid.
nx^D=ixMhi^D-ixMlo^D+1;
nxC^D=nx^D+1;
nc={nx^D*}
np={nxC^D*}

! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
do level=levmin,levmax
  if (writelevel(level)) then
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (node(plevel_,igrid)/=level) cycle
      block=>ps(igrid)
      ! only output a grid when fully within clipped region selected
      ! by writespshift array
      if(({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
           *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
          <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})) then
        call calc_x(igrid,xC,xCC)
        call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                      ixC^L,ixCC^L,.true.)
        select case(convert_type)
        case('vtu')
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') &
             '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<PointData>'
          do iw=1,nw+nwauxio
            if(iw<=nw) then 
              if(.not.w_write(iw)) cycle
            end if
            write(qunit,'(a,a,a)')&
            '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') {(|}wC_TMP(ix^D,iw)*normconv(iw),{ix^D=ixCmin^D,ixCmax^D)}
            write(qunit,'(a)')'</DataArray>'
          end do
          write(qunit,'(a)')'</PointData>'
          write(qunit,'(a)')'<Points>'
          write(qunit,'(a)')'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          ! write cell corner coordinates in a backward dimensional loop, always 3D output
          {do ix^DB=ixCmin^DB,ixCmax^DB \}
             x_VTK(1:3)=zero;
             x_VTK(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0);
             write(qunit,'(3(1pe14.6))') x_VTK
          {end do \}
          write(qunit,'(a)')'</DataArray>'
          write(qunit,'(a)')'</Points>'
        case('vtuCC')
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') &
             '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<CellData>'
          do iw=1,nw+nwauxio
            if(iw<=nw) then 
              if(.not.w_write(iw)) cycle
            end if
            write(qunit,'(a,a,a)')&
            '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') {(|}wCC_TMP(ix^D,iw)*normconv(iw),{ix^D=ixCCmin^D,ixCCmax^D)}
            write(qunit,'(a)')'</DataArray>'
          end do
          write(qunit,'(a)')'</CellData>'
          write(qunit,'(a)')'<Points>'
          write(qunit,'(a)')'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          ! write cell corner coordinates in a backward dimensional loop, always 3D output
          {do ix^DB=ixCmin^DB,ixCmax^DB \}
             x_VTK(1:3)=zero;
             x_VTK(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0);
             write(qunit,'(3(1pe14.6))') x_VTK
          {end do \}
          write(qunit,'(a)')'</DataArray>'
          write(qunit,'(a)')'</Points>'
        end select

        write(qunit,'(a)')'<Cells>'
        ! connectivity part
        write(qunit,'(a)')'<DataArray type="Int32" Name="connectivity" format="ascii">'
        call save_connvtk(qunit,igrid)
        write(qunit,'(a)')'</DataArray>'

        ! offsets data array
        write(qunit,'(a)')'<DataArray type="Int32" Name="offsets" format="ascii">'
        do icel=1,nc
          write(qunit,'(i7)') icel*(2**^ND)
        end do
        write(qunit,'(a)')'</DataArray>'

        ! VTK cell type data array
        write(qunit,'(a)')'<DataArray type="Int32" Name="types" format="ascii">'
        ! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax
        {^IFONED VTK_type=3 \}
        {^IFTWOD VTK_type=8 \}
        {^IFTHREED VTK_type=11 \}
        do icel=1,nc
          write(qunit,'(i2)') VTK_type
        end do
        write(qunit,'(a)')'</DataArray>'

        write(qunit,'(a)')'</Cells>'

        write(qunit,'(a)')'</Piece>'
      end if
    end do
  end if
end do

write(qunit,'(a)')'</UnstructuredGrid>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)

end subroutine unstructuredvtk
!====================================================================================
subroutine unstructuredvtkB(qunit)

! output for vtu format to paraview, binary version output
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors
use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters
use mod_calculate_xw

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio):: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)  :: wCC_TMP
double precision :: normconv(0:nw+nwauxio)

integer, allocatable :: intstatus(:,:)
integer :: itag,ipe,igrid,level,icel,ixC^L,ixCC^L,Morton_no,Morton_length
integer :: nx^D,nxC^D,nc,np,VTK_type,ix^D,filenr
integer*8 :: offset

integer::  k,iw
integer::  length,lengthcc,length_coords,length_conn,length_offsets
character::  buf
character(len=80)::  filename
character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical ::   fileopen,cell_corner=.false.
logical, allocatable :: Morton_aim(:),Morton_aim_p(:)
!-----------------------------------------------------------------------------

normconv=one
Morton_length=Morton_stop(npe-1)-Morton_start(0)+1
allocate(Morton_aim(Morton_start(0):Morton_stop(npe-1)))
allocate(Morton_aim_p(Morton_start(0):Morton_stop(npe-1)))
Morton_aim=.false.
Morton_aim_p=.false.
do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  level=node(plevel_,igrid)
  ! we can clip parts of the grid away, select variables, levels etc.
  if(writelevel(level)) then
    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    if(({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
          *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
         <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})) then
      Morton_aim_p(Morton_no)=.true.
    end if
  end if
end do
call MPI_ALLREDUCE(Morton_aim_p,Morton_aim,Morton_length,MPI_LOGICAL,MPI_LOR,&
                         icomm,ierrmpi)
select case(convert_type)
 case('vtuB','vtuBmpi')
   cell_corner=.true.
 case('vtuBCC','vtuBCCmpi')
   cell_corner=.false.
end select
if (mype /= 0) then
  do Morton_no=Morton_start(mype),Morton_stop(mype)
    if(.not. Morton_aim(Morton_no)) cycle
    igrid=sfc_to_igrid(Morton_no)
    call calc_x(igrid,xC,xCC)
    call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                      ixC^L,ixCC^L,.true.)
    itag=Morton_no
    call MPI_SEND(xC_TMP,1,type_block_xc_io, 0,itag,icomm,ierrmpi)
    if(cell_corner) then
      call MPI_SEND(wC_TMP,1,type_block_wc_io, 0,itag,icomm,ierrmpi)
    else
      call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
    end if
  end do

else
  ! mype==0
  offset=0
  inquire(qunit,opened=fileopen)
  if(.not.fileopen)then
    ! generate filename 
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
    write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".vtu"
    ! Open the file for the header part
    open(qunit,file=filename,status='replace')
  end if
  call getheadernames(wnamei,xandwnamei,outfilehead)
  ! generate xml header
  write(qunit,'(a)')'<?xml version="1.0"?>'
  write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
  write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
  write(qunit,'(a)')'<UnstructuredGrid>'
  write(qunit,'(a)')'<FieldData>'
  write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                     'NumberOfTuples="1" format="ascii">'
  write(qunit,*) real(global_time*time_convert_factor)
  write(qunit,'(a)')'</DataArray>'
  write(qunit,'(a)')'</FieldData>'

  ! number of cells, number of corner points, per grid.
  nx^D=ixMhi^D-ixMlo^D+1;
  nxC^D=nx^D+1;
  nc={nx^D*}
  np={nxC^D*}
  length=np*size_real
  lengthcc=nc*size_real
  length_coords=3*length
  length_conn=2**^ND*size_int*nc
  length_offsets=nc*size_int

  ! Note: using the w_write, writelevel, writespshift
  do Morton_no=Morton_start(0),Morton_stop(0)
    if(.not. Morton_aim(Morton_no)) cycle
    if(cell_corner) then
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') &
         '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<PointData>'
      do iw=1,nw+nwauxio
        if(iw<=nw) then 
          if(.not.w_write(iw)) cycle
        endif

        write(qunit,'(a,a,a,i16,a)')&
            '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
            '" format="appended" offset="',offset,'">'
        write(qunit,'(a)')'</DataArray>'
        offset=offset+length+size_int
      end do
      write(qunit,'(a)')'</PointData>'
      write(qunit,'(a)')'<Points>'
      write(qunit,'(a,i16,a)') &
      '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      offset=offset+length_coords+size_int
      write(qunit,'(a)')'</Points>'
    else
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') &
         '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<CellData>'
      do iw=1,nw+nwauxio
        if(iw<=nw) then 
           if(.not.w_write(iw)) cycle
        end if
        write(qunit,'(a,a,a,i16,a)')&
            '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
            '" format="appended" offset="',offset,'">'
        write(qunit,'(a)')'</DataArray>'
        offset=offset+lengthcc+size_int
      end do
      write(qunit,'(a)')'</CellData>'
      write(qunit,'(a)')'<Points>'
      write(qunit,'(a,i16,a)') &
      '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      offset=offset+length_coords+size_int
      write(qunit,'(a)')'</Points>'
    end if
    write(qunit,'(a)')'<Cells>'

    ! connectivity part
    write(qunit,'(a,i16,a)')&
      '<DataArray type="Int32" Name="connectivity" format="appended" offset="',offset,'"/>'
    offset=offset+length_conn+size_int    

    ! offsets data array
    write(qunit,'(a,i16,a)') &
      '<DataArray type="Int32" Name="offsets" format="appended" offset="',offset,'"/>'
    offset=offset+length_offsets+size_int    

    ! VTK cell type data array
    write(qunit,'(a,i16,a)') &
      '<DataArray type="Int32" Name="types" format="appended" offset="',offset,'"/>' 
    offset=offset+size_int+nc*size_int

    write(qunit,'(a)')'</Cells>'

    write(qunit,'(a)')'</Piece>'
  end do
  ! write metadata communicated from other processors
  if(npe>1)then
    do ipe=1, npe-1
      do Morton_no=Morton_start(ipe),Morton_stop(ipe)
        if(.not. Morton_aim(Morton_no)) cycle
        if(cell_corner) then
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') &
             '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<PointData>'
          do iw=1,nw+nwauxio
            if(iw<=nw) then 
              if(.not.w_write(iw)) cycle
            end if
            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+length+size_int
          end do
          write(qunit,'(a)')'</PointData>'
          write(qunit,'(a)')'<Points>'
          write(qunit,'(a,i16,a)') &
          '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
          ! write cell corner coordinates in a backward dimensional loop, always 3D output
          offset=offset+length_coords+size_int
          write(qunit,'(a)')'</Points>'
        else
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') &
             '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<CellData>'
          do iw=1,nw+nwauxio
            if(iw<=nw) then 
              if(.not.w_write(iw)) cycle
            end if
            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+lengthcc+size_int
          end do
          write(qunit,'(a)')'</CellData>'
          write(qunit,'(a)')'<Points>'
          write(qunit,'(a,i16,a)') &
          '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
          ! write cell corner coordinates in a backward dimensional loop, always 3D output
          offset=offset+length_coords+size_int
          write(qunit,'(a)')'</Points>'
        end if
        write(qunit,'(a)')'<Cells>'
        ! connectivity part
        write(qunit,'(a,i16,a)')&
          '<DataArray type="Int32" Name="connectivity" format="appended" offset="',offset,'"/>'
        offset=offset+length_conn+size_int    

        ! offsets data array
        write(qunit,'(a,i16,a)') &
          '<DataArray type="Int32" Name="offsets" format="appended" offset="',offset,'"/>'
        offset=offset+length_offsets+size_int    

        ! VTK cell type data array
        write(qunit,'(a,i16,a)') &
          '<DataArray type="Int32" Name="types" format="appended" offset="',offset,'"/>' 
        offset=offset+size_int+nc*size_int

        write(qunit,'(a)')'</Cells>'

        write(qunit,'(a)')'</Piece>'
      end do
    end do
  end if

  write(qunit,'(a)')'</UnstructuredGrid>'
  write(qunit,'(a)')'<AppendedData encoding="raw">'
  close(qunit)
  open(qunit,file=filename,access='stream',form='unformatted',position='append')
  buf='_'
  write(qunit) TRIM(buf)

  do Morton_no=Morton_start(0),Morton_stop(0)
    if(.not. Morton_aim(Morton_no)) cycle
    igrid=sfc_to_igrid(Morton_no)
    call calc_x(igrid,xC,xCC)
    call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                   ixC^L,ixCC^L,.true.)
    do iw=1,nw+nwauxio
      if(iw<=nw) then 
        if(.not.w_write(iw)) cycle
      end if
      if(cell_corner) then
        write(qunit) length
        write(qunit) {(|}real(wC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCmin^D,ixCmax^D)}
      else
        write(qunit) lengthcc
        write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
      end if
    end do

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
    end do
  end do
  allocate(intstatus(MPI_STATUS_SIZE,1))
  if(npe>1)then
    ixCCmin^D=ixMlo^D; ixCCmax^D=ixMhi^D;
    ixCmin^D=ixMlo^D-1; ixCmax^D=ixMhi^D;
    do ipe=1, npe-1
      do Morton_no=Morton_start(ipe),Morton_stop(ipe)
        if(.not. Morton_aim(Morton_no)) cycle
        itag=Morton_no
        call MPI_RECV(xC_TMP,1,type_block_xc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
        if(cell_corner) then
          call MPI_RECV(wC_TMP,1,type_block_wc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
        else
          call MPI_RECV(wCC_TMP,1,type_block_wcc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
        end if
        do iw=1,nw+nwauxio
          if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
          end if
          if(cell_corner) then
            write(qunit) length
            write(qunit) {(|}real(wC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCmin^D,ixCmax^D)}
          else
            write(qunit) lengthcc
            write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
          end if
        end do

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
        end do
      end do
    end do
  end if
  close(qunit)
  open(qunit,file=filename,status='unknown',form='formatted',position='append')
  write(qunit,'(a)')'</AppendedData>'
  write(qunit,'(a)')'</VTKFile>'
  close(qunit)
  deallocate(intstatus)
end if

deallocate(Morton_aim,Morton_aim_p)
if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
end if

end subroutine unstructuredvtkB
!====================================================================================
subroutine unstructuredvtkB64(qunit)

! output for vtu format to paraview, binary version output
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors
use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters
use mod_calculate_xw

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio):: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)  :: wCC_TMP
double precision :: normconv(0:nw+nwauxio)

integer, allocatable :: intstatus(:,:)
integer :: itag,ipe,igrid,level,icel,ixC^L,ixCC^L,Morton_no,Morton_length
integer :: nx^D,nxC^D,nc,np,VTK_type,ix^D,filenr
integer*8 :: offset

integer::  k,iw
integer::  length,lengthcc,length_coords,length_conn,length_offsets
character::  buf
character(len=80)::  filename
character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical ::   fileopen,cell_corner=.false.
logical, allocatable :: Morton_aim(:),Morton_aim_p(:)
!-----------------------------------------------------------------------------

normconv=one
Morton_length=Morton_stop(npe-1)-Morton_start(0)+1
allocate(Morton_aim(Morton_start(0):Morton_stop(npe-1)))
allocate(Morton_aim_p(Morton_start(0):Morton_stop(npe-1)))
Morton_aim=.false.
Morton_aim_p=.false.
do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  level=node(plevel_,igrid)
  ! we can clip parts of the grid away, select variables, levels etc.
  if(writelevel(level)) then
    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    if(({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
          *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
         <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})) then
      Morton_aim_p(Morton_no)=.true.
    end if
  end if
end do
call MPI_ALLREDUCE(Morton_aim_p,Morton_aim,Morton_length,MPI_LOGICAL,MPI_LOR,&
                         icomm,ierrmpi)
select case(convert_type)
 case('vtuB64','vtuBmpi64')
   cell_corner=.true.
 case('vtuBCC64','vtuBCCmpi64')
   cell_corner=.false.
end select
if (mype /= 0) then
  do Morton_no=Morton_start(mype),Morton_stop(mype)
    if(.not. Morton_aim(Morton_no)) cycle
    igrid=sfc_to_igrid(Morton_no)
    call calc_x(igrid,xC,xCC)
    call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                      ixC^L,ixCC^L,.true.)
    itag=Morton_no
    call MPI_SEND(xC_TMP,1,type_block_xc_io, 0,itag,icomm,ierrmpi)
    if(cell_corner) then
      call MPI_SEND(wC_TMP,1,type_block_wc_io, 0,itag,icomm,ierrmpi)
    else
      call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
    end if
  end do

else
  ! mype==0
  offset=0
  inquire(qunit,opened=fileopen)
  if(.not.fileopen)then
    ! generate filename 
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
    write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".vtu"
    ! Open the file for the header part
    open(qunit,file=filename,status='replace')
  end if
  call getheadernames(wnamei,xandwnamei,outfilehead)
  ! generate xml header
  write(qunit,'(a)')'<?xml version="1.0"?>'
  write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
  write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
  write(qunit,'(a)')'<UnstructuredGrid>'
  write(qunit,'(a)')'<FieldData>'
  write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                     'NumberOfTuples="1" format="ascii">'
  write(qunit,*) real(global_time*time_convert_factor)
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
  do Morton_no=Morton_start(0),Morton_stop(0)
    if(.not. Morton_aim(Morton_no)) cycle
    if(cell_corner) then
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') &
         '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<PointData>'
      do iw=1,nw+nwauxio
        if(iw<=nw) then 
          if(.not.w_write(iw)) cycle
        endif

        write(qunit,'(a,a,a,i16,a)')&
            '<DataArray type="Float64" Name="',TRIM(wnamei(iw)), &
            '" format="appended" offset="',offset,'">'
        write(qunit,'(a)')'</DataArray>'
        offset=offset+length+size_int
      end do
      write(qunit,'(a)')'</PointData>'
      write(qunit,'(a)')'<Points>'
      write(qunit,'(a,i16,a)') &
      '<DataArray type="Float64" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      offset=offset+length_coords+size_int
      write(qunit,'(a)')'</Points>'
    else
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') &
         '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<CellData>'
      do iw=1,nw+nwauxio
        if(iw<=nw) then 
           if(.not.w_write(iw)) cycle
        end if
        write(qunit,'(a,a,a,i16,a)')&
            '<DataArray type="Float64" Name="',TRIM(wnamei(iw)), &
            '" format="appended" offset="',offset,'">'
        write(qunit,'(a)')'</DataArray>'
        offset=offset+lengthcc+size_int
      end do
      write(qunit,'(a)')'</CellData>'
      write(qunit,'(a)')'<Points>'
      write(qunit,'(a,i16,a)') &
      '<DataArray type="Float64" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      offset=offset+length_coords+size_int
      write(qunit,'(a)')'</Points>'
    end if
    write(qunit,'(a)')'<Cells>'

    ! connectivity part
    write(qunit,'(a,i16,a)')&
      '<DataArray type="Int32" Name="connectivity" format="appended" offset="',offset,'"/>'
    offset=offset+length_conn+size_int    

    ! offsets data array
    write(qunit,'(a,i16,a)') &
      '<DataArray type="Int32" Name="offsets" format="appended" offset="',offset,'"/>'
    offset=offset+length_offsets+size_int    

    ! VTK cell type data array
    write(qunit,'(a,i16,a)') &
      '<DataArray type="Int32" Name="types" format="appended" offset="',offset,'"/>' 
    offset=offset+size_int+nc*size_int

    write(qunit,'(a)')'</Cells>'

    write(qunit,'(a)')'</Piece>'
  end do
  ! write metadata communicated from other processors
  if(npe>1)then
    do ipe=1, npe-1
      do Morton_no=Morton_start(ipe),Morton_stop(ipe)
        if(.not. Morton_aim(Morton_no)) cycle
        if(cell_corner) then
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') &
             '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<PointData>'
          do iw=1,nw+nwauxio
            if(iw<=nw) then 
              if(.not.w_write(iw)) cycle
            end if
            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float64" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+length+size_int
          end do
          write(qunit,'(a)')'</PointData>'
          write(qunit,'(a)')'<Points>'
          write(qunit,'(a,i16,a)') &
          '<DataArray type="Float64" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
          ! write cell corner coordinates in a backward dimensional loop, always 3D output
          offset=offset+length_coords+size_int
          write(qunit,'(a)')'</Points>'
        else
          ! we write out every grid as one VTK PIECE
          write(qunit,'(a,i7,a,i7,a)') &
             '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<CellData>'
          do iw=1,nw+nwauxio
            if(iw<=nw) then 
              if(.not.w_write(iw)) cycle
            end if
            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float64" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+lengthcc+size_int
          end do
          write(qunit,'(a)')'</CellData>'
          write(qunit,'(a)')'<Points>'
          write(qunit,'(a,i16,a)') &
          '<DataArray type="Float64" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
          ! write cell corner coordinates in a backward dimensional loop, always 3D output
          offset=offset+length_coords+size_int
          write(qunit,'(a)')'</Points>'
        end if
        write(qunit,'(a)')'<Cells>'
        ! connectivity part
        write(qunit,'(a,i16,a)')&
          '<DataArray type="Int32" Name="connectivity" format="appended" offset="',offset,'"/>'
        offset=offset+length_conn+size_int    

        ! offsets data array
        write(qunit,'(a,i16,a)') &
          '<DataArray type="Int32" Name="offsets" format="appended" offset="',offset,'"/>'
        offset=offset+length_offsets+size_int    

        ! VTK cell type data array
        write(qunit,'(a,i16,a)') &
          '<DataArray type="Int32" Name="types" format="appended" offset="',offset,'"/>' 
        offset=offset+size_int+nc*size_int

        write(qunit,'(a)')'</Cells>'

        write(qunit,'(a)')'</Piece>'
      end do
    end do
  end if

  write(qunit,'(a)')'</UnstructuredGrid>'
  write(qunit,'(a)')'<AppendedData encoding="raw">'
  close(qunit)
  open(qunit,file=filename,access='stream',form='unformatted',position='append')
  buf='_'
  write(qunit) TRIM(buf)

  do Morton_no=Morton_start(0),Morton_stop(0)
    if(.not. Morton_aim(Morton_no)) cycle
    igrid=sfc_to_igrid(Morton_no)
    call calc_x(igrid,xC,xCC)
    call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                   ixC^L,ixCC^L,.true.)
    do iw=1,nw+nwauxio
      if(iw<=nw) then 
        if(.not.w_write(iw)) cycle
      end if
      if(cell_corner) then
        write(qunit) length
        write(qunit) {(|}wC_TMP(ix^D,iw)*normconv(iw),{ix^D=ixCmin^D,ixCmax^D)}
      else
        write(qunit) lengthcc
        write(qunit) {(|}wCC_TMP(ix^D,iw)*normconv(iw),{ix^D=ixCCmin^D,ixCCmax^D)}
      end if
    end do

    write(qunit) length_coords
    {do ix^DB=ixCmin^DB,ixCmax^DB \}
      x_VTK(1:3)=zero;
      x_VTK(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0);
      do k=1,3
        write(qunit) x_VTK(k)
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
    end do
  end do
  allocate(intstatus(MPI_STATUS_SIZE,1))
  if(npe>1)then
    ixCCmin^D=ixMlo^D; ixCCmax^D=ixMhi^D;
    ixCmin^D=ixMlo^D-1; ixCmax^D=ixMhi^D;
    do ipe=1, npe-1
      do Morton_no=Morton_start(ipe),Morton_stop(ipe)
        if(.not. Morton_aim(Morton_no)) cycle
        itag=Morton_no
        call MPI_RECV(xC_TMP,1,type_block_xc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
        if(cell_corner) then
          call MPI_RECV(wC_TMP,1,type_block_wc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
        else
          call MPI_RECV(wCC_TMP,1,type_block_wcc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
        end if
        do iw=1,nw+nwauxio
          if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
          end if
          if(cell_corner) then
            write(qunit) length
            write(qunit) {(|}wC_TMP(ix^D,iw)*normconv(iw),{ix^D=ixCmin^D,ixCmax^D)}
          else
            write(qunit) lengthcc
            write(qunit) {(|}wCC_TMP(ix^D,iw)*normconv(iw),{ix^D=ixCCmin^D,ixCCmax^D)}
          end if
        end do

        write(qunit) length_coords
        {do ix^DB=ixCmin^DB,ixCmax^DB \}
          x_VTK(1:3)=zero;
          x_VTK(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0);
          do k=1,3
           write(qunit) x_VTK(k)
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
        end do
      end do
    end do
  end if
  close(qunit)
  open(qunit,file=filename,status='unknown',form='formatted',position='append')
  write(qunit,'(a)')'</AppendedData>'
  write(qunit,'(a)')'</VTKFile>'
  close(qunit)
  deallocate(intstatus)
end if

deallocate(Morton_aim,Morton_aim_p)
if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
end if

end subroutine unstructuredvtkB64
!====================================================================================
subroutine save_connvtk(qunit,igrid)

! this saves the basic line, pixel and voxel connectivity,
! as used by VTK file outputs for unstructured grid

use mod_global_parameters

integer, intent(in) :: qunit, igrid

integer :: nx^D, nxC^D, ix^D
!-----------------------------------------------------------------------------
nx^D=ixMhi^D-ixMlo^D+1;
nxC^D=nx^D+1;

{do ix^DB=1,nx^DB\}
        {^IFONED write(qunit,'(2(i7,1x))')ix1-1,ix1 \}
        {^IFTWOD
        write(qunit,'(4(i7,1x))')(ix2-1)*nxC1+ix1-1, &
               (ix2-1)*nxC1+ix1,ix2*nxC1+ix1-1,ix2*nxC1+ix1
        \}
        {^IFTHREED
        write(qunit,'(8(i7,1x))')&
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

end subroutine save_connvtk
!=============================================================================
subroutine ImageDataVtk_mpi(qunit)

! output for vti format to paraview, non-binary version output
! parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors
! allows skipping of w_write selected variables

! implementation such that length of ASCII output is identical when 
! run on 1 versus multiple CPUs (however, the order of the vtu pieces can differ)

use mod_forest, only: Morton_start, Morton_stop, tree_node_ptr, igrid_to_node, sfc_to_igrid
use mod_global_parameters
use mod_calculate_xw

integer, intent(in) ::    qunit

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP,xC_TMP_recv
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP,xCC_TMP_recv
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP,wC_TMP_recv
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP,wCC_TMP_recv
double precision, dimension(0:nw+nwauxio)                   :: normconv
integer::               igrid,iigrid,level,ixC^L,ixCC^L
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx^D

character(len=80)::  filename
integer ::           filenr

integer, allocatable :: intstatus(:,:)
logical, allocatable :: Morton_aim(:),Morton_aim_p(:)

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical :: fileopen
integer :: itag,ipe,Morton_no,Morton_length
integer :: ixrvC^L, ixrvCC^L, siz_ind, ind_send(5*^ND), ind_recv(5*^ND)
double precision    :: origin(1:3), spacing(1:3)
integer :: wholeExtent(1:6), ig^D
type(tree_node_ptr) :: tree
!-----------------------------------------------------------------------------
if(levmin/=levmax) call mpistop('ImageData can only be used when levmin=levmax')

normconv(0) = length_convert_factor
normconv(1:nw) = w_convert_factor
siz_ind=5*^ND
Morton_length=Morton_stop(npe-1)-Morton_start(0)+1
allocate(Morton_aim(Morton_start(0):Morton_stop(npe-1)))
allocate(Morton_aim_p(Morton_start(0):Morton_stop(npe-1)))
Morton_aim=.false.
Morton_aim_p=.false.
do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  level=node(plevel_,igrid)
  ! we can clip parts of the grid away, select variables, levels etc.
  if(writelevel(level)) then
   ! only output a grid when fully within clipped region selected
   ! by writespshift array
   if(({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
         *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
        <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})) then
     Morton_aim_p(Morton_no)=.true.
   end if
  end if
end do
call MPI_ALLREDUCE(Morton_aim_p,Morton_aim,Morton_length,MPI_LOGICAL,MPI_LOR,&
                         icomm,ierrmpi)


if (mype /= 0) then
 do Morton_no=Morton_start(mype),Morton_stop(mype)
   if(.not. Morton_aim(Morton_no)) cycle
   igrid=sfc_to_igrid(Morton_no)
   call calc_x(igrid,xC,xCC)
   call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                     ixC^L,ixCC^L,.true.)
   tree%node => igrid_to_node(igrid, mype)%node
   {^D& ig^D = tree%node%ig^D; }
   itag=Morton_no
   ind_send=(/ ixC^L,ixCC^L, ig^D /)
   call MPI_SEND(ind_send,siz_ind,MPI_INTEGER, 0,itag,icomm,ierrmpi)
   call MPI_SEND(wC_TMP,1,type_block_wc_io, 0,itag,icomm,ierrmpi)
   call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
 end do

else

 inquire(qunit,opened=fileopen)
 if(.not.fileopen)then
    ! generate filename 
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
    write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".vti"
   ! Open the file for the header part
   open(qunit,file=filename,status='unknown',form='formatted')
 endif

call getheadernames(wnamei,xandwnamei,outfilehead)

! number of cells per grid.
nx^D=ixMhi^D-ixMlo^D+1;

origin      = 0
{^D& origin(^D) = xprobmin^D*normconv(0); }
spacing     = zero
{^D&spacing(^D) = dxlevel(^D)*normconv(0); }

wholeExtent = 0
! if we use writespshift, the whole extent has to be calculated:
{^D&wholeExtent(^D*2-1) = nx^D * ceiling(((xprobmax^D-xprobmin^D)*writespshift(^D,1)) &
     /(nx^D*dxlevel(^D))) \}
{^D&wholeExtent(^D*2)   = nx^D * floor(((xprobmax^D-xprobmin^D)*(1.0d0-writespshift(^D,2))) &
     /(nx^D*dxlevel(^D))) \}

! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="ImageData"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a,3(1pe14.6),a,6(i10),a,3(1pe14.6),a)')'  <ImageData Origin="',&
     origin,'" WholeExtent="',wholeExtent,'" Spacing="',spacing,'">'
 write(qunit,'(a)')'<FieldData>'
 write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                    'NumberOfTuples="1" format="ascii">'
 write(qunit,*) real(global_time*time_convert_factor)
 write(qunit,'(a)')'</DataArray>'
 write(qunit,'(a)')'</FieldData>'

! write the data from proc 0
do Morton_no=Morton_start(0),Morton_stop(0)
   if(.not. Morton_aim(Morton_no)) cycle
   igrid=sfc_to_igrid(Morton_no)
   tree%node => igrid_to_node(igrid, 0)%node
   {^D& ig^D = tree%node%ig^D; }
   call calc_x(igrid,xC,xCC)
   call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
        ixC^L,ixCC^L,.true.)
   call write_vti(qunit,ixG^LL,ixC^L,ixCC^L,ig^D,&
        nx^D,normconv,wnamei,wC_TMP,wCC_TMP)   
end do

if(npe>1)then
   allocate(intstatus(MPI_STATUS_SIZE,1))
   do ipe=1, npe-1
      do Morton_no=Morton_start(ipe),Morton_stop(ipe)
         if(.not. Morton_aim(Morton_no)) cycle
         itag=Morton_no
         call MPI_RECV(ind_recv,siz_ind, MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
         ixrvCmin^D=ind_recv(^D);ixrvCmax^D=ind_recv(^ND+^D);
         ixrvCCmin^D=ind_recv(2*^ND+^D);ixrvCCmax^D=ind_recv(3*^ND+^D);
         ig^D=ind_recv(4*^ND+^D);
         call MPI_RECV(wC_TMP,1,type_block_wc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
         call MPI_RECV(wCC_TMP,1,type_block_wcc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
         call write_vti(qunit,ixG^LL,ixrvC^L,ixrvCC^L,ig^D,&
              nx^D,normconv,wnamei,wC_TMP,wCC_TMP)   
      end do
   end do
end if

write(qunit,'(a)')'</ImageData>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)
if(npe>1) deallocate(intstatus)
endif

deallocate(Morton_aim,Morton_aim_p)
if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
endif

end subroutine ImageDataVtk_mpi
!============================================================================
subroutine punstructuredvtk_mpi(qunit)

! Write one pvtu and vtu files for each processor
! Otherwise like unstructuredvtk_mpi

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters
use mod_calculate_xw

integer, intent(in) ::    qunit
!
double precision, dimension(0:nw+nwauxio)                   :: normconv
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim)         :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)           :: xCC_TMP
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim)         :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)           :: xCC
double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP
character(len=name_len)   :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead
integer             :: nx^D,nxC^D,nc,np, igrid,ixC^L,ixCC^L,level,Morton_no
character(len=80)   :: pfilename
integer             :: filenr
logical             :: fileopen,conv_grid
!----------------------------------------------------------------------------

! Write pvtu-file:
if (mype==0) then
   call write_pvtu(qunit)
endif
! Now write the Source files:

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
   ! generate filename 
   filenr=snapshotini
   if (autoconvert) filenr=snapshotnext
   ! Open the file for the header part
   write(pfilename,'(a,i4.4,a,i4.4,a)') TRIM(base_filename),filenr,"p",mype,".vtu"
   open(qunit,file=pfilename,status='unknown',form='formatted')
endif
! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'  <UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(qunit,*) real(global_time*time_convert_factor)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'

call getheadernames(wnamei,xandwnamei,outfilehead)

! number of cells, number of corner points, per grid.
nx^D=ixMhi^D-ixMlo^D+1;
nxC^D=nx^D+1;
nc={nx^D*}
np={nxC^D*}

! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
do level=levmin,levmax
   if (.not.writelevel(level)) cycle
   do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    if (node(plevel_,igrid)/=level) cycle

    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    conv_grid=({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
          *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
         <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})
    if (.not.conv_grid) cycle

    call calc_x(igrid,xC,xCC)
    call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
         ixC^L,ixCC^L,.true.)

    call write_vtk(qunit,ixG^LL,ixC^L,ixCC^L,igrid,nc,np,nx^D,nxC^D,&
         normconv,wnamei,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP)
   end do ! Morton_no loop
end do ! level loop

 write(qunit,'(a)')'  </UnstructuredGrid>'
 write(qunit,'(a)')'</VTKFile>'
 close(qunit)

if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
endif
end subroutine punstructuredvtk_mpi
!============================================================================
subroutine unstructuredvtk_mpi(qunit)

! output for vtu format to paraview, non-binary version output
! parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors
! allows skipping of w_write selected variables

! implementation such that length of ASCII output is identical when 
! run on 1 versus multiple CPUs (however, the order of the vtu pieces can differ)

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters
use mod_calculate_xw

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP,xC_TMP_recv
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP,xCC_TMP_recv
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP,wC_TMP_recv
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP,wCC_TMP_recv
double precision, dimension(0:nw+nwauxio)                   :: normconv
integer::               igrid,iigrid,level,ixC^L,ixCC^L
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx^D,nxC^D,nodesonlevel,elemsonlevel,nc,np,ix^D

character(len=80)::  filename
integer ::           filenr

integer, allocatable :: intstatus(:,:)

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical :: fileopen,conv_grid,cond_grid_recv
integer :: itag,ipe,Morton_no,siz_ind
integer :: ind_send(4*^ND),ind_recv(4*^ND)
integer :: levmin_recv,levmax_recv,level_recv,igrid_recv,ixrvC^L,ixrvCC^L
!-----------------------------------------------------------------------------
if (mype==0) then
 inquire(qunit,opened=fileopen)
 if(.not.fileopen)then
    ! generate filename 
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
    write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".vtu"
   ! Open the file for the header part
   open(qunit,file=filename,status='unknown',form='formatted')
 endif
 ! generate xml header
 write(qunit,'(a)')'<?xml version="1.0"?>'
 write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
 write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
 write(qunit,'(a)')'<UnstructuredGrid>'
 write(qunit,'(a)')'<FieldData>'
 write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                    'NumberOfTuples="1" format="ascii">'
 write(qunit,*) real(global_time*time_convert_factor)
 write(qunit,'(a)')'</DataArray>'
 write(qunit,'(a)')'</FieldData>'
end if

call getheadernames(wnamei,xandwnamei,outfilehead)
! number of cells, number of corner points, per grid.
nx^D=ixMhi^D-ixMlo^D+1;
nxC^D=nx^D+1;
nc={nx^D*}
np={nxC^D*}

! all slave processors send their minmal/maximal levels
if  (mype/=0) then
 if (Morton_stop(mype)==0) call mpistop("nultag")
 itag=1000*Morton_stop(mype)
 !print *,'ype,itag for levmin=',mype,itag,levmin
 call MPI_SEND(levmin,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
 itag=2000*Morton_stop(mype)
 !print *,'mype,itag for levmax=',mype,itag,levmax
 call MPI_SEND(levmax,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
end if


! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
do level=levmin,levmax
   if (.not.writelevel(level)) cycle
   do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    if (mype/=0)then
      itag=Morton_no
      call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
      itag=igrid
      call MPI_SEND(node(plevel_,igrid),1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
    end if
    if (node(plevel_,igrid)/=level) cycle

    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    conv_grid=({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
          *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
         <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})
    if (mype/=0)then
      call MPI_SEND(conv_grid,1,MPI_LOGICAL,0,itag,icomm,ierrmpi)
    end if
    if (.not.conv_grid) cycle

    call calc_x(igrid,xC,xCC)
    call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                     ixC^L,ixCC^L,.true.)

    if (mype/=0) then
       itag=Morton_no
       ind_send=(/ ixC^L,ixCC^L /)
       siz_ind=4*^ND
       call MPI_SEND(ind_send,siz_ind,MPI_INTEGER, 0,itag,icomm,ierrmpi)
       call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION, 0,itag,icomm,ierrmpi)

       call MPI_SEND(wC_TMP,1,type_block_wc_io, 0,itag,icomm,ierrmpi)
       call MPI_SEND(xC_TMP,1,type_block_xc_io, 0,itag,icomm,ierrmpi)
       itag=igrid
       call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
       call MPI_SEND(xCC_TMP,1,type_block_xcc_io, 0,itag,icomm,ierrmpi)
    else
       call write_vtk(qunit,ixG^LL,ixC^L,ixCC^L,igrid,nc,np,nx^D,nxC^D,&
                          normconv,wnamei,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP)
    end if
   end do ! Morton_no loop
end do ! level loop


if (mype==0) then
 allocate(intstatus(MPI_STATUS_SIZE,1))
 if(npe>1)then
  do ipe=1,npe-1
   itag=1000*Morton_stop(ipe)
   call MPI_RECV(levmin_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
   !!print *,'mype RECEIVES,itag for levmin=',mype,itag,levmin_recv
   itag=2000*Morton_stop(ipe)
   call MPI_RECV(levmax_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
   !!print *,'mype RECEIVES itag for levmax=',mype,itag,levmax_recv
   do level=levmin_recv,levmax_recv
    if (.not.writelevel(level)) cycle
    do  Morton_no=Morton_start(ipe),Morton_stop(ipe)
     itag=Morton_no
     call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
     itag=igrid_recv
     call MPI_RECV(level_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
     if (level_recv/=level) cycle

     call MPI_RECV(cond_grid_recv,1,MPI_LOGICAL, ipe,itag,icomm,intstatus(:,1),ierrmpi)
     if(.not.cond_grid_recv)cycle

     itag=Morton_no
     siz_ind=4*^ND
     call MPI_RECV(ind_recv,siz_ind, MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
     ixrvCmin^D=ind_recv(^D);ixrvCmax^D=ind_recv(^ND+^D);
     ixrvCCmin^D=ind_recv(2*^ND+^D);ixrvCCmax^D=ind_recv(3*^ND+^D);
     call MPI_RECV(normconv,nw+nwauxio+1, MPI_DOUBLE_PRECISION,ipe,itag&
                  ,icomm,intstatus(:,1),ierrmpi)

     call MPI_RECV(wC_TMP_recv,1,type_block_wc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
     call MPI_RECV(xC_TMP_recv,1,type_block_xc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)

     itag=igrid_recv
     call MPI_RECV(wCC_TMP_recv,1,type_block_wcc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
     call MPI_RECV(xCC_TMP_recv,1,type_block_xcc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
     call write_vtk(qunit,ixG^LL,ixrvC^L,ixrvCC^L,igrid_recv,&
                    nc,np,nx^D,nxC^D,normconv,wnamei,&
                    xC_TMP_recv,xCC_TMP_recv,wC_TMP_recv,wCC_TMP_recv)
    enddo ! Morton_no loop
   enddo ! level loop
  enddo ! processor loop
 endif ! multiple processors
 write(qunit,'(a)')'</UnstructuredGrid>'
 write(qunit,'(a)')'</VTKFile>'
 close(qunit)
endif

if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
  if(mype==0)deallocate(intstatus)
endif

end subroutine unstructuredvtk_mpi
!============================================================================
subroutine write_vtk(qunit,ixI^L,ixC^L,ixCC^L,igrid,nc,np,nx^D,nxC^D,&
                     normconv,wnamei,xC,xCC,wC,wCC)

use mod_global_parameters

integer, intent(in) :: qunit
integer, intent(in) :: ixI^L,ixC^L,ixCC^L
integer, intent(in) :: igrid,nc,np,nx^D,nxC^D
double precision, intent(in) :: normconv(0:nw+nwauxio) 
character(len=name_len), intent(in)::  wnamei(1:nw+nwauxio)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC

double precision ::  x_VTK(1:3)
integer :: iw,ix^D,icel,VTK_type
!----------------------------------------------------------------------------

select case(convert_type)
    case('vtumpi','pvtumpi')
         ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') &
            '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<PointData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
         endif

            write(qunit,'(a,a,a)')&
          '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') {(|}wC(ix^D,iw)*normconv(iw),{ix^D=ixCmin^D,ixCmax^D)}
            write(qunit,'(a)')'</DataArray>'
      enddo
      write(qunit,'(a)')'</PointData>'

      write(qunit,'(a)')'<Points>'
      write(qunit,'(a)')'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
      {do ix^DB=ixCmin^DB,ixCmax^DB \}
            x_VTK(1:3)=zero;
            x_VTK(1:ndim)=xC(ix^D,1:ndim)*normconv(0);
            write(qunit,'(3(1pe14.6))') x_VTK
      {end do \}
      write(qunit,'(a)')'</DataArray>'
      write(qunit,'(a)')'</Points>'

    case('vtuCCmpi','pvtuCCmpi')
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') &
         '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<CellData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
         endif
            write(qunit,'(a,a,a)')&
          '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') {(|}wCC(ix^D,iw)*normconv(iw),{ix^D=ixCCmin^D,ixCCmax^D)}
            write(qunit,'(a)')'</DataArray>'
      enddo
      write(qunit,'(a)')'</CellData>'

      write(qunit,'(a)')'<Points>'
      write(qunit,'(a)')'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      {do ix^DB=ixCmin^DB,ixCmax^DB \}
            x_VTK(1:3)=zero;
            x_VTK(1:ndim)=xC(ix^D,1:ndim)*normconv(0);
            write(qunit,'(3(1pe14.6))') x_VTK
      {end do \}
      write(qunit,'(a)')'</DataArray>'
      write(qunit,'(a)')'</Points>'
end select

write(qunit,'(a)')'<Cells>'

! connectivity part
write(qunit,'(a)')'<DataArray type="Int32" Name="connectivity" format="ascii">'
call save_connvtk(qunit,igrid)
write(qunit,'(a)')'</DataArray>'

! offsets data array
write(qunit,'(a)')'<DataArray type="Int32" Name="offsets" format="ascii">'
do icel=1,nc
    write(qunit,'(i7)') icel*(2**^ND)
end do
write(qunit,'(a)')'</DataArray>'

! VTK cell type data array
write(qunit,'(a)')'<DataArray type="Int32" Name="types" format="ascii">'
! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax
{^IFONED VTK_type=3 \}
{^IFTWOD VTK_type=8 \}
{^IFTHREED VTK_type=11 \}
do icel=1,nc
   write(qunit,'(i2)') VTK_type
enddo
write(qunit,'(a)')'</DataArray>'

write(qunit,'(a)')'</Cells>'

write(qunit,'(a)')'</Piece>'

end subroutine write_vtk
!============================================================================
subroutine write_vti(qunit,ixI^L,ixC^L,ixCC^L,ig^D,nx^D,&
                     normconv,wnamei,wC,wCC)
use mod_global_parameters

integer, intent(in) :: qunit
integer, intent(in) :: ixI^L,ixC^L,ixCC^L
integer, intent(in) :: ig^D,nx^D
double precision, intent(in) :: normconv(0:nw+nwauxio) 
character(len=name_len), intent(in)::  wnamei(1:nw+nwauxio)

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC

integer :: iw,ix^D
integer :: extent(1:6)
!----------------------------------------------------------------------------

extent = 0
{^D& extent(^D*2-1) = (ig^D-1) * nx^D; }
{^D& extent(^D*2)   = (ig^D)   * nx^D; }


select case(convert_type)
    case('vtimpi','pvtimpi')
         ! we write out every grid as one VTK PIECE
      write(qunit,'(a,6(i10),a)') &
            '<Piece Extent="',extent,'">'
      write(qunit,'(a)')'<PointData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
         endif

            write(qunit,'(a,a,a)')&
          '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe20.12))') {(|}wC(ix^D,iw)*normconv(iw),{ix^D=ixCmin^D,ixCmax^D)}
            write(qunit,'(a)')'</DataArray>'
      enddo
      write(qunit,'(a)')'</PointData>'

    case('vtiCCmpi','pvtiCCmpi')
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,6(i10),a)') &
            '<Piece Extent="',extent,'">'
      write(qunit,'(a)')'<CellData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
         endif
            write(qunit,'(a,a,a)')&
          '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe20.12))') {(|}wCC(ix^D,iw)*normconv(iw),{ix^D=ixCCmin^D,ixCCmax^D)}
            write(qunit,'(a)')'</DataArray>'
      enddo
      write(qunit,'(a)')'</CellData>'
end select

write(qunit,'(a)')'</Piece>'

end subroutine write_vti
!=============================================================================
subroutine write_pvtu(qunit)

use mod_global_parameters
use mod_calculate_xw

integer, intent(in) :: qunit

character(len=name_len)   :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio),outtype
character(len=1024) :: outfilehead
character(len=80)   :: filename,pfilename
integer             :: filenr,iw,ipe,iscalars
logical             :: fileopen

select case(convert_type)
case('pvtumpi','pvtuBmpi')
   outtype="PPointData"
case('pvtuCCmpi','pvtuBCCmpi')
   outtype="PCellData"
end select
inquire(qunit,opened=fileopen)
if(.not.fileopen)then
   ! generate filename 
   filenr=snapshotini
   if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".pvtu"
   ! Open the file
   open(qunit,file=filename,status='unknown',form='formatted')
endif

call getheadernames(wnamei,xandwnamei,outfilehead)

! Get the default selection:
iscalars=1
do iw=nw,1, -1
   if (w_write(iw)) iscalars=iw
end do


! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="PUnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'  <PUnstructuredGrid GhostLevel="0">'
! Either celldata or pointdata:
write(qunit,'(a,a,a,a,a)')&
     '    <',TRIM(outtype),' Scalars="',TRIM(wnamei(iscalars))//'">'
do iw=1,nw
   if(.not.w_write(iw))cycle
   write(qunit,'(a,a,a)')&
        '      <PDataArray type="Float32" Name="',TRIM(wnamei(iw)),'"/>'
end do
do iw=nw+1,nw+nwauxio
   write(qunit,'(a,a,a)')&
        '      <PDataArray type="Float32" Name="',TRIM(wnamei(iw)),'"/>'
end do
write(qunit,'(a,a,a)')'    </',TRIM(outtype),'>'

write(qunit,'(a)')'    <PPoints>'
write(qunit,'(a)')'      <PDataArray type="Float32" NumberOfComponents="3"/>'
write(qunit,'(a)')'    </PPoints>'

do ipe=0,npe-1
   write(pfilename,'(a,i4.4,a,i4.4,a)') TRIM(base_filename(&
        INDEX (base_filename, '/', BACK = .TRUE.)+1:&
        LEN(base_filename))),filenr,"p",&
        ipe,".vtu"
   write(qunit,'(a,a,a)')'    <Piece Source="',TRIM(pfilename),'"/>'
end do
write(qunit,'(a)')'  </PUnstructuredGrid>'
write(qunit,'(a)')'</VTKFile>'

close(qunit)

end subroutine write_pvtu
!=============================================================================
subroutine tecplot_mpi(qunit)

! output for tecplot (ASCII format)
! parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors

! the current implementation is such that tecplotmpi and tecplotCCmpi will 
! create different length output ASCII files when used on 1 versus multiple CPUs
! in fact, on 1 CPU, there will be as many zones as there are levels
! on multiple CPUs, there will be a number of zones up to the number of
! levels times the number of CPUs (can be less, when some level not on a CPU)

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters
use mod_calculate_xw

integer, intent(in) :: qunit

integer::               igrid,iigrid,level,igonlevel,iw,idim,ix^D
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx^D,nxC^D,nodesonlevel,elemsonlevel,ixC^L,ixCC^L
integer :: nodesonlevelmype,elemsonlevelmype

integer ::              nodes, elems

integer, allocatable :: intstatus(:,:)

double precision :: x_TEC(ndim), w_TEC(nw+nwauxio)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP,xC_TMP_recv
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP,xCC_TMP_recv
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP,wC_TMP_recv
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP,wCC_TMP_recv
double precision, dimension(0:nw+nwauxio)                   :: normconv
logical :: fileopen,first
integer :: itag,Morton_no,ipe,levmin_recv,levmax_recv,igrid_recv,level_recv
integer :: ixrvC^L,ixrvCC^L
integer :: ind_send(2*^ND),ind_recv(2*^ND),siz_ind,igonlevel_recv
integer :: NumGridsOnLevel_mype(1:nlevelshi,0:npe-1)
character(len=80) :: filename
integer ::           filenr
character(len=1024) :: tecplothead

character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead
!-----------------------------------------------------------------------------
if(nw/=count(w_write(1:nw)))then
 if(mype==0) PRINT *,'tecplot_mpi does not use w_write=F'
 call mpistop('w_write, tecplot')
end if

if(nocartesian)then
 if(mype==0) PRINT *,'tecplot_mpi with nocartesian'
endif

Master_cpu_open : if (mype == 0) then
 inquire(qunit,opened=fileopen)
 if (.not.fileopen) then
   ! generate filename
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".plt"
   open(qunit,file=filename,status='unknown')
 end if

 call getheadernames(wnamei,xandwnamei,outfilehead)

 write(tecplothead,'(a)') "VARIABLES = "//TRIM(outfilehead)
 write(qunit,'(a)') tecplothead(1:len_trim(tecplothead))
end if  Master_cpu_open


! determine overall number of grids per level, and the same info per CPU
NumGridsOnLevel(1:nlevelshi)=0
do level=levmin,levmax
   NumGridsOnLevel(level)=0
   do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid = sfc_to_igrid(Morton_no)
      if (node(plevel_,igrid)/=level) cycle
      NumGridsOnLevel(level)=NumGridsOnLevel(level)+1
   end do
   NumGridsOnLevel_mype(level,0:npe-1)=0
   NumGridsOnLevel_mype(level,mype) = NumGridsOnLevel(level)
   call MPI_ALLREDUCE(MPI_IN_PLACE,NumGridsOnLevel_mype(level,0:npe-1),npe,MPI_INTEGER,&
                 MPI_MAX,icomm,ierrmpi)
   call MPI_ALLREDUCE(MPI_IN_PLACE,NumGridsOnLevel(level),1,MPI_INTEGER,MPI_SUM, &
                   icomm,ierrmpi)
end do


!!do level=levmin,levmax
!!  print *,'mype, level en NumGridsOnLevel_mype(level,0:npe-1)=', &
!!     mype,level,NumGridsOnLevel_mype(level,0:npe-1)
!!  print *,'mype, level en NumGridsOnLevel(level)=', &
!!     mype,level,NumGridsOnLevel(level)
!!enddo


nx^D=ixMhi^D-ixMlo^D+1;
nxC^D=nx^D+1;

if(mype==0.and.npe>1) allocate(intstatus(MPI_STATUS_SIZE,1))

{^IFONED
if(convert_type=='teclinempi') then
   nodes=0
   elems=0
   do level=levmin,levmax
      nodes=nodes + NumGridsOnLevel(level)*{nxC^D*}
      elems=elems + NumGridsOnLevel(level)*{nx^D*}
   enddo

   if (mype==0) write(qunit,"(a,i7,a,1pe12.5,a)") &
         'ZONE T="all levels", I=',elems, &
         ', SOLUTIONTIME=',global_time*time_convert_factor,', F=POINT' 

   igonlevel=0
   do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid = sfc_to_igrid(Morton_no)
      call calc_x(igrid,xC,xCC)
      call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,ixC^L,ixCC^L,.true.)
      if (mype==0) then
       {do ix^DB=ixCCmin^DB,ixCCmax^DB\}
            x_TEC(1:ndim)=xCC_TMP(ix^D,1:ndim)*normconv(0)
            w_TEC(1:nw+nwauxio)=wCC_TMP(ix^D,1:nw+nwauxio)*normconv(1:nw+nwauxio)
           write(qunit,fmt="(100(e14.6))") x_TEC, w_TEC
       {end do\}
       else if (mype/=0) then
        itag=Morton_no
        call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
        call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION,0,itag,icomm,ierrmpi)
        call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
        call MPI_SEND(xCC_TMP,1,type_block_xcc_io, 0,itag,icomm,ierrmpi)
       end if
    enddo
    if (mype==0) then
       do ipe=1,npe-1
        do  Morton_no=Morton_start(ipe),Morton_stop(ipe)
          itag=Morton_no
          call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
          call MPI_RECV(normconv,nw+nwauxio+1, MPI_DOUBLE_PRECISION,ipe,&
                        itag,icomm,intstatus(:,1),ierrmpi)
          call MPI_RECV(wCC_TMP_recv,1,type_block_wcc_io, ipe,itag,&
                             icomm,intstatus(:,1),ierrmpi)
          call MPI_RECV(xCC_TMP_recv,1,type_block_xcc_io, ipe,itag,&
                             icomm,intstatus(:,1),ierrmpi)
         {do ix^DB=ixCCmin^DB,ixCCmax^DB\}
             x_TEC(1:ndim)=xCC_TMP_recv(ix^D,1:ndim)*normconv(0)
             w_TEC(1:nw+nwauxio)=wCC_TMP_recv(ix^D,1:nw+nwauxio)*normconv(1:nw+nwauxio)
             write(qunit,fmt="(100(e14.6))") x_TEC, w_TEC
         {end do\}
        end do 
       end do
       close(qunit)
    end if
else
}

if  (mype/=0) then
 itag=1000*Morton_stop(mype)
 call MPI_SEND(levmin,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
 itag=2000*Morton_stop(mype)
 call MPI_SEND(levmax,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
end if

do level=levmin,levmax
   nodesonlevelmype=NumGridsOnLevel_mype(level,mype)*{nxC^D*}
   elemsonlevelmype=NumGridsOnLevel_mype(level,mype)*{nx^D*}
   nodesonlevel=NumGridsOnLevel(level)*{nxC^D*}
   elemsonlevel=NumGridsOnLevel(level)*{nx^D*}
   ! for all tecplot variants coded up here, we let the TECPLOT ZONES coincide
   ! with the AMR grid LEVEL. Other options would be
   !    let each grid define a zone: inefficient for TECPLOT internal workings
   !       hence not implemented
   !    let entire octree define 1 zone: no difference in interpolation 
   !       properties across TECPLOT zones detected as yet, hence not done
   select case(convert_type)
     case('tecplotmpi')
       ! in this option, we store the corner coordinates, as well as the corner
       ! values of all variables (obtained by averaging). This allows POINT packaging, 
       ! and thus we can save full grid info by using one call to calc_grid
       if (mype==0.and.(nodesonlevelmype>0.and.elemsonlevelmype>0))&
        write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,a)") &
             'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
             ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=POINT, ZONETYPE=', &
          {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
      do Morton_no=Morton_start(mype),Morton_stop(mype)
         igrid = sfc_to_igrid(Morton_no)
         if (mype/=0)then
           itag=Morton_no
           call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
           itag=igrid
           call MPI_SEND(node(plevel_,igrid),1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
         end if
         if (node(plevel_,igrid)/=level) cycle
         call calc_x(igrid,xC,xCC)
         call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                        ixC^L,ixCC^L,.true.)
         if (mype/=0) then
            itag=Morton_no
            ind_send=(/ ixC^L /)
            siz_ind=2*^ND
            call MPI_SEND(ind_send,siz_ind, MPI_INTEGER, 0,itag,icomm,ierrmpi)
            call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION, 0,itag,icomm,ierrmpi)

            call MPI_SEND(wC_TMP,1,type_block_wc_io, 0,itag,icomm,ierrmpi)
            call MPI_SEND(xC_TMP,1,type_block_xc_io, 0,itag,icomm,ierrmpi)
         else  
           {do ix^DB=ixCmin^DB,ixCmax^DB\}
              x_TEC(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0)
              w_TEC(1:nw+nwauxio)=wC_TMP(ix^D,1:nw+nwauxio)*normconv(1:nw+nwauxio)
              write(qunit,fmt="(100(e14.6))") x_TEC, w_TEC
           {end do\}
         end if
       enddo

     case('tecplotCCmpi')
       ! in this option, we store the corner coordinates, and the cell center
       ! values of all variables. Due to this mix of corner/cell center, we must 
       ! use BLOCK packaging, and thus we have enormous overhead by using 
       ! calc_grid repeatedly to merely fill values of cell corner coordinates 
       ! and cell center values per dimension, per variable
       if(ndim+nw+nwauxio>99) call mpistop("adjust format specification in writeout")
       if(nw+nwauxio==1)then
         ! to make tecplot happy: avoid [ndim+1-ndim+1] in varlocation varset
         ! and just set [ndim+1]
         if (mype==0.and.(nodesonlevelmype>0.and.elemsonlevelmype>0))&
          write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
            ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
       else
        if(ndim+nw+nwauxio<10) then
         ! difference only in length of integer format specification for ndim+nw+nwauxio
         if (mype==0.and.(nodesonlevelmype>0.and.elemsonlevelmype>0))&
          write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i1,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
            ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,'-',ndim+nw+nwauxio,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
        else
         if (mype==0.and.(nodesonlevelmype>0.and.elemsonlevelmype>0))&
          write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i2,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
            ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,'-',ndim+nw+nwauxio,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
        endif
       endif
       
       do idim=1,ndim
         first=(idim==1)
         do Morton_no=Morton_start(mype),Morton_stop(mype)
          igrid = sfc_to_igrid(Morton_no)
          if (mype/=0)then
           itag=Morton_no*idim
           call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
           itag=igrid*idim
           call MPI_SEND(node(plevel_,igrid),1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
          end if
          if (node(plevel_,igrid)/=level) cycle

          call calc_x(igrid,xC,xCC)
          call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                           ixC^L,ixCC^L,first)
          if (mype/=0)then
            ind_send=(/ ixC^L /)
            siz_ind=2*^ND
            itag=igrid*idim
            call MPI_SEND(ind_send,siz_ind, MPI_INTEGER, 0,itag,icomm,ierrmpi)
            call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION, 0,itag,icomm,ierrmpi)
            call MPI_SEND(xC_TMP,1,type_block_xc_io, 0,itag,icomm,ierrmpi)
          else
            write(qunit,fmt="(100(e14.6))") xC_TMP(ixC^S,idim)*normconv(0)
          end if
         enddo
       enddo
      
       do iw=1,nw+nwauxio
        do Morton_no=Morton_start(mype),Morton_stop(mype)
         igrid = sfc_to_igrid(Morton_no)
         if (mype/=0)then
           itag=Morton_no*(ndim+iw)
           call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
           itag=igrid*(ndim+iw)
           call MPI_SEND(node(plevel_,igrid),1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
         end if
         if (node(plevel_,igrid)/=level) cycle

         call calc_x(igrid,xC,xCC)
         call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                           ixC^L,ixCC^L,.true.)
            
         if (mype/=0)then
            ind_send=(/ ixCC^L /)
            siz_ind=2*^ND
            itag=igrid*(ndim+iw)
            call MPI_SEND(ind_send,siz_ind, MPI_INTEGER, 0,itag,icomm,ierrmpi)
            call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION, 0,itag,icomm,ierrmpi)
            call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
         else
            write(qunit,fmt="(100(e14.6))") wCC_TMP(ixCC^S,iw)*normconv(iw)
         endif
        enddo
       enddo
     case default
       call mpistop('no such tecplot type')
   end select
 

   igonlevel=0
   do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid = sfc_to_igrid(Morton_no)
      if (mype/=0)then
          itag=Morton_no
          call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
          itag=igrid
          call MPI_SEND(node(plevel_,igrid),1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
      end if
      if (node(plevel_,igrid)/=level) cycle

      igonlevel=igonlevel+1
      if (mype/=0)then
          itag=igrid
          call MPI_SEND(igonlevel,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
      end if
      if(mype==0)then
        call save_conntec(qunit,igrid,igonlevel)
      endif
   end do
end do

if (mype==0) then
 if (npe>1) then
  do ipe=1,npe-1
   itag=1000*Morton_stop(ipe)
   call MPI_RECV(levmin_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
   itag=2000*Morton_stop(ipe)
   call MPI_RECV(levmax_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
   do level=levmin_recv,levmax_recv
    nodesonlevelmype=NumGridsOnLevel_mype(level,ipe)*{nxC^D*}
    elemsonlevelmype=NumGridsOnLevel_mype(level,ipe)*{nx^D*}
    nodesonlevel=NumGridsOnLevel(level)*{nxC^D*}
    elemsonlevel=NumGridsOnLevel(level)*{nx^D*}
    select case(convert_type)
     case('tecplotmpi')
        ! in this option, we store the corner coordinates, as well as the corner
        ! values of all variables (obtained by averaging). This allows POINT packaging, 
        ! and thus we can save full grid info by using one call to calc_grid
        if(nodesonlevelmype>0.and.elemsonlevelmype>0) &
        write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,a)") &
             'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
             ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=POINT, ZONETYPE=', &
          {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
        do  Morton_no=Morton_start(ipe),Morton_stop(ipe)
         itag=Morton_no
         call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
         itag=igrid_recv
         call MPI_RECV(level_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
         if (level_recv/=level) cycle

         itag=Morton_no
         siz_ind=2*^ND
         call MPI_RECV(ind_recv,siz_ind, MPI_INTEGER, ipe,itag,&
                        icomm,intstatus(:,1),ierrmpi)
         ixrvCmin^D=ind_recv(^D);ixrvCmax^D=ind_recv(^ND+^D);
         call MPI_RECV(normconv,nw+nwauxio+1, MPI_DOUBLE_PRECISION,ipe,itag&
                  ,icomm,intstatus(:,1),ierrmpi)
     
         call MPI_RECV(wC_TMP_recv,1,type_block_wc_io, ipe,itag,&
                        icomm,intstatus(:,1),ierrmpi)
         call MPI_RECV(xC_TMP_recv,1,type_block_xc_io, ipe,itag,&
                        icomm,intstatus(:,1),ierrmpi)

         {do ix^DB=ixrvCmin^DB,ixrvCmax^DB\}
              x_TEC(1:ndim)=xC_TMP_RECV(ix^D,1:ndim)*normconv(0)
              w_TEC(1:nw+nwauxio)=wC_TMP_RECV(ix^D,1:nw+nwauxio)*normconv(1:nw+nwauxio)
              write(qunit,fmt="(100(e14.6))") x_TEC, w_TEC
         {end do\}
        end do
     case('tecplotCCmpi')
       ! in this option, we store the corner coordinates, and the cell center
       ! values of all variables. Due to this mix of corner/cell center, we must 
       ! use BLOCK packaging, and thus we have enormous overhead by using 
       ! calc_grid repeatedly to merely fill values of cell corner coordinates 
       ! and cell center values per dimension, per variable
       if(ndim+nw+nwauxio>99) call mpistop("adjust format specification in writeout")
       if(nw+nwauxio==1)then
         ! to make tecplot happy: avoid [ndim+1-ndim+1] in varlocation varset
         ! and just set [ndim+1]
         if(nodesonlevelmype>0.and.elemsonlevelmype>0) &
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
            ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
       else
        if(ndim+nw+nwauxio<10) then
         ! difference only in length of integer format specification for ndim+nw+nwauxio
         if(nodesonlevelmype>0.and.elemsonlevelmype>0) &
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i1,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
            ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,'-',ndim+nw+nwauxio,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
        else
         if(nodesonlevelmype>0.and.elemsonlevelmype>0) &
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i2,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
            ', SOLUTIONTIME=',global_time*time_convert_factor,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,'-',ndim+nw+nwauxio,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
        endif
       endif

       do idim=1,ndim
         do  Morton_no=Morton_start(ipe),Morton_stop(ipe)
           itag=Morton_no*idim
           call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
           itag=igrid_recv*idim
           call MPI_RECV(level_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
           if (level_recv/=level) cycle
           
           siz_ind=2*^ND
           itag=igrid_recv*idim
           call MPI_RECV(ind_recv,siz_ind, MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
           ixrvCmin^D=ind_recv(^D);ixrvCmax^D=ind_recv(^ND+^D);     
           call MPI_RECV(normconv,nw+nwauxio+1, MPI_DOUBLE_PRECISION,ipe,itag&
                  ,icomm,intstatus(:,1),ierrmpi)
           call MPI_RECV(xC_TMP_recv,1,type_block_xc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
           write(qunit,fmt="(100(e14.6))") xC_TMP_recv(ixrvC^S,idim)*normconv(0)
         end do
       end do
    
       do iw=1,nw+nwauxio
        do Morton_no=Morton_start(ipe),Morton_stop(ipe)
           itag=Morton_no*(ndim+iw)
           call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
           itag=igrid_recv*(ndim+iw)
           call MPI_RECV(level_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
           if (level_recv/=level) cycle

           siz_ind=2*^ND
           itag=igrid_recv*(ndim+iw)
           call MPI_RECV(ind_recv,siz_ind, MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
           ixrvCCmin^D=ind_recv(^D);ixrvCCmax^D=ind_recv(^ND+^D);
           call MPI_RECV(normconv,nw+nwauxio+1, MPI_DOUBLE_PRECISION,ipe,itag&
                  ,icomm,intstatus(:,1),ierrmpi)
           call MPI_RECV(wCC_TMP_recv,1,type_block_wcc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
           write(qunit,fmt="(100(e14.6))") wCC_TMP_recv(ixrvCC^S,iw)*normconv(iw)
        enddo
       enddo
     case default
       call mpistop('no such tecplot type')
    end select

    do Morton_no=Morton_start(ipe),Morton_stop(ipe)
      itag=Morton_no
      call MPI_RECV(igrid_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
      itag=igrid_recv
      call MPI_RECV(level_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
      if (level_recv/=level) cycle

      itag=igrid_recv
      call MPI_RECV(igonlevel_recv,1,MPI_INTEGER, ipe,itag,icomm,intstatus(:,1),ierrmpi)
      call save_conntec(qunit,igrid_recv,igonlevel_recv)
    end do ! morton loop
   end do ! level loop
  end do ! ipe loop
 end if ! npe>1 if
end if ! mype=0 if
{^IFONED endif}

if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
  if(mype==0)deallocate(intstatus)
endif

end subroutine tecplot_mpi
!=============================================================================
subroutine punstructuredvtkB_mpi(qunit)

! Write one pvtu and vtu files for each processor
! Otherwise like unstructuredvtk_mpi
! output for vtu format to paraview, binary version output
! uses calc_grid to compute nwauxio variables
! allows renormalizing using convert factors

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters
use mod_calculate_xw

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP

integer :: igrid,iigrid,level,igonlevel,icel,ixC^L,ixCC^L,Morton_no
integer ::               NumGridsOnLevel(1:nlevelshi)
integer :: nx^D,nxC^D,nodesonlevel,elemsonlevel,nc,np,VTK_type,ix^D
double precision :: normconv(0:nw+nwauxio)
character(len=80) :: pfilename
character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

integer*8 :: offset
integer::  recsep,k,iw,filenr
integer::  length,lengthcc,offset_points,offset_cells, &
           length_coords,length_conn,length_offsets
character::  buf
character(len=6)::  bufform

logical ::   fileopen
!-----------------------------------------------------------------------------

! Write pvtu-file:
if (mype==0) then
   call write_pvtu(qunit)
endif
! Now write the Source files:
inquire(qunit,opened=fileopen)
if(.not.fileopen)then
   ! generate filename 
   filenr=snapshotnext-1
   if (autoconvert) filenr=snapshotnext
   ! Open the file for the header part
   write(pfilename,'(a,i4.4,a,i4.4,a)') TRIM(base_filename),filenr,"p",mype,".vtu"
   open(qunit,file=pfilename,status='unknown',form='formatted')
endif
! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
write(qunit,'(a)')'  <UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(qunit,*) real(global_time*time_convert_factor)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'


offset=0
recsep=4

call getheadernames(wnamei,xandwnamei,outfilehead)

! number of cells, number of corner points, per grid.
nx^D=ixMhi^D-ixMlo^D+1;
nxC^D=nx^D+1;
nc={nx^D*}
np={nxC^D*}

length=np*size_real
lengthcc=nc*size_real

length_coords=3*length
length_conn=2**^ND*size_int*nc
length_offsets=nc*size_int

! Note: using the w_write, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
do level=levmin,levmax
 if (writelevel(level)) then
   do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    if (node(plevel_,igrid)/=level) cycle
    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    if (({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
          *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
         <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})) then
      select case(convert_type)
       case('pvtuBmpi')
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
            offset=offset+length+size_int
         enddo
         do iw=nw+1,nw+nwauxio

            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+length+size_int
         enddo
         write(qunit,'(a)')'</PointData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
     '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_int
         write(qunit,'(a)')'</Points>'
       case('pvtuBCCmpi')
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
            offset=offset+lengthcc+size_int
         enddo
         do iw=nw+1,nw+nwauxio

            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+lengthcc+size_int
         enddo
         write(qunit,'(a)')'</CellData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
     '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_int
         write(qunit,'(a)')'</Points>'
      end select

   
      write(qunit,'(a)')'<Cells>'

      ! connectivity part
      write(qunit,'(a,i16,a)')&
        '<DataArray type="Int32" Name="connectivity" format="appended" offset="',offset,'"/>'
      offset=offset+length_conn+size_int    

      ! offsets data array
      write(qunit,'(a,i16,a)') &
        '<DataArray type="Int32" Name="offsets" format="appended" offset="',offset,'"/>'
      offset=offset+length_offsets+size_int    

      ! VTK cell type data array
      write(qunit,'(a,i16,a)') &
        '<DataArray type="Int32" Name="types" format="appended" offset="',offset,'"/>' 
      offset=offset+size_int+nc*size_int

      write(qunit,'(a)')'</Cells>'

      write(qunit,'(a)')'</Piece>'
    endif
   enddo
 endif
enddo

write(qunit,'(a)')'</UnstructuredGrid>'
write(qunit,'(a)')'<AppendedData encoding="raw">'

close(qunit)
! next to make gfortran compiler happy, as it does not know
! form='binary' and produces error on compilation
!bufform='binary'
!open(qunit,file=pfilename,form=bufform,position='append')
!This should in principle do also for gfortran (tested with gfortran 4.6.0 and Intel 11.1):
open(qunit,file=pfilename,access='stream',form='unformatted',position='append')
buf='_'
write(qunit) TRIM(buf)

do level=levmin,levmax
 if (writelevel(level)) then
   do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    if (node(plevel_,igrid)/=level) cycle
      ! only output a grid when fully within clipped region selected
      ! by writespshift array
      if (({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
          *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
         <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})) then
        call calc_x(igrid,xC,xCC)
        call calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                       ixC^L,ixCC^L,.true.)
        do iw=1,nw
          if(.not.w_write(iw))cycle
          select case(convert_type)
            case('pvtuBmpi')
              write(qunit) length
              write(qunit) {(|}real(wC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCmin^D,ixCmax^D)}
            case('pvtuBCCmpi')
              write(qunit) lengthcc
              write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
          end select 
        enddo
        do iw=nw+1,nw+nwauxio
          select case(convert_type)
            case('pvtuBmpi')
              write(qunit) length
              write(qunit) {(|}real(wC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCmin^D,ixCmax^D)}
            case('pvtuBCCmpi')
              write(qunit) lengthcc
              write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
          end select 
        enddo

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
    endif
  end do
 endif
end do

close(qunit)
open(qunit,file=pfilename,status='unknown',form='formatted',position='append')

write(qunit,'(a)')'</AppendedData>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)

end subroutine punstructuredvtkB_mpi
{^IFTWOD
! subroutines to convert 2.5D data to 3D data
!====================================================================================
subroutine unstructuredvtkB23(qunit)

! output for vtu format to paraview, binary version output
! not parallel, uses calc_grid to compute nwauxio variables

use mod_global_parameters
use mod_physics
use mod_calculate_xw

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
character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:3+nw+nwauxio)
character(len=1024) :: outfilehead
integer*8 :: offset
integer :: size_length,recsep,k,iw
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
size_length=4

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
  ! generate filename 
  write(filename,'(a,a,i4.4,a)') TRIM(base_filename),"3D",snapshotini,".vtu"
  ! Open the file for the header part
  open(qunit,file=filename,status='replace')
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

length=np*size_real
lengthcc=nc*size_real

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
      n3grid=nint(zlength/d3grid)
      do i3grid=1,n3grid !subcycles
      select case(convert_type)
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
            offset=offset+length+size_int
         enddo
         if(nwauxio>0)then
          do iw=nw+1,nw+nwauxio
             write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
                TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+length+size_int
          enddo
         endif
         write(qunit,'(a)')'</PointData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
            offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_int
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
            offset=offset+lengthcc+size_int
         enddo
         if(nwauxio>0)then
          do iw=nw+1,nw+nwauxio
             write(qunit,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
                TRIM(wnamei(iw)), '" format="appended" offset="',offset,'">'
             write(qunit,'(a)')'</DataArray>'
             offset=offset+lengthcc+size_int
          enddo
         endif
         write(qunit,'(a)')'</CellData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
'<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',&
            offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_int
         write(qunit,'(a)')'</Points>'
      end select

      write(qunit,'(a)')'<Cells>'

      ! connectivity part
      write(qunit,'(a,i16,a)')&
'<DataArray type="Int32" Name="connectivity" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_conn+size_int

      ! offsets data array
      write(qunit,'(a,i16,a)') &
         '<DataArray type="Int32" Name="offsets" format="appended" offset="',&
         offset,'"/>'
      offset=offset+length_offsets+size_int

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
       n3grid=nint(zlength/d3grid)
       ! In case primitives to be saved: use primitive subroutine
       !  extra layer around mesh only needed when storing corner values and averaging
       if(saveprim) then
         call phys_to_primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
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
          select case(convert_type)
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
           select case(convert_type)
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
         end do
        end if
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
        end do
       end do !subcycles
      end if
   end do
 end if
end do

close(qunit)
open(qunit,file=filename,status='unknown',form='formatted',position='append')

write(qunit,'(a)')'</AppendedData>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)
end subroutine unstructuredvtkB23
!====================================================================================
subroutine unstructuredvtkBsym23(qunit)

! output for vtu format to paraview, binary version output
! not parallel, uses calc_grid to compute nwauxio variables
! use this subroutine  when the physical domain is symmetric/asymmetric about (0,y,z) 
! plane, xprobmin1=0 and the computational domain is a half of the physical domain

use mod_global_parameters
use mod_physics
use mod_calculate_xw

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
character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:3+nw+nwauxio)
character(len=1024) :: outfilehead
integer*8 :: offset
integer :: size_length,recsep,k,iw
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
size_length=4

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

length=np*size_real
lengthcc=nc*size_real

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
     n3grid=nint(zlength/d3grid)
     do i3grid=1,n3grid !subcycles
      !! original domain ----------------------------------start 
      select case(convert_type)
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
      select case(convert_type)
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
       n3grid=nint(zlength/d3grid)
       ! In case primitives to be saved: use primitive subroutine
       !  extra layer around mesh only needed when storing corner values and averaging
       if(saveprim) then
         call phys_to_primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
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
          select case(convert_type)
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
           select case(convert_type)
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
        end do
        !! original domain ----------------------------------end 
        !! symetric/asymetric mirror domain -----------------start
        do iw=1,nw
          if(.not.w_write(iw))cycle
          if(iw==2 .or. iw==4 .or. iw==7) then
           wC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)=&
           -wC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
           wCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,iw)=&
           -wCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,iw)
          end if
          select case(convert_type)
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
           select case(convert_type)
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
         end do
        end if
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
        end do
        !! symetric/asymetric mirror domain -----------------end
       end do !subcycles
      end if
   end do
 end if
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
logical, save :: subfirst=.true.
!-----------------------------------------------------------------------------
! following only for allowing compiler to go through with debug on

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
  ! default (no) normalization for auxiliary variables
  normconv(nw+1:nw+nwauxio)=one
  ! maybe need for restriction to ixG^LL^LSUB1 
  call specialvar_output23(ixGlo1,ixGlo2,ixGlo1,ixGhi1,ixGhi2,ixGhi1,ixGlo1&
     +1,ixGlo2+1,ixGlo1+1,ixGhi1-1,ixGhi2-1,ixGhi1-1,w,xCC,normconv)
endif

! compute the cell-center values for w first
!===========================================
! cell center values obtained from mere copy, while B0+B1 split handled here

wCC(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,:)=w(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,:)
if(B0field) then
  do ix3=ixCCmin3,ixCCmax3
    do ix2=ixCCmin2,ixCCmax2
      do ix1=ixCCmin1,ixCCmax1
        wCC(ix1,ix2,ix3,iw_mag(:))=wCC(ix1,ix2,ix3,iw_mag(:))+ps(igrid)%B0(ix1,ix2,&
            :,0)
      end do
    end do
  end do
end if

if(.not.saveprim .and. B0field .and. iw_e>0) then
   do ix3=ixCCmin3,ixCCmax3
   do ix2=ixCCmin2,ixCCmax2
   do ix1=ixCCmin1,ixCCmax1
       wCC(ix1,ix2,ix3,iw_e)=w(ix1,ix2,ix3,iw_e) +half*sum(ps(igrid)%B0(ix1,&
          ix2,:,0)**2 ) + sum(w(ix1,ix2,ix3,&
          iw_mag(:))*ps(igrid)%B0(ix1,ix2,:,0))
   end do
   end do
   end do
end if
! compute the corner values for w now by averaging
!=================================================
if(slab_uniform)then
   ! for slab symmetry: no geometrical info required
   do iw=1,nw+nwauxio
      if (B0field.and.iw>iw_mag(1)-1.and.iw<=iw_mag(ndir)) then
         idir=iw-iw_mag(1)+1
         do ix3=ixCmin3,ixCmax3
         do ix2=ixCmin2,ixCmax2
         do ix1=ixCmin1,ixCmax1
           wC(ix1,ix2,ix3,iw)=sum(w(ix1:ix1+1,ix2:ix2+1,ix3,iw) &
              +ps(igrid)%B0(ix1:ix1+1,ix2:ix2+1&
              ,idir,0))/dble(2**3)+&
              sum(w(ix1:ix1+1,ix2:ix2+1,ix3+1,iw) &
              +ps(igrid)%B0(ix1:ix1+1,ix2:ix2+1&
              ,idir,0))/dble(2**3)
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
   if(.not.saveprim .and. B0field .and. iw_e>0) then
      do ix3=ixCmin3,ixCmax3
      do ix2=ixCmin2,ixCmax2
      do ix1=ixCmin1,ixCmax1
         wC(ix1,ix2,ix3,iw_e)=sum( w(ix1:ix1+1,ix2:ix2+1,ix3,iw_e) &
            +half*sum(ps(igrid)%B0(ix1:ix1+1,ix2:ix2+1&
            ,:,0)**2,dim=ndim+1) + sum( w(ix1:ix1+1,ix2:ix2+1,ix3&
            ,iw_mag(:))*ps(igrid)%B0(ix1:ix1+1,ix2:ix2+1&
            ,:,0),dim=ndim+1) ) /dble(2**3)+&
            sum( w(ix1:ix1+1,ix2:ix2+1,ix3+1,iw_e) &
            +half*sum(ps(igrid)%B0(ix1:ix1+1,ix2:ix2+1&
            ,:,0)**2,dim=ndim+1) + sum( w(ix1:ix1+1,ix2:ix2+1,ix3&
            +1,iw_mag(:))*ps(igrid)%B0(ix1:ix1+1,ix2:ix2+1&
            ,:,0),dim=ndim+1) ) /dble(2**3)
      end do
      end do
      end do
   end if
end if
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
!if(saveprim)then
!  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1)=w(ixOmin1:ixOmax1,&
!   ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_e)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
!   ixOmin3:ixOmax3,iw_rho)
!endif
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
end subroutine specialvar_output23 
\}
!=============================================================================
