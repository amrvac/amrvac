!=============================================================================
subroutine generate_plotfile

{#IFDEF PARTICLES
use mod_timing, only: tpartc, tpartc0
}
include 'amrvacdef.f'
!-----------------------------------------------------------------------------

if(mype==0.and.level_io>0)write(unitterm,*)'reset tree to fixed level=',level_io
if(level_io>0 .or. level_io_min.ne.1 .or. level_io_max.ne.nlevelshi) then 
   call resettree_convert
end if

call getbc(t,ixG^LL,pw,pwCoarse,pgeo,pgeoCoarse,.false.,0,nwflux+nwaux)

!!!call Global_useroutput !compute at user level any global variable over all grids

select case(convert_type)
  case('idl','idlCC')
   call valout_idl(unitconvert)
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
  case('pvtumpi','pvtuCCmpi')
   call punstructuredvtk_mpi(unitconvert)
  case('pvtuBmpi','pvtuBCCmpi')
   call punstructuredvtkB_mpi(unitconvert)
  case('vtimpi','vtiCCmpi')
   call ImageDataVtk_mpi(unitconvert)
  case('dx')
   call valout_dx(unitconvert)
  case('onegrid','onegridmpi')
   call onegrid(unitconvert)
  case('oneblock','oneblockB')
   call oneblock(unitconvert)
{#IFDEF PARTICLES
  case('particles', 'particlesmpi')
     tpartc0 = MPI_WTIME()
     call handle_particles()
     tpartc = tpartc + (MPI_WTIME() - tpartc0)
}
{#IFDEF UCONVERT
  case('user','usermpi')
   call userspecialconvert(unitconvert)
}
  case default
   call mpistop("Error in generate_plotfile: Unknown convert_type")
end select

end subroutine generate_plotfile
!=============================================================================
subroutine getheadernames(wnamei,xandwnamei,outfilehead)

! this collects all variables names in the wnamei character array, getting the info from
! the primnames/wnames strings (depending on saveprim). It combines this info with names
! for the dimensional directions in the xandwnamei array. In the outfilehead, it collects
! the dimensional names, and only those names from the nw variables for output (through writew)
! together with the added names for nwauxio variables

include 'amrvacdef.f'

character(len=10)   :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

integer::  space_position,iw
character(len=10)::  wname
character(len=len(primnames))::  scanstring

logical, save:: first=.true.
!-----------------------------------------------------------------------------

! in case additional variables are computed and stored for output, adjust 
! the wnames and primnames string
if(nwauxio>0 .and. first) call specialvarnames_output

! --- part to provide variable names from primnames/varnames strings
if(saveprim) then
   scanstring=TRIM(primnames)
else
   scanstring=TRIM(wnames)
endif

space_position=index(scanstring,' ')
do iw=1,nw+nwauxio
   do while (space_position==1)
     scanstring=scanstring(2:)
     space_position=index(scanstring,' ')
   enddo
   wname=scanstring(:space_position-1)
   scanstring=scanstring(space_position+1:)
   space_position=index(scanstring,' ')

   ! fill all names, even those that we will not write away from the first nw
   wnamei(iw)=TRIM(wname)
enddo
! --- end of part to provide variable names 

select case (typeaxial)
   case( "spherical" )
      xandwnamei(1)="r";{^NOONED xandwnamei(2)="Theta"};{^IFTHREED xandwnamei(3)="Phi"}
   case( "cylindrical" )
      xandwnamei(1)="R";
      {^NOONED
      if( ^PHI == 2 )then
         xandwnamei(2)="Phi"
      else
         xandwnamei(2)="Z"
      endif}
      {^IFTHREED
      if( ^PHI == 2 )then
         xandwnamei(3)="Z"
      else
         xandwnamei(3)="Phi"
      endif}
   case default
      xandwnamei(1)="X";{^NOONED xandwnamei(2)="Y"};{^IFTHREED xandwnamei(3)="Z"}
end select

xandwnamei(ndim+1:ndim+nw+nwauxio)=wnamei(1:nw+nwauxio)

! in outfilehead, collect the dimensional names, and all output variable names
! first all dimensions
write(outfilehead,'(a)') TRIM(xandwnamei(1))
{^NOONED
do iw=2,ndim
   wname=xandwnamei(iw)
write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(wname)
enddo
}
! then all nw variables, with writew control for inclusion
do iw=ndim+1,ndim+nw
   wname=xandwnamei(iw)
   if(writew(iw-ndim)) then
write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(wname)
   endif
enddo
! then all nwauxio variables
if(nwauxio>0) then
  do iw=ndim+nw+1,ndim+nw+nwauxio
     wname=xandwnamei(iw)
write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(wname)
  enddo
endif

if(first.and.mype==0)then
  print*,'-----------------------------------------------------------------------------'
  write(unitterm,*)'Saving visual data. Coordinate directions and variable names are:'
  do iw=1,ndim
    print *,iw,xandwnamei(iw)
  enddo
  do iw=ndim+1,ndim+nw+nwauxio
    print *,iw,wnamei(iw-ndim),xandwnamei(iw)
  enddo
  write(unitterm,*)'time =', t
  print*,'-----------------------------------------------------------------------------'
  first=.false.
endif

end subroutine getheadernames
!=============================================================================
subroutine oneblock(qunit)

! this is for turning an AMR run into a single block
! the data will be all on selected level level_io

! this version should work for any dimension
! only writes writew selected 1:nw variables, also nwauxio
! may use saveprim to switch to primitives
! this version can not work on multiple CPUs
! does not renormalize variables

! header info differs from onegrid below 

! ASCII or binary output

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid, igrid_to_node
include 'amrvacdef.f'
integer, intent(in) :: qunit

integer             :: Morton_no,igrid,ix^D,ig^D,level
integer, pointer    :: ig_to_igrid(:^D&,:)
logical             :: fileopen
character(len=80)   :: filename
integer             :: filenr

character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
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

if(npe>1)then
 if(mype==0) PRINT *,'ONEBLOCK as yet to be parallelized'
 call mpistop('npe>1, oneblock')
end if

! only variables selected by writew will be written out
normconv(0:nw+nwauxio)=one
normconv(0:nw)=normvar(0:nw)
writenw=count(writew(1:nw))+nwauxio
iiw=0
do iw =1,nw
 if (.not.writew(iw))cycle
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
  call primitive(ixG^LL,ixG^LL^LSUB1,pw(igrid)%w,px(igrid)%x)
 end do
else
 if (nwaux>0) then
  do iigrid=1,igridstail; igrid=igrids(iigrid)
   call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixG^LL^LSUB1,"oneblock")
  end do
 end if
end if


Master_cpu_open : if (mype == 0) then
 inquire(qunit,opened=fileopen)
 if (.not.fileopen) then
   ! generate filename
    filenr=snapshotini
    if (.not.convert) filenr=snapshot-1
   write(filename,'(a,i4.4,a)') TRIM(filenameout),filenr,".blk"
   select case(convert_type)
    case("oneblock")
     open(qunit,file=filename,status='unknown')
     write(qunit,*) TRIM(outfilehead)
     write(qunit,*)( {^D&(ixMhi^D-ixMlo^D+1)*})*(Morton_stop(npe-1)-Morton_start(0)+1),&
                    {ng^D(level_io)*(ixMhi^D-ixMlo^D+1)|,}
     write(qunit,*)t*normt
    case("oneblockB")
     open(qunit,file=filename,form='unformatted',status='unknown')
     write(qunit) outfilehead
     write(qunit)  ( {^D&(ixMhi^D-ixMlo^D+1)*})*(Morton_stop(npe-1)-Morton_start(0)+1),&
                    {ng^D(level_io)*(ixMhi^D-ixMlo^D+1)|,}
     write(qunit)t*normt
   end select
 end if
end if Master_cpu_open

{^IFTHREED 
do ig3=1,ng3(level_io)}
   {^NOONED
   do ig2=1,ng2(level_io)}
       do ig1=1,ng1(level_io)
         igrid=ig_to_igrid(ig^D,mype)
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
         allocate(pwio(igrid)%w(ixG^T,1:nw+nwauxio))
         pwio(igrid)%w(ixG^T,1:nw)=pw(igrid)%w(ixG^T,1:nw)
         if(nwauxio>=1)then
            call specialvar_output(ixG^LL,ixM^LL^LADD1,pwio(igrid)%w,px(igrid)%x,normconv)
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

       do ig1=1,ng1(level_io)
         do ix1=ixMlo1,ixMhi1
           igrid=ig_to_igrid(ig^D,mype)
           if (B0field) then
             myB0_cell => pB0_cell(igrid)
           end if
           Master_write : if(mype==0) then
           {^IFMHD
             if(B0field) then
               ^C&pwio(igrid)%w(ix^D,b^C_)=pwio(igrid)%w(ix^D,b^C_)+myB0_cell%w(ix^D,^C);\
             end if
           }
             select case(convert_type)
               case("oneblock")
                 write(qunit,fmt="(100(e14.6))") &
                  px(igrid)%x(ix^D,1:ndim)*normconv(0),&
                  (pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw)),iw=1,writenw)
               case("oneblockB")
                 write(qunit) real(px(igrid)%x(ix^D,1:ndim)*normconv(0)),&
                  (real(pwio(igrid)%w(ix^D,iwrite(iw))*normconv(iwrite(iw))),iw=1,writenw)
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
  call conserve(ixG^LL,ixG^LL^LSUB1,pw(igrid)%w,px(igrid)%x,patchw)
 end do
endif

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
include 'amrvacdef.f'
integer, intent(in) :: qunit

integer             :: Morton_no,igrid,ix^D,iw
logical             :: fileopen
character(len=80)   :: filename
integer             :: filenr

character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
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
    if (.not.convert) filenr=snapshot-1
   write(filename,'(a,i4.4,a)') TRIM(filenameout),filenr,".blk"
   open(qunit,file=filename,status='unknown')
 end if
 write(qunit,"(a)")outfilehead
 write(qunit,"(i7)") ( {^D&(ixMhi^D-ixMlo^D+1)*} )*(Morton_stop(npe-1)-Morton_start(0)+1)
end if Master_cpu_open

do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  if(saveprim) call primitive(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x)
  if (mype/=0)then
      itag=Morton_no
      call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
      call MPI_SEND(px(igrid)%x,1,type_block_xcc_io, 0,itag,icomm,ierrmpi)
      itag=igrid
      call MPI_SEND(pw(igrid)%w,1,type_block_io, 0,itag,icomm,ierrmpi)
  else
   {do ix^DB=ixMlo^DB,ixMhi^DB\}
      do iw=1,nw
        if( dabs(pw(igrid)%w(ix^D,iw)) < 1.0d-32 ) pw(igrid)%w(ix^D,iw) = zero
      enddo
       write(qunit,fmt="(100(e14.6))") px(igrid)%x(ix^D,1:ndim)&
                                     ,pw(igrid)%w(ix^D,1:nw)
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
              if( dabs(pw(igrid)%w(ix^D,iw)) < smalldouble ) pw(igrid)%w(ix^D,iw) = zero
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
subroutine valout_idl(qunit)

! output for idl macros from (amr)vac
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using normt and normvar-array

! binary output format

use mod_forest, only: nleafs
include 'amrvacdef.f'

integer, intent(in) :: qunit

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP

character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical :: fileopen
integer :: iigrid,igrid,nx^D,nxC^D,ixC^L,ixCC^L,iw
character(len=80) :: filename
integer :: filenr

!!! length mismatch possible
character(len=80) :: tmpnames

double precision:: rnode_IDL(rnodehi), normconv(0:nw+nwauxio)
!-----------------------------------------------------------------------------

if(npe>1)then
 if(mype==0) PRINT *,'valoutidl as yet to be parallelized'
 call mpistop('npe>1, valoutidl')
end if

if(nw/=count(writew(1:nw)))then
 if(mype==0) PRINT *,'valoutidl does not use writew=F'
 call mpistop('writew, valoutidl')
end if

inquire(qunit,opened=fileopen)
if (.not.fileopen) then
   ! generate filename
    filenr=snapshotini
    if (.not.convert) filenr=snapshot-1
   write(filename,'(a,i4.4,a)') TRIM(filenameout),filenr,".out"
   open(qunit,file=filename,status='unknown',form='unformatted')
end if

write(qunit)fileheadout
write(qunit)it,t*normt,ndim,neqpar+nspecialpar,nw+nwauxio

nx^D=ixMhi^D-ixMlo^D+1;
select case(convert_type)
  case('idl')
    ! store cell corner quantities, involves averaging to corners
    nxC^D=nx^D+1;
  case('idlCC')
    ! store cell center quantities, do not average to corners
    nxC^D=nx^D;
  case default
   call mpistop("Error in valout_idl: Unknown convert_type")
end select

call getheadernames(wnamei,xandwnamei,outfilehead)

! for idl output: add the eqparnames, note the length mismatch!!!
tmpnames=TRIM(outfilehead)//' '//TRIM(eqparname)//' '//TRIM(specialparname)

! use -nleafs to indicate amr grid
if (nleafs==1 .and. mxnest==1) then
  write(qunit) nxC^D
  write(qunit)eqpar
  write(qunit)tmpnames
else
  write(qunit) ^D&-nleafs
  write(qunit) eqpar
  write(qunit)tmpnames
  ! write out individual grid sizes, grid level, and corners
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     write(qunit) nxC^D
     write(qunit) node(plevel_,igrid)
     select case (typeaxial)
      case ("slab","slabtest")
      {rnode_IDL(rpxmin1_:rpxmin^ND_)=rnode(rpxmin1_:rpxmin^ND_,igrid)};
      {rnode_IDL(rpxmax1_:rpxmax^ND_)=rnode(rpxmax1_:rpxmax^ND_,igrid)}; 
      case ("cylindrical")
      {^IFONED
       rnode_IDL(rpxmin1_)=rnode(rpxmin1_,igrid)
       rnode_IDL(rpxmax1_)=rnode(rpxmax1_,igrid)}
      {^IFTWOD
       select case (phi_)
          case (2) 
        rnode_IDL(rpxmin1_)=rnode(rpxmin1_,igrid)*cos(rnode(rpxmin2_,igrid))
        rnode_IDL(rpxmax1_)=rnode(rpxmax1_,igrid)*cos(rnode(rpxmax2_,igrid))
        rnode_IDL(rpxmin2_)=rnode(rpxmin1_,igrid)*sin(rnode(rpxmin2_,igrid))
        rnode_IDL(rpxmax2_)=rnode(rpxmax1_,igrid)*sin(rnode(rpxmax2_,igrid))
        ! next avoids that IDL clips off part of the domain at circle edges
        if(sin(rnode(rpxmin2_,igrid))*sin(rnode(rpxmax2_,igrid))< zero)then
         if(cos(rnode(rpxmin2_,igrid))>zero)then
           rnode_IDL(rpxmin1_)=rnode(rpxmin1_,igrid)
           rnode_IDL(rpxmax1_)=rnode(rpxmax1_,igrid)
         else
           rnode_IDL(rpxmin1_)=-rnode(rpxmax1_,igrid)
           rnode_IDL(rpxmax1_)=rnode(rpxmin1_,igrid)
         endif
        endif
        if(cos(rnode(rpxmin2_,igrid))*cos(rnode(rpxmax2_,igrid))< zero)then
         if(sin(rnode(rpxmin2_,igrid))>zero)then
           rnode_IDL(rpxmin2_)=rnode(rpxmin1_,igrid)
           rnode_IDL(rpxmax2_)=rnode(rpxmax1_,igrid)
         else
           rnode_IDL(rpxmin2_)=-rnode(rpxmax1_,igrid)
           rnode_IDL(rpxmax2_)=rnode(rpxmin1_,igrid)
         endif
        endif
          case default
      {rnode_IDL(rpxmin1_:rpxmin^ND_)=rnode(rpxmin1_:rpxmin^ND_,igrid)};
      {rnode_IDL(rpxmax1_:rpxmax^ND_)=rnode(rpxmax1_:rpxmax^ND_,igrid)}; 
       end select}
      {^IFTHREED
        rnode_IDL(rpxmin1_)=rnode(rpxmin1_,igrid)*cos(rnode(rpxmin2_,igrid))
        rnode_IDL(rpxmax1_)=rnode(rpxmax1_,igrid)*cos(rnode(rpxmax2_,igrid))
        rnode_IDL(rpxmin2_)=rnode(rpxmin1_,igrid)*sin(rnode(rpxmin2_,igrid))
        rnode_IDL(rpxmax2_)=rnode(rpxmax1_,igrid)*sin(rnode(rpxmax2_,igrid))
        rnode_IDL(rpxmin3_)=rnode(rpxmin3_,igrid)
        rnode_IDL(rpxmax3_)=rnode(rpxmax3_,igrid)
        if(sin(rnode(rpxmin2_,igrid))*sin(rnode(rpxmax2_,igrid))< zero)then
         if(cos(rnode(rpxmin2_,igrid))>zero)then
           rnode_IDL(rpxmin1_)=rnode(rpxmin1_,igrid)
           rnode_IDL(rpxmax1_)=rnode(rpxmax1_,igrid)
         else
           rnode_IDL(rpxmin1_)=-rnode(rpxmax1_,igrid)
           rnode_IDL(rpxmax1_)=rnode(rpxmin1_,igrid)
         endif
        endif
        if(cos(rnode(rpxmin2_,igrid))*cos(rnode(rpxmax2_,igrid))< zero)then
         if(sin(rnode(rpxmin2_,igrid))>zero)then
           rnode_IDL(rpxmin2_)=rnode(rpxmin1_,igrid)
           rnode_IDL(rpxmax2_)=rnode(rpxmax1_,igrid)
         else
           rnode_IDL(rpxmin2_)=-rnode(rpxmax1_,igrid)
           rnode_IDL(rpxmax2_)=rnode(rpxmin1_,igrid)
         endif
        endif}
      case ("spherical")
       rnode_IDL(rpxmin1_)=rnode(rpxmin1_,igrid) &
         {^NOONED*sin(rnode(rpxmin2_,igrid))}{^IFTHREED*cos(rnode(rpxmin3_,igrid))}
       rnode_IDL(rpxmax1_)=rnode(rpxmax1_,igrid) &
         {^NOONED*sin(rnode(rpxmax2_,igrid))}{^IFTHREED*cos(rnode(rpxmax3_,igrid))}
       {^IFTWOD
       rnode_IDL(rpxmin2_)=rnode(rpxmin1_,igrid)*cos(rnode(rpxmin2_,igrid))
       rnode_IDL(rpxmax2_)=rnode(rpxmax1_,igrid)*cos(rnode(rpxmax2_,igrid))}
       {^IFTHREED
       rnode_IDL(rpxmin2_)=rnode(rpxmin1_,igrid) &
                *sin(rnode(rpxmin2_,igrid))*sin(rnode(rpxmin3_,igrid))
       rnode_IDL(rpxmax2_)=rnode(rpxmax1_,igrid) &
                *sin(rnode(rpxmax2_,igrid))*sin(rnode(rpxmax3_,igrid))
       rnode_IDL(rpxmin3_)=rnode(rpxmin1_,igrid)*cos(rnode(rpxmin2_,igrid))
       rnode_IDL(rpxmax3_)=rnode(rpxmax1_,igrid)*cos(rnode(rpxmax2_,igrid))}
     end select 
     normconv(0)=normvar(0)
     write(qunit) {rnode_IDL(rpxmin^D_)*normconv(0)},{rnode_IDL(rpxmax^D_)*normconv(0)}
  end do
end if

! write out variable values
do iigrid=1,igridstail; igrid=igrids(iigrid);
  call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,ixC^L,ixCC^L,.true.)
  select case(convert_type)
  case('idl')
    ! write out corner coordinates and (averaged) corner values
    write(qunit) xC_TMP(ixC^S,1:ndim)*normconv(0)
    do iw=1,nw+nwauxio
      write(qunit) wC_TMP(ixC^S,iw)*normconv(iw)
    end do
  case('idlCC')
    ! write out cell center coordinates and cell center values
    write(qunit) xCC_TMP(ixCC^S,1:ndim)*normconv(0)
    do iw=1,nw+nwauxio
      write(qunit) wCC_TMP(ixCC^S,iw)*normconv(iw)
    end do
  end select 
end do

end subroutine valout_idl
!=============================================================================
subroutine tecplot(qunit)

! output for tecplot (ASCII format)
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using normt and normvar-array

include 'amrvacdef.f'

integer, intent(in) :: qunit

integer::               igrid,iigrid,level,igonlevel,iw,idim,ix^D
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx^D,nxC^D,nodesonlevel,elemsonlevel,ixC^L,ixCC^L

integer ::              nodes, elems


double precision :: x_TEC(ndim), w_TEC(nw+nwauxio)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP
double precision, dimension(0:nw+nwauxio)                   :: normconv
logical :: fileopen,first
character(len=80) :: filename
integer  :: filenr

!!! possible length conflict
character(len=1024) :: tecplothead

character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead
!-----------------------------------------------------------------------------

if(npe>1)then
 if(mype==0) PRINT *,'tecplot not parallel, use tecplotmpi'
 call mpistop('npe>1, tecplot')
end if

if(nw/=count(writew(1:nw)))then
 if(mype==0) PRINT *,'tecplot does not use writew=F'
 call mpistop('writew, tecplot')
end if

if(nocartesian)then
 if(mype==0) PRINT *,'tecplot with nocartesian and typeaxial=',typeaxial
endif

inquire(qunit,opened=fileopen)
if (.not.fileopen) then
   ! generate filename    
   filenr=snapshotini
   if (.not.convert) filenr=snapshot-1
   write(filename,'(a,i4.4,a)') TRIM(filenameout),filenr,".plt"
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
         ', SOLUTIONTIME=',t*normt,', F=POINT' 

   igonlevel=0
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,ixC^L,ixCC^L,.true.)
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
            ', SOLUTIONTIME=',t*normt,', DATAPACKING=POINT, ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
       do iigrid=1,igridstail; igrid=igrids(iigrid);
         if (node(plevel_,igrid)/=level) cycle
         call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
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
            ', SOLUTIONTIME=',t*normt,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
       else
        if(ndim+nw+nwauxio<10) then
         ! difference only in length of integer format specification for ndim+nw+nwauxio
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i1,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevel,', E=',elemsonlevel, &
            ', SOLUTIONTIME=',t*normt,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,'-',ndim+nw+nwauxio,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
        else
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i2,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevel,', E=',elemsonlevel, &
            ', SOLUTIONTIME=',t*normt,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,'-',ndim+nw+nwauxio,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
        endif
       endif
       do idim=1,ndim
         first=(idim==1) 
         do iigrid=1,igridstail; igrid=igrids(iigrid);
            if (node(plevel_,igrid)/=level) cycle
            call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                           ixC^L,ixCC^L,first)
            write(qunit,fmt="(100(e14.6))") xC_TMP(ixC^S,idim)*normconv(0)
         enddo
       enddo
       do iw=1,nw+nwauxio
         do iigrid=1,igridstail; igrid=igrids(iigrid);
            if (node(plevel_,igrid)/=level) cycle
            call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
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
      igonlevel=igonlevel+1
      call save_conntec(qunit,igrid,igonlevel)
   end do
end do
{^IFONED endif}

end subroutine tecplot
!=============================================================================
subroutine save_conntec(qunit,igrid,igonlevel)

! this saves the basic line, quad and brick connectivity,
! as used by TECPLOT file outputs for unstructured grid

include 'amrvacdef.f'

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
subroutine calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                     ixC^L,ixCC^L,first)

! this subroutine computes both corner as well as cell-centered values
! it handles how we do the center to corner averaging, as well as 
! whether we switch to cartesian or want primitive or conservative output,
! handling the addition of B0 in B0+B1 cases, ...
!
! the normconv is passed on to specialvar_output for extending with
! possible normalization values for the nw+1:nw+nwauxio entries

include 'amrvacdef.f'

integer, intent(in) :: qunit, igrid
logical, intent(in) :: first

integer :: nx^D, nxC^D, ix^D, ix, iw, level, idir
integer :: ixC^L,ixCC^L,nxCC^D
double precision :: dx^D

integer :: idims,jxC^L
double precision :: ldw(ixG^T), dwC(ixG^T)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP
double precision, dimension(ixG^T,1:nw+nwauxio)   :: w

double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv 
integer ::iwe,iwb^C
logical, save :: subfirst=.true.
!-----------------------------------------------------------------------------
! following only for allowing compiler to go through with debug on
{#IFDEF ENERGY
iwe=e_
}
iwb^C=b^C_;


nx^D=ixMhi^D-ixMlo^D+1;
level=node(plevel_,igrid)
dx^D=dx(^D,level);

! for normalization within the code
if(saveprim) then
  normconv(0:nw)=normvar(0:nw)
else
  normconv(0)=normvar(0)
  ! assuming density
  normconv(1)=normvar(1)
  ! assuming momentum=density*velocity
  if (nw>=2) normconv(2:1+^NC)=normvar(1)*normvar(2:1+^NC)
  ! assuming energy/pressure and magnetic field
  if (nw>=2+^NC) normconv(2+^NC:nw)=normvar(2+^NC:nw)
  if (typephys=='hdmdust') then
  ! energy followed by dust density and momentum
     normconv(2+^NC)=normvar(2+^NC) !energy
     normconv(3+^NC:2+^NC+^NDS)=normvar(3+^NC:2+^NC+^NDS) !dust density
     {^DS&{^C&normconv(2+^NC+^NDS+^C+(^DS-1)*^NC)= &
            normvar(2+^NC+^NDS+^C+(^DS-1)*^NC)* &
            normvar(2+^NC+^DS);}\} !dust momentum     
  end if     
end if

! coordinates of cell centers
nxCC^D=nx^D;
ixCCmin^D=ixMlo^D; ixCCmax^D=ixMhi^D;
{do ix=ixCCmin^D,ixCCmax^D
    xCC(ix^D%ixCC^S,^D)=rnode(rpxmin^D_,igrid)+(dble(ix-ixCCmin^D)+half)*dx^D
end do\}

! coordinates of cell corners
nxC^D=nx^D+1;
ixCmin^D=ixMlo^D-1; ixCmax^D=ixMhi^D;
{do ix=ixCmin^D,ixCmax^D
    xC(ix^D%ixC^S,^D)=rnode(rpxmin^D_,igrid)+dble(ix-ixCmin^D)*dx^D
end do\}

! In case primitives to be saved: use primitive subroutine
!  extra layer around mesh only needed when storing corner values and averaging
w(ixG^T,1:nw)=pw(igrid)%w(ixG^T,1:nw)
if(saveprim.and.first) call primitive(ixG^LL,ixM^LL^LADD1,w,px(igrid)%x)

if (nwextra>0) then
 ! here we actually fill the ghost layers for the nwextra variables using 
 ! continuous extrapolation (as these values do not exist normally in ghost cells)
 do idims=1,ndim
  select case(idims)
   {case(^D)
     jxCmin^DD=ixGhi^D+1-dixB^D%jxCmin^DD=ixGlo^DD;
     jxCmax^DD=ixGhi^DD;
     do ix^D=jxCmin^D,jxCmax^D
         w(ix^D^D%jxC^S,nw-nwextra+1:nw) = w(jxCmin^D-1^D%jxC^S,nw-nwextra+1:nw)
     end do 
     jxCmin^DD=ixGlo^DD;
     jxCmax^DD=ixGlo^D-1+dixB^D%jxCmax^DD=ixGhi^DD;
     do ix^D=jxCmin^D,jxCmax^D
         w(ix^D^D%jxC^S,nw-nwextra+1:nw) = w(jxCmax^D+1^D%jxC^S,nw-nwextra+1:nw)
     end do \}
  end select
 end do
end if

! next lines needed when specialvar_output uses gradients
! and later on when dwlimiter2 is used 
typelimiter=typelimiter1(node(plevel_,igrid))
typegradlimiter=typegradlimiter1(node(plevel_,igrid))
^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
if(nwauxio>0)then
  ! auxiliary io variables can be computed and added by user
  ! next few lines ensure correct usage of routines like divvector etc
  if (.not.slab) mygeo => pgeo(igrid)
  if (B0field) then
    myB0_cell => pB0_cell(igrid)
    myB0      => pB0_cell(igrid)
    {^D&myB0_face^D => pB0_face^D(igrid)\}
  end if
  ! default (no) normalization for auxiliary variables
  normconv(nw+1:nw+nwauxio)=one
  ! maybe need for restriction to ixG^LL^LSUB1 ??
  !call specialvar_output(ixG^LL,ixG^LL,w,px(igrid)%x,normconv)
  !call specialvar_output(ixG^LL,ixG^LL^LSUB1,w,px(igrid)%x,normconv)
  call specialvar_output(ixG^LL,ixM^LL^LADD1,w,px(igrid)%x,normconv)
endif

! compute the cell-center values for w first
!===========================================
! cell center values obtained from mere copy, while B0+B1 split handled here
do iw=1,nw+nwauxio
   if (B0field.and.iw>b0_.and.iw<=b0_+ndir) then
      idir=iw-b0_
      {do ix^DB=ixCCmin^DB,ixCCmax^DB\}
         wCC(ix^D,iw)=w(ix^D,iw)+pB0_cell(igrid)%w(ix^D,idir)
      {end do\}
   else
      {do ix^DB=ixCCmin^DB,ixCCmax^DB\}
          wCC(ix^D,iw)=w(ix^D,iw)
      {end do\}
   end if
end do
{#IFDEF ENERGY
if((.not.saveprim) .and. B0field) then
   {do ix^DB=ixCCmin^DB,ixCCmax^DB\}
       wCC(ix^D,iwe)=w(ix^D,iwe) &
           +half*( ^C&pB0_cell(igrid)%w(ix^D,iwb^C-b0_)**2+ ) &
           + ( ^C&w(ix^D,iwb^C)*pB0_cell(igrid)%w(ix^D,iwb^C-b0_)+ )
   {end do\}
endif
}
! compute the corner values for w now by averaging
!=================================================

if(typeaxial=='slab')then
   ! for slab symmetry: no geometrical info required
   do iw=1,nw+nwauxio
      if (B0field.and.iw>b0_.and.iw<=b0_+ndir) then
         idir=iw-b0_
         {do ix^DB=ixCmin^DB,ixCmax^DB\}
           wC(ix^D,iw)=sum(w(ix^D:ix^D+1,iw) &
                           +pB0_cell(igrid)%w(ix^D:ix^D+1,idir))/dble(2**ndim)
         {end do\}
      else
        if(uselimiter)then
           if(ndim>1)call mpistop("to be corrected for multi-D")
           do idims =1,ndim
              jxC^L=ixC^L+kr(idims,^D);
              dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)
              call dwlimiter2(dwC,ixG^LL,ixC^L,iw,idims,ldw,dxlevel(idims))
              wC(ixC^S,iw)=w(ixC^S,iw)+half*ldw(ixC^S)
           end do
        else
          {do ix^DB=ixCmin^DB,ixCmax^DB\}
             wC(ix^D,iw)=sum(w(ix^D:ix^D+1,iw))/dble(2**ndim)
          {end do\}
       end if
      end if
   end do
{#IFDEF ENERGY
   if((.not.saveprim) .and. B0field) then
      {do ix^DB=ixCmin^DB,ixCmax^DB\}
         wC(ix^D,iwe)=sum(w(ix^D:ix^D+1,iwe) &
           +half*( ^C&pB0_cell(igrid)%w(ix^D:ix^D+1,iwb^C-b0_)**2+ ) &
           + ( ^C&w(ix^D:ix^D+1,iwb^C)*pB0_cell(igrid)%w(ix^D:ix^D+1,iwb^C-b0_)+ ) ) &
            /dble(2**ndim)
      {end do\}
   endif
}
else
   do iw=1,nw+nwauxio
      if (B0field.and.iw>b0_.and.iw<=b0_+ndir) then
         idir=iw-b0_
         {do ix^DB=ixCmin^DB,ixCmax^DB\}
           wC(ix^D,iw)= sum((w(ix^D:ix^D+1,iw)+pB0_cell(igrid)%w(ix^D:ix^D+1,idir)) &
                            *pgeo(igrid)%dvolume(ix^D:ix^D+1))    &
                    /sum(pgeo(igrid)%dvolume(ix^D:ix^D+1))
         {end do\}
      else
         {do ix^DB=ixCmin^DB,ixCmax^DB\}
           wC(ix^D,iw)=sum(w(ix^D:ix^D+1,iw)*pgeo(igrid)%dvolume(ix^D:ix^D+1)) &
                    /sum(pgeo(igrid)%dvolume(ix^D:ix^D+1))
         {end do\}
      end if
   end do
{#IFDEF ENERGY
   if((.not.saveprim) .and. B0field) then
      {do ix^DB=ixCmin^DB,ixCmax^DB\}
         wC(ix^D,iwe)=sum((w(ix^D:ix^D+1,iwe) &
           +half*( ^C&pB0_cell(igrid)%w(ix^D:ix^D+1,iwb^C-b0_)**2+ ) &
           + ( ^C&w(ix^D:ix^D+1,iwb^C)*pB0_cell(igrid)%w(ix^D:ix^D+1,iwb^C-b0_)+ ) ) &
                            *pgeo(igrid)%dvolume(ix^D:ix^D+1))    &
                    /sum(pgeo(igrid)%dvolume(ix^D:ix^D+1))
      {end do\}
   endif
}
endif

if(nocartesian) then
  ! keep the coordinate and vector components
  xC_TMP(ixC^S,1:ndim)          = xC(ixC^S,1:ndim)
  wC_TMP(ixC^S,1:nw+nwauxio)    = wC(ixC^S,1:nw+nwauxio)
  xCC_TMP(ixCC^S,1:ndim)        = xCC(ixCC^S,1:ndim)
  wCC_TMP(ixCC^S,1:nw+nwauxio)  = wCC(ixCC^S,1:nw+nwauxio)
else
  ! do all conversions to cartesian coordinates and vector components
  ! start for the corner values
  call cartesian(xC_TMP,wC_TMP,ixC^L,xC,wC)
  ! then cell center values
  call cartesian(xCC_TMP,wCC_TMP,ixCC^L,xCC,wCC)
endif

! Warning: differentiate between idl/idlCC/tecplot/tecplotCC/vtu(B)/vtu(B)CC
if(nwaux>0 .and. mype==0 .and. first.and.subfirst) then
  ! when corner values are computed and auxiliaries present: warn!
  if(convert_type=='idl'.or.convert_type=='tecplot' &
   .or.convert_type=='vtu'.or.convert_type=='vtuB') &
      write(*,*) 'Warning: also averaged auxiliaries within calc_grid'
  subfirst=.false.
endif

end subroutine calc_grid
!=============================================================================
subroutine cartesian(x_TMP,w_TMP,ix^L,xC,wC)

! conversion of coordinate and vector components from cylindrical/spherical
! to cartesian coordinates and components done here
! Also: nullifying values lower than smalldouble

include 'amrvacdef.f'

integer :: ix^L, ix^D, idim, iw, ivector, iw0
integer, dimension(nw) :: vectoriw
double precision :: x_TEC(ndim), w_TEC(nw+nwauxio)
double precision, dimension(ndim,ndim) :: normal

double precision, dimension(ix^S,ndim) :: xC
double precision, dimension(ix^S,nw+nwauxio)   :: wC

double precision, dimension(ix^S,ndim) :: x_TMP
double precision, dimension(ix^S,nw+nwauxio)   :: w_TMP
!-----------------------------------------------------------------------------

vectoriw=-1
if(nvector>0) then
  do ivector=1,nvector
     do idim=1,ndim
        vectoriw(iw_vector(ivector)+idim)=iw_vector(ivector)
     end do
  end do
endif

{do ix^DB=ixmin^DB,ixmax^DB\}
   select case (typeaxial)
   case ("slab","slabtest")
      x_TEC(1:ndim)=xC(ix^D,1:ndim)
      w_TEC(1:nw+nwauxio)=wC(ix^D,1:nw+nwauxio)
   case ("cylindrical")
      {^IFONED
      x_TEC(1)=xC(ix^D,1)}
      {^IFTWOD
      select case (phi_)
      case (2) 
         x_TEC(1)=xC(ix^D,1)*cos(xC(ix^D,2))
         x_TEC(2)=xC(ix^D,1)*sin(xC(ix^D,2))
      case default
         x_TEC(1)=xC(ix^D,1)
         x_TEC(2)=xC(ix^D,2)
      end select}
      {^IFTHREED
      x_TEC(1)=xC(ix^D,1)*cos(xC(ix^D,pphi_))
      x_TEC(2)=xC(ix^D,1)*sin(xC(ix^D,pphi_))
      x_TEC(3)=xC(ix^D,zz_)}

      if (nvector>0) then
         {^IFONED normal(1,1)=one}

         {^IFTWOD
         select case (phi_)
         case (2) 
            normal(1,1)=cos(xC(ix^D,2))
            normal(1,2)=-sin(xC(ix^D,2))
            normal(2,1)=sin(xC(ix^D,2))
            normal(2,2)=cos(xC(ix^D,2))
         case default
            normal(1,1)=one
            normal(1,2)=zero
            normal(2,1)=zero
            normal(2,2)=one
         end select}

         {^IFTHREED
         normal(1,1)=cos(xC(ix^D,pphi_))
         normal(1,pphi_)=-sin(xC(ix^D,pphi_))
         normal(1,zz_)=zero

         normal(2,1)=sin(xC(ix^D,pphi_))
         normal(2,pphi_)=cos(xC(ix^D,pphi_))
         normal(2,zz_)=zero

         normal(3,1)=zero
         normal(3,pphi_)=zero
         normal(3,zz_)=one}
      end if
      do iw=1,nw+nwauxio
         if (iw<=nw) iw0=vectoriw(iw)
         if (iw0>=0.and.iw<=iw0+ndim.and.iw<=nw) then
            idim=iw-iw0
            w_TEC(iw0+idim)={^D&wC(ix^DD,iw0+^D)*normal(idim,^D)+}
         else
            w_TEC(iw)=wC(ix^D,iw)
         end if
      end do
   case ("spherical")
      x_TEC(1)=xC(ix^D,1){^NOONED*sin(xC(ix^D,2))}{^IFTHREED*cos(xC(ix^D,3))}
      {^IFTWOD
      x_TEC(2)=xC(ix^D,1)*cos(xC(ix^D,2))}
      {^IFTHREED
      x_TEC(2)=xC(ix^D,1)*sin(xC(ix^D,2))*sin(xC(ix^D,3))
      x_TEC(3)=xC(ix^D,1)*cos(xC(ix^D,2))}

      if (nvector>0) then
         {^IFONED normal(1,1)=one}
         {^NOONED
         normal(1,1)=sin(xC(ix^D,2)){^IFTHREED*cos(xC(ix^D,3))}
         normal(1,2)=cos(xC(ix^D,2)){^IFTHREED*cos(xC(ix^D,3))
         normal(1,3)=-sin(xC(ix^D,3))}}

         {^IFTWOD
         normal(2,1)=cos(xC(ix^D,2))
         normal(2,2)=-sin(xC(ix^D,2))}
         {^IFTHREED
         normal(2,1)=sin(xC(ix^D,2))*sin(xC(ix^D,3))
         normal(2,2)=cos(xC(ix^D,2))*sin(xC(ix^D,3))
         normal(2,3)=cos(xC(ix^D,3))

         normal(3,1)=cos(xC(ix^D,2))
         normal(3,2)=-sin(xC(ix^D,2))
         normal(3,3)=zero}
      end if
      do iw=1,nw+nwauxio
         if (iw<=nw) iw0=vectoriw(iw)
         if (iw0>=0.and.iw<=iw0+ndim.and.iw<=nw) then
            idim=iw-iw0
            w_TEC(iw0+idim)={^D&wC(ix^DD,iw0+^D)*normal(idim,^D)+}
         else
            w_TEC(iw)=wC(ix^D,iw)
         end if
      end do
   case default
      write(*,*) "No converter for typeaxial=",typeaxial
   end select
   x_TMP(ix^D,1:ndim)=x_TEC(1:ndim)
   w_TMP(ix^D,1:nw+nwauxio)=w_TEC(1:nw+nwauxio)
   ! Be aware that small values are nullified here!!!
   where(dabs(w_TMP(ix^D,1:nw+nwauxio))<smalldouble)
         w_TMP(ix^D,1:nw+nwauxio)=zero
   endwhere
{end do\}

end subroutine cartesian
!=============================================================================
subroutine unstructuredvtk(qunit)

! output for vtu format to paraview
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using normvar-array

include 'amrvacdef.f'

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP
double precision, dimension(0:nw+nwauxio)                   :: normconv
integer::               igrid,iigrid,level,igonlevel,icel,ixC^L,ixCC^L,iw
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx^D,nxC^D,nodesonlevel,elemsonlevel,nc,np,VTK_type,ix^D

character(len=80)::  filename
integer          :: filenr

character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
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
   if (.not.convert) filenr=snapshot-1
  write(filename,'(a,i4.4,a)') TRIM(filenameout),filenr,".vtu"
  ! Open the file for the header part
  open(qunit,file=filename,status='unknown')
endif

call getheadernames(wnamei,xandwnamei,outfilehead)

! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
{#IFDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="BigEndian">'}
{#IFNDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'}
write(qunit,'(a)')'<UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(qunit,*) dble(dble(t)*normt)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'

! number of cells, number of corner points, per grid.
nx^D=ixMhi^D-ixMlo^D+1;
nxC^D=nx^D+1;
nc={nx^D*}
np={nxC^D*}

! Note: using the writew, writelevel, writespshift
! we can clip parts of the grid away, select variables, levels etc.
do level=levmin,levmax
 if (writelevel(level)) then
   do iigrid=1,igridstail; igrid=igrids(iigrid);
    if (node(plevel_,igrid)/=level) cycle
    ! only output a grid when fully within clipped region selected
    ! by writespshift array
    if (({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
          *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
         <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})) then
      call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                     ixC^L,ixCC^L,.true.)
      select case(convert_type)
       case('vtu')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') &
            '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<PointData>'
         do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.writew(iw)) cycle
         endif

            write(qunit,'(a,a,a)')&
          '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') {(|}wC_TMP(ix^D,iw)*normconv(iw),{ix^D=ixCmin^D,ixCmax^D)}
            write(qunit,'(a)')'</DataArray>'
         enddo
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
            if(.not.writew(iw)) cycle
         endif

            write(qunit,'(a,a,a)')&
          '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') {(|}wCC_TMP(ix^D,iw)*normconv(iw),{ix^D=ixCCmin^D,ixCCmax^D)}
            write(qunit,'(a)')'</DataArray>'
         enddo
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
      enddo
      write(qunit,'(a)')'</DataArray>'

      write(qunit,'(a)')'</Cells>'

      write(qunit,'(a)')'</Piece>'
    endif
   enddo
 endif
enddo

write(qunit,'(a)')'</UnstructuredGrid>'
write(qunit,'(a)')'</VTKFile>'
close(qunit)

end subroutine unstructuredvtk
!====================================================================================
subroutine unstructuredvtkB(qunit)

! output for vtu format to paraview, binary version output
! not parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using normvar-array
use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
include 'amrvacdef.f'

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio):: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)  :: wCC_TMP
double precision :: normconv(0:nw+nwauxio)

integer, allocatable :: intstatus(:,:)
integer :: ipe,igrid,level,icel,ixC^L,ixCC^L,Morton_no,Morton_length
integer :: nx^D,nxC^D,nc,np,VTK_type,ix^D,filenr
integer*8 :: offset

integer::  size_int,size_double,size_length,k,iw
integer::  length,lengthcc,length_coords,length_conn,length_offsets
character::  buf
character(len=80)::  filename
character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical ::   fileopen,cell_corner
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
   call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                     ixC^L,ixCC^L,.true.)
   itag=Morton_no
   call MPI_SEND(xC_TMP,1,type_block_xc_io, 0,itag,icomm,ierrmpi)
   if(cell_corner) then
     call MPI_SEND(wC_TMP,1,type_block_wc_io, 0,itag,icomm,ierrmpi)
   else
     call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
   endif
 end do

else
 ! mype==0
 offset=0
 !size_double=8
 size_double=4
 size_length=4
 size_int=size_length
 
 inquire(qunit,opened=fileopen)
 if(.not.fileopen)then
   ! generate filename 
    filenr=snapshotini
    if (.not.convert) filenr=snapshot-1
   write(filename,'(a,i4.4,a)') TRIM(filenameout),filenr,".vtu"
   ! Open the file for the header part
   open(qunit,file=filename,status='replace')
 endif
 
 call getheadernames(wnamei,xandwnamei,outfilehead)
 
 ! generate xml header
 write(qunit,'(a)')'<?xml version="1.0"?>'
 write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
 {#IFDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="BigEndian">'}
 {#IFNDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'}
 write(qunit,'(a)')'<UnstructuredGrid>'
 write(qunit,'(a)')'<FieldData>'
 write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                    'NumberOfTuples="1" format="ascii">'
 write(qunit,*) real(t*normt)
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

 ! Note: using the writew, writelevel, writespshift
 do Morton_no=Morton_start(0),Morton_stop(0)
    if(.not. Morton_aim(Morton_no)) cycle
    if(cell_corner) then
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') &
         '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<PointData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.writew(iw)) cycle
         endif

         write(qunit,'(a,a,a,i16,a)')&
             '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
             '" format="appended" offset="',offset,'">'
         write(qunit,'(a)')'</DataArray>'
         offset=offset+length+size_length
      enddo
      write(qunit,'(a)')'</PointData>'

      write(qunit,'(a)')'<Points>'
      write(qunit,'(a,i16,a)') &
   '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      offset=offset+length_coords+size_length
      write(qunit,'(a)')'</Points>'
    else
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') &
         '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<CellData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.writew(iw)) cycle
         endif

         write(qunit,'(a,a,a,i16,a)')&
             '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
             '" format="appended" offset="',offset,'">'
         write(qunit,'(a)')'</DataArray>'
         offset=offset+lengthcc+size_length
      enddo
      write(qunit,'(a)')'</CellData>'

      write(qunit,'(a)')'<Points>'
      write(qunit,'(a,i16,a)') &
   '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      offset=offset+length_coords+size_length
      write(qunit,'(a)')'</Points>'
    end if
   
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
            if(.not.writew(iw)) cycle
         endif

           write(qunit,'(a,a,a,i16,a)')&
               '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
               '" format="appended" offset="',offset,'">'
           write(qunit,'(a)')'</DataArray>'
           offset=offset+length+size_length
        enddo
        write(qunit,'(a)')'</PointData>'

        write(qunit,'(a)')'<Points>'
        write(qunit,'(a,i16,a)') &
     '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
        ! write cell corner coordinates in a backward dimensional loop, always 3D output
        offset=offset+length_coords+size_length
        write(qunit,'(a)')'</Points>'
      else
        ! we write out every grid as one VTK PIECE
        write(qunit,'(a,i7,a,i7,a)') &
           '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
        write(qunit,'(a)')'<CellData>'
        do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.writew(iw)) cycle
         endif

           write(qunit,'(a,a,a,i16,a)')&
               '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
               '" format="appended" offset="',offset,'">'
           write(qunit,'(a)')'</DataArray>'
           offset=offset+lengthcc+size_length
        enddo
        write(qunit,'(a)')'</CellData>'

        write(qunit,'(a)')'<Points>'
        write(qunit,'(a,i16,a)') &
     '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
        ! write cell corner coordinates in a backward dimensional loop, always 3D output
        offset=offset+length_coords+size_length
        write(qunit,'(a)')'</Points>'
      end if
     
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
   call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                  ixC^L,ixCC^L,.true.)
   do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.writew(iw)) cycle
         endif
     if(cell_corner) then
       write(qunit) length
       write(qunit) {(|}real(wC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCmin^D,ixCmax^D)}
     else
       write(qunit) lengthcc
       write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
     end if
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
            if(.not.writew(iw)) cycle
         endif
        if(cell_corner) then
          write(qunit) length
          write(qunit) {(|}real(wC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCmin^D,ixCmax^D)}
        else
          write(qunit) lengthcc
          write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
        end if
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
endif

end subroutine unstructuredvtkB
!====================================================================================
subroutine save_connvtk(qunit,igrid)

! this saves the basic line, pixel and voxel connectivity,
! as used by VTK file outputs for unstructured grid

include 'amrvacdef.f'

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
subroutine valout_dx(qunit)

!
!  Numberings in DX start at zero.
!  Array ordering becomes row-major (C/DX style).
!
include 'amrvacdef.f'

integer, intent(in) :: qunit

integer :: iigrid, igrid, level, ngrids, nx^D

integer,parameter::     size_double = 8
integer,parameter::     size_byte   = 1
integer,parameter::     size_recsep = 4
character(len=5) ::     byteorder

character(len=80)::     filename
character(len=80)::     name,physics,scanstring,wname
integer          ::     filenr

integer::               underscore_position
integer::               iw, space_position, max_name_len
integer::               offset
integer::               nummeshpoints
integer::               firstgridonlevel,lastgridonlevel
integer::               NumGridsOnLevel(1:nlevelshi)

integer,parameter ::    byte=selected_int_kind(1)

character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

character(LEN=10)    :: dummy_date,dummy_time,dummy_zone
integer,dimension(8) :: DateAndTime
!-----------------------------------------------------------------------------
if(npe>1)then
 if(mype==0) PRINT *,'valoutdx not parallel'
 call mpistop('npe>1, valoutdx')
end if

if(nw/=count(writew(1:nw)))then
 if(mype==0) PRINT *,'valoutdx does not use writew=F'
 call mpistop('writew, valoutdx')
end if

if(saveprim)then
 if(mype==0) PRINT *,'valoutdx does not use saveprim=T'
 call mpistop('saveprim, valoutdx')
end if

if(nwauxio>0)then
 if(mype==0) PRINT *,'valoutdx does not use nwauxio>0'
 call mpistop('nwauxio>0, valoutdx')
end if

if (B0field) call mpistop("No B0 field implemented in dx plotfile")

nx^D=ixGhi^D-2*dixB;

byteorder = ' '//TRIM(dxfiletype)//' '
   ! generate filename    
   filenr=snapshotini
   if (.not.convert) filenr=snapshot-1
   write(filename,'(a,i4.4,a)') TRIM(filenameout),filenr,".plt"

call date_and_time(dummy_date,dummy_time,dummy_zone,DateAndTime)

! Open the file for the header part, ascii data
open(qunit,file=filename,status='unknown',form='formatted')

write(qunit,'(2a)') '### AMRVAC datafile for simulation ',TRIM(fileheadout)
write(qunit,'(a,i02,a,i02,a,i4,a,i02,a,i02)') '### Generated on ', &
     DateAndTime(3),'/',DateAndTime(2),'/',DateAndTime(1), &
     ' at ',DateAndTime(5),'h',DateAndTime(6)

call getheadernames(wnamei,xandwnamei,outfilehead)

offset = 0
ngrids = 0
do level=levmin,levmax
   NumGridsOnLevel(level)=0
   firstgridonlevel = ngrids
   write(qunit,'(a,i3.3)') '# start level',level
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (node(plevel_,igrid)/=level) cycle
      NumGridsOnLevel(level)=NumGridsOnLevel(level)+1
      ngrids = ngrids+1
      write(qunit,'(a,i5.5)') '# start grid ',ngrids-1
      !
      ! write positions object
      !
      write(qunit,'(a,i5.5,a,3(x,i11))') 'object "pos',ngrids-1, &
                  '" class gridpositions counts ',{nx^D+1}
      write(qunit,'(a,3(x,f25.16))') '  origin ',{rnode(rpxmin^D_,igrid)*normvar(0)}
      {^IFONED 
      write(qunit,'(a,x,f25.16)') '  delta  ',rnode(rpdx1_,igrid)*normvar(0)
      \}
      {^IFTWOD
      write(qunit,'(a,2(x,f25.16))') '  delta  ',rnode(rpdx1_,igrid)*normvar(0),0.0d0
      write(qunit,'(a,2(x,f25.16))') '  delta  ',0.0d0,rnode(rpdx2_,igrid)*normvar(0)
      \}
      {^IFTHREED 
      write(qunit,'(a,3(x,f25.16))') '  delta  ',rnode(rpdx1_,igrid)*normvar(0),0.0d0,0.0d0
      write(qunit,'(a,3(x,f25.16))') '  delta  ',0.0d0,rnode(rpdx2_,igrid)*normvar(0),0.0d0
      write(qunit,'(a,3(x,f25.16))') '  delta  ',0.0d0,0.0d0,rnode(rpdx3_,igrid)*normvar(0)
      \}
      write(qunit,'(a)') 'attribute "dep" string "positions" '
      write(qunit,'(a)') '#'
      !
      ! write connections object
      !
      write(qunit,'(a,i5.5,a,3(x,i11))') 'object "con',ngrids-1, &
                  '" class gridconnections counts ',{nx^D+1}
      {^IFONED 
        write(qunit,'(a)') 'attribute "element type" string "lines" '
      \}
      {^IFTWOD 
        write(qunit,'(a)') 'attribute "element type" string "quads" '
      \}
      {^IFTHREED 
        write(qunit,'(a)') 'attribute "element type" string "cubes" '
      \}
      write(qunit,'(a)') 'attribute "ref" string "positions" '
      write(qunit,'(a)') '#'
      !      
      ! write data object (header info)
      !
      nummeshpoints={nx^D*}
      offset = offset + size_recsep

      write(qunit,'(a,i5.5,a,i3,a,x,i11,2a,x,i11)') &
           'object "dat',                            ngrids-1, &
           '" class array type double rank 1 shape ',nw+nwauxio, &
           ' items ',                                nummeshpoints, &
           byteorder, 'binary data ',                offset
      offset = offset + (nw+nwauxio)*nummeshpoints*size_double + size_recsep
      write(qunit,'(a)') 'attribute "dep" string "connections" '
      write(qunit,'(a)') '#'
      !
      ! write field object
      !
      write(qunit,'(a,i5.5,a)') 'object ',ngrids-1,' class field'
      write(qunit,'(a,i5.5,a)') '  component "positions" value "pos', &
                                ngrids-1,'"'
      write(qunit,'(a,i5.5,a)') '  component "connections" value "con', &
                                ngrids-1,'"'
      write(qunit,'(a,i5.5,a)') '  component "data" value "dat', &
                                ngrids-1,'"'
      write(qunit,'(a,i5.5)') '  attribute "refinement level" number ',level
      write(qunit,'(a,i5.5)') '# end grid ',ngrids-1
   end do

   lastgridonlevel=ngrids-1

   write(qunit,'(a)') '#'
   write(qunit,'(a,i3.3,a,i5.5,a)') '# end level',level, &
                           ' (',NumGridsOnLevel(level),' grids)'
end do
write(qunit,'(a)') '#'
!
! eqpar array
!
write(qunit,'(a,x,i11,a)') &
  'object "eqpararray" class array type float items ', &
  neqpar+nspecialpar,' data follows'
write(qunit,'(f24.12)') eqpar
write(qunit,'(a)') '#'
!
! # grids on level
!
write(qunit,'(a,x,i11,a)') &
  'object "ngridsonlevarray" class array type int items ', &
  levmax-levmin+1,' data follows'
write(qunit,'(100(x,i11))') NumGridsOnLevel(levmin:levmax)
write(qunit,'(a)') '#'
!
! wnames array
!
write(qunit,'(a,x,i11,a)') &
  'object "wnamesarray" class array type string rank 1 shape 80 items ', &
  nw+nwauxio,' data follows'
do iw=1,nw+nwauxio
   write(qunit,'(a,a,a)', advance='no') ' "',TRIM(wnamei(iw)),'"'
enddo
write(qunit,'(a)') ' '
write(qunit,'(a)') '#'

! Separate name and physics from fileheadout
underscore_position = index(TRIM(fileheadout),'_',.true.)
if (underscore_position == 0) then
   name=fileheadout
   physics='unknown'
else
   name=fileheadout(:underscore_position-1)
   physics=fileheadout(underscore_position+1:)
endif
!
! Top level group with all attributes
!
write(qunit,'(a)') 'object "default" class multigrid'
do igrid=1,ngrids
   write(qunit,'(a,i5,a,i5.5,a)') 'member ',igrid-1,' value ',igrid-1
end do
write(qunit,'(3a)')  'attribute "simulationname"     string "',TRIM(name),'"'
write(qunit,'(3a)')  'attribute "physics"  string "',TRIM(physics),'"'
write(qunit,'(a,x,i11)') 'attribute "ndim"     number ',ndim
write(qunit,'(a,x,i11)') 'attribute "ndir"     number ',ndir
write(qunit,'(a,x,i11)') 'attribute "nw"       number ',nw+nwauxio
write(qunit,'(a,x,i11)') 'attribute "timestep" number ',it
write(qunit,'(a,f25.16)')'attribute "time"     number ',t *normt
write(qunit,'(a)')   'attribute "eqpar"    value "eqpararray"'
write(qunit,'(a)')   'attribute "ngrids"   value "ngridsonlevarray"'
write(qunit,'(a)')   'attribute "wnames"   value "wnamesarray"'
write(qunit,'(a)') '#'

! denote the end of the header section
write(qunit,'(a,i5.5)') '# end header section'
write(qunit,'(a)') 'end'

close(qunit)
! now for the binary part...
open(qunit,file=filename,status='unknown', &
           form='unformatted',position='append')

ngrids=0
offset = 0
do level=levmin,levmax
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (node(plevel_,igrid)/=level) cycle
      ngrids = ngrids+1
      nummeshpoints={nx^D*}
      offset = offset + size_recsep
      ! write data array
      call varout_dx_condep(qunit,pw(igrid)%w,ixG^LL)
      offset = offset + (nw+nwauxio)*nummeshpoints*size_double + size_recsep
   end do
end do

close(qunit)

end subroutine valout_dx
!============================================================================
subroutine varout_dx_condep(qunit,w,ixG^L)

include 'amrvacdef.f'

integer, intent(in) :: qunit, ixG^L
double precision, intent(in) :: w(ixG^S,1:nw)

integer :: ixM^L, ix^D, iw
!-----------------------------------------------------------------------------
ixM^L=ixG^L^LSUBdixB;

! We write the arrays in row-major order (C/DX style) for the spatial indices
write(qunit) ({(|}w(ix^D,iw)*normvar(iw),iw=1,nw),{ix^DB=ixMmin^DB,ixMmax^DB)}

end subroutine varout_dx_condep


!============================================================================
subroutine ImageDataVtk_mpi(qunit)

! output for vti format to paraview, non-binary version output
! parallel, uses calc_grid to compute nwauxio variables
! allows renormalizing using normvar-array
! allows skipping of writew selected variables

! implementation such that length of ASCII output is identical when 
! run on 1 versus multiple CPUs (however, the order of the vtu pieces can differ)

use mod_forest, only: Morton_start, Morton_stop, tree_node_ptr, igrid_to_node, sfc_to_igrid
include 'amrvacdef.f'

integer, intent(in) ::    qunit

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP,xC_TMP_recv
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP,xCC_TMP_recv

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

character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical :: fileopen
integer :: ipe,Morton_no,Morton_length
integer :: ixrvC^L, ixrvCC^L, siz_ind, ind_send(5*^ND), ind_recv(5*^ND)
double precision    :: origin(1:3), spacing(1:3)
integer :: wholeExtent(1:6), ig^D
type(tree_node_ptr) :: tree
!-----------------------------------------------------------------------------
if(levmin/=levmax) call mpistop('ImageData can only be used when levmin=levmax')

normconv(0:nw)=normvar(0:nw)
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
   call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
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
    if (.not.convert) filenr=snapshot-1
    write(filename,'(a,i4.4,a)') TRIM(filenameout),filenr,".vti"
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
{#IFDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="BigEndian">'}
{#IFNDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'}
write(qunit,'(a,3(1pe14.6),a,6(i10),a,3(1pe14.6),a)')'  <ImageData Origin="',origin,'" WholeExtent="',wholeExtent,'" Spacing="',spacing,'">'
 write(qunit,'(a)')'<FieldData>'
 write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                    'NumberOfTuples="1" format="ascii">'
 write(qunit,*) real(t*normt)
 write(qunit,'(a)')'</DataArray>'
 write(qunit,'(a)')'</FieldData>'

! write the data from proc 0
do Morton_no=Morton_start(0),Morton_stop(0)
   if(.not. Morton_aim(Morton_no)) cycle
   igrid=sfc_to_igrid(Morton_no)
   tree%node => igrid_to_node(igrid, 0)%node
   {^D& ig^D = tree%node%ig^D; }
   call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
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
include 'amrvacdef.f'

integer, intent(in) ::    qunit
!
double precision, dimension(0:nw+nwauxio)                   :: normconv
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim)         :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)           :: xCC_TMP
double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP
character(len=10)   :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
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
   if (.not.convert) filenr=snapshot-1
   ! Open the file for the header part
   write(pfilename,'(a,i4.4,a,i4.4,a)') TRIM(filenameout),filenr,"p",mype,".vtu"
   open(qunit,file=pfilename,status='unknown',form='formatted')
endif
! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
{#IFDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="BigEndian">'}
{#IFNDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'}
write(qunit,'(a)')'  <UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(qunit,*) real(t*normt)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'

call getheadernames(wnamei,xandwnamei,outfilehead)

! number of cells, number of corner points, per grid.
nx^D=ixMhi^D-ixMlo^D+1;
nxC^D=nx^D+1;
nc={nx^D*}
np={nxC^D*}

! Note: using the writew, writelevel, writespshift
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

    call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
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
! allows renormalizing using normvar-array
! allows skipping of writew selected variables

! implementation such that length of ASCII output is identical when 
! run on 1 versus multiple CPUs (however, the order of the vtu pieces can differ)

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
include 'amrvacdef.f'

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP,xC_TMP_recv
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP,xCC_TMP_recv

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP,wC_TMP_recv
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP,wCC_TMP_recv
double precision, dimension(0:nw+nwauxio)                   :: normconv
integer::               igrid,iigrid,level,ixC^L,ixCC^L
integer::               NumGridsOnLevel(1:nlevelshi)
integer :: nx^D,nxC^D,nodesonlevel,elemsonlevel,nc,np,ix^D

character(len=80)::  filename
integer ::           filenr

integer, allocatable :: intstatus(:,:)

character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical :: fileopen,conv_grid,cond_grid_recv
integer :: ipe,Morton_no,siz_ind
integer :: ind_send(4*^ND),ind_recv(4*^ND)
integer :: levmin_recv,levmax_recv,level_recv,igrid_recv,ixrvC^L,ixrvCC^L
!-----------------------------------------------------------------------------
if (mype==0) then
 inquire(qunit,opened=fileopen)
 if(.not.fileopen)then
    ! generate filename 
    filenr=snapshotini
    if (.not.convert) filenr=snapshot-1
    write(filename,'(a,i4.4,a)') TRIM(filenameout),filenr,".vtu"
   ! Open the file for the header part
   open(qunit,file=filename,status='unknown',form='formatted')
 endif
 ! generate xml header
 write(qunit,'(a)')'<?xml version="1.0"?>'
 write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
{#IFDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="BigEndian">'}
{#IFNDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'}
 write(qunit,'(a)')'<UnstructuredGrid>'
 write(qunit,'(a)')'<FieldData>'
 write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                    'NumberOfTuples="1" format="ascii">'
 write(qunit,*) real(t*normt)
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


! Note: using the writew, writelevel, writespshift
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

    call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
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

include 'amrvacdef.f'

integer, intent(in) :: qunit
integer, intent(in) :: ixI^L,ixC^L,ixCC^L
integer, intent(in) :: igrid,nc,np,nx^D,nxC^D
double precision, intent(in) :: normconv(0:nw+nwauxio) 
character(len=10), intent(in)::  wnamei(1:nw+nwauxio)

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
            if(.not.writew(iw)) cycle
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
            if(.not.writew(iw)) cycle
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
include 'amrvacdef.f'

integer, intent(in) :: qunit
integer, intent(in) :: ixI^L,ixC^L,ixCC^L
integer, intent(in) :: ig^D,nx^D
double precision, intent(in) :: normconv(0:nw+nwauxio) 
character(len=10), intent(in)::  wnamei(1:nw+nwauxio)

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
            if(.not.writew(iw)) cycle
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
            if(.not.writew(iw)) cycle
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

include 'amrvacdef.f'

integer, intent(in) :: qunit

character(len=10)   :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio),outtype
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
   if (.not.convert) filenr=snapshot-1
   write(filename,'(a,i4.4,a)') TRIM(filenameout),filenr,".pvtu"
   ! Open the file
   open(qunit,file=filename,status='unknown',form='formatted')
endif

call getheadernames(wnamei,xandwnamei,outfilehead)

! Get the default selection:
iscalars=1
do iw=nw,1, -1
   if (writew(iw)) iscalars=iw
end do


! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="PUnstructuredGrid"'
{#IFDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="BigEndian">'}
{#IFNDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'}
write(qunit,'(a)')'  <PUnstructuredGrid GhostLevel="0">'
! Either celldata or pointdata:
write(qunit,'(a,a,a,a,a)')&
     '    <',TRIM(outtype),' Scalars="',TRIM(wnamei(iscalars))//'">'
do iw=1,nw
   if(.not.writew(iw))cycle
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
   write(pfilename,'(a,i4.4,a,i4.4,a)') TRIM(filenameout(&
        INDEX (filenameout, '/', BACK = .TRUE.)+1:&
        LEN(filenameout))),filenr,"p",&
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
! allows renormalizing using normt and normvar-array

! the current implementation is such that tecplotmpi and tecplotCCmpi will 
! create different length output ASCII files when used on 1 versus multiple CPUs
! in fact, on 1 CPU, there will be as many zones as there are levels
! on multiple CPUs, there will be a number of zones up to the number of
! levels times the number of CPUs (can be less, when some level not on a CPU)

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
include 'amrvacdef.f'

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

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP,wC_TMP_recv
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP,wCC_TMP_recv
double precision, dimension(0:nw+nwauxio)                   :: normconv
logical :: fileopen,first
integer :: Morton_no,ipe,levmin_recv,levmax_recv,igrid_recv,level_recv
integer :: ixrvC^L,ixrvCC^L
integer :: ind_send(2*^ND),ind_recv(2*^ND),siz_ind,igonlevel_recv
integer :: NumGridsOnLevel_mype(1:nlevelshi,0:npe-1)
character(len=80) :: filename
integer ::           filenr
character(len=1024) :: tecplothead

character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead
!-----------------------------------------------------------------------------
if(nw/=count(writew(1:nw)))then
 if(mype==0) PRINT *,'tecplot_mpi does not use writew=F'
 call mpistop('writew, tecplot')
end if

if(nocartesian)then
 if(mype==0) PRINT *,'tecplot_mpi with nocartesian and typeaxial=',typeaxial
endif

Master_cpu_open : if (mype == 0) then
 inquire(qunit,opened=fileopen)
 if (.not.fileopen) then
   ! generate filename
    filenr=snapshotini
    if (.not.convert) filenr=snapshot-1
   write(filename,'(a,i4.4,a)') TRIM(filenameout),filenr,".plt"
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
         ', SOLUTIONTIME=',t*normt,', F=POINT' 

   igonlevel=0
   do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid = sfc_to_igrid(Morton_no)
      call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,ixC^L,ixCC^L,.true.)
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
             ', SOLUTIONTIME=',t*normt,', DATAPACKING=POINT, ZONETYPE=', &
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
         call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
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
            ', SOLUTIONTIME=',t*normt,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
       else
        if(ndim+nw+nwauxio<10) then
         ! difference only in length of integer format specification for ndim+nw+nwauxio
         if (mype==0.and.(nodesonlevelmype>0.and.elemsonlevelmype>0))&
          write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i1,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
            ', SOLUTIONTIME=',t*normt,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,'-',ndim+nw+nwauxio,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
        else
         if (mype==0.and.(nodesonlevelmype>0.and.elemsonlevelmype>0))&
          write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i2,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
            ', SOLUTIONTIME=',t*normt,', DATAPACKING=BLOCK, VARLOCATION=([', &
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

          call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
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

         call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
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
             ', SOLUTIONTIME=',t*normt,', DATAPACKING=POINT, ZONETYPE=', &
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
            ', SOLUTIONTIME=',t*normt,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
       else
        if(ndim+nw+nwauxio<10) then
         ! difference only in length of integer format specification for ndim+nw+nwauxio
         if(nodesonlevelmype>0.and.elemsonlevelmype>0) &
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i1,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
            ', SOLUTIONTIME=',t*normt,', DATAPACKING=BLOCK, VARLOCATION=([', &
            ndim+1,'-',ndim+nw+nwauxio,']=CELLCENTERED), ZONETYPE=', &
         {^IFONED 'FELINESEG'}{^IFTWOD 'FEQUADRILATERAL'}{^IFTHREED 'FEBRICK'}
        else
         if(nodesonlevelmype>0.and.elemsonlevelmype>0) &
         write(qunit,"(a,i7,a,a,i7,a,i7,a,f25.16,a,i1,a,i2,a,a)") &
            'ZONE T="',level,'"',', N=',nodesonlevelmype,', E=',elemsonlevelmype, &
            ', SOLUTIONTIME=',t*normt,', DATAPACKING=BLOCK, VARLOCATION=([', &
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
! allows renormalizing using normvar-array

use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
include 'amrvacdef.f'

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP

integer :: igrid,iigrid,level,igonlevel,icel,ixC^L,ixCC^L,Morton_no
integer ::               NumGridsOnLevel(1:nlevelshi)
integer :: nx^D,nxC^D,nodesonlevel,elemsonlevel,nc,np,VTK_type,ix^D
double precision :: normconv(0:nw+nwauxio)
character(len=80) :: pfilename
character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

integer*8 :: offset
integer::  size_int,size_double,size_length,recsep,k,iw,filenr
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
   filenr=snapshotini
   if (.not.convert) filenr=snapshot-1
   ! Open the file for the header part
   write(pfilename,'(a,i4.4,a,i4.4,a)') TRIM(filenameout),filenr,"p",mype,".vtu"
   open(qunit,file=pfilename,status='unknown',form='formatted')
endif
! generate xml header
write(qunit,'(a)')'<?xml version="1.0"?>'
write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
{#IFDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="BigEndian">'}
{#IFNDEF BIGENDIAN write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'}
write(qunit,'(a)')'  <UnstructuredGrid>'
write(qunit,'(a)')'<FieldData>'
write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                   'NumberOfTuples="1" format="ascii">'
write(qunit,*) real(t*normt)
write(qunit,'(a)')'</DataArray>'
write(qunit,'(a)')'</FieldData>'


offset=0
recsep=4
size_double=4
size_length=4
size_int=size_length

call getheadernames(wnamei,xandwnamei,outfilehead)

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

! Note: using the writew, writelevel, writespshift
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
            if(.not.writew(iw))cycle

            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+length+size_length
         enddo
         do iw=nw+1,nw+nwauxio

            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+length+size_length
         enddo
         write(qunit,'(a)')'</PointData>'

         write(qunit,'(a)')'<Points>'
         write(qunit,'(a,i16,a)') &
     '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
         ! write cell corner coordinates in a backward dimensional loop, always 3D output
         offset=offset+length_coords+size_length
         write(qunit,'(a)')'</Points>'
       case('pvtuBCCmpi')
         ! we write out every grid as one VTK PIECE
         write(qunit,'(a,i7,a,i7,a)') &
            '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
         write(qunit,'(a)')'<CellData>'
         do iw=1,nw
            if(.not.writew(iw))cycle

            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+lengthcc+size_length
         enddo
         do iw=nw+1,nw+nwauxio

            write(qunit,'(a,a,a,i16,a)')&
                '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
                '" format="appended" offset="',offset,'">'
            write(qunit,'(a)')'</DataArray>'
            offset=offset+lengthcc+size_length
         enddo
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
        call calc_grid(qunit,igrid,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                       ixC^L,ixCC^L,.true.)
        do iw=1,nw
          if(.not.writew(iw))cycle
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
!=============================================================================
