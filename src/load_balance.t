!=============================================================================
subroutine load_balance
use mod_forest
use mod_global_parameters

integer :: Morton_no, recv_igrid, recv_ipe, send_igrid, send_ipe, igrid, ipe

integer, external :: getnode

call get_Morton_range_active

if (npe==1) then
   sfc_to_igrid(:)=sfc(1,Morton_start(mype):Morton_stop(mype))
   return
end if

irecv=0
isend=0

do ipe=0,npe-1; do Morton_no=Morton_start(ipe),Morton_stop(ipe)
   recv_ipe=ipe

   send_igrid=sfc(1,Morton_no)
   send_ipe=sfc(2,Morton_no)

   if (recv_ipe/=send_ipe) then
      recv_igrid=getnode(recv_ipe)
      call change_ipe_tree_leaf(recv_igrid,recv_ipe,send_igrid,send_ipe)
      if (recv_ipe==mype) call lb_recv
      if (send_ipe==mype) call lb_send
   end if
   if (recv_ipe==mype) then
      if (recv_ipe==send_ipe) then
         sfc_to_igrid(Morton_no)=send_igrid
      else
         sfc_to_igrid(Morton_no)=recv_igrid
      end if
   end if
end do; end do

if (irecv>0) call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
if (isend>0) call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)

! post processing
do ipe=0,npe-1; do Morton_no=Morton_start(ipe),Morton_stop(ipe)
   recv_ipe=ipe

   send_igrid=sfc(1,Morton_no)
   send_ipe=sfc(2,Morton_no)

   if (recv_ipe/=send_ipe) then
      if (send_ipe==mype) call dealloc_node(send_igrid)
      call putnode(send_igrid,send_ipe)
   end if
end do; end do
{#IFDEF EVOLVINGBOUNDARY
! mark physical-boundary blocks on space-filling curve
do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   if (phyboundblock(igrid)) sfc_phybound(Morton_no)=1
end do
call MPI_ALLREDUCE(MPI_IN_PLACE,sfc_phybound,nleafs,MPI_INTEGER,&
                   MPI_SUM,icomm,ierrmpi)
}

! Update sfc array: igrid and ipe info in space filling curve
call amr_Morton_order()

!!$if (nwaux>0) then
!!$   do Morton_no=Morton_start(mype),Morton_stop(mype)
!!$      if (sfc(2,Morton_no)/=mype) then
!!$         igrid=sfc_to_igrid(Morton_no)
!!$         saveigrid=igrid
!!$         call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL,"load_balance")
!!$      end if
!!$   end do
!!$end if

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine lb_recv
!-----------------------------------------------------------------------------
call alloc_node(recv_igrid)

itag=recv_igrid
irecv=irecv+1
{#IFDEF EVOLVINGBOUNDARY
if (phyboundblock(recv_igrid)) then
   call MPI_IRECV(pw(recv_igrid)%w,1,type_block,send_ipe,itag, &
                  icomm,recvrequest(irecv),ierrmpi)
else
   call MPI_IRECV(pw(recv_igrid)%w,1,type_block_io,send_ipe,itag, &
                  icomm,recvrequest(irecv),ierrmpi)
end if
}{#IFNDEF EVOLVINGBOUNDARY
call MPI_IRECV(pw(recv_igrid)%w,1,type_block_io,send_ipe,itag, &
               icomm,recvrequest(irecv),ierrmpi)
}

end subroutine lb_recv
!=============================================================================
subroutine lb_send
!-----------------------------------------------------------------------------
itag=recv_igrid
isend=isend+1
{#IFDEF EVOLVINGBOUNDARY
if (phyboundblock(send_igrid)) then
   call MPI_ISEND(pw(send_igrid)%w,1,type_block,recv_ipe,itag, &
                  icomm,sendrequest(isend),ierrmpi)
else
   call MPI_ISEND(pw(send_igrid)%w,1,type_block_io,recv_ipe,itag, &
                  icomm,sendrequest(isend),ierrmpi)
end if
}{#IFNDEF EVOLVINGBOUNDARY
call MPI_ISEND(pw(send_igrid)%w,1,type_block_io,recv_ipe,itag, &
               icomm,sendrequest(isend),ierrmpi)
}

end subroutine lb_send
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine load_balance
!=============================================================================
subroutine level1_Morton_order
! use Morton curve to connect level-1 grid blocks
use mod_forest
use mod_global_parameters

integer, allocatable :: gsq_sfc(:^D&)
integer :: ig^D, ngsq^D, Morton_no
integer(kind=8), external :: mortonEncode
!-----------------------------------------------------------------------------
! use the smallest square/cube to cover the full domain 
ngsq^D=2**ceiling(log(real(ng^D(1)))/log(2.0));
{^NOONED
{ngsq^D=max(ngsq^DD) \}
}
allocate(gsq_sfc(ngsq^D))
! get Morton-order numbers in the square/cube
{do ig^DB=1,ngsq^DB\}
   gsq_sfc(ig^D)=int(mortonEncode(ig^D-1,ndim))+1
{end do\}
! delete Morton blocks that are out of the domain
{do ig^DB=1,ngsq^DB\}
   if (ig^D>ng^D(1)|.or.) then
      where(gsq_sfc>=gsq_sfc(ig^D))
         gsq_sfc=gsq_sfc-1
      end where 
   end if  
{end do\}
! copy the modified Morton numbers to the blocks in the domain
allocate(iglevel1_sfc(ng^D(1)))
allocate(sfc_iglevel1(ndim,nglev1))
{do ig^DB=1,ng^DB(1)\}
   iglevel1_sfc(ig^D)=gsq_sfc(ig^D)
   {sfc_iglevel1(^D,iglevel1_sfc(ig^DD))=ig^D \}
{end do\}
!do Morton_no=1,nglev1
!   ig^D=sfc_iglevel1(^D,Morton_no)\ 
!   print*,'Morton',Morton_no,'ig',ig^D
!end do
!stop

end subroutine level1_Morton_order
!=============================================================================
integer(kind=8) function mortonEncode(ig^D,ndim)
use iso_fortran_env, only : int64
implicit none
integer(kind=4) :: i,ig^D,ndim
integer(kind=8) :: answer
answer = 0
do i=0,64/ndim
  {^IFONED answer=ig1}
  {^IFTWOD
   answer=ior(answer,ior(ishft(iand(ig1,ishft(1_int64,i)),i),&
          ishft(iand(ig2,ishft(1_int64,i)),i+1)))
  \}
  {^IFTHREED
   answer=ior(answer,ior(ishft(iand(ig1,ishft(1_int64,i)),2*i),ior(ishft(&
   iand(ig2,ishft(1_int64,i)),2*i+1),ishft(iand(ig3,ishft(1_int64,i)),2*i+2))))
  \}
end do
mortonEncode=answer
return
end function mortonEncode
!=============================================================================
subroutine amr_Morton_order
! Construct Morton-order as a global recursive lexicographic ordering.

use mod_forest
use mod_global_parameters

integer :: ig^D, Morton_no, isfc
!-----------------------------------------------------------------------------
Morton_no=0
nglev1={ng^D(1)*}
do isfc=1,nglev1
   ig^D=sfc_iglevel1(^D,isfc)\ 
   call get_Morton_number(tree_root(ig^D))
end do

if (Morton_no/=nleafs) then
   call mpistop("error in amr_Morton_order: Morton_no/=nleafs")
end if

contains
!=============================================================================
! internal procedures
!=============================================================================
recursive subroutine get_Morton_number(tree)

type(tree_node_ptr) :: tree

integer :: ic^D
!-----------------------------------------------------------------------------
if (tree%node%leaf) then
   Morton_no=Morton_no+1
   sfc(1,Morton_no)=tree%node%igrid
   sfc(2,Morton_no)=tree%node%ipe
   if (tree%node%active) then 
      sfc(3,Morton_no)=1 
   else 
      sfc(3,Morton_no)=0 
   end if
   igrid_to_sfc(tree%node%igrid)=Morton_no
else
   {do ic^DB=1,2\}
      call get_Morton_number(tree%node%child(ic^D))
   {end do\}
end if

end subroutine get_Morton_number
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine amr_Morton_order
!=============================================================================
subroutine get_Morton_range
use mod_forest
use mod_global_parameters

integer :: ipe
!-----------------------------------------------------------------------------
if (allocated(sfc_to_igrid)) deallocate(sfc_to_igrid)
{#IFDEF EVOLVINGBOUNDARY
if (allocated(sfc_phybound)) deallocate(sfc_phybound)
}
do ipe=0,npe-1
   Morton_start(ipe)=1+ipe*int(nleafs/npe)+min(ipe,mod(nleafs,npe))
   Morton_stop(ipe)=(ipe+1)*int(nleafs/npe)+min(ipe+1,mod(nleafs,npe))
end do
if (Morton_stop(mype)>=Morton_start(mype)) then
   allocate(sfc_to_igrid(Morton_start(mype):Morton_stop(mype)))
{#IFDEF EVOLVINGBOUNDARY
   allocate(sfc_phybound(nleafs))
   sfc_phybound=0
}
end if

end subroutine get_Morton_range
!=============================================================================

!=============================================================================
subroutine get_Morton_range_active
use mod_forest
use mod_global_parameters

! Cut the sfc based on weighted decision.  
! Oliver Porth, 02.02.2012

!!!Here you choose the weithts:!!
integer, parameter :: wa=3, wp=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! wp : Weight for passive block
! wa : Weight for active block
! wp=0 : balance load (active blocks) exactly and 
! don't care about memory imbalance
! wp=wa : balance memory exactly and don't care about load
! Maximum possible memory imbalance is X=wa/wp.

! If you run into memory issues, decrease this ratio.
! Scaling should be better if you allow for higher ratio.  
! Note also that passive cells still do regridding and boundary swap.  
! I have best results with a ratio 2:1, but it is problem dependent.
! It can't do magic though...
! Best to make sure that the sfc is properly aligned with your problem. 

integer :: ipe, Morton_no
integer :: Mtot, Mstop, Mcurr
! For debugging: 
integer :: nactive(0:npe-1),npassive(0:npe-1)
!double precision, save :: ptasum=0
!-----------------------------------------------------------------------------
if (allocated(sfc_to_igrid)) deallocate(sfc_to_igrid)
{#IFDEF EVOLVINGBOUNDARY
if (allocated(sfc_phybound)) deallocate(sfc_phybound)
}

Mtot  = nleafs_active*wa+(nleafs-nleafs_active)*wp
ipe = 0 
Mcurr = 0

nactive=0
npassive=0

Morton_start(0) = 1
do Morton_no=1,nleafs
   ! This is where we ideally would like to make the cuts:
   Mstop  = (ipe+1)*int(Mtot/npe)+min(ipe+1,mod(Mtot,npe))
   ! Build up mass:
   Mcurr = Mcurr + (wa*sfc(3,Morton_no)+wp*(1-sfc(3,Morton_no)))

   if (sfc(3,Morton_no)==1) then 
      nactive(ipe) = nactive(ipe) +1
      else
         npassive(ipe) = npassive(ipe) +1
      end if

   if (Mcurr >= Mstop) then 
      Morton_stop(ipe) = Morton_no
      ipe = ipe +1
      if (ipe>=npe) exit
      Morton_start(ipe) = Morton_no + 1
   end if
end do

Xmemory=dble(maxval(npassive+nactive))/&
     dble(minval(npassive+nactive))
Xload=dble(maxval(nactive))/&
     dble(minval(nactive))

!ptasum = ptasum +dble(nleafs-nleafs_active)/dble(nleafs_active)

!if (mype == 0) print*, 'nleafs_passive:',nleafs-nleafs_active, 'nleafs_active:',nleafs_active,'ratio:',dble(nleafs-nleafs_active)/dble(nleafs_active),'mean ratio:',ptasum/it



if (Morton_stop(mype)>=Morton_start(mype)) then
   allocate(sfc_to_igrid(Morton_start(mype):Morton_stop(mype)))
{#IFDEF EVOLVINGBOUNDARY
   allocate(sfc_phybound(nleafs))
   sfc_phybound=0
}
end if

end subroutine get_Morton_range_active
!=============================================================================
