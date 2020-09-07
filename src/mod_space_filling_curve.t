module mod_space_filling_curve
  implicit none

contains

  !> build Morton space filling curve for level 1 grid
  subroutine level1_Morton_order
    ! use Morton curve to connect level 1 grid blocks
    use mod_forest
    use mod_global_parameters

    integer, allocatable :: gsq_sfc(:^D&),seq_sfc(:),seq_ig^D(:)
    integer :: ig^D, ngsq^D, isq, total_number
    !integer(kind=8), external :: mortonEncode
    logical, allocatable :: in_domain(:)

    ! use the smallest square/cube to cover the full domain 
    ngsq^D=2**ceiling(log(real(ng^D(1)))/log(2.0));
    {^NOONED
    {ngsq^D=max(ngsq^DD) \}
    }
    total_number={ngsq^D|*}
    ! Morton number acquired by block numbers
    allocate(gsq_sfc(ngsq^D))
    ! Morton number in sequence
    allocate(seq_sfc(total_number))
    ! block numbers in sequence
    {allocate(seq_ig^D(total_number))\}
    allocate(in_domain(total_number))
    in_domain=.true.
    ! get Morton-order numbers in the square/cube
    {do ig^DB=1,ngsq^DB\}
       gsq_sfc(ig^D)=int(mortonEncode(ig^D-1,ndim))+1
       seq_sfc(gsq_sfc(ig^D))=gsq_sfc(ig^D)
       {seq_ig^D(gsq_sfc(ig^DD))=ig^D \}
    {end do\}
    ! mark blocks that are out of the domain and change Morton number
    do isq=1,total_number
       if (seq_ig^D(isq)>ng^D(1)|.or.) then
         seq_sfc(isq:total_number)=seq_sfc(isq:total_number)-1
         in_domain(isq)=.false.
       end if
    end do
    ! copy the modified Morton numbers to the blocks in the domain
    allocate(iglevel1_sfc(ng^D(1)))
    allocate(sfc_iglevel1(ndim,nglev1))
    do isq=1,total_number
      if(in_domain(isq)) then
        iglevel1_sfc(seq_ig^D(isq))=seq_sfc(isq)
        {sfc_iglevel1(^D,seq_sfc(isq))=seq_ig^D(isq) \}
      end if
    end do

    deallocate(gsq_sfc,seq_sfc,seq_ig^D,in_domain)

  end subroutine level1_Morton_order

  integer(kind=8) function mortonEncode(ig^D,ndim)
    use iso_fortran_env, only : int64
    implicit none
    integer(kind=4), intent(in) :: ig^D,ndim
    integer(kind=4) :: i
    integer(kind=8) :: answer, lg^D

    ! Create a 64-bit version of ig^D
    lg^D=ig^D;
    answer = 0

    do i=0,64/ndim
      {^IFONED answer=ig1}
      {^IFTWOD
       answer=ior(answer,ior(ishft(iand(lg1,ishft(1_int64,i)),i),&
              ishft(iand(lg2,ishft(1_int64,i)),i+1)))
      \}
      {^IFTHREED
       answer=ior(answer,ior(ishft(iand(lg1,ishft(1_int64,i)),2*i),ior(ishft(&
       iand(lg2,ishft(1_int64,i)),2*i+1),ishft(iand(lg3,ishft(1_int64,i)),2*i+2))))
      \}
    end do
    mortonEncode=answer
    return
  end function mortonEncode

  !> Construct Morton-order as a global recursive lexicographic ordering.
  subroutine amr_Morton_order
    use mod_forest
    use mod_global_parameters

    integer :: ig^D, Morton_no, isfc

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

    recursive subroutine get_Morton_number(tree)

      type(tree_node_ptr) :: tree

      integer :: ic^D

      if (tree%node%leaf) then
         Morton_no=Morton_no+1
         sfc(1,Morton_no)=tree%node%igrid
         sfc(2,Morton_no)=tree%node%ipe
         if (tree%node%active) then 
            sfc(3,Morton_no)=1 
         else 
            sfc(3,Morton_no)=0 
         end if
         if(tree%node%ipe==mype) igrid_to_sfc(tree%node%igrid)=Morton_no
      else
         {do ic^DB=1,2\}
            call get_Morton_number(tree%node%child(ic^D))
         {end do\}
      end if

    end subroutine get_Morton_number

  end subroutine amr_Morton_order

  !> Set the Morton range for each processor
  subroutine get_Morton_range
    use mod_forest
    use mod_global_parameters

    integer :: ipe, blocks_left, procs_left, num_blocks

    if (allocated(sfc_to_igrid)) deallocate(sfc_to_igrid)
    {#IFDEF EVOLVINGBOUNDARY
    if (allocated(sfc_phybound)) deallocate(sfc_phybound)
    }

    blocks_left = nleafs

    do ipe = 0, npe-1
      if (ipe == 0) then
        Morton_start(ipe) = 1
      else
        Morton_start(ipe) = Morton_stop(ipe-1) + 1
      end if

      ! Compute how many blocks this cpu should take
      procs_left = npe - ipe
      num_blocks = ceiling(blocks_left / dble(procs_left))
      Morton_stop(ipe) = Morton_start(ipe) + num_blocks - 1
      blocks_left = blocks_left - num_blocks
    end do

    allocate(sfc_to_igrid(Morton_start(mype):Morton_stop(mype)))
    {#IFDEF EVOLVINGBOUNDARY
    allocate(sfc_phybound(nleafs))
    sfc_phybound=0
    }

  end subroutine get_Morton_range

  subroutine get_Morton_range_active
    use mod_forest
    use mod_global_parameters

    ! Cut the sfc based on weighted decision.  
    ! Oliver Porth, 02.02.2012

    !!!Here you choose the weithts:!!
    integer, parameter :: wa=3, wp=1
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

end module mod_space_filling_curve
