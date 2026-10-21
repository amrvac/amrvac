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
    if(.not. allocated(iglevel1_sfc)) allocate(iglevel1_sfc(ng^D(1)))
    if(.not. allocated(sfc_iglevel1)) allocate(sfc_iglevel1(ndim,nglev1))
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
    integer(kind=8) :: answer, lg^D
    integer(kind=4), intent(in) :: ig^D,ndim
    integer(kind=4) :: i

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
    use mod_comm_lib, only: mpistop

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

  !> Cost-weighted SFC partition. Reduces the per-rank block_cost array
  !> into a global Morton-indexed cost, EWMA-blends with the persistent
  !> costlist, and partitions the Morton order so each rank's cumulative
  !> cost is balanced. Falls back to equal-block-count if all costs are
  !> zero (cold start before the first measured step).
  subroutine get_Morton_range_costed
    use mod_forest
    use mod_global_parameters

    integer :: ipe, Morton_no, igrid, ix, nseen, ncut, nmax, nrem
    integer :: ncut_lo, ncut_hi
    integer, save :: last_report = -1
    double precision :: cost_total, cost_target, cost_mean
    double precision :: maxcount, mincount, maxcost, meancost
    double precision, allocatable :: cost_local(:), cost_cumul(:)
    integer :: nblocks_per(0:npe-1)
    logical :: partition_ok
    logical, save :: warned_fallback = .false.

    if (allocated(sfc_to_igrid)) deallocate(sfc_to_igrid)
    {#IFDEF EVOLVINGBOUNDARY
    if (allocated(sfc_phybound)) deallocate(sfc_phybound)
    }

    ! 1. Each rank fills its measured per-step block cost into the
    !    Morton-indexed cost_local at its own slots, zeros elsewhere.
    !    Morton numbering is invariant across load_balance migration, so
    !    costlist values built up via EWMA below stay meaningful when
    !    blocks change rank. Refinement events insert and remove Morton indices;
    !    slots created that way are seeded at step 3b and then converge on their
    !    own measurements.
    allocate(cost_local(nleafs))
    cost_local = 0.0d0
    do Morton_no = 1, nleafs
      if (sfc(2,Morton_no) == mype) then
        igrid = sfc(1,Morton_no)
        if (igrid > 0 .and. igrid <= max_blocks) then
          cost_local(Morton_no) = block_cost(igrid)
        end if
      end if
    end do

    ! 2. Reduce so every rank has the global per-Morton cost vector. Since
    !    each Morton index is owned by exactly one rank in the pre-balance
    !    state, the SUM is just the value contributed by that one rank.
    if (npe > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, cost_local, nleafs, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierrmpi)
    end if

    ! 3. EWMA blend the measurement into the persistent global costlist.
    !    Skip slots with zero measurement (block did not run advect this
    !    cycle for some reason) and keep the prior estimate. This blend at
    !    the GLOBAL Morton level (not per-igrid) is what makes the
    !    partition migration-safe.
    do Morton_no = 1, nleafs
      if (cost_local(Morton_no) > 0.0d0) then
        if (costlist_seen(Morton_no)) then
          costlist(Morton_no) = lb_alpha * costlist(Morton_no) &
                              + (1.0d0 - lb_alpha) * cost_local(Morton_no)
        else
          ! First measurement for this slot: take it rather than blend it
          ! against a value that was never measured.
          costlist(Morton_no) = cost_local(Morton_no)
          costlist_seen(Morton_no) = .true.
        end if
      end if
    end do
    deallocate(cost_local)

    ! 3b. Give never-measured slots the mean of the measured ones. costlist is
    !     a wall time, so any fixed seed is a units error and would let a single
    !     untimed leaf dominate cost_total.
    nseen = count(costlist_seen(1:nleafs))
    if (nseen > 0) then
      cost_mean = sum(costlist(1:nleafs), mask=costlist_seen(1:nleafs)) / dble(nseen)
      do Morton_no = 1, nleafs
        if (.not. costlist_seen(Morton_no)) costlist(Morton_no) = cost_mean
      end do
    end if

    cost_total = sum(costlist(1:nleafs))
    if (cost_total <= 0.0d0 .or. nseen == 0 .or. nleafs < npe) then
      ! Fully cold, or fewer leaves than ranks (no cost-based cut can give
      ! every rank work): fall back to the equal-block partition.
      call get_Morton_range
      return
    end if

    ! 4. Cumulative-cost cut along the Morton-sorted leaves. Each rank is closed
    !    explicitly, in rank order, so the ranges are always a contiguous cover
    !    of 1..nleafs. Morton_start/Morton_stop persist between calls, so a rank
    !    left unassigned here would silently retain a stale range.
    allocate(cost_cumul(0:nleafs))
    cost_cumul(0) = 0.0d0
    do Morton_no = 1, nleafs
      cost_cumul(Morton_no) = cost_cumul(Morton_no-1) + costlist(Morton_no)
    end do

    ! Upper bound on the blocks one rank may hold, from lb_max_block_ratio.
    nmax = max(1, ceiling(lb_max_block_ratio * dble(nleafs) / dble(npe)))

    nblocks_per = 0
    ix = 0                                  ! last leaf handed out so far
    Morton_start(0) = 1
    do ipe = 0, npe-2
      nrem = npe-1-ipe                      ! ranks still to be assigned
      cost_target = (dble(ipe) + 1.0d0) * cost_total / dble(npe)
      ! smallest leaf whose cumulative cost reaches this rank's target
      ncut = ix + 1
      do while (ncut < nleafs .and. cost_cumul(ncut) < cost_target)
        ncut = ncut + 1
      end do
      ! Feasible window: lo keeps at least one leaf here and leaves the
      ! remaining ranks no more than nmax each; hi respects nmax for this rank
      ! and leaves at least one leaf per remaining rank.
      ncut_lo = max(ix + 1,     nleafs - nmax*nrem)
      ncut_hi = min(ix + nmax,  nleafs - nrem)
      ncut = min(max(ncut, ncut_lo), ncut_hi)
      Morton_stop(ipe)    = ncut
      Morton_start(ipe+1) = ncut + 1
      nblocks_per(ipe)    = ncut - ix
      ix = ncut
    end do
    Morton_stop(npe-1) = nleafs
    nblocks_per(npe-1) = nleafs - ix
    deallocate(cost_cumul)

    ! 4b. Verify the cover before returning: load_balance migrates blocks from
    !     these ranges and a malformed partition corrupts the forest rather than
    !     failing cleanly. Fall back to the equal-block cut if it is wrong.
    partition_ok = (Morton_start(0) == 1 .and. Morton_stop(npe-1) == nleafs)
    do ipe = 1, npe-1
      if (Morton_start(ipe) /= Morton_stop(ipe-1) + 1) partition_ok = .false.
    end do
    do ipe = 0, npe-1
      if (Morton_stop(ipe) < Morton_start(ipe)) partition_ok = .false.
      if (Morton_start(ipe) < 1 .or. Morton_stop(ipe) > nleafs) partition_ok = .false.
    end do
    if (.not. partition_ok) then
      if (mype == 0 .and. .not. warned_fallback) then
        write(*,*) 'get_Morton_range_costed: partition failed validation at it=', it
        write(*,*) '  nleafs=', nleafs, ' npe=', npe, ' cost_total=', cost_total
        write(*,*) '  falling back to equal-block partition (warned once)'
        warned_fallback = .true.
      end if
      call get_Morton_range
      return
    end if

    ! 5. Diagnostics (also drive the existing log columns Xload/Xmemory).
    maxcount = dble(maxval(nblocks_per))
    mincount = dble(max(1,minval(nblocks_per)))
    Xmemory = maxcount / mincount
    maxcost  = 0.0d0
    do ix = 0, npe-1
      maxcost = max(maxcost, sum(costlist(Morton_start(ix):Morton_stop(ix))))
    end do
    meancost = cost_total / dble(npe)
    Xload    = maxcost / max(meancost, tiny(1.0d0))

    ! Report the achieved block spread and load imbalance, rank 0, at most once
    ! per 500 iterations.
    if (mype == 0 .and. it >= last_report + 500) then
      write(*,'(a,i9,a,i6,a,i5,a,i5,a,i5,a,f6.2)') &
        ' [lb] it=', it, ' nleafs=', nleafs, &
        '  blocks/rank min=', minval(nblocks_per), ' max=', maxval(nblocks_per), &
        ' cap=', nmax, '  Xload=', Xload
      last_report = it
    end if

    if (Morton_stop(mype) >= Morton_start(mype)) then
      allocate(sfc_to_igrid(Morton_start(mype):Morton_stop(mype)))
      {#IFDEF EVOLVINGBOUNDARY
      allocate(sfc_phybound(nleafs))
      sfc_phybound = 0
      }
    end if

  end subroutine get_Morton_range_costed

end module mod_space_filling_curve
