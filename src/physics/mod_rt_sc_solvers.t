!=============================================================================
!> Solvers of the short-characteristics (SC) radiative transfer and their dispatch
!>
!> rt_sc_solve runs once per step. It sweeps every ray through every band on rt_sc_amr_level and
!> leaves the heating Q in i_qrad and the vertical flux Fy in i_frad. The pipelined KBA wavefront
!> is the default. Its sweep corners run at once, each in its own intensity slot, and the rays of
!> a corner follow one another through the wavefront. The ray-by-ray wavefront passes one ray at
!> a time, agrees with it to round-off, and is the simplest to read. In both, each rank
!> sweeps only its own band blocks, one diagonal of blocks at a time, with a getbc of the ray
!> intensity between diagonals, so that no rank holds the whole band. The replicated gather
!> assembles chi and B of the whole band on every rank and is kept as the reference. It differs
!> from the wavefront only where rays cross a periodic seam
!>
!> The loops over block cells are written along (hdim, updim) through hv, so that one loop
!> serves either stratification axis
!>
!> Structure:
!>   A. dispatch      rt_sc_solve and the validation rt_sc_validate_pipe;
!>   B. gather        rt_sc_solve_gather and rt_grid_index;
!>   C. band blocks   rt_band_blocks, which lists the band blocks and caches chi and B, with
!>                    rt_tau_tops and rt_block_tau for tau from the band top, which
!>                    rt_sc_emergent_tau1 shares;
!>   D. wave stages   rt_wave_begin, rt_wave_band_open, rt_wave_band_close, which forms Q, and
!>                    rt_wave_end, with rt_add_moments and rt_exchange;
!>   E. block work    rt_seam_read and rt_seam_write for the periodic seam, rt_inflow_ghosts for
!>                    the band-end and edge inflows, and rt_sweep_one;
!>   F. schedules     rt_sc_solve_wave and rt_sc_solve_wave_pipe, with wave_stage
!>
!> Jack Jenkins (10.2026); More or less finished in v1 state for amrvac4.0
!=============================================================================
module mod_rt_sc_solvers
  use mod_rt_sc_state
  implicit none

  public :: rt_sc_solve

contains

  !> Collective, once per step. For each band on rt_sc_amr_level it sweeps all rays and forms the
  !> heating Q, the thin/thick blend with the tau^2 taper, in i_qrad, which is zero elsewhere,
  !> together with Fy, and Iout on the gather path. rt_sc_add_source applies the energy every
  !> stage. The gather path is for comparison only, see rt_sc_solve_gather
  subroutine rt_sc_solve()
    use mod_global_parameters
    use mod_timing, only: time_rt, timer_open, timer_close
    double precision :: tw(3)
{^IFTWOD
    if(.not.rt_sc_active) return
    if(rt_sc_valpipe .and. rt_sc_distributed) then
      call rt_sc_validate_pipe()
      return
    end if
    call timer_open(tw)
    if(rt_sc_distributed) then
      if(rt_sc_pipeline) then
        call rt_sc_solve_wave_pipe()
      else
        call rt_sc_solve_wave()
      end if
      rt_wt_wave = rt_wt_wave + (MPI_WTIME()-tw(1))
      rt_nsolve = rt_nsolve + 1
      if(rt_sc_debug .and. mype==0 .and. mod(rt_nsolve,10)==0) write(*,'(a,i5,a,f9.3,a)') &
        ' RT wave nsolve=',rt_nsolve,'  cumulative wave-solve time=',rt_wt_wave,' s'
    else
      call rt_sc_solve_gather()
    end if
    time_rt = time_rt + timer_close(tw)
    ! The per-rank wall time, which includes the pipeline waits that show the imbalance
    if(lb_diagnose) lb_rt_accum = MPI_WTIME() - tw(1)
}
  end subroutine rt_sc_solve

  !> Replicated gather, the original implementation, kept for comparison and not for production.
  !> Every rank holds chi and B of the whole band, the rays are shared out over the node-local
  !> ranks, and the moments are reduced within the node. Its memory grows with the whole band on
  !> every rank, and its periodic seam takes a vacuum inflow for rays swept column by column
  !> (sc_sweep_2d). Use KBA, pipelined or not
  subroutine rt_sc_solve_gather()
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    {^IFTWOD
    double precision :: chi_l(ixG^T), b_l(ixG^T), tw
    double precision, allocatable :: fill_g(:,:)
    integer :: iigrid, igrid, ix1, ix2, ig, jg, k, ierr, ib
    if(.not.rt_sc_active) return
    if(.not.rt_sc_ready) call mpistop("rt_sc_solve_gather: rt_sc_setup not called")
    if(i_qrad < 0)       call mpistop("rt_sc_solve_gather: i_qrad not registered")
    if(.not.associated(rt_sc_chiB)) call mpistop("rt_sc_solve_gather: opacity callback not set")

    ! The node-local and leader communicators are built on the first call. Every rank joins its
    ! node communicator, and only the node leaders join the leader communicator
    if(rt_node_comm == -1) then
      call MPI_Comm_split_type(icomm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, rt_node_comm, ierr)
      call MPI_Comm_rank(rt_node_comm, rt_node_rank, ierr)
      call MPI_Comm_size(rt_node_comm, rt_node_size, ierr)
      if(rt_node_rank == 0) then
        call MPI_Comm_split(icomm, 0, mype, rt_lead_comm, ierr)
      else
        call MPI_Comm_split(icomm, MPI_UNDEFINED, mype, rt_lead_comm, ierr)
      end if
    end if

    nxg = rg_nh; dxg = rg_dh; dyg = rg_dv

    ! Q is zeroed everywhere first, so the cells outside the bands receive no heating
    do iigrid=1,igridstail; igrid=igrids(iigrid)
      ps(igrid)%w(ixG^T,i_qrad) = 0.0d0
    end do

    do ib = 1, rt_sc_nband
      ! Each band has its own incoming conditions and its own number of rows
      bc_lo_thermal = (trim(rt_sc_bc_min(ib)) == 'thermal')
      bc_hi_thermal = (trim(rt_sc_bc_max(ib)) == 'thermal')
      nyg = band_nv(ib)
      call rt_grid_alloc()
      ! fill_g counts how many cells land on each grid point, which checks the AMR pin below
      allocate(fill_g(nxg,nyg)); fill_g = 0.0d0
      chi_g = 0.0d0; b_g = 0.0d0
      tw = MPI_WTIME()

      ! Each rank writes chi and B of its own level-L band cells into the global grid
      do iigrid=1,igridstail; igrid=igrids(iigrid)
        if(node(plevel_,igrid) /= rt_sc_amr_level) cycle
        call rt_sc_chiB(ixG^LL, ixM^LL, ps(igrid)%w, ps(igrid)%x, chi_l, b_l)
        do ix2=ixMlo2,ixMhi2; do ix1=ixMlo1,ixMhi1
          if(.not.rt_grid_index(ps(igrid)%x(ix1,ix2,hdim), ps(igrid)%x(ix1,ix2,updim), ib, ig, jg)) cycle
          chi_g(ig,jg) = chi_l(ix1,ix2); b_g(ig,jg) = b_l(ix1,ix2); fill_g(ig,jg) = 1.0d0
        end do; end do
      end do
      ! Each point is written by one rank and is zero on all the others, so summing within the
      ! node, then over the leaders, and broadcasting within the node gives every rank the field
      call rt_allreduce_chunk(chi_g,  nxg, nyg, rt_node_comm)
      call rt_allreduce_chunk(b_g,    nxg, nyg, rt_node_comm)
      call rt_allreduce_chunk(fill_g, nxg, nyg, rt_node_comm)
      if(rt_lead_comm /= MPI_COMM_NULL) then
        call rt_allreduce_chunk(chi_g,  nxg, nyg, rt_lead_comm)
        call rt_allreduce_chunk(b_g,    nxg, nyg, rt_lead_comm)
        call rt_allreduce_chunk(fill_g, nxg, nyg, rt_lead_comm)
      end if
      if(rt_node_size > 1) then
        call MPI_Bcast(chi_g,  nxg*nyg, MPI_DOUBLE_PRECISION, 0, rt_node_comm, ierr)
        call MPI_Bcast(b_g,    nxg*nyg, MPI_DOUBLE_PRECISION, 0, rt_node_comm, ierr)
        call MPI_Bcast(fill_g, nxg*nyg, MPI_DOUBLE_PRECISION, 0, rt_node_comm, ierr)
      end if
      ! Every grid point must have been filled exactly once. Anything else means the band is not
      ! one complete level and the AMR pin has failed
      if(minval(fill_g) /= 1.0d0 .or. maxval(fill_g) /= 1.0d0) then
        if(mype==0) print *,'rt_sc_solve_gather: band',ib,' not a complete single level ', &
             rt_sc_amr_level,' (minfill',minval(fill_g),' maxfill',maxval(fill_g),')'
        call mpistop("rt_sc_solve_gather: RT band is not a complete single AMR level")
      end if
      deallocate(fill_g)
      call sc_pad_grid()
      rt_wt_gath = rt_wt_gath + (MPI_WTIME()-tw); tw = MPI_WTIME()

      ! Without scattering each ray is one independent sweep, so the rays are shared out over the
      ! node-local ranks and the moments summed within the node. Each sweep starts from zero
      ! intensity, so the sum depends on the order and the number of ranks only through round-off
      jm_g = 0.0d0; fx_g = 0.0d0; fy_g = 0.0d0
      do k = rt_node_rank+1, rt_sc_nray, rt_node_size
        call sc_sweep_2d(mu_h(k), mu_v(k))
        ! J = sum_k w_k I_k,  Fx = sum_k w_k mu_h,k I_k,  Fy = sum_k w_k mu_v,k I_k
        jm_g = jm_g + wmu(k)*ii_g
        fx_g = fx_g + wmu(k)*mu_h(k)*ii_g
        fy_g = fy_g + wmu(k)*mu_v(k)*ii_g
      end do
      rt_wt_swp = rt_wt_swp + (MPI_WTIME()-tw); tw = MPI_WTIME()
      if(rt_node_size>1) then
        call rt_allreduce_chunk(jm_g, nxg, nyg, rt_node_comm)
        call rt_allreduce_chunk(fx_g, nxg, nyg, rt_node_comm)
        call rt_allreduce_chunk(fy_g, nxg, nyg, rt_node_comm)
      end if
      rt_wt_red = rt_wt_red + (MPI_WTIME()-tw); tw = MPI_WTIME()
      ! F = c sum_k w_k mu_k I_k
      fx_g = c_norm*fx_g; fy_g = c_norm*fy_g
      if(rt_sc_flux_blend) then
        call rt_sc_heat_blend()
      else
        ! Q = c chi (J - B)
        q_g = c_norm*chi_g*(jm_g - b_g)
      end if
      ! A vertical ray swept afterwards gives the emergent intensity Iout
      call sc_sweep_2d(0.0d0, 1.0d0)

      ! Q, Fy and Iout go back to the band's level-L cells through the same index map
      do iigrid=1,igridstail; igrid=igrids(iigrid)
        if(node(plevel_,igrid) /= rt_sc_amr_level) cycle
        do ix2=ixMlo2,ixMhi2; do ix1=ixMlo1,ixMhi1
          if(.not.rt_grid_index(ps(igrid)%x(ix1,ix2,hdim), ps(igrid)%x(ix1,ix2,updim), ib, ig, jg)) cycle
          ps(igrid)%w(ix1,ix2,i_qrad) = q_g(ig,jg)
          if(i_frad>0) ps(igrid)%w(ix1,ix2,i_frad) = fy_g(ig,jg)
          if(i_iout>0) ps(igrid)%w(ix1,ix2,i_iout) = ii_g(ig,jg)
        end do; end do
      end do
      rt_wt_scat = rt_wt_scat + (MPI_WTIME()-tw)
    end do
    rt_nsolve = rt_nsolve + 1
    if(rt_sc_debug .and. mype==0 .and. mod(rt_nsolve,10)==0) &
      write(*,'(a,i5,4(a,f8.3),a,f8.3,a)') ' RT timing nsolve=',rt_nsolve, &
        ' | gather+fieldreduce=',rt_wt_gath,' sweep=',rt_wt_swp,' momentreduce=',rt_wt_red, &
        ' blend+scatter=',rt_wt_scat,' | total=',rt_wt_gath+rt_wt_swp+rt_wt_red+rt_wt_scat,' s (cumulative)'
    }
  end subroutine rt_sc_solve_gather

  !> Indices (ig, jg) on the replicated grid of band ib of the cell centred at xh along hdim and xv
  !> along updim, false when the cell lies outside the band
  logical function rt_grid_index(xh, xv, ib, ig, jg) result(inside)
    double precision, intent(in) :: xh, xv
    integer, intent(in)          :: ib
    integer, intent(out)         :: ig, jg
    inside = .false.
    if(xv <= band_lo_c(ib) .or. xv >= band_hi_c(ib)) return
    ! ig = nint((x - x_lo)/dx + 1/2),  jg = nint((y - y_band)/dy + 1/2)
    ig = nint((xh - rg_hlo)/rg_dh + 0.5d0)
    jg = nint((xv - band_lo_c(ib))/rg_dv + 0.5d0)
    inside = ig >= 1 .and. ig <= rg_nh .and. jg >= 1 .and. jg <= band_nv(ib)
  end function rt_grid_index

  {^IFTWOD
  !> Validation of the pipelined KBA. The ray-by-ray Q is copied into i_iout and the pipelined Q
  !> left in i_qrad, and the largest difference is printed relative to the largest Q, both over
  !> the whole band and away from the seam blocks, the outermost x-block on either side. Both
  !> solves start from the same seam inflow, since the first one hands its own on to the second
  subroutine rt_sc_validate_pipe()
    use mod_global_parameters
    integer :: iigrid, igrid, ierr, ix1, ix2
    double precision :: dd, dmx, qmx, dmi, xedge, gbuf(3), gout(3)
    double precision, allocatable :: seamL0(:,:), seamR0(:,:)
    logical :: had_seam
    had_seam = allocated(rt_sc_seamL_prev)
    if(had_seam) then
      seamL0 = rt_sc_seamL_prev; seamR0 = rt_sc_seamR_prev
    end if
    call rt_sc_solve_wave()
    if(allocated(rt_sc_seamL_prev)) then
      if(had_seam) then
        rt_sc_seamL_prev = seamL0; rt_sc_seamR_prev = seamR0
      else
        rt_sc_seamL_prev = 0.0d0; rt_sc_seamR_prev = 0.0d0
      end if
    end if
    do iigrid=1,igridstail; igrid=igrids(iigrid)
      ps(igrid)%w(ixG^T,i_iout) = ps(igrid)%w(ixG^T,i_qrad)
    end do
    call rt_sc_solve_wave_pipe()
    dmx=0.0d0; qmx=0.0d0; dmi=0.0d0
    xedge = block_nx1*(xprobmax1-xprobmin1)/dble(domain_nx1)
    do iigrid=1,igridstail; igrid=igrids(iigrid)
      qmx=max(qmx, maxval(dabs(ps(igrid)%w(ixM^T,i_qrad))))
      do ix2=ixMlo2,ixMhi2; do ix1=ixMlo1,ixMhi1
        dd=dabs(ps(igrid)%w(ix1,ix2,i_qrad)-ps(igrid)%w(ix1,ix2,i_iout))
        dmx=max(dmx,dd)
        if(ps(igrid)%x(ix1,ix2,1) > xprobmin1+xedge .and. &
           ps(igrid)%x(ix1,ix2,1) < xprobmax1-xedge) dmi=max(dmi,dd)
      end do; end do
    end do
    ! rel = max abs(Q_pipe - Q_wave)/max abs(Q)
    gbuf = (/ dmx, qmx, dmi /)
    call MPI_ALLREDUCE(gbuf, gout, 3, MPI_DOUBLE_PRECISION, MPI_MAX, icomm, ierr)
    if(mype==0) write(*,'(a,es12.5,a,es12.5,a,es12.5,a,es12.5)') &
      ' rt_valpipe rel=',gout(1)/max(gout(2),1.0d-300),'  rel(x-interior,no seam)=', &
      gout(3)/max(gout(2),1.0d-300),'  max|dQ|=',gout(1),'  max|Q|=',gout(2)
  end subroutine rt_sc_validate_pipe

  !> List the blocks of this rank that lie in band ib, with their block position (bi, bj), cache chi
  !> and B on them, ghosts included, and check that the band is one complete level
  subroutine rt_band_blocks(ib)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    integer, intent(in) :: ib
    integer :: iigrid, igrid, iblk, ierr
    integer, allocatable :: cover(:,:)
    ! The band holds whole blocks, since rt_sc_setup snapped band_hi_c to a block face, and its
    ! longest diagonal p+q is rg_nbh + wv_nbxv
    wv_nbxv = band_nbv(ib)
    wv_dmax = rg_nbh + wv_nbxv
    if(allocated(wv_ig)) deallocate(wv_ig,wv_bi,wv_bj)
    allocate(wv_ig(max(igridstail,1)),wv_bi(max(igridstail,1)),wv_bj(max(igridstail,1)))
    wv_nblk = 0
    do iigrid=1,igridstail; igrid=igrids(iigrid)
      if(node(plevel_,igrid)/=rt_sc_amr_level) cycle
      if(rnode(rpxmin0_+updim,igrid) < band_lo_c(ib)-0.5d0*rg_dv) cycle
      if(rnode(rpxmin0_+updim,igrid)+rg_bv*rg_dv > band_hi_c(ib)+0.5d0*rg_dv) cycle
      wv_nblk = wv_nblk + 1
      wv_ig(wv_nblk) = igrid
      wv_bi(wv_nblk) = nint((rnode(rpxmin0_+hdim,igrid) -rg_hlo)/(rg_bh*rg_dh)) + 1
      wv_bj(wv_nblk) = nint((rnode(rpxmin0_+updim,igrid)-band_lo_c(ib))/(rg_bv*rg_dv)) + 1
    end do
    ! The caches only grow and are otherwise reused
    if(allocated(rtw_chib)) then
      if(size(rtw_chib,3) < max(wv_nblk,1)) deallocate(rtw_chib,rtw_bb,rtw_fx)
    end if
    if(.not.allocated(rtw_chib)) allocate(rtw_chib(ixGlo1:ixGhi1,ixGlo2:ixGhi2,max(wv_nblk,1)), &
                                          rtw_bb  (ixGlo1:ixGhi1,ixGlo2:ixGhi2,max(wv_nblk,1)), &
                                          rtw_fx  (ixGlo1:ixGhi1,ixGlo2:ixGhi2,max(wv_nblk,1)))
    ! Every (bi, bj) of the band must be owned exactly once, or the AMR pin has failed
    allocate(cover(rg_nbh, wv_nbxv)); cover = 0
    do iblk = 1, wv_nblk
      igrid = wv_ig(iblk)
      call rt_sc_chiB(ixG^LL, ixG^LL, ps(igrid)%w, ps(igrid)%x, rtw_chib(:,:,iblk), rtw_bb(:,:,iblk))
      ! Beyond the band ends chi and B take the edge values, as in the gather, so that the deep
      ! row thermalises. A ghost extrapolation would turn the large deep chi into spurious heating
      if(wv_bj(iblk)==1)       call copy_edge_row(iGlo_v, iMlo_v-1, iMlo_v)
      if(wv_bj(iblk)==wv_nbxv) call copy_edge_row(iMhi_v+1, iGhi_v, iMhi_v)
      cover(wv_bi(iblk), wv_bj(iblk)) = cover(wv_bi(iblk), wv_bj(iblk)) + 1
    end do
    ! The coverage is summed over the ranks before it is checked
    call MPI_ALLREDUCE(MPI_IN_PLACE, cover, rg_nbh*wv_nbxv, MPI_INTEGER, MPI_SUM, icomm, ierr)
    if(minval(cover) /= 1 .or. maxval(cover) /= 1) then
      if(mype==0) print *,'rt_band_blocks: band',ib,' not a complete single level ', &
           rt_sc_amr_level,' (mincover',minval(cover),' maxcover',maxval(cover),')'
      call mpistop("rt_band_blocks: RT band is not a complete single AMR level")
    end if
  contains
    !> Rows iv_lo to iv_hi of the chi and B caches of block iblk take row iv_src
    subroutine copy_edge_row(iv_lo, iv_hi, iv_src)
      integer, intent(in) :: iv_lo, iv_hi, iv_src
      integer :: ih, iv, i1, i2, j1, j2
      do iv = iv_lo, iv_hi
        do ih = iGlo_h, iGhi_h
          call hv(ih, iv, i1, i2)
          call hv(ih, iv_src, j1, j2)
          rtw_chib(i1,i2,iblk) = rtw_chib(j1,j2,iblk)
          rtw_bb(i1,i2,iblk)   = rtw_bb(j1,j2,iblk)
        end do
      end do
    end subroutine copy_edge_row
  end subroutine rt_band_blocks

  !> Tau at the top face of every block row of the current band, per global column along hdim.
  !> Each block adds the optical thickness of its columns, one allreduce assembles the table, and
  !> the top of a row lies below all the rows above it
  subroutine rt_tau_tops(tautop)
    use mod_global_parameters
    double precision, intent(out) :: tautop(rg_nh,wv_nbxv)
    double precision :: colthk(rg_nh,wv_nbxv), s
    integer :: iblk, ih, iv, i1, i2, ig, bj, jb, ierr
    colthk = 0.0d0
    do iblk = 1, wv_nblk
      do ih = iMlo_h, iMhi_h
        ig = (wv_bi(iblk)-1)*rg_bh + (ih - iMlo_h + 1)
        ! colthk = sum_j chi_j dy down the block column
        s = 0.0d0
        do iv = iMlo_v, iMhi_v
          call hv(ih, iv, i1, i2)
          s = s + rtw_chib(i1,i2,iblk)
        end do
        colthk(ig,wv_bj(iblk)) = s*rg_dv
      end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE, colthk, rg_nh*wv_nbxv, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr)
    do bj = 1, wv_nbxv
      ! tau_top(bj) = sum of colthk(jb) over the rows jb > bj
      tautop(:,bj) = 0.0d0
      do jb = bj+1, wv_nbxv; tautop(:,bj) = tautop(:,bj) + colthk(:,jb); end do
    end do
  end subroutine rt_tau_tops

  !> Tau from the band top to every cell centre of block iblk, integrated down from the top face
  !> of the block
  subroutine rt_block_tau(iblk, tautop, taub)
    use mod_global_parameters
    integer, intent(in)           :: iblk
    double precision, intent(in)  :: tautop(rg_nh,wv_nbxv)
    double precision, intent(out) :: taub(ixG^T)
    integer :: ih, iv, i1, i2, j1, j2, ig
    do ih = iMlo_h, iMhi_h
      ig = (wv_bi(iblk)-1)*rg_bh + (ih - iMlo_h + 1)
      ! tau = tau_top + chi dy/2 in the top cell,  tau_j = tau_(j+1) + (chi_(j+1) + chi_j) dy/2
      call hv(ih, iMhi_v, i1, i2)
      taub(i1,i2) = tautop(ig,wv_bj(iblk)) + 0.5d0*rtw_chib(i1,i2,iblk)*rg_dv
      do iv = iMhi_v-1, iMlo_v, -1
        call hv(ih, iv, i1, i2)
        call hv(ih, iv+1, j1, j2)
        taub(i1,i2) = taub(j1,j2) + 0.5d0*(rtw_chib(j1,j2,iblk)+rtw_chib(i1,i2,iblk))*rg_dv
      end do
    end do
  end subroutine rt_block_tau

  !> Prologue of both KBA schedules on level rt_sc_amr_level. It checks the setup, sets the
  !> x-distribution mode, refreshes the hydro ghosts, allocates the intensity state and the seam
  !> buffers, and builds the nslot wavefront datatype. Each rank sweeps only its own band blocks,
  !> ordered per ray by the diagonal p+q counted from the ray's upwind corner, with a getbc between
  !> diagonals, so every cell is swept once per ray
  subroutine rt_wave_begin(nslot)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    use mod_ghostcells_update, only: getbc, create_bc_mpi_datatype
    use mod_variables, only: iwstart, nwgc
    integer, intent(in) :: nslot
    if(.not.rt_sc_ready) call mpistop("rt_wave_begin: rt_sc_setup not called")
    if(i_isw < 0)        call mpistop("rt_wave_begin: wavefront slot not reserved")
    if(i_qrad < 0 .or. i_frad < 0) call mpistop("rt_wave_begin: the wave needs i_qrad and i_frad")
    ! With one block across x, sc_sweep_block wraps in place. With several, the x-upwind values
    ! come from the neighbour's ghosts, the diagonal widens, and the periodic seam takes the
    ! lagged wrap
    wv_xwrap = (rg_nbh == 1)
    if(.not.wv_xwrap .and. mype==0 .and. .not.rt_sc_seam_warned) then
      write(*,*) 'rt_wave_begin: periodic seam closed by the lagged wrap'
      rt_sc_seam_warned = .true.
    end if
    ! The hydro ghosts are refreshed so that chi and B can be evaluated on the ghost ring
    call getbc(global_time, 0.0d0, ps, iwstart, nwgc)
    call rt_sc_alloc_state()
    ! A non-periodic hdim has inflow edges and therefore no seam
    wv_seam = rt_sc_seam_on .and. .not.wv_xwrap .and. bc_h_periodic
    if(wv_seam) call rt_sc_seam_alloc(rg_nv, rt_sc_nray)
    ! From here on getbc of ps_rt exchanges only the nslot intensity columns
    call create_bc_mpi_datatype(i_isw, nslot)
  end subroutine rt_wave_begin

  !> Open band ib. It lists and caches the band blocks and zeroes the accumulators
  subroutine rt_wave_band_open(ib)
    use mod_global_parameters
    integer, intent(in) :: ib
    integer :: iblk, igrid
    ! Each band has its own incoming conditions at the ends
    bc_lo_thermal = (trim(rt_sc_bc_min(ib)) == 'thermal')
    bc_hi_thermal = (trim(rt_sc_bc_max(ib)) == 'thermal')
    call rt_band_blocks(ib)
    do iblk = 1, wv_nblk
      igrid = wv_ig(iblk)
      ! i_qrad accumulates J and later holds Q, while i_frad and rtw_fx accumulate Fy and Fx
      ps(igrid)%w(ixG^T,i_qrad) = 0.0d0
      ps(igrid)%w(ixG^T,i_frad) = 0.0d0
      rtw_fx(:,:,iblk) = 0.0d0
    end do
  end subroutine rt_wave_band_open

  !> Close the band. The fluxes are scaled by c_norm, Q is formed as the imbalance or as the -div F
  !> blend with the tau^2 taper, and the nslot wavefront datatype is restored
  subroutine rt_wave_band_close(nslot)
    use mod_global_parameters
    use mod_ghostcells_update, only: getbc, create_bc_mpi_datatype
    integer, intent(in) :: nslot
    integer :: iblk, igrid, ix1, ix2, ih, iv, i1, i2
    double precision :: taub(ixG^T), dfxdx, dfydy
    double precision, allocatable :: tautop(:,:)
    ! Fy = c sum_k w_k mu_v,k I_k,  Fx = c sum_k w_k mu_h,k I_k
    do iblk = 1, wv_nblk
      igrid = wv_ig(iblk)
      ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,i_frad) = c_norm*ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,i_frad)
      rtw_fx(ixMlo1:ixMhi1,ixMlo2:ixMhi2,iblk) = c_norm*rtw_fx(ixMlo1:ixMhi1,ixMlo2:ixMhi2,iblk)
    end do

    if(.not.rt_sc_flux_blend) then
      ! Without the blend Q is the imbalance, which is right in the thin limit but becomes stiff
      ! and cancelling deep down
      do iblk = 1, wv_nblk
        igrid = wv_ig(iblk)
        do ix2=ixMlo2,ixMhi2; do ix1=ixMlo1,ixMhi1
          ! Q = c chi (J - B)
          ps(igrid)%w(ix1,ix2,i_qrad) = c_norm*rtw_chib(ix1,ix2,iblk)* &
               (ps(igrid)%w(ix1,ix2,i_qrad) - rtw_bb(ix1,ix2,iblk))
        end do; end do
      end do
      return
    end if

    ! With the blend, Q = w c chi (J - B) + (1 - w)(-div F), where w = exp(-k tau) with tau
    ! measured from the top, and the taper multiplies the result. dFx/dx needs Fx in the x
    ! ghosts, which it reaches through the i_isw slot of ps_rt, and dFy/dy needs Fy in the y
    ! ghosts from a getbc of i_frad
    do iblk = 1, wv_nblk
      ps_rt(wv_ig(iblk))%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,i_isw) = rtw_fx(ixMlo1:ixMhi1,ixMlo2:ixMhi2,iblk)
    end do
    call create_bc_mpi_datatype(i_isw, 1)
    call rt_exchange(1)
    call create_bc_mpi_datatype(i_frad, 1)
    call getbc(global_time, 0.0d0, ps, i_frad, 1)
    allocate(tautop(rg_nh,wv_nbxv))
    call rt_tau_tops(tautop)
    do iblk = 1, wv_nblk
      igrid = wv_ig(iblk)
      call rt_block_tau(iblk, tautop, taub)
      do iv = iMlo_v, iMhi_v
        do ih = iMlo_h, iMhi_h
          call hv(ih, iv, i1, i2)
          ! dFx/dx = (Fx_(i+1) - Fx_(i-1))/(2 dx), centred, since the ghosts hold the neighbour
          ! block or the periodic image
          dfxdx = (ps_rt(igrid)%w(i1+eh(1),i2+eh(2),i_isw) &
                  -ps_rt(igrid)%w(i1-eh(1),i2-eh(2),i_isw))/(2.0d0*rg_dh)
          ! dFy/dy = (Fy_(j+1) - Fy_(j-1))/(2 dy), one-sided in the bottom and top rows of the
          ! band, where the ghosts lie outside it
          if(wv_bj(iblk)==1 .and. iv==iMlo_v) then
            dfydy = (ps(igrid)%w(i1+ev(1),i2+ev(2),i_frad)-ps(igrid)%w(i1,i2,i_frad))/rg_dv
          else if(wv_bj(iblk)==wv_nbxv .and. iv==iMhi_v) then
            dfydy = (ps(igrid)%w(i1,i2,i_frad)-ps(igrid)%w(i1-ev(1),i2-ev(2),i_frad))/rg_dv
          else
            dfydy = (ps(igrid)%w(i1+ev(1),i2+ev(2),i_frad) &
                    -ps(igrid)%w(i1-ev(1),i2-ev(2),i_frad))/(2.0d0*rg_dv)
          end if
          ps(igrid)%w(i1,i2,i_qrad) = blend_heating(rtw_chib(i1,i2,iblk), &
               ps(igrid)%w(i1,i2,i_qrad), rtw_bb(i1,i2,iblk), dfxdx+dfydy, taub(i1,i2))
        end do
      end do
    end do
    deallocate(tautop)
    ! The datatype goes back to the nslot intensity columns for the next band
    call create_bc_mpi_datatype(i_isw, nslot)
  end subroutine rt_wave_band_close

  !> Add the weighted intensity iib of ray k on block iblk to J in i_qrad, Fy in i_frad and Fx in
  !> rtw_fx
  subroutine rt_add_moments(iblk, k, iib)
    use mod_global_parameters
    integer, intent(in) :: iblk, k
    double precision, intent(in) :: iib(ixG^T)
    integer :: igrid
    igrid = wv_ig(iblk)
    ! J += w_k I_k,  Fy += w_k mu_v,k I_k,  Fx += w_k mu_h,k I_k
    ps(igrid)%w(ixM^T,i_qrad) = ps(igrid)%w(ixM^T,i_qrad) + wmu(k)*iib(ixM^T)
    ps(igrid)%w(ixM^T,i_frad) = ps(igrid)%w(ixM^T,i_frad) + wmu(k)*mu_v(k)*iib(ixM^T)
    rtw_fx(ixM^T,iblk) = rtw_fx(ixM^T,iblk) + wmu(k)*mu_h(k)*iib(ixM^T)
  end subroutine rt_add_moments

  !> Ghost exchange of the nslot intensity columns of ps_rt from i_isw, with the physical
  !> boundaries left alone, since the inflows are written separately
  subroutine rt_exchange(nslot)
    use mod_global_parameters
    use mod_ghostcells_update, only: getbc, bcphys
    integer, intent(in) :: nslot
    bcphys = .false.
    call getbc(global_time, 0.0d0, ps_rt, i_isw, nslot)
    bcphys = .true.
  end subroutine rt_exchange

  !> Epilogue. It hands the seam edges over to the next solve, restores the hydro datatype and
  !> releases the block list
  subroutine rt_wave_end()
    use mod_global_parameters
    use mod_ghostcells_update, only: create_bc_mpi_datatype
    use mod_variables, only: iwstart, nwgc
    if(wv_seam) call rt_sc_seam_finalize(rg_nv, rt_sc_nray)
    call create_bc_mpi_datatype(iwstart, nwgc)
    if(allocated(wv_ig)) deallocate(wv_ig,wv_bi,wv_bj)
  end subroutine rt_wave_end

  !> Seam read. At the domain edge where ray k enters, the wrapped upwind ghost column of the block
  !> takes the outflow of the far edge from the last solve. A ray moving in +x enters through the
  !> first block column and reads seamR, and one moving in -x enters through the last and reads
  !> seamL
  subroutine rt_seam_read(iib, igrid, ibi, sx, k)
    use mod_global_parameters
    double precision, intent(inout) :: iib(ixG^T)
    integer, intent(in) :: igrid, ibi, sx, k
    integer :: iv, i1, i2, gy0
    ! gy0 is the global row offset of the block along updim
    gy0 = nint((rnode(rpxmin0_+updim,igrid)-rg_vlo)/rg_dv)
    if(ibi==1 .and. sx>0) then
      do iv = iMlo_v, iMhi_v
        call hv(iMlo_h-1, iv, i1, i2); iib(i1,i2) = rt_sc_seamR_prev(gy0+(iv-iMlo_v)+1,k)
      end do
    else if(ibi==rg_nbh .and. sx<0) then
      do iv = iMlo_v, iMhi_v
        call hv(iMhi_h+1, iv, i1, i2); iib(i1,i2) = rt_sc_seamL_prev(gy0+(iv-iMlo_v)+1,k)
      end do
    end if
  end subroutine rt_seam_read

  !> Seam write. At the far domain edge the outflow of ray k is kept for the next solve. A ray
  !> moving in +x leaves through the last block column into seamR, and one moving in -x leaves
  !> through the first into seamL
  subroutine rt_seam_write(iib, igrid, ibi, sx, k)
    use mod_global_parameters
    double precision, intent(in) :: iib(ixG^T)
    integer, intent(in) :: igrid, ibi, sx, k
    integer :: iv, i1, i2, gy0
    gy0 = nint((rnode(rpxmin0_+updim,igrid)-rg_vlo)/rg_dv)
    if(ibi==rg_nbh .and. sx>0) then
      do iv = iMlo_v, iMhi_v
        call hv(iMhi_h, iv, i1, i2); rt_sc_seamR_cur(gy0+(iv-iMlo_v)+1,k) = iib(i1,i2)
      end do
    else if(ibi==1 .and. sx<0) then
      do iv = iMlo_v, iMhi_v
        call hv(iMlo_h, iv, i1, i2); rt_sc_seamL_cur(gy0+(iv-iMlo_v)+1,k) = iib(i1,i2)
      end do
    end if
  end subroutine rt_seam_write

  !> Inflow intensity in slot islot. The ghosts beyond the band ends, and beyond the horizontal
  !> edges where hdim is not periodic, take B of the edge cell when that end is thermal and zero
  !> when it is open. The band ends are interior rows of the domain, so getbc never fills them
  subroutine rt_inflow_ghosts(islot)
    use mod_global_parameters
    integer, intent(in) :: islot
    integer :: iblk, igrid
    do iblk = 1, wv_nblk
      igrid = wv_ig(iblk)
      if(wv_bj(iblk)==1)       call fill_rows(iGlo_v, iMlo_v-1, iMlo_v, bc_lo_thermal)
      if(wv_bj(iblk)==wv_nbxv) call fill_rows(iMhi_v+1, iGhi_v, iMhi_v, bc_hi_thermal)
    end do
    if(bc_h_periodic) return
    do iblk = 1, wv_nblk
      igrid = wv_ig(iblk)
      if(wv_bi(iblk)==1)      call fill_columns(iGlo_h, iMlo_h-1, iMlo_h, bc_hlo_thermal)
      if(wv_bi(iblk)==rg_nbh) call fill_columns(iMhi_h+1, iGhi_h, iMhi_h, bc_hhi_thermal)
    end do
  contains
    !> Ghost rows iv_lo to iv_hi, from the edge row iv_src
    subroutine fill_rows(iv_lo, iv_hi, iv_src, thermal)
      integer, intent(in) :: iv_lo, iv_hi, iv_src
      logical, intent(in) :: thermal
      integer :: ih, iv, i1, i2, j1, j2
      do iv = iv_lo, iv_hi
        do ih = iGlo_h, iGhi_h
          call hv(ih, iv, i1, i2)
          call hv(ih, iv_src, j1, j2)
          ps_rt(igrid)%w(i1,i2,islot) = merge(rtw_bb(j1,j2,iblk), 0.0d0, thermal)
        end do
      end do
    end subroutine fill_rows
    !> Ghost columns ih_lo to ih_hi, from the edge column ih_src
    subroutine fill_columns(ih_lo, ih_hi, ih_src, thermal)
      integer, intent(in) :: ih_lo, ih_hi, ih_src
      logical, intent(in) :: thermal
      integer :: ih, iv, i1, i2, j1, j2
      do ih = ih_lo, ih_hi
        do iv = iGlo_v, iGhi_v
          call hv(ih, iv, i1, i2)
          call hv(ih_src, iv, j1, j2)
          ps_rt(igrid)%w(i1,i2,islot) = merge(rtw_bb(j1,j2,iblk), 0.0d0, thermal)
        end do
      end do
    end subroutine fill_columns
  end subroutine rt_inflow_ghosts

  !> Sweep block iblk for ray k in place on iib. With updim = 1 the block is transposed around the
  !> kernel, which expects its first index to be horizontal. The sweep time of each block is
  !> accumulated for the load balancer
  subroutine rt_sweep_one(iblk, k, iib)
    use mod_global_parameters
    integer, intent(in) :: iblk, k
    double precision, intent(inout) :: iib(ixG^T)
    double precision :: tii(ixGlo2:ixGhi2,ixGlo1:ixGhi1), t0
    integer :: igrid
    if(lb_automatic) t0 = MPI_WTIME()
    igrid = wv_ig(iblk)
    if(updim==2) then
      call sc_sweep_block(mu_h(k), mu_v(k), rg_dh, rg_dv, wv_xwrap, ixG^LL, ixM^LL, &
                          rtw_chib(:,:,iblk), rtw_bb(:,:,iblk), iib)
    else
      tii = transpose(iib)
      call sc_sweep_block(mu_h(k), mu_v(k), rg_dh, rg_dv, wv_xwrap, &
                          ixGlo2,ixGlo1,ixGhi2,ixGhi1, ixMlo2,ixMlo1,ixMhi2,ixMhi1, &
                          transpose(rtw_chib(:,:,iblk)), transpose(rtw_bb(:,:,iblk)), tii)
      iib = transpose(tii)
    end if
    if(lb_automatic) block_cost_rt(igrid) = block_cost_rt(igrid) + (MPI_WTIME() - t0)
  end subroutine rt_sweep_one

  !> Position of a block counted from the upwind end of a row of n blocks, for a ray whose cosine
  !> along that row has sign s. With s = 0 no block depends on its neighbour along the row, so
  !> every block sits at the first position
  pure integer function wave_stage(b, n, s)
    integer, intent(in) :: b, n, s
    wave_stage = merge(b, merge(n+1-b, 1, s<0), s>0)
  end function wave_stage

  !> Ray-by-ray KBA, with one ray at a time passing through the diagonal wavefront
  subroutine rt_sc_solve_wave()
    use mod_global_parameters
    integer :: ib, k, sx, sy, d, iblk, igrid
    double precision :: iib(ixG^T)
    call rt_wave_begin(1)
    do ib = 1, rt_sc_nband
      call rt_wave_band_open(ib)
      do k = 1, rt_sc_nray
        ! The signs of the cosines fix the upwind corner the wavefront starts from
        sx = ray_sign(mu_h(k)); sy = ray_sign(mu_v(k))
        ! The intensity is zeroed for every ray, as in the gather, so the rays stay independent
        do iblk = 1, wv_nblk
          ps_rt(wv_ig(iblk))%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,i_isw) = 0.0d0
        end do
        ! Diagonal d holds the blocks with p + q = d, starting from the upwind corner at d = 2.
        ! Before each diagonal the getbc hands the outflow faces of the last one to their
        ! downstream neighbours, and the band-end and edge inflows are written over its ghosts
        do d = 2, wv_dmax
          call rt_exchange(1)
          call rt_inflow_ghosts(i_isw)
          do iblk = 1, wv_nblk
            if(wave_stage(wv_bi(iblk),rg_nbh,sx) + wave_stage(wv_bj(iblk),wv_nbxv,sy) /= d) cycle
            igrid = wv_ig(iblk)
            ! The block is swept on a copy that holds its ghosts, with the seam inflow added at the
            ! domain edge, and its interior is written back
            iib = ps_rt(igrid)%w(:,:,i_isw)
            if(wv_seam) call rt_seam_read(iib, igrid, wv_bi(iblk), sx, k)
            call rt_sweep_one(iblk, k, iib)
            ps_rt(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,i_isw) = iib(ixMlo1:ixMhi1,ixMlo2:ixMhi2)
            if(wv_seam) call rt_seam_write(iib, igrid, wv_bi(iblk), sx, k)
          end do
        end do
        ! Once the ray has crossed the whole band, its weighted intensity adds to J, Fy and Fx
        do iblk = 1, wv_nblk
          call rt_add_moments(iblk, k, ps_rt(wv_ig(iblk))%w(:,:,i_isw))
        end do
      end do
      call rt_wave_band_close(1)
    end do
    call rt_wave_end()
  end subroutine rt_sc_solve_wave

  !> Pipelined KBA. The direction-sign corners run together in batches of nisw slots, and within
  !> each corner the rays follow one another through its wavefront. The result agrees with that
  !> of rt_sc_solve_wave to round-off, since the moments are summed in another order
  subroutine rt_sc_solve_wave_pipe()
    use mod_global_parameters
    use mod_ghostcells_update, only: create_bc_mpi_datatype
    integer :: ib, k, sx, sy, iblk, igrid, islot, nslot
    integer :: ng, gidx, g, jj, tt_clk, tmax, gbase, gend, gg
    integer :: grp_sx(9), grp_sy(9), grp_n(9)
    integer, allocatable :: grp_ray(:,:)
    double precision :: iib(ixG^T)

    ! The rays are grouped by the signs of (mu_h, mu_v), since the rays of one group share an
    ! upwind corner and therefore one schedule. Counting the zero sign there are nine sign pairs,
    ! so no more groups can occur
    allocate(grp_ray(9, max(rt_sc_nray,1)))
    ng = 0
    do k = 1, rt_sc_nray
      sx = ray_sign(mu_h(k)); sy = ray_sign(mu_v(k))
      gidx = 0
      do g = 1, ng
        if(grp_sx(g)==sx .and. grp_sy(g)==sy) then; gidx = g; exit; end if
      end do
      if(gidx==0) then
        ng = ng + 1
        grp_sx(ng)=sx; grp_sy(ng)=sy; grp_n(ng)=0; gidx=ng
      end if
      grp_n(gidx) = grp_n(gidx) + 1
      grp_ray(gidx, grp_n(gidx)) = k
    end do
    call rt_wave_begin(nisw)
    do ib = 1, rt_sc_nband
      call rt_wave_band_open(ib)
      ! In group gg a block sits at stage p+q, and at clock tt_clk it sweeps ray
      ! jj = tt_clk-(p+q)+2 of the group into slot i_isw+(gg-gbase). One getbc per clock
      ! advances every front of the batch by one stage. The accumulators were zeroed when the band
      ! opened, so they sum over the batches
      do gbase = 1, ng, nisw
        gend  = min(gbase+nisw-1, ng)
        nslot = gend-gbase+1
        ! The slots of the batch start from zero intensity
        do iblk = 1, wv_nblk
          ps_rt(wv_ig(iblk))%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,i_isw:i_isw+nisw-1) = 0.0d0
        end do
        ! The last ray of a group reaches the last diagonal at clock grp_n + wv_dmax - 2, and the
        ! batch runs until its longest group has finished
        tmax = wv_dmax - 2
        do gg = gbase, gend; tmax = max(tmax, grp_n(gg) + wv_dmax - 2); end do
        ! The datatype is built once per batch, for its number of slots
        call create_bc_mpi_datatype(i_isw, nslot)
        do tt_clk = 1, tmax
          ! One getbc moves every front, and the band-end and edge inflows are then written into
          ! each slot
          call rt_exchange(nslot)
          do gg = gbase, gend
            call rt_inflow_ghosts(i_isw+(gg-gbase))
          end do
          do gg = gbase, gend
            islot = i_isw+(gg-gbase)
            do iblk = 1, wv_nblk
              ! jj = t - (p + q) + 2, the ray of the group that reaches this block at clock t. A
              ! block sweeps at this clock only when that is a ray of the group
              jj = tt_clk - (wave_stage(wv_bi(iblk),rg_nbh,grp_sx(gg)) &
                             + wave_stage(wv_bj(iblk),wv_nbxv,grp_sy(gg))) + 2
              if(jj < 1 .or. jj > grp_n(gg)) cycle
              k = grp_ray(gg, jj)
              igrid = wv_ig(iblk)
              ! The copy keeps the upwind ghosts from getbc, and its interior is cleared for ray k
              iib = ps_rt(igrid)%w(:,:,islot)
              iib(ixMlo1:ixMhi1,ixMlo2:ixMhi2) = 0.0d0
              if(wv_seam) call rt_seam_read(iib, igrid, wv_bi(iblk), grp_sx(gg), k)
              call rt_sweep_one(iblk, k, iib)
              ! The interior goes back into the slot, the seam edge is kept, and the weighted
              ! intensity adds to J, Fy and Fx straight away
              ps_rt(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,islot) = iib(ixMlo1:ixMhi1,ixMlo2:ixMhi2)
              if(wv_seam) call rt_seam_write(iib, igrid, wv_bi(iblk), grp_sx(gg), k)
              call rt_add_moments(iblk, k, iib)
            end do
          end do
        end do
      end do
      call rt_wave_band_close(nisw)
    end do
    call rt_wave_end()
  end subroutine rt_sc_solve_wave_pipe
  }
end module mod_rt_sc_solvers
