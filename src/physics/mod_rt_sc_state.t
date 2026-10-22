!=============================================================================
!> Parameters, resolved state and lifecycle of the short-characteristics (SC) radiative
!> transfer, the Olson-Kunasz linear RT of AMRVAC
!>
!> The &rt_list parameters are physical, with heights in cm, pressures in dyn/cm^2 and
!> temperatures in K. rt_sc_params_read reads them, and rt_sc_setup, run once the units are
!> fixed, converts them to code units, snaps the top of each band to a block face, sets the
!> horizontal boundary conditions, binds the opacity callback and builds the quadrature.
!> rt_sc_activate, called by the physics module when hd_radiative_transfer or
!> mhd_radiative_transfer is set, runs the reading and reserves the intensity slots that the
!> distributed wavefront needs
!>
!> The module also holds the state that the solvers share. This is the ps_rt intensity state,
!> kept out of the hydro w so that it skips the stage copies, the prolongation and the output,
!> together with the block caches of chi and B, the KBA wave state of the current solve, the
!> buffers of the periodic seam and the node communicators of the gather. It uses
!> mod_rt_sc_opacity and mod_rt_sc_kernel, and mod_rt_sc_solvers uses it in turn
!>
!> Jack Jenkins (10.2026); More or less finished in v1 state for amrvac4.0
!=============================================================================
module mod_rt_sc_state
  use mod_global_parameters, only: bigdouble, smalldouble, std_len, max_blocks
  use mod_physicaldata, only: state
  use mod_rt_sc_opacity
  use mod_rt_sc_kernel
  implicit none

  !> True when the physics module has switched the RT on, with hd_radiative_transfer or
  !> mhd_radiative_transfer, through rt_sc_activate
  logical, public :: rt_sc_active = .false.
  !> Axis of the stratification, 1 for x and 2 for y. The default of -1 takes the last axis
  integer, public :: rt_sc_updim = -1
  !> The AMR level the RT runs on. The bands are held at this level so that the RT sees a uniform
  !> grid, and it cannot exceed refine_max_level
  integer, public :: rt_sc_amr_level = 1
  !> Centre and half-width in K of the temperature gate in rt_sc_thin_weight, which divides the
  !> radiative losses between the SC transfer and the thin cooling
  double precision, public :: rt_sc_blend_T0 = 2.0d4
  double precision, public :: rt_sc_blend_dT = 5.0d3

  !> The RT bands along updim, of which a bar with a photosphere at each end needs two. The edges
  !> are in cm, with -1e30 and +1e30 standing for the domain ends, and each end has its own
  !> incoming-intensity condition, 'thermal' or 'open'
  integer, parameter :: rt_sc_maxband = 8
  integer, public :: rt_sc_nband = 1
  double precision, public :: rt_sc_band_min(rt_sc_maxband) = -1.0d30
  double precision, public :: rt_sc_band_max(rt_sc_maxband) =  1.0d30
  character(len=10), public :: rt_sc_bc_min(rt_sc_maxband) = 'thermal'
  character(len=10), public :: rt_sc_bc_max(rt_sc_maxband) = 'open'
  !> Incoming intensity at the horizontal edges where they are not periodic, 'open' for I = 0 or
  !> 'thermal' for I = B of the edge column
  character(len=10), public :: rt_sc_bc_hmin = 'open', rt_sc_bc_hmax = 'open'
  !> rt_sc_flux_blend blends the thin heating c chi (J - B) with the thick -div F. rt_sc_cfl limits
  !> the heating so that one stage changes the gas energy by at most that fraction of itself
  logical, public :: rt_sc_flux_blend = .true.
  double precision, public :: rt_sc_cfl = 0.5d0
  !> Prints the cumulative RT timing every 10 solves
  logical, public :: rt_sc_debug = .false.
  !> Not read by the RT itself. The user files set it from eos%eint_rho_FI_threshold for their own
  !> SC source
  double precision, public :: rt_sc_fi_thresh = -1.0d0
  !> Extra-variable indices of the heating Q, the vertical flux F and the emergent intensity Iout,
  !> left at -1 when they are not registered
  integer, public :: i_qrad = -1, i_frad = -1, i_iout = -1
  !> Horizontal mean of the flux on the tau = rt_sc_tau_emergent surface, from rt_sc_emergent_tau1.
  !> A flux controller should read this one, because Frad still grows by 10-16% above tau = 1 and
  !> the flux at the band top therefore reads high
  double precision, public :: rt_sc_Ftau1 = 0.0d0
  double precision, public :: rt_sc_tau_emergent = 1.0d0
  !> The distributed KBA wavefront. Each rank sweeps only its own band blocks, one diagonal at a
  !> time, with a getbc of the ray intensity between diagonals, so no rank ever holds the whole
  !> slab. Its intensity slots are reserved by rt_sc_activate. Setting it false selects the
  !> replicated gather, the original implementation, kept for comparison only: it holds the whole
  !> band on every rank and its periodic seam is not exact. Production runs should keep KBA
  logical, public :: rt_sc_distributed = .true.
  !> The pipelined KBA. The four sweep corners run at once, each with its own intensity slot, and
  !> the rays of a corner follow one another through the wavefront. It agrees with the ray-by-ray
  !> wave to round-off and is 1.4-1.8 times faster
  logical, public :: rt_sc_pipeline = .true.
  !> For validation, runs the ray-by-ray and pipelined waves on the same state and prints max|dQ|
  logical, public :: rt_sc_valpipe = .false.
  !> First column of the ray intensity in ps_rt, and -1 until reserved. It is the only field getbc
  !> exchanges during a sweep, because chi and B are cached per block while J and Fy accumulate in
  !> i_qrad and i_frad, perhaps confusingly but it works and prevents size creep in ps_rt
  integer :: i_isw = -1
  !> Number of intensity slots, one for the ray-by-ray wave and up to four when pipelined
  integer :: nisw = 1
  !> The RT intensity state. It is kept out of the hydro w so that it skips the stage copies, the
  !> prolongation and the output, and getbc exchanges it with bcphys off. rt_sc_alloc_state
  !> allocates it per block
  type(state), dimension(:), allocatable, target :: ps_rt
  ! Set once the seam notice has been printed, so that it appears only once
  logical :: rt_sc_seam_warned = .false.
  !> The periodic seam of the single-pass KBA. A ray crossing the seam takes as its inflow the
  !> outflow at the far edge from the previous solve, which is exact at steady state. The buffers
  !> are indexed by global row along updim and by ray. seamR holds the right-hand outflow of the
  !> last block column, which feeds the +x wrap into the first, and seamL the reverse. The _prev
  !> buffers are read and the _cur buffers written, the allreduce at the end of each solve makes
  !> them global, and both start at zero. With rt_sc_seam_on false the seam inflow stays zero
  logical :: rt_sc_seam_on = .true.
  double precision, allocatable :: rt_sc_seamL_prev(:,:), rt_sc_seamR_prev(:,:)
  double precision, allocatable :: rt_sc_seamL_cur(:,:),  rt_sc_seamR_cur(:,:)
  !> Node-local communicator for the moment reduction of the gather. The ranks reach that
  !> reduction far out of step because of the ray imbalance, and inter-node collectives under
  !> such a stagger fail on the tcp/ofi provider. Each node therefore forms the full moments from
  !> its own copy of the field
  integer :: rt_node_comm = -1, rt_node_rank = 0, rt_node_size = 1
  !> Communicator of one leader rank per node. The field gather sums within each node, then over
  !> the leaders, and then broadcasts within the node, so that no collective spans all the ranks.
  !> The tcp provider refuses those at 46 Mm
  integer :: rt_lead_comm = -1
  !> Cumulative wall time of the field gather, the sweep, the moment reduction and the scatter,
  !> and of the wave, together with the number of solves
  double precision, public :: rt_wt_gath=0.d0, rt_wt_swp=0.d0, rt_wt_red=0.d0, rt_wt_scat=0.d0
  double precision, public :: rt_wt_wave=0.d0
  integer, public :: rt_nsolve=0
  !> Block chi and B including the ghosts, rebuilt every solve, and the block horizontal flux Fx
  !> for the -div F blend. Fy accumulates in i_frad instead, and Fx is copied into ps_rt for its
  !> ghost exchange when the band closes
  double precision, allocatable :: rtw_chib(:,:,:), rtw_bb(:,:,:), rtw_fx(:,:,:)
  !> The band blocks of this rank with their block position (bi, bj), the block rows of the
  !> current band and its longest diagonal, and the wrap and seam modes of the KBA
  integer :: wv_nblk, wv_nbxv, wv_dmax
  integer, allocatable :: wv_ig(:), wv_bi(:), wv_bj(:)
  logical :: wv_xwrap
  ! True when the lagged wrap applies, which needs the distributed x-sweep and a periodic domain
  logical :: wv_seam
  !> Ray cosines along the RT grid axes (hdim, updim). They are (mux, muy) for updim = 2 and the
  !> swapped pair for updim = 1
  double precision, allocatable :: mu_h(:), mu_v(:)

  ! The axes, band edges and readiness flag, all resolved in code units by rt_sc_setup
  integer :: updim = 2, hdim = 1
  double precision :: band_lo_c(rt_sc_maxband), band_hi_c(rt_sc_maxband)
  logical :: rt_sc_ready = .false.
  !> The RT grid on rt_sc_amr_level along (hdim, updim): cells across the domain, cells per block,
  !> blocks across the width, cell sizes and the lower domain edges, with the cell and block
  !> rows of each band
  integer :: rg_nh, rg_nv, rg_bh, rg_bv, rg_nbh
  double precision :: rg_dh, rg_dv, rg_hlo, rg_vlo
  integer :: band_nv(rt_sc_maxband), band_nbv(rt_sc_maxband)
  !> Unit steps of hdim and updim in block index space and the block index ranges along them, so
  !> that a loop written in (h, v) serves either stratification axis
  integer :: eh(2), ev(2)
  integer :: iMlo_h, iMhi_h, iMlo_v, iMhi_v, iGlo_h, iGhi_h, iGlo_v, iGhi_v

  namelist /rt_list/ rt_sc_updim, rt_sc_amr_level, rt_sc_nband, &
       rt_sc_band_min, rt_sc_band_max, rt_sc_bc_min, rt_sc_bc_max, &
       rt_sc_bc_hmin, rt_sc_bc_hmax, &
       rt_sc_opacity, rt_sc_kappa_table, rt_sc_Tcap, rt_sc_dTcap, rt_sc_p_TR, rt_sc_dp_TR, &
       rt_sc_quad, rt_sc_ndir, rt_sc_nmu, rt_sc_flux_blend, rt_sc_blend_k, rt_sc_tau_min, rt_sc_cfl, &
       rt_sc_blend_T0, rt_sc_blend_dT, rt_sc_debug, &
       rt_sc_distributed, rt_sc_pipeline, rt_sc_valpipe, rt_sc_seam_on

  public :: rt_sc_params_read
  public :: rt_sc_setup
  public :: rt_sc_reserve_slots
  public :: rt_sc_init

contains

  !> Read &rt_list from each parameter file. The values stay physical until rt_sc_setup, and a
  !> file without the namelist is passed over
  subroutine rt_sc_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer :: n
    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status="old")
      read(unitpar, rt_list, end=222)
222   close(unitpar)
    end do
  end subroutine rt_sc_params_read

  !> Check the geometry, convert the parameters to code units, resolve the bands and the boundary
  !> conditions, and bind rt_sc_chiB to the chosen opacity. Runs after the unit system is fixed
  subroutine rt_sc_setup()
    use mod_global_parameters
    use mod_geometry, only: coordinate, Cartesian
    use mod_usr_methods, only: usr_rt_sc_chiB
    use mod_comm_lib, only: mpistop
    integer :: ib, fac
    integer :: dnx(ndim), bnx(ndim), mlo(ndim), mhi(ndim), glo(ndim), ghi(ndim)
    double precision :: dlo(ndim), dhi(ndim), hi_req
    if(.not. rt_sc_active) return
    if(rt_sc_updim < 1) rt_sc_updim = ndim
    ! The RT grid arithmetic assumes uniform Cartesian cells. All else is rejected for now
    if(coordinate /= Cartesian) then
      call mpistop("rt_sc_setup: uniform Cartesian only (no stretched grids)")
    end if
    if(ndim /= 2) call mpistop("rt_sc_setup: 2D only")
    if(rt_sc_updim > ndim) call mpistop("rt_sc_setup: bad rt_sc_updim")
    ! The horizontal axis is whichever of the two is not updim
    updim = rt_sc_updim
    hdim  = 3 - updim
    if(rt_sc_amr_level < 1) rt_sc_amr_level = 1
    if(rt_sc_amr_level > refine_max_level) call mpistop("rt_sc_setup: rt_sc_amr_level above refine_max_level")
    if(rt_sc_nband < 1) rt_sc_nband = 1
    if(rt_sc_nband > rt_sc_maxband) call mpistop("rt_sc_setup: rt_sc_nband above rt_sc_maxband")
    ! The distributed wave keeps its intensity in slots that rt_sc_activate reserves
    if(rt_sc_distributed .and. i_isw < 0) call mpistop("rt_sc_setup: call rt_sc_activate first")
    ! The RT grid along (hdim, updim). The horizontal direction splits into whole blocks, since
    ! AMRVAC requires domain_nx to be a multiple of block_nx
    {dlo(^D)=xprobmin^D; dhi(^D)=xprobmax^D; dnx(^D)=domain_nx^D; bnx(^D)=block_nx^D;}
    fac    = 2**(rt_sc_amr_level-1)
    rg_nh  = dnx(hdim)*fac;  rg_nv = dnx(updim)*fac
    rg_bh  = bnx(hdim);      rg_bv = bnx(updim)
    rg_nbh = rg_nh/rg_bh
    rg_hlo = dlo(hdim);      rg_vlo = dlo(updim)
    rg_dh  = (dhi(hdim) -dlo(hdim)) /dble(rg_nh)
    rg_dv  = (dhi(updim)-dlo(updim))/dble(rg_nv)
    {mlo(^D)=ixMlo^D; mhi(^D)=ixMhi^D; glo(^D)=ixGlo^D; ghi(^D)=ixGhi^D;}
    eh = 0; eh(hdim) = 1
    ev = 0; ev(updim) = 1
    iMlo_h = mlo(hdim);  iMhi_h = mhi(hdim);  iGlo_h = glo(hdim);  iGhi_h = ghi(hdim)
    iMlo_v = mlo(updim); iMhi_v = mhi(updim); iGlo_v = glo(updim); iGhi_v = ghi(updim)
    ! The band edges go to code units, and the 1e30 sentinels resolve to the domain ends
    do ib = 1, rt_sc_nband
      band_lo_c(ib) = merge(dlo(updim), rt_sc_band_min(ib)/unit_length, rt_sc_band_min(ib) <= -1.0d29)
      band_hi_c(ib) = merge(dhi(updim), rt_sc_band_max(ib)/unit_length, rt_sc_band_max(ib) >=  1.0d29)
      if(band_hi_c(ib) <= band_lo_c(ib)) call mpistop("rt_sc_setup: band top not above its bottom")
    end do
    ! The KBA needs every band to tile into whole blocks along updim, so its top is snapped down to
    ! the last block edge below the requested height, and everything downstream reads band_hi_c.
    ! It never moves up, which could leave the domain or reach blocks the refinement pin has not
    ! refined
    do ib = 1, rt_sc_nband
      hi_req = band_hi_c(ib)
      ! Integer division counts only the whole blocks that fit in the requested band
      ! nbv = floor(ncells/bv),  band_hi = band_lo + nbv bv dv
      band_nbv(ib) = nint((band_hi_c(ib)-band_lo_c(ib))/rg_dv)/rg_bv
      if(band_nbv(ib) < 1) call mpistop("rt_sc_setup: band thinner than one block along updim")
      band_nv(ib)   = band_nbv(ib)*rg_bv
      band_hi_c(ib) = band_lo_c(ib) + dble(band_nv(ib))*rg_dv
      if(mype==0 .and. dabs(hi_req-band_hi_c(ib)) > 1.0d-10*dabs(rg_dv)) then
        write(*,'(a,i2,a,es13.6,a,es13.6,a,i0,a,i0,a)') &
          ' rt_sc_setup: band ', ib, ' top snapped ', hi_req*unit_length, &
          ' -> ', band_hi_c(ib)*unit_length, ' cm  (', band_nbv(ib), ' blocks of ', rg_bv, ' cells)'
      end if
    end do
    ! The switch pressure is held as a code-unit logarithm, so the cells need no conversion. The
    ! opacity callback is bound once here, and the table is read only when it is needed
    rt_sc_lp_TR = dlog10(rt_sc_p_TR/unit_pressure)
    select case(trim(rt_sc_opacity))
    case('Rosseland')
      call read_opacity_table()
      rt_sc_chiB => rt_sc_chiB_rosseland
    case('Kramers')
      rt_sc_chiB => rt_sc_chiB_kramers
    case('usr')
      if(.not. associated(usr_rt_sc_chiB)) &
        call mpistop("rt_sc_setup: rt_sc_opacity='usr' but usr_rt_sc_chiB not set in mod_usr.t")
      rt_sc_chiB => usr_rt_sc_chiB
    case default
      call mpistop("rt_sc_setup: unknown rt_sc_opacity (Rosseland, Kramers or usr)")
    end select
    ! The horizontal edges wrap only where the domain is periodic, and otherwise take the inflow
    ! conditions given in &rt_list
    bc_h_periodic  = periodB(hdim)
    bc_hlo_thermal = (trim(rt_sc_bc_hmin) == 'thermal')
    bc_hhi_thermal = (trim(rt_sc_bc_hmax) == 'thermal')
    if(mype==0 .and. .not.bc_h_periodic) print *, &
      'rt_sc_setup: hdim not periodic, inflow BC ', trim(rt_sc_bc_hmin),'/',trim(rt_sc_bc_hmax)

    ! The GL polar nodes are placed about the stratification axis, and the cosines are then mapped
    ! onto the RT grid axes
    rt_sc_polar_axis = updim
    call rt_sc_init()
    call rt_sc_axis_cosines()
    rt_sc_ready = .true.
    if(mype==0) print *, "rt_sc_setup: done. updim=",updim," level=",rt_sc_amr_level, &
         " nband=",rt_sc_nband," opacity=",trim(rt_sc_opacity)
  end subroutine rt_sc_setup

  !> Map the quadrature onto the RT grid axes (hdim, updim). Called after rt_sc_init
  subroutine rt_sc_axis_cosines()
    if(allocated(mu_h)) deallocate(mu_h, mu_v)
    allocate(mu_h(rt_sc_nray), mu_v(rt_sc_nray))
    if(hdim == 1) then
      mu_h = mux; mu_v = muy
    else
      mu_h = muy; mu_v = mux
    end if
  end subroutine rt_sc_axis_cosines

  !> Block indices (i1, i2) of the cell at position ih along hdim and iv along updim
  pure subroutine hv(ih, iv, i1, i2)
    integer, intent(in)  :: ih, iv
    integer, intent(out) :: i1, i2
    i1 = eh(1)*ih + ev(1)*iv
    i2 = eh(2)*ih + ev(2)*iv
  end subroutine hv

  !> Reserve the wavefront intensity slots in ps_rt, four in all, one per sweep corner, of which
  !> the ray-by-ray wave uses only the first. Called from rt_sc_activate, and does nothing
  !> unless rt_sc_distributed is set
  subroutine rt_sc_reserve_slots()
    if(.not. rt_sc_distributed) return
    ! A second call leaves the reservation as it is
    if(i_isw > 0) return
    ! The slots are the first columns of ps_rt, which is separate from the hydro state
    nisw  = 4
    i_isw = 1
  end subroutine rt_sc_reserve_slots

  !> Allocate ps_rt(igrid)%w on every active block and zero it, sharing the geometry of ps, which
  !> getbc reads. With bcphys off only %w is exchanged. Blocks created by a regrid get their state
  !> at the next solve, and blocks that go inactive keep it until they are reused
  subroutine rt_sc_alloc_state()
    use mod_global_parameters
    integer :: iigrid, igrid
    if(.not.allocated(ps_rt)) allocate(ps_rt(max_blocks))
    do iigrid=1,igridstail; igrid=igrids(iigrid)
      ! A block whose state was sized for a different nw is allocated afresh
      if(allocated(ps_rt(igrid)%w)) then
        if(size(ps_rt(igrid)%w,ndim+1) /= nw) deallocate(ps_rt(igrid)%w)
      end if
      if(.not.allocated(ps_rt(igrid)%w)) allocate(ps_rt(igrid)%w(ixG^T,1:nw))
      ! The state is zeroed every solve, because getbc with bcphys off never fills the
      ! physical-boundary ghosts, which would otherwise keep stale values
      ps_rt(igrid)%w = 0.0d0
      ps_rt(igrid)%igrid = igrid
      ps_rt(igrid)%level = node(plevel_,igrid)
      ! The geometry points at the hydro block rather than being copied
      {ps_rt(igrid)%ixGmin^D = ixGlo^D; ps_rt(igrid)%ixGmax^D = ixGhi^D;}
      ps_rt(igrid)%x                   => ps(igrid)%x
      ps_rt(igrid)%dvolume             => ps(igrid)%dvolume
      ps_rt(igrid)%is_physical_boundary=> ps(igrid)%is_physical_boundary
    end do
  end subroutine rt_sc_alloc_state

  !> Allocate the seam buffers once, or again when the row or ray count changes, and zero the
  !> buffers written this solve
  subroutine rt_sc_seam_alloc(dnxv, nray)
    integer, intent(in) :: dnxv, nray
    if(allocated(rt_sc_seamL_prev)) then
      if(size(rt_sc_seamL_prev,1)/=dnxv .or. size(rt_sc_seamL_prev,2)/=nray) &
        deallocate(rt_sc_seamL_prev,rt_sc_seamR_prev,rt_sc_seamL_cur,rt_sc_seamR_cur)
    end if
    ! The previous edges start at zero, so the first solve sees no seam inflow
    if(.not.allocated(rt_sc_seamL_prev)) then
      allocate(rt_sc_seamL_prev(dnxv,nray), rt_sc_seamR_prev(dnxv,nray), &
               rt_sc_seamL_cur (dnxv,nray), rt_sc_seamR_cur (dnxv,nray))
      rt_sc_seamL_prev = 0.0d0; rt_sc_seamR_prev = 0.0d0
    end if
    rt_sc_seamL_cur = 0.0d0; rt_sc_seamR_cur = 0.0d0
  end subroutine rt_sc_seam_alloc

  !> Sum the seam edges of this solve over the ranks, where each edge has a single owner, and keep
  !> them for the next solve
  subroutine rt_sc_seam_finalize(dnxv, nray)
    use mod_global_parameters
    integer, intent(in) :: dnxv, nray
    integer :: ierr
    call MPI_ALLREDUCE(MPI_IN_PLACE, rt_sc_seamL_cur, dnxv*nray, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, rt_sc_seamR_cur, dnxv*nray, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr)
    rt_sc_seamL_prev = rt_sc_seamL_cur
    rt_sc_seamR_prev = rt_sc_seamR_cur
  end subroutine rt_sc_seam_finalize

  !> In-place sum allreduce of a 2D array in chunks of whole columns of about 1 MB
  subroutine rt_allreduce_chunk(a, n1, n2, comm)
    use mod_global_parameters
    integer, intent(in)             :: n1, n2, comm
    double precision, intent(inout) :: a(n1,n2)
    integer :: j0, j1, nc, ierr
    ! 131072 doubles make 1 MB, so nc whole columns fit in one chunk
    nc = max(1, 131072/max(n1,1))
    j0 = 1
    do while(j0 <= n2)
      j1 = min(j0+nc-1, n2)
      call MPI_ALLREDUCE(MPI_IN_PLACE, a(1,j0), n1*(j1-j0+1), MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      j0 = j1 + 1
    end do
  end subroutine rt_allreduce_chunk
end module mod_rt_sc_state
