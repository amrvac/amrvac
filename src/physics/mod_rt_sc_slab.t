!=============================================================================
!> Public face of the short-characteristics (SC) radiative transfer and its coupling to the
!> (M)HD
!>
!> The rest of AMRVAC (amrvac.t, mod_errest and mod_radiative_cooling) and the user files reach
!> the RT through this module alone. It re-exports mod_rt_sc_state and mod_rt_sc_solvers and adds
!> the routines that couple the transfer to the gas. rt_sc_add_source applies the heating every
!> stage, limited so that one stage changes the gas energy by at most rt_sc_cfl of itself.
!> rt_sc_thin_weight divides the radiative losses between SC and the thin cooling.
!> rt_sc_pin_refine holds the bands on rt_sc_amr_level, and rt_sc_emergent_tau1 measures the flux
!> on the tau = 1 surface for a flux controller, and rt_sc_selftest checks the sweep on an analytic
!> slab
!>
!> The modules form a chain without cycles. mod_rt_sc_opacity and mod_rt_sc_kernel are leaves,
!> mod_rt_sc_state uses both, mod_rt_sc_solvers uses the state, and this module uses the solvers
!>
!> Jack Jenkins (10.2026); More or less finished in v1 state for amrvac4.0
!=============================================================================
module mod_rt_sc_slab
  use mod_rt_sc_state
  use mod_rt_sc_solvers
  implicit none

  public :: rt_sc_activate
  public :: rt_sc_thin_weight
  public :: rt_sc_add_source
  public :: rt_sc_selftest
  public :: rt_sc_pin_refine
  public :: rt_sc_emergent_tau1

contains

  !> Switch the RT on. The physics module calls this from its initialisation when
  !> hd_radiative_transfer or mhd_radiative_transfer is set. It reads &rt_list, reserves the
  !> intensity slots of the wavefront and registers the output variables Qrad and Frad, with Iout
  !> as well for the replicated solver and the validation mode, which fill it. rt_sc_setup follows
  !> once the units and the domain are known
  subroutine rt_sc_activate()
    use mod_global_parameters, only: par_files
    use mod_variables, only: var_set_extravar
    call rt_sc_params_read(par_files)
    rt_sc_active = .true.
    call rt_sc_reserve_slots()
    i_qrad = var_set_extravar('Qrad', 'Qrad')
    i_frad = var_set_extravar('Frad', 'Frad')
    if(.not.rt_sc_distributed .or. rt_sc_valpipe) i_iout = var_set_extravar('Iout', 'Iout')
  end subroutine rt_sc_activate

  !> Hold every block overlapping an RT band on rt_sc_amr_level. Called from forcedrefine_grid
  !> after usr_refine_grid, so it overrides both the estimator and the user
  subroutine rt_sc_pin_refine(igrid, level, my_refine, my_coarsen)
    use mod_global_parameters
    integer, intent(in)    :: igrid, level
    integer, intent(inout) :: my_refine, my_coarsen
    integer :: ib
    double precision :: blo, bhi
    if(.not. rt_sc_active) return
    ! Only the extent along updim matters, since a band spans the full width of the domain
    blo = rnode(rpxmin0_+updim,igrid)
    bhi = rnode(rpxmax0_+updim,igrid)
    do ib = 1, rt_sc_nband
      if(bhi > band_lo_c(ib) .and. blo < band_hi_c(ib)) then
        if(level < rt_sc_amr_level) then
          my_refine = 1;  my_coarsen = 0
        else if(level == rt_sc_amr_level) then
          my_refine = -1; my_coarsen = -1
        else
          my_refine = 0;  my_coarsen = 1
        end if
        return
      end if
    end do
  end subroutine rt_sc_pin_refine

  !> Self-test of the real sweep on an analytic slab. J must approach B
  !> where the slab is optically thick, stay finite in the tenuous top, where the gray FLD
  !> coefficient D = c lambda/(kappa rho) would diverge, and stay bounded by the source everywhere.
  !> The slab is horizontally infinite, so the test sets thermal side inflows, which are exact for
  !> it, a thermal bottom and an open top, and restores the run's boundaries afterwards
  subroutine rt_sc_selftest()
    use mod_global_parameters
    {^IFTWOD
    integer :: j, k
    double precision :: tau, taumax, B0, ap, jb_err_deep, jmax, jmin, chimax, chimin
    logical :: run_bc(5)
    run_bc = [bc_h_periodic, bc_hlo_thermal, bc_hhi_thermal, bc_lo_thermal, bc_hi_thermal]
    bc_h_periodic = .false.; bc_hlo_thermal = .true.; bc_hhi_thermal = .true.
    bc_lo_thermal = .true.;  bc_hi_thermal = .false.
    ! The source is uniform, S = B = B0, and chi grows geometrically downwards over six decades
    nxg = 4; nyg = 200
    dxg = 1.0d0; dyg = 1.0d0
    call rt_grid_alloc()
    ! The quadrature is built here when no solve has built it yet
    if(.not.allocated(mux)) call rt_sc_init()
    B0 = 1.0d0
    taumax = 1.0d3
    do j = 1, nyg
      ! ap runs from one in the bottom row to zero in the top row
      ap = dble(nyg-j)/dble(nyg-1)
      ! chi_j = (taumax/nyg) 10^(6 ap - 3), which falls from 1e3 to 1e-3 times taumax/nyg
      chi_g(:,j) = (taumax/dble(nyg)) * 10.0d0**(6.0d0*ap - 3.0d0)
      b_g(:,j)   = B0
    end do
    call sc_pad_grid()
    ! J is as always the weighted sum of the intensities of all rays
    jm_g = 0.0d0
    do k = 1, rt_sc_nray
      call sc_sweep_2d(mux(k), muy(k))
      ! J = sum_k w_k I_k
      jm_g(:,:) = jm_g(:,:) + wmu(k)*ii_g(:,:)
    end do
    jmax = maxval(jm_g); jmin = minval(jm_g)
    chimax = maxval(chi_g); chimin = minval(chi_g)
    ! The optical depth is accumulated down the first column. Below tau ~ 10 the relative
    ! departure of J from B measures how well the deep slab thermalises
    tau = 0.0d0; jb_err_deep = 0.0d0
    do j = nyg, 1, -1
      tau = tau + chi_g(1,j)*dyg
      ! err = max abs(J - B)/B over the cells with tau > 10
      if(tau > 10.0d0) jb_err_deep = max(jb_err_deep, dabs(jm_g(1,j)-b_g(1,j))/b_g(1,j))
    end do
    ! The test passes when J is finite and non-negative everywhere and within 5% of B deep down
    if(mype==0) then
      print *, "=== rt_sc_selftest (uniform-source slab, B0=1) ==="
      print *, " grid nxg,nyg =", nxg, nyg, "  rays =", rt_sc_nray
      print *, " chi per-cell range :", chimin, chimax
      print *, " J range            :", jmin, jmax, "  (should be finite, in [0, ~B0])"
      print *, " max |J-B|/B for tau>10 (deep thermalisation, want <~1e-2):", jb_err_deep
      print *, " J at top cell (thin), col1                   :", jm_g(1,nyg)
      print *, " J at bottom cell (deep), col1               :", jm_g(1,1)
      if(jmax < bigdouble .and. jmin >= -1.0d-10 .and. jb_err_deep < 5.0d-2) then
        print *, " result: pass (finite, non-negative, thermalises deep)"
      else
        print *, " result: fail"
      end if
    end if
    bc_h_periodic = run_bc(1); bc_hlo_thermal = run_bc(2); bc_hhi_thermal = run_bc(3)
    bc_lo_thermal = run_bc(4); bc_hi_thermal = run_bc(5)
    }
  end subroutine rt_sc_selftest

  !> Horizontal mean of Fy on the tau = rt_sc_tau_emergent surface, the photosphere, into
  !> rt_sc_Ftau1. The flux controller needs this value, since Frad is not conserved above tau = 1
  !> and the flux at the band top reads larger. Only the first band is searched, since it is the
  !> one that holds the photosphere, and tau comes from the band top as in the KBA heat blend
  subroutine rt_sc_emergent_tau1()
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    {^IFTWOD
    double precision :: taub(ixG^T), tau_up, fy_up, ww, gbuf(2), gout(2)
    double precision, allocatable :: tautop(:,:)
    integer :: igrid, iblk, ih, iv, i1, i2, ig, ierr
    if(.not.rt_sc_active) return
    if(.not.rt_sc_ready) call mpistop("rt_sc_emergent_tau1: setup not called")
    if(i_frad < 0)       call mpistop("rt_sc_emergent_tau1: i_frad slot not registered")
    if(.not.associated(rt_sc_chiB)) call mpistop("rt_sc_emergent_tau1: opacity callback not set")
    call rt_band_blocks(1)
    allocate(tautop(rg_nh,wv_nbxv))
    call rt_tau_tops(tautop)
    ! Each column walks down from the top face of its block, and Fy is interpolated linearly in
    ! tau where the column crosses the surface. At the top face the flux takes its value in the
    ! top cell. gbuf sums Fy at the crossings and counts the columns that found one
    gbuf = 0.0d0
    do iblk = 1, wv_nblk
      igrid = wv_ig(iblk)
      call rt_block_tau(iblk, tautop, taub)
      do ih = iMlo_h, iMhi_h
        ig = (wv_bi(iblk)-1)*rg_bh + (ih - iMlo_h + 1)
        tau_up = tautop(ig,wv_bj(iblk))
        ! When the top face already lies below the surface, the crossing is in a block above
        if(tau_up >= rt_sc_tau_emergent) cycle
        call hv(ih, iMhi_v, i1, i2)
        fy_up = ps(igrid)%w(i1,i2,i_frad)
        do iv = iMhi_v, iMlo_v, -1
          call hv(ih, iv, i1, i2)
          if(taub(i1,i2) >= rt_sc_tau_emergent) then
            ! Fy(tau_e) = Fy_up + ww (Fy - Fy_up),  ww = (tau_e - tau_up)/(tau_c - tau_up)
            ww = (rt_sc_tau_emergent-tau_up)/max(taub(i1,i2)-tau_up,1.0d-30)
            gbuf(1) = gbuf(1) + fy_up + ww*(ps(igrid)%w(i1,i2,i_frad)-fy_up)
            gbuf(2) = gbuf(2) + 1.0d0
            exit
          end if
          tau_up = taub(i1,i2)
          fy_up  = ps(igrid)%w(i1,i2,i_frad)
        end do
      end do
    end do
    ! Ftau1 = sum of Fy(tau_e) over the crossings / number of crossings
    call MPI_ALLREDUCE(gbuf, gout, 2, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr)
    rt_sc_Ftau1 = gout(1)/max(gout(2),1.0d0)
    deallocate(tautop)
    }
  end subroutine rt_sc_emergent_tau1

  !> Unsplit heating source of one block and stage, e += qdt Q, with Q read from i_qrad, which is
  !> zero outside the bands. The heating is limited so that one stage changes the gas energy by
  !> at most rt_sc_cfl of itself
  subroutine rt_sc_add_source(qdt, ixI^L, ixO^L, wCT, w, x, thin_share, e_internal)
    use mod_global_parameters
    use mod_eos, only: eos
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    !> True when the thin cooling is also on. The SC heating then takes only the share
    !> 1 - w_thin, so that no loss is counted twice
    logical, intent(in)             :: thin_share
    !> True when the energy variable is the internal energy, as with mhd_internal_e
    logical, intent(in)             :: e_internal
    {^IFTWOD
    integer :: ix1, ix2
    double precision :: egas, ke, em, qlim, dn, qc, Te(ixI^S), pth(ixI^S)
    if(.not.rt_sc_active .or. i_qrad < 0) return
    ! The thin-cooling share needs the temperature and the gas pressure, both from the EoS
    if(thin_share) then
      call eos%get_Te(wCT, x, ixI^L, ixO^L, Te)
      call eos%get_thermal_pressure(wCT, x, ixI^L, ixO^L, pth)
    end if
    do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        ! Cells without heating, which include every cell outside the bands, are skipped
        qc = wCT(ix1,ix2,i_qrad)
        if(qc == 0.0d0) cycle
        ! The gas energy is the total less the kinetic and magnetic parts, taken from wCT like Q,
        ! or the energy variable itself when that is already the internal energy
        ! e_gas = e - rho v^2/2 - B^2/2
        if(e_internal) then
          egas = wCT(ix1,ix2,iw_e)
        else
          ke = 0.5d0*sum(wCT(ix1,ix2,iw_mom(:))**2)/wCT(ix1,ix2,iw_rho)
          em = 0.0d0
          if(allocated(iw_mag)) em = 0.5d0*sum(wCT(ix1,ix2,iw_mag(:))**2)
          egas = wCT(ix1,ix2,iw_e) - ke - em
        end if
        ! dn is the largest rate that changes egas by rt_sc_cfl of itself over qdt. The limiter
        ! leaves smaller rates untouched and caps larger ones at dn, keeping their sign
        ! dn = cfl e_gas/qdt,  q_lim = q dn/max(dn, abs(q)),  e = e + qdt w q_lim
        dn   = rt_sc_cfl*max(egas,smalldouble)/max(qdt,smalldouble)
        qlim = qc*dn/max(dn,dabs(qc))
        if(thin_share) qlim = qlim*(1.0d0 - rt_sc_thin_weight(Te(ix1,ix2), x(ix1,ix2,updim), pth(ix1,ix2)))
        w(ix1,ix2,iw_e) = w(ix1,ix2,iw_e) + qdt*qlim
      end do
    end do
    }
  end subroutine rt_sc_add_source

  !> Share of the radiative losses taken by the optically thin treatment, with SC taking the rest,
  !> 1 - w_thin. The thin cooling and the user SC source both call this, so they share one
  !> definition. Inside a band w_thin is the product of a temperature gate and a pressure gate, so
  !> a cell goes to SC when it is cool or dense, as in the opacity switch. The pressure gate
  !> max(0, 1 - (p/p_TR)^2) of MURaM's Qthin is exactly zero in the photosphere. Outside every
  !> band, or without RT, w_thin = 1, because there is no SC there to take over. The step to one
  !> at the band edge marks where SC stops, and the gas there is coronal, with the temperature
  !> gate already close to one
  pure elemental function rt_sc_thin_weight(Te_code, y_code, p_code) result(w_thin)
    use mod_global_parameters, only: unit_temperature, unit_pressure
    !> Temperature in code units
    double precision, intent(in) :: Te_code
    !> Height along updim in code units
    double precision, intent(in) :: y_code
    !> Thermal pressure in code units
    double precision, intent(in) :: p_code
    double precision :: w_thin, pr
    integer :: ib

    w_thin = 1.0d0
    if(.not. rt_sc_active) return
    do ib = 1, rt_sc_nband
      ! The band edges are those after the block snap, which is where SC acts
      if(y_code >= band_lo_c(ib) .and. y_code <= band_hi_c(ib)) then
        ! The gate parameters are in K and cgs, so the comparison is made in physical units
        ! w_thin = (1 + tanh((T - T0)/dT))/2 max(0, 1 - (p/p_TR)^2)
        pr = p_code*unit_pressure/rt_sc_p_TR
        w_thin = 0.5d0*(1.0d0 + tanh((Te_code*unit_temperature - rt_sc_blend_T0)/rt_sc_blend_dT)) &
               * max(0.0d0, 1.0d0 - pr*pr)
        return
      end if
    end do
  end function rt_sc_thin_weight

end module mod_rt_sc_slab
