!> Well-balanced atmosphere module.
!>
!> Builds grid-aligned hydrostatic equilibrium tables at the finest AMR
!> resolution from a user-provided fine 1D reference profile. Provides
!> routines for setting the initial condition and boundary conditions
!> from these tables.
!>
!> The module is EoS-independent: it receives p(s) and rho(s) and never
!> calls EoS routines. The multiplicative HSE recurrence uses T = p/rho
!> which is valid for any equation of state.
!>
!> Usage:
!>   1. User integrates HSE to produce fine pa(:), ra(:) arrays
!>   2. call wb_atmos_init(pa, ra, grav, jmax, ds, gzone, grav_dir)
!>   3. call wb_atmos_set_w(ixI^L, ixO^L, w, x) in IC and BC routines
!>   4. call phys_to_conserved(...) after set_w (user's responsibility)
module mod_wb_atmosphere
  use mod_global_parameters
  use mod_comm_lib, only: mpistop

  implicit none
  private

  !> Variable indices for density and pressure — set during wb_atmos_init
  !> to avoid depending on a specific physics module.
  integer :: wb_iw_rho = -1
  integer :: wb_iw_p   = -1

  !> Grid-aligned discrete HSE (built at finest AMR resolution)
  double precision, allocatable, public :: wb_pa(:)   !< equilibrium pressure
  double precision, allocatable, public :: wb_ra(:)   !< equilibrium density
  integer, public :: wb_n_grid = 0                    !< table size

  !> Grid parameters
  integer          :: wb_grav_dir = 1   !< gravity direction (1, 2, or 3)
  double precision :: wb_dx_grid        !< cell spacing at Lmax
  double precision :: wb_x_lo           !< coordinate of first cell centre
  integer          :: wb_n_ghost        !< ghost cells at Lmax resolution
  logical          :: wb_initialised = .false.

  public :: wb_atmos_init
  public :: wb_atmos_set_w
  public :: wb_atmos_get_eq

contains

  !> Build grid-aligned HSE tables from a fine 1D reference profile.
  !>
  !> The fine table (pa, ra) spans [xmin - gzone, xmax + gzone] along the
  !> gravity direction, evenly spaced at ds. The gravity array grav(jmax)
  !> gives the gravitational acceleration at each fine table point.
  !>
  !> The grid-aligned table is built at the finest AMR resolution using the
  !> multiplicative trapezoidal recurrence, which matches the discrete HSE
  !> that the WB reconstruction preserves.
  subroutine wb_atmos_init(pa, ra, grav, jmax, ds, gzone, grav_dir, irho, ip)
    double precision, intent(in) :: pa(:)    !< fine pressure table
    double precision, intent(in) :: ra(:)    !< fine density table
    double precision, intent(in) :: grav(:)  !< fine gravity table
    integer, intent(in) :: jmax             !< fine table size
    double precision, intent(in) :: ds       !< fine table spacing
    double precision, intent(in) :: gzone    !< ghost zone width (code units)
    integer, intent(in) :: grav_dir          !< gravity direction (1, 2, or 3)
    integer, intent(in) :: irho              !< w-array index for density
    integer, intent(in) :: ip                !< w-array index for pressure

    double precision :: x_g, T_lo, T_hi, g_lo, g_hi, res
    double precision :: alpha_prev, alpha_curr
    double precision :: x_min, x_max, p_dev, max_dev
    integer :: ig, na, n_cells, level_factor

    wb_grav_dir = grav_dir
    wb_iw_rho = irho
    wb_iw_p   = ip

    ! Grid parameters at finest AMR level
    level_factor = 2**(refine_max_level - 1)
    select case(grav_dir)
    case(1)
      n_cells = domain_nx1 * level_factor
      x_min = xprobmin1; x_max = xprobmax1
    {^IFTWOD
    case(2)
      n_cells = domain_nx2 * level_factor
      x_min = xprobmin2; x_max = xprobmax2
    }
    {^IFTHREED
    case(2)
      n_cells = domain_nx2 * level_factor
      x_min = xprobmin2; x_max = xprobmax2
    case(3)
      n_cells = domain_nx3 * level_factor
      x_min = xprobmin3; x_max = xprobmax3
    }
    case default
      call mpistop('wb_atmos_init: invalid grav_dir')
    end select

    wb_n_ghost = nghostcells * level_factor
    wb_n_grid  = n_cells + 2 * wb_n_ghost
    wb_dx_grid = (x_max - x_min) / dble(n_cells)
    wb_x_lo    = x_min - (dble(wb_n_ghost) - 0.5d0) * wb_dx_grid

    allocate(wb_pa(wb_n_grid), wb_ra(wb_n_grid))

    ! Interpolate T = p/rho and g from fine table at each grid cell centre,
    ! then build grid-aligned table via multiplicative HSE recurrence.
    !
    ! Seed the first cell: interpolate p and T from fine table.
    x_g = wb_x_lo
    call interp_fine(x_g, x_min, gzone, ds, jmax, pa, ra, grav, T_lo, g_lo)
    ! Use same indexing for seed pressure
    na = floor((x_g - (x_min - gzone)) / ds + 0.5d0)
    na = max(1, min(na, jmax - 1))
    res = (x_g - (x_min - gzone)) / ds - (dble(na) - 0.5d0)
    res = max(0.0d0, min(1.0d0, res))
    wb_pa(1) = pa(na) + res * (pa(na + 1) - pa(na))
    wb_ra(1) = wb_pa(1) / T_lo

    ! Forward recurrence
    do ig = 2, wb_n_grid
      x_g = wb_x_lo + dble(ig - 1) * wb_dx_grid
      call interp_fine(x_g, x_min, gzone, ds, jmax, pa, ra, grav, T_hi, g_hi)

      alpha_prev = 0.5d0 * wb_dx_grid * g_lo / T_lo
      alpha_curr = 0.5d0 * wb_dx_grid * g_hi / T_hi

      wb_pa(ig) = wb_pa(ig - 1) * (1.0d0 + alpha_prev) &
                  / (1.0d0 - alpha_curr)
      wb_ra(ig) = wb_pa(ig) / T_hi

      T_lo = T_hi
      g_lo = g_hi
    end do

    ! Diagnostic: max deviation from fine table
    max_dev = 0.0d0
    do ig = 1, wb_n_grid
      x_g = wb_x_lo + dble(ig - 1) * wb_dx_grid
      na = max(1, min(nint((x_g - (x_min - gzone))/ds) + 1, jmax))
      if (pa(na) > 0.0d0) then
        p_dev = dabs(wb_pa(ig) - pa(na)) / pa(na)
        max_dev = max(max_dev, p_dev)
      end if
    end do

    if (mype == 0) then
      write(*,'(A,I7,A,ES10.3,A,I1)') &
        ' WB atmosphere: ', wb_n_grid, ' cells, max p deviation = ', &
        max_dev, ', grav_dir = ', wb_grav_dir
    end if

    wb_initialised = .true.

  end subroutine wb_atmos_init

  !> Set w to the equilibrium state (rho, v=0, p) at each cell position.
  !>
  !> Point-value sample of the trap-rule-built grid table at cell centres.
  !> Consistent with the legacy_trap WB transform which expects point-value
  !> IC matching its own trap-rule recurrence at L_max.
  !>
  !> On exit: w(rho_) = rho_eq, w(mom(:)) = 0, w(p_) = p_eq.
  !> The caller must call phys_to_conserved afterwards.
  subroutine wb_atmos_set_w(ixI^L, ixO^L, w, x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision :: s
    integer :: ig
    integer :: ix^D

    if (.not. wb_initialised) &
      call mpistop('wb_atmos_set_w: call wb_atmos_init first')

    ! Zero all conserved/primitive variables first
    w(ixO^S, 1:nw) = 0.0d0

    ! Set density and pressure from the equilibrium table
    {do ix^DB = ixOmin^DB, ixOmax^DB\}
      s = x(ix^D, wb_grav_dir)
      ig = nint((s - wb_x_lo) / wb_dx_grid) + 1
      ig = max(1, min(ig, wb_n_grid))

      w(ix^D, wb_iw_rho) = wb_ra(ig)
      w(ix^D, wb_iw_p)   = wb_pa(ig)
    {end do\}

  end subroutine wb_atmos_set_w

  !> Look up equilibrium (rho, p) at an arbitrary coordinate along the
  !> gravity direction. Uses nearest-cell lookup from the grid-aligned table.
  subroutine wb_atmos_get_eq(x_grav, rho_eq, p_eq)
    double precision, intent(in)  :: x_grav
    double precision, intent(out) :: rho_eq, p_eq

    integer :: ig

    if (.not. wb_initialised) &
      call mpistop('wb_atmos_get_eq: call wb_atmos_init first')

    ig = nint((x_grav - wb_x_lo) / wb_dx_grid) + 1
    ig = max(1, min(ig, wb_n_grid))

    rho_eq = wb_ra(ig)
    p_eq   = wb_pa(ig)

  end subroutine wb_atmos_get_eq

  ! ── Private helpers ──────────────────────────────────────────────

  !> Interpolate T = p/rho and g from the fine table at position x_g.
  subroutine interp_fine(x_g, x_min, gzone, ds, jmax, pa, ra, grav, T_out, g_out)
    double precision, intent(in)  :: x_g, x_min, gzone, ds
    integer, intent(in)           :: jmax
    double precision, intent(in)  :: pa(:), ra(:), grav(:)
    double precision, intent(out) :: T_out, g_out

    integer :: na
    double precision :: frac

    ! Cell-centred indexing matching inithdstatic convention:
    ! na = floor((x + gzone)/ds + 0.5), res = remainder within cell
    na = floor((x_g - (x_min - gzone)) / ds + 0.5d0)
    na = max(1, min(na, jmax - 1))
    frac = (x_g - (x_min - gzone)) / ds - (dble(na) - 0.5d0)
    frac = frac / 1.0d0  ! already normalised to ds
    frac = max(0.0d0, min(1.0d0, frac))

    ! Interpolate T = p/rho ratio directly (NOT p and rho separately).
    T_out = pa(na)/ra(na) + frac * (pa(na+1)/ra(na+1) - pa(na)/ra(na))

    g_out = grav(na) + frac * (grav(na + 1) - grav(na))

  end subroutine interp_fine

end module mod_wb_atmosphere
