!=============================================================================
!> Entropy-method LTE EoS: every query is ONE bicubic Hermite evaluation.
!>
!> Each quantity is stored as four tables -- the value plus its three
!> derivatives (d/dx, d/dy, d2/dxdy) at every node -- which fully determine
!> the bicubic Hermite polynomial in each cell: exact at nodes, O(h^4) inside,
!> no closure. No Newton, bisection, or iteration on the hot path; the
!> iterative work is done offline against the analytic Saha solver (see
!> entropy128/generate_entropy_tables.py).
!>
!> Forward tables, axes (log10 nH, log10 eint/nH):
!>     Tfwd -> T    pfwd -> p    neOnH -> ne/nH
!> Inverse tables:
!>     eintP (log nH, log p/nH) -> log10(eint/nH)
!>     eintT (log nH, log T)    -> log10(eint/nH)
!>     g1p   (log nH, log p/nH) -> Gamma_1
!>
!> Gamma_1 is tabulated directly (g1p) rather than recovered from second
!> derivatives, so the lookup stays smooth; values are Maxwell-consistent with
!> (p, T) by construction.
!=============================================================================
module mod_eos_lte_entropy
    use mod_global_parameters
    use mod_comm_lib,      only: mpistop
    use mod_eos_container, only: eos_table_container, eos
    use mod_eos_lte_tables, only: load_tables_LTE, precompute_FI_bypass_constants, &
         ensure_axis_nodes, shift_axis_to_code, shift_axis_to_code_T, &
         eos_build_guards, eos_validate_table
    use mod_eos_interp, only: precompute_step_inv
    implicit none
    private

    double precision, parameter :: LN10 = 2.302585092994046d0

    !> Unit conversions fixed at finalise: log10 of the code-to-CGS factor of
    !> eint/nH (= p/nH), and 1/unit_pressure
    double precision :: ent_shift_pe = 0.0d0, ent_inv_up = 1.0d0

    !> Value and derivative tables of each quantity interleaved as
    !> (f, fx, fy, fxy) per node, built at finalise; axes stay in eos%<quantity>
    double precision, allocatable :: il_neOnH(:,:,:), il_Tfwd(:,:,:), il_pfwd(:,:,:)
    double precision, allocatable :: il_eintP(:,:,:), il_g1p(:,:,:), il_eintT(:,:,:)

    public :: load_entropy_LTE, finalise_entropy_LTE

    !> Forward (log nH, log eint/nH) -> thermodynamic state
    public :: entropy_T_from_nH_eint
    public :: entropy_p_nH_from_eint, entropy_p_from_eint
    public :: entropy_y_from_nH_eint
    public :: entropy_T_and_y_from_nH_eint
    !> Inverse (log nH, log p/nH) -> eint/p ratio, Gamma_1
    public :: entropy_eint_from_p_bisect
    public :: entropy_eint_from_nH_p
    public :: entropy_gamma1_from_nH_p
    !> Inverse (log nH, log T) -> log10(eint/nH)
    public :: entropy_eint_from_nH_T

contains

    !> Entropy-method load: six quantities, each a value table plus its three
    !> derivative tables, on the forward (nH, eint/nH), inverse-p (nH, p/nH) and
    !> inverse-T (nH, T) grids. Every runtime query is one bicubic-Hermite
    !> evaluation (see generate_entropy_tables.py).
    subroutine load_entropy_LTE()
        integer :: iq, id
        character(len=8), parameter :: quantity(6) = &
            [character(8):: 'neOnH', 'Tfwd', 'pfwd', 'eintP', 'g1p', 'eintT']
        character(len=3), parameter :: deriv(4) = &
            [character(3):: '', '_x', '_y', '_xy']

        if (mype == 0) write(*,*) &
            'EoS method: entropy (bicubic Hermite, 4-derivative corners, no closure)'
        do iq = 1, size(quantity)
            do id = 1, size(deriv)
                call load_tables_LTE(trim(quantity(iq))//trim(deriv(id)))
            end do
        end do
    end subroutine load_entropy_LTE

    !> Entropy-method finalise: shift each loaded table's axes from CGS to code
    !> units and run the standard prepare. No aux builds -- every runtime quantity
    !> is a single Hermite lookup. neOnH/Tfwd/pfwd/eintP/g1p are on (log nH,
    !> log eint/nH) (axis2_is_T=.false.); eintT is on (log nH, log T)
    !> (axis2_is_T=.true.). Each quantity ships value + three derivative tables.
    subroutine finalise_entropy_LTE()
        call entropy_shift_prepare(eos%neOnH,    .false.)
        call entropy_shift_prepare(eos%neOnH_x,  .false.)
        call entropy_shift_prepare(eos%neOnH_y,  .false.)
        call entropy_shift_prepare(eos%neOnH_xy, .false.)
        call entropy_shift_prepare(eos%Tfwd,     .false.)
        call entropy_shift_prepare(eos%Tfwd_x,   .false.)
        call entropy_shift_prepare(eos%Tfwd_y,   .false.)
        call entropy_shift_prepare(eos%Tfwd_xy,  .false.)
        call entropy_shift_prepare(eos%pfwd,     .false.)
        call entropy_shift_prepare(eos%pfwd_x,   .false.)
        call entropy_shift_prepare(eos%pfwd_y,   .false.)
        call entropy_shift_prepare(eos%pfwd_xy,  .false.)
        call entropy_shift_prepare(eos%eintP,    .false.)
        call entropy_shift_prepare(eos%eintP_x,  .false.)
        call entropy_shift_prepare(eos%eintP_y,  .false.)
        call entropy_shift_prepare(eos%eintP_xy, .false.)
        call entropy_shift_prepare(eos%g1p,      .false.)
        call entropy_shift_prepare(eos%g1p_x,    .false.)
        call entropy_shift_prepare(eos%g1p_y,    .false.)
        call entropy_shift_prepare(eos%g1p_xy,   .false.)
        call entropy_shift_prepare(eos%eintT,    .true.)
        call entropy_shift_prepare(eos%eintT_x,  .true.)
        call entropy_shift_prepare(eos%eintT_y,  .true.)
        call entropy_shift_prepare(eos%eintT_xy, .true.)
        call precompute_FI_bypass_constants()
        ent_shift_pe = dlog10(unit_pressure / unit_numberdensity)
        ent_inv_up   = 1.0d0 / unit_pressure
        call interleave(eos%neOnH, eos%neOnH_x, eos%neOnH_y, eos%neOnH_xy, il_neOnH)
        call interleave(eos%Tfwd,  eos%Tfwd_x,  eos%Tfwd_y,  eos%Tfwd_xy,  il_Tfwd)
        call interleave(eos%pfwd,  eos%pfwd_x,  eos%pfwd_y,  eos%pfwd_xy,  il_pfwd)
        call interleave(eos%eintP, eos%eintP_x, eos%eintP_y, eos%eintP_xy, il_eintP)
        call interleave(eos%g1p,   eos%g1p_x,   eos%g1p_y,   eos%g1p_xy,   il_g1p)
        call interleave(eos%eintT, eos%eintT_x, eos%eintT_y, eos%eintT_xy, il_eintT)
        !> The entropy method never loads eos%T, but several callers read
        !> eos%T%var{1,2}_{min,max} (the cold-cell eint floor, FI-fallback
        !> checks). Source them from eos%Tfwd, on the same forward grid.
        eos%T%var2_min = eos%Tfwd%var2_min
        eos%eint_nH_min = 10.0d0**eos%T%var2_min
        eos%T%var2_max = eos%Tfwd%var2_max
        eos%T%var1_min = eos%Tfwd%var1_min
        eos%T%var1_max = eos%Tfwd%var1_max
        if (mype == 0) write(*,'(A)') &
            ' EoS entropy: forward + inverse tables loaded; no aux builds.'
    end subroutine finalise_entropy_LTE

    !> Index location: uniform and non-uniform variants.
    pure subroutine locate_idx_uniform(n, vmin, vmax, val, ix, t, h)
        integer,          intent(in)  :: n
        double precision, intent(in)  :: vmin, vmax, val
        integer,          intent(out) :: ix
        double precision, intent(out) :: t, h
        h = (vmax - vmin) / dble(n - 1)
        if (val <= vmin) then
            ix = 1;     t = 0.0d0;  return
        end if
        if (val >= vmax) then
            ix = n - 1; t = 1.0d0;  return
        end if
        ix = 1 + int((val - vmin) / h)
        if (ix < 1)     ix = 1
        if (ix > n - 1) ix = n - 1
        t  = (val - (vmin + (ix - 1) * h)) / h
        if (t < 0.0d0) t = 0.0d0
        if (t > 1.0d0) t = 1.0d0
    end subroutine locate_idx_uniform

    pure subroutine locate_idx_nodes(n, nodes, val, ix, t, h)
        integer,          intent(in)  :: n
        double precision, intent(in)  :: nodes(n), val
        integer,          intent(out) :: ix
        double precision, intent(out) :: t, h
        integer :: lo, mid, hi
        if (val <= nodes(1)) then
            ix = 1; t = 0.0d0; h = nodes(2) - nodes(1); return
        end if
        if (val >= nodes(n)) then
            ix = n - 1; t = 1.0d0; h = nodes(n) - nodes(n - 1); return
        end if
        ! Binary search for cell containing val.
        lo = 1; hi = n
        do
            if (hi - lo <= 1) exit
            mid = (lo + hi) / 2
            if (nodes(mid) <= val) then
                lo = mid
            else
                hi = mid
            end if
        end do
        ix = lo
        h = nodes(ix + 1) - nodes(ix)
        t = (val - nodes(ix)) / h
        if (t < 0.0d0) t = 0.0d0
        if (t > 1.0d0) t = 1.0d0
    end subroutine locate_idx_nodes

    !> The bracket of locate_idx_nodes (largest node index lo <= n-1 with
    !> nodes(lo) <= val) in O(1) from the guard buckets: guard(k) is the first node
    !> at or above the start of bucket k, corrected forward to the first node
    !> >= val and back across a bucket edge that rounding put above val
    pure subroutine locate_idx_guard(n, nodes, guard, M, scale, val, ix, t, h)
        integer,          intent(in)  :: n, M
        double precision, intent(in)  :: nodes(n), scale, val
        integer,          intent(in)  :: guard(M)
        integer,          intent(out) :: ix
        double precision, intent(out) :: t, h
        integer :: k, j
        if (val <= nodes(1)) then
            ix = 1; t = 0.0d0; h = nodes(2) - nodes(1); return
        end if
        if (val >= nodes(n)) then
            ix = n - 1; t = 1.0d0; h = nodes(n) - nodes(n - 1); return
        end if
        k = 1 + min(M-1, max(0, int((val - nodes(1)) * scale)))
        j = guard(k)
        do while (j < n .and. nodes(j) < val)
            j = j + 1
        end do
        do while (j > 1 .and. nodes(j-1) >= val)
            j = j - 1
        end do
        if (nodes(j) > val) j = j - 1
        ix = min(max(j, 1), n - 1)
        h = nodes(ix + 1) - nodes(ix)
        t = (val - nodes(ix)) / h
        if (t < 0.0d0) t = 0.0d0
        if (t > 1.0d0) t = 1.0d0
    end subroutine locate_idx_guard

    pure subroutine locate_idx_axis1(tab, val, ix, t, h)
        type(eos_table_container), intent(in) :: tab
        double precision, intent(in)  :: val
        integer,          intent(out) :: ix
        double precision, intent(out) :: t, h
        if (tab%is_uniform) then
            call locate_idx_uniform(tab%dim1, tab%var1_min, tab%var1_max, val, ix, t, h)
        else if (tab%guard_M_1 > 0) then
            call locate_idx_guard(tab%dim1, tab%var1_nodes, tab%guard_1, tab%guard_M_1, &
                                  tab%guard_scale_1, val, ix, t, h)
        else
            call locate_idx_nodes(tab%dim1, tab%var1_nodes, val, ix, t, h)
        end if
    end subroutine locate_idx_axis1

    pure subroutine locate_idx_axis2(tab, val, ix, t, h)
        type(eos_table_container), intent(in) :: tab
        double precision, intent(in)  :: val
        integer,          intent(out) :: ix
        double precision, intent(out) :: t, h
        if (tab%is_uniform) then
            call locate_idx_uniform(tab%dim2, tab%var2_min, tab%var2_max, val, ix, t, h)
        else if (tab%guard_M_2 > 0) then
            call locate_idx_guard(tab%dim2, tab%var2_nodes, tab%guard_2, tab%guard_M_2, &
                                  tab%guard_scale_2, val, ix, t, h)
        else
            call locate_idx_nodes(tab%dim2, tab%var2_nodes, val, ix, t, h)
        end if
    end subroutine locate_idx_axis2

    !> 1D cubic Hermite basis on t in [0, 1] (value weights).
    !> H0 = 2t^3 - 3t^2 + 1   value at t=0      H2 = -2t^3 + 3t^2   value at t=1
    !> H1 = t^3 - 2t^2 + t    slope at t=0      H3 = t^3 - t^2      slope at t=1
    pure subroutine cubic_basis(t, H0, H1, H2, H3)
        double precision, intent(in)  :: t
        double precision, intent(out) :: H0, H1, H2, H3
        double precision :: t2, t3
        t2 = t * t
        t3 = t2 * t
        H0 = 2.0d0*t3 - 3.0d0*t2 + 1.0d0
        H1 = t3 - 2.0d0*t2 + t
        H2 = -2.0d0*t3 + 3.0d0*t2
        H3 = t3 - t2
    end subroutine cubic_basis

    !> Bicubic Hermite contraction at a located cell. c(1:4,i,j) holds
    !> (f, df/dx, df/dy, d2f/dxdy) at node (i,j); the derivatives are physical-axis
    !> derivatives, scaled by the local cell widths to unit-interval coordinates.
    !> 4 stored values x 4 corners = 16 conditions, which fully determine the
    !> degree-(3,3) polynomial in the cell
    pure double precision function contract_il(c, ix, iy, dx, dy, &
                                               H0x, H1x, H2x, H3x, &
                                               H0y, H1y, H2y, H3y) result(val)
        double precision, intent(in) :: c(:,:,:)
        integer, intent(in) :: ix, iy
        double precision, intent(in) :: dx, dy
        double precision, intent(in) :: H0x, H1x, H2x, H3x
        double precision, intent(in) :: H0y, H1y, H2y, H3y
        double precision :: f00, f10, f01, f11
        double precision :: fx00, fx10, fx01, fx11
        double precision :: fy00, fy10, fy01, fy11
        double precision :: fxy00, fxy10, fxy01, fxy11
        f00 = c(1,ix,  iy  );    f10 = c(1,ix+1,iy  )
        f01 = c(1,ix,  iy+1);    f11 = c(1,ix+1,iy+1)
        fx00 = c(2,ix,  iy  ) * dx
        fx10 = c(2,ix+1,iy  ) * dx
        fx01 = c(2,ix,  iy+1) * dx
        fx11 = c(2,ix+1,iy+1) * dx
        fy00 = c(3,ix,  iy  ) * dy
        fy10 = c(3,ix+1,iy  ) * dy
        fy01 = c(3,ix,  iy+1) * dy
        fy11 = c(3,ix+1,iy+1) * dy
        fxy00 = c(4,ix,  iy  ) * dx * dy
        fxy10 = c(4,ix+1,iy  ) * dx * dy
        fxy01 = c(4,ix,  iy+1) * dx * dy
        fxy11 = c(4,ix+1,iy+1) * dx * dy
        val =  H0x*H0y*f00 + H2x*H0y*f10 + H0x*H2y*f01 + H2x*H2y*f11 &
             + H1x*H0y*fx00 + H3x*H0y*fx10 + H1x*H2y*fx01 + H3x*H2y*fx11 &
             + H0x*H1y*fy00 + H2x*H1y*fy10 + H0x*H3y*fy01 + H2x*H3y*fy11 &
             + H1x*H1y*fxy00 + H3x*H1y*fxy10 + H1x*H3y*fxy01 + H3x*H3y*fxy11
    end function contract_il

    !> One bicubic Hermite evaluation of quantity c on the axes of container ax
    pure double precision function hermite_eval(ax, c, x, y) result(val)
        type(eos_table_container), intent(in) :: ax
        double precision, intent(in) :: c(:,:,:)
        double precision, intent(in) :: x, y
        integer :: ix, iy
        double precision :: tx, ty, dx, dy
        double precision :: H0x, H1x, H2x, H3x
        double precision :: H0y, H1y, H2y, H3y
        call locate_idx_axis1(ax, x, ix, tx, dx)
        call locate_idx_axis2(ax, y, iy, ty, dy)
        call cubic_basis(tx, H0x, H1x, H2x, H3x)
        call cubic_basis(ty, H0y, H1y, H2y, H3y)
        val = contract_il(c, ix, iy, dx, dy, H0x, H1x, H2x, H3x, H0y, H1y, H2y, H3y)
    end function hermite_eval

    !> Interleave a quantity's value and three derivative tables node by node,
    !> so the two corners of a cell row share one cache line
    subroutine interleave(f_t, fx_t, fy_t, fxy_t, c)
        type(eos_table_container), intent(in) :: f_t, fx_t, fy_t, fxy_t
        double precision, allocatable, intent(out) :: c(:,:,:)
        allocate(c(4, size(f_t%table,1), size(f_t%table,2)))
        c(1,:,:) = f_t%table
        c(2,:,:) = fx_t%table
        c(3,:,:) = fy_t%table
        c(4,:,:) = fxy_t%table
    end subroutine interleave

    !> Public lookups, code units in and out. Forward tables are on
    !> (log nH, log eint/nH), eintP and g1p on (log nH, log p/nH), eintT on
    !> (log nH, log T)
    double precision function entropy_T_from_nH_eint(log_nH_code, log_e_nh_code) result(T_code)
        double precision, intent(in) :: log_nH_code, log_e_nh_code
        T_code = hermite_eval(eos%Tfwd, il_Tfwd, log_nH_code, log_e_nh_code) / unit_temperature
    end function entropy_T_from_nH_eint

    pure double precision function entropy_p_nH_from_eint(log_nH_code, log_e_nh_code) result(p_nH_code)
        double precision, intent(in) :: log_nH_code, log_e_nh_code
        p_nH_code = hermite_eval(eos%pfwd, il_pfwd, log_nH_code, log_e_nh_code) &
            * ent_inv_up / 10.0d0**log_nH_code
    end function entropy_p_nH_from_eint

    !> Pressure in code units from (log nH, log eint/nH); the pfwd table holds p itself
    pure double precision function entropy_p_from_eint(log_nH_code, log_e_nh_code) result(p_code)
        double precision, intent(in) :: log_nH_code, log_e_nh_code
        p_code = hermite_eval(eos%pfwd, il_pfwd, log_nH_code, log_e_nh_code) * ent_inv_up
    end function entropy_p_from_eint

    double precision function entropy_y_from_nH_eint(log_nH_code, log_e_nh_code) result(y)
        double precision, intent(in) :: log_nH_code, log_e_nh_code
        y = hermite_eval(eos%neOnH, il_neOnH, log_nH_code, log_e_nh_code)
    end function entropy_y_from_nH_eint

    !> T and ne/nH at one point: Tfwd and neOnH share the forward grid, so the
    !> cell location and the basis are computed once
    subroutine entropy_T_and_y_from_nH_eint(log_nH_code, log_e_nh_code, T_code, y_out)
        double precision, intent(in)  :: log_nH_code, log_e_nh_code
        double precision, intent(out) :: T_code, y_out
        integer :: ix, iy
        double precision :: tx, ty, dx, dy
        double precision :: H0x, H1x, H2x, H3x
        double precision :: H0y, H1y, H2y, H3y
        call locate_idx_axis1(eos%Tfwd, log_nH_code, ix, tx, dx)
        call locate_idx_axis2(eos%Tfwd, log_e_nh_code, iy, ty, dy)
        call cubic_basis(tx, H0x, H1x, H2x, H3x)
        call cubic_basis(ty, H0y, H1y, H2y, H3y)
        T_code = contract_il(il_Tfwd, ix, iy, dx, dy, H0x, H1x, H2x, H3x, H0y, H1y, H2y, H3y) &
            / unit_temperature
        y_out  = contract_il(il_neOnH, ix, iy, dx, dy, H0x, H1x, H2x, H3x, H0y, H1y, H2y, H3y)
    end subroutine entropy_T_and_y_from_nH_eint

    !> eint/p in code units from (log nH, log p/nH). The table holds
    !> log10(eint/nH) in CGS
    double precision function entropy_eint_from_nH_p(log_nH_code, log_p_nH_code) result(ratio)
        double precision, intent(in) :: log_nH_code, log_p_nH_code
        double precision :: log_e_nh_code
        log_e_nh_code = hermite_eval(eos%eintP, il_eintP, log_nH_code, log_p_nH_code) - ent_shift_pe
        ratio = 10.0d0**(log_e_nh_code - log_p_nH_code)
    end function entropy_eint_from_nH_p

    !> Bisection inverse: log10(eint/nH) [code] such that the forward pressure
    !> matches p exactly; the eintP table only supplies the initial guess
    subroutine entropy_eint_from_p_bisect(log_nH_code, log_p_nH_code, log_eint_nH_code)
        double precision, intent(in)  :: log_nH_code, log_p_nH_code
        double precision, intent(out) :: log_eint_nH_code
        double precision :: ratio, guess, lo, hi, mid, f_lo, f_hi, f_mid
        double precision :: lim_lo, lim_hi, log_p_code
        integer :: iter

        ! eint-axis limits of the forward table, uniform or adaptive
        if (eos%pfwd%is_uniform) then
            lim_lo = eos%pfwd%var2_min
            lim_hi = eos%pfwd%var2_max
        else
            lim_lo = eos%pfwd%var2_nodes(1)
            lim_hi = eos%pfwd%var2_nodes(eos%pfwd%dim2)
        end if

        ratio = entropy_eint_from_nH_p(log_nH_code, log_p_nH_code)
        guess = dlog10(max(ratio, 1.0d-300)) + log_p_nH_code
        guess = max(lim_lo, min(lim_hi, guess))

        ! Bracket, expanding until the target is straddled
        log_p_code = log_nH_code + log_p_nH_code
        lo = max(lim_lo, guess - 2.0d-2)
        hi = min(lim_hi, guess + 2.0d-2)
        f_lo = resid(lo)
        f_hi = resid(hi)
        do iter = 1, 12
            if (f_lo*f_hi <= 0.0d0) exit
            lo = max(lim_lo, lo - 1.0d-1)
            hi = min(lim_hi, hi + 1.0d-1)
            f_lo = resid(lo)
            f_hi = resid(hi)
            if (lo <= lim_lo .and. hi >= lim_hi) exit
        end do

        ! No bracket, the target lies outside the table, so the table guess stands
        if (f_lo*f_hi > 0.0d0) then
            log_eint_nH_code = guess
            return
        end if

        do iter = 1, 40
            mid = 0.5d0*(lo + hi)
            f_mid = resid(mid)
            if (f_mid*f_lo <= 0.0d0) then
                hi = mid; f_hi = f_mid
            else
                lo = mid; f_lo = f_mid
            end if
            if (hi - lo < 1.0d-12) exit
        end do
        log_eint_nH_code = 0.5d0*(lo + hi)

    contains

        !> Residual in log10 p of the forward table at log10(eint/nH) = x
        double precision function resid(x)
            double precision, intent(in) :: x
            resid = dlog10(max(entropy_p_from_eint(log_nH_code, x), 1.0d-300)) - log_p_code
        end function resid

    end subroutine entropy_eint_from_p_bisect

    double precision function entropy_gamma1_from_nH_p(log_nH_code, log_p_nH_code) result(g1)
        double precision, intent(in) :: log_nH_code, log_p_nH_code
        g1 = hermite_eval(eos%g1p, il_g1p, log_nH_code, log_p_nH_code)
    end function entropy_gamma1_from_nH_p

    !> log10(eint/nH) in code units from (log nH, log T); the table holds it in CGS
    double precision function entropy_eint_from_nH_T(log_nH_code, log_T_code) result(log_e_nh_code)
        double precision, intent(in) :: log_nH_code, log_T_code
        log_e_nh_code = hermite_eval(eos%eintT, il_eintT, log_nH_code, log_T_code) - ent_shift_pe
    end function entropy_eint_from_nH_T


    !> Per-container code-unit shift + prepare for loaded entropy tables
    !> (moved from mod_eos_lte_tables; uses the shared infra there).
    subroutine entropy_table_prepare(tc)
        !> Bundle ensure_axis_nodes + precompute_step_inv + build_guards + validate,
        !> the identical four-step setup each loaded entropy table container needs.
        type(eos_table_container), intent(inout) :: tc
        call ensure_axis_nodes(tc)
        call precompute_step_inv(tc)
        call eos_build_guards(tc)
        call eos_validate_table(tc, trim(tc%filename))
    end subroutine entropy_table_prepare


    subroutine entropy_shift_prepare(tc, axis2_is_T)
        type(eos_table_container), intent(inout) :: tc
        logical, intent(in) :: axis2_is_T
        if (axis2_is_T) then
            call shift_axis_to_code_T(tc)
        else
            call shift_axis_to_code(tc)
        end if
        call entropy_table_prepare(tc)
    end subroutine entropy_shift_prepare

end module mod_eos_lte_entropy
!> Needs a line after to pass the preprocessor
