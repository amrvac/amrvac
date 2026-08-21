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
module mod_eos_LTE_entropy
    use mod_global_parameters
    use mod_comm_lib,      only: mpistop
    use mod_eos_container, only: eos_table_container, eos
    use mod_eos_LTE_tables, only: load_tables_LTE, precompute_FI_bypass_constants, &
         ensure_axis_nodes, shift_axis_to_code, shift_axis_to_code_T, &
         eos_build_guards, eos_validate_table
    use mod_eos_interp, only: precompute_step_inv
    implicit none
    private

    double precision, parameter :: LN10 = 2.302585092994046d0

    public :: load_entropy_LTE, finalise_entropy_LTE

    !> Forward (log nH, log eint/nH) -> thermodynamic state
    public :: entropy_T_from_nH_eint
    public :: entropy_p_nH_from_eint
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
        !> The entropy method never loads eos%T, but several callers read
        !> eos%T%var{1,2}_{min,max} (the cold-cell eint floor, FI-fallback
        !> checks). Source them from eos%Tfwd, on the same forward grid.
        eos%T%var2_min = eos%Tfwd%var2_min
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

    pure subroutine locate_idx_axis1(tab, val, ix, t, h)
        type(eos_table_container), intent(in) :: tab
        double precision, intent(in)  :: val
        integer,          intent(out) :: ix
        double precision, intent(out) :: t, h
        if (tab%is_uniform) then
            call locate_idx_uniform(tab%dim1, tab%var1_min, tab%var1_max, val, ix, t, h)
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

    !> Bicubic Hermite evaluator: value only.
    !
    !> 4 stored derivatives per corner (f, fx, fy, fxy) x 4 corners = 16
    !> conditions; fully determines the 16-coefficient degree-(3,3)
    !> polynomial. The stored fx, fy, fxy are physical-axis derivatives;
    !> we multiply by the local cell widths to convert to unit-interval
    !> basis coordinates.
    pure subroutine bicubic_hermite_eval(f_t, fx_t, fy_t, fxy_t, x, y, val)
        type(eos_table_container), intent(in) :: f_t, fx_t, fy_t, fxy_t
        double precision, intent(in)  :: x, y
        double precision, intent(out) :: val
        integer :: ix, iy
        double precision :: tx, ty, dx, dy
        double precision :: H0x, H1x, H2x, H3x
        double precision :: H0y, H1y, H2y, H3y
        double precision :: f00, f10, f01, f11
        double precision :: fx00, fx10, fx01, fx11
        double precision :: fy00, fy10, fy01, fy11
        double precision :: fxy00, fxy10, fxy01, fxy11

        call locate_idx_axis1(f_t, x, ix, tx, dx)
        call locate_idx_axis2(f_t, y, iy, ty, dy)

        call cubic_basis(tx, H0x, H1x, H2x, H3x)
        call cubic_basis(ty, H0y, H1y, H2y, H3y)

        f00 = f_t%table(ix,   iy  );    f10 = f_t%table(ix+1, iy  )
        f01 = f_t%table(ix,   iy+1);    f11 = f_t%table(ix+1, iy+1)
        fx00 = fx_t%table(ix,   iy  ) * dx
        fx10 = fx_t%table(ix+1, iy  ) * dx
        fx01 = fx_t%table(ix,   iy+1) * dx
        fx11 = fx_t%table(ix+1, iy+1) * dx
        fy00 = fy_t%table(ix,   iy  ) * dy
        fy10 = fy_t%table(ix+1, iy  ) * dy
        fy01 = fy_t%table(ix,   iy+1) * dy
        fy11 = fy_t%table(ix+1, iy+1) * dy
        fxy00 = fxy_t%table(ix,   iy  ) * dx * dy
        fxy10 = fxy_t%table(ix+1, iy  ) * dx * dy
        fxy01 = fxy_t%table(ix,   iy+1) * dx * dy
        fxy11 = fxy_t%table(ix+1, iy+1) * dx * dy

        val =  H0x*H0y*f00 + H2x*H0y*f10 + H0x*H2y*f01 + H2x*H2y*f11 &
             + H1x*H0y*fx00 + H3x*H0y*fx10 + H1x*H2y*fx01 + H3x*H2y*fx11 &
             + H0x*H1y*fy00 + H2x*H1y*fy10 + H0x*H3y*fy01 + H2x*H3y*fy11 &
             + H1x*H1y*fxy00 + H3x*H1y*fxy10 + H1x*H3y*fxy01 + H3x*H3y*fxy11
    end subroutine bicubic_hermite_eval

    !> Public wrappers -- code-unit out.
    double precision function entropy_T_from_nH_eint(Tfwd, Tfwd_x, Tfwd_y, Tfwd_xy, &
                                                           log_nH_code, log_e_nh_code) &
                                                           result(T_code)
        type(eos_table_container), intent(in) :: Tfwd, Tfwd_x, Tfwd_y, Tfwd_xy
        double precision, intent(in) :: log_nH_code, log_e_nh_code
        double precision :: T_cgs
        call bicubic_hermite_eval(Tfwd, Tfwd_x, Tfwd_y, Tfwd_xy, &
                                    log_nH_code, log_e_nh_code, T_cgs)
        T_code = T_cgs / unit_temperature
    end function entropy_T_from_nH_eint

    pure double precision function entropy_p_nH_from_eint(pfwd, pfwd_x, pfwd_y, pfwd_xy, &
                                                       log_nH_code, log_e_nh_code) &
                                                       result(p_nH_code)
        type(eos_table_container), intent(in) :: pfwd, pfwd_x, pfwd_y, pfwd_xy
        double precision, intent(in) :: log_nH_code, log_e_nh_code
        double precision :: p_cgs, nH_cgs
        call bicubic_hermite_eval(pfwd, pfwd_x, pfwd_y, pfwd_xy, &
                                    log_nH_code, log_e_nh_code, p_cgs)
        nH_cgs    = 10.0d0**(log_nH_code + dlog10(unit_numberdensity))
        p_nH_code = (p_cgs / nH_cgs) * unit_numberdensity / unit_pressure
    end function entropy_p_nH_from_eint

    !> Forward: ne/nH from (log nH, log eint/nH). Single bicubic Hermite
    !> eval of the neOnH table.
    double precision function entropy_y_from_nH_eint(neOnH, neOnH_x, neOnH_y, &
                                                          neOnH_xy, log_nH_code, &
                                                          log_e_nh_code) result(y)
        type(eos_table_container), intent(in) :: neOnH, neOnH_x, neOnH_y, neOnH_xy
        double precision, intent(in) :: log_nH_code, log_e_nh_code
        call bicubic_hermite_eval(neOnH, neOnH_x, neOnH_y, neOnH_xy, &
                                    log_nH_code, log_e_nh_code, y)
    end function entropy_y_from_nH_eint

    !> Forward: T and y from (log nH, log eint/nH) -- combined for the
    !> update_eos_LTE hot path that needs both at once.
    subroutine entropy_T_and_y_from_nH_eint(Tfwd, Tfwd_x, Tfwd_y, Tfwd_xy, &
                                                  neOnH, neOnH_x, neOnH_y, neOnH_xy, &
                                                  log_nH_code, log_e_nh_code, &
                                                  T_code, y_out)
        !> Fused (T, y) lookup: shares ONE cell-location call between the two
        !> bicubic Hermite evaluations. Tfwd and neOnH live on the same
        !> adaptive (lr, le) grid (same axis nodes), so the (ix, iy, tx, ty,
        !> dx, dy) coordinates and the cubic-basis values are identical for
        !> both. We compute them once and re-use, saving the locate work
        !> (binary search on adaptive axes) and the basis evaluations.
        type(eos_table_container), intent(in) :: Tfwd, Tfwd_x, Tfwd_y, Tfwd_xy
        type(eos_table_container), intent(in) :: neOnH, neOnH_x, neOnH_y, neOnH_xy
        double precision, intent(in)  :: log_nH_code, log_e_nh_code
        double precision, intent(out) :: T_code, y_out
        double precision :: T_cgs
        integer :: ix, iy
        double precision :: tx, ty, dx, dy
        double precision :: H0x, H1x, H2x, H3x
        double precision :: H0y, H1y, H2y, H3y

        call locate_idx_axis1(Tfwd, log_nH_code, ix, tx, dx)
        call locate_idx_axis2(Tfwd, log_e_nh_code, iy, ty, dy)
        call cubic_basis(tx, H0x, H1x, H2x, H3x)
        call cubic_basis(ty, H0y, H1y, H2y, H3y)

        T_cgs = contract_value(Tfwd, Tfwd_x, Tfwd_y, Tfwd_xy, ix, iy, &
                                 dx, dy, H0x, H1x, H2x, H3x, H0y, H1y, H2y, H3y)
        y_out = contract_value(neOnH, neOnH_x, neOnH_y, neOnH_xy, ix, iy, &
                                 dx, dy, H0x, H1x, H2x, H3x, H0y, H1y, H2y, H3y)
        T_code = T_cgs / unit_temperature
    end subroutine entropy_T_and_y_from_nH_eint

    !> Bicubic Hermite contraction given pre-computed cell coordinates and
    !> basis values. Lets multiple quantities at the SAME query point share
    !> the cell-location and basis-evaluation work. Returns value only.
    pure double precision function contract_value(f_t, fx_t, fy_t, fxy_t, &
                                                   ix, iy, dx, dy, &
                                                   H0x, H1x, H2x, H3x, &
                                                   H0y, H1y, H2y, H3y) result(val)
        !> Shared-locate bicubic-Hermite contraction for fused multi-quantity
        !> lookups at the same query point.
        type(eos_table_container), intent(in) :: f_t, fx_t, fy_t, fxy_t
        integer, intent(in) :: ix, iy
        double precision, intent(in) :: dx, dy
        double precision, intent(in) :: H0x, H1x, H2x, H3x
        double precision, intent(in) :: H0y, H1y, H2y, H3y
        double precision :: f00, f10, f01, f11
        double precision :: fx00, fx10, fx01, fx11
        double precision :: fy00, fy10, fy01, fy11
        double precision :: fxy00, fxy10, fxy01, fxy11
        f00 = f_t%table(ix,   iy  );    f10 = f_t%table(ix+1, iy  )
        f01 = f_t%table(ix,   iy+1);    f11 = f_t%table(ix+1, iy+1)
        fx00 = fx_t%table(ix,   iy  ) * dx
        fx10 = fx_t%table(ix+1, iy  ) * dx
        fx01 = fx_t%table(ix,   iy+1) * dx
        fx11 = fx_t%table(ix+1, iy+1) * dx
        fy00 = fy_t%table(ix,   iy  ) * dy
        fy10 = fy_t%table(ix+1, iy  ) * dy
        fy01 = fy_t%table(ix,   iy+1) * dy
        fy11 = fy_t%table(ix+1, iy+1) * dy
        fxy00 = fxy_t%table(ix,   iy  ) * dx * dy
        fxy10 = fxy_t%table(ix+1, iy  ) * dx * dy
        fxy01 = fxy_t%table(ix,   iy+1) * dx * dy
        fxy11 = fxy_t%table(ix+1, iy+1) * dx * dy
        val =  H0x*H0y*f00 + H2x*H0y*f10 + H0x*H2y*f01 + H2x*H2y*f11 &
             + H1x*H0y*fx00 + H3x*H0y*fx10 + H1x*H2y*fx01 + H3x*H2y*fx11 &
             + H0x*H1y*fy00 + H2x*H1y*fy10 + H0x*H3y*fy01 + H2x*H3y*fy11 &
             + H1x*H1y*fxy00 + H3x*H1y*fxy10 + H1x*H3y*fxy01 + H3x*H3y*fxy11
    end function contract_value

    !> Inverse: log10(eint/nH) from (log nH, log p/nH). Single bicubic
    !> Hermite eval of the eintP table. Replaces p2eint bisection.
    !
    !> Returns the eint/p ratio for drop-in compatibility with the
    !> existing p2eint_from_nH_p call sites which use it as
    !> eint = p * (returned ratio).
    double precision function entropy_eint_from_nH_p(eintP, eintP_x, eintP_y, eintP_xy, &
                                                       log_nH_code, log_p_nH_code) result(ratio)
        !> Returns eint/p in dimensionless code units. The stored table value is
        !> log10(eint/nH) in CGS; we convert to code (axis-2 shift = the same
        !> log10(unit_pressure/unit_numberdensity) as for log(p/nH) and log(eint/nH))
        !> before forming the ratio with log_p_nH_code, otherwise the units mix.
        type(eos_table_container), intent(in) :: eintP, eintP_x, eintP_y, eintP_xy
        double precision, intent(in) :: log_nH_code, log_p_nH_code
        double precision :: log_e_nh_cgs, log_e_nh_code
        call bicubic_hermite_eval(eintP, eintP_x, eintP_y, eintP_xy, &
                                    log_nH_code, log_p_nH_code, log_e_nh_cgs)
        log_e_nh_code = log_e_nh_cgs - dlog10(unit_pressure / unit_numberdensity)
        ratio = 10.0d0**(log_e_nh_code - log_p_nH_code)
    end function entropy_eint_from_nH_p

    !> Bisection inverse for the entropy table set: find log10(eint/nH) [code] such that
    !> p(nH, eint) = p_target exactly, using the forward pfwd table. The eintP table only
    !> supplies the initial guess. This is the entropy-mode counterpart of
    !> eint_from_p_bisect (mod_eos_LTE), which is bound to eos%log_p - a table that is not
    !> loaded when eos_method='entropy' (its use there segfaulted on an unallocated array).
    subroutine entropy_eint_from_p_bisect(pfwd, pfwd_x, pfwd_y, pfwd_xy, &
                                          eintP, eintP_x, eintP_y, eintP_xy, &
                                          log_nH_code, log_p_nH_code, log_eint_nH_code)
        type(eos_table_container), intent(in) :: pfwd, pfwd_x, pfwd_y, pfwd_xy
        type(eos_table_container), intent(in) :: eintP, eintP_x, eintP_y, eintP_xy
        double precision, intent(in)  :: log_nH_code, log_p_nH_code
        double precision, intent(out) :: log_eint_nH_code
        double precision :: ratio, guess, lo, hi, mid, f_lo, f_hi, f_mid
        double precision :: lim_lo, lim_hi
        integer :: iter

        !> eint-axis limits of the forward table (uniform vs adaptive grid)
        if (pfwd%is_uniform) then
            lim_lo = pfwd%var2_min
            lim_hi = pfwd%var2_max
        else
            lim_lo = pfwd%var2_nodes(1)
            lim_hi = pfwd%var2_nodes(pfwd%dim2)
        end if

        !> initial guess from the inverse table
        ratio = entropy_eint_from_nH_p(eintP, eintP_x, eintP_y, eintP_xy, &
                                       log_nH_code, log_p_nH_code)
        guess = dlog10(max(ratio, 1.0d-300)) + log_p_nH_code
        guess = max(lim_lo, min(lim_hi, guess))

        !> bracket, expanding until the target is straddled
        lo = max(lim_lo, guess - 2.0d-2)
        hi = min(lim_hi, guess + 2.0d-2)
        f_lo = dlog10(max(entropy_p_nH_from_eint(pfwd, pfwd_x, pfwd_y, pfwd_xy, &
                          log_nH_code, lo), 1.0d-300)) - log_p_nH_code
        f_hi = dlog10(max(entropy_p_nH_from_eint(pfwd, pfwd_x, pfwd_y, pfwd_xy, &
                          log_nH_code, hi), 1.0d-300)) - log_p_nH_code
        do iter = 1, 12
            if (f_lo*f_hi <= 0.0d0) exit
            lo = max(lim_lo, lo - 1.0d-1)
            hi = min(lim_hi, hi + 1.0d-1)
            f_lo = dlog10(max(entropy_p_nH_from_eint(pfwd, pfwd_x, pfwd_y, pfwd_xy, &
                              log_nH_code, lo), 1.0d-300)) - log_p_nH_code
            f_hi = dlog10(max(entropy_p_nH_from_eint(pfwd, pfwd_x, pfwd_y, pfwd_xy, &
                              log_nH_code, hi), 1.0d-300)) - log_p_nH_code
            if (lo <= lim_lo .and. hi >= lim_hi) exit
        end do

        !> no bracket (target outside the table range) -> keep the table guess
        if (f_lo*f_hi > 0.0d0) then
            log_eint_nH_code = guess
            return
        end if

        do iter = 1, 40
            mid = 0.5d0*(lo + hi)
            f_mid = dlog10(max(entropy_p_nH_from_eint(pfwd, pfwd_x, pfwd_y, pfwd_xy, &
                               log_nH_code, mid), 1.0d-300)) - log_p_nH_code
            if (f_mid*f_lo <= 0.0d0) then
                hi = mid; f_hi = f_mid
            else
                lo = mid; f_lo = f_mid
            end if
            if (hi - lo < 1.0d-12) exit
        end do
        log_eint_nH_code = 0.5d0*(lo + hi)
    end subroutine entropy_eint_from_p_bisect

    !> Inverse: Gamma_1 from (log nH, log p/nH). Bicubic Hermite of g1p.
    double precision function entropy_gamma1_from_nH_p(g1p, g1p_x, g1p_y, g1p_xy, &
                                                         log_nH_code, log_p_nH_code) &
                                                         result(g1)
        type(eos_table_container), intent(in) :: g1p, g1p_x, g1p_y, g1p_xy
        double precision, intent(in) :: log_nH_code, log_p_nH_code
        call bicubic_hermite_eval(g1p, g1p_x, g1p_y, g1p_xy, &
                                    log_nH_code, log_p_nH_code, g1)
    end function entropy_gamma1_from_nH_p

    !> Inverse: log10(eint/nH) from (log nH, log T). Single bicubic
    !> Hermite eval of the eintT table.
    double precision function entropy_eint_from_nH_T(eintT, eintT_x, eintT_y, eintT_xy, &
                                                       log_nH_code, log_T_code) &
                                                       result(log_e_nh_code)
        !> Stored table value is log10(eint/nH) in CGS; convert to code units
        !> before returning so the caller can use it as a code-unit log directly.
        type(eos_table_container), intent(in) :: eintT, eintT_x, eintT_y, eintT_xy
        double precision, intent(in) :: log_nH_code, log_T_code
        double precision :: log_e_nh_cgs
        call bicubic_hermite_eval(eintT, eintT_x, eintT_y, eintT_xy, &
                                    log_nH_code, log_T_code, log_e_nh_cgs)
        log_e_nh_code = log_e_nh_cgs - dlog10(unit_pressure / unit_numberdensity)
    end function entropy_eint_from_nH_T


    !> Per-container code-unit shift + prepare for loaded entropy tables
    !> (moved from mod_eos_LTE_tables; uses the shared infra there).
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

end module mod_eos_LTE_entropy
!> Needs a line after to pass the preprocessor
