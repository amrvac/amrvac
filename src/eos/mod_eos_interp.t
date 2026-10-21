!=============================================================================
!> Interpolation kernels for the EoS tables (pure math; no EoS state).
!>
!> Every routine operates on an eos_table_container passed by the caller:
!>   bicubic_lookup / bilinear_lookup  -- dispatch on tc%is_uniform and return
!>       the interpolated value (monotone bicubic PCHIP, or bilinear).
!>   interp_pchip_interleaved(_nu)      -- evaluate several quantities sharing
!>       one (nH, eint/nH) grid in a single cache-friendly pass.
!>   find_index_guard                   -- O(1) bracket lookup on a non-uniform
!>       axis (falls back to binary search if no guard array is built).
!> Uniform grids use an affine index; non-uniform (_nu) grids use explicit
!> node arrays. The PCHIP kernel is monotone (no overshoot on each 1D slice).
!=============================================================================
module mod_eos_interp
    use mod_eos_container, only: eos_table_container
    implicit none
    private

    public :: bicubic_lookup, bilinear_lookup
    public :: interp_pchip_interleaved, interp_pchip_interleaved_nu
    public :: find_index_guard, interp_clamped_bilinear_table
    public :: precompute_step_inv

contains

    !> The single monotone cubic Hermite (PCHIP-style) interval kernel for this
    !> module. Interpolates on [p1,p2] (t in [0,1]) from 4 uniformly spaced
    !> samples p0..p3: secant slopes d0,d1,d2; interior derivatives by harmonic
    !> mean with a monotonicity sign-clamp and a Hyman magnitude clamp (no
    !> overshoot on the interval). EVERY PCHIP slice in this module -- the
    !> uniform/non-uniform bicubic routines and both interleaved evaluators --
    !> routes through here, so the limiter exists in exactly one place.
    pure double precision function pchip_interval(p0, p1, p2, p3, t) result(v)
        double precision, intent(in) :: p0, p1, p2, p3, t
        double precision :: d0, d1, d2, m1, m2
        double precision :: tt, ttt, h00, h10, h01, h11
        double precision :: s, a1, a2, lim

        d0 = p1 - p0
        d1 = p2 - p1
        d2 = p3 - p2

        if (d1 == 0.0d0) then
            m1 = 0.0d0
            m2 = 0.0d0
        else
            if (d0*d1 <= 0.0d0) then
                m1 = 0.0d0
            else
                m1 = 2.0d0*d0*d1 / (d0 + d1)
            end if
            if (d1*d2 <= 0.0d0) then
                m2 = 0.0d0
            else
                m2 = 2.0d0*d1*d2 / (d1 + d2)
            end if
            s = sign(1.0d0, d1)
            a1 = s*m1
            a2 = s*m2
            if (a1 < 0.0d0) a1 = 0.0d0
            if (a2 < 0.0d0) a2 = 0.0d0
            lim = 3.0d0*abs(d1)
            if (a1 > lim) a1 = lim
            if (a2 > lim) a2 = lim
            m1 = s*a1
            m2 = s*a2
        end if

        tt  = t*t
        ttt = tt*t
        h00 = 2.0d0*ttt - 3.0d0*tt + 1.0d0
        h10 = ttt - 2.0d0*tt + t
        h01 = -2.0d0*ttt + 3.0d0*tt
        h11 = ttt - tt

        v = h00*p1 + h10*m1 + h01*p2 + h11*m2
    end function pchip_interval

    !> Fast adaptive index lookup. Returns the smallest 1-based ii such that
    !> nodes(ii) >= val. Falls back to binary search if the guard array is
    !> not built (M = 0). Equivalent in result to find_index_bsearch but O(1)
    !> when the guard is present.
    pure integer function find_index_guard(nodes, n, val, guard, M, scale_M) result(ix)
        use mod_lookup_table, only: find_index_bsearch
        double precision, intent(in) :: nodes(n), val, scale_M
        integer, intent(in) :: n, M
        integer, intent(in) :: guard(M)
        integer :: k
        if (M > 0) then
            k = 1 + min(M-1, max(0, int(floor((val - nodes(1)) * scale_M))))
            ix = guard(k)
            do while (ix < n .and. nodes(ix) < val)
                ix = ix + 1
            end do
        else
            ix = find_index_bsearch(nodes(1:n), val)
        end if
    end function find_index_guard

    pure double precision function interp_clamped_bilinear_table(vary, varx, table, nx, ny, vmin_y, vmax_y, vmin_x, vmax_x) result(z)
        double precision, intent(in) :: vary, varx
        integer, intent(in)  :: nx, ny
        double precision, intent(in) :: table(ny, nx)
        double precision, intent(in) :: vmin_y, vmax_y, vmin_x, vmax_x

        double precision :: fx, fy, tx, ty, rx, ry, xstep_inv, ystep_inv
        integer  :: ix, iy, ix1, iy1

        xstep_inv = dble(nx-1) / (vmax_x - vmin_x)
        ystep_inv = dble(ny-1) / (vmax_y - vmin_y)

        !> fractional positions
        fx = (varx - vmin_x) * xstep_inv
        fy = (vary - vmin_y) * ystep_inv

        !> clamp to [0, nx-1] and [0, ny-1]
        rx = max(0.0d0, min(fx, dble(nx-1)))
        ry = max(0.0d0, min(fy, dble(ny-1)))

        ix = int(rx)
        iy = int(ry)

        ix1 = min(ix+1, nx-1)
        iy1 = min(iy+1, ny-1)

        tx = rx - dble(ix)
        ty = ry - dble(iy)

        z = (1.0d0-ty) * ((1.0d0-tx)*table(iy+1, ix+1) + tx*table(iy+1, ix1+1)) &
        +        ty  * ((1.0d0-tx)*table(iy1+1,ix+1) + tx*table(iy1+1,ix1+1))
    end function interp_clamped_bilinear_table

    pure double precision function interp_clamped_monotone_bicubic_table( &
    vary, varx, table, nx, ny, vmin_y, vmax_y, vmin_x, vmax_x) result(z)

        !> Monotone bicubic (practical tensor-product variant):
        !>   - Do a *monotone* cubic Hermite (PCHIP-style) interpolation in x on 4 points
        !>     for each of 4 surrounding y-rows (j-1..j+2) -> gives 4 intermediate values.
        !>   - Then do the same monotone cubic Hermite interpolation in y on those 4 values.
        !>
        !> This construction is monotone along x-lines and y-lines (no overshoot on each 1D slice)
        !> https://jacobwilliams.github.io/PCHIP/
        !> https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html

        double precision, intent(in) :: vary, varx
        integer, intent(in)          :: nx, ny
        double precision, intent(in) :: table(ny, nx)
        double precision, intent(in) :: vmin_y, vmax_y, vmin_x, vmax_x

        double precision :: fx, fy, rx, ry, tx, ty, xstep_inv, ystep_inv
        integer :: ix, iy

        double precision :: g0, g1, g2, g3

        xstep_inv = dble(nx-1) / (vmax_x - vmin_x)
        ystep_inv = dble(ny-1) / (vmax_y - vmin_y)

        fx = (varx - vmin_x) * xstep_inv
        fy = (vary - vmin_y) * ystep_inv

        rx = max(0.0d0, min(fx, dble(nx-1)))
        ry = max(0.0d0, min(fy, dble(ny-1)))

        ix = int(rx)
        iy = int(ry)

        tx = rx - dble(ix)
        ty = ry - dble(iy)

        !> Four rows around iy: (iy-1, iy, iy+1, iy+2)
        g0 = x_interp_row(iy-1, ix, tx)
        g1 = x_interp_row(iy  , ix, tx)
        g2 = x_interp_row(iy+1, ix, tx)
        g3 = x_interp_row(iy+2, ix, tx)

        !> Now monotone cubic in y between g1 and g2 with neighbors g0,g3
        z = y_interp_from4(g0, g1, g2, g3, ty)

    contains

        pure integer function clampi(i, lo, hi) result(o)
            integer, intent(in) :: i, lo, hi
            o = max(lo, min(hi, i))
        end function clampi

        pure double precision function x_interp_row(j, i, t) result(v)
            !> Interpolate in x for fixed row j at cell starting index i (0-based),
            !> using points i-1, i, i+1, i+2 with clamped indexing.
            integer, intent(in) :: j, i
            double precision, intent(in) :: t
            integer :: i0, i1, i2, i3, jj
            double precision :: p0, p1, p2, p3

            jj = clampi(j, 0, ny-1)
            i0 = clampi(i-1, 0, nx-1)
            i1 = clampi(i  , 0, nx-1)
            i2 = clampi(i+1, 0, nx-1)
            i3 = clampi(i+2, 0, nx-1)

            p0 = table(jj+1, i0+1)
            p1 = table(jj+1, i1+1)
            p2 = table(jj+1, i2+1)
            p3 = table(jj+1, i3+1)

            v = pchip_interval(p0, p1, p2, p3, t)
        end function x_interp_row

        pure double precision function y_interp_from4(v0, v1, v2, v3, t) result(v)
            !> Interpolate in y between v1 and v2 using surrounding v0..v3, t in [0,1]
            double precision, intent(in) :: v0, v1, v2, v3, t
            v = pchip_interval(v0, v1, v2, v3, t)
        end function y_interp_from4

    end function interp_clamped_monotone_bicubic_table

    !> Dispatch wrapper: bicubic lookup that branches on the table's
    !> is_uniform flag. Hides the uniform-vs-adaptive distinction from
    !> call sites -- they pass (var1, var2, eos%<table>) and get the
    !> interpolated value, regardless of grid type.
    pure double precision function bicubic_lookup(var1, var2, tc) result(z)
        double precision, intent(in) :: var1, var2
        type(eos_table_container), intent(in) :: tc
        if (tc%is_uniform) then
            !> Existing uniform routine signature has (vary, varx, ...) where
            !> vary is the FIRST table axis (var1). Caller var1 -> vary,
            !> caller var2 -> varx. The existing argument ordering passes
            !> dim1, dim2 and the var1 / var2 bounds; the routine's internal
            !> "y/x" labelling is consistent with that ordering for square
            !> tables (which all our tables are).
            z = interp_clamped_monotone_bicubic_table(var1, var2, tc%table, &
                tc%dim1, tc%dim2, &
                tc%var1_min, tc%var1_max, tc%var2_min, tc%var2_max)
        else
            z = interp_clamped_monotone_bicubic_table_nu(var1, var2, tc%table, &
                tc%dim1, tc%dim2, tc%var1_nodes, tc%var2_nodes, &
                tc%guard_1, tc%guard_M_1, tc%guard_scale_1, &
                tc%guard_2, tc%guard_M_2, tc%guard_scale_2)
        end if
    end function bicubic_lookup

    !> Dispatch wrapper: bilinear lookup that branches on is_uniform.
    pure double precision function bilinear_lookup(var1, var2, tc) result(z)
        double precision, intent(in) :: var1, var2
        type(eos_table_container), intent(in) :: tc
        if (tc%is_uniform) then
            z = interp_clamped_bilinear_table(var1, var2, tc%table, &
                tc%dim1, tc%dim2, &
                tc%var1_min, tc%var1_max, tc%var2_min, tc%var2_max)
        else
            z = interp_clamped_bilinear_table_nu(var1, var2, tc%table, &
                tc%dim1, tc%dim2, tc%var1_nodes, tc%var2_nodes, &
                tc%guard_1, tc%guard_M_1, tc%guard_scale_1, &
                tc%guard_2, tc%guard_M_2, tc%guard_scale_2)
        end if
    end function bilinear_lookup

    !> Non-uniform-grid bilinear lookup. Same semantics as
    !> interp_clamped_bilinear_table but with explicit per-axis node arrays.
    !> Use binary search to locate the enclosing cell, then compute the
    !> local fractional coordinate from the actual node positions. Falls
    !> back to the boundary node value outside the grid (clamped).
    !>
    !> Convention (matches struct fields var1_nodes / var2_nodes):
    !>   var1 has dim1 nodes, var2 has dim2 nodes.
    !>   Table is stored as table(dim1, dim2) -- Fortran column-major,
    !>   matching how generate_state_tables.py writes it (data.T.tofile).
    !>   Element (i1, i2) of the table corresponds to node positions
    !>   (var1_nodes(i1), var2_nodes(i2)).
    pure double precision function interp_clamped_bilinear_table_nu( &
        var1, var2, table, dim1, dim2, var1_nodes, var2_nodes, &
        guard_1, M_1, scale_1, guard_2, M_2, scale_2) result(z)
        double precision, intent(in) :: var1, var2
        integer, intent(in)          :: dim1, dim2
        double precision, intent(in) :: table(dim1, dim2)
        double precision, intent(in) :: var1_nodes(dim1), var2_nodes(dim2)
        integer, intent(in)          :: M_1, M_2
        double precision, intent(in) :: scale_1, scale_2
        integer, intent(in)          :: guard_1(M_1), guard_2(M_2)

        integer          :: i1, i2, i1p, i2p, ii
        double precision :: t1, t2

        !> Axis-1 cell location
        if (var1 <= var1_nodes(1)) then
            i1 = 0; t1 = 0.0d0
        else if (var1 >= var1_nodes(dim1)) then
            i1 = dim1 - 2; t1 = 1.0d0
        else
            ii = find_index_guard(var1_nodes, dim1, var1, guard_1, M_1, scale_1)
            i1 = max(0, min(ii - 2, dim1 - 2))            !> 0-based cell index
            t1 = (var1 - var1_nodes(i1+1)) / (var1_nodes(i1+2) - var1_nodes(i1+1))
            t1 = max(0.0d0, min(t1, 1.0d0))
        end if

        !> Axis-2 cell location
        if (var2 <= var2_nodes(1)) then
            i2 = 0; t2 = 0.0d0
        else if (var2 >= var2_nodes(dim2)) then
            i2 = dim2 - 2; t2 = 1.0d0
        else
            ii = find_index_guard(var2_nodes, dim2, var2, guard_2, M_2, scale_2)
            i2 = max(0, min(ii - 2, dim2 - 2))
            t2 = (var2 - var2_nodes(i2+1)) / (var2_nodes(i2+2) - var2_nodes(i2+1))
            t2 = max(0.0d0, min(t2, 1.0d0))
        end if

        i1p = min(i1+1, dim1-1)
        i2p = min(i2+1, dim2-1)

        z = (1.0d0-t2) * ((1.0d0-t1)*table(i1+1, i2+1) + t1*table(i1p+1, i2+1)) &
        +        t2  * ((1.0d0-t1)*table(i1+1, i2p+1) + t1*table(i1p+1, i2p+1))
    end function interp_clamped_bilinear_table_nu

    !> Non-uniform-grid monotone bicubic lookup. Same semantics as
    !> interp_clamped_monotone_bicubic_table but with explicit per-axis
    !> node arrays. The PCHIP kernel itself is unchanged: it uses the
    !> "uniform" secant-slope formulation operating on the local 4-tuple
    !> of values (this is the standard pragmatic approach for adaptive
    !> grids with rectangular node structure -- strictly suboptimal vs.
    !> proper non-uniform PCHIP but empirically near-optimal for grids
    !> designed by curvature equidistribution).
    pure double precision function interp_clamped_monotone_bicubic_table_nu( &
        var1, var2, table, dim1, dim2, var1_nodes, var2_nodes, &
        guard_1, M_1, scale_1, guard_2, M_2, scale_2) result(z)
        double precision, intent(in) :: var1, var2
        integer, intent(in)          :: dim1, dim2
        double precision, intent(in) :: table(dim1, dim2)
        double precision, intent(in) :: var1_nodes(dim1), var2_nodes(dim2)
        integer, intent(in)          :: M_1, M_2
        double precision, intent(in) :: scale_1, scale_2
        integer, intent(in)          :: guard_1(M_1), guard_2(M_2)

        integer          :: i1, i2, ii
        double precision :: t1, t2
        double precision :: g0, g1, g2, g3

        !> Axis-1 cell location
        if (var1 <= var1_nodes(1)) then
            i1 = 0; t1 = 0.0d0
        else if (var1 >= var1_nodes(dim1)) then
            i1 = dim1 - 2; t1 = 1.0d0
        else
            ii = find_index_guard(var1_nodes, dim1, var1, guard_1, M_1, scale_1)
            i1 = max(0, min(ii - 2, dim1 - 2))
            t1 = (var1 - var1_nodes(i1+1)) / (var1_nodes(i1+2) - var1_nodes(i1+1))
            t1 = max(0.0d0, min(t1, 1.0d0))
        end if

        !> Axis-2 cell location
        if (var2 <= var2_nodes(1)) then
            i2 = 0; t2 = 0.0d0
        else if (var2 >= var2_nodes(dim2)) then
            i2 = dim2 - 2; t2 = 1.0d0
        else
            ii = find_index_guard(var2_nodes, dim2, var2, guard_2, M_2, scale_2)
            i2 = max(0, min(ii - 2, dim2 - 2))
            t2 = (var2 - var2_nodes(i2+1)) / (var2_nodes(i2+2) - var2_nodes(i2+1))
            t2 = max(0.0d0, min(t2, 1.0d0))
        end if

        !> Four columns around i2: (i2-1, i2, i2+1, i2+2), interpolate
        !> in axis-1 along each fixed-i2 column, then pchip in axis-2.
        g0 = ax1_interp_col_nu(i2-1, i1, t1)
        g1 = ax1_interp_col_nu(i2  , i1, t1)
        g2 = ax1_interp_col_nu(i2+1, i1, t1)
        g3 = ax1_interp_col_nu(i2+2, i1, t1)

        z = pchip_interval(g0, g1, g2, g3, t2)

    contains

        pure integer function clampi_nu(i, lo, hi) result(o)
            integer, intent(in) :: i, lo, hi
            o = max(lo, min(hi, i))
        end function clampi_nu

        pure double precision function ax1_interp_col_nu(j2, i1_cell, t) result(v)
            !> Interpolate in axis-1 at fixed axis-2 column j2 (0-based,
            !> clamped) using axis-1 cell starting index i1_cell (0-based),
            !> with the 4-point stencil i1_cell-1, i1_cell, i1_cell+1, i1_cell+2.
            integer, intent(in) :: j2, i1_cell
            double precision, intent(in) :: t
            integer :: a0, a1, a2, a3, b
            double precision :: p0, p1, p2, p3

            b  = clampi_nu(j2, 0, dim2-1)
            a0 = clampi_nu(i1_cell-1, 0, dim1-1)
            a1 = clampi_nu(i1_cell  , 0, dim1-1)
            a2 = clampi_nu(i1_cell+1, 0, dim1-1)
            a3 = clampi_nu(i1_cell+2, 0, dim1-1)

            p0 = table(a0+1, b+1)
            p1 = table(a1+1, b+1)
            p2 = table(a2+1, b+1)
            p3 = table(a3+1, b+1)

            v = pchip_interval(p0, p1, p2, p3, t)
        end function ax1_interp_col_nu

    end function interp_clamped_monotone_bicubic_table_nu

    !> Interleaved PCHIP: evaluate N quantities at the same (vary, varx) point.
    !> Table layout: table_il(nq, ny, nx) where nq quantities share the same grid.
    !> All nq values at each grid point are contiguous in memory (cache-optimal).
    subroutine interp_pchip_interleaved(vary, varx, table_il, nq, nx, ny, &
        vmin_y, vmax_y, vmin_x, vmax_x, results)
        double precision, intent(in) :: vary, varx
        integer, intent(in)          :: nq, nx, ny
        double precision, intent(in) :: table_il(nq, ny, nx)
        double precision, intent(in) :: vmin_y, vmax_y, vmin_x, vmax_x
        double precision, intent(out) :: results(nq)

        double precision :: fx, fy, rx, ry, tx, ty, xstep_inv, ystep_inv
        integer :: ix, iy, q
        integer :: i0, i1, i2, i3, j0, j1, j2, j3
        double precision :: g0(nq), g1(nq), g2(nq), g3(nq)
        double precision :: p0, p1, p2, p3

        xstep_inv = dble(nx-1) / (vmax_x - vmin_x)
        ystep_inv = dble(ny-1) / (vmax_y - vmin_y)

        fx = (varx - vmin_x) * xstep_inv
        fy = (vary - vmin_y) * ystep_inv

        rx = max(0.0d0, min(fx, dble(nx-1)))
        ry = max(0.0d0, min(fy, dble(ny-1)))

        ix = int(rx); iy = int(ry)
        tx = rx - dble(ix); ty = ry - dble(iy)

        ! Clamped column indices (shared across all rows and quantities)
        i0 = max(0, min(nx-1, ix-1))
        i1 = max(0, min(nx-1, ix))
        i2 = max(0, min(nx-1, ix+1))
        i3 = max(0, min(nx-1, ix+2))

        ! For each of 4 y-rows, interpolate in x for ALL quantities
        do j0 = 0, 3
            j1 = max(0, min(ny-1, iy - 1 + j0))
            do q = 1, nq
                p0 = table_il(q, j1+1, i0+1)
                p1 = table_il(q, j1+1, i1+1)
                p2 = table_il(q, j1+1, i2+1)
                p3 = table_il(q, j1+1, i3+1)
                select case(j0)
                case(0); g0(q) = pchip_interval(p0, p1, p2, p3, tx)
                case(1); g1(q) = pchip_interval(p0, p1, p2, p3, tx)
                case(2); g2(q) = pchip_interval(p0, p1, p2, p3, tx)
                case(3); g3(q) = pchip_interval(p0, p1, p2, p3, tx)
                end select
            end do
        end do

        ! Final PCHIP interpolation in y for each quantity
        do q = 1, nq
            results(q) = pchip_interval(g0(q), g1(q), g2(q), g3(q), ty)
        end do

    end subroutine interp_pchip_interleaved

    !> Adaptive-grid sibling of interp_pchip_interleaved.  Replaces the
    !> affine index calculation `(x - x_min) * step_inv` with binary
    !> searches on the explicit node arrays var{1,2}_nodes (Q axis-1
    !> binary searches: O(log_2 N) comparisons each, total ~16 cmps for
    !> N = 256).  The PCHIP kernel itself is unchanged: it operates on
    !> the local 4-tuple of values regardless of grid spacing.  The
    !> stride-1 access pattern over the q-slot dimension is preserved,
    !> so cache behaviour matches the uniform variant for the inner kernel.
    subroutine interp_pchip_interleaved_nu(vary, varx, table_il, nq, nx, ny, &
        var1_nodes, var2_nodes, &
        guard_1, M_1, scale_1, guard_2, M_2, scale_2, results)
        double precision, intent(in) :: vary, varx
        integer, intent(in)          :: nq, nx, ny
        double precision, intent(in) :: table_il(nq, ny, nx)
        double precision, intent(in) :: var1_nodes(ny), var2_nodes(nx)
        integer, intent(in)          :: M_1, M_2
        double precision, intent(in) :: scale_1, scale_2
        integer, intent(in)          :: guard_1(M_1), guard_2(M_2)
        double precision, intent(out) :: results(nq)

        double precision :: tx, ty
        integer :: ix, iy, q, ii
        integer :: i0, i1, i2, i3, j0, j1, j2, j3
        double precision :: g0(nq), g1(nq), g2(nq), g3(nq)
        double precision :: p0, p1, p2, p3

        ! Axis-1 (vary) cell index + local fractional coordinate
        if (vary <= var1_nodes(1)) then
            iy = 0; ty = 0.0d0
        else if (vary >= var1_nodes(ny)) then
            iy = ny - 2; ty = 1.0d0
        else
            ii = find_index_guard(var1_nodes, ny, vary, guard_1, M_1, scale_1)
            iy = max(0, min(ii - 2, ny - 2))
            ty = (vary - var1_nodes(iy+1)) &
               / (var1_nodes(iy+2) - var1_nodes(iy+1))
            ty = max(0.0d0, min(ty, 1.0d0))
        end if

        ! Axis-2 (varx) cell index + local fractional coordinate
        if (varx <= var2_nodes(1)) then
            ix = 0; tx = 0.0d0
        else if (varx >= var2_nodes(nx)) then
            ix = nx - 2; tx = 1.0d0
        else
            ii = find_index_guard(var2_nodes, nx, varx, guard_2, M_2, scale_2)
            ix = max(0, min(ii - 2, nx - 2))
            tx = (varx - var2_nodes(ix+1)) &
               / (var2_nodes(ix+2) - var2_nodes(ix+1))
            tx = max(0.0d0, min(tx, 1.0d0))
        end if

        i0 = max(0, min(nx-1, ix-1))
        i1 = max(0, min(nx-1, ix))
        i2 = max(0, min(nx-1, ix+1))
        i3 = max(0, min(nx-1, ix+2))

        ! Same kernel as the uniform variant -- only the index calculation differs
        do j0 = 0, 3
            j1 = max(0, min(ny-1, iy - 1 + j0))
            do q = 1, nq
                p0 = table_il(q, j1+1, i0+1)
                p1 = table_il(q, j1+1, i1+1)
                p2 = table_il(q, j1+1, i2+1)
                p3 = table_il(q, j1+1, i3+1)
                select case(j0)
                case(0); g0(q) = pchip_interval(p0, p1, p2, p3, tx)
                case(1); g1(q) = pchip_interval(p0, p1, p2, p3, tx)
                case(2); g2(q) = pchip_interval(p0, p1, p2, p3, tx)
                case(3); g3(q) = pchip_interval(p0, p1, p2, p3, tx)
                end select
            end do
        end do

        do q = 1, nq
            results(q) = pchip_interval(g0(q), g1(q), g2(q), g3(q), ty)
        end do

    end subroutine interp_pchip_interleaved_nu

    subroutine precompute_step_inv(tc)
        type(eos_table_container), intent(inout) :: tc
        if (allocated(tc%table)) then
            tc%step_inv_1 = dble(tc%dim1-1) / (tc%var1_max - tc%var1_min)
            tc%step_inv_2 = dble(tc%dim2-1) / (tc%var2_max - tc%var2_min)
        end if
    end subroutine precompute_step_inv

end module mod_eos_interp
!> Needs a line after to pass the preprocessor
