!=============================================================================
!> 'state' method (eos_method='state', legacy 'tables') of the LTE EoS family.
!>
!> The direct per-quantity tabulation of the Saha ionisation-equilibrium state:
!> loads the T/neOnH (and ionE p2eint) tables, shifts them to code units, and
!> builds the derived tables the runtime needs (gamma1, log_p, p_over_nH, the
!> interleaved fast-path block, the pressure-indexed gamma1, and the bisected
!> p2eint / eint-from-T inverses). The shared binary-file I/O, axis/guard/validate
!> machinery and the FI-bypass constants live in mod_eos_lte_tables (used by all
!> LTE methods); this module owns only the state-method build + finalise. Runs
!> once, during eos_init_LTE / eos_finalise_LTE.
!=============================================================================
module mod_eos_lte_state
    use mod_global_parameters
    use mod_eos_container
    use mod_eos_interp
    use mod_eos_lte_tables, only: load_tables_LTE, try_load_tables_LTE, &
         ensure_axis_nodes, eos_build_guards, eos_validate_table, &
         precompute_FI_bypass_constants
    use mod_comm_lib, only: mpistop
    implicit none
    private

    !> state arms of the eos_init / eos_finalise dispatchers (mod_eos_lte)
    public :: load_state_LTE, finalise_state_LTE

contains

    !> 'tables' method load: T + neOnH, plus the ionE inverse/derived tables
    !> (p2eint loaded; gamma1 / eint_from_T loaded if available, else built later).
    !> Called from eos_init_LTE for eos_method='tables'.
    subroutine load_state_LTE()
        call load_tables_LTE("T")
        call load_tables_LTE("neOnH")
        if (eos%ionE) then
            call load_tables_LTE("p2eint")
            call try_load_tables_LTE("gamma1")
            call try_load_tables_LTE("eint_from_T")
        end if
    end subroutine load_state_LTE

    !> Table-based method finalise: shift the loaded T/neOnH (and ionE-derived)
    !> tables to code units, build the derived tables not present on disk, then
    !> precompute step/guard arrays and validate. shift_table_to_code: axis2='eint'
    !> uses the log(p/nH) axis-2 unit, value_shift subtracted from table VALUES
    !> (T stores log T -> shift by log(unit_temperature); neOnH dimensionless -> 0);
    !> it also ensures node arrays exist so the build_*_table routines see a single
    !> var{1,2}_nodes path. Called from eos_finalise_LTE for eos_method='tables'.
    subroutine finalise_state_LTE()
        call shift_table_to_code(eos%T,     .false., dlog10(unit_temperature))
        call shift_table_to_code(eos%neOnH, .false., 0.d0)

        if (eos%ionE) then
            if (allocated(eos%p2eint%table)) then
                call shift_table_to_code(eos%p2eint, .false., 0.d0)
            end if

            !> Build derived tables from loaded T and neOnH. For any table that
            !> has a loaded version on disk (uniform or adaptive), KEEP the loaded
            !> one -- this enables the hybrid mode (T+neOnH uniform for the
            !> interleaved fast path; gamma1/p2eint/eint_from_T adaptive for the
            !> gamma_1 peak / tighter closure). Loaded versions came from the same
            !> Saha solver as T/neOnH so they are physically consistent. Skipping
            !> rebuild also avoids the one-time runtime cost but sacrifices accuracy
            if (.not. allocated(eos%gamma1%table)) then
                call build_gamma1_table()
            else
                if (mype == 0) write(*,'(A)') " Using loaded gamma1 table (skip rebuild)"
                !> gamma1 on (log_nH, log_eint/nH) CGS log10; axes -> code units,
                !> dimensionless values (no value shift).
                call shift_table_to_code(eos%gamma1, .false., 0.d0)
            end if

            !> Order matters: build_gamma1_p_table reads eos%p2eint%var*_nodes,
            !> so p2eint must be allocated first (loaded from disk above, or built
            !> here for the entropy-sourced case).
            if (.not. allocated(eos%p2eint%table)) then
                call build_p2eint_table()
            else
                if (mype == 0) write(*,'(A)') " Using loaded p2eint table (skip rebuild)"
                call ensure_axis_nodes(eos%p2eint)
            end if

            call build_gamma1_p_table()
            call build_log_p_table()
            call build_p_over_nH_table()

            if (.not. allocated(eos%eint_from_T%table)) then
                call build_eint_from_T_table()
            else
                if (mype == 0) write(*,'(A)') " Using loaded eint_from_T table (skip rebuild)"
                !> eint_from_T on (log_nH, log_T) CGS log10, values log10(eint/nH)
                !> CGS. axis2='T'; value shift log(p/nH).
                call shift_table_to_code(eos%eint_from_T, .true., &
                     dlog10(unit_pressure/unit_numberdensity))
            end if
        endif

        call precompute_step_inv(eos%T)
        call precompute_step_inv(eos%neOnH)
        call eos_build_guards(eos%T)
        call eos_build_guards(eos%neOnH)
        if (eos%ionE) then
            call precompute_step_inv(eos%p2eint)
            call precompute_step_inv(eos%gamma1)
            call precompute_step_inv(eos%gamma1_p)
            call precompute_step_inv(eos%eint_from_T)
            call precompute_step_inv(eos%log_p)
            call precompute_step_inv(eos%p_over_nH)
            call eos_build_guards(eos%p2eint)
            call eos_build_guards(eos%gamma1)
            call eos_build_guards(eos%gamma1_p)
            call eos_build_guards(eos%eint_from_T)
            call eos_build_guards(eos%log_p)
            call eos_build_guards(eos%p_over_nH)

            ! Build interleaved table: T, neOnH, p_over_nH
            call build_interleaved_eint_table()
        end if

        ! Defensive check: invariants on every populated table. Catches
        ! silent-failure modes (forgotten guard, mismatched node-array length,
        ! non-monotonic nodes) at startup rather than at first lookup.
        call eos_validate_table(eos%T,           "eos%T")
        call eos_validate_table(eos%neOnH,       "eos%neOnH")
        if (eos%ionE) then
            call eos_validate_table(eos%p2eint,      "eos%p2eint")
            call eos_validate_table(eos%gamma1,      "eos%gamma1")
            call eos_validate_table(eos%gamma1_p,    "eos%gamma1_p")
            call eos_validate_table(eos%eint_from_T, "eos%eint_from_T")
            call eos_validate_table(eos%log_p,       "eos%log_p")
            call eos_validate_table(eos%p_over_nH,   "eos%p_over_nH")
        end if

        call precompute_FI_bypass_constants()

        if (eos%ionE) call verify_eos_round_trips()
        !> Crucial to not end up with interpolation drift...
    end subroutine finalise_state_LTE

    !> 'tables' method: shift one loaded table from CGS log10 storage to code
    !> units and ensure its node arrays exist. axis2_kind selects the axis-2
    !> unit ('T' -> log(unit_temperature); else -> log(unit_pressure/
    !> unit_numberdensity)). value_shift is subtracted from the table VALUES
    !> (0 = dimensionless table, no value shift; e.g. unit_temperature for the
    !> T table, unit_pressure/unit_numberdensity for eint_from_T). Replaces the
    !> per-table copy-pasted shift idiom (the subtle per-table differences here
    !> were a past source of unit bugs).
    subroutine shift_table_to_code(tc, axis2_is_T, value_shift)
        type(eos_table_container), intent(inout) :: tc
        logical, intent(in) :: axis2_is_T 
        double precision, intent(in) :: value_shift
        double precision :: a1, a2
        a1 = dlog10(unit_numberdensity)
        if (axis2_is_T) then
            a2 = dlog10(unit_temperature)
        else
            a2 = dlog10(unit_pressure/unit_numberdensity)
        end if
        tc%var1_min = tc%var1_min - a1
        tc%var1_max = tc%var1_max - a1
        tc%var2_min = tc%var2_min - a2
        tc%var2_max = tc%var2_max - a2
        if (value_shift /= 0.d0) tc%table = tc%table - value_shift
        if (.not. tc%is_uniform) then
            tc%var1_nodes = tc%var1_nodes - a1
            tc%var2_nodes = tc%var2_nodes - a2
        end if
        call ensure_axis_nodes(tc)
    end subroutine shift_table_to_code

    !> Gamma_1 from (log nH, log eint/nH). Used only by build_gamma1_p_table,
    !> i.e. the 'tables' method; the entropy method reads its own g1p table.
    double precision function gamma1_from_nH_eint(log_nH, log_eint_nH) result(g1)
        double precision, intent(in) :: log_nH, log_eint_nH
        g1 = bicubic_lookup(log_nH, log_eint_nH, eos%gamma1)
    end function gamma1_from_nH_eint

    !> Build the first adiabatic index Gamma_1 table from T and neOnH tables.
    !> Called during eos_finalise after unit conversion, when tables are in code units.
    !>
    !> Gamma_1 = (rho/p) * (dp/drho)_s = a + b * (p / eint_vol)
    !>   where a = d(log10 p)/d(log10 nH) at constant eint/nH
    !>         b = d(log10 p)/d(log10 eint/nH) at constant nH
    !>   and p = nH * T * (1 + He + y), eint_vol = nH * (eint/nH)
    !>
    !> Derivation: from the thermodynamic identity
    !>   a^2 = (dp/drho)_{eint} + ((eint+p)/rho) * (dp/deint)_rho
    !> expressed in table coordinates (log10 nH, log10 eint/nH).
    !> For ideal gas: a=1, b=1, p/eint=gamma-1, so Gamma_1 = 1+(gamma-1) = gamma.
    subroutine build_gamma1_table()
        integer :: n1, n2, i, j, im, ip, jm, jp
        double precision :: x1, x2
        double precision, allocatable :: log_p(:,:)
        double precision :: T_val, y_val, p_val, eint_vol
        double precision :: a_val, b_val, g1_val
        double precision :: g1_min, g1_max

        n1 = eos%T%dim1
        n2 = eos%T%dim2

        !> Gamma1 table shares the same grid as T and neOnH (inherits
        !> uniform-or-adaptive node layout from the source).
        eos%gamma1%dim1 = n1
        eos%gamma1%dim2 = n2
        eos%gamma1%var1_min = eos%T%var1_min
        eos%gamma1%var1_max = eos%T%var1_max
        eos%gamma1%var2_min = eos%T%var2_min
        eos%gamma1%var2_max = eos%T%var2_max
        eos%gamma1%filename = 'computed_gamma1'
        eos%gamma1%is_uniform = eos%T%is_uniform
        if (allocated(eos%gamma1%var1_nodes)) deallocate(eos%gamma1%var1_nodes)
        if (allocated(eos%gamma1%var2_nodes)) deallocate(eos%gamma1%var2_nodes)
        allocate(eos%gamma1%var1_nodes(n1)); eos%gamma1%var1_nodes = eos%T%var1_nodes
        allocate(eos%gamma1%var2_nodes(n2)); eos%gamma1%var2_nodes = eos%T%var2_nodes

        if (allocated(eos%gamma1%table)) deallocate(eos%gamma1%table)
        allocate(eos%gamma1%table(n1, n2))
        allocate(log_p(n1, n2))

        !> Compute log10(p) at each grid point (using actual node positions)
        !> In code units: p = nH * T * (1 + He + y), kB absorbed
        do j = 1, n2
            x2 = eos%T%var2_nodes(j)
            do i = 1, n1
                x1 = eos%T%var1_nodes(i)
                T_val = 10.0d0**eos%T%table(i, j)         !> T in code units
                y_val = 10.0d0**eos%neOnH%table(i, j)     !> Ne/nH (dimensionless)
                !> log10(p) = log10(nH) + log10(T) + log10(1 + He + y)
                log_p(i, j) = x1 + eos%T%table(i, j) + dlog10(1.0d0 + eos%He_abundance + y_val)
            end do
        end do

        !> Compute Gamma_1 from finite differences of log_p using
        !> actual node spacings (correct for both uniform and adaptive grids).
        g1_min = 1.0d30
        g1_max = -1.0d30

        do j = 1, n2
            x2 = eos%T%var2_nodes(j)
            do i = 1, n1
                x1 = eos%T%var1_nodes(i)

                !> Centered finite differences with one-sided at boundaries
                im = max(i - 1, 1)
                ip = min(i + 1, n1)
                jm = max(j - 1, 1)
                jp = min(j + 1, n2)

                !> p/eint_vol = T * (1+He+y) / (eint/nH)
                T_val = 10.0d0**eos%T%table(i, j)
                y_val = 10.0d0**eos%neOnH%table(i, j)
                eint_vol = 10.0d0**x2  !> eint/nH in code units
                p_val = T_val * (1.0d0 + eos%He_abundance + y_val)

                if (eos%gamma1_method == 'effective') then
                    !> Approximate: gamma_eff = 1 + p/eint_vol
                    g1_val = 1.0d0 + p_val / eint_vol
                else
                    !> Full formula: Gamma_1 = a + b * (p / eint_vol)
                    !> a = d(log10 p) / d(log10 nH) at constant eint/nH
                    a_val = (log_p(ip, j) - log_p(im, j)) &
                          / (eos%T%var1_nodes(ip) - eos%T%var1_nodes(im))
                    !> b = d(log10 p) / d(log10 eint/nH) at constant nH
                    b_val = (log_p(i, jp) - log_p(i, jm)) &
                          / (eos%T%var2_nodes(jp) - eos%T%var2_nodes(jm))
                    g1_val = a_val + b_val * (p_val / eint_vol)
                endif

                !> Clamp to physical range [1, 5/3] (monatomic H+He gas)
                g1_val = max(1.001d0, min(g1_val, 5.0d0/3.0d0))

                eos%gamma1%table(i, j) = g1_val

                g1_min = min(g1_min, g1_val)
                g1_max = max(g1_max, g1_val)
            end do
        end do

        deallocate(log_p)

        if (mype == 0) then
            write(*, '(A,A,A,F8.4,A,F8.4)') &
                ' Gamma1 table built (', trim(eos%gamma1_method), '): min = ', &
                g1_min, ', max = ', g1_max
        end if

    end subroutine build_gamma1_table

    !> Build merged log10(p/nH) table: log10(p/nH)(nH, eint/nH).
    !> Same grid as T and neOnH tables. Each entry is:
    !>   log10(p/nH) = log10( 10^T * (1 + He + 10^y) )
    !> where T and y are from the forward tables at the same grid point.
    !> Used by the WB bisection to replace two separate PCHIP lookups
    !> (T + neOnH) with a single lookup per iteration.
    subroutine build_log_p_table()
        integer :: n1, n2, i, j
        double precision :: log_nH, log_eint_nH, T_val, y_val

        n1 = eos%T%dim1
        n2 = eos%T%dim2

        eos%log_p%dim1 = n1
        eos%log_p%dim2 = n2
        eos%log_p%var1_min = eos%T%var1_min
        eos%log_p%var1_max = eos%T%var1_max
        eos%log_p%var2_min = eos%T%var2_min
        eos%log_p%var2_max = eos%T%var2_max
        eos%log_p%is_uniform = eos%T%is_uniform
        if (allocated(eos%log_p%var1_nodes)) deallocate(eos%log_p%var1_nodes)
        if (allocated(eos%log_p%var2_nodes)) deallocate(eos%log_p%var2_nodes)
        allocate(eos%log_p%var1_nodes(n1)); eos%log_p%var1_nodes = eos%T%var1_nodes
        allocate(eos%log_p%var2_nodes(n2)); eos%log_p%var2_nodes = eos%T%var2_nodes

        if (allocated(eos%log_p%table)) deallocate(eos%log_p%table)
        allocate(eos%log_p%table(n1, n2))

        do j = 1, n2
            do i = 1, n1
                T_val = eos%T%table(i, j)
                y_val = eos%neOnH%table(i, j)
                eos%log_p%table(i, j) = dlog10(10.0d0**T_val &
                    * (1.0d0 + eos%He_abundance + 10.0d0**y_val))
            end do
        end do

        if (mype == 0) then
            write(*, '(A,I4,A,I4)') &
                ' Merged log_p table built: ', n1, ' x ', n2
        end if

    end subroutine build_log_p_table

    !> Build p/nH table: (1+He+y)*T at each (nH, eint/nH) grid point.
    !> Shares axes with the T and neOnH tables.
    !> Used by to_primitive_LTE for single-lookup pressure computation.
    subroutine build_p_over_nH_table()
        integer :: n1, n2, i, j
        double precision :: T_val, y_val, p_nH_val
        double precision :: p_min, p_max

        n1 = eos%T%dim1
        n2 = eos%T%dim2

        eos%p_over_nH%dim1 = n1
        eos%p_over_nH%dim2 = n2
        eos%p_over_nH%var1_min = eos%T%var1_min
        eos%p_over_nH%var1_max = eos%T%var1_max
        eos%p_over_nH%var2_min = eos%T%var2_min
        eos%p_over_nH%var2_max = eos%T%var2_max
        eos%p_over_nH%filename = 'computed_p_over_nH'
        eos%p_over_nH%is_uniform = eos%T%is_uniform
        if (allocated(eos%p_over_nH%var1_nodes)) deallocate(eos%p_over_nH%var1_nodes)
        if (allocated(eos%p_over_nH%var2_nodes)) deallocate(eos%p_over_nH%var2_nodes)
        allocate(eos%p_over_nH%var1_nodes(n1)); eos%p_over_nH%var1_nodes = eos%T%var1_nodes
        allocate(eos%p_over_nH%var2_nodes(n2)); eos%p_over_nH%var2_nodes = eos%T%var2_nodes

        if (allocated(eos%p_over_nH%table)) deallocate(eos%p_over_nH%table)
        allocate(eos%p_over_nH%table(n1, n2))

        p_min = 1.0d30
        p_max = -1.0d30

        do j = 1, n2
            do i = 1, n1
                !> T table stores log10(T_code), neOnH stores log10(y)
                T_val = 10.0d0**eos%T%table(i, j)
                y_val = 10.0d0**eos%neOnH%table(i, j)
                !> p/nH = (1 + He + y) * T in code units
                p_nH_val = (1.0d0 + eos%He_abundance + y_val) * T_val
                !> Store as log10 for PCHIP interpolation consistency
                eos%p_over_nH%table(i, j) = dlog10(p_nH_val)
                p_min = min(p_min, p_nH_val)
                p_max = max(p_max, p_nH_val)
            end do
        end do

        call precompute_step_inv(eos%p_over_nH)

        if (mype == 0) then
            write(*, '(A,ES10.3,A,ES10.3)') &
                ' p/nH table built: min = ', p_min, ', max = ', p_max
        end if

    end subroutine build_p_over_nH_table

    !> Build interleaved table from existing T, neOnH, p_over_nH tables.
    !> Layout: table_eint_il(3, ny, nx) where:
    !>   slot 1 = T (log10), slot 2 = neOnH (log10), slot 3 = p_over_nH (log10)
    !> All three are on the same (nH, eint/nH) axes.
    subroutine build_interleaved_eint_table()
        integer :: n1, n2, i, j

        n1 = eos%T%dim1
        n2 = eos%T%dim2

        if (allocated(eos%table_eint_il)) deallocate(eos%table_eint_il)
        allocate(eos%table_eint_il(3, n1, n2))

        do j = 1, n2
            do i = 1, n1
                eos%table_eint_il(1, i, j) = eos%T%table(i, j)
                eos%table_eint_il(2, i, j) = eos%neOnH%table(i, j)
                eos%table_eint_il(3, i, j) = eos%p_over_nH%table(i, j)
            end do
        end do

        if (mype == 0) then
            write(*, '(A,I0,A,I0,A,I0,A,F6.1,A)') &
                ' Interleaved table built: ', 3, ' x ', n1, ' x ', n2, &
                ' (', 3.0d0*n1*n2*8.0d0/1024.0d0, ' KB)'
        end if

    end subroutine build_interleaved_eint_table

    !> Build pressure-indexed Gamma_1 table: Gamma_1(nH, p/nH).
    !> Re-indexes the eint-based gamma1 table into pressure space,
    !> sharing the same grid axes as the p2eint table.
    !> This eliminates the intermediate p2eint lookup from hd_get_csound2_LTE.
    subroutine build_gamma1_p_table()
        integer :: n1, n2, i, j
        double precision :: log_nH, log_p_nH
        double precision :: p2eint_ratio, log_eint_nH
        double precision :: g1_val, g1_min, g1_max

        n1 = eos%p2eint%dim1
        n2 = eos%p2eint%dim2

        ! gamma1_p shares the same axis layout as p2eint (inherits adaptive
        ! v1 nodes from T via p2eint; v2 is uniform over the log_p/nH range).
        eos%gamma1_p%dim1 = n1
        eos%gamma1_p%dim2 = n2
        eos%gamma1_p%var1_min = eos%p2eint%var1_min
        eos%gamma1_p%var1_max = eos%p2eint%var1_max
        eos%gamma1_p%var2_min = eos%p2eint%var2_min
        eos%gamma1_p%var2_max = eos%p2eint%var2_max
        eos%gamma1_p%filename = 'computed_gamma1_p'
        eos%gamma1_p%is_uniform = eos%p2eint%is_uniform
        if (allocated(eos%gamma1_p%var1_nodes)) deallocate(eos%gamma1_p%var1_nodes)
        if (allocated(eos%gamma1_p%var2_nodes)) deallocate(eos%gamma1_p%var2_nodes)
        allocate(eos%gamma1_p%var1_nodes(n1)); eos%gamma1_p%var1_nodes = eos%p2eint%var1_nodes
        allocate(eos%gamma1_p%var2_nodes(n2)); eos%gamma1_p%var2_nodes = eos%p2eint%var2_nodes

        if (allocated(eos%gamma1_p%table)) deallocate(eos%gamma1_p%table)
        allocate(eos%gamma1_p%table(n1, n2))

        g1_min = 1.0d30
        g1_max = -1.0d30

        do j = 1, n2
            log_p_nH = eos%p2eint%var2_nodes(j)
            do i = 1, n1
                log_nH = eos%p2eint%var1_nodes(i)

                ! Convert from (nH, p/nH) to (nH, eint/nH) via the p2eint table
                ! p2eint table stores the ratio eint/(p), so:
                ! log10(eint/nH) = log10(p/nH) + log10(p2eint_ratio)
                p2eint_ratio = eos%p2eint%table(i, j)
                log_eint_nH = log_p_nH + dlog10(p2eint_ratio)

                ! Look up Gamma_1 from the eint-indexed table
                g1_val = gamma1_from_nH_eint(log_nH, log_eint_nH)

                eos%gamma1_p%table(i, j) = g1_val

                g1_min = min(g1_min, g1_val)
                g1_max = max(g1_max, g1_val)
            end do
        end do

        if (mype == 0) then
            write(*, '(A,F8.4,A,F8.4)') &
                ' Gamma1_p table built (pressure-indexed): min = ', g1_min, ', max = ', g1_max
        end if

    end subroutine build_gamma1_p_table

    !> Build p2eint table at runtime by inverting the forward (T, neOnH) tables.
    !> For each (nH, p/nH) grid point, find eint/nH by bisection such that
    !> the PCHIP-interpolated forward tables give p(nH, eint/nH) = p_target.
    !> This guarantees exact round-trip: to_conserved(p) -> eint -> to_primitive -> p.
    !>
    !> The p/nH grid bounds are computed from the forward tables (not hardcoded),
    !> ensuring coverage of reasonable pressures.
    subroutine build_p2eint_table()
        integer :: n1, n2_fwd, n2_p, i, j, iter
        double precision :: log_nH_i, log_eint_lo, log_eint_hi, log_eint_mid
        double precision :: log_p_target, log_p_eval
        double precision :: T_val, y_val
        double precision :: p_global_min, p_global_max
        double precision :: dx_p, log_p_nH_ij
        double precision :: max_err, err
        double precision :: log_p_lo_i, log_p_hi_i
        integer :: i_worst, j_worst

        if (mype == 0) write(*,*) 'Building p2eint table from forward tables...'

        n1 = eos%T%dim1
        n2_fwd = eos%T%dim2

        !> Compute p/nH at every forward grid node to find the pressure range
        p_global_min =  1.0d30
        p_global_max = -1.0d30
        do j = 1, n2_fwd
            do i = 1, n1
                T_val = 10.0d0**eos%T%table(i, j)
                y_val = 10.0d0**eos%neOnH%table(i, j)
                log_p_nH_ij = dlog10(T_val * (1.0d0 + eos%He_abundance + y_val))
                p_global_min = min(p_global_min, log_p_nH_ij)
                p_global_max = max(p_global_max, log_p_nH_ij)
            end do
        end do

        !> Pad pressure range to avoid edge clamping at off-grid PCHIP evaluations.
        !> The PCHIP polynomial can overshoot grid-node extrema,
        !> which for steep ionisation-zone gradients needs generous padding.
        p_global_min = p_global_min - 0.2d0
        p_global_max = p_global_max + 0.2d0

        n2_p = n2_fwd  !> Same resolution as forward table
        if (allocated(eos%p2eint%table)) deallocate(eos%p2eint%table)
        eos%p2eint%dim1 = n1
        eos%p2eint%dim2 = n2_p
        eos%p2eint%var1_min = eos%T%var1_min
        eos%p2eint%var1_max = eos%T%var1_max
        eos%p2eint%var2_min = p_global_min
        eos%p2eint%var2_max = p_global_max
        eos%p2eint%filename = 'computed_p2eint'
        !> v1 axis inherited from T (adaptive if T is); v2 axis is uniform
        !> over the (different) log_p/nH range. When T is adaptive, we still
        !> store explicit var2_nodes so lookups dispatch via _nu uniformly.
        eos%p2eint%is_uniform = eos%T%is_uniform
        if (allocated(eos%p2eint%var1_nodes)) deallocate(eos%p2eint%var1_nodes)
        if (allocated(eos%p2eint%var2_nodes)) deallocate(eos%p2eint%var2_nodes)
        allocate(eos%p2eint%var1_nodes(n1)); eos%p2eint%var1_nodes = eos%T%var1_nodes
        allocate(eos%p2eint%var2_nodes(n2_p))
        do j = 1, n2_p
            eos%p2eint%var2_nodes(j) = p_global_min &
                + (j - 1) * (p_global_max - p_global_min) / dble(n2_p - 1)
        end do
        allocate(eos%p2eint%table(n1, n2_p))

        dx_p = (p_global_max - p_global_min) / dble(n2_p - 1)
        max_err = 0.0d0

        !> For each (nH_i, p_j/nH), bisect on log10(eint/nH) to find the
        !> value where the PCHIP-interpolated forward tables give this pressure.
        i_worst = 1
        j_worst = 1
        do i = 1, n1
            log_nH_i = eos%T%var1_nodes(i)

            !> Compute achievable p range for this nH from forward table endpoints
            T_val = bicubic_lookup(log_nH_i, eos%T%var2_min, eos%T)
            y_val = bicubic_lookup(log_nH_i, eos%T%var2_min, eos%neOnH)
            log_p_lo_i = dlog10(10.0d0**T_val &
                * (1.0d0 + eos%He_abundance + 10.0d0**y_val))
            T_val = bicubic_lookup(log_nH_i, eos%T%var2_max, eos%T)
            y_val = bicubic_lookup(log_nH_i, eos%T%var2_max, eos%neOnH)
            log_p_hi_i = dlog10(10.0d0**T_val &
                * (1.0d0 + eos%He_abundance + 10.0d0**y_val))

            do j = 1, n2_p
                log_p_target = p_global_min + (j - 1) * dx_p

                !> Clamp target to achievable range -- outside this range,
                !> map to the table boundary, should be fine
                if (log_p_target <= log_p_lo_i) then
                    eos%p2eint%table(i, j) = 10.0d0**(eos%T%var2_min - log_p_target)
                    cycle
                else if (log_p_target >= log_p_hi_i) then
                    eos%p2eint%table(i, j) = 10.0d0**(eos%T%var2_max - log_p_target)
                    cycle
                end if

                !> Bisection on log10(eint/nH)
                log_eint_lo = eos%T%var2_min
                log_eint_hi = eos%T%var2_max

                do iter = 1, 52  !> 52 bisections -> 2^-52 ~ 10^-15.7 precision
                    log_eint_mid = 0.5d0 * (log_eint_lo + log_eint_hi)

                    !> Evaluate p from forward tables at (nH_i, eint_mid)
                    T_val = bicubic_lookup(log_nH_i, log_eint_mid, eos%T)
                    y_val = bicubic_lookup(log_nH_i, log_eint_mid, eos%neOnH)
                    log_p_eval = dlog10(10.0d0**T_val &
                        * (1.0d0 + eos%He_abundance + 10.0d0**y_val))

                    if (log_p_eval < log_p_target) then
                        log_eint_lo = log_eint_mid
                    else
                        log_eint_hi = log_eint_mid
                    end if

                    if (dabs(log_eint_hi - log_eint_lo) < 1.0d-14) exit
                end do

                eos%p2eint%table(i, j) = 10.0d0**(log_eint_mid - log_p_target)

                !> Track round-trip error (only for in-range points)
                err = dabs(log_p_eval - log_p_target)
                if (err > max_err) then
                    max_err = err
                    i_worst = i
                    j_worst = j
                end if
            end do
        end do

        if (mype == 0) then
            write(*, '(A,ES10.3,A,I4,A,I4)') &
                ' p2eint table built: max bisection err = ', max_err, &
                ' at (i,j)=(', i_worst, ',', j_worst, ')'
            write(*, '(A,F8.4,A,F8.4)') &
                '   log10(p/nH) range = ', p_global_min, ' to ', p_global_max
        end if

    end subroutine build_p2eint_table

    !> Build inverse T table: given (nH, T), return log10(eint/nH).
    !> Inverts the forward T(nH, eint/nH) table by bisection using the SAME
    !> PCHIP interpolation kernel as T_from_nH_eint. This guarantees that
    !> T_from_nH_eint(nH, eint_from_T(nH, T)) = T to machine precision.
    !> Called during eos_finalise after unit conversion.
    subroutine build_eint_from_T_table()
        integer :: n1, n2_fwd, n2_inv, i, j, iter
        double precision :: dx2_inv
        double precision :: T_global_min, T_global_max, T_FI_threshold
        double precision :: log_nH_i, log_eint_lo, log_eint_hi, log_eint_mid
        double precision :: T_target, T_eval, T_max_at_nH, eint_nH_FI
        double precision :: eint_nH_min, eint_nH_max, max_err, err
        integer :: n_FI_filled

        n1 = eos%T%dim1
        n2_fwd = eos%T%dim2

        !> Find global T range from the forward T table
        T_global_min =  1.0d30
        T_global_max = -1.0d30
        do j = 1, n2_fwd
            do i = 1, n1
                T_global_min = min(T_global_min, eos%T%table(i, j))
                T_global_max = max(T_global_max, eos%T%table(i, j))
            end do
        end do

        !> FI threshold temperature in code units (200,000 K)
        T_FI_threshold = eos%p_rho_FI_threshold &
            * (1.0d0 + 4.0d0*eos%He_abundance) / eos%n_per_nH_FI

        !> Extend T axis to cover the FI regime up to 10 MK (or T_global_max if larger)
        T_global_max = max(T_global_max, dlog10(1.0d7/unit_temperature))

        !> Set up inverse table with same nH axis, extended T axis
        n2_inv = n2_fwd  !> Same resolution as forward table
        if (allocated(eos%eint_from_T%table)) deallocate(eos%eint_from_T%table)
        eos%eint_from_T%dim1 = n1
        eos%eint_from_T%dim2 = n2_inv
        eos%eint_from_T%var1_min = eos%T%var1_min
        eos%eint_from_T%var1_max = eos%T%var1_max
        eos%eint_from_T%var2_min = T_global_min
        eos%eint_from_T%var2_max = T_global_max
        eos%eint_from_T%filename = 'computed_eint_from_T'
        !> Inherit v1 axis from T (adaptive if T is); uniform v2 axis covers
        !> [T_global_min, T_global_max] with explicit nodes.
        eos%eint_from_T%is_uniform = eos%T%is_uniform
        if (allocated(eos%eint_from_T%var1_nodes)) deallocate(eos%eint_from_T%var1_nodes)
        if (allocated(eos%eint_from_T%var2_nodes)) deallocate(eos%eint_from_T%var2_nodes)
        allocate(eos%eint_from_T%var1_nodes(n1)); eos%eint_from_T%var1_nodes = eos%T%var1_nodes
        allocate(eos%eint_from_T%var2_nodes(n2_inv))
        do j = 1, n2_inv
            eos%eint_from_T%var2_nodes(j) = T_global_min &
                + (j - 1) * (T_global_max - T_global_min) / dble(n2_inv - 1)
        end do

        allocate(eos%eint_from_T%table(n1, n2_inv))

        dx2_inv = (T_global_max - T_global_min) / dble(n2_inv - 1)

        eint_nH_min =  1.0d30
        eint_nH_max = -1.0d30
        max_err = 0.0d0
        n_FI_filled = 0

        !> For each (nH_i, T_j): use bisection on the forward table where it
        !> has coverage, and the analytic FI formula above the table range.
        do i = 1, n1
            log_nH_i = eos%T%var1_nodes(i)

            !> Find the maximum T in the forward table at this nH
            T_max_at_nH = eos%T%table(i, n2_fwd)

            do j = 1, n2_inv
                T_target = T_global_min + (j - 1) * dx2_inv

                if (10.0d0**T_target > T_FI_threshold .or. &
                    T_target > T_max_at_nH) then
                    !> Above table coverage or FI regime: use analytic formula
                    !> eint/nH = n_per_nH_FI / (gamma-1) * T + neOnH_FI * eion_per_nH
                    eint_nH_FI = eos%n_per_nH_FI * eos%inv_gamma_minus_1 &
                        * 10.0d0**T_target
                    if (eos%ionE) eint_nH_FI = eint_nH_FI &
                        + eos%neOnH_FI * eos%eion_per_nH
                    eos%eint_from_T%table(i, j) = dlog10(eint_nH_FI)
                    n_FI_filled = n_FI_filled + 1
                else
                    !> Bisection on log10(eint/nH) using the forward T table
                    log_eint_lo = eos%T%var2_min
                    log_eint_hi = eos%T%var2_max

                    do iter = 1, 52
                        log_eint_mid = 0.5d0 * (log_eint_lo + log_eint_hi)

                        T_eval = bicubic_lookup(log_nH_i, log_eint_mid, eos%T)

                        if (T_eval < T_target) then
                            log_eint_lo = log_eint_mid
                        else
                            log_eint_hi = log_eint_mid
                        end if

                        if (dabs(log_eint_hi - log_eint_lo) < 1.0d-14) exit
                    end do

                    eos%eint_from_T%table(i, j) = log_eint_mid

                    err = dabs(T_eval - T_target)
                    max_err = max(max_err, err)
                end if

                eint_nH_min = min(eint_nH_min, eos%eint_from_T%table(i, j))
                eint_nH_max = max(eint_nH_max, eos%eint_from_T%table(i, j))
            end do
        end do

        if (mype == 0) then
            write(*, '(A,ES10.3,A,I0,A,I0)') &
                ' Inverse T table built: max bisection err = ', max_err, &
                ', FI-filled = ', n_FI_filled, ' / ', n1*n2_inv
            write(*, '(A,F8.4,A,F8.4)') &
                '   log10(T) range = ', T_global_min, ' to ', T_global_max
            write(*, '(A,F8.4,A,F8.4)') &
                '   log10(eint/nH) range = ', eint_nH_min, ' to ', eint_nH_max
        end if

    end subroutine build_eint_from_T_table

    !> Verify round-trip consistency of ALL EoS table pathways.
    !> Tests at off-grid random points to measure real interpolation error.
    !> Called after all tables are built.
    subroutine verify_eos_round_trips()
        integer :: n1, n2, i, j, n_tested
        double precision :: log_nH, log_eint_nH
        double precision :: dx1, dx2
        double precision :: T_val, y_val, log_p_nH, p2eint_val, log_eint_recovered
        double precision :: T_recovered, p_from_Ty, log_p_recovered
        double precision :: eint_from_T_val, log_eint_from_T_recovered
        double precision :: err_p, err_eint_T
        double precision :: max_err_p, max_err_eint_T, mean_err_p, mean_err_eint_T
        double precision :: log_nH_worst_p, log_eint_worst_p

        n1 = eos%T%dim1
        n2 = eos%T%dim2
        dx1 = (eos%T%var1_max - eos%T%var1_min) / dble(n1 - 1)
        dx2 = (eos%T%var2_max - eos%T%var2_min) / dble(n2 - 1)

        max_err_p = 0.0d0
        max_err_eint_T = 0.0d0
        mean_err_p = 0.0d0
        mean_err_eint_T = 0.0d0
        n_tested = 0

        !> Test at cell-centre-like points (offset by 0.5*dx from grid nodes)
        !> to exercise the PCHIP interpolation between nodes.
        do i = 2, n1 - 1
            log_nH = eos%T%var1_min + (dble(i) - 0.5d0) * dx1
            do j = 2, n2 - 1
                log_eint_nH = eos%T%var2_min + (dble(j) - 0.5d0) * dx2
                n_tested = n_tested + 1

                !> === Round-trip 1: eint -> T,y -> p -> p2eint -> eint' ===
                !> Forward: get T and y from eint
                T_val = bicubic_lookup(log_nH, log_eint_nH, eos%T)
                y_val = bicubic_lookup(log_nH, log_eint_nH, eos%neOnH)
                !> Compute p/nH
                log_p_nH = dlog10(10.0d0**T_val &
                    * (1.0d0 + eos%He_abundance + 10.0d0**y_val))
                !> Inverse: get eint from p via p2eint table
                p2eint_val = bicubic_lookup(log_nH, log_p_nH, eos%p2eint)
                log_eint_recovered = log_p_nH + dlog10(p2eint_val)

                err_p = dabs(log_eint_recovered - log_eint_nH)
                mean_err_p = mean_err_p + err_p
                if (err_p > max_err_p) then
                    max_err_p = err_p
                    log_nH_worst_p = log_nH
                    log_eint_worst_p = log_eint_nH
                end if

                !> === Round-trip 2: eint -> T -> eint_from_T -> eint' ===
                eint_from_T_val = bicubic_lookup(log_nH, T_val, eos%eint_from_T)

                err_eint_T = dabs(eint_from_T_val - log_eint_nH)
                mean_err_eint_T = mean_err_eint_T + err_eint_T
            end do
        end do

        mean_err_p = mean_err_p / dble(max(n_tested, 1))
        mean_err_eint_T = mean_err_eint_T / dble(max(n_tested, 1))

        if (mype == 0) then
            write(*, '(A)') ' === EoS round-trip verification (off-grid points) ==='
            write(*, '(A,I6,A)') '   Tested ', n_tested, ' points (half-cell offsets)'
            write(*, '(A,ES10.3,A,ES10.3)') &
                '   p round-trip:      max err = ', max_err_p, &
                '  mean err = ', mean_err_p
            write(*, '(A,F8.4,A,F8.4)') &
                '     worst at log_nH = ', log_nH_worst_p, &
                ', log_eint = ', log_eint_worst_p
            write(*, '(A,ES10.3,A,ES10.3)') &
                '   T round-trip:      max err = ', max_err_eint_T, &
                '  mean err = ', mean_err_eint_T
            write(*, '(A)') ' ====================================================='
        end if

    end subroutine verify_eos_round_trips

end module mod_eos_lte_state
!> Needs a line after to pass the preprocessor
