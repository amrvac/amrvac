!=============================================================================
!> Shared LTE lookup-table infrastructure (build/IO; not on the hot path).
!>
!> The table-file machinery used by ALL LTE methods (state, entropy, analytic):
!> binary-file readers (read_eos_from_file / load_tables_LTE / try_load_tables_LTE),
!> per-container preparation (ensure_axis_nodes, code-unit axis shifts, O(1) guard
!> arrays via eos_build_guards, eos_validate_table) and the fully-ionised bypass
!> constants. NOT a method: 'tables' here means the literal lookup tables, not the
!> 'state' method (which lives in mod_eos_LTE_state and consumes this infra).
!> Everything here runs once, during eos_init_LTE / eos_finalise_LTE.
!=============================================================================
module mod_eos_LTE_tables
    use mod_global_parameters
    use mod_eos_container
    use mod_eos_interp
    use mod_comm_lib, only: mpistop
    implicit none
    private

    !> Shared table-file I/O + per-container prep, consumed by mod_eos_LTE_state
    !> (state method) and mod_eos_LTE_entropy (entropy method).
    public :: load_tables_LTE, try_load_tables_LTE
    public :: ensure_axis_nodes, shift_axis_to_code, shift_axis_to_code_T
    public :: eos_build_guards, eos_validate_table, precompute_FI_bypass_constants

    character(len=std_len), parameter :: eos_table_prefix = "LTEeos_"

contains

    !> Ensure tc%var1_nodes / tc%var2_nodes are allocated and populated so
    !> that build_*_table routines and other code can reference axis positions
    !> uniformly (no is_uniform branching). For uniform tables we generate
    !> evenly-spaced node arrays from var{1,2}_min/max.
    !>
    !> Lookup performance is unaffected: bicubic_lookup / bilinear_lookup
    !> still branch on tc%is_uniform -- uniform tables go through the affine
    !> fast path; only adaptive tables consult var{1,2}_nodes.
    subroutine ensure_axis_nodes(tc)
        type(eos_table_container), intent(inout) :: tc
        integer :: k
        double precision :: dx
        if (.not. allocated(tc%table)) return
        if (.not. allocated(tc%var1_nodes)) then
            allocate(tc%var1_nodes(tc%dim1))
            dx = (tc%var1_max - tc%var1_min) / dble(tc%dim1 - 1)
            do k = 1, tc%dim1
                tc%var1_nodes(k) = tc%var1_min + (k - 1) * dx
            end do
        end if
        if (.not. allocated(tc%var2_nodes)) then
            allocate(tc%var2_nodes(tc%dim2))
            dx = (tc%var2_max - tc%var2_min) / dble(tc%dim2 - 1)
            do k = 1, tc%dim2
                tc%var2_nodes(k) = tc%var2_min + (k - 1) * dx
            end do
        end if
    end subroutine ensure_axis_nodes

    !> Shift a table's (log_nH, log_eint/nH) axes -- and adaptive node arrays
    !> if present -- from CGS storage to code units.  Used by the entropy
    !> method so eos_finalise can apply the shift AFTER unit_* are populated.
    subroutine shift_axis_to_code(tc)
        !> CGS -> code shift assuming axis-1 = log(nH) and axis-2 = log(eint/nH)
        !> (or log(p/nH); both share the same erg/particle unit conversion).
        type(eos_table_container), intent(inout) :: tc
        tc%var1_min = tc%var1_min - dlog10(unit_numberdensity)
        tc%var1_max = tc%var1_max - dlog10(unit_numberdensity)
        tc%var2_min = tc%var2_min - dlog10(unit_pressure/unit_numberdensity)
        tc%var2_max = tc%var2_max - dlog10(unit_pressure/unit_numberdensity)
        if (.not. tc%is_uniform) then
            tc%var1_nodes = tc%var1_nodes - dlog10(unit_numberdensity)
            tc%var2_nodes = tc%var2_nodes - dlog10(unit_pressure/unit_numberdensity)
        end if
    end subroutine shift_axis_to_code

    subroutine shift_axis_to_code_T(tc)
        !> CGS -> code shift for tables whose axis-2 = log(T) (K).
        !> Axis-1 = log(nH) shift as in the standard case.
        type(eos_table_container), intent(inout) :: tc
        tc%var1_min = tc%var1_min - dlog10(unit_numberdensity)
        tc%var1_max = tc%var1_max - dlog10(unit_numberdensity)
        tc%var2_min = tc%var2_min - dlog10(unit_temperature)
        tc%var2_max = tc%var2_max - dlog10(unit_temperature)
        if (.not. tc%is_uniform) then
            tc%var1_nodes = tc%var1_nodes - dlog10(unit_numberdensity)
            tc%var2_nodes = tc%var2_nodes - dlog10(unit_temperature)
        end if
    end subroutine shift_axis_to_code_T

    !> Build guard (bucket) arrays so that adaptive index lookup is O(1).
    !> No-op for uniform tables. Safe to call repeatedly.
    !> Idea here is to map nonuniform table nodes to a uniform search grid
    subroutine eos_build_guards(tc)
        type(eos_table_container), intent(inout) :: tc
        if (tc%is_uniform) return
        if (.not. allocated(tc%var1_nodes)) return
        if (.not. allocated(tc%var2_nodes)) return
        call build_guard_for_axis(tc%var1_nodes, tc%dim1, &
            tc%guard_1, tc%guard_M_1, tc%guard_scale_1)
        call build_guard_for_axis(tc%var2_nodes, tc%dim2, &
            tc%guard_2, tc%guard_M_2, tc%guard_scale_2)
    end subroutine eos_build_guards

    subroutine build_guard_for_axis(nodes, n, guard, M_out, scale_out)
        double precision, intent(in) :: nodes(:)
        integer, intent(in) :: n
        integer, allocatable, intent(inout) :: guard(:)
        integer, intent(out) :: M_out
        double precision, intent(out) :: scale_out

        integer, parameter :: SAFETY = 4
        integer :: k, ii, M
        double precision :: xmin, xmax, span, min_gap, bk

        xmin = nodes(1); xmax = nodes(n)
        span = xmax - xmin
        if (span <= 0.0d0 .or. n < 2) then
            M_out = 0; scale_out = 0.0d0
            if (allocated(guard)) deallocate(guard)
            return
        end if
        min_gap = huge(1.0d0)
        do ii = 1, n-1
            if (nodes(ii+1) - nodes(ii) < min_gap) min_gap = nodes(ii+1) - nodes(ii)
        end do
        if (min_gap <= 0.0d0) then
            !> Degenerate node array -- fall back to bsearch by leaving M=0
            M_out = 0; scale_out = 0.0d0
            if (allocated(guard)) deallocate(guard)
            return
        end if
        M = max(2*n, SAFETY * int(ceiling(span / min_gap)))
        if (allocated(guard)) deallocate(guard)
        allocate(guard(M))
        scale_out = dble(M) / span
        M_out = M

        ii = 1
        do k = 1, M
            bk = xmin + dble(k-1) * span / dble(M)
            do while (ii < n .and. nodes(ii) < bk)
                ii = ii + 1
            end do
            guard(k) = ii
        end do
    end subroutine build_guard_for_axis

    !> Defensive consistency check for one table after init. Aborts with a
    !> diagnostic message if any silent-failure invariant is violated.
    subroutine eos_validate_table(tc, name)
        type(eos_table_container), intent(in) :: tc
        character(len=*), intent(in) :: name
        integer :: i

        if (.not. allocated(tc%table)) then
            if (mype == 0) write(*,'(A,A)') " EOS validate FAIL: table not allocated for ", trim(name)
            call mpistop("eos_validate_table: missing table data")
        end if
        if (size(tc%table, 1) /= tc%dim1 .or. size(tc%table, 2) /= tc%dim2) then
            if (mype == 0) write(*,'(A,A)') " EOS validate FAIL: dim mismatch for ", trim(name)
            call mpistop("eos_validate_table: dim mismatch")
        end if
        !> Value finiteness: a NaN (e.g. from dlog10 of a non-positive stored
        !> value, or a truncated/corrupt file read) otherwise passes the
        !> structural checks and only surfaces as a runtime FPE deep in a lookup.
        !> The x /= x idiom is NaN-true without needing the ieee_arithmetic module.
        if (any(tc%table /= tc%table)) then
            if (mype == 0) write(*,'(A,A)') " EOS validate FAIL: NaN in table ", trim(name)
            call mpistop("eos_validate_table: NaN in table data")
        end if
        if (.not. tc%is_uniform) then
            if (.not. allocated(tc%var1_nodes) .or. .not. allocated(tc%var2_nodes)) then
                if (mype == 0) write(*,'(A,A)') &
                    " EOS validate FAIL: adaptive table missing node arrays for ", trim(name)
                call mpistop("eos_validate_table: adaptive table without nodes")
            end if
            if (size(tc%var1_nodes) /= tc%dim1 .or. size(tc%var2_nodes) /= tc%dim2) then
                if (mype == 0) write(*,'(A,A)') &
                    " EOS validate FAIL: node array length mismatch for ", trim(name)
                call mpistop("eos_validate_table: node length mismatch")
            end if
            do i = 1, tc%dim1 - 1
                if (tc%var1_nodes(i+1) <= tc%var1_nodes(i)) then
                    if (mype == 0) write(*,'(A,A,A,I0)') &
                        " EOS validate FAIL: var1_nodes not strictly monotonic for ", &
                        trim(name), " at index ", i
                    call mpistop("eos_validate_table: non-monotonic nodes")
                end if
            end do
            do i = 1, tc%dim2 - 1
                if (tc%var2_nodes(i+1) <= tc%var2_nodes(i)) then
                    if (mype == 0) write(*,'(A,A,A,I0)') &
                        " EOS validate FAIL: var2_nodes not strictly monotonic for ", &
                        trim(name), " at index ", i
                    call mpistop("eos_validate_table: non-monotonic nodes")
                end if
            end do
            if (.not. allocated(tc%guard_1) .or. .not. allocated(tc%guard_2) .or. &
                tc%guard_M_1 <= 0 .or. tc%guard_M_2 <= 0 .or. &
                tc%guard_scale_1 <= 0.0d0 .or. tc%guard_scale_2 <= 0.0d0) then
                if (mype == 0) write(*,'(A,A,A)') &
                    " EOS validate WARN: adaptive table without guard for ", trim(name), &
                    " - lookup will fall back to bsearch"
            end if
        end if
        if (mype == 0) write(*,'(A,A)') " EOS validate OK: ", trim(name)
    end subroutine eos_validate_table

    subroutine precompute_FI_bypass_constants()
        !> Precompute constants for bypassing table lookups when gas is fully ionised.
        !> Physics: above T ~ 50 kK, H is >99.99% ionised and He is >90% ionised.
        !> The EoS reduces to ideal gas with constant mu, Gamma_1 = 5/3, ne/nH = 1+2*A_He.
        !> We check eint/rho > threshold (1 division + 1 comparison) to skip all table lookups.
        double precision :: chi_H, chi_HeI, chi_HeII
        double precision :: T_thresh, eint_nH_thresh, p_nH_thresh

        !> Ionisation energies in CGS (erg)
        chi_H   = 13.6d0  * 1.602d-12
        chi_HeI = 24.6d0  * 1.602d-12
        chi_HeII = 54.4d0 * 1.602d-12

        !> Fully-ionised particle counts and electron fraction
        eos%n_per_nH_FI = 2.0d0 + 3.0d0 * eos%He_abundance
        eos%neOnH_FI    = 1.0d0 + 2.0d0 * eos%He_abundance

        !> Total ionisation energy per hydrogen atom (CGS), converted to code units
        !> eion_per_nH [code] = eion_per_nH [erg] / (unit_pressure / unit_numberdensity)
        eos%eion_per_nH = (chi_H + eos%He_abundance * (chi_HeI + chi_HeII)) &
                        * unit_numberdensity / unit_pressure
        !> No-energy mode (.not. ionE): the ionisation energy is NOT folded into
        !> eint, so it must not appear anywhere in the thermodynamics. Zeroing it
        !> here collapses both the FI threshold below and every FI-zone temperature
        !> path (update_eos_LTE, get_temperature_from_eint[_fast]_LTE, the analytic
        !> FI bypasses) to their correct no-energy form from one place -- the
        !> from-eint getters subtract eion unconditionally, so this is what keeps
        !> them consistent with update_eos_LTE's ionE-gated branch.
        if (.not. eos%ionE) eos%eion_per_nH = 0.0d0

        !> Bypass threshold: eint/rho value corresponding to T = 200,000 K fully ionised.
        !> Raised from 50 kK to 200 kK so that blended LTE tables (which converge to
        !> FI by 200 kK) match the bypass formula exactly at the threshold -- no
        !> discontinuity in ne/nH, T, Gamma_1, or p when cells cross the boundary.
        !> eint/nH = 1.5 * n_per_nH * kB * T + eion_per_nH  (all in code units)
        !> eint/rho = eint/nH / nH2rhoFactor
        T_thresh = 200000.0d0 / unit_temperature
        eint_nH_thresh = 1.5d0 * eos%n_per_nH_FI * T_thresh + eos%eion_per_nH
        eos%eint_rho_FI_threshold = eint_nH_thresh / eos%nH2rhoFactor

        !> Pressure-based threshold for primitive-variable checks:
        !> For FI gas, T = p / (rho * Rfactor_FI), so p/rho > Rfactor_FI * T_thresh
        !> where Rfactor_FI = n_per_nH_FI / (1 + 4*A_He) in code units
        eos%p_rho_FI_threshold = eos%n_per_nH_FI &
            / (1.0d0 + 4.0d0*eos%He_abundance) * T_thresh

        if (mype == 0) then
            write(*,'(a,es12.4)') ' FI bypass: eion_per_nH (code) = ', eos%eion_per_nH
            write(*,'(a,es12.4)') ' FI bypass: eint/rho threshold = ', eos%eint_rho_FI_threshold
            write(*,'(a,f8.4)')   ' FI bypass: n_per_nH_FI        = ', eos%n_per_nH_FI
            write(*,'(a,f8.4)')   ' FI bypass: neOnH_FI           = ', eos%neOnH_FI
        end if

    end subroutine precompute_FI_bypass_constants

    !> Read one named table file into its eos% container. The filename encodes
    !> the composition (H or HHe) and whether ionisation energy is included.
    subroutine load_tables_LTE(fieldname)
        character(len=*), intent(in) :: fieldname
        character(len=std_len) :: subname, tablename, filename, filepath
        subname = "H" !> Always consider at least hydrogen

        if (eos%He_abundance > 0.0d0) then
            write(subname,"(A,A)") trim(subname), "He"
        endif

        if (eos%ionE) then !> Consider the ionisation energy?
            write(subname,"(A,A)") trim(subname), "_IonE"
        else
            write(subname,"(A,A)") trim(subname), "_NoIonE"
        endif

        write(tablename, "(A,A,A)") trim(eos_table_prefix), fieldname, "_"
        write(filename,"(A,A,A)") trim(tablename), trim(subname), '.bin'

        if (mype==0) then
            print*, "Reading EoS tables from: ", trim(filename)
        endif

        write(filepath,"(A,A)") trim(eos%table_location), trim(filename)

        select case (fieldname)
        case("T")
            eos%T%filename = filename
            call read_eos_from_file(trim(filepath), eos%T, .true.)
        case("neOnH")
            eos%neOnH%filename = filename
            ! Legacy tables (method='tables'/'analytic') store log10(neOnH); the
            ! entropy method's tables store raw (linear) neOnH because they go
            ! through bicubic Hermite (which interpolates the value directly,
            ! no log10 transform needed and would break stored derivatives).
            call read_eos_from_file(trim(filepath), eos%neOnH, &
                                     eos%method_id /= EOS_ENTROPY)
        case("p2eint")
            eos%p2eint%filename = filename
            call read_eos_from_file(trim(filepath), eos%p2eint, .false.)
        case("gamma1")
            eos%gamma1%filename = filename
            call read_eos_from_file(trim(filepath), eos%gamma1, .false.)
        case("eint_from_T")
            eos%eint_from_T%filename = filename
            call read_eos_from_file(trim(filepath), eos%eint_from_T, .false.)
            !> Entropy-method tables. Unit shifts MUST happen in eos_finalise (where
            !> unit_* are set), not here.
        case("neOnH_x")
            eos%neOnH_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%neOnH_x, .false.)
        case("neOnH_y")
            eos%neOnH_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%neOnH_y, .false.)
        case("neOnH_xy")
            eos%neOnH_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%neOnH_xy, .false.)
        case("eintP")
            eos%eintP%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintP, .false.)
        case("eintP_x")
            eos%eintP_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintP_x, .false.)
        case("eintP_y")
            eos%eintP_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintP_y, .false.)
        case("eintP_xy")
            eos%eintP_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintP_xy, .false.)
        case("g1p")
            eos%g1p%filename = filename
            call read_eos_from_file(trim(filepath), eos%g1p, .false.)
        case("g1p_x")
            eos%g1p_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%g1p_x, .false.)
        case("g1p_y")
            eos%g1p_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%g1p_y, .false.)
        case("g1p_xy")
            eos%g1p_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%g1p_xy, .false.)
        case("eintT")
            eos%eintT%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintT, .false.)
        case("eintT_x")
            eos%eintT_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintT_x, .false.)
        case("eintT_y")
            eos%eintT_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintT_y, .false.)
        case("eintT_xy")
            eos%eintT_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%eintT_xy, .false.)
            ! 'yT' tables dropped -- zero runtime references.
        case("Tfwd")
            eos%Tfwd%filename = filename
            call read_eos_from_file(trim(filepath), eos%Tfwd, .false.)
        case("Tfwd_x")
            eos%Tfwd_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%Tfwd_x, .false.)
        case("Tfwd_y")
            eos%Tfwd_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%Tfwd_y, .false.)
        case("Tfwd_xy")
            eos%Tfwd_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%Tfwd_xy, .false.)
        case("pfwd")
            eos%pfwd%filename = filename
            call read_eos_from_file(trim(filepath), eos%pfwd, .false.)
        case("pfwd_x")
            eos%pfwd_x%filename = filename
            call read_eos_from_file(trim(filepath), eos%pfwd_x, .false.)
        case("pfwd_y")
            eos%pfwd_y%filename = filename
            call read_eos_from_file(trim(filepath), eos%pfwd_y, .false.)
        case("pfwd_xy")
            eos%pfwd_xy%filename = filename
            call read_eos_from_file(trim(filepath), eos%pfwd_xy, .false.)
        case default
            call mpistop("eos table name "//trim(fieldname)//" not recognised in load_tables_LTE")
        end select

    end subroutine load_tables_LTE

    !> Attempt to load pre-computed table from binary file.
    !> If the file does not exist, silently skip and the table will be built at runtime.
    subroutine try_load_tables_LTE(fieldname)
        character(len=*), intent(in) :: fieldname
        character(len=std_len) :: subname, tablename, filename, filepath
        logical :: file_exists

        subname = "H"
        if (eos%He_abundance > 0.0d0) then
            write(subname,"(A,A)") trim(subname), "He"
        endif
        if (eos%ionE) then
            write(subname,"(A,A)") trim(subname), "_IonE"
        else
            write(subname,"(A,A)") trim(subname), "_NoIonE"
        endif
        write(tablename, "(A,A,A)") trim(eos_table_prefix), fieldname, "_"
        write(filename,"(A,A,A)") trim(tablename), trim(subname), '.bin'
        write(filepath,"(A,A)") trim(eos%table_location), trim(filename)

        inquire(file=trim(filepath), exist=file_exists)
        if (file_exists) then
            call load_tables_LTE(fieldname)
        else
            if (mype == 0) write(*,*) &
                trim(fieldname)//' table not found, will build at runtime: '//trim(filename)
        endif
    end subroutine try_load_tables_LTE

    !> Read an EoS table from a binary file written by generate_state_tables.py.
    !>
    !> File format (uniform -- legacy):
    !>   [2 x int32]   : (dim1, dim2) = (dimy, dimx)   in Fortran column-major
    !>   [4 x float64] : (var1_min, var1_max, var2_min, var2_max)
    !>   [dim1*dim2 x float64] : table values
    !>
    !> File format (adaptive -- extension): same prefix, then optional trailer
    !>   [1 x int32]      : trailer_flag (= 1)
    !>   [dim1 x float64] : x_nodes (axis-1 node positions, log10 space)
    !>   [dim2 x float64] : y_nodes (axis-2 node positions, log10 space)
    !>
    !> The trailer is detected via end-of-file on the trailer_flag read:
    !> uniform tables (no trailer) leave is_uniform = .true. and produce
    !> ios /= 0 on the read; adaptive tables set is_uniform = .false. and
    !> populate x_nodes, y_nodes.
    subroutine read_eos_from_file(filename, tc, logtable)
        character(len=*), intent(in) :: filename
        type(eos_table_container), intent(inout) :: tc
        logical, intent(in) :: logtable

        integer :: ios, trailer_flag, lun

        !> The mandatory header/payload reads are iostat-guarded so a missing,
        !> truncated, or wrong-format file fails with a clear mpistop rather than
        !> an opaque "Fortran runtime error: end of file". (Only the OPTIONAL
        !> adaptive trailer below is allowed to hit EOF cleanly.)
        open(newunit=lun, file=filename, access='stream', form='unformatted', &
             action='read', iostat=ios)
        if (ios /= 0) call mpistop( &
            "read_eos_from_file: cannot open EoS table "//trim(filename))

        !> First read the header information, size from python file, and allocate
        read(lun, iostat=ios) tc%dim1, tc%dim2
        if (ios /= 0) call mpistop( &
            "read_eos_from_file: cannot read dims from "//trim(filename))

        !> Then read the header information, bounds from python file
        read(lun, iostat=ios) tc%var1_min, tc%var1_max, tc%var2_min, tc%var2_max
        if (ios /= 0) call mpistop( &
            "read_eos_from_file: cannot read axis bounds from "//trim(filename))

        if (allocated(tc%table)) deallocate(tc%table)
        allocate(tc%table(tc%dim1, tc%dim2))

        !> Then read the 2D cube from the python file
        read(lun, iostat=ios) tc%table
        if (ios /= 0) call mpistop( &
            "read_eos_from_file: cannot read table payload from "//trim(filename))

        if (logtable) then
            tc%table = dlog10(tc%table) !> Store tables in log10 space for interpolation
        endif

        !> Optional adaptive-grid trailer. Uniform (legacy) tables hit EOF here
        !> and the iostat-guarded read keeps is_uniform = .true.
        tc%is_uniform = .true.
        if (allocated(tc%var1_nodes)) deallocate(tc%var1_nodes)
        if (allocated(tc%var2_nodes)) deallocate(tc%var2_nodes)

        read(lun, iostat=ios) trailer_flag
        if (ios == 0) then
            if (trailer_flag == 1) then
                tc%is_uniform = .false.
                allocate(tc%var1_nodes(tc%dim1))
                allocate(tc%var2_nodes(tc%dim2))
                read(lun) tc%var1_nodes
                read(lun) tc%var2_nodes
                if (mype == 0) then
                    print*, trim(filename), " (adaptive grid trailer detected)"
                endif
            else
                if (mype == 0) then
                    print*, "WARNING: unrecognised trailer flag ", trailer_flag, &
                        " in ", trim(filename), " - treating as uniform"
                endif
            endif
        end if

        close(lun)

    end subroutine read_eos_from_file

end module mod_eos_LTE_tables
!> Needs a line after to pass the preprocessor
