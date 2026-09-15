!> Equation of state for AMRVAC, handled through a single eos_container object.
!>
!> Two gas models (eos_type):
!>   FI  -- fully ionised ideal gas (constant gamma, mu).
!>   LTE -- local thermodynamic equilibrium with partial ionisation of H (+He),
!>          via one of three interchangeable methods (eos_method):
!>            state    -- per-quantity (T,p,ne/nH,eint/p) Saha-equilibrium tables,
!>                        PCHIP/bilinear interpolated ('tables' accepted as legacy alias).
!>            entropy  -- bicubic-Hermite reconstruction (see mod_eos_lte_entropy).
!>            analytic -- on-the-fly H-only Saha solve (see saha_* routines).
!>          Tables are generated offline and read as binary; their axes are
!>          (log10 nH, log10 eint/nH) in code units after the unit shift applied
!>          in eos_finalise. ionE optionally folds ionisation energy into eint.
!>
!> The code-unit normalisation is kept so an FI atmosphere recovers R = 1;
!> update_eos_LTE relies on this. update_eos refreshes the cached EoS fields
!> (Te, ne) once per RK substep, and for boundary cells via update_eos_4_bc;
!> call eos%update_eos manually before any extra output.
!>
!> Jack Jenkins (27.01.2026); legacy EoS type 1: uniform tabulated (by Chris Osborne) H(+He) LTE EoS, in
!> consultation with Damien Przybylski (MURaM). Routed through legacy usr_Rfactor approach.
!> Jack Jenkins (05.05.2026); EoS type 2 finished, now separate module that each hd/mhd physics module 
!>                            from upstream/amrvac3.3 route through.
!> Jack Jenkins (11.05.2026); Benchmarked curvature-density-hybrid-tabulated EoS tables.
!> Jack Jenkins (25.05.2026); Single curvature-density-tabulated entropy potential EoS prototyped. Validated
!>                            gamma(T) against Ballester et al. 2001 to ~ 1.1% .
!> Jack Jenkins (16.06.2026); Folded in fld, ffhd, and new partial_ionisation modules to be EoS compliant 
!>                            and refactored for readability.

module mod_eos
    use mod_global_parameters
    use mod_eos_container
    use mod_eos_shared_functions
    use mod_eos_interp
    use mod_eos_lte_saha
    use mod_eos_lte_tables
    use mod_eos_fi
    use mod_eos_lte
    use mod_eos_pi, only: eos_init_PI, eos_finalise_PI
    use mod_timing
    use mod_comm_lib, only: mpistop

    implicit none
    private

    character(len=std_len) :: AMRVAC_DIR

    !> The EoS state object (eos) and the cached-log10(nH) wextra index
    !> (iw_log_nH) now live in mod_eos_container so every EoS sub-module can
    !> reach them; re-export here so existing `use mod_eos` callers still see them.
    public :: eos, iw_log_nH

    !> Lifecycle (called from amrvac.t)
    public :: eos_init, eos_finalise, prepare_eos_w_fields
    !> Type-agnostic helpers, safe in any eos_type. The mode-specific scalar
    !> kernels (LTE T/y/p/eint/gamma1, saha_*, the PI state/eint/csound2 backend
    !> and fl-port shims) are deliberately NOT re-exported here: the hd/mhd/ffhd
    !> seams reach them via their sub-modules (mod_eos_lte / mod_eos_lte_saha /
    !> mod_eos_pi), and each kernel mpistops if called under the wrong
    !> eos_type/method. So `use mod_eos` no longer exposes mode-specific routines.
    public :: get_ne_nH, eos_get_log_T_floor
    !> Ideal-gas Gamma_1 (returns eos%gamma); bound to phys_get_gamma1 for FI and
    !> no-energy PI, valid in any ideal-gas context, so kept on the facade.
    public :: get_gamma1_FI

contains

    !> Lifecycle: read and validate parameters, allocate, load and finalise
    !> the chosen backend. Called from amrvac.t.
    !> Read this module"s parameters from a file
    subroutine eos_read_params(files)
        character(len=*), intent(in) :: files(:)
        integer                      :: n

        !> Namelist mirror variables: the parfile uses these simple names
        !> (e.g. eos_type='FI'); they are committed to eos% in one block after
        !> the read. Derived eos% quantities (gamma_minus_1, method_id, ...) are
        !> NOT mirrored here -- they are computed from these inputs below.
        character(len=std_len) :: eos_type, table_location
        character(len=20)      :: eos_method, gamma1_method, p2eint_method, pi_table
        double precision       :: He_abundance, gamma
        logical                :: ionE

        namelist /eos_list/ eos_type, table_location, ionE, He_abundance, gamma, &
            eos_method, gamma1_method, p2eint_method, pi_table

        !> defaults
        eos_type          = 'FI'
        eos_method        = 'state'
        pi_table          = 'chromosphere'
        gamma1_method     = 'exact'
        p2eint_method     = 'table'
        ionE              = .false.
        He_abundance      = 0.1d0
        gamma             = 5.0d0/3.0d0
        call get_environment_variable("AMRVAC_DIR", AMRVAC_DIR)
        table_location    = trim(AMRVAC_DIR)//"/src/tables/eos_tables/uniform256/"

        !> read parfile(s): each successive file overrides the previous
        do n = 1, size(files)
            open(unitpar, file=trim(files(n)), status="old")
            read(unitpar, eos_list, end=111)
            111     close(unitpar)
        end do

        !> commit raw inputs to eos% (plain copies, no logic)
        eos%eos_type          = eos_type
        eos%method            = eos_method
        eos%pi_table          = pi_table
        eos%gamma1_method     = gamma1_method
        eos%p2eint_method     = p2eint_method
        eos%table_location    = table_location
        eos%He_abundance      = He_abundance
        eos%gamma             = gamma
        eos%ionE              = ionE

        !> derive secondary quantities from the raw inputs
        eos%gamma_minus_1 = eos%gamma - 1.0d0
        eos%inv_gamma     = 1.0d0 / eos%gamma
        if (abs(eos%gamma_minus_1) > 1.0d-12) then
            eos%inv_gamma_minus_1 = 1.0d0 / eos%gamma_minus_1
        else
            eos%inv_gamma_minus_1 = 0.0d0
        end if
        eos%nH2rhoFactor = 1.0d0       !> FI default; (m)hd_physical_units may reset
        select case (eos%method)
        case ('analytic'); eos%method_id = EOS_ANALYTIC
        case ('entropy');  eos%method_id = EOS_ENTROPY
        case default;      eos%method_id = EOS_STATE
        end select
        select case (eos%eos_type)
        case ('LTE'); eos%type_id = EOS_TYPE_LTE
        case ('PI');  eos%type_id = EOS_TYPE_PI
        case default; eos%type_id = EOS_TYPE_FI
        end select

        !> report
        if (mype == 0) then
            write(*,*) "EoS type: "//trim(eos%eos_type)
            if (eos%eos_type == 'LTE') &
                write(*,*) "EoS table location: "//trim(eos%table_location)
        end if

        !> validate / coerce (all consistency guards live in one place)
        call eos_validate_params()

    end subroutine eos_read_params

    !> Validate and coerce the parsed eos_list against the chosen eos_type /
    !> method. Grouped by concern so the rules are readable in one block (rather
    !> than interleaved with the field copies in eos_read_params). May coerce
    !> eos% fields (e.g. FI forces ionE=.false.) and mpistop on illegal combos.
    subroutine eos_validate_params()

        !> eos_type x ionE
        if (eos%eos_type == 'FI' .and. eos%ionE) then
            eos%ionE = .false.
            if(mype==0) write(*,*) "WARNING: eos_type = 'FI' requires ionE = .false."
        end if

        !> PI ('partial approximations') table selector
        if (eos%eos_type == 'PI') then
            select case (trim(eos%pi_table))
            case ('chromosphere','flare','prominence')
                ! valid
            case default
                call mpistop("eos_type='PI': pi_table must be "// &
                     "'chromosphere', 'flare', or 'prominence'")
            end select
            !> ionisation-energy mode (ionE) is supported only for the T-only
            !> tables; the prominence table is (T,p)-dependent -> no-energy only.
            !> This also gates the only PI configuration that cannot model He
            !> ionisation energy: the chromosphere/flare energy paths now fold
            !> the He ionisation energy into eint (mod_eos_pi_tables,
            !> ionization_eps_ion_of_degrees), consistent with the He electron
            !> count already in the R-factor, so He>0 energy mode is allowed there.
            if (eos%ionE .and. trim(eos%pi_table) == 'prominence') then
                call mpistop("eos_type='PI' with ionE=T requires a T-only table "// &
                     "(chromosphere or flare); prominence is no-energy only")
            end if
            if(mype==0) write(*,*) "EoS PI table: "//trim(eos%pi_table)// &
                 ", ionE=", eos%ionE
        end if

        !> eos_method validity and eos_type x method ('tables' kept as a legacy
        !> alias for 'state' -- both decode to EOS_STATE via the case default).
        if (eos%method /= 'analytic' .and. eos%method /= 'state' &
            .and. eos%method /= 'tables' .and. eos%method_id /= EOS_ENTROPY) then
            call mpistop("Unknown eos_method: "//trim(eos%method)// &
                " (expected 'state', 'analytic', or 'entropy')")
        end if
        if (eos%eos_type == 'FI' .and. eos%method_id == EOS_ANALYTIC) then
            call mpistop("eos_method='analytic' requires eos_type='LTE'")
        end if

        !> gamma1_method validity
        if (eos%gamma1_method /= 'exact' .and. eos%gamma1_method /= 'effective' &
            .and. eos%gamma1_method /= 'constant') then
            call mpistop("Unknown gamma1_method: "//trim(eos%gamma1_method)// &
                " (expected 'exact', 'effective', or 'constant')")
        end if

        !> analytic-method caveats (H-only Saha)
        if (eos%method_id == EOS_ANALYTIC) then
            if (eos%He_abundance > 0.0d0 .and. mype == 0) write(*,*) &
                "WARNING: eos_method='analytic' is H-only. He_abundance=", &
                eos%He_abundance, " used in pressure law but not in Saha ionization."
            ! effective gamma1 needs T,neOnH tables -> not available analytically
            if (eos%gamma1_method == 'effective') then
                if (mype == 0) write(*,*) &
                    "WARNING: gamma1_method='effective' not applicable with analytic EoS."
                if (mype == 0) write(*,*) &
                    "         Falling back to gamma1_method='exact' (analytic 2D table)."
                eos%gamma1_method = 'exact'
            end if
            if (.not. eos%ionE .and. mype == 0) then
                write(*,*) "WARNING: eos_method='analytic' without ionE"
                write(*,*) "  reduces to ideal gas. Consider eos_type='FI'."
            end if
        end if

        !> LTE+ionE effective-gamma1 note
        if (eos%eos_type == 'LTE' .and. eos%ionE .and. &
            eos%gamma1_method == 'effective') then
            if (mype == 0) write(*,*) &
                "NOTE: gamma1_method='effective' builds gamma_eff = 1+p/eint table."
        end if

    end subroutine eos_validate_params

    !> Phase 'create' (before units are known): allocate the EoS object, read
    !> its &eos_list parameters, and load the raw tables for the chosen method.
    !> The getters wired here must be live before (m)hd_link_eos runs.
    subroutine eos_init()
        allocate(eos)
        call eos_read_params(par_files)

        !> Type-agnostic getters: always live (eos% pointer targets).
        eos%get_rho   => get_rho
        eos%get_nH    => get_nH
        eos%get_ne_nH => get_ne_nH

        !> Per-type create-phase init (each owns its raw table loads / pointer
        !> wiring that must be live before (m)hd_link_eos runs), in the matching
        !> mod_eos_<TYPE> module.
        select case (eos%eos_type)
        case ('FI');  call eos_init_FI()
        case ('LTE'); call eos_init_LTE()
        case ('PI');  call eos_init_PI()
        case default
            call mpistop('eos_init: eos_type '//trim(eos%eos_type)//' not recognised')
        end select
    end subroutine eos_init

    !> Phase 'commit' (after units are known): finalise the dispatch for the
    !> loaded physics -- wire the method's runtime pointers and finalise its
    !> tables. The cross-module wiring that follows (into the ghost-cell update
    !> and the physics source terms) is done by the caller in amrvac.t, so the
    !> EoS never reaches back into those modules.
    !> Finalise the EoS once units are set: dispatch on eos_type to the per-type
    !> finaliser (each owns its pointer wiring + backend table loading, living in
    !> mod_eos_fi / mod_eos_lte / mod_eos_pi), then wire the shared total-energy
    !> temperature getter. Called from amrvac.t after phys init.
    subroutine eos_finalise()
        if (.not. allocated(eos)) &
            call mpistop('eos_finalise: eos not allocated')

        select case (eos%eos_type)
        case ('LTE')
            call eos_finalise_LTE()
        case ('FI')
            call eos_finalise_FI()
        case ('PI')
            call eos_finalise_PI()
        case default
            call mpistop('eos_finalise: eos_type '//trim(eos%eos_type)//' not recognised')
        end select

        eos%get_temperature_from_etot => get_temperature_from_etot
    end subroutine eos_finalise

    subroutine prepare_eos_w_fields()
        integer :: iigrid, igrid

        ! fill eos space of all grids, inc boundaries
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
            call eos%update_eos(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
        end do

    end subroutine prepare_eos_w_fields

end module mod_eos
!> Needs a line after to pass the preprocesor