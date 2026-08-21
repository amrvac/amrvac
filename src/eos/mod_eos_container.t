!=============================================================================
!> EoS state container -- the single thermodynamic authority for AMRVAC.
!>
!> Defines the eos_container derived type: its scalar fields (gamma, He_abundance,
!> eos_type, method_id, ionE, the unit-normalisation factors) and its
!> procedure-pointer slots -- to_conserved/to_primitive, p_to_e,
!> get_thermal_pressure, get_csound2, get_Rfactor, get_temperature_from_{eint,
!> etot,pressure}, get_Te, get_ne_nH, get_rho/get_nH, update_eos -- which ARE the
!> EoS public API the physics modules call. Also defines eos_table_container (one
!> lookup table + its axes/guards) and the global `eos` instance plus the
!> cached-log10(nH) wextra index iw_log_nH.
!>
!> Pure declarations -- no logic. The routines bound to these pointers live in
!> the EoS sub-modules (mod_eos*) and the physics seams (mod_(m)hd_eos).
!=============================================================================
module mod_eos_container
    use mod_global_parameters
    implicit none
    private

    public :: eos_container, eos_table_container

    !> LTE method selector, decoded once from eos%method into eos%method_id so
    !> the runtime kernels dispatch on an integer rather than a string compare.
    integer, parameter, public :: EOS_STATE = 1, EOS_ANALYTIC = 2, EOS_ENTROPY = 3

    !> EoS-type selector, decoded once from eos%eos_type into eos%type_id. Used by
    !> the mode-specific kernels to mpistop if called under the wrong eos_type
    !> (e.g. an LTE table kernel invoked while running FI or PI).
    integer, parameter, public :: EOS_TYPE_FI = 1, EOS_TYPE_LTE = 2, EOS_TYPE_PI = 3

    abstract interface
        subroutine convert_condition(ixI^L, ixO^L, w, x)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
        end subroutine convert_condition

        subroutine update_eos_spaces(ixI^L, ixO^L, w, x)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(inout) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
        end subroutine update_eos_spaces

        subroutine get_eos_variable(ixI^L, ixO^L, w, x, res)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(in) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
            double precision, intent(out)   :: res(ixI^S)
        end subroutine get_eos_variable

        subroutine get_eos_variable_alt(w, x, ixI^L, ixO^L,  res)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(in) :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
            double precision, intent(out)   :: res(ixI^S)
        end subroutine get_eos_variable_alt

        subroutine get_ne_nH_iface(ixI^L, ixO^L, w, x, ne, nH)
            use mod_global_parameters
            integer, intent(in)             :: ixI^L, ixO^L
            double precision, intent(in)    :: w(ixI^S, nw)
            double precision, intent(in)    :: x(ixI^S, 1:ndim)
            double precision, intent(out)   :: ne(ixI^S), nH(ixI^S)
        end subroutine get_ne_nH_iface
    end interface

    type eos_table_container

        character(len=std_len) :: filename
        double precision, allocatable :: table(:,:)
        integer :: dim1, dim2
        double precision :: var1_min, var1_max, var2_min, var2_max
        double precision :: step_inv_1, step_inv_2  !> precomputed (dim-1)/(max-min) for each axis

        !> Curvature-adaptive (non-uniform) grid support.
        !> When .true. (default), the grid is uniform in log10 space and
        !> the existing affine index calculation applies. When .false.,
        !> the table file carries an explicit pair of axis node arrays
        !> read from the binary trailer; index lookup uses binary search
        !> and the local fractional coordinate becomes
        !>   t1 = (var1 - var1_nodes(i)) / (var1_nodes(i+1) - var1_nodes(i))
        !> The PCHIP/bicubic kernels themselves are unchanged.
        !>
        !> Naming follows the existing var1/var2 convention:
        !>   var1_nodes has length dim1, var2_nodes has length dim2.
        logical :: is_uniform = .true.
        double precision, allocatable :: var1_nodes(:), var2_nodes(:)

        !> Guard (bucket) array for O(1) adaptive index lookup. Built when
        !> is_uniform = .false. via eos_build_guards. For each axis:
        !>   guard_M    : number of buckets covering [v(1), v(N)]
        !>   guard_scale: M / (v(N) - v(1))
        !>   guard(:)   : 1-based array, size M; guard(k) is the smallest
        !>                1-based node index ii such that v(ii) >= left
        !>                edge of bucket k (i.e., v(1) + (k-1)*span/M).
        !> Lookup:
        !>   k  = 1 + min(M-1, max(0, floor((val - v(1)) * scale)))
        !>   ii = guard(k); while (ii < N .and. v(ii) < val) ii = ii + 1
        !> M is sized at build time so the safety while-loop iterates at
        !> most once for typical curvature-equidistributed grids.
        integer :: guard_M_1 = 0, guard_M_2 = 0
        double precision :: guard_scale_1 = 0.0d0, guard_scale_2 = 0.0d0
        integer, allocatable :: guard_1(:), guard_2(:)

    end type eos_table_container

    type eos_container

        character(len=std_len) :: eos_type
        character(len=20)     :: method = 'state'         !> 'state', 'analytic' or 'entropy' ('tables' = legacy alias for 'state')
        integer               :: method_id = EOS_STATE   !> integer form of method (set at init)
        integer               :: type_id   = EOS_TYPE_FI  !> integer form of eos_type (set at init)
        character(len=20)     :: gamma1_method = 'exact'   !> 'exact' or 'effective'
        character(len=20)     :: p2eint_method = 'table'   !> 'table' (fast) or 'bisect' (accurate)
        character(len=20)     :: pi_table = 'chromosphere'  !> PI backend table: 'chromosphere'|'flare'|'prominence'
        logical :: ionE
        character(len=std_len) :: table_location
        double precision :: He_abundance
        double precision :: gamma
        double precision :: gamma_minus_1
        double precision :: inv_gamma
        double precision :: inv_gamma_minus_1
        double precision :: inv_squared_c0
        double precision :: inv_squared_c
        double precision :: nH2rhoFactor

        !> Fully-ionised regime bypass constants (precomputed in eos_finalise)
        double precision :: eion_per_nH          !> Total ionisation energy per nH [code units]
        double precision :: eint_rho_FI_threshold !> eint/rho above which gas is fully ionised [code units]
        double precision :: p_rho_FI_threshold    !> p/rho above which gas is fully ionised [code units]
        double precision :: n_per_nH_FI          !> Total particles per nH when FI = 2 + 3*A_He
        double precision :: neOnH_FI             !> ne/nH when fully ionised = 1 + 2*A_He

        !>Leaving a comment here to remind me that it's a bad idea to have anything outside
        !> of the w() array if one wants to /ensure/ the fields are consistent
        !> given OMP directives and load balancing

        ! Expected components of the EoS
        type(eos_table_container) :: neOnH
        type(eos_table_container) :: T
        type(eos_table_container) :: p2eint
        type(eos_table_container) :: gamma1    !> First adiabatic index Gamma_1(nH, eint/nH)
        type(eos_table_container) :: gamma1_p  !> First adiabatic index Gamma_1(nH, p/nH) for fast csound2
        type(eos_table_container) :: eint_from_T  !> Inverse T table: eint/nH(nH, T)
        type(eos_table_container) :: log_p        !> Merged log10(p/nH)(nH, eint/nH) for WB bisection
        type(eos_table_container) :: p_over_nH   !> (1+He+y)*T for direct p lookup in to_primitive

        !> Entropy-method tables (eos%method == 'entropy'). Each quantity Q is
        !> stored as four containers -- Q, Q_x, Q_y, Q_xy (value plus the two
        !> first derivatives and the cross derivative) -- fully determining the
        !> bicubic Hermite polynomial in each cell (mod_eos_LTE_entropy). Five
        !> quantities:
        !>   Forward (log_nH, log_eint/nH): Tfwd, pfwd, neOnH.
        !>   Inverse (log_nH, log_p/nH):    eintP, g1p (= Gamma_1).
        !>   Inverse (log_nH, log_T):       eintT.
        !> The inverse tables are bisected against the analytic Saha solver at
        !> build time and frozen to disk, so every runtime query is one lookup.
        !
        !> Forward (lr, le): ionisation fraction. The neOnH value container
        !> is the existing neOnH field above; add three derivative tables.
        type(eos_table_container) :: neOnH_x  !> d(ne/nH)/dx
        type(eos_table_container) :: neOnH_y  !> d(ne/nH)/dy
        type(eos_table_container) :: neOnH_xy !> d2(ne/nH)/(dx dy)
        !
        !> Forward (lr, le): direct T and p tables. We do NOT derive T from
        !> ds/dy of the s polynomial at runtime; the bicubic Hermite gives
        !> O(h^3) error in derivatives vs O(h^4) in values. Storing T and p
        !> as their own bicubic Hermite tables (built from the same Saha
        !> call as s) keeps O(h^4) accuracy and preserves Maxwell consistency
        !> at every node (T, p, s, gamma1 all from one Saha state).
        type(eos_table_container) :: Tfwd     !> T [K] (CGS-stored) at (lr, le)
        type(eos_table_container) :: Tfwd_x   !> dT/dx
        type(eos_table_container) :: Tfwd_y   !> dT/dy
        type(eos_table_container) :: Tfwd_xy  !> d2T/(dx dy)
        type(eos_table_container) :: pfwd     !> p [erg/cm^3] (CGS-stored) at (lr, le)
        type(eos_table_container) :: pfwd_x   !> dp/dx
        type(eos_table_container) :: pfwd_y   !> dp/dy
        type(eos_table_container) :: pfwd_xy  !> d2p/(dx dy)
        !
        !> Inverse (lr, lp): eint from (rho, p), replaces p2eint bisection
        !> at runtime. Stored as log10(eint/nH), CGS axis units before shift.
        type(eos_table_container) :: eintP    !> log10(eint/nH) [CGS log10]
        type(eos_table_container) :: eintP_x  !> d(eintP)/dx    x = log10(nH)
        type(eos_table_container) :: eintP_y  !> d(eintP)/dy    y = log10(p/nH)
        type(eos_table_container) :: eintP_xy !> d2(eintP)/(dx dy)
        !
        !> Inverse (lr, lp): gamma_1 at (rho, p) for csound2-from-pressure paths
        type(eos_table_container) :: g1p      !> Gamma_1(rho, p) dimensionless
        type(eos_table_container) :: g1p_x    !> dGamma_1/dx
        type(eos_table_container) :: g1p_y    !> dGamma_1/dy    y = log10(p/nH)
        type(eos_table_container) :: g1p_xy   !> d2Gamma_1/(dx dy)
        !
        !> Inverse (lr, lT): eint from (rho, T). Used by IC, HSE BC, cooling.
        type(eos_table_container) :: eintT    !> log10(eint/nH) [CGS log10]
        type(eos_table_container) :: eintT_x  !> d(eintT)/dx   x = log10(nH)
        type(eos_table_container) :: eintT_y  !> d(eintT)/dy   y = log10(T)
        type(eos_table_container) :: eintT_xy !> d2(eintT)/(dx dy)
        !> Interleaved Group A table: T, neOnH, p_over_nH at same (nH, eint/nH) grid
        double precision, allocatable :: table_eint_il(:,:,:)  !> (3, dim1, dim2)
        integer :: il_nq = 3  !> number of interleaved quantities

        procedure (convert_condition), pointer, nopass :: to_conserved => null()
        procedure (convert_condition), pointer, nopass :: to_primitive => null()
        procedure (convert_condition), pointer, nopass :: p_to_e => null()
        procedure (update_eos_spaces), pointer, nopass :: update_eos => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_thermal_pressure => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_temperature_from_eint => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_temperature_from_etot => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_temperature_from_pressure => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_Rfactor => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_rho => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_nH => null()
        procedure (get_ne_nH_iface), pointer, nopass :: get_ne_nH => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_Te => null()
        procedure (get_eos_variable_alt), pointer, nopass :: get_csound2 => null()

    contains

        !> Internal setter for gamma: updates gamma + the three derived cached
        !> values (gamma_minus_1, inv_gamma, inv_gamma_minus_1) atomically.
        !> PRIVATE on purpose: gamma is set ONLY via the parfile (&eos_list gamma=),
        !> read once at eos_init, so the derived constants can never go stale.
        !> User code must not set gamma; assigning eos%gamma directly is caught by
        !> the consistency check in (m)hd_check_params.
        procedure, private :: set_gamma => eos_set_gamma

    end type eos_container

    !> The single EoS state object, allocated in eos_init and shared (read-mostly)
    !> across all EoS sub-modules and the physics layer.
    type(eos_container), allocatable, public :: eos

    !> wextra index for the cached log10(nH) used during STS substeps
    !> (-1 = not allocated). Set by the physics module that owns the wextra slot.
    integer, public :: iw_log_nH = -1

contains

    subroutine eos_set_gamma(this, gamma_new)
        class(eos_container), intent(inout) :: this
        double precision, intent(in)        :: gamma_new

        this%gamma = gamma_new
        this%gamma_minus_1 = this%gamma - 1.0d0
        this%inv_gamma = 1.0d0 / this%gamma
        if (abs(this%gamma_minus_1) > 1.0d-12) then
            this%inv_gamma_minus_1 = 1.0d0 / this%gamma_minus_1
        else
            this%inv_gamma_minus_1 = 0.0d0
        end if
    end subroutine eos_set_gamma

end module mod_eos_container
!> Needs a line after to pass the preprocesor