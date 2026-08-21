!=============================================================================
!> LTE (Saha-table) EoS kernels and finalise for the eos% family.
!>
!> Carved out of mod_eos.t (the per-type split). Owns everything LTE: the
!> per-substep update_eos_LTE (caches Te_/Ne_), the cached-Te/Rfactor getters,
!> the scalar table kernels (T, y, p2eint, gamma1, p/nH, eint, bisections), the
!> per-method table finalisers, and eos_finalise_LTE (the LTE arm of the
!> eos_finalise dispatcher). Depends only on the EoS leaf modules + shared
!> accessors, never on mod_eos -> no circular use. mod_eos re-exports the public
!> kernels so existing `use mod_eos` callers are unaffected.
!=============================================================================
module mod_eos_LTE
    use mod_global_parameters
    use mod_eos_container
    use mod_eos_shared_functions
    use mod_eos_interp
    use mod_eos_LTE_saha
    use mod_eos_LTE_entropy
    use mod_eos_LTE_tables
    use mod_eos_LTE_state
    use mod_timing
    use mod_comm_lib, only: mpistop

    implicit none
    private

    !> update_eos_LTE, get_Te_LTE, get_temperature_from_eint_LTE and
    !> Rfactor_from_LTE are PRIVATE: all bound to eos% pointers inside
    !> eos_finalise_LTE, so no external module names them.
    !> get_temperature_from_eint_fast_LTE stays public -- the seam binds it to
    !> tc_fl%get_temperature_from_eint (fast TC path).
    public :: get_temperature_from_eint_fast_LTE
    !> Scalar table kernels (re-exported by mod_eos for the hd/mhd/ffhd seams).
    !> T_and_y_from_nH_eint is internal-only (entropy update path).
    public :: y_from_nH_eint, T_from_nH_eint
    public :: p2eint_from_nH_p, gamma1_from_nH_p, p_nH_from_eint
    public :: eint_nH_from_T, eint_from_p_bisect, eos_get_eintT_grid
    !> LTE arms of the eos_init / eos_finalise dispatchers
    public :: eos_init_LTE, eos_finalise_LTE

contains

    !> LTE arm of eos_init (before units are known): dispatch on eos_method to
    !> the chosen method's loader, each owned by its method module (analytic loads
    !> nothing). The LTE temperature-from-pressure getter is a later pass, so
    !> unlike FI/PI it is left unset here.
    subroutine eos_init_LTE()
        select case (eos%method_id)
        case (EOS_ANALYTIC); call load_analytic_LTE()
        case (EOS_ENTROPY);  call load_entropy_LTE()
        case default;        call load_state_LTE()   ! EOS_STATE
        end select
    end subroutine eos_init_LTE

    !> LTE arm of eos_finalise: wire the LTE runtime pointer targets, then
    !> dispatch on eos_method to the per-method finaliser, each owned by its
    !> method module (analytic->mod_eos_LTE_saha, entropy->mod_eos_LTE_entropy,
    !> tables->mod_eos_LTE_tables).
    subroutine eos_finalise_LTE()
        eos%update_eos => update_eos_LTE
        ! eos%get_thermal_pressure is set by (m)hd_link_eos
        eos%get_temperature_from_eint => get_temperature_from_eint_LTE
        eos%get_Te => get_Te_LTE
        !> R-factor: physics-independent (cached-Ne lookup), bound here rather
        !> than in the hd/mhd/ffhd link arms so the target stays private.
        eos%get_Rfactor => Rfactor_from_LTE
        select case (eos%method_id)
        case (EOS_ANALYTIC); call finalise_analytic_LTE()
        case (EOS_ENTROPY);  call finalise_entropy_LTE()
        case default;        call finalise_state_LTE()   ! EOS_STATE
        end select
    end subroutine eos_finalise_LTE


    subroutine update_eos_LTE(ixI^L, ixO^L, w, x)
        use mod_physics
        !> This routine is called before each RK substep
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw)

        double precision :: wlocal(ixI^S,1:nw)
        double precision :: deltaEfactor(ixI^S)
        double precision :: prev_y(ixI^S), new_y(ixI^S)
        double precision :: pth(ixI^S)
        double precision :: nH_in(ixI^S), nH(ixI^S), eint_in(ixI^S)
        double precision :: Rfactor_FI
        double precision :: yy, eint_nH_floor
        double precision :: time0
        integer :: ix^D

        timeeos0 = MPI_WTIME() !> For monitoring cost of eos module

        wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
        call phys_e_to_ei(ixI^L,ixO^L,wlocal,x) !> wlocal now contains internal energy, NOT total energy

        pth(ixO^S) = eos%gamma_minus_1 * wlocal(ixO^S,iw_e) !> pressure from internal energy only - SHOULD ONLY BE USED WHEN IonE IS UNIMPORTANT

        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        ! Enforce internal energy floor for EoS lookup only.
        ! Prevents NaN from dlog10(eint<0) after strong rarefactions where
        ! kinetic energy can numerically exceed total energy.
        ! The conserved energy w(iw_e) is NOT modified: the physical fluxes
        ! based on the small positive pressure will naturally restore the cell.
        if (eos%method_id == EOS_ANALYTIC) then
            ! Analytic: floor to ~100 K equivalent (no tables loaded)
            ! In code units: eint/nH = T_code / (gamma-1) for neutral gas
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                eint_nH_floor = nH(ix^D) * eos%inv_gamma_minus_1 * 100.0d0 / unit_temperature
                wlocal(ix^D,iw_e) = max(wlocal(ix^D,iw_e), eint_nH_floor)
            {end do\}
        else
            nH_in(ixO^S) = dlog10(nH(ixO^S))
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                !> eos%T%var2_min holds the forward log(eint/nH) lower bound for
                !> all table methods (entropy sources it from Tfwd in eos_finalise).
                eint_nH_floor = nH(ix^D) * 10.0d0**eos%T%var2_min
                if (wlocal(ix^D,iw_e) < eint_nH_floor) then
                    wlocal(ix^D,iw_e) = eint_nH_floor
                end if
            {end do\}
            eint_in(ixO^S) = dlog10(wlocal(ixO^S,iw_e)) - nH_in(ixO^S)
        end if

        !> Constant FI Rfactor for the fully-ionised fast path.
        !> Must NOT use eos%get_Rfactor which reads iw_ne (uninitialised at IC).
        Rfactor_FI = eos%n_per_nH_FI / eos%nH2rhoFactor

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if (wlocal(ix^D,iw_e) / w(ix^D,iw_rho) > eos%eint_rho_FI_threshold) then
                !> Fully ionised: analytical formulae, no table lookups
                new_y(ix^D) = eos%neOnH_FI
                if (eos%ionE) then
                    !> T = (eint - eion*nH) * (gamma-1) / (Rfactor_FI * rho)
                    w(ix^D,iw_te) = eos%gamma_minus_1 &
                        * (wlocal(ix^D,iw_e) - eos%eion_per_nH * nH(ix^D)) &
                        / (Rfactor_FI * w(ix^D,iw_rho))
                else
                    w(ix^D,iw_te) = pth(ix^D) / (nH(ix^D) * (1.0d0 + eos%He_abundance + new_y(ix^D)))
                end if
            else
                !> Ionisation zone
                if (eos%method_id == EOS_ANALYTIC) then
                    !> Analytical Saha: solve quadratic for T and y
                    call saha_T_from_nH_eint(nH(ix^D), &
                        wlocal(ix^D,iw_e) / nH(ix^D), &
                        w(ix^D,iw_te), yy)
                    new_y(ix^D) = yy
                else
                    !> Table lookups: fused T+y for shared index computation
                    if (.not. eos%ionE) then
                        new_y(ix^D) = y_from_nH_eint(nH_in(ix^D),eint_in(ix^D))
                        w(ix^D,iw_te) = pth(ix^D) / (nH(ix^D) * (1.0d0 + eos%He_abundance + new_y(ix^D)))
                    else if (eos%method_id == EOS_ENTROPY) then
                        !> Entropy method dispatches T_and_y_from_nH_eint to the
                        !> biquintic-Hermite forward at any eint; no FI fallback
                        !> needed (and eos%T isn't loaded so var2_max is 0 here).
                        call T_and_y_from_nH_eint(nH_in(ix^D), eint_in(ix^D), &
                            w(ix^D,iw_te), yy)
                        new_y(ix^D) = yy
                    else
                        if (eint_in(ix^D) < eos%T%var2_max) then
                            call T_and_y_from_nH_eint(nH_in(ix^D), eint_in(ix^D), &
                                w(ix^D,iw_te), yy)
                            new_y(ix^D) = yy
                        else
                            !> Above-table fallback with ionisation energy correction
                            new_y(ix^D) = eos%neOnH_FI
                            w(ix^D,iw_te) = eos%gamma_minus_1 &
                                * (wlocal(ix^D,iw_e) - eos%eion_per_nH * nH(ix^D)) &
                                / (Rfactor_FI * w(ix^D,iw_rho))
                        end if
                    end if
                end if
            end if
        {end do\}

        w(ixO^S,iw_ne) = new_y(ixO^S) * nH(ixO^S)

        timeeos_update=timeeos_update+(MPI_WTIME()-timeeos0)

    end subroutine update_eos_LTE

    !> Thermodynamic getters shared verbatim by the hd/mhd seams (eos%/phys
    !> pointer targets). Only routines free of physics-module symbols live here:
    !> Rfactor_from_LTE uses the cached Ne_ field; get_gamma1_FI returns
    !> eos%gamma. The FI/LTE variants that read primitive indices (rho_,p_) or
    !> the physics RR stay in each seam.
    !> Rfactor from LTE EOS: R = (1 + He + ne/nH) / nH2rhoFactor
    !> Uses stored Ne_ wextra field from update_eos_LTE.
    subroutine Rfactor_from_LTE(w,x,ixI^L,ixO^L,Rfactor)
        use mod_global_parameters
        integer, intent(in) :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: Rfactor(ixI^S)

        double precision :: y_nH
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            y_nH = w(ix^D, iw_ne) * eos%nH2rhoFactor / w(ix^D, iw_rho)
            Rfactor(ix^D) = (1.0d0 + eos%He_abundance + y_nH) / eos%nH2rhoFactor
        {end do\}

    end subroutine Rfactor_from_LTE

    !>#######################################################################
    !> Temperature getters: cached Te accessors and temperature from the
    !> energy/pressure variable. The from-eint variants are split by eos_type
    !> (needed by the STS path: dt sees conserved, set_source sees eint -- see
    !> mod_thermal_conduction).
    !>#######################################################################

    subroutine get_Te_LTE(w,x,ixI^L,ixO^L,T)
        use mod_global_parameters
        integer, intent(in)           :: ixI^L, ixO^L
        double precision, intent(in)  :: w(ixI^S,1:nw)
        double precision, intent(in)  :: x(ixI^S,1:ndim)
        double precision, intent(out) :: T(ixI^S)

        T(ixO^S) = w(ixO^S,iw_te)

    end subroutine get_Te_LTE

    subroutine get_temperature_from_eint_LTE(w, x, ixI^L, ixO^L, res)
        !> Assumes input energy is internal energy.
        !> Includes FI bypass: cells with eint/rho above threshold skip table lookups.
        use mod_physics
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(out)   :: res(ixI^S)

        double precision :: nH(ixI^S),nH_in(ixI^S), eint_in(ixI^S)
        double precision :: Rfactor(ixI^S), Rfactor_FI, T_loc, y_loc
        integer :: ix^D

        timeeos0 = MPI_WTIME()

        Rfactor_FI = eos%n_per_nH_FI / (1.0d0 + 4.0d0*eos%He_abundance)

        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        nH_in(ixO^S) = dlog10(nH(ixO^S))
        eint_in(ixO^S) = dlog10(w(ixO^S,iw_e)) - nH_in(ixO^S)

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if (w(ix^D,iw_e) / w(ix^D,iw_rho) > eos%eint_rho_FI_threshold) then
                !> Fully ionised: T = (eint - eion*nH) * (gamma-1) / (Rfactor_FI * rho)
                res(ix^D) = eos%gamma_minus_1 &
                    * (w(ix^D,iw_e) - eos%eion_per_nH * nH(ix^D)) &
                    / (Rfactor_FI * w(ix^D,iw_rho))
            else if (eos%method_id == EOS_ANALYTIC) then
                call saha_T_from_nH_eint(nH(ix^D), &
                    w(ix^D,iw_e) / nH(ix^D), T_loc, y_loc)
                res(ix^D) = T_loc
            else
                res(ix^D) = T_from_nH_eint(nH_in(ix^D),eint_in(ix^D))
            endif
        {end do\}

        timeeos_Tfromei=timeeos_Tfromei+(MPI_WTIME()-timeeos0)
    end subroutine get_temperature_from_eint_LTE

    subroutine get_temperature_from_eint_fast_LTE(w, x, ixI^L, ixO^L, res)
        !> Fast TC variant: two-pass regime-aware bypass.
        !>
        !> Pass 1 (vectorised): compute FI formula for ALL cells as array ops.
        !>   The compiler vectorises this with SVML (no branches, pure arithmetic).
        !>
        !> Pass 2 (scalar): overwrite ionisation-zone cells with
        !>   bilinear table lookup + dexp.  Threshold check uses multiply
        !>   (eint <= threshold * rho) to avoid a 14-cycle scalar division.
        use mod_physics
        use mod_eos_LTE_entropy, only: entropy_T_from_nH_eint
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(in)    :: w(ixI^S,1:nw)
        double precision, intent(out)   :: res(ixI^S)

        double precision :: inv_Rfactor_FI, eion_rho_inv
        double precision :: log_nH_val, log_eint_nH_val
        double precision :: fx, fy, rx, ry, nH_loc, T_loc, y_loc
        integer :: jx, jy, jx1, jy1, ix^D
        double precision, parameter :: ln10 = 2.302585092994046d0

        if (eos%type_id /= EOS_TYPE_LTE) call mpistop("get_temperature_from_eint_fast_LTE called outside its eos_type (LTE)")
        timeeos0 = MPI_WTIME()

        !> Precompute scalar constants (avoid per-cell divisions)
        inv_Rfactor_FI = (1.0d0 + 4.0d0*eos%He_abundance) / eos%n_per_nH_FI
        eion_rho_inv = eos%eion_per_nH / eos%nH2rhoFactor

        !> Pass 1: FI formula for ALL cells (vectorisable array operations).
        !> T = (gamma-1) * (eint - eion_per_rho * rho) / (Rfactor_FI * rho)
        res(ixO^S) = eos%gamma_minus_1 * inv_Rfactor_FI &
            * (w(ixO^S,iw_e) - eion_rho_inv * w(ixO^S,iw_rho)) &
            / w(ixO^S,iw_rho)

        !> Pass 2: overwrite ionisation-zone cells.
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if (w(ix^D,iw_e) <= eos%eint_rho_FI_threshold * w(ix^D,iw_rho)) then
                if (eos%method_id == EOS_ANALYTIC) then
                    nH_loc = w(ix^D, iw_rho) / eos%nH2rhoFactor
                    call saha_T_from_nH_eint(nH_loc, &
                        w(ix^D,iw_e) / nH_loc, T_loc, y_loc)
                    res(ix^D) = T_loc
                else if (eos%method_id == EOS_ENTROPY) then
                    if (iw_log_nH > 0) then
                        log_nH_val = block%wextra(ix^D, iw_log_nH)
                    else
                        log_nH_val = dlog10(w(ix^D, iw_rho) / eos%nH2rhoFactor)
                    end if
                    log_eint_nH_val = dlog10(w(ix^D, iw_e)) - log_nH_val
                    res(ix^D) = entropy_T_from_nH_eint(eos%Tfwd, eos%Tfwd_x, &
                        eos%Tfwd_y, eos%Tfwd_xy, log_nH_val, log_eint_nH_val)
                else
                    if (iw_log_nH > 0) then
                        log_nH_val = block%wextra(ix^D, iw_log_nH)
                    else
                        log_nH_val = dlog10(w(ix^D, iw_rho) / eos%nH2rhoFactor)
                    end if
                    log_eint_nH_val = dlog10(w(ix^D, iw_e)) - log_nH_val

                    if (eos%T%is_uniform) then
                        !> Inlined uniform bilinear (hot path: no function call).
                        ry = max(0.0d0, min((log_nH_val - eos%T%var1_min) * eos%T%step_inv_1, &
                                             dble(eos%T%dim1-1)))
                        rx = max(0.0d0, min((log_eint_nH_val - eos%T%var2_min) * eos%T%step_inv_2, &
                                             dble(eos%T%dim2-1)))
                        jy = int(ry); jx = int(rx)
                        jy1 = min(jy+1, eos%T%dim1-1)
                        jx1 = min(jx+1, eos%T%dim2-1)
                        fy = ry - dble(jy); fx = rx - dble(jx)
                        res(ix^D) = dexp(ln10 * ( &
                            (1.0d0-fy)*((1.0d0-fx)*eos%T%table(jy+1,jx+1)  &
                                              + fx *eos%T%table(jy+1,jx1+1)) &
                           +      fy *((1.0d0-fx)*eos%T%table(jy1+1,jx+1)   &
                                              + fx *eos%T%table(jy1+1,jx1+1))))
                    else
                        !> Non-uniform: dispatch to the binary-search bilinear.
                        !> Slightly slower per cell (binary search adds ~8 cmps per axis)
                        !> but unavoidable when the grid is non-uniform.
                        res(ix^D) = dexp(ln10 * bilinear_lookup( &
                            log_nH_val, log_eint_nH_val, eos%T))
                    end if
                end if
            end if
        {end do\}

        timeeos_Tfromei = timeeos_Tfromei + (MPI_WTIME()-timeeos0)
    end subroutine get_temperature_from_eint_fast_LTE

    !>#######################################################################
    !> Scalar EoS kernels: pointwise lookups on (log10 nH, log10 eint/nH) or
    !> (log10 nH, log10 p/nH) in code units, called by the hd/mhd EoS layer.
    !> Each dispatches analytic / entropy / table. Ionisation fraction and
    !> temperature first, then pressure / Gamma_1 / eint.
    !>#######################################################################

    !> Ionization fraction from (log10 nH, log10 eint/nH) in code units.
    !> Dispatches: analytic -> Saha quadratic, tables -> PCHIP interpolation.
    double precision function y_from_nH_eint(nH, eint_nh) result(result_val)
        use mod_eos_LTE_entropy, only: entropy_y_from_nH_eint
        double precision, intent(in) :: nH, eint_nh
        double precision, parameter :: ln10 = 2.302585092994046d0
        double precision :: T_loc, y_loc, eint_rho

        if (eos%type_id /= EOS_TYPE_LTE) call mpistop("y_from_nH_eint called outside its eos_type (LTE)")
        if (eos%method_id == EOS_ANALYTIC) then
            ! FI bypass: skip Saha solve for fully ionized cells
            eint_rho = 10.0d0**eint_nh / eos%nH2rhoFactor
            if (eint_rho > eos%eint_rho_FI_threshold) then
                result_val = eos%neOnH_FI
                return
            end if
            call saha_T_from_nH_eint(10.0d0**nH, 10.0d0**eint_nh, T_loc, y_loc)
            result_val = y_loc
        else if (eos%method_id == EOS_ENTROPY) then
            ! neOnH stored linearly (not log10); single bicubic Hermite lookup.
            result_val = entropy_y_from_nH_eint(eos%neOnH, eos%neOnH_x, &
                eos%neOnH_y, eos%neOnH_xy, nH, eint_nh)
        else
            result_val = dexp(ln10 * bicubic_lookup(nH, eint_nh, eos%neOnH))
        end if
    end function y_from_nH_eint

    !> Temperature from (log10 nH, log10 eint/nH) in code units.
    !> Dispatches: analytic -> Saha bisection/Newton, tables -> PCHIP interpolation.
    double precision function T_from_nH_eint(nH, eint_nh) result(result_val)
        use mod_eos_LTE_entropy, only: entropy_T_from_nH_eint
        double precision, intent(in) :: nH, eint_nh
        double precision, parameter :: ln10 = 2.302585092994046d0
        double precision :: T_loc, y_loc, eint_rho, Rfactor_FI

        if (eos%type_id /= EOS_TYPE_LTE) call mpistop("T_from_nH_eint called outside its eos_type (LTE)")
        if (eos%method_id == EOS_ANALYTIC) then
            ! FI bypass: skip Saha solve for fully ionized cells
            eint_rho = 10.0d0**eint_nh / eos%nH2rhoFactor
            if (eint_rho > eos%eint_rho_FI_threshold) then
                Rfactor_FI = eos%n_per_nH_FI / (1.0d0 + 4.0d0*eos%He_abundance)
                result_val = eos%gamma_minus_1 &
                    * (10.0d0**eint_nh - eos%eion_per_nH) / Rfactor_FI
                return
            end if
            call saha_T_from_nH_eint(10.0d0**nH, 10.0d0**eint_nh, T_loc, y_loc)
            result_val = T_loc
        else if (eos%method_id == EOS_ENTROPY) then
            result_val = entropy_T_from_nH_eint(eos%Tfwd, eos%Tfwd_x, &
                eos%Tfwd_y, eos%Tfwd_xy, nH, eint_nh)
        else
            result_val = dexp(ln10 * bicubic_lookup(nH, eint_nh, eos%T))
        end if
    end function T_from_nH_eint

    !> Fused T+y lookup from (log10 nH, log10 eint/nH) in code units.
    !> Computes grid indices once, evaluates both T and y tables.
    !> Saves one index computation + better cache utilisation vs separate calls.
    subroutine T_and_y_from_nH_eint(log_nH, log_eint_nH, T_out, y_out)
        use mod_eos_LTE_entropy, only: entropy_T_and_y_from_nH_eint
        double precision, intent(in) :: log_nH, log_eint_nH
        double precision, intent(out) :: T_out, y_out
        double precision, parameter :: ln10 = 2.302585092994046d0
        double precision :: results(3)

        if (eos%method_id == EOS_ANALYTIC) then
            call saha_T_from_nH_eint(10.0d0**log_nH, 10.0d0**log_eint_nH, T_out, y_out)
        else if (eos%method_id == EOS_ENTROPY) then
            call entropy_T_and_y_from_nH_eint(eos%Tfwd, eos%Tfwd_x, &
                eos%Tfwd_y, eos%Tfwd_xy, &
                eos%neOnH, eos%neOnH_x, eos%neOnH_y, eos%neOnH_xy, &
                log_nH, log_eint_nH, T_out, y_out)
        else if (allocated(eos%table_eint_il)) then
            ! Interleaved PCHIP fast-path.
            ! (T, n_e/n_H, p/n_H), so cache behaviour is preserved.
            if (eos%T%is_uniform) then
                call interp_pchip_interleaved(log_nH, log_eint_nH, &
                    eos%table_eint_il, 3, eos%T%dim2, eos%T%dim1, &
                    eos%T%var1_min, eos%T%var1_max, &
                    eos%T%var2_min, eos%T%var2_max, results)
            else
                call interp_pchip_interleaved_nu(log_nH, log_eint_nH, &
                    eos%table_eint_il, 3, eos%T%dim2, eos%T%dim1, &
                    eos%T%var1_nodes, eos%T%var2_nodes, &
                    eos%T%guard_1, eos%T%guard_M_1, eos%T%guard_scale_1, &
                    eos%T%guard_2, eos%T%guard_M_2, eos%T%guard_scale_2, &
                    results)
            end if
            T_out = dexp(ln10 * results(1))
            y_out = dexp(ln10 * results(2))
        else
            ! No interleaved table built -- separate dispatcher calls
            T_out = dexp(ln10 * bicubic_lookup(log_nH, log_eint_nH, eos%T))
            y_out = dexp(ln10 * bicubic_lookup(log_nH, log_eint_nH, eos%neOnH))
        end if
    end subroutine T_and_y_from_nH_eint

    !> Pressure-to-eint ratio from (log10 nH, log10 p/nH) in code units.
    !> Dispatches: analytic -> Saha solve for eint/p, tables -> PCHIP interpolation.
    double precision function p2eint_from_nH_p(nH, ponH) result(result_val)
        use mod_eos_LTE_entropy, only: entropy_eint_from_nH_p
        double precision, intent(in) :: nH, ponH
        double precision :: nH_code, p_code, T_loc, y_loc, eint_nH_loc, p_rho

        if (eos%type_id /= EOS_TYPE_LTE) call mpistop("p2eint_from_nH_p called outside its eos_type (LTE)")
        if (eos%method_id == EOS_ANALYTIC) then
            nH_code = 10.0d0**nH
            p_code = nH_code * 10.0d0**ponH
            ! FI bypass: skip Saha solve for fully ionized cells
            p_rho = p_code / (nH_code * eos%nH2rhoFactor)
            if (p_rho > eos%p_rho_FI_threshold) then
                result_val = eos%inv_gamma_minus_1 &
                    + eos%eion_per_nH * nH_code / p_code
                return
            end if
            call saha_state_from_nH_p(nH_code, p_code, T_loc, y_loc, eint_nH_loc)
            result_val = eint_nH_loc * nH_code / p_code
        else if (eos%method_id == EOS_ENTROPY) then
            ! One bicubic-Hermite lookup on the eintP table.
            result_val = entropy_eint_from_nH_p(eos%eintP, eos%eintP_x, &
                eos%eintP_y, eos%eintP_xy, nH, ponH)
        else
            result_val = bicubic_lookup(nH, ponH, eos%p2eint)
        end if
    end function p2eint_from_nH_p

    !> Gamma_1 from pressure-indexed table: (log10 nH, log10 p/nH) -> Gamma_1.
    !> For 'entropy' the conversion p -> eint -> Gamma_1 via formula is intended
    !> to keep Maxwell consistency with the forward at runtime. The p->eint inverse 
    !> is one bisection per cell -- non-trivial cost; this function is in the hot path
    !> via hd_get_csound2_LTE.
    double precision function gamma1_from_nH_p(log_nH, log_p_nH) result(g1)
        use mod_eos_LTE_entropy, only: entropy_gamma1_from_nH_p
        double precision, intent(in) :: log_nH, log_p_nH
        if (eos%type_id /= EOS_TYPE_LTE) call mpistop("gamma1_from_nH_p called outside its eos_type (LTE)")
        if (eos%method_id == EOS_ENTROPY) then
            ! Single bicubic Hermite lookup on the pre-built Gamma_1(rho, p) table.
            ! ZERO runtime iterations; no p->eint intermediate inversion.
            g1 = entropy_gamma1_from_nH_p(eos%g1p, eos%g1p_x, eos%g1p_y, &
                eos%g1p_xy, log_nH, log_p_nH)
        else
            g1 = bicubic_lookup(log_nH, log_p_nH, eos%gamma1_p)
        end if
    end function gamma1_from_nH_p

    !> Merged log10(p/nH) lookup: (log10 nH, log10 eint/nH) -> log10(p/nH)
    !> Single PCHIP evaluation replacing separate T + neOnH lookups.
    double precision function log_p_from_nH_eint(log_nH, log_eint_nH) result(lp)
        double precision, intent(in) :: log_nH, log_eint_nH
        lp = bicubic_lookup(log_nH, log_eint_nH, eos%log_p)
    end function log_p_from_nH_eint

    !> p/nH from (log10 nH, log10 eint/nH) in code units.
    !> Returns (1+He+y)*T directly -- single lookup replaces T + y lookups.
    double precision function p_nH_from_eint(log_nH, log_eint_nH) result(p_nH)
        use mod_eos_LTE_entropy, only: entropy_p_nH_from_eint
        double precision, intent(in) :: log_nH, log_eint_nH
        double precision, parameter :: ln10 = 2.302585092994046d0
        double precision :: T_loc, y_loc

        if (eos%type_id /= EOS_TYPE_LTE) call mpistop("p_nH_from_eint called outside its eos_type (LTE)")
        if (eos%method_id == EOS_ANALYTIC) then
            call saha_T_from_nH_eint(10.0d0**log_nH, 10.0d0**log_eint_nH, T_loc, y_loc)
            p_nH = (1.0d0 + eos%He_abundance + y_loc) * T_loc
        else if (eos%method_id == EOS_ENTROPY) then
            p_nH = entropy_p_nH_from_eint(eos%pfwd, eos%pfwd_x, eos%pfwd_y, &
                                           eos%pfwd_xy, log_nH, log_eint_nH)
        else
            p_nH = dexp(ln10 * bicubic_lookup(log_nH, log_eint_nH, eos%p_over_nH))
        end if
    end function p_nH_from_eint

    !> Internal energy per nH from (log10 nH, log10 T) in code units.
    !> Uses the bisection-built inverse table (H+He, machine precision).
    !> Fallback: H-only Saha if table not built.
    double precision function eint_nH_from_T(log_nH, log_T) result(eint_nH)
        use mod_eos_LTE_entropy, only: entropy_eint_from_nH_T
        double precision, intent(in) :: log_nH, log_T
        double precision, parameter :: ln10 = 2.302585092994046d0
        double precision :: log_e_nh

        if (eos%type_id /= EOS_TYPE_LTE) call mpistop("eint_nH_from_T called outside its eos_type (LTE)")
        if (eos%method_id == EOS_ENTROPY) then
            ! Single bicubic Hermite lookup on the pre-built eint(rho, T)
            ! inverse table. ZERO runtime iterations.
            log_e_nh = entropy_eint_from_nH_T(eos%eintT, eos%eintT_x, &
                eos%eintT_y, eos%eintT_xy, log_nH, log_T)
            eint_nH = 10.0d0**log_e_nh
        else if (allocated(eos%eint_from_T%table)) then
            eint_nH = dexp(ln10 * bicubic_lookup(log_nH, log_T, eos%eint_from_T))
        else
            eint_nH = saha_eint_from_nH_T(10.0d0**log_nH, 10.0d0**log_T)
        end if
    end function eint_nH_from_T

    !>#######################################################################
    !> Iterative solvers: table-guessed bisection for the WB pressure->eint
    !> inversion since the method needs high accuracy to not self-seed perturbations.
    !>#######################################################################

    !> Cached bisection on the log_p table for WB p->eint inversion.
    !> Precomputes nH-direction indices and table values once, then
    !> bisects using only cheap PCHIP evaluations with varying tx.
    !> Expects a narrow initial bracket [lo, hi] (e.g. from p2eint guess).
    !>
    !> Adaptive-grid support: when eos%log_p%is_uniform is .false., the
    !> y- and x-direction index calculations switch to binary search on
    !> the explicit node arrays. The cached 4x4 stencil mechanism is
    !> unchanged; only the index-to-cell mapping changes.
    subroutine log_p_bisect_cached(log_nH, log_p_target, &
        log_eint_lo, log_eint_hi, max_iter, log_eint_result)
        use mod_lookup_table, only: find_index_bsearch
        double precision, intent(in)    :: log_nH, log_p_target
        double precision, intent(inout) :: log_eint_lo, log_eint_hi
        integer, intent(in)             :: max_iter
        double precision, intent(out)   :: log_eint_result

        integer :: nx, ny, ix, iy, iter, ii
        integer :: i0, i1, i2, i3, j0, j1, j2, j3
        double precision :: tx, ty, rx, ry
        double precision :: xstep_inv, ystep_inv
        double precision :: vmin_x, vmax_x, vmin_y, vmax_y
        double precision :: log_eint_mid, log_p_eval
        logical :: is_unif

        !> Table values: tv(y_row, x_col) for 4 y-rows x 4 x-cols
        double precision :: tv(4,4)
        !> Cached ix for detecting grid cell change
        integer :: ix_cached

        nx = eos%log_p%dim2   ! eint axis (varx)
        ny = eos%log_p%dim1   ! nH axis (vary)
        is_unif = eos%log_p%is_uniform

        if (is_unif) then
            vmin_x = eos%log_p%var2_min
            vmax_x = eos%log_p%var2_max
            vmin_y = eos%log_p%var1_min
            vmax_y = eos%log_p%var1_max
            xstep_inv = dble(nx-1) / (vmax_x - vmin_x)
            ystep_inv = dble(ny-1) / (vmax_y - vmin_y)
        end if

        !> Precompute nH indices (constant across all iterations)
        if (is_unif) then
            ry = (log_nH - vmin_y) * ystep_inv
            ry = max(0.0d0, min(ry, dble(ny-1)))
            iy = int(ry)
            ty = ry - dble(iy)
        else
            !> Non-uniform: binary search on var1_nodes (length ny = dim1)
            if (log_nH <= eos%log_p%var1_nodes(1)) then
                iy = 0; ty = 0.0d0
            else if (log_nH >= eos%log_p%var1_nodes(ny)) then
                iy = ny - 2; ty = 1.0d0
            else
                ii = find_index_guard(eos%log_p%var1_nodes, ny, log_nH, &
                    eos%log_p%guard_1, eos%log_p%guard_M_1, eos%log_p%guard_scale_1)
                iy = max(0, min(ii - 2, ny - 2))
                ty = (log_nH - eos%log_p%var1_nodes(iy+1)) &
                   / (eos%log_p%var1_nodes(iy+2) - eos%log_p%var1_nodes(iy+1))
                ty = max(0.0d0, min(ty, 1.0d0))
            end if
        end if

        !> Clamped y-row indices
        j0 = max(0, min(ny-1, iy-1))
        j1 = max(0, min(ny-1, iy  ))
        j2 = max(0, min(ny-1, iy+1))
        j3 = max(0, min(ny-1, iy+2))

        ix_cached = -1

        do iter = 1, max_iter
            log_eint_mid = 0.5d0 * (log_eint_lo + log_eint_hi)

            !> Compute eint-direction index
            if (is_unif) then
                rx = (log_eint_mid - vmin_x) * xstep_inv
                rx = max(0.0d0, min(rx, dble(nx-1)))
                ix = int(rx)
                tx = rx - dble(ix)
            else
                !> Non-uniform: binary search on var2_nodes (length nx = dim2)
                if (log_eint_mid <= eos%log_p%var2_nodes(1)) then
                    ix = 0; tx = 0.0d0
                else if (log_eint_mid >= eos%log_p%var2_nodes(nx)) then
                    ix = nx - 2; tx = 1.0d0
                else
                    ii = find_index_guard(eos%log_p%var2_nodes, nx, log_eint_mid, &
                        eos%log_p%guard_2, eos%log_p%guard_M_2, eos%log_p%guard_scale_2)
                    ix = max(0, min(ii - 2, nx - 2))
                    tx = (log_eint_mid - eos%log_p%var2_nodes(ix+1)) &
                       / (eos%log_p%var2_nodes(ix+2) - eos%log_p%var2_nodes(ix+1))
                    tx = max(0.0d0, min(tx, 1.0d0))
                end if
            end if

            !> Load 16 table values only when ix changes
            if (ix /= ix_cached) then
                i0 = max(0, min(nx-1, ix-1))
                i1 = max(0, min(nx-1, ix  ))
                i2 = max(0, min(nx-1, ix+1))
                i3 = max(0, min(nx-1, ix+2))

                tv(1,1) = eos%log_p%table(j0+1, i0+1)
                tv(1,2) = eos%log_p%table(j0+1, i1+1)
                tv(1,3) = eos%log_p%table(j0+1, i2+1)
                tv(1,4) = eos%log_p%table(j0+1, i3+1)
                tv(2,1) = eos%log_p%table(j1+1, i0+1)
                tv(2,2) = eos%log_p%table(j1+1, i1+1)
                tv(2,3) = eos%log_p%table(j1+1, i2+1)
                tv(2,4) = eos%log_p%table(j1+1, i3+1)
                tv(3,1) = eos%log_p%table(j2+1, i0+1)
                tv(3,2) = eos%log_p%table(j2+1, i1+1)
                tv(3,3) = eos%log_p%table(j2+1, i2+1)
                tv(3,4) = eos%log_p%table(j2+1, i3+1)
                tv(4,1) = eos%log_p%table(j3+1, i0+1)
                tv(4,2) = eos%log_p%table(j3+1, i1+1)
                tv(4,3) = eos%log_p%table(j3+1, i2+1)
                tv(4,4) = eos%log_p%table(j3+1, i3+1)

                ix_cached = ix
            end if

            !> Evaluate PCHIP: 4 x-rows then 1 y-interp
            log_p_eval = pchip_2d_from_cache(tv, tx, ty)

            if (log_p_eval < log_p_target) then
                log_eint_lo = log_eint_mid
            else
                log_eint_hi = log_eint_mid
            end if

            if (dabs(log_eint_hi - log_eint_lo) < 1.0d-14) exit
        end do

        log_eint_result = 0.5d0 * (log_eint_lo + log_eint_hi)

    contains

        pure double precision function pchip_1d(p0, p1, p2, p3, t) result(v)
            double precision, intent(in) :: p0, p1, p2, p3, t
            double precision :: d0, d1, d2, m1, m2, s, a1, a2, lim
            double precision :: tt, ttt, h00, h10, h01, h11

            d0 = p1 - p0;  d1 = p2 - p1;  d2 = p3 - p2

            if (d1 == 0.0d0) then
                m1 = 0.0d0;  m2 = 0.0d0
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
                a1 = s*m1;  a2 = s*m2
                if (a1 < 0.0d0) a1 = 0.0d0
                if (a2 < 0.0d0) a2 = 0.0d0
                lim = 3.0d0*abs(d1)
                if (a1 > lim) a1 = lim
                if (a2 > lim) a2 = lim
                m1 = s*a1;  m2 = s*a2
            end if

            tt = t*t;  ttt = tt*t
            h00 = 2.0d0*ttt - 3.0d0*tt + 1.0d0
            h10 = ttt - 2.0d0*tt + t
            h01 = -2.0d0*ttt + 3.0d0*tt
            h11 = ttt - tt
            v = h00*p1 + h10*m1 + h01*p2 + h11*m2
        end function pchip_1d

        pure double precision function pchip_2d_from_cache(c, tx, ty) result(z)
            double precision, intent(in) :: c(4,4), tx, ty
            double precision :: g0, g1, g2, g3
            !> 4 x-direction interpolations (one per y-row)
            g0 = pchip_1d(c(1,1), c(1,2), c(1,3), c(1,4), tx)
            g1 = pchip_1d(c(2,1), c(2,2), c(2,3), c(2,4), tx)
            g2 = pchip_1d(c(3,1), c(3,2), c(3,3), c(3,4), tx)
            g3 = pchip_1d(c(4,1), c(4,2), c(4,3), c(4,4), tx)
            !> 1 y-direction interpolation
            z = pchip_1d(g0, g1, g2, g3, ty)
        end function pchip_2d_from_cache

    end subroutine log_p_bisect_cached

    !> Given log10(nH) and log10(p), find log10(eint/nH) by table-guessed
    !> bisection on the forward pressure table.  Wraps the bracketing logic
    !> and log_p_bisect_cached into a single call for use by both HD and MHD
    !> from_conserved routines.
    subroutine eint_from_p_bisect(log_nH_val, log_p_val, log_eint_nH_out)
        double precision, intent(in)  :: log_nH_val, log_p_val
        double precision, intent(out) :: log_eint_nH_out

        double precision :: log_p_target, p2eint_ratio, log_eint_guess
        double precision :: log_p_at_guess, log_eint_lo, log_eint_hi
        double precision :: f_bracket, margin
        double precision :: eint_lim_lo, eint_lim_hi
        integer :: max_iter

        if (eos%type_id /= EOS_TYPE_LTE) call mpistop("eint_from_p_bisect called outside its eos_type (LTE)")
        log_p_target = log_p_val - log_nH_val

        !> eint-axis limits. var2_min/var2_max are populated for uniform tables only; with an
        !> adaptive (non-uniform) log_p table they are meaningless, so the guess/bracket clamps
        !> below let log_eint run out of range -> out-of-bounds stencil load in bicubic_lookup
        !> (SIGSEGV at init). log_p_bisect_cached already branches on is_uniform; this wrapper did not.
        if (eos%log_p%is_uniform) then
            eint_lim_lo = eos%log_p%var2_min
            eint_lim_hi = eos%log_p%var2_max
        else
            eint_lim_lo = eos%log_p%var2_nodes(1)
            eint_lim_hi = eos%log_p%var2_nodes(eos%log_p%dim2)
        end if

        ! Initial guess from p2eint table
        p2eint_ratio = p2eint_from_nH_p(log_nH_val, log_p_target)
        log_eint_guess = dlog10(p2eint_ratio) + log_p_target

        log_eint_guess = max(log_eint_guess, eint_lim_lo)
        log_eint_guess = min(log_eint_guess, eint_lim_hi)

        log_p_at_guess = log_p_from_nH_eint(log_nH_val, log_eint_guess)

        ! Establish bracket around the guess
        margin = 5.0d-4
        if (log_p_at_guess < log_p_target) then
            log_eint_lo = log_eint_guess
            log_eint_hi = min(log_eint_guess + margin, eint_lim_hi)
            f_bracket = log_p_from_nH_eint(log_nH_val, log_eint_hi) - log_p_target
        else
            log_eint_lo = max(log_eint_guess - margin, eint_lim_lo)
            log_eint_hi = log_eint_guess
            f_bracket = -(log_p_from_nH_eint(log_nH_val, log_eint_lo) - log_p_target)
        end if

        if (f_bracket >= 0.0d0) then
            max_iter = 8
        else
            log_eint_lo = eint_lim_lo
            log_eint_hi = eint_lim_hi
            max_iter = 20
        end if

        call log_p_bisect_cached(log_nH_val, log_p_target, &
            log_eint_lo, log_eint_hi, max_iter, log_eint_nH_out)

    end subroutine eint_from_p_bisect

    !>#######################################################################
    !> Table-metadata helpers: expose inverse-table grid bounds without
    !> reaching into EoS table internals.
    !>#######################################################################

    !> log_nH grid metadata of the (log_nH, log_T) inverse table (eint from T),
    !> choosing the container by method ('tables'->eint_from_T, 'entropy'->eintT).
    !> n_nH=0 if no such table (analytic/FI). Lets cooling's build_Y_mod_table
    !> get its grid without reaching into EoS table internals.
    subroutine eos_get_eintT_grid(n_nH, lg_nH_min, lg_nH_max)
        integer, intent(out)          :: n_nH
        double precision, intent(out) :: lg_nH_min, lg_nH_max

        if (eos%type_id /= EOS_TYPE_LTE) call mpistop("eos_get_eintT_grid called outside its eos_type (LTE)")
        if (allocated(eos%eint_from_T%table)) then
            n_nH = eos%eint_from_T%dim1
            lg_nH_min = eos%eint_from_T%var1_min
            lg_nH_max = eos%eint_from_T%var1_max
        else if (allocated(eos%eintT%table)) then
            n_nH = eos%eintT%dim1
            lg_nH_min = eos%eintT%var1_min
            lg_nH_max = eos%eintT%var1_max
        else
            n_nH = 0
            lg_nH_min = 0.0d0
            lg_nH_max = 0.0d0
        end if
    end subroutine eos_get_eintT_grid

end module mod_eos_LTE
