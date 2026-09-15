!=============================================================================
!> FFHD <-> EoS seam: binds the eos% authority into force-free hydrodynamics.
!>
!> Mirrors mod_hd_eos.t, adapted to FFHD's structure: a single field-aligned
!> momentum component mom(1) and NO magnetic energy (force-free, frozen B0), so
!> every kinetic-energy term is half*rho*mom(1)**2. The B0 projection lives in
!> the flux/cmax/velocity layer (mod_ffhd_phys.t) and is orthogonal to the EoS
!> seam, so the thermodynamic conversions are HD-like with one momentum
!> component.
!>
!> ffhd_link_eos        -- wires eos%/phys_ pointers per eos_type (FI now;
!>                         LTE/PI guarded until milestones 2/3).
!> ffhd_bind_eos_to_source -- wires the thermal-conduction (tc_fl), radiative-
!>                         cooling (rc_fl) and thermal-emission (te_fl_ffhd)
!>                         fluid-port callbacks from eos%. This is what the
!>                         pre-eos% hand-wiring was missing (rc_fl%get_Te /
!>                         get_ne_nH), which caused the cooling SIGSEGV.
!=============================================================================
module mod_ffhd_eos
    use mod_global_parameters
    use mod_physics
    use mod_eos
    !> Mode-specific kernels come from their sub-modules (the facade no longer
    !> re-exports them); each mpistops if called under the wrong eos_type/method.
    use mod_eos_lte
    use mod_eos_lte_saha
    use mod_eos_pi
    use mod_eos_container
    use mod_ffhd_phys
    use mod_radiative_cooling, only: build_Y_mod_table

    use mod_comm_lib, only: mpistop

    implicit none
    private

    public :: ffhd_link_eos

contains

    !> Link the appropriate EoS conversion/closure routines per eos_type.
    !> Called from ffhd_activate after ffhd_phys_init (mirrors hd_activate).
    subroutine ffhd_link_eos()
        use mod_usr_methods

        ! Isothermal FFHD (no energy equation) has no thermodynamic state for
        ! the EoS to close -- only FI is meaningful there.
        if (.not. ffhd_energy .and. eos%eos_type /= 'FI') &
            call mpistop("FFHD: ffhd_energy=.false. requires eos_type='FI'")

        ! FI / PI base conversions = ffhd's origin routines (single-component
        ! KE). PI energy mode overrides them below (milestone 3). LTE routes
        ! primitive<->conserved through ffhd_p_to_e / table inversion.
        if (eos%eos_type == 'FI' .or. eos%eos_type == 'PI') then
            eos%to_conserved => ffhd_to_conserved_origin
            eos%to_primitive => ffhd_to_primitive_origin
        else if (eos%eos_type == 'LTE') then
            eos%to_conserved => ffhd_to_conserved_LTE
            eos%to_primitive => ffhd_to_primitive_LTE
        else
            call mpistop("Unknown FFHD EOS type: " // trim(eos%eos_type))
        end if

        phys_to_primitive => eos%to_primitive
        phys_to_conserved => eos%to_conserved
        ffhd_to_primitive => eos%to_primitive
        ffhd_to_conserved => eos%to_conserved
        phys_get_rho      => eos%get_rho

        ! suitable for both FI (origin inlines it) and LTE (conversions call it)
        eos%p_to_e => ffhd_p_to_e

        ! pthermal: the isothermal (no-energy) pointer is wired in phys_init and
        ! left untouched; for the energy equation route through eos%. LTE reads
        ! the stored Te_/Ne_ aux state; FI/PI use (gamma-1)*(e-KE).
        if (ffhd_energy) then
            if (eos%eos_type == 'LTE') then
                phys_get_pthermal        => ffhd_get_pthermal_LTE
                ffhd_get_pthermal        => ffhd_get_pthermal_LTE
                eos%get_thermal_pressure => ffhd_get_pthermal_LTE
            else
                phys_get_pthermal        => ffhd_get_pthermal_origin
                ffhd_get_pthermal        => ffhd_get_pthermal_origin
                eos%get_thermal_pressure => ffhd_get_pthermal_origin
            end if
        end if

        phys_bind_eos_to_source => ffhd_bind_eos_to_source

        ! Energy<->internal heads: ffhd previously handed these only to the TC
        ! STS via set_conversion_methods_to_head; the LTE fast-T kernel and
        ! phys_update_temperature also need the module-level phys_ pointers.
        phys_e_to_ei => ffhd_e_to_ei
        phys_ei_to_e => ffhd_ei_to_e
        phys_get_ei  => ffhd_get_ei   ! LTE+ionE cooling recovers eint via this

        ! sound speed / gamma1. LTE+ionE: Gamma_1(nH,p) from the EoS tables and
        ! EoS-aware AMR prolongation (interpolate in (rho,v,T) to avoid Jensen's
        ! inequality across the ionisation plateau). Otherwise FI: cs2=gamma*p/rho.
        if (eos%eos_type == 'LTE' .and. eos%ionE) then
            eos%get_csound2   => ffhd_get_csound2_LTE
            phys_get_gamma1   => ffhd_get_gamma1_LTE
            phys_to_prolong   => ffhd_to_prolong_LTE
            phys_from_prolong => ffhd_from_prolong_LTE
        else
            eos%get_csound2 => ffhd_get_csound2
            phys_get_gamma1 => get_gamma1_FI
        end if

        ! ffhd_get_temperature (phys pointer) by eos_type: PI reads the cached
        ! Te_, the others recompute from etot.
        if (eos%eos_type == 'PI') then
            ffhd_get_temperature => ffhd_get_temperature_from_Te
        else
            ffhd_get_temperature => ffhd_get_temperature_from_etot
        end if
        ! Rfactor: only the FI case is physics-dependent (usr_Rfactor). LTE and PI
        ! bind pure eos% routines in eos_finalise_{LTE,PI} (run after this, before
        ! ffhd_bind_eos_to_source -> targets stay private). usr_Rfactor thus applies
        ! to FI only, consistent with mhd/hd.
        if (associated(usr_Rfactor)) then
            eos%get_Rfactor => usr_Rfactor
        else
            eos%get_Rfactor => Rfactor_from_constant_ionization
        end if
        ! ffhd_get_Rfactor (runtime alias) captured in ffhd_bind_eos_to_source,
        ! AFTER eos_finalise has set eos%get_Rfactor for LTE/PI.

        !> PI energy-mode overrides (eos_type='PI', ionE=.true.). Energy mode
        !> differs from no-energy PI only in the eint<->p relation (eint carries
        !> the ionisation energy), so override exactly the routines touching it,
        !> all via the portable mod_eos_pi backend. get_Rfactor stays the
        !> Te_-addressed routine above (Te_ refreshed each step by
        !> ffhd_update_temperature_PI).
        if (eos%eos_type == 'PI' .and. eos%ionE) then
            if (.not. ffhd_energy) &
                call mpistop('FFHD PI energy EoS requires ffhd_energy=.true.')
            eos%to_conserved  => ffhd_to_conserved_PI
            eos%to_primitive  => ffhd_to_primitive_PI
            phys_to_primitive => eos%to_primitive
            phys_to_conserved => eos%to_conserved
            ffhd_to_primitive => eos%to_primitive
            ffhd_to_conserved => eos%to_conserved
            eos%p_to_e        => ffhd_p_to_e_PI
            phys_get_pthermal        => ffhd_get_pthermal_PI
            ffhd_get_pthermal        => ffhd_get_pthermal_PI
            eos%get_thermal_pressure => ffhd_get_pthermal_PI
            !> eos%get_csound2 => get_csound2_PI set in eos_finalise_PI (private target)
            phys_get_gamma1   => get_gamma1_PI
        end if

    end subroutine ffhd_link_eos

    !> Wire the source-term fluid ports from eos%. Called from eos_finalise via
    !> phys_bind_eos_to_source, AFTER eos_finalise has set eos%get_Te,
    !> get_ne_nH, eion_per_nH, n_per_nH_FI/neOnH_FI, get_temperature_from_*.
    subroutine ffhd_bind_eos_to_source()

        !> Runtime R-factor alias: eos%get_Rfactor is fully set by now
        !> (eos_finalise ran before phys_bind_eos_to_source), incl. the LTE/PI
        !> bindings done in eos_finalise_{LTE,PI}.
        ffhd_get_Rfactor => eos%get_Rfactor

        if (allocated(tc_fl)) then
            tc_fl%get_temperature_from_conserved => eos%get_temperature_from_etot
            if (eos%eos_type == 'LTE' .and. eos%ionE) then
                tc_fl%get_temperature_from_eint => get_temperature_from_eint_fast_LTE
            else
                tc_fl%get_temperature_from_eint => eos%get_temperature_from_eint
            end if
            tc_fl%get_rho           => eos%get_rho
            tc_fl%get_ne_nH         => eos%get_ne_nH
            tc_fl%get_var_Rfactor   => eos%get_Rfactor
            tc_fl%inv_gamma_minus_1 =  eos%inv_gamma_minus_1
            tc_fl%nH2rhoFactor      =  eos%nH2rhoFactor
            tc_fl%log_T_floor       =  eos_get_log_T_floor()
            tc_fl%eint_from_T       => eint_nH_from_T
        end if

        if (allocated(rc_fl)) then
            rc_fl%get_rho           => eos%get_rho
            rc_fl%get_pthermal      => eos%get_thermal_pressure
            rc_fl%get_var_Rfactor   => eos%get_Rfactor
            rc_fl%get_Te            => eos%get_Te
            rc_fl%get_ne_nH         => eos%get_ne_nH
            rc_fl%ionE              =  eos%ionE
            rc_fl%method            =  eos%method
            rc_fl%inv_gamma_minus_1 =  eos%inv_gamma_minus_1
            rc_fl%nH2rhoFactor      =  eos%nH2rhoFactor
            rc_fl%eion_per_nH       =  eos%eion_per_nH
            rc_fl%eint_from_T       => eint_nH_from_T
            rc_fl%p2eint            => p2eint_from_nH_p
            rc_fl%T_from_eint       => T_from_nH_eint
            rc_fl%y_from_eint       => y_from_nH_eint
            if (eos%ionE .and. eos%eos_type == 'LTE') then
                call eos_get_eintT_grid(rc_fl%Y_mod_n_nH, &
                     rc_fl%Y_mod_lg_nH_min, rc_fl%Y_mod_lg_nH_max)
                call build_Y_mod_table(rc_fl)
            end if
        end if

        !> PI energy-mode: repoint eint<->p kernels at the portable backend.
        if (eos%eos_type == 'PI' .and. eos%ionE) then
            if (allocated(tc_fl)) tc_fl%eint_from_T => eint_from_T_PI
            if (allocated(rc_fl)) then
                rc_fl%eint_from_T => eint_from_T_PI
                rc_fl%p2eint      => p2eint_PI
                rc_fl%T_from_eint => T_from_eint_PI
                rc_fl%y_from_eint => y_from_eint_PI
            end if
        end if

        {^IFTHREED
        if (allocated(te_fl_ffhd)) then
            te_fl_ffhd%get_rho         => eos%get_rho
            te_fl_ffhd%get_pthermal    => eos%get_thermal_pressure
            te_fl_ffhd%get_var_Rfactor => eos%get_Rfactor
            te_fl_ffhd%get_ne_nH       => eos%get_ne_nH
        end if
        }

    end subroutine ffhd_bind_eos_to_source

    !=========================================================================
    !> LTE (Saha table) thermodynamics for FFHD. Mirrors mod_hd_eos.t adapted
    !> to FFHD's single field-aligned momentum mom(1) and total-energy-only
    !> ("origin") formulation; no dust, no magnetic energy.
    !=========================================================================

    !> LTE primitive -> conserved. On entry rho/v/p; on exit rho/mom/E.
    subroutine ffhd_to_conserved_LTE(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        call ffhd_p_to_e(ixI^L, ixO^L, w, x)
        w(ixO^S,mom(1)) = w(ixO^S,rho_)*w(ixO^S,mom(1))
    end subroutine ffhd_to_conserved_LTE

    !> Convert pressure to total energy: E = eint(rho, p) + KE.
    !> On entry rho/v/p; on exit w(e_) = total energy (rho, mom unchanged).
    !> Four ionE paths mirror hd p_to_e: FI bypass, analytic Saha, p-bisect,
    !> p2eint table. Non-ionE LTE uses the constant inverse-gamma relation.
    subroutine ffhd_p_to_e(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        integer :: ix^D
        double precision :: p_to_eint, p_over_rho
        double precision :: nH(ixI^S), nH_in(ixI^S), p_in(ixI^S)
        double precision :: log_eint_mid, eint_total
        double precision :: T_solve, y_solve, eint_nH_solve

        if (.not. ffhd_energy) return

        if (eos%ionE) then
            call eos%get_nH(w, x, ixI^L, ixO^L, nH)
            nH_in(ixO^S) = dlog10(nH(ixO^S))
            if (eos%p2eint_method /= 'bisect') then
                p_in(ixO^S) = dlog10(w(ixO^S,p_)) - nH_in(ixO^S)
            end if
        endif

        p_to_eint = eos%inv_gamma_minus_1
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if (eos%ionE) then
                p_over_rho = w(ix^D,p_) / w(ix^D,rho_)
                if (p_over_rho > eos%p_rho_FI_threshold) then
                    p_to_eint = eos%inv_gamma_minus_1 &
                        + eos%eion_per_nH * nH(ix^D) / w(ix^D,p_)
                    w(ix^D,e_) = w(ix^D,p_)*p_to_eint + &
                        half*w(ix^D,mom(1))**2*w(ix^D,rho_)
                else if (eos%method == 'analytic') then
                    call saha_state_from_nH_p(nH(ix^D), w(ix^D,p_), &
                        T_solve, y_solve, eint_nH_solve)
                    w(ix^D,e_) = eint_nH_solve * nH(ix^D) + &
                        half*w(ix^D,mom(1))**2*w(ix^D,rho_)
                else if (eos%p2eint_method == 'bisect') then
                    call eint_from_p_bisect(nH_in(ix^D), &
                        dlog10(w(ix^D,p_)), log_eint_mid)
                    eint_total = nH(ix^D) * 10.0d0**log_eint_mid
                    eint_total = max(eint_total, &
                        nH(ix^D) * 10.0d0**eos%T%var2_min)
                    w(ix^D,e_) = eint_total + &
                        half*w(ix^D,mom(1))**2*w(ix^D,rho_)
                else
                    p_to_eint = p2eint_from_nH_p(nH_in(ix^D), p_in(ix^D))
                    w(ix^D,e_) = w(ix^D,p_)*p_to_eint + &
                        half*w(ix^D,mom(1))**2*w(ix^D,rho_)
                end if
            else
                w(ix^D,e_) = w(ix^D,p_)*p_to_eint + &
                    half*w(ix^D,mom(1))**2*w(ix^D,rho_)
            end if
        {end do\}
    end subroutine ffhd_p_to_e

    !> LTE conserved -> primitive. On entry rho/mom/E; on exit rho/v/p.
    !> Pressure is energy-consistent from actual eint (stored Ne_/Te_ may be
    !> stale after AMR prolong/coarsen on the nonlinear EoS).
    subroutine ffhd_to_primitive_LTE(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision                :: inv_rho
        double precision                :: nH(ixI^S), log_nH(ixI^S)
        double precision                :: eint_val, eint_in
        integer :: ix^D

        if (fix_small_values) then
            call ffhd_handle_small_values(.false., w, x, ixI^L, ixO^L, &
                'ffhd_to_primitive_LTE')
        end if

        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        if (eos%ionE) log_nH(ixO^S) = dlog10(nH(ixO^S))

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0/w(ix^D,rho_)
            w(ix^D,mom(1)) = w(ix^D,mom(1))*inv_rho
            if (ffhd_energy) then
                if (eos%ionE) then
                    eint_val = w(ix^D,e_) - half*w(ix^D,rho_)*w(ix^D,mom(1))**2
                    if (eos%method /= 'analytic') then
                        eint_val = max(eint_val, nH(ix^D)*10.0d0**eos%T%var2_min)
                    end if
                    eint_val = max(eint_val, smalldouble)
                    if (eint_val * inv_rho > eos%eint_rho_FI_threshold) then
                        w(ix^D,p_) = eos%gamma_minus_1 &
                            * (eint_val - eos%eion_per_nH * nH(ix^D))
                    else
                        eint_in = dlog10(eint_val) - log_nH(ix^D)
                        w(ix^D,p_) = nH(ix^D) &
                            * p_nH_from_eint(log_nH(ix^D), eint_in)
                    end if
                else
                    w(ix^D,p_) = eos%gamma_minus_1 &
                        * (w(ix^D,e_) - half*w(ix^D,rho_)*w(ix^D,mom(1))**2)
                end if
            end if
        {end do\}
    end subroutine ffhd_to_primitive_LTE

    !> LTE thermal pressure from the stored aux state: p = nH(1+He+ne/nH)*Te.
    !> Mirrors the LTE branch of hd_get_pthermal (small-value guards from
    !> ffhd_get_pthermal_origin).
    subroutine ffhd_get_pthermal_LTE(w, x, ixI^L, ixO^L, pth)
        use mod_global_parameters
        use mod_small_values, only: trace_small_values
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: pth(ixI^S)
        integer                      :: iw, ix^D
        double precision             :: nH(ixI^S)

        if (ffhd_energy) then
            call eos%get_nH(w, x, ixI^L, ixO^L, nH)
            pth(ixO^S) = nH(ixO^S) * (1.0d0 + eos%He_abundance &
                + (w(ixO^S,iw_ne) / nH(ixO^S))) * w(ixO^S,iw_te)
        else
            pth(ixO^S) = ffhd_adiab * w(ixO^S, rho_)**eos%gamma
        end if

        if (fix_small_values) then
            {do ix^DB= ixO^LIM^DB\}
               if(pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
            {end do^D&\}
        else if (check_small_values) then
            {do ix^DB= ixO^LIM^DB\}
               if(pth(ix^D)<small_pressure) then
                 write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                      " encountered when call ffhd_get_pthermal_LTE"
                 write(*,*) "Iteration: ", it, " Time: ", global_time
                 write(*,*) "Location: ", x(ix^D,:)
                 write(*,*) "Cell number: ", ix^D
                 do iw=1,nw
                   write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
                 end do
                 if(trace_small_values) write(*,*) dsqrt(pth(ix^D)-bigdouble)
                 write(*,*) "Saving status at the previous time step"
                 crash=.true.
               end if
            {end do^D&\}
        end if
    end subroutine ffhd_get_pthermal_LTE

    !> cs2 = Gamma_1(LTE) * p/rho.
    subroutine ffhd_get_csound2_LTE(w, x, ixI^L, ixO^L, cs2)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision, intent(out)   :: cs2(ixI^S)
        integer :: ix^D

        call ffhd_get_gamma1_LTE(w, x, ixI^L, ixO^L, cs2)
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cs2(ix^D) = cs2(ix^D) * w(ix^D, p_) / w(ix^D, rho_)
        {end do\}
    end subroutine ffhd_get_csound2_LTE

    !> Effective Gamma_1 for LTE+IonE (dispatch on gamma1_method / eos%method).
    subroutine ffhd_get_gamma1_LTE(w, x, ixI^L, ixO^L, gamma1)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision, intent(out)   :: gamma1(ixI^S)
        double precision :: nH_val, p_over_rho
        integer :: ix^D

        if (eos%gamma1_method == 'constant') then
            gamma1(ixO^S) = eos%gamma
            return
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            p_over_rho = w(ix^D, p_) / w(ix^D, rho_)
            if (p_over_rho > eos%p_rho_FI_threshold) then
                gamma1(ix^D) = eos%gamma
            else
                nH_val = w(ix^D, rho_) / eos%nH2rhoFactor
                if (eos%method == 'analytic') then
                    if (iw_te > 0 .and. w(ix^D,iw_te) > 0.0d0) then
                        gamma1(ix^D) = saha_gamma1_from_nH_T(nH_val, w(ix^D,iw_te))
                    else
                        gamma1(ix^D) = eos%gamma
                    end if
                else
                    gamma1(ix^D) = gamma1_from_nH_p(dlog10(nH_val), &
                        dlog10(w(ix^D, p_) / nH_val))
                end if
            end if
        {end do\}
    end subroutine ffhd_get_gamma1_LTE

    !> Conserved (rho, rho*v, E) -> prolong form (rho, v, T); T in the e_ slot.
    !> Interpolating in (rho, v, T) avoids Jensen's inequality across the
    !> ionisation plateau. Mirrors hd_to_prolong_LTE.
    subroutine ffhd_to_prolong_LTE(ixI^L, ixO^L, w, x)
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: inv_rho, eint_val, nH_val, log_nH, T_loc, y_loc
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0 / w(ix^D, rho_)
            nH_val = w(ix^D, rho_) / eos%nH2rhoFactor
            log_nH = dlog10(nH_val)

            w(ix^D,mom(1)) = w(ix^D,mom(1))*inv_rho

            eint_val = w(ix^D, e_) - half*w(ix^D, rho_)*w(ix^D,mom(1))**2
            if (eos%method /= 'analytic') then
                eint_val = max(eint_val, nH_val * 10.0d0**eos%T%var2_min)
            end if
            eint_val = max(eint_val, smalldouble)

            if (eint_val * inv_rho > eos%eint_rho_FI_threshold) then
                w(ix^D, e_) = eos%gamma_minus_1 &
                    * (eint_val - eos%eion_per_nH * nH_val) &
                    / (nH_val * eos%n_per_nH_FI)
            else if (eos%method == 'analytic') then
                call saha_T_from_nH_eint(nH_val, eint_val / nH_val, T_loc, y_loc)
                w(ix^D, e_) = T_loc
            else
                w(ix^D, e_) = T_from_nH_eint(log_nH, &
                    dlog10(eint_val) - log_nH)
            end if
        {end do\}
    end subroutine ffhd_to_prolong_LTE

    !> Prolong form (rho, v, T) -> conserved (rho, rho*v, E); T read from e_.
    !> Mirrors hd_from_prolong_LTE.
    subroutine ffhd_from_prolong_LTE(ixI^L, ixO^L, w, x)
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: T_val, eint_val, T_FI, nH_val, log_nH, log_T_min
        integer :: ix^D

        T_FI = (eos%eint_rho_FI_threshold &
            * eos%nH2rhoFactor - eos%eion_per_nH) &
            * eos%gamma_minus_1 / eos%n_per_nH_FI

        if (eos%method == 'entropy') then
            log_T_min = eos%eintT%var2_min
        else
            log_T_min = eos%eint_from_T%var2_min
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            T_val = w(ix^D, e_)
            nH_val = w(ix^D, rho_) / eos%nH2rhoFactor
            log_nH = dlog10(nH_val)

            if (T_val > T_FI) then
                eint_val = nH_val &
                    * (eos%n_per_nH_FI * T_val * eos%inv_gamma_minus_1 &
                    + eos%eion_per_nH)
            else if (eos%method == 'analytic') then
                eint_val = saha_eint_from_nH_T(nH_val, T_val) * nH_val
            else
                eint_val = eint_nH_from_T(log_nH, &
                    dlog10(max(T_val, 10.0d0**log_T_min))) * nH_val
            end if

            w(ix^D, e_) = eint_val + half*w(ix^D, rho_)*w(ix^D,mom(1))**2
            w(ix^D,mom(1)) = w(ix^D,rho_)*w(ix^D,mom(1))
        {end do\}
    end subroutine ffhd_from_prolong_LTE

    !=========================================================================
    !> PI energy-mode (eos_type='PI', ionE=.true.) thermodynamics for FFHD.
    !> eint carries the ionisation-energy term, so p=(gamma-1)*eint no longer
    !> holds; the eint<->p relation is delegated to the portable scalar backend
    !> (mod_eos_pi). Single momentum, total-energy only. Mirrors the
    !> hd PI-energy family.
    !=========================================================================

    !> Primitive pressure -> total energy via backend (mom(1) is velocity here).
    subroutine ffhd_p_to_e_PI(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision :: eint
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if (ffhd_energy) then
                call eint_from_rho_p_PI(w(ix^D,rho_), w(ix^D,p_), eint)
                w(ix^D,e_) = eint + half*w(ix^D,mom(1))**2*w(ix^D,rho_)
            end if
        {end do\}
    end subroutine ffhd_p_to_e_PI

    !> Primitive -> conserved (PI energy).
    subroutine ffhd_to_conserved_PI(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        call ffhd_p_to_e_PI(ixI^L, ixO^L, w, x)
        w(ixO^S,mom(1)) = w(ixO^S,rho_)*w(ixO^S,mom(1))
    end subroutine ffhd_to_conserved_PI

    !> Conserved -> primitive (PI energy): KE removed -> eint, backend -> p.
    subroutine ffhd_to_primitive_PI(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision :: inv_rho, eint_val, T, Rfac
        integer :: ix^D

        if (fix_small_values) then
            call ffhd_handle_small_values(.false., w, x, ixI^L, ixO^L, &
                'ffhd_to_primitive_PI')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0/w(ix^D,rho_)
            w(ix^D,mom(1)) = w(ix^D,mom(1))*inv_rho
            if (ffhd_energy) then
                eint_val = w(ix^D,e_) - half*w(ix^D,rho_)*w(ix^D,mom(1))**2
                eint_val = max(eint_val, smalldouble)
                call state_from_eint_PI(w(ix^D,rho_), eint_val, &
                    T, w(ix^D,p_), Rfac)
            end if
        {end do\}
    end subroutine ffhd_to_primitive_PI

    !> PI-energy thermal pressure: eint (= e-KE) -> p via the backend.
    subroutine ffhd_get_pthermal_PI(w, x, ixI^L, ixO^L, pth)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
        double precision, intent(out):: pth(ixI^S)
        double precision :: ei(ixI^S), Tpi, Rpi
        integer :: ix^D

        ei(ixO^S) = ffhd_get_ei(w, ixI^L, ixO^L)
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            call state_from_eint_PI(w(ix^D,rho_), ei(ix^D), Tpi, pth(ix^D), Rpi)
        {end do\}
    end subroutine ffhd_get_pthermal_PI

    !> Adiabatic sound speed^2 from primitive (rho, p) via the backend.

    !> Effective Gamma_1 = cs2 * rho / p (PI energy).

    !> PI-energy Te_ refresh: eint (= e-KE) -> T via the backend.

    !> Prominence (T,p table) R-factor from (rho, pth) -- no-energy only, so
    !> pth=(gamma-1)*eint and R is independent of it. Mirrors hd_Rfactor_prominence.

    !> Prominence (T,p) Te_ update from (rho, pth).

end module mod_ffhd_eos
