!=============================================================================
!> MHD <-> EoS seam: binds the eos% authority into magnetohydrodynamics.
!>
!> mhd_link_eos wires the eos%/phys_ procedure pointers (conversions, pthermal,
!> csound2/gamma1, get_Rfactor, temperature, update_temperature) per eos_type
!> (FI / LTE / PI, with PI energy- and prominence-mode overrides) and per energy
!> formulation (origin / internal_e / hydrodynamic_e / semirelativistic).
!> bind_eos_to_source wires the thermal-conduction / radiative-cooling / thermal-
!> emission / FLD fluid-port callbacks from eos%.
!>
!> Holds the MHD block routines that wrap the shared thermodynamics with MHD's
!> mechanical-energy bookkeeping (kinetic + MAGNETIC energy): the FI/LTE/PI
!> conversions, p_to_e, the pthermal variants, csound2/gamma1, and the PI energy
!> + prominence R-factor / temperature kernels.
!=============================================================================
module mod_mhd_eos
    use mod_global_parameters
    use mod_physics
    use mod_eos
    !> The mode-specific scalar kernels are reached from their sub-modules, not the
    !> mod_eos facade (which no longer re-exports them). Each kernel mpistops if
    !> called under the wrong eos_type/method.
    use mod_eos_LTE
    use mod_eos_LTE_saha
    use mod_eos_PI
    use mod_eos_container
    use mod_mhd_phys
    use mod_timing
    use mod_radiative_cooling, only: build_Y_mod_table

    use mod_comm_lib, only: mpistop

    implicit none
    private

    !> Thermal pressure pointer: set by mhd_link_eos based on energy formulation.
    !> Internal to mod_mhd_eos — external callers use eos%get_thermal_pressure.
    procedure(sub_get_pthermal), pointer, public :: mhd_get_pthermal  => null()

    !> Temperature pointer: set by mhd_link_eos based on EoS type and energy formulation.
    procedure(sub_get_pthermal), pointer, public :: mhd_get_temperature => null()
    !> use habitual name of converting to primitive
    procedure(sub_convert), pointer, public :: mhd_to_primitive  => null()
    !> use habitual name of converting to conserved
    procedure(sub_convert), pointer, public :: mhd_to_conserved  => null()

    public :: mhd_link_eos
contains

    !> Link the appropriate EOS conversion routines based on the selected EoS type
    subroutine mhd_link_eos()
        use mod_usr_methods

        !> Primitive <-> Conserved routing
        ! PI (partial approximations) takes the FI conversions as its BASE: in
        ! no-energy mode this is exact (e_int = p/(gamma-1); the variable R-factor
        ! enters only via eos%get_Rfactor, below). In energy mode (eos%ionE) the
        ! PI override block at the end of this routine repoints to_conserved/
        ! to_primitive (and p_to_e/pthermal/csound) at the PI-energy variants, so
        ! the FI assignment here is just a harmless base that gets replaced.
        if (eos%eos_type == 'FI' .or. eos%eos_type == 'PI') then
            ! PI takes its electron count from a variable R factor evaluated on
            ! the rho slot, which under splitting holds only the perturbation.
            if (eos%eos_type == 'PI' .and. has_equi_rho_and_p) &
              call mpistop('PI EoS not supported with equilibrium splitting')
            if(mhd_hydrodynamic_e) then
                eos%to_primitive        => mhd_to_primitive_hde
                eos%to_conserved        => mhd_to_conserved_hde
            else if(mhd_semirelativistic) then
                if(mhd_energy) then
                    eos%to_primitive        => mhd_to_primitive_semirelati
                    eos%to_conserved        => mhd_to_conserved_semirelati
                else
                    eos%to_primitive        => mhd_to_primitive_semirelati_noe
                    eos%to_conserved        => mhd_to_conserved_semirelati_noe
                end if
            else
                if(has_equi_rho_and_p) then
                    eos%to_primitive        => mhd_to_primitive_split_rho
                    eos%to_conserved        => mhd_to_conserved_split_rho
                else if(mhd_internal_e) then
                    eos%to_primitive        => mhd_to_primitive_inte
                    eos%to_conserved        => mhd_to_conserved_inte
                else if(mhd_energy) then
                    eos%to_primitive         => mhd_to_primitive_origin
                    eos%to_conserved         => mhd_to_conserved_origin
                else
                    eos%to_primitive         => mhd_to_primitive_origin_noe
                    eos%to_conserved         => mhd_to_conserved_origin_noe
                end if
            end if
        else if (eos%eos_type == 'LTE') then
            ! Guard unsupported combinations
            if (mhd_hydrodynamic_e) &
              call mpistop('LTE EoS not supported with mhd_hydrodynamic_e')
            if (mhd_semirelativistic) &
              call mpistop('LTE EoS not supported with mhd_semirelativistic')
            if (has_equi_rho_and_p) &
              call mpistop('LTE EoS not supported with equilibrium splitting')
            if (.not. mhd_energy) &
              call mpistop('LTE EoS requires mhd_energy=.true.')
            ! Route LTE conversion based on energy formulation
            if (mhd_internal_e) then
                eos%to_primitive  => mhd_to_primitive_inte_LTE
                eos%to_conserved  => mhd_to_conserved_inte_LTE
            else
                eos%to_primitive  => mhd_to_primitive_origin_LTE
                eos%to_conserved  => mhd_to_conserved_origin_LTE
            end if
        else
            call mpistop('Error: Unknown MHD EOS type: ' // trim(eos%eos_type))
        end if

        phys_to_primitive       => eos%to_primitive
        phys_to_conserved       => eos%to_conserved
        phys_bind_eos_to_source => bind_eos_to_source

        !> p_to_e (pressure -> total energy for origin, or eint for inte)
        if (mhd_internal_e) then
            eos%p_to_e => mhd_p_to_eint
        else
            eos%p_to_e => mhd_p_to_e
        end if

        !> Sound speed and Gamma1
        if (eos%eos_type == 'LTE' .and. eos%ionE) then
            eos%get_csound2 => mhd_get_csound2_LTE
            phys_get_gamma1 => mhd_get_gamma1_LTE
        else
            eos%get_csound2 => mhd_get_csound2_FI
            phys_get_gamma1 => get_gamma1_FI
        end if

        !> Prolongation
        if (eos%eos_type == 'LTE' .and. eos%ionE) then
            phys_to_prolong   => mhd_to_prolong_LTE
            phys_from_prolong => mhd_from_prolong_LTE
        end if

        !> Rfactor: only the FI case is physics-dependent (the usr_Rfactor user
        !> hook). The LTE and PI bindings are pure eos% routines, so they are set
        !> in eos_finalise_{LTE,PI} (which run after this and before
        !> bind_eos_to_source, so they win) -- keeping those targets private.
        if(associated(usr_Rfactor)) then
            eos%get_Rfactor => usr_Rfactor
        else
            eos%get_Rfactor => Rfactor_from_constant_ionization
        end if


        !> Thermal pressure routing
        if(mhd_internal_e) then
            phys_get_pthermal        => mhd_get_pthermal_inte
            mhd_get_pthermal         => mhd_get_pthermal_inte
        else if(mhd_hydrodynamic_e) then
            phys_get_pthermal        => mhd_get_pthermal_hde
            mhd_get_pthermal         => mhd_get_pthermal_hde
        else if(mhd_semirelativistic) then
            phys_get_pthermal        => mhd_get_pthermal_semirelati
            mhd_get_pthermal         => mhd_get_pthermal_semirelati
        else if(mhd_energy) then
            phys_get_pthermal        => mhd_get_pthermal_origin
            mhd_get_pthermal         => mhd_get_pthermal_origin
        else
            phys_get_pthermal        => mhd_get_pthermal_noe
            mhd_get_pthermal         => mhd_get_pthermal_noe
        end if

        ! For LTE, override with version that delegates to eos%get_thermal_pressure
        if (eos%eos_type == 'LTE') then
            phys_get_pthermal => mhd_get_pthermal_LTE
            mhd_get_pthermal  => mhd_get_pthermal_LTE
        end if
        ! Single entry point: all callers use eos%get_thermal_pressure
        eos%get_thermal_pressure => mhd_get_pthermal

        !> Temperature routing
        if(eos%eos_type == 'PI' .or. eos%eos_type == 'LTE') then
            mhd_get_temperature => mhd_get_temperature_from_Te
        else
            if(mhd_internal_e) then
                if(has_equi_rho_and_p) then
                    mhd_get_temperature => mhd_get_temperature_from_eint_with_equi
                else
                    mhd_get_temperature => mhd_get_temperature_from_eint
                end if
            else
                if(has_equi_rho_and_p) then
                    mhd_get_temperature => mhd_get_temperature_from_etot_with_equi
                else
                    mhd_get_temperature => mhd_get_temperature_from_etot
                end if
            end if
        end if

        !> Internal energy extraction
        if(mhd_internal_e) then
            phys_get_ei => mhd_get_ei_inte
        else if(mhd_energy) then
            phys_get_ei => mhd_get_ei_origin
        end if

        !> PI energy-mode overrides (eos_type='PI', ionE=.true.)
        ! PI shares FI's normalisation and no-energy path verbatim; energy mode
        ! differs ONLY in the eint<->p relation (eint carries the ionisation
        ! energy). Override exactly the routines that touch that relation, all
        ! delegating to the portable scalar backend in mod_eos_PI. FI
        ! routines stay pristine (no flag branch in the hot path).
        if (eos%eos_type == 'PI' .and. eos%ionE) then
            if (.not. mhd_energy) &
              call mpistop('PI energy EoS requires mhd_energy=.true.')
            if (mhd_hydrodynamic_e) &
              call mpistop('PI energy EoS not supported with mhd_hydrodynamic_e')
            if (mhd_semirelativistic) &
              call mpistop('PI energy EoS not supported with mhd_semirelativistic')
            if (has_equi_rho_and_p) &
              call mpistop('PI energy EoS not supported with equilibrium splitting')

            if (mhd_internal_e) then
                eos%to_conserved  => mhd_to_conserved_inte_PI
                eos%to_primitive  => mhd_to_primitive_inte_PI
                eos%p_to_e        => mhd_p_to_eint_PI
                phys_get_pthermal => mhd_get_pthermal_inte_PI
                mhd_get_pthermal  => mhd_get_pthermal_inte_PI
            else
                eos%to_conserved  => mhd_to_conserved_origin_PI
                eos%to_primitive  => mhd_to_primitive_origin_PI
                eos%p_to_e        => mhd_p_to_e_PI
                phys_get_pthermal => mhd_get_pthermal_origin_PI
                mhd_get_pthermal  => mhd_get_pthermal_origin_PI
            end if
            phys_to_primitive        => eos%to_primitive
            phys_to_conserved        => eos%to_conserved
            eos%get_thermal_pressure => mhd_get_pthermal

            !> eos%get_csound2 => get_csound2_PI is physics-independent: set in
            !> eos_finalise_PI (ionE arm) so the target stays private.
            phys_get_gamma1 => get_gamma1_PI

            ! Te_ refreshed from eint each substep (direct inversion, no lag).
            ! Module-specific (uses MHD's Te_ index): PI registers Te via
            ! var_set_auxvar, which does NOT set the generic iw_te, so the
            ! temperature aux must be addressed through the module's own Te_.
        end if

        ! phys_e_to_ei / phys_ei_to_e assigned in mhd_phys_init (mod_mhd_phys.t)
        mhd_to_primitive => eos%to_primitive
        mhd_to_conserved => eos%to_conserved

    end subroutine mhd_link_eos

    !> Called by eos_finalise via phys_bind_eos_to_source to link
    !> TC and RC modules to EoS-aware function pointers.
    subroutine bind_eos_to_source()

        ! Override eos%get_temperature_from_{etot,eint} per (eos_type, internal_e, equi).
        ! Runs AFTER eos_finalise's unconditional generic-helper assignment.
        ! Two semantically distinct hooks:
        !   - from_etot: caller's w(:,e_) is total energy; subtract KE+ME before T.
        !   - from_eint: caller's w(:,e_) is gas internal energy; compute T directly.
        ! In mhd_internal_e mode w(:,e_) IS already eint, so the etot hook is bound
        ! to the same routine as the eint hook (no subtraction). Required because
        ! eos_finalise's generic helper calls phys_e_to_ei, which is left unbound
        ! by mhd_phys_init when mhd_internal_e=.true.
        if (eos%eos_type == 'LTE') then
            if (mhd_internal_e) then
                ! w(:,e_) is eint; both hooks do the LTE table lookup directly.
                eos%get_temperature_from_etot => eos%get_temperature_from_eint
            else
                ! Total-energy: subtract KE+ME (via phys_e_to_ei) then LTE table lookup.
                eos%get_temperature_from_etot => mhd_get_temperature_from_etot_LTE
            end if
        else  ! FI
            if (mhd_internal_e) then
                ! w(:,e_) is eint; both hooks share the eint variant.
                if (has_equi_rho_and_p) then
                    eos%get_temperature_from_etot => mhd_get_temperature_from_eint_with_equi
                    eos%get_temperature_from_eint => mhd_get_temperature_from_eint_with_equi
                else
                    eos%get_temperature_from_etot => mhd_get_temperature_from_eint
                    eos%get_temperature_from_eint => mhd_get_temperature_from_eint
                end if
            else
                ! Total-energy: distinct etot/eint paths.
                if (has_equi_rho_and_p) then
                    eos%get_temperature_from_etot => mhd_get_temperature_from_etot_with_equi
                    eos%get_temperature_from_eint => mhd_get_temperature_from_eint_with_equi
                else
                    eos%get_temperature_from_etot => mhd_get_temperature_from_etot
                    eos%get_temperature_from_eint => mhd_get_temperature_from_eint
                end if
            end if
        end if

        if (allocated(tc_fl)) then
            tc_fl%get_temperature_from_conserved => eos%get_temperature_from_etot
            if (eos%eos_type == 'LTE' .and. eos%ionE) then
                tc_fl%get_temperature_from_eint => get_temperature_from_eint_fast_LTE
            else
                tc_fl%get_temperature_from_eint => eos%get_temperature_from_eint
            end if
            tc_fl%get_rho => eos%get_rho
            tc_fl%get_ne_nH => eos%get_ne_nH
            tc_fl%get_var_Rfactor => eos%get_Rfactor
            tc_fl%inv_gamma_minus_1 =  eos%inv_gamma_minus_1
            tc_fl%nH2rhoFactor      =  eos%nH2rhoFactor
            tc_fl%log_T_floor       =  eos_get_log_T_floor()
            tc_fl%eint_from_T       => eint_nH_from_T
            ! Equilibrium-specific pointers
            if(has_equi_rho_and_p .and. mhd_equi_thermal) then
                tc_fl%subtract_equi = .true.
                tc_fl%get_temperature_equi => mhd_get_temperature_equi
                tc_fl%get_rho_equi => mhd_get_rho_equi
            else
                tc_fl%subtract_equi = .false.
            end if
        end if

        if (allocated(rc_fl)) then
            rc_fl%get_rho          => eos%get_rho
            rc_fl%get_pthermal     => eos%get_thermal_pressure
            rc_fl%get_var_Rfactor  => eos%get_Rfactor
            rc_fl%get_Te           => eos%get_Te
            rc_fl%get_ne_nH        => eos%get_ne_nH
            nullify(rc_fl%get_rho2_factor)
            if(mhd_uawsom) rc_fl%get_rho2_factor => mhd_uawsom_rho2_factor
            rc_fl%ionE              =  eos%ionE
            rc_fl%method            =  eos%method
            rc_fl%inv_gamma_minus_1 =  eos%inv_gamma_minus_1
            rc_fl%nH2rhoFactor      =  eos%nH2rhoFactor
            rc_fl%eion_per_nH       =  eos%eion_per_nH
            rc_fl%eint_from_T       => eint_nH_from_T
            rc_fl%p2eint            => p2eint_from_nH_p
            rc_fl%T_from_eint       => T_from_nH_eint
            rc_fl%y_from_eint       => y_from_nH_eint
            ! Equilibrium-specific pointers
            if(has_equi_rho_and_p .and. mhd_equi_thermal) then
                rc_fl%subtract_equi = .true.
                rc_fl%get_rho_equi => mhd_get_rho_equi
                rc_fl%get_pthermal_equi => mhd_get_pe_equi
                rc_fl%get_ne_nH_equi => mhd_get_ne_nH_equi
                rc_fl%get_temperature_equi => mhd_get_temperature_equi
            else
                rc_fl%subtract_equi = .false.
            end if
            !> Build the variable-c_V Townsend Y_mod table now that all
            !> EoS tables (eint_from_T, T, neOnH) are in code units.
            !> build_Y_mod_table checks coolmethod=='exact' and .not.isPPL
            !> internally and early-returns otherwise. Feed it the inverse-table nH
            !> grid via the port so it never touches eos% directly. LTE only:
            !> eos_get_eintT_grid reads the LTE eintT grid, and PI (T-only) cooling
            !> runs the classical Townsend path (Y_mod-for-PI is a later refinement).
            if (eos%ionE .and. eos%eos_type == 'LTE') then
                call eos_get_eintT_grid(rc_fl%Y_mod_n_nH, &
                     rc_fl%Y_mod_lg_nH_min, rc_fl%Y_mod_lg_nH_max)
                call build_Y_mod_table(rc_fl)
            end if
        end if

        !> PI energy mode: the eint<->T scalar callbacks wired above are the
        !> LTE-table lookups, which PI does not load. Repoint cooling/conduction
        !> at the ionisation backend (rho=nH=1 per H). Classical Townsend; Y_mod
        !> is not built for PI (guard above), so the classical branch is used.
        if (eos%eos_type == 'PI' .and. eos%ionE) then
            if (allocated(tc_fl)) tc_fl%eint_from_T => eint_from_T_PI
            if (allocated(rc_fl)) then
                rc_fl%eint_from_T => eint_from_T_PI
                rc_fl%p2eint      => p2eint_PI
                rc_fl%T_from_eint => T_from_eint_PI
                rc_fl%y_from_eint => y_from_eint_PI
            end if
        end if

        if (allocated(te_fl_mhd)) then
            te_fl_mhd%get_rho          => eos%get_rho
            te_fl_mhd%get_pthermal     => eos%get_thermal_pressure
            te_fl_mhd%get_var_Rfactor  => eos%get_Rfactor
            te_fl_mhd%get_ne_nH        => eos%get_ne_nH
        end if

        if (allocated(fld_fl)) then
            !> Radiation (FLD) fluid: gas-EoS callbacks. get_temperature_from_pressure
            !> is the FI T=p/(R*rho) routine (the LTE variant is a later pass).
            fld_fl%gamma       =  eos%gamma
            fld_fl%get_tgas    => eos%get_temperature_from_pressure
            fld_fl%get_Rfactor => eos%get_Rfactor
        end if

    end subroutine bind_eos_to_source

    !> FI conversion routines (signatures updated for convert_condition interface)

    !> Transform primitive variables into conservative ones (origin, total energy)
    subroutine mhd_to_conserved_origin(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            ! Calculate total energy from pressure, kinetic and magnetic energy
            w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1&
                      +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                      +(^C&w(ix^D,b^C_)**2+))
            if(mhd_uawsom) w(ix^D,e_)=w(ix^D,e_)+&
                 w(ix^D,wAplus_)+w(ix^D,wAminus_)+w(ix^D,wkplus_)+w(ix^D,wkminus_)
            ! Convert velocity to momentum
            ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
        {end do\}

    end subroutine mhd_to_conserved_origin

    !> Transform primitive variables into conservative ones (no energy)
    subroutine mhd_to_conserved_origin_noe(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            ! Convert velocity to momentum
            ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
        {end do\}

    end subroutine mhd_to_conserved_origin_noe

    !> Transform primitive variables into conservative ones (hydrodynamic energy)
    subroutine mhd_to_conserved_hde(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            ! Calculate total energy from pressure, kinetic and magnetic energy
            w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1&
                      +half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
            ! Convert velocity to momentum
            ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
        {end do\}

    end subroutine mhd_to_conserved_hde

    !> Transform primitive variables into conservative ones (internal energy)
    subroutine mhd_to_conserved_inte(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            ! Calculate internal energy from pressure
            w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1
            ! Convert velocity to momentum
            ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
        {end do\}

    end subroutine mhd_to_conserved_inte

    !> Transform primitive variables into conservative ones (split rho)
    subroutine mhd_to_conserved_split_rho(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: rho
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            rho=w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i)
            ! Calculate total energy from pressure, kinetic and magnetic energy
            w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1&
                      +half*((^C&w(ix^D,m^C_)**2+)*rho&
                            +(^C&w(ix^D,b^C_)**2+))
            ! Convert velocity to momentum
            ^C&w(ix^D,m^C_)=rho*w(ix^D,m^C_)\
        {end do\}

    end subroutine mhd_to_conserved_split_rho

    !> Transform primitive variables into conservative ones (semirelativistic)
    subroutine mhd_to_conserved_semirelati(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: E(ixO^S,1:ndir), S(ixO^S,1:ndir)
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            {^IFTHREEC
                E(ix^D,1)=w(ix^D,b2_)*w(ix^D,m3_)-w(ix^D,b3_)*w(ix^D,m2_)
                E(ix^D,2)=w(ix^D,b3_)*w(ix^D,m1_)-w(ix^D,b1_)*w(ix^D,m3_)
                E(ix^D,3)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
                S(ix^D,1)=E(ix^D,2)*w(ix^D,b3_)-E(ix^D,3)*w(ix^D,b2_)
                S(ix^D,2)=E(ix^D,3)*w(ix^D,b1_)-E(ix^D,1)*w(ix^D,b3_)
                S(ix^D,3)=E(ix^D,1)*w(ix^D,b2_)-E(ix^D,2)*w(ix^D,b1_)
            }
            {^IFTWOC
                E(ix^D,1)=zero
                ! switch 3 with 2 to add 3 when ^C from 1 to 2
                E(ix^D,2)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
                S(ix^D,1)=-E(ix^D,2)*w(ix^D,b2_)
                S(ix^D,2)=E(ix^D,2)*w(ix^D,b1_)
            }
            {^IFONEC
                E(ix^D,1)=zero
                S(ix^D,1)=zero
            }
            if(mhd_internal_e) then
                ! internal energy
                w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1
            else
                ! equation (9)
                ! Calculate total energy from internal, kinetic and magnetic energy
                w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1&
                          +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                          +(^C&w(ix^D,b^C_)**2+)&
                          +(^C&e(ix^D,^C)**2+)*eos%inv_squared_c)
            end if

            ! Convert velocity to momentum, equation (9)
            ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)+S(ix^D,^C)*eos%inv_squared_c\

        {end do\}

    end subroutine mhd_to_conserved_semirelati

    subroutine mhd_to_conserved_semirelati_noe(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: E(ixO^S,1:ndir), S(ixO^S,1:ndir)
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            {^IFTHREEC
                E(ix^D,1)=w(ix^D,b2_)*w(ix^D,m3_)-w(ix^D,b3_)*w(ix^D,m2_)
                E(ix^D,2)=w(ix^D,b3_)*w(ix^D,m1_)-w(ix^D,b1_)*w(ix^D,m3_)
                E(ix^D,3)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
                S(ix^D,1)=E(ix^D,2)*w(ix^D,b3_)-E(ix^D,3)*w(ix^D,b2_)
                S(ix^D,2)=E(ix^D,3)*w(ix^D,b1_)-E(ix^D,1)*w(ix^D,b3_)
                S(ix^D,3)=E(ix^D,1)*w(ix^D,b2_)-E(ix^D,2)*w(ix^D,b1_)
            }
            {^IFTWOC
                E(ix^D,1)=zero
                ! switch 3 with 2 to add 3 when ^C from 1 to 2
                E(ix^D,2)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
                S(ix^D,1)=-E(ix^D,2)*w(ix^D,b2_)
                S(ix^D,2)=E(ix^D,2)*w(ix^D,b1_)
            }
            {^IFONEC
                S(ix^D,1)=zero
            }
            ! Convert velocity to momentum, equation (9)
            ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)+S(ix^D,^C)*eos%inv_squared_c\

        {end do\}

    end subroutine mhd_to_conserved_semirelati_noe

    !> Transform conservative variables into primitive ones (origin)
    subroutine mhd_to_primitive_origin(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision                :: inv_rho
        integer :: ix^D

        if (fix_small_values) then
            call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_origin')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0/w(ix^D,rho_)
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
            ! Calculate pressure = (gamma-1) * (e-ek-eb)
            w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)&
                      -half*(w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+)&
                        +(^C&w(ix^D,b^C_)**2+))&
                      -mhd_uawsom_wave_energy_cell(w(ix^D,:)))
        {end do\}

    end subroutine mhd_to_primitive_origin

    !> Transform conservative variables into primitive ones (no energy)
    subroutine mhd_to_primitive_origin_noe(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision                :: inv_rho
        integer :: ix^D

        if (fix_small_values) then
            call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_origin_noe')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0/w(ix^D,rho_)
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
        {end do\}

    end subroutine mhd_to_primitive_origin_noe

    !> Transform conservative variables into primitive ones (hde)
    subroutine mhd_to_primitive_hde(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision                :: inv_rho
        integer                         :: ix^D

        if (fix_small_values) then
            call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_hde')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0/w(ix^D,rho_)
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
            ! Calculate pressure = (gamma-1) * (e-ek)
            w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)-half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+))
        {end do\}

    end subroutine mhd_to_primitive_hde

    !> Transform conservative variables into primitive ones (internal energy)
    subroutine mhd_to_primitive_inte(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision                :: inv_rho
        integer                         :: ix^D

        if (fix_small_values) then
            call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_inte')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            ! Calculate pressure = (gamma-1) * e_internal
            w(ix^D,p_)=w(ix^D,e_)*eos%gamma_minus_1
            ! Convert momentum to velocity
            inv_rho = 1.d0/w(ix^D,rho_)
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
        {end do\}

    end subroutine mhd_to_primitive_inte

    !> Transform conservative variables into primitive ones (split rho)
    subroutine mhd_to_primitive_split_rho(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: inv_rho
        integer :: ix^D

        if (fix_small_values) then
            call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_split_rho')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho=1.d0/(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
            ! Calculate pressure = (gamma-1) * (e-ek-eb)
            w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)&
                        -half*((w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,b0i))*&
                        (^C&w(ix^D,m^C_)**2+)+(^C&w(ix^D,b^C_)**2+)))
        {end do\}

    end subroutine mhd_to_primitive_split_rho

    !> Transform conservative variables into primitive ones (semirelativistic)
    subroutine mhd_to_primitive_semirelati(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: b(ixO^S,1:ndir), tmp, b2, gamma2, inv_rho
        integer :: ix^D

        if (fix_small_values) then
            call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_semirelati')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            b2=(^C&w(ix^D,b^C_)**2+)
            if(b2>smalldouble) then
                tmp=1.d0/sqrt(b2)
            else
                tmp=0.d0
            end if
            ^C&b(ix^D,^C)=w(ix^D,b^C_)*tmp\
            tmp=(^C&b(ix^D,^C)*w(ix^D,m^C_)+)

            inv_rho=1.d0/w(ix^D,rho_)
            ! Va^2/c^2
            b2=b2*inv_rho*eos%inv_squared_c
            ! equation (15)
            gamma2=1.d0/(1.d0+b2)
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=gamma2*(w(ix^D,m^C_)+b2*b(ix^D,^C)*tmp)*inv_rho\

            if(mhd_internal_e) then
                ! internal energy to pressure
                w(ix^D,p_)=eos%gamma_minus_1*w(ix^D,e_)
            else
                ! E=Bxv
                {^IFTHREEC
                    b(ix^D,1)=w(ix^D,b2_)*w(ix^D,m3_)-w(ix^D,b3_)*w(ix^D,m2_)
                    b(ix^D,2)=w(ix^D,b3_)*w(ix^D,m1_)-w(ix^D,b1_)*w(ix^D,m3_)
                    b(ix^D,3)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
                }
                {^IFTWOC
                    b(ix^D,1)=zero
                    b(ix^D,2)=w(ix^D,b1_)*w(ix^D,m2_)-w(ix^D,b2_)*w(ix^D,m1_)
                }
                {^IFONEC
                    b(ix^D,1)=zero
                }
                ! Calculate pressure = (gamma-1) * (e-eK-eB-eE)
                w(ix^D,p_)=eos%gamma_minus_1*(w(ix^D,e_)&
                          -half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                          +(^C&w(ix^D,b^C_)**2+)&
                          +(^C&b(ix^D,^C)**2+)*eos%inv_squared_c))
            end if
        {end do\}

    end subroutine mhd_to_primitive_semirelati

    !> Transform conservative variables into primitive ones (semirelativistic noe)
    subroutine mhd_to_primitive_semirelati_noe(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: b(ixO^S,1:ndir),tmp,b2,gamma2,inv_rho
        integer :: ix^D, idir

        if (fix_small_values) then
            call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, 'mhd_to_primitive_semirelati_noe')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            b2=(^C&w(ix^D,b^C_)**2+)
            if(b2>smalldouble) then
                tmp=1.d0/sqrt(b2)
            else
                tmp=0.d0
            end if
            ^C&b(ix^D,^C)=w(ix^D,b^C_)*tmp\
            tmp=(^C&b(ix^D,^C)*w(ix^D,m^C_)+)

            inv_rho=1.d0/w(ix^D,rho_)
            ! Va^2/c^2
            b2=b2*inv_rho*eos%inv_squared_c
            ! equation (15)
            gamma2=1.d0/(1.d0+b2)
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=gamma2*(w(ix^D,m^C_)+b2*b(ix^D,^C)*tmp)*inv_rho\
        {end do\}

    end subroutine mhd_to_primitive_semirelati_noe

    !> LTE conversion routines

    !> LTE: Primitive (rho, v, p, B) -> Conserved (rho, rho*v, E_total, B)
    !> E_total = eint + 0.5*rho*v^2 + 0.5*B^2
    subroutine mhd_to_conserved_origin_LTE(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        timeeos0 = MPI_WTIME()

        ! Convert p -> E_total (eint + eK + eB) via EoS table
        call mhd_p_to_e(ixI^L, ixO^L, w, x)
        ! Convert velocity to momentum
        ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,m^C_)\

        timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine mhd_to_conserved_origin_LTE

    !> LTE: Primitive (rho, v, p, B) -> Conserved (rho, rho*v, eint, B)
    !> Internal energy formulation: e_ stores eint only
    subroutine mhd_to_conserved_inte_LTE(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        timeeos0 = MPI_WTIME()

        ! Convert p -> eint via EoS table (no mechanical energy)
        call mhd_p_to_eint(ixI^L, ixO^L, w, x)
        ! Convert velocity to momentum
        ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,m^C_)\

        timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine mhd_to_conserved_inte_LTE

    !> Convert pressure to total energy (eint + eK + eB) for origin formulation.
    !> Uses p2eint table for LTE ionE, with FI bypass for hot cells.
    subroutine mhd_p_to_e(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        integer :: ix^D
        double precision :: p_to_eint, p_over_rho
        double precision :: nH(ixI^S), nH_in(ixI^S), p_in(ixI^S)
        double precision :: T_solve, y_solve, eint_nH_solve
        double precision :: log_eint_mid, eint_total

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
                    w(ix^D,e_)=w(ix^D,p_)*p_to_eint&
                              +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                              +(^C&w(ix^D,b^C_)**2+))
                else if (eos%method == 'analytic') then
                    call saha_state_from_nH_p(nH(ix^D), w(ix^D,p_), &
                        T_solve, y_solve, eint_nH_solve)
                    w(ix^D,e_) = eint_nH_solve * nH(ix^D) &
                        +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                        +(^C&w(ix^D,b^C_)**2+))
                else if (eos%p2eint_method == 'bisect') then
                    call eint_from_p_bisect(nH_in(ix^D), &
                        dlog10(w(ix^D,p_)), log_eint_mid)
                    eint_total = nH(ix^D) * 10.0d0**log_eint_mid
                    eint_total = max(eint_total, &
                        nH(ix^D) * 10.0d0**eos%T%var2_min)
                    w(ix^D,e_) = eint_total + &
                        half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                        +(^C&w(ix^D,b^C_)**2+))
                else
                    p_to_eint = p2eint_from_nH_p(nH_in(ix^D), p_in(ix^D))
                    w(ix^D,e_)=w(ix^D,p_)*p_to_eint&
                              +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                              +(^C&w(ix^D,b^C_)**2+))
                end if
            else
                w(ix^D,e_)=w(ix^D,p_)*p_to_eint&
                          +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                          +(^C&w(ix^D,b^C_)**2+))
            end if
            if(mhd_uawsom) w(ix^D,e_)=w(ix^D,e_)+&
                 w(ix^D,wAplus_)+w(ix^D,wAminus_)+w(ix^D,wkplus_)+w(ix^D,wkminus_)
        {end do\}

    end subroutine mhd_p_to_e

    !> Convert pressure to internal energy only (for internal_e formulation).
    !> Uses p2eint table for LTE ionE, with FI bypass for hot cells.
    subroutine mhd_p_to_eint(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        integer :: ix^D
        double precision :: p_to_eint, p_over_rho
        double precision :: nH(ixI^S), nH_in(ixI^S), p_in(ixI^S)
        double precision :: T_solve, y_solve, eint_nH_solve
        double precision :: log_eint_mid, eint_total

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
                    w(ix^D,e_) = w(ix^D,p_) * p_to_eint
                else if (eos%method == 'analytic') then
                    call saha_state_from_nH_p(nH(ix^D), w(ix^D,p_), &
                        T_solve, y_solve, eint_nH_solve)
                    w(ix^D,e_) = eint_nH_solve * nH(ix^D)
                else if (eos%p2eint_method == 'bisect') then
                    call eint_from_p_bisect(nH_in(ix^D), &
                        dlog10(w(ix^D,p_)), log_eint_mid)
                    eint_total = nH(ix^D) * 10.0d0**log_eint_mid
                    eint_total = max(eint_total, &
                        nH(ix^D) * 10.0d0**eos%T%var2_min)
                    w(ix^D,e_) = eint_total
                else
                    p_to_eint = p2eint_from_nH_p(nH_in(ix^D), p_in(ix^D))
                    w(ix^D,e_) = w(ix^D,p_) * p_to_eint
                end if
            else
                w(ix^D,e_) = w(ix^D,p_) * p_to_eint
            end if
        {end do\}

    end subroutine mhd_p_to_eint

    !> LTE: Conserved (rho, rho*v, E_total, B) -> Primitive (rho, v, p, B)
    !> Uses T and neOnH tables with FI bypass for hot cells.
    subroutine mhd_to_primitive_origin_LTE(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: inv_rho, eint_val, eint_in
        double precision :: nH(ixI^S), log_nH(ixI^S)
        integer :: ix^D

        timeeos0 = MPI_WTIME()

        if (fix_small_values) then
            call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, &
                'mhd_to_primitive_origin_LTE')
        end if

        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        if (eos%ionE) then
            log_nH(ixO^S) = dlog10(nH(ixO^S))
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0/w(ix^D,rho_)
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
            ! Extract internal energy: eint = E - 0.5*rho*v^2 - 0.5*B^2
            eint_val = w(ix^D,e_) &
                      - half*(w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+) &
                      + (^C&w(ix^D,b^C_)**2+))
            ! Floor eint to prevent unphysical values
            if (eos%method /= 'analytic') then
                eint_val = max(eint_val, nH(ix^D) * 10.0d0**eos%T%var2_min)
            end if
            eint_val = max(eint_val, smalldouble)

            if (eos%ionE) then
                if (eint_val * inv_rho > eos%eint_rho_FI_threshold) then
                    ! FI bypass: p = (gamma-1) * (eint - eion*nH)
                    w(ix^D,p_) = eos%gamma_minus_1 &
                        * (eint_val - eos%eion_per_nH * nH(ix^D))
                else
                    ! Ionisation zone: single p/nH lookup (replaces separate T + y)
                    eint_in = dlog10(eint_val) - log_nH(ix^D)
                    w(ix^D,p_) = nH(ix^D) * p_nH_from_eint(log_nH(ix^D), eint_in)
                end if
            else
                ! FI without ionE tables
                w(ix^D,p_) = eos%gamma_minus_1 * eint_val
            end if
        {end do\}

        timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine mhd_to_primitive_origin_LTE

    !> LTE: Conserved (rho, rho*v, eint, B) -> Primitive (rho, v, p, B)
    !> Internal energy formulation: eint = w(e_) directly.
    subroutine mhd_to_primitive_inte_LTE(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: inv_rho, eint_val, eint_in, T_loc, y_loc
        double precision :: nH(ixI^S), log_nH(ixI^S)
        integer :: ix^D

        timeeos0 = MPI_WTIME()

        if (fix_small_values) then
            call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, &
                'mhd_to_primitive_inte_LTE')
        end if

        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        if (eos%ionE) then
            log_nH(ixO^S) = dlog10(nH(ixO^S))
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            ! eint = w(e_) directly for internal energy formulation
            eint_val = w(ix^D,e_)
            eint_val = max(eint_val, nH(ix^D) * 10.0d0**eos%T%var2_min)

            if (eos%ionE) then
                if (eint_val / w(ix^D,rho_) > eos%eint_rho_FI_threshold) then
                    w(ix^D,p_) = eos%gamma_minus_1 &
                        * (eint_val - eos%eion_per_nH * nH(ix^D))
                else
                    eint_in = dlog10(eint_val) - log_nH(ix^D)
                    w(ix^D,p_) = nH(ix^D) * p_nH_from_eint(log_nH(ix^D), eint_in)
                end if
            else
                w(ix^D,p_) = eos%gamma_minus_1 * eint_val
            end if

            ! Convert momentum to velocity
            inv_rho = 1.d0/w(ix^D,rho_)
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
        {end do\}

        timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine mhd_to_primitive_inte_LTE

    !> LTE prolongation routines

    !> Convert conserved (rho, rho*v, E_total, B) to prolongation form (rho, v, T, B).
    !> Temperature is stored in the p_ slot for AMR interpolation.
    !> T is smoothest through the ionisation zone (conduction-dominated, linear gradient).
    subroutine mhd_to_prolong_LTE(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: inv_rho, eint_val, T_loc, y_loc
        double precision :: nH(ixI^S), log_nH(ixI^S)
        integer :: ix^D

        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        log_nH(ixO^S) = dlog10(nH(ixO^S))

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0/w(ix^D,rho_)
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
            if (mhd_internal_e) then
                eint_val = w(ix^D,e_)
            else
                ! Compute eint = E - 0.5*rho*v^2 - 0.5*B^2
                eint_val = w(ix^D,e_) &
                    - half*(w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+) &
                    + (^C&w(ix^D,b^C_)**2+))
            end if
            ! Floor eint to prevent unphysical values
            if (eos%method /= 'analytic') then
                eint_val = max(eint_val, nH(ix^D) * 10.0d0**eos%T%var2_min)
            end if
            eint_val = max(eint_val, smalldouble)

            if (eint_val * inv_rho > eos%eint_rho_FI_threshold) then
                ! FI: T = (gamma-1)*(eint - eion*nH) / (nH * n_per_nH_FI)
                w(ix^D,p_) = eos%gamma_minus_1 &
                    * (eint_val - eos%eion_per_nH * nH(ix^D)) &
                    / (nH(ix^D) * eos%n_per_nH_FI)
            else if (eos%method == 'analytic') then
                ! Analytical Saha: solve for T from eint
                call saha_T_from_nH_eint(nH(ix^D), &
                    eint_val / nH(ix^D), T_loc, y_loc)
                w(ix^D,p_) = T_loc
            else
                ! Ionisation zone: T from table
                w(ix^D,p_) = T_from_nH_eint( &
                    log_nH(ix^D), &
                    dlog10(eint_val) - log_nH(ix^D))
            end if
        {end do\}

    end subroutine mhd_to_prolong_LTE

    !> Convert prolongation form (rho, v, T, B) to conserved (rho, rho*v, E_total, B).
    !> T is read from the p_ slot. Uses eint_nH_from_T table for back-conversion.
    subroutine mhd_from_prolong_LTE(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: T_val, eint_val, T_FI, log_T_min
        double precision :: nH(ixI^S), log_nH(ixI^S)
        integer :: ix^D

        ! Temperature above which gas is fully ionised
        T_FI = (eos%eint_rho_FI_threshold &
            * eos%nH2rhoFactor - eos%eion_per_nH) &
            * eos%gamma_minus_1 / eos%n_per_nH_FI

        ! Floor for log_T into the (rho, T) inverse table; entropy method uses
        ! eos%eintT, legacy 'tables' method uses eos%eint_from_T. Picking the
        ! wrong container leaves var2_min = 0 → floors T at 10^6 K.
        if (eos%method == 'entropy') then
            log_T_min = eos%eintT%var2_min
        else
            log_T_min = eos%eint_from_T%var2_min
        end if

        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        log_nH(ixO^S) = dlog10(nH(ixO^S))

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            T_val = w(ix^D,p_)  ! T stored in p_ slot
            if (T_val > T_FI) then
                ! FI: eint = nH*(n_per_nH*T/(gamma-1) + eion)
                eint_val = nH(ix^D) &
                    * (eos%n_per_nH_FI * T_val * eos%inv_gamma_minus_1 &
                    + eos%eion_per_nH)
            else if (eos%method == 'analytic') then
                ! Analytical Saha: eint from T directly
                eint_val = saha_eint_from_nH_T(nH(ix^D), T_val) * nH(ix^D)
            else
                ! Ionisation zone: eint/nH from T table
                eint_val = eint_nH_from_T( &
                    log_nH(ix^D), &
                    dlog10(max(T_val, 10.0d0**log_T_min))) &
                    * nH(ix^D)
            end if
            if (mhd_internal_e) then
                w(ix^D,e_) = eint_val
            else
                ! E = eint + 0.5*rho*v^2 + 0.5*B^2
                w(ix^D,e_) = eint_val &
                    + half*(w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+) &
                    + (^C&w(ix^D,b^C_)**2+))
            end if
            ! Convert velocity to momentum
            ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
        {end do\}

    end subroutine mhd_from_prolong_LTE

    !> Sound speed routines

    !> Acoustic sound speed squared for FI (constant gamma) EoS.
    !> Expects w in primitive form: w(p_) = pressure, w(rho_) = density.
    subroutine mhd_get_csound2_FI(w, x, ixI^L, ixO^L, cs2)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision, intent(out)   :: cs2(ixI^S)

        double precision :: rho(ixI^S), pth(ixI^S)

        timeeos0 = MPI_WTIME()

        !> Both slots hold perturbations under equilibrium splitting, so the
        !> background is restored before the ratio is formed.
        call eos%get_rho(w, x, ixI^L, ixO^L, rho)
        pth(ixO^S) = w(ixO^S, p_)
        if (iw_equi_p > 0) pth(ixO^S) = pth(ixO^S) &
             + block%equi_vars(ixO^S, iw_equi_p, b0i)
        cs2(ixO^S) = eos%gamma * pth(ixO^S) / rho(ixO^S)

        timeeos_csound = timeeos_csound + (MPI_WTIME()-timeeos0)

    end subroutine mhd_get_csound2_FI

    !> Acoustic sound speed squared for LTE+IonE EoS using Gamma_1 table.
    !> Expects w in primitive form. Uses FI bypass for hot cells.
    subroutine mhd_get_csound2_LTE(w, x, ixI^L, ixO^L, cs2)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision, intent(out)   :: cs2(ixI^S)

        double precision :: nH_val, log_nH, log_p_nH, g1, p_over_rho
        integer :: ix^D

        timeeos0 = MPI_WTIME()

        if (eos%gamma1_method == 'constant') then
            cs2(ixO^S) = eos%gamma * w(ixO^S, p_) / w(ixO^S, rho_)
        else
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                p_over_rho = w(ix^D, p_) / w(ix^D, rho_)
                if (p_over_rho > eos%p_rho_FI_threshold) then
                    cs2(ix^D) = eos%gamma * p_over_rho
                else
                    nH_val = w(ix^D, rho_) / eos%nH2rhoFactor
                    if (eos%method == 'analytic') then
                        if (iw_te > 0 .and. w(ix^D,iw_te) > 0.0d0) then
                            g1 = saha_gamma1_from_nH_T(nH_val, w(ix^D,iw_te))
                        else
                            g1 = eos%gamma
                        end if
                    else
                        log_nH = dlog10(nH_val)
                        log_p_nH = dlog10(w(ix^D, p_) / nH_val)
                        g1 = gamma1_from_nH_p(log_nH, log_p_nH)
                    end if
                    cs2(ix^D) = g1 * p_over_rho
                end if
            {end do\}
        end if

        timeeos_csound = timeeos_csound + (MPI_WTIME()-timeeos0)

    end subroutine mhd_get_csound2_LTE

    subroutine mhd_get_gamma1_LTE(w, x, ixI^L, ixO^L, gamma1)
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

    end subroutine mhd_get_gamma1_LTE

    !> Rfactor routines (matching HD: EoS functions live in EoS module)

    !> Rfactor = p/(rho*T) for constant ionisation degree (FI/PI no-energy).
    !> Stays in the seam: RR is the physics module's gas-constant factor, not
    !> visible to mod_eos. Rfactor_from_LTE lives in mod_eos (no RR).
    subroutine Rfactor_from_constant_ionization(w,x,ixI^L,ixO^L,Rfactor)
        use mod_global_parameters
        integer, intent(in) :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: Rfactor(ixI^S)

        Rfactor(ixO^S)=RR

    end subroutine Rfactor_from_constant_ionization

    !> Rfactor from partial ionisation temperature lookup (Leenaarts et al. 2012)

    !> Thermal pressure routines (moved from mod_mhd_phys.t)

    !> Calculate isothermal thermal pressure
    subroutine mhd_get_pthermal_noe(w,x,ixI^L,ixO^L,pth)
        use mod_global_parameters

        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: pth(ixI^S)

        if(has_equi_rho_and_p) then
            pth(ixO^S)=mhd_adiab*(w(ixO^S,rho_)+block%equi_vars(ixO^S,equi_rho0_,0))**eos%gamma
        else
            pth(ixO^S)=mhd_adiab*w(ixO^S,rho_)**eos%gamma
        end if

    end subroutine mhd_get_pthermal_noe

    !> Calculate thermal pressure from internal energy
    subroutine mhd_get_pthermal_inte(w,x,ixI^L,ixO^L,pth)
        use mod_global_parameters
        use mod_small_values, only: trace_small_values

        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: pth(ixI^S)

        integer :: iw, ix^D

        {do ix^DB= ixOmin^DB,ixOmax^DB\}
            if(has_equi_rho_and_p) then
                pth(ix^D)=eos%gamma_minus_1*w(ix^D,e_)+block%equi_vars(ix^D,equi_pe0_,0)
            else
                pth(ix^D)=eos%gamma_minus_1*w(ix^D,e_)
            end if
            if(fix_small_values.and.pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
        {end do\}

        if(check_small_values.and..not.fix_small_values) then
            {do ix^DB= ixOmin^DB,ixOmax^DB\}
                if(pth(ix^D)<small_pressure) then
                    write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                         " encountered when call mhd_get_pthermal_inte"
                    write(*,*) "Iteration: ", it, " Time: ", global_time
                    write(*,*) "Location: ", x(ix^D,:)
                    write(*,*) "Cell number: ", ix^D
                    do iw=1,nw
                        write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
                    end do
                    if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
                    write(*,*) "Saving status at the previous time step"
                    crash=.true.
                end if
            {end do\}
        end if

    end subroutine mhd_get_pthermal_inte

    !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
    subroutine mhd_get_pthermal_origin(w,x,ixI^L,ixO^L,pth)
        use mod_global_parameters
        use mod_small_values, only: trace_small_values

        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: pth(ixI^S)

        integer :: iw, ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if(has_equi_rho_and_p) then
                pth(ix^D)=eos%gamma_minus_1*(w(ix^D,e_)-half*((^C&w(ix^D,m^C_)**2+)/(w(ix^D,rho_)+block%equi_vars(ix^D,equi_rho0_,0))&
                     +(^C&w(ix^D,b^C_)**2+)))+block%equi_vars(ix^D,equi_pe0_,0)
            else
                pth(ix^D)=eos%gamma_minus_1*(w(ix^D,e_)-half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)&
                     +(^C&w(ix^D,b^C_)**2+))-mhd_uawsom_wave_energy_cell(w(ix^D,:)))
            end if
            if(fix_small_values.and.pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
        {end do\}

        if(check_small_values.and..not.fix_small_values) then
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                if(pth(ix^D)<small_pressure) then
                    write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                         " encountered when call mhd_get_pthermal"
                    write(*,*) "Iteration: ", it, " Time: ", global_time
                    write(*,*) "Location: ", x(ix^D,:)
                    write(*,*) "Cell number: ", ix^D
                    do iw=1,nw
                        write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
                    end do
                    if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
                    write(*,*) "Saving status at the previous time step"
                    crash=.true.
                end if
            {end do\}
        end if

    end subroutine mhd_get_pthermal_origin

    !> Calculate thermal pressure for LTE EoS (delegates to eos%get_thermal_pressure)
    subroutine mhd_get_pthermal_LTE(w,x,ixI^L,ixO^L,pth)
        use mod_global_parameters
        use mod_small_values, only: trace_small_values

        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: pth(ixI^S)
        double precision :: nH(ixI^S)

        integer :: iw, ix^D

        ! LTE: p = nH * (1 + He + ne/nH) * T from stored state variables
        call eos%get_nH(w, x, ixI^L, ixO^L, nH)
        pth(ixO^S) = nH(ixO^S) * (1.0d0 + eos%He_abundance &
            + (w(ixO^S,Ne_) / nH(ixO^S))) * w(ixO^S,Te_)

        if(fix_small_values) then
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                if(pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
            {end do\}
        else if(check_small_values) then
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                if(pth(ix^D)<small_pressure) then
                    write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                         " encountered when call mhd_get_pthermal_LTE"
                    write(*,*) "Iteration: ", it, " Time: ", global_time
                    write(*,*) "Location: ", x(ix^D,:)
                    write(*,*) "Cell number: ", ix^D
                    do iw=1,nw
                        write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
                    end do
                    if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
                    write(*,*) "Saving status at the previous time step"
                    crash=.true.
                end if
            {end do\}
        end if

    end subroutine mhd_get_pthermal_LTE

    !> Calculate thermal pressure for semirelativistic MHD
    subroutine mhd_get_pthermal_semirelati(w,x,ixI^L,ixO^L,pth)
        use mod_global_parameters
        use mod_small_values, only: trace_small_values

        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: pth(ixI^S)

        double precision :: b(ixO^S,1:ndir), v(ixO^S,1:ndir), tmp, b2, gamma2, inv_rho
        integer :: iw, ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            b2=(^C&w(ix^D,b^C_)**2+)
            if(b2>smalldouble) then
                tmp=1.d0/sqrt(b2)
            else
                tmp=0.d0
            end if
            ^C&b(ix^D,^C)=w(ix^D,b^C_)*tmp\
            tmp=(^C&b(ix^D,^C)*w(ix^D,m^C_)+)

            inv_rho=1.d0/w(ix^D,rho_)
            ! Va^2/c^2
            b2=b2*inv_rho*eos%inv_squared_c
            ! equation (15)
            gamma2=1.d0/(1.d0+b2)
            ! Convert momentum to velocity
            ^C&v(ix^D,^C)=gamma2*(w(ix^D,m^C_)+b2*b(ix^D,^C)*tmp)*inv_rho\

            ! E=Bxv
            {^IFTHREEC
                b(ix^D,1)=w(ix^D,b2_)*v(ix^D,3)-w(ix^D,b3_)*v(ix^D,2)
                b(ix^D,2)=w(ix^D,b3_)*v(ix^D,1)-w(ix^D,b1_)*v(ix^D,3)
                b(ix^D,3)=w(ix^D,b1_)*v(ix^D,2)-w(ix^D,b2_)*v(ix^D,1)
            }
            {^IFTWOC
                b(ix^D,1)=zero
                b(ix^D,2)=w(ix^D,b1_)*v(ix^D,2)-w(ix^D,b2_)*v(ix^D,1)
            }
            {^IFONEC
                b(ix^D,1)=zero
            }
            ! Calculate pressure = (gamma-1) * (e-eK-eB-eE)
            pth(ix^D)=eos%gamma_minus_1*(w(ix^D,e_)&
                       -half*((^C&v(ix^D,^C)**2+)*w(ix^D,rho_)&
                       +(^C&w(ix^D,b^C_)**2+)&
                       +(^C&b(ix^D,^C)**2+)*eos%inv_squared_c))
            if(fix_small_values.and.pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
        {end do\}

        if(check_small_values.and..not.fix_small_values) then
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                if(pth(ix^D)<small_pressure) then
                    write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                         " encountered when call mhd_get_pthermal_semirelati"
                    write(*,*) "Iteration: ", it, " Time: ", global_time
                    write(*,*) "Location: ", x(ix^D,:)
                    write(*,*) "Cell number: ", ix^D
                    do iw=1,nw
                        write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
                    end do
                    if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
                    write(*,*) "Saving status at the previous time step"
                    crash=.true.
                end if
            {end do\}
        end if

    end subroutine mhd_get_pthermal_semirelati

    !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L
    subroutine mhd_get_pthermal_hde(w,x,ixI^L,ixO^L,pth)
        use mod_global_parameters
        use mod_small_values, only: trace_small_values

        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: pth(ixI^S)

        integer :: iw, ix^D

        {do ix^DB= ixOmin^DB,ixOmax^DB\}
            pth(ix^D)=eos%gamma_minus_1*(w(ix^D,e_)-half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_)))
            if(fix_small_values.and.pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
        {end do\}
        if(check_small_values.and..not.fix_small_values) then
            {do ix^DB= ixOmin^DB,ixOmax^DB\}
                if(pth(ix^D)<small_pressure) then
                    write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                         " encountered when call mhd_get_pthermal_hde"
                    write(*,*) "Iteration: ", it, " Time: ", global_time
                    write(*,*) "Location: ", x(ix^D,:)
                    write(*,*) "Cell number: ", ix^D
                    do iw=1,nw
                        write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
                    end do
                    if(trace_small_values) write(*,*) sqrt(pth(ix^D)-bigdouble)
                    write(*,*) "Saving status at the previous time step"
                    crash=.true.
                end if
            {end do\}
        end if

    end subroutine mhd_get_pthermal_hde

    !> Temperature routines (moved from mod_mhd_phys.t)

    !> Copy temperature from stored Te variable
    subroutine mhd_get_temperature_from_Te(w, x, ixI^L, ixO^L, res)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: res(ixI^S)
        res(ixO^S) = w(ixO^S, Te_)
    end subroutine mhd_get_temperature_from_Te

    !> Calculate temperature=p/rho when in e_ the internal energy is stored
    subroutine mhd_get_temperature_from_eint(w, x, ixI^L, ixO^L, res)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: res(ixI^S)

        double precision :: R(ixI^S)

        call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
        res(ixO^S) = eos%gamma_minus_1 * w(ixO^S, e_)/(w(ixO^S,rho_)*R(ixO^S))
    end subroutine mhd_get_temperature_from_eint

    !> Calculate temperature=p/rho when in e_ the total energy is stored
    subroutine mhd_get_temperature_from_etot(w, x, ixI^L, ixO^L, res)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: res(ixI^S)

        double precision :: R(ixI^S)

        call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
        call mhd_get_pthermal(w,x,ixI^L,ixO^L,res)
        res(ixO^S)=res(ixO^S)/(R(ixO^S)*w(ixO^S,rho_))

    end subroutine mhd_get_temperature_from_etot

    subroutine mhd_get_temperature_from_etot_with_equi(w, x, ixI^L, ixO^L, res)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: res(ixI^S)

        double precision :: R(ixI^S)

        call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
        call mhd_get_pthermal(w,x,ixI^L,ixO^L,res)
        res(ixO^S)=res(ixO^S)/(R(ixO^S)*(w(ixO^S,rho_)+block%equi_vars(ixO^S,equi_rho0_,b0i)))

    end subroutine mhd_get_temperature_from_etot_with_equi

    subroutine mhd_get_temperature_from_eint_with_equi(w, x, ixI^L, ixO^L, res)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: res(ixI^S)

        double precision :: R(ixI^S)

        call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
        res(ixO^S) = (eos%gamma_minus_1 * w(ixO^S, e_) + block%equi_vars(ixO^S,equi_pe0_,b0i)) /&
                    ((w(ixO^S,rho_) +block%equi_vars(ixO^S,equi_rho0_,b0i))*R(ixO^S))

    end subroutine mhd_get_temperature_from_eint_with_equi

    !> LTE+MHD: T from total energy. Subtract KE+ME (via the MHD-aware
    !> phys_e_to_ei dispatcher) to get eint, then look up T in the LTE table
    !> through eos%get_temperature_from_eint (bound to get_temperature_from_eint_LTE
    !> by eos_finalise). Always fresh; does not use the iw_te cache.
    !> Should NOT be bound when mhd_internal_e=.true. (w(:,e_) is already eint
    !> and phys_e_to_ei is unbound) - use eos%get_temperature_from_eint directly.
    subroutine mhd_get_temperature_from_etot_LTE(w, x, ixI^L, ixO^L, res)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: res(ixI^S)

        double precision :: wlocal(ixI^S, 1:nw)

        wlocal(ixI^S, 1:nw) = w(ixI^S, 1:nw)
        ! Subtract KE+ME from w(:,e_) via the MHD-aware dispatcher
        call phys_e_to_ei(ixI^L, ixO^L, wlocal, x)
        ! Now wlocal(:,e_) is the gas internal energy; LTE table lookup
        call eos%get_temperature_from_eint(wlocal, x, ixI^L, ixO^L, res)

    end subroutine mhd_get_temperature_from_etot_LTE

    subroutine mhd_get_temperature_equi(w,x, ixI^L, ixO^L, res)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: res(ixI^S)

        double precision :: R(ixI^S)

        call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
        res(ixO^S)= block%equi_vars(ixO^S,equi_pe0_,b0i)/(block%equi_vars(ixO^S,equi_rho0_,b0i)*R(ixO^S))

    end subroutine mhd_get_temperature_equi

    !> Equilibrium extraction routines (moved from mod_mhd_phys.t)

    subroutine mhd_get_rho_equi(w, x, ixI^L, ixO^L, res)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: res(ixI^S)
        res(ixO^S) = block%equi_vars(ixO^S,equi_rho0_,b0i)
    end subroutine mhd_get_rho_equi

    subroutine mhd_get_pe_equi(w,x, ixI^L, ixO^L, res)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: res(ixI^S)
        res(ixO^S) = block%equi_vars(ixO^S,equi_pe0_,b0i)
    end subroutine mhd_get_pe_equi

    !> Electron and hydrogen number densities of the background alone, for the
    !> equilibrium rate that subtract_equi removes. Splitting is restricted to
    !> the FI EoS, so the electron count is the fully ionised one.
    subroutine mhd_get_ne_nH_equi(ixI^L, ixO^L, w, x, ne, nH)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: ne(ixI^S), nH(ixI^S)
        nH(ixO^S) = block%equi_vars(ixO^S,equi_rho0_,b0i) / eos%nH2rhoFactor
        ne(ixO^S) = nH(ixO^S) * eos%neOnH_FI
    end subroutine mhd_get_ne_nH_equi

    !> Internal energy extraction + small value handling (moved from mod_mhd_phys.t)

    !> Get internal energy from conserved state (origin formulation)
    function mhd_get_ei_origin(w, ixI^L, ixO^L) result(ei)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, nw)
        double precision             :: ei(ixO^S)

        ! ei = E_total - e_kinetic - e_magnetic
        ei(ixO^S) = w(ixO^S,e_) - half*((^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_) &
             + (^C&w(ixO^S,b^C_)**2+))
    end function mhd_get_ei_origin

    !> Get internal energy from conserved state (internal energy formulation)
    function mhd_get_ei_inte(w, ixI^L, ixO^L) result(ei)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, nw)
        double precision             :: ei(ixO^S)

        ! e_ already stores internal energy
        ei(ixO^S) = w(ixO^S,e_)
    end function mhd_get_ei_inte

    ! mhd_handle_small_ei stays in mod_mhd_phys.t (callers are there)

    !> Update temperature Te_ from pressure using partial ionization

    !> PI energy-mode routines (eos_type='PI', ionE=.true.)
    !> Structural mirror of the FI conversion/pthermal/csound family, sharing
    !> FI's eq_state_units / RR=1 normalisation; the only difference is that the
    !> eint<->p relation is delegated to the portable scalar backend
    !> (mod_eos_PI), so eint carries the ionisation-energy term.
    !> Origin: w(e_) is total energy; inte: w(e_) is gas internal energy.

    !> Primitive pressure -> total energy (origin). m^C_ is still velocity here
    !> (momentum conversion follows in to_conserved), so KE = half*rho*v^2.
    subroutine mhd_p_to_e_PI(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: eint
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            call eint_from_rho_p_PI(w(ix^D,rho_), w(ix^D,p_), eint)
            w(ix^D,e_)=eint &
                      +half*((^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)&
                      +(^C&w(ix^D,b^C_)**2+))
        {end do\}

    end subroutine mhd_p_to_e_PI

    !> Primitive pressure -> gas internal energy (inte formulation)
    subroutine mhd_p_to_eint_PI(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: eint
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            call eint_from_rho_p_PI(w(ix^D,rho_), w(ix^D,p_), eint)
            w(ix^D,e_)=eint
        {end do\}

    end subroutine mhd_p_to_eint_PI

    !> Primitive -> conserved (origin, total energy)
    subroutine mhd_to_conserved_origin_PI(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        call mhd_p_to_e_PI(ixI^L, ixO^L, w, x)
        ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,m^C_)\

    end subroutine mhd_to_conserved_origin_PI

    !> Primitive -> conserved (internal energy)
    subroutine mhd_to_conserved_inte_PI(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        call mhd_p_to_eint_PI(ixI^L, ixO^L, w, x)
        ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,m^C_)\

    end subroutine mhd_to_conserved_inte_PI

    !> Conserved -> primitive (origin). Subtract KE+ME -> eint, invert -> p.
    subroutine mhd_to_primitive_origin_PI(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: inv_rho, eint_val, T, Rfac
        integer :: ix^D

        if (fix_small_values) then
            call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, &
                'mhd_to_primitive_origin_PI')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0/w(ix^D,rho_)
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
            ! eint = E - 0.5*rho*v^2 - 0.5*B^2
            eint_val = w(ix^D,e_) &
                      - half*(w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+) &
                      + (^C&w(ix^D,b^C_)**2+))
            eint_val = max(eint_val, smalldouble)
            call state_from_eint_PI(w(ix^D,rho_), eint_val, T, w(ix^D,p_), Rfac)
        {end do\}

    end subroutine mhd_to_primitive_origin_PI

    !> Conserved -> primitive (internal energy). e_ already holds eint.
    subroutine mhd_to_primitive_inte_PI(ixI^L,ixO^L,w,x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: inv_rho, eint_val, T, Rfac
        integer :: ix^D

        if (fix_small_values) then
            call phys_handle_small_values(.false., w, x, ixI^L, ixO^L, &
                'mhd_to_primitive_inte_PI')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            eint_val = max(w(ix^D,e_), smalldouble)
            call state_from_eint_PI(w(ix^D,rho_), eint_val, T, w(ix^D,p_), Rfac)
            ! Convert momentum to velocity
            inv_rho = 1.d0/w(ix^D,rho_)
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
        {end do\}

    end subroutine mhd_to_primitive_inte_PI

    !> Thermal pressure (origin): subtract KE+ME, invert eint -> p.
    subroutine mhd_get_pthermal_origin_PI(w,x,ixI^L,ixO^L,pth)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: pth(ixI^S)

        double precision :: eint_val, T, Rfac
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            eint_val = w(ix^D,e_) &
                      - half*((^C&w(ix^D,m^C_)**2+)/w(ix^D,rho_) &
                      + (^C&w(ix^D,b^C_)**2+))
            eint_val = max(eint_val, smalldouble)
            call state_from_eint_PI(w(ix^D,rho_), eint_val, T, pth(ix^D), Rfac)
            if(fix_small_values.and.pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
        {end do\}

    end subroutine mhd_get_pthermal_origin_PI

    !> Thermal pressure (inte): e_ holds eint, invert eint -> p.
    subroutine mhd_get_pthermal_inte_PI(w,x,ixI^L,ixO^L,pth)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: pth(ixI^S)

        double precision :: eint_val, T, Rfac
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            eint_val = max(w(ix^D,e_), smalldouble)
            call state_from_eint_PI(w(ix^D,rho_), eint_val, T, pth(ix^D), Rfac)
            if(fix_small_values.and.pth(ix^D)<small_pressure) pth(ix^D)=small_pressure
        {end do\}

    end subroutine mhd_get_pthermal_inte_PI

    !> Adiabatic sound speed squared from primitive (rho, p).

    !> Effective Gamma1 = cs2 * rho / p (for the same primitive state).

    !> PI energy mode: refresh Te_ from gas internal energy via the backend
    !> eint->T inversion (direct, so no lagged wCT -- unlike no-energy
    !> mhd_update_temperature which lags iz_H via wCT(Te_)). Uses MHD's Te_ index
    !> (PI registers Te via var_set_auxvar, so the generic iw_te is unset).

    !> PI prominence (T,p) routines (eos_type='PI', pi_table='prominence').
    !> Prominence ionisation depends on (T,p), not T alone, so the chromosphere
    !> T-only seam does not apply. Following the collaborator's design, R and T
    !> are recomputed from (rho, pth) via ionization_get_state, which dispatches
    !> the prominence (pressure-bisection) inversion internally. No-energy only
    !> (guarded in eos_validate_params), so pth = (gamma-1)*eint is R-independent
    !> and obtained from phys_get_ei (conserved w; KE/ME removed).

    !> R-factor for the prominence (T,p) table: recompute from (rho, pth).

    !> Refresh stored Te_ for the prominence table from (rho, pth).

end module mod_mhd_eos
!> Needs a line after to pass the preprocesor
