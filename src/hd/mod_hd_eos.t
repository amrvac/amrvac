!=============================================================================
!> HD <-> EoS seam: binds the eos% authority into hydrodynamics.
!>
!> hd_link_eos wires the eos%/phys_ procedure pointers (conversions, pthermal,
!> csound2/gamma1, get_Rfactor, temperature, update_temperature) per eos_type
!> (FI / LTE / PI, with PI energy- and prominence-mode overrides).
!> bind_eos_to_source wires the thermal-conduction / radiative-cooling / FLD
!> fluid-port callbacks from eos%.
!>
!> Holds the HD block routines that wrap the shared thermodynamics with HD's
!> mechanical-energy bookkeeping -- total-energy ("origin") formulation,
!> dust-aware: FI/LTE/PI conversions, p_to_e, pthermal, csound2/gamma1, and the
!> PI energy + prominence R-factor / temperature kernels.
!=============================================================================
module mod_hd_eos
    use mod_global_parameters
    use mod_physics
    use mod_eos
    !> Mode-specific kernels come from their sub-modules (the facade no longer
    !> re-exports them); each mpistops if called under the wrong eos_type/method.
    use mod_eos_LTE
    use mod_eos_LTE_entropy, only: entropy_eint_from_p_bisect
    use mod_eos_LTE_saha
    use mod_eos_PI
    use mod_eos_container
    use mod_hd_phys
    use mod_timing
    use mod_radiative_cooling, only: build_Y_mod_table

    use mod_comm_lib, only: mpistop

    implicit none
    private

    procedure(sub_convert), pointer, public :: hd_to_primitive  => null()
    procedure(sub_convert), pointer, public :: hd_to_conserved  => null()

    public :: hd_link_eos, hd_get_pthermal

contains

    !> Link the appropriate EOS conversion routines based on the selected EoS type
    subroutine hd_link_eos()
        use mod_usr_methods

        ! PI (partial approximations) takes the FI ideal-gas conversions as
        ! its BASE (e_int=p/(gamma-1); variable R via get_Rfactor). Energy
        ! mode (ionE) overrides them with the PI-energy variants at the end
        ! of this routine; FI conversions stay pristine.
        if (eos%eos_type == 'FI' .or. eos%eos_type == 'PI') then
            eos%to_conserved => hd_to_conserved_origin
            eos%to_primitive => hd_to_primitive_origin
        else if (eos%eos_type == 'LTE') then
            eos%to_conserved => hd_to_conserved_LTE
            eos%to_primitive => hd_to_primitive_LTE
        else
            call mpistop('Error: Unknown HD EOS type: ' // trim(eos%eos_type))
        end if

        phys_to_primitive => eos%to_primitive
        phys_to_conserved => eos%to_conserved
        phys_get_rho      => eos%get_rho
        phys_get_pthermal => hd_get_pthermal
        eos%get_thermal_pressure => hd_get_pthermal
        phys_bind_eos_to_source  => bind_eos_to_source

        eos%p_to_e => p_to_e !> suitable for both FI and LTE

        ! Link sound speed and gamma1 computation
        if (eos%eos_type == 'LTE' .and. eos%ionE) then
            eos%get_csound2 => hd_get_csound2_LTE
            phys_get_gamma1 => hd_get_gamma1_LTE
            ! EoS-aware prolongation: interpolate in (rho, v, T) space
            phys_to_prolong   => hd_to_prolong_LTE
            phys_from_prolong => hd_from_prolong_LTE
        else
            eos%get_csound2 => hd_get_csound2_FI
            phys_get_gamma1 => get_gamma1_FI
        end if

        ! Rfactor: only the FI case is physics-dependent (the usr_Rfactor user
        ! hook). LTE and PI bind pure eos% routines in eos_finalise_{LTE,PI},
        ! which run after this and before bind_eos_to_source (so they win and
        ! their targets stay private). Note: this makes usr_Rfactor apply only to
        ! FI for HD -- consistent with MHD, where the EoS owns R for LTE/PI.
        if(associated(usr_Rfactor)) then
            eos%get_Rfactor=>usr_Rfactor
        else
            eos%get_Rfactor=>Rfactor_from_constant_ionization
        end if

        !> PI energy-mode overrides (eos_type='PI', ionE=.true.)
        ! Energy mode differs from no-energy PI only in the eint<->p relation,
        ! so override exactly the routines that touch it (conversions, p_to_e,
        ! csound2/gamma1, Te update), all via the portable backend. pthermal
        ! is unchanged: hd_get_pthermal already carries a PI-energy branch.
        if (eos%eos_type == 'PI' .and. eos%ionE) then
            if (.not. hd_energy) &
                call mpistop('PI energy EoS requires hd_energy=.true.')
            eos%to_conserved  => hd_to_conserved_PI
            eos%to_primitive  => hd_to_primitive_PI
            phys_to_primitive => eos%to_primitive
            phys_to_conserved => eos%to_conserved
            eos%p_to_e        => p_to_e_PI
            !> eos%get_csound2 => get_csound2_PI set in eos_finalise_PI (private target)
            phys_get_gamma1   => get_gamma1_PI
            ! Te_ from eint each substep (uses HD's Te_ index: PI registers
            ! Te via var_set_auxvar, so the generic iw_te is unset).
        end if
        hd_to_primitive          => eos%to_primitive
        hd_to_conserved          => eos%to_conserved

    end subroutine hd_link_eos

    subroutine bind_eos_to_source() !> this is called in eos_finalise through mod_physics procedure linking
        if (allocated(tc_fl)) then
            tc_fl%get_temperature_from_conserved => eos%get_temperature_from_etot
            !> Use fast bilinear T lookup for TC STS substeps (density fixed,
            !> TC flux dominated by corona where T(eint) is smooth)
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
        end if

        if (allocated(rc_fl)) then
            rc_fl%get_rho => eos%get_rho
            rc_fl%get_pthermal => eos%get_thermal_pressure
            rc_fl%get_var_Rfactor => eos%get_Rfactor
            rc_fl%get_Te => eos%get_Te
            rc_fl%get_ne_nH => eos%get_ne_nH
            rc_fl%ionE              =  eos%ionE
            rc_fl%method            =  eos%method
            rc_fl%inv_gamma_minus_1 =  eos%inv_gamma_minus_1
            rc_fl%nH2rhoFactor      =  eos%nH2rhoFactor
            rc_fl%eion_per_nH       =  eos%eion_per_nH
            rc_fl%eint_from_T       => eint_nH_from_T
            rc_fl%p2eint            => p2eint_from_nH_p
            rc_fl%T_from_eint       => T_from_nH_eint
            rc_fl%y_from_eint       => y_from_nH_eint
            !> Build the variable-c_V Townsend Y_mod table now that all
            !> EoS tables (eint_from_T, T, neOnH) are in code units.
            !> build_Y_mod_table checks coolmethod=='exact' and .not.isPPL
            !> internally and early-returns otherwise. LTE only:
            !> eos_get_eintT_grid reads the LTE eintT grid; PI runs classical
            !> Townsend (Y_mod-for-PI is a later refinement).
            if (eos%ionE .and. eos%eos_type == 'LTE') then
                call eos_get_eintT_grid(rc_fl%Y_mod_n_nH, &
                     rc_fl%Y_mod_lg_nH_min, rc_fl%Y_mod_lg_nH_max)
                call build_Y_mod_table(rc_fl)
            end if
        end if

        !> PI energy mode: repoint cooling/conduction eint<->T callbacks at
        !> the ionisation backend (the LTE-table versions above assume LTE
        !> tables PI does not load). Classical Townsend; Y_mod not built (above).
        if (eos%eos_type == 'PI' .and. eos%ionE) then
            if (allocated(tc_fl)) tc_fl%eint_from_T => eint_from_T_PI
            if (allocated(rc_fl)) then
                rc_fl%eint_from_T => eint_from_T_PI
                rc_fl%p2eint      => p2eint_PI
                rc_fl%T_from_eint => T_from_eint_PI
                rc_fl%y_from_eint => y_from_eint_PI
            end if
        end if

        if (allocated(fld_fl)) then
            !> Radiation (FLD) fluid: gas-EoS callbacks. get_temperature_from_pressure
            !> is the FI T=p/(R*rho) routine (the LTE variant is a later pass).
            fld_fl%gamma       =  eos%gamma
            fld_fl%get_tgas    => eos%get_temperature_from_pressure
            fld_fl%get_Rfactor => eos%get_Rfactor
        end if
    end subroutine bind_eos_to_source

    !> Transform primitive variables into conservative ones
    subroutine hd_to_conserved_origin(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        use mod_dust, only: dust_to_conserved
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        integer :: ix^D

        timeeos0 = MPI_WTIME() !> For monitoring cost of eos module

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if (hd_energy) then
                ! Calculate total energy from pressure and kinetic energy
                w(ix^D,e_)=w(ix^D,p_)*eos%inv_gamma_minus_1+&
                    half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
            end if
            ! Convert velocity to momentum
            ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
        {end do\}

        if (hd_dust) then
            call dust_to_conserved(ixI^L, ixO^L, w, x)
        end if

        timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine hd_to_conserved_origin

    !> Transform conservative variables into primitive ones
    subroutine hd_to_primitive_origin(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        use mod_dust, only: dust_to_primitive
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision                :: inv_rho
        integer :: ix^D

        timeeos0 = MPI_WTIME() !> For monitoring cost of eos module

        if (fix_small_values) then
            call hd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'hd_to_primitive')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0/w(ix^D,rho_)
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
            ! Calculate pressure = (gamma-1) * (e-ek)
            if(hd_energy) then
                ! Compute pressure
                w(ix^D,p_)=(eos%gamma_minus_1)*(w(ix^D,e_)&
                    -half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+))
            end if
        {end do\}

        ! Convert dust momentum to dust velocity
        if (hd_dust) then
            call dust_to_primitive(ixI^L, ixO^L, w, x)
        end if

        timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine hd_to_primitive_origin

    !> LTE primitive -> conserved conversion.
    !>
    !> On entry: rho_ = density, m_ = velocity, p_ = pressure.
    !> On exit:  rho_ = density (unchanged), m_ = momentum, e_ = total energy.
    subroutine hd_to_conserved_LTE(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        use mod_dust, only: dust_to_conserved
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        timeeos0 = MPI_WTIME()

        call p_to_e(ixI^L, ixO^L, w, x)

        ! Convert velocity to momentum
        ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,m^C_)\

        if (hd_dust) then
            call dust_to_conserved(ixI^L, ixO^L, w, x)
        end if

        timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine hd_to_conserved_LTE

    subroutine p_to_e(ixI^L, ixO^L, w, x)
        !> Convert pressure to total energy: E = eint(rho, p) + KE.
        !>
        !> On entry: w(rho_) = density, w(m_) = velocity, w(p_) = pressure.
        !> On exit:  w(rho_) unchanged, w(m_) unchanged, w(e_) = total energy.
        !>
        !> Four paths for ionE (LTE with ionisation):
        !>   1. FI bypass (p/rho > threshold): analytic eint = p/(gamma-1) + eion*nH
        !>   2. Analytical Saha (eos_method='analytic'): direct Saha solve for T,y from p
        !>   3. WB mode (hd_well_balanced): cached bisection on forward tables
        !>   4. Standard: p2eint inverse table lookup (fast, ~0.01% round-trip error)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        integer :: ix^D
        double precision :: p_to_eint, p_over_rho
        double precision :: nH(ixI^S), nH_in(ixI^S), p_in(ixI^S)
        double precision :: log_eint_mid, eint_total
        double precision :: T_solve, y_solve, eint_nH_solve

        if (eos%ionE) then
            call eos%get_nH(w, x, ixI^L, ixO^L, nH)
            nH_in(ixO^S) = dlog10(nH(ixO^S))
            if (eos%p2eint_method /= 'bisect') then
                p_in(ixO^S) = dlog10(w(ixO^S,p_)) - nH_in(ixO^S)
            end if
        endif

        p_to_eint = eos%inv_gamma_minus_1
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if (hd_energy) then
                if (eos%ionE) then
                    p_over_rho = w(ix^D,p_) / w(ix^D,rho_)
                    if (p_over_rho > eos%p_rho_FI_threshold) then
                        !> FI bypass: exact inverse of to_primitive
                        p_to_eint = eos%inv_gamma_minus_1 &
                            + eos%eion_per_nH * nH(ix^D) / w(ix^D,p_)
                        w(ix^D,e_) = w(ix^D,p_)*p_to_eint + &
                            half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                    else if (eos%method == 'analytic') then
                        !> Analytical Saha: bisect for T from p, return eint directly
                        call saha_state_from_nH_p(nH(ix^D), w(ix^D,p_), &
                            T_solve, y_solve, eint_nH_solve)
                        w(ix^D,e_) = eint_nH_solve * nH(ix^D) + &
                            half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                    else if (eos%p2eint_method == 'bisect') then
                        !> Bisection on forward p table: finds eint such that
                        !> p(nH, eint) = p_target exactly. More accurate than
                        !> the p2eint table in the ionisation zone.
                        !> Entropy mode uses the pfwd/eintP pair: eos%log_p (used by
                        !> eint_from_p_bisect) is not loaded there.
                        if (eos%method == 'entropy') then
                            call entropy_eint_from_p_bisect( &
                                eos%pfwd, eos%pfwd_x, eos%pfwd_y, eos%pfwd_xy, &
                                eos%eintP, eos%eintP_x, eos%eintP_y, eos%eintP_xy, &
                                nH_in(ix^D), dlog10(w(ix^D,p_)) - nH_in(ix^D), &
                                log_eint_mid)
                            eint_total = nH(ix^D) * 10.0d0**log_eint_mid
                        else
                            call eint_from_p_bisect(nH_in(ix^D), &
                                dlog10(w(ix^D,p_)), log_eint_mid)
                            eint_total = nH(ix^D) * 10.0d0**log_eint_mid
                            eint_total = max(eint_total, &
                                nH(ix^D) * 10.0d0**eos%T%var2_min)
                        end if
                        w(ix^D,e_) = eint_total + &
                            half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                    else
                        !> Standard table lookup (fast, default)
                        p_to_eint = p2eint_from_nH_p(nH_in(ix^D), p_in(ix^D))
                        w(ix^D,e_) = w(ix^D,p_)*p_to_eint + &
                            half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                    end if
                else
                    w(ix^D,e_) = w(ix^D,p_)*p_to_eint + &
                        half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
                end if
            end if
        {end do\}

    end subroutine p_to_e

    !> LTE conserved -> primitive conversion.
    !>
    !> On entry: rho_ = density, m_ = momentum, e_ = total energy.
    !> On exit:  rho_ = density (unchanged), m_ = velocity, p_ = pressure.
    !>
    !> Pressure is computed energy-consistently from actual eint via EoS
    !> table lookups (T and ne/nH). Cannot use stored Ne_/Te_ because
    !> they may be stale after AMR prolongation/coarsening.
    subroutine hd_to_primitive_LTE(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        use mod_dust, only: dust_to_primitive
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision                :: inv_rho
        double precision                :: nH(ixI^S)
        double precision                :: log_nH(ixI^S)
        double precision                :: eint_val, eint_in, T_loc, y_loc
        integer :: ix^D

        timeeos0 = MPI_WTIME()

        if (fix_small_values) then
            call hd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'hd_to_primitive_LTE')
        end if

        call eos%get_nH(w, x, ixI^L, ixO^L, nH)

        ! Cache log10(nH) for all cells (used by table lookups for IonE)
        if (eos%ionE) then
            log_nH(ixO^S) = dlog10(nH(ixO^S))
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0/w(ix^D,rho_)
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
            ! Calculate pressure
            if(hd_energy) then
                if (eos%ionE) then
                    ! Energy-consistent pressure from actual eint via EoS tables.
                    ! Cannot use stored Ne_/Te_ because they may be stale after
                    ! AMR prolongation/coarsening (nonlinear EoS breaks averaging).
                    eint_val = w(ix^D,e_) - half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+)
                    ! Floor eint to prevent unphysical values
                    if (eos%method /= 'analytic') then
                        eint_val = max(eint_val, nH(ix^D) * 10.0d0**eos%T%var2_min)
                    end if
                    eint_val = max(eint_val, smalldouble)
                    if (eint_val * inv_rho > eos%eint_rho_FI_threshold) then
                        ! FI bypass: p = (gamma-1)*(eint - eion*nH)
                        w(ix^D,p_) = eos%gamma_minus_1 &
                            * (eint_val - eos%eion_per_nH * nH(ix^D))
                    else
                        ! Ionisation zone: single p/nH lookup
                        eint_in = dlog10(eint_val) - log_nH(ix^D)
                        w(ix^D,p_) = nH(ix^D) &
                            * p_nH_from_eint(log_nH(ix^D), eint_in)
                    end if
                else
                    w(ix^D,p_)=(eos%gamma_minus_1)*(w(ix^D,e_)&
                        -half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+))
                end if
            end if
        {end do\}

        ! Convert dust momentum to dust velocity
        if (hd_dust) then
            call dust_to_primitive(ixI^L, ixO^L, w, x)
        end if

        timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)

    end subroutine hd_to_primitive_LTE

    !> Calculate thermal pressure within ixO^L.
    !> For energy runs delegates to eos%get_thermal_pressure; for no-energy
    !> (isothermal) uses the adiabatic relation p = adiab * rho^gamma.
    subroutine hd_get_pthermal(w, x, ixI^L, ixO^L, pth)
        use mod_global_parameters
        use mod_physics, only: phys_get_ei
        use mod_usr_methods, only: usr_set_pthermal
        use mod_small_values, only: trace_small_values

        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, 1:nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: pth(ixI^S)
        integer                      :: iw, ix^D
        double precision :: nH(ixI^S), ei(ixI^S), Tpi, Rpi

        if (hd_energy) then
            if (eos%eos_type == 'LTE') then
                ! LTE: p = nH * (1 + He + ne/nH) * T from stored state
                call eos%get_nH(w, x, ixI^L, ixO^L, nH)
                pth(ixO^S) = nH(ixO^S) * (1.0d0 + eos%He_abundance &
                    + (w(ixO^S,iw_ne) / nH(ixO^S))) * w(ixO^S,iw_te)
            else if (eos%eos_type == 'PI' .and. eos%ionE) then
                ! PI energy: invert eint (= e - KE, via phys_get_ei) -> p with
                ! the ionisation backend (eint carries the ionisation energy).
                ei(ixO^S) = phys_get_ei(w, ixI^L, ixO^L)
                {do ix^DB=ixOmin^DB,ixOmax^DB\}
                    call state_from_eint_PI(w(ix^D,rho_), ei(ix^D), &
                        Tpi, pth(ix^D), Rpi)
                {end do\}
            else
                ! FI: p = (gamma-1) * (e - KE)
                pth(ixO^S) = eos%gamma_minus_1 * phys_get_ei(w, ixI^L, ixO^L)
            end if
        else
            if (.not. associated(usr_set_pthermal)) then
                pth(ixO^S) = hd_adiab * w(ixO^S, rho_)**eos%gamma
            else
                call usr_set_pthermal(w,x,ixI^L,ixO^L,pth)
            end if
        end if

        if (fix_small_values) then
            {do ix^DB= ixO^LIM^DB\}
                if(pth(ix^D)<small_pressure) then
                    pth(ix^D)=small_pressure
                endif
            {enddo^D&\}
        else if (check_small_values) then
            {do ix^DB= ixO^LIM^DB\}
                if(pth(ix^D)<small_pressure) then
                    write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                        " encountered when call hd_get_pthermal"
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
            {enddo^D&\}
        end if

    end subroutine hd_get_pthermal

    !> Sound speed squared for FI (fully ionized / constant gamma) EoS.
    !> Expects w in primitive form: w(p_) = pressure, w(rho_) = density.
    subroutine hd_get_csound2_FI(w, x, ixI^L, ixO^L, cs2)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision, intent(out)   :: cs2(ixI^S)

        timeeos0 = MPI_WTIME()

        cs2(ixO^S) = eos%gamma * w(ixO^S, p_) / w(ixO^S, rho_)

        timeeos_csound = timeeos_csound + (MPI_WTIME()-timeeos0)

    end subroutine hd_get_csound2_FI

    !> Sound speed squared for LTE+IonE EoS.
    !> Delegates Gamma_1 computation to hd_get_gamma1_LTE, then cs2 = Gamma_1 * p/rho.
    subroutine hd_get_csound2_LTE(w, x, ixI^L, ixO^L, cs2)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision, intent(out)   :: cs2(ixI^S)
        integer :: ix^D

        timeeos0 = MPI_WTIME()

        call hd_get_gamma1_LTE(w, x, ixI^L, ixO^L, cs2)
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cs2(ix^D) = cs2(ix^D) * w(ix^D, p_) / w(ix^D, rho_)
        {end do\}

        timeeos_csound = timeeos_csound + (MPI_WTIME()-timeeos0)

    end subroutine hd_get_csound2_LTE

    !> PI energy-mode routines (eos_type='PI', ionE=.true.)
    !> Mirror the HD FI/LTE family; the eint<->p relation is delegated to the
    !> portable scalar backend (mod_eos_PI). Same eq_state_units /
    !> RR=1 normalisation as FI. Origin (total-energy) only, as in HD. eint
    !> carries the ionisation-energy term, so p=(gamma-1)*eint no longer holds.

    !> Primitive pressure -> total energy via backend (m_ is velocity here).
    subroutine p_to_e_PI(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision :: eint
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if (hd_energy) then
                call eint_from_rho_p_PI(w(ix^D,rho_), w(ix^D,p_), eint)
                w(ix^D,e_) = eint + half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
            end if
        {end do\}
    end subroutine p_to_e_PI

    !> Primitive -> conserved (PI energy)
    subroutine hd_to_conserved_PI(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        use mod_dust, only: dust_to_conserved
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        timeeos0 = MPI_WTIME()
        call p_to_e_PI(ixI^L, ixO^L, w, x)
        ! Convert velocity to momentum
        ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,m^C_)\
        if (hd_dust) call dust_to_conserved(ixI^L, ixO^L, w, x)
        timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)
    end subroutine hd_to_conserved_PI

    !> Conserved -> primitive (PI energy): KE removed -> eint, backend -> p.
    subroutine hd_to_primitive_PI(ixI^L, ixO^L, w, x)
        use mod_global_parameters
        use mod_dust, only: dust_to_primitive
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision :: inv_rho, eint_val, T, Rfac
        integer :: ix^D

        timeeos0 = MPI_WTIME()

        if (fix_small_values) then
            call hd_handle_small_values(.false., w, x, ixI^L, ixO^L, &
                'hd_to_primitive_PI')
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0/w(ix^D,rho_)
            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
            if (hd_energy) then
                eint_val = w(ix^D,e_) &
                    - half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+)
                eint_val = max(eint_val, smalldouble)
                call state_from_eint_PI(w(ix^D,rho_), eint_val, &
                    T, w(ix^D,p_), Rfac)
            end if
        {end do\}

        if (hd_dust) call dust_to_primitive(ixI^L, ixO^L, w, x)
        timeeos_conv=timeeos_conv+(MPI_WTIME()-timeeos0)
    end subroutine hd_to_primitive_PI

    !> Adiabatic sound speed squared from primitive (rho, p).

    !> Effective Gamma1 = cs2 * rho / p.

    !> PI energy mode: refresh Te_ from gas internal energy via the backend
    !> eint->T inversion. Uses HD's Te_ index (PI registers Te via
    !> var_set_auxvar, so the generic iw_te is unset). phys_get_ei supplies
    !> eint (KE removed).

    !> PI no-energy R-factor from the stored Te_ (HD's index). Mirrors MHD's
    !> Rfactor_from_temperature_ionization. The generic mod_eos version uses
    !> the global iw_te, which PI leaves unset (var_set_auxvar) -- so HD needs
    !> its own Te_-addressed routine, exactly as MHD does.

    !> PI no-energy Te_ update (pth/(rho*R) with lagged iz_H from wCT(Te_)).
    !> HD's Te_-addressed mirror of MHD's mhd_update_temperature.

    !> Prominence (T,p) R-factor: recompute from (rho, pth) via the module's
    !> prominence inversion (no-energy only -> pth=(gamma-1)*eint, R-indep).

    !> Prominence (T,p) Te_ update from (rho, pth).

    !> Return effective adiabatic index for LTE+IonE EoS.
    !> Dispatches on gamma1_method and eos%method.
    subroutine hd_get_gamma1_LTE(w, x, ixI^L, ixO^L, gamma1)
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

    end subroutine hd_get_gamma1_LTE

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

    !> Convert conserved (rho, rho*v, E) to prolong form (rho, v, T).
    !> T is stored in the e_ slot. Interpolation in this space avoids
    !> Jensen's inequality across the ionisation plateau.
    subroutine hd_to_prolong_LTE(ixI^L, ixO^L, w, x)
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: inv_rho, eint_val, nH_val, log_nH, T_loc, y_loc
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            inv_rho = 1.d0 / w(ix^D, rho_)
            nH_val = w(ix^D, rho_) / eos%nH2rhoFactor
            log_nH = dlog10(nH_val)

            ! Convert momentum to velocity
            ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\

            ! Compute eint = E - 0.5*rho*v^2
            eint_val = w(ix^D, e_) &
                - half * w(ix^D, rho_) * (^C&w(ix^D,m^C_)**2+)
            ! Floor eint to prevent unphysical values
            if (eos%method /= 'analytic') then
                eint_val = max(eint_val, nH_val * 10.0d0**eos%T%var2_min)
            end if
            eint_val = max(eint_val, smalldouble)

            ! Convert eint to T via EoS
            if (eint_val * inv_rho > eos%eint_rho_FI_threshold) then
                ! FI bypass: T = (gamma-1)*(eint - eion*nH) / (nH * n_per_nH_FI)
                w(ix^D, e_) = eos%gamma_minus_1 &
                    * (eint_val - eos%eion_per_nH * nH_val) &
                    / (nH_val * eos%n_per_nH_FI)
            else if (eos%method == 'analytic') then
                ! Analytical Saha: solve for T from eint
                call saha_T_from_nH_eint(nH_val, &
                    eint_val / nH_val, T_loc, y_loc)
                w(ix^D, e_) = T_loc
            else
                ! Ionisation zone: T from table
                w(ix^D, e_) = T_from_nH_eint(log_nH, &
                    dlog10(eint_val) - log_nH)
            end if
        {end do\}

    end subroutine hd_to_prolong_LTE

    !> Convert prolong form (rho, v, T) back to conserved (rho, rho*v, E).
    !> T is read from the e_ slot. Uses eint_nH_from_T for back-conversion.
    subroutine hd_from_prolong_LTE(ixI^L, ixO^L, w, x)
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(inout) :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)

        double precision :: T_val, eint_val, T_FI, nH_val, log_nH, log_T_min
        integer :: ix^D

        T_FI = (eos%eint_rho_FI_threshold &
            * eos%nH2rhoFactor - eos%eion_per_nH) &
            * eos%gamma_minus_1 / eos%n_per_nH_FI

        ! Floor for log_T when calling the (rho, T) inverse table.
        ! Legacy 'tables' method populates eos%eint_from_T; entropy method
        ! populates eos%eintT. Picking the wrong container leaves var2_min = 0
        ! (uninitialised), which floors T at 10^0 = 1 code unit (= 10^6 K) —
        ! that clobbers any cold cell going through AMR prolongation.
        if (eos%method == 'entropy') then
            log_T_min = eos%eintT%var2_min
        else
            log_T_min = eos%eint_from_T%var2_min
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            T_val = w(ix^D, e_)  ! T stored in e_ slot
            nH_val = w(ix^D, rho_) / eos%nH2rhoFactor
            log_nH = dlog10(nH_val)

            if (T_val > T_FI) then
                ! FI: eint = nH*(n_per_nH*T/(gamma-1) + eion)
                eint_val = nH_val &
                    * (eos%n_per_nH_FI * T_val * eos%inv_gamma_minus_1 &
                    + eos%eion_per_nH)
            else if (eos%method == 'analytic') then
                ! Analytical Saha: eint from T directly
                eint_val = saha_eint_from_nH_T(nH_val, T_val) * nH_val
            else
                ! Ionisation zone: eint/nH from T table
                eint_val = eint_nH_from_T(log_nH, &
                    dlog10(max(T_val, 10.0d0**log_T_min))) &
                    * nH_val
            end if

            ! E = eint + 0.5*rho*v^2
            w(ix^D, e_) = eint_val &
                + half * w(ix^D, rho_) * (^C&w(ix^D,m^C_)**2+)
            ! Convert velocity to momentum
            ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
        {end do\}

    end subroutine hd_from_prolong_LTE

end module mod_hd_eos
!> Needs a line after to pass the preprocesor
