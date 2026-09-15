!=============================================================================
!> PI (partial-ionisation, eos_type='PI') arm of the eos% family.
!>
!> The adapter between the ionisation backend (mod_eos_pi_tables) and the eos%
!> authority. PI rides the FI fully-ionised ideal-gas closure (RR=1
!> normalisation, p=(gamma-1) eint when ionE=F) and differs only through the
!> variable mean-mass that enters via eos%get_Rfactor and -- in energy mode --
!> the ionisation energy folded into eint. mod_eos_pi_tables has no mod_eos
!> dependency, so eos% -> ionisation is acyclic.
!>
!> Structure mirrors mod_eos_fi / mod_eos_lte:
!>   A. init / finalise               -- wire the eos% pointer targets
!>   B. generic block getters         -- the eos%-interface routines bound for PI
!>                                       (csound2, gamma1, Rfactor, Te update);
!>                                       module-agnostic via iw_rho/iw_e/iw_te,
!>                                       so the hd/mhd/ffhd seams share ONE copy
!>                                       instead of three. The seam keeps only
!>                                       the B/KE-specific conversions.
!>   C. scalar combined-solve backend -- one ionisation solve yields T,p,R(,iz)
!>                                       together; called per-cell by the seam
!>                                       conversions and by section D.
!>   D. fl-port scalar shims          -- f(log_nH,log_x) callbacks for the
!>                                       cooling / conduction port objects.
!> mod_eos re-exports the public names so the seams reach them through the single
!> `use mod_eos` facade, exactly as for the FI/LTE kernels.
!=============================================================================
module mod_eos_pi
    use mod_global_parameters
    use mod_eos_container, only: eos, EOS_TYPE_PI
    use mod_eos_fi, only: get_temperature_from_eint_FI, &
         get_temperature_from_pressure_FI
    use mod_eos_pi_tables
    use mod_timing
    use mod_comm_lib, only: mpistop
    implicit none
    private

    !> Lifecycle (PI arms of the eos_init / eos_finalise dispatchers). The eos%
    !> block getters (update_eos_PI, get_Te_PI, get_csound2_PI, get_Rfactor_*_PI)
    !> are PRIVATE: they are bound to eos% pointers inside eos_finalise_PI, so no
    !> external module needs to name them.
    public :: eos_init_PI, eos_finalise_PI

    !> get_gamma1_PI is bound to the phys_get_gamma1 pointer in the seam (the
    !> phys layer owns that pointer), so it must stay public.
    public :: get_gamma1_PI

    !> Scalar combined-solve backend, called per-cell by the seam conversions
    !> (rho = nH = 1 per H under the FI normalisation; eint is the GAS internal
    !> energy, carrying the ionisation energy in energy mode).
    public :: state_from_eint_PI   !> rho,eint -> T, p, R   (energy mode)
    public :: p_eint_from_rho_T_PI !> rho,T    -> p, eint, R
    public :: eint_from_rho_p_PI   !> rho,p    -> eint      (prim -> cons)
    public :: csound2_prim_PI      !> rho,p    -> csound2

    !> fl-port scalar callbacks (signature f(log_nH, log_x) -> scalar) for the
    !> cooling / conduction port objects. PI's T-only ionisation is
    !> nH-independent, so log_nH is ignored and quantities are evaluated per H.
    public :: eint_from_T_PI       !> (log_nH, log_T)      -> eint/nH
    public :: p2eint_PI            !> (log_nH, log_p_nH)   -> eint/p factor
    public :: T_from_eint_PI       !> (log_nH, log_eint_nH)-> T
    public :: y_from_eint_PI       !> (log_nH, log_eint_nH)-> ne/nH

contains

    !=========================================================================
    ! A. Lifecycle
    !=========================================================================

    !> PI arm of eos_init (before units are known): PI shares the FI
    !> temperature-from-pressure getter -- T=p/(R*rho) routed through
    !> eos%get_Rfactor, which the seam points at the ionisation R-factor.
    subroutine eos_init_PI()
        eos%get_temperature_from_pressure => get_temperature_from_pressure_FI
    end subroutine eos_init_PI

    !> PI arm of eos_finalise (after unit_* are set): ride FI's ideal-gas getter
    !> set and RR=1 normalisation (the variable mean-mass enters only via
    !> eos%get_Rfactor, wired by the seam link arm), then initialise the
    !> ionisation backend. In energy mode eint carries the ionisation energy, so
    !> swap in the backend eint->T inversion.
    subroutine eos_finalise_PI()
        use mod_global_parameters, only: mype
        double precision :: Rfactor_norm

        eos%update_eos => update_eos_PI   !> direct Te cache refresh each substep
        eos%get_temperature_from_eint => get_temperature_from_eint_FI
        eos%get_Te => get_Te_PI           !> read the cached Te (refreshed above)

        !> R-factor: dispatch on pi_table. Physics-independent, so bound here
        !> (was in the hd/mhd/ffhd link arms) -> these targets stay private.
        if (eos%pi_table == 'prominence') then
            eos%get_Rfactor => get_Rfactor_prominence_PI
        else
            eos%get_Rfactor => get_Rfactor_tonly_PI
        end if
        !> Fully-ionised particle counts (2 + 3*A_He per H; ne/nH = 1 + 2*A_He)
        eos%n_per_nH_FI = 2.0d0 + 3.0d0 * eos%He_abundance
        eos%neOnH_FI    = 1.0d0 + 2.0d0 * eos%He_abundance

        !> Initialise the ionisation backend. Under the FI normalisation the
        !> (a,b) factors are absorbed into unit_*, RR=1, and b=(2+3 A_He) is the
        !> full-ionisation reference in unit_pressure; the module's
        !> R = (1+iz_H+A_He(...))/Rfactor_norm must reduce to 1 at full ionisation,
        !> so Rfactor_norm = 2+3 A_He (the LTE a=b=1 value 1+4 A_He would be wrong).
        Rfactor_norm = 2.0d0 + 3.0d0 * eos%He_abundance
        call ionization_degree_init(eos%He_abundance, Rfactor_norm, &
             trim(eos%pi_table), include_energy=eos%ionE)
        if (mype == 0) write(*,*) "[eos PI] ionisation backend initialised: table=", &
             trim(eos%pi_table), " ionE=", eos%ionE, " Rfactor_norm=", Rfactor_norm

        !> Energy mode: eint carries the ionisation energy, so swap in the
        !> backend eint->T inversion and the energy-mode csound2 (both
        !> physics-independent; get_csound2_PI was in the link arms).
        if (eos%ionE) then
            eos%get_temperature_from_eint => get_temperature_from_eint_PI
            eos%get_csound2 => get_csound2_PI
        end if
    end subroutine eos_finalise_PI

    !=========================================================================
    ! B. Generic block getters (eos%-interface shape).
    !    Module-agnostic: the primitive pressure is w(iw_e) (e_ aliases p_ in
    !    primitive state) and the cached temperature is w(iw_te), published by
    !    the phys module. One copy serves hd, mhd and ffhd.
    !=========================================================================

    !> Adiabatic sound speed squared from primitive (rho, p). PI energy mode.
    subroutine get_csound2_PI(w, x, ixI^L, ixO^L, cs2)
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: cs2(ixI^S)

        integer :: ix^D

        timeeos0 = MPI_WTIME()
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            call csound2_prim_PI(w(ix^D,iw_rho), w(ix^D,iw_e), cs2(ix^D))
        {end do\}
        timeeos_csound = timeeos_csound + (MPI_WTIME()-timeeos0)
    end subroutine get_csound2_PI

    !> Effective Gamma1 = cs2 * rho / p for the same primitive state.
    subroutine get_gamma1_PI(w, x, ixI^L, ixO^L, gamma1)
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S, nw)
        double precision, intent(in) :: x(ixI^S, 1:ndim)
        double precision, intent(out):: gamma1(ixI^S)

        double precision :: cs2
        integer :: ix^D

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            call csound2_prim_PI(w(ix^D,iw_rho), w(ix^D,iw_e), cs2)
            gamma1(ix^D) = cs2 * w(ix^D,iw_rho) / w(ix^D,iw_e)
        {end do\}
    end subroutine get_gamma1_PI

    !> No-energy R-factor (chromosphere/flare T-only tables) from the stored Te_.
    !> Denominator (2+3 He) = b (full-ion reference) absorbed into unit_pressure
    !> under the FI normalisation, RR=1 -> R->1 at full ionisation.
    subroutine get_Rfactor_tonly_PI(w, x, ixI^L, ixO^L, Rfactor)
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: Rfactor(ixI^S)
        double precision :: iz_H(ixO^S), iz_He(ixO^S)

        call ionization_degree_from_temperature(ixI^L,ixO^L,w(ixI^S,iw_te),iz_H,iz_He)
        Rfactor(ixO^S) = (1.d0 + iz_H(ixO^S) + eos%He_abundance &
             * (1.d0 + iz_He(ixO^S)*(1.d0+iz_He(ixO^S)))) / (2.d0 + 3.d0*eos%He_abundance)
    end subroutine get_Rfactor_tonly_PI

    !> Prominence (T,p) R-factor: recompute from (rho, pth) via the backend's
    !> prominence inversion (no-energy only -> pth=(gamma-1)*eint, R-independent).
    subroutine get_Rfactor_prominence_PI(w, x, ixI^L, ixO^L, Rfactor)
        use mod_physics, only: phys_get_ei
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: Rfactor(ixI^S)
        double precision :: ei(ixI^S), pth(ixI^S), rho(ixI^S), T(ixI^S)

        ei(ixO^S)  = phys_get_ei(w, ixI^L, ixO^L)
        pth(ixO^S) = eos%gamma_minus_1 * ei(ixO^S)
        rho(ixO^S) = w(ixO^S,iw_rho)
        call ionization_get_state(ixI^L, ixO^L, rho, pth, T, Rfactor)
    end subroutine get_Rfactor_prominence_PI

    !> Direct (no-lag) temperature refresh: recompute the cached w(iw_te) from the
    !> CURRENT conserved state, mirroring update_eos_LTE. Energy mode inverts the
    !> gas internal energy (incl. ionisation energy) by Newton; no-energy modes
    !> form pth=(gamma-1)*eint and invert (rho,pth)->T via the backend
    !> (y-bisection for chromosphere/flare, p-bisection for prominence). Shared by
    !> update_eos_PI (eos%update_eos) and the gated source/STS call sites.
    subroutine refresh_Te_PI(ixI^L, ixO^L, w, x)
        use mod_physics, only: phys_get_ei
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw)
        double precision :: ei(ixI^S), pth(ixI^S), rho(ixI^S), T(ixI^S), Rf(ixI^S)
        double precision :: Tsc, psc, Rsc
        integer :: ix^D

        ei(ixO^S) = phys_get_ei(w, ixI^L, ixO^L)
        if (eos%ionE) then
            !> energy mode (chromosphere/flare): eint carries the ionisation
            !> energy, invert eint->T by Newton per cell.
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
                call state_from_eint_PI(w(ix^D,iw_rho), ei(ix^D), Tsc, psc, Rsc)
                w(ix^D,iw_te) = Tsc
            {end do\}
        else
            !> no-energy: pth=(gamma-1)*eint; ionization_get_state inverts
            !> (rho,pth)->T (T-only y-bisection or prominence p-bisection).
            pth(ixO^S) = eos%gamma_minus_1 * ei(ixO^S)
            rho(ixO^S) = w(ixO^S,iw_rho)
            call ionization_get_state(ixI^L, ixO^L, rho, pth, T, Rf)
            w(ixO^S,iw_te) = T(ixO^S)
        end if
    end subroutine refresh_Te_PI

    !> eos%update_eos arm: per-RK-substep cache refresh (called from mod_advance,
    !> exactly as update_eos_LTE). Direct, so the cached Te never lags.
    subroutine update_eos_PI(ixI^L, ixO^L, w, x)
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw)
        call refresh_Te_PI(ixI^L, ixO^L, w, x)
    end subroutine update_eos_PI

    !> eos%get_Te arm: read the cached temperature (refreshed by update_eos_PI /
    !> the gated source/STS sites), mirroring get_Te_LTE.
    subroutine get_Te_PI(w, x, ixI^L, ixO^L, T)
        integer, intent(in)           :: ixI^L, ixO^L
        double precision, intent(in)  :: w(ixI^S,1:nw)
        double precision, intent(in)  :: x(ixI^S,1:ndim)
        double precision, intent(out) :: T(ixI^S)
        T(ixO^S) = w(ixO^S,iw_te)
    end subroutine get_Te_PI

    !> PI energy mode: temperature from GAS internal energy via the backend
    !> eint->T inversion (eint carries the ionisation-energy term, so the FI
    !> p=(gamma-1)*eint relation does not hold). w(iw_e) is assumed to be eint
    !> already (the etot hook subtracts KE/ME via phys_e_to_ei before calling).
    subroutine get_temperature_from_eint_PI(w, x, ixI^L, ixO^L, res)
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(out):: res(ixI^S)

        double precision :: pth, Rfac
        integer :: ix^D

        timeeos0 = MPI_WTIME()
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            call state_from_eint_PI(w(ix^D,iw_rho), w(ix^D,iw_e), &
                 res(ix^D), pth, Rfac)
        {end do\}
        timeeos_Tfromei=timeeos_Tfromei+(MPI_WTIME()-timeeos0)
    end subroutine get_temperature_from_eint_PI

    !=========================================================================
    ! C. Scalar combined-solve backend (eos%inv_gamma_minus_1 throughout,
    !    matching the FI ideal-gas split; the ionisation-energy term is folded
    !    in by the ionisation module when ionE=.true.).
    !=========================================================================

    !> Conserved gas internal energy -> temperature, thermal pressure, R.
    !> Inverse of p_eint_from_rho_T_PI. eint is the GAS internal energy
    !> (mechanical energy already removed by the caller).
    subroutine state_from_eint_PI(rho, eint, T, p, Rfactor, iz_H, iz_He)
        double precision, intent(in)  :: rho, eint
        double precision, intent(out) :: T, p, Rfactor
        double precision, intent(out), optional :: iz_H, iz_He

        if (eos%type_id /= EOS_TYPE_PI) call mpistop( &
            "state_from_eint_PI called outside eos_type='PI'")
        call ionization_get_state_from_eint(rho, eint, eos%inv_gamma_minus_1, &
             T, p, Rfactor, iz_H, iz_He)
    end subroutine state_from_eint_PI

    !> Temperature -> thermal pressure and GAS internal energy (incl. ion energy)
    subroutine p_eint_from_rho_T_PI(rho, T, p, eint, Rfactor, iz_H, iz_He)
        double precision, intent(in)  :: rho, T
        double precision, intent(out) :: p, eint, Rfactor
        double precision, intent(out), optional :: iz_H, iz_He

        if (eos%type_id /= EOS_TYPE_PI) call mpistop( &
            "p_eint_from_rho_T_PI called outside eos_type='PI'")
        call ionization_get_p_eint_from_rho_T(rho, T, eos%inv_gamma_minus_1, &
             p, eint, Rfactor, iz_H, iz_He)
    end subroutine p_eint_from_rho_T_PI

    !> Primitive pressure -> GAS internal energy (prim -> conserved direction):
    !> invert (rho,p)->T then map T->eint, so the forward/inverse pair stay consistent.
    subroutine eint_from_rho_p_PI(rho, p, eint)
        double precision, intent(in)  :: rho, p
        double precision, intent(out) :: eint
        double precision :: T, pcheck, Rfactor

        if (eos%type_id /= EOS_TYPE_PI) call mpistop( &
            "eint_from_rho_p_PI called outside eos_type='PI'")
        call ionization_get_state_scalar(rho, p, T, Rfactor)
        call ionization_get_p_eint_from_rho_T(rho, T, eos%inv_gamma_minus_1, &
             pcheck, eint, Rfactor)
    end subroutine eint_from_rho_p_PI

    !> Adiabatic sound speed squared from primitive (rho, p).
    subroutine csound2_prim_PI(rho, p, csound2)
        double precision, intent(in)  :: rho, p
        double precision, intent(out) :: csound2
        double precision :: T, Rfactor

        if (eos%type_id /= EOS_TYPE_PI) call mpistop( &
            "csound2_prim_PI called outside eos_type='PI'")
        call ionization_get_state_scalar(rho, p, T, Rfactor)
        call ionization_get_csound2_T(T, eos%inv_gamma_minus_1, csound2)
    end subroutine csound2_prim_PI

    !=========================================================================
    ! D. fl-port scalar callbacks (rho = nH = 1 per H). Thin shims over the
    !    section-C backend matching the cooling/conduction port signatures.
    !=========================================================================

    !> Internal energy per H from temperature: eint/nH(T).
    double precision function eint_from_T_PI(log_nH, log_T) result(eint)
        double precision, intent(in) :: log_nH, log_T
        double precision :: p, Rfactor
        call p_eint_from_rho_T_PI(1.d0, 10.d0**log_T, p, eint, Rfactor)
    end function eint_from_T_PI

    !> eint/p factor from pressure per H: maps p -> eint = p * (this).
    double precision function p2eint_PI(log_nH, log_p_nH) result(factor)
        double precision, intent(in) :: log_nH, log_p_nH
        double precision :: p, eint
        p = 10.d0**log_p_nH
        call eint_from_rho_p_PI(1.d0, p, eint)
        factor = eint / p
    end function p2eint_PI

    !> Temperature from internal energy per H.
    double precision function T_from_eint_PI(log_nH, log_eint_nH) result(T)
        double precision, intent(in) :: log_nH, log_eint_nH
        double precision :: p, Rfactor
        call state_from_eint_PI(1.d0, 10.d0**log_eint_nH, T, p, Rfactor)
    end function T_from_eint_PI

    !> Electron-to-hydrogen ratio ne/nH from internal energy per H.
    !> ne/nH = iz_H + A_He*iz_He*(1+iz_He) (matches the R-factor numerator's
    !> electron count; second He ionisation assumed equal to the first).
    double precision function y_from_eint_PI(log_nH, log_eint_nH) result(y)
        double precision, intent(in) :: log_nH, log_eint_nH
        double precision :: T, p, Rfactor, iz_H, iz_He
        call state_from_eint_PI(1.d0, 10.d0**log_eint_nH, T, p, Rfactor, iz_H, iz_He)
        y = iz_H + eos%He_abundance*iz_He*(1.d0+iz_He)
    end function y_from_eint_PI

end module mod_eos_pi
