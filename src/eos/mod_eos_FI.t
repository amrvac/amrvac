!=============================================================================
!> FI (fully-ionised, constant-gamma) EoS kernels for the eos% family.
!>
!> Carved out of mod_eos.t (the per-type split). The fully-ionised ideal gas:
!> p = (gamma-1) eint, constant Gamma_1 = gamma, temperature T = p/(R rho) with
!> R the constant ionisation R-factor. update_eos is a no-op (nothing cached).
!> mod_eos re-exports the public names so existing `use mod_eos` callers are
!> unaffected; eos_finalise wires these as the FI pointer targets.
!=============================================================================
module mod_eos_FI
    use mod_global_parameters
    use mod_eos_container
    use mod_eos_shared_functions, only: get_rho
    use mod_timing

    implicit none
    private

    !> get_Te_FI is PRIVATE: bound to eos%get_Te inside eos_finalise_FI only.
    !> The rest stay public: update_eos_FI (ffhd_phys), get_gamma1_FI (the seam
    !> phys_get_gamma1), get_temperature_from_{eint,pressure}_FI (reused by
    !> mod_eos_PI), and the init/finalise dispatcher arms.
    public :: update_eos_FI
    public :: get_gamma1_FI
    public :: get_temperature_from_eint_FI
    public :: get_temperature_from_pressure_FI
    public :: eos_init_FI
    public :: eos_finalise_FI

contains

    !> FI arm of eos_init (before units are known): wire the FI
    !> temperature-from-pressure getter, which must be live before
    !> (m)hd_link_eos captures it into the radiation fluid object.
    subroutine eos_init_FI()
        eos%get_temperature_from_pressure => get_temperature_from_pressure_FI
    end subroutine eos_init_FI

    !> FI arm of eos_finalise: wire the ideal-gas getter set and the fully-ionised
    !> particle counts (2 + 3*A_He per H; ne/nH = 1 + 2*A_He). No tables to load.
    subroutine eos_finalise_FI()
        eos%update_eos => update_eos_FI   !> nothing cached for ideal gas
        eos%get_temperature_from_eint => get_temperature_from_eint_FI
        eos%get_Te => get_Te_FI
        eos%n_per_nH_FI = 2.0d0 + 3.0d0 * eos%He_abundance
        eos%neOnH_FI    = 1.0d0 + 2.0d0 * eos%He_abundance
    end subroutine eos_finalise_FI

    !> FI update_eos: nothing to cache for an ideal gas. Called each RK substep.
    subroutine update_eos_FI(ixI^L, ixO^L, w, x)
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(inout) :: w(ixI^S,1:nw)
    end subroutine update_eos_FI

    !> Gamma_1 for the fully-ionised ideal gas: constant eos%gamma.
    subroutine get_gamma1_FI(w, x, ixI^L, ixO^L, gamma1)
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: w(ixI^S, nw)
        double precision, intent(in)    :: x(ixI^S, 1:ndim)
        double precision, intent(out)   :: gamma1(ixI^S)

        gamma1(ixO^S) = eos%gamma

    end subroutine get_gamma1_FI

    !> Cached-Te accessor for FI: T = pth/(rho*R) (no stored Te_; recompute).
    subroutine get_Te_FI(w,x,ixI^L,ixO^L,T)
        integer, intent(in)           :: ixI^L, ixO^L
        double precision, intent(in)  :: w(ixI^S,1:nw)
        double precision, intent(in)  :: x(ixI^S,1:ndim)
        double precision, intent(out) :: T(ixI^S)
        double precision :: Rfactor(ixI^S), pth(ixI^S), rho(ixI^S)

        !> No timing here: get_thermal_pressure is already timed.
        !> The division and Rfactor call are trivial.
        !> get_thermal_pressure already returns the total pressure under
        !> equilibrium splitting, so only the density needs the background.
        call eos%get_thermal_pressure(w, x, ixI^L, ixO^L, pth)
        call eos%get_Rfactor(w,x,ixI^L,ixO^L,Rfactor)
        call get_rho(w,x,ixI^L,ixO^L,rho)
        T(ixO^S) = pth(ixO^S) / (rho(ixO^S) * Rfactor(ixO^S))

    end subroutine get_Te_FI

    !> Temperature from internal energy, fully-ionised ideal gas: T = pth/(rho*R).
    subroutine get_temperature_from_eint_FI(w, x, ixI^L, ixO^L, res)
        !> Assumes input energy is internal energy
        use mod_physics
        integer, intent(in)             :: ixI^L,ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(in) :: w(ixI^S,1:nw)
        double precision, intent(out)   :: res(ixI^S)

        double precision :: Rfactor(ixI^S), rho(ixI^S), pth(ixI^S)

        timeeos0 = MPI_WTIME()

        call eos%get_Rfactor(w,x,ixI^L,ixO^L,Rfactor)
        call get_rho(w,x,ixI^L,ixO^L,rho)
        !> Under equilibrium splitting the e slot holds the perturbation internal
        !> energy, so the background pressure is added back to recover the total.
        pth(ixO^S) = eos%gamma_minus_1 * w(ixO^S,iw_e)
        if (iw_equi_p > 0) pth(ixO^S) = pth(ixO^S) &
             + block%equi_vars(ixO^S,iw_equi_p,b0i)
        res(ixO^S) = pth(ixO^S) / (Rfactor(ixO^S) * rho(ixO^S))

        timeeos_Tfromei=timeeos_Tfromei+(MPI_WTIME()-timeeos0)
    end subroutine get_temperature_from_eint_FI

    !> FI temperature from primitive pressure: T = p / (R * rho).
    !> w is PRIMITIVE here, so iw_e holds the thermal pressure. Lives in the EoS
    !> so hd and mhd share one routine (replaces the orphaned, never-called
    !> hd_get_temperature_from_prim that used to sit in mod_hd_phys).
    subroutine get_temperature_from_pressure_FI(w, x, ixI^L, ixO^L, res)
        integer, intent(in)             :: ixI^L, ixO^L
        double precision, intent(in)    :: x(ixI^S,1:ndim)
        double precision, intent(in)    :: w(ixI^S,1:nw)
        double precision, intent(out)   :: res(ixI^S)

        double precision :: Rfactor(ixI^S)

        call eos%get_Rfactor(w, x, ixI^L, ixO^L, Rfactor)
        res(ixO^S) = w(ixO^S,iw_e) / (Rfactor(ixO^S) * w(ixO^S,iw_rho))
    end subroutine get_temperature_from_pressure_FI

end module mod_eos_FI
