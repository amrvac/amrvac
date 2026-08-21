"""Pure-numpy H+He Saha (no JAX) — fast and side-effect-free.

For building the table in seconds. The JAX path in saha.py is kept for
sanity-checking derivatives via autodiff later.
"""
import numpy as np

kB = 1.380649e-16; mp = 1.6726e-24; me = 9.1094e-28
h  = 6.6261e-27; eV = 1.602e-12
chi_H    = 13.598 * eV
chi_HeI  = 24.587 * eV
chi_HeII = 54.418 * eV
HE = 0.1

# Atomic partition functions, Barklem & Collet (2016) table 8 -- the same source
# saha_multi.py reads from data/bc16_table8_atomic_partfn.dat. Over 3-25 kK these are
# constant to <1e-4 (U(H_I)=2.0000-2.0002), so they are held fixed here to keep this
# module dependency-free; saha_multi carries the full T-dependent tables.
U_HI, U_HII             = 2.0, 1.0
U_HeI, U_HeII, U_HeIII  = 1.0, 2.0, 1.0

def K_saha(T, chi, U_up=1.0, U_lo=1.0):
    """Saha ratio  n_{j+1} n_e / n_j = 2 (U_{j+1}/U_j) (2 pi m_e k T/h^2)^{3/2} e^{-chi/kT}.

    The leading 2 is the free-electron spin degeneracy (the electron's own g is applied
    separately in eos_state's Sackur-Tetrode term). The U ratio was previously OMITTED,
    i.e. implicitly 1 for every species, which made K exactly 2x too large for hydrogen
    (U_HI=2) and 2x too small for He I->II (U_HeII/U_HeI=2). Verified against saha_multi
    (independent implementation + B&C16 data): agreement 0.07% over x_H = 1e-4 .. 0.55."""
    return 2.0 * (U_up/U_lo) * (2.0*np.pi*me*kB*T/h**2)**1.5 * np.exp(-chi/(kB*T))

def saha_HHe_solve(rho, T):
    nH = rho / ((1.0 + 4.0*HE) * mp)
    K_H    = K_saha(T, chi_H,    U_HII,   U_HI)
    K_HeI  = K_saha(T, chi_HeI,  U_HeII,  U_HeI)
    K_HeII = K_saha(T, chi_HeII, U_HeIII, U_HeII)
    ne = nH * 0.5
    for _ in range(200):
        y_H = K_H / (K_H + ne)
        r1 = K_HeI / max(ne, 1e-30)
        r2 = K_HeII / max(ne, 1e-30)
        f_HeI = 1.0 / (1.0 + r1 + r1*r2)
        y_HeII = r1 * f_HeI
        y_HeIII = r1 * r2 * f_HeI
        ne_new = y_H*nH + (y_HeII + 2.0*y_HeIII) * HE * nH
        if abs(ne_new - ne) < 1e-14 * max(ne_new, 1e-30):
            ne = ne_new; break
        ne = 0.5*(ne + ne_new)
    # Final populations at converged ne
    y_H = K_H / (K_H + ne)
    r1 = K_HeI / ne; r2 = K_HeII / ne
    f_HeI = 1.0 / (1.0 + r1 + r1*r2)
    y_HeII = r1*f_HeI
    y_HeIII = r1*r2*f_HeI
    return y_H, y_HeII, y_HeIII, ne/nH

def eos_state(rho, T):
    nH = rho / ((1.0 + 4.0*HE)*mp)
    nHe = HE*nH
    y_H, y_HeII, y_HeIII, ne_over_nH = saha_HHe_solve(rho, T)
    ne = ne_over_nH*nH
    nHI   = (1.0 - y_H) * nH
    nHII  = y_H * nH
    nHeI  = (1.0 - y_HeII - y_HeIII) * nHe
    nHeII = y_HeII * nHe
    nHeIII= y_HeIII * nHe
    n_tot = nHI + nHII + nHeI + nHeII + nHeIII + ne
    p = n_tot * kB * T
    e_int = (1.5 * n_tot * kB * T
             + chi_H * nHII
             + chi_HeI * nHeII
             + (chi_HeI + chi_HeII) * nHeIII)
    # Per-species Sackur-Tetrode entropy density
    def stran(n, mass, g=2.0):
        n_safe = max(n, 1e-30)
        return n_safe * kB * (2.5 + np.log(
            g * (2.0*np.pi*mass*kB*T/h**2)**1.5 / n_safe))
    # Internal degeneracies must be the SAME partition functions the Saha ratios use,
    # or the EoS is not thermodynamically self-consistent. g_HI was 4.0 here while the
    # Saha implied U_HI=1 -- inconsistent under either convention (electronic-only is
    # (H_I,H_II)=(2,1); including the proton's nuclear spin would be (4,2), never (4,1)).
    # Nuclear spin cancels in the Saha ratio and, applied consistently, contributes only
    # a constant n_H,tot k ln2; applied to H I alone it varies with ionisation and
    # corrupts the adiabat.
    s = (stran(nHI,   mp,   g=U_HI)   + stran(nHII,  mp,   g=U_HII) +
         stran(nHeI,  4*mp, g=U_HeI)  + stran(nHeII, 4*mp, g=U_HeII) +
         stran(nHeIII,4*mp, g=U_HeIII)+ stran(ne,    me,   g=2.0))
    return p, e_int, s

def eos_at_rho_eint(rho, e_int_target):
    """Bisection on log T such that eos_state(rho, T)[1] == e_int_target."""
    logT_lo, logT_hi = np.log(100.0), np.log(2e6)
    for _ in range(80):
        logT_mid = 0.5*(logT_lo + logT_hi)
        T_mid = np.exp(logT_mid)
        _, e_int, _ = eos_state(rho, T_mid)
        if e_int < e_int_target: logT_lo = logT_mid
        else:                    logT_hi = logT_mid
    T = np.exp(0.5*(logT_lo + logT_hi))
    p, e_int, s = eos_state(rho, T)
    return T, p, s
