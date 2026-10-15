"""H+He gas in Saha equilibrium, in plain numpy, with exact first derivatives

Hydrogen has two ionisation stages and helium three. The state at (rho, T) comes from the
charge balance, and its first derivatives from implicit differentiation of that balance, so
generate_entropy_tables.py needs no finite differences for them. CGS units throughout
"""
import numpy as np

K_B = 1.380649e-16       # Boltzmann constant (erg/K)
M_P = 1.6726e-24         # proton mass (g)
M_E = 9.1094e-28         # electron mass (g)
H_PLANCK = 6.6261e-27    # Planck constant (erg s)
EV = 1.602e-12           # electron volt (erg)
CHI_H = 13.598 * EV      # ionisation energy of H I
CHI_HEI = 24.587 * EV    # ionisation energy of He I
CHI_HEII = 54.418 * EV   # ionisation energy of He II
A_HE = 0.1               # default helium abundance, He atoms per H atom

# Atomic partition functions, held at their ground-state values. Over 3-25 kK the tabulated ones
# are constant to 1e-4 (U(H I) = 2.0000-2.0002), so this module needs no data files
U_HI, U_HII = 2.0, 1.0
U_HEI, U_HEII, U_HEIII = 1.0, 2.0, 1.0


def saha_ratio(T, chi, u_up=1.0, u_lo=1.0):
    """Saha ratio n_{j+1} n_e / n_j = 2 (U_{j+1}/U_j) (2 pi m_e k T/h^2)^{3/2} exp(-chi/kT)

    The leading 2 is the spin degeneracy of the free electron
    """
    return (2.0 * (u_up/u_lo) * (2.0*np.pi*M_E*K_B*T/H_PLANCK**2)**1.5
            * np.exp(-chi/(K_B*T)))


def _translational_entropy(n, mass, g, T):
    """Sackur-Tetrode entropy per volume of n particles of one species with degeneracy g

    An empty stage contributes nothing, since n ln(1/n) vanishes as n goes to zero. The floor on n
    only keeps the logarithm finite
    """
    n = max(n, 1e-30)
    return n * K_B * (2.5 + np.log(g * (2.0*np.pi*mass*K_B*T/H_PLANCK**2)**1.5 / n))


class SahaHHe:
    """Saha equation of state for a mixture of a_he helium atoms per hydrogen atom"""

    def __init__(self, a_he=A_HE):
        self.a_he = a_he

    def n_h(self, rho):
        """Number density of hydrogen nuclei"""
        return rho / ((1.0 + 4.0*self.a_he) * M_P)

    def ionisation(self, rho, T):
        """Ionisation state (y_H, y_HeII, y_HeIII, ne/nH) at (rho, T)

        Newton iteration in u = ln ne on the charge balance ne = nH Z(ne), started from the
        hydrogen-only root. The residual u - ln(nH Z) has slope 1 + var/Z >= 1, so the iteration
        converges from either side
        """
        n_h = self.n_h(rho)
        k_h = saha_ratio(T, CHI_H, U_HII, U_HI)
        k_hei = saha_ratio(T, CHI_HEI, U_HEII, U_HEI)
        k_heii = saha_ratio(T, CHI_HEII, U_HEIII, U_HEII)

        def stages(ne):
            # Fractions of H II, He II and He III at electron density ne
            r1 = k_hei / ne
            r2 = k_heii / ne
            f0 = 1.0 / (1.0 + r1 + r1*r2)
            return k_h / (k_h + ne), r1*f0, r1*r2*f0

        ne = 2.0*n_h / (1.0 + np.sqrt(1.0 + 4.0*n_h/k_h))
        for _ in range(100):
            y, q1, q2 = stages(ne)
            q_he = q1 + 2.0*q2
            z = y + self.a_he*q_he
            du = (np.log(ne / (n_h*z))
                  / (1.0 + (y*(1.0 - y) + self.a_he*(q1 + 4.0*q2 - q_he**2)) / z))
            ne *= np.exp(-du)
            if abs(du) < 1.0e-14:
                break
        else:
            raise RuntimeError('SahaHHe.ionisation: no convergence at rho=%g T=%g' % (rho, T))
        y, q1, q2 = stages(ne)
        return y, q1, q2, ne/n_h

    def thermo(self, rho, T):
        """State at (rho, T) and its exact first derivatives

        With u = ln ne, a stage of charge j has a Saha weight proportional to ne^-j at fixed T,
        and d ln(weight)/d ln T = a_j at fixed ne. Over the stages of each element this gives
        d<X>/du = -cov(j, X) and d<X>/d ln T = cov(a, X). The charge balance ne = nH Z(u, T)
        then gives du/d ln rho = 1/D and du/d ln T = Z_T/(Z D), with D = 1 - Z_u/Z

        Keys of the returned dict, all derivatives taken per ln rho (_r) and ln T (_T):
          T, p, e, s, nH      temperature, pressure, internal energy and entropy per volume,
                              and the hydrogen density
          ne_nH               electrons per hydrogen nucleus, with ne_nH_r and ne_nH_T
          chi_r, chi_T        d ln p / d ln rho and d ln p / d ln T
          E_r, E_T            de/d ln rho and de/d ln T
          G1                  first adiabatic exponent
        """
        n_h = self.n_h(rho)
        kT = K_B*T
        y, q1, q2, z = self.ionisation(rho, T)
        # a_j = d ln(Saha weight)/d ln T of each upper stage, and the energy stored in each
        a_h = 1.5 + CHI_H/kT
        a1 = 1.5 + CHI_HEI/kT
        a2 = a1 + 1.5 + CHI_HEII/kT
        e1 = CHI_HEI
        e2 = CHI_HEI + CHI_HEII
        # Helium moments: mean charge, mean a and mean stored energy per helium atom
        q_he = q1 + 2.0*q2
        a_he_mean = q1*a1 + q2*a2
        e_he = q1*e1 + q2*e2
        v_h = y*(1.0 - y)
        # Derivatives of the electrons per H (Z) and the ionisation energy per H (W) along u, ln T
        z_u = -(v_h + self.a_he*(q1 + 4.0*q2 - q_he**2))
        z_t = v_h*a_h + self.a_he*(q1*a1 + 2.0*q2*a2 - q_he*a_he_mean)
        w_u = -(CHI_H*v_h + self.a_he*(q1*e1 + 2.0*q2*e2 - q_he*e_he))
        w_t = CHI_H*v_h*a_h + self.a_he*(q1*a1*e1 + q2*a2*e2 - a_he_mean*e_he)
        # Response of u = ln ne to ln rho and ln T through the charge balance
        d = 1.0 - z_u/z
        u_r = 1.0/d
        u_t = z_t/(z*d)
        n_tot = 1.0 + self.a_he + z
        p = n_h*n_tot*kT
        e = n_h*(1.5*n_tot*kT + CHI_H*y + self.a_he*e_he)
        e_r = e + n_h*(1.5*kT*z*(u_r - 1.0) + w_u*u_r)
        e_t = n_h*(1.5*kT*(n_tot + z*u_t) + w_t + w_u*u_t)
        chi_r = 1.0 + z*(u_r - 1.0)/n_tot
        chi_t = 1.0 + z*u_t/n_tot
        return dict(T=T, p=p, e=e, s=self._entropy(n_h, T, y, q1, q2, z), ne_nH=z, nH=n_h,
                    chi_r=chi_r, chi_T=chi_t, E_r=e_r, E_T=e_t,
                    ne_nH_r=z*(u_r - 1.0), ne_nH_T=z*u_t, G1=chi_r + chi_t**2*p/e_t)

    def _entropy(self, n_h, T, y, q1, q2, z):
        """Entropy per volume, summed over every stage and the electrons"""
        n_he = self.a_he*n_h
        return (_translational_entropy((1.0 - y)*n_h, M_P, U_HI, T)
                + _translational_entropy(y*n_h, M_P, U_HII, T)
                + _translational_entropy((1.0 - q1 - q2)*n_he, 4*M_P, U_HEI, T)
                + _translational_entropy(q1*n_he, 4*M_P, U_HEII, T)
                + _translational_entropy(q2*n_he, 4*M_P, U_HEIII, T)
                + _translational_entropy(z*n_h, M_E, 2.0, T))

    def T_from_e(self, rho, e_target):
        """Temperature and state at which the internal energy per volume is e_target"""
        return self._root_lnT(rho, lambda th: (th['e'] - e_target, th['E_T']))

    def T_from_p(self, rho, p_target):
        """Temperature and state at which the pressure is p_target"""
        return self._root_lnT(rho, lambda th: (np.log(th['p']/p_target), th['chi_T']))

    def _root_lnT(self, rho, f, lo=np.log(300.0), hi=np.log(1e8)):
        """Root in ln T of f(state) -> (residual, d residual/d ln T), the residual rising with T

        Newton kept inside the bracket [lo, hi], which shrinks every step. A step that leaves the
        bracket, or that does not halve the step before last, becomes a bisection. A root at the
        bracket edge is an error rather than a clamp
        """
        lo0, hi0 = lo, hi
        lt = 0.5*(lo + hi)
        dx = dxold = hi - lo
        for _ in range(200):
            th = self.thermo(rho, np.exp(lt))
            r, dr = f(th)
            if r == 0.0:
                break
            if r > 0:
                hi = lt
            else:
                lo = lt
            nl = lt - r/dr
            if not (lo < nl < hi) or abs(nl - lt) > 0.5*abs(dxold):
                nl = 0.5*(lo + hi)
            dxold = dx
            dx = nl - lt
            lt = nl
            if abs(dx) <= 4e-16*abs(lt) or hi - lo <= 4e-16*abs(lt):
                break
        else:
            raise RuntimeError('SahaHHe: ln T root not converged at rho=%g' % rho)
        if not lo0 + 1e-12 < lt < hi0 - 1e-12:
            raise RuntimeError('SahaHHe: ln T root at the bracket edge, rho=%g' % rho)
        T = np.exp(lt)
        return T, self.thermo(rho, T)
