"""
Snow two-gamma Rankine-Hugoniot jump condition solver.

Solves the quadratic for the compression ratio r = rho2/rho1 across a shock
where the effective adiabatic index differs upstream (gamma1) and downstream
(gamma2).  Derived from the HD conservation laws in the shock frame with
piecewise-constant gamma.
"""

import numpy as np


class SnowJumpSolver:
    """Solve the two-gamma Rankine-Hugoniot quadratic (Eqns 10-12)."""

    def solve(self, M1, gamma1, gamma2):
        """Solve the quadratic for compression ratio r = rho2/rho1.

        Parameters
        ----------
        M1 : float
            Upstream Mach number (v1 / c1, where c1^2 = gamma1 * P1 / rho1).
        gamma1 : float
            Upstream adiabatic index.
        gamma2 : float
            Downstream adiabatic index.

        Returns
        -------
        r : float
            Physical compression ratio (r >= 1 for shocks).
        """
        M2 = M1**2

        # Quadratic coefficients (Eqns 10-12)
        A = (gamma1 / gamma2) * (gamma2 - 1) * (M2 / 2.0 + 1.0 / (gamma1 - 1))
        B = -(1.0 + gamma1 * M2)
        C = gamma1 * M2 - (gamma1 / gamma2) * (gamma2 - 1) * M2 / 2.0

        disc = B**2 - 4 * A * C
        if disc < 0:
            return np.nan

        r1 = (-B + np.sqrt(disc)) / (2 * A)
        r2 = (-B - np.sqrt(disc)) / (2 * A)

        # Physical root: r > 1 for shocks.  One root should be ~1 (trivial
        # when gamma1=gamma2) or < 1.  Take the larger root.
        r = max(r1, r2)
        return r

    def both_roots(self, M1, gamma1, gamma2):
        """Return both roots of the quadratic."""
        M2 = M1**2
        A = (gamma1 / gamma2) * (gamma2 - 1) * (M2 / 2.0 + 1.0 / (gamma1 - 1))
        B = -(1.0 + gamma1 * M2)
        C = gamma1 * M2 - (gamma1 / gamma2) * (gamma2 - 1) * M2 / 2.0
        disc = B**2 - 4 * A * C
        sq = np.sqrt(max(disc, 0))
        return (-B + sq) / (2 * A), (-B - sq) / (2 * A)

    def post_shock_state(self, rho1, v1, P1, gamma1, gamma2):
        """Compute full post-shock state from pre-shock primitives.

        Parameters
        ----------
        rho1, v1, P1 : float
            Pre-shock density, velocity (in shock frame), pressure.
        gamma1, gamma2 : float
            Upstream and downstream adiabatic indices.

        Returns
        -------
        dict with keys: rho2, v2, P2, r, M1
        """
        c1 = np.sqrt(gamma1 * P1 / rho1)
        M1 = abs(v1) / c1
        r = self.solve(M1, gamma1, gamma2)

        rho2 = rho1 * r
        v2 = v1 / r  # mass conservation

        # P2 from Eqn 4:
        # P2 = (g2-1)/g2 * rho2 * [0.5*(v1^2 - v2^2) + g1/(g1-1) * P1/rho1]
        P2 = ((gamma2 - 1) / gamma2) * rho2 * (
            0.5 * (v1**2 - v2**2) + gamma1 / (gamma1 - 1) * P1 / rho1
        )

        return dict(rho2=rho2, v2=v2, P2=P2, r=r, M1=M1)

    def strong_shock_limit(self, gamma2):
        """Maximum compression ratio in the strong shock limit (M1 -> inf).

        r_max = (gamma2 + 1) / (gamma2 - 1)
        """
        return (gamma2 + 1) / (gamma2 - 1)

    def standard_rh(self, M1, gamma):
        """Standard single-gamma Rankine-Hugoniot compression ratio.

        r = (gamma+1)*M1^2 / ((gamma-1)*M1^2 + 2)
        """
        M2 = M1**2
        return (gamma + 1) * M2 / ((gamma - 1) * M2 + 2)


# -- Verification ------------------------------------------------------

def verify_trivial():
    """Verify r=1 is a root when gamma1 = gamma2."""
    solver = SnowJumpSolver()
    for gamma in [1.1, 1.4, 5.0/3.0]:
        for M1 in [1.5, 3.0, 10.0]:
            r1, r2 = solver.both_roots(M1, gamma, gamma)
            # One root should be 1.0
            if not (abs(r1 - 1.0) < 1e-12 or abs(r2 - 1.0) < 1e-12):
                return False, f"gamma={gamma}, M1={M1}: roots = {r1}, {r2}"
    return True, "r=1 is always a root when gamma1=gamma2"


def verify_standard_rh():
    """Verify non-trivial root matches standard RH when gamma1 = gamma2."""
    solver = SnowJumpSolver()
    max_err = 0.0
    for gamma in [1.1, 1.4, 5.0/3.0]:
        for M1 in [1.01, 1.5, 2.0, 5.0, 10.0, 50.0]:
            r_snow = solver.solve(M1, gamma, gamma)
            r_std = solver.standard_rh(M1, gamma)
            err = abs(r_snow - r_std) / r_std
            max_err = max(max_err, err)
            if err > 1e-10:
                return False, f"gamma={gamma}, M1={M1}: r_snow={r_snow}, r_std={r_std}, err={err}"
    return True, f"max relative error = {max_err:.2e}"


def verify_strong_shock():
    """Verify r -> (gamma2+1)/(gamma2-1) as M1 -> infinity."""
    solver = SnowJumpSolver()
    max_err = 0.0
    for gamma2 in [1.06, 1.2, 1.4, 5.0/3.0]:
        gamma1 = 5.0 / 3.0  # upstream always fully neutral/ionised
        r_inf = solver.solve(1e6, gamma1, gamma2)
        r_lim = solver.strong_shock_limit(gamma2)
        err = abs(r_inf - r_lim) / r_lim
        max_err = max(max_err, err)
        if err > 1e-4:
            return False, f"gamma2={gamma2}: r(M=1e6)={r_inf}, r_lim={r_lim}, err={err}"
    return True, f"max relative error = {max_err:.2e}"


def run_verification():
    """Run all verification checks and print results."""
    checks = [
        ("Trivial (r=1 when g1=g2)", verify_trivial),
        ("Standard RH recovery",     verify_standard_rh),
        ("Strong shock limit",       verify_strong_shock),
    ]
    all_pass = True
    for name, fn in checks:
        passed, msg = fn()
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {name}: {msg}")
        if not passed:
            all_pass = False
    return all_pass


if __name__ == "__main__":
    print("Snow two-gamma jump condition verification:")
    ok = run_verification()
    print()
    if ok:
        print("All checks passed.")
    else:
        print("SOME CHECKS FAILED.")

    # Print a few example solutions
    solver = SnowJumpSolver()
    print("\nExample: gamma1=5/3, gamma2=1.2, M1=3.0")
    r = solver.solve(3.0, 5.0/3.0, 1.2)
    r_std = solver.standard_rh(3.0, 5.0/3.0)
    print(f"  r (two-gamma) = {r:.6f}")
    print(f"  r (standard, gamma=5/3) = {r_std:.6f}")
    print(f"  r_max(gamma2=1.2) = {solver.strong_shock_limit(1.2):.2f}")

    print("\nExample: gamma1=5/3, gamma2=1.06, M1=5.0")
    r = solver.solve(5.0, 5.0/3.0, 1.06)
    r_std = solver.standard_rh(5.0, 5.0/3.0)
    print(f"  r (two-gamma) = {r:.6f}")
    print(f"  r (standard, gamma=5/3) = {r_std:.6f}")
    print(f"  r_max(gamma2=1.06) = {solver.strong_shock_limit(1.06):.2f}")
