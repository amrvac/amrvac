"""Coleman (2020) Riemann test initial conditions in code units.

Single source of truth for the regression suite. Each entry contains
the (rho_L, v_L, T_L) and (rho_R, v_R, T_R) primitive states plus the
final time of the comparison snapshot.

These values match the par files in riemann_tests/coleman_test{N}.par
and the TESTS dicts in riemann_tests/snow_gamma_check/snow_*.py.
"""

TESTS = {
    1: dict(
        rho_L=1.5136914e1, rho_R=1.8921142e0,
        v_L=0.0, v_R=0.0,
        T_L=2.3670611e-2, T_R=9.7838523e-3,
        t_final=0.62933,
        shock_sides=("R",),
        description="Sod-like through He ionisation onset",
    ),
    2: dict(
        rho_L=6.0547655e2, rho_R=6.0547655e0,
        v_L=0.0, v_R=0.0,
        T_L=1.8936488e-2, T_R=2.9982773e-3,
        t_final=0.75520,
        shock_sides=("R",),
        description="Strong shock through H ionisation zone",
    ),
    3: dict(
        rho_L=1.2109531e2, rho_R=6.0547655e1,
        v_L=4.3697016e-1, v_R=-6.7531753e-1,
        T_L=9.4682442e-4, T_R=9.4682442e-4,
        # Coleman (2020) uses t_final = 3.776 for this test, but at that
        # time BOTH outgoing shocks (s_L=-0.18, s_R=+0.20) have left the
        # [0,1] domain (left at t=2.78 and 2.48 respectively). The
        # analytical solution is then a constant star state, and any
        # comparison is dominated by continuous-BC pollution injecting
        # the original L/R IC back in. Capping at 1.5 keeps both shocks
        # inside the domain (x_L=0.23, x_R=0.80) so the test exercises
        # the actual wave structure.
        t_final=1.50000,
        shock_sides=("L", "R"),
        description="Asymmetric collision (double shock)",
    ),
    4: dict(
        rho_L=7.5684569e1, rho_R=6.0547655e1,
        v_L=5.9586841e-1, v_R=-7.1504209e-1,
        T_L=9.4682442e-4, T_R=9.4682442e-4,
        # Same rationale as test 3: shocks at s_L=-0.24, s_R=+0.22 leave
        # at t=2.06 and 2.24; t_final=1.0 keeps them at x_L=0.26, x_R=0.72.
        t_final=1.00000,
        shock_sides=("L", "R"),
        description="Asymmetric collision (double shock)",
    ),
    5: dict(
        rho_L=1.2109531e4, rho_R=1.2109531e4,
        v_L=-3.1779648e-1, v_R=3.1779648e-1,
        T_L=1.4991387e-2, T_R=1.4991387e-2,
        t_final=0.62933,
        # Empty -- this is a double rarefaction; no shocks form, so the
        # Snow two-gamma jump check is not meaningful. The plotter skips
        # markers on tests with empty shock_sides.
        shock_sides=(),
        description="Symmetric outflow (double rarefaction)",
    ),
    6: dict(
        rho_L=9.0821483e3, rho_R=1.2109531e4,
        v_L=-1.9862280e-1, v_R=3.5752104e-1,
        T_L=1.4991387e-2, T_R=1.4991387e-2,
        t_final=0.62933,
        shock_sides=(),
        description="Asymmetric outflow (double rarefaction)",
    ),
}

# Variant suffix -> par-file naming convention used by riemann_tests/
VARIANTS = {
    "default": "",
    "geff":    "_geff",
    "hires":   "_hires",
    "pureH":   "_pureH",
}
