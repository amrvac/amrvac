#!/usr/bin/env python3
"""
Exact Riemann solver for general (tabulated) equation of state.

Implements the algorithm from:
  - Toro (2009), Ch. 4: general framework
  - Zingale & Katz (2015), ApJS 216:31: stellar EoS implementation
  - Coleman (2020), ApJS 248:7: hydrogen ionization tests
  - Chen et al. (2019), JCoPh 388:490: general EoS Riemann solver

The solver uses AMRVAC's LTE EoS tables (T, Ne/nH, p2eint, Gamma1)
to compute pressure, internal energy, and sound speed as functions
of density and internal energy per particle.

All quantities are in AMRVAC code units unless otherwise noted.
"""

import os
import struct
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


# 
# EoS TABLE READER AND INTERPOLATOR
# 

class EosTable:
    """A single 2D EoS table with PCHIP interpolation.

    Supports both uniform (legacy) and curvature-adaptive grids: a binary
    file with the new optional trailer carries explicit per-axis node
    arrays which override the uniform-spacing assumption. The Fortran
    side reads the same trailer (see read_eos_from_file in mod_eos.t)."""

    def __init__(self, dim1, dim2, v1min, v1max, v2min, v2max, table,
                 v1_nodes=None, v2_nodes=None):
        self.dim1 = dim1      # axis-1 (nH) length
        self.dim2 = dim2      # axis-2 (eint/nH or p/nH or T) length
        self.v1min = v1min    # log10 min/max bounds (still meaningful for
        self.v1max = v1max    # adaptive tables -- boundary nodes preserved)
        self.v2min = v2min
        self.v2max = v2max
        self.table = table    # (dim1, dim2) array
        # Adaptive trailer (None => uniform grid)
        self.v1_nodes = v1_nodes
        self.v2_nodes = v2_nodes
        self.is_uniform = (v1_nodes is None) and (v2_nodes is None)

    @classmethod
    def from_file(cls, fname):
        """Read AMRVAC binary EoS table (stream access, no record markers).

        Detects the optional adaptive trailer at end-of-file: trailer flag
        (int32 = 1) followed by dim1 + dim2 doubles giving explicit node
        coordinates. Uniform tables produce no trailer (EOF after the
        table values) and yield is_uniform = True."""
        with open(fname, 'rb') as f:
            dim1, dim2 = struct.unpack('<ii', f.read(8))
            v1min, v1max, v2min, v2max = struct.unpack('<dddd', f.read(32))
            data = np.frombuffer(f.read(dim1 * dim2 * 8), dtype='<f8')
            # Fortran column-major to C row-major
            table = data.reshape((dim2, dim1)).T.copy()

            # Optional adaptive-grid trailer
            v1_nodes = None
            v2_nodes = None
            trailer_buf = f.read(4)
            if len(trailer_buf) == 4:
                trailer_flag = struct.unpack('<i', trailer_buf)[0]
                if trailer_flag == 1:
                    v1_nodes = np.frombuffer(f.read(dim1 * 8), dtype='<f8').copy()
                    v2_nodes = np.frombuffer(f.read(dim2 * 8), dtype='<f8').copy()
        return cls(dim1, dim2, v1min, v1max, v2min, v2max, table,
                   v1_nodes=v1_nodes, v2_nodes=v2_nodes)

    def _locate_cell(self, v, nodes, vmin, vmax, n):
        """Return (i, t) -- 0-based cell index and local fractional offset.
        For uniform: affine map. For non-uniform: binary search on nodes."""
        if nodes is None:
            # Uniform branch
            r = (v - vmin) / (vmax - vmin) * (n - 1)
            r = max(0.0, min(r, n - 1))
            i = int(np.floor(r))
            t = r - i
            return i, t
        # Non-uniform branch
        if v <= nodes[0]:
            return 0, 0.0
        if v >= nodes[n - 1]:
            return n - 2, 1.0
        # np.searchsorted returns smallest i such that nodes[i] >= v
        ii = int(np.searchsorted(nodes, v, side='left'))
        i = max(0, min(ii - 1, n - 2))
        t = (v - nodes[i]) / (nodes[i + 1] - nodes[i])
        t = max(0.0, min(t, 1.0))
        return i, t

    def interp(self, v1, v2):
        """Monotone bicubic (PCHIP) interpolation. Branches on adaptive vs
        uniform grid via the cached node arrays."""
        ix, tx = self._locate_cell(v1, self.v1_nodes, self.v1min, self.v1max,
                                    self.dim1)
        iy, ty = self._locate_cell(v2, self.v2_nodes, self.v2min, self.v2max,
                                    self.dim2)

        # Gather 4 rows for PCHIP in y, each interpolated in x via PCHIP
        vals = []
        for di in [-1, 0, 1, 2]:
            ri = np.clip(ix + di, 0, self.dim1 - 1)
            # PCHIP along y-axis (dim2) for this row
            row = self.table[ri, :]
            v = self._pchip_1d(row, iy, ty)
            vals.append(v)

        # PCHIP along x-axis on the 4 intermediate values
        return self._pchip_4pt(vals, tx)

    @staticmethod
    def _pchip_1d(data, i, t):
        """PCHIP interpolation at fractional position i+t in a 1D array."""
        n = len(data)
        i0 = max(0, min(i - 1, n - 1))
        i1 = max(0, min(i, n - 1))
        i2 = max(0, min(i + 1, n - 1))
        i3 = max(0, min(i + 2, n - 1))

        p0, p1, p2, p3 = data[i0], data[i1], data[i2], data[i3]
        return EosTable._pchip_4pt([p0, p1, p2, p3], t)

    @staticmethod
    def _pchip_4pt(p, t):
        """Monotone cubic Hermite interpolation on 4 points p[0..3] at parameter t in [0,1].
        Points are at positions -1, 0, 1, 2; interpolation between p[1] and p[2]."""
        # Secant slopes
        d0 = p[1] - p[0]
        d1 = p[2] - p[1]
        d2 = p[3] - p[2]

        # PCHIP derivatives at p[1] and p[2]
        m1 = EosTable._pchip_slope(d0, d1)
        m2 = EosTable._pchip_slope(d1, d2)

        # Hermite basis
        h00 = 2*t**3 - 3*t**2 + 1
        h10 = t**3 - 2*t**2 + t
        h01 = -2*t**3 + 3*t**2
        h11 = t**3 - t**2

        return h00*p[1] + h10*m1 + h01*p[2] + h11*m2

    @staticmethod
    def _pchip_slope(d0, d1):
        """Monotone PCHIP derivative from two secants."""
        if d0 * d1 <= 0:
            return 0.0
        m = 2.0 * d0 * d1 / (d0 + d1)  # harmonic mean
        # Hyman monotonicity clamp
        if abs(m) > 3 * abs(d1):
            m = 3 * d1 * np.sign(m) / abs(np.sign(m)) if np.sign(m) != 0 else 0.0
            m = 3 * d1
        return m


class LTEEoS:
    """
    LTE equation of state with ionization energy, using AMRVAC tables.

    All public methods work in AMRVAC code units.
    """

    def __init__(self, tables_dir, composition='HHe', He_abundance=0.1,
                 unit_temperature=1e6, unit_numberdensity=1e15, unit_length=1e9,
                 gamma1_approx=False):

        self.He_abundance = He_abundance
        self.nH2rhoFactor = 1.0 + 4.0 * He_abundance

        # CGS constants
        mp_cgs = 1.672621777e-24
        kB_cgs = 1.3806488e-16

        # Unit system (a=1, b=1 for LTE in AMRVAC)
        self.unit_temperature = unit_temperature
        self.unit_numberdensity = unit_numberdensity
        self.unit_density = mp_cgs * unit_numberdensity         # a=1
        self.unit_pressure = kB_cgs * unit_numberdensity * unit_temperature  # b=1
        self.unit_velocity = np.sqrt(self.unit_pressure / self.unit_density)
        self.unit_length = unit_length
        self.unit_time = unit_length / self.unit_velocity

        log_nH_shift = np.log10(unit_numberdensity)
        log_eint_nH_shift = np.log10(self.unit_pressure / unit_numberdensity)
        log_T_shift = np.log10(unit_temperature)

        # Helper: shift both bounds and (if present) adaptive node arrays
        # together. Adaptive tables carry explicit node coords on disk in
        # CGS log-space; they must be converted to code-unit log-space
        # alongside v{1,2}min/max -- otherwise lookups silently clamp to
        # the boundary node and return garbage.
        def _shift_axes(tab, dx1, dx2):
            tab.v1min -= dx1; tab.v1max -= dx1
            tab.v2min -= dx2; tab.v2max -= dx2
            if tab.v1_nodes is not None:
                tab.v1_nodes = tab.v1_nodes - dx1
            if tab.v2_nodes is not None:
                tab.v2_nodes = tab.v2_nodes - dx2

        # --- Load T table ---
        T_raw = EosTable.from_file(f'{tables_dir}/LTEeos_T_{composition}_IonE.bin')
        # Raw values are T_cgs; convert to log10(T_code)
        T_raw.table = np.log10(T_raw.table) - log_T_shift
        _shift_axes(T_raw, log_nH_shift, log_eint_nH_shift)
        self.T_tab = T_raw

        # --- Load neOnH table ---
        ne_raw = EosTable.from_file(f'{tables_dir}/LTEeos_neOnH_{composition}_IonE.bin')
        # Raw values are Ne/nH (dimensionless); store as log10
        ne_raw.table = np.log10(np.maximum(ne_raw.table, 1e-30))
        _shift_axes(ne_raw, log_nH_shift, log_eint_nH_shift)
        self.neOnH_tab = ne_raw

        # --- Load p2eint table ---
        p2e_raw = EosTable.from_file(f'{tables_dir}/LTEeos_p2eint_{composition}_IonE.bin')
        # Raw values are eint/p (dimensionless, linear); axes are (nH, p/nH) in CGS
        log_p_nH_shift = log_eint_nH_shift  # same: pressure/nH has same units as eint/nH
        _shift_axes(p2e_raw, log_nH_shift, log_p_nH_shift)
        self.p2eint_tab = p2e_raw

        # Gamma_1 table: prefer loading the autodiff-exact one written by
        # generate_lte_tables.py (LTEeos_gamma1_*_IonE.bin) over the
        # finite-difference rebuild. Falls back to FD only if the binary
        # is missing.
        self.gamma1_approx = gamma1_approx
        gamma1_path = f'{tables_dir}/LTEeos_gamma1_{composition}_IonE.bin'
        if (not gamma1_approx) and os.path.exists(gamma1_path):
            g1_raw = EosTable.from_file(gamma1_path)
            # Gamma_1 is dimensionless; only the (log_nH, log_eint_nH) axes
            # need code-unit shifting.
            _shift_axes(g1_raw, log_nH_shift, log_eint_nH_shift)
            self.gamma1_tab = g1_raw
            self._gamma1_source = 'autodiff_binary'
        else:
            self._build_gamma1_table()
            self._gamma1_source = 'finite_difference_rebuild'

    # -- Thermodynamic functions (code units) --

    def pressure(self, log_nH, log_eint_nH):
        """Pressure from (nH, eint/nH). Returns p in code units."""
        log_T = self.T_tab.interp(log_nH, log_eint_nH)      # log10(T_code)
        log_y = self.neOnH_tab.interp(log_nH, log_eint_nH)   # log10(Ne/nH)
        T = 10.0**log_T
        y = 10.0**log_y
        nH = 10.0**log_nH
        return nH * T * (1.0 + self.He_abundance + y)

    def eint_from_p(self, log_nH, p):
        """Internal energy density from (nH, p). Returns eint_vol in code units."""
        nH = 10.0**log_nH
        log_p_nH = np.log10(p / nH)
        eint_over_p = self.p2eint_tab.interp(log_nH, log_p_nH)  # linear, no antilog
        return p * eint_over_p  # eint_vol = p * (eint/p)

    def temperature(self, log_nH, log_eint_nH):
        """Temperature in code units."""
        return 10.0**self.T_tab.interp(log_nH, log_eint_nH)

    def gamma1(self, log_nH, log_eint_nH):
        """First adiabatic index Gamma_1."""
        return self.gamma1_tab.interp(log_nH, log_eint_nH)

    def sound_speed(self, log_nH, log_eint_nH):
        """Sound speed in code units."""
        g1 = self.gamma1(log_nH, log_eint_nH)
        p = self.pressure(log_nH, log_eint_nH)
        rho = 10.0**log_nH * self.nH2rhoFactor
        return np.sqrt(g1 * p / rho)

    def _build_gamma1_table(self):
        """Build Gamma_1 table from T and neOnH tables via finite differences.

        Gamma_1 = a + b * (p / eint_vol), where:
            a = d(log p) / d(log nH) at constant eint/nH
            b = d(log p) / d(log eint/nH) at constant nH
        """
        dim1 = self.T_tab.dim1
        dim2 = self.T_tab.dim2
        v1min = self.T_tab.v1min
        v1max = self.T_tab.v1max
        v2min = self.T_tab.v2min
        v2max = self.T_tab.v2max

        # Build log10(p) on the grid
        log_nH_arr = np.linspace(v1min, v1max, dim1)
        log_eint_nH_arr = np.linspace(v2min, v2max, dim2)
        d1 = (v1max - v1min) / (dim1 - 1)
        d2 = (v2max - v2min) / (dim2 - 1)

        log_p = np.zeros((dim1, dim2))
        for i in range(dim1):
            for j in range(dim2):
                p = self.pressure(log_nH_arr[i], log_eint_nH_arr[j])
                log_p[i, j] = np.log10(max(p, 1e-30))

        # Compute partial derivatives via central differences
        g1_table = np.zeros((dim1, dim2))
        for i in range(dim1):
            for j in range(dim2):
                # a = d(log p)/d(log nH)
                if i == 0:
                    a = (log_p[i+1, j] - log_p[i, j]) / d1
                elif i == dim1 - 1:
                    a = (log_p[i, j] - log_p[i-1, j]) / d1
                else:
                    a = (log_p[i+1, j] - log_p[i-1, j]) / (2 * d1)

                # b = d(log p)/d(log eint/nH)
                if j == 0:
                    b = (log_p[i, j+1] - log_p[i, j]) / d2
                elif j == dim2 - 1:
                    b = (log_p[i, j] - log_p[i, j-1]) / d2
                else:
                    b = (log_p[i, j+1] - log_p[i, j-1]) / (2 * d2)

                nH = 10.0**log_nH_arr[i]
                eint_nH = 10.0**(log_eint_nH_arr[j])
                eint_vol = nH * eint_nH
                p_val = 10.0**log_p[i, j]

                if self.gamma1_approx:
                    g1_val = 1.0 + p_val / eint_vol
                else:
                    g1_val = a + b * (p_val / eint_vol)
                g1_table[i, j] = np.clip(g1_val, 1.0, 5.0/3.0)

        self.gamma1_tab = EosTable(dim1, dim2, v1min, v1max, v2min, v2max, g1_table)


# 
# PURE HYDROGEN EoS (Coleman 2020, analytic Saha)
# 

class PureHydrogenEoS:
    """
    Analytic pure-hydrogen EoS using the Saha equation.

    Implements Coleman (2020, ApJS, 248, 7) equations exactly:
      - Saha ionization: x^2/(1-x) = (2/rho) T^{3/2} exp(-1/T)
      - Pressure:  p = nH T (1+x)
      - Energy:    e_int = nH [x T_ion + 1.5 T (1+x)]
      - Gamma_1:   from analytic c_s^2 (Coleman Appendix B)

    All public methods use AMRVAC code units, matching the LTEEoS interface.
    """

    def __init__(self, unit_temperature=1e6, unit_numberdensity=1e15, unit_length=1e9):
        mp_cgs = 1.672621777e-24
        kB_cgs = 1.3806488e-16
        me_cgs = 9.10938e-28
        h_cgs  = 6.62607e-27

        self.nH2rhoFactor = 1.0   # pure H, no He

        # Unit system (same as AMRVAC with eq_state_units=.false.)
        self.unit_temperature = unit_temperature
        self.unit_numberdensity = unit_numberdensity
        self.unit_density  = mp_cgs * unit_numberdensity
        self.unit_pressure = kB_cgs * unit_numberdensity * unit_temperature
        self.unit_velocity = np.sqrt(self.unit_pressure / self.unit_density)
        self.unit_length = unit_length
        self.unit_time = unit_length / self.unit_velocity

        # Coleman reference scales: E_ion = alpha^2 m_e c^2 / 2
        alpha = 7.2973525693e-3
        c_cgs = 2.99792458e10
        chi_H = 0.5 * alpha**2 * me_cgs * c_cgs**2                      # 2.1799e-11 erg
        self.T_ion = chi_H / kB_cgs                                      # 157888 K
        self.n_q = (2*np.pi * me_cgs * kB_cgs * self.T_ion
                    / h_cgs**2)**1.5                                      # 1.515e23 cm^-3

        # Dimensionless conversion factors
        self.T_ion_code = self.T_ion / unit_temperature      # T_ion in code units
        self.nq_code    = self.n_q / unit_numberdensity      # n_q in code nH units

    # -- Saha ionization fraction --

    def _ionization_fraction(self, nH_code, T_code):
        """Saha ionization fraction x for pure H (scalar or array)."""
        T_dim = np.asarray(T_code) / self.T_ion_code        # T / T_ion
        rho_dim = np.asarray(nH_code) / self.nq_code         # nH / n_q

        # K = (1/rho) T^{3/2} exp(-1/T)   [g_e g_ion/g_atom = 2x1/2 = 1]
        with np.errstate(over='ignore', divide='ignore'):
            exponent = np.where(T_dim > 1e-4, -1.0 / T_dim, -1e4)
            K = (1.0 / rho_dim) * T_dim**1.5 * np.exp(exponent)

        # x^2 + Kx - K = 0  ->  x = [-K + sqrt(K^2 + 4K)] / 2
        x = (-K + np.sqrt(K**2 + 4.0 * K)) / 2.0
        return np.clip(x, 0.0, 1.0)

    # -- Temperature  eint inversion --

    def _T_from_eint_nH(self, nH_code, eint_nH_code):
        """Invert e_int/nH -> T using Brent's method (scalar)."""
        target = float(eint_nH_code)
        nH = float(nH_code)

        def residual(T_trial):
            x = float(self._ionization_fraction(nH, T_trial))
            return x * self.T_ion_code + 1.5 * T_trial * (1.0 + x) - target

        # Bracket: at T->0, eint/nH->0; eint/nH increases with T
        T_lo = 1e-10
        T_hi = max(target / 1.5, 1.0)
        while residual(T_hi) < 0:
            T_hi *= 2.0

        return brentq(residual, T_lo, T_hi, xtol=1e-14, rtol=1e-14)

    def eint_nH_from_T(self, log_nH, T_code):
        """Direct: (nH, T) -> eint/nH in code units (no table inversion needed)."""
        nH = 10.0**log_nH
        x = self._ionization_fraction(nH, T_code)
        return float(x) * self.T_ion_code + 1.5 * float(T_code) * (1.0 + float(x))

    # -- Thermodynamic functions (code units, same interface as LTEEoS) --

    def pressure(self, log_nH, log_eint_nH):
        """p = nH T (1+x)."""
        nH = 10.0**log_nH
        eint_nH = 10.0**log_eint_nH
        T = self._T_from_eint_nH(nH, eint_nH)
        x = float(self._ionization_fraction(nH, T))
        return nH * T * (1.0 + x)

    def eint_from_p(self, log_nH, p):
        """eint_vol from (nH, p) by inverting p = nH T (1+x)."""
        nH = 10.0**log_nH

        def residual(T_trial):
            x = float(self._ionization_fraction(nH, T_trial))
            return nH * T_trial * (1.0 + x) - p

        T_lo = 1e-10
        T_hi = p / nH   # upper bound: 1+x >= 1
        while residual(T_hi) < 0:
            T_hi *= 2.0

        T = brentq(residual, T_lo, T_hi, xtol=1e-14, rtol=1e-14)
        x = float(self._ionization_fraction(nH, T))
        return nH * (x * self.T_ion_code + 1.5 * T * (1.0 + x))

    def temperature(self, log_nH, log_eint_nH):
        """T in code units."""
        nH = 10.0**log_nH
        eint_nH = 10.0**log_eint_nH
        return self._T_from_eint_nH(nH, eint_nH)

    def gamma1(self, log_nH, log_eint_nH):
        """Gamma_1 from Coleman (2020) analytic formula."""
        nH = 10.0**log_nH
        eint_nH = 10.0**log_eint_nH
        T = self._T_from_eint_nH(nH, eint_nH)
        x = float(self._ionization_fraction(nH, T))
        return self._gamma1_analytic(T, x)

    def _gamma1_analytic(self, T_code, x):
        """
        Coleman (2020) analytic Gamma_1 for pure hydrogen.

        c_s^2 / (p/rho) = Gamma_1 = (1+x)T [ 1 + (2c + 2T(4+3T)(1-x)x)
                                          / (3c + (2+3T)^2(1-x)x) ]  /  [(1+x)T]
        where c = 2 T^2 (2 + x - x^2),  T = T/T_ion.
        """
        T_dim = T_code / self.T_ion_code
        c = 2.0 * T_dim**2 * (2.0 + x - x**2)
        denom = 3.0 * c + (2.0 + 3.0 * T_dim)**2 * (1.0 - x) * x
        if denom < 1e-30:
            return 5.0 / 3.0
        numer = 2.0 * c + 2.0 * T_dim * (4.0 + 3.0 * T_dim) * (1.0 - x) * x
        a_sq = (1.0 + x) * T_dim * (1.0 + numer / denom)
        p_over_rho = T_dim * (1.0 + x)   # p/(nH m_H) in dimensionless
        return a_sq / p_over_rho

    def sound_speed(self, log_nH, log_eint_nH):
        """Sound speed in code units."""
        g1 = self.gamma1(log_nH, log_eint_nH)
        p  = self.pressure(log_nH, log_eint_nH)
        rho = 10.0**log_nH * self.nH2rhoFactor
        return np.sqrt(g1 * p / rho)


# 
# EXACT RIEMANN SOLVER FOR GENERAL EoS
# 

class GeneralEoSRiemann:
    """
    Exact Riemann solver for a general (tabulated) equation of state.

    Algorithm:
    1. For each trial p*, compute the wave function on each side:
       - Shock (p* > p_k): solve Hugoniot for rho*, then f = velocity jump
       - Rarefaction (p* < p_k): integrate isentrope, f = velocity integral
    2. Newton-Raphson/Brent iteration on p* such that:
       f_L(p*) + f_R(p*) + (u_R - u_L) = 0
    3. Sample the solution at each x for a given time t.
    """

    def __init__(self, eos, tol=1e-10, max_iter=100):
        self.eos = eos
        self.tol = tol
        self.max_iter = max_iter

    def solve(self, rho_L, u_L, p_L, eint_L,
                    rho_R, u_R, p_R, eint_R):
        """
        Solve the Riemann problem.

        Parameters: primitive state on each side (code units).
        eint_L/R is internal energy DENSITY (not specific).

        Returns dict with star-region quantities and wave structure.
        """
        F = self.eos.nH2rhoFactor

        # Precompute EoS quantities for L and R states
        nH_L = rho_L / F
        nH_R = rho_R / F
        log_nH_L = np.log10(nH_L)
        log_nH_R = np.log10(nH_R)
        log_eint_nH_L = np.log10(eint_L / nH_L)
        log_eint_nH_R = np.log10(eint_R / nH_R)

        c_L = self.eos.sound_speed(log_nH_L, log_eint_nH_L)
        c_R = self.eos.sound_speed(log_nH_R, log_eint_nH_R)

        # Initial guess: two-rarefaction approximation (linearized)
        p_guess = 0.5 * (p_L + p_R) - 0.125 * (u_R - u_L) * (rho_L + rho_R) * (c_L + c_R)
        p_guess = max(p_guess, 1e-10 * min(p_L, p_R))

        # Store states for wave functions
        self._state_L = (rho_L, u_L, p_L, eint_L, nH_L, log_nH_L, log_eint_nH_L, c_L)
        self._state_R = (rho_R, u_R, p_R, eint_R, nH_R, log_nH_R, log_eint_nH_R, c_R)

        du = u_R - u_L

        def residual(p_try):
            fL = self._wave_function_L(p_try)
            fR = self._wave_function_R(p_try)
            return fL + fR + du

        # Bracket the root: residual < 0 at low p (rarefaction dominates),
        # > 0 at high p (shock dominates). Use geometric search.
        p_lo = p_R * 0.5
        p_hi = p_L * 1.5
        r_lo = residual(p_lo)
        r_hi = residual(p_hi)

        # Expand bracket if needed
        for _ in range(30):
            if r_lo * r_hi < 0:
                break
            if r_lo > 0:
                p_lo *= 0.1
                r_lo = residual(p_lo)
            if r_hi < 0:
                p_hi *= 3.0
                r_hi = residual(p_hi)

        # Use Brent's method (robust, guaranteed convergence with bracket)
        try:
            p_star = brentq(residual, p_lo, p_hi, xtol=1e-12, rtol=self.tol)
        except ValueError:
            # Fallback: Newton from the initial guess
            p_star = p_guess
            for it in range(self.max_iter):
                res = residual(p_star)
                dfL = self._wave_function_deriv_L(p_star)
                dfR = self._wave_function_deriv_R(p_star)
                dres = dfL + dfR
                if abs(dres) < 1e-30:
                    break
                dp = -res / dres
                p_new = max(p_star + dp, 1e-10 * min(p_L, p_R))
                if abs(dp) / (0.5 * (p_star + p_new) + 1e-30) < self.tol:
                    p_star = p_new
                    break
                p_star = p_new

        # Velocity in star region: u* = u_L - f_L (Toro convention)
        fL = self._wave_function_L(p_star)
        u_star = u_L - fL

        # Density in star region (both sides of contact)
        rho_starL = self._star_density_L(p_star)
        rho_starR = self._star_density_R(p_star)
        eint_starL = self.eos.eint_from_p(np.log10(rho_starL / F), p_star)
        eint_starR = self.eos.eint_from_p(np.log10(rho_starR / F), p_star)

        return {
            'p_star': p_star, 'u_star': u_star,
            'rho_starL': rho_starL, 'rho_starR': rho_starR,
            'eint_starL': eint_starL, 'eint_starR': eint_starR,
            'c_L': c_L, 'c_R': c_R,
        }

    def sample(self, sol, x, t, x0=0.5,
               rho_L=None, u_L=None, p_L=None, eint_L=None,
               rho_R=None, u_R=None, p_R=None, eint_R=None):
        """Sample the exact solution at positions x and time t."""
        F = self.eos.nH2rhoFactor
        p_star = sol['p_star']
        u_star = sol['u_star']
        rho_sL = sol['rho_starL']
        rho_sR = sol['rho_starR']
        eint_sL = sol['eint_starL']
        eint_sR = sol['eint_starR']

        nH_sL = rho_sL / F
        nH_sR = rho_sR / F
        c_sL = self.eos.sound_speed(np.log10(nH_sL), np.log10(eint_sL / nH_sL))
        c_sR = self.eos.sound_speed(np.log10(nH_sR), np.log10(eint_sR / nH_sR))

        c_L = sol['c_L']
        c_R = sol['c_R']

        rho = np.empty_like(x, dtype=float)
        u = np.empty_like(x, dtype=float)
        p = np.empty_like(x, dtype=float)
        eint = np.empty_like(x, dtype=float)

        for i, xi in enumerate(x):
            S = (xi - x0) / t if t > 0 else 0.0

            if S < u_star:
                # Left of contact
                if p_star > p_L:
                    # Left shock
                    nH_L_val = rho_L / F
                    j_L = np.sqrt((p_star - p_L) / (1.0/rho_L - 1.0/rho_sL + 1e-30))
                    s_L = u_L - j_L / rho_L
                    if S < s_L:
                        rho[i] = rho_L; u[i] = u_L; p[i] = p_L; eint[i] = eint_L
                    else:
                        rho[i] = rho_sL; u[i] = u_star; p[i] = p_star; eint[i] = eint_sL
                else:
                    # Left rarefaction
                    s_HL = u_L - c_L        # head speed
                    s_TL = u_star - c_sL    # tail speed
                    if S < s_HL:
                        rho[i] = rho_L; u[i] = u_L; p[i] = p_L; eint[i] = eint_L
                    elif S > s_TL:
                        rho[i] = rho_sL; u[i] = u_star; p[i] = p_star; eint[i] = eint_sL
                    else:
                        # Inside the rarefaction fan
                        rho_f, u_f, p_f, eint_f = self._rarefaction_fan_L(
                            S, rho_L, u_L, p_L, eint_L, rho_sL, u_star, p_star, eint_sL)
                        rho[i] = rho_f; u[i] = u_f; p[i] = p_f; eint[i] = eint_f
            else:
                # Right of contact
                if p_star > p_R:
                    # Right shock
                    j_R = np.sqrt((p_star - p_R) / (1.0/rho_R - 1.0/rho_sR + 1e-30))
                    s_R = u_R + j_R / rho_R
                    if S > s_R:
                        rho[i] = rho_R; u[i] = u_R; p[i] = p_R; eint[i] = eint_R
                    else:
                        rho[i] = rho_sR; u[i] = u_star; p[i] = p_star; eint[i] = eint_sR
                else:
                    # Right rarefaction
                    s_HR = u_R + c_R
                    s_TR = u_star + c_sR
                    if S > s_HR:
                        rho[i] = rho_R; u[i] = u_R; p[i] = p_R; eint[i] = eint_R
                    elif S < s_TR:
                        rho[i] = rho_sR; u[i] = u_star; p[i] = p_star; eint[i] = eint_sR
                    else:
                        rho_f, u_f, p_f, eint_f = self._rarefaction_fan_R(
                            S, rho_R, u_R, p_R, eint_R, rho_sR, u_star, p_star, eint_sR)
                        rho[i] = rho_f; u[i] = u_f; p[i] = p_f; eint[i] = eint_f

        return rho, u, p, eint

    # -- Wave functions --

    def _wave_function_L(self, p_star):
        """Left wave function f_L (Toro convention).

        f_L + f_R + (u_R - u_L) = 0,  u* = u_L - f_L.

        Left shock (p* > p_L):  f_L = (p* - p_L)/j > 0
        Left rarefaction (p* < p_L): f_L = u_L - u* < 0
        """
        rho_L, u_L, p_L, eint_L, nH_L, log_nH_L, log_eint_nH_L, c_L = self._state_L
        if p_star >= p_L:
            # Shock: f_L = +(p* - p_L)/j  (sign=+1, not -1)
            return self._shock_f(p_star, rho_L, p_L, eint_L, sign=+1)
        else:
            # Rarefaction integral returns u* - u_L; negate for f_L = u_L - u*
            return -self._rarefaction_f_L(p_star)

    def _wave_function_R(self, p_star):
        """Right wave function f_R (Toro convention).

        f_L + f_R + (u_R - u_L) = 0,  u* = u_R + f_R.

        Right shock (p* > p_R):  f_R = (p* - p_R)/j > 0
        Right rarefaction (p* < p_R): f_R = u* - u_R < 0
        """
        rho_R, u_R, p_R, eint_R, nH_R, log_nH_R, log_eint_nH_R, c_R = self._state_R
        if p_star >= p_R:
            return self._shock_f(p_star, rho_R, p_R, eint_R, sign=+1)
        else:
            return self._rarefaction_f_R(p_star)

    def _wave_function_deriv_L(self, p_star):
        """Numerical derivative of left wave function."""
        dp = max(abs(p_star) * 1e-6, 1e-12)
        return (self._wave_function_L(p_star + dp) - self._wave_function_L(p_star - dp)) / (2 * dp)

    def _wave_function_deriv_R(self, p_star):
        """Numerical derivative of right wave function."""
        dp = max(abs(p_star) * 1e-6, 1e-12)
        return (self._wave_function_R(p_star + dp) - self._wave_function_R(p_star - dp)) / (2 * dp)

    def _shock_f(self, p_star, rho_0, p_0, eint_0, sign):
        """Shock wave function: velocity jump magnitude across a shock.

        Returns sign * (p* - p_0) / j where j is the mass flux.
        With sign=+1, this gives (p*-p_0)/j > 0, matching Toro's convention
        for both left and right shocks.
        """
        F = self.eos.nH2rhoFactor

        # Solve Hugoniot for rho_star using SPECIFIC internal energy:
        #   epsilon(rho*, p*) = epsilon_0 + 0.5*(p* + p_0)*(1/rho_0 - 1/rho*)
        # where epsilon = eint_vol / rho (Toro eq 4.8, Zeldovich & Raizer).
        epsilon_0 = eint_0 / rho_0

        def hugoniot_residual(rho_star):
            nH_star = rho_star / F
            eint_star_vol = self.eos.eint_from_p(np.log10(nH_star), p_star)
            epsilon_star = eint_star_vol / rho_star
            rhs = epsilon_0 + 0.5 * (p_star + p_0) * (1.0/rho_0 - 1.0/rho_star)
            return epsilon_star - rhs

        # For a compressive shock: rho* > rho_0
        rho_lo = rho_0 * 1.001
        rho_hi = rho_0 * 10.0

        # Expand bracket if needed
        for _ in range(20):
            try:
                f_lo = hugoniot_residual(rho_lo)
                f_hi = hugoniot_residual(rho_hi)
                if f_lo * f_hi < 0:
                    break
            except:
                pass
            rho_hi *= 2.0
        else:
            # Fallback: use ideal gas estimate
            gamma = 5.0/3.0
            gp1 = gamma + 1; gm1 = gamma - 1
            rho_star = rho_0 * ((p_star/p_0 + gm1/gp1) / (gm1/gp1 * p_star/p_0 + 1.0))
            j = np.sqrt(rho_0 * (0.5*gp1*p_star + 0.5*gm1*p_0))
            return sign * (p_star - p_0) / j

        try:
            rho_star = brentq(hugoniot_residual, rho_lo, rho_hi, xtol=1e-12, rtol=1e-12)
        except:
            rho_star = 0.5 * (rho_lo + rho_hi)

        # Mass flux and velocity jump
        dv = 1.0/rho_0 - 1.0/rho_star
        if abs(dv) < 1e-30:
            return 0.0
        j_sq = (p_star - p_0) / dv
        j = np.sqrt(abs(j_sq))
        return sign * (p_star - p_0) / j

    def _rarefaction_f_L(self, p_star):
        """Left rarefaction wave function via isentrope integration.

        Integrates de/drho = p/rho along the isentrope from state_L down to p*.
        Returns u* - u_L = -integral(c/rho, drho).
        """
        rho_L, u_L, p_L, eint_L, nH_L, log_nH_L, log_eint_nH_L, c_L = self._state_L
        return self._rarefaction_integral(rho_L, p_L, eint_L, p_star, sign=-1)

    def _rarefaction_f_R(self, p_star):
        """Right rarefaction: integrates from state_R to p*."""
        rho_R, u_R, p_R, eint_R, nH_R, log_nH_R, log_eint_nH_R, c_R = self._state_R
        return self._rarefaction_integral(rho_R, p_R, eint_R, p_star, sign=+1)

    def _rarefaction_integral(self, rho_0, p_0, eint_0, p_target, sign):
        """Integrate along the isentrope from (rho_0, eint_0) until p = p_target.

        The isentrope ODE (ds=0, first law Tds = de + pd(1/rho) = 0):
            d(epsilon)/d(rho) = p/rho^2      (specific internal energy)
            d(eint_vol)/d(rho) = (eint_vol + p)/rho   (volumetric, via product rule)

        Uses scipy.integrate.solve_ivp with adaptive RK45 for high accuracy.

        sign=-1 for left (du = -c/rho * drho), +1 for right (du = +c/rho * drho).

        Returns the velocity change u* - u_0.
        """
        F = self.eos.nH2rhoFactor

        # State: y = [eint_vol, u_change]
        def rhs(rho, y):
            eint_vol = y[0]
            nH = rho / F
            if nH <= 0 or eint_vol <= 0:
                return [0.0, 0.0]
            log_nH = np.log10(nH)
            eint_nH = eint_vol / nH
            if eint_nH <= 0:
                return [0.0, 0.0]
            log_eint_nH = np.log10(eint_nH)
            p_val = self.eos.pressure(log_nH, log_eint_nH)
            c_val = self.eos.sound_speed(log_nH, log_eint_nH)
            return [(eint_vol + p_val) / rho, sign * c_val / rho]

        # Event: pressure crosses target (terminal event to stop integration)
        def pressure_event(rho, y):
            eint_vol = y[0]
            nH = rho / F
            if nH <= 0 or eint_vol <= 0:
                return -1.0
            log_nH = np.log10(nH)
            eint_nH = eint_vol / nH
            if eint_nH <= 0:
                return -1.0
            p_val = self.eos.pressure(log_nH, np.log10(eint_nH))
            return p_val - p_target
        pressure_event.terminal = True
        pressure_event.direction = -1  # p decreasing through target

        # Integration range: rho_0 down to a conservative lower bound
        gamma_eff = 5.0 / 3.0
        rho_est = rho_0 * (p_target / p_0) ** (1.0 / gamma_eff)
        rho_end = max(rho_est * 0.1, rho_0 * 1e-4)

        sol = solve_ivp(rhs, [rho_0, rho_end], [eint_0, 0.0],
                        method='RK45', events=pressure_event,
                        rtol=1e-10, atol=1e-12, max_step=abs(rho_0 - rho_end) / 100)

        if sol.t_events[0].size > 0:
            # Event triggered: return u_change at the event
            return sol.y_events[0][0][1]
        else:
            # Event not triggered: return last value
            return sol.y[1, -1]

    def _star_density_L(self, p_star):
        """Find rho* on the left side of the contact."""
        rho_L, u_L, p_L, eint_L, nH_L, log_nH_L, log_eint_nH_L, c_L = self._state_L
        if p_star >= p_L:
            return self._hugoniot_density(rho_L, p_L, eint_L, p_star)
        else:
            return self._isentrope_density(rho_L, p_L, eint_L, p_star)

    def _star_density_R(self, p_star):
        """Find rho* on the right side of the contact."""
        rho_R, u_R, p_R, eint_R, nH_R, log_nH_R, log_eint_nH_R, c_R = self._state_R
        if p_star >= p_R:
            return self._hugoniot_density(rho_R, p_R, eint_R, p_star)
        else:
            return self._isentrope_density(rho_R, p_R, eint_R, p_star)

    def _hugoniot_density(self, rho_0, p_0, eint_0, p_star):
        """Solve Hugoniot for post-shock density using specific internal energy."""
        F = self.eos.nH2rhoFactor
        epsilon_0 = eint_0 / rho_0

        def residual(rho_s):
            nH_s = rho_s / F
            eint_s_vol = self.eos.eint_from_p(np.log10(nH_s), p_star)
            epsilon_s = eint_s_vol / rho_s
            rhs = epsilon_0 + 0.5 * (p_star + p_0) * (1.0/rho_0 - 1.0/rho_s)
            return epsilon_s - rhs

        rho_lo = rho_0 * 1.001
        rho_hi = rho_0 * 10.0
        for _ in range(20):
            try:
                if residual(rho_lo) * residual(rho_hi) < 0:
                    break
            except:
                pass
            rho_hi *= 2.0

        try:
            return brentq(residual, rho_lo, rho_hi, xtol=1e-12, rtol=1e-12)
        except:
            # Ideal gas fallback
            gamma = 5.0/3.0
            gp1 = gamma + 1; gm1 = gamma - 1
            return rho_0 * ((p_star/p_0 + gm1/gp1) / (gm1/gp1 * p_star/p_0 + 1.0))

    def _isentrope_density(self, rho_0, p_0, eint_0, p_target):
        """Follow isentrope from (rho_0, eint_0) to p_target, return rho."""
        F = self.eos.nH2rhoFactor

        def rhs(rho, y):
            eint_vol = y[0]
            nH = rho / F
            if nH <= 0 or eint_vol <= 0:
                return [0.0]
            p_val = self.eos.pressure(np.log10(nH), np.log10(eint_vol / nH))
            return [(eint_vol + p_val) / rho]

        def pressure_event(rho, y):
            eint_vol = y[0]
            nH = rho / F
            if nH <= 0 or eint_vol <= 0:
                return -1.0
            p_val = self.eos.pressure(np.log10(nH), np.log10(eint_vol / nH))
            return p_val - p_target
        pressure_event.terminal = True
        pressure_event.direction = -1

        gamma_eff = 5.0 / 3.0
        rho_est = rho_0 * (p_target / p_0) ** (1.0 / gamma_eff)
        rho_end = max(rho_est * 0.1, rho_0 * 1e-4)

        sol = solve_ivp(rhs, [rho_0, rho_end], [eint_0],
                        method='RK45', events=pressure_event,
                        rtol=1e-10, atol=1e-12, max_step=abs(rho_0 - rho_end) / 100)

        if sol.t_events[0].size > 0:
            return sol.t_events[0][0]  # rho at the event
        else:
            return sol.t[-1]

    def _rarefaction_fan_L(self, S, rho_L, u_L, p_L, eint_L,
                            rho_sL, u_star, p_star, eint_sL):
        """Sample inside a left rarefaction fan at speed S = (x-x0)/t.

        For a left-going 1-wave: u - c = S defines the fan state.
        Integrate the isentrope from rho_L to rho_sL, tracking u, and
        find rho where u(rho) - c(rho) = S.
        """
        F = self.eos.nH2rhoFactor

        def rhs(rho, y):
            eint_vol, u_val = y
            nH = rho / F
            if nH <= 0 or eint_vol <= 0:
                return [0.0, 0.0]
            log_nH = np.log10(nH)
            eint_nH = eint_vol / nH
            if eint_nH <= 0:
                return [0.0, 0.0]
            p_val = self.eos.pressure(log_nH, np.log10(eint_nH))
            c_val = self.eos.sound_speed(log_nH, np.log10(eint_nH))
            return [(eint_vol + p_val) / rho, -c_val / rho]

        rho_end = max(rho_sL, rho_L * 0.01)
        sol = solve_ivp(rhs, [rho_L, rho_end], [eint_L, u_L],
                        method='RK45', rtol=1e-8, atol=1e-10, dense_output=True,
                        max_step=abs(rho_L - rho_end) / 50)

        # Find rho where u - c = S using the dense output
        best_err = float('inf')
        best_state = (rho_L, u_L, p_L, eint_L)

        for rho_k in np.linspace(rho_L, rho_end, 500):
            y_k = sol.sol(rho_k)
            eint_vol_k, u_k = y_k
            nH = rho_k / F
            if nH <= 0 or eint_vol_k <= 0:
                continue
            log_nH = np.log10(nH)
            eint_nH = eint_vol_k / nH
            if eint_nH <= 0:
                continue
            c_k = self.eos.sound_speed(log_nH, np.log10(eint_nH))
            p_k = self.eos.pressure(log_nH, np.log10(eint_nH))
            err = abs(u_k - c_k - S)
            if err < best_err:
                best_err = err
                best_state = (rho_k, u_k, p_k, eint_vol_k)

        return best_state

    def _rarefaction_fan_R(self, S, rho_R, u_R, p_R, eint_R,
                            rho_sR, u_star, p_star, eint_sR):
        """Sample inside a right rarefaction fan at speed S = (x-x0)/t.

        For a right-going 3-wave: u + c = S defines the fan state.
        """
        F = self.eos.nH2rhoFactor

        def rhs(rho, y):
            eint_vol, u_val = y
            nH = rho / F
            if nH <= 0 or eint_vol <= 0:
                return [0.0, 0.0]
            log_nH = np.log10(nH)
            eint_nH = eint_vol / nH
            if eint_nH <= 0:
                return [0.0, 0.0]
            p_val = self.eos.pressure(log_nH, np.log10(eint_nH))
            c_val = self.eos.sound_speed(log_nH, np.log10(eint_nH))
            return [(eint_vol + p_val) / rho, c_val / rho]

        rho_end = max(rho_sR, rho_R * 0.01)
        sol = solve_ivp(rhs, [rho_R, rho_end], [eint_R, u_R],
                        method='RK45', rtol=1e-8, atol=1e-10, dense_output=True,
                        max_step=abs(rho_R - rho_end) / 50)

        best_err = float('inf')
        best_state = (rho_R, u_R, p_R, eint_R)

        for rho_k in np.linspace(rho_R, rho_end, 500):
            y_k = sol.sol(rho_k)
            eint_vol_k, u_k = y_k
            nH = rho_k / F
            if nH <= 0 or eint_vol_k <= 0:
                continue
            log_nH = np.log10(nH)
            eint_nH = eint_vol_k / nH
            if eint_nH <= 0:
                continue
            c_k = self.eos.sound_speed(log_nH, np.log10(eint_nH))
            p_k = self.eos.pressure(log_nH, np.log10(eint_nH))
            err = abs(u_k + c_k - S)
            if err < best_err:
                best_err = err
                best_state = (rho_k, u_k, p_k, eint_vol_k)

        return best_state


# 
# CONVENIENCE: solve a Coleman-style test problem
# 

def solve_coleman_test(eos, rho_L, rho_R, v_L, v_R, T_L, T_R, x, t, x0=0.5):
    """
    Solve a Coleman-style Riemann problem given (rho, v, T) on each side.

    Uses the EoS to convert T -> eint, then solves and samples.
    Returns: rho, v, p, eint arrays at positions x and time t.
    """
    F = eos.nH2rhoFactor

    # Convert T -> eint using the EoS
    nH_L = rho_L / F
    nH_R = rho_R / F

    # Use the eint_from_T inverse table: given (nH, T) -> eint/nH
    # For now, invert by finding eint/nH such that T(nH, eint/nH) = T
    log_nH_L = np.log10(nH_L)
    log_nH_R = np.log10(nH_R)

    eint_nH_L = _find_eint_from_T(eos, log_nH_L, T_L)
    eint_nH_R = _find_eint_from_T(eos, log_nH_R, T_R)

    eint_L = nH_L * eint_nH_L  # internal energy density
    eint_R = nH_R * eint_nH_R

    # Compute pressures
    p_L = eos.pressure(log_nH_L, np.log10(eint_nH_L))
    p_R = eos.pressure(log_nH_R, np.log10(eint_nH_R))

    print(f"    L state: rho={rho_L:.4f}, v={v_L:.4f}, T={T_L:.6f}, "
          f"p={p_L:.6f}, eint={eint_L:.6f}")
    print(f"    R state: rho={rho_R:.4f}, v={v_R:.4f}, T={T_R:.6f}, "
          f"p={p_R:.6f}, eint={eint_R:.6f}")

    solver = GeneralEoSRiemann(eos)
    sol = solver.solve(rho_L, v_L, p_L, eint_L, rho_R, v_R, p_R, eint_R)

    print(f"    Star: p*={sol['p_star']:.6f}, u*={sol['u_star']:.6f}, "
          f"rho*_L={sol['rho_starL']:.4f}, rho*_R={sol['rho_starR']:.4f}")

    rho, v, p, eint = solver.sample(sol, x, t, x0,
                                     rho_L, v_L, p_L, eint_L,
                                     rho_R, v_R, p_R, eint_R)
    return rho, v, p, eint


def _find_eint_from_T(eos, log_nH, T_target):
    """Find eint/nH such that T(nH, eint/nH) = T_target."""
    # Analytic EoS: direct computation (no table inversion needed)
    if hasattr(eos, 'eint_nH_from_T'):
        return eos.eint_nH_from_T(log_nH, T_target)

    # Table-based EoS: bisect on the T table
    log_T_target = np.log10(T_target)
    lo = eos.T_tab.v2min
    hi = eos.T_tab.v2max

    def residual(log_eint_nH):
        log_T = eos.T_tab.interp(log_nH, log_eint_nH)
        return log_T - log_T_target

    try:
        log_eint_nH = brentq(residual, lo, hi, xtol=1e-12, rtol=1e-12)
    except ValueError:
        # If T is outside table range, clamp
        if residual(lo) > 0:
            log_eint_nH = lo
        else:
            log_eint_nH = hi

    return 10.0**log_eint_nH
