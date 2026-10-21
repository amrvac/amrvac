"""LTE EoS Python wrapper for the *entropy method* table set.

Same public interface as `eos_riemann.LTEEoS` (pressure, temperature, gamma1,
eint_from_p, sound_speed, He_abundance, nH2rhoFactor), so it can drop in
under `_compare.make_lte_eos`. Reads the bicubic-Hermite tables produced by
`src/eos/entropy/generate_all_tables.py` (28 .bin files per composition,
4 per quantity: value + x + y + xy) and evaluates them with the same
kernel AMRVAC's mod_eos_entropy.t uses.

Forward queries (log_nH, log_e/nH):
  - Tfwd  -> T_cgs
  - neOnH -> y = n_e/n_H
  - pfwd  -> p_cgs
  - g1e   -> Gamma_1 at (rho, eps)

Inverse queries (log_nH, log_p/nH):
  - eintP -> log10(e_int/nH)_cgs
  - g1p   -> Gamma_1 at (rho, p)
"""

from __future__ import annotations

import os
import struct
import numpy as np


# -- Bicubic-Hermite kernel -- mirrors AMRVAC's mod_eos_entropy.t ----------

def _cubic_basis(t):
    t2 = t*t; t3 = t*t*t
    return (2*t3 - 3*t2 + 1.0,
            t3 - 2*t2 + t,
            -2*t3 + 3*t2,
            t3 - t2)


def _locate(nodes, val):
    """Find cell index ix such that nodes[ix] <= val <= nodes[ix+1]."""
    n = len(nodes)
    if val <= nodes[0]:
        return 0, 0.0, nodes[1] - nodes[0]
    if val >= nodes[-1]:
        return n - 2, 1.0, nodes[-1] - nodes[-2]
    ix = int(np.searchsorted(nodes, val) - 1)
    ix = max(0, min(n - 2, ix))
    h = nodes[ix + 1] - nodes[ix]
    t = (val - nodes[ix]) / h
    return ix, t, h


class _BicubicTable:
    """Holds f, fx, fy, fxy on a possibly-adaptive 2D grid (log axes)."""

    @staticmethod
    def _read_one(fname):
        with open(fname, 'rb') as fh:
            dim1, dim2 = struct.unpack('<ii', fh.read(8))
            v1_min, v1_max, v2_min, v2_max = struct.unpack('<4d', fh.read(32))
            tbl = np.frombuffer(fh.read(dim1 * dim2 * 8), dtype='<f8').reshape(dim2, dim1).T
            try:
                flag = struct.unpack('<i', fh.read(4))[0]
                if flag == 1:
                    v1n = np.frombuffer(fh.read(dim1 * 8), dtype='<f8').copy()
                    v2n = np.frombuffer(fh.read(dim2 * 8), dtype='<f8').copy()
                else:
                    v1n = np.linspace(v1_min, v1_max, dim1)
                    v2n = np.linspace(v2_min, v2_max, dim2)
            except struct.error:
                v1n = np.linspace(v1_min, v1_max, dim1)
                v2n = np.linspace(v2_min, v2_max, dim2)
        return tbl, v1n, v2n, (v1_min, v1_max, v2_min, v2_max)

    def __init__(self, tables_dir, stem, composition='HHe'):
        f_path  = f'{tables_dir}/LTEeos_{stem}_{composition}_IonE.bin'
        fx_path = f'{tables_dir}/LTEeos_{stem}_x_{composition}_IonE.bin'
        fy_path = f'{tables_dir}/LTEeos_{stem}_y_{composition}_IonE.bin'
        fxy_path = f'{tables_dir}/LTEeos_{stem}_xy_{composition}_IonE.bin'
        self.f,   self.v1n, self.v2n, bounds = self._read_one(f_path)
        self.fx,  _, _, _ = self._read_one(fx_path)
        self.fy,  _, _, _ = self._read_one(fy_path)
        self.fxy, _, _, _ = self._read_one(fxy_path)
        self.v1_min, self.v1_max, self.v2_min, self.v2_max = bounds

    def shift(self, dx1, dx2):
        """Subtract dx1/dx2 from axes so subsequent evals use code-unit coords."""
        self.v1n = self.v1n - dx1
        self.v2n = self.v2n - dx2
        self.v1_min -= dx1; self.v1_max -= dx1
        self.v2_min -= dx2; self.v2_max -= dx2

    def eval(self, x, y):
        ix, tx, dx = _locate(self.v1n, x)
        iy, ty, dy = _locate(self.v2n, y)
        H0x, H1x, H2x, H3x = _cubic_basis(tx)
        H0y, H1y, H2y, H3y = _cubic_basis(ty)
        f, fx, fy, fxy = self.f, self.fx, self.fy, self.fxy
        return (H0x*H0y*f[ix,iy] + H2x*H0y*f[ix+1,iy] + H0x*H2y*f[ix,iy+1] + H2x*H2y*f[ix+1,iy+1]
              + H1x*H0y*fx[ix,iy]*dx + H3x*H0y*fx[ix+1,iy]*dx + H1x*H2y*fx[ix,iy+1]*dx + H3x*H2y*fx[ix+1,iy+1]*dx
              + H0x*H1y*fy[ix,iy]*dy + H2x*H1y*fy[ix+1,iy]*dy + H0x*H3y*fy[ix,iy+1]*dy + H2x*H3y*fy[ix+1,iy+1]*dy
              + H1x*H1y*fxy[ix,iy]*dx*dy + H3x*H1y*fxy[ix+1,iy]*dx*dy + H1x*H3y*fxy[ix,iy+1]*dx*dy + H3x*H3y*fxy[ix+1,iy+1]*dx*dy)


# -- Public class with the same interface as eos_riemann.LTEEoS ----------

class EntropyEoS:
    """Drop-in replacement for `eos_riemann.LTEEoS` using entropy bicubic
    Hermite tables. All public methods work in AMRVAC code units, identical
    to `LTEEoS`.
    """

    def __init__(self, tables_dir, composition='HHe', He_abundance=0.1,
                 unit_temperature=1e6, unit_numberdensity=1e15, unit_length=1e9,
                 gamma1_approx=False):
        self.He_abundance = He_abundance
        self.nH2rhoFactor = 1.0 + 4.0 * He_abundance
        self.gamma1_approx = False

        mp_cgs = 1.672621777e-24
        kB_cgs = 1.3806488e-16
        self.unit_temperature = unit_temperature
        self.unit_numberdensity = unit_numberdensity
        self.unit_density = mp_cgs * unit_numberdensity
        self.unit_pressure = kB_cgs * unit_numberdensity * unit_temperature
        self.unit_velocity = np.sqrt(self.unit_pressure / self.unit_density)
        self.unit_length = unit_length
        self.unit_time = unit_length / self.unit_velocity

        # Code-unit log shifts: subtract these from CGS axes so the loaded
        # tables can be queried directly at code-unit log coordinates.
        log_nH_shift = np.log10(unit_numberdensity)
        # axis 2 of forward tables: log10(eint/nH) in CGS. The CGS unit of
        # eint/nH is erg/H = unit_pressure/unit_numberdensity in code terms.
        log_eint_nH_shift = np.log10(self.unit_pressure / unit_numberdensity)
        # axis 2 of inverse-p tables: log10(p/nH) -- same erg/H unit.
        log_p_nH_shift = log_eint_nH_shift
        log_T_shift = np.log10(unit_temperature)

        # -- Forward tables on adaptive (log_nH, log_eint/nH) grid --
        self.Tfwd_tab  = _BicubicTable(tables_dir, 'Tfwd',  composition)
        self.pfwd_tab  = _BicubicTable(tables_dir, 'pfwd',  composition)
        self.neOnH_tab = _BicubicTable(tables_dir, 'neOnH', composition)
        self.g1e_tab   = _BicubicTable(tables_dir, 'g1e',   composition)
        for tab in (self.Tfwd_tab, self.pfwd_tab, self.neOnH_tab, self.g1e_tab):
            tab.shift(log_nH_shift, log_eint_nH_shift)

        # -- Inverse-p tables on regular (log_nH, log_p/nH) grid --
        self.eintP_tab = _BicubicTable(tables_dir, 'eintP', composition)
        self.g1p_tab   = _BicubicTable(tables_dir, 'g1p',   composition)
        for tab in (self.eintP_tab, self.g1p_tab):
            tab.shift(log_nH_shift, log_p_nH_shift)

        # -- Inverse-T table on regular (log_nH, log_T) grid --
        self.eintT_tab = _BicubicTable(tables_dir, 'eintT', composition)
        self.eintT_tab.shift(log_nH_shift, log_T_shift)

        self._gamma1_source = 'entropy_bicubic_g1e'

    # -- Thermodynamic functions (code units) --

    def pressure(self, log_nH, log_eint_nH):
        """p in code units, from (rho, eint). Uses direct pfwd bicubic, which
        matches the runtime hd_get_pthermal LTE path exactly (since runtime
        also returns nH(1+He+y)T evaluated from the same node values)."""
        p_cgs = self.pfwd_tab.eval(log_nH, log_eint_nH)
        return p_cgs / self.unit_pressure

    def eint_from_p(self, log_nH, p):
        """eint_vol in code units from (rho, p)."""
        nH = 10.0 ** log_nH
        # log_p/nH in code units = log10(p_code / nH_code).
        log_p_nH = np.log10(max(p / nH, 1e-300))
        log_e_nh_cgs = self.eintP_tab.eval(log_nH, log_p_nH)
        # eintP stores log10(eint/nH) in CGS. Convert to code units.
        log_e_nh_code = log_e_nh_cgs - np.log10(self.unit_pressure / self.unit_numberdensity)
        e_nh_code = 10.0 ** log_e_nh_code
        return nH * e_nh_code

    def temperature(self, log_nH, log_eint_nH):
        """T in code units."""
        T_cgs = self.Tfwd_tab.eval(log_nH, log_eint_nH)
        return T_cgs / self.unit_temperature

    def gamma1(self, log_nH, log_eint_nH):
        """Gamma_1 at (rho, eint).

        Routes through pfwd -> g1p (the inverse-axis Gamma_1 table) so this method
        returns the SAME value AMRVAC's set_output_vars stores at the same
        state. Both AMRVAC's `gamma1_from_nH_p` (used to write the `g1` dat
        column) and this method now query the g1p table -- eliminating the
        forward/inverse interpolation-noise mismatch that made the Gamma_1 panel
        in plot_combined show non-overlapping sim and exact curves.
        """
        # p in code units at this (rho, eint) -- uses pfwd table.
        p_code = self.pressure(log_nH, log_eint_nH)
        nH_code = 10.0 ** log_nH
        log_p_nH_code = np.log10(p_code / nH_code)
        return self.g1p_tab.eval(log_nH, log_p_nH_code)

    def gamma1_from_g1e(self, log_nH, log_eint_nH):
        """Direct forward Gamma_1(rho,eps) via the g1e table -- kept for diagnostics
        (e.g. comparing g1p vs g1e interpolation noise at the same state)."""
        return self.g1e_tab.eval(log_nH, log_eint_nH)

    def sound_speed(self, log_nH, log_eint_nH):
        g1 = self.gamma1(log_nH, log_eint_nH)
        p = self.pressure(log_nH, log_eint_nH)
        rho = 10.0 ** log_nH * self.nH2rhoFactor
        return np.sqrt(g1 * p / rho)

    def eint_nH_from_T(self, log_nH, T_code):
        """eint/nH in code units, from (rho, T). Direct eintT bicubic lookup
        -- no bisection. Matches the analytic-EoS branch of _find_eint_from_T
        so the exact tabular Riemann solver invokes this directly instead of
        bisecting on T_tab (which entropy_eos doesn't expose)."""
        log_T = np.log10(T_code)
        log_e_nh_cgs = self.eintT_tab.eval(log_nH, log_T)
        log_e_nh_code = log_e_nh_cgs - np.log10(self.unit_pressure / self.unit_numberdensity)
        return 10.0 ** log_e_nh_code
