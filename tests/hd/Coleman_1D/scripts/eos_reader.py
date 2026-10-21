"""Slim LTE EoS reader for the Coleman_1D autotest.

Reads the entropy128 bicubic-Hermite table set produced by
`src/eos/entropy/generate_all_tables.py` and exposes the three queries
needed by the autotest's Snow two-gamma jump check:

  - eint_nH_from_T(log_nH, T_code)  : eint/nH at given (rho, T)
  - gamma1(log_nH, log_eint_nH)     : Gamma_1 at given (rho, eint)
  - temperature(log_nH, log_eint_nH): T at given (rho, eint)

The bicubic-Hermite kernel mirrors AMRVAC's mod_eos_entropy.t exactly,
so AMRVAC and this reader sample the same polynomial at table cells --
any mismatch is in the Fortran cons->prim path, not the table.

Use IdealGasEoS for the FI variants (gamma=5/3, no table lookups).
"""

from __future__ import annotations

import os
import struct
import numpy as np


# -- Constants matching AMRVAC's mod_global_parameters defaults ---------
MP_CGS = 1.672621777e-24
KB_CGS = 1.3806488e-16


# -- Bicubic-Hermite kernel ---------------------------------------------

def _cubic_basis(t):
    t2 = t * t
    t3 = t2 * t
    return (
        2 * t3 - 3 * t2 + 1.0,
        t3 - 2 * t2 + t,
        -2 * t3 + 3 * t2,
        t3 - t2,
    )


def _locate(nodes, val):
    n = len(nodes)
    if val <= nodes[0]:
        return 0, 0.0, nodes[1] - nodes[0]
    if val >= nodes[-1]:
        return n - 2, 1.0, nodes[-1] - nodes[-2]
    ix = int(np.searchsorted(nodes, val) - 1)
    ix = max(0, min(n - 2, ix))
    h = nodes[ix + 1] - nodes[ix]
    return ix, (val - nodes[ix]) / h, h


class _BicubicTable:
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

    def __init__(self, tables_dir, stem, composition):
        paths = [
            f'{tables_dir}/LTEeos_{stem}_{composition}_IonE.bin',
            f'{tables_dir}/LTEeos_{stem}_x_{composition}_IonE.bin',
            f'{tables_dir}/LTEeos_{stem}_y_{composition}_IonE.bin',
            f'{tables_dir}/LTEeos_{stem}_xy_{composition}_IonE.bin',
        ]
        self.f,   self.v1n, self.v2n, bounds = self._read_one(paths[0])
        self.fx,  _, _, _ = self._read_one(paths[1])
        self.fy,  _, _, _ = self._read_one(paths[2])
        self.fxy, _, _, _ = self._read_one(paths[3])
        self.v1_min, self.v1_max, self.v2_min, self.v2_max = bounds

    def shift(self, dx1, dx2):
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
        return (
            H0x*H0y*f[ix,iy]     + H2x*H0y*f[ix+1,iy]     + H0x*H2y*f[ix,iy+1]     + H2x*H2y*f[ix+1,iy+1]
          + H1x*H0y*fx[ix,iy]*dx + H3x*H0y*fx[ix+1,iy]*dx + H1x*H2y*fx[ix,iy+1]*dx + H3x*H2y*fx[ix+1,iy+1]*dx
          + H0x*H1y*fy[ix,iy]*dy + H2x*H1y*fy[ix+1,iy]*dy + H0x*H3y*fy[ix,iy+1]*dy + H2x*H3y*fy[ix+1,iy+1]*dy
          + H1x*H1y*fxy[ix,iy]*dx*dy + H3x*H1y*fxy[ix+1,iy]*dx*dy
          + H1x*H3y*fxy[ix,iy+1]*dx*dy + H3x*H3y*fxy[ix+1,iy+1]*dx*dy
        )


# -- Public EoS classes -------------------------------------------------

class IdealGasEoS:
    """gamma=5/3 ideal gas. No table lookups; matches mod_eos's FI dispatch."""

    def __init__(self, He_abundance=0.1, gamma=5.0/3.0,
                 unit_temperature=1e6, unit_numberdensity=1e15, unit_length=1e9):
        self.gamma = gamma
        self.He_abundance = He_abundance
        self.nH2rhoFactor = 1.0 + 4.0 * He_abundance
        self.unit_temperature   = unit_temperature
        self.unit_numberdensity = unit_numberdensity
        self.unit_density       = MP_CGS * unit_numberdensity * self.nH2rhoFactor
        self.unit_pressure      = KB_CGS * unit_numberdensity * unit_temperature
        self.unit_velocity      = np.sqrt(self.unit_pressure / self.unit_density)
        self.unit_length        = unit_length
        self.unit_time          = unit_length / self.unit_velocity

    def eint_nH_from_T(self, log_nH, T_code):
        # p = rho*T (eq_state_units=false). eint = p/(gamma-1) = rho*T/(gamma-1).
        # eint/nH = nH2rhoFactor * T / (gamma-1)
        return self.nH2rhoFactor * T_code / (self.gamma - 1.0)

    def gamma1(self, log_nH, log_eint_nH):
        return self.gamma

    def temperature(self, log_nH, log_eint_nH):
        # eint/nH = nH2rhoFactor * T / (gamma-1)  ->  T = (eint/nH) * (gamma-1) / nH2rhoFactor
        return (10.0 ** log_eint_nH) * (self.gamma - 1.0) / self.nH2rhoFactor

    def pressure(self, log_nH, log_eint_nH):
        # p = (gamma-1) * eint_vol = (gamma-1) * nH * (eint/nH)
        return (self.gamma - 1.0) * (10.0 ** log_nH) * (10.0 ** log_eint_nH)

    def sound_speed(self, log_nH, log_eint_nH):
        p = self.pressure(log_nH, log_eint_nH)
        rho = (10.0 ** log_nH) * self.nH2rhoFactor
        return np.sqrt(self.gamma * p / rho)

    def eint_from_p(self, log_nH, p):
        # eint_vol = p / (gamma-1)
        return p / (self.gamma - 1.0)


class EntropyEoS:
    """Entropy-method bicubic-Hermite LTE EoS reader. Reads tables in
    src/tables/eos_tables/entropy128/ via the same kernel that
    mod_eos_entropy.t uses at runtime."""

    def __init__(self, tables_dir, composition='HHe', He_abundance=0.1,
                 unit_temperature=1e6, unit_numberdensity=1e15, unit_length=1e9):
        self.He_abundance = He_abundance
        self.nH2rhoFactor = 1.0 + 4.0 * He_abundance
        self.unit_temperature   = unit_temperature
        self.unit_numberdensity = unit_numberdensity
        self.unit_density       = MP_CGS * unit_numberdensity * self.nH2rhoFactor
        self.unit_pressure      = KB_CGS * unit_numberdensity * unit_temperature
        self.unit_velocity      = np.sqrt(self.unit_pressure / self.unit_density)
        self.unit_length        = unit_length
        self.unit_time          = unit_length / self.unit_velocity

        # Code-unit log shifts. Subtract from CGS axes so subsequent queries
        # use code-unit log coordinates directly.
        log_nH_shift     = np.log10(unit_numberdensity)
        log_eint_nH_shift = np.log10(self.unit_pressure / unit_numberdensity)
        log_p_nH_shift    = log_eint_nH_shift
        log_T_shift       = np.log10(unit_temperature)

        # Forward tables on adaptive (log_nH, log_eint/nH) grid
        self.Tfwd  = _BicubicTable(tables_dir, 'Tfwd',  composition)
        self.pfwd  = _BicubicTable(tables_dir, 'pfwd',  composition)
        self.neOnH = _BicubicTable(tables_dir, 'neOnH', composition)
        self.g1e   = _BicubicTable(tables_dir, 'g1e',   composition)
        for tab in (self.Tfwd, self.pfwd, self.neOnH, self.g1e):
            tab.shift(log_nH_shift, log_eint_nH_shift)

        # Inverse-p tables on regular (log_nH, log_p/nH) grid
        self.eintP = _BicubicTable(tables_dir, 'eintP', composition)
        self.g1p   = _BicubicTable(tables_dir, 'g1p',   composition)
        for tab in (self.eintP, self.g1p):
            tab.shift(log_nH_shift, log_p_nH_shift)

        # Inverse-T table on regular (log_nH, log_T) grid
        self.eintT = _BicubicTable(tables_dir, 'eintT', composition)
        self.eintT.shift(log_nH_shift, log_T_shift)

    def eint_nH_from_T(self, log_nH, T_code):
        log_T = np.log10(T_code)
        log_e_nh_cgs = self.eintT.eval(log_nH, log_T)
        return 10.0 ** (log_e_nh_cgs - np.log10(self.unit_pressure / self.unit_numberdensity))

    def pressure(self, log_nH, log_eint_nH):
        p_cgs = self.pfwd.eval(log_nH, log_eint_nH)
        return p_cgs / self.unit_pressure

    def temperature(self, log_nH, log_eint_nH):
        return self.Tfwd.eval(log_nH, log_eint_nH) / self.unit_temperature

    def gamma1(self, log_nH, log_eint_nH):
        # Route through (p,nH) -> g1p, matching AMRVAC's set_output_vars
        p_code = self.pressure(log_nH, log_eint_nH)
        log_p_nH_code = np.log10(p_code / 10.0 ** log_nH)
        return self.g1p.eval(log_nH, log_p_nH_code)


# -- Factory ------------------------------------------------------------

def make_eos(eos_type: str, tables_dir: str, He_abundance: float):
    """Return an EoS instance for the given variant.

    eos_type in {'FI', 'LTE'}; tables_dir only matters for LTE."""
    if eos_type == 'FI':
        return IdealGasEoS(He_abundance=He_abundance)
    composition = 'HHe' if He_abundance > 0.0 else 'H'
    return EntropyEoS(tables_dir, composition=composition, He_abundance=He_abundance)
