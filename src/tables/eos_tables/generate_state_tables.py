"""
Generate LTE equation-of-state tables for AMRVAC.

Solves coupled Saha ionization equilibrium + energy conservation for H+He plasma
using a damped Newton root finder (JAX + optimistix). Produces binary table files
for four variants:
  - H_NoIonE   : pure H, no ionization energy in eint
  - H_IonE     : pure H, ionization energy included
  - HHe_NoIonE : H + 10% He (by number), no ionization energy
  - HHe_IonE   : H + 10% He, ionization energy included

Each variant produces 2-4 binary files:
  - LTEeos_neOnH_{comp}.bin  : ne/nH(nH, eint/nH)    [dimensionless]
  - LTEeos_T_{comp}.bin      : T(nH, eint/nH)         [K]
  - LTEeos_P_{comp}.bin      : log10(p/nH)(nH, eint/nH)  [log10(erg) CGS or log10(J) SI]
  - LTEeos_p2eint_{comp}.bin : eint/p(nH, p/nH)       [dimensionless, IonE only]

Binary format per file (matches AMRVAC read_eos_from_file):
  [2 x int32]   : (dim_nH, dim_var2)     -- grid dimensions
  [4 x float64] : (nH_min, nH_max, var2_min, var2_max)  -- log10 bounds
  [dim_nH x dim_var2 x float64] : table data in Fortran column-major order

Grid: 128x128 in log10 space (SI).
  nH:      [1e11, 1e25] m^-3
  eint/nH: [1e-19.8, 1e-16] J
  p/nH:    [1e-20.2, 1e-16.4] Pa m^3

Unit conventions for bounds:
  CGS (default): nH in cm^-3, eint/nH and p/nH in erg (= dyn cm)
  SI:            nH in m^-3,  eint/nH in J, p/nH in Pa m^3

Based on jax.ipynb by Chris Osborne.
"""

import argparse
import os
import struct
import sys

# Force CPU backend (avoids CUDA driver version mismatches and launch overhead)
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

import astropy.constants as const
import astropy.units as u
from tqdm import tqdm
import optimistix
from optimistix import root_find

from optimistix._solver.newton_chord import (
    _AbstractNewtonChord,
    _NewtonChordState,
    _small,
    _diverged,
    _converged,
)
import lineax as lx
from optimistix._custom_types import Aux, Fn, Out, Y
from optimistix._solution import RESULTS
from optimistix._misc import cauchy_termination
from equinox.internal import ω
import jax.tree_util as jtu
from jax import lax
from typing import Any
from jaxtyping import PyTree

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class DampedNewton(_AbstractNewtonChord[Y, Out, Aux]):
    """Newton's method for root finding. Also sometimes known as Newton--Raphson.

    Unlike the SciPy implementation of Newton's method, the Optimistix version also
    works for vector-valued (or PyTree-valued) `y`.

    This solver optionally accepts the following `options`:

    - `lower`: The lower bound on the hypercube which contains the root. Should be a
        PyTree of arrays each broadcastable to the corresponding element of `y`. The
        iterates of `y` will be clipped to this hypercube.
    - `upper`: The upper bound on the hypercube which contains the root. Should be a
        PyTree of arrays each broadcastable to the corresponding element of `y`. The
        iterates of `y` will be clipped to this hypercube.
    """

    _is_newton = True
    cauchy_termination = True

    def step(self, fn, y, args, options, state, tags):
        lower = options.get("lower")
        upper = options.get("upper")
        del options
        if self._is_newton:
            fx, lin_fn, aux = jax.linearize(
                lambda _y: fn(_y, args), y, has_aux=True
            )
            jac = lx.FunctionLinearOperator(
                lin_fn, jax.eval_shape(lambda: y), tags=tags
            )
            sol = lx.linear_solve(jac, fx, self.linear_solver, throw=False)
        else:
            fx, aux = fn(y, args)
            jac, linear_state = state.linear_state   # pyright: ignore
            linear_state = lax.stop_gradient(linear_state)
            sol = lx.linear_solve(
                jac, fx, self.linear_solver, state=linear_state, throw=False
            )
        diff = sol.value
        new_y = (y**ω - 0.1 * diff**ω).ω
        if lower is not None:
            new_y = jtu.tree_map(lambda a, b: jnp.clip(a, min=b), new_y, lower)
        if upper is not None:
            new_y = jtu.tree_map(lambda a, b: jnp.clip(a, max=b), new_y, upper)
        if lower is not None or upper is not None:
            diff = (y**ω - new_y**ω).ω
        scale = (self.atol + self.rtol * ω(new_y).call(jnp.abs)).ω
        with jax.numpy_dtype_promotion("standard"):
            diffsize = self.norm((diff**ω / scale**ω).ω)
        f_val = fx if self.cauchy_termination else None
        new_state = _NewtonChordState(
            f=f_val,
            linear_state=state.linear_state,
            diff=diff,
            diffsize=jnp.asarray(diffsize, dtype=state.diffsize.dtype),
            diffsize_prev=state.diffsize,
            result=RESULTS.promote(sol.result),
            step=state.step + 1,
        )
        return new_y, new_state, aux

    def terminate(self, fn, y, args, options, state, tags):
        del fn, args, options
        if self.cauchy_termination:
            # NOTE(cmo): As we struggle with numerical precision/boundaries in a
            # few spots, do a classic Cauchy termination rather than just
            # comparing f against 0 as done in the base Newton root finder
            # ... but maybe we should terminate on either
            terminate_a = cauchy_termination(
                self.rtol, self.atol, self.norm, y, state.diff, state.f, state.diff
            )
            terminate_b = cauchy_termination(
                self.rtol, self.atol, self.norm, y, state.diff, state.f, state.diff
            )
            terminate = terminate_a | terminate_b
            terminate_result = RESULTS.successful
        else:
            # TODO(kidger): perform only one iteration when solving a linear system!
            at_least_two = state.step >= 2
            rate = state.diffsize / state.diffsize_prev
            factor = state.diffsize * rate / (1 - rate)
            small = _small(state.diffsize)
            diverged = _diverged(rate)
            converged = _converged(factor, self.kappa)
            terminate = at_least_two & (small | diverged | converged)
            terminate_result = RESULTS.where(
                jnp.invert(small) & (diverged | jnp.invert(converged)),
                RESULTS.nonlinear_divergence,
                RESULTS.successful,
            )
        linsolve_fail = state.result != RESULTS.successful
        result = RESULTS.where(linsolve_fail, state.result, terminate_result)
        terminate = linsolve_fail | terminate
        return terminate, result


chi_H = (13.595 << u.eV).to("J").value  # H ionization energy
chi_HeII = (24.48 << u.eV).to("J").value  # He first ionization energy (He0 -> He+)
chi_HeIII = (54.17 << u.eV).to("J").value  # He second ionization energy (He+ -> He2+)

kB = const.k_B.value
m_e = const.m_e.value
h_planck = const.h.value

# SI -> CGS log10 offsets
# nH: 1 m^-3 = 1e-6 cm^-3  =>  log10(nH_cgs) = log10(nH_SI) - 6
# energy/nH: 1 J = 1e7 erg  =>  log10(E_cgs) = log10(E_SI) + 7
LOG10_NH_SI_TO_CGS = -6.0
LOG10_ENERGY_SI_TO_CGS = 7.0

NH_BOUNDS = [11, 25]  # log10(nH / m^-3)
EINT_NH_BOUNDS = [-19.8, -16]  # log10(eint_per_nH / J)
PRES_NH_BOUNDS = [-20.2, -16.4]  # log10(pres_per_nH / Pa m^3)
N_GRID = 256

# FI blending: smooth table convergence to fully-ionised ideal gas above T_blend.
# Eliminates the discontinuity at the bypass threshold in AMRVAC.
FI_BLEND_T_LO = 50e3    # K: blend starts (f=0, pure Saha)
FI_BLEND_T_HI = 200e3   # K: blend ends   (f=1, pure FI)


def fi_blend_fraction(T_K, T_lo=FI_BLEND_T_LO, T_hi=FI_BLEND_T_HI):
    """Smooth [0,1] blend fraction: 0 below T_lo, 1 above T_hi.
    Hermite smoothstep (C1) in log10(T) space."""
    t = np.clip((np.log10(np.maximum(T_K, 1.0)) - np.log10(T_lo))
                / (np.log10(T_hi) - np.log10(T_lo)), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def fi_reference_values(A_He):
    """Return FI analytical constants for given helium abundance."""
    n_per_nH = 2.0 + 3.0 * A_He       # total particles per nH
    neOnH = 1.0 + 2.0 * A_He           # electron fraction
    eion_per_nH = chi_H + A_He * (chi_HeII + chi_HeIII)  # ionisation energy per nH [J]
    return n_per_nH, neOnH, eion_per_nH


def saha_rhs_H(T):
    """Saha RHS for H: n_HII * n_e / n_HI.  Statistical weight ratio 2*g_HII/g_HI = 1."""
    saha_const = ((2 * jnp.pi * m_e * kB) / h_planck**2) ** 1.5
    return saha_const * T**1.5 * jnp.exp(-chi_H / (kB * T))


def saha_rhs_HeII(T):
    """Saha RHS for He I -> He II: n_HeII * n_e / n_HeI.  Weight ratio 2*g_HeII/g_HeI = 4."""
    saha_const = ((2 * jnp.pi * m_e * kB) / h_planck**2) ** 1.5
    return 4 * saha_const * T**1.5 * jnp.exp(-chi_HeII / (kB * T))


def saha_rhs_HeIII(T):
    """Saha RHS for He II -> He III: n_HeIII * n_e / n_HeII.  Weight ratio 2*g_HeIII/g_HeII = 1."""
    saha_const = ((2 * jnp.pi * m_e * kB) / h_planck**2) ** 1.5
    return saha_const * T**1.5 * jnp.exp(-chi_HeIII / (kB * T))


# NOTE(cmo): When everything balances this returns [0, 0] - we want to differentiate it and do an NR
def ion_eq_H(ion_frac, temperature, nhtot, eint, include_ion_e=True):
    """Pure-H Saha + energy balance.  2 equations in (ne/nH, T)."""
    saha_const = saha_rhs_H(temperature) / nhtot
    ne = ion_frac * nhtot
    ion_mult = 1.0 if include_ion_e else 0.0

    saha_term = jnp.where(
        temperature < 5e3,
        ion_frac**2 + saha_const * ion_frac - saha_const,
        ion_frac + saha_const * (1.0 - 1.0 / ion_frac),
    )
    return jnp.array(
        [
            saha_term,
            1.5 * nhtot * (1.0 + ion_frac) * kB * temperature
            + ion_mult * ne * chi_H
            - eint,
        ]
    )


def ion_eq_HHe(ion_frac, temperature, hei_frac, heii_frac, nhtot, eint, include_ion_e=True):
    """H+He reduced system.  4 equations in (ne/nH, T, He_I_frac, He_II_frac).

    He III fraction is derived: heiii_frac = 1 - hei_frac - heii_frac.
    H fractions are derived from H Saha given ne.
    """
    saha_const_H = saha_rhs_H(temperature)
    saha_const_HeII = saha_rhs_HeII(temperature)
    saha_const_HeIII = saha_rhs_HeIII(temperature)
    ne = ion_frac * nhtot

    hii_hi = saha_const_H / ne
    hi_frac = 1.0 / (1.0 + hii_hi)
    hii_frac = 1.0 - hi_frac
    heiii_frac = 1.0 - hei_frac - heii_frac

    ion_mult = 1.0 if include_ion_e else 0.0

    return jnp.array(
        [
            heii_frac / hei_frac - saha_const_HeII / ne,
            heiii_frac / heii_frac - saha_const_HeIII / ne,
            1.5 * nhtot * (1.1 + hii_frac + 0.1 * (heii_frac + 2 * heiii_frac)) * kB * temperature
            + ion_mult * (
                nhtot * hii_frac * chi_H
                + 0.1 * nhtot * (heii_frac * chi_HeII + heiii_frac * (chi_HeII + chi_HeIII))
            )
            - eint,
            hii_frac + 0.1 * (heii_frac + 2 * heiii_frac) - ion_frac,
        ]
    )


def ion_eq_H_pres(ion_frac, temperature, nhtot, pres):
    """Pure-H Saha + ideal-gas pressure.  2 equations in (ne/nH, T)."""
    saha_const = saha_rhs_H(temperature) / nhtot

    saha_term = jnp.where(
        temperature < 5e3,
        ion_frac**2 + saha_const * ion_frac - saha_const,
        ion_frac + saha_const * (1.0 - 1.0 / ion_frac),
    )
    return jnp.array(
        [
            saha_term,
            nhtot * (1.0 + ion_frac) * kB * temperature - pres,
        ]
    )


def ion_eq_HHe_pres(ion_frac, temperature, hei_frac, heii_frac, nhtot, pres):
    """H+He reduced system with pressure constraint (no ionization energy in p)."""
    saha_const_H = saha_rhs_H(temperature)
    saha_const_HeII = saha_rhs_HeII(temperature)
    saha_const_HeIII = saha_rhs_HeIII(temperature)
    ne = ion_frac * nhtot

    hii_hi = saha_const_H / ne
    hi_frac = 1.0 / (1.0 + hii_hi)
    hii_frac = 1.0 - hi_frac
    heiii_frac = 1.0 - hei_frac - heii_frac

    return jnp.array(
        [
            heii_frac / hei_frac - saha_const_HeII / ne,
            heiii_frac / heii_frac - saha_const_HeIII / ne,
            nhtot * (1.1 + hii_frac + 0.1 * (heii_frac + 2 * heiii_frac)) * kB * temperature - pres,
            hii_frac + 0.1 * (heii_frac + 2 * heiii_frac) - ion_frac,
        ]
    )


INCLUDE_ION_E = True


def fn_eint_H(y, args):
    nhtot, eint = args
    return ion_eq_H(y[0], y[1], nhtot, eint, include_ion_e=INCLUDE_ION_E)


def fn_eint_HHe(y, args):
    nhtot, eint = args
    return ion_eq_HHe(y[0], y[1], y[2], y[3], nhtot, eint, include_ion_e=INCLUDE_ION_E)


def fn_pres_H(y, args):
    nhtot, pres = args
    return ion_eq_H_pres(y[0], y[1], nhtot, pres)


def fn_pres_HHe(y, args):
    nhtot, pres = args
    return ion_eq_HHe_pres(y[0], y[1], y[2], y[3], nhtot, pres)


_solver = DampedNewton(rtol=1e-8, atol=1e-12)


def solve_eint(nhtot, eint, H_only=False):
    """Solve ionization equilibrium given (nH, eint).

    Returns (ne/nH, T) for H-only, or (ne/nH, T, He_I_frac, He_II_frac) for H+He.
    """
    # Pass 1: H-only
    sol_H = root_find(
        fn_eint_H,
        solver=_solver,
        y0=jnp.array([0.5, 8e3]),
        args=(nhtot, eint),
        throw=False,
        options=dict(lower=jnp.array([1e-30, 100.0]), upper=jnp.array([1.0, 1e8])),
        max_steps=2000,
    )
    y, T = sol_H.value

    if H_only:
        return sol_H.value

    # Pass 2: H+He using H-only solution as initial guess
    heii_hei = float(saha_rhs_HeII(T)) / (y * nhtot)
    heiii_heii = float(saha_rhs_HeIII(T)) / (y * nhtot)
    hei_frac = 1.0 / (1.0 + heii_hei + heii_hei * heiii_heii)
    heii_frac = heii_hei * hei_frac
    heiii_frac = heiii_heii * heii_frac

    y0 = np.array([y + 0.1 * heii_frac + 0.2 * heiii_frac, T, hei_frac, heii_frac])
    sol = root_find(
        fn_eint_HHe,
        solver=_solver,
        y0=y0,
        args=(nhtot, eint),
        throw=False,
        options=dict(
            lower=jnp.array([1e-30, 100.0, 1e-30, 1e-30]),
            upper=jnp.array([1.2, 1e8, 1.0, 1.0]),
        ),
        max_steps=1000,
    )
    return sol.value


def solve_pres(nhtot, pres, H_only=False):
    """Solve ionization equilibrium given (nH, p).  Returns eint/p ratio."""
    # Pass 1: H-only
    sol_H = root_find(
        fn_pres_H,
        solver=_solver,
        y0=jnp.array([0.5, 8e3]),
        args=(nhtot, pres),
        throw=False,
        options=dict(lower=jnp.array([1e-30, 100.0]), upper=jnp.array([1.0, 1e8])),
        max_steps=2000,
    )
    y, T = sol_H.value

    if H_only:
        eint = 1.5 * nhtot * (1.0 + y) * kB * T + nhtot * y * chi_H
        return eint / pres

    # Pass 2: H+He
    heii_hei = float(saha_rhs_HeII(T)) / (y * nhtot)
    heiii_heii = float(saha_rhs_HeIII(T)) / (y * nhtot)
    hei_frac = 1.0 / (1.0 + heii_hei + heii_hei * heiii_heii)
    heii_frac = heii_hei * hei_frac
    heiii_frac = heiii_heii * heii_frac

    y0 = np.array([y + 0.1 * heii_frac + 0.2 * heiii_frac, T, hei_frac, heii_frac])
    sol = root_find(
        fn_pres_HHe,
        solver=_solver,
        y0=y0,
        args=(nhtot, pres),
        throw=False,
        options=dict(
            lower=jnp.array([1e-30, 100.0, 1e-30, 1e-30]),
            upper=jnp.array([1.2, 1e8, 1.0, 1.0]),
        ),
        max_steps=1000,
    )

    y_sol = sol.value[0]
    T_sol = sol.value[1]
    hei_frac = sol.value[2]
    heii_frac = sol.value[3]

    ne = y_sol * nhtot
    hii_hi = saha_rhs_H(T_sol) / ne
    hi_frac = 1.0 / (1.0 + hii_hi)
    hii_frac = 1.0 - hi_frac
    heiii_frac = 1.0 - hei_frac - heii_frac

    eint = (
        1.5 * nhtot * (1.1 + hii_frac + 0.1 * (heii_frac + 2 * heiii_frac)) * kB * T_sol
        + nhtot * hii_frac * chi_H
        + 0.1 * nhtot * (heii_frac * chi_HeII + heiii_frac * (chi_HeII + chi_HeIII))
    )
    return eint / pres


def solve_at_T(nH, T, A_He, include_ion_e):
    """Forward Saha solve at fixed (nH, T). Returns eint/nH in SI (J).

    For H-only: closed-form quadratic for ionization fraction.
    For H+He: fixed-point iteration on coupled Saha equations (~5 iterations).
    """
    K_H = float(saha_rhs_H(jnp.float64(T)))

    if A_He == 0.0:  # H-only
        S = K_H / nH
        y = 2.0 * S / (S + np.sqrt(S**2 + 4.0 * S))
        eint_nH = 1.5 * (1.0 + y) * kB * T
        if include_ion_e:
            eint_nH += y * chi_H
        return eint_nH

    # H+He: fixed-point iteration
    K_HeII = float(saha_rhs_HeII(jnp.float64(T)))
    K_HeIII = float(saha_rhs_HeIII(jnp.float64(T)))

    # Initial guess from H-only Saha
    S = K_H / nH
    y = (-S + np.sqrt(S**2 + 4.0 * S)) / 2.0
    ne = max(y * nH, 1e-30)

    for _ in range(20):
        hii_hi = K_H / ne
        hi_frac = 1.0 / (1.0 + hii_hi)
        hii_frac = 1.0 - hi_frac

        heii_hei = K_HeII / ne
        heiii_heii = K_HeIII / ne
        hei_frac = 1.0 / (1.0 + heii_hei * (1.0 + heiii_heii))
        heii_frac = heii_hei * hei_frac
        heiii_frac = heiii_heii * heii_frac

        y_new = hii_frac + 0.1 * (heii_frac + 2.0 * heiii_frac)
        ne_new = max(y_new * nH, 1e-30)

        if abs(ne_new - ne) / ne < 1e-12:
            ne = ne_new
            break
        ne = ne_new

    # Recompute final fractions at converged ne
    hii_hi = K_H / ne
    hi_frac = 1.0 / (1.0 + hii_hi)
    hii_frac = 1.0 - hi_frac
    heii_hei = K_HeII / ne
    heiii_heii = K_HeIII / ne
    hei_frac = 1.0 / (1.0 + heii_hei * (1.0 + heiii_heii))
    heii_frac = heii_hei * hei_frac
    heiii_frac = heiii_heii * heii_frac

    # eint/nH: kinetic + ionization energy
    eint_nH = 1.5 * (1.1 + hii_frac + 0.1 * (heii_frac + 2.0 * heiii_frac)) * kB * T
    if include_ion_e:
        eint_nH += (hii_frac * chi_H
                    + 0.1 * (heii_frac * chi_HeII
                             + heiii_frac * (chi_HeII + chi_HeIII)))
    return eint_nH


def _log10_p_H_func(log10_nH, log10_eint_nH, y0):
    """Differentiable log10(p) for H-only via root_find + IFT."""
    nH = 10.0**log10_nH
    eint_nH = 10.0**log10_eint_nH
    eint = eint_nH * nH

    y0_sg = jax.lax.stop_gradient(y0)
    sol = root_find(
        fn_eint_H, solver=_solver, y0=y0_sg, args=(nH, eint),
        throw=False,
        options=dict(lower=jnp.array([1e-30, 100.0]),
                     upper=jnp.array([1.0, 1e8])),
        max_steps=100,
    )
    y, T = sol.value[0], sol.value[1]
    p = nH * (1.0 + y) * kB * T
    return jnp.log10(p)


def _log10_p_HHe_func(log10_nH, log10_eint_nH, y0):
    """Differentiable log10(p) for H+He via root_find + IFT."""
    nH = 10.0**log10_nH
    eint_nH = 10.0**log10_eint_nH
    eint = eint_nH * nH

    y0_sg = jax.lax.stop_gradient(y0)
    sol = root_find(
        fn_eint_HHe, solver=_solver, y0=y0_sg, args=(nH, eint),
        throw=False,
        options=dict(lower=jnp.array([1e-30, 100.0, 1e-30, 1e-30]),
                     upper=jnp.array([1.2, 1e8, 1.0, 1.0])),
        max_steps=100,
    )
    y = sol.value[0]
    T = sol.value[1]
    # p = nH * (1 + A_He + ne/nH) * kB * T; A_He = 0.1 for H+He
    p = nH * (1.1 + y) * kB * T
    return jnp.log10(p)


_grad_log10_p_H = jax.jit(jax.grad(_log10_p_H_func, argnums=(0, 1)))
_grad_log10_p_HHe = jax.jit(jax.grad(_log10_p_HHe_func, argnums=(0, 1)))



def generate_eint_tables(n_grid, nh_grid, eint_nh_grid, H_only=False, blend=True):
    """Generate (ne/nH, T) tables on the (nH, eint/nH) grid.

    If blend=True and IonE is enabled, smoothly blend toward FI analytical
    values above T_blend, so the tables converge to ideal gas at the hot end.
    """
    A_He = 0.0 if H_only else 0.1
    n_per_nH_FI, neOnH_FI, eion_per_nH = fi_reference_values(A_He)
    do_blend = blend and INCLUDE_ION_E

    y_table = np.zeros((n_grid, n_grid))
    T_table = np.zeros((n_grid, n_grid))
    n_blended = 0
    for inh, n_h in enumerate(tqdm(nh_grid, desc=f"eint table ({'H' if H_only else 'HHe'})")):
        for ie, eint_nh in enumerate(eint_nh_grid):
            result = solve_eint(n_h, eint_nh * n_h, H_only=H_only)
            y_saha = result[0]
            T_saha = result[1]

            if do_blend:
                f = fi_blend_fraction(T_saha)
                if f > 0:
                    n_blended += 1
                    T_FI = max((eint_nh - eion_per_nH) / (1.5 * n_per_nH_FI * kB), 100.0)
                    # Geometric interpolation (log-space blend)
                    T_table[inh, ie] = T_saha**(1-f) * T_FI**f
                    y_table[inh, ie] = y_saha**(1-f) * neOnH_FI**f
                else:
                    T_table[inh, ie] = T_saha
                    y_table[inh, ie] = y_saha
            else:
                T_table[inh, ie] = T_saha
                y_table[inh, ie] = y_saha

    if do_blend:
        print(f"  FI blend applied to {n_blended}/{n_grid*n_grid} grid points")
    return y_table, T_table


def generate_pres_table_independent(n_grid, nh_grid, pres_nh_grid,
                                    H_only=False, blend=True):
    """Generate eint/p table by independently solving Saha on the (nH, p/nH) grid.

    This is the LEGACY method. The round-trip e->p->e is not exact because the
    Saha solver produces slightly different (T, y) values than the forward tables,
    and the tables sample different grids with independent interpolation errors.

    If blend=True, smoothly blend toward FI eint/p = 1.5 + eion/p above T_blend.
    """
    A_He = 0.0 if H_only else 0.1
    n_per_nH_FI, neOnH_FI, eion_per_nH = fi_reference_values(A_He)
    do_blend = blend and INCLUDE_ION_E

    eint_pres_table = np.zeros((n_grid, n_grid))
    n_blended = 0
    for inh, n_h in enumerate(tqdm(nh_grid, desc=f"pres table ({'H' if H_only else 'HHe'})")):
        for ip, pres_nh in enumerate(pres_nh_grid):
            ratio = solve_pres(n_h, pres_nh * n_h, H_only=H_only)
            ratio_saha = float(ratio) if np.ndim(ratio) == 0 else float(ratio[0])

            if do_blend:
                # Estimate T from FI formula for blend fraction
                T_approx = pres_nh / (n_per_nH_FI * kB)
                f = fi_blend_fraction(T_approx)
                if f > 0:
                    n_blended += 1
                    ratio_FI = 1.5 + eion_per_nH / pres_nh
                    ratio_saha = (1.0 - f) * ratio_saha + f * ratio_FI
            eint_pres_table[inh, ip] = ratio_saha

    eint_pres_table[np.isnan(eint_pres_table)] = 1.5
    if do_blend:
        print(f"  FI blend (pres): {n_blended}/{n_grid*n_grid} grid points")
    return eint_pres_table


def generate_pres_table(y_table, T_table, n_grid, nh_grid, eint_nh_grid,
                        pres_nh_grid, H_only=False, blend=True):
    """Derive eint/p table from the forward (T, neOnH) tables.

    Instead of independently solving Saha on the pressure grid, numerically
    invert the forward tables using PCHIP interpolation. This ensures the
    round-trip e->p->e is exact at forward grid nodes and O(h^4) between them.

    The forward tables already have FI blending applied, so the derived p
    values automatically incorporate the blend — no separate blend step needed.

    Falls back to solve_pres() for pressure values outside the range of the
    forward table's pressure coverage (typically only a few deep-chromosphere
    points at the cold end).
    """
    from scipy.interpolate import PchipInterpolator

    A_He = 0.0 if H_only else 0.1
    eint_pres_table = np.zeros((n_grid, n_grid))
    log_eint_nh = np.log10(eint_nh_grid)
    log_pres_nh = np.log10(pres_nh_grid)
    n_fallback = 0
    n_monotone_fix = 0

    for i in tqdm(range(n_grid), desc=f"p2eint ({'H' if H_only else 'HHe'})"):
        # Compute p/nH at each forward grid node from the stored T and y.
        # This is EXACTLY the same computation as the P table (line 799),
        # so the values are guaranteed consistent with the forward tables.
        p_over_nH = (1.0 + A_He + y_table[i, :]) * kB * T_table[i, :]
        log_p_nH = np.log10(np.maximum(p_over_nH, 1e-100))

        # Check and enforce monotonicity (physics guarantees p increases with
        # eint at fixed nH, but numerical noise in the Saha solver can create
        # tiny non-monotonic blips)
        diffs = np.diff(log_p_nH)
        if np.any(diffs <= 0):
            n_bad = np.sum(diffs <= 0)
            n_monotone_fix += n_bad
            # Forward accumulate: ensure strictly increasing
            for j in range(1, len(log_p_nH)):
                if log_p_nH[j] <= log_p_nH[j-1]:
                    log_p_nH[j] = log_p_nH[j-1] + 1e-12

        # Build inverse mapping: log(p/nH) -> log(eint/nH) via PCHIP
        inv_interp = PchipInterpolator(log_p_nH, log_eint_nh)

        # Range covered by forward table for this nH slice
        p_min, p_max = log_p_nH[0], log_p_nH[-1]

        for k in range(n_grid):
            if p_min <= log_pres_nh[k] <= p_max:
                # In range: use PCHIP inversion of forward table
                log_eint_k = float(inv_interp(log_pres_nh[k]))
                eint_pres_table[i, k] = 10.0**log_eint_k / pres_nh_grid[k]
            else:
                # Out of range: fall back to independent Saha solve
                ratio = solve_pres(nh_grid[i], pres_nh_grid[k] * nh_grid[i],
                                   H_only=H_only)
                ratio_val = float(ratio) if np.ndim(ratio) == 0 \
                            else float(ratio[0])
                # Apply FI blend if the fallback point is in the hot regime
                if blend and INCLUDE_ION_E:
                    _, _, eion_per_nH = fi_reference_values(A_He)
                    n_per_nH_FI = 2.3 if A_He > 0 else 2.0
                    T_approx = pres_nh_grid[k] / (n_per_nH_FI * kB)
                    f = fi_blend_fraction(T_approx)
                    if f > 0:
                        ratio_FI = 1.5 + eion_per_nH / pres_nh_grid[k]
                        ratio_val = (1.0 - f) * ratio_val + f * ratio_FI
                eint_pres_table[i, k] = ratio_val
                n_fallback += 1

    eint_pres_table[np.isnan(eint_pres_table)] = 1.5

    print(f"  p2eint derived from forward tables (round-trip consistent)")
    if n_fallback > 0:
        print(f"  Saha fallback: {n_fallback}/{n_grid*n_grid} points "
              f"({100*n_fallback/(n_grid*n_grid):.1f}%)")
    if n_monotone_fix > 0:
        print(f"  Monotonicity fixes: {n_monotone_fix} intervals")
    return eint_pres_table


def compute_gamma1(y_table, T_table, nh_grid, eint_nh_grid, A_He, H_only,
                   blend=True):
    """Compute Gamma1 table using JAX autodiff (implicit function theorem).

    Uses jax.grad through optimistix.root_find to get exact partial derivatives
    of log10(p) w.r.t. (log10 nH, log10 eint/nH). Zero finite-difference error.

    Gamma1 = a + b * (p / eint_vol) where:
      a = d(log10 p) / d(log10 nH)       at constant eint/nH
      b = d(log10 p) / d(log10 eint/nH)  at constant nH

    If blend=True, smoothly blend gamma1 toward 5/3 above T_blend.
    """
    n_nH, n_eint = y_table.shape
    gamma1_table = np.zeros((n_nH, n_eint))
    grad_func = _grad_log10_p_H if H_only else _grad_log10_p_HHe
    g1_min, g1_max = 1e30, -1e30
    gamma_mono = 5.0 / 3.0  # monatomic ideal gas

    for i in tqdm(range(n_nH), desc=f"gamma1 ({'H' if H_only else 'HHe'})"):
        log10_nH = jnp.float64(np.log10(nh_grid[i]))
        for j in range(n_eint):
            log10_eint_nH = jnp.float64(np.log10(eint_nh_grid[j]))

            if H_only:
                y0 = jnp.array([y_table[i, j], T_table[i, j]])
            else:
                y_val = y_table[i, j]
                T_val = T_table[i, j]
                ne = max(y_val * nh_grid[i], 1e-30)
                K_HeII = float(saha_rhs_HeII(jnp.float64(T_val)))
                K_HeIII = float(saha_rhs_HeIII(jnp.float64(T_val)))
                heii_hei = K_HeII / ne
                heiii_heii = K_HeIII / ne
                hei_frac = 1.0 / (1.0 + heii_hei * (1.0 + heiii_heii))
                heii_frac = heii_hei * hei_frac
                y0 = jnp.array([y_val, T_val, hei_frac, heii_frac])

            a, b = grad_func(log10_nH, log10_eint_nH, y0)
            p_over_eint = ((1.0 + A_He + y_table[i, j]) * kB * T_table[i, j]
                           / eint_nh_grid[j])
            g1 = float(a) + float(b) * p_over_eint

            # Blend toward 5/3 at hot end (Gamma1 cannot exceed 5/3 for H+He)
            if blend:
                f = fi_blend_fraction(T_table[i, j])
                g1 = (1.0 - f) * g1 + f * gamma_mono

            g1 = max(1.001, min(g1, gamma_mono))
            gamma1_table[i, j] = g1
            g1_min = min(g1_min, g1)
            g1_max = max(g1_max, g1)

    print(f"  Gamma1 range: [{g1_min:.4f}, {g1_max:.4f}]")
    return gamma1_table


def generate_eint_from_T_table_independent(n_grid, nh_grid, T_table, A_He,
                                            include_ion_e, blend=True):
    """Generate eint/nH(nH, T) table using independent forward Saha solve.

    Legacy method: solves Saha independently at each (nH, T) grid point.
    Not round-trip consistent with the forward T(nH, eint/nH) table because
    the forward and inverse tables are computed independently.

    Returns (table, T_grid) where:
      - table: log10(eint/nH) in SI (J), shape (n_nH, n_T)
      - T_grid: temperature grid in Kelvin

    If blend=True and include_ion_e, blend toward FI eint above T_blend.
    """
    n_per_nH_FI, neOnH_FI, eion_per_nH = fi_reference_values(A_He)
    do_blend = blend and include_ion_e

    T_min = np.min(T_table)
    T_max = np.max(T_table)

    T_grid = np.logspace(np.log10(T_min), np.log10(T_max), n_grid)
    eint_from_T = np.zeros((len(nh_grid), n_grid))

    for i, nH in enumerate(tqdm(nh_grid, desc="eint_from_T (independent)")):
        for j, T in enumerate(T_grid):
            eint_saha = solve_at_T(nH, T, A_He, include_ion_e)

            if do_blend:
                f = fi_blend_fraction(T)
                if f > 0:
                    eint_FI = 1.5 * n_per_nH_FI * kB * T + eion_per_nH
                    eint_blended = np.exp(np.log(max(eint_saha, 1e-100))*(1-f)
                                         + np.log(eint_FI)*f)
                    eint_from_T[i, j] = np.log10(max(eint_blended, 1e-100))
                else:
                    eint_from_T[i, j] = np.log10(max(eint_saha, 1e-100))
            else:
                eint_from_T[i, j] = np.log10(max(eint_saha, 1e-100))

    print(f"  eint_from_T via independent Saha solve")
    print(f"  log10(T/K) range: [{np.log10(T_min):.4f}, {np.log10(T_max):.4f}]")
    print(f"  log10(eint/nH) range: [{np.min(eint_from_T):.4f}, {np.max(eint_from_T):.4f}]")
    return eint_from_T, T_grid


def generate_eint_from_T_table(n_grid, nh_grid, T_table, eint_nh_grid,
                                blend=True):
    """Generate eint/nH(nH, T) table by PCHIP inversion of the forward T table.

    Numerically inverts the forward T(nH, eint/nH) table per nH slice using
    PCHIP interpolation.  This ensures the round-trip eint -> T -> eint is
    exact at forward grid nodes and O(h^4) between them.

    The forward T table already has FI blending applied, so the derived
    eint_from_T values automatically inherit the blend.

    Returns (table, T_grid) where:
      - table: log10(eint/nH) in SI (J), shape (n_nH, n_T)
      - T_grid: temperature grid in Kelvin
    """
    from scipy.interpolate import PchipInterpolator

    log_eint_nh = np.log10(eint_nh_grid)

    # Global T range from the forward table
    T_min = np.min(T_table)
    T_max = np.max(T_table)
    T_grid = np.logspace(np.log10(T_min), np.log10(T_max), n_grid)
    log_T_grid = np.log10(T_grid)

    eint_from_T = np.zeros((len(nh_grid), n_grid))
    n_monotone_fix = 0
    n_clamp = 0

    for i in tqdm(range(len(nh_grid)), desc="eint_from_T (derived)"):
        # Forward T column: log10(T) at each eint/nH grid point for this nH
        log_T_fwd = np.log10(np.maximum(T_table[i, :], 1e-100))

        # Ensure strict monotonicity (T must increase with eint/nH)
        for j in range(1, len(log_T_fwd)):
            if log_T_fwd[j] <= log_T_fwd[j - 1]:
                log_T_fwd[j] = log_T_fwd[j - 1] + 1e-12
                n_monotone_fix += 1

        # Build inverse: log10(T) -> log10(eint/nH) via PCHIP
        inv_interp = PchipInterpolator(log_T_fwd, log_eint_nh)

        # T range covered by the forward table for this nH slice
        T_lo, T_hi = log_T_fwd[0], log_T_fwd[-1]

        for j in range(n_grid):
            if T_lo <= log_T_grid[j] <= T_hi:
                eint_from_T[i, j] = float(inv_interp(log_T_grid[j]))
            elif log_T_grid[j] < T_lo:
                eint_from_T[i, j] = log_eint_nh[0]
                n_clamp += 1
            else:
                eint_from_T[i, j] = log_eint_nh[-1]
                n_clamp += 1

    print(f"  eint_from_T derived from forward T table (round-trip consistent)")
    print(f"  log10(T/K) range: [{np.log10(T_min):.4f}, {np.log10(T_max):.4f}]")
    print(f"  log10(eint/nH) range: [{np.min(eint_from_T):.4f}, "
          f"{np.max(eint_from_T):.4f}]")
    if n_monotone_fix > 0:
        print(f"  Monotonicity fixes: {n_monotone_fix} intervals")
    if n_clamp > 0:
        print(f"  Clamped: {n_clamp}/{len(nh_grid) * n_grid} points "
              f"({100 * n_clamp / (len(nh_grid) * n_grid):.1f}%)")
    return eint_from_T, T_grid



def write_bin_table(filepath, data, nH_bounds_log10, var2_bounds_log10):
    """Write a single table in AMRVAC binary format.

    Format:
      [2 x int32]   : (dim1, dim2) = data.shape
      [4 x float64] : (nH_min, nH_max, var2_min, var2_max) in log10
      [dim1 x dim2 x float64] : data in Fortran column-major order

    data is stored transposed (.T) so that numpy C-order write produces
    Fortran column-major layout matching read(10) quantity_in(dimy, dimx).
    """
    with open(filepath, "wb") as f:
        f.write(struct.pack("2i", data.shape[0], data.shape[1]))
        f.write(struct.pack("4d", *nH_bounds_log10, *var2_bounds_log10))
        data.T.tofile(f)


def save_tables_bin(outdir, comp_str, y_table, T_table, eint_pres_table,
                    nh_grid, eint_nh_grid, pres_nh_grid, A_He, units,
                    gamma1_table=None, eint_from_T_table=None, T_grid=None):
    """Write all binary table files for one variant.

    Args:
        outdir: output directory
        comp_str: composition string for filenames, e.g. "HHe_IonE"
        y_table: ne/nH table (n_nH, n_eint)
        T_table: temperature table (n_nH, n_eint) in K
        eint_pres_table: eint/p table (n_nH, n_pres) or None
        nh_grid: nH grid in SI (m^-3)
        eint_nh_grid: eint/nH grid in SI (J)
        pres_nh_grid: pres/nH grid in SI (Pa m^3)
        A_He: helium abundance (0.0 or 0.1)
        units: "cgs" or "si"
        gamma1_table: Gamma1 table (n_nH, n_eint), dimensionless, or None
        eint_from_T_table: log10(eint/nH) table (n_nH, n_T) in SI (J), or None
        T_grid: temperature grid in Kelvin for eint_from_T, or None
    """
    
    if units == "cgs":
        nH_offset = LOG10_NH_SI_TO_CGS      # -6: m^-3 -> cm^-3
        energy_offset = LOG10_ENERGY_SI_TO_CGS  # +7: J -> erg (also Pa m^3 -> dyn cm)
    else:
        nH_offset = 0.0
        energy_offset = 0.0

    # Bounds in log10 with unit conversion
    nH_bounds = (
        np.log10(nh_grid[0]) + nH_offset,
        np.log10(nh_grid[-1]) + nH_offset,
    )
    eint_bounds = (
        np.log10(eint_nh_grid[0]) + energy_offset,
        np.log10(eint_nh_grid[-1]) + energy_offset,
    )

    prefix = f"{outdir}/LTEeos_"

    # 1. ne/nH table (dimensionless — no unit conversion on values)
    path = f"{prefix}neOnH_{comp_str}.bin"
    write_bin_table(path, y_table, nH_bounds, eint_bounds)
    print(f"  -> {path}")

    # 2. T table (K — no unit conversion on values)
    path = f"{prefix}T_{comp_str}.bin"
    write_bin_table(path, T_table, nH_bounds, eint_bounds)
    print(f"  -> {path}")

    # 3. P table: log10(p/nH) derived from y and T
    #    p/nH = (1 + A_He + ne/nH) * kB * T  [SI: J, CGS: erg]
    P_table = np.log10((1.0 + A_He + y_table) * kB * T_table) + energy_offset
    path = f"{prefix}P_{comp_str}.bin"
    write_bin_table(path, P_table, nH_bounds, eint_bounds)
    print(f"  -> {path}")

    # 4. p2eint table (dimensionless ratio — no unit conversion on values)
    if eint_pres_table is not None:
        pres_bounds = (
            np.log10(pres_nh_grid[0]) + energy_offset,
            np.log10(pres_nh_grid[-1]) + energy_offset,
        )
        path = f"{prefix}p2eint_{comp_str}.bin"
        write_bin_table(path, eint_pres_table, nH_bounds, pres_bounds)
        print(f"  -> {path}")

    # 5. Gamma1 table (dimensionless — same (nH, eint/nH) grid as T/neOnH)
    if gamma1_table is not None:
        path = f"{prefix}gamma1_{comp_str}.bin"
        write_bin_table(path, gamma1_table, nH_bounds, eint_bounds)
        print(f"  -> {path}")

    # 6. eint_from_T table: log10(eint/nH) on (nH, T) grid
    if eint_from_T_table is not None and T_grid is not None:
        T_bounds = (
            np.log10(T_grid[0]),    # log10(T_min / K), no CGS/SI offset
            np.log10(T_grid[-1]),   # log10(T_max / K)
        )
        # Values are log10(eint/nH) in SI; apply energy offset for CGS
        eint_from_T_out = eint_from_T_table + energy_offset
        path = f"{prefix}eint_from_T_{comp_str}.bin"
        write_bin_table(path, eint_from_T_out, nH_bounds, T_bounds)
        print(f"  -> {path}")


def verify_round_trip(y_table, T_table, eint_pres_table, nh_grid,
                      eint_nh_grid, pres_nh_grid, A_He, name):
    """Verify round-trip consistency: eint -> p -> eint' and report max error.

    Simulates the AMRVAC round-trip:
      Forward: eint -> (T, y) via PCHIP on forward tables -> p = nH*(1+He+y)*kB*T
      Backward: p -> eint' = p * p2eint via PCHIP on p2eint table

    Reports the relative error |eint' - eint| / eint at random test points.
    """
    from scipy.interpolate import PchipInterpolator

    n_test = 500  # test points per nH slice (total = n_nH * n_test)
    n_nH = len(nh_grid)
    log_eint_nh = np.log10(eint_nh_grid)
    log_pres_nh = np.log10(pres_nh_grid)
    n_pres = len(pres_nh_grid)

    all_errors = []
    all_T = []

    rng = np.random.default_rng(42)

    for i in range(n_nH):
        nH = nh_grid[i]
        log_nH = np.log10(nH)

        # Forward table data for this nH slice
        T_row = T_table[i, :]      # T(eint/nH)
        y_row = y_table[i, :]      # ne/nH(eint/nH)

        # Build PCHIP interpolators for forward direction (mimics Fortran)
        T_interp = PchipInterpolator(log_eint_nh, np.log10(T_row))
        y_interp = PchipInterpolator(log_eint_nh, np.log10(y_row))

        # p2eint table data for this nH slice
        p2eint_row = eint_pres_table[i, :]

        # Build PCHIP interpolator for backward direction
        p2eint_interp = PchipInterpolator(log_pres_nh, p2eint_row)

        # Random test points in the eint/nH range
        log_eint_test = rng.uniform(log_eint_nh[1], log_eint_nh[-2], n_test)

        for log_e in log_eint_test:
            eint_nH = 10.0**log_e

            # Forward: eint -> T, y -> p
            log_T = float(T_interp(log_e))
            log_y = float(y_interp(log_e))
            T_val = 10.0**log_T
            y_val = 10.0**log_y
            p_nH = (1.0 + A_He + y_val) * kB * T_val
            log_p_nH = np.log10(p_nH)

            # Check p/nH is within the p2eint table range
            if log_p_nH < log_pres_nh[0] or log_p_nH > log_pres_nh[-1]:
                continue

            # Backward: p -> eint' = p * p2eint
            p2eint_val = float(p2eint_interp(log_p_nH))
            eint_nH_prime = p_nH * p2eint_val

            # Relative error
            rel_err = (eint_nH_prime - eint_nH) / eint_nH
            all_errors.append(rel_err)
            all_T.append(T_val)

    errors = np.array(all_errors)
    T_arr = np.array(all_T)

    print(f"\n  Round-trip verification ({name}):")
    print(f"    {len(errors)} test points")
    print(f"    max |error|:    {np.max(np.abs(errors)):.3e}")
    print(f"    mean error:     {np.mean(errors):+.3e}")
    print(f"    median |error|: {np.median(np.abs(errors)):.3e}")
    print(f"    95th pct:       {np.percentile(np.abs(errors), 95):.3e}")
    print(f"    99th pct:       {np.percentile(np.abs(errors), 99):.3e}")

    # Breakdown by temperature range
    T_bins = [(6e3, 50e3, "ionisation zone (6-50 kK)"),
              (50e3, 200e3, "TR (50-200 kK)"),
              (200e3, 1e7, "corona (>200 kK)")]
    for T_lo, T_hi, label in T_bins:
        mask = (T_arr >= T_lo) & (T_arr < T_hi)
        if np.any(mask):
            print(f"    {label}: max |e|={np.max(np.abs(errors[mask])):.3e}, "
                  f"mean={np.mean(errors[mask]):+.3e} ({np.sum(mask)} pts)")


def verify_T_round_trip(T_table, eint_from_T_table, nh_grid, eint_nh_grid,
                        T_grid, name):
    """Verify round-trip consistency: eint -> T -> eint' and report max error.

    Simulates the AMRVAC T-based prolongation round-trip:
      Forward: eint -> T via PCHIP on forward T table
      Backward: T -> eint' via PCHIP on eint_from_T table

    Reports the relative error |eint' - eint| / eint at random test points.
    """
    from scipy.interpolate import PchipInterpolator

    n_test = 500
    n_nH = len(nh_grid)
    log_eint_nh = np.log10(eint_nh_grid)
    log_T_grid = np.log10(T_grid)

    all_errors = []
    all_T = []

    rng = np.random.default_rng(42)

    for i in range(n_nH):
        # Forward table: T(eint/nH) for this nH slice
        log_T_fwd = np.log10(np.maximum(T_table[i, :], 1e-100))

        # Build PCHIP: log10(eint/nH) -> log10(T) (forward)
        T_interp = PchipInterpolator(log_eint_nh, log_T_fwd)

        # Inverse table: eint/nH(T) for this nH slice
        eint_row = eint_from_T_table[i, :]  # log10(eint/nH) in SI

        # Build PCHIP: log10(T) -> log10(eint/nH) (backward)
        inv_interp = PchipInterpolator(log_T_grid, eint_row)

        # Random test points in the eint/nH range
        log_eint_test = rng.uniform(log_eint_nh[1], log_eint_nh[-2], n_test)

        for log_e in log_eint_test:
            eint_nH = 10.0**log_e

            # Forward: eint -> T
            log_T = float(T_interp(log_e))
            T_val = 10.0**log_T

            # Check T is within the inverse table range
            if log_T < log_T_grid[0] or log_T > log_T_grid[-1]:
                continue

            # Backward: T -> eint'
            log_eint_prime = float(inv_interp(log_T))
            eint_nH_prime = 10.0**log_eint_prime

            rel_err = (eint_nH_prime - eint_nH) / eint_nH
            all_errors.append(rel_err)
            all_T.append(T_val)

    errors = np.array(all_errors)
    T_arr = np.array(all_T)

    print(f"\n  T-pathway round-trip verification ({name}):")
    print(f"    {len(errors)} test points")
    print(f"    max |error|:    {np.max(np.abs(errors)):.3e}")
    print(f"    mean error:     {np.mean(errors):+.3e}")
    print(f"    median |error|: {np.median(np.abs(errors)):.3e}")
    print(f"    95th pct:       {np.percentile(np.abs(errors), 95):.3e}")
    print(f"    99th pct:       {np.percentile(np.abs(errors), 99):.3e}")

    T_bins = [(6e3, 50e3, "ionisation zone (6-50 kK)"),
              (50e3, 200e3, "TR (50-200 kK)"),
              (200e3, 1e7, "corona (>200 kK)")]
    for T_lo, T_hi, label in T_bins:
        mask = (T_arr >= T_lo) & (T_arr < T_hi)
        if np.any(mask):
            print(f"    {label}: max |e|={np.max(np.abs(errors[mask])):.3e}, "
                  f"mean={np.mean(errors[mask]):+.3e} ({np.sum(mask)} pts)")


def plot_summary(results, units, outdir, filename="lte_eos_tables_summary.png"):

    if units == "cgs":
        nH_off = LOG10_NH_SI_TO_CGS
        e_off = LOG10_ENERGY_SI_TO_CGS
        nH_unit = r"cm$^{-3}$"
        e_unit = "erg"
        p_unit = "dyn cm"
    else:
        nH_off = 0.0
        e_off = 0.0
        nH_unit = r"m$^{-3}$"
        e_unit = "J"
        p_unit = r"Pa m$^3$"

    nH_lbl = rf"$\log_{{10}}(n_\mathrm{{H}}$ / {nH_unit})"
    eint_lbl = rf"$\log_{{10}}(e_\mathrm{{int}}/n_\mathrm{{H}}$ / {e_unit})"
    pres_lbl = rf"$\log_{{10}}(p/n_\mathrm{{H}}$ / {p_unit})"
    T_lbl = r"$\log_{10}(T$ / K)"

    # 6 panels per variant, laid out as 2 rows x 3 cols per variant
    panels_per_variant = 6
    n_rows = 2 * len(results)
    n_cols = 3
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(5 * n_cols, 4 * n_rows),
                              squeeze=False)

    for vi, res in enumerate(results):
        name = res['name']
        nh = res['nh_grid']
        eint_nh = res['eint_nh_grid']
        pres_nh = res['pres_nh_grid']
        y_tab = res['y_table']
        T_tab = res['T_table']
        A_He = res['A_He']

        log_nH = np.log10(nh) + nH_off
        log_eint = np.log10(eint_nh) + e_off
        log_pres = np.log10(pres_nh) + e_off

        r0 = 2 * vi      # first row for this variant
        r1 = 2 * vi + 1  # second row

        # Row 0, Col 0: ne/nH
        ax = axes[r0, 0]
        im = ax.pcolormesh(log_eint, log_nH,
                           np.log10(np.clip(y_tab, 1e-30, None)),
                           shading='auto', cmap='viridis')
        ax.set_xlabel(eint_lbl, fontsize=8)
        ax.set_ylabel(nH_lbl, fontsize=8)
        ax.set_title(f"{name}\nlog$_{{10}}$(ne/nH)", fontsize=9)
        ax.tick_params(labelsize=7)
        fig.colorbar(im, ax=ax)

        # Row 0, Col 1: T
        ax = axes[r0, 1]
        im = ax.pcolormesh(log_eint, log_nH, np.log10(T_tab),
                           shading='auto', cmap='inferno')
        ax.set_xlabel(eint_lbl, fontsize=8)
        ax.set_ylabel(nH_lbl, fontsize=8)
        ax.set_title("log$_{10}$(T / K)", fontsize=9)
        ax.tick_params(labelsize=7)
        fig.colorbar(im, ax=ax)

        # Row 0, Col 2: P = log10(p/nH)
        ax = axes[r0, 2]
        P_tab = np.log10((1.0 + A_He + y_tab) * kB * T_tab) + e_off
        im = ax.pcolormesh(log_eint, log_nH, P_tab,
                           shading='auto', cmap='viridis')
        ax.set_xlabel(eint_lbl, fontsize=8)
        ax.set_ylabel(nH_lbl, fontsize=8)
        ax.set_title(f"log$_{{10}}$(p/nH / {e_unit})", fontsize=9)
        ax.tick_params(labelsize=7)
        fig.colorbar(im, ax=ax)

        # Row 1, Col 0: p2eint (IonE only)
        ax = axes[r1, 0]
        if res.get('eint_pres_table') is not None:
            im = ax.pcolormesh(log_pres, log_nH, res['eint_pres_table'],
                               shading='auto', cmap='viridis')
            ax.set_xlabel(pres_lbl, fontsize=8)
            ax.set_ylabel(nH_lbl, fontsize=8)
            ax.set_title("eint/p", fontsize=9)
            ax.tick_params(labelsize=7)
            fig.colorbar(im, ax=ax)
        else:
            ax.set_visible(False)

        # Row 1, Col 1: gamma1 (IonE only)
        ax = axes[r1, 1]
        if res.get('gamma1_table') is not None:
            im = ax.pcolormesh(log_eint, log_nH, res['gamma1_table'],
                               shading='auto', cmap='RdBu_r',
                               vmin=1.0, vmax=5.0/3.0 + 0.02)
            ax.set_xlabel(eint_lbl, fontsize=8)
            ax.set_ylabel(nH_lbl, fontsize=8)
            ax.set_title(r"$\Gamma_1$", fontsize=9)
            ax.tick_params(labelsize=7)
            fig.colorbar(im, ax=ax)
        else:
            ax.set_visible(False)

        # Row 1, Col 2: eint_from_T (IonE only)
        ax = axes[r1, 2]
        if (res.get('eint_from_T_table') is not None
                and res.get('T_grid') is not None):
            log_T = np.log10(res['T_grid'])
            eint_vals = res['eint_from_T_table'] + e_off
            im = ax.pcolormesh(log_T, log_nH, eint_vals,
                               shading='auto', cmap='viridis')
            ax.set_xlabel(T_lbl, fontsize=8)
            ax.set_ylabel(nH_lbl, fontsize=8)
            ax.set_title(f"log$_{{10}}$(eint/nH / {e_unit})", fontsize=9)
            ax.tick_params(labelsize=7)
            fig.colorbar(im, ax=ax)
        else:
            ax.set_visible(False)

    fig.tight_layout()
    outpath = f"{outdir}/{filename}"
    fig.savefig(outpath, dpi=400, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSummary plot: {outpath}")


# -- Main ----------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--tables",
        nargs="+",
        choices=["H_NoIonE", "H_IonE", "HHe_NoIonE", "HHe_IonE", "all"],
        default=["all"],
        help="Which table variants to generate (default: all)",
    )
    parser.add_argument(
        "--grid", type=int, default=N_GRID, help=f"Grid size (default: {N_GRID})"
    )
    parser.add_argument(
        "--units",
        choices=["cgs", "si"],
        default="cgs",
        help="Unit system for table bounds: cgs (default, for AMRVAC) or si",
    )
    parser.add_argument(
        "--no-plot", action="store_true",
        help="Skip generating summary plot",
    )
    parser.add_argument(
        "--no-blend", action="store_true",
        help="Disable FI blending at hot end (pure Saha tables)",
    )
    parser.add_argument(
        "--p2eint-method",
        choices=["derived", "independent"],
        default="derived",
        help="p2eint table method: 'derived' (default) inverts forward tables "
             "for round-trip consistency; 'independent' solves Saha on pressure "
             "grid separately (legacy, has ~0.9%% round-trip error in ionisation zone)",
    )
    parser.add_argument(
        "--eint-from-T-method",
        choices=["derived", "independent"],
        default="derived",
        help="eint_from_T table method: 'derived' (default) inverts forward T "
             "table via PCHIP for round-trip consistency; 'independent' solves "
             "Saha independently on the T grid (legacy)",
    )
    args = parser.parse_args()

    if "all" in args.tables:
        tables = ["H_NoIonE", "H_IonE", "HHe_NoIonE", "HHe_IonE"]
    else:
        tables = args.tables

    n = args.grid
    nh_grid = np.logspace(NH_BOUNDS[0], NH_BOUNDS[1], n)
    eint_nh_grid = np.logspace(EINT_NH_BOUNDS[0], EINT_NH_BOUNDS[1], n)
    pres_nh_grid = np.logspace(PRES_NH_BOUNDS[0], PRES_NH_BOUNDS[1], n)

    outdir = "."

    cache = {}
    all_results = []

    for table_name in tables:
        print(f"\n=== Generating {table_name} ({args.units.upper()}) ===")

        H_only = table_name.startswith("H_")
        has_IonE = table_name.endswith("IonE")
        A_He = 0.0 if H_only else 0.1

        # Set the global flag for ionization energy inclusion
        global INCLUDE_ION_E
        INCLUDE_ION_E = has_IonE

        do_blend = not args.no_blend
        cache_key = (H_only, has_IonE, do_blend)
        if cache_key not in cache:
            y_tab, T_tab = generate_eint_tables(n, nh_grid, eint_nh_grid,
                                                 H_only=H_only, blend=do_blend)
            cache[cache_key] = (y_tab, T_tab)
        else:
            y_tab, T_tab = cache[cache_key]

        eint_pres_tab = None
        gamma1_tab = None
        eint_from_T_tab = None
        T_grid_tab = None

        if has_IonE:
            if args.p2eint_method == "derived":
                eint_pres_tab = generate_pres_table(
                    y_tab, T_tab, n, nh_grid, eint_nh_grid, pres_nh_grid,
                    H_only=H_only, blend=do_blend)
            else:
                eint_pres_tab = generate_pres_table_independent(
                    n, nh_grid, pres_nh_grid,
                    H_only=H_only, blend=do_blend)

            # Gamma1 via JAX autodiff (exact derivatives via IFT)
            gamma1_tab = compute_gamma1(
                y_tab, T_tab, nh_grid, eint_nh_grid, A_He, H_only,
                blend=do_blend)

            # eint_from_T: PCHIP inversion of forward T table (default) or
            # independent Saha solve (legacy)
            if args.eint_from_T_method == "derived":
                eint_from_T_tab, T_grid_tab = generate_eint_from_T_table(
                    n, nh_grid, T_tab, eint_nh_grid, blend=do_blend)
            else:
                eint_from_T_tab, T_grid_tab = \
                    generate_eint_from_T_table_independent(
                        n, nh_grid, T_tab, A_He, has_IonE, blend=do_blend)

        # Round-trip verification for IonE tables
        if has_IonE and eint_pres_tab is not None:
            verify_round_trip(y_tab, T_tab, eint_pres_tab, nh_grid,
                              eint_nh_grid, pres_nh_grid, A_He, table_name)

        # T-pathway round-trip verification
        if has_IonE and eint_from_T_tab is not None:
            verify_T_round_trip(T_tab, eint_from_T_tab, nh_grid,
                                eint_nh_grid, T_grid_tab, table_name)

        save_tables_bin(
            outdir, table_name, y_tab, T_tab, eint_pres_tab,
            nh_grid, eint_nh_grid, pres_nh_grid, A_He, args.units,
            gamma1_table=gamma1_tab,
            eint_from_T_table=eint_from_T_tab,
            T_grid=T_grid_tab,
        )

        all_results.append({
            'name': table_name,
            'nh_grid': nh_grid,
            'eint_nh_grid': eint_nh_grid,
            'pres_nh_grid': pres_nh_grid,
            'y_table': y_tab,
            'T_table': T_tab,
            'A_He': A_He,
            'eint_pres_table': eint_pres_tab,
            'gamma1_table': gamma1_tab,
            'eint_from_T_table': eint_from_T_tab,
            'T_grid': T_grid_tab,
        })

    if not args.no_plot and all_results:
        h_results = [r for r in all_results if r['name'].startswith('H_')]
        hhe_results = [r for r in all_results if r['name'].startswith('HHe_')]
        if h_results:
            plot_summary(h_results, args.units, outdir, "lte_eos_tables_H.png")
        if hhe_results:
            plot_summary(hhe_results, args.units, outdir, "lte_eos_tables_HHe.png")

    print("\nDone.")


if __name__ == "__main__":
    main()
