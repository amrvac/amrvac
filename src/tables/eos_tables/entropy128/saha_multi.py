"""Multi-element + molecular LTE chemical-equilibrium solver (CGS).

Drop-in replacement for saha_np's eos_state / eos_at_rho_eint, generalised from
H+He to an arbitrary composition with metals and molecules, using the vetted
B&C16 / AGSS09 data in eos_data.py.

Scope (v1, matches MURaM Wittmann scope + B&C16 data):
  - Atomic ionisation: every element in the composition, all stages for which an
    ionisation energy is available (B&C16 table4 gives 3 stages), with
    T-dependent partition functions (B&C16 table8).
  - Molecules: H2 (+ optionally H2+, H-) -- the bulk-EoS molecules. Heavier
    diatomics (CO/OH/CN/TiO) are Tier-2 and added later via the same machinery.
  - Charge neutrality solved by Newton/bisection on n_e.

Returned thermodynamics use the same conventions as saha_np:
  p     = (sum of all heavy particles + n_e) kB T
  e_int = 1.5 (sum particles) kB T + ionisation energy stored + dissociation
          energy stored (relative to fully-dissociated neutral atoms)
  s     = sum of per-species Sackur-Tetrode translational entropy
"""
import numpy as np
import eos_data as D

kB, eV, mp, me, h, amu = D.kB, D.eV, D.mp, D.me, D.h, D.amu

# ---------- preload vetted data once ----------
_Tpf, _PF = D.load_atomic_partition_functions()      # {'El_I': Q[42]}, Tgrid
_TKp, _KP = D.load_molecular_logKp()                 # {'H2': log10pK[42]}, Tgrid
_IE = D.load_ionization_energies_bc16()              # {'El': [IE1,IE2,IE3]} eV
_lnTpf = np.log(_Tpf)
_lnTKp = np.log(_TKp)
_DE = D.load_dissociation_energies()                 # {molid: De_eV}

# Diatomic molecules carried (Tier-1 H2 + Tier-2 cool-end set). Each maps to its
# two constituent elements. Only those with a B&C16 Kp and constituents present
# in the composition are activated. (H2O etc. are triatomic -> not in B&C16.)
DEFAULT_MOLECULES = {
    "H2": ("H", "H"), "CO": ("C", "O"), "OH": ("O", "H"), "CN": ("C", "N"),
    "C2": ("C", "C"), "N2": ("N", "N"), "NO": ("N", "O"), "NH": ("N", "H"),
    "CH": ("C", "H"), "TiO": ("Ti", "O"), "SiO": ("Si", "O"), "SiH": ("Si", "H"),
    "MgH": ("Mg", "H"), "CaH": ("Ca", "H"), "AlH": ("Al", "H"), "FeH": ("Fe", "H"),
}


def _pf(species, T):
    """Atomic partition function Q(T), log-log interpolated, clamped to grid."""
    q = _PF.get(species)
    if q is None:
        return 1.0
    lt = np.log(min(max(T, _Tpf[0]), _Tpf[-1]))
    return float(np.exp(np.interp(lt, _lnTpf, np.log(q))))


def _logKp(mol, T):
    """log10 of the dissociation equilibrium constant pK(T) [dyn/cm^2]."""
    k = _KP.get(mol)
    if k is None:
        return None
    lt = np.log(min(max(T, _TKp[0]), _TKp[-1]))
    return float(np.interp(lt, _lnTKp, k))


def make_composition(abundance_set="AGSS09", elements=None, molecules=None):
    """Build a composition record. `molecules` is an iterable of molecule names
    (defaults to DEFAULT_MOLECULES); each is activated only if it has a B&C16 Kp
    and both constituents are in the composition. Each activated molecule is
    stored as (name, (elA, elB), De_eV)."""
    ab = D.load_abundances(abundance_set)
    if elements is None:
        elements = [e for e in D.ELEMENTS if ab.get(e, 0.0) > 0.0]
    elset = set(elements)
    if molecules is None:
        molecules = list(DEFAULT_MOLECULES.keys())
    mols = []
    for m in molecules:
        cons = DEFAULT_MOLECULES.get(m)
        if cons is None or _logKp(m, 5000.0) is None:
            continue
        if cons[0] in elset and cons[1] in elset:
            mols.append((m, cons, _DE.get(m, 0.0)))
    comp = {
        "set": abundance_set,
        "elements": elements,
        "A": {e: ab[e] for e in elements},                 # rel. to H
        "mass": {e: D.ATOMIC_MASS_AMU[e] * amu for e in elements},
        "IE": {e: _IE.get(e, []) for e in elements},       # eV per stage
        "molecules": mols,                                 # [(name,(elA,elB),De)]
    }
    comp["mass_per_H"] = sum(comp["A"][e] * comp["mass"][e] for e in elements)
    # Energy zero-point: reference e_int to the cold fully-molecular neutral ground
    # state so e_int > 0 everywhere a molecular table needs log10(eint). A pure
    # per-nucleus constant -> changes no physics (Gamma_1, latent heat, p, T all
    # untouched); only shifts the absolute zero. Value = the most-negative raw
    # eint/nH over the cold molecular regime (i.e. the full-molecularisation
    # binding per H), found by a quick cold scan, +5% margin.
    comp["eint_offset"] = 0.0
    _lr = np.linspace(-9.0, -5.0, 8); _T = np.array([300., 500., 800., 1200., 1800.])
    _RHO, _TT = np.meshgrid(10.0 ** _lr, _T, indexing="ij")
    _st = solve_grid(_RHO.ravel(), _TT.ravel(), comp)   # offset still 0 here
    _mn = float(np.min(_st["e_int"] / _st["nH"]))
    comp["eint_offset"] = max(0.0, -_mn) * 1.05
    return comp


def _saha_S(e, T, sahafac):
    """Saha ratios S_{j} = n_{j+1} n_e / n_j for element e (length = #stages)."""
    ies = _IE.get(e, [])
    S = []
    for j, chi in enumerate(ies):
        Uj = _pf(f"{e}_{D.ROMAN[j]}", T)
        Uj1 = _pf(f"{e}_{D.ROMAN[j+1]}", T) if j + 1 < len(D.ROMAN) else 1.0
        S.append(2.0 * (Uj1 / max(Uj, 1e-300)) * sahafac
                 * np.exp(-chi * eV / (kB * T)))
    return S


def _element_state(e, T, ne, sahafac):
    """Given n_e, return (frac[stages], mean_charge, neutral_fraction) for element e
    via Saha (atomic only; H molecules handled separately in solve)."""
    S = _saha_S(e, T, sahafac)
    phi = [1.0]
    for Sj in S:
        phi.append(phi[-1] * Sj / max(ne, 1e-300))
    tot = sum(phi)
    frac = [p / tot for p in phi]
    qbar = sum(j * frac[j] for j in range(len(frac)))
    return frac, qbar


def _ion_factors(S, ne):
    """phi[stages], Phi=sum phi, qfac=sum j*phi (electrons per neutral atom)
    from the ne-independent Saha ratios S."""
    phi = [1.0]
    for Sj in S:
        phi.append(phi[-1] * Sj / ne)
    Phi = 0.0; qfac = 0.0
    for j, p in enumerate(phi):
        Phi += p; qfac += j * p
    return phi, Phi, qfac


def _activities(elements, S_all, N, mols, Cm, ne):
    """Coupled neutral-atom activities a_X (= n_{X,I}) at fixed n_e, by Gauss-Seidel
    on element conservation  N_X = a_X*Phi_X + sum_m nu_{X,m} n_m  with diatomic
    molecules  n_m = C_m a_A a_B.  Returns (a, phi_all, qfac)."""
    phi_all = {}; Phi = {}; qfac = {}
    for e in elements:
        phi_all[e], Phi[e], qfac[e] = _ion_factors(S_all[e], ne)
    a = {e: N[e] / Phi[e] for e in elements}        # init: no molecules
    if mols:
        for _ in range(60):
            maxch = 0.0
            for e in elements:
                lin = Phi[e]; quad = 0.0
                for (name, (A, B), De), C in zip(mols, Cm):
                    if C <= 0.0:
                        continue
                    if A == e and B == e:           # homonuclear: 2 of e per mol
                        quad += 2.0 * C
                    elif A == e:
                        lin += C * a[B]
                    elif B == e:
                        lin += C * a[A]
                ratio = (4.0 * quad / lin) * (N[e] / lin)   # = 4 q N / lin^2 (overflow-safe)
                anew = 2.0 * N[e] / (lin * (1.0 + (1.0 + ratio) ** 0.5))
                if anew > 0.0:
                    maxch = max(maxch, abs(anew - a[e]) / anew)
                a[e] = anew
            if maxch < 1e-11:
                break
    return a, phi_all, qfac


def solve(rho, T, comp):
    """Solve chemical equilibrium at (rho, T). Returns number densities + summary.

    Molecular thermodynamics: the B&C16 Kp already carries the full internal
    partition functions, so the molecular *populations* are exact. Internal
    (rotational/vibrational) energy and entropy of molecules are not yet added
    (Tier-2 refinement; matches the MURaM Wittmann port which zeroed them) -- the
    dominant molecular EoS effect, the dissociation latent heat, IS included."""
    nH = rho / comp["mass_per_H"]
    sahafac = (2.0 * np.pi * me * kB * T / h**2) ** 1.5
    elements = comp["elements"]
    N = {e: comp["A"][e] * nH for e in elements}
    S_all = {e: _saha_S(e, T, sahafac) for e in elements}      # ne-independent
    # molecules negligible above the B&C16 Kp grid max (~1e4 K); disable (do not
    # clamp Kp there: clamping while kT uses the true T fabricates spurious mols).
    mols = comp["molecules"] if T <= _TKp[-1] else []
    Cm = [kB * T / (10.0 ** _logKp(name, T)) for (name, c, De) in mols]

    # n_e by bisection in log n_e; upper bound = fully-ionised electron budget
    zmax = sum(comp["A"][e] * max(len(S_all[e]), 1) for e in elements)
    lo, hi = np.log(nH * 1e-18), np.log(nH * (zmax + 1.0))
    for _ in range(64):
        ne = np.exp(0.5 * (lo + hi))
        a, phi_all, qfac = _activities(elements, S_all, N, mols, Cm, ne)
        ne_from = sum(a[e] * qfac[e] for e in elements)
        if ne_from > ne:
            lo = 0.5 * (lo + hi)
        else:
            hi = 0.5 * (lo + hi)
    ne = np.exp(0.5 * (lo + hi))
    a, phi_all, qfac = _activities(elements, S_all, N, mols, Cm, ne)

    species = []; n_heavy = 0.0; e_ion = 0.0; e_diss = 0.0
    for e in elements:
        ies = _IE.get(e, []); cum = 0.0; phi = phi_all[e]
        for j in range(len(phi)):
            nj = a[e] * phi[j]
            n_heavy += nj
            species.append((f"{e}_{D.ROMAN[j]}", nj, comp["mass"][e], _pf(f"{e}_{D.ROMAN[j]}", T)))
            if j >= 1:
                cum += ies[j - 1] * eV
                e_ion += nj * cum
    nmol = {}
    for (name, (A, B), De), C in zip(mols, Cm):
        nm = C * a[A] * a[B]
        nmol[name] = nm
        n_heavy += nm
        species.append((name, nm, comp["mass"][A] + comp["mass"][B], 1.0))
        e_diss += nm * (-De * eV)              # dissociation latent heat (bound = -De)
    n_tot = n_heavy + ne
    p = n_tot * kB * T
    e_int = 1.5 * n_tot * kB * T + e_ion + e_diss + comp.get('eint_offset', 0.0) * nH

    def stran(n, mass, g):
        if n <= 0.0:
            return 0.0
        logterm = (np.log(max(g, 1e-30))
                   + 1.5 * np.log(2.0 * np.pi * mass * kB * T / h**2) - np.log(n))
        return n * kB * (2.5 + logterm)
    s = stran(ne, me, 2.0)
    for name, n, mass, g in species:
        s += stran(n, mass, g)

    phiH = phi_all["H"]
    nHI = a["H"] * phiH[0]
    nHII = a["H"] * phiH[1] if len(phiH) > 1 else 0.0
    return dict(nH=nH, ne=ne, ne_over_nH=ne / nH, nHI=nHI, nHII=nHII,
                nH2=nmol.get("H2", 0.0), nmol=nmol, p=p, e_int=e_int, s=s,
                n_tot=n_tot, a=a)


# ===================================================================
# Vectorised solver: same physics as solve(), evaluated over a whole
# (rho, T) grid with numpy array ops (loops only over elements / molecules
# / sweeps). 100-1000x faster than per-point solve() for table generation.
# ===================================================================
def _pf_grid(species, Tarr):
    q = _PF.get(species)
    if q is None:
        return np.ones_like(Tarr)
    lt = np.log(np.clip(Tarr, _Tpf[0], _Tpf[-1]))
    return np.exp(np.interp(lt, _lnTpf, np.log(q)))


def _logKp_grid(mol, Tarr):
    lt = np.log(np.clip(Tarr, _TKp[0], _TKp[-1]))
    return np.interp(lt, _lnTKp, _KP[mol])


def _saha_S_grid(e, Tarr, sahafac):
    ies = _IE.get(e, [])
    S = []
    for j, chi in enumerate(ies):
        Uj = _pf_grid(f"{e}_{D.ROMAN[j]}", Tarr)
        Uj1 = _pf_grid(f"{e}_{D.ROMAN[j+1]}", Tarr) if j + 1 < len(D.ROMAN) else np.ones_like(Tarr)
        S.append(2.0 * (Uj1 / np.maximum(Uj, 1e-300)) * sahafac
                 * np.exp(-chi * eV / (kB * Tarr)))
    return S


def _activities_grid(elements, S_all, N, mols, Cm, ne):
    phi_all = {}; Phi = {}; qfac = {}
    for e in elements:
        phi = [np.ones_like(ne)]
        for Sj in S_all[e]:
            phi.append(phi[-1] * Sj / ne)
        Phi[e] = sum(phi)
        qfac[e] = sum(j * phi[j] for j in range(len(phi)))
        phi_all[e] = phi
    a = {e: N[e] / Phi[e] for e in elements}
    if mols:
      with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for _ in range(60):
            maxch = 0.0
            for e in elements:
                lin = Phi[e].copy(); quad = np.zeros_like(ne)
                for (name, (A, B), De), C in zip(mols, Cm):
                    if A == e and B == e:
                        quad = quad + 2.0 * C
                    elif A == e:
                        lin = lin + C * a[B]
                    elif B == e:
                        lin = lin + C * a[A]
                # stable positive root of quad*a^2 + lin*a - N = 0 (lin,quad,N >= 0;
                # lin = Phi >= 1 > 0). Avoids -lin+sqrt cancellation and lin^2 overflow;
                # reduces to N/lin when quad=0.
                ratio = (4.0 * quad / lin) * (N[e] / lin)   # = 4 q N / lin^2
                anew = 2.0 * N[e] / (lin * (1.0 + np.sqrt(1.0 + ratio)))
                maxch = max(maxch, float(np.max(np.abs(anew - a[e]) / np.maximum(anew, 1e-300))))
                a[e] = anew
            if maxch < 1e-11:
                break
    return a, phi_all, qfac


def solve_grid(rho, T, comp):
    """Vectorised chemical equilibrium over flat arrays rho[K], T[K]. Returns a
    dict of K-length arrays (p, e_int, s, ne_over_nH, nH2, n_tot, ...)."""
    rho = np.asarray(rho, float); T = np.asarray(T, float)
    nH = rho / comp["mass_per_H"]
    sahafac = (2.0 * np.pi * me * kB * T / h**2) ** 1.5
    elements = comp["elements"]
    N = {e: comp["A"][e] * nH for e in elements}
    S_all = {e: _saha_S_grid(e, T, sahafac) for e in elements}
    mols = comp["molecules"]
    mask = (T <= _TKp[-1])
    # C = kT/pK = kT*10^(-logKp); clip exponent so a deeply-bound molecule at very
    # low T (logKp -> large negative) gives a huge-but-finite C (-> all atoms bound)
    # rather than overflowing to inf.
    Cm = [np.where(mask, kB * T * 10.0 ** (-np.maximum(_logKp_grid(name, T), -290.0)), 0.0)
          for (name, c, De) in mols]

    zmax = sum(comp["A"][e] * max(len(S_all[e]), 1) for e in elements)
    lo = np.log(nH * 1e-18); hi = np.log(nH * (zmax + 1.0))
    for _ in range(64):
        ne = np.exp(0.5 * (lo + hi))
        a, phi_all, qfac = _activities_grid(elements, S_all, N, mols, Cm, ne)
        ne_from = sum(a[e] * qfac[e] for e in elements)
        up = ne_from > ne
        mid = 0.5 * (lo + hi)
        lo = np.where(up, mid, lo)
        hi = np.where(up, hi, mid)
    ne = np.exp(0.5 * (lo + hi))
    a, phi_all, qfac = _activities_grid(elements, S_all, N, mols, Cm, ne)

    n_heavy = np.zeros_like(nH); e_ion = np.zeros_like(nH); e_diss = np.zeros_like(nH)
    s = np.zeros_like(nH)

    def stran(n, mass, g):
        nn = np.maximum(n, 1e-300)
        val = nn * kB * (2.5 + np.log(np.maximum(g, 1e-30))
                         + 1.5 * np.log(2.0 * np.pi * mass * kB * T / h**2) - np.log(nn))
        return np.where(n > 0.0, val, 0.0)

    s = s + stran(ne, me, 2.0)
    for e in elements:
        ies = _IE.get(e, []); cum = np.zeros_like(nH); phi = phi_all[e]
        for j in range(len(phi)):
            nj = a[e] * phi[j]
            n_heavy = n_heavy + nj
            s = s + stran(nj, comp["mass"][e], _pf_grid(f"{e}_{D.ROMAN[j]}", T))
            if j >= 1:
                cum = cum + ies[j - 1] * eV
                e_ion = e_ion + nj * cum
    nmol = {}
    for (name, (A, B), De), C in zip(mols, Cm):
        nm = C * a[A] * a[B]
        nmol[name] = nm
        n_heavy = n_heavy + nm
        s = s + stran(nm, comp["mass"][A] + comp["mass"][B], 1.0)
        e_diss = e_diss + nm * (-De * eV)
    n_tot = n_heavy + ne
    p = n_tot * kB * T
    e_int = 1.5 * n_tot * kB * T + e_ion + e_diss + comp.get('eint_offset', 0.0) * nH
    phiH = phi_all["H"]
    nHI = a["H"] * phiH[0]
    nHII = a["H"] * phiH[1] if len(phiH) > 1 else np.zeros_like(nH)
    return dict(nH=nH, ne=ne, ne_over_nH=ne / nH, nHI=nHI, nHII=nHII,
                nH2=nmol.get("H2", np.zeros_like(nH)), nmol=nmol,
                p=p, e_int=e_int, s=s, n_tot=n_tot)


def eos_state(rho, T, comp):
    st = solve(rho, T, comp)
    return st["p"], st["e_int"], st["s"]


def eos_at_rho_eint(rho, e_int_target, comp):
    lo, hi = np.log(100.0), np.log(2e6)
    for _ in range(60):
        m = 0.5 * (lo + hi)
        _, e_int, _ = eos_state(rho, np.exp(m), comp)
        if e_int < e_int_target:
            lo = m
        else:
            hi = m
    T = np.exp(0.5 * (lo + hi))
    p, e_int, s = eos_state(rho, T, comp)
    return T, p, s


if __name__ == "__main__":
    comp = make_composition("AGSS09")   # full default molecule set
    print(f"composition: {len(comp['elements'])} elements; "
          f"molecules: {[m[0] for m in comp['molecules']]}")
    print(f"mass_per_H={comp['mass_per_H']:.4e} g")
    print("\ncool dense (rho=1e-7) -- molecular regime:")
    for T in [2500., 3000., 3500., 4000., 5000.]:
        st = solve(1e-7, T, comp)
        frac = {k: v / st['nH'] for k, v in st['nmol'].items()}
        top = sorted(frac.items(), key=lambda kv: -kv[1])[:4]
        print(f"  T={T:6.0f} ne/nH={st['ne_over_nH']:.2e}  top mol/nH: "
              + "  ".join(f"{k}={v:.2e}" for k, v in top))
    print("\nhot dilute (rho=1e-9) -- ionisation regime:")
    for T in [8000., 1e4, 3e4, 1e6]:
        st = solve(1e-9, T, comp)
        print(f"  T={T:8.0f} ne/nH={st['ne_over_nH']:.4f}")
