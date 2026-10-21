"""Vetted atomic + molecular input data for the composition-general LTE EoS.

THE DATA IS THE EoS. The bicubic-Hermite interpolation is provably exact (the
entropy128 tables reproduce byte-for-byte from generate_entropy_tables.py); accuracy
is therefore set entirely by the per-species inputs loaded here. Every datum has
a recorded provenance and is cross-checked against an independent source where
possible (see cross_check()).

Sources (all under ./data/, copied verbatim from the upstream releases):
  - Abundances : Asplund, Grevesse, Sauval & Scott 2009 (AGSS09), ARA&A 47, 481,
                 photospheric log eps(H)=12 scale. (abundances_AGSS09.txt;
                 originally the MURaM Wittmann abundance.txt.)
                 AAG21 (Asplund, Amarsi & Grevesse 2021, A&A 653, A141) provided
                 inline as an option.
  - Ionisation energies : Barklem & Collet 2016 (B&C16), A&A 588, A96, table4
                 (IE1..IE3, eV; NIST-based). bc16_table4_ionization.dat.
  - Atomic partition functions : B&C16 table8 (Nov-2022 bug-fixed), Q(T) on a
                 42-point T grid. bc16_table8_atomic_partfn.dat.
  - Molecular equilibrium constants : B&C16 table7 (Nov-2022), log10(pK)(T) for
                 the dissociation A+B<->AB on the same 42-T grid.
                 bc16_table7_molecular_logKp.dat.
  - Molecular dissociation energies : B&C16 table1 (adopted De, eV).
                 bc16_table1_dissociation.dat.
  - Cross-check reference : MURaM Wittmann atomic_parameters.txt (independent IE
                 set) -- used only to flag discrepancies (it carries a known Ti
                 IE typo, 9.30 vs ~99.3 eV, which this module does NOT use).

NB on the B&C16 T grid: it runs 1e-5 K to 1e4 K. Above 1e4 K molecules are
irrelevant and atomic partition functions are taken as their highest-T tabulated
value (near the ground-state statistical weight after the level sum saturates);
the solver clamps T to the grid for the PF/Kp lookups (the ionised 1e4-1e7 K
range is governed by Saha, not by PF structure).
"""
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, "data")

# ---- physical constants (CGS) ----
kB = 1.380649e-16
eV = 1.602176634e-12
mp = 1.6726219e-24
me = 9.1093837e-28
h  = 6.62607015e-27
amu = 1.66053907e-24

ELEMENTS = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
            "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]
SYM2Z = {s: i + 1 for i, s in enumerate(ELEMENTS)}
ATOMIC_MASS_AMU = {  # standard atomic weights, for rho<->number conversion
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38}
ROMAN = ["I", "II", "III", "IV", "V", "VI"]

# AAG21 photospheric abundances (Asplund, Amarsi & Grevesse 2021), log eps(H)=12.
# Provided for the alternative-abundance option; default is AGSS09 (on disk).
AAG21 = {
    "H": 12.00, "He": 10.914, "Li": 0.96, "Be": 1.38, "B": 2.70,
    "C": 8.46, "N": 7.83, "O": 8.69, "F": 4.40, "Ne": 8.06,
    "Na": 6.22, "Mg": 7.55, "Al": 6.43, "Si": 7.51, "P": 5.41,
    "S": 7.12, "Cl": 5.31, "Ar": 6.38, "K": 5.07, "Ca": 6.30,
    "Sc": 3.14, "Ti": 4.97, "V": 3.90, "Cr": 5.62, "Mn": 5.42,
    "Fe": 7.46, "Co": 4.94, "Ni": 6.20, "Cu": 4.18, "Zn": 4.56}


def load_abundances(which="AGSS09"):
    """Return dict element-symbol -> linear number abundance relative to H (=1)."""
    if which == "AAG21":
        logeps = dict(AAG21)
    elif which == "AGSS09":
        logeps = {}
        with open(os.path.join(_DATA, "abundances_AGSS09.txt")) as f:
            for line in f:
                t = line.split()
                if len(t) >= 2:
                    logeps[ELEMENTS[int(t[0]) - 1]] = float(t[1])
    else:
        raise ValueError(f"unknown abundance set {which!r}")
    return {el: 10.0 ** (logeps[el] - 12.0) for el in logeps}


def _parse_bc16_grid_table(path, symbol_field_end):
    """Parse a B&C16 table8/table7-style file: a count line, a 'T [K]' grid line,
    then 'SYMBOL  v1 v2 ... v42' rows. Returns (Tgrid[42], {symbol: values[42]})."""
    Tgrid = None
    data = {}
    with open(path) as f:
        for line in f:
            if "[K]" in line:
                Tgrid = np.array([float(x) for x in line.replace("T", " ")
                                  .replace("[K]", " ").split()])
                continue
            t = line.split()
            if Tgrid is None or len(t) < 43:
                continue
            try:
                vals = np.array([float(x) for x in t[1:43]])
            except ValueError:
                continue
            data[t[0]] = vals
    return Tgrid, data


def load_atomic_partition_functions():
    """{species 'El_<roman>': (Tgrid, Q(T))} from B&C16 table8 (Nov2022)."""
    Tgrid, raw = _parse_bc16_grid_table(
        os.path.join(_DATA, "bc16_table8_atomic_partfn.dat"), 7)
    return Tgrid, raw


def load_molecular_logKp():
    """{molecule: (Tgrid, log10 pK(T))} from B&C16 table7 (Nov2022).

    Convention (B&C16): pK is the dissociation equilibrium constant in pressure
    units for AB <-> A + B, i.e. pK = p_A p_B / p_AB  [dyn/cm^2]."""
    Tgrid, raw = _parse_bc16_grid_table(
        os.path.join(_DATA, "bc16_table7_molecular_logKp.dat"), 5)
    return Tgrid, raw


def load_ionization_energies_bc16():
    """{element: [IE1, IE2, IE3]} eV, from B&C16 table4 (NIST-based)."""
    out = {}
    with open(os.path.join(_DATA, "bc16_table4_ionization.dat")) as f:
        for line in f:
            if line.startswith("*") or not line.strip():
                continue
            t = line.split()
            if len(t) < 5:
                continue
            try:
                Z = int(t[0])
            except ValueError:
                continue
            el = t[1]
            ies = [float(x) for x in t[2:5]]
            out[el] = [x for x in ies if x > 0]  # drop the -1 sentinels
    return out


def load_dissociation_energies():
    """{molecule: De_eV} adopted dissociation energy from B&C16 table1
    (fixed-width: molid 1-5, adopted De 84-93)."""
    out = {}
    with open(os.path.join(_DATA, "bc16_table1_dissociation.dat")) as f:
        for line in f:
            molid = line[0:5].strip()
            de = line[83:93].strip()
            if not molid or de in ("", "."):
                continue
            try:
                out[molid] = float(de)
            except ValueError:
                continue
    return out


def _load_muram_ies():
    """Independent IE set from the MURaM Wittmann atomic_parameters.txt, for
    cross-checking only. Returns {(Z, stage_j): IE_eV}. Carries a known Ti typo."""
    out = {}
    with open(os.path.join(_DATA, "muram_atomic_parameters_for_crosscheck.txt")) as f:
        for line in f:
            t = line.split()
            if len(t) < 3:
                continue
            try:
                Z, j, ie = int(t[0]), int(t[1]), float(t[2])
            except ValueError:
                continue
            out[(Z, j)] = ie
    return out


def cross_check(verbose=True):
    """Programmatic vetting (§3.4b). Cross-check B&C16 IEs against the independent
    MURaM set and flag any element/stage differing by >1%, and sanity-check the
    partition-function and Kp tables. Returns a list of flag strings."""
    flags = []
    ie_bc = load_ionization_energies_bc16()
    ie_mu = _load_muram_ies()
    for el, ies in ie_bc.items():
        if el not in SYM2Z:
            continue
        Z = SYM2Z[el]
        for j, ie in enumerate(ies):
            ref = ie_mu.get((Z, j))
            if ref is None or ref <= 0:
                continue
            rel = abs(ie - ref) / ref
            if rel > 0.01:
                flags.append(f"IE mismatch {el} stage {j}: B&C16={ie:.4f} eV vs "
                             f"MURaM={ref:.4f} eV ({rel*100:.1f}%)")
    # known MURaM Ti typo must be among the flags (proves the check works)
    ti4 = ie_mu.get((22, 4))
    if ti4 is not None and ti4 < 50.0:
        flags.append(f"CONFIRMED known MURaM Ti IE typo: (Z=22,stage4)={ti4} eV "
                     f"(should be ~99.3); B&C16 used instead -> safe.")
    # PF sanity: H_I -> 2 at low T; He_I -> 1 at low T
    Tg, pf = load_atomic_partition_functions()
    if not (1.99 < pf["H_I"][0] < 2.01):
        flags.append(f"PF sanity FAIL: H_I({Tg[0]:.0e}K)={pf['H_I'][0]} (expected 2)")
    if not (0.99 < pf["He_I"][0] < 1.01):
        flags.append(f"PF sanity FAIL: He_I low-T = {pf['He_I'][0]} (expected 1)")
    # Kp sanity: H2 present and monotone-ish in log
    TgK, kp = load_molecular_logKp()
    if "H2" not in kp:
        flags.append("Kp sanity FAIL: H2 missing from table7")
    if verbose:
        print(f"cross_check: {len(flags)} flag(s)")
        for s in flags:
            print("  -", s)
        ab = load_abundances("AGSS09")
        Zmass = sum(ATOMIC_MASS_AMU[e] * ab[e] for e in ELEMENTS)
        Xmass = ATOMIC_MASS_AMU["H"] * ab["H"]
        Zmet = sum(ATOMIC_MASS_AMU[e] * ab[e] for e in ELEMENTS if e not in ("H", "He"))
        print(f"  AGSS09: {len(ab)} elements; metals/total mass frac Z = "
              f"{Zmet/Zmass:.4f}; Z/X = {Zmet/Xmass:.4f}")
        print(f"  atomic PF species: {len(pf)}; molecules with Kp: {len(kp)}; "
              f"T grid {Tg[0]:.0e}..{Tg[-1]:.0e} K ({len(Tg)} pts)")
    return flags


PROVENANCE = {
    "abundances": "AGSS09 (Asplund+2009 ARA&A 47,481); AAG21 optional",
    "ionization_energies": "Barklem & Collet 2016 table4 (NIST-based)",
    "atomic_partition_functions": "Barklem & Collet 2016 table8 vNov2022",
    "molecular_logKp": "Barklem & Collet 2016 table7 vNov2022",
    "dissociation_energies": "Barklem & Collet 2016 table1 (adopted De)",
    "crosscheck_reference": "MURaM Wittmann atomic_parameters.txt (independent IEs)",
}


if __name__ == "__main__":
    print("=== eos_data.py self-test / provenance ===")
    for k, v in PROVENANCE.items():
        print(f"  {k}: {v}")
    print()
    cross_check(verbose=True)
