# Equation of state in MPI-AMRVAC

[TOC]

The HD, MHD and FFHD modules do not compute pressure, temperature or sound speed
from a fixed ideal-gas law. They ask a single equation-of-state object, `eos`,
which is chosen at run time in the `&eos_list` namelist. The default is the fully
ionised ideal gas that MPI-AMRVAC has always used, so existing setups run
unchanged. The other two models, LTE and PI, describe partially ionised gas, as
found in the solar photosphere, chromosphere and prominences. This page describes
the models, how to select them, and how the module is organised. Set the
[normalization units](amrvac_units.md) in `usr_init()` before activating the
physics module, since the tables are converted to code units with them.

# Physics {#eos-physics}

Three gas models are available through `eos_type`.

**FI, fully ionised.** The gas is an ideal gas with a constant ratio of specific
heats \f$\gamma\f$ and a constant mean molecular weight:
\f[
 p=(\gamma-1)\,e_{\rm int},\qquad T=\frac{p}{R\,\rho},
\f]
where \f$R\f$ is the constant ionisation factor of a fully ionised H and He
plasma. Nothing is cached, and this is the classical behaviour of the code.

**LTE, local thermodynamic equilibrium.** Hydrogen, and helium if
`He_abundance` is positive, ionise according to the Saha equations at the local
temperature and density. The pressure then depends on the electron fraction
\f$y=n_e/n_{\rm H}\f$,
\f[
 p=(1+A_{\rm He}+y)\,n_{\rm H}k_{\rm B}T,
\f]
with \f$A_{\rm He}\f$ the helium abundance by number. With `ionE=.true.` the
internal energy also holds the energy spent on ionisation, so that
\f$e_{\rm int}=\tfrac{3}{2}p+n_{\rm H}E_{\rm ion}\f$, where \f$E_{\rm ion}\f$ is the
ionisation energy per hydrogen nucleus at the local state. The temperature and
the adiabatic index \f$\Gamma_1\f$ then vary through the partially ionised layers,
and \f$\Gamma_1\f$ falls well below 5/3 where hydrogen is half ionised. Without
`ionE`, only the mean molecular weight varies.

**PI, partial approximations.** The gas keeps the ideal-gas closure of FI, but
its mean mass follows an ionisation table chosen with `pi_table`. The
chromosphere and flare tables depend on temperature alone, and with `ionE=.true.`
they also fold the hydrogen and helium ionisation energy into \f$e_{\rm int}\f$. The
prominence table depends on temperature and pressure, and is only available
without the ionisation energy.

Hot cells take a shortcut in LTE. Above a threshold in internal energy per unit
mass the gas is fully ionised to within the table accuracy, and the FI formulas
are used without a table lookup.

# Practical use {#eos-practical}

## Choosing the gas model {#eos-namelist}

The model and its options are set in `&eos_list`. Without this namelist the run
uses FI with \f$\gamma=5/3\f$.

name | type | default | description
---|---|---|---
eos_type | string | 'FI' | Gas model, 'FI', 'LTE' or 'PI'.
eos_method | string | 'state' | How LTE evaluates the equilibrium. <br>'state' interpolates tabulated quantities (T, p, ne/nH, eint/p) on a grid of \f$\log_{10}n_{\rm H}\f$ and \f$\log_{10}(e_{\rm int}/n_{\rm H})\f$. <br>'entropy' reconstructs every quantity from one tabulated entropy potential with bicubic Hermite interpolation, so that they stay thermodynamically consistent. <br>'analytic' solves the hydrogen Saha equation on the fly and keeps helium neutral. <br>'tables' is accepted as the former name of 'state'.
ionE | logical | F | Include the ionisation energy in the internal energy. Only meaningful for LTE and PI, and forced to false for FI.
He_abundance | double precision | 0.1 | Helium abundance by number relative to hydrogen. With zero the hydrogen-only tables are used.
gamma | double precision | 5/3 | Ratio of specific heats of FI, and of fully ionised cells in LTE and PI.
gamma1_method | string | 'exact' | How LTE obtains \f$\Gamma_1\f$. <br>'exact' differentiates the pressure table. <br>'effective' uses \f$1+p/e_{\rm int}\f$, which is cheaper but only approximate. <br>'constant' returns `gamma` everywhere.
p2eint_method | string | 'table' | How LTE converts pressure to internal energy, 'table' (fast) or 'bisect' (accurate, and forced by the well-balanced HD setup).
table_location | string | `$AMRVAC_DIR/src/tables/eos_tables/uniform256/` | Directory that holds the LTE tables.
pi_table | string | 'chromosphere' | PI ionisation table, 'chromosphere', 'flare' or 'prominence'.

An LTE run with ionisation energy and helium, for example, uses

    &eos_list
       eos_type = 'LTE'
       eos_method = 'entropy'
       ionE = .true.
       He_abundance = 0.1d0
       table_location = '/path/to/amrvac/src/tables/eos_tables/entropy128/'
    /

and a PI prominence run uses

    &eos_list
       eos_type = 'PI'
       pi_table = 'prominence'
    /

## Tables {#eos-tables}

The LTE tables ship with the code in `src/tables/eos_tables`, as binary files on
axes of \f$\log_{10}n_{\rm H}\f$ and \f$\log_{10}(e_{\rm int}/n_{\rm H})\f$ in cgs.
There are three sets.

- `uniform256` holds the tables of the 'state' method on a uniform grid, for
  hydrogen only and for hydrogen and helium, with and without the ionisation
  energy. It is the default `table_location`.
- `hybrid128` holds the same tables on a non-uniform grid, with nodes gathered
  where the ionisation changes quickly. They are read by the 'state' method in
  the same way.
- `entropy128` holds the tables of the 'entropy' method, which include the
  ionisation energy.

The tables are shifted to code units once, when the run starts, so the unit
system must be set before the physics module is activated. The scripts
`generate_state_tables.py` and `entropy128/generate_entropy_tables.py`
regenerate them.

Every rank reads the tables at start-up. On a cluster it can pay to copy them to
node-local storage and point `table_location` there.

## Cached temperature and electron density {#eos-cache}

LTE adds two variables to each block, the temperature `Te` and the electron
density `Ne`. They are refreshed on every grid, ghost cells included, before each
substep of the time integration, and in the ghost cells after every boundary
update, so that the fluxes and sources read a consistent state without repeating
the table lookups. Both are written to the data files with the other variables.
A user routine that needs the temperature should call `eos%get_Te`, which returns
the cached field in LTE and computes \f$p/(R\rho)\f$ otherwise.

## Coupling to other physics {#eos-coupling}

Radiative cooling and thermal conduction take the temperature and the electron
density from the EoS, so their rates follow the chosen model. Within one
conduction step the ionisation map is frozen by default (`tc_frozen_ion=.true.` in
`&tc_list`), which avoids a table lookup at every super-timestepping stage. AMR
prolongation interpolates the temperature, not the internal energy, when LTE is
active, so that a refined block does not inherit a spurious ionisation front.

# Numerics and implementation {#eos-implementation}

The EoS lives in `src/eos`. `mod_eos` is the module the rest of the code uses. It
reads `&eos_list`, and its two lifecycle routines are called from `amrvac.t`:
`eos_init` before `usr_init`, so that the object exists when the physics module
links to it, and `eos_finalise` after the physics initialisation, once the units
are known and the tables can be converted. The object itself, with its parameters
and procedure pointers, is defined in `mod_eos_container`.

Each gas model has its own module. `mod_eos_fi` holds the ideal-gas kernels,
`mod_eos_lte` the LTE kernels and the per-substep update of `Te` and `Ne`, and
`mod_eos_pi` the PI adapter to its ionisation tables in `mod_eos_pi_tables`. The
LTE methods are implemented in `mod_eos_lte_state`, `mod_eos_lte_entropy` and
`mod_eos_lte_saha`, with the table handling in `mod_eos_lte_tables` and the
interpolation in `mod_eos_interp`. The physics modules connect through their own
EoS files, `mod_hd_eos`, `mod_mhd_eos` and `mod_ffhd_eos`, which point the
`phys_` routines for pressure, sound speed, temperature and the conversion
between conserved and primitive variables at the chosen model.

For LTE and PI runs the time spent in the EoS appears as separate categories in
the timing summary at the end of a run, so that its cost can be compared with
that of the fluxes and sources.
