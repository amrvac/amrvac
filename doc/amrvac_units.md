# Units

Describes the
- `hd_physical_units` subroutine in `hd/mod_hd_phys.t`.
- `mhd_physical_units` subroutine in `mhd/mod_mhd_phys.t`.

Consider two physical quantities \f$ A \f$ and \f$ B \f$ related by \f$ c \f$, so that \f$A=cB\f$. Write \f$A\f$ in units of \f$A_0\f$ as \f$A=A_0\widetilde A\f$, with \f$A\f$ variable and dimensional, \f$A_0\f$ constant and dimensional, and \f$\widetilde A\f$ variable and dimensionless. If we choose \f$\widetilde A=\widetilde B\f$, we have \f$A_0=cB_0\f$. The user specifies initial/boundary conditions in \f$\widetilde X\f$-quantities or scaled \f$\widetilde X\f$-quantities, and MPI-AMRVAC solves \f$\widetilde X\f$-quantities.  
More generally of course \f$\widetilde X\f$-quantities are related, e.g. if \f$A=cBD^2\f$, then \f$\widetilde A=\widetilde B\widetilde D^2\f$ and \f$A_0=cB_0D_0^2\f$.

MPI-AMRVAC constructs a *dependent* set of units \f$\{X_0\}\f$ for the set of physical quantities \f$\{X\}\f$ based on a set of *independent but complete* units \f$\{Y_0\}\f$ the user provides, various combinations are possible. If a unit \f$X_0\f$ cannot be determined since the user does not provide a fully complete set of independent units \f$\{Y_0\}\f$, \f$X_0=1\f$ is used.  
Constructing a set of units isn't always necessary since the user may simply provide \f$\widetilde X\f$-data, but the user may actually want their data to represent \f$a\widetilde A\f$ with some scale factor \f$a\f$, or MPI-AMRVAC may need to use data from tables which must be converted appropriately to the correct set of units, e.g. for optically thin radiative cooling from `physics/mod_radiative_cooling.t` used in e.g. `hd/mod_hd_phys.t`, specified in tables of \f$\log_{10}\text{ K}\f$. In such cases, the user should specify a complete set of units \f$\{Y_0\}\f$ in the `usr_init` subroutine of the `mod_usr.t` file.  
The constants \f$c\f$ are in general dimensional. Their dimensions are set either in SI units or Gaussian CGS units by MPI-AMRVAC (the user chooses), and these constants may depend also on additional dimensionless quantities the user provides, e.g. ionization levels.

MPI-AMRVAC's internal set of units for hydrodynamics/magnetohydrodynamics are:
- mass (\f$m_0\f$)
- length (\f$L_0\f$)
- time (\f$t_0\f$)
- velocity (\f$v_0\f$)
- number density (\f$n_0\f$)
- density (\f$\rho_0\f$)
- temperature (\f$T_0\f$)
- pressure (\f$P_0\f$)
- magnetic field (\f$B_0\f$)
- electric charge (\f$q_0\f$)

They are related by:
- \f$v_0=L_0/t_0\f$
- \f$m_0=\rho_0L_0^3\f$
- \f$P_0=\rho_0v_0^2\f$
- \f$\rho_0=am_{\text p} n_0\f$
- \f$P_0=bk_{\text B}n_0T_0\f$
- \f$P_0=B_0^2/\mu_0\f$
- \f$q_0=(B_0L_0^2)/(v_0\mu_0)\f$ (SI)
- \f$q_0=(B_0cL_0^2)/(v_0^2\mu_0)\f$ (Gaussian CGS, \f$c\f$ is speed of light)

Some information on the constants used.  
The proton mass `mp` is needed to relate number density and mass density, and the Boltzmann constant `kB` is needed to relate temperature, density, and pressure. They may be specified in either SI or Gaussian CGS units by MPI-AMRVAC.  
Helium abundance and ionization levels also affect these relationships.  
Electrons are assumed to be massless, neutrons are assumed to have proton mass, and binding energies (atomic and nuclear) are ignored.  
The constant \f$\mu_0\f$ in Gaussian CGS units is \f$1\f$, but electromagnetic equations in Gaussian CGS units are off by a factor of \f$4\pi\f$. MPI-AMRVAC defines \f$\mu_0\f$ to be \f$4\pi\f$ to ensure consistent conversions.  
The electric charge units are necessary for particle motion only.

Option 1: `eq_state_units`  
With `eq_state_units`, two molecular scaling factors are used, `a` for mass and `b` for the 'degrees of freedom per mass'.  
MPI-AMRVAC understands number density units in this case as **H number density units** (ionized or not). So `a` is the ratio of fluid mass to H component mass, `a = 1+4*He_abundance`. And this is used as `unit_density = a*mp*unit_numberdensity`.  
So we have \f$\widetilde\rho=\widetilde n\f$ and \f$n_0=am_pn_0\f$.  
`b` similarly associates number density units to pressure units. It is the number of degrees of freedom that come with one H proton, the unit pressure per H proton per unit temperature, using the equipartition theorem. Each H proton contributes `H_ion_fr` electrons, as well as `He_abundance` He atoms, which each contribute one alpha particle as well as `He_ion_fr` 1st electrons and `He_ion_fr*He_ion_fr2` 2nd electrons. So `b = 1+H_ion_fr+He_abundance*(He_ion_fr*(He_ion_fr2+1)+1)` in general (if `hd_partial_ionization`). If fully ionized (not `hd_partial_ionization`) this simplifies to `b = 2+3*He_abundance`. Now `unit_pressure = b*kB*unit_numberdensity*unit_temperature`, so \f$P_0=(b/a)(k_B/m_p)\rho_0T_0\f$ and \f$\widetilde P=\widetilde n\widetilde T=\widetilde\rho\widetilde T\f$.

Option 2: not `eq_state_units`  
In this mode `a` and `b` are kept fixed to `1`, so that number density units count **baryon number** (any form of mass). Instead of the units being affected by ionization, ionization affects the \f$\widetilde X\f$-equations directly. `b/a` is replaced by `RR`, and \f$\widetilde P=\text{RR}\widetilde\rho\widetilde T\f$. `kB` and `mp` are still used to relate \f$\rho_0\f$ to \f$n_0\f$, and \f$P_0\f$ to \f$n_0\f$ and \f$T_0\f$.
