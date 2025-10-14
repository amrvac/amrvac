
# Overview of directory files

The base directory of the source code contains a few example LTE line-opacity tables for astrophysical applications. The following tables are included by default:

- `Y02400`: tables for a hydrogen-rich atmosphere appropriate for wind outflows from main-sequence OB-type stars (hydrogen mass fraction $X=0.74$, helium mass fraction $Y=0.24$, metal mass fraction $Z=0.02$).

- `Y09800`: tables for a pure helium atmosphere appropriate for wind outflows from classical Wolf-Rayet stars (hydrogen mass fraction $X=0$, helium mass fraction $Y=0.98$, metal mass fraction $Z=0.02$).

All these tables run in temperature from $\log_{10}(T) = 4.1...4.7$ (by 20 points) and mass density $\log_{10}(\rho) = -8...-22$ (by 50 points). Within each directory a given filename corresponds to the CAK line-opacity parameter tabulated as a function of mass density and temperature. For legacy reasons these directories are named after the helium mass fraction of the tables contained within the directory.

# Generating new tables

The line opacity tables use the parameters ($\alpha$, $\bar{Q}$, $Q_0$) introduced by Gayley (1995) together with the Thomson (free-electron) opacity $\kappa_e$. These opacities are appropriate to model wind outflows using the line-driven wind formalism of Castor, Abbott, & Klein (CAK). New tabulations for different mass fractions and/or metalicities can be made by following the method outlined in [Poniatowski et al. (2021, A&A 667)](https://ui.adsabs.harvard.edu/abs/2022A%26A...667A.113P/abstract).

It is important to follow the file structure and naming convention of files, as in the example directories, in order to avoid problems with reading in the tables. Because four files need to be given for a given metalicity, corresponding to the three line opacity parameters and the free-electron opacity, the files are expected to be grouped in a directory. Additionally, the user can provide a custom directory, different from this AMRVAC directory, where such table files are stored. Important is that such custom table files also have a matching data format of 20 (columns) and 50 (rows) for temperature and density, respectively.

# How-to

The usage of the CAK line-opacity tables is akin to that of the OPAL_tables (see README of OPAL tables). Except that the initialisation call is made to the subroutine `init_cak_table` and access to line opacity values is achieved by calling the subroutine `set_cak_opacity`, possibly also pointing to a different directory with table files than default here.
