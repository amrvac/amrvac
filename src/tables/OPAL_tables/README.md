
# Overview of directory files

The base directory of the source code contains the following text files by default:

- `GN93hz`: gives an overview of all currently available Rosseland-mean opacity tabulations from the OPAL project (126 tables in total). All tabulations assume solar abundances from [Grevesse & Noels (1993)](https://ui.adsabs.harvard.edu/abs/1993oee..conf...15G/abstract).
- `Y09800`: example OPAL table for a pure helium atmosphere appropriate to model a wind outflows from classical Wolf-Rayet stars (hydrogen mass fraction $X=0$, helium mass fraction $Y=0.98$, metal mass fraction $Z=0.02$).

# Generating new tables

Different Rosseland-mean opacity tabulations can be used in a radiation-hydrodynamic simulation. The procedure for using a different Rosseland-mean opacity tabulation than the provided `Y09800` table is straightforward. First look up in the overview file `GN93hz` a suitable table, then copy its content into a new file, and save it accordingly in the current directory.

For example, when opacities from Table 92 are required, create a file with the following header in the current directory:

```
TABLE # 92     $G&N'93 Solar$     X=0.9000 Y=0.1000 Z=0.0000 dXc=0.0000 dXo=0.0000

                                                  log R

logT  -8.0   -7.5   -7.0   -6.5   -6.0   -5.5   -5.0   -4.5   -4.0   -3.5   -3.0   -2.5   -2.0   -1.5   -1.0   -0.5    0.0    0.5    1.0

3.75 -0.700 -0.870 -1.069 -1.281 -1.498 -1.711 -1.911 -2.078 -2.179 -2.182 -2.092 -1.936 -1.746 -1.539 -1.321 -1.097 -0.866 -0.629 -0.389
3.80 -0.473 -0.519 -0.606 -0.735 -0.889 -1.056 -1.217 -1.346 -1.422 -1.419 -1.337 -1.203 -1.026 -0.839 -0.639 -0.431 -0.214  0.013  0.244
3.85 -0.446 -0.445 -0.443 -0.445 -0.466 -0.514 -0.575 -0.630 -0.659 -0.639 -0.573 -0.462 -0.323 -0.164  0.008  0.197  0.391  0.602  0.820
3.90 -0.444 -0.440 -0.426 -0.390 -0.319 -0.217 -0.112 -0.020  0.054  0.120  0.192  0.281  0.391  0.508  0.645  0.804  0.972  1.156  1.351
......
```

It is important to not change the original file structure in order to avoid problems with reading in the table. In essence this means the newly created file containing a table should count 70 lines with different $log_{10}(T)$ values, as is standard for OPAL tables, and 19 columns with different $log_{10}(\rho)$ values. Additionally, the user can provide a custom directory, different from this AMRVAC directory, where such table files are stored (in the same format).

# How-to

There are essentially two ways to use a (new) OPAL opacity table in a simulation.

1. When using the FLD module of the code, the filename containing the table can be passed as a string to the `fld_opal_table` variable within the `fld_list` namelist. By default the FLD module will use the already included `Y09800` table.

2. When requiring the tables inside `mod_usr.t` the user needs to import some subroutines from `mod_opal_opacity.t`. Firstly, the table needs to be initialised by calling the subroutine `init_opal_table`. This should be done once, for example, inside the subroutine `initglobaldata_usr`. Secondly, to apply the values from the table to a given problem (e.g. boundary condition, source term computation), include a call to the subroutine `set_opal_opacity`.
See `mod_opal_opacity.t` for what arguments these subroutines require in their calls. An example of how to call OPAL tables inside the `mod_usr.t` file is given by the Wolf-Rayet stellar wind test problems of the rhd module. (Note: these test problems also rely on the FLD module, but this is not a prerequisite to use the tables.)
