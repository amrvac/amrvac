# Small values options in MPI-AMRVAC

# MPI-AMRVAC and negative pressure (or density)

When doing full hydrodynamic (HD) or MHD simulations, the density on the mesh `w(ixM^S,rho_)` must obviously always be strictly positive, and the same is true for the total energy density `w(ixM^S,e_)` (or the internal energy density, if used). The corresponding primitive variables, density and pressure, must similarly be strictly positive, always [note: When HD is extended to include pressureless dust, the density for any of the added dust species can be locally equal to zero, so their density is positive, but not strictly positive. Such pressureless dust, by definition, has zero pressure also. In MHD, one can handle zero temperature (hence zero pressure) in a zero-beta simulation, by using the isothermal closure with zero isothermal sound speed. In that case the simulation has no energy variable at all: `mhd_energy=.false.`]. The following discussion explains the use of positivity fixes, which may need to be activated to stably compute full (M)HD evolutions over extended timescales, at very high resolutions.


## Default settings: fix nothing, check nothing, and stop when negative pressure error occurs.

For HD and MHD, the default behavior to handle negative or small density/pressure values is controlled by the settings


```fortran
&methodlist
  check_small_values = .true.
  fix_small_values = .false.
  small_values_method='error'
  small_pressure=0.0d0
  small_density=0.0d0
  small_temperature=0.0d0
  small_values_fix_iw(1:nw)=.true.
  small_values_daverage=1
  trace_small_values=.false.
/
```

These default settings are initialized in various places, e.g. `check_small_values=T` and `fix_small_values=F` are in `mod_global_parameters.t`, while `small_values_method` is set in the module `mod_small_values.t`. Some defaults are also set in `amrvacio/mod_input_output.t`. They can all be overwritten in the `&methodlist` namelist. Note that not all of these defaults are always relevant, or are not even independent of eachother: e.g. a positive `small_temperature` could be used in combination with `small_density` to deduce `small_pressure=small_density*small_temperature` (in dimensionless code units). Similarly, an internal energy density `small_e` value can be used in the physics module, deduced from these and the EOS in use.

This particular default setting means that we NEVER do any artificial fix of positivity (since `fix_small_values=F`), assuming everything works as physically expected, and we only do a check on getting a positive pressure deduced each time we use the `phys_get_pthermal` subroutine (i.e. hd_get_pthermal or mhd_get_pthermal). Pressure positivity (actually p<small_pressure) is ALWAYS checked in that subroutine (thanks to `check_small_values=T`), and if we encounter p<small_pressure, an error is thrown with some info on its location. The code is then halted and forced to save the last timestep available, for possible debug analysis.

It is to be noted that this default setting is usually sufficient for most practical simulations. It also helps signaling possible implementation errors in user-coded initial conditions: the idea is obviously that you initialize with physically correct conservative values in all grid cells, for which the corresponding primitive variables are all positive. The same must be true for the ghost cells. If the code signals a negative pressure error within the first timestep, or after a few tens of timesteps, first scrutinize your own initial and boundary conditions for mistakes. If they are physically meaningful, and the code can set several 1000 or so timesteps without triggering any error, you are probably doing fine, and may proceed to activate possible non-default settings as described next.

## Level one bootstrapping: artificially fix/replace small pressure, only when needed.

If the code eventually (after thousands of steps without problem) throws an error, you may continue the run with a simply changed setting like

```fortran
&methodlist
  check_small_values = .false.
  fix_small_values = .true.
  small_values_method='error'
  small_pressure=1.0d-14
  small_density=1.0d-12
/
```
This setting will replace the deduced pressure in `phys_get_pthermal` by the small but positive `small_pressure` whenever the deduced pressure would end up below it. This sometimes suffices to fix and continue a run without further issue. Keeping the `small_values_method='error'` in combination with the `fix_small_values=T` implies that now, in various places throughout the code, positivity checks will become activated (all using the `small_pressure` and `small_density` values), and again errors will be thrown to stop the code if needed. 

## Level two bootstrapping: activate a small-values treatment, and hope for the best.

If the previous setting did not yet fully solve your problem, and you are certain that no coding error is made, you may proceed to a higher level of bootstrapping, eg. by changing to

```fortran
&methodlist
  check_small_values = .false.
  fix_small_values = .true.
  small_values_method='replace'
  small_pressure=1.0d-14
  small_density=1.0d-12
/
```

This activates a (much) more aggressive fixing of negative pressure (density) instances, and in reality will replace pressure and density values with their input floor values `small_density` and `small_pressure`. Due to `fix_small_values=T`, this now happens in several locations, so not only in the `phys_get_pthermal` subroutine. E.g. it happens when converting conservative to primitive variables, at the end of the `phys_to_primitive` subroutine. What exactly is done to fix things is controlled in the `phys_handle_small_values` subroutine, which has an HD and an MHD version. If `small_values_method\='ignore'`, we check and flag erronous entries, to be 'corrected' or reported further on. The actually implemented checks are found in the `hd_check_w` or `mhd_check_w` instances, using `small_pressure` and `small_density`.

The `small_values_method` can be one of 'replace' or 'average', and the latter 'average' variant uses a kind of buffer zone (with a width controlled by `small_values_daverage=1`) to average info from surrounding cells that are physically ok, to replace the faulty cell. This averaging approach works best on primitive variables directly.

In the `replace` variant, different behaviors can be enforced by the logicals `small_values_fix_iw(1:nw)`: these are normally all set to T, meaning that a small density will also nullify momenta and set the pressure to `small_pressure` (hence internal energy density to `small_e`). You are sort of personally responsible for finding out which (if any) combination and fixing variant works best.

# Debug options

For those familiar with using debuggers (like gdb), you may want to force the filling of corner (in 2D or 3D) ghost cells when trying to find out when and where a bug or problem happens. This is because in many instances we do not need the corner cells at all, and this means in practice that `phys_req_diagonal=.false.` unless the stencil (from specific added source terms) requires it. We added the `hd_force_diagonal` switch to the `&hd_list` for this purpose only: you can set it to T and then corners will be filled even though the stencil does not need it, but then no divisions by zero density values in corners will be reported, allowing to find out actual errors. You will get info on where a problem in a HD run occurs by a setting like

```fortran
&methodlist
  check_small_values = .true.
  fix_small_values = .true.
  small_values_method='error'
  small_pressure=1.0d-14
  small_density=1.0d-12
/

&hd_list
  hd_force_diagonal=.true.
/
```
