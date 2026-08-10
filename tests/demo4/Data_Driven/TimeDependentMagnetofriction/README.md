# Time-dependent magnetofriction

This template uses the standalone `mf` physics module, not the legacy
magnetofriction loop embedded in MHD.  Its six evolved variables are three
artificial magnetofriction velocities and `B1:B3`; there is no density.

It restarts from either a PotentialField or legacy
MagnetofrictionalRelaxation snapshot. `usr_transform_w` imports only the three
magnetic components and initializes the artificial velocity to zero. It then
drives the lower boundary with the B-only sequence
`B_0001.dat`, `B_0002.dat`, ... .  The frame time is stored in each file;
`driving_time_scale` in `usr_list` controls observational seconds advanced per
simulated second and defaults to 12. No physical velocity boundary files are
read. The observed B is imposed in the innermost lower ghost layer, while
farther B ghost layers and the artificial mf velocity are extrapolated.
Only the two magnetic frames bracketing the current observation time are kept
in memory; the small array of frame timestamps is retained for lookup.

The Python `DataDriven.ipynb` workflow copies this template and writes
`time_dependent_mf.par` with portable project-relative paths.
