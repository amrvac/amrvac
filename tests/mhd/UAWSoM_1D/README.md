# UAWSoM one-dimensional regression tests

This directory exercises the optional four-equation UAWSoM extension in the
ordinary MHD module.  The four parameter combinations test:

1. propagation direction and periodic conservation without effective sources;
2. conservative `wAplus`/`wAminus` exchange from the Eq. 34 reflection term;
3. nonlinear Alfvén and kink damping, with lost wave energy recovered as heat
   because total energy is conserved; and
4. repair of a negative wave energy by the standard small-value mechanism.

Run all short cases with `make -f test.make`.  Regression logs use the common
MPI-AMRVAC tolerances (relative `1e-5`, absolute `1e-8`).
