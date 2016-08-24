# Test stressedIsland

\test stressedIsland
\todo Does not compile:

    comm_lib.f:63.10:
    size_block=nxG1*nxG2*nw*size_double
          1
    Error: Symbol 'size_block' at (1) has no IMPLICIT type

# Setup instructions

Setup this test case with

    $AMRVAC_DIR/setup.pl -d=23 -phi=0 -z=0 -g=12,12 -p=ff -eos=default -nf=0 -ndust=0 -u=nul -arch=default

# Description

This setup explores the explosive formation of a current sheet from an X-point
configuration. We solve the equations of force-free magnetodynamics, valid for
highly magnetised plasma. The overall evolution takes around one lightcrossing
time. You can explore the impact of resistivity `(eqpar(kpar_))` on the
reconnection rate and the width of the current sheet. There is quite a bit of
run-time analysis coded up in `analysis.t`, e.g. it measures the reconnection
rate and gives maximal values for the current and \f$E^2/B^2\f$.

Along with Maxwells equations, we push test-particles, each advanced due to the
Lorentz force given by the electric and magnetic fields of the simulation. This
particle-treatment is activated by defining in definitions.h:

    #define PARTICLES
    #define PARTICLES_LORENTZ

and by including the corresponding subroutines:

    INCLUDE:amrvacmodules/handle_particles.t
    INCLUDE:amrvacmodules/integrate_particles.t

In Paraview, you can load particle *.csv data and visualize them on top of the
fluid data. Use the filter `TableToPoints` to convert the data first. Try to
identify particles that show interesting evolution and follow them in the
simulation. This is done by setting `follow(index) = .true.` in subroutine
`init_particles` and re-running the simulation.


