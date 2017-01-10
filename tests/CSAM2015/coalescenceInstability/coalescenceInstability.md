# Test coalescenceInstability

\test coalescenceInstability

# Setup instructions

Setup this test case with

    $AMRVAC_DIR/setup.pl -d=23 -phi=0 -z=0 -g=12,12 -p=ff -eos=default -nf=0 -ndust=0 -u=nul -arch=default

# Description

Similarly to stressedIsland, this setup explores the formation of a
current-sheet. Two flux tubes are slightly pushed together, starting the
reconnection process at their interface. Since the currents in the tubes are
parallel, the tubes begin to attract each other which speeds up the formation of
the current-sheet, and the tubes merge into a single one. This setup again uses
force-free magnetodynamics, simulating infinitely magnetised plasma. Tearing
instability sets in at Lundquist numbers of \f$\sim2000\f$ and you can see the
formation of plasmoids and saturation of the reconnection rate.

You can explore how the evolution depends on the perturbation `(eqpar(vcoll_))`
and on the resistivity `(eqpar(kpar_))`. As in stressedIsland, 100 000
test-particles are advanced with the simulation and meaningful quantifications
are done during runtime in `analysis.t`.


