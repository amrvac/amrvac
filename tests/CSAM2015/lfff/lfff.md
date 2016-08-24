# Test lfff

\test lfff
\todo Does not compile:

    eigensystem.f:430.68:
    call gradient(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,p_),ixImin1,&
                                                                    1
    Error: Symbol 'p_' at (1) has no IMPLICIT type

# Setup instructions

Setup this test case with

    $AMRVAC_DIR/setup.pl -d=33 -phi=0 -z=0 -g=16,16,16 -p=mhd -eos=iso

# Description

This is a test to demonstrate how to generate a potential or linear force-free
field extrapolation from a given (HMI) magnetogram. One can download such
magnetograms from the SDO/HMI website, you then get these in FITS format. An IDL
script `converthmi.pro` is provided, that assumes you have an IDL license and
access to the library solar software (SSW) for SDO data handling. This routine
just converts the FITS file into a (binary, uniform Cartesian grid) dat file,
containing the vertical magnetic field component only. For convenience, an
already converted magnetogram (in *dat) format is provided as well.

The user file shows how we use the LFFF module (see the
`INCLUDE::amrvacmodules/lfff.t` statement) to then generate a potential or
linear force-free (look for the `alpha` parameter) magnetic field extrapolation.
You should set the code to a 3D isothermal MHD run (use the `-eos=iso` switch),
and the parfile is set up to stop when the initial condition for a possible
dynamical evolution is generated (i.e., the only thing that is computed here is
the LFFF extrapolation, you should use Paraview to visualize the magnetogram and
some selected fieldlines).

Note: to see fieldlines properly in Paraview, we need corner-valued data. This
can be done directly on converting to vtu, by selecting e.g.
`convert\_type='vtuBmpi'`. You can also let the code dump the actual
cell-centered values `convert\_type='vtuBCCmpi'`, but then need to use a
CelldatatoPointdata filter, before you can start drawing fieldlines. Another
quirk of paraview is that it needs the Interpolator type to be set to
`Interpolator with Cell locator`, before the Streamline tracer will be able to
handle field lines across resolution changes (which we have in this 3-level AMR
extrapolated field).

A description of this functionality and how the LFFF extrapolation works
technically is found in (the references of) *MPI-AMRVAC for solar and
astrophysics*, O. Porth, C. Xia, T. Hendrix, S.P. Moschou, \& R. Keppens, 2014,
ApJS 214, 4 (26pp) `doi:10.1088/0067-0049/214/1/4`.


