# Test doubleGEM

\test doubleGEM

# Setup instructions

Setup this test case with

    $AMRVAC_DIR/setup.pl -d=23 -phi=0 -z=0 -g=18,18 -p=mhd -eos=default -nf=0 -ndust=0 -u=nul -arch=default

# Description

This test is meant to illustrate the use of resistive to Hall-MHD evolutions. It
simulates the double periodic, double GEM test (this test generalizes the
standard Geospace Environment Modeling Challenge, a standard reconnection setup,
to a setup where energy conservation can be verified easily, see *Resistive MHD
reconnection: resolving long-term, chaotic dynamics*, R. Keppens, O. Porth, K.
Galsgaard, J.T. Frederiksen, A.L. Restante, G. Lapenta, \& C. Parnell, 2013,
Phys. of Plasmas 20, 092109 (17pp). `doi: 10.1063/1.4820946`).

The setup here illustrates the use of Hall MHD, where one additionally needs to
activate Hall-physics in the `definitions.h` file. Also, it shows how to use the
GLM approach to monopole control, which also requires defining in
`definitions.h`. Look for

    #define GLM
    #define HALL

and spot the corresponding code parts in the user file (amrvacusr.t) or the
source code (once compiled, these are in `F90sources/`, where you can inspect
the `amrvacphys.f` file, compared to the `amrvacphys.t` file from the
`\$AMRVAC_DIR/src/mhd` directory).

Note that our (explicit) treatment of the Hall term makes the Hall run
(`doublegemmhdrunA`) much more challenging than the resistive MHD run
(`doublegemmhdrunB`).

See how we use the `iprob` problem switch in the amrvacusr.t and the parfiles,
to use the same compiled code for running different problems.

We set up the problems to use fourth order finite differences here (using MP5
limiting) and strong stability preserving Runge-Kutta schemes. Note the use of
extra ghost cells (`dixB`).

The test is done in 2.5D (invariance in the ignored z-direction), note that the
resistive run keeps \f$B_z\f$ and \f$v_z\f$ zero at all times, while Hall
physics causes the generation of out-of-plane vector components.

A visco-resistive variant of this test was discussed in *MPI-AMRVAC for solar
and astrophysics*, O. Porth, C. Xia, T. Hendrix, S.P. Moschou, \& R. Keppens,
2014, ApJS 214, 4 (26pp) `doi:10.1088/0067-0049/214/1/4`.


