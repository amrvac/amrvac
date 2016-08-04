# Test promRTI

\test promRTI
\todo Describe what this test does
\todo The test does not compile

    amrvacusr.f:888.26:
    call random_seed(put=randarr)
                          1
    Error: Size of 'put' argument of 'random_seed' intrinsic at (1) too small (2/12)

# Setup instructions

Setup this test case in 2D with

    $AMRVAC_DIR/setup.pl -d=23 -phi=0 -z=0 -g=16,16 -p=mhd

# Description

This test has no description yet.


