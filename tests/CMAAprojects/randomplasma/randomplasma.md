# Test randomplasma

\test randomplasma
\todo Describe what this test does
\todo The test does not compile

    eigensystem.f:323.4:
    jcobi(ix1,1:nwprim,1:nwprim)=matmul(ritmp(ix1,1:nwwave,1:nwprim),&
    1
    Error: Different shape for array assignment at (1) on dimension 1 (6 and 8)

# Setup instructions

Setup this test case with

    $AMRVAC_DIR/setup.pl -d=12 -phi=0 -z=0 -g=14 -p=mhd -eos=default -nf=0 -ndust=0 -u=nul -arch=default

# Description

This test has no description yet.


