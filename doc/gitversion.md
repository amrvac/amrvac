# How to migrate to the git version

## Step 1

Go through the page [Getting started](getting_started.md)

## Step 2

New switches for setup.pl (former setamrvac): For MHD and HD, the equation of
state is controlled via: -eos=gamma/iso where gamma is the ideal gas with
adiabatic index eqpar(gamma_) iso inherits the functionality from the hdadiab or
mhdiso modules. Currently -eos=default will result in -eos=gamma (but this could
change in the relativistic modules). This is still to be implemented for the
relativistic modules sr(m)hd and sr(m)hdeos -nf=0(default) :: number of tracers
which move with fluid element adds variables tr1_ to trn_.Select regions that
you want to trace and let w(ix^D,tr1_)=one (primitive variable) or
w(ix^D,Dtr1_)=w(ix^D,rho_) (conservative variable) there. In .par file, add
names for tr1_ in primnames and wnames and a suitable boundary type for tr1_ in
parameter typeB. -arch=default/macbook/debug/...:: select a compiling mode
defined in arch/*.defs. you can easily add more machine specific definition
files here. -ndust=0(default) :: number of dust species to use in hd physics
module

Make a folder for your own problem anywhere you want, for example,
~/simulations/CME/. In this folder type

    $AMRVAC_DIR/setup -d=33 -z=0 -phi=0 -g=16,16,16 -p=hd -u=CME -s

to set up the compiling environment.
Your original usr files _$AMRVAC_DIR/src/usr/amrvacusr.t.CME_ and
_amrvacusrpar.t.CME_ will be copied into this folder ~/simulation/CME/ as
amrvacusr.t and amrvacusrpar.t.
From then on you just need to modify these two files to change your code.
Besides, you need to create a parameter file named as amrvac.par (or a
symbolic link file linked to your parameter file with another name) in the
folder.

Then compile your code:

    make

and run on 2 processors:

    mpirun -np 2 amrvac

If you changed the compiling environment, you should make clean before
compile. You may also create a folder datamr/ in ~/simulations/CME/ to output
data in it. By default, the output will be stored in dataxxxx.dat in the
simulation directory and since we recommend to run each simulation in is own
folder (~/simulations/CME1/, ~/simulations/CME2/ ...) creating a datamr is not
longer required.

## Step 3

Modify your amrvacusr.t The way to include modules at the beginning of
amrvacusr.t is changed a little:

    INCLUDE:amrvacusr.gravity.t => INCLUDE:amrvacmodules/gravity.t
    INCLUDE:amrvacnul.specialbound.t => INCLUDE:amrvacnul/specialbound.t

et al. A new module must be included:

    INCLUDE:amrvacnul/usrflags.t

Variable iw can not be used in subroutine specialbound_usr anymore. Just
delete "select case(iw) .. case(rho_) .. end select". And an new subroutine
bc_int must be included. See an example in

    amrvacnul/specialbound.t

Add x argument in these subroutines:

    call conserve(ixG^L,ix^L,w,patchw) => call conserve(ixG^L,ix^L,w,x,patchw)
    call primitive(ixI^L,ixI^L,w) => call primitive(ixI^L,ixI^L,w,x)
    call getpthermal(w,ixG^L,ixO^L,pth) => call getpthermal(w,x,ixG^L,ixO^L,pth)
    call smallvalues(w,ixI^L,ixO^L,"aa") => call smallvalues(w,x,ixI^L,ixO^L,"aa")

## Step 4

Modify your amrvac.par The way adding split and unsplit source terms is
improved. Parameters "sourcesplit" and "sourceunsplit" are removed. The sources,
if any, can be added in a split or unsplit way according to the logical
parameters "ssplitdivb", "ssplitdust", "ssplitresis", and "ssplituser" which
correspond to divb source to maintain divergence free of magnetic field, dust
effect, resistivity, and other sources added by user, respectively. The way of
cleaning the divergence of magnetic field divB is improved and rearranged. The
parameter "divfix" is removed.
_typedivbfix='powel'(default)/'janhunen'/'linde'/'none'/'glm1'(default)/'glm2'/'glm3'_
If choose 'linde', two extra parameters will work: _divbdiff=0.5d0(default)_ ::
coefficient to control diffusion speed _typedivbdiff='all'(default)/'ind'_ ::
add/exclude divB diffusive term in energy equation. If you choose
_'glm1'/'glm2'/'glm3'_, a line,

    #define GLM

must exits in definitions.h, a new variable _psi_ must exist in _amrvac.par_
after magnetic field _b2_ (if 2 directions) or _b3_ (if 3 directions), and its
corresponding boundary type _typeB_ must be set up (as _'cont'_ or
_'periodic'_). If choose other than glm, _#define GLM_ in definitions.h must
be changed, e.g, into _#undefine GLM_ or deleted.

## Step 5

Check `definitions.h`. This is used to define additional switches. Currently
there are the following switches: switch on Hall MHD (_#define HALL_), switch on
the GLM treatment (_#define GLM_), add synchrotron cooling to srmhd and srhdeos
physics modules (_#define EPSINF_) or switch the binary vtu output to big endian
(#define BIGENDIAN). This list will grow. An empty definitions.h should behave
as it used to. If definitions.h is changed, you need to recompile the code to
make it effective.

## Step 6

A new parameter _ditregrid_ is introduced to reconstruct the whole AMR grids once every ditregrid iteration(s) instead of regridding every iteration.

## Step 7

subroutine specialbound_usr for special boundaries in amrvacusr.t should be modified, as the select case(iw) has been removed. A correct example is given as follows:

    subroutine specialbound_usr(qt,ixI^L,ixO^L,iw,iB,w,x)

    ! special boundary types, user defined
    ! user must assign conservative variables in bounderies

    include 'amrvacdef.f'

    integer, intent(in) :: ixI^L, ixO^L, iw, iB
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    !----------------------------------------------------------------------------

    ! just to give an example for 3D MHD :
    select case(iB)
     case(1)
       ! min boundary in the 1st dimension
       w(ixO^S,rho_)=1.d0
       w(ixO^S,v1_)=0.d0
       w(ixO^S,v2_)=0.d0
       w(ixO^S,v3_)=0.d0
       w(ixO^S,p_)=2.d0
       w(ixO^S,b1_)=1.d0
       w(ixO^S,b2_)=0.d0
       w(ixO^S,b3_)=0.d0
       call conserve(ixI^L,ixO^L,w,x,patchfalse)
     case(2)
       ! max boundary in the 1st dimension
       w(ixO^S,rho_)=1.d0
       w(ixO^S,v1_)=0.d0
       w(ixO^S,v2_)=0.d0
       w(ixO^S,v3_)=0.d0
       w(ixO^S,p_)=2.d0
       w(ixO^S,b1_)=1.d0
       w(ixO^S,b2_)=0.d0
       w(ixO^S,b3_)=0.d0
       call conserve(ixI^L,ixO^L,w,x,patchfalse)
     ...
    end select
    end subroutine specialbound_usr

Last modified: Fri Oct 30 10:50:00 CEST 2012
