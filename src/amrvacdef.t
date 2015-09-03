!=============================================================================
! include file amrvacdef.t

use mod_indices
use mod_physicaldata
use mod_connectivity

IMPLICIT NONE

! DEFINITIONS OF GLOBAL PARAMETERS AND VARIABLES
! Parameters:

! Indices for cylindrical coordinates FOR TESTS, negative value when not used:
INTEGER,PARAMETER:: r_=1, phi_=^PHI, z_=^Z

! Indices for cylindrical coordinates FOR INDEXING, always positive
INTEGER,PARAMETER:: pphi_=^PPHI, zz_=^ZZ

include 'amrvacpar.f'

INTEGER,PARAMETER:: ndim=^ND, ndir=^NC

include 'amrvacsettings.f'

INTEGER,PARAMETER:: filelog_=1,fileout_=2,fileslice_=3,filecollapse_=4,fileanalysis_=5,nfile=5 ! outputfiles
INTEGER,PARAMETER:: nslicemax=1000

INTEGER,PARAMETER:: unitstdin=5,unitterm=6,uniterr=6 ! Unit names.

! Units reserved for files:
INTEGER,PARAMETER:: unitpar=9
INTEGER,PARAMETER:: unitconvert=10
INTEGER,PARAMETER:: unitslice=11
INTEGER,PARAMETER:: unitsnapshot=12
INTEGER,PARAMETER:: unitcollapse=13
INTEGER,PARAMETER:: unitanalysis=14

INTEGER,PARAMETER:: biginteger=10000000

! Note: smalldouble must be above machine precision 
DOUBLE PRECISION,PARAMETER:: smalldouble=1.D-12, bigdouble=1.D+99
DOUBLE PRECISION,PARAMETER:: zero=0D0,one=1D0,two=2D0,half=0.5D0,quarter=0.25D0,third=0.33333333333333333333d0
DOUBLE PRECISION,PARAMETER:: dpi=3.141592653589793238462643383279502884197169399375105d0

! Physical scaling parameters:
COMMON, DOUBLE PRECISION:: UNIT_LENGTH, UNIT_DENSITY, UNIT_VELOCITY

include 'amrvacusrpar.f'

! For transform variables and save selected data
COMMON, INTEGER :: nwtf
COMMON, INTEGER :: neqpartf

!Kronecker delta and Levi-Civita tensors
COMMON, INTEGER:: kr(3,3),lvc(3,3,3)

!Equation and method parameters
COMMON, DOUBLE PRECISION:: eqpar(neqpar+nspecialpar)

! Time step control parameters
COMMON, DOUBLE PRECISION :: courantpar, dtpar, dtdiffpar, dtTCpar{#IFDEF MAGNETOFRICTION ,cmf_c,cmf_y,cmf_divb}
COMMON, CHARACTER*131 :: typecourant,typeresid
COMMON, LOGICAL :: time_accurate, addmpibarrier

!Time parameters
INTEGER,PARAMETER:: nsavehi=100       ! maximum No. saves into outputfiles
                                      ! defined by arrays of tsave or itsave
COMMON, DOUBLE PRECISION:: t,tmax,dtmin,residmin,residmax,residual
COMMON, DOUBLE PRECISION:: tfixgrid
COMMON, DOUBLE PRECISION:: tsave(nsavehi,nfile),tsavelast(nfile),dtsave(nfile),slicecoord(nslicemax)
COMMON, LOGICAL:: tmaxexact,treset,itreset,firstprocess,resetgrid,fixprocess,changeglobals,collapse(ndim)
COMMON, INTEGER:: it,itmax,itmin,slowsteps{#IFDEF MAGNETOFRICTION , mfitmax}
COMMON, INTEGER:: itsave(nsavehi,nfile),itsavelast(nfile),ditsave(nfile)
COMMON, INTEGER:: isavet(nfile),isaveit(nfile), nslices, slicedir(nslicemax), collapseLevel
COMMON, INTEGER:: n_saves(1:nfile)
COMMON, INTEGER:: typeparIO
COMMON, INTEGER:: itfixgrid,ditregrid
COMMON, INTEGER:: nwauxio
COMMON, INTEGER:: istep, nstep

!Method switches
COMMON, CHARACTER*131 :: typeadvance
COMMON, CHARACTER*131 :: typelow1(nlevelshi),typelimited,typesourcesplit
COMMON, CHARACTER*131 :: typefull1(nlevelshi), typepred1(nlevelshi)
COMMON, CHARACTER*131 :: typelimiter1(nlevelshi),typegradlimiter1(nlevelshi)
COMMON, CHARACTER*131 :: typelimiter,typegradlimiter,typeprolonglimit
COMMON, CHARACTER*131 :: typeentropy(nw),typetvd,typetvdlf,typeaverage
COMMON, CHARACTER*131 :: typedimsplit,typeaxial,typecoord,typepoly
COMMON, INTEGER:: errorestimate,nxdiffusehllc,typespherical,ncyclemax
COMMON, DOUBLE PRECISION:: entropycoef(nw)
COMMON, DOUBLE PRECISION:: tvdlfeps, mcbeta, parastsnu, TCphi
COMMON, LOGICAL:: sourceparasts,sourceimpl,sourceimplcycle,conduction,TCsaturate,bcphys
COMMON, LOGICAL:: loglimit(nw),logflag(nw),flathllc,flatcd,flatsh,flatppm
COMMON, LOGICAL:: ssplitdust,ssplitdivb,ssplitresis,ssplituser,useprimitive,dimsplit
COMMON, LOGICAL:: restrictprimitive,prolongprimitive, &
                  coarsenprimitive,useprimitiveRel, &
                  amrentropy

COMMON, LOGICAL:: divbwave,compactres,BnormLF
COMMON, DOUBLE PRECISION:: divbdiff,smallT,smallp,smallrho,amr_wavefilter(nlevelshi)
COMMON, CHARACTER*131 :: typedivbdiff,typedivbfix,typediv,typegrad

COMMON, LOGICAL:: fixsmall,strictnr,strictsmall,strictzero,strictgetaux
COMMON, DOUBLE PRECISION::dmaxvel,tolernr,absaccnr,tlow
COMMON, INTEGER:: maxitnr,nflatgetaux

COMMON, LOGICAL:: nocartesian
COMMON, LOGICAL:: writew(nw),writelevel(nlevelshi)
COMMON, DOUBLE PRECISION:: writespshift(ndim,2)
COMMON, INTEGER:: level_io, level_io_min, level_io_max

! cooling related parameters
COMMON, INTEGER:: ncool, cmulti
COMMON, CHARACTER*131 :: coolcurve,coolmethod
COMMON, DOUBLE PRECISION :: cfrac
COMMON, LOGICAL :: Tfix

! dust related paramters
COMMON, LOGICAL  :: dustzero
COMMON, DOUBLE PRECISION :: smallrhod
COMMON, CHARACTER*131 :: dustmethod,dustspecies,dusttemp

! local and global fastest wave speed (computed in setdt):
COMMON, DOUBLE PRECISION :: cmax_mype, cmax_global

!Gravity related parameters
COMMON, DOUBLE PRECISION ::  x1ptms,x2ptms,x3ptms,ptmass

!Boundary region parameters
INTEGER,PARAMETER:: nhiB=2*ndim         ! maximum No. boundary sections
COMMON, LOGICAL:: periodB(ndim), poleB(2,ndim), aperiodB(ndim)
COMMON, CHARACTER*131 :: typeB(nw,nhiB)
COMMON, CHARACTER*131 :: typeghostfill,typegridfill
COMMON, DOUBLE PRECISION::ratebdflux
COMMON, LOGICAL:: internalboundary

!File parameters
COMMON, CHARACTER*131 :: inifile,filenameout,filenameini,filenamelog
COMMON, CHARACTER*131 :: fileheadout
COMMON, CHARACTER*1024 :: wnames,primnames,wnameslog
COMMON, CHARACTER*131 :: typefilelog
COMMON, INTEGER :: snapshotini
COMMON, LOGICAL :: sliceascii

!Convert parameters
COMMON, LOGICAL :: convert,autoconvert,saveprim,uselimiter,endian_swap
COMMON, CHARACTER*131 :: convert_type, dxfiletype, collapse_type
COMMON, DOUBLE PRECISION :: normvar(0:nw),normt
! --------------------------------------------
!Test parameters
COMMON, CHARACTER*131 :: teststr
COMMON, INTEGER :: ixtest1,ixtest2,ixtest3,iwtest,idimtest
COMMON, INTEGER:: saveigrid
! Stores the memory and load imbalance, to be used in printlog:
COMMON, DOUBLE PRECISION :: Xload, Xmemory

LOGICAL:: oktest    !This is a local variable for all subroutines and functions

COMMON, DOUBLE PRECISION:: time_bc

integer,parameter:: nodehi=^ND+1
integer,parameter:: plevel_=1
integer,parameter:: pig^D_=plevel_+^D

integer,parameter:: rnodehi=3*^ND
integer,parameter:: rpxmin0_=0
integer,parameter:: rpxmin^D_=rpxmin0_+^D 
integer,parameter:: rpxmax0_=^ND
integer,parameter:: rpxmax^D_=rpxmax0_+^D 
integer,parameter:: rpdx^D_=2*^ND+^D

! parameters for bc_phys
integer,parameter:: ismin^D=-1+2*^D
integer,parameter:: ismax^D=2*^D

include 'mpif.h'

!-----------------------------------------------------------------------------
! common block variables
!-----------------------------------------------------------------------------

! nodal info per grid: 
! rnode:            corner coordinates,dx(idim)
! node:             dimensions and pointers 

! mxnest:           maximal number of levels

! tol:              error tolerance used for refinement decision
! nbufferx^D:       number of cells as buffer zone

! dixB:             number of ghost cells surrounding a grid for boundary cndt.

double precision :: tol, tolratio
double precision :: rnode, rnode_sub, dx, dt, dtimpl, dt_grid, dxlevel
integer          :: mxnest, dixB, nbufferx^D
integer          :: node, node_sub, levmin, levmax, levmax_sub
logical          :: skipfinestep, time_advance

common /nodalr/     rnode(rnodehi,ngridshi), rnode_sub(rnodehi,ngridshi), tol(nlevelshi), tolratio(nlevelshi), &
                    dx(ndim,nlevelshi), dt, dtimpl, dt_grid(ngridshi), dxlevel(ndim)
common /nodali/     node(nodehi,ngridshi), node_sub(nodehi,ngridshi), &
                    nbufferx^D, mxnest, dixB, levmin, levmax, levmax_sub
common /nodall/     skipfinestep, time_advance

! end include file amrvacdef.t
!============================================================================
