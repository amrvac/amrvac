!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_datatypes
!! NAME
!! defs_datatypes
!!
!! FUNCTION
!! This module contains definitions of all structured datatypes for the
!! ABINIT package.
!! If you want to add one new datatype, please, examine first whether
!! another datatype might meet your need (e.g. adding some records to it).
!! Then, if you are sure your new structured datatype is needed,
!! write it here, and DOCUMENT it properly (not all datastructure here are
!! well documented, it is a shame ...).
!! Do not forget : you will likely be the major winner if you document
!! properly.
!! Proper documentation of a structured datatype means :
!!  (1) Mention it in the list just below
!!  (2) Describe it in the NOTES section
!!  (3) Put it in alphabetical order in the the main section of this module
!!  (4) Document each of its records, except if they are described elsewhere
!!      (this exception is typically the case of the dataset associated with
!!      input variables, for which there is a help file)
!!
!! List of datatypes :
!! * aim_dataset_type : the "dataset" for aim
!! * anaddb_dataset_type : the "dataset" for anaddb
!! * bandstructure_type : different information about the band structure
!! * bcp_type : a "bonding critical point" for aim
!! * dataset_type : the "dataset" for the main abinit code
!! * datafil_type : the data (units,filenames) related to files
!! * dens_sym_operator_type : information for symmetrizing the density
!! * efield_type : First-principles calculations in a finite electric field
!! * electronic_structure : for GW part of ABINIT, energies, occupations,
!!   wavefunctions in real and reciprocal space (big set of data !)
!! * epsilonm1_parameters : for GW part of ABINIT, parameters for epsilon-1
!! * epsilonm1_results : for GW part of ABINIT, results of screening
!! * gs_hamiltonian_type : datastructure describing an Hamiltonian
!! * hdr_type   : the header of wf, den and pot files
!! * MPI_type : the data related to MPI parallelization
!! * pawang_type : for PAW, ANGular mesh discretization and related data
!! * pawfgr_type : for PAW, Fine rectangular GRid parameters and related data
!! * pawfgrtab_type : for PAW, various arrays giving data related to fine grid for a given atom
!! * pawrad_type : for PAW, RADial mesh discretization and related data
!! * pawtab_type : for PAW, TABulated data initialized at start
!! * paw_an_type : for PAW, various arrays given on ANgular mesh or ANgular moments
!! * paw_ij_type : for PAW, various arrays given on (i,j) (partial waves)
!!   channels
!! * pseudopotential_type : for norm-conserving pseudopotential, all the
!!   information
!! * pspheader_paw_type : for PAW, the header of the atomic file
!! * pspheader_type : for norm-conserving pseudopotentials, the header of
!!   the file
!! * results_gs_type : contains the results of a GS calculation
!! * results_out_type : contains a subset of the results, for internal
!!   tests and timing analysis
!! * sigma_parameters : for GW part of ABINIT, parameters for sigma
!! * sigma_results : for GW part of ABINIT, results of sigma
!! * wffile_type : a handler for dealing with the IO of a wavefunction file
!!
!! COPYRIGHT
!! Copyright (C) 2001-2006 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! PAW structured datatypes to be described ...
!! * pawang_type : ANGular mesh discretization and related data
!! * pawfgr_type : Fine rectangular GRid parameters and related data
!! * pawfgrtab_type : various arrays giving data related to fine grid for a given atom
!! * pawrad_type :  RADial mesh discretization and related data
!! * pawtab_type : TABulated data initialized at start
!! * paw_an_type : various arrays given on ANgular mesh or
!! * paw_ij_type : various arrays given on (i,j) (partial waves) channels
!! * pspheader_paw_type: the header of the atomic file
!!
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_datatypes

 use defs_basis

 implicit none

!Structures

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/aim_dataset_type
!! NAME
!! aim_dataset_type
!!
!! FUNCTION
!! The aim_dataset_type structured datatype
!! gathers all the input variables for the aim code
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type aim_dataset_type

! Since all these input variables are described in the aim_help.html
! file, they are not described in length here ...

! Integer
  integer :: crit,denout,dltyp,gpsurf,irho,ivol,lapout,nsa,nsb,nsc
  integer :: ngrid(3)
  integer :: batom  !! Warning : corresponds to the input variable atom
  integer :: foll   !! Warning : corresponds to the input variable follow
  integer :: isurf  !! Warning : corresponds to the input variable surf
  integer :: irsur  !! Warning : corresponds to the input variable rsurf
  integer :: nph    !! Warning : corresponds to the input variable nphi
  integer :: npt    !! Warning : corresponds to the input variable inpt
  integer :: nth    !! Warning : corresponds to the input variable ntheta
  integer :: plden  !! Warning : not documented in help file ?!

! Real
  real(dp) :: atrad,coff1,coff2,dpclim,folstp,lgrad,lgrad2,lstep,lstep2,&
&  maxatd,maxcpd,phimax,phimin
  real(dp) :: foldep(3),scal(3),vpts(3,4)
  real(dp) :: dr0    !! Warning : correspond to the input variable radstp
  real(dp) :: phi0   !! Warning : correspond to the input variable rsurdir(2)
  real(dp) :: rmin   !! Warning : correspond to the input variable ratmin
  real(dp) :: th0    !! Warning : correspond to the input variable rsurdir(1)
  real(dp) :: themax !! Warning : correspond to the input variable thetamax
  real(dp) :: themin !! Warning : correspond to the input variable thetamin

 end type aim_dataset_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/anaddb_dataset_type
!! NAME
!! anaddb_dataset_type
!!
!! FUNCTION
!! The anaddb_dataset_type structured datatype
!! gather all the input variables for the anaddb code.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type anaddb_dataset_type

! Since all these input variables are described in the anaddb_help.html
! file, they are not described in length here ...
! Integer
integer :: alphon,asr,brav,chneut,dieflag,dipdip,doscalprod,eivec,elaflag,elphflag,enunit
  integer :: gkk2exist,gkk2write,gkk_rptexist,gkk_rptwrite,gkqexist,gkqwrite
  integer :: ifcana,ifcflag,ifcout,instrflag,natfix,natifc,natom
  integer :: nchan,nfreq,ngrids,nlflag,nph1l,nph2l,nqpath
  integer :: nqshft,nsphere,nstrfix,ntemper,nwchan
  integer :: phfrqexist,phfrqwrite,piezoflag,polflag,prtmbm,prtfsurf,prtnest,ramansr
  integer :: relaxat,relaxstr,rfmeth,selectz,symdynmat,telphint,tkeepbands,thmflag
  integer :: ngqpt(9)             ! ngqpt(9) instead of ngqpt(3) is needed in wght9.f
  integer :: istrfix(6),ng2qpt(3),kptrlatt(3,3)

! Real(dp)
  real(dp) :: a2fsmear,dostol,elphsmear,elph_fermie,frmax,frmin,temperinc,tempermin,thmtol,mustar,rifcsph
  real(dp) :: q1shft(3,4),q2shft(3),targetpol(3)

! Integer pointers
  integer, pointer :: atifc(:)    ! atifc(natom) WARNING : there is a transformation
                                  ! of this input variable, in chkin9
                                  ! This should be changed ...
  integer, pointer :: iatfix(:)   ! iatfix(natom)

! Real pointers
  real(dp), pointer :: qnrml1(:)  ! qnrml1(nph1l)
  real(dp), pointer :: qnrml2(:)  ! qnrml2(nph2l)
  real(dp), pointer :: qpath(:,:) ! qpath(3,nqpath)
  real(dp), pointer :: qph1l(:,:) ! qph1l(3,nph1l)
  real(dp), pointer :: qph2l(:,:) ! qph2l(3,nph2l)

 end type anaddb_dataset_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/bandstructure_type
!! NAME
!! bandstructure_type
!!
!! FUNCTION
!! It contains different information about the band structure
!! (eigenenergies, residuals, derivative of occupation number
!!  vs energy in case of metallic occupations)
!! and Brillouin zone according to the context : k points, occupation numbers,
!! storage mode of wavefunctions, weights ...
!! For example, the initial Brillouin zone, set up in the dataset, will be treated
!! in the response function part of the code, to give a reduced
!! Brillouin zone different from the original one, due to the
!! breaking of the symmetries related to the existence of a wavevector,
!! or the lack of time-reversal invariance
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type bandstructure_type

  integer :: bantot                  ! total number of bands (sum(nband(:))
  integer :: nkpt                    ! number of k points
  integer :: nsppol                  ! number of spin-polarizations
  integer, pointer :: istwfk(:)      ! istwfk(nkpt) storage mode at each k point
  integer, pointer :: nband(:)       ! nband(nkpt*nsppol) number of bands
                                     !    at each k point and spin-polarisation
  integer, pointer :: npwarr(:)      ! npwarr(nkpt) number of plane waves at each k point
  real(dp), pointer :: kptns(:,:)    ! kptns(3,nkpt)  k-point vectors
  real(dp), pointer :: eig(:)        ! eig(bantot)  eigenvalues of each band
  real(dp), pointer :: occ(:)        ! occ(bantot)  occupation of each band
  real(dp), pointer :: doccde(:)     ! doccde(bantot)  derivative of the
                                     !    occupation of each band wrt energy (needed for RF)
  real(dp), pointer :: wtk(:)        ! wtk(nkpt)  weight of each k point

 end type bandstructure_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/bcp_type
!! NAME
!! bcp_type
!!
!! FUNCTION
!! a "bonding critical point" for aim
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type bcp_type

! Integer
  integer :: iat     !! number of the bonding atom inside a primitive cell
  integer :: ipos    !! number of the primitive cell of the bonding atom

! Real
  real(dp) :: chg     !! charge at the critical point
  real(dp) :: diff(3) !! three distances : AT-CP,BAT-CP,AT-BAT
  real(dp) :: ev(3)   !! eigenvalues of the Hessian
  real(dp) :: pom(3)  !! position of the bonding atom
  real(dp) :: rr(3)   !! position of the bcp
  real(dp) :: vec(3,3)!! eigenvectors of the Hessian
  real(dp) :: vv(3)   !! position of the bcp relative to the central atom

 end type bcp_type

!----------------------------------------------------------------------

!!****t* defs_datatypes/datafiles_type
!! NAME
!! datafiles_type
!!
!! FUNCTION
!! The datafiles_type structures datatype gather all the variables related
!! to files, such as filename, and file units.
!! For one dataset, it is initialized in driver.f, and will not change
!! at all during the treatment of the dataset.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type datafiles_type

  integer :: ireadden
   !   ireadden non-zero  if the den file must be read
  integer :: ireadwf
   ! if(optdriver/=1), that is, no response-function computation,
   !   ireadwf non-zero  if the wffk file must be read
   !   (if irdwfk non-zero or getwfk non-zero)
   ! if(optdriver==1), that is, response-function computation,
   !   ireadwf non-zero  if the wff1 file must be read
   !   (if ird1wf non-zero or get1wf non-zero)
  integer :: unddb   ! unit number for Derivative DataBase
  integer :: unddk   ! unit number for ddk 1WF file
  integer :: unkg    ! unit number for k+G data
  integer :: unkgq   ! unit number for k+G+q data
  integer :: unkg1   ! unit number for first-order k+G+q data
  integer :: unwff1  ! unit number for wavefunctions, number one
  integer :: unwff2  ! unit number for wavefunctions, number two
  integer :: unwffgs ! unit number for ground-state wavefunctions
  integer :: unwffkq ! unit number for k+q ground-state wavefunctions
  integer :: unwft1  ! unit number for wavefunctions, temporary one
  integer :: unwft2  ! unit number for wavefunctions, temporary two
  integer :: unwftgs ! unit number for ground-state wavefunctions, temporary
  integer :: unwftkq ! unit number for k+q ground-state wavefunctions, temporary
  integer :: unylm   ! unit number for Ylm(k) data
  integer :: unylm1  ! unit number for first-order Ylm(k+q) data
  integer :: unpaw   ! unit number for temporary PAW data (for ex. rhoij_nk) (Paw only)
  integer :: ungsc1  ! unit number for <g|S|c> data (Paw only), temporary one
  integer :: ungsc2  ! unit number for <g|S|c> data (Paw only), temporary two
  integer :: unpos   ! unit number for restart molecular dynamics

  character(len=fnlen) :: filnam_ds(5)
   ! if no dataset mode, the five names from the standard input :
   !   ab_in, ab_out, abi, abo, tmp
   ! if dataset mode, the same 5 filenames, appended with //'_DS'//trim(jdtset)

  character(len=fnlen) :: fildensin
   ! if no dataset mode             : abi//'DEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'DEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'DEN'

  character(len=fnlen) :: filvhain
   ! if no dataset mode             : abi//'VHA'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'VHA'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'VHA'

  character(len=fnlen) :: filkss
   ! if no dataset mode             : abi//'KSS'
   ! if dataset mode, and getkss==0 : abi//'_DS'//trim(jdtset)//'KSS'
   ! if dataset mode, and getkss/=0 : abo//'_DS'//trim(jgetkss)//'KSS'

  character(len=fnlen) :: filqps
   ! if no dataset mode             : abi//'SCR'
   ! if dataset mode, and getqps==0 : abi//'_DS'//trim(jdtset)//'SCR'
   ! if dataset mode, and getqps/=0 : abo//'_DS'//trim(jgetqps)//'SCR'

  character(len=fnlen) :: filscr
   ! if no dataset mode             : abi//'SCR'
   ! if dataset mode, and getscr==0 : abi//'_DS'//trim(jdtset)//'SCR'
   ! if dataset mode, and getscr/=0 : abo//'_DS'//trim(jgetscr)//'SCR'

! character(len=fnlen) :: filpsp(ntypat)
   ! the filenames of the pseudopotential files, from the standard input.

  character(len=fnlen) :: filstat
   ! tmp//'_STATUS'

  character(len=fnlen) :: fnamewffk
   ! the name of the ground-state wavefunction file to be read (see driver.f)

  character(len=fnlen) :: fnamewffq
   ! the name of the k+q ground-state wavefunction file to be read (see driver.f)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffddk
   ! the generic name of the ddk response wavefunction file(s) to be read (see driver.f)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewff1
   ! the generic name of the first-order wavefunction file(s) to be read (see driver.f)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fildens1in   ! to be described by MVeithen

 end type datafiles_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/dataset_type
!! NAME
!! dataset_type
!!
!! FUNCTION
!! The dataset_type structured datatype gather all the input variables,
!! except those that are labelled NOT INTERNAL.
!! For one dataset, it is initialized in driver.f, and will not change
!! at all during the treatment of the dataset.
!! The "evolving" input variables are also stored, with their
!! name appended with _orig, to make clear that this is the original
!! value, decided by the user, and not a possibly modified, intermediate value.
!! The following input variables are NOT INTERNAL, that is, they
!! are input variables used to determine other input variables,
!! after suitable processing, and do not appear anymore afterwards
!! (so, they do not appear as components of a dataset_type variable) :
!! cpuh,cpum(but cpus is present),fband,kptbounds,ndivk,nobj,
!! objaat,objbat,objaax,objbax,objan,objbn,objarf,objbrf,objaro,objbro
!! objatr,objbtr,vaclst,vacuum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type dataset_type
! Since all these input variables are described in the abinis_help.html
! file, they are not described in length here ...
! Integer
  integer :: accesswff,berryopt,brvltt,ceksph,chkexit,chkprim,&
&  delayperm,enunit,exchn2n3,fft_opt_lob,freqremax,freqspmax,frzfermi,getacfd,&
&  getcell,getddk,getden,getkss,getocc,getqps,getscr,getvel,getwfk,&
&  getwfq,getxcart,getxred,get1den,get1wf,gpara,gwcalctyp,iboxcut,&
&  icoultrtmt,idyson,ikhxc,inclvkb,intexact,intxc,ionmov,&
&  iprcch,iprcel,iprcfc,irdddk,irdkss,irdqps,irdscr,irdwfk,irdwfq,ird1wf,&
&  iscf,isecur,istatr,istatshft,ixc,ixcpositron,&
!  jdtset contains the actual number of the dataset
&  jdtset,kpara,kptopt,kssform,ldgapp,localrdwf,lofwrite,mband,mffmem,mgfft,mgfftdg,&
&  mkmem,mkqmem,mk1mem,mpw,mqgrid,mqgriddg,natom,natrd,natsph,nbandsus,nbdblock,nbdbuf,&
&  nberry,nbandkss,ncenter,nconeq,nctime,ndtset,ndyson,&
&  nfft,nfftdg,nfreqim,nfreqre,nfreqsp,nfreqsus,ngroup_rf,nkptgw,nkpt,nline,&
&  nnsclo,nomegasrd,norb,npack,npara,npband,npfft,npsp,npspalch,npulayit,&
&  npweps,npwkss,npwsigx,npwwfn,nqpt,nqptdm,nscforder,&
&  nsheps,nshiftk,nshsigx,nshwfn,nspden,nspinor,nsppol,nstep,nsym,ntime,&
&  ntypalch,ntypat,ntyppure,occopt,optcell,optdriver,&
&  optforces,optfreqsus,optnlxccc,optstress,ortalg,&
&  outputXML,paral_rf,parareel,&
&  pawlcutd,pawlmix,pawmixdg,pawnphi,pawntheta,pawnzlm,pawoptmix,pawprtvol,pawxcdev,&
&  positron,ppmodel,prepanl,prtacfd,prtbbb,prtcml,&
&  prtden,prtdos,prteig,prtfsurf,prtgeo,prtgkk,prtkpt,prtnabla,prtpot,prtstm,prtvha,prtvhxc,prtvol,prtvxc,&
&  prtwant,prtwf,prt1dm,ptgroupma,restartxf,rfasr,rfelfd,rfmeth,rfphon,rfstrs,rfthrd,&
&  rfuser,rf1elfd,rf1phon,rf2elfd,rf2phon,rf3elfd,rf3phon,&
&  signperm,spgaxor,spgorig,spgroup,splitsigc,suskxcrs,symmorphi,td_mexcit,tfkinfunc,timopt,usepaw,&
&  usepawu,useria,userib,useric,userid,userie,useylm,vacnum,wfoptalg
! Integer arrays
  integer :: bdberry(4),dsifkpt(3),kptrlatt(3,3),ngfft(18),ngfftdg(18),nloalg(5),&
&  qprtrb(3),rfatpol(2),rfdir(3),rf1atpol(2),rf1dir(3),&
&  rf2atpol(2),rf2dir(3),rf3atpol(2),rf3dir(3),supercell(3)
! Integer pointers
  integer, pointer ::  algalch(:)    ! algalch(ntypalch)
  integer, pointer ::  bdgw(:,:)     ! bdgw(2,nkptgw)
  integer, pointer ::  iatfix(:,:)   ! iatfix(3,natom)
  integer, pointer ::  iatsph(:)     ! iatsph(natsph)
  integer, pointer ::  istwfk(:)     ! istwfk(nkpt)
  integer, pointer ::  kberry(:,:)   ! kberry(3,nberry)
  integer, pointer ::  lpawu(:)      ! lpawu(ntypat)
  integer, pointer ::  ltypeorb(:)   ! ltypeorb(norb)
  integer, pointer ::  nband(:)      ! nband(nkpt*nsppol)
  integer, pointer ::  numorb(:)     ! numorb(ncenter)
  integer, pointer ::  so_typat(:)   ! so_typat(ntypat)
  integer, pointer ::  symafm(:)     ! symafm(nsym)
  integer, pointer ::  symrel(:,:,:) ! symrel(3,3,nsym)
  integer, pointer ::  typat(:)      ! typat(natom)
! Real
  real(dp) :: alpha,boxcutmin,bxctmindg,charge,cpus,dedlnn,diecut,diegap,dielam,&
&  dielng,diemac,diemix,dilatmx,dosdeltae,dtion,&
&  ecut,ecuteps,ecutsigx,ecutsm,ecutwfn,effmass,&
&  eshift,fband,fixmom,freqsusin,freqsuslo,friction,kptnrm,kptrlen,mdftemp,&
&  mditemp,mdwall,nelect,noseinert,omegasrdmax,pawecutdg,pawovlp,pawsphmix,ppmfrq,qptnrm,ratsph,&
&  sciss,soenergy,stmbias,strfact,strprecon,td_maxene,tfnewton,toldfe,toldff,&
&  tolmxf,tolvrs,tolwfr,tphysel,tsmear,userra,userrb,userrc,userrd,&
&  userre,vacwidth,vis,zcut
! Real arrays
  real(dp) :: acell_orig(3),angdeg_orig(3),boxcenter(3),&
&  efield(3),genafm(3),qpt(3),qptn(3),rprim_orig(3,3),&
&  rprimd_orig(3,3),strtarget(6),vprtrb(2)
! Real pointers
  real(dp), pointer :: amu(:)         ! amu(ntypat)
  real(dp), pointer :: densty(:,:)    ! densty(ntypat,4)
  real(dp), pointer :: jpaw(:)        ! jpaw(ntypat)
  real(dp), pointer :: kpt(:,:)       ! kpt(3,nkpt)
  real(dp), pointer :: kptgw(:,:)     ! kptgw(3,nkptgw)
  real(dp), pointer :: kptns(:,:)     ! kptns(3,nkpt)
  real(dp), pointer :: mixalch(:,:)   ! mixalch(npspalch,ntypalch)
  real(dp), pointer :: occ_orig(:)    ! occ_orig(mband*nkpt*nsppol)
  real(dp), pointer :: qptdm(:,:)     ! qptdm(3,nqptdm)
  real(dp), pointer :: rcoord(:,:)    ! rcoord(3,ncenter)
  real(dp), pointer :: rtheta(:,:)    ! rtheta(3,norb)
  real(dp), pointer :: shiftk(:,:)    ! shiftk(3,nshiftk)
  real(dp), pointer :: spinat(:,:)    ! spinat(3,natom)
  real(dp), pointer :: tnons(:,:)     ! tnons(3,nsym)
  real(dp), pointer :: upaw(:)        ! upaw(ntypat)
  real(dp), pointer :: vel_orig(:,:)  ! vel_orig(3,natom)
  real(dp), pointer :: wtatcon(:,:,:) ! wtatcon(3,natom,nconeq)
  real(dp), pointer :: wtk(:)         ! wtk(nkpt)
  real(dp), pointer :: xred_orig(:,:) ! xred_orig(3,natom)
  real(dp), pointer :: ziontypat(:)   ! ziontypat(ntypat)
  real(dp), pointer :: znucl(:)       ! znucl(npsp)
 end type dataset_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/dens_sym_operator_type
!! NAME
!! dens_sym_operator_type
!!
!! FUNCTION
!! Information for symmetrizing the density
!! This datastructure contains the information needed for symmetrizing
!! the density : number of symmetry operations, location of related
!! points in reciprocal space, phases, etc ...
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type dens_sym_operator_type

! Integer scalar

  integer :: flagdensymop
   ! if 1, the density symmetrization operator is to be applied,
   ! if 0, do not apply it

  integer :: nfft
   ! number of FFT grid points

  integer :: nspdensymop
   ! number of spin-density components of irrzon and phnons

  integer :: nsym
   ! number of symmetries

! Integer arrays

! integer, pointer :: irrzon(:,:,:)
   ! irrzon(nfft*flagdensymop,2,nspdensymop)
   ! contains the locations of related
   ! grid points and the repetition number for each symmetry class.

! integer, pointer :: symafm(:)
   ! symafm(nsym)
   ! anti-ferromagnetic character of the symmetry operations (+ if the
   ! magnetisation is not conserved, -1 if the magnetisation is reversed)

! Real (double precision) arrays

! real(dp), pointer :: phnons(:,:,:)
   ! phnons(2,nfft*flagdensymop,nspdensymop)
   ! phases associated with nonsymmorphic translations

 end type dens_sym_operator_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/efield_type
!! NAME
!! efield_type
!!
!! FUNCTION
!! First-principles calculations in a finite electric field
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type efield_type

! Integer variables
  integer :: fmkmem              ! number of k-points in the FBZ per cpu
  integer :: fmkmem_max          ! max of fmkmem
  integer :: fnkpt               ! number of k-points in the FBZ
  integer :: maxnstr             ! max number of strings along idir=1,2,3
  integer :: maxnkstr            ! max number of k-points per string
  integer :: mkmem_max           ! max of mkmem
  integer :: nband_occ           ! number of occupied bands
                                 ! this number must be the same for every k

! Integer arrays
  integer :: nstr(3)             ! nstr(idir) = number of strings along idir
  integer :: nkstr(3)            ! nkstr(idir) = number of k-points per string

! Real(dp) scalars
  real(dp) :: sdeg               ! spin degeneracy: sdeg = 2 if nsppol = 1
                                 !                         1 if nsppol = 2

! Real(dp) arrays
  real(dp) :: dkvecs(3,3)        ! dkvec(:,idir) = vector between a k-poinit
                                 ! and its nearest neighbour along idir
  real(dp) :: efield_dot(3)      ! reciprocal lattice coordinates of the
                                 ! electric field

! Integer pointers
  integer, pointer :: cgindex(:,:)    ! cgindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cg array
  integer, pointer :: cgqindex(:,:,:) ! cgqindex(3,6,nkpt*nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cgq and pwnsfacq
                                      ! arrays
                                      ! (see vtorho.f and initberry.f)
  integer, pointer :: ikpt_dk(:,:,:)  ! ikpt_dk(nkpt,2,3)
                                      ! ikpt_dp(ikpt,ii,idir) = index of the
                                      ! k-point at k+dk (ii=1) and k-dk (ii=2)
  integer, pointer :: idxkstr(:,:,:)  ! idxkstr(maxnkstr,maxnstr,3)
                                      ! idxkstr(ikstr,istr,idir) index (ikpt) of
                                      ! k-point ikstr on string istr along idir
  integer, pointer :: indkk_f2ibz(:,:)   ! indkk_f2ibz(1:dtefield%fnkpt,1:6)
                                         ! information needed to fold a
                                         ! k-point in the FBZ into the IBZ;
                                         ! the second index (1:6)
                                         ! is as described in listkk
  integer, pointer :: i2fbz(:)           ! i2fbz(1:nkpt) gives index of IBZ
                                         ! k-points in the FBZ k-point list
  integer, pointer :: nneigh(:)          ! nneigh(nkpt)
                                         ! for each k-point, nneigh stores
                                         ! the number of its nearest neighbours
                                         ! that are not related by symmetry
  integer, pointer :: kgindex(:)      ! kgind(nkpt)
                                      ! kgind(ikpt) = ikg
  integer, pointer :: fkgindex(:)     ! same as kgindex, but defined
                                      ! for the FBZ and intended to use
                                      ! with pwindf
  integer, pointer :: sflag(:,:,:,:)  ! sflag(nband_occ,nkpt*nsppol,2,3)
                                      ! sflag = 0 : compute the whole row of
                                      !             smat
                                      ! sflag = 1 : the row is up to date

! Real(dp) pointers
  real(dp), pointer :: fkptns(:,:)       ! fkptns(3,1:dtefield%fnkpt)
                                         ! k-points in FBZ

  real(dp), pointer :: smat(:,:,:,:,:,:)
! smat(2,nband_occ,nband_occ,nkpt*nsppol,2,3)
! Overlap matrix for every k-point. In an electric field calculation,
! smat is updated at every iteration.


 end type efield_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/electronic_structure
!! NAME
!! electronic_structure
!!
!! FUNCTION
!! For the GW part of ABINIT, the electronic_structure structured datatype
!! gather all the energies, occupations, wavefunctions
!!    in real and reciprocal space (big set of data !)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type electronic_structure
  integer nb,nk,ng,nr
  real(dp), pointer :: en(:,:)       ! en(n,k) Bloch order
  real(dp), pointer :: oc(:,:)       ! same
  real(dp), pointer :: wfg(:,:,:)    ! wfg(g,n,k) order
  real(dp), pointer :: wfr(:,:,:)    ! wfr(r,n,k) order
  real(dp) :: fermie
  real(dp) :: etotal, residm, ecut_eff
  integer :: xc
 end type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/epsilonm1_parameters
!! NAME
!! epsilonm1_parameters
!!
!! FUNCTION
!! For the GW part of ABINIT, the  epsilonm1_parameters structured datatype
!! gather different parameters that characterize the inverse dielectric
!! matrices.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type epsilonm1_parameters
  integer :: gwcalctyp
  integer :: npwwfn,npwe
  integer :: nb
  integer :: nk,nkbz,nop
  integer :: nq,nomega
  integer :: nomegaer,nomegaei
  real(dp) :: soenergy
  real(dp) :: omegaermax
 end type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/ epsilonm1_results
!! NAME
!! epsilonm1_results
!!
!! FUNCTION
!! For the GW part of ABINIT, the  epsilonm1_results structured datatype
!! gather the results of screening : the inverse dielectric
!! matrix, and the omega matrices .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type epsilonm1_results
  integer :: nq,nomega,npwe
  complex, pointer :: omega(:)
  complex, pointer :: epsm1(:,:,:,:)
  complex, pointer :: bigomegatwsq(:,:,:),omegatw(:,:,:)
 end type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/gs_hamiltonian_type
!! NAME
!! gs_hamiltonian_type
!!
!! FUNCTION
!! This datastructure contains the information about one Hamiltonian,
!! needed in the "getghc" routine, that apply the Hamiltonian
!! on a wavefunction.
!! About the non-local part of the Hamiltonian
!! The operator Onl has the following general form:
!! $Onl=sum_{R,lmn,l''m''n''} {|P_{Rlmn}> Enl^{R}_{lmn,l''m''n''} <P_{Rl''m''n''}|}$
!! Operator Onl is -- in the typical case -- the nonlocal potential.
!! - In a classical plane-wave calculation, $Enl^{R}_{lmn,l''m''n''}$ is the
!!   Kleinmann-Bylander energy $Ekb^{R}_{ln}$.
!! - In a PAW calculation, $Enl^{R}_{lmn,l''m''n''}$ can either be the nonlocal
!!   contribution to the Hamiltonian $D_{ij}$ or the overlap matrix $S_{ij}$.
!! - The |P_{Rlmn}> are the projector functions.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type gs_hamiltonian_type

! Integer scalar

  integer :: dimekb1
   ! First dimension of Ekb (see ekb in this file)
   ! Same as psps%dimekb
   ! ->Norm conserving : Max. number of Kleinman-Bylander energies
   !                     for each atom type
   !                     dimekb1=lnmax
   ! ->PAW : Max. number of Dij coefficients connecting projectors
   !                     for each atom
   !                     dimekb1=lmnmax*(lmnmax+1)/2

  integer :: dimekb2
   ! Second dimension of Ekb (see ekb in this file)
   ! ->Norm conserving psps: dimekb2=ntypat
   ! ->PAW                 : dimekb2=natom

  integer :: istwf_k
   ! option parameter that describes the storage of wfs

  integer :: lmnmax
   ! Maximum number of different l,m,n components over all types of psps.
   ! same as dtset%lmnmax

  integer :: matblk
   ! dimension of the array ph3d

  integer :: mgfft
   ! maximum size for 1D FFTs (same as dtset%mgfft)

  integer :: mproj  ! TO BE SUPPRESSED LATER
   ! Maximum number of non-local projectors over all angular momenta
   !  and type of psps
   ! 0 only if all psps are local
   ! same as psps%mproj

  integer :: mpsang
   ! Highest angular momentum of non-local projectors over all type of psps.
   ! shifted by 1 : for all local psps, mpsang=0; for largest s, mpsang=1,
   ! for largest p, mpsang=2; for largest d, mpsang=3; for largest f, mpsang=4
   ! This gives also the number of non-local "channels"
   ! same as psps%mpsang

  integer :: mpssoang
   ! Maximum number of channels, including those for treating the spin-orbit coupling
   ! when mpspso=1, mpssoang=mpsang
   ! when mpspso=2, mpssoang=2*mpsang-1
   ! same as psps%mpssoang

  integer :: natom
   ! The number of atoms for this dataset
   ! same as dtset%natom

  integer :: nfft
   ! number of FFT grid points
   ! same as dtset%nfft

  integer :: npw
   ! number of plane waves

  integer :: ntypat
   ! Number of types of pseudopotentials
   ! same as dtset%ntypat

  integer :: nvloc
   ! Number of components of vloc
   ! usually, nvloc=1, except in the non-collinear magnetism case, where nvloc=4

  integer :: n4,n5,n6
   ! same as ngfft(4:6)

  integer :: usepaw
   ! if usepaw=0 , use norm-conserving psps part of the code
   ! is usepaw=1 , use paw part of the code

  integer :: useylm
   ! governs the way the nonlocal operator is to be applied:
   !   1=using Ylm, 0=using Legendre polynomials

! Integer arrays

  integer, pointer :: atindx(:)
   ! atindx(natom)
   ! index table for atoms (see scfcv.f)

  integer, pointer :: atindx1(:)
   ! atindx1(natom)
   ! index table for atoms, inverse of atindx (see scfcv.f)

  integer, pointer :: gbound(:,:)
   ! gbound(2*mgfft+8,2)
   ! G sphere boundary

  integer, pointer :: indlmn(:,:,:)
   ! indlmn(6,lmnmax,ntypat)
   ! For each type of psp,
   ! array giving l,m,n,ln,lm,spin for i=ln  (if useylm=0)
   !                                or i=lmn (if useylm=1)

! integer, pointer :: indpw_k(:,:)
   ! indpw_k(4,npw)
   ! array which gives fft box index for given basis sphere
   ! This component was taken away : CPU time problem !

! integer, pointer :: kg_k(:,:)
   ! kg_k(3,npw)
   ! G vector coordinates with respect to reciprocal lattice translations
   ! This component was taken away : CPU time problem !

  integer, pointer :: nattyp(:)
   ! nattyp(ntypat)
   ! # of atoms of each type

  integer :: ngfft(18)
   ! ngfft(1:3)=integer fft box dimensions
   ! ngfft(4:6)=integer fft box dimensions, might be augmented for CPU speed
   ! ngfft(7)=fftalg
   ! ngfft(8)=fftalg
   ! same as dtset%ngfft

  integer :: nloalg(5)
   ! governs the choice of the algorithm for non-local operator
   ! same as dtset%nloalg

  integer, pointer :: pspso(:)
   ! pspso(ntypat)
   ! For each type of psp, 1 if no spin-orbit component is taken
   ! into account, 2 if a spin-orbit component is used

! Real (double precision) scalar

  real(dp) :: ucvol
   ! unit cell volume (Bohr**3)

! Real (double precision) arrays

  real(dp), pointer :: ekb(:,:)
   ! ekb(dimekb1,dimekb2)
   !  ->Norm conserving : (Real) Kleinman-Bylander energies (hartree)
   !          for number of basis functions (l,n) (lnmax)
   !          and number of atom types (ntypat)
   !          dimekb1=lnmax ; dimekb2=ntypat
   !  ->PAW : (Real, symmetric) Frozen part of Dij coefficients
   !                            to connect projectors
   !          for number of basis functions (l,m,n) (lmnmax)
   !          and number of atom (natom)
   !          dimekb1=lmnmax*(lmnmax+1)/2 ; dimekb2=natom
   ! NOTE (MT) : ekb (norm-conserving) is now diagonal (one dimension
   !             lnmax); it would be easy to give it a second
   !             (symmetric) dimension by putting
   !             dimekb1=lnmax*(lnmax+1)/2
   !             in the place of dimekb1=lnmax.

  real(dp), pointer :: sij(:,:)
   ! sij(dimekb1,ntypat*usepaw) = overlap matrix for paw calculation

! real(dp), pointer :: ffnl(:,:,:,:)
   ! ffnl(npw,2,lmnmax,ntypat)
   ! nonlocal form factors
   ! This component was taken away : CPU time problem !

  real(dp) :: gmet(3,3)
   ! reciprocal space metric tensor in Bohr**-2

  real(dp) :: gprimd(3,3)
   ! dimensional reciprocal space primitive translations (Bohr^-1)

! real(dp), pointer :: kinpw(:)
   ! kinpw(npw)
   ! (modified) kinetic energy for each plane wave (Hartree)
   ! This component was taken away : CPU time problem !

  real(dp) :: kpoint(3)
   ! dimensionless k point coordinates wrt reciprocal lattice vectors

  real(dp), pointer :: phkxred(:,:)
   ! phkxred(2,natom)
   ! phase factors exp(2 pi k.xred)

  real(dp), pointer :: ph1d(:,:)
   ! ph1d(2,3*(2*mgfft+1)*natom)
   ! 1-dim phase arrays for structure factor (see getph.f).

! real(dp), pointer :: ph3d(:,:,:)
   ! ph3d(2,npw,matblk)
   ! 3-dim structure factors, for each atom and plane wave
   ! This component was taken away : CPU time problem !

! real(dp), pointer :: vlocal(:,:,:,:)
   ! vlocal(n4,n5,n6,nvloc)
   ! local potential in real space, on the augmented fft grid
   ! This component was taken away : CPU time problem !

  real(dp),pointer :: xred(:,:)
   ! xred(3,natom)
   ! reduced coordinates of atoms (dimensionless)

! real(dp),pointer :: ylm(:,:)
   ! ylm(npw,mpsang*mpsang*useylm)
   ! Real spherical harmonics for each G
   ! This component was taken away : CPU time problem !

 end type gs_hamiltonian_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/hdr_type
!! NAME
!! hdr_type
!!
!! FUNCTION
!! It contains all the information needed to write a header for a
!! wf, den or pot file.
!! The structure of the header is explained in the abinis_help.html file.
!! The datatype is considered as an object, to which are attached a whole
!! set of "methods", actually, different subroutines.
!! A few of these subroutines are : hdr_init, hdr_update, hdr_clean,
!! hdr_check, hdr_io, hdr_skip.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type hdr_type
  integer :: bantot        ! total number of bands (sum of nband on all kpts and spins)
  integer :: date          ! starting date
  integer :: headform      ! format of the header
  integer :: intxc,ixc,natom,nkpt,npsp,nspden        ! input variables
  integer :: nspinor,nsppol,nsym,ntypat,occopt        ! input variables
  integer :: pertcase      ! the index of the perturbation, 0 if GS calculation
  integer :: usepaw        ! input variable (0=norm-conserving psps, 1=paw)
  integer :: ngfft(3)      ! input variable

! This record is not a part of the hdr_type, although it is present in the
! header of the files. This is because it depends on the kind of file
! that is written, while all other information does not depend on it.
! It was preferred to let it be initialized or defined outside of hdr_type.
! integer :: fform         ! file descriptor (or file format)

  integer, pointer :: istwfk(:)    ! input variable istwfk(nkpt)
  integer, pointer :: lmn_size(:)  ! lmn_size(npsp) from psps
  integer, pointer :: nband(:)     ! input variable nband(nkpt*nsppol)
  integer, pointer :: npwarr(:)    ! npwarr(nkpt) array holding npw for each k point
  integer, pointer :: pspcod(:)    ! pscod(npsp) from psps
  integer, pointer :: pspdat(:)    ! psdat(npsp) from psps
  integer, pointer :: pspso(:)     ! pspso(npsp) from psps
  integer, pointer :: pspxc(:)     ! pspxc(npsp) from psps
  integer, pointer :: so_typat(:)  ! input variable so_typat(ntypat)
  integer, pointer :: symafm(:)    ! input variable symafm(nsym)
  integer, pointer :: symrel(:,:,:)! input variable symrel(3,3,nsym)
  integer, pointer :: typat(:)     ! input variable typat(natom)

  real(dp) :: ecut                  ! input variable
  real(dp) :: ecutdg                ! input variable (ecut for NC psps, pawecutdg for paw)
  real(dp) :: ecutsm                ! input variable
  real(dp) :: ecut_eff              ! ecut*dilatmx**2 (dilatmx is an input variable)
  real(dp) :: etot,fermie,residm    ! EVOLVING variables
  real(dp) :: qptn(3)               ! the wavevector, in case of a perturbation
  real(dp) :: rprimd(3,3)           ! EVOLVING variables
  real(dp) :: stmbias               ! input variable
  real(dp) :: tphysel               ! input variable
  real(dp) :: tsmear                ! input variable
  real(dp), pointer :: kptns(:,:)   ! input variable kptns(3,nkpt)
  real(dp), pointer :: occ(:)       ! EVOLVING variable occ(bantot)
  real(dp), pointer :: tnons(:,:)   ! input variable tnons(3,nsym)
  real(dp), pointer :: xred(:,:)    ! EVOLVING variable xred(3,natom)
  real(dp), pointer :: zionpsp(:)   ! zionpsp(npsp) from psps
  real(dp), pointer :: znuclpsp(:)  ! znuclpsp(npsp) from psps
                                    ! Note the difference between znucl and znuclpsp !!
  real(dp), pointer :: znucltypat(:)! znucltypat(ntypat) from alchemy

  character(len=6) :: codvsn              ! version of the code
  character(len=132), pointer :: title(:) ! title(npsp) from psps

  type(atmrhoij_type), pointer :: atmrhoij(:) ! EVOLVING variable paw_ij(natom)%rhoij(lmn2_size,nspden), only for paw

!Should make a list of supplementary infos

 end type hdr_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/MPI_type
!! NAME
!! MPI_type
!!
!! FUNCTION
!! The MPI_type structured datatype gather different information
!! about the MPI parallelisation : number of processors,
!! the index of my processor, the different groups of processors, etc ...
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type MPI_type

! Integer scalar

  integer :: paral_compil_kpt
   ! paral_compil_kpt =0 : no -DMPI flag was activated in the compiling procedure
   ! paral_compil_kpt =1 : the -DMPI flag was activated in the compiling procedure

  integer :: paral_compil_fft
   ! paral_compil_fft =0 : no -DMPIFFT flag was activated in the compiling procedure
   ! paral_compil_fft =1 : the -DMPIFFT flag was activated in the compiling procedure

  integer :: paral_compil_mpio
   ! paral_compil_mpio =0 : no -DMPIO flag was activated in the compiling procedure
   ! paral_compil_mpio =1 : the -DMPIO flag was activated in the compiling procedure

  integer :: paral_level
   ! level of parallelization at a moment in the code
   ! level = 1 : level parareel
   ! level = 2 : level nkpt
   ! level = 3 : level FFT

  integer :: paralbd
   ! relevant only if paral_compil_kpt=1 . So, in addition to the kpt parallelization :
   ! paralbd=0 : (no //ization on bands)
   ! paralbd=1 : (//ization on bands)
   ! paralbd>1 : (//ization on blocks of bands)

  integer :: me               ! number of my processor in the group of all processors
  integer :: nproc            ! number of processors
  integer :: me_group         ! number of my processor in my group of kpt
  integer :: nproc_group      ! number of processors in my group of kpt
  integer :: me_fft           ! number of my processor in my group of FFT
  integer :: me_band           ! number of my processor in my group of bands
  integer :: nproc_fft        ! number of processors in my group of FFT
  integer :: master_fft       ! number of master of my fft group (in the world_group)
  integer :: paral_fft        ! set to 1 if the FFT parallelisation is active
  integer :: me_g0            ! if set to 1, means that the current processor is taking care of the G(0 0 0) planewave.
  integer :: num_group_fft    ! number of FFT group of my processor. 0 if my processor is not in a group
  integer :: num_group        ! number of group of my processor. 0 if my processor is not in a group
  integer :: nproc_per_kpt    ! number of processors per kpt
  integer :: world_group      ! number of the group of processor (for MPI communicator)

  integer :: fft_master_group
   ! fft_master_group
   ! group of processors of fft_master_comm
   ! exists only when paral_fft = 1

  integer :: fft_master_comm
   ! fft_master_comm
   ! communicator on master processors
   ! (one processor per fft_group or all processors when paral_fft = 0)

integer :: fft_option_lob
   ! fft_option_lob
   ! option for lob
   ! fft_option_lob=1 : old version of lob
   ! fft_option_lob=2 : new version of lob
   ! exists only when paral_fft = 1


! Integer arrays

  integer, pointer :: fft_group(:)
   ! fft_group(nkpt*nsppol)
   ! tab of groups of processors which treat ffts
   ! exists only when paral_fft = 1

  integer, pointer :: fft_comm(:)
   ! fft_comm(nkpt*nsppol)
   ! tab of communicators of processors which treat ffts of a kpt
   ! exists only when paral_fft = 1

  integer, pointer :: proc_distrb(:,:,:)
   ! proc_distrb(nkpt,mband,nsppol)
   ! number of the processor that will treat
   ! each band in each k point.

  integer, pointer :: kpt_group(:)
   ! kpt_group(nproc_per_kpt)
   ! tab of groups of processors which treat one nkpt/nsppol
   ! exists only when paralbd > 1

  integer, pointer :: kpt_comm(:)
   ! kpt_comm(nproc_per_kpt)
   ! tab of communicators of processors which treat one nkpt/nsppol
   ! exists only when paralbd > 1

  integer, pointer :: kptdstrb(:,:,:)
   ! kptdstrb(me,ineigh,ikptloc)
   ! tab of processors required for mv_3dte.f and berryphase_new.f

  integer, pointer :: kptdstrbi(:,:,:)
   ! same as kptdstrb, but for k-points in the iBZ
   ! required for MPI // of the finite electric field (see vtorho.f)

  integer, pointer :: nplanes_fft(:)
   ! nplanes_fft(nkpt)
   ! number of planes for my proc me_fft
   ! exists only if mpi_enreg%paral_compil_fft==1

   integer, pointer :: ind_fft_planes(:,:)
   ! ind_fft_planes(nkpt,nplanes_fft)
   ! indice of planes for each kpoint for my proc me_fft
   ! exists only if mpi_enreg%paral_compil_fft==1

!BEGIN TF_CHANGES
! Adds for parallelization over perturbations
  integer :: paral_compil_respfn
   ! paral_compil_respfn =0 : no -DMPI flag was activated in the compiling procedure
   ! paral_compil_respfn =1 : the -DMPI flag was activated in the compiling procedure

  integer :: me_respfn           ! number of my processor in my group of perturbations
  integer :: nproc_respfn        ! number of processors in my group of perturbations
  integer :: my_respfn_group     ! my group for calculating perturbations
  integer :: my_respfn_comm      ! my communicator of my_respfn_group
  integer :: respfn_master_group ! groups for masters of respfn_groups
  integer :: respfn_master_comm  ! communicator for masters of respfn_groups
  integer :: ngroup_respfn       ! number of groups for calculating perturbations
  integer :: spaceComm           ! communicator for calculating responsefunction
                                 ! default is MPI_COMM_WORLD but may be changed in 08seqpar/loper3.F90

  integer, pointer :: respfn_group(:) ! groups for calculating perturbations
  integer, pointer :: respfn_comm(:)  ! communicators for respfn_group
!END TF_CHANGES

  integer :: gpara
  ! describes if parallelization over G is selected in the inputfile
  integer :: mgblk
  ! maximal block size for blocks of G
  integer :: gmin
  ! describes the start-G for each processor
  integer :: gmax
  ! describes the end-G for each processor
  integer, pointer :: gmpigroup(:)
  ! one processor group for each k-point
  integer, pointer :: gmpicomm(:)
  ! one communicator for each k-point (each group)
  integer :: ggroup
  ! my group (each k-point is processed in one group)
  integer :: gindex
  ! my group index
  integer :: gmaster
  ! my group master (reserved for future use)
  integer :: gngroup
  ! number of processors in my group

!This is for the bandFFT case
   character :: mode_para
   !If mode_para=='bandFFT', we are in bandFFT mode
   integer :: commcart
   !This is the communicator for the full cartesian array
   integer :: comm_band, comm_fft
   !The communicators over bands and fft respectively
   integer :: me_cart
   !This is the rank of the proc in the full cartesian array
   integer :: dimcart
   !This is the dimension of the cartesian array (2 for 2-dim)
   integer :: nproc_band
   !This is the number of procs on which we distribute bands
   integer, pointer :: sizecart(:)
   !The first dimension is the number of fft processors, the second the number of bands
   integer, pointer :: coords(:)
   !The coordinate of the proc in the cartesian array

! Adds for parareel
  integer :: parareel
   ! parareel = 0 default
   ! parareel = 1 if treats parareel case

! All the following data exist only in the parareel=1 case
  integer :: npara                 ! number of loops on gstate
  integer :: ipara                 ! number of actual internal loop on gstate
  integer :: jpara                 ! number of actual external loop on gstate
  integer :: me_group_para         ! number of my processor in my group of para
  integer :: nproc_group_para      ! number of processors in my group of para
  integer :: num_group_para        ! number of group of my processor. 0 if my processor is not in a group
  integer :: nproc_per_para        ! number of processors per para
  integer :: master_group_para     ! number of the master processor (in the world group) of my group of para

  integer, pointer :: proc_distrb_para(:,:)
   ! proc_distrb_para(npara,nkpt)
   ! exists only when parareel = 1
   ! number of the processor that will treat
   ! each kpt in each para.

  integer, pointer :: kpt_group_para(:)
   ! kpt_group_para(npara)
   ! tab of groups of processors which treat one npara
   ! exists only when parareel = 1

  integer, pointer :: kpt_comm_para(:)
   ! kpt_comm_para(npara)
   ! tab of communicators of processors which treat one npara
   ! exists only when parareel = 1

 end type MPI_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pawang_type
!! NAME
!! pawang_type
!!
!! FUNCTION
!! For PAW, ANGular mesh discretization and related data
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type pawang_type

!Integer scalars

  integer :: angl_size
   ! Dimension of paw angular mesh (angl_size=ntheta*nphi)

  integer :: l_max
   ! Maximum value of angular momentum l+1

  integer :: l_size_max
   ! Maximum value of angular momentum l_size=2*l_max-1

  integer :: lcutd
   !  1+max. index (l) of the moments used in the develop. of densities

  integer :: ngnt
   ! Number of non zero Gaunt coefficients

  integer :: ntheta, nphi
   ! Dimensions of paw angular mesh

  integer :: nsym
   ! Number of symmetry elements in space group

!Integer arrays

  integer, pointer :: gntselect(:,:)
   ! gntselect(l_size_max**2,l_max**2*(l_max**2+1)/2)
   ! Selection rules for Gaunt coefficients
   ! (if gntselect>0, Gaunt coeff. is non-zero)

!Real (double precision) arrays

  real(dp), pointer :: anginit(:,:)
   ! anginit(3,angl_size)
   ! For each point of the angular mesh, gives the coordinates
   ! of the corresponding point on an unitary sphere

  real(dp), pointer :: angwgth(:)
   ! angwgth(angl_size)
   ! For each point of the angular mesh, gives the weight
   ! of the corresponding point on an unitary sphere

  real(dp), pointer :: realgnt(:)
   ! realgnt(ngnt)
   ! Non zero real Gaunt coefficients

  real(dp), pointer :: ylmr(:,:)
   ! ylmr(l_size_max**2,angl_size)
   ! Real Ylm calculated in real space

  real(dp), pointer :: zarot(:,:,:,:)
   !  zarot(l_size_max,l_size_max,l_max,nsym)
   !  Coeffs of the transformation of real spherical
   !  harmonics under the symmetry operations

 end type pawang_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pawfgr_type
!! NAME
!! pawfgr_type
!!
!! FUNCTION
!! For PAW, Fine rectangular GRid parameters and related data
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type pawfgr_type

!Integer scalars

  integer :: mgfft, nfft
   ! Values of mffft and nfft for the fine rectangular grid:
   !   mgfft= max(ngfft(i)) [max. size of 1D FFT grid]
   !   nfft=ngfft1*ngfft2*ngfft3 [number of pts in the FFT box]

  integer :: mgfftc, nfftc
   ! Values of mffft and nfft for the COARSE rectangular grid:
   !   mgfftc= max(ngfftc(i)) [max. size of 1D FFT grid]
   !   nfftc=ngfftc1*ngfftc2*ngfftc3 [number of pts in the FFT box]

  integer :: usefinegrid
   ! Flag: =1 if a double-grid is used to convert spherical data
   !       to Fourier grid. =0 otherwise

  integer :: natom
   ! Number of atoms in the unit cell

!Integer arrays

  integer, pointer :: coatofin(:)
   ! coatofin(nfftc)
   ! Index of the points of the coarse grid on the fine grid

  integer, pointer :: fintocoa(:)
   ! fintocoa(nfft)
   ! Index of the points of the fine grid on the coarse grid
   !  (=0 if the point of the fine grid does not belong to the coarse grid)

  integer :: ngfft(18)
   ! ngfft(1:18)=integer array with FFT box dimensions and other
   ! information on FFTs, for the fine rectangular grid

  integer :: ngfftc(18)
   ! ngfft(1:18)=integer array with FFT box dimensions and other
   ! information on FFTs, for the COARSE rectangular grid

!Real (double precision)

  real(dp) :: gsqcut
   ! Fourier cutoff on G^2 for "large sphere" of radius double
   ! that of the basis sphere corresponding to paw_ecutdg
   ! (concerns the fine rectangular grid)

 end type pawfgr_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pawfgrtab_type
!! NAME
!! pawfgrtab_type
!!
!! FUNCTION
!! For PAW, various arrays giving data related to fine grid for a given atom
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type pawfgrtab_type

!Integer scalars

  integer :: l_size
   ! 1+maximum value of l leading to non zero Gaunt coeffs
   ! for the considered atom type

  integer :: nfgd
   ! Number of Fine rectangular GriD points
   ! in the paw sphere around considered atom

!Integer arrays

  integer, pointer :: ifftsph(:)
   ! ifftsph(nfgd)
   ! Array giving the FFT index (fine grid) of a point in the paw
   ! sphere around considered atom (ifftsph=ix+n1*(iy-1+n2*(iz-1))

!Real (double precision) arrays

  real(dp), pointer :: gylm(:,:)
   ! gylm(nfgd,l_size*l_size)
   ! Gives g(r)*Ylm(r) on the fine rectangular grid
   ! around considered atom

  real(dp), pointer :: gylmgr(:,:,:)
   ! gylmgr(3,nfgd,l_size*l_size)
   ! Gives the gradient of g(r)*Ylm(r) wrt cart. coordinates
   ! on the fine rectangular grid around considered atom

  real(dp), pointer :: rfgd(:,:)
   ! r(3,nfgd)
   ! Gives all R vectors (r-r_atom) on the Fine rectangular GriD
   ! around considered atom

 end type pawfgrtab_type

!!***

!-------------------------------------------------------------------------

!!****t* defs_datatypes/pawrad_type
!! NAME
!! pawrad_type
!!
!! FUNCTION
!! For PAW, RADial mesh discretization and related data
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type pawrad_type

!Integer scalars

  integer :: mesh_size
   ! Dimension of radial mesh

  integer :: mesh_type
   ! Type of mesh
   !     1=regular grid: r(i)=(i-1)*AA
   !     2=logarithmic grid: r(i)=AA*(exp[BB*(i-1)]-1)
   !     3=logarithmic grid: r(i>1)=AA*exp[BB*(i-1)] and r(1)=0
   !     4=logarithmic grid: r(i)=-AA*ln[1-BB*(i-1)] with BB=1/n

!Real (double precision) scalars

  real(dp) :: lstep
   ! Exponential step of the mesh (BB parameter above)
   ! Defined only if mesh type is logarithmic

  real(dp) :: rmax
   ! Max. value of r = rad(mesh_size)

  real(dp) :: rstep
   ! Radial step of the mesh (AA parameter above)

  real(dp) :: stepint
   ! Radial step used to convert any function from the
   ! present grid onto a regular grid in order to
   ! integrate it using trapeze method

!Real (double precision) arrays

  real(dp), pointer :: rad(:)
   ! rad(mesh_size)
   ! Coordinates of all the points of the mesh

  real(dp), pointer :: radfact(:)
   ! radfact(mesh_size)
   ! Factor used to compute radial integrals
   ! Before being integrated on the present mesh,
   ! any function is multiplied by this factor

  real(dp), pointer :: simfact(:)
   ! radfact(mesh_size)
   ! Factor used to compute radial integrals by the a Simpson scheme
   ! Integral[f] = Sum_i [simfact(i)*f(i)]

 end type pawrad_type

!!***

!!****t* defs_datatypes/pawtab_type
!! NAME
!! pawtab_type
!!
!! FUNCTION
!! For PAW, TABulated data initialized at start
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type pawtab_type

!Integer scalars

  integer :: basis_size
   ! Number of elements for the paw nl basis on the considered atom type

!lda+u
  integer :: ij_proj
   ! Number of (i,j) elements for the orbitals on which U acts (PAW+U only)
   ! on the considered atom type (ij_proj=1 (1 projector), 3 (2 projectors)...)

  integer :: ij_size
   ! Number of (i,j) elements for the symetric paw basis
   ! on the considered atom type (ij_size=basis_size*(basis_size+1)/2)

  integer :: l_size
   ! Maximum value of l+1 leading to non zero Gaunt coeffs
   ! (l_size=2*l_max+1)

  integer :: lmn_size
   ! Number of (l,m,n) elements for the paw basis

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for the paw basis

  integer :: lmnmix_sz
   ! lmnmix_sz=number of klmn=(lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix

  integer :: lpawu
   ! lpawu gives l on which U is applied for a given type of atom.

  integer :: nproju
   ! nproju is the number of projectors for orbitals on which paw+u acts.

  integer :: mesh_size
   ! Dimension of radial mesh

  integer :: shape_lambda
   ! Lambda parameter in gaussian shapefunction (shape_type=2)

  integer :: shape_type
   ! Radial shape function type
   ! shape_type=-1 ; g(r)=numeric (read from psp file)
   ! shape_type= 1 ; g(r)=[sin(pi*r/rshp)/(pi*r/rshp)]**2 if r<=rshp, zero if r>rshp
   ! shape_type= 2 ; g(r)=exp[-(r/sigma)**lambda]
   ! shape_type= 3 ; gl(r)=Alpha(1,l)*jl(q(1,l)*r)+Alpha(2,l)*jl(q(2,l)*r) for each l

!lda+u
  integer :: usepawu
   ! usepawu=0 ; do not use PAW+U formalism
   ! usepawu=1 ; use PAW+U formalism (Full localized limit)
   ! usepawu=2 ; use PAW+U formalism (Around Mean Field)

  integer :: usetcore
   ! Flag controling use of pseudized core density (0 if tncore=zero)

  integer :: vlocopt
   ! 0 if Vloc in atomic data is Vbare    (Blochl s formulation)
   ! 1 if Vloc in atomic data is VH(tnzc) (Kresse s formulation)

!Real (double precision) scalars

  real(dp) :: exccore
   ! Exchange-correlation energy for the core density

 real(dp) :: jpaw
   ! jpaw
   ! Value of J parameter for paw+u for a given type.

  real(dp) :: rpaw
   ! Radius of PAW sphere

  real(dp) :: rshp
   ! Compensation charge radius (if r>rshp, g(r)=zero)

  real(dp) :: shape_sigma
   ! Sigma parameter in gaussian shapefunction (shape_type=2)

  real(dp) :: upaw
   ! upaw
   ! Value of U parameter for paw+u for a given type.


!Integer arrays

  integer, pointer :: indklmn(:,:)
   ! indklmn(4,lmn2_size)
   ! Array giving klm, kln, abs(il-jl) and (il+jl) for each klmn=(ilmn,jlmn)
   ! Note: ilmn=(il,im,in) and ilmn<=jlmn

  integer, pointer :: klmntomn(:,:)
   ! klmntomn(4,lmn2_size)
   ! Array giving im, jm ,in, and jn for each klmn=(ilmn,jlmn)
   ! Note: ilmn=(il,im,in) and ilmn<=jlmn
   ! NB: klmntomn is an application and not a bijection

  integer, pointer :: kmix(:)
   ! kmix(lmnmix_sz)
   ! Indirect array selecting the klmn=(lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix

  integer, pointer :: lnproju(:)
   ! lnproju(nproju) gives ln (index for phi) for each projectors on which U acts (PAW+U only)
   ! nproju is 1 or 2 and  is the number of projectors for correlated orbitals

!Real (double precision) arrays

  real(dp), pointer :: coredens(:)
   ! coredens(mesh_size_max)
   ! Gives the core density of the atom

  real(dp), pointer :: dij0(:)
   ! dij0(lmn2_size)
   ! Part of the Dij term (non-local operator) completely
   ! calculated in the atomic data part

  real(dp), pointer :: dltij(:)
   ! dltij(lmn2_size)
   ! Factor used to compute sums over klmn=(ilmn,jlmn)
   ! ((ilmn,ilmn) term has to be added once)
   ! dltij(klmn)=1 if ilmn=jlmn, else dltij(klmn)=2

  real(dp), pointer :: eijkl(:,:)
   ! eijkl(lmn2_size,lmn2_size)
   ! Part of the Dij term (non-local operator) that depends only from
   ! the projected occupation coeffs in the self-consistent loop

  real(dp), pointer :: gnorm(:)
   ! gnorm(l_size)
   ! Give the the normalization factor of each radial shape function

  real(dp), pointer :: phi(:,:)
   ! phi(mesh_size_max, basis_size)
   ! Gives, on the radial grid, the paw all electron wavefunctions

  real(dp), pointer :: phiphj(:,:)
   ! phiphj(mesh_size,ij_size)
   ! Useful product Phi(:,i)*Phi(:,j)

  real(dp), pointer :: phiphjint(:)
   ! phiphjint(ij_proj)
   ! Integration of Phi(:,i)*Phi(:,j) for LDA+U occupation matrix

  real(dp), pointer :: qijl(:,:)
   ! qijl(l_size**2,lmn2_size)
   ! The qijl are the moments of the charge density difference between
   ! the AE and PS partial wave for each channel (i,j). They take part
   ! to the building of the compensation charge

  real(dp), pointer :: rhoij0(:)
   ! rhoij0(lmn2_size)
   ! Initial guess for rhoij

  real(dp), pointer :: shape_alpha(:,:)
   ! shape_alpha(2,l_size)
   ! Alpha_i parameters in Bessel shapefunctions (shape_type=3)

  real(dp), pointer :: shape_q(:,:)
   ! shape_q(2,l_size)
   ! Q_i parameters in Bessel shapefunctions (shape_type=3)

  real(dp), pointer :: shapefunc(:,:)
   ! shapefunc(mesh_size_max,l_size)
   ! Gives the normalized radial shape function for each l component

  real(dp), pointer :: sij(:)
   ! sij(lmn2_size)
   ! Nonlocal part of the overlap operator

  real(dp), pointer :: tcoredens(:)
   ! tcoredens(mesh_size_max)
   ! Gives the pseudo core density of the atom

  real(dp), pointer :: tphi(:,:)
   ! tphi(mesh_size_max,basis_size)
   ! Gives, on the radial grid, the paw atomic pseudowavefunctions

  real(dp), pointer :: tphitphj(:,:)
   ! tphitphj(mesh_size,ij_size)
   ! Useful product tPhi(:,i)*tPhi(:,j)

! lda+u
  real(dp), pointer :: Vee(:,:,:,:)
   ! Screened interaction matrix  Deduced from U and J parameters
   ! computed on the basis of orbitals on which U acts.

 end type pawtab_type

!!***

!-------------------------------------------------------------------------

!!****t* defs_datatypes/paw_an_type
!! NAME
!! paw_an_type
!!
!! FUNCTION
!! For PAW, various arrays given on ANgular mesh or ANgular moments
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type paw_an_type

!Integer scalars

  integer :: angl_size
   ! Dimension of paw angular mesh (angl_size=ntheta*nphi)

  integer :: lm_size
   ! lm_size=(l_size)**2
   ! l is Maximum value of l+1 leading to non zero Gaunt coeffs (l_size=2*l_max+1)

  integer :: mesh_size
   ! Dimension of radial mesh

  integer :: nspden
   ! Number of spin-density components

!Integer arrays

  integer, pointer :: lmselect(:,:)
   ! lmselect(lm_size,nspden)
   ! lmselect(ilm,ispden)=select the non-zero LM-moments of spherical densities n1 and tn1

!Real (double precision) arrays

  real(dp), pointer :: vxc1 (:,:,:)
   ! vxc1(mesh_size,angl_size,nspden)
   ! Gives xc potential inside the sphere
   ! Only if dtset%pawxcdev=0

  real(dp), pointer :: vxc1m (:,:,:)
   ! vxc1m(mesh_size,lm_size,nspden)
   ! Gives (l,m) moments of xc potential inside the sphere
   ! Only if dtset%pawxcdev/=0

  real(dp), pointer :: vxct1 (:,:,:)
   ! vxct1(mesh_size,angl_size,nspden)
   ! Gives xc tild potential inside the sphere
   ! Only if dtset%pawxcdev=0

   real(dp), pointer :: vxct1m (:,:,:)
   ! vxct1m(mesh_size,lm_size,nspden)
   ! Gives (l,m) moments of xc tild potential inside the sphere
   ! Only if dtset%pawxcdev/=0

 end type paw_an_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/paw_ij_type
!! NAME
!! paw_ij_type
!!
!! FUNCTION
!! For PAW, various arrays given on (i,j) (partial waves) channels
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type paw_ij_type

!Integer scalars

  integer :: lmn_size
   ! Number of (l,m,n) elements for the paw basis

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for the paw basis

  integer :: lmnmix_sz
   ! lmnmix_sz=number of (lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix
   !           i.e. number of rhoij elements being mixed during SCF cycle

  integer :: ngrhoij
   ! First dimension of array grhoij

  integer :: nspden
   ! Number of spin-density components

!Integer arrays

  integer, pointer :: nrhoijsel(:)
   ! nrhoijsel(nspden)
   ! Number of non-zero value of rhoij for given spin component
   ! This is the size of rhoijp(:,:) (see below in this datastructure)

  integer, pointer :: kpawmix(:)
   ! kpawmix(lmnmix_sz)
   ! Indirect array selecting the elements of spherical part (rhoij or dij) being mixed during SCF cycle

  integer, pointer :: rhoijselect(:,:)
   ! rhoijselect(lmn2_size,nspden)
   ! Indirect array selecting the non-zero elements of rhoij:
   ! rhoijselect(isel,ispden)=klmn if rhoij(klmn,ispden) is non-zero

!Real (double precision) arrays

  real(dp), pointer :: dij(:,:)
   ! dij(lmn2_size,nspden)
   ! Dij term (non-local operator)

  real(dp), pointer :: grhoij (:,:,:)
   ! grhoij(ngrhoij,lmn2_size,nspden)
   ! Symetrized gradients of Rho_ij wrt xred, strains, ...

  real(dp), pointer :: noccmmp(:,:,:)
   ! noccmmp(2*lpawu+1,2*lpawu+1,nspden)
   ! gives occupation matrix for lda+u (computed in pawdenpot)

  real(dp), pointer :: nocctot(:)
   ! nocctot(nspden)
   ! gives trace of occupation matrix for lda+u (computed in pawdenpot)
   ! for each value of ispden (1 or 2)

  real(dp), pointer :: rhoij (:,:)
   ! rhoij(lmn2_size,nspden)
   ! Symetrized (augmentation) waves occupancies Rho_ij

  real(dp), pointer :: rhoij_ (:,:)
   ! rhoij_(lmn2_size,nspden)
   ! NON symetrized (augmentation) waves occupancies Rho_ij - used only for printing

  real(dp), pointer :: rhoijp (:,:)
   ! rhoijp(lmn2_size,nspden)
   ! Symetrized (augmentation) waves occupancies Rho_ij
   ! in PACKED STORAGE (only non-zero elements are stored)

  real(dp), pointer :: rhoijres (:,:)
   ! rhoijres(lmn2_size,nspden)
   ! Rho_ij residuals

  real(dp), pointer :: veij(:)
   ! veij(lmn2_size)
   ! (i,j) channel that enters into the calculation of
   !  the hartree energy term inside the spheres for
   !  the density n1-ntild1-nhat

 end type paw_ij_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/atmrhoij_type
!! NAME
!! atmrhoij_type
!!
!! FUNCTION
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type atmrhoij_type
  integer :: lmn2_size                 ! First dimension of rhoij
  integer :: nspden                    ! Number of spin-density components
  integer, pointer :: nrhoijsel(:)     ! Number of non-zero values of rhoij(:,ispden)
  integer, pointer :: rhoijselect(:,:) ! (i,j) indexes of non-zero elements of rhoij
  real(dp), pointer :: rhoij(:,:)      ! rhoij(lmn2_size,nspden) augmentation waves occupancies
 end type atmrhoij_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pseudopotential_type
!! NAME
!! pseudopotential_type
!!
!! FUNCTION
!! This structured datatype contains all the information about one
!! norm-conserving pseudopotential, including the description of the local
!! and non-local parts, the different projectors, the non-linear core
!! correction ...
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type pseudopotential_type

! Integer scalars

  integer :: dimekb
   ! Dimension of Ekb
   ! ->Norm conserving : Max. number of Kleinman-Bylander energies
   !                     for each atom type
   !                     dimekb=lnmax (lnmax: see this file)
   ! ->PAW : Max. number of Dij coefficients connecting projectors
   !                     for each atom type
   !                     dimekb=lmnmax*(lmnmax+1)/2 (lmnmax: see this file)

  integer :: lmnmax
   !  If useylm=0, max number of (l,m,n) comp. over all type of psps (lnproj)
   !  If useylm=1, max number of (l,n)   comp. over all type of psps (lmnproj)
   !  If mpspso is 2, lmnmax takes into account the spin-orbit projectors,
   !  so, it is equal to the max of lmnprojso or lnprojso, see pspheader_type

  integer :: lnmax
   !  Max. number of (l,n) components over all type of psps
   !  If mpspso is 2, lmnmax takes into account the spin-orbit projectors,
   !  so, it is equal to the max of lnprojso, see pspheader_type

  integer :: mproj    ! TO BE SUPPRESSED
   ! Maximum number of non-local projectors over all angular momenta
   !  and type of psps
   ! 0 only if all psps are local

  integer :: mpsang
   ! Highest angular momentum of non-local projectors over all type of psps.
   ! shifted by 1 : for all local psps, mpsang=0; for largest s, mpsang=1,
   ! for largest p, mpsang=2; for largest d, mpsang=3; for largest f, mpsang=4
   ! This gives also the number of non-local "channels"

  integer :: mpspso
   ! mpspso is set to 1 if none of the psps is used with a spin-orbit part (that
   !  is, if the user input variable so_typat (new name of pspso) is not equal
   !  to 1 in at least one case
   ! otherwise, it is set to 2

  integer :: mpssoang
   ! Maximum number of channels, including those for treating the spin-orbit coupling
   ! when mpspso=1, mpssoang=mpsang
   ! when mpspso=2, mpssoang=2*mpsang-1

  integer :: mqgrid_ff
   ! Number of points in the reciprocal space grid on which
   ! the radial functions ffspl are specified

  integer :: mqgrid_vl
   ! Number of points in the reciprocal space grid on which
   ! the radial functions vlspl are specified

  integer :: mtypalch
   ! Maximum number of alchemical pseudo atoms. If non-zero,
   ! the mechanism to generate mixing of pseudopotentials is activated

  integer :: npsp
   ! Number of types of pseudopotentials

  integer :: npspalch
   ! Number of types of pseudopotentials use for alchemical purposes

  integer :: ntypat
   ! Number of types of atoms (might be alchemy wrt pseudopotentials)

  integer :: ntypalch
   ! Number of types of alchemical pseudoatoms

  integer :: ntyppure
   ! Number of types of pure pseudoatoms

  integer :: n1xccc
   ! Number of radial points for the description of the pseudo-core charge
   ! (in the framework of the non-linear XC core correction)

  integer :: optnlxccc
   ! Option for the choice of non-linear XC core correction treatment (see the input variable)

  integer :: positron
   ! Option for the choice of type of GS calculation (electron or positron)

  integer :: usepaw
   ! if usepaw=0 , use norm-conserving psps part of the code
   ! is usepaw=1 , use paw part of the code

  integer :: useylm
   ! governs the way the nonlocal operator is to be applied:
   !   1=using Ylm, 0=using Legendre polynomials

! Logical scalars

  logical :: vlspl_recipSpace
   ! governs if vlspl is compute in reciprocal space or in real
   ! space (when available).

! Integer arrays

  integer, pointer :: algalch(:)   ! algalch(ntypalch)
   ! For each type of pseudo atom, the algorithm to mix the pseudopotentials

  integer, pointer :: indlmn(:,:,:)
   ! indlmn(6,lmnmax,ntypat)
   ! For each type of psp,
   ! array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
   !                                or i=lmn (if useylm=1)

  integer, pointer :: pspdat(:)
   ! pspdat(ntypat)
   ! For each type of psp, the date of psp generation, as given by the psp file

  integer, pointer :: pspcod(:)
   ! pspcod(ntypat)
   ! For each type of psp, the format -or code- of psp generation,
   !  as given by the psp file

  integer, pointer :: pspso(:)
   ! pspso(ntypat)
   ! For each type of psp, 1 if no spin-orbit component is taken
   ! into account, 2 if a spin-orbit component is used

  integer, pointer :: pspxc(:)
   ! pspxc(ntypat)
   ! For each type of psp, the XC functional that was used to generate it,
   ! as given by the psp file

! Real (double precision) arrays

  real(dp), pointer :: dnqdq0(:)
   ! dnqdq0(ntypat)
   ! Gives 1/q d(tNcore(q))/dq for q=0 , for each type of PAW psp.
   ! (tNcore(q) = FT of pseudo core density)

  real(dp), pointer :: ekb(:,:)
   ! ekb(dimekb,ntypat*(1-usepaw))
   !  ->NORM-CONSERVING PSPS ONLY:
   !    (Real) Kleinman-Bylander energies (hartree)
   !           for number of basis functions (l,n) (lnmax)
   !           and number of atom types (ntypat)
   ! NOTE (MT) : ekb (norm-conserving) is now diagonal (one dimension
   !             lnmax); it would be easy to give it a second
   !             (symmetric) dimension by putting
   !             dimekb=lnmax*(lnmax+1)/2
   !             in the place of dimekb=lmnmax.

  real(dp), pointer :: ffspl(:,:,:,:)
   ! ffspl(mqgrid_ff,2,lnmax,ntypat)
   ! Gives, on the radial grid, the different non-local projectors,
   ! in both the norm-conserving case, and the PAW case

  real(dp), pointer :: mixalch(:,:)
   ! mixalch(npspalch,ntypalch)
   ! Mixing coefficients to generate alchemical pseudo atoms

  real(dp), pointer :: qgrid_ff(:)
   ! qgrid_ff(mqgrid_ff)
   ! The coordinates of all the points of the radial grid for the nl form factors

  real(dp), pointer :: qgrid_vl(:)
   ! qgrid_vl(mqgrid_vl)
   ! The coordinates of all the points of the radial grid for the local part of psp

  real(dp), pointer :: vlspl(:,:,:)
   ! vlspl(mqgrid_vl,2,ntypat)
   ! Gives, on the radial grid, the local part of each type of psp.

  real(dp), pointer :: dvlspl(:,:,:)
   ! dvlspl(mqgrid_vl,2,ntypat)
   ! Gives, on the radial grid, the first derivative of the local
   ! part of each type of psp (computed when the flag 'vlspl_recipSpace'
   ! is true).

  real(dp), pointer :: ncspl(:,:,:)
   ! ncspl(mqgrid_vl,2,ntypat*usepaw)
   ! PAW only (use xccc1d in NC)
   ! Gives, on the radial grid, the tncore part of each type of psp.

  real(dp), pointer :: xcccrc(:)
   ! xcccrc(ntypat)
   ! Gives the maximum radius of the pseudo-core charge, for each type of psp.

  real(dp), pointer :: xccc1d(:,:,:)
   ! xccc1d(n1xccc*(1-usepaw),6,ntypat)
   ! Norm-conserving only (use ncspl in PAW)
   ! The component xccc1d(n1xccc,1,ntypat) is the pseudo-core charge
   ! for each type of atom, on the radial grid. The components
   ! xccc1d(n1xccc,ideriv,ntypat) give the ideriv-th derivative of the
   ! pseudo-core charge with respect to the radial distance.

  real(dp), pointer :: zionpsp(:)
   ! zionpsp(npsp)
   ! For each pseudopotential, the ionic pseudo-charge
   ! (giving raise to a long-range coulomb potential)

  real(dp), pointer :: ziontypat(:)
   ! ziontypat(ntypat)
   !  For each type of atom (might be alchemy wrt psps), the ionic pseudo-charge
   ! (giving raise to a long-range coulomb potential)

  real(dp), pointer :: znuclpsp(:)
   ! znuclpsp(npsp)
   ! The atomic number of each pseudopotential

  real(dp), pointer :: znucltypat(:)
   ! znucltypat(ntypat)
   ! The atomic number of each type of atom (might be alchemy wrt psps)

! Character arrays

  character(len=fnlen), pointer :: filpsp(:)
   ! filpsp(ntypat)
   ! The filename of the pseudopotential

  character(len=fnlen), pointer :: title(:)
   ! title(ntypat)
   ! The content of first line read from the psp file

 end type pseudopotential_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pspheader_paw_type
!! NAME
!! pspheader_paw_type
!!
!! FUNCTION
!! The pspheader_paw_type structured datatype gather additional information
!! about a PAW pseudopotential file, from its header.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type pspheader_paw_type
  integer :: basis_size    ! Number of elements of the wf basis ((l,n) quantum numbers)
  integer :: l_size        ! Maximum value of l+1 leading to a non zero Gaunt coefficient
  integer :: lmn_size      ! Number of elements of the paw basis
  integer :: mesh_size     ! Dimension of (main) radial mesh
  integer :: pawver        ! Version number of paw psp format
 end type pspheader_paw_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pspheader_type
!! NAME
!! pspheader_type
!!
!! FUNCTION
!! The pspheader_type structured datatype gather different information
!! about a pseudopotential file, from its header.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type pspheader_type
  integer :: nproj(0:3) ! number of scalar projectors for each angular momentum
  integer :: nprojso(3) ! number of spin-orbit projectors for each angular momentum
  integer :: lmax       ! maximum l quantum number (-1 if only local)
                        ! Example : s only       -> lmax=0
                        !           s and p      -> lmax=1
                        !           d only       -> lmax=2
  integer :: pspcod     ! code number of the pseudopotential
  integer :: pspdat     ! date of generation of the pseudopotential
  integer :: pspxc      ! exchange-correlation functional
  integer :: pspso      ! spin-orbit characteristics
  integer :: xccc       ! =0 if no XC core correction, non-zero if XC core correction
  real(dp) :: zionpsp     ! charge of the ion made of core electrons only
  real(dp) :: znuclpsp    ! atomic number of the nuclei
  character(len=fnlen) :: filpsp   ! name of the psp file
  character(len=fnlen) :: title    ! content of first line read from the psp file
  type(pspheader_paw_type) :: pawheader ! only for PAW psps. See above
 end type pspheader_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/results_gs_type
!! NAME
!! results_gs_type
!!
!! FUNCTION
!! This structured datatype contains the results of a GS calculation :
!! energy and its decomposition, forces and their decompositions, stresses
!! and their decompositions
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type results_gs_type

! Integer scalar

  integer :: natom
   ! The number of atoms for this dataset

! Real (double precision) scalars

! All the energies are in Hartree, obtained "per unit cell".
  real(dp) :: eei      ! local pseudopotential energy (Hartree)
  real(dp) :: eeig     ! sum of eigenvalue energy (Hartree)
  real(dp) :: eew      ! Ewald energy (Hartree)
  real(dp) :: ehart    ! Hartree part of total energy (Hartree)
  real(dp) :: eii      ! pseudopotential core-core energy
  real(dp) :: ek       ! kinetic energy (Hartree)
  real(dp) :: enefield ! the term of the energy functional that depends
                       ! explicitely on the electric field
                       ! enefield = -ucvol*E*P
  real(dp) :: enl      ! nonlocal pseudopotential energy (Hartree)
  real(dp) :: entropy  ! entropy (Hartree)
  real(dp) :: enxc     ! exchange-correlation energy (Hartree)
  real(dp) :: enxcdc   ! exchange-correlation double-counting energy (Hartree)
  real(dp) :: epaw     ! PAW spherical energy (Hartree)
  real(dp) :: epawdc   ! PAW spherical double-counting energy (Hartree)
  real(dp) :: etotal   ! total energy (Hartree)
                       ! for fixed occupation numbers (occopt==0,1,or 2):
                       !   etotal=ek+ehart+enxc+eei+eew+eii+enl+PAW_spherical_part
                       ! for varying occupation numbers (occopt>=3):
                       !   etotal=ek+ehart+enxc+eei+eew+eii+enl - tsmear*entropy +PAW_spherical_part
  real(dp) :: fermie   ! Fermi energy (Hartree)
  real(dp) :: residm   ! maximum value for the residual over all bands, all k points,
                       !   and all spins (Hartree or Hartree**2, to be checked !)
  real(dp) :: vxcavg   ! Average of the exchange-correlation energy. The average
                       ! of the local psp pot and the Hartree pot is set to zero (due
                       ! to the usual problem at G=0 for Coulombic system, so vxcavg
                       ! is also the average of the local part of the Hamiltonian

! Real (double precision) arrays

  real(dp), pointer :: fcart(:,:)
   ! fcart(3,natom)
   ! Cartesian forces (Hartree/Bohr)

  real(dp), pointer :: fred(:,:)
   ! fred(3,natom)
   ! Forces in reduced coordinates (Hartree)
   ! Actually, gradient of the total energy with respect
   ! to change of reduced coordinates

  real(dp), pointer :: gresid(:,:)
   ! gresid(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the residual
   ! of the potential

  real(dp), pointer :: grewtn(:,:)
   ! grewtn(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the Ewald energy

  real(dp), pointer :: grhat(:,:)
   ! grhat(3,natom)
   ! PAW ONLY
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the compensation charge

  real(dp), pointer :: grxc(:,:)
   ! grxc(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the XC energy

  real(dp) :: pel(3)
   ! ucvol times the electronic polarization in reduced coordinates

  real(dp) :: strten(6)
   ! Stress tensor in cartesian coordinates (Hartree/Bohr^3)
   ! 6 unique components of this symmetric 3x3 tensor:
   ! Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).

  real(dp), pointer :: synlgr(:,:)
   ! synlgr(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the non-local energy
   ! The "sy" prefix refer to the fact that this gradient has been
   ! symmetrized.

 end type results_gs_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/results_out_type
!! NAME
!! results_out_type
!!
!! FUNCTION
!! This structured datatype contains a subset of the results of a GS
!! calculation, needed to perform the so-called "internal tests", and
!! to perform the timing analysis
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type results_out_type

! Integer scalar

  integer :: natom ! The number of atoms for this dataset

! Integer arrays

  integer, pointer :: npwtot(:)      ! npw(mxnkpt) Full number of plane waves for each
                                     ! k point, computed with the "true" rprimd
                                     ! Not taking into account the decrease due to istwfk
                                     ! Not taking into account the spread of pws on different procs
! Real (double precision) scalars

! All the energies are in Hartree, obtained "per unit cell".
  real(dp) :: etotal  ! total energy (Hartree)

! Real (double precision) arrays

  real(dp) :: acell(3),rprim(3,3),rprimd(3,3),strten(6)
  real(dp), pointer :: fcart(:,:) ! fcart(3,natom) Cartesian forces (Hartree/Bohr)
  real(dp), pointer :: fred(:,:)  ! fred(3,natom)
   ! Forces in reduced coordinates (Hartree)
   ! Actually, gradient of the total energy with respect
   ! to change of reduced coordinates
  real(dp), pointer :: occ(:)     ! occ(mxmband_upper*mxnkpt*mxnsppol)
  real(dp), pointer :: vel(:,:)   ! vel(3,natom)
  real(dp), pointer :: xred(:,:)  ! xred(3,natom)

 end type results_out_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/sigma_parameters
!! NAME
!! sigma_parameters
!!
!! FUNCTION
!! For the GW part of ABINIT, the sigma_parameters structured datatype
!! gather different parameters that characterize the self-energy operator.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type sigma_parameters
  integer :: gwcalctyp
  integer :: npwwfn,npwx,npwc
  integer :: nb
  integer :: nk,nq
  integer :: nkbz,nqbz
  integer :: nop
  integer :: splitsigc
  integer :: ppmodel
  real(dp) :: zcut, deltae
  integer :: nomegasr,nomegasrd,nomegasi
  real(dp) :: omegasrmax,omegasrdmax,omegasimax,omegasimin
  integer :: nkcalc
  integer, pointer :: kcalc(:), minbnd(:), maxbnd(:)
  real(dp), pointer :: xkcalc(:,:)
 end type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/sigma_results
!! NAME
!! sigma_results
!!
!! FUNCTION
!! For the GW part of ABINIT, the sigma_results structured datatype
!! gather the results of sigma (to be explained in more details)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type sigma_results
  integer nb,nk
  integer :: nomega,nomegasrd,nomegasi
  real(dp), pointer :: e0(:,:)
  real(dp), pointer :: en_qp_diago(:,:)
  real(dp), pointer :: vxcme(:,:)
  real(dp), pointer :: sigxme(:,:)
  complex, pointer :: hhartree(:,:,:)
  complex, pointer :: eigvec_qp(:,:,:)
  complex, pointer :: sigcmee0(:,:), dsigmee0(:,:), ze0(:,:)
  complex, pointer :: sigmee(:,:), degw(:,:), egw(:,:)
  real(dp), pointer :: e0gap(:), degwgap(:), egwgap(:)
  complex(dp), pointer :: omegasrd(:,:,:)
  complex, pointer :: sigcmesrd(:,:,:), sigxcmesrd(:,:,:)
  complex(dp), pointer :: omegasi(:)
  complex, pointer :: sigcmesi(:,:,:), sigxcmesi(:,:,:)
  complex(dp), pointer :: omega(:)
  complex, pointer :: sigcme(:,:,:), sigxcme(:,:,:)
  real(dp), pointer :: ame(:,:,:), ak(:,:)
 end type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/vardims_type
!! NAME
!!  vardims_type
!!
!! FUNCTION
!!  Stores dimensions of dataset variables.
!!
!! SOURCE

 type vardims_type

  integer :: mband,ntypat,natom,natsph,nkpt,nkptgw,nshiftk,nsppol,nberry,&
&            nsym,npsp,nconeq,ntypalch,npspalch,nfft,nspden,wfs_dim1,wfs_dim2,&
&            nfreqsus,npw_tiny,nqptdm,norb,ncenter

 end type vardims_type

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/wffile_type
!! NAME
!! wffile_type
!!
!! FUNCTION
!! This structure datatype is a handler for dealing with the IO of a
!! wavefunction file.
!! It contains, among other things, the method of access to the file
!! (standard F90 read/write, or NetCDF call, or MPI IO), the unit number
!! if applicable, the filename, the information on the
!! parallelism, etc ...
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 type wffile_type

! Integer scalar
  integer :: unwff
   ! unwff  unit number of unformatted wavefunction disk file
  integer :: accesswff
   ! Method to access the wavefunction file
   ! =0 if usual Fortran IO routines
   ! =1 if MPI/IO routines (this access method is only available in parallel)
   ! =2 if NetCDF routines (not used yet)
   ! =-1 if usual Fortran IO routines, but only the master node in the parallel case
  integer :: formwff
   ! formwff=format of the eigenvalues
   !   -1 => not used
   !    0 => vector of eigenvalues
   !    1 => hermitian matrix of eigenvalues
  integer ::  kgwff
   ! kgwff  if 1 , read or write kg_k ; if 0, do not care about kg_k
  character(len=fnlen) :: fname
   ! filename (if available)

! In case of MPI parallel use
  integer :: master
   ! master = number of the processor master of the IO procedure when the WffOpen call is issued
  integer :: me
   ! me = number of my processor
  integer :: nproc
   ! nproc = number of processors that will have access to the file
  integer :: spaceComm
   ! spaceComm = space communicator of the IO procedure when the WffOpen call is issued

! In case of MPI/IO : additional information
  integer :: fhwff
   ! fhwff  file handle of unformatted wavefunction disk file (use in MPI/IO only)
  integer :: nbOct_int,nbOct_dp,nbOct_ch,lght_recs
   ! nbOct_int octet number of int value
   ! nbOct_dp octet number of dp value
   ! nbOct_ch octet number of character value
   ! lght_recs length of record

  integer(abinit_offset)  :: offwff,off_recs
   ! offwff  offset position of unformatted wavefunction disk file
   ! off_recs  offset position of start record
   !             (use in parallel)

 end type wffile_type

!----------------------------------------------------------------------

end module defs_datatypes
!!***
