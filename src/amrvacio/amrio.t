!=============================================================================
!> Read the command line arguments passed to amrvac
subroutine readcommandline

use M_kracken
include 'amrvacdef.f'
integer           :: len, ier
logical           :: help

!----------------------------------------------------------------------------
!defaults and usage:
call kracken('cmd','-i amrvac.par -if data -restart -1 -slice 0 -collapse 0 '//&
  '--help .false. -convert .false.')

! Getting the filename
      call retrev('cmd_i',inifile,len,ier)
      call retrev('cmd_if',filenameini,len,ier)
      snapshotini = iget('cmd_restart')
      snapshotnext= snapshotini+1
      slicenext   = iget('cmd_slice')
      collapseNext   = iget('cmd_collapse')
      help = lget('cmd_-help')                    ! --help present?
      convert = lget('cmd_convert')               ! -convert present?

if (mype==0) then
   print*,'-----------------------------------------------------------------------------'
   print*,'-----------------------------------------------------------------------------'
   print*,'|         __  __ ____ ___        _    __  __ ______     ___    ____         |'
   print*,'|        |  \/  |  _ \_ _|      / \  |  \/  |  _ \ \   / / \  / ___|        |'
   print*,'|        | |\/| | |_) | |_____ / _ \ | |\/| | |_) \ \ / / _ \| |            |'
   print*,'|        | |  | |  __/| |_____/ ___ \| |  | |  _ < \ V / ___ \ |___         |'
   print*,'|        |_|  |_|_|  |___|   /_/   \_\_|  |_|_| \_\ \_/_/   \_\____|        |'
   print*,'-----------------------------------------------------------------------------'
   print*,'-----------------------------------------------------------------------------'
end if


if (help) then
   if (mype==0) then 
      print*,'calling example:                                         '
      print*,'./amrvac -i parameterfile -restart 100 -convert -if datamr/data'
      print*,'default parameterfile is amrvac.par                            '
      print*,'Note that parameterfile parameters overwrite the commandline   '
      print*,'-----------------------------------------------------------------------------'
      print*,'                                                         '
   endif
   call comm_finalize
   STOP
endif

end subroutine readcommandline
!=============================================================================
!> Read in the user-supplied parameter-file
subroutine readparameters

include 'amrvacdef.f'

logical :: fileopen
integer :: i, j, k, ifile, iB, isave, iw, level, idim, islice
integer :: nxlone^D
double precision :: dxlone^D

namelist /testlist/   teststr,ixtest1,ixtest2,ixtest3,iwtest,idimtest
namelist /filelist/   filenameout,filenameini,filenamelog, &
                      snapshotini,typefilelog,firstprocess,resetgrid,changeglobals,snapshotnext, &
                      convert,convert_type,dxfiletype,saveprim,primnames, &
                      typeparIO,uselimiter,nwauxio,nocartesian,addmpibarrier, &
                      writew,writelevel,writespshift,endian_swap, &
                      normvar,normt,level_io,level_io_min,level_io_max, &
                      autoconvert,sliceascii,slicenext,collapseNext,collapse_type
namelist /savelist/   tsave,itsave,dtsave,ditsave,nslices,slicedir,slicecoord,collapse,collapseLevel{#IFDEF MAGNETOFRICTION , ditsavemf}
namelist /stoplist/   itmax,tmax,tmaxexact,dtmin,t,it,treset,itreset,residmin,residmax,typeresid{#IFDEF MAGNETOFRICTION , itmaxmf}
namelist /methodlist/ wnames,fileheadout,typeadvance, &
                      ssplitdust,ssplitdivb,ssplitresis,ssplituser,typesourcesplit,&
                      sourceimpl,sourceimplcycle,conduction,TCsaturate,TCphi,ncyclemax,sourceparasts,parastsnu,&
                      dimsplit,typedimsplit,typeaxial,typecoord,&
                      typefull1,typepred1,typelow1,&
                      typelimiter1,mcbeta,typegradlimiter1,&
                      flatcd,flatsh,flatppm,&
                      loglimit,typelimited,useprimitive,typetvdlf, &
                      typetvd,typeentropy,entropycoef,typeaverage, &
                      B0field,Bdip,Bquad,Boct,Busr,divbwave,&
                      typedivbfix,divbdiff,typedivbdiff,compactres,&
                      useprimitiveRel,maxitnr,dmaxvel,typepoly,&
                      tvdlfeps,BnormLF,&
                      smallT,smallp,smallrho,typegrad,typediv,tolernr,absaccnr,&
                      strictnr,fixsmall,strictsmall,strictgetaux,nflatgetaux,&
                      strictzero,nxdiffusehllc,typespherical,&
                      fixprocess,flathllc, &
                      ncool,cmulti,coolmethod,coolcurve,Tfix, &
                      smallrhod, dustzero, dustmethod,dustspecies,dusttemp, &
                      x1ptms,x2ptms,x3ptms,ptmass,tlow,nwtf,neqpartf
namelist /boundlist/  dixB,typeB,typeghostfill,typegridfill,ratebdflux,&
                      internalboundary
namelist /amrlist/    mxnest,nbufferx^D,tol,tolratio,errorestimate, &
                      amr_wavefilter,nxlone^D,dxlone^D,iprob,xprob^L, &
                      skipfinestep,wflags,flags,&
                      restrictprimitive,prolongprimitive,coarsenprimitive, &
                      typeprolonglimit, &
                      amrentropy,logflag,tfixgrid,itfixgrid,ditregrid{#IFDEF STRETCHGRID ,qst}
namelist /paramlist/  time_accurate, courantpar, dtpar, dtdiffpar, dtTCpar,&
                      typecourant, slowsteps, cfrac{#IFDEF MAGNETOFRICTION , cmf_c, cmf_y, cmf_divb}
!----------------------------------------------------------------------------

! defaults for boundary treatments
ratebdflux=one
typeghostfill='linear' 
dixB=2
typeB(1:nw,1:nhiB)='cont'
internalboundary=.false.

! code behavior and testing defaults
addmpibarrier=.false.
teststr=' '
ixtest1=1
ixtest2=1
ixtest3=1
iwtest=1
idimtest=1

! defaults for parameters for optional cooling module (van Marle)
ncool      = 100
cmulti     = 1
coolcurve  = 'DM'
coolmethod = 'explicit2'
cfrac      = 0.1d0
Tfix       = .false.

! defaults for dust
dustzero = .false.
dustmethod = 'Kwok'
dustspecies = 'graphite'
dusttemp = 'constant'

! defaults for parameters for optional pointgrav module (van Marle)
! --> set here mass to zero, coordinates to zero
x1ptms=zero
x2ptms=zero
x3ptms=zero
ptmass=zero

! defaults for specific options
fixprocess=.false.
typegrad='central'
typediv='central'
smallT=-one
smallp=-one
smallrho=-one
smallrhod=-one

! relativistic module defaults
useprimitiveRel=.true.
strictnr=.true.
fixsmall=.false.
strictsmall=.true.
strictgetaux=.false.
nflatgetaux=1
strictzero=.true.
typepoly='meliani'
tolernr=1.0d-13
absaccnr=1.0d-13
maxitnr=100
dmaxvel=1.0d-8
tlow=zero

! defaults for convert behavior
nwauxio=0
nocartesian=.false.
saveprim=.false.
autoconvert=.false.
endian_swap=.false.
convert_type='vtuBCCmpi'
collapse_type='vti'
dxfiletype='lsb'
writew(1:nw)=.true.
writelevel(1:nlevelshi)=.true.
writespshift(1:ndim,1:2)=zero
level_io=-1
level_io_min=1
level_io_max=nlevelshi

! normalization of primitive variables: only for output
! note that normvar(0) is for length
! this scaling is optional, and must be set consistently if used
normvar(0:nw)=one
normt=one

! residual defaults
residmin=-1.0d0
residmax=bigdouble
typeresid='relative'

! AMR related defaults
mxnest=1
nbufferx^D=0;
tol(1:nlevelshi)=0.1d0
tolratio(1:nlevelshi)=1.0d0/8.0d0
typegridfill='linear'
amrentropy=.false.
restrictprimitive=.false.
coarsenprimitive=.false.
prolongprimitive=.false.
typeprolonglimit='default'
errorestimate=3
flags(1:nflag_)=0
wflags(1:nflag_)=zero
flags(nflag_)=1
flags(1)=1
wflags(1)=one
logflag(1:nw)=.false.
amr_wavefilter(1:nlevelshi)=1.0d-2
skipfinestep=.false.
tfixgrid=bigdouble
itfixgrid=biginteger
ditregrid=1
{#IFDEF STRETCHGRID
qst=bigdouble
}

! MHD specific defaults
B0field=.false.
Bdip=zero
Bquad=zero
Boct=zero
Busr=zero
{#IFNDEF GLM
typedivbfix='linde'\}
{#IFDEF GLM
typedivbfix='glm1'\}
divbdiff=0.5d0
typedivbdiff='all'
compactres=.false.
divbwave=.true.


! IO defaults
itmax=biginteger
{#IFDEF MAGNETOFRICTION
itmaxmf=0
ditsavemf=20000
}
tmax=bigdouble
tmaxexact=.true.
dtmin=1.0d-10
typeparIO=0
nslices=0
collapse=.false.
collapseLevel=1
sliceascii=.false.
do ifile=1,nfile
   do isave=1,nsavehi
      tsave(isave,ifile)=bigdouble   ! t  of saves into the output files
      itsave(isave,ifile)=biginteger ! it of saves into the output files
   end do
   dtsave(ifile)=bigdouble           ! time between saves
   ditsave(ifile)=biginteger         ! timesteps between saves
   isavet(ifile)=1                   ! index for saves by t
   isaveit(ifile)=1                  ! index for saves by it
end do
typefilelog='default'
fileheadout = 'AMRVAC'
nwtf=0
neqpartf=0

! defaults for input 
firstprocess=.false.
resetgrid=.false.
changeglobals=.false.
treset=.false.
itreset=.false.
filenameout='data'
filenamelog='amrvac'

! Defaults for discretization methods
typeaverage='default'
tvdlfeps=one
{#IFDEF FCT
BnormLF=.false.
}
{#IFNDEF FCT
BnormLF=.true.
}
nxdiffusehllc=0
flathllc=.false.
typeaxial='slab'
typecoord='default'
typespherical=1
slowsteps=-1
courantpar=0.8d0
typecourant='minimum'
dimsplit=.false.
typedimsplit='default'
typelimited='predictor'
mcbeta=1.4d0
useprimitive=.true.
typetvd='roe'
typetvdlf='cmaxmean'
sourceimpl=.false.
sourceimplcycle=.false.
conduction=.false.
TCsaturate=.false.
TCphi=1.d0
bcphys=.true.
ncyclemax=1000
sourceparasts=.false.
parastsnu=0.001d0
ssplitdust=.false.
ssplitdivb=.false.
{^IFMHDPHYS
ssplitdivb=.true.
}
ssplitresis=.false.
ssplituser=.false.
typeadvance='twostep'
do level=1,nlevelshi
   typefull1(level)='tvdlf'
   typepred1(level)='default'
   typelow1(level)='default'
   typelimiter1(level)='minmod'
   typegradlimiter1(level)='minmod'
enddo
flatcd=.false.
flatsh=.false.
flatppm=.true.
typesourcesplit='sfs'
loglimit(1:nw)=.false.
do iw=1,nw
   typeentropy(iw)='nul'      ! Entropy fix type
end do
dtdiffpar=0.5d0
dtTCpar=0.9d0
dtpar=-1.d0
time_accurate=.true.
{#IFDEF MAGNETOFRICTION
cmf_c=0.3d0
cmf_y=0.2d0
cmf_divb=0.1d0
}

! problem setup defaults
dxlone^D=zero;
nxlone^D=0;
iprob=1

! end defaults

! Initialize Kronecker delta, and Levi-Civita tensor
do i=1,3
   do j=1,3
      if(i==j)then
         kr(i,j)=1
      else
         kr(i,j)=0
      endif
      do k=1,3
         if(i==j.or.j==k.or.k==i)then
            lvc(i,j,k)=0
         else if(i+1==j.or.i-2==j)then
            lvc(i,j,k)=1
         else
            lvc(i,j,k)=-1
         endif
      enddo
   enddo
enddo

! MPI reads from a file
open(unitpar,file=inifile,status='old')

! Start reading from standard input
read(unitpar,testlist)

oktest=index(teststr,'readparameters')>=1
if(oktest) write(unitterm,*)'ReadParameters'
if(oktest) write(unitterm,testlist)

primnames='default'

read(unitpar,filelist)

if(TRIM(primnames)=='default'.and.mype==0) write(uniterr,*) &
   'Warning in ReadParameters: primnames not given!'

if(firstprocess .and. snapshotini<0) &
  call mpistop("Please restart from a snapshot when firstprocess=T")
if(convert .and. snapshotini<0) then
  convert=.false.
  write(uniterr,*) 'Warning in ReadParameters: ',&
        'Please change convert to .false. when start a new run!'
end if
if(convert) autoconvert=.false.

read(unitpar,savelist)
do ifile=1,nfile
   if(dtsave(ifile)<bigdouble/2.and.oktest) &
      write(unitterm,'(" DTSAVE  for file",i2," =",g12.5)') &
         ifile,dtsave(ifile)
   if(ditsave(ifile)<biginteger.and.oktest) &
      write(unitterm,'(" DITSAVE for file",i2," =",i10)') &
         ifile,ditsave(ifile)
   if(tsave(1,ifile)==bigdouble.and.itsave(1,ifile)==biginteger.and. &
       dtsave(ifile)==bigdouble.and.ditsave(ifile)==biginteger.and.mype==0) &
       write(uniterr,*)'Warning in ReadParameters: ', &
                       'No save condition for file ',ifile
enddo

do islice=1,nslices
   if(slicedir(islice) > ndim) &
   write(uniterr,*)'Warning in ReadParameters: ', &
        'Slice ', islice,' direction',slicedir(islice),'larger than ndim=',ndim
   if(slicedir(islice) < 1) &
   write(uniterr,*)'Warning in ReadParameters: ', &
        'Slice ', islice,' direction',slicedir(islice),'too small, should be [',1,ndim,']'
end do

read(unitpar,stoplist)
if(oktest)then
   if(itmax<biginteger)       write(unitterm,*) 'ITMAX=',itmax
   if(tmax<bigdouble)         write(unitterm,*) 'TMAX=',tmax
   write(unitterm,*)'tmaxexact=',tmaxexact
   if(dtmin>smalldouble)      write(unitterm,*) 'DTMIN=',dtmin
endif

if(itmax==biginteger .and. tmax==bigdouble.and.mype==0) write(uniterr,*) &
   'Warning in ReadParameters: itmax or tmax not given!'

if(residmin>=zero) then
   if (mype==0) write(unitterm,*)"SS computation with input value residmin"
   if(residmin<=smalldouble) call mpistop("Provide value for residual above smalldouble")
end if

wnames='default'

read(unitpar,methodlist)

if(compactres)then
 if(typeaxial/='slab') call mpistop("compactres for MHD only in cartesian case")
endif

if(TRIM(wnames)=='default') call mpistop("Provide wnames and restart code")
wnameslog=wnames

do level=1,nlevelshi
   !if(typefull1(level)=='tvdlf1'.and.typeadvance=='twostep') &
   !   call mpistop(" tvdlf1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='hll1'.and.typeadvance=='twostep') &
   !   call mpistop(" hll1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='hllc1'.and.typeadvance=='twostep') &
   !   call mpistop(" hllc1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='hllcd1'.and.typeadvance=='twostep') &
   !   call mpistop(" hllcd1 is onestep method, reset typeadvance=onestep!")
   !if(typefull1(level)=='tvdmu1'.and.typeadvance=='twostep') &
   !   call mpistop(" tvdmu1 is onestep method, reset typeadvance=onestep!")
   if(typefull1(level)=='tvd'.and.typeadvance=='twostep') &
      call mpistop(" tvd is onestep method, reset typeadvance=onestep!")
   if(typefull1(level)=='tvd1'.and.typeadvance=='twostep') &
      call mpistop(" tvd1 is onestep method, reset typeadvance=onestep!")
   if(typefull1(level)=='tvd'.or.typefull1(level)=='tvd1')then 
      if(mype==0.and.(.not.dimsplit)) write(unitterm,*) &
         'Warning: setting dimsplit=T for tvd, as used for level=',level
      dimsplit=.true.
   endif

   if (typepred1(level)=='default') then
      select case (typefull1(level))
       case ('cd')
         typepred1(level)='cd'
      case ('cd4')
         typepred1(level)='cd4' 
      case ('fd')
         typepred1(level)='fd'
      case ('tvdlf','tvdmu')
         typepred1(level)='hancock'
      case ('hll')
         typepred1(level)='hll'
      case ('hllc')
         typepred1(level)='hllc'
      case ('hllcd')
         typepred1(level)='hllcd'
      case ('hlld')
         typepred1(level)='hlld'
      case ('hlldd')
         typepred1(level)='hlldd'
      case ('tvdlf1','tvdmu1','tvd1','tvd','hll1','hllc1', &
            'hlld1','hllcd1','hlldd1','nul','source')
         typepred1(level)='nul'
      case default
         call mpistop("No default predictor for this full step")
      end select
   end if
end do

select case (typeadvance)
case ("onestep")
   nstep=1
case ("twostep")
   nstep=2
case ("threestep")
   nstep=3
case ("fourstep","rk4","jameson","ssprk43")
   nstep=4
case ("ssprk54")
   nstep=5
case default
   call mpistop("Unknown typeadvance")
end select


do level=1,nlevelshi
  if (typelow1(level)=='default') then
   select case (typefull1(level))
   case ('cd')
      typelow1(level)='cd'
   case ('cd4')
      typelow1(level)='cd4'
   case ('fd')
      typelow1(level)='fd'
   case ('hancock')
      typelow1(level)='hancock'
   case ('tvdlf','tvdlf1','tvdmu','tvdmu1','tvd1','tvd')
      typelow1(level)='tvdlf1'
   case ('hll','hll1')
      typelow1(level)='hll1'
   case ('hllc','hllc1')
      typelow1(level)='hllc1'
   case ('hllcd','hllcd1')
      typelow1(level)='hllcd1'
   case ('hlld','hlld1')
      typelow1(level)='hlld1'
   case ('hlldd','hlldd1')
      typelow1(level)='hlldd1'
   case ('nul')
      typelow1(level)='nul'
   case ('source')
      typelow1(level)='source'
   case default
      call mpistop("No default typelow for this full step")
   end select
  end if
enddo

! Harmonize the parameters for dimensional splitting and source splitting
if(typedimsplit   =='default'.and.     dimsplit)   typedimsplit='xyyx'
if(typedimsplit   =='default'.and..not.dimsplit)   typedimsplit='unsplit'
dimsplit   = typedimsplit   /='unsplit'


if (typeaxial=="slab") then
   slab=.true.
else
   slab=.false.
end if

if (typeaxial=='spherical') then
   if (dimsplit) then
      if(mype==0)print *,'Warning: spherical symmetry needs dimsplit=F, resetting'
      dimsplit=.false.
   end if
end if

if (typecoord=='default') then
   typecoord = typeaxial
end if

if (ndim==1) dimsplit=.false.
if (.not.dimsplit.and.ndim>1) then
   select case (typeadvance)
   case ("ssprk54","ssprk43","fourstep", "rk4", "threestep", "twostep")
      ! Runge-Kutta needs predictor
      typelimited="predictor"
      if(mype==0)print *,'typelimited to predictor for RK'
   end select
end if

{#IFDEF GLM
if (typedivbfix/='glm1' .and. typedivbfix/='glm2' .and. typedivbfix/='glm3') &
     call mpistop('using GLM, so typedivbfix should be either glm1 or glm2 (or glm3 with no additional sources)')
if (ssplitdivb .eqv. .false.) &
     call mpistop('GLM needs ssplitdivb = .true.')
\}
{#IFNDEF GLM
if(typedivbfix=='glm1' .or. typedivbfix=='glm2' .or. typedivbfix=='glm3') then
  call mpistop('Not using GLM, to enable GLM, change definitions.h. '//&
      'Alternatives for typedivbfix are powel, janhunen, linde or none')
end if
\}
if (B0field) then
   if(mype==0)print *,'B0+B1 split for MHD'
   if (.not.typephys=='mhd') call mpistop("B0+B1 split for MHD only")
end if

if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    .and.(flatsh.and.typephys=='rho')) then
    call mpistop(" PPM with flatsh=.true. can not be used with typephys='rho'!")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    .and.(flatsh.and.typephys=='hdadiab')) then
     call mpistop(" PPM with flatsh=.true. can not be used with typephys='hdadiab'!")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    .and.(flatcd.and.typephys=='hdadiab')) then
     call mpistop(" PPM with flatcd=.true. can not be used with typephys='hdadiab'!")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    .and.(flatsh.and..not.useprimitive)) then
     call mpistop(" PPM with flatsh=.true. needs useprimitive=T!")
end if
if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    .and.(flatcd.and..not.useprimitive)) then
     call mpistop(" PPM with flatcd=.true. needs useprimitive=T!")
end if

if(oktest) write(unitterm,methodlist)

read(unitpar,boundlist)
if(oktest) write(unitterm,boundlist)
do idim=1,ndim
   periodB(idim)=(any(typeB(:,2*idim-1:2*idim)=='periodic'))
   aperiodB(idim)=(any(typeB(:,2*idim-1:2*idim)=='aperiodic'))
   if (periodB(idim).or.aperiodB(idim)) then
      do iw=1,nw
         if (typeB(iw,2*idim-1) .ne. typeB(iw,2*idim)) &
              call mpistop("Wrong counterpart in periodic boundary")
         if (typeB(iw,2*idim-1) /= 'periodic' .and. typeB(iw,2*idim-1) /= 'aperiodic') &
              call mpistop("Each dimension should either have all &
              or no variables periodic, some can be aperiodic")
      end do
   end if
end do

if (any(typelimiter1(1:nlevelshi)=='ppm').and.(dixB<4)) then
    call mpistop(" PPM works only with dixB>=4 !")
end if

if (any(typelimiter1(1:nlevelshi)=='mp5') .and. (dixB<3)) then
 call mpistop("mp5 needs at at least 3 ghost cells! Set dixB=3 in boundlist.")
end if

if(oktest) write(unitterm,*)"periodB: ",periodB(1:ndim)

read(unitpar,amrlist)

select case (typeaxial)
{^NOONED
case ("spherical")
   xprob^LIM^DE=xprob^LIM^DE*two*dpi;
\}
case ("cylindrical")
   {if (^D==^PHI) then
      xprob^LIM^D=xprob^LIM^D*two*dpi;
   end if\}
end select

{#IFDEF STRETCHGRID
!if (mxnest>1) call mpistop("No refinement possible with a loggrid")
if (typeaxial=='slab') call mpistop("Cartesian log grid not implemented")
if (xprobmin1<=0) call mpistop("xprobmin1 must be positive in a stretched grid")
if (qst/=bigdouble) then
   xprobmax1=xprobmin1*qst**nxlone1
   logG=2.d0*(qst-1.d0)/(qst+1.d0)
   if(mype==0) write(*,*) 'xprobmax1 is computed for given nxlone1 and qst:', xprobmax1
else if (qst==bigdouble .and. xprobmax1/=bigdouble) then
   qst=(xprobmax1/xprobmin1)**(1.d0/dble(nxlone1))
   logG=2.d0*(qst-1.d0)/(qst+1.d0)
   if(mype==0) write(*,*) 'logG and qst computed from xprobmax1: ', logG, qst
end if
if(mype==0) write(unitterm,*)'Stretched grid level one minimal dx1:',xprobmin1/(one-half*logG)*logG,&
         ' maximal dx1:',xprobmax1/qst/(one-half*logG)*logG
if(mype==0) write(unitterm,*)'first cell center at x1:',xprobmin1*qst**(-dixB)/(one-half*logG)
}

{
if(nxlone^D>1 .and. mod(nxlone^D,2)==0)then
   dxlone^D=(xprobmax^D-xprobmin^D)/dble(nxlone^D)
   if (mype==0) then
      write(unitterm,*)'Using ',nxlone^D,' cells in dimension ',^D
      write(unitterm,*)'level one dx(',^D,')=',dxlone^D
   end if
end if
\}
if(oktest) write(unitterm,amrlist)
if({dxlone^D*}<smalldouble)then
   write(unitterm,*)'Wrong value(s) for level one dx:',dxlone^D
   call mpistop("Reset nxlone or dxlone!")
endif
^D&dx(^D,1)=dxlone^D;
if(mxnest>nlevelshi.or.mxnest<1)then
   write(unitterm,*)'Error: mxnest',mxnest,'>nlevelshi ',nlevelshi
   call mpistop("Reset nlevelshi and recompile!")
endif

if (flags(nflag_)>nw) then
   write(unitterm,*)'Error: flags(nw+1)=',flags(nw+1),'>nw ',nw
   call mpistop("Reset flags(nw+1)!")
end if
if (flags(nflag_)==0) errorestimate=0
if (flags(nflag_)<0) then
   if (mype==0) then
      write(unitterm,*) "flags(",nflag_,") can not be negative"
      call mpistop("")
   end if
end if
select case (errorestimate)
case (0)
   if (mype==0) write(unitterm,*)"Error estimation is user defined"
case (1)
   if (mype==0) write(unitterm,*)"Error estimation is richardson procedure"
case (2)
   if (mype==0) write(unitterm,*)"Error estimation is relative error"
case (3)
   if (mype==0) write(unitterm,*)"Error estimation is Lohner's scheme"
case (4)
   if (mype==0) write(unitterm,*)"Error estimation is Lohner's original scheme"
case default
   call mpistop("Unknown error estimator, change errorestimate")
end select
if (B0field.and.errorestimate==1) then
   call mpistop("No Richardson procedure in combination with B0")
end if

if (tfixgrid<bigdouble/2.0d0) then
   if(mype==0)print*,'Warning, at time=',tfixgrid,'the grid will be fixed'
end if
if (itfixgrid<biginteger/2) then
   if(mype==0)print*,'Warning, at iteration=',itfixgrid,'the grid will be fixed'
end if
if (ditregrid>1) then
   if(mype==0)print*,'Note, Grid is reconstructed once every',ditregrid,'iterations'
end if


do islice=1,nslices
select case(slicedir(islice))
{case(^D)
   if(slicecoord(islice)<xprobmin^D.or.slicecoord(islice)>xprobmax^D) &
   write(uniterr,*)'Warning in ReadParameters: ', &
        'Slice ', islice, ' coordinate',slicecoord(islice),'out of bounds for dimension ',slicedir(islice)
\}
end select
end do

read(unitpar,paramlist)
if(oktest) write(unitterm,paramlist)
if (dtpar>zero) time_accurate=.true.

if(.not.time_accurate) then
  if(residmin<=smalldouble .or. residmax==bigdouble) then
   if (mype==0) write(unitterm,*)"Non time_accurate SS computation needs values residmin and residmax"
   call mpistop("Provide values for residual bounds in stoplist")
  end if
end if

! Warn when too few blocks at start of simulation 
if (mype.eq.0 .and. snapshotini.eq.-1 .and. {^D& floor(dble(nxlone^D)/(dble(ixGhi^D)-2.0d0*dble(dixB))) |*} .lt. npe) then
   call mpistop('Need at least as many blocks on level 1 as cores to initialize!')
end if


close(unitpar)

if (mype==0) then
   print*,'Reading from inifile: ', trim(inifile)
   print*,'snapshotini         : ', snapshotini
   print*,'slicenext           : ', slicenext
   print*,'collapsenext        : ', collapsenext
   print*,'Filenameini         : ', trim(filenameini)
   print*,'Converting?         : ', convert
   print*,'                                                                '
endif

if(oktest)write(unitterm,*)'End readparameters'

end subroutine readparameters
!=============================================================================
subroutine saveamrfile(ifile)

! following specific for Intel compiler and use on VIC3 with MPT
!DEC$ ATTRIBUTES NOINLINE :: write_snapshot

include 'amrvacdef.f'
integer:: ifile
!-----------------------------------------------------------------------------
select case (ifile)
case (fileout_)
   if(endian_swap) typeparIO=-1
   if (typeparIO==1)then
     call write_snapshot
   else if(typeparIO==0) then
     call write_snapshot_nopar
   else if(typeparIO==-1) then
     call write_snapshot_noparf
   endif
{#IFDEF BOUNDARYDRIVER
   call write_boundary
}
!opedit: now we can also convert directly and will when autoconvert is set in inifile: 
   if (autoconvert) call generate_plotfile
{#IFDEF PARTICLES
   call write_particles_snapshot
}
case (fileslice_)
   call write_slice
case (filecollapse_)
   call write_collapsed
case (filelog_)
   select case (typefilelog)
   case ('default')
      call printlog_default
   case ('special')
      call printlog_special
   case default
      call mpistop("Error in SaveFile: Unknown typefilelog")
   end select
case (fileanalysis_)
  call write_analysis
case default
   write(*,*) 'No save method is defined for ifile=',ifile
   call mpistop("")
end select

! opedit: Flush stdout and stderr from time to time.
flush(unit=unitterm)

end subroutine saveamrfile
!=============================================================================
subroutine write_snapshot
use mod_forest
include 'amrvacdef.f'

{#IFDEF TRANSFORMW
double precision, allocatable :: wtf(:^D&,:)
double precision :: eqpar_tf(neqpartf)
integer :: file_handle_tf
character(len=80) :: filenametf
}
integer :: file_handle, amode, igrid, Morton_no, iwrite
integer :: nx^D
integer(kind=MPI_OFFSET_KIND) :: offset
integer, dimension(MPI_STATUS_SIZE) :: status
character(len=80) :: filename, line
logical, save :: firstsnapshot=.true.
!-----------------------------------------------------------------------------
if (firstsnapshot) then
   snapshot=snapshotnext
   firstsnapshot=.false.
end if

if (snapshot >= 10000) then
   if (mype==0) then
      write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
      write(*,*) "overwriting first frames"
   end if
   snapshot=0
end if

! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"

if(mype==0) then
   open(unit=unitsnapshot,file=filename,status='replace')
   close(unit=unitsnapshot, status='delete')
end if
call MPI_BARRIER(icomm,ierrmpi)

amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
call MPI_FILE_OPEN(icomm,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)

{#IFDEF TRANSFORMW
if(nwtf>0 .and. neqpartf>0) then
  write(filenametf,"(a,i4.4,a)") TRIM(filenameout),snapshot,"tf.dat"
  if(mype==0) then
     open(unit=unitsnapshot,file=filenametf,status='replace')
     close(unit=unitsnapshot)
  end if
  amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
  call MPI_FILE_OPEN(icomm,filenametf,amode,MPI_INFO_NULL,file_handle_tf,ierrmpi)
  allocate(wtf(ixG^T,1:nwtf))
endif
}

iwrite=0
do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine, 
      ! which might be used in getaux
      saveigrid=igrid
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      if (.not.slab) mygeo => pgeo(igrid)
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
   endif
   iwrite=iwrite+1
{#IFDEF TRANSFORMW
   if(nwtf>0 .and. neqpartf>0) then
     call transformw_usr(pw(igrid)%w,wtf,eqpar_tf,ixG^LL,ixM^LL)
     offset=int(size_block_io_tf,kind=MPI_OFFSET_KIND) &
            *int(Morton_no-1,kind=MPI_OFFSET_KIND)
     call MPI_FILE_WRITE_AT(file_handle_tf,offset,wtf,1,type_block_io_tf, &
                             status,ierrmpi)     
   endif
}
{#IFDEF EVOLVINGBOUNDARY
   nphyboundblock=sum(sfc_phybound(1:Morton_no-1))
   offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
          *int(Morton_no-1-nphyboundblock,kind=MPI_OFFSET_KIND) + &
          int(size_block,kind=MPI_OFFSET_KIND) &
          *int(nphyboundblock,kind=MPI_OFFSET_KIND)
   if (sfc_phybound(Morton_no)==1) then
      call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,type_block, &
                              status,ierrmpi)
   else
      call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,&
                              type_block_io,status,ierrmpi)
   end if
}{#IFNDEF EVOLVINGBOUNDARY
   offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
          *int(Morton_no-1,kind=MPI_OFFSET_KIND)
   call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,type_block_io, &
                           status,ierrmpi)
}
end do


call MPI_FILE_CLOSE(file_handle,ierrmpi)
{#IFDEF TRANSFORMW
if(nwtf>0 .and. neqpartf>0) then
  call MPI_FILE_CLOSE(file_handle_tf,ierrmpi)
endif
}
if (mype==0) then
   amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
   call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL, &
                      file_handle,ierrmpi)

   call write_forest(file_handle)

   {nx^D=ixMhi^D-ixMlo^D+1
   call MPI_FILE_WRITE(file_handle,nx^D,1,MPI_INTEGER,status,ierrmpi)\}
   call MPI_FILE_WRITE(file_handle,eqpar,neqpar+nspecialpar, &
                       MPI_DOUBLE_PRECISION,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nleafs,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,levmax,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndim,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndir,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nw,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,neqpar+nspecialpar,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,it,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)
{#IFDEF EVOLVINGBOUNDARY
   nphyboundblock=sum(sfc_phybound)
   call MPI_FILE_WRITE(file_handle,nphyboundblock,1,MPI_INTEGER,status,ierrmpi)
}
   call MPI_FILE_CLOSE(file_handle,ierrmpi)
end if
{#IFDEF TRANSFORMW
if (mype==0 .and. nwtf>0 .and. neqpartf>0) then
   amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
   call MPI_FILE_OPEN(MPI_COMM_SELF,filenametf,amode,MPI_INFO_NULL, &
                      file_handle_tf,ierrmpi)

   call write_forest(file_handle_tf)

   {nx^D=ixMhi^D-ixMlo^D+1
   call MPI_FILE_WRITE(file_handle_tf,nx^D,1,MPI_INTEGER,status,ierrmpi)\}
   call MPI_FILE_WRITE(file_handle_tf,eqpar_tf,neqpartf, &
                       MPI_DOUBLE_PRECISION,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,nleafs,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,levmax,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,ndim,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,ndir,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,nwtf,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,neqpartf,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,it,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)

   call MPI_FILE_CLOSE(file_handle_tf,ierrmpi)
endif
}
snapshot=snapshot+1

end subroutine write_snapshot
!=============================================================================
subroutine write_snapshot_nopar
use mod_forest
include 'amrvacdef.f'

{#IFDEF TRANSFORMW
double precision, allocatable :: wtf(:^D&,:)
double precision :: eqpar_tf(neqpartf)
integer :: file_handle_tf
character(len=80) :: filenametf
}
integer :: file_handle, amode, igrid, Morton_no, iwrite
integer :: nx^D

integer(kind=MPI_OFFSET_KIND) :: offset

integer, allocatable :: iorecvstatus(:,:),ioastatus(:,:)
integer, allocatable :: igrecvstatus(:,:)
integer, allocatable :: igrid_recv(:) 

integer, dimension(MPI_STATUS_SIZE) :: status

integer  :: ipe,insend,inrecv,nrecv,nwrite
character(len=80) :: filename, line
logical, save :: firstsnapshot=.true.
!-----------------------------------------------------------------------------
call MPI_BARRIER(icomm,ierrmpi)

if (firstsnapshot) then
   snapshot=snapshotnext
   firstsnapshot=.false.
end if

if (snapshot >= 10000) then
   if (mype==0) then
      write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
      write(*,*) "overwriting first frames"
   end if
   snapshot=0
end if

nrecv=0
inrecv=0
nwrite=0
insend=0
iwrite=0

if (mype /= 0) then
 do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    itag=Morton_no
    insend=insend+1
    if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine, 
      ! which might be used in getaux
      saveigrid=igrid
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      mygeo =>pgeo(igrid)
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
    endif
    call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
    call MPI_SEND(pw(igrid)%w,1,type_block_io, 0,itag,icomm,ierrmpi)
 end do
else 
 ! mype==0
 nwrite=(Morton_stop(0)-Morton_start(0)+1)

 ! master processor writes out
 write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"

 open(unit=unitsnapshot,file=filename,status='replace')
 close(unitsnapshot, status='delete')

 amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
 call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)

{#IFDEF TRANSFORMW
 if(nwtf>0 .and. neqpartf>0) then
   write(filenametf,"(a,i4.4,a)") TRIM(filenameout),snapshot,"tf.dat"
   call MPI_FILE_OPEN(MPI_COMM_SELF,filenametf,amode,MPI_INFO_NULL,file_handle_tf,ierrmpi)
   allocate(wtf(ixG^T,1:nwtf))

   open(unit=unitsnapshot,file=filenametf,status='replace')
   close(unit=unitsnapshot)
 endif
}

 ! writing his local data first
 do Morton_no=Morton_start(0),Morton_stop(0)
   igrid=sfc_to_igrid(Morton_no)
   iwrite=iwrite+1
   if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine,
      ! which might be used in getaux
      saveigrid=igrid
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      mygeo =>pgeo(igrid)
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
   endif
{#IFDEF TRANSFORMW
   if(nwtf>0 .and. neqpartf>0) then
     call transformw_usr(pw(igrid)%w,wtf,eqpar_tf,ixG^LL,ixM^LL)
     offset=int(size_block_io_tf,kind=MPI_OFFSET_KIND) &
            *int(Morton_no-1,kind=MPI_OFFSET_KIND)
     call MPI_FILE_WRITE_AT(file_handle_tf,offset,wtf,1,type_block_io_tf, &
                             status,ierrmpi)
   endif
}
   offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
          *int(Morton_no-1,kind=MPI_OFFSET_KIND)
   call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,type_block_io, &
                           status,ierrmpi)
 end do
 ! write data communicated from other processors
 if(npe>1)then
  nrecv=(Morton_stop(npe-1)-Morton_start(1)+1)
  inrecv=0
  allocate(igrid_recv(nrecv))
  allocate(igrecvstatus(MPI_STATUS_SIZE,nrecv),iorecvstatus(MPI_STATUS_SIZE,nrecv))
  allocate(ioastatus(MPI_STATUS_SIZE,nrecv))

  do ipe =1, npe-1
   do Morton_no=Morton_start(ipe),Morton_stop(ipe)
     iwrite=iwrite+1
     itag=Morton_no
     inrecv=inrecv+1
     call MPI_RECV(igrid_recv(inrecv),1,MPI_INTEGER, ipe,itag,icomm,&
                   igrecvstatus(:,inrecv),ierrmpi)

     allocate(pwio(igrid_recv(inrecv))%w(ixG^T,1:nw))
     call MPI_RECV(pwio(igrid_recv(inrecv))%w,1,type_block_io,ipe,itag,icomm,&
                   iorecvstatus(:,inrecv),ierrmpi)
{#IFDEF TRANSFORMW
     if(nwtf>0 .and. neqpartf>0) then
       call transformw_usr(pwio(igrid_recv(inrecv))%w,wtf,eqpar_tf,ixG^LL,ixM^LL)
       offset=int(size_block_io_tf,kind=MPI_OFFSET_KIND) &
              *int(Morton_no-1,kind=MPI_OFFSET_KIND)
       call MPI_FILE_WRITE_AT(file_handle_tf,offset,wtf,1,type_block_io_tf,&
                               ioastatus(:,inrecv),ierrmpi)
     endif
}
     offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
          *int(Morton_no-1,kind=MPI_OFFSET_KIND)
     call MPI_FILE_WRITE_AT(file_handle,offset,pwio(igrid_recv(inrecv))%w,1,&
                          type_block_io,ioastatus(:,inrecv),ierrmpi)
     deallocate(pwio(igrid_recv(inrecv))%w)
   end do
  end do
  deallocate(igrecvstatus,iorecvstatus,ioastatus,igrid_recv)
 end if
end if

if(mype==0) call MPI_FILE_CLOSE(file_handle,ierrmpi)

{#IFDEF TRANSFORMW
if(nwtf>0 .and. neqpartf>0) then
  if (nwrite>0) then
    if(mype==0)deallocate(wtf)
  endif
  call MPI_FILE_CLOSE(file_handle_tf,ierrmpi)
endif
}

if (mype==0) then
   amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
   call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL, &
                      file_handle,ierrmpi)

   call write_forest(file_handle)

   {nx^D=ixMhi^D-ixMlo^D+1
   call MPI_FILE_WRITE(file_handle,nx^D,1,MPI_INTEGER,status,ierrmpi)\}
   call MPI_FILE_WRITE(file_handle,eqpar,neqpar+nspecialpar, &
                       MPI_DOUBLE_PRECISION,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nleafs,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,levmax,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndim,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,ndir,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,nw,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,neqpar+nspecialpar,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,it,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)

   call MPI_FILE_CLOSE(file_handle,ierrmpi)
end if
{#IFDEF TRANSFORMW
if (mype==0 .and. nwtf>0 .and. neqpartf>0) then
   amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
   call MPI_FILE_OPEN(MPI_COMM_SELF,filenametf,amode,MPI_INFO_NULL, &
                      file_handle_tf,ierrmpi)

   call write_forest(file_handle_tf)

   {nx^D=ixMhi^D-ixMlo^D+1
   call MPI_FILE_WRITE(file_handle_tf,nx^D,1,MPI_INTEGER,status,ierrmpi)\}
   call MPI_FILE_WRITE(file_handle_tf,eqpar_tf,neqpartf, &
                       MPI_DOUBLE_PRECISION,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,nleafs,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,levmax,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,ndim,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,ndir,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,nwtf,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,neqpartf,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,it,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(file_handle_tf,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)

   call MPI_FILE_CLOSE(file_handle_tf,ierrmpi)
endif
}
snapshot=snapshot+1

call MPI_BARRIER(icomm,ierrmpi)

end subroutine write_snapshot_nopar
!=============================================================================
subroutine write_snapshot_noparf
use mod_forest
include 'amrvacdef.f'

integer :: igrid, Morton_no
integer :: nx^D

integer, allocatable :: iostatus(:,:),iorecvstatus(:,:),ioastatus(:,:)
integer, allocatable :: igrecvstatus(:,:)
integer, allocatable :: iorequest(:),igrid_recv(:) 

integer  :: ipe,insend,inrecv,nrecv,nwrite
character(len=80) :: filename, line
logical, save :: firstsnapshot=.true.
!-----------------------------------------------------------------------------
call MPI_BARRIER(icomm,ierrmpi)

if (firstsnapshot) then
   snapshot=snapshotnext
   firstsnapshot=.false.
end if

if (snapshot >= 10000) then
   if (mype==0) then
      write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
      write(*,*) "overwriting first frames"
   end if
   snapshot=0
end if

nrecv=0
inrecv=0
nwrite=0
insend=0

if (mype /= 0) then
 do Morton_no=Morton_start(mype),Morton_stop(mype)
    igrid=sfc_to_igrid(Morton_no)
    itag=Morton_no
    insend=insend+1
    if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine, 
      ! which might be used in getaux
      saveigrid=igrid
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      mygeo =>pgeo(igrid)
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
    endif
    call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
    call MPI_SEND(pw(igrid)%w,1,type_block_io, 0,itag,icomm,ierrmpi)
 end do
else 
 ! mype==0
 nwrite=(Morton_stop(0)-Morton_start(0)+1)
 allocate(iorequest(nwrite),iostatus(MPI_STATUS_SIZE,nwrite))
 iorequest=MPI_REQUEST_NULL

 ! master processor writes out
 write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"
 if(endian_swap) then
  {#IFNDEF BIGENDIAN
   open(unit=unitsnapshot,file=filename,form='unformatted',access='stream',&
        status='replace',convert='BIG_ENDIAN')
  }
  {#IFDEF BIGENDIAN
   open(unit=unitsnapshot,file=filename,form='unformatted',access='stream',&
        status='replace',convert='LITTLE_ENDIAN')
  }
 else
   open(unit=unitsnapshot,file=filename,form='unformatted',access='stream',&
        status='replace')
 end if
 ! writing his local data first
 do Morton_no=Morton_start(0),Morton_stop(0)
   igrid=sfc_to_igrid(Morton_no)
   if (nwaux>0) then
      ! extra layer around mesh only for later averaging in convert
      ! set dxlevel value for use in gradient subroutine,
      ! which might be used in getaux
      saveigrid=igrid
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      mygeo =>pgeo(igrid)
      if (B0field) then
         myB0_cell => pB0_cell(igrid)
         {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      call getaux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
   endif
   write(unitsnapshot) pw(igrid)%w(ixM^T,1:nw)
 end do
 ! write data communicated from other processors
 if(npe>1)then
  nrecv=(Morton_stop(npe-1)-Morton_start(1)+1)
  inrecv=0
  allocate(igrid_recv(nrecv))
  allocate(igrecvstatus(MPI_STATUS_SIZE,nrecv),iorecvstatus(MPI_STATUS_SIZE,nrecv))
  allocate(ioastatus(MPI_STATUS_SIZE,nrecv))

  do ipe =1, npe-1
   do Morton_no=Morton_start(ipe),Morton_stop(ipe)
     itag=Morton_no
     inrecv=inrecv+1
     call MPI_RECV(igrid_recv(inrecv),1,MPI_INTEGER, ipe,itag,icomm,&
                   igrecvstatus(:,inrecv),ierrmpi)

     allocate(pwio(igrid_recv(inrecv))%w(ixG^T,1:nw))
     call MPI_RECV(pwio(igrid_recv(inrecv))%w,1,type_block_io,ipe,itag,icomm,&
                   iorecvstatus(:,inrecv),ierrmpi)

     write(unitsnapshot) pwio(igrid_recv(inrecv))%w(ixM^T,1:nw)
     deallocate(pwio(igrid_recv(inrecv))%w)
   end do
  end do
  deallocate(igrecvstatus,iorecvstatus,ioastatus,igrid_recv)
 end if
end if

if(nwrite>0)then
  call MPI_WAITALL(nwrite,iorequest,iostatus,ierrmpi) 
  if(mype==0)deallocate(iorequest,iostatus)
end if

if(mype==0) then
  call write_forest(unitsnapshot)
  {nx^D=ixMhi^D-ixMlo^D+1
  write(unitsnapshot) nx^D\}
  write(unitsnapshot) eqpar
  write(unitsnapshot) nleafs
  write(unitsnapshot) levmax 
  write(unitsnapshot) ndim 
  write(unitsnapshot) ndir
  write(unitsnapshot) nw
  write(unitsnapshot) neqpar+nspecialpar
  write(unitsnapshot) it
  write(unitsnapshot) t
  close(unitsnapshot)
end if

snapshot=snapshot+1

call MPI_BARRIER(icomm,ierrmpi)

end subroutine write_snapshot_noparf
!=============================================================================
subroutine read_snapshot
use mod_forest
include 'amrvacdef.f'

integer :: file_handle, amode, igrid, Morton_no, iread
integer :: levmaxini, ndimini, ndirini, nwini, neqparini, nxini^D

{^IFMPT integer :: size_double, size_int,lb}
{^IFNOMPT  integer(kind=MPI_ADDRESS_KIND) :: size_double, size_int, lb}

integer(kind=MPI_OFFSET_KIND) :: offset
integer, dimension(MPI_STATUS_SIZE,ngridshi) :: iostatus
integer, dimension(MPI_STATUS_SIZE) :: status
character(len=80) :: filename
logical :: fexist
!-----------------------------------------------------------------------------
!!!call MPI_BARRIER(icomm,ierrmpi)
! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameini),snapshotini,".dat"

if(mype==0) then
  inquire(file=filename,exist=fexist)
  if(.not.fexist) call mpistop(filename//"as an input snapshot file is not found!")
endif

amode=MPI_MODE_RDONLY
call MPI_FILE_OPEN(icomm,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)


!call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,size_double,ierrmpi)
!call MPI_TYPE_EXTENT(MPI_INTEGER,size_int,ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,size_int,ierrmpi)

{#IFDEF EVOLVINGBOUNDARY
offset=-int(8*size_int+size_double,kind=MPI_OFFSET_KIND)
}{#IFNDEF EVOLVINGBOUNDARY
offset=-int(7*size_int+size_double,kind=MPI_OFFSET_KIND)
}
!call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

call MPI_FILE_READ_ALL(file_handle,nleafs,1,MPI_INTEGER,status,ierrmpi)
nleafs_active = nleafs
call MPI_FILE_READ_ALL(file_handle,levmaxini,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,ndimini,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,ndirini,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,nwini,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,neqparini,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,it,1,MPI_INTEGER,status,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)
{#IFDEF EVOLVINGBOUNDARY
call MPI_FILE_READ_ALL(file_handle,nphyboundblock,1,MPI_INTEGER,status,ierrmpi)
}

! check if settings are suitable for restart
if (levmaxini>mxnest) then
   if (mype==0) write(*,*) "number of levels in restart file = ",levmaxini
   if (mype==0) write(*,*) "mxnest = ",mxnest
   call mpistop("mxnest should be at least number of levels in restart file")
end if
if (ndimini/=ndim) then
   if (mype==0) write(*,*) "ndim in restart file = ",ndimini
   if (mype==0) write(*,*) "ndim = ",ndim
   call mpistop("reset ndim to ndim in restart file")
end if
if (ndirini/=ndir) then
   if (mype==0) write(*,*) "ndir in restart file = ",ndirini
   if (mype==0) write(*,*) "ndir = ",ndir
   call mpistop("reset ndir to ndir in restart file")
end if
if (nw/=nwini) then
   if (mype==0) write(*,*) "nw=",nw," and nw in restart file=",nwini
   call mpistop("currently, changing nw at restart is not allowed")
end if

offset=offset-int(ndimini*size_int+neqparini*size_double,kind=MPI_OFFSET_KIND)
!call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

{call MPI_FILE_READ_ALL(file_handle,nxini^D,1,MPI_INTEGER,status,ierrmpi)\}
if (ixGhi^D/=nxini^D+2*dixB|.or.) then
   if (mype==0) write(*,*) "Error: reset resolution to ",nxini^D+2*dixB
   call mpistop("change with setup.pl")
end if
neqparini=min(neqparini,neqpar+nspecialpar)
call MPI_FILE_READ_ALL(file_handle,eqpar,neqparini, &
                       MPI_DOUBLE_PRECISION,status,ierrmpi)



call read_forest(file_handle)

iread=0
{#IFDEF EVOLVINGBOUNDARY
! mark physical-boundary blocks on space-filling curve
do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   call alloc_node(igrid)
   if (phyboundblock(igrid)) sfc_phybound(Morton_no)=1
end do
call MPI_ALLREDUCE(MPI_IN_PLACE,sfc_phybound,nleafs,MPI_INTEGER,&
                   MPI_SUM,icomm,ierrmpi)

do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   iread=iread+1
   nphyboundblock=sum(sfc_phybound(1:Morton_no-1))
   offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
          *int(Morton_no-1-nphyboundblock,kind=MPI_OFFSET_KIND) + &
          int(size_block,kind=MPI_OFFSET_KIND) &
          *int(nphyboundblock,kind=MPI_OFFSET_KIND)
   if (sfc_phybound(Morton_no)==1) then
      call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1,type_block, &
                             status,ierrmpi)
   else
      call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1,type_block_io, &
                             status,ierrmpi)
   end if
end do
}{#IFNDEF EVOLVINGBOUNDARY
do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   call alloc_node(igrid)
   iread=iread+1
   offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
          *int(Morton_no-1,kind=MPI_OFFSET_KIND)
   call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1,type_block_io, &
                          status,ierrmpi)
end do
}

call MPI_FILE_CLOSE(file_handle,ierrmpi)

!!!call MPI_BARRIER(icomm,ierrmpi)
end subroutine read_snapshot
!=============================================================================
subroutine read_snapshotnopar
use mod_forest
include 'amrvacdef.f'

double precision :: wio(ixG^T,1:nw)
integer :: file_handle, amode, igrid, Morton_no, iread
integer :: levmaxini, ndimini, ndirini, nwini, neqparini, nxini^D

{^IFMPT integer :: size_double, size_int,lb}
{^IFNOMPT  integer(kind=MPI_ADDRESS_KIND) :: size_double, size_int, lb}

integer(kind=MPI_OFFSET_KIND) :: offset
integer, dimension(MPI_STATUS_SIZE) :: status
integer, dimension(MPI_STATUS_SIZE) :: iostatus

integer, allocatable :: iorecvstatus(:,:)
integer :: ipe,inrecv,nrecv
integer :: sendini(7+^ND)
character(len=80) :: filename
logical :: fexist
!-----------------------------------------------------------------------------
!!!call MPI_BARRIER(icomm,ierrmpi)
! generate filename
write(filename,"(a,i4.4,a)") TRIM(filenameini),snapshotini,".dat"
if (mype==0) then
 inquire(file=filename,exist=fexist)
 if(.not.fexist) call mpistop(filename//"as an input snapshot file is not found!")
 amode=MPI_MODE_RDONLY
 call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)

 !call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,size_double,ierrmpi)
 !call MPI_TYPE_EXTENT(MPI_INTEGER,size_int,ierrmpi)
 call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
 call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,size_int,ierrmpi)

 offset=-int(7*size_int+size_double,kind=MPI_OFFSET_KIND)
 !call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
 call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

 call MPI_FILE_READ(file_handle,nleafs,1,MPI_INTEGER,status,ierrmpi)
 nleafs_active = nleafs
 call MPI_FILE_READ(file_handle,levmaxini,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,ndimini,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,ndirini,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,nwini,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,neqparini,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,it,1,MPI_INTEGER,status,ierrmpi)
 call MPI_FILE_READ(file_handle,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)
 ! check if settings are suitable for restart
 if (levmaxini>mxnest) then
    write(*,*) "number of levels in restart file = ",levmaxini
    write(*,*) "mxnest = ",mxnest
    call mpistop("mxnest should be at least number of levels in restart file")
 end if
 if (ndimini/=ndim) then
    write(*,*) "ndim in restart file = ",ndimini
    write(*,*) "ndim = ",ndim
    call mpistop("reset ndim to ndim in restart file")
 end if
 if (ndirini/=ndir) then
    write(*,*) "ndir in restart file = ",ndirini
    write(*,*) "ndir = ",ndir
    call mpistop("reset ndir to ndir in restart file")
 end if
 if (nw/=nwini) then
    write(*,*) "nw=",nw," and nw in restart file=",nwini
    call mpistop("currently, changing nw at restart is not allowed")
 end if


 offset=offset-int(ndimini*size_int+neqparini*size_double,kind=MPI_OFFSET_KIND)
 !call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
 call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

 {call MPI_FILE_READ(file_handle,nxini^D,1,MPI_INTEGER,status,ierrmpi)\}
 if (ixGhi^D/=nxini^D+2*dixB|.or.) then
    write(*,*) "Error: reset resolution to ",nxini^D+2*dixB
    call mpistop("change with setamrvac")
 end if
 neqparini=min(neqparini,neqpar+nspecialpar)
 call MPI_FILE_READ(file_handle,eqpar,neqparini, &
                       MPI_DOUBLE_PRECISION,status,ierrmpi)
end if

! broadcast the global parameters first
if (npe>1) then
  if (mype==0) then
     sendini=(/nleafs,levmaxini,ndimini,ndirini,nwini,neqparini,it ,^D&nxini^D /)
  end if
  call MPI_BCAST(sendini,7+^ND,MPI_INTEGER,0,icomm,ierrmpi)
  nleafs=sendini(1);levmaxini=sendini(2);ndimini=sendini(3);
  ndirini=sendini(4);nwini=sendini(5);
  neqparini=sendini(6);it=sendini(7);
  nxini^D=sendini(7+^D);
  nleafs_active = nleafs
  call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
  call MPI_BCAST(eqpar,neqparini,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
end if

call read_forest(file_handle)

if (mype==0)then
   iread=0
   do Morton_no=Morton_start(0),Morton_stop(0)
      igrid=sfc_to_igrid(Morton_no)
      call alloc_node(igrid)
      iread=iread+1
      offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
             *int(Morton_no-1,kind=MPI_OFFSET_KIND)
      call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1,type_block_io, &
                            iostatus,ierrmpi)
   end do
   if (npe>1) then
    do ipe=1,npe-1
     do Morton_no=Morton_start(ipe),Morton_stop(ipe)
       iread=iread+1
       itag=Morton_no
       offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
              *int(Morton_no-1,kind=MPI_OFFSET_KIND)
       call MPI_FILE_READ_AT(file_handle,offset,wio,1,type_block_io, &
                            iostatus,ierrmpi)
       call MPI_SEND(wio,1,type_block_io, ipe,itag,icomm,ierrmpi)
     end do
    end do
   end if
   call MPI_FILE_CLOSE(file_handle,ierrmpi)
else
   nrecv=(Morton_stop(mype)-Morton_start(mype)+1)
   allocate(iorecvstatus(MPI_STATUS_SIZE,nrecv))
   inrecv=0
   do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid=sfc_to_igrid(Morton_no)
      itag=Morton_no
      call alloc_node(igrid)
      inrecv=inrecv+1
      call MPI_RECV(pw(igrid)%w,1,type_block_io,0,itag,icomm,&
                    iorecvstatus(:,inrecv),ierrmpi)
   end do
   deallocate(iorecvstatus)
end if

call MPI_BARRIER(icomm,ierrmpi)

end subroutine read_snapshotnopar
!=============================================================================
subroutine printlog_default

! printlog: calculates volume averaged mean values 
use mod_timing
use mod_forest,only:nleafs,nleafs_active,nleafs_level
include 'amrvacdef.f'

logical          :: fileopen
integer          :: iigrid, igrid, level, iw, i
double precision :: wmean(1:nw), volume(1:nlevelshi), volprob, voltotal
double precision :: dvolume(ixG^T), volumeflat(1:nlevelshi)
integer          :: numlevels, nx^D, nc, ncells, dit
double precision :: dtTimeLast, now, cellupdatesPerSecond, activeBlocksPerCore, wctPerCodeTime, timeToFinish
integer, dimension(1:nlevelshi) :: isum_send, isum_recv
double precision, dimension(1:nw+1+nlevelshi) :: dsum_send, dsum_recv
character(len=80) :: filename
character(len=2048) :: line
logical, save :: opened=.false.
integer :: amode, status(MPI_STATUS_SIZE)
!-----------------------------------------------------------------------------

volume(1:mxnest)=zero
volumeflat(1:mxnest)=zero
wmean(1:nw)= zero

do iigrid=1,igridstail; igrid=igrids(iigrid);
   level=node(plevel_,igrid)
   volumeflat(level)=volumeflat(level)+ &
          {(rnode(rpxmax^D_,igrid)-rnode(rpxmin^D_,igrid))|*}
   if (slab) then
      dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*}
   else
      dvolume(ixM^T)=pgeo(igrid)%dvolume(ixM^T)
      volume(level)=volume(level)+sum(dvolume(ixM^T))
   end if
   do iw=1,nw
      wmean(iw)=wmean(iw)+sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,iw))
   end do
end do
if (slab) volume(levmin:levmax)=volumeflat(levmin:levmax)

voltotal=sum(volume(levmin:levmax))

numlevels=levmax-levmin+1
dsum_send(1:nw)=wmean(1:nw)
dsum_send(nw+1)=voltotal
dsum_send(nw+2:nw+1+numlevels)=volumeflat(levmin:levmax)
call MPI_REDUCE(dsum_send,dsum_recv,nw+1+numlevels,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,icomm,ierrmpi)

if (mype==0) then

! To compute cell updates per second, we do the following:
nx^D=ixMhi^D-ixMlo^D+1;
nc={nx^D*}
ncells = nc * nleafs_active
! assumes the number of active leafs haven't changed since last compute.
now        = MPI_WTIME()
dit        = it - itTimeLast
dtTimeLast = now - timeLast
itTimeLast = it
timeLast   = now
cellupdatesPerSecond = dble(ncells) * dble(nstep) * dble(dit) / (dtTimeLast * dble(npe))
! blocks per core:
activeBlocksPerCore = dble(nleafs_active) / dble(npe)

! Wall clock time per code time unit in seconds:
wctPerCodeTime = dtTimeLast / max(dit * dt, epsilon(1.0d0))

! Wall clock time to finish in hours:
timeToFinish = (tmax - t) * wctPerCodeTime / 3600.0d0

   wmean(1:nw)=dsum_recv(1:nw)
   voltotal=dsum_recv(nw+1)
   volumeflat(levmin:levmax)=dsum_recv(nw+2:nw+1+numlevels)

   wmean=wmean/voltotal

   ! determine coverage in coordinate space
   volprob={(xprobmax^D-xprobmin^D)|*}
   volumeflat(levmin:levmax)=volumeflat(levmin:levmax)/volprob

   if (.not.opened) then
      ! generate filename
      write(filename,"(a,a)") TRIM(filenamelog),".log"

      amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      amode=ior(amode,MPI_MODE_APPEND)
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode, &
                         MPI_INFO_NULL,log_fh,ierrmpi)
      opened=.true.

      call MPI_FILE_WRITE(log_fh,fileheadout,len_trim(fileheadout), &
                          MPI_CHARACTER,status,ierrmpi)
      !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
      call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

      i=len_trim(wnameslog)-1
      do level=1,mxnest
          i=i+3
          if(level<10) then
            if (i+1<1024) write(wnameslog(i:i+1),"(a,i1)") "c",level
          else
            if (i+2<1024) write(wnameslog(i:i+2),"(a,i2)") "c",level
          endif
      end do

      do level=1,mxnest
          i=i+3
          if(level<10) then
            if (i+1<1024) write(wnameslog(i:i+1),"(a,i1)") "n",level
          else
            if (i+2<1024) write(wnameslog(i:i+2),"(a,i2)") "n",level
          endif
      end do
      if (time_accurate) then
         if(residmin>smalldouble) then
           write(line,'(a15,a1024)')"it   t  dt res ",wnameslog
         else
           write(line,'(a15,a1024)')"it   t   dt    ",wnameslog
         endif
      else
         if(residmin>smalldouble) then
           write(line,'(a7,a1024)')"it res ",wnameslog
         else
           write(line,'(a7,a1024)')"it     ",wnameslog
         endif
      end if

      line=trim(line)//"| Xload Xmemory 'Cell_Updates /second/core'"
      line=trim(line)//" 'Active_Blocks/Core' 'Wct Per Code Time [s]' 'TimeToFinish [hrs]'"

      
      call MPI_FILE_WRITE(log_fh,line,len_trim(line),MPI_CHARACTER, &
                          status,ierrmpi)
   end if
   !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
   call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

   if (time_accurate) then
      if(residmin>smalldouble) then
         write(line,'(i7,3(es12.4))')it,t,dt,residual
      else
         write(line,'(i7,2(es12.4))')it,t,dt
      endif
   else
      if(residmin>smalldouble) then
         write(line,'(i7,1(es12.4))')it,residual
      else
         write(line,'(i7)')it
      endif
   end if
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                       MPI_CHARACTER,status,ierrmpi)
   do iw=1,nw
      write(line,'(es12.4)')wmean(iw)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do
   do level=1,mxnest
      write(line,'(es12.4)')volumeflat(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do

   do level=1,mxnest
      write(line,'(i8)') nleafs_level(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do

   write(line,'(a3,6(es10.2))') ' | ', Xload, Xmemory, cellupdatesPerSecond, &
        activeBlocksPerCore, wctPerCodeTime, timeToFinish
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)

end if

end subroutine printlog_default

!> Compute mean(w**power) over the leaves of the grid. The first mode
!> (power=1) corresponds to to the mean, the second to the mean squared values
!> and so on.
subroutine get_volume_average(power, mode, volume)
  include 'amrvacdef.f'

  integer, intent(in)           :: power     !< Which mode to compute
  double precision, intent(out) :: mode(nw)  !< The computed mode
  double precision, intent(out) :: volume    !< The total grid volume
  integer                       :: iigrid, igrid, iw
  double precision              :: wsum(nw+1)
  double precision              :: dvolume(ixG^T)
  double precision              :: dsum_recv(1:nw+1)

  wsum(:) = 0

  ! Loop over all the grids
  do iigrid = 1, igridstail
     igrid = igrids(iigrid)

     ! Determine the volume of the grid cells
     if (slab) then
        dvolume(ixM^T) = {rnode(rpdx^D_,igrid)|*}
     else
        dvolume(ixM^T) = pgeo(igrid)%dvolume(ixM^T)
     end if

     ! Store total volume in last element
     wsum(nw+1) = wsum(nw+1) + sum(dvolume(ixM^T))

     ! Compute the modes of the cell-centered variables, weighted by volume
     do iw = 1, nw
        wsum(iw) = wsum(iw) + &
             sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,iw)**power)
     end do
  end do

  ! Make the information available on all tasks
  call MPI_ALLREDUCE(wsum, dsum_recv, nw+1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, icomm, ierrmpi)

  ! Set the volume and the average
  volume = dsum_recv(nw+1)
  mode   = dsum_recv(1:nw) / volume

end subroutine get_volume_average

!> Compute the volume average of func(w) over the leaves of the grid.
subroutine get_volume_average_func(func, f_avg, volume)
  include 'amrvacdef.f'

  interface
     pure function func(w_vec, w_size) result(val)
       integer, intent(in)          :: w_size
       double precision, intent(in) :: w_vec(w_size)
       double precision             :: val
     end function func
  end interface
  double precision, intent(out) :: f_avg  !< The volume average of func
  double precision, intent(out) :: volume    !< The total grid volume
  integer                       :: iigrid, igrid, i^D
  double precision              :: wsum(2)
  double precision              :: dvolume(ixG^T)
  double precision              :: dsum_recv(2)

  wsum(:) = 0

  ! Loop over all the grids
  do iigrid = 1, igridstail
     igrid = igrids(iigrid)

     ! Determine the volume of the grid cells
     if (slab) then
        dvolume(ixM^T) = {rnode(rpdx^D_,igrid)|*}
     else
        dvolume(ixM^T) = pgeo(igrid)%dvolume(ixM^T)
     end if

     ! Store total volume in last element
     wsum(2) = wsum(2) + sum(dvolume(ixM^T))

     ! Compute the modes of the cell-centered variables, weighted by volume
     {do i^D = ixMlo^D, ixMhi^D\}
     wsum(1) = wsum(1) + dvolume(i^D) * &
          func(pw(igrid)%w(i^D, :), nw)
     {end do\}
  end do

  ! Make the information available on all tasks
  call MPI_ALLREDUCE(wsum, dsum_recv, 2, MPI_DOUBLE_PRECISION, &
       MPI_SUM, icomm, ierrmpi)

  ! Set the volume and the average
  volume = dsum_recv(2)
  f_avg  = dsum_recv(1) / volume

end subroutine get_volume_average_func
!=============================================================================
