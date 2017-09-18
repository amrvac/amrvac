!##############################################################################
! include amrvacusrpar - binaryhdWR98a

INCLUDE:amrvacmodules/coolingpar.t


INTEGER,PARAMETER:: MassA_=neqpar+1,MassB_=MassA_+1,jumpsh_=MassB_+1, &
                    vinfA_=jumpsh_+1,vinfB_=vinfA_+1,rinfA_=vinfB_+1,rinfB_=rinfA_+1, &
                    rhoinfA_=rinfB_+1,rhoinfB_=rhoinfA_+1,Porb_=rhoinfB_+1,ecc_=Porb_+1, &
                    semimajor_=ecc_+1,rhoism_=semimajor_+1,pism_=rhoism_+1,&
                    TwindA_=pism_+1,TwindB_=TwindA_+1,Tscale_=TwindB_+1,Lscale_=Tscale_+1, &
                    Mue_=Lscale_+1,DfracA_=Mue_+1,DfracB_=DfracA_+1, &
                    min_ar_=DfracB_+1,max_ar_=min_ar_+1,tfix_=max_ar_+1,&
                    Macc_=tfix_+1,nspecialpar=25
CHARACTER*74,PARAMETER:: specialparname= &
'mA mB js vA vB rA rB hA hB P e a rhs pi TA TB Ts Ls Mu DA DB ma mi tf Macc'

! constants in cgs units here
double precision, parameter :: kboltz = 1.38065D-16
double precision, parameter :: mhydro = 1.6733D-24
double precision, parameter :: msolar = 1.989D+33
double precision, parameter :: rsolar = 7.0D+10
double precision, parameter :: sinyear= 3.1536D+7
double precision, parameter :: ggrav  = 6.671D-8

INTEGER, PARAMETER:: jmax=5000

integer :: nxmapA,nxmapB
double precision :: rmapA,logrhomapA,vmapA
double precision :: rmapB,logrhomapB,vmapB
double precision :: rfracA1,rfracA2,rfracB1,rfracB2

common /fracmap  / rfracA1,rfracA2,rfracB1,rfracB2
common /iwindmap / nxmapA, nxmapB
common /windmaps / rmapA(jmax),rmapB(jmax), &
              logrhomapA(jmax),logrhomapB(jmax),vmapA(jmax),vmapB(jmax)

! end include amrvacusrpar - binaryhdWR98a
!##############################################################################
