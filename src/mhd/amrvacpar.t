!===============================================================================
! include file amrvacpar.t.mhd

CHARACTER*3,PARAMETER:: typephys='mhd'            ! VACPHYS module name

{#IFDEF GLM 
{#IFDEF ENERGY
CHARACTER*11,PARAMETER:: eqparname='gamma eta etah etahyper Cr'    ! Equation parameter names
}{#IFDEF ISO
CHARACTER*11,PARAMETER:: eqparname='gamma eta etah etahyper Cr adiab'    ! Equation parameter names
}}
{#IFNDEF GLM 
{#IFDEF ENERGY
CHARACTER*11,PARAMETER:: eqparname='gamma eta etah etahyper'    ! Equation parameter names
}{#IFDEF ISO
CHARACTER*11,PARAMETER:: eqparname='gamma eta etah etahyper adiab'    ! Equation parameter names
}}

INTEGER,PARAMETER:: rho_=1,m0_=rho_,m^C_=m0_+^C
{#IFDEF GLM 
{#IFDEF ENERGY
INTEGER,PARAMETER:: e_=m^NC_+1
INTEGER,PARAMETER:: b0_=e_,b^C_=b0_+^C,psi_=b^NC_+1  ! flow variables
INTEGER,PARAMETER:: nwflux=3+2*^NC+^NFL
}{#IFNDEF ENERGY
INTEGER,PARAMETER:: b0_=m^NC_,b^C_=b0_+^C,psi_=b^NC_+1  ! flow variables
INTEGER,PARAMETER:: nwflux=2+2*^NC+^NFL
}
{#IFDEF TRACER
INTEGER,PARAMETER:: Dtr^FL_=psi_+^FL
}}
{#IFNDEF GLM
{#IFDEF ENERGY
INTEGER,PARAMETER:: e_=m^NC_+1
INTEGER,PARAMETER:: b0_=e_,b^C_=b0_+^C  ! flow variables
INTEGER,PARAMETER:: nwflux=2+2*^NC+^NFL
}{#IFNDEF ENERGY
INTEGER,PARAMETER:: b0_=m^NC_,b^C_=b0_+^C  ! flow variables
INTEGER,PARAMETER:: nwflux=1+2*^NC+^NFL
}
{#IFDEF TRACER
INTEGER,PARAMETER:: Dtr^FL_=b^NC_+^FL
}}

integer,parameter:: nwaux=0
integer,parameter:: nwextra=0
integer,parameter:: nw=nwflux+nwaux+nwextra

{#IFDEF ENERGY
integer,parameter:: ee_=e_
integer,parameter:: rhos_=e_
INTEGER, PARAMETER:: p_=e_, pp_=ee_    ! Primitive variables
}{#IFDEF TRACER
integer,parameter:: tr^FL_=Dtr^FL_
}
INTEGER,PARAMETER:: v0_=m0_, v^C_=m^C_

INTEGER,PARAMETER:: mr_=m0_+r_,mphi_=m0_+phi_,mz_=m0_+z_  ! Polar var. names
INTEGER,PARAMETER:: br_=b0_+r_,bphi_=b0_+phi_,bz_=b0_+z_

integer, parameter :: nvector=2                             ! No. vector vars
integer, dimension(nvector), parameter :: iw_vector=(/ m0_, b0_ /)

INTEGER,PARAMETER:: fastRW_=1,fastLW_=2,slowRW_=3,slowLW_=4 ! Characteristic
INTEGER,PARAMETER:: entroW_=5,diverW_=6,alfvRW_=7,alfvLW_=8 ! waves
INTEGER,PARAMETER:: nworkroe=15


{#IFDEF GLM 
{#IFDEF ENERGY
INTEGER,PARAMETER:: gamma_=1,eta_=2,etah_=3,etahyper_=4,Cr_=5,neqpar=5     ! equation params
}{#IFDEF ISO
INTEGER,PARAMETER:: gamma_=1,eta_=2,etah_=2,etahyper_=4,Cr_=5,adiab_=6,neqpar=6     ! equation params
}}
{#IFNDEF GLM 
{#IFDEF ENERGY
INTEGER,PARAMETER:: gamma_=1,eta_=2,etah_=3,etahyper_=4,neqpar=4     ! equation params
}{#IFDEF ISO
INTEGER,PARAMETER:: gamma_=1,eta_=2,etah_=3,etahyper_=4,adiab_=5,neqpar=5     ! equation params
}}

COMMON, DOUBLE PRECISION::smalle,minrho,minp
INTEGER,PARAMETER:: nflag_=nw+1
COMMON, INTEGER:: flags(nflag_)
COMMON, DOUBLE PRECISION:: wflags(nflag_)

! xprob: problem box; iprob: problem
COMMON, INTEGER:: iprob
COMMON, DOUBLE PRECISION:: xprob^L

! end include file amrvacpar.t.mhd
!===============================================================================
