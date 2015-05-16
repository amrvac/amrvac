!==============================================================================
! include file amrvacpar.t.srhdeos
!------------------------------------------------------------------------------

CHARACTER*7,PARAMETER:: typephys='srhdeos'           ! VACPHYS module name
CHARACTER*5,PARAMETER:: eqparname='gamma'        ! Equation parameter names

integer,parameter:: b0_=-1 ! No magnetic field
integer,parameter:: b^C_=-1 ! No magnetic field
integer,parameter:: e_=-1 ! for compilation of convert routine

integer,parameter:: d_=1
integer,parameter:: s0_=d_
integer,parameter:: s^C_=s0_+^C
integer,parameter:: tau_=s^NC_+1
{#IFDEF TRACER
! Scalar tracers
integer,parameter:: Dtr^FL_=tau_+^FL
{#IFDEF EPSINF
!======Cutoff energy ========!
integer,parameter:: Depsinf_=Dtr^NFL_+1
!======Conserved electron density=====!
integer,parameter::Dne_  = Depsinf_+1
integer,parameter::Drho1_  = Dne_             ! for compatibility with the srmhd nomenclature (used in cooling)
!======Advected electron density======!
integer,parameter::Dne0_ = Depsinf_+2
! Number of variables
integer,parameter:: nwflux=Dne0_
}{#IFNDEF EPSINF
! Number of variables
integer,parameter:: nwflux=Dtr^NFL_
}}
{#IFNDEF TRACER
{#IFDEF EPSINF
!======Cutoff energy ========!
integer,parameter:: Depsinf_=tau_+1
!======Conserved electron density=====!
integer,parameter::Dne_  = Depsinf_+1
integer,parameter::Drho1_  = Dne_             ! for compatibility with the srmhd nomenclature (used in cooling)
!======Advected electron density======!
integer,parameter::Dne0_ = Depsinf_+2
! Number of variables
integer,parameter:: nwflux=Dne0_
}{#IFNDEF EPSINF
! Number of variables
integer,parameter:: nwflux=tau_
}}
! auxiliary variables
integer,parameter:: lfac_=nwflux+1
integer,parameter:: p_=lfac_+1

integer,parameter:: nwaux=2
integer,parameter:: nwextra=0
integer,parameter:: nw=nwflux+nwaux+nwextra


!Primitive variables
integer,parameter:: rho_=d_
integer,parameter:: v0_=s0_
integer,parameter:: v^C_=v0_+^C
integer,parameter:: u0_=s0_
integer,parameter:: u^C_=u0_+^C
integer,parameter:: pp_=tau_
{#IFDEF TRACER
integer,parameter:: tr^FL_=Dtr^FL_
}
{#IFDEF EPSINF
integer,parameter:: epsinf_=Depsinf_
integer,parameter:: ne_=Dne_
integer,parameter:: rho1_  = Dne_             ! for compatibility with the srmhd nomenclature (used in cooling)
integer,parameter:: ne0_=Dne0_
}
! polar variable names
integer,parameter:: sr_=s0_+r_
integer,parameter:: sphi_=s0_+phi_
integer,parameter:: sz_=s0_+z_
integer,parameter:: vr_=v0_+r_
integer,parameter:: vphi_=v0_+phi_
integer,parameter:: vz_=v0_+z_

integer, parameter :: nvector=1                             ! No. vector vars
integer, dimension(nvector), parameter :: iw_vector=(/ s0_ /)

                                                            ! Characteristic
INTEGER,PARAMETER:: soundRW_=1,soundLW_=2,entropW_=3,shearW0_=3      ! waves
INTEGER,PARAMETER:: nworkroe=3      

{#IFNDEF EPSINF
INTEGER,PARAMETER:: gamma_=1,neqpar=1                     ! equation parameters
}{#IFDEF EPSINF
INTEGER,PARAMETER:: gamma_=1, epsfloor_=2, rho0floor_=3, rho1e_=4, neqpar=4     ! equation params
}
COMMON, DOUBLE PRECISION::smalltau,smallxi,minrho,minp

INTEGER,PARAMETER:: nflag_=nw+1
COMMON, INTEGER:: flags(nflag_)
COMMON, DOUBLE PRECISION:: wflags(nflag_)

! xprob: problem box; iprob: problem
COMMON, INTEGER:: iprob
COMMON, DOUBLE PRECISION:: xprob^L

! end include file amrvacpar.t.srhdeos
!==============================================================================
