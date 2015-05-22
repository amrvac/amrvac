!==============================================================================
! include file amrvacpar.t.hd

{#IFDEF ENERGY
{#IFNDEF DUST
CHARACTER*2,PARAMETER:: typephys='hd'            ! physics id for idl
}
{#IFDEF DUST
CHARACTER*7,PARAMETER:: typephys='hdmdust'            ! physics id for idl
}
CHARACTER*8,PARAMETER:: eqparname='gamma mu'        ! Equation parameter names
integer,parameter:: nwflux=^NC*(1+^NDS)+2+^NDS+^NFL
}

{#IFDEF ISO
{#IFNDEF DUST
CHARACTER*7,PARAMETER:: typephys='hdadiab'            ! physics id for idl
}
{#IFDEF DUST
CHARACTER*11,PARAMETER:: typephys='hdadiabdust'            ! physics id for idl
}
CHARACTER*14,PARAMETER:: eqparname='gamma adiab mu'        ! Equation parameter names
integer,parameter:: nwflux=^NC*(1+^NDS)+1+^NDS+^NFL
}

integer,parameter:: nwaux=0
integer,parameter:: nwextra=0
integer,parameter:: nw=nwflux+nwaux+nwextra

integer,parameter:: rho_=1
integer,parameter:: m0_=rho_
integer,parameter:: m^C_=m0_+^C

{#IFDEF ENERGY
integer,parameter:: e_=m^NC_+1
integer,parameter:: ee_=e_
integer,parameter:: rhos_=e_
{#IFDEF TRACER
INTEGER,PARAMETER:: Dtr^FL_=e_+^FL
{#IFDEF DUST
INTEGER,PARAMETER:: rhod0_=e_+^NFL
}}
{#IFNDEF TRACER
{#IFDEF DUST
INTEGER,PARAMETER:: rhod0_=e_
}}}

{#IFNDEF ENERGY
{#IFDEF TRACER
INTEGER,PARAMETER:: Dtr^FL_=m^NC_+^FL
{#IFDEF DUST
INTEGER,PARAMETER:: rhod0_=m^NC_+^NFL
}}
{#IFNDEF TRACER
{#IFDEF DUST
INTEGER,PARAMETER:: rhod0_=m^NC_
}}}

{#IFDEF DUST
integer,parameter:: rhod^DS_=rhod0_+^DS
integer,parameter:: m0d^DS_=rhod^DS_
integer,parameter:: m1d^DS_=rhod^NDS_+^DS
{^NOONEC
integer,parameter:: m2d^DS_=m1d^NDS_+^DS }
{^IFTHREEC
integer,parameter:: m3d^DS_=m2d^NDS_+^DS }
}

integer,parameter:: b0_=-1 ! No magnetic field
integer,parameter:: b^C_=-1 ! No magnetic field

{#IFDEF ENERGY
INTEGER,PARAMETER:: p_=e_, pp_=ee_    ! Primitive variables
}

INTEGER,PARAMETER:: v0_=m0_, v^C_=m^C_

{#IFDEF TRACER
integer,parameter:: tr^FL_=Dtr^FL_
}

INTEGER,PARAMETER:: mr_=m0_+r_,mphi_=m0_+phi_,mz_=m0_+z_  ! Polar var. names
{#IFDEF DUST
integer,parameter:: v1d^DS_=rhod^NDS_+^DS
{^NOONEC
integer,parameter:: v2d^DS_=m1d^NDS_+^DS }
{^IFTHREEC
integer,parameter:: v3d^DS_=m2d^NDS_+^DS }
INTEGER,PARAMETER:: mrd^DS_=m0d^DS_+(^NDS*r_)              ! Polar var. names
INTEGER,PARAMETER:: mphid^DS_=m0d^DS_+(^NDS*phi_)          ! Polar var. names
INTEGER,PARAMETER:: mzd^DS_=m0d^DS_+(^NDS*z_)              ! Polar var. names
}

! Note: dust related velocity vectors not handled here
integer, parameter :: nvector=1                             ! No. vector vars
integer, dimension(nvector), parameter :: iw_vector=(/ m0_ /)

{#IFDEF ENERGY
                                                            ! Characteristic
INTEGER,PARAMETER:: soundRW_=1,soundLW_=2,entropW_=3,shearW0_=3      ! waves
INTEGER,PARAMETER:: nworkroe=3      
}

{#IFDEF ISO
                                                            ! Characteristic
INTEGER,PARAMETER:: soundRW_=1,soundLW_=2,shearW0_=2        ! waves
INTEGER,PARAMETER:: nworkroe=1
}

{#IFDEF ENERGY
INTEGER,PARAMETER:: gamma_=1,mu_=2,neqpar=2                     ! equation parameters
}
{#IFDEF ISO
INTEGER,PARAMETER:: gamma_=1,adiab_=2,mu_=3,neqpar=3                     ! equation parameters
}
INTEGER,PARAMETER:: nflag_=nw+1
COMMON, INTEGER:: flags(nflag_)
COMMON, DOUBLE PRECISION:: wflags(nflag_)

! xprob: problem box; iprob: problem
COMMON, INTEGER:: iprob
COMMON, DOUBLE PRECISION:: xprob^L

COMMON, DOUBLE PRECISION::{#IFDEF ENERGY smalle, } minrho, minp 
{#IFDEF DUST 
DOUBLE PRECISION:: minrhod, sdust(1:^NDS),dsdust(1:^NDS),rhodust(1:^NDS),mhcgspar,kbcgspar
DOUBLE PRECISION :: Lstar, Tdust
COMMON /DPDUST/ minrhod,sdust,dsdust,rhodust,mhcgspar,kbcgspar,Tdust,Lstar
}

! end include file amrvacpar.t.hd
!==============================================================================
