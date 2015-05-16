!=============================================================================
! include file amrvacpar.t.srmhdiso  -- Inspired from original SRMHD
! Designed by Fabien Casse 04/12/2013 
!

character*8,parameter:: typephys='srmhdiso'            ! VACPHYS module name
{#IFDEF EPSINF
character*13,parameter:: eqparname='gamma csound2 epsfloor rho0floor rho1e Cr'      ! Equation parameter names
}{#IFNDEF EPSINF
character*13,parameter:: eqparname='gamma csound2 Cr'      ! Equation parameter names
}
! flow variables
!=====Conserve variables=====!
integer,parameter:: d_=1
integer,parameter:: s0_=d_
integer,parameter:: s^C_=s0_+^C
integer,parameter:: b0_=s^NC_
integer,parameter:: b^C_=b0_+^C
{#IFDEF GLM
integer,parameter:: psi_=b^NC_+1
integer,parameter:: nwmhd = psi_
}
{#IFNDEF GLM
integer,parameter:: nwmhd=b^NC_
}
!=============================!

{#IFDEF TRACER
!======Scalar tracers=========!
integer,parameter:: Dtr^FL_=nwmhd+^FL
!=============================!
! Number of variables with tracer
integer,parameter:: nwfluxtr=Dtr^NFL_
}
{#IFNDEF TRACER
! Number of variables with tracer
integer,parameter:: nwfluxtr=nwmhd
}
{#IFDEF EPSINF
!======Cutoff energy ========!
integer,parameter:: Depsinf_=nwfluxtr+1
!======Advected electron density ========!
integer,parameter:: Drho0_=nwfluxtr+2
!======Conserved electron density ========!
integer,parameter:: Drho1_=nwfluxtr+3
!======Advected electron density v2======!
integer,parameter::Dn0_ = nwfluxtr+4
!======Conserved electron density v2=====!
integer,parameter::Dn_  = nwfluxtr+5
! Number of variables
integer,parameter:: nwflux=Dn_
!=============================!
}
{#IFNDEF EPSINF
! Number of variables
integer,parameter:: nwflux=nwfluxtr
}

!==== auxiliary variables====!
integer,parameter:: e_=nwflux+1
integer,parameter:: rhos_=e_
integer,parameter:: tau_=e_
integer,parameter:: lfac_=e_+1     ! Lorentz factor
integer,parameter:: xi_=lfac_+1       
!=============================!

! Number of variables
integer,parameter:: nwaux=3 ! tau lfac xi
integer,parameter:: nwextra=0! 
integer,parameter:: nw=nwflux+nwaux+nwextra

!=====Primitive variables=====!
integer,parameter:: rho_=d_
integer,parameter:: pp_=e_
  !---- 3-velocities ----!
integer,parameter:: v0_=s0_
integer,parameter:: v^C_=v0_+^C
  !---- 4-velocities ----!
integer,parameter:: u0_=s0_
integer,parameter:: u^C_=u0_+^C
{#IFDEF TRACER
integer,parameter:: tr^FL_=Dtr^FL_
}
{#IFDEF EPSINF
integer,parameter:: epsinf_=Depsinf_
integer,parameter:: rho0_=Drho0_
integer,parameter:: rho1_=Drho1_
integer,parameter:: n0_=Dn0_
integer,parameter:: n_=Dn_
}

!=============================!
! polar variable names
integer,parameter:: sr_=s0_+r_
integer,parameter:: sphi_=s0_+phi_
integer,parameter:: sz_=s0_+z_
integer,parameter:: vr_=v0_+r_
integer,parameter:: vphi_=v0_+phi_
integer,parameter:: vz_=v0_+z_
integer,parameter:: uz_=v0_+z_
integer,parameter:: ur_=v0_+r_
integer,parameter:: uphi_=v0_+phi_
integer,parameter:: br_=b0_+r_
integer,parameter:: bphi_=b0_+phi_
integer,parameter:: bz_=b0_+z_
integer,parameter:: ee_=e_

integer, parameter :: nvector=2                             ! No. vector vars
integer, dimension(nvector), parameter :: iw_vector=(/ s0_, b0_ /)

integer,parameter:: fastRW_=1,fastLW_=2,slowRW_=3,slowLW_=4 ! Characteristic
integer,parameter:: entroW_=5,diverW_=6,alfvRW_=7,alfvLW_=8 ! waves
integer,parameter:: nworkroe=15

{#IFNDEF EPSINF
integer,parameter:: gamma_=1, csound2_=2, Cr_=3, neqpar=3              ! equation params
}{#IFDEF EPSINF
INTEGER,PARAMETER:: gamma_=1, csound2_=2, epsfloor_=3, rho0floor_=4, rho1e_=5, Cr_=6, neqpar=6     ! equation params
}

INTEGER,PARAMETER:: nflag_=nw+1
COMMON, INTEGER:: flags(nflag_)
COMMON, DOUBLE PRECISION:: wflags(nflag_)

COMMON, DOUBLE PRECISION::minp,minrho,smallxi,smalltau,govergminone
COMMON, DOUBLE PRECISION::limitvalue

! xprob: problem box; iprob: problem
COMMON, INTEGER:: iprob
COMMON, DOUBLE PRECISION:: xprob^L
{#IFDEF GLM COMMON, DOUBLE PRECISION:: storeeqparch }

! end include file amrvacpar.t.srmhd
!=============================================================================
