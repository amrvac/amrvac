!=============================================================================
! include file amrvacpar.t.srmhd

character*5,parameter:: typephys='srrmhd'            ! VACPHYS module name



CHARACTER*9,PARAMETER:: eqparname='gamma adiab kappa eta' !Equation parameter names



! flow variables
!=====Conserve variables=====!
integer,parameter:: d_=1
integer,parameter:: s0_=d_
integer,parameter:: s1_=s0_+1,s2_=s0_+2,s3_=s0_+3
integer,parameter:: e_=s3_+1
integer,parameter:: tau_=e_
integer,parameter:: b0_=e_
integer,parameter:: b1_=b0_+1,b2_=b0_+2,b3_=b0_+3
integer,parameter:: e0_=b3_
integer,parameter:: e1_=e0_+1,e2_=e0_+2,e3_=e0_+3
integer,parameter:: q_=e3_+1
integer,parameter:: phib_=q_+1
integer,parameter:: psi_=phib_+1
integer,parameter:: nwmhd=psi_


! Number of variables with tracer
integer,parameter:: nwfluxtr=nwmhd



! Number of variables
integer,parameter:: nwflux=nwfluxtr

integer,parameter:: lfac_=nwflux+1     ! Lorentz factor
integer,parameter:: xi_=lfac_+1       ! lfac2 Enthalpy
!=============================!

integer,parameter:: nwaux=2
integer,parameter:: nwextra=0
integer,parameter:: nw=nwflux+nwaux+nwextra

!=====Primitive variables=====!
integer,parameter:: rho_=d_
  !---- 3-velocities ----!
integer,parameter:: v0_=s0_
integer,parameter:: v1_=v0_+1,v2_=v0_+2,v3_=v0_+3
  !---- 4-velocities ----!
integer,parameter:: u0_=s0_
integer,parameter:: u1_=u0_+1,u2_=u0_+2,u3_=u0_+3


integer,parameter:: pp_=e_
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
integer,parameter:: er_=e0_+r_
integer,parameter:: ephi_=e0_+phi_
integer,parameter:: ez_=e0_+z_
integer,parameter:: ee_=e_

integer, parameter :: nvector=3                      ! No. vector vars
integer, dimension(nvector), parameter :: iw_vector=(/ s0_, b0_, e0_ /)

integer,parameter:: fastRW_=1,fastLW_=2,slowRW_=3,slowLW_=4 ! Characteristic
integer,parameter:: entroW_=5,diverW_=6,alfvRW_=7,alfvLW_=8 ! waves
integer,parameter:: nworkroe=15



INTEGER,PARAMETER:: gamma_=1, adiab_=2, kappa_=3, eta_=4, neqpar=4 !equation params


INTEGER,PARAMETER:: nflag_=nw+1
INTEGER:: flags(nflag_)
DOUBLE PRECISION:: wflags(nflag_)

DOUBLE PRECISION::minp,minrho,smallxi,smalltau , govergminone
DOUBLE PRECISION::limitvalue

! xprob: problem box; iprob: problem
INTEGER:: iprob
DOUBLE PRECISION:: xprobmin1,xprobmax1

! end include file amrvacpar.t.srmhd
!=============================================================================
COMMON /DOUB/ wflags,minp,minrho,smallxi,smalltau , govergminone,limitvalue,&
   xprobmin1,xprobmax1
COMMON /INTE/ flags,iprob
