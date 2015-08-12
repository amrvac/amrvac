!=============================================================================
! include file amrvacpar.t.ff

character*5,parameter:: typephys='ff'            ! VACPHYS module name


CHARACTER*9,PARAMETER:: eqparname='kpar kperp kappa'    ! Equation parameter names

! flow variables
!=====Conserve variables=====!
integer,parameter:: b0_=0
integer,parameter:: b^C_=b0_+^C
integer,parameter:: e0_=b^NC_
integer,parameter:: e^C_=e0_+^C
integer,parameter:: phib_=e^NC_+1

! Number of variables
integer,parameter:: nwflux=phib_

integer,parameter:: nwaux=0
integer,parameter:: nwextra=0
integer,parameter:: nw=nwflux+nwaux+nwextra

!=====Primitive variables=====!
! polar variable names
integer,parameter:: br_=b0_+r_
integer,parameter:: bphi_=b0_+phi_
integer,parameter:: bz_=b0_+z_
integer,parameter:: er_=e0_+r_
integer,parameter:: ephi_=e0_+phi_
integer,parameter:: ez_=e0_+z_

!====== Dummies ====!
integer,parameter:: e_=1

integer, parameter :: nvector=2                      ! No. vector vars
integer, dimension(nvector), parameter :: iw_vector=(/ b0_, e0_ /)

integer,parameter:: fastRW_=1,fastLW_=2,slowRW_=3,slowLW_=4 ! Characteristic
integer,parameter:: entroW_=5,diverW_=6,alfvRW_=7,alfvLW_=8 ! waves
integer,parameter:: nworkroe=15

INTEGER,PARAMETER:: kpar_=1, kperp_=2, kappa_=3, neqpar=3     ! equation params

INTEGER,PARAMETER:: nflag_=nw+1
COMMON, INTEGER:: flags(nflag_)
COMMON, DOUBLE PRECISION:: wflags(nflag_)

! xprob: problem box; iprob: problem
COMMON, INTEGER:: iprob
COMMON, DOUBLE PRECISION:: xprob^L

! end include file amrvacpar.t.srmhd
!=============================================================================
