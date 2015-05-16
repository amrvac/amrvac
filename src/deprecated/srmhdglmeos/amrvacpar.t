!=============================================================================
! include file amrvacpar.t.srmhdglmeos

character*8,parameter:: typephys='srmhdglmeos'    ! VACPHYS module name
character*9,parameter:: eqparname='gamma'      ! Equation parameter names

! flow variables
!=====Conserve variables=====!
integer,parameter:: d_=1
integer,parameter:: s0_=d_
integer,parameter:: s^C_=s0_+^C
integer,parameter:: e_=s^NC_+1
integer,parameter:: tau_=e_
integer,parameter:: b0_=e_
integer,parameter:: b^C_=b0_+^C
integer,parameter:: rhos_=e_
integer,parameter:: psi_=b^NC_+1
{^IFMLTFLUID
!======Scalar tracers=========!
integer,parameter:: Dtr^FL_=psi_+^FL
!=============================!
! Number of variables
integer,parameter:: nwflux=Dtr^NFL_
integer,parameter:: lfac_=Dtr^NFL_+1     ! Lorentz factor
}
{^IFONEFLUID
! Number of variables
integer,parameter:: nwflux=psi_
integer,parameter:: lfac_=psi_+1     ! Lorentz factor
}
integer,parameter:: xi_=lfac_+1       ! lfac^2 Enthalpy
!=============================!

integer,parameter:: nwaux=2
integer,parameter:: nwextra=0
integer,parameter:: nw=nwflux+nwaux+nwextra

!=====Primitive variables=====!
integer,parameter:: rho_=d_
  !---- 3-velocities ----!
integer,parameter:: v0_=s0_
integer,parameter:: v^C_=v0_+^C
  !---- 4-velocities ----!
integer,parameter:: u0_=s0_
integer,parameter:: u^C_=u0_+^C

{^IFMLTFLUID
!======Primitive tracers=========!
integer,parameter:: tr^FL_=Dtr^FL_
!=============================!
}

integer,parameter:: pp_=e_,p_=e_
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

! --- to allow compilation of read_relax
integer,parameter:: m0_=s0_,m^C_=s^C_

character(len=10),dimension(1:nw),parameter :: wnamear=(/'d    '&
,'s1   ',{^NOONEC 's2   ',|}{^IFTHREEC 's3   ',|} 'e    '&
,'b1   ',{^NOONEC 'b2   ',|}{^IFTHREEC 'b3   ',|} 'psi  ',&
{^IFMLTFLUID 'Dtr1 ',|}&
{^IFTOFLUID 'Dtr2 ',|}&
{^IFTRFLUID 'Dtr3 ',|}&
{^IFFRFLUID 'Dtr4 ',|}&
{^IFFVFLUID 'Dtr5 ',|}&
{^IFSIFLUID 'Dtr6 ',|}&
{^IFSEFLUID 'Dtr7 ',|}&
{^IFHEFLUID 'Dtr8 ',|}&
{^IFNIFLUID 'Dtr9 ',|}&
{^IFTEFLUID 'Dtr10',|}&
'lfac ','xi   '/)
character(len=10),dimension(1:nw),parameter::primwnamear=(/'rho  '&
,'v1   ',{^NOONEC 'v2   ',|}{^IFTHREEC 'v3   ',|} 'p    '&
,'b1   ',{^NOONEC 'b2   ',|}{^IFTHREEC 'b3   ',|} 'psi  ',&
{^IFMLTFLUID 'tr1  ',|}&
{^IFTOFLUID 'tr2  ',|}&
{^IFTRFLUID 'tr3  ',|}&
{^IFFRFLUID 'tr4  ',|}&
{^IFFVFLUID 'tr5  ',|}&
{^IFSIFLUID 'tr6  ',|}&
{^IFSEFLUID 'tr7  ',|}&
{^IFHEFLUID 'tr8  ',|}&
{^IFNIFLUID 'tr9  ',|}&
{^IFTEFLUID 'tr10 ',|}&
'lfac ','xi   '/)

integer, parameter :: nvector=2                             ! No. vector vars
integer, dimension(nvector), parameter :: iw_vector=(/ s0_, b0_ /)


integer, parameter :: nvectorw=2 *^NC                            ! No. vector vars
integer, dimension(nvectorw), parameter :: iwvector=(/ s^C_, b^C_ /)
COMMON,logical                          :: jac_iw(1:nw)

integer,parameter:: fastRW_=1,fastLW_=2,slowRW_=3,slowLW_=4 ! Characteristic
integer,parameter:: entroW_=5,diverW_=6,alfvRW_=7,alfvLW_=8 ! waves
integer,parameter:: nworkroe=15

integer,parameter:: gamma_=1,Cr_=2,neqpar=2           ! equation params
COMMON, DOUBLE PRECISION:: storeeqparch

COMMON, DOUBLE PRECISION::smalltau,smallxi,minrho,minp
COMMON, DOUBLE PRECISION::limitvalue
INTEGER,PARAMETER:: nflag_=nw+1
COMMON, INTEGER:: flags(nflag_)
COMMON, DOUBLE PRECISION:: wflags(nflag_)

! xprob: problem box; iprob: problem
COMMON, INTEGER:: iprob
COMMON, DOUBLE PRECISION:: xprob^L

! end include file amrvacpar.t.srmhdglmeos
!=============================================================================
