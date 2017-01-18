!===============================================================================
! include file amrvacpar.t.nonlinear

CHARACTER*3,PARAMETER:: typephys='rho'            ! VACPHYS module name

CHARACTER*8,PARAMETER:: eqparname='fluxtype'

integer,parameter:: nwflux=1
integer,parameter:: nwaux=0
integer,parameter:: nwextra=0
integer,parameter:: nw=nwflux+nwaux+nwextra

integer,parameter:: rho_=1
integer,parameter:: b0_=-1 ! No magnetic field
integer,parameter:: b^C_=-1 ! No magnetic field
integer,parameter:: e_=-1 ! No energy (compilation of convert)

INTEGER,PARAMETER:: fluxtype_=1,neqpar=1   ! equation parameters

integer, parameter :: nvector=0                        ! No. vector vars.
integer, dimension(:), allocatable :: iw_vector

integer,parameter:: nworkroe=1

INTEGER,PARAMETER:: nflag_=nw+1
COMMON, INTEGER:: flags(nflag_)
COMMON, DOUBLE PRECISION:: wflags(nflag_)

! xprob: problem box; iprob: problem
COMMON, INTEGER:: iprob
COMMON, DOUBLE PRECISION:: xprob^L
 
! end include file amrvacpar.t.nonlinear
!===============================================================================
