!===============================================================================
! include file amrvacpar.t.nonlinear

CHARACTER*3,PARAMETER:: typephys='nonlinear+rho'            ! VACPHYS module name

{^IFONED   CHARACTER*10,PARAMETER:: eqparname='v1 fluxtype'}   ! Equation parameter names
{^IFTWOD   CHARACTER*13,PARAMETER:: eqparname='v1 v2 fluxtype'}
{^IFTHREED CHARACTER*16,PARAMETER:: eqparname='v1 v2 v3 fluxtype'}

integer,parameter:: nwflux=1
integer,parameter:: nwaux=0
integer,parameter:: nwextra=0
integer,parameter:: nw=nwflux+nwaux+nwextra

integer,parameter:: rho_=1
integer,parameter:: b0_=-1 ! No magnetic field
integer,parameter:: b^C_=-1 ! No magnetic field
integer,parameter:: e_=-1 ! No energy (compilation of convert)

INTEGER,PARAMETER:: fluxtype_=1   ! equation parameters
INTEGER,PARAMETER:: v0_=fluxtype_,v^D_=v0_+^D,neqpar=v^ND_   ! equation parameters

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
