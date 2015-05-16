!##############################################################################
! amrvacusrpar.t.movstar


INCLUDE:amrvacmodules/coolingpar.t


INTEGER,PARAMETER :: nspecialpar=11
CHARACTER*40,PARAMETER ::  &
   specialparname='Mdot1  vw1  T1  &
                   rhoISM vISM TISM Rstar Rwind &
                   Tsc Lsc Mue'

!
! INPUT parameters
!

INTEGER,PARAMETER :: Mdot1_  = neqpar+1
INTEGER,PARAMETER :: vwind1_ = neqpar+2
INTEGER,PARAMETER :: Twind1_ = neqpar+3

INTEGER,PARAMETER :: rhoISM_ = neqpar+4
INTEGER,PARAMETER :: vISM_   = neqpar+5
INTEGER,PARAMETER :: TISM_   = neqpar+6

INTEGER,PARAMETER :: Rstar_ = neqpar+7
INTEGER,PARAMETER :: Rwind_ = neqpar+8

! Scaling

INTEGER,PARAMETER :: Tscale_= neqpar+9
INTEGER,PARAMETER :: Lscale_= neqpar+10
INTEGER,PARAMETER :: Mue_   = neqpar+11



double precision, parameter :: msol   = 1.989D33
double precision, parameter :: year   = 3.15576D7
double precision, parameter :: kboltz = 1.38065D-16
double precision, parameter :: mhydro = 1.6733D-24

! amrvacusrpar.t.movstar
!##############################################################################
