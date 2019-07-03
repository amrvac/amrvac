!##############################################################################
! amrvacusrpar.t.testrad


INCLUDE:amrvacmodules/coolingpar.t


INTEGER,PARAMETER :: nspecialpar=11
CHARACTER*40,PARAMETER ::  &
   specialparname='Tsc Lsc Mue  rhoISM TISM  Fs alphah md1 md2 v1 v2'



!
! SCaling factors
!
INTEGER,PARAMETER :: Tscale_ = neqpar+1
INTEGER,PARAMETER :: Lscale_ = neqpar+2
INTEGER,PARAMETER :: Mue_    = neqpar+3


!
! INPUT parameters
!


INTEGER,PARAMETER :: rhoISM_  = neqpar+4
INTEGER,PARAMETER :: TISM_    = neqpar+5
INTEGER,PARAMETER :: Fstar_   = neqpar+6
INTEGER,PARAMETER :: alphah_  = neqpar+7

INTEGER,PARAMETER  :: Mdot1_  = neqpar+8
INTEGER,PARAMETER :: Mdot2_   = neqpar+9
INTEGER,PARAMETER :: vel1_    = neqpar+10
INTEGER,PARAMETER :: vel2_    = neqpar+11



double precision, parameter :: msol = 1.989D33
double precision, parameter :: year = 3.15576D7
double precision, parameter :: kboltz = 1.38065D-16
double precision, parameter :: mhydro = 1.6733D-24



! amrvacusrpar.t.windwind
!##############################################################################
