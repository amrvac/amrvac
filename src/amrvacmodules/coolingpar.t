!==============================================================================
! include file coolingpar.t
! this can be used for radiative losses in HD and MHD
! 02.09.2012: moved to amrvacmodules folder and renamed from amrvacusrpar.cooling.t
! by Oliver Porth
!-----------------------------------------------------------------------------
!
!(1)write in your usr file the line
!
!INCLUDE:amrvacmodules/cooling.t
!
! and put a call to addsource_cooling in specialsource (amrvacnul/specialsource.t)
!         a call to getdt_cooling in getdt_special     (amrvacnul/specialsource.t)
!         a call to coolinit in initglobaldata_usr     (amrvacnul/specialini.t)
!
!(2) write in your usrpar file the line which includes this file
!
!INCLUDE:amrvacmodules/coolingpar.t
!
!(3) edit in your usrpar file the specialpar and specialparname 
!    entries to include for sure:
!
!INTEGER,PARAMETER :: Tscale_=neqpar+1                   ! dimension: kboltz/mhydro *1/velocity^2 
!               
!INTEGER,PARAMETER :: Lscale_=Tscale_+1                  ! dimension: (1/mhydro^2) * (mass*time^3)/(length^5)
!
!INTEGER,PARAMETER :: Mue_ = Lscale_+1,nspecialpar=3     ! dimensionless
!
!CHARACTER*11,PARAMETER :: specialparname='Tsc Lsc Mue'
!-----------------------------------------------------------------------------


! parameters used for implicit cooling source calculations
INTEGER, PARAMETER          :: maxiter = 100
DOUBLE PRECISION, PARAMETER :: e_error = 1.0D-6

! Radiation related variables
integer, parameter :: ncoolmax=100000

COMMON,DOUBLE PRECISION :: tcool(1:ncoolmax), Lcool(1:ncoolmax), dLdtcool(1:ncoolmax)
COMMON,DOUBLE PRECISION :: tcoolmin,tcoolmax
COMMON,DOUBLE PRECISION :: tcmulti(1:3,1:ncoolmax), Lcmulti(1:3,1:ncoolmax)
COMMON,DOUBLE PRECISION :: dLcmulti(1:3,1:ncoolmax)

COMMON,DOUBLE PRECISION :: Yc(1:ncoolmax), invYc(1:ncoolmax), Tref, Lref

! end include file coolingpar.t
!==============================================================================
