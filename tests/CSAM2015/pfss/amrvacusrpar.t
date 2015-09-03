!##############################################################################
! include amrvacusrpar - pfss 

! This file should contain tha number of PROBLEM dependent equation parameters,
! the index names for them with values neqpar+1..neqpar+nspecialpar,
! and a string giving the names for the file header. For example:
!
! INTEGER,PARAMETER :: mass=neqpar+1, nspecialpar=1
! CHARACTER*4,PARAMETER :: specialparname='mass'
!
! By default there are no special parameters

!INTEGER,PARAMETER:: grav0_=neqpar, grav^D_=grav0_+^D, mu_=grav^ND_+1,&
!vdr_=mu_+1,nspecialpar=2+^ND
!CHARACTER*20,PARAMETER:: specialparname='grav1 grav2 grav3 mu vdr'
!COMMON, double precision ::
!Lunit,Teunit,nHunit,runit,Bunit,mHunit,vunit,tunit,punit
!COMMON, double precision :: SRadius,rhob,Tiso,delb2


integer, parameter :: grav0_=neqpar, grav^D_=grav0_+^D, vescape_=grav0_+^ND+1,&
                   nspecialpar=^ND+1
character*80, parameter :: specialparname='g1 g2 g3 vesc'

COMMON, double precision :: Lunit,Teunit,nHunit,runit,punit,mHunit,k_B,miu0,vunit,tunit,Bunit
COMMON, double precision :: rhob,Tiso
! end include amrvacusrpar - pfss
!##############################################################################
