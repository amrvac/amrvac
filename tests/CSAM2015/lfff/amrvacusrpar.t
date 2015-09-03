!##############################################################################
! include amrvacusrpar - arfff

INTEGER,PARAMETER:: grav0_=neqpar, grav^D_=grav0_+^D, nspecialpar=^ND
CHARACTER*20,PARAMETER:: specialparname='grav1 grav2 grav3'
COMMON, double precision :: Lunit,Teunit,nHunit,runit,Bunit,mHunit,vunit,tunit,punit
COMMON, double precision :: SRadius,rhob,Tiso
! end include amrvacusrpar - arfff
!##############################################################################
