!##############################################################################
! include amrvacusrpar - prombbs

INTEGER,PARAMETER:: grav0_=neqpar, grav^D_=grav0_+^D, nspecialpar=^ND
CHARACTER*20,PARAMETER:: specialparname='grav1 grav2'
INTEGER, PARAMETER:: jmax=80000
COMMON, double precision :: pa(jmax),ra(jmax),ya(jmax)
COMMON, double precision :: Lunit,Teunit,nHunit,runit,Bunit,mHunit,k_B,miu0,vunit,tunit,punit,heatunit
COMMON, double precision :: SRadius,dr,gzone,rho0,Tch,Tco,B0,bsca,htra2
! end include amrvacusrpar - prombbs
!##############################################################################
