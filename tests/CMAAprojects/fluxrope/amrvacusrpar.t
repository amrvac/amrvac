!##############################################################################
! include amrvacusrpar - fluxrope23

INTEGER,PARAMETER:: grav0_=neqpar, grav^D_=grav0_+^D,LB0_=neqpar+^ND+1,nspecialpar=^ND+1
CHARACTER*28,PARAMETER:: specialparname='grav1 grav2 LB0'
INTEGER, PARAMETER:: jmax=80000
COMMON, double precision :: pa(jmax),ra(jmax),ya(jmax),dya
COMMON, double precision :: Lunit,Teunit,nHunit,runit,punit,mHunit,k_B,miu0,vunit,tunit,heatunit,Bunit
COMMON, double precision :: Tpho,rpho,Ttop,gzone,B0,theta,tdstop
! end include amrvacusrpar - fluxrope23
!##############################################################################
