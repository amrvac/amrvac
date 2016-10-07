!##############################################################################
! include amrvacusrpar - solaratmosphere23

INTEGER,PARAMETER:: grav0_=neqpar, grav^D_=grav0_+^D, kappa_=neqpar+^ND+1,&
    Tscale_=kappa_+1,Lscale_=Tscale_+1,Mue_=Lscale_+1,nspecialpar=^ND+4
CHARACTER*29,PARAMETER:: specialparname='grav1 grav2 kappa Tsc Lsc Mue'
INTEGER, PARAMETER:: jmax=80000
COMMON, double precision :: pa(jmax),ra(jmax),ya(jmax),dya
COMMON, double precision :: Lunit,Teunit,nHunit,runit,punit,mHunit,k_B,miu0,vunit,tunit,heatunit,Bunit
COMMON, double precision :: gzone,B0,theta,SRadius,kx,ly,bQ0
! end include amrvacusrpar - solaratmosphere23
!##############################################################################
