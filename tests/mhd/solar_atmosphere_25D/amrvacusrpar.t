!##############################################################################
! include amrvacusrpar - solaratmosphere23

INCLUDE:amrvacmodules/coolingpar.t

INTEGER,PARAMETER:: grav0_=neqpar, grav^D_=grav0_+^D, Tscale_=neqpar+^ND+1,&
    Lscale_=Tscale_+1,Mue_=Lscale_+1,nspecialpar=^ND+3
CHARACTER*29,PARAMETER:: specialparname='grav1 grav2 Tsc Lsc Mue'
INTEGER, PARAMETER:: jmax=80000
COMMON, double precision :: pa(jmax),ra(jmax),ya(jmax),dya
COMMON, double precision :: heatunit,gzone,B0,theta,SRadius,kx,ly,bQ0
! end include amrvacusrpar - solaratmosphere23
!##############################################################################
