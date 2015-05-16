!##############################################################################
! include amrvacusrpar - randomplasma

INTEGER,PARAMETER:: delrho_=neqpar+1,nxmodes_=delrho_+1, &
        Lx0_=nxmodes_+1,a0_=Lx0_+1,sig0_=a0_+1,x0_=sig0_+1,beta_=x0_+1,  &
        nspecialpar=7
CHARACTER*30,PARAMETER:: specialparname='delrho nxm Lx0 a0 sig0 x0 beta'
COMMON, DOUBLE PRECISION:: randphasey(1000),randampliy(1000),rho_std0,rho_std,rho_mean0
COMMON, LOGICAL:: norm_rho 
COMMON, double precision :: Lunit,runit,punit,k_B,miu0,vunit,tunit,Bunit

! end include amrvacusrpar - randomplasma
!##############################################################################
