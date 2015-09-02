!##############################################################################
! include amrvacusrpar - promRTideal

INTEGER,PARAMETER:: grav0_=neqpar, grav^D_=grav0_+^D, &
                    eps_=grav0_+^ND+1, nxmodes_=eps_+1, &
                    BB1_=nxmodes_+1,BB3_=BB1_+1,nspecialpar=^ND+4
{^IFTWOD
CHARACTER*21,PARAMETER:: specialparname='g1 g2 eps nxm bb1 bb3'
}
{^IFTHREED
CHARACTER*24,PARAMETER:: specialparname='g1 g2 g3 eps nxm bb1 bb3'
}

INTEGER, PARAMETER:: jmax=80000
COMMON, double precision :: pa(jmax),ra(jmax),ya(jmax),raext(jmax),paext(jmax)
COMMON, double precision :: Lunit,Teunit,nHunit,runit,Bunit,mHunit,k_B,miu0,vunit,tunit,punit,heatunit
COMMON, double precision :: SRadius,dr,gzone,rho0,Tch,Tco,bsca,htra1,htra2,htra3, &
                            ybot,ytop,Tpromax,Tpromin,pwidth,Tcoext

COMMON, DOUBLE PRECISION:: randphase(1000)

! end include amrvacusrpar - promRTideal
!##############################################################################
