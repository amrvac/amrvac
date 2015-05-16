!##############################################################################
! include amrvacusrpar - mhdconvection

INTEGER,PARAMETER:: grav0_=neqpar, grav^D_=grav0_+^D, kappa_=neqpar+^ND+1,&
    mu_=kappa_+1,bstr_=mu_+1,temptop_=bstr_+1,nspecialpar=^ND+4

{^IFTWOD
CHARACTER*24,PARAMETER:: specialparname='g1 g2 kappa mu Bstr Ttop'
}
{^IFTHREED
CHARACTER*27,PARAMETER:: specialparname='g1 g2 g3 kappa mu Bstr Ttop'
}

! end include amrvacusrpar - mhdconvection
!##############################################################################
