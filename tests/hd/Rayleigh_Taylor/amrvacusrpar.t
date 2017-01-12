!##############################################################################
! include amrvacusrpar - gravity

INTEGER,PARAMETER:: grav0_=neqpar, grav^D_=neqpar+^D, nspecialpar=^ND
{^IFONED   CHARACTER*5 ,PARAMETER:: specialparname='grav1'}
{^IFTWOD   CHARACTER*11,PARAMETER:: specialparname='grav1 grav2'}
{^IFTHREED CHARACTER*17,PARAMETER:: specialparname='grav1 grav2 grav3'}

! end include amrvacusrpar - gravity
!##############################################################################
