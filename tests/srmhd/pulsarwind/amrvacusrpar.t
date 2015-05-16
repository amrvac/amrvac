!##############################################################################
! include amrvacusrpar - nul

! This file should contain the number of PROBLEM dependent equation parameters,
! the index names for them with values neqpar+1..neqpar+nspecialpar,
! and a string giving the names for the file header. For example:
!
! INTEGER,PARAMETER:: mass_=neqpar+1, nspecialpar=1
! CHARACTER*4,PARAMETER:: specialparname='mass'
!
! By default there are no special parameters

INTEGER,PARAMETER:: sig0_=neqpar+1, lor_=neqpar+2, &
     b_=neqpar+3, eta_=neqpar+4, ve_=neqpar+5, rhoe_=neqpar+6, ri_=neqpar+7, &
     re_=neqpar+8, &
     theta0_=neqpar+9, flfac_=neqpar+10, flfacNeighbor_=neqpar+11, &
     m_=neqpar+12, &
     DeltaRho_=neqpar+13, &
     nspecialpar=13
COMMON, INTEGER      :: addlevel, addlevel_origin
CHARACTER*1,PARAMETER:: specialparname=' '

! end include amrvacusrpar - nul
!##############################################################################
