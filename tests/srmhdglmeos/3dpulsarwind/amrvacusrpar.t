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

INTEGER,PARAMETER:: f0_=neqpar+1, sig0_=neqpar+2, lor_=neqpar+3, &
      xi0_=neqpar+4, ve_=neqpar+5, rhoe_=neqpar+6, r0_=neqpar+7, & 
      rc_=neqpar+8, nspecialpar=8
CHARACTER*1,PARAMETER:: specialparname=' '

! end include amrvacusrpar - nul
!##############################################################################
