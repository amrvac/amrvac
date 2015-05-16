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

INTEGER,PARAMETER:: beta_=neqpar+1,eta_=beta_+1,ca_=eta_+1,Ma_=ca_+1,nspecialpar=4
CHARACTER*14,PARAMETER:: specialparname='beta eta ca Ma'

! end include amrvacusrpar - nul
!##############################################################################
