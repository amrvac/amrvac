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

INTEGER,PARAMETER:: beta0_ = neqpar+1, rho0_= neqpar+2, vcoll_=neqpar+3, nspecialpar=3
CHARACTER*1,PARAMETER:: specialparname=' '

! end include amrvacusrpar - nul
!##############################################################################
