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

INTEGER,PARAMETER:: rho1_=neqpar+1,vel1_=rho1_+1,T1_=vel1_+1, &
                    min_ar_=T1_+1,max_ar_=min_ar_+1,      &
                    nspecialpar=5
CHARACTER*99,PARAMETER:: specialparname='rho1 vel1 T1 min_ar max_ar'

! end include amrvacusrpar - nul
!##############################################################################
