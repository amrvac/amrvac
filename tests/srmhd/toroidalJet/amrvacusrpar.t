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

INTEGER,PARAMETER:: bm_=neqpar+1, rm_=neqpar+2, rj_=neqpar+3, lor_=neqpar+4, &
     rhoe_=neqpar+5, betam_=neqpar+6, nspecialpar=6
CHARACTER*1,PARAMETER:: specialparname=' '

! end include amrvacusrpar - nul
!##############################################################################
