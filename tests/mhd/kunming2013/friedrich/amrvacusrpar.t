!##############################################################################
! include amrvacusrpar - friedrich

! This file should contain the number of PROBLEM dependent equation parameters,
! the index names for them with values neqpar+1..neqpar+nspecialpar,
! and a string giving the names for the file header. 

INTEGER,PARAMETER:: bxini_=neqpar+1, dvz_=neqpar+2, dp_=neqpar+3, nspecialpar=3
CHARACTER*13,PARAMETER:: specialparname='bx delvz delp'

! end include amrvacusrpar - friedrich
!##############################################################################
