!##############################################################################
! include amrvacusrpar - jetCloudInteraction

! This file should contain the number of PROBLEM dependent equation parameters,
! the index names for them with values neqpar+1..neqpar+nspecialpar,
! and a string giving the names for the file header. 

INTEGER,PARAMETER:: beta_=neqpar+1,eta_=beta_+1,ca_=eta_+1,Ma_=ca_+1, &
                    rc_=Ma_+1,nspecialpar=5
CHARACTER*17,PARAMETER:: specialparname='beta eta ca Ma rc'

! end include amrvacusrpar - jetCloudInteraction
!##############################################################################
