!##############################################################################
! include amrvacusrpar - srmhdeos23JetLaun

! This file should contain the number of PROBLEM dependent equation parameters,
! the index names for them with values neqpar+1..neqpar+nspecialpar,
! and a string giving the names for the file header. For example:
!
! INTEGER,PARAMETER:: mass_=neqpar+1, nspecialpar=1
! CHARACTER*4,PARAMETER:: specialparname='mass'
!
! By default there are no special parameters

INTEGER,PARAMETER:: rhoin_=neqpar+1,beta_=neqpar+2,vphim_=neqpar+3,Bmax_=neqpar+4,rhoe_=neqpar+5,&
                    rs_=neqpar+6,r1_=neqpar+7,r2_=neqpar+8,nspecialpar=8
CHARACTER*1,PARAMETER:: specialparname=' '

! end include amrvacusrpar - srmhdeos23JetLaun 
!##############################################################################
