!##############################################################################
! include amrvacusrpar - shockcloud

INTEGER,PARAMETER:: chi_=neqpar+1,nval_=chi_+1,rcore_=nval_+1,rbound_=rcore_+1, & 
            xshock_=rbound_+1,beta_=xshock_+1,machs_=beta_+1,nspecialpar=7

CHARACTER*29,PARAMETER:: specialparname='chi nval rc rb xsh beta machs'

! end include amrvacusrpar - shockcloud
!##############################################################################
