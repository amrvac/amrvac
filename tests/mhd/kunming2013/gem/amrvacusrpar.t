!##############################################################################
! include amrvacusrpar - gemmhd

INTEGER,PARAMETER:: sheetl_=neqpar+1,T0_=sheetl_+1,rhorat_=T0_+1, &
                    llx_=rhorat_+1,lly_=llx_+1,psi0_=lly_+1, &
                    nspecialpar=6
CHARACTER*17,PARAMETER:: specialparname='L T0 Rr lx ly psi'

! end include amrvacusrpar - gemmhd
!##############################################################################
