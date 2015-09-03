!##############################################################################
! include amrvacusrpar - doublegemmhd

INTEGER,PARAMETER:: sheetl_=neqpar+1,BB0_=sheetl_+1,rhorat_=BB0_+1, &
                    llx_=rhorat_+1,lly_=llx_+1,psi0top_=lly_+1,psi0bot_=psi0top_+1, &
                    nspecialpar=7
CHARACTER*23,PARAMETER:: specialparname='L B0 Rr lx ly psiT psiB'

! end include amrvacusrpar - doublegemmhd
!##############################################################################
