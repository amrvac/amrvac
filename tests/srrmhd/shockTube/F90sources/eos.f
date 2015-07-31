!###########################################################################
! module amrvacphys - srmhdeos version april 2009
! This module is developed using the paper Meliani et al 2004
!===========================================================================
subroutine Enthalpy(w,ixImin1,ixImax1,ixOmin1,ixOmax1,patchw,rhoh)

!================== IMPORTANT ==================!
!This subroutine is used only with primitive variables on input w
!===============================================!

include 'amrvacdef.f'

integer, intent(in)                                    :: ixImin1,ixImax1,&
   ixOmin1,ixOmax1
double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)    :: w
logical,          dimension(ixGlo1:ixGhi1), intent(in)         :: patchw
double precision, dimension(ixGlo1:ixGhi1), intent(out)        :: rhoh

!--------------------------------------! 




where(.not.patchw(ixOmin1:ixOmax1))
 rhoh(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,rho_) + govergminone &
    * w(ixOmin1:ixOmax1,pp_)
end where




end subroutine Enthalpy
!=============================================================================
subroutine Pressuren(w,ixImin1,ixImax1,ixOmin1,ixOmax1,varconserve,p,patchw)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!   both cases assume updated auxiliary variables xi_ en lfac_
!===============================================!

include 'amrvacdef.f'

integer,intent(in)                                :: ixImin1,ixImax1,ixOmin1,&
   ixOmax1
double precision, intent(in),dimension(ixImin1:ixImax1,1:nw):: w
logical,intent(in)                                :: varconserve
logical,intent(in),dimension(ixGlo1:ixGhi1)               :: patchw
double precision,intent(out), dimension(ixGlo1:ixGhi1)    :: p

!--------------------------------------! 



if (varconserve) then
 where(.not.patchw(ixOmin1:ixOmax1))
    p(ixOmin1:ixOmax1) = ( w(ixOmin1:ixOmax1,xi_)/(w(ixOmin1:ixOmax1,lfac_)&
       **2.0d0) - w(ixOmin1:ixOmax1,d_)/w(ixOmin1:ixOmax1,lfac_) ) &
       / govergminone
 end where
else
 where(.not.patchw(ixOmin1:ixOmax1))
    p(ixOmin1:ixOmax1) = ( w(ixOmin1:ixOmax1,xi_)/(w(ixOmin1:ixOmax1,lfac_)&
       **2.0d0) - w(ixOmin1:ixOmax1,rho_) ) / govergminone
 end where




end if

return
end subroutine Pressuren
!============================================================================
subroutine Pressure(w,ixImin1,ixImax1,ixOmin1,ixOmax1,varconserve,p)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!   both cases assume updated auxiliary variables xi_ en lfac_
!===============================================!

include 'amrvacdef.f'

integer,intent(in)                                :: ixImin1,ixImax1,ixOmin1,&
   ixOmax1
double precision, intent(in),dimension(ixImin1:ixImax1,1:nw):: w
logical,intent(in)                                :: varconserve
double precision,intent(out), dimension(ixGlo1:ixGhi1)    :: p

!--------------------------------------! 



if (varconserve) then
    p(ixOmin1:ixOmax1) = (w(ixOmin1:ixOmax1,xi_)/(w(ixOmin1:ixOmax1,lfac_)&
       **2.0d0)-(w(ixOmin1:ixOmax1,d_)/w(ixOmin1:ixOmax1,lfac_)))/govergminone
else
    p(ixOmin1:ixOmax1) = (w(ixOmin1:ixOmax1,xi_)/(w(ixOmin1:ixOmax1,lfac_)&
       **2.0d0)-w(ixOmin1:ixOmax1,rho_))/govergminone
end if




return
end subroutine Pressure
!============================================================================
subroutine getinvcsound2(w,ixImin1,ixImax1,ixOmin1,ixOmax1,invcsound2)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w
!===============================================!

include 'amrvacdef.f'

integer,intent(in)                                    :: ixImin1,ixImax1,&
   ixOmin1,ixOmax1
double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)   :: w
double precision, dimension(ixGlo1:ixGhi1), intent(out)       :: invcsound2
double precision, dimension(ixGlo1:ixGhi1)                    :: h

!--------------------------------------!  

h(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,xi_)/(w(ixOmin1:ixOmax1,lfac_)**2.0d0)



invcsound2(ixOmin1:ixOmax1) =  (w(ixOmin1:ixOmax1,xi_)/(w(ixOmin1:ixOmax1,&
   xi_) - w(ixOmin1:ixOmax1,lfac_)*w(ixOmin1:ixOmax1,d_))) &
               /(eqpar(gamma_)-one)




end subroutine getinvcsound2
!=============================================================================
subroutine getcsound2(w,ixImin1,ixImax1,ixOmin1,ixOmax1,varconserve,csound2)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!===============================================!

include 'amrvacdef.f'

integer,intent(in)                                   :: ixImin1,ixImax1,&
   ixOmin1,ixOmax1
double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)  :: w
logical, intent(in)                                  :: varconserve
double precision, dimension(ixGlo1:ixGhi1),intent(out)       :: csound2
double precision, dimension(ixGlo1:ixGhi1)                   :: h, P

!--------------------------------------!

h(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,xi_)/(w(ixOmin1:ixOmax1,lfac_)**2.0d0)


call Pressure(w,ixImin1,ixImax1,ixOmin1,ixOmax1,varconserve,P)
csound2(ixOmin1:ixOmax1) = (eqpar(gamma_) * P(ixOmin1:ixOmax1)) &
   / h(ixOmin1:ixOmax1)



end subroutine getcsound2
!=============================================================================
subroutine FuncPressure(xicurrent,lfac,d,ssqr,tau,dlfacdxi,p,dpdxi)

! compute pointwise value for pressure p and dpdxi

include 'amrvacdef.f'
  
double precision, intent(in)         :: xicurrent,lfac,d,ssqr,tau,dlfacdxi
double precision, intent(out)        :: p,dpdxi
! .. local ..
double precision                     :: rho,h

double precision                     :: dpdchi, dpdrho 


!_______________________________________________________________!


h=xicurrent/(lfac**2.0d0)
rho=d/lfac




p = (h - rho)/govergminone
dpdchi = one/govergminone
dpdrho = zero





! This is quite general.  Whenever you can give expressions 
 !for dpdchi and dpdrho from your EOS, this will do the inversion for your EOS.
dpdxi = dpdchi * (one/lfac**2 + ((d*lfac-2.0d0*xicurrent)/lfac&
   **3.0d0) * dlfacdxi) &
      - dpdrho * d/lfac**2 * dlfacdxi

return
end subroutine FuncPressure
!=============================================================================!
subroutine smallvaluesEOS
! This is the smallvalues for Synge gas to be precise.  
include 'amrvacdef.f'

double precision::Lsmallrho,Lsmallp,LsmallE
!--------------------------------------------------
Lsmallrho=(1.0d0 + 10.0d0 * minrho) * minrho
Lsmallp=(1.0d0 + 10.0d0 * minp) * minp
LsmallE=Lsmallp/(eqpar(gamma_)-one)+    dsqrt((Lsmallp/(eqpar(gamma_)&
   -one))**2.0d0+Lsmallrho**2.0d0)

smalltau=half*((eqpar(gamma_)+one)*LsmallE-(eqpar(gamma_)-one)*Lsmallrho&
   **2.0d0/LsmallE)-Lsmallp-Lsmallrho
! may need to replace by smallxi above

smallxi=half*((eqpar(gamma_)+one)*LsmallE-(eqpar(gamma_)-one)*Lsmallrho&
   **2.0d0/LsmallE)

end Subroutine smallvaluesEOS
!=============================================================================!
subroutine xinoFlow(xi,tau,d,bb)

include 'amrvacdef.f'

double precision:: xi,tau,d,bb
!_______________________________________________________________!



xi = eqpar(gamma_) *(tau - half * bb)+d


end subroutine xinoFlow

!=============================================================================
subroutine FuncEnthalpy(pcurrent,lfac2inv,d,s1,s2,s3,tau,sqrs,xicurrent,&
   dv2d2p,h,dhdp,ierror)

include 'amrvacdef.f'
  
double precision, intent(in)             :: pcurrent,lfac2inv,d,s1,s2,s3,tau,&
   sqrs,xicurrent,dv2d2p
double precision, intent(out)            :: h,dhdp
integer, intent(inout)                   :: ierror

double precision                         :: rho

!_______________________________________________________________!
rho=d*dsqrt(lfac2inv)




h = rho + govergminone * pcurrent
dhdp = govergminone + d/sqrt(lfac2inv)*sqrs/xicurrent**3.0d0



return
end subroutine FuncEnthalpy
!=========================================================================
subroutine Bisection_Enthalpy(pcurrent,lfac2inv,d,s1,s2,s3,tau,sqrs,xicurrent,&
   h,ierror)

include 'amrvacdef.f'
  
integer:: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision:: pcurrent,lfac2inv,d,s1,s2,s3,tau,sqrs,xicurrent
double precision:: h
integer::ierror
double precision:: rho

!_______________________________________________________________!
rho=d*dsqrt(lfac2inv)



h = rho + govergminone * pcurrent


return
end subroutine Bisection_Enthalpy

!=============================================================================
subroutine Calcule_Geffect(w,ixImin1,ixImax1,ixOmin1,ixOmax1,varconserve,&
   Geffect)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!   both cases assume updated auxiliary variables xi_ en lfac_
!===============================================!

include 'amrvacdef.f'

integer:: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision,intent(in):: w(ixImin1:ixImax1,1:nw)
logical,intent(in)         :: varconserve
double precision,intent(out):: Geffect(ixGlo1:ixGhi1)

!_______________________________________________________________!


Geffect(ixOmin1:ixOmax1) = eqpar(gamma_)



end subroutine Calcule_Geffect
!=============================================================================
! end module amrvacphys - srmhdeos
!##############################################################################
