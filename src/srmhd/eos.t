!###########################################################################
! module amrvacphys - srmhdeos version april 2009
! This module is developed using the paper Meliani et al 2004
!===========================================================================
subroutine Enthalpy(w,ixI^L,ixO^L,patchw,rhoh)

!================== IMPORTANT ==================!
!This subroutine is used only with primitive variables on input w
!===============================================!

include 'amrvacdef.f'

integer, intent(in)                                    :: ixI^L,ixO^L
double precision, dimension(ixI^S,1:nw), intent(in)    :: w
logical,          dimension(ixG^T), intent(in)         :: patchw
double precision, dimension(ixG^T), intent(out)        :: rhoh
{#IFDEF SYNGE
double precision, dimension(ixG^T)                     :: E_Th,E }
!--------------------------------------! 

{#IFDEF SYNGE
where(.not.patchw(ixO^S))
 E_Th(ixO^S) = w(ixO^S,pp_)/(eqpar(gamma_)-one) 
 E(ixO^S) = E_Th(ixO^S) + dsqrt(E_Th(ixO^S)**2.0d0+w(ixO^S,rho_)**2.0d0)
 rhoh(ixO^S) = half*((eqpar(gamma_)+one) * E(ixO^S)-&
	       (eqpar(gamma_)-one)* w(ixO^S,rho_)*(w(ixO^S,rho_)/E(ixO^S)))
end where
}

{#IFDEF GAMMA
where(.not.patchw(ixO^S))
 rhoh(ixO^S) = w(ixO^S,rho_) + govergminone * w(ixO^S,pp_)
end where
}

{#IFDEF ISO
where(.not.patchw(ixO^S))
 rhoh(ixO^S) = w(ixO^S,rho_) + govergminone * eqpar(adiab_)*w(ixO^S,rho_)**eqpar(gamma_)
end where
}

end subroutine Enthalpy
!=============================================================================
subroutine Pressuren(w,ixI^L,ixO^L,varconserve,p,patchw)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!   both cases assume updated auxiliary variables xi_ en lfac_
!===============================================!

include 'amrvacdef.f'

integer,intent(in)                                :: ixI^L,ixO^L
double precision, intent(in),dimension(ixI^S,1:nw):: w
logical,intent(in)                                :: varconserve
logical,intent(in),dimension(ixG^T)               :: patchw
double precision,intent(out), dimension(ixG^T)    :: p
{#IFDEF SYNGE
double precision, dimension(ixG^T)     :: h,E }
!--------------------------------------! 
{#IFDEF SYNGE
if (varconserve) then
 where(.not.patchw(ixO^S))
  h(ixO^S) = w(ixO^S,xi_)/(w(ixO^S,lfac_)**2.0d0)
  E(ixO^S) = (h(ixO^S)+dsqrt(h(ixO^S)**2.0d0+(eqpar(gamma_)**2.0d0-one)&
             *(w(ixO^S,d_)/w(ixO^S,lfac_))**2.0d0))/(eqpar(gamma_)+one)
  p(ixO^S) = (eqpar(gamma_)-one)/2.0d0&
             * (E(ixO^S)-((w(ixO^S,d_)/w(ixO^S,lfac_))**2.0d0/E(ixO^S)))
 end where
else
 where(.not.patchw(ixO^S))
  h(ixO^S) = w(ixO^S,xi_)/(w(ixO^S,lfac_)**2.0d0)
  E(ixO^S) = (h(ixO^S)+dsqrt(h(ixO^S)**2.0d0+(eqpar(gamma_)**2.0d0-one)&
             *w(ixO^S,rho_)**2.0d0))/(eqpar(gamma_)+one)
  p(ixO^S) = (eqpar(gamma_)-one)/2.0d0&
             * (E(ixO^S)-(w(ixO^S,rho_)**2.0d0/E(ixO^S)))
 end where
}

{#IFDEF GAMMA
if (varconserve) then
 where(.not.patchw(ixO^S))
    p(ixO^S) = ( w(ixO^S,xi_)/(w(ixO^S,lfac_)**2.0d0) - w(ixO^S,d_)/w(ixO^S,lfac_) ) / govergminone
 end where
else
 where(.not.patchw(ixO^S))
    p(ixO^S) = ( w(ixO^S,xi_)/(w(ixO^S,lfac_)**2.0d0) - w(ixO^S,rho_) ) / govergminone
 end where
}

{#IFDEF ISO
if (varconserve) then
 where(.not.patchw(ixO^S))
    p(ixO^S) = eqpar(adiab_) * (w(ixO^S,d_)/w(ixO^S,lfac_))**eqpar(gamma_)
 end where
else
 where(.not.patchw(ixO^S))
    p(ixO^S) = eqpar(adiab_) * w(ixO^S,rho_)**eqpar(gamma_)
 end where
}

end if

return
end subroutine Pressuren
!============================================================================
subroutine Pressure(w,ixI^L,ixO^L,varconserve,p)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!   both cases assume updated auxiliary variables xi_ en lfac_
!===============================================!

include 'amrvacdef.f'

integer,intent(in)                                :: ixI^L,ixO^L
double precision, intent(in),dimension(ixI^S,1:nw):: w
logical,intent(in)                                :: varconserve
double precision,intent(out), dimension(ixG^T)    :: p
{#IFDEF SYNGE
double precision, dimension(ixG^T)     :: h,E }
!--------------------------------------! 
{#IFDEF SYNGE
if (varconserve) then
  h(ixO^S) = w(ixO^S,xi_)/(w(ixO^S,lfac_)**2.0d0)
  E(ixO^S) = (h(ixO^S)+dsqrt(h(ixO^S)**2.0d0+(eqpar(gamma_)**2.0d0-one)&
             *(w(ixO^S,d_)/w(ixO^S,lfac_))**2.0d0))/(eqpar(gamma_)+one)
  p(ixO^S) = (eqpar(gamma_)-one)/2.0d0&
             * (E(ixO^S)-((w(ixO^S,d_)/w(ixO^S,lfac_))**2.0d0/E(ixO^S)))
else
  h(ixO^S) = w(ixO^S,xi_)/(w(ixO^S,lfac_)**2.0d0)
  E(ixO^S) = (h(ixO^S)+dsqrt(h(ixO^S)**2.0d0+(eqpar(gamma_)**2.0d0-one)&
             *w(ixO^S,rho_)**2.0d0))/(eqpar(gamma_)+one)
  p(ixO^S) = (eqpar(gamma_)-one)/2.0d0&
             * (E(ixO^S)-(w(ixO^S,rho_)**2.0d0/E(ixO^S)))
end if
}

{#IFDEF GAMMA
if (varconserve) then
    p(ixO^S) = (w(ixO^S,xi_)/(w(ixO^S,lfac_)**2.0d0)-(w(ixO^S,d_)/w(ixO^S,lfac_)))/govergminone
else
    p(ixO^S) = (w(ixO^S,xi_)/(w(ixO^S,lfac_)**2.0d0)-w(ixO^S,rho_))/govergminone
end if
}

{#IFDEF ISO
if (varconserve) then
    p(ixO^S) = eqpar(adiab_) * (w(ixO^S,d_)/w(ixO^S,lfac_))**eqpar(gamma_)
else
    p(ixO^S) = eqpar(adiab_) * w(ixO^S,rho_)**eqpar(gamma_)
end if
}

return
end subroutine Pressure
!============================================================================
subroutine getinvcsound2(w,ixI^L,ixO^L,invcsound2)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w
!===============================================!

include 'amrvacdef.f'

integer,intent(in)                                    :: ixI^L,ixO^L
double precision, dimension(ixI^S,1:nw), intent(in)   :: w
double precision, dimension(ixG^T), intent(out)       :: invcsound2
double precision, dimension(ixG^T)                    :: h
{#IFDEF SYNGE
double precision, dimension(ixG^T)                    :: rho,E,P }
!--------------------------------------!  

h(ixO^S)=w(ixO^S,xi_)/(w(ixO^S,lfac_)**2.0d0)
{#IFDEF SYNGE 
rho(ixO^S)=w(ixO^S,d_)/w(ixO^S,lfac_)
E(ixO^S)=(h(ixO^S)+dsqrt(h(ixO^S)**2.0d0+(eqpar(gamma_)**2.0d0-one)&
          *rho(ixO^S)**2.0d0))/(eqpar(gamma_)+one)
P(ixO^S) = (eqpar(gamma_)-one)/2.0d0&
	 * (E(ixO^S)-rho(ixO^S)**2.0d0/E(ixO^S))

invcsound2(ixO^S)=2.0d0*h(ixO^S)/(P(ixO^S)&
             *((eqpar(gamma_)+one)&
             +(eqpar(gamma_)-one)*(rho(ixO^S)/E(ixO^S))**2.0d0))
}

{#IFDEF GAMMA
invcsound2(ixO^S) =  (w(ixO^S,xi_)/(w(ixO^S,xi_) - w(ixO^S,lfac_)*w(ixO^S,d_))) &
               /(eqpar(gamma_)-one)
}

{#IFDEF ISO
invcsound2(ixO^S) = h(ixO^S) / (eqpar(gamma_) * eqpar(adiab_) * ( w(ixO^S,d_)/w(ixO^S,lfac_) )**eqpar(gamma_))
}

end subroutine getinvcsound2
!=============================================================================
subroutine getcsound2(w,ixI^L,ixO^L,varconserve,csound2)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!===============================================!

include 'amrvacdef.f'

integer,intent(in)                                   :: ixI^L,ixO^L
double precision, dimension(ixI^S,1:nw), intent(in)  :: w
logical, intent(in)                                  :: varconserve
double precision, dimension(ixG^T),intent(out)       :: csound2
double precision, dimension(ixG^T)                   :: h, P
{#IFDEF SYNGE
double precision, dimension(ixG^T)                   :: rho,E }
!--------------------------------------!

h(ixO^S)=w(ixO^S,xi_)/(w(ixO^S,lfac_)**2.0d0)
{#IFDEF SYNGE 
if (varconserve) then
   rho(ixO^S)=w(ixO^S,d_)/w(ixO^S,lfac_)
   E(ixO^S)=(h(ixO^S)+dsqrt(h(ixO^S)**2.0d0+(eqpar(gamma_)**2.0d0-one)&
             *rho(ixO^S)**2.0d0))/(eqpar(gamma_)+one)
   P(ixO^S) = (eqpar(gamma_)-one)/2.0d0&
              * (E(ixO^S)-rho(ixO^S)**2.0d0/E(ixO^S))

   csound2(ixO^S)=(P(ixO^S)&
                   *((eqpar(gamma_)+one)&
                   +(eqpar(gamma_)-one)*(rho(ixO^S)/E(ixO^S))**2.0d0))&
                   /(2.0d0*h(ixO^S))
else 
   E(ixO^S)=(h(ixO^S)+dsqrt(h(ixO^S)**2.0d0+(eqpar(gamma_)**2.0d0-one)&
             *w(ixO^S,rho_)**2.0d0))/(eqpar(gamma_)+one)
   P(ixO^S) = (eqpar(gamma_)-one)/2.0d0&
              * (E(ixO^S)-w(ixO^S,rho_)**2.0d0/E(ixO^S))

   csound2(ixO^S)=(P(ixO^S)&
                   *((eqpar(gamma_)+one)&
                   +(eqpar(gamma_)-one)*(w(ixO^S,rho_)/E(ixO^S))**2.0d0))&
                   /(2.0*h(ixO^S))
end if
}
{#IFDEF GAMMA
call Pressure(w,ixI^L,ixO^L,varconserve,P)
csound2(ixO^S) = (eqpar(gamma_) * P(ixO^S)) / h(ixO^S)
}
{#IFDEF ISO
call Pressure(w,ixI^L,ixO^L,varconserve,P)
csound2(ixO^S) = (eqpar(gamma_) * P(ixO^S)) / h(ixO^S)
}

end subroutine getcsound2
!=============================================================================
subroutine FuncPressure(xicurrent,lfac,d,ssqr,tau,dlfacdxi,p,dpdxi)

! compute pointwise value for pressure p and dpdxi

include 'amrvacdef.f'
  
double precision, intent(in)         :: xicurrent,lfac,d,ssqr,tau,dlfacdxi
double precision, intent(out)        :: p,dpdxi
! .. local ..
double precision                     :: rho,h
{#IFNDEF SYNGE
double precision                     :: dpdchi, dpdrho }
{#IFDEF SYNGE
double precision:: E
double precision:: dhdxi,dEdxi  }

!_______________________________________________________________!


h=xicurrent/(lfac**2.0d0)
rho=d/lfac

{#IFDEF SYNGE
E =(h+dsqrt(h**2.0d0+(eqpar(gamma_)**2.0d0-one)*rho**2.0d0))&
	      /(eqpar(gamma_)+one)

! output pressure
p=(eqpar(gamma_)-one)/2.0d0* (E-rho**2.0d0/E)

dhdxi = one/(lfac**2.0d0)-2.0d0*xicurrent/(lfac**2.0d0)*dlfacdxi/lfac

dEdxi=(dhdxi+(h*dhdxi-(eqpar(gamma_)**2.0d0-one)*rho**2.0d0*dlfacdxi/lfac)&
	/dsqrt(h**2.0d0+(eqpar(gamma_)**2.0d0-one)*rho**2.0d0))&
	/(eqpar(gamma_)+one) 

! output pressure derivative to xi
dpdxi=half*(eqpar(gamma_)-one)*(2.0d0*(rho**2.0d0/E)*dlfacdxi/lfac+&
	(one+(rho/E)**2.0)*dEdxi)

}

{#IFDEF GAMMA
p = (h - rho)/govergminone
dpdchi = one/govergminone
dpdrho = zero
}

{#IFDEF ISO
! This should work in principle. However we are using 
! the inversion based on momentum equation (con2primIso.t)
! So this will never occur.  
p = eqpar(adiab_) * rho**eqpar(gamma_)
dpdchi = zero
dpdrho = eqpar(adiab_) * eqpar(gamma_) * rho**(eqpar(gamma_)-one)
}

{#IFNDEF SYNGE
! This is quite general.  Whenever you can give expressions 
! for dpdchi and dpdrho from your EOS, this will do the inversion for your EOS.
dpdxi = dpdchi * (one/lfac**2 + ((d*lfac-2.0d0*xicurrent)/lfac**3.0d0) * dlfacdxi) &
      - dpdrho * d/lfac**2 * dlfacdxi
}
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
LsmallE=Lsmallp/(eqpar(gamma_)-one)+&
	dsqrt((Lsmallp/(eqpar(gamma_)-one))**2.0d0+Lsmallrho**2.0d0)

smalltau=half*((eqpar(gamma_)+one)*LsmallE-&
    (eqpar(gamma_)-one)*Lsmallrho**2.0d0/LsmallE)-Lsmallp-Lsmallrho
! may need to replace by smallxi above

smallxi=half*((eqpar(gamma_)+one)*LsmallE-&
   (eqpar(gamma_)-one)*Lsmallrho**2.0d0/LsmallE)

end Subroutine smallvaluesEOS
!=============================================================================!
subroutine xinoFlow(xi,tau,d,bb)

include 'amrvacdef.f'

double precision:: xi,tau,d,bb
!_______________________________________________________________!

{#IFDEF SYNGE
xi = half*((eqpar(gamma_)+one)*(tau+d-half*bb)-(eqpar(gamma_)-one)*d**2/(tau+d-half*bb))
}
{#IFDEF GAMMA
xi = eqpar(gamma_) *(tau - half * bb)+d
}
{#IFDEF ISO
xi = d + govergminone * eqpar(adiab_) * d**eqpar(gamma_)
}
end subroutine xinoFlow
{#IFDEF ENERGY
!=============================================================================
subroutine FuncEnthalpy(pcurrent,lfac2inv,d,s^C,tau,sqrs,xicurrent,dv2d2p,h,dhdp,ierror)

include 'amrvacdef.f'
  
double precision, intent(in)             :: pcurrent,lfac2inv,d,s^C,tau,sqrs,xicurrent,dv2d2p
double precision, intent(out)            :: h,dhdp
integer, intent(inout)                   :: ierror

double precision                         :: rho
{#IFDEF SYNGE
double precision:: E_th,E
double precision:: dE_thdp,dEdp
}
!_______________________________________________________________!
rho=d*dsqrt(lfac2inv)

{#IFDEF SYNGE
!== The Classical definition of the thermal energy ==!
E_th = pcurrent/(eqpar(gamma_)-one)
E = (E_th + dsqrt(E_th**2.0d0+rho**2.0d0))
!== Enthalpy ==!
h = half *((eqpar(gamma_)+one) * E-(eqpar(gamma_)-one) * rho*(rho/E))


!=== Derivative of thermal energy ===!
dE_thdp = (one/(eqpar(gamma_)-one))

!=== Derivative of internal energy ===!
dEdp = dE_thdp * (one+E_th/dsqrt(E_th**2.0d0+rho**2.0d0))&
	+  d**2.0d0*dv2d2p/dsqrt(E_th**2.0d0+rho**2.0d0)

!====== Derivative of Enthalpy ======!
dhdp = half*((eqpar(gamma_)+one)*dEdp + &
 (eqpar(gamma_)-one)*(rho**2.0d0/E)*(-2.0d0*dv2d2p/lfac2inv+dEdp/E)) 
}

{#IFDEF GAMMA
h = rho + govergminone * pcurrent
dhdp = govergminone + d/sqrt(lfac2inv)*sqrs/xicurrent**3.0d0
}


return
end subroutine FuncEnthalpy
!=========================================================================
subroutine Bisection_Enthalpy(pcurrent,lfac2inv,d,s^C,tau,sqrs,xicurrent,h,ierror)

include 'amrvacdef.f'
  
integer:: ixI^L,ixO^L
double precision:: pcurrent,lfac2inv,d,s^C,tau,sqrs,xicurrent
double precision:: h
integer::ierror
double precision:: rho
{#IFDEF SYNGE
double precision:: E_th,E
}
!_______________________________________________________________!
rho=d*dsqrt(lfac2inv)
{#IFDEF SYNGE
E_th = (one/(eqpar(gamma_)-one) * pcurrent)
E = (E_th + dsqrt(E_th**2.0+rho**2.0d0))
!== Enthalpy ==!
h = half *((eqpar(gamma_)+one) * E-(eqpar(gamma_)-one) * rho*(rho/E))
}

{#IFDEF GAMMA
h = rho + govergminone * pcurrent
}

return
end subroutine Bisection_Enthalpy
}
!=============================================================================
subroutine Calcule_Geffect(w,ixI^L,ixO^L,varconserve,Geffect)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!   both cases assume updated auxiliary variables xi_ en lfac_
!===============================================!

include 'amrvacdef.f'

integer:: ixI^L,ixO^L
double precision,intent(in):: w(ixI^S,1:nw)
logical,intent(in)         :: varconserve
double precision,intent(out):: Geffect(ixG^T)
{#IFDEF SYNGE
double precision:: p(ixG^T)
double precision:: rho(ixG^T)
}
!_______________________________________________________________!
{#IFDEF SYNGE
call Pressure(w,ixI^L,ixO^L,varconserve,p)
if (varconserve) then
  rho(ixO^S)=w(ixO^S,d_)/w(ixO^S,lfac_)
  Geffect(ixO^S) = eqpar(gamma_)-half*(eqpar(gamma_)-one)*(one-(rho(ixO^S)/&
    (p(ixO^S)/(eqpar(gamma_)-one)+&
    dsqrt((p(ixO^S)/(eqpar(gamma_)-one))**2.0d0+rho(ixO^S)**2.0d0)))**2.0d0)
else
  Geffect(ixO^S) = eqpar(gamma_)-half*(eqpar(gamma_)-one)*(one-(w(ixO^S,rho_)/&
    (p(ixO^S)/(eqpar(gamma_)-one)+&
    dsqrt((p(ixO^S)/(eqpar(gamma_)-one))**2.0d0+w(ixO^S,rho_)**2.0d0)))**2.0d0)
end if
}
{#IFDEF GAMMA
Geffect(ixO^S) = eqpar(gamma_)
}
{#IFDEF ISO
Geffect(ixO^S) = eqpar(gamma_)
}

end subroutine Calcule_Geffect
!=============================================================================
! end module amrvacphys - srmhdeos
!##############################################################################
