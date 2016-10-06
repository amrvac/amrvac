!=============================================================================
! amrvacusr.t.nul
!=============================================================================
!INCLUDE:amrvacnul/specialini.t
!INCLUDE:amrvacnul/speciallog.t
!INCLUDE:amrvacnul/specialbound.t
!INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
!=============================================================================
! amrvacusr.t.nul
!=============================================================================

subroutine initglobaldata_usr
use constants
use mod_global_parameters

double precision    :: Me, Ee, L, Ldash
!----------------------------------------------------------------------------

!!!! Special refinement !!!!!!!!
addlevel=1
addlevel_origin=7
! addlevel is adaptively modified in printlog_special (see below).
! set the limiter on shock to minmod:
!typelimiter1(mxnest-addlevel_origin+addlevel:mxnest) = 'minmod'
! But user has to reset this value at restart!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Parameters in cgs units:
L     = 5.d38              ! erg s^-1; Spin down luminosity
Me    = 3.0d0*CONST_MSun   ! g       ; Mass of shell
Ee    = 1.0d51             ! erg     ; Supernova energy in shell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dynamical parameters: 
eqpar(gamma_) = 4.0d0/3.0d0
eqpar(b_)     = 0.03d0
eqpar(eta_)   = 45.0d0*dpi/180.
eqpar(theta0_)= 10.0d0*dpi/180.
eqpar(sig0_)  = 0.01d0
eqpar(lor_)   = 10.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Perturbations of the SNR:
eqpar(m_) = 0.0d0
eqpar(DeltaRho_) = 0.0d0
if (typeaxial .ne. 'spherical') then
   eqpar(DeltaRho_) = zero
endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Power injected into domain: 
Ldash = dpi * (8.0d0/3.0d0+4.0d0*eqpar(b_))    !   ; Spin down luminosity in code units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Perform the physical scaling: 
UNIT_VELOCITY  = CONST_c     ! cm s^-1
UNIT_LENGTH    = 1.D18       ! cm
UNIT_DENSITY   = L/Ldash * (UNIT_LENGTH/UNIT_VELOCITY)**3.0D0 * UNIT_LENGTH**(-5.0D0)  ! g cm^-3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Adopt parameters for Crab: 
eqpar(ri_)    = 1.0               !    ; initial inner shell radius in code units
eqpar(re_)    = 5.0               !    ; initial outer shell radius in code units
eqpar(ve_)    = sqrt(10.0d0/3.0d0*Ee/Me)/UNIT_VELOCITY    !    ; velocity at outer shell radius
eqpar(rhoe_)  = Me*(4.0d0*dpi/3.0d0*(eqpar(re_)**3.0d0-eqpar(ri_)**3.0d0))**(-1.0d0)    &
     /(UNIT_DENSITY * UNIT_LENGTH**3.0d0)                 !    ; density of shell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Fill the scaling arrays for output: 
normvar(0)       = UNIT_LENGTH
normvar(rho_)    = UNIT_DENSITY
{^C&normvar(v^C_)   = UNIT_VELOCITY \}
normvar(pp_)     = UNIT_VELOCITY**2.0d0 * UNIT_DENSITY
{^C&normvar(b^C_)   = sqrt(4.0d0*dpi*normvar(pp_)) \}
{#IFDEF EPSINF
normvar(epsinf_) = CONST_Peta*CONST_eV
normvar(rho0_)   = UNIT_DENSITY
normvar(rho1_)   = UNIT_DENSITY
normvar(n0_)     = UNIT_DENSITY
normvar(n_)      = UNIT_DENSITY
}
normvar(xi_)     = CONST_c**2.0d0
{^FL& normvar(tr^FL_)      = one\}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{#IFDEF EPSINF
! Electron density in shell:
eqpar(rho1e_)    = 1.d-5
! Floor value for epsinf evolution:
eqpar(epsfloor_) = 1.d-9 !PeV
! Floor value for advected rho0 (variable rho0_):
eqpar(rho0floor_)= eqpar(rho1e_)*1.d3
}
! threshold for Lorentz factor:
eqpar(flfac_)=0.9d0
eqpar(flfacNeighbor_)=0.7d0


if (mype == 0) then
   print*, '---------------------------------------------------------'
   print*, 'UNIT_LENGTH: ',     UNIT_LENGTH, 'cm'
   print*, 'UNIT_VELOCITY: ',   UNIT_VELOCITY, 'cm s^-1'
   print*, 'UNIT_DENSITY: ',    UNIT_DENSITY, 'g cm^-3'
   print*, '---------------------------------------------------------'
   print*, 'Density of ejecta: ', eqpar(rhoe_), 'that is', eqpar(rhoe_)*UNIT_DENSITY, 'g/cm^3'
   print*, 'Velocity of ejecta: ', eqpar(ve_), 'that is', eqpar(ve_)*UNIT_VELOCITY, 'cm/s'
   print*, 'Convergence time:' ,  eqpar(re_)*UNIT_LENGTH/(UNIT_VELOCITY*eqpar(ve_))/CONST_years, 'years'
   print*, 'One timeunit corresponds to ', UNIT_LENGTH/UNIT_VELOCITY/CONST_years, 'years'
   print*, 'Perturbarion of SNR, m: ', eqpar(m_)
   print*, 'Perturbarion of SNR, DeltaRho: ', eqpar(DeltaRho_)
   print*, '---------------------------------------------------------'
end if

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision  :: ftot(ixG^S), &
     fm(ixG^S), fk(ixG^S), sigma(ixG^S), delta(ixG^S), &
     rcyl(ixG^S), rsphere(ixG^S), costheta(ixG^S), sintheta(ixG^S), theta(ixG^S),&
     vR(ixG^S), rhoh(ixG^S), bphi(ixG^S)
{^IFTHREED double precision  :: sinphi(ixG^S), cosphi(ixG^S), phi(ixG^S)}
logical :: patchw(ixG^T)
!----------------------------------------------------------------------------
patchw(ixG^S)=.false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do some geometry: 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{^IFTHREED
if (typeaxial .ne. 'spherical') then
   rsphere(ixG^S)  = dsqrt({^D&dabs(x(ixG^S,^D))**2.0d0+})
   rcyl(ixG^S)     = dsqrt(dabs(x(ixG^S,1))**2.0d0 + dabs(x(ixG^S,2))**2.0d0)
   costheta(ixG^S) = x(ixG^S,3)/rsphere(ixG^S)
   cosphi(ixG^S)   = x(ixG^S,1)/rcyl(ixG^S)
   sinphi(ixG^S)   = x(ixG^S,2)/rcyl(ixG^S)
   theta(ixG^S)    = acos(costheta(ixG^S))  
   phi(ixG^S)      = acos(cosphi(ixG^S))  
else 
   rsphere(ixG^S)  = dabs(x(ixG^S,1))
   theta(ixG^S)    = x(ixG^S,2)
   phi(ixG^S)      = x(ixG^S,3)
   sintheta(ixG^S) = sin(x(ixG^S,2))
   costheta(ixG^S) = cos(x(ixG^S,2))
   cosphi(ixG^S)   = cos(x(ixG^S,3))
   sinphi(ixG^S)   = sin(x(ixG^S,3))
   rcyl(ixG^S)     = rsphere(ixG^S) * sintheta(ixG^S)
end if
}
{^IFTWOD
if (typeaxial .ne. 'spherical') then
   rsphere(ixG^S)  = dsqrt({^D&dabs(x(ixG^S,^D))**2.0d0+})
   rcyl(ixG^S)     = dsqrt(dabs(x(ixG^S,1))**2.0d0)
   costheta(ixG^S) = x(ixG^S,2)/rsphere(ixG^S)
   theta(ixG^S)    = acos(costheta(ixG^S))  
else 
   rsphere(ixG^S)  = dabs(x(ixG^S,1))
   theta(ixG^S)    = x(ixG^S,2)
   sintheta(ixG^S) = sin(x(ixG^S,2))
   costheta(ixG^S) = cos(x(ixG^S,2))
   rcyl(ixG^S)     = rsphere(ixG^S) * sintheta(ixG^S)
end if
}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get sigma:
call get_sigma(ixG^L,ixG^L,w,x,sigma)
! Get delta:
call get_delta(ixG^L,ixG^L,w,x,delta)

! Energy fluxes:
ftot(ixG^S) = 1./dabs(rsphere(ixG^S))**2.0d0*&
     (dabs(dsin(theta(ixG^S)))**2.d0+eqpar(b_))
fm(ixG^S) = ftot(ixG^S)/(1.+1./(sigma(ixG^S)))
fk(ixG^S) = ftot(ixG^S)/(1.+sigma(ixG^S))

! Set Lorentz factor and radial velocity:
w(ixG^S,lfac_) = eqpar(lor_)
vR(ixG^S)     = dsqrt(one-one/w(ixG^S,lfac_)**2.0d0)

! Obtain density:
w(ixG^S,rho_) = fk(ixG^S)/(eqpar(lor_)**2.0d0*vR(ixG^S))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fill the wind: 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{^IFTHREED
if (typeaxial .ne. 'spherical') then
   bphi(ixG^S)       = dsqrt(fm(ixG^S)/vR(ixG^S))*x(ixG^S,3)/dabs(x(ixG^S,3))
   w(ixG^S,b1_)      = - bphi(ixG^S) * sinphi(ixG^S)
   w(ixG^S,b2_)      = bphi(ixG^S) * cosphi(ixG^S)
   {^D&w(ixG^S,v^D_) = vR(ixG^S) * x(ixG^S,^D)/rsphere(ixG^S)\}
else
   bphi(ixG^S)  = dsqrt(fm(ixG^S)/vR(ixG^S)) * costheta(ixG^S)/dabs(costheta(ixG^S))
   w(ixG^S,b3_) = bphi(ixG^S)
   w(ixG^S,v1_) = vR(ixG^S)
end if
}
{^IFTWOD
if (typeaxial .ne. 'spherical') then
   bphi(ixG^S)       = dsqrt(fm(ixG^S)/vR(ixG^S))*x(ixG^S,2)/dabs(x(ixG^S,2))
   w(ixG^S,bphi_)    = bphi(ixG^S)
   {^D&w(ixG^S,v^D_) = vR(ixG^S) * x(ixG^S,^D)/rsphere(ixG^S)\}
else
   bphi(ixG^S)  = dsqrt(fm(ixG^S)/vR(ixG^S)) * costheta(ixG^S)/dabs(costheta(ixG^S))
   w(ixG^S,b3_) = bphi(ixG^S)
   w(ixG^S,v1_) = vR(ixG^S)
end if
}

w(ixG^S,pp_)   = w(ixG^S,rho_)/100.d0
rhoh(ixG^S) = w(ixG^S,rho_) + w(ixG^S,pp_)*eqpar(gamma_)/(eqpar(gamma_)-one)
w(ixG^S,xi_) = w(ixG^S,lfac_)**2.0d0*rhoh(ixG^S)

w(ixG^S,flg_)   = one
{#IFDEF EPSINF
w(ixG^S,epsinf_)= one         ! PeV
w(ixG^S,rho0_)  = 1./dabs(rsphere(ixG^S))**2.0d0
w(ixG^S,rho1_)  = w(ixG^S,rho0_)
w(ixG^S,n_)     = (1.-delta(ixG^S))*w(ixG^S,rho0_) + delta(ixG^S) * 0.001* w(ixG^S,rho0_)
w(ixG^S,n0_)    = w(ixG^S,n_)
}
w(ixG^S,tr1_)   = one
{#IFDEF GLM
w(ixG^S,psi_)   = zero
}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fill external medium: 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (.not. time_advance) then 
   where(rsphere(ixG^S)>eqpar(ri_))
      {^C&w(ixG^S,b^C_)=zero \}
      {^C&w(ixG^S,v^C_)=zero \}
{^IFTWOD
      w(ixG^S,rho_)    = eqpar(rhoe_) &
           + eqpar(DeltaRho_)*eqpar(rhoe_)*sin(2.0d0*dpi*theta(ixG^S) * eqpar(m_)/(xprobmax2-xprobmin2))
}{^IFTHREED
      w(ixG^S,rho_)    = eqpar(rhoe_) &
           + eqpar(DeltaRho_)*eqpar(rhoe_)*sin(2.0d0*dpi*theta(ixG^S) * eqpar(m_)/(xprobmax2-xprobmin2)) &
           * sin(2.0d0*dpi*phi(ixG^S) * eqpar(m_)/(xprobmax3-xprobmin3))
}
      w(ixG^S,pp_)     = eqpar(rhoe_)*tlow
      w(ixG^S,flg_)    = zero
{#IFDEF EPSINF
      w(ixG^S,rho0_)   = eqpar(rho0floor_)
      w(ixG^S,rho1_)   = eqpar(rho1e_)
      w(ixG^S,n0_)     = eqpar(rho0floor_)
      w(ixG^S,n_)      = eqpar(rho1e_)
      w(ixG^S,epsinf_) = eqpar(epsfloor_)
}
      w(ixG^S,tr1_)    = zero
   end where

   if (typeaxial .ne. 'spherical') then
      where(rsphere(ixG^S)>eqpar(ri_))
         {^D&w(ixG^S,v^D_)= eqpar(ve_) * x(ixG^S,^D)/eqpar(re_)\}
      end where
   else
      where(rsphere(ixG^S)>eqpar(ri_))
         w(ixG^S,v1_)     = eqpar(ve_) * rsphere(ixG^S)/eqpar(re_)
      end where
   end if ! spherical
end if ! time_advance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (useprimitiveRel) then
   w(ixG^S,lfac_)=one/dsqrt(one-({^C&w(ixG^S,v^C_)**2.0d0+}))
   {^C&w(ixG^S,u^C_)=w(ixG^S,lfac_)*w(ixG^S,v^C_);\}
endif
call conserve(ixG^L,ix^L,w,x,patchw)

end subroutine initonegrid_usr
!=============================================================================
subroutine flag_grid_usr(qt,ixI^L,ixO^L,w,x,flag)
use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L
integer, intent(inout)          :: flag
double precision, intent(in)    :: qt
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)

double precision, dimension(ixI^S) :: lfac_pw

! flag=0 : Treat as normal domain
! flag=1 : Treat as passive, but reduce by safety belt
! flag=2 : Treat as passive and dont care about safety belt

!-----------------------------------------------------------------------------
flag = 0

call get_lfac(ixI^L,ixO^L,w,x,lfac_pw)
if ( .not.any( w(ixO^S,lfac_)/lfac_pw(ixO^S) < eqpar(flfac_)) ) &
     flag=1
if ({^D& minval(x(ixI^S,^D))*maxval(x(ixI^S,^D)) < zero |.and.}) &
     flag=3

end subroutine flag_grid_usr
!=============================================================================
subroutine fixp_usr(ixI^L,ixO^L,w,x)
use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(inout)    :: w(ixI^S,1:nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)

logical :: patchp(ixI^S)
!----------------------------------------------------------------------------
patchp(ixO^S) = .true.

where(w(ixO^S,pp_)/w(ixO^S,rho_) < tlow) 
   patchp(ixO^S) = .false.
end where

if (any(.not.patchp(ixO^S))) then
   where(.not.patchp(ixO^S))
      w(ixO^S,pp_) = tlow*w(ixO^S,rho_)
   end where
   call conserven(ixI^L,ixO^L,w,patchp)
   call primitiven(ixI^L,ixO^L,w,patchp)
endif

end subroutine fixp_usr
!=============================================================================
subroutine correctaux_usr(ixI^L,ixO^L,w,x,patchierror,subname)

use mod_global_parameters

! This is for TM eos:
! E_Th = p/(eqpar(gamma_)-one) 
! E = E_Th + dsqrt(E_Th**2.0d0+rho**2.0d0)
! rhoh = half*((eqpar(gamma_)+one) * E-&
!	       (eqpar(gamma_)-one)* rho*(rho/E))
!
! This is for adiabatic eos:
!rhoh = rho + p*eqpar(gamma_)/(eqpar(gamma_)-one)

integer, intent(in)            :: ixI^L, ixO^L
integer, intent(inout)         :: patchierror(ixG^T)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixI^S,1:nw)
double precision, intent(in)   :: x(ixI^S,1:ndim)

! .. local ..
integer                        :: ix^D
double precision               :: rhoh(ixI^S), vR(ixI^S)
double precision, dimension(ixI^S) :: SdotB, sqrB, absu,rsphere
double precision, dimension(ixI^S) :: lfac_pw, delta, lfac_save
logical           :: patchw(ixG^T)
integer           :: idirmin
!----------------------------------------------------------------------------

  patchw=.true.

! Obtain pw Lorentz factor:
call get_lfac(ixI^L,ixO^L,w,x,lfac_pw)
call get_delta(ixI^L,ixO^L,w,x,delta)

if (typeaxial .ne. 'spherical') then
   rsphere(ixO^S)  = dsqrt({^D&dabs(x(ixO^S,^D))**2.0d0+})
else
   rsphere(ixO^S)  = dabs(x(ixO^S,1))
endif
! Take density and bfield from solution, set vel perp. to B:
   where (patchierror(ixO^S) /=0 .and. w(ixO^S,d_) > zero .and. &
      w(ixO^S,lfac_)/lfac_pw(ixO^S) > eqpar(flfac_) .and. &
      w(ixO^S,flg_) /= 0 )

      patchw(ixO^S)=.false.
      patchierror(ixO^S)=0
      lfac_save(ixO^S) = w(ixO^S,lfac_)
   
!      SdotB(ixO^S)= {^C&w(ixO^S,s^C_)*w(ixO^S,b^C_)+}
! Assume flow perpendicular to B like pw solution:
      SdotB(ixO^S)= zero
      sqrB(ixO^S) = {^C&w(ixO^S,b^C_)**2.0d0+}  

      w(ixO^S,rho_) = w(ixO^S,d_)/lfac_pw(ixO^S)
{#IFDEF EPSINF
      w(ixO^S,rho1_)= one/dabs(rsphere(ixO^S))**2.0d0
      w(ixO^S,rho0_)= w(ixO^S,rho1_)
      w(ixO^S,n_)     = (1.-delta(ixO^S))*w(ixO^S,rho0_) + delta(ixO^S) * 0.001* w(ixO^S,rho0_)
      w(ixO^S,n0_)    = w(ixO^S,n_)
      w(ixO^S,epsinf_)= one         ! PeV
}
      w(ixO^S,tr1_) = one
      w(ixO^S,pp_)   = w(ixO^S,rho_)/100.d0      
      rhoh(ixO^S) = w(ixO^S,rho_) + w(ixO^S,pp_)*eqpar(gamma_)/(eqpar(gamma_)-one)
! Recovering lfac and xi:
      w(ixO^S,lfac_) = lfac_pw(ixO^S)
      w(ixO^S,xi_) = w(ixO^S,lfac_)**2.0d0*rhoh(ixO^S)

! This gets the vel-dir perp. to B.

     {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*(w(ixO^S,s^C_)&
                      +SdotB(ixO^S)*w(ixO^S,b^C_)/w(ixO^S,xi_))/ &
                      (w(ixO^S,xi_)+sqrB(ixO^S))\}

     absu(ixO^S) = dsqrt({^C&w(ixO^S,u^C_)**2.0d0+})
     vR(ixO^S)   = dsqrt(one-one/lfac_pw(ixO^S)**2.0d0)

     {^C&w(ixO^S,u^C_)= w(ixO^S,u^C_)/absu(ixO^S)&
          *lfac_pw(ixO^S)*vR(ixO^S)
\}

  end where

  {do ix^D= ixO^LIM^D\}
  if (.not. patchw(ix^D)) then
     print*,mype, ': ','it:',it,', repaired getaux, rho=',w(ix^D,d_),'lfac=',lfac_save(ix^D)
     print*,mype, ': ','igrid:',saveigrid
     print*,mype, ': ','position  ', px(saveigrid)%x(ix^D, 1:ndim)
     print*,mype, ': ','called from ', trim(subname)
  end if
  {enddo^D&\}

if ( any(w(ixO^S,d_) <= zero) ) then
   print*,mype,': ', 'rhos:',w(ixO^S,d_),&
        'ps:',w(ixO^S,pp_),&
        'v1:',w(ixO^S,v1_),&
        'v2:',w(ixO^S,v2_),&
        'v3:',w(ixO^S,v3_),&
        'b1:',w(ixO^S,b1_),&
        'b2:',w(ixO^S,b2_),&
        'b3:',w(ixO^S,b3_)
!     call mpistop('negative density encountered.')
end if

call conserven(ixI^L,ixO^L,w,patchw)

end subroutine correctaux_usr
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

! you must set consistent values for integers refine/coarsen:

! refine = -1 enforce to not refine
! refine =  0 doesn't enforce anything
! refine =  1 enforce refinement

! coarsen = -1 enforce to not coarsen
! coarsen =  0 doesn't enforce anything
! coarsen =  1 enforce coarsen

use mod_global_parameters

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
! .. local ..
double precision, dimension(ixG^T) :: lfac_pw
double precision, parameter :: &
     lfacpw = 0.7, lfacext = 0.5
integer                     :: ixONE^L
!----------------------------------------------------------------------------
! Add one ghost zone (two would cause and oscialltory behaviour)
ixONE^L=ix^L^LADD1;

! Obtain pw Lorentz factor:
call get_lfac(ixG^L,ixG^L,w,x,lfac_pw)

! refine the blocks with shock to addlevel:

   if (level >= mxnest-(addlevel_origin-1)) coarsen = 1
   if (level >= mxnest-(addlevel_origin)) refine  = -1

if (level <mxnest-(addlevel_origin-addlevel)) then
   if (any(w(ixONE^S,lfac_)/lfac_pw(ixONE^S)>lfacpw) .and.&
    any(w(ixONE^S,lfac_)/lfac_pw(ixONE^S)<lfacext)) &
      refine = 1
end if
if (level <=mxnest-(addlevel_origin-addlevel)) then
   if (any(w(ixG^S,lfac_)/lfac_pw(ixG^S)>lfacpw) .and.&
    any(w(ixG^S,lfac_)/lfac_pw(ixG^S)<lfacext)) &
        coarsen = -1
end if

if (typeaxial .ne. 'spherical') then
   ! Check if at origin and refine to maximum: 
   if (({^D& minval(x(ixG^S,^D))*maxval(x(ixG^S,^D)) < zero |.and.})) then
      refine = 1
      coarsen = -1
   end if
end if

! Use tracers not to refine into SN remnant:
  if (.not.any(abs(w(ixG^S,Dtr1_)/w(ixG^S,d_))>=smalldouble)) then
     coarsen = 1
  end if
  if (.not.any(abs(w(ixONE^S,Dtr1_)/w(ixG^S,d_))>=smalldouble)) then
     refine = -1
  end if

end subroutine specialrefine_grid
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

use mod_global_parameters

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)
! .. local ..
double precision  :: ftot(ixI^S), fm(ixI^S), fk(ixI^S),sigma(ixI^S),&
     rcyl(ixI^S), &
     costheta(ixI^S), sintheta(ixI^S), vR(ixI^S), rsphere(ixI^S), theta(ixI^S), rhoh(ixI^S), bphi(ixI^S),&
        lfac_pw(ixI^S),delta(ixI^S)
{^IFTHREED double precision  :: sinphi(ixI^S), cosphi(ixI^S)}
logical :: patchw(ixI^S)
integer :: idirmin, i^D, iadd^D
integer :: isafety
!-----------------------------------------------------------------------------

patchw = .true.
isafety = 2
if (level < mxnest) isafety = 1

   where (w(ixO^S,flg_) < two )
      w(ixO^S,flg_) =   zero
   end where

if (any(w(ixO^S,lfac_) > eqpar(flfac_)*eqpar(lor_))) then


! Do some geometry: 

{^IFTHREED
if (typeaxial .ne. 'spherical') then
   rsphere(ixI^S)  = dsqrt({^D&dabs(x(ixI^S,^D))**2.0d0+})
   rcyl(ixI^S)     = dsqrt(dabs(x(ixI^S,1))**2.0d0 + dabs(x(ixI^S,2))**2.0d0)
   costheta(ixI^S) = x(ixI^S,3)/rsphere(ixI^S)
   cosphi(ixI^S)   = x(ixI^S,1)/rcyl(ixI^S)
   sinphi(ixI^S)   = x(ixI^S,2)/rcyl(ixI^S)
   theta(ixI^S)    = acos(costheta(ixI^S))  
else 
   rsphere(ixI^S)  = dabs(x(ixI^S,1))
   theta(ixI^S)    = x(ixI^S,2)
   sintheta(ixI^S) = sin(x(ixI^S,2))
   costheta(ixI^S) = cos(x(ixI^S,2))
   cosphi(ixI^S)   = cos(x(ixI^S,3))
   sinphi(ixI^S)   = sin(x(ixI^S,3))
   rcyl(ixI^S)     = rsphere(ixI^S) * sintheta(ixI^S)
end if
}
{^IFTWOD
if (typeaxial .ne. 'spherical') then
   rsphere(ixI^S)  = dsqrt({^D&dabs(x(ixI^S,^D))**2.0d0+})
   rcyl(ixI^S)     = dsqrt(dabs(x(ixI^S,1))**2.0d0)
   costheta(ixI^S) = x(ixI^S,2)/rsphere(ixI^S)
   theta(ixI^S)    = acos(costheta(ixI^S))  
else 
   rsphere(ixI^S)  = dabs(x(ixI^S,1))
   theta(ixI^S)    = x(ixI^S,2)
   sintheta(ixI^S) = sin(x(ixI^S,2))
   costheta(ixI^S) = cos(x(ixI^S,2))
   rcyl(ixI^S)     = rsphere(ixI^S) * sintheta(ixI^S)
end if
}

! Get sigma:
    call get_sigma(ixI^L,ixO^L,w,x,sigma)
    call get_delta(ixI^L,ixO^L,w,x,delta)

! Get Lorentz factor:
  lfac_pw(ixI^S) = eqpar(lor_)
 
! Construct array range to be processed:
{do i^D=ixOmin^D,ixOmax^D \}
if (w(i^D,lfac_)/lfac_pw(i^D) > eqpar(flfac_))&
    patchw(i^D) = .false.
{do iadd^D=-isafety,isafety \}
if (w(i^D+iadd^D,lfac_)/lfac_pw(i^D)<=eqpar(flfacNeighbor_)) &
   patchw(i^D) = .true.
{end do\}
{end do\}

where (.not.patchw(ixO^S) )
! Energy fluxes:
   ftot(ixO^S) = 1./dabs(rsphere(ixO^S))**2.0d0*&
        (dabs(dsin(theta(ixO^S)))**2.d0+eqpar(b_))
   fm(ixO^S) = ftot(ixO^S)/(1.+1./(sigma(ixO^S)))
   fk(ixO^S) = ftot(ixO^S)/(1.+sigma(ixO^S))    

! Obtain density from energy flux: 
   vR(ixO^S)     = dsqrt(one-one/eqpar(lor_)**2.0d0)
   w(ixO^S,lfac_) = lfac_pw(ixO^S)
! set scalars:
   w(ixO^S,rho_) = fk(ixO^S)/(eqpar(lor_)**2.0d0*vR(ixO^S))
   w(ixO^S,pp_)   = w(ixO^S,rho_)/100.d0
   rhoh(ixO^S) = w(ixO^S,rho_) + w(ixO^S,pp_)*eqpar(gamma_)/(eqpar(gamma_)-one)
   w(ixO^S,xi_) = w(ixO^S,lfac_)**2.0d0*rhoh(ixO^S)

{#IFDEF EPSINF
  w(ixO^S,rho1_)= one/dabs(rsphere(ixO^S))**2.0d0
   w(ixO^S,rho0_)= w(ixO^S,rho1_)
   w(ixO^S,n_)     = (1.-delta(ixO^S))*w(ixO^S,rho0_) + delta(ixO^S) * 0.001* w(ixO^S,rho0_)
   w(ixO^S,n0_)    = w(ixO^S,n_)
   w(ixO^S,epsinf_)= one         ! PeV
}
   w(ixO^S,tr1_)   = one

{#IFDEF GLM
   w(ixO^S,psi_) = 0.0d0
}
end where

if (typeaxial .ne. 'spherical') then
   where (.not.patchw(ixO^S) )
      ! set vectors:
      {^IFTHREED
      bphi(ixO^S) = dsqrt(fm(ixO^S)/vR(ixO^S))*x(ixO^S,3)/dabs(x(ixO^S,3))
      w(ixO^S,b1_) = - bphi(ixO^S) * sinphi(ixO^S)
      w(ixO^S,b2_) = bphi(ixO^S) * cosphi(ixO^S)
      w(ixO^S,b3_) = zero
      }
      {^IFTWOD
      bphi(ixO^S) = dsqrt(fm(ixO^S)/vR(ixO^S))*x(ixO^S,2)/dabs(x(ixO^S,2))
      w(ixO^S,bphi_) = bphi(ixO^S)
      }
      {^D&w(ixO^S,v^D_) = vR(ixO^S) * x(ixO^S,^D)/rsphere(ixO^S)*w(ixO^S,lfac_)\}
   end where
else
   where (.not.patchw(ixO^S) )
      ! set vectors:
      {^D&w(ixO^S,b^D_) = zero\}
      bphi(ixO^S)  = dsqrt(fm(ixO^S)/vR(ixO^S)) * costheta(ixO^S)/dabs(costheta(ixO^S))
      w(ixO^S,b3_) = bphi(ixO^S)
      {^D&w(ixO^S,v^D_) = zero\}
      w(ixO^S,v1_) = vR(ixO^S)*w(ixO^S,lfac_)
   end where
end if ! spherical

where (.not.patchw(ixO^S) .and. w(ixO^S,flg_) < two )
   w(ixO^S,flg_) = one
end where

call conserven(ixI^L,ixO^L,w,patchw)

end if

! Reset pressure in the shell, needed to reduce errors occuring at low temperatures:
if (any(abs(w(ixO^S,Dtr1_)/(w(ixO^S,d_)))<=smalldouble) .or.&
     any(abs(w(ixO^S,Dtr1_)/(w(ixO^S,d_)))<=0.5d0).and.level<mxnest-addlevel_origin) then 
patchw=.true.
   
where (abs(w(ixO^S,Dtr1_)/(w(ixO^S,d_)))<=smalldouble .or.&
     abs(w(ixO^S,Dtr1_)/(w(ixO^S,d_)))<=0.5d0.and.level<mxnest-addlevel_origin)
   patchw(ixO^S) = .false.
end where

call primitiven(ixI^L,ixO^L,w,patchw)

where (.not.patchw(ixO^S))
   w(ixO^S,pp_)   = w(ixO^S,rho_)*tlow
end where
call conserven(ixI^L,ixO^L,w,patchw)

end if

{#IFDEF EPSINF
! reset rho1 to floor:
where (w(ixO^S,Drho1_)<w(ixO^S,lfac_)*eqpar(rho1e_))
   w(ixO^S,Drho1_)=w(ixO^S,lfac_)*eqpar(rho1e_)
end where
! reset rho0 to floor:
where (w(ixO^S,Drho0_)<w(ixO^S,Drho1_)*eqpar(rho0floor_))
   w(ixO^S,Drho0_)=w(ixO^S,Drho1_)*eqpar(rho0floor_)
end where
! reset n to floor:
where (w(ixO^S,Dn_)<w(ixO^S,lfac_)*eqpar(rho1e_))
   w(ixO^S,Dn_)=w(ixO^S,lfac_)*eqpar(rho1e_)
end where
! reset n0 to floor:
where (w(ixO^S,Dn0_)<w(ixO^S,Dn_)*eqpar(rho0floor_))
   w(ixO^S,Dn0_)=w(ixO^S,Dn_)*eqpar(rho0floor_)
end where
! reset eps_inf to floor: 
where (w(ixO^S,Depsinf_)<w(ixO^S,lfac_)**(1.0d0/3.0d0)*eqpar(epsfloor_)*w(ixO^S,rho1_)**(2.0d0/3.0d0))
   w(ixO^S,Depsinf_) = w(ixO^S,lfac_)**(1.0d0/3.0d0)*eqpar(epsfloor_)*w(ixO^S,rho1_)**(2.0d0/3.0d0)
end where
! try to catch NaNs (or infs just to be save):
where ((w(ixO^S,Depsinf_) .ne. w(ixO^S,Depsinf_)) .or.(w(ixO^S,Depsinf_)*zero .ne. zero))
   w(ixO^S,Depsinf_) = eqpar(epsfloor_)*eqpar(rho1e_)**(2.0d0/3.0d0)
end where
}
end subroutine process_grid_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
! .. local ..
{#IFDEF EPSINF
logical          :: patchw(ixG^T)
double precision :: c2tilde, epstocode
double precision :: bdash2(ixG^T), e2(ixG^T), b2(ixG^T), epsinf(ixG^T)
!-----------------------------------------------------------------------------
patchw = .true.
patchw(ixO^S) = .false.

! Primitive call because we need velocities, Lorentz factor and primitive epsinf:
call primitive(ixI^L,ixO^L,w,x)
! Check first if we need to cool at all:
if (.not.any(w(ixO^S,epsinf_)>eqpar(epsfloor_))) then
   call conserven(ixI^L,ixO^L,w,patchw)
   return
end if

! Synchrotron constant c2tilde = 2/3*c2
c2tilde = 1.578d-3     ! g^-2 cm^-1 s^3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Bring to code units:
c2tilde   = c2tilde   / (UNIT_DENSITY**(-2.0D0)*UNIT_LENGTH**(-4.0D0)*UNIT_VELOCITY**(-3.0D0)) 
epstocode = normvar(epsinf_)  / (UNIT_DENSITY*UNIT_LENGTH**3.0D0*UNIT_VELOCITY**2.0D0) * 4.0d0*dpi


! Obtain co-moving magnetic field:
! using the invariant Bdash^2 = B^2-E^2
b2(ixO^S)     = {^C&w(ixO^S,b^C_)**2.0d0|+}
e2(ixO^S) = (w(ixO^S,b2_)*w(ixO^S,u3_)/w(ixO^S,lfac_) - w(ixO^S,b3_)*w(ixO^S,u2_)/w(ixO^S,lfac_))**2.0d0 &
     + (w(ixO^S,b3_)*w(ixO^S,u1_)/w(ixO^S,lfac_) - w(ixO^S,b1_)*w(ixO^S,u3_)/w(ixO^S,lfac_))**2.0d0 &
     + (w(ixO^S,b1_)*w(ixO^S,u2_)/w(ixO^S,lfac_) - w(ixO^S,b2_)*w(ixO^S,u1_)/w(ixO^S,lfac_))**2.0d0
bdash2(ixO^S) = b2(ixO^S) - e2(ixO^S)

! Store the primitive cutoff energy:
epsinf(ixO^S) = w(ixO^S,epsinf_)

call conserve(ixI^L,ixO^L,w,x,patchw)

where (epsinf(ixO^S)>eqpar(epsfloor_))
   w(ixO^S,Depsinf_) = w(ixO^S,Depsinf_) - qdt*c2tilde &
        * abs(epsinf(ixO^S))**2.0d0*bdash2(ixO^S)*abs(w(ixO^S,rho1_)/w(ixO^S,lfac_))**(2.0d0/3.0d0) &
        * epstocode
end where

! reset to floor: 
where (w(ixO^S,Depsinf_)<w(ixO^S,lfac_)**(1.0d0/3.0d0)*eqpar(epsfloor_)*abs(w(ixO^S,rho1_))**(2.0d0/3.0d0))
   w(ixO^S,Depsinf_) = w(ixO^S,lfac_)**(1.0d0/3.0d0)*eqpar(epsfloor_)*abs(w(ixO^S,rho1_))**(2.0d0/3.0d0)
end where
}
end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!.. local ..
{#IFDEF EPSINF
logical          :: patchw(ixG^T)
double precision :: epstocode, epsinf(ixG^T)
double precision :: bdash2(ixG^T), e2(ixG^T), b2(ixG^T)
double precision :: c2tilde
!
double precision, parameter :: dtfac=0.3d0
!-----------------------------------------------------------------------------
! Check first if we need to cool further:
if (.not.any(w(ix^S,Depsinf_)>w(ix^S,lfac_)**(1.0d0/3.0d0)*eqpar(epsfloor_)*w(ix^S,Drho1_)**(2.0d0/3.0d0))) then 
   dtnew=bigdouble
   return
end if

patchw = .true.
patchw(ix^S) = .false.

! Synchrotron constant c2tilde = 2/3*c2
c2tilde = 1.578d-3     ! g^-2 cm^-1 s^3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Bring to code units:
c2tilde   = c2tilde   / (UNIT_DENSITY**(-2.0D0)*UNIT_LENGTH**(-4.0D0)*UNIT_VELOCITY**(-3.0D0)) 
epstocode = normvar(epsinf_)  / (UNIT_DENSITY*UNIT_LENGTH**3.0D0*UNIT_VELOCITY**2.0D0) * 4.0d0*dpi

call primitiven(ixG^L,ix^L,w,patchw)

! Obtain co-moving magnetic field:
! using the invariant Bdash^2 = B^2-E^2
b2(ix^S)     = {^C&w(ix^S,b^C_)**2.0d0|+}
e2(ix^S) = (w(ix^S,b2_)*w(ix^S,u3_)/w(ix^S,lfac_) - w(ix^S,b3_)*w(ix^S,u2_)/w(ix^S,lfac_))**2.0d0 &
     + (w(ix^S,b3_)*w(ix^S,u1_)/w(ix^S,lfac_) - w(ix^S,b1_)*w(ix^S,u3_)/w(ix^S,lfac_))**2.0d0 &
     + (w(ix^S,b1_)*w(ix^S,u2_)/w(ix^S,lfac_) - w(ix^S,b2_)*w(ix^S,u1_)/w(ix^S,lfac_))**2.0d0
bdash2(ix^S) = b2(ix^S) - e2(ix^S)

! Store the primitive cutoff energy:
epsinf(ix^S) = w(ix^S,epsinf_)
call conserven(ixG^L,ix^L,w,patchw)


dtnew=dtfac/c2tilde/epstocode * min(abs(minval(w(ix^S,lfac_)/epsinf(ix^S)/bdash2(ix^S))),bigdouble)
}
end subroutine getdt_special
!=============================================================================
subroutine specialbound_usr(qt,ixI^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw, iB
double precision, intent(in) :: qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)
! .. local ..
double precision  :: ftot(ixG^T), fm(ixG^T), fk(ixG^T), sigma(ixG^T), &
     vR(ixG^T), rsphere(ixG^T), theta(ixG^T), rhoh(ixG^T), bphi(ixG^T),&
     lfac_pw(ixG^T), delta(ixG^T)
logical :: patchw(ixG^T)
!----------------------------------------------------------------------------
patchw(ixG^T)         = .true.

if (iB == 1) then
   if (typeaxial .eq. 'spherical') then
      patchw(ixO^S)   = .false.
      ! Do some geometry:
      rsphere(ixO^S)  = abs(x(ixO^S,1))
      theta(ixO^S)    = x(ixO^S,2)

      if (any(rsphere(ixO^S).lt.zero)) call mpistop('negative radius in injection boundary')

      ! Get sigma:
      call get_sigma(ixI^L,ixO^L,w,x,sigma)
      call get_delta(ixI^L,ixO^L,w,x,delta)
      ! Get Lorentz factor:
      call get_lfac(ixI^L,ixO^L,w,x,lfac_pw)

      ! Energy fluxes:
      ftot(ixO^S)     = 1./rsphere(ixO^S)**2.0d0*&
           (dabs(dsin(theta(ixO^S)))**2.d0+eqpar(b_))
      fm(ixO^S)       = ftot(ixO^S)/(1.+1./(sigma(ixO^S)))
      fk(ixO^S)       = ftot(ixO^S)/(1.+sigma(ixO^S))    

      vR(ixO^S)       = dsqrt(one-one/lfac_pw(ixO^S)**2.0d0)
      w(ixO^S,lfac_)  = lfac_pw(ixO^S)

      ! set scalars:
      w(ixO^S,rho_)   = fk(ixO^S)/(eqpar(lor_)**2.0d0*vR(ixO^S))
      w(ixO^S,pp_)    = w(ixO^S,rho_)/100.d0
      rhoh(ixO^S)     = w(ixO^S,rho_) + w(ixO^S,pp_)*eqpar(gamma_)/(eqpar(gamma_)-one)
      w(ixO^S,xi_)    = w(ixO^S,lfac_)**2.0d0*rhoh(ixO^S)

{#IFDEF EPSINF
      w(ixO^S,rho1_)  = one/rsphere(ixO^S)**2.0d0
      w(ixO^S,rho0_)  = w(ixO^S,rho1_)
      w(ixO^S,n_)     = (1.-delta(ixO^S))*w(ixO^S,rho0_) + delta(ixO^S) * 0.001* w(ixO^S,rho0_)
      w(ixO^S,n0_)    = w(ixO^S,n_)
      w(ixO^S,epsinf_)= one         ! PeV
}
      w(ixO^S,tr1_)   = one
      w(ixO^S,flg_)   = one

{#IFDEF GLM
      w(ixO^S,psi_)   = 0.0d0
}

      ! set vectors:
      {^C&w(ixO^S,b^C_) = zero\}
      bphi(ixO^S)     = dsqrt(fm(ixO^S)/vR(ixO^S)) * cos(theta(ixO^S))/dabs(cos(theta(ixO^S)))
      w(ixO^S,b3_)    = bphi(ixO^S)
      {^C&w(ixO^S,v^C_) = zero\}
      w(ixO^S,v1_)    = vR(ixO^S)*w(ixO^S,lfac_)

      call conserven(ixI^L,ixO^L,w,patchw)
   end if !spherical
end if !iB

end subroutine specialbound_usr
!=============================================================================
subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
! corresponding normalization values (default value 1)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)
{^IFTHREED
! .. local ..
double precision:: divb(ixG^T)
!-----------------------------------------------------------------------------

!call mpistop("special output file undefined")
call getdivb(w,ixI^L,ixO^L,divb)
 w(ixO^S,nw+1)=divb(ixO^S)
}
end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

!call mpistop("special varnames and primnames undefined")
{^IFTHREED
 primnames= TRIM(primnames)//' '//'divb'
 wnames=TRIM(wnames)//' '//'divb'
}
end subroutine specialvarnames_output
!=============================================================================
subroutine printlog_special

! printlog: calculates volume averaged mean values 
use mod_forest,only:nleafs,nleafs_active,nleafs_level
use mod_timing
use mod_global_parameters

logical :: fileopen
integer :: iigrid, igrid, level, iw, i
double precision :: wmean(1:nw), volume(1:nlevelshi), volprob, voltotal
double precision :: dvolume(ixG^T), volumeflat(1:nlevelshi)
integer :: numlevels, nx^D, nc, ncells, dit
double precision :: dtTimeLast, now, cellupdatesPerSecond, activeBlocksPerCore, wctPerCodeTime, timeToFinish
integer, dimension(1:nlevelshi) :: isum_send
double precision, dimension(1:nw+1+nlevelshi) :: dsum_send, dsum_recv
character(len=80) :: filename
character(len=1024) :: line
logical, save :: opened=.false.
integer :: amode, status(MPI_STATUS_SIZE)
{^IFTHREED integer, parameter :: nmin=100, nmax=10000}
{^IFTWOD integer, parameter :: nmin=100, nmax=600}
!-----------------------------------------------------------------------------
volume(1:mxnest)=zero
volumeflat(1:mxnest)=zero
wmean(1:nw)= zero

do iigrid=1,igridstail; igrid=igrids(iigrid);
   level=node(plevel_,igrid)
   volumeflat(level)=volumeflat(level)+ &
          {(rnode(rpxmax^D_,igrid)-rnode(rpxmin^D_,igrid))|*}
   if (slab) then
      dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*}
   else
      dvolume(ixM^T)=pgeo(igrid)%dvolume(ixM^T)
      volume(level)=volume(level)+sum(dvolume(ixM^T))
   end if
   do iw=1,nw
      wmean(iw)=wmean(iw)+sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,iw))
   end do
end do
if (slab) volume(levmin:levmax)=volumeflat(levmin:levmax)

voltotal=sum(volume(levmin:levmax))

numlevels=levmax-levmin+1
dsum_send(1:nw)=wmean(1:nw)
dsum_send(nw+1)=voltotal
dsum_send(nw+2:nw+1+numlevels)=volumeflat(levmin:levmax)
call MPI_REDUCE(dsum_send,dsum_recv,nw+1+numlevels,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,icomm,ierrmpi)
!------------------------------------------------------------------
! Try to be smart: increase levels on shock when under some threshold:
if (nleafs_level(mxnest-addlevel_origin+addlevel)<nmin) then
        addlevel=addlevel+1
        if (mype==0) print*, "Increased refinement of shock to &
level", mxnest-addlevel_origin+addlevel
        typelimiter1(mxnest-addlevel_origin+addlevel-1) = 'minmod'
end if
if (nleafs_level(mxnest-addlevel_origin+addlevel)>nmax) then
        addlevel=addlevel-1
        if (mype==0) print*, "Decreased refinement of shock to &
level", mxnest-addlevel_origin+addlevel
        typelimiter1(mxnest-addlevel_origin+addlevel+1) = 'minmod'
end if
! set the limiter on shock to minmod:
typelimiter1(mxnest-addlevel_origin+addlevel:mxnest) = 'minmod'
!------------------------------------------------------------------

if (mype==0) then
! To compute cell updates per second, we do the following:
nx^D=ixMhi^D-ixMlo^D+1;
nc={nx^D*}
ncells = nc * nleafs_active
! assumes the number of active leafs haven't changed since last compute.
now        = MPI_WTIME()
dit        = it - itTimeLast
dtTimeLast = now - timeLast
itTimeLast = it
timeLast   = now
cellupdatesPerSecond = dble(ncells * nstep * dit) / (dtTimeLast * npe)
! blocks per core:
activeBlocksPerCore = nleafs_active / npe
! Wall clock time per code time unit in seconds:
wctPerCodeTime = dtTimeLast / (dble(dit) * dt)
! Wall clock time to finish in hours:
timeToFinish = (tmax - t) * wctPerCodeTime / 3600.0d0

   wmean(1:nw)=dsum_recv(1:nw)
   voltotal=dsum_recv(nw+1)
   volumeflat(levmin:levmax)=dsum_recv(nw+2:nw+1+numlevels)

   wmean=wmean/voltotal

   ! determine coverage in coordinate space
   volprob={(xprobmax^D-xprobmin^D)|*}
   volumeflat(levmin:levmax)=volumeflat(levmin:levmax)/volprob

   if (.not.opened) then
      ! generate filename
      write(filename,"(a,a)") TRIM(filenamelog),".log"

      amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      amode=ior(amode,MPI_MODE_APPEND)
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode, &
                         MPI_INFO_NULL,log_fh,ierrmpi)
      opened=.true.

      call MPI_FILE_WRITE(log_fh,fileheadout,len_trim(fileheadout), &
                          MPI_CHARACTER,status,ierrmpi)
      call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

      i=len_trim(wnames)-1
      do level=1,mxnest
          i=i+3
          if(level<10) then
            if (i+1<79) write(wnames(i:i+1),"(a,i1)") "c",level
          else
            if (i+2<79) write(wnames(i:i+2),"(a,i2)") "c",level
          endif
      end do

      do level=1,mxnest
          i=i+3
          if(level<10) then
            if (i+1<79) write(wnames(i:i+1),"(a,i1)") "n",level
          else
            if (i+2<79) write(wnames(i:i+2),"(a,i2)") "n",level
          endif
      end do

      if(residmin>smalldouble) then
        write(line,'(a15,a79)')"it   t  dt res ",wnames
      else
        write(line,'(a15,a79)')"it   t   dt    ",wnames
      endif

      line=trim(line)//"| Xload Xmemory 'Cell_Updates /second/core'"
      line=trim(line)//" 'Active_Blocks/Core' 'Wct Per Code Time [s]' 'TimeToFinish [hrs]'"

      call MPI_FILE_WRITE(log_fh,line,len_trim(line),MPI_CHARACTER, &
                          status,ierrmpi)
   end if
   call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

   if(residmin>smalldouble) then
      write(line,'(i7,3(es12.4))')it,t,dt,residual
   else
      write(line,'(i7,2(es12.4))')it,t,dt
   endif

   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                       MPI_CHARACTER,status,ierrmpi)
   do iw=1,nw
      write(line,'(es12.4)')wmean(iw)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do
   do level=1,mxnest
      write(line,'(es12.4)')volumeflat(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do

   do level=1,mxnest
      write(line,'(i8)') nleafs_level(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do

   write(line,'(a3,6(es10.2))') ' | ', Xload, Xmemory, cellupdatesPerSecond, &
        activeBlocksPerCore, wctPerCodeTime, timeToFinish
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)

end if

end subroutine printlog_special
!=============================================================================
subroutine userspecialconvert(qunitconvert)

! Allow user to use their own data-converting proceures

use constants
use mod_global_parameters
integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------
logical :: fileopen
integer :: iigrid, igrid, level, iw
double precision :: volume(1:nlevelshi),  voltotal
double precision :: re, ke, me, te
double precision :: dvolume(ixG^T), inBubble(ixG^T), lfac_pw(ixG^T), volumeflat(1:nlevelshi)
integer :: numlevels
integer, dimension(1:nlevelshi) :: isum_send
double precision, dimension(1:4) :: dsum_send, dsum_recv
character(len=80) :: filename
character(len=1024) :: line
logical, save :: opened=.false.
logical, save :: file_exists=.false.
integer :: amode, status(MPI_STATUS_SIZE)
double precision :: trcut
!-----------------------------------------------------------------------------

!!!Selects the bubble volume: !!!!!!!!!!
trcut = 1.0d-3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

volume(1:mxnest)=zero
volumeflat(1:mxnest)=zero
re = 0.0d0
ke = 0.0d0
me = 0.0d0
te = 0.0d0


do iigrid=1,igridstail; igrid=igrids(iigrid);

call primitive(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x)
! First normalize to cgs units:
   {^D&px(igrid)%x(ixM^T,^D) = normvar(0) * px(igrid)%x(ixM^T,^D)\}
   do iw = 1, nw
      pw(igrid)%w(ixM^T,iw) = normvar(iw) * pw(igrid)%w(ixM^T,iw)
   end do
! Mask out bubble volume:
   call get_lfac(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,lfac_pw)
   where (pw(igrid)%w(ixM^T,lfac_) .le. eqpar(flfac_)*lfac_pw(ixM^T) .and. &
       pw(igrid)%w(ixM^T,tr1_) .ge. trcut )
      inBubble(ixM^T) = one
   elsewhere
      inBubble(ixM^T) = zero
   end where

   level=node(plevel_,igrid)
   volumeflat(level)=volumeflat(level)+ &
          {(rnode(rpxmax^D_,igrid)-rnode(rpxmin^D_,igrid))|*} &
          * normvar(0)**3.0d0
   if (slab) then
      dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*} &
           * normvar(0)**3.0d0
   else
      dvolume(ixM^T)=pgeo(igrid)%dvolume(ixM^T) &
           * normvar(0)**3.0d0 * 2.0d0* dpi
      volume(level)=volume(level)+sum(dvolume(ixM^T))
   end if
! Now compute the energies:

   re = re + sum(pw(igrid)%w(ixM^T,lfac_)*pw(igrid)%w(ixM^T,rho_)*CONST_c**2.0d0 & 
        * dvolume(ixM^T) * inBubble(ixM^T))
   ke = ke + sum((pw(igrid)%w(ixM^T,lfac_) - 1.0d0) &
        * pw(igrid)%w(ixM^T,lfac_)*pw(igrid)%w(ixM^T,rho_)*CONST_c**2.0d0 &
        * dvolume(ixM^T) * inBubble(ixM^T))
   me = me + sum((({^C& pw(igrid)%w(ixM^T,b^C_)**2.0d0 |+})/4.0d0/dpi &
        - one/8.0d0/dpi/pw(igrid)%w(ixM^T,lfac_)**2.0d0 * &
        (({^C& pw(igrid)%w(ixM^T,b^C_)**2.0d0 |+}) &
        + ({^C& pw(igrid)%w(ixM^T,u^C_)/CONST_c * pw(igrid)%w(ixM^T,b^C_) |+})**2.0d0 )) &
        * dvolume(ixM^T) * inBubble(ixM^T))
   te = te + sum((4.0d0*pw(igrid)%w(ixM^T,lfac_)**2.0d0*pw(igrid)%w(ixM^T,pp_)  &
        - pw(igrid)%w(ixM^T,pp_) ) &
        * dvolume(ixM^T) * inBubble(ixM^T))

end do
if (slab) volume(levmin:levmax)=volumeflat(levmin:levmax)

voltotal=sum(volume(levmin:levmax))


dsum_send(1)=re
dsum_send(2)=ke
dsum_send(3)=me
dsum_send(4)=te

call MPI_REDUCE(dsum_send,dsum_recv,4,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,icomm,ierrmpi)


if (mype==0) then

   re = dsum_recv(1)
   ke = dsum_recv(2)
   me = dsum_recv(3)
   te = dsum_recv(4)

   if (.not.opened) then
      ! generate filename
      write(filename,"(a,a)") TRIM("energetics"),".csv"
      INQUIRE(FILE=filename, EXIST=file_exists)

      if (.not. file_exists) then
         open(unit=qunitconvert,file=filename,status='unknown',access='append')
         write(qunitconvert,"(a)") trim('# t [years] re [erg] ke [erg] me [erg] te [erg]')
      else
         open(unit=qunitconvert,file=filename,status='unknown',access='append')
      end if
      opened=.true.
   end if

   write(qunitconvert,'(5(es14.6))')t*UNIT_LENGTH/UNIT_VELOCITY/CONST_years,re,ke,me,te

end if


end subroutine userspecialconvert
!=============================================================================
subroutine get_lfac(ixI^L,ixO^L,w,x,lfac_pw)

      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: w(ixI^S,1:nw)
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      double precision, intent(out)   :: lfac_pw(ixI^S)

 lfac_pw(ixO^S) = eqpar(lor_)
 
end subroutine get_lfac
!=============================================================================
subroutine get_delta(ixI^L,ixO^L,w,x,delta)

      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: w(ixI^S,1:nw)
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      double precision, intent(out)   :: delta(ixI^S)
! ..local..
      double precision, dimension(ixI^S) :: rsphere, costheta,theta

if (typeaxial .ne. 'spherical') then
   rsphere(ixI^S)  = dsqrt({^D&dabs(x(ixI^S,^D))**2.0d0+})
{^IFTHREED
   costheta(ixI^S) = x(ixI^S,3)/rsphere(ixI^S)
}
{^IFTWOD
   costheta(ixI^S) = x(ixI^S,2)/rsphere(ixI^S)
}
   theta(ixI^S)    = acos(costheta(ixI^S))  
else
   rsphere(ixI^S)  = dabs(x(ixI^S,1))
   costheta(ixI^S) = cos(x(ixI^S,2))
   theta(ixI^S)    = x(ixI^S,2)
end if

    where (theta(ixO^S) < dpi/2.0d0-eqpar(eta_)&
    .or. theta(ixO^S) > dpi/2.0d0+eqpar(eta_))
        delta(ixO^S) = one
        elsewhere
    delta(ixO^S) = dabs(2.0d0/dpi*acos(-1.0d0/tan(theta(ixO^S))*1.0d0/tan(eqpar(eta_)) )-1.0d0)**2.0d0 
    end where

end subroutine get_delta
!=============================================================================
subroutine get_sigma(ixI^L,ixO^L,w,x,sigma)

      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: w(ixI^S,1:nw)
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      double precision, intent(out)   :: sigma(ixI^S)
! ..local..
      double precision, dimension(ixI^S) :: rsphere, costheta,theta,delta

if (typeaxial .ne. 'spherical') then
   rsphere(ixI^S)  = dsqrt({^D&dabs(x(ixI^S,^D))**2.0d0+})
{^IFTHREED
   costheta(ixI^S) = x(ixI^S,3)/rsphere(ixI^S)
}
{^IFTWOD
   costheta(ixI^S) = x(ixI^S,2)/rsphere(ixI^S)
}
   theta(ixI^S)    = acos(costheta(ixI^S))  
else
   rsphere(ixI^S)  = dabs(x(ixI^S,1))
   costheta(ixI^S) = cos(x(ixI^S,2))
   theta(ixI^S)    = x(ixI^S,2)
end if

   call get_delta(ixI^L,ixO^L,w,x,delta)

where (theta(ixO^S)>dpi/2.-eqpar(eta_) .and. theta(ixO^S)<dpi/2.+eqpar(eta_))
    sigma(ixO^S)= eqpar(sig0_)*delta(ixO^S)/&
    ( 1. + eqpar(sig0_)*(1.-delta(ixO^S)))
    elsewhere (theta(ixO^S)<= eqpar(theta0_))
       sigma(ixO^S)=eqpar(sig0_)*(theta(ixO^S)/eqpar(theta0_))**2.0d0
    elsewhere (theta(ixO^S)>= dpi-eqpar(theta0_))
       sigma(ixO^S)=eqpar(sig0_)*((dpi-theta(ixO^S))/eqpar(theta0_))**2.0d0
    elsewhere
       sigma(ixO^S)=eqpar(sig0_)
end where

end subroutine get_sigma
!=============================================================================
! The rest is just dummies
!=============================================================================









































!=============================================================================
subroutine bc_int(level,qt,ixI^L,ixO^L,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(in) :: x(ixI^S,1:ndim)

! .. local ..
!-----------------------------------------------------------------------------

end subroutine bc_int
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

use mod_global_parameters

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in &
      parfile or in user file')

var(ixI^S) = zero 

end subroutine specialvarforerrest
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0
!=============================================================================
!=============================================================================
