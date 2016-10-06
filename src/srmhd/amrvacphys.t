!#############################################################################
! Module amrvacphys/- srmhdeos version april 2009
! added positions to fluxes etc. march 2012

INCLUDE:amrvacnul/getdt.t
INCLUDE:amrvacnul/roe.t
!=============================================================================
subroutine checkglobaldata
use mod_global_parameters
!-----------------------------------------------------------------------------
minrho = max(zero,smallrho)
{#IFDEF ISO
minp=eqpar(adiab_)*minrho**eqpar(gamma_)
govergminone =eqpar(gamma_)/(eqpar(gamma_)-one)
}
{#IFDEF GAMMA
govergminone =eqpar(gamma_)/(eqpar(gamma_)-one)
minp  = max(zero,smallp)
smallxi=minrho+minp*govergminone
smalltau = minp/(eqpar(gamma_)-one)
}
{#IFDEF SYNGE
call smallvaluesEOS
}
end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

! place to set entropy fixes etc, absent for now

use mod_global_parameters
!-----------------------------------------------------------------------------
eqpar(adiab_)=1.0d0

{#IFNDEF SYNGE
eqpar(gamma_)=4.0d0/3.0d0
}
{#IFDEF SYNGE
eqpar(gamma_)=5.d0/3.d0
}

{#IFDEF GLM
eqpar(Cr_)= 0.18d0
}

if(strictzero)limitvalue=zero
if(.not.strictzero)limitvalue=smalldouble**2.0d0

end subroutine initglobaldata
!=============================================================================
subroutine checkw(checkprimitive,ixI^L,ixO^L,w,flag)

use mod_global_parameters
  
logical, intent(in)          :: checkprimitive
integer, intent(in)          :: ixI^L,ixO^L
double precision, intent(in) :: w(ixI^S,1:nw)
logical, intent(out)         :: flag(ixG^T)
!-----------------------------------------------------------------------------

flag(ixG^T)=.true.

{#IFDEF ENERGY
if (checkprimitive) then
  if(useprimitiveRel)then
     ! check   rho>=0, p>=smallp
     flag(ixO^S) = (w(ixO^S,rho_) > minrho).and. &
                   (w(ixO^S,pp_)  >=minp)
  else
    ! check  v^2 < 1, rho>=0, p>=minp
     ! v^2 < 1, rho>0, p>0
     flag(ixO^S) = (({^C&w(ixO^S,v^C_)**2.0d0+})< one).and. &
                   (w(ixO^S,rho_) > minrho).and. &
                   (w(ixO^S,pp_) >=minp)
  endif
else
  ! checks on the conservative variables
  flag(ixO^S)= (w(ixO^S,d_)   > minrho).and. &
               (w(ixO^S,tau_) > smalltau)
end if
}
{#IFNDEF ENERGY
if (checkprimitive) then
  if(useprimitiveRel)then
     ! check   rho>=0
     flag(ixO^S) = (w(ixO^S,rho_) > minrho)
  else
    ! check  v^2 < 1, rho>=0
     ! v^2 < 1, rho>0, p>0
     flag(ixO^S) = (({^C&w(ixO^S,v^C_)**2.0d0+})< one).and. &
                   (w(ixO^S,rho_) > minrho)
  endif
else
  ! checks on the conservative variables
  flag(ixO^S)= (w(ixO^S,d_)   > minrho)
end if
}

end subroutine checkw
!=============================================================================
subroutine conserve(ixI^L,ixO^L,w,x,patchw)

! Transform primitive variables into conservative ones
! (rho,v,p,B) ---> (D,S,tau,B,lfac,xi)
! call to smallvalues
! --> latter only used for correcting procedure in correctaux
! --> input array patchw for spatially selective transformation

use mod_global_parameters

integer, intent(in)               :: ixI^L,ixO^L
double precision, intent(inout)   :: w(ixI^S,nw)
logical, intent(in)               :: patchw(ixG^T)
double precision, intent(in)      :: x(ixI^S,1:ndim)

double precision, dimension(ixG^T):: sqrV,sqrU,sqrB,VdotB,rhoh
!-----------------------------------------------------------------------------

if(useprimitiveRel)then
  ! assumes four velocity computed in primitive (rho u p B) with u=lfac*v
  where(.not.patchw(ixO^S))
     sqrU(ixO^S) ={^C&w(ixO^S,u^C_)**2.0d0+} 
     sqrV(ixO^S) =sqrU(ixO^S)/(one+sqrU(ixO^S))
     sqrB(ixO^S) ={^C&w(ixO^S,b^C_)**2.0d0+}
  endwhere
else
  ! assumes velocity in primitive (rho v p B) 
  where(.not.patchw(ixO^S))
     sqrV(ixO^S) ={^C&w(ixO^S,v^C_)**2.0d0+}
     sqrB(ixO^S) ={^C&w(ixO^S,b^C_)**2.0d0+} 
  endwhere
endif

! fill the auxilary variable lfac (lorentz factor)
if(useprimitiveRel)then
   where(.not.patchw(ixO^S))
      w(ixO^S,lfac_)=dsqrt(one+sqrU(ixO^S)) 
      VdotB(ixO^S)  =({^C&w(ixO^S,b^C_)*w(ixO^S,u^C_)+})/w(ixO^S,lfac_)
   endwhere
else
   where(.not.patchw(ixO^S))
      w(ixO^S,lfac_)=one/dsqrt(one-sqrV(ixO^S))
      VdotB(ixO^S)  ={^C&w(ixO^S,b^C_)*w(ixO^S,v^C_)+}
   endwhere
endif

! fill the auxilary variable xi and density D
! with enthalpy w: xi= lfac^2 rhoh
! density: d = lfac * rho
call Enthalpy(w,ixI^L,ixO^L,patchw,rhoh)
where(.not.patchw(ixO^S))
   w(ixO^S,xi_)=w(ixO^S,lfac_)*w(ixO^S,lfac_)* rhoh(ixO^S)
   w(ixO^S,d_) =w(ixO^S,rho_)*w(ixO^S,lfac_)
   ! recycle sqrU array for storing temporary positive 
   ! array for use in energy variable
   sqrU(ixO^S)=sqrB(ixO^S)*sqrV(ixO^S)-VdotB(ixO^S)**2.0d0
endwhere
where(.not.patchw(ixO^S).and.sqrU(ixO^S)<zero)
   sqrU(ixO^S)=zero
endwhere

{#IFDEF TRACER
! We got D, now we can get the conserved tracers:
where(.not.patchw(ixO^S))
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,d_)*w(ixO^S,tr^FL_)\}
endwhere
}

{#IFDEF EPSINF
where(.not.patchw(ixO^S))
   w(ixO^S,Dn_) =w(ixO^S,n_)*w(ixO^S,lfac_)
   w(ixO^S,Dn0_) =w(ixO^S,Dn_)*w(ixO^S,n0_)
   w(ixO^S,Depsinf_) = w(ixO^S,epsinf_)*w(ixO^S,Dn_)**(2.0d0/3.0d0) &
        *w(ixO^S,lfac_)**(1.0d0/3.0d0)
   w(ixO^S,Depslow_) = w(ixO^S,epslow_)*w(ixO^S,Dn_)**(2.0d0/3.0d0) &
        *w(ixO^S,lfac_)**(1.0d0/3.0d0)
endwhere
}
! fill the vector S
! s= (xi + B^2) * v - (v.B) * B
if(useprimitiveRel)then
   where(.not.patchw(ixO^S))
       {^C&w(ixO^S,s^C_)=(w(ixO^S,xi_)+sqrB(ixO^S))&
               *w(ixO^S,u^C_)/w(ixO^S,lfac_) - &
                VdotB(ixO^S)*w(ixO^S,b^C_);}
   endwhere
else
   where(.not.patchw(ixO^S))
       {^C&w(ixO^S,s^C_)=(w(ixO^S,xi_)+sqrB(ixO^S))*w(ixO^S,v^C_) - &
                       VdotB(ixO^S)*w(ixO^S,b^C_);}
   endwhere
endif

!{#IFDEF ENERGY
! E = xi - p +B^2/2 + (v^2 B^2 - (v.B)^2)/2 
! instead of E use tau= E - D
where(.not.patchw(ixO^S))
    w(ixO^S,tau_)=w(ixO^S,xi_) - w(ixO^S,pp_) - w(ixO^S,d_) +& 
         half*(sqrB(ixO^S) + sqrU(ixO^S))
endwhere
!}
if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,patchw(ixO^S),'conserve')

end subroutine conserve
!=============================================================================
subroutine conserven(ixI^L,ixO^L,w,patchw)

! Transform primitive variables into conservative ones
! (rho,v,p,B) ---> (D,S,tau,B,lfac,xi)
! no call to smallvalues
! --> latter only used for correcting procedure in correctaux
! --> input array patchw for spatially selective transformation

use mod_global_parameters

integer, intent(in)               :: ixI^L,ixO^L
double precision, intent(inout)   :: w(ixI^S,nw)
logical, intent(in)               :: patchw(ixG^T)

double precision, dimension(ixG^T):: sqrV,sqrU,sqrB,VdotB,rhoh
!-----------------------------------------------------------------------------

if(useprimitiveRel)then
  ! assumes four velocity computed in primitive (rho u p B) with u=lfac*v
  where(.not.patchw(ixO^S))
     sqrU(ixO^S) ={^C&w(ixO^S,u^C_)**2.0d0+} 
     sqrV(ixO^S) =sqrU(ixO^S)/(one+sqrU(ixO^S))
     sqrB(ixO^S) ={^C&w(ixO^S,b^C_)**2.0d0+}
  endwhere
else
  ! assumes velocity in primitive (rho v p B) 
  where(.not.patchw(ixO^S))
     sqrV(ixO^S) ={^C&w(ixO^S,v^C_)**2.0d0+}
     sqrB(ixO^S) ={^C&w(ixO^S,b^C_)**2.0d0+} 
  endwhere
endif

! fill the auxilary variable lfac (lorentz factor)
if(useprimitiveRel)then
   where(.not.patchw(ixO^S))
      w(ixO^S,lfac_)=dsqrt(one+sqrU(ixO^S)) 
      VdotB(ixO^S)  =({^C&w(ixO^S,b^C_)*w(ixO^S,u^C_)+})/w(ixO^S,lfac_)
   endwhere
else
   where(.not.patchw(ixO^S))
      w(ixO^S,lfac_)=one/dsqrt(one-sqrV(ixO^S))
      VdotB(ixO^S)  ={^C&w(ixO^S,b^C_)*w(ixO^S,v^C_)+}
   endwhere
endif

! fill the auxilary variable xi and density D
! with enthalpy w: xi= lfac^2 rhoh
! density: d = lfac * rho
call Enthalpy(w,ixI^L,ixO^L,patchw,rhoh)
where(.not.patchw(ixO^S))
   w(ixO^S,xi_)=w(ixO^S,lfac_)*w(ixO^S,lfac_)* rhoh(ixO^S)
   w(ixO^S,d_) =w(ixO^S,rho_)*w(ixO^S,lfac_)
   ! recycle sqrU array for storing temporary positive 
   ! array for use in energy variable
   sqrU(ixO^S)=sqrB(ixO^S)*sqrV(ixO^S)-VdotB(ixO^S)**2.0d0
endwhere
where(.not.patchw(ixO^S).and.sqrU(ixO^S)<zero)
   sqrU(ixO^S)=zero
endwhere

{#IFDEF TRACER
! We got D, now we can get the conserved tracers:
where(.not.patchw(ixO^S))
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,d_)*w(ixO^S,tr^FL_)
\}
endwhere
}{#IFDEF EPSINF
where(.not.patchw(ixO^S))
   w(ixO^S,Dn_) =w(ixO^S,n_)*w(ixO^S,lfac_)
   w(ixO^S,Dn0_) =w(ixO^S,Dn_)*w(ixO^S,n0_)
   w(ixO^S,Depsinf_) = w(ixO^S,epsinf_)*w(ixO^S,Dn_)**(2.0d0/3.0d0) &
        *w(ixO^S,lfac_)**(1.0d0/3.0d0)
   w(ixO^S,Depslow_) = w(ixO^S,epslow_)*w(ixO^S,Dn_)**(2.0d0/3.0d0) &
        *w(ixO^S,lfac_)**(1.0d0/3.0d0)
endwhere
}
! fill the vector S
! s= (xi + B^2) * v - (v.B) * B
if(useprimitiveRel)then
   where(.not.patchw(ixO^S))
       {^C&w(ixO^S,s^C_)=(w(ixO^S,xi_)+sqrB(ixO^S))&
               *w(ixO^S,u^C_)/w(ixO^S,lfac_) - &
                VdotB(ixO^S)*w(ixO^S,b^C_);}
   endwhere
else
   where(.not.patchw(ixO^S))
       {^C&w(ixO^S,s^C_)=(w(ixO^S,xi_)+sqrB(ixO^S))*w(ixO^S,v^C_) - &
                       VdotB(ixO^S)*w(ixO^S,b^C_);}
   endwhere
endif

!{#IFDEF ENERGY
! E = xi - p +B^2/2 + (v^2 B^2 - (v.B)^2)/2 
! instead of E use tau= E - D
where(.not.patchw(ixO^S))
    w(ixO^S,tau_)=w(ixO^S,xi_) - w(ixO^S,pp_) - w(ixO^S,d_) +& 
         half*(sqrB(ixO^S) + sqrU(ixO^S))
endwhere
!}
end subroutine conserven
!=============================================================================
subroutine primitive(ixI^L,ixO^L,w,x)

! Transform conservative variables into primitive ones
! (D,S,tau,B)-->(rho,v,p,B,lfac,xi)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(inout)    :: w(ixI^S,1:nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
double precision, dimension(ixG^T) :: SdotB, sqrB,tmpP
!-----------------------------------------------------------------------------

! calculate lorentz factor and xi from conservatives only
call getaux(.true.,w,x,ixI^L,ixO^L,'primitive')

!{#IFDEF ENERGY
call Pressure(w,ixI^L,ixO^L,.true.,tmpP)
w(ixO^S,pp_)=tmpP(ixO^S)
!}

w(ixO^S,rho_)=w(ixO^S,d_)/w(ixO^S,lfac_)

SdotB(ixO^S)= {^C&w(ixO^S,s^C_)*w(ixO^S,b^C_)+}
sqrB(ixO^S) = {^C&w(ixO^S,b^C_)**2.0d0+}  

if(useprimitiveRel)then   
  where (w(ixO^S,xi_)<smallxi)
     {^C&w(ixO^S,u^C_)=zero\}
  elsewhere
     {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*(w(ixO^S,s^C_)&
                      +SdotB(ixO^S)*w(ixO^S,b^C_)/w(ixO^S,xi_))/ &
                      (w(ixO^S,xi_)+sqrB(ixO^S))\}
  end where
else
  where (w(ixO^S,xi_)<smallxi)
     {^C&w(ixO^S,v^C_)=zero\}
  elsewhere
     {^C&w(ixO^S,v^C_)=(w(ixO^S,s^C_)+SdotB(ixO^S)*w(ixO^S,b^C_)/w(ixO^S,xi_))&
                      /(w(ixO^S,xi_)+sqrB(ixO^S))\}
  endwhere
endif

{#IFDEF TRACER
! We got lor, rho, Dtr, now we can get the tracers:
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,Dtr^FL_)/w(ixO^S,lfac_)/w(ixO^S,rho_)\}
}

{#IFDEF EPSINF
w(ixO^S,n_)   = w(ixO^S,Dn_)/w(ixO^S,lfac_)
w(ixO^S,n0_)   = w(ixO^S,Dn0_)/w(ixO^S,lfac_)/w(ixO^S,n_)
w(ixO^S,epsinf_) = w(ixO^S,Depsinf_)/(w(ixO^S,lfac_) &
     *w(ixO^S,n_)**(2.0d0/3.0d0))
w(ixO^S,epslow_) = w(ixO^S,Depslow_)/(w(ixO^S,lfac_) &
     *w(ixO^S,n_)**(2.0d0/3.0d0))
}

{#IFDEF ENERGY
if (tlow>zero) call fixp_usr(ixI^L,ixO^L,w,x)
}

end subroutine primitive
!==============================================================================
subroutine primitiven(ixI^L,ixO^L,w,patchw)

! same as primitive, but with padding by patchw
!  --> needed in correctaux, avoiding recursive calling by smallvalues
! assumes updated xi and lfac, and does not use smallxi cutoff
! Transform conservative variables into primitive ones
! (D,S,tau,B)-->(rho,v,p,B,lfac,xi)

use mod_global_parameters

integer, intent(in)                    :: ixI^L,ixO^L
double precision, intent(inout)        :: w(ixI^S,1:nw)
logical, intent(in),dimension(ixG^T)   :: patchw
double precision, dimension(ixG^T) :: SdotB, sqrB,tmpP
!-----------------------------------------------------------------------------

!{#IFDEF ENERGY
call Pressuren(w,ixI^L,ixO^L,.true.,tmpP,patchw)
where(.not.patchw(ixO^S))
   w(ixO^S,pp_)=tmpP(ixO^S)
end where
!}

where(.not.patchw(ixO^S))
   w(ixO^S,rho_)=w(ixO^S,d_)/w(ixO^S,lfac_)
   SdotB(ixO^S)= {^C&w(ixO^S,s^C_)*w(ixO^S,b^C_)+}
   sqrB(ixO^S) = {^C&w(ixO^S,b^C_)**2.0d0+}
end where


if(useprimitiveRel)then   
  where(.not.patchw(ixO^S))
     {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*(w(ixO^S,s^C_)&
                      +SdotB(ixO^S)*w(ixO^S,b^C_)/w(ixO^S,xi_))/ &
                      (w(ixO^S,xi_)+sqrB(ixO^S))\}
  end where
else
  where(.not.patchw(ixO^S))
     {^C&w(ixO^S,v^C_)=(w(ixO^S,s^C_)+SdotB(ixO^S)*w(ixO^S,b^C_)/w(ixO^S,xi_))&
                      /(w(ixO^S,xi_)+sqrB(ixO^S))\}
  endwhere
endif

{#IFDEF TRACER
! We got lor, rho, Dtr, now we can get the tracers:
where(.not.patchw(ixO^S))
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,Dtr^FL_)/w(ixO^S,lfac_)/w(ixO^S,rho_)\}
endwhere
}
{#IFDEF EPSINF
where(.not.patchw(ixO^S))
w(ixO^S,n_)   = w(ixO^S,Dn_)/w(ixO^S,lfac_)
w(ixO^S,n0_)   = w(ixO^S,Dn0_)/w(ixO^S,lfac_)/w(ixO^S,n_)
w(ixO^S,epsinf_) = w(ixO^S,Depsinf_)/(w(ixO^S,lfac_) &
     *w(ixO^S,n_)**(2.0d0/3.0d0))
w(ixO^S,epslow_) = w(ixO^S,Depslow_)/(w(ixO^S,lfac_) &
     *w(ixO^S,n_)**(2.0d0/3.0d0))
end where
}
end subroutine primitiven
!==============================================================================
subroutine e_to_rhos(ixI^L,ixO^L,w,x)

use mod_global_parameters
integer:: ixI^L,ixO^L
double precision:: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-------------------------------------------------------------------------

call mpistop("e to rhos for SRMHDEOS unavailable")

end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixI^L,ixO^L,w,x)

use mod_global_parameters

integer:: ixI^L,ixO^L
double precision:: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

call mpistop("e to rhos for SRMHDEOS unavailable")

end subroutine rhos_to_e
!=============================================================================
subroutine ppmflatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L,ixL^L,ixR^L
double precision, intent(in)  :: w(ixI^S,nw),d2w(ixG^T,1:nwflux)

double precision, intent(inout) :: drho(ixG^T),dp(ixG^T)

double precision     :: ptot(ixG^T)
!-----------------------------------------------------------------------------

{#IFDEF ENERGY
if(useprimitive)then
   ! note that ixI range is ultimately ixG full range
   if (useprimitiveRel) then
      ptot(ixI^S) = w(ixI^S,pp_) &
                  + half*(({^C&w(ixI^S,u^C_)*w(ixI^S,b^C_)+})/w(ixI^S,lfac_)  &
                  +({^C&w(ixI^S,b^C_)**2.0d0+})/w(ixI^S,lfac_)**2.0d0)
   else
      ptot(ixI^S) = w(ixI^S,pp_)+half*(({^C&w(ixI^S,v^C_)*w(ixI^S,b^C_)+})  &
                  +({^C&w(ixI^S,b^C_)**2.0d0+})/w(ixI^S,lfac_)**2.0d0)
   end if
   call Calcule_Geffect(w,ixI^L,ixO^L,.false.,dp)
   drho(ixO^S) =dp(ixO^S)*dabs(d2w(ixO^S,rho_))&
              /min(w(ixL^S,rho_),w(ixR^S,rho_))
   dp(ixO^S) = dabs(ptot(ixR^S)-ptot(ixL^S))/min(ptot(ixL^S),ptot(ixR^S))
end if
}
{#IFDEF ISO
call mpistop("PPM with flatsh=.true. can not be used with eos=iso !")
}

end subroutine ppmflatcd
!=============================================================================
subroutine ppmflatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
integer, intent(in)           :: idims
double precision, intent(in)  :: w(ixI^S,nw)

double precision, intent(inout) :: drho(ixG^T),dp(ixG^T),dv(ixG^T)

double precision :: v(ixG^T),ptot(ixG^T)
!-----------------------------------------------------------------------------
{#IFDEF ENERGY
if(useprimitive)then
   ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
   ! note that ixI range is ultimately ixG full range
   if (useprimitiveRel) then
      ptot(ixI^S) = w(ixI^S,pp_) &
                  + half*(({^C&w(ixI^S,u^C_)*w(ixI^S,b^C_)+})/w(ixI^S,lfac_)  &
                  +({^C&w(ixI^S,b^C_)**2.0d0+})/w(ixI^S,lfac_)**2.0d0)
   else
      ptot(ixI^S) = w(ixI^S,pp_)+half*(({^C&w(ixI^S,v^C_)*w(ixI^S,b^C_)+})  &
                  +({^C&w(ixI^S,b^C_)**2.0d0+})/w(ixI^S,lfac_)**2.0d0)
   end if
   where (dabs(ptot(ixRR^S)-ptot(ixLL^S))>smalldouble)
      drho(ixO^S) = dabs((ptot(ixR^S)-ptot(ixL^S))&
                   /(ptot(ixRR^S)-ptot(ixLL^S)))
   elsewhere
      drho(ixO^S) = zero
   end where

   !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26 
   !  use "dp" to save sound speed squared
   call getcsound2(w,ixI^L,ixO^L,.false.,dp) 
   dp(ixO^S) = dabs(ptot(ixR^S)-ptot(ixL^S))&
                /(w(ixO^S,rho_)*dp(ixO^S))
   if (useprimitiveRel) then
     v(ixI^S)  = w(ixI^S,u0_+idims)/w(ixI^S,lfac_)
   else
     v(ixI^S)  = w(ixI^S,v0_+idims)
   end if
   call gradient(v,ixI^L,ixO^L,idims,dv)
end if
}
{#IFDEF ISO
call mpistop("PPM with flatsh=.true. can not be used with eos=iso !")
}
end subroutine ppmflatsh
!=============================================================================
subroutine getv(w,x,ixI^L,ixO^L,idims,v)

! Calculate v_idim=m_idim/rho within ixO^L

use mod_global_parameters
  
integer, intent(in)                              :: ixI^L,ixO^L,idims
double precision, intent(in)                     :: w(ixI^S,1:nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
double precision, dimension(ixG^T), intent(out)  :: v

double precision, dimension(ixG^T)               ::  VdotB, sqrB
!-----------------------------------------------------------------------------

VdotB(ixO^S)= ({^C&w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/w(ixO^S,xi_)
sqrB(ixO^S) = {^C&w(ixO^S,b^C_)**2.0d0+}  
where(w(ixO^S,xi_)+sqrB(ixO^S)<=smallxi)
   v(ixO^S)=zero
elsewhere
v(ixO^S)=(w(ixO^S,s0_+idims)+VdotB(ixO^S)*w(ixO^S,b0_+idims))/ &
         (w(ixO^S,xi_)+sqrB(ixO^S))
endwhere

end subroutine getv
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixI^L,ixO^L,idims,cmax,cmin,needcmin)
  
! Calculate cmax_idim within ixO^L
! new_cmax is not used
  
use mod_global_parameters
  
logical, intent(in)           :: new_cmax,needcmin
integer, intent(in)           :: ixI^L,ixO^L,idims
double precision, intent(in)  :: w(ixI^S,1:nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
double precision, intent(out) :: cmax(ixG^T),cmin(ixG^T)


integer,parameter:: order=4
integer          :: ix^D,iroot
double precision :: VdotB,sqrlfac,dedp,vi,onemvi,root,xscaleinv,lfac
double precision :: poly(0:order),bi,bivdotb,biolfac2
double precision :: alpha,beta,delta,zeta

double precision, dimension(1:order):: roots,lambdaroots
double precision                    :: lambdamax,lambdamin
double precision, dimension(ixG^T)  :: sqrB,rhoh,csound2
double precision, dimension(ixG^T)  :: vidim,vidim2,v2
double precision, dimension(ixG^T)  :: A,B,C,Mag
logical, dimension(ixG^T)           :: Hydro,s2
logical, dimension(ixG^T,5)         :: Domain
!-----------------------------------------------------------------------------


roots(1:order)=zero
lambdaroots(1:order)=zero
sqrB(ixO^S)  = {^C&w(ixO^S,b^C_)**2.0d0+} 


!-------------!
Hydro(ixO^S) = (sqrB(ixO^S)<=limitvalue)
s2(ixO^S)    = ({^C&w(ixO^S,s^C_)**2.0d0+})<=limitvalue

! the square of sound speed using EOS in eos.t
call getcsound2(w,ixI^L,ixO^L,.true.,csound2)


!---- Various domains ----!
Domain(ixO^S,1) = .not.Hydro(ixO^S) .and. .not.s2(ixO^S) &
.and. w(ixO^S,b0_+idims) /= zero
Domain(ixO^S,2) = .not.Hydro(ixO^S) .and.      s2(ixO^S) &
.and. w(ixO^S,b0_+idims) /= zero
Domain(ixO^S,3) = .not.Hydro(ixO^S) .and. .not.s2(ixO^S) &
 .and. w(ixO^S,b0_+idims) == zero 
Domain(ixO^S,4) = .not.Hydro(ixO^S) .and.      s2(ixO^S) &
 .and. w(ixO^S,b0_+idims) == zero 
Domain(ixO^S,5) =      Hydro(ixO^S)
!-------------------------!


if(ANY(Domain(ixO^S,1)))then
 select case(typepoly)
 case('original')
  {do ix^D= ixO^LIM^D\}
   if(Domain(ix^D,1))then
      VdotB       = ({^C&w(ix^D,s^C_)*w(ix^D,b^C_)+})/w(ix^D,xi_)
      sqrlfac     = w(ix^D,lfac_)**2.0d0
        
      dedp =  one/csound2(ix^D)

      delta = sqrlfac*w(ix^D,xi_)*(dedp-one)
      alpha  = w(ix^D,xi_)+VdotB*VdotB*sqrlfac*(dedp-one)+sqrB(ix^D)*dedp
        
      vi = (w(ix^D,s0_+idims)+VdotB*w(ix^D,b0_+idims))/(w(ix^D,xi_)+sqrB(ix^D))
      onemvi= one-vi
      bi = w(ix^D,b0_+idims)
      bivdotb=bi*VdotB
      biolfac2=bi**2/sqrlfac

      xscaleinv=one/(alpha+delta)
     ! coefficients of the non-transformed polynomial 
     ! multiplied by xscaleinv
      poly(0)= xscaleinv*(vi**4*delta-vi**2*alpha+two*vi*bivdotb+biolfac2)
      poly(1)= xscaleinv*(-4.0d0*delta*vi**3+two*vi*alpha-two*bivdotb)
      poly(2)= xscaleinv*(6.0d0*vi**2*delta+alpha*(vi**2-one)-two*vi*bivdotb&
               -biolfac2)
      poly(3)= xscaleinv*(-4.0d0*vi*delta-two*vi*alpha+two*bivdotb)
      poly(4)= one

      call zroots4(poly,roots)
      lambdaroots=roots ! transform back to lambda-space
      if(needcmin)then
        lambdamax = maxval(lambdaroots) 
        lambdamin = minval(lambdaroots) 
        cmax(ix^D)=min(max(lambdamax,zero),one)
        cmin(ix^D)=max(min(lambdamin,zero),-one)
      else
        lambdamax =maxval(dabs(lambdaroots)) 
        cmax(ix^D)=min(one,lambdamax)
      end if
 
      !if (lambdamax>one) then 
      !   print*,'WARNING: getcmax found eigenvalue >= 1 !'
      !   print*,ix^D
      !   print*,'roots in mu: ',roots
      !   do iroot=1,order
      !     root=roots(iroot)
      !     print*,iroot,root,poly(0)+poly(1)*root+poly(2)*root**2&
      !            +poly(3)*root**3+poly(4)*root**4
      !   enddo
      !   print*,'eigenvalues: ',lambdaroots
      !   print*,'lambdamax=',lambdamax
         !write(*,'(a1,5(f20.12,a1))') &
         !    '{',poly(0),',',poly(1),',',poly(2),',',poly(3),',',poly(4),'}'
         !print*,'d=',w(ix^D,d_),'s=',{^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
         !print*,'B=',{^C&w(ix^D,b^C_)},'lfac=',w(ix^D,lfac_),'xi=',w(ix^D,xi_)
      !   call mpistop("== getcmax ==")
      !end if
        

    endif ! Domain 1 original
  {enddo^D&\}
 case('bergmans')
  {do ix^D= ixO^LIM^D\}
   if(Domain(ix^D,1))then
      VdotB = ({^C&w(ix^D,s^C_)*w(ix^D,b^C_)+})/w(ix^D,xi_)
      sqrlfac = w(ix^D,lfac_)**2.0d0
      dedp =  one/csound2(ix^D)
      delta = sqrlfac*w(ix^D,xi_)*(dedp-one)
      alpha  = w(ix^D,xi_)+VdotB*VdotB*sqrlfac*(dedp-one)+sqrB(ix^D)*dedp      
      vi = (w(ix^D,s0_+idims)+VdotB*w(ix^D,b0_+idims))/(w(ix^D,xi_)+sqrB(ix^D))
      onemvi= one-vi
      bi = w(ix^D,b0_+idims)
      bivdotb=bi*VdotB
      biolfac2=bi**2/sqrlfac

      xscaleinv=one/(alpha+delta)
     ! coefficients of the transformed polynomial with 
     ! lambda -> 1-1/y , multiplied by y^4 and xscaleinv
      poly(0)= one
      poly(1)= xscaleinv*((-4.0d0*delta-two*alpha)*onemvi-two*(alpha+bivdotb))
      poly(2)= xscaleinv*(6.0d0*onemvi*(delta*onemvi+alpha) &
                 -alpha*(one-vi**2)+6.0d0*bivdotb-two*vi*bivdotb-biolfac2)
      poly(3)= xscaleinv*(-two*onemvi &
              *(two*delta*onemvi**2+alpha*onemvi+two*bivdotb)+two*biolfac2)
      poly(4)= xscaleinv*delta*onemvi**4

      call zroots4(poly,roots)
      lambdaroots=one-one/roots ! transform back to lambda-space
      if(needcmin)then
        lambdamax = maxval(lambdaroots) 
        lambdamin = minval(lambdaroots) 
        cmax(ix^D)=min(max(lambdamax,zero),one)
        cmin(ix^D)=max(min(lambdamin,zero),-one)
      else
        lambdamax = maxval(dabs(lambdaroots)) 
        cmax(ix^D)=min(one,lambdamax)
      end if
 
      !if (lambdamax>one) then 
      !   print*,'WARNING: getcmax found eigenvalue >= 1 !'
      !   print*,ix^D
      !   print*,'roots in mu: ',roots
      !   do iroot=1,order
      !     root=roots(iroot)
      !     print*,iroot,root,poly(0)+poly(1)*root+poly(2)*root**2&
      !            +poly(3)*root**3+poly(4)*root**4
      !   enddo
      !   print*,'eigenvalues: ',lambdaroots
      !   print*,'lambdamax=',lambdamax
         !write(*,'(a1,5(f20.12,a1))') &
         !    '{',poly(0),',',poly(1),',',poly(2),',',poly(3),',',poly(4),'}'
         !print*,'d=',w(ix^D,d_),'s=',{^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
         !print*,'B=',{^C&w(ix^D,b^C_)},'lfac=',w(ix^D,lfac_),'xi=',w(ix^D,xi_)
      !   call mpistop("== getcmax ==")
      !end if

   endif ! Domain 1 bergmans
  {enddo^D&\}
 case('meliani')
  {do ix^D= ixO^LIM^D\}
   if(Domain(ix^D,1))then
    VdotB = ({^C&w(ix^D,s^C_)*w(ix^D,b^C_)+})/w(ix^D,xi_)
    sqrlfac = w(ix^D,lfac_)**2.0d0
    dedp =  one/csound2(ix^D)
    vi = (w(ix^D,s0_+idims)+VdotB*w(ix^D,b0_+idims))/(w(ix^D,xi_)+sqrB(ix^D))
    onemvi= one-vi    
    bi = w(ix^D,b0_+idims)
    bivdotb=bi*VdotB
    biolfac2=bi**2/sqrlfac

    alpha = sqrlfac*w(ix^D,xi_)*(dedp-one)
    beta  = w(ix^D,xi_)+VdotB*VdotB*sqrlfac*(dedp-one)+sqrB(ix^D)*dedp
    delta = two*w(ix^D,lfac_)* bivdotb
    zeta  = - bi**2.0d0        
    xscaleinv=one/(alpha+beta)

   !*********************************************************!
   ! coefficients of the transformed polynomial with 
   ! lambda -> lfac (y-vi) , multiplied by y^4 and xscaleinv
    poly(4) = one
    poly(3) = xscaleinv*(two*vi*w(ix^D,lfac_)* beta+delta)
    poly(2) = xscaleinv*(two*vi*w(ix^D,lfac_)* delta+zeta&
             -sqrlfac*beta*(one-vi**2))
    poly(1) = xscaleinv*(two*vi*w(ix^D,lfac_)*zeta&
             -sqrlfac*delta*(one-vi**2))
    poly(0) = xscaleinv*(bi**2.0d0*sqrlfac*(one-vi**2))
   !*********************************************************!
    call zroots4(poly,roots)
    lambdaroots=roots/w(ix^D,lfac_)+vi
    if(needcmin)then
      lambdamax = maxval(lambdaroots) 
      lambdamin = minval(lambdaroots) 
      cmax(ix^D)=min(max(lambdamax,zero),one)
      cmin(ix^D)=max(min(lambdamin,zero),-one)
    else
      lambdamax = maxval(dabs(lambdaroots)) 
      cmax(ix^D)=min(one,lambdamax)
    end if

      !if (lambdamax>one) then 
      !   print*,'WARNING: getcmax found eigenvalue >= 1 !'
      !   print*,ix^D
      !   print*,'roots in mu: ',roots
      !   do iroot=1,order
      !     root=roots(iroot)
      !     print*,iroot,root,poly(0)+poly(1)*root+poly(2)*root**2&
      !            +poly(3)*root**3+poly(4)*root**4
      !   enddo
      !   print*,'eigenvalues: ',lambdaroots
      !   print*,'lambdamax=',lambdamax
         !write(*,'(a1,5(f20.12,a1))') &
         !    '{',poly(0),',',poly(1),',',poly(2),',',poly(3),',',poly(4),'}'
         !print*,'d=',w(ix^D,d_),'s=',{^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
         !print*,'B=',{^C&w(ix^D,b^C_)},'lfac=',w(ix^D,lfac_),'xi=',w(ix^D,xi_)
      !   call mpistop("== getcmax ==")
      !end if

   endif ! Domain 1 meliani
  {enddo^D&\}
 case('gammie')

   where(Domain(ixO^S,1))
     rhoh(ixO^S)   = w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
     ! squared Alfven speed
     B(ixO^S)  = sqrB(ixO^S)/(sqrB(ixO^S)+rhoh(ixO^S))
     ! equation 72 Del zanna et al. 2007A&A...473...11D
     A(ixO^S)  = csound2(ixO^S) + B(ixO^S)- csound2(ixO^S)*B(ixO^S)
     ! C= VdotB
     C(ixO^S)= (({^C&w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/w(ixO^S,xi_))
   endwhere

   needming : if(.not.needcmin)then

     where(Domain(ixO^S,1))
        vidim2(ixO^S)  = ((w(ixO^S,s0_+idims)+C(ixO^S)*w(ixO^S,b0_+idims))&
                      /(w(ixO^S,xi_)+sqrB(ixO^S)))**2.0d0
     end where
     directg: if (ndir==1) then
        where(Domain(ixO^S,1))
          cmax(ixO^S)= (dsqrt(vidim2(ixO^S))+dsqrt(A(ixO^S)))/ & 
                (one+dsqrt(A(ixO^S)*vidim2(ixO^S)))
        endwhere
     else ! directg
        where(Domain(ixO^S,1))
           v2(ixO^S)  = ({^C&(w(ixO^S,s^C_)+C(ixO^S)*w(ixO^S,b^C_))**2.0d0+})&
                     /(w(ixO^S,xi_)+sqrB(ixO^S))**2.0d0

           cmax(ixO^S)=min(one,max((dsqrt(vidim2(ixO^S))*(one-A(ixO^S)) + &
                dsqrt(A(ixO^S)*(one-v2(ixO^S))*((one-v2(ixO^S)&
                *A(ixO^S))-vidim2(ixO^S)*(one-A(ixO^S))))) &
                /(one-v2(ixO^S)*A(ixO^S)),zero))
        endwhere
     endif directg
 
   else ! needming

     where(Domain(ixO^S,1))
       vidim(ixO^S) =  (w(ixO^S,s0_+idims)+C(ixO^S)*w(ixO^S,b0_+idims))&
                      /(w(ixO^S,xi_)+sqrB(ixO^S))
     end where

     ! characteristic speed Eq. 76 Del zanna et al. 2007A&A...473...11D
     if (ndir==1) then
        where(Domain(ixO^S,1))
          cmax(ixO^S)= min( one,max((vidim(ixO^S)+dsqrt(A(ixO^S)))/ & 
                  (one+dsqrt(A(ixO^S))*vidim(ixO^S)),zero))
          cmin(ixO^S)= max(-one,min((vidim(ixO^S)-dsqrt(A(ixO^S)))/ & 
                  (one-dsqrt(A(ixO^S))*vidim(ixO^S)),zero))
        endwhere
     else ! ndir
        where(Domain(ixO^S,1))
          v2(ixO^S)    = ({^C&(w(ixO^S,s^C_)+C(ixO^S)*w(ixO^S,b^C_))**2.0d0+})&
                     /(w(ixO^S,xi_)+sqrB(ixO^S))**2.0d0
          cmax(ixO^S)=min( one,max((vidim(ixO^S)*(one-A(ixO^S))  &
                 +dsqrt(A(ixO^S)*(one-v2(ixO^S))&
                 *(one-v2(ixO^S)*A(ixO^S)-vidim(ixO^S)**2.0d0&
                 *(one-A(ixO^S)))))/(one-v2(ixO^S)*A(ixO^S)),zero))
          cmin(ixO^S)=max(-one,min((vidim(ixO^S)*(one-A(ixO^S))  &
                 -dsqrt(A(ixO^S)*(one-v2(ixO^S))&
                 *(one-v2(ixO^S)*A(ixO^S)-vidim(ixO^S)**2.0d0&
                 *(one-A(ixO^S)))))/(one-v2(ixO^S)*A(ixO^S)),zero))
        endwhere
     endif ! ndir
   end if needming
 case default
   call mpistop("Error in getcmax: Unknown polynomial type")
 end select
endif ! any domaine 1 

!--------------------------------------------------------------------------!
!Vanishing total velocity
! biquadratic equation: A cmax^4-B cmax^2 +C=0
s2zero : if(any(Domain(ixO^S,2)))then
 where(Domain(ixO^S,2))
  rhoh(ixO^S) = w(ixO^S,xi_)

  A(ixO^S) = rhoh(ixO^S)+ sqrB(ixO^S)
  B(ixO^S) = sqrB(ixO^S)+(rhoh(ixO^S)+w(ixO^S,b0_+idims)**2.0d0)*csound2(ixO^S)
  C(ixO^S) = csound2(ixO^S) * w(ixO^S,b0_+idims)**2.0d0
  cmax(ixO^S)= min(dabs(dsqrt((B(ixO^S)+&
        dsqrt(B(ixO^S)**2.0d0-4.0d0*A(ixO^S)*C(ixO^S)))/(2.0d0*A(ixO^S)))),one)
 endwhere
 if(needcmin)then
  where(Domain(ixO^S,2))
  cmin(ixO^S)= max(min(dsqrt((B(ixO^S)-&
               dsqrt(B(ixO^S)**2.0d0-4.0d0*A(ixO^S)*C(ixO^S)))&
               /(2.0d0*A(ixO^S))),zero),-one)
  endwhere
 endif 
endif s2zero
!--------------------------------------------------------------------------!
! Vanishing normal component
! quadratic equation  A cmax^2 - B cmax + C = 0
bxzero : if(any(Domain(ixO^S,3)))then
 if(needcmin)then
  where(Domain(ixO^S,3))
   rhoh(ixO^S)   = w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
   vidim(ixO^S)  = (w(ixO^S,s0_+idims)/(w(ixO^S,xi_)+sqrB(ixO^S)))

   Mag(ixO^S) = sqrB(ixO^S)/w(ixO^S,lfac_)**2.0d0 &
              +(({^C&w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/w(ixO^S,xi_))**2.0d0&
              *(one-csound2(ixO^S))

   A(ixO^S) = rhoh(ixO^S) *(csound2(ixO^S)&
             +w(ixO^S,lfac_)**2.0d0*(one-csound2(ixO^S))) + Mag(ixO^S)
   B(ixO^S) =   w(ixO^S,xi_)* vidim(ixO^S)* (one-csound2(ixO^S))
   C(ixO^S) = rhoh(ixO^S) * (-csound2(ixO^S) &
             +w(ixO^S,lfac_)**2.0d0*vidim(ixO^S)**2.0d0*(one-csound2(ixO^S))) &
             -Mag(ixO^S)

   cmax(ixO^S)= min(max((B(ixO^S)+dsqrt(B(ixO^S)**2.0d0-A(ixO^S)*C(ixO^S)))/&
                A(ixO^S),zero),  one)
   cmin(ixO^S)= max(min((B(ixO^S)-dsqrt(B(ixO^S)**2.0d0-A(ixO^S)*C(ixO^S)))/&
                A(ixO^S),zero), -one)
  endwhere
 else
  where(Domain(ixO^S,3))
   rhoh(ixO^S)   = w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
   vidim2(ixO^S)  = (w(ixO^S,s0_+idims)/(w(ixO^S,xi_)+sqrB(ixO^S)))**2.0d0


   Mag(ixO^S) = sqrB(ixO^S)/w(ixO^S,lfac_)**2.0d0 &
               +(({^C&w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/w(ixO^S,xi_))**2.0d0&
               *(one-csound2(ixO^S))

   A(ixO^S)   = rhoh(ixO^S) *(csound2(ixO^S)&
               +w(ixO^S,lfac_)**2.0d0*(one-csound2(ixO^S))) + Mag(ixO^S)
   B(ixO^S)   = w(ixO^S,xi_)* Dsqrt(vidim2(ixO^S))* (one-csound2(ixO^S))
   C(ixO^S)   = rhoh(ixO^S) * (-csound2(ixO^S) +w(ixO^S,lfac_)**2.0d0&
              *vidim2(ixO^S)*(one-csound2(ixO^S))) - Mag(ixO^S)

   cmax(ixO^S)= min(max((B(ixO^S)+dsqrt(B(ixO^S)**2.0d0-A(ixO^S)*C(ixO^S)))&
              /A(ixO^S),zero), one)

  endwhere
 end if
endif bxzero
!--------------------------------------------------------------------------!
! Vanishing normal component and speed
! A cmax - C =0
bxs2zero :if(Any(Domain(ixO^S,4)))then
 where(Domain(ixO^S,4))
  rhoh(ixO^S)   = w(ixO^S,xi_)

  A(ixO^S) = rhoh(ixO^S) +sqrB(ixO^S)
  C(ixO^S) = rhoh(ixO^S) * csound2(ixO^S)+sqrB(ixO^S)

  cmax(ixO^S)= min(dsqrt(C(ixO^S)/A(ixO^S)), one)
 endwhere
 if(needcmin)then
  where(Domain(ixO^S,4))
     cmin(ixO^S)= max(-dsqrt(C(ixO^S)/A(ixO^S)), -one)
  end where
 end if
endif bxs2zero
!--------------------------------------------------------------------------!
! Vanishing magnetic field
bbzero : if(Any(Domain(ixO^S,5)))Then
 needmin : if(.not.needcmin)then
 
  where(Domain(ixO^S,5))
   vidim2(ixO^S)  = (w(ixO^S,s0_+idims)/w(ixO^S,xi_))**2.0d0
   v2(ixO^S)      = ({^C&w(ixO^S,s^C_)**2.0d0+})/w(ixO^S,xi_)**2.0d0
  endwhere

  direct: if (ndir==1) then
   where(Domain(ixO^S,5))
    cmax(ixO^S)= (dsqrt(v2(ixO^S))+dsqrt(csound2(ixO^S)))/ & 
                (one+dsqrt(csound2(ixO^S)*v2(ixO^S)))
   endwhere
  else ! direct
   where(Domain(ixO^S,5))
    cmax(ixO^S)=max(( Dsqrt(vidim2(ixO^S))*(one-csound2(ixO^S)) + &
                dsqrt(csound2(ixO^S)*(one-v2(ixO^S))*(one-v2(ixO^S)&
                *csound2(ixO^S)-vidim2(ixO^S)*(one-csound2(ixO^S))))) &
                /(one-v2(ixO^S)*csound2(ixO^S)),zero) 
   endwhere
  endif direct
 else
   where(Domain(ixO^S,5))
    vidim(ixO^S)   = w(ixO^S,s0_+idims)/w(ixO^S,xi_)
    v2(ixO^S)      = ({^C&w(ixO^S,s^C_)**2.0d0+})/w(ixO^S,xi_)**2.0d0
  endwhere


  if (ndir==1) then
     where(Domain(ixO^S,5))
      cmax(ixO^S)= Max((vidim(ixO^S)+dsqrt(csound2(ixO^S)))/ & 
                  (one+dsqrt(csound2(ixO^S))*vidim(ixO^S)),zero)
      cmin(ixO^S)= Min((vidim(ixO^S)-dsqrt(csound2(ixO^S)))/ & 
                  (one-dsqrt(csound2(ixO^S))*vidim(ixO^S)),zero)
     endwhere
  else ! ndir
     where(Domain(ixO^S,5))

      cmax(ixO^S)=max((vidim(ixO^S)*(one-csound2(ixO^S))  &
                 +dsqrt(csound2(ixO^S)*(one-v2(ixO^S))*( &
                  one-v2(ixO^S)*csound2(ixO^S)-vidim(ixO^S)**2.0d0&
                 *(one-csound2(ixO^S)))))/(one-v2(ixO^S)*csound2(ixO^S)),zero)

      cmin(ixO^S)=min((vidim(ixO^S)*(one-csound2(ixO^S))  &
                 -dsqrt(csound2(ixO^S)*(one-v2(ixO^S))&
                 *(one-v2(ixO^S)*csound2(ixO^S)-vidim(ixO^S)**2.0d0&
                 *(one-csound2(ixO^S)))))/ (one-v2(ixO^S)*csound2(ixO^S)),zero)
     endwhere
  endif ! ndir
 end if needmin
endif  bbzero

end subroutine getcmax
!=============================================================================
subroutine getflux(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.
use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L,iw,idims
double precision, intent(in)       :: w(ixI^S,nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(out)      :: f(ixG^T)
logical, intent(out)               :: transport

double precision, dimension(ixG^T) :: VdotB, sqrB, ptot
!-----------------------------------------------------------------------------
transport=.true.

if (iw==d_) then
   f(ixO^S)=zero
{#IFDEF TRACER
{else if (iw==tr^FL_) then 
      f(ixO^S)=zero\}
}
{#IFDEF EPSINF
else if (iw==epsinf_) then 
      f(ixO^S)=zero
else if (iw==epslow_) then 
      f(ixO^S)=zero
else if (iw==n0_) then 
      f(ixO^S)=zero
else if (iw==n_) then 
      f(ixO^S)=zero
}
{#IFDEF GLM !f_i[psi]=Ch^2*b_{i}
! Eq. 24e Dedner et al 2002 JCP, 175, 645
else if(iw==psi_) then
{#IFNDEF FCT
       f(ixO^S)=cmax_global**2*w(ixO^S,b0_+idims)
}{#IFDEF FCT
       f(ixO^S)=zero
}
      transport=.false. 
}
else
   VdotB(ixO^S) = ({^C&w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/w(ixO^S,xi_)
   sqrB(ixO^S)  = {^C&w(ixO^S,b^C_)**2.0d0+}
   {if (iw==b^C_) then    !!! F^i[B_c] = B_c*v_i - v_c*B_i 
        if (idims==^C) then
           ! f_i[b_i] should be exactly 0, so we do not use the transport flux
           {#IFNDEF GLM f(ixO^S)=zero}
           {#IFDEF GLM   f(ixO^S)=w(ixO^S,psi_)}
           transport=.false.            
        else
           f(ixO^S)=-(w(ixO^S,s^C_)+VdotB(ixO^S)*w(ixO^S,b^C_))* &
                      w(ixO^S,b0_+idims)/(w(ixO^S,xi_)+sqrB(ixO^S))
        end if
     end if\}
    {if (iw==s^C_) then
        f(ixO^S)=-w(ixO^S,b0_+idims)*(VdotB(ixO^S)*&
                   (w(ixO^S,s^C_)+VdotB(ixO^S)*w(ixO^S,b^C_))/&
                   (w(ixO^S,xi_)+sqrB(ixO^S)) + &
                    w(ixO^S,b^C_)/w(ixO^S,lfac_)**2.0d0)
        if (idims==^C) then !!! add ptot
           call Pressure(w,ixI^L,ixO^L,.true.,ptot)
 
           ptot(ixO^S) = half*(VdotB(ixO^S)**2.0d0  & 
                       + sqrB(ixO^S)/w(ixO^S,lfac_)**2.0d0)  &
                       + ptot(ixO^S)
           f(ixO^S)=f(ixO^S) + ptot(ixO^S)
        end if
     end if\}
     if (iw==e_) then
        call Pressure(w,ixI^L,ixO^L,.true.,ptot)
        where(w(ixO^S,xi_)<smallxi)
           f(ixO^S)=zero
        elsewhere
          ptot(ixO^S) = half*(VdotB(ixO^S)**2.0d0  & 
                       + sqrB(ixO^S)/(w(ixO^S,lfac_)**2.0d0))  &
                       + ptot(ixO^S)
           ! f=ptot*v_i 
          f(ixO^S)=ptot(ixO^S)* &
                (w(ixO^S,s0_+idims)+VdotB(ixO^S)*w(ixO^S,b0_+idims))/ &
                (w(ixO^S,xi_)+sqrB(ixO^S))
           ! f=ptot*v_i - v.B*B_i
           f(ixO^S)=f(ixO^S) - VdotB(ixO^S)*w(ixO^S,b0_+idims)
        end where
     end if
  end if

end subroutine getflux
!=============================================================================
subroutine getfluxforhllc(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

use mod_global_parameters

integer, intent(in)::          ixI^L,ixO^L,iw,idims
double precision, intent(in)  :: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
double precision, intent(out) :: f(ixG^T,1:nwflux)
logical, intent(out)          :: transport

double precision, dimension(ixG^T):: VdotB, sqrB, ptot
!-----------------------------------------------------------------------------



transport=.true.
if (iw==d_) then
   f(ixO^S,iw)=zero
{#IFDEF TRACER
{else if (iw==tr^FL_) then 
      f(ixO^S,iw)=zero\}
}
{#IFDEF EPSINF
else if (iw==epsinf_) then 
      f(ixO^S,iw)=zero
else if (iw==epslow_) then 
      f(ixO^S,iw)=zero
else if (iw==n0_) then 
      f(ixO^S,iw)=zero
else if (iw==n_) then 
      f(ixO^S,iw)=zero
}
{#IFDEF GLM !f_i[psi]=Ch^2*b_{i}
! Eq. 24e Dedner et al 2002 JCP, 175, 645
else if(iw==psi_) then
{#IFNDEF FCT
       f(ixO^S,iw)=cmax_global**2*w(ixO^S,b0_+idims)
}{#IFDEF FCT
       f(ixO^S,iw)=zero
}
      transport=.false. 
}
else
   VdotB(ixO^S) = ({^C&w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/w(ixO^S,xi_)
   sqrB(ixO^S) = {^C&w(ixO^S,b^C_)**2.0d0+}
   {if (iw==b^C_) then    !!! F^i[B_c] = B_c*v_i - v_c*B_i 
        if (idims==^C) then
           ! f_i[b_i] should be exactly 0, so we do not use the transport flux
{#IFNDEF GLM f(ixO^S,iw)=zero }
{#IFDEF GLM   f(ixO^S,iw)=w(ixO^S,psi_) }
           transport=.false.            
        else
           f(ixO^S,iw)=-(w(ixO^S,s^C_)+VdotB(ixO^S)*w(ixO^S,b^C_))* &
                      w(ixO^S,b0_+idims)/(w(ixO^S,xi_)+sqrB(ixO^S))
        end if
     end if\}
    {if (iw==s^C_) then
        f(ixO^S,iw)=-w(ixO^S,b0_+idims)*(VdotB(ixO^S)*&
                   (w(ixO^S,s^C_)+VdotB(ixO^S)*w(ixO^S,b^C_))/&
                   (w(ixO^S,xi_)+sqrB(ixO^S)) + &
                    w(ixO^S,b^C_)/w(ixO^S,lfac_)**2.0d0)
        if (idims==^C) then !!! add ptot
           call Pressure(w,ixI^L,ixO^L,.true.,ptot)
           ptot(ixO^S) = half*(VdotB(ixO^S)**2.0d0  & 
                       + sqrB(ixO^S)/(w(ixO^S,lfac_)**2.0d0))  &
                       + ptot(ixO^S)
           f(ixO^S,iw)=f(ixO^S,iw) + ptot(ixO^S)
        end if
     end if\}
     if (iw==e_) then
        call Pressure(w,ixI^L,ixO^L,.true.,ptot)
        where(w(ixO^S,xi_)<smallxi)
           f(ixO^S,iw)=zero
        elsewhere
           ptot(ixO^S) = half*(VdotB(ixO^S)**2.0d0  & 
                       + sqrB(ixO^S)/(w(ixO^S,lfac_)**2.0d0))  &
                       + ptot(ixO^S)
           ! f=ptot*v_i 
           f(ixO^S,iw)=ptot(ixO^S)* &
                (w(ixO^S,s0_+idims)+VdotB(ixO^S)*w(ixO^S,b0_+idims))/ &
                (w(ixO^S,xi_)+sqrB(ixO^S))
           ! f=ptot*v_i - v.B*B_i
           f(ixO^S,iw)=f(ixO^S,iw) - VdotB(ixO^S)*w(ixO^S,b0_+idims)
        end where
     end if
  end if

end subroutine getfluxforhllc
!=============================================================================
subroutine addgeometry(qdt,ixI^L,ixO^L,wCT,w,x)

! Add geometrical source terms to w

use mod_global_parameters

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(in)       :: qdt
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(inout)    :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

integer                            :: iw,ix,idir, h1x^L{^NOONED, h2x^L}
double precision, dimension(ixG^T) :: tmp, sqrVdotB, sqrB, P
logical                            :: angmomfix=.false.
!-----------------------------------------------------------------------------

if(typeaxial /= 'slab')then
  call getaux(.true.,wCT,x,ixI^L,ixO^L,'addgeometry')
  sqrVdotB(ixO^S) = (({^C&wCT(ixO^S,s^C_)*wCT(ixO^S,b^C_)+})/wCT(ixO^S,xi_))**2
  sqrB(ixO^S) = {^C&wCT(ixO^S,b^C_)**2.0d0+}
endif

select case (typeaxial)
case ('slab')
   ! No source terms in slab symmetry
case ('cylindrical')
   do iw=1,nwflux
      select case (iw)
      ! source[sr]=(sphi*vphi-Bphi^2/lfac^2-VdotB*Bphi*vphi)/radius
      !             + ptot/radius
      case (sr_)
         !==== Calculate the pressure ====!
         call Pressure(wCT,ixI^L,ixO^L,.true.,P)
         tmp(ixO^S) = half*(sqrVdotB(ixO^S) + &
                sqrB(ixO^S)/(wCT(ixO^S,lfac_)**2)) + P(ixO^S)
         w(ixO^S,sr_)=w(ixO^S,sr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
         tmp(ixO^S)=zero
{^IFPHI
         ! now the rest.
         tmp(ixO^S)=(wCT(ixO^S,sphi_)**2-sqrVdotB(ixO^S)*wCT(ixO^S,bphi_)**2)/ &
              (wCT(ixO^S,xi_)+sqrB(ixO^S))  - &
              wCT(ixO^S,bphi_)**2/(wCT(ixO^S,lfac_)**2)
      ! source[sphi]=-(sphi*vr-Bphi*Br/lfac^2-VdotB*Bphi*vr)/radius
      case (sphi_)
         tmp(ixO^S)=-((wCT(ixO^S,sr_)*wCT(ixO^S,sphi_) - &
                    sqrVdotB(ixO^S)*wCT(ixO^S,br_)*wCT(ixO^S,bphi_))/ &
                    (wCT(ixO^S,xi_)+sqrB(ixO^S))  - &
                    wCT(ixO^S,br_)*wCT(ixO^S,bphi_)/(wCT(ixO^S,lfac_)**2) )
      case (bphi_)
          tmp(ixO^S)=-(wCT(ixO^S,sphi_)*wCT(ixO^S,br_) - &
                       wCT(ixO^S,sr_)*wCT(ixO^S,bphi_)) / &
                  (wCT(ixO^S,xi_)+sqrB(ixO^S))
}
{#IFDEF GLM      ! s[br]=psi/radius
      case (br_)
         tmp(ixO^S)=wCT(ixO^S,psi_)
}
      end select

      ! Divide by radius and add to w
      if (iw==sr_ {#IFDEF GLM .or.iw==br_} {^IFPHI .or.iw==sphi_.or.iw==bphi_}) then
           w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
      end if
   end do

   case('spherical')
      h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
      call Pressure(wCT,ixI^L,ixO^L,.true.,P)
      do iw=1,nwflux
         select case(iw)
         ! s[s1]=((Stheta**2+Sphi**2-(vdotB)**2*(Btheta**2+Bphi**2))/(xi+b2)
         !       +2*ptotal-(btheta**2+bphi**2)/lfac**2)/r
         case(s1_)
           ! total pressure
           tmp(ixO^S) = half*(sqrVdotB(ixO^S) + &
                sqrB(ixO^S)/wCT(ixO^S,lfac_)**2.0d0) + P(ixO^S)
           tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
                    *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
                    /mygeo%dvolume(ixO^S){&^CE&
            +(wCT(ixO^S,s^CE_)**2.0d0-sqrVdotB(ixO^S)*wCT(ixO^S,b^CE_)**2)/ &
             (wCT(ixO^S,xi_)+sqrB(ixO^S))-&
             wCT(ixO^S,b^CE_)**2.0d0/wCT(ixO^S,lfac_)**2.0 }
{^NOONEC
         ! s[s2]=-(Sr*Stheta-(vdotB)**2*Br*Btheta)/(xi+b2)/r+(Br*Btheta/lfac**2)/r
         !       + cot(theta)*(Sphi*Sphi-(vdotB)**2*Bphi*Bphi)/(xi+b2)/r
         !       + cot(theta)*(ptot- Bphi**2/lfac**2)/r
         case(s2_)
}
{^NOONED
           ! total pressure
           tmp(ixO^S) = half*(sqrVdotB(ixO^S) + &
                sqrB(ixO^S)/wCT(ixO^S,lfac_)**2.0d0) + P(ixO^S)
           w(ixO^S,iw)=w(ixO^S,iw) &
            +qdt*tmp(ixO^S)*(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
                      /mygeo%dvolume(ixO^S)
}
{^NOONEC
           tmp(ixO^S)=-((wCT(ixO^S,s1_)*wCT(ixO^S,s2_) - &
                         sqrVdotB(ixO^S)*wCT(ixO^S,b1_)*wCT(ixO^S,b2_))/ &
                       (wCT(ixO^S,xi_)+sqrB(ixO^S)) &
                        -wCT(ixO^S,b1_)*wCT(ixO^S,b2_)/wCT(ixO^S,lfac_)**2)
}
{^IFTHREEC
{^NOONED
            tmp(ixO^S)= tmp(ixO^S)+ ((wCT(ixO^S,s3_)**2.0d0&
                      - sqrVdotB(ixO^S)*wCT(ixO^S,b3_)**2.0d0)/ &
                        (wCT(ixO^S,xi_)+sqrB(ixO^S))&
                      - wCT(ixO^S,b3_)**2.0d0/wCT(ixO^S,lfac_)**2.0d0)&
                      * dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
}
         ! s[s3]=-(Sphi*Sr-(vdotB)**2*Bphi*Br)/(xi+b2)/r+(Bphi*Br/lfac**2)/r
         !       -cot(theta)*(Sphi*Stheta-(vdotB)**2*Bphi*Btheta)/(xi+b2)/r
         !       +cot(theta)*Bphi*Btheta/lfac**2/r
         case(s3_)
            if(.not.angmomfix) &
            tmp(ixO^S)=-((wCT(ixO^S,s3_)*wCT(ixO^S,s1_) - &
                         sqrVdotB(ixO^S)*wCT(ixO^S,b3_)*wCT(ixO^S,b1_))/ &
                       (wCT(ixO^S,xi_)+sqrB(ixO^S)) &
                        -wCT(ixO^S,b3_)*wCT(ixO^S,b1_)/wCT(ixO^S,lfac_)**2.0d0)
 {^NOONED 
            tmp(ixO^S)= tmp(ixO^S)-((wCT(ixO^S,s2_)*wCT(ixO^S,s3_) - &
                         sqrVdotB(ixO^S)*wCT(ixO^S,b2_)*wCT(ixO^S,b3_))/ &
                       (wCT(ixO^S,xi_)+sqrB(ixO^S)) &
                        -wCT(ixO^S,b2_)*wCT(ixO^S,b3_)/wCT(ixO^S,lfac_)**2.0d0) &
                      *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) 
}
}
{#IFDEF GLM
      ! s[b1]=2*psi/r
      case (b1_)
         tmp(ixO^S)=2.0d0*wCT(ixO^S,psi_)
}
{^NOONEC
         ! s[b2]=(sr*Btheta-stheta*Br)/(xi+B^2)/r
         !       + cot(theta)*psi/r
         case(b2_)

            tmp(ixO^S)=(wCT(ixO^S,s1_)*wCT(ixO^S,b2_) &
                       -wCT(ixO^S,s2_)*wCT(ixO^S,b1_))&
                       /(wCT(ixO^S,xi_)+sqrB(ixO^S))
}
{#IFDEF GLM
{^NOONED
            tmp(ixO^S)=tmp(ixO^S) &
                 + dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) &
                 * wCT(ixO^S,psi_)
}
}

{^IFTHREEC
         ! s[b3]=(sr*Bphi-sphi*Br)/(xi+B^2)/r
         !       -cot(theta)*(sphi*Btheta-stheta*Bphi)/(xi+B^2)/r
         case(b3_)
            tmp(ixO^S)=(wCT(ixO^S,s1_)*wCT(ixO^S,b3_) &
                   -wCT(ixO^S,s3_)*wCT(ixO^S,b1_))&
                   /(wCT(ixO^S,xi_)+sqrB(ixO^S)) {^NOONED &
                   -(wCT(ixO^S,s3_)*wCT(ixO^S,b2_) &
                   -wCT(ixO^S,s2_)*wCT(ixO^S,b3_))*dcos(x(ixO^S,2)) &
                   /((wCT(ixO^S,xi_)+sqrB(ixO^S))*dsin(x(ixO^S,2))) }
}
         end select
         ! Divide by radius and add to wnew
         if(iw==s1_{#IFDEF GLM .or.iw==b1_}{^NOONEC.or.iw==s2_.or.iw==b2_}&
              {^IFTHREEC .or.iw==b3_ .or.(iw==s3_.and..not.angmomfix)}) &
            w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
      end do
end select

end subroutine addgeometry
!=============================================================================
subroutine addsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(in)    :: wCT(ixI^S,1:nw)
double precision, intent(inout) :: w(ixI^S,1:nw)
logical, intent(in)             :: qsourcesplit

double precision :: dx^D
!-----------------------------------------------------------------------------

dx^D=dxlevel(^D);

! Sources related to div B
{^NOONED
if(qsourcesplit .eqv. ssplitdivb) then
select case (typedivbfix)
{#IFDEF GLM
case ('glm1')
   call addsource_glm1(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
case ('glm2')
   call addsource_glm2(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)\}
case ('powel')
   call addsource_powel(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
case ('janhunen')
   call addsource_janhunen(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
case ('linde')
   call addsource_linde(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
end select

! now update the auxiliaries to the new state
if (typedivbfix=='powel' .or. typedivbfix=='janhunen' .or. typedivbfix=='linde' .or. typedivbfix=='glm2') &
call getaux(.true.,w,x,ixI^L,ixO^L,'addsource')
endif
}

end subroutine addsource
!=============================================================================
subroutine getcurrent(w,ixI^L,ix^L,idirmin,current)

! Calculate idirmin and the idirmin:3 components of the common current array

use mod_global_parameters

integer, parameter:: idirmin0=7-2*ndir
integer :: ix^L, idirmin, ixI^L
double precision :: w(ixI^S,1:nw)

! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
double precision :: current(ixG^T,7-2*ndir:3),bvec(ixG^T,1:ndir)
!-----------------------------------------------------------------------------

bvec(ixI^S,1:ndir)=w(ixI^S,b0_+1:b0_+ndir)
call curlvector(bvec,ixI^L,ix^L,current,idirmin,idirmin0,ndir)

end subroutine getcurrent
!=============================================================================
subroutine laguer(a,m,x)
! slightly adapted from Num. recipes
implicit none

integer:: m
double complex:: a(1:m+1),x

integer:: maxit,mr,mt
double precision:: epss
parameter (epss=2.0d-15,mr=8,mt=10,maxit=mt*mr)

integer:: iter,j
double precision:: abx,abp,abm,err,frac(1:mr)
double complex:: dx,x1,b,d,f,g,h,sq,gp,gm,g2
!-----------------------------------------------------------------------------

frac(1)=0.5d0
frac(2)=0.25d0
frac(3)=0.75d0
frac(4)=0.13d0
frac(5)=0.38d0
frac(6)=0.62d0
frac(7)=0.88d0
frac(8)=1.0d0

do iter=1,maxit
   b=a(m+1)
   err=cdabs(b)
   d=dcmplx(0.0d0,0.0d0)
   f=dcmplx(0.0d0,0.0d0)
   abx=cdabs(x)
   do j=m,1,-1
      f=x*f+d
      d=x*d+b
      b=x*b+a(j)
      err=cdabs(b)+abx*err
   end do
   err=epss*err
   if(cdabs(b)<=err) then
      return
   else
      g=d/b
      g2=g*g
      h=g2-2.0d0*f/b
      sq=cdsqrt((dble(m)-1.0d0)*(dble(m)*h-g2))
      gp=g+sq
      gm=g-sq
      abp=cdabs(gp)
      abm=cdabs(gm)
      if(abp<abm) gp=gm
      if (max(abp,abm)>0.0d0) then
         dx=dble(m)/gp
      else
         dx=cdexp(dcmplx(dlog(1.0d0+abx),dble(iter)))
      endif
   endif
   x1=x-dx
   if(x==x1)return
   if (mod(iter,mt)/=0) then
      x=x1
   else
      x=x-dx*frac(iter/mt)
   endif
end do

!print *,'WARNING: too many iterations in laguer'
!print *,'         coefficients:'
!print *,m,a(1:m+1)
call mpistop("laguer problem")

return
end subroutine laguer
!=============================================================================
subroutine zroots4(a,roots)
! adapted version from Num. recipes (no sorting)
use mod_global_parameters

  
double precision,parameter :: eps=1.0d-14
double precision,parameter :: accroot=1.0d-8
double precision,parameter :: relaccroot=1.0d-8
integer,parameter          :: maxm=5
integer,parameter          :: m=4
integer,parameter          :: NRcycles=100

double precision:: a(1:maxm),roots(1:m)

integer :: iroot,j,jj
double precision:: p1,p,root,valuepoly,root1,valuepoly1
double complex:: croots(1:m),ad(1:maxm),x,b,c
!-----------------------------------------------------------------------------

ad(1:m+1)=dcmplx(a(1:m+1),0.0d0)

do j=m,1,-1
   x=dcmplx(0.0d0,0.0d0)
   call laguer(ad,j,x)
   if (dabs(dimag(x)) <= 2.0d0*eps**2*dabs(dble(x))) x=dcmplx(dble(x),0.0d0)
   croots(j)=x
   b=ad(j+1)
   do jj=j,1,-1
      c=ad(jj)
      ad(jj)=b
      b=x*b+c
   end do
end do

roots(1:m) = dble(croots(1:m))

! polish the real roots using Newton-raphson
do iroot=1,m
    root=roots(iroot)
    root1=root
    valuepoly=a(1)+root*(a(2)+root*(a(3)+root*(a(4)+root*a(5))))
    valuepoly1=valuepoly
    !if(dabs(valuepoly)>accroot) then
    !   print *,' root =',root,' value polynomial=',valuepoly
    !endif
    ! next test only for bergmans rescaling, inactivated
    if(typepoly=='bergmans' .and. root<0.5d0)then
        print *,'error in rootfinder'
        print *,roots(1:m),iroot
        call mpistop("zroots problem")
    endif
    if(dabs(valuepoly)>accroot)then
      !print *,'above ',accroot,' doing NR'
      do j=1,NRcycles
         p  = a(m+1)*roots(iroot) + a(m)
         p1 = a(m+1)
         do jj=m-1,1,-1
            p1 = p + p1*roots(iroot)
            p = a(jj) + p*roots(iroot)
         end do
         if (dabs(p)<accroot) exit
         if (p1==0.0d0) then
            print *,'error in NR'
            call mpistop("zroots error")
         end if
         roots(iroot)=roots(iroot) - p/p1
         if (dabs(p/p1/roots(iroot))<relaccroot) exit
         !print *,'NR cycle ',j,' root=',roots(iroot), &
         !    ' p=',p,' p1=',p1,' p/p1=',p/p1
      end do
      root=roots(iroot)
      valuepoly=a(1)+root*(a(2)+root*(a(3)+root*(a(4)+root*a(5))))
      if(typepoly=='bergmans' .and. root<0.5d0)then
        print *,'error in rootfinder AFTER NR'
        print *,roots(1:m),iroot
        call mpistop("zroots problem")
      endif
      if(dabs(valuepoly)>accroot) then
        print *,'AFTER NR root =',root,' value polynomial=',valuepoly
        print *,'before NR root:',root1,' value=',valuepoly1
        !call mpistop("")
      endif
    endif
enddo

end subroutine zroots4
!=============================================================================
! end module amrvacphys/- srmhdeos
!#############################################################################
