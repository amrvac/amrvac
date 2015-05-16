!#############################################################################
! module amrvacphys/- srmhd  version april 2009

INCLUDE:amrvacnul/getdt.t
INCLUDE:amrvacnul/roe.t
!=============================================================================
subroutine checkglobaldata
include 'amrvacdef.f'
!-----------------------------------------------------------------------------
govergminone =eqpar(gamma_)/(eqpar(gamma_)-one)
minp   = max(zero,smallp)
minrho = max(zero,smallrho)
smallxi=minrho+minp*govergminone
smalltau = minp/(eqpar(gamma_)-one)

end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

! place to set entropy fixes etc, absent for now

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
eqpar(gamma_)=4./3.

{#IFDEF GLM
eqpar(Cr_)=0.18
}

if(strictzero)limitvalue=zero
if(.not.strictzero)limitvalue=smalldouble**2.0d0

end subroutine initglobaldata
!=============================================================================
subroutine checkw(checkprimitive,ixI^L,ixO^L,w,flag)

include 'amrvacdef.f'
  
logical, intent(in)          :: checkprimitive
integer, intent(in)          :: ixI^L,ixO^L
double precision, intent(in) :: w(ixI^S,1:nw)
logical, intent(out)         :: flag(ixG^T)
!-----------------------------------------------------------------------------

flag(ixG^T)=.true.

if (checkprimitive) then
  if(useprimitiveRel)then
     ! check   rho>=0, p>=smallp
     flag(ixO^S) = (w(ixO^S,rho_) > minrho).and. &
                   (w(ixO^S,pp_)  >=minp&
{#IFDEF EPSINF         .and. w(ixO^S,rho1_) > minrho})
               
  else
    ! check  v^2 < 1, rho>=0, p>=minp
     ! v^2 < 1, rho>0, p>0
     flag(ixO^S) = ({^C&w(ixO^S,v^C_)**2.0d0+}< one).and. &
                   (w(ixO^S,rho_) > minrho).and. &
                   (w(ixO^S,pp_) >=minp&
{#IFDEF EPSINF         .and. w(ixO^S,rho1_) > minrho&
                    .and. w(ixO^S,n_) > minrho})
  endif
else
  ! checks on the conservative variables
  flag(ixO^S)= (w(ixO^S,d_)   > minrho).and. &
               (w(ixO^S,tau_) > smalltau&
{#IFDEF EPSINF         .and. w(ixO^S,Drho1_) > minrho&
                    .and. w(ixO^S,Dn_) > minrho})
end if
  
!if(any(.not.flag(ixO^S))) call mpistop("bad check found")

end subroutine checkw
!=============================================================================
subroutine conserve(ixI^L,ixO^L,w,x,patchw)

! Transform primitive variables into conservative ones
! (rho,v,p,B,tr) ---> (D,S,tau,B,Dtr,lfac,xi)
! call to smallvalues
! --> latter only used for correcting procedure in correctaux
! --> input array patchw for spatially selective transformation
! opedit: Added tracer scalars tr that are advected.  Set nf>1.  
! I do that by solving the equivalent conservation law for D*tr
! It is d_t (D tr) + div (D tr \vec{v})= 0
! <=>   d_t (tr) + \vec{v} div (tr)= 0
! since d_t (D) + div (D \vec{v})= 0.  

include 'amrvacdef.f'

integer, intent(in)               :: ixI^L,ixO^L
double precision, intent(inout)   :: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
logical, intent(in)               :: patchw(ixG^T)

double precision, dimension(ixG^T):: sqrV,sqrU,sqrB,VdotB,rhoh
!-----------------------------------------------------------------------------

if(useprimitiveRel)then
  ! assumes four velocity computed in primitive (rho u p B) with u=lfac*v
  where(.not.patchw(ixO^S))
     sqrU(ixO^S) ={^C&w(ixO^S,u^C_)**2.0d0+} 
     sqrV(ixO^S)   =sqrU(ixO^S)/(one+sqrU(ixO^S))
     sqrB(ixO^S) ={^C&w(ixO^S,b^C_)**2.0d0+}
  endwhere
else
  ! assumes velocity in primitive (rho v p B) 
  where(.not.patchw(ixO^S))
     sqrV(ixO^S) ={^C&w(ixO^S,v^C_)**2.0d0+}
     sqrB(ixO^S) ={^C&w(ixO^S,b^C_)**2.0d0+} 
  endwhere
endif

! fill the auxiliary variable lfac (lorentz factor)
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

! fill the auxiliary variable xi and density D
! with enthalpy w: xi= lfac^2 rhoh
! density: d = lfac * rho
where(.not.patchw(ixO^S))
   w(ixO^S,xi_)=w(ixO^S,lfac_)*w(ixO^S,lfac_) &
                 *(w(ixO^S,rho_)+govergminone*w(ixO^S,pp_)) 
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
   w(ixO^S,Drho1_) =w(ixO^S,rho1_)*w(ixO^S,lfac_)
   w(ixO^S,Drho0_) =w(ixO^S,Drho1_)*w(ixO^S,rho0_)
   w(ixO^S,Dn_) =w(ixO^S,n_)*w(ixO^S,lfac_)
   w(ixO^S,Dn0_) =w(ixO^S,Dn_)*w(ixO^S,n0_)
   w(ixO^S,Depsinf_) = w(ixO^S,epsinf_)*w(ixO^S,Drho1_)**(2.0d0/3.0d0) &
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

! E = xi - p +B^2/2 + (v^2 B^2 - (v.B)^2)/2 
! instead of E use tau= E - D
where(.not.patchw(ixO^S))
    w(ixO^S,tau_)=w(ixO^S,xi_) - w(ixO^S,pp_) - w(ixO^S,d_) +& 
         half*(sqrB(ixO^S) + sqrU(ixO^S))
endwhere

if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,patchw,'conserve')

end subroutine conserve
!=============================================================================
subroutine conserven(ixI^L,ixO^L,w,patchw)

! Transform primitive variables into conservative ones
! (rho,v,p,B) ---> (D,S,tau,B,lfac,xi)
! no call to smallvalues
! --> latter only used for correcting procedure in correctaux
! --> input array patchw for spatially selective transformation

include 'amrvacdef.f'

integer, intent(in)               :: ixI^L,ixO^L
double precision, intent(inout)   :: w(ixI^S,nw)
logical, intent(in)               :: patchw(ixG^T)

double precision, dimension(ixG^T):: sqrV,sqrU,sqrB,VdotB,rhoh
!-----------------------------------------------------------------------------

if(useprimitiveRel)then
  ! assumes four velocity computed in primitive (rho u p B) with u=lfac*v
  where(.not.patchw(ixO^S))
     sqrU(ixO^S) ={^C&w(ixO^S,u^C_)**2.0d0+} 
     sqrV(ixO^S)   =sqrU(ixO^S)/(one+sqrU(ixO^S))
     sqrB(ixO^S) ={^C&w(ixO^S,b^C_)**2.0d0+}
  endwhere
else
  ! assumes velocity in primitive (rho v p B) 
  where(.not.patchw(ixO^S))
     sqrV(ixO^S) ={^C&w(ixO^S,v^C_)**2.0d0+}
     sqrB(ixO^S) ={^C&w(ixO^S,b^C_)**2.0d0+} 
  endwhere
endif

! fill the auxiliary variable lfac (lorentz factor)
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

! fill the auxiliary variable xi and density D
! with enthalpy w: xi= lfac^2 rhoh
! density: d = lfac * rho
where(.not.patchw(ixO^S))
   w(ixO^S,xi_)=w(ixO^S,lfac_)*w(ixO^S,lfac_) &
                 *(w(ixO^S,rho_)+govergminone*w(ixO^S,pp_)) 
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
   w(ixO^S,Drho1_) =w(ixO^S,rho1_)*w(ixO^S,lfac_)
   w(ixO^S,Drho0_) =w(ixO^S,Drho1_)*w(ixO^S,rho0_)
   w(ixO^S,Dn_) =w(ixO^S,n_)*w(ixO^S,lfac_)
   w(ixO^S,Dn0_) =w(ixO^S,Dn_)*w(ixO^S,n0_)
   w(ixO^S,Depsinf_) = w(ixO^S,epsinf_)*w(ixO^S,Drho1_)**(2.0d0/3.0d0) &
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

! E = xi - p +B^2/2 + (v^2 B^2 - (v.B)^2)/2 
! instead of E use tau= E - D
where(.not.patchw(ixO^S))
    w(ixO^S,tau_)=w(ixO^S,xi_) - w(ixO^S,pp_) - w(ixO^S,d_) +& 
         half*(sqrB(ixO^S) + sqrU(ixO^S))
endwhere

end subroutine conserven
!=============================================================================
subroutine primitive(ixI^L,ixO^L,w,x)

! Transform conservative variables into primitive ones
! (D,S,tau,B)-->(rho,v,p,B,lfac,xi)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(inout)    :: w(ixI^S,1:nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, dimension(ixG^T) :: SdotB, sqrB
!-----------------------------------------------------------------------------

! calculate lorentz factor and xi from conservatives only
call getaux(.true.,w,x,ixI^L,ixO^L,'primitive')

w(ixO^S,rho_)=w(ixO^S,d_)/w(ixO^S,lfac_)
w(ixO^S,pp_)=(w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0-w(ixO^S,rho_))/govergminone

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
w(ixO^S,rho1_)   = w(ixO^S,Drho1_)/w(ixO^S,lfac_)
w(ixO^S,rho0_)   = w(ixO^S,Drho0_)/w(ixO^S,lfac_)/w(ixO^S,rho1_)
w(ixO^S,n_)   = w(ixO^S,Dn_)/w(ixO^S,lfac_)
w(ixO^S,n0_)   = w(ixO^S,Dn0_)/w(ixO^S,lfac_)/w(ixO^S,n_)
w(ixO^S,epsinf_) = w(ixO^S,Depsinf_)/(w(ixO^S,lfac_) &
     *w(ixO^S,rho1_)**(2.0d0/3.0d0))
}
if (tlow>zero) call fixp_usr(ixI^L,ixO^L,w,x)

end subroutine primitive
!==============================================================================
subroutine primitiven(ixI^L,ixO^L,w,patchw)

! same as primitive, but with padding, no getaux call
! --> needed in correctaux, assumes updated xi and lfac, no smallxi fix
! Transform conservative variables into primitive ones
! (D,S,tau,B)-->(rho,v,p,B,lfac,xi)

include 'amrvacdef.f'

integer, intent(in)                    :: ixI^L,ixO^L
double precision, intent(inout)        :: w(ixI^S,1:nw)
logical, intent(in),dimension(ixG^T)   :: patchw

double precision, dimension(ixG^T)     :: SdotB, sqrB
!-----------------------------------------------------------------------------

where(.not.patchw(ixO^S))
   w(ixO^S,rho_)=w(ixO^S,d_)/w(ixO^S,lfac_)
   w(ixO^S,pp_)=(w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0-w(ixO^S,rho_))/govergminone
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
  end where
endif

{#IFDEF TRACER
! We got lor, rho, Dtr, now we can get the tracers:
where(.not.patchw(ixO^S))
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,Dtr^FL_)/w(ixO^S,lfac_)/w(ixO^S,rho_)\}
endwhere
}
{#IFDEF EPSINF
where(.not.patchw(ixO^S))
w(ixO^S,rho1_)   = w(ixO^S,Drho1_)/w(ixO^S,lfac_)
w(ixO^S,rho0_)   = w(ixO^S,Drho0_)/w(ixO^S,lfac_)/w(ixO^S,rho1_)
w(ixO^S,n_)      = w(ixO^S,Dn_)/w(ixO^S,lfac_)
w(ixO^S,n0_)     = w(ixO^S,Dn0_)/w(ixO^S,lfac_)/w(ixO^S,n_)
w(ixO^S,epsinf_) = w(ixO^S,Depsinf_)/(w(ixO^S,lfac_) &
     *w(ixO^S,rho1_)**(2.0d0/3.0d0))
end where
}

end subroutine primitiven
!=============================================================================
subroutine e_to_rhos(ixI^L,ixO^L,w,x)

include 'amrvacdef.f'
integer:: ixI^L,ixO^L
double precision:: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-------------------------------------------------------------------------

call mpistop("e to rhos for SRMHD unavailable")

end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixI^L,ixO^L,w,x)

include 'amrvacdef.f'

integer:: ixI^L,ixO^L
double precision:: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

call mpistop("e to rhos for SRMHD unavailable")

end subroutine rhos_to_e
!=============================================================================
subroutine ppmflatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L,ixL^L,ixR^L
double precision, intent(in)  :: w(ixI^S,nw),d2w(ixG^T,1:nwflux)

double precision, intent(out) :: drho(ixG^T),dp(ixG^T)

double precision     :: ptot(ixG^T)
!-----------------------------------------------------------------------------

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
   drho(ixO^S) =eqpar(gamma_)*dabs(d2w(ixO^S,rho_))&
              /min(w(ixL^S,rho_),w(ixR^S,rho_))
   dp(ixO^S) = dabs(ptot(ixR^S)-ptot(ixL^S))/min(ptot(ixL^S),ptot(ixR^S))
end if

end subroutine ppmflatcd
!=============================================================================
subroutine ppmflatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
integer, intent(in)           :: idims
double precision, intent(in)  :: w(ixI^S,nw)

double precision, intent(out) :: drho(ixG^T),dp(ixG^T),dv(ixG^T)

double precision :: v(ixG^T),ptot(ixG^T)
!-----------------------------------------------------------------------------

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
   where (dabs(ptot(ixRR^S)-ptot(ixLL^S))>smalldouble)
      drho(ixO^S) = dabs((ptot(ixR^S)-ptot(ixL^S))&
                   /(ptot(ixRR^S)-ptot(ixLL^S)))
   elsewhere
      drho(ixO^S) = zero
   end where

   !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26 
   !  use "dp" to save sound speed squared
   dp(ixO^S) = (eqpar(gamma_)-one) &
             *(w(ixO^S,xi_)-w(ixO^S,rho_)*w(ixO^S,lfac_)**2.0d0)/w(ixO^S,xi_)
   dp(ixO^S) = dabs(ptot(ixR^S)-ptot(ixL^S))&
                /(w(ixO^S,rho_)*dp(ixO^S))
   if (useprimitiveRel) then
     v(ixI^S)  = w(ixI^S,u0_+idims)/w(ixI^S,lfac_)
   else
     v(ixI^S)  = w(ixI^S,v0_+idims)
   end if
   call gradient(v,ixO^L,idims,dv)
end if

end subroutine ppmflatsh
!=============================================================================
subroutine getv(w,x,ixI^L,ixO^L,idims,v)

! Calculate v_idim=m_idim/rho within ixO^L

include 'amrvacdef.f'
  
integer, intent(in)                              :: ixI^L,ixO^L,idims
double precision, intent(in)                     :: w(ixI^S,1:nw)
double precision, intent(in)                     :: x(ixI^S,1:ndim)
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
  
include 'amrvacdef.f'
  
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
        
      dedp =  (w(ix^D,xi_)/(w(ix^D,xi_) - w(ix^D,lfac_)*w(ix^D,d_))) &
               /(eqpar(gamma_)-one)

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
        lambdamax = maxval(dabs(lambdaroots)) 
        cmax(ix^D)=min(one,lambdamax)
      end if
 
      !if (lambdamax>=one) then 
         !print*,'WARNING: getcmax found eigenvalue >= 1 !'
         !print*,ix^D
         !print*,'roots in mu: ',roots
         !do iroot=1,order
         !  root=roots(iroot)
         !  print*,iroot,root,poly(0)+poly(1)*root+poly(2)*root**2&
         !         +poly(3)*root**3+poly(4)*root**4
         !enddo
         !print*,'eigenvalues: ',lambdaroots
         !print*,'lambdamax=',lambdamax
         !write(*,'(a1,5(f20.12,a1))') &
         !    '{',poly(0),',',poly(1),',',poly(2),',',poly(3),',',poly(4),'}'
         !print*,'d=',w(ix^D,d_),'s=',{^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
         !print*,'B=',{^C&w(ix^D,b^C_)},'lfac=',w(ix^D,lfac_),'xi=',w(ix^D,xi_)
      !   call mpistop("== getcmax ==")
      !end if
        

    endif ! Domain 1
  {enddo^D&\}
 case('bergmans')
  {do ix^D= ixO^LIM^D\}
   if(Domain(ix^D,1))then
      VdotB = ({^C&w(ix^D,s^C_)*w(ix^D,b^C_)+})/w(ix^D,xi_)
      sqrlfac = w(ix^D,lfac_)**2.0d0
      dedp =  (w(ix^D,xi_)/(w(ix^D,xi_) - w(ix^D,lfac_)*w(ix^D,d_))) &
               /(eqpar(gamma_)-one)
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
 
      !if (lambdamax>=one) then 
         !print*,'WARNING: getcmax found eigenvalue >= 1 !'
         !print*,ix^D
         !print*,'roots in mu: ',roots
         !do iroot=1,order
         !  root=roots(iroot)
         !  print*,iroot,root,poly(0)+poly(1)*root+poly(2)*root**2&
         !         +poly(3)*root**3+poly(4)*root**4
         !enddo
         !print*,'eigenvalues: ',lambdaroots
         !print*,'lambdamax=',lambdamax
         !write(*,'(a1,5(f20.12,a1))') &
         !     '{',poly(0),',',poly(1),',',poly(2),',',poly(3),',',poly(4),'}'
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
    dedp =  (w(ix^D,xi_)/(w(ix^D,xi_) - w(ix^D,lfac_)*w(ix^D,d_))) &
               /(eqpar(gamma_)-one)
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

    !if (lambdamax>=one) then 
         !print*,'WARNING: getcmax found eigenvalue >= 1 !'
         !print*,ix^D
         !print*,'roots in mu: ',roots
         !do iroot=1,order
         !  root=roots(iroot)
         !  print*,iroot,root,poly(0)+poly(1)*root+poly(2)*root**2&
         !         +poly(3)*root**3+poly(4)*root**4
         !enddo
         !print*,'eigenvalues: ',lambdaroots
         !print*,'lambdamax=',lambdamax
         !write(*,'(a1,5(f20.12,a1))') &
         !     '{',poly(0),',',poly(1),',',poly(2),',',poly(3),',',poly(4),'}'
         !print*,'d=',w(ix^D,d_),'s=',{^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
         !print*,'B=',{^C&w(ix^D,b^C_)},'lfac=',w(ix^D,lfac_),'xi=',w(ix^D,xi_)
    !     call mpistop("== getcmax ==")
    !end if

   endif ! Domain 1 meliani
  {enddo^D&\}
 case('gammie')

   where(Domain(ixO^S,1))
     rhoh(ixO^S)   = w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
     ! squared Alfven speed
     B(ixO^S)  = sqrB(ixO^S)/(sqrB(ixO^S)+rhoh(ixO^S))
     ! squared sound speed
     csound2(ixO^S)= (eqpar(gamma_)-one) &
              *(w(ixO^S,xi_)-w(ixO^S,d_)*w(ixO^S,lfac_))/w(ixO^S,xi_)
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
  csound2(ixO^S)=(eqpar(gamma_)-one)*(rhoh(ixO^S)-w(ixO^S,d_))/rhoh(ixO^S)

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
   csound2(ixO^S)= (eqpar(gamma_)-one)*(w(ixO^S,xi_)-w(ixO^S,d_)*w(ixO^S,lfac_))/w(ixO^S,xi_)
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
   csound2(ixO^S)=(eqpar(gamma_)-one)*(w(ixO^S,xi_)-w(ixO^S,d_)*w(ixO^S,lfac_))/w(ixO^S,xi_)
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
  csound2(ixO^S)= (eqpar(gamma_)-one)*(w(ixO^S,xi_)-w(ixO^S,d_))/w(ixO^S,xi_)

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
   csound2(ixO^S) = (eqpar(gamma_)-one)*(w(ixO^S,xi_)-w(ixO^S,d_)*w(ixO^S,lfac_))/w(ixO^S,xi_) 
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
    csound2(ixO^S) = (eqpar(gamma_)-one)*(w(ixO^S,xi_)-w(ixO^S,d_)*w(ixO^S,lfac_))/w(ixO^S,xi_) 
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
!--------------------------------------------------------------------------!

end subroutine getcmax
!=============================================================================
subroutine getflux(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

include 'amrvacdef.f'

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
else if (iw==rho0_) then 
      f(ixO^S)=zero
else if (iw==rho1_) then 
      f(ixO^S)=zero
else if (iw==n0_) then 
      f(ixO^S)=zero
else if (iw==n_) then 
      f(ixO^S)=zero
}
{#IFDEF GLM !f_i[psi]=Ch^2*b_{i}
   ! Eq. 24e Dedner et al 2002 JCP, 175, 645
else if (iw==psi_) then
     f(ixO^S)=storeeqparch**2.0d0*w(ixO^S,b0_+idims)
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
           ptot(ixO^S) = half*(VdotB(ixO^S)**2.0d0  & 
                       + sqrB(ixO^S)/w(ixO^S,lfac_)**2.0d0)  &
                       + (eqpar(gamma_)-one) &
                       * (w(ixO^S,xi_)-w(ixO^S,lfac_)*w(ixO^S,d_)) &
                       / (eqpar(gamma_)*w(ixO^S,lfac_)**2.0d0)
           f(ixO^S)=f(ixO^S) + ptot(ixO^S)
        end if
     end if\}
     if (iw==e_) then
        where(w(ixO^S,xi_)<smallxi)
           f(ixO^S)=zero
        elsewhere
          ptot(ixO^S) = half*(VdotB(ixO^S)**2.0d0  & 
                       + sqrB(ixO^S)/(w(ixO^S,lfac_)**2.0d0))  &
                       + (eqpar(gamma_)-one) &
                       * (w(ixO^S,xi_)-w(ixO^S,lfac_)*w(ixO^S,d_)) &
                       / (eqpar(gamma_)*w(ixO^S,lfac_)**2.0d0)
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

include 'amrvacdef.f'

integer, intent(in)::          ixI^L,ixO^L,iw,idims
double precision, intent(in)  :: w(ixI^S,nw)
double precision, intent(in)  :: x(ixI^S,1:ndim)
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
else if (iw==rho0_) then 
      f(ixO^S,iw)=zero
else if (iw==rho1_) then 
      f(ixO^S,iw)=zero
else if (iw==n0_) then 
      f(ixO^S,iw)=zero
else if (iw==n_) then 
      f(ixO^S,iw)=zero
}
{#IFDEF GLM !f_i[psi]=Ch^2*b_{i}
! Eq. 24e Dedner et al 2002 JCP, 175, 645
else if(iw==psi_) then
      f(ixO^S,iw)=storeeqparch**2.0d0*w(ixO^S,b0_+idims)
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
           ptot(ixO^S) = half*(VdotB(ixO^S)**2.0d0  & 
                       + sqrB(ixO^S)/(w(ixO^S,lfac_)**2.0d0))  &
                       + (eqpar(gamma_)-one) &
                       * (w(ixO^S,xi_)-w(ixO^S,lfac_)*w(ixO^S,d_)) &
                       / (eqpar(gamma_)*w(ixO^S,lfac_)**2.0d0)
           f(ixO^S,iw)=f(ixO^S,iw) + ptot(ixO^S)
        end if
     end if\}
     if (iw==e_) then
        where(w(ixO^S,xi_)<smallxi)
           f(ixO^S,iw)=zero
        elsewhere
           ptot(ixO^S) = half*(VdotB(ixO^S)**2.0d0  & 
                       + sqrB(ixO^S)/(w(ixO^S,lfac_)**2.0d0))  &
                       + (eqpar(gamma_)-one) &
                       * (w(ixO^S,xi_)-w(ixO^S,lfac_)*w(ixO^S,d_)) &
                       / (eqpar(gamma_)*w(ixO^S,lfac_)**2.0d0)
           ! f=ptot*v_i 
           f(ixO^S,iw)=ptot(ixO^S)* &
                (w(ixO^S,s0_+idims)+VdotB(ixO^S)*w(ixO^S,b0_+idims))/ &
                (w(ixO^S,xi_)+sqrB(ixO^S))
           ! f=ptot*v_i - v.B*B_i
           f(ixO^S,iw)=f(ixO^S,iw) - VdotB(ixO^S)*w(ixO^S,b0_+idims)
        end where
     end if ! iw==e_
  end if

end subroutine getfluxforhllc
!=============================================================================
subroutine con2prim(lfac,xi,d,s^C,tau,b^C,ierror)

! (D,S,tau,B) --> compute auxiliaries lfac and xi

include 'amrvacdef.f'

double precision:: lfac,xi
double precision:: d,s^C,tau,b^C
integer         :: ierror
      
double precision:: sdotb,ssqr,bsqr
double precision:: xi1,xi2,xii,xih,xil,dxi,f,df
double precision:: temp,vsqr,fl,fh,lfacl,lfach,lfaci
double precision:: er,er1,xiprev,dplus,lfacmax
double precision:: fv,dfv,v^C,pressure
double precision:: lastxiok,lastlfacok,lastf
integer:: ni,nit,niiter

logical:: finished,tests,testsl,testsh
!-----------------------------------------------------------------------------

ierror=0

! ierror=1 : error on entry: must have D>=0, tau>=smallp/(gamma-1)
! ierror=2 : maxitnr reached without convergence
! ierror=3 : no range for solution was find
! ierror=4 : final xi <smallxi
! ierror=5 : lower bound non-negative f(xi) value
! ierror=6 : v^2>1
! ierror=7 :  solution out from the range?
! ierror=8 : nonmonotonic function f?

if(d<minrho .or. tau<smalltau)then
  ierror=1
  return
endif

bsqr={b^C**2.0d0+}
ssqr={s^C**2.0d0+}

! handle the case of no flow first
if(ssqr<=limitvalue)then
  xi = eqpar(gamma_) *(tau - half * bsqr)+d
  lfac=one
  return
endif

! Hydro case: handle exactly as before (using p to iterate on)
if(bsqr<=limitvalue)then
  call con2primHydro(lfac,xi,d,s^C,tau,ierror)
  return 
end if ! Hydro

temp={s^C*b^C+}
if(temp<zero)then
  sdotb=-min(dsqrt(ssqr*bsqr),-temp)
else
  sdotb=min(dsqrt(ssqr*bsqr),temp)
endif

! the maximal allowed Lorentz factor
lfacmax=one/dsqrt(one-(one-dmaxvel)**2.0d0)
! starting point for NR on xi
! xi1 is meant to be the lower bound for the bracket to be used in NR
dplus=d+eqpar(gamma_)*minp/(eqpar(gamma_)-one)
xi1=dplus

! compute v^2(xi1)
vsqr = ((ssqr + (xi1*two*sdotb*sdotb + bsqr*sdotb*sdotb)/(xi1*xi1))&
      /((xi1+bsqr)*(xi1+bsqr)))

!=== find new xi1, in the case v2(xi1) > maxvel^2 = (1-dmaxvel)^2 
! locate xi corresponding to maximal velocity allowed, namely 1-dmaxvel
niiter=-1
if(vsqr > (one-dmaxvel)**2.0d0 )then 
   er1=one
   xiprev=xi1
LoopVmax:  do ni = 1,maxitnr

      if(ni>maxitnr/2)then
        xi1=half*(xi1+xiprev)
        er1=10.0d0*er1
      endif

      ! v^2(xi1) - maxvel^2
      fv= ((ssqr + (xi1*two*sdotb*sdotb + bsqr*sdotb*sdotb)/(xi1*xi1))&
        /((xi1+bsqr)*(xi1+bsqr)))-(one-dmaxvel )**2.0d0
      ! d(v^2(xi)-maxvel^2)/dxi
      dfv= -two * (sdotb*sdotb*(3.0d0*xi1*(xi1+bsqr)+bsqr*bsqr)+ssqr*xi1**3)/ &
                     ((xi1*(xi1+bsqr))**3.0d0)

      if(fv*dfv==zero) then
         if(fv==zero)then
            exit LoopVmax
         else
            !print *,'stop: dfv becomes zero, non-monotonic function of xi!!'
            ierror=8
            return
         endif
      else
        xiprev=xi1
        xi1   =xi1 -fv/dfv
        if(fv*dfv>zero)then
          ! xi-iterate decreased
          ! restrict to left
          xi1=max(xi1,dplus)
        else ! fv*dfv <0 
          ! xi-iterate increased
          ! restrict to right
          xi1=min(xi1,tau+d)
        endif
      endif
      er=dabs(fv/dfv)/xi1
      if((er<tolernr*er1).or.(dabs(fv/dfv)<absaccnr))exit LoopVmax
      niiter=ni
   enddo LoopVmax
endif 
 !=============================================!

if(niiter==maxitnr)then
 ! could not find consistent value of lower bound for xi compliant with maxvel
 print *,' could not find value of lower bound for xi compliant with maxvel'
 print *,'xi1=',xi1,'dplus=',dplus,' tau+d=',tau+d
 ierror=2
 return
endif

! we now compute f(xi1) and lfac(xi1)

call funcd(xi1,fl,df,lfacl,d,ssqr,tau,bsqr,sdotb,ierror)
if(ierror /=0)return
testsl=(xi1>=dplus*lfacl.and.lfacl<=lfacmax)
 
if(fl>zero)then
  !print *,'warning: lower bound non-negative f(xi) value!'
  !print *,'xi1=',xi1,' vs dplus=',dplus,'dplus*lfacl=',dplus*lfacl&
  !        ,'lfacl=',lfacl,'lfacmax=',lfacmax
  !print *,'fl=',fl,' niiter=',niiter,' vsqr(dplus)=',vsqr
  ierror=5
  return
endif
 
!--------------------------------------------------------!
! xi2 is meant to be the maximal bound for the bracket to be used in NR
! for this we take the value of E=tau+d, increased by smallp

xi2= max(tau+d+minp - half *bsqr,10.0d0*xi1)
niiter=-1

LoopxiMax : do ni=1,maxitnr
    ! we now compute f(xi2) and lfac(xi2)
    call funcd(xi2,fh,df,lfach,d,ssqr,tau,bsqr,sdotb,ierror)
    if(ierror /=0)return
    testsh=(xi2>=dplus*lfach.and.lfach<=lfacmax)

    ! maximal bound found when f(xi1) opposite sign of f(xi2)
    !, enforce consistency tests on xi2
    if (testsh.and.fh *fl <=zero) exit LoopxiMax
    !==== Zak 17/05 fast convergence ====!
    xi1=xi2
    fl=fh
    lfacl=lfach
    testsl=testsh
    !====================================!
    xi2=two*xi2
    niiter=ni
end do   LoopxiMax
!--------------------------------------------------------!

if(niiter == maxitnr .or. (fh*fl>zero))then
  ! could not find upper bound for NR on xi
  print *,'could not find upper bound for NR on xi'
  print *,'niiter=',niiter,' versus maxitnr=',maxitnr
  print *,'xi1=',xi1,' fl=',fl,' lfacl=',lfacl,' vs        dplus=',dplus
  print *,'xi2=',xi2,' fh=',fh,' lfach=',lfach,' vs tau+d+smallp=',tau+d+smallp
  ierror = 3
  return
end if

finished=.false.
if(fl==zero)then
  xii=xi1
  lfaci=lfacl
  finished=(xi1>=dplus*lfacl.and.lfacl<=lfacmax)
endif

if(fh==zero)then
  xii=xi2
  lfaci=lfach
  finished=(xi2>=dplus*lfach.and.lfach<=lfacmax)
endif

if(finished)then
  xi=xii
  lfac=lfaci
  return
end if


xil=xi1
xih=xi2
xii=half*(xih+xil)    !Initialize the guess for rootfinder
call funcd(xii,f,df,lfaci,d,ssqr,tau,bsqr,sdotb,ierror)
if(ierror /=0)return

er1 = one
nit = 0
niiter=-1
xiprev=xii
lastxiok=-one
lastlfacok=-one
lastf=-one

!--- Start iteration ---!
LoopNRRMHD :  do ni=1,maxitnr 
      nit = nit + 1
      if(nit>maxitnr/2)then
        ! mix the last  value for convergence
        xii=half*(xii+xiprev)
        ! relax accuracy requirement
        er1=10.0d0*er1
        ! following avoids decrease of accuracy requirement 
        ! *every* iteration step beyond maxitnr/2
        nit = nit - maxitnr/10
      endif
      call funcd(xii,f,df,lfaci,d,ssqr,tau,bsqr,sdotb,ierror) 
      if(ierror /=0)return
      tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
      if(tests)then
        lastxiok=xii
        lastlfacok=lfaci
        lastf=f
      endif
      if(f*df==zero) then
         if(f==zero) then
            tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
            if(tests)then
               exit ! zero found and consistency checks fullfilled
            else
               !print*,mype,' error 7: ', xii,dplus,lfaci,lfacmax
              ierror=7
              return
            endif
         else
            print *,'stop: df becomes zero, non-monotonic function of xi!!'
            ierror=8
            return
         endif
      else
        xiprev=xii
        xii   =xii -f/df
        if(f*df>zero)then
          ! xi-iterate decreased
          ! restrict to left
          xii=max(xii,xil)
        else
         ! xi-iterate increased
         ! restrict to right
          xii=min(xii,xih)
        endif
      endif
      er=dabs(f/df)/xii
      if((er<tolernr*er1).or.(dabs(f/df)<absaccnr))then
            call funcd(xii,f,df,lfaci,d,ssqr,tau,bsqr,sdotb,ierror) 
            tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
            if(tests)then
               exit LoopNRRMHD ! converged solution with tests ensured
            else
               !print*,mype,' error 7: ', xii,dplus,lfaci,lfacmax
              ierror=7
              return
            endif
      endif
      niiter=ni

enddo LoopNRRMHD

if(niiter==maxitnr) then
     ! no convergence of NR for xi, although zero bracketed
     print *,'no convergence of NR for xi, although zero bracketed'
     print *,'er=',er,'tolernr=',tolernr,'er1=',er1,'df=',df,'absaccnr=',absaccnr
     print *,'xii=',xii,' f=',f
     print *,'brackets xil=',xil,' and xih=',xih,' with fl fh=',fl,fh
     print *,'lastxiok=',lastxiok,'lastlfacok=',lastlfacok,'lastf=',lastf
     ierror=2
     return
endif ! niiter==maxitnr

 !===============================!
 ! final values for auxiliary variables are now passed to w-array
  xi=xii
  lfac=lfaci
 !===============================!
 ! we now perform some additional consistency checks
  
  if(xi<smallxi)then
    print *,'xi smaller than smallxi!!! '
    ierror=4
    return
  endif

 ! compute v*(xi+B^2)
if(ssqr .ne. zero)then
  temp=dsqrt({(s^C + b^C*sdotb/xi)*(s^C + b^C*sdotb/xi)+})
   {v^C=dsqrt(one-one/(lfac**2))*(s^C + b^C*sdotb/xi)/temp\}
  if(({v^C**2.0d0+})>one) then
    print *,'bsqr=',bsqr,'v=',v^C,'v^2=',{v^C**2+}
    print *,'v(xi+B^2)=',temp
    print *,'xi=',xi,'lfac=',lfac
    print *,' v^2>1!!! '
    ierror=6
    return
  end if
end if

contains
!=============================================================================
!DEC$ ATTRIBUTES INLINE :: funcd
subroutine funcd(xi,F,dF,lfac,d,ssqr,tau,bsqr,sdotb,ierror)

include 'amrvacdef.f'

double precision, intent(in)  :: xi,d,ssqr,tau,bsqr,sdotb
double precision, intent(out) :: F,dF,lfac
integer, intent(inout)        :: ierror
  
double precision  :: dlfac,sb2
double precision  :: vsqr,g1
!-----------------------------------------------------------------------------

sb2 = sdotb*sdotb
vsqr = ((ssqr + (xi*two*sb2 + bsqr*sb2)/(xi*xi))/((xi+bsqr)*(xi+bsqr)))
g1 = (eqpar(gamma_)-one)/eqpar(gamma_)


if (vsqr<one) then
   lfac = one/dsqrt(one-vsqr)
   dlfac = -lfac**3.0d0*(sb2*(3.0d0*xi*(xi+bsqr)+bsqr*bsqr)+ssqr*xi**3.0d0) &
           /((xi*(xi+bsqr))**3.0d0)

   F  = xi-tau-d+half*bsqr*(one+vsqr)-half*sb2/xi**2.0d0 &
       -g1*(xi - d*lfac)/lfac**2.0d0
   dF = one + sb2/xi**3.0d0 - g1/lfac**2.0d0  &
       +dlfac*((two*xi - lfac*d)*g1 + bsqr)/lfac**3.0d0

else 
   print *, mype,' Warning: erroneous input to funcd &
      since vsrq=',vsqr,' >=1'
   print *,mype, ' input values d, ssqr, tau, bsqr, sdotb:',d,ssqr,tau,bsqr,sdotb
   ierror =6
   return
end if

end subroutine funcd
!=============================================================================
end subroutine con2prim
!=============================================================================
subroutine addgeometry(qdt,ixI^L,ixO^L,wCT,w,x)

! Add geometrical source terms to w

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(in)       :: qdt
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(inout)    :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

integer                            :: iw,ix,idir, h1x^L{^NOONED, h2x^L}
double precision, dimension(ixG^T) :: tmp, sqrVdotB, sqrB
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
         tmp(ixO^S) = half*(sqrVdotB(ixO^S) + &
                sqrB(ixO^S)/(wCT(ixO^S,lfac_)**2)) +(eqpar(gamma_)-one)* &
                (wCT(ixO^S,xi_)-wCT(ixO^S,lfac_)*wCT(ixO^S,d_))/ &
                (eqpar(gamma_)*wCT(ixO^S,lfac_)**2)
 
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
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*wCT(ixO^S,psi_)/x(ixO^S,1)
}
      end select
      ! Divide by radius and add to w
      if (iw==sr_ .or.(iw==sphi_.or.iw==bphi_)) then
           w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
      end if
   end do

   case('spherical')
      h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
      do iw=1,nwflux
         select case(iw)
         ! s[s1]=((Stheta**2+Sphi**2-(vdotB)**2*(Btheta**2+Bphi**2))/(xi+b2)
         !       +2*ptotal-(btheta**2+bphi**2)/lfac**2)/r
         case(s1_)
           ! total pressure
           tmp(ixO^S) = half*(sqrVdotB(ixO^S) + &
                sqrB(ixO^S)/wCT(ixO^S,lfac_)**2.0d0) + (eqpar(gamma_)-one)* &
                (wCT(ixO^S,xi_)-wCT(ixO^S,lfac_)*wCT(ixO^S,d_))/ &
                (eqpar(gamma_)*wCT(ixO^S,lfac_)**2)

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
                sqrB(ixO^S)/wCT(ixO^S,lfac_)**2.0d0) + (eqpar(gamma_)-one)* &
                (wCT(ixO^S,xi_)-wCT(ixO^S,lfac_)*wCT(ixO^S,d_))/ &
                (eqpar(gamma_)*wCT(ixO^S,lfac_)**2)

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
         tmp(ixO^S)=wCT(ixO^S,psi_)*x(ixO^S,1) &
                 *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
                 /mygeo%dvolume(ixO^S)
}
{^NOONEC
         ! s[b2]=(sr*Btheta-stheta*Br)/(xi+B^2)/r
         !       + cot(theta)*psi/r
         case(b2_)
{#IFDEF GLM
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*wCT(ixO^S,psi_) &
                     *(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
                     /mygeo%dvolume(ixO^S)
}
            tmp(ixO^S)=(wCT(ixO^S,s1_)*wCT(ixO^S,b2_) &
                       -wCT(ixO^S,s2_)*wCT(ixO^S,b1_))&
                       /(wCT(ixO^S,xi_)+sqrB(ixO^S))
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
         if(iw==s1_{^NOONEC.or.iw==s2_.or.iw==b2_}{^IFTHREEC .or.iw==b3_ &
               .or.(iw==s3_.and..not.angmomfix)}) &
            w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
      end do
end select

end subroutine addgeometry
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
include 'amrvacdef.f'

  
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
      !!if(dabs(valuepoly)>accroot) then
        !!print *,'AFTER NR root =',root,' value polynomial=',valuepoly
        !!print *,'before NR root:',root1,' value=',valuepoly1
        !call mpistop("")
      !!endif
    endif
enddo

end subroutine zroots4
!=============================================================================
subroutine addsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

include 'amrvacdef.f'

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

include 'amrvacdef.f'

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

!=============================================================================
subroutine con2primHydro(lfac,xi,d,s^C,tau,ierror)
!use ieee_arithmetic
include 'amrvacdef.f'

! this is a copy of the HD iteration, where we solve for p via NR, and modified
! to give xi on output

double precision :: xi,lfac
double precision :: d,s^C,tau
integer          :: ierror

integer          :: ni,niiter
double precision :: pcurrent,pnew,pL,pR
double precision :: er,er1,ff,df,dp,v^C
double precision :: pmin,lfac2inv,pLabs,pRabs,pprev
double precision :: s2overcubeG2rh,sqrs
double precision :: xicurrent
double precision :: oldff1,oldff2
double precision :: Nff
double precision :: pleft,pright,pnewi
integer          ::nit,n2it,ni2,ni3
double precision :: h,dhdp
!-----------------------------------------------------------------------------
ierror=0
! ierror=0 : ok
!
! ierror<>0
!
! ierror=1 : error on entry: must have D>=minrho, tau>=smalltau
! ierror=2 : maxitnr reached without convergence
! ierror=3 : final pressure value < smallp or xi<smallxi during iteration
! ierror=4 : final v^2=1 hence problem as lfac=1/0
! ierror=5 : nonmonotonic function f?
! ierror=7 : stop due to strictnr violation

if(d<minrho .or. tau<smalltau) then
  ierror=1
  return
endif

! incase input pressure is not available or random value: replace by smallp

sqrs={s^C**2.0d0+}


! left and right brackets for p-range
pmin=dsqrt(sqrs)/(one-dmaxvel)-tau-d
pLabs=max(minp,pmin)
pRabs=1.0d99
! start value from input
pcurrent=pLabs

er1=one
pprev=pcurrent

! Fudge Parameters
oldff1=1.0d7  ! High number
oldff2=1.0d9  ! High number bigger then oldff1
n2it = 0
nit  = 0



LoopNR:  do ni=1,maxitnr
     nit = nit + 1
     !============= Controle ~1~=============!
     if(nit>maxitnr/4)then
        !print *,'ni,er,p',ni,er,pcurrent
        ! mix pressure value for convergence
        pcurrent=half*(pcurrent+pprev)
        ! relax accuracy requirement
        er1=10.*er1
        nit = nit - maxitnr/10
     endif
     !=======================================!

     niiter=ni  
     xicurrent=tau+d+pcurrent

     if(xicurrent<smallxi) then
       print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
       print *,'stop: too small xi iterate:',xicurrent
       print *,'for pressure iterate p',pcurrent
       print *,'pressure bracket pLabs pRabs',pLabs,pRabs
       print *,'iteration number:',ni
       print *,'values for d,s,tau,s^2:',d,s^C,tau,sqrs
       ierror=3 
       return
     endif

     {v^C=s^C/xicurrent\}
     lfac2inv=one - ({v^C**2.0d0+})
     if(lfac2inv>zero) then
       lfac=one/dsqrt(lfac2inv)
     else
       print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
       print *,'stop: negative or zero factor 1-v^2:',lfac2inv
       print *,'for pressure iterate p',pcurrent
       print *,'pressure bracket pL pR',pL,pR
       print *,'absolute pressure bracket pLabs pRabs',pLabs,pRabs
       print *,'iteration number:',ni
       print *,'values for d,s,tau,s^2:',d,s^C,tau,sqrs
       print *,'values for v,xi:',v^C,xicurrent
       ierror=4
       return
     endif
       
     s2overcubeG2rh=sqrs/(xicurrent**3.0d0)
     ff=(d*(lfac-one)-pcurrent-tau)*lfac2inv + pcurrent*govergminone
     df=(d*lfac - 2.0d0*xicurrent)*s2overcubeG2rh  &
          + govergminone - lfac2inv

     if (ff*df==zero) then
        if (ff==zero) then
            exit ! zero found
        else
            print *,'stop: df becomes zero, non-monotonic f(p)!!'
            ierror=5
            return
        endif
     else 
        pnew=pcurrent-ff/df
        if (ff*df>zero) then
            ! pressure iterate has decreased
            ! restrict to left 
            pnew=max(pnew,pLabs)
        else  ! ff*df<0
            ! pressure iterate has increased
            ! restrict to right 
            pnew=min(pnew,pRabs)
        endif
     endif
        

     ! handle special case where NR incorrectly believes in convergence
     if(pnew == pLabs .and. pcurrent==pnew .and. &
        abs(ff)> absaccnr .and. sqrs > zero)then
        pnewi=pnew
        ! try 2 higher pressure values to locate a sign change for f(p)
LoopCor:  do ni2=1,2
	   !=====================!
	   pcurrent=pnewi*500.0d0
  	   xicurrent=tau+d+pcurrent
     	   {v^C=s^C/xicurrent\}
     	   lfac2inv=one - ({v^C**2.0d0+})
 	   !=====================!
        
	   !=====================!
    	   if(lfac2inv>zero)then
	      lfac=one/dsqrt(lfac2inv)
	   else
              ierror=4
              return
	   endif
 	   !=====================!

	   !==== Calculate enthalpy and derivative ====!
           Nff=(d*(lfac-one)-pcurrent-tau)*lfac2inv + pcurrent*govergminone

	   !== Save old value of pressure ==!
	   pnewi=pcurrent
	   !================================!

	   !== find the interval where is the root ==!
	   if(Nff * ff <=zero)then
	      pnew=pcurrent
	      exit LoopCor
	   endif
	   !=========================================!
        enddo LoopCor

        !== No possible solution, correct all including the conservatives ==!
        if( Nff*ff>zero)then
     
           ! following is in accord with trick done in smallvalues
           d   = 2.0d0*(one + 10.0d0 * minrho) * minrho
           tau = 2.0d0*(one + 10.0d0 * smalltau) * smalltau
           {^C&s^C =zero;}
           pcurrent     = (eqpar(gamma_)-one)*tau
           xi=tau+d+pcurrent
           lfac = one

           if(strictnr)ierror=7
           ! leave the do loop here
           return
        endif
     endif
     !===============================================!
     dp=pcurrent-pnew
     er=two*dabs(dp)/(pnew+pcurrent)
     if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
     !===============================================!

     ! For very small values of pressure, NR algorithm is not efficient to
     ! find root, use Euler algorithm to find precise value of pressure
     if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.&
           ff * oldff1 < zero    .and.  dabs(ff)>absaccnr)then

       n2it=n2it+1
       if(n2it<=3) pcurrent=half*(pnew+pcurrent)
       if(n2it>3)then
         pright =pcurrent
         pleft=pprev
         pcurrent=half*(pleft+pright)
 Dicho:  do ni3=1,maxitnr
           !===================!
           xicurrent=tau+d+pcurrent
           {v^C=s^C/xicurrent\}
           lfac2inv=one - ({v^C**2.0d0+})
           if(lfac2inv>zero)then
             lfac=one/dsqrt(lfac2inv)
           else
             ierror=4
             return
           endif
           !===================!

           Nff=(d*(lfac-one)-pcurrent-tau)*lfac2inv + pcurrent*govergminone

	   !==== Iterate ====!
	   if(ff * Nff < zero)then
     		pleft=pcurrent
  	   else
		pright=pcurrent
    	   endif

	   pcurrent=half*(pleft+pright)
	   !==================!

	   !=== The iteration converge ===!
           if(2.0d0*dabs(pleft-pright)/(pleft+pright)< absaccnr &
	      .or. dabs(ff)<absaccnr)then
              pnew=pcurrent
	      exit LoopNR
            endif
	    !==============================!

	    !=== conserve the last value of Nff ===!
	    ff=Nff
	    !======================================!
         enddo    Dicho
       endif

     else
       !====== There is no problems, continue the NR iteration ======!
       pprev=pcurrent
       pcurrent=pnew
       !=============================================================!
     endif 

 
     !=== keep the values of the 2 last ff ===!
     oldff2=oldff1
     oldff1=ff
     !========================================!
enddo LoopNR
  
if(niiter==maxitnr)then
   ierror=2
   return
endif

if(pcurrent<minp) then
   ierror=3
   return
endif

!------------------------------!
xi=tau+d+pcurrent
{v^C = s^C/xi\}
lfac2inv=one - ({v^C**2.0d0+})
if(lfac2inv>zero) then
   lfac=one/dsqrt(lfac2inv)
else
   ierror=4
   return
endif

end subroutine con2primHydro
!=============================================================================
! end module amrvacphys/- srmhd
!#############################################################################
