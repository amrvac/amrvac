!#############################################################################
! module amrvacphys/- mhd

!=============================================================================
subroutine checkglobaldata
include 'amrvacdef.f'
!-----------------------------------------------------------------------------
minrho= max(zero,smallrho)
{#IFDEF ISO
minp=eqpar(adiab_)*minrho**eqpar(gamma_)
}
{#IFDEF GAMMA
minp  = max(zero,smallp)
smalle= minp/(eqpar(gamma_)-one)
}
end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

! set default values for entropy fixes

include 'amrvacdef.f'

integer :: il
!-----------------------------------------------------------------------------
{#IFDEF GAMMA
eqpar(gamma_)=5.d0/3.d0
}
{#IFDEF ISO
eqpar(gamma_)=1.d0
eqpar(adiab_)=1.d0
}

{#IFDEF GLM
eqpar(Cr_)=0.18d0
}

eqpar(eta_)=0.d0
eqpar(etahyper_)=0.d0
eqpar(etah_)=0.d0

do il=1,nw
   select case(il)
   case(fastRW_,fastLW_,slowRW_,slowLW_)
      entropycoef(il)=0.2d0
   case(alfvRW_,alfvLW_)
      entropycoef(il)=0.4d0
   case default
      entropycoef(il)= -one
   end select
end do

end subroutine initglobaldata
!=============================================================================
subroutine getaux(clipping,w,x,ixI^L,ixO^L,subname)

! Calculate auxilary variables ixO^L from non-auxiliary entries in w
! clipping can be set to .true. to e.g. correct unphysical pressures,
! densities, v>c,  etc.

include 'amrvacdef.f'

logical                :: clipping
integer                :: ixI^L, ixO^L
double precision       :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
character(len=*)       :: subname
!-----------------------------------------------------------------------------

end subroutine getaux
!=============================================================================
subroutine checkw(checkprimitive,ixI^L,ixO^L,w,flag)

include 'amrvacdef.f'

logical :: checkprimitive
integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
logical :: flag(ixI^S)

double precision :: tmp1(ixI^S)
!-----------------------------------------------------------------------------
flag(ixI^S)=.true.

{#IFDEF ENERGY
if(checkprimitive)then
   tmp1(ixO^S)=w(ixO^S,p_)
else
   ! First calculate kinetic energy*2=m**2/rho
   tmp1(ixO^S)=( ^C&w(ixO^S,m^C_)**2+ )/w(ixO^S,rho_)
   ! Add magnetic energy*2=b**2
   tmp1(ixO^S)=tmp1(ixO^S)+ ^C&w(ixO^S,b^C_)**2+
   ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
   tmp1(ixO^S)=(eqpar(gamma_)-one)*(w(ixO^S,e_)-half*tmp1(ixO^S))
endif

flag(ixO^S)=(tmp1(ixO^S)>=minp .and. w(ixO^S,rho_)>=minrho)
}
{#IFDEF ISO
flag(ixO^S)=(w(ixO^S,rho_)>=minrho)
}
end subroutine checkw
!=============================================================================
subroutine conserve(ixI^L,ixO^L,w,x,patchw)

! Transform primitive variables into conservative ones

include 'amrvacdef.f'

integer, intent(in)    :: ixI^L, ixO^L
double precision       :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
logical                :: patchw(ixI^S)
double precision       :: invgam
!-----------------------------------------------------------------------------

invgam=1.d0/(eqpar(gamma_)-one)
where(.not.patchw(ixO^S))
{#IFDEF GAMMA
   ! Calculate total energy from pressure, kinetic and magnetic energy
   w(ixO^S,e_)=w(ixO^S,p_)*invgam+&
      half*(w(ixO^S,rho_)*(^C&w(ixO^S,v^C_)**2+)+(^C&w(ixO^S,b^C_)**2+))
}
   ! Convert velocity to momentum
   ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,v^C_);
{#IFDEF TRACER
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,rho_)*w(ixO^S,tr^FL_)\}
}
end where

if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"conserve")

end subroutine conserve
!=============================================================================
subroutine conserven(ixI^L,ixO^L,w,patchw)

! Transform primitive variables into conservative ones

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(inout) :: w(ixI^S,nw)
logical, intent(in)             :: patchw(ixI^S)
double precision       :: invgam
!-----------------------------------------------------------------------------

invgam=1.d0/(eqpar(gamma_)-one)
where(.not.patchw(ixO^S))
{#IFDEF GAMMA
   ! Calculate total energy from pressure, kinetic and magnetic energy
   w(ixO^S,e_)=w(ixO^S,p_)*invgam+&
      half*(w(ixO^S,rho_)*(^C&w(ixO^S,v^C_)**2+)+(^C&w(ixO^S,b^C_)**2+))
}
   ! Convert velocity to momentum
   ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,v^C_);
{#IFDEF TRACER
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,rho_)*w(ixO^S,tr^FL_)\}
}
end where

end subroutine conserven
!=============================================================================
subroutine primitive(ixI^L,ixO^L,w,x)

! Transform conservative variables into primitive ones

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(inout) :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
! .. local ..
{#IFDEF GAMMA
integer, dimension(ixI^S)       :: patchierror
}
integer, dimension(ndim)       :: lowpindex
!-----------------------------------------------------------------------------
if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"primitive")

! Convert momentum to velocity
^C&w(ixO^S,v^C_)=w(ixO^S,m^C_)/w(ixO^S,rho_);
{#IFDEF ENERGY
! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
w(ixO^S,p_)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
       half*(({^C&w(ixO^S,v^C_)**2+})*w(ixO^S,rho_)&
       +{ ^C&w(ixO^S,b^C_)**2+}))
if(strictsmall) then
  if(any(w(ixO^S,p_)<minp)) then
    lowpindex=minloc(w(ixO^S,p_))
    ^D&lowpindex(^D)=lowpindex(^D)+ixOmin^D-1;
    write(*,*)'too small pressure = ',minval(w(ixO^S,p_)),' with limit=',&
    minp,' at x=',x(^D&lowpindex(^D),1:ndim),' array index=',lowpindex,&
    ' where E_k=',half*(^C&w(^D&lowpindex(^D),v^C_)**2+)*&
    w(^D&lowpindex(^D),rho_),&
    ' E_B=',half*(^C&w(^D&lowpindex(^D),b^C_)**2+),' E_total=',&
    w(^D&lowpindex(^D),p_)/(eqpar(gamma_)-one)+half*&
    (^C&w(^D&lowpindex(^D),v^C_)**2+)*w(^D&lowpindex(^D),rho_)+&
    half*(^C&w(^D&lowpindex(^D),b^C_)**2+),' w(1:nwflux)=',&
    w(^D&lowpindex(^D),1:nwflux),' when t=',t,' it=',it
    call mpistop("=== primitive pressure problem===")
  end if
else
  if (strictgetaux) then
     where(w(ixO^S,p_)<minp)
        w(ixO^S,p_)=minp
     endwhere
  else
     where(w(ixO^S,p_)<minp)
       patchierror(ixO^S) = 1
     elsewhere
       patchierror(ixO^S) = 0
     end where
     if (any(patchierror(ixO^S)/=0)) &
   call correctaux(ixI^L,ixO^L,w,x,patchierror,'primitive')
 end if
end if
}
{#IFDEF TRACER
! We got rho, Dtr, now we can get the tracers:
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,Dtr^FL_)/w(ixO^S,rho_)\}
}
end subroutine primitive
!=============================================================================
subroutine primitiven(ixI^L,ixO^L,w,patchw)

! Transform conservative variables into primitive ones

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(inout) :: w(ixI^S,nw)
logical, intent(in),dimension(ixI^S)   :: patchw
!-----------------------------------------------------------------------------

where(.not.patchw(ixO^S))
  ! Convert momentum to velocity
  ^C&w(ixO^S,v^C_)=w(ixO^S,m^C_)/w(ixO^S,rho_);
{#IFDEF GAMMA
  ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
  w(ixO^S,p_)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
         half*(({^C&w(ixO^S,v^C_)**2+})*w(ixO^S,rho_)&
         +{ ^C&w(ixO^S,b^C_)**2+}))
}
{#IFDEF TRACER
! We got rho, Dtr, now we can get the tracers:
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,Dtr^FL_)/w(ixO^S,rho_)\}
}
end where

end subroutine primitiven
!=============================================================================
subroutine e_to_rhos(ixI^L,ixO^L,w,x)

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision,intent(inout)  :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------
{#IFDEF GAMMA
w(ixO^S,rhos_)=(eqpar(gamma_)-one)*w(ixO^S,rho_)**(one-eqpar(gamma_)) &
               *(w(ixO^S,e_)-half*((^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_) &
                                   +(^C&w(ixO^S,b^C_)**2+)))
}
{#IFDEF ISO
call mpistop("e_to_rhos can not be used with eos=iso !")
}
end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixI^L,ixO^L,w,x)

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------
{#IFDEF GAMMA
w(ixO^S,e_)=(one/(eqpar(gamma_)-one))*w(ixO^S,rho_)**(eqpar(gamma_)-one)*&
             w(ixO^S,rhos_)+half*((^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_)+&
             (^C&w(ixO^S,b^C_)**2+))
}
{#IFDEF ISO
call mpistop("rhos_to_e can not be used with eos=iso !")
}
end subroutine rhos_to_e
!=============================================================================
subroutine internalenergy(ixI^L,ixO^L,w,x,ie)

! get internal energy

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision                :: ie(ixI^S)
!-----------------------------------------------------------------------------
{#IFDEF GAMMA
ie(ixO^S)=( ^C&w(ixO^S,m^C_)**2+ )/w(ixO^S,rho_)+ ^C&w(ixO^S,b^C_)**2+
! internal energy=e-0.5*(2ek+2eb)
ie(ixO^S)=w(ixO^S,e_)-half*ie(ixO^S)
}
{#IFDEF ISO
call mpistop("internalenergy can not be used with eos=iso !")
}

end subroutine internalenergy
!=============================================================================
subroutine getv(w,x,ixI^L,ixO^L,idims,v)

! Calculate v_idim=m_idim/rho within ixO^L

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, idims
double precision, intent(in)  :: w(ixI^S,nw)
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(out) :: v(ixI^S)
!-----------------------------------------------------------------------------

v(ixO^S)=w(ixO^S,m0_+idims)/w(ixO^S,rho_)

end subroutine getv
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixI^L,ixO^L,idims,cmax,cmin,needcmin)

! Calculate cmax_idim=csound+abs(v_idim) within ixO^L

include 'amrvacdef.f'

logical :: new_cmax,needcmin
integer, intent(in) :: ixI^L, ixO^L, idims
double precision :: w(ixI^S,nw), cmax(ixI^S), cmin(ixI^S)
double precision, intent(in)    :: x(ixI^S,1:ndim)

double precision :: csound2(ixI^S), cfast2(ixI^S), AvMinCs2(ixI^S), tmp(ixI^S)
{#IFDEF HALL double precision ::  bmag2(ixI^S), kmax}
!, vh(ixI^S), va2(ixI^S), Omegai(ixI^S)
!-----------------------------------------------------------------------------

!Direction independent part of getcmax:
call getcsound2(w,x,ixI^L,ixO^L,csound2)

if (B0field) then
   cfast2(ixO^S)=( ^C&(w(ixO^S,b^C_)+myB0%w(ixO^S,^C))**2+ ) &
                 /w(ixO^S,rho_)+csound2(ixO^S)
   AvMinCs2(ixO^S)=cfast2(ixO^S)**2-4.0d0*csound2(ixO^S)&
                *((myB0%w(ixO^S,idims)+w(ixO^S,b0_+idims))**2) &
                 /w(ixO^S,rho_)
else
   cfast2(ixO^S)=( ^C&w(ixO^S,b^C_)**2+ )/w(ixO^S,rho_)+csound2(ixO^S)
   AvMinCs2(ixO^S)=cfast2(ixO^S)**2-4.0d0*csound2(ixO^S)&
                *((w(ixO^S,b0_+idims))**2)/w(ixO^S,rho_)
endif

where(AvMinCs2(ixO^S)<zero)
   AvMinCs2(ixO^S)=zero
end where

AvMinCs2(ixO^S)=dsqrt(AvMinCs2(ixO^S))

{#IFNDEF HALL
tmp(ixO^S) = dsqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
}
{#IFDEF HALL
if (.not. B0field) then
   bmag2(ixO^S)={^C& w(ixO^S,b^C_)**2 +}
else
   bmag2(ixO^S)={^C& (w(ixO^S,b^C_)+myB0%w(ixO^S,^C))**2 +}
end if
! take the Hall velocity into account:
! most simple estimate, high k limit:
! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
kmax = dpi/min({dxlevel(^D)},bigdouble)*half
tmp(ixO^S) = max(dsqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
     eqpar(etah_) * sqrt(bmag2(ixO^S))/w(ixO^S,rho_)*kmax)

! Phase speed of the fast mode in direction of B can be estimated to be smaller than vh, (this is more strict):
!Omegai(ixO^S) = sqrt(bmag2(ixO^S))/eqpar(etah_)
!va2(ixO^S) = bmag2(ixO^S) / w(ixO^S,rho_)
!vh(ixO^S) = one/sqrt(3.0d0)*&
!     sqrt(csound2(ixO^S)+two*va2(ixO^S)+two*sqrt( (csound2(ixO^S)-va2(ixO^S))**2 &
!     + kmax**4 * va2(ixO^S)**4/Omegai(ixO^S)**4 & 
!     + kmax**2*va2(ixO^S)**2 * (4.0d0*va2(ixO^S)-csound2(ixO^S))/Omegai(ixO^S)**2) &
!     + kmax**2*va2(ixO^S)**2/Omegai(ixO^S)**2)
!tmp(ixO^S) = max(dsqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
!     vh(ixO^S)) 
}

if(needcmin)then
  cmax(ixO^S)=max(+tmp(ixO^S) &
      +(w(ixO^S,m0_+idims)/w(ixO^S,rho_)),zero)
  cmin(ixO^S)=min(-tmp(ixO^S) &
      +(w(ixO^S,m0_+idims)/w(ixO^S,rho_)),zero)
else
  cmax(ixO^S)=tmp(ixO^S) &
      +dabs(w(ixO^S,m0_+idims)/w(ixO^S,rho_))
end if

end subroutine getcmax
!=============================================================================
subroutine getpthermal(w,x,ixI^L,ixO^L,p)

! Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw), p(ixI^S)
double precision, intent(in)    :: x(ixI^S,1:ndim)
integer, dimension(ixI^S)       :: patchierror
integer, dimension(ndim)       :: lowpindex
!-----------------------------------------------------------------------------
{#IFDEF GAMMA
if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,'getpthermal')

! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
p(ixO^S)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
       half*(({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_)&
       +{^C&w(ixO^S,b^C_)**2+}))

! Clip off negative pressure if smallp is set
if(strictsmall) then
  if(any(p(ixO^S)<minp)) then
    lowpindex=minloc(p(ixO^S))
    ^D&lowpindex(^D)=lowpindex(^D)+ixOmin^D-1;
    write(*,*)'too small pressure = ',minval(p(ixO^S)),' with limit=',minp,&
    ' at x=',x(^D&lowpindex(^D),1:ndim),' array index=',lowpindex,&
    ' where E_k=',half*(^C&w(^D&lowpindex(^D),m^C_)**2+)/&
    w(^D&lowpindex(^D),rho_),&
    ' E_B=',half*(^C&w(^D&lowpindex(^D),b^C_)**2+),' E_total=',w(^D&lowpindex(^D),e_),&
    ' w(1:nwflux)=',w(^D&lowpindex(^D),1:nwflux),&
    ' when t=',t,' it=',it
    call mpistop("=== strictsmall in getpthermal ===")
  end if
else
  if (strictgetaux) then
     where(p(ixO^S)<minp)
        p(ixO^S)=minp
     endwhere
  else
     where(p(ixO^S)<minp)
       patchierror(ixO^S) = 1
     elsewhere
       patchierror(ixO^S) = 0
     end where
     if (any(patchierror(ixO^S)/=0))then
       call correctaux(ixI^L,ixO^L,w,x,patchierror,'getpthermal')
       where(patchierror(ixO^S)/=0)
         p(ixO^S)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
        half*(({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_)+{^C&w(ixO^S,b^C_)**2+}))
       end where
     end if
 end if
end if
}
{#IFDEF ISO
p(ixO^S)=eqpar(adiab_)*w(ixO^S,rho_)**eqpar(gamma_)
}
end subroutine getpthermal
!=============================================================================
subroutine getcsound2prim(w,x,ixI^L,ixO^L,csound2)

! Calculate the square of the thermal sound speed csound2 within ixO^L
! from the primitive variables in w.
! csound2=gamma*p/rho

include 'amrvacdef.f'

integer, intent(in)             :: ixO^L, ixI^L
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(in)    :: w(ixI^S,nw)
double precision, intent(out)   :: csound2(ixI^S)
!-----------------------------------------------------------------------------
{#IFDEF GAMMA
csound2(ixO^S)=eqpar(gamma_)*w(ixO^S,p_)/w(ixO^S,rho_)
}
{#IFDEF ISO
csound2(ixO^S)=eqpar(gamma_)*eqpar(adiab_)*w(ixO^S,rho_)**(eqpar(gamma_)-one)
}
end subroutine getcsound2prim
!=============================================================================
subroutine getcsound2(w,x,ixI^L,ixO^L,csound2)

! Calculate the square of the thermal sound speed csound2 within ixO^L.
! csound2=gamma*p/rho

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(out)   :: csound2(ixI^S)
!-----------------------------------------------------------------------------
{#IFDEF GAMMA
call getpthermal(w,x,ixI^L,ixO^L,csound2)
csound2(ixO^S)=eqpar(gamma_)*csound2(ixO^S)/w(ixO^S,rho_)
}
{#IFDEF ISO
csound2(ixO^S)=eqpar(gamma_)*eqpar(adiab_)*w(ixO^S,rho_)**(eqpar(gamma_)-one)
}
end subroutine getcsound2
!=============================================================================
subroutine getptotal(w,x,ixI^L,ixO^L,p)

! Calculate total pressure within ixO^L including magnetic pressure
! p=(g-1)*e-0.5*(g-1)*m**2/rho+(1-0.5*g)*b**2

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(out)   :: p(ixI^S)
!.. local ..
{#IFDEF GAMMA
double precision :: gamma
}
integer, dimension(ixI^S)       :: patchierror
integer, dimension(ndim)       :: lowpindex
!-----------------------------------------------------------------------------
if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,'getptotal')

{#IFDEF GAMMA
gamma=eqpar(gamma_)
p(ixO^S)=(one-half*gamma)*( ^C&w(ixO^S,b^C_)**2+ )+(gamma-one)*&
  (w(ixO^S,e_)-half*(^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_))

if(strictsmall) then
  if(any(p(ixO^S)<minp)) then
    lowpindex=minloc(p(ixO^S))
    ^D&lowpindex(^D)=lowpindex(^D)+ixOmin^D-1;
    write(*,*)'too small pressure = ',minval(p(ixO^S)),' at x=',&
    x(^D&lowpindex(^D),1:ndim),lowpindex,' with limit=',minp,' where E_k=',&
    half*(^C&w(^D&lowpindex(^D),m^C_)**2+)/w(^D&lowpindex(^D),rho_),' E_B=',&
    half*(^C&w(^D&lowpindex(^D),b^C_)**2+),'E_total=',w(^D&lowpindex(^D),e_),&
    ' w(1:nwflux)=',w(^D&lowpindex(^D),1:nwflux),' when t=',t,' it=',it
    call mpistop("=== strictsmall in getptotal ===")
  end if
else
  if (strictgetaux) then
     where(p(ixO^S)<minp)
        p(ixO^S)=minp
     endwhere
  else
     where(p(ixO^S)<minp)
       patchierror(ixO^S) = 1
     elsewhere
       patchierror(ixO^S) = 0
     end where
     if (any(patchierror(ixO^S)/=0))then
       call correctaux(ixI^L,ixO^L,w,x,patchierror,'getptotal')
       where(patchierror(ixO^S)/=0)
          p(ixO^S)=(one-half*gamma)*( ^C&w(ixO^S,b^C_)**2+ )+(gamma-one)*&
            (w(ixO^S,e_)-half*(^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_))
       end where
     end if
 end if
end if
}
{#IFDEF ISO
p(ixO^S)=eqpar(adiab_)*w(ixO^S,rho_)**eqpar(gamma_)+(^C&w(ixO^S,b^C_)**2+)*half
}


end subroutine getptotal
!=============================================================================
subroutine getfluxforhllc(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L, iw, idims
double precision, intent(in)    :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(out)   :: f(ixI^S,1:nwflux)
!.. local ..
logical :: transport
integer          :: idir
double precision :: tmp(ixI^S),tmp2(ixI^S){#IFDEF HALL , vh(ixI^S,1:3)}
!-----------------------------------------------------------------------------
transport=.true.

if (B0field) then
   if (iw==m0_+idims{#IFDEF ENERGY .or. iw==e_}) tmp(ixO^S)={^C&myB0%w(ixO^S,^C)*w(ixO^S,b^C_)+}
end if

select case (iw)
   ! f_i[rho]=v_i*rho
   case (rho_)
      f(ixO^S,iw)=zero
{#IFDEF TRACER
{  case (tr^FL_)
      f(ixO^S,iw)=zero\}
}
   ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
   {case (m^C_)
      if (idims==^C) then
         call getptotal(w,x,ixI^L,ixO^L,tmp2)
         f(ixO^S,iw)=tmp2(ixO^S)-w(ixO^S,b0_+idims)*w(ixO^S,b^C_)
         if (B0field) f(ixO^S,iw)=f(ixO^S,iw)+tmp(ixO^S)
      else
         f(ixO^S,iw)= -w(ixO^S,b^C_)*w(ixO^S,b0_+idims)
      end if
      if (B0field) then
         f(ixO^S,iw)=f(ixO^S,iw)-myB0%w(ixO^S,idims)*w(ixO^S,b^C_) &
                                -w(ixO^S,b0_+idims)*myB0%w(ixO^S,^C)
      end if\}
{#IFDEF ENERGY
   ! f_i[e]=v_i*e+(m_i*ptotal-b_i*(b_k*m_k))/rho
   case (e_)
      call getptotal(w,x,ixI^L,ixO^L,tmp2)
      f(ixO^S,iw)=(w(ixO^S,m0_+idims)*tmp2(ixO^S)- &
          w(ixO^S,b0_+idims)*( ^C&w(ixO^S,b^C_)*w(ixO^S,m^C_)+ ))/w(ixO^S,rho_)
{#IFDEF HALL
   ! f_i[e]= f_i[e] + vh_i*(b_k*b_k) - b_i*(vh_k*b_k)
      if (eqpar(etah_)>zero) then
         call getvh(w,x,ixI^L,ixO^L,vh)
         f(ixO^S,iw) = f(ixO^S,iw) + vh(ixO^S,idims)* &
           (^C&w(ixO^S,b^C_)*w(ixO^S,b^C_)+ ) &
           - w(ixO^S,b0_+idims) * (^C&vh(ixO^S,^C)*w(ixO^S,b^C_)+ )
         if (B0field) then
            f(ixO^S,iw) = f(ixO^S,iw) &
                 + vh(ixO^S,idims) * tmp(ixO^S) &
                 - (^C&vh(ixO^S,^C)*w(ixO^S,b^C_)+ ) * myB0%w(ixO^S,idims)
         end if
      end if
}
      if (B0field) then
         f(ixO^S,iw)=f(ixO^S,iw)+ tmp(ixO^S) &
                          *w(ixO^S,m0_+idims)/w(ixO^S,rho_) &
                          -( ^C&w(ixO^S,m^C_)*w(ixO^S,b^C_)+ )/w(ixO^S,rho_) &
                          *myB0%w(ixO^S,idims)
      end if
}
   ! f_i[b_k]=v_i*b_k-m_k/rho*b_i
   {case (b^C_)
      if (idims==^C) then
         ! f_i[b_i] should be exactly 0, so we do not use the transport flux
{#IFNDEF GLM f(ixO^S,iw)=zero }
{#IFDEF GLM   f(ixO^S,iw)=w(ixO^S,psi_)}
         transport=.false.
      else
         f(ixO^S,iw)= -w(ixO^S,b0_+idims)*w(ixO^S,m^C_)/w(ixO^S,rho_)
         if (B0field) then
            f(ixO^S,iw)=f(ixO^S,iw) &
                     +w(ixO^S,m0_+idims)/w(ixO^S,rho_)*myB0%w(ixO^S,^C) &
                     -myB0%w(ixO^S,idims)*w(ixO^S,m^C_)/w(ixO^S,rho_)
         end if
{#IFDEF HALL
   ! Hall MHD:
   ! f_i[b_k] = f_i[b_k] + vh_i*b_k - vh_k*b_i
         if (eqpar(etah_)>zero) then
            call getvh(w,x,ixI^L,ixO^L,vh)
            if (B0field) then
               f(ixO^S,iw) = f(ixO^S,iw) &
                    - vh(ixO^S,^C)*(w(ixO^S,b0_+idims)+myB0%w(ixO^S,idims)) &
                    + vh(ixO^S,idims)*(w(ixO^S,b^C_)+myB0%w(ixO^S,^C))
            else 
               f(ixO^S,iw) = f(ixO^S,iw) &
                    - vh(ixO^S,^C)*w(ixO^S,b0_+idims) &
                    + vh(ixO^S,idims)*w(ixO^S,b^C_)
            end if
         end if
}
      end if\}
{#IFDEF GLM 
!f_i[psi]=Ch^2*b_{i}
! Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
   case(psi_)
{#IFNDEF FCT
       f(ixO^S,iw)=cmax_global**2*w(ixO^S,b0_+idims)
}{#IFDEF FCT
       f(ixO^S,iw)=zero
}
       transport=.false. }
   case default
      call mpistop("Error in getflux: unknown flow variable!")
end select

end subroutine getfluxforhllc
!=============================================================================
subroutine getflux(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L, iw, idims
double precision, intent(in)    :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision,intent(out)    :: f(ixI^S)
!.. local ..
logical :: transport
double precision :: tmp(ixI^S){#IFDEF HALL , vh(ixI^S,1:3)}
integer          :: idirmin, idir
!-----------------------------------------------------------------------------
transport=.true.

if (B0field) then
   if (iw==m0_+idims{#IFDEF ENERGY .or. iw==e_}) tmp(ixO^S)={^C&myB0%w(ixO^S,^C)*w(ixO^S,b^C_)+}
end if

select case (iw)
   ! f_i[rho]=v_i*rho
   case (rho_)
      f(ixO^S)=zero
{#IFDEF TRACER
{  case (tr^FL_)
      f(ixO^S)=zero\}
}
   ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
   {case (m^C_)
      if (idims==^C) then
         call getptotal(w,x,ixI^L,ixO^L,f)
         f(ixO^S)=f(ixO^S)-w(ixO^S,b0_+idims)*w(ixO^S,b^C_)
         if (B0field) f(ixO^S)=f(ixO^S)+tmp(ixO^S)
      else
         f(ixO^S)= -w(ixO^S,b^C_)*w(ixO^S,b0_+idims)
      end if
      if (B0field) then
         f(ixO^S)=f(ixO^S)-myB0%w(ixO^S,idims)*w(ixO^S,b^C_) &
                          -w(ixO^S,b0_+idims)*myB0%w(ixO^S,^C)
      end if\}
{#IFDEF ENERGY
   ! f_i[e]=v_i*e+(m_i*ptotal-b_i*(b_k*m_k))/rho
   case (e_)
      call getptotal(w,x,ixI^L,ixO^L,f)
      f(ixO^S)=(w(ixO^S,m0_+idims)*f(ixO^S)- &
          w(ixO^S,b0_+idims)*( ^C&w(ixO^S,b^C_)*w(ixO^S,m^C_)+ ))/w(ixO^S,rho_)
{#IFDEF HALL
   ! f_i[e]= f_i[e] + vh_i*(b_k*b_k) - b_i*(vh_k*b_k)
      if (eqpar(etah_)>zero) then
         call getvh(w,x,ixI^L,ixO^L,vh)
         f(ixO^S) = f(ixO^S) + vh(ixO^S,idims)*&
              (^C&w(ixO^S,b^C_)*w(ixO^S,b^C_)+ )&
              - w(ixO^S,b0_+idims) * (^C&vh(ixO^S,^C)*w(ixO^S,b^C_)+ )
         if (B0field) then
            f(ixO^S) = f(ixO^S) &
                 + vh(ixO^S,idims) * tmp(ixO^S) &
                 - (^C&vh(ixO^S,^C)*w(ixO^S,b^C_)+ ) * myB0%w(ixO^S,idims)
         end if
      end if
}
      if (B0field) then
         f(ixO^S)=f(ixO^S)+ tmp(ixO^S) &
                          *w(ixO^S,m0_+idims)/w(ixO^S,rho_) &
                          -( ^C&w(ixO^S,m^C_)*w(ixO^S,b^C_)+ )/w(ixO^S,rho_) &
                          *myB0%w(ixO^S,idims)
      end if
}
   ! f_i[b_k]=v_i*b_k-m_k/rho*b_i
   {case (b^C_)
      if (idims==^C) then
         ! f_i[b_i] should be exactly 0, so we do not use the transport flux
{#IFNDEF GLM f(ixO^S)=zero}
{#IFDEF GLM   f(ixO^S)=w(ixO^S,psi_)}
         transport=.false.
      else
         f(ixO^S)= -w(ixO^S,b0_+idims)*w(ixO^S,m^C_)/w(ixO^S,rho_)
         if (B0field) then
            f(ixO^S)=f(ixO^S) &
                     +w(ixO^S,m0_+idims)/w(ixO^S,rho_)*myB0%w(ixO^S,^C) &
                     -myB0%w(ixO^S,idims)*w(ixO^S,m^C_)/w(ixO^S,rho_)
         end if
{#IFDEF HALL
         ! Hall MHD:
         ! f_i[b_k] = f_i[b_k] + vh_i*b_k - vh_k*b_i 
         if (eqpar(etah_)>zero) then
            call getvh(w,x,ixI^L,ixO^L,vh)
            if (B0field) then
               f(ixO^S) = f(ixO^S) &
                    - vh(ixO^S,^C)*(w(ixO^S,b0_+idims)+myB0%w(ixO^S,idims)) &
                    + vh(ixO^S,idims)*(w(ixO^S,b^C_)+myB0%w(ixO^S,^C))
            else 
               f(ixO^S) = f(ixO^S) &
                    - vh(ixO^S,^C)*w(ixO^S,b0_+idims) &
                    + vh(ixO^S,idims)*w(ixO^S,b^C_)
            end if
         end if
}
      end if\}
{#IFDEF GLM 
!f_i[psi]=Ch^2*b_{i}
! Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
   case(psi_)
{#IFNDEF FCT
       f(ixO^S)=cmax_global**2*w(ixO^S,b0_+idims)
}{#IFDEF FCT
       f(ixO^S)=zero
}
     transport=.false.\}
   case default
      call mpistop("Error in getflux: unknown flow variable!")
end select

end subroutine getflux
!=============================================================================
subroutine addsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)
logical, intent(in)             :: qsourcesplit
!.. local ..
double precision:: dx^D
!-----------------------------------------------------------------------------

dx^D=dxlevel(^D);
if(qsourcesplit .eqv. ssplitresis) then
! Sources for resistivity in eqs. for e, B1, B2 and B3
if(dabs(eqpar(eta_))>smalldouble)then
   if (.not.slab) call mpistop("no resistivity in non-slab geometry")
   if(compactres)then
      call addsource_res1(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
   else
      call addsource_res2(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
   endif
endif

if (eqpar(etahyper_)>0.d0)then
   call addsource_hyperres(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
end if
endif

{^NOONED
if(qsourcesplit .eqv. ssplitdivb) then
! Sources related to div B
select case (typedivbfix)
{#IFDEF GLM
case ('glm1')
   call addsource_glm1(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
case ('glm2')
   call addsource_glm2(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D) 
case ('glm3')
   call addsource_glm3(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)\}
case ('powel')
   call addsource_powel(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
case ('janhunen')
   call addsource_janhunen(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
case ('linde')
   call addsource_linde(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
case ('lindejanhunen')
   call addsource_linde(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
   call addsource_janhunen(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
case ('lindepowel')
   call addsource_linde(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
   call addsource_powel(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
end select
endif
}
end subroutine addsource
!=============================================================================
subroutine addsource_res1(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add resistive source to w within ixO 
! Uses 3 point stencil (1 neighbour) in each direction, non-conservative
! If the fourthorder precompiler flag is set, uses fourth order central difference for the laplacian. Then the stencil is 5 (2 neighbours).  
include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in)    :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
!.. local ..
integer :: ix^L,idir,jdir,kdir,idirmin,iw,idims,jxO^L,hxO^L,ix
{#IFDEF FOURTHORDER  integer :: lxO^L, kxO^L}

double precision :: tmp(ixI^S),tmp2(ixI^S)

! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S)
double precision :: gradeta(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

{#IFNDEF FOURTHORDER
! Calculating resistive sources involve one extra layer
ix^L=ixO^L^LADD1;
}{#IFDEF FOURTHORDER
ix^L=ixO^L^LADD2;

}
if(ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
   call mpistop("Error in addsource_res1: Non-conforming input limits")

! Calculate current density and idirmin
call getcurrent(wCT,ixI^L,ixO^L,idirmin,current)

if(eqpar(eta_)>zero)then
   eta(ix^S)=eqpar(eta_)
   gradeta(ixO^S,1:ndim)=zero
else
   call specialeta(wCT,ixI^L,ix^L,idirmin,x,current,eta)
   ! assumes that eta is not function of current?
   do idims=1,ndim
      call gradient(eta,ixI^L,ixO^L,idims,tmp)
      gradeta(ixO^S,idims)=tmp(ixO^S)
   enddo
endif

do idir=1,ndir

   ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
{#IFNDEF FOURTHORDER
   tmp(ixO^S)=zero
   tmp2(ixI^S)=wCT(ixI^S,b0_+idir)
   do idims=1,ndim
      jxO^L=ixO^L+kr(idims,^D);
      hxO^L=ixO^L-kr(idims,^D);
      tmp(ixO^S)=tmp(ixO^S)+&
           (tmp2(jxO^S)-2.0d0*tmp2(ixO^S)+tmp2(hxO^S))/dxlevel(idims)**2
   enddo
}
{#IFDEF FOURTHORDER
   tmp(ixO^S)=zero
   tmp2(ixI^S)=wCT(ixI^S,b0_+idir)
   do idims=1,ndim
      lxO^L=ixO^L+2*kr(idims,^D);
      jxO^L=ixO^L+kr(idims,^D);
      hxO^L=ixO^L-kr(idims,^D);
      kxO^L=ixO^L-2*kr(idims,^D);
      tmp(ixO^S)=tmp(ixO^S)+&
           (-tmp2(lxO^S)+16.0d0*tmp2(jxO^S)-30.0d0*tmp2(ixO^S)+16.0d0*tmp2(hxO^S)-tmp2(kxO^S)) &
           /(12.0d0 * dxlevel(idims)**2)
   enddo
}

   ! Multiply by eta
   tmp(ixO^S)=tmp(ixO^S)*eta(ixO^S)

   ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
   if(eqpar(eta_)<zero)then
      do jdir=1,ndim; do kdir=idirmin,3
         if(lvc(idir,jdir,kdir)/=0)then
            if(lvc(idir,jdir,kdir)==1)then
               tmp(ixO^S)=tmp(ixO^S)-gradeta(ixO^S,jdir)*current(ixO^S,kdir)
            else
               tmp(ixO^S)=tmp(ixO^S)+gradeta(ixO^S,jdir)*current(ixO^S,kdir)
            endif
         endif
      enddo; enddo
   endif

   ! Add sources related to eta*laplB-grad(eta) x J to B and e
   do iw=1,nw
      if(iw==b0_+idir)then
         ! dB_idir/dt+=tmp
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)
{#IFDEF ENERGY
      else if(iw==e_)then
         ! de/dt+=B.tmp
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)*wCT(ixO^S,b0_+idir)
}
      endif
   end do  ! iw
enddo ! idir


{#IFDEF ENERGY
! de/dt+=eta*J**2
tmp(ixO^S)=zero
do idir=idirmin,3
   tmp(ixO^S)=tmp(ixO^S)+current(ixO^S,idir)**2
enddo
w(ixO^S,e_)=w(ixO^S,e_)+qdt*eta(ixO^S)*tmp(ixO^S)

if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"addsource_res1")
}
end subroutine addsource_res1
!=============================================================================
subroutine addsource_res2(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add resistive source to w within ixO 
! Uses 5 point stencil (2 neighbours) in each direction, conservative

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in)    :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
!.. local ..
integer :: ix^L,idir,jdir,kdir,idirmin,iw,idims,idirmin1

double precision :: tmp(ixI^S),tmp2(ixI^S)

! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S),curlj(ixI^S,1:3)
double precision :: tmpvec(ixI^S,1:3),tmpvec2(ixI^S,1:ndir)
!-----------------------------------------------------------------------------

ix^L=ixO^L^LADD2;

if(ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
   call mpistop("Error in addsource_res2: Non-conforming input limits")

ix^L=ixO^L^LADD1;
! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
! Determine exact value of idirmin while doing the loop.
call getcurrent(wCT,ixI^L,ix^L,idirmin,current)

if(eqpar(eta_)>zero)then
   eta(ix^S)=eqpar(eta_)
else
   call specialeta(wCT,ixI^L,ix^L,idirmin,x,current,eta)
endif

! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
tmpvec(ix^S,1:ndir)=zero
do jdir=idirmin,3
   tmpvec(ix^S,jdir)=current(ix^S,jdir)*eta(ix^S)*qdt
enddo
call curlvector(tmpvec,ix^L,ixO^L,curlj,idirmin1,1,3)
{^C& w(ixO^S,b^C_) = w(ixO^S,b^C_)-curlj(ixO^S,b^C_-b0_)\}

{#IFDEF ENERGY
! de/dt= +div(B x Jeta)
tmpvec2(ix^S,1:ndir)=zero
do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
   if(lvc(idir,jdir,kdir)/=0)then
      tmp(ix^S)=wCT(ix^S,b0_+jdir)*current(ix^S,kdir)*eta(ix^S)*qdt
      if(lvc(idir,jdir,kdir)==1)then
         tmpvec2(ix^S,idir)=tmpvec2(ix^S,idir)+tmp(ix^S)
      else
         tmpvec2(ix^S,idir)=tmpvec2(ix^S,idir)-tmp(ix^S)
      endif
   endif
enddo; enddo; enddo
!select case(typediv)
!case("central")
   call divvector(tmpvec2,ixI^L,ixO^L,tmp)
!case("limited")
!   call divvectorS(tmpvec2,ixI^L,ixO^L,tmp)
!end select
w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)

if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"addsource_res2")
}
end subroutine addsource_res2
!=============================================================================
subroutine addsource_hyperres(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add Hyper-resistive source to w within ixO 
! Uses 9 point stencil (4 neighbours) in each direction.

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in)    :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
!.. local ..
double precision                :: current(ixI^S,7-2*ndir:3)
double precision                :: tmpvec(ixI^S,1:3),tmpvec2(ixI^S,1:3),tmp(ixI^S),ehyper(ixI^S,1:3)
integer                         :: ix^L,idir,jdir,kdir,idirmin,idirmin1
!-----------------------------------------------------------------------------
ix^L=ixO^L^LADD3;
if(ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
   call mpistop("Error in addsource_hyperres: Non-conforming input limits")

call getcurrent(wCT,ixI^L,ix^L,idirmin,current)
tmpvec(ix^S,1:ndir)=zero
do jdir=idirmin,3
   tmpvec(ix^S,jdir)=current(ix^S,jdir)
enddo

ix^L=ixO^L^LADD2;
call curlvector(tmpvec,ix^L^LADD1,ix^L,tmpvec2,idirmin1,1,3)

ix^L=ixO^L^LADD1;
tmpvec(ix^S,1:ndir)=zero
call curlvector(tmpvec2,ix^L^LADD1,ix^L,tmpvec,idirmin1,1,3)
ehyper(ix^S,1:ndir) = - tmpvec(ix^S,1:ndir)*eqpar(etahyper_)

ix^L=ixO^L;
tmpvec2(ix^S,1:ndir)=zero
call curlvector(ehyper,ix^L^LADD1,ix^L,tmpvec2,idirmin1,1,3)

{^C& w(ixO^S,b^C_) = w(ixO^S,b^C_)-tmpvec2(ixO^S,^C)*qdt\}


{#IFDEF ENERGY
! de/dt= +div(B x Ehyper)
ix^L=ixO^L^LADD1;
tmpvec2(ix^S,1:ndir)=zero
do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
   tmpvec2(ix^S,idir) = tmpvec(ix^S,idir)&
        + lvc(idir,jdir,kdir)*wCT(ix^S,b0_+jdir)*ehyper(ix^S,kdir)
enddo; enddo; enddo
tmp(ixO^S)=zero
call divvector(tmpvec2,ixI^L,ixO^L,tmp)
w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)*qdt

if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"addsource_hyperres")
}
end subroutine addsource_hyperres
!=============================================================================
subroutine getcurrent(w,ixI^L,ix^L,idirmin,current)

! Calculate idirmin and the idirmin:3 components of the common current array
! make sure that dxlevel(^D) is set correctly.
include 'amrvacdef.f'

integer, parameter:: idirmin0=7-2*ndir
integer :: ix^L, idirmin, ixI^L
double precision :: w(ixI^S,1:nw)

! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
double precision :: current(ixI^S,7-2*ndir:3),bvec(ixI^S,1:ndir)
!-----------------------------------------------------------------------------

if(B0field) then
 ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
else
 ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
endif
!bvec(ixI^S,1:ndir)=w(ixI^S,b0_+1:b0_+ndir)
call curlvector(bvec,ixI^L,ix^L,current,idirmin,idirmin0,ndir)

end subroutine getcurrent
!=============================================================================
subroutine getdt(w,ixI^L,ix^L,dtnew,dx^D,x)

! If resistivity is not zero, check diffusion time limit for dt

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ix^L
double precision, intent(out)   :: dtnew
double precision, intent(in)    :: dx^D
double precision, intent(in)    :: w(ixI^S,1:nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!.. local ..
integer :: idirmin,idims
double precision :: dxarr(ndim)
double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S) 
!double precision :: {#IFDEF HALL , dthall }
!-----------------------------------------------------------------------------
dtnew=bigdouble

^D&dxarr(^D)=dx^D;
^D&dxlevel(^D)=dx^D;
if(eqpar(eta_)>zero)then
   dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/eqpar(eta_)
else if(eqpar(eta_)<zero)then
   call getcurrent(w,ixI^L,ix^L,idirmin,current)
   call specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)
   dtnew=bigdouble
   do idims=1,ndim
      dtnew=min(dtnew,&
                dtdiffpar/(smalldouble+maxval(eta(ix^S)/dxarr(idims)**2)))
   enddo
endif

if(eqpar(etahyper_)>zero)then
   dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/eqpar(etahyper_),dtnew)
end if 

!{#IFDEF HALL
! This is now covered in cmax, so no need.
!if(eqpar(etah_)>zero)then
!   call getdthall(w,x,ixI^L,ix^L,dx^D,dthall)
!   dtnew=min(dtnew,dthall)
!end if
!}
end subroutine getdt
!=============================================================================
subroutine ppmflatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L,ixL^L,ixR^L
double precision, intent(in)  :: w(ixI^S,nw),d2w(ixI^S,1:nwflux)
double precision, intent(inout) :: drho(ixI^S),dp(ixI^S)
!-----------------------------------------------------------------------------
{#IFDEF GAMMA
if(useprimitive)then
 drho(ixO^S) =eqpar(gamma_)*dabs(d2w(ixO^S,rho_))&
              /min(w(ixL^S,rho_),w(ixR^S,rho_))
 dp(ixO^S) = dabs(d2w(ixO^S,p_))/min(w(ixL^S,p_),w(ixR^S,p_))
end if
}
{#IFDEF ISO
call mpistop("PPM with flatcd=.true. can not be used with eos=iso !")
}
end subroutine ppmflatcd
!=============================================================================
subroutine ppmflatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)

include 'amrvacdef.f'

! based on Mignone and Miller and Collela 2002
! PPM flattening at shocks: we use total pressure and not thermal pressure 

integer, intent(in)           :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
integer, intent(in)           :: idims
double precision, intent(in)  :: w(ixI^S,nw)

double precision, intent(inout) :: drho(ixI^S),dp(ixI^S),dv(ixI^S)
double precision :: ptot(ixI^S)
!-----------------------------------------------------------------------------

{#IFDEF GAMMA
if(useprimitive)then
   ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
   ptot(ixO^S)=w(ixO^S,p_)+half*( ^C&w(ixO^S,b^C_)**2+ )
   where (dabs(ptot(ixRR^S)-ptot(ixLL^S))>smalldouble)
      drho(ixO^S) = dabs((ptot(ixR^S)-ptot(ixL^S))&
                        /(ptot(ixRR^S)-ptot(ixLL^S)))
   elsewhere
      drho(ixO^S) = zero
   end where

   !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26 
   !  use "dp" to save squared sound speed, assume primitive in w
   dp(ixO^S)=(eqpar(gamma_)*w(ixO^S,p_)/w(ixO^S,rho_))

   dp(ixO^S)  = dabs(ptot(ixR^S)-ptot(ixL^S))&
                /(w(ixO^S,rho_)*dp(ixO^S))
   ! recycle ptot to store v
   ptot(ixI^S)= w(ixI^S,v0_+idims)
   call gradient(ptot,ixI^L,ixO^L,idims,dv)
end if
}
{#IFDEF ISO
call mpistop("PPM with flatsh=.true. can not be used with eos=iso !")
}

end subroutine ppmflatsh
!=============================================================================
subroutine addgeometry(qdt,ixI^L,ixO^L,wCT,w,x)

! Add geometrical source terms to w

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
!.. local ..
integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
double precision :: tmp(ixI^S)
logical          :: angmomfix=.false.
!-----------------------------------------------------------------------------

select case (typeaxial)
case ('slab')
   ! No source terms in slab symmetry
case ('cylindrical')
   do iw=1,nwflux
      select case (iw)
      ! s[mr]=(ptotal-Bphi**2+mphi**2/rho)/radius
      case (mr_)
         call getptotal(wCT,x,ixI^L,ixO^L,tmp)
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
         tmp(ixO^S)=zero
{^IFPHI
         tmp(ixO^S)= &
            -wCT(ixO^S,bphi_)**2+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_)

      ! s[mphi]=(-mphi*mr/rho+Bphi*Br)/radius
      case (mphi_)
        tmp(ixO^S)= &
           -wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_) &
           +wCT(ixO^S,bphi_)*wCT(ixO^S,br_)

      ! s[Bphi]=((Bphi*mr-Br*mphi)/rho)/radius
      case (bphi_)
        tmp(ixO^S)=(wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
                   -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
               /wCT(ixO^S,rho_)
}
{#IFDEF GLM      ! s[br]=psi/radius
      case (br_)
         tmp(ixO^S)=wCT(ixO^S,psi_)
}
      end select

      ! Divide by radius and add to w
      if (iw==mr_{#IFDEF GLM .or.iw==br_}{^IFPHI .or.iw==mphi_.or.iw==bphi_}) then
          w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
      endif

   end do
case ('spherical')
   h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
   do iw=1,nwflux
      select case (iw)
      ! s[m1]=((mtheta**2+mphi**2)/rho+2*ptotal-(Btheta**2+Bphi**2))/r
      case (m1_)
         call getptotal(wCT,x,ixI^L,ixO^L,tmp)
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)+{^C&myB0_cell%w(ixO^S,^C)*wCT(ixO^S,b^C_)+}
         end if
         ! For nonuniform Cartesian grid this provides hydrostatic equil.
         tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
                 *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
                 /mygeo%dvolume(ixO^S){&^CE&
               +wCT(ixO^S,m^CE_)**2/wCT(ixO^S,rho_)-wCT(ixO^S,b^CE_)**2 }
         if (B0field.and.ndir>1) then
            tmp(ixO^S)=tmp(ixO^S){^CE&-2.0d0*myB0_cell%w(ixO^S,^CE) &
                                            *wCT(ixO^S,b^CE_)|}
         end if
{^NOONEC
      ! s[m2]=-(mr*mtheta/rho-Br*Btheta)/r
      !       + cot(theta)*(mphi**2/rho+(p+0.5*B**2)-Bphi**2)/r
      case (m2_)
}
{^NOONED
         call getptotal(wCT,x,ixI^L,ixO^L,tmp)
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)+{^C&myB0_cell%w(ixO^S,^C)*wCT(ixO^S,b^C_)+}
         end if
         ! This will make hydrostatic p=const an exact solution
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S) &
                     *(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
                     /mygeo%dvolume(ixO^S)
}
{^NOONEC
         tmp(ixO^S)=-(wCT(ixO^S,m1_)*wCT(ixO^S,m2_)/wCT(ixO^S,rho_) &
                     -wCT(ixO^S,b1_)*wCT(ixO^S,b2_))
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)+myB0_cell%w(ixO^S,1)*wCT(ixO^S,b2_) &
                                 +wCT(ixO^S,b1_)*myB0_cell%w(ixO^S,2)
         end if
}
{^IFTHREEC
{^NOONED
         tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,m3_)**2/wCT(ixO^S,rho_) &
                        -wCT(ixO^S,b3_)**2)*dcos(x(ixO^S,2)) &
                                           /dsin(x(ixO^S,2))
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)-2.0d0*myB0_cell%w(ixO^S,3)*wCT(ixO^S,b3_)&
                       *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
         end if
}
      ! s[m3]=-(mphi*mr/rho-Bphi*Br)/r
      !       -cot(theta)*(mtheta*mphi/rho-Btheta*Bphi)/r
      case (m3_)
         if (.not.angmomfix) then
            tmp(ixO^S)=-(wCT(ixO^S,m3_)*wCT(ixO^S,m1_)/wCT(ixO^S,rho_) &
                     -wCT(ixO^S,b3_)*wCT(ixO^S,b1_)) {^NOONED &
                   -(wCT(ixO^S,m2_)*wCT(ixO^S,m3_)/wCT(ixO^S,rho_) &
                     -wCT(ixO^S,b2_)*wCT(ixO^S,b3_)) &
                   *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
            if (B0field) then
               tmp(ixO^S)=tmp(ixO^S)+myB0_cell%w(ixO^S,1)*wCT(ixO^S,b3_) &
                          +wCT(ixO^S,b1_)*myB0_cell%w(ixO^S,3) {^NOONED &
                          +(myB0_cell%w(ixO^S,2)*wCT(ixO^S,b3_) &
                            +wCT(ixO^S,b2_)*myB0_cell%w(ixO^S,3)) &
                   *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
            end if
         end if
}
{#IFDEF GLM
      ! s[b1]=2*psi/r
      case (b1_)
         tmp(ixO^S)=2.0d0*wCT(ixO^S,psi_)
}
{^NOONEC
      ! s[b2]=(mr*Btheta-mtheta*Br)/rho/r
      !       + cot(theta)*psi/r
      case (b2_)
         tmp(ixO^S)=(wCT(ixO^S,m1_)*wCT(ixO^S,b2_) &
                    -wCT(ixO^S,m2_)*wCT(ixO^S,b1_))/wCT(ixO^S,rho_)
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,m1_)*myB0_cell%w(ixO^S,2) &
                       -wCT(ixO^S,m2_)*myB0_cell%w(ixO^S,1))/wCT(ixO^S,rho_)
         end if
{#IFDEF GLM
            tmp(ixO^S)=tmp(ixO^S) &
                 + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
}
}
{^IFTHREEC
      ! s[b3]=(mr*Bphi-mphi*Br)/rho/r
      !       -cot(theta)*(mphi*Btheta-mtheta*Bphi)/rho/r
      case (b3_)
         tmp(ixO^S)=(wCT(ixO^S,m1_)*wCT(ixO^S,b3_) &
                 -wCT(ixO^S,m3_)*wCT(ixO^S,b1_))/wCT(ixO^S,rho_) {^NOONED &
                -(wCT(ixO^S,m3_)*wCT(ixO^S,b2_) &
                 -wCT(ixO^S,m2_)*wCT(ixO^S,b3_))*dcos(x(ixO^S,2)) &
                               /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,m1_)*myB0_cell%w(ixO^S,3) &
               -wCT(ixO^S,m3_)*myB0_cell%w(ixO^S,1))/wCT(ixO^S,rho_){^NOONED &
               -(wCT(ixO^S,m3_)*myB0_cell%w(ixO^S,2) &
                -wCT(ixO^S,m2_)*myB0_cell%w(ixO^S,3))*dcos(x(ixO^S,2)) &
                               /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
         end if
}
      end select
      ! Divide by radius and add to w
      if (iw==m1_{#IFDEF GLM .or.iw==b1_}{^NOONEC.or.iw==m2_.or.iw==b2_}&
        {^IFTHREEC .or.iw==b3_ .or.(iw==m3_.and..not.angmomfix)}) &
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
   end do
end select

end subroutine addgeometry
!=============================================================================
! end module amrvacphys/- mhd
!#############################################################################
