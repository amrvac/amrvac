!##############################################################################
! module amrvacphys/- hd

{#IFNDEF DUST
INCLUDE:amrvacnul/addsource.t
}
{#IFNDEF DUST
INCLUDE:amrvacnul/getdt.t
}
!=============================================================================
subroutine checkglobaldata

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
{#IFDEF ISO
if(eqpar(gamma_)<=zero) call mpistop ("gamma negative not ok")
if(eqpar(adiab_)<zero) call mpistop ("adiab strict negative not ok")
minrho= max(zero,smallrho)
minp=eqpar(adiab_)*minrho**eqpar(gamma_)
}
{#IFDEF GAMMA
if(eqpar(gamma_)<=zero.or.eqpar(gamma_)==one) call mpistop ("gamma negative or 1 not ok")
minp  = max(zero,smallp)
minrho= max(zero,smallrho)
smalle= minp/(eqpar(gamma_)-one)
}
{#IFDEF DUST
if(eqpar(mu_)<=zero) call mpistop ("mu (molecular weight) negative not ok")
minrhod= max(zero,smallrhod)
}

end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

! set default values for entropy fixes for 'yee' type

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
{#IFDEF DUST
eqpar(mu_) = one
mhcgspar=1.6733D-24
kbcgspar=1.38065D-16
}

do il=1,nw
   select case(il)
   case(soundRW_,soundLW_)
      entropycoef(il)= 0.2d0
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
integer,intent(in)              :: ixI^L, ixO^L
double precision, intent(inout) :: w(ixI^S,nw)
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
logical :: flag(ixG^T)

double precision :: tmp(ixG^T)
!-----------------------------------------------------------------------------
flag(ixG^T)=.true.

{#IFDEF GAMMA
if(checkprimitive)then
  flag(ixO^S)=(w(ixO^S,p_)>=minp .and. w(ixO^S,rho_)>=minrho)
else
  tmp(ixO^S)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
         half*( ^C&w(ixO^S,m^C_)**2+ )/w(ixO^S,rho_))
  flag(ixO^S)=(tmp(ixO^S)>=minp .and. w(ixO^S,rho_)>=minrho)
endif
}
{#IFDEF ISO
if(eqpar(adiab_)/=zero)then
  flag(ixO^S)=(w(ixO^S,rho_)>=zero)
endif
}

end subroutine checkw
!=============================================================================
subroutine conserve(ixI^L,ixO^L,w,x,patchw)

! Transform primitive variables into conservative ones

include 'amrvacdef.f'

integer, intent(in)    :: ixI^L, ixO^L
double precision       :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
logical                :: patchw(ixG^T)
double precision       :: invgam
!-----------------------------------------------------------------------------

invgam=1.d0/(eqpar(gamma_)-one)
where(.not.patchw(ixO^S))
{#IFDEF GAMMA
   ! Calculate total energy from pressure and kinetic energy
   w(ixO^S,e_)=w(ixO^S,p_)*invgam+ &
           half*w(ixO^S,rho_)*(^C&w(ixO^S,v^C_)**2+)
}
   ! Convert velocity to momentum
   ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,v^C_);
{#IFDEF TRACER
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,rho_)*w(ixO^S,tr^FL_)\}
}
{#IFDEF DUST
  ! Convert dust velocity to dust momentum
   {^DS&{^C&w(ixO^S,m^Cd^DS_)=w(ixO^S,rhod^DS_)*w(ixO^S,v^Cd^DS_);}\}
}
end where

if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"conserve")

end subroutine conserve
!=============================================================================
subroutine conserven(ixI^L,ixO^L,w,patchw)

! Transform primitive variables into conservative ones
! Idem to conserve, no smallvalues call

include 'amrvacdef.f'

integer, intent(in)    :: ixI^L, ixO^L
double precision       :: w(ixI^S,nw)
logical                :: patchw(ixG^T)
double precision       :: invgam
!-----------------------------------------------------------------------------

invgam=1.d0/(eqpar(gamma_)-one)
where(.not.patchw(ixO^S))
{#IFDEF GAMMA
   ! Calculate total energy from pressure and kinetic energy
   w(ixO^S,e_)=w(ixO^S,p_)*invgam+ &
           half*w(ixO^S,rho_)*(^C&w(ixO^S,v^C_)**2+)
}
   ! Convert velocity to momentum
   ^C&w(ixO^S,m^C_)=w(ixO^S,rho_)*w(ixO^S,v^C_);
{#IFDEF TRACER
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,rho_)*w(ixO^S,tr^FL_)\}
}
{#IFDEF DUST
   ! Convert dust velocity to dust momentum
   {^DS&{^C&w(ixO^S,m^Cd^DS_)=w(ixO^S,rhod^DS_)*w(ixO^S,v^Cd^DS_);}\}
}
end where

end subroutine conserven
!=============================================================================
subroutine primitive(ixI^L,ixO^L,w,x)

! Transform conservative variables into primitive ones

include 'amrvacdef.f'

integer, intent(in)    :: ixI^L, ixO^L
double precision       :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)

integer, dimension(ixG^T)       :: patchierror
integer, dimension(ndim)        :: lowpindex
!-----------------------------------------------------------------------------
if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"primitive")
{#IFDEF GAMMA
! compute pressure 
w(ixO^S,p_)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
         half*( ^C&w(ixO^S,m^C_)**2+ )/w(ixO^S,rho_))
}
{#IFDEF DUST
! Convert dust momentum to dust velocity
{^DS&where(w(ixO^S,rhod^DS_)>minrhod)
  {^C&w(ixO^S,v^Cd^DS_)=w(ixO^S,m^Cd^DS_)/w(ixO^S,rhod^DS_);}
else where
  {^C&w(ixO^S,v^Cd^DS_)=zero;}
end where\}
}
{#IFDEF ENERGY
! Convert momentum to velocity
^C&w(ixO^S,v^C_)=w(ixO^S,m^C_)/w(ixO^S,rho_);
if(strictsmall) then
  if(any(w(ixO^S,p_)<minp)) then
    lowpindex=minloc(w(ixO^S,p_))
    ^D&lowpindex(^D)=lowpindex(^D)+ixOmin^D-1;
    write(*,*)'too small pressure = ',minval(w(ixO^S,p_)),' with limit=',minp,&
    ' at x=',x(^D&lowpindex(^D),1:ndim),lowpindex,' where E_k=',&
    half*(^C&w(^D&lowpindex(^D),m^C_)**2+)/w(^D&lowpindex(^D),rho_),&
    ' E_total=',w(^D&lowpindex(^D),e_),&
    ' w(1:nwflux)=',w(^D&lowpindex(^D),1:nwflux),' when t=',t,' it=',it
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
     else where
       patchierror(ixO^S) = 0
     end where
     if (any(patchierror(ixO^S)/=0)) &
   call correctaux(ixI^L,ixO^L,w,x,patchierror,'primitive')
 end if
end if
}

{#IFDEF ISO
if(eqpar(adiab_)>zero)then
  ! Convert momentum to velocity
  ^C&w(ixO^S,v^C_)=w(ixO^S,m^C_)/w(ixO^S,rho_);
else
  ! case of zero temperature: allow zero density
  where(w(ixO^S,rho_)/=zero)
        ^C&w(ixO^S,v^C_)=w(ixO^S,m^C_)/w(ixO^S,rho_);
  elsewhere
        ^C&w(ixO^S,v^C_)=zero;
  endwhere
endif
}

{#IFDEF TRACER
! We got rho, Dtr, now we can get the tracers:
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,Dtr^FL_)/w(ixO^S,rho_)\}
}

end subroutine primitive
!=============================================================================
subroutine primitiven(ixI^L,ixO^L,w,patchw)

! Transform conservative variables into primitive ones
! Idem to primitive, no smallvalues call

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
logical, intent(in),dimension(ixG^T)   :: patchw
!-----------------------------------------------------------------------------
{#IFDEF ISO
if(eqpar(adiab_)>zero)then
  where(.not.patchw(ixO^S))
  ! Convert momentum to velocity
    ^C&w(ixO^S,v^C_)=w(ixO^S,m^C_)/w(ixO^S,rho_);
  endwhere
else
  ! case of zero temperature: allow zero density
  where(w(ixO^S,rho_)/=zero.and..not.patchw(ixO^S))
        ^C&w(ixO^S,v^C_)=w(ixO^S,m^C_)/w(ixO^S,rho_);
  elsewhere(w(ixO^S,rho_)==zero.and..not.patchw(ixO^S))
        ^C&w(ixO^S,v^C_)=zero;
  endwhere
endif
}
{#IFDEF GAMMA
where(.not.patchw(ixO^S))
  ! compute pressure 
  w(ixO^S,p_)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
               half*( ^C&w(ixO^S,m^C_)**2+ )/w(ixO^S,rho_))
  ! Convert momentum to velocity
  {^C&w(ixO^S,v^C_)=w(ixO^S,m^C_)/w(ixO^S,rho_);}
end where
}
{#IFDEF TRACER
where(.not.patchw(ixO^S))
! We got rho, Dtr, now we can get the tracers:
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,Dtr^FL_)/w(ixO^S,rho_)\}
end where
}
{#IFDEF DUST
where(.not.patchw(ixO^S))
! Convert dust momentum to dust velocity
    {^DS&where(w(ixO^S,rhod^DS_)>minrhod)
        {^C&w(ixO^S,v^Cd^DS_)=w(ixO^S,m^Cd^DS_)/w(ixO^S,rhod^DS_);}
    else where
        {^C&w(ixO^S,v^Cd^DS_)=zero;}
    end where
    \}
end where
}

end subroutine primitiven
!=============================================================================
subroutine e_to_rhos(ixI^L,ixO^L,w,x)

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------
{#IFDEF GAMMA
w(ixO^S,rhos_)=(eqpar(gamma_)-one)*w(ixO^S,rho_)**(one-eqpar(gamma_)) &
               *(w(ixO^S,e_)-half*(^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_))
}
{#IFDEF ISO
call mpistop("energy from entropy can not be used with -eos=iso !")
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
w(ixO^S,e_)=w(ixO^S,rho_)**(eqpar(gamma_)-one)*w(ixO^S,rhos_) &
            /(eqpar(gamma_)-one) +half*(^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_)
}
{#IFDEF ISO
call mpistop("entropy from energy can not be used with -eos=iso !")
}
end subroutine rhos_to_e
!=============================================================================
subroutine ppmflatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L,ixL^L,ixR^L
double precision, intent(in)  :: w(ixI^S,nw),d2w(ixG^T,1:nwflux)

double precision, intent(inout) :: drho(ixG^T),dp(ixG^T)
!-----------------------------------------------------------------------------
{#IFDEF GAMMA
if(useprimitive)then
 drho(ixO^S) =eqpar(gamma_)*dabs(d2w(ixO^S,rho_))&
              /min(w(ixL^S,rho_),w(ixR^S,rho_))
 dp(ixO^S) = dabs(d2w(ixO^S,p_))/min(w(ixL^S,p_),w(ixR^S,p_))
end if
}
{#IFDEF ISO
call mpistop("PPM with flatcd=.true. can not be used with -eos=iso !")
}
end subroutine ppmflatcd
!=============================================================================
subroutine ppmflatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
integer, intent(in)           :: idims
double precision, intent(in)  :: w(ixI^S,nw)

double precision, intent(inout) :: drho(ixG^T),dp(ixG^T),dv(ixG^T)
double precision :: v(ixG^T)
!-----------------------------------------------------------------------------
{#IFDEF GAMMA
if(useprimitive)then
   ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
   where (dabs(w(ixRR^S,p_)-w(ixLL^S,p_))>smalldouble)
      drho(ixO^S) = dabs((w(ixR^S,p_)-w(ixL^S,p_))&
                        /(w(ixRR^S,p_)-w(ixLL^S,p_)))
   else where
      drho(ixO^S) = zero
   end where

   !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26 
   !  use "dp" to save squared sound speed, assuming primitives
   dp(ixO^S)=(eqpar(gamma_)*w(ixO^S,p_)/w(ixO^S,rho_))

   dp(ixO^S) = dabs(w(ixR^S,p_)-w(ixL^S,p_))&
                /(w(ixO^S,rho_)*dp(ixO^S))
   v(ixI^S)  = w(ixI^S,v0_+idims)
   call gradient(v,ixO^L,idims,dv)
end if
}
{#IFDEF ISO
call mpistop("PPM with flatsh=.true. can not be used with -eos=iso !")
}
end subroutine ppmflatsh
!=============================================================================
subroutine getv(w,x,ixI^L,ixO^L,idims,v)

! Calculate v_idim=m_idim/rho within ixO^L

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, idims
double precision :: w(ixI^S,nw), v(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------
{#IFNDEF ISO
v(ixO^S)=w(ixO^S,m0_+idims)/w(ixO^S,rho_)
}
{#IFDEF ISO
if(eqpar(adiab_)>zero)then
   v(ixO^S)=w(ixO^S,m0_+idims)/w(ixO^S,rho_)
else
   ! case of zero temperature: allow zero density
   where(w(ixO^S,rho_)/=zero)
        v(ixO^S)=w(ixO^S,m0_+idims)/w(ixO^S,rho_)
   elsewhere
        v(ixO^S)=zero
   endwhere
endif
}
end subroutine getv
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixI^L,ixO^L,idims,cmax,cmin,needcmin)

! Calculate cmax_idim=csound+abs(v_idim) within ixO^L

include 'amrvacdef.f'

logical :: new_cmax,needcmin
integer, intent(in) :: ixI^L, ixO^L, idims
double precision :: w(ixI^S,nw), cmax(ixG^T), cmin(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)

double precision :: csound(ixG^T){#IFDEF DUST ,speeddust(ixG^T,1:^NDS)}
!-----------------------------------------------------------------------------
{#IFDEF DUST
select case(idims)
{case(^C)
    {^DS&where(w(ixO^S,rhod^DS_)>minrhod)
        speeddust(ixO^S,^DS)=w(ixO^S,m^Cd^DS_)/w(ixO^S,rhod^DS_);
    else where
        speeddust(ixO^S,^DS)=zero;
    end where\}
\}
end select
}
{#IFDEF GAMMA
call getpthermal(w,x,ixI^L,ixO^L,csound)
csound(ixO^S)=sqrt(eqpar(gamma_)*csound(ixO^S)/w(ixO^S,rho_))
if(needcmin)then
  cmax(ixO^S)=max(w(ixO^S,m0_+idims)/w(ixO^S,rho_)+csound(ixO^S), &
    {#IFDEF DUST {^DS&speeddust(ixO^S,^DS)}, } zero)
  cmin(ixO^S)=min(w(ixO^S,m0_+idims)/w(ixO^S,rho_)-csound(ixO^S),&
    {#IFDEF DUST {^DS&speeddust(ixO^S,^DS)}, } zero)
else
  cmax(ixO^S)=max(csound(ixO^S)+dabs(w(ixO^S,m0_+idims)/w(ixO^S,rho_)), &
    {#IFDEF DUST {^DS&dabs(speeddust(ixO^S,^DS))}, } zero)
endif
}

{#IFDEF ISO
if(eqpar(adiab_)>zero)then
  call getpthermal(w,x,ixI^L,ixO^L,csound)
  csound(ixO^S)=sqrt(eqpar(gamma_)*csound(ixO^S)/w(ixO^S,rho_))
  if(needcmin)then
    cmax(ixO^S)=max(w(ixO^S,m0_+idims)/w(ixO^S,rho_)+csound(ixO^S),&
    {#IFDEF DUST {^DS&speeddust(ixO^S,^DS)}, } zero)
    cmin(ixO^S)=min(w(ixO^S,m0_+idims)/w(ixO^S,rho_)-csound(ixO^S),&
    {#IFDEF DUST {^DS&speeddust(ixO^S,^DS)}, } zero)
  else
    cmax(ixO^S)=max(csound(ixO^S)+abs(w(ixO^S,m0_+idims)/w(ixO^S,rho_)), &
    {#IFDEF DUST {^DS&dabs(speeddust(ixO^S,^DS))}, } zero)
  endif
else
  ! case of zero temperature: allow zero density
  if(needcmin)then
    where(w(ixO^S,rho_)/=zero)
       cmax(ixO^S)= &
    max(w(ixO^S,m0_+idims)/w(ixO^S,rho_),{#IFDEF DUST {^DS&speeddust(ixO^S,^DS)}, } zero)
       cmin(ixO^S)= &
    min(w(ixO^S,m0_+idims)/w(ixO^S,rho_),{#IFDEF DUST {^DS&speeddust(ixO^S,^DS)}, } zero)
    elsewhere
      cmax(ixO^S)={#IFDEF DUST max({^DS&speeddust(ixO^S,^DS)}, } zero{#IFDEF DUST )}
      cmin(ixO^S)={#IFDEF DUST min({^DS&speeddust(ixO^S,^DS)}, } zero{#IFDEF DUST )}
    endwhere
  else
    where(w(ixO^S,rho_)/=zero)
      cmax(ixO^S)={#IFDEF DUST max({^DS&dabs(speeddust(ixO^S,^DS))}, } abs(w(ixO^S,m0_+idims)/w(ixO^S,rho_)){#IFDEF DUST )}
    elsewhere
      cmax(ixO^S)={#IFDEF DUST max({^DS&dabs(speeddust(ixO^S,^DS))}, } zero{#IFDEF DUST )}
    endwhere
  endif
endif
}
end subroutine getcmax

!=============================================================================
subroutine getpthermal(w,x,ixI^L,ixO^L,p)

! Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision                :: w(ixI^S,nw), p(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
integer, dimension(ixG^T)       :: patchierror
integer, dimension(ndim)        :: lowpindex
!-----------------------------------------------------------------------------
if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"getpthermal")

{#IFDEF GAMMA
p(ixO^S)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
         half*({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_))
}
{#IFDEF ISO
p(ixO^S)=eqpar(adiab_)*w(ixO^S,rho_)**eqpar(gamma_)
}

if(any(p(ixO^S)<minp)) then
  if(strictsmall) then
    lowpindex=minloc(p(ixO^S))
    ^D&lowpindex(^D)=lowpindex(^D)+ixOmin^D-1;
    write(*,*)'too small pressure = ',minval(p(ixO^S)),' with limit=',minp,&
    ' at x=',x(^D&lowpindex(^D),1:ndim),lowpindex,' where E_k=',&
    half*(^C&w(^D&lowpindex(^D),m^C_)**2+)/w(^D&lowpindex(^D),rho_)
{#IFDEF GAMMA
    write(*,*)' E_total=',w(^D&lowpindex(^D),e_)
}
    write(*,*)' w(1:nwflux)=',w(^D&lowpindex(^D),1:nwflux),' when t=',t,' it=',it
    call mpistop("=== strictsmall in getpthermal ===")
  else 
    if(strictgetaux) then
      where(p(ixO^S)<minp)
        p(ixO^S)=minp
{#IFDEF ISO
        w(ixO^S,rho_)=minrho
        {^C&w(ixO^S,m^C_)=zero;}
}
      endwhere
    else
      where(p(ixO^S)<minp)
        patchierror(ixO^S) = 1
      else where
        patchierror(ixO^S) = 0
      end where
      if(any(patchierror(ixO^S)/=0)) then
        call correctaux(ixI^L,ixO^L,w,x,patchierror,'getpthermal')
        where(patchierror(ixO^S)/=0)
{#IFDEF GAMMA
          p(ixO^S)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
          half*({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_))
}
{#IFDEF ISO
          p(ixO^S)=eqpar(adiab_)*w(ixO^S,rho_)**eqpar(gamma_)
}
        end where
      end if
    end if
  end if
end if

end subroutine getpthermal
!=============================================================================
subroutine getfluxforhllc(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw, idims
double precision :: w(ixI^S,nw),f(ixG^T,1:nwflux),tmp(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
logical :: transport
!-----------------------------------------------------------------------------
transport=.true.

if(iw==m0_+idims)then
   ! f_i[m_i]=v_i*m_i+p
   call getpthermal(w,x,ixI^L,ixO^L,tmp)
   f(ixO^S,iw)=tmp(ixO^S)
{#IFDEF ENERGY
else if(iw==e_)then
   ! f_i[e]=v_i*e+m_i/rho*p
   call getpthermal(w,x,ixI^L,ixO^L,tmp)
   f(ixO^S,iw)=w(ixO^S,m0_+idims)/w(ixO^S,rho_)*tmp(ixO^S)
}
{#IFDEF DUST
{^DS&else if(iw==rhod^DS_) then
    where(w(ixO^S,rhod^DS_)>minrhod)
        f(ixO^S,iw)=w(ixO^S,(rhod^DS_)+idims*^NDS);
    else where
        f(ixO^S,iw) = zero
    end where
    transport=.false.
\}

{^DS&{else if (iw==m^Cd^DS_) then
    ! use tmp for speeddust here
    where(w(ixO^S,rhod^DS_)>minrhod)
        tmp(ixO^S)=w(ixO^S,(rhod^DS_)+idims*^NDS)/w(ixO^S,rhod^DS_);
    else where
        tmp(ixO^S)=zero;
    end where
    f(ixO^S,iw)=w(ixO^S,iw)*tmp(ixO^S)
    transport=.false.
\}
\}

}
else
   f(ixO^S,iw)=zero
endif


end subroutine getfluxforhllc
!=============================================================================
subroutine getflux(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw, idims
double precision :: w(ixI^S,nw), f(ixG^T), tmp(ixG^T)
double precision, intent(in)    :: x(ixI^S,1:ndim)
logical :: transport
!-----------------------------------------------------------------------------
transport=.true.

if(iw==m0_+idims)then
   ! f_i[m_i]=v_i*m_i+p
   call getpthermal(w,x,ixI^L,ixO^L,f)
{#IFDEF ENERGY
else if(iw==e_)then
   ! f_i[e]=v_i*e+m_i/rho*p
   call getpthermal(w,x,ixI^L,ixO^L,f)
   f(ixO^S)=w(ixO^S,m0_+idims)/w(ixO^S,rho_)*f(ixO^S)
}
{#IFDEF DUST
{^DS&else if(iw==rhod^DS_) then
    where(w(ixO^S,rhod^DS_)>minrhod)
        f(ixO^S)=w(ixO^S,(rhod^DS_)+idims*^NDS);
    else where
        f(ixO^S) = zero
    end where
    transport=.false.
\}

{^DS&{else if (iw==m^Cd^DS_) then
    ! use tmp for speeddust here
    where(w(ixO^S,rhod^DS_)>minrhod)
        tmp(ixO^S)=w(ixO^S,(rhod^DS_)+idims*^NDS)/w(ixO^S,rhod^DS_);
    else where
        tmp(ixO^S)=zero;
    end where
    f(ixO^S)=w(ixO^S,iw)*tmp(ixO^S)
    transport=.false.
\}
\}

}
else
   f(ixO^S)=zero
endif


end subroutine getflux
!=============================================================================
subroutine addgeometry(qdt,ixI^L,ixO^L,wCT,w,x)

! Add geometrical source terms to w

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: qdt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

integer :: iw,idir, h1x^L{^NOONED, h2x^L}
double precision :: tmp(ixG^T){#IFDEF DUST ,speeddust(ixG^T,2:3,1:^NDS)}
logical :: angmomfix=.false.
!-----------------------------------------------------------------------------

select case (typeaxial)
case ('slab')
   ! No source terms in slab symmetry
case ('cylindrical')
{#IFDEF DUST
{^IFPHI
   {^DS&where(wCT(ixO^S,rhod^DS_)>minrhod)
     speeddust(ixO^S,3,^DS)=wCT(ixO^S,mphid^DS_)/wCT(ixO^S,rhod^DS_);
    else where
     speeddust(ixO^S,3,^DS)=zero;
   end where\}
}}
   do iw=1,nwflux
      select case(iw)
      ! s[mr]=(pthermal+mphi**2/rho)/radius
      case (mr_)
         call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
         tmp(ixO^S)=zero
{^IFPHI
         tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_)
      ! s[mphi]=(-mphi*mr/rho)/radius
      case (mphi_)
         tmp(ixO^S)=-wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_)
}
{#IFDEF DUST
{^IFPHI
{^DS&case (mrd^DS_)
        tmp(ixO^S)=speeddust(ixO^S,3,^DS)*wCT(ixO^S,mphid^DS_)
        ! s[mphi]=(-mphi*mr/rho)/radius
      case (mphid^DS_)
         tmp(ixO^S)=-speeddust(ixO^S,3,^DS)*wCT(ixO^S,mrd^DS_)
      \}
}}
      end select
      ! Divide by radius and add to w
      if (iw==mr_.or.iw==mphi_) then
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
      end if
{#IFDEF DUST
     {^DS& if ( iw==mrd^DS_.or.iw==mphid^DS_ ) then
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,r_)
      end if\}
}
   end do
case ('spherical')
   h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
{#IFDEF DUST
  {^NOONED
   {^DS&where(wCT(ixO^S,rhod^DS_)>minrhod)
     speeddust(ixO^S,2,^DS)=wCT(ixO^S,m2d^DS_)/wCT(ixO^S,rhod^DS_);
   else where
     speeddust(ixO^S,2,^DS)=zero;
   end where\}
   }
   {^IFTHREED
   {^DS&where(wCT(ixO^S,rhod^DS_)>minrhod)
     speeddust(ixO^S,3,^DS)=wCT(ixO^S,m3d^DS_)/wCT(ixO^S,rhod^DS_);
   else where
     speeddust(ixO^S,3,^DS)=zero;
   end where\}
   }
}
   do iw=1,nwflux
      select case (iw)
      ! s[m1]=((mtheta**2+mphi**2)/rho+2*p)/r
      case (m1_)
         call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
         ! For nonuniform Cartesian grid this provides hydrostatic equil.
         tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
                 *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
                 /mygeo%dvolume(ixO^S){&^CE&
               +wCT(ixO^S,m^CE_)**2/wCT(ixO^S,rho_) }
{^NOONEC
      ! s[m2]=-(mr*mtheta/rho)/r
      !       + cot(theta)*(mphi**2/rho+p)/r
      case (m2_)
}
{^NOONED
         call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
         ! This will make hydrostatic p=const an exact solution
         w(ixO^S,iw)=w(ixO^S,iw) &
            +qdt*tmp(ixO^S)*(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
            /mygeo%dvolume(ixO^S)
}
{^NOONEC
         tmp(ixO^S)=-(wCT(ixO^S,m1_)*wCT(ixO^S,m2_)/wCT(ixO^S,rho_))
}
{^IFTHREEC
{^NOONED
         tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,m3_)**2/wCT(ixO^S,rho_)) &
                        *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
}
      ! s[m3]=-(mphi*mr/rho)/r
      !       -cot(theta)*(mtheta*mphi/rho)/r
      case (m3_)
         if (.not.angmomfix) &
         tmp(ixO^S)=-(wCT(ixO^S,m3_)*wCT(ixO^S,m1_)/wCT(ixO^S,rho_)) {^NOONED &
                   -(wCT(ixO^S,m2_)*wCT(ixO^S,m3_)/wCT(ixO^S,rho_)) &
                   *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
}
{#IFDEF DUST
      ! s[md1]=((mdtheta**2+mdphi**2)/rho)/r
      {^DS&case (m1d^DS_)

         ! For nonuniform Cartesian grid this provides hydrostatic equil.
         tmp(ixO^S)=0{^CE&+wCT(ixO^S,m^CEd^DS_)*speeddust(ixO^S,^CE,^DS) }
       \}
{^NOONEC
      ! s[m2]=-(mr*mtheta/rho)/r
      !       + cot(theta)*(mphi**2/rho)/r

      {^DS&case (m2d^DS_)

         tmp(ixO^S)=-(wCT(ixO^S,m1d^DS_)*speeddust(ixO^S,2,^DS))
         if(ndim>1) then
{^IFTHREEC
             tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,m3d^DS_)*speeddust(ixO^S,3,^DS)) &
                        *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
}
            endif
      \}
}

{^IFTHREEC
      ! s[m3]=-(mphi*mr/rho)/r
      !       -cot(theta)*(mtheta*mphi/rho)/r
      {^DS&case (m3d^DS_)
         if (.not.angmomfix) &
         tmp(ixO^S)=-(speeddust(ixO^S,3,^DS)*wCT(ixO^S,m1d^DS_)) {^NOONED &
                   -(wCT(ixO^S,m2d^DS_)*speeddust(ixO^S,3,^DS)) &
                   *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
      \}
}
}
      end select
      ! Divide by radius and add to w
      if (iw==m1_{^NOONEC.or.iw==m2_}{^IFTHREEC  &
            .or.(iw==m3_.and..not.angmomfix)}) &
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
{#IFDEF DUST
      {^DS&if (iw==m1d^DS_{^NOONEC.or.iw==m2d^DS_}{^IFTHREEC  &
            .or.(iw==m3d^DS_.and..not.angmomfix)}) &
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
      \}
}
   end do
end select

if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"addgeometry")

end subroutine addgeometry
!=============================================================================
{#IFDEF DUST
subroutine addsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

! w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
logical, intent(in) :: qsourcesplit

! .. local ..
double precision, dimension(ixG^T,1:^NC,1:^NDS)  :: fdrag
integer :: idir
!-----------------------------------------------------------------------------

select case( TRIM(dustmethod) )
   case( 'none' ) 
	!do nothing here
   case default !all regular dust methods here
	if(qsourcesplit .eqv. ssplitdust) then
		call get_3d_dragforce(ixI^L,ixO^L,wCT,x,fdrag)
		do idir=1,ndir
		    fdrag(ixO^S,idir,1:^NDS)=fdrag(ixO^S,idir,1:^NDS)*qdt
		    {^DS&
		    !!RK!where( w(ixO^S,rhod^DS_)>minrhod )
		       w(ixO^S,m0_+idir)  = w(ixO^S,m0_+idir)  + fdrag(ixO^S,idir,^DS)
		{#IFDEF ENERGY
		       w(ixO^S,e_)        = w(ixO^S,e_)        + (wCT(ixO^S,m0_+idir)/ &
				            wCT(ixO^S,rho_))*fdrag(ixO^S,idir,^DS)
		}
		       w(ixO^S,(rhod^DS_)+idir*^NDS) = w(ixO^S,(rhod^DS_)+idir*^NDS) - fdrag(ixO^S,idir,^DS)
		    !!RK!end where
		    \}
		end do

		if( dustzero ) call set_dusttozero(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
	endif
end select

end subroutine addsource
!=============================================================================
subroutine set_dusttozero(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
!
!  Force dust density to zero if rho_dust <= minrhod
!
include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

{^DS&where( w(ixO^S,rhod^DS_)<=minrhod )
    w(ixO^S,rhod^DS_) = zero
    {^C&w(ixO^S,m^Cd^DS_) = zero;}
end where
\}

end subroutine set_dusttozero
!=============================================================================
subroutine get_3d_dragforce(ixI^L,ixO^L,w,x,fdrag)
!
! Calculate drag force based on Epstein's or Stokes' law
! From Kwok 1975,page 584 (between eqn 8 and 9)
!

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)

! local
double precision, dimension(ixG^T) :: vt2,deltav,fd,ptherm,vdust,vgas,Tgas
double precision, dimension(ixG^T,1:^NC,1:^NDS) :: fdrag
double precision, dimension(ixG^T,1:^NDS) :: alpha_T 

integer :: idir
double precision :: K
!-----------------------------------------------------------------------------

call getpthermal(w,x,ixI^L,ixO^L,ptherm)

vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S,rho_)

select case( TRIM(dustmethod) )
 case( 'Kwok' ) ! assume sticking coefficient equals 0.25

  do idir=1,^NC
    call getv(w,x,ixI^L,ixO^L,idir,vgas)

    {^DS&where(w(ixO^S,rhod^DS_)>minrhod)
       vdust(ixO^S)  = w(ixO^S,(rhod^DS_)+idir*^NDS)/w(ixO^S,rhod^DS_)
       deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))
       ! 0.75 from sticking coefficient
       fd(ixO^S)     = 0.75d0*w(ixO^S,rhod^DS_)*w(ixO^S,rho_)*deltav(ixO^S) &
                     / (rhodust(^DS)*sdust(^DS))
       ! 0.75 from spherical grainvolume
       fd(ixO^S)     = -fd(ixO^S)*0.75d0*dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)
    else where
       fd(ixO^S) = zero
    end where
    fdrag(ixO^S,idir,^DS) = fd(ixO^S)
    \}
  enddo
{#IFDEF ENERGY
 case( 'sticking' ) ! Calculate sticking coefficient based on the gas and dust temperatures
!
!  Equation from Decin et al. 2006
!
  Tgas(ixO^S) = ( ptherm(ixO^S)*normvar(p_)*mhcgspar)/(w(ixO^S,rho_)*normvar(rho_)*kbcgspar)
  call get_tdust(w,x,ixI^L,ixO^L,alpha_T)
  
  do idir=1,^NC
    call getv(w,x,ixI^L,ixO^L,idir,vgas)

    {^DS&where(w(ixO^S,rhod^DS_)>minrhod)
       alpha_T(ixO^S,^DS) = max(0.35d0*dexp(-dsqrt((Tgas(ixO^S)+alpha_T(ixO^S,^DS))/5.0d2))+0.1d0,smalldouble)
       vdust(ixO^S)  = w(ixO^S,(rhod^DS_)+idir*^NDS)/w(ixO^S,rhod^DS_)
       deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))
       fd(ixO^S)     = (one-alpha_T(ixO^S,^DS))*w(ixO^S,rhod^DS_)*w(ixO^S,rho_)* &
                        deltav(ixO^S) / (rhodust(^DS)*sdust(^DS))
       fd(ixO^S)     = -fd(ixO^S)*0.75d0*dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)       
    else where
       fd(ixO^S) = zero
    end where
    fdrag(ixO^S,idir,^DS) = fd(ixO^S)
    \}
  enddo 
}
 case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)   
  K = 3.4d5/^NDS
  do idir=1,^NC
    call getv(w,x,ixI^L,ixO^L,idir,vgas)

    {^DS&where(w(ixO^S,rhod^DS_)>minrhod)
       vdust(ixO^S)  = w(ixO^S,(rhod^DS_)+idir*^NDS)/w(ixO^S,rhod^DS_)
       deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))
       
       fd(ixO^S)     = -K*deltav(ixO^S)
    else where
       fd(ixO^S) = zero
    end where
    fdrag(ixO^S,idir,^DS) = fd(ixO^S)
    \}
  enddo
 case('none')  
    {^DS&fdrag(ixO^S,idir,^DS) = zero\}

 case default
      call mpistop( "===This dust method has not been implemented===" )
end select

return
end subroutine get_3d_dragforce
!=============================================================================

subroutine get_tdust(w,x,ixI^L,ixO^L,Td)
!
! Returns dust temperature (in K), either as constant or based 
! on equ. 5.41,5.42 and 5.44 from Tielens (2005)
!
!

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(out)   :: Td(ixG^T,1:^NDS)

double precision :: G0(ixO^S)

!---------------------------------------------------------------------------

select case( TRIM(dusttemp) )
case( 'constant' )
  Td(ixO^S,1:^NDS) = Tdust
case( 'ism' )
{^DS&select case( TRIM(dustspecies) )
  case( 'graphite' )
    Td(ixO^S,^DS) = 15.8d0*((1.0d0/(sdust(^DS)*normvar(0)))**0.06d0)
  case( 'silicate' )
    Td(ixO^S,^DS) = 13.6d0*((1.0d0/(sdust(^DS)*normvar(0)))**0.06d0)
  case default
      call mpistop( "===Dust species undetermined===" )
  end select
  \}  
case( 'stellar' )
  select case( TRIM(typeaxial) )
  case( 'spherical' )
    G0(ixO^S) = max(x(ixO^S,1)*normvar(0),smalldouble)
  case( 'cylindrical' )
    G0(ixO^S) = max(dsqrt(x(ixO^S,1)**2 + x(ixO^S,2)**2)*normvar(0),smalldouble)
  case( 'slab' )
{^IFTHREED
    G0(ixO^S) = max(dsqrt((x(ixO^S,1)-xptmass1)**2 + (x(ixO^S,2)-xptmass2)**2  & 
                + (x(ixO^S,3)-xptmass3)**2)*normvar(0),smalldouble)
}
  end select
  G0(ixO^S)= 2.1d4*(Lstar/1.0d8)*((3.0857d17/G0(ixO^S))**2)
{^DS&select case( TRIM(dustspecies) )
  case( 'graphite' )
    Td(ixO^S,^DS) = 61.0d0*((1.0d0/(sdust(^DS)*normvar(0)))**0.06d0) &
              *(G0(ixO^S)**(one/5.8d0))
  case( 'silicate' )
    Td(ixO^S,^DS) = 50.0d0*((1.0d0/(sdust(^DS)*normvar(0)))**0.06d0) &
              *(G0(ixO^S)**(one/6.0d0))
  case default
      call mpistop( "===Dust species undetermined===" )
  end select
  \}
case default
    call mpistop( "===Dust temperature undetermined===" )
end select    
   
return
end subroutine get_tdust
!=============================================================================

subroutine getdt(w,ixI^L,ixO^L,dtnew,dx^D,x)

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
! .. local ..
integer                                             :: idims,idust,idir
double precision                                    :: dxinv(^ND)
double precision, dimension(ixG^T,1:^ND)            :: agas
double precision, dimension(ixG^T,1:^ND,1:^NDS)     :: adust
double precision                                    :: dtgas
double precision, dimension(1:^NDS)                 :: dtdust
double precision, dimension(ixG^T,1:^NC,1:^NDS)     :: fdrag
double precision, dimension(ixG^T) :: vt2,deltav,tstop,ptherm,vdust,vgas
double precision :: K
!-----------------------------------------------------------------------------

dtnew=bigdouble
dtgas=bigdouble
dtdust=bigdouble


!-----------------------------
!get dt related to dust and gas stopping time (Laibe 2011)
!-----------------------------


select case( TRIM(dustmethod) )
  
  case( 'Kwok' ) ! assume sticking coefficient equals 0.25
    call mpistop( "===This dust method getdt has not been implemented===" )
    
  
  case( 'sticking' ) ! Calculate sticking coefficient based on the gas temperature
    call get_3d_dragforce(ixI^L,ixO^L,w,x,fdrag)
    dtdust(1:^NDS) = bigdouble

    call getpthermal(w,x,ixI^L,ixO^L,ptherm)
    vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S,rho_)


    !Tgas, mu=mean molecular weight
      ptherm(ixO^S) = ( ptherm(ixO^S)*normvar(p_)*mhcgspar*eqpar(mu_))/(w(ixO^S,rho_)*normvar(rho_)*kbcgspar)
      ptherm(ixO^S) = max(0.35d0*dexp(-dsqrt(ptherm(ixO^S)*2.d-3))+0.1d0,smalldouble)

    do idir=1,^NC
    call getv(w,x,ixI^L,ixO^L,idir,vgas)

    {^DS&where(w(ixO^S,rhod^DS_)>minrhod)
       vdust(ixO^S)  = w(ixO^S,(rhod^DS_)+idir*^NDS)/w(ixO^S,rhod^DS_)
       deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))
       tstop(ixO^S)  = 4.0d0*(rhodust(^DS)*sdust(^DS))/ &
                     (3.0d0*(one-ptherm(ixO^S))*dsqrt(vt2(ixO^S) + &
                     deltav(ixO^S)**2)*(w(ixO^S,rhod^DS_) + &
                     w(ixO^S,rho_)))
    else where
       tstop(ixO^S) = bigdouble
    end where


    dtdust(^DS) = min(minval(tstop(ixO^S)),dtdust(^DS))
    \}
    enddo

    dtnew = min(minval(dtdiffpar*dtdust(1:^NDS)),dtnew)
 
 
 
 case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)   
    K = 3.4d5/^NDS
    dtdust(1:^NDS) = bigdouble

    {^DS&where(w(ixO^S,rhod^DS_)>minrhod)
       tstop(ixO^S)  = (w(ixO^S,rhod^DS_)*w(ixO^S,rho_))/ &
                     (K*(w(ixO^S,rhod^DS_) + w(ixO^S,rho_)))
    else where
       tstop(ixO^S) = bigdouble
    end where


    dtdust(^DS) = min(minval(tstop(ixO^S)),dtdust(^DS))
    \} 
    
    dtnew = min(minval(dtdiffpar*dtdust(1:^NDS)),dtnew)
 case('none')
    ! no dust timestep
 case default
      call mpistop( "===This dust method has not been implemented===" )
end select


if(dtnew<dtmin)then
   write(unitterm,*)"-------------------------------------"
   write(unitterm,*)"Warning: found DUST related time step too small! dtnew=",dtnew
   write(unitterm,*)"on grid with index:", saveigrid," grid level=",node(plevel_,saveigrid)
   write(unitterm,*)"grid corners are=",{^D&rnode(rpxmin^D_,saveigrid), rnode(rpxmax^D_,saveigrid)}
   write(unitterm,*)"dtgas=",dtgas, " dtdust =",dtdust(1:^NDS)
   write(unitterm,*)"on processor:", mype
   write(unitterm,*)"-------------------------------------"
endif

end subroutine getdt
}
!=============================================================================
! end module amrvacphys/- hd
!##############################################################################
