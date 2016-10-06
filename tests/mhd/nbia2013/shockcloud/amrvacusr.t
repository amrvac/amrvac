!=============================================================================
! amrvacusr.t.shockcloud
!=============================================================================
!INCLUDE:amrvacnul/specialini.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/speciallog.t
!INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------
eqpar(gamma_)=1.66666667d0
eqpar(eta_)=zero
eqpar(etah_)=zero

eqpar(chi_)=10.0d0
eqpar(nval_)=8.0d0
eqpar(rcore_)=0.62d0
eqpar(rbound_)=1.77d0
eqpar(xshock_)=-2.66d0
eqpar(machs_)=10.d0

select case(iprob)
  case(1,11)
    eqpar(beta_)=1.0d0
  case(2,12)
    eqpar(beta_)=10.0d0
  case(3,13)
    eqpar(beta_)=0.5d0
  case default
    call mpistop('iprob to implement')
endselect

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: rval(ixG^T),pleft,rholeft,vleft,Prat,alfa,cleft
double precision:: pright,rhoright,vright,deltax

logical :: patchw(ixG^T)
logical, save :: first=.true.
!----------------------------------------------------------------------------

if(typephys/='mhd')call mpistop("this is an MHD problem: set typephys!")


! uniform cloud surroundings
pleft=one
rholeft=one+(eqpar(chi_)-one)/(one+(eqpar(rbound_)/eqpar(rcore_))**eqpar(nval_))
vleft=zero

! set uniform parallel B field
w(ix^S,b1_)=dsqrt(two*pleft/eqpar(beta_))
w(ix^S,b2_)=zero
{^IFTHREED
w(ix^S,b3_)=zero
}

! set trivial velocity components
w(ix^S,v2_)=zero
{^IFTHREED
w(ix^S,v3_)=zero
}

! compute the RH related states
Prat=one/(one+(eqpar(machs_)**2-one)*two*eqpar(gamma_)/(eqpar(gamma_)+one))
alfa=(eqpar(gamma_)+1)/(eqpar(gamma_)-one)
cleft=dsqrt(eqpar(gamma_)*pleft/rholeft)
rhoright=rholeft*(alfa+Prat)/(alfa*Prat+one)
pright=pleft/Prat
vright=cleft*eqpar(machs_)*(one-(alfa*Prat+one)/(alfa+Prat))

deltax=2.0d-2
select case(iprob)
  case(11,12,13)
    w(ix^S,rho_)=rhoright+(rholeft-rhoright) &
            *half*(dtanh((x(ix^S,1)-eqpar(xshock_))/deltax)+one)
    w(ix^S,v1_) =vright  +(vleft-vright)     &
            *half*(dtanh((x(ix^S,1)-eqpar(xshock_))/deltax)+one)
    w(ix^S,p_)  =pright  +(pleft-pright)     &
            *half*(dtanh((x(ix^S,1)-eqpar(xshock_))/deltax)+one)
  case(1,2,3)
    where(x(ix^S,1)<eqpar(xshock_))
      w(ix^S,rho_)=rhoright
      w(ix^S,v1_) =vright
      w(ix^S,p_)  =pright
    elsewhere
      w(ix^S,rho_)=rholeft
      w(ix^S,v1_) =vleft
      w(ix^S,p_)  =pleft
    endwhere
  case default
    call mpistop('iprob to implement')
end select

! overwrite cloud region density, set tracer values
rval(ix^S)=dsqrt(^D&x(ix^S,^D)**2+)
where(rval(ix^S)<eqpar(rbound_))
  w(ix^S,rho_)=one+(eqpar(chi_)-one)/(one+(rval(ix^S)/eqpar(rcore_))**eqpar(nval_))
  w(ix^S,tr1_)=one
elsewhere
  w(ix^S,tr1_)=zero
endwhere

if(mype==0.and.first)then
   write(*,*)'Doing shock-cloud challenge, ideal MHD'
   write(*,*)'density is=',rholeft,' pressure is=',pleft
   write(*,*)'plasma beta is set to=',eqpar(beta_)
   write(*,*)'shock Mach is=',eqpar(machs_)
   first=.false.
endif

patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined
! user must assign conservative variables in bounderies

use mod_global_parameters

integer, intent(in) :: ixG^L, ixO^L, iw, iB
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: pleft,rholeft,Prat,alfa,cleft
!----------------------------------------------------------------------------

pleft=one
rholeft=one+(eqpar(chi_)-one)/(one+(eqpar(rbound_)/eqpar(rcore_))**eqpar(nval_))
Prat=one/(one+(eqpar(machs_)**2-one)*two*eqpar(gamma_)/(eqpar(gamma_)+one))
alfa=(eqpar(gamma_)+1)/(eqpar(gamma_)-one)
cleft=dsqrt(eqpar(gamma_)*pleft/rholeft)

select case(iB)
 case(1)
  ! fix the postshock values
  w(ixO^S,rho_)=rholeft*(alfa+Prat)/(alfa*Prat+one)
  w(ixO^S,p_)=pleft/Prat
  w(ixO^S,b1_)=dsqrt(two*pleft/eqpar(beta_))
  w(ixO^S,b2_)=zero
  {^IFTHREED
  w(ixO^S,b3_)=zero
  }
  w(ixO^S,v1_)=cleft*eqpar(machs_)*(one-(alfa*Prat+one)/(alfa+Prat))
  w(ixO^S,v2_)=zero
  {^IFTHREED
  w(ixO^S,v3_)=zero
  }
  w(ixO^S,tr1_)=zero
  call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case default
   call mpistop('boundary not defined')
end select

end subroutine specialbound_usr
!=============================================================================
subroutine bc_int(level,qt,ixG^L,ixO^L,w,x)

! internal boundary, user defined
!
! This subroutine can be used to artificially overwrite ALL conservative 
! variables in a user-selected region of the mesh, and thereby act as
! an internal boundary region. It is called just before external (ghost cell)
! boundary regions will be set by the BC selection. Here, you could e.g. 
! want to introduce an extra variable (nwextra, to be distinguished from nwaux)
! which can be used to identify the internal boundary region location.
! Its effect should always be local as it acts on the mesh.
!

use mod_global_parameters

integer, intent(in) :: ixG^L,ixO^L,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

end subroutine bc_int
!=============================================================================
! amrvacusr.t.shockcloud
!=============================================================================
