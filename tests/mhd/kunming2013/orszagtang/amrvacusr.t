!=============================================================================
! amrvacusr.t.mhdOT

! INCLUDE:amrvacnul/specialini.t
! INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------
eqpar(gamma_)=5.0d0/3.0d0
eqpar(eta_)=zero
{#IFDEF ISO
eqpar(gamma_)=1.0d0
eqpar(adiab_)=1.641192d0
}
end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: rho0,p0,b0

logical, save :: first=.true.
logical:: patchw(ixG^T)
!----------------------------------------------------------------------------

rho0=25.0d0/9.0d0
p0=5.0d0/3.0d0
b0=1.0d0

w(ix^S,rho_)=rho0
w(ix^S,v1_)=-sin(x(ix^S,2)){^IFTHREED*dcos(x(ixG^S,3))}
w(ix^S,v2_)= sin(x(ix^S,1)){^IFTHREED*dcos(x(ixG^S,3))}
{^IFTHREED
w(ix^S,v3_)= zero
}
{#IFDEF ENERGY
w(ix^S,pp_)=p0
}
w(ix^S,b1_)=-b0*sin(x(ix^S,2)){^IFTHREED*dcos(x(ixG^S,3))}
w(ix^S,b2_)= b0*sin(two*x(ix^S,1)){^IFTHREED*dcos(x(ixG^S,3))}
{^IFTHREED
w(ix^S,b3_)= zero
}

patchw(ix^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

if(first)then
      if(mype==0)then
       write(*,*)'Doing ideal MHD, Orszag Tang problem'
       write(*,*)'rho - p - b - gamma?:',rho0,p0,b0,eqpar(gamma_)
      endif
      first=.false.
endif


return
end subroutine initonegrid_usr
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
! .. local ..
double precision :: current(ixG^T,7-2*ndir:3)
integer          :: idirmin
double precision :: wloc(ixI^S,1:nw)
double precision :: pth(ixG^T), rho(ixG^T)
!-----------------------------------------------------------------------------

wloc(ixI^S,1:nw)=w(ixI^S,1:nw)
if(saveprim)then
  {#IFDEF ENERGY
  pth(ixO^S)=wloc(ixO^S,p_)
  }
  {#IFDEF ISO
  pth(ixO^S)=eqpar(adiab_)*wloc(ixO^S,rho_)**eqpar(gamma_)
  }
else
  call getpthermal(wloc,x,ixI^L,ixO^L,pth)
endif
rho(ixO^S)=wloc(ixO^S,rho_)
w(ixO^S,nw+1)=pth(ixO^S)/rho(ixO^S)


call getcurrent(wloc,ixI^L,ixO^L,idirmin,current)
w(ixO^S,nw+2)=current(ixO^S,3)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'T'
primnames= TRIM(primnames)//' '//'jz'
wnames=TRIM(wnames)//' '//'T'
wnames=TRIM(wnames)//' '//'jz'


end subroutine specialvarnames_output
!=============================================================================
subroutine printlog_special

use mod_global_parameters
!-----------------------------------------------------------------------------

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

use mod_global_parameters

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
!=============================================================================
! amrvacusr.t.mhdOT
!=============================================================================
