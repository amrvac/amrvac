!=============================================================================
! amrvacusr.t.srmhdOT

! INCLUDE:amrvacnul/specialini.t
! INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
INCLUDE:amrvacmodules/handle_particles.t
INCLUDE:amrvacmodules/integrate_particles.t
!=============================================================================
subroutine initglobaldata_usr

use constants
include 'amrvacdef.f'
!-----------------------------------------------------------------------------

UNIT_LENGTH = 1.0d18
UNIT_VELOCITY = CONST_c
UNIT_DENSITY = 1.0d-36

normvar(0) = UNIT_LENGTH
normvar(rho_) = UNIT_DENSITY
{^C&normvar(v^C_)   = UNIT_VELOCITY \}
normvar(pp_)     = UNIT_VELOCITY**2.0d0 * UNIT_DENSITY
{^C&normvar(b^C_)   = sqrt(4.0d0*dpi*normvar(pp_)) \}
normvar(xi_)     = UNIT_DENSITY * CONST_c**2
normt = UNIT_LENGTH/UNIT_VELOCITY

eqpar(gamma_)=4.0d0/3.0d0
end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: rho0,p0,vmax

logical, save :: first=.true.
logical:: patchw(ixG^T)
!----------------------------------------------------------------------------
oktest = index(teststr,'initonegrid_usr')>=1
if (oktest) write(unitterm,*) ' === initonegrid_usr  (in ) : ', &
      'ixG^L : ',ixG^L

rho0=one
p0=10.0d0
vmax=0.99d0

vmax=vmax/dsqrt(two)
w(ix^S,rho_)=rho0
w(ix^S,v1_)=-vmax*sin(x(ix^S,2))
w(ix^S,v2_)= vmax*sin(x(ix^S,1))
w(ix^S,pp_)=p0
w(ix^S,b1_)=-sin(x(ix^S,2))
w(ix^S,b2_)= sin(two*x(ix^S,1))

{#IFDEF C3
w(ix^S,v3_)= zero
w(ix^S,b3_)= one
}

{#IFDEF EPSINF
w(ix^S,epsinf_)=one
w(ix^S,rho0_)=one
w(ix^S,rho1_)=one
w(ix^S,n_)=one
w(ix^S,n0_)=one
}

{#IFDEF GLM
w(ix^S,psi_)=zero
}

  
w(ix^S,lfac_)=one/dsqrt(one-({^C&w(ix^S,v^C_)**2+}))
if(useprimitiveRel)then
   {^C&w(ix^S,u^C_)=w(ix^S,lfac_)*w(ix^S,v^C_)\}
endif
  
patchw(ix^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

if(first)then
      write(*,*)'Doing 2D ideal SRMHD, Orszag Tang problem'
      write(*,*)'rho - p - gamma - primRel?:',rho0,p0,eqpar(gamma_),useprimitiveRel
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

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)
! .. local ..
double precision:: divb(ixG^T)
double precision :: current(ixG^T,7-2*ndir:3)
integer          :: idirmin
!-----------------------------------------------------------------------------

!call mpistop("special output file undefined")
call getdivb(w,ixI^L,ixO^L,divb)
w(ixO^S,nw+1)=divb(ixO^S)
call getcurrent(w,ixI^L,ixO^L,idirmin,current)
w(ixO^S,nw+2)=current(ixO^S,3)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

!call mpistop("special varnames and primnames undefined")
primnames= TRIM(primnames)//' '//'divb'
primnames= TRIM(primnames)//' '//'jz'
wnames=TRIM(wnames)//' '//'divb'
wnames=TRIM(wnames)//' '//'jz'

end subroutine specialvarnames_output
!=============================================================================
subroutine printlog_special

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

include 'amrvacdef.f'

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
!=============================================================================
! amrvacusr.t.srmhdOT
!=============================================================================
