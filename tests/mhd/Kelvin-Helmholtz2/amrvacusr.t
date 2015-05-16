!=============================================================================
! amrvacusr.t.KH

! INCLUDE:amrvacnul/specialini.t
! INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
eqpar(gamma_)=5.0d0/3.0d0
{#IFDEF ISO
eqpar(adiab_)=0.0d0
}
{#IFDEF GLM
eqpar(Cr_)=0.18d0
}
end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

! .. local ..
double precision, parameter :: y0=1./20., M=1., theta=dpi/3., ca=0.1, rho0=1., vy0=1.0d-2, sigma=0.1
double precision            :: p0

logical:: patchw(ixG^T)
!----------------------------------------------------------------------------
p0 = 1.0d0/eqpar(gamma_)

patchw(ix^S)  = .false.

w(ix^S,v1_)  = M/2. * tanh(x(ix^S,2)/y0)
! Perturbation in vy:
w(ix^S,v2_)  = vy0 * sin(2.0d0*dpi*x(ix^S,1)) * exp(-(x(ix^S,2)/sigma)**2.0d0)

w(ix^S,b1_)  = ca * sqrt(rho0) * cos(theta)
w(ix^S,b3_)  = ca * sqrt(rho0) * sin(theta)

w(ix^S,rho_) = rho0
w(ix^S,p_)   = p0

call conserven(ixG^L,ix^L,w,patchw)

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
