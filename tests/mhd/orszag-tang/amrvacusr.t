!=============================================================================
! amrvacusr.t.srmhdOT

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
{#IFDEF ISO
eqpar(adiab_)=1.641192d0
}
{#IFDEF GLM
eqpar(Cr_) = 0.18d0
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

rho0=25.0d0/36.0d0/dpi
p0=5.0d0/12.0d0/dpi
b0=1.0d0/sqrt(4.0d0*dpi)

w(ix^S,rho_)=rho0
w(ix^S,v1_)=-sin(2.0d0*dpi*x(ix^S,2))
w(ix^S,v2_)= sin(2.0d0*dpi*x(ix^S,1))
{#IFDEF ENERGY
w(ix^S,pp_)=p0
}
w(ix^S,b1_)=-b0*sin(2.0d0*dpi*x(ix^S,2))
w(ix^S,b2_)= b0*sin(4.0d0*dpi*x(ix^S,1))

{#IFDEF C3
w(ix^S,b3_)= b0
w(ix^S,v3_)= 0.0d0
}

patchw(ix^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

if(first)then
      write(*,*)'Doing 2D ideal MHD, Orszag Tang problem'
      write(*,*)'rho - p - b - gamma?:',rho0,p0,b0,eqpar(gamma_)
      first=.false.
endif


return
end subroutine initonegrid_usr
!=============================================================================
{#IFDEF FCT
subroutine initvecpot_usr(ixI^L, ixC^L, xC, A)

! initialize the vectorpotential on the corners
! used by b_from_vectorpotential()


use mod_global_parameters

integer, intent(in)                :: ixI^L, ixC^L
double precision, intent(in)       :: xC(ixI^S,1:ndim)
double precision, intent(out)      :: A(ixI^S,1:ndir)

!-----------------------------------------------------------------------------

A=zero

end subroutine initvecpot_usr
!=============================================================================
}
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
double precision:: divb(ixG^T)
double precision :: current(ixG^T,7-2*ndir:3)
integer          :: idirmin
!-----------------------------------------------------------------------------

call getdivb(w,ixI^L,ixO^L,divb)
w(ixO^S,nw+1)=divb(ixO^S)

call getcurrent(w,ixI^L,ixO^L,idirmin,current)
w(ixO^S,nw+2)=current(ixO^S,3)
end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

!call mpistop("special varnames and primnames undefined")
primnames= TRIM(primnames)//' '//'divb'
primnames= TRIM(primnames)//' '//'jz'
wnames=TRIM(wnames)//' '//'divb'
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
! amrvacusr.t.srmhdOT
!=============================================================================
