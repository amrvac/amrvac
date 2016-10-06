!=============================================================================
! amrvacusr.t.nul
!=============================================================================
INCLUDE:amrvacnul/speciallog.t
!INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------
{#IFDEF GLM
eqpar(Cr_)    = 0.18d0
}

{#IFDEF SYNGE
eqpar(gamma_) = 5.0d0/3.0d0
}{#IFNDEF SYNGE
eqpar(gamma_) = 4.0d0/3.0d0
}

eqpar(bm_) = one
eqpar(rm_) = 0.37d0
eqpar(rj_) = one
eqpar(lor_) = 10.0d0
eqpar(rhoe_) = 1.0d3
eqpar(betam_) = 0.34d0


end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

! initialize one grid within ixO^L

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)
logical patchw(ixG^T)
double precision                   :: pe
double precision, dimension(ixG^T) :: r
!-----------------------------------------------------------------------------
patchw(ixO^S) = .false.

!------------------------------ Geometry: ------------------------------
r(ixI^S) = abs(x(ixI^S,1))
!-----------------------------------------------------------------------

{#IFDEF GLM
w(ixO^S,psi_) = zero
}

pe = eqpar(betam_)*eqpar(bm_)**2/2.0d0

w(ixO^S,b3_) = zero
w(ixO^S,rho_)= eqpar(rhoe_)
w(ixO^S,pp_) = pe
w(ixO^S,v2_) = zero


w(ixO^S,lfac_)=one/dsqrt(one-({^C&w(ixO^S,v^C_)**2.0d0+}))
if (useprimitiveRel)  then 
   {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,v^C_);\}
end if
call conserve(ixI^L,ixO^L,w,x,patchw)
end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixI^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined
! user must assign conservative variables in bounderies

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw, iB
double precision, intent(in) :: qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)
logical patchw(ixG^T)
double precision                   :: pe, alpha
double precision, dimension(ixG^T) :: r
!-----------------------------------------------------------------------------
select case(iB)
 case(3)
!   ! min boundary in the 2nd dimension 
patchw(ixO^S) = .false.

!------------------------------ Geometry: ------------------------------
r(ixI^S) = abs(x(ixI^S,1))
!-----------------------------------------------------------------------

{#IFDEF GLM
w(ixO^S,psi_) = zero
}

{^C& w(ixO^S,v^C_)=zero\}
{^C& w(ixO^S,b^C_)=zero\}

pe = eqpar(betam_)*eqpar(bm_)**2/2.0d0
alpha = one - one/eqpar(betam_)*(eqpar(rm_)/eqpar(rj_))**2

where (r(ixO^S) .lt. eqpar(rm_))
   w(ixO^S,b3_) = - eqpar(lor_) * eqpar(bm_) * (r(ixO^S)/eqpar(rm_))
   w(ixO^S,rho_) = one
   w(ixO^S,pp_) = pe * (alpha + 2.0d0/eqpar(betam_)*(one - (r(ixO^S)/eqpar(rm_))**2))
   w(ixO^S,v2_) = sqrt(one - one/eqpar(lor_)**2)
elsewhere(eqpar(rm_) .le. r(ixO^S) .and. r(ixO^S) .lt. eqpar(rj_))
   w(ixO^S,b3_) = - eqpar(lor_) * eqpar(bm_) * (eqpar(rm_)/r(ixO^S))
   w(ixO^S,rho_) = one
   w(ixO^S,pp_) = alpha * pe
   w(ixO^S,v2_) = sqrt(one - one/eqpar(lor_)**2)
elsewhere
   w(ixO^S,b3_) = zero
   w(ixO^S,rho_)= eqpar(rhoe_)
   w(ixO^S,pp_) = pe
   w(ixO^S,v2_) = zero
end where


w(ixO^S,lfac_)=one/dsqrt(one-({^C&w(ixO^S,v^C_)**2.0d0+}))
if (useprimitiveRel)  then 
   {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,v^C_);\}
end if
call conserve(ixI^L,ixO^L,w,x,patchw)
end select

end subroutine specialbound_usr
!=============================================================================






























!=============================================================================
subroutine bc_int(level,qt,ixI^L,ixO^L,w,x)

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

integer, intent(in) :: ixI^L,ixO^L,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(in) :: x(ixI^S,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

! just to give an example for relativistic MHD
!  -----------------------------------------
!patchw(ixO^S)=.true.
!where (({^D&x(ixO^S,^D)**2+})<half**2.0d0) 
!    patchw(ixO^S) = .false.
!  ^C&w(ixO^S,v^C_)=zero;
!  ^C&w(ixO^S,b^C_)=zero;
!    w(ixO^S,b3_) = one
!    w(ixO^S,v1_) = 0.99
!    w(ixO^S,rho_) = 1.d0
!    w(ixO^S,pp_)  = 2.0d0
!    w(ixO^S,lfac_)=one/dsqrt(one-({^C&w(ixO^S,v^C_)**2.0d0+}))
!end where
!!if (useprimitiveRel) then
!!  where (({^D&x(ixO^S,^D)**2+})<half**2.0d0) 
!!  {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,v^C_);\}
!!  end where
!!endif
!call conserve(ixI^L,ixO^L,w,patchw)

end subroutine bc_int
!=============================================================================
! amrvacusr.t.nul
!=============================================================================
