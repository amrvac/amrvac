!=============================================================================
! amrvacusr.t.nul
!=============================================================================
INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------

eqpar(gamma_) = 2.0d0

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid within ix^L

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
logical patchw(ixG^T)
patchw(ix^S) = .false.
!-----------------------------------------------------------------------------

{#IFDEF PSI
   w(ix^S,psi_) = 0.0d0
}

where(x(ix^S,1) .lt. 0.0d0) 
   w(ix^S,rho_) = 1.0d0
   w(ix^S,pp_)  = 1.0d0

   w(ix^S,u1_)  = 0.0d0
   w(ix^S,u2_)  = 0.0d0
   w(ix^S,u3_)  = 0.0d0 

   w(ix^S,b1_)  = 0.5d0
   w(ix^S,b2_)  = 1.0d0
   w(ix^S,b3_)  = 0.0d0
elsewhere
   w(ix^S,rho_) = 0.125d0
   w(ix^S,pp_)  = 0.1d0

   w(ix^S,u1_)  = 0.0d0
   w(ix^S,u2_)  = 0.0d0
   w(ix^S,u3_)  = 0.0d0 

   w(ix^S,b1_)  = 0.5d0
   w(ix^S,b2_)  =-1.0d0
   w(ix^S,b3_)  = 0.0d0
   
end where


call conserve(ixG^L,ix^L,w,x,patchw)
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

! Not needed for this setup.

A(ixC^S,1:ndir) = zero


end subroutine initvecpot_usr
}
!=============================================================================
! amrvacusr.t.nul
!=============================================================================
