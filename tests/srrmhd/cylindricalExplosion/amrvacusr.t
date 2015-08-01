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

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

eqpar(gamma_) = 4.0d0/3.0d0
eqpar(eta_)   = 0.018d0
eqpar(kappa_) = 1.0d0/0.18d0

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid within ix^L

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
! .. local ..
double precision, dimension(ixG^T)  :: r
!-----------------------------------------------------------------------------

! =====================
! Geometry:
r(ix^S) = sqrt({^D& x(ix^S,^D)**2|+})
! =====================

w(ix^S,psi_)  = 0.0d0
w(ix^S,phib_) = 0.0d0
w(ix^S,q_)    = 0.0d0

w(ix^S,u1_)  = 0.0d0
w(ix^S,u2_)  = 0.0d0
w(ix^S,u3_)  = 0.0d0 

w(ix^S,b1_)  = 0.1d0
w(ix^S,b2_)  = 0.0d0
w(ix^S,b3_)  = 0.0d0

where (r(ix^S) .lt. 0.8d0)
   w(ix^S,rho_) = 0.01d0
   w(ix^S,pp_)  = 1.0d0
elsewhere (r(ix^S) .ge. 0.8d0 .and. r(ix^S) .lt. one)
   w(ix^S,rho_) = 0.01d0 * exp(-(r(ix^S)-0.8d0))
   w(ix^S,pp_)  = 1.0d0 * exp(-(r(ix^S)-0.8d0))
elsewhere
   w(ix^S,rho_) = 0.001d0
   w(ix^S,pp_)  = 0.001d0
end where

w(ix^S,lfac_)= sqrt(one + {^C& w(ix^S,u^C_)**2|+})
w(ix^S,e1_)  = w(ix^S,b2_)*w(ix^S,u3_) - w(ix^S,b3_)*w(ix^S,u2_)
w(ix^S,e2_)  = w(ix^S,b3_)*w(ix^S,u1_) - w(ix^S,b1_)*w(ix^S,u3_)
w(ix^S,e3_)  = w(ix^S,b1_)*w(ix^S,u2_) - w(ix^S,b2_)*w(ix^S,u1_)
{^C& w(ix^S,e^C_)  = w(ix^S,e^C_)/w(ix^S,lfac_)\}

call conserve(ixG^L,ix^L,w,x,patchfalse)

end subroutine initonegrid_usr
!=============================================================================
{#IFDEF FCT
subroutine initvecpot_usr(ixI^L, ixC^L, xC, A)

! initialize the vectorpotential on the corners
! used by b_from_vectorpotential()


include 'amrvacdef.f'

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
