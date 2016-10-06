!=============================================================================
! amrvacusr.t.testrho

! INCLUDE:amrvacnul.specialini.t
INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!----------------------------------------------------------------------------

{^IFONED   eqpar(v1_)=one }
{^IFTWOD   call mpistop("just a 1D test") }
{^IFTHREED call mpistop("just a 1D test") }

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
!----------------------------------------------------------------------------
select case (iprob)
case (1)
   ! Advection of \sin(\pi x)**4
   w(ix^S,rho_) = sin(dpi*x(ix^S,1))**4.0d0
case default
   call mpistop("iprob not available!")
end select

end subroutine initonegrid_usr
!=============================================================================
! amrvacusr.t.testrho
!=============================================================================
