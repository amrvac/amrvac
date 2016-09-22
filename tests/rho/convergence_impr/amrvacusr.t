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
   call test_solution_1(ix^L, x(ix^S, :), 0.0d0, w(ix^S,rho_))
case default
   call mpistop("iprob not available!")
end select

end subroutine initonegrid_usr

subroutine test_solution_1(ix^L, x, time, val)
  use mod_global_parameters
  integer, intent(in)          :: ix^L
  double precision, intent(in) :: x(ix^S, ndim), time
  double precision             :: xrel(ix^S, ndim), val(ix^S)
  double precision             :: vel(ndim)
  integer                      :: idim

  vel = eqpar(v1_:v1_+ndim-1)
  val = 1

  do idim = 1, ndim
     xrel(ix^S, idim) = x(ix^S, idim) - vel(idim) * time
     val = val * sin(dpi * xrel(ix^S, idim))**4
  end do

end subroutine test_solution_1

!=============================================================================
! amrvacusr.t.testrho
!=============================================================================
