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

! Set velocity to one
eqpar(v1_:v1_+ndim-1) = 1.0d0

end subroutine initglobaldata_usr

! initialize one grid
subroutine initonegrid_usr(ixG^L,ix^L,w,x)
  use mod_global_parameters

  integer, intent(in)             :: ixG^L, ix^L
  double precision, intent(in)    :: x(ixG^S,1:ndim)
  double precision, intent(inout) :: w(ixG^S,1:nw)
  double precision                :: x_vec(ndim)
  integer                         :: i^D

   {do i^D = ixmin^D, ixmax^D\}
   x_vec = x(i^D, :)
   call test_solution(x_vec, 0.0d0, w(i^D, rho_))
   {end do\}

end subroutine initonegrid_usr

subroutine test_solution(x_vec, time, val)
  use mod_global_parameters

  double precision, intent(in)  :: x_vec(ndim), time
  double precision, intent(out) :: val
  double precision              :: xrel(ndim)
  double precision              :: vel(ndim)

  vel  = eqpar(v1_:v1_+ndim-1)
  xrel = x_vec - vel * time

  select case (iprob)
  case (1)
     val = product(sin(dpi * xrel))**4
  case default
     call mpistop("iprob not available!")
  end select
end subroutine test_solution

!=============================================================================
! amrvacusr.t.testrho
!=============================================================================
