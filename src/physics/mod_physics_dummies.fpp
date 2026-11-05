#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute
! Dummy routines which can be overwritten by a physics-dependent implementation
! On non-Cray compilers, these call mpistop with a descriptive message.
! On Cray, STOP cannot be inlined into OpenACC kernels, so dummies just return -1. 

#:def estimate_speeds_minmax()
subroutine estimate_speeds_minmax(uL, uR, xC, flux_dim, wL, wR)
#ifndef _CRAYFTN
  use mod_comm_lib, only: mpistop
#endif
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: wL, wR

#ifndef _CRAYFTN
  call mpistop("estimate_speeds_minmax not implemented for this physics module")
#endif
  wL = -1._dp
  wR = -1._dp

end subroutine estimate_speeds_minmax
#:enddef


#:def estimate_speeds_hllc()
subroutine estimate_speeds_hllc(uL, uR, xC, flux_dim, sL, sR)
#ifndef _CRAYFTN
  use mod_comm_lib, only: mpistop
#endif
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer,  intent(in)  :: flux_dim
  real(dp), intent(out) :: sL, sR

#ifndef _CRAYFTN
  call mpistop("estimate_speeds_hllc not implemented for this physics module")
#endif
  sL = -1._dp
  sR = -1._dp

end subroutine estimate_speeds_hllc
#:enddef


#:def hllc_flux()
subroutine hllc_flux(uL, uR, flux_l, flux_r, sL, sR, flux_dim, F, degenerate)
#ifndef _CRAYFTN
  use mod_comm_lib, only: mpistop
#endif
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(inout) :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)    :: flux_l(nw_flux), flux_r(nw_flux)
  real(dp), intent(in)    :: sL, sR
  integer,  intent(in)    :: flux_dim
  real(dp), intent(out)   :: F(nw_flux)
  logical,  intent(out)   :: degenerate

#ifndef _CRAYFTN
  call mpistop("hllc_flux not implemented for this physics module")
#endif
  F = 0._dp
  degenerate = .true.

end subroutine hllc_flux
#:enddef


#:def addsource_nonlocal()
subroutine addsource_nonlocal(qdt, dtfactor, qtC, wCTprim, qt, wnew, x, dx, idir, &
     qsourcesplit)
  ${GPU_ROUTINE_SEQ()}$

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCTprim(nw_phys,5)
  real(dp), intent(in)     :: x(1:ndim), dx(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  integer, intent(in)      :: idir
  logical, intent(in)      :: qsourcesplit


end subroutine addsource_nonlocal
#:enddef

#:def addsource_compact()
subroutine addsource_compact(qdt, dtfactor, qtC, wCTprim1, wCTprim2, wCTprim3, qt, wnew, x, dx, &
     qsourcesplit)
  ${GPU_ROUTINE_SEQ()}$

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCTprim1(nw_phys,3),wCTprim2(nw_phys,3),wCTprim3(nw_phys,3)
  real(dp), intent(in)     :: x(1:ndim), dx(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit


end subroutine addsource_compact
#:enddef

#:def addsource_geometry()
subroutine addsource_geometry(qdt, wprim, wnew, x, dAdV)
#ifndef _CRAYFTN
  use mod_comm_lib, only: mpistop
#endif
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)     :: qdt
  real(dp), intent(in)     :: wprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim)
  !> (upper minus lower face area) / cell volume, per direction
  real(dp), intent(in)     :: dAdV(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)

#ifndef _CRAYFTN
  call mpistop("curvilinear geometry not implemented for this physics module")
#endif

end subroutine addsource_geometry
#:enddef

#:def fix_prim_state()
!> Make a reconstructed primitive face state admissible before it is handed to
!> the Riemann solver. Genuinely empty by default: a physics module without
!> auxiliary variables and without a floor of its own needs nothing here, and
!> the call is inlined away.
pure subroutine fix_prim_state(u)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(inout) :: u(nw_phys)

end subroutine fix_prim_state
#:enddef

! Physics modules with auxiliary variables (nwaux > 0) that are derived from
! the conserved state override these. get_aux has a physics-specific signature
! and is only called by that physics' own to_primitive and set_aux, so
! everyone else needs no default for it.
#:def get_aux()
#:enddef

#:def set_aux()
!> Overwrite the auxiliary slots of a conserved state with the values
!> consistent with its other components. Genuinely empty by default: it is
!> called by set_aux_all, prolong_2nd and coarsen_grid for every physics, and
!> a physics module without auxiliaries has nothing to refresh.
pure subroutine set_aux(u)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(inout) :: u(nw_phys)

end subroutine set_aux
#:enddef
