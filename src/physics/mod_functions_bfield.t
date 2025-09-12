module mod_functions_bfield

  implicit none
  private

  public :: get_divb

  !> Indices of the magnetic field
  integer, allocatable, public :: mag(:)

contains

  !> Calculate div B within ixO
  subroutine get_divb(w,ixI^L,ixO^L,divb,nth_in)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: divb(ixI^S)
    integer, intent(in), optional   :: nth_in
    integer                         :: ixC^L, idir, nth

    if(present(nth_in)) then
      nth=nth_in
    else
      nth=1
    endif

    if(stagger_grid) then
      divb(ixO^S)=0.d0
      do idir=1,ndim
        ixC^L=ixO^L-kr(idir,^D);
        divb(ixO^S)=divb(ixO^S)+block%ws(ixO^S,idir)*block%surfaceC(ixO^S,idir)-&
                                block%ws(ixC^S,idir)*block%surfaceC(ixC^S,idir)
      end do
      divb(ixO^S)=divb(ixO^S)/block%dvolume(ixO^S)
    else
      select case(typediv)
      case("central")
        call divvector(w(ixI^S,mag(1:ndir)),ixI^L,ixO^L,divb,nth)
      case("limited")
        call divvectorS(w(ixI^S,mag(1:ndir)),ixI^L,ixO^L,divb)
      end select
    end if
  end subroutine get_divb


end module mod_functions_bfield
