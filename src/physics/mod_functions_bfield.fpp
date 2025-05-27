module mod_functions_bfield

  implicit none
  private


  public :: get_divb

  !> Indices of the magnetic field
  integer, allocatable, public :: mag(:)

contains


  !> Calculate div B within ixO
  subroutine get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb, fourthorder)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    logical, intent(in), optional   :: fourthorder

    integer                            :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
       ixCmax2,ixCmax3, idir

    if(stagger_grid) then
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=0.d0
      do idir=1,ndim
        ixCmin1=ixOmin1-kr(idir,1);ixCmin2=ixOmin2-kr(idir,2)
        ixCmin3=ixOmin3-kr(idir,3);ixCmax1=ixOmax1-kr(idir,1)
        ixCmax2=ixOmax2-kr(idir,2);ixCmax3=ixOmax3-kr(idir,3);
        divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+block%ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idir)*block%surfaceC(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)-block%ws(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)*block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,idir)
      end do
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else
      select case(typediv)
      case("central")
        call divvector(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           mag(1:ndir)),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb,fourthorder)
      case("limited")
        call divvectorS(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           mag(1:ndir)),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
      end select
    end if

  end subroutine get_divb


end module mod_functions_bfield
