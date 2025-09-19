module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

contains

  subroutine usr_init()

    call set_coordinate_system("Cartesian_3D")

    usr_init_one_grid => initonegrid_usr

    call phys_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    double precision :: rbs,xc1,xc2,xc3
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          print *,'MHD blast wave in Cartesian coordinate'
       end if
       first=.false.
    end if

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,rho_)=1.d0
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,p_)=1.d0
    rbs=0.2d0
    xc1=0.5d0*(xprobmin1+xprobmax1)
    xc2=0.5d0*(xprobmin2+xprobmax2)
    xc3=0.5d0*(xprobmin3+xprobmax3)

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1)-xc1)**2+&
     (x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,2)-xc2)**2+&
     (x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,3)-xc3)**2 < rbs**2)
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,p_)=100.d0
    endwhere

   
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mag(1))=1.d0
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mag(2))=1.d0
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mag(3))=1.d0

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  !> analytical fomula for the unit vectors along B
  pure real(dp) function bfield(x, idim) result(field)
    !$acc routine seq
    real(dp), intent(in)    :: x(1:ndim)
    integer, value, intent(in)     :: idim

    if (idim == 1) field =  1.0_dp
    if (idim == 2) field =  1.0_dp
    if (idim == 3) field =  1.0_dp

  end function bfield

end module mod_usr
