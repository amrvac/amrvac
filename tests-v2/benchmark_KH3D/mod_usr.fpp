module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision :: rhodens  = 10.0d0
  double precision :: rholight = 1.0d0
  double precision :: kx       = 4.0d0 * acos(-1.0d0)
  double precision :: kz       = 32d0*atan(1.0d0)
  double precision :: pint     = 2.5d0
  double precision :: dlin     = 0.025d0
  double precision :: sigma    = 0.00125d0

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
    double precision                :: vextra

    select case(iprob)
    case(1)
      vextra=0.0d0
    case(2)
      vextra = 10.0d0
    case(3)
      vextra=100.0d0
    case default
      error stop "Unknown iprob"
    endselect

    where(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)<0.25d0.or.x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>0.75d0)
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,rho_)=rholight
    elsewhere
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,rho_)=rhodens
    endwhere

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,p_)=pint

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)<=(0.25d0-dlin)).or.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)>=(0.75d0+dlin)))
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(1))=vextra-0.5d0
    endwhere

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>=(0.25d0+dlin)).and.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)<=(0.75d0-dlin)))
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(1))=vextra+0.5d0
    endwhere

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>(0.25d0-dlin)).and.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)<(0.25d0+dlin)))
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
         mom(1))=(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
         2)-(0.25d0-dlin))* ((vextra+0.5d0)-(vextra-0.5d0))/(2.0d0*dlin)+&
         (vextra-0.5d0)
    endwhere

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>(0.75d0-dlin)).and.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)<(0.75d0+dlin)))
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
         mom(1))=(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
         2)-(0.75d0-dlin))* ((vextra-0.5d0)-(vextra+0.5d0))/(2.0d0*dlin)+&
         (vextra+0.5d0)
    endwhere

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       mom(2))=0.01d0*dsin(kx*x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1))* (dexp(-0.5d0*(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2,ixGmin3:ixGmax3,2)-0.75d0)**2/sigma))
    
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       mom(3))=0.1d0*dsin(kz*x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       3))* (dexp(-0.5d0*(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)-0.75d0)**2/sigma))
   

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

end module mod_usr
