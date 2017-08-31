module mod_usr
  use mod_hd

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
    use mod_usr_methods

    usr_init_one_grid => initonegrid_usr
    call set_coordinate_system("Cartesian_3D")

    call hd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
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

    where(x(ixG^S,2)<0.25d0.or.x(ixG^S,2)>0.75d0)
      w(ixG^S,rho_)=rholight
    elsewhere
      w(ixG^S,rho_)=rhodens
    endwhere

    w(ixG^S,p_)=pint

    where((x(ixG^S,2)<=(0.25d0-dlin)).or.(x(ixG^S,2)>=(0.75d0+dlin)))
      w(ixG^S,mom(1))=vextra-0.5d0
    endwhere

    where((x(ixG^S,2)>=(0.25d0+dlin)).and.(x(ixG^S,2)<=(0.75d0-dlin)))
      w(ixG^S,mom(1))=vextra+0.5d0
    endwhere

    where((x(ixG^S,2)>(0.25d0-dlin)).and.(x(ixG^S,2)<(0.25d0+dlin)))
      w(ixG^S,mom(1))=(x(ixG^S,2)-(0.25d0-dlin))* &
           ((vextra+0.5d0)-(vextra-0.5d0))/(2.0d0*dlin)+(vextra-0.5d0)
    endwhere

    where((x(ixG^S,2)>(0.75d0-dlin)).and.(x(ixG^S,2)<(0.75d0+dlin)))
      w(ixG^S,mom(1))=(x(ixG^S,2)-(0.75d0-dlin))* &
           ((vextra-0.5d0)-(vextra+0.5d0))/(2.0d0*dlin)+(vextra+0.5d0)
    endwhere

    w(ixG^S,mom(2))=0.01d0*dsin(kx*x(ixG^S,1))* &
         (dexp(-0.5d0*(x(ixG^S,2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixG^S,2)-0.75d0)**2/sigma))
    {^IFTHREED
    w(ixG^S,mom(3))=0.1d0*dsin(kz*x(ixG^S,3))* &
         (dexp(-0.5d0*(x(ixG^S,2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixG^S,2)-0.75d0)**2/sigma))
    }

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

end module mod_usr
