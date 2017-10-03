! setup.pl -d=3
! test thermal conduction in circular magnetic field loops
module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    usr_init_one_grid => initonegrid_usr

    call set_coordinate_system("Cartesian")
    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: r(ixI^S), theta(ixI^S), B(ixI^S), al,ab
    double precision :: xp(ixI^S,1:ndim),xpp(ixI^S,1:ndim)

    al=0.25d0*dpi
    ab=0.25d0*dpi
    ! first rotate around z axis
    xp(ixO^S,1)=dcos(al)*x(ixO^S,1)-dsin(al)*x(ixO^S,2)
    xp(ixO^S,2)=dsin(al)*x(ixO^S,1)+dcos(al)*x(ixO^S,2)
    xp(ixO^S,3)=x(ixO^S,3)
    ! then rotate around x axis
    xpp(ixO^S,1)=xp(ixO^S,1)
    xpp(ixO^S,2)=dcos(ab)*xp(ixO^S,2)-dsin(ab)*xp(ixO^S,3)
    xpp(ixO^S,3)=dsin(ab)*xp(ixO^S,2)+dcos(ab)*xp(ixO^S,3)
    ! set all velocity to zero
    w(ixO^S, mom(:)) = zero
    ! uniform density
    w(ixO^S,rho_) =1.d0
    r(ixO^S)=dsqrt(xpp(ixO^S,1)**2+xpp(ixO^S,2)**2)
    where(xpp(ixO^S,1)>0.d0)
      theta(ixO^S)=atan(xpp(ixO^S,2)/xpp(ixO^S,1))
    elsewhere(xpp(ixO^S,1)<0.d0)
      theta(ixO^S)=dpi-atan(xpp(ixO^S,2)/abs(xpp(ixO^S,1)))
    elsewhere
      theta(ixO^S)=0.d0
    endwhere
    ! hot central circular spot with uniform pressure
    where(xpp(ixO^S,3)>-0.2d0 .and. xpp(ixO^S,3)<0.2d0 .and. r(ixO^S)>0.5d0 .and. r(ixO^S)<0.7d0 .and. &
        theta(ixO^S)>11.d0/12.d0*dpi .and. theta(ixO^S)<13.d0/12.d0*dpi)
      w(ixO^S,p_)=12.d0
    elsewhere
      w(ixO^S,p_)=10.d0
    endwhere
    ! straight line current
    B(ixO^S)=Busr/r(ixO^S)
    w(ixO^S,mag(1))=B(ixO^S)*dcos(theta(ixO^S)+0.5*dpi)
    w(ixO^S,mag(2))=B(ixO^S)*dsin(theta(ixO^S)+0.5*dpi)
    w(ixO^S,mag(3))=0.d0
    ! first rotate around a axis -ab rad
    xp(ixO^S,1)=w(ixO^S,mag(1))
    xp(ixO^S,2)=dcos(-ab)*w(ixO^S,mag(2))-dsin(-ab)*w(ixO^S,mag(3))
    xp(ixO^S,3)=dsin(-ab)*w(ixO^S,mag(2))+dcos(-ab)*w(ixO^S,mag(3))
    ! then rotate around z axis -al rad
    w(ixO^S,mag(1))=dcos(-al)*xp(ixO^S,1)-dsin(-al)*xp(ixO^S,2)
    w(ixO^S,mag(2))=dsin(-al)*xp(ixO^S,1)+dcos(-al)*xp(ixO^S,2)
    w(ixO^S,mag(3))=xp(ixO^S,3)

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

end module mod_usr
