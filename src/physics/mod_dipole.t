!> sets up a magnetic dipole in a 3D cartesian box
!  parameters for dipole moment: 
!      strength mu_dipole (factor 4pi included)
!      polar angle in degrees theta_dipole [0 to 90]
!      azimuthal angle in degrees phi_dipole [0 to 360]
!  center of dipole at x1_dip, x2_dip, x3_dip
module mod_dipole
  implicit none
  double precision :: mu_dipole,theta_dipole,phi_dipole,x1_dip,x2_dip,x3_dip
  double precision :: mx, my, mz

contains

{^IFTHREED
  subroutine dipole(ixI^L,ixO^L,x,B)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: B(ixI^S,1:ndir)

    double precision :: xsincos(ixI^S),xsinsin(ixI^S),xcos(ixI^S)
    double precision :: rr0(ixI^S),mdotr(ixI^S)
 
    B=0.0d0
    mx=mu_dipole*dsin(theta_dipole*dpi/180.0d0)*dcos(phi_dipole*dpi/180.d0)
    my=mu_dipole*dsin(theta_dipole*dpi/180.0d0)*dsin(phi_dipole*dpi/180.0d0)
    mz=mu_dipole*dcos(theta_dipole*dpi/180.0d0)
    rr0(ixO^S)=dsqrt((x(ixO^S,1)-x1_dip)**2 &
                    +(x(ixO^S,2)-x2_dip)**2 &
                    +(x(ixO^S,3)-x3_dip)**2)
    if(any(rr0(ixO^S)<smalldouble)) call mpistop("dipole center at cell center")
    !where(rr0(ixO^S)<smalldouble) 
    !   rr0(ixO^S)=bigdouble
    !endwhere
    xsincos(ixO^S)=(x(ixO^S,1)-x1_dip)/rr0(ixO^S)
    xsinsin(ixO^S)=(x(ixO^S,2)-x2_dip)/rr0(ixO^S)
    xcos(ixO^S)   =(x(ixO^S,3)-x3_dip)/rr0(ixO^S)
    mdotr(ixO^S)=mx*xsincos(ixO^S)+my*xsinsin(ixO^S)+mz*xcos(ixO^S)
    B(ixO^S,1)=(3.0d0*xsincos(ixO^S)*mdotr(ixO^S)-mx)/rr0(ixO^S)**3
    B(ixO^S,2)=(3.0d0*xsinsin(ixO^S)*mdotr(ixO^S)-my)/rr0(ixO^S)**3
    B(ixO^S,3)=(3.0d0*xcos(ixO^S)   *mdotr(ixO^S)-mz)/rr0(ixO^S)**3

  end subroutine dipole
}
end module mod_dipole
