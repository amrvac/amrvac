!> test isotropic thermal conduction from a central hot ball
module mod_usr
  use mod_hd
  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    usr_init_one_grid => initonegrid_usr

    call set_coordinate_system("Cartesian_2D")

    call hd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! set all velocity to zero
    w(ixO^S, mom(:)) = zero
    ! uniform pressure
    w(ixO^S,p_)   =2.d0
    ! hot central circular spot with uniform pressure
    where((^D&x(ixO^S,^D)**2+) .lt. 0.25d0**2)
     w(ixO^S,rho_) =1.d0
    elsewhere
     w(ixO^S,rho_) =10.d0
    endwhere

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

end module mod_usr
