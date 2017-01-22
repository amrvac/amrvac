! setup.pl -d=33
! test radiative cooling in a central hot ball
module mod_usr
  use mod_hd
  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    unit_length        = 1.d9   ! cm
    unit_temperature   = 1.d6   ! K
    unit_numberdensity = 1.d9   ! cm^-3

    usr_init_one_grid => initonegrid_usr

    call hd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_global_parameters
    use mod_physics

    integer, intent(in)             :: ixG^L,ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    ! set all velocity to zero
    w(ixG^S, mom(:)) = zero
    ! uniform pressure
    w(ixG^S,p_)   =2.d0
    ! hot central circular spot with uniform pressure
    where((^D&x(ixG^S,^D)**2+) .lt. 0.25d0**2)
      w(ixG^S,rho_) =1.d0
    elsewhere
      w(ixG^S,rho_) =10.d0
    endwhere

    call phys_to_conserved(ixG^L,ixG^L,w,x)

  end subroutine initonegrid_usr

end module mod_usr
