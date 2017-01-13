! setup.pl -d=33
! test thermal conduction in a central hot ball
module mod_usr
  use mod_hd
  implicit none

contains

  subroutine usr_init()
    use mod_usr_methods

    call initglobaldata_usr

    usr_init_one_grid => initonegrid_usr

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters

    double precision :: munit,kB

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3
    munit              = 1.67262d-24                                  ! g
    unit_density       = 1.4d0*munit*unit_numberdensity               ! 2.34167d-15 g*cm^-3
    kB                 = 1.3806d-16                                   ! erg*K^-1
    unit_velocity      = dsqrt(unit_temperature*kB/munit*2.3d0/1.4d0) ! 1.16449d7 cm/s
    unit_time          = unit_length/unit_velocity                    ! 85.87461 s

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_global_parameters
    use mod_physics

    integer, intent(in)             :: ixG^L,ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision                :: radius, xcircle^D

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
