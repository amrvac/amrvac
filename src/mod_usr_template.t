!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_mhd

  implicit none

  ! Custom global variables can be defined here
  ! ...

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2.5D")

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see src/mod_usr_methods.t
    ! ...

    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call mhd_activate()

  end subroutine usr_init

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! index for do loop and cell faces
    integer :: idir, ixC^L

    ! Set initial values for w
    ! 1.d0 and 0.d0 are just examples should be adjusted to user's problem

    ! cell-center density
    w(ixO^S, rho_) = 1.d0
    ! cell-center velocity 
    w(ixO^S, mom(1)) = 0.d0
    w(ixO^S, mom(2)) = 0.d0
    ! cell-center pressure
    w(ixO^S, e_) = 1.d0
    if(stagger_grid) then
      ! set cell-face magnetic field B using CT method for divB control

      ! set B from vector potential (zero divB guaranteed) given 
      ! usr_init_vector_potential pointing to  a subroutine to set vector potential
      call b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block%ws,x)

      ! or directly set cell-face B (divB maybe non-zero) as following:
      do idir=1,ndim
        ixCmin^D=ixImin^D;
        ixCmax^D=ixImax^D-kr(idir,^D);
        ! cell-face B_idir
        block%ws(ixC^S,idir)=1.d0
      end do
      ! update cell-center B from cell-face B
      call mhd_face_to_center(ixO^L,block)
    else
      ! cell-center magnetic field
      w(ixO^S,mag(1)) = 1.d0
      w(ixO^S,mag(2)) = 1.d0
    end if

    ! convert primitive variables to conservative variables
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

  ! Extra routines can be placed here
  ! ...

end module mod_usr
