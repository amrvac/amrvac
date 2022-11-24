module mod_usr
  use mod_hd
  implicit none

contains

  subroutine usr_init()
    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr

    call hd_activate()
    call set_coordinate_system('Cartesian_1D')
  end subroutine usr_init

  subroutine initglobaldata_usr
    hd_gamma = 1.4d0
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    where(x(ix^S,1) .lt. -4)
      w(ix^S,rho_)=3.857143d0
      w(ix^S,mom(1))=2.629369d0
      w(ix^S,p_)=10.33333d0
    else where
      w(ix^S,rho_)=1.d0+0.2d0*dsin(5.d0*x(ix^S,1))
      w(ix^S,mom(1))=zero
      w(ix^S,p_)=1.d0
    end where
    call hd_to_conserved(ixG^L,ix^L,w,x)
  end subroutine initonegrid_usr
end module mod_usr
