! Gombosi et al. 2002 JCP 177, 176, section 7.2 Alfven wave
module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr

    call set_coordinate_system("Cartesian_1.75D")
    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          print *,'semirelativisitic Alfven wave in Cartesian coordinate'
       end if
       first=.false.
    end if
    w(ixO^S,rho_)=1.d0
    w(ixO^S,mom(1))=0.5d0
    w(ixO^S,mom(2))=0.d0
    where(x(ixO^S,1)<0.d0)
      w(ixO^S,mom(3))=0.01d0
    else where
      w(ixO^S,mom(3))=-0.01d0
    end where
    w(ixO^S,p_)=0.1d0
    w(ixO^S,mag(1))=1.d0
    w(ixO^S,mag(2))=1.d0
    w(ixO^S,mag(3))=0.d0

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

end module mod_usr
