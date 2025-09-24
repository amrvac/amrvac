! setup.pl -d=2
! test thermal conduction in a ring in ffhd module
module mod_usr
  use mod_ffhd
  implicit none

contains

  subroutine usr_init()
    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    usr_init_one_grid   => initonegrid_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output 
    usr_set_B0          => specialset_B0

    call set_coordinate_system("Cartesian")
    call ffhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: r(ixI^S), theta(ixI^S), B(ixI^S)

    ! set all velocity to zero
    w(ixO^S, mom(:)) = zero
    ! uniform pressure
    w(ixO^S,rho_) =1.d0
    r(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
    where(x(ixO^S,1)>0.d0)
      theta(ixO^S)=atan(x(ixO^S,2)/x(ixO^S,1))
    elsewhere(x(ixO^S,1)<0.d0)
      theta(ixO^S)=dpi-atan(x(ixO^S,2)/abs(x(ixO^S,1)))
    elsewhere
      theta(ixO^S)=0.d0
    endwhere
    ! hot central circular spot with uniform pressure
    where(r(ixO^S)>0.5d0 .and. r(ixO^S)<0.7d0 .and. theta(ixO^S)>11.d0/12.d0*dpi .and. theta(ixO^S)<13.d0/12.d0*dpi)
      w(ixO^S,p_)=4.d0
    elsewhere
      w(ixO^S,p_)=1.d0
    endwhere

    call ffhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    double precision :: r(ixI^S), theta(ixI^S), B(ixI^S)

    r(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
    where(x(ixO^S,1)>0.d0)
      theta(ixO^S)=atan(x(ixO^S,2)/x(ixO^S,1))
    elsewhere(x(ixO^S,1)<0.d0)
      theta(ixO^S)=dpi-atan(x(ixO^S,2)/abs(x(ixO^S,1)))
    elsewhere
      theta(ixO^S)=0.d0
    endwhere
    ! straight line current
    B(ixO^S)=Busr/r(ixO^S)
    wB0(ixO^S,1)=B(ixO^S)*dcos(theta(ixO^S)+0.5*dpi)
    wB0(ixO^S,2)=B(ixO^S)*dsin(theta(ixO^S)+0.5*dpi)
    B(ixO^S)=dsqrt(wB0(ixO^S,1)**2+wB0(ixO^S,2)**2)+1.d-319
    wB0(ixO^S,1)=wB0(ixO^S,1)/B(ixO^S)
    wB0(ixO^S,2)=wB0(ixO^S,2)/B(ixO^S)
  end subroutine specialset_B0

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S)

    call ffhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te'
  end subroutine specialvarnames_output

end module mod_usr
