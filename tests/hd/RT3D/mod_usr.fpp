module mod_usr
  use mod_physics
  use mod_amrvac

  implicit none
  
contains

  subroutine usr_init()

    call set_coordinate_system("Cartesian_3D")

    usr_init_one_grid => initonegrid_usr

    call phys_activate()

  end subroutine usr_init

  pure real(dp) function gravity_field(wCT, x, idim) result(field)
    !$acc routine seq
    real(dp), intent(in)    :: wCT(nw_phys)
    real(dp), intent(in)    :: x(1:ndim)
    integer, value, intent(in)     :: idim
    ! real(dp)                :: field

    if (idim == 1) field =  0.0_dp
    if (idim == 2) field =  0.0_dp
    if (idim == 3) field = -1.0_dp

  end function gravity_field

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    integer :: ix1,ix2,ix3

    double precision:: z0,epsilon,rhodens,rholight,kx,ky,pint
    logical::          first
    data first/.true./


    ! the location of demarcation line`
    z0=0.8d0

    ! density of two types
    rhodens=one
    rholight=0.1d0

    ! setup the perturbation
    epsilon=0.05d0
    ! kx=2 pi
    kx=8.0d0*atan(one)
    ! ky=8 pi
    ky=4.0d0*8.0d0*atan(one)

    ! print out the info
    if (first) then
       if (mype==0) then
          print *,'HD Rayleigh Taylor problem'
          print *,'  --assuming z ranging from 0-1!'
          print *,'  --interface z0-epsilon:',z0,epsilon
          print *,'  --density ratio:',rhodens/rholight
          print *,'  --kx:',kx
          print *,'  --ky:',ky
       end if
       first=.false.
    end if

    ! initialize the density
    where(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)>z0+epsilon*sin(kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1))*sin(ky*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=rhodens
    elsewhere
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=rholight
    endwhere

    ! set all velocity to zero
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(:)) = zero

    ! pressure at interface
    pint=one
    if(hd_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_)=pint-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)-z0)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)/(hd_gamma-one)
    end if
  end subroutine initonegrid_usr

end module mod_usr
