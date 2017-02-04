module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()
    use mod_usr_methods

    usr_init_one_grid => initonegrid_usr
    usr_gravity=> gravity

    call hd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ix^D

    double precision:: y0,epsilon,rhodens,rholight,kx,kz,pint,dely
    logical::          first
    data first/.true./


    ! the location of demarcation line`
    y0=0.8d0

    ! density of two types
    rhodens=one
    rholight=0.1d0

    ! setup the perturbation
    epsilon=0.05d0
    ! kx=2 pi
    kx=8.0d0*atan(one)
    ! kz=8 pi
    kz=32d0*atan(one)

    ! print out the info
    if (first) then
       if (mype==0) then
          print *,'HD Rayleigh Taylor problem'
          print *,'  --assuming y ranging from 0-1!'
          print *,'  --interface y0-epsilon:',y0,epsilon
          print *,'  --density ratio:',rhodens/rholight
          print *,'  --kx:',kx
          print *,'  --kz:',kz
       end if
       first=.false.
    end if

    ! initialize the density
    if(kx*kz/=zero)then
      where(x(ixO^S,2)>y0+epsilon*sin(kx*x(ixO^S,1))*sin(kz*x(ixO^S,3)))
        w(ixO^S,rho_)=rhodens
      elsewhere
        w(ixO^S,rho_)=rholight
      endwhere
    else
      if(kx==0.d0) then
        where(x(ixO^S,2)>y0+epsilon*sin(kx*x(ixO^S,3)))
          w(ixO^S,rho_)=rhodens
        elsewhere
          w(ixO^S,rho_)=rholight
        endwhere
      else
        where(x(ixO^S,2)>y0+epsilon*sin(kx*x(ixO^S,1)))
          w(ixO^S,rho_)=rhodens
        elsewhere
          w(ixO^S,rho_)=rholight
        endwhere
      end if
    endif

    ! set all velocity to zero
    w(ixO^S, mom(:)) = zero

    ! pressure at interface
    pint=one
    if(hd_energy) then
      w(ixO^S,e_)=pint-w(ixO^S,rho_)*(x(ixO^S,2)-y0)
      w(ixO^S,e_)=w(ixO^S,e_)/(hd_gamma-one)
    end if
  end subroutine initonegrid_usr

  ! Calculate gravitational acceleration in each dimension
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:)=0.d0
    gravity_field(ixO^S,2)=-1.d0

  end subroutine gravity

end module mod_usr
