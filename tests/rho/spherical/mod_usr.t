module mod_usr
  use mod_rho

  implicit none

  integer :: i_sol,i_err

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_variables

    usr_init_one_grid => initonegrid_usr
    usr_process_grid => store_sol_err
    usr_print_log => print_min_max

    call set_coordinate_system("spherical_3D")
    call rho_activate()

    i_sol = var_set_extravar("solution", "solution")
    i_err = var_set_extravar("error", "error")
  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters
    use mod_physics

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rhoprofile(ixI^S)
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          select case(iprob)
             case(1)
                print *,'uniform advection in spherical coord'
             case(2)
                print *,'3D spheres advected in sphere: spherical coord'
             case(3)
                print *,'gaussian advected in spherical coordinates'
             case default
                call mpistop("iprob not available: edit mod_usr.t")
          end select
       end if
       first=.false.
    end if
    call set_density_profile(ixI^L,ixO^L,0.0d0,x,rhoprofile)
    w(ixO^S,rho_)=rhoprofile(ixO^S)

  end subroutine initonegrid_usr

  subroutine set_density_profile(ixI^L,ixO^L,qt,x,rhoprofile)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qt,x(ixI^S,1:ndim)
    double precision, intent(out) :: rhoprofile(ixI^S)

    double precision :: rad,xcc^D
    double precision :: xcart(ixI^S,1:ndim),rr2(ixI^S)
    logical :: mask(ixI^S)

    rhoprofile(ixO^S)=1.d0
    if(iprob==2) then
       ! radius of the spheres
       rad=0.2d0*(xprobmax1-xprobmin1)*0.5d0
       ! cartesian coordinates from r, theta, phi
       xcart(ixO^S,1)=x(ixO^S,1)*dsin(x(ixO^S,2))*dcos(x(ixO^S,3))
       xcart(ixO^S,2)=x(ixO^S,1)*dsin(x(ixO^S,2))*dsin(x(ixO^S,3))
       xcart(ixO^S,3)=x(ixO^S,1)*dcos(x(ixO^S,2))
       ! set up sphere positioned on positive x-axis
       xcc1=xprobmax1*0.5d0+rho_v(1)*qt
       xcc2=rho_v(2)*qt
       xcc3=rho_v(3)*qt
       where({^D&(xcart(ixO^S,^D)-xcc^D)**2+}<rad**2)
         rhoprofile(ixO^S)=100.d0
       endwhere
       ! set up sphere positioned on positive y-axis
       xcc1=rho_v(1)*qt
       xcc2=xprobmax1*0.5d0+rho_v(2)*qt
       xcc3=rho_v(3)*qt
       where({^D&(xcart(ixO^S,^D)-xcc^D)**2+}<rad**2)
         rhoprofile(ixO^S)=100.d0
       endwhere
       ! set up sphere positioned on positive z-axis
       xcc1=rho_v(1)*qt
       xcc2=rho_v(2)*qt
       xcc3=xprobmax1*0.5d0+rho_v(3)*qt
       where({^D&(xcart(ixO^S,^D)-xcc^D)**2+}<rad**2)
         rhoprofile(ixO^S)=100.d0
       endwhere
       ! set up sphere positioned on negative z-axis
       xcc1=rho_v(1)*qt
       xcc2=rho_v(2)*qt
       xcc3=-xprobmax1*0.5d0+rho_v(3)*qt
       where({^D&(xcart(ixO^S,^D)-xcc^D)**2+}<rad**2)
         rhoprofile(ixO^S)=100.d0
       endwhere
    end if
    if(iprob==3)then
       rad=0.2d0*(xprobmax1-xprobmin1)*0.5d0
       ! cartesian coordinates from R, Z, varphi
       xcart(ixO^S,1)=x(ixO^S,1)*dsin(x(ixO^S,2))*dcos(x(ixO^S,3))
       xcart(ixO^S,2)=x(ixO^S,1)*dsin(x(ixO^S,2))*dsin(x(ixO^S,3))
       xcart(ixO^S,3)=x(ixO^S,1)*dcos(x(ixO^S,2))
       ! set up gaussian bell in sphere positioned on positive x-axis
       xcc1=(xprobmax1-xprobmin1)*0.5d0+rho_v(1)*qt
       xcc2=rho_v(2)*qt
       xcc3=rho_v(3)*qt
       mask(ixO^S)=({^D&(xcart(ixO^S,^D)-xcc^D)**2+}<=rad**2)
       rr2(ixO^S)={^D&(xcart(ixO^S,^D)-xcc^D)**2+}
       where(mask(ixO^S))
         rhoprofile(ixO^S) = dexp(-rr2(ixO^S)/rad**2)
       elsewhere
         rhoprofile(ixO^S) = dexp(-1.0d0)
       endwhere
    end if

  end subroutine set_density_profile

  subroutine store_sol_err(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rhoprofile(ixI^S)

    call set_density_profile(ixI^L,ixO^L,qt,x,rhoprofile)
    w(ixO^S,i_sol) = rhoprofile(ixO^S)
    w(ixO^S,i_err) = dabs(w(ixO^S,rho_) - w(ixO^S,i_sol))
  end subroutine store_sol_err

  subroutine print_min_max()
    use mod_global_parameters
    use mod_input_output, only: get_global_minima, get_global_maxima, &
                                get_volume_average
    double precision   :: minvals(nw),maxvals(nw)

    integer :: iw
    double precision :: modes(nw,2), volume

    call get_global_minima(minvals)
    call get_global_maxima(maxvals)
    call get_volume_average(1,modes(:,1),volume)
    call get_volume_average(2,modes(:,2),volume)

    if (mype == 0) then
       write(*, "(A,4E16.8,A,3E16.8)") " time + rho min-max-tot:", &
              global_time, minvals(rho_), maxvals(rho_), modes(rho_,1), &
             " Error: Linf-L1-L2:", &
              maxvals(i_err),modes(i_err,1),dsqrt(modes(i_err,2))
    end if
  end subroutine print_min_max

end module mod_usr

