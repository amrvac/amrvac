module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_gravity

    usr_init_one_grid => initonegrid_usr
    usr_source => specialsource
    usr_get_dt => getdt_special

    gravity_field = [0.0d0, -1.0d0, 0.0d0]

    call hd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    integer :: ix^D

    double precision:: y0,epsilon,rhodens,rholight,kx,kz,pint,dely
    logical::          first
    data first/.true./
    !----------------------------------------------------------------------------

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
       where(x(ixG^S,2)>y0+epsilon*sin(kx*x(ixG^S,1))*sin(kz*x(ixG^S,3)))
          w(ixG^S,rho_)=rhodens
       elsewhere
          w(ixG^S,rho_)=rholight
       endwhere
    else
       if(kx==zero)then
          where(x(ixG^S,2)>y0+epsilon*sin(kz*x(ixG^S,3)))
             w(ixG^S,rho_)=rhodens
          elsewhere
             w(ixG^S,rho_)=rholight
          endwhere
       else
          where(x(ixG^S,2)>y0+epsilon*sin(kx*x(ixG^S,1)))
             w(ixG^S,rho_)=rhodens
          elsewhere
             w(ixG^S,rho_)=rholight
          endwhere
       endif
    endif

    ! set all velocity to zero
    w(ixG^S, mom(:)) = zero

    ! pressure at interface
    pint=one
    if (hd_energy) then
       w(ixG^S,e_)=pint-w(ixG^S,rho_)*(x(ixG^S,2)-y0)
       w(ixG^S,e_)=w(ixG^S,e_)/(hd_gamma-one)
    end if
  end subroutine initonegrid_usr

  ! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  ! iw=iwmin...iwmax.  wCT is at time qCT
  subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    use mod_gravity

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    call gravity_add_source(qdt,ixI^L,ixO^L,iw^LIM, rho_, mom, e_, &
         qtC,wCT,qt,w,x)

  end subroutine specialsource

  ! Limit "dt" further if necessary, e.g. due to the special source terms.
  ! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
  ! module have already been called.
  subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_gravity

    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: dx^D, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

    call gravity_get_dt(w,ixG^L,ix^L,dtnew,dx^D,x)

  end subroutine getdt_special

end module mod_usr
