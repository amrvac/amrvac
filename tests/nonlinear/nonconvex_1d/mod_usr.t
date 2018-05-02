module mod_usr
  use mod_nonlinear

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer :: i_err
  integer :: i_sol

contains

  subroutine usr_init()
    use mod_variables

    usr_init_one_grid => initonegrid_usr
    usr_process_grid => set_error
    usr_print_log => print_error

    call nonlinear_activate()

    i_err = var_set_extravar("error", "error")
    i_sol = var_set_extravar("solution", "solution")

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    w(ix^S,rho_) = solution(x(ix^S, 1), 0.0_dp)

  end subroutine initonegrid_usr

  elemental function solution(x, t) result(val)
    real(dp), intent(in) :: x, t
    real(dp)             :: val, xrel
    real(dp), parameter  :: pi = acos(-1.0_dp)

    real(dp) :: uleft,uright,umid
    integer :: icase

    select case(iprob)
       case(1)
        ! compound case
        uleft=1.0d0
        uright=-1.0d0
       case(2)
        ! shock case
        uleft=1.0d0
        uright=0.1d0
       case(3)
        ! rarefaction case
        uleft=0.1d0
        uright=1.0d0
       case default
        uleft=1.0d0
        uright=1.0d0
    end select

    ! we start the discontinuity at t=0 always at x=0
    ! the following assumes uleft>0 positive always
    ! rarefaction case
    if(uleft<uright) icase=1
    ! shock case
    if(uright>-uleft/2.0d0 .and. uright<=uleft) icase=2 
    ! compound case
    if(uright<=-uleft/2.0d0) icase=3
        
    select case(icase)
      case(1) ! rarefaction
        if(x<=3.0d0*uleft**2*t ) val=uleft
        if(x>=3.0d0*uright**2*t) val=uright
        if(3.0d0*uleft**2*t<x .and. x<3.0d0*uright**2*t) val=dsqrt(x/(3.0d0*t))
      case(2) ! shock
        if(x<=(uleft**2+uleft*uright+uright**2)*t) val=uleft
        if(x>(uleft**2+uleft*uright+uright**2)*t) val=uright
      case(3) ! compound
        umid=-uleft/2.0d0
        if(x<=3.0d0*umid**2*t) val=uleft
        if(x>=3.0d0*uright**2*t) val=uright
        if(3.0d0*umid**2*t<x .and. x<3.0d0*uright**2*t) val=-dsqrt(x/(3.0d0*t))
    end select

  end function solution

  subroutine print_error()
    use mod_input_output, only: get_volume_average, get_global_maxima
    double precision   :: modes(nw, 2), volume
    double precision   :: maxvals(nw)

    call get_global_maxima(maxvals)
    call get_volume_average(1, modes(:, 1), volume)
    call get_volume_average(2, modes(:, 2), volume)

    if (mype == 0) then
       write(*, "(A,4E14.6)") " CONV_TEST (t,err_1,err_2,err_inf):", global_time, &
            modes(i_err, 1), sqrt(modes(i_err, 2)), maxvals(i_err)
    end if
  end subroutine print_error

  subroutine set_error(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,i_sol) = solution(x(ixO^S, 1), qt)
    w(ixO^S,i_err) = abs(w(ixO^S,rho_) - w(ixO^S,i_sol))
  end subroutine set_error

end module mod_usr

