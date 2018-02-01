! Smooth vortex test, e.g. in Duffell ApJS 226, 2016
module mod_usr
  use mod_hd

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer :: i_err_p, i_err_r
  integer :: i_sol

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_variables

    usr_init_one_grid => SV_init_one_grid
    usr_process_grid => set_error
    usr_print_log => print_error
    usr_refine_grid    => specialrefine_grid

    call set_coordinate_system("polar_2D")
    call hd_activate()

    i_err_r = var_set_extravar("rho_error", "rho_error")
    i_err_p = var_set_extravar("p_error", "p_error")
    i_sol = var_set_extravar("p_solution", "p_solution")

  end subroutine usr_init

  ! Initialize one grid
  subroutine SV_init_one_grid(ixG^L,ix^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    
    {^IFONED   call mpistop("This is a 2D HD problem") }
    {^IFTHREED call mpistop("This is a 2D HD problem") }

    w(ix^S,rho_) = 1.0d0
    {^IFTWOD
    w(ix^S,e_)   = p_solution(x(ix^S, 1))
    w(ix^S,mom(1)) = 0.0d0
    w(ix^S,mom(2)) = vphi_solution(x(ix^S, 1))
    where(x(ix^S,2)<dpi) 
      w(ix^S,tracer(1))=1.0
    elsewhere
      w(ix^S,tracer(1))=-1.0
    endwhere
    }

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine SV_init_one_grid

  elemental function p_solution(rad) result(val)
    use mod_global_parameters
    real(dp), intent(in) :: rad
    real(dp)             :: val
    real(dp) :: p0

    val=1.0d0-0.5d0*dexp(-rad**2)

  end function p_solution

  elemental function vphi_solution(rad) result(val)
    use mod_global_parameters
    real(dp), intent(in) :: rad
    real(dp)             :: val

    val=rad*dexp(-0.5d0*rad**2)

  end function vphi_solution

  subroutine set_error(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,i_sol) = p_solution(x(ixO^S, 1))
    call hd_to_primitive(ixI^L,ixO^L,w,x)
    w(ixO^S,i_err_p) = abs(w(ixO^S,e_) - w(ixO^S,i_sol))
    call hd_to_conserved(ixI^L,ixO^L,w,x)
    w(ixO^S,i_err_r) = abs(w(ixO^S,rho_) - 1.0d0)
  end subroutine set_error

  subroutine print_error()
    use mod_global_parameters
    use mod_input_output, only: get_volume_average, get_global_maxima

    double precision   :: modes(nw, 2), volume
    double precision   :: maxvals(nw)

    call get_global_maxima(maxvals)
    call get_volume_average(1, modes(:, 1), volume)
    call get_volume_average(2, modes(:, 2), volume)

    if (mype == 0) then
       write(*, "(A,7E14.6)") " CONVTEST (t rho_1 rho_2 rho_inf p_1 p_2 p_inf):", global_time, &
            modes(i_err_r, 1), sqrt(modes(i_err_r, 2)), maxvals(i_err_r), &
            modes(i_err_p, 1), sqrt(modes(i_err_p, 2)), maxvals(i_err_p)
    end if
  end subroutine print_error

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! test with different levels of refinement enforced
    !if (any((x(ix^S,1)<0.3d0).and.(x(ix^S,1)>0.1d0))) then
    !   refine=1
    !endif

  end subroutine specialrefine_grid

end module mod_usr
