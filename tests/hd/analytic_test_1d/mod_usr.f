module mod_usr
  use mod_hd
  use mod_riemann_exact

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer            :: i_err_rho
  integer            :: i_sol_rho
  type(riemann_t)    :: rp

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_initialize
    use mod_variables

    usr_init_one_grid => rm1d_init_one_grid
    usr_process_global => set_error
    usr_print_log => print_error

    call set_coordinate_system("Cartesian")
    call hd_activate()

    i_err_rho = var_set_extravar("err_rho", "err_rho")
    i_sol_rho = var_set_extravar("sol_rho", "sol_rho")

    call initialize_amrvac()
    call select_riemann_problem()

  end subroutine usr_init

  subroutine select_riemann_problem()
    use mod_global_parameters

    select case (iprob)
    case (1)
       ! state at left of discontinuity
       rp%rhol = 10.0
       rp%pl   = 100.0
       rp%ul   = 0.

       ! state at right of discontinuity
       rp%rhor = 1.0
       rp%pr   = 1.0
       rp%ur   = 0.

       !  equation of state
       rp%gamma = hd_gamma

       ! spatial interval over which to compute solution
       rp%xl = xprobmin1 + 0.5_dp * dxlevel(1)
       rp%xr = xprobmax1 - 0.5_dp * dxlevel(1)

       ! location of discontinuity at t = 0
       rp%xi = 0.5_dp * (rp%xl + rp%xr)
    case (2)
       ! Sod shock tube
       rp%rhol = 1.0
       rp%pl   = 1.0
       rp%ul   = 0.
       rp%rhor = 0.125
       rp%pr   = 0.1
       rp%ur   = 0.

       !  equation of state
       rp%gamma = hd_gamma

       ! spatial interval over which to compute solution
       rp%xl = xprobmin1 + 0.5_dp * dxlevel(1)
       rp%xr = xprobmax1 - 0.5_dp * dxlevel(1)

       ! location of discontinuity at t = 0
       rp%xi = 0.5_dp * (rp%xl + rp%xr)
    case default
       error stop "Unknown iprob"
    end select

  end subroutine select_riemann_problem

  ! Initialize one grid
  subroutine rm1d_init_one_grid(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)

    where (abs(x(ixmin1:ixmax1,1))<rp%xi)
       w(ixmin1:ixmax1,rho_)   = rp%rhol
       w(ixmin1:ixmax1,mom(1)) = rp%ul
       w(ixmin1:ixmax1,e_)     = rp%pl
    elsewhere
       w(ixmin1:ixmax1,rho_)   = rp%rhor
       w(ixmin1:ixmax1,mom(1)) = rp%ur
       w(ixmin1:ixmax1,e_)     = rp%pr
    end where

    call hd_to_conserved(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)

  end subroutine rm1d_init_one_grid

  subroutine set_error(iit, time)
    use mod_global_parameters

    integer, intent(in) :: iit
    real(dp), intent(in) :: time

    double precision :: sol_x(domain_nx1)
    double precision :: sol_rho(domain_nx1)
    double precision :: sol_u(domain_nx1)
    double precision :: sol_p(domain_nx1)
    integer :: iigrid, igrid
    integer :: ix1, ix2

    rp%t = time
    call riemann_solve(rp, domain_nx1, sol_x, sol_rho, sol_u, sol_p)

    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ix1 = ceiling((pw(igrid)%x(ixMlo1, 1) - xprobmin1) / dxlevel(1))
       ix2 = ix1 + block_nx1 - 1
       pw(igrid)%w(ixMlo1:ixMhi1, i_sol_rho) = sol_rho(ix1:ix2)
       pw(igrid)%w(ixMlo1:ixMhi1, i_err_rho) = &
            pw(igrid)%w(ixMlo1:ixMhi1, rho_) - &
            pw(igrid)%w(ixMlo1:ixMhi1, i_sol_rho)
    end do

  end subroutine set_error

  subroutine print_error()
    use mod_global_parameters
    use mod_input_output, only: get_volume_average
    double precision   :: modes(nw, 2), volume

    call get_volume_average(1, modes(:, 1), volume)
    call get_volume_average(2, modes(:, 2), volume)

    if (mype == 0) then
       write(*, "(A,2E16.8)") " time -- RMSE:", global_time, sqrt(modes(i_err_rho, 2))
    end if
  end subroutine print_error

end module mod_usr
