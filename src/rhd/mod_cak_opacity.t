module mod_cak_opacity
    implicit NONE

    !> min and max indices for R,T-range in opacity table
    integer, parameter :: Dmin = 2
    integer, parameter :: Dmax = 21
    integer, parameter :: Tmin = 2
    integer, parameter :: Tmax = 21

    !> The opacity tables are read once and stored globally in Kappa_vals
    double precision, public :: alpha_vals(Dmin:Dmax,Tmin:Tmax)
    double precision, public :: Qbar_vals(Dmin:Dmax,Tmin:Tmax)
    double precision, public :: Q0_vals(Dmin:Dmax,Tmin:Tmax)
    double precision, public :: kappa_e_vals(Dmin:Dmax,Tmin:Tmax)

    double precision, public :: Log_D_list(Dmin:Dmax)
    double precision, public :: Log_T_list(Tmin:Tmax)

    character(255), public :: AMRVAC_DIR
    character(255), public :: fileplace

    public :: init_cak
    public :: set_cak_opacity

  contains

!> This routine is called when the fld radiation module is initialised.
!> Here, the tables for different He Abndcs are read and interpolated
subroutine init_cak(tablename)
  ! use mod_global_parameters
  character(6), intent(in) :: tablename

  CALL get_environment_variable("AMRVAC_DIR", AMRVAC_DIR)
  fileplace = TRIM(AMRVAC_DIR)//"/src/rhd/CAK_tables/"

  call read_table(Log_D_list, Log_T_list, alpha_vals, TRIM(tablename)//"/al_TD")
  call read_table(Log_D_list, Log_T_list, Qbar_vals, TRIM(tablename)//"/Qb_TD")
  call read_table(Log_D_list, Log_T_list, Q0_vals, TRIM(tablename)//"/Q0_TD")
  call read_table(Log_D_list, Log_T_list, kappa_e_vals, TRIM(tablename)//"/Ke_TD")

  ! print*, "Read Luka's tables"

end subroutine init_cak

!> This subroutine calculates the opacity for
!> a given temperature-density structure.
!> The opacities are read from a table that has the initialised metalicity
subroutine set_cak_opacity(rho,temp, gradv,alpha_output, Qbar_output, Q0_output, kappa_e_output)
  double precision, intent(in) :: rho, temp, gradv
  double precision, intent(out) :: alpha_output, Qbar_output, Q0_output, kappa_e_output

  double precision, PARAMETER :: const_c     = 2.99792458d10   ! cm S^-1           ; Speed of light

  double precision :: D_input, T_input

  double precision :: kappa_cak
  double precision :: tau, M_t

  D_input = dlog10(rho)
  T_input = dlog10(temp)

  D_input = min(-7.d0-1.d-5, D_input)
  D_input = max(-16.d0+1.d-5, D_input)
  T_input = min(5d0-1.d-5, T_input)
  T_input = max(4d0+1.d-5, T_input)

  call get_val_comb(alpha_vals,Qbar_vals,Q0_vals,kappa_e_vals, &
                    Log_D_list, Log_T_list, D_input, T_input, &
                    alpha_output, Qbar_output, Q0_output, kappa_e_output)

  ! alpha = alpha_output
  ! Qbar = Qbar_output
  ! Q0 = Q0_output
  ! kappa_e = kappa_e_output

end subroutine set_cak_opacity

!> This routine reads out values and arguments from a table
subroutine read_table(D, T, K, filename)
    !> This routine reads in the the values for log kappa, and the values for log T and log R on the x and y axis

    double precision, intent(out) :: K(Dmin:Dmax,Tmin:Tmax), D(Dmin:Dmax), T(Tmin:Tmax)
    character(*), intent(in) :: filename

    character :: dum
    integer :: row, col

    OPEN(1,status = 'old', FILE=TRIM(fileplace)//filename)

    !> Read logT
    READ(1,*) dum,T(Tmin:Tmax)

    !> Read T and K
    do row = Dmin,Dmax !> NOT READING ENTIRE TABLE
      READ(1,*) D(row), K(row,Tmin:Tmax)
    enddo

    CLOSE(1)

end subroutine read_table

!>This subroutine looks in the table for the four couples (T,R)
!surrounding a given input for T and R
subroutine get_val_comb(K1_vals,K2_vals,K3_vals,K4_vals, &
                   Log_D_list, Log_T_list, D, T, &
                   K1, K2, K3, K4)

    double precision, intent(in) :: K1_vals(Dmin:Dmax,Tmin:Tmax)
    double precision, intent(in) :: K2_vals(Dmin:Dmax,Tmin:Tmax)
    double precision, intent(in) :: K3_vals(Dmin:Dmax,Tmin:Tmax)
    double precision, intent(in) :: K4_vals(Dmin:Dmax,Tmin:Tmax)

    double precision, intent(in) :: Log_D_list(Dmin:Dmax)
    double precision, intent(in) :: Log_T_list(Tmin:Tmax)
    double precision, intent(in) :: D, T

    double precision, intent(out) :: K1
    double precision, intent(out) :: K2
    double precision, intent(out) :: K3
    double precision, intent(out) :: K4

    integer :: low_r_index, up_r_index
    integer :: low_t_index, up_t_index

    if (D .gt. maxval(Log_D_list)) then
        ! print*, 'Extrapolating in logR'
        low_r_index = Dmax-1
        up_r_index = Dmax
    elseif (D .lt. minval(Log_D_list)) then
        ! print*, 'Extrapolating in logR'
        low_r_index = Dmin
        up_r_index = Dmin+1
    else
        call get_low_up_index(D, Log_D_list, Dmin, Dmax, low_r_index, up_r_index)
    endif

    if (T .gt. maxval(Log_T_list)) then
        ! print*, 'Extrapolating in logT'
        low_t_index = Tmax-1
        up_t_index = Tmax
    elseif ( T .lt. minval(Log_T_list)) then
        ! print*, 'Extrapolating in logT'
        low_t_index = Tmin
        up_t_index = Tmin+1
    else
        call get_low_up_index(T, Log_T_list, Tmin, Tmax, low_t_index, up_t_index)
    endif

    call interpolate_KRT(low_r_index, up_r_index, low_t_index, up_t_index, &
                       Log_D_list, Log_T_list, K1_vals, D, T, K1)

    call interpolate_KRT(low_r_index, up_r_index, low_t_index, up_t_index, &
                       Log_D_list, Log_T_list, K2_vals, D, T, K2)

    call interpolate_KRT(low_r_index, up_r_index, low_t_index, up_t_index, &
                       Log_D_list, Log_T_list, K3_vals, D, T, K3)

    call interpolate_KRT(low_r_index, up_r_index, low_t_index, up_t_index, &
                       Log_D_list, Log_T_list, K4_vals, D, T, K4)
end subroutine get_val_comb


!> this subroutine finds the indexes in R and T arrays of the two values surrounding the input R and T
subroutine get_low_up_index(x, x_list, imin, imax, low_i, up_i)

    integer, intent(in) :: imin, imax
    double precision, intent(in) :: x
    double precision, intent(in) :: x_list(imin:imax)

    integer, intent(out) :: low_i, up_i

    double precision :: low_val, up_val
    double precision :: res(imin:imax)

    res = x_list - x

    up_val = minval(res, MASK = res .ge. 0) + x
    low_val = maxval(res, MASK = res .le. 0) + x

    up_i = minloc(abs(x_list - up_val),1) + imin -1
    low_i = minloc(abs(x_list - low_val),1) + imin -1


    if (up_i .eq. low_i) low_i = low_i - 1

end subroutine get_low_up_index

!> This subroutine does a bilinear interpolation in the R,T-plane
subroutine interpolate_KRT(low_r, up_r, low_t, up_t, Log_D_list, Log_T_list, Kappa_vals, D, T, k_interp)

    integer, intent(in) :: low_r, up_r, low_t, up_t
    double precision, intent(in) :: Kappa_vals(Dmin:Dmax,Tmin:Tmax)
    double precision, intent(in) :: Log_D_list(Dmin:Dmax)
    double precision, intent(in) :: Log_T_list(Tmin:Tmax)
    double precision, intent(in) :: D,T
    double precision, intent(out) :: k_interp

    double precision :: r1,r2,t1,t2
    double precision :: k1, k2, k3, k4
    double precision :: ka, kb

    !Cool ascii drawing of interpolation scheme: first interpolate twice in the T coord to get
    !ka and kb, then interpolate in the R coord to get ki

!   r_1    R        r_2
!     |                |
!     |                |
! ----k1--------ka-----k2----- t_1
!     |          |     |
!     |          |     |
!   T |          |     |
!     |          |     |
!     |          ki    |
!     |          |     |
! ----k3--------kb-----k4----- t_2
!     |                |
!     |                |


    r1 = Log_D_list(low_r)
    r2 = Log_D_list(up_r)
    t1 = Log_T_list(low_t)
    t2 = Log_T_list(up_t)

    ! k1 = Kappa_vals(low_t, low_r)
    ! k2 = Kappa_vals(low_t, up_r)
    ! k3 = Kappa_vals(up_t, low_r)
    ! k4 = Kappa_vals(up_t, up_r)

    k1 = Kappa_vals(low_r, low_t)
    k2 = Kappa_vals(low_r, up_t)
    k3 = Kappa_vals(up_r, low_t)
    k4 = Kappa_vals(up_r, up_t)

    ! print*, 'surounding indexes'
    ! print*, low_r, up_r, low_t, up_t
    ! print*, 'surounding input values'
    ! print*, r1,r2,t1,t2
    ! print*, 'surounding table values'
    ! print*, k1, k2, k3, k4

    call interpolate1D(r1,r2,D,k1,k2,ka)
    call interpolate1D(r1,r2,D,k3,k4,kb)
    call interpolate1D(t1,t2,T,ka,kb,k_interp)

    ! print*, 'interpolated value'
    ! print*, ka, kb, k_interp

end subroutine interpolate_KRT


!> Interpolation in one dimension
subroutine interpolate1D(x1, x2, x, y1, y2, y)

    double precision, intent(in) :: x, x1, x2
    double precision, intent(in) :: y1, y2
    double precision, intent(out) :: y

    y = y1 + (x-x1)*(y2-y1)/(x2-x1)

end subroutine interpolate1D

end module mod_cak_opacity

! !> Interpolation on logarithmic scale
! subroutine log_interpolate1D(x1, x2, x, y1, y2, y)
!
!    double precision, intent(in) :: x, x1, x2
!    double precision, intent(in) :: y1, y2
!    double precision, intent(out) :: y
!
!    double precision :: expx, expx1, expx2
!    double precision :: expy1, expy2
!
!    expx = 10**x
!    expx1 = 10**x1
!    expx2 = 10**x2
!
!    expy1 = 10**y1
!    expy2 = 10**y2
!
!    y = expy1 + (expx-expx1)*(expy2-expy1)/(expx2-expx1)
!    y = log10(y)
!
! end subroutine log_interpolate1D
