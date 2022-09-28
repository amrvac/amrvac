!> This module reads in opacities from opal tables.

module mod_opal_opacity
    implicit NONE

    !> min and max indices for R,T-range in opacity table
    integer, parameter :: rmin = 2
    integer, parameter :: rmax = 20
    integer, parameter :: tmin = 7
    integer, parameter :: tmax = 76

    !> The opacity tables are read once and stored globally in Kappa_vals
    double precision, public :: Kappa_vals(7:76,2:20)
    double precision, public :: Kappa_vals1(7:76,2:20)
    double precision, public :: Kappa_vals2(7:76,2:20)

    double precision, public :: Log_R_list(2:20)
    double precision, public :: Log_T_list(7:76)

    character(255), public :: AMRVAC_DIR
    character(255), public :: fileplace
    ! character(*), parameter, public :: AMRVAC_DIR = "/STER/luka/codes/Nico_amrvac/" ! use call getenv("AMRVAC_DIR", AMRVAC_DIR)
    ! character(*), parameter, public :: fileplace = AMRVAC_DIR//"src/rhd/Opacity_tables/"

    ! character(len=:), allocatable :: AMRVAC_DIR != "/STER/luka/codes/Nico_amrvac" ! use call getenv("AMRVAC_DIR", AMRVAC_DIR)
    ! character, allocatable :: fileplace != AMRVAC_DIR//"src/rhd/Opacity_tables/"


    public :: init_opal
    public :: set_opal_opacity

  contains

!> This routine is called when the fld radiation module is initialised.
!> Here, the tables for different He Abndcs are read and interpolated
subroutine init_opal(He_abundance,tablename)
  ! use mod_global_parameters

  double precision, intent(in) :: He_abundance
  character(6), intent(in) :: tablename

  !> Y1 actually 1.0, NOT 0.1???
  double precision :: Y1 = 0.1000
  double precision :: Y2 = 0.0999


  CALL get_environment_variable("AMRVAC_DIR", AMRVAC_DIR)
  fileplace = TRIM(AMRVAC_DIR)//"/src/rhd/Opacity_tables/"


  ! character(len=std_len) :: DIR_tmp
  ! character(len=std_len) :: FILE_tmp
  !
  !
  ! call GET_ENVIRONMENT_VARIABLE("AMRVAC_DIR", DIR_tmp)
  ! ! allocate(AMRVAC_DIR(len=len_trim(DIR_tmp)))
  ! print*, DIR_tmp
  ! print*, TRIM(DIR_tmp)
  !
  ! allocate(character(len=len_trim(DIR_tmp)) :: AMRVAC_DIR)
  !
  ! print*, AMRVAC_DIR
  !
  ! FILE_tmp = AMRVAC_DIR//"src/rhd/Opacity_tables/"
  ! print*, FILE_tmp
  ! ! allocate(fileplace(len=len_trim(FILE_tmp)))
  ! fileplace = trim(FILE_tmp)
  !
  ! print*, '---',fileplace,'---'
  ! stop

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! if (Y1 .gt. Y2) then
  !   if (He_abundance .gt. Y1) call mpistop('OPAL table not covered')
  !   if (He_abundance .lt. Y2) call mpistop('OPAL table not covered')
  ! else
  !   if (He_abundance .lt. Y1) call mpistop('OPAL table not covered')
  !   if (He_abundance .gt. Y2) call mpistop('OPAL table not covered')
  ! endif

  !> WR
  ! call read_table(Log_R_list, Log_T_list, Kappa_vals1,'Y09800')

  !> RSG
  ! call read_table(Log_R_list, Log_T_list, Kappa_vals1,'Y00999')

  !> General
  call read_table(Log_R_list, Log_T_list, Kappa_vals1,tablename)

  ! call read_table(Log_R_list, Log_T_list, Kappa_vals1,'Y01000')
  ! call read_table(Log_R_list, Log_T_list, Kappa_vals2,'Y00999')

  ! if (He_abundance .eq. Y1) then
  !   Kappa_vals = Kappa_vals1
  !   print*, 'Using correct table'
  !   stop
  ! elseif (He_abundance .eq. Y2) then
  !   Kappa_vals = Kappa_vals2
  ! else
  !   call interpolate_two_tables(Y1,Y2, He_abundance, Kappa_vals1, Kappa_vals2, Kappa_vals)
  ! endif

  Kappa_vals = Kappa_vals1

end subroutine init_opal

!> This subroutine calculates the opacity for
!> a given temperature-density structure.
!> The opacities are read from a table that has the initialised metalicity
subroutine set_opal_opacity(rho,temp,kappa)
  double precision, intent(in) :: rho, temp
  double precision, intent(out) :: kappa

  double precision :: R_input, T_input, K_output

  R_input = rho/(temp*1d-6)**3
  T_input = temp

  R_input = dlog10(R_input)
  T_input = dlog10(T_input)

  call get_kappa(Kappa_vals, Log_R_list, Log_T_list, R_input, T_input, K_output)

  !> If the outcome is 9.999, look right in the table
  do while (K_output .gt. 9.0d0)
      ! print*, 'R,T datapoint out of opal table'
      R_input = R_input + 0.5
      call get_kappa(Kappa_vals, Log_R_list, Log_T_list, R_input, T_input, K_output)
  enddo

  !> If the outcome is NaN, look left in the table
  do while (K_output .eq. 0.0d0)
      ! print*, 'R,T datapoint out of opal table'
      R_input = R_input - 0.5d0
      call get_kappa(Kappa_vals, Log_R_list, Log_T_list, R_input, T_input, K_output)
  enddo

  kappa = 10d0**K_output

end subroutine set_opal_opacity

!> This routine reads out values and arguments from an opacity table
subroutine read_table(R, T, K, filename)
    !> This routine reads in the the values for log kappa, and the values for log T and log R on the x and y axis

    double precision, intent(out) :: K(7:76,2:20), R(2:20), T(7:76)
    character(*), intent(in) :: filename

    character :: dum
    integer :: row, col

    OPEN(1,status = 'old', FILE=TRIM(fileplace)//filename)

    !> Skip first 4 lines
    do row = 1,4
        READ(1,*)
    enddo

    !> Read R
    READ(1,*) dum,R(2:20)

    READ(1,*)

    !> Read T and K
    do row = 7,76 !> NOT READING ENTIRE TABLE
        READ(1,'(f4.2,19f7.3)') T(row), K(row,2:20)
    enddo

    CLOSE(1)

end subroutine read_table


!> This subroutine creates a new table for a given He abundance,
! by interpolating to known tables at every point in the R,T plane
subroutine interpolate_two_tables(Y1, Y2, Y_in, K1, K2, K_interp)

    double precision, intent(in) :: K1(7:76,2:20), K2(7:76,2:20)
    double precision, intent(in) :: Y1, Y2, Y_in
    double precision, intent(out) :: K_interp(7:76,2:20)

    integer row, colum

    do colum=2,20
    do row=7,76
        call interpolate1D(Y1,Y2,Y_in,K1(row,colum),K2(row,colum),K_interp(row,colum))
    enddo
    enddo

end subroutine interpolate_two_tables


!>This subroutine looks in the table for the four couples (T,R)
!surrounding a given input for T and R
subroutine get_kappa(Kappa_vals, Log_R_list, Log_T_list, R, T, K)

    double precision, intent(in) :: Kappa_vals(7:76,2:20)
    double precision, intent(in) :: Log_R_list(2:20)
    double precision, intent(in) :: Log_T_list(7:76)

    double precision, intent(in) :: R, T
    double precision, intent(out) :: K

    integer :: low_r_index, up_r_index
    integer :: low_t_index, up_t_index

    if (R .gt. maxval(Log_R_list)) then
        ! print*, 'Extrapolating in logR'
        low_r_index = 19
        up_r_index = 20
    elseif (R .lt. minval(Log_R_list)) then
        ! print*, 'Extrapolating in logR'
        low_r_index = 2
        up_r_index = 3
    else
        call get_low_up_index(R, Log_R_list, 2, 20, low_r_index, up_r_index)
    endif

    if (T .gt. maxval(Log_T_list)) then
        ! print*, 'Extrapolating in logT'
        low_t_index = 75
        up_t_index = 76
    elseif ( T .lt. minval(Log_T_list)) then
        ! print*, 'Extrapolating in logT'
        low_t_index = 7
        up_t_index = 8
    else
        call get_low_up_index(T, Log_T_list, 7, 76, low_t_index, up_t_index)
    endif

    call interpolate_KRT(low_r_index, up_r_index, low_t_index, up_t_index, Log_R_list, Log_T_list, Kappa_vals, R, T, K)

end subroutine get_kappa


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
subroutine interpolate_KRT(low_r, up_r, low_t, up_t, Log_R_list, Log_T_list, Kappa_vals, R, T, k_interp)

    integer, intent(in) :: low_r, up_r, low_t, up_t
    double precision, intent(in) :: Kappa_vals(7:76,2:20)
    double precision, intent(in) :: Log_R_list(2:20)
    double precision, intent(in) :: Log_T_list(7:76)
    double precision, intent(in) :: R,T
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


    r1 = Log_R_list(low_r)
    r2 = Log_R_list(up_r)
    t1 = Log_T_list(low_t)
    t2 = Log_T_list(up_t)

    k1 = Kappa_vals(low_t, low_r)
    k2 = Kappa_vals(low_t, up_r)
    k3 = Kappa_vals(up_t, low_r)
    k4 = Kappa_vals(up_t, up_r)

    call interpolate1D(r1,r2,R,k1,k2,ka)
    call interpolate1D(r1,r2,R,k3,k4,kb)
    call interpolate1D(t1,t2,T,ka,kb,k_interp)

end subroutine interpolate_KRT


!> Interpolation in one dimension
subroutine interpolate1D(x1, x2, x, y1, y2, y)

    double precision, intent(in) :: x, x1, x2
    double precision, intent(in) :: y1, y2
    double precision, intent(out) :: y

    y = y1 + (x-x1)*(y2-y1)/(x2-x1)

end subroutine interpolate1D

end module mod_opal_opacity

!> Interpolation on logarithmic scale
subroutine log_interpolate1D(x1, x2, x, y1, y2, y)

   double precision, intent(in) :: x, x1, x2
   double precision, intent(in) :: y1, y2
   double precision, intent(out) :: y

   double precision :: expx, expx1, expx2
   double precision :: expy1, expy2

   expx = 10**x
   expx1 = 10**x1
   expx2 = 10**x2

   expy1 = 10**y1
   expy2 = 10**y2

   y = expy1 + (expx-expx1)*(expy2-expy1)/(expx2-expx1)
   y = log10(y)

end subroutine log_interpolate1D
