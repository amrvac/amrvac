!> This module reads in Rosseland-mean opacities from OPAL tables. Table opacity
!> values are given in base 10 logarithm and are a function of mass density (R)
!> and temperature (T), which are also both given in base 10 logarithm.

module mod_opal_opacity
  
  implicit none

  !> min and max indices for R,T-range in opacity table
  integer, parameter :: iRmin = 2, iRmax = 20
  integer, parameter :: iTmin = 7, iTmax = 76

  !> The opacity tables are read once and stored globally in Kappa_vals
  double precision, public :: Kappa_vals(iTmin:iTmax,iRmin:iRmax)
  double precision, public :: logR_list(iRmin:iRmax), logT_list(iTmin:iTmax)

  public :: init_opal_table
  public :: set_opal_opacity

contains

  !> This subroutine is called when the FLD radiation module is initialised.
  !> The OPAL tables for different helium abundances are read and interpolated.
  subroutine init_opal_table(He_abundance, tablename)
    
    double precision, intent(in) :: He_abundance
    character(len=*), intent(in) :: tablename

    ! Local variables
    character(len=:), allocatable :: AMRVAC_DIR, path_table_dir
    
    call get_environment_variable("AMRVAC_DIR", AMRVAC_DIR)
        
    ! Absolute path to OPAL table location
    path_table_dir = trim(AMRVAC_DIR)//"/src/rhd/OPAL_tables/"//trim(tablename)
    
    ! Read in the OPAL table
    call read_table(logR_list, logT_list, Kappa_vals, path_table_dir)

  end subroutine init_opal_table

  !> This subroutine calculates the opacity for a given temperature-density
  !> structure. Opacities are read from a table with given metalicity.
  subroutine set_opal_opacity(rho, temp, kappa)
        
    double precision, intent(in)  :: rho, temp
    double precision, intent(out) :: kappa

    ! Local variables
    double precision :: logR_in, logT_in, logKappa_out

    ! Convert input density and temperature to OPAL table convention
    logR_in = log10(rho/(temp*1d-6)**3)
    logT_in = log10(temp)

    call get_kappa(Kappa_vals, logR_list, logT_list, logR_in, logT_in, logKappa_out)

    ! If the outcome is 9.999 (no table entry), look right in the table
    do while (logKappa_out .gt. 9.0d0)
        logR_in = logR_in + 0.5d0
        call get_kappa(Kappa_vals, logR_list, logT_list, logR_in, logT_in, logKappa_out)
    enddo

    ! If the outcome is NaN, look left in the table
    do while (logKappa_out .eq. 0.0d0)
        logR_in = logR_in - 0.5d0
        call get_kappa(Kappa_vals, logR_list, logT_list, logR_in, logT_in, logKappa_out)
    enddo

    ! Convert OPAL opacity for output
    kappa = 10d0**logKappa_out

  end subroutine set_opal_opacity

  !> This routine reads in 1-D (rho,T) values from an opacity table and gives
  !> back as output the 1-D (rho,T) values and 2-D opacity
  subroutine read_table(R, T, K, filename)
        
    character(*), intent(in)      :: filename
    double precision, intent(out) :: K(iTmin:iTmax,iRmin:iRmax), R(iRmin:iRmax), T(iTmin:iTmax)

    ! Local variables
    character :: dum
    integer   :: row, col

    open(unit=1, status='old', file=trim(filename))

    ! Skip first 4 lines
    do row = 1,4
        read(1,*)
    enddo

    ! Read rho
    read(1,*) dum, R(iRmin:iRmax)
    read(1,*)

    ! Read T and kappa
    do row = iTmin,iTmax
        read(1,'(f4.2,19f7.3)') T(row), K(row,iRmin:iRmax)
    enddo

    close(1)

  end subroutine read_table

  !> This subroutine looks in the table for the four couples (T,rho) surrounding
  !> a given input T and rho
  subroutine get_kappa(Kappa_vals, logR_list, logT_list, R, T, K)
    
    double precision, intent(in)  :: Kappa_vals(iTmin:iTmax,iRmin:iRmax)
    double precision, intent(in)  :: logR_list(iRmin:iRmax), logT_list(iTmin:iTmax)
    double precision, intent(in)  :: R, T
    double precision, intent(out) :: K

    ! Local variables
    integer :: low_R_index, up_R_index
    integer :: low_T_index, up_T_index

    if (R .gt. maxval(logR_list)) then
        ! print*, 'Extrapolating in logR'
        low_R_index = iRmax-1
        up_R_index = iRmax
    elseif (R .lt. minval(logR_list)) then
        ! print*, 'Extrapolating in logR'
        low_R_index = iRmin
        up_R_index = iRmin+1
    else
        call get_low_up_index(R, logR_list, iRmin, iRmax, low_R_index, up_R_index)
    endif

    if (T .gt. maxval(logT_list)) then
        ! print*, 'Extrapolating in logT'
        low_T_index = iTmax-1
        up_T_index = iTmax
    elseif ( T .lt. minval(logT_list)) then
        ! print*, 'Extrapolating in logT'
        low_T_index = iTmin
        up_T_index = iTmin+1
    else
        call get_low_up_index(T, logT_list, iTmin, iTmax, low_T_index, up_T_index)
    endif

    call interpolate_KRT(low_R_index, up_R_index, low_T_index, up_T_index, logR_list, logT_list, Kappa_vals, R, T, K)

  end subroutine get_kappa

  !> This subroutine finds the indices in rho and T arrays of the two values
  !> surrounding the input rho and T
  subroutine get_low_up_index(var_in, var_list, imin, imax, low_i, up_i)
    
    integer, intent(in)          :: imin, imax
    double precision, intent(in) :: var_in, var_list(imin:imax)
    integer, intent(out)         :: low_i, up_i

    ! Local variables
    double precision :: low_val, up_val, res(imin:imax)

    res = var_list - var_in
    
    ! Find all bounding values for given input
    up_val = minval(res, MASK = res .ge. 0) + var_in
    low_val = maxval(res, MASK = res .le. 0) + var_in

    ! Find all bounding indices
    up_i = minloc(abs(var_list - up_val),1) + imin -1
    low_i = minloc(abs(var_list - low_val),1) + imin -1

    if (up_i .eq. low_i) low_i = low_i - 1

  end subroutine get_low_up_index

  !> This subroutine does a bilinear interpolation in the R,T-plane
  subroutine interpolate_KRT(low_r, up_r, low_t, up_t, logR_list, logT_list, Kappa_vals, R, T, kappa_interp)
    
    integer, intent(in)           :: low_r, up_r, low_t, up_t
    double precision, intent(in)  :: Kappa_vals(iTmin:iTmax,iRmin:iRmax)
    double precision, intent(in)  :: logR_list(iRmin:iRmax), logT_list(iTmin:iTmax)
    double precision, intent(in)  :: R, T
    double precision, intent(out) :: kappa_interp

    ! Local variables
    double precision :: r1, r2, t1, t2, k1, k2, k3, k4, ka, kb

    ! First interpolate twice in the T coord to get ka and kb, then interpolate
    ! in the R coord to get ki

    !    r_1    rho       r_2
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

    r1 = logR_list(low_r)
    r2 = logR_list(up_r)
    t1 = logT_list(low_t)
    t2 = logT_list(up_t)

    k1 = Kappa_vals(low_t, low_r)
    k2 = Kappa_vals(low_t, up_r)
    k3 = Kappa_vals(up_t, low_r)
    k4 = Kappa_vals(up_t, up_r)

    call interpolate1D(r1,r2,R,k1,k2,ka)
    call interpolate1D(r1,r2,R,k3,k4,kb)
    call interpolate1D(t1,t2,T,ka,kb,kappa_interp)

  end subroutine interpolate_KRT

  !> Interpolation in one dimension
  subroutine interpolate1D(x1, x2, x, y1, y2, y)
        
    double precision, intent(in)  :: x, x1, x2, y1, y2
    double precision, intent(out) :: y

    y = y1 + (x-x1)*(y2-y1)/(x2-x1)

  end subroutine interpolate1D

end module mod_opal_opacity
