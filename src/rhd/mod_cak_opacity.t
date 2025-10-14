!> This module reads in CAK line opacities in the Gayley (1995) notation
!> (alpha, Qbar, Q0, kappae) from corresponding tables. Tabulated values assume
!> LTE conditions and are a function of mass density (D) and temperature (T),
!> which are both given in base 10 logarithm.
!> The construction of the tables is outlined in Poniatowski+ (2021), A&A, 667.

module mod_cak_opacity

  use mod_comm_lib, only: mpistop

  implicit none
  private

  !> min and max indices for R,T-range in opacity table
  integer, parameter :: iDmin = 2, iDmax = 51
  integer, parameter :: iTmin = 2, iTmax = 21

  !> If user wants to read tables from different location than AMRVAC source
  !> NOTE: the tables should have same size as AMRVAC tables
  logical :: set_custom_tabledir = .false.

  !> The opacity tables are read once and stored globally
  double precision, public :: alpha_vals(iDmin:iDmax,iTmin:iTmax)
  double precision, public :: Qbar_vals(iDmin:iDmax,iTmin:iTmax)
  double precision, public :: Q0_vals(iDmin:iDmax,iTmin:iTmax)
  double precision, public :: kappae_vals(iDmin:iDmax,iTmin:iTmax)
  double precision, public :: logD_list(iDmin:iDmax), logT_list(iTmin:iTmax)

  public :: init_cak_table
  public :: set_cak_opacity

contains

  !> This routine is called when the FLD radiation module is initialised.
  subroutine init_cak_table(tabledir,set_custom_tabledir)

    character(len=*),  intent(in) :: tabledir
    logical, optional, intent(in) :: set_custom_tabledir

    ! Local variables
    character(len=256) :: AMRVAC_DIR, path_table_dir

    if (present(set_custom_tabledir) .and. set_custom_tabledir) then
      path_table_dir = trim(tabledir)
    else
      call get_environment_variable("AMRVAC_DIR", AMRVAC_DIR)
      path_table_dir = trim(AMRVAC_DIR)//"/src/tables/CAK_tables/"//trim(tabledir)
    endif

    call read_table(logD_list, logT_list, alpha_vals, trim(path_table_dir)//"/al_TD")
    call read_table(logD_list, logT_list, Qbar_vals, trim(path_table_dir)//"/Qb_TD")
    call read_table(logD_list, logT_list, Q0_vals, trim(path_table_dir)//"/Q0_TD")
    call read_table(logD_list, logT_list, kappae_vals, trim(path_table_dir)//"/Ke_TD")

  end subroutine init_cak_table

  !> This subroutine calculates the opacity for a given temperature-density
  !> structure. Opacities are read from a table with given metalicity.
  subroutine set_cak_opacity(rho, temp, alpha_output, Qbar_output, Q0_output, kappae_output)

    double precision, intent(in)  :: rho, temp
    double precision, intent(out) :: alpha_output, Qbar_output, Q0_output, kappae_output

    ! Local variables
    double precision :: D_input, T_input, kappa_cak

    if (temp <= 0.0d0) call mpistop('Input temperature < 0 in set_cak_opacity.')
    if (rho <= 0.0d0) call mpistop('Input density < 0 in set_cak_opacity.')

    D_input = log10(rho)
    T_input = log10(temp)

    call get_val_comb(alpha_vals,Qbar_vals,Q0_vals,kappae_vals, &
                      logD_list, logT_list, D_input, T_input, &
                      alpha_output, Qbar_output, Q0_output, kappae_output)

  end subroutine set_cak_opacity

  !> This routine reads in 1-D (T,rho) values from a line opacity table and gives
  !> back as output the 1-D (T,rho) values and 2-D line opacity value
  subroutine read_table(D, T, K, filename)

    character(*), intent(in)      :: filename
    double precision, intent(out) :: K(iDmin:iDmax,iTmin:iTmax), D(iDmin:iDmax), T(iTmin:iTmax)

    ! Local variables
    character :: dum
    integer   :: row, col, funit=99
    logical   :: alive

    inquire(file=trim(filename), exist=alive)

    if (alive) then
      open(funit, file=trim(filename), status='unknown')
    else
      call mpistop('Table file you want to use cannot be found: '//filename)
    endif

    ! Read temperature
    read(funit,*) dum, T(iTmin:iTmax)

    ! Read rho and kappa
    do row = iDmin,iDmax
      read(funit,*) D(row), K(row,iTmin:iTmax)
    enddo

    close(funit)

  end subroutine read_table

  !> This subroutine looks in the table for the four couples (T,rho) surrounding
  !> a given input T and rho
  subroutine get_val_comb(K1_vals,K2_vals,K3_vals,K4_vals, &
                     logD_list, logT_list, D, T, K1, K2, K3, K4)

    double precision, intent(in)  :: K1_vals(iDmin:iDmax,iTmin:iTmax)
    double precision, intent(in)  :: K2_vals(iDmin:iDmax,iTmin:iTmax)
    double precision, intent(in)  :: K3_vals(iDmin:iDmax,iTmin:iTmax)
    double precision, intent(in)  :: K4_vals(iDmin:iDmax,iTmin:iTmax)
    double precision, intent(in)  :: logD_list(iDmin:iDmax), logT_list(iTmin:iTmax)
    double precision, intent(in)  :: D, T
    double precision, intent(out) :: K1, K2, K3, K4

    ! Local variables
    integer :: low_D_index, up_D_index, low_T_index, up_T_index

    if (D > maxval(logD_list)) then
      ! print*, 'Extrapolating in logR'
      low_D_index = iDmax-1
      up_D_index = iDmax
    elseif (D < minval(logD_list)) then
      ! print*, 'Extrapolating in logR'
      low_D_index = iDmin
      up_D_index = iDmin+1
    else
      call get_low_up_index(D, logD_list, iDmin, iDmax, low_D_index, up_D_index)
    endif

    if (T > maxval(logT_list)) then
      ! print*, 'Extrapolating in logT'
      low_T_index = iTmax-1
      up_T_index = iTmax
    elseif (T < minval(logT_list)) then
      ! print*, 'Extrapolating in logT'
      low_T_index = iTmin
      up_T_index = iTmin+1
    else
      call get_low_up_index(T, logT_list, iTmin, iTmax, low_T_index, up_T_index)
    endif

    call interpolate_KRT(low_D_index, up_D_index, low_T_index, up_T_index, &
                       logD_list, logT_list, K1_vals, D, T, K1)

    call interpolate_KRT(low_D_index, up_D_index, low_T_index, up_T_index, &
                       logD_list, logT_list, K2_vals, D, T, K2)

    call interpolate_KRT(low_D_index, up_D_index, low_T_index, up_T_index, &
                       logD_list, logT_list, K3_vals, D, T, K3)

    call interpolate_KRT(low_D_index, up_D_index, low_T_index, up_T_index, &
                       logD_list, logT_list, K4_vals, D, T, K4)

  end subroutine get_val_comb

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
    up_val = minval(res, MASK = res >= 0) + var_in
    low_val = maxval(res, MASK = res <= 0) + var_in

    ! Find all bounding indices
    up_i = minloc(abs(var_list - up_val),1) + imin -1
    low_i = minloc(abs(var_list - low_val),1) + imin -1

    if (up_i == low_i) low_i = low_i - 1

  end subroutine get_low_up_index

  !> This subroutine does a bilinear interpolation in the R,T-plane
  subroutine interpolate_KRT(low_r, up_r, low_t, up_t, logD_list, logT_list, Kappa_vals, D, T, kappa_interp)

    integer, intent(in)           :: low_r, up_r, low_t, up_t
    double precision, intent(in)  :: Kappa_vals(iDmin:iDmax,iTmin:iTmax)
    double precision, intent(in)  :: logD_list(iDmin:iDmax), logT_list(iTmin:iTmax)
    double precision, intent(in)  :: D, T
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

    r1 = logD_list(low_r)
    r2 = logD_list(up_r)
    t1 = logT_list(low_t)
    t2 = logT_list(up_t)

    k1 = Kappa_vals(low_r, low_t)
    k2 = Kappa_vals(low_r, up_t)
    k3 = Kappa_vals(up_r, low_t)
    k4 = Kappa_vals(up_r, up_t)

    call interpolate1D(r1,r2,D,k1,k2,ka)
    call interpolate1D(r1,r2,D,k3,k4,kb)
    call interpolate1D(t1,t2,T,ka,kb,kappa_interp)

  end subroutine interpolate_KRT

  !> Interpolation in one dimension
  subroutine interpolate1D(x1, x2, x, y1, y2, y)
        
    double precision, intent(in)  :: x, x1, x2, y1, y2
    double precision, intent(out) :: y

    y = y1 + (x-x1)*(y2-y1)/(x2-x1)

  end subroutine interpolate1D

end module mod_cak_opacity
