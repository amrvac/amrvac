! Semirelativistic wave tests of Gombosi et al. 2002 JCoPh 177, 176, section 7,
! set up so that each one can be run with or without splitting off an analytic
! background field (section 6, eqs 91-97)
!
! The initial state is prescribed in terms of the total magnetic field. When
! B0field is on, the analytic background is subtracted so that the total field
! is identical in every splitting mode, and the solution must be too
module mod_usr
  use mod_mhd
  implicit none

  ! uniform background state, defaults are the section 7.2 Alfven wave test
  double precision :: rho_ini = 1.0d0
  double precision :: p_ini   = 0.1d0
  double precision :: u_ini(3) = (/ 0.5d0, 0.0d0, 0.0d0 /)
  double precision :: b_ini(3) = (/ 1.0d0, 1.0d0, 0.0d0 /)

  ! velocity jump, added for x<0 and subtracted for x>0
  double precision :: du_jump(3) = (/ 0.0d0, 0.0d0, 0.01d0 /)

  ! analytic background field: 0 leaves it uniform so that the split field and
  ! the background current both vanish, 1 shears the y component so that
  ! neither does and the j0xB0 and E.j0 source terms are exercised
  integer          :: b0_profile = 0
  double precision :: b0_shear   = 0.5d0
  double precision :: b0_width   = 50.0d0

  ! unit system, which fixes the true light speed relative to the reduced one (see gombosi.par)
  double precision :: unit_v = 1.0d0, unit_nd = 1.0d0

contains

  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ rho_ini, p_ini, u_ini, b_ini, du_jump, &
                        b0_profile, b0_shear, b0_width, unit_v, unit_nd

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine usr_init()

    call usr_params_read(par_files)

    unit_velocity      = unit_v
    unit_numberdensity = unit_nd

    usr_init_one_grid => initonegrid_usr
    usr_set_B0        => specialset_B0

    call set_coordinate_system("Cartesian_1.75D")
    call mhd_activate()

  end subroutine usr_init

  !> Analytic background field, uniform unless b0_profile selects the shear
  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    integer :: idir

    do idir = 1, ndir
      wB0(ixO^S,idir) = b_ini(idir)
    end do

    if (b0_profile == 1) then
      wB0(ixO^S,2) = wB0(ixO^S,2) + b0_shear*dtanh(x(ixO^S,1)/b0_width)
    end if

  end subroutine specialset_B0

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: wB0loc(ixI^S,1:ndir)
    integer          :: idir

    w(ixO^S,rho_) = rho_ini
    w(ixO^S,p_)   = p_ini

    do idir = 1, ndir
      where (x(ixO^S,1) < 0.0d0)
        w(ixO^S,mom(idir)) = u_ini(idir) + du_jump(idir)
      else where
        w(ixO^S,mom(idir)) = u_ini(idir) - du_jump(idir)
      end where
    end do

    ! total field, reduced to the split field where the background is carried
    ! separately, so that every splitting mode starts from the same state
    do idir = 1, ndir
      w(ixO^S,mag(idir)) = b_ini(idir)
    end do
    if (B0field) then
      call specialset_B0(ixI^L,ixO^L,x,wB0loc)
      do idir = 1, ndir
        w(ixO^S,mag(idir)) = w(ixO^S,mag(idir)) - wB0loc(ixO^S,idir)
      end do
    end if

    call eos%to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

end module mod_usr
