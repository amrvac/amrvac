!> thermal instability test problem copied from amrvac/tests/demo/thermal_instability_HD/
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision :: scale, radius

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ scale, radius

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read


  subroutine usr_init()

    call set_coordinate_system("Cartesian_3D")
    call usr_params_read(par_files)

    unit_length = 1.d9 !< cm
    unit_temperature   = 1.d6 !< K
    unit_numberdensity = 1.d9 !< cm^-3

    usr_init_one_grid => initonegrid_usr

    call phys_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
      ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
      ixImin3:ixImax3,1:nw)

    double precision :: xmid1, xmid2, xmid3
    logical :: mask(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    logical::          first
    data first/.true./


    ! print out the info
    if (first .and. mype==0) then
      print *,'ffHD test problem thermal instability'
      print *, "density contrast ",scale," input"
      print *,"    within radius ",radius," input"
      print *,"derived units: ",unit_pressure,unit_density,unit_velocity
      first=.false.
    end if

    xmid1 = xprobmin1 + 0.5d0 * (xprobmax1 - xprobmin1)
    xmid2 = xprobmin2 + 0.5d0 * (xprobmax2 - xprobmin2)
    xmid3 = xprobmin3 + 0.5d0 * (xprobmax3 - xprobmin3)

    mask(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
     ((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)-xmid1)**2 + &
      (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)-xmid2)**2 + &
      (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)-xmid3)**2) < radius**2

    where(mask(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)  =1.0d0+scale
    elsewhere
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)  =1.0d0
    endwhere
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1)) = zero
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)    = 1.0d0/phys_gamma  

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, q_) = zero

    call phys_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

    if(mype==0.and.first)then
       write(*,*)'Doing TI setup with gamma=',phys_gamma
       first=.false.
    endif

  end subroutine initonegrid_usr


  pure real(dp) function gravity_field(wCT, x, idim) result(field)
    !$acc routine seq
    real(dp), intent(in)    :: wCT(nw_phys)
    real(dp), intent(in)    :: x(1:ndim)
    integer, value, intent(in)     :: idim

    if (idim == 1) field =  0.0_dp
    if (idim == 2) field =  0.0_dp
    if (idim == 3) field = -1.0_dp

  end function gravity_field

  !> analytical fomula for the unit vectors along B
  pure real(dp) function bfield(x, idim) result(field)
    !>$acc routine seq
    real(dp), intent(in)    :: x(1:ndim)
    integer, value, intent(in)     :: idim

    if (idim == 1) field =  0.0_dp
    if (idim == 2) field =  0.0_dp
    if (idim == 3) field = -1.0_dp

  end function bfield

  subroutine addsource_usr(qdt, dtfactor, qtC, wCT, wCTprim, qt, wnew, x, qsourcesplit)
    #:if defined('COOLING')
    use mod_physics, only: rc_fl
    use mod_radiative_cooling, only: getvar_cooling_exact
    #:endif

    double precision, intent(in) :: qdt, dtfactor, qtC, qt
    double precision, intent(in) :: wCT(nw_phys), wCTprim(nw_phys) 
    double precision, intent(in) :: x(ndim)
    double precision, intent(inout) :: wnew(nw_phys)
    logical, intent(in) :: qsourcesplit

    integer             :: idim
    double precision :: coolrate
    double precision :: winit(nw_phys)

    #:if defined('COOLING')
    if (scale /= 0.d0) return
    winit(iw_rho) = 1.d0
    winit(iw_mom(1)) = 0.d0
    winit(iw_e) = 1.d0/phys_gamma/(phys_gamma-1.d0)
    call getvar_cooling_exact(qdt, winit, winit, x, coolrate, rc_fl)
    wnew(iw_e) = wnew(iw_e) + coolrate*qdt
    #:endif
  end subroutine addsource_usr

end module mod_usr
