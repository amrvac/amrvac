
! This example showcases the use of the following subroutines:
!   - usr_print_log
!   - get_volume_average_func

! Example case based on the 2D version of test problem 'radiative_cooling_3D'
! A circular, overdense spot is thermally unstable due to radiative losses.
! In a custom _c.log file, the average kinetic energy and maximum and minimum temperature are written as a time series,
!   at the same rate as the usual .log file. The original .log file is also saved.

module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    unit_length        = 1.d9   ! cm
    unit_temperature   = 1.d6   ! K
    unit_numberdensity = 1.d9   ! cm^-3

    usr_init_one_grid => initonegrid_usr
    usr_source        => special_source
    usr_print_log     => energies_log

    call set_coordinate_system("Cartesian_2D")

    call mhd_activate()

  end subroutine usr_init


  ! Set up a circular density perturbation to excite thermal instability
  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

    use mod_physics
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! set all velocity to zero
    w(ixO^S, mom(:)) = zero
    ! uniform pressure
    w(ixO^S,p_)   = 1.d0
    ! hot central circular spot with uniform pressure
    where((^D&x(ixO^S,^D)**2+) .lt. 0.25d0**2)
      w(ixO^S,rho_) = 2.d0*w(ixO^S,p_)
    elsewhere
      w(ixO^S,rho_) = 1.0d0
    endwhere
    ! magnetic field diagonally upwards
    w(ixO^S,mag(:)) = 1.d0

    call phys_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr


  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

    !use mod_radiative_cooling
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: bQgrid(ixO^S)

    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+(qdt*bQgrid(ixO^S))

  end subroutine special_source


  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)

    ! calculate background heating bQ
    use mod_radiative_cooling
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent (out) :: bQgrid(ixO^S)

    double precision  :: T_cor,l_cor,rho_e

    ! set background density and temperature as in initonegrid_usr
    T_cor = 1.0d6/unit_temperature
    rho_e = 1.0d0

    ! get cooling curve value
    call findL(T_cor,l_cor)

    ! approximate thermal equilibrium: matching the initial losses due to density, except in the perturbation
    bQgrid(ixO^S) = rho_e**2.d0 * l_cor

  end subroutine getbQ


  ! OUTPUT CUSTOM LOG FILE
  !__________________________________________________________

  subroutine energies_log

   use mod_global_parameters
   use mod_input_output, only: printlog_default, get_volume_average_func

   integer, parameter :: my_unit = 123
   character(len=80)  :: fmt_string = '(5(es12.4))', filename
   logical, save      :: alive, visited=.false.
   double precision   :: v_avg, Tmax, Tmin, volume

   ! First output the standard log file:
   call printlog_default

   ! Now make the custom _c.log file:
   ! Calculate average temperature:
   call get_volume_average_func(kinetic, v_avg, volume)
   ! Calculate minimum and maximum temperature
   call get_minmax_temperature(Tmax, Tmin)

   filename = trim(base_filename) // "_c.log"

   if (.not. visited) then
     ! Delete the log when not doing a restart run
     if (restart_from_file == undefined) then
        open(unit=my_unit,file=trim(filename),form='formatted',status='replace')
        write(my_unit,'(a)') '#Global_time mean(v) Tmax Tmin'
     end if
     visited = .true.
  end if

  if (mype == 0) then
     write(filename,"(a)") filename
     inquire(file=filename,exist=alive)
     if(alive) then
       open(unit=my_unit,file=filename,form='formatted',status='old',access='append')
     else
       open(unit=my_unit,file=filename,form='formatted',status='new')
     endif

     ! if number of output doubles is increased, don't forget to change the fmt_string above
      write(my_unit, fmt_string) global_time, v_avg, Tmax, Tmin
      close(my_unit)
    end if

  end subroutine energies_log


  ! Function that calculates kinetic energy, to be used in get_volume_average_func
  pure function kinetic(w_vec, w_size) result(kin_energy)
    integer, intent(in)          :: w_size
    double precision, intent(in) :: w_vec(w_size)
    double precision             :: kin_energy

    kin_energy = 0.5*(w_vec(mom(1))**2 + w_vec(mom(2))**2) / w_vec(rho_) ! divide by rho for conserved
  end function kinetic


  ! Calculates both min and max temperature over all grids in one go by mimicking subroutine get_global_maxima/minima
  subroutine get_minmax_temperature(Tmax, Tmin)

   use mod_global_parameters

   double precision, intent(out) :: Tmax, Tmin

   integer                       :: iigrid, igrid
   double precision              :: Tmax_mype,Tmax_recv,Tmin_mype,Tmin_recv
   double precision              :: w(ixG^T,1:nw),x(ixG^T,1:ndim),wlocal(ixG^T,1:nw),xlocal(ixG^T,1:ndim),Te(ixG^T),pth(ixG^T)

   !!! ixG are indices of grid including ghost cells; ixM is mesh (grid without ghost cells)
   !!! Need to use ^T and ^LL instead of the ^S and ^L you would use for ixO, ixI

   Tmax_mype = -bigdouble
   Tmin_mype = bigdouble

  ! Loop over all the grids and keep track of the temporary min/max
   do iigrid = 1, igridstail
      igrid = igrids(iigrid)

      wlocal(ixG^T,1:nw) = ps(igrid)%w(ixG^T,1:nw)
      xlocal(ixG^T,1:ndim) = ps(igrid)%x(ixG^T,1:ndim)
      ! for some reason, get_pthermal only works with (ixG,ixG) not with (ixM,ixG)
      call mhd_get_pthermal(wlocal,xlocal,ixG^LL,ixG^LL,pth)
      Te(ixG^T) = pth(ixG^T)/wlocal(ixG^T,rho_)

      ! Compare values on current grid to temporary max/min
      Tmax_mype = max(Tmax_mype,maxval(Te(ixM^T)))
      Tmin_mype = min(Tmin_mype,minval(Te(ixM^T)))
   end do

   ! Make the information available on all tasks
   call mpi_allreduce(Tmax_mype, Tmax_recv, 1, mpi_double_precision, &
        mpi_max, icomm, ierrmpi)
   call mpi_allreduce(Tmin_mype, Tmin_recv, 1, mpi_double_precision, &
        mpi_min, icomm, ierrmpi)

   Tmax = Tmax_recv
   Tmin = Tmin_recv

  end subroutine get_minmax_temperature

end module mod_usr
