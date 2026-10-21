!> Analytic TD99 boundary generator and weighted-optimization NLFFF benchmark.
module mod_usr
  use mod_bfield
  use mod_functions_bfield, only: mag
  use mod_nlfff_optimization
  use mod_tdfluxrope
  implicit none

  character(len=32) :: run_mode
  character(len=256) :: boundary_filename,boundary_export_filename
  integer :: fft_padding_factor,nlfff_buffer_cells
  integer :: nlfff_max_iterations,nlfff_log_interval
  integer :: nlfff_plateau_interval,nlfff_plateau_window
  character(len=16) :: fft_top_boundary,lfff_flux_treatment
  character(len=24) :: nlfff_update_preconditioner
  double precision :: i0_effective_ta,lfff_max_flux_imbalance
  double precision :: nlfff_initial_step_scale,nlfff_plateau_tolerance
  logical :: nlfff_write_detailed_history,nlfff_plateau_enabled
  type(nlfff_optimization_result), save :: optimization_result
  logical, save :: firstusrglobaldata=.true.

contains

  subroutine usr_init()
    use mod_geometry, only: set_coordinate_system
    use mod_global_parameters, only: unit_length,unit_temperature,&
         unit_numberdensity
    use mod_usr_methods

    unit_length        = 1.d9
    unit_temperature   = 1.d6
    unit_numberdensity = 1.d10

    call set_coordinate_system('Cartesian_3D')
    usr_set_parameters            => initglobaldata_usr
    usr_init_one_grid             => initonegrid_usr
    usr_improve_initial_condition => improve_initial_condition_usr
    call bfield_activate()
  end subroutine usr_init

  subroutine usr_params_read(files)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    character(len=*), intent(in) :: files(:)
    integer :: n
    namelist /usr_list/ run_mode,boundary_filename,&
         boundary_export_filename,i0_effective_ta,fft_padding_factor,&
         fft_top_boundary,lfff_flux_treatment,lfff_max_flux_imbalance,&
         nlfff_buffer_cells,nlfff_max_iterations,&
         nlfff_initial_step_scale,nlfff_log_interval,&
         nlfff_update_preconditioner,nlfff_write_detailed_history,&
         nlfff_plateau_enabled,nlfff_plateau_interval,&
         nlfff_plateau_window,nlfff_plateau_tolerance

    run_mode='optimization'
    boundary_filename=''
    boundary_export_filename=''
    i0_effective_ta=13.d0
    fft_padding_factor=-1
    fft_top_boundary=''
    lfff_flux_treatment=''
    lfff_max_flux_imbalance=-1.d0
    nlfff_buffer_cells=-1
    nlfff_max_iterations=-1
    nlfff_initial_step_scale=-1.d0
    nlfff_log_interval=-1
    nlfff_update_preconditioner=''
    nlfff_write_detailed_history=.false.
    nlfff_plateau_enabled=.true.
    nlfff_plateau_interval=10
    nlfff_plateau_window=10
    nlfff_plateau_tolerance=1.d-4
    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,usr_list,end=111)
111   close(unitpar)
    end do

    select case(trim(adjustl(run_mode)))
    case('generate_boundary')
      if(len_trim(boundary_export_filename)==0) &
           call mpistop('TD99 boundary_export_filename must be set')
    case('optimization')
      if(len_trim(boundary_filename)==0) &
           call mpistop('TD99 Optimization boundary_filename must be set')
    case default
      call mpistop("run_mode must be 'generate_boundary' or 'optimization'")
    end select
  end subroutine usr_params_read

  subroutine initglobaldata_usr()
    use mod_global_parameters

    call usr_params_read(par_files)
    select case(trim(adjustl(run_mode)))
    case('generate_boundary')
      call configure_td99()
    case('optimization')
      if(firstusrglobaldata) then
        call init_nlfff_optimization_boundary(trim(boundary_filename),&
             unit_length,unit_magneticfield)
        firstusrglobaldata=.false.
      end if
    end select
  end subroutine initglobaldata_usr

  subroutine configure_td99()
    use mod_global_parameters

    double precision :: itube,nt_td99

    Li_TD99    = 0.5d0
    p_Bt_ratio = 1.d0
    d_TD99     = 6.d9/unit_length
    L_TD99     = 6.d9/unit_length
    R_TD99     = 1.d10/unit_length
    a_TD99     = 3.5d9/unit_length
    q_TD99     = 80.d0*1.d20/(unit_magneticfield*unit_length**2)
    Izero_TD99 = (-i0_effective_ta*1.d12)*2.99792456d9/&
         (unit_magneticfield*unit_length*const_c)

    itube=2.d0*q_TD99*L_TD99*R_TD99/&
         ((L_TD99**2+R_TD99**2)**1.5d0*&
          (log(8.d0*R_TD99/a_TD99)-1.5d0+0.5d0*Li_TD99))
    nt_td99=abs(itube/Izero_TD99)*R_TD99**2/a_TD99**2
    if(mype==0) then
      write(*,*) 'TD99 analytic boundary generator'
      write(*,*) 'R,a,L,d [Mm]: ',R_TD99*unit_length/1.d8,&
           a_TD99*unit_length/1.d8,L_TD99*unit_length/1.d8,&
           d_TD99*unit_length/1.d8
      write(*,*) 'q, I0, Nt: ',q_TD99,Izero_TD99,nt_td99
    end if
  end subroutine configure_td99

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,:)=zero
  end subroutine initonegrid_usr

  subroutine improve_initial_condition_usr()
    use mod_global_parameters

    type(nlfff_optimization_config) :: config

    select case(trim(adjustl(run_mode)))
    case('generate_boundary')
      call export_td99_boundary()
    case('optimization')
      config%fft_padding_factor=fft_padding_factor
      config%fft_top_boundary=trim(fft_top_boundary)
      config%flux_treatment=trim(lfff_flux_treatment)
      config%max_flux_imbalance=lfff_max_flux_imbalance
      config%buffer_cells=nlfff_buffer_cells
      config%max_iterations=nlfff_max_iterations
      config%initial_step_scale=nlfff_initial_step_scale
      config%log_interval=nlfff_log_interval
      config%update_preconditioner=trim(nlfff_update_preconditioner)
      config%write_detailed_history=nlfff_write_detailed_history
      config%plateau_enabled=nlfff_plateau_enabled
      config%plateau_interval=nlfff_plateau_interval
      config%plateau_window=nlfff_plateau_window
      config%plateau_tolerance=nlfff_plateau_tolerance
      call extrapolate_nlfff_optimization(mag(:),config,optimization_result)
    end select
  end subroutine improve_initial_condition_usr

  ! The 504x504 map contains two padding cells on each horizontal side.  Its
  ! central 500x500 values match the Optimization domain cell centers.
  subroutine export_td99_boundary()
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    integer :: ixImin1,ixImax1,ixImin2,ixImax2,ixImin3,ixImax3
    integer :: ixOmin1,ixOmax1,ixOmin2,ixOmax2,ixOmin3,ixOmax3
    integer :: ix1,ix2,nxb,nyb,stream_unit
    double precision :: dx_km,dy_km
    double precision, allocatable :: xplane(:,:,:,:),bplane(:,:,:,:)

    if(mype/=0) return
    if(domain_nx1<2*nghostcells+2 .or. domain_nx2<2*nghostcells+2) &
         call mpistop('TD99 V1 exporter mesh is too small for padding')

    nxb=domain_nx1
    nyb=domain_nx2
    allocate(xplane(1:nxb,1:nyb,1:1,1:ndim))
    allocate(bplane(1:nxb,1:nyb,1:1,1:ndir))
    do ix2=1,nyb
      do ix1=1,nxb
        xplane(ix1,ix2,1,1)=xprobmin1+&
             (dble(ix1-nghostcells)-0.5d0)*dx(1,1)
        xplane(ix1,ix2,1,2)=xprobmin2+&
             (dble(ix2-nghostcells)-0.5d0)*dx(2,1)
        xplane(ix1,ix2,1,3)=xprobmin3-0.5d0*dx(3,1)
      end do
    end do

    ixImin1=1; ixImax1=nxb
    ixImin2=1; ixImax2=nyb
    ixImin3=1; ixImax3=1
    ixOmin1=ixImin1; ixOmax1=ixImax1
    ixOmin2=ixImin2; ixOmax2=ixImax2
    ixOmin3=ixImin3; ixOmax3=ixImax3
    call TD99(ixI^L,ixO^L,xplane,bplane)

    dx_km=dx(1,1)*unit_length/1.d5
    dy_km=dx(2,1)*unit_length/1.d5
    open(newunit=stream_unit,file=trim(boundary_export_filename),&
         status='replace',access='stream',form='unformatted',action='write')
    write(stream_unit) zero,nxb,nyb,dx_km,dy_km
    write(stream_unit) bplane(:,:,1,:)*unit_magneticfield
    close(stream_unit)
    write(*,*) 'Exported TD99 V1 boundary: ',&
         trim(boundary_export_filename),nxb,nyb
    deallocate(xplane,bplane)
  end subroutine export_td99_boundary

end module mod_usr
