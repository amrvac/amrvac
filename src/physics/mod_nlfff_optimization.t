!> Nonlinear force-free extrapolation by the weighted optimization method.
!>
!> The implementation follows Equations (10), (13), (15), and (17) of
!> Wiegelmann (2004), Solar Physics 219, 87, using fixed nodal physical
!> boundaries and explicit adaptive-step controls.
module mod_nlfff_optimization
  use mod_physicaldata, only: state
  implicit none
  private

  double precision, parameter :: nlfff_step_grow=1.01d0
  double precision, parameter :: nlfff_step_shrink=0.5d0
  ! Thomas's reference controller uses mue > 1e-7 * dx**2 as its
  ! absolute pseudo-time-step floor.  Keep this tied to the mesh spacing,
  ! rather than to a user-selected initial step scale.
  double precision, parameter :: nlfff_step_floor_dx2=1.d-7

  type, public :: nlfff_optimization_config
    integer :: fft_padding_factor=-1
    character(len=16) :: fft_top_boundary=''
    character(len=16) :: flux_treatment=''
    double precision :: max_flux_imbalance=-1.d0
    integer :: buffer_cells=-1
    integer :: max_iterations=-1
    double precision :: initial_step_scale=-1.d0
    integer :: log_interval=-1
    character(len=24) :: update_preconditioner=''
    character(len=16) :: initialization_mode='potential'
    logical :: write_detailed_history=.false.
    logical :: plateau_enabled=.true.
    integer :: plateau_interval=10
    integer :: plateau_window=10
    double precision :: plateau_tolerance=1.d-4
  end type nlfff_optimization_config

  type, public :: nlfff_optimization_result
    integer :: attempts=0
    integer :: accepted_steps=0
    integer :: rejected_steps=0
    integer :: functional_evaluations=0
    character(len=32) :: stop_reason='not_started'
    double precision :: initial_L=0.d0
    double precision :: initial_L_force=0.d0
    double precision :: initial_L_div=0.d0
    double precision :: initial_B2_integral=0.d0
    double precision :: initial_epsilon_force=0.d0
    double precision :: initial_epsilon_div=0.d0
    double precision :: step_floor=0.d0
    double precision :: final_L=0.d0
    double precision :: final_L_force=0.d0
    double precision :: final_L_div=0.d0
    double precision :: final_B2_integral=0.d0
    double precision :: final_epsilon_force=0.d0
    double precision :: final_epsilon_div=0.d0
    double precision :: final_step=0.d0
    integer :: plateau_count=0
    double precision :: final_relative_functional_change=1.d0
  end type nlfff_optimization_result

{^IFTHREED
  ! Internal performance counters.  They are reported as rank maxima after a
  ! run and are deliberately kept out of the public result and CSV schemas.
  type :: nlfff_optimization_timing
    double precision :: update_kernel=0.d0
    double precision :: binomial_halo_x=0.d0
    double precision :: binomial_halo_y=0.d0
    double precision :: trial_b_exchange=0.d0
    double precision :: qs_exchange=0.d0
    double precision :: functional_local=0.d0
    double precision :: functional_reductions=0.d0
    double precision :: unified_diagnostics=0.d0
  end type nlfff_optimization_timing

  type(nlfff_optimization_timing), save :: optimization_timing
  logical, save :: optimization_timing_active=.false.

  double precision, allocatable, save :: bottom_b(:,:,:)
  ! Static Cartesian-grid weights are reused by every trial and functional
  ! evaluation.  Cache one-dimensional factors rather than a block-sized
  ! three-dimensional array, so the cache remains negligible on HMI cases.
  double precision, allocatable, save :: nlfff_weight_x(:)
  double precision, allocatable, save :: nlfff_weight_y(:)
  double precision, allocatable, save :: nlfff_weight_z(:)
  ! Q=(Omega_a x B) and s=(Omega_b dot B) require four communication slots,
  ! independent of the persistent physics width.  Keeping accepted and trial
  ! values separate lets rejected trials be discarded without rollback.
  type(state), allocatable, target, save :: accepted_qs_state(:)
  type(state), allocatable, target, save :: trial_qs_state(:)

  public :: init_nlfff_optimization_boundary
  public :: extrapolate_nlfff_optimization
  public :: nlfff_cosine_side_weight
  public :: nlfff_cosine_top_weight
  public :: nlfff_weight_product
  public :: nlfff_fixed_active_index
  public :: nlfff_boundary_derivative
  public :: nlfff_compute_omega_cell
  public :: nlfff_compose_update_cell
  public :: nlfff_dimensionless_metrics
}

contains

{^IFTHREED
  !> Read one Python V1 data-driven vector magnetogram.  The same read also
  !> initializes mod_lfff's normalized Bz data for the potential FFT.
  subroutine init_nlfff_optimization_boundary(filename,qLunit,qBunit,qxc1,qxc2)
    use mod_lfff, only: init_b_fff_data_driven_boundary

    character(len=*), intent(in) :: filename
    double precision, intent(in) :: qLunit,qBunit
    double precision, intent(in), optional :: qxc1,qxc2

    if(allocated(bottom_b)) deallocate(bottom_b)
    call init_b_fff_data_driven_boundary(filename,qLunit,qBunit,qxc1,qxc2,bottom_b)
  end subroutine init_nlfff_optimization_boundary

  !> Perform a one-shot, fixed-grid weighted optimization extrapolation.
  subroutine extrapolate_nlfff_optimization(iw_b,config,result)
    use mpi
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_lfff, only: extrapolate_potential_fft
    use mod_nlfff_diagnostics, only: nlfff_physical_metrics,&
       evaluate_nlfff_metrics_amrvac,write_nlfff_metrics_header,&
       write_nlfff_metrics_row

    integer, intent(in) :: iw_b(3)
    type(nlfff_optimization_config), intent(in) :: config
    type(nlfff_optimization_result), intent(out) :: result

    double precision :: step,step_floor,Lold,Lforce_old,Ldiv_old,B2old
    double precision :: Ltrial,Lforce_trial,Ldiv_trial,B2trial
    double precision :: epsilon_force,epsilon_div,characteristic_spacing
    double precision :: plateau_reference_L,plateau_relative_change
    double precision :: timing_start
    type(nlfff_physical_metrics) :: physical_metrics
    integer :: log_unit,metrics_unit
    integer :: plateau_count
    logical :: accepted_trial,plateau_stop

    result=nlfff_optimization_result()
    optimization_timing=nlfff_optimization_timing()
    optimization_timing_active=.true.
    call validate_configuration(iw_b,config)
    call initialize_nlfff_weight_cache(config%buffer_cells)
    call validate_exchange_state()
    call allocate_qs_states()
    call prepare_partial_exchange_types(iw_b)

    if(trim(config%initialization_mode)=='potential') then
      ! The observed magnetogram is centered in the adjacent lower ghost layer.
      call extrapolate_potential_fft(iw_b,config%fft_padding_factor,&
         0.5d0*dx(3,1),0.d0,config%fft_top_boundary,&
         config%flux_treatment,config%max_flux_imbalance)
    end if

    call apply_physical_boundaries(ps,iw_b)
    call exchange_b_internal(iw_b)
    call evaluate_functional(ps,iw_b,accepted_qs_state,config%buffer_cells,&
       Lold,Lforce_old,Ldiv_old,B2old)
    result%functional_evaluations=1

    step=config%initial_step_scale*minval(dx(:,1))**2
    step_floor=nlfff_step_floor_dx2*minval(dx(:,1))**2
    result%step_floor=step_floor
    plateau_count=0
    plateau_reference_L=Lold
    plateau_relative_change=1.d0
    plateau_stop=.false.
    result%final_relative_functional_change=plateau_relative_change
    characteristic_spacing=product(dx(:,1))**(1.d0/3.d0)
    call nlfff_dimensionless_metrics(Lforce_old,Ldiv_old,B2old,&
       characteristic_spacing,epsilon_force,epsilon_div)
    result%initial_L=Lold
    result%initial_L_force=Lforce_old
    result%initial_L_div=Ldiv_old
    result%initial_B2_integral=B2old
    result%initial_epsilon_force=epsilon_force
    result%initial_epsilon_div=epsilon_div
    result%final_L=Lold
    result%final_L_force=Lforce_old
    result%final_L_div=Ldiv_old
    result%final_B2_integral=B2old
    result%final_epsilon_force=epsilon_force
    result%final_epsilon_div=epsilon_div
    result%final_step=step

    log_unit=-1
    metrics_unit=-1
    accepted_trial=.false.
    if(mype==0) then
      if(config%write_detailed_history) then
        open(newunit=log_unit,file=trim(base_filename)//'_nlfff_opt.csv',&
           status='replace',action='write')
        write(log_unit,'(a)') 'attempt,accepted,rejected,step,L,L_force,L_div,'//&
           'B2_integral,epsilon_force,epsilon_div,functional_evaluations,'//&
           'step_floor,plateau_count,plateau_relative_change,accepted_trial'
        call write_log_row(log_unit,result,step,Lold,Lforce_old,Ldiv_old,B2old,&
           epsilon_force,epsilon_div,accepted_trial)
      end if
      open(newunit=metrics_unit,file=trim(base_filename)//'_nlfff_metrics.csv',&
         status='replace',action='write')
      call write_nlfff_metrics_header(metrics_unit)
      write(*,*) 'Weighted NLFFF optimization initial L:',Lold
    end if
    timing_start=MPI_WTIME()
    call evaluate_nlfff_metrics_amrvac(iw_b,physical_metrics)
    call accumulate_timing(optimization_timing%unified_diagnostics,timing_start)
    if(mype==0) call write_nlfff_metrics_row(metrics_unit,0,physical_metrics)

    do while(result%accepted_steps<config%max_iterations .and. step>=step_floor .and. &
       .not.plateau_stop)
      result%attempts=result%attempts+1
      call form_trial_field(iw_b,config%buffer_cells,&
         config%update_preconditioner,step)
      call apply_physical_boundaries(ps1,iw_b)
      timing_start=MPI_WTIME()
      call exchange_trial_b_internal(iw_b)
      call accumulate_timing(optimization_timing%trial_b_exchange,timing_start)
      call evaluate_functional(ps1,iw_b,trial_qs_state,config%buffer_cells,&
         Ltrial,Lforce_trial,Ldiv_trial,B2trial)
      result%functional_evaluations=result%functional_evaluations+1

      accepted_trial=(Ltrial<Lold)
      if(accepted_trial) then
        result%accepted_steps=result%accepted_steps+1
        call accept_trial_field(iw_b)
        call swap_qs_states()
        Lold=Ltrial
        Lforce_old=Lforce_trial
        Ldiv_old=Ldiv_trial
        B2old=B2trial
        step=step*nlfff_step_grow
      else
        result%rejected_steps=result%rejected_steps+1
        step=step*nlfff_step_shrink
      end if

      result%final_L=Lold
      result%final_L_force=Lforce_old
      result%final_L_div=Ldiv_old
      result%final_B2_integral=B2old
      call nlfff_dimensionless_metrics(Lforce_old,Ldiv_old,B2old,&
         characteristic_spacing,epsilon_force,epsilon_div)
      result%final_epsilon_force=epsilon_force
      result%final_epsilon_div=epsilon_div
      result%final_step=step

      ! Thomas checks the relative functional change every ten accepted
      ! iterations and stops after ten consecutive plateau checks.  Rejected
      ! trials do not advance this diagnostic, matching Thomas's it=it-1
      ! rollback of rejected iterations.
      if(config%plateau_enabled .and. accepted_trial .and. &
         result%accepted_steps>=config%plateau_interval .and. &
         mod(result%accepted_steps,config%plateau_interval)==0) then
        plateau_relative_change=abs((Lold-plateau_reference_L)/&
           max(abs(Lold),tiny(1.d0)))
        if(plateau_relative_change<config%plateau_tolerance) then
          plateau_count=plateau_count+1
        else
          plateau_count=0
        end if
        plateau_reference_L=Lold
        result%plateau_count=plateau_count
        result%final_relative_functional_change=plateau_relative_change
        if(plateau_count>=config%plateau_window) plateau_stop=.true.
      end if
      if(mod(result%attempts,config%log_interval)==0) then
        timing_start=MPI_WTIME()
        call evaluate_nlfff_metrics_amrvac(iw_b,physical_metrics)
        call accumulate_timing(optimization_timing%unified_diagnostics,timing_start)
        if(mype==0) then
          if(config%write_detailed_history) call write_log_row(log_unit,result,&
             step,Lold,Lforce_old,Ldiv_old,B2old,epsilon_force,epsilon_div,&
             accepted_trial)
          call write_nlfff_metrics_row(metrics_unit,result%attempts,&
             physical_metrics)
          write(*,*) 'NLFFF optimization:',result%accepted_steps,&
             result%rejected_steps,step,Lold
        end if
      end if
    end do

    if(plateau_stop) then
      result%stop_reason='plateau'
    else if(result%accepted_steps>=config%max_iterations) then
      result%stop_reason='max_iterations'
    else
      result%stop_reason='step_floor'
    end if
    result%final_step=step
    if(mod(result%attempts,config%log_interval)/=0) then
      timing_start=MPI_WTIME()
      call evaluate_nlfff_metrics_amrvac(iw_b,physical_metrics)
      call accumulate_timing(optimization_timing%unified_diagnostics,timing_start)
      if(mype==0) then
        if(config%write_detailed_history) call write_log_row(log_unit,result,&
           step,Lold,Lforce_old,Ldiv_old,B2old,epsilon_force,epsilon_div,&
           accepted_trial)
        call write_nlfff_metrics_row(metrics_unit,result%attempts,&
           physical_metrics)
      end if
    end if
    if(mype==0) then
      if(config%write_detailed_history) close(log_unit)
      close(metrics_unit)
      write(*,*) 'Weighted NLFFF optimization stopped:',trim(result%stop_reason)
      write(*,*) 'accepted, rejected, final L:',result%accepted_steps,&
         result%rejected_steps,result%final_L
      write(*,*) 'step floor, plateau count, relative dL/L:',&
         result%step_floor,result%plateau_count,&
         result%final_relative_functional_change
      write(*,*) 'functional evaluations:',result%functional_evaluations
    end if
    call report_nlfff_timing()
    optimization_timing_active=.false.
    call deallocate_nlfff_weight_cache()
    call restore_full_exchange_types()
    call deallocate_qs_states()
  end subroutine extrapolate_nlfff_optimization

  subroutine allocate_qs_states()
    use mod_amr_solution_node, only: alloc_state
    use mod_global_parameters

    integer :: iigrid,igrid

    if(allocated(accepted_qs_state)) deallocate(accepted_qs_state)
    if(allocated(trial_qs_state)) deallocate(trial_qs_state)
    allocate(accepted_qs_state(max_blocks),trial_qs_state(max_blocks))
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      call alloc_state(igrid,accepted_qs_state(igrid),ixG^LL,ixG^LL,.false.)
      call alloc_state(igrid,trial_qs_state(igrid),ixG^LL,ixG^LL,.false.)
      deallocate(accepted_qs_state(igrid)%w,trial_qs_state(igrid)%w)
      allocate(accepted_qs_state(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:4))
      allocate(trial_qs_state(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:4))
      accepted_qs_state(igrid)%w=0.d0
      trial_qs_state(igrid)%w=0.d0
    end do
  end subroutine allocate_qs_states

  subroutine deallocate_qs_states()
    if(allocated(accepted_qs_state)) deallocate(accepted_qs_state)
    if(allocated(trial_qs_state)) deallocate(trial_qs_state)
  end subroutine deallocate_qs_states

  subroutine initialize_nlfff_weight_cache(buffer_cells)
    use mod_global_parameters

    integer, intent(in) :: buffer_cells
    integer :: ig1,ig2,ig3

    if(allocated(nlfff_weight_x)) deallocate(nlfff_weight_x)
    if(allocated(nlfff_weight_y)) deallocate(nlfff_weight_y)
    if(allocated(nlfff_weight_z)) deallocate(nlfff_weight_z)
    allocate(nlfff_weight_x(1-nghostcells:domain_nx1+nghostcells))
    allocate(nlfff_weight_y(1-nghostcells:domain_nx2+nghostcells))
    allocate(nlfff_weight_z(1-nghostcells:domain_nx3+nghostcells))
    do ig1=lbound(nlfff_weight_x,1),ubound(nlfff_weight_x,1)
      nlfff_weight_x(ig1)=nlfff_cosine_side_weight(ig1,domain_nx1,buffer_cells)
    end do
    do ig2=lbound(nlfff_weight_y,1),ubound(nlfff_weight_y,1)
      nlfff_weight_y(ig2)=nlfff_cosine_side_weight(ig2,domain_nx2,buffer_cells)
    end do
    do ig3=lbound(nlfff_weight_z,1),ubound(nlfff_weight_z,1)
      nlfff_weight_z(ig3)=nlfff_cosine_top_weight(ig3,domain_nx3,buffer_cells)
    end do
  end subroutine initialize_nlfff_weight_cache

  subroutine deallocate_nlfff_weight_cache()
    if(allocated(nlfff_weight_x)) deallocate(nlfff_weight_x)
    if(allocated(nlfff_weight_y)) deallocate(nlfff_weight_y)
    if(allocated(nlfff_weight_z)) deallocate(nlfff_weight_z)
  end subroutine deallocate_nlfff_weight_cache

  logical function cached_fixed_active_cell(igrid,ix1,ix2,ix3)
    use mod_global_parameters, only: block_nx1,block_nx2,block_nx3,nghostcells,node,&
       pig1_,pig2_,pig3_,domain_nx1,domain_nx2,domain_nx3

    integer, intent(in) :: igrid,ix1,ix2,ix3
    integer :: ig1,ig2,ig3

    ig1=(node(pig1_,igrid)-1)*block_nx1+ix1-nghostcells
    ig2=(node(pig2_,igrid)-1)*block_nx2+ix2-nghostcells
    ig3=(node(pig3_,igrid)-1)*block_nx3+ix3-nghostcells
    cached_fixed_active_cell=nlfff_fixed_active_index(ig1,ig2,ig3,&
       domain_nx1,domain_nx2,domain_nx3)
  end function cached_fixed_active_cell

  double precision function cached_cell_weight(igrid,ix1,ix2,ix3)
    use mod_global_parameters, only: block_nx1,block_nx2,block_nx3,nghostcells,node,&
       pig1_,pig2_,pig3_

    integer, intent(in) :: igrid,ix1,ix2,ix3
    integer :: ig1,ig2,ig3

    ig1=(node(pig1_,igrid)-1)*block_nx1+ix1-nghostcells
    ig2=(node(pig2_,igrid)-1)*block_nx2+ix2-nghostcells
    ig3=(node(pig3_,igrid)-1)*block_nx3+ix3-nghostcells
    cached_cell_weight=nlfff_weight_x(ig1)*nlfff_weight_y(ig2)*nlfff_weight_z(ig3)
  end function cached_cell_weight

  subroutine report_nlfff_timing()
    use mpi
    use mod_global_parameters, only: icomm,ierrmpi,mype

    double precision :: local_values(8),global_values(8)

    local_values=(/optimization_timing%update_kernel,&
       optimization_timing%binomial_halo_x,&
       optimization_timing%binomial_halo_y,&
       optimization_timing%trial_b_exchange,&
       optimization_timing%qs_exchange,&
       optimization_timing%functional_local,&
       optimization_timing%functional_reductions,&
       optimization_timing%unified_diagnostics/)
    call MPI_ALLREDUCE(local_values,global_values,8,MPI_DOUBLE_PRECISION,&
       MPI_MAX,icomm,ierrmpi)
    if(mype==0) then
      write(*,'(a,8(1x,es16.8))') 'NLFFF timing rank-max [s]:',global_values
      write(*,'(a)') '  update_kernel binomial_halo_x binomial_halo_y '//&
         'trial_b_exchange qs_exchange functional_local '//&
         'functional_reductions unified_diagnostics'
    end if
  end subroutine report_nlfff_timing

  subroutine accumulate_timing(counter,timing_start)
    use mpi

    double precision, intent(inout) :: counter
    double precision, intent(in) :: timing_start

    if(optimization_timing_active) counter=counter+MPI_WTIME()-timing_start
  end subroutine accumulate_timing

  subroutine validate_configuration(iw_b,config)
    use mod_comm_lib, only: mpistop
    use mod_geometry, only: Cartesian,coordinate
    use mod_global_parameters
    use mod_lfff, only: xa1,xa2
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    integer, intent(in) :: iw_b(3)
    type(nlfff_optimization_config), intent(in) :: config
    double precision :: dx1m,dx2m,tol
    integer :: pad1,pad2

    if(ndim/=3 .or. ndir/=3) call mpistop('NLFFF optimization requires 3D vectors')
    if(coordinate/=Cartesian) call mpistop('NLFFF optimization requires Cartesian coordinates')
    if(any(stretched_dim)) call mpistop('NLFFF optimization requires a uniform mesh')
    if(refine_max_level/=1 .or. levmax/=1) &
       call mpistop('NLFFF optimization requires refine_max_level=1')
    if(stagger_grid) call mpistop('NLFFF optimization requires cell-centered B')
    if(B0field) call mpistop('NLFFF optimization does not support B0 splitting')
    if(any(periodB)) call mpistop('NLFFF optimization requires physical side and top boundaries')
    if(nghostcells<2) call mpistop('NLFFF optimization requires at least two ghost cells')
    if(block_nx1<3 .or. block_nx2<3 .or. block_nx3<3) &
       call mpistop('NLFFF one-sided boundary stencil requires at least three cells per block')
    if(nw<3) call mpistop('NLFFF optimization requires three magnetic variables')
    if(iw_b(1)<1 .or. iw_b(3)>nw) &
       call mpistop('NLFFF magnetic variable indices are outside the state')
    if(any(iw_b/=(/iw_b(1),iw_b(1)+1,iw_b(1)+2/))) &
       call mpistop('NLFFF optimization requires contiguous B components')
    if(.not.allocated(bottom_b)) call mpistop('NLFFF vector boundary has not been initialized')
    if(.not.all(ieee_is_finite(bottom_b))) call mpistop('NLFFF vector boundary is non-finite')

    if(config%fft_padding_factor<1) call mpistop('NLFFF fft_padding_factor must be set')
    if(trim(config%fft_top_boundary)/='open' .and. &
       trim(config%fft_top_boundary)/='closed') &
       call mpistop("NLFFF fft_top_boundary must be 'open' or 'closed'")
    if(trim(config%flux_treatment)/='strict' .and. &
       trim(config%flux_treatment)/='subtract_mean') &
       call mpistop("NLFFF flux_treatment must be 'strict' or 'subtract_mean'")
    if(.not.ieee_is_finite(config%max_flux_imbalance) .or. &
       config%max_flux_imbalance<0.d0 .or. config%max_flux_imbalance>1.d0) &
       call mpistop('NLFFF max_flux_imbalance must be explicitly set in [0,1]')
    if(config%buffer_cells<2) call mpistop('NLFFF buffer_cells must be at least two')
    if(2*config%buffer_cells>=domain_nx1 .or. &
       2*config%buffer_cells>=domain_nx2 .or. &
       config%buffer_cells>=domain_nx3) &
       call mpistop('NLFFF cosine buffer leaves no interior physical region')
    if(config%max_iterations<1) call mpistop('NLFFF max_iterations must be positive')
    if(.not.ieee_is_finite(config%initial_step_scale) .or. &
       config%initial_step_scale<=0.d0) &
       call mpistop('NLFFF initial_step_scale must be positive')
    if(config%log_interval<1) call mpistop('NLFFF log_interval must be positive')
    if(config%plateau_interval<1) &
       call mpistop('NLFFF plateau_interval must be positive')
    if(config%plateau_window<1) &
       call mpistop('NLFFF plateau_window must be positive')
    if(.not.ieee_is_finite(config%plateau_tolerance) .or. &
       config%plateau_tolerance<0.d0) &
       call mpistop('NLFFF plateau_tolerance must be non-negative')
    if(trim(config%update_preconditioner)/='none' .and. &
       trim(config%update_preconditioner)/='binomial_xy') &
       call mpistop("NLFFF update_preconditioner must be 'none' or 'binomial_xy'")
    if(trim(config%initialization_mode)/='potential' .and. &
       trim(config%initialization_mode)/='current_state') &
       call mpistop("NLFFF initialization_mode must be 'potential' or 'current_state'")

    if(size(bottom_b,1)<domain_nx1 .or. size(bottom_b,2)<domain_nx2) &
       call mpistop('NLFFF vector boundary is smaller than the physical domain')
    if(mod(size(bottom_b,1)-domain_nx1,2)/=0 .or. &
       mod(size(bottom_b,2)-domain_nx2,2)/=0) &
       call mpistop('NLFFF vector boundary padding must be symmetric')
    pad1=(size(bottom_b,1)-domain_nx1)/2
    pad2=(size(bottom_b,2)-domain_nx2)/2
    if(pad1<nghostcells .or. pad2<nghostcells) &
       call mpistop('NLFFF vector boundary lacks horizontal ghost padding')
    if(size(xa1)<2 .or. size(xa2)<2) call mpistop('NLFFF boundary coordinates are incomplete')
    dx1m=xa1(2)-xa1(1)
    dx2m=xa2(2)-xa2(1)
    tol=1.d-10*max(1.d0,abs(dx(1,1)),abs(dx(2,1)))
    if(abs(dx1m-dx(1,1))>tol .or. abs(dx2m-dx(2,1))>tol) &
       call mpistop('NLFFF vector boundary spacing does not match the grid')
    if(abs(xa1(pad1+1)-(xprobmin1+0.5d0*dx(1,1)))>tol .or. &
       abs(xa1(pad1+domain_nx1)-(xprobmax1-0.5d0*dx(1,1)))>tol .or. &
       abs(xa2(pad2+1)-(xprobmin2+0.5d0*dx(2,1)))>tol .or. &
       abs(xa2(pad2+domain_nx2)-(xprobmax2-0.5d0*dx(2,1)))>tol) &
       call mpistop('NLFFF vector boundary center does not match the grid')
  end subroutine validate_configuration

  !> The initialization hook normally enters with the full-state datatype.
  !> Refuse a nested partial exchange because its MPI datatype targets would
  !> otherwise be overwritten and could not be reconstructed here.
  subroutine validate_exchange_state()
    use mod_comm_lib, only: mpistop
    use mod_ghostcells_update

    if(.not.associated(type_send_srl,type_send_srl_f) .or. &
       .not.associated(type_recv_srl,type_recv_srl_f) .or. &
       .not.associated(type_send_r,type_send_r_f) .or. &
       .not.associated(type_recv_r,type_recv_r_f) .or. &
       .not.associated(type_send_p,type_send_p_f) .or. &
       .not.associated(type_recv_p,type_recv_p_f) .or. .not.bcphys) &
       call mpistop('NLFFF optimization requires full-state ghost exchange on entry')
  end subroutine validate_exchange_state

  subroutine prepare_partial_exchange_types(iw_b)
    use mod_ghostcells_update

    integer, intent(in) :: iw_b(3)

    type_send_srl=>type_send_srl_p1
    type_recv_srl=>type_recv_srl_p1
    type_send_r=>type_send_r_p1
    type_recv_r=>type_recv_r_p1
    type_send_p=>type_send_p_p1
    type_recv_p=>type_recv_p_p1
    call create_bc_mpi_datatype(iw_b(1),3)

    type_send_srl=>type_send_srl_p2
    type_recv_srl=>type_recv_srl_p2
    type_send_r=>type_send_r_p2
    type_recv_r=>type_recv_r_p2
    type_send_p=>type_send_p_p2
    type_recv_p=>type_recv_p_p2
    call create_bc_mpi_datatype(1,4,4)
    call restore_full_exchange_types()
  end subroutine prepare_partial_exchange_types

  subroutine restore_full_exchange_types()
    use mod_ghostcells_update

    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f
    bcphys=.true.
  end subroutine restore_full_exchange_types

  subroutine exchange_b_internal(iw_b)
    use mod_global_parameters
    use mod_ghostcells_update

    integer, intent(in) :: iw_b(3)

    type_send_srl=>type_send_srl_p1
    type_recv_srl=>type_recv_srl_p1
    type_send_r=>type_send_r_p1
    type_recv_r=>type_recv_r_p1
    type_send_p=>type_send_p_p1
    type_recv_p=>type_recv_p_p1
    bcphys=.false.
    call getbc(global_time,0.d0,ps,iw_b(1),3)
    call restore_full_exchange_types()
  end subroutine exchange_b_internal

  subroutine exchange_qs_internal(qs_state)
    use mod_global_parameters
    use mod_ghostcells_update

    type(state), target, intent(inout) :: qs_state(max_blocks)

    type_send_srl=>type_send_srl_p2
    type_recv_srl=>type_recv_srl_p2
    type_send_r=>type_send_r_p2
    type_recv_r=>type_recv_r_p2
    type_send_p=>type_send_p_p2
    type_recv_p=>type_recv_p_p2
    bcphys=.false.
    call getbc(global_time,0.d0,qs_state,1,4)
    call restore_full_exchange_types()
  end subroutine exchange_qs_internal

  subroutine exchange_trial_b_internal(iw_b)
    use mod_global_parameters
    use mod_ghostcells_update

    integer, intent(in) :: iw_b(3)

    type_send_srl=>type_send_srl_p1
    type_recv_srl=>type_recv_srl_p1
    type_send_r=>type_send_r_p1
    type_recv_r=>type_recv_r_p1
    type_send_p=>type_send_p_p1
    type_recv_p=>type_recv_p_p1
    bcphys=.false.
    call getbc(global_time,0.d0,ps1,iw_b(1),3)
    call restore_full_exchange_types()
  end subroutine exchange_trial_b_internal

  !> Fill physical ghosts without invoking user MHD boundary callbacks.
  subroutine apply_physical_boundaries(field_state,iw_b)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    use mod_lfff, only: xa1,xa2

    integer, intent(in) :: iw_b(3)
    type(state), target, intent(inout) :: field_state(max_blocks)
    integer :: iigrid,igrid,idir,ix1,ix2,ix3,ib1,ib2
    double precision :: dx1m,dx2m

    dx1m=xa1(2)-xa1(1)
    dx2m=xa2(2)-xa2(1)
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      if(ps(igrid)%is_physical_boundary(1)) then
        do ix1=ixGlo1,ixMlo1-1
          do idir=1,3
            field_state(igrid)%w(ix1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_b(idir))=&
               field_state(igrid)%w(ixMlo1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_b(idir))
          end do
        end do
      end if
      if(ps(igrid)%is_physical_boundary(2)) then
        do ix1=ixMhi1+1,ixGhi1
          do idir=1,3
            field_state(igrid)%w(ix1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_b(idir))=&
               field_state(igrid)%w(ixMhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_b(idir))
          end do
        end do
      end if
      if(ps(igrid)%is_physical_boundary(3)) then
        do ix2=ixGlo2,ixMlo2-1
          do idir=1,3
            field_state(igrid)%w(ixGlo1:ixGhi1,ix2,ixGlo3:ixGhi3,iw_b(idir))=&
               field_state(igrid)%w(ixGlo1:ixGhi1,ixMlo2,ixGlo3:ixGhi3,iw_b(idir))
          end do
        end do
      end if
      if(ps(igrid)%is_physical_boundary(4)) then
        do ix2=ixMhi2+1,ixGhi2
          do idir=1,3
            field_state(igrid)%w(ixGlo1:ixGhi1,ix2,ixGlo3:ixGhi3,iw_b(idir))=&
               field_state(igrid)%w(ixGlo1:ixGhi1,ixMhi2,ixGlo3:ixGhi3,iw_b(idir))
          end do
        end do
      end if
      if(ps(igrid)%is_physical_boundary(6)) then
        do ix3=ixMhi3+1,ixGhi3
          do idir=1,3
            field_state(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix3,iw_b(idir))=&
               field_state(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixMhi3,iw_b(idir))
          end do
        end do
      end if

      ! Apply the measured lower boundary last, so it owns lower-side corners.
      if(ps(igrid)%is_physical_boundary(5)) then
        do ix3=ixGlo3,ixMlo3-1
          do ix2=ixGlo2,ixGhi2
            do ix1=ixGlo1,ixGhi1
              ib1=nint((ps(igrid)%x(ix1,ix2,ix3,1)-xa1(1))/dx1m)+1
              ib2=nint((ps(igrid)%x(ix1,ix2,ix3,2)-xa2(1))/dx2m)+1
              if(ib1<1 .or. ib1>size(bottom_b,1) .or. &
                 ib2<1 .or. ib2>size(bottom_b,2)) then
                write(*,*) 'NLFFF boundary pixel outside frame:',mype,igrid,ib1,ib2
                call mpistop('NLFFF lower-boundary mapping failed')
              end if
              do idir=1,3
                field_state(igrid)%w(ix1,ix2,ix3,iw_b(idir))=&
                   bottom_b(ib1,ib2,idir)
              end do
            end do
          end do
        end do
      end if
    end do
  end subroutine apply_physical_boundaries

  !> Evaluate the unmodified weighted functional and prepare Q,s in one pass.
  subroutine evaluate_functional(field_state,iw_b,qs_state,buffer_cells,L,&
     Lforce,Ldiv,B2integral)
    use mpi
    use mod_global_parameters
    use mod_geometry, only: curlvector,divvector
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    integer, intent(in) :: iw_b(3),buffer_cells
    type(state), target, intent(in) :: field_state(max_blocks)
    type(state), target, intent(inout) :: qs_state(max_blocks)
    double precision, intent(out) :: L,Lforce,Ldiv,B2integral

    double precision :: local_force,local_div,local_B2
    double precision :: part_force,part_div,part_B2
    double precision :: cell_force,cell_div,diagnostic(3)
    double precision :: local_values(3),global_values(3)
    double precision :: timing_start
    double precision :: bvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:3)
    double precision :: current(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:3)
    double precision :: divb_array(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)
    double precision :: oa(3),ob(3),b(3),j(3),divb,b2,wcell,cell_volume
    integer :: iigrid,igrid,ix1,ix2,ix3,ixB3,idirmin

    local_force=0.d0
    local_div=0.d0
    local_B2=0.d0
    timing_start=MPI_WTIME()
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      qs_state(igrid)%w(ixG^T,1:4)=0.d0
    end do

    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      block=>ps(igrid)
      dxlevel(:)=dx(:,1)
      do ix3=ixGlo3,ixGhi3
        do ix2=ixGlo2,ixGhi2
          do ix1=ixGlo1,ixGhi1
            bvec(ix1,ix2,ix3,:)=&
               field_state(igrid)%w(ix1,ix2,ix3,iw_b(:))
          end do
        end do
      end do
      current=0.d0
      divb_array=0.d0
      idirmin=1
      call curlvector(bvec,ixG^LL,ixM^LL,current,idirmin,1,3)
      call divvector(bvec,ixG^LL,ixM^LL,divb_array,1)
      call override_fixed_boundary_curl_div(bvec,igrid,current,divb_array)
      part_force=0.d0
      part_div=0.d0
      part_B2=0.d0
      do ix3=ixMlo3,ixMhi3
        do ix2=ixMlo2,ixMhi2
          do ix1=ixMlo1,ixMhi1
            b=bvec(ix1,ix2,ix3,:)
            j=current(ix1,ix2,ix3,:)
            divb=divb_array(ix1,ix2,ix3)
            b2=dot_product(b,b)
            call nlfff_compute_omega_cell(b,j,divb,oa,ob)
            qs_state(igrid)%w(ix1,ix2,ix3,1:3)=cross3(oa,b)
            qs_state(igrid)%w(ix1,ix2,ix3,4)=dot_product(ob,b)
            wcell=cached_cell_weight(igrid,ix1,ix2,ix3)
            cell_volume=ps(igrid)%dvolume(ix1,ix2,ix3)
            cell_force=0.d0
            cell_div=0.d0
            if(b2>0.d0) then
              cell_force=wcell*b2*dot_product(oa,oa)*cell_volume
              cell_div=wcell*b2*dot_product(ob,ob)*cell_volume
            end if
            if(.not.all(ieee_is_finite(b)) .or. &
               .not.all(ieee_is_finite(j)) .or. &
               .not.ieee_is_finite(divb) .or. &
               .not.all(ieee_is_finite(oa)) .or. &
               .not.all(ieee_is_finite(ob)) .or. &
               .not.ieee_is_finite(cell_force) .or. &
               .not.ieee_is_finite(cell_div)) then
              diagnostic=(/wcell,cell_force,cell_div/)
              call report_nonfinite('functional',igrid,ix1,ix2,ix3,&
                 b,j,divb,oa,ob,diagnostic)
            end if
            part_force=part_force+cell_force
            part_div=part_div+cell_div
            part_B2=part_B2+wcell*b2*cell_volume
          end do
        end do
      end do

      ! Match the nodal optimization discretization at the photosphere: the
      ! fixed magnetogram node contributes one full uniform-cell volume to L.
      if(ps(igrid)%is_physical_boundary(5)) then
        ixB3=ixMlo3-1
        do ix2=ixMlo2,ixMhi2
          do ix1=ixMlo1,ixMhi1
            b=bvec(ix1,ix2,ixB3,:)
            call compute_curl_div_node(bvec,igrid,ix1,ix2,ixB3,.true.,j,divb)
            b2=dot_product(b,b)
            call nlfff_compute_omega_cell(b,j,divb,oa,ob)
            qs_state(igrid)%w(ix1,ix2,ixB3,1:3)=cross3(oa,b)
            qs_state(igrid)%w(ix1,ix2,ixB3,4)=dot_product(ob,b)
            wcell=cached_cell_weight(igrid,ix1,ix2,ixB3)
            cell_volume=ps(igrid)%dvolume(ix1,ix2,ixMlo3)
            cell_force=0.d0
            cell_div=0.d0
            if(b2>0.d0) then
              cell_force=wcell*b2*dot_product(oa,oa)*cell_volume
              cell_div=wcell*b2*dot_product(ob,ob)*cell_volume
            end if
            if(.not.all(ieee_is_finite(b)) .or. &
               .not.all(ieee_is_finite(j)) .or. &
               .not.ieee_is_finite(divb) .or. &
               .not.all(ieee_is_finite(oa)) .or. &
               .not.all(ieee_is_finite(ob)) .or. &
               .not.ieee_is_finite(cell_force) .or. &
               .not.ieee_is_finite(cell_div)) then
              diagnostic=(/wcell,cell_force,cell_div/)
              call report_nonfinite('lower-boundary functional',igrid,ix1,&
                 ix2,ixB3,b,j,divb,oa,ob,diagnostic)
            end if
            part_force=part_force+cell_force
            part_div=part_div+cell_div
            part_B2=part_B2+wcell*b2*cell_volume
          end do
        end do
      end if
      local_force=local_force+part_force
      local_div=local_div+part_div
      local_B2=local_B2+part_B2
    end do
    call accumulate_timing(optimization_timing%functional_local,timing_start)
    timing_start=MPI_WTIME()
    local_values=(/local_force,local_div,local_B2/)
    call MPI_ALLREDUCE(local_values,global_values,3,MPI_DOUBLE_PRECISION,&
       MPI_SUM,icomm,ierrmpi)
    Lforce=global_values(1)
    Ldiv=global_values(2)
    B2integral=global_values(3)
    call accumulate_timing(optimization_timing%functional_reductions,timing_start)
    L=Lforce+Ldiv
    timing_start=MPI_WTIME()
    call exchange_qs_internal(qs_state)
    call accumulate_timing(optimization_timing%qs_exchange,timing_start)
  end subroutine evaluate_functional

  subroutine form_trial_field(iw_b,buffer_cells,update_preconditioner,step)
    use mpi
    use mod_global_parameters
    use mod_geometry, only: curlvector,divvector,gradient
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    integer, intent(in) :: iw_b(3),buffer_cells
    character(len=*), intent(in) :: update_preconditioner
    double precision, intent(in) :: step
    double precision :: bvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:3)
    double precision :: current(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:3)
    double precision :: divb_array(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)
    double precision :: qvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:3)
    double precision :: scalar_s(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision :: curlq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:3)
    double precision :: weight(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision :: grad_s_array(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:3)
    double precision :: grad_w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:3)
    double precision :: b(3),j(3),oa(3),ob(3),ftilde(3),gw(3)
    double precision :: q(3),curl_q(3),grad_scalar_s(3),divb,wcell
    double precision :: timing_start
    integer :: iigrid,igrid,ix1,ix2,ix3,idir,idirmin

    timing_start=MPI_WTIME()
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      ps1(igrid)%w(ixG^T,iw_b(:))=0.d0
    end do
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      block=>ps(igrid)
      dxlevel(:)=dx(:,1)
      do ix3=ixGlo3,ixGhi3
        do ix2=ixGlo2,ixGhi2
          do ix1=ixGlo1,ixGhi1
            bvec(ix1,ix2,ix3,:)=ps(igrid)%w(ix1,ix2,ix3,iw_b(:))
            qvec(ix1,ix2,ix3,:)=accepted_qs_state(igrid)%w(ix1,ix2,ix3,1:3)
            scalar_s(ix1,ix2,ix3)=accepted_qs_state(igrid)%w(ix1,ix2,ix3,4)
            weight(ix1,ix2,ix3)=cached_cell_weight(igrid,ix1,ix2,ix3)
          end do
        end do
      end do
      current=0.d0
      divb_array=0.d0
      curlq=0.d0
      grad_s_array=0.d0
      grad_w=0.d0
      idirmin=1
      call curlvector(bvec,ixG^LL,ixM^LL,current,idirmin,1,3)
      call divvector(bvec,ixG^LL,ixM^LL,divb_array,1)
      idirmin=1
      call curlvector(qvec,ixG^LL,ixM^LL,curlq,idirmin,1,3)
      do idir=1,3
        call gradient(scalar_s,ixG^LL,ixM^LL,idir,&
           grad_s_array(ixG^T,idir),1)
        call gradient(weight,ixG^LL,ixM^LL,idir,grad_w(ixG^T,idir),1)
      end do
      do ix3=ixMlo3,ixMhi3
        do ix2=ixMlo2,ixMhi2
          do ix1=ixMlo1,ixMhi1
            if(cached_fixed_active_cell(igrid,ix1,ix2,ix3)) cycle
            b=bvec(ix1,ix2,ix3,:)
            j=current(ix1,ix2,ix3,:)
            divb=divb_array(ix1,ix2,ix3)
            q=qvec(ix1,ix2,ix3,:)
            curl_q=curlq(ix1,ix2,ix3,:)
            grad_scalar_s=grad_s_array(ix1,ix2,ix3,:)
            gw=grad_w(ix1,ix2,ix3,:)
            call nlfff_compute_omega_cell(b,j,divb,oa,ob)
            call nlfff_compose_update_cell(b,j,divb,oa,ob,q,&
               scalar_s(ix1,ix2,ix3),curl_q,grad_scalar_s,gw,&
               weight(ix1,ix2,ix3),ftilde)
            if(.not.all(ieee_is_finite(b)) .or. &
               .not.all(ieee_is_finite(j)) .or. &
               .not.ieee_is_finite(divb) .or. &
               .not.all(ieee_is_finite(oa)) .or. &
               .not.all(ieee_is_finite(ob)) .or. &
               .not.all(ieee_is_finite(ftilde))) &
               call report_nonfinite('update force',igrid,ix1,ix2,ix3,&
                  b,j,divb,oa,ob,ftilde)
            ps1(igrid)%w(ix1,ix2,ix3,iw_b(:))=ftilde
          end do
        end do
      end do
    end do
    call accumulate_timing(optimization_timing%update_kernel,timing_start)

    if(trim(update_preconditioner)=='binomial_xy') then
      timing_start=MPI_WTIME()
      call exchange_trial_b_internal(iw_b)
      call accumulate_timing(optimization_timing%binomial_halo_x,timing_start)
      do iigrid=1,igridstail
        igrid=igrids(iigrid)
        trial_qs_state(igrid)%w(ixG^T,1:3)=0.d0
        do ix3=ixMlo3,ixMhi3
          do ix2=ixMlo2,ixMhi2
            do ix1=ixMlo1,ixMhi1
              trial_qs_state(igrid)%w(ix1,ix2,ix3,1:3)=&
                 0.25d0*ps1(igrid)%w(ix1-1,ix2,ix3,iw_b(:))+&
                 0.5d0*ps1(igrid)%w(ix1,ix2,ix3,iw_b(:))+&
                 0.25d0*ps1(igrid)%w(ix1+1,ix2,ix3,iw_b(:))
            end do
          end do
        end do
        do ix3=ixMlo3,ixMhi3
          do ix2=ixMlo2,ixMhi2
            do ix1=ixMlo1,ixMhi1
              ps1(igrid)%w(ix1,ix2,ix3,iw_b(:))=&
                 trial_qs_state(igrid)%w(ix1,ix2,ix3,1:3)
            end do
          end do
        end do
      end do
      timing_start=MPI_WTIME()
      call exchange_trial_b_internal(iw_b)
      call accumulate_timing(optimization_timing%binomial_halo_y,timing_start)
      do iigrid=1,igridstail
        igrid=igrids(iigrid)
        do ix3=ixMlo3,ixMhi3
          do ix2=ixMlo2,ixMhi2
            do ix1=ixMlo1,ixMhi1
              trial_qs_state(igrid)%w(ix1,ix2,ix3,1:3)=&
                 0.25d0*ps1(igrid)%w(ix1,ix2-1,ix3,iw_b(:))+&
                 0.5d0*ps1(igrid)%w(ix1,ix2,ix3,iw_b(:))+&
                 0.25d0*ps1(igrid)%w(ix1,ix2+1,ix3,iw_b(:))
            end do
          end do
        end do
        do ix3=ixMlo3,ixMhi3
          do ix2=ixMlo2,ixMhi2
            do ix1=ixMlo1,ixMhi1
              ps1(igrid)%w(ix1,ix2,ix3,iw_b(:))=&
                 trial_qs_state(igrid)%w(ix1,ix2,ix3,1:3)
            end do
          end do
        end do
      end do
    end if

    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      do ix3=ixMlo3,ixMhi3
        do ix2=ixMlo2,ixMhi2
          do ix1=ixMlo1,ixMhi1
            if(cached_fixed_active_cell(igrid,ix1,ix2,ix3)) then
              ps1(igrid)%w(ix1,ix2,ix3,iw_b(:))=&
                 ps(igrid)%w(ix1,ix2,ix3,iw_b(:))
            else
              ps1(igrid)%w(ix1,ix2,ix3,iw_b(:))=&
                 ps(igrid)%w(ix1,ix2,ix3,iw_b(:))+&
                 step*ps1(igrid)%w(ix1,ix2,ix3,iw_b(:))
            end if
          end do
        end do
      end do
    end do
  end subroutine form_trial_field

  subroutine accept_trial_field(iw_b)
    use mod_global_parameters

    integer, intent(in) :: iw_b(3)
    integer :: iigrid,igrid

    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      ps(igrid)%w(ixG^T,iw_b(:))=ps1(igrid)%w(ixG^T,iw_b(:))
    end do
  end subroutine accept_trial_field

  subroutine swap_qs_states()
    type(state), allocatable :: temporary(:)

    call move_alloc(accepted_qs_state,temporary)
    call move_alloc(trial_qs_state,accepted_qs_state)
    call move_alloc(temporary,trial_qs_state)
  end subroutine swap_qs_states

  !> Replace centered physical-boundary derivatives by second-order inward
  !> one-sided derivatives.  Block and MPI interfaces remain centered.
  subroutine override_fixed_boundary_curl_div(bvec,igrid,current,divb)
    use mod_global_parameters

    integer, intent(in) :: igrid
    double precision, intent(in) :: bvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:3)
    double precision, intent(inout) :: current(ixGlo1:ixGhi1,&
       ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:3)
    double precision, intent(inout) :: divb(ixGlo1:ixGhi1,&
       ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision :: j(3),div
    integer :: ix1,ix2,ix3

    do ix3=ixMlo3,ixMhi3
      do ix2=ixMlo2,ixMhi2
        do ix1=ixMlo1,ixMhi1
          if((ps(igrid)%is_physical_boundary(1) .and. ix1==ixMlo1) .or. &
             (ps(igrid)%is_physical_boundary(2) .and. ix1==ixMhi1) .or. &
             (ps(igrid)%is_physical_boundary(3) .and. ix2==ixMlo2) .or. &
             (ps(igrid)%is_physical_boundary(4) .and. ix2==ixMhi2) .or. &
             (ps(igrid)%is_physical_boundary(6) .and. ix3==ixMhi3)) then
            call compute_curl_div_node(bvec,igrid,ix1,ix2,ix3,.false.,j,div)
            current(ix1,ix2,ix3,:)=j
            divb(ix1,ix2,ix3)=div
          end if
        end do
      end do
    end do
  end subroutine override_fixed_boundary_curl_div

  subroutine compute_curl_div_node(bvec,igrid,ix1,ix2,ix3,bottom_node,&
     current,divb)
    use mod_global_parameters

    integer, intent(in) :: igrid,ix1,ix2,ix3
    logical, intent(in) :: bottom_node
    double precision, intent(in) :: bvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:3)
    double precision, intent(out) :: current(3),divb
    double precision :: db(3,3)
    integer :: ic,idir

    do idir=1,3
      do ic=1,3
        call derivative_at_node(bvec,igrid,ix1,ix2,ix3,ic,idir,&
           bottom_node,db(ic,idir))
      end do
    end do
    current(1)=db(3,2)-db(2,3)
    current(2)=db(1,3)-db(3,1)
    current(3)=db(2,1)-db(1,2)
    divb=db(1,1)+db(2,2)+db(3,3)
  end subroutine compute_curl_div_node

  subroutine derivative_at_node(bvec,igrid,ix1,ix2,ix3,ic,idir,&
     bottom_node,derivative)
    use mod_global_parameters

    integer, intent(in) :: igrid,ix1,ix2,ix3,ic,idir
    logical, intent(in) :: bottom_node
    double precision, intent(in) :: bvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:3)
    double precision, intent(out) :: derivative

    select case(idir)
    case(1)
      if(ps(igrid)%is_physical_boundary(1) .and. ix1==ixMlo1) then
        derivative=nlfff_boundary_derivative(bvec(ix1,ix2,ix3,ic),&
           bvec(ix1+1,ix2,ix3,ic),bvec(ix1+2,ix2,ix3,ic),dx(1,1),.true.)
      else if(ps(igrid)%is_physical_boundary(2) .and. ix1==ixMhi1) then
        derivative=nlfff_boundary_derivative(bvec(ix1,ix2,ix3,ic),&
           bvec(ix1-1,ix2,ix3,ic),bvec(ix1-2,ix2,ix3,ic),dx(1,1),.false.)
      else
        derivative=(bvec(ix1+1,ix2,ix3,ic)-&
                    bvec(ix1-1,ix2,ix3,ic))/(2.d0*dx(1,1))
      end if
    case(2)
      if(ps(igrid)%is_physical_boundary(3) .and. ix2==ixMlo2) then
        derivative=nlfff_boundary_derivative(bvec(ix1,ix2,ix3,ic),&
           bvec(ix1,ix2+1,ix3,ic),bvec(ix1,ix2+2,ix3,ic),dx(2,1),.true.)
      else if(ps(igrid)%is_physical_boundary(4) .and. ix2==ixMhi2) then
        derivative=nlfff_boundary_derivative(bvec(ix1,ix2,ix3,ic),&
           bvec(ix1,ix2-1,ix3,ic),bvec(ix1,ix2-2,ix3,ic),dx(2,1),.false.)
      else
        derivative=(bvec(ix1,ix2+1,ix3,ic)-&
                    bvec(ix1,ix2-1,ix3,ic))/(2.d0*dx(2,1))
      end if
    case(3)
      if(bottom_node) then
        derivative=nlfff_boundary_derivative(bvec(ix1,ix2,ix3,ic),&
           bvec(ix1,ix2,ix3+1,ic),bvec(ix1,ix2,ix3+2,ic),dx(3,1),.true.)
      else if(ps(igrid)%is_physical_boundary(6) .and. ix3==ixMhi3) then
        derivative=nlfff_boundary_derivative(bvec(ix1,ix2,ix3,ic),&
           bvec(ix1,ix2,ix3-1,ic),bvec(ix1,ix2,ix3-2,ic),dx(3,1),.false.)
      else
        derivative=(bvec(ix1,ix2,ix3+1,ic)-&
                    bvec(ix1,ix2,ix3-1,ic))/(2.d0*dx(3,1))
      end if
    end select
  end subroutine derivative_at_node

  pure double precision function nlfff_boundary_derivative(f0,f1,f2,h,at_low)
    double precision, intent(in) :: f0,f1,f2,h
    logical, intent(in) :: at_low

    if(at_low) then
      nlfff_boundary_derivative=(-3.d0*f0+4.d0*f1-f2)/(2.d0*h)
    else
      nlfff_boundary_derivative=(3.d0*f0-4.d0*f1+f2)/(2.d0*h)
    end if
  end function nlfff_boundary_derivative

  pure subroutine nlfff_compute_omega_cell(b,j,divb,omega_a,omega_b)
    double precision, intent(in) :: b(3),j(3),divb
    double precision, intent(out) :: omega_a(3),omega_b(3)
    double precision :: b2

    b2=dot_product(b,b)
    if(b2>0.d0) then
      omega_a=cross3(j,b)/b2
      omega_b=divb*b/b2
    else
      omega_a=0.d0
      omega_b=0.d0
    end if
  end subroutine nlfff_compute_omega_cell

  !> Pointwise terms in Equations (13) and (15) of Wiegelmann (2004).
  !> Derivatives and Q=Omega_a x B, s=Omega_b dot B are supplied by the
  !> caller so this kernel stays independent of a particular stencil.
  pure subroutine nlfff_compose_update_cell(b,j,divb,omega_a,omega_b,q,s,&
     curl_q,grad_s,grad_w,weight,ftilde)
    double precision, intent(in) :: b(3),j(3),divb
    double precision, intent(in) :: omega_a(3),omega_b(3),q(3),s
    double precision, intent(in) :: curl_q(3),grad_s(3),grad_w(3),weight
    double precision, intent(out) :: ftilde(3)
    double precision :: f(3)

    f=curl_q-cross3(omega_a,j)+grad_s-omega_b*divb+&
       (dot_product(omega_a,omega_a)+dot_product(omega_b,omega_b))*b
    ftilde=weight*f+cross3(q,grad_w)+s*grad_w
  end subroutine nlfff_compose_update_cell

  pure subroutine nlfff_dimensionless_metrics(Lforce,Ldiv,B2integral,&
     characteristic_spacing,epsilon_force,epsilon_div)
    double precision, intent(in) :: Lforce,Ldiv,B2integral
    double precision, intent(in) :: characteristic_spacing
    double precision, intent(out) :: epsilon_force,epsilon_div

    if(B2integral>0.d0) then
      epsilon_force=characteristic_spacing*sqrt(max(0.d0,Lforce)/B2integral)
      epsilon_div=characteristic_spacing*sqrt(max(0.d0,Ldiv)/B2integral)
    else
      epsilon_force=0.d0
      epsilon_div=0.d0
    end if
  end subroutine nlfff_dimensionless_metrics

  !> One-dimensional binomial factor used by the transverse update
  !> preconditioner.  Exposing the three-point kernel keeps its numerical
  !> properties testable without an initialized AMRVAC mesh.
  pure function nlfff_binomial_triplet(left,center,right) result(filtered)
    double precision, intent(in) :: left(3),center(3),right(3)
    double precision :: filtered(3)

    filtered=0.25d0*left+0.5d0*center+0.25d0*right
  end function nlfff_binomial_triplet

  pure function cross3(a,b) result(c)
    double precision, intent(in) :: a(3),b(3)
    double precision :: c(3)

    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
  end function cross3

  pure double precision function nlfff_cosine_side_weight(index,ncell,buffer_cells)
    integer, intent(in) :: index,ncell,buffer_cells
    double precision :: s,pi

    pi=4.d0*datan(1.d0)
    if(index<1 .or. index>ncell) then
      nlfff_cosine_side_weight=0.d0
    else if(index<=buffer_cells) then
      s=dble(index-1)/dble(buffer_cells-1)
      nlfff_cosine_side_weight=0.5d0*(1.d0-dcos(pi*s))
    else if(index>=ncell-buffer_cells+1) then
      s=dble(ncell-index)/dble(buffer_cells-1)
      nlfff_cosine_side_weight=0.5d0*(1.d0-dcos(pi*s))
    else
      nlfff_cosine_side_weight=1.d0
    end if
  end function nlfff_cosine_side_weight

  pure double precision function nlfff_cosine_top_weight(index,ncell,buffer_cells)
    integer, intent(in) :: index,ncell,buffer_cells
    double precision :: s,pi

    pi=4.d0*datan(1.d0)
    if(index>ncell) then
      nlfff_cosine_top_weight=0.d0
    else if(index<1) then
      nlfff_cosine_top_weight=1.d0
    else if(index>=ncell-buffer_cells+1) then
      s=dble(ncell-index)/dble(buffer_cells-1)
      nlfff_cosine_top_weight=0.5d0*(1.d0-dcos(pi*s))
    else
      nlfff_cosine_top_weight=1.d0
    end if
  end function nlfff_cosine_top_weight

  pure double precision function nlfff_weight_product(index1,index2,index3,&
     ncell1,ncell2,ncell3,buffer_cells)
    integer, intent(in) :: index1,index2,index3
    integer, intent(in) :: ncell1,ncell2,ncell3,buffer_cells

    nlfff_weight_product=&
       nlfff_cosine_side_weight(index1,ncell1,buffer_cells)*&
       nlfff_cosine_side_weight(index2,ncell2,buffer_cells)*&
       nlfff_cosine_top_weight(index3,ncell3,buffer_cells)
  end function nlfff_weight_product

  pure logical function nlfff_fixed_active_index(index1,index2,index3,&
     ncell1,ncell2,ncell3)
    integer, intent(in) :: index1,index2,index3,ncell1,ncell2,ncell3

    nlfff_fixed_active_index=(index3==ncell3 .or. index1==1 .or. &
       index1==ncell1 .or. index2==1 .or. index2==ncell2)
  end function nlfff_fixed_active_index

  pure double precision function cell_weight(x1,x2,x3,buffer_cells)
    use mod_global_parameters, only: xprobmin1,xprobmin2,xprobmin3,dx,&
       domain_nx1,domain_nx2,domain_nx3

    double precision, intent(in) :: x1,x2,x3
    integer, intent(in) :: buffer_cells
    integer :: ig1,ig2,ig3

    ig1=nint((x1-xprobmin1)/dx(1,1)+0.5d0)
    ig2=nint((x2-xprobmin2)/dx(2,1)+0.5d0)
    ig3=nint((x3-xprobmin3)/dx(3,1)+0.5d0)
    cell_weight=nlfff_weight_product(ig1,ig2,ig3,domain_nx1,domain_nx2,&
       domain_nx3,buffer_cells)
  end function cell_weight

  pure logical function fixed_active_cell(x1,x2,x3)
    use mod_global_parameters, only: xprobmin1,xprobmin2,xprobmin3,dx,&
       domain_nx1,domain_nx2,domain_nx3

    double precision, intent(in) :: x1,x2,x3
    integer :: ig1,ig2,ig3

    ig1=nint((x1-xprobmin1)/dx(1,1)+0.5d0)
    ig2=nint((x2-xprobmin2)/dx(2,1)+0.5d0)
    ig3=nint((x3-xprobmin3)/dx(3,1)+0.5d0)
    fixed_active_cell=nlfff_fixed_active_index(ig1,ig2,ig3,domain_nx1,&
       domain_nx2,domain_nx3)
  end function fixed_active_cell

  subroutine report_nonfinite(quantity,igrid,ix1,ix2,ix3,b,j,divb,oa,ob,extra)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters, only: mype

    character(len=*), intent(in) :: quantity
    integer, intent(in) :: igrid,ix1,ix2,ix3
    double precision, intent(in) :: b(3),j(3),divb,oa(3),ob(3),extra(3)

    write(*,*) 'Non-finite NLFFF ',trim(quantity),' at rank/grid/cell:',&
       mype,igrid,ix1,ix2,ix3
    write(*,*) 'B:',b
    write(*,*) 'curl(B):',j
    write(*,*) 'div(B):',divb
    write(*,*) 'Omega_a:',oa
    write(*,*) 'Omega_b:',ob
    write(*,*) 'stage values:',extra
    call mpistop('non-finite value in NLFFF optimization')
  end subroutine report_nonfinite

  subroutine write_log_row(unit,result,step,L,Lforce,Ldiv,B2integral,&
     epsilon_force,epsilon_div,accepted_trial)
    integer, intent(in) :: unit
    type(nlfff_optimization_result), intent(in) :: result
    double precision, intent(in) :: step,L,Lforce,Ldiv,B2integral
    double precision, intent(in) :: epsilon_force,epsilon_div
    logical, intent(in) :: accepted_trial

    write(unit,'(i0,a,i0,a,i0,7(a,es24.16),a,i0,a,es24.16,a,i0,a,es24.16,a,l1)') result%attempts,',',&
       result%accepted_steps,',',result%rejected_steps,',',step,',',L,',',&
       Lforce,',',Ldiv,',',B2integral,',',epsilon_force,',',epsilon_div,',',&
       result%functional_evaluations,',',result%step_floor,',',&
       result%plateau_count,',',result%final_relative_functional_change,',',&
       accepted_trial
    flush(unit)
  end subroutine write_log_row
}

end module mod_nlfff_optimization
