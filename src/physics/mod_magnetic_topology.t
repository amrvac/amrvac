module mod_magnetic_topology
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan, &
       ieee_is_finite
  use mod_global_parameters, only: ndim,npe,slab_uniform,periodB, &
       xprobmin1,xprobmax1,global_time,ps,iwstart,nwgc,par_files
  {^IFTHREED
  use mod_global_parameters, only: xprobmin2,xprobmax2,xprobmin3,xprobmax3
  }
  use mod_comm_lib, only: mpistop
  use mod_geometry, only: geo_coordinate => coordinate, &
       geo_cartesian => Cartesian, geo_spherical => spherical, &
       geo_cartesian_stretched => Cartesian_stretched
  use mod_trace_field, only: trace_length_result,trace_twist_result, &
       trace_mapping_result,trace_qperp_result,trace_topology_result, &
       trace_field_length_single,trace_field_length_multi, &
       trace_field_twist_single,trace_field_twist_multi, &
       trace_field_mapping_single,trace_field_mapping_multi, &
       trace_field_topology_multi, &
       trace_field_qperp_single,trace_field_qperp_multi, &
       trace_field_rk2_short_boundary_q_multi, &
       trace_field_spherical_rmin_q_multi, &
       trace_field_spherical_rmin_q_qperp_multi, &
       trace_field_spherical_qperp_multi, &
       trace_spherical_curl_cache_build, &
       trace_spherical_curl_cache_clear, &
       trace_spherical_profile_set, &
       trace_spherical_profile_reset, &
       trace_spherical_profile_count_seeds, &
       trace_spherical_profile_report, &
       trace_cartesian_global_min_cell_size, &
       trace_spherical_global_min_cell_size, &
       trace_set_step_control, &
       trace_set_integrator, &
       trace_debug_cartesian_rk45_tangent_q0_multi, &
       trace_rk2_stats_set_enabled, &
       trace_rk2_stats_reset, &
       trace_rk2_stats_report, &
       trace_rk45_stats_set_enabled, &
       trace_rk45_stats_reset, &
       trace_rk45_stats_report, &
       trace_status_active, &
       trace_status_boundary, &
       trace_status_unsupported_geometry, &
       trace_face_none,trace_face_xmin,trace_face_xmax, &
       trace_face_ymin,trace_face_ymax,trace_face_zmin, &
       trace_face_zmax,trace_face_ambiguous
  implicit none
  private

  public :: mt_length_single,mt_length_seeds,mt_length_plane_xy
  public :: mt_length_plane_xz,mt_length_plane_yz
  public :: mt_fieldline_products_seeds
  public :: mt_twist_single,mt_twist_seeds,mt_twist_plane_xy
  public :: mt_twist_plane_xz,mt_twist_plane_yz
  public :: mt_mapping_plane_xy
  public :: mt_qperp_plane_xy,mt_qperp_plane_xz,mt_qperp_plane_yz
  public :: mt_qperp_plane_arbitrary
  public :: mt_fieldline_products_plane_arbitrary
  public :: mt_topology_plane_xy,mt_topology_plane_xz,mt_topology_plane_yz
  public :: mt_qsl_plane_vtu_xy,mt_qsl_plane_vtu_xz,mt_qsl_plane_vtu_yz
  public :: mt_write_cartesian_vti_pointdata
  public :: mt_fieldline_products_volume_vti
  public :: mt_params_read,mt_run_topology_task

  integer, parameter :: mt_vti_kind_float64 = 1
  integer, parameter :: mt_vti_kind_int32   = 2
  integer, parameter :: mt_vti_name_len     = 64
  integer, parameter :: mt_task_name_len    = 128
  double precision, parameter :: mt_unset_real = huge(1.d0)

  character(len=mt_task_name_len) :: mt_mode = ''
  character(len=mt_task_name_len) :: mt_output_file = ''
  character(len=mt_task_name_len) :: mt_output_prefix = ''
  character(len=mt_task_name_len) :: mt_seed_file = ''
  character(len=mt_task_name_len) :: mt_seed_surface = ''
  character(len=mt_task_name_len) :: mt_seed_layout = 'endpoint'
  character(len=mt_task_name_len) :: mt_vtk_detail = 'minimal'
  character(len=mt_task_name_len) :: mt_step_control = 'global_cell_fraction'
  character(len=mt_task_name_len) :: mt_trace_integrator = 'rk2'
  double precision :: mt_dL = -1.d0
  double precision :: mt_step_fraction = 0.25d0
  double precision :: mt_dL_min = 0.d0
  double precision :: mt_rk45_atol = 1.d-8
  double precision :: mt_rk45_rtol = 1.d-6
  double precision :: mt_rk45_safety = 0.9d0
  double precision :: mt_rk45_min_shrink = 0.2d0
  double precision :: mt_rk45_max_grow = 5.d0
  double precision :: mt_rk45_tangent_floor = 1.d0
  double precision :: mt_rk45_tangent_rtol = 2.d-5
  integer :: mt_max_steps = -1
  double precision :: mt_max_steps_factor = 2.d0
  double precision :: mt_b_min = -1.d0
  double precision :: mt_seed_coord = mt_unset_real
  double precision :: mt_seed_theta0 = mt_unset_real
  double precision :: mt_seed_phi0 = mt_unset_real
  double precision :: mt_seed_alpha = mt_unset_real
  logical :: mt_compute_length = .true.
  logical :: mt_compute_twist = .false.
  logical :: mt_compute_q = .false.
  logical :: mt_compute_qperp = .false.
  logical :: mt_rk2_step_diagnostic = .false.
  logical :: mt_rk45_step_diagnostic = .false.
  logical :: mt_rk45_tangent_diagnostic = .false.
  logical :: mt_rk2_fusion_diagnostic = .false.
  logical :: mt_write_csv = .false.
  character(len=mt_task_name_len) :: mt_plane = ''
  double precision :: mt_origin(3) = mt_unset_real
  double precision :: mt_e1(3) = mt_unset_real
  double precision :: mt_e2(3) = mt_unset_real
  double precision :: mt_xmin = mt_unset_real
  double precision :: mt_xmax = mt_unset_real
  double precision :: mt_ymin = mt_unset_real
  double precision :: mt_ymax = mt_unset_real
  double precision :: mt_zmin = mt_unset_real
  double precision :: mt_zmax = mt_unset_real
  double precision :: mt_x0 = mt_unset_real
  double precision :: mt_y0 = mt_unset_real
  double precision :: mt_z0 = mt_unset_real
  integer :: mt_nx = -1
  integer :: mt_ny = -1
  integer :: mt_nz = -1
  double precision :: mt_s1min = mt_unset_real
  double precision :: mt_s1max = mt_unset_real
  double precision :: mt_s2min = mt_unset_real
  double precision :: mt_s2max = mt_unset_real
  double precision :: mt_s3min = mt_unset_real
  double precision :: mt_s3max = mt_unset_real
  integer :: mt_n1 = -1
  integer :: mt_n2 = -1
  integer :: mt_n3 = -1
  integer :: mt_chunk_nz = -1
  logical :: mt_profile_spherical = .false.
  logical :: mt_params_loaded = .false.

  type, private :: mt_vti_array_desc
    character(len=mt_vti_name_len) :: name
    integer :: kind
    integer(kind=8) :: nbytes
    integer(kind=8) :: offset
  end type mt_vti_array_desc

  type, private :: mt_volume_products
    double precision, allocatable :: length_total(:)
    double precision, allocatable :: length_backward(:)
    double precision, allocatable :: length_forward(:)
    integer, allocatable :: nstep_backward_length(:)
    integer, allocatable :: nstep_forward_length(:)
    integer, allocatable :: status_backward_length(:)
    integer, allocatable :: status_forward_length(:)
    double precision, allocatable :: twist_total(:)
    double precision, allocatable :: twist_backward(:)
    double precision, allocatable :: twist_forward(:)
    integer, allocatable :: nstep_backward_twist(:)
    integer, allocatable :: nstep_forward_twist(:)
    integer, allocatable :: status_backward_twist(:)
    integer, allocatable :: status_forward_twist(:)
    double precision, allocatable :: q(:)
    double precision, allocatable :: logq(:)
    double precision, allocatable :: N2_q(:)
    double precision, allocatable :: bfactor_q(:)
    double precision, allocatable :: length_forward_q(:)
    double precision, allocatable :: length_backward_q(:)
    double precision, allocatable :: Bseed_norm_q(:)
    double precision, allocatable :: Bf_norm_q(:)
    double precision, allocatable :: Bb_norm_q(:)
    integer, allocatable :: valid_q(:)
    integer, allocatable :: status_q(:)
    integer, allocatable :: face_forward_q(:)
    integer, allocatable :: face_backward_q(:)
    integer, allocatable :: status_forward_q(:)
    integer, allocatable :: status_backward_q(:)
    double precision, allocatable :: qperp(:)
    double precision, allocatable :: logqperp(:)
    double precision, allocatable :: N2(:)
    double precision, allocatable :: bfactor(:)
    double precision, allocatable :: length_forward_qperp(:)
    double precision, allocatable :: length_backward_qperp(:)
    double precision, allocatable :: Bseed_norm(:)
    double precision, allocatable :: Bf_norm(:)
    double precision, allocatable :: Bb_norm(:)
    integer, allocatable :: valid_qperp(:)
    integer, allocatable :: status_qperp(:)
    integer, allocatable :: face_forward_qperp(:)
    integer, allocatable :: face_backward_qperp(:)
    integer, allocatable :: status_forward_qperp(:)
    integer, allocatable :: status_backward_qperp(:)
  end type mt_volume_products

contains

  subroutine mt_params_read(files)
    ! Read the one-task magnetic-topology namelist.
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)

    integer :: n

    namelist /magnetic_topology_list/ mt_mode,mt_output_file, &
         mt_output_prefix,mt_seed_file,mt_vtk_detail, &
         mt_seed_surface,mt_seed_layout,mt_seed_coord, &
         mt_seed_theta0,mt_seed_phi0,mt_seed_alpha, &
         mt_step_control,mt_dL,mt_step_fraction,mt_dL_min, &
         mt_trace_integrator,mt_rk45_atol,mt_rk45_rtol, &
         mt_rk45_safety,mt_rk45_min_shrink,mt_rk45_max_grow, &
         mt_rk45_tangent_floor,mt_rk45_tangent_rtol, &
         mt_max_steps,mt_max_steps_factor,mt_b_min, &
         mt_compute_length,mt_compute_twist,mt_compute_q,mt_compute_qperp, &
         mt_rk2_step_diagnostic,mt_rk45_step_diagnostic, &
         mt_rk45_tangent_diagnostic,mt_rk2_fusion_diagnostic, &
         mt_write_csv,mt_plane, &
         mt_origin,mt_e1,mt_e2,mt_xmin,mt_xmax,mt_ymin,mt_ymax, &
         mt_zmin,mt_zmax,mt_nx,mt_ny,mt_nz,mt_x0,mt_y0,mt_z0, &
         mt_s1min,mt_s1max,mt_n1,mt_s2min,mt_s2max,mt_n2, &
         mt_s3min,mt_s3max,mt_n3, &
         mt_chunk_nz,mt_profile_spherical

    call mt_set_default_params()
    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,magnetic_topology_list,end=111)
111   close(unitpar)
    enddo
    mt_params_loaded=.true.
  end subroutine mt_params_read

  subroutine mt_run_topology_task()
    ! Dispatch the namelist-selected topology/QSL postprocessing task.
    ! This backend is deliberately single-MPI-rank and OpenMP capable.
    use mod_comm_lib, only: mpistop
    use mod_ghostcells_update, only: getbc
    character(len=mt_task_name_len) :: mode
    logical :: report_rk2,report_rk45

    if (npe/=1) call mpistop(&
         'magnetic-topology conversion requires npe=1; use OpenMP threads')
    if (.not.mt_params_loaded) call mt_params_read(par_files)
    call getbc(global_time,0.d0,ps,iwstart,nwgc)

    mode=mt_lowercase(trim(mt_mode))
    call mt_validate_common_params(mode)
    call mt_validate_trace_integrator(mode)
    call mt_apply_trace_step_control()
    call mt_apply_auto_max_steps()
    call trace_set_integrator(mt_trace_integrator,mt_rk45_atol, &
         mt_rk45_rtol,mt_rk45_safety,mt_rk45_min_shrink, &
         mt_rk45_max_grow,mt_rk45_tangent_floor,mt_rk45_tangent_rtol)
    report_rk2=mt_lowercase(trim(mt_trace_integrator))=='rk2' .and. &
         mt_rk2_step_diagnostic
    report_rk45=(mt_lowercase(trim(mt_trace_integrator))=='rk45_cartesian' &
         .or. mt_lowercase(trim(mt_trace_integrator))=='rk45_spherical') &
         .and. mt_rk45_step_diagnostic
    call trace_rk2_stats_set_enabled(report_rk2)
    call trace_rk45_stats_set_enabled(report_rk45)
    if (report_rk2) call trace_rk2_stats_reset()
    if (report_rk45) call trace_rk45_stats_reset()
    if (mt_profile_spherical .and. &
         trim(mode)/='spherical_surface_products' .and. &
         trim(mode)/='spherical_cloud_products') then
      write(*,'(a)') 'mt_run_topology_task: mt_profile_spherical applies '// &
           'only to spherical topology modes'
    endif

    select case (trim(mode))
    case ('axis_plane_full_vtu')
      call mt_run_axis_plane_full_vtu_task()
    case ('volume_vti')
      call mt_run_volume_vti_task()
    case ('arbitrary_plane_products')
      call mt_run_arbitrary_plane_products_task()
    case ('seed_products')
      call mt_run_seed_products_task()
    case ('axis_plane_csv')
      call mt_run_axis_plane_csv_task()
    case ('spherical_surface_products')
      call mt_run_spherical_surface_products_task()
    case ('spherical_cloud_products')
      call mt_run_spherical_cloud_products_task()
    case default
      call mpistop('mt_run_topology_task: unknown mt_mode='//trim(mt_mode))
    end select
    if (report_rk2) call trace_rk2_stats_report(trim(mode))
    call trace_rk2_stats_set_enabled(.false.)
    if (report_rk45) call trace_rk45_stats_report(trim(mode))
    call trace_rk45_stats_set_enabled(.false.)
  end subroutine mt_run_topology_task

  subroutine mt_set_default_params()
    mt_params_loaded=.false.
    mt_mode=''
    mt_output_file=''
    mt_output_prefix=''
    mt_seed_file=''
    mt_seed_surface=''
    mt_seed_layout='endpoint'
    mt_vtk_detail='minimal'
    mt_step_control='global_cell_fraction'
    mt_trace_integrator='rk2'
    mt_dL=-1.d0
    mt_step_fraction=0.25d0
    mt_dL_min=0.d0
    mt_rk45_atol=1.d-8
    mt_rk45_rtol=1.d-6
    mt_rk45_safety=0.9d0
    mt_rk45_min_shrink=0.2d0
    mt_rk45_max_grow=5.d0
    mt_rk45_tangent_floor=1.d0
    mt_rk45_tangent_rtol=2.d-5
    mt_max_steps=-1
    mt_max_steps_factor=2.d0
    mt_b_min=-1.d0
    mt_seed_coord=mt_unset_real
    mt_seed_theta0=mt_unset_real
    mt_seed_phi0=mt_unset_real
    mt_seed_alpha=mt_unset_real
    mt_compute_length=.true.
    mt_compute_twist=.false.
    mt_compute_q=.false.
    mt_compute_qperp=.false.
    mt_rk2_step_diagnostic=.false.
    mt_rk45_step_diagnostic=.false.
    mt_rk45_tangent_diagnostic=.false.
    mt_rk2_fusion_diagnostic=.false.
    mt_write_csv=.false.
    mt_plane=''
    mt_origin=mt_unset_real
    mt_e1=mt_unset_real
    mt_e2=mt_unset_real
    mt_xmin=mt_unset_real
    mt_xmax=mt_unset_real
    mt_ymin=mt_unset_real
    mt_ymax=mt_unset_real
    mt_zmin=mt_unset_real
    mt_zmax=mt_unset_real
    mt_x0=mt_unset_real
    mt_y0=mt_unset_real
    mt_z0=mt_unset_real
    mt_nx=-1
    mt_ny=-1
    mt_nz=-1
    mt_s1min=mt_unset_real
    mt_s1max=mt_unset_real
    mt_s2min=mt_unset_real
    mt_s2max=mt_unset_real
    mt_s3min=mt_unset_real
    mt_s3max=mt_unset_real
    mt_n1=-1
    mt_n2=-1
    mt_n3=-1
    mt_chunk_nz=-1
    mt_profile_spherical=.false.
  end subroutine mt_set_default_params

  subroutine mt_validate_common_params(mode)
    character(len=*), intent(in) :: mode

    if (len_trim(mode)==0) then
      call mpistop('mt_run_topology_task requires mt_mode')
    endif
    select case (trim(mode))
    case ('spherical_surface_products','spherical_cloud_products')
      if (ndim/=3) then
        call mpistop(trim(mode)//' requires 3D spherical geometry')
      endif
      {^IFTHREED
      if (geo_coordinate/=geo_spherical) then
        call mpistop(trim(mode)//' requires 3D spherical geometry')
      endif
      if (periodB(3)) then
        call mpistop(trim(mode)//' does not yet support periodic phi')
      endif
      }
    case ('seed_products')
      if (ndim/=3) then
        call mpistop('seed_products requires 3D Cartesian or spherical geometry')
      endif
      {^IFTHREED
      select case (geo_coordinate)
      case (geo_cartesian,geo_cartesian_stretched,geo_spherical)
      case default
        call mpistop('seed_products requires Cartesian or spherical geometry')
      end select
      if (geo_coordinate==geo_spherical .and. periodB(3)) then
        call mpistop('spherical seed_products does not yet support periodic phi')
      endif
      }
    case ('volume_vti')
      if (ndim/=3) then
        call mpistop('volume_vti requires 3D Cartesian geometry')
      endif
      {^IFTHREED
      select case (geo_coordinate)
      case (geo_cartesian)
        ! The VTI sampling grid is user-defined; tracing may interpolate
        ! through uniform or AMR Cartesian simulation grids.
      case (geo_cartesian_stretched)
        ! The VTI sampling grid is uniform Cartesian; tracing may interpolate
        ! through a stretched Cartesian simulation grid.
      case default
        call mpistop('volume_vti requires Cartesian geometry')
      end select
      }
    case ('axis_plane_full_vtu','axis_plane_csv')
      if (ndim/=3) then
        call mpistop(trim(mode)//' requires 3D Cartesian geometry')
      endif
      {^IFTHREED
      select case (geo_coordinate)
      case (geo_cartesian)
        ! Axis-plane output is a seed-surface product; RK2 tracing supports
        ! both slab-uniform and AMR Cartesian simulation grids.
      case (geo_cartesian_stretched)
      case default
        call mpistop(trim(mode)//' requires Cartesian geometry')
      end select
      }
    case ('arbitrary_plane_products')
      if (ndim/=3) then
        call mpistop('arbitrary_plane_products requires 3D Cartesian geometry')
      endif
      {^IFTHREED
      select case (geo_coordinate)
      case (geo_cartesian,geo_cartesian_stretched)
      case default
        call mpistop('arbitrary_plane_products requires Cartesian geometry')
      end select
      }
    case default
      if (ndim/=3 .or. geo_coordinate/=geo_cartesian .or. &
           .not.slab_uniform) then
        call mpistop('mt_run_topology_task mode requires 3D uniform Cartesian geometry')
      endif
    end select
    if (len_trim(mt_output_file)==0 .and. &
        len_trim(mt_output_prefix)==0) then
      call mpistop('mt_run_topology_task requires mt_output_file or mt_output_prefix')
    endif
    select case (mt_lowercase(trim(mt_step_control)))
    case ('fixed')
      if (mt_dL<=0.d0) then
        call mpistop('mt_step_control=fixed requires mt_dL > 0')
      endif
    case ('cell_fraction')
      select case (geo_coordinate)
      case (geo_cartesian,geo_cartesian_stretched)
      case (geo_spherical)
        select case (trim(mode))
        case ('seed_products','spherical_surface_products', &
              'spherical_cloud_products')
        case default
          call mpistop('mt_step_control=cell_fraction is only supported '// &
               'for spherical topology modes')
        end select
      case default
        call mpistop('mt_step_control=cell_fraction requires Cartesian or '// &
             'supported spherical geometry')
      end select
      if (mt_step_fraction<=0.d0) then
        call mpistop('mt_step_control=cell_fraction requires '// &
             'mt_step_fraction > 0')
      endif
      if (mt_dL_min<0.d0) then
        call mpistop('mt_step_control=cell_fraction requires mt_dL_min >= 0')
      endif
    case ('global_cell_fraction')
      if (mt_step_fraction<=0.d0) then
        call mpistop('mt_step_control=global_cell_fraction requires '// &
             'mt_step_fraction > 0')
      endif
      if (mt_dL_min<0.d0) then
        call mpistop('mt_step_control=global_cell_fraction requires mt_dL_min >= 0')
      endif
    case default
      call mpistop('mt_run_topology_task requires mt_step_control=fixed '// &
           'or cell_fraction/global_cell_fraction')
    end select
    if (mt_max_steps<=0 .and. mt_max_steps_factor<=0.d0) then
      call mpistop('mt_max_steps auto requires mt_max_steps_factor > 0')
    endif
    select case (mt_lowercase(trim(mt_vtk_detail)))
    case ('minimal','full')
    case default
      call mpistop('mt_run_topology_task requires mt_vtk_detail=minimal or full')
    end select
  end subroutine mt_validate_common_params

  subroutine mt_validate_trace_integrator(mode)
    character(len=*), intent(in) :: mode

    character(len=mt_task_name_len) :: integrator

    integrator=mt_lowercase(trim(mt_trace_integrator))
    select case (trim(integrator))
    case ('rk2')
      ! RK2 is the default tracing path; product-specific guards live below.
    case ('rk45_cartesian')
      if (trim(mode)/='seed_products' .and. trim(mode)/='volume_vti' .and. &
           trim(mode)/='arbitrary_plane_products' .and. &
           trim(mode)/='axis_plane_full_vtu' .and. &
           trim(mode)/='axis_plane_csv') then
        write(*,'(a)') 'mt_trace_integrator=rk45_cartesian supports '// &
             'only Cartesian seed_products, volume_vti, arbitrary_plane_products, '// &
             'or axis-plane products'
        flush(6)
        call mpistop('mt_trace_integrator=rk45_cartesian supports '// &
             'only seed_products, volume_vti, arbitrary_plane_products, '// &
             'or axis-plane products')
      endif
      if (ndim/=3) then
        call mpistop('mt_trace_integrator=rk45_cartesian requires 3D Cartesian geometry')
      endif
      {^IFTHREED
      select case (geo_coordinate)
      case (geo_cartesian,geo_cartesian_stretched)
      case default
        write(*,'(a)') 'mt_trace_integrator=rk45_cartesian is Cartesian-only'
        flush(6)
        call mpistop('mt_trace_integrator=rk45_cartesian is Cartesian-only')
      end select
      }
      if (mt_rk45_atol<0.d0 .or. mt_rk45_rtol<0.d0 .or. &
           mt_rk45_atol+mt_rk45_rtol<=0.d0) then
        call mpistop('mt_trace_integrator=rk45_cartesian requires '// &
             'non-negative tolerances with atol+rtol > 0')
      endif
      if (mt_rk45_safety<=0.d0 .or. mt_rk45_safety>1.d0) then
        call mpistop('mt_rk45_safety must be in (0,1]')
      endif
      if (mt_rk45_min_shrink<=0.d0 .or. mt_rk45_min_shrink>1.d0) then
        call mpistop('mt_rk45_min_shrink must be in (0,1]')
      endif
      if (mt_rk45_max_grow<1.d0) then
        call mpistop('mt_rk45_max_grow must be >= 1')
      endif
      if (mt_rk45_tangent_floor<=0.d0) then
        call mpistop('mt_rk45_tangent_floor must be > 0')
      endif
      if (mt_rk45_tangent_rtol<=0.d0) then
        call mpistop('mt_rk45_tangent_rtol must be > 0')
      endif
    case ('rk45_spherical')
      if (trim(mode)/='seed_products' .and. &
           trim(mode)/='spherical_cloud_products' .and. &
           trim(mode)/='spherical_surface_products') then
        write(*,'(a)') 'mt_trace_integrator=rk45_spherical supports '// &
             'only spherical seed_products, spherical_cloud_products, '// &
             'or spherical_surface_products'
        flush(6)
        call mpistop('mt_trace_integrator=rk45_spherical supports '// &
             'only seed_products, spherical_cloud_products, '// &
             'or spherical_surface_products')
      endif
      if (ndim/=3) then
        call mpistop('mt_trace_integrator=rk45_spherical requires 3D spherical geometry')
      endif
      {^IFTHREED
      if (geo_coordinate/=geo_spherical) then
        write(*,'(a)') 'mt_trace_integrator=rk45_spherical is spherical-only'
        flush(6)
        call mpistop('mt_trace_integrator=rk45_spherical is spherical-only')
      endif
      if (periodB(3)) then
        call mpistop('mt_trace_integrator=rk45_spherical does not yet support periodic phi')
      endif
      }
      if (mt_compute_q) then
        {^IFTHREED
        if (geo_coordinate/=geo_spherical) then
          call mpistop('mt_trace_integrator=rk45_spherical standard '// &
               'logQ currently requires spherical geometry')
        endif
        }
      endif
      if (mt_rk45_atol<0.d0 .or. mt_rk45_rtol<0.d0 .or. &
           mt_rk45_atol+mt_rk45_rtol<=0.d0) then
        call mpistop('mt_trace_integrator=rk45_spherical requires '// &
             'non-negative tolerances with atol+rtol > 0')
      endif
      if (mt_rk45_safety<=0.d0 .or. mt_rk45_safety>1.d0) then
        call mpistop('mt_rk45_safety must be in (0,1]')
      endif
      if (mt_rk45_min_shrink<=0.d0 .or. mt_rk45_min_shrink>1.d0) then
        call mpistop('mt_rk45_min_shrink must be in (0,1]')
      endif
      if (mt_rk45_max_grow<1.d0) then
        call mpistop('mt_rk45_max_grow must be >= 1')
      endif
      if (mt_rk45_tangent_floor<=0.d0) then
        call mpistop('mt_rk45_tangent_floor must be > 0')
      endif
      if (mt_rk45_tangent_rtol<=0.d0) then
        call mpistop('mt_rk45_tangent_rtol must be > 0')
      endif
    case default
      call mpistop('mt_trace_integrator must be rk2, rk45_cartesian, or rk45_spherical')
    end select
    if (mt_rk45_tangent_diagnostic) then
      if (trim(mode)/='seed_products') then
        call mpistop('mt_rk45_tangent_diagnostic requires seed_products')
      endif
      if (trim(integrator)/='rk45_cartesian') then
        call mpistop('mt_rk45_tangent_diagnostic requires '// &
             'mt_trace_integrator=rk45_cartesian')
      endif
      if (mt_compute_q .or. mt_compute_qperp) then
        call mpistop('mt_rk45_tangent_diagnostic is q0 diagnostic-only; '// &
             'disable public Q and Qperp')
      endif
      if (ndim/=3 .or. .not.slab_uniform) then
        call mpistop('mt_rk45_tangent_diagnostic requires 3D uniform Cartesian geometry')
      endif
      {^IFTHREED
      if (geo_coordinate/=geo_cartesian) then
        call mpistop('mt_rk45_tangent_diagnostic requires Cartesian geometry')
      endif
      }
    endif
    if (mt_rk2_fusion_diagnostic) then
      if (trim(mode)/='seed_products') then
        call mpistop('mt_rk2_fusion_diagnostic requires seed_products')
      endif
      if (trim(integrator)/='rk2') then
        call mpistop('mt_rk2_fusion_diagnostic requires '// &
             'mt_trace_integrator=rk2')
      endif
      if (.not.mt_compute_twist) then
        call mpistop('mt_rk2_fusion_diagnostic requires twist diagnostics')
      endif
      if (mt_compute_q .or. mt_compute_qperp) then
        call mpistop('mt_rk2_fusion_diagnostic is q0 diagnostic-only; '// &
             'disable public Q and Qperp')
      endif
    endif
  end subroutine mt_validate_trace_integrator

  subroutine mt_apply_trace_step_control()
    character(len=mt_task_name_len) :: step_mode
    double precision :: hmin,dL_eff,dL_cap
    integer :: status

    step_mode=mt_lowercase(trim(mt_step_control))
    select case (trim(step_mode))
    case ('cell_fraction')
      if ((geo_coordinate==geo_cartesian .or. &
           geo_coordinate==geo_cartesian_stretched .or. &
           geo_coordinate==geo_spherical) .and. mt_dL<=0.d0) then
        dL_cap=mt_domain_diagonal()
        if (dL_cap<=0.d0) then
          call mpistop('mt_step_control=cell_fraction could not determine '// &
               'domain step cap')
        endif
        mt_dL=dL_cap
        write(*,'(a,es16.8,a,es16.8,a,es16.8)') &
             'mt_step_control=cell_fraction: mt_step_fraction=', &
             mt_step_fraction,', automatic dL_cap=',mt_dL, &
             ', mt_dL_min=',mt_dL_min
      else
        write(*,'(a,es16.8,a,es16.8,a,es16.8)') &
             'mt_step_control=cell_fraction: mt_step_fraction=', &
             mt_step_fraction,', mt_dL_cap=',mt_dL, &
             ', mt_dL_min=',mt_dL_min
      endif
      call trace_set_step_control('cell_fraction',mt_step_fraction,mt_dL_min)
    case ('global_cell_fraction')
      select case (geo_coordinate)
      case (geo_cartesian,geo_cartesian_stretched)
        call trace_cartesian_global_min_cell_size(hmin,status)
      case (geo_spherical)
        call trace_spherical_global_min_cell_size(hmin,status)
      case default
        status=trace_status_unsupported_geometry
        hmin=-1.d0
      end select
      if (status/=trace_status_active .or. hmin<=0.d0) then
        call mpistop('mt_step_control=global_cell_fraction could not determine h_global_min')
      endif
      dL_eff=mt_step_fraction*hmin
      if (mt_dL>0.d0) dL_eff=min(dL_eff,mt_dL)
      if (mt_dL_min>0.d0) then
        if (mt_dL>0.d0) then
          dL_eff=max(dL_eff,min(mt_dL_min,mt_dL))
        else
          dL_eff=max(dL_eff,mt_dL_min)
        endif
      endif
      if (dL_eff<=0.d0) then
        call mpistop('mt_step_control=global_cell_fraction produced non-positive mt_dL')
      endif
      if (mt_dL>0.d0) then
        write(*,'(a,es16.8,a,es16.8,a,es16.8,a,es16.8)') &
             'mt_step_control=global_cell_fraction: h_global_min=',hmin, &
             ', mt_step_fraction=',mt_step_fraction,', mt_dL_cap=',mt_dL, &
             ', effective_mt_dL=',dL_eff
      else
        write(*,'(a,es16.8,a,es16.8,a,es16.8)') &
             'mt_step_control=global_cell_fraction: h_global_min=',hmin, &
             ', mt_step_fraction=',mt_step_fraction,', effective_mt_dL=',dL_eff
      endif
      mt_dL=dL_eff
      call trace_set_step_control('fixed',mt_step_fraction,mt_dL_min)
    case default
      call trace_set_step_control('fixed',mt_step_fraction,mt_dL_min)
    end select
  end subroutine mt_apply_trace_step_control

  subroutine mt_apply_auto_max_steps()
    character(len=mt_task_name_len) :: step_mode
    double precision :: hmin,step_est,Ldiag,nsteps_real
    integer :: status

    if (mt_max_steps>0) return

    Ldiag=mt_domain_diagonal()
    if (Ldiag<=0.d0) then
      call mpistop('mt_max_steps auto could not determine domain diagonal')
    endif

    step_mode=mt_lowercase(trim(mt_step_control))
    select case (trim(step_mode))
    case ('fixed')
      step_est=mt_dL
    case ('cell_fraction','global_cell_fraction')
      select case (geo_coordinate)
      case (geo_cartesian,geo_cartesian_stretched)
        call trace_cartesian_global_min_cell_size(hmin,status)
      case (geo_spherical)
        call trace_spherical_global_min_cell_size(hmin,status)
      case default
        status=trace_status_unsupported_geometry
        hmin=-1.d0
      end select
      if (status/=trace_status_active .or. hmin<=0.d0) then
        call mpistop('mt_max_steps auto could not determine h_global_min')
      endif
      step_est=mt_step_fraction*hmin
      if (mt_dL>0.d0) step_est=min(step_est,mt_dL)
    case default
      step_est=mt_dL
    end select

    if (step_est<=0.d0) then
      call mpistop('mt_max_steps auto produced non-positive step estimate')
    endif
    nsteps_real=mt_max_steps_factor*Ldiag/step_est
    if (nsteps_real>dble(huge(mt_max_steps)-1)) then
      call mpistop('mt_max_steps auto estimate exceeds integer range')
    endif
    mt_max_steps=max(1,ceiling(nsteps_real))
    write(*,'(a,es16.8,a,es16.8,a,es16.8,a,es16.8,a,i0)') &
         'mt_max_steps auto: Ldiag=',Ldiag,', step_est=',step_est, &
         ', factor=',mt_max_steps_factor,', nsteps_real=',nsteps_real, &
         ', mt_max_steps=',mt_max_steps
  end subroutine mt_apply_auto_max_steps

  double precision function mt_domain_diagonal() result(Ldiag)
    integer :: i,j,k,n,m
    double precision :: r,theta,phi,dist
    double precision :: corners(8,3),dx(3)

    Ldiag=abs(xprobmax1-xprobmin1)
    {^IFTHREED
    if (geo_coordinate==geo_spherical) then
      n=0
      do i=0,1
        if (i==0) then
          r=xprobmin1
        else
          r=xprobmax1
        endif
        do j=0,1
          if (j==0) then
            theta=xprobmin2
          else
            theta=xprobmax2
          endif
          do k=0,1
            if (k==0) then
              phi=xprobmin3
            else
              phi=xprobmax3
            endif
            n=n+1
            corners(n,1)=r*dsin(theta)*dcos(phi)
            corners(n,2)=r*dsin(theta)*dsin(phi)
            corners(n,3)=r*dcos(theta)
          end do
        end do
      end do
      Ldiag=0.d0
      do i=1,n-1
        do j=i+1,n
          do m=1,3
            dx(m)=corners(i,m)-corners(j,m)
          end do
          dist=dsqrt(sum(dx*dx))
          Ldiag=max(Ldiag,dist)
        end do
      end do
    else
      Ldiag=dsqrt((xprobmax1-xprobmin1)**2+ &
           (xprobmax2-xprobmin2)**2+(xprobmax3-xprobmin3)**2)
    endif
    }
  end function mt_domain_diagonal

  subroutine mt_run_axis_plane_full_vtu_task()
    character(len=mt_task_name_len) :: plane
    character(len=mt_task_name_len) :: output_file
    logical :: minimal_output

    plane=mt_lowercase(trim(mt_plane))
    minimal_output=.not.mt_vtk_detail_is_full()
    if (minimal_output) then
      call mt_resolve_output_file('axis_plane_full_vtu','.vti', &
           '_'//trim(plane)//'_minimal.vti',output_file)
    else
      call mt_resolve_output_file('axis_plane_full_vtu','.vtu', &
           '_'//trim(plane)//'_full.vtu',output_file)
    endif
    if (.not.minimal_output .and. &
         (mt_compute_twist .or. mt_compute_q .or. mt_compute_qperp)) then
      write(*,'(a)') 'mt_run_topology_task: axis_plane_full_vtu ignores compute flags'
      write(*,'(a)') 'mt_run_topology_task: it always writes full plane products'
    endif

    select case (trim(plane))
    case ('xy')
      call mt_require_real('mt_xmin',mt_xmin)
      call mt_require_real('mt_xmax',mt_xmax)
      call mt_require_real('mt_ymin',mt_ymin)
      call mt_require_real('mt_ymax',mt_ymax)
      call mt_require_real('mt_z0',mt_z0)
      call mt_require_positive_int('mt_nx',mt_nx)
      call mt_require_positive_int('mt_ny',mt_ny)
      call mt_require_ordered('mt_xmin','mt_xmax',mt_xmin,mt_xmax)
      call mt_require_ordered('mt_ymin','mt_ymax',mt_ymin,mt_ymax)
      if (minimal_output) then
        write(*,'(a)') 'mt_run_topology_task: writing axis-plane minimal VTI '//trim(output_file)
        if (mt_b_min>0.d0) then
          call mt_qsl_plane_vti_xy(mt_xmin,mt_xmax,mt_nx,mt_ymin,mt_ymax, &
               mt_ny,mt_z0,mt_dL,mt_max_steps,trim(output_file), &
               mt_compute_length,mt_compute_twist,mt_compute_q, &
               mt_compute_qperp,b_min=mt_b_min)
        else
          call mt_qsl_plane_vti_xy(mt_xmin,mt_xmax,mt_nx,mt_ymin,mt_ymax, &
               mt_ny,mt_z0,mt_dL,mt_max_steps,trim(output_file), &
               mt_compute_length,mt_compute_twist,mt_compute_q, &
               mt_compute_qperp)
        endif
      else
        write(*,'(a)') 'mt_run_topology_task: writing axis-plane full VTU '//trim(output_file)
        if (mt_b_min>0.d0) then
          call mt_qsl_plane_vtu_xy(mt_xmin,mt_xmax,mt_nx,mt_ymin,mt_ymax, &
               mt_ny,mt_z0,mt_dL,mt_max_steps,trim(output_file), &
               b_min=mt_b_min)
        else
          call mt_qsl_plane_vtu_xy(mt_xmin,mt_xmax,mt_nx,mt_ymin,mt_ymax, &
               mt_ny,mt_z0,mt_dL,mt_max_steps,trim(output_file))
        endif
      endif
    case ('xz')
      call mt_require_real('mt_xmin',mt_xmin)
      call mt_require_real('mt_xmax',mt_xmax)
      call mt_require_real('mt_zmin',mt_zmin)
      call mt_require_real('mt_zmax',mt_zmax)
      call mt_require_real('mt_y0',mt_y0)
      call mt_require_positive_int('mt_nx',mt_nx)
      call mt_require_positive_int('mt_nz',mt_nz)
      call mt_require_ordered('mt_xmin','mt_xmax',mt_xmin,mt_xmax)
      call mt_require_ordered('mt_zmin','mt_zmax',mt_zmin,mt_zmax)
      if (minimal_output) then
        write(*,'(a)') 'mt_run_topology_task: writing axis-plane minimal VTI '//trim(output_file)
        if (mt_b_min>0.d0) then
          call mt_qsl_plane_vti_xz(mt_xmin,mt_xmax,mt_nx,mt_zmin,mt_zmax, &
               mt_nz,mt_y0,mt_dL,mt_max_steps,trim(output_file), &
               mt_compute_length,mt_compute_twist,mt_compute_q, &
               mt_compute_qperp,b_min=mt_b_min)
        else
          call mt_qsl_plane_vti_xz(mt_xmin,mt_xmax,mt_nx,mt_zmin,mt_zmax, &
               mt_nz,mt_y0,mt_dL,mt_max_steps,trim(output_file), &
               mt_compute_length,mt_compute_twist,mt_compute_q, &
               mt_compute_qperp)
        endif
      else
        write(*,'(a)') 'mt_run_topology_task: writing axis-plane full VTU '//trim(output_file)
        if (mt_b_min>0.d0) then
          call mt_qsl_plane_vtu_xz(mt_xmin,mt_xmax,mt_nx,mt_zmin,mt_zmax, &
               mt_nz,mt_y0,mt_dL,mt_max_steps,trim(output_file), &
               b_min=mt_b_min)
        else
          call mt_qsl_plane_vtu_xz(mt_xmin,mt_xmax,mt_nx,mt_zmin,mt_zmax, &
               mt_nz,mt_y0,mt_dL,mt_max_steps,trim(output_file))
        endif
      endif
    case ('yz')
      call mt_require_real('mt_ymin',mt_ymin)
      call mt_require_real('mt_ymax',mt_ymax)
      call mt_require_real('mt_zmin',mt_zmin)
      call mt_require_real('mt_zmax',mt_zmax)
      call mt_require_real('mt_x0',mt_x0)
      call mt_require_positive_int('mt_ny',mt_ny)
      call mt_require_positive_int('mt_nz',mt_nz)
      call mt_require_ordered('mt_ymin','mt_ymax',mt_ymin,mt_ymax)
      call mt_require_ordered('mt_zmin','mt_zmax',mt_zmin,mt_zmax)
      if (minimal_output) then
        write(*,'(a)') 'mt_run_topology_task: writing axis-plane minimal VTI '//trim(output_file)
        if (mt_b_min>0.d0) then
          call mt_qsl_plane_vti_yz(mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax, &
               mt_nz,mt_x0,mt_dL,mt_max_steps,trim(output_file), &
               mt_compute_length,mt_compute_twist,mt_compute_q, &
               mt_compute_qperp,b_min=mt_b_min)
        else
          call mt_qsl_plane_vti_yz(mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax, &
               mt_nz,mt_x0,mt_dL,mt_max_steps,trim(output_file), &
               mt_compute_length,mt_compute_twist,mt_compute_q, &
               mt_compute_qperp)
        endif
      else
        write(*,'(a)') 'mt_run_topology_task: writing axis-plane full VTU '//trim(output_file)
        if (mt_b_min>0.d0) then
          call mt_qsl_plane_vtu_yz(mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax, &
               mt_nz,mt_x0,mt_dL,mt_max_steps,trim(output_file), &
               b_min=mt_b_min)
        else
          call mt_qsl_plane_vtu_yz(mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax, &
               mt_nz,mt_x0,mt_dL,mt_max_steps,trim(output_file))
        endif
      endif
    case default
      call mpistop('mt_run_topology_task: axis_plane_full_vtu requires mt_plane=xy, xz, or yz')
    end select
  end subroutine mt_run_axis_plane_full_vtu_task

  subroutine mt_run_axis_plane_csv_task()
    character(len=mt_task_name_len) :: plane
    character(len=mt_task_name_len) :: length_csv,twist_csv
    character(len=mt_task_name_len) :: q_csv
    character(len=mt_task_name_len) :: qperp_csv

    if (len_trim(mt_output_prefix)==0) then
      call mpistop('axis_plane_csv requires mt_output_prefix')
    endif
    if (len_trim(mt_output_file)>0) then
      write(*,'(a)') 'mt_run_topology_task: axis_plane_csv ignores mt_output_file'
      write(*,'(a)') 'mt_run_topology_task: axis_plane_csv uses mt_output_prefix for CSV outputs'
    endif

    plane=mt_lowercase(trim(mt_plane))
    call mt_resolve_prefix_file('_'//trim(plane)//'_length.csv', &
         length_csv,'axis_plane_csv length requires mt_output_prefix')
    twist_csv=''
    q_csv=''
    qperp_csv=''
    if (mt_compute_twist) then
      call mt_resolve_prefix_file('_'//trim(plane)//'_twist.csv', &
           twist_csv,'axis_plane_csv twist requires mt_output_prefix')
    endif
    if (mt_compute_q) then
      call mt_resolve_prefix_file('_'//trim(plane)//'_q.csv', &
           q_csv,'axis_plane_csv Q requires mt_output_prefix')
    endif
    if (mt_compute_qperp) then
      call mt_resolve_prefix_file('_'//trim(plane)//'_qperp_method2.csv', &
           qperp_csv,'axis_plane_csv Qperp requires mt_output_prefix')
    endif

    select case (trim(plane))
    case ('xy')
      call mt_require_real('mt_xmin',mt_xmin)
      call mt_require_real('mt_xmax',mt_xmax)
      call mt_require_real('mt_ymin',mt_ymin)
      call mt_require_real('mt_ymax',mt_ymax)
      call mt_require_real('mt_z0',mt_z0)
      call mt_require_positive_int('mt_nx',mt_nx)
      call mt_require_positive_int('mt_ny',mt_ny)
      call mt_require_ordered('mt_xmin','mt_xmax',mt_xmin,mt_xmax)
      call mt_require_ordered('mt_ymin','mt_ymax',mt_ymin,mt_ymax)
      if (mt_b_min>0.d0) then
        call mt_axis_plane_products_csv_axis(mt_xmin,mt_xmax,mt_nx, &
             mt_ymin,mt_ymax,mt_ny,mt_z0,1,2,3,mt_dL,mt_max_steps, &
             trim(length_csv),trim(twist_csv),trim(q_csv), &
             trim(qperp_csv),'mt_axis_plane_products_csv_xy','ix,iy', &
             'ix,iy',b_min=mt_b_min)
      else
        call mt_axis_plane_products_csv_axis(mt_xmin,mt_xmax,mt_nx, &
             mt_ymin,mt_ymax,mt_ny,mt_z0,1,2,3,mt_dL,mt_max_steps, &
             trim(length_csv),trim(twist_csv),trim(q_csv), &
             trim(qperp_csv),'mt_axis_plane_products_csv_xy','ix,iy', &
             'ix,iy')
      endif
    case ('xz')
      call mt_require_real('mt_xmin',mt_xmin)
      call mt_require_real('mt_xmax',mt_xmax)
      call mt_require_real('mt_zmin',mt_zmin)
      call mt_require_real('mt_zmax',mt_zmax)
      call mt_require_real('mt_y0',mt_y0)
      call mt_require_positive_int('mt_nx',mt_nx)
      call mt_require_positive_int('mt_nz',mt_nz)
      call mt_require_ordered('mt_xmin','mt_xmax',mt_xmin,mt_xmax)
      call mt_require_ordered('mt_zmin','mt_zmax',mt_zmin,mt_zmax)
      if (mt_b_min>0.d0) then
        call mt_axis_plane_products_csv_axis(mt_xmin,mt_xmax,mt_nx, &
             mt_zmin,mt_zmax,mt_nz,mt_y0,1,3,2,mt_dL,mt_max_steps, &
             trim(length_csv),trim(twist_csv),trim(q_csv), &
             trim(qperp_csv),'mt_axis_plane_products_csv_xz','ix,iz', &
             'i,j',b_min=mt_b_min)
      else
        call mt_axis_plane_products_csv_axis(mt_xmin,mt_xmax,mt_nx, &
             mt_zmin,mt_zmax,mt_nz,mt_y0,1,3,2,mt_dL,mt_max_steps, &
             trim(length_csv),trim(twist_csv),trim(q_csv), &
             trim(qperp_csv),'mt_axis_plane_products_csv_xz','ix,iz', &
             'i,j')
      endif
    case ('yz')
      call mt_require_real('mt_ymin',mt_ymin)
      call mt_require_real('mt_ymax',mt_ymax)
      call mt_require_real('mt_zmin',mt_zmin)
      call mt_require_real('mt_zmax',mt_zmax)
      call mt_require_real('mt_x0',mt_x0)
      call mt_require_positive_int('mt_ny',mt_ny)
      call mt_require_positive_int('mt_nz',mt_nz)
      call mt_require_ordered('mt_ymin','mt_ymax',mt_ymin,mt_ymax)
      call mt_require_ordered('mt_zmin','mt_zmax',mt_zmin,mt_zmax)
      if (mt_b_min>0.d0) then
        call mt_axis_plane_products_csv_axis(mt_ymin,mt_ymax,mt_ny, &
             mt_zmin,mt_zmax,mt_nz,mt_x0,2,3,1,mt_dL,mt_max_steps, &
             trim(length_csv),trim(twist_csv),trim(q_csv), &
             trim(qperp_csv),'mt_axis_plane_products_csv_yz','iy,iz', &
             'i,j',b_min=mt_b_min)
      else
        call mt_axis_plane_products_csv_axis(mt_ymin,mt_ymax,mt_ny, &
             mt_zmin,mt_zmax,mt_nz,mt_x0,2,3,1,mt_dL,mt_max_steps, &
             trim(length_csv),trim(twist_csv),trim(q_csv), &
             trim(qperp_csv),'mt_axis_plane_products_csv_yz','iy,iz', &
             'i,j')
      endif
    case default
      call mpistop('mt_run_topology_task: axis_plane_csv requires mt_plane=xy, xz, or yz')
    end select

    write(*,'(a)') 'mt_run_topology_task: wrote axis-plane CSV length '// &
         trim(length_csv)
    if (mt_compute_twist) write(*,'(a)') &
         'mt_run_topology_task: wrote axis-plane CSV twist '//trim(twist_csv)
    if (mt_compute_q) write(*,'(a)') &
         'mt_run_topology_task: wrote axis-plane CSV Q '//trim(q_csv)
    if (mt_compute_qperp) write(*,'(a)') &
         'mt_run_topology_task: wrote axis-plane CSV Qperp '//trim(qperp_csv)
  end subroutine mt_run_axis_plane_csv_task

  subroutine mt_run_volume_vti_task()
    character(len=mt_task_name_len) :: output_file

    call mt_resolve_output_file('volume_vti','.vti','_volume.vti', &
         output_file)
    call mt_require_real('mt_xmin',mt_xmin)
    call mt_require_real('mt_xmax',mt_xmax)
    call mt_require_real('mt_ymin',mt_ymin)
    call mt_require_real('mt_ymax',mt_ymax)
    call mt_require_real('mt_zmin',mt_zmin)
    call mt_require_real('mt_zmax',mt_zmax)
    call mt_require_positive_int('mt_nx',mt_nx)
    call mt_require_positive_int('mt_ny',mt_ny)
    call mt_require_positive_int('mt_nz',mt_nz)
    call mt_require_ordered('mt_xmin','mt_xmax',mt_xmin,mt_xmax)
    call mt_require_ordered('mt_ymin','mt_ymax',mt_ymin,mt_ymax)
    call mt_require_ordered('mt_zmin','mt_zmax',mt_zmin,mt_zmax)
    write(*,'(a)') 'mt_run_topology_task: writing volume VTI '//trim(output_file)
    if (mt_b_min>0.d0) then
      if (mt_chunk_nz>0) then
        call mt_fieldline_products_volume_vti(mt_xmin,mt_xmax,mt_nx, &
             mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax,mt_nz, &
             mt_dL,mt_max_steps,trim(output_file),b_min=mt_b_min, &
           compute_twist=mt_compute_twist, &
             compute_length=mt_compute_length, &
             compute_q=mt_compute_q,compute_qperp=mt_compute_qperp, &
             chunk_nz=mt_chunk_nz)
      else
        call mt_fieldline_products_volume_vti(mt_xmin,mt_xmax,mt_nx, &
             mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax,mt_nz, &
             mt_dL,mt_max_steps,trim(output_file),b_min=mt_b_min, &
             compute_twist=mt_compute_twist, &
             compute_length=mt_compute_length, &
             compute_q=mt_compute_q,compute_qperp=mt_compute_qperp)
      endif
    else
      if (mt_chunk_nz>0) then
        call mt_fieldline_products_volume_vti(mt_xmin,mt_xmax,mt_nx, &
             mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax,mt_nz, &
             mt_dL,mt_max_steps,trim(output_file), &
             compute_twist=mt_compute_twist, &
             compute_length=mt_compute_length, &
             compute_q=mt_compute_q,compute_qperp=mt_compute_qperp, &
             chunk_nz=mt_chunk_nz)
      else
        call mt_fieldline_products_volume_vti(mt_xmin,mt_xmax,mt_nx, &
             mt_ymin,mt_ymax,mt_ny,mt_zmin,mt_zmax,mt_nz, &
             mt_dL,mt_max_steps,trim(output_file), &
             compute_twist=mt_compute_twist, &
             compute_length=mt_compute_length, &
             compute_q=mt_compute_q,compute_qperp=mt_compute_qperp)
      endif
    endif
  end subroutine mt_run_volume_vti_task

  subroutine mt_run_arbitrary_plane_products_task()
    character(len=mt_task_name_len) :: vtu_file,csv_file

    call mt_resolve_output_file('arbitrary_plane_products','.vtu', &
         '_arbitrary_plane_products.vtu',vtu_file)
    csv_file=''
    if (mt_write_csv) then
      call mt_resolve_prefix_file('_arbitrary_plane_products.csv',csv_file, &
           'arbitrary_plane_products CSV requires mt_output_prefix')
    endif
    call mt_require_real('mt_origin(1)',mt_origin(1))
    call mt_require_real('mt_origin(2)',mt_origin(2))
    call mt_require_real('mt_origin(3)',mt_origin(3))
    call mt_require_real('mt_e1(1)',mt_e1(1))
    call mt_require_real('mt_e1(2)',mt_e1(2))
    call mt_require_real('mt_e1(3)',mt_e1(3))
    call mt_require_real('mt_e2(1)',mt_e2(1))
    call mt_require_real('mt_e2(2)',mt_e2(2))
    call mt_require_real('mt_e2(3)',mt_e2(3))
    call mt_require_real('mt_s1min',mt_s1min)
    call mt_require_real('mt_s1max',mt_s1max)
    call mt_require_real('mt_s2min',mt_s2min)
    call mt_require_real('mt_s2max',mt_s2max)
    call mt_require_positive_int('mt_n1',mt_n1)
    call mt_require_positive_int('mt_n2',mt_n2)
    call mt_require_ordered('mt_s1min','mt_s1max',mt_s1min,mt_s1max)
    call mt_require_ordered('mt_s2min','mt_s2max',mt_s2min,mt_s2max)

    write(*,'(a)') 'mt_run_topology_task: writing arbitrary-plane VTU '// &
         trim(vtu_file)
    if (mt_write_csv) then
      write(*,'(a)') 'mt_run_topology_task: writing arbitrary-plane CSV '// &
           trim(csv_file)
    endif

    if (mt_b_min>0.d0) then
      call mt_fieldline_products_plane_arbitrary(mt_origin(1:ndim), &
           mt_e1(1:ndim),mt_e2(1:ndim),mt_s1min,mt_s1max,mt_n1, &
           mt_s2min,mt_s2max,mt_n2,mt_dL,mt_max_steps,trim(csv_file), &
           b_min=mt_b_min,compute_length=mt_compute_length, &
           compute_twist=mt_compute_twist, &
           compute_q=mt_compute_q,compute_qperp=mt_compute_qperp, &
           vtu_file=trim(vtu_file),write_csv=mt_write_csv)
    else
      call mt_fieldline_products_plane_arbitrary(mt_origin(1:ndim), &
           mt_e1(1:ndim),mt_e2(1:ndim),mt_s1min,mt_s1max,mt_n1, &
           mt_s2min,mt_s2max,mt_n2,mt_dL,mt_max_steps,trim(csv_file), &
           compute_length=mt_compute_length, &
           compute_twist=mt_compute_twist, &
           compute_q=mt_compute_q,compute_qperp=mt_compute_qperp, &
           vtu_file=trim(vtu_file),write_csv=mt_write_csv)
    endif
  end subroutine mt_run_arbitrary_plane_products_task

  subroutine mt_run_seed_products_task()
    character(len=mt_task_name_len) :: csv_file,vtu_file,diag_file
    character(len=mt_task_name_len) :: rk2_diag_file
    double precision, allocatable :: seeds(:,:)
    integer :: nseed

    if (len_trim(mt_output_prefix)==0) then
      call mpistop('seed_products requires mt_output_prefix')
    endif
    if (len_trim(mt_output_file)>0) then
      write(*,'(a)') 'mt_run_topology_task: seed_products ignores mt_output_file'
      write(*,'(a)') 'mt_run_topology_task: seed_products uses mt_output_prefix for CSV and VTU'
    endif
    if (len_trim(mt_seed_file)==0) then
      call mpistop('seed_products requires mt_seed_file')
    endif

    csv_file=''
    if (mt_write_csv) then
      call mt_resolve_prefix_file('_seed_products.csv',csv_file, &
           'seed_products CSV requires mt_output_prefix')
    endif
    call mt_resolve_prefix_file('_seed_products.vtu',vtu_file, &
         'seed_products VTU requires mt_output_prefix')
    if (mt_rk45_tangent_diagnostic) then
      call mt_resolve_prefix_file('_rk45_tangent_diag.csv',diag_file, &
           'seed_products RK45 tangent diagnostic requires mt_output_prefix')
    endif
    if (mt_rk2_fusion_diagnostic) then
      call mt_resolve_prefix_file('_rk2_fusion_diag.csv',rk2_diag_file, &
           'seed_products RK2 fusion diagnostic requires mt_output_prefix')
    endif
    call mt_read_seed_file(trim(mt_seed_file),seeds,nseed)

    if (mt_write_csv) write(*,'(a)') &
         'mt_run_topology_task: writing seed-products CSV '//trim(csv_file)
    write(*,'(a)') 'mt_run_topology_task: writing seed-products VTU '// &
         trim(vtu_file)

    if (mt_b_min>0.d0) then
      call mt_fieldline_products_seeds(seeds,nseed,mt_dL,mt_max_steps, &
           trim(csv_file),b_min=mt_b_min,compute_length=mt_compute_length, &
           compute_twist=mt_compute_twist, &
           compute_q=mt_compute_q,compute_qperp=mt_compute_qperp, &
           vtu_file=trim(vtu_file))
      if (mt_rk45_tangent_diagnostic) then
        write(*,'(a)') 'mt_run_topology_task: writing RK45 tangent '// &
             'diagnostic CSV '//trim(diag_file)
        call mt_rk45_tangent_diagnostic_seeds(seeds,nseed,mt_dL, &
             mt_max_steps,trim(diag_file),b_min=mt_b_min)
      endif
      if (mt_rk2_fusion_diagnostic) then
        write(*,'(a)') 'mt_run_topology_task: writing RK2 fusion '// &
             'diagnostic CSV '//trim(rk2_diag_file)
        call mt_rk2_fusion_diagnostic_seeds(seeds,nseed,mt_dL, &
             mt_max_steps,trim(rk2_diag_file),b_min=mt_b_min)
      endif
    else
      call mt_fieldline_products_seeds(seeds,nseed,mt_dL,mt_max_steps, &
           trim(csv_file),compute_length=mt_compute_length, &
           compute_twist=mt_compute_twist, &
           compute_q=mt_compute_q,compute_qperp=mt_compute_qperp, &
           vtu_file=trim(vtu_file))
      if (mt_rk45_tangent_diagnostic) then
        write(*,'(a)') 'mt_run_topology_task: writing RK45 tangent '// &
             'diagnostic CSV '//trim(diag_file)
        call mt_rk45_tangent_diagnostic_seeds(seeds,nseed,mt_dL, &
             mt_max_steps,trim(diag_file))
      endif
      if (mt_rk2_fusion_diagnostic) then
        write(*,'(a)') 'mt_run_topology_task: writing RK2 fusion '// &
             'diagnostic CSV '//trim(rk2_diag_file)
        call mt_rk2_fusion_diagnostic_seeds(seeds,nseed,mt_dL, &
             mt_max_steps,trim(rk2_diag_file))
      endif
    endif

    deallocate(seeds)
  end subroutine mt_run_seed_products_task

  subroutine mt_run_spherical_surface_products_task()
    character(len=mt_task_name_len) :: surface,csv_file,vtu_file,suffix
    double precision :: seed_coord,s1_min,s1_max,s2_min,s2_max
    double precision :: seed_theta0,seed_phi0,seed_alpha

    surface=mt_lowercase(trim(mt_seed_surface))
    if (len_trim(surface)==0) surface='rmin'
    select case (trim(surface))
    case ('rmin')
      seed_coord=xprobmin1
    case ('rconst','r_const')
      surface='rconst'
      call mt_require_real('mt_seed_coord',mt_seed_coord)
      seed_coord=mt_seed_coord
      if (seed_coord<xprobmin1 .or. seed_coord>xprobmax1) then
        call mpistop('spherical_surface_products requires rconst '// &
             'mt_seed_coord inside domain')
      endif
    {^IFTHREED
    case ('theta_const','thetaconst')
      surface='theta_const'
      call mt_require_real('mt_seed_coord',mt_seed_coord)
      seed_coord=mt_seed_coord
      if (seed_coord<xprobmin2 .or. seed_coord>xprobmax2) then
        call mpistop('spherical_surface_products requires theta_const '// &
             'mt_seed_coord inside domain')
      endif
    case ('phi_const','phiconst')
      surface='phi_const'
      call mt_require_real('mt_seed_coord',mt_seed_coord)
      seed_coord=mt_seed_coord
      if (seed_coord<xprobmin3 .or. seed_coord>xprobmax3) then
        call mpistop('spherical_surface_products requires phi_const '// &
             'mt_seed_coord inside domain')
      endif
    case ('radial_plane','radialplane')
      surface='radial_plane'
      seed_coord=0.d0
      call mt_require_real('mt_seed_theta0',mt_seed_theta0)
      call mt_require_real('mt_seed_phi0',mt_seed_phi0)
      call mt_require_real('mt_seed_alpha',mt_seed_alpha)
      if (mt_seed_theta0<xprobmin2 .or. mt_seed_theta0>xprobmax2) then
        call mpistop('spherical_surface_products requires radial_plane '// &
             'mt_seed_theta0 inside domain')
      endif
      if (mt_seed_phi0<xprobmin3 .or. mt_seed_phi0>xprobmax3) then
        call mpistop('spherical_surface_products requires radial_plane '// &
             'mt_seed_phi0 inside domain')
      endif
      if (abs(dsin(mt_seed_theta0))<=1.d-12) then
        call mpistop('spherical_surface_products radial_plane requires '// &
             'mt_seed_theta0 away from the polar singularity')
      endif
    }
    case default
      call mpistop('spherical_surface_products supports '// &
           'mt_seed_surface=rmin, rconst, theta_const, phi_const, or radial_plane')
    end select
    call mt_require_positive_int('mt_n1',mt_n1)
    call mt_require_positive_int('mt_n2',mt_n2)
    s1_min=0.d0
    s1_max=0.d0
    s2_min=0.d0
    s2_max=0.d0
    {^IFTHREED
    select case (trim(surface))
    case ('rmin','rconst')
      s1_min=xprobmin2
      s1_max=xprobmax2
      s2_min=xprobmin3
      s2_max=xprobmax3
    case ('theta_const')
      s1_min=xprobmin1
      s1_max=xprobmax1
      s2_min=xprobmin3
      s2_max=xprobmax3
    case ('phi_const')
      s1_min=xprobmin1
      s1_max=xprobmax1
      s2_min=xprobmin2
      s2_max=xprobmax2
    case ('radial_plane')
      s1_min=xprobmin1
      s1_max=xprobmax1
      s2_min=-0.5d0*min(xprobmax2-xprobmin2,xprobmax3-xprobmin3)
      s2_max= 0.5d0*min(xprobmax2-xprobmin2,xprobmax3-xprobmin3)
    end select
    }
    if (mt_is_set_real(mt_s1min)) s1_min=mt_s1min
    if (mt_is_set_real(mt_s1max)) s1_max=mt_s1max
    if (mt_is_set_real(mt_s2min)) s2_min=mt_s2min
    if (mt_is_set_real(mt_s2max)) s2_max=mt_s2max
    call mt_require_ordered('mt_s1min','mt_s1max',s1_min,s1_max)
    call mt_require_ordered('mt_s2min','mt_s2max',s2_min,s2_max)

    csv_file=''
    if (mt_write_csv) then
      if (trim(surface)=='rmin') then
        suffix='_spherical_rmin_products.csv'
      else if (trim(surface)=='radial_plane') then
        suffix='_spherical_radial_plane_products.csv'
      else
        suffix='_spherical_'//trim(surface)//'_products.csv'
      endif
      call mt_resolve_prefix_file(trim(suffix),csv_file, &
           'spherical_surface_products CSV requires mt_output_prefix')
    endif
    if (trim(surface)=='rmin') then
      suffix='_spherical_rmin_products.vtu'
    else if (trim(surface)=='radial_plane') then
      suffix='_spherical_radial_plane_products.vtu'
    else
      suffix='_spherical_'//trim(surface)//'_products.vtu'
    endif
    call mt_resolve_output_file('spherical_surface_products','.vtu', &
         trim(suffix),vtu_file)

    if (mt_write_csv) then
      write(*,'(a)') 'mt_run_topology_task: writing spherical '// &
           trim(surface)//' CSV '//trim(csv_file)
    endif
    write(*,'(a)') 'mt_run_topology_task: writing spherical '// &
         trim(surface)//' VTU '//trim(vtu_file)

    if (mt_profile_spherical) then
      call trace_spherical_profile_reset()
      call trace_spherical_profile_set(.true.)
    endif

    seed_theta0=mt_seed_theta0
    seed_phi0=mt_seed_phi0
    seed_alpha=mt_seed_alpha
    if (mt_b_min>0.d0) then
      call mt_fieldline_products_spherical_surface(trim(surface),seed_coord, &
           s1_min,s1_max,mt_n1,s2_min,s2_max,mt_n2,mt_dL,mt_max_steps, &
           trim(csv_file),trim(vtu_file),trim(mt_seed_layout), &
           seed_theta0,seed_phi0,seed_alpha, &
           compute_twist=mt_compute_twist,compute_q=mt_compute_q, &
           compute_qperp=mt_compute_qperp,b_min=mt_b_min)
    else
      call mt_fieldline_products_spherical_surface(trim(surface),seed_coord, &
           s1_min,s1_max,mt_n1,s2_min,s2_max,mt_n2,mt_dL,mt_max_steps, &
           trim(csv_file),trim(vtu_file),trim(mt_seed_layout), &
           seed_theta0,seed_phi0,seed_alpha, &
           compute_twist=mt_compute_twist,compute_q=mt_compute_q, &
           compute_qperp=mt_compute_qperp)
    endif
    if (mt_profile_spherical) then
      call trace_spherical_profile_report('spherical_surface_products '// &
           trim(surface))
      call trace_spherical_profile_set(.false.)
    endif
  end subroutine mt_run_spherical_surface_products_task

  subroutine mt_run_spherical_cloud_products_task()
    character(len=mt_task_name_len) :: vtu_file
    double precision :: s1_min,s1_max,s2_min,s2_max,s3_min,s3_max

    if (mt_write_csv) then
      write(*,'(a)') 'mt_run_topology_task: spherical_cloud_products '// &
           'does not write CSV; use seed_products for selected-point diagnostics'
    endif
    call mt_require_positive_int('mt_n1',mt_n1)
    call mt_require_positive_int('mt_n2',mt_n2)
    call mt_require_positive_int('mt_n3',mt_n3)
    s1_min=xprobmin1
    s1_max=xprobmax1
    s2_min=0.d0
    s2_max=0.d0
    s3_min=0.d0
    s3_max=0.d0
    {^IFTHREED
    s2_min=xprobmin2
    s2_max=xprobmax2
    s3_min=xprobmin3
    s3_max=xprobmax3
    }
    if (mt_is_set_real(mt_s1min)) s1_min=mt_s1min
    if (mt_is_set_real(mt_s1max)) s1_max=mt_s1max
    if (mt_is_set_real(mt_s2min)) s2_min=mt_s2min
    if (mt_is_set_real(mt_s2max)) s2_max=mt_s2max
    if (mt_is_set_real(mt_s3min)) s3_min=mt_s3min
    if (mt_is_set_real(mt_s3max)) s3_max=mt_s3max
    call mt_require_ordered('mt_s1min','mt_s1max',s1_min,s1_max)
    call mt_require_ordered('mt_s2min','mt_s2max',s2_min,s2_max)
    call mt_require_ordered('mt_s3min','mt_s3max',s3_min,s3_max)

    call mt_resolve_output_file('spherical_cloud_products','.vtu', &
         '_spherical_cloud_products.vtu',vtu_file)

    write(*,'(a)') 'mt_run_topology_task: writing spherical cloud VTU '// &
         trim(vtu_file)

    if (mt_profile_spherical) then
      call trace_spherical_profile_reset()
      call trace_spherical_profile_set(.true.)
    endif

    if (mt_b_min>0.d0) then
      call mt_fieldline_products_spherical_cloud(s1_min,s1_max,mt_n1, &
           s2_min,s2_max,mt_n2,s3_min,s3_max,mt_n3,mt_dL,mt_max_steps, &
           trim(vtu_file),trim(mt_seed_layout), &
           compute_twist=mt_compute_twist,compute_q=mt_compute_q, &
           compute_qperp=mt_compute_qperp,b_min=mt_b_min)
    else
      call mt_fieldline_products_spherical_cloud(s1_min,s1_max,mt_n1, &
           s2_min,s2_max,mt_n2,s3_min,s3_max,mt_n3,mt_dL,mt_max_steps, &
           trim(vtu_file),trim(mt_seed_layout), &
           compute_twist=mt_compute_twist,compute_q=mt_compute_q, &
           compute_qperp=mt_compute_qperp)
    endif
    if (mt_profile_spherical) then
      call trace_spherical_profile_report('spherical_cloud_products')
      call trace_spherical_profile_set(.false.)
    endif
  end subroutine mt_run_spherical_cloud_products_task

  subroutine mt_resolve_output_file(mode,extension,suffix,output_file)
    character(len=*), intent(in) :: mode,extension,suffix
    character(len=*), intent(out) :: output_file
    integer :: required_len

    output_file=''
    if (len_trim(mt_output_file)>0) then
      output_file=trim(mt_output_file)
      if (len_trim(mt_output_prefix)>0) then
        write(*,'(a)') 'mt_run_topology_task: mt_output_file takes precedence over mt_output_prefix'
      endif
      call mt_warn_extension(output_file,extension,mode)
    else
      required_len=len_trim(mt_output_prefix)+len_trim(suffix)
      if (required_len>len(output_file)) then
        call mpistop('mt_run_topology_task output filename exceeds internal length')
      endif
      output_file=trim(mt_output_prefix)//trim(suffix)
    endif
  end subroutine mt_resolve_output_file

  subroutine mt_resolve_prefix_file(suffix,output_file,missing_message)
    character(len=*), intent(in) :: suffix,missing_message
    character(len=*), intent(out) :: output_file
    integer :: required_len

    if (len_trim(mt_output_prefix)==0) then
      call mpistop(trim(missing_message))
    endif
    required_len=len_trim(mt_output_prefix)+len_trim(suffix)
    if (required_len>len(output_file)) then
      call mpistop('mt_run_topology_task output filename exceeds internal length')
    endif
    output_file=trim(mt_output_prefix)//trim(suffix)
  end subroutine mt_resolve_prefix_file

  logical function mt_is_set_real(value) result(is_set)
    double precision, intent(in) :: value

    is_set=(value/=mt_unset_real)
  end function mt_is_set_real

  subroutine mt_require_real(name,value)
    character(len=*), intent(in) :: name
    double precision, intent(in) :: value

    if (.not.mt_is_set_real(value)) then
      call mpistop('mt_run_topology_task requires '//trim(name))
    endif
  end subroutine mt_require_real

  subroutine mt_require_positive_int(name,value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: value

    if (value<=0) then
      call mpistop('mt_run_topology_task requires '//trim(name)//' > 0')
    endif
  end subroutine mt_require_positive_int

  subroutine mt_require_ordered(name_min,name_max,value_min,value_max)
    character(len=*), intent(in) :: name_min,name_max
    double precision, intent(in) :: value_min,value_max

    if (value_max<value_min) then
      call mpistop('mt_run_topology_task requires ordered '// &
           trim(name_min)//'/'//trim(name_max))
    endif
  end subroutine mt_require_ordered

  subroutine mt_require_requested_science(do_length,do_twist,do_q,do_qperp, &
       caller)
    logical, intent(in) :: do_length,do_twist,do_q,do_qperp
    character(len=*), intent(in) :: caller

    if (.not.(do_length .or. do_twist .or. do_q .or. do_qperp)) then
      call mpistop(trim(caller)//' requires at least one requested science array')
    endif
  end subroutine mt_require_requested_science

  subroutine mt_read_seed_file(filename,seeds,nseed)
    character(len=*), intent(in) :: filename
    double precision, allocatable, intent(out) :: seeds(:,:)
    integer, intent(out) :: nseed

    character(len=1024) :: line
    double precision :: xyz(3)
    integer :: io,ios,line_no,iseed

    nseed=0
    open(newunit=io,file=trim(filename),status='old',action='read', &
         form='formatted',iostat=ios)
    if (ios/=0) then
      call mpistop('mt_read_seed_file could not open '//trim(filename))
    endif
    line_no=0
    do
      read(io,'(a)',iostat=ios) line
      if (ios/=0) exit
      line_no=line_no+1
      if (.not.mt_seed_line_is_data(line)) cycle
      call mt_validate_seed_line(filename,line,line_no)
      nseed=nseed+1
    enddo
    close(io)
    if (nseed<=0) then
      call mpistop('mt_read_seed_file found no seeds in '//trim(filename))
    endif

    allocate(seeds(nseed,ndim))
    seeds=0.d0
    open(newunit=io,file=trim(filename),status='old',action='read', &
         form='formatted',iostat=ios)
    if (ios/=0) then
      call mpistop('mt_read_seed_file could not reopen '//trim(filename))
    endif
    line_no=0
    iseed=0
    do
      read(io,'(a)',iostat=ios) line
      if (ios/=0) exit
      line_no=line_no+1
      if (.not.mt_seed_line_is_data(line)) cycle
      iseed=iseed+1
      read(line,*,iostat=ios) xyz
      if (ios/=0) then
        call mt_fail_seed_line(filename,line_no)
      endif
      seeds(iseed,1:ndim)=xyz(1:ndim)
    enddo
    close(io)
  end subroutine mt_read_seed_file

  logical function mt_seed_line_is_data(line) result(is_data)
    character(len=*), intent(in) :: line
    character(len=len(line)) :: trimmed

    trimmed=adjustl(line)
    is_data=len_trim(trimmed)>0
    if (is_data) is_data=trimmed(1:1)/='#'
  end function mt_seed_line_is_data

  subroutine mt_validate_seed_line(filename,line,line_no)
    character(len=*), intent(in) :: filename,line
    integer, intent(in) :: line_no

    double precision :: x,y,z
    integer :: ios

    if (index(line,',')>0) then
      call mt_fail_seed_line(filename,line_no)
    endif
    if (mt_count_tokens(line)/=3) then
      call mt_fail_seed_line(filename,line_no)
    endif
    read(line,*,iostat=ios) x,y,z
    if (ios/=0) then
      call mt_fail_seed_line(filename,line_no)
    endif
  end subroutine mt_validate_seed_line

  integer function mt_count_tokens(line) result(ntoken)
    character(len=*), intent(in) :: line
    integer :: i
    logical :: in_token

    ntoken=0
    in_token=.false.
    do i=1,len_trim(line)
      if (line(i:i)==' ' .or. line(i:i)==achar(9)) then
        in_token=.false.
      else
        if (.not.in_token) then
          ntoken=ntoken+1
          in_token=.true.
        endif
      endif
    enddo
  end function mt_count_tokens

  subroutine mt_fail_seed_line(filename,line_no)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: line_no

    call mpistop('mt_read_seed_file malformed line '// &
         trim(mt_int_to_string(line_no))//' in '//trim(filename))
  end subroutine mt_fail_seed_line

  function mt_int_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text

    write(text,'(i0)') value
  end function mt_int_to_string

  subroutine mt_warn_extension(filename,extension,mode)
    character(len=*), intent(in) :: filename,extension,mode

    if (.not.mt_has_extension(filename,extension)) then
      write(*,'(a)') 'mt_run_topology_task warning: '//trim(mode)// &
           ' output file does not end in '//trim(extension)
    endif
  end subroutine mt_warn_extension

  logical function mt_has_extension(filename,extension) result(has_ext)
    character(len=*), intent(in) :: filename,extension
    character(len=mt_task_name_len) :: fname,ext
    integer :: lf,le

    fname=mt_lowercase(trim(filename))
    ext=mt_lowercase(trim(extension))
    lf=len_trim(fname)
    le=len_trim(ext)
    has_ext=.false.
    if (lf>=le) has_ext=(fname(lf-le+1:lf)==ext(1:le))
  end function mt_has_extension

  function mt_lowercase(input) result(output)
    character(len=*), intent(in) :: input
    character(len=len(input)) :: output
    integer :: i,code

    output=input
    do i=1,len(input)
      code=iachar(output(i:i))
      if (code>=iachar('A') .and. code<=iachar('Z')) then
        output(i:i)=achar(code+iachar('a')-iachar('A'))
      endif
    enddo
  end function mt_lowercase

  logical function mt_vtk_detail_is_full() result(is_full)
    is_full=(mt_lowercase(trim(mt_vtk_detail))=='full')
  end function mt_vtk_detail_is_full

  double precision function mt_visual_float(value) result(out_value)
    double precision, intent(in) :: value

    if (ieee_is_finite(value)) then
      out_value=value
    else
      out_value=0.d0
    endif
  end function mt_visual_float

  double precision function mt_visual_valid_float(value,is_valid) &
       result(out_value)
    double precision, intent(in) :: value
    logical, intent(in) :: is_valid

    if (is_valid .and. ieee_is_finite(value)) then
      out_value=value
    else
      out_value=0.d0
    endif
  end function mt_visual_valid_float

  double precision function mt_vc(vec,icomp) result(value)
    ! Return a physical coordinate component without creating fixed
    ! out-of-bounds references in 1D/2D source expansions.
    double precision, intent(in) :: vec(ndim)
    integer, intent(in) :: icomp

    if (icomp>=1 .and. icomp<=ndim) then
      value=vec(icomp)
    else
      value=0.d0
    endif
  end function mt_vc

  subroutine mt_length_single(seed,dL,max_steps,csv_file,b_min)
    ! Trace one magnetic field line and write its summary to a CSV file.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_length_result) :: result
    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status

    if (npe/=1) then
      call mpistop('mt_length_single currently requires npe=1')
    endif

    if (present(b_min)) then
      call trace_field_length_single(seed,dL,max_steps,result,b_min)
    else
      call trace_field_length_single(seed,dL,max_steps,result)
    endif

    seed_xyz=0.d0
    seed_xyz(1:ndim)=result%seed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop('mt_length_single could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_length_single could not write CSV header')
    endif

    write(csv_unit,'(es24.16,5(",",es24.16),4(",",i0))',iostat=io_status) &
         seed_xyz,result%total_length, &
         result%backward_length,result%forward_length, &
         result%backward_nstep,result%forward_nstep, &
         result%backward_status,result%forward_status
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_length_single could not write CSV data')
    endif

    close(csv_unit)
  end subroutine mt_length_single

  subroutine mt_length_seeds(seeds,nseed,dL,max_steps,csv_file,b_min)
    ! Trace multiple magnetic field lines and write one summary row per seed.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_length_result), allocatable :: results(:)
    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status,iseed

    if (npe/=1) then
      call mpistop('mt_length_seeds currently requires npe=1')
    endif
    if (nseed<0) then
      call mpistop('mt_length_seeds requires nseed>=0')
    endif

    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_length_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_length_multi(seeds,nseed,dL,max_steps,results)
    endif

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      deallocate(results)
      call mpistop('mt_length_seeds could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'seed_id,seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      deallocate(results)
      call mpistop('mt_length_seeds could not write CSV header')
    endif

    do iseed=1,nseed
      seed_xyz=0.d0
      seed_xyz(1:ndim)=results(iseed)%seed
      write(csv_unit,'(i0,6(",",es24.16),4(",",i0))',iostat=io_status) &
           iseed,seed_xyz,results(iseed)%total_length, &
           results(iseed)%backward_length,results(iseed)%forward_length, &
           results(iseed)%backward_nstep,results(iseed)%forward_nstep, &
           results(iseed)%backward_status,results(iseed)%forward_status
      if (io_status/=0) then
        close(csv_unit)
        deallocate(results)
        call mpistop('mt_length_seeds could not write CSV data')
      endif
    enddo

    close(csv_unit)
    deallocate(results)
  end subroutine mt_length_seeds

  subroutine mt_fieldline_products_seeds(seeds,nseed,dL,max_steps,csv_file, &
       b_min,compute_length,compute_twist,compute_q,compute_qperp,vtu_file)
    ! Write selected per-field-line diagnostics for an arbitrary seed set.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min
    logical, intent(in), optional :: compute_length,compute_twist
    logical, intent(in), optional :: compute_q,compute_qperp
    character(len=*), intent(in), optional :: vtu_file

    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_qperp_result), allocatable :: q_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    logical :: do_length,do_twist,do_q,do_qperp
    character(len=mt_task_name_len) :: integrator

    do_length=.true.
    do_twist=.false.
    do_q=.false.
    do_qperp=.false.
    if (present(compute_length)) do_length=compute_length
    if (present(compute_twist)) do_twist=compute_twist
    if (present(compute_q)) do_q=compute_q
    if (present(compute_qperp)) do_qperp=compute_qperp
    integrator=mt_lowercase(trim(mt_trace_integrator))

    if (npe/=1) then
      call mpistop('mt_fieldline_products_seeds currently requires npe=1')
    endif
    if (ndim/=3) then
      call mpistop('mt_fieldline_products_seeds requires 3D Cartesian or spherical geometry')
    endif
    {^IFTHREED
    select case (geo_coordinate)
    case (geo_cartesian,geo_cartesian_stretched)
    case (geo_spherical)
      if (periodB(3)) then
        call mpistop('spherical seed_products does not yet support periodic phi')
      endif
    case default
      call mpistop('mt_fieldline_products_seeds requires Cartesian or spherical geometry')
    end select
    }
    if (do_qperp) then
      select case (trim(integrator))
      case ('rk2')
        ! RK2 Qperp is the mature tangent-transport path.
      case ('rk45_cartesian')
        {^IFTHREED
        if (geo_coordinate==geo_spherical) then
          call mpistop('seed_products logQperp with spherical geometry '// &
               'requires mt_trace_integrator=rk45_spherical')
        endif
        }
      case ('rk45_spherical')
        {^IFTHREED
        if (geo_coordinate/=geo_spherical) then
          call mpistop('seed_products logQperp with Cartesian geometry '// &
               'requires mt_trace_integrator=rk45_cartesian')
        endif
        }
      case default
        call mpistop('seed_products logQperp requires mt_trace_integrator='// &
             'rk2, rk45_cartesian, or rk45_spherical')
      end select
    endif
    if (do_q) then
      select case (trim(integrator))
      case ('rk2')
        ! RK2 standard logQ is the mature tangent-transport path.
      case ('rk45_cartesian')
        {^IFTHREED
        if (geo_coordinate==geo_spherical) then
          call mpistop('seed_products logQ with spherical geometry '// &
               'requires mt_trace_integrator=rk45_spherical')
        endif
        }
      case ('rk45_spherical')
        {^IFTHREED
        if (geo_coordinate/=geo_spherical) then
          call mpistop('seed_products logQ with Cartesian geometry '// &
               'requires mt_trace_integrator=rk45_cartesian')
        endif
        }
      case default
        call mpistop('seed_products logQ requires mt_trace_integrator='// &
             'rk2, rk45_cartesian, or rk45_spherical')
      end select
    endif
    if (nseed<0) then
      call mpistop('mt_fieldline_products_seeds requires nseed>=0')
    endif
    call mt_require_requested_science(do_length,do_twist,do_q,do_qperp, &
         'mt_fieldline_products_seeds')

    allocate(length_results(nseed))
    if (do_twist) then
      allocate(twist_results(nseed))
    else
      allocate(twist_results(0))
    endif
    if (do_q) then
      allocate(q_results(nseed))
    else
      allocate(q_results(0))
    endif
    if (do_qperp) then
      allocate(qperp_results(nseed))
    else
      allocate(qperp_results(0))
    endif

    if (present(b_min)) then
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,q_results,qperp_results,do_twist, &
           do_q,do_qperp,b_min)
    else
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,q_results,qperp_results,do_twist, &
           do_q,do_qperp)
    endif

    if (len_trim(csv_file)>0) then
      call mt_write_fieldline_products_seeds_csv(length_results, &
           twist_results,q_results,qperp_results,nseed,csv_file,do_twist, &
           do_q,do_qperp)
    endif
    if (present(vtu_file)) then
      if (len_trim(vtu_file)>0) then
        call mt_write_fieldline_products_vtu_vertices(vtu_file, &
             length_results,twist_results,q_results,qperp_results,nseed, &
             do_length,do_twist,do_q,do_qperp,'mt_fieldline_products_seeds')
      endif
    endif

    if (allocated(qperp_results)) deallocate(qperp_results)
    if (allocated(q_results)) deallocate(q_results)
    if (allocated(twist_results)) deallocate(twist_results)
    deallocate(length_results)
  end subroutine mt_fieldline_products_seeds

  subroutine mt_rk45_tangent_diagnostic_seeds(seeds,nseed,dL,max_steps, &
       csv_file,b_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_qperp_result), allocatable :: rk2_results(:)
    type(trace_qperp_result), allocatable :: rk45_results(:)

    if (npe/=1) then
      call mpistop('mt_rk45_tangent_diagnostic_seeds currently requires npe=1')
    endif
    if (nseed<0) then
      call mpistop('mt_rk45_tangent_diagnostic_seeds requires nseed>=0')
    endif

    allocate(rk2_results(nseed),rk45_results(nseed))
    if (present(b_min)) then
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,rk2_results, &
           b_min)
      call trace_debug_cartesian_rk45_tangent_q0_multi(seeds,nseed,dL, &
           max_steps,rk45_results,b_min)
    else
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,rk2_results)
      call trace_debug_cartesian_rk45_tangent_q0_multi(seeds,nseed,dL, &
           max_steps,rk45_results)
    endif

    call mt_write_rk45_tangent_diagnostic_csv(rk2_results,rk45_results, &
         nseed,csv_file)
    deallocate(rk45_results,rk2_results)
  end subroutine mt_rk45_tangent_diagnostic_seeds

  subroutine mt_write_rk45_tangent_diagnostic_csv(rk2_results,rk45_results, &
       nseed,csv_file)
    integer, intent(in) :: nseed
    type(trace_qperp_result), intent(in) :: rk2_results(nseed)
    type(trace_qperp_result), intent(in) :: rk45_results(nseed)
    character(len=*), intent(in) :: csv_file

    double precision :: seed_xyz(3),logq_diff,q_diff,length_diff
    double precision :: f_endpoint_diff,b_endpoint_diff
    double precision :: uf_diff,vf_diff,ub_diff,vb_diff
    integer :: csv_unit,io_status,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop('mt_rk45_tangent_diagnostic could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'seed_id,seed_x,seed_y,seed_z,'// &
         'rk2_status_q0,rk45_status_q0,'// &
         'rk2_forward_status,rk45_forward_status,'// &
         'rk2_backward_status,rk45_backward_status,'// &
         'rk2_logq0,rk45_logq0,abs_diff_logq0,'// &
         'rk2_q0,rk45_q0,abs_diff_q0,'// &
         'rk2_length_total,rk45_length_total,abs_diff_length_total,'// &
         'endpoint_diff_forward,endpoint_diff_backward,'// &
         'u_forward_perp_diff,v_forward_perp_diff,'// &
         'u_backward_perp_diff,v_backward_perp_diff'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_rk45_tangent_diagnostic could not write CSV header')
    endif

    do iseed=1,nseed
      seed_xyz=0.d0
      seed_xyz(1:ndim)=rk2_results(iseed)%seed
      logq_diff=mt_absdiff_or_nan(rk2_results(iseed)%logq0, &
           rk45_results(iseed)%logq0)
      q_diff=mt_absdiff_or_nan(rk2_results(iseed)%q0, &
           rk45_results(iseed)%q0)
      length_diff=abs((rk2_results(iseed)%forward_length+ &
           rk2_results(iseed)%backward_length)- &
           (rk45_results(iseed)%forward_length+ &
           rk45_results(iseed)%backward_length))
      f_endpoint_diff=dsqrt(sum((rk2_results(iseed)%forward_endpoint- &
           rk45_results(iseed)%forward_endpoint)**2))
      b_endpoint_diff=dsqrt(sum((rk2_results(iseed)%backward_endpoint- &
           rk45_results(iseed)%backward_endpoint)**2))
      uf_diff=dsqrt(sum((rk2_results(iseed)%u_forward_perp- &
           rk45_results(iseed)%u_forward_perp)**2))
      vf_diff=dsqrt(sum((rk2_results(iseed)%v_forward_perp- &
           rk45_results(iseed)%v_forward_perp)**2))
      ub_diff=dsqrt(sum((rk2_results(iseed)%u_backward_perp- &
           rk45_results(iseed)%u_backward_perp)**2))
      vb_diff=dsqrt(sum((rk2_results(iseed)%v_backward_perp- &
           rk45_results(iseed)%v_backward_perp)**2))

      write(csv_unit,'(i0,3(",",es24.16),6(",",i0),15(",",es24.16))', &
           iostat=io_status) &
           iseed,seed_xyz, &
           rk2_results(iseed)%status_q0, &
           rk45_results(iseed)%status_q0, &
           rk2_results(iseed)%forward_status, &
           rk45_results(iseed)%forward_status, &
           rk2_results(iseed)%backward_status, &
           rk45_results(iseed)%backward_status, &
           rk2_results(iseed)%logq0,rk45_results(iseed)%logq0,logq_diff, &
           rk2_results(iseed)%q0,rk45_results(iseed)%q0,q_diff, &
           rk2_results(iseed)%forward_length+ &
           rk2_results(iseed)%backward_length, &
           rk45_results(iseed)%forward_length+ &
           rk45_results(iseed)%backward_length, &
           length_diff,f_endpoint_diff,b_endpoint_diff, &
           uf_diff,vf_diff,ub_diff,vb_diff
      if (io_status/=0) then
        close(csv_unit)
        call mpistop('mt_rk45_tangent_diagnostic could not write CSV data')
      endif
    enddo

    close(csv_unit)
  end subroutine mt_write_rk45_tangent_diagnostic_csv

  subroutine mt_rk2_fusion_diagnostic_seeds(seeds,nseed,dL,max_steps, &
       csv_file,b_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_topology_result), allocatable :: summary(:)
    type(trace_qperp_result), allocatable :: q_trace(:)
    type(trace_qperp_result), allocatable :: q_short(:)
    type(trace_twist_result), allocatable :: q_twist(:)
    type(trace_twist_result), allocatable :: q_short_twist(:)
    integer :: cache_status
    logical :: use_spherical_cache

    if (npe/=1) then
      call mpistop('mt_rk2_fusion_diagnostic currently requires npe=1')
    endif
    if (nseed<0) then
      call mpistop('mt_rk2_fusion_diagnostic requires nseed>=0')
    endif
    if (mt_lowercase(trim(mt_trace_integrator))/='rk2') then
      call mpistop('mt_rk2_fusion_diagnostic requires rk2 tracing')
    endif

    allocate(summary(nseed),q_trace(nseed),q_short(nseed), &
         q_twist(nseed),q_short_twist(nseed))
    use_spherical_cache=.false.
    {^IFTHREED
    if (geo_coordinate==geo_spherical) then
      call trace_spherical_curl_cache_build(cache_status)
      if (cache_status/=trace_status_active) then
        call mpistop('mt_rk2_fusion_diagnostic failed to build '// &
             'spherical curl cache')
      endif
      use_spherical_cache=.true.
    endif
    }

    if (present(b_min)) then
      call trace_field_topology_multi(seeds,nseed,dL,max_steps,summary, &
           need_twist=.true.,need_mapping=.false.,b_min=b_min)
      if (geo_coordinate==geo_spherical) then
        call trace_field_spherical_rmin_q_multi(seeds,nseed,dL,max_steps, &
             q_trace,b_min,twist_results=q_twist)
        call trace_field_rk2_short_boundary_q_multi(seeds,nseed,dL, &
             max_steps,q_short,b_min,twist_results=q_short_twist)
      else
        call trace_field_qperp_multi(seeds,nseed,dL,max_steps,q_trace, &
             b_min,twist_results=q_twist)
        call trace_field_rk2_short_boundary_q_multi(seeds,nseed,dL, &
             max_steps,q_short,b_min,twist_results=q_short_twist)
      endif
    else
      call trace_field_topology_multi(seeds,nseed,dL,max_steps,summary, &
           need_twist=.true.,need_mapping=.false.)
      if (geo_coordinate==geo_spherical) then
        call trace_field_spherical_rmin_q_multi(seeds,nseed,dL,max_steps, &
             q_trace,twist_results=q_twist)
        call trace_field_rk2_short_boundary_q_multi(seeds,nseed,dL, &
             max_steps,q_short,twist_results=q_short_twist)
      else
        call trace_field_qperp_multi(seeds,nseed,dL,max_steps,q_trace, &
             twist_results=q_twist)
        call trace_field_rk2_short_boundary_q_multi(seeds,nseed,dL, &
             max_steps,q_short,twist_results=q_short_twist)
      endif
    endif

    if (use_spherical_cache) call trace_spherical_curl_cache_clear()
    call mt_write_rk2_fusion_diagnostic_csv(summary,q_trace,q_twist, &
         q_short,q_short_twist,nseed,csv_file)
    deallocate(q_short_twist,q_twist,q_short,q_trace,summary)
  end subroutine mt_rk2_fusion_diagnostic_seeds

  subroutine mt_write_rk2_fusion_diagnostic_csv(summary,q_trace,q_twist, &
       q_short,q_short_twist,nseed,csv_file)
    integer, intent(in) :: nseed
    type(trace_topology_result), intent(in) :: summary(nseed)
    type(trace_qperp_result), intent(in) :: q_trace(nseed)
    type(trace_twist_result), intent(in) :: q_twist(nseed)
    type(trace_qperp_result), intent(in) :: q_short(nseed)
    type(trace_twist_result), intent(in) :: q_short_twist(nseed)
    character(len=*), intent(in) :: csv_file

    double precision :: seed_xyz(3),length_summary,length_q,length_diff
    double precision :: length_short,length_diff_short,length_q_short_diff
    double precision :: twist_summary,twist_q,twist_diff
    double precision :: twist_short,twist_diff_short,twist_q_short_diff
    double precision :: max_length_diff,max_twist_diff
    double precision :: max_length_short_diff,max_twist_short_diff
    integer :: csv_unit,io_status,iseed
    integer :: status_mismatch,face_mismatch,twist_status_mismatch
    integer :: short_status_mismatch,short_face_mismatch
    integer :: short_twist_status_mismatch,short_valid_q_mismatch

    max_length_diff=0.d0
    max_twist_diff=0.d0
    max_length_short_diff=0.d0
    max_twist_short_diff=0.d0
    status_mismatch=0
    face_mismatch=0
    twist_status_mismatch=0
    short_status_mismatch=0
    short_face_mismatch=0
    short_twist_status_mismatch=0
    short_valid_q_mismatch=0

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop('mt_rk2_fusion_diagnostic could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'seed_id,seed_x,seed_y,seed_z,'// &
         'length_summary,length_qtrace,dlength_summary_minus_qtrace,'// &
         'length_forward_summary,length_forward_qtrace,'// &
         'length_backward_summary,length_backward_qtrace,'// &
         'twist_summary,twist_qtrace,dtwist_summary_minus_qtrace,'// &
         'twist_forward_summary,twist_forward_qtrace,'// &
         'twist_backward_summary,twist_backward_qtrace,'// &
         'length_short,dlength_summary_minus_short,'// &
         'dlength_qtrace_minus_short,'// &
         'length_forward_short,length_backward_short,'// &
         'twist_short,dtwist_summary_minus_short,'// &
         'dtwist_qtrace_minus_short,'// &
         'twist_forward_short,twist_backward_short,'// &
         'status_forward_summary,status_forward_qtrace,status_forward_short,'// &
         'status_backward_summary,status_backward_qtrace,status_backward_short,'// &
         'face_forward_summary,face_forward_qtrace,face_forward_short,'// &
         'face_backward_summary,face_backward_qtrace,face_backward_short,'// &
         'status_twist_summary,status_twist_qtrace,status_twist_short,'// &
         'logQ,valid_Q,status_Q,logQ_short,valid_Q_short,status_Q_short'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_rk2_fusion_diagnostic could not write CSV header')
    endif

    do iseed=1,nseed
      seed_xyz=0.d0
      seed_xyz(1:ndim)=summary(iseed)%seed
      length_summary=summary(iseed)%length_total
      length_q=q_trace(iseed)%forward_length+q_trace(iseed)%backward_length
      length_short=q_short(iseed)%forward_length+ &
           q_short(iseed)%backward_length
      length_diff=length_summary-length_q
      length_diff_short=length_summary-length_short
      length_q_short_diff=length_q-length_short
      twist_summary=summary(iseed)%twist_total
      twist_q=q_twist(iseed)%total_twist
      twist_short=q_short_twist(iseed)%total_twist
      twist_diff=twist_summary-twist_q
      twist_diff_short=twist_summary-twist_short
      twist_q_short_diff=twist_q-twist_short
      max_length_diff=max(max_length_diff,abs(length_diff))
      max_twist_diff=max(max_twist_diff,abs(twist_diff))
      max_length_short_diff=max(max_length_short_diff, &
           abs(length_diff_short))
      max_twist_short_diff=max(max_twist_short_diff, &
           abs(twist_diff_short))
      if (summary(iseed)%forward_status/=q_trace(iseed)%forward_status .or. &
           summary(iseed)%backward_status/=q_trace(iseed)%backward_status) &
           status_mismatch=status_mismatch+1
      if (summary(iseed)%forward_face/=q_trace(iseed)%forward_face .or. &
           summary(iseed)%backward_face/=q_trace(iseed)%backward_face) &
           face_mismatch=face_mismatch+1
      if (summary(iseed)%status_twist/=q_twist(iseed)%status_twist) &
           twist_status_mismatch=twist_status_mismatch+1
      if (summary(iseed)%forward_status/=q_short(iseed)%forward_status .or. &
           summary(iseed)%backward_status/=q_short(iseed)%backward_status) &
           short_status_mismatch=short_status_mismatch+1
      if (summary(iseed)%forward_face/=q_short(iseed)%forward_face .or. &
           summary(iseed)%backward_face/=q_short(iseed)%backward_face) &
           short_face_mismatch=short_face_mismatch+1
      if (summary(iseed)%status_twist/=q_short_twist(iseed)%status_twist) &
           short_twist_status_mismatch=short_twist_status_mismatch+1
      if (q_trace(iseed)%valid_q0 .neqv. q_short(iseed)%valid_q0) &
           short_valid_q_mismatch=short_valid_q_mismatch+1

      write(csv_unit, &
           '(i0,3(",",es24.16),24(",",es24.16),15(",",i0),'// &
           '",",es24.16,",",l1,",",i0,",",es24.16,",",l1,",",i0)', &
           iostat=io_status) &
           iseed,seed_xyz, &
           length_summary,length_q,length_diff, &
           summary(iseed)%length_forward,q_trace(iseed)%forward_length, &
           summary(iseed)%length_backward,q_trace(iseed)%backward_length, &
           twist_summary,twist_q,twist_diff, &
           summary(iseed)%twist_forward,q_twist(iseed)%forward_twist, &
           summary(iseed)%twist_backward,q_twist(iseed)%backward_twist, &
           length_short,length_diff_short,length_q_short_diff, &
           q_short(iseed)%forward_length,q_short(iseed)%backward_length, &
           twist_short,twist_diff_short,twist_q_short_diff, &
           q_short_twist(iseed)%forward_twist, &
           q_short_twist(iseed)%backward_twist, &
           summary(iseed)%forward_status,q_trace(iseed)%forward_status, &
           q_short(iseed)%forward_status, &
           summary(iseed)%backward_status,q_trace(iseed)%backward_status, &
           q_short(iseed)%backward_status, &
           summary(iseed)%forward_face,q_trace(iseed)%forward_face, &
           q_short(iseed)%forward_face, &
           summary(iseed)%backward_face,q_trace(iseed)%backward_face, &
           q_short(iseed)%backward_face, &
           summary(iseed)%status_twist,q_twist(iseed)%status_twist, &
           q_short_twist(iseed)%status_twist, &
           q_trace(iseed)%logq0,q_trace(iseed)%valid_q0, &
           q_trace(iseed)%status_q0, &
           q_short(iseed)%logq0,q_short(iseed)%valid_q0, &
           q_short(iseed)%status_q0
      if (io_status/=0) then
        close(csv_unit)
        call mpistop('mt_rk2_fusion_diagnostic could not write CSV data')
      endif
    enddo

    close(csv_unit)
    write(*,'(a,es12.4)') 'mt_rk2_fusion_diagnostic max_abs_length_diff: ', &
         max_length_diff
    write(*,'(a,es12.4)') 'mt_rk2_fusion_diagnostic max_abs_twist_diff: ', &
         max_twist_diff
    write(*,'(a,es12.4)') &
         'mt_rk2_fusion_diagnostic max_abs_length_diff_short: ', &
         max_length_short_diff
    write(*,'(a,es12.4)') &
         'mt_rk2_fusion_diagnostic max_abs_twist_diff_short: ', &
         max_twist_short_diff
    write(*,'(a,i0)') 'mt_rk2_fusion_diagnostic status_mismatch_count: ', &
         status_mismatch
    write(*,'(a,i0)') 'mt_rk2_fusion_diagnostic face_mismatch_count: ', &
         face_mismatch
    write(*,'(a,i0)') &
         'mt_rk2_fusion_diagnostic twist_status_mismatch_count: ', &
         twist_status_mismatch
    write(*,'(a,i0)') &
         'mt_rk2_fusion_diagnostic short_status_mismatch_count: ', &
         short_status_mismatch
    write(*,'(a,i0)') &
         'mt_rk2_fusion_diagnostic short_face_mismatch_count: ', &
         short_face_mismatch
    write(*,'(a,i0)') &
         'mt_rk2_fusion_diagnostic short_twist_status_mismatch_count: ', &
         short_twist_status_mismatch
    write(*,'(a,i0)') &
         'mt_rk2_fusion_diagnostic short_valid_q_mismatch_count: ', &
         short_valid_q_mismatch
  end subroutine mt_write_rk2_fusion_diagnostic_csv

  double precision function mt_absdiff_or_nan(a,b) result(diff)
    double precision, intent(in) :: a,b

    if (ieee_is_finite(a) .and. ieee_is_finite(b)) then
      diff=abs(a-b)
    else
      diff=ieee_value(0.d0,ieee_quiet_nan)
    endif
  end function mt_absdiff_or_nan

  subroutine mt_fieldline_products_spherical_surface(surface,seed_coord, &
       s1_min,s1_max,n1,s2_min,s2_max,n2,dL,max_steps,csv_file, &
       vtu_file,seed_layout,seed_theta0,seed_phi0,seed_alpha, &
       compute_twist,compute_q,compute_qperp,b_min)
    character(len=*), intent(in) :: surface,csv_file,vtu_file,seed_layout
    integer, intent(in) :: n1,n2,max_steps
    double precision, intent(in) :: seed_coord,s1_min,s1_max,s2_min,s2_max,dL
    double precision, intent(in) :: seed_theta0,seed_phi0,seed_alpha
    logical, intent(in), optional :: compute_twist,compute_q,compute_qperp
    double precision, intent(in), optional :: b_min

    type(trace_topology_result), allocatable :: topology(:)
    type(trace_qperp_result), allocatable :: q_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed
    logical :: do_twist,do_q,do_qperp

    if (npe/=1) then
      call mpistop('mt_fieldline_products_spherical_surface currently '// &
           'requires npe=1')
    endif
    if (ndim/=3 .or. geo_coordinate/=geo_spherical) then
      call mpistop('mt_fieldline_products_spherical_surface requires '// &
           '3D spherical geometry')
    endif
    {^IFTHREED
    if (periodB(3)) then
      call mpistop('mt_fieldline_products_spherical_surface does not '// &
           'yet support periodic phi')
    endif
    }
    if (n1<1 .or. n2<1) then
      call mpistop('mt_fieldline_products_spherical_surface requires '// &
           'sample counts >=1')
    endif
    if (s1_max<s1_min .or. s2_max<s2_min) then
      call mpistop('mt_fieldline_products_spherical_surface requires ordered bounds')
    endif
    do_twist=.false.
    do_q=.false.
    do_qperp=.false.
    if (present(compute_twist)) do_twist=compute_twist
    if (present(compute_q)) do_q=compute_q
    if (present(compute_qperp)) do_qperp=compute_qperp
    call mt_build_spherical_surface_seeds(surface,seed_coord,s1_min,s1_max, &
         n1,s2_min,s2_max,n2,seed_layout,seed_theta0,seed_phi0, &
         seed_alpha,seeds)
    nseed=n1*n2
    if (mt_profile_spherical) call trace_spherical_profile_count_seeds(nseed)
    allocate(topology(nseed))
    if (do_q) then
      allocate(q_results(nseed))
    else
      allocate(q_results(0))
    endif
    if (do_qperp) then
      allocate(qperp_results(nseed))
    else
      allocate(qperp_results(0))
    endif

    if (present(b_min)) then
      call mt_trace_spherical_surface_products(seeds,nseed,dL,max_steps, &
           topology,q_results,qperp_results,do_twist,do_q,do_qperp,b_min)
    else
      call mt_trace_spherical_surface_products(seeds,nseed,dL,max_steps, &
           topology,q_results,qperp_results,do_twist,do_q,do_qperp)
    endif

    if (len_trim(csv_file)>0) then
      call mt_write_spherical_rmin_csv(topology,n1,n2,csv_file, &
           'mt_fieldline_products_spherical_surface',do_twist,do_q, &
           q_results,do_qperp,qperp_results)
    endif
    call mt_write_spherical_topology_vtu(vtu_file,topology,n1,n2, &
         'mt_fieldline_products_spherical_surface',do_twist,do_q,q_results, &
         do_qperp,qperp_results)

    deallocate(qperp_results,q_results,topology,seeds)
  end subroutine mt_fieldline_products_spherical_surface

  subroutine mt_fieldline_products_spherical_cloud(s1_min,s1_max,n1, &
       s2_min,s2_max,n2,s3_min,s3_max,n3,dL,max_steps,vtu_file, &
       seed_layout,compute_twist,compute_q,compute_qperp,b_min)
    character(len=*), intent(in) :: vtu_file,seed_layout
    integer, intent(in) :: n1,n2,n3,max_steps
    double precision, intent(in) :: s1_min,s1_max,s2_min,s2_max
    double precision, intent(in) :: s3_min,s3_max,dL
    logical, intent(in), optional :: compute_twist,compute_q,compute_qperp
    double precision, intent(in), optional :: b_min

    type(trace_topology_result), allocatable :: topology(:)
    type(trace_qperp_result), allocatable :: q_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed
    logical :: do_twist,do_q,do_qperp
    character(len=mt_task_name_len) :: integrator

    if (npe/=1) then
      call mpistop('mt_fieldline_products_spherical_cloud currently '// &
           'requires npe=1')
    endif
    if (ndim/=3 .or. geo_coordinate/=geo_spherical) then
      call mpistop('mt_fieldline_products_spherical_cloud requires '// &
           '3D spherical geometry')
    endif
    {^IFTHREED
    if (periodB(3)) then
      call mpistop('mt_fieldline_products_spherical_cloud does not '// &
           'yet support periodic phi')
    endif
    }
    if (n1<1 .or. n2<1 .or. n3<1) then
      call mpistop('mt_fieldline_products_spherical_cloud requires '// &
           'sample counts >=1')
    endif
    if (s1_max<s1_min .or. s2_max<s2_min .or. s3_max<s3_min) then
      call mpistop('mt_fieldline_products_spherical_cloud requires ordered bounds')
    endif

    do_twist=.false.
    do_q=.false.
    do_qperp=.false.
    if (present(compute_twist)) do_twist=compute_twist
    if (present(compute_q)) do_q=compute_q
    if (present(compute_qperp)) do_qperp=compute_qperp
    integrator=mt_lowercase(trim(mt_trace_integrator))

    call mt_build_spherical_cloud_seeds(s1_min,s1_max,n1,s2_min,s2_max, &
         n2,s3_min,s3_max,n3,seed_layout,seeds)
    nseed=n1*n2*n3
    if (mt_profile_spherical) call trace_spherical_profile_count_seeds(nseed)
    allocate(topology(nseed))
    if (do_q) then
      allocate(q_results(nseed))
    else
      allocate(q_results(0))
    endif
    if (do_qperp) then
      allocate(qperp_results(nseed))
    else
      allocate(qperp_results(0))
    endif

    if (present(b_min)) then
      call mt_trace_spherical_cloud_products(seeds,nseed,dL,max_steps, &
           topology,q_results,qperp_results,do_twist,do_q,do_qperp,b_min)
    else
      call mt_trace_spherical_cloud_products(seeds,nseed,dL,max_steps, &
           topology,q_results,qperp_results,do_twist,do_q,do_qperp)
    endif

    call mt_write_spherical_cloud_vtu(vtu_file,topology,nseed, &
         'mt_fieldline_products_spherical_cloud',do_twist,do_q,q_results, &
         do_qperp,qperp_results)

    deallocate(qperp_results,q_results,topology,seeds)
  end subroutine mt_fieldline_products_spherical_cloud

  subroutine mt_trace_spherical_surface_products(seeds,nseed,dL,max_steps, &
       topology,q_results,qperp_results,do_twist,do_q,do_qperp,b_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_topology_result), intent(out) :: topology(nseed)
    type(trace_qperp_result), intent(out) :: q_results(:)
    type(trace_qperp_result), intent(out) :: qperp_results(:)
    logical, intent(in) :: do_twist,do_q,do_qperp
    double precision, intent(in), optional :: b_min

    type(trace_qperp_result) :: q_local,qperp_local
    double precision :: seed_local(ndim)
    integer :: iseed,cache_status

    if (nseed<=0) return

    if (do_twist) then
      call trace_spherical_curl_cache_build(cache_status)
      if (cache_status/=trace_status_active) then
        call mpistop('mt_trace_spherical_surface_products failed to '// &
             'build spherical curl cache')
      endif
    endif

    if (present(b_min)) then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local,q_local,qperp_local) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call mt_trace_spherical_surface_seed(seed_local,dL,max_steps, &
             topology(iseed),q_local,qperp_local,do_twist,do_q, &
             do_qperp,b_min)
        if (do_q) q_results(iseed)=q_local
        if (do_qperp) qperp_results(iseed)=qperp_local
      enddo
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local,q_local,qperp_local) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call mt_trace_spherical_surface_seed(seed_local,dL,max_steps, &
             topology(iseed),q_local,qperp_local,do_twist,do_q, &
             do_qperp)
        if (do_q) q_results(iseed)=q_local
        if (do_qperp) qperp_results(iseed)=qperp_local
      enddo
      !$OMP END PARALLEL DO
    endif

    if (do_twist) call trace_spherical_curl_cache_clear()
  end subroutine mt_trace_spherical_surface_products

  subroutine mt_trace_spherical_surface_seed(seed,dL,max_steps,topology, &
       q_result,qperp_result,do_twist,do_q,do_qperp,b_min)
    integer, intent(in) :: max_steps
    double precision, intent(in) :: seed(ndim),dL
    type(trace_topology_result), intent(out) :: topology
    type(trace_qperp_result), intent(out) :: q_result,qperp_result
    logical, intent(in) :: do_twist,do_q,do_qperp
    double precision, intent(in), optional :: b_min

    type(trace_topology_result) :: topology_one(1)
    type(trace_qperp_result) :: q_one(1),qperp_one(1)
    type(trace_twist_result) :: twist_one(1)
    double precision :: seed_one(1,ndim),source_normal(3)

    seed_one(1,:)=seed
    source_normal=0.d0
    source_normal(1)=-1.d0

    if (do_q .and. do_qperp) then
      if (present(b_min)) then
        if (do_twist) then
          call trace_field_spherical_rmin_q_qperp_multi(seed_one,1,dL, &
               max_steps,q_one,qperp_one,b_min,twist_results=twist_one)
        else
          call trace_field_spherical_rmin_q_qperp_multi(seed_one,1,dL, &
               max_steps,q_one,qperp_one,b_min)
        endif
      else
        if (do_twist) then
          call trace_field_spherical_rmin_q_qperp_multi(seed_one,1,dL, &
               max_steps,q_one,qperp_one,twist_results=twist_one)
        else
          call trace_field_spherical_rmin_q_qperp_multi(seed_one,1,dL, &
               max_steps,q_one,qperp_one)
        endif
      endif
      call mt_q0_trace_to_topology(q_one,twist_one,1,topology_one,do_twist)
      topology=topology_one(1)
      q_result=q_one(1)
      qperp_result=qperp_one(1)
      return
    endif

    if (do_q) then
      if (present(b_min)) then
        if (do_twist) then
          call trace_field_spherical_rmin_q_multi(seed_one,1,dL,max_steps, &
               q_one,b_min,twist_results=twist_one)
        else
          call trace_field_spherical_rmin_q_multi(seed_one,1,dL,max_steps, &
               q_one,b_min)
        endif
      else
        if (do_twist) then
          call trace_field_spherical_rmin_q_multi(seed_one,1,dL,max_steps, &
               q_one,twist_results=twist_one)
        else
          call trace_field_spherical_rmin_q_multi(seed_one,1,dL,max_steps, &
               q_one)
        endif
      endif
      call mt_q0_trace_to_topology(q_one,twist_one,1,topology_one,do_twist)
      topology=topology_one(1)
      q_result=q_one(1)
      return
    endif

    if (present(b_min)) then
      call trace_field_topology_multi(seed_one,1,dL,max_steps,topology_one, &
           need_twist=do_twist,need_mapping=.false.,b_min=b_min, &
           source_normal=source_normal)
    else
      call trace_field_topology_multi(seed_one,1,dL,max_steps,topology_one, &
           need_twist=do_twist,need_mapping=.false., &
           source_normal=source_normal)
    endif
    topology=topology_one(1)

    if (do_qperp) then
      if (present(b_min)) then
        call trace_field_spherical_qperp_multi(seed_one,1,dL,max_steps, &
             qperp_one,b_min)
      else
        call trace_field_spherical_qperp_multi(seed_one,1,dL,max_steps, &
             qperp_one)
      endif
      qperp_result=qperp_one(1)
    endif
  end subroutine mt_trace_spherical_surface_seed

  subroutine mt_trace_spherical_cloud_products(seeds,nseed,dL,max_steps, &
       topology,q_results,qperp_results,do_twist,do_q,do_qperp,b_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_topology_result), intent(out) :: topology(nseed)
    type(trace_qperp_result), intent(out) :: q_results(:)
    type(trace_qperp_result), intent(out) :: qperp_results(:)
    logical, intent(in) :: do_twist,do_q,do_qperp
    double precision, intent(in), optional :: b_min

    type(trace_topology_result) :: topology_one(1)
    type(trace_twist_result), allocatable :: twist_results(:)
    double precision :: seed_local(ndim),seed_one(1,ndim)
    integer :: iseed,cache_status
    character(len=mt_task_name_len) :: integrator

    if (nseed<=0) return

    integrator=mt_lowercase(trim(mt_trace_integrator))
    if (do_q .and. do_qperp .and. &
         (integrator=='rk2' .or. integrator=='rk45_spherical')) then
      ! The spherical grouped RK45 driver scans all active states by grid.
      ! For volume clouds this is much slower and much more memory-hungry than
      ! the surface product's seed-wise OpenMP schedule.
      if (present(b_min)) then
        call mt_trace_spherical_surface_products(seeds,nseed,dL,max_steps, &
             topology,q_results,qperp_results,do_twist,do_q,do_qperp,b_min)
      else
        call mt_trace_spherical_surface_products(seeds,nseed,dL,max_steps, &
             topology,q_results,qperp_results,do_twist,do_q,do_qperp)
      endif
      return
    endif

    if (do_q .and. .not.do_qperp) then
      if (do_twist) then
        allocate(twist_results(nseed))
        if (present(b_min)) then
          call trace_field_spherical_rmin_q_multi(seeds,nseed,dL,max_steps, &
               q_results,b_min,twist_results=twist_results)
        else
          call trace_field_spherical_rmin_q_multi(seeds,nseed,dL,max_steps, &
               q_results,twist_results=twist_results)
        endif
        call mt_q0_trace_to_topology(q_results,twist_results,nseed,topology, &
             do_twist)
        deallocate(twist_results)
      else
        allocate(twist_results(0))
        if (present(b_min)) then
          call trace_field_spherical_rmin_q_multi(seeds,nseed,dL,max_steps, &
               q_results,b_min)
        else
          call trace_field_spherical_rmin_q_multi(seeds,nseed,dL,max_steps, &
               q_results)
        endif
        call mt_q0_trace_to_topology(q_results,twist_results,nseed,topology, &
             do_twist)
        deallocate(twist_results)
      endif
      return
    endif

    if (do_qperp) then
      if (do_twist) then
        allocate(twist_results(nseed))
        if (present(b_min)) then
          call trace_field_spherical_qperp_multi(seeds,nseed,dL,max_steps, &
               qperp_results,b_min,twist_results=twist_results)
        else
          call trace_field_spherical_qperp_multi(seeds,nseed,dL,max_steps, &
               qperp_results,twist_results=twist_results)
        endif
        call mt_qperp_trace_to_topology(qperp_results,twist_results,nseed, &
             topology,do_twist)
        deallocate(twist_results)
      else
        allocate(twist_results(0))
        if (present(b_min)) then
          call trace_field_spherical_qperp_multi(seeds,nseed,dL,max_steps, &
               qperp_results,b_min)
        else
          call trace_field_spherical_qperp_multi(seeds,nseed,dL,max_steps, &
               qperp_results)
        endif
        call mt_qperp_trace_to_topology(qperp_results,twist_results,nseed, &
             topology,do_twist)
        deallocate(twist_results)
      endif
      return
    endif

    if (do_twist) then
      call trace_spherical_curl_cache_build(cache_status)
      if (cache_status/=trace_status_active) then
        call mpistop('mt_trace_spherical_cloud_products failed to '// &
             'build spherical curl cache')
      endif
    endif

    if (present(b_min)) then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local,seed_one,topology_one) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        seed_one(1,:)=seed_local
        call trace_field_topology_multi(seed_one,1,dL,max_steps, &
             topology_one,need_twist=do_twist,need_mapping=.false., &
             b_min=b_min)
        topology(iseed)=topology_one(1)
      enddo
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local,seed_one,topology_one) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        seed_one(1,:)=seed_local
        call trace_field_topology_multi(seed_one,1,dL,max_steps, &
             topology_one,need_twist=do_twist,need_mapping=.false.)
        topology(iseed)=topology_one(1)
      enddo
      !$OMP END PARALLEL DO
    endif

    if (do_twist) call trace_spherical_curl_cache_clear()
  end subroutine mt_trace_spherical_cloud_products

  subroutine mt_length_plane_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,csv_file,b_min)
    ! Trace a uniform seed grid on a constant-z Cartesian plane.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_length_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_length_plane_xy','ix,iy',b_min)
    else
      call mt_length_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_length_plane_xy','ix,iy')
    endif
  end subroutine mt_length_plane_xy

  subroutine mt_mapping_plane_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,csv_file,b_min)
    ! Trace a constant-z seed grid and write endpoint metadata for Q mapping.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_mapping_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    double precision :: source_normal(3)
    integer :: nseed

    call mt_validate_axis_plane(xmin,xmax,nx,ymin,ymax,ny, &
         'mt_mapping_plane_xy')
    call mt_build_axis_plane_seeds(xmin,xmax,nx,ymin,ymax,ny,z0, &
         1,2,3,seeds)

    source_normal=0.d0
    source_normal(3)=1.d0
    nseed=nx*ny
    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_mapping_multi(seeds,nseed,dL,max_steps,results, &
           b_min=b_min,source_normal=source_normal)
    else
      call trace_field_mapping_multi(seeds,nseed,dL,max_steps,results, &
           source_normal=source_normal)
    endif
    call mt_write_mapping_plane_xy_csv(results,nx,ny,csv_file)

    deallocate(seeds,results)
  end subroutine mt_mapping_plane_xy

  subroutine mt_q_plane_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,csv_file,b_min)
    ! Compute standard Cartesian logQ diagnostics on a constant-z seed plane.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_q_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_q_plane_xy','ix,iy',b_min)
    else
      call mt_q_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_q_plane_xy','ix,iy')
    endif
  end subroutine mt_q_plane_xy

  subroutine mt_q_plane_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,csv_file,b_min)
    ! Compute standard Cartesian logQ diagnostics on a constant-y seed plane.
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_q_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_q_plane_xz','ix,iz',b_min)
    else
      call mt_q_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_q_plane_xz','ix,iz')
    endif
  end subroutine mt_q_plane_xz

  subroutine mt_q_plane_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,csv_file,b_min)
    ! Compute standard Cartesian logQ diagnostics on a constant-x seed plane.
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_q_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_q_plane_yz','iy,iz',b_min)
    else
      call mt_q_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_q_plane_yz','iy,iz')
    endif
  end subroutine mt_q_plane_yz

  subroutine mt_qperp_plane_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,csv_file,b_min)
    ! Compute Method-II Qperp diagnostics on a constant-z source plane.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qperp_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_qperp_plane_xy','ix,iy',b_min)
    else
      call mt_qperp_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_qperp_plane_xy','ix,iy')
    endif
  end subroutine mt_qperp_plane_xy

  subroutine mt_qperp_plane_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,csv_file,b_min)
    ! Compute Method-II Qperp diagnostics on a constant-y source plane.
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qperp_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_qperp_plane_xz','ix,iz',b_min)
    else
      call mt_qperp_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_qperp_plane_xz','ix,iz')
    endif
  end subroutine mt_qperp_plane_xz

  subroutine mt_qperp_plane_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,csv_file,b_min)
    ! Compute Method-II Qperp diagnostics on a constant-x source plane.
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qperp_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_qperp_plane_yz','iy,iz',b_min)
    else
      call mt_qperp_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_qperp_plane_yz','iy,iz')
    endif
  end subroutine mt_qperp_plane_yz

  subroutine mt_qperp_plane_arbitrary(origin,e1,e2,s1min,s1max,n1, &
       s2min,s2max,n2,dL,max_steps,csv_file,b_min)
    ! Compute Method-II Qperp diagnostics on an orthonormal arbitrary plane.
    integer, intent(in) :: n1,n2,max_steps
    double precision, intent(in) :: origin(ndim),e1(ndim),e2(ndim)
    double precision, intent(in) :: s1min,s1max,s2min,s2max,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_qperp_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:),s1(:),s2(:)
    integer :: nseed

    if (.not.mt_validate_arbitrary_plane_basis(e1,e2,s1min,s1max,n1, &
         s2min,s2max,n2,'mt_qperp_plane_arbitrary')) then
      call mt_write_qperp_arbitrary_header(csv_file, &
           'mt_qperp_plane_arbitrary')
      return
    endif

    call mt_build_arbitrary_plane_seeds(origin,e1,e2,s1min,s1max,n1, &
         s2min,s2max,n2,seeds,s1,s2)

    nseed=n1*n2
    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results)
    endif
    call mt_write_qperp_arbitrary_csv(results,s1,s2,n1,n2,csv_file, &
         'mt_qperp_plane_arbitrary')

    deallocate(seeds,s1,s2,results)
  end subroutine mt_qperp_plane_arbitrary

  subroutine mt_fieldline_products_plane_arbitrary(origin,e1,e2,s1min, &
       s1max,n1,s2min,s2max,n2,dL,max_steps,csv_file,b_min, &
       compute_length,compute_twist,compute_q,compute_qperp,vtu_file, &
       write_csv)
    ! Write selected per-field-line diagnostics on an arbitrary seed plane.
    integer, intent(in) :: n1,n2,max_steps
    double precision, intent(in) :: origin(ndim),e1(ndim),e2(ndim)
    double precision, intent(in) :: s1min,s1max,s2min,s2max,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min
    logical, intent(in), optional :: compute_length,compute_twist
    logical, intent(in), optional :: compute_q,compute_qperp
    character(len=*), intent(in), optional :: vtu_file
    logical, intent(in), optional :: write_csv

    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_qperp_result), allocatable :: q_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    double precision, allocatable :: seeds(:,:),s1(:),s2(:)
    integer :: nseed
    logical :: do_length,do_twist,do_q,do_qperp,do_csv

    do_length=.true.
    do_twist=.false.
    do_q=.false.
    do_qperp=.false.
    do_csv=.true.
    if (present(compute_length)) do_length=compute_length
    if (present(compute_twist)) do_twist=compute_twist
    if (present(compute_q)) do_q=compute_q
    if (present(compute_qperp)) do_qperp=compute_qperp
    if (present(write_csv)) do_csv=write_csv

    if (.not.mt_validate_arbitrary_plane_basis(e1,e2,s1min,s1max,n1, &
         s2min,s2max,n2,'mt_fieldline_products_plane_arbitrary')) then
      if (do_csv) then
        call mt_write_fieldline_products_plane_arbitrary_header(csv_file, &
             'mt_fieldline_products_plane_arbitrary',do_twist,do_q,do_qperp)
      endif
      if (present(vtu_file)) then
        if (len_trim(vtu_file)>0) then
          call mt_write_fieldline_products_vtu_empty(vtu_file,do_twist, &
               do_q,do_qperp,'mt_fieldline_products_plane_arbitrary', &
               do_length=do_length)
        endif
      endif
      return
    endif
    call mt_require_requested_science(do_length,do_twist,do_q,do_qperp, &
         'mt_fieldline_products_plane_arbitrary')

    call mt_build_arbitrary_plane_seeds(origin,e1,e2,s1min,s1max,n1, &
         s2min,s2max,n2,seeds,s1,s2)

    nseed=n1*n2
    allocate(length_results(nseed))
    if (do_twist) then
      allocate(twist_results(nseed))
    else
      allocate(twist_results(0))
    endif
    if (do_q) then
      allocate(q_results(nseed))
    else
      allocate(q_results(0))
    endif
    if (do_qperp) then
      allocate(qperp_results(nseed))
    else
      allocate(qperp_results(0))
    endif

    if (present(b_min)) then
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,q_results,qperp_results,do_twist, &
           do_q,do_qperp,b_min)
    else
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,q_results,qperp_results,do_twist, &
           do_q,do_qperp)
    endif

    if (do_csv) then
      call mt_write_fieldline_products_plane_arbitrary_csv(length_results, &
           twist_results,q_results,qperp_results,s1,s2,n1,n2,csv_file, &
           do_twist,do_q,do_qperp,'mt_fieldline_products_plane_arbitrary')
    endif
    if (present(vtu_file)) then
      if (len_trim(vtu_file)>0) then
        call mt_write_fieldline_products_vtu_plane(vtu_file, &
             length_results,twist_results,q_results,qperp_results,n1,n2, &
             do_twist,do_q,do_qperp,'mt_fieldline_products_plane_arbitrary', &
             do_length=do_length)
      endif
    endif

    deallocate(qperp_results,q_results,twist_results,length_results,seeds,s1,s2)
  end subroutine mt_fieldline_products_plane_arbitrary

  subroutine mt_topology_plane_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,length_csv,twist_csv,mapping_csv,b_min)
    ! Trace a constant-z seed plane once and write selected topology products.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: length_csv,twist_csv,mapping_csv
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_topology_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_xy','ix,iy','ix,iy',b_min)
    else
      call mt_topology_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_xy','ix,iy','ix,iy')
    endif
  end subroutine mt_topology_plane_xy

  subroutine mt_topology_plane_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,length_csv,twist_csv,mapping_csv,b_min)
    ! Trace a constant-y seed plane once and write selected topology products.
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: length_csv,twist_csv,mapping_csv
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_topology_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_xz','ix,iz','i,j',b_min)
    else
      call mt_topology_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_xz','ix,iz','i,j')
    endif
  end subroutine mt_topology_plane_xz

  subroutine mt_topology_plane_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,length_csv,twist_csv,mapping_csv,b_min)
    ! Trace a constant-x seed plane once and write selected topology products.
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: length_csv,twist_csv,mapping_csv
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_topology_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_yz','iy,iz','i,j',b_min)
    else
      call mt_topology_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,length_csv,twist_csv,mapping_csv, &
           'mt_topology_plane_yz','iy,iz','i,j')
    endif
  end subroutine mt_topology_plane_yz

  subroutine mt_qsl_plane_vtu_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,vtu_file,b_min)
    ! Write a full topology/QSL visualization file for a constant-z plane.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: vtu_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qsl_plane_vtu_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_xy',b_min)
    else
      call mt_qsl_plane_vtu_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_xy')
    endif
  end subroutine mt_qsl_plane_vtu_xy

  subroutine mt_qsl_plane_vtu_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,vtu_file,b_min)
    ! Write a full topology/QSL visualization file for a constant-y plane.
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: vtu_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qsl_plane_vtu_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_xz',b_min)
    else
      call mt_qsl_plane_vtu_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_xz')
    endif
  end subroutine mt_qsl_plane_vtu_xz

  subroutine mt_qsl_plane_vtu_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,vtu_file,b_min)
    ! Write a full topology/QSL visualization file for a constant-x plane.
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: vtu_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qsl_plane_vtu_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_yz',b_min)
    else
      call mt_qsl_plane_vtu_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,vtu_file,'mt_qsl_plane_vtu_yz')
    endif
  end subroutine mt_qsl_plane_vtu_yz

  subroutine mt_qsl_plane_vti_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,vti_file,do_length,do_twist,do_q,do_qperp,b_min)
    ! Write requested minimal science arrays for a constant-z plane as VTI.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: vti_file
    logical, intent(in) :: do_length,do_twist,do_q,do_qperp
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qsl_plane_vti_axis(xmin,xmax,nx,ymin,ymax,ny,z0, &
           1,2,3,dL,max_steps,vti_file,do_twist,do_q,do_qperp, &
           'mt_qsl_plane_vti_xy',b_min,do_length=do_length)
    else
      call mt_qsl_plane_vti_axis(xmin,xmax,nx,ymin,ymax,ny,z0, &
           1,2,3,dL,max_steps,vti_file,do_twist,do_q,do_qperp, &
           'mt_qsl_plane_vti_xy',do_length=do_length)
    endif
  end subroutine mt_qsl_plane_vti_xy

  subroutine mt_qsl_plane_vti_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,vti_file,do_length,do_twist,do_q,do_qperp,b_min)
    ! Write requested minimal science arrays for a constant-y plane as VTI.
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: vti_file
    logical, intent(in) :: do_length,do_twist,do_q,do_qperp
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qsl_plane_vti_axis(xmin,xmax,nx,zmin,zmax,nz,y0, &
           1,3,2,dL,max_steps,vti_file,do_twist,do_q,do_qperp, &
           'mt_qsl_plane_vti_xz',b_min,do_length=do_length)
    else
      call mt_qsl_plane_vti_axis(xmin,xmax,nx,zmin,zmax,nz,y0, &
           1,3,2,dL,max_steps,vti_file,do_twist,do_q,do_qperp, &
           'mt_qsl_plane_vti_xz',do_length=do_length)
    endif
  end subroutine mt_qsl_plane_vti_xz

  subroutine mt_qsl_plane_vti_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,vti_file,do_length,do_twist,do_q,do_qperp,b_min)
    ! Write requested minimal science arrays for a constant-x plane as VTI.
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: vti_file
    logical, intent(in) :: do_length,do_twist,do_q,do_qperp
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_qsl_plane_vti_axis(ymin,ymax,ny,zmin,zmax,nz,x0, &
           2,3,1,dL,max_steps,vti_file,do_twist,do_q,do_qperp, &
           'mt_qsl_plane_vti_yz',b_min,do_length=do_length)
    else
      call mt_qsl_plane_vti_axis(ymin,ymax,ny,zmin,zmax,nz,x0, &
           2,3,1,dL,max_steps,vti_file,do_twist,do_q,do_qperp, &
           'mt_qsl_plane_vti_yz',do_length=do_length)
    endif
  end subroutine mt_qsl_plane_vti_yz

  subroutine mt_write_cartesian_vti_pointdata(vti_file,xmin,xmax,nx, &
       ymin,ymax,ny,zmin,zmax,nz,length_total,qperp,status)
    ! Low-level appended-binary VTI writer for future uniform seed-volume
    ! products. PointData values use i-fastest ordering:
    ! ipoint = (k-1)*nx*ny + (j-1)*nx + i.
    character(len=*), intent(in) :: vti_file
    integer, intent(in) :: nx,ny,nz
    double precision, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
    double precision, intent(in) :: length_total(:),qperp(:)
    integer, intent(in) :: status(:)

    integer :: npoint
    double precision :: origin(3),spacing(3)

    if (len_trim(vti_file)==0) then
      call mpistop('mt_write_cartesian_vti_pointdata requires a VTI file')
    endif
    if (nx<=0 .or. ny<=0 .or. nz<=0) then
      call mpistop('mt_write_cartesian_vti_pointdata requires positive dimensions')
    endif

    npoint=nx*ny*nz
    if (size(length_total)/=npoint .or. size(qperp)/=npoint .or. &
         size(status)/=npoint) then
      call mpistop('mt_write_cartesian_vti_pointdata array size mismatch')
    endif

    origin=(/ xmin,ymin,zmin /)
    spacing(1)=mt_vti_axis_spacing(xmin,xmax,nx)
    spacing(2)=mt_vti_axis_spacing(ymin,ymax,ny)
    spacing(3)=mt_vti_axis_spacing(zmin,zmax,nz)
    call mt_write_vti_pointdata_fixed(vti_file,origin,spacing,nx,ny,nz, &
         length_total,qperp,status,'mt_write_cartesian_vti_pointdata')
  end subroutine mt_write_cartesian_vti_pointdata

  subroutine mt_fieldline_products_volume_vti(xmin,xmax,nx, &
       ymin,ymax,ny,zmin,zmax,nz,dL,max_steps,vti_file,b_min, &
       compute_length,compute_twist,compute_q,compute_qperp,chunk_nz)
    ! Compute per-seed field-line products on a user-defined uniform
    ! Cartesian sampling volume and write PointData to appended-binary VTI.
    integer, intent(in) :: nx,ny,nz,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax,dL
    character(len=*), intent(in) :: vti_file
    double precision, intent(in), optional :: b_min
    logical, intent(in), optional :: compute_length,compute_twist
    logical, intent(in), optional :: compute_q,compute_qperp
    integer, intent(in), optional :: chunk_nz

    type(trace_length_result), allocatable :: length_slab(:)
    type(trace_twist_result), allocatable :: twist_slab(:)
    type(trace_qperp_result), allocatable :: q_slab(:)
    type(trace_qperp_result), allocatable :: qperp_slab(:)
    type(mt_volume_products) :: products
    double precision, allocatable :: seeds_slab(:,:)
    double precision :: origin(3),spacing(3)
    integer :: k_start,k_end,slab_nz,slab_nseed,nseed
    integer :: chunk_nz_eff
    integer(kind=8) :: nseed64
    logical :: do_length,do_twist,do_q,do_qperp

    if (len_trim(vti_file)==0) then
      call mpistop('mt_fieldline_products_volume_vti requires a VTI file')
    endif
    if (npe/=1) then
      call mpistop('mt_fieldline_products_volume_vti currently requires npe=1')
    endif
    if (ndim/=3) then
      call mpistop('mt_fieldline_products_volume_vti requires 3D Cartesian geometry')
    endif
    {^IFTHREED
    select case (geo_coordinate)
    case (geo_cartesian)
      ! The output VTI is a user-defined sampling lattice; RK2 tracing can
      ! interpolate through either slab-uniform or AMR Cartesian grids.
    case (geo_cartesian_stretched)
    case default
      call mpistop('volume_vti requires Cartesian geometry')
    end select
    }
    if (nx<=0 .or. ny<=0 .or. nz<=0) then
      call mpistop('mt_fieldline_products_volume_vti requires positive dimensions')
    endif
    if (xmax<xmin .or. ymax<ymin .or. zmax<zmin) then
      call mpistop('mt_fieldline_products_volume_vti requires ordered bounds')
    endif

    nseed64=int(nx,kind=8)*int(ny,kind=8)*int(nz,kind=8)
    if (nseed64>int(huge(nseed),kind=8)) then
      call mpistop('mt_fieldline_products_volume_vti seed count overflows integer')
    endif
    if (nseed64>1000000_8) then
      write(*,'(a,i0,a)') &
           'mt_fieldline_products_volume_vti warning: volume has ', &
           nseed64,' seeds; final VTI arrays still scale with nseed'
    endif
    nseed=int(nseed64)

    do_length=.true.
    do_twist=.false.
    do_q=.false.
    do_qperp=.false.
    if (present(compute_length)) do_length=compute_length
    if (present(compute_twist)) do_twist=compute_twist
    if (present(compute_q)) do_q=compute_q
    if (present(compute_qperp)) do_qperp=compute_qperp
    call mt_require_requested_science(do_length,do_twist,do_q,do_qperp, &
         'mt_fieldline_products_volume_vti')

    chunk_nz_eff=nz
    if (present(chunk_nz)) chunk_nz_eff=max(1,min(nz,chunk_nz))

    origin=(/ xmin,ymin,zmin /)
    spacing(1)=mt_vti_axis_spacing(xmin,xmax,nx)
    spacing(2)=mt_vti_axis_spacing(ymin,ymax,ny)
    spacing(3)=mt_vti_axis_spacing(zmin,zmax,nz)

    call mt_allocate_volume_products(products,nseed,do_twist,do_q,do_qperp)

    do k_start=1,nz,chunk_nz_eff
      k_end=min(nz,k_start+chunk_nz_eff-1)
      slab_nz=k_end-k_start+1
      slab_nseed=nx*ny*slab_nz

      allocate(seeds_slab(slab_nseed,ndim))
      call mt_build_volume_slab_seeds(seeds_slab,xmin,ymin,zmin, &
           spacing,nx,ny,k_start,slab_nz)

      allocate(length_slab(slab_nseed))
      if (do_twist) then
        allocate(twist_slab(slab_nseed))
      else
        allocate(twist_slab(0))
      endif
      if (do_q) then
        allocate(q_slab(slab_nseed))
      else
        allocate(q_slab(0))
      endif

      if (do_qperp) then
        allocate(qperp_slab(slab_nseed))
      else
        allocate(qperp_slab(0))
      endif

      if (present(b_min)) then
        call mt_trace_fieldline_products_seedset(seeds_slab,slab_nseed, &
             dL,max_steps,length_slab,twist_slab,q_slab,qperp_slab, &
             do_twist,do_q,do_qperp,b_min)
      else
        call mt_trace_fieldline_products_seedset(seeds_slab,slab_nseed, &
             dL,max_steps,length_slab,twist_slab,q_slab,qperp_slab, &
             do_twist,do_q,do_qperp)
      endif

      call mt_copy_volume_length_slab(products,length_slab,nx,ny, &
           k_start,slab_nz)
      deallocate(length_slab)

      if (do_twist) then
        call mt_copy_volume_twist_slab(products,twist_slab,nx,ny, &
             k_start,slab_nz)
      endif
      deallocate(twist_slab)

      if (do_q) then
        call mt_copy_volume_q_slab(products,q_slab,nx,ny,k_start,slab_nz)
      endif
      deallocate(q_slab)

      if (do_qperp) then
        call mt_copy_volume_qperp_slab(products,qperp_slab,nx,ny, &
             k_start,slab_nz)
      endif
      deallocate(qperp_slab)

      deallocate(seeds_slab)
    enddo

    call mt_write_fieldline_products_volume_vti(vti_file,origin,spacing, &
         nx,ny,nz,products, &
         do_length,do_twist,do_q,do_qperp,'mt_fieldline_products_volume_vti')
    call mt_deallocate_volume_products(products)
  end subroutine mt_fieldline_products_volume_vti

  subroutine mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
       length_results,twist_results,q_results,qperp_results,do_twist, &
       do_q,do_qperp,b_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_length_result), intent(out) :: length_results(nseed)
    type(trace_twist_result), intent(out) :: twist_results(:)
    type(trace_qperp_result), intent(out) :: q_results(:)
    type(trace_qperp_result), intent(out) :: qperp_results(:)
    logical, intent(in) :: do_twist,do_q,do_qperp
    double precision, intent(in), optional :: b_min

    type(trace_topology_result) :: topology_one(1)
    double precision :: seed_local(ndim),seed_one(1,ndim)
    integer :: iseed
    character(len=mt_task_name_len) :: integrator

    integrator=mt_lowercase(trim(mt_trace_integrator))
    if (do_q .and. do_qperp .and. geo_coordinate==geo_spherical .and. &
         (integrator=='rk2' .or. integrator=='rk45_spherical')) then
      if (do_twist) then
        if (present(b_min)) then
          call trace_field_spherical_rmin_q_qperp_multi(seeds,nseed,dL, &
               max_steps,q_results,qperp_results,b_min, &
               twist_results=twist_results)
        else
          call trace_field_spherical_rmin_q_qperp_multi(seeds,nseed,dL, &
               max_steps,q_results,qperp_results, &
               twist_results=twist_results)
        endif
      else
        if (present(b_min)) then
          call trace_field_spherical_rmin_q_qperp_multi(seeds,nseed,dL, &
               max_steps,q_results,qperp_results,b_min)
        else
          call trace_field_spherical_rmin_q_qperp_multi(seeds,nseed,dL, &
               max_steps,q_results,qperp_results)
        endif
      endif
      call mt_q0_trace_to_length(q_results,nseed,length_results)
      return
    endif

    if (do_q .and. do_qperp .and. geo_coordinate/=geo_spherical .and. &
         (integrator=='rk2' .or. integrator=='rk45_cartesian')) then
      if (do_twist) then
        if (present(b_min)) then
          call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
               qperp_results,b_min,twist_results=twist_results)
        else
          call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
               qperp_results,twist_results=twist_results)
        endif
      else
        if (present(b_min)) then
          call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
               qperp_results,b_min)
        else
          call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
               qperp_results)
        endif
      endif
      q_results=qperp_results
      call mt_qperp_trace_to_length(qperp_results,nseed,length_results)
      return
    endif

    if (do_qperp .and. .not.do_q .and. &
         (integrator=='rk2' .or. &
         (integrator=='rk45_cartesian' .and. geo_coordinate/=geo_spherical) &
         .or. (integrator=='rk45_spherical' .and. &
         geo_coordinate==geo_spherical))) then
      select case (geo_coordinate)
      case (geo_spherical)
        if (do_twist) then
          if (present(b_min)) then
            call trace_field_spherical_qperp_multi(seeds,nseed,dL,max_steps, &
                 qperp_results,b_min,twist_results=twist_results)
          else
            call trace_field_spherical_qperp_multi(seeds,nseed,dL,max_steps, &
                 qperp_results,twist_results=twist_results)
          endif
        else
          if (present(b_min)) then
            call trace_field_spherical_qperp_multi(seeds,nseed,dL,max_steps, &
                 qperp_results,b_min)
          else
            call trace_field_spherical_qperp_multi(seeds,nseed,dL,max_steps, &
                 qperp_results)
          endif
        endif
      case default
        if (do_twist) then
          if (present(b_min)) then
            call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
                 qperp_results,b_min,twist_results=twist_results)
          else
            call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
                 qperp_results,twist_results=twist_results)
          endif
        else
          if (present(b_min)) then
            call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
                 qperp_results,b_min)
          else
            call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
                 qperp_results)
          endif
        endif
      end select
      call mt_qperp_trace_to_length(qperp_results,nseed,length_results)
      return
    endif

    if (do_q .and. .not.do_qperp .and. geo_coordinate==geo_spherical .and. &
         (integrator=='rk2' .or. integrator=='rk45_spherical')) then
      if (do_twist) then
        if (present(b_min)) then
          call trace_field_spherical_rmin_q_multi(seeds,nseed,dL,max_steps, &
               q_results,b_min,twist_results=twist_results)
        else
          call trace_field_spherical_rmin_q_multi(seeds,nseed,dL,max_steps, &
               q_results,twist_results=twist_results)
        endif
      else
        if (present(b_min)) then
          call trace_field_spherical_rmin_q_multi(seeds,nseed,dL,max_steps, &
               q_results,b_min)
        else
          call trace_field_spherical_rmin_q_multi(seeds,nseed,dL,max_steps, &
               q_results)
        endif
      endif
      call mt_q0_trace_to_length(q_results,nseed,length_results)
      return
    endif

    if (do_q .and. .not.do_qperp .and. &
         (integrator=='rk2' .or. &
         (integrator=='rk45_cartesian' .and. geo_coordinate/=geo_spherical) &
         .or. (integrator=='rk45_spherical' .and. &
         geo_coordinate==geo_spherical))) then
      if (integrator=='rk2') then
        if (geo_coordinate==geo_cartesian_stretched .or. &
             (geo_coordinate==geo_cartesian .and. .not.slab_uniform)) then
          if (do_twist) then
            if (present(b_min)) then
              call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
                   q_results,b_min,twist_results=twist_results)
            else
              call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
                   q_results,twist_results=twist_results)
            endif
          else
            if (present(b_min)) then
              call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
                   q_results,b_min)
            else
              call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
                   q_results)
            endif
          endif
        else
          if (do_twist) then
            if (present(b_min)) then
              call trace_field_rk2_short_boundary_q_multi(seeds,nseed,dL, &
                   max_steps,q_results,b_min,twist_results=twist_results)
            else
              call trace_field_rk2_short_boundary_q_multi(seeds,nseed,dL, &
                   max_steps,q_results,twist_results=twist_results)
            endif
          else
            if (present(b_min)) then
              call trace_field_rk2_short_boundary_q_multi(seeds,nseed,dL, &
                   max_steps,q_results,b_min)
            else
              call trace_field_rk2_short_boundary_q_multi(seeds,nseed,dL, &
                   max_steps,q_results)
            endif
          endif
        endif
      else if (integrator=='rk45_cartesian') then
        if (do_twist) then
          if (present(b_min)) then
            call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
                 q_results,b_min,twist_results=twist_results)
          else
            call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
                 q_results,twist_results=twist_results)
          endif
        else
          if (present(b_min)) then
            call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
                 q_results,b_min)
          else
            call trace_field_qperp_multi(seeds,nseed,dL,max_steps,q_results)
          endif
        endif
      else
        if (do_twist) then
          if (present(b_min)) then
            call trace_field_spherical_rmin_q_multi(seeds,nseed,dL, &
                 max_steps,q_results,b_min,twist_results=twist_results)
          else
            call trace_field_spherical_rmin_q_multi(seeds,nseed,dL, &
                 max_steps,q_results,twist_results=twist_results)
          endif
        else
          if (present(b_min)) then
            call trace_field_spherical_rmin_q_multi(seeds,nseed,dL, &
                 max_steps,q_results,b_min)
          else
            call trace_field_spherical_rmin_q_multi(seeds,nseed,dL, &
                 max_steps,q_results)
          endif
        endif
      endif
      call mt_q0_trace_to_length(q_results,nseed,length_results)
      return
    endif

    if (present(b_min)) then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local,seed_one,topology_one) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        if (do_twist) then
          seed_one(1,:)=seed_local
          call trace_field_topology_multi(seed_one,1,dL,max_steps, &
               topology_one,need_twist=.true.,need_mapping=.false., &
               b_min=b_min)
          call mt_topology_to_length(topology_one,1, &
               length_results(iseed:iseed))
          call mt_topology_to_twist(topology_one,1, &
               twist_results(iseed:iseed))
        else
          call trace_field_length_single(seed_local,dL,max_steps, &
               length_results(iseed),b_min)
        endif
        if (do_qperp) then
          call trace_field_qperp_single(seed_local,dL,max_steps, &
               qperp_results(iseed),b_min)
        endif
      enddo
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local,seed_one,topology_one) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        if (do_twist) then
          seed_one(1,:)=seed_local
          call trace_field_topology_multi(seed_one,1,dL,max_steps, &
               topology_one,need_twist=.true.,need_mapping=.false.)
          call mt_topology_to_length(topology_one,1, &
               length_results(iseed:iseed))
          call mt_topology_to_twist(topology_one,1, &
               twist_results(iseed:iseed))
        else
          call trace_field_length_single(seed_local,dL,max_steps, &
               length_results(iseed))
        endif
        if (do_qperp) then
          call trace_field_qperp_single(seed_local,dL,max_steps, &
               qperp_results(iseed))
        endif
      enddo
      !$OMP END PARALLEL DO
    endif

    if (do_q) then
      if (mt_lowercase(trim(mt_trace_integrator))=='rk45_cartesian') then
        if (present(b_min)) then
          call trace_field_qperp_multi(seeds,nseed,dL,max_steps,q_results, &
               b_min)
        else
          call trace_field_qperp_multi(seeds,nseed,dL,max_steps,q_results)
        endif
      else if (geo_coordinate==geo_spherical) then
        if (present(b_min)) then
          call trace_field_spherical_rmin_q_multi(seeds,nseed,dL, &
               max_steps,q_results,b_min)
        else
          call trace_field_spherical_rmin_q_multi(seeds,nseed,dL, &
               max_steps,q_results)
        endif
      else if (do_qperp) then
        q_results=qperp_results
      else
        if (present(b_min)) then
          call trace_field_qperp_multi(seeds,nseed,dL,max_steps,q_results, &
               b_min)
        else
          call trace_field_qperp_multi(seeds,nseed,dL,max_steps,q_results)
        endif
      endif
    endif
  end subroutine mt_trace_fieldline_products_seedset

  subroutine mt_twist_single(seed,dL,max_steps,csv_file,b_min)
    ! Trace one magnetic field line and write its length and twist summary.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_twist_result) :: result
    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status

    if (npe/=1) then
      call mpistop('mt_twist_single currently requires npe=1')
    endif
    if (ndim/=3 .or. .not.slab_uniform) then
      call mpistop('mt_twist_single requires 3D uniform Cartesian geometry')
    endif

    if (present(b_min)) then
      call trace_field_twist_single(seed,dL,max_steps,result,b_min)
    else
      call trace_field_twist_single(seed,dL,max_steps,result)
    endif

    seed_xyz=0.d0
    seed_xyz(1:ndim)=result%line%seed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop('mt_twist_single could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'seed_id,seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'twist_total,twist_backward,twist_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_twist_single could not write CSV header')
    endif

    write(csv_unit,'(i0,9(",",es24.16),4(",",i0))',iostat=io_status) &
         1,seed_xyz,result%line%total_length, &
         result%line%backward_length,result%line%forward_length, &
         result%total_twist,result%backward_twist,result%forward_twist, &
         result%line%backward_nstep,result%line%forward_nstep, &
         result%line%backward_status,result%line%forward_status
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_twist_single could not write CSV data')
    endif

    close(csv_unit)
  end subroutine mt_twist_single

  subroutine mt_twist_seeds(seeds,nseed,dL,max_steps,csv_file,b_min)
    ! Trace multiple magnetic field lines and write length and twist summaries.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    type(trace_twist_result), allocatable :: results(:)
    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status,iseed

    if (npe/=1) then
      call mpistop('mt_twist_seeds currently requires npe=1')
    endif
    if (ndim/=3 .or. .not.slab_uniform) then
      call mpistop('mt_twist_seeds requires 3D uniform Cartesian geometry')
    endif
    if (nseed<0) then
      call mpistop('mt_twist_seeds requires nseed>=0')
    endif

    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_twist_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_twist_multi(seeds,nseed,dL,max_steps,results)
    endif

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      deallocate(results)
      call mpistop('mt_twist_seeds could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'seed_id,seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'twist_total,twist_backward,twist_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      deallocate(results)
      call mpistop('mt_twist_seeds could not write CSV header')
    endif

    do iseed=1,nseed
      seed_xyz=0.d0
      seed_xyz(1:ndim)=results(iseed)%line%seed
      write(csv_unit,'(i0,9(",",es24.16),4(",",i0))',iostat=io_status) &
           iseed,seed_xyz,results(iseed)%line%total_length, &
           results(iseed)%line%backward_length, &
           results(iseed)%line%forward_length, &
           results(iseed)%total_twist,results(iseed)%backward_twist, &
           results(iseed)%forward_twist, &
           results(iseed)%line%backward_nstep, &
           results(iseed)%line%forward_nstep, &
           results(iseed)%line%backward_status, &
           results(iseed)%line%forward_status
      if (io_status/=0) then
        close(csv_unit)
        deallocate(results)
        call mpistop('mt_twist_seeds could not write CSV data')
      endif
    enddo

    close(csv_unit)
    deallocate(results)
  end subroutine mt_twist_seeds

  subroutine mt_twist_plane_xy(xmin,xmax,nx,ymin,ymax,ny,z0,dL, &
       max_steps,csv_file,b_min)
    ! Trace a uniform seed grid on a constant-z plane and write twist summaries.
    integer, intent(in) :: nx,ny,max_steps
    double precision, intent(in) :: xmin,xmax,ymin,ymax,z0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_twist_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_twist_plane_xy','ix,iy',b_min)
    else
      call mt_twist_plane_axis(xmin,xmax,nx,ymin,ymax,ny,z0,1,2,3, &
           dL,max_steps,csv_file,'mt_twist_plane_xy','ix,iy')
    endif
  end subroutine mt_twist_plane_xy

  subroutine mt_length_plane_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,csv_file,b_min)
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_length_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_length_plane_xz','i,j',b_min)
    else
      call mt_length_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_length_plane_xz','i,j')
    endif
  end subroutine mt_length_plane_xz

  subroutine mt_length_plane_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,csv_file,b_min)
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_length_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_length_plane_yz','i,j',b_min)
    else
      call mt_length_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_length_plane_yz','i,j')
    endif
  end subroutine mt_length_plane_yz

  subroutine mt_twist_plane_xz(xmin,xmax,nx,zmin,zmax,nz,y0,dL, &
       max_steps,csv_file,b_min)
    integer, intent(in) :: nx,nz,max_steps
    double precision, intent(in) :: xmin,xmax,zmin,zmax,y0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_twist_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_twist_plane_xz','i,j',b_min)
    else
      call mt_twist_plane_axis(xmin,xmax,nx,zmin,zmax,nz,y0,1,3,2, &
           dL,max_steps,csv_file,'mt_twist_plane_xz','i,j')
    endif
  end subroutine mt_twist_plane_xz

  subroutine mt_twist_plane_yz(ymin,ymax,ny,zmin,zmax,nz,x0,dL, &
       max_steps,csv_file,b_min)
    integer, intent(in) :: ny,nz,max_steps
    double precision, intent(in) :: ymin,ymax,zmin,zmax,x0,dL
    character(len=*), intent(in) :: csv_file
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call mt_twist_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_twist_plane_yz','i,j',b_min)
    else
      call mt_twist_plane_axis(ymin,ymax,ny,zmin,zmax,nz,x0,2,3,1, &
           dL,max_steps,csv_file,'mt_twist_plane_yz','i,j')
    endif
  end subroutine mt_twist_plane_yz

  subroutine mt_length_plane_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,csv_file, &
       caller,index_header,b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: csv_file,caller,index_header
    double precision, intent(in), optional :: b_min

    type(trace_length_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)

    nseed=n1*n2
    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_length_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_length_multi(seeds,nseed,dL,max_steps,results)
    endif
    call mt_write_length_plane_csv(results,n1,n2,csv_file,caller, &
         index_header)

    deallocate(seeds,results)
  end subroutine mt_length_plane_axis

  subroutine mt_twist_plane_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,csv_file, &
       caller,index_header,b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: csv_file,caller,index_header
    double precision, intent(in), optional :: b_min

    type(trace_twist_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)

    nseed=n1*n2
    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_twist_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_twist_multi(seeds,nseed,dL,max_steps,results)
    endif
    call mt_write_twist_plane_csv(results,n1,n2,csv_file,caller, &
         index_header)

    deallocate(seeds,results)
  end subroutine mt_twist_plane_axis

  subroutine mt_q_plane_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,csv_file,caller, &
       index_header,b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: csv_file,caller,index_header
    double precision, intent(in), optional :: b_min

    type(trace_qperp_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)

    nseed=n1*n2
    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results)
    endif
    call mt_write_q_plane_csv(results,n1,n2,csv_file,caller,index_header)

    deallocate(seeds,results)
  end subroutine mt_q_plane_axis

  subroutine mt_qperp_plane_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,csv_file, &
       caller,index_header,b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: csv_file,caller,index_header
    double precision, intent(in), optional :: b_min

    type(trace_qperp_result), allocatable :: results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)

    nseed=n1*n2
    allocate(results(nseed))
    if (present(b_min)) then
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results,b_min)
    else
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps,results)
    endif
    call mt_write_qperp_plane_csv(results,n1,n2,csv_file,caller, &
         index_header)

    deallocate(seeds,results)
  end subroutine mt_qperp_plane_axis

  subroutine mt_axis_plane_products_csv_axis(c1min,c1max,n1,c2min,c2max, &
       n2,fixed_value,axis1,axis2,fixed_axis,dL,max_steps,length_csv, &
       twist_csv,q_csv,qperp_csv,caller,index_header,length_index_header, &
       b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: length_csv,twist_csv,q_csv,qperp_csv
    character(len=*), intent(in) :: caller,index_header,length_index_header
    double precision, intent(in), optional :: b_min

    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_qperp_result), allocatable :: q_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    double precision, allocatable :: seeds(:,:)
    integer :: nseed
    logical :: do_twist,do_q,do_qperp

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    if (len_trim(length_csv)==0) then
      call mpistop(trim(caller)//' requires a length CSV file')
    endif

    do_twist=len_trim(twist_csv)>0
    do_q=len_trim(q_csv)>0
    do_qperp=len_trim(qperp_csv)>0

    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)
    nseed=n1*n2
    allocate(length_results(nseed))
    if (do_twist) then
      allocate(twist_results(nseed))
    else
      allocate(twist_results(0))
    endif
    if (do_q) then
      allocate(q_results(nseed))
    else
      allocate(q_results(0))
    endif
    if (do_qperp) then
      allocate(qperp_results(nseed))
    else
      allocate(qperp_results(0))
    endif

    if (present(b_min)) then
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,q_results,qperp_results, &
           do_twist,do_q,do_qperp,b_min)
    else
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,q_results,qperp_results, &
           do_twist,do_q,do_qperp)
    endif

    call mt_write_length_plane_csv(length_results,n1,n2,length_csv, &
         caller,length_index_header)
    if (do_twist) then
      call mt_write_twist_plane_csv(twist_results,n1,n2,twist_csv, &
           caller,length_index_header)
    endif
    if (do_q) then
      call mt_write_q_plane_csv(q_results,n1,n2,q_csv,caller, &
           index_header)
    endif
    if (do_qperp) then
      call mt_write_qperp_plane_csv(qperp_results,n1,n2,qperp_csv, &
           caller,index_header)
    endif

    deallocate(qperp_results,q_results,twist_results,length_results,seeds)
  end subroutine mt_axis_plane_products_csv_axis

  subroutine mt_topology_plane_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,length_csv, &
       twist_csv,mapping_csv,caller,index_header,length_index_header, &
       b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: length_csv,twist_csv,mapping_csv
    character(len=*), intent(in) :: caller,index_header,length_index_header
    double precision, intent(in), optional :: b_min

    type(trace_topology_result), allocatable :: topology(:)
    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_mapping_result), allocatable :: mapping_results(:)
    double precision, allocatable :: seeds(:,:)
    double precision :: source_normal(3)
    integer :: nseed
    logical :: write_twist,write_mapping,need_mapping

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    if (len_trim(length_csv)==0) then
      call mpistop(trim(caller)//' requires a length CSV file')
    endif

    write_twist=len_trim(twist_csv)>0
    write_mapping=len_trim(mapping_csv)>0
    need_mapping=write_mapping

    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)
    nseed=n1*n2
    allocate(topology(nseed),length_results(nseed))
    source_normal=0.d0
    source_normal(fixed_axis)=1.d0
    if (present(b_min)) then
      call trace_field_topology_multi(seeds,nseed,dL,max_steps,topology, &
           need_twist=write_twist,need_mapping=need_mapping,b_min=b_min, &
           source_normal=source_normal)
    else
      call trace_field_topology_multi(seeds,nseed,dL,max_steps,topology, &
           need_twist=write_twist,need_mapping=need_mapping, &
           source_normal=source_normal)
    endif

    call mt_topology_to_length(topology,nseed,length_results)
    call mt_write_length_plane_csv(length_results,n1,n2,length_csv, &
         caller,length_index_header)

    if (write_twist) then
      allocate(twist_results(nseed))
      call mt_topology_to_twist(topology,nseed,twist_results)
      call mt_write_twist_plane_csv(twist_results,n1,n2,twist_csv, &
           caller,length_index_header)
      deallocate(twist_results)
    endif

    if (need_mapping) then
      allocate(mapping_results(nseed))
      call mt_topology_to_mapping(topology,nseed,mapping_results)
      if (write_mapping) then
        call mt_write_mapping_plane_csv(mapping_results,n1,n2,mapping_csv, &
             caller,index_header)
      endif
      deallocate(mapping_results)
    endif

    deallocate(length_results,topology,seeds)
  end subroutine mt_topology_plane_axis

  subroutine mt_qsl_plane_vtu_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,vtu_file,caller, &
       b_min)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: vtu_file,caller
    double precision, intent(in), optional :: b_min

    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_mapping_result), allocatable :: mapping_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    double precision, allocatable :: seeds(:,:)
    double precision :: source_normal(3)
    integer :: nseed
    logical :: need_mapping

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    if (len_trim(vtu_file)==0) then
      call mpistop(trim(caller)//' requires a VTU file')
    endif

    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)
    nseed=n1*n2
    allocate(length_results(nseed),twist_results(nseed), &
         mapping_results(nseed),qperp_results(nseed))

    source_normal=0.d0
    source_normal(fixed_axis)=1.d0
    need_mapping=mt_vtk_detail_is_full()
    if (need_mapping) then
      if (present(b_min)) then
        call mt_trace_qsl_plane_full(seeds,nseed,dL,max_steps,source_normal, &
             qperp_results,twist_results,mapping_results,b_min)
      else
        call mt_trace_qsl_plane_full(seeds,nseed,dL,max_steps,source_normal, &
             qperp_results,twist_results,mapping_results)
      endif
      call mt_qperp_trace_to_length(qperp_results,nseed,length_results)
    else
      if (present(b_min)) then
        call mt_trace_qsl_plane_minimal(seeds,nseed,dL,max_steps, &
             qperp_results,twist_results,b_min)
      else
        call mt_trace_qsl_plane_minimal(seeds,nseed,dL,max_steps, &
             qperp_results,twist_results)
      endif
      call mt_qperp_trace_to_length(qperp_results,nseed,length_results)
    endif

    call mt_write_qsl_plane_vtu(vtu_file,length_results,twist_results, &
         mapping_results,qperp_results,n1,n2,caller)

    deallocate(qperp_results,mapping_results,twist_results, &
         length_results,seeds)
  end subroutine mt_qsl_plane_vtu_axis

  subroutine mt_qsl_plane_vti_axis(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,dL,max_steps,vti_file, &
       do_twist,do_q,do_qperp,caller,b_min,do_length)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis,max_steps
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    double precision, intent(in) :: fixed_value,dL
    character(len=*), intent(in) :: vti_file,caller
    logical, intent(in) :: do_twist,do_q,do_qperp
    double precision, intent(in), optional :: b_min
    logical, intent(in), optional :: do_length

    type(trace_length_result), allocatable :: length_results(:)
    type(trace_twist_result), allocatable :: twist_results(:)
    type(trace_qperp_result), allocatable :: q_results(:)
    type(trace_qperp_result), allocatable :: qperp_results(:)
    type(mt_volume_products) :: products
    double precision, allocatable :: seeds(:,:)
    double precision :: origin(3),spacing(3)
    integer :: nseed,nx_vti,ny_vti,nz_vti
    logical :: do_length_eff

    call mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    if (len_trim(vti_file)==0) then
      call mpistop(trim(caller)//' requires a VTI file')
    endif
    do_length_eff=.true.
    if (present(do_length)) do_length_eff=do_length
    call mt_require_requested_science(do_length_eff,do_twist,do_q,do_qperp, &
         caller)

    call mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
         fixed_value,axis1,axis2,fixed_axis,seeds)
    nseed=n1*n2

    allocate(length_results(nseed))
    if (do_twist) then
      allocate(twist_results(nseed))
    else
      allocate(twist_results(0))
    endif
    if (do_q) then
      allocate(q_results(nseed))
    else
      allocate(q_results(0))
    endif
    if (do_qperp) then
      allocate(qperp_results(nseed))
    else
      allocate(qperp_results(0))
    endif

    if (present(b_min)) then
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,q_results,qperp_results,do_twist, &
           do_q,do_qperp,b_min)
    else
      call mt_trace_fieldline_products_seedset(seeds,nseed,dL,max_steps, &
           length_results,twist_results,q_results,qperp_results,do_twist, &
           do_q,do_qperp)
    endif

    call mt_allocate_volume_products(products,nseed,do_twist,do_q,do_qperp)
    call mt_copy_volume_length_slab(products,length_results,nseed,1,1,1)
    if (do_twist) call mt_copy_volume_twist_slab(products,twist_results, &
         nseed,1,1,1)
    if (do_q) call mt_copy_volume_q_slab(products,q_results,nseed,1,1,1)
    if (do_qperp) call mt_copy_volume_qperp_slab(products,qperp_results, &
         nseed,1,1,1)

    origin=0.d0
    spacing=0.d0
    nx_vti=1
    ny_vti=1
    nz_vti=1
    origin(axis1)=c1min
    origin(axis2)=c2min
    origin(fixed_axis)=fixed_value
    spacing(axis1)=mt_vti_axis_spacing(c1min,c1max,n1)
    spacing(axis2)=mt_vti_axis_spacing(c2min,c2max,n2)
    select case (axis1)
    case (1)
      nx_vti=n1
    case (2)
      ny_vti=n1
    case (3)
      nz_vti=n1
    end select
    select case (axis2)
    case (1)
      nx_vti=n2
    case (2)
      ny_vti=n2
    case (3)
      nz_vti=n2
    end select

    call mt_write_fieldline_products_volume_vti(vti_file,origin,spacing, &
         nx_vti,ny_vti,nz_vti,products,do_length_eff,do_twist,do_q, &
         do_qperp,caller)

    call mt_deallocate_volume_products(products)
    deallocate(qperp_results,q_results,twist_results,length_results,seeds)
  end subroutine mt_qsl_plane_vti_axis

  subroutine mt_trace_qsl_plane_full(seeds,nseed,dL,max_steps,source_normal, &
       qperp_results,twist_results,mapping_results,b_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL,source_normal(3)
    type(trace_qperp_result), intent(out) :: qperp_results(nseed)
    type(trace_twist_result), intent(out) :: twist_results(nseed)
    type(trace_mapping_result), intent(out) :: mapping_results(nseed)
    double precision, intent(in), optional :: b_min

    double precision :: seed_local(ndim)
    integer :: iseed

    if (present(b_min)) then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_field_mapping_single(seed_local,dL,max_steps, &
             mapping_results(iseed),b_min,source_normal)
      enddo
      !$OMP END PARALLEL DO
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
           qperp_results,b_min,twist_results=twist_results)
    else
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iseed,seed_local) SCHEDULE(DYNAMIC,16)
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_field_mapping_single(seed_local,dL,max_steps, &
             mapping_results(iseed),source_normal=source_normal)
      enddo
      !$OMP END PARALLEL DO
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
           qperp_results,twist_results=twist_results)
    endif
  end subroutine mt_trace_qsl_plane_full

  subroutine mt_trace_qsl_plane_minimal(seeds,nseed,dL,max_steps, &
       qperp_results,twist_results,b_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_qperp_result), intent(out) :: qperp_results(nseed)
    type(trace_twist_result), intent(out) :: twist_results(nseed)
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
           qperp_results,b_min,twist_results=twist_results)
    else
      call trace_field_qperp_multi(seeds,nseed,dL,max_steps, &
           qperp_results,twist_results=twist_results)
    endif
  end subroutine mt_trace_qsl_plane_minimal

  subroutine mt_qperp_to_length(qperp_results,nseed,results)
    integer, intent(in) :: nseed
    type(trace_qperp_result), intent(in) :: qperp_results(nseed)
    type(trace_length_result), intent(out) :: results(nseed)

    integer :: iseed

    do iseed=1,nseed
      results(iseed)%seed=qperp_results(iseed)%seed
      results(iseed)%forward_footpoint=qperp_results(iseed)%forward_endpoint
      results(iseed)%backward_footpoint=qperp_results(iseed)%backward_endpoint
      results(iseed)%forward_status=qperp_results(iseed)%forward_status
      results(iseed)%backward_status=qperp_results(iseed)%backward_status
      results(iseed)%forward_nstep=0
      results(iseed)%backward_nstep=0
      if (qperp_results(iseed)%valid_q0) then
        results(iseed)%forward_length=qperp_results(iseed)%forward_length
        results(iseed)%backward_length=qperp_results(iseed)%backward_length
        results(iseed)%total_length=results(iseed)%forward_length &
             +results(iseed)%backward_length
      else
        results(iseed)%forward_length=0.d0
        results(iseed)%backward_length=0.d0
        results(iseed)%total_length=0.d0
      endif
    enddo
  end subroutine mt_qperp_to_length

  subroutine mt_qperp_trace_to_length(qperp_results,nseed,results)
    integer, intent(in) :: nseed
    type(trace_qperp_result), intent(in) :: qperp_results(nseed)
    type(trace_length_result), intent(out) :: results(nseed)

    integer :: iseed

    do iseed=1,nseed
      results(iseed)%seed=qperp_results(iseed)%seed
      results(iseed)%forward_footpoint=qperp_results(iseed)%forward_endpoint
      results(iseed)%backward_footpoint=qperp_results(iseed)%backward_endpoint
      results(iseed)%forward_length=qperp_results(iseed)%forward_length
      results(iseed)%backward_length=qperp_results(iseed)%backward_length
      results(iseed)%total_length=results(iseed)%forward_length+ &
           results(iseed)%backward_length
      results(iseed)%forward_nstep=qperp_results(iseed)%forward_nstep
      results(iseed)%backward_nstep=qperp_results(iseed)%backward_nstep
      results(iseed)%forward_status=qperp_results(iseed)%forward_status
      results(iseed)%backward_status=qperp_results(iseed)%backward_status
    enddo
  end subroutine mt_qperp_trace_to_length

  subroutine mt_qperp_trace_to_topology(qperp_results,twist_results,nseed, &
       topology,do_twist)
    integer, intent(in) :: nseed
    type(trace_qperp_result), intent(in) :: qperp_results(nseed)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_topology_result), intent(out) :: topology(nseed)
    logical, intent(in) :: do_twist

    integer :: iseed

    do iseed=1,nseed
      topology(iseed)%seed=qperp_results(iseed)%seed
      topology(iseed)%length_forward=qperp_results(iseed)%forward_length
      topology(iseed)%length_backward=qperp_results(iseed)%backward_length
      topology(iseed)%length_total=topology(iseed)%length_forward &
           +topology(iseed)%length_backward
      topology(iseed)%twist_forward=0.d0
      topology(iseed)%twist_backward=0.d0
      topology(iseed)%twist_total=0.d0
      if (do_twist) then
        topology(iseed)%twist_forward=twist_results(iseed)%forward_twist
        topology(iseed)%twist_backward=twist_results(iseed)%backward_twist
        topology(iseed)%twist_total=twist_results(iseed)%total_twist
      endif
      topology(iseed)%forward_endpoint=qperp_results(iseed)%forward_endpoint
      topology(iseed)%backward_endpoint=qperp_results(iseed)%backward_endpoint
      topology(iseed)%forward_nstep=qperp_results(iseed)%forward_nstep
      topology(iseed)%backward_nstep=qperp_results(iseed)%backward_nstep
      topology(iseed)%forward_face=qperp_results(iseed)%forward_face
      topology(iseed)%backward_face=qperp_results(iseed)%backward_face
      topology(iseed)%forward_status=qperp_results(iseed)%forward_status
      topology(iseed)%backward_status=qperp_results(iseed)%backward_status
      topology(iseed)%map_forward_endpoint=qperp_results(iseed)%forward_endpoint
      topology(iseed)%map_backward_endpoint=qperp_results(iseed)%backward_endpoint
      topology(iseed)%map_forward_length=qperp_results(iseed)%forward_length
      topology(iseed)%map_backward_length=qperp_results(iseed)%backward_length
      topology(iseed)%map_forward_face=qperp_results(iseed)%forward_face
      topology(iseed)%map_backward_face=qperp_results(iseed)%backward_face
      topology(iseed)%map_forward_status=qperp_results(iseed)%forward_status
      topology(iseed)%map_backward_status=qperp_results(iseed)%backward_status
      topology(iseed)%source_B=0.d0
      topology(iseed)%forward_B=0.d0
      topology(iseed)%backward_B=0.d0
      topology(iseed)%source_B(1:ndim)=qperp_results(iseed)%B_seed
      topology(iseed)%forward_B(1:ndim)=qperp_results(iseed)%forward_B
      topology(iseed)%backward_B(1:ndim)=qperp_results(iseed)%backward_B
      topology(iseed)%source_Bn=0.d0
      topology(iseed)%forward_Bn=0.d0
      topology(iseed)%backward_Bn=0.d0
      topology(iseed)%has_twist=do_twist
      topology(iseed)%has_mapping=.false.
      topology(iseed)%valid_twist=.false.
      topology(iseed)%status_twist=qperp_results(iseed)%status
      if (do_twist) then
        topology(iseed)%valid_twist=twist_results(iseed)%valid_twist
        topology(iseed)%status_twist=twist_results(iseed)%status_twist
      endif
      topology(iseed)%valid=qperp_results(iseed)%status==trace_status_boundary
      topology(iseed)%status=qperp_results(iseed)%status
    enddo
  end subroutine mt_qperp_trace_to_topology

  subroutine mt_q0_trace_to_length(q_results,nseed,results)
    integer, intent(in) :: nseed
    type(trace_qperp_result), intent(in) :: q_results(nseed)
    type(trace_length_result), intent(out) :: results(nseed)

    integer :: iseed

    do iseed=1,nseed
      results(iseed)%seed=q_results(iseed)%seed
      results(iseed)%forward_footpoint=q_results(iseed)%forward_endpoint
      results(iseed)%backward_footpoint=q_results(iseed)%backward_endpoint
      results(iseed)%forward_length=q_results(iseed)%forward_length
      results(iseed)%backward_length=q_results(iseed)%backward_length
      results(iseed)%total_length=results(iseed)%forward_length &
           +results(iseed)%backward_length
      results(iseed)%forward_nstep=q_results(iseed)%forward_nstep
      results(iseed)%backward_nstep=q_results(iseed)%backward_nstep
      results(iseed)%forward_status=q_results(iseed)%forward_status
      results(iseed)%backward_status=q_results(iseed)%backward_status
    enddo
  end subroutine mt_q0_trace_to_length

  subroutine mt_q0_trace_to_topology(q_results,twist_results,nseed, &
       topology,do_twist)
    integer, intent(in) :: nseed
    type(trace_qperp_result), intent(in) :: q_results(nseed)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_topology_result), intent(out) :: topology(nseed)
    logical, intent(in) :: do_twist

    integer :: iseed

    do iseed=1,nseed
      topology(iseed)%seed=q_results(iseed)%seed
      topology(iseed)%length_forward=q_results(iseed)%forward_length
      topology(iseed)%length_backward=q_results(iseed)%backward_length
      topology(iseed)%length_total=topology(iseed)%length_forward &
           +topology(iseed)%length_backward
      topology(iseed)%twist_forward=0.d0
      topology(iseed)%twist_backward=0.d0
      topology(iseed)%twist_total=0.d0
      if (do_twist) then
        topology(iseed)%twist_forward=twist_results(iseed)%forward_twist
        topology(iseed)%twist_backward=twist_results(iseed)%backward_twist
        topology(iseed)%twist_total=twist_results(iseed)%total_twist
      endif
      topology(iseed)%forward_endpoint=q_results(iseed)%forward_endpoint
      topology(iseed)%backward_endpoint=q_results(iseed)%backward_endpoint
      topology(iseed)%forward_nstep=q_results(iseed)%forward_nstep
      topology(iseed)%backward_nstep=q_results(iseed)%backward_nstep
      topology(iseed)%forward_face=q_results(iseed)%forward_face
      topology(iseed)%backward_face=q_results(iseed)%backward_face
      topology(iseed)%forward_status=q_results(iseed)%forward_status
      topology(iseed)%backward_status=q_results(iseed)%backward_status
      topology(iseed)%map_forward_endpoint=q_results(iseed)%forward_endpoint
      topology(iseed)%map_backward_endpoint=q_results(iseed)%backward_endpoint
      topology(iseed)%map_forward_length=q_results(iseed)%forward_length
      topology(iseed)%map_backward_length=q_results(iseed)%backward_length
      topology(iseed)%map_forward_face=q_results(iseed)%forward_face
      topology(iseed)%map_backward_face=q_results(iseed)%backward_face
      topology(iseed)%map_forward_status=q_results(iseed)%forward_status
      topology(iseed)%map_backward_status=q_results(iseed)%backward_status
      topology(iseed)%source_B=0.d0
      topology(iseed)%forward_B=0.d0
      topology(iseed)%backward_B=0.d0
      topology(iseed)%source_B(1:ndim)=q_results(iseed)%B_seed
      topology(iseed)%forward_B(1:ndim)=q_results(iseed)%forward_B
      topology(iseed)%backward_B(1:ndim)=q_results(iseed)%backward_B
      topology(iseed)%source_Bn=0.d0
      topology(iseed)%forward_Bn=q_results(iseed)%forward_Bn_q0
      topology(iseed)%backward_Bn=q_results(iseed)%backward_Bn_q0
      topology(iseed)%has_twist=do_twist
      topology(iseed)%has_mapping=.false.
      topology(iseed)%valid_twist=.false.
      topology(iseed)%status_twist=q_results(iseed)%status
      if (do_twist) then
        topology(iseed)%valid_twist=twist_results(iseed)%valid_twist
        topology(iseed)%status_twist=twist_results(iseed)%status_twist
      endif
      topology(iseed)%valid=q_results(iseed)%status==trace_status_boundary
      topology(iseed)%status=q_results(iseed)%status
    enddo
  end subroutine mt_q0_trace_to_topology

  subroutine mt_topology_to_length(topology,nseed,results)
    integer, intent(in) :: nseed
    type(trace_topology_result), intent(in) :: topology(nseed)
    type(trace_length_result), intent(out) :: results(nseed)

    integer :: iseed

    do iseed=1,nseed
      results(iseed)%seed=topology(iseed)%seed
      results(iseed)%forward_footpoint=topology(iseed)%forward_endpoint
      results(iseed)%backward_footpoint=topology(iseed)%backward_endpoint
      results(iseed)%forward_length=topology(iseed)%length_forward
      results(iseed)%backward_length=topology(iseed)%length_backward
      results(iseed)%total_length=topology(iseed)%length_total
      results(iseed)%forward_nstep=topology(iseed)%forward_nstep
      results(iseed)%backward_nstep=topology(iseed)%backward_nstep
      results(iseed)%forward_status=topology(iseed)%forward_status
      results(iseed)%backward_status=topology(iseed)%backward_status
    enddo
  end subroutine mt_topology_to_length

  subroutine mt_topology_to_twist(topology,nseed,results)
    integer, intent(in) :: nseed
    type(trace_topology_result), intent(in) :: topology(nseed)
    type(trace_twist_result), intent(out) :: results(nseed)

    integer :: iseed

    do iseed=1,nseed
      results(iseed)%line%seed=topology(iseed)%seed
      results(iseed)%line%forward_footpoint=topology(iseed)%forward_endpoint
      results(iseed)%line%backward_footpoint=topology(iseed)%backward_endpoint
      results(iseed)%line%forward_length=topology(iseed)%length_forward
      results(iseed)%line%backward_length=topology(iseed)%length_backward
      results(iseed)%line%total_length=topology(iseed)%length_total
      results(iseed)%line%forward_nstep=topology(iseed)%forward_nstep
      results(iseed)%line%backward_nstep=topology(iseed)%backward_nstep
      results(iseed)%line%forward_status=topology(iseed)%forward_status
      results(iseed)%line%backward_status=topology(iseed)%backward_status
      results(iseed)%forward_twist=topology(iseed)%twist_forward
      results(iseed)%backward_twist=topology(iseed)%twist_backward
      results(iseed)%total_twist=topology(iseed)%twist_total
      results(iseed)%valid_twist=topology(iseed)%valid_twist
      results(iseed)%status_twist=topology(iseed)%status_twist
    enddo
  end subroutine mt_topology_to_twist

  subroutine mt_topology_to_mapping(topology,nseed,results)
    integer, intent(in) :: nseed
    type(trace_topology_result), intent(in) :: topology(nseed)
    type(trace_mapping_result), intent(out) :: results(nseed)

    integer :: iseed

    do iseed=1,nseed
      results(iseed)%seed=topology(iseed)%seed
      results(iseed)%source_B=topology(iseed)%source_B
      results(iseed)%forward_footpoint=topology(iseed)%map_forward_endpoint
      results(iseed)%backward_footpoint=topology(iseed)%map_backward_endpoint
      results(iseed)%forward_B=topology(iseed)%forward_B
      results(iseed)%backward_B=topology(iseed)%backward_B
      results(iseed)%forward_length=topology(iseed)%map_forward_length
      results(iseed)%backward_length=topology(iseed)%map_backward_length
      results(iseed)%source_Bn=topology(iseed)%source_Bn
      results(iseed)%forward_Bn=topology(iseed)%forward_Bn
      results(iseed)%backward_Bn=topology(iseed)%backward_Bn
      results(iseed)%forward_face=topology(iseed)%map_forward_face
      results(iseed)%backward_face=topology(iseed)%map_backward_face
      results(iseed)%forward_status=topology(iseed)%map_forward_status
      results(iseed)%backward_status=topology(iseed)%map_backward_status
      results(iseed)%valid=topology(iseed)%has_mapping .and. &
           topology(iseed)%valid
    enddo
  end subroutine mt_topology_to_mapping

  subroutine mt_validate_axis_plane(c1min,c1max,n1,c2min,c2max,n2,caller)
    integer, intent(in) :: n1,n2
    double precision, intent(in) :: c1min,c1max,c2min,c2max
    character(len=*), intent(in) :: caller

    if (npe/=1) then
      call mpistop(trim(caller)//' currently requires npe=1')
    endif
    if (ndim/=3) then
      call mpistop(trim(caller)//' requires ndim=3')
      return
    endif
    {^IFTHREED
    select case (geo_coordinate)
    case (geo_cartesian)
      ! Coordinate-plane products are sampled seed surfaces. RK2 tracing can
      ! use either slab-uniform or AMR Cartesian simulation grids.
    case (geo_cartesian_stretched)
    case default
      call mpistop(trim(caller)//' requires Cartesian geometry')
    end select
    }
    if (n1<1 .or. n2<1) then
      call mpistop(trim(caller)//' requires both sample counts >=1')
    endif
    if (c1max<c1min .or. c2max<c2min) then
      call mpistop(trim(caller)//' requires ordered plane bounds')
    endif
  end subroutine mt_validate_axis_plane

  subroutine mt_build_axis_plane_seeds(c1min,c1max,n1,c2min,c2max,n2, &
       fixed_value,axis1,axis2,fixed_axis,seeds)
    integer, intent(in) :: n1,n2,axis1,axis2,fixed_axis
    double precision, intent(in) :: c1min,c1max,c2min,c2max,fixed_value
    double precision, allocatable, intent(out) :: seeds(:,:)

    double precision :: dc1,dc2
    integer :: i,j,iseed

    dc1=0.d0
    dc2=0.d0
    if (n1>1) dc1=(c1max-c1min)/dble(n1-1)
    if (n2>1) dc2=(c2max-c2min)/dble(n2-1)

    allocate(seeds(n1*n2,ndim))
    seeds=0.d0
    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        seeds(iseed,axis1)=c1min+dble(i-1)*dc1
        seeds(iseed,axis2)=c2min+dble(j-1)*dc2
        seeds(iseed,fixed_axis)=fixed_value
      enddo
    enddo
  end subroutine mt_build_axis_plane_seeds

  subroutine mt_build_spherical_surface_seeds(surface,seed_coord,s1_min, &
       s1_max,n1,s2_min,s2_max,n2,seed_layout,seed_theta0,seed_phi0, &
       seed_alpha,seeds)
    character(len=*), intent(in) :: surface,seed_layout
    integer, intent(in) :: n1,n2
    double precision, intent(in) :: seed_coord,s1_min,s1_max,s2_min,s2_max
    double precision, intent(in) :: seed_theta0,seed_phi0,seed_alpha
    double precision, allocatable, intent(out) :: seeds(:,:)

    character(len=mt_task_name_len) :: layout
    double precision :: s1,s2
    integer :: i,j,iseed,status

    layout=mt_lowercase(trim(seed_layout))
    if (len_trim(layout)==0) layout='endpoint'
    select case (trim(layout))
    case ('endpoint','endpoints','boundary')
    case ('cell_centered','cell-centered','centered','center')
    case default
      call mpistop('spherical_surface_products requires '// &
           'mt_seed_layout=endpoint or cell_centered')
    end select

    allocate(seeds(n1*n2,ndim))
    seeds=0.d0
    do j=1,n2
      s2=mt_seed_axis_coord(s2_min,s2_max,n2,j,layout)
      do i=1,n1
        s1=mt_seed_axis_coord(s1_min,s1_max,n1,i,layout)
        iseed=(j-1)*n1+i
        select case (trim(surface))
        case ('rmin','rconst')
          seeds(iseed,1)=seed_coord
          seeds(iseed,2)=s1
          seeds(iseed,3)=s2
        case ('theta_const')
          seeds(iseed,1)=s1
          seeds(iseed,2)=seed_coord
          seeds(iseed,3)=s2
        case ('phi_const')
          seeds(iseed,1)=s1
          seeds(iseed,2)=s2
          seeds(iseed,3)=seed_coord
        case ('radial_plane')
          call mt_spherical_radial_plane_seed(s1,s2,seed_theta0,seed_phi0, &
               seed_alpha,seeds(iseed,:),status)
          if (status/=0) then
            call mpistop('mt_build_spherical_surface_seeds radial_plane '// &
                 'seed lies outside the spherical domain; adjust '// &
                 'mt_seed_theta0/mt_seed_phi0/mt_seed_alpha or mt_s2 bounds')
          endif
        case default
          call mpistop('mt_build_spherical_surface_seeds got '// &
               'unsupported surface')
        end select
      enddo
    enddo
  end subroutine mt_build_spherical_surface_seeds

  subroutine mt_spherical_radial_plane_seed(radius,u,theta0,phi0,alpha,seed, &
       status)
    double precision, intent(in) :: radius,u,theta0,phi0,alpha
    double precision, intent(out) :: seed(ndim)
    integer, intent(out) :: status

    double precision :: sin_theta0,cos_theta0,sin_phi0,cos_phi0
    double precision :: n0(3),etheta(3),ephi(3),tangent(3),nhat(3)
    double precision :: norm_n,theta,phi,cos_u,sin_u

    seed=0.d0
    status=1
    {^IFTHREED
    if (radius<xprobmin1 .or. radius>xprobmax1) return
    if (theta0<xprobmin2 .or. theta0>xprobmax2) return
    if (phi0<xprobmin3 .or. phi0>xprobmax3) return

    sin_theta0=dsin(theta0)
    cos_theta0=dcos(theta0)
    if (abs(sin_theta0)<=1.d-12) return
    sin_phi0=dsin(phi0)
    cos_phi0=dcos(phi0)

    n0=(/sin_theta0*cos_phi0,sin_theta0*sin_phi0,cos_theta0/)
    etheta=(/cos_theta0*cos_phi0,cos_theta0*sin_phi0,-sin_theta0/)
    ephi=(/-sin_phi0,cos_phi0,0.d0/)
    tangent=dcos(alpha)*etheta+dsin(alpha)*ephi

    cos_u=dcos(u)
    sin_u=dsin(u)
    nhat=cos_u*n0+sin_u*tangent
    norm_n=dsqrt(sum(nhat**2))
    if (norm_n<=0.d0) return
    nhat=nhat/norm_n

    theta=dacos(max(-1.d0,min(1.d0,nhat(3))))
    phi=datan2(nhat(2),nhat(1))
    if (theta<xprobmin2 .or. theta>xprobmax2) return
    if (phi<xprobmin3 .or. phi>xprobmax3) return

    seed(1)=radius
    seed(2)=theta
    seed(3)=phi
    status=0
    }
  end subroutine mt_spherical_radial_plane_seed

  subroutine mt_build_spherical_cloud_seeds(s1_min,s1_max,n1,s2_min, &
       s2_max,n2,s3_min,s3_max,n3,seed_layout,seeds)
    character(len=*), intent(in) :: seed_layout
    integer, intent(in) :: n1,n2,n3
    double precision, intent(in) :: s1_min,s1_max,s2_min,s2_max
    double precision, intent(in) :: s3_min,s3_max
    double precision, allocatable, intent(out) :: seeds(:,:)

    character(len=mt_task_name_len) :: layout
    double precision :: s1,s2,s3
    integer :: i,j,k,iseed

    layout=mt_lowercase(trim(seed_layout))
    if (len_trim(layout)==0) layout='endpoint'
    select case (trim(layout))
    case ('endpoint','endpoints','boundary')
    case ('cell_centered','cell-centered','centered','center')
    case default
      call mpistop('spherical_cloud_products requires '// &
           'mt_seed_layout=endpoint or cell_centered')
    end select

    allocate(seeds(n1*n2*n3,ndim))
    seeds=0.d0
    do k=1,n3
      s3=mt_seed_axis_coord(s3_min,s3_max,n3,k,layout)
      do j=1,n2
        s2=mt_seed_axis_coord(s2_min,s2_max,n2,j,layout)
        do i=1,n1
          s1=mt_seed_axis_coord(s1_min,s1_max,n1,i,layout)
          iseed=(k-1)*n1*n2+(j-1)*n1+i
          seeds(iseed,1)=s1
          seeds(iseed,2)=s2
          seeds(iseed,3)=s3
        enddo
      enddo
    enddo
  end subroutine mt_build_spherical_cloud_seeds

  double precision function mt_seed_axis_coord(cmin,cmax,n,i,layout) &
       result(coord)
    double precision, intent(in) :: cmin,cmax
    integer, intent(in) :: n,i
    character(len=*), intent(in) :: layout

    select case (trim(layout))
    case ('cell_centered','cell-centered','centered','center')
      coord=cmin+(dble(i)-0.5d0)*(cmax-cmin)/dble(n)
    case default
      if (n>1) then
        coord=cmin+dble(i-1)*(cmax-cmin)/dble(n-1)
      else
        coord=cmin
      endif
    end select
  end function mt_seed_axis_coord

  integer function mt_spherical_connection_type(face_b,face_f,is_valid) &
       result(connection_type)
    integer, intent(in) :: face_b,face_f
    logical, intent(in) :: is_valid

    logical :: b_rmin,f_rmin,b_rmax,f_rmax,b_side,f_side

    connection_type=0
    if (.not.is_valid) return
    if (.not.mt_q_face_valid(face_b)) return
    if (.not.mt_q_face_valid(face_f)) return

    b_rmin=face_b==trace_face_xmin
    f_rmin=face_f==trace_face_xmin
    b_rmax=face_b==trace_face_xmax
    f_rmax=face_f==trace_face_xmax
    b_side=.not.(b_rmin .or. b_rmax)
    f_side=.not.(f_rmin .or. f_rmax)

    if (b_rmin .and. f_rmin) then
      connection_type=1
    else if ((b_rmin .and. f_rmax) .or. (b_rmax .and. f_rmin)) then
      connection_type=2
    else if ((b_rmin .and. f_side) .or. (f_rmin .and. b_side)) then
      connection_type=3
    else if (b_rmax .and. f_rmax) then
      connection_type=4
    else if (b_side .and. f_side .and. face_b==face_f) then
      connection_type=5
    else
      connection_type=6
    endif
  end function mt_spherical_connection_type

  logical function mt_validate_arbitrary_plane_basis(e1,e2,s1min,s1max, &
       n1,s2min,s2max,n2,caller) result(valid)
    integer, intent(in) :: n1,n2
    double precision, intent(in) :: e1(ndim),e2(ndim)
    double precision, intent(in) :: s1min,s1max,s2min,s2max
    character(len=*), intent(in) :: caller

    double precision, parameter :: basis_tol=1.d-10
    double precision :: e1_norm,e2_norm,e12_dot

    valid=.false.
    if (npe/=1) then
      write(*,'(a)') trim(caller)//' currently requires npe=1'
      return
    endif
    if (ndim/=3) then
      write(*,'(a)') trim(caller)//' requires ndim=3'
      return
    endif
    {^IFTHREED
    select case (geo_coordinate)
    case (geo_cartesian,geo_cartesian_stretched)
    case default
      write(*,'(a)') trim(caller)//' requires Cartesian geometry'
      return
    end select
    }
    if (n1<1 .or. n2<1) then
      write(*,'(a)') trim(caller)//' requires both sample counts >=1'
      return
    endif
    if (s1max<s1min .or. s2max<s2min) then
      write(*,'(a)') trim(caller)//' requires ordered plane bounds'
      return
    endif
    e1_norm=dsqrt(sum(e1**2))
    e2_norm=dsqrt(sum(e2**2))
    e12_dot=sum(e1*e2)
    if (abs(e1_norm-1.d0)>basis_tol .or. &
         abs(e2_norm-1.d0)>basis_tol .or. &
         abs(e12_dot)>basis_tol) then
      write(*,'(a)') trim(caller)//' requires orthonormal e1/e2'
      return
    endif

    valid=.true.
  end function mt_validate_arbitrary_plane_basis

  subroutine mt_build_arbitrary_plane_seeds(origin,e1,e2,s1min,s1max,n1, &
       s2min,s2max,n2,seeds,s1,s2)
    integer, intent(in) :: n1,n2
    double precision, intent(in) :: origin(ndim),e1(ndim),e2(ndim)
    double precision, intent(in) :: s1min,s1max,s2min,s2max
    double precision, allocatable, intent(out) :: seeds(:,:),s1(:),s2(:)

    double precision :: ds1,ds2
    integer :: i,j,iseed

    ds1=0.d0
    ds2=0.d0
    if (n1>1) ds1=(s1max-s1min)/dble(n1-1)
    if (n2>1) ds2=(s2max-s2min)/dble(n2-1)

    allocate(seeds(n1*n2,ndim),s1(n1*n2),s2(n1*n2))
    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        s1(iseed)=s1min+dble(i-1)*ds1
        s2(iseed)=s2min+dble(j-1)*ds2
        seeds(iseed,:)=origin+s1(iseed)*e1+s2(iseed)*e2
      enddo
    enddo
  end subroutine mt_build_arbitrary_plane_seeds

  subroutine mt_write_length_plane_csv(results,n1,n2,csv_file,caller, &
       index_header)
    integer, intent(in) :: n1,n2
    type(trace_length_result), intent(in) :: results(n1*n2)
    character(len=*), intent(in) :: csv_file,caller,index_header

    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         trim(index_header)//',seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        seed_xyz=0.d0
        seed_xyz(1:ndim)=results(iseed)%seed
        write(csv_unit,'(i0,",",i0,6(",",es24.16),4(",",i0))', &
             iostat=io_status) i,j,seed_xyz,results(iseed)%total_length, &
             results(iseed)%backward_length, &
             results(iseed)%forward_length, &
             results(iseed)%backward_nstep, &
             results(iseed)%forward_nstep, &
             results(iseed)%backward_status, &
             results(iseed)%forward_status
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_length_plane_csv

  subroutine mt_write_spherical_rmin_csv(topology,n1,n2,csv_file,caller, &
       do_twist,do_q,q_results,do_qperp,qperp_results)
    integer, intent(in) :: n1,n2
    type(trace_topology_result), intent(in) :: topology(n1*n2)
    character(len=*), intent(in) :: csv_file,caller
    logical, intent(in) :: do_twist,do_q,do_qperp
    type(trace_qperp_result), intent(in) :: q_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)

    integer :: csv_unit,io_status,i,j,iseed,connection_type
    logical :: valid
    character(len=2048) :: header

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    header='i,j,seed_r,seed_theta,seed_phi,'// &
         'length_total,length_backward,length_forward'
    if (do_twist) then
      header=trim(header)//','// &
           'twist_total,twist_backward,twist_forward'
    endif
    if (do_q) then
      header=trim(header)//','// &
           'logQ,valid_Q,status_Q'
    endif
    if (do_qperp) then
      header=trim(header)//','// &
           'logQperp,valid_Qperp,status_Qperp'
    endif
    header=trim(header)//','// &
         'r_b,theta_b,phi_b,r_f,theta_f,phi_f,'// &
         'face_backward,face_forward,connection_type_spherical,'// &
         'nstep_backward,nstep_forward,status_backward,status_forward'
    if (do_twist) then
      header=trim(header)//','// &
           'status_twist,valid_twist'
    endif
    header=trim(header)//',valid'
    write(csv_unit,'(a)',iostat=io_status) trim(header)
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        valid=topology(iseed)%valid
        connection_type=mt_spherical_connection_type( &
             topology(iseed)%backward_face,topology(iseed)%forward_face, &
             valid)
        write(csv_unit,'(i0,",",i0,6(",",es24.16))',advance='no', &
             iostat=io_status) i,j,topology(iseed)%seed, &
             topology(iseed)%length_total, &
             topology(iseed)%length_backward, &
             topology(iseed)%length_forward
        if (io_status/=0) exit
        if (do_twist) then
          write(csv_unit,'(3(",",es24.16))',advance='no', &
               iostat=io_status) topology(iseed)%twist_total, &
               topology(iseed)%twist_backward, &
               topology(iseed)%twist_forward
          if (io_status/=0) exit
        endif
        if (do_q) then
          write(csv_unit,'(",",es24.16,",",l1,",",i0)',advance='no', &
               iostat=io_status) q_results(iseed)%logq0, &
               q_results(iseed)%valid_q0,q_results(iseed)%status_q0
          if (io_status/=0) exit
        endif
        if (do_qperp) then
          write(csv_unit,'(",",es24.16,",",l1,",",i0)',advance='no', &
               iostat=io_status) qperp_results(iseed)%logqperp, &
               qperp_results(iseed)%valid,qperp_results(iseed)%status
          if (io_status/=0) exit
        endif
        write(csv_unit,'(6(",",es24.16),7(",",i0))',advance='no', &
             iostat=io_status) topology(iseed)%backward_endpoint, &
             topology(iseed)%forward_endpoint, &
             topology(iseed)%backward_face,topology(iseed)%forward_face, &
             connection_type,topology(iseed)%backward_nstep, &
             topology(iseed)%forward_nstep, &
             topology(iseed)%backward_status, &
             topology(iseed)%forward_status
        if (io_status/=0) exit
        if (do_twist) then
          write(csv_unit,'(",",i0,",",l1)',advance='no',iostat=io_status) &
               topology(iseed)%status_twist,topology(iseed)%valid_twist
          if (io_status/=0) exit
        endif
        write(csv_unit,'(",",l1)',iostat=io_status) valid
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_spherical_rmin_csv

  subroutine mt_write_spherical_topology_vtu(vtu_file,topology,n1,n2, &
       caller,do_twist,do_q,q_results,do_qperp,qperp_results)
    character(len=*), intent(in) :: vtu_file,caller
    integer, intent(in) :: n1,n2
    type(trace_topology_result), intent(in) :: topology(n1*n2)
    logical, intent(in) :: do_twist,do_q,do_qperp
    type(trace_qperp_result), intent(in) :: q_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)

    integer :: vtu_unit,io_status,npoint,ncell

    npoint=n1*n2
    ncell=0
    if (n1>1 .and. n2>1) ncell=(n1-1)*(n2-1)
    if (ncell==0) ncell=npoint

    open(newunit=vtu_unit,file=trim(vtu_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTU file')
    endif

    call mt_write_vtu_file_header(vtu_unit,npoint,ncell)
    call mt_write_vtu_spherical_topology_pointdata(vtu_unit,topology, &
         npoint,do_twist,do_q,q_results,do_qperp,qperp_results,io_status)
    if (io_status==0) call mt_write_vtu_topology_points(vtu_unit, &
         topology,npoint,io_status)
    if (io_status==0) then
      if (n1>1 .and. n2>1) then
        call mt_write_vtu_quad_cells(vtu_unit,n1,n2,io_status)
      else
        call mt_write_vtu_vertex_cells(vtu_unit,npoint,io_status)
      endif
    endif
    if (io_status==0) call mt_write_vtu_file_footer(vtu_unit,io_status)
    close(vtu_unit)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not write VTU file')
    endif
  end subroutine mt_write_spherical_topology_vtu

  subroutine mt_write_spherical_cloud_csv(topology,n1,n2,n3,csv_file, &
       caller,do_twist,do_qperp,qperp_results)
    integer, intent(in) :: n1,n2,n3
    type(trace_topology_result), intent(in) :: topology(n1*n2*n3)
    character(len=*), intent(in) :: csv_file,caller
    logical, intent(in) :: do_twist,do_qperp
    type(trace_qperp_result), intent(in) :: qperp_results(:)

    integer :: csv_unit,io_status,i,j,k,iseed,connection_type
    logical :: valid
    character(len=2048) :: header

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    header='i,j,k,seed_r,seed_theta,seed_phi,'// &
         'length_total,length_backward,length_forward'
    if (do_twist) then
      header=trim(header)//','// &
           'twist_total,twist_backward,twist_forward'
    endif
    if (do_qperp) then
      header=trim(header)//','// &
           'logQperp,valid_Qperp,status_Qperp'
    endif
    header=trim(header)//','// &
         'r_b,theta_b,phi_b,r_f,theta_f,phi_f,'// &
         'face_backward,face_forward,connection_type_spherical,'// &
         'nstep_backward,nstep_forward,status_backward,status_forward'
    if (do_twist) then
      header=trim(header)//','// &
           'status_twist,valid_twist'
    endif
    header=trim(header)//',valid'
    write(csv_unit,'(a)',iostat=io_status) trim(header)
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do k=1,n3
      do j=1,n2
        do i=1,n1
          iseed=(k-1)*n1*n2+(j-1)*n1+i
          valid=topology(iseed)%valid
          connection_type=mt_spherical_connection_type( &
               topology(iseed)%backward_face,topology(iseed)%forward_face, &
               valid)
          write(csv_unit,'(i0,",",i0,",",i0,6(",",es24.16))', &
               advance='no',iostat=io_status) i,j,k,topology(iseed)%seed, &
               topology(iseed)%length_total, &
               topology(iseed)%length_backward, &
               topology(iseed)%length_forward
          if (io_status/=0) exit
          if (do_twist) then
            write(csv_unit,'(3(",",es24.16))',advance='no', &
                 iostat=io_status) topology(iseed)%twist_total, &
                 topology(iseed)%twist_backward, &
                 topology(iseed)%twist_forward
            if (io_status/=0) exit
          endif
          if (do_qperp) then
            write(csv_unit,'(",",es24.16,",",l1,",",i0)',advance='no', &
                 iostat=io_status) qperp_results(iseed)%logqperp, &
                 qperp_results(iseed)%valid,qperp_results(iseed)%status
            if (io_status/=0) exit
          endif
          write(csv_unit,'(6(",",es24.16),7(",",i0))',advance='no', &
               iostat=io_status) topology(iseed)%backward_endpoint, &
               topology(iseed)%forward_endpoint, &
               topology(iseed)%backward_face,topology(iseed)%forward_face, &
               connection_type,topology(iseed)%backward_nstep, &
               topology(iseed)%forward_nstep, &
               topology(iseed)%backward_status, &
               topology(iseed)%forward_status
          if (io_status/=0) exit
          if (do_twist) then
            write(csv_unit,'(",",i0,",",l1)',advance='no', &
                 iostat=io_status) topology(iseed)%status_twist, &
                 topology(iseed)%valid_twist
            if (io_status/=0) exit
          endif
          write(csv_unit,'(",",l1)',iostat=io_status) valid
          if (io_status/=0) exit
        enddo
        if (io_status/=0) exit
      enddo
      if (io_status/=0) exit
    enddo

    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV data')
    endif
    close(csv_unit)
  end subroutine mt_write_spherical_cloud_csv

  subroutine mt_write_spherical_cloud_vtu(vtu_file,topology,npoint,caller, &
       do_twist,do_q,q_results,do_qperp,qperp_results)
    character(len=*), intent(in) :: vtu_file,caller
    integer, intent(in) :: npoint
    type(trace_topology_result), intent(in) :: topology(npoint)
    logical, intent(in) :: do_twist,do_q,do_qperp
    type(trace_qperp_result), intent(in) :: q_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)

    integer :: vtu_unit,io_status

    open(newunit=vtu_unit,file=trim(vtu_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTU file')
    endif

    call mt_write_vtu_file_header(vtu_unit,npoint,npoint)
    call mt_write_vtu_spherical_topology_pointdata(vtu_unit,topology, &
         npoint,do_twist,do_q,q_results,do_qperp,qperp_results,io_status)
    if (io_status==0) call mt_write_vtu_topology_points(vtu_unit, &
         topology,npoint,io_status)
    if (io_status==0) call mt_write_vtu_vertex_cells(vtu_unit,npoint, &
         io_status)
    if (io_status==0) call mt_write_vtu_file_footer(vtu_unit,io_status)
    close(vtu_unit)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not write VTU file')
    endif
  end subroutine mt_write_spherical_cloud_vtu

  subroutine mt_write_fieldline_products_seeds_csv(length_results, &
       twist_results,q_results,qperp_results,nseed,csv_file,do_twist, &
       do_q,do_qperp)
    integer, intent(in) :: nseed
    type(trace_length_result), intent(in) :: length_results(nseed)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: q_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    character(len=*), intent(in) :: csv_file
    logical, intent(in) :: do_twist,do_q,do_qperp

    double precision :: seed_xyz(3)
    double precision :: Bseed_norm,Bf_norm,Bb_norm
    character(len=2048) :: header
    integer :: csv_unit,io_status,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop('mt_fieldline_products_seeds could not open CSV file')
    endif

    header='seed_id,seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'nstep_backward_length,nstep_forward_length,'// &
         'status_backward_length,status_forward_length'
    if (do_twist) then
      header=trim(header)//','// &
           'twist_total,twist_backward,twist_forward,'// &
           'nstep_backward_twist,nstep_forward_twist,'// &
           'status_backward_twist,status_forward_twist'
    endif
    if (do_q) then
      header=trim(header)//','// &
           'logQ,valid_Q,status_Q,'// &
           'face_forward_Q,face_backward_Q,'// &
           'status_forward_Q,status_backward_Q,'// &
           'length_forward_Q,length_backward_Q,'// &
           'x_f_Q,y_f_Q,z_f_Q,x_b_Q,y_b_Q,z_b_Q'
    endif
    if (do_qperp) then
      header=trim(header)//','// &
           'qperp,logqperp,valid_qperp,status_qperp,'// &
           'N2,bfactor,'// &
           'face_forward_qperp,face_backward_qperp,'// &
           'status_forward_qperp,status_backward_qperp,'// &
           'length_forward_qperp,length_backward_qperp,'// &
           'Bseed_norm,Bf_norm,Bb_norm,'// &
           'x_f_qperp,y_f_qperp,z_f_qperp,'// &
           'x_b_qperp,y_b_qperp,z_b_qperp'
    endif
    write(csv_unit,'(a)',iostat=io_status) trim(header)
    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_fieldline_products_seeds could not write CSV header')
    endif

    do iseed=1,nseed
      seed_xyz=0.d0
      seed_xyz(1:ndim)=length_results(iseed)%seed
      write(csv_unit,'(i0,6(",",es24.16),4(",",i0))', &
           advance='no',iostat=io_status) &
           iseed,seed_xyz,length_results(iseed)%total_length, &
           length_results(iseed)%backward_length, &
           length_results(iseed)%forward_length, &
           length_results(iseed)%backward_nstep, &
           length_results(iseed)%forward_nstep, &
           length_results(iseed)%backward_status, &
           length_results(iseed)%forward_status
      if (io_status/=0) exit

      if (do_twist) then
        write(csv_unit,'(3(",",es24.16),4(",",i0))', &
             advance='no',iostat=io_status) &
             twist_results(iseed)%total_twist, &
             twist_results(iseed)%backward_twist, &
             twist_results(iseed)%forward_twist, &
             twist_results(iseed)%line%backward_nstep, &
             twist_results(iseed)%line%forward_nstep, &
             twist_results(iseed)%line%backward_status, &
             twist_results(iseed)%line%forward_status
        if (io_status/=0) exit
      endif

      if (do_q) then
        write(csv_unit, &
             '(",",es24.16,",",l1,",",i0,4(",",i0),8(",",es24.16))', &
             advance='no',iostat=io_status) &
             q_results(iseed)%logq0, &
             q_results(iseed)%valid_q0,q_results(iseed)%status_q0, &
             q_results(iseed)%forward_face, &
             q_results(iseed)%backward_face, &
             q_results(iseed)%forward_status, &
             q_results(iseed)%backward_status, &
             q_results(iseed)%forward_length, &
             q_results(iseed)%backward_length, &
             q_results(iseed)%forward_endpoint, &
             q_results(iseed)%backward_endpoint
        if (io_status/=0) exit
      endif

      if (do_qperp) then
        Bseed_norm=dsqrt(sum(qperp_results(iseed)%B_seed**2))
        Bf_norm=dsqrt(sum(qperp_results(iseed)%forward_B**2))
        Bb_norm=dsqrt(sum(qperp_results(iseed)%backward_B**2))
        write(csv_unit, &
             '(2(",",es24.16),",",l1,",",i0,2(",",es24.16),'// &
             '4(",",i0),11(",",es24.16))', &
             advance='no',iostat=io_status) &
             qperp_results(iseed)%qperp, &
             qperp_results(iseed)%logqperp, &
             qperp_results(iseed)%valid,qperp_results(iseed)%status, &
             qperp_results(iseed)%N2,qperp_results(iseed)%bfactor, &
             qperp_results(iseed)%forward_face, &
             qperp_results(iseed)%backward_face, &
             qperp_results(iseed)%forward_status, &
             qperp_results(iseed)%backward_status, &
             qperp_results(iseed)%forward_length, &
             qperp_results(iseed)%backward_length, &
             Bseed_norm,Bf_norm,Bb_norm, &
             qperp_results(iseed)%forward_endpoint, &
             qperp_results(iseed)%backward_endpoint
        if (io_status/=0) exit
      endif

      write(csv_unit,'()',iostat=io_status)
      if (io_status/=0) exit
    enddo

    if (io_status/=0) then
      close(csv_unit)
      call mpistop('mt_fieldline_products_seeds could not write CSV data')
    endif
    close(csv_unit)
  end subroutine mt_write_fieldline_products_seeds_csv

  subroutine mt_fieldline_products_append_header(header,do_twist,do_q,do_qperp)
    character(len=*), intent(inout) :: header
    logical, intent(in) :: do_twist,do_q,do_qperp

    header=trim(header)//','// &
         'length_total,length_backward,length_forward,'// &
         'nstep_backward_length,nstep_forward_length,'// &
         'status_backward_length,status_forward_length'
    if (do_twist) then
      header=trim(header)//','// &
           'twist_total,twist_backward,twist_forward,'// &
           'nstep_backward_twist,nstep_forward_twist,'// &
           'status_backward_twist,status_forward_twist'
    endif
    if (do_q) then
      header=trim(header)//','// &
           'logQ,valid_Q,status_Q,'// &
           'face_forward_Q,face_backward_Q,'// &
           'status_forward_Q,status_backward_Q,'// &
           'length_forward_Q,length_backward_Q,'// &
           'x_f_Q,y_f_Q,z_f_Q,x_b_Q,y_b_Q,z_b_Q'
    endif
    if (do_qperp) then
      header=trim(header)//','// &
           'qperp,logqperp,valid_qperp,status_qperp,'// &
           'N2,bfactor,'// &
           'face_forward_qperp,face_backward_qperp,'// &
           'status_forward_qperp,status_backward_qperp,'// &
           'length_forward_qperp,length_backward_qperp,'// &
           'Bseed_norm,Bf_norm,Bb_norm,'// &
           'x_f_qperp,y_f_qperp,z_f_qperp,'// &
           'x_b_qperp,y_b_qperp,z_b_qperp'
    endif
  end subroutine mt_fieldline_products_append_header

  subroutine mt_write_fieldline_products_plane_arbitrary_header(csv_file, &
       caller,do_twist,do_q,do_qperp)
    character(len=*), intent(in) :: csv_file,caller
    logical, intent(in) :: do_twist,do_q,do_qperp

    character(len=2048) :: header
    integer :: csv_unit,io_status

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    header='i,j,s1,s2,seed_x,seed_y,seed_z'
    call mt_fieldline_products_append_header(header,do_twist,do_q,do_qperp)
    write(csv_unit,'(a)',iostat=io_status) trim(header)
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    close(csv_unit)
  end subroutine mt_write_fieldline_products_plane_arbitrary_header

  subroutine mt_write_fieldline_products_plane_arbitrary_csv( &
       length_results,twist_results,q_results,qperp_results,s1,s2,n1,n2, &
       csv_file,do_twist,do_q,do_qperp,caller)
    integer, intent(in) :: n1,n2
    type(trace_length_result), intent(in) :: length_results(n1*n2)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: q_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    double precision, intent(in) :: s1(n1*n2),s2(n1*n2)
    character(len=*), intent(in) :: csv_file,caller
    logical, intent(in) :: do_twist,do_q,do_qperp

    double precision :: seed_xyz(3)
    double precision :: Bseed_norm,Bf_norm,Bb_norm
    character(len=2048) :: header
    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    header='i,j,s1,s2,seed_x,seed_y,seed_z'
    call mt_fieldline_products_append_header(header,do_twist,do_q,do_qperp)
    write(csv_unit,'(a)',iostat=io_status) trim(header)
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        seed_xyz=0.d0
        seed_xyz(1:ndim)=length_results(iseed)%seed
        write(csv_unit,'(i0,",",i0,8(",",es24.16),4(",",i0))', &
             advance='no',iostat=io_status) &
             i,j,s1(iseed),s2(iseed),seed_xyz, &
             length_results(iseed)%total_length, &
             length_results(iseed)%backward_length, &
             length_results(iseed)%forward_length, &
             length_results(iseed)%backward_nstep, &
             length_results(iseed)%forward_nstep, &
             length_results(iseed)%backward_status, &
             length_results(iseed)%forward_status
        if (io_status/=0) exit

        if (do_twist) then
          write(csv_unit,'(3(",",es24.16),4(",",i0))', &
               advance='no',iostat=io_status) &
               twist_results(iseed)%total_twist, &
               twist_results(iseed)%backward_twist, &
               twist_results(iseed)%forward_twist, &
               twist_results(iseed)%line%backward_nstep, &
               twist_results(iseed)%line%forward_nstep, &
               twist_results(iseed)%line%backward_status, &
               twist_results(iseed)%line%forward_status
          if (io_status/=0) exit
        endif

        if (do_q) then
          write(csv_unit, &
               '(",",es24.16,",",l1,",",i0,4(",",i0),8(",",es24.16))', &
               advance='no',iostat=io_status) &
               q_results(iseed)%logq0, &
               q_results(iseed)%valid_q0,q_results(iseed)%status_q0, &
               q_results(iseed)%forward_face, &
               q_results(iseed)%backward_face, &
               q_results(iseed)%forward_status, &
               q_results(iseed)%backward_status, &
               q_results(iseed)%forward_length, &
               q_results(iseed)%backward_length, &
               q_results(iseed)%forward_endpoint, &
               q_results(iseed)%backward_endpoint
          if (io_status/=0) exit
        endif

        if (do_qperp) then
          Bseed_norm=dsqrt(sum(qperp_results(iseed)%B_seed**2))
          Bf_norm=dsqrt(sum(qperp_results(iseed)%forward_B**2))
          Bb_norm=dsqrt(sum(qperp_results(iseed)%backward_B**2))
          write(csv_unit, &
               '(2(",",es24.16),",",l1,",",i0,2(",",es24.16),'// &
               '4(",",i0),11(",",es24.16))', &
               advance='no',iostat=io_status) &
               qperp_results(iseed)%qperp, &
               qperp_results(iseed)%logqperp, &
               qperp_results(iseed)%valid,qperp_results(iseed)%status, &
               qperp_results(iseed)%N2,qperp_results(iseed)%bfactor, &
               qperp_results(iseed)%forward_face, &
               qperp_results(iseed)%backward_face, &
               qperp_results(iseed)%forward_status, &
               qperp_results(iseed)%backward_status, &
               qperp_results(iseed)%forward_length, &
               qperp_results(iseed)%backward_length, &
               Bseed_norm,Bf_norm,Bb_norm, &
               qperp_results(iseed)%forward_endpoint, &
               qperp_results(iseed)%backward_endpoint
          if (io_status/=0) exit
        endif

        write(csv_unit,'()',iostat=io_status)
        if (io_status/=0) exit
      enddo
      if (io_status/=0) exit
    enddo

    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV data')
    endif
    close(csv_unit)
  end subroutine mt_write_fieldline_products_plane_arbitrary_csv

  subroutine mt_write_fieldline_products_vtu_empty(vtu_file,do_twist, &
       do_q,do_qperp,caller,do_length)
    character(len=*), intent(in) :: vtu_file,caller
    logical, intent(in) :: do_twist,do_q,do_qperp
    logical, intent(in), optional :: do_length

    integer :: vtu_unit,io_status
    logical :: do_length_eff

    open(newunit=vtu_unit,file=trim(vtu_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTU file')
    endif

    call mt_write_vtu_file_header(vtu_unit,0,0)
    do_length_eff=.true.
    if (present(do_length)) do_length_eff=do_length
    call mt_write_vtu_empty_pointdata(vtu_unit,do_length_eff,do_twist, &
         do_q,do_qperp)
    write(vtu_unit,'(a)',iostat=io_status) '<CellData>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</CellData>'
    if (io_status==0) then
      write(vtu_unit,'(a)',iostat=io_status) '<Points>'
    endif
    if (io_status==0) then
      write(vtu_unit,'(a)',iostat=io_status) &
           '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    endif
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</Points>'
    if (io_status==0) call mt_write_vtu_cells_empty(vtu_unit,io_status)
    if (io_status==0) call mt_write_vtu_file_footer(vtu_unit,io_status)
    if (io_status/=0) then
      close(vtu_unit)
      call mpistop(trim(caller)//' could not write VTU file')
    endif

    close(vtu_unit)
  end subroutine mt_write_fieldline_products_vtu_empty

  subroutine mt_write_fieldline_products_vtu_vertices(vtu_file, &
       length_results,twist_results,q_results,qperp_results,npoint, &
       do_length,do_twist,do_q,do_qperp,caller)
    integer, intent(in) :: npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: q_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    character(len=*), intent(in) :: vtu_file,caller
    logical, intent(in) :: do_length,do_twist,do_q,do_qperp

    integer :: vtu_unit,io_status

    open(newunit=vtu_unit,file=trim(vtu_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTU file')
    endif

    call mt_write_vtu_file_header(vtu_unit,npoint,npoint)
    call mt_write_vtu_product_pointdata(vtu_unit,length_results, &
         twist_results,q_results,qperp_results,npoint,do_twist,do_q, &
         do_qperp,io_status,do_length=do_length)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '<CellData>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</CellData>'
    if (io_status==0) call mt_write_vtu_points(vtu_unit,length_results, &
         npoint,io_status)
    if (io_status==0) call mt_write_vtu_vertex_cells(vtu_unit,npoint, &
         io_status)
    if (io_status==0) call mt_write_vtu_file_footer(vtu_unit,io_status)
    if (io_status/=0) then
      close(vtu_unit)
      call mpistop(trim(caller)//' could not write VTU file')
    endif

    close(vtu_unit)
  end subroutine mt_write_fieldline_products_vtu_vertices

  subroutine mt_write_fieldline_products_vtu_plane(vtu_file, &
       length_results,twist_results,q_results,qperp_results,n1,n2,do_twist, &
       do_q,do_qperp,caller,do_length)
    integer, intent(in) :: n1,n2
    type(trace_length_result), intent(in) :: length_results(n1*n2)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: q_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    character(len=*), intent(in) :: vtu_file,caller
    logical, intent(in) :: do_twist,do_q,do_qperp
    logical, intent(in), optional :: do_length

    integer :: npoint,ncell,vtu_unit,io_status
    logical :: do_length_eff

    npoint=n1*n2
    if (n1>1 .and. n2>1) then
      ncell=(n1-1)*(n2-1)
    else
      ncell=npoint
    endif

    open(newunit=vtu_unit,file=trim(vtu_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTU file')
    endif

    call mt_write_vtu_file_header(vtu_unit,npoint,ncell)
    do_length_eff=.true.
    if (present(do_length)) do_length_eff=do_length
    call mt_write_vtu_product_pointdata(vtu_unit,length_results, &
         twist_results,q_results,qperp_results,npoint,do_twist,do_q, &
         do_qperp,io_status,do_length=do_length_eff)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '<CellData>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</CellData>'
    if (io_status==0) call mt_write_vtu_points(vtu_unit,length_results, &
         npoint,io_status)
    if (io_status==0) then
      if (n1>1 .and. n2>1) then
        call mt_write_vtu_quad_cells(vtu_unit,n1,n2,io_status)
      else
        call mt_write_vtu_vertex_cells(vtu_unit,npoint,io_status)
      endif
    endif
    if (io_status==0) call mt_write_vtu_file_footer(vtu_unit,io_status)
    if (io_status/=0) then
      close(vtu_unit)
      call mpistop(trim(caller)//' could not write VTU file')
    endif

    close(vtu_unit)
  end subroutine mt_write_fieldline_products_vtu_plane

  subroutine mt_write_qsl_plane_vtu(vtu_file,length_results,twist_results, &
       mapping_results,qperp_results,n1,n2,caller)
    integer, intent(in) :: n1,n2
    type(trace_length_result), intent(in) :: length_results(n1*n2)
    type(trace_twist_result), intent(in) :: twist_results(n1*n2)
    type(trace_mapping_result), intent(in) :: mapping_results(n1*n2)
    type(trace_qperp_result), intent(in) :: qperp_results(n1*n2)
    character(len=*), intent(in) :: vtu_file,caller

    integer :: npoint,ncell,vtu_unit,io_status

    npoint=n1*n2
    if (n1>1 .and. n2>1) then
      ncell=(n1-1)*(n2-1)
    else
      ncell=npoint
    endif

    open(newunit=vtu_unit,file=trim(vtu_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTU file')
    endif

    call mt_write_vtu_file_header(vtu_unit,npoint,ncell)
    call mt_write_qsl_plane_vtu_pointdata(vtu_unit,length_results, &
         twist_results,mapping_results,qperp_results,npoint,io_status)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '<CellData>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</CellData>'
    if (io_status==0) call mt_write_vtu_points(vtu_unit,length_results, &
         npoint,io_status)
    if (io_status==0) then
      if (n1>1 .and. n2>1) then
        call mt_write_vtu_quad_cells(vtu_unit,n1,n2,io_status)
      else
        call mt_write_vtu_vertex_cells(vtu_unit,npoint,io_status)
      endif
    endif
    if (io_status==0) call mt_write_vtu_file_footer(vtu_unit,io_status)
    if (io_status/=0) then
      close(vtu_unit)
      call mpistop(trim(caller)//' could not write VTU file')
    endif

    close(vtu_unit)
  end subroutine mt_write_qsl_plane_vtu

  subroutine mt_write_qsl_plane_vtu_pointdata(vtu_unit,length_results, &
       twist_results,mapping_results,qperp_results,npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    type(trace_twist_result), intent(in) :: twist_results(npoint)
    type(trace_mapping_result), intent(in) :: mapping_results(npoint)
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    write(vtu_unit,'(a)',iostat=io_status) '<PointData>'
    if (io_status/=0) return

    if (mt_vtk_detail_is_full()) then
      call mt_write_vtu_length_pointdata(vtu_unit,length_results,npoint, &
           io_status)
      call mt_write_vtu_twist_pointdata(vtu_unit,twist_results,npoint, &
           io_status)
      call mt_write_vtu_mapping_pointdata(vtu_unit,mapping_results,npoint, &
           io_status)
      call mt_write_vtu_q_product_pointdata(vtu_unit,qperp_results, &
           npoint,.true.,io_status)
      call mt_write_vtu_qperp_public_pointdata(vtu_unit,qperp_results, &
           npoint,io_status)
      call mt_write_vtu_qperp_method2_pointdata(vtu_unit,qperp_results, &
           npoint,io_status)
    else
      call mt_write_vtu_qsl_minimal_pointdata(vtu_unit,length_results, &
           twist_results,qperp_results,npoint,io_status)
    endif

    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</PointData>'
  end subroutine mt_write_qsl_plane_vtu_pointdata

  subroutine mt_write_vtu_qsl_minimal_pointdata(vtu_unit,length_results, &
       twist_results,qperp_results,npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    type(trace_twist_result), intent(in) :: twist_results(npoint)
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'length_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_float(length_results(ipoint)%total_length), &
         ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'twist_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_float(twist_results(ipoint)%total_twist), &
         ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'logQ',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_valid_float(qperp_results(ipoint)%logq0, &
         qperp_results(ipoint)%valid_q0),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'logQperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_valid_float(qperp_results(ipoint)%logqperp, &
         qperp_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

  end subroutine mt_write_vtu_qsl_minimal_pointdata

  double precision function mt_vti_axis_spacing(xmin,xmax,n)
    integer, intent(in) :: n
    double precision, intent(in) :: xmin,xmax

    if (n>1) then
      mt_vti_axis_spacing=(xmax-xmin)/dble(n-1)
    else
      mt_vti_axis_spacing=0.d0
    endif
  end function mt_vti_axis_spacing

  subroutine mt_allocate_volume_products(products,nseed,do_twist,do_q, &
       do_qperp)
    type(mt_volume_products), intent(inout) :: products
    integer, intent(in) :: nseed
    logical, intent(in) :: do_twist,do_q,do_qperp

    allocate(products%length_total(nseed))
    allocate(products%length_backward(nseed))
    allocate(products%length_forward(nseed))
    allocate(products%nstep_backward_length(nseed))
    allocate(products%nstep_forward_length(nseed))
    allocate(products%status_backward_length(nseed))
    allocate(products%status_forward_length(nseed))

    if (do_twist) then
      allocate(products%twist_total(nseed))
      allocate(products%twist_backward(nseed))
      allocate(products%twist_forward(nseed))
      allocate(products%nstep_backward_twist(nseed))
      allocate(products%nstep_forward_twist(nseed))
      allocate(products%status_backward_twist(nseed))
      allocate(products%status_forward_twist(nseed))
    endif

    if (do_q) then
      allocate(products%q(nseed))
      allocate(products%logq(nseed))
      allocate(products%N2_q(nseed))
      allocate(products%bfactor_q(nseed))
      allocate(products%length_forward_q(nseed))
      allocate(products%length_backward_q(nseed))
      allocate(products%Bseed_norm_q(nseed))
      allocate(products%Bf_norm_q(nseed))
      allocate(products%Bb_norm_q(nseed))
      allocate(products%valid_q(nseed))
      allocate(products%status_q(nseed))
      allocate(products%face_forward_q(nseed))
      allocate(products%face_backward_q(nseed))
      allocate(products%status_forward_q(nseed))
      allocate(products%status_backward_q(nseed))
    endif

    if (do_qperp) then
      allocate(products%qperp(nseed))
      allocate(products%logqperp(nseed))
      allocate(products%N2(nseed))
      allocate(products%bfactor(nseed))
      allocate(products%length_forward_qperp(nseed))
      allocate(products%length_backward_qperp(nseed))
      allocate(products%Bseed_norm(nseed))
      allocate(products%Bf_norm(nseed))
      allocate(products%Bb_norm(nseed))
      allocate(products%valid_qperp(nseed))
      allocate(products%status_qperp(nseed))
      allocate(products%face_forward_qperp(nseed))
      allocate(products%face_backward_qperp(nseed))
      allocate(products%status_forward_qperp(nseed))
      allocate(products%status_backward_qperp(nseed))
    endif
  end subroutine mt_allocate_volume_products

  subroutine mt_deallocate_volume_products(products)
    type(mt_volume_products), intent(inout) :: products

    if (allocated(products%length_total)) deallocate(products%length_total)
    if (allocated(products%length_backward)) deallocate(products%length_backward)
    if (allocated(products%length_forward)) deallocate(products%length_forward)
    if (allocated(products%nstep_backward_length)) &
         deallocate(products%nstep_backward_length)
    if (allocated(products%nstep_forward_length)) &
         deallocate(products%nstep_forward_length)
    if (allocated(products%status_backward_length)) &
         deallocate(products%status_backward_length)
    if (allocated(products%status_forward_length)) &
         deallocate(products%status_forward_length)

    if (allocated(products%twist_total)) deallocate(products%twist_total)
    if (allocated(products%twist_backward)) deallocate(products%twist_backward)
    if (allocated(products%twist_forward)) deallocate(products%twist_forward)
    if (allocated(products%nstep_backward_twist)) &
         deallocate(products%nstep_backward_twist)
    if (allocated(products%nstep_forward_twist)) &
         deallocate(products%nstep_forward_twist)
    if (allocated(products%status_backward_twist)) &
         deallocate(products%status_backward_twist)
    if (allocated(products%status_forward_twist)) &
         deallocate(products%status_forward_twist)

    if (allocated(products%q)) deallocate(products%q)
    if (allocated(products%logq)) deallocate(products%logq)
    if (allocated(products%N2_q)) deallocate(products%N2_q)
    if (allocated(products%bfactor_q)) deallocate(products%bfactor_q)
    if (allocated(products%length_forward_q)) &
         deallocate(products%length_forward_q)
    if (allocated(products%length_backward_q)) &
         deallocate(products%length_backward_q)
    if (allocated(products%Bseed_norm_q)) deallocate(products%Bseed_norm_q)
    if (allocated(products%Bf_norm_q)) deallocate(products%Bf_norm_q)
    if (allocated(products%Bb_norm_q)) deallocate(products%Bb_norm_q)
    if (allocated(products%valid_q)) deallocate(products%valid_q)
    if (allocated(products%status_q)) deallocate(products%status_q)
    if (allocated(products%face_forward_q)) &
         deallocate(products%face_forward_q)
    if (allocated(products%face_backward_q)) &
         deallocate(products%face_backward_q)
    if (allocated(products%status_forward_q)) &
         deallocate(products%status_forward_q)
    if (allocated(products%status_backward_q)) &
         deallocate(products%status_backward_q)

    if (allocated(products%qperp)) deallocate(products%qperp)
    if (allocated(products%logqperp)) deallocate(products%logqperp)
    if (allocated(products%N2)) deallocate(products%N2)
    if (allocated(products%bfactor)) deallocate(products%bfactor)
    if (allocated(products%length_forward_qperp)) &
         deallocate(products%length_forward_qperp)
    if (allocated(products%length_backward_qperp)) &
         deallocate(products%length_backward_qperp)
    if (allocated(products%Bseed_norm)) deallocate(products%Bseed_norm)
    if (allocated(products%Bf_norm)) deallocate(products%Bf_norm)
    if (allocated(products%Bb_norm)) deallocate(products%Bb_norm)
    if (allocated(products%valid_qperp)) deallocate(products%valid_qperp)
    if (allocated(products%status_qperp)) deallocate(products%status_qperp)
    if (allocated(products%face_forward_qperp)) &
         deallocate(products%face_forward_qperp)
    if (allocated(products%face_backward_qperp)) &
         deallocate(products%face_backward_qperp)
    if (allocated(products%status_forward_qperp)) &
         deallocate(products%status_forward_qperp)
    if (allocated(products%status_backward_qperp)) &
         deallocate(products%status_backward_qperp)
  end subroutine mt_deallocate_volume_products

  subroutine mt_build_volume_slab_seeds(seeds,xmin,ymin,zmin,spacing, &
       nx,ny,k_start,slab_nz)
    double precision, intent(out) :: seeds(:,:)
    double precision, intent(in) :: xmin,ymin,zmin,spacing(3)
    integer, intent(in) :: nx,ny,k_start,slab_nz

    integer :: i,j,k,kk,iseed

    do kk=1,slab_nz
      k=k_start+kk-1
      do j=1,ny
        do i=1,nx
          iseed=(kk-1)*nx*ny+(j-1)*nx+i
          seeds(iseed,1)=xmin+dble(i-1)*spacing(1)
          seeds(iseed,2)=ymin+dble(j-1)*spacing(2)
          seeds(iseed,3)=zmin+dble(k-1)*spacing(3)
        enddo
      enddo
    enddo
  end subroutine mt_build_volume_slab_seeds

  subroutine mt_copy_volume_length_slab(products,results,nx,ny,k_start, &
       slab_nz)
    type(mt_volume_products), intent(inout) :: products
    type(trace_length_result), intent(in) :: results(:)
    integer, intent(in) :: nx,ny,k_start,slab_nz

    integer :: i,j,kk,k,ilocal,iglobal

    do kk=1,slab_nz
      k=k_start+kk-1
      do j=1,ny
        do i=1,nx
          ilocal=(kk-1)*nx*ny+(j-1)*nx+i
          iglobal=(k-1)*nx*ny+(j-1)*nx+i
          products%length_total(iglobal)=results(ilocal)%total_length
          products%length_backward(iglobal)=results(ilocal)%backward_length
          products%length_forward(iglobal)=results(ilocal)%forward_length
          products%nstep_backward_length(iglobal)= &
               results(ilocal)%backward_nstep
          products%nstep_forward_length(iglobal)= &
               results(ilocal)%forward_nstep
          products%status_backward_length(iglobal)= &
               results(ilocal)%backward_status
          products%status_forward_length(iglobal)= &
               results(ilocal)%forward_status
        enddo
      enddo
    enddo
  end subroutine mt_copy_volume_length_slab

  subroutine mt_copy_volume_twist_slab(products,results,nx,ny,k_start, &
       slab_nz)
    type(mt_volume_products), intent(inout) :: products
    type(trace_twist_result), intent(in) :: results(:)
    integer, intent(in) :: nx,ny,k_start,slab_nz

    integer :: i,j,kk,k,ilocal,iglobal

    do kk=1,slab_nz
      k=k_start+kk-1
      do j=1,ny
        do i=1,nx
          ilocal=(kk-1)*nx*ny+(j-1)*nx+i
          iglobal=(k-1)*nx*ny+(j-1)*nx+i
          products%twist_total(iglobal)=results(ilocal)%total_twist
          products%twist_backward(iglobal)=results(ilocal)%backward_twist
          products%twist_forward(iglobal)=results(ilocal)%forward_twist
          products%nstep_backward_twist(iglobal)= &
               results(ilocal)%line%backward_nstep
          products%nstep_forward_twist(iglobal)= &
               results(ilocal)%line%forward_nstep
          products%status_backward_twist(iglobal)= &
               results(ilocal)%line%backward_status
          products%status_forward_twist(iglobal)= &
               results(ilocal)%line%forward_status
        enddo
      enddo
    enddo
  end subroutine mt_copy_volume_twist_slab

  subroutine mt_copy_volume_q_slab(products,results,nx,ny,k_start,slab_nz)
    type(mt_volume_products), intent(inout) :: products
    type(trace_qperp_result), intent(in) :: results(:)
    integer, intent(in) :: nx,ny,k_start,slab_nz

    integer :: i,j,kk,k,ilocal,iglobal

    do kk=1,slab_nz
      k=k_start+kk-1
      do j=1,ny
        do i=1,nx
          ilocal=(kk-1)*nx*ny+(j-1)*nx+i
          iglobal=(k-1)*nx*ny+(j-1)*nx+i
          products%q(iglobal)=results(ilocal)%q0
          products%logq(iglobal)=results(ilocal)%logq0
          products%N2_q(iglobal)=results(ilocal)%N2_qperp0
          products%bfactor_q(iglobal)=results(ilocal)%bfactor_qperp0
          products%length_forward_q(iglobal)=results(ilocal)%forward_length
          products%length_backward_q(iglobal)=results(ilocal)%backward_length
          products%Bseed_norm_q(iglobal)=dsqrt(sum(results(ilocal)%B_seed**2))
          products%Bf_norm_q(iglobal)=dsqrt(sum(results(ilocal)%forward_B**2))
          products%Bb_norm_q(iglobal)=dsqrt(sum(results(ilocal)%backward_B**2))
          products%valid_q(iglobal)=merge(1,0,results(ilocal)%valid_q0)
          products%status_q(iglobal)=results(ilocal)%status_q0
          products%face_forward_q(iglobal)=results(ilocal)%forward_face
          products%face_backward_q(iglobal)=results(ilocal)%backward_face
          products%status_forward_q(iglobal)=results(ilocal)%forward_status
          products%status_backward_q(iglobal)=results(ilocal)%backward_status
        enddo
      enddo
    enddo
  end subroutine mt_copy_volume_q_slab

  subroutine mt_copy_volume_qperp_slab(products,results,nx,ny,k_start, &
       slab_nz)
    type(mt_volume_products), intent(inout) :: products
    type(trace_qperp_result), intent(in) :: results(:)
    integer, intent(in) :: nx,ny,k_start,slab_nz

    integer :: i,j,kk,k,ilocal,iglobal

    do kk=1,slab_nz
      k=k_start+kk-1
      do j=1,ny
        do i=1,nx
          ilocal=(kk-1)*nx*ny+(j-1)*nx+i
          iglobal=(k-1)*nx*ny+(j-1)*nx+i
          products%qperp(iglobal)=results(ilocal)%qperp
          products%logqperp(iglobal)=results(ilocal)%logqperp
          products%N2(iglobal)=results(ilocal)%N2
          products%bfactor(iglobal)=results(ilocal)%bfactor
          products%length_forward_qperp(iglobal)= &
               results(ilocal)%forward_length
          products%length_backward_qperp(iglobal)= &
               results(ilocal)%backward_length
          products%Bseed_norm(iglobal)=dsqrt(sum(results(ilocal)%B_seed**2))
          products%Bf_norm(iglobal)=dsqrt(sum(results(ilocal)%forward_B**2))
          products%Bb_norm(iglobal)=dsqrt(sum(results(ilocal)%backward_B**2))
          products%valid_qperp(iglobal)=merge(1,0,results(ilocal)%valid)
          products%status_qperp(iglobal)=results(ilocal)%status
          products%face_forward_qperp(iglobal)=results(ilocal)%forward_face
          products%face_backward_qperp(iglobal)=results(ilocal)%backward_face
          products%status_forward_qperp(iglobal)= &
               results(ilocal)%forward_status
          products%status_backward_qperp(iglobal)= &
               results(ilocal)%backward_status
        enddo
      enddo
    enddo
  end subroutine mt_copy_volume_qperp_slab

  subroutine mt_build_vti_desc(desc,name,array_kind,npoint)
    type(mt_vti_array_desc), intent(out) :: desc
    character(len=*), intent(in) :: name
    integer, intent(in) :: array_kind,npoint

    desc%name=''
    desc%name=trim(name)
    desc%kind=array_kind
    select case (array_kind)
    case (mt_vti_kind_float64)
      desc%nbytes=int(npoint,kind=8)*8_8
    case (mt_vti_kind_int32)
      desc%nbytes=int(npoint,kind=8)*4_8
    case default
      desc%nbytes=0_8
    end select
    desc%offset=0_8
  end subroutine mt_build_vti_desc

  subroutine mt_finalize_vti_desc_offsets(descs,ndesc)
    integer, intent(in) :: ndesc
    type(mt_vti_array_desc), intent(inout) :: descs(ndesc)

    integer :: idesc
    integer(kind=8) :: offset

    offset=0_8
    do idesc=1,ndesc
      descs(idesc)%offset=offset
      offset=offset+descs(idesc)%nbytes+4_8
    enddo
  end subroutine mt_finalize_vti_desc_offsets

  subroutine mt_append_vti_desc(descs,idesc,name,array_kind,npoint)
    type(mt_vti_array_desc), intent(inout) :: descs(:)
    integer, intent(inout) :: idesc
    character(len=*), intent(in) :: name
    integer, intent(in) :: array_kind,npoint

    idesc=idesc+1
    call mt_build_vti_desc(descs(idesc),name,array_kind,npoint)
  end subroutine mt_append_vti_desc

  subroutine mt_build_cartesian_vti_pointdata_descs(npoint,descs,ndesc)
    integer, intent(in) :: npoint
    type(mt_vti_array_desc), allocatable, intent(out) :: descs(:)
    integer, intent(out) :: ndesc

    ndesc=3
    allocate(descs(ndesc))
    call mt_build_vti_desc(descs(1),'length_total',mt_vti_kind_float64, &
         npoint)
    call mt_build_vti_desc(descs(2),'qperp',mt_vti_kind_float64,npoint)
    call mt_build_vti_desc(descs(3),'status',mt_vti_kind_int32,npoint)
    call mt_finalize_vti_desc_offsets(descs,ndesc)
  end subroutine mt_build_cartesian_vti_pointdata_descs

  subroutine mt_build_volume_vti_descs(npoint,do_length,do_twist,do_q, &
       do_qperp,descs,ndesc)
    integer, intent(in) :: npoint
    logical, intent(in) :: do_length,do_twist,do_q,do_qperp
    type(mt_vti_array_desc), allocatable, intent(out) :: descs(:)
    integer, intent(out) :: ndesc

    integer :: idesc,max_desc

    if (.not.mt_vtk_detail_is_full()) then
      max_desc=0
      if (do_length) max_desc=max_desc+1
      if (do_twist) max_desc=max_desc+1
      if (do_q) max_desc=max_desc+1
      if (do_qperp) max_desc=max_desc+1
    else
      max_desc=0
      if (do_length) max_desc=max_desc+7
      if (do_twist) max_desc=max_desc+7
      if (do_q) max_desc=max_desc+15
      if (do_qperp) max_desc=max_desc+15
    endif
    allocate(descs(max_desc))

    idesc=0
    if (do_length) then
      call mt_append_vti_desc(descs,idesc,'length_total', &
           mt_vti_kind_float64,npoint)
    endif
    if (.not.mt_vtk_detail_is_full()) then
      if (do_twist) then
        call mt_append_vti_desc(descs,idesc,'twist_total', &
             mt_vti_kind_float64,npoint)
      endif
      if (do_q) then
        call mt_append_vti_desc(descs,idesc,'logQ', &
             mt_vti_kind_float64,npoint)
      endif
      if (do_qperp) then
        call mt_append_vti_desc(descs,idesc,'logQperp', &
             mt_vti_kind_float64,npoint)
      endif
      ndesc=idesc
      call mt_finalize_vti_desc_offsets(descs,ndesc)
      return
    endif

    if (do_length) then
      call mt_append_vti_desc(descs,idesc,'length_backward', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'length_forward', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'nstep_backward_length', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'nstep_forward_length', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_backward_length', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_forward_length', &
           mt_vti_kind_int32,npoint)
    endif

    if (do_twist) then
      call mt_append_vti_desc(descs,idesc,'twist_total', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'twist_backward', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'twist_forward', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'nstep_backward_twist', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'nstep_forward_twist', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_backward_twist', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_forward_twist', &
           mt_vti_kind_int32,npoint)
    endif

    if (do_q) then
      call mt_append_vti_desc(descs,idesc,'q',mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'logq', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'N2_q',mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'bfactor_q', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'length_forward_q', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'length_backward_q', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'Bseed_norm_q', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'Bf_norm_q', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'Bb_norm_q', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'valid_Q', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_Q', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'face_forward_Q', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'face_backward_Q', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_forward_Q', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_backward_Q', &
           mt_vti_kind_int32,npoint)
    endif

    if (do_qperp) then
      call mt_append_vti_desc(descs,idesc,'qperp', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'logqperp', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'N2',mt_vti_kind_float64, &
           npoint)
      call mt_append_vti_desc(descs,idesc,'bfactor', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'length_forward_qperp', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'length_backward_qperp', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'Bseed_norm', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'Bf_norm', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'Bb_norm', &
           mt_vti_kind_float64,npoint)
      call mt_append_vti_desc(descs,idesc,'valid_qperp', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_qperp', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'face_forward_qperp', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'face_backward_qperp', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_forward_qperp', &
           mt_vti_kind_int32,npoint)
      call mt_append_vti_desc(descs,idesc,'status_backward_qperp', &
           mt_vti_kind_int32,npoint)
    endif

    ndesc=idesc
    call mt_finalize_vti_desc_offsets(descs,ndesc)
  end subroutine mt_build_volume_vti_descs

  character(len=8) function mt_vti_type_name(array_kind)
    integer, intent(in) :: array_kind

    select case (array_kind)
    case (mt_vti_kind_float64)
      mt_vti_type_name='Float64'
    case (mt_vti_kind_int32)
      mt_vti_type_name='Int32'
    case default
      mt_vti_type_name='Unknown'
    end select
  end function mt_vti_type_name

  subroutine mt_write_vti_image_header(vti_unit,origin,spacing,nx,ny,nz, &
       descs,ndesc,io_status)
    integer, intent(in) :: vti_unit,nx,ny,nz,ndesc
    double precision, intent(in) :: origin(3),spacing(3)
    type(mt_vti_array_desc), intent(in) :: descs(ndesc)
    integer, intent(out) :: io_status

    integer :: extent(6),idesc

    extent=(/ 0,nx-1,0,ny-1,0,nz-1 /)
    write(vti_unit,'(a)',iostat=io_status) '<?xml version="1.0"?>'
    if (io_status/=0) return
    write(vti_unit,'(a)',iostat=io_status) &
         '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'
    if (io_status/=0) return
    write(vti_unit,'(a,3(1pe24.16),a,6(i0,1x),a,3(1pe24.16),a)', &
         iostat=io_status) '  <ImageData Origin="',origin, &
         '" WholeExtent="',extent,'" Spacing="',spacing,'">'
    if (io_status/=0) return
    write(vti_unit,'(a,6(i0,1x),a)',iostat=io_status) &
         '    <Piece Extent="',extent,'">'
    if (io_status/=0) return
    write(vti_unit,'(a)',iostat=io_status) '      <PointData>'
    if (io_status/=0) return

    do idesc=1,ndesc
      call mt_write_vti_appended_array(vti_unit, &
           mt_vti_type_name(descs(idesc)%kind),trim(descs(idesc)%name), &
           descs(idesc)%offset,io_status)
      if (io_status/=0) return
    enddo

    write(vti_unit,'(a)',iostat=io_status) '      </PointData>'
    if (io_status/=0) return
    write(vti_unit,'(a)',iostat=io_status) '    </Piece>'
    if (io_status/=0) return
    write(vti_unit,'(a)',iostat=io_status) '  </ImageData>'
    if (io_status/=0) return
    write(vti_unit,'(a)',iostat=io_status) '<AppendedData encoding="raw">'
  end subroutine mt_write_vti_image_header

  subroutine mt_write_fieldline_products_volume_vti(vti_file,origin, &
       spacing,nx,ny,nz,products,do_length,do_twist,do_q,do_qperp,caller)
    integer, intent(in) :: nx,ny,nz
    double precision, intent(in) :: origin(3),spacing(3)
    type(mt_volume_products), intent(in) :: products
    logical, intent(in) :: do_length,do_twist,do_q,do_qperp
    character(len=*), intent(in) :: vti_file,caller

    integer :: vti_unit,io_status,npoint
    character(len=1) :: marker
    type(mt_vti_array_desc), allocatable :: descs(:)
    integer :: ndesc

    npoint=nx*ny*nz
    call mt_check_vti_byte_count(npoint,caller)
    call mt_build_volume_vti_descs(npoint,do_length,do_twist,do_q, &
         do_qperp,descs,ndesc)

    open(newunit=vti_unit,file=trim(vti_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTI file')
    endif

    call mt_write_vti_image_header(vti_unit,origin,spacing,nx,ny,nz, &
         descs,ndesc,io_status)
    if (io_status/=0) then
      close(vti_unit)
      call mpistop(trim(caller)//' could not write VTI header')
    endif
    close(vti_unit)

    open(newunit=vti_unit,file=trim(vti_file),access='stream', &
         form='unformatted',status='old',position='append', &
         action='write',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not append VTI payload')
    endif

    marker='_'
    write(vti_unit,iostat=io_status) marker
    if (io_status==0) then
      call mt_write_volume_vti_payload(vti_unit,products,npoint,descs, &
           ndesc,io_status)
    endif
    close(vti_unit)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not write VTI payload')
    endif

    open(newunit=vti_unit,file=trim(vti_file),status='old', &
         action='write',form='formatted',position='append',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not append VTI footer')
    endif
    write(vti_unit,'(a)',iostat=io_status) '</AppendedData>'
    if (io_status==0) write(vti_unit,'(a)',iostat=io_status) '</VTKFile>'
    close(vti_unit)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not write VTI footer')
    endif
    deallocate(descs)
  end subroutine mt_write_fieldline_products_volume_vti

  subroutine mt_write_volume_vti_payload(vti_unit,products,npoint,descs, &
       ndesc,io_status)
    integer, intent(in) :: vti_unit,npoint,ndesc
    type(mt_volume_products), intent(in) :: products
    type(mt_vti_array_desc), intent(in) :: descs(ndesc)
    integer, intent(inout) :: io_status

    integer :: idesc

    do idesc=1,ndesc
      select case (trim(descs(idesc)%name))
      case ('length_total')
        if (mt_vtk_detail_is_full()) then
          call mt_write_vti_payload_float64(vti_unit,products%length_total, &
               npoint,io_status)
        else
          call mt_write_vti_payload_float64_visual(vti_unit, &
               products%length_total,npoint,io_status)
        endif
      case ('length_backward')
        call mt_write_vti_payload_float64(vti_unit,products%length_backward, &
             npoint,io_status)
      case ('length_forward')
        call mt_write_vti_payload_float64(vti_unit,products%length_forward, &
             npoint,io_status)
      case ('nstep_backward_length')
        call mt_write_vti_payload_int32(vti_unit, &
             products%nstep_backward_length,npoint,io_status)
      case ('nstep_forward_length')
        call mt_write_vti_payload_int32(vti_unit, &
             products%nstep_forward_length,npoint,io_status)
      case ('status_backward_length')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_backward_length,npoint,io_status)
      case ('status_forward_length')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_forward_length,npoint,io_status)
      case ('twist_total')
        if (mt_vtk_detail_is_full()) then
          call mt_write_vti_payload_float64(vti_unit,products%twist_total, &
               npoint,io_status)
        else
          call mt_write_vti_payload_float64_visual(vti_unit, &
               products%twist_total,npoint,io_status)
        endif
      case ('twist_backward')
        call mt_write_vti_payload_float64(vti_unit,products%twist_backward, &
             npoint,io_status)
      case ('twist_forward')
        call mt_write_vti_payload_float64(vti_unit,products%twist_forward, &
             npoint,io_status)
      case ('nstep_backward_twist')
        call mt_write_vti_payload_int32(vti_unit, &
             products%nstep_backward_twist,npoint,io_status)
      case ('nstep_forward_twist')
        call mt_write_vti_payload_int32(vti_unit, &
             products%nstep_forward_twist,npoint,io_status)
      case ('status_backward_twist')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_backward_twist,npoint,io_status)
      case ('status_forward_twist')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_forward_twist,npoint,io_status)
      case ('q')
        call mt_write_vti_payload_float64(vti_unit,products%q,npoint, &
             io_status)
      case ('logq')
        call mt_write_vti_payload_float64(vti_unit,products%logq,npoint, &
             io_status)
      case ('logQ')
        call mt_write_vti_payload_logq_visual(vti_unit,products,npoint, &
             io_status)
      case ('N2_q')
        call mt_write_vti_payload_float64(vti_unit,products%N2_q,npoint, &
             io_status)
      case ('bfactor_q')
        call mt_write_vti_payload_float64(vti_unit,products%bfactor_q, &
             npoint,io_status)
      case ('length_forward_q')
        call mt_write_vti_payload_float64(vti_unit, &
             products%length_forward_q,npoint,io_status)
      case ('length_backward_q')
        call mt_write_vti_payload_float64(vti_unit, &
             products%length_backward_q,npoint,io_status)
      case ('Bseed_norm_q')
        call mt_write_vti_payload_float64(vti_unit,products%Bseed_norm_q, &
             npoint,io_status)
      case ('Bf_norm_q')
        call mt_write_vti_payload_float64(vti_unit,products%Bf_norm_q,npoint, &
             io_status)
      case ('Bb_norm_q')
        call mt_write_vti_payload_float64(vti_unit,products%Bb_norm_q,npoint, &
             io_status)
      case ('valid_Q')
        call mt_write_vti_payload_int32(vti_unit,products%valid_q,npoint, &
             io_status)
      case ('status_Q')
        call mt_write_vti_payload_int32(vti_unit,products%status_q,npoint, &
             io_status)
      case ('face_forward_Q')
        call mt_write_vti_payload_int32(vti_unit,products%face_forward_q, &
             npoint,io_status)
      case ('face_backward_Q')
        call mt_write_vti_payload_int32(vti_unit,products%face_backward_q, &
             npoint,io_status)
      case ('status_forward_Q')
        call mt_write_vti_payload_int32(vti_unit,products%status_forward_q, &
             npoint,io_status)
      case ('status_backward_Q')
        call mt_write_vti_payload_int32(vti_unit,products%status_backward_q, &
             npoint,io_status)
      case ('qperp')
        call mt_write_vti_payload_float64(vti_unit,products%qperp,npoint, &
             io_status)
      case ('logqperp')
        call mt_write_vti_payload_float64(vti_unit,products%logqperp, &
             npoint,io_status)
      case ('logQperp')
        call mt_write_vti_payload_logqperp_visual(vti_unit,products,npoint, &
             io_status)
      case ('N2')
        call mt_write_vti_payload_float64(vti_unit,products%N2,npoint, &
             io_status)
      case ('bfactor')
        call mt_write_vti_payload_float64(vti_unit,products%bfactor,npoint, &
             io_status)
      case ('length_forward_qperp')
        call mt_write_vti_payload_float64(vti_unit, &
             products%length_forward_qperp,npoint,io_status)
      case ('length_backward_qperp')
        call mt_write_vti_payload_float64(vti_unit, &
             products%length_backward_qperp,npoint,io_status)
      case ('Bseed_norm')
        call mt_write_vti_payload_float64(vti_unit,products%Bseed_norm, &
             npoint,io_status)
      case ('Bf_norm')
        call mt_write_vti_payload_float64(vti_unit,products%Bf_norm,npoint, &
             io_status)
      case ('Bb_norm')
        call mt_write_vti_payload_float64(vti_unit,products%Bb_norm,npoint, &
             io_status)
      case ('valid_qperp')
        call mt_write_vti_payload_int32(vti_unit,products%valid_qperp, &
             npoint,io_status)
      case ('valid_Qperp')
        call mt_write_vti_payload_int32(vti_unit,products%valid_qperp, &
             npoint,io_status)
      case ('status_qperp')
        call mt_write_vti_payload_int32(vti_unit,products%status_qperp, &
             npoint,io_status)
      case ('status_Qperp')
        call mt_write_vti_payload_int32(vti_unit,products%status_qperp, &
             npoint,io_status)
      case ('face_forward_qperp')
        call mt_write_vti_payload_int32(vti_unit, &
             products%face_forward_qperp,npoint,io_status)
      case ('face_backward_qperp')
        call mt_write_vti_payload_int32(vti_unit, &
             products%face_backward_qperp,npoint,io_status)
      case ('status_forward_qperp')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_forward_qperp,npoint,io_status)
      case ('status_backward_qperp')
        call mt_write_vti_payload_int32(vti_unit, &
             products%status_backward_qperp,npoint,io_status)
      case default
        io_status=1
      end select
      if (io_status/=0) exit
    enddo
  end subroutine mt_write_volume_vti_payload

  subroutine mt_write_vti_payload_logq_visual(vti_unit,products,npoint, &
       io_status)
    integer, intent(in) :: vti_unit,npoint
    type(mt_volume_products), intent(in) :: products
    integer, intent(inout) :: io_status

    double precision, allocatable :: visual_values(:)
    integer :: ipoint

    allocate(visual_values(npoint))
    do ipoint=1,npoint
      visual_values(ipoint)=mt_visual_valid_float(products%logq(ipoint), &
           products%valid_q(ipoint)==1)
    enddo
    call mt_write_vti_payload_float64(vti_unit,visual_values,npoint, &
         io_status)
    deallocate(visual_values)
  end subroutine mt_write_vti_payload_logq_visual

  subroutine mt_write_vti_payload_float64_visual(vti_unit,values,npoint, &
       io_status)
    integer, intent(in) :: vti_unit,npoint
    double precision, intent(in) :: values(npoint)
    integer, intent(inout) :: io_status

    double precision, allocatable :: visual_values(:)
    integer :: ipoint

    if (io_status/=0) return
    allocate(visual_values(npoint))
    do ipoint=1,npoint
      visual_values(ipoint)=mt_visual_float(values(ipoint))
    enddo
    call mt_write_vti_payload_float64(vti_unit,visual_values,npoint, &
         io_status)
    deallocate(visual_values)
  end subroutine mt_write_vti_payload_float64_visual

  subroutine mt_write_vti_payload_logqperp_visual(vti_unit,products,npoint, &
       io_status)
    integer, intent(in) :: vti_unit,npoint
    type(mt_volume_products), intent(in) :: products
    integer, intent(inout) :: io_status

    double precision, allocatable :: visual_values(:)
    integer :: ipoint

    if (io_status/=0) return
    allocate(visual_values(npoint))
    do ipoint=1,npoint
      visual_values(ipoint)=mt_visual_valid_float(products%logqperp(ipoint), &
           products%valid_qperp(ipoint)==1)
    enddo
    call mt_write_vti_payload_float64(vti_unit,visual_values,npoint, &
         io_status)
    deallocate(visual_values)
  end subroutine mt_write_vti_payload_logqperp_visual

  subroutine mt_check_vti_byte_count(npoint,caller)
    integer, intent(in) :: npoint
    character(len=*), intent(in) :: caller

    if (int(npoint,kind=8)*8_8>int(huge(0),kind=8)) then
      call mpistop(trim(caller)//' exceeds 32-bit VTI byte count')
    endif
  end subroutine mt_check_vti_byte_count

  subroutine mt_write_vti_payload_float64(vti_unit,values,npoint,io_status)
    integer, intent(in) :: vti_unit,npoint
    double precision, intent(in) :: values(npoint)
    integer, intent(inout) :: io_status

    integer :: byte_count

    if (io_status/=0) return
    byte_count=npoint*8
    write(vti_unit,iostat=io_status) byte_count
    if (io_status==0) write(vti_unit,iostat=io_status) values
  end subroutine mt_write_vti_payload_float64

  subroutine mt_write_vti_payload_int32(vti_unit,values,npoint,io_status)
    integer, intent(in) :: vti_unit,npoint
    integer, intent(in) :: values(npoint)
    integer, intent(inout) :: io_status

    integer :: byte_count

    if (io_status/=0) return
    byte_count=npoint*4
    write(vti_unit,iostat=io_status) byte_count
    if (io_status==0) write(vti_unit,iostat=io_status) values
  end subroutine mt_write_vti_payload_int32

  subroutine mt_write_vti_pointdata_fixed(vti_file,origin,spacing,nx,ny,nz, &
       length_total,qperp,status,caller)
    character(len=*), intent(in) :: vti_file,caller
    integer, intent(in) :: nx,ny,nz
    double precision, intent(in) :: origin(3),spacing(3)
    double precision, intent(in) :: length_total(nx*ny*nz),qperp(nx*ny*nz)
    integer, intent(in) :: status(nx*ny*nz)

    integer :: vti_unit,io_status,npoint
    integer :: idesc,ndesc
    character(len=1) :: marker
    type(mt_vti_array_desc), allocatable :: descs(:)

    npoint=nx*ny*nz
    call mt_check_vti_byte_count(npoint,caller)
    call mt_build_cartesian_vti_pointdata_descs(npoint,descs,ndesc)

    open(newunit=vti_unit,file=trim(vti_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open VTI file')
    endif

    call mt_write_vti_image_header(vti_unit,origin,spacing,nx,ny,nz, &
         descs,ndesc,io_status)
    if (io_status/=0) then
      close(vti_unit)
      call mpistop(trim(caller)//' could not write VTI header')
    endif
    close(vti_unit)

    open(newunit=vti_unit,file=trim(vti_file),access='stream', &
         form='unformatted',status='old',position='append', &
         action='write',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not append VTI payload')
    endif

    marker='_'
    write(vti_unit,iostat=io_status) marker
    do idesc=1,ndesc
      if (io_status/=0) exit
      select case (trim(descs(idesc)%name))
      case ('length_total')
        call mt_write_vti_payload_float64(vti_unit,length_total,npoint, &
             io_status)
      case ('qperp')
        call mt_write_vti_payload_float64(vti_unit,qperp,npoint,io_status)
      case ('status')
        call mt_write_vti_payload_int32(vti_unit,status,npoint,io_status)
      case default
        io_status=1
      end select
    enddo
    close(vti_unit)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not write VTI payload')
    endif

    open(newunit=vti_unit,file=trim(vti_file),status='old', &
         action='write',form='formatted',position='append',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not append VTI footer')
    endif
    write(vti_unit,'(a)',iostat=io_status) '</AppendedData>'
    if (io_status==0) write(vti_unit,'(a)',iostat=io_status) '</VTKFile>'
    close(vti_unit)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not write VTI footer')
    endif
    deallocate(descs)
  end subroutine mt_write_vti_pointdata_fixed

  subroutine mt_write_vti_appended_array(vti_unit,vtk_type,name,offset, &
       io_status)
    integer, intent(in) :: vti_unit
    character(len=*), intent(in) :: vtk_type,name
    integer(kind=8), intent(in) :: offset
    integer, intent(out) :: io_status

    write(vti_unit,'(a)',advance='no',iostat=io_status) &
         '        <DataArray type="'
    if (io_status==0) write(vti_unit,'(a)',advance='no', &
         iostat=io_status) trim(vtk_type)
    if (io_status==0) write(vti_unit,'(a)',advance='no', &
         iostat=io_status) '" Name="'
    if (io_status==0) write(vti_unit,'(a)',advance='no', &
         iostat=io_status) trim(name)
    if (io_status==0) write(vti_unit,'(a)',advance='no', &
         iostat=io_status) '" format="appended" offset="'
    if (io_status==0) write(vti_unit,'(i0)',advance='no', &
         iostat=io_status) offset
    if (io_status==0) write(vti_unit,'(a)',iostat=io_status) '"/>'
  end subroutine mt_write_vti_appended_array

  subroutine mt_write_vtu_length_pointdata(vtu_unit,length_results,npoint, &
       io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'length_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (length_results(ipoint)%total_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_backward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (length_results(ipoint)%backward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_forward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (length_results(ipoint)%forward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_int_array_start(vtu_unit,'nstep_backward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%backward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'nstep_forward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%forward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward_length', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (length_results(ipoint)%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_length_pointdata

  subroutine mt_write_vtu_mapping_pointdata(vtu_unit,mapping_results, &
       npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_mapping_result), intent(in) :: mapping_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'x_f_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(mapping_results(ipoint)%forward_footpoint,1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_f_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(mapping_results(ipoint)%forward_footpoint,2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_f_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(mapping_results(ipoint)%forward_footpoint,3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'x_b_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(mapping_results(ipoint)%backward_footpoint,1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_b_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(mapping_results(ipoint)%backward_footpoint,2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_b_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(mapping_results(ipoint)%backward_footpoint,3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'source_Bn_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%source_Bn,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'forward_Bn_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%forward_Bn,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'backward_Bn_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mapping_results(ipoint)%backward_Bn,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_int_array_start(vtu_unit,'face_forward_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mapping_results(ipoint)%forward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_backward_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mapping_results(ipoint)%backward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mapping_results(ipoint)%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward_mapping', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mapping_results(ipoint)%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'valid_mapping',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (merge(1,0,mapping_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_mapping_pointdata

  subroutine mt_write_vtu_qperp_public_pointdata(vtu_unit,qperp_results, &
       npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'logQperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_valid_float(qperp_results(ipoint)%logqperp, &
         qperp_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'valid_Qperp',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (merge(1,0,qperp_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_Qperp',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_qperp_public_pointdata

  subroutine mt_write_vtu_qperp_method2_pointdata(vtu_unit,qperp_results, &
       npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'qperp_method2',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%qperp,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'logqperp_method2', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%logqperp,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'N2_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%N2,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'bfactor_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%bfactor,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%forward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%backward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'Bseed_norm_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%B_seed**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'Bf_norm_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%forward_B**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'Bb_norm_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%backward_B**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'x_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%forward_endpoint,1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%forward_endpoint,2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%forward_endpoint,3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'x_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%backward_endpoint,1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%backward_endpoint,2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%backward_endpoint,3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_int_array_start(vtu_unit,'valid_qperp_method2', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (merge(1,0,qperp_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_qperp_method2', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%forward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%backward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_qperp_method2_pointdata

  subroutine mt_write_vtu_q_product_pointdata(vtu_unit,qperp_results, &
       npoint,write_raw_q,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    logical, intent(in) :: write_raw_q
    integer, intent(inout) :: io_status

    integer :: ipoint

    if (write_raw_q) then
      call mt_write_vtu_float_array_start(vtu_unit,'q',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (qperp_results(ipoint)%q0,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    call mt_write_vtu_float_array_start(vtu_unit,'logQ',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_valid_float(qperp_results(ipoint)%logq0, &
         qperp_results(ipoint)%valid_q0),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'valid_Q',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (merge(1,0,qperp_results(ipoint)%valid_q0),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_Q',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%status_q0,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_pair_Q',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mt_q_result_face_pair(qperp_results(ipoint)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'connection_type_Q', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mt_q_result_connection_type(qperp_results(ipoint)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_q_product_pointdata

  subroutine mt_write_vtu_spherical_topology_pointdata(vtu_unit,topology, &
       npoint,do_twist,do_q,q_results,do_qperp,qperp_results,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_topology_result), intent(in) :: topology(npoint)
    logical, intent(in) :: do_twist,do_q,do_qperp
    type(trace_qperp_result), intent(in) :: q_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    integer, intent(inout) :: io_status

    integer :: ipoint

    write(vtu_unit,'(a)',iostat=io_status) '<PointData>'
    if (io_status/=0) return

    if (.not.mt_vtk_detail_is_full()) then
      call mt_write_vtu_spherical_minimal_pointdata(vtu_unit,topology, &
           npoint,do_twist,do_q,q_results,do_qperp,qperp_results, &
           io_status)
      if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
           '</PointData>'
      return
    endif

    call mt_write_vtu_float_array_start(vtu_unit,'length_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_float(topology(ipoint)%length_total),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_backward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (topology(ipoint)%length_backward,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_forward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (topology(ipoint)%length_forward,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    if (do_twist) then
      call mt_write_vtu_float_array_start(vtu_unit,'twist_total',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (mt_visual_valid_float(topology(ipoint)%twist_total, &
           topology(ipoint)%valid_twist),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_float_array_start(vtu_unit,'twist_backward', &
           io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (mt_visual_valid_float(topology(ipoint)%twist_backward, &
           topology(ipoint)%valid_twist),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_float_array_start(vtu_unit,'twist_forward', &
           io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (mt_visual_valid_float(topology(ipoint)%twist_forward, &
           topology(ipoint)%valid_twist),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_int_array_start(vtu_unit,'valid_twist',io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (merge(1,0,topology(ipoint)%valid_twist),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_int_array_start(vtu_unit,'status_twist',io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (topology(ipoint)%status_twist,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    if (do_q) then
      call mt_write_vtu_float_array_start(vtu_unit,'logQ',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (mt_visual_valid_float(q_results(ipoint)%logq0, &
           q_results(ipoint)%valid_q0),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_int_array_start(vtu_unit,'valid_Q',io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (merge(1,0,q_results(ipoint)%valid_q0),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_int_array_start(vtu_unit,'status_Q',io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (q_results(ipoint)%status_q0,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    if (do_qperp) then
      call mt_write_vtu_float_array_start(vtu_unit,'logQperp',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (mt_visual_valid_float(qperp_results(ipoint)%logqperp, &
           qperp_results(ipoint)%valid),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_int_array_start(vtu_unit,'valid_Qperp',io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (merge(1,0,qperp_results(ipoint)%valid),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_int_array_start(vtu_unit,'status_Qperp',io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (qperp_results(ipoint)%status,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    call mt_write_vtu_float_array_start(vtu_unit,'r_b',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (topology(ipoint)%backward_endpoint(1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'theta_b',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(topology(ipoint)%backward_endpoint,2), &
         ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'phi_b',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(topology(ipoint)%backward_endpoint,3), &
         ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'r_f',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (topology(ipoint)%forward_endpoint(1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'theta_f',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(topology(ipoint)%forward_endpoint,2), &
         ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'phi_f',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(topology(ipoint)%forward_endpoint,3), &
         ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_int_array_start(vtu_unit,'face_backward',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (topology(ipoint)%backward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_forward',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (topology(ipoint)%forward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'connection_type_spherical', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (mt_spherical_connection_type(topology(ipoint)%backward_face, &
         topology(ipoint)%forward_face,topology(ipoint)%valid), &
         ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'nstep_backward',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (topology(ipoint)%backward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'nstep_forward',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (topology(ipoint)%forward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (topology(ipoint)%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (topology(ipoint)%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'valid',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (merge(1,0,topology(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</PointData>'
  end subroutine mt_write_vtu_spherical_topology_pointdata

  subroutine mt_write_vtu_spherical_minimal_pointdata(vtu_unit,topology, &
       npoint,do_twist,do_q,q_results,do_qperp,qperp_results,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_topology_result), intent(in) :: topology(npoint)
    logical, intent(in) :: do_twist,do_q,do_qperp
    type(trace_qperp_result), intent(in) :: q_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'length_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_visual_float(topology(ipoint)%length_total),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    if (do_twist) then
      call mt_write_vtu_float_array_start(vtu_unit,'twist_total',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (mt_visual_valid_float(topology(ipoint)%twist_total, &
           topology(ipoint)%valid_twist),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    if (do_q) then
      call mt_write_vtu_float_array_start(vtu_unit,'logQ',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (mt_visual_valid_float(q_results(ipoint)%logq0, &
           q_results(ipoint)%valid_q0),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    if (do_qperp) then
      call mt_write_vtu_float_array_start(vtu_unit,'logQperp',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (mt_visual_valid_float(qperp_results(ipoint)%logqperp, &
           qperp_results(ipoint)%valid),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif
  end subroutine mt_write_vtu_spherical_minimal_pointdata

  subroutine mt_write_vtu_file_header(vtu_unit,npoint,ncell)
    integer, intent(in) :: vtu_unit,npoint,ncell

    write(vtu_unit,'(a)') '<?xml version="1.0"?>'
    write(vtu_unit,'(a)') &
         '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(vtu_unit,'(a)') '<UnstructuredGrid>'
    write(vtu_unit,'(a,i0,a,i0,a)') '<Piece NumberOfPoints="',npoint, &
         '" NumberOfCells="',ncell,'">'
  end subroutine mt_write_vtu_file_header

  subroutine mt_write_vtu_file_footer(vtu_unit,io_status)
    integer, intent(in) :: vtu_unit
    integer, intent(inout) :: io_status

    write(vtu_unit,'(a)',iostat=io_status) '</Piece>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</UnstructuredGrid>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</VTKFile>'
  end subroutine mt_write_vtu_file_footer

  subroutine mt_write_vtu_product_pointdata(vtu_unit,length_results, &
       twist_results,q_results,qperp_results,npoint,do_twist,do_q, &
       do_qperp,io_status,do_length)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: q_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    logical, intent(in) :: do_twist,do_q,do_qperp
    integer, intent(inout) :: io_status
    logical, intent(in), optional :: do_length

    integer :: ipoint
    logical :: do_length_eff

    write(vtu_unit,'(a)',iostat=io_status) '<PointData>'
    if (io_status/=0) return
    do_length_eff=.true.
    if (present(do_length)) do_length_eff=do_length

    if (.not.mt_vtk_detail_is_full()) then
      call mt_write_vtu_product_minimal_pointdata(vtu_unit,length_results, &
           twist_results,q_results,qperp_results,npoint,do_twist,do_q, &
           do_qperp,io_status,do_length=do_length_eff)
      if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
           '</PointData>'
      return
    endif

    if (do_length_eff) then
      call mt_write_vtu_float_array_start(vtu_unit,'length_total',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (length_results(ipoint)%total_length,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_float_array_start(vtu_unit,'length_backward', &
           io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (length_results(ipoint)%backward_length,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_float_array_start(vtu_unit,'length_forward', &
           io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (length_results(ipoint)%forward_length,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)

      call mt_write_vtu_int_array_start(vtu_unit,'nstep_backward_length', &
           io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (length_results(ipoint)%backward_nstep,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_int_array_start(vtu_unit,'nstep_forward_length', &
           io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (length_results(ipoint)%forward_nstep,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_int_array_start(vtu_unit,'status_backward_length', &
           io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (length_results(ipoint)%backward_status,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
      call mt_write_vtu_int_array_start(vtu_unit,'status_forward_length', &
           io_status)
      if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
           (length_results(ipoint)%forward_status,ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    if (do_twist) then
      call mt_write_vtu_twist_pointdata(vtu_unit,twist_results,npoint, &
           io_status)
    endif
    if (do_q) then
      call mt_write_vtu_q_product_pointdata(vtu_unit,q_results,npoint, &
           .false.,io_status)
    endif
    if (do_qperp) then
      call mt_write_vtu_qperp_pointdata(vtu_unit,qperp_results,npoint, &
           io_status)
    endif

    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</PointData>'
  end subroutine mt_write_vtu_product_pointdata

  subroutine mt_write_vtu_product_minimal_pointdata(vtu_unit, &
       length_results,twist_results,q_results,qperp_results,npoint, &
       do_twist,do_q,do_qperp,io_status,do_length)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    type(trace_twist_result), intent(in) :: twist_results(:)
    type(trace_qperp_result), intent(in) :: q_results(:)
    type(trace_qperp_result), intent(in) :: qperp_results(:)
    logical, intent(in) :: do_twist,do_q,do_qperp
    integer, intent(inout) :: io_status
    logical, intent(in), optional :: do_length

    integer :: ipoint
    logical :: do_length_eff

    do_length_eff=.true.
    if (present(do_length)) do_length_eff=do_length

    if (do_length_eff) then
      call mt_write_vtu_float_array_start(vtu_unit,'length_total',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
           (mt_visual_float(length_results(ipoint)%total_length), &
           ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    if (do_twist) then
      call mt_write_vtu_float_array_start(vtu_unit,'twist_total',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))', &
           iostat=io_status) &
           (mt_visual_float(twist_results(ipoint)%total_twist), &
           ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    if (do_q) then
      call mt_write_vtu_float_array_start(vtu_unit,'logQ',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))', &
           iostat=io_status) &
           (mt_visual_valid_float(q_results(ipoint)%logq0, &
           q_results(ipoint)%valid_q0),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif

    if (do_qperp) then
      call mt_write_vtu_float_array_start(vtu_unit,'logQperp',io_status)
      if (io_status==0) write(vtu_unit,'(4(1pe24.16))', &
           iostat=io_status) &
           (mt_visual_valid_float(qperp_results(ipoint)%logqperp, &
           qperp_results(ipoint)%valid),ipoint=1,npoint)
      call mt_write_vtu_data_array_end(vtu_unit,io_status)
    endif
  end subroutine mt_write_vtu_product_minimal_pointdata

  subroutine mt_write_vtu_twist_pointdata(vtu_unit,twist_results,npoint, &
       io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_twist_result), intent(in) :: twist_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'twist_total',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (twist_results(ipoint)%total_twist,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'twist_backward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (twist_results(ipoint)%backward_twist,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'twist_forward', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (twist_results(ipoint)%forward_twist,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_int_array_start(vtu_unit,'nstep_backward_twist', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (twist_results(ipoint)%line%backward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'nstep_forward_twist', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (twist_results(ipoint)%line%forward_nstep,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward_twist', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (twist_results(ipoint)%line%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward_twist', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (twist_results(ipoint)%line%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_twist_pointdata

  subroutine mt_write_vtu_qperp_pointdata(vtu_unit,qperp_results,npoint, &
       io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%qperp,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'logqperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%logqperp,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'N2',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%N2,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'bfactor',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%bfactor,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%forward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'length_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (qperp_results(ipoint)%backward_length,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_float_array_start(vtu_unit,'Bseed_norm',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%B_seed**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'Bf_norm',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%forward_B**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'Bb_norm',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (dsqrt(sum(qperp_results(ipoint)%backward_B**2)),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)

    call mt_write_vtu_qperp_endpoint_pointdata(vtu_unit,qperp_results, &
         npoint,io_status)
    call mt_write_vtu_qperp_int_pointdata(vtu_unit,qperp_results,npoint, &
         io_status)
  end subroutine mt_write_vtu_qperp_pointdata

  subroutine mt_write_vtu_qperp_endpoint_pointdata(vtu_unit,qperp_results, &
       npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_float_array_start(vtu_unit,'x_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%forward_endpoint,1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%forward_endpoint,2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_f_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%forward_endpoint,3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'x_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%backward_endpoint,1),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'y_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%backward_endpoint,2),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_float_array_start(vtu_unit,'z_b_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(4(1pe24.16))',iostat=io_status) &
         (mt_vc(qperp_results(ipoint)%backward_endpoint,3),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_qperp_endpoint_pointdata

  subroutine mt_write_vtu_qperp_int_pointdata(vtu_unit,qperp_results, &
       npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_qperp_result), intent(in) :: qperp_results(npoint)
    integer, intent(inout) :: io_status

    integer :: ipoint

    call mt_write_vtu_int_array_start(vtu_unit,'valid_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (merge(1,0,qperp_results(ipoint)%valid),ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_qperp',io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%forward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'face_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%backward_face,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_forward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%forward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
    call mt_write_vtu_int_array_start(vtu_unit,'status_backward_qperp', &
         io_status)
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (qperp_results(ipoint)%backward_status,ipoint=1,npoint)
    call mt_write_vtu_data_array_end(vtu_unit,io_status)
  end subroutine mt_write_vtu_qperp_int_pointdata

  subroutine mt_write_vtu_empty_pointdata(vtu_unit,do_length,do_twist, &
       do_q,do_qperp)
    integer, intent(in) :: vtu_unit
    logical, intent(in) :: do_length,do_twist,do_q,do_qperp

    write(vtu_unit,'(a)') '<PointData>'
    if (.not.mt_vtk_detail_is_full()) then
      if (do_length) then
        call mt_write_vtu_empty_float_array(vtu_unit,'length_total')
      endif
      if (do_twist) then
        call mt_write_vtu_empty_float_array(vtu_unit,'twist_total')
      endif
      if (do_q) then
        call mt_write_vtu_empty_float_array(vtu_unit,'logQ')
      endif
      if (do_qperp) then
        call mt_write_vtu_empty_float_array(vtu_unit,'logQperp')
      endif
      write(vtu_unit,'(a)') '</PointData>'
      return
    endif
    if (do_length) then
      call mt_write_vtu_empty_float_array(vtu_unit,'length_total')
      call mt_write_vtu_empty_float_array(vtu_unit,'length_backward')
      call mt_write_vtu_empty_float_array(vtu_unit,'length_forward')
      call mt_write_vtu_empty_int_array(vtu_unit,'nstep_backward_length')
      call mt_write_vtu_empty_int_array(vtu_unit,'nstep_forward_length')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_backward_length')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_forward_length')
    endif
    if (do_twist) then
      call mt_write_vtu_empty_float_array(vtu_unit,'twist_total')
      call mt_write_vtu_empty_float_array(vtu_unit,'twist_backward')
      call mt_write_vtu_empty_float_array(vtu_unit,'twist_forward')
      call mt_write_vtu_empty_int_array(vtu_unit,'nstep_backward_twist')
      call mt_write_vtu_empty_int_array(vtu_unit,'nstep_forward_twist')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_backward_twist')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_forward_twist')
    endif
    if (do_q) then
      call mt_write_vtu_empty_float_array(vtu_unit,'q')
      call mt_write_vtu_empty_float_array(vtu_unit,'logQ')
      call mt_write_vtu_empty_float_array(vtu_unit,'N2_q')
      call mt_write_vtu_empty_float_array(vtu_unit,'bfactor_q')
      call mt_write_vtu_empty_float_array(vtu_unit,'length_forward_q')
      call mt_write_vtu_empty_float_array(vtu_unit,'length_backward_q')
      call mt_write_vtu_empty_float_array(vtu_unit,'Bseed_norm_q')
      call mt_write_vtu_empty_float_array(vtu_unit,'Bf_norm_q')
      call mt_write_vtu_empty_float_array(vtu_unit,'Bb_norm_q')
      call mt_write_vtu_empty_int_array(vtu_unit,'valid_Q')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_Q')
      call mt_write_vtu_empty_int_array(vtu_unit,'face_forward_Q')
      call mt_write_vtu_empty_int_array(vtu_unit,'face_backward_Q')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_forward_Q')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_backward_Q')
    endif
    if (do_qperp) then
      call mt_write_vtu_empty_float_array(vtu_unit,'qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'logqperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'N2')
      call mt_write_vtu_empty_float_array(vtu_unit,'bfactor')
      call mt_write_vtu_empty_float_array(vtu_unit,'length_forward_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'length_backward_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'Bseed_norm')
      call mt_write_vtu_empty_float_array(vtu_unit,'Bf_norm')
      call mt_write_vtu_empty_float_array(vtu_unit,'Bb_norm')
      call mt_write_vtu_empty_float_array(vtu_unit,'x_f_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'y_f_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'z_f_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'x_b_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'y_b_qperp')
      call mt_write_vtu_empty_float_array(vtu_unit,'z_b_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'valid_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'face_forward_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'face_backward_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_forward_qperp')
      call mt_write_vtu_empty_int_array(vtu_unit,'status_backward_qperp')
    endif
    write(vtu_unit,'(a)') '</PointData>'
  end subroutine mt_write_vtu_empty_pointdata

  subroutine mt_write_vtu_empty_float_array(vtu_unit,name)
    integer, intent(in) :: vtu_unit
    character(len=*), intent(in) :: name

    write(vtu_unit,'(a,a,a)') '<DataArray type="Float64" Name="', &
         trim(name),'" format="ascii">'
    write(vtu_unit,'(a)') '</DataArray>'
  end subroutine mt_write_vtu_empty_float_array

  subroutine mt_write_vtu_empty_int_array(vtu_unit,name)
    integer, intent(in) :: vtu_unit
    character(len=*), intent(in) :: name

    write(vtu_unit,'(a,a,a)') '<DataArray type="Int32" Name="', &
         trim(name),'" format="ascii">'
    write(vtu_unit,'(a)') '</DataArray>'
  end subroutine mt_write_vtu_empty_int_array

  subroutine mt_write_vtu_float_array_start(vtu_unit,name,io_status)
    integer, intent(in) :: vtu_unit
    character(len=*), intent(in) :: name
    integer, intent(inout) :: io_status

    if (io_status/=0) return
    write(vtu_unit,'(a,a,a)',iostat=io_status) &
         '<DataArray type="Float64" Name="',trim(name),'" format="ascii">'
  end subroutine mt_write_vtu_float_array_start

  subroutine mt_write_vtu_int_array_start(vtu_unit,name,io_status)
    integer, intent(in) :: vtu_unit
    character(len=*), intent(in) :: name
    integer, intent(inout) :: io_status

    if (io_status/=0) return
    write(vtu_unit,'(a,a,a)',iostat=io_status) &
         '<DataArray type="Int32" Name="',trim(name),'" format="ascii">'
  end subroutine mt_write_vtu_int_array_start

  subroutine mt_write_vtu_data_array_end(vtu_unit,io_status)
    integer, intent(in) :: vtu_unit
    integer, intent(inout) :: io_status

    if (io_status/=0) return
    write(vtu_unit,'(a)',iostat=io_status) '</DataArray>'
  end subroutine mt_write_vtu_data_array_end

  subroutine mt_write_vtu_points(vtu_unit,length_results,npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_length_result), intent(in) :: length_results(npoint)
    integer, intent(inout) :: io_status

    double precision :: seed_xyz(3)
    integer :: ipoint

    write(vtu_unit,'(a)',iostat=io_status) '<Points>'
    if (io_status==0) then
      write(vtu_unit,'(a)',iostat=io_status) &
           '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    endif
    do ipoint=1,npoint
      call mt_vtu_point_from_coord(length_results(ipoint)%seed,seed_xyz)
      if (io_status==0) write(vtu_unit,'(3(1pe24.16))', &
           iostat=io_status) seed_xyz
    enddo
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</Points>'
  end subroutine mt_write_vtu_points

  subroutine mt_write_vtu_topology_points(vtu_unit,topology,npoint, &
       io_status)
    integer, intent(in) :: vtu_unit,npoint
    type(trace_topology_result), intent(in) :: topology(npoint)
    integer, intent(inout) :: io_status

    double precision :: seed_xyz(3)
    integer :: ipoint

    write(vtu_unit,'(a)',iostat=io_status) '<Points>'
    if (io_status==0) then
      write(vtu_unit,'(a)',iostat=io_status) &
           '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    endif
    do ipoint=1,npoint
      call mt_vtu_point_from_coord(topology(ipoint)%seed,seed_xyz)
      if (io_status==0) write(vtu_unit,'(3(1pe24.16))', &
           iostat=io_status) seed_xyz
    enddo
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</Points>'
  end subroutine mt_write_vtu_topology_points

  subroutine mt_vtu_point_from_coord(coord,point_xyz)
    double precision, intent(in) :: coord(ndim)
    double precision, intent(out) :: point_xyz(3)

    double precision :: r,theta,phi,sin_theta

    point_xyz=0.d0
    if (geo_coordinate==geo_spherical .and. ndim==3) then
      {^IFTHREED
      r=coord(1)
      theta=coord(2)
      phi=coord(3)
      sin_theta=dsin(theta)
      point_xyz(1)=r*sin_theta*dcos(phi)
      point_xyz(2)=r*sin_theta*dsin(phi)
      point_xyz(3)=r*dcos(theta)
      }
    else
      point_xyz(1)=mt_vc(coord,1)
      point_xyz(2)=mt_vc(coord,2)
      point_xyz(3)=mt_vc(coord,3)
    endif
  end subroutine mt_vtu_point_from_coord

  subroutine mt_write_vtu_vertex_cells(vtu_unit,npoint,io_status)
    integer, intent(in) :: vtu_unit,npoint
    integer, intent(inout) :: io_status

    integer :: ipoint

    write(vtu_unit,'(a)',iostat=io_status) '<Cells>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="connectivity" format="ascii">'
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (ipoint-1,ipoint=1,npoint)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="offsets" format="ascii">'
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (ipoint,ipoint=1,npoint)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="UInt8" Name="types" format="ascii">'
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (1,ipoint=1,npoint)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</Cells>'
  end subroutine mt_write_vtu_vertex_cells

  subroutine mt_write_vtu_quad_cells(vtu_unit,n1,n2,io_status)
    integer, intent(in) :: vtu_unit,n1,n2
    integer, intent(inout) :: io_status

    integer :: i,j,icell,p,ncell

    ncell=(n1-1)*(n2-1)
    write(vtu_unit,'(a)',iostat=io_status) '<Cells>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="connectivity" format="ascii">'
    do j=1,n2-1
      do i=1,n1-1
        p=(j-1)*n1+(i-1)
        if (io_status==0) write(vtu_unit,'(4(i0,1x))',iostat=io_status) &
             p,p+1,p+1+n1,p+n1
      enddo
    enddo
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="offsets" format="ascii">'
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (4*icell,icell=1,ncell)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="UInt8" Name="types" format="ascii">'
    if (io_status==0) write(vtu_unit,'(12(i0,1x))',iostat=io_status) &
         (9,icell=1,ncell)
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</Cells>'
  end subroutine mt_write_vtu_quad_cells

  subroutine mt_write_vtu_cells_empty(vtu_unit,io_status)
    integer, intent(in) :: vtu_unit
    integer, intent(inout) :: io_status

    write(vtu_unit,'(a)',iostat=io_status) '<Cells>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="connectivity" format="ascii">'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="Int32" Name="offsets" format="ascii">'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '<DataArray type="UInt8" Name="types" format="ascii">'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) &
         '</DataArray>'
    if (io_status==0) write(vtu_unit,'(a)',iostat=io_status) '</Cells>'
  end subroutine mt_write_vtu_cells_empty

  subroutine mt_write_twist_plane_csv(results,n1,n2,csv_file,caller, &
       index_header)
    integer, intent(in) :: n1,n2
    type(trace_twist_result), intent(in) :: results(n1*n2)
    character(len=*), intent(in) :: csv_file,caller,index_header

    double precision :: seed_xyz(3)
    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         trim(index_header)//',seed_x,seed_y,seed_z,'// &
         'length_total,length_backward,length_forward,'// &
         'twist_total,twist_backward,twist_forward,'// &
         'nstep_backward,nstep_forward,'// &
         'status_backward,status_forward'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        seed_xyz=0.d0
        seed_xyz(1:ndim)=results(iseed)%line%seed
        write(csv_unit,'(i0,",",i0,9(",",es24.16),4(",",i0))', &
             iostat=io_status) i,j,seed_xyz, &
             results(iseed)%line%total_length, &
             results(iseed)%line%backward_length, &
             results(iseed)%line%forward_length, &
             results(iseed)%total_twist, &
             results(iseed)%backward_twist, &
             results(iseed)%forward_twist, &
             results(iseed)%line%backward_nstep, &
             results(iseed)%line%forward_nstep, &
             results(iseed)%line%backward_status, &
             results(iseed)%line%forward_status
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_twist_plane_csv

  subroutine mt_write_q_plane_csv(results,n1,n2,csv_file,caller, &
       index_header)
    integer, intent(in) :: n1,n2
    type(trace_qperp_result), intent(in) :: results(n1*n2)
    character(len=*), intent(in) :: csv_file,caller,index_header

    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         trim(index_header)//',seed_x,seed_y,seed_z,'// &
         'logQ,valid_Q,status_Q,'// &
         'face_forward_Q,face_backward_Q,'// &
         'status_forward_Q,status_backward_Q,'// &
         'length_forward_Q,length_backward_Q,'// &
         'x_f_Q,y_f_Q,z_f_Q,x_b_Q,y_b_Q,z_b_Q'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        write(csv_unit, &
             '(i0,",",i0,4(",",es24.16),",",l1,5(",",i0),'// &
             '8(",",es24.16))',iostat=io_status) &
             i,j,results(iseed)%seed, &
             results(iseed)%logq0, &
             results(iseed)%valid_q0,results(iseed)%status_q0, &
             results(iseed)%forward_face,results(iseed)%backward_face, &
             results(iseed)%forward_status,results(iseed)%backward_status, &
             results(iseed)%forward_length,results(iseed)%backward_length, &
             results(iseed)%forward_endpoint, &
             results(iseed)%backward_endpoint
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_q_plane_csv

  subroutine mt_write_qperp_plane_csv(results,n1,n2,csv_file,caller, &
       index_header)
    integer, intent(in) :: n1,n2
    type(trace_qperp_result), intent(in) :: results(n1*n2)
    character(len=*), intent(in) :: csv_file,caller,index_header

    double precision :: Bseed_norm,Bf_norm,Bb_norm
    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         trim(index_header)//',seed_x,seed_y,seed_z,'// &
         'qperp,logqperp,valid,status,'// &
         'N2,bfactor,'// &
         'face_forward,face_backward,'// &
         'status_forward,status_backward,'// &
         'length_forward,length_backward,'// &
         'Bseed_norm,Bf_norm,Bb_norm,'// &
         'x_f,y_f,z_f,x_b,y_b,z_b'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        Bseed_norm=dsqrt(sum(results(iseed)%B_seed**2))
        Bf_norm=dsqrt(sum(results(iseed)%forward_B**2))
        Bb_norm=dsqrt(sum(results(iseed)%backward_B**2))
        write(csv_unit, &
             '(i0,",",i0,5(",",es24.16),",",l1,",",i0,'// &
             '2(",",es24.16),4(",",i0),11(",",es24.16))', &
             iostat=io_status) &
             i,j,results(iseed)%seed, &
             results(iseed)%qperp,results(iseed)%logqperp, &
             results(iseed)%valid,results(iseed)%status, &
             results(iseed)%N2,results(iseed)%bfactor, &
             results(iseed)%forward_face,results(iseed)%backward_face, &
             results(iseed)%forward_status,results(iseed)%backward_status, &
             results(iseed)%forward_length,results(iseed)%backward_length, &
             Bseed_norm,Bf_norm,Bb_norm, &
             results(iseed)%forward_endpoint, &
             results(iseed)%backward_endpoint
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_qperp_plane_csv

  subroutine mt_write_qperp_arbitrary_header(csv_file,caller)
    character(len=*), intent(in) :: csv_file,caller

    integer :: csv_unit,io_status

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'i,j,s1,s2,seed_x,seed_y,seed_z,'// &
         'qperp,logqperp,valid,status,'// &
         'N2,bfactor,'// &
         'face_forward,face_backward,'// &
         'status_forward,status_backward,'// &
         'length_forward,length_backward,'// &
         'Bseed_norm,Bf_norm,Bb_norm,'// &
         'x_f,y_f,z_f,x_b,y_b,z_b'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    close(csv_unit)
  end subroutine mt_write_qperp_arbitrary_header

  subroutine mt_write_qperp_arbitrary_csv(results,s1,s2,n1,n2,csv_file, &
       caller)
    integer, intent(in) :: n1,n2
    type(trace_qperp_result), intent(in) :: results(n1*n2)
    double precision, intent(in) :: s1(n1*n2),s2(n1*n2)
    character(len=*), intent(in) :: csv_file,caller

    double precision :: Bseed_norm,Bf_norm,Bb_norm
    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         'i,j,s1,s2,seed_x,seed_y,seed_z,'// &
         'qperp,logqperp,valid,status,'// &
         'N2,bfactor,'// &
         'face_forward,face_backward,'// &
         'status_forward,status_backward,'// &
         'length_forward,length_backward,'// &
         'Bseed_norm,Bf_norm,Bb_norm,'// &
         'x_f,y_f,z_f,x_b,y_b,z_b'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        Bseed_norm=dsqrt(sum(results(iseed)%B_seed**2))
        Bf_norm=dsqrt(sum(results(iseed)%forward_B**2))
        Bb_norm=dsqrt(sum(results(iseed)%backward_B**2))
        write(csv_unit, &
             '(i0,",",i0,7(",",es24.16),",",l1,",",i0,'// &
             '2(",",es24.16),4(",",i0),11(",",es24.16))', &
             iostat=io_status) &
             i,j,s1(iseed),s2(iseed),results(iseed)%seed, &
             results(iseed)%qperp,results(iseed)%logqperp, &
             results(iseed)%valid,results(iseed)%status, &
             results(iseed)%N2,results(iseed)%bfactor, &
             results(iseed)%forward_face,results(iseed)%backward_face, &
             results(iseed)%forward_status,results(iseed)%backward_status, &
             results(iseed)%forward_length,results(iseed)%backward_length, &
             Bseed_norm,Bf_norm,Bb_norm, &
             results(iseed)%forward_endpoint, &
             results(iseed)%backward_endpoint
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_qperp_arbitrary_csv

  subroutine mt_write_mapping_plane_xy_csv(results,nx,ny,csv_file)
    integer, intent(in) :: nx,ny
    type(trace_mapping_result), intent(in) :: results(nx*ny)
    character(len=*), intent(in) :: csv_file

    call mt_write_mapping_plane_csv(results,nx,ny,csv_file, &
         'mt_mapping_plane_xy','ix,iy')
  end subroutine mt_write_mapping_plane_xy_csv

  subroutine mt_write_mapping_plane_csv(results,n1,n2,csv_file,caller, &
       index_header)
    integer, intent(in) :: n1,n2
    type(trace_mapping_result), intent(in) :: results(n1*n2)
    character(len=*), intent(in) :: csv_file,caller,index_header

    integer :: csv_unit,io_status,i,j,iseed

    open(newunit=csv_unit,file=trim(csv_file),status='replace', &
         action='write',form='formatted',iostat=io_status)
    if (io_status/=0) then
      call mpistop(trim(caller)//' could not open CSV file')
    endif

    write(csv_unit,'(a)',iostat=io_status) &
         trim(index_header)//','// &
         'seed_x,seed_y,seed_z,'// &
         'source_Bx,source_By,source_Bz,source_Bn,'// &
         'backward_x,backward_y,backward_z,'// &
         'backward_Bx,backward_By,backward_Bz,backward_Bn,'// &
         'backward_face,backward_length,backward_status,'// &
         'forward_x,forward_y,forward_z,'// &
         'forward_Bx,forward_By,forward_Bz,forward_Bn,'// &
         'forward_face,forward_length,forward_status,valid'
    if (io_status/=0) then
      close(csv_unit)
      call mpistop(trim(caller)//' could not write CSV header')
    endif

    do j=1,n2
      do i=1,n1
        iseed=(j-1)*n1+i
        write(csv_unit, &
             '(i0,",",i0,14(",",es24.16),",",i0,",",es24.16,'// &
             '",",i0,7(",",es24.16),",",i0,",",es24.16,'// &
             '",",i0,",",l1)',iostat=io_status) &
             i,j, &
             results(iseed)%seed, &
             results(iseed)%source_B,results(iseed)%source_Bn, &
             results(iseed)%backward_footpoint, &
             results(iseed)%backward_B,results(iseed)%backward_Bn, &
             results(iseed)%backward_face, &
             results(iseed)%backward_length, &
             results(iseed)%backward_status, &
             results(iseed)%forward_footpoint, &
             results(iseed)%forward_B,results(iseed)%forward_Bn, &
             results(iseed)%forward_face, &
             results(iseed)%forward_length, &
             results(iseed)%forward_status, &
             results(iseed)%valid
        if (io_status/=0) then
          close(csv_unit)
          call mpistop(trim(caller)//' could not write CSV data')
        endif
      enddo
    enddo

    close(csv_unit)
  end subroutine mt_write_mapping_plane_csv

  integer function mt_q_result_face_pair(qperp_result) result(face_pair)
    type(trace_qperp_result), intent(in) :: qperp_result

    face_pair=0
    if (.not.qperp_result%valid_q0) return
    face_pair=mt_q_face_pair_from_faces(qperp_result%backward_face, &
         qperp_result%forward_face)
  end function mt_q_result_face_pair

  integer function mt_q_result_connection_type(qperp_result) &
       result(connection_type)
    type(trace_qperp_result), intent(in) :: qperp_result

    connection_type=0
    if (.not.qperp_result%valid_q0) return
    connection_type=mt_q_connection_type_from_faces( &
         qperp_result%backward_face,qperp_result%forward_face)
  end function mt_q_result_connection_type

  integer function mt_q_face_pair_from_faces(face_b,face_f) result(face_pair)
    integer, intent(in) :: face_b,face_f

    face_pair=0
    if (.not.mt_q_face_valid(face_b)) return
    if (.not.mt_q_face_valid(face_f)) return
    face_pair=10*face_b+face_f
  end function mt_q_face_pair_from_faces

  integer function mt_q_connection_type_from_faces(face_b,face_f) &
       result(connection_type)
    integer, intent(in) :: face_b,face_f

    connection_type=0
    if (.not.mt_q_face_valid(face_b)) return
    if (.not.mt_q_face_valid(face_f)) return
    if (face_b==trace_face_zmin .and. face_f==trace_face_zmin) then
      connection_type=1
    else if (face_b==trace_face_zmax .and. face_f==trace_face_zmax) then
      connection_type=2
    else if ((face_b==trace_face_zmin .and. face_f==trace_face_zmax) .or. &
             (face_b==trace_face_zmax .and. face_f==trace_face_zmin)) then
      connection_type=3
    else if (face_b==face_f) then
      connection_type=4
    else
      connection_type=5
    endif
  end function mt_q_connection_type_from_faces

  logical function mt_q_face_valid(face)
    integer, intent(in) :: face

    mt_q_face_valid=face>=trace_face_xmin .and. face<=trace_face_zmax
  end function mt_q_face_valid

end module mod_magnetic_topology
