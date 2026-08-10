module mod_trace_field
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan, &
       ieee_is_finite
  use mod_global_parameters
  use mod_geometry, only: geo_coordinate => coordinate, &
       geo_cartesian => Cartesian, geo_spherical => spherical, &
       geo_cartesian_stretched => Cartesian_stretched, &
       curlvector
  use mod_physics
  use mod_connectivity, only: igrids,igridstail
  implicit none

  integer, parameter, public :: trace_status_active               = 0
  integer, parameter, public :: trace_status_boundary             = 1
  integer, parameter, public :: trace_status_weak_field           = 2
  integer, parameter, public :: trace_status_max_steps            = 3
  integer, parameter, public :: trace_status_seed_outside         = 4
  integer, parameter, public :: trace_status_out_of_domain        = 5
  integer, parameter, public :: trace_status_invalid_input        = 6
  integer, parameter, public :: trace_status_unsupported_geometry = 7
  integer, parameter, public :: trace_status_trac_stop            = 8
  integer, parameter, public :: trace_status_mpi_unsupported      = 9
  integer, parameter, public :: trace_status_bad_curl_stencil     = 10
  integer, parameter, public :: trace_status_bad_grad_stencil     = 11
  integer, parameter, public :: trace_status_bad_face_limit_sample= 12
  integer, parameter, public :: trace_status_singular_q           = 13
  integer, parameter, public :: trace_status_bad_q_bound          = 14

  ! Both standard and perpendicular squashing factors have the analytic
  ! lower bound Q >= 2.  Cartesian tangent transport can undershoot this
  ! bound slightly because of integration error, so valid outputs are floored
  ! before taking the logarithm.
  double precision, parameter, private :: trace_q_min=2.d0

  integer, parameter, public :: trace_face_none      = 0
  integer, parameter, public :: trace_face_xmin      = 1
  integer, parameter, public :: trace_face_xmax      = 2
  integer, parameter, public :: trace_face_ymin      = 3
  integer, parameter, public :: trace_face_ymax      = 4
  integer, parameter, public :: trace_face_zmin      = 5
  integer, parameter, public :: trace_face_zmax      = 6
  integer, parameter, public :: trace_face_ambiguous = 7

  public :: trace_field_length_single,trace_field_length_multi
  public :: trace_field_twist_single,trace_field_twist_multi
  public :: trace_field_mapping_single,trace_field_mapping_multi
  public :: trace_field_topology_multi
  public :: trace_field_qperp_single,trace_field_qperp_multi
  public :: trace_field_spherical_rmin_q_multi
  public :: trace_field_rk2_short_boundary_q_multi
  public :: trace_field_spherical_rmin_q_qperp_multi
  public :: trace_field_spherical_qperp_multi
  public :: trace_spherical_curl_cache_build
  public :: trace_spherical_curl_cache_clear
  public :: trace_cartesian_global_min_cell_size
  public :: trace_spherical_global_min_cell_size
  public :: trace_set_step_control
  public :: trace_set_integrator
  public :: trace_rk2_stats_set_enabled
  public :: trace_rk2_stats_reset
  public :: trace_rk2_stats_report
  public :: trace_rk45_stats_set_enabled
  public :: trace_rk45_stats_reset
  public :: trace_rk45_stats_report
  public :: trace_spherical_profile_set
  public :: trace_spherical_profile_reset
  public :: trace_spherical_profile_count_seeds
  public :: trace_spherical_profile_report
  public :: trace_debug_sample_bhat_gradbhat
  public :: trace_debug_compare_gradbhat
  public :: trace_debug_compare_gradbhat_methods
  public :: trace_debug_transport_tangents
  public :: trace_debug_transport_tangents_to_boundary
  public :: trace_debug_sample_endpoint_B_face_limit
  public :: trace_debug_qperp_single
  public :: trace_field_q0_multi_rk45_cartesian
  public :: trace_debug_cartesian_rk45_tangent_q0_multi
  private :: trace_summary_seed,trace_summary_multi
  private :: trace_summary_twist_seed,trace_summary_twist_multi
  private :: trace_summary_mapping_seed,trace_summary_mapping_multi
  private :: trace_summary_topology_multi
  private :: trace_summary_validate,trace_summary_init_result
  private :: trace_summary_init_twist_result
  private :: trace_summary_init_mapping_result
  private :: trace_summary_init_topology_result
  private :: trace_summary_locate_seed,trace_summary_trace_seed
  private :: trace_summary_trace_twist_seed
  private :: trace_summary_init_state,trace_summary_advance_state
  private :: trace_summary_trace_state_to_end,trace_summary_trace_states_grouped
  private :: trace_summary_advance_state_in_grid
  private :: trace_summary_rk2_step,trace_summary_rk45_cartesian_step
  private :: trace_summary_rk45_spherical_step
  private :: trace_intersect_domain
  private :: trace_summary_accumulate_twist,trace_twist_density_at_point
  private :: trace_tangent_accumulate_twist,trace_tangent_fill_twist_result
  private :: sample_B_curlB_at_point
  private :: sample_B_curlB_spherical_at_point
  private :: sample_B_at_point,trace_summary_fill_mapping
  private :: trace_total_B_at_cell,trace_bhat_at_cell
  private :: sample_bhat_gradbhat_at_point
  private :: sample_bhat_gradbhat_cellfd_at_point
  private :: sample_bhat_gradbhat_xeps_at_point
  private :: sample_bhat_gradbhat_interpderiv_at_point
  private :: trace_debug_locate_point,trace_locate_point_with_hint
  private :: trace_tangent_rhs
  private :: trace_advance_tangent_state_rk2,trace_project_to_perp_bhat
  private :: trace_tangent_rk2_trial,trace_endpoint_B_bhat
  private :: trace_tangent_rk2_trial_from_rhs
  private :: trace_tangent_state_advance_in_grid_rk45_cartesian
  private :: trace_tangent_rk45_cartesian_step
  private :: trace_tangent_rhs_cartesian_located
  private :: trace_tangent_rk45_cartesian_finish_boundary
  private :: trace_tangent_trace_states_grouped_rk45_spherical
  private :: trace_tangent_state_advance_in_grid_rk45_spherical
  private :: trace_tangent_rk45_spherical_step
  private :: trace_tangent_rk45_spherical_finish_boundary
  private :: trace_q0_finalize_from_states
  private :: trace_sample_B_on_domain_face_limit
  private :: trace_locate_face_limit_grid
  private :: trace_face_normal,trace_project_to_boundary_face
  private :: trace_make_perp_basis,trace_debug_nan
  private :: trace_face_from_dim_side,trace_face_Bn,trace_face_is_boundary
  private :: trace_boundary_face_at_point
  private :: trace_is_supported_summary_geometry
  private :: trace_effective_step,trace_spherical_metric_ok
  private :: trace_segment_length
  private :: trace_summary_twist_status_from_states
  private :: trace_init_qperp_result,trace_qperp_single_core
  private :: trace_qperp_compute_q0_scalars
  private :: trace_qperp_project_to_boundary_face
  private :: trace_tangent_state_init
  private :: trace_tangent_trace_states_grouped
  private :: trace_tangent_state_advance_in_grid
  private :: trace_tangent_state_finalize_boundary
  private :: trace_qperp_prepare_seed_result
  private :: trace_qperp_compute_scalars
  private :: trace_qperp_finalize_from_states
  private :: trace_spherical_rmin_q_prepare_seed_result
  private :: trace_spherical_rmin_q_finalize_from_states
  private :: trace_spherical_radial_q_finalize_from_states
  private :: trace_spherical_radial_q_face_pair_admitted
  private :: trace_spherical_radial_endpoint_matrix
  private :: trace_spherical_radial_project_to_surface
  private :: trace_spherical_qperp_prepare_seed_result
  private :: trace_spherical_qperp_finalize_from_states
  private :: trace_spherical_qperp_compute_from_states
  private :: trace_spherical_q_from_endpoint_matrices
  private :: trace_spherical_basis,trace_spherical_coord_to_cart
  private :: trace_cart_to_spherical_coord
  private :: trace_spherical_sample_B_bhat_cart
  private :: trace_spherical_sample_bhat_gradbhat_covariant
  private :: trace_spherical_sample_bhat_gradbhat_cartfd
  private :: trace_spherical_bhat_cart_to_rhs
  private :: trace_spherical_endpoint_B_bhat
  private :: trace_spherical_face_probe_coord
  private :: trace_spherical_physical_cell_scale
  private :: trace_spherical_cell_scale_from_widths
  private :: trace_spherical_interp_ctx_build
  private :: trace_spherical_interp_ctx_build_cached
  private :: trace_spherical_effective_step_ctx
  private :: trace_spherical_sample_Bsph_ctx
  private :: trace_spherical_sample_B_bhat_cart_ctx
  private :: trace_spherical_sample_cached_curlB_ctx
  private :: trace_spherical_sample_bhat_gradbhat_covariant_ctx
  private :: get_K_spherical_ctx
  private :: trace_spherical_interp_ctx_fill
  private :: trace_spherical_interp_ctx_load_Bcorners
  private :: trace_interp_weights_block,trace_interp_index_3d
  private :: trace_interp_weights_block_near,trace_interp_index_3d_near
  private :: trace_spherical_local_cell_widths
  private :: trace_spherical_profile_add_count
  private :: trace_spherical_profile_add_time
  private :: trace_spherical_profile_note_trace_steps
  private :: trace_spherical_profile_time
  private :: trace_cartesian_like_geometry,trace_rk45_position_integrator
  private :: trace_cartesian_rhs_bhat
  private :: trace_cartesian_local_cell_size
  private :: trace_rk45_try_boundary_finish
  private :: trace_rk2_stats_note_direction
  private :: trace_rk2_stats_note_step
  private :: trace_rk2_stats_note_completion
  private :: trace_rk2_stats_note_step_limit
  private :: trace_rk2_stats_note_group_iteration
  private :: trace_rk2_stats_note_grid_group
  private :: trace_rk2_stats_note_rhs
  private :: trace_rk2_stats_add_grouped_time
  private :: trace_rk45_stats_note_direction
  private :: trace_rk45_stats_note_attempt
  private :: trace_rk45_stats_note_tangent_error

  type, private :: trace_sph_interp_ctx
    logical :: valid=.false.
    integer :: igrid=-1
    integer :: ixbl(3)=0
    double precision :: x(3)=zero
    double precision :: xd(3)=zero
    double precision :: dxc(3)=zero
    double precision :: w(0:1,3)=zero
    double precision :: dxloc(3)=zero
    double precision :: r=zero,theta=zero,phi=zero
    double precision :: sin_theta=zero,cos_theta=zero
    double precision :: sin_phi=zero,cos_phi=zero
    double precision :: rsin_theta=zero
    double precision :: h_local=zero
    logical :: bcorner_valid=.false.
    double precision :: bcorner(0:1,0:1,0:1,3)=zero
  end type trace_sph_interp_ctx

  type, private :: trace_spherical_curl_cache_block
    double precision, dimension(:^D&,:), allocatable :: current
    logical :: ready=.false.
  end type trace_spherical_curl_cache_block

  type(trace_spherical_curl_cache_block), allocatable, private :: &
       trace_spherical_curl_cache(:)
  logical, private :: trace_spherical_curl_cache_ready=.false.

  integer, parameter, private :: trace_step_control_fixed=0
  integer, parameter, private :: trace_step_control_cell_fraction=1
  integer, private :: trace_step_control_mode=trace_step_control_fixed
  double precision, private :: trace_step_fraction=0.1d0
  double precision, private :: trace_step_min=0.d0

  integer, parameter, private :: trace_integrator_rk2=0
  integer, parameter, private :: trace_integrator_rk45_cartesian=1
  integer, parameter, private :: trace_integrator_rk45_spherical=2
  integer, parameter, private :: trace_tangent_group_rk2_general=1
  integer, parameter, private :: trace_tangent_group_rk2_short_boundary=2
  integer, parameter, private :: trace_tangent_group_rk45_cartesian=3
  integer, parameter, private :: trace_rk45_reject_error=1
  integer, parameter, private :: trace_rk45_reject_stage_failure=2
  integer, parameter, private :: trace_rk45_reject_stage_outside=3
  integer, parameter, private :: trace_rk45_reject_boundary=4
  integer, parameter, private :: trace_rk45_reject_tangent_error=5
  integer, private :: trace_integrator_mode=trace_integrator_rk2
  double precision, private :: trace_rk45_atol=1.d-8
  double precision, private :: trace_rk45_rtol=1.d-6
  double precision, private :: trace_rk45_safety=0.9d0
  double precision, private :: trace_rk45_min_shrink=0.2d0
  double precision, private :: trace_rk45_max_grow=5.d0
  double precision, private :: trace_rk45_tangent_floor=1.d0
  double precision, private :: trace_rk45_tangent_rtol=2.d-5
  logical, private :: trace_rk2_stats_enabled=.false.
  integer(kind=8), private :: trace_rk2_directions=0_8
  integer(kind=8), private :: trace_rk2_completed_directions=0_8
  integer(kind=8), private :: trace_rk2_total_steps=0_8
  integer(kind=8), private :: trace_rk2_boundary_completions=0_8
  integer(kind=8), private :: trace_rk2_max_steps_failures=0_8
  integer(kind=8), private :: trace_rk2_domain_limited_steps=0_8
  integer(kind=8), private :: trace_rk2_grid_transition_steps=0_8
  integer(kind=8), private :: trace_rk2_cell_fraction_limited=0_8
  integer(kind=8), private :: trace_rk2_dL_limited=0_8
  integer(kind=8), private :: trace_rk2_dL_min_limited=0_8
  integer(kind=8), private :: trace_rk2_group_iterations=0_8
  integer(kind=8), private :: trace_rk2_grid_groups=0_8
  integer(kind=8), private :: trace_rk2_group_state_sum=0_8
  integer(kind=8), private :: trace_rk2_group_state_max=0_8
  integer(kind=8), private :: trace_rk2_rhs_calls=0_8
  double precision, private :: trace_rk2_sum_step_length=0.d0
  double precision, private :: trace_rk2_min_step_length=huge(one)
  double precision, private :: trace_rk2_max_step_length=0.d0
  double precision, private :: trace_rk2_grouped_time=0.d0
  integer(kind=8), private :: trace_rk2_max_steps_per_direction=0_8
  logical, private :: trace_rk45_stats_enabled=.false.
  integer(kind=8), private :: trace_rk45_directions=0_8
  integer(kind=8), private :: trace_rk45_attempts=0_8
  integer(kind=8), private :: trace_rk45_accepted=0_8
  integer(kind=8), private :: trace_rk45_rejected=0_8
  integer(kind=8), private :: trace_rk45_rejected_error=0_8
  integer(kind=8), private :: trace_rk45_rejected_position_error=0_8
  integer(kind=8), private :: trace_rk45_rejected_tangent_error=0_8
  integer(kind=8), private :: trace_rk45_rejected_stage_failure=0_8
  integer(kind=8), private :: trace_rk45_rejected_stage_outside=0_8
  integer(kind=8), private :: trace_rk45_rejected_boundary=0_8
  integer(kind=8), private :: trace_rk45_boundary_limited=0_8
  double precision, private :: trace_rk45_sum_accepted_h=0.d0
  double precision, private :: trace_rk45_min_accepted_h=huge(one)
  double precision, private :: trace_rk45_max_accepted_h=0.d0
  integer(kind=8), private :: trace_rk45_tangent_error_samples=0_8
  integer(kind=8), private :: trace_rk45_tangent_would_reject=0_8
  double precision, private :: trace_rk45_sum_tangent_error_ratio=0.d0
  double precision, private :: trace_rk45_max_tangent_error_ratio=0.d0
  double precision, private :: trace_rk45_max_tangent_u_error_ratio=0.d0
  double precision, private :: trace_rk45_max_tangent_v_error_ratio=0.d0
  double precision, private :: trace_rk45_max_position_error_ratio=0.d0

  logical, private :: trace_spherical_profile_enabled=.false.
  logical, parameter, private :: trace_spherical_profile_timers=.false.
  integer(kind=8), private :: trace_profile_seeds=0_8
  integer(kind=8), private :: trace_profile_directions=0_8
  integer(kind=8), private :: trace_profile_steps=0_8
  integer(kind=8), private :: trace_profile_boundary_events=0_8
  integer(kind=8), private :: trace_profile_endpoint_finalizations=0_8
  integer(kind=8), private :: trace_profile_context_requests=0_8
  integer(kind=8), private :: trace_profile_same_cell_hits=0_8
  integer(kind=8), private :: trace_profile_same_grid_hits=0_8
  integer(kind=8), private :: trace_profile_full_context_builds=0_8
  integer(kind=8), private :: trace_profile_context_failures=0_8
  integer(kind=8), private :: trace_profile_cache_invalidations=0_8
  integer(kind=8), private :: trace_profile_b_samples=0_8
  integer(kind=8), private :: trace_profile_curl_samples=0_8
  integer(kind=8), private :: trace_profile_grad_samples=0_8
  integer(kind=8), private :: trace_profile_grad_success=0_8
  integer(kind=8), private :: trace_profile_grad_fallbacks=0_8
  integer(kind=8), private :: trace_profile_b_corner_loads=0_8
  integer(kind=8), private :: trace_profile_b_corner_hits=0_8
  integer(kind=8), private :: trace_profile_hlocal_evals=0_8
  integer(kind=8), private :: trace_profile_face_probes=0_8
  integer(kind=8), private :: trace_profile_tangent_systems=0_8
  integer(kind=8), private :: trace_profile_q_tangent_systems=0_8
  integer(kind=8), private :: trace_profile_qperp_tangent_systems=0_8
  integer(kind=8), private :: trace_profile_combined_tangent_systems=0_8
  integer(kind=8), private :: trace_profile_max_steps_per_trace=0_8
  double precision, private :: trace_profile_time_context=0.d0
  double precision, private :: trace_profile_time_b_interp=0.d0
  double precision, private :: trace_profile_time_curl_interp=0.d0
  double precision, private :: trace_profile_time_grad_bhat=0.d0
  double precision, private :: trace_profile_time_step=0.d0
  double precision, private :: trace_profile_time_tangent=0.d0

  type, private :: trace_summary_state
    double precision :: x(ndim)
    double precision :: footpoint(ndim)
    double precision :: length
    double precision :: twist
    integer :: nstep
    integer :: status
    integer :: twist_status
    integer :: igrid
    integer :: seed_id
    integer :: face
    logical :: forward
    logical :: active
    logical :: accumulate_twist
    double precision :: rk45_h
    type(trace_sph_interp_ctx) :: sph_cache
  end type trace_summary_state

  type, private :: trace_tangent_state
    integer :: seed_id
    integer :: direction
    logical :: active
    logical :: complete
    logical :: has_extra
    double precision :: x(ndim)
    double precision :: u(ndim)
    double precision :: v(ndim)
    double precision :: p(ndim)
    double precision :: q(ndim)
    double precision :: length
    double precision :: twist
    integer :: nstep
    integer :: igrid
    integer :: status
    integer :: twist_status
    integer :: face_id
    double precision :: endpoint(ndim)
    double precision :: endpoint_B(ndim)
    double precision :: endpoint_bhat(ndim)
    double precision :: u_perp(ndim)
    double precision :: v_perp(ndim)
    double precision :: p_perp(ndim)
    double precision :: q_perp(ndim)
    double precision :: rk45_h
    logical :: accumulate_twist
    type(trace_sph_interp_ctx) :: sph_cache
  end type trace_tangent_state

  type, public :: trace_length_result
    double precision :: seed(ndim)
    double precision :: forward_footpoint(ndim)
    double precision :: backward_footpoint(ndim)
    double precision :: forward_length
    double precision :: backward_length
    double precision :: total_length
    integer :: forward_nstep
    integer :: backward_nstep
    integer :: forward_status
    integer :: backward_status
  end type trace_length_result

  type, public :: trace_twist_result
    type(trace_length_result) :: line
    double precision :: forward_twist
    double precision :: backward_twist
    double precision :: total_twist
    logical :: valid_twist
    integer :: status_twist
  end type trace_twist_result

  type, public :: trace_mapping_result
    double precision :: seed(ndim)
    double precision :: source_B(3)
    double precision :: forward_footpoint(ndim)
    double precision :: backward_footpoint(ndim)
    double precision :: forward_B(3)
    double precision :: backward_B(3)
    double precision :: forward_length
    double precision :: backward_length
    double precision :: source_Bn
    double precision :: forward_Bn
    double precision :: backward_Bn
    integer :: forward_face
    integer :: backward_face
    integer :: forward_status
    integer :: backward_status
    logical :: valid
  end type trace_mapping_result

  type, public :: trace_qperp_result
    double precision :: seed(ndim)
    double precision :: qperp
    double precision :: logqperp
    double precision :: N2
    double precision :: bfactor
    double precision :: q0
    double precision :: logq0
    double precision :: qperp0
    double precision :: logqperp0
    double precision :: N2_qperp0
    double precision :: bfactor_qperp0
    double precision :: forward_Bn_q0
    double precision :: backward_Bn_q0
    double precision :: B_seed(ndim)
    double precision :: bhat_seed(ndim)
    double precision :: forward_endpoint(ndim)
    double precision :: backward_endpoint(ndim)
    double precision :: forward_B(ndim)
    double precision :: backward_B(ndim)
    double precision :: forward_bhat(ndim)
    double precision :: backward_bhat(ndim)
    double precision :: forward_length
    double precision :: backward_length
    integer :: forward_nstep
    integer :: backward_nstep
    double precision :: u0(ndim)
    double precision :: v0(ndim)
    double precision :: u_forward_perp(ndim)
    double precision :: v_forward_perp(ndim)
    double precision :: u_backward_perp(ndim)
    double precision :: v_backward_perp(ndim)
    integer :: forward_face
    integer :: backward_face
    integer :: forward_status
    integer :: backward_status
    integer :: status_q0
    integer :: status_qperp0
    integer :: status
    logical :: valid
    logical :: valid_q0
    logical :: valid_qperp0
  end type trace_qperp_result

  type, public :: trace_topology_result
    double precision :: seed(ndim)
    double precision :: length_forward
    double precision :: length_backward
    double precision :: length_total
    double precision :: twist_forward
    double precision :: twist_backward
    double precision :: twist_total
    double precision :: forward_endpoint(ndim)
    double precision :: backward_endpoint(ndim)
    integer :: forward_nstep
    integer :: backward_nstep
    integer :: forward_face
    integer :: backward_face
    integer :: forward_status
    integer :: backward_status
    double precision :: map_forward_endpoint(ndim)
    double precision :: map_backward_endpoint(ndim)
    double precision :: map_forward_length
    double precision :: map_backward_length
    integer :: map_forward_face
    integer :: map_backward_face
    integer :: map_forward_status
    integer :: map_backward_status
    double precision :: source_B(3)
    double precision :: forward_B(3)
    double precision :: backward_B(3)
    double precision :: source_Bn
    double precision :: forward_Bn
    double precision :: backward_Bn
    logical :: has_twist
    logical :: has_mapping
    logical :: valid_twist
    logical :: valid
    integer :: status_twist
    integer :: status
  end type trace_topology_result

contains

  subroutine trace_set_step_control(mode,step_fraction,step_min)
    character(len=*), intent(in) :: mode
    double precision, intent(in) :: step_fraction,step_min

    character(len=len(mode)) :: mode_lc
    integer :: i,code

    mode_lc=mode
    do i=1,len(mode_lc)
      code=iachar(mode_lc(i:i))
      if (code>=iachar('A') .and. code<=iachar('Z')) then
        mode_lc(i:i)=achar(code+iachar('a')-iachar('A'))
      endif
    enddo

    select case (trim(mode_lc))
    case ('cell_fraction')
      trace_step_control_mode=trace_step_control_cell_fraction
    case default
      trace_step_control_mode=trace_step_control_fixed
    end select
    trace_step_fraction=step_fraction
    trace_step_min=max(zero,step_min)
  end subroutine trace_set_step_control

  subroutine trace_set_integrator(mode,atol,rtol,safety,min_shrink, &
       max_grow,tangent_floor,tangent_rtol)
    character(len=*), intent(in) :: mode
    double precision, intent(in) :: atol,rtol,safety,min_shrink,max_grow
    double precision, intent(in) :: tangent_floor,tangent_rtol

    character(len=len(mode)) :: mode_lc
    integer :: i,code

    mode_lc=mode
    do i=1,len(mode_lc)
      code=iachar(mode_lc(i:i))
      if (code>=iachar('A') .and. code<=iachar('Z')) then
        mode_lc(i:i)=achar(code+iachar('a')-iachar('A'))
      endif
    enddo

    select case (trim(mode_lc))
    case ('rk45_cartesian')
      trace_integrator_mode=trace_integrator_rk45_cartesian
    case ('rk45_spherical')
      trace_integrator_mode=trace_integrator_rk45_spherical
    case default
      trace_integrator_mode=trace_integrator_rk2
    end select
    trace_rk45_atol=max(atol,zero)
    trace_rk45_rtol=max(rtol,zero)
    trace_rk45_safety=max(0.1d0,min(safety,one))
    trace_rk45_min_shrink=max(0.01d0,min(min_shrink,one))
    trace_rk45_max_grow=max(one,max_grow)
    trace_rk45_tangent_floor=max(tangent_floor,smalldouble)
    trace_rk45_tangent_rtol=max(tangent_rtol,smalldouble)
  end subroutine trace_set_integrator

  subroutine trace_rk2_stats_set_enabled(enabled)
    logical, intent(in) :: enabled

    trace_rk2_stats_enabled=enabled
  end subroutine trace_rk2_stats_set_enabled

  subroutine trace_rk2_stats_reset()
    trace_rk2_directions=0_8
    trace_rk2_completed_directions=0_8
    trace_rk2_total_steps=0_8
    trace_rk2_boundary_completions=0_8
    trace_rk2_max_steps_failures=0_8
    trace_rk2_domain_limited_steps=0_8
    trace_rk2_grid_transition_steps=0_8
    trace_rk2_cell_fraction_limited=0_8
    trace_rk2_dL_limited=0_8
    trace_rk2_dL_min_limited=0_8
    trace_rk2_group_iterations=0_8
    trace_rk2_grid_groups=0_8
    trace_rk2_group_state_sum=0_8
    trace_rk2_group_state_max=0_8
    trace_rk2_rhs_calls=0_8
    trace_rk2_sum_step_length=0.d0
    trace_rk2_min_step_length=huge(one)
    trace_rk2_max_step_length=0.d0
    trace_rk2_grouped_time=0.d0
    trace_rk2_max_steps_per_direction=0_8
  end subroutine trace_rk2_stats_reset

  subroutine trace_rk2_stats_report(label)
    character(len=*), intent(in) :: label

    double precision :: mean_group_size,mean_h,mean_steps_per_direction

    if (.not.trace_rk2_stats_enabled) return

    mean_h=zero
    if (trace_rk2_total_steps>0_8) then
      mean_h=trace_rk2_sum_step_length/dble(trace_rk2_total_steps)
    endif
    mean_steps_per_direction=zero
    if (trace_rk2_completed_directions>0_8) then
      mean_steps_per_direction=dble(trace_rk2_total_steps)/ &
           dble(trace_rk2_completed_directions)
    endif
    mean_group_size=zero
    if (trace_rk2_grid_groups>0_8) then
      mean_group_size=dble(trace_rk2_group_state_sum)/ &
           dble(trace_rk2_grid_groups)
    endif

    write(*,'(a)') '[mt_rk2] Cartesian field-line RK2 summary: '// &
         trim(label)
    write(*,'(a,i0)') '  traced_directions: ',trace_rk2_directions
    write(*,'(a,i0)') '  completed_directions: ', &
         trace_rk2_completed_directions
    write(*,'(a,i0)') '  total_steps: ',trace_rk2_total_steps
    write(*,'(a,es12.4)') '  mean_step_length: ',mean_h
    if (trace_rk2_total_steps>0_8) then
      write(*,'(a,es12.4)') '  min_step_length: ', &
           trace_rk2_min_step_length
      write(*,'(a,es12.4)') '  max_step_length: ', &
           trace_rk2_max_step_length
    endif
    write(*,'(a,es12.4)') '  mean_steps_per_direction: ', &
         mean_steps_per_direction
    write(*,'(a,i0)') '  max_steps_per_direction: ', &
         trace_rk2_max_steps_per_direction
    write(*,'(a,i0)') '  boundary_completions: ', &
         trace_rk2_boundary_completions
    write(*,'(a,i0)') '  max_steps_failures: ', &
         trace_rk2_max_steps_failures
    write(*,'(a,i0)') '  domain_boundary_limited_steps: ', &
         trace_rk2_domain_limited_steps
    write(*,'(a,i0)') '  block_grid_clipped_steps: 0'
    write(*,'(a,i0)') '  grid_transition_steps: ', &
         trace_rk2_grid_transition_steps
    write(*,'(a,i0)') '  cell_fraction_limited_evals: ', &
         trace_rk2_cell_fraction_limited
    write(*,'(a,i0)') '  dL_cap_limited_evals: ',trace_rk2_dL_limited
    write(*,'(a,i0)') '  dL_min_limited_evals: ',trace_rk2_dL_min_limited
    write(*,'(a,i0)') '  grouped_driver_outer_iterations: ', &
         trace_rk2_group_iterations
    write(*,'(a,i0)') '  grouped_driver_parallel_batches: ', &
         trace_rk2_grid_groups
    write(*,'(a,es12.4)') '  grouped_driver_mean_batch_size: ', &
         mean_group_size
    write(*,'(a,i0)') '  grouped_driver_max_batch_size: ', &
         trace_rk2_group_state_max
    write(*,'(a,es12.4)') '  grouped_driver_cpu_time: ', &
         trace_rk2_grouped_time
    write(*,'(a,i0)') '  trace_tangent_rhs_calls: ',trace_rk2_rhs_calls
    write(*,'(a)') '  note: Cartesian RK2 grouped tracing does not clip '// &
         'steps to AMR block/grid boundaries; grid changes are counted after '// &
         'accepted full steps.'
  end subroutine trace_rk2_stats_report

  subroutine trace_rk2_stats_note_direction()
    if (.not.trace_rk2_stats_enabled) return
    !$OMP CRITICAL(trace_rk2_stats)
    trace_rk2_directions=trace_rk2_directions+1_8
    !$OMP END CRITICAL(trace_rk2_stats)
  end subroutine trace_rk2_stats_note_direction

  subroutine trace_rk2_stats_note_step(h,domain_limited,grid_transition)
    double precision, intent(in) :: h
    logical, intent(in) :: domain_limited,grid_transition

    if (.not.trace_rk2_stats_enabled) return
    !$OMP CRITICAL(trace_rk2_stats)
    trace_rk2_total_steps=trace_rk2_total_steps+1_8
    trace_rk2_sum_step_length=trace_rk2_sum_step_length+abs(h)
    trace_rk2_min_step_length=min(trace_rk2_min_step_length,abs(h))
    trace_rk2_max_step_length=max(trace_rk2_max_step_length,abs(h))
    if (domain_limited) trace_rk2_domain_limited_steps= &
         trace_rk2_domain_limited_steps+1_8
    if (grid_transition) trace_rk2_grid_transition_steps= &
         trace_rk2_grid_transition_steps+1_8
    !$OMP END CRITICAL(trace_rk2_stats)
  end subroutine trace_rk2_stats_note_step

  subroutine trace_rk2_stats_note_completion(status,nstep)
    integer, intent(in) :: status,nstep

    if (.not.trace_rk2_stats_enabled) return
    !$OMP CRITICAL(trace_rk2_stats)
    trace_rk2_completed_directions=trace_rk2_completed_directions+1_8
    trace_rk2_max_steps_per_direction=max( &
         trace_rk2_max_steps_per_direction,int(nstep,kind=8))
    select case (status)
    case (trace_status_boundary)
      trace_rk2_boundary_completions=trace_rk2_boundary_completions+1_8
    case (trace_status_max_steps)
      trace_rk2_max_steps_failures=trace_rk2_max_steps_failures+1_8
    end select
    !$OMP END CRITICAL(trace_rk2_stats)
  end subroutine trace_rk2_stats_note_completion

  subroutine trace_rk2_stats_note_step_limit(step_cap,cell_step,actual_step)
    double precision, intent(in) :: step_cap,cell_step,actual_step

    double precision :: tol

    if (.not.trace_rk2_stats_enabled) return
    tol=100.d0*epsilon(one)*max(one,step_cap,cell_step,actual_step)
    !$OMP CRITICAL(trace_rk2_stats)
    if (actual_step>=step_cap-tol) then
      trace_rk2_dL_limited=trace_rk2_dL_limited+1_8
    else if (trace_step_min>zero .and. actual_step>=trace_step_min-tol .and. &
         actual_step>cell_step+tol) then
      trace_rk2_dL_min_limited=trace_rk2_dL_min_limited+1_8
    else
      trace_rk2_cell_fraction_limited=trace_rk2_cell_fraction_limited+1_8
    endif
    !$OMP END CRITICAL(trace_rk2_stats)
  end subroutine trace_rk2_stats_note_step_limit

  subroutine trace_rk2_stats_note_group_iteration()
    if (.not.trace_rk2_stats_enabled) return
    !$OMP CRITICAL(trace_rk2_stats)
    trace_rk2_group_iterations=trace_rk2_group_iterations+1_8
    !$OMP END CRITICAL(trace_rk2_stats)
  end subroutine trace_rk2_stats_note_group_iteration

  subroutine trace_rk2_stats_note_grid_group(group_size)
    integer, intent(in) :: group_size

    if (.not.trace_rk2_stats_enabled) return
    !$OMP CRITICAL(trace_rk2_stats)
    trace_rk2_grid_groups=trace_rk2_grid_groups+1_8
    trace_rk2_group_state_sum=trace_rk2_group_state_sum+ &
         int(group_size,kind=8)
    trace_rk2_group_state_max=max(trace_rk2_group_state_max, &
         int(group_size,kind=8))
    !$OMP END CRITICAL(trace_rk2_stats)
  end subroutine trace_rk2_stats_note_grid_group

  subroutine trace_rk2_stats_note_rhs()
    if (.not.trace_rk2_stats_enabled) return
    !$OMP ATOMIC
    trace_rk2_rhs_calls=trace_rk2_rhs_calls+1_8
  end subroutine trace_rk2_stats_note_rhs

  subroutine trace_rk2_stats_add_grouped_time(dt)
    double precision, intent(in) :: dt

    if (.not.trace_rk2_stats_enabled) return
    !$OMP CRITICAL(trace_rk2_stats)
    trace_rk2_grouped_time=trace_rk2_grouped_time+dt
    !$OMP END CRITICAL(trace_rk2_stats)
  end subroutine trace_rk2_stats_add_grouped_time

  subroutine trace_rk45_stats_set_enabled(enabled)
    logical, intent(in) :: enabled

    trace_rk45_stats_enabled=enabled
  end subroutine trace_rk45_stats_set_enabled

  subroutine trace_rk45_stats_reset()
    trace_rk45_directions=0_8
    trace_rk45_attempts=0_8
    trace_rk45_accepted=0_8
    trace_rk45_rejected=0_8
    trace_rk45_rejected_error=0_8
    trace_rk45_rejected_position_error=0_8
    trace_rk45_rejected_tangent_error=0_8
    trace_rk45_rejected_stage_failure=0_8
    trace_rk45_rejected_stage_outside=0_8
    trace_rk45_rejected_boundary=0_8
    trace_rk45_boundary_limited=0_8
    trace_rk45_sum_accepted_h=0.d0
    trace_rk45_min_accepted_h=huge(one)
    trace_rk45_max_accepted_h=0.d0
    trace_rk45_tangent_error_samples=0_8
    trace_rk45_tangent_would_reject=0_8
    trace_rk45_sum_tangent_error_ratio=0.d0
    trace_rk45_max_tangent_error_ratio=0.d0
    trace_rk45_max_tangent_u_error_ratio=0.d0
    trace_rk45_max_tangent_v_error_ratio=0.d0
    trace_rk45_max_position_error_ratio=0.d0
  end subroutine trace_rk45_stats_reset

  subroutine trace_rk45_stats_report(label)
    character(len=*), intent(in) :: label

    double precision :: mean_h,mean_reject_per_trace,mean_tangent_error

    if (.not.trace_rk45_stats_enabled) return
    if (trace_integrator_mode/=trace_integrator_rk45_cartesian .and. &
         trace_integrator_mode/=trace_integrator_rk45_spherical) return

    mean_h=zero
    if (trace_rk45_accepted>0_8) then
      mean_h=trace_rk45_sum_accepted_h/dble(trace_rk45_accepted)
    endif
    mean_reject_per_trace=zero
    if (trace_rk45_directions>0_8) then
      mean_reject_per_trace=dble(trace_rk45_rejected)/dble(trace_rk45_directions)
    endif

    select case (trace_integrator_mode)
    case (trace_integrator_rk45_cartesian)
      write(*,'(a)') '[mt_rk45] Cartesian field-line RK45 summary: '// &
           trim(label)
    case (trace_integrator_rk45_spherical)
      write(*,'(a)') '[mt_rk45] Spherical field-line RK45 summary: '// &
           trim(label)
    end select
    write(*,'(a,i0)') '  traced_directions: ',trace_rk45_directions
    write(*,'(a,i0)') '  attempts: ',trace_rk45_attempts
    write(*,'(a,i0)') '  accepted_steps: ',trace_rk45_accepted
    write(*,'(a,i0)') '  rejected_steps: ',trace_rk45_rejected
    write(*,'(a,i0)') '  rejected_error_steps: ', &
         trace_rk45_rejected_error
    write(*,'(a,i0)') '  rejected_position_error_steps: ', &
         trace_rk45_rejected_position_error
    write(*,'(a,i0)') '  rejected_tangent_error_steps: ', &
         trace_rk45_rejected_tangent_error
    write(*,'(a,i0)') '  rejected_stage_failure_steps: ', &
         trace_rk45_rejected_stage_failure
    write(*,'(a,i0)') '  rejected_stage_outside_steps: ', &
         trace_rk45_rejected_stage_outside
    write(*,'(a,i0)') '  rejected_boundary_steps: ', &
         trace_rk45_rejected_boundary
    write(*,'(a,es12.4)') '  rejected_steps_per_direction: ', &
         mean_reject_per_trace
    write(*,'(a,i0)') '  boundary_limited_steps: ', &
         trace_rk45_boundary_limited
    write(*,'(a,es12.4)') '  mean_accepted_step: ',mean_h
    if (trace_rk45_accepted>0_8) then
      write(*,'(a,es12.4)') '  min_accepted_step: ', &
           trace_rk45_min_accepted_h
      write(*,'(a,es12.4)') '  max_accepted_step: ', &
           trace_rk45_max_accepted_h
    endif
    if (trace_rk45_tangent_error_samples>0_8) then
      write(*,'(a)') '  error_control_mode: position_tangent'
      mean_tangent_error=trace_rk45_sum_tangent_error_ratio/ &
           dble(trace_rk45_tangent_error_samples)
      write(*,'(a,i0)') '  tangent_error_samples: ', &
           trace_rk45_tangent_error_samples
      write(*,'(a,i0)') '  tangent_would_reject_steps: ', &
           trace_rk45_tangent_would_reject
      write(*,'(a,es12.4)') '  max_position_error_ratio: ', &
           trace_rk45_max_position_error_ratio
      write(*,'(a,es12.4)') '  mean_tangent_error_ratio: ', &
           mean_tangent_error
      write(*,'(a,es12.4)') '  max_tangent_error_ratio: ', &
           trace_rk45_max_tangent_error_ratio
      write(*,'(a,es12.4)') '  max_tangent_u_error_ratio: ', &
           trace_rk45_max_tangent_u_error_ratio
      write(*,'(a,es12.4)') '  max_tangent_v_error_ratio: ', &
           trace_rk45_max_tangent_v_error_ratio
      write(*,'(a,es12.4)') '  tangent_error_rtol: ', &
           trace_rk45_tangent_rtol
      write(*,'(a,es12.4)') '  tangent_error_floor: ', &
           trace_rk45_tangent_floor
      write(*,'(a)') '  tangent_error_note: standard-logQ tangent '// &
           'traces use position+tangent RK45 accept/reject.'
    else
      write(*,'(a)') '  error_control_mode: position'
    endif
    write(*,'(a)') '  note: this experimental RK45 mode traces '// &
         'position/length/twist and tangent-transport Q products where enabled.'
  end subroutine trace_rk45_stats_report

  subroutine trace_rk45_stats_note_direction()
    if (.not.trace_rk45_stats_enabled) return
    if (trace_integrator_mode/=trace_integrator_rk45_cartesian .and. &
         trace_integrator_mode/=trace_integrator_rk45_spherical) return
    !$OMP CRITICAL(trace_rk45_stats)
    trace_rk45_directions=trace_rk45_directions+1_8
    !$OMP END CRITICAL(trace_rk45_stats)
  end subroutine trace_rk45_stats_note_direction

  subroutine trace_rk45_stats_note_attempt(accepted,boundary_limited,h, &
       reject_reason)
    logical, intent(in) :: accepted,boundary_limited
    double precision, intent(in) :: h
    integer, intent(in), optional :: reject_reason

    if (.not.trace_rk45_stats_enabled) return
    if (trace_integrator_mode/=trace_integrator_rk45_cartesian .and. &
         trace_integrator_mode/=trace_integrator_rk45_spherical) return
    !$OMP CRITICAL(trace_rk45_stats)
    trace_rk45_attempts=trace_rk45_attempts+1_8
    if (accepted) then
      trace_rk45_accepted=trace_rk45_accepted+1_8
      trace_rk45_sum_accepted_h=trace_rk45_sum_accepted_h+abs(h)
      trace_rk45_min_accepted_h=min(trace_rk45_min_accepted_h,abs(h))
      trace_rk45_max_accepted_h=max(trace_rk45_max_accepted_h,abs(h))
      if (boundary_limited) trace_rk45_boundary_limited= &
           trace_rk45_boundary_limited+1_8
    else
      trace_rk45_rejected=trace_rk45_rejected+1_8
      if (present(reject_reason)) then
        select case (reject_reason)
        case (trace_rk45_reject_error)
          trace_rk45_rejected_error=trace_rk45_rejected_error+1_8
          trace_rk45_rejected_position_error= &
               trace_rk45_rejected_position_error+1_8
        case (trace_rk45_reject_tangent_error)
          trace_rk45_rejected_error=trace_rk45_rejected_error+1_8
          trace_rk45_rejected_tangent_error= &
               trace_rk45_rejected_tangent_error+1_8
        case (trace_rk45_reject_stage_failure)
          trace_rk45_rejected_stage_failure= &
               trace_rk45_rejected_stage_failure+1_8
        case (trace_rk45_reject_stage_outside)
          trace_rk45_rejected_stage_outside= &
               trace_rk45_rejected_stage_outside+1_8
        case (trace_rk45_reject_boundary)
          trace_rk45_rejected_boundary=trace_rk45_rejected_boundary+1_8
        end select
      endif
    endif
    !$OMP END CRITICAL(trace_rk45_stats)
  end subroutine trace_rk45_stats_note_attempt

  subroutine trace_rk45_stats_note_tangent_error(position_error_ratio, &
       u_error_ratio,v_error_ratio)
    double precision, intent(in) :: position_error_ratio
    double precision, intent(in) :: u_error_ratio,v_error_ratio

    double precision :: tangent_error_ratio

    if (.not.trace_rk45_stats_enabled) return
    if (trace_integrator_mode/=trace_integrator_rk45_cartesian .and. &
         trace_integrator_mode/=trace_integrator_rk45_spherical) return
    tangent_error_ratio=max(u_error_ratio,v_error_ratio)
    !$OMP CRITICAL(trace_rk45_stats)
    trace_rk45_tangent_error_samples= &
         trace_rk45_tangent_error_samples+1_8
    trace_rk45_sum_tangent_error_ratio= &
         trace_rk45_sum_tangent_error_ratio+tangent_error_ratio
    trace_rk45_max_tangent_error_ratio= &
         max(trace_rk45_max_tangent_error_ratio,tangent_error_ratio)
    trace_rk45_max_tangent_u_error_ratio= &
         max(trace_rk45_max_tangent_u_error_ratio,u_error_ratio)
    trace_rk45_max_tangent_v_error_ratio= &
         max(trace_rk45_max_tangent_v_error_ratio,v_error_ratio)
    trace_rk45_max_position_error_ratio= &
         max(trace_rk45_max_position_error_ratio,position_error_ratio)
    if (tangent_error_ratio>trace_rk45_tangent_rtol) &
         trace_rk45_tangent_would_reject= &
         trace_rk45_tangent_would_reject+1_8
    !$OMP END CRITICAL(trace_rk45_stats)
  end subroutine trace_rk45_stats_note_tangent_error

  subroutine trace_spherical_profile_set(enabled)
    logical, intent(in) :: enabled

    trace_spherical_profile_enabled=enabled
  end subroutine trace_spherical_profile_set

  subroutine trace_spherical_profile_reset()
    trace_profile_seeds=0_8
    trace_profile_directions=0_8
    trace_profile_steps=0_8
    trace_profile_boundary_events=0_8
    trace_profile_endpoint_finalizations=0_8
    trace_profile_context_requests=0_8
    trace_profile_same_cell_hits=0_8
    trace_profile_same_grid_hits=0_8
    trace_profile_full_context_builds=0_8
    trace_profile_context_failures=0_8
    trace_profile_cache_invalidations=0_8
    trace_profile_b_samples=0_8
    trace_profile_curl_samples=0_8
    trace_profile_grad_samples=0_8
    trace_profile_grad_success=0_8
    trace_profile_grad_fallbacks=0_8
    trace_profile_b_corner_loads=0_8
    trace_profile_b_corner_hits=0_8
    trace_profile_hlocal_evals=0_8
    trace_profile_face_probes=0_8
    trace_profile_tangent_systems=0_8
    trace_profile_q_tangent_systems=0_8
    trace_profile_qperp_tangent_systems=0_8
    trace_profile_combined_tangent_systems=0_8
    trace_profile_max_steps_per_trace=0_8
    trace_profile_time_context=0.d0
    trace_profile_time_b_interp=0.d0
    trace_profile_time_curl_interp=0.d0
    trace_profile_time_grad_bhat=0.d0
    trace_profile_time_step=0.d0
    trace_profile_time_tangent=0.d0
  end subroutine trace_spherical_profile_reset

  subroutine trace_spherical_profile_count_seeds(nseed)
    integer, intent(in) :: nseed

    if (.not.trace_spherical_profile_enabled) return
    call trace_spherical_profile_add_count(trace_profile_seeds, &
         int(max(nseed,0),kind=8))
  end subroutine trace_spherical_profile_count_seeds

  subroutine trace_spherical_profile_report(label)
    character(len=*), intent(in) :: label

    double precision :: avg_steps,denom,hit_denom

    if (.not.trace_spherical_profile_enabled) return
    if (mype/=0) return

    denom=max(one,dble(trace_profile_directions))
    avg_steps=dble(trace_profile_steps)/denom
    hit_denom=max(one,dble(trace_profile_context_requests))

    write(*,'(a)') '[mt_profile] spherical topology profiling summary '// &
         trim(label)
    write(*,'(a,i0)') '  seeds: ',trace_profile_seeds
    write(*,'(a,i0)') '  traced_directions: ',trace_profile_directions
    write(*,'(a,i0)') '  accepted_steps_all_integrators: ', &
         trace_profile_steps
    write(*,'(a,es12.4)') '  avg_steps_per_counted_trace: ',avg_steps
    write(*,'(a,i0)') '  max_steps_per_trace: ', &
         trace_profile_max_steps_per_trace
    write(*,'(a)') '  note: accepted steps include counted spherical '// &
         'summary and tangent integrations; shared Q/Qperp tangent '// &
         'transport is counted once.'
    write(*,'(a)') '  note: profiling counters are diagnostic and may '// &
         'significantly slow hot loops; do not use profiling-on wall time '// &
         'as production timing.'
    write(*,'(a,i0)') '  boundary_events: ',trace_profile_boundary_events
    write(*,'(a,i0)') '  endpoint_finalizations: ', &
         trace_profile_endpoint_finalizations
    write(*,'(a,i0)') '  context_requests: ',trace_profile_context_requests
    write(*,'(a,i0,a,f6.2,a)') '  same_cell_hits: ', &
         trace_profile_same_cell_hits,' (', &
         100.d0*dble(trace_profile_same_cell_hits)/hit_denom,'%)'
    write(*,'(a,i0,a,f6.2,a)') '  same_grid_hits: ', &
         trace_profile_same_grid_hits,' (', &
         100.d0*dble(trace_profile_same_grid_hits)/hit_denom,'%)'
    write(*,'(a,i0,a,f6.2,a)') '  full_context_builds: ', &
         trace_profile_full_context_builds,' (', &
         100.d0*dble(trace_profile_full_context_builds)/hit_denom,'%)'
    write(*,'(a,i0)') '  context_failures: ',trace_profile_context_failures
    write(*,'(a,i0)') '  cache_invalidations: ', &
         trace_profile_cache_invalidations
    write(*,'(a,i0)') '  B_samples: ',trace_profile_b_samples
    write(*,'(a,i0)') '  curl_samples: ',trace_profile_curl_samples
    write(*,'(a,i0)') '  grad_bhat_samples: ',trace_profile_grad_samples
    write(*,'(a,i0)') '  grad_bhat_success: ',trace_profile_grad_success
    write(*,'(a,i0)') '  grad_bhat_fallbacks: ', &
         trace_profile_grad_fallbacks
    write(*,'(a,i0)') '  B_corner_loads: ',trace_profile_b_corner_loads
    write(*,'(a,i0)') '  B_corner_cache_hits: ', &
         trace_profile_b_corner_hits
    write(*,'(a,i0)') '  h_local_evals: ',trace_profile_hlocal_evals
    write(*,'(a,i0)') '  boundary_face_probes: ', &
         trace_profile_face_probes
    write(*,'(a,i0)') '  tangent_systems: ', &
         trace_profile_tangent_systems
    write(*,'(a,i0)') '  q_tangent_systems: ', &
         trace_profile_q_tangent_systems
    write(*,'(a,i0)') '  qperp_tangent_systems: ', &
         trace_profile_qperp_tangent_systems
    write(*,'(a,i0)') '  combined_q_qperp_tangent_systems: ', &
         trace_profile_combined_tangent_systems
    write(*,'(a)') '  timing_cpu_seconds:'
    if (.not.trace_spherical_profile_timers) then
      write(*,'(a)') '    disabled (hot-loop timers distort profiling runs)'
      return
    endif
    write(*,'(a,es12.4)') '    context: ',trace_profile_time_context
    write(*,'(a,es12.4)') '    B_interp: ',trace_profile_time_b_interp
    write(*,'(a,es12.4)') '    curl_interp: ', &
         trace_profile_time_curl_interp
    write(*,'(a,es12.4)') '    grad_bhat: ',trace_profile_time_grad_bhat
    write(*,'(a,es12.4)') '    step_size: ',trace_profile_time_step
    write(*,'(a,es12.4)') '    tangent_transport: ', &
         trace_profile_time_tangent
  end subroutine trace_spherical_profile_report

  subroutine trace_spherical_profile_add_count(counter,delta)
    integer(kind=8), intent(inout) :: counter
    integer(kind=8), intent(in) :: delta

    if (.not.trace_spherical_profile_enabled) return
    !$OMP ATOMIC UPDATE
    counter=counter+delta
  end subroutine trace_spherical_profile_add_count

  subroutine trace_spherical_profile_add_time(timer,dt)
    double precision, intent(inout) :: timer
    double precision, intent(in) :: dt

    if (.not.trace_spherical_profile_enabled) return
    if (.not.trace_spherical_profile_timers) return
    !$OMP ATOMIC UPDATE
    timer=timer+dt
  end subroutine trace_spherical_profile_add_time

  subroutine trace_spherical_profile_note_trace_steps(nstep)
    integer, intent(in) :: nstep

    if (.not.trace_spherical_profile_enabled) return
    !$OMP CRITICAL(trace_spherical_profile_max_steps)
    trace_profile_max_steps_per_trace=max(trace_profile_max_steps_per_trace, &
         int(nstep,kind=8))
    !$OMP END CRITICAL(trace_spherical_profile_max_steps)
  end subroutine trace_spherical_profile_note_trace_steps

  double precision function trace_spherical_profile_time() result(tnow)
    if (.not.trace_spherical_profile_timers) then
      tnow=0.d0
      return
    endif
    call cpu_time(tnow)
  end function trace_spherical_profile_time

  subroutine trace_spherical_curl_cache_build(status)
    integer, intent(out) :: status

    double precision :: bvec(ixG^T,1:3)
    integer :: iigrid,igrid,idirmin
    integer :: ixI^L,ixO^L,ixA^L

    status=trace_status_active
    call trace_spherical_curl_cache_clear()

    if (ndim/=3 .or. geo_coordinate/=geo_spherical) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (.not.B0field .and. .not.allocated(iw_mag)) then
      status=trace_status_out_of_domain
      return
    endif

    {^IFTHREED
    allocate(trace_spherical_curl_cache(max_blocks))
    ixI^L=ixG^LL;
    ixO^L=ixM^LL^LADD1;
    ixA^L=ixO^L^LADD1;
    if (ixImin1>ixAmin1 .or. ixImax1<ixAmax1 .or. &
         ixImin2>ixAmin2 .or. ixImax2<ixAmax2 .or. &
         ixImin3>ixAmin3 .or. ixImax3<ixAmax3) then
      call trace_spherical_curl_cache_clear()
      status=trace_status_bad_curl_stencil
      return
    endif

    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      allocate(trace_spherical_curl_cache(igrid)%current( &
           lbound(ps(igrid)%w,1):ubound(ps(igrid)%w,1), &
           lbound(ps(igrid)%w,2):ubound(ps(igrid)%w,2), &
           lbound(ps(igrid)%w,3):ubound(ps(igrid)%w,3),1:3))
      trace_spherical_curl_cache(igrid)%current=zero
      bvec=zero
      if (allocated(iw_mag)) then
        bvec(:,:,:,1)=ps(igrid)%w(:,:,:,iw_mag(1))
        bvec(:,:,:,2)=ps(igrid)%w(:,:,:,iw_mag(2))
        bvec(:,:,:,3)=ps(igrid)%w(:,:,:,iw_mag(3))
      endif

      block=>ps(igrid)
      call curlvector(bvec,ixI^L,ixO^L, &
           trace_spherical_curl_cache(igrid)%current,idirmin,1,3)
      if (B0field) then
        trace_spherical_curl_cache(igrid)%current(ixO^S,1:3)= &
             trace_spherical_curl_cache(igrid)%current(ixO^S,1:3) &
             +ps(igrid)%J0(ixO^S,1:3)
      endif
      trace_spherical_curl_cache(igrid)%ready=.true.
    enddo
    trace_spherical_curl_cache_ready=.true.
    }
  end subroutine trace_spherical_curl_cache_build

  subroutine trace_spherical_curl_cache_clear()
    integer :: igrid

    if (allocated(trace_spherical_curl_cache)) then
      do igrid=1,size(trace_spherical_curl_cache)
        if (allocated(trace_spherical_curl_cache(igrid)%current)) then
          deallocate(trace_spherical_curl_cache(igrid)%current)
        endif
        trace_spherical_curl_cache(igrid)%ready=.false.
      enddo
      deallocate(trace_spherical_curl_cache)
    endif
    trace_spherical_curl_cache_ready=.false.
  end subroutine trace_spherical_curl_cache_clear

  subroutine trace_field_length_single(seed,dL,max_steps,result,b_min)
    ! Trace one magnetic field line in both directions without storing its path.
    ! This first implementation supports uniform Cartesian grids and npe=1.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_length_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call trace_summary_seed(seed,dL,max_steps,result,b_min)
    else
      call trace_summary_seed(seed,dL,max_steps,result)
    endif
  end subroutine trace_field_length_single

  subroutine trace_field_twist_single(seed,dL,max_steps,result,b_min)
    ! Trace one magnetic field line and integrate twist without storing its path.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_twist_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call trace_summary_twist_seed(seed,dL,max_steps,result,b_min)
    else
      call trace_summary_twist_seed(seed,dL,max_steps,result)
    endif
  end subroutine trace_field_twist_single

  subroutine trace_field_mapping_single(seed,dL,max_steps,result,b_min, &
       source_normal)
    ! Trace one seed and return endpoint metadata for future mapping/Q work.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_mapping_result), intent(out) :: result
    double precision, intent(in), optional :: b_min
    double precision, intent(in), optional :: source_normal(3)

    if (present(b_min)) then
      if (present(source_normal)) then
        call trace_summary_mapping_seed(seed,dL,max_steps,result,b_min, &
             source_normal)
      else
        call trace_summary_mapping_seed(seed,dL,max_steps,result,b_min)
      endif
    else
      if (present(source_normal)) then
        call trace_summary_mapping_seed(seed,dL,max_steps,result, &
             source_normal=source_normal)
      else
        call trace_summary_mapping_seed(seed,dL,max_steps,result)
      endif
    endif
  end subroutine trace_field_mapping_single

  subroutine trace_field_qperp_single(seed,dL,max_steps,result,b_min)
    ! Trace one seed with Method-II tangent transport and return Q_perp.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_qperp_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call trace_qperp_single_core(seed,dL,max_steps,result,b_min)
    else
      call trace_qperp_single_core(seed,dL,max_steps,result)
    endif
  end subroutine trace_field_qperp_single

  subroutine trace_debug_sample_bhat_gradbhat(seed,bhat,grad_bhat,status,b_min)
    ! Debug/validation hook for the future Method-II sampler.
    double precision, intent(in) :: seed(ndim)
    double precision, intent(out) :: bhat(3),grad_bhat(3,3)
    integer, intent(out) :: status
    double precision, intent(in), optional :: b_min

    type(trace_length_result) :: located
    double precision :: field_min
    integer :: igrid

    bhat=zero
    grad_bhat=zero
    status=trace_status_active
    if (npe/=1) then
      status=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status=trace_status_unsupported_geometry
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_summary_locate_seed(seed,located,igrid)
    status=located%forward_status
    if (status/=trace_status_active) return

    call sample_bhat_gradbhat_at_point(seed,igrid,field_min,bhat, &
         grad_bhat,status)
  end subroutine trace_debug_sample_bhat_gradbhat

  subroutine trace_debug_sample_endpoint_B_face_limit(xhit,face_id,B,bhat, &
       status,b_min)
    ! Debug/validation hook for the Qperp endpoint face-limit B sampler.
    double precision, intent(in) :: xhit(ndim)
    integer, intent(in) :: face_id
    double precision, intent(out) :: B(ndim),bhat(ndim)
    integer, intent(out) :: status
    double precision, intent(in), optional :: b_min

    double precision :: field_min,xprobe(ndim),B3(3),bhat3(3)
    double precision :: xface,span
    integer :: normal_dim,side,igrid

    B=zero
    bhat=zero
    status=trace_status_active
    if (npe/=1) then
      status=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (.not.trace_face_is_boundary(face_id)) then
      status=trace_status_bad_face_limit_sample
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    {^IFTHREED
    select case(face_id)
    case(trace_face_xmin)
      normal_dim=1
      side=-1
      xface=xprobmin1
    case(trace_face_xmax)
      normal_dim=1
      side=1
      xface=xprobmax1
    case(trace_face_ymin)
      normal_dim=2
      side=-1
      xface=xprobmin2
    case(trace_face_ymax)
      normal_dim=2
      side=1
      xface=xprobmax2
    case(trace_face_zmin)
      normal_dim=3
      side=-1
      xface=xprobmin3
    case(trace_face_zmax)
      normal_dim=3
      side=1
      xface=xprobmax3
    case default
      status=trace_status_bad_face_limit_sample
      return
    end select

    span=max(abs(xprobmax1-xprobmin1), &
         max(abs(xprobmax2-xprobmin2),abs(xprobmax3-xprobmin3)))
    xprobe=xhit
    xprobe(normal_dim)=xface-side*1.d-6*max(one,span)
    call trace_debug_locate_point(xprobe,igrid,status)
    if (status/=trace_status_active) return

    call trace_sample_B_on_domain_face_limit(xhit,face_id,igrid,field_min, &
         B3,bhat3,status)
    if (status/=trace_status_active) return
    B=B3(1:ndim)
    bhat=bhat3(1:ndim)
    status=trace_status_boundary
    }
  end subroutine trace_debug_sample_endpoint_B_face_limit

  subroutine trace_debug_compare_gradbhat(seed,eps,bhat_current, &
       grad_current,status_current,bhat_xeps,grad_xeps,status_xeps,b_min)
    ! Debug/validation hook comparing production and x+-eps grad(bhat).
    double precision, intent(in) :: seed(ndim),eps
    double precision, intent(out) :: bhat_current(3),grad_current(3,3)
    double precision, intent(out) :: bhat_xeps(3),grad_xeps(3,3)
    integer, intent(out) :: status_current,status_xeps
    double precision, intent(in), optional :: b_min

    type(trace_length_result) :: located
    double precision :: field_min
    integer :: igrid

    bhat_current=zero
    grad_current=zero
    bhat_xeps=zero
    grad_xeps=zero
    status_current=trace_status_active
    status_xeps=trace_status_active
    if (npe/=1) then
      status_current=trace_status_mpi_unsupported
      status_xeps=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status_current=trace_status_unsupported_geometry
      status_xeps=trace_status_unsupported_geometry
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_summary_locate_seed(seed,located,igrid)
    status_current=located%forward_status
    if (status_current==trace_status_active) then
      call sample_bhat_gradbhat_at_point(seed,igrid,field_min, &
           bhat_current,grad_current,status_current)
    endif

    call sample_bhat_gradbhat_xeps_at_point(seed,eps,field_min,bhat_xeps, &
         grad_xeps,status_xeps)
  end subroutine trace_debug_compare_gradbhat

  subroutine trace_debug_compare_gradbhat_methods(seed,eps,bhat_current, &
       grad_current,status_current,bhat_xeps,grad_xeps,status_xeps, &
       bhat_interp,grad_interp,status_interp,b_min)
    ! Debug/validation hook comparing all grad(bhat) sampler candidates.
    double precision, intent(in) :: seed(ndim),eps
    double precision, intent(out) :: bhat_current(3),grad_current(3,3)
    double precision, intent(out) :: bhat_xeps(3),grad_xeps(3,3)
    double precision, intent(out) :: bhat_interp(3),grad_interp(3,3)
    integer, intent(out) :: status_current,status_xeps,status_interp
    double precision, intent(in), optional :: b_min

    type(trace_length_result) :: located
    double precision :: field_min
    integer :: igrid

    bhat_current=zero
    grad_current=zero
    bhat_xeps=zero
    grad_xeps=zero
    bhat_interp=zero
    grad_interp=zero
    status_current=trace_status_active
    status_xeps=trace_status_active
    status_interp=trace_status_active
    if (npe/=1) then
      status_current=trace_status_mpi_unsupported
      status_xeps=trace_status_mpi_unsupported
      status_interp=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status_current=trace_status_unsupported_geometry
      status_xeps=trace_status_unsupported_geometry
      status_interp=trace_status_unsupported_geometry
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_summary_locate_seed(seed,located,igrid)
    status_current=located%forward_status
    status_interp=located%forward_status
    if (status_current==trace_status_active) then
      call sample_bhat_gradbhat_cellfd_at_point(seed,igrid,field_min, &
           bhat_current,grad_current,status_current)
    endif
    if (status_interp==trace_status_active) then
      call sample_bhat_gradbhat_interpderiv_at_point(seed,igrid,field_min, &
           bhat_interp,grad_interp,status_interp)
    endif

    call sample_bhat_gradbhat_xeps_at_point(seed,eps,field_min,bhat_xeps, &
         grad_xeps,status_xeps)
  end subroutine trace_debug_compare_gradbhat_methods

  subroutine trace_debug_transport_tangents(seed,u0,v0,dL,nstep,forward, &
       x_end,u_end,v_end,bhat_end,u_perp_end,v_perp_end,status,b_min)
    ! Debug/validation hook for future Method-II tangent transport.
    double precision, intent(in) :: seed(ndim),u0(ndim),v0(ndim),dL
    integer, intent(in) :: nstep
    logical, intent(in) :: forward
    double precision, intent(out) :: x_end(ndim),u_end(ndim),v_end(ndim)
    double precision, intent(out) :: bhat_end(ndim)
    double precision, intent(out) :: u_perp_end(ndim),v_perp_end(ndim)
    integer, intent(out) :: status
    double precision, intent(in), optional :: b_min

    double precision :: field_min,h,bhat3(3),grad_end(3,3)
    integer :: igrid,istep

    x_end=seed
    u_end=u0
    v_end=v0
    bhat_end=zero
    u_perp_end=zero
    v_perp_end=zero

    status=trace_status_active
    if (npe/=1) then
      status=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status=trace_status_unsupported_geometry
      return
    else if (dL<=zero .or. nstep<0) then
      status=trace_status_invalid_input
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    h=dL
    if (.not.forward) h=-dL

    call trace_debug_locate_point(x_end,igrid,status)
    if (status/=trace_status_active) return

    do istep=1,nstep
      call trace_advance_tangent_state_rk2(x_end,u_end,v_end,igrid,h, &
           field_min,status)
      if (status/=trace_status_active) exit
    enddo

    if (status==trace_status_active) then
      {^IFTHREED
      call sample_bhat_gradbhat_at_point(x_end,igrid,field_min,bhat3, &
           grad_end,status)
      if (status==trace_status_active) bhat_end=bhat3(1:ndim)
      }
    endif
    if (status/=trace_status_active) return

    call trace_project_to_perp_bhat(u_end,bhat_end,u_perp_end)
    call trace_project_to_perp_bhat(v_end,bhat_end,v_perp_end)
  end subroutine trace_debug_transport_tangents

  subroutine trace_debug_transport_tangents_to_boundary(seed,u0,v0,dL, &
       max_steps,forward,x_end,u_end,v_end,bhat_end,b_end,u_perp_end, &
       v_perp_end,u_face_end,v_face_end,length,face_id,status,b_min)
    ! Debug/validation hook for Method-II style tangent transport to boundary.
    double precision, intent(in) :: seed(ndim),u0(ndim),v0(ndim),dL
    integer, intent(in) :: max_steps
    logical, intent(in) :: forward
    double precision, intent(out) :: x_end(ndim),u_end(ndim),v_end(ndim)
    double precision, intent(out) :: bhat_end(ndim),b_end(ndim)
    double precision, intent(out) :: u_perp_end(ndim),v_perp_end(ndim)
    double precision, intent(out) :: u_face_end(ndim),v_face_end(ndim)
    double precision, intent(out) :: length
    integer, intent(out) :: face_id,status
    double precision, intent(in), optional :: b_min

    double precision :: field_min,h,h_partial,partial_length
    double precision :: xprobe(ndim),xtrial(ndim),utrial(ndim),vtrial(ndim)
    double precision :: xhit(ndim),xold(ndim),uold(ndim),vold(ndim)
    double precision :: kx1(ndim),ku1(ndim),kv1(ndim)
    integer :: igrid,igrid_trial,istep,point_domain
    logical :: hit_ok

    x_end=seed
    u_end=u0
    v_end=v0
    bhat_end=zero
    b_end=zero
    u_perp_end=zero
    v_perp_end=zero
    u_face_end=zero
    v_face_end=zero
    length=zero
    face_id=trace_face_none

    status=trace_status_active
    if (npe/=1) then
      status=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.slab_uniform) then
      status=trace_status_unsupported_geometry
      return
    else if (dL<=zero .or. max_steps<0) then
      status=trace_status_invalid_input
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    h=dL
    if (.not.forward) h=-dL

    call trace_debug_locate_point(x_end,igrid,status)
    if (status/=trace_status_active) return

    do istep=1,max_steps
      xold=x_end
      uold=u_end
      vold=v_end

      call trace_tangent_rhs(xold,uold,vold,igrid,field_min,kx1,ku1,kv1, &
           status)
      if (status/=trace_status_active) return
      xprobe=xold+h*kx1

      point_domain=0
      {if (xprobe(^DB)>=xprobmin^DB .and. xprobe(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        call trace_intersect_domain(xold,xprobe,xhit,hit_ok,face_id)
        if (.not.hit_ok) then
          call trace_boundary_face_at_point(xold,face_id,hit_ok)
          if (.not.hit_ok) then
            status=trace_status_boundary
            x_end=xold
            return
          endif
          x_end=xold
          u_end=uold
          v_end=vold
          status=trace_status_boundary
          exit
        endif
        partial_length=dsqrt(sum((xhit-xold)**2))
        if (partial_length<=100.d0*epsilon(one)*max(one,abs(h))) then
          x_end=xhit
          u_end=uold
          v_end=vold
          status=trace_status_boundary
          exit
        endif
        h_partial=sign(partial_length,h)
        call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,igrid, &
             h_partial,field_min,kx1,ku1,kv1,xtrial,utrial,vtrial,status)
        if (status/=trace_status_active) return
        x_end=xhit
        u_end=utrial
        v_end=vtrial
        length=length+partial_length
        status=trace_status_boundary
        exit
      endif

      call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,igrid,h, &
           field_min,kx1,ku1,kv1,xtrial,utrial,vtrial,status)
      if (status/=trace_status_active) return

      point_domain=0
      {if (xtrial(^DB)>=xprobmin^DB .and. xtrial(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain==ndim) then
        call trace_locate_point_with_hint(xtrial,igrid,igrid_trial,status)
        if (status/=trace_status_active) return
        x_end=xtrial
        u_end=utrial
        v_end=vtrial
        igrid=igrid_trial
        length=length+abs(h)
      else
        call trace_intersect_domain(xold,xtrial,xhit,hit_ok,face_id)
        if (.not.hit_ok) then
          call trace_boundary_face_at_point(xold,face_id,hit_ok)
          if (.not.hit_ok) then
            status=trace_status_boundary
            x_end=xold
            return
          endif
          x_end=xold
          u_end=uold
          v_end=vold
          status=trace_status_boundary
          exit
        endif
        partial_length=dsqrt(sum((xhit-xold)**2))
        if (partial_length<=100.d0*epsilon(one)*max(one,abs(h))) then
          x_end=xhit
          u_end=uold
          v_end=vold
          status=trace_status_boundary
          exit
        endif
        h_partial=sign(partial_length,h)
        call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,igrid, &
             h_partial,field_min,kx1,ku1,kv1,xtrial,utrial,vtrial,status)
        if (status/=trace_status_active) return
        x_end=xhit
        u_end=utrial
        v_end=vtrial
        length=length+partial_length
        status=trace_status_boundary
        exit
      endif
    enddo

    if (status==trace_status_active) then
      status=trace_status_max_steps
      return
    endif
    if (status/=trace_status_boundary) return

    call trace_endpoint_B_bhat(x_end,face_id,igrid,field_min,b_end, &
         bhat_end,status)
    if (status/=trace_status_boundary) return

    call trace_project_to_perp_bhat(u_end,bhat_end,u_perp_end)
    call trace_project_to_perp_bhat(v_end,bhat_end,v_perp_end)
    call trace_project_to_boundary_face(u_end,bhat_end,face_id,u_face_end)
    call trace_project_to_boundary_face(v_end,bhat_end,face_id,v_face_end)
  end subroutine trace_debug_transport_tangents_to_boundary

  subroutine trace_debug_qperp_single(seed,dL,max_steps,qperp,logqperp, &
       x_f,x_b,u_f_perp,v_f_perp,u_b_perp,v_b_perp,b_seed,b_f,b_b, &
       length_f,length_b,face_f,face_b,status,N2,bfactor,b_min)
    ! Debug/validation hook for single-seed Method-II Q_perp diagnostics.
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    double precision, intent(out) :: qperp,logqperp,N2,bfactor
    double precision, intent(out) :: x_f(ndim),x_b(ndim)
    double precision, intent(out) :: u_f_perp(ndim),v_f_perp(ndim)
    double precision, intent(out) :: u_b_perp(ndim),v_b_perp(ndim)
    double precision, intent(out) :: b_seed(ndim),b_f(ndim),b_b(ndim)
    double precision, intent(out) :: length_f,length_b
    integer, intent(out) :: face_f,face_b,status
    double precision, intent(in), optional :: b_min

    type(trace_qperp_result) :: result

    if (present(b_min)) then
      call trace_field_qperp_single(seed,dL,max_steps,result,b_min)
    else
      call trace_field_qperp_single(seed,dL,max_steps,result)
    endif

    qperp=result%qperp
    logqperp=result%logqperp
    N2=result%N2
    bfactor=result%bfactor
    x_f=result%forward_endpoint
    x_b=result%backward_endpoint
    u_f_perp=result%u_forward_perp
    v_f_perp=result%v_forward_perp
    u_b_perp=result%u_backward_perp
    v_b_perp=result%v_backward_perp
    b_seed=result%B_seed
    b_f=result%forward_B
    b_b=result%backward_B
    length_f=result%forward_length
    length_b=result%backward_length
    face_f=result%forward_face
    face_b=result%backward_face
    status=result%status
  end subroutine trace_debug_qperp_single

  subroutine trace_init_qperp_result(seed,result)
    double precision, intent(in) :: seed(ndim)
    type(trace_qperp_result), intent(out) :: result

    result%seed=seed
    result%qperp=trace_debug_nan()
    result%logqperp=trace_debug_nan()
    result%N2=trace_debug_nan()
    result%bfactor=trace_debug_nan()
    result%q0=trace_debug_nan()
    result%logq0=trace_debug_nan()
    result%qperp0=trace_debug_nan()
    result%logqperp0=trace_debug_nan()
    result%N2_qperp0=trace_debug_nan()
    result%bfactor_qperp0=trace_debug_nan()
    result%forward_Bn_q0=trace_debug_nan()
    result%backward_Bn_q0=trace_debug_nan()
    result%B_seed=zero
    result%bhat_seed=zero
    result%forward_endpoint=seed
    result%backward_endpoint=seed
    result%forward_B=zero
    result%backward_B=zero
    result%forward_bhat=zero
    result%backward_bhat=zero
    result%forward_length=zero
    result%backward_length=zero
    result%forward_nstep=0
    result%backward_nstep=0
    result%u0=zero
    result%v0=zero
    result%u_forward_perp=zero
    result%v_forward_perp=zero
    result%u_backward_perp=zero
    result%v_backward_perp=zero
    result%forward_face=trace_face_none
    result%backward_face=trace_face_none
    result%forward_status=trace_status_active
    result%backward_status=trace_status_active
    result%status_q0=trace_status_active
    result%status_qperp0=trace_status_active
    result%status=trace_status_active
    result%valid=.false.
    result%valid_q0=.false.
    result%valid_qperp0=.false.
  end subroutine trace_init_qperp_result

  subroutine trace_qperp_single_core(seed,dL,max_steps,result,b_min)
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_qperp_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    double precision :: field_min,B3(3),Bseed_norm
    double precision :: u_end(ndim),v_end(ndim),u_face(ndim),v_face(ndim)
    integer :: igrid,status,indomain

    call trace_init_qperp_result(seed,result)

    if (npe/=1) then
      result%status=trace_status_mpi_unsupported
      return
    else if (ndim/=3 .or. .not.trace_cartesian_like_geometry()) then
      result%status=trace_status_unsupported_geometry
      return
    else if (dL<=zero .or. max_steps<=0) then
      result%status=trace_status_invalid_input
      return
    endif

    indomain=0
    {if (seed(^DB)>=xprobmin^DB .and. seed(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain/=ndim) then
      result%status=trace_status_seed_outside
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_debug_locate_point(seed,igrid,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif
    call sample_B_at_point(seed,igrid,field_min,B3,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif
    result%B_seed=B3(1:ndim)
    Bseed_norm=dsqrt(sum(result%B_seed**2))
    if (Bseed_norm<=zero .or. Bseed_norm<field_min) then
      result%status=trace_status_weak_field
      return
    endif
    result%bhat_seed=result%B_seed/Bseed_norm
    call trace_make_perp_basis(B3/Bseed_norm,result%u0,result%v0,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif

    call trace_debug_transport_tangents_to_boundary(seed,result%u0, &
         result%v0,dL,max_steps,.true.,result%forward_endpoint,u_end, &
         v_end,result%forward_bhat,result%forward_B,result%u_forward_perp, &
         result%v_forward_perp,u_face,v_face,result%forward_length, &
         result%forward_face,status,b_min=field_min)
    result%forward_status=status
    if (status/=trace_status_boundary) then
      result%status=status
      return
    endif

    call trace_debug_transport_tangents_to_boundary(seed,result%u0, &
         result%v0,dL,max_steps,.false.,result%backward_endpoint,u_end, &
         v_end,result%backward_bhat,result%backward_B,result%u_backward_perp, &
         result%v_backward_perp,u_face,v_face,result%backward_length, &
         result%backward_face,status,b_min=field_min)
    result%backward_status=status
    if (status/=trace_status_boundary) then
      result%status=status
      return
    endif

    call trace_qperp_compute_scalars(result,field_min)
  end subroutine trace_qperp_single_core

  subroutine trace_tangent_state_init(seed,u0,v0,igrid,seed_id,direction, &
       state,accumulate_twist)
    double precision, intent(in) :: seed(ndim),u0(ndim),v0(ndim)
    integer, intent(in) :: igrid,seed_id,direction
    type(trace_tangent_state), intent(out) :: state
    logical, intent(in), optional :: accumulate_twist

    state%seed_id=seed_id
    state%direction=direction
    state%active=.true.
    state%complete=.false.
    state%has_extra=.false.
    state%x=seed
    state%u=u0
    state%v=v0
    state%p=zero
    state%q=zero
    state%length=zero
    state%twist=zero
    state%nstep=0
    state%igrid=igrid
    state%status=trace_status_active
    state%twist_status=trace_status_active
    state%face_id=trace_face_none
    state%endpoint=seed
    state%endpoint_B=zero
    state%endpoint_bhat=zero
    state%u_perp=zero
    state%v_perp=zero
    state%p_perp=zero
    state%q_perp=zero
    state%rk45_h=zero
    state%accumulate_twist=.false.
    if (present(accumulate_twist)) state%accumulate_twist=accumulate_twist
    state%sph_cache%valid=.false.
    if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) then
      call trace_spherical_profile_add_count(trace_profile_directions,1_8)
      call trace_spherical_profile_add_count(trace_profile_tangent_systems,1_8)
    endif
  end subroutine trace_tangent_state_init

  subroutine trace_tangent_trace_states_grouped(states,nstate,dL,max_steps, &
       threshold,trace_mode)
    ! One grouped-state/OpenMP driver for RK2, RK2 short-boundary, and RK45
    ! Cartesian tangent tracing.  The mode selects only the single-state
    ! advance kernel, preserving the distinct numerical semantics.
    integer, intent(in) :: nstate,max_steps
    integer, intent(in), optional :: trace_mode
    type(trace_tangent_state), intent(inout) :: states(nstate)
    double precision, intent(in) :: dL,threshold

    logical, allocatable :: processed(:)
    double precision :: t0,tgroup0,tgroup1
    integer :: active_count,group_size,istate,jstate,target_grid,group_mode
    logical :: any_active
    logical :: profile_tangent

    group_mode=trace_tangent_group_rk2_general
    if (present(trace_mode)) group_mode=trace_mode
    profile_tangent=(group_mode/=trace_tangent_group_rk45_cartesian .and. &
         trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical)
    if (profile_tangent) t0=trace_spherical_profile_time()

    if (trace_cartesian_like_geometry()) then
      ! Cartesian interpolation/finalization is thread-safe across grids, and
      ! each state advance stops when that state leaves its current grid.  A
      ! single active-state parallel loop avoids thousands of tiny grid-group
      ! scans/parallel-region launches on AMR products.
      if (trace_rk2_stats_enabled) call cpu_time(tgroup0)
      do
        any_active=.false.
        active_count=0
        do istate=1,nstate
          if (states(istate)%active) then
            any_active=.true.
            active_count=active_count+1
          endif
        enddo
        if (.not.any_active) exit
        call trace_rk2_stats_note_group_iteration()
        call trace_rk2_stats_note_grid_group(active_count)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jstate) SCHEDULE(DYNAMIC,16)
        do jstate=1,nstate
          if (.not.states(jstate)%active) cycle
          select case (group_mode)
          case (trace_tangent_group_rk2_general)
            call trace_tangent_state_advance_in_grid(states(jstate), &
                 states(jstate)%igrid,dL,max_steps,threshold)
          case (trace_tangent_group_rk2_short_boundary)
            call trace_tangent_state_advance_in_grid(states(jstate), &
                 states(jstate)%igrid,dL,max_steps,threshold, &
                 short_boundary=.true.)
          case (trace_tangent_group_rk45_cartesian)
            call trace_tangent_state_advance_in_grid_rk45_cartesian( &
                 states(jstate),states(jstate)%igrid,dL,max_steps,threshold)
          case default
            states(jstate)%status=trace_status_invalid_input
            states(jstate)%active=.false.
            states(jstate)%complete=.true.
          end select
        enddo
        !$OMP END PARALLEL DO
      enddo
      if (trace_rk2_stats_enabled) then
        call cpu_time(tgroup1)
        call trace_rk2_stats_add_grouped_time(tgroup1-tgroup0)
      endif
      return
    endif

    if (trace_rk2_stats_enabled) call cpu_time(tgroup0)
    allocate(processed(nstate))
    do
      any_active=.false.
      do istate=1,nstate
        if (states(istate)%active) then
          any_active=.true.
          exit
        endif
      enddo
      if (.not.any_active) exit
      call trace_rk2_stats_note_group_iteration()

      processed=.false.
      do istate=1,nstate
        if (.not.states(istate)%active .or. processed(istate)) cycle
        target_grid=states(istate)%igrid
        if (trace_rk2_stats_enabled) then
          group_size=0
          do jstate=1,nstate
            if (states(jstate)%active .and. .not.processed(jstate) .and. &
                 states(jstate)%igrid==target_grid) group_size=group_size+1
          enddo
          call trace_rk2_stats_note_grid_group(group_size)
        endif
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jstate) SCHEDULE(DYNAMIC,16)
        do jstate=1,nstate
          if (.not.states(jstate)%active .or. processed(jstate)) cycle
          if (states(jstate)%igrid/=target_grid) cycle
          select case (group_mode)
          case (trace_tangent_group_rk2_general)
            call trace_tangent_state_advance_in_grid(states(jstate), &
                 target_grid,dL,max_steps,threshold)
          case (trace_tangent_group_rk2_short_boundary)
            call trace_tangent_state_advance_in_grid(states(jstate), &
                 target_grid,dL,max_steps,threshold,short_boundary=.true.)
          case (trace_tangent_group_rk45_cartesian)
            call trace_tangent_state_advance_in_grid_rk45_cartesian( &
                 states(jstate),target_grid,dL,max_steps,threshold)
          case default
            states(jstate)%status=trace_status_invalid_input
            states(jstate)%active=.false.
            states(jstate)%complete=.true.
          end select
          processed(jstate)=.true.
        enddo
        !$OMP END PARALLEL DO
      enddo
    enddo
    deallocate(processed)
    if (trace_rk2_stats_enabled) then
      call cpu_time(tgroup1)
      call trace_rk2_stats_add_grouped_time(tgroup1-tgroup0)
    endif
    if (profile_tangent) call trace_spherical_profile_add_time( &
         trace_profile_time_tangent,trace_spherical_profile_time()-t0)
  end subroutine trace_tangent_trace_states_grouped

  subroutine trace_tangent_state_advance_in_grid_rk45_cartesian(state, &
       igrid,dL,max_steps,threshold)
    type(trace_tangent_state), intent(inout) :: state
    integer, intent(in) :: igrid,max_steps
    double precision, intent(in) :: dL,threshold

    integer :: status

    if (.not.state%active) return

    do while(state%active .and. state%igrid==igrid .and. &
         state%nstep<max_steps)
      call trace_tangent_rk45_cartesian_step(state,dL,threshold,status)
      if (status/=trace_status_active) then
        state%status=status
        state%active=.false.
        state%complete=.true.
        return
      endif
    enddo

    if (state%active .and. state%nstep>=max_steps) then
      state%status=trace_status_max_steps
      state%active=.false.
      state%complete=.true.
    endif
  end subroutine trace_tangent_state_advance_in_grid_rk45_cartesian

  subroutine trace_tangent_rk45_cartesian_step(state,dL,threshold,status)
    type(trace_tangent_state), intent(inout) :: state
    double precision, intent(in) :: dL,threshold
    integer, intent(out) :: status

    double precision, parameter :: b21=1.d0/5.d0
    double precision, parameter :: b31=3.d0/40.d0,b32=9.d0/40.d0
    double precision, parameter :: b41=3.d0/10.d0,b42=-9.d0/10.d0, &
         b43=6.d0/5.d0
    double precision, parameter :: b51=-11.d0/54.d0,b52=5.d0/2.d0, &
         b53=-70.d0/27.d0,b54=35.d0/27.d0
    double precision, parameter :: b61=1631.d0/55296.d0, &
         b62=175.d0/512.d0,b63=575.d0/13824.d0, &
         b64=44275.d0/110592.d0,b65=253.d0/4096.d0
    double precision, parameter :: c1=37.d0/378.d0,c3=250.d0/621.d0, &
         c4=125.d0/594.d0,c6=512.d0/1771.d0
    double precision, parameter :: cs1=2825.d0/27648.d0, &
         cs3=18575.d0/48384.d0,cs4=13525.d0/55296.d0, &
         cs5=277.d0/14336.d0,cs6=one/4.d0

    double precision :: xold(ndim),uold(ndim),vold(ndim)
    double precision :: kx1(ndim),kx2(ndim),kx3(ndim),kx4(ndim)
    double precision :: kx5(ndim),kx6(ndim)
    double precision :: ku1(ndim),ku2(ndim),ku3(ndim),ku4(ndim)
    double precision :: ku5(ndim),ku6(ndim)
    double precision :: kv1(ndim),kv2(ndim),kv3(ndim),kv4(ndim)
    double precision :: kv5(ndim),kv6(ndim)
    double precision :: x2(ndim),x3(ndim),x4(ndim),x5s(ndim),x6(ndim)
    double precision :: u2(ndim),u3(ndim),u4(ndim),u5s(ndim),u6(ndim)
    double precision :: v2(ndim),v3(ndim),v4(ndim),v5s(ndim),v6(ndim)
    double precision :: x5sol(ndim),x4err(ndim)
    double precision :: u5sol(ndim),v5sol(ndim)
    double precision :: u4err(ndim),v4err(ndim)
    double precision :: h,hmax,hfloor,hcell,tol,err,grow,hsigned
    double precision :: pos_err_ratio,u_err_ratio,v_err_ratio
    double precision :: tan_err_ratio,control_ratio
    double precision :: tw1,tw3,tw4,tw6,twist_increment
    double precision :: dxb^D
    integer :: igrid1,igrid2,igrid3,igrid4,igrid5,igrid6
    integer :: iter,point_domain,reject_reason,twist_status
    logical :: boundary_finished,pos_ok,tan_ok

    status=trace_status_active
    if (.not.trace_cartesian_like_geometry()) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (dL<=zero) then
      status=trace_status_invalid_input
      return
    endif

    xold=state%x
    uold=state%u
    vold=state%v
    ^D&dxb^D=rnode(rpdx^D_,state%igrid);
    hmax=abs(trace_effective_step(xold,dL,state%igrid,dxb^D))
    if (hmax<=zero) then
      status=trace_status_invalid_input
      return
    endif
    h=state%rk45_h
    if (h<=zero) h=hmax
    h=min(h,hmax)
    hfloor=max(trace_step_min,100.d0*epsilon(one)*max(one,hmax))
    h=max(h,hfloor)

    do iter=1,100
      h=max(hfloor,min(h,hmax))
      hsigned=h
      if (state%direction<0) hsigned=-h

      call trace_tangent_rhs_cartesian_located(xold,uold,vold, &
           state%igrid,threshold,kx1,ku1,kv1,igrid1,status)
      if (status/=trace_status_active) return

      x2=xold+hsigned*b21*kx1
      u2=uold+hsigned*b21*ku1
      v2=vold+hsigned*b21*kv1
      call trace_tangent_rhs_cartesian_located(x2,u2,v2,igrid1, &
           threshold,kx2,ku2,kv2,igrid2,status)
      if (status/=trace_status_active) then
        if (status==trace_status_weak_field) return
        call trace_tangent_rk45_cartesian_finish_boundary(state,xold,uold, &
             vold,x2,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished)
        if (boundary_finished) return
        reject_reason=trace_rk45_reject_stage_failure
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        if (trace_rk45_stats_enabled) &
             call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif

      x3=xold+hsigned*(b31*kx1+b32*kx2)
      u3=uold+hsigned*(b31*ku1+b32*ku2)
      v3=vold+hsigned*(b31*kv1+b32*kv2)
      call trace_tangent_rhs_cartesian_located(x3,u3,v3,igrid2, &
           threshold,kx3,ku3,kv3,igrid3,status)
      if (status/=trace_status_active) then
        if (status==trace_status_weak_field) return
        call trace_tangent_rk45_cartesian_finish_boundary(state,xold,uold, &
             vold,x3,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished)
        if (boundary_finished) return
        reject_reason=trace_rk45_reject_stage_failure
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        if (trace_rk45_stats_enabled) &
             call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif

      x4=xold+hsigned*(b41*kx1+b42*kx2+b43*kx3)
      u4=uold+hsigned*(b41*ku1+b42*ku2+b43*ku3)
      v4=vold+hsigned*(b41*kv1+b42*kv2+b43*kv3)
      call trace_tangent_rhs_cartesian_located(x4,u4,v4,igrid3, &
           threshold,kx4,ku4,kv4,igrid4,status)
      if (status/=trace_status_active) then
        if (status==trace_status_weak_field) return
        call trace_tangent_rk45_cartesian_finish_boundary(state,xold,uold, &
             vold,x4,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished)
        if (boundary_finished) return
        reject_reason=trace_rk45_reject_stage_failure
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        if (trace_rk45_stats_enabled) &
             call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif

      x5s=xold+hsigned*(b51*kx1+b52*kx2+b53*kx3+b54*kx4)
      u5s=uold+hsigned*(b51*ku1+b52*ku2+b53*ku3+b54*ku4)
      v5s=vold+hsigned*(b51*kv1+b52*kv2+b53*kv3+b54*kv4)
      call trace_tangent_rhs_cartesian_located(x5s,u5s,v5s,igrid4, &
           threshold,kx5,ku5,kv5,igrid5,status)
      if (status/=trace_status_active) then
        if (status==trace_status_weak_field) return
        call trace_tangent_rk45_cartesian_finish_boundary(state,xold,uold, &
             vold,x5s,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished)
        if (boundary_finished) return
        reject_reason=trace_rk45_reject_stage_failure
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        if (trace_rk45_stats_enabled) &
             call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif

      x6=xold+hsigned*(b61*kx1+b62*kx2+b63*kx3+b64*kx4+b65*kx5)
      u6=uold+hsigned*(b61*ku1+b62*ku2+b63*ku3+b64*ku4+b65*ku5)
      v6=vold+hsigned*(b61*kv1+b62*kv2+b63*kv3+b64*kv4+b65*kv5)
      call trace_tangent_rhs_cartesian_located(x6,u6,v6,igrid5, &
           threshold,kx6,ku6,kv6,igrid6,status)
      if (status/=trace_status_active) then
        if (status==trace_status_weak_field) return
        call trace_tangent_rk45_cartesian_finish_boundary(state,xold,uold, &
             vold,x6,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished)
        if (boundary_finished) return
        reject_reason=trace_rk45_reject_stage_failure
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        if (trace_rk45_stats_enabled) &
             call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif

      x5sol=xold+hsigned*(c1*kx1+c3*kx3+c4*kx4+c6*kx6)
      u5sol=uold+hsigned*(c1*ku1+c3*ku3+c4*ku4+c6*ku6)
      v5sol=vold+hsigned*(c1*kv1+c3*kv3+c4*kv4+c6*kv6)
      x4err=xold+hsigned*(cs1*kx1+cs3*kx3+cs4*kx4+ &
           cs5*kx5+cs6*kx6)
      u4err=uold+hsigned*(cs1*ku1+cs3*ku3+cs4*ku4+ &
           cs5*ku5+cs6*ku6)
      v4err=vold+hsigned*(cs1*kv1+cs3*kv3+cs4*kv4+ &
           cs5*kv5+cs6*kv6)

      point_domain=0
      {if (x5sol(^DB)>=xprobmin^DB .and. x5sol(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        call trace_tangent_rk45_cartesian_finish_boundary(state,xold,uold, &
             vold,x5sol,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished)
        if (boundary_finished) return
        if (trace_rk45_stats_enabled) &
             call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             trace_rk45_reject_boundary)
        status=trace_status_out_of_domain
        return
      endif

      call trace_cartesian_local_cell_size(xold,igrid1,hcell,status)
      if (status/=trace_status_active .or. hcell<=zero) hcell=hmax
      err=dsqrt(sum((x5sol-x4err)**2))
      tol=trace_rk45_atol+trace_rk45_rtol*max(h,hcell)
      if (tol>zero) then
        pos_err_ratio=err/tol
      else
        pos_err_ratio=zero
      endif
      u_err_ratio=dsqrt(sum((u5sol-u4err)**2))/ &
           max(dsqrt(sum(u5sol**2)),trace_rk45_tangent_floor)
      v_err_ratio=dsqrt(sum((v5sol-v4err)**2))/ &
           max(dsqrt(sum(v5sol**2)),trace_rk45_tangent_floor)
      tan_err_ratio=max(u_err_ratio,v_err_ratio)
      if (trace_rk45_stats_enabled) &
           call trace_rk45_stats_note_tangent_error(pos_err_ratio, &
           u_err_ratio,v_err_ratio)
      pos_ok=(err<=tol)
      tan_ok=(tan_err_ratio<=trace_rk45_tangent_rtol)
      if ((pos_ok .and. tan_ok) .or. h<=hfloor*(one+epsilon(one))) then
        if (state%accumulate_twist .and. &
             state%twist_status==trace_status_active) then
          twist_status=trace_status_active
          call trace_twist_density_at_point(xold,igrid1,threshold,tw1, &
               twist_status)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x3,igrid3,threshold,tw3, &
               twist_status)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x4,igrid4,threshold,tw4, &
               twist_status)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x6,igrid6,threshold,tw6, &
               twist_status)
          if (twist_status==trace_status_active) then
            twist_increment=h*(c1*tw1+c3*tw3+c4*tw4+c6*tw6)
            state%twist=state%twist+twist_increment
          else
            state%twist_status=twist_status
          endif
        endif
        state%x=x5sol
        state%u=u5sol
        state%v=v5sol
        state%length=state%length+h
        state%nstep=state%nstep+1
        call trace_locate_point_with_hint(state%x,igrid6,state%igrid,status)
        if (status/=trace_status_active) return
        control_ratio=max(pos_err_ratio, &
             tan_err_ratio/trace_rk45_tangent_rtol)
        if (control_ratio>zero) then
          grow=trace_rk45_safety*(one/control_ratio)**0.2d0
          grow=max(trace_rk45_min_shrink,min(trace_rk45_max_grow,grow))
        else
          grow=trace_rk45_max_grow
        endif
        state%rk45_h=min(hmax,max(hfloor,h*grow))
        if (trace_rk45_stats_enabled) &
             call trace_rk45_stats_note_attempt(.true.,.false.,h)
        status=trace_status_active
        return
      endif

      if (.not.pos_ok) then
        grow=trace_rk45_safety*(tol/max(err,smalldouble))**0.25d0
        if (trace_rk45_stats_enabled) &
             call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             trace_rk45_reject_error)
      else
        grow=trace_rk45_safety*(trace_rk45_tangent_rtol/ &
             max(tan_err_ratio,smalldouble))**0.25d0
        if (trace_rk45_stats_enabled) &
             call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             trace_rk45_reject_tangent_error)
      endif
      grow=max(trace_rk45_min_shrink,min(one,grow))
      h=max(hfloor,h*grow)
    enddo

    call trace_tangent_rk45_cartesian_finish_boundary(state,xold,uold, &
         vold,xold+hsigned*kx1,hsigned,threshold,kx1,ku1,kv1,status, &
         boundary_finished)
    if (.not.boundary_finished) status=trace_status_out_of_domain
  end subroutine trace_tangent_rk45_cartesian_step

  subroutine trace_tangent_rhs_cartesian_located(x,u,v,igrid_hint, &
       threshold,kx,ku,kv,igrid,status)
    double precision, intent(in) :: x(ndim),u(ndim),v(ndim),threshold
    integer, intent(in) :: igrid_hint
    double precision, intent(out) :: kx(ndim),ku(ndim),kv(ndim)
    integer, intent(out) :: igrid,status

    kx=zero
    ku=zero
    kv=zero
    igrid=-1
    status=trace_status_unsupported_geometry
    if (.not.trace_cartesian_like_geometry()) return

    call trace_locate_point_with_hint(x,igrid_hint,igrid,status)
    if (status/=trace_status_active) return
    call trace_tangent_rhs(x,u,v,igrid,threshold,kx,ku,kv,status)
  end subroutine trace_tangent_rhs_cartesian_located

  subroutine trace_tangent_rk45_cartesian_finish_boundary(state,xold,uold, &
       vold,xstage,hsigned,threshold,kx1,ku1,kv1,status,finished)
    type(trace_tangent_state), intent(inout) :: state
    double precision, intent(in) :: xold(ndim),uold(ndim),vold(ndim)
    double precision, intent(in) :: xstage(ndim),hsigned,threshold
    double precision, intent(in) :: kx1(ndim),ku1(ndim),kv1(ndim)
    integer, intent(out) :: status
    logical, intent(out) :: finished

    double precision :: xhit(ndim),utrial(ndim),vtrial(ndim)
    double precision :: xtrial(ndim),partial_length,h_partial
    logical :: hit_ok

    finished=.false.
    call trace_intersect_domain(xold,xstage,xhit,hit_ok,state%face_id)
    if (.not.hit_ok) return

    partial_length=dsqrt(sum((xhit-xold)**2))
    if (partial_length>100.d0*epsilon(one)*max(one,abs(hsigned))) then
      h_partial=sign(partial_length,hsigned)
      call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid, &
           h_partial,threshold,kx1,ku1,kv1,xtrial,utrial,vtrial,status)
      if (status/=trace_status_active) return
      call trace_tangent_accumulate_twist(state,xhit,partial_length, &
           threshold)
      state%u=utrial
      state%v=vtrial
      state%length=state%length+partial_length
      state%nstep=state%nstep+1
    else
      state%u=uold
      state%v=vold
    endif
    state%x=xhit
    state%endpoint=xhit
    state%status=trace_status_boundary
    state%active=.false.
    state%complete=.true.
    state%rk45_h=zero
    call trace_tangent_state_finalize_boundary(state,threshold)
    call trace_rk45_stats_note_attempt(.true.,.true.,partial_length)
    status=trace_status_active
    finished=.true.
  end subroutine trace_tangent_rk45_cartesian_finish_boundary

  subroutine trace_tangent_trace_states_grouped_rk45_spherical(states, &
       nstate,dL,max_steps,threshold)
    integer, intent(in) :: nstate,max_steps
    type(trace_tangent_state), intent(inout) :: states(nstate)
    double precision, intent(in) :: dL,threshold

    logical, allocatable :: processed(:)
    integer :: istate,jstate,target_grid
    logical :: any_active

    allocate(processed(nstate))
    do
      any_active=.false.
      do istate=1,nstate
        if (states(istate)%active) then
          any_active=.true.
          exit
        endif
      enddo
      if (.not.any_active) exit

      processed=.false.
      do istate=1,nstate
        if (.not.states(istate)%active .or. processed(istate)) cycle
        target_grid=states(istate)%igrid
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jstate) SCHEDULE(DYNAMIC,16)
        do jstate=1,nstate
          if (.not.states(jstate)%active .or. processed(jstate)) cycle
          if (states(jstate)%igrid/=target_grid) cycle
          call trace_tangent_state_advance_in_grid_rk45_spherical( &
               states(jstate),target_grid,dL,max_steps,threshold)
          processed(jstate)=.true.
        enddo
        !$OMP END PARALLEL DO
      enddo
    enddo
    deallocate(processed)
  end subroutine trace_tangent_trace_states_grouped_rk45_spherical

  subroutine trace_tangent_state_advance_in_grid_rk45_spherical(state, &
       igrid,dL,max_steps,threshold)
    type(trace_tangent_state), intent(inout) :: state
    integer, intent(in) :: igrid,max_steps
    double precision, intent(in) :: dL,threshold

    integer :: status

    if (.not.state%active) return

    do while(state%active .and. state%igrid==igrid .and. &
         state%nstep<max_steps)
      call trace_tangent_rk45_spherical_step(state,dL,threshold,status)
      if (status/=trace_status_active) then
        state%status=status
        state%active=.false.
        state%complete=.true.
        return
      endif
    enddo

    if (state%active .and. state%nstep>=max_steps) then
      state%status=trace_status_max_steps
      state%active=.false.
      state%complete=.true.
    endif
    if ((.not.state%active) .and. trace_spherical_profile_enabled .and. &
         geo_coordinate==geo_spherical) then
      call trace_spherical_profile_note_trace_steps(state%nstep)
    endif
  end subroutine trace_tangent_state_advance_in_grid_rk45_spherical

  subroutine trace_tangent_rk45_spherical_step(state,dL,threshold,status)
    type(trace_tangent_state), intent(inout) :: state
    double precision, intent(in) :: dL,threshold
    integer, intent(out) :: status

    double precision, parameter :: b21=1.d0/5.d0
    double precision, parameter :: b31=3.d0/40.d0,b32=9.d0/40.d0
    double precision, parameter :: b41=3.d0/10.d0,b42=-9.d0/10.d0, &
         b43=6.d0/5.d0
    double precision, parameter :: b51=-11.d0/54.d0,b52=5.d0/2.d0, &
         b53=-70.d0/27.d0,b54=35.d0/27.d0
    double precision, parameter :: b61=1631.d0/55296.d0, &
         b62=175.d0/512.d0,b63=575.d0/13824.d0, &
         b64=44275.d0/110592.d0,b65=253.d0/4096.d0
    double precision, parameter :: c1=37.d0/378.d0,c3=250.d0/621.d0, &
         c4=125.d0/594.d0,c6=512.d0/1771.d0
    double precision, parameter :: cs1=2825.d0/27648.d0, &
         cs3=18575.d0/48384.d0,cs4=13525.d0/55296.d0, &
         cs5=277.d0/14336.d0,cs6=one/4.d0

    double precision :: xold(ndim),uold(ndim),vold(ndim)
    double precision :: pold(ndim),qold(ndim)
    double precision :: kx1(ndim),kx2(ndim),kx3(ndim),kx4(ndim)
    double precision :: kx5(ndim),kx6(ndim)
    double precision :: ku1(ndim),ku2(ndim),ku3(ndim),ku4(ndim)
    double precision :: ku5(ndim),ku6(ndim)
    double precision :: kv1(ndim),kv2(ndim),kv3(ndim),kv4(ndim)
    double precision :: kv5(ndim),kv6(ndim)
    double precision :: kp1(ndim),kp2(ndim),kp3(ndim),kp4(ndim)
    double precision :: kp5(ndim),kp6(ndim)
    double precision :: kq1(ndim),kq2(ndim),kq3(ndim),kq4(ndim)
    double precision :: kq5(ndim),kq6(ndim)
    double precision :: x2(ndim),x3(ndim),x4(ndim),x5s(ndim),x6(ndim)
    double precision :: u2(ndim),u3(ndim),u4(ndim),u5s(ndim),u6(ndim)
    double precision :: v2(ndim),v3(ndim),v4(ndim),v5s(ndim),v6(ndim)
    double precision :: p2(ndim),p3(ndim),p4(ndim),p5s(ndim),p6(ndim)
    double precision :: q2(ndim),q3(ndim),q4(ndim),q5s(ndim),q6(ndim)
    double precision :: x5sol(ndim),x4err(ndim)
    double precision :: u5sol(ndim),v5sol(ndim)
    double precision :: u4err(ndim),v4err(ndim)
    double precision :: p5sol(ndim),q5sol(ndim)
    double precision :: p4err(ndim),q4err(ndim)
    double precision :: h,hmax,hfloor,hcell,tol,err,grow,hsigned
    double precision :: pos_err_ratio,u_err_ratio,v_err_ratio
    double precision :: p_err_ratio,q_err_ratio
    double precision :: tan_err_ratio,control_ratio
    double precision :: tw1,tw3,tw4,tw6,twist_increment
    double precision :: r_metric,sin_theta_metric
    double precision :: domain_min(ndim),domain_max(ndim)
    type(trace_sph_interp_ctx) :: ctx1,ctx2,ctx3,ctx4,ctx5,ctx6
    integer :: igrid1,igrid2,igrid3,igrid4,igrid5,igrid6
    integer :: iter,point_domain,ctx_status,reject_reason,twist_status
    logical :: boundary_finished,pos_ok,tan_ok

    status=trace_status_active
    if (geo_coordinate/=geo_spherical) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (dL<=zero) then
      status=trace_status_invalid_input
      return
    endif

    xold=state%x
    uold=state%u
    vold=state%v
    pold=state%p
    qold=state%q
    call trace_spherical_interp_ctx_build_cached(xold,state%igrid, &
         state%sph_cache,ctx1,ctx_status)
    if (ctx_status/=trace_status_active) then
      status=ctx_status
      return
    endif
    hmax=abs(trace_spherical_effective_step_ctx(dL,ctx1))
    hcell=ctx1%h_local
    if (hmax<=zero) then
      status=trace_status_invalid_input
      return
    endif
    h=state%rk45_h
    if (h<=zero) h=hmax
    h=min(h,hmax)
    hfloor=max(trace_step_min,100.d0*epsilon(one)*max(one,hmax))
    h=max(h,hfloor)
    ^D&domain_min(^D)=xprobmin^D;
    ^D&domain_max(^D)=xprobmax^D;

    do iter=1,100
      h=max(hfloor,min(h,hmax))
      hsigned=h
      if (state%direction<0) hsigned=-h

      if (state%has_extra) then
        call trace_tangent_rhs(xold,uold,vold,ctx1%igrid,threshold,kx1, &
             ku1,kv1,status,pold,qold,kp1,kq1,sph_ctx=ctx1)
      else
        call trace_tangent_rhs(xold,uold,vold,ctx1%igrid,threshold,kx1, &
             ku1,kv1,status,sph_ctx=ctx1)
      endif
      if (status/=trace_status_active) return
      igrid1=ctx1%igrid

      x2=xold+hsigned*b21*kx1
      u2=uold+hsigned*b21*ku1
      v2=vold+hsigned*b21*kv1
      if (state%has_extra) then
        p2=pold+hsigned*b21*kp1
        q2=qold+hsigned*b21*kq1
      endif
      call trace_spherical_interp_ctx_build_cached(x2,igrid1, &
           state%sph_cache,ctx2,ctx_status)
      if (ctx_status/=trace_status_active) then
        call trace_tangent_rk45_spherical_finish_boundary(state,xold,uold, &
             vold,x2,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished,pold,qold,kp1,kq1)
        if (boundary_finished) return
        reject_reason=trace_rk45_reject_stage_failure
        if (ctx_status==trace_status_boundary .or. &
             ctx_status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h,reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif
      if (state%has_extra) then
        call trace_tangent_rhs(x2,u2,v2,ctx2%igrid,threshold,kx2,ku2, &
             kv2,status,p2,q2,kp2,kq2,sph_ctx=ctx2)
      else
        call trace_tangent_rhs(x2,u2,v2,ctx2%igrid,threshold,kx2,ku2, &
             kv2,status,sph_ctx=ctx2)
      endif
      if (status/=trace_status_active) return
      igrid2=ctx2%igrid

      x3=xold+hsigned*(b31*kx1+b32*kx2)
      u3=uold+hsigned*(b31*ku1+b32*ku2)
      v3=vold+hsigned*(b31*kv1+b32*kv2)
      if (state%has_extra) then
        p3=pold+hsigned*(b31*kp1+b32*kp2)
        q3=qold+hsigned*(b31*kq1+b32*kq2)
      endif
      call trace_spherical_interp_ctx_build_cached(x3,igrid2, &
           state%sph_cache,ctx3,ctx_status)
      if (ctx_status/=trace_status_active) then
        call trace_tangent_rk45_spherical_finish_boundary(state,xold,uold, &
             vold,x3,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished,pold,qold,kp1,kq1)
        if (boundary_finished) return
        reject_reason=trace_rk45_reject_stage_failure
        if (ctx_status==trace_status_boundary .or. &
             ctx_status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h,reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif
      if (state%has_extra) then
        call trace_tangent_rhs(x3,u3,v3,ctx3%igrid,threshold,kx3,ku3, &
             kv3,status,p3,q3,kp3,kq3,sph_ctx=ctx3)
      else
        call trace_tangent_rhs(x3,u3,v3,ctx3%igrid,threshold,kx3,ku3, &
             kv3,status,sph_ctx=ctx3)
      endif
      if (status/=trace_status_active) return
      igrid3=ctx3%igrid

      x4=xold+hsigned*(b41*kx1+b42*kx2+b43*kx3)
      u4=uold+hsigned*(b41*ku1+b42*ku2+b43*ku3)
      v4=vold+hsigned*(b41*kv1+b42*kv2+b43*kv3)
      if (state%has_extra) then
        p4=pold+hsigned*(b41*kp1+b42*kp2+b43*kp3)
        q4=qold+hsigned*(b41*kq1+b42*kq2+b43*kq3)
      endif
      call trace_spherical_interp_ctx_build_cached(x4,igrid3, &
           state%sph_cache,ctx4,ctx_status)
      if (ctx_status/=trace_status_active) then
        call trace_tangent_rk45_spherical_finish_boundary(state,xold,uold, &
             vold,x4,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished,pold,qold,kp1,kq1)
        if (boundary_finished) return
        reject_reason=trace_rk45_reject_stage_failure
        if (ctx_status==trace_status_boundary .or. &
             ctx_status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h,reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif
      if (state%has_extra) then
        call trace_tangent_rhs(x4,u4,v4,ctx4%igrid,threshold,kx4,ku4, &
             kv4,status,p4,q4,kp4,kq4,sph_ctx=ctx4)
      else
        call trace_tangent_rhs(x4,u4,v4,ctx4%igrid,threshold,kx4,ku4, &
             kv4,status,sph_ctx=ctx4)
      endif
      if (status/=trace_status_active) return
      igrid4=ctx4%igrid

      x5s=xold+hsigned*(b51*kx1+b52*kx2+b53*kx3+b54*kx4)
      u5s=uold+hsigned*(b51*ku1+b52*ku2+b53*ku3+b54*ku4)
      v5s=vold+hsigned*(b51*kv1+b52*kv2+b53*kv3+b54*kv4)
      if (state%has_extra) then
        p5s=pold+hsigned*(b51*kp1+b52*kp2+b53*kp3+b54*kp4)
        q5s=qold+hsigned*(b51*kq1+b52*kq2+b53*kq3+b54*kq4)
      endif
      call trace_spherical_interp_ctx_build_cached(x5s,igrid4, &
           state%sph_cache,ctx5,ctx_status)
      if (ctx_status/=trace_status_active) then
        call trace_tangent_rk45_spherical_finish_boundary(state,xold,uold, &
             vold,x5s,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished,pold,qold,kp1,kq1)
        if (boundary_finished) return
        reject_reason=trace_rk45_reject_stage_failure
        if (ctx_status==trace_status_boundary .or. &
             ctx_status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h,reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif
      if (state%has_extra) then
        call trace_tangent_rhs(x5s,u5s,v5s,ctx5%igrid,threshold,kx5,ku5, &
             kv5,status,p5s,q5s,kp5,kq5,sph_ctx=ctx5)
      else
        call trace_tangent_rhs(x5s,u5s,v5s,ctx5%igrid,threshold,kx5,ku5, &
             kv5,status,sph_ctx=ctx5)
      endif
      if (status/=trace_status_active) return
      igrid5=ctx5%igrid

      x6=xold+hsigned*(b61*kx1+b62*kx2+b63*kx3+b64*kx4+b65*kx5)
      u6=uold+hsigned*(b61*ku1+b62*ku2+b63*ku3+b64*ku4+b65*ku5)
      v6=vold+hsigned*(b61*kv1+b62*kv2+b63*kv3+b64*kv4+b65*kv5)
      if (state%has_extra) then
        p6=pold+hsigned*(b61*kp1+b62*kp2+b63*kp3+b64*kp4+b65*kp5)
        q6=qold+hsigned*(b61*kq1+b62*kq2+b63*kq3+b64*kq4+b65*kq5)
      endif
      call trace_spherical_interp_ctx_build_cached(x6,igrid5, &
           state%sph_cache,ctx6,ctx_status)
      if (ctx_status/=trace_status_active) then
        call trace_tangent_rk45_spherical_finish_boundary(state,xold,uold, &
             vold,x6,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished,pold,qold,kp1,kq1)
        if (boundary_finished) return
        reject_reason=trace_rk45_reject_stage_failure
        if (ctx_status==trace_status_boundary .or. &
             ctx_status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h,reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif
      if (state%has_extra) then
        call trace_tangent_rhs(x6,u6,v6,ctx6%igrid,threshold,kx6,ku6, &
             kv6,status,p6,q6,kp6,kq6,sph_ctx=ctx6)
      else
        call trace_tangent_rhs(x6,u6,v6,ctx6%igrid,threshold,kx6,ku6, &
             kv6,status,sph_ctx=ctx6)
      endif
      if (status/=trace_status_active) return
      igrid6=ctx6%igrid

      x5sol=xold+hsigned*(c1*kx1+c3*kx3+c4*kx4+c6*kx6)
      u5sol=uold+hsigned*(c1*ku1+c3*ku3+c4*ku4+c6*ku6)
      v5sol=vold+hsigned*(c1*kv1+c3*kv3+c4*kv4+c6*kv6)
      x4err=xold+hsigned*(cs1*kx1+cs3*kx3+cs4*kx4+ &
           cs5*kx5+cs6*kx6)
      u4err=uold+hsigned*(cs1*ku1+cs3*ku3+cs4*ku4+ &
           cs5*ku5+cs6*ku6)
      v4err=vold+hsigned*(cs1*kv1+cs3*kv3+cs4*kv4+ &
           cs5*kv5+cs6*kv6)
      if (state%has_extra) then
        p5sol=pold+hsigned*(c1*kp1+c3*kp3+c4*kp4+c6*kp6)
        q5sol=qold+hsigned*(c1*kq1+c3*kq3+c4*kq4+c6*kq6)
        p4err=pold+hsigned*(cs1*kp1+cs3*kp3+cs4*kp4+ &
             cs5*kp5+cs6*kp6)
        q4err=qold+hsigned*(cs1*kq1+cs3*kq3+cs4*kq4+ &
             cs5*kq5+cs6*kq6)
      endif

      point_domain=0
      {if (x5sol(^DB)>=domain_min(^DB) .and. x5sol(^DB)<domain_max(^DB)) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        call trace_tangent_rk45_spherical_finish_boundary(state,xold,uold, &
             vold,x5sol,hsigned,threshold,kx1,ku1,kv1,status, &
             boundary_finished,pold,qold,kp1,kq1)
        if (boundary_finished) return
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             trace_rk45_reject_boundary)
        status=trace_status_out_of_domain
        return
      endif

      err=zero
      {^IFTHREED
      r_metric=max(abs(xold(1)),smalldouble)
      sin_theta_metric=max(abs(dsin(xold(2))),smalldouble)
      err=dsqrt((x5sol(1)-x4err(1))**2+ &
           (r_metric*(x5sol(2)-x4err(2)))**2+ &
           (r_metric*sin_theta_metric*(x5sol(3)-x4err(3)))**2)
      }
      if (hcell<=zero) hcell=hmax
      tol=trace_rk45_atol+trace_rk45_rtol*max(h,hcell)
      if (tol>zero) then
        pos_err_ratio=err/tol
      else
        pos_err_ratio=zero
      endif
      u_err_ratio=dsqrt(sum((u5sol-u4err)**2))/ &
           max(dsqrt(sum(u5sol**2)),trace_rk45_tangent_floor)
      v_err_ratio=dsqrt(sum((v5sol-v4err)**2))/ &
           max(dsqrt(sum(v5sol**2)),trace_rk45_tangent_floor)
      if (state%has_extra) then
        p_err_ratio=dsqrt(sum((p5sol-p4err)**2))/ &
             max(dsqrt(sum(p5sol**2)),trace_rk45_tangent_floor)
        q_err_ratio=dsqrt(sum((q5sol-q4err)**2))/ &
             max(dsqrt(sum(q5sol**2)),trace_rk45_tangent_floor)
        u_err_ratio=max(u_err_ratio,p_err_ratio)
        v_err_ratio=max(v_err_ratio,q_err_ratio)
      endif
      tan_err_ratio=max(u_err_ratio,v_err_ratio)
      call trace_rk45_stats_note_tangent_error(pos_err_ratio, &
           u_err_ratio,v_err_ratio)
      pos_ok=(err<=tol)
      tan_ok=(tan_err_ratio<=trace_rk45_tangent_rtol)
      if ((pos_ok .and. tan_ok) .or. h<=hfloor*(one+epsilon(one))) then
        if (state%accumulate_twist .and. &
             state%twist_status==trace_status_active) then
          twist_status=trace_status_active
          call trace_twist_density_at_point(xold,ctx1%igrid,threshold,tw1, &
               twist_status,sph_cache=ctx1)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x3,ctx3%igrid,threshold, &
               tw3,twist_status,sph_cache=ctx3)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x4,ctx4%igrid,threshold, &
               tw4,twist_status,sph_cache=ctx4)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x6,ctx6%igrid,threshold, &
               tw6,twist_status,sph_cache=ctx6)
          if (twist_status==trace_status_active) then
            twist_increment=h*(c1*tw1+c3*tw3+c4*tw4+c6*tw6)
            state%twist=state%twist+twist_increment
          else
            state%twist_status=twist_status
          endif
        endif
        state%x=x5sol
        state%u=u5sol
        state%v=v5sol
        if (state%has_extra) then
          state%p=p5sol
          state%q=q5sol
        endif
        state%length=state%length+h
        state%nstep=state%nstep+1
        call trace_locate_point_with_hint(state%x,igrid6,state%igrid,status)
        if (status/=trace_status_active) return
        state%sph_cache%valid=.false.
        if (trace_spherical_profile_enabled) &
             call trace_spherical_profile_add_count(trace_profile_steps,1_8)
        control_ratio=max(pos_err_ratio, &
             tan_err_ratio/trace_rk45_tangent_rtol)
        if (control_ratio>zero) then
          grow=trace_rk45_safety*(one/control_ratio)**0.2d0
          grow=max(trace_rk45_min_shrink,min(trace_rk45_max_grow,grow))
        else
          grow=trace_rk45_max_grow
        endif
        state%rk45_h=min(hmax,max(hfloor,h*grow))
        call trace_rk45_stats_note_attempt(.true.,.false.,h)
        status=trace_status_active
        return
      endif

      if (.not.pos_ok) then
        grow=trace_rk45_safety*(tol/max(err,smalldouble))**0.25d0
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             trace_rk45_reject_error)
      else
        grow=trace_rk45_safety*(trace_rk45_tangent_rtol/ &
             max(tan_err_ratio,smalldouble))**0.25d0
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             trace_rk45_reject_tangent_error)
      endif
      grow=max(trace_rk45_min_shrink,min(one,grow))
      h=max(hfloor,h*grow)
    enddo

    call trace_tangent_rk45_spherical_finish_boundary(state,xold,uold, &
         vold,xold+hsigned*kx1,hsigned,threshold,kx1,ku1,kv1,status, &
         boundary_finished,pold,qold,kp1,kq1)
    if (.not.boundary_finished) status=trace_status_out_of_domain
  end subroutine trace_tangent_rk45_spherical_step

  subroutine trace_tangent_rk45_spherical_finish_boundary(state,xold,uold, &
       vold,xstage,hsigned,threshold,kx1,ku1,kv1,status,finished, &
       pold,qold,kp1,kq1)
    type(trace_tangent_state), intent(inout) :: state
    double precision, intent(in) :: xold(ndim),uold(ndim),vold(ndim)
    double precision, intent(in) :: xstage(ndim),hsigned,threshold
    double precision, intent(in) :: kx1(ndim),ku1(ndim),kv1(ndim)
    integer, intent(out) :: status
    logical, intent(out) :: finished
    double precision, intent(in), optional :: pold(ndim),qold(ndim)
    double precision, intent(in), optional :: kp1(ndim),kq1(ndim)

    double precision :: xhit(ndim),utrial(ndim),vtrial(ndim)
    double precision :: ptrial(ndim),qtrial(ndim)
    double precision :: xtrial(ndim),partial_length,h_partial,alpha_hit
    logical :: hit_ok

    finished=.false.
    call trace_intersect_domain(xold,xstage,xhit,hit_ok,state%face_id, &
         alpha_hit)
    if (.not.hit_ok) return

    partial_length=trace_segment_length(xold,xhit,hsigned,alpha_hit)
    if (partial_length>100.d0*epsilon(one)*max(one,abs(hsigned))) then
      h_partial=sign(partial_length,hsigned)
      if (state%has_extra) then
        if (.not.(present(pold) .and. present(qold) .and. present(kp1) &
             .and. present(kq1))) then
          status=trace_status_invalid_input
          return
        endif
        call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid, &
             h_partial,threshold,kx1,ku1,kv1,xtrial,utrial,vtrial,status, &
             pold,qold,kp1,kq1,ptrial,qtrial,sph_cache=state%sph_cache)
      else
        call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid, &
             h_partial,threshold,kx1,ku1,kv1,xtrial,utrial,vtrial,status, &
             sph_cache=state%sph_cache)
      endif
      if (status/=trace_status_active) return
      call trace_tangent_accumulate_twist(state,xhit,partial_length, &
           threshold)
      state%u=utrial
      state%v=vtrial
      if (state%has_extra) then
        state%p=ptrial
        state%q=qtrial
      endif
      state%length=state%length+partial_length
      state%nstep=state%nstep+1
      if (trace_spherical_profile_enabled) &
           call trace_spherical_profile_add_count(trace_profile_steps,1_8)
    else
      state%u=uold
      state%v=vold
      if (state%has_extra) then
        if (.not.(present(pold) .and. present(qold))) then
          status=trace_status_invalid_input
          return
        endif
        state%p=pold
        state%q=qold
      endif
    endif
    state%x=xhit
    state%endpoint=xhit
    state%status=trace_status_boundary
    state%active=.false.
    state%complete=.true.
    state%rk45_h=zero
    state%sph_cache%valid=.false.
    if (trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_count(trace_profile_boundary_events,1_8)
    endif
    call trace_tangent_state_finalize_boundary(state,threshold)
    call trace_rk45_stats_note_attempt(.true.,.true.,partial_length)
    status=trace_status_active
    finished=.true.
  end subroutine trace_tangent_rk45_spherical_finish_boundary

  subroutine trace_tangent_state_advance_in_grid(state,igrid,dL,max_steps, &
       threshold,short_boundary)
    type(trace_tangent_state), intent(inout) :: state
    integer, intent(in) :: igrid,max_steps
    double precision, intent(in) :: dL,threshold
    logical, intent(in), optional :: short_boundary

    integer, parameter :: max_bisect=50
    double precision :: h,h_partial,partial_length,alpha_hit
    double precision :: h_abs_low,h_abs_high,h_abs_mid,h_abs_final,hsign
    double precision :: xprobe(ndim),xtrial(ndim),utrial(ndim),vtrial(ndim)
    double precision :: ptrial(ndim),qtrial(ndim)
    double precision :: xhit(ndim),xold(ndim),uold(ndim),vold(ndim)
    double precision :: pold(ndim),qold(ndim)
    double precision :: xlow(ndim),ulow(ndim),vlow(ndim)
    double precision :: plow(ndim),qlow(ndim),xhigh(ndim)
    double precision :: kx1(ndim),ku1(ndim),kv1(ndim),kp1(ndim),kq1(ndim)
    double precision :: dxb^D
    type(trace_sph_interp_ctx) :: ctx_old
    integer :: igrid_old,igrid_trial,point_domain,iter
    integer :: ctx_status,trial_status
    logical :: hit_ok,use_short_boundary

    if (.not.state%active) return
    use_short_boundary=.false.
    if (present(short_boundary)) use_short_boundary=short_boundary

    do while(state%active .and. state%igrid==igrid .and. &
         state%nstep<max_steps)
      xold=state%x
      uold=state%u
      vold=state%v
      pold=state%p
      qold=state%q
      ^D&dxb^D=rnode(rpdx^D_,state%igrid);
      ctx_status=trace_status_unsupported_geometry
      if (geo_coordinate==geo_spherical) then
        call trace_spherical_interp_ctx_build_cached(xold,state%igrid, &
             state%sph_cache,ctx_old,ctx_status)
      endif
      if (ctx_status==trace_status_active) then
        h=trace_spherical_effective_step_ctx(dL,ctx_old)
      else
        h=trace_effective_step(xold,dL,state%igrid,dxb^D)
      endif
      if (state%direction<0) h=-h

      if (state%has_extra) then
        if (ctx_status==trace_status_active) then
          call trace_tangent_rhs(xold,uold,vold,state%igrid,threshold,kx1, &
               ku1,kv1,state%status,pold,qold,kp1,kq1,sph_ctx=ctx_old)
        else
          call trace_tangent_rhs(xold,uold,vold,state%igrid,threshold,kx1, &
               ku1,kv1,state%status,pold,qold,kp1,kq1)
        endif
      else
        if (ctx_status==trace_status_active) then
          call trace_tangent_rhs(xold,uold,vold,state%igrid,threshold,kx1, &
               ku1,kv1,state%status,sph_ctx=ctx_old)
        else
          call trace_tangent_rhs(xold,uold,vold,state%igrid,threshold,kx1, &
               ku1,kv1,state%status)
        endif
      endif
      if (state%status/=trace_status_active) then
        state%active=.false.
        state%complete=.true.
        if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
        return
      endif
      if (use_short_boundary) then
        ! The q0/logQ short-boundary path keeps the same RK2 RHS but clamps
        ! the trial stage and bisects the final step to the domain boundary.
        hsign=sign(one,h)
        if (state%has_extra) then
          call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid, &
               h,threshold,kx1,ku1,kv1,xtrial,utrial,vtrial,state%status, &
               pold,qold,kp1,kq1,ptrial,qtrial,sph_cache=state%sph_cache, &
               clamp_stage=.true.)
        else
          call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid, &
               h,threshold,kx1,ku1,kv1,xtrial,utrial,vtrial,state%status, &
               sph_cache=state%sph_cache,clamp_stage=.true.)
        endif
        if (state%status/=trace_status_active) then
          state%active=.false.
          state%complete=.true.
          if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
          return
        endif

        point_domain=0
        {if (xtrial(^DB)>=xprobmin^DB .and. xtrial(^DB)<xprobmax^DB) point_domain=point_domain+1\}
        if (point_domain==ndim) then
          igrid_old=state%igrid
          call trace_locate_point_with_hint(xtrial,state%igrid,igrid_trial, &
               state%status)
          if (state%status/=trace_status_active) then
            state%active=.false.
            state%complete=.true.
            if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
            return
          endif
          call trace_tangent_accumulate_twist(state,xtrial,abs(h),threshold)
          state%x=xtrial
          state%u=utrial
          state%v=vtrial
          if (state%has_extra) then
            state%p=ptrial
            state%q=qtrial
          endif
          state%igrid=igrid_trial
          state%length=state%length+abs(h)
          state%nstep=state%nstep+1
          if (trace_rk2_stats_enabled) call trace_rk2_stats_note_step(abs(h),.false., &
               igrid_trial/=igrid_old)
          if (trace_spherical_profile_enabled .and. &
               geo_coordinate==geo_spherical) then
            call trace_spherical_profile_add_count(trace_profile_steps,1_8)
          endif
        else
          h_abs_low=zero
          h_abs_high=abs(h)
          xlow=xold
          ulow=uold
          vlow=vold
          plow=pold
          qlow=qold
          xhigh=xtrial
          do iter=1,max_bisect
            h_abs_mid=half*(h_abs_low+h_abs_high)
            if (h_abs_high-h_abs_low<=100.d0*epsilon(one)*max(one,abs(h))) &
                 exit
            trial_status=trace_status_active
            if (state%has_extra) then
              call trace_tangent_rk2_trial_from_rhs(xold,uold,vold, &
                   state%igrid,hsign*h_abs_mid,threshold,kx1,ku1,kv1, &
                   xtrial,utrial,vtrial,trial_status,pold,qold,kp1,kq1, &
                   ptrial,qtrial,sph_cache=state%sph_cache, &
                   clamp_stage=.true.)
            else
              call trace_tangent_rk2_trial_from_rhs(xold,uold,vold, &
                   state%igrid,hsign*h_abs_mid,threshold,kx1,ku1,kv1, &
                   xtrial,utrial,vtrial,trial_status, &
                   sph_cache=state%sph_cache,clamp_stage=.true.)
            endif
            if (trial_status/=trace_status_active) then
              h_abs_high=h_abs_mid
              cycle
            endif
            point_domain=0
            {if (xtrial(^DB)>=xprobmin^DB .and. xtrial(^DB)<xprobmax^DB) point_domain=point_domain+1\}
            if (point_domain==ndim) then
              h_abs_low=h_abs_mid
              xlow=xtrial
              ulow=utrial
              vlow=vtrial
              if (state%has_extra) then
                plow=ptrial
                qlow=qtrial
              endif
            else
              h_abs_high=h_abs_mid
              xhigh=xtrial
            endif
          enddo

          call trace_intersect_domain(xlow,xhigh,xhit,hit_ok,state%face_id, &
               alpha_hit)
          if (hit_ok) then
            h_abs_final=h_abs_low+max(zero,min(one,alpha_hit))* &
                 (h_abs_high-h_abs_low)
          else
            xhit=xlow
            h_abs_final=h_abs_low
            call trace_boundary_face_at_point(xhit,state%face_id,hit_ok)
          endif
          if (.not.hit_ok) then
            state%status=trace_status_out_of_domain
            state%active=.false.
            state%complete=.true.
            if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
            return
          endif

          partial_length=h_abs_final
          if (partial_length>100.d0*epsilon(one)*max(one,abs(h))) then
            trial_status=trace_status_active
            if (state%has_extra) then
              call trace_tangent_rk2_trial_from_rhs(xold,uold,vold, &
                   state%igrid,hsign*h_abs_final,threshold,kx1,ku1,kv1, &
                   xtrial,utrial,vtrial,trial_status,pold,qold,kp1,kq1, &
                   ptrial,qtrial,sph_cache=state%sph_cache, &
                   clamp_stage=.true.)
            else
              call trace_tangent_rk2_trial_from_rhs(xold,uold,vold, &
                   state%igrid,hsign*h_abs_final,threshold,kx1,ku1,kv1, &
                   xtrial,utrial,vtrial,trial_status, &
                   sph_cache=state%sph_cache,clamp_stage=.true.)
            endif
            if (trial_status/=trace_status_active) then
              state%status=trial_status
              state%active=.false.
              state%complete=.true.
              if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
              return
            endif
            call trace_tangent_accumulate_twist(state,xhit,partial_length, &
                 threshold)
            state%u=utrial
            state%v=vtrial
            if (state%has_extra) then
              state%p=ptrial
              state%q=qtrial
            endif
            state%length=state%length+partial_length
            state%nstep=state%nstep+1
            if (trace_rk2_stats_enabled) call trace_rk2_stats_note_step(partial_length,.true.,.false.)
            if (trace_spherical_profile_enabled .and. &
                 geo_coordinate==geo_spherical) then
              call trace_spherical_profile_add_count(trace_profile_steps,1_8)
            endif
          else
            state%u=ulow
            state%v=vlow
            if (state%has_extra) then
              state%p=plow
              state%q=qlow
            endif
          endif
          state%x=xhit
          state%endpoint=xhit
          state%status=trace_status_boundary
          state%active=.false.
          if (trace_spherical_profile_enabled .and. &
               geo_coordinate==geo_spherical) then
            call trace_spherical_profile_add_count( &
                 trace_profile_boundary_events,1_8)
          endif
          call trace_tangent_state_finalize_boundary(state,threshold)
          if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
          return
        endif
        cycle
      endif
      xprobe=xold+h*kx1

      point_domain=0
      {if (xprobe(^DB)>=xprobmin^DB .and. xprobe(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        call trace_intersect_domain(xold,xprobe,xhit,hit_ok,state%face_id, &
             alpha_hit)
        if (.not.hit_ok) then
          call trace_boundary_face_at_point(xold,state%face_id,hit_ok)
          state%status=trace_status_boundary
          state%active=.false.
          state%x=xold
          state%u=uold
          state%v=vold
          if (state%has_extra) then
            state%p=pold
            state%q=qold
          endif
          state%endpoint=xold
          if (hit_ok) call trace_tangent_state_finalize_boundary(state, &
               threshold)
          if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
          return
        endif
        partial_length=trace_segment_length(xold,xhit,h,alpha_hit)
        if (partial_length<=100.d0*epsilon(one)*max(one,abs(h))) then
          state%x=xhit
          state%u=uold
          state%v=vold
          if (state%has_extra) then
            state%p=pold
            state%q=qold
          endif
          state%endpoint=xhit
          state%status=trace_status_boundary
          state%active=.false.
          call trace_tangent_state_finalize_boundary(state,threshold)
          if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
          return
        endif
        h_partial=sign(partial_length,h)
        if (state%has_extra) then
          call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid, &
               h_partial,threshold,kx1,ku1,kv1,xtrial,utrial,vtrial, &
               state%status,pold,qold,kp1,kq1,ptrial,qtrial, &
               sph_cache=state%sph_cache)
        else
          call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid, &
               h_partial,threshold,kx1,ku1,kv1,xtrial,utrial,vtrial, &
               state%status,sph_cache=state%sph_cache)
        endif
        if (state%status/=trace_status_active) then
          state%active=.false.
          state%complete=.true.
          if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
          return
        endif
        call trace_tangent_accumulate_twist(state,xhit,partial_length, &
             threshold)
        state%x=xhit
        state%u=utrial
        state%v=vtrial
        if (state%has_extra) then
          state%p=ptrial
          state%q=qtrial
        endif
        state%length=state%length+partial_length
        state%endpoint=xhit
        state%status=trace_status_boundary
        state%active=.false.
        if (trace_rk2_stats_enabled) call trace_rk2_stats_note_step(partial_length,.true.,.false.)
        call trace_tangent_state_finalize_boundary(state,threshold)
        if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
        return
      endif

      if (state%has_extra) then
        call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid,h, &
             threshold,kx1,ku1,kv1,xtrial,utrial,vtrial,state%status, &
             pold,qold,kp1,kq1,ptrial,qtrial,sph_cache=state%sph_cache)
      else
        call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid,h, &
             threshold,kx1,ku1,kv1,xtrial,utrial,vtrial,state%status, &
             sph_cache=state%sph_cache)
      endif
      if (state%status/=trace_status_active) then
        state%active=.false.
        state%complete=.true.
        if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
        return
      endif

      point_domain=0
      {if (xtrial(^DB)>=xprobmin^DB .and. xtrial(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain==ndim) then
        igrid_old=state%igrid
        call trace_locate_point_with_hint(xtrial,state%igrid,igrid_trial, &
             state%status)
        if (state%status/=trace_status_active) then
          state%active=.false.
          state%complete=.true.
          if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
          return
        endif
        call trace_tangent_accumulate_twist(state,xtrial,abs(h),threshold)
        state%x=xtrial
        state%u=utrial
        state%v=vtrial
        if (state%has_extra) then
          state%p=ptrial
          state%q=qtrial
        endif
        state%igrid=igrid_trial
        state%length=state%length+abs(h)
        state%nstep=state%nstep+1
        if (trace_rk2_stats_enabled) call trace_rk2_stats_note_step(abs(h),.false., &
             igrid_trial/=igrid_old)
        if (trace_spherical_profile_enabled .and. &
             geo_coordinate==geo_spherical) then
          call trace_spherical_profile_add_count(trace_profile_steps,1_8)
        endif
      else
        call trace_intersect_domain(xold,xtrial,xhit,hit_ok,state%face_id, &
             alpha_hit)
        if (.not.hit_ok) then
          call trace_boundary_face_at_point(xold,state%face_id,hit_ok)
          state%status=trace_status_boundary
          state%active=.false.
          state%x=xold
          state%u=uold
          state%v=vold
          if (state%has_extra) then
            state%p=pold
            state%q=qold
          endif
          state%endpoint=xold
          if (hit_ok) call trace_tangent_state_finalize_boundary(state, &
               threshold)
          if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
          return
        endif
        partial_length=trace_segment_length(xold,xhit,h,alpha_hit)
        if (partial_length<=100.d0*epsilon(one)*max(one,abs(h))) then
          state%x=xhit
          state%u=uold
          state%v=vold
          if (state%has_extra) then
            state%p=pold
            state%q=qold
          endif
          state%endpoint=xhit
          state%status=trace_status_boundary
          state%active=.false.
          call trace_tangent_state_finalize_boundary(state,threshold)
          if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
          return
        endif
        h_partial=sign(partial_length,h)
        if (state%has_extra) then
          call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid, &
               h_partial,threshold,kx1,ku1,kv1,xtrial,utrial,vtrial, &
               state%status,pold,qold,kp1,kq1,ptrial,qtrial, &
               sph_cache=state%sph_cache)
        else
          call trace_tangent_rk2_trial_from_rhs(xold,uold,vold,state%igrid, &
               h_partial,threshold,kx1,ku1,kv1,xtrial,utrial,vtrial, &
               state%status,sph_cache=state%sph_cache)
        endif
        if (state%status/=trace_status_active) then
          state%active=.false.
          state%complete=.true.
          if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
          return
        endif
        call trace_tangent_accumulate_twist(state,xhit,partial_length, &
             threshold)
        state%x=xhit
        state%u=utrial
        state%v=vtrial
        if (state%has_extra) then
          state%p=ptrial
          state%q=qtrial
        endif
        state%length=state%length+partial_length
        state%endpoint=xhit
        state%status=trace_status_boundary
        state%active=.false.
        if (trace_rk2_stats_enabled) call trace_rk2_stats_note_step(partial_length,.true.,.false.)
        if (trace_spherical_profile_enabled .and. &
             geo_coordinate==geo_spherical) then
          call trace_spherical_profile_add_count(trace_profile_boundary_events, &
               1_8)
        endif
        call trace_tangent_state_finalize_boundary(state,threshold)
        if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
        return
      endif
    enddo

    if (state%active .and. state%nstep>=max_steps) then
      state%status=trace_status_max_steps
      state%active=.false.
      state%complete=.true.
      if (trace_rk2_stats_enabled) call trace_rk2_stats_note_completion(state%status,state%nstep)
    endif
    if ((.not.state%active) .and. trace_spherical_profile_enabled .and. &
         geo_coordinate==geo_spherical) then
      call trace_spherical_profile_note_trace_steps(state%nstep)
    endif
  end subroutine trace_tangent_state_advance_in_grid

  subroutine trace_tangent_state_finalize_boundary(state,threshold)
    type(trace_tangent_state), intent(inout) :: state
    double precision, intent(in) :: threshold

    if (state%status/=trace_status_boundary) return
    if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
         call trace_spherical_profile_add_count( &
         trace_profile_endpoint_finalizations,1_8)

    call trace_endpoint_B_bhat(state%endpoint,state%face_id,state%igrid, &
         threshold,state%endpoint_B,state%endpoint_bhat,state%status)
    if (state%status/=trace_status_boundary) return

    call trace_project_to_perp_bhat(state%u,state%endpoint_bhat, &
         state%u_perp)
    call trace_project_to_perp_bhat(state%v,state%endpoint_bhat, &
         state%v_perp)
    if (state%has_extra) then
      call trace_project_to_perp_bhat(state%p,state%endpoint_bhat, &
           state%p_perp)
      call trace_project_to_perp_bhat(state%q,state%endpoint_bhat, &
           state%q_perp)
    endif
    state%complete=.true.
  end subroutine trace_tangent_state_finalize_boundary

  subroutine trace_qperp_prepare_seed_result(seed,field_min,result,igrid, &
       status)
    double precision, intent(in) :: seed(ndim),field_min
    type(trace_qperp_result), intent(out) :: result
    integer, intent(out) :: igrid,status

    double precision :: B3(3),Bseed_norm
    integer :: indomain

    call trace_init_qperp_result(seed,result)
    igrid=-1
    status=trace_status_active

    indomain=0
    {if (seed(^DB)>=xprobmin^DB .and. seed(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain/=ndim) then
      status=trace_status_seed_outside
      result%status=status
      return
    endif

    call trace_debug_locate_point(seed,igrid,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif

    call sample_B_at_point(seed,igrid,field_min,B3,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif
    result%B_seed=B3(1:ndim)
    Bseed_norm=dsqrt(sum(result%B_seed**2))
    if (Bseed_norm<=zero .or. Bseed_norm<field_min) then
      status=trace_status_weak_field
      result%status=status
      return
    endif

    result%bhat_seed=result%B_seed/Bseed_norm
    call trace_make_perp_basis(B3/Bseed_norm,result%u0,result%v0,status)
    if (status/=trace_status_active) result%status=status
  end subroutine trace_qperp_prepare_seed_result

  subroutine trace_qperp_compute_scalars(result,field_min)
    type(trace_qperp_result), intent(inout) :: result
    double precision, intent(in) :: field_min

    double precision :: Bseed_norm,Bf_norm,Bb_norm,dot_f,dot_b,qtmp

    Bseed_norm=dsqrt(sum(result%B_seed**2))
    Bf_norm=dsqrt(sum(result%forward_B**2))
    Bb_norm=dsqrt(sum(result%backward_B**2))
    if (Bseed_norm<=zero .or. Bf_norm<=zero .or. Bb_norm<=zero .or. &
         Bseed_norm<field_min .or. Bf_norm<field_min .or. &
         Bb_norm<field_min) then
      result%status=trace_status_weak_field
      return
    endif

    dot_f=sum(result%u_forward_perp*result%v_forward_perp)
    dot_b=sum(result%u_backward_perp*result%v_backward_perp)
    result%N2=sum(result%u_forward_perp**2)* &
         sum(result%v_backward_perp**2)+ &
         sum(result%u_backward_perp**2)* &
         sum(result%v_forward_perp**2)-2.d0*dot_f*dot_b
    result%bfactor=abs(Bf_norm*Bb_norm)/(Bseed_norm**2)
    qtmp=result%N2*result%bfactor
    if (qtmp>zero .and. ieee_is_finite(qtmp)) then
      qtmp=max(qtmp,trace_q_min)
      result%qperp=qtmp
      result%logqperp=dlog10(qtmp)
      result%valid=.true.
      result%status=trace_status_boundary
      result%qperp0=qtmp
      result%logqperp0=result%logqperp
      result%valid_qperp0=.true.
      result%status_qperp0=trace_status_boundary
    else
      result%N2=trace_debug_nan()
      result%bfactor=trace_debug_nan()
      result%qperp=trace_debug_nan()
      result%logqperp=trace_debug_nan()
      result%status=trace_status_invalid_input
      result%valid=.false.
      result%qperp0=trace_debug_nan()
      result%logqperp0=trace_debug_nan()
      result%valid_qperp0=.false.
      result%status_qperp0=trace_status_invalid_input
    endif
    call trace_qperp_compute_q0_scalars(result,field_min)
  end subroutine trace_qperp_compute_scalars

  subroutine trace_qperp_compute_q0_scalars(result,field_min)
    type(trace_qperp_result), intent(inout) :: result
    double precision, intent(in) :: field_min

    double precision :: Bseed_norm,Bnf,Bnb,N2bnd,qtmp
    double precision :: uf_face(ndim),vf_face(ndim)
    double precision :: ub_face(ndim),vb_face(ndim)
    integer :: status_f,status_b

    result%q0=trace_debug_nan()
    result%logq0=trace_debug_nan()
    result%N2_qperp0=trace_debug_nan()
    result%bfactor_qperp0=trace_debug_nan()
    result%forward_Bn_q0=trace_debug_nan()
    result%backward_Bn_q0=trace_debug_nan()
    result%valid_q0=.false.

    if (result%forward_status/=trace_status_boundary) then
      result%status_q0=result%forward_status
      return
    endif
    if (result%backward_status/=trace_status_boundary) then
      result%status_q0=result%backward_status
      return
    endif

    Bseed_norm=dsqrt(sum(result%B_seed**2))
    if (Bseed_norm<=zero .or. Bseed_norm<field_min) then
      result%status_q0=trace_status_weak_field
      return
    endif

    call trace_qperp_project_to_boundary_face(result%u_forward_perp, &
         result%forward_B,result%forward_face,uf_face,Bnf,status_f)
    call trace_qperp_project_to_boundary_face(result%v_forward_perp, &
         result%forward_B,result%forward_face,vf_face,Bnf,status_f)
    call trace_qperp_project_to_boundary_face(result%u_backward_perp, &
         result%backward_B,result%backward_face,ub_face,Bnb,status_b)
    call trace_qperp_project_to_boundary_face(result%v_backward_perp, &
         result%backward_B,result%backward_face,vb_face,Bnb,status_b)
    if (status_f/=trace_status_active) then
      result%status_q0=status_f
      return
    endif
    if (status_b/=trace_status_active) then
      result%status_q0=status_b
      return
    endif
    if (abs(Bnf)<field_min .or. abs(Bnb)<field_min) then
      result%status_q0=trace_status_weak_field
      return
    endif

    N2bnd=sum(uf_face**2)*sum(vb_face**2)+ &
         sum(ub_face**2)*sum(vf_face**2)- &
         2.d0*sum(uf_face*vf_face)*sum(ub_face*vb_face)
    result%N2_qperp0=abs(N2bnd)
    result%bfactor_qperp0=abs(Bnf*Bnb)/(Bseed_norm**2)
    result%forward_Bn_q0=Bnf
    result%backward_Bn_q0=Bnb
    qtmp=result%N2_qperp0*result%bfactor_qperp0
    if (qtmp>zero .and. ieee_is_finite(qtmp)) then
      qtmp=max(qtmp,trace_q_min)
      result%q0=qtmp
      result%logq0=dlog10(qtmp)
      result%valid_q0=.true.
      result%status_q0=trace_status_boundary
    else
      result%q0=trace_debug_nan()
      result%logq0=trace_debug_nan()
      result%valid_q0=.false.
      result%status_q0=trace_status_invalid_input
    endif
  end subroutine trace_qperp_compute_q0_scalars

  subroutine trace_qperp_project_to_boundary_face(vec,B,face_id,vec_face, &
       Bn,status)
    double precision, intent(in) :: vec(ndim),B(ndim)
    integer, intent(in) :: face_id
    double precision, intent(out) :: vec_face(ndim),Bn
    integer, intent(out) :: status

    double precision :: normal(ndim),vecn
    logical :: ok

    vec_face=zero
    Bn=zero
    status=trace_status_invalid_input
    call trace_face_normal(face_id,normal,ok)
    if (.not.ok) return
    Bn=sum(B*normal)
    if (Bn==zero) then
      status=trace_status_weak_field
      return
    endif
    vecn=sum(vec*normal)
    vec_face=vec-vecn/Bn*B
    status=trace_status_active
  end subroutine trace_qperp_project_to_boundary_face

  subroutine trace_qperp_finalize_from_states(forward_state,backward_state, &
       result,field_min)
    type(trace_tangent_state), intent(in) :: forward_state,backward_state
    type(trace_qperp_result), intent(inout) :: result
    double precision, intent(in) :: field_min

    result%forward_endpoint=forward_state%endpoint
    result%forward_B=forward_state%endpoint_B
    result%forward_bhat=forward_state%endpoint_bhat
    result%forward_length=forward_state%length
    result%forward_nstep=forward_state%nstep
    result%forward_face=forward_state%face_id
    result%forward_status=forward_state%status
    result%u_forward_perp=forward_state%u_perp
    result%v_forward_perp=forward_state%v_perp

    result%backward_endpoint=backward_state%endpoint
    result%backward_B=backward_state%endpoint_B
    result%backward_bhat=backward_state%endpoint_bhat
    result%backward_length=backward_state%length
    result%backward_nstep=backward_state%nstep
    result%backward_face=backward_state%face_id
    result%backward_status=backward_state%status
    result%u_backward_perp=backward_state%u_perp
    result%v_backward_perp=backward_state%v_perp

    if (forward_state%status/=trace_status_boundary) then
      result%status=forward_state%status
      return
    endif
    if (backward_state%status/=trace_status_boundary) then
      result%status=backward_state%status
      return
    endif

    call trace_qperp_compute_scalars(result,field_min)
  end subroutine trace_qperp_finalize_from_states

  subroutine trace_q0_finalize_from_states(forward_state,backward_state, &
       result,field_min)
    type(trace_tangent_state), intent(in) :: forward_state,backward_state
    type(trace_qperp_result), intent(inout) :: result
    double precision, intent(in) :: field_min

    result%forward_endpoint=forward_state%endpoint
    result%forward_B=forward_state%endpoint_B
    result%forward_bhat=forward_state%endpoint_bhat
    result%forward_length=forward_state%length
    result%forward_nstep=forward_state%nstep
    result%forward_face=forward_state%face_id
    result%forward_status=forward_state%status
    result%u_forward_perp=forward_state%u_perp
    result%v_forward_perp=forward_state%v_perp

    result%backward_endpoint=backward_state%endpoint
    result%backward_B=backward_state%endpoint_B
    result%backward_bhat=backward_state%endpoint_bhat
    result%backward_length=backward_state%length
    result%backward_nstep=backward_state%nstep
    result%backward_face=backward_state%face_id
    result%backward_status=backward_state%status
    result%u_backward_perp=backward_state%u_perp
    result%v_backward_perp=backward_state%v_perp

    if (forward_state%status/=trace_status_boundary) then
      result%status=forward_state%status
      result%status_q0=forward_state%status
      return
    endif
    if (backward_state%status/=trace_status_boundary) then
      result%status=backward_state%status
      result%status_q0=backward_state%status
      return
    endif

    result%status=trace_status_boundary
    call trace_qperp_compute_q0_scalars(result,field_min)
  end subroutine trace_q0_finalize_from_states

  subroutine trace_spherical_rmin_q_prepare_seed_result(seed,field_min, &
       result,igrid,status)
    double precision, intent(in) :: seed(ndim),field_min
    type(trace_qperp_result), intent(out) :: result
    integer, intent(out) :: igrid,status

    double precision :: Bcart(3),bhat_cart(3),er(3),etheta(3),ephi(3)
    double precision :: radial_tol
    integer :: indomain

    call trace_init_qperp_result(seed,result)
    igrid=-1
    status=trace_status_active

    if (ndim/=3 .or. geo_coordinate/=geo_spherical) then
      status=trace_status_unsupported_geometry
      result%status=status
      result%status_q0=status
      return
    endif

    indomain=0
    {if (seed(^DB)>=xprobmin^DB .and. seed(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain/=ndim) then
      status=trace_status_seed_outside
      result%status=status
      result%status_q0=status
      return
    endif

    call trace_debug_locate_point(seed,igrid,status)
    if (status/=trace_status_active) then
      result%status=status
      result%status_q0=status
      return
    endif

    {^IFTHREED
    call trace_spherical_sample_B_bhat_cart(seed,igrid,field_min,Bcart, &
         bhat_cart,status)
    if (status/=trace_status_active) then
      result%status=status
      result%status_q0=status
      return
    endif
    call trace_spherical_basis(seed,er,etheta,ephi,status)
    if (status/=trace_status_active) then
      result%status=status
      result%status_q0=status
      return
    endif
    result%B_seed=Bcart(1:ndim)
    result%bhat_seed=bhat_cart(1:ndim)
    radial_tol=100.d0*epsilon(one)*max(one,abs(xprobmin1),abs(xprobmax1))
    ! Boundary seeds use the radial-surface basis; interior sampled Q uses a
    ! B-transverse basis to avoid a nearly singular initial tangent plane.
    if (abs(seed(1)-xprobmin1)<=radial_tol .or. &
         abs(seed(1)-xprobmax1)<=radial_tol) then
      result%u0=etheta(1:ndim)
      result%v0=ephi(1:ndim)
    else
      call trace_make_perp_basis(bhat_cart,result%u0,result%v0,status)
      if (status/=trace_status_active) then
        result%status=status
        result%status_q0=status
        return
      endif
    endif
    }
  end subroutine trace_spherical_rmin_q_prepare_seed_result

  subroutine trace_spherical_rmin_q_finalize_from_states(forward_state, &
       backward_state,result,field_min)
    type(trace_tangent_state), intent(in) :: forward_state,backward_state
    type(trace_qperp_result), intent(inout) :: result
    double precision, intent(in) :: field_min

    call trace_spherical_radial_q_finalize_from_states(forward_state, &
         backward_state,result,field_min)

  end subroutine trace_spherical_rmin_q_finalize_from_states

  subroutine trace_spherical_radial_q_finalize_from_states(forward_state, &
       backward_state,result,field_min)
    type(trace_tangent_state), intent(in) :: forward_state,backward_state
    type(trace_qperp_result), intent(inout) :: result
    double precision, intent(in) :: field_min

    double precision :: Bnf,Bnb,frob2,qtmp,logqtmp,bfactor
    double precision :: Af(2,2),Ab(2,2)
    integer :: status_f,status_b,status_q
    logical :: valid_q

    result%forward_endpoint=forward_state%endpoint
    result%forward_B=forward_state%endpoint_B
    result%forward_bhat=forward_state%endpoint_bhat
    result%forward_length=forward_state%length
    result%forward_nstep=forward_state%nstep
    result%forward_face=forward_state%face_id
    result%forward_status=forward_state%status
    result%u_forward_perp=forward_state%u
    result%v_forward_perp=forward_state%v

    result%backward_endpoint=backward_state%endpoint
    result%backward_B=backward_state%endpoint_B
    result%backward_bhat=backward_state%endpoint_bhat
    result%backward_length=backward_state%length
    result%backward_nstep=backward_state%nstep
    result%backward_face=backward_state%face_id
    result%backward_status=backward_state%status
    result%u_backward_perp=backward_state%u
    result%v_backward_perp=backward_state%v

    result%q0=trace_debug_nan()
    result%logq0=trace_debug_nan()
    result%N2_qperp0=trace_debug_nan()
    result%bfactor_qperp0=trace_debug_nan()
    result%forward_Bn_q0=trace_debug_nan()
    result%backward_Bn_q0=trace_debug_nan()
    result%valid_q0=.false.

    if (forward_state%status/=trace_status_boundary) then
      result%status_q0=forward_state%status
      result%status=forward_state%status
      return
    endif
    if (backward_state%status/=trace_status_boundary) then
      result%status_q0=backward_state%status
      result%status=backward_state%status
      return
    endif
    result%status=trace_status_boundary

    if (.not.trace_spherical_radial_q_face_pair_admitted( &
         forward_state%face_id,backward_state%face_id)) then
      result%status_q0=trace_status_invalid_input
      return
    endif

    if (dsqrt(sum(result%B_seed**2))<=zero .or. &
         dsqrt(sum(result%B_seed**2))<field_min) then
      result%status_q0=trace_status_weak_field
      return
    endif

    call trace_spherical_radial_endpoint_matrix(forward_state,Af,Bnf,status_f)
    result%forward_Bn_q0=Bnf
    if (status_f/=trace_status_active) then
      result%status_q0=status_f
      return
    endif

    call trace_spherical_radial_endpoint_matrix(backward_state,Ab,Bnb,status_b)
    result%backward_Bn_q0=Bnb
    if (status_b/=trace_status_active) then
      result%status_q0=status_b
      return
    endif

    if (abs(Bnf)<1.d-10 .or. abs(Bnb)<1.d-10) then
      result%status_q0=trace_status_weak_field
      return
    endif

    call trace_spherical_q_from_endpoint_matrices(Af,Ab,qtmp,logqtmp, &
         frob2,bfactor,valid_q,status_q)
    result%N2_qperp0=frob2
    result%bfactor_qperp0=bfactor
    if (.not.valid_q) then
      result%status_q0=status_q
      return
    endif

    result%q0=qtmp
    result%logq0=logqtmp
    result%valid_q0=.true.
    result%status_q0=trace_status_boundary
  end subroutine trace_spherical_radial_q_finalize_from_states

  logical function trace_spherical_radial_q_face_pair_admitted(face_forward, &
       face_backward) result(is_admitted)
    integer, intent(in) :: face_forward,face_backward

    is_admitted=.false.
    select case (face_forward)
    case (trace_face_xmin)
      is_admitted=(face_backward==trace_face_xmin .or. &
           face_backward==trace_face_xmax)
    case (trace_face_xmax)
      is_admitted=(face_backward==trace_face_xmin)
    end select

  end function trace_spherical_radial_q_face_pair_admitted

  subroutine trace_spherical_radial_endpoint_matrix(state,Amat,Bn,status)
    type(trace_tangent_state), intent(in) :: state
    double precision, intent(out) :: Amat(2,2),Bn
    integer, intent(out) :: status

    double precision :: u_face(ndim),v_face(ndim)
    double precision :: er(3),etheta(3),ephi(3)
    double precision :: Bn_u,Bn_v
    integer :: status_u,status_v,basis_status

    Amat=zero
    Bn=zero
    status=trace_status_active

    call trace_spherical_radial_project_to_surface(state%u,state%endpoint_B, &
         state%endpoint,state%face_id,u_face,Bn_u,status_u)
    if (status_u/=trace_status_active) then
      status=status_u
      return
    endif

    call trace_spherical_radial_project_to_surface(state%v,state%endpoint_B, &
         state%endpoint,state%face_id,v_face,Bn_v,status_v)
    if (status_v/=trace_status_active) then
      status=status_v
      return
    endif

    Bn=Bn_u

    {^IFTHREED
    call trace_spherical_basis(state%endpoint,er,etheta,ephi,basis_status)
    if (basis_status/=trace_status_active) then
      status=basis_status
      return
    endif
    Amat(1,1)=sum(u_face*etheta(1:ndim))
    Amat(2,1)=sum(u_face*ephi(1:ndim))
    Amat(1,2)=sum(v_face*etheta(1:ndim))
    Amat(2,2)=sum(v_face*ephi(1:ndim))
    }

  end subroutine trace_spherical_radial_endpoint_matrix

  subroutine trace_spherical_radial_project_to_surface(vec,B,x,face_id, &
       vec_face,Bn,status)
    double precision, intent(in) :: vec(ndim),B(ndim),x(ndim)
    integer, intent(in) :: face_id
    double precision, intent(out) :: vec_face(ndim),Bn
    integer, intent(out) :: status

    double precision :: er(3),etheta(3),ephi(3),normal(ndim),vecn
    double precision :: bhat(ndim),Bnorm

    vec_face=zero
    Bn=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical) return

    {^IFTHREED
    call trace_spherical_basis(x,er,etheta,ephi,status)
    if (status/=trace_status_active) return
    select case (face_id)
    case (trace_face_xmin)
      normal=-er(1:ndim)
    case (trace_face_xmax)
      normal=er(1:ndim)
    case default
      status=trace_status_invalid_input
      return
    end select
    Bnorm=dsqrt(sum(B**2))
    if (Bnorm<=smalldouble) then
      status=trace_status_weak_field
      return
    endif
    bhat=B/Bnorm
    Bn=sum(bhat*normal)
    if (abs(Bn)<=smalldouble) then
      status=trace_status_weak_field
      return
    endif
    vecn=sum(vec*normal)
    vec_face=vec-vecn/Bn*bhat
    status=trace_status_active
    }
  end subroutine trace_spherical_radial_project_to_surface

  subroutine trace_spherical_qperp_prepare_seed_result(seed,field_min, &
       result,igrid,status)
    double precision, intent(in) :: seed(ndim),field_min
    type(trace_qperp_result), intent(out) :: result
    integer, intent(out) :: igrid,status

    double precision :: Bcart(3),bhat_cart(3)
    integer :: indomain

    call trace_init_qperp_result(seed,result)
    igrid=-1
    status=trace_status_active

    if (ndim/=3 .or. geo_coordinate/=geo_spherical) then
      status=trace_status_unsupported_geometry
      result%status=status
      return
    endif

    indomain=0
    {if (seed(^DB)>=xprobmin^DB .and. seed(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain/=ndim) then
      status=trace_status_seed_outside
      result%status=status
      return
    endif

    call trace_debug_locate_point(seed,igrid,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif

    {^IFTHREED
    call trace_spherical_sample_B_bhat_cart(seed,igrid,field_min,Bcart, &
         bhat_cart,status)
    if (status/=trace_status_active) then
      result%status=status
      return
    endif
    result%B_seed=Bcart(1:ndim)
    result%bhat_seed=bhat_cart(1:ndim)
    call trace_make_perp_basis(bhat_cart,result%u0,result%v0,status)
    if (status/=trace_status_active) result%status=status
    }
  end subroutine trace_spherical_qperp_prepare_seed_result

  subroutine trace_spherical_qperp_finalize_from_states(forward_state, &
       backward_state,result,field_min)
    type(trace_tangent_state), intent(in) :: forward_state,backward_state
    type(trace_qperp_result), intent(inout) :: result
    double precision, intent(in) :: field_min

    result%forward_endpoint=forward_state%endpoint
    result%forward_B=forward_state%endpoint_B
    result%forward_bhat=forward_state%endpoint_bhat
    result%forward_length=forward_state%length
    result%forward_nstep=forward_state%nstep
    result%forward_face=forward_state%face_id
    result%forward_status=forward_state%status
    result%u_forward_perp=forward_state%u_perp
    result%v_forward_perp=forward_state%v_perp

    result%backward_endpoint=backward_state%endpoint
    result%backward_B=backward_state%endpoint_B
    result%backward_bhat=backward_state%endpoint_bhat
    result%backward_length=backward_state%length
    result%backward_nstep=backward_state%nstep
    result%backward_face=backward_state%face_id
    result%backward_status=backward_state%status
    result%u_backward_perp=backward_state%u_perp
    result%v_backward_perp=backward_state%v_perp

    if (forward_state%status/=trace_status_boundary) then
      result%status=forward_state%status
      result%status_qperp0=forward_state%status
      return
    endif
    if (backward_state%status/=trace_status_boundary) then
      result%status=backward_state%status
      result%status_qperp0=backward_state%status
      return
    endif

    call trace_spherical_qperp_compute_from_states(result,field_min)
  end subroutine trace_spherical_qperp_finalize_from_states

  subroutine trace_spherical_qperp_compute_from_states(result,field_min)
    type(trace_qperp_result), intent(inout) :: result
    double precision, intent(in) :: field_min

    double precision :: f1(ndim),f2(ndim),b1(ndim),b2(ndim)
    double precision :: forward_bhat3(3),backward_bhat3(3)
    double precision :: Af(2,2),Ab(2,2)
    double precision :: Bseed_norm,Bf_norm,Bb_norm,frob2,qtmp,logqtmp,bfactor
    integer :: status_f,status_b,status_q
    logical :: valid_q

    result%qperp=trace_debug_nan()
    result%logqperp=trace_debug_nan()
    result%N2=trace_debug_nan()
    result%bfactor=trace_debug_nan()
    result%qperp0=trace_debug_nan()
    result%logqperp0=trace_debug_nan()
    result%N2_qperp0=trace_debug_nan()
    result%bfactor_qperp0=trace_debug_nan()
    result%valid=.false.
    result%valid_qperp0=.false.

    if (ndim/=3) then
      result%status=trace_status_unsupported_geometry
      result%status_qperp0=trace_status_unsupported_geometry
      return
    endif

    Bseed_norm=dsqrt(sum(result%B_seed**2))
    Bf_norm=dsqrt(sum(result%forward_B**2))
    Bb_norm=dsqrt(sum(result%backward_B**2))
    if (Bseed_norm<=zero .or. Bf_norm<=zero .or. Bb_norm<=zero .or. &
         Bseed_norm<field_min .or. Bf_norm<field_min .or. &
         Bb_norm<field_min) then
      result%status=trace_status_weak_field
      result%status_qperp0=trace_status_weak_field
      return
    endif

    forward_bhat3=zero
    backward_bhat3=zero
    forward_bhat3(1:ndim)=result%forward_bhat
    backward_bhat3(1:ndim)=result%backward_bhat
    call trace_make_perp_basis(forward_bhat3,f1,f2,status_f)
    call trace_make_perp_basis(backward_bhat3,b1,b2,status_b)
    if (status_f/=trace_status_active) then
      result%status=status_f
      result%status_qperp0=status_f
      return
    endif
    if (status_b/=trace_status_active) then
      result%status=status_b
      result%status_qperp0=status_b
      return
    endif

    Af(1,1)=sum(result%u_forward_perp*f1)
    Af(2,1)=sum(result%u_forward_perp*f2)
    Af(1,2)=sum(result%v_forward_perp*f1)
    Af(2,2)=sum(result%v_forward_perp*f2)
    Ab(1,1)=sum(result%u_backward_perp*b1)
    Ab(2,1)=sum(result%u_backward_perp*b2)
    Ab(1,2)=sum(result%v_backward_perp*b1)
    Ab(2,2)=sum(result%v_backward_perp*b2)

    call trace_spherical_q_from_endpoint_matrices(Af,Ab,qtmp,logqtmp, &
         frob2,bfactor,valid_q,status_q)
    result%N2=frob2
    result%bfactor=bfactor
    result%N2_qperp0=frob2
    result%bfactor_qperp0=bfactor
    if (.not.valid_q) then
      result%status=status_q
      result%status_qperp0=status_q
      return
    endif

    result%qperp=qtmp
    result%logqperp=logqtmp
    result%valid=.true.
    result%status=trace_status_boundary
    result%qperp0=qtmp
    result%logqperp0=result%logqperp
    result%valid_qperp0=.true.
    result%status_qperp0=trace_status_boundary
  end subroutine trace_spherical_qperp_compute_from_states

  subroutine trace_spherical_q_from_endpoint_matrices(Af,Ab,qval,logq, &
       frob2,bfactor,valid,status)
    double precision, intent(in) :: Af(2,2),Ab(2,2)
    double precision, intent(out) :: qval,logq,frob2,bfactor
    logical, intent(out) :: valid
    integer, intent(out) :: status

    double precision :: Dmat(2,2),detAb,detD,q_tol,det_tol

    qval=trace_debug_nan()
    logq=trace_debug_nan()
    frob2=trace_debug_nan()
    bfactor=trace_debug_nan()
    valid=.false.
    status=trace_status_invalid_input

    detAb=Ab(1,1)*Ab(2,2)-Ab(1,2)*Ab(2,1)
    det_tol=1.d-12
    if (.not.ieee_is_finite(detAb) .or. abs(detAb)<=det_tol) then
      status=trace_status_singular_q
      return
    endif

    Dmat(1,1)=(Af(1,1)*Ab(2,2)-Af(1,2)*Ab(2,1))/detAb
    Dmat(1,2)=(-Af(1,1)*Ab(1,2)+Af(1,2)*Ab(1,1))/detAb
    Dmat(2,1)=(Af(2,1)*Ab(2,2)-Af(2,2)*Ab(2,1))/detAb
    Dmat(2,2)=(-Af(2,1)*Ab(1,2)+Af(2,2)*Ab(1,1))/detAb
    detD=Dmat(1,1)*Dmat(2,2)-Dmat(1,2)*Dmat(2,1)
    frob2=sum(Dmat**2)
    bfactor=abs(detD)
    if (.not.ieee_is_finite(detD) .or. abs(detD)<=det_tol) then
      status=trace_status_singular_q
      return
    endif

    qval=frob2/abs(detD)
    q_tol=1.d-8
    if (qval>zero .and. ieee_is_finite(qval)) then
      if (qval<trace_q_min-q_tol) then
        status=trace_status_bad_q_bound
        return
      endif
      qval=max(qval,trace_q_min)
      logq=dlog10(qval)
      valid=.true.
      status=trace_status_boundary
    else
      qval=trace_debug_nan()
      logq=trace_debug_nan()
      valid=.false.
      status=trace_status_invalid_input
    endif

  end subroutine trace_spherical_q_from_endpoint_matrices

  subroutine trace_spherical_basis(x,er,etheta,ephi,status)
    double precision, intent(in) :: x(ndim)
    double precision, intent(out) :: er(3),etheta(3),ephi(3)
    integer, intent(out) :: status

    double precision :: theta,phi,sint,cost,sinp,cosp

    er=zero
    etheta=zero
    ephi=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical) return

    {^IFTHREED
    if (.not.trace_spherical_metric_ok(x)) then
      status=trace_status_invalid_input
      return
    endif
    theta=x(2)
    phi=x(3)
    sint=dsin(theta)
    cost=dcos(theta)
    sinp=dsin(phi)
    cosp=dcos(phi)
    er=(/ sint*cosp,sint*sinp,cost /)
    etheta=(/ cost*cosp,cost*sinp,-sint /)
    ephi=(/ -sinp,cosp,zero /)
    status=trace_status_active
    }
  end subroutine trace_spherical_basis

  subroutine trace_spherical_coord_to_cart(x,xcart,status)
    double precision, intent(in) :: x(ndim)
    double precision, intent(out) :: xcart(3)
    integer, intent(out) :: status

    double precision :: r,theta,phi,sint

    xcart=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical) return

    {^IFTHREED
    if (.not.trace_spherical_metric_ok(x)) then
      status=trace_status_invalid_input
      return
    endif
    r=x(1)
    theta=x(2)
    phi=x(3)
    sint=dsin(theta)
    xcart(1)=r*sint*dcos(phi)
    xcart(2)=r*sint*dsin(phi)
    xcart(3)=r*dcos(theta)
    status=trace_status_active
    }
  end subroutine trace_spherical_coord_to_cart

  subroutine trace_cart_to_spherical_coord(xcart,x,status)
    double precision, intent(in) :: xcart(3)
    double precision, intent(out) :: x(ndim)
    integer, intent(out) :: status

    double precision :: r,rho,phi,twopi

    x=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical) return

    {^IFTHREED
    r=dsqrt(sum(xcart**2))
    if (r<=100.d0*epsilon(one)) then
      status=trace_status_invalid_input
      return
    endif
    rho=dsqrt(xcart(1)**2+xcart(2)**2)
    phi=datan2(xcart(2),xcart(1))
    twopi=2.d0*dacos(-one)
    if (phi<xprobmin3 .and. phi+twopi<xprobmax3) phi=phi+twopi
    if (phi>=xprobmax3 .and. phi-twopi>=xprobmin3) phi=phi-twopi
    x(1)=r
    x(2)=datan2(rho,xcart(3))
    x(3)=phi
    status=trace_status_active
    }
  end subroutine trace_cart_to_spherical_coord

  subroutine trace_spherical_interp_ctx_build(x,igrid,ctx,status)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: igrid
    type(trace_sph_interp_ctx), intent(out) :: ctx
    integer, intent(out) :: status

    double precision :: xd^D,dxc^D
    double precision :: t0
    integer :: ixI^L,ixbl^D

    if (trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_count(trace_profile_context_requests, &
           1_8)
      call trace_spherical_profile_add_count( &
           trace_profile_full_context_builds,1_8)
      t0=trace_spherical_profile_time()
    endif

    ctx%valid=.false.
    ctx%igrid=-1
    ctx%ixbl=0
    ctx%x=zero
    ctx%xd=zero
    ctx%dxc=zero
    ctx%w=zero
    ctx%dxloc=zero
    ctx%r=zero
    ctx%theta=zero
    ctx%phi=zero
    ctx%sin_theta=zero
    ctx%cos_theta=zero
    ctx%sin_phi=zero
    ctx%cos_phi=zero
    ctx%rsin_theta=zero
    ctx%h_local=zero
    ctx%bcorner_valid=.false.
    ctx%bcorner=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical .or. igrid<0) then
      if (trace_spherical_profile_enabled) then
        call trace_spherical_profile_add_count(trace_profile_context_failures, &
             1_8)
        call trace_spherical_profile_add_time(trace_profile_time_context, &
             trace_spherical_profile_time()-t0)
      endif
      return
    endif

    {^IFTHREED
    status=trace_status_out_of_domain
    if (.not.trace_spherical_metric_ok(x)) then
      status=trace_status_invalid_input
      if (trace_spherical_profile_enabled) then
        call trace_spherical_profile_add_count(trace_profile_context_failures, &
             1_8)
        call trace_spherical_profile_add_time(trace_profile_time_context, &
             trace_spherical_profile_time()-t0)
      endif
      return
    endif

    ixI^L=ixG^LL;
    call trace_interp_weights_block(x,igrid,ixI^L,ixbl^D,xd^D,dxc^D,status)
    if (status/=trace_status_active) then
      if (trace_spherical_profile_enabled) then
        call trace_spherical_profile_add_count(trace_profile_context_failures, &
             1_8)
        call trace_spherical_profile_add_time(trace_profile_time_context, &
             trace_spherical_profile_time()-t0)
      endif
      return
    endif

    call trace_spherical_interp_ctx_fill(x,igrid,ixbl^D,xd^D,dxc^D,ctx, &
         status)
    if (status==trace_status_active) then
      call trace_spherical_interp_ctx_load_Bcorners(ctx,status)
    endif
    if (trace_spherical_profile_enabled) then
      if (status/=trace_status_active) then
        call trace_spherical_profile_add_count(trace_profile_context_failures, &
             1_8)
      endif
      call trace_spherical_profile_add_time(trace_profile_time_context, &
           trace_spherical_profile_time()-t0)
    endif
    }
  end subroutine trace_spherical_interp_ctx_build

  subroutine trace_spherical_interp_ctx_build_cached(x,igrid,cache,ctx,status)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: igrid
    type(trace_sph_interp_ctx), intent(inout) :: cache
    type(trace_sph_interp_ctx), intent(out) :: ctx
    integer, intent(out) :: status

    double precision :: xd^D,dxc^D
    double precision :: t0
    integer :: ixI^L,ixbl^D
    logical :: same_cell

    if (trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_count(trace_profile_context_requests, &
           1_8)
      t0=trace_spherical_profile_time()
    endif

    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical .or. igrid<0) then
      if (trace_spherical_profile_enabled) then
        call trace_spherical_profile_add_count(trace_profile_context_failures, &
             1_8)
        call trace_spherical_profile_add_time(trace_profile_time_context, &
             trace_spherical_profile_time()-t0)
      endif
      call trace_spherical_interp_ctx_build(x,igrid,ctx,status)
      return
    endif

    {^IFTHREED
    if (cache%valid .and. cache%igrid==igrid) then
      ixI^L=ixG^LL;
      call trace_interp_weights_block_near(x,igrid,ixI^L,cache%ixbl, &
           ixbl^D,xd^D,dxc^D,status)
      if (status==trace_status_active) then
        call trace_spherical_interp_ctx_fill(x,igrid,ixbl^D,xd^D,dxc^D, &
             ctx,status)
        if (status==trace_status_active) then
          same_cell=(ixbl1==cache%ixbl(1) .and. ixbl2==cache%ixbl(2) .and. &
               ixbl3==cache%ixbl(3))
          if (same_cell .and. cache%bcorner_valid) then
            ctx%bcorner_valid=.true.
            ctx%bcorner=cache%bcorner
            if (trace_spherical_profile_enabled) then
              call trace_spherical_profile_add_count( &
                   trace_profile_b_corner_hits,1_8)
            endif
          else
            call trace_spherical_interp_ctx_load_Bcorners(ctx,status)
          endif
        endif
        if (status==trace_status_active) then
          cache=ctx
          if (trace_spherical_profile_enabled) then
            if (same_cell) then
              call trace_spherical_profile_add_count( &
                   trace_profile_same_cell_hits,1_8)
            else
              call trace_spherical_profile_add_count( &
                   trace_profile_same_grid_hits,1_8)
            endif
            call trace_spherical_profile_add_time(trace_profile_time_context, &
                 trace_spherical_profile_time()-t0)
          endif
        else if (trace_spherical_profile_enabled) then
          call trace_spherical_profile_add_count(trace_profile_context_failures, &
               1_8)
          call trace_spherical_profile_add_time(trace_profile_time_context, &
               trace_spherical_profile_time()-t0)
        endif
        return
      endif
    else if (cache%valid .and. trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_count(trace_profile_cache_invalidations, &
           1_8)
    endif

    if (trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_time(trace_profile_time_context, &
           trace_spherical_profile_time()-t0)
    endif
    call trace_spherical_interp_ctx_build(x,igrid,ctx,status)
    if (status==trace_status_active) cache=ctx
    }
  end subroutine trace_spherical_interp_ctx_build_cached

  subroutine trace_spherical_interp_ctx_fill(x,igrid,ixbl^D,xd^D,dxc^D,ctx, &
       status)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: igrid,ixbl^D
    double precision, intent(in) :: xd^D,dxc^D
    type(trace_sph_interp_ctx), intent(out) :: ctx
    integer, intent(out) :: status

    if (trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_count(trace_profile_hlocal_evals,1_8)
    endif

    ctx%valid=.false.
    ctx%igrid=-1
    ctx%ixbl=0
    ctx%x=zero
    ctx%xd=zero
    ctx%dxc=zero
    ctx%w=zero
    ctx%dxloc=zero
    ctx%r=zero
    ctx%theta=zero
    ctx%phi=zero
    ctx%sin_theta=zero
    ctx%cos_theta=zero
    ctx%sin_phi=zero
    ctx%cos_phi=zero
    ctx%rsin_theta=zero
    ctx%h_local=zero
    ctx%bcorner_valid=.false.
    ctx%bcorner=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical .or. igrid<0) return

    {^IFTHREED
    ctx%igrid=igrid
    ctx%x(1:3)=x(1:3)
    ctx%ixbl=(/ixbl1,ixbl2,ixbl3/)
    ctx%xd=(/xd1,xd2,xd3/)
    ctx%dxc=(/dxc1,dxc2,dxc3/)
    ctx%w(0,1)=one-xd1
    ctx%w(1,1)=xd1
    ctx%w(0,2)=one-xd2
    ctx%w(1,2)=xd2
    ctx%w(0,3)=one-xd3
    ctx%w(1,3)=xd3

    ctx%dxloc(1)=min(abs(ps(igrid)%dx(ixbl1,ixbl2,ixbl3,1)), &
         abs(ps(igrid)%dx(ixbl1+1,ixbl2,ixbl3,1)))
    ctx%dxloc(2)=min(abs(ps(igrid)%dx(ixbl1,ixbl2,ixbl3,2)), &
         abs(ps(igrid)%dx(ixbl1,ixbl2+1,ixbl3,2)))
    ctx%dxloc(3)=min(abs(ps(igrid)%dx(ixbl1,ixbl2,ixbl3,3)), &
         abs(ps(igrid)%dx(ixbl1,ixbl2,ixbl3+1,3)))
    ctx%dxloc=max(ctx%dxloc,smalldouble)

    ctx%r=x(1)
    ctx%theta=x(2)
    ctx%phi=x(3)
    ctx%sin_theta=dsin(ctx%theta)
    ctx%cos_theta=dcos(ctx%theta)
    ctx%sin_phi=dsin(ctx%phi)
    ctx%cos_phi=dcos(ctx%phi)
    ctx%rsin_theta=ctx%r*ctx%sin_theta
    if (ctx%r<=smalldouble .or. abs(ctx%sin_theta)<=smalldouble) then
      status=trace_status_invalid_input
      return
    endif

    ctx%h_local=trace_spherical_cell_scale_from_widths(ctx%r, &
         ctx%sin_theta,ctx%dxloc)
    ctx%h_local=max(ctx%h_local,smalldouble)
    ctx%valid=.true.
    status=trace_status_active
    }
  end subroutine trace_spherical_interp_ctx_fill

  subroutine trace_spherical_interp_ctx_load_Bcorners(ctx,status)
    type(trace_sph_interp_ctx), intent(inout) :: ctx
    integer, intent(out) :: status

    integer :: i1,i2,i3,igrid,ix1,ix2,ix3

    status=trace_status_out_of_domain
    if (.not.ctx%valid) return
    if (ctx%bcorner_valid) then
      if (trace_spherical_profile_enabled) then
        call trace_spherical_profile_add_count(trace_profile_b_corner_hits, &
             1_8)
      endif
      status=trace_status_active
      return
    endif
    if (.not.B0field .and. .not.allocated(iw_mag)) return

    {^IFTHREED
    igrid=ctx%igrid
    ctx%bcorner=zero
    do i3=0,1
      do i2=0,1
        do i1=0,1
          ix1=ctx%ixbl(1)+i1
          ix2=ctx%ixbl(2)+i2
          ix3=ctx%ixbl(3)+i3
          if (allocated(iw_mag)) then
            ctx%bcorner(i1,i2,i3,1:3)= &
                 ps(igrid)%w(ix1,ix2,ix3,iw_mag(1:3))
          endif
          if (B0field) then
            ctx%bcorner(i1,i2,i3,1:3)= &
                 ctx%bcorner(i1,i2,i3,1:3) &
                 +ps(igrid)%B0(ix1,ix2,ix3,1:3,0)
          endif
        enddo
      enddo
    enddo
    ctx%bcorner_valid=.true.
    if (trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_count(trace_profile_b_corner_loads, &
           1_8)
    endif
    status=trace_status_active
    }
  end subroutine trace_spherical_interp_ctx_load_Bcorners

  double precision function trace_spherical_effective_step_ctx(ds,ctx) &
       result(ds_eff)
    double precision, intent(in) :: ds
    type(trace_sph_interp_ctx), intent(in) :: ctx

    double precision :: step_cap,t0

    if (trace_spherical_profile_enabled) t0=trace_spherical_profile_time()

    ds_eff=ds
    if (.not.ctx%valid) then
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_step,trace_spherical_profile_time()-t0)
      return
    endif
    step_cap=abs(ds)
    if (ctx%h_local<=zero) then
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_step,trace_spherical_profile_time()-t0)
      return
    endif
    select case (trace_step_control_mode)
    case (trace_step_control_cell_fraction)
      ds_eff=min(step_cap,trace_step_fraction*ctx%h_local)
      if (trace_step_min>zero) ds_eff=max(ds_eff,min(trace_step_min,step_cap))
    case default
      ds_eff=min(step_cap,ctx%h_local)
    end select
    if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
         trace_profile_time_step,trace_spherical_profile_time()-t0)
  end function trace_spherical_effective_step_ctx

  subroutine trace_spherical_sample_Bsph_ctx(ctx,threshold,Bsph,status)
    type(trace_sph_interp_ctx), intent(in) :: ctx
    double precision, intent(in) :: threshold
    double precision, intent(out) :: Bsph(3)
    integer, intent(out) :: status

    double precision :: B2,weight
    double precision :: t0
    integer :: i1,i2,i3,j

    if (trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_count(trace_profile_b_samples,1_8)
      t0=trace_spherical_profile_time()
    endif

    Bsph=zero
    status=trace_status_out_of_domain
    if (.not.ctx%valid) then
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_b_interp,trace_spherical_profile_time()-t0)
      return
    endif
    if (.not.ctx%bcorner_valid) then
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_b_interp,trace_spherical_profile_time()-t0)
      return
    endif

    {^IFTHREED
    do i3=0,1
      do i2=0,1
        do i1=0,1
          weight=ctx%w(i1,1)*ctx%w(i2,2)*ctx%w(i3,3)
          do j=1,3
            Bsph(j)=Bsph(j)+ctx%bcorner(i1,i2,i3,j)*weight
          enddo
        enddo
      enddo
    enddo
    B2=sum(Bsph**2)
    if (B2<=zero .or. dsqrt(B2)<threshold) then
      status=trace_status_weak_field
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_b_interp,trace_spherical_profile_time()-t0)
      return
    endif
    status=trace_status_active
    if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
         trace_profile_time_b_interp,trace_spherical_profile_time()-t0)
    }
  end subroutine trace_spherical_sample_Bsph_ctx

  subroutine trace_spherical_sample_B_bhat_cart_ctx(ctx,threshold,Bcart, &
       bhat_cart,status)
    type(trace_sph_interp_ctx), intent(in) :: ctx
    double precision, intent(in) :: threshold
    double precision, intent(out) :: Bcart(3),bhat_cart(3)
    integer, intent(out) :: status

    double precision :: Bsph(3),Bnorm

    Bcart=zero
    bhat_cart=zero
    call trace_spherical_sample_Bsph_ctx(ctx,threshold,Bsph,status)
    if (status/=trace_status_active) return

    Bcart(1)=Bsph(1)*ctx%sin_theta*ctx%cos_phi &
         +Bsph(2)*ctx%cos_theta*ctx%cos_phi &
         -Bsph(3)*ctx%sin_phi
    Bcart(2)=Bsph(1)*ctx%sin_theta*ctx%sin_phi &
         +Bsph(2)*ctx%cos_theta*ctx%sin_phi &
         +Bsph(3)*ctx%cos_phi
    Bcart(3)=Bsph(1)*ctx%cos_theta-Bsph(2)*ctx%sin_theta
    Bnorm=dsqrt(sum(Bcart**2))
    if (Bnorm<=zero .or. Bnorm<threshold) then
      status=trace_status_weak_field
      return
    endif
    bhat_cart=Bcart/Bnorm
    status=trace_status_active
  end subroutine trace_spherical_sample_B_bhat_cart_ctx

  subroutine trace_spherical_sample_cached_curlB_ctx(ctx,curlB,status)
    type(trace_sph_interp_ctx), intent(in) :: ctx
    double precision, intent(out) :: curlB(3)
    integer, intent(out) :: status

    double precision :: weight
    double precision :: t0
    integer :: i1,i2,i3,j,igrid,ix1,ix2,ix3

    if (trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_count(trace_profile_curl_samples,1_8)
      t0=trace_spherical_profile_time()
    endif

    curlB=zero
    status=trace_status_bad_curl_stencil
    if (.not.ctx%valid) then
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_curl_interp,trace_spherical_profile_time()-t0)
      return
    endif
    if (.not.trace_spherical_curl_cache_ready) then
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_curl_interp,trace_spherical_profile_time()-t0)
      return
    endif
    if (.not.allocated(trace_spherical_curl_cache)) then
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_curl_interp,trace_spherical_profile_time()-t0)
      return
    endif
    igrid=ctx%igrid
    if (igrid>size(trace_spherical_curl_cache)) then
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_curl_interp,trace_spherical_profile_time()-t0)
      return
    endif
    if (.not.trace_spherical_curl_cache(igrid)%ready) then
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_curl_interp,trace_spherical_profile_time()-t0)
      return
    endif
    if (.not.allocated(trace_spherical_curl_cache(igrid)%current)) then
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_curl_interp,trace_spherical_profile_time()-t0)
      return
    endif

    {^IFTHREED
    do i3=0,1
      do i2=0,1
        do i1=0,1
          ix1=ctx%ixbl(1)+i1
          ix2=ctx%ixbl(2)+i2
          ix3=ctx%ixbl(3)+i3
          weight=ctx%w(i1,1)*ctx%w(i2,2)*ctx%w(i3,3)
          do j=1,3
            curlB(j)=curlB(j)+ &
                 trace_spherical_curl_cache(igrid)%current(ix1,ix2,ix3,j) &
                 *weight
          enddo
        enddo
      enddo
    enddo
    status=trace_status_active
    if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
         trace_profile_time_curl_interp,trace_spherical_profile_time()-t0)
    }
  end subroutine trace_spherical_sample_cached_curlB_ctx

  subroutine trace_spherical_sample_bhat_gradbhat_covariant_ctx(ctx, &
       threshold,bhat_cart,grad_bhat_cart,status)
    type(trace_sph_interp_ctx), intent(in) :: ctx
    double precision, intent(in) :: threshold
    double precision, intent(out) :: bhat_cart(3),grad_bhat_cart(3,3)
    integer, intent(out) :: status

    double precision :: dw(0:1,3),weight,Bnorm,projected
    double precision :: Bsph(3),bhat_sph(3),dBdq(3,3),dbhatdq(3,3)
    double precision :: A_local(3,3),E(3,3)
    double precision :: t0
    integer :: i1,i2,i3,j

    if (trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_count(trace_profile_grad_samples,1_8)
      t0=trace_spherical_profile_time()
    endif

    bhat_cart=zero
    grad_bhat_cart=zero
    status=trace_status_bad_grad_stencil
    if (.not.ctx%valid) then
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_grad_bhat,trace_spherical_profile_time()-t0)
      return
    endif
    if (.not.ctx%bcorner_valid) then
      status=trace_status_out_of_domain
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_grad_bhat,trace_spherical_profile_time()-t0)
      return
    endif

    {^IFTHREED
    dw(0,1)=-one/ctx%dxc(1)
    dw(1,1)= one/ctx%dxc(1)
    dw(0,2)=-one/ctx%dxc(2)
    dw(1,2)= one/ctx%dxc(2)
    dw(0,3)=-one/ctx%dxc(3)
    dw(1,3)= one/ctx%dxc(3)

    Bsph=zero
    dBdq=zero
    do i3=0,1
      do i2=0,1
        do i1=0,1
          weight=ctx%w(i1,1)*ctx%w(i2,2)*ctx%w(i3,3)
          Bsph(1:3)=Bsph(1:3)+ctx%bcorner(i1,i2,i3,1:3)*weight
          dBdq(1:3,1)=dBdq(1:3,1)+ctx%bcorner(i1,i2,i3,1:3)* &
               dw(i1,1)*ctx%w(i2,2)*ctx%w(i3,3)
          dBdq(1:3,2)=dBdq(1:3,2)+ctx%bcorner(i1,i2,i3,1:3)* &
               ctx%w(i1,1)*dw(i2,2)*ctx%w(i3,3)
          dBdq(1:3,3)=dBdq(1:3,3)+ctx%bcorner(i1,i2,i3,1:3)* &
               ctx%w(i1,1)*ctx%w(i2,2)*dw(i3,3)
        enddo
      enddo
    enddo

    Bnorm=dsqrt(sum(Bsph**2))
    if (Bnorm<=zero .or. Bnorm<threshold) then
      status=trace_status_weak_field
      if (trace_spherical_profile_enabled) call trace_spherical_profile_add_time( &
           trace_profile_time_grad_bhat,trace_spherical_profile_time()-t0)
      return
    endif

    bhat_sph=Bsph/Bnorm
    do j=1,3
      projected=sum(bhat_sph*dBdq(1:3,j))
      dbhatdq(1:3,j)=(dBdq(1:3,j)-bhat_sph(1:3)*projected)/Bnorm
    enddo

    E(1:3,1)=(/ctx%sin_theta*ctx%cos_phi, &
         ctx%sin_theta*ctx%sin_phi,ctx%cos_theta/)
    E(1:3,2)=(/ctx%cos_theta*ctx%cos_phi, &
         ctx%cos_theta*ctx%sin_phi,-ctx%sin_theta/)
    E(1:3,3)=(/-ctx%sin_phi,ctx%cos_phi,zero/)
    bhat_cart=bhat_sph(1)*E(1:3,1)+bhat_sph(2)*E(1:3,2) &
         +bhat_sph(3)*E(1:3,3)

    A_local(1:3,1)=dbhatdq(1:3,1)
    A_local(1,2)=(dbhatdq(1,2)-bhat_sph(2))/ctx%r
    A_local(2,2)=(dbhatdq(2,2)+bhat_sph(1))/ctx%r
    A_local(3,2)= dbhatdq(3,2)/ctx%r
    A_local(1,3)=(dbhatdq(1,3)-ctx%sin_theta*bhat_sph(3))/ctx%rsin_theta
    A_local(2,3)=(dbhatdq(2,3)-ctx%cos_theta*bhat_sph(3))/ctx%rsin_theta
    A_local(3,3)=(dbhatdq(3,3)+ctx%sin_theta*bhat_sph(1) &
         +ctx%cos_theta*bhat_sph(2))/ctx%rsin_theta

    grad_bhat_cart=matmul(E,matmul(A_local,transpose(E)))
    status=trace_status_active
    if (trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_count(trace_profile_grad_success,1_8)
      call trace_spherical_profile_add_time(trace_profile_time_grad_bhat, &
           trace_spherical_profile_time()-t0)
    endif
    }
  end subroutine trace_spherical_sample_bhat_gradbhat_covariant_ctx

  subroutine get_K_spherical_ctx(ctx,K,ftype,b_min,field_ok)
    type(trace_sph_interp_ctx), intent(in) :: ctx
    double precision, intent(out) :: K(ndim)
    character(len=std_len), intent(in) :: ftype
    double precision, intent(in), optional :: b_min
    logical, intent(out), optional :: field_ok

    double precision :: B(3),Ftotal,field_min
    logical :: valid_field
    integer :: status

    K=zero
    valid_field=.false.
    if (ftype/='Bfield') then
      if (present(field_ok)) field_ok=.false.
      return
    endif
    call trace_spherical_sample_Bsph_ctx(ctx,zero,B,status)
    if (status/=trace_status_active) then
      if (present(field_ok)) field_ok=.false.
      return
    endif
    Ftotal=dsqrt(sum(B**2))
    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    valid_field=Ftotal>zero .and. Ftotal>=field_min
    if (present(field_ok)) field_ok=valid_field
    if (.not.valid_field) return

    {^IFTHREED
    K(1)=B(1)/Ftotal
    K(2)=B(2)/(ctx%r*Ftotal)
    K(3)=B(3)/(ctx%rsin_theta*Ftotal)
    }
  end subroutine get_K_spherical_ctx

  subroutine trace_spherical_sample_B_bhat_cart(x,igrid,threshold,Bcart, &
       bhat_cart,status)
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: Bcart(3),bhat_cart(3)
    integer, intent(out) :: status

    type(trace_sph_interp_ctx) :: ctx

    Bcart=zero
    bhat_cart=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical) return

    {^IFTHREED
    call trace_spherical_interp_ctx_build(x,igrid,ctx,status)
    if (status/=trace_status_active) return
    call trace_spherical_sample_B_bhat_cart_ctx(ctx,threshold,Bcart, &
         bhat_cart,status)
    }
  end subroutine trace_spherical_sample_B_bhat_cart

  subroutine trace_spherical_sample_bhat_gradbhat_covariant(x,igrid, &
       threshold,bhat_cart,grad_bhat_cart,status)
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: bhat_cart(3),grad_bhat_cart(3,3)
    integer, intent(out) :: status

    type(trace_sph_interp_ctx) :: ctx

    bhat_cart=zero
    grad_bhat_cart=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical) return

    {^IFTHREED
    call trace_spherical_interp_ctx_build(x,igrid,ctx,status)
    if (status/=trace_status_active) return
    call trace_spherical_sample_bhat_gradbhat_covariant_ctx(ctx,threshold, &
         bhat_cart,grad_bhat_cart,status)
    }
  end subroutine trace_spherical_sample_bhat_gradbhat_covariant

  subroutine trace_spherical_sample_bhat_gradbhat_cartfd(x,igrid,threshold, &
       bhat,grad_bhat,status)
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: bhat(3),grad_bhat(3,3)
    integer, intent(out) :: status

    double precision :: xcart(3),xplus_cart(3),xminus_cart(3)
    double precision :: xplus(ndim),xminus(ndim),bplus(3),bminus(3)
    double precision :: Btmp(3),eps
    integer :: j,igrid_plus,igrid_minus,status_plus,status_minus

    bhat=zero
    grad_bhat=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical) return

    {^IFTHREED
    call trace_spherical_sample_B_bhat_cart(x,igrid,threshold,Btmp,bhat, &
         status)
    if (status/=trace_status_active) return
    call trace_spherical_coord_to_cart(x,xcart,status)
    if (status/=trace_status_active) return

    eps=max(1.d-8*max(one,abs(x(1))), &
         0.01d0*trace_spherical_physical_cell_scale(x,igrid))
    do j=1,3
      xplus_cart=xcart
      xminus_cart=xcart
      xplus_cart(j)=xplus_cart(j)+eps
      xminus_cart(j)=xminus_cart(j)-eps

      status_plus=trace_status_bad_grad_stencil
      call trace_cart_to_spherical_coord(xplus_cart,xplus,status_plus)
      if (status_plus==trace_status_active) then
        call trace_locate_point_with_hint(xplus,igrid,igrid_plus,status_plus)
      endif
      if (status_plus==trace_status_active) then
        call trace_spherical_sample_B_bhat_cart(xplus,igrid_plus,threshold, &
             Btmp,bplus,status_plus)
      endif

      status_minus=trace_status_bad_grad_stencil
      call trace_cart_to_spherical_coord(xminus_cart,xminus,status_minus)
      if (status_minus==trace_status_active) then
        call trace_locate_point_with_hint(xminus,igrid,igrid_minus, &
             status_minus)
      endif
      if (status_minus==trace_status_active) then
        call trace_spherical_sample_B_bhat_cart(xminus,igrid_minus,threshold, &
             Btmp,bminus,status_minus)
      endif

      if (status_plus==trace_status_active .and. &
           status_minus==trace_status_active) then
        grad_bhat(:,j)=(bplus-bminus)/(2.d0*eps)
      else if (status_plus==trace_status_active) then
        grad_bhat(:,j)=(bplus-bhat)/eps
      else if (status_minus==trace_status_active) then
        grad_bhat(:,j)=(bhat-bminus)/eps
      else
        status=trace_status_bad_grad_stencil
        return
      endif
    enddo
    status=trace_status_active
    }
  end subroutine trace_spherical_sample_bhat_gradbhat_cartfd

  subroutine trace_spherical_bhat_cart_to_rhs(x,bhat_cart,kx,status)
    double precision, intent(in) :: x(ndim),bhat_cart(3)
    double precision, intent(out) :: kx(ndim)
    integer, intent(out) :: status

    double precision :: er(3),etheta(3),ephi(3),r,sint

    kx=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical) return

    {^IFTHREED
    call trace_spherical_basis(x,er,etheta,ephi,status)
    if (status/=trace_status_active) return
    r=x(1)
    sint=dsin(x(2))
    if (r<=smalldouble .or. abs(sint)<=smalldouble) then
      status=trace_status_invalid_input
      return
    endif
    kx(1)=sum(bhat_cart*er)
    kx(2)=sum(bhat_cart*etheta)/r
    kx(3)=sum(bhat_cart*ephi)/(r*sint)
    status=trace_status_active
    }
  end subroutine trace_spherical_bhat_cart_to_rhs

  subroutine trace_spherical_endpoint_B_bhat(x,face_id,igrid,threshold,B, &
       bhat,status)
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: face_id,igrid
    double precision, intent(out) :: B(ndim),bhat(ndim)
    integer, intent(out) :: status

    double precision :: xprobe(ndim),Bcart(3),bhat_cart(3)
    integer :: normal_dim,sample_igrid

    B=zero
    bhat=zero
    status=trace_status_bad_face_limit_sample
    if (ndim/=3 .or. geo_coordinate/=geo_spherical .or. igrid<0) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (.not.trace_face_is_boundary(face_id)) return

    {^IFTHREED
    select case(face_id)
    case(trace_face_xmin)
      normal_dim=1
    case(trace_face_xmax)
      normal_dim=1
    case(trace_face_ymin)
      normal_dim=2
    case(trace_face_ymax)
      normal_dim=2
    case(trace_face_zmin)
      normal_dim=3
    case(trace_face_zmax)
      normal_dim=3
    case default
      return
    end select

    xprobe=x
    call trace_spherical_face_probe_coord(face_id,igrid, &
         xprobe(normal_dim),status)
    if (status/=trace_status_active) return
    xprobe(1)=max(xprobmin1,min(xprobe(1), &
         xprobmax1-100.d0*epsilon(one)*max(one,abs(xprobmax1))))
    xprobe(2)=max(xprobmin2,min(xprobe(2), &
         xprobmax2-100.d0*epsilon(one)*max(one,abs(xprobmax2))))
    xprobe(3)=max(xprobmin3,min(xprobe(3), &
         xprobmax3-100.d0*epsilon(one)*max(one,abs(xprobmax3))))
    call trace_debug_locate_point(xprobe,sample_igrid,status)
    if (status/=trace_status_active) return
    call trace_spherical_sample_B_bhat_cart(xprobe,sample_igrid,threshold, &
         Bcart,bhat_cart,status)
    if (status/=trace_status_active) return
    B=Bcart(1:ndim)
    bhat=bhat_cart(1:ndim)
    status=trace_status_boundary
    }
  end subroutine trace_spherical_endpoint_B_bhat

  subroutine trace_spherical_face_probe_coord(face_id,igrid,xprobe,status)
    integer, intent(in) :: face_id,igrid
    double precision, intent(out) :: xprobe
    integer, intent(out) :: status

    integer :: i1,i2,i3,normal_dim,side

    if (trace_spherical_profile_enabled) then
      call trace_spherical_profile_add_count(trace_profile_face_probes,1_8)
    endif

    xprobe=zero
    status=trace_status_bad_face_limit_sample
    if (ndim/=3 .or. geo_coordinate/=geo_spherical .or. igrid<0) then
      status=trace_status_unsupported_geometry
      return
    endif

    {^IFTHREED
    i1=ixMlo1
    i2=ixMlo2
    i3=ixMlo3
    select case(face_id)
    case(trace_face_xmin)
      normal_dim=1
      side=-1
      i1=ixMlo1
    case(trace_face_xmax)
      normal_dim=1
      side=1
      i1=ixMhi1
    case(trace_face_ymin)
      normal_dim=2
      side=-1
      i2=ixMlo2
    case(trace_face_ymax)
      normal_dim=2
      side=1
      i2=ixMhi2
    case(trace_face_zmin)
      normal_dim=3
      side=-1
      i3=ixMlo3
    case(trace_face_zmax)
      normal_dim=3
      side=1
      i3=ixMhi3
    case default
      return
    end select

    xprobe=ps(igrid)%x(i1,i2,i3,normal_dim)
    if (side<0) then
      if (xprobe<xprobmin1 .and. normal_dim==1) return
      if (xprobe<xprobmin2 .and. normal_dim==2) return
      if (xprobe<xprobmin3 .and. normal_dim==3) return
    else
      if (xprobe>xprobmax1 .and. normal_dim==1) return
      if (xprobe>xprobmax2 .and. normal_dim==2) return
      if (xprobe>xprobmax3 .and. normal_dim==3) return
    endif
    status=trace_status_active
    }
  end subroutine trace_spherical_face_probe_coord

  double precision function trace_spherical_physical_cell_scale(x,igrid) &
       result(cell_scale)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: igrid

    double precision :: r,sint,dxloc(ndim)
    integer :: status

    cell_scale=one
    if (ndim/=3 .or. geo_coordinate/=geo_spherical .or. igrid<0) return

    {^IFTHREED
    call trace_spherical_local_cell_widths(x,igrid,dxloc,status)
    if (status/=trace_status_active) return
    r=max(abs(x(1)),smalldouble)
    sint=max(abs(dsin(x(2))),smalldouble)
    cell_scale=trace_spherical_cell_scale_from_widths(r,sint,dxloc)
    }
  end function trace_spherical_physical_cell_scale

  double precision function trace_spherical_cell_scale_from_widths(r,sint, &
       dxloc) result(cell_scale)
    double precision, intent(in) :: r,sint,dxloc(ndim)

    cell_scale=one
    {^IFTHREED
    cell_scale=min(abs(dxloc(1)),abs(r)*abs(dxloc(2)))
    cell_scale=min(cell_scale,abs(r)*max(abs(sint),smalldouble)* &
         abs(dxloc(3)))
    cell_scale=max(cell_scale,smalldouble)
    }
  end function trace_spherical_cell_scale_from_widths

  subroutine trace_spherical_global_min_cell_size(hmin,status)
    double precision, intent(out) :: hmin
    integer, intent(out) :: status

    double precision :: dxloc(ndim),hcell,r,sint
    integer :: iigrid,igrid,ix1,ix2,ix3

    hmin=huge(one)
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. geo_coordinate/=geo_spherical) return

    {^IFTHREED
    status=trace_status_out_of_domain
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      do ix3=ixMlo3,ixMhi3
        do ix2=ixMlo2,ixMhi2
          do ix1=ixMlo1,ixMhi1
            r=ps(igrid)%x(ix1,ix2,ix3,1)
            sint=dsin(ps(igrid)%x(ix1,ix2,ix3,2))
            if (r<=smalldouble .or. abs(sint)<=smalldouble) cycle
            dxloc(1)=abs(ps(igrid)%dx(ix1,ix2,ix3,1))
            dxloc(2)=abs(ps(igrid)%dx(ix1,ix2,ix3,2))
            dxloc(3)=abs(ps(igrid)%dx(ix1,ix2,ix3,3))
            hcell=trace_spherical_cell_scale_from_widths(r,sint,dxloc)
            hmin=min(hmin,hcell)
          enddo
        enddo
      enddo
    enddo
    if (hmin<huge(one)) status=trace_status_active
    }
  end subroutine trace_spherical_global_min_cell_size

  subroutine trace_cartesian_global_min_cell_size(hmin,status)
    double precision, intent(out) :: hmin
    integer, intent(out) :: status

    double precision :: hlocal,hcell
    integer :: iigrid,igrid,ix1,ix2,ix3

    hmin=huge(one)
    hlocal=huge(one)
    status=trace_status_unsupported_geometry
    if (ndim/=3 .or. .not.trace_cartesian_like_geometry()) return

    {^IFTHREED
    status=trace_status_out_of_domain
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      do ix3=ixMlo3,ixMhi3
        do ix2=ixMlo2,ixMhi2
          do ix1=ixMlo1,ixMhi1
            hcell=min(abs(ps(igrid)%dx(ix1,ix2,ix3,1)), &
                 abs(ps(igrid)%dx(ix1,ix2,ix3,2)))
            hcell=min(hcell,abs(ps(igrid)%dx(ix1,ix2,ix3,3)))
            if (hcell>zero) hlocal=min(hlocal,hcell)
          enddo
        enddo
      enddo
    enddo
    if (npe>1) then
      call MPI_ALLREDUCE(hlocal,hmin,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
           icomm,ierrmpi)
    else
      hmin=hlocal
    endif
    if (hmin<huge(one)) status=trace_status_active
    }
  end subroutine trace_cartesian_global_min_cell_size

  logical function trace_cartesian_like_geometry() result(is_cart_like)
    is_cart_like=.false.
    if (ndim/=3) return
    {^IFTHREED
    select case (geo_coordinate)
    case (geo_cartesian,geo_cartesian_stretched)
      is_cart_like=.true.
    end select
    }
  end function trace_cartesian_like_geometry

  logical function trace_rk45_position_integrator() result(is_rk45)
    is_rk45=.false.
    if (ndim/=3) return
    {^IFTHREED
    select case (trace_integrator_mode)
    case (trace_integrator_rk45_cartesian)
      is_rk45=trace_cartesian_like_geometry()
    case (trace_integrator_rk45_spherical)
      is_rk45=geo_coordinate==geo_spherical
    end select
    }
  end function trace_rk45_position_integrator

  subroutine trace_interp_weights_block(x,igrid,ixI^L,ixbl^D,xd^D,dxc^D, &
       status)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: igrid,ixI^L
    integer, intent(out) :: ixbl^D,status
    double precision, intent(out) :: xd^D,dxc^D

    status=trace_status_out_of_domain
    ^D&ixbl^D=ixImin^D;
    ^D&xd^D=zero;
    ^D&dxc^D=zero;
    if (igrid<0) return

    {^IFTHREED
    call trace_interp_index_3d(x(1),igrid,1,ixI^L,ixbl1,xd1,dxc1,status)
    if (status/=trace_status_active) return
    call trace_interp_index_3d(x(2),igrid,2,ixI^L,ixbl2,xd2,dxc2,status)
    if (status/=trace_status_active) return
    call trace_interp_index_3d(x(3),igrid,3,ixI^L,ixbl3,xd3,dxc3,status)
    }
  end subroutine trace_interp_weights_block

  subroutine trace_interp_weights_block_near(x,igrid,ixI^L,ixstart,ixbl^D, &
       xd^D,dxc^D,status)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: igrid,ixI^L
    integer, intent(in) :: ixstart(3)
    integer, intent(out) :: ixbl^D,status
    double precision, intent(out) :: xd^D,dxc^D

    status=trace_status_out_of_domain
    ^D&ixbl^D=ixImin^D;
    ^D&xd^D=zero;
    ^D&dxc^D=zero;
    if (igrid<0) return

    {^IFTHREED
    call trace_interp_index_3d_near(x(1),igrid,1,ixI^L,ixstart(1), &
         ixbl1,xd1,dxc1,status)
    if (status/=trace_status_active) return
    call trace_interp_index_3d_near(x(2),igrid,2,ixI^L,ixstart(2), &
         ixbl2,xd2,dxc2,status)
    if (status/=trace_status_active) return
    call trace_interp_index_3d_near(x(3),igrid,3,ixI^L,ixstart(3), &
         ixbl3,xd3,dxc3,status)
    }
  end subroutine trace_interp_weights_block_near

  subroutine trace_interp_index_3d(xval,igrid,idim,ixI^L,ixlo,frac,dxc, &
       status)
    double precision, intent(in) :: xval
    integer, intent(in) :: igrid,idim,ixI^L
    integer, intent(out) :: ixlo,status
    double precision, intent(out) :: frac,dxc

    double precision :: xlo,xhi,tol
    integer :: i

    ixlo=ixImin1
    frac=zero
    dxc=zero
    status=trace_status_out_of_domain
    if (igrid<0) return

    {^IFTHREED
    select case(idim)
    case(1)
      do i=ixImin1,ixImax1-1
        xlo=ps(igrid)%x(i,ixImin2,ixImin3,1)
        xhi=ps(igrid)%x(i+1,ixImin2,ixImin3,1)
        dxc=xhi-xlo
        tol=100.d0*epsilon(one)*max(one,max(abs(xlo),abs(xhi)))
        if (abs(dxc)<=smalldouble) cycle
        if (xval>=min(xlo,xhi)-tol .and. &
             xval<=max(xlo,xhi)+tol) then
          ixlo=i
          frac=max(zero,min(one,(xval-xlo)/dxc))
          status=trace_status_active
          return
        endif
      enddo
    case(2)
      do i=ixImin2,ixImax2-1
        xlo=ps(igrid)%x(ixImin1,i,ixImin3,2)
        xhi=ps(igrid)%x(ixImin1,i+1,ixImin3,2)
        dxc=xhi-xlo
        tol=100.d0*epsilon(one)*max(one,max(abs(xlo),abs(xhi)))
        if (abs(dxc)<=smalldouble) cycle
        if (xval>=min(xlo,xhi)-tol .and. &
             xval<=max(xlo,xhi)+tol) then
          ixlo=i
          frac=max(zero,min(one,(xval-xlo)/dxc))
          status=trace_status_active
          return
        endif
      enddo
    case(3)
      do i=ixImin3,ixImax3-1
        xlo=ps(igrid)%x(ixImin1,ixImin2,i,3)
        xhi=ps(igrid)%x(ixImin1,ixImin2,i+1,3)
        dxc=xhi-xlo
        tol=100.d0*epsilon(one)*max(one,max(abs(xlo),abs(xhi)))
        if (abs(dxc)<=smalldouble) cycle
        if (xval>=min(xlo,xhi)-tol .and. &
             xval<=max(xlo,xhi)+tol) then
          ixlo=i
          frac=max(zero,min(one,(xval-xlo)/dxc))
          status=trace_status_active
          return
        endif
      enddo
    end select
    }
  end subroutine trace_interp_index_3d

  subroutine trace_interp_index_3d_near(xval,igrid,idim,ixI^L,ixstart, &
       ixlo,frac,dxc,status)
    double precision, intent(in) :: xval
    integer, intent(in) :: igrid,idim,ixI^L,ixstart
    integer, intent(out) :: ixlo,status
    double precision, intent(out) :: frac,dxc

    double precision :: xlo,xhi,tol
    integer :: i,ilo,ihi,idir

    ixlo=ixImin1
    frac=zero
    dxc=zero
    status=trace_status_out_of_domain
    if (igrid<0) return

    {^IFTHREED
    select case(idim)
    case(1)
      ilo=ixImin1
      ihi=ixImax1-1
    case(2)
      ilo=ixImin2
      ihi=ixImax2-1
    case(3)
      ilo=ixImin3
      ihi=ixImax3-1
    case default
      return
    end select

    i=max(ilo,min(ixstart,ihi))
    do
      select case(idim)
      case(1)
        xlo=ps(igrid)%x(i,ixImin2,ixImin3,1)
        xhi=ps(igrid)%x(i+1,ixImin2,ixImin3,1)
      case(2)
        xlo=ps(igrid)%x(ixImin1,i,ixImin3,2)
        xhi=ps(igrid)%x(ixImin1,i+1,ixImin3,2)
      case(3)
        xlo=ps(igrid)%x(ixImin1,ixImin2,i,3)
        xhi=ps(igrid)%x(ixImin1,ixImin2,i+1,3)
      end select
      dxc=xhi-xlo
      tol=100.d0*epsilon(one)*max(one,max(abs(xlo),abs(xhi)))
      if (abs(dxc)>smalldouble .and. &
           xval>=min(xlo,xhi)-tol .and. &
           xval<=max(xlo,xhi)+tol) then
        ixlo=i
        frac=max(zero,min(one,(xval-xlo)/dxc))
        status=trace_status_active
        return
      endif

      if (abs(dxc)<=smalldouble) then
        idir=1
      else if (xval<min(xlo,xhi)-tol) then
        idir=-1
      else
        idir=1
      endif

      i=i+idir
      if (i<ilo .or. i>ihi) exit
    enddo
    }
  end subroutine trace_interp_index_3d_near

  subroutine trace_spherical_local_cell_widths(x,igrid,dxloc,status)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: igrid
    double precision, intent(out) :: dxloc(ndim)
    integer, intent(out) :: status

    type(trace_sph_interp_ctx) :: ctx

    dxloc=one
    status=trace_status_out_of_domain
    if (ndim/=3 .or. geo_coordinate/=geo_spherical .or. igrid<0) return

    {^IFTHREED
    call trace_spherical_interp_ctx_build(x,igrid,ctx,status)
    if (status/=trace_status_active) return
    dxloc(1:3)=ctx%dxloc(1:3)
    }
  end subroutine trace_spherical_local_cell_widths

  subroutine trace_debug_locate_point(x,igrid,status)
    use mod_particle_base, only: find_particle_ipe

    double precision, intent(in) :: x(ndim)
    integer, intent(out) :: igrid,status

    double precision :: x3d(3)
    integer :: indomain,ipe,j

    igrid=-1
    status=trace_status_active
    indomain=0
    {if (x(^DB)>=xprobmin^DB .and. x(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain/=ndim) then
      status=trace_status_boundary
      return
    endif

    x3d=zero
    do j=1,ndim
      x3d(j)=x(j)
    enddo
    call find_particle_ipe(x3d,igrid,ipe)
    if (igrid<0 .or. ipe/=mype) then
      status=trace_status_out_of_domain
    endif
  end subroutine trace_debug_locate_point

  subroutine trace_locate_point_with_hint(x,igrid_hint,igrid,status)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: igrid_hint
    integer, intent(out) :: igrid,status

    double precision :: xbmin^D,xbmax^D
    integer :: inblock

    igrid=-1
    status=trace_status_active
    if (igrid_hint>=0) then
      ^D&xbmin^D=rnode(rpxmin^D_,igrid_hint);
      ^D&xbmax^D=rnode(rpxmax^D_,igrid_hint);
      inblock=0
      {if (x(^DB)>xbmin^DB .and. x(^DB)<xbmax^DB) inblock=inblock+1\}
      if (inblock==ndim) then
        igrid=igrid_hint
        return
      endif
    endif

    call trace_debug_locate_point(x,igrid,status)
  end subroutine trace_locate_point_with_hint

  subroutine trace_tangent_rhs(x,u,v,igrid,threshold,kx,ku,kv,status, &
       p,q,kp,kq,sph_ctx)
    double precision, intent(in) :: x(ndim),u(ndim),v(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: kx(ndim),ku(ndim),kv(ndim)
    integer, intent(out) :: status
    double precision, intent(in), optional :: p(ndim),q(ndim)
    double precision, intent(out), optional :: kp(ndim),kq(ndim)
    type(trace_sph_interp_ctx), intent(in), optional :: sph_ctx

    double precision :: bhat(3),grad_bhat(3,3)
    double precision :: u3(3),v3(3),p3(3),q3(3)
    double precision :: ku3(3),kv3(3),kp3(3),kq3(3)
    logical :: use_ctx

    kx=zero
    ku=zero
    kv=zero
    if (present(kp)) kp=zero
    if (present(kq)) kq=zero
    status=trace_status_unsupported_geometry
    if (ndim/=3) return
    if (trace_rk2_stats_enabled) call trace_rk2_stats_note_rhs()

    {^IFTHREED
    if (geo_coordinate==geo_spherical) then
      use_ctx=.false.
      if (present(sph_ctx)) then
        use_ctx=sph_ctx%valid .and. sph_ctx%igrid==igrid .and. &
             maxval(abs(sph_ctx%x(1:3)-x(1:3))) <= &
             100.d0*epsilon(one)*max(one,maxval(abs(x(1:3))))
      endif
      if (use_ctx) then
        call trace_spherical_sample_bhat_gradbhat_covariant_ctx(sph_ctx, &
             threshold,bhat,grad_bhat,status)
      else
        call trace_spherical_sample_bhat_gradbhat_covariant(x,igrid, &
             threshold,bhat,grad_bhat,status)
      endif
      if (status==trace_status_bad_grad_stencil) then
        if (trace_spherical_profile_enabled) then
          call trace_spherical_profile_add_count(trace_profile_grad_fallbacks, &
               1_8)
        endif
        call trace_spherical_sample_bhat_gradbhat_cartfd(x,igrid,threshold, &
           bhat,grad_bhat,status)
      endif
      if (status/=trace_status_active) return
      call trace_spherical_bhat_cart_to_rhs(x,bhat,kx,status)
      if (status/=trace_status_active) return
    else
      call sample_bhat_gradbhat_at_point(x,igrid,threshold,bhat,grad_bhat,status)
      if (status/=trace_status_active) return
      kx=bhat(1:ndim)
    endif
    if (status/=trace_status_active) return

    u3=zero
    v3=zero
    u3(1:ndim)=u
    v3(1:ndim)=v
    ku3=matmul(grad_bhat,u3)
    kv3=matmul(grad_bhat,v3)
    ku=ku3(1:ndim)
    kv=kv3(1:ndim)
    if (present(p) .and. present(q) .and. present(kp) .and. &
         present(kq)) then
      p3=zero
      q3=zero
      p3(1:ndim)=p
      q3(1:ndim)=q
      kp3=matmul(grad_bhat,p3)
      kq3=matmul(grad_bhat,q3)
      kp=kp3(1:ndim)
      kq=kq3(1:ndim)
    endif
    }
  end subroutine trace_tangent_rhs

  subroutine trace_advance_tangent_state_rk2(x,u,v,igrid,h,threshold,status)
    double precision, intent(inout) :: x(ndim),u(ndim),v(ndim)
    integer, intent(inout) :: igrid
    double precision, intent(in) :: h,threshold
    integer, intent(out) :: status

    double precision :: xnew(ndim),unew(ndim),vnew(ndim)
    integer :: igrid_mid,igrid_new

    call trace_tangent_rk2_trial(x,u,v,igrid,h,threshold,xnew,unew,vnew, &
         status)
    if (status/=trace_status_active) return

    call trace_debug_locate_point(xnew,igrid_new,status)
    if (status/=trace_status_active) return

    x=xnew
    u=unew
    v=vnew
    igrid=igrid_new
  end subroutine trace_advance_tangent_state_rk2

  subroutine trace_tangent_rk2_trial(x,u,v,igrid,h,threshold,xnew,unew, &
       vnew,status)
    double precision, intent(in) :: x(ndim),u(ndim),v(ndim),h,threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: xnew(ndim),unew(ndim),vnew(ndim)
    integer, intent(out) :: status

    double precision :: kx1(ndim),ku1(ndim),kv1(ndim)

    call trace_tangent_rhs(x,u,v,igrid,threshold,kx1,ku1,kv1,status)
    if (status/=trace_status_active) return

    call trace_tangent_rk2_trial_from_rhs(x,u,v,igrid,h,threshold, &
         kx1,ku1,kv1,xnew,unew,vnew,status)
  end subroutine trace_tangent_rk2_trial

  subroutine trace_tangent_rk2_trial_from_rhs(x,u,v,igrid,h,threshold, &
       kx1,ku1,kv1,xnew,unew,vnew,status,p,q,kp1,kq1,pnew,qnew, &
       sph_cache,clamp_stage)
    double precision, intent(in) :: x(ndim),u(ndim),v(ndim),h,threshold
    double precision, intent(in) :: kx1(ndim),ku1(ndim),kv1(ndim)
    integer, intent(in) :: igrid
    double precision, intent(out) :: xnew(ndim),unew(ndim),vnew(ndim)
    integer, intent(out) :: status
    double precision, intent(in), optional :: p(ndim),q(ndim)
    double precision, intent(in), optional :: kp1(ndim),kq1(ndim)
    double precision, intent(out), optional :: pnew(ndim),qnew(ndim)
    type(trace_sph_interp_ctx), intent(inout), optional :: sph_cache
    logical, intent(in), optional :: clamp_stage

    double precision :: kx2(ndim),ku2(ndim),kv2(ndim),kp2(ndim),kq2(ndim)
    double precision :: xmid(ndim),xmid_sample(ndim)
    double precision :: umid(ndim),vmid(ndim),pmid(ndim),qmid(ndim)
    type(trace_sph_interp_ctx) :: ctx_mid
    integer :: igrid_mid
    integer :: ctx_status
    logical :: use_clamped_stage

    xmid=x+half*h*kx1
    umid=u+half*h*ku1
    vmid=v+half*h*kv1
    use_clamped_stage=.false.
    if (present(clamp_stage)) use_clamped_stage=clamp_stage
    xmid_sample=xmid
    if (use_clamped_stage) then
      ^D&xmid_sample(^D)=max(xprobmin^D,min(xmid_sample(^D), &
           xprobmax^D-100.d0*epsilon(one)*max(one,abs(xprobmax^D))));
    endif
    call trace_locate_point_with_hint(xmid_sample,igrid,igrid_mid,status)
    if (status/=trace_status_active) return
    ctx_status=trace_status_unsupported_geometry
    if (geo_coordinate==geo_spherical) then
      if (present(sph_cache)) then
        call trace_spherical_interp_ctx_build_cached(xmid_sample,igrid_mid, &
             sph_cache,ctx_mid,ctx_status)
      else
        call trace_spherical_interp_ctx_build(xmid_sample,igrid_mid,ctx_mid, &
             ctx_status)
      endif
    endif

    if (present(p) .and. present(q) .and. present(kp1) .and. &
         present(kq1) .and. present(pnew) .and. present(qnew)) then
      pmid=p+half*h*kp1
      qmid=q+half*h*kq1
      if (ctx_status==trace_status_active) then
        call trace_tangent_rhs(xmid_sample,umid,vmid,igrid_mid,threshold,kx2,ku2, &
             kv2,status,pmid,qmid,kp2,kq2,sph_ctx=ctx_mid)
      else
        call trace_tangent_rhs(xmid_sample,umid,vmid,igrid_mid,threshold,kx2,ku2, &
             kv2,status,pmid,qmid,kp2,kq2)
      endif
    else
      if (ctx_status==trace_status_active) then
        call trace_tangent_rhs(xmid_sample,umid,vmid,igrid_mid,threshold,kx2,ku2, &
             kv2,status,sph_ctx=ctx_mid)
      else
        call trace_tangent_rhs(xmid_sample,umid,vmid,igrid_mid,threshold,kx2,ku2, &
             kv2,status)
      endif
    endif
    if (status/=trace_status_active) return

    xnew=x+h*kx2
    unew=u+h*ku2
    vnew=v+h*kv2
    if (present(pnew) .and. present(qnew)) then
      pnew=p+h*kp2
      qnew=q+h*kq2
    endif
  end subroutine trace_tangent_rk2_trial_from_rhs

  subroutine trace_project_to_perp_bhat(vec,bhat,vec_perp)
    double precision, intent(in) :: vec(ndim),bhat(ndim)
    double precision, intent(out) :: vec_perp(ndim)

    vec_perp=vec-sum(vec*bhat)*bhat
  end subroutine trace_project_to_perp_bhat

  subroutine trace_endpoint_B_bhat(x,face_id,igrid,threshold,B,bhat,status)
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: face_id,igrid
    double precision, intent(out) :: B(ndim),bhat(ndim)
    integer, intent(out) :: status

    double precision :: B3(3),bhat3(3)
    integer :: sample_igrid

    B=zero
    bhat=zero
    if (geo_coordinate==geo_spherical) then
      call trace_spherical_endpoint_B_bhat(x,face_id,igrid,threshold,B,bhat, &
           status)
      return
    endif
    call trace_locate_face_limit_grid(x,face_id,igrid,sample_igrid,status)
    if (status/=trace_status_active) return
    call trace_sample_B_on_domain_face_limit(x,face_id,sample_igrid,threshold, &
         B3,bhat3,status)
    if (status/=trace_status_active) return
    B=B3(1:ndim)
    bhat=bhat3(1:ndim)
    status=trace_status_boundary
  end subroutine trace_endpoint_B_bhat

  subroutine trace_locate_face_limit_grid(x,face_id,igrid_in,igrid_out, &
       status)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: face_id,igrid_in
    integer, intent(out) :: igrid_out,status

    double precision :: xprobe(ndim),xface
    integer :: ixO^L
    integer :: normal_dim,side

    igrid_out=-1
    status=trace_status_bad_face_limit_sample
    if (ndim/=3 .or. .not.trace_cartesian_like_geometry() .or. &
         igrid_in<0) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (.not.trace_face_is_boundary(face_id)) return

    {^IFTHREED
    select case(face_id)
    case(trace_face_xmin)
      normal_dim=1
      side=-1
      xface=xprobmin1
    case(trace_face_xmax)
      normal_dim=1
      side=1
      xface=xprobmax1
    case(trace_face_ymin)
      normal_dim=2
      side=-1
      xface=xprobmin2
    case(trace_face_ymax)
      normal_dim=2
      side=1
      xface=xprobmax2
    case(trace_face_zmin)
      normal_dim=3
      side=-1
      xface=xprobmin3
    case(trace_face_zmax)
      normal_dim=3
      side=1
      xface=xprobmax3
    case default
      return
    end select

    ixO^L=ixM^LL;
    xprobe=x
    select case(normal_dim)
    case(1)
      if (side<0) then
        xprobe(normal_dim)=ps(igrid_in)%x(ixOmin1,ixOmin2,ixOmin3,1)
      else
        xprobe(normal_dim)=ps(igrid_in)%x(ixOmax1,ixOmin2,ixOmin3,1)
      endif
    case(2)
      if (side<0) then
        xprobe(normal_dim)=ps(igrid_in)%x(ixOmin1,ixOmin2,ixOmin3,2)
      else
        xprobe(normal_dim)=ps(igrid_in)%x(ixOmin1,ixOmax2,ixOmin3,2)
      endif
    case(3)
      if (side<0) then
        xprobe(normal_dim)=ps(igrid_in)%x(ixOmin1,ixOmin2,ixOmin3,3)
      else
        xprobe(normal_dim)=ps(igrid_in)%x(ixOmin1,ixOmin2,ixOmax3,3)
      endif
    end select
    call trace_debug_locate_point(xprobe,igrid_out,status)
    }
  end subroutine trace_locate_face_limit_grid

  subroutine trace_sample_B_on_domain_face_limit(xhit,face_id,igrid, &
       threshold,B,bhat,status)
    ! Sample total B at a rectangular domain face from the inside-side limit.
    ! The normal direction is linearly extrapolated from the two nearest
    ! interior cell centers; tangential directions use bilinear interpolation.
    double precision, intent(in) :: xhit(ndim),threshold
    integer, intent(in) :: face_id,igrid
    double precision, intent(out) :: B(3),bhat(3)
    integer, intent(out) :: status

    double precision :: Bnear(3),Binner(3)
    double precision :: xface,xnear,xinner
    double precision :: xnear_pt(ndim),xinner_pt(ndim)
    double precision :: Bnorm
    integer :: ixO^L
    integer :: normal_dim,side
    integer :: igrid_near,igrid_inner

    B=zero
    bhat=zero
    status=trace_status_bad_face_limit_sample
    if (ndim/=3 .or. .not.trace_cartesian_like_geometry() .or. igrid<0) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (.not.trace_face_is_boundary(face_id)) return

    {^IFTHREED
    select case(face_id)
    case(trace_face_xmin)
      normal_dim=1
      side=-1
      xface=xprobmin1
    case(trace_face_xmax)
      normal_dim=1
      side=1
      xface=xprobmax1
    case(trace_face_ymin)
      normal_dim=2
      side=-1
      xface=xprobmin2
    case(trace_face_ymax)
      normal_dim=2
      side=1
      xface=xprobmax2
    case(trace_face_zmin)
      normal_dim=3
      side=-1
      xface=xprobmin3
    case(trace_face_zmax)
      normal_dim=3
      side=1
      xface=xprobmax3
    case default
      return
    end select

    ixO^L=ixM^LL;
    select case(normal_dim)
    case(1)
      if (side<0) then
        xnear=ps(igrid)%x(ixOmin1,ixOmin2,ixOmin3,1)
        xinner=ps(igrid)%x(ixOmin1+1,ixOmin2,ixOmin3,1)
      else
        xnear=ps(igrid)%x(ixOmax1,ixOmin2,ixOmin3,1)
        xinner=ps(igrid)%x(ixOmax1-1,ixOmin2,ixOmin3,1)
      endif
    case(2)
      if (side<0) then
        xnear=ps(igrid)%x(ixOmin1,ixOmin2,ixOmin3,2)
        xinner=ps(igrid)%x(ixOmin1,ixOmin2+1,ixOmin3,2)
      else
        xnear=ps(igrid)%x(ixOmin1,ixOmax2,ixOmin3,2)
        xinner=ps(igrid)%x(ixOmin1,ixOmax2-1,ixOmin3,2)
      endif
    case(3)
      if (side<0) then
        xnear=ps(igrid)%x(ixOmin1,ixOmin2,ixOmin3,3)
        xinner=ps(igrid)%x(ixOmin1,ixOmin2,ixOmin3+1,3)
      else
        xnear=ps(igrid)%x(ixOmin1,ixOmin2,ixOmax3,3)
        xinner=ps(igrid)%x(ixOmin1,ixOmin2,ixOmax3-1,3)
      endif
    end select
    xnear_pt=xhit
    xinner_pt=xhit
    xnear_pt(normal_dim)=xnear
    xinner_pt(normal_dim)=xinner

    call trace_debug_locate_point(xnear_pt,igrid_near,status)
    if (status/=trace_status_active) return
    call sample_B_at_point(xnear_pt,igrid_near,threshold,Bnear,status)
    if (status/=trace_status_active) return
    call trace_debug_locate_point(xinner_pt,igrid_inner,status)
    if (status/=trace_status_active) return
    call sample_B_at_point(xinner_pt,igrid_inner,threshold,Binner,status)
    if (status/=trace_status_active) return

    B=Bnear+(Bnear-Binner)/(xnear-xinner)*(xface-xnear)
    Bnorm=dsqrt(sum(B**2))
    if (Bnorm<=zero .or. Bnorm<threshold) then
      status=trace_status_weak_field
      return
    endif
    bhat=B/Bnorm
    status=trace_status_active
    }
  end subroutine trace_sample_B_on_domain_face_limit

  subroutine trace_face_normal(face_id,normal,ok)
    integer, intent(in) :: face_id
    double precision, intent(out) :: normal(ndim)
    logical, intent(out) :: ok

    normal=zero
    ok=.true.
    select case(face_id)
    case(trace_face_xmin)
      normal(1)=-one
    case(trace_face_xmax)
      normal(1)=one
    {^IFONED
    case default
      ok=.false.
    }
    {^IFTWOD
    case(trace_face_ymin)
      normal(2)=-one
    case(trace_face_ymax)
      normal(2)=one
    case default
      ok=.false.
    }
    {^IFTHREED
    case(trace_face_ymin)
      normal(2)=-one
    case(trace_face_ymax)
      normal(2)=one
    case(trace_face_zmin)
      normal(3)=-one
    case(trace_face_zmax)
      normal(3)=one
    case default
      ok=.false.
    }
    end select
  end subroutine trace_face_normal

  subroutine trace_project_to_boundary_face(vec,bhat,face_id,vec_face)
    double precision, intent(in) :: vec(ndim),bhat(ndim)
    integer, intent(in) :: face_id
    double precision, intent(out) :: vec_face(ndim)

    double precision :: normal(ndim),bn
    logical :: ok

    vec_face=zero
    call trace_face_normal(face_id,normal,ok)
    if (.not.ok) return
    bn=sum(bhat*normal)
    if (abs(bn)<=smalldouble) return
    vec_face=vec-sum(vec*normal)/bn*bhat
  end subroutine trace_project_to_boundary_face

  subroutine trace_make_perp_basis(bhat,u0,v0,status)
    double precision, intent(in) :: bhat(3)
    double precision, intent(out) :: u0(ndim),v0(ndim)
    integer, intent(out) :: status

    double precision :: ref(3),u3(3),v3(3),b3(3),unorm,bnorm
    integer :: iref

    u0=zero
    v0=zero
    status=trace_status_active
    if (ndim/=3) then
      status=trace_status_unsupported_geometry
      return
    endif

    b3=bhat
    bnorm=dsqrt(sum(b3**2))
    if (bnorm<=zero) then
      status=trace_status_weak_field
      return
    endif
    b3=b3/bnorm

    iref=1
    if (abs(b3(2))<abs(b3(iref))) iref=2
    if (abs(b3(3))<abs(b3(iref))) iref=3
    ref=zero
    ref(iref)=one

    u3=ref-sum(ref*b3)*b3
    unorm=dsqrt(sum(u3**2))
    if (unorm<=smalldouble) then
      status=trace_status_weak_field
      return
    endif
    u3=u3/unorm

    v3(1)=b3(2)*u3(3)-b3(3)*u3(2)
    v3(2)=b3(3)*u3(1)-b3(1)*u3(3)
    v3(3)=b3(1)*u3(2)-b3(2)*u3(1)
    v3=v3/dsqrt(sum(v3**2))

    u0=u3(1:ndim)
    v0=v3(1:ndim)
  end subroutine trace_make_perp_basis

  double precision function trace_debug_nan()
    trace_debug_nan=ieee_value(0.d0,ieee_quiet_nan)
  end function trace_debug_nan

  subroutine trace_summary_seed(seed,dL,max_steps,result,b_min)
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_length_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    double precision :: field_min
    integer :: common_status,igrid

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status/=trace_status_active) then
      call trace_summary_init_result(seed,result,common_status)
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_summary_locate_seed(seed,result,igrid)
    if (result%forward_status/=trace_status_active) return
    call trace_summary_trace_seed(igrid,dL,max_steps,field_min,result)
  end subroutine trace_summary_seed

  subroutine trace_summary_twist_seed(seed,dL,max_steps,result,b_min)
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_twist_result), intent(out) :: result
    double precision, intent(in), optional :: b_min

    double precision :: field_min
    integer :: common_status,igrid

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status==trace_status_active .and. ndim/=3) then
      common_status=trace_status_unsupported_geometry
    endif
    if (common_status/=trace_status_active) then
      call trace_summary_init_twist_result(seed,result,common_status)
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)

    call trace_summary_locate_seed(seed,result%line,igrid)
    result%forward_twist=zero
    result%backward_twist=zero
    result%total_twist=zero
    result%valid_twist=.false.
    result%status_twist=result%line%forward_status
    if (result%line%forward_status/=trace_status_active) return
    call trace_summary_trace_twist_seed(igrid,dL,max_steps,field_min,result)
  end subroutine trace_summary_twist_seed

  subroutine trace_summary_mapping_seed(seed,dL,max_steps,result,b_min, &
       source_normal)
    double precision, intent(in) :: seed(ndim),dL
    integer, intent(in) :: max_steps
    type(trace_mapping_result), intent(out) :: result
    double precision, intent(in), optional :: b_min
    double precision, intent(in), optional :: source_normal(3)

    type(trace_summary_state) :: forward_state,backward_state
    type(trace_length_result) :: located
    double precision :: field_min,normal(3)
    integer :: common_status,igrid
    logical :: have_normal

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status==trace_status_active .and. ndim/=3) then
      common_status=trace_status_unsupported_geometry
    endif
    if (common_status/=trace_status_active) then
      call trace_summary_init_mapping_result(seed,result,common_status)
      return
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    normal=zero
    have_normal=present(source_normal)
    if (have_normal) normal=source_normal

    call trace_summary_locate_seed(seed,located,igrid)
    call trace_summary_init_mapping_result(seed,result,located%forward_status)
    if (located%forward_status/=trace_status_active) return

    call trace_summary_init_state(seed,igrid,.true.,1,.false.,forward_state)
    call trace_summary_trace_state_to_end(forward_state,dL,max_steps,field_min)
    call trace_summary_init_state(seed,igrid,.false.,1,.false.,backward_state)
    call trace_summary_trace_state_to_end(backward_state,dL,max_steps,field_min)
    call trace_summary_fill_mapping(seed,igrid,forward_state,backward_state, &
         field_min,have_normal,normal,result)
  end subroutine trace_summary_mapping_seed

  subroutine trace_summary_multi(seeds,nseed,dL,max_steps,results,field_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL,field_min
    type(trace_length_result), intent(out) :: results(nseed)

    type(trace_summary_state), allocatable :: states(:)
    double precision :: seed_local(ndim)
    integer :: common_status,iseed,igrid,iforward,ibackward

    if (nseed<=0) return

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_summary_init_result(seed_local,results(iseed),common_status)
      enddo
      return
    endif

    allocate(states(2*nseed))
    do iseed=1,nseed
      seed_local=seeds(iseed,:)
      call trace_summary_locate_seed(seed_local,results(iseed),igrid)
      iforward=2*iseed-1
      ibackward=2*iseed
      call trace_summary_init_state(seed_local,igrid,.true.,iseed,.false., &
           states(iforward))
      call trace_summary_init_state(seed_local,igrid,.false.,iseed,.false., &
           states(ibackward))
      if (results(iseed)%forward_status/=trace_status_active) then
        states(iforward)%status=results(iseed)%forward_status
        states(iforward)%active=.false.
        states(ibackward)%status=results(iseed)%backward_status
        states(ibackward)%active=.false.
      endif
    enddo

    call trace_summary_trace_states_grouped(states,2*nseed,dL,max_steps, &
         field_min)

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      results(iseed)%forward_footpoint=states(iforward)%footpoint
      results(iseed)%forward_length=states(iforward)%length
      results(iseed)%forward_nstep=states(iforward)%nstep
      results(iseed)%forward_status=states(iforward)%status
      results(iseed)%backward_footpoint=states(ibackward)%footpoint
      results(iseed)%backward_length=states(ibackward)%length
      results(iseed)%backward_nstep=states(ibackward)%nstep
      results(iseed)%backward_status=states(ibackward)%status
      results(iseed)%total_length=results(iseed)%forward_length &
           +results(iseed)%backward_length
    enddo
    deallocate(states)
  end subroutine trace_summary_multi

  subroutine trace_summary_twist_multi(seeds,nseed,dL,max_steps,results, &
       field_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL,field_min
    type(trace_twist_result), intent(out) :: results(nseed)

    type(trace_summary_state), allocatable :: states(:)
    double precision :: seed_local(ndim)
    integer :: common_status,iseed,igrid,iforward,ibackward

    if (nseed<=0) return

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status==trace_status_active .and. ndim/=3) then
      common_status=trace_status_unsupported_geometry
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_summary_init_twist_result(seed_local,results(iseed), &
             common_status)
      enddo
      return
    endif

    allocate(states(2*nseed))
    do iseed=1,nseed
      seed_local=seeds(iseed,:)
      call trace_summary_locate_seed(seed_local,results(iseed)%line,igrid)
      results(iseed)%forward_twist=zero
      results(iseed)%backward_twist=zero
      results(iseed)%total_twist=zero
      results(iseed)%valid_twist=.false.
      results(iseed)%status_twist=results(iseed)%line%forward_status
      iforward=2*iseed-1
      ibackward=2*iseed
      call trace_summary_init_state(seed_local,igrid,.true.,iseed,.true., &
           states(iforward))
      call trace_summary_init_state(seed_local,igrid,.false.,iseed,.true., &
           states(ibackward))
      if (results(iseed)%line%forward_status/=trace_status_active) then
        states(iforward)%status=results(iseed)%line%forward_status
        states(iforward)%active=.false.
        states(ibackward)%status=results(iseed)%line%backward_status
        states(ibackward)%active=.false.
      endif
    enddo

    call trace_summary_trace_states_grouped(states,2*nseed,dL,max_steps, &
         field_min)

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      results(iseed)%line%forward_footpoint=states(iforward)%footpoint
      results(iseed)%line%forward_length=states(iforward)%length
      results(iseed)%line%forward_nstep=states(iforward)%nstep
      results(iseed)%line%forward_status=states(iforward)%status
      results(iseed)%line%backward_footpoint=states(ibackward)%footpoint
      results(iseed)%line%backward_length=states(ibackward)%length
      results(iseed)%line%backward_nstep=states(ibackward)%nstep
      results(iseed)%line%backward_status=states(ibackward)%status
      results(iseed)%line%total_length=results(iseed)%line%forward_length &
           +results(iseed)%line%backward_length
      results(iseed)%forward_twist=states(iforward)%twist
      results(iseed)%backward_twist=states(ibackward)%twist
      results(iseed)%total_twist=results(iseed)%forward_twist &
           +results(iseed)%backward_twist
      results(iseed)%status_twist=trace_summary_twist_status_from_states( &
           states(iforward),states(ibackward))
      results(iseed)%valid_twist=results(iseed)%status_twist &
           ==trace_status_boundary
    enddo
    deallocate(states)
  end subroutine trace_summary_twist_multi

  subroutine trace_summary_mapping_multi(seeds,nseed,dL,max_steps,results, &
       field_min,have_source_normal,source_normal)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL,field_min
    logical, intent(in) :: have_source_normal
    double precision, intent(in) :: source_normal(3)
    type(trace_mapping_result), intent(out) :: results(nseed)

    type(trace_summary_state), allocatable :: states(:)
    integer, allocatable :: source_igrids(:)
    type(trace_length_result) :: located
    double precision :: seed_local(ndim)
    integer :: common_status,iseed,igrid,iforward,ibackward

    if (nseed<=0) return

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status==trace_status_active .and. ndim/=3) then
      common_status=trace_status_unsupported_geometry
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        seed_local=seeds(iseed,:)
        call trace_summary_init_mapping_result(seed_local,results(iseed), &
             common_status)
      enddo
      return
    endif

    allocate(states(2*nseed))
    allocate(source_igrids(nseed))
    source_igrids=-1
    do iseed=1,nseed
      seed_local=seeds(iseed,:)
      call trace_summary_locate_seed(seed_local,located,igrid)
      source_igrids(iseed)=igrid
      call trace_summary_init_mapping_result(seed_local,results(iseed), &
           located%forward_status)
      iforward=2*iseed-1
      ibackward=2*iseed
      call trace_summary_init_state(seed_local,igrid,.true.,iseed,.false., &
           states(iforward))
      call trace_summary_init_state(seed_local,igrid,.false.,iseed,.false., &
           states(ibackward))
      if (located%forward_status/=trace_status_active) then
        states(iforward)%status=located%forward_status
        states(iforward)%active=.false.
        states(ibackward)%status=located%backward_status
        states(ibackward)%active=.false.
      endif
    enddo

    call trace_summary_trace_states_grouped(states,2*nseed,dL,max_steps, &
         field_min)

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      if (source_igrids(iseed)>=0) then
        call trace_summary_fill_mapping(results(iseed)%seed, &
             source_igrids(iseed),states(iforward),states(ibackward), &
             field_min,have_source_normal,source_normal,results(iseed))
      endif
    enddo
    deallocate(source_igrids,states)
  end subroutine trace_summary_mapping_multi

  subroutine trace_summary_topology_multi(seeds,nseed,dL,max_steps,results, &
       field_min,need_twist,need_mapping,have_source_normal,source_normal)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL,field_min
    logical, intent(in) :: need_twist,need_mapping,have_source_normal
    double precision, intent(in) :: source_normal(3)
    type(trace_topology_result), intent(out) :: results(nseed)

    type(trace_summary_state), allocatable :: states(:)
    type(trace_summary_state), allocatable :: map_states(:)
    integer, allocatable :: source_igrids(:)
    type(trace_length_result) :: located
    type(trace_mapping_result) :: mapping
    double precision :: seed_local(ndim)
    integer :: common_status,iseed,igrid,iforward,ibackward,idim

    if (nseed<=0) return

    call trace_summary_validate(dL,max_steps,common_status)
    if (common_status==trace_status_active .and. &
         (need_twist .or. need_mapping) .and. ndim/=3) then
      common_status=trace_status_unsupported_geometry
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        do idim=1,ndim
          seed_local(idim)=seeds(iseed,idim)
        enddo
        call trace_summary_init_topology_result(seed_local,results(iseed), &
             common_status,need_twist,need_mapping)
      enddo
      return
    endif

    allocate(states(2*nseed))
    ! Mapping summaries need only B/Bn; keep them independent of
    ! curl-stencil failures that can stop twist accumulation near boundaries.
    if (need_mapping .and. need_twist) allocate(map_states(2*nseed))
    allocate(source_igrids(nseed))
    source_igrids=-1
    do iseed=1,nseed
      do idim=1,ndim
        seed_local(idim)=seeds(iseed,idim)
      enddo
      call trace_summary_locate_seed(seed_local,located,igrid)
      source_igrids(iseed)=igrid
      call trace_summary_init_topology_result(seed_local,results(iseed), &
           located%forward_status,need_twist,need_mapping)
      iforward=2*iseed-1
      ibackward=2*iseed
      call trace_summary_init_state(seed_local,igrid,.true.,iseed, &
           need_twist,states(iforward))
      call trace_summary_init_state(seed_local,igrid,.false.,iseed, &
           need_twist,states(ibackward))
      if (allocated(map_states)) then
        call trace_summary_init_state(seed_local,igrid,.true.,iseed, &
             .false.,map_states(iforward))
        call trace_summary_init_state(seed_local,igrid,.false.,iseed, &
             .false.,map_states(ibackward))
      endif
      if (located%forward_status/=trace_status_active) then
        states(iforward)%status=located%forward_status
        states(iforward)%active=.false.
        states(ibackward)%status=located%backward_status
        states(ibackward)%active=.false.
        if (allocated(map_states)) then
          map_states(iforward)%status=located%forward_status
          map_states(iforward)%active=.false.
          map_states(ibackward)%status=located%backward_status
          map_states(ibackward)%active=.false.
        endif
      endif
    enddo

    call trace_summary_trace_states_grouped(states,2*nseed,dL,max_steps, &
         field_min)
    if (allocated(map_states)) then
      call trace_summary_trace_states_grouped(map_states,2*nseed,dL, &
           max_steps,field_min)
    endif

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      results(iseed)%forward_endpoint=states(iforward)%footpoint
      results(iseed)%backward_endpoint=states(ibackward)%footpoint
      results(iseed)%length_forward=states(iforward)%length
      results(iseed)%length_backward=states(ibackward)%length
      results(iseed)%length_total=results(iseed)%length_forward &
           +results(iseed)%length_backward
      results(iseed)%forward_nstep=states(iforward)%nstep
      results(iseed)%backward_nstep=states(ibackward)%nstep
      results(iseed)%forward_face=states(iforward)%face
      results(iseed)%backward_face=states(ibackward)%face
      results(iseed)%forward_status=states(iforward)%status
      results(iseed)%backward_status=states(ibackward)%status
      if (need_twist) then
        results(iseed)%twist_forward=states(iforward)%twist
        results(iseed)%twist_backward=states(ibackward)%twist
        results(iseed)%twist_total=results(iseed)%twist_forward &
             +results(iseed)%twist_backward
        results(iseed)%status_twist=trace_summary_twist_status_from_states( &
             states(iforward),states(ibackward))
        results(iseed)%valid_twist=results(iseed)%status_twist &
             ==trace_status_boundary
      endif

      if (states(iforward)%status==trace_status_boundary .and. &
           states(ibackward)%status==trace_status_boundary) then
        results(iseed)%status=trace_status_boundary
        results(iseed)%valid=.true.
      else
        results(iseed)%valid=.false.
        if (states(iforward)%status/=trace_status_boundary) then
          results(iseed)%status=states(iforward)%status
        else
          results(iseed)%status=states(ibackward)%status
        endif
      endif

      if (need_mapping .and. source_igrids(iseed)>=0) then
        if (allocated(map_states)) then
          call trace_summary_fill_mapping(results(iseed)%seed, &
               source_igrids(iseed),map_states(iforward), &
               map_states(ibackward),field_min,have_source_normal, &
               source_normal,mapping)
        else
          call trace_summary_fill_mapping(results(iseed)%seed, &
               source_igrids(iseed),states(iforward),states(ibackward), &
               field_min,have_source_normal,source_normal,mapping)
        endif
        results(iseed)%source_B=mapping%source_B
        results(iseed)%forward_B=mapping%forward_B
        results(iseed)%backward_B=mapping%backward_B
        results(iseed)%source_Bn=mapping%source_Bn
        results(iseed)%forward_Bn=mapping%forward_Bn
        results(iseed)%backward_Bn=mapping%backward_Bn
        results(iseed)%map_forward_endpoint=mapping%forward_footpoint
        results(iseed)%map_backward_endpoint=mapping%backward_footpoint
        results(iseed)%map_forward_length=mapping%forward_length
        results(iseed)%map_backward_length=mapping%backward_length
        results(iseed)%map_forward_face=mapping%forward_face
        results(iseed)%map_backward_face=mapping%backward_face
        results(iseed)%map_forward_status=mapping%forward_status
        results(iseed)%map_backward_status=mapping%backward_status
        results(iseed)%forward_status=mapping%forward_status
        results(iseed)%backward_status=mapping%backward_status
        results(iseed)%valid=mapping%valid
        if (mapping%valid) then
          results(iseed)%status=trace_status_boundary
        else if (mapping%forward_status/=trace_status_boundary) then
          results(iseed)%status=mapping%forward_status
        else
          results(iseed)%status=mapping%backward_status
        endif
      endif
    enddo
    if (allocated(map_states)) deallocate(map_states)
    deallocate(source_igrids,states)
  end subroutine trace_summary_topology_multi

  subroutine trace_summary_validate(dL,max_steps,status)
    double precision, intent(in) :: dL
    integer, intent(in) :: max_steps
    integer, intent(out) :: status

    status=trace_status_active
    if (npe/=1) then
      status=trace_status_mpi_unsupported
    else if (dL<=zero .or. max_steps<0) then
      status=trace_status_invalid_input
    else if (.not.trace_is_supported_summary_geometry()) then
      status=trace_status_unsupported_geometry
    endif
  end subroutine trace_summary_validate

  logical function trace_is_supported_summary_geometry() result(is_supported)
    is_supported=.false.
    if (ndim/=3) return
    {^IFTHREED
    select case (geo_coordinate)
    case (geo_cartesian,geo_cartesian_stretched)
      is_supported=.true.
    case (geo_spherical)
      ! Phase-1 spherical tracing uses native coordinates and assumes a
      ! non-periodic phi interval. Periodic wrapping is a later extension.
      is_supported=.not.periodB(3)
    end select
    }
  end function trace_is_supported_summary_geometry

  subroutine trace_summary_init_result(seed,result,status)
    double precision, intent(in) :: seed(ndim)
    type(trace_length_result), intent(out) :: result
    integer, intent(in) :: status

    result%seed=seed
    result%forward_footpoint=seed
    result%backward_footpoint=seed
    result%forward_length=zero
    result%backward_length=zero
    result%total_length=zero
    result%forward_nstep=0
    result%backward_nstep=0
    result%forward_status=status
    result%backward_status=status
  end subroutine trace_summary_init_result

  subroutine trace_summary_init_twist_result(seed,result,status)
    double precision, intent(in) :: seed(ndim)
    type(trace_twist_result), intent(out) :: result
    integer, intent(in) :: status

    call trace_summary_init_result(seed,result%line,status)
    result%forward_twist=zero
    result%backward_twist=zero
    result%total_twist=zero
    result%valid_twist=.false.
    result%status_twist=status
  end subroutine trace_summary_init_twist_result

  subroutine trace_summary_init_mapping_result(seed,result,status)
    double precision, intent(in) :: seed(ndim)
    type(trace_mapping_result), intent(out) :: result
    integer, intent(in) :: status

    result%seed=seed
    result%source_B=zero
    result%forward_footpoint=seed
    result%backward_footpoint=seed
    result%forward_B=zero
    result%backward_B=zero
    result%forward_length=zero
    result%backward_length=zero
    result%source_Bn=zero
    result%forward_Bn=zero
    result%backward_Bn=zero
    result%forward_face=trace_face_none
    result%backward_face=trace_face_none
    result%forward_status=status
    result%backward_status=status
    result%valid=.false.
  end subroutine trace_summary_init_mapping_result

  subroutine trace_summary_init_topology_result(seed,result,status, &
       need_twist,need_mapping)
    double precision, intent(in) :: seed(ndim)
    type(trace_topology_result), intent(out) :: result
    integer, intent(in) :: status
    logical, intent(in) :: need_twist,need_mapping

    result%seed=seed
    result%length_forward=zero
    result%length_backward=zero
    result%length_total=zero
    result%twist_forward=zero
    result%twist_backward=zero
    result%twist_total=zero
    result%forward_endpoint=seed
    result%backward_endpoint=seed
    result%forward_nstep=0
    result%backward_nstep=0
    result%forward_face=trace_face_none
    result%backward_face=trace_face_none
    result%forward_status=status
    result%backward_status=status
    result%map_forward_endpoint=seed
    result%map_backward_endpoint=seed
    result%map_forward_length=zero
    result%map_backward_length=zero
    result%map_forward_face=trace_face_none
    result%map_backward_face=trace_face_none
    result%map_forward_status=status
    result%map_backward_status=status
    result%source_B=zero
    result%forward_B=zero
    result%backward_B=zero
    result%source_Bn=zero
    result%forward_Bn=zero
    result%backward_Bn=zero
    result%has_twist=need_twist
    result%has_mapping=need_mapping
    result%valid_twist=.false.
    result%valid=.false.
    result%status_twist=status
    result%status=status
  end subroutine trace_summary_init_topology_result

  subroutine trace_summary_locate_seed(seed,result,igrid)
    use mod_particle_base, only: find_particle_ipe

    double precision, intent(in) :: seed(ndim)
    type(trace_length_result), intent(out) :: result
    integer, intent(out) :: igrid

    double precision :: x3d(3),domain_min(ndim),domain_max(ndim),tol
    integer :: indomain,ipe,j

    call trace_summary_init_result(seed,result,trace_status_active)
    igrid=-1
    indomain=0
    ^D&domain_min(^D)=xprobmin^D;
    ^D&domain_max(^D)=xprobmax^D;
    tol=100.d0*epsilon(one)*max(one,maxval(abs(domain_max-domain_min)))
    do j=1,ndim
      if (seed(j)>=domain_min(j)-tol .and. &
           seed(j)<=domain_max(j)+tol) indomain=indomain+1
    enddo
    if (indomain/=ndim) then
      result%forward_status=trace_status_seed_outside
      result%backward_status=trace_status_seed_outside
      return
    endif

    x3d=zero
    do j=1,ndim
      x3d(j)=min(max(seed(j),domain_min(j)+tol),domain_max(j)-tol)
    enddo
    call find_particle_ipe(x3d,igrid,ipe)
    if (igrid<0 .or. ipe/=mype) then
      result%forward_status=trace_status_out_of_domain
      result%backward_status=trace_status_out_of_domain
      return
    endif
  end subroutine trace_summary_locate_seed

  subroutine trace_summary_trace_seed(igrid,dL,max_steps,field_min,result)
    integer, intent(in) :: igrid,max_steps
    double precision, intent(in) :: dL,field_min
    type(trace_length_result), intent(inout) :: result

    type(trace_summary_state) :: forward_state,backward_state

    call trace_summary_init_state(result%seed,igrid,.true.,1,.false.,forward_state)
    call trace_summary_trace_state_to_end(forward_state,dL,max_steps,field_min)
    call trace_summary_init_state(result%seed,igrid,.false.,1,.false.,backward_state)
    call trace_summary_trace_state_to_end(backward_state,dL,max_steps,field_min)

    result%forward_footpoint=forward_state%footpoint
    result%forward_length=forward_state%length
    result%forward_nstep=forward_state%nstep
    result%forward_status=forward_state%status
    result%backward_footpoint=backward_state%footpoint
    result%backward_length=backward_state%length
    result%backward_nstep=backward_state%nstep
    result%backward_status=backward_state%status
    result%total_length=result%forward_length+result%backward_length
  end subroutine trace_summary_trace_seed

  subroutine trace_summary_trace_twist_seed(igrid,dL,max_steps,field_min, &
       result)
    integer, intent(in) :: igrid,max_steps
    double precision, intent(in) :: dL,field_min
    type(trace_twist_result), intent(inout) :: result

    type(trace_summary_state) :: forward_state,backward_state

    call trace_summary_init_state(result%line%seed,igrid,.true.,1,.true., &
         forward_state)
    call trace_summary_trace_state_to_end(forward_state,dL,max_steps,field_min)
    call trace_summary_init_state(result%line%seed,igrid,.false.,1,.true., &
         backward_state)
    call trace_summary_trace_state_to_end(backward_state,dL,max_steps,field_min)

    result%line%forward_footpoint=forward_state%footpoint
    result%line%forward_length=forward_state%length
    result%line%forward_nstep=forward_state%nstep
    result%line%forward_status=forward_state%status
    result%line%backward_footpoint=backward_state%footpoint
    result%line%backward_length=backward_state%length
    result%line%backward_nstep=backward_state%nstep
    result%line%backward_status=backward_state%status
    result%line%total_length=result%line%forward_length &
         +result%line%backward_length
    result%forward_twist=forward_state%twist
    result%backward_twist=backward_state%twist
    result%total_twist=result%forward_twist+result%backward_twist
    result%status_twist=trace_summary_twist_status_from_states( &
         forward_state,backward_state)
    result%valid_twist=result%status_twist==trace_status_boundary
  end subroutine trace_summary_trace_twist_seed

  subroutine trace_summary_init_state(xseed,igrid,forward,seed_id, &
       accumulate_twist,state)
    double precision, intent(in) :: xseed(ndim)
    integer, intent(in) :: igrid,seed_id
    logical, intent(in) :: forward,accumulate_twist
    type(trace_summary_state), intent(out) :: state

    state%x=xseed
    state%footpoint=xseed
    state%length=zero
    state%twist=zero
    state%nstep=0
    state%status=trace_status_active
    state%twist_status=trace_status_active
    state%igrid=igrid
    state%seed_id=seed_id
    state%face=trace_face_none
    state%forward=forward
    state%active=.true.
    state%accumulate_twist=accumulate_twist
    state%rk45_h=zero
    state%sph_cache%valid=.false.
    if (trace_rk45_position_integrator()) then
      call trace_rk45_stats_note_direction()
    endif
    if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
         call trace_spherical_profile_add_count(trace_profile_directions,1_8)
  end subroutine trace_summary_init_state

  subroutine trace_summary_trace_state_to_end(state,ds,max_nstep,threshold)
    type(trace_summary_state), intent(inout) :: state
    double precision, intent(in) :: ds,threshold
    integer, intent(in) :: max_nstep

    do while(state%active .and. state%nstep<max_nstep)
      call trace_summary_advance_state(state,ds,threshold)
    enddo

    if (state%active) then
      state%status=trace_status_max_steps
      state%active=.false.
    endif
    state%footpoint=state%x
    if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
         call trace_spherical_profile_note_trace_steps(state%nstep)
  end subroutine trace_summary_trace_state_to_end

  subroutine trace_summary_trace_states_grouped(states,nstate,ds,max_nstep,threshold)
    integer, intent(in) :: nstate,max_nstep
    type(trace_summary_state), intent(inout) :: states(nstate)
    double precision, intent(in) :: ds,threshold

    logical, allocatable :: processed(:)
    integer :: istate,jstate,target_grid
    logical :: any_active

    allocate(processed(nstate))
    do
      any_active=.false.
      do istate=1,nstate
        if (states(istate)%active) then
          any_active=.true.
          exit
        endif
      enddo
      if (.not.any_active) exit

      processed=.false.
      do istate=1,nstate
        if (.not.states(istate)%active .or. processed(istate)) cycle
        target_grid=states(istate)%igrid
        do jstate=1,nstate
          if (.not.states(jstate)%active .or. processed(jstate)) cycle
          if (states(jstate)%igrid/=target_grid) cycle
          call trace_summary_advance_state_in_grid(states(jstate),target_grid, &
               ds,max_nstep,threshold)
          processed(jstate)=.true.
        enddo
      enddo
    enddo
    deallocate(processed)
  end subroutine trace_summary_trace_states_grouped

  subroutine trace_summary_advance_state_in_grid(state,igrid,ds,max_nstep,threshold)
    type(trace_summary_state), intent(inout) :: state
    integer, intent(in) :: igrid,max_nstep
    double precision, intent(in) :: ds,threshold

    do while(state%active .and. state%igrid==igrid .and. &
         state%nstep<max_nstep)
      call trace_summary_advance_state(state,ds,threshold)
    enddo

    if (state%active .and. state%nstep>=max_nstep) then
      state%status=trace_status_max_steps
      state%active=.false.
      state%footpoint=state%x
      if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
           call trace_spherical_profile_note_trace_steps(state%nstep)
    endif
  end subroutine trace_summary_advance_state_in_grid

  subroutine trace_summary_advance_state(state,ds,threshold)
    type(trace_summary_state), intent(inout) :: state
    double precision, intent(in) :: ds,threshold

    double precision :: xnext(ndim),xhit(ndim),ds_actual,segment_length
    double precision :: alpha_hit
    double precision :: xbmin^D,xbmax^D
    double precision :: domain_min(ndim),domain_max(ndim)
    double precision :: near_tol,best_dist,face_dist
    integer :: igrid_next,ipe_next,inblock,step_status,face_id
    integer :: idim,best_dim,best_side
    integer :: point_domain
    double precision :: rk45_twist_increment
    integer :: rk45_twist_status
    logical :: newpe,stop_trace,hit_ok
    logical :: rk45_step,rk45_twist_integrated

    if (.not.state%active) return

    rk45_step=.false.
    rk45_twist_increment=zero
    rk45_twist_status=trace_status_active
    rk45_twist_integrated=.false.
    if (trace_integrator_mode==trace_integrator_rk45_cartesian .and. &
         trace_cartesian_like_geometry()) then
      rk45_step=.true.
      call trace_summary_rk45_cartesian_step(state%x,state%igrid,ds, &
           state%forward,threshold,state%rk45_h,xnext,ds_actual, &
           step_status,state%accumulate_twist,rk45_twist_increment, &
           rk45_twist_status,rk45_twist_integrated)
    else if (trace_integrator_mode==trace_integrator_rk45_spherical .and. &
         geo_coordinate==geo_spherical) then
      rk45_step=.true.
      call trace_summary_rk45_spherical_step(state%x,state%igrid,ds, &
           state%forward,threshold,state%rk45_h,state%sph_cache,xnext, &
           ds_actual,step_status,state%accumulate_twist, &
           rk45_twist_increment,rk45_twist_status,rk45_twist_integrated)
    else
      call trace_summary_rk2_step(state%x,state%igrid,ds,state%forward, &
           threshold,xnext,ds_actual,step_status,sph_cache=state%sph_cache)
    endif
    if (step_status/=trace_status_active) then
      if (step_status==trace_status_out_of_domain .and. &
           trace_rk45_position_integrator()) then
        call trace_boundary_face_at_point(state%x,state%face,hit_ok)
        if (.not.hit_ok) then
          ^D&domain_min(^D)=xprobmin^D;
          ^D&domain_max(^D)=xprobmax^D;
          near_tol=max(100.d0*epsilon(one)* &
               max(one,maxval(abs(domain_max-domain_min))), &
               1.d-5*max(one,abs(ds)))
          best_dist=huge(one)
          best_dim=0
          best_side=0
          do idim=1,ndim
            face_dist=abs(state%x(idim)-domain_min(idim))
            if (face_dist<best_dist) then
              best_dist=face_dist
              best_dim=idim
              best_side=-1
            endif
            face_dist=abs(state%x(idim)-domain_max(idim))
            if (face_dist<best_dist) then
              best_dist=face_dist
              best_dim=idim
              best_side=1
            endif
          enddo
          if (best_dim>0 .and. best_dist<=near_tol) then
            if (best_side<0) then
              state%x(best_dim)=domain_min(best_dim)
            else
              state%x(best_dim)=domain_max(best_dim)
            endif
            state%footpoint=state%x
            state%face=trace_face_from_dim_side(best_dim,best_side)
            hit_ok=.true.
          endif
        endif
        if (hit_ok) then
          state%status=trace_status_boundary
          state%active=.false.
          return
        endif
      endif
      state%status=step_status
      state%active=.false.
      return
    endif

    point_domain=0
    {if (xnext(^DB)>=xprobmin^DB .and. xnext(^DB)<xprobmax^DB) point_domain=point_domain+1\}
    if (point_domain/=ndim) then
      call trace_intersect_domain(state%x,xnext,xhit,hit_ok,face_id, &
           alpha_hit)
      if (hit_ok) then
        if (trace_rk45_position_integrator() .and. &
             rk45_step) then
          segment_length=abs(ds_actual)
        else
          segment_length=trace_segment_length(state%x,xhit,ds_actual, &
               alpha_hit)
        endif
        if (state%accumulate_twist .and. &
             state%twist_status==trace_status_active) then
          if (rk45_step .and. rk45_twist_integrated) then
            state%twist=state%twist+rk45_twist_increment
          else if (rk45_step .and. rk45_twist_status/=trace_status_active) then
            state%twist_status=rk45_twist_status
          else
            call trace_summary_accumulate_twist(state,xhit,segment_length, &
                 threshold,step_status)
            if (step_status/=trace_status_active) then
              state%twist_status=step_status
            endif
          endif
        endif
        state%length=state%length+segment_length
        ! Count the nonzero terminal boundary segment as one completed step.
        state%nstep=state%nstep+1
        if (trace_spherical_profile_enabled .and. &
             geo_coordinate==geo_spherical) then
          call trace_spherical_profile_add_count(trace_profile_steps,1_8)
        endif
        state%x=xhit
        state%footpoint=state%x
        state%face=face_id
      endif
      if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
           call trace_spherical_profile_add_count(trace_profile_boundary_events, &
           1_8)
      state%status=trace_status_boundary
      state%active=.false.
      if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
           call trace_spherical_profile_note_trace_steps(state%nstep)
      return
    endif

    if (state%accumulate_twist .and. &
         state%twist_status==trace_status_active) then
      if (rk45_step .and. rk45_twist_integrated) then
        state%twist=state%twist+rk45_twist_increment
      else if (rk45_step .and. rk45_twist_status/=trace_status_active) then
        state%twist_status=rk45_twist_status
      else
        call trace_summary_accumulate_twist(state,xnext,ds_actual,threshold, &
             step_status)
        if (step_status/=trace_status_active) then
          state%twist_status=step_status
        endif
      endif
    endif
    if (rk45_step) then
      state%length=state%length+abs(ds_actual)
    else
      state%length=state%length+trace_segment_length(state%x,xnext, &
           ds_actual)
    endif
    state%nstep=state%nstep+1
    if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
         call trace_spherical_profile_add_count(trace_profile_steps,1_8)
    state%x=xnext
    state%footpoint=state%x

    ^D&xbmin^D=rnode(rpxmin^D_,state%igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,state%igrid)\
    inblock=0
    {if (state%x(^DB)>=xbmin^DB .and. state%x(^DB)<xbmax^DB) inblock=inblock+1\}
    if (inblock==ndim) return

    igrid_next=state%igrid
    ipe_next=mype
    newpe=.false.
    stop_trace=.false.
    call find_next_grid(state%igrid,igrid_next,ipe_next,state%x,newpe,stop_trace)
    if (stop_trace) then
      state%status=trace_status_boundary
      state%active=.false.
      if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
           call trace_spherical_profile_add_count(trace_profile_boundary_events, &
           1_8)
      if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
           call trace_spherical_profile_note_trace_steps(state%nstep)
      return
    endif
    if (newpe .or. ipe_next/=mype) then
      if (trace_rk45_position_integrator()) then
        call trace_boundary_face_at_point(state%x,state%face,hit_ok)
        if (.not.hit_ok) then
          ^D&domain_min(^D)=xprobmin^D;
          ^D&domain_max(^D)=xprobmax^D;
          near_tol=max(100.d0*epsilon(one)* &
               max(one,maxval(abs(domain_max-domain_min))), &
               1.d-6*max(one,abs(ds_actual)))
          best_dist=huge(one)
          best_dim=0
          best_side=0
          do idim=1,ndim
            face_dist=abs(state%x(idim)-domain_min(idim))
            if (face_dist<best_dist) then
              best_dist=face_dist
              best_dim=idim
              best_side=-1
            endif
            face_dist=abs(state%x(idim)-domain_max(idim))
            if (face_dist<best_dist) then
              best_dist=face_dist
              best_dim=idim
              best_side=1
            endif
          enddo
          if (best_dim>0 .and. best_dist<=near_tol) then
            if (best_side<0) then
              state%x(best_dim)=domain_min(best_dim)
            else
              state%x(best_dim)=domain_max(best_dim)
            endif
            state%footpoint=state%x
            state%face=trace_face_from_dim_side(best_dim,best_side)
            hit_ok=.true.
          endif
        endif
        if (hit_ok) then
          state%status=trace_status_boundary
          state%active=.false.
          return
        endif
      endif
      state%status=trace_status_out_of_domain
      state%active=.false.
      if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
           call trace_spherical_profile_note_trace_steps(state%nstep)
      return
    endif
    state%igrid=igrid_next
  end subroutine trace_summary_advance_state

  subroutine trace_summary_accumulate_twist(state,xend,segment_length, &
       threshold,status)
    type(trace_summary_state), intent(inout) :: state
    double precision, intent(in) :: xend(ndim),segment_length,threshold
    integer, intent(out) :: status

    double precision :: xmid(ndim),B(3),curlB(3),B2,twist_density

    status=trace_status_active
    if (segment_length<=zero) return

    xmid=half*(state%x+xend)
    call sample_B_curlB_at_point(xmid,state%igrid,threshold,B,curlB,status, &
         sph_cache=state%sph_cache)
    if (status/=trace_status_active) return

    B2=sum(B**2)
    twist_density=sum(curlB*B)/(4.d0*dpi*B2)
    ! Both directions contribute with positive arc length to the full-line Tw.
    state%twist=state%twist+twist_density*segment_length
  end subroutine trace_summary_accumulate_twist

  subroutine trace_twist_density_at_point(x,igrid,threshold,twist_density, &
       status,sph_cache)
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: twist_density
    integer, intent(out) :: status
    type(trace_sph_interp_ctx), intent(inout), optional :: sph_cache

    double precision :: B(3),curlB(3),B2

    twist_density=zero
    call sample_B_curlB_at_point(x,igrid,threshold,B,curlB,status, &
         sph_cache=sph_cache)
    if (status/=trace_status_active) return

    B2=sum(B**2)
    if (B2<=zero) then
      status=trace_status_weak_field
      return
    endif
    twist_density=sum(curlB*B)/(4.d0*dpi*B2)
  end subroutine trace_twist_density_at_point

  subroutine trace_tangent_accumulate_twist(state,xend,segment_length, &
       threshold)
    type(trace_tangent_state), intent(inout) :: state
    double precision, intent(in) :: xend(ndim),segment_length,threshold

    double precision :: xmid(ndim),twist_density

    if (.not.state%accumulate_twist) return
    if (state%twist_status/=trace_status_active) return
    if (segment_length<=zero) return
    xmid=half*(state%x+xend)
    call trace_twist_density_at_point(xmid,state%igrid,threshold, &
         twist_density,state%twist_status,sph_cache=state%sph_cache)
    if (state%twist_status/=trace_status_active) return
    state%twist=state%twist+twist_density*segment_length
  end subroutine trace_tangent_accumulate_twist

  integer function trace_summary_twist_status_from_states(forward_state, &
       backward_state) result(status_twist)
    type(trace_summary_state), intent(in) :: forward_state,backward_state

    status_twist=trace_status_boundary
    if (forward_state%status/=trace_status_boundary) then
      status_twist=forward_state%status
    else if (backward_state%status/=trace_status_boundary) then
      status_twist=backward_state%status
    else if (forward_state%twist_status/=trace_status_active) then
      status_twist=forward_state%twist_status
    else if (backward_state%twist_status/=trace_status_active) then
      status_twist=backward_state%twist_status
    endif
  end function trace_summary_twist_status_from_states

  subroutine trace_tangent_fill_twist_result(seed,forward_state, &
       backward_state,result)
    double precision, intent(in) :: seed(ndim)
    type(trace_tangent_state), intent(in) :: forward_state,backward_state
    type(trace_twist_result), intent(out) :: result

    result%line%seed=seed
    result%line%forward_footpoint=forward_state%endpoint
    result%line%backward_footpoint=backward_state%endpoint
    result%line%forward_length=forward_state%length
    result%line%backward_length=backward_state%length
    result%line%total_length=forward_state%length+backward_state%length
    result%line%forward_nstep=forward_state%nstep
    result%line%backward_nstep=backward_state%nstep
    result%line%forward_status=forward_state%status
    result%line%backward_status=backward_state%status
    result%forward_twist=forward_state%twist
    result%backward_twist=backward_state%twist
    result%total_twist=forward_state%twist+backward_state%twist
    result%status_twist=trace_status_boundary
    if (forward_state%status/=trace_status_boundary) then
      result%status_twist=forward_state%status
    else if (backward_state%status/=trace_status_boundary) then
      result%status_twist=backward_state%status
    else if (forward_state%twist_status/=trace_status_active) then
      result%status_twist=forward_state%twist_status
    else if (backward_state%twist_status/=trace_status_active) then
      result%status_twist=backward_state%twist_status
    endif
    result%valid_twist=result%status_twist==trace_status_boundary
  end subroutine trace_tangent_fill_twist_result

  subroutine trace_total_B_at_cell(igrid,ix1,ix2,ix3,B,status)
    integer, intent(in) :: igrid,ix1,ix2,ix3
    double precision, intent(out) :: B(3)
    integer, intent(out) :: status

    B=zero
    status=trace_status_out_of_domain
    if (ndim/=3 .or. igrid<0) return

    {^IFTHREED
    if (.not.B0field .and. .not.allocated(iw_mag)) return
    if (allocated(iw_mag)) then
      B(1:3)=ps(igrid)%w(ix1,ix2,ix3,iw_mag(1:3))
    endif
    if (B0field) then
      B(1:3)=B(1:3)+ps(igrid)%B0(ix1,ix2,ix3,1:3,0)
    endif
    status=trace_status_active
    }
  end subroutine trace_total_B_at_cell

  subroutine trace_bhat_at_cell(igrid,ix1,ix2,ix3,threshold,bhat,status)
    integer, intent(in) :: igrid,ix1,ix2,ix3
    double precision, intent(in) :: threshold
    double precision, intent(out) :: bhat(3)
    integer, intent(out) :: status

    double precision :: B(3),Bnorm

    bhat=zero
    call trace_total_B_at_cell(igrid,ix1,ix2,ix3,B,status)
    if (status/=trace_status_active) return

    Bnorm=dsqrt(sum(B**2))
    if (Bnorm<=zero .or. Bnorm<threshold) then
      status=trace_status_weak_field
      return
    endif
    bhat=B/Bnorm
  end subroutine trace_bhat_at_cell

  subroutine sample_bhat_gradbhat_at_point(x,igrid,threshold,bhat, &
       grad_bhat,status)
    ! Production Method-II sampler: trilinear total-B derivative + chain rule.
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: bhat(3),grad_bhat(3,3)
    integer, intent(out) :: status

    call sample_bhat_gradbhat_interpderiv_at_point(x,igrid,threshold,bhat, &
         grad_bhat,status)
  end subroutine sample_bhat_gradbhat_at_point

  subroutine sample_bhat_gradbhat_cellfd_at_point(x,igrid,threshold,bhat, &
       grad_bhat,status)
    ! Legacy diagnostic sampler: interpolate cell-centered finite-difference grad(bhat).
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: bhat(3),grad_bhat(3,3)
    integer, intent(out) :: status

    double precision :: dxb^D,dxc^D,xd^D
    double precision :: bhat_cell(0:1^D&,3)
    double precision :: grad_cell(0:1^D&,3,3)
    double precision :: factor(0:1^D&)
    double precision :: bhat_center(3),bhat_plus(3),bhat_minus(3)
    integer :: ixI^L,ixbl^D,ix^D,j,k

    bhat=zero
    grad_bhat=zero
    status=trace_status_bad_grad_stencil
    if (ndim/=3 .or. igrid<0) return

    {^IFTHREED
    ixI^L=ixG^LL;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    ^D&ixbl^D=floor((x(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
    ^D&xd^D=(x(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D;

    if (ixbl1-1<ixImin1 .or. ixbl1+2>ixImax1 .or. &
         ixbl2-1<ixImin2 .or. ixbl2+2>ixImax2 .or. &
         ixbl3-1<ixImin3 .or. ixbl3+2>ixImax3) return

    bhat_cell=zero
    grad_cell=zero
    {do ix^D=0,1\}
      call trace_bhat_at_cell(igrid,ixbl1+ix1,ixbl2+ix2, &
           ixbl3+ix3,threshold,bhat_center,status)
      if (status/=trace_status_active) return
      bhat_cell(ix^D,1:3)=bhat_center

      call trace_bhat_at_cell(igrid,ixbl1+ix1+1,ixbl2+ix2, &
           ixbl3+ix3,threshold,bhat_plus,status)
      if (status/=trace_status_active) return
      call trace_bhat_at_cell(igrid,ixbl1+ix1-1,ixbl2+ix2, &
           ixbl3+ix3,threshold,bhat_minus,status)
      if (status/=trace_status_active) return
      grad_cell(ix^D,1:3,1)=(bhat_plus-bhat_minus)/(2.d0*dxb1)

      call trace_bhat_at_cell(igrid,ixbl1+ix1,ixbl2+ix2+1, &
           ixbl3+ix3,threshold,bhat_plus,status)
      if (status/=trace_status_active) return
      call trace_bhat_at_cell(igrid,ixbl1+ix1,ixbl2+ix2-1, &
           ixbl3+ix3,threshold,bhat_minus,status)
      if (status/=trace_status_active) return
      grad_cell(ix^D,1:3,2)=(bhat_plus-bhat_minus)/(2.d0*dxb2)

      call trace_bhat_at_cell(igrid,ixbl1+ix1,ixbl2+ix2, &
           ixbl3+ix3+1,threshold,bhat_plus,status)
      if (status/=trace_status_active) return
      call trace_bhat_at_cell(igrid,ixbl1+ix1,ixbl2+ix2, &
           ixbl3+ix3-1,threshold,bhat_minus,status)
      if (status/=trace_status_active) return
      grad_cell(ix^D,1:3,3)=(bhat_plus-bhat_minus)/(2.d0*dxb3)
    {enddo\}

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}
    {do ix^D=0,1\}
      do j=1,3
        bhat(j)=bhat(j)+bhat_cell(ix^D,j)*factor(ix^D)
        do k=1,3
          grad_bhat(j,k)=grad_bhat(j,k)+grad_cell(ix^D,j,k)*factor(ix^D)
        enddo
      enddo
    {enddo\}

    status=trace_status_active
    }
  end subroutine sample_bhat_gradbhat_cellfd_at_point

  subroutine sample_bhat_gradbhat_xeps_at_point(x,eps,threshold,bhat, &
       grad_bhat,status)
    ! Reference-only grad(bhat) by central differences of interpolated bhat.
    double precision, intent(in) :: x(ndim),eps,threshold
    double precision, intent(out) :: bhat(3),grad_bhat(3,3)
    integer, intent(out) :: status

    double precision :: xp(ndim),xm(ndim)
    double precision :: B(3),Bp(3),Bm(3),bhat_p(3),bhat_m(3)
    double precision :: Bnorm,Bpnorm,Bmnorm
    integer :: igrid,idir,point_domain

    bhat=zero
    grad_bhat=zero
    status=trace_status_bad_grad_stencil
    if (ndim/=3 .or. .not.slab_uniform) then
      status=trace_status_unsupported_geometry
      return
    endif
    if (eps<=zero) then
      status=trace_status_invalid_input
      return
    endif

    call trace_debug_locate_point(x,igrid,status)
    if (status/=trace_status_active) return
    call sample_B_at_point(x,igrid,threshold,B,status)
    if (status/=trace_status_active) return
    Bnorm=dsqrt(sum(B**2))
    if (Bnorm<=zero .or. Bnorm<threshold) then
      status=trace_status_weak_field
      return
    endif
    bhat=B/Bnorm

    do idir=1,ndim
      xp=x
      xm=x
      xp(idir)=xp(idir)+eps
      xm(idir)=xm(idir)-eps

      point_domain=0
      {if (xp(^DB)>=xprobmin^DB .and. xp(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        status=trace_status_bad_grad_stencil
        return
      endif
      call trace_debug_locate_point(xp,igrid,status)
      if (status/=trace_status_active) return
      call sample_B_at_point(xp,igrid,threshold,Bp,status)
      if (status/=trace_status_active) return
      Bpnorm=dsqrt(sum(Bp**2))
      if (Bpnorm<=zero .or. Bpnorm<threshold) then
        status=trace_status_weak_field
        return
      endif
      bhat_p=Bp/Bpnorm

      point_domain=0
      {if (xm(^DB)>=xprobmin^DB .and. xm(^DB)<xprobmax^DB) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        status=trace_status_bad_grad_stencil
        return
      endif
      call trace_debug_locate_point(xm,igrid,status)
      if (status/=trace_status_active) return
      call sample_B_at_point(xm,igrid,threshold,Bm,status)
      if (status/=trace_status_active) return
      Bmnorm=dsqrt(sum(Bm**2))
      if (Bmnorm<=zero .or. Bmnorm<threshold) then
        status=trace_status_weak_field
        return
      endif
      bhat_m=Bm/Bmnorm

      grad_bhat(1:3,idir)=(bhat_p-bhat_m)/(2.d0*eps)
    enddo

    status=trace_status_active
  end subroutine sample_bhat_gradbhat_xeps_at_point

  subroutine sample_bhat_gradbhat_interpderiv_at_point(x,igrid,threshold, &
       bhat,grad_bhat,status)
    ! Candidate sampler: differentiate trilinear total-B interpolation.
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: bhat(3),grad_bhat(3,3)
    integer, intent(out) :: status

    double precision :: dxb^D,dxc^D,xd^D
    double precision :: field(0:1^D&,3)
    double precision :: B(3),dBdx(3,3)
    double precision :: wx(0:1),wy(0:1),wz(0:1)
    double precision :: dwx(0:1),dwy(0:1),dwz(0:1)
    double precision :: weight,Bnorm,projected
    integer :: ixI^L,ixbl^D,ix^D,j

    bhat=zero
    grad_bhat=zero
    status=trace_status_bad_grad_stencil
    if (ndim/=3 .or. .not.trace_cartesian_like_geometry() .or. igrid<0) then
      status=trace_status_unsupported_geometry
      return
    endif

    {^IFTHREED
    ixI^L=ixG^LL;
    call trace_interp_weights_block(x,igrid,ixI^L,ixbl^D,xd^D,dxc^D, &
         status)
    if (status/=trace_status_active) return

    if (ixbl1<ixImin1 .or. ixbl1+1>ixImax1 .or. &
         ixbl2<ixImin2 .or. ixbl2+1>ixImax2 .or. &
         ixbl3<ixImin3 .or. ixbl3+1>ixImax3) return

    status=trace_status_out_of_domain
    if (.not.B0field .and. .not.allocated(iw_mag)) return

    field=zero
    {do ix^D=0,1\}
      if (allocated(iw_mag)) then
        field(ix^D,1:3)=ps(igrid)%w(ixbl1+ix1,ixbl2+ix2, &
             ixbl3+ix3,iw_mag(1:3))
      endif
      if (B0field) then
        field(ix^D,1:3)=field(ix^D,1:3) &
             +ps(igrid)%B0(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3,1:3,0)
      endif
    {enddo\}

    wx(0)=one-xd1
    wx(1)=xd1
    wy(0)=one-xd2
    wy(1)=xd2
    wz(0)=one-xd3
    wz(1)=xd3
    dwx(0)=-one/dxc1
    dwx(1)= one/dxc1
    dwy(0)=-one/dxc2
    dwy(1)= one/dxc2
    dwz(0)=-one/dxc3
    dwz(1)= one/dxc3

    B=zero
    dBdx=zero
    {do ix^D=0,1\}
      weight=wx(ix1)*wy(ix2)*wz(ix3)
      B(1:3)=B(1:3)+field(ix^D,1:3)*weight
      dBdx(1:3,1)=dBdx(1:3,1)+field(ix^D,1:3)* &
           dwx(ix1)*wy(ix2)*wz(ix3)
      dBdx(1:3,2)=dBdx(1:3,2)+field(ix^D,1:3)* &
           wx(ix1)*dwy(ix2)*wz(ix3)
      dBdx(1:3,3)=dBdx(1:3,3)+field(ix^D,1:3)* &
           wx(ix1)*wy(ix2)*dwz(ix3)
    {enddo\}

    Bnorm=dsqrt(sum(B**2))
    if (Bnorm<=zero .or. Bnorm<threshold) then
      status=trace_status_weak_field
      return
    endif

    bhat=B/Bnorm
    do j=1,3
      projected=sum(bhat*dBdx(1:3,j))
      grad_bhat(1:3,j)=(dBdx(1:3,j)-bhat(1:3)*projected)/Bnorm
    enddo

    status=trace_status_active
    }
  end subroutine sample_bhat_gradbhat_interpderiv_at_point

  subroutine sample_B_curlB_at_point(x,igrid,threshold,B,curlB,status, &
       sph_cache)
    ! Interpolate total B and its cell-centered second-order curl locally.
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: B(3),curlB(3)
    integer, intent(out) :: status
    type(trace_sph_interp_ctx), intent(inout), optional :: sph_cache

    double precision :: dxb^D,dxc^D,xd^D
    double precision :: field(0:1^D&,3),current(0:1^D&,3)
    double precision :: factor(0:1^D&),B2
    double precision :: dx1c,dx2c,dx3c
    integer :: ixI^L,ixJ^L,ixbl^D,ix^D,j

    B=zero
    curlB=zero
    status=trace_status_bad_curl_stencil
    if (ndim/=3) return

    if (geo_coordinate==geo_spherical) then
      if (present(sph_cache)) then
        call sample_B_curlB_spherical_at_point(x,igrid,threshold,B,curlB, &
             status,sph_cache=sph_cache)
      else
        call sample_B_curlB_spherical_at_point(x,igrid,threshold,B,curlB, &
             status)
      endif
      return
    endif

    {^IFTHREED
    ixI^L=ixG^LL;
    ixJ^L=ixM^LL^LADD1;
    call trace_interp_weights_block(x,igrid,ixI^L,ixbl^D,xd^D,dxc^D, &
         status)
    if (status/=trace_status_active) return

    if (ixbl1-1<ixImin1 .or. ixbl1+2>ixImax1 .or. &
         ixbl2-1<ixImin2 .or. ixbl2+2>ixImax2 .or. &
         ixbl3-1<ixImin3 .or. ixbl3+2>ixImax3) return
    if (B0field) then
      if (ixbl1<ixJmin1 .or. ixbl1+1>ixJmax1 .or. &
           ixbl2<ixJmin2 .or. ixbl2+1>ixJmax2 .or. &
           ixbl3<ixJmin3 .or. ixbl3+1>ixJmax3) return
    endif
    status=trace_status_out_of_domain
    if (.not.B0field .and. .not.allocated(iw_mag)) return

    field=zero
    current=zero
    {do ix^D=0,1\}
      if (allocated(iw_mag)) then
        field(ix^D,1:3)=ps(igrid)%w(ixbl1+ix1,ixbl2+ix2, &
             ixbl3+ix3,iw_mag(1:3))
      endif
      if (B0field) then
        field(ix^D,1:3)=field(ix^D,1:3) &
             +ps(igrid)%B0(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3,1:3,0)
      endif
      if (allocated(iw_mag)) then
        dx1c=ps(igrid)%x(ixbl1+ix1+1,ixbl2+ix2,ixbl3+ix3,1) &
             -ps(igrid)%x(ixbl1+ix1-1,ixbl2+ix2,ixbl3+ix3,1)
        dx2c=ps(igrid)%x(ixbl1+ix1,ixbl2+ix2+1,ixbl3+ix3,2) &
             -ps(igrid)%x(ixbl1+ix1,ixbl2+ix2-1,ixbl3+ix3,2)
        dx3c=ps(igrid)%x(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3+1,3) &
             -ps(igrid)%x(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3-1,3)
        if (abs(dx1c)<=smalldouble .or. abs(dx2c)<=smalldouble .or. &
             abs(dx3c)<=smalldouble) then
          status=trace_status_bad_curl_stencil
          return
        endif
        current(ix^D,1)=half*( &
             (ps(igrid)%w(ixbl1+ix1,ixbl2+ix2+1,ixbl3+ix3,iw_mag(3)) &
             -ps(igrid)%w(ixbl1+ix1,ixbl2+ix2-1,ixbl3+ix3,iw_mag(3)))/dx2c &
             -(ps(igrid)%w(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3+1,iw_mag(2)) &
             -ps(igrid)%w(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3-1,iw_mag(2)))/dx3c)
        current(ix^D,2)=half*( &
             (ps(igrid)%w(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3+1,iw_mag(1)) &
             -ps(igrid)%w(ixbl1+ix1,ixbl2+ix2,ixbl3+ix3-1,iw_mag(1)))/dx3c &
             -(ps(igrid)%w(ixbl1+ix1+1,ixbl2+ix2,ixbl3+ix3,iw_mag(3)) &
             -ps(igrid)%w(ixbl1+ix1-1,ixbl2+ix2,ixbl3+ix3,iw_mag(3)))/dx1c)
        current(ix^D,3)=half*( &
             (ps(igrid)%w(ixbl1+ix1+1,ixbl2+ix2,ixbl3+ix3,iw_mag(2)) &
             -ps(igrid)%w(ixbl1+ix1-1,ixbl2+ix2,ixbl3+ix3,iw_mag(2)))/dx1c &
             -(ps(igrid)%w(ixbl1+ix1,ixbl2+ix2+1,ixbl3+ix3,iw_mag(1)) &
             -ps(igrid)%w(ixbl1+ix1,ixbl2+ix2-1,ixbl3+ix3,iw_mag(1)))/dx2c)
      endif
      if (B0field) then
        current(ix^D,1:3)=current(ix^D,1:3) &
             +ps(igrid)%J0(ixbl^D+ix^D,1:3)
      endif
    {enddo\}

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}
    {do ix^D=0,1\}
      do j=1,3
        B(j)=B(j)+field(ix^D,j)*factor(ix^D)
        curlB(j)=curlB(j)+current(ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

    B2=sum(B**2)
    if (B2<=zero .or. dsqrt(B2)<threshold) then
      status=trace_status_weak_field
      return
    endif
    status=trace_status_active
    }
  end subroutine sample_B_curlB_at_point

  subroutine sample_B_curlB_spherical_at_point(x,igrid,threshold,B,curlB, &
       status,sph_cache)
    ! Reuse AMRVAC's geometry-aware curl operator for local spherical
    ! physical components. B0/J0 handling matches the MHD current path.
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: B(3),curlB(3)
    integer, intent(out) :: status
    type(trace_sph_interp_ctx), intent(inout), optional :: sph_cache

    double precision :: dxb^D,dxc^D,xd^D
    double precision :: bvec(ixG^T,1:3),current(ixG^T,1:3)
    double precision :: factor(0:1^D&),B2
    type(trace_sph_interp_ctx) :: ctx
    integer :: ixI^L,ixO^L,ixA^L,ixJ^L,ixbl^D,ix^D,j,idirmin

    B=zero
    curlB=zero
    status=trace_status_bad_curl_stencil
    if (ndim/=3 .or. geo_coordinate/=geo_spherical .or. igrid<0) return

    {^IFTHREED
    if (present(sph_cache)) then
      call trace_spherical_interp_ctx_build_cached(x,igrid,sph_cache,ctx, &
           status)
    else
      call trace_spherical_interp_ctx_build(x,igrid,ctx,status)
    endif
    if (status/=trace_status_active) return
    call trace_spherical_sample_Bsph_ctx(ctx,threshold,B,status)
    if (status/=trace_status_active) return
    if (trace_spherical_curl_cache_ready) then
      call trace_spherical_sample_cached_curlB_ctx(ctx,curlB,status)
      if (status/=trace_status_active) return
      B2=sum(B**2)
      if (B2<=zero .or. dsqrt(B2)<threshold) then
        status=trace_status_weak_field
        return
      endif
      status=trace_status_active
      return
    endif

    ixI^L=ixG^LL;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    ^D&dxlevel(^D)=dxb^D;
    ^D&dxc^D=ctx%dxc(^D);
    ^D&xd^D=ctx%xd(^D);
    ^D&ixbl^D=ctx%ixbl(^D);
    ^D&ixOmin^D=ixbl^D;
    ^D&ixOmax^D=ixbl^D+1;
    ixA^L=ixO^L^LADD1;
    ixJ^L=ixM^LL^LADD1;

    if (ixImin1>ixAmin1 .or. ixImax1<ixAmax1 .or. &
         ixImin2>ixAmin2 .or. ixImax2<ixAmax2 .or. &
         ixImin3>ixAmin3 .or. ixImax3<ixAmax3) return
    if (B0field) then
      if (ixJmin1>ixOmin1 .or. ixJmax1<ixOmax1 .or. &
           ixJmin2>ixOmin2 .or. ixJmax2<ixOmax2 .or. &
           ixJmin3>ixOmin3 .or. ixJmax3<ixOmax3) return
    endif
    status=trace_status_out_of_domain
    if (.not.B0field .and. .not.allocated(iw_mag)) return

    bvec=zero
    current=zero
    if (allocated(iw_mag)) then
      bvec(ixA^S,1:3)=ps(igrid)%w(ixA^S,iw_mag(1:3))
    endif

    block=>ps(igrid)
    call curlvector(bvec,ixI^L,ixO^L,current,idirmin,1,3)
    if (B0field) then
      current(ixO^S,1:3)=current(ixO^S,1:3)+ps(igrid)%J0(ixO^S,1:3)
    endif

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}
    {do ix^D=0,1\}
      do j=1,3
        curlB(j)=curlB(j)+current(ixbl^D+ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

    B2=sum(B**2)
    if (B2<=zero .or. dsqrt(B2)<threshold) then
      status=trace_status_weak_field
      return
    endif
    status=trace_status_active
    }
  end subroutine sample_B_curlB_spherical_at_point

  subroutine sample_B_at_point(x,igrid,threshold,B,status)
    ! Interpolate total magnetic field at x using the local block stencil.
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid
    double precision, intent(out) :: B(3)
    integer, intent(out) :: status

    double precision :: dxb^D,dxc^D,xd^D
    double precision :: field(0:1^D&,3),factor(0:1^D&),B2
    double precision :: Bcell(3)
    type(trace_sph_interp_ctx) :: ctx
    integer :: ixI^L,ixbl^D,ix^D,j

    B=zero
    status=trace_status_out_of_domain
    if (ndim/=3 .or. igrid<0) return

    {^IFTHREED
    ixI^L=ixG^LL;
    if (geo_coordinate==geo_spherical) then
      call trace_spherical_interp_ctx_build(x,igrid,ctx,status)
      if (status/=trace_status_active) return
      call trace_spherical_sample_Bsph_ctx(ctx,threshold,B,status)
      return
    else
      call trace_interp_weights_block(x,igrid,ixI^L,ixbl^D,xd^D,dxc^D, &
           status)
      if (status/=trace_status_active) return
    endif

    if (ixbl1<ixImin1 .or. ixbl1+1>ixImax1 .or. &
         ixbl2<ixImin2 .or. ixbl2+1>ixImax2 .or. &
         ixbl3<ixImin3 .or. ixbl3+1>ixImax3) return
    field=zero
    {do ix^D=0,1\}
      call trace_total_B_at_cell(igrid,ixbl1+ix1,ixbl2+ix2, &
           ixbl3+ix3,Bcell,status)
      if (status/=trace_status_active) return
      field(ix^D,1:3)=Bcell
    {enddo\}

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}
    {do ix^D=0,1\}
      do j=1,3
        B(j)=B(j)+field(ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

    B2=sum(B**2)
    if (B2<=zero .or. dsqrt(B2)<threshold) then
      status=trace_status_weak_field
      return
    endif
    status=trace_status_active
    }
  end subroutine sample_B_at_point

  subroutine trace_summary_fill_mapping(seed,source_igrid,forward_state, &
       backward_state,threshold,have_source_normal,source_normal,result)
    double precision, intent(in) :: seed(ndim),threshold,source_normal(3)
    integer, intent(in) :: source_igrid
    logical, intent(in) :: have_source_normal
    type(trace_summary_state), intent(in) :: forward_state,backward_state
    type(trace_mapping_result), intent(inout) :: result

    integer :: status

    call trace_summary_init_mapping_result(seed,result,trace_status_active)

    result%forward_footpoint=forward_state%footpoint
    result%backward_footpoint=backward_state%footpoint
    result%forward_length=forward_state%length
    result%backward_length=backward_state%length
    result%forward_face=forward_state%face
    result%backward_face=backward_state%face
    result%forward_status=forward_state%status
    result%backward_status=backward_state%status

    call sample_B_at_point(seed,source_igrid,threshold,result%source_B,status)
    if (status/=trace_status_active) then
      result%forward_status=status
      result%backward_status=status
      result%valid=.false.
      return
    endif
    if (have_source_normal) then
      result%source_Bn=sum(result%source_B*source_normal)
    endif

    if (forward_state%status==trace_status_boundary) then
      call sample_B_at_point(forward_state%footpoint,forward_state%igrid, &
           threshold,result%forward_B,status)
      if (status/=trace_status_active) result%forward_status=status
    endif
    if (backward_state%status==trace_status_boundary) then
      call sample_B_at_point(backward_state%footpoint,backward_state%igrid, &
           threshold,result%backward_B,status)
      if (status/=trace_status_active) result%backward_status=status
    endif

    result%forward_Bn=trace_face_Bn(result%forward_face,result%forward_B)
    result%backward_Bn=trace_face_Bn(result%backward_face,result%backward_B)
    result%valid=result%forward_status==trace_status_boundary .and. &
         result%backward_status==trace_status_boundary .and. &
         trace_face_is_boundary(result%forward_face) .and. &
         trace_face_is_boundary(result%backward_face)
  end subroutine trace_summary_fill_mapping

  double precision function trace_face_Bn(face_id,B) result(Bn)
    integer, intent(in) :: face_id
    double precision, intent(in) :: B(3)

    Bn=zero
    select case(face_id)
    case(trace_face_xmin)
      Bn=-B(1)
    case(trace_face_xmax)
      Bn=B(1)
    case(trace_face_ymin)
      Bn=-B(2)
    case(trace_face_ymax)
      Bn=B(2)
    case(trace_face_zmin)
      Bn=-B(3)
    case(trace_face_zmax)
      Bn=B(3)
    end select
  end function trace_face_Bn

  logical function trace_face_is_boundary(face_id) result(is_boundary)
    integer, intent(in) :: face_id

    is_boundary=face_id>=trace_face_xmin .and. face_id<=trace_face_zmax
  end function trace_face_is_boundary

  subroutine trace_boundary_face_at_point(x,face_id,on_boundary)
    double precision, intent(in) :: x(ndim)
    integer, intent(out) :: face_id
    logical, intent(out) :: on_boundary

    double precision :: domain_min(ndim),domain_max(ndim),tol
    integer :: idim,nface,candidate_face

    ^D&domain_min(^D)=xprobmin^D;
    ^D&domain_max(^D)=xprobmax^D;
    tol=100.d0*epsilon(one)*max(one,maxval(abs(domain_max-domain_min)))
    face_id=trace_face_none
    nface=0
    do idim=1,ndim
      if (abs(x(idim)-domain_min(idim))<=tol) then
        candidate_face=trace_face_from_dim_side(idim,-1)
        nface=nface+1
        if (nface==1) then
          face_id=candidate_face
        else if (candidate_face/=face_id) then
          face_id=trace_face_ambiguous
        endif
      endif
      if (abs(x(idim)-domain_max(idim))<=tol) then
        candidate_face=trace_face_from_dim_side(idim,1)
        nface=nface+1
        if (nface==1) then
          face_id=candidate_face
        else if (candidate_face/=face_id) then
          face_id=trace_face_ambiguous
        endif
      endif
    enddo
    on_boundary=nface>0
  end subroutine trace_boundary_face_at_point

  subroutine trace_intersect_domain(xinside,xoutside,xhit,hit_ok,face_id, &
       alpha_out)
    double precision, intent(in) :: xinside(ndim),xoutside(ndim)
    double precision, intent(out) :: xhit(ndim)
    logical, intent(out) :: hit_ok
    integer, intent(out) :: face_id
    double precision, intent(out), optional :: alpha_out

    double precision :: alpha,alpha_hit,alpha_tol,delta
    double precision :: domain_min(ndim),domain_max(ndim)
    integer :: idim,hit_dim,hit_side,candidate_face

    ^D&domain_min(^D)=xprobmin^D;
    ^D&domain_max(^D)=xprobmax^D;
    alpha_hit=huge(one)
    alpha_tol=100.d0*epsilon(one)
    hit_dim=0
    hit_side=0
    face_id=trace_face_none
    do idim=1,ndim
      delta=xoutside(idim)-xinside(idim)
      if (xoutside(idim)<domain_min(idim) .and. abs(delta)>smalldouble) then
        alpha=(domain_min(idim)-xinside(idim))/delta
        if (alpha>=zero .and. alpha<=one) then
          candidate_face=trace_face_from_dim_side(idim,-1)
          if (alpha<alpha_hit-alpha_tol) then
            alpha_hit=alpha
            hit_dim=idim
            hit_side=-1
            face_id=candidate_face
          else if (abs(alpha-alpha_hit)<=alpha_tol .and. &
               candidate_face/=face_id) then
            face_id=trace_face_ambiguous
          endif
        endif
      else if (xoutside(idim)>=domain_max(idim) .and. &
           abs(delta)>smalldouble) then
        alpha=(domain_max(idim)-xinside(idim))/delta
        if (alpha>=zero .and. alpha<=one) then
          candidate_face=trace_face_from_dim_side(idim,1)
          if (alpha<alpha_hit-alpha_tol) then
            alpha_hit=alpha
            hit_dim=idim
            hit_side=1
            face_id=candidate_face
          else if (abs(alpha-alpha_hit)<=alpha_tol .and. &
               candidate_face/=face_id) then
            face_id=trace_face_ambiguous
          endif
        endif
      endif
    enddo

    hit_ok=hit_dim>0
    xhit=xinside
    if (present(alpha_out)) alpha_out=zero
    if (.not.hit_ok) return

    xhit=xinside+alpha_hit*(xoutside-xinside)
    do idim=1,ndim
      xhit(idim)=min(max(xhit(idim),domain_min(idim)),domain_max(idim))
    enddo
    if (hit_side<0) then
      xhit(hit_dim)=domain_min(hit_dim)
    else
      xhit(hit_dim)=domain_max(hit_dim)
    endif
    if (present(alpha_out)) alpha_out=alpha_hit
  end subroutine trace_intersect_domain

  integer function trace_face_from_dim_side(idim,side) result(face_id)
    integer, intent(in) :: idim,side

    face_id=trace_face_none
    select case(idim)
    case(1)
      if (side<0) then
        face_id=trace_face_xmin
      else
        face_id=trace_face_xmax
      endif
    case(2)
      if (side<0) then
        face_id=trace_face_ymin
      else
        face_id=trace_face_ymax
      endif
    case(3)
      if (side<0) then
        face_id=trace_face_zmin
      else
        face_id=trace_face_zmax
      endif
    end select
  end function trace_face_from_dim_side

  subroutine trace_cartesian_local_cell_size(x,igrid,hcell,status)
    double precision, intent(in) :: x(ndim)
    integer, intent(in) :: igrid
    double precision, intent(out) :: hcell
    integer, intent(out) :: status

    double precision :: dxc^D,xd^D
    integer :: ixI^L,ixbl^D

    hcell=zero
    status=trace_status_unsupported_geometry
    if (.not.trace_cartesian_like_geometry() .or. igrid<0) return

    {^IFTHREED
    ixI^L=ixG^LL;
    call trace_interp_weights_block(x,igrid,ixI^L,ixbl^D,xd^D,dxc^D, &
         status)
    if (status/=trace_status_active) return
    hcell=min(abs(dxc1),abs(dxc2))
    hcell=min(hcell,abs(dxc3))
    if (hcell<=zero) status=trace_status_out_of_domain
    }
  end subroutine trace_cartesian_local_cell_size

  subroutine trace_cartesian_rhs_bhat(x,igrid_hint,threshold,bhat,igrid, &
       status)
    double precision, intent(in) :: x(ndim),threshold
    integer, intent(in) :: igrid_hint
    double precision, intent(out) :: bhat(ndim)
    integer, intent(out) :: igrid,status

    double precision :: dxc^D,xd^D
    double precision :: field(0:1^D&,3),factor(0:1^D&)
    double precision :: Bcell(3),B(3),Bnorm
    integer :: ixI^L,ixbl^D,ix^D,j

    bhat=zero
    igrid=-1
    status=trace_status_unsupported_geometry
    if (.not.trace_cartesian_like_geometry()) return

    {^IFTHREED
    call trace_locate_point_with_hint(x,igrid_hint,igrid,status)
    if (status/=trace_status_active) return

    ixI^L=ixG^LL;
    call trace_interp_weights_block(x,igrid,ixI^L,ixbl^D,xd^D,dxc^D, &
         status)
    if (status/=trace_status_active) return

    field=zero
    {do ix^D=0,1\}
      call trace_total_B_at_cell(igrid,ixbl1+ix1,ixbl2+ix2, &
           ixbl3+ix3,Bcell,status)
      if (status/=trace_status_active) return
      field(ix^D,1:3)=Bcell
    {enddo\}

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}

    B=zero
    {do ix^D=0,1\}
      do j=1,3
        B(j)=B(j)+field(ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

    Bnorm=dsqrt(sum(B**2))
    if (Bnorm<=zero .or. Bnorm<threshold) then
      status=trace_status_weak_field
      return
    endif
    bhat=B(1:ndim)/Bnorm
    status=trace_status_active
    }
  end subroutine trace_cartesian_rhs_bhat

  subroutine trace_rk45_try_boundary_finish(xnow,xstage,hfloor,xnext, &
       ds_actual,h_suggest,status,finished)
    double precision, intent(in) :: xnow(ndim),xstage(ndim),hfloor
    double precision, intent(out) :: xnext(ndim),ds_actual
    double precision, intent(inout) :: h_suggest
    integer, intent(out) :: status
    logical, intent(out) :: finished

    double precision :: xhit(ndim)
    integer :: face_id
    logical :: hit_ok

    finished=.false.
    call trace_intersect_domain(xnow,xstage,xhit,hit_ok,face_id)
    if (.not.hit_ok) return

    xnext=xhit
    ds_actual=dsqrt(sum((xnext-xnow)**2))
    h_suggest=hfloor
    status=trace_status_active
    call trace_rk45_stats_note_attempt(.true.,.true.,ds_actual)
    finished=.true.
  end subroutine trace_rk45_try_boundary_finish

  subroutine trace_summary_rk45_cartesian_step(xnow,igrid,ds,forward, &
       threshold,h_suggest,xnext,ds_actual,status,accumulate_twist, &
       twist_increment,twist_status,twist_integrated)
    double precision, intent(in) :: xnow(ndim),ds,threshold
    integer, intent(in) :: igrid
    logical, intent(in) :: forward
    double precision, intent(inout) :: h_suggest
    double precision, intent(out) :: xnext(ndim),ds_actual
    integer, intent(out) :: status
    logical, intent(in) :: accumulate_twist
    double precision, intent(out) :: twist_increment
    integer, intent(out) :: twist_status
    logical, intent(out) :: twist_integrated

    double precision, parameter :: b21=1.d0/5.d0
    double precision, parameter :: b31=3.d0/40.d0,b32=9.d0/40.d0
    double precision, parameter :: b41=3.d0/10.d0,b42=-9.d0/10.d0, &
         b43=6.d0/5.d0
    double precision, parameter :: b51=-11.d0/54.d0,b52=5.d0/2.d0, &
         b53=-70.d0/27.d0,b54=35.d0/27.d0
    double precision, parameter :: b61=1631.d0/55296.d0, &
         b62=175.d0/512.d0,b63=575.d0/13824.d0, &
         b64=44275.d0/110592.d0,b65=253.d0/4096.d0
    double precision, parameter :: c1=37.d0/378.d0,c3=250.d0/621.d0, &
         c4=125.d0/594.d0,c6=512.d0/1771.d0
    double precision, parameter :: cs1=2825.d0/27648.d0, &
         cs3=18575.d0/48384.d0,cs4=13525.d0/55296.d0, &
         cs5=277.d0/14336.d0,cs6=one/4.d0

    double precision :: k1(ndim),k2(ndim),k3(ndim),k4(ndim),k5(ndim),k6(ndim)
    double precision :: x2(ndim),x3(ndim),x4(ndim),x5s(ndim),x6(ndim)
    double precision :: x5(ndim),x4err(ndim),xhit(ndim)
    double precision :: h,hmax,hfloor,hcell,tol,err,grow,sgn
    double precision :: tw1,tw3,tw4,tw6
    double precision :: domain_min(ndim),domain_max(ndim)
    integer :: igrid1,igrid2,igrid3,igrid4,igrid5,igrid6
    integer :: iter,point_domain,face_id
    integer :: reject_reason
    logical :: hit_ok,boundary_limited,boundary_finished

    xnext=xnow
    ds_actual=zero
    status=trace_status_active
    twist_increment=zero
    twist_status=trace_status_active
    twist_integrated=.false.
    if (.not.trace_cartesian_like_geometry()) then
      status=trace_status_unsupported_geometry
      return
    endif

    hmax=abs(ds)
    if (hmax<=zero) then
      status=trace_status_invalid_input
      return
    endif
    h=h_suggest
    if (h<=zero) h=hmax
    h=min(h,hmax)
    hfloor=max(trace_step_min,100.d0*epsilon(one)*max(one,hmax))
    h=max(h,hfloor)
    sgn=one
    if (.not.forward) sgn=-one
    ^D&domain_min(^D)=xprobmin^D;
    ^D&domain_max(^D)=xprobmax^D;

    do iter=1,100
      h=max(hfloor,min(h,hmax))
      call trace_cartesian_rhs_bhat(xnow,igrid,threshold,k1,igrid1,status)
      if (status/=trace_status_active) return

      x2=xnow+sgn*h*b21*k1
      call trace_cartesian_rhs_bhat(x2,igrid1,threshold,k2,igrid2,status)
      if (status/=trace_status_active) then
        if (status==trace_status_weak_field) return
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) then
          call trace_rk45_try_boundary_finish(xnow,x2,hfloor,xnext, &
               ds_actual,h_suggest,status,boundary_finished)
          if (boundary_finished) return
        endif
        reject_reason=trace_rk45_reject_stage_failure
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif

      x3=xnow+sgn*h*(b31*k1+b32*k2)
      call trace_cartesian_rhs_bhat(x3,igrid2,threshold,k3,igrid3,status)
      if (status/=trace_status_active) then
        if (status==trace_status_weak_field) return
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) then
          call trace_rk45_try_boundary_finish(xnow,x3,hfloor,xnext, &
               ds_actual,h_suggest,status,boundary_finished)
          if (boundary_finished) return
        endif
        reject_reason=trace_rk45_reject_stage_failure
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif

      x4=xnow+sgn*h*(b41*k1+b42*k2+b43*k3)
      call trace_cartesian_rhs_bhat(x4,igrid3,threshold,k4,igrid4,status)
      if (status/=trace_status_active) then
        if (status==trace_status_weak_field) return
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) then
          call trace_rk45_try_boundary_finish(xnow,x4,hfloor,xnext, &
               ds_actual,h_suggest,status,boundary_finished)
          if (boundary_finished) return
        endif
        reject_reason=trace_rk45_reject_stage_failure
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif

      x5s=xnow+sgn*h*(b51*k1+b52*k2+b53*k3+b54*k4)
      call trace_cartesian_rhs_bhat(x5s,igrid4,threshold,k5,igrid5,status)
      if (status/=trace_status_active) then
        if (status==trace_status_weak_field) return
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) then
          call trace_rk45_try_boundary_finish(xnow,x5s,hfloor,xnext, &
               ds_actual,h_suggest,status,boundary_finished)
          if (boundary_finished) return
        endif
        reject_reason=trace_rk45_reject_stage_failure
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif

      x6=xnow+sgn*h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5)
      call trace_cartesian_rhs_bhat(x6,igrid5,threshold,k6,igrid6,status)
      if (status/=trace_status_active) then
        if (status==trace_status_weak_field) return
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) then
          call trace_rk45_try_boundary_finish(xnow,x6,hfloor,xnext, &
               ds_actual,h_suggest,status,boundary_finished)
          if (boundary_finished) return
        endif
        reject_reason=trace_rk45_reject_stage_failure
        if (status==trace_status_boundary .or. &
             status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif

      x5=xnow+sgn*h*(c1*k1+c3*k3+c4*k4+c6*k6)
      x4err=xnow+sgn*h*(cs1*k1+cs3*k3+cs4*k4+cs5*k5+cs6*k6)

      point_domain=0
      {if (x5(^DB)>=domain_min(^DB) .and. x5(^DB)<domain_max(^DB)) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        call trace_intersect_domain(xnow,x5,xhit,hit_ok,face_id)
        if (hit_ok) then
          xnext=xhit
          ds_actual=dsqrt(sum((xnext-xnow)**2))
          h_suggest=hfloor
          status=trace_status_active
          call trace_rk45_stats_note_attempt(.true.,.true.,ds_actual)
          return
        endif
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             trace_rk45_reject_boundary)
        status=trace_status_out_of_domain
        return
      endif

      call trace_cartesian_local_cell_size(xnow,igrid1,hcell,status)
      if (status/=trace_status_active .or. hcell<=zero) hcell=hmax
      err=dsqrt(sum((x5-x4err)**2))
      tol=trace_rk45_atol+trace_rk45_rtol*max(h,hcell)
      if (err<=tol .or. h<=hfloor*(one+epsilon(one))) then
        xnext=x5
        ds_actual=h
        if (accumulate_twist) then
          call trace_twist_density_at_point(xnow,igrid1,threshold,tw1, &
               twist_status)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x3,igrid3,threshold,tw3, &
               twist_status)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x4,igrid4,threshold,tw4, &
               twist_status)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x6,igrid6,threshold,tw6, &
               twist_status)
          if (twist_status==trace_status_active) then
            twist_increment=h*(c1*tw1+c3*tw3+c4*tw4+c6*tw6)
            twist_integrated=.true.
          endif
        endif
        if (err>zero .and. tol>zero) then
          grow=trace_rk45_safety*(tol/err)**0.2d0
          grow=max(trace_rk45_min_shrink,min(trace_rk45_max_grow,grow))
        else
          grow=trace_rk45_max_grow
        endif
        h_suggest=min(hmax,max(hfloor,h*grow))
        call trace_rk45_stats_note_attempt(.true.,.false.,h)
        status=trace_status_active
        return
      endif

      grow=trace_rk45_safety*(tol/max(err,smalldouble))**0.25d0
      grow=max(trace_rk45_min_shrink,min(one,grow))
      call trace_rk45_stats_note_attempt(.false.,.false.,h, &
           trace_rk45_reject_error)
      h=max(hfloor,h*grow)
    enddo

    call trace_intersect_domain(xnow,xnow+sgn*h*k1,xhit,hit_ok,face_id)
    if (hit_ok) then
      boundary_limited=.true.
      xnext=xhit
      ds_actual=dsqrt(sum((xnext-xnow)**2))
      h_suggest=hfloor
      call trace_rk45_stats_note_attempt(.true.,boundary_limited,ds_actual)
      status=trace_status_active
    else
      status=trace_status_out_of_domain
    endif
  end subroutine trace_summary_rk45_cartesian_step

  subroutine trace_summary_rk45_spherical_step(xnow,igrid,ds,forward, &
       threshold,h_suggest,sph_cache,xnext,ds_actual,status, &
       accumulate_twist,twist_increment,twist_status,twist_integrated)
    double precision, intent(in) :: xnow(ndim),ds,threshold
    integer, intent(in) :: igrid
    logical, intent(in) :: forward
    double precision, intent(inout) :: h_suggest
    type(trace_sph_interp_ctx), intent(inout) :: sph_cache
    double precision, intent(out) :: xnext(ndim),ds_actual
    integer, intent(out) :: status
    logical, intent(in) :: accumulate_twist
    double precision, intent(out) :: twist_increment
    integer, intent(out) :: twist_status
    logical, intent(out) :: twist_integrated

    double precision, parameter :: b21=1.d0/5.d0
    double precision, parameter :: b31=3.d0/40.d0,b32=9.d0/40.d0
    double precision, parameter :: b41=3.d0/10.d0,b42=-9.d0/10.d0, &
         b43=6.d0/5.d0
    double precision, parameter :: b51=-11.d0/54.d0,b52=5.d0/2.d0, &
         b53=-70.d0/27.d0,b54=35.d0/27.d0
    double precision, parameter :: b61=1631.d0/55296.d0, &
         b62=175.d0/512.d0,b63=575.d0/13824.d0, &
         b64=44275.d0/110592.d0,b65=253.d0/4096.d0
    double precision, parameter :: c1=37.d0/378.d0,c3=250.d0/621.d0, &
         c4=125.d0/594.d0,c6=512.d0/1771.d0
    double precision, parameter :: cs1=2825.d0/27648.d0, &
         cs3=18575.d0/48384.d0,cs4=13525.d0/55296.d0, &
         cs5=277.d0/14336.d0,cs6=one/4.d0

    double precision :: k1(ndim),k2(ndim),k3(ndim),k4(ndim),k5(ndim),k6(ndim)
    double precision :: x2(ndim),x3(ndim),x4(ndim),x5s(ndim),x6(ndim)
    double precision :: x5(ndim),x4err(ndim),xhit(ndim)
    double precision :: h,hmax,hfloor,hcell,tol,err,grow,sgn
    double precision :: tw1,tw3,tw4,tw6
    double precision :: r_metric,sin_theta_metric,alpha_hit
    double precision :: domain_min(ndim),domain_max(ndim)
    type(trace_sph_interp_ctx) :: ctx1,ctx2,ctx3,ctx4,ctx5,ctx6
    integer :: igrid1,igrid2,igrid3,igrid4,igrid5,igrid6
    integer :: iter,point_domain,face_id,ctx_status
    integer :: reject_reason
    logical :: hit_ok,field_ok
    character(len=std_len) :: field_type

    xnext=xnow
    ds_actual=zero
    status=trace_status_active
    twist_increment=zero
    twist_status=trace_status_active
    twist_integrated=.false.
    if (geo_coordinate/=geo_spherical) then
      status=trace_status_unsupported_geometry
      return
    endif

    hmax=abs(ds)
    if (hmax<=zero) then
      status=trace_status_invalid_input
      return
    endif
    h=h_suggest
    if (h<=zero) h=hmax
    h=min(h,hmax)
    hfloor=max(trace_step_min,100.d0*epsilon(one)*max(one,hmax))
    h=max(h,hfloor)
    sgn=one
    if (.not.forward) sgn=-one
    field_type='Bfield'
    ^D&domain_min(^D)=xprobmin^D;
    ^D&domain_max(^D)=xprobmax^D;

    do iter=1,100
      h=max(hfloor,min(h,hmax))
      call trace_spherical_interp_ctx_build_cached(xnow,igrid,sph_cache, &
           ctx1,ctx_status)
      if (ctx_status/=trace_status_active) then
        status=ctx_status
        return
      endif
      call get_K_spherical_ctx(ctx1,k1,field_type,threshold,field_ok)
      if (.not.field_ok) then
        status=trace_status_weak_field
        return
      endif
      igrid1=ctx1%igrid
      hcell=ctx1%h_local

      x2=xnow+sgn*h*b21*k1
      call trace_spherical_interp_ctx_build_cached(x2,igrid1,sph_cache, &
           ctx2,ctx_status)
      if (ctx_status/=trace_status_active) then
        reject_reason=trace_rk45_reject_stage_failure
        if (ctx_status==trace_status_boundary .or. &
             ctx_status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif
      call get_K_spherical_ctx(ctx2,k2,field_type,threshold,field_ok)
      if (.not.field_ok) then
        status=trace_status_weak_field
        return
      endif
      igrid2=ctx2%igrid

      x3=xnow+sgn*h*(b31*k1+b32*k2)
      call trace_spherical_interp_ctx_build_cached(x3,igrid2,sph_cache, &
           ctx3,ctx_status)
      if (ctx_status/=trace_status_active) then
        reject_reason=trace_rk45_reject_stage_failure
        if (ctx_status==trace_status_boundary .or. &
             ctx_status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif
      call get_K_spherical_ctx(ctx3,k3,field_type,threshold,field_ok)
      if (.not.field_ok) then
        status=trace_status_weak_field
        return
      endif
      igrid3=ctx3%igrid

      x4=xnow+sgn*h*(b41*k1+b42*k2+b43*k3)
      call trace_spherical_interp_ctx_build_cached(x4,igrid3,sph_cache, &
           ctx4,ctx_status)
      if (ctx_status/=trace_status_active) then
        reject_reason=trace_rk45_reject_stage_failure
        if (ctx_status==trace_status_boundary .or. &
             ctx_status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif
      call get_K_spherical_ctx(ctx4,k4,field_type,threshold,field_ok)
      if (.not.field_ok) then
        status=trace_status_weak_field
        return
      endif
      igrid4=ctx4%igrid

      x5s=xnow+sgn*h*(b51*k1+b52*k2+b53*k3+b54*k4)
      call trace_spherical_interp_ctx_build_cached(x5s,igrid4,sph_cache, &
           ctx5,ctx_status)
      if (ctx_status/=trace_status_active) then
        reject_reason=trace_rk45_reject_stage_failure
        if (ctx_status==trace_status_boundary .or. &
             ctx_status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif
      call get_K_spherical_ctx(ctx5,k5,field_type,threshold,field_ok)
      if (.not.field_ok) then
        status=trace_status_weak_field
        return
      endif
      igrid5=ctx5%igrid

      x6=xnow+sgn*h*(b61*k1+b62*k2+b63*k3+b64*k4+b65*k5)
      call trace_spherical_interp_ctx_build_cached(x6,igrid5,sph_cache, &
           ctx6,ctx_status)
      if (ctx_status/=trace_status_active) then
        reject_reason=trace_rk45_reject_stage_failure
        if (ctx_status==trace_status_boundary .or. &
             ctx_status==trace_status_out_of_domain) &
             reject_reason=trace_rk45_reject_stage_outside
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             reject_reason)
        if (h<=hfloor*(one+epsilon(one))) exit
        h=max(hfloor,half*h)
        cycle
      endif
      call get_K_spherical_ctx(ctx6,k6,field_type,threshold,field_ok)
      if (.not.field_ok) then
        status=trace_status_weak_field
        return
      endif

      x5=xnow+sgn*h*(c1*k1+c3*k3+c4*k4+c6*k6)
      x4err=xnow+sgn*h*(cs1*k1+cs3*k3+cs4*k4+cs5*k5+cs6*k6)

      point_domain=0
      {if (x5(^DB)>=domain_min(^DB) .and. x5(^DB)<domain_max(^DB)) point_domain=point_domain+1\}
      if (point_domain/=ndim) then
        call trace_intersect_domain(xnow,x5,xhit,hit_ok,face_id, &
             alpha_hit)
        if (hit_ok) then
          xnext=xhit
          ds_actual=max(zero,min(one,alpha_hit))*h
          h_suggest=hfloor
          status=trace_status_active
          call trace_rk45_stats_note_attempt(.true.,.true.,ds_actual)
          return
        endif
        call trace_rk45_stats_note_attempt(.false.,.false.,h, &
             trace_rk45_reject_boundary)
        status=trace_status_out_of_domain
        return
      endif

      err=zero
      {^IFTHREED
      r_metric=max(abs(xnow(1)),smalldouble)
      sin_theta_metric=max(abs(dsin(xnow(2))),smalldouble)
      err=dsqrt((x5(1)-x4err(1))**2+ &
           (r_metric*(x5(2)-x4err(2)))**2+ &
           (r_metric*sin_theta_metric*(x5(3)-x4err(3)))**2)
      }
      if (hcell<=zero) hcell=hmax
      tol=trace_rk45_atol+trace_rk45_rtol*max(h,hcell)
      if (err<=tol .or. h<=hfloor*(one+epsilon(one))) then
        xnext=x5
        ds_actual=h
        if (accumulate_twist) then
          call trace_twist_density_at_point(xnow,ctx1%igrid,threshold,tw1, &
               twist_status,sph_cache)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x3,ctx3%igrid,threshold,tw3, &
               twist_status,sph_cache)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x4,ctx4%igrid,threshold,tw4, &
               twist_status,sph_cache)
          if (twist_status==trace_status_active) &
               call trace_twist_density_at_point(x6,ctx6%igrid,threshold,tw6, &
               twist_status,sph_cache)
          if (twist_status==trace_status_active) then
            twist_increment=h*(c1*tw1+c3*tw3+c4*tw4+c6*tw6)
            twist_integrated=.true.
          endif
        endif
        if (err>zero .and. tol>zero) then
          grow=trace_rk45_safety*(tol/err)**0.2d0
          grow=max(trace_rk45_min_shrink,min(trace_rk45_max_grow,grow))
        else
          grow=trace_rk45_max_grow
        endif
        h_suggest=min(hmax,max(hfloor,h*grow))
        call trace_rk45_stats_note_attempt(.true.,.false.,h)
        status=trace_status_active
        return
      endif

      grow=trace_rk45_safety*(tol/max(err,smalldouble))**0.25d0
      grow=max(trace_rk45_min_shrink,min(one,grow))
      call trace_rk45_stats_note_attempt(.false.,.false.,h, &
           trace_rk45_reject_error)
      h=max(hfloor,h*grow)
    enddo

    status=trace_status_out_of_domain
  end subroutine trace_summary_rk45_spherical_step

  subroutine trace_summary_rk2_step(xnow,igrid,ds,forward,threshold, &
       xnext,ds_actual,status,sph_cache)
    double precision, intent(in) :: xnow(ndim),ds,threshold
    integer, intent(in) :: igrid
    logical, intent(in) :: forward
    double precision, intent(out) :: xnext(ndim),ds_actual
    integer, intent(out) :: status
    type(trace_sph_interp_ctx), intent(inout), optional :: sph_cache

    double precision :: xstage(ndim),xstage_sample(ndim),K1(ndim),K2(ndim)
    double precision :: ds_step
    double precision :: dxb^D
    type(trace_sph_interp_ctx) :: ctx_now,ctx_stage
    integer :: ixI^L
    integer :: igrid_stage
    integer :: ctx_status
    logical :: field_ok
    character(len=std_len) :: field_type

    ixI^L=ixG^LL;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    field_type='Bfield'
    xnext=xnow
    ds_actual=zero
    status=trace_status_active

    if (geo_coordinate==geo_spherical) then
      if (present(sph_cache)) then
        call trace_spherical_interp_ctx_build_cached(xnow,igrid,sph_cache, &
             ctx_now,ctx_status)
      else
        call trace_spherical_interp_ctx_build(xnow,igrid,ctx_now,ctx_status)
      endif
      if (ctx_status/=trace_status_active) then
        status=ctx_status
        return
      endif
      ds_step=trace_spherical_effective_step_ctx(ds,ctx_now)
      call get_K_spherical_ctx(ctx_now,K1,field_type,threshold,field_ok)
    else
      ds_step=trace_effective_step(xnow,ds,igrid,dxb^D)
      call get_K(xnow,igrid,K1,ixI^L,dxb^D,field_type,threshold,field_ok)
    endif
    if (.not.field_ok) then
      status=trace_status_weak_field
      return
    endif

    if (forward) then
      xstage=xnow+ds_step*K1
    else
      xstage=xnow-ds_step*K1
    endif
    xstage_sample=xstage
    ^D&xstage_sample(^D)=max(xprobmin^D,min(xstage_sample(^D), &
         xprobmax^D-100.d0*epsilon(one)*max(one,abs(xprobmax^D))));
    call trace_locate_point_with_hint(xstage_sample,igrid,igrid_stage,status)
    if (status/=trace_status_active) return
    ^D&dxb^D=rnode(rpdx^D_,igrid_stage);
    if (geo_coordinate==geo_spherical) then
      if (present(sph_cache)) then
        call trace_spherical_interp_ctx_build_cached(xstage_sample, &
             igrid_stage,sph_cache,ctx_stage,ctx_status)
      else
        call trace_spherical_interp_ctx_build(xstage_sample,igrid_stage, &
             ctx_stage,ctx_status)
      endif
      if (ctx_status/=trace_status_active) then
        status=ctx_status
        return
      endif
      call get_K_spherical_ctx(ctx_stage,K2,field_type,threshold,field_ok)
    else
      call get_K(xstage_sample,igrid_stage,K2,ixI^L,dxb^D,field_type, &
           threshold,field_ok)
    endif
    if (.not.field_ok) then
      status=trace_status_weak_field
      return
    endif

    if (forward) then
      xnext=xnow+ds_step*(half*K1+half*K2)
    else
      xnext=xnow-ds_step*(half*K1+half*K2)
    endif
    if (geo_coordinate==geo_spherical) then
      ds_actual=abs(ds_step)
    else
      ds_actual=dsqrt(sum((xnext-xnow)**2))
    endif
  end subroutine trace_summary_rk2_step

  double precision function trace_effective_step(x,ds,igrid,dxb^D) &
       result(ds_eff)
    double precision, intent(in) :: x(ndim),ds
    integer, intent(in) :: igrid
    double precision, intent(in) :: dxb^D

    double precision :: r,sin_theta,cell_scale,dxloc(ndim),step_cap,hcell
    type(trace_sph_interp_ctx) :: ctx
    integer :: status

    ds_eff=ds
    if (trace_cartesian_like_geometry()) then
      if (trace_step_control_mode==trace_step_control_cell_fraction) then
        step_cap=abs(ds)
        call trace_cartesian_local_cell_size(x,igrid,hcell,status)
        if (status==trace_status_active .and. hcell>zero) then
          ds_eff=dsign(min(step_cap,trace_step_fraction*hcell),ds)
          if (trace_step_min>zero) then
            ds_eff=dsign(max(abs(ds_eff),min(trace_step_min,step_cap)),ds)
          endif
          if (trace_rk2_stats_enabled) call trace_rk2_stats_note_step_limit(step_cap, &
               trace_step_fraction*hcell,abs(ds_eff))
        endif
      endif
      return
    endif
    if (geo_coordinate/=geo_spherical) return

    {^IFTHREED
    call trace_spherical_interp_ctx_build(x,igrid,ctx,status)
    if (status==trace_status_active) then
      ds_eff=trace_spherical_effective_step_ctx(ds,ctx)
      return
    endif

    step_cap=abs(ds)
    call trace_spherical_local_cell_widths(x,igrid,dxloc,status)
    if (status/=trace_status_active) return
    r=max(abs(x(1)),smalldouble)
    sin_theta=abs(dsin(x(2)))
    cell_scale=min(abs(dxloc(1)),r*abs(dxloc(2)))
    cell_scale=min(cell_scale,r*max(sin_theta,smalldouble)*abs(dxloc(3)))
    if (cell_scale>zero) then
      select case (trace_step_control_mode)
      case (trace_step_control_cell_fraction)
        ds_eff=min(step_cap,trace_step_fraction*cell_scale)
        if (trace_step_min>zero) ds_eff=max(ds_eff,min(trace_step_min,step_cap))
      case default
        ds_eff=min(step_cap,cell_scale)
      end select
    endif
    }
  end function trace_effective_step

  logical function trace_spherical_metric_ok(x) result(is_ok)
    double precision, intent(in) :: x(ndim)

    double precision :: metric_tol

    is_ok=.false.
    if (ndim/=3) return
    {^IFTHREED
    metric_tol=100.d0*epsilon(one)
    is_ok=x(1)>metric_tol .and. abs(dsin(x(2)))>metric_tol
    }
  end function trace_spherical_metric_ok

  double precision function trace_segment_length(xa,xb,ds_full,alpha) &
       result(length)
    double precision, intent(in) :: xa(ndim),xb(ndim),ds_full
    double precision, intent(in), optional :: alpha

    if (geo_coordinate==geo_spherical) then
      if (present(alpha)) then
        length=max(zero,min(one,alpha))*abs(ds_full)
      else
        length=abs(ds_full)
      endif
    else
      length=dsqrt(sum((xb-xa)**2))
    endif
  end function trace_segment_length

  subroutine trace_field_length_multi(seeds,nseed,dL,max_steps,results,b_min)
    ! Trace multiple seeds without storing field-line paths.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_length_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min

    double precision :: field_min

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    call trace_summary_multi(seeds,nseed,dL,max_steps,results,field_min)
  end subroutine trace_field_length_multi

  subroutine trace_field_twist_multi(seeds,nseed,dL,max_steps,results,b_min)
    ! Trace multiple seeds and integrate twist without storing field-line paths.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_twist_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min

    double precision :: field_min

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    call trace_summary_twist_multi(seeds,nseed,dL,max_steps,results,field_min)
  end subroutine trace_field_twist_multi

  subroutine trace_field_mapping_multi(seeds,nseed,dL,max_steps,results, &
       b_min,source_normal)
    ! Trace multiple seeds and return endpoint metadata without storing paths.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_mapping_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min
    double precision, intent(in), optional :: source_normal(3)

    double precision :: field_min,normal(3)
    logical :: have_normal

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    normal=zero
    have_normal=present(source_normal)
    if (have_normal) normal=source_normal
    call trace_summary_mapping_multi(seeds,nseed,dL,max_steps,results, &
         field_min,have_normal,normal)
  end subroutine trace_field_mapping_multi

  subroutine trace_field_topology_multi(seeds,nseed,dL,max_steps,results, &
       need_twist,need_mapping,b_min,source_normal)
    ! Trace multiple seeds once and return length plus optional Tw/mapping.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_topology_result), intent(out) :: results(nseed)
    logical, intent(in), optional :: need_twist,need_mapping
    double precision, intent(in), optional :: b_min
    double precision, intent(in), optional :: source_normal(3)

    double precision :: field_min,normal(3)
    logical :: do_twist,do_mapping,have_normal

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    do_twist=.false.
    do_mapping=.false.
    if (present(need_twist)) do_twist=need_twist
    if (present(need_mapping)) do_mapping=need_mapping
    normal=zero
    have_normal=present(source_normal)
    if (have_normal) normal=source_normal

    call trace_summary_topology_multi(seeds,nseed,dL,max_steps,results, &
         field_min,do_twist,do_mapping,have_normal,normal)
  end subroutine trace_field_topology_multi

  subroutine trace_field_qperp_multi(seeds,nseed,dL,max_steps,results,b_min, &
       twist_results)
    ! Block-grouped tangent-state tracing for multiple Method-II/Q_perp seeds.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_qperp_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min
    type(trace_twist_result), intent(out), optional :: twist_results(nseed)

    type(trace_tangent_state), allocatable :: states(:)
    double precision :: field_min
    double precision :: seed_local(ndim),zero_vec(ndim)
    integer :: common_status,iseed,idim,igrid,iforward,ibackward
    integer :: seed_status
    logical :: do_twist

    if (nseed<=0) return

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    zero_vec=zero
    do_twist=present(twist_results)

    common_status=trace_status_active
    if (npe/=1) then
      common_status=trace_status_mpi_unsupported
    else if (ndim/=3 .or. .not.trace_cartesian_like_geometry()) then
      common_status=trace_status_unsupported_geometry
    else if (dL<=zero .or. max_steps<=0) then
      common_status=trace_status_invalid_input
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        do idim=1,ndim
          seed_local(idim)=seeds(iseed,idim)
        enddo
        call trace_init_qperp_result(seed_local,results(iseed))
        results(iseed)%status=common_status
        if (do_twist) call trace_summary_init_twist_result(seed_local, &
             twist_results(iseed),common_status)
      enddo
      return
    endif

    allocate(states(2*nseed))
    do iseed=1,nseed
      do idim=1,ndim
        seed_local(idim)=seeds(iseed,idim)
      enddo
      call trace_qperp_prepare_seed_result(seed_local,field_min, &
           results(iseed),igrid,seed_status)
      iforward=2*iseed-1
      ibackward=2*iseed
      if (seed_status==trace_status_active) then
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,1,states(iforward), &
             accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,-1,states(ibackward), &
             accumulate_twist=do_twist)
        if (trace_integrator_mode==trace_integrator_rk45_cartesian .or. &
             trace_integrator_mode==trace_integrator_rk45_spherical) then
          call trace_rk45_stats_note_direction()
          call trace_rk45_stats_note_direction()
        else if (trace_integrator_mode==trace_integrator_rk2) then
          call trace_rk2_stats_note_direction()
          call trace_rk2_stats_note_direction()
        endif
      else
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,1,states(iforward),accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,-1,states(ibackward),accumulate_twist=do_twist)
        states(iforward)%active=.false.
        states(ibackward)%active=.false.
        states(iforward)%status=seed_status
        states(ibackward)%status=seed_status
      endif
    enddo

    select case (trace_integrator_mode)
    case (trace_integrator_rk45_cartesian)
      call trace_tangent_trace_states_grouped(states,2*nseed,dL,max_steps, &
           field_min,trace_mode=trace_tangent_group_rk45_cartesian)
    case (trace_integrator_rk45_spherical)
      call trace_tangent_trace_states_grouped_rk45_spherical(states, &
           2*nseed,dL,max_steps,field_min)
    case default
      call trace_tangent_trace_states_grouped(states,2*nseed,dL,max_steps, &
           field_min)
    end select

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      if (results(iseed)%status==trace_status_active) then
        call trace_qperp_finalize_from_states(states(iforward), &
             states(ibackward),results(iseed),field_min)
      endif
      if (do_twist) call trace_tangent_fill_twist_result(seeds(iseed,:), &
           states(iforward),states(ibackward),twist_results(iseed))
    enddo
    deallocate(states)
  end subroutine trace_field_qperp_multi

  subroutine trace_field_q0_multi_rk45_cartesian(seeds,nseed,dL, &
       max_steps,results,b_min,twist_results)
    ! Cartesian standard-q0 tangent transport using RK45 stages.
    ! The final q0 projection is valid for Cartesian-like grids because it
    ! uses domain-face endpoints, endpoint face-limit B, and transported
    ! tangent vectors from the same RK45 tangent states as Qperp.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_qperp_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min
    type(trace_twist_result), intent(out), optional :: twist_results(nseed)

    type(trace_tangent_state), allocatable :: states(:)
    double precision :: field_min
    double precision :: seed_local(ndim),zero_vec(ndim)
    integer :: common_status,iseed,idim,igrid,iforward,ibackward
    integer :: seed_status
    logical :: do_twist

    if (nseed<=0) return

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    zero_vec=zero
    do_twist=present(twist_results)

    common_status=trace_status_active
    if (npe/=1) then
      common_status=trace_status_mpi_unsupported
    else if (ndim/=3 .or. .not.trace_cartesian_like_geometry()) then
      common_status=trace_status_unsupported_geometry
    else if (dL<=zero .or. max_steps<=0) then
      common_status=trace_status_invalid_input
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        do idim=1,ndim
          seed_local(idim)=seeds(iseed,idim)
        enddo
        call trace_init_qperp_result(seed_local,results(iseed))
        results(iseed)%status=common_status
        results(iseed)%status_q0=common_status
        if (do_twist) call trace_summary_init_twist_result(seed_local, &
             twist_results(iseed),common_status)
      enddo
      return
    endif

    allocate(states(2*nseed))
    do iseed=1,nseed
      do idim=1,ndim
        seed_local(idim)=seeds(iseed,idim)
      enddo
      call trace_qperp_prepare_seed_result(seed_local,field_min, &
           results(iseed),igrid,seed_status)
      iforward=2*iseed-1
      ibackward=2*iseed
      if (seed_status==trace_status_active) then
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,1,states(iforward), &
             accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,-1,states(ibackward), &
             accumulate_twist=do_twist)
        call trace_rk45_stats_note_direction()
        call trace_rk45_stats_note_direction()
      else
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,1,states(iforward),accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,-1,states(ibackward),accumulate_twist=do_twist)
        states(iforward)%active=.false.
        states(ibackward)%active=.false.
        states(iforward)%status=seed_status
        states(ibackward)%status=seed_status
        results(iseed)%status_q0=seed_status
      endif
    enddo

    call trace_tangent_trace_states_grouped(states,2*nseed,dL,max_steps, &
         field_min,trace_mode=trace_tangent_group_rk45_cartesian)

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      if (results(iseed)%status==trace_status_active) then
        call trace_q0_finalize_from_states(states(iforward), &
             states(ibackward),results(iseed),field_min)
      endif
      if (do_twist) call trace_tangent_fill_twist_result(seeds(iseed,:), &
           states(iforward),states(ibackward),twist_results(iseed))
    enddo
    deallocate(states)
  end subroutine trace_field_q0_multi_rk45_cartesian

  subroutine trace_debug_cartesian_rk45_tangent_q0_multi(seeds,nseed,dL, &
       max_steps,results,b_min)
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_qperp_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min

    if (present(b_min)) then
      call trace_field_q0_multi_rk45_cartesian(seeds,nseed,dL,max_steps, &
           results,b_min)
    else
      call trace_field_q0_multi_rk45_cartesian(seeds,nseed,dL,max_steps, &
           results)
    endif
  end subroutine trace_debug_cartesian_rk45_tangent_q0_multi

  subroutine trace_field_spherical_rmin_q_multi(seeds,nseed,dL,max_steps, &
       results,b_min,twist_results)
    ! Legacy-named radial-boundary spherical Q from tangent transport.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_qperp_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min
    type(trace_twist_result), intent(out), optional :: twist_results(nseed)

    type(trace_tangent_state), allocatable :: states(:)
    double precision :: field_min
    double precision :: seed_local(ndim),zero_vec(ndim)
    integer :: common_status,iseed,idim,igrid,iforward,ibackward
    integer :: seed_status
    logical :: do_twist

    if (nseed<=0) return
    if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
         call trace_spherical_profile_add_count(trace_profile_q_tangent_systems, &
         int(2*nseed,kind=8))

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    zero_vec=zero
    do_twist=present(twist_results)

    common_status=trace_status_active
    if (npe/=1) then
      common_status=trace_status_mpi_unsupported
    else if (ndim/=3 .or. geo_coordinate/=geo_spherical) then
      common_status=trace_status_unsupported_geometry
    else if (dL<=zero .or. max_steps<=0) then
      common_status=trace_status_invalid_input
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        do idim=1,ndim
          seed_local(idim)=seeds(iseed,idim)
        enddo
        call trace_init_qperp_result(seed_local,results(iseed))
        results(iseed)%status=common_status
        results(iseed)%status_q0=common_status
        if (do_twist) call trace_summary_init_twist_result(seed_local, &
             twist_results(iseed),common_status)
      enddo
      return
    endif

    allocate(states(2*nseed))
    do iseed=1,nseed
      do idim=1,ndim
        seed_local(idim)=seeds(iseed,idim)
      enddo
      call trace_spherical_rmin_q_prepare_seed_result(seed_local,field_min, &
           results(iseed),igrid,seed_status)
      iforward=2*iseed-1
      ibackward=2*iseed
      if (seed_status==trace_status_active) then
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,1,states(iforward), &
             accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,-1,states(ibackward), &
             accumulate_twist=do_twist)
        if (trace_integrator_mode==trace_integrator_rk45_spherical) then
          call trace_rk45_stats_note_direction()
          call trace_rk45_stats_note_direction()
        else if (trace_integrator_mode==trace_integrator_rk2) then
          call trace_rk2_stats_note_direction()
          call trace_rk2_stats_note_direction()
        endif
      else
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,1,states(iforward),accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,-1,states(ibackward),accumulate_twist=do_twist)
        states(iforward)%active=.false.
        states(ibackward)%active=.false.
        states(iforward)%status=seed_status
        states(ibackward)%status=seed_status
      endif
    enddo

    if (trace_integrator_mode==trace_integrator_rk45_spherical) then
      call trace_tangent_trace_states_grouped_rk45_spherical(states, &
           2*nseed,dL,max_steps,field_min)
    else
      call trace_tangent_trace_states_grouped(states,2*nseed,dL,max_steps, &
           field_min)
    endif

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      if (results(iseed)%status==trace_status_active) then
        call trace_spherical_rmin_q_finalize_from_states(states(iforward), &
             states(ibackward),results(iseed),field_min)
      endif
      if (do_twist) call trace_tangent_fill_twist_result(seeds(iseed,:), &
           states(iforward),states(ibackward),twist_results(iseed))
    enddo
    deallocate(states)
  end subroutine trace_field_spherical_rmin_q_multi

  subroutine trace_field_rk2_short_boundary_q_multi(seeds,nseed,dL, &
       max_steps,results,b_min,twist_results)
    ! RK2 q0/logQ tangent tracing with a shortened final boundary step.
    ! Seed-products use this as the unified non-Qperp length/Tw/logQ trace.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_qperp_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min
    type(trace_twist_result), intent(out), optional :: twist_results(nseed)

    type(trace_tangent_state), allocatable :: states(:)
    double precision :: field_min
    double precision :: seed_local(ndim),zero_vec(ndim)
    integer :: common_status,iseed,idim,igrid,iforward,ibackward
    integer :: seed_status
    logical :: do_twist,is_spherical

    if (nseed<=0) return

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    zero_vec=zero
    do_twist=present(twist_results)
    is_spherical=(geo_coordinate==geo_spherical)

    common_status=trace_status_active
    if (npe/=1) then
      common_status=trace_status_mpi_unsupported
    else if (ndim/=3) then
      common_status=trace_status_unsupported_geometry
    else if (geo_coordinate==geo_cartesian .and. .not.slab_uniform) then
      common_status=trace_status_unsupported_geometry
    else if (geo_coordinate/=geo_cartesian .and. &
         geo_coordinate/=geo_spherical) then
      common_status=trace_status_unsupported_geometry
    else if (dL<=zero .or. max_steps<=0) then
      common_status=trace_status_invalid_input
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        do idim=1,ndim
          seed_local(idim)=seeds(iseed,idim)
        enddo
        call trace_init_qperp_result(seed_local,results(iseed))
        results(iseed)%status=common_status
        results(iseed)%status_q0=common_status
        if (do_twist) call trace_summary_init_twist_result(seed_local, &
             twist_results(iseed),common_status)
      enddo
      return
    endif

    allocate(states(2*nseed))
    do iseed=1,nseed
      do idim=1,ndim
        seed_local(idim)=seeds(iseed,idim)
      enddo
      if (is_spherical) then
        call trace_spherical_rmin_q_prepare_seed_result(seed_local,field_min, &
             results(iseed),igrid,seed_status)
      else
        call trace_qperp_prepare_seed_result(seed_local,field_min, &
             results(iseed),igrid,seed_status)
      endif
      iforward=2*iseed-1
      ibackward=2*iseed
      if (seed_status==trace_status_active) then
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,1,states(iforward), &
             accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,-1,states(ibackward), &
             accumulate_twist=do_twist)
        if (trace_integrator_mode==trace_integrator_rk45_spherical) then
          call trace_rk45_stats_note_direction()
          call trace_rk45_stats_note_direction()
        else if (trace_integrator_mode==trace_integrator_rk2) then
          call trace_rk2_stats_note_direction()
          call trace_rk2_stats_note_direction()
        endif
      else
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,1,states(iforward),accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,-1,states(ibackward),accumulate_twist=do_twist)
        states(iforward)%active=.false.
        states(ibackward)%active=.false.
        states(iforward)%status=seed_status
        states(ibackward)%status=seed_status
        results(iseed)%status_q0=seed_status
      endif
    enddo

    call trace_tangent_trace_states_grouped(states,2*nseed,dL,max_steps, &
         field_min,trace_mode=trace_tangent_group_rk2_short_boundary)

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      if (results(iseed)%status==trace_status_active) then
        if (is_spherical) then
          call trace_spherical_rmin_q_finalize_from_states(states(iforward), &
               states(ibackward),results(iseed),field_min)
        else
          call trace_q0_finalize_from_states(states(iforward), &
               states(ibackward),results(iseed),field_min)
        endif
      endif
      if (do_twist) call trace_tangent_fill_twist_result(seeds(iseed,:), &
           states(iforward),states(ibackward),twist_results(iseed))
    enddo
    deallocate(states)
  end subroutine trace_field_rk2_short_boundary_q_multi

  subroutine trace_field_spherical_qperp_multi(seeds,nseed,dL,max_steps, &
       results,b_min,twist_results)
    ! Spherical Q_perp from physical Cartesian tangent transport.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_qperp_result), intent(out) :: results(nseed)
    double precision, intent(in), optional :: b_min
    type(trace_twist_result), intent(out), optional :: twist_results(nseed)

    type(trace_tangent_state), allocatable :: states(:)
    double precision :: field_min
    double precision :: seed_local(ndim),zero_vec(ndim)
    integer :: common_status,iseed,idim,igrid,iforward,ibackward
    integer :: seed_status
    integer :: cache_status
    logical :: do_twist

    if (nseed<=0) return
    if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) &
         call trace_spherical_profile_add_count( &
         trace_profile_qperp_tangent_systems,int(2*nseed,kind=8))

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    zero_vec=zero
    do_twist=present(twist_results)

    common_status=trace_status_active
    if (npe/=1) then
      common_status=trace_status_mpi_unsupported
    else if (ndim/=3 .or. geo_coordinate/=geo_spherical) then
      common_status=trace_status_unsupported_geometry
    else if (dL<=zero .or. max_steps<=0) then
      common_status=trace_status_invalid_input
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        do idim=1,ndim
          seed_local(idim)=seeds(iseed,idim)
        enddo
        call trace_init_qperp_result(seed_local,results(iseed))
        results(iseed)%status=common_status
        results(iseed)%status_qperp0=common_status
        if (do_twist) call trace_summary_init_twist_result(seed_local, &
             twist_results(iseed),common_status)
      enddo
      return
    endif

    if (do_twist) then
      call trace_spherical_curl_cache_build(cache_status)
      if (cache_status/=trace_status_active) then
        do iseed=1,nseed
          do idim=1,ndim
            seed_local(idim)=seeds(iseed,idim)
          enddo
          call trace_init_qperp_result(seed_local,results(iseed))
          results(iseed)%status=cache_status
          results(iseed)%status_qperp0=cache_status
          call trace_summary_init_twist_result(seed_local, &
               twist_results(iseed),cache_status)
        enddo
        return
      endif
    endif

    allocate(states(2*nseed))
    do iseed=1,nseed
      do idim=1,ndim
        seed_local(idim)=seeds(iseed,idim)
      enddo
      call trace_spherical_qperp_prepare_seed_result(seed_local,field_min, &
           results(iseed),igrid,seed_status)
      iforward=2*iseed-1
      ibackward=2*iseed
      if (seed_status==trace_status_active) then
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,1,states(iforward), &
             accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,results(iseed)%u0, &
             results(iseed)%v0,igrid,iseed,-1,states(ibackward), &
             accumulate_twist=do_twist)
        if (trace_integrator_mode==trace_integrator_rk45_spherical) then
          call trace_rk45_stats_note_direction()
          call trace_rk45_stats_note_direction()
        else if (trace_integrator_mode==trace_integrator_rk2) then
          call trace_rk2_stats_note_direction()
          call trace_rk2_stats_note_direction()
        endif
      else
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,1,states(iforward),accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid, &
             iseed,-1,states(ibackward),accumulate_twist=do_twist)
        states(iforward)%active=.false.
        states(ibackward)%active=.false.
        states(iforward)%status=seed_status
        states(ibackward)%status=seed_status
      endif
    enddo

    if (trace_integrator_mode==trace_integrator_rk45_spherical) then
      call trace_tangent_trace_states_grouped_rk45_spherical(states, &
           2*nseed,dL,max_steps,field_min)
    else
      call trace_tangent_trace_states_grouped(states,2*nseed,dL,max_steps, &
           field_min)
    endif

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      if (results(iseed)%status==trace_status_active) then
        call trace_spherical_qperp_finalize_from_states(states(iforward), &
             states(ibackward),results(iseed),field_min)
      endif
      if (do_twist) call trace_tangent_fill_twist_result(seeds(iseed,:), &
           states(iforward),states(ibackward),twist_results(iseed))
    enddo
    deallocate(states)
    if (do_twist) call trace_spherical_curl_cache_clear()
  end subroutine trace_field_spherical_qperp_multi

  subroutine trace_field_spherical_rmin_q_qperp_multi(seeds,nseed,dL, &
       max_steps,q_results,qperp_results,b_min,twist_results)
    ! Shared spherical radial-boundary Q and Q_perp tangent transport.
    integer, intent(in) :: nseed,max_steps
    double precision, intent(in) :: seeds(nseed,ndim),dL
    type(trace_qperp_result), intent(out) :: q_results(nseed)
    type(trace_qperp_result), intent(out) :: qperp_results(nseed)
    double precision, intent(in), optional :: b_min
    type(trace_twist_result), intent(out), optional :: twist_results(nseed)

    type(trace_tangent_state), allocatable :: states(:)
    type(trace_tangent_state) :: qperp_forward_state,qperp_backward_state
    double precision :: field_min
    double precision :: seed_local(ndim),zero_vec(ndim)
    integer :: common_status,iseed,idim,igrid_q,igrid_qperp
    integer :: iforward,ibackward,seed_status_q,seed_status_qperp
    logical :: do_twist

    if (nseed<=0) return
    if (trace_spherical_profile_enabled .and. geo_coordinate==geo_spherical) then
      call trace_spherical_profile_add_count(trace_profile_q_tangent_systems, &
           int(2*nseed,kind=8))
      call trace_spherical_profile_add_count( &
           trace_profile_qperp_tangent_systems,int(2*nseed,kind=8))
      call trace_spherical_profile_add_count( &
           trace_profile_combined_tangent_systems,int(2*nseed,kind=8))
    endif

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    zero_vec=zero
    do_twist=present(twist_results)

    common_status=trace_status_active
    if (npe/=1) then
      common_status=trace_status_mpi_unsupported
    else if (ndim/=3 .or. geo_coordinate/=geo_spherical) then
      common_status=trace_status_unsupported_geometry
    else if (dL<=zero .or. max_steps<=0) then
      common_status=trace_status_invalid_input
    endif
    if (common_status/=trace_status_active) then
      do iseed=1,nseed
        do idim=1,ndim
          seed_local(idim)=seeds(iseed,idim)
        enddo
        call trace_init_qperp_result(seed_local,q_results(iseed))
        q_results(iseed)%status=common_status
        q_results(iseed)%status_q0=common_status
        call trace_init_qperp_result(seed_local,qperp_results(iseed))
        qperp_results(iseed)%status=common_status
        qperp_results(iseed)%status_qperp0=common_status
        if (do_twist) call trace_summary_init_twist_result(seed_local, &
             twist_results(iseed),common_status)
      enddo
      return
    endif

    allocate(states(2*nseed))
    do iseed=1,nseed
      do idim=1,ndim
        seed_local(idim)=seeds(iseed,idim)
      enddo
      call trace_spherical_rmin_q_prepare_seed_result(seed_local,field_min, &
           q_results(iseed),igrid_q,seed_status_q)
      call trace_spherical_qperp_prepare_seed_result(seed_local,field_min, &
           qperp_results(iseed),igrid_qperp,seed_status_qperp)

      iforward=2*iseed-1
      ibackward=2*iseed
      if (seed_status_q==trace_status_active .and. &
           seed_status_qperp==trace_status_active) then
        call trace_tangent_state_init(seed_local,q_results(iseed)%u0, &
             q_results(iseed)%v0,igrid_q,iseed,1,states(iforward), &
             accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,q_results(iseed)%u0, &
             q_results(iseed)%v0,igrid_q,iseed,-1,states(ibackward), &
             accumulate_twist=do_twist)
        states(iforward)%has_extra=.true.
        states(ibackward)%has_extra=.true.
        states(iforward)%p=qperp_results(iseed)%u0
        states(iforward)%q=qperp_results(iseed)%v0
        states(ibackward)%p=qperp_results(iseed)%u0
        states(ibackward)%q=qperp_results(iseed)%v0
        if (trace_integrator_mode==trace_integrator_rk45_spherical) then
          call trace_rk45_stats_note_direction()
          call trace_rk45_stats_note_direction()
        else if (trace_integrator_mode==trace_integrator_rk2) then
          call trace_rk2_stats_note_direction()
          call trace_rk2_stats_note_direction()
        endif
      else
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid_q, &
             iseed,1,states(iforward),accumulate_twist=do_twist)
        call trace_tangent_state_init(seed_local,zero_vec,zero_vec,igrid_q, &
             iseed,-1,states(ibackward),accumulate_twist=do_twist)
        states(iforward)%active=.false.
        states(ibackward)%active=.false.
        if (seed_status_q/=trace_status_active) then
          states(iforward)%status=seed_status_q
          states(ibackward)%status=seed_status_q
        else
          states(iforward)%status=seed_status_qperp
          states(ibackward)%status=seed_status_qperp
        endif
      endif
    enddo

    if (trace_integrator_mode==trace_integrator_rk45_spherical) then
      call trace_tangent_trace_states_grouped_rk45_spherical(states, &
           2*nseed,dL,max_steps,field_min)
    else
      call trace_tangent_trace_states_grouped(states,2*nseed,dL,max_steps, &
           field_min)
    endif

    do iseed=1,nseed
      iforward=2*iseed-1
      ibackward=2*iseed
      if (q_results(iseed)%status==trace_status_active) then
        call trace_spherical_rmin_q_finalize_from_states(states(iforward), &
             states(ibackward),q_results(iseed),field_min)
      endif
      if (qperp_results(iseed)%status==trace_status_active) then
        qperp_forward_state=states(iforward)
        qperp_backward_state=states(ibackward)
        qperp_forward_state%u_perp=states(iforward)%p_perp
        qperp_forward_state%v_perp=states(iforward)%q_perp
        qperp_backward_state%u_perp=states(ibackward)%p_perp
        qperp_backward_state%v_perp=states(ibackward)%q_perp
        call trace_spherical_qperp_finalize_from_states(qperp_forward_state, &
             qperp_backward_state,qperp_results(iseed),field_min)
      endif
      if (do_twist) call trace_tangent_fill_twist_result(seeds(iseed,:), &
           states(iforward),states(ibackward),twist_results(iseed))
    enddo
    deallocate(states)
  end subroutine trace_field_spherical_rmin_q_qperp_multi

  subroutine trace_field_multi(xfm,wPm,wLm,dL,numL,numP,nwP,nwL,forwardm,ftype,tcondi)
    ! trace multiple field lines
    ! xfm: locations of points at the field lines. User should provide xfm(1:numL,1,1:ndim)
    !   as seed points, then field line wills be traced from the seed points. xfm(1:numL,1,:) 
    !   are user given and other points are given by the subroutine as feedback.
    ! numL: number of field lines user wants to trace. user given
    ! numP: maximum number of points at the field line. User defined. Note that not every
    !   point of numP is valid, as the tracing will stop when the field line leave the
    !   simulation box. The number of valid point is given in wLm(1)
    ! wPm: point variables, which have values at each point of xfm. The way to get wP is
    !   user defined, with help of IO given by subroutine set_field_w in mod_usr_methods.
    !   User can calculate the density/temperature at the point and then store the values
    !   into wPm, or do something else.
    ! nwP: number of point variables. user given
    ! wLm: line variables, variables for the lines rather than for each point of the lines.
    !   For example, wLm(1) stores the number of valid points at the field lines. The way to 
    !   get wLm is also user defined with set_field_w, the same as wPm. User can calculate
    !   the maximum density/temperature at the field lines and stores it in wLm
    ! nwL: number of line variables. user given
    ! dL: length step for field line tracing. user given
    ! forwardm: true--trace field forward; false--trace field line backward. user given
    ! ftype: type of field user wants to trace. User can trace velocity field by setting
    !   ftype='Vfield' or trace magnetic field by setting ftype='Bfield'. It is possible 
    !   to trace other fields, e.g. electric field, where user can define the field with
    !   IO given by subroutine set_field in mod_usr_methods. user given
    ! tcondi: user given
    use mod_particle_base

    integer, intent(in) :: numL,numP,nwP,nwL
    double precision, intent(inout) :: xfm(numL,numP,ndim),wPm(numL,numP,nwP),wLm(numL,1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forwardm(numL)
    character(len=std_len), intent(in) :: ftype,tcondi

    double precision :: x3d(3),statusF(4+ndim),statusL(numL,4+ndim+nwL),statusS(numL,4+ndim+nwL)
    double precision :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, allocatable :: data_send(:,:,:),data_recv(:,:,:)
    integer :: indomain,ipe_now,igrid_now,igrid,j,iL
    integer :: ipoint_in,ipoint_out,numSend,nRT,nRTmax
    integer :: ipointm(numL),igridm(numL)
    logical :: continueL(numL),myL(numL)
    logical :: stopT,forward

    if (tcondi/='TRAC') then
      wPm=zero
    else
      wPm=-1
    endif
    wLm=zero
    xfm(1:numL,2:numP,:)=zero
    stopT=.TRUE.
    myL=.FALSE.
    xf=zero
    wP=zero

    ! find the pe and igrid for the first point
    do iL=1,numL
      indomain=0
      wLm(iL,1)=0
      {if (xfm(iL,1,^DB)>=xprobmin^DB .and. xfm(iL,1,^DB)<xprobmax^DB) indomain=indomain+1\}
      if (indomain==ndim) then
        if (tcondi/='TRAC') wLm(iL,1)=1
        continueL(iL)=.TRUE.
        ! find pe and igrid
        x3d=0.d0
        do j=1,ndim
          x3d(j)=xfm(iL,1,j)
        enddo
        call find_particle_ipe(x3d,igrid_now,ipe_now)
        stopT=.FALSE.
        ipointm(iL)=1
        igridm(iL)=igrid_now
        if (mype==ipe_now) then 
          myL(iL)=.TRUE.
        else
          xfm(iL,1,:)=zero
        endif
      else
        continueL(iL)=.FALSE.
        wLm(iL,1)=zero
      endif
    enddo

    do while(stopT .eqv. .FALSE.)
      ! tracing multiple field lines inside pe
      statusS=zero
      do iL=1,numL
        if (myL(iL) .and. continueL(iL)) then
          igrid=igridm(iL)
          ipoint_in=ipointm(iL)
          xf(ipoint_in,:)=xfm(iL,ipoint_in,:)
          wL(:)=wLm(iL,:)
          forward=forwardm(iL)
          statusF=zero
          call find_points_in_pe(igrid,ipoint_in,xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi,statusF)
          ipoint_out=int(statusF(1))
          xfm(iL,ipoint_in:ipoint_out-1,:)=xf(ipoint_in:ipoint_out-1,:)
          wPm(iL,ipoint_in:ipoint_out-1,:)=wP(ipoint_in:ipoint_out-1,:)
          ! status for each field line
          ! 1: index of next point
          ! 2: ipe of next point
          ! 3: igrid of next point
          ! 4: trace_status_active -> continue; otherwise stop tracing
          ! 5:4+ndim: coordinate of next point
          ! 4+ndim+1:4+ndim+nwL: wL(2:1+nwL)
          ! for TRAC nwL=2 -> wL(2): current Tcoff; wL(3): Tmax
          ! for TRAC nwP=2 -> wP(:,1):ipe; wP(:,2):igrid
          statusS(iL,1:4+ndim)=statusF(1:4+ndim)
          statusS(iL,4+ndim+1:4+ndim+nwL)=wL(2:1+nwL)
          if (tcondi=='TRAC') wLm(iL,1)=ipoint_out-1
        endif
      enddo

      ! comunicating tracing results
      numSend=numL*(4+ndim+nwL)
      call MPI_ALLREDUCE(statusS,statusL,numSend,MPI_DOUBLE_PRECISION,&
                       MPI_SUM,icomm,ierrmpi)

      ! for next step
      stopT=.TRUE.
      myL=.FALSE.
      do iL=1,numL
        if (continueL(iL)) then
          ipointm(iL)=int(statusL(iL,1))
          if (mype==int(statusL(iL,2))) myL(iL)=.TRUE.
          igridm(iL)=int(statusL(iL,3))
          if (int(statusL(iL,4))==trace_status_active) then
            stopT=.FALSE.
          else
            continueL(iL)=.FALSE.
          endif
          if (myL(iL)) xfm(iL,ipointm(iL),1:ndim)=statusL(iL,4+1:4+ndim)
          if (tcondi/='TRAC') then
            if (int(statusL(iL,4))==trace_status_weak_field) then
              wLm(iL,1)=ipointm(iL)
            else
              wLm(iL,1)=ipointm(iL)-1
            endif
          endif
          wLm(iL,2:1+nwL)=statusL(iL,4+ndim+1:4+ndim+nwL)
        endif
      enddo
    enddo

    ! communication after tracing
    if (tcondi/='TRAC') then
      nRTmax=0
      do iL=1,numL
        if (nRTmax<int(wLm(iL,1))) nRTmax=int(wLm(iL,1))
      enddo
      numSend=numL*nRTmax*(ndim+nwP)

      allocate(data_send(numL,nRTmax,ndim+nwP),data_recv(numL,nRTmax,ndim+nwP))
      data_send(:,:,:)=zero
      do iL=1,numL
        nRT=int(wLm(iL,1))
        data_send(iL,1:nRT,1:ndim)=xfm(iL,1:nRT,1:ndim)
        if (nwP>0) data_send(iL,1:nRT,1+ndim:ndim+nwP)=wPm(iL,1:nRT,1:nwP)
      enddo
      call MPI_ALLREDUCE(data_send,data_recv,numSend,MPI_DOUBLE_PRECISION,&
                         MPI_SUM,icomm,ierrmpi)
      do iL=1,numL
        nRT=int(wLm(iL,1))
        xfm(iL,1:nRT,1:ndim)=data_recv(iL,1:nRT,1:ndim)
        if (nwP>0) wPm(iL,1:nRT,1:nwP)=data_recv(iL,1:nRT,1+ndim:ndim+nwP)
      enddo
      deallocate(data_send,data_recv)
    endif

  end subroutine trace_field_multi

  subroutine trace_field_single(xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi)
    ! trace a field line
    ! xf: locations of points at the field line. User should provide xf(1,1:ndim)
    !   as seed point, then field line will be traced from the seed point. xf(1,:) is 
    !   user given and xf(2:wL(1),:) are given by the subroutine as feedback.
    ! numP: maximum number of points at the field line. User defined. Note that not every
    !   point of numP is valid, as the tracing will stop when the field line leave the
    !   simulation box. The number of valid point is given in wL(1)
    ! wP: point variables, which have values at each point of xf. The way to get wP is
    !   user defined, with help of IO given by subroutine set_field_w in mod_usr_methods.
    !   User can calculate the density/temperature at the point and then store the values
    !   into wP, or do something else.
    ! nwP: number of point variables. user given
    ! wL: line variables, variables for the line rather than for each point of the line.
    !   For example, wL(1) stores the number of valid points at the field line. The way to 
    !   get wL is also user defined with set_field_w, the same as wP. User can calculate
    !   the maximum density/temperature at the field line and stores it in wL
    ! nwL: number of line variables. user given
    ! dL: length step for field line tracing. user given
    ! forward: true--trace field forward; false--trace field line backward. user given
    ! ftype: type of field user wants to trace. User can trace velocity field by setting
    !   ftype='Vfield' or trace magnetic field by setting ftype='Bfield'. It is possible 
    !   to trace other fields, e.g. electric field, where user can define the field with
    !   IO given by subroutine set_field in mod_usr_methods. user given
    ! tcondi: user given
    use mod_usr_methods
    use mod_particle_base

    integer, intent(in) :: numP,nwP,nwL
    double precision, intent(inout) :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forward
    character(len=std_len), intent(in) :: ftype,tcondi

    double precision :: x3d(3),statusF(4+ndim),status_bcast(4+ndim+nwL)
    double precision, allocatable :: data_send(:,:),data_recv(:,:)
    integer :: indomain,ipoint_in,ipe_now,igrid_now,igrid,j
    integer :: ipoint_out,ipe_next,igrid_next,numRT
    logical :: stopT

    wP=zero
    wL=zero
    xf(2:numP,:)=zero

    ! check whether or the first point is inside simulation box. if yes, find
    ! the pe and igrid for the point
    indomain=0
    wL(1)=0
    {if (xf(1,^DB)>=xprobmin^DB .and. xf(1,^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain==ndim) then
      wL(1)=1

       ! find pe and igrid
       x3d=0.d0
       do j=1,ndim
         x3d(j)=xf(1,j)
       enddo
      call find_particle_ipe(x3d,igrid_now,ipe_now)
      stopT=.FALSE.
      ipoint_in=1
      if (mype/=ipe_now) xf(1,:)=zero
    else
      if (mype==0) then
        call MPISTOP('Field tracing error: given point is not in simulation box!')
      endif
    endif


    ! other points in field line
    do while(stopT .eqv. .FALSE.)

      if (mype==ipe_now) then
        igrid=igrid_now
        ! looking for points in one pe
        call find_points_in_pe(igrid,ipoint_in,xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi,statusF)
        status_bcast(1:4+ndim)=statusF(1:4+ndim)
        status_bcast(4+ndim+1:4+ndim+nwL)=wL(2:1+nwL)
      endif
      ! comunication
      call MPI_BCAST(status_bcast,4+ndim+nwL,MPI_DOUBLE_PRECISION,ipe_now,icomm,ierrmpi)
      statusF(1:4+ndim)=status_bcast(1:4+ndim)
      wL(2:1+nwL)=status_bcast(4+ndim+1:4+ndim+nwL)

      ! prepare for next step
      ipoint_out=int(statusF(1))
      ipe_next=int(statusF(2))
      igrid_next=int(statusF(3))
      if (int(statusF(4))/=trace_status_active) then
        stopT=.TRUE.
        if (int(statusF(4))==trace_status_weak_field) then
          wL(1)=ipoint_out
        else
          wL(1)=ipoint_out-1
        endif
      endif
      if (mype==ipe_next) then
        do j=1,ndim
          xf(ipoint_out,j)=statusF(4+j)
        enddo
      else
        xf(ipoint_out,:)=zero
      endif

      ! pe and grid of next point
      ipe_now=ipe_next
      igrid_now=igrid_next
      ipoint_in=ipoint_out
    enddo

    if (tcondi/='TRAC') then
      numRT=int(wL(1))
      allocate(data_send(numRT,ndim+nwP),data_recv(numRT,ndim+nwP))
      data_send(:,:)=zero
      data_recv(:,:)=zero
      data_send(1:numRT,1:ndim)=xf(1:numRT,1:ndim)
      if (nwP>0) data_send(1:numRT,1+ndim:ndim+nwP)=wP(1:numRT,1:nwP)
      call MPI_ALLREDUCE(data_send,data_recv,numRT*(ndim+nwP),MPI_DOUBLE_PRECISION,&
                         MPI_SUM,icomm,ierrmpi)
      xf(1:numRT,1:ndim)=data_recv(1:numRT,1:ndim)
      if (nwP>0) wP(1:numRT,1:nwP)=data_recv(1:numRT,1+ndim:ndim+nwP)
      deallocate(data_send,data_recv)
    endif

  end subroutine trace_field_single

  subroutine find_points_in_pe(igrid,ipoint_in,xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi,statusF)

    integer, intent(inout) :: igrid
    integer, intent(in) :: ipoint_in,numP,nwP,nwL
    double precision, intent(inout) :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forward
    character(len=std_len), intent(in) :: ftype,tcondi
    double precision, intent(inout) :: statusF(4+ndim)

    double precision :: xfout(ndim)
    integer :: ipe_next,igrid_next,ip_in,ip_out,j,indomain,trace_status
    logical :: newpe,stopT

    ip_in=ipoint_in
    newpe=.FALSE.
    ipe_next=mype
    igrid_next=igrid
    statusF=zero

    do while(newpe .eqv. .FALSE.)
      ! looking for points in given grid    
      call find_points_interp(igrid,ip_in,ip_out,xf,wP,wL,numP,nwP,nwL, &
           dL,forward,ftype,tcondi,trace_status)
      ip_in=ip_out

      ! when next point is out of given grid, find next grid  
      if (trace_status/=trace_status_active) then
        newpe=.TRUE.
        stopT=.TRUE.
      else
        indomain=0
        {if (xf(ip_out,^DB)>=xprobmin^DB .and. xf(ip_out,^DB)<=xprobmax^DB) indomain=indomain+1\}
        if (ip_out<numP .and. indomain==ndim) then
          if (tcondi/='TRAC') then
            stopT=.FALSE.
            xfout=xf(ip_out,:)
            call find_next_grid(igrid,igrid_next,ipe_next,xfout,newpe,stopT)
            if (stopT) trace_status=trace_status_boundary
          else
            if (xf(ip_out,ndim)>phys_trac_mask) then
              newpe=.TRUE.
              stopT=.TRUE.
              trace_status=trace_status_trac_stop
            else
              stopT=.FALSE.
              xfout=xf(ip_out,:)
              call find_next_grid(igrid,igrid_next,ipe_next,xfout,newpe,stopT)
              if (stopT) trace_status=trace_status_boundary
            endif
          endif
        else
          newpe=.TRUE.
          stopT=.TRUE.
          if (ip_out>=numP) then
            trace_status=trace_status_max_steps
          else
            trace_status=trace_status_boundary
          endif
        endif
      endif

      if (newpe) then
        statusF(1)=ip_out
        statusF(2)=ipe_next
        statusF(3)=igrid_next
        statusF(4)=trace_status_active
        if (stopT) statusF(4)=trace_status
        do j=1,ndim
          statusF(4+j)=xf(ip_out,j)
        enddo
      endif

      if (newpe .eqv. .FALSE.) igrid=igrid_next
    enddo

  end subroutine find_points_in_pe

  subroutine find_next_grid(igrid,igrid_next,ipe_next,xf1,newpe,stopT)
    ! check the grid and pe of next point
    use mod_usr_methods
    use mod_global_parameters
    use mod_forest

    integer, intent(inout) :: igrid,igrid_next,ipe_next
    double precision, intent(in) :: xf1(ndim)
    logical, intent(inout) :: newpe,stopT

    double precision :: dxb^D,xb^L,xbmid^D
    double precision :: xbn^L
    integer :: idn^D,my_neighbor_type,inblock
    integer :: ic^D,inc^D,ipe_neighbor,igrid_neighbor

    igrid_next=igrid
    ipe_next=mype

    ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
    inblock=0

    ! direction of next grid
    idn^D=0\
    {if (xf1(^D)<=xbmin^D) idn^D=-1\}
    {if (xf1(^D)>=xbmax^D) idn^D=1\}
    my_neighbor_type=neighbor_type(idn^D,igrid)
    igrid_neighbor=neighbor(1,idn^D,igrid)
    ipe_neighbor=neighbor(2,idn^D,igrid)

    ! ipe and igrid of next grid
    select case(my_neighbor_type)
    case (neighbor_boundary)
      ! next point is not in simulation box
      newpe=.TRUE.
      stopT=.TRUE.

    case(neighbor_coarse)
      ! neighbor grid has lower refinement level      
      igrid_next=igrid_neighbor
      ipe_next=ipe_neighbor
      if (mype==ipe_neighbor) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif

    case(neighbor_sibling)
      ! neighbor grid has lower refinement level 
      igrid_next=igrid_neighbor
      ipe_next=ipe_neighbor
      if (mype==ipe_neighbor) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif

    case(neighbor_fine)
      ! neighbor grid has higher refinement level 
      {xbmid^D=(xbmin^D+xbmax^D)/2.d0\}
      ^D&inc^D=1\
      {if (xf1(^D)<=xbmin^D) inc^D=0\}
      {if (xf1(^D)>xbmin^D .and. xf1(^D)<=xbmid^D) inc^D=1\}
      {if (xf1(^D)>xbmid^D .and. xf1(^D)<xbmax^D) inc^D=2\}
      {if (xf1(^D)>=xbmax^D) inc^D=3\}
      ipe_next=neighbor_child(2,inc^D,igrid)
      igrid_next=neighbor_child(1,inc^D,igrid)
      if (mype==ipe_next) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif
    end select

  end subroutine find_next_grid

  subroutine find_points_interp(igrid,ip_in,ip_out,xf,wP,wL,numP,nwP,nwL, &
       dL,forward,ftype,tcondi,trace_status)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in) :: igrid,ip_in,numP,nwP,nwL
    integer, intent(inout) :: ip_out
    double precision, intent(inout) :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forward
    character(len=std_len), intent(in) :: ftype,tcondi
    integer, intent(out) :: trace_status

    double precision :: dxb^D,xb^L
    double precision :: field(ixg^T,ndir)
    double precision :: xs1(ndim),xs2(ndim),K1(ndim),K2(ndim)
    double precision :: xfpre(ndim),xfnow(ndim),xfnext(ndim)
    double precision :: Tpre,Tnow,Tnext,dTds,Lt,Lr,ds,T_bott,trac_delta
    integer          :: ip,inblock,ixI^L,ixO^L,j
    logical :: field_ok

    ixI^L=ixG^LL;
    ixO^L=ixM^LL;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    ^D&xbmin^D=rnode(rpxmin^D_,igrid);
    ^D&xbmax^D=rnode(rpxmax^D_,igrid);
    ip_out=ip_in
    trace_status=trace_status_active

    if (tcondi/='TRAC') then
      ds=dL
    else
      ds=dxb^ND
      T_bott=2.d4/unit_temperature
      trac_delta=0.25d0
    endif

    ! main loop
    MAINLOOP: do ip=ip_in,numP-1

      ! integrate magnetic field with Runge-Kutta method
      xs1(:)=xf(ip,:)
      call get_K(xs1,igrid,K1,ixI^L,dxb^D,ftype,smalldouble,field_ok)
      if (.not.field_ok) then
        ip_out=ip
        trace_status=trace_status_weak_field
        return
      endif
      if (forward) then
        xs2(:)=xf(ip,:)+ds*K1(:)
      else
        xs2(:)=xf(ip,:)-ds*K1(:)
      endif
      call get_K(xs2,igrid,K2,ixI^L,dxb^D,ftype,smalldouble,field_ok)
      if (.not.field_ok) then
        ip_out=ip
        trace_status=trace_status_weak_field
        return
      endif
      if (forward) then
        xf(ip+1,:)=xf(ip,:)+ds*(0.5*K1(:)+0.5*K2(:))
      else
        xf(ip+1,:)=xf(ip,:)-ds*(0.5*K1(:)+0.5*K2(:))
      endif
      ip_out=ip+1

      ! get local values for variable via interpolation
      if (tcondi/='TRAC') then
        if (associated(usr_set_field_w)) then 
          call usr_set_field_w(igrid,ip,xf,wP,wL,numP,nwP,nwL,dL,forward,ftype,tcondi)
        endif
      else  ! get TRAC Tcoff
        wP(ip,1)=mype
        wP(ip,2)=igrid
        if (ip==ip_in) then
          if (forward) then
            xfpre(:)=xf(ip,:)-ds*K1(:)
          else
            xfpre(:)=xf(ip,:)+ds*K1(:)
          endif
          xfnow(:)=xf(ip,:)
          xfnext(:)=xf(ip+1,:)
          call get_T_loc_TRAC(igrid,xfpre,Tpre,ixI^L,dxb^D)
          call get_T_loc_TRAC(igrid,xfnow,Tnow,ixI^L,dxb^D)
          call get_T_loc_TRAC(igrid,xfnext,Tnext,ixI^L,dxb^D)
        else
          xfpre=xf(ip-1,:)
          xfnow(:)=xf(ip,:)
          xfnext(:)=xf(ip+1,:)
          Tpre=Tnow
          Tnow=Tnext
          call get_T_loc_TRAC(igrid,xfnext,Tnext,ixI^L,dxb^D)
        endif
        dTds=abs(Tnext-Tpre)/(2*ds)
        if (ip==1) then
          Lt=0.d0
          wL(2)=T_bott   ! current Tcofl
          wL(3)=Tnow     ! current Tlmax
        else
          Lt=0.d0
          if (dTds>0.d0) then
            Lt=Tnow/dTds
            Lr=ds
            ! renew cutoff temperature
            if(Lr>trac_delta*Lt) then
              if (Tnow>wL(2)) wL(2)=Tnow
            endif
          endif
          if (Tnow>wL(3)) wL(3)=Tnow
        endif
      endif

      ! exit loop if next point is not in this block
      inblock=0
      {if (xf(ip+1,^DB)>=xbmin^DB .and. xf(ip+1,^DB)<xbmax^DB) inblock=inblock+1\}
      if (tcondi=='TRAC' .and. xf(ip+1,ndim)>phys_trac_mask) inblock=0
      if (inblock/=ndim) exit MAINLOOP

    enddo MAINLOOP

  end subroutine find_points_interp

  subroutine get_K(xfn,igrid,K,ixI^L,dxb^D,ftype,b_min,field_ok)
    use mod_usr_methods

    integer :: ixI^L,igrid
    double precision :: dxb^D
    double precision :: xfn(ndim),K(ndim)
    character(len=std_len) :: ftype
    double precision, intent(in), optional :: b_min
    logical, intent(out), optional :: field_ok

    double precision :: dxc^D,xd^D
    double precision :: field(0:1^D&,ndir),Fx(ndim),factor(0:1^D&)
    double precision :: Ftotal,field_min
    logical :: valid_field
    integer          :: ixb^D,ix^D,ixbl^D,j,status_interp

    if (geo_coordinate==geo_spherical) then
      call trace_interp_weights_block(xfn,igrid,ixI^L,ixbl^D,xd^D,dxc^D, &
           status_interp)
      if (status_interp/=trace_status_active) then
        if (present(field_ok)) field_ok=.false.
        K=zero
        return
      endif
    else
      ^D&ixbl^D=floor((xfn(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
      ^D&xd^D=(xfn(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D;
    endif

    field=zero
    if(ftype=='Bfield') then
      if(B0field) then
        if(allocated(iw_mag)) then
          {do ix^D=0,1\}
            field(ix^D,1:ndir)=ps(igrid)%w(ixbl^D+ix^D,iw_mag(1:ndir)) &
                 +ps(igrid)%B0(ixbl^D+ix^D,1:ndir,0)
          {enddo\}
        else
          {do ix^D=0,1\}
            field(ix^D,1:ndir)=ps(igrid)%B0(ixbl^D+ix^D,1:ndir,0)
          {enddo\}
        endif
      else
        {do ix^D=0,1\}
          field(ix^D,1:ndir)=ps(igrid)%w(ixbl^D+ix^D,iw_mag(1:ndir))
        {enddo\}
      endif
    else if (ftype=='Vfield') then
      block
        double precision :: w_stencil(ixbl^D:ixbl^D+1^D&,1:nw)
        double precision :: x_stencil(ixbl^D:ixbl^D+1^D&,1:ndim)
        double precision :: vector_stencil(ixbl^D:ixbl^D+1^D&,1:ndir)
        integer :: ixS^L

        ^D&ixSmin^D=ixbl^D;
        ^D&ixSmax^D=ixbl^D+1;
        w_stencil(ixS^S,1:nw)=ps(igrid)%w(ixS^S,1:nw)
        x_stencil(ixS^S,1:ndim)=ps(igrid)%x(ixS^S,1:ndim)
        call phys_get_v(w_stencil,x_stencil,ixS^L,ixS^L,vector_stencil)
        {do ix^D=0,1\}
          field(ix^D,1:ndir)=vector_stencil(ixbl^D+ix^D,1:ndir)
        {enddo\}
      end block
    endif
    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}

    if (ftype=='Bfield' .or. ftype=='Vfield') then
      Fx=0.d0
      {do ix^DB=0,1\}
        do j=1,ndim
          Fx(j)=Fx(j)+field(ix^D,j)*factor(ix^D)
        enddo
      {enddo\}
    else if (associated(usr_set_field)) then
      call usr_set_field(xfn,igrid,Fx,ftype)
    else
      call MPISTOP('Field tracing error: wrong field type!')
    endif

    Ftotal=zero
    do j=1,ndim
      Ftotal=Ftotal+(Fx(j))**2
    enddo
    Ftotal=dsqrt(Ftotal)

    field_min=smalldouble
    if (present(b_min)) field_min=max(b_min,zero)
    valid_field=Ftotal>zero .and. Ftotal>=field_min
    if (present(field_ok)) field_ok=valid_field

    K=zero
    if (valid_field) then
      select case (geo_coordinate)
      case (geo_spherical)
        {^IFTHREED
        if (.not.trace_spherical_metric_ok(xfn)) then
          valid_field=.false.
          if (present(field_ok)) field_ok=.false.
          return
        endif
        K(1)=Fx(1)/Ftotal
        K(2)=Fx(2)/(xfn(1)*Ftotal)
        K(3)=Fx(3)/(xfn(1)*dsin(xfn(2))*Ftotal)
        }
      case default
        K(1:ndim)=Fx(1:ndim)/Ftotal
      end select
    endif

  end subroutine get_K

  subroutine get_T_loc_TRAC(igrid,xloc,Tloc,ixI^L,dxb^D)
    ! grid T has been calculated and stored in wextra(ixI^S,Tcoff_)
    ! see mod_trac
    integer, intent(in) :: igrid,ixI^L
    double precision, intent(inout) :: xloc(ndim)
    double precision, intent(inout) :: Tloc
    double precision, intent(in) :: dxb^D

    double precision :: xd^D
    double precision :: factor(0:1^D&),Tnear(0:1^D&)
    integer          :: ixb^D,ix^D,ixbl^D,j,ixO^L

    ^D&ixbl^D=floor((xloc(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
    ^D&xd^D=(xloc(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D;
    ^D&ixOmin^D=ixbl^D;
    ^D&ixOmax^D=ixOmin^D+1;

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
      Tnear(ix^D)=ps(igrid)%wextra(ixbl^D+ix^D,iw_tcoff)
    {enddo\}

    Tloc=0.d0
    ! interpolation
    {do ix^DB=0,1\}
      Tloc=Tloc+Tnear(ix^D)*factor(ix^D)
    {enddo\}

  end subroutine get_T_loc_TRAC

end module mod_trace_field
