!> Minimal spherical-coordinate TDm/RBSL flux-rope initial condition.
!> The analytic TDm field is evaluated in Cartesian coordinates and transformed
!> back to AMRVAC spherical components (r, theta, phi; theta is colatitude).
!>
!> In AMRVAC spherical geometry, magnetic variables are stored internally as
!> local physical components mag(:)=(B_r,B_theta,B_phi). Default VTU conversion
!> writes vector components in Cartesian coordinates while preserving names such
!> as b1, b2, b3.
module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  use mod_global_parameters, only: par_files,mype,refine_max_level,&
      restart_from_file,undefined,firstprocess
  use mod_magnetic_topology, only: mt_params_read, mt_run_topology_task
  use mod_usr_methods, only: usr_special_convert,usr_refine_grid
  implicit none

  integer :: np

  ! TDm / background-field parameters in normalized units after init.
  double precision :: q_para, d_para, L_para
  double precision :: a0, F_flx
  double precision :: radial_center0

  ! Physical parameters read from amrvac.par.
  double precision :: minor_radius_cm
  double precision :: rho0, pressure0
  integer :: n_axis_points

  double precision, parameter :: solar_radius_cm = 6.961d10

  ! Cartesian coordinates of the flux-rope axis used by the RBSL integral.
  double precision, allocatable :: x_axis(:,:)
  logical :: tdm_setup_ready = .false.

contains

  subroutine usr_init()

    ! Local-box length scale used to make raw Cartesian/spherical snapshots
    ! directly comparable: 1 code unit = 1e9 cm = 10 Mm.
    unit_length        = 1.d9 ! cm
    unit_numberdensity = 1.d9     ! cm-3,cm-3,cm-3
    unit_temperature   = 1.d6     ! K

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_init_vector_potential => initvecpot_usr
    usr_special_convert => spherical_topology_special_convert
    usr_refine_grid => spherical_tdm_refine_grid

    call set_coordinate_system('spherical_3D')
    call mhd_activate()

  end subroutine usr_init

  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /usr_list/ minor_radius_cm, n_axis_points, rho0, pressure0

    minor_radius_cm = 2.25d9
    n_axis_points = 400
    rho0 = 1.d0
    pressure0 = 1.d0

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status='old')
      read(unitpar, usr_list, end=111)
111   close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine spherical_topology_special_convert(qunitconvert)
    integer, intent(in) :: qunitconvert
    double precision :: t0,t1

    call mt_params_read(par_files)
    call spherical_topology_print_amr_summary()
    t0 = MPI_WTIME()
    call mt_run_topology_task()
    t1 = MPI_WTIME()
    if (mype==0) write(*,'(a,es14.6)') 'TOPOLOGY_BENCH_WALL seconds=',t1-t0
  end subroutine spherical_topology_special_convert

  subroutine spherical_topology_print_amr_summary()
    use mod_forest, only: nleafs,nleafs_active,nleafs_level
    integer :: level

    if (mype==0) then
      write(*,'(a,i0,a,i0,a,i0)') 'TOPOLOGY_BENCH_AMR nleafs=',nleafs,&
          ' active=',nleafs_active,' refine_max_level=',refine_max_level
      do level=1,refine_max_level
        write(*,'(a,i0,a,i0)') 'TOPOLOGY_BENCH_AMR_LEVEL level=',level,&
            ' nleafs=',nleafs_level(level)
      end do
    end if
  end subroutine spherical_topology_print_amr_summary

  subroutine initglobaldata_usr()
    double precision :: r0
    double precision :: b_perp_apex, shafranov_factor, mu0I_equilibrium
    double precision :: source_depth_cm, source_halfsep_cm
    double precision :: axis_half_width_cm, axis_top_height_cm,&
        axis_bottom_height_cm
    double precision :: path_minor_cm

    call usr_params_read(par_files)

    r0 = 8.d0
    a0 = minor_radius_cm/unit_length

    if (a0 <= zero .or. a0 >= r0) then
      call mpistop('TDm requires 0 < minor_radius_cm < 8e9 cm')
    end if

    ! Exact circular TDm axis in the local Cartesian frame.
    radial_center0 = solar_radius_cm/unit_length

    ! A normal restart reads the magnetic field from the snapshot.  Keep only
    ! the lightweight scalars needed by optional AMR/output callbacks and skip
    ! all analytic TDm/RBSL axis construction.  firstprocess=.true. is the
    ! explicit exception because AMRVAC then reruns usr_init_one_grid.
    tdm_setup_ready = .false.
    if(restart_from_file/=undefined .and. .not.firstprocess) then
      if(allocated(x_axis)) deallocate(x_axis)
      if(mype==0) print *,&
          'Restart mode: skipping spherical TDm/RBSL field setup'
      return
    end if

    d_para = 4.5d0
    L_para = 3.d0
    ! Cartesian bipoB sums nb0=3 coincident source pairs. This spherical
    ! setup uses one vector-potential source pair, whose curl has the opposite
    ! sign relative to the direct Cartesian source-field expression.
    q_para = -600.d0/sqrt(4.d0*dpi)
    ! At the torus apex, the bipolar background field is perpendicular to the
    ! torus plane. Titov et al. (2014), Equations (7) and (11), give the
    ! Shafranov equilibrium current and flux for the parabolic current profile.
    ! The same flux-current relation is Equation (14) of Titov et al. (2021).
    b_perp_apex = abs(2.d0*L_para*q_para/&
       (r0**2+L_para**2)**1.5d0)
    shafranov_factor = log(8.d0*r0/a0)-25.d0/24.d0
    if (shafranov_factor <= zero) then
      call mpistop('Invalid TDm Shafranov equilibrium factor')
    end if
    mu0I_equilibrium = 4.d0*dpi*r0*b_perp_apex/shafranov_factor
    F_flx = 3.d0*mu0I_equilibrium*a0/(5.d0*sqrt(2.d0))

    np = max(16, n_axis_points)
    if(allocated(x_axis)) deallocate(x_axis)
    allocate(x_axis(np, ndim))

    call calc_cartesian_tdm_axis(x_axis, np, radial_center0)
    tdm_setup_ready = .true.

    source_depth_cm = d_para*unit_length
    source_halfsep_cm = L_para*unit_length
    axis_half_width_cm = maxval(abs(x_axis(:,3)))*unit_length
    axis_top_height_cm = (maxval(x_axis(:,1)) - radial_center0)*unit_length
    axis_bottom_height_cm = (minval(x_axis(:,1)) - radial_center0)*unit_length
    path_minor_cm = a0*unit_length

    if(mype == 0) then
      print *, 'Minimal spherical TDm initial condition'
      print *, 'Using one complete circular TDm/RBSL integration path'
      print *, 'unit_length [cm]: ', unit_length
      print *, 'solar radius [code units]: ', radial_center0
      print *, 'q_para normalized: ', q_para
      print *, 'apex strapping-field magnitude: ', b_perp_apex
      print *, 'equilibrium mu0*I normalized: ', mu0I_equilibrium
      print *, 'F_flx normalized: ', F_flx
      print *, 'source depth [cm]: ', source_depth_cm
      print *, 'source half separation [cm]: ', source_halfsep_cm
      print *, 'major radius [cm]: ', r0*unit_length
      print *, 'axis transverse half-width [cm]: ', axis_half_width_cm
      print *, 'axis radial height range [cm]: ', axis_bottom_height_cm, axis_top_height_cm
      print *, 'minor radius [cm]: ', path_minor_cm
      print *, 'axis integration points: ', np
    end if

  end subroutine initglobaldata_usr

  subroutine calc_cartesian_tdm_axis(xs,npo,radial_center)
    integer, intent(in) :: npo
    double precision, intent(in) :: radial_center
    double precision, dimension(npo,3), intent(inout) :: xs

    integer :: i
    double precision :: theta_start, theta, xref, zref
    double precision, parameter :: arc_center_x = 0.d0
    double precision, parameter :: arc_center_z = -4.5d0
    double precision, parameter :: arc_radius = 8.d0

    ! The TDm field is analytic, so evaluate the RBSL line integral on the
    ! complete torus circle. The endpoint is not duplicated because the RBSL
    ! quadrature closes the path periodically between points npo and 1.
    theta_start = asin(-arc_center_z/arc_radius)
    do i = 1, npo
      theta = theta_start + 2.d0*dpi*dble(i-1)/dble(npo)
      xref = arc_center_x + arc_radius*cos(theta)
      zref = arc_center_z + arc_radius*sin(theta)
      xs(i,1) = radial_center + zref
      xs(i,2) = zero
      xs(i,3) = xref
    end do
  end subroutine calc_cartesian_tdm_axis

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: Bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    logical, save :: first = .true.

    if(.not.tdm_setup_ready) then
      call mpistop('TDm field initialization called after restart setup skip')
    end if

    if(first) then
      if(mype == 0) print *, 'Constructing TDm/RBSL magnetic field'
      first = .false.
    end if

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) = rho0
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(:)) = 0.d0
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_) = pressure0

    if(stagger_grid) then
      call b_from_vector_potential(block%ixGsmin1,block%ixGsmin2,&
         block%ixGsmin3,block%ixGsmax1,block%ixGsmax2,block%ixGsmax3, ixImin1,&
         ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
         ixOmax1,ixOmax2,ixOmax3, block%ws, x)
      call mhd_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
          block)
    else
      ! Fast cell-centered path used by the init-data pars. This mirrors the
      ! Cartesian demo: direct B-field accumulation in Cartesian coordinates,
      ! followed by conversion to local spherical physical components.
      call tdm_magnetic_field_spherical(ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, x,&
          Bvec)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:)) = Bvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)
    end if

    if(mhd_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_) = 0.d0

    call eos%to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)

  end subroutine initonegrid_usr

  subroutine spherical_tdm_refine_grid(igrid, level, ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      qt, w, x, refine, coarsen)
    ! Deterministic AMR benchmark criterion in the same local Cartesian frame
    ! used by the spherical TDm/RBSL analytic field evaluation.
    integer, intent(in) :: igrid, level, ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qt, w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    integer, intent(inout) :: refine, coarsen

    double precision :: xcart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: xloc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        yloc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        zloc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    logical :: broad_region, rope_region, core_region

    xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2))*cos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2))*sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*cos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))

    zloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1) - radial_center0
    yloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)
    xloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)

    broad_region = any(abs(xloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) <= 8.5d0 .and. abs(yloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)) <= 4.5d0 .and. zloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) >= 0.d0 .and. zloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) <= 7.5d0)
    rope_region = any(abs(xloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) <= 7.5d0 .and. abs(yloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)) <= 2.5d0 .and. zloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) >= 0.d0 .and. zloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) <= 5.5d0) .or. &
       any(abs(xloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) <= 9.d0 .and. abs(yloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)) <= 4.d0 .and. zloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) >= 0.d0 .and. zloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) <= 2.d0)
    core_region = any(abs(xloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) <= 5.5d0 .and. abs(yloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)) <= 1.8d0 .and. zloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) >= 0.5d0 .and. zloc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) <= 4.5d0)

    if (level==1 .and. broad_region) then
      refine = 1
      coarsen = -1
    else if (level==2 .and. rope_region) then
      refine = 1
      coarsen = -1
    else if (level==3 .and. core_region) then
      refine = 1
      coarsen = -1
    end if
  end subroutine spherical_tdm_refine_grid

  subroutine initvecpot_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, xC, A, idir)
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idir
    double precision, intent(in) :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: A(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: Avec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)

    if(.not.tdm_setup_ready) then
      call mpistop('TDm vector potential called after restart setup skip')
    end if

    call tdm_vector_potential_spherical(ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, xC,&
        Avec)
    A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = Avec(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)

  end subroutine initvecpot_usr

  subroutine tdm_magnetic_field_spherical(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, xsph,&
      Bsph)
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: xsph(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: Bsph(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    double precision :: xcart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: xlocal(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: B_bg(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim), B_rope(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision :: B_local(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim), B_cart(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision :: axis_local(np,1:ndim)
    double precision :: q_direct

    B_bg = 0.d0
    B_rope = 0.d0

    ! Local Cartesian coordinates matching the Cartesian TDm/RBSL reference:
    ! xlocal=(x_ref, y_ref, z_ref), with z_ref radial-outward from R_sun.
    xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) = xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*sin(xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2))*cos(xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2) = xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*sin(xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2))*sin(xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3) = xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*cos(xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
    xlocal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) = xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
    xlocal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2) = xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
    xlocal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3) = xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) - radial_center0

    axis_local(:,1) = x_axis(:,3)
    axis_local(:,2) = x_axis(:,2)
    axis_local(:,3) = x_axis(:,1) - radial_center0

    q_direct = -q_para/3.d0
    call bipolar_field_direct_B(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, xlocal,&
        L_para, d_para, q_direct, 0.d0, 3, B_bg)
    call RBSL_flux_rope_direct_B(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, np, a0, F_flx,&
        .false., xlocal, axis_local, B_rope)

    B_local = B_bg + B_rope

    ! Convert local reference components back to the standard Cartesian basis
    ! used by AMRVAC's spherical vector transform.
    B_cart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) = B_local(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
    B_cart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2) = B_local(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
    B_cart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3) = B_local(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
    call Cart2SphereVector(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, xsph, B_cart, Bsph)

  end subroutine tdm_magnetic_field_spherical

  subroutine tdm_vector_potential_spherical(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, xsph,&
      Asph)
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: xsph(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: Asph(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    double precision :: xcart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: A_bg(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim), A_rope(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision :: A_cart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)

    ! Standard AMRVAC spherical coordinates: x=(r, theta, phi), theta colatitude.
    xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) = xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*sin(xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2))*cos(xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2) = xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*sin(xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2))*sin(xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    xcart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3) = xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*cos(xsph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))

    A_bg = 0.d0
    A_rope = 0.d0
    call bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, xcart, A_bg)
    call RBSL_flux_rope(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, np, a0, F_flx, .true.,&
        xcart, x_axis, A_rope)

    A_cart = A_bg + A_rope
    call Cart2SphereVector(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, xsph, A_cart, Asph)

  end subroutine tdm_vector_potential_spherical


  subroutine Cart2SphereVector(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,A_in,A_out)

    integer :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision :: x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: A_in(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: A_out(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)

    real*8 :: lon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       lat(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: bxCart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       byCart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       bzCart(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: br(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       bth(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       bph(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: a11(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a12(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a13(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: a21(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a22(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a23(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: a31(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a32(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       a33(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real*8 :: latc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       lonc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       pAng(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    bxCart=0.0d0
    byCart=0.0d0
    bzCart=0.0d0
    lon=0.0d0
    lat=0.0d0
    bph=0.0d0
    bth=0.0d0
    br=0.0d0
    A_out=0.0d0

    bxCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = A_in(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)
    byCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = A_in(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)
    bzCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = A_in(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)
    lon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) =  x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)
    lat(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) =  0.5d0*dpi - &
       x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
    latc = 0.0d0
    lonc = 0.0d0
    pAng = 0.0d0

    a11 = -sin(latc) * sin(pAng) * sin(lon - lonc) + cos(pAng) * cos(lon - &
       lonc)
    a12 =  sin(latc) * cos(pAng) * sin(lon - lonc) + sin(pAng) * cos(lon - &
       lonc)
    a13 = -cos(latc) * sin(lon - lonc)
    a21 = -sin(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * &
       sin(lon - lonc)) - cos(lat) * cos(latc) * sin(pAng)
    a22 =  sin(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * &
       sin(lon - lonc)) + cos(lat) * cos(latc) * cos(pAng)
    a23 = -cos(latc) * sin(lat) * cos(lon - lonc) + sin(latc) * cos(lat)
    a31 =  cos(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * &
       sin(lon - lonc)) - sin(lat) * cos(latc) * sin(pAng)
    a32 = -cos(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * &
       sin(lon - lonc)) + sin(lat) * cos(latc) * cos(pAng)
    a33 =  cos(lat) * cos(latc) * cos(lon - lonc) + sin(lat) * sin(latc)

    bph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = a11(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * bxCart(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) +a12(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * byCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) +a13(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * bzCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    bth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = a21(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * bxCart(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) +a22(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * byCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) +a23(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * bzCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    br(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = a31(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * bxCart(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) +a32(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * byCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) +a33(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * bzCart(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    A_out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) = br(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    A_out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2) = -1.0d0*bth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    A_out(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3) = bph(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine Cart2SphereVector

  subroutine bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,A,Bbp)

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: A(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, optional, intent(out) :: Bbp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)

    double precision :: Aphi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    Aphi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)= q_para*(L_para-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2))/(sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)+d_para-radial_center0)**2)*sqrt(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)**2+(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)+d_para-radial_center0)**2+&
       (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)-L_para)**2))+q_para*(L_para+x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2))/(sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)+d_para-radial_center0)**2)*sqrt(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)**2+(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)+d_para-radial_center0)**2+&
       (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)+L_para)**2))

    A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)=-Aphi(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)/sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)+d_para-radial_center0)**2)
    A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)= 0.d0
    A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)= Aphi(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)+d_para-radial_center0)/sqrt(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)**2+(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)+d_para-radial_center0)**2)

    if(present(Bbp)) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)+d_para-radial_center0)**2+(x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)+L_para)**2)**3
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)=       (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)+L_para)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)=                x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)= (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)+d_para-radial_center0)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)+d_para-radial_center0)**2+(x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)-L_para)**2)**3
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)=      -(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)-L_para)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)=               -x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)=-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)+d_para-radial_center0)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)=q_para*Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)
    end if

  end subroutine bipolar_field

  subroutine bipolar_field_direct_B(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,L_cha,d_cha,&
     q_cha,x_cha,nb,Bout)
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,nb
    double precision, intent(in) :: L_cha,d_cha,q_cha,x_cha
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out) :: Bout(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    integer :: i
    double precision :: xpos
    double precision :: rpv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        rmv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: rplus(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim), rminus(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)

    Bout = 0.d0
    do i=1,nb
      xpos = 2.d0*dble(i-1)*x_cha/dble(max(1,nb-1))-x_cha
      rplus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)-xpos
      rminus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)-xpos
      rplus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)-L_cha
      rminus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)+L_cha
      rplus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)+d_cha
      rminus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)+d_cha
      rpv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sqrt(rplus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)**2+rplus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)**2+rplus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)**2)
      rmv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sqrt(rminus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)**2+rminus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)**2+rminus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)**2)
      Bout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1) = Bout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)+q_cha*(rplus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)/rpv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**3-rminus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)/rmv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**3)
      Bout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2) = Bout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)+q_cha*(rplus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)/rpv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**3-rminus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)/rmv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**3)
      Bout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3) = Bout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)+q_cha*(rplus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)/rpv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**3-rminus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)/rmv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**3)
    end do

  end subroutine bipolar_field_direct_B

  subroutine RBSL_flux_rope(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,np,a,F_flx,&
     positive_helicity,x,x_axis,Atotal,Bfr)

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,np
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: x_axis(np,1:ndim)
    double precision, intent(in)    :: a
    double precision, intent(in)    :: F_flx
    logical, intent(in)             :: positive_helicity
    double precision, intent(out)   :: Atotal(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, optional, intent(out) :: Bfr(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)

    double precision :: I_cur
    double precision :: AIx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: AFx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    integer :: ix1,ix2,ix3, ixp, idirmin, ixMmin1,ixMmin2,ixMmin3,ixMmax1,&
       ixMmax2,ixMmax3
    double precision :: r_mag, KIr, KFr, re_pi, sqrt1r, f52r, fsqrt6
    double precision :: Rpl(1:ndim), r_vec(1:ndim), Rcr(1:ndim)
    double precision :: axis_element(np,1:ndim)

    if(positive_helicity) then
      I_cur = 5.d0*sqrt(2.d0)*F_flx/(3.d0*4.d0*dpi*a)
    else
      I_cur =-5.d0*sqrt(2.d0)*F_flx/(3.d0*4.d0*dpi*a)
    end if

    re_pi=1.d0/dpi
    AIx = 0.d0
    AFx = 0.d0
    fsqrt6=1.d0/sqrt(6.d0)
    do ixp=1,np
      if (ixp == 1) then
        axis_element(ixp,:) = 0.5d0*(x_axis(ixp+1,:)-x_axis(np,:))
      else if (ixp == np) then
        axis_element(ixp,:) = 0.5d0*(x_axis(1,:)-x_axis(ixp-1,:))
      else
        axis_element(ixp,:) = 0.5d0*(x_axis(ixp+1,:)-x_axis(ixp-1,:))
      end if
    end do
    if(present(Bfr)) then
      ixMmin1=ixImin1;ixMmin2=ixImin2;ixMmin3=ixImin3;ixMmax1=ixImax1
      ixMmax2=ixImax2;ixMmax3=ixImax3;
    else
      ixMmin1=ixOmin1;ixMmin2=ixOmin2;ixMmin3=ixOmin3;ixMmax1=ixOmax1
      ixMmax2=ixOmax2;ixMmax3=ixOmax3;
    end if
    do ix3=ixMmin3,ixMmax3
    do ix2=ixMmin2,ixMmax2
    do ix1=ixMmin1,ixMmax1
      do ixp=1,np
        r_vec(:) = (x(ix1,ix2,ix3,:) - x_axis(ixp,:))/a
        r_mag = sqrt(sum(r_vec(:)**2))
        Rpl(:) = axis_element(ixp,:)
        Rcr(1) = Rpl(2)*r_vec(3) - Rpl(3)*r_vec(2)
        Rcr(2) = Rpl(3)*r_vec(1) - Rpl(1)*r_vec(3)
        Rcr(3) = Rpl(1)*r_vec(2) - Rpl(2)*r_vec(1)
        if (r_mag <= 1.d-3) then
          ! Analytic r -> 0 limits from Titov et al. (2021), Equations (6, 11).
          KIr = 16.d0/(3.d0*dpi)
          KFr = 10.d0/(3.d0*dpi) + 5.d0/(2.d0*sqrt(6.d0)*dpi)*&
             (dpi-2.d0*asin(0.2d0))
        else if (r_mag < 1.d0) then
          sqrt1r=sqrt(1.d0-r_mag**2)
          f52r=5.d0-2.d0*r_mag**2
          KIr = 2.d0*re_pi*(asin(r_mag)/r_mag + f52r*third*sqrt1r)
          KFr = 2.d0*re_pi/r_mag**2*(asin(r_mag)/r_mag-sqrt1r) + &
             2.d0*re_pi*sqrt1r + f52r*0.5d0*fsqrt6*(1.d0 - &
             2.d0*re_pi*asin((1.d0+2.d0*r_mag**2)/f52r))
        else
          KIr = 1.d0/r_mag
          KFr = KIr**3
        endif
        AIx(ix1,ix2,ix3,:) = AIx(ix1,ix2,ix3,:) + KIr*Rpl(:)
        AFx(ix1,ix2,ix3,:) = AFx(ix1,ix2,ix3,:) + KFr*Rcr(:)
      end do
      AIx(ix1,ix2,ix3,:) = AIx(ix1,ix2,ix3,:)*I_cur/a
      AFx(ix1,ix2,ix3,:) = AFx(ix1,ix2,ix3,:)*F_flx*0.25d0*re_pi/a**2
    end do
    end do
    end do
    Atotal=AIx+AFx
    if(present(Bfr)) call curlvector(Atotal,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Bfr,&
       idirmin,1,ndir)

  end subroutine RBSL_flux_rope

  subroutine RBSL_flux_rope_direct_B(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,np,a,F_flx,&
     positive_helicity,x,x_axis,Btotal)
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,np
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: x_axis(np,1:ndim)
    double precision, intent(in)    :: a
    double precision, intent(in)    :: F_flx
    logical, intent(in)             :: positive_helicity
    double precision, intent(out)   :: Btotal(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    double precision :: I_cur
    double precision :: BIx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: BFx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    integer :: ix1,ix2,ix3, ixp
    double precision :: r_mag, KIr, KFr1, KFr2, Rdr
    double precision :: asr, or2, asrr, re_pi
    double precision :: Rpl(1:ndim), r_vec(1:ndim), Rcr(1:ndim)
    double precision :: axis_element(np,1:ndim)

    if(positive_helicity) then
      I_cur = 5.d0*sqrt(2.d0)*F_flx/(3.d0*a)
    else
      I_cur =-5.d0*sqrt(2.d0)*F_flx/(3.d0*a)
    end if

    re_pi=1.d0/dpi
    BIx = 0.d0
    BFx = 0.d0

    do ixp=1,np
      if (ixp == 1) then
        axis_element(ixp,:) = 0.5d0*(x_axis(ixp+1,:)-x_axis(np,:))
      else if (ixp == np) then
        axis_element(ixp,:) = 0.5d0*(x_axis(1,:)-x_axis(ixp-1,:))
      else
        axis_element(ixp,:) = 0.5d0*(x_axis(ixp+1,:)-x_axis(ixp-1,:))
      end if
    end do

    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
      do ixp=1,np
        r_vec(:) = (x(ix1,ix2,ix3,:) - x_axis(ixp,:))/a
        r_mag = sqrt(sum(r_vec(:)**2))
        Rpl(:) = axis_element(ixp,:)
        Rcr(1) = Rpl(2)*r_vec(3) - Rpl(3)*r_vec(2)
        Rcr(2) = Rpl(3)*r_vec(1) - Rpl(1)*r_vec(3)
        Rcr(3) = Rpl(1)*r_vec(2) - Rpl(2)*r_vec(1)
        Rdr = sum(Rpl(:)*r_vec(:))
        if (r_mag <= 1.d-3) then
          ! Analytic r -> 0 limits from Titov et al. (2021), Equations (A2,
          ! A7, and A8), avoiding cancellation in the finite-r expressions.
          KIr = 16.d0/(3.d0*dpi)
          KFr1 = 5.d0/sqrt(6.d0) + 10.d0/dpi*&
             (2.d0/3.d0-asin(0.2d0)/sqrt(6.d0))
          KFr2 = sqrt(6.d0)/3.d0 + 2.d0/(15.d0*dpi)*&
             (24.d0-5.d0*sqrt(6.d0)*asin(0.2d0))
        else if (r_mag <= 1.d0) then
          asr=asin(r_mag)/r_mag
          or2=sqrt(1.d0-r_mag**2)
          asrr=asin((1.d0+2.d0*r_mag*r_mag)/(5.d0-2.d0*r_mag*r_mag))
          KIr=2.d0*re_pi*((asr-or2)/r_mag**2+2.d0*or2)
          KFr1=2.d0*re_pi/r_mag**2*(or2-asr)+8.d0*re_pi*or2+&
             (5.d0-4.d0*r_mag**2)/sqrt(6.d0)*(1.d0-2.d0*re_pi*asrr)
          KFr2=2.d0*re_pi/r_mag**4*(3.d0*asr-(3.d0+2.d0*r_mag**2)*or2)+&
             2.d0/sqrt(6.d0)*(1.d0-2.d0*re_pi*asrr)
        else
          KIr=1.d0/r_mag**3
          KFr1=-1.d0/r_mag**3
          KFr2=3.d0/r_mag**5
        end if
        BIx(ix1,ix2,ix3,:) = BIx(ix1,ix2,ix3,&
           :) + I_cur*0.25d0*re_pi*KIr*Rcr(:)/a**2
        BFx(ix1,ix2,ix3,:) = BFx(ix1,ix2,ix3,&
           :) + F_flx*0.25d0*re_pi*(KFr1*Rpl(:)+KFr2*Rdr*r_vec(:))/a**3
      end do
    end do
    end do
    end do
    Btotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       :) = BIx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       :) + BFx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)

  end subroutine RBSL_flux_rope_direct_B

end module mod_usr
