!> Spherical TDm/RBSL radiation-synthesis demonstration.
!> The magnetic field and AMR regions are identical to MagneticTopologyQSL_Sph.
!> The atmosphere follows an AL-C7 hydrostatic profile and prominence material
!> is placed in radial magnetic dips inside the analytic flux rope.
!>
!> In AMRVAC spherical geometry, magnetic variables are stored internally as
!> local physical components mag(:)=(B_r,B_theta,B_phi). Default VTU conversion
!> writes vector components in Cartesian coordinates while preserving names such
!> as b1, b2, b3.
module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  use mod_solar_atmosphere, only: get_atm_para
  use mod_global_parameters, only: par_files,mype,refine_max_level,&
      restart_from_file,undefined,firstprocess
  use mod_usr_methods, only: usr_set_parameters,usr_init_one_grid,&
      usr_init_vector_potential,usr_refine_grid,usr_aux_output,&
      usr_add_aux_names
  implicit none

  integer, parameter :: n_atmosphere=4096
  integer :: np

  ! TDm / background-field parameters in normalized units after init.
  double precision :: q_para, d_para, L_para
  double precision :: a0, F_flx
  double precision :: radial_center0

  ! Physical parameters read from amrvac.par.
  double precision :: minor_radius_cm,rho_ref_height_cm,rho_ref
  double precision :: prom_density_factor,dip_bz_max,prom_axis_fraction
  integer :: n_axis_points
  character(len=16) :: atmosphere_curve
  double precision :: atmosphere_hmin,atmosphere_dh

  double precision, parameter :: solar_radius_cm = 6.961d10

  ! Cartesian coordinates of the flux-rope axis used by the RBSL integral.
  double precision, allocatable :: x_axis(:,:)
  double precision, allocatable :: atmosphere_rho(:),atmosphere_p(:),&
      atmosphere_T(:)
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
    usr_refine_grid => spherical_tdm_refine_grid
    usr_aux_output => specialvar_output
    usr_add_aux_names => specialvarnames_output

    call set_coordinate_system('spherical_3D')
    call mhd_activate()

  end subroutine usr_init

  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /usr_list/ minor_radius_cm,n_axis_points,atmosphere_curve,&
        rho_ref_height_cm,rho_ref,prom_density_factor,dip_bz_max,&
        prom_axis_fraction

    minor_radius_cm = 2.25d9
    n_axis_points = 400
    atmosphere_curve = 'AL-C7'
    rho_ref_height_cm = 1.d9
    rho_ref = 0.5d0
    prom_density_factor = 100.d0
    dip_bz_max = 0.4d0
    prom_axis_fraction = 1.d0

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status='old')
      read(unitpar, usr_list, end=111)
111   close(unitpar)
    end do

  end subroutine usr_params_read

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
    if(rho_ref<=zero) call mpistop('rho_ref must be positive')
    if(prom_density_factor<one) then
      call mpistop('prom_density_factor must be at least one')
    end if
    if(dip_bz_max<zero) call mpistop('dip_bz_max must be non-negative')
    if(prom_axis_fraction<=zero) then
      call mpistop('prom_axis_fraction must be positive')
    end if

    ! Exact circular TDm axis in the local Cartesian frame.
    radial_center0 = solar_radius_cm/unit_length

    ! Conversion/diagnostic restarts consume the state stored in the snapshot.
    ! Preserve only scalars needed by optional output/refinement callbacks and
    ! skip the analytic axis and atmosphere table.  firstprocess=.true. remains
    ! the explicit route for deliberately regenerating the initial condition.
    tdm_setup_ready = .false.
    if(restart_from_file/=undefined .and. .not.firstprocess) then
      if(allocated(x_axis)) deallocate(x_axis)
      if(allocated(atmosphere_rho)) deallocate(atmosphere_rho)
      if(allocated(atmosphere_p)) deallocate(atmosphere_p)
      if(allocated(atmosphere_T)) deallocate(atmosphere_T)
      if(mype==0) print *,&
          'Restart mode: skipping spherical TDm/atmosphere setup'
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
    call initialize_solar_atmosphere()
    tdm_setup_ready = .true.

    source_depth_cm = d_para*unit_length
    source_halfsep_cm = L_para*unit_length
    axis_half_width_cm = maxval(abs(x_axis(:,3)))*unit_length
    axis_top_height_cm = (maxval(x_axis(:,1)) - radial_center0)*unit_length
    axis_bottom_height_cm = (minval(x_axis(:,1)) - radial_center0)*unit_length
    path_minor_cm = a0*unit_length

    if(mype == 0) then
      print *, 'Spherical TDm/RBSL radiation-synthesis initial condition'
      print *, 'Using the MagneticTopologyQSL_Sph analytic magnetic field'
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
      print *, 'atmosphere curve: ',trim(atmosphere_curve)
      print *, 'reference height [Mm]: ',rho_ref_height_cm/1.d8
      print *, 'reference number density [cm^-3]: ',&
          rho_ref*unit_numberdensity
      print *, 'prominence density factor: ',prom_density_factor
      print *, 'dip |Br| threshold: ',dip_bz_max
      print *, 'dip axis-distance limit [code units]: ',&
          prom_axis_fraction*a0
    end if

  end subroutine initglobaldata_usr

  subroutine initialize_solar_atmosphere()
    double precision :: h(n_atmosphere),grav(n_atmosphere)
    double precision :: atmosphere_hmax,gravity0,reference_height,padding
    integer :: j

    padding = 2.d8/unit_length
    atmosphere_hmin = xprobmin1-radial_center0-padding
    atmosphere_hmax = xprobmax1-radial_center0+padding
    atmosphere_dh = (atmosphere_hmax-atmosphere_hmin)/&
        dble(n_atmosphere-1)
    gravity0 = -2.74d4*unit_length/unit_velocity**2
    reference_height = rho_ref_height_cm/unit_length

    do j=1,n_atmosphere
      h(j) = atmosphere_hmin+dble(j-1)*atmosphere_dh
      grav(j) = gravity0*(radial_center0/(radial_center0+h(j)))**2
    end do

    if(allocated(atmosphere_rho)) deallocate(atmosphere_rho)
    if(allocated(atmosphere_p)) deallocate(atmosphere_p)
    if(allocated(atmosphere_T)) deallocate(atmosphere_T)
    allocate(atmosphere_rho(n_atmosphere),atmosphere_p(n_atmosphere),&
        atmosphere_T(n_atmosphere))
    call get_atm_para(h,atmosphere_rho,atmosphere_p,grav,n_atmosphere,&
        trim(atmosphere_curve),reference_height,rho_ref,&
        Tem=atmosphere_T,clamp_low_T=.true.)
  end subroutine initialize_solar_atmosphere

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

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bvec(ixI^S,1:ndim)
    logical :: prom(ixI^S)

    if(.not.tdm_setup_ready) then
      call mpistop('TDm field initialization called after restart setup skip')
    end if

    if(stagger_grid) then
      call mpistop('RadiationSynthesis_Sph dip mask requires cell-centered B')
    end if

    ! Use exactly the direct cell-centered magnetic-field path of the topology
    ! case, including guard cells required by the magnetic-dip derivative.
    call tdm_magnetic_field_spherical(ixI^L,ixI^L,x,Bvec)
    w(ixI^S,mag(:)) = Bvec(ixI^S,:)

    call set_atmosphere_on_grid(ixI^L,ixO^L,x,w)
    w(ixO^S,mom(:)) = zero
    if(mhd_glm) w(ixO^S,psi_) = zero

    prom = .false.
    call get_prominence_mask(ixI^L,ixO^L,w,x,prom)
    where(prom(ixO^S))
      w(ixO^S,rho_) = prom_density_factor*w(ixO^S,rho_)
    end where

    call eos%to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine set_atmosphere_on_grid(ixI^L,ixO^L,x,w)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: position,fraction,height
    integer :: ix^D,j

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      height = x(ix^D,1)-radial_center0
      position = (height-atmosphere_hmin)/atmosphere_dh
      if(position<=zero) then
        j = 1
        fraction = zero
      else if(position>=dble(n_atmosphere-1)) then
        j = n_atmosphere-1
        fraction = one
      else
        j = int(floor(position))+1
        fraction = position-dble(j-1)
      end if
      w(ix^D,rho_) = atmosphere_rho(j)+fraction*&
          (atmosphere_rho(j+1)-atmosphere_rho(j))
      w(ix^D,p_) = atmosphere_p(j)+fraction*&
          (atmosphere_p(j+1)-atmosphere_p(j))
    {end do\}
  end subroutine set_atmosphere_on_grid

  subroutine get_prominence_mask(ixI^L,ixO^L,w,x,prom)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: w(ixI^S,1:*)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    logical, intent(out) :: prom(ixI^S)
    double precision :: grad_br(ixI^S),dip_metric(ixI^S)
    double precision :: xcart(ixI^S,1:ndim)
    double precision :: xloc(ixI^S),yloc(ixI^S),zloc(ixI^S)
    double precision :: radial_distance(ixI^S),axis_distance(ixI^S)
    integer :: idir

    dip_metric = zero
    do idir=1,ndir
      call gradient(w(ixI^S,mag(1)),ixI^L,ixO^L,idir,grad_br)
      dip_metric(ixO^S) = dip_metric(ixO^S)+&
          w(ixO^S,mag(idir))*grad_br(ixO^S)
    end do

    call spherical_to_local_reference(ixI^L,ixO^L,x,xcart,xloc,yloc,zloc)
    radial_distance(ixO^S) = sqrt(xloc(ixO^S)**2+&
        (zloc(ixO^S)+4.5d0)**2)
    axis_distance(ixO^S) = sqrt(yloc(ixO^S)**2+&
        (radial_distance(ixO^S)-8.d0)**2)

    prom = .false.
    where(axis_distance(ixO^S)<=prom_axis_fraction*a0 .and.&
        abs(w(ixO^S,mag(1)))<=dip_bz_max .and.&
        dip_metric(ixO^S)>=zero)
      prom(ixO^S) = .true.
    end where
  end subroutine get_prominence_mask

  subroutine spherical_to_local_reference(ixI^L,ixO^L,x,xcart,xloc,yloc,&
      zloc)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: xcart(ixI^S,1:ndim)
    double precision, intent(out) :: xloc(ixI^S),yloc(ixI^S),zloc(ixI^S)

    xcart(ixO^S,1) = x(ixO^S,1)*sin(x(ixO^S,2))*cos(x(ixO^S,3))
    xcart(ixO^S,2) = x(ixO^S,1)*sin(x(ixO^S,2))*sin(x(ixO^S,3))
    xcart(ixO^S,3) = x(ixO^S,1)*cos(x(ixO^S,2))
    xloc(ixO^S) = xcart(ixO^S,3)
    yloc(ixO^S) = xcart(ixO^S,2)
    zloc(ixO^S) = xcart(ixO^S,1)-radial_center0
  end subroutine spherical_to_local_reference

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: w(ixI^S,nw+nwauxio)
    double precision :: normconv(0:nw+nwauxio)
    double precision :: pth(ixI^S)
    logical :: prom(ixI^S)

    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)
    call get_prominence_mask(ixI^L,ixO^L,w,x,prom)
    w(ixO^S,nw+1) = pth(ixO^S)/w(ixO^S,rho_)
    w(ixO^S,nw+2) = zero
    where(prom(ixO^S)) w(ixO^S,nw+2) = one
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames = 'T prominence_mask'
  end subroutine specialvarnames_output

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
