!> Cartesian TDm/RBSL radiation-synthesis demonstration.
!> The magnetic field and AMR regions are identical to MagneticTopologyQSL_Cart.
!> The thermodynamic atmosphere is an AL-C7 hydrostatic atmosphere, and cool
!> prominence material is placed only in magnetic dips inside the flux rope.
module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  use mod_solar_atmosphere, only: get_atm_para
  use mod_global_parameters, only: par_files,mype,refine_max_level,&
     restart_from_file,undefined,firstprocess
  use mod_usr_methods, only: usr_set_parameters,usr_init_one_grid,&
     usr_refine_grid,usr_aux_output,usr_add_aux_names
  implicit none

  integer, parameter :: n_atmosphere=4096
  integer :: np
  double precision :: q_para,d_para,L_para
  double precision :: a0,F_flx
  double precision :: minor_radius_cm,rho_ref_height_cm,rho_ref
  double precision :: prom_density_factor,dip_bz_max,prom_axis_fraction
  integer :: n_axis_points
  character(len=16) :: atmosphere_curve
  double precision :: atmosphere_hmin,atmosphere_dh
  double precision, allocatable :: x_axis(:,:)
  double precision, allocatable :: atmosphere_rho(:),atmosphere_p(:),&
     atmosphere_T(:)
  logical :: tdm_setup_ready=.false.

contains

  subroutine usr_init()
    unit_length        = 1.d9
    unit_numberdensity = 1.d9
    unit_temperature   = 1.d6

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_refine_grid => tdm_bench_refine_grid
    usr_aux_output => specialvar_output
    usr_add_aux_names => specialvarnames_output

    call set_coordinate_system('Cartesian_3D')
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

    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,usr_list,end=111)
111   close(unitpar)
    end do
  end subroutine usr_params_read

  subroutine initglobaldata_usr()
    double precision, parameter :: major_radius=8.d0
    double precision :: b_perp_apex,shafranov_factor,mu0I_equilibrium

    call usr_params_read(par_files)

    a0 = minor_radius_cm/unit_length
    if(a0<=zero .or. a0>=major_radius) then
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

    ! Conversion/diagnostic restarts consume the magnetic and thermodynamic
    ! state stored in the snapshot.  Avoid constructing the analytic TDm axis
    ! and the hydrostatic atmosphere table unless firstprocess explicitly asks
    ! AMRVAC to regenerate the initial condition.
    tdm_setup_ready = .false.
    if(restart_from_file/=undefined .and. .not.firstprocess) then
      if(allocated(x_axis)) deallocate(x_axis)
      if(allocated(atmosphere_rho)) deallocate(atmosphere_rho)
      if(allocated(atmosphere_p)) deallocate(atmosphere_p)
      if(allocated(atmosphere_T)) deallocate(atmosphere_T)
      if(mype==0) print *,&
         'Restart mode: skipping Cartesian TDm/atmosphere setup'
      return
    end if

    d_para = 4.5d0
    L_para = 3.d0
    q_para = -600.d0/sqrt(4.d0*dpi)

    b_perp_apex = abs(2.d0*L_para*q_para/&
       (major_radius**2+L_para**2)**1.5d0)
    shafranov_factor = log(8.d0*major_radius/a0)-25.d0/24.d0
    if(shafranov_factor<=zero) then
      call mpistop('Invalid TDm Shafranov equilibrium factor')
    end if
    mu0I_equilibrium = 4.d0*dpi*major_radius*b_perp_apex/shafranov_factor
    F_flx = 3.d0*mu0I_equilibrium*a0/(5.d0*sqrt(2.d0))

    np = max(16,n_axis_points)
    if(allocated(x_axis)) deallocate(x_axis)
    allocate(x_axis(np,ndim))
    call calc_cartesian_tdm_axis(x_axis,np)
    call initialize_solar_atmosphere()
    tdm_setup_ready = .true.

    if(mype==0) then
      print *, 'Cartesian TDm/RBSL radiation-synthesis initial condition'
      print *, 'Using the MagneticTopologyQSL_Cart analytic magnetic field'
      print *, 'unit_length [cm]: ',unit_length
      print *, 'q_para normalized: ',q_para
      print *, 'apex strapping-field magnitude: ',b_perp_apex
      print *, 'equilibrium mu0*I normalized: ',mu0I_equilibrium
      print *, 'F_flx normalized: ',F_flx
      print *, 'major radius [code units]: ',major_radius
      print *, 'minor radius [code units]: ',a0
      print *, 'axis integration points: ',np
      print *, 'atmosphere curve: ',trim(atmosphere_curve)
      print *, 'reference height [Mm]: ',rho_ref_height_cm/1.d8
      print *, 'reference number density [cm^-3]: ',rho_ref*unit_numberdensity
      print *, 'prominence density factor: ',prom_density_factor
      print *, 'dip |Bz| threshold: ',dip_bz_max
      print *, 'dip axis-distance limit [code units]: ',prom_axis_fraction*a0
    end if
  end subroutine initglobaldata_usr

  subroutine initialize_solar_atmosphere()
    double precision, parameter :: solar_radius_cm=6.961d10
    double precision :: h(n_atmosphere),grav(n_atmosphere)
    double precision :: atmosphere_hmax,gravity0,solar_radius
    double precision :: reference_height,padding
    integer :: j

    padding = 2.d8/unit_length
    atmosphere_hmin = min(zero,xprobmin3)-padding
    atmosphere_hmax = xprobmax3+padding
    atmosphere_dh = (atmosphere_hmax-atmosphere_hmin)/&
       dble(n_atmosphere-1)
    gravity0 = -2.74d4*unit_length/unit_velocity**2
    solar_radius = solar_radius_cm/unit_length
    reference_height = rho_ref_height_cm/unit_length

    do j=1,n_atmosphere
      h(j) = atmosphere_hmin+dble(j-1)*atmosphere_dh
      grav(j) = gravity0*(solar_radius/(solar_radius+h(j)))**2
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

  subroutine calc_cartesian_tdm_axis(xs,npo)
    integer, intent(in) :: npo
    double precision, intent(out) :: xs(npo,3)
    integer :: i
    double precision :: theta_start,theta
    double precision, parameter :: center_x=0.d0,center_z=-4.5d0
    double precision, parameter :: major_radius=8.d0

    theta_start = asin(-center_z/major_radius)
    do i=1,npo
      theta = theta_start+2.d0*dpi*dble(i-1)/dble(npo)
      xs(i,1) = center_x+major_radius*cos(theta)
      xs(i,2) = zero
      xs(i,3) = center_z+major_radius*sin(theta)
    end do
  end subroutine calc_cartesian_tdm_axis

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: B_bg(ixI^S,1:ndim),B_rope(ixI^S,1:ndim)
    double precision :: q_direct
    logical :: prom(ixI^S)

    if(.not.tdm_setup_ready) then
      call mpistop('TDm field initialization called after restart setup skip')
    end if

    B_bg = zero
    B_rope = zero
    q_direct = -q_para/3.d0

    ! Evaluate analytic B also in guard cells so the dip derivative is local
    ! and independent of neighboring-block initialization order.
    call bipolar_field_direct_B(ixI^L,ixI^L,x,L_para,d_para,q_direct,zero,3,&
       B_bg)
    call rbsl_flux_rope_direct_B(ixI^L,ixI^L,np,a0,F_flx,.false.,x,x_axis,&
       B_rope)
    w(ixI^S,mag(:)) = B_bg(ixI^S,:)+B_rope(ixI^S,:)

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
    double precision :: position,fraction
    integer :: ix^D,j

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      position = (x(ix^D,3)-atmosphere_hmin)/atmosphere_dh
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

  subroutine get_prominence_mask(ixI^L,ixO^L,w,x,prom,dip_metric_out,&
     axis_distance_out)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: w(ixI^S,1:*)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    logical, intent(out) :: prom(ixI^S)
    double precision, optional, intent(out) :: dip_metric_out(ixI^S)
    double precision, optional, intent(out) :: axis_distance_out(ixI^S)
    double precision :: grad_bz(ixI^S),dip_metric(ixI^S)
    double precision :: radial_distance(ixI^S),axis_distance(ixI^S)
    integer :: idir

    dip_metric = zero
    do idir=1,ndir
      call gradient(w(ixI^S,mag(3)),ixI^L,ixO^L,idir,grad_bz)
      dip_metric(ixO^S) = dip_metric(ixO^S)+&
         w(ixO^S,mag(idir))*grad_bz(ixO^S)
    end do

    radial_distance(ixO^S) = sqrt(x(ixO^S,1)**2+&
       (x(ixO^S,3)+4.5d0)**2)
    axis_distance(ixO^S) = sqrt(x(ixO^S,2)**2+&
       (radial_distance(ixO^S)-8.d0)**2)

    prom = .false.
    where(axis_distance(ixO^S)<=prom_axis_fraction*a0 .and.&
       abs(w(ixO^S,mag(3)))<=dip_bz_max .and.&
       dip_metric(ixO^S)>=zero)
      prom(ixO^S) = .true.
    end where

    if(present(dip_metric_out)) then
      dip_metric_out(ixO^S) = dip_metric(ixO^S)
    end if
    if(present(axis_distance_out)) then
      axis_distance_out(ixO^S) = axis_distance(ixO^S)
    end if
  end subroutine get_prominence_mask

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

  subroutine tdm_bench_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,&
     coarsen)
    integer, intent(in) :: igrid,level,ixI^L,ixO^L
    double precision, intent(in) :: qt,w(ixI^S,1:nw),x(ixI^S,1:ndim)
    integer, intent(inout) :: refine,coarsen
    logical :: broad_region,rope_region,core_region

    broad_region = any(abs(x(ixO^S,1))<=8.5d0 .and.&
       abs(x(ixO^S,2))<=4.5d0 .and. x(ixO^S,3)>=zero .and.&
       x(ixO^S,3)<=7.5d0)
    rope_region = any(abs(x(ixO^S,1))<=7.5d0 .and.&
       abs(x(ixO^S,2))<=2.5d0 .and. x(ixO^S,3)>=zero .and.&
       x(ixO^S,3)<=5.5d0) .or.&
       any(abs(x(ixO^S,1))<=9.d0 .and. abs(x(ixO^S,2))<=4.d0 .and.&
       x(ixO^S,3)>=zero .and. x(ixO^S,3)<=2.d0)
    core_region = any(abs(x(ixO^S,1))<=5.5d0 .and.&
       abs(x(ixO^S,2))<=1.8d0 .and. x(ixO^S,3)>=0.5d0 .and.&
       x(ixO^S,3)<=4.5d0)

    if(level==1 .and. broad_region) then
      refine = 1
      coarsen = -1
    else if(level==2 .and. rope_region) then
      refine = 1
      coarsen = -1
    else if(level==3 .and. core_region) then
      refine = 1
      coarsen = -1
    end if
  end subroutine tdm_bench_refine_grid

  subroutine bipolar_field_direct_B(ixI^L,ixO^L,x,L_cha,d_cha,q_cha,x_cha,&
     nb,Bout)
    integer, intent(in) :: ixI^L,ixO^L,nb
    double precision, intent(in) :: L_cha,d_cha,q_cha,x_cha
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: Bout(ixI^S,1:ndim)
    integer :: i,idir
    double precision :: xpos
    double precision :: rpv(ixI^S),rmv(ixI^S)
    double precision :: rplus(ixI^S,1:ndim),rminus(ixI^S,1:ndim)

    Bout = zero
    do i=1,nb
      xpos = 2.d0*dble(i-1)*x_cha/dble(max(1,nb-1))-x_cha
      rplus(ixO^S,1) = x(ixO^S,1)-xpos
      rminus(ixO^S,1) = x(ixO^S,1)-xpos
      rplus(ixO^S,2) = x(ixO^S,2)-L_cha
      rminus(ixO^S,2) = x(ixO^S,2)+L_cha
      rplus(ixO^S,3) = x(ixO^S,3)+d_cha
      rminus(ixO^S,3) = x(ixO^S,3)+d_cha
      rpv(ixO^S) = sqrt(sum(rplus(ixO^S,:)**2,dim=ndim+1))
      rmv(ixO^S) = sqrt(sum(rminus(ixO^S,:)**2,dim=ndim+1))
      do idir=1,ndim
        Bout(ixO^S,idir) = Bout(ixO^S,idir)+q_cha*&
           (rplus(ixO^S,idir)/rpv(ixO^S)**3-&
            rminus(ixO^S,idir)/rmv(ixO^S)**3)
      end do
    end do
  end subroutine bipolar_field_direct_B

  subroutine rbsl_flux_rope_direct_B(ixI^L,ixO^L,np,a,F_flx,&
     positive_helicity,x,x_axis,Btotal)
    integer, intent(in) :: ixI^L,ixO^L,np
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(in) :: x_axis(np,1:ndim),a,F_flx
    logical, intent(in) :: positive_helicity
    double precision, intent(out) :: Btotal(ixI^S,1:ndim)
    double precision :: I_cur
    double precision :: BIx(ixI^S,1:ndim),BFx(ixI^S,1:ndim)
    double precision :: axis_element(np,1:ndim)
    double precision :: r_mag,KIr,KFr1,KFr2,Rdr
    double precision :: asr,or2,asrr,re_pi
    double precision :: Rpl(1:ndim),r_vec(1:ndim),Rcr(1:ndim)
    integer :: ix^D,ixp

    if(positive_helicity) then
      I_cur = 5.d0*sqrt(2.d0)*F_flx/(3.d0*a)
    else
      I_cur = -5.d0*sqrt(2.d0)*F_flx/(3.d0*a)
    end if

    re_pi = 1.d0/dpi
    BIx = zero
    BFx = zero
    do ixp=1,np
      if(ixp==1) then
        axis_element(ixp,:) = 0.5d0*(x_axis(ixp+1,:)-x_axis(np,:))
      else if(ixp==np) then
        axis_element(ixp,:) = 0.5d0*(x_axis(1,:)-x_axis(ixp-1,:))
      else
        axis_element(ixp,:) = 0.5d0*(x_axis(ixp+1,:)-x_axis(ixp-1,:))
      end if
    end do

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      do ixp=1,np
        r_vec(:) = (x(ix^D,:)-x_axis(ixp,:))/a
        r_mag = sqrt(sum(r_vec(:)**2))
        Rpl(:) = axis_element(ixp,:)
        Rcr(1) = Rpl(2)*r_vec(3)-Rpl(3)*r_vec(2)
        Rcr(2) = Rpl(3)*r_vec(1)-Rpl(1)*r_vec(3)
        Rcr(3) = Rpl(1)*r_vec(2)-Rpl(2)*r_vec(1)
        Rdr = sum(Rpl(:)*r_vec(:))
        if(r_mag<=1.d-3) then
          KIr = 16.d0/(3.d0*dpi)
          KFr1 = 5.d0/sqrt(6.d0)+10.d0/dpi*&
             (2.d0/3.d0-asin(0.2d0)/sqrt(6.d0))
          KFr2 = sqrt(6.d0)/3.d0+2.d0/(15.d0*dpi)*&
             (24.d0-5.d0*sqrt(6.d0)*asin(0.2d0))
        else if(r_mag<=one) then
          asr = asin(r_mag)/r_mag
          or2 = sqrt(one-r_mag**2)
          asrr = asin((one+two*r_mag*r_mag)/(5.d0-two*r_mag*r_mag))
          KIr = two*re_pi*((asr-or2)/r_mag**2+two*or2)
          KFr1 = two*re_pi/r_mag**2*(or2-asr)+8.d0*re_pi*or2+&
             (5.d0-4.d0*r_mag**2)/sqrt(6.d0)*(one-two*re_pi*asrr)
          KFr2 = two*re_pi/r_mag**4*&
             (3.d0*asr-(3.d0+two*r_mag**2)*or2)+&
             two/sqrt(6.d0)*(one-two*re_pi*asrr)
        else
          KIr = one/r_mag**3
          KFr1 = -one/r_mag**3
          KFr2 = 3.d0/r_mag**5
        end if
        BIx(ix^D,:) = BIx(ix^D,:)+I_cur*0.25d0*re_pi*KIr*Rcr(:)/a**2
        BFx(ix^D,:) = BFx(ix^D,:)+F_flx*0.25d0*re_pi*&
           (KFr1*Rpl(:)+KFr2*Rdr*r_vec(:))/a**3
      end do
    {end do\}

    Btotal(ixO^S,:) = BIx(ixO^S,:)+BFx(ixO^S,:)
  end subroutine rbsl_flux_rope_direct_B

end module mod_usr
