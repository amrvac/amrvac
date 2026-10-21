!> Cartesian TDm/RBSL field used by the finite-volume magnetic-helicity demo.
!> Both cases use the same local Cartesian magnetic field, analytic circular
!> axis, normalization, background bipolar field, and RBSL kernels.
module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  use mod_global_parameters, only: par_files,mype,&
     restart_from_file,undefined,firstprocess
  use mod_usr_methods, only: usr_refine_grid
  implicit none

  integer :: np
  double precision :: q_para,d_para,L_para
  double precision :: a0,F_flx
  double precision :: minor_radius_cm,rho0,pressure0
  integer :: n_axis_points
  double precision, allocatable :: x_axis(:,:)
  logical :: tdm_setup_ready=.false.

contains

  subroutine usr_init()
    unit_length        = 1.d9
    unit_numberdensity = 1.d9
    unit_temperature   = 1.d6

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_refine_grid => tdm_bench_refine_grid

    call set_coordinate_system('Cartesian_3D')
    call mhd_activate()
  end subroutine usr_init

  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /usr_list/ minor_radius_cm,n_axis_points,rho0,pressure0

    minor_radius_cm = 2.25d9
    n_axis_points = 400
    rho0 = 1.d0
    pressure0 = 1.d0

    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,usr_list,end=111)
111   close(unitpar)
    end do
  end subroutine usr_params_read

  subroutine initglobaldata_usr()
    double precision, parameter :: major_radius = 8.d0
    double precision :: b_perp_apex,shafranov_factor,mu0I_equilibrium

    call usr_params_read(par_files)

    a0 = minor_radius_cm/unit_length
    if(a0<=zero .or. a0>=major_radius) then
      call mpistop('TDm requires 0 < minor_radius_cm < 8e9 cm')
    end if

    ! A normal restart reads the magnetic field from the snapshot.  Do not
    ! construct the analytic axis or any TDm/RBSL setup data unless the user
    ! explicitly requests firstprocess, which reruns usr_init_one_grid.
    tdm_setup_ready = .false.
    if(restart_from_file/=undefined .and. .not.firstprocess) then
      if(allocated(x_axis)) deallocate(x_axis)
      if(mype==0) print *,&
         'Restart mode: skipping Cartesian TDm/RBSL field setup'
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
    tdm_setup_ready = .true.

    if(mype==0) then
      print *, 'Cartesian TDm/RBSL initial condition'
      print *, 'Using one complete circular TDm/RBSL integration path'
      print *, 'unit_length [cm]: ',unit_length
      print *, 'q_para normalized: ',q_para
      print *, 'apex strapping-field magnitude: ',b_perp_apex
      print *, 'equilibrium mu0*I normalized: ',mu0I_equilibrium
      print *, 'F_flx normalized: ',F_flx
      print *, 'source depth [code units]: ',d_para
      print *, 'source half separation [code units]: ',L_para
      print *, 'major radius [code units]: ',major_radius
      print *, 'minor radius [code units]: ',a0
      print *, 'axis integration points: ',np
    end if
  end subroutine initglobaldata_usr

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

    if(.not.tdm_setup_ready) then
      call mpistop('TDm field initialization called after restart setup skip')
    end if

    B_bg = zero
    B_rope = zero
    q_direct = -q_para/3.d0

    call bipolar_field_direct_B(ixI^L,ixO^L,x,L_para,d_para,q_direct,zero,3,&
       B_bg)
    call rbsl_flux_rope_direct_B(ixI^L,ixO^L,np,a0,F_flx,.false.,x,x_axis,&
       B_rope)

    w(ixO^S,rho_) = rho0
    w(ixO^S,mom(:)) = zero
    w(ixO^S,p_) = pressure0
    w(ixO^S,mag(:)) = B_bg(ixO^S,:)+B_rope(ixO^S,:)
    if(mhd_glm) w(ixO^S,psi_) = zero

    call eos%to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

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
