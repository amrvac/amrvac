!> Demcsak et al. (2020) TD99 boundary generator and Grad--Rubin benchmark.
module mod_usr
  use mod_mhd
  use mod_nlfff_grad_rubin
  use mod_tdfluxrope
  implicit none

  character(len=32) :: run_mode
  character(len=256) :: boundary_filename,gr_alpha_filename
  character(len=256) :: boundary_export_filename,alpha_export_filename
  character(len=24) :: gr_alpha_source
  integer :: gr_fft_padding_factor,gr_polarity,gr_fieldline_max_steps
  integer :: gr_max_iterations,gr_convergence_streak,gr_log_interval
  integer :: gr_self_consistency_cycles,gr_self_consistency_streak
  character(len=16) :: gr_flux_treatment
  double precision :: i0_effective_ta,alpha_step_scale
  double precision :: gr_max_flux_imbalance,gr_bz_taper_zero,gr_bz_taper_full
  double precision :: gr_relaxation_factor,gr_fieldline_step_fraction
  double precision :: gr_field_change_tolerance,gr_energy_change_tolerance
  double precision :: gr_self_consistency_field_tolerance
  double precision :: gr_self_consistency_alpha_tolerance
  double precision :: gr_memory_limit_mb
  logical :: gr_write_detailed_history
  type(nlfff_grad_rubin_result), save :: gr_result
  logical, save :: first_parameters=.true.

contains

  subroutine usr_init()
    use mod_usr_methods

    unit_length        = 5.d9
    unit_temperature   = 1.d6
    unit_numberdensity = 1.d10
    call set_coordinate_system('Cartesian_3D')
    usr_set_parameters            => initglobaldata_usr
    usr_init_one_grid             => initonegrid_usr
    usr_improve_initial_condition => improve_initial_condition_usr
    call mhd_activate()
  end subroutine usr_init

  subroutine usr_params_read(files)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    character(len=*), intent(in) :: files(:)
    integer :: n
    namelist /usr_list/ run_mode,boundary_filename,gr_alpha_filename,&
         boundary_export_filename,alpha_export_filename,&
         i0_effective_ta,alpha_step_scale,gr_alpha_source,&
         gr_fft_padding_factor,&
         gr_flux_treatment,gr_max_flux_imbalance,gr_polarity,&
         gr_bz_taper_zero,gr_bz_taper_full,gr_relaxation_factor,&
         gr_fieldline_step_fraction,gr_fieldline_max_steps,&
         gr_max_iterations,gr_convergence_streak,&
         gr_field_change_tolerance,gr_energy_change_tolerance,&
         gr_self_consistency_cycles,gr_self_consistency_streak,&
         gr_self_consistency_field_tolerance,&
         gr_self_consistency_alpha_tolerance,&
         gr_log_interval,gr_memory_limit_mb,gr_write_detailed_history

    run_mode='grad_rubin'
    boundary_filename=''
    gr_alpha_filename=''
    boundary_export_filename=''
    alpha_export_filename=''
    i0_effective_ta=13.d0
    alpha_step_scale=0.25d0
    gr_alpha_source='external'
    gr_fft_padding_factor=1
    gr_flux_treatment='strict'
    gr_max_flux_imbalance=0.1d0
    gr_polarity=-1
    gr_bz_taper_zero=0.01d0
    gr_bz_taper_full=0.02d0
    gr_relaxation_factor=0.5d0
    gr_fieldline_step_fraction=0.5d0
    gr_fieldline_max_steps=10000
    gr_max_iterations=16
    gr_convergence_streak=3
    gr_field_change_tolerance=0.d0
    gr_energy_change_tolerance=0.d0
    gr_self_consistency_cycles=0
    gr_self_consistency_streak=1
    gr_self_consistency_field_tolerance=1.d-3
    gr_self_consistency_alpha_tolerance=1.d-3
    gr_log_interval=1
    gr_memory_limit_mb=2048.d0
    gr_write_detailed_history=.true.
    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,usr_list,end=111)
111   close(unitpar)
    end do

    select case(trim(adjustl(run_mode)))
    case('generate_boundary')
      if(len_trim(boundary_export_filename)==0 .or. &
         len_trim(alpha_export_filename)==0) &
         call mpistop('TD99 boundary and alpha export filenames must be set')
      if(alpha_step_scale<=0.d0) &
         call mpistop('TD99 alpha_step_scale must be positive')
    case('grad_rubin')
      if(len_trim(boundary_filename)==0) &
         call mpistop('Grad-Rubin boundary filename must be set')
      if(trim(adjustl(gr_alpha_source))/='vector_magnetogram' .and. &
         trim(adjustl(gr_alpha_source))/='external') &
         call mpistop("gr_alpha_source must be 'vector_magnetogram' or 'external'")
      if(trim(adjustl(gr_alpha_source))=='external' .and. &
         len_trim(gr_alpha_filename)==0) &
         call mpistop('external Grad-Rubin alpha filename must be set')
    case default
      call mpistop("run_mode must be 'generate_boundary' or 'grad_rubin'")
    end select
  end subroutine usr_params_read

  subroutine initglobaldata_usr()
    use mod_global_parameters

    call usr_params_read(par_files)
    if(len_trim(base_filename)==0) then
      if(trim(adjustl(run_mode))=='generate_boundary') then
        base_filename='grad_rubin_td99_boundary'
      else
        base_filename='grad_rubin_td99'
      end if
    end if
    call configure_td99()
    if(trim(adjustl(run_mode))=='grad_rubin' .and. first_parameters) then
      call init_nlfff_grad_rubin_boundary(trim(boundary_filename),&
           unit_length,unit_magneticfield)
      if(trim(adjustl(gr_alpha_source))=='external') &
         call init_nlfff_grad_rubin_external_alpha(trim(gr_alpha_filename),&
              unit_length,unit_magneticfield)
      first_parameters=.false.
    end if
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
      write(*,*) 'Demcsak et al. TD99 parameters R,a,L,d [Mm]:',&
           R_TD99*unit_length/1.d8,a_TD99*unit_length/1.d8,&
           L_TD99*unit_length/1.d8,d_TD99*unit_length/1.d8
      write(*,*) 'I0 [TA], surface twist:',i0_effective_ta,nt_td99
    end if
  end subroutine configure_td99

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,:)=zero
    w(ixO^S,rho_)=one
  end subroutine initonegrid_usr

  subroutine improve_initial_condition_usr()
    use mod_global_parameters
    use mod_ghostcells_update, only: getbc

    type(nlfff_grad_rubin_config) :: config

    if(trim(adjustl(run_mode))=='generate_boundary') then
      call export_td99_products()
      return
    end if

    config%fft_padding_factor=gr_fft_padding_factor
    config%flux_treatment=trim(gr_flux_treatment)
    config%max_flux_imbalance=gr_max_flux_imbalance
    config%polarity=gr_polarity
    config%alpha_source=trim(adjustl(gr_alpha_source))
    config%bz_taper_zero=gr_bz_taper_zero
    config%bz_taper_full=gr_bz_taper_full
    config%relaxation_factor=gr_relaxation_factor
    config%fieldline_step_fraction=gr_fieldline_step_fraction
    config%fieldline_max_steps=gr_fieldline_max_steps
    config%max_iterations=gr_max_iterations
    config%convergence_streak=gr_convergence_streak
    config%field_change_tolerance=gr_field_change_tolerance
    config%energy_change_tolerance=gr_energy_change_tolerance
    config%self_consistency_cycles=gr_self_consistency_cycles
    config%self_consistency_streak=gr_self_consistency_streak
    config%self_consistency_field_tolerance=&
         gr_self_consistency_field_tolerance
    config%self_consistency_alpha_tolerance=&
         gr_self_consistency_alpha_tolerance
    config%log_interval=gr_log_interval
    config%memory_limit_mb=gr_memory_limit_mb
    config%write_detailed_history=gr_write_detailed_history
    call extrapolate_nlfff_grad_rubin(mag(:),config,gr_result)
    if(mype==0) then
      write(*,*) 'Grad-Rubin stop reason:',trim(gr_result%stop_reason)
      write(*,*) 'Grad-Rubin iterations/converged:',gr_result%iterations,&
           gr_result%converged
    end if
    call getbc(global_time,0.d0,ps,iwstart,nwgc)
  end subroutine improve_initial_condition_usr

  subroutine export_td99_products()
    use mod_global_parameters

    integer :: ixI^L,ixO^L,ix1,ix2,nx,ny,iu
    double precision :: dx_km,dy_km
    double precision, allocatable :: xplane(:,:,:,:),bplane(:,:,:,:)
    double precision, allocatable :: jz(:,:),alpha_raw(:,:),alpha(:,:)
    double precision, allocatable :: weight(:,:),xcoord(:),ycoord(:)
    integer, allocatable :: pil_mask(:,:),valid_mask(:,:),polarity_mask(:,:)
    logical, allocatable :: rope_mask(:,:)

    if(mype/=0) return
    nx=domain_nx1
    ny=domain_nx2
    allocate(xplane(1:nx,1:ny,1:1,1:ndim))
    allocate(bplane(1:nx,1:ny,1:1,1:ndir))
    allocate(jz(nx,ny),alpha_raw(nx,ny),alpha(nx,ny),weight(nx,ny))
    allocate(xcoord(nx),ycoord(ny))
    allocate(pil_mask(nx,ny),valid_mask(nx,ny),polarity_mask(nx,ny))
    allocate(rope_mask(nx,ny))
    do ix2=1,ny
      ycoord(ix2)=xprobmin2+(dble(ix2)-0.5d0)*dx(2,1)
      do ix1=1,nx
        xcoord(ix1)=xprobmin1+(dble(ix1)-0.5d0)*dx(1,1)
        xplane(ix1,ix2,1,1)=xcoord(ix1)
        xplane(ix1,ix2,1,2)=ycoord(ix2)
        xplane(ix1,ix2,1,3)=0.d0
      end do
    end do
    ixImin1=1; ixImax1=nx
    ixImin2=1; ixImax2=ny
    ixImin3=1; ixImax3=1
    ixOmin1=ixImin1; ixOmax1=ixImax1
    ixOmin2=ixImin2; ixOmax2=ixImax2
    ixOmin3=ixImin3; ixOmax3=ixImax3
    call TD99(ixI^L,ixO^L,xplane,bplane)

    dx_km=dx(1,1)*unit_length/1.d5
    dy_km=dx(2,1)*unit_length/1.d5
    open(newunit=iu,file=trim(boundary_export_filename),status='replace',&
         access='stream',form='unformatted',action='write')
    write(iu) zero,nx,ny,dx_km,dy_km
    write(iu) bplane(:,:,1,:)*unit_magneticfield
    close(iu)

    call compute_td99_jz(nx,ny,alpha_step_scale*dx(1,1),jz,rope_mask)
    alpha_raw=zero
    alpha=zero
    weight=zero
    pil_mask=0
    valid_mask=0
    polarity_mask=0
    do ix2=1,ny
      do ix1=1,nx
        if(abs(bplane(ix1,ix2,1,3))>epsilon(one)) then
          alpha_raw(ix1,ix2)=jz(ix1,ix2)/bplane(ix1,ix2,1,3)
        end if
        if(rope_mask(ix1,ix2) .and. &
           abs(bplane(ix1,ix2,1,3))>epsilon(one)) then
          alpha(ix1,ix2)=alpha_raw(ix1,ix2)
          weight(ix1,ix2)=one
          pil_mask(ix1,ix2)=1
          valid_mask(ix1,ix2)=1
          if(bplane(ix1,ix2,1,3)>zero) then
            polarity_mask(ix1,ix2)=1
          else
            polarity_mask(ix1,ix2)=-1
          end if
        end if
      end do
    end do
    call write_external_alpha(nx,ny,xcoord,ycoord,alpha_raw,alpha,weight,&
         pil_mask,valid_mask,polarity_mask)
    write(*,*) 'Exported TD99 boundary/alpha:',trim(boundary_export_filename),&
         trim(alpha_export_filename)
    write(*,*) 'TD99 alpha support/min/max:',count(valid_mask/=0),&
         minval(alpha,mask=valid_mask/=0),maxval(alpha,mask=valid_mask/=0)

    deallocate(xplane,bplane,jz,alpha_raw,alpha,weight,xcoord,ycoord)
    deallocate(pil_mask,valid_mask,polarity_mask,rope_mask)
  end subroutine export_td99_products

  subroutine write_external_alpha(nx,ny,xcoord,ycoord,alpha_raw,alpha,&
       weight,pil_mask,valid_mask,polarity_mask)
    use mod_global_parameters

    integer, intent(in) :: nx,ny
    double precision, intent(in) :: xcoord(nx),ycoord(ny)
    double precision, intent(in) :: alpha_raw(nx,ny),alpha(nx,ny)
    double precision, intent(in) :: weight(nx,ny)
    integer, intent(in) :: pil_mask(nx,ny),valid_mask(nx,ny)
    integer, intent(in) :: polarity_mask(nx,ny)
    character(len=32) :: magic
    character(len=4096) :: metadata_blob
    integer :: iu,version,reserved,alpha_code,mask_code,metadata_length
    double precision :: xc,yc

    magic='AMRVAC_EXTERNAL_ALPHA_V1'
    version=1
    reserved=0
    alpha_code=1
    mask_code=1
    metadata_length=0
    metadata_blob=''
    xc=0.5d0*(xcoord(1)+xcoord(nx))
    yc=0.5d0*(ycoord(1)+ycoord(ny))
    open(newunit=iu,file=trim(alpha_export_filename),status='replace',&
         access='stream',form='unformatted',action='write')
    write(iu) magic,version,nx,ny,reserved,unit_length,unit_magneticfield,&
         dx(1,1),dx(2,1),xc,yc,alpha_code,mask_code
    write(iu) xcoord,ycoord
    write(iu) alpha_raw,alpha,weight,pil_mask,valid_mask,polarity_mask
    write(iu) metadata_length,metadata_blob
    close(iu)
  end subroutine write_external_alpha

  subroutine compute_td99_jz(nx,ny,hnom,jz,rope_mask)
    use mod_global_parameters

    integer, intent(in) :: nx,ny
    double precision, intent(in) :: hnom
    double precision, intent(out) :: jz(nx,ny)
    logical, intent(out) :: rope_mask(nx,ny)
    integer :: m,ix1,ix2
    double precision, allocatable :: xst(:,:,:,:),bst(:,:,:,:)
    double precision, allocatable :: xoffset(:,:),yoffset(:,:)
    double precision, allocatable :: s0(:,:),s1(:,:),s2(:,:),s3(:,:),s4(:,:)
    double precision, allocatable :: dbydx(:,:),dbxdy(:,:)
    double precision, allocatable :: hx(:,:),hy(:,:)
    integer, allocatable :: stx(:,:),sty(:,:)

    allocate(xst(1:nx,1:ny,1:1,1:ndim),bst(1:nx,1:ny,1:1,1:ndir))
    allocate(xoffset(nx,ny),yoffset(nx,ny))
    allocate(s0(nx,ny),s1(nx,ny),s2(nx,ny),s3(nx,ny),s4(nx,ny))
    allocate(dbydx(nx,ny),dbxdy(nx,ny),hx(nx,ny),hy(nx,ny))
    allocate(stx(nx,ny),sty(nx,ny))
    call build_rope_mask(nx,ny,rope_mask)
    call select_stencils(nx,ny,hnom,rope_mask,hx,hy,stx,sty)

    do m=0,4
      call fill_offset_grid(nx,ny,1,m,hx,stx,xoffset,yoffset)
      call eval_td99_grid(nx,ny,xoffset,yoffset,xst,bst)
      select case(m)
      case(0); s0=bst(:,:,1,2)
      case(1); s1=bst(:,:,1,2)
      case(2); s2=bst(:,:,1,2)
      case(3); s3=bst(:,:,1,2)
      case(4); s4=bst(:,:,1,2)
      end select
    end do
    do ix2=1,ny
      do ix1=1,nx
        call derivative_value(stx(ix1,ix2),hx(ix1,ix2),&
             s0(ix1,ix2),s1(ix1,ix2),s2(ix1,ix2),s3(ix1,ix2),&
             s4(ix1,ix2),dbydx(ix1,ix2))
      end do
    end do
    do m=0,4
      call fill_offset_grid(nx,ny,2,m,hy,sty,xoffset,yoffset)
      call eval_td99_grid(nx,ny,xoffset,yoffset,xst,bst)
      select case(m)
      case(0); s0=bst(:,:,1,1)
      case(1); s1=bst(:,:,1,1)
      case(2); s2=bst(:,:,1,1)
      case(3); s3=bst(:,:,1,1)
      case(4); s4=bst(:,:,1,1)
      end select
    end do
    do ix2=1,ny
      do ix1=1,nx
        call derivative_value(sty(ix1,ix2),hy(ix1,ix2),&
             s0(ix1,ix2),s1(ix1,ix2),s2(ix1,ix2),s3(ix1,ix2),&
             s4(ix1,ix2),dbxdy(ix1,ix2))
      end do
    end do
    jz=dbydx-dbxdy
    deallocate(xst,bst,xoffset,yoffset,s0,s1,s2,s3,s4,dbydx,dbxdy)
    deallocate(hx,hy,stx,sty)
  end subroutine compute_td99_jz

  subroutine derivative_value(stype,h,s0,s1,s2,s3,s4,value)
    integer, intent(in) :: stype
    double precision, intent(in) :: h,s0,s1,s2,s3,s4
    double precision, intent(out) :: value

    select case(stype)
    case(0)
      value=(-s4+8.d0*s3-8.d0*s1+s0)/(12.d0*h)
    case(1)
      value=(-25.d0*s0+48.d0*s1-36.d0*s2+16.d0*s3-3.d0*s4)/(12.d0*h)
    case(-1)
      value=(25.d0*s4-48.d0*s3+36.d0*s2-16.d0*s1+3.d0*s0)/(12.d0*h)
    end select
  end subroutine derivative_value

  subroutine build_rope_mask(nx,ny,rope_mask)
    use mod_global_parameters

    integer, intent(in) :: nx,ny
    logical, intent(out) :: rope_mask(nx,ny)
    integer :: ix1,ix2
    double precision :: xc,yc

    do ix2=1,ny
      yc=xprobmin2+(dble(ix2)-0.5d0)*dx(2,1)
      do ix1=1,nx
        xc=xprobmin1+(dble(ix1)-0.5d0)*dx(1,1)
        rope_mask(ix1,ix2)=td99_rho_z0(xc,yc)<a_TD99
      end do
    end do
  end subroutine build_rope_mask

  function td99_rho_z0(xc,yc) result(rho)
    double precision, intent(in) :: xc,yc
    double precision :: rho,rvertical

    rvertical=sqrt(yc*yc+d_TD99*d_TD99)
    rho=sqrt(xc*xc+(rvertical-R_TD99)**2)
  end function td99_rho_z0

  logical function stencil_inside(xc,yc,direction,stype,h)
    double precision, intent(in) :: xc,yc,h
    integer, intent(in) :: direction,stype
    integer :: m,mmin,mmax
    double precision :: xt,yt,offset

    select case(stype)
    case(0);  mmin=-2; mmax=2
    case(1);  mmin=0;  mmax=4
    case(-1); mmin=-4; mmax=0
    case default
      stencil_inside=.false.
      return
    end select
    stencil_inside=.true.
    do m=mmin,mmax
      offset=dble(m)*h
      xt=xc
      yt=yc
      if(direction==1) xt=xt+offset
      if(direction==2) yt=yt+offset
      if(td99_rho_z0(xt,yt)>=a_TD99) then
        stencil_inside=.false.
        return
      end if
    end do
  end function stencil_inside

  subroutine choose_stencil(xc,yc,hnom,direction,h,stype)
    double precision, intent(in) :: xc,yc,hnom
    integer, intent(in) :: direction
    double precision, intent(out) :: h
    integer, intent(out) :: stype
    integer :: attempt

    h=hnom
    stype=99
    do attempt=1,40
      if(stencil_inside(xc,yc,direction,0,h)) then
        stype=0
      else if(stencil_inside(xc,yc,direction,1,h)) then
        stype=1
      else if(stencil_inside(xc,yc,direction,-1,h)) then
        stype=-1
      end if
      if(stype/=99) return
      h=0.5d0*h
    end do
    stop 'adaptive TD99 stencil did not fit inside rho<a'
  end subroutine choose_stencil

  subroutine select_stencils(nx,ny,hnom,rope_mask,hx,hy,stx,sty)
    use mod_global_parameters

    integer, intent(in) :: nx,ny
    double precision, intent(in) :: hnom
    logical, intent(in) :: rope_mask(nx,ny)
    double precision, intent(out) :: hx(nx,ny),hy(nx,ny)
    integer, intent(out) :: stx(nx,ny),sty(nx,ny)
    integer :: ix1,ix2
    double precision :: xc,yc

    do ix2=1,ny
      yc=xprobmin2+(dble(ix2)-0.5d0)*dx(2,1)
      do ix1=1,nx
        xc=xprobmin1+(dble(ix1)-0.5d0)*dx(1,1)
        if(rope_mask(ix1,ix2)) then
          call choose_stencil(xc,yc,hnom,1,hx(ix1,ix2),stx(ix1,ix2))
          call choose_stencil(xc,yc,hnom,2,hy(ix1,ix2),sty(ix1,ix2))
        else
          hx(ix1,ix2)=hnom
          hy(ix1,ix2)=hnom
          stx(ix1,ix2)=0
          sty(ix1,ix2)=0
        end if
      end do
    end do
  end subroutine select_stencils

  subroutine fill_offset_grid(nx,ny,direction,m,h,stype,xoffset,yoffset)
    integer, intent(in) :: nx,ny,direction,m
    double precision, intent(in) :: h(nx,ny)
    integer, intent(in) :: stype(nx,ny)
    double precision, intent(out) :: xoffset(nx,ny),yoffset(nx,ny)
    integer :: ix1,ix2
    double precision :: offset

    xoffset=zero
    yoffset=zero
    do ix2=1,ny
      do ix1=1,nx
        select case(stype(ix1,ix2))
        case(0);  offset=(dble(m)-2.d0)*h(ix1,ix2)
        case(1);  offset=dble(m)*h(ix1,ix2)
        case(-1); offset=(dble(m)-4.d0)*h(ix1,ix2)
        end select
        if(direction==1) xoffset(ix1,ix2)=offset
        if(direction==2) yoffset(ix1,ix2)=offset
      end do
    end do
  end subroutine fill_offset_grid

  subroutine eval_td99_grid(nx,ny,xoffset,yoffset,xst,bst)
    use mod_global_parameters

    integer, intent(in) :: nx,ny
    double precision, intent(in) :: xoffset(nx,ny),yoffset(nx,ny)
    double precision, intent(out) :: xst(1:nx,1:ny,1:1,1:ndim)
    double precision, intent(out) :: bst(1:nx,1:ny,1:1,1:ndir)
    integer :: ix1,ix2

    do ix2=1,ny
      do ix1=1,nx
        xst(ix1,ix2,1,1)=xprobmin1+(dble(ix1)-0.5d0)*dx(1,1)+xoffset(ix1,ix2)
        xst(ix1,ix2,1,2)=xprobmin2+(dble(ix2)-0.5d0)*dx(2,1)+yoffset(ix1,ix2)
        xst(ix1,ix2,1,3)=0.d0
      end do
    end do
    call TD99(1,1,1,nx,ny,1,1,1,1,nx,ny,1,xst,bst)
  end subroutine eval_td99_grid

end module mod_usr
