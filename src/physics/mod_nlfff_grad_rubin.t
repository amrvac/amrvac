!> MPI-parallel fixed-grid Grad--Rubin NLFFF extrapolation.
!>
!> This is an independent implementation of the current-field iteration in
!> Wheatland (2006, 2007).  The dense-grid representation is deliberate: it
!> gives every rank deterministic read-only access during field-line tracing.
module mod_nlfff_grad_rubin
  implicit none
  private

  integer, parameter, public :: gr_trace_bottom=1
  integer, parameter, public :: gr_trace_open=2
  integer, parameter, public :: gr_trace_weak=3
  integer, parameter, public :: gr_trace_max_steps=4

  type, public :: nlfff_grad_rubin_config
    integer :: fft_padding_factor=2
    character(len=16) :: flux_treatment='strict'
    double precision :: max_flux_imbalance=0.1d0
    integer :: polarity=1
    character(len=24) :: alpha_source='vector_magnetogram'
    double precision :: bz_taper_zero=0.01d0
    double precision :: bz_taper_full=0.02d0
    double precision :: relaxation_factor=0.5d0
    double precision :: fieldline_step_fraction=0.5d0
    integer :: fieldline_max_steps=10000
    integer :: max_iterations=50
    integer :: convergence_streak=3
    double precision :: field_change_tolerance=1.d-5
    double precision :: energy_change_tolerance=1.d-6
    integer :: self_consistency_cycles=0
    integer :: self_consistency_streak=1
    double precision :: self_consistency_field_tolerance=1.d-3
    double precision :: self_consistency_alpha_tolerance=1.d-3
    integer :: log_interval=1
    double precision :: memory_limit_mb=2048.d0
    logical :: write_detailed_history=.false.
  end type nlfff_grad_rubin_config

  type, public :: nlfff_grad_rubin_result
    integer :: polarity=0
    integer :: iterations=0
    logical :: converged=.false.
    character(len=32) :: stop_reason='not_started'
    double precision :: initial_energy=0.d0
    double precision :: final_energy=0.d0
    double precision :: rms_field_change=0.d0
    double precision :: relative_energy_change=0.d0
    double precision :: current_weighted_theta=0.d0
    double precision :: epsilon_div=0.d0
    integer :: closed_fieldlines=0
    integer :: open_fieldlines=0
    integer :: weak_fieldlines=0
    integer :: max_step_fieldlines=0
    integer :: valid_alpha_pixels=0
    double precision :: alpha_min=0.d0
    double precision :: alpha_max=0.d0
    double precision :: memory_required_mb=0.d0
    integer :: self_consistency_cycles_completed=0
    double precision :: positive_negative_rms_difference=0.d0
    double precision :: boundary_alpha_rms_change=0.d0
  end type nlfff_grad_rubin_result

{^IFTHREED
  double precision, allocatable, save :: boundary_b(:,:,:)
  double precision, allocatable, save :: external_alpha_raw_full(:,:)
  double precision, allocatable, save :: external_alpha_full(:,:)
  double precision, allocatable, save :: external_alpha_weight_full(:,:)
  integer, allocatable, save :: external_alpha_pil_full(:,:)
  integer, allocatable, save :: external_alpha_valid_full(:,:)
  integer, allocatable, save :: external_alpha_polarity_full(:,:)
  double precision, allocatable, save :: external_alpha_x(:),external_alpha_y(:)
  character(len=1024), save :: external_alpha_filename=''
  double precision, save :: external_alpha_unit_length_cm=0.d0
  double precision, save :: external_alpha_unit_magneticfield_g=0.d0

  public :: init_nlfff_grad_rubin_boundary
  public :: init_nlfff_grad_rubin_external_alpha
  public :: extrapolate_nlfff_grad_rubin
  public :: gr_taper_weight
  public :: gr_compute_alpha_boundary
  public :: gr_trilinear_sample
  public :: gr_bilinear_sample
  public :: gr_combine_alpha_unweighted
  public :: gr_solve_current_field_mpi
}

contains

{^IFTHREED
  subroutine init_nlfff_grad_rubin_boundary(filename,unit_length,&
     unit_magneticfield,qxc1,qxc2,alpha_filename)
    use mod_lfff, only: init_b_fff_data_driven_boundary

    character(len=*), intent(in) :: filename
    double precision, intent(in) :: unit_length,unit_magneticfield
    double precision, intent(in), optional :: qxc1,qxc2
    character(len=*), intent(in), optional :: alpha_filename

    if(allocated(boundary_b)) deallocate(boundary_b)
    call init_b_fff_data_driven_boundary(filename,unit_length,&
       unit_magneticfield,qxc1,qxc2,boundary_b)
    if(present(alpha_filename)) then
      if(len_trim(alpha_filename)>0) call init_nlfff_grad_rubin_external_alpha(&
         trim(alpha_filename),unit_length,unit_magneticfield)
    end if
  end subroutine init_nlfff_grad_rubin_boundary

  subroutine init_nlfff_grad_rubin_external_alpha(filename,unit_length,&
     unit_magneticfield)
    use mod_comm_lib, only: mpistop
    character(len=*), intent(in) :: filename
    double precision, intent(in) :: unit_length,unit_magneticfield

    if(.not.allocated(boundary_b)) call mpistop('external alpha requires an initialized vector boundary')
    call read_external_alpha_product(trim(filename),unit_length,unit_magneticfield)
  end subroutine init_nlfff_grad_rubin_external_alpha

  subroutine read_external_alpha_product(filename,unit_length,unit_magneticfield)
    use mpi
    use mod_global_parameters, only: icomm,mype
    use mod_lfff, only: xa1,xa2,nx1,nx2
    use mod_comm_lib, only: mpistop
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    character(len=*), intent(in) :: filename
    double precision, intent(in) :: unit_length,unit_magneticfield
    character(len=32) :: magic
    character(len=4096) :: metadata_blob
    integer :: iu,ios,i,j,version,nx,ny,reserved,alpha_code,mask_code
    integer :: metadata_length,ierrmpi
    double precision :: product_length,product_bfield,dx_product,dy_product
    double precision :: xc_product,yc_product,tol
    logical :: exists
    character(len=32), parameter :: expected_magic='AMRVAC_EXTERNAL_ALPHA_V1'

    if(len_trim(filename)>len(external_alpha_filename)) &
       call mpistop('external alpha filename is too long')
    if(mype==0) then
      inquire(file=trim(filename),exist=exists)
      if(.not.exists) call mpistop('missing external alpha product')
      open(newunit=iu,file=trim(filename),status='old',access='stream',&
         form='unformatted',action='read',iostat=ios)
      if(ios/=0) call mpistop('cannot open external alpha product')
      read(iu,iostat=ios) magic,version,nx,ny,reserved,product_length,&
         product_bfield,dx_product,dy_product,xc_product,yc_product,&
         alpha_code,mask_code
      if(ios/=0) call mpistop('cannot read external alpha header')
      if(magic/=expected_magic .or. version/=1 .or. nx<1 .or. ny<1) &
         call mpistop('invalid external alpha magic/version/grid')
      if(alpha_code/=1 .or. mask_code/=1) &
         call mpistop('unsupported external alpha unit/mask code')
      allocate(external_alpha_x(nx),external_alpha_y(ny))
      allocate(external_alpha_raw_full(nx,ny),external_alpha_full(nx,ny))
      allocate(external_alpha_weight_full(nx,ny))
      allocate(external_alpha_pil_full(nx,ny),external_alpha_valid_full(nx,ny))
      allocate(external_alpha_polarity_full(nx,ny))
      read(iu,iostat=ios) external_alpha_x,external_alpha_y
      read(iu,iostat=ios) external_alpha_raw_full,external_alpha_full,&
         external_alpha_weight_full,external_alpha_pil_full,&
         external_alpha_valid_full,external_alpha_polarity_full
      read(iu,iostat=ios) metadata_length,metadata_blob
      close(iu)
      if(ios/=0) call mpistop('cannot read external alpha payload')
      if(metadata_length<0 .or. metadata_length>len(metadata_blob)) &
         call mpistop('invalid external alpha metadata length')
    else
      nx=0; ny=0; product_length=0.d0; product_bfield=0.d0
      dx_product=0.d0; dy_product=0.d0; xc_product=0.d0; yc_product=0.d0
      alpha_code=0; mask_code=0
    end if
    call MPI_BCAST(nx,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(ny,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(product_length,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(product_bfield,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(dx_product,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(dy_product,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(xc_product,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(yc_product,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    if(mype/=0) then
      allocate(external_alpha_x(nx),external_alpha_y(ny))
      allocate(external_alpha_raw_full(nx,ny),external_alpha_full(nx,ny))
      allocate(external_alpha_weight_full(nx,ny))
      allocate(external_alpha_pil_full(nx,ny),external_alpha_valid_full(nx,ny))
      allocate(external_alpha_polarity_full(nx,ny))
    end if
    call MPI_BCAST(external_alpha_x,nx,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(external_alpha_y,ny,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(external_alpha_raw_full,nx*ny,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(external_alpha_full,nx*ny,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(external_alpha_weight_full,nx*ny,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(external_alpha_pil_full,nx*ny,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(external_alpha_valid_full,nx*ny,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(external_alpha_polarity_full,nx*ny,MPI_INTEGER,0,icomm,ierrmpi)

    if(.not.ieee_is_finite(product_length) .or. .not.ieee_is_finite(product_bfield) .or. &
       product_length<=0.d0 .or. product_bfield<=0.d0) &
       call mpistop('external alpha contains invalid unit metadata')
    tol=1.d-12*max(1.d0,dabs(unit_length),dabs(product_length))
    if(dabs(product_length-unit_length)>tol) &
       call mpistop('external alpha unit_length does not match AMRVAC')
    tol=1.d-12*max(1.d0,dabs(unit_magneticfield),dabs(product_bfield))
    if(dabs(product_bfield-unit_magneticfield)>tol) &
       call mpistop('external alpha magnetic unit does not match AMRVAC')
    if(nx/=nx1 .or. ny/=nx2) then
      if(mype==0) write(*,*) 'external alpha grid/product:',nx,ny,' magnetogram:',nx1,nx2
      call mpistop('external alpha grid shape does not match magnetogram')
    end if
    tol=1.d-12*max(1.d0,maxval(dabs(external_alpha_x)),maxval(dabs(xa1)))
    if(maxval(dabs(external_alpha_x-xa1))>tol) &
       call mpistop('external alpha x coordinates do not match magnetogram')
    tol=1.d-12*max(1.d0,maxval(dabs(external_alpha_y)),maxval(dabs(xa2)))
    if(maxval(dabs(external_alpha_y-xa2))>tol) &
       call mpistop('external alpha y coordinates do not match magnetogram')
    if(dx_product<=0.d0 .or. dy_product<=0.d0) &
       call mpistop('external alpha spacing must be positive')
    if(.not.all(ieee_is_finite(external_alpha_raw_full)) .or. &
       .not.all(ieee_is_finite(external_alpha_full)) .or. &
       .not.all(ieee_is_finite(external_alpha_weight_full))) &
       call mpistop('external alpha contains non-finite arrays')
    if(any(external_alpha_weight_full<0)) &
       call mpistop('external alpha weight contains negative values')
    do j=1,ny
      do i=1,nx
        if(external_alpha_pil_full(i,j)/=-1 .and. external_alpha_pil_full(i,j)/=0 .and. &
           external_alpha_pil_full(i,j)/=1) call mpistop('invalid external alpha PIL mask')
        if(external_alpha_valid_full(i,j)/=0 .and. external_alpha_valid_full(i,j)/=1) &
           call mpistop('invalid external alpha valid mask')
        if(external_alpha_polarity_full(i,j)/=-1 .and. &
           external_alpha_polarity_full(i,j)/=0 .and. external_alpha_polarity_full(i,j)/=1) &
           call mpistop('invalid external alpha polarity mask')
        if(external_alpha_valid_full(i,j)/=0 .and. &
           external_alpha_polarity_full(i,j)/=0) then
          if(external_alpha_polarity_full(i,j)*boundary_b(i,j,3)<=0.d0) &
             call mpistop('external alpha polarity disagrees with Bz')
        end if
        if(external_alpha_valid_full(i,j)/=0 .and. &
           external_alpha_polarity_full(i,j)==0) &
           call mpistop('external alpha valid support requires nonzero polarity')
      end do
    end do
    external_alpha_filename=trim(filename)
    external_alpha_unit_length_cm=product_length
    external_alpha_unit_magneticfield_g=product_bfield
    if(mype==0) then
      write(*,*) 'external alpha product:',trim(filename)
      write(*,*) 'external alpha support:',count(external_alpha_valid_full/=0)
    end if
  end subroutine read_external_alpha_product

  subroutine extrapolate_nlfff_grad_rubin(iw_b,config,result)
    use mpi
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    use mod_lfff, only: extrapolate_potential_fft_dense
    use mod_nlfff_diagnostics, only: nlfff_physical_metrics,&
       evaluate_nlfff_metrics_dense,write_nlfff_metrics_header,&
       write_nlfff_metrics_row
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    integer, intent(in) :: iw_b(3)
    type(nlfff_grad_rubin_config), intent(in) :: config
    type(nlfff_grad_rubin_result), intent(out) :: result

    double precision, allocatable :: bcore(:,:,:),alpha0(:,:),alpha_raw(:,:),alpha_weight(:,:)
    double precision, allocatable :: potential(:,:,:,:)
    double precision, allocatable :: b(:,:,:,:),bnew(:,:,:,:),alpha(:,:,:)
    double precision, allocatable :: current(:,:,:,:),bc(:,:,:,:)
    logical, allocatable :: alpha_mask(:,:)
    integer, allocatable :: alpha_polarity(:,:)
    double precision :: dx1,dx2,dx3,energy_old,energy_new,rms_change,energy_change
    double precision :: memory_mb,flux_before,flux_after,mean_correction
    double precision :: theta,epsdiv
    type(nlfff_physical_metrics) :: physical_metrics
    integer :: iter,stable,log_unit,metrics_unit,flux_status,ic
    integer :: counts(4),valid_pixels
    character(len=16) :: polarity_name
    logical :: log_open

    result=nlfff_grad_rubin_result()
    result%polarity=config%polarity
    call validate_configuration(iw_b,config)
    if(config%self_consistency_cycles>0) then
      call extrapolate_self_consistent(iw_b,config,result)
      return
    end if
    memory_mb=estimate_memory_mb(config%fft_padding_factor)
    result%memory_required_mb=memory_mb
    if(memory_mb>config%memory_limit_mb) then
      if(mype==0) write(*,*) 'Grad-Rubin required/allowed MiB per rank:',&
         memory_mb,config%memory_limit_mb
      call mpistop('Grad-Rubin dense solver exceeds memory_limit_mb')
    end if
    call extract_boundary_core(bcore)
    dx1=dx(1,1)
    dx2=dx(2,1)
    dx3=dx(3,1)
    call balance_normal_field(bcore(:,:,3),config%flux_treatment,&
       config%max_flux_imbalance,flux_before,flux_after,mean_correction,&
       flux_status)
    if(flux_status/=0 .and. flux_status/=1) then
      if(mype==0) write(*,*) 'Grad-Rubin flux imbalance:',flux_before
      call mpistop('Grad-Rubin bottom normal field failed flux-balance policy')
    end if

    allocate(alpha0(domain_nx1,domain_nx2),alpha_raw(domain_nx1,domain_nx2))
    allocate(alpha_weight(domain_nx1,domain_nx2),alpha_polarity(domain_nx1,domain_nx2))
    allocate(alpha_mask(domain_nx1,domain_nx2))
    if(trim(adjustl(config%alpha_source))=='external') then
      call extract_external_alpha_core(alpha_raw,alpha0,alpha_weight,&
         alpha_mask,alpha_polarity)
      valid_pixels=count(alpha_mask)
    else
      call gr_compute_alpha_boundary(bcore,dx1,dx2,config%bz_taper_zero,&
         config%bz_taper_full,alpha0,valid_pixels,alpha_raw)
      alpha_weight=0.d0
      where(alpha0/=0.d0) alpha_weight=1.d0
      alpha_polarity=0
      where(bcore(:,:,3)>0.d0) alpha_polarity=1
      where(bcore(:,:,3)<0.d0) alpha_polarity=-1
      alpha_mask=config%polarity*bcore(:,:,3)>0.d0 .and. &
         dabs(bcore(:,:,3))/max(maxval(dabs(bcore(:,:,3))),tiny(1.d0))>&
         config%bz_taper_zero
    end if
    if(trim(adjustl(config%alpha_source))=='external') then
      alpha_mask=alpha_mask .and. config%polarity*bcore(:,:,3)>0.d0
    end if
    where(.not.alpha_mask) alpha0=0.d0
    result%valid_alpha_pixels=count(alpha_mask)
    if(result%valid_alpha_pixels>0) then
      result%alpha_min=minval(alpha0,mask=alpha_mask)
      result%alpha_max=maxval(alpha0,mask=alpha_mask)
    end if
    if(config%write_detailed_history .and. mype==0) then
      call write_alpha_diagnostics(trim(base_filename),config%polarity,&
         trim(adjustl(config%alpha_source)),alpha_raw,alpha0,alpha_weight,alpha_mask,&
         alpha_polarity)
    end if

    call extrapolate_potential_fft_dense(bcore(:,:,3),dx1,dx2,dx3,&
       config%fft_padding_factor,potential)
    allocate(b,source=potential)
    allocate(bnew,source=potential)
    allocate(alpha(domain_nx1,domain_nx2,domain_nx3+1))
    allocate(current(domain_nx1,domain_nx2,domain_nx3+1,3))
    allocate(bc(domain_nx1,domain_nx2,domain_nx3+1,3))

    energy_old=gr_energy(b,dx1,dx2,dx3)
    result%initial_energy=energy_old
    stable=0
    log_open=.false.
    log_unit=-1
    metrics_unit=-1
    if(mype==0) then
      if(config%polarity>0) then
        polarity_name='positive'
      else
        polarity_name='negative'
      end if
      if(config%write_detailed_history) then
        open(newunit=log_unit,file=trim(base_filename)//'_grad_rubin_'//&
           trim(polarity_name)//'.csv',status='replace',action='write')
        log_open=.true.
        write(log_unit,'(a)') 'iteration,energy,rms_field_change,'//&
           'relative_energy_change,theta_j_deg,epsilon_div,closed,open,weak,max_steps'
      end if
      open(newunit=metrics_unit,file=trim(base_filename)//'_nlfff_metrics.csv',&
         status='replace',action='write')
      call write_nlfff_metrics_header(metrics_unit)
    end if
    call evaluate_nlfff_metrics_dense(b,dx1,dx2,dx3,physical_metrics)
    if(mype==0) call write_nlfff_metrics_row(metrics_unit,0,physical_metrics)

    do iter=1,config%max_iterations
      call transport_alpha_mpi(b,alpha0,bcore(:,:,3),config,alpha,counts)
      do concurrent(ic=1:3)
        current(:,:,:,ic)=alpha*b(:,:,:,ic)
      end do
      call gr_solve_current_field_mpi(current,config%fft_padding_factor,dx1,dx2,&
         dx3,bc)
      bnew=(1.d0-config%relaxation_factor)*b+&
         config%relaxation_factor*(potential+bc)
      bnew(:,:,1,3)=bcore(:,:,3)
      if(.not.all(ieee_is_finite(bnew))) &
         call mpistop('non-finite field in Grad-Rubin iteration')

      energy_new=gr_energy(bnew,dx1,dx2,dx3)
      rms_change=dsqrt(sum((bnew-b)**2)/max(sum(b**2),tiny(1.d0)))
      energy_change=dabs(energy_new-energy_old)/max(dabs(energy_old),tiny(1.d0))
      call field_diagnostics(bnew,dx1,dx2,dx3,theta,epsdiv)
      b=bnew
      energy_old=energy_new
      result%iterations=iter
      result%final_energy=energy_new
      result%rms_field_change=rms_change
      result%relative_energy_change=energy_change
      result%current_weighted_theta=theta
      result%epsilon_div=epsdiv
      result%closed_fieldlines=counts(1)
      result%open_fieldlines=counts(2)
      result%weak_fieldlines=counts(3)
      result%max_step_fieldlines=counts(4)

      if(mod(iter,config%log_interval)==0 .or. iter==1) then
        call evaluate_nlfff_metrics_dense(b,dx1,dx2,dx3,physical_metrics)
        if(mype==0) then
          if(config%write_detailed_history) then
            write(log_unit,'(i0,5(a,es24.16),4(a,i0))') iter,',',energy_new,',',&
               rms_change,',',energy_change,',',theta,',',epsdiv,',',counts(1),&
               ',',counts(2),',',counts(3),',',counts(4)
            flush(log_unit)
          end if
          call write_nlfff_metrics_row(metrics_unit,iter,physical_metrics)
        end if
      end if

      if(rms_change<=config%field_change_tolerance .and. &
         energy_change<=config%energy_change_tolerance) then
        stable=stable+1
      else
        stable=0
      end if
      if(stable>=config%convergence_streak) then
        result%converged=.true.
        result%stop_reason='converged'
        exit
      end if
    end do
    if(.not.result%converged) result%stop_reason='max_iterations'
    if(mype==0) then
      if(log_open) close(log_unit)
      close(metrics_unit)
    end if
    call scatter_dense_field(b,iw_b)

    deallocate(bcore,alpha0,alpha_raw,alpha_weight,alpha_polarity,alpha_mask,&
       potential,b,bnew,alpha,current,bc)
  end subroutine extrapolate_nlfff_grad_rubin

  subroutine extrapolate_self_consistent(iw_b,config,result)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    use mod_lfff, only: extrapolate_potential_fft_dense
    use mod_nlfff_diagnostics, only: nlfff_physical_metrics,&
       evaluate_nlfff_metrics_dense,write_nlfff_metrics_header,&
       write_nlfff_metrics_row
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    integer, intent(in) :: iw_b(3)
    type(nlfff_grad_rubin_config), intent(in) :: config
    type(nlfff_grad_rubin_result), intent(out) :: result

    type(nlfff_grad_rubin_config) :: polarity_config
    type(nlfff_grad_rubin_result) :: positive_result,negative_result
    double precision, allocatable :: bcore(:,:,:),alpha0(:,:),alpha_old(:,:)
    double precision, allocatable :: alpha_raw(:,:),alpha_weight(:,:)
    double precision, allocatable :: alpha_positive(:,:,:),alpha_negative(:,:,:)
    double precision, allocatable :: potential(:,:,:,:),bpositive(:,:,:,:)
    double precision, allocatable :: bnegative(:,:,:,:),bfinal(:,:,:,:)
    logical, allocatable :: valid_mask(:,:),alpha_mask(:,:)
    integer, allocatable :: alpha_polarity(:,:)
    double precision :: dx1,dx2,dx3,memory_mb,flux_before,flux_after
    double precision :: mean_correction,pn_difference,alpha_change,denom
    type(nlfff_physical_metrics) :: physical_metrics
    integer :: cycle,stable,flux_status,valid_pixels,log_unit,metrics_unit
    integer :: total_iterations
    integer :: counts_positive(4),counts_negative(4)
    logical :: log_open

    result=nlfff_grad_rubin_result()
    result%polarity=0
    memory_mb=estimate_self_consistent_memory_mb(config%fft_padding_factor)
    result%memory_required_mb=memory_mb
    if(memory_mb>config%memory_limit_mb) then
      if(mype==0) write(*,*) 'Grad-Rubin self-consistent required/allowed MiB per rank:',&
         memory_mb,config%memory_limit_mb
      call mpistop('Grad-Rubin self-consistency exceeds memory_limit_mb')
    end if

    call extract_boundary_core(bcore)
    dx1=dx(1,1); dx2=dx(2,1); dx3=dx(3,1)
    call balance_normal_field(bcore(:,:,3),config%flux_treatment,&
       config%max_flux_imbalance,flux_before,flux_after,mean_correction,&
       flux_status)
    if(flux_status/=0 .and. flux_status/=1) &
       call mpistop('Grad-Rubin bottom normal field failed flux-balance policy')

    allocate(alpha0(domain_nx1,domain_nx2),alpha_old(domain_nx1,domain_nx2))
    allocate(alpha_raw(domain_nx1,domain_nx2),alpha_weight(domain_nx1,domain_nx2))
    allocate(alpha_polarity(domain_nx1,domain_nx2))
    allocate(valid_mask(domain_nx1,domain_nx2),alpha_mask(domain_nx1,domain_nx2))
    if(trim(adjustl(config%alpha_source))=='external') then
      call extract_external_alpha_core(alpha_raw,alpha0,alpha_weight,&
         alpha_mask,alpha_polarity)
      valid_mask=alpha_mask
    else
      call gr_compute_alpha_boundary(bcore,dx1,dx2,config%bz_taper_zero,&
         config%bz_taper_full,alpha0,valid_pixels,alpha_raw)
      valid_mask=dabs(bcore(:,:,3))/max(maxval(dabs(bcore(:,:,3))),tiny(1.d0))>&
         config%bz_taper_zero
      alpha_weight=0.d0
      where(alpha0/=0.d0) alpha_weight=1.d0
      alpha_polarity=0
      where(bcore(:,:,3)>0.d0) alpha_polarity=1
      where(bcore(:,:,3)<0.d0) alpha_polarity=-1
    end if
    where(.not.valid_mask) alpha0=0.d0
    result%valid_alpha_pixels=count(valid_mask)
    if(result%valid_alpha_pixels>0) then
      result%alpha_min=minval(alpha0,mask=valid_mask)
      result%alpha_max=maxval(alpha0,mask=valid_mask)
    end if

    call extrapolate_potential_fft_dense(bcore(:,:,3),dx1,dx2,dx3,&
       config%fft_padding_factor,potential)
    allocate(bpositive,source=potential)
    allocate(bnegative,source=potential)
    allocate(bfinal,source=potential)
    allocate(alpha_positive(domain_nx1,domain_nx2,domain_nx3+1))
    allocate(alpha_negative(domain_nx1,domain_nx2,domain_nx3+1))
    result%initial_energy=gr_energy(potential,dx1,dx2,dx3)
    stable=0
    total_iterations=0
    log_open=.false.; log_unit=-1; metrics_unit=-1
    if(mype==0) then
      if(config%write_detailed_history) then
        open(newunit=log_unit,file=trim(base_filename)//&
           '_grad_rubin_self_consistent.csv',status='replace',action='write')
        log_open=.true.
        write(log_unit,'(a)') 'cycle,positive_iterations,negative_iterations,'//&
           'positive_converged,negative_converged,'//&
           'positive_negative_rms_difference,boundary_alpha_rms_change,'//&
           'positive_energy,negative_energy,positive_closed,negative_closed,'//&
           'positive_open,negative_open,positive_weak,negative_weak,'//&
           'positive_max_steps,negative_max_steps'
      end if
      open(newunit=metrics_unit,file=trim(base_filename)//'_nlfff_metrics.csv',&
         status='replace',action='write')
      call write_nlfff_metrics_header(metrics_unit)
    end if
    call evaluate_nlfff_metrics_dense(potential,dx1,dx2,dx3,physical_metrics)
    if(mype==0) call write_nlfff_metrics_row(metrics_unit,0,physical_metrics)

    polarity_config=config
    polarity_config%self_consistency_cycles=0
    do cycle=1,config%self_consistency_cycles
      alpha_old=alpha0
      polarity_config%polarity=1
      call solve_polarity(potential,bcore(:,:,3),alpha0,polarity_config,&
         bpositive,alpha_positive,positive_result,counts_positive)
      polarity_config%polarity=-1
      call solve_polarity(potential,bcore(:,:,3),alpha0,polarity_config,&
         bnegative,alpha_negative,negative_result,counts_negative)
      total_iterations=total_iterations+positive_result%iterations+&
         negative_result%iterations

      call gr_combine_alpha_unweighted(alpha_positive(:,:,1),&
         alpha_negative(:,:,1),alpha0)
      where(.not.valid_mask) alpha0=0.d0
      denom=max(sum(alpha_old**2),tiny(1.d0))
      alpha_change=dsqrt(sum((alpha0-alpha_old)**2)/denom)
      pn_difference=dsqrt(sum((bpositive-bnegative)**2)/&
         max(0.5d0*(sum(bpositive**2)+sum(bnegative**2)),tiny(1.d0)))
      if(.not.ieee_is_finite(alpha_change) .or. &
         .not.ieee_is_finite(pn_difference)) &
         call mpistop('non-finite Grad-Rubin self-consistency diagnostic')

      result%self_consistency_cycles_completed=cycle
      result%positive_negative_rms_difference=pn_difference
      result%boundary_alpha_rms_change=alpha_change
      result%iterations=total_iterations
      result%final_energy=negative_result%final_energy
      result%rms_field_change=negative_result%rms_field_change
      result%relative_energy_change=negative_result%relative_energy_change
      result%closed_fieldlines=counts_positive(1)+counts_negative(1)
      result%open_fieldlines=counts_positive(2)+counts_negative(2)
      result%weak_fieldlines=counts_positive(3)+counts_negative(3)
      result%max_step_fieldlines=counts_positive(4)+counts_negative(4)
      call evaluate_nlfff_metrics_dense(bnegative,dx1,dx2,dx3,physical_metrics)

      if(mype==0) then
        if(config%write_detailed_history) then
          write(log_unit,'(i0,2(",",i0),2(",",l1),4(",",es24.16),'//&
             '8(",",i0))') cycle,positive_result%iterations,&
             negative_result%iterations,positive_result%converged,&
             negative_result%converged,pn_difference,alpha_change,&
             positive_result%final_energy,negative_result%final_energy,&
             counts_positive(1),counts_negative(1),counts_positive(2),&
             counts_negative(2),counts_positive(3),counts_negative(3),&
             counts_positive(4),counts_negative(4)
          flush(log_unit)
        end if
        call write_nlfff_metrics_row(metrics_unit,cycle,physical_metrics)
      end if

      if(pn_difference<=config%self_consistency_field_tolerance .and. &
         alpha_change<=config%self_consistency_alpha_tolerance) then
        stable=stable+1
      else
        stable=0
      end if
      if(stable>=config%self_consistency_streak) then
        result%converged=.true.
        result%stop_reason='self_consistent'
        exit
      end if
    end do
    if(.not.result%converged) result%stop_reason='self_consistency_cycles'
    if(mype==0) then
      if(log_open) close(log_unit)
      close(metrics_unit)
    end if
    if(result%valid_alpha_pixels>0) then
      result%alpha_min=minval(alpha0,mask=valid_mask)
      result%alpha_max=maxval(alpha0,mask=valid_mask)
    end if

    ! CFIT alternates complete P and N solutions and leaves the last polarity
    ! field as its output.  Use the final N solution because every AMRVAC cycle
    ! above is a complete P/N pair; the P/N discrepancy is returned explicitly.
    bfinal=bnegative
    call field_diagnostics(bfinal,dx1,dx2,dx3,result%current_weighted_theta,&
       result%epsilon_div)
    call scatter_dense_field(bfinal,iw_b)

    deallocate(bcore,alpha0,alpha_old,alpha_raw,alpha_weight,alpha_polarity,&
       valid_mask,alpha_mask,alpha_positive,alpha_negative,potential,bpositive,&
       bnegative,bfinal)
  end subroutine extrapolate_self_consistent

  subroutine solve_polarity(potential,bz0,alpha0,config,b,alpha,result,counts)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters, only: dx
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    double precision, intent(in) :: potential(:,:,:,:),bz0(:,:),alpha0(:,:)
    type(nlfff_grad_rubin_config), intent(in) :: config
    double precision, intent(out) :: b(:,:,:,:),alpha(:,:,:)
    type(nlfff_grad_rubin_result), intent(out) :: result
    integer, intent(out) :: counts(4)
    double precision, allocatable :: bnew(:,:,:,:),current(:,:,:,:),bc(:,:,:,:)
    double precision :: energy_old,energy_new,rms_change,energy_change,theta,epsdiv
    integer :: iter,stable,ic

    result=nlfff_grad_rubin_result()
    result%polarity=config%polarity
    allocate(bnew,source=potential)
    allocate(current(size(potential,1),size(potential,2),size(potential,3),3))
    allocate(bc(size(potential,1),size(potential,2),size(potential,3),3))
    b=potential
    energy_old=gr_energy(b,dx(1,1),dx(2,1),dx(3,1))
    result%initial_energy=energy_old
    stable=0; counts=0
    do iter=1,config%max_iterations
      call transport_alpha_mpi(b,alpha0,bz0,config,alpha,counts)
      do concurrent(ic=1:3)
        current(:,:,:,ic)=alpha*b(:,:,:,ic)
      end do
      call gr_solve_current_field_mpi(current,config%fft_padding_factor,&
         dx(1,1),dx(2,1),dx(3,1),bc)
      bnew=(1.d0-config%relaxation_factor)*b+&
         config%relaxation_factor*(potential+bc)
      bnew(:,:,1,3)=bz0
      if(.not.all(ieee_is_finite(bnew))) &
         call mpistop('non-finite field in Grad-Rubin polarity solve')
      energy_new=gr_energy(bnew,dx(1,1),dx(2,1),dx(3,1))
      rms_change=dsqrt(sum((bnew-b)**2)/max(sum(b**2),tiny(1.d0)))
      energy_change=dabs(energy_new-energy_old)/max(dabs(energy_old),tiny(1.d0))
      b=bnew; energy_old=energy_new
      result%iterations=iter
      result%final_energy=energy_new
      result%rms_field_change=rms_change
      result%relative_energy_change=energy_change
      if(rms_change<=config%field_change_tolerance .and. &
         energy_change<=config%energy_change_tolerance) then
        stable=stable+1
      else
        stable=0
      end if
      if(stable>=config%convergence_streak) then
        result%converged=.true.
        result%stop_reason='converged'
        exit
      end if
    end do
    if(.not.result%converged) result%stop_reason='max_iterations'
    call field_diagnostics(b,dx(1,1),dx(2,1),dx(3,1),theta,epsdiv)
    result%current_weighted_theta=theta; result%epsilon_div=epsdiv
    result%closed_fieldlines=counts(1); result%open_fieldlines=counts(2)
    result%weak_fieldlines=counts(3); result%max_step_fieldlines=counts(4)
    deallocate(bnew,current,bc)
  end subroutine solve_polarity

  pure subroutine gr_combine_alpha_unweighted(alpha_positive,alpha_negative,&
     alpha_combined)
    double precision, intent(in) :: alpha_positive(:,:),alpha_negative(:,:)
    double precision, intent(out) :: alpha_combined(size(alpha_positive,1),&
       size(alpha_positive,2))
    alpha_combined=0.5d0*(alpha_positive+alpha_negative)
  end subroutine gr_combine_alpha_unweighted

  subroutine validate_configuration(iw_b,config)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    use mod_geometry, only: coordinate,Cartesian
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    integer, intent(in) :: iw_b(3)
    type(nlfff_grad_rubin_config), intent(in) :: config

    if(.not.allocated(boundary_b)) call mpistop('Grad-Rubin boundary is not initialized')
    if(.not.all(ieee_is_finite(boundary_b))) &
       call mpistop('Grad-Rubin boundary contains non-finite values')
    if(ndim/=3 .or. coordinate/=Cartesian) &
       call mpistop('Grad-Rubin v1 requires Cartesian 3D')
    if(any(stretched_dim) .or. refine_max_level/=1 .or. levmax/=1) &
       call mpistop('Grad-Rubin v1 requires a single uniform level')
    if(stagger_grid) call mpistop('Grad-Rubin v1 does not support stagger_grid')
    if(any(iw_b<1) .or. any(iw_b>nw)) call mpistop('invalid Grad-Rubin B indices')
    if(config%polarity/=1 .and. config%polarity/=-1) &
       call mpistop('Grad-Rubin polarity must be +1 or -1')
    select case(trim(adjustl(config%alpha_source)))
    case('vector_magnetogram')
      continue
    case('external')
      if(.not.allocated(external_alpha_full)) &
         call mpistop('alpha_source=external requires an external alpha product')
    case default
      call mpistop('Grad-Rubin alpha_source must be vector_magnetogram or external')
    end select
    if(config%fft_padding_factor<1) call mpistop('Grad-Rubin padding must be positive')
    if(domain_nx1<3 .or. domain_nx2<3 .or. domain_nx3<2) &
       call mpistop('Grad-Rubin grid must be at least 3 x 3 x 2')
    if(mod(domain_nx1,block_nx1)/=0 .or. mod(domain_nx2,block_nx2)/=0 .or. &
       mod(domain_nx3,block_nx3)/=0) &
       call mpistop('Grad-Rubin domain sizes must be divisible by block sizes')
    if(config%max_flux_imbalance<0.d0 .or. config%max_flux_imbalance>1.d0) &
       call mpistop('Grad-Rubin max_flux_imbalance must be in [0,1]')
    if(config%bz_taper_zero<0.d0 .or. &
       config%bz_taper_full<=config%bz_taper_zero .or. &
       config%bz_taper_full>1.d0) &
       call mpistop('invalid Grad-Rubin Bz taper thresholds')
    if(config%relaxation_factor<=0.d0 .or. config%relaxation_factor>1.d0) &
       call mpistop('Grad-Rubin relaxation factor must be in (0,1]')
    if(config%fieldline_step_fraction<=0.d0 .or. &
       config%fieldline_max_steps<1 .or. config%max_iterations<1 .or. &
       config%convergence_streak<1 .or. config%log_interval<1) &
       call mpistop('invalid Grad-Rubin iteration controls')
    if(config%field_change_tolerance<0.d0 .or. &
       config%energy_change_tolerance<0.d0 .or. config%memory_limit_mb<=0.d0) &
       call mpistop('invalid Grad-Rubin tolerance or memory limit')
    if(config%self_consistency_cycles<0 .or. &
       config%self_consistency_streak<1 .or. &
       config%self_consistency_field_tolerance<0.d0 .or. &
       config%self_consistency_alpha_tolerance<0.d0) &
       call mpistop('invalid Grad-Rubin self-consistency controls')
    if(.not.ieee_is_finite(config%max_flux_imbalance) .or. &
       .not.ieee_is_finite(config%bz_taper_zero) .or. &
       .not.ieee_is_finite(config%bz_taper_full) .or. &
       .not.ieee_is_finite(config%relaxation_factor) .or. &
       .not.ieee_is_finite(config%fieldline_step_fraction) .or. &
       .not.ieee_is_finite(config%field_change_tolerance) .or. &
       .not.ieee_is_finite(config%energy_change_tolerance) .or. &
       .not.ieee_is_finite(config%self_consistency_field_tolerance) .or. &
       .not.ieee_is_finite(config%self_consistency_alpha_tolerance) .or. &
       .not.ieee_is_finite(config%memory_limit_mb)) &
       call mpistop('non-finite Grad-Rubin configuration')
    select case(trim(adjustl(config%flux_treatment)))
    case('strict','subtract_mean')
      continue
    case default
      call mpistop('Grad-Rubin flux_treatment must be strict or subtract_mean')
    end select
  end subroutine validate_configuration

  subroutine extract_external_alpha_core(alpha_raw,alpha_clean,alpha_weight,&
     alpha_mask,alpha_polarity)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    use mod_lfff, only: xa1,xa2

    double precision, intent(out) :: alpha_raw(domain_nx1,domain_nx2)
    double precision, intent(out) :: alpha_clean(domain_nx1,domain_nx2)
    double precision, intent(out) :: alpha_weight(domain_nx1,domain_nx2)
    logical, intent(out) :: alpha_mask(domain_nx1,domain_nx2)
    integer, intent(out) :: alpha_polarity(domain_nx1,domain_nx2)
    double precision :: err,tol
    integer :: i,j,i0,j0,s

    i0=0
    tol=1.d-8*max(1.d0,dabs(xprobmin1),dabs(xprobmax1),dx(1,1))
    do s=1,size(external_alpha_x)-domain_nx1+1
      err=0.d0
      do i=1,domain_nx1
        err=max(err,dabs(external_alpha_x(s+i-1)-&
           (xprobmin1+(dble(i)-0.5d0)*dx(1,1))))
      end do
      if(err<=tol) then
        i0=s
        exit
      end if
    end do
    j0=0
    tol=1.d-8*max(1.d0,dabs(xprobmin2),dabs(xprobmax2),dx(2,1))
    do s=1,size(external_alpha_y)-domain_nx2+1
      err=0.d0
      do j=1,domain_nx2
        err=max(err,dabs(external_alpha_y(s+j-1)-&
           (xprobmin2+(dble(j)-0.5d0)*dx(2,1))))
      end do
      if(err<=tol) then
        j0=s
        exit
      end if
    end do
    if(i0==0 .or. j0==0) call mpistop('external alpha coordinates do not match active grid')
    alpha_raw=external_alpha_raw_full(i0:i0+domain_nx1-1,j0:j0+domain_nx2-1)
    alpha_clean=external_alpha_full(i0:i0+domain_nx1-1,j0:j0+domain_nx2-1)
    alpha_weight=external_alpha_weight_full(i0:i0+domain_nx1-1,j0:j0+domain_nx2-1)
    alpha_mask=external_alpha_valid_full(i0:i0+domain_nx1-1,j0:j0+domain_nx2-1)/=0
    alpha_polarity=external_alpha_polarity_full(i0:i0+domain_nx1-1,j0:j0+domain_nx2-1)
  end subroutine extract_external_alpha_core

  subroutine write_alpha_diagnostics(base,polarity,source,alpha_raw,alpha_selected,&
     alpha_weight,alpha_mask,alpha_polarity)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters, only: mype

    character(len=*), intent(in) :: base,source
    integer, intent(in) :: polarity
    double precision, intent(in) :: alpha_raw(:,:),alpha_selected(:,:),alpha_weight(:,:)
    logical, intent(in) :: alpha_mask(:,:)
    integer, intent(in) :: alpha_polarity(:,:)
    character(len=1024) :: filename,polarity_name
    integer :: iu,i,j,ios

    if(mype/=0) return
    if(polarity>0) then
      polarity_name='positive'
    else
      polarity_name='negative'
    end if
    filename=trim(base)//'_grad_rubin_'//trim(polarity_name)//'_alpha.csv'
    open(newunit=iu,file=trim(filename),status='replace',action='write',iostat=ios)
    if(ios/=0) call mpistop('cannot write Grad-Rubin alpha diagnostics')
    write(iu,'(a)') 'i,j,alpha_raw,alpha_selected,weight,valid,polarity,source'
    do j=1,size(alpha_raw,2)
      do i=1,size(alpha_raw,1)
        write(iu,'(2(i0,","),3(es24.16,","),i0,",",i0,",",a)') i,j,&
           alpha_raw(i,j),alpha_selected(i,j),alpha_weight(i,j),merge(1,0,alpha_mask(i,j)),&
           alpha_polarity(i,j),trim(source)
      end do
    end do
    close(iu)
  end subroutine write_alpha_diagnostics

  subroutine extract_boundary_core(core)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    use mod_lfff, only: xa1,xa2

    double precision, allocatable, intent(out) :: core(:,:,:)
    double precision :: err,tol
    integer :: i,j,i0,j0,s

    i0=0
    tol=1.d-8*max(1.d0,dabs(xprobmin1),dabs(xprobmax1),dx(1,1))
    do s=1,size(boundary_b,1)-domain_nx1+1
      err=0.d0
      do i=1,domain_nx1
        err=max(err,dabs(xa1(s+i-1)-(xprobmin1+(dble(i)-0.5d0)*dx(1,1))))
      end do
      if(err<=tol) then
        i0=s
        exit
      end if
    end do
    j0=0
    tol=1.d-8*max(1.d0,dabs(xprobmin2),dabs(xprobmax2),dx(2,1))
    do s=1,size(boundary_b,2)-domain_nx2+1
      err=0.d0
      do j=1,domain_nx2
        err=max(err,dabs(xa2(s+j-1)-(xprobmin2+(dble(j)-0.5d0)*dx(2,1))))
      end do
      if(err<=tol) then
        j0=s
        exit
      end if
    end do
    if(i0==0 .or. j0==0) call mpistop('Grad-Rubin boundary does not match grid centres')
    allocate(core(domain_nx1,domain_nx2,3))
    core=boundary_b(i0:i0+domain_nx1-1,j0:j0+domain_nx2-1,:)
  end subroutine extract_boundary_core

  subroutine balance_normal_field(bz,treatment,max_imbalance,before,after,&
     correction,status)
    double precision, intent(inout) :: bz(:,:)
    character(len=*), intent(in) :: treatment
    double precision, intent(in) :: max_imbalance
    double precision, intent(out) :: before,after,correction
    integer, intent(out) :: status
    double precision :: unsigned

    unsigned=sum(dabs(bz))
    before=dabs(sum(bz))/max(unsigned,tiny(1.d0))
    after=before
    correction=0.d0
    status=0
    select case(trim(adjustl(treatment)))
    case('strict')
      if(before>1.d-8) status=2
    case('subtract_mean')
      if(before>max_imbalance) then
        status=3
      else if(before>1.d-8) then
        correction=sum(bz)/dble(size(bz))
        bz=bz-correction
        after=dabs(sum(bz))/max(sum(dabs(bz)),tiny(1.d0))
        status=1
      end if
    case default
      status=4
    end select
  end subroutine balance_normal_field

  pure double precision function gr_taper_weight(r,zero_threshold,full_threshold)
    double precision, intent(in) :: r,zero_threshold,full_threshold
    double precision :: pi,s

    pi=4.d0*datan(1.d0)
    if(r<=zero_threshold) then
      gr_taper_weight=0.d0
    else if(r>=full_threshold) then
      gr_taper_weight=1.d0
    else
      s=(r-zero_threshold)/(full_threshold-zero_threshold)
      gr_taper_weight=0.5d0*(1.d0-dcos(pi*s))
    end if
  end function gr_taper_weight

  subroutine gr_compute_alpha_boundary(b,dx1,dx2,zero_threshold,full_threshold,&
     alpha,valid_pixels,alpha_raw)
    double precision, intent(in) :: b(:,:,:),dx1,dx2
    double precision, intent(in) :: zero_threshold,full_threshold
    double precision, intent(out) :: alpha(size(b,1),size(b,2))
    integer, intent(out) :: valid_pixels
    double precision, intent(out), optional :: alpha_raw(size(b,1),size(b,2))
    double precision :: dbxdy,dbydx,bzmax,r,w
    integer :: i,j,nx,ny

    nx=size(b,1)
    ny=size(b,2)
    bzmax=maxval(dabs(b(:,:,3)))
    alpha=0.d0
    if(present(alpha_raw)) alpha_raw=0.d0
    valid_pixels=0
    if(bzmax<=0.d0) return
    do j=1,ny
      do i=1,nx
        dbydx=derivative_1d(b(:,j,2),i,dx1)
        dbxdy=derivative_1d(b(i,:,1),j,dx2)
        r=dabs(b(i,j,3))/bzmax
        w=gr_taper_weight(r,zero_threshold,full_threshold)
        if(w>0.d0 .and. b(i,j,3)/=0.d0) then
          if(present(alpha_raw)) alpha_raw(i,j)=(dbydx-dbxdy)/b(i,j,3)
          alpha(i,j)=w*(dbydx-dbxdy)/b(i,j,3)
          valid_pixels=valid_pixels+1
        end if
      end do
    end do
  end subroutine gr_compute_alpha_boundary

  pure double precision function derivative_1d(f,i,h)
    double precision, intent(in) :: f(:),h
    integer, intent(in) :: i
    integer :: n

    n=size(f)
    if(n<3) then
      derivative_1d=0.d0
    else if(i==1) then
      derivative_1d=(-3.d0*f(1)+4.d0*f(2)-f(3))/(2.d0*h)
    else if(i==n) then
      derivative_1d=(3.d0*f(n)-4.d0*f(n-1)+f(n-2))/(2.d0*h)
    else
      derivative_1d=(f(i+1)-f(i-1))/(2.d0*h)
    end if
  end function derivative_1d

  subroutine transport_alpha_mpi(b,alpha0,bz0,config,alpha,global_counts)
    use mpi
    use mod_global_parameters

    double precision, intent(in) :: b(:,:,:,:),alpha0(:,:),bz0(:,:)
    type(nlfff_grad_rubin_config), intent(in) :: config
    double precision, intent(out) :: alpha(:,:,:)
    integer, intent(out) :: global_counts(4)
    double precision, allocatable :: local_alpha(:),gathered_alpha(:)
    double precision :: x(3),end_forward(3),end_backward(3),value,bzf,bzb
    double precision :: z0,step,bfloor
    integer, allocatable :: recvcounts(:),displs(:)
    integer :: i,j,k,index,first_index,last_index,local_index,nseed
    integer :: statusf,statusb,local_counts(4),chosen,rank

    nseed=size(alpha)
    allocate(recvcounts(npe),displs(npe))
    do rank=0,npe-1
      displs(rank+1)=nseed*rank/npe
      recvcounts(rank+1)=nseed*(rank+1)/npe-displs(rank+1)
    end do
    first_index=displs(mype+1)+1
    last_index=displs(mype+1)+recvcounts(mype+1)
    allocate(local_alpha(recvcounts(mype+1)),gathered_alpha(nseed))
    local_alpha=0.d0
    local_counts=0
    z0=xprobmin3-0.5d0*dx(3,1)
    step=config%fieldline_step_fraction*minval(dx(:,1))
    bfloor=1.d-12*maxval(dsqrt(sum(b**2,dim=4)))
    do index=first_index,last_index
          k=(index-1)/(size(alpha,1)*size(alpha,2))+1
          j=mod(index-1,size(alpha,1)*size(alpha,2))/size(alpha,1)+1
          i=mod(index-1,size(alpha,1))+1
          local_index=index-first_index+1
          x=(/xprobmin1+(dble(i)-0.5d0)*dx(1,1),&
             xprobmin2+(dble(j)-0.5d0)*dx(2,1),z0+dble(k-1)*dx(3,1)/)
          call trace_to_boundary(b,x,1.d0,step,config%fieldline_max_steps,&
             bfloor,end_forward,statusf)
          call trace_to_boundary(b,x,-1.d0,step,config%fieldline_max_steps,&
             bfloor,end_backward,statusb)
          if(statusf==gr_trace_bottom .and. statusb==gr_trace_bottom) then
            local_counts(1)=local_counts(1)+1
            bzf=gr_bilinear_sample(bz0,end_forward(1),end_forward(2),&
               xprobmin1,xprobmin2,dx(1,1),dx(2,1))
            bzb=gr_bilinear_sample(bz0,end_backward(1),end_backward(2),&
               xprobmin1,xprobmin2,dx(1,1),dx(2,1))
            chosen=0
            if(config%polarity*bzf>0.d0) chosen=1
            if(config%polarity*bzb>0.d0) chosen=2
            if(chosen==1) then
              value=gr_bilinear_sample(alpha0,end_forward(1),end_forward(2),&
                 xprobmin1,xprobmin2,dx(1,1),dx(2,1))
              local_alpha(local_index)=value
            else if(chosen==2) then
              value=gr_bilinear_sample(alpha0,end_backward(1),end_backward(2),&
                 xprobmin1,xprobmin2,dx(1,1),dx(2,1))
              local_alpha(local_index)=value
            end if
          else if(statusf==gr_trace_max_steps .or. statusb==gr_trace_max_steps) then
            local_counts(4)=local_counts(4)+1
          else if(statusf==gr_trace_weak .or. statusb==gr_trace_weak) then
            local_counts(3)=local_counts(3)+1
          else
            local_counts(2)=local_counts(2)+1
          end if
    end do
    call MPI_ALLGATHERV(local_alpha,size(local_alpha),MPI_DOUBLE_PRECISION,&
       gathered_alpha,recvcounts,displs,MPI_DOUBLE_PRECISION,icomm,ierrmpi)
    alpha=reshape(gathered_alpha,shape(alpha))
    call MPI_ALLREDUCE(local_counts,global_counts,4,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
    deallocate(local_alpha,gathered_alpha,recvcounts,displs)
  end subroutine transport_alpha_mpi

  subroutine trace_to_boundary(b,start,direction,step,max_steps,bfloor,endpoint,status)
    use mod_global_parameters, only: xprobmin1,xprobmax1,xprobmin2,xprobmax2,&
       xprobmax3,dx,xprobmin3

    double precision, intent(in) :: b(:,:,:,:),start(3),direction,step,bfloor
    integer, intent(in) :: max_steps
    double precision, intent(out) :: endpoint(3)
    integer, intent(out) :: status
    double precision :: x(3),xnew(3),k1(3),k2(3),k3(3),k4(3),z0,t
    logical :: ok
    integer :: n

    z0=xprobmin3-0.5d0*dx(3,1)
    x=start
    do n=1,max_steps
      call field_direction(b,x,direction,bfloor,k1,ok)
      if(.not.ok) then
        endpoint=x
        status=gr_trace_weak
        return
      end if
      if(point_is_open(x+0.5d0*step*k1)) then
        endpoint=x+0.5d0*step*k1
        status=gr_trace_open
        return
      end if
      if(x(3)+0.5d0*step*k1(3)<=z0) then
        t=(z0-x(3))/(0.5d0*step*k1(3))
        endpoint=x+t*(0.5d0*step*k1)
        endpoint(3)=z0
        status=gr_trace_bottom
        return
      end if
      call field_direction(b,x+0.5d0*step*k1,direction,bfloor,k2,ok)
      if(.not.ok) then
        endpoint=x
        status=gr_trace_weak
        return
      end if
      if(point_is_open(x+0.5d0*step*k2)) then
        endpoint=x+0.5d0*step*k2
        status=gr_trace_open
        return
      end if
      if(x(3)+0.5d0*step*k2(3)<=z0) then
        t=(z0-x(3))/(0.5d0*step*k2(3))
        endpoint=x+t*(0.5d0*step*k2)
        endpoint(3)=z0
        status=gr_trace_bottom
        return
      end if
      call field_direction(b,x+0.5d0*step*k2,direction,bfloor,k3,ok)
      if(.not.ok) then
        endpoint=x
        status=gr_trace_weak
        return
      end if
      if(point_is_open(x+step*k3)) then
        endpoint=x+step*k3
        status=gr_trace_open
        return
      end if
      if(x(3)+step*k3(3)<=z0) then
        t=(z0-x(3))/(step*k3(3))
        endpoint=x+t*(step*k3)
        endpoint(3)=z0
        status=gr_trace_bottom
        return
      end if
      call field_direction(b,x+step*k3,direction,bfloor,k4,ok)
      if(.not.ok) then
        endpoint=x
        status=gr_trace_weak
        return
      end if
      xnew=x+step*(k1+2.d0*k2+2.d0*k3+k4)/6.d0
      if(xnew(3)<=z0) then
        t=(z0-x(3))/min(xnew(3)-x(3),-tiny(1.d0))
        endpoint=x+t*(xnew-x)
        endpoint(3)=z0
        status=gr_trace_bottom
        return
      end if
      if(xnew(1)<xprobmin1 .or. xnew(1)>xprobmax1 .or. &
         xnew(2)<xprobmin2 .or. xnew(2)>xprobmax2 .or. xnew(3)>xprobmax3) then
        endpoint=xnew
        status=gr_trace_open
        return
      end if
      x=xnew
    end do
    endpoint=x
    status=gr_trace_max_steps
  end subroutine trace_to_boundary

  logical function point_is_open(x)
    use mod_global_parameters, only: xprobmin1,xprobmax1,xprobmin2,xprobmax2,&
       xprobmax3
    double precision, intent(in) :: x(3)
    point_is_open=x(1)<xprobmin1 .or. x(1)>xprobmax1 .or. &
       x(2)<xprobmin2 .or. x(2)>xprobmax2 .or. x(3)>xprobmax3
  end function point_is_open

  subroutine field_direction(b,x,direction,bfloor,bhat,ok)
    use mod_global_parameters, only: xprobmin1,xprobmax1,xprobmin2,xprobmax2,&
       xprobmin3,xprobmax3,dx
    double precision, intent(in) :: b(:,:,:,:),x(3),direction,bfloor
    double precision, intent(out) :: bhat(3)
    logical, intent(out) :: ok
    double precision :: value(3),norm,z0

    z0=xprobmin3-0.5d0*dx(3,1)
    ! RK stages may lie just across a boundary within the current step. The
    ! sampler clamps those stages; trace_to_boundary classifies the completed
    ! step and computes the exact lower-plane intersection.
    call gr_trilinear_sample(b,x,xprobmin1,xprobmin2,z0,dx(1,1),dx(2,1),&
       dx(3,1),value)
    norm=dsqrt(dot_product(value,value))
    ok=(norm>bfloor)
    if(ok) then
      bhat=direction*value/norm
    else
      bhat=0.d0
    end if
  end subroutine field_direction

  subroutine gr_trilinear_sample(field,x,xmin,ymin,zmin,dx1,dx2,dx3,value)
    double precision, intent(in) :: field(:,:,:,:),x(3),xmin,ymin,zmin,dx1,dx2,dx3
    double precision, intent(out) :: value(size(field,4))
    double precision :: tx,ty,tz,fx,fy,fz
    integer :: i0,j0,k0,i1,j1,k1

    fx=(x(1)-(xmin+0.5d0*dx1))/dx1+1.d0
    fy=(x(2)-(ymin+0.5d0*dx2))/dx2+1.d0
    fz=(x(3)-zmin)/dx3+1.d0
    fx=max(1.d0,min(dble(size(field,1)),fx))
    fy=max(1.d0,min(dble(size(field,2)),fy))
    fz=max(1.d0,min(dble(size(field,3)),fz))
    i0=min(int(floor(fx)),size(field,1)-1); i1=i0+1; tx=fx-dble(i0)
    j0=min(int(floor(fy)),size(field,2)-1); j1=j0+1; ty=fy-dble(j0)
    k0=min(int(floor(fz)),size(field,3)-1); k1=k0+1; tz=fz-dble(k0)
    if(size(field,1)==1) then; i0=1; i1=1; tx=0.d0; end if
    if(size(field,2)==1) then; j0=1; j1=1; ty=0.d0; end if
    if(size(field,3)==1) then; k0=1; k1=1; tz=0.d0; end if
    value=(1.d0-tz)*((1.d0-ty)*((1.d0-tx)*field(i0,j0,k0,:)+&
       tx*field(i1,j0,k0,:))+ty*((1.d0-tx)*field(i0,j1,k0,:)+&
       tx*field(i1,j1,k0,:)))+tz*((1.d0-ty)*((1.d0-tx)*&
       field(i0,j0,k1,:)+tx*field(i1,j0,k1,:))+ty*((1.d0-tx)*&
       field(i0,j1,k1,:)+tx*field(i1,j1,k1,:)))
  end subroutine gr_trilinear_sample

  pure double precision function gr_bilinear_sample(field,x,y,xmin,ymin,dx1,dx2)
    double precision, intent(in) :: field(:,:),x,y,xmin,ymin,dx1,dx2
    double precision :: fx,fy,tx,ty
    integer :: i0,i1,j0,j1

    fx=max(1.d0,min(dble(size(field,1)),(x-(xmin+0.5d0*dx1))/dx1+1.d0))
    fy=max(1.d0,min(dble(size(field,2)),(y-(ymin+0.5d0*dx2))/dx2+1.d0))
    i0=min(int(floor(fx)),size(field,1)-1); i1=i0+1; tx=fx-dble(i0)
    j0=min(int(floor(fy)),size(field,2)-1); j1=j0+1; ty=fy-dble(j0)
    gr_bilinear_sample=(1.d0-ty)*((1.d0-tx)*field(i0,j0)+tx*field(i1,j0))+&
       ty*((1.d0-tx)*field(i0,j1)+tx*field(i1,j1))
  end function gr_bilinear_sample

  subroutine gr_solve_current_field_mpi(current,padding,dx1,dx2,dx3,bc)
    use mpi
    use mod_global_parameters, only: mype,npe,icomm,ierrmpi,dpi
    use mod_fft, only: fft_2d

    double precision, intent(in) :: current(:,:,:,:),dx1,dx2,dx3
    integer, intent(in) :: padding
    double precision, intent(out) :: bc(:,:,:,:)
    double complex, allocatable :: fj(:,:,:,:),fb(:,:,:,:),plane(:,:)
    double complex :: ia(3),ib(3),im(3),sumz,diffx,diffy
    double precision, allocatable :: kx(:),ky(:),z(:),bc_global(:,:,:,:)
    double precision :: q
    integer :: nx,ny,nz,npx,npy,ip0,jp0,i,j,k,s,ic,mode

    nx=size(current,1); ny=size(current,2); nz=size(current,3)
    npx=padding*nx; npy=padding*ny
    ip0=(npx-nx)/2+1; jp0=(npy-ny)/2+1
    allocate(fj(npx,npy,nz,3),fb(npx,npy,nz,3),plane(npx,npy))
    allocate(kx(npx),ky(npy),z(nz))
    fj=(0.d0,0.d0); fb=(0.d0,0.d0)
    do i=1,npx
      mode=i-1; if(mode>npx/2) mode=mode-npx
      kx(i)=2.d0*dpi*dble(mode)/(dble(npx)*dx1)
    end do
    do j=1,npy
      mode=j-1; if(mode>npy/2) mode=mode-npy
      ky(j)=2.d0*dpi*dble(mode)/(dble(npy)*dx2)
    end do
    do k=1,nz; z(k)=dble(k-1)*dx3; end do

    do k=1,nz
      if(mod(k-1,npe)/=mype) cycle
      do ic=1,3
        plane=(0.d0,0.d0)
        plane(ip0:ip0+nx-1,jp0:jp0+ny-1)=&
           dcmplx(current(:,:,k,ic),0.d0)
        call fft_2d(plane,.false.)
        fj(:,:,k,ic)=plane
      end do
    end do
    if(npe>1) then
      call MPI_ALLREDUCE(fj,fb,size(fj),MPI_DOUBLE_COMPLEX,MPI_SUM,icomm,ierrmpi)
      fj=fb
      fb=(0.d0,0.d0)
    end if

    do k=1,nz
      if(mod(k-1,npe)/=mype) cycle
      do j=1,npy
        do i=1,npx
          q=dsqrt(kx(i)**2+ky(j)**2)
          if(q<=0.d0) then
            ia=(0.d0,0.d0)
            if(k<nz) then
              do s=k,nz-1
                ia=ia+0.5d0*dx3*(fj(i,j,s,:)+fj(i,j,s+1,:))
              end do
            end if
            fb(i,j,k,1)=-ia(2)
            fb(i,j,k,2)= ia(1)
            fb(i,j,k,3)=(0.d0,0.d0)
          else
            ia=(0.d0,0.d0); ib=(0.d0,0.d0); im=(0.d0,0.d0)
            if(k<nz) then
              do s=k,nz-1
                ia=ia+0.5d0*dx3*(dexp(-q*(z(s)-z(k)))*fj(i,j,s,:)+&
                   dexp(-q*(z(s+1)-z(k)))*fj(i,j,s+1,:))
              end do
            end if
            if(k>1) then
              do s=1,k-1
                ib=ib+0.5d0*dx3*(dexp(-q*(z(k)-z(s)))*fj(i,j,s,:)+&
                   dexp(-q*(z(k)-z(s+1)))*fj(i,j,s+1,:))
              end do
            end if
            do s=1,nz-1
              im=im+0.5d0*dx3*(dexp(-q*(z(k)+z(s)))*fj(i,j,s,:)+&
                 dexp(-q*(z(k)+z(s+1)))*fj(i,j,s+1,:))
            end do
            sumz=ia(3)+ib(3)+im(3)
            diffx=ia(1)-ib(1)+im(1)
            diffy=ia(2)-ib(2)+im(2)
            ! With mod_fft's exp(-i k.x) forward convention, curl(A) uses
            ! +i ky Az-dAy/dz, dAx/dz-i kx Az, and
            ! +i kx Ay-i ky Ax.  These signs are also the ones used by the
            ! Wheatland/CFIT open-half-space solution.
            fb(i,j,k,1)= (0.d0,1.d0)*ky(j)*sumz/(2.d0*q)-0.5d0*diffy
            fb(i,j,k,2)= 0.5d0*diffx-(0.d0,1.d0)*kx(i)*sumz/(2.d0*q)
            fb(i,j,k,3)=((0.d0,1.d0)*kx(i)*(ia(2)+ib(2)-im(2))-&
               (0.d0,1.d0)*ky(j)*(ia(1)+ib(1)-im(1)))/(2.d0*q)
          end if
        end do
      end do
    end do
    if(npe>1) then
      call MPI_ALLREDUCE(fb,fj,size(fb),MPI_DOUBLE_COMPLEX,MPI_SUM,icomm,ierrmpi)
      fb=fj
    end if
    bc=0.d0
    do k=1,nz
      if(mod(k-1,npe)/=mype) cycle
      do ic=1,3
        plane=fb(:,:,k,ic)
        call fft_2d(plane,.true.)
        bc(:,:,k,ic)=dble(plane(ip0:ip0+nx-1,jp0:jp0+ny-1))
      end do
    end do
    if(npe>1) then
      allocate(bc_global(nx,ny,nz,3))
      call MPI_ALLREDUCE(bc,bc_global,size(bc),MPI_DOUBLE_PRECISION,MPI_SUM,&
         icomm,ierrmpi)
      bc=bc_global
      deallocate(bc_global)
    end if
    bc(:,:,1,3)=0.d0
    deallocate(fj,fb,plane,kx,ky,z)
  end subroutine gr_solve_current_field_mpi

  subroutine scatter_dense_field(b,iw_b)
    use mod_global_parameters
    use mod_forest, only: tree_root
    double precision, intent(in) :: b(:,:,:,:)
    integer, intent(in) :: iw_b(3)
    integer :: ig1,ig2,ig3,igrid,i,j,k,ic,gi,gj,gk

    do ig3=1,domain_nx3/block_nx3
      do ig2=1,domain_nx2/block_nx2
        do ig1=1,domain_nx1/block_nx1
          if(tree_root(ig1,ig2,ig3)%node%ipe/=mype) cycle
          igrid=tree_root(ig1,ig2,ig3)%node%igrid
          do k=1,block_nx3; gk=(ig3-1)*block_nx3+k+1
            do j=1,block_nx2; gj=(ig2-1)*block_nx2+j
              do i=1,block_nx1; gi=(ig1-1)*block_nx1+i
                do ic=1,3
                  ps(igrid)%w(ixMlo1+i-1,ixMlo2+j-1,ixMlo3+k-1,iw_b(ic))=&
                     b(gi,gj,gk,ic)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine scatter_dense_field

  pure double precision function gr_energy(b,dx1,dx2,dx3)
    double precision, intent(in) :: b(:,:,:,:),dx1,dx2,dx3
    gr_energy=0.5d0*sum(b(:,:,2:,:)**2)*dx1*dx2*dx3
  end function gr_energy

  subroutine field_diagnostics(b,dx1,dx2,dx3,theta_deg,epsilon_div)
    double precision, intent(in) :: b(:,:,:,:),dx1,dx2,dx3
    double precision, intent(out) :: theta_deg,epsilon_div
    double precision :: db(3,3),jvec(3),bv(3),cross(3),divb,b2
    double precision :: sum_cross,sum_j,sum_div2,sum_b2,pi,current_scale,h
    double precision :: ncells
    integer :: i,j,k,ic

    sum_cross=0.d0; sum_j=0.d0; sum_div2=0.d0; sum_b2=0.d0
    do k=2,size(b,3)
      do j=1,size(b,2)
        do i=1,size(b,1)
          do ic=1,3
            db(ic,1)=derivative_1d(b(:,j,k,ic),i,dx1)
            db(ic,2)=derivative_1d(b(i,:,k,ic),j,dx2)
            db(ic,3)=derivative_1d(b(i,j,:,ic),k,dx3)
          end do
          jvec=(/db(3,2)-db(2,3),db(1,3)-db(3,1),db(2,1)-db(1,2)/)
          bv=b(i,j,k,:)
          cross=(/jvec(2)*bv(3)-jvec(3)*bv(2),&
             jvec(3)*bv(1)-jvec(1)*bv(3),jvec(1)*bv(2)-jvec(2)*bv(1)/)
          b2=dot_product(bv,bv)
          if(b2>0.d0) sum_cross=sum_cross+dsqrt(dot_product(cross,cross)/b2)
          sum_j=sum_j+dsqrt(dot_product(jvec,jvec))
          divb=db(1,1)+db(2,2)+db(3,3)
          sum_div2=sum_div2+divb**2
          sum_b2=sum_b2+b2
        end do
      end do
    end do
    pi=4.d0*datan(1.d0)
    ncells=dble(size(b,1)*size(b,2)*(size(b,3)-1))
    h=(dx1*dx2*dx3)**(1.d0/3.d0)
    current_scale=dsqrt(max(sum_b2*ncells,tiny(1.d0)))/h
    if(sum_j<=1.d-12*current_scale) then
      theta_deg=0.d0
    else
      theta_deg=dasin(min(1.d0,sum_cross/sum_j))*180.d0/pi
    end if
    epsilon_div=h*dsqrt(sum_div2/max(sum_b2,tiny(1.d0)))
  end subroutine field_diagnostics

  double precision function estimate_memory_mb(padding)
    use mod_global_parameters, only: domain_nx1,domain_nx2,domain_nx3,npe
    integer, intent(in) :: padding
    double precision :: n,nxy,nfft,nxyfft,boundary_n,base_bytes
    double precision :: trace_bytes,fft_bytes,potential_bytes

    n=dble(domain_nx1)*dble(domain_nx2)*dble(domain_nx3+1)
    nxy=dble(domain_nx1)*dble(domain_nx2)
    nxyfft=dble(padding*domain_nx1)*dble(padding*domain_nx2)
    nfft=dble(padding*domain_nx1)*dble(padding*domain_nx2)*&
       dble(domain_nx3+1)
    boundary_n=0.d0
    if(allocated(boundary_b)) boundary_n=dble(size(boundary_b))
    ! Persistent arrays: input boundary, extracted vector boundary, alpha0,
    ! five dense vector fields, and dense alpha.  The two alternatives add
    ! either Allgatherv buffers or the padded complex Poisson work arrays.
    base_bytes=8.d0*(boundary_n+5.d0*nxy+16.d0*n)
    trace_bytes=8.d0*(2.d0*n+dble(ceiling(n/dble(npe))))+&
       4.d0*dble(2*npe)
    fft_bytes=16.d0*(6.d0*nfft+nxyfft)+max(8.d0*3.d0*n,16.d0*nxyfft)+&
       8.d0*dble(padding*domain_nx1+padding*domain_nx2+domain_nx3+1)
    ! During the initial potential solve only potential and four padded real
    ! planes coexist with the boundary arrays.
    potential_bytes=8.d0*(boundary_n+5.d0*nxy+3.d0*n+4.d0*nxyfft+&
       dble(padding*domain_nx1+padding*domain_nx2))
    estimate_memory_mb=max(base_bytes+trace_bytes,base_bytes+fft_bytes,&
       potential_bytes)/(1024.d0**2)
  end function estimate_memory_mb

  double precision function estimate_self_consistent_memory_mb(padding)
    use mod_global_parameters, only: domain_nx1,domain_nx2,domain_nx3
    integer, intent(in) :: padding
    double precision :: n,nxy,extra_bytes

    n=dble(domain_nx1)*dble(domain_nx2)*dble(domain_nx3+1)
    nxy=dble(domain_nx1)*dble(domain_nx2)
    ! In addition to a single-polarity solve, self-consistency retains the
    ! potential field, both P/N fields and both transported-alpha volumes.
    extra_bytes=8.d0*(9.d0*n+3.d0*nxy)
    estimate_self_consistent_memory_mb=estimate_memory_mb(padding)+&
       extra_bytes/(1024.d0**2)
  end function estimate_self_consistent_memory_mb
}

end module mod_nlfff_grad_rubin
