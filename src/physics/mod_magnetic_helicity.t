!> Finite-volume relative magnetic helicity on a uniform Cartesian mesh.
!>
!> The implementation uses a DeVore gauge and streams one plane at a time.
!> Plane data are reduced to a deterministic owner; consecutive plane owners
!> exchange only the O(N_perp**2) integration state. No magnetic-field-line
!> tracing or replicated three-dimensional field is used.
module mod_magnetic_helicity
  use mod_magnetic_reference_fv, only: magnetic_reference_config,&
     magnetic_reference_result,magnetic_reference_field,&
     solve_magnetic_reference_fv,free_magnetic_reference_field
  implicit none
  private

  integer, parameter :: mh_name_len=256
  integer, parameter :: mh_nstate=24
  integer, parameter :: mh_state_tag=28431
  integer, parameter :: mh_debug_tag=28432
  integer, parameter :: mh_debug_nvalue=12
  integer, parameter :: mh_debug_narray=6
  ! Plane-pipeline channel starts: two integration levels for B and Bp,
  ! two magnetic planes for centered derivatives, one Bp plane, and top b.
  integer, parameter :: mh_q1_b=1,mh_q2_b=4
  integer, parameter :: mh_q1_bp=7,mh_q2_bp=10
  integer, parameter :: mh_b1=13,mh_b2=16,mh_bp1=19,mh_btop=22

  type, public :: magnetic_helicity_config
    integer :: gauge_axis=3
    double precision :: ratio_tolerance=1.d-12
    logical :: write_debug_vti=.false.
    logical :: write_mg_timing=.false.
    character(len=mh_name_len) :: debug_vti_file=''
    type(magnetic_reference_config) :: reference
  end type magnetic_helicity_config

  type, public :: magnetic_helicity_result
    double precision :: Hm=0.d0,HJ=0.d0,HPJ=0.d0
    double precision :: abs_HJ_over_abs_Hm=0.d0
    double precision :: energy=0.d0,potential_energy=0.d0,free_energy=0.d0
    double precision :: Hm_physical=0.d0,HJ_physical=0.d0
    double precision :: HPJ_physical=0.d0
    double precision :: energy_physical=0.d0
    double precision :: potential_energy_physical=0.d0
    double precision :: free_energy_physical=0.d0
    double precision :: helicity_unit=1.d0,energy_unit=1.d0
    double precision :: curl_A_error=0.d0,curl_Ap_error=0.d0
    double precision :: decomposition_error=0.d0
    double precision :: epsilon_div_B=0.d0
    double precision :: net_flux_imbalance=0.d0
    double precision :: boundary_normal_error=0.d0
    double precision :: mg_residual=0.d0
    integer :: mg_cycles=0
    logical :: ratio_is_valid=.false.
  end type magnetic_helicity_result

  type, private :: mh_debug_vti_state
    logical :: enabled=.false.
    integer :: unit=0
    integer :: nxyz(3)=0
    integer(kind=8) :: data_pos=0_8
    integer(kind=8) :: offsets(mh_debug_narray)=0_8
    integer(kind=8) :: nbytes(mh_debug_narray)=0_8
  end type mh_debug_vti_state

  character(len=mh_name_len), public :: mh_output_file=''
  integer, public :: mh_gauge_axis=3
  double precision, public :: mh_mg_tolerance=1.d-4
  integer, public :: mh_mg_max_cycles=50
  double precision, public :: mh_max_flux_imbalance=1.d-6
  logical, public :: mh_write_debug_vti=.false.
  logical, public :: mh_write_mg_timing=.false.
  character(len=mh_name_len), public :: mh_debug_vti_file=''
  type(magnetic_helicity_result), public, save :: mh_last_result

  public :: mh_params_read,mh_run_task,compute_magnetic_helicity_fv

contains

  subroutine mh_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /magnetic_helicity_list/ mh_output_file,&
       mh_gauge_axis,mh_mg_tolerance,mh_mg_max_cycles,&
       mh_max_flux_imbalance,mh_write_debug_vti,mh_write_mg_timing,&
       mh_debug_vti_file

    mh_output_file=''
    mh_gauge_axis=3
    mh_mg_tolerance=1.d-4
    mh_mg_max_cycles=50
    mh_max_flux_imbalance=1.d-6
    mh_write_debug_vti=.false.
    mh_write_mg_timing=.false.
    mh_debug_vti_file=''
    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,magnetic_helicity_list,end=111)
111   close(unitpar)
    end do
  end subroutine mh_params_read

  subroutine mh_run_task()
    use mod_comm_lib, only: mpistop
    use mod_global_parameters, only: par_files,convert,level_io,global_time,&
         ps,iwstart,nwgc
    use mod_ghostcells_update, only: getbc
    type(magnetic_helicity_config) :: config

    if(.not.convert) call mpistop(&
         'magnetic-helicity conversion requires convert=.true.')
    if(level_io<1) call mpistop(&
         'magnetic-helicity conversion requires level_io > 0')
    call mh_params_read(par_files)
    call getbc(global_time,0.d0,ps,iwstart,nwgc)
    config%gauge_axis=mh_gauge_axis
    config%reference%residual_tolerance=mh_mg_tolerance
    config%reference%max_cycles=mh_mg_max_cycles
    config%reference%max_flux_imbalance=mh_max_flux_imbalance
    config%reference%write_timing=mh_write_mg_timing
    config%write_debug_vti=mh_write_debug_vti
    config%debug_vti_file=mh_debug_vti_file
    call compute_magnetic_helicity_fv(config,mh_last_result)
    call mh_write_result_csv(mh_last_result)
  end subroutine mh_run_task

  subroutine mh_write_result_csv(result)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters, only: mype,base_filename,snapshotini,&
         global_time,time_convert_factor,SI_unit
    type(magnetic_helicity_result), intent(in) :: result
    character(len=mh_name_len) :: filename
    character(len=5) :: helicity_label
    character(len=3) :: energy_label
    character(len=5) :: ratio_label
    logical :: exists
    integer :: csv_unit,io_status
    integer(kind=8) :: file_size

    if(mype/=0) return
    if(len_trim(mh_output_file)>0) then
      filename=trim(mh_output_file)
    else
      filename=trim(base_filename)//'_helicity.csv'
    end if
    if(SI_unit) then
      helicity_label='Wb2'
      energy_label='J'
    else
      helicity_label='Mx2'
      energy_label='erg'
    end if
    ratio_label=merge('true ','false',result%ratio_is_valid)

    file_size=0_8
    inquire(file=trim(filename),exist=exists)
    if(exists) inquire(file=trim(filename),size=file_size)
    open(newunit=csv_unit,file=trim(filename),status='unknown',&
         action='write',position='append',iostat=io_status)
    if(io_status/=0) call mpistop('could not open magnetic-helicity CSV')
    if(.not.exists .or. file_size==0) then
      write(csv_unit,'(a)') 'snapshot,code_time,physical_time_s,'//&
           'Hm_code,HJ_code,HPJ_code,abs_HJ_over_abs_Hm,'//&
           'Hm_'//trim(helicity_label)//',HJ_'//trim(helicity_label)//&
           ',HPJ_'//trim(helicity_label)//',E_code,Ep_code,Efree_code,'//&
           'E_'//trim(energy_label)//',Ep_'//trim(energy_label)//&
           ',Efree_'//trim(energy_label)//',net_flux_imbalance,'//&
           'epsilon_div_B,boundary_normal_error,curl_A_error,'//&
           'curl_Ap_error,decomposition_error,mg_residual,mg_cycles,'//&
           'ratio_is_valid'
    end if

    write(csv_unit,'(i0)',advance='no') snapshotini
    call mh_csv_real(csv_unit,global_time)
    call mh_csv_real(csv_unit,global_time*time_convert_factor)
    call mh_csv_real(csv_unit,result%Hm)
    call mh_csv_real(csv_unit,result%HJ)
    call mh_csv_real(csv_unit,result%HPJ)
    if(result%ratio_is_valid) then
      call mh_csv_real(csv_unit,result%abs_HJ_over_abs_Hm)
    else
      write(csv_unit,'(",")',advance='no')
    end if
    call mh_csv_real(csv_unit,result%Hm_physical)
    call mh_csv_real(csv_unit,result%HJ_physical)
    call mh_csv_real(csv_unit,result%HPJ_physical)
    call mh_csv_real(csv_unit,result%energy)
    call mh_csv_real(csv_unit,result%potential_energy)
    call mh_csv_real(csv_unit,result%free_energy)
    call mh_csv_real(csv_unit,result%energy_physical)
    call mh_csv_real(csv_unit,result%potential_energy_physical)
    call mh_csv_real(csv_unit,result%free_energy_physical)
    call mh_csv_real(csv_unit,result%net_flux_imbalance)
    call mh_csv_real(csv_unit,result%epsilon_div_B)
    call mh_csv_real(csv_unit,result%boundary_normal_error)
    call mh_csv_real(csv_unit,result%curl_A_error)
    call mh_csv_real(csv_unit,result%curl_Ap_error)
    call mh_csv_real(csv_unit,result%decomposition_error)
    call mh_csv_real(csv_unit,result%mg_residual)
    write(csv_unit,'(",",i0,",",a)') result%mg_cycles,trim(ratio_label)
    close(csv_unit)
  end subroutine mh_write_result_csv

  subroutine mh_csv_real(csv_unit,value)
    integer, intent(in) :: csv_unit
    double precision, intent(in) :: value
    write(csv_unit,'(",",es24.16)',advance='no') value
  end subroutine mh_csv_real

  subroutine compute_magnetic_helicity_fv(config,result)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    use mod_geometry, only: Cartesian,coordinate
    type(magnetic_helicity_config), intent(in) :: config
    type(magnetic_helicity_result), intent(out) :: result

    result=magnetic_helicity_result()
    if(ndim/=3) call mpistop('finite-volume magnetic helicity requires three dimensions')
    if(coordinate/=Cartesian) &
       call mpistop('finite-volume magnetic helicity requires Cartesian coordinates')
    if(any(stretched_dim)) &
       call mpistop('finite-volume magnetic helicity requires an unstretched mesh')
    if(config%gauge_axis<1 .or. config%gauge_axis>3) &
       call mpistop('mh_gauge_axis must be 1, 2, or 3')
    if(config%ratio_tolerance<=0.d0) &
       call mpistop('magnetic-helicity ratio tolerance must be positive')

    {^IFTHREED
    call mh_compute_3d(config,result)
    }
  end subroutine compute_magnetic_helicity_fv

{^IFTHREED
  subroutine mh_compute_3d(config,result)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    type(magnetic_helicity_config), intent(in) :: config
    type(magnetic_helicity_result), intent(inout) :: result

    type(magnetic_reference_result) :: ref_result
    type(magnetic_reference_field) :: bp
    type(mh_debug_vti_state) :: debug_vti
    double precision, allocatable :: send_b(:,:,:),recv_b(:,:,:)
    double precision, allocatable :: send_bp(:,:,:),recv_bp(:,:,:)
    double precision, allocatable :: avec(:,:,:),apvec(:,:,:)
    double precision, allocatable :: q_b(:,:,:),q_bp(:,:,:)
    double precision, allocatable :: bdiff(:,:,:)
    double precision, allocatable :: plane_state(:,:,:)
    double precision, allocatable :: debug_plane(:,:)
    double precision, allocatable :: face_send(:,:),face_recv(:,:)
    integer, allocatable :: count_send(:,:),count_recv(:,:)
    integer, allocatable :: face_count_send(:,:),face_count_recv(:,:)
    double precision :: local_values(11),global_values(11)
    double precision :: cell_dx(3),du,dv,da,dvolume,denom
    integer :: nxyz(3),axis,uaxis,vaxis,orientation
    integer :: nu,nv,na,k,owner,previous_owner,next_owner
    integer :: request,status(MPI_STATUS_SIZE),state_count

    call mh_uniform_grid(config%gauge_axis,nxyz,cell_dx,uaxis,vaxis,&
       orientation)
    axis=config%gauge_axis
    nu=nxyz(uaxis)
    nv=nxyz(vaxis)
    na=nxyz(axis)
    du=cell_dx(uaxis)
    dv=cell_dx(vaxis)
    da=cell_dx(axis)
    dvolume=cell_dx(1)*cell_dx(2)*cell_dx(3)
    if(min(nu,nv,na)<3) &
       call mpistop('magnetic helicity requires at least three cells per direction')

    call mh_debug_vti_begin(debug_vti,config,nxyz,cell_dx)

    call solve_magnetic_reference_fv(config%reference,bp,ref_result)

    allocate(send_b(nu,nv,3),recv_b(nu,nv,3))
    allocate(send_bp(nu,nv,3),recv_bp(nu,nv,3))
    allocate(avec(nu,nv,3),apvec(nu,nv,3))
    allocate(q_b(nu,nv,3),q_bp(nu,nv,3))
    allocate(bdiff(nu,nv,3))
    allocate(plane_state(nu,nv,mh_nstate))
    allocate(face_send(nu,nv),face_recv(nu,nv))
    allocate(count_send(nu,nv),count_recv(nu,nv))
    allocate(face_count_send(nu,nv),face_count_recv(nu,nv))
    if(debug_vti%enabled) allocate(debug_plane(mh_debug_nvalue,nu*nv))

    plane_state=0.d0
    local_values=0.d0
    state_count=nu*nv*mh_nstate
    do k=na,1,-1
      owner=mh_plane_owner(k)
      request=MPI_REQUEST_NULL
      if(k<na .and. mype==owner) then
        previous_owner=mh_plane_owner(k+1)
        if(previous_owner/=owner) call MPI_IRECV(plane_state,state_count,&
           MPI_DOUBLE_PRECISION,previous_owner,mh_state_tag,icomm,request,&
           ierrmpi)
      end if

      call mh_assemble_plane(axis,uaxis,vaxis,k,nxyz,cell_dx,bp,send_b,&
         send_bp,count_send,face_send,face_count_send)
      call MPI_REDUCE(send_b,recv_b,nu*nv*3,MPI_DOUBLE_PRECISION,&
         MPI_SUM,owner,icomm,ierrmpi)
      call MPI_REDUCE(send_bp,recv_bp,nu*nv*3,MPI_DOUBLE_PRECISION,&
         MPI_SUM,owner,icomm,ierrmpi)
      call MPI_REDUCE(count_send,count_recv,nu*nv,MPI_INTEGER,&
         MPI_SUM,owner,icomm,ierrmpi)
      if(k==na) then
        call MPI_REDUCE(face_send,face_recv,nu*nv,MPI_DOUBLE_PRECISION,&
           MPI_SUM,owner,icomm,ierrmpi)
        call MPI_REDUCE(face_count_send,face_count_recv,nu*nv,MPI_INTEGER,&
           MPI_SUM,owner,icomm,ierrmpi)
      end if

      if(mype==owner) then
        if(request/=MPI_REQUEST_NULL) call MPI_WAIT(request,status,ierrmpi)
        if(any(count_recv/=1)) &
           call mpistop('magnetic-helicity plane assembly is not one-to-one')
        if(k==na) then
          if(any(face_count_recv/=1)) call mpistop(&
             'magnetic-helicity top boundary assembly is not one-to-one')
          plane_state=0.d0
          call mh_build_top_b(face_recv,du,dv,uaxis,vaxis,axis,&
             orientation,plane_state(:,:,mh_btop:mh_btop+2))
          q_b=0.5d0*da*recv_b
          q_bp=0.5d0*da*recv_bp
        else if(k==na-1) then
          ! Seed the second plane by a trapezoid between cell centers.
          q_b=plane_state(:,:,mh_q1_b:mh_q1_b+2)+0.5d0*da*(&
             plane_state(:,:,mh_b1:mh_b1+2)+recv_b)
          q_bp=plane_state(:,:,mh_q1_bp:mh_q1_bp+2)+0.5d0*da*(&
             plane_state(:,:,mh_bp1:mh_bp1+2)+recv_bp)
        else
          ! I(k)=I(k+2)+2*ds*B(k+1), so the centered derivative of I
          ! is exactly -B at every interior plane.
          q_b=plane_state(:,:,mh_q2_b:mh_q2_b+2)+2.d0*da*&
             plane_state(:,:,mh_b1:mh_b1+2)
          q_bp=plane_state(:,:,mh_q2_bp:mh_q2_bp+2)+2.d0*da*&
             plane_state(:,:,mh_bp1:mh_bp1+2)
        end if

        call mh_make_vector_potential(&
           plane_state(:,:,mh_btop:mh_btop+2),q_b,uaxis,&
           vaxis,axis,orientation,avec)
        call mh_make_vector_potential(&
           plane_state(:,:,mh_btop:mh_btop+2),q_bp,uaxis,&
           vaxis,axis,orientation,apvec)
        if(k<=na-2) call mh_plane_quality(plane_state,q_b,q_bp,recv_b,&
           uaxis,vaxis,axis,orientation,du,dv,da,dvolume,local_values)
        call mh_plane_integrals(avec,apvec,recv_b,recv_bp,bdiff,&
           dvolume,local_values)
        if(debug_vti%enabled) call mh_debug_pack_plane(recv_bp,avec,apvec,&
             bdiff,debug_plane)

        plane_state(:,:,mh_q2_b:mh_q2_b+2)=&
           plane_state(:,:,mh_q1_b:mh_q1_b+2)
        plane_state(:,:,mh_q1_b:mh_q1_b+2)=q_b
        plane_state(:,:,mh_q2_bp:mh_q2_bp+2)=&
           plane_state(:,:,mh_q1_bp:mh_q1_bp+2)
        plane_state(:,:,mh_q1_bp:mh_q1_bp+2)=q_bp
        plane_state(:,:,mh_b2:mh_b2+2)=plane_state(:,:,mh_b1:mh_b1+2)
        plane_state(:,:,mh_b1:mh_b1+2)=recv_b
        plane_state(:,:,mh_bp1:mh_bp1+2)=recv_bp

      end if
      if(debug_vti%enabled) call mh_debug_vti_write_plane(debug_vti,owner,&
           k,uaxis,vaxis,axis,debug_plane)
      if(mype==owner .and. k>1) then
        next_owner=mh_plane_owner(k-1)
        if(next_owner/=owner) call MPI_SEND(plane_state,state_count,&
             MPI_DOUBLE_PRECISION,next_owner,mh_state_tag,icomm,ierrmpi)
      end if
    end do
    call mh_debug_vti_end(debug_vti)

    call MPI_ALLREDUCE(local_values,global_values,size(local_values),&
       MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    result%Hm=global_values(1)
    result%HJ=global_values(2)
    result%HPJ=global_values(3)
    result%energy=global_values(4)
    result%potential_energy=ref_result%magnetic_energy
    result%free_energy=result%energy-result%potential_energy
    result%curl_A_error=dsqrt(global_values(6)/&
       max(global_values(7),tiny(1.d0)))
    result%curl_Ap_error=dsqrt(global_values(8)/&
       max(global_values(9),tiny(1.d0)))
    result%epsilon_div_B=dsqrt(global_values(10)/&
       max(global_values(11),tiny(1.d0)))
    denom=max(abs(result%Hm),abs(result%HJ)+abs(result%HPJ),tiny(1.d0))
    result%decomposition_error=abs(result%Hm-result%HJ-result%HPJ)/denom
    result%ratio_is_valid=abs(result%Hm)>&
       config%ratio_tolerance*max(global_values(5),tiny(1.d0))
    if(result%ratio_is_valid) &
       result%abs_HJ_over_abs_Hm=abs(result%HJ)/abs(result%Hm)

    result%net_flux_imbalance=ref_result%flux_imbalance
    result%boundary_normal_error=ref_result%boundary_normal_error
    result%mg_residual=ref_result%residual
    result%mg_cycles=ref_result%cycles
    result%helicity_unit=unit_magneticfield**2*unit_length**4
    result%energy_unit=unit_pressure*unit_length**3
    result%Hm_physical=result%Hm*result%helicity_unit
    result%HJ_physical=result%HJ*result%helicity_unit
    result%HPJ_physical=result%HPJ*result%helicity_unit
    result%energy_physical=result%energy*result%energy_unit
    result%potential_energy_physical=result%potential_energy*result%energy_unit
    result%free_energy_physical=result%free_energy*result%energy_unit

    call free_magnetic_reference_field(bp)
    deallocate(send_b,recv_b,send_bp,recv_bp,avec,apvec,q_b,q_bp,&
       bdiff,plane_state)
    deallocate(face_send,face_recv,count_send,count_recv)
    deallocate(face_count_send,face_count_recv)
    if(allocated(debug_plane)) deallocate(debug_plane)
  end subroutine mh_compute_3d

  integer function mh_plane_owner(k)
    use mod_global_parameters, only: npe
    integer, intent(in) :: k
    mh_plane_owner=mod(k-1,npe)
  end function mh_plane_owner

  subroutine mh_uniform_grid(axis,nxyz,cell_dx,uaxis,vaxis,orientation)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    integer, intent(in) :: axis
    integer, intent(out) :: nxyz(3),uaxis,vaxis,orientation
    double precision, intent(out) :: cell_dx(3)
    integer :: iigrid,igrid,local_min,local_max,grid_min,grid_max
    double precision :: local_error,global_error

    local_min=huge(1)
    local_max=-huge(1)
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      local_min=min(local_min,node(plevel_,igrid))
      local_max=max(local_max,node(plevel_,igrid))
    end do
    call MPI_ALLREDUCE(local_min,grid_min,1,MPI_INTEGER,MPI_MIN,icomm,ierrmpi)
    call MPI_ALLREDUCE(local_max,grid_max,1,MPI_INTEGER,MPI_MAX,icomm,ierrmpi)
    if(grid_min/=grid_max) call mpistop(&
       'magnetic helicity requires one uniform AMR level; set level_io')
    nxyz=(/domain_nx1,domain_nx2,domain_nx3/)*2**(grid_min-1)
    cell_dx(1)=(xprobmax1-xprobmin1)/dble(nxyz(1))
    cell_dx(2)=(xprobmax2-xprobmin2)/dble(nxyz(2))
    cell_dx(3)=(xprobmax3-xprobmin3)/dble(nxyz(3))
    local_error=0.d0
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      local_error=max(local_error,abs(rnode(rpdx1_,igrid)-cell_dx(1))/&
         cell_dx(1))
      local_error=max(local_error,abs(rnode(rpdx2_,igrid)-cell_dx(2))/&
         cell_dx(2))
      local_error=max(local_error,abs(rnode(rpdx3_,igrid)-cell_dx(3))/&
         cell_dx(3))
    end do
    call MPI_ALLREDUCE(local_error,global_error,1,MPI_DOUBLE_PRECISION,&
       MPI_MAX,icomm,ierrmpi)
    if(global_error>1.d-12) call mpistop(&
       'magnetic helicity requires a uniform Cartesian cell size')

    select case(axis)
    case(1)
      uaxis=2; vaxis=3; orientation=1
    case(2)
      uaxis=1; vaxis=3; orientation=-1
    case(3)
      uaxis=1; vaxis=2; orientation=1
    end select
  end subroutine mh_uniform_grid

  subroutine mh_assemble_plane(axis,uaxis,vaxis,kplane,nxyz,cell_dx,bp,&
       bsend,bpsend,counts,face_send,face_counts)
    use mod_global_parameters
    integer, intent(in) :: axis,uaxis,vaxis,kplane,nxyz(3)
    double precision, intent(in) :: cell_dx(3)
    type(magnetic_reference_field), intent(in) :: bp
    double precision, intent(out) :: bsend(:,:,:),bpsend(:,:,:)
    integer, intent(out) :: counts(:,:)
    double precision, intent(out) :: face_send(:,:)
    integer, intent(out) :: face_counts(:,:)
    integer :: iigrid,igrid,i,j,l,idir,ig,jg,lg,iloc,jloc,lloc
    integer :: start(3),block_extent(3),gu,gv

    bsend=0.d0
    bpsend=0.d0
    counts=0
    face_send=0.d0
    face_counts=0
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      start(1)=nint((ps(igrid)%x(ixMlo1,ixMlo2,ixMlo3,1)-&
         xprobmin1)/cell_dx(1)+0.5d0)
      start(2)=nint((ps(igrid)%x(ixMlo1,ixMlo2,ixMlo3,2)-&
         xprobmin2)/cell_dx(2)+0.5d0)
      start(3)=nint((ps(igrid)%x(ixMlo1,ixMlo2,ixMlo3,3)-&
         xprobmin3)/cell_dx(3)+0.5d0)
      block_extent=(/ixMhi1-ixMlo1,ixMhi2-ixMlo2,ixMhi3-ixMlo3/)
      if(kplane<start(axis) .or. &
         kplane>start(axis)+block_extent(axis)) cycle
      select case(axis)
      case(1)
        iloc=ixMlo1+kplane-start(1)
        do l=ixMlo3,ixMhi3; lg=start(3)+l-ixMlo3
          do j=ixMlo2,ixMhi2; jg=start(2)+j-ixMlo2
            gu=jg; gv=lg
            do idir=1,3
              bsend(gu,gv,idir)=mh_total_b(igrid,iloc,j,l,idir)
              bpsend(gu,gv,idir)=bp%blocks(igrid)%b(iloc,j,l,idir)
            end do
            counts(gu,gv)=1
            if(kplane==nxyz(axis) .and.&
               ps(igrid)%is_physical_boundary(2*axis)) then
              face_send(gu,gv)=mh_boundary_normal(igrid,iloc,j,l,axis)
              face_counts(gu,gv)=1
            end if
          end do
        end do
      case(2)
        jloc=ixMlo2+kplane-start(2)
        do l=ixMlo3,ixMhi3; lg=start(3)+l-ixMlo3
          do i=ixMlo1,ixMhi1; ig=start(1)+i-ixMlo1
            gu=ig; gv=lg
            do idir=1,3
              bsend(gu,gv,idir)=mh_total_b(igrid,i,jloc,l,idir)
              bpsend(gu,gv,idir)=bp%blocks(igrid)%b(i,jloc,l,idir)
            end do
            counts(gu,gv)=1
            if(kplane==nxyz(axis) .and.&
               ps(igrid)%is_physical_boundary(2*axis)) then
              face_send(gu,gv)=mh_boundary_normal(igrid,i,jloc,l,axis)
              face_counts(gu,gv)=1
            end if
          end do
        end do
      case(3)
        lloc=ixMlo3+kplane-start(3)
        do j=ixMlo2,ixMhi2; jg=start(2)+j-ixMlo2
          do i=ixMlo1,ixMhi1; ig=start(1)+i-ixMlo1
            gu=ig; gv=jg
            do idir=1,3
              bsend(gu,gv,idir)=mh_total_b(igrid,i,j,lloc,idir)
              bpsend(gu,gv,idir)=bp%blocks(igrid)%b(i,j,lloc,idir)
            end do
            counts(gu,gv)=1
            if(kplane==nxyz(axis) .and.&
               ps(igrid)%is_physical_boundary(2*axis)) then
              face_send(gu,gv)=mh_boundary_normal(igrid,i,j,lloc,axis)
              face_counts(gu,gv)=1
            end if
          end do
        end do
      end select
    end do
  end subroutine mh_assemble_plane

  double precision function mh_total_b(igrid,i,j,k,idir)
    use mod_global_parameters
    integer, intent(in) :: igrid,i,j,k,idir
    mh_total_b=ps(igrid)%w(i,j,k,iw_mag(idir))
    if(B0field) mh_total_b=mh_total_b+ps(igrid)%B0(i,j,k,idir,0)
  end function mh_total_b

  double precision function mh_boundary_normal(igrid,i,j,k,axis)
    use mod_global_parameters
    integer, intent(in) :: igrid,i,j,k,axis
    integer :: ip,jp,kp

    if(stagger_grid) then
      mh_boundary_normal=ps(igrid)%ws(i,j,k,axis)
      if(B0field) mh_boundary_normal=mh_boundary_normal+&
         ps(igrid)%B0(i,j,k,axis,axis)
    else
      ip=i; jp=j; kp=k
      select case(axis)
      case(1); ip=i+1
      case(2); jp=j+1
      case(3); kp=k+1
      end select
      mh_boundary_normal=0.5d0*(ps(igrid)%w(i,j,k,iw_mag(axis))+&
         ps(igrid)%w(ip,jp,kp,iw_mag(axis)))
      if(B0field) mh_boundary_normal=mh_boundary_normal+0.5d0*(&
         ps(igrid)%B0(i,j,k,axis,0)+ps(igrid)%B0(ip,jp,kp,axis,0))
    end if
  end function mh_boundary_normal

  subroutine mh_build_top_b(bnormal,du,dv,uaxis,vaxis,axis,&
       orientation,btop)
    double precision, intent(in) :: bnormal(:,:),du,dv
    integer, intent(in) :: uaxis,vaxis,axis,orientation
    double precision, intent(out) :: btop(:,:,:)
    integer :: i,j,nu,nv

    nu=size(bnormal,1); nv=size(bnormal,2)
    btop=0.d0
    ! The same odd/even recurrence makes the centered transverse curl of b
    ! equal to the supplied top-face normal field at interior points.
    if(nu>=2) then
      do j=1,nv
        btop(2,j,vaxis)=0.25d0*dble(orientation)*du*&
           (bnormal(1,j)+bnormal(2,j))
        do i=3,nu
          btop(i,j,vaxis)=btop(i-2,j,vaxis)+&
             dble(orientation)*du*bnormal(i-1,j)
        end do
      end do
    end if
    if(nv>=2) then
      do i=1,nu
        btop(i,2,uaxis)=-0.25d0*dble(orientation)*dv*&
           (bnormal(i,1)+bnormal(i,2))
        do j=3,nv
          btop(i,j,uaxis)=btop(i,j-2,uaxis)-&
             dble(orientation)*dv*bnormal(i,j-1)
        end do
      end do
    end if
    btop(:,:,axis)=0.d0
  end subroutine mh_build_top_b

  subroutine mh_make_vector_potential(btop,q,uaxis,vaxis,axis,&
       orientation,avec)
    double precision, intent(in) :: btop(:,:,:),q(:,:,:)
    integer, intent(in) :: uaxis,vaxis,axis,orientation
    double precision, intent(out) :: avec(:,:,:)

    avec=btop
    avec(:,:,uaxis)=avec(:,:,uaxis)-dble(orientation)*q(:,:,vaxis)
    avec(:,:,vaxis)=avec(:,:,vaxis)+dble(orientation)*q(:,:,uaxis)
    avec(:,:,axis)=0.d0
  end subroutine mh_make_vector_potential

  subroutine mh_plane_integrals(avec,apvec,b,bp,bdiff,dvolume,values)
    double precision, intent(in) :: avec(:,:,:),apvec(:,:,:)
    double precision, intent(in) :: b(:,:,:),bp(:,:,:),dvolume
    double precision, intent(out) :: bdiff(:,:,:)
    double precision, intent(inout) :: values(11)

    bdiff=b-bp
    values(1)=values(1)+sum((avec+apvec)*bdiff)*dvolume
    values(2)=values(2)+sum((avec-apvec)*bdiff)*dvolume
    values(3)=values(3)+2.d0*sum(apvec*bdiff)*dvolume
    values(4)=values(4)+0.5d0*sum(b**2)*dvolume
    ! Use the full-field A*B scale for the near-zero Hm decision. A scale
    ! based only on B-Bp collapses together with a potential field and would
    ! make a roundoff-level Hm look like a valid denominator.
    values(5)=values(5)+sum((dsqrt(sum(avec**2,dim=3))+&
       dsqrt(sum(apvec**2,dim=3)))*(dsqrt(sum(b**2,dim=3))+&
       dsqrt(sum(bp**2,dim=3))))*dvolume
  end subroutine mh_plane_integrals

  subroutine mh_plane_quality(state,qcur,qcurp,bcur,uaxis,vaxis,axis,&
       orientation,du,dv,da,dvolume,values)
    double precision, intent(in) :: state(:,:,:),qcur(:,:,:),qcurp(:,:,:)
    double precision, intent(in) :: bcur(:,:,:),du,dv,da,dvolume
    integer, intent(in) :: uaxis,vaxis,axis,orientation
    double precision, intent(inout) :: values(11)
    double precision :: curl_a(3),curl_ap(3),btarget(3),bptarget(3),divb
    double precision :: s,du_av,dv_au,du_apv,dv_apu
    integer :: i,j

    s=dble(orientation)
    do j=2,size(state,2)-1
      do i=2,size(state,1)-1
        curl_a=0.d0
        curl_ap=0.d0
        curl_a(uaxis)=-(state(i,j,mh_q2_b-1+uaxis)-&
           qcur(i,j,uaxis))/(2.d0*da)
        curl_a(vaxis)=-(state(i,j,mh_q2_b-1+vaxis)-&
           qcur(i,j,vaxis))/(2.d0*da)
        curl_ap(uaxis)=-(state(i,j,mh_q2_bp-1+uaxis)-&
           qcurp(i,j,uaxis))/(2.d0*da)
        curl_ap(vaxis)=-(state(i,j,mh_q2_bp-1+vaxis)-&
           qcurp(i,j,vaxis))/(2.d0*da)
        du_av=(state(i+1,j,mh_btop-1+vaxis)+&
               s*state(i+1,j,mh_q1_b-1+uaxis)-&
               state(i-1,j,mh_btop-1+vaxis)-&
               s*state(i-1,j,mh_q1_b-1+uaxis))/(2.d0*du)
        dv_au=(state(i,j+1,mh_btop-1+uaxis)-&
               s*state(i,j+1,mh_q1_b-1+vaxis)-&
               state(i,j-1,mh_btop-1+uaxis)+&
               s*state(i,j-1,mh_q1_b-1+vaxis))/(2.d0*dv)
        du_apv=(state(i+1,j,mh_btop-1+vaxis)+&
                s*state(i+1,j,mh_q1_bp-1+uaxis)-&
                state(i-1,j,mh_btop-1+vaxis)-&
                s*state(i-1,j,mh_q1_bp-1+uaxis))/(2.d0*du)
        dv_apu=(state(i,j+1,mh_btop-1+uaxis)-&
                s*state(i,j+1,mh_q1_bp-1+vaxis)-&
                state(i,j-1,mh_btop-1+uaxis)+&
                s*state(i,j-1,mh_q1_bp-1+vaxis))/(2.d0*dv)
        curl_a(axis)=s*(du_av-dv_au)
        curl_ap(axis)=s*(du_apv-dv_apu)
        btarget=state(i,j,mh_b1:mh_b1+2)
        bptarget=state(i,j,mh_bp1:mh_bp1+2)
        values(6)=values(6)+sum((curl_a-btarget)**2)*dvolume
        values(7)=values(7)+sum(btarget**2)*dvolume
        values(8)=values(8)+sum((curl_ap-bptarget)**2)*dvolume
        values(9)=values(9)+sum(bptarget**2)*dvolume
        divb=(state(i+1,j,mh_b1-1+uaxis)-&
              state(i-1,j,mh_b1-1+uaxis))/(2.d0*du)+&
             (state(i,j+1,mh_b1-1+vaxis)-&
              state(i,j-1,mh_b1-1+vaxis))/(2.d0*dv)+&
             (state(i,j,mh_b2-1+axis)-bcur(i,j,axis))/(2.d0*da)
        values(10)=values(10)+divb**2*dvolume
        values(11)=values(11)+sum(btarget**2)*dvolume
      end do
    end do
  end subroutine mh_plane_quality

  subroutine mh_debug_pack_plane(bp,avec,apvec,bdiff,packed)
    double precision, intent(in) :: bp(:,:,:),avec(:,:,:),apvec(:,:,:)
    double precision, intent(in) :: bdiff(:,:,:)
    double precision, intent(out) :: packed(:,:)
    integer :: i,j,ip

    ip=0
    do j=1,size(bp,2)
      do i=1,size(bp,1)
        ip=ip+1
        packed(1:3,ip)=bp(i,j,1:3)
        packed(4:6,ip)=avec(i,j,1:3)
        packed(7:9,ip)=apvec(i,j,1:3)
        packed(10,ip)=sum((avec(i,j,1:3)+apvec(i,j,1:3))*&
             bdiff(i,j,1:3))
        packed(11,ip)=sum((avec(i,j,1:3)-apvec(i,j,1:3))*&
             bdiff(i,j,1:3))
        packed(12,ip)=2.d0*sum(apvec(i,j,1:3)*bdiff(i,j,1:3))
      end do
    end do
  end subroutine mh_debug_pack_plane

  subroutine mh_debug_vti_begin(writer,config,nxyz,cell_dx)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters, only: mype,base_filename,snapshotini,&
         xprobmin1,xprobmin2,xprobmin3,type_endian
    type(mh_debug_vti_state), intent(out) :: writer
    type(magnetic_helicity_config), intent(in) :: config
    integer, intent(in) :: nxyz(3)
    double precision, intent(in) :: cell_dx(3)
    character(len=mh_name_len) :: filename
    character(len=1024) :: line
    character(len=12) :: byte_order
    integer :: io_status,i
    integer(kind=8) :: ncell

    writer=mh_debug_vti_state()
    writer%enabled=config%write_debug_vti
    writer%nxyz=nxyz
    if(.not.writer%enabled) return
    if(len_trim(config%debug_vti_file)>0) then
      filename=trim(config%debug_vti_file)
    else
      write(filename,'(a,"_helicity_debug",i4.4,".vti")') &
           trim(base_filename),snapshotini
    end if

    ncell=int(nxyz(1),8)*int(nxyz(2),8)*int(nxyz(3),8)
    writer%nbytes(1:3)=3_8*ncell*8_8
    writer%nbytes(4:6)=ncell*8_8
    writer%offsets(1)=0_8
    do i=2,mh_debug_narray
      writer%offsets(i)=writer%offsets(i-1)+8_8+writer%nbytes(i-1)
    end do

    if(mype/=0) return
    byte_order=merge('LittleEndian','BigEndian   ',type_endian==1)
    open(newunit=writer%unit,file=trim(filename),status='replace',&
         action='write',access='stream',form='unformatted',iostat=io_status)
    if(io_status/=0) call mpistop('could not open magnetic-helicity debug VTI')
    call mh_stream_line(writer%unit,'<?xml version="1.0"?>')
    call mh_stream_line(writer%unit,'<VTKFile type="ImageData" version="1.0" '//&
         'byte_order="'//trim(byte_order)//'" header_type="UInt64">')
    write(line,'(a,6(i0,1x),a,3(es24.16,1x),a,3(es24.16,1x),a)') &
         '  <ImageData WholeExtent="',0,nxyz(1),0,nxyz(2),0,nxyz(3),&
         '" Origin="',xprobmin1,xprobmin2,xprobmin3,'" Spacing="',&
         cell_dx,'">'
    call mh_stream_line(writer%unit,trim(line))
    write(line,'(a,6(i0,1x),a)') '    <Piece Extent="',0,nxyz(1),&
         0,nxyz(2),0,nxyz(3),'">'
    call mh_stream_line(writer%unit,trim(line))
    call mh_stream_line(writer%unit,'      <CellData>')
    call mh_debug_vti_array_line(writer%unit,'Bp',3,writer%offsets(1))
    call mh_debug_vti_array_line(writer%unit,'A',3,writer%offsets(2))
    call mh_debug_vti_array_line(writer%unit,'Ap',3,writer%offsets(3))
    call mh_debug_vti_array_line(writer%unit,'h_m',1,writer%offsets(4))
    call mh_debug_vti_array_line(writer%unit,'h_J',1,writer%offsets(5))
    call mh_debug_vti_array_line(writer%unit,'h_PJ',1,writer%offsets(6))
    call mh_stream_line(writer%unit,'      </CellData>')
    call mh_stream_line(writer%unit,'    </Piece>')
    call mh_stream_line(writer%unit,'  </ImageData>')
    write(writer%unit) '<AppendedData encoding="raw">_'
    inquire(unit=writer%unit,pos=writer%data_pos)
    do i=1,mh_debug_narray
      write(writer%unit,pos=writer%data_pos+writer%offsets(i),&
           iostat=io_status) writer%nbytes(i)
      if(io_status/=0) call mpistop('could not initialize debug VTI payload')
    end do
  end subroutine mh_debug_vti_begin

  subroutine mh_debug_vti_array_line(vti_unit,name,ncomponent,offset)
    integer, intent(in) :: vti_unit,ncomponent
    character(len=*), intent(in) :: name
    integer(kind=8), intent(in) :: offset
    character(len=512) :: line

    write(line,'(a,a,a,i0,a,i0,a)') &
         '        <DataArray type="Float64" Name="',trim(name),&
         '" NumberOfComponents="',ncomponent,&
         '" format="appended" offset="',offset,'"/>'
    call mh_stream_line(vti_unit,trim(line))
  end subroutine mh_debug_vti_array_line

  subroutine mh_stream_line(vti_unit,line)
    integer, intent(in) :: vti_unit
    character(len=*), intent(in) :: line
    write(vti_unit) trim(line)//new_line('a')
  end subroutine mh_stream_line

  subroutine mh_debug_vti_write_plane(writer,owner,kplane,uaxis,vaxis,&
       axis,packed)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    type(mh_debug_vti_state), intent(inout) :: writer
    integer, intent(in) :: owner,kplane,uaxis,vaxis,axis
    double precision, intent(inout) :: packed(:,:)
    integer :: status(MPI_STATUS_SIZE),i,j,ip,ixyz(3),io_status
    integer(kind=8) :: global_index,position

    if(.not.writer%enabled) return
    if(owner/=0) then
      if(mype==owner) call MPI_SEND(packed,size(packed),&
           MPI_DOUBLE_PRECISION,0,mh_debug_tag,icomm,ierrmpi)
      if(mype==0) call MPI_RECV(packed,size(packed),MPI_DOUBLE_PRECISION,&
           owner,mh_debug_tag,icomm,status,ierrmpi)
    end if
    if(mype/=0) return

    ip=0
    do j=1,writer%nxyz(vaxis)
      do i=1,writer%nxyz(uaxis)
        ip=ip+1
        ixyz=1
        ixyz(axis)=kplane
        ixyz(uaxis)=i
        ixyz(vaxis)=j
        global_index=int(ixyz(1),8)+int(writer%nxyz(1),8)*(&
             int(ixyz(2)-1,8)+int(writer%nxyz(2),8)*int(ixyz(3)-1,8))
        position=writer%data_pos+writer%offsets(1)+8_8+&
             (global_index-1_8)*24_8
        write(writer%unit,pos=position,iostat=io_status) packed(1:3,ip)
        if(io_status/=0) call mpistop('could not write Bp debug VTI payload')
        position=writer%data_pos+writer%offsets(2)+8_8+&
             (global_index-1_8)*24_8
        write(writer%unit,pos=position,iostat=io_status) packed(4:6,ip)
        if(io_status/=0) call mpistop('could not write A debug VTI payload')
        position=writer%data_pos+writer%offsets(3)+8_8+&
             (global_index-1_8)*24_8
        write(writer%unit,pos=position,iostat=io_status) packed(7:9,ip)
        if(io_status/=0) call mpistop('could not write Ap debug VTI payload')
        position=writer%data_pos+writer%offsets(4)+8_8+&
             (global_index-1_8)*8_8
        write(writer%unit,pos=position,iostat=io_status) packed(10,ip)
        if(io_status/=0) call mpistop('could not write hm debug VTI payload')
        position=writer%data_pos+writer%offsets(5)+8_8+&
             (global_index-1_8)*8_8
        write(writer%unit,pos=position,iostat=io_status) packed(11,ip)
        if(io_status/=0) call mpistop('could not write hJ debug VTI payload')
        position=writer%data_pos+writer%offsets(6)+8_8+&
             (global_index-1_8)*8_8
        write(writer%unit,pos=position,iostat=io_status) packed(12,ip)
        if(io_status/=0) call mpistop('could not write hPJ debug VTI payload')
      end do
    end do
  end subroutine mh_debug_vti_write_plane

  subroutine mh_debug_vti_end(writer)
    use mod_global_parameters, only: mype,icomm,ierrmpi
    type(mh_debug_vti_state), intent(inout) :: writer
    integer(kind=8) :: footer_pos

    if(.not.writer%enabled) return
    if(mype==0) then
      footer_pos=writer%data_pos+writer%offsets(mh_debug_narray)+8_8+&
           writer%nbytes(mh_debug_narray)
      write(writer%unit,pos=footer_pos) new_line('a')//&
           '</AppendedData>'//new_line('a')//'</VTKFile>'//new_line('a')
      close(writer%unit)
    end if
    call MPI_BARRIER(icomm,ierrmpi)
  end subroutine mh_debug_vti_end
}

end module mod_magnetic_helicity
