!> Utilities for reading AMRVAC data-driven magnetic boundary frames.
!>
!> Python V1 writes each frame as a stream binary file with layout:
!> snapshot_time, nx, ny, dx, dy, Bx, By, Bz.
!>
!> The spacing dx/dy is stored in km and magnetic-field components are stored
!> in Gauss. The magnetic data are ordered as a Fortran array (nx, ny, 3).
module mod_data_driven_boundary
  implicit none

  type data_driven_boundary_series
    integer :: nframe=0,nx=0,ny=0
    double precision :: dx_km=0.d0,dy_km=0.d0
    character(len=1024) :: directory=''
    character(len=64) :: prefix='B_'
    double precision, allocatable :: times(:)
    integer :: cache_indices(2)=0
    ! Only the two frames bracketing the current observation time are cached.
    double precision, allocatable :: values(:,:,:,:)
  end type data_driven_boundary_series

contains

  subroutine read_data_driven_boundary_frame(filename,snapshot_time,nx,ny,dx_km,dy_km,bframe)
    character(len=*), intent(in) :: filename
    double precision, intent(out) :: snapshot_time,dx_km,dy_km
    integer, intent(out) :: nx,ny
    double precision, allocatable, intent(out) :: bframe(:,:,:)

    call read_data_driven_boundary_header(filename,snapshot_time,nx,ny,dx_km,dy_km)
    allocate(bframe(nx,ny,3))
    call read_data_driven_boundary_frame_into(filename,snapshot_time,nx,ny,dx_km,dy_km,bframe)
  end subroutine read_data_driven_boundary_frame

  subroutine read_data_driven_boundary_header(filename,snapshot_time,nx,ny,dx_km,dy_km)
    use mod_comm_lib, only: mpistop
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    character(len=*), intent(in) :: filename
    double precision, intent(out) :: snapshot_time,dx_km,dy_km
    integer, intent(out) :: nx,ny

    integer :: iu,ios
    integer(kind=8) :: file_size,expected_size
    logical :: exists

    inquire(file=filename,exist=exists)
    if(.not.exists) call mpistop('missing data-driven boundary frame')
    open(newunit=iu,file=trim(filename),status='old',access='stream', &
      form='unformatted',action='read',iostat=ios)
    if(ios/=0) call mpistop('cannot open data-driven boundary frame')
    read(iu,iostat=ios) snapshot_time
    if(ios==0) read(iu,iostat=ios) nx
    if(ios==0) read(iu,iostat=ios) ny
    if(ios==0) read(iu,iostat=ios) dx_km
    if(ios==0) read(iu,iostat=ios) dy_km
    inquire(unit=iu,size=file_size)
    close(iu)
    if(ios/=0) call mpistop('cannot read data-driven boundary header')
    if(nx<=0 .or. ny<=0) call mpistop('invalid nx/ny in data-driven boundary frame')
    if(.not.ieee_is_finite(snapshot_time) .or. .not.ieee_is_finite(dx_km) .or. &
       .not.ieee_is_finite(dy_km) .or. dx_km<=0.d0 .or. dy_km<=0.d0) &
      call mpistop('invalid time or spacing in data-driven boundary frame')
    expected_size=32_8+int(nx,kind=8)*int(ny,kind=8)*24_8
    if(file_size/=expected_size) &
      call mpistop('data-driven boundary frame size does not match its header')
  end subroutine read_data_driven_boundary_header

  subroutine read_data_driven_boundary_frame_into(filename,snapshot_time,nx,ny, &
      dx_km,dy_km,bframe)
    use mod_comm_lib, only: mpistop

    character(len=*), intent(in) :: filename
    double precision, intent(out) :: snapshot_time,dx_km,dy_km
    integer, intent(out) :: nx,ny
    double precision, intent(out) :: bframe(:,:,:)

    integer :: iu,ios

    call read_data_driven_boundary_header(filename,snapshot_time,nx,ny,dx_km,dy_km)
    if(size(bframe,1)/=nx .or. size(bframe,2)/=ny .or. size(bframe,3)/=3) &
      call mpistop('destination array does not match data-driven boundary frame')
    open(newunit=iu,file=trim(filename),status='old',access='stream', &
      form='unformatted',action='read',iostat=ios)
    if(ios/=0) call mpistop('cannot reopen data-driven boundary frame')
    read(iu,iostat=ios) snapshot_time
    if(ios==0) read(iu,iostat=ios) nx
    if(ios==0) read(iu,iostat=ios) ny
    if(ios==0) read(iu,iostat=ios) dx_km
    if(ios==0) read(iu,iostat=ios) dy_km
    if(ios==0) read(iu,iostat=ios) bframe
    close(iu)
    if(ios/=0) call mpistop('cannot read data-driven boundary payload')
  end subroutine read_data_driven_boundary_frame_into

  subroutine read_data_driven_boundary_series(directory,series,prefix,expected_nframe)
    use mod_comm_lib, only: mpistop

    character(len=*), intent(in) :: directory
    type(data_driven_boundary_series), intent(inout) :: series
    character(len=*), intent(in), optional :: prefix
    integer, intent(in), optional :: expected_nframe

    character(len=1024) :: filename
    character(len=64) :: frame_prefix
    double precision :: snapshot_time,dx_km,dy_km
    integer :: iframe,nx,ny,left_slot,right_slot
    logical :: exists

    frame_prefix='B_'
    if(present(prefix)) frame_prefix=trim(prefix)
    if(len_trim(directory)>len(series%directory)) &
      call mpistop('data-driven boundary directory path is too long')
    if(len_trim(frame_prefix)>len(series%prefix)) &
      call mpistop('data-driven boundary prefix is too long')
    series%directory=trim(directory)
    series%prefix=trim(frame_prefix)
    if(present(expected_nframe)) then
      if(expected_nframe<2) call mpistop('expected boundary frame count must be at least two')
      series%nframe=expected_nframe
      do iframe=1,series%nframe
        write(filename,'(a,"/",a,i4.4,".dat")') trim(directory),trim(frame_prefix),iframe
        inquire(file=trim(filename),exist=exists)
        if(.not.exists) call mpistop('missing frame in data-driven boundary series')
      end do
      write(filename,'(a,"/",a,i4.4,".dat")') trim(directory),trim(frame_prefix),series%nframe+1
      inquire(file=trim(filename),exist=exists)
      if(exists) call mpistop('boundary series contains more frames than declared')
    else
      iframe=1
      do
        write(filename,'(a,"/",a,i4.4,".dat")') trim(directory),trim(frame_prefix),iframe
        inquire(file=trim(filename),exist=exists)
        if(.not.exists) exit
        iframe=iframe+1
      end do
      series%nframe=iframe-1
    end if
    if(series%nframe<2) call mpistop('data-driven boundary series requires at least two frames')

    if(allocated(series%times)) deallocate(series%times)
    if(allocated(series%values)) deallocate(series%values)
    allocate(series%times(series%nframe))

    do iframe=1,series%nframe
      write(filename,'(a,"/",a,i4.4,".dat")') trim(directory),trim(frame_prefix),iframe
      call read_data_driven_boundary_header(trim(filename),snapshot_time,nx,ny,dx_km,dy_km)
      if(iframe==1) then
        series%nx=nx
        series%ny=ny
        series%dx_km=dx_km
        series%dy_km=dy_km
      else
        if(nx/=series%nx .or. ny/=series%ny) &
          call mpistop('inconsistent nx/ny in data-driven boundary series')
        if(abs(dx_km-series%dx_km)>1.d-10*max(1.d0,abs(series%dx_km))) &
          call mpistop('inconsistent dx in data-driven boundary series')
        if(abs(dy_km-series%dy_km)>1.d-10*max(1.d0,abs(series%dy_km))) &
          call mpistop('inconsistent dy in data-driven boundary series')
        if(snapshot_time<=series%times(iframe-1)) &
          call mpistop('data-driven boundary times must be strictly increasing')
      end if
      series%times(iframe)=snapshot_time
    end do
    if(abs(series%times(1))>1.d-10) &
      call mpistop('first data-driven boundary snapshot_time must be zero')
    allocate(series%values(series%nx,series%ny,3,2))
    series%cache_indices=0
    call ensure_data_driven_boundary_pair(series,1,left_slot,right_slot)
  end subroutine read_data_driven_boundary_series

  subroutine interpolate_data_driven_boundary(series,time_seconds,bframe,left_index,weight)
    type(data_driven_boundary_series), intent(inout) :: series
    double precision, intent(in) :: time_seconds
    double precision, intent(out) :: bframe(series%nx,series%ny,3)
    integer, intent(out), optional :: left_index
    double precision, intent(out), optional :: weight

    integer :: lo,left_slot,right_slot
    double precision :: alpha,denominator

    if(time_seconds<=series%times(1)) then
      lo=1
      alpha=0.d0
    else if(time_seconds>=series%times(series%nframe)) then
      lo=series%nframe-1
      alpha=1.d0
    else
      lo=1
      do while(lo<series%nframe-1 .and. time_seconds>=series%times(lo+1))
        lo=lo+1
      end do
      denominator=series%times(lo+1)-series%times(lo)
      alpha=(time_seconds-series%times(lo))/denominator
    end if
    call ensure_data_driven_boundary_pair(series,lo,left_slot,right_slot)
    bframe=(1.d0-alpha)*series%values(:,:,:,left_slot)+ &
      alpha*series%values(:,:,:,right_slot)
    if(present(left_index)) left_index=lo
    if(present(weight)) weight=alpha
  end subroutine interpolate_data_driven_boundary

  subroutine ensure_data_driven_boundary_pair(series,lo,left_slot,right_slot)
    type(data_driven_boundary_series), intent(inout) :: series
    integer, intent(in) :: lo
    integer, intent(out) :: left_slot,right_slot
    integer :: slot

    left_slot=0
    right_slot=0
    do slot=1,2
      if(series%cache_indices(slot)==lo) left_slot=slot
      if(series%cache_indices(slot)==lo+1) right_slot=slot
    end do
    if(left_slot>0 .and. right_slot>0) return
    if(left_slot>0) then
      right_slot=3-left_slot
      call load_data_driven_boundary_cache_slot(series,lo+1,right_slot)
    else if(right_slot>0) then
      left_slot=3-right_slot
      call load_data_driven_boundary_cache_slot(series,lo,left_slot)
    else
      left_slot=1
      right_slot=2
      call load_data_driven_boundary_cache_slot(series,lo,left_slot)
      call load_data_driven_boundary_cache_slot(series,lo+1,right_slot)
    end if
  end subroutine ensure_data_driven_boundary_pair

  subroutine load_data_driven_boundary_cache_slot(series,iframe,slot)
    use mod_comm_lib, only: mpistop

    type(data_driven_boundary_series), intent(inout) :: series
    integer, intent(in) :: iframe,slot
    character(len=1024) :: filename
    double precision :: snapshot_time,dx_km,dy_km
    double precision, parameter :: reltol=1.d-10
    integer :: nx,ny

    if(iframe<1 .or. iframe>series%nframe) &
      call mpistop('data-driven boundary cache index is out of range')
    write(filename,'(a,"/",a,i4.4,".dat")') trim(series%directory), &
      trim(series%prefix),iframe
    call read_data_driven_boundary_frame_into(trim(filename),snapshot_time,nx,ny, &
      dx_km,dy_km,series%values(:,:,:,slot))
    if(nx/=series%nx .or. ny/=series%ny .or. &
       abs(dx_km-series%dx_km)>reltol*max(1.d0,abs(series%dx_km)) .or. &
       abs(dy_km-series%dy_km)>reltol*max(1.d0,abs(series%dy_km)) .or. &
       abs(snapshot_time-series%times(iframe))>reltol*max(1.d0,abs(snapshot_time))) &
      call mpistop('data-driven boundary frame changed after header scan')
    series%cache_indices(slot)=iframe
  end subroutine load_data_driven_boundary_cache_slot

  subroutine interpolate_data_driven_boundary_scaled(series,time_code,unit_time_seconds, &
      driving_time_scale,bframe,left_index,weight)
    !> Interpolate using AMRVAC code time. ``driving_time_scale`` is the
    !> observational-time advance per simulated second; values greater than
    !> one therefore accelerate the observed boundary evolution.
    use mod_comm_lib, only: mpistop
    type(data_driven_boundary_series), intent(inout) :: series
    double precision, intent(in) :: time_code,unit_time_seconds,driving_time_scale
    double precision, intent(out) :: bframe(series%nx,series%ny,3)
    integer, intent(out), optional :: left_index
    double precision, intent(out), optional :: weight
    double precision :: observation_time

    if(unit_time_seconds<=0.d0) call mpistop('unit_time_seconds must be positive')
    if(driving_time_scale<=0.d0) call mpistop('driving_time_scale must be positive')
    observation_time=time_code*unit_time_seconds*driving_time_scale
    call interpolate_data_driven_boundary(series,observation_time,bframe,left_index,weight)
  end subroutine interpolate_data_driven_boundary_scaled

  subroutine destroy_data_driven_boundary_series(series)
    type(data_driven_boundary_series), intent(inout) :: series
    if(allocated(series%times)) deallocate(series%times)
    if(allocated(series%values)) deallocate(series%values)
    series%nframe=0
    series%nx=0
    series%ny=0
    series%directory=''
    series%prefix='B_'
    series%cache_indices=0
  end subroutine destroy_data_driven_boundary_series

end module mod_data_driven_boundary
