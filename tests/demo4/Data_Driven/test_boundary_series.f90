program test_boundary_series
  use mpi
  use mod_global_parameters, only: icomm,mype,npe,ierrmpi
  use mod_data_driven_boundary
  implicit none

  type(data_driven_boundary_series) :: series
  double precision, allocatable :: frame(:,:,:)
  character(len=1024) :: directory
  integer :: left
  double precision :: weight

  call MPI_INIT(ierrmpi)
  icomm=MPI_COMM_WORLD
  call MPI_COMM_RANK(icomm,mype,ierrmpi)
  call MPI_COMM_SIZE(icomm,npe,ierrmpi)
  call get_command_argument(1,directory)
  if(len_trim(directory)==0) error stop 'usage: test_boundary_series DIRECTORY'

  call read_data_driven_boundary_series(trim(directory),series,expected_nframe=3)
  allocate(frame(series%nx,series%ny,3))
  call interpolate_data_driven_boundary(series,0.d0,frame,left,weight)
  if(left/=1 .or. abs(weight)>1.d-12) error stop 'wrong first-frame endpoint'
  if(maxval(abs(frame))>1.d-12) error stop 'wrong field at first-frame endpoint'
  call interpolate_data_driven_boundary(series,5.d0,frame,left,weight)
  if(series%nframe/=3) error stop 'wrong frame count'
  if(size(series%values,4)/=2) error stop 'series cache must contain exactly two frames'
  if(left/=1 .or. abs(weight-0.5d0)>1.d-12) error stop 'wrong interpolation weight'
  if(maxval(abs(frame-5.d0))>1.d-12) error stop 'wrong interpolated field'
  call interpolate_data_driven_boundary(series,10.d0,frame,left,weight)
  if(left/=2 .or. abs(weight)>1.d-12) error stop 'wrong interior-frame endpoint'
  if(maxval(abs(frame-10.d0))>1.d-12) error stop 'wrong field at interior-frame endpoint'
  call interpolate_data_driven_boundary(series,15.d0,frame,left,weight)
  if(left/=2 .or. abs(weight-0.5d0)>1.d-12) error stop 'wrong sliding-cache weight'
  if(maxval(abs(frame-15.d0))>1.d-12) error stop 'wrong sliding-cache field'
  call interpolate_data_driven_boundary(series,20.d0,frame,left,weight)
  if(left/=2 .or. abs(weight-1.d0)>1.d-12) error stop 'wrong last-frame endpoint'
  if(maxval(abs(frame-20.d0))>1.d-12) error stop 'wrong field at last-frame endpoint'
  if(mype==0) write(*,*) 'boundary series interpolation: PASS'
  call destroy_data_driven_boundary_series(series)
  deallocate(frame)
  call MPI_FINALIZE(ierrmpi)
end program test_boundary_series
