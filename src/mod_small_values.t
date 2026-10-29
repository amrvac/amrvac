!> Module for handling problematic values in simulations, such as negative
!> pressures
module mod_small_values

  use iso_fortran_env, only: int64

  implicit none
  private

  !> How to handle small values
  character(len=20), public :: small_values_method = "error"

  !> Average over this many cells in each direction
  integer, public :: small_values_daverage = 1

  !> trace small values in the source file using traceback flag of compiler
  logical, public :: trace_small_values=.false.

  !> Whether to apply small value fixes to certain variables
  logical, public, allocatable :: small_values_fix_iw(:)

  public :: small_values_error
  public :: small_values_average

  ! Per-process state; access only through report_repair_warning's critical region.
  type, public :: repair_warning_state
    integer(int64) :: events = 0_int64, next_report = 10_int64
    integer :: details = 0
  end type repair_warning_state
  public :: report_repair_warning

contains

  ! Counts repair events, not distinct cells. No communication in this local path.
  subroutine report_repair_warning(state,label,events,edges,iw,xfirst,before,after)
    use mod_global_parameters, only: mype, it
    type(repair_warning_state), intent(inout) :: state
    character(len=*), intent(in) :: label
    integer(int64), intent(in) :: events, edges
    integer, intent(in) :: iw
    double precision, intent(in) :: xfirst(:), before, after
    logical :: milestone

    if(events<=0_int64) return
    !$omp critical(amrvac_repair_warning)
    ! Saturate rather than wrap on exceptionally long runs.
    state%events=state%events+min(events,huge(state%events)-state%events)
    milestone=.false.
    if(state%next_report>0_int64) then
      milestone=state%events>=state%next_report
      do while(state%events>=state%next_report)
        if(state%next_report>huge(state%next_report)/10_int64) then
          state%next_report=0_int64
          exit
        end if
        state%next_report=state%next_report*10_int64
      end do
    end if
    if(state%details<3) then
      state%details=state%details+1
      write(*,*) 'WARNING ',label,': rank=',mype,' it=',it, &
           ' events=',events,' edges=',edges,' cumulative_events=',state%events
      write(*,*) '  first: variable=',iw,' x=',xfirst,' before=',before,' after=',after
      if(state%details==3) write(*,*) '  Further detailed reports suppressed; cumulative milestones only.'
    else if(milestone) then
      write(*,*) 'WARNING ',label,': rank=',mype,' it=',it,' cumulative_events=',state%events
    end if
    !$omp end critical(amrvac_repair_warning)
  end subroutine report_repair_warning

  subroutine small_values_error(wprim, x, ixI^L, ixO^L, w_flag, subname)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: wprim(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    logical, intent(in)          :: w_flag(ixI^S,1:nw)
    character(len=*), intent(in) :: subname

    double precision :: x^D
    integer :: iw,iiw,ix^D

    if (.not.crash) then
      do iw=1,nw
       {do ix^DB= ixO^LIM^DB\}
          if(w_flag(ix^D,iw)) then
            write(*,*) "Error: small value of ", trim(prim_wnames(iw)),wprim(ix^D,iw),&
                 " encountered when call ", subname
            write(*,*) "Iteration: ", it, " Time: ", global_time, "Processor: ",mype
            write(*,*) "Location: ", x(ix^D,:)
            select case (coordinate)
              case (Cartesian,Cartesian_stretched,Cartesian_expansion)
              case (cylindrical)
                {^IFONED
                x1=x(ix^D,1)}
                {^IFTWOD
                select case (phi_)
                case (2)
                   x1=x(ix^D,1)*cos(x(ix^D,2))
                   x2=x(ix^D,1)*sin(x(ix^D,2))
                case default
                   x1=x(ix^D,1)
                   x2=x(ix^D,2)
                end select}
                {^IFTHREED
                x1=x(ix^D,1)*cos(x(ix^D,phi_))
                x2=x(ix^D,1)*sin(x(ix^D,phi_))
                x3=x(ix^D,z_)}
                write(*,*) "Cartesian location: ", x^D
              case (spherical)
                x1=x(ix^D,1){^NOONED*sin(x(ix^D,2))}{^IFTHREED*cos(x(ix^D,3))}
                {^IFTWOD
                x2=x(ix^D,1)*cos(x(ix^D,2))}
                {^IFTHREED
                x2=x(ix^D,1)*sin(x(ix^D,2))*sin(x(ix^D,3))
                x3=x(ix^D,1)*cos(x(ix^D,2))}
                write(*,*) "Cartesian location: ", x^D
              case default
                write(*,*) "No converter for coordinate=",coordinate
            end select
            write(*,*) "Cell number: ", ix^D
            do iiw=1,nw
              write(*,*) trim(prim_wnames(iiw)),": ",wprim(ix^D,iiw)
            end do
            ! use erroneous arithmetic operation to crash the run
            if(trace_small_values) write(*,*) sqrt(wprim(ix^D,iw)-bigdouble)
            write(*,*) "Saving status at the previous time step"
            crash=.true.
          end if
       {end do^D&\}
      end do
    end if
  end subroutine small_values_error

  subroutine small_values_average(ixI^L, ixO^L, w, x, w_flag, windex)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    logical, intent(in)             :: w_flag(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer, optional, intent(in)   :: windex
    integer                         :: iw, kxO^L, ix^D, i, nwstart, nwend

    type(repair_warning_state), save :: warning_state
    integer(int64) :: fallback_events
    integer :: first_iw
    double precision :: first_x(ndim), first_before, first_after

    fallback_events=0_int64
    if(present(windex)) then
      nwstart=windex
      nwend=windex
    else
      nwstart=1
      nwend=nw
    end if

    do iw = nwstart, nwend
     {do ix^DB= ixO^LIM^DB\}
      ! point with local failure identified by w_flag
        if (w_flag(ix^D,iw)) then
          ! verify in cube with border width small_values_daverage the presence of
          ! cells where all went ok
          do i = 1, max(small_values_daverage, 1)
            {kxOmin^D= max(ix^D-i, ixImin^D);
            kxOmax^D= min(ix^D+i, ixImax^D);\}
            ! in case cells are fine within smaller cube than 
            ! the userset small_values_daverage: use that smaller cube
            if(any(w_flag(kxO^S,iw) .eqv. .false.)) exit
          end do

          if(any(w_flag(kxO^S,iw) .eqv. .false.)) then
            ! within surrounding cube, cells without problem were found

            ! faulty cells are corrected by averaging here
            ! only average those which were ok and replace faulty cells
            if(small_values_fix_iw(iw)) then
              w(ix^D, iw) = sum(w(kxO^S, iw), w_flag(kxO^S,iw) .eqv. .false.)&
                   / count(w_flag(kxO^S,iw) .eqv. .false.)
            end if
          else
            fallback_events=fallback_events+1_int64
            if(fallback_events==1_int64) then
              first_iw=iw
              first_x=x(ix^D,1:ndim)
              first_before=w(ix^D,iw)
            end if
            if(iw==iw_e) w(ix^D, iw)=small_pressure
            if(iw==iw_rho) w(ix^D, iw)=small_density
            if(fallback_events==1_int64) first_after=w(ix^D,iw)
          end if
        end if
     {enddo^D&\}
    end do
    ! Other variables can remain unchanged in this existing fallback path.
    if(fallback_events>0_int64) call report_repair_warning(warning_state, &
         'small_values_average: no valid neighbours',fallback_events,0_int64, &
         first_iw,first_x,first_before,first_after)
  end subroutine small_values_average

end module mod_small_values
