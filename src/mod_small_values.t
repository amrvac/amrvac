!> Module for handling problematic values in simulations, such as negative
!> pressures
module mod_small_values

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

contains

  subroutine small_values_error(wprim, x, ixI^L, ixO^L, w_flag, subname)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: wprim(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    logical, intent(in)          :: w_flag(ixI^S,1:nw)
    character(len=*), intent(in) :: subname
    integer                      :: iw,iiw,ix^D

    if (.not.crash) then
      do iw=1,nw
       {do ix^DB= ixO^LIM^DB\}
          if(w_flag(ix^D,iw)) then
            write(*,*) "Error: small value of ", trim(prim_wnames(iw)),wprim(ix^D,iw),&
                 " encountered when call ", subname
            write(*,*) "Iteration: ", it, " Time: ", global_time
            write(*,*) "Location: ", x(ix^D,:)
            write(*,*) "Cell number: ", ix^D
            do iiw=1,nw
              write(*,*) trim(prim_wnames(iiw)),": ",wprim(ix^D,iiw)
            end do
            ! use erroneous arithmetic operation to crash the run
            if(trace_small_values) write(*,*) sqrt(wprim(ix^D,iw)-bigdouble)
            write(*,*) "Saving status at the previous time step"
            crash=.true.
          end if
       {enddo^D&\}
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
            {kxOmin^D= max(ix^D-i, ixOmin^D);
            kxOmax^D= min(ix^D+i, ixOmax^D);\}
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
            write(*,*) "no cells without error were found in cube of size", & 
                 small_values_daverage
            write(*,*) "at location:", x(ix^D, 1:ndim)
            write(*,*) "at index:", ix^D
            write(*,*) "w numer:", iw
            write(*,*) "Saving status at the previous time step"
            crash=.true.
          end if
        end if
     {enddo^D&\}
    end do

  end subroutine small_values_average

end module mod_small_values
