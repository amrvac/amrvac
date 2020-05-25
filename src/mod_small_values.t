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

  public :: small_values_error
  public :: small_values_average

contains

  subroutine small_values_error(w, x, ixI^L, ixO^L, w_flag, subname, smallw)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(in) :: smallw(1:nw)
    integer, intent(in)          :: w_flag(ixI^S)
    integer                      :: ix_bad(ndim), iw
    character(len=*), intent(in) :: subname

    ix_bad = maxloc(w_flag(ixO^S)) + [ ixOmin^D-1 ]

    if (.not.crash) then
      write(*,*) "Error: small value of ", trim(prim_wnames(maxval(w_flag(ixO^S)))), &
           " encountered when call ", subname
      write(*,*) "Iteration: ", it, " Time: ", global_time
      write(*,*) "Location: ", x({ix_bad(^D)}, :)
      write(*,*) "Cell number: ", ix_bad(:)
      do iw = 1, nw
         write(*, '(A20,A,E14.7)') trim(cons_wnames(iw)), ": ", &
              w({ix_bad(^D)}, iw)
      end do
      write(*,*) 'small value:',smallw(maxval(w_flag(ixO^S)))
      ! use erroneous arithmetic operation to crash the run
      if(trace_small_values) write(*,*) sqrt(smallw(1)-bigdouble)
      write(*,*) "Saving status at the previous time step"
      crash=.true.
    end if
  end subroutine small_values_error

  subroutine small_values_average(ixI^L, ixO^L, w, x, w_flag)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    integer, intent(in)             :: w_flag(ixI^S)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: iw, kxO^L, ix^D, i

    {do ix^DB= ixO^LIM^DB\}

    ! point with local failure identified by w_flag
    if (w_flag(ix^D) /= 0) then
      ! verify in cube with border width small_values_daverage the presence of
      ! cells where all went ok
      do i = 1, max(small_values_daverage, 1)
        {kxOmin^D= max(ix^D-i, ixOmin^D);
        kxOmax^D= min(ix^D+i, ixOmax^D);\}

        ! in case cells are fine within smaller cube than 
        ! the userset small_values_daverage: use that smaller cube
        if (any(w_flag(kxO^S) == 0)) exit
      end do

      if (any(w_flag(kxO^S) == 0)) then
        ! within surrounding cube, cells without problem were found

        ! faulty cells are corrected by averaging here
        ! only average those which were ok and replace faulty cells
        do iw = 1, nw
          if (small_values_fix_iw(iw)) then
            w(ix^D, iw) = sum(w(kxO^S, iw), w_flag(kxO^S) == 0)&
                 / count(w_flag(kxO^S) == 0)
          end if
        end do
      else
        write(*,*) "no cells without error were found in cube of size", & 
             small_values_daverage
        write(*,*) "at location:", x(ix^D, 1:ndim)
        write(*,*) "at index:", ix^D
        write(*,*) "w_flag(ix^D):", w_flag(ix^D)
        write(*,*) "Saving status at the previous time step"
        crash=.true.
      end if
    end if
    {enddo^D&\}

  end subroutine small_values_average

end module mod_small_values
