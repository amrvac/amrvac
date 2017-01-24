module mod_small_values

  implicit none
  private

  integer :: daverage = 1

  character(len=20), public :: small_values_method = "error"

  public :: small_values_init
  public :: small_values_error
  public :: small_values_average

contains

  subroutine small_values_init()
    use mod_global_parameters

    character(len=20) :: method

    call CFG_add_get(cfg, "small_value%method", small_values_method, &
         "How to handle small values. Options: error, replace, average")

    call CFG_add_get(cfg, "small_value%daverage", daverage, &
         "Distance (in cells) used for averaging small values")

  end subroutine small_values_init

  subroutine small_values_error(w, x, ixI^L, ixO^L, w_flag)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    integer, intent(in)          :: w_flag(ixI^S)
    integer                      :: ix_bad(ndim)

    ix_bad = maxloc(w_flag(ixO^S)) + [ ixOmin^D-1 ]

    write(*, *) "Error: small value encountered"
    write(*, *) "Iteration: ", it, " Time: ", t
    write(*, *) "Location: ", x({ix_bad(^D)}, :)
    write(*, *) "w(1:nw): ", w({ix_bad(^D)}, 1:nw)

    call mpistop("Error: small value encountered")
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

       ! verify in cube with border width daverage the presence
       ! of cells where all went ok
       do i = 1, daverage
          {kxOmin^D= max(ix^D-i, ixOmin^D);
          kxOmax^D= min(ix^D+i, ixOmax^D);\}

          ! in case cells are fine within smaller cube than 
          ! the userset daverage: use that smaller cube
          if (any(w_flag(kxO^S) == 0)) exit
       end do

       if (any(w_flag(kxO^S) == 0)) then
          ! within surrounding cube, cells without problem were found

          ! faulty cells are corrected by averaging here
          ! only average those which were ok and replace faulty cells
          do iw = 1, nw
             w(ix^D, iw) = sum(w(kxO^S, iw), w_flag(kxO^S) == 0)&
                  / count(w_flag(kxO^S) == 0)
          end do
       else
          ! no cells without error were found in cube of size daverage
          ! --> point of no recovery
          print*,"Getaux error:", w_flag(ix^D),"ix^D=", ix^D
          !print*,"New ","rho=", w(ix^D, rho_),"m=", &
          !        {^C&w(ix^D, m^C_)},"e=", w(ix^D, e_)
          !print*,"position  ", px(saveigrid)%x(ix^D, 1:ndim)
          if (w_flag(ix^D)<0) then
             call mpistop("-small_values_average from smallvalues-----")
          else
             call mpistop("-small_values_average from primitive or getpthermal--")
          end if
       end if
    end if
    {enddo^D&\}

  end subroutine small_values_average

end module mod_small_values
