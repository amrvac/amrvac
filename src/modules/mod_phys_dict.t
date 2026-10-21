!> Quick module that mimics the behaviour of python dictionaries for passing attributes around 
!> without explicit dependencies

!> Jack Jenkins Dec 2025

module mod_phys_dict
    use mod_global_parameters

    implicit none
    private

    type phys_dict
        character(len=std_len) :: physics_type
        integer :: n = 0
        character(len=std_len), allocatable :: keys(:)
        logical, allocatable          :: values(:)
    end type phys_dict

    public :: phys_dict, set_dict_flag, get_dict_flag

contains

    subroutine set_dict_flag(d, key, val)
        type(phys_dict), intent(inout) :: d
        character(len=*), intent(in)   :: key
        logical, intent(in)            :: val

        character(len=std_len), allocatable :: new_keys(:)
        logical, allocatable          :: new_values(:)

        integer :: i

        !> If key exists, update it
        do i = 1, d%n
        if (trim(d%keys(i)) == trim(key)) then
            d%values(i) = val
            return
        end if
        end do

        !> Otherwise, append
        if (.not. allocated(d%keys)) then
            allocate(d%keys(1))
            allocate(d%values(1))
            d%n        = 1
            d%keys(1)  = trim(key)
            d%values(1)= val
        else
            !> Extend by 1

            allocate(new_keys(d%n+1))
            allocate(new_values(d%n+1))

            new_keys(1:d%n)   = d%keys
            new_values(1:d%n) = d%values

            new_keys(d%n+1)   = trim(key)
            new_values(d%n+1) = val

            call move_alloc(new_keys,   d%keys)
            call move_alloc(new_values, d%values)

            d%n = d%n + 1
        end if
    end subroutine set_dict_flag


    function get_dict_flag(d, key) result(val)
        type(phys_dict), intent(in)       :: d
        character(len=*), intent(in)      :: key
        logical                           :: val

        integer :: i

        do i = 1, d%n
        if (trim(d%keys(i)) == trim(key)) then
            val = d%values(i)
            return
        end if
        end do
        !> If key not found, return false
        val = .false.
    end function get_dict_flag

end module mod_phys_dict
