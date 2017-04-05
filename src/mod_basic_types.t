!> Module with basic data types used in amrvac
module mod_basic_types
  implicit none
  public

    !> Size (in bytes) of single precision real
  integer, parameter :: size_real = 4

  !> Size (in bytes) of double precision real
  integer, parameter :: size_double = 8

  !> Size (in bytes) of default integer
  integer, parameter :: size_int = 4

  !> Size (in bytes) of default logical
  integer, parameter :: size_logical = 4

    !> Default length for strings
  integer, parameter :: std_len = 131

  !> Default length for names (of e.g. variables)
  integer, parameter :: name_len = 16

end module mod_basic_types
