!> This module contains global integers and indices
!> @todo Why are there logicals/reals here?
!> @todo Clean this list up
module mod_indices

  implicit none

  !> The maximum number of levels in the grid refinement
  !> @todo Don't use a fixed upper bound
  integer, parameter :: nlevelshi = 13

end module mod_indices
