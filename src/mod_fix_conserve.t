module mod_fix_conserve
   use mod_indices, only: ngridshi
   implicit none

   type fluxalloc
      double precision, dimension(:^D&,:), pointer:: flux
   end type fluxalloc
   type(fluxalloc), dimension(2,^ND,ngridshi), save :: pflux

   integer, save :: nrecv, nsend
   double precision, dimension(:), allocatable, save :: recvbuffer
   integer, dimension(^ND), save :: isize

end module mod_fix_conserve
