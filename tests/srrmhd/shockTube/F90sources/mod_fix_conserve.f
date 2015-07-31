module mod_fix_conserve
   use mod_indices, only: ngridshi
   implicit none

   type fluxalloc
      double precision, dimension(:,:), pointer:: flux
   end type fluxalloc
   type(fluxalloc), dimension(2,1,ngridshi), save :: pflux

   integer, save :: nrecv, nsend
   double precision, dimension(:), allocatable, save :: recvbuffer
   integer, dimension(1), save :: isize

end module mod_fix_conserve
