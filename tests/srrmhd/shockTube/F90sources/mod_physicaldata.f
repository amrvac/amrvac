module mod_physicaldata
  use mod_indices, only: ngridshi
   implicit none
   save

   type walloc
      double precision, dimension(:,:), pointer :: w
   end type walloc



   type walloc_sub
      double precision, dimension(:), pointer :: w
   end type walloc_sub

   type(walloc), dimension(ngridshi) :: pw, pwold, pw1, pw2, pw3, pw4, pwres
   type(walloc), dimension(ngridshi) :: pwCoarse, pwCoCo
   type(walloc), dimension(ngridshi) :: pwio
   type(walloc), dimension(ngridshi), target :: pB0_cell,  pB0_face1
   type(walloc), pointer :: myB0_cell, myB0_face1, myB0
   type(walloc_sub), dimension(ngridshi) :: pw_sub

   double precision, dimension(:), allocatable :: collapsedData



   type xalloc
      double precision, dimension(:,:), pointer :: x
   end type xalloc



   type xalloc_sub
      double precision, dimension(:), pointer :: x
   end type xalloc_sub

   type(xalloc), dimension(ngridshi), target :: px, pxCoarse
   type(xalloc_sub), dimension(ngridshi) :: px_sub
   !!! type(xalloc), pointer :: myx

   type geoalloc
      double precision, dimension(:), pointer :: dvolume
      double precision, dimension(:), pointer :: surfaceC1,surface1
      double precision, dimension(:,:), pointer :: dx
   end type geoalloc

   type(geoalloc), dimension(ngridshi), target :: pgeo, pgeoCoarse, pgeoCoCo
   type(geoalloc), pointer                     :: mygeo
   
end module mod_physicaldata
