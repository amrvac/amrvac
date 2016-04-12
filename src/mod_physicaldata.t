module mod_physicaldata
  use mod_indices, only: ngridshi
   implicit none
   save

   type walloc
      double precision, dimension(:^D&,:), pointer :: w
   end type walloc

{^NOONED
   type walloc_sub
      double precision, dimension(:^DE&,:), pointer :: w
   end type walloc_sub
}
{^IFONED
   type walloc_sub
      double precision, dimension(:), pointer :: w
   end type walloc_sub
}
   type(walloc), dimension(ngridshi) :: pw, pwold, pw1, pw2, pw3, pw4, pwres
   type(walloc), dimension(ngridshi) :: pwCoarse, pwCoCo
   type(walloc), dimension(ngridshi) :: pwio
   type(walloc), dimension(ngridshi), target :: pB0_cell,  pB0_face^D
{#IFDEF BOUNDARYDRIVER
   type(walloc), dimension(2*^ND), target :: pB0bc_cell,  pB0bc_face^D
}
   type(walloc), pointer :: myB0_cell, myB0_face^D, myB0
   type(walloc_sub), dimension(ngridshi) :: pw_sub
{^IFONED
   double precision, dimension(:), allocatable :: collapsedData
}
{^NOONED
   double precision, dimension(:^DE&,:), allocatable :: collapsedData
}

   type xalloc
      double precision, dimension(:^D&,:), pointer :: x
   end type xalloc

{^NOONED
   type xalloc_sub
      double precision, dimension(:^DE&,:), pointer :: x
   end type xalloc_sub
}
{^IFONED
   type xalloc_sub
      double precision, dimension(:), pointer :: x
   end type xalloc_sub
}
   type(xalloc), dimension(ngridshi), target :: px, pxCoarse
   type(xalloc_sub), dimension(ngridshi) :: px_sub
   !!! type(xalloc), pointer :: myx

   type geoalloc
      double precision, dimension(:^D&), pointer :: dvolume
      double precision, dimension(:^D&), pointer :: surfaceC^D,surface^D
      double precision, dimension(:^D&,:), pointer :: dx,x
   end type geoalloc

   type(geoalloc), dimension(ngridshi), target :: pgeo, pgeoCoarse, pgeoCoCo
   type(geoalloc), pointer                     :: mygeo
!$OMP THREADPRIVATE(myB0_cell,myB0_face^D,myB0,mygeo)   
end module mod_physicaldata
