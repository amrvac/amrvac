module mod_physicaldata
   implicit none
   save

   type walloc
      !> ID of a grid block
      integer :: igrid=-1
      !> location of w-array, 0: cell center, ^D : cell interface in dimension ^D
      integer :: iwpos=0
      !> Is w in primitive state or not
      logical :: w_is_primitive=.false.
      !> Variables, normally center
      double precision, dimension(:^D&,:), allocatable :: w
      !> Background fixed components of w at cell center
      double precision, dimension(:^D&,:), allocatable :: w0
      !> Background fixed components of w at cell interface
      double precision, dimension(:^D&,:,:), allocatable :: w0face
      !> Cell-center positions
      double precision, dimension(:^D&,:), allocatable :: x
      !> Cell-center positions for non-Cartesian coordinates
      double precision, dimension(:^D&,:), allocatable :: xg
      !> Cell sizes
      double precision, dimension(:^D&,:), allocatable :: dx
      !> Volumes of a cell
      double precision, dimension(:^D&), allocatable :: dvolume
      !> Areas of cell-center surfaces
      double precision, dimension(:^D&), allocatable :: surface^D
      !> Areas of cell-face surfaces
      double precision, dimension(:^D&), allocatable :: surfaceC^D
   end type walloc

{^NOONED
   type walloc_sub
      !> ID of a grid block
      integer :: igrid=-1
      !> location of w-array, 0: cell center, ^D : cell interface in dimension ^D
      integer :: iwpos=0
      !> Is w in primitive state or not
      logical :: w_is_primitive=.false.
      !> Variables, normally center
      double precision, dimension(:^DE&,:), allocatable :: w
      !> Background fixed components of w at cell center
      double precision, dimension(:^DE&,:), allocatable :: w0
      !> Background fixed components of w at cell interface
      double precision, dimension(:^DE&,:,:), allocatable :: w0face
      !> Cell-center positions
      double precision, dimension(:^DE&,:), allocatable :: x
      !> Cell-center positions for non-Cartesian coordinates
      double precision, dimension(:^DE&,:), allocatable :: xg
      !> Cell sizes
      double precision, dimension(:^DE&,:), allocatable :: dx
      !> Volumes of a cell
      double precision, dimension(:^DE&), allocatable :: dvolume
      !> Areas of cell-center surfaces 
      double precision, dimension(:^DE&), allocatable :: surface^DE
      !> Areas of cell-face surfaces
      double precision, dimension(:^DE&), allocatable :: surfaceC^DE
   end type walloc_sub
}
{^IFONED
   type walloc_sub
      !> ID of a grid block
      integer :: igrid=-1
      !> location of w-array, 0: cell center, ^D : cell interface in dimension ^D
      integer :: iwpos=0
      !> Is w in primitive state or not
      logical :: w_is_primitive=.false.
      !> Variables, normally center
      double precision, dimension(:), allocatable :: w
      !> Background fixed components of w at cell center
      double precision, dimension(:), allocatable :: w0
      !> Background fixed components of w at cell interface
      double precision, dimension(:), allocatable :: w0face
      !> Cell-center positions
      double precision, dimension(:), allocatable :: x
   end type walloc_sub
}
   ! array of physical variables
   type(walloc), dimension(:), allocatable :: pw, pwold, pw1, pw2, pw3, pw4, &
                                              pwCoarse, pwio
   type(walloc), dimension(:), allocatable, target :: pB0_cell,  pB0_face^D
{#IFDEF BOUNDARYDRIVER
   type(walloc), dimension(2*^ND), target :: pB0bc_cell,  pB0bc_face^D
}
   type(walloc), pointer :: myB0_cell, myB0_face^D, myB0
   type(walloc_sub), dimension(:), allocatable :: pw_sub
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
   ! array of spacial coordinates
   type(xalloc), dimension(:), allocatable, target :: px, pxCoarse
   type(xalloc_sub), dimension(:), allocatable :: px_sub

   type geoalloc
      double precision, dimension(:^D&), pointer :: dvolume
      double precision, dimension(:^D&), pointer :: surfaceC^D,surface^D
      double precision, dimension(:^D&,:), pointer :: dx,x
   end type geoalloc
   ! array of geometric information
   type(geoalloc), dimension(:), allocatable, target :: pgeo, pgeoCoarse
   type(geoalloc), pointer                     :: mygeo
!$OMP THREADPRIVATE(myB0_cell,myB0_face^D,myB0,mygeo)   
end module mod_physicaldata
