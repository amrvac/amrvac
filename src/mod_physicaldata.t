module mod_physicaldata
   implicit none
   save

   type walloc
      !> ID of a grid block
      integer :: igrid=-1
      !> location of w-array, 0: cell center, ^D : cell interface in dimension ^D
      integer :: iwpos=0
      !> location of w0-array, 0: cell center, ^D : cell interface in dimension ^D
      integer :: iw0=0
      !> Is w in primitive state or not
      logical :: w_is_primitive=.false.
      !> Variables, normally center
      double precision, dimension(:^D&,:), allocatable :: w
      !> Variables, normally center, temperary state for multi-step scheme
      double precision, dimension(:^D&,:), allocatable :: w1
      !> Variables, normally center, temperary state for multi-step scheme
      double precision, dimension(:^D&,:), allocatable :: w2
      !> Variables, normally center, temperary state for multi-step scheme
      double precision, dimension(:^D&,:), allocatable :: w3
      !> Variables, normally center, temperary state for multi-step scheme
      double precision, dimension(:^D&,:), allocatable :: w4
      !> Variables, normally center, pointer of reference state
      double precision, dimension(:^D&,:), pointer :: wa => Null()
      !> Variables, normally center, pointer of updated state
      double precision, dimension(:^D&,:), pointer :: wb => Null()
      !> Variables, normally center, initial state at the beginning of iteration
      double precision, dimension(:^D&,:), allocatable :: wold
      !> Variables, normally center, for visualization data
      double precision, dimension(:^D&,:), allocatable :: wio
      !> Variables, normally center, one level coarser representative
      double precision, dimension(:^D&,:), allocatable :: wcoarse
      !> Background fixed components of w at cell center and cell interface
      double precision, dimension(:^D&,:,:), allocatable :: w0
      !> Cell-center positions
      double precision, dimension(:^D&,:), allocatable :: x
      !> Cell-center positions, one level coarser representative
      double precision, dimension(:^D&,:), allocatable :: xcoarse
      !> Cell-center positions for non-Cartesian coordinates
      double precision, dimension(:^D&,:), allocatable :: xg
      !> Cell sizes
      double precision, dimension(:^D&,:), allocatable :: dx
      !> Volumes of a cell
      double precision, dimension(:^D&), allocatable :: dvolume
      !> Volumes of a cell, one level coarser representative
      double precision, dimension(:^D&), allocatable :: dvolumecoarse
      !> Volumes of a cell, pointer
      double precision, dimension(:^D&), pointer :: dvolumep => Null()
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
      !> Variables, normally center, one level coarser representative
      double precision, dimension(:^DE&,:), allocatable :: wcoarse
      !> Cell-center positions
      double precision, dimension(:^DE&,:), allocatable :: x
      !> Cell-center positions, one level coarser representative
      double precision, dimension(:^DE&,:), allocatable :: xcoarse
      !> Cell-center positions for non-Cartesian coordinates
      double precision, dimension(:^DE&,:), allocatable :: xg
      !> Cell sizes
      double precision, dimension(:^DE&,:), allocatable :: dx
      !> Volumes of a cell
      double precision, dimension(:^DE&), allocatable :: dvolume
      !> Volumes of a cell, one level coarser representative
      double precision, dimension(:^DE&), allocatable :: dvolumecoarse
      !> Volumes of a cell, pointer
      double precision, dimension(:^DE&), pointer :: dvolumep => Null()
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
      !> Variables, normally center, one level coarser representative
      double precision, dimension(:), allocatable :: wcoarse
      !> Cell-center positions
      double precision, dimension(:), allocatable :: x
      !> Cell-center positions, one level coarser representative
      double precision, dimension(:), allocatable :: xcoarse
   end type walloc_sub
}
   !> Block pointer for using current block
   type(walloc), pointer :: block

   !> array of physical blocks
   type(walloc), dimension(:), allocatable, target :: pw

   !> array of physical blocks in reduced dimension
   type(walloc_sub), dimension(:), allocatable, target :: pw_sub

{^IFONED
   double precision, dimension(:), allocatable :: collapsedData
}
{^NOONED
   double precision, dimension(:^DE&,:), allocatable :: collapsedData
}

end module mod_physicaldata
