!> Module with basic grid data structures
module mod_forest
   implicit none

   !> Pointer to a tree_node
   type tree_node_ptr
      type(tree_node), pointer :: node
   end type tree_node_ptr

   !> The data structure that contains information about a tree node/grid block
   type tree_node
      integer :: ig^D !< Spatial indices of the grid block
      integer :: level !< Refinement level
      integer :: igrid !< Index of grid on processor's grid-array
      integer :: ipe   !< On which processor the grid is stored
      logical :: leaf !< Is the grid a leaf (no further refinement)
      logical :: active
      type(tree_node_ptr) :: parent !< Pointer to parent grid
      type(tree_node_ptr) :: child(2^D&) !< Pointer to children
      type(tree_node_ptr) :: neighbor(2,^ND) !< Pointer to neighbors
      type(tree_node_ptr) :: next !< Next node at refinement level (linked list)
      type(tree_node_ptr) :: prev !< Previous node at refinement level (linked list)
   end type tree_node

   !> Pointers to the coarse grid
   type(tree_node_ptr), dimension(:^D&), allocatable, save :: tree_root

   !> Array to go from an [igrid, ipe] index to a node pointer
   type(tree_node_ptr), dimension(:,:), allocatable, save :: igrid_to_node

   !> The head pointer of the linked list per refinement level
   type(tree_node_ptr), dimension(:), allocatable, save :: level_head

   !> The tail pointer of the linked list per refinement level
   type(tree_node_ptr), dimension(:), allocatable, save :: level_tail

   !> Array to go from a Morton number to an igrid and processor index. Sfc(1:3,
   !> MN) contains [igrid, ipe, active], where MN is a morton number and active
   !> is 0 or 1.
   integer, dimension(:,:), allocatable, save :: sfc

   !> Space filling curve for level 1 grid. sfc_iglevel1(^D, MN) gives ig^D (the
   !> spatial index of the grid) 
   integer, dimension(:,:), allocatable, save :: sfc_iglevel1

   !> iglevel1_sfc(ig^D) gives the Morton number for grid ig^D
   integer, dimension(:^D&), allocatable, save :: iglevel1_sfc

   !> Go from a Morton number to an igrid index (for a single processor)
   integer, dimension(:), allocatable, save :: sfc_to_igrid

   !> Go from a grid index to Morton number (for a single processor)
   integer, dimension(:), allocatable, save :: igrid_to_sfc

   !> Space filling curve used for physical boundary blocks
   integer, dimension(:), allocatable, save :: sfc_phybound

   !> First Morton number per processor
   integer, dimension(:), allocatable, save :: Morton_start

   !> Last Morton number per processor
   integer, dimension(:), allocatable, save :: Morton_stop

   integer, dimension(:), allocatable, save :: Morton_sub_start, Morton_sub_stop

   !> AMR flags and grids-in-use identifier per processor (igrid,ipe) 
   logical, dimension(:,:), allocatable, save :: coarsen, refine, buffer, igrid_inuse

   !> Number of parent blocks
   integer, save :: nparents

   !> Number of leaf block
   integer, save :: nleafs

   integer :: nleafs_active, nglev1{#IFDEF EVOLVINGBOUNDARY , nphyboundblock}

   !> How many leaves are present per refinement level
   integer, dimension(:), allocatable, save :: nleafs_level

end module mod_forest
