module mod_physicaldata
   implicit none
   save

   type mesh_t
      !> ID of a grid block
      integer :: igrid=-1
      !> index range of block array in cell centers
      integer :: ixG^L
      !> index range of block array in cell faces
      integer :: ixGs^L
      !> level of AMR
      integer :: level
      !> If it face a physical boundary
      logical, dimension(:), allocatable :: is_physical_boundary(:)

      !> Cell-center positions
      double precision, dimension(:^D&,:),     allocatable :: x
      !> Barycenter positions
      double precision, dimension(:^D&,:),     allocatable :: xbar
      !> Cell sizes in coordinate units
      double precision, dimension(:^D&,:),     allocatable :: dx
      !> Cell sizes at cell center in length unit
      double precision, dimension(:^D&,:),     allocatable :: ds
      !> Cell sizes at cell face in length unit
      double precision, dimension(:^D&,:),     allocatable :: dsC
      !> line integral \sqrt{gamma} dx at cell edge in length unit 
      !fixme: seems no need
      double precision, dimension(:^D&,:),     allocatable :: dlE
      !> Volumes of a cell
      double precision, dimension(:^D&),       allocatable :: dvolume
      !> Areas of cell-center surfaces
      double precision, dimension(:^D&,:),     allocatable :: surface
      !> Areas of cell-face surfaces
      double precision, dimension(:^D&,:),     allocatable :: surfaceC

      !> Christoffel symbols of a cell 
      double precision, dimension(:^D&,:), allocatable :: christoffel
   end type mesh_t

   type metric_t
      !> metric quantities
      ! this array is used in MPI communication,
      ! which includes: alpha, betai, vecX, h_ij
      ! fixme: remove vecX
      double precision, dimension(:^D&,:), allocatable :: vars
   end type metric_t

   type state
      logical :: is_prim
      !> Variables, cell center values
      double precision, dimension(:^D&,:), allocatable :: w
      !> Variables, cell face values
      double precision, dimension(:^D&,:), allocatable :: ws
      !> extra variables do not need ghost cell and equation flux
      double precision, dimension(:^D&,:), allocatable :: wextra
      !> metric varialbes
      type(metric_t), pointer :: metric => Null()
      !> mesh varialbes
      type(mesh_t),   pointer :: mesh   => Null()
   end type state

   ! fixme: this is for output, move this to io module
{^NOONED
   type state_sub
      !> ID of a grid block
      integer :: igrid=-1
      !> Variables, normally center
      double precision, dimension(:^DE&,:), allocatable :: w 
      !> Variables for the cornerpositions on the slice 
      double precision, dimension(:^DE&,:), allocatable :: wC
      !> Variables, normally center, one level coarser representative
      double precision, dimension(:^DE&,:), allocatable :: wcoarse
      !> Cell-center positions
      double precision, dimension(:^DE&,:), allocatable :: x
      !> Corner positions on the slice
      double precision, dimension(:^DE&,:), allocatable :: xC
      !> Cell-center positions, one level coarser representative
      double precision, dimension(:^DE&,:), allocatable :: xcoarse
      !> Cell sizes
      double precision, dimension(:^DE&,:), allocatable :: dx
      !> Cell sizes, one level coarser
      double precision, dimension(:^D&,:), allocatable :: dxcoarse
      !> Cell sizes in length unit
      double precision, dimension(:^D&,:), allocatable :: ds
      !> Volumes of a cell
      double precision, dimension(:^DE&), allocatable :: dvolume
      !> Volumes of a cell, one level coarser representative
      double precision, dimension(:^DE&), allocatable :: dvolumecoarse
      !> Areas of cell-center surfaces 
      double precision, dimension(:^DE&,:), allocatable :: surface
      !> Areas of cell-face surfaces
      double precision, dimension(:^DE&,:), allocatable :: surfaceC
   end type state_sub
}
{^IFONED
   type state_sub
      !> ID of a grid block
      integer :: igrid=-1
      !> Variables, normally center
      double precision, dimension(:), allocatable :: w
      !> Variables for the cornerpositions on the slice 
      double precision, dimension(:), allocatable :: wC
      !> Variables, normally center, one level coarser representative
      double precision, dimension(:), allocatable :: wcoarse
      !> Cell-center positions
      double precision, dimension(:), allocatable :: x
      !> Corner positions on the slice
      double precision, dimension(:), allocatable :: xC
      !> Cell-center positions, one level coarser representative
      double precision, dimension(:), allocatable :: xcoarse
   end type state_sub
}

   !> array of meshes for all blocks on my processor
   type(mesh_t), dimension(:), allocatable, target :: mesh
   !> array of physical blocks, one level coarser representative
   type(mesh_t), dimension(:), allocatable, target :: mesh_co

   !> array of meshes for all blocks on my processor
   type(metric_t), dimension(:), allocatable, target :: metric
   !> array of physical blocks, one level coarser representative
   type(metric_t), dimension(:), allocatable, target :: metric_co

   !> buffer for pole boundary
   type(state) :: pole_buf
   !> array of physical states for all blocks on my processor
   type(state), dimension(:), allocatable, target :: ps
   !> array of physical states, temp 1 for multi-step time integrator
   type(state), dimension(:), allocatable, target :: ps1
   !> array of physical states, temp 2 for multi-step time integrator
   type(state), dimension(:), allocatable, target :: ps2
   !> array of physical states, temp 3 for multi-step time integrator
   type(state), dimension(:), allocatable, target :: ps3
   !> array of physical states, temp 4 for multi-step time integrator
   type(state), dimension(:), allocatable, target :: ps4
   !> array of physical states, at the beginning of each iteration
   type(state), dimension(:), allocatable, target :: pso
   !> array of physical blocks, one level coarser representative
   type(state), dimension(:), allocatable, target :: psc

   !> array of physical blocks in reduced dimension
   type(state_sub), dimension(:), allocatable, target :: ps_sub

   {^IFONED
   double precision, dimension(:), allocatable :: collapsedData
   }
   {^NOONED
   double precision, dimension(:^DE&,:), allocatable :: collapsedData
   }
  !> velocities store for constrained transport
  type ct_velocity
    double precision, dimension(:^D&,:), allocatable :: vnorm,cbarmin,cbarmax
    double precision, dimension(:^D&,:,:), allocatable :: vbarC,vbarLC,vbarRC
    double precision, dimension(:^D&,:,:), allocatable :: betaC
    double precision, dimension(:^D&,:), allocatable :: EhatLC, EhatRC
  end type ct_velocity

end module mod_physicaldata
