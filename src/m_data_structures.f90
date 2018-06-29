#include "cpp_macros.h"
module m_data_structures

  implicit none
  public

  !> Type of reals
  integer, parameter :: dp = kind(0.0d0)

  !> Type for 64-bit integers
  integer, parameter :: i8 = selected_int_kind(18)

  !> Indicates a standard Laplacian
  integer, parameter :: mg_laplacian = 1

  !> Indicates a variable-coefficient Laplacian
  integer, parameter :: mg_vlaplacian = 2

  !> Indicates a constant-coefficient Helmholtz equation
  integer, parameter :: mg_helmholtz = 3

  !> Indicates a variable-coefficient Helmholtz equation
  integer, parameter :: mg_vhelmholtz = 4

  integer, parameter :: mg_cartesian   = 1 !< Cartesian coordinate system
  integer, parameter :: mg_cylindrical = 2 !< Cylindrical coordinate system
  integer, parameter :: mg_spherical   = 3 !< Spherical coordinate system

  integer, parameter :: mg_smoother_gs     = 1
  integer, parameter :: mg_smoother_gsrb   = 2
  integer, parameter :: mg_smoother_jacobi = 3

  !> Problem dimension
  integer, parameter :: mg_ndim = NDIM

  !> Number of predefined multigrid variables
  integer, parameter :: mg_num_vars = 4
  !> Maximum number of variables
  integer, parameter :: mg_max_num_vars = 10
  !> Index of solution
  integer, parameter :: mg_iphi = 1
  !> Index of right-hand side
  integer, parameter :: mg_irhs = 2
  !> Index of previous solution (used for correction)
  integer, parameter :: mg_iold = 3
  !> Index of residual
  integer, parameter :: mg_ires = 4

  !> Index of the variable coefficient (at cell centers)
  integer, parameter :: mg_iveps = 5

  !> Minimum allowed grid level
  integer, parameter :: mg_lvl_lo = -20
  !> Maximum allowed grid level
  integer, parameter :: mg_lvl_hi = 20

  !> Value to indicate a Dirichlet boundary condition
  integer, parameter :: mg_bc_dirichlet = -10

  !> Value to indicate a Neumann boundary condition
  integer, parameter :: mg_bc_neumann = -11

  !> Value to indicate a continuous boundary condition
  integer, parameter :: mg_bc_continuous = -12

  !> Special value that indicates there is no box
  integer, parameter :: mg_no_box = 0
  !> Special value that indicates there is a physical boundary
  integer, parameter :: mg_physical_boundary = -1

  !> Maximum number of timers to use
  integer, parameter :: mg_max_timers = 20

#if NDIM == 2
  ! Numbering of children (same location as **corners**)
  integer, parameter :: mg_num_children = 4

  ! Index offset for each child
  integer, parameter :: mg_child_dix(2, 4) = reshape([0,0,1,0,0,1,1,1], [2,4])
  ! Reverse child index in each direction
  integer, parameter :: mg_child_rev(4, 2) = reshape([2,1,4,3,3,4,1,2], [4,2])
  ! Children adjacent to a neighbor
  integer, parameter :: mg_child_adj_nb(2, 4) = reshape([1,3,2,4,1,2,3,4], [2,4])
  ! Which children have a low index per dimension
  logical, parameter :: mg_child_low(2, 4) = reshape([.true., .true., &
       .false., .true., .true., .false., .false., .false.], [2, 4])

  ! Neighbor topology information
  integer, parameter :: mg_num_neighbors = 4
  integer, parameter :: mg_neighb_lowx = 1
  integer, parameter :: mg_neighb_highx = 2
  integer, parameter :: mg_neighb_lowy = 3
  integer, parameter :: mg_neighb_highy = 4

  ! Index offsets of neighbors
  integer, parameter :: mg_neighb_dix(2, 4) = reshape([-1,0,1,0,0,-1,0,1], [2,4])
  ! Which neighbors have a lower index
  logical, parameter :: mg_neighb_low(4) = [.true., .false., .true., .false.]
  ! Opposite of nb_low, but now as -1,1 integers
  integer, parameter :: mg_neighb_high_pm(4) = [-1, 1, -1, 1]

  ! Reverse neighbors
  integer, parameter :: mg_neighb_rev(4) = [2, 1, 4, 3]
  ! Direction (dimension) for a neighbor
  integer, parameter :: mg_neighb_dim(4) = [1, 1, 2, 2]
#elif NDIM == 3
    ! Numbering of children (same location as **corners**)
  integer, parameter :: mg_num_children = 8

  ! Index offset for each child
  integer, parameter :: mg_child_dix(3, 8) = reshape( &
       [0,0,0, 1,0,0, 0,1,0, 1,1,0, &
       0,0,1, 1,0,1, 0,1,1, 1,1,1], [3,8])
  ! Reverse child index in each direction
  integer, parameter :: mg_child_rev(8, 3) = reshape( &
       [2,1,4,3,6,5,8,7, 3,4,1,2,7,8,5,6, 5,6,7,8,1,2,3,4], [8,3])
  ! Children adjacent to a neighbor
  integer, parameter :: mg_child_adj_nb(4, 6) = reshape( &
       [1,3,5,7, 2,4,6,8, 1,2,5,6, 3,4,7,8, 1,2,3,4, 5,6,7,8], [4,6])
  ! Which children have a low index per dimension
  logical, parameter :: mg_child_low(3, 8) = reshape([ &
       .true., .true., .true., .false., .true., .true., &
       .true., .false., .true., .false., .false., .true., &
       .true., .true., .false., .false., .true., .false., &
       .true., .false., .false., .false., .false., .false.], [3, 8])

  ! Neighbor topology information
  integer, parameter :: mg_num_neighbors = 6
  integer, parameter :: mg_neighb_lowx = 1
  integer, parameter :: mg_neighb_highx = 2
  integer, parameter :: mg_neighb_lowy = 3
  integer, parameter :: mg_neighb_highy = 4
  integer, parameter :: mg_neighb_lowz = 5
  integer, parameter :: mg_neighb_highz = 6
  ! Index offsets of neighbors
  integer, parameter :: mg_neighb_dix(3, 6) = reshape( &
       [-1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1], [3,6])
  ! Which neighbors have a lower index
  logical, parameter :: mg_neighb_low(6) = &
       [.true., .false., .true., .false., .true., .false.]
  ! Opposite of nb_low, but now as -1,1 integers
  integer, parameter :: mg_neighb_high_pm(6) = [-1, 1, -1, 1, -1, 1]
  ! Reverse neighbors
  integer, parameter :: mg_neighb_rev(6) = [2, 1, 4, 3, 6, 5]
  ! Direction (dimension) for a neighbor
  integer, parameter :: mg_neighb_dim(6) = [1, 1, 2, 2, 3, 3]
#endif

  !> Lists of blocks per refinement level
  type mg_lvl_t
     integer, allocatable :: leaves(:)
     integer, allocatable :: parents(:)
     integer, allocatable :: ref_bnds(:)
     integer, allocatable :: ids(:)
     integer, allocatable :: my_leaves(:)
     integer, allocatable :: my_parents(:)
     integer, allocatable :: my_ref_bnds(:)
     integer, allocatable :: my_ids(:)
  end type mg_lvl_t

  !> Box data structure
  type mg_box_t
     integer  :: rank              !< Which process owns this box
     integer  :: id                !< Box id (index in boxes(:) array)
     integer  :: lvl               !< Refinement level
     integer  :: ix(NDIM)          !< Spatial index
     integer  :: parent            !< Id of parent
     integer  :: children(2**NDIM) !< Ids of children
     integer  :: neighbors(2*NDIM) !< Ids of neighbors
     real(dp) :: r_min(NDIM)       !< Minimum coordinate
     real(dp) :: dr(NDIM)          !< Grid spacing
     !> Cell-centered data
     real(dp), allocatable :: cc(DTIMES(:), :)
  end type mg_box_t

  !> Buffer type (one is used for each pair of communicating processes)
  type mg_buf_t
     integer               :: i_send !< Index in send array
     integer               :: i_recv
     integer               :: i_ix
     integer, allocatable  :: ix(:) !< Will be used to sort the data
     real(dp), allocatable :: send(:)
     real(dp), allocatable :: recv(:)
  end type mg_buf_t

  type mg_comm_t
     integer, allocatable :: n_send(:, :)
     integer, allocatable :: n_recv(:, :)
  end type mg_comm_t

  type mg_bc_t
     integer  :: bc_type  = mg_bc_dirichlet !< Type of boundary condition
     real(dp) :: bc_value = 0.0_dp       !< Value (for e.g. Dirichlet or Neumann)
     !> To set user-defined boundary conditions (overrides bc(:))
     procedure(mg_subr_bc), pointer, nopass :: boundary_cond  => null()
     !> To set a user-defined refinement boundary method
     procedure(mg_subr_rb), pointer, nopass :: refinement_bnd => null()
  end type mg_bc_t

  type mg_timer_t
     character(len=20) :: name
     real(dp)          :: t = 0.0_dp
     real(dp)          :: t0
  end type mg_timer_t

  type mg_t
     !> Whether the multigrid tree structure has been created
     logical                  :: tree_created     = .false.
     !> Whether storage has been allocated for boxes
     logical                  :: is_allocated     = .false.
     !> Number of extra cell-centered variable (e.g., for coefficients)
     integer                  :: n_extra_vars      = 0
     !> MPI communicator
     integer                  :: comm             = -1
     !> Number of MPI tasks
     integer                  :: n_cpu            = -1
     !> MPI rank of this task
     integer                  :: my_rank          = -1
     !> Size of boxes in cells, be equal for all dimensions
     integer                  :: box_size         = -1
     !> Highest grid level in the tree
     integer                  :: highest_lvl      = -1
     !> Lowest grid level in the tree
     integer                  :: lowest_lvl       = -1
     !> First normal level of the quadtree/octree, at coarser levels parents
     !> have only one child
     integer                  :: first_normal_lvl = -1
     !> Total number of boxes in the tree (among all processes)
     integer                  :: n_boxes          = 0
     !> Size of boxes per level (differs for coarsest levels)
     integer                  :: box_size_lvl(mg_lvl_lo:mg_lvl_hi)
     !> Size of domain per level (if uniformly refined)
     integer                  :: domain_size_lvl(NDIM, mg_lvl_lo:mg_lvl_hi)
     !> Grid spacing per level
     real(dp)                 :: dr(NDIM, mg_lvl_lo:mg_lvl_hi)
     !> Minimum coordinates
     real(dp)                 :: r_min(NDIM)
     !> List of all levels
     type(mg_lvl_t)              :: lvls(mg_lvl_lo:mg_lvl_hi)
     !> Array with all boxes in the tree. Only boxes owned by this task are
     !> actually allocated
     type(mg_box_t), allocatable :: boxes(:)
     !> Buffer for communication with each other task
     type(mg_buf_t), allocatable :: buf(:)

     !> Communication info for restriction
     type(mg_comm_t)             :: comm_restrict
     !> Communication info for prolongation
     type(mg_comm_t)             :: comm_prolong
     !> Communication info for ghost cell filling
     type(mg_comm_t)             :: comm_ghostcell

     !> Whether boundary condition data has been stored for mg solution
     logical :: phi_bc_data_stored = .false.

     !> To store pre-defined boundary conditions per direction per variable
     type(mg_bc_t) :: bc(mg_num_neighbors, mg_max_num_vars)

     !> Type of operator
     integer :: operator_type = mg_laplacian
     !> Type of grid geometry
     integer :: geometry_type = mg_cartesian

     !> Whether the mean has to be subtracted from the multigrid solution
     logical  :: subtract_mean       = .false.
     !> Type of multigrid smoother
     integer  :: smoother_type       = mg_smoother_gs
     !> Number of substeps for the smoother (for GSRB this is 2)
     integer  :: n_smoother_substeps = 1
     !> Number of cycles when doing downwards relaxation
     integer  :: n_cycle_down        = 2
     !> Number of cycles when doing upwards relaxation
     integer  :: n_cycle_up          = 2
     !> Maximum number of cycles on the coarse grid
     integer  :: max_coarse_cycles   = 1000
     integer  :: coarsest_grid(NDIM) = 2
     !> Stop coarse grid when max. residual is smaller than this
     real(dp) :: residual_coarse_abs = 1e-8_dp
     !> Stop coarse grid when residual has been reduced by this factor
     real(dp) :: residual_coarse_rel = 1e-8_dp

     !> Multigrid operator (e.g., Laplacian)
     procedure(mg_box_op), pointer, nopass   :: box_op => null()

     !> Multigrid smoother
     procedure(mg_box_gsrb), pointer, nopass :: box_smoother => null()

     !> Multigrid prolongation method
     procedure(mg_box_prolong), pointer, nopass :: box_prolong => null()

     !> Number of timers
     integer       :: n_timers = 0
     !> Values for the timers
     type(mg_timer_t) :: timers(mg_max_timers)
  end type mg_t

  interface
     !> To fill ghost cells near physical boundaries
     subroutine mg_subr_bc(box, nc, iv, nb, bc_type, bc)
       import
       type(mg_box_t), intent(in) :: box
       integer, intent(in)     :: nc
       integer, intent(in)     :: iv      !< Index of variable
       integer, intent(in)     :: nb      !< Direction
       integer, intent(out)    :: bc_type !< Type of b.c.
       !> Boundary values
#if NDIM == 2
       real(dp), intent(out)   :: bc(nc)
#elif NDIM == 3
       real(dp), intent(out)   :: bc(nc, nc)
#endif
     end subroutine mg_subr_bc

     !> To fill ghost cells near refinement boundaries
     subroutine mg_subr_rb(box, nc, iv, nb, cgc)
       import
       type(mg_box_t), intent(inout) :: box
       integer, intent(in)        :: nc
       integer, intent(in)        :: iv !< Index of variable
       integer, intent(in)        :: nb !< Direction
       !> Coarse data
#if NDIM == 2
       real(dp), intent(in)       :: cgc(nc)
#elif NDIM == 3
       real(dp), intent(in)       :: cgc(nc, nc)
#endif
     end subroutine mg_subr_rb

     !> Subroutine that performs A * cc(..., i_in) = cc(..., i_out)
     subroutine mg_box_op(mg, id, nc, i_out)
       import
       type(mg_t), intent(inout) :: mg
       integer, intent(in)       :: id
       integer, intent(in)       :: nc
       integer, intent(in)       :: i_out
     end subroutine mg_box_op

     !> Subroutine that performs Gauss-Seidel relaxation
     subroutine mg_box_gsrb(mg, id, nc, redblack_cntr)
       import
       type(mg_t), intent(inout) :: mg
       integer, intent(in)       :: id
       integer, intent(in)       :: nc
       integer, intent(in)       :: redblack_cntr
     end subroutine mg_box_gsrb

     !> Subroutine that performs prolongation to a single child
     subroutine mg_box_prolong(mg, p_id, dix, nc, iv, fine)
       import
       type(mg_t), intent(inout) :: mg
       integer, intent(in)       :: p_id             !< Id of parent
       integer, intent(in)       :: dix(NDIM)        !< Offset of child in parent grid
       integer, intent(in)       :: nc               !< Child grid size
       integer, intent(in)       :: iv               !< Prolong from this variable
       real(dp), intent(out)     :: fine(DTIMES(nc)) !< Prolonged values
     end subroutine mg_box_prolong
  end interface

contains

  !> Return .true. if a box has children
  elemental logical function mg_has_children(box)
    type(mg_box_t), intent(in) :: box

    ! Boxes are either fully refined or not, so we only need to check one of the
    ! children
    mg_has_children = (box%children(1) /= mg_no_box)
  end function mg_has_children

  !> Compute the 'child index' for a box with spatial index ix. With 'child
  !> index' we mean the index in the children(:) array of its parent.
  integer function mg_ix_to_ichild(ix)
    integer, intent(in) :: ix(NDIM) !< Spatial index of the box
    ! The index can range from 1 (all ix odd) and 2**$D (all ix even)
#if NDIM == 2
    mg_ix_to_ichild = 4 - 2 * iand(ix(2), 1) - iand(ix(1), 1)
#elif NDIM == 3
    mg_ix_to_ichild = 8 - 4 * iand(ix(3), 1) - &
         2 * iand(ix(2), 1) - iand(ix(1), 1)
#endif
  end function mg_ix_to_ichild

  !> Get the offset of a box with respect to its parent (e.g. in 2d, there can
  !> be a child at offset 0,0, one at n_cell/2,0, one at 0,n_cell/2 and one at
  !> n_cell/2, n_cell/2)
  pure function mg_get_child_offset(mg, id) result(ix_offset)
    type(mg_t), intent(in) :: mg
    integer, intent(in)    :: id
    integer                :: ix_offset(NDIM)

    if (mg%boxes(id)%lvl <= mg%first_normal_lvl) then
       ix_offset(:) = 0
    else
       ix_offset = iand(mg%boxes(id)%ix-1, 1) * &
            ishft(mg%box_size, -1) ! * n_cell / 2
    end if
  end function mg_get_child_offset

  pure function mg_highest_uniform_lvl(mg) result(lvl)
    type(mg_t), intent(in) :: mg
    integer                :: lvl

    do lvl = mg%first_normal_lvl, mg%highest_lvl-1
       ! Exit if a grid is partially refined
       if (size(mg%lvls(lvl)%leaves) /= 0 .and. &
           size(mg%lvls(lvl)%parents) /= 0) exit
    end do
    ! If the loop did not exit, we get lvl equals mg%highest_lvl
  end function mg_highest_uniform_lvl

  !> Determine total number of unknowns (on leaves)
  function mg_number_of_unknowns(mg) result(n_unknowns)
    type(mg_t), intent(in) :: mg
    integer                :: lvl
    integer(i8)            :: n_unknowns

    n_unknowns = 0
    do lvl = mg%first_normal_lvl, mg%highest_lvl
       n_unknowns = n_unknowns + size(mg%lvls(lvl)%leaves)
    end do
    n_unknowns = n_unknowns * int(mg%box_size**3, i8)
  end function mg_number_of_unknowns

  !> Get coordinates at the face of a box
  subroutine mg_get_face_coords(box, nb, nc, x)
    type(mg_box_t), intent(in) :: box
    integer, intent(in)        :: nb
    integer, intent(in)        :: nc
#if NDIM == 2
    real(dp), intent(out)      :: x(nc, 2)
    integer                    :: i
#elif NDIM == 3
    real(dp), intent(out)      :: x(nc, nc, 3)
    integer                    :: i, j
#endif
    integer                    :: nb_dim, ixs(NDIM-1)
    real(dp)                   :: rmin(NDIM)

    ! Determine directions perpendicular to neighbor
    nb_dim = mg_neighb_dim(nb)
    ixs                     = [(i, i = 1, NDIM-1)]
    ixs(nb_dim:) = ixs(nb_dim:) + 1

    rmin = box%r_min
    if (.not. mg_neighb_low(nb)) then
       rmin(nb_dim) = rmin(nb_dim) + box%dr(nb_dim) * nc
    end if

#if NDIM == 2
    do i = 1, nc
       x(i, :) = rmin
       x(i, ixs(1)) = x(i, ixs(1)) + (i-0.5d0) * box%dr(ixs(1))
    end do
#elif NDIM == 3
    do j = 1, nc
       do i = 1, nc
          x(i, j, :) = rmin
          x(i, j, ixs) = x(i, j, ixs) + ([i, j] - 0.5d0) * box%dr(ixs)
       end do
    end do
#endif
  end subroutine mg_get_face_coords

  integer function mg_add_timer(mg, name)
    type(mg_t), intent(inout) :: mg
    character(len=*), intent(in) :: name

    mg%n_timers                  = mg%n_timers + 1
    mg_add_timer                 = mg%n_timers
    mg%timers(mg_add_timer)%name = name
  end function mg_add_timer

  subroutine mg_timer_start(timer)
    use mpi
    type(mg_timer_t), intent(inout) :: timer
    timer%t0 = mpi_wtime()
  end subroutine mg_timer_start

  subroutine mg_timer_end(timer)
    use mpi
    type(mg_timer_t), intent(inout) :: timer
    timer%t = timer%t + mpi_wtime() - timer%t0
  end subroutine mg_timer_end

  subroutine mg_timers_show(mg)
    use mpi
    type(mg_t), intent(in) :: mg
    integer                :: n, ierr
    real(dp)               :: tmin(mg%n_timers)
    real(dp)               :: tmax(mg%n_timers)

    call mpi_reduce(mg%timers(1:mg%n_timers)%t, tmin, mg%n_timers, &
         mpi_double, mpi_min, 0, mg%comm, ierr)
    call mpi_reduce(mg%timers(1:mg%n_timers)%t, tmax, mg%n_timers, &
         mpi_double, mpi_max, 0, mg%comm, ierr)

    if (mg%my_rank == 0) then
       write(*, "(A20,2A16)") "name                ", "min(s)", "max(s)"
       do n = 1, mg%n_timers
          write(*, "(A20,2F16.6)") mg%timers(n)%name, &
               tmin(n), tmax(n)
       end do
    end if
  end subroutine mg_timers_show

end module m_data_structures
