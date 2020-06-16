

! Single module generated from the octree-mg sources.
! This file can be easier to include in existing projects.
!
! Notes:
! 1. The module name is here extended by _1d, _2d or _3d
! 2. The free space Poisson solver is not included here.
! 3. It is best to make changes in the original repository at
!    https://github.com/jannisteunissen/octree-mg
! 4. This file can be generated as follows:
!    cd octree-mg/single_module
!    ./to_single_module.sh
!
! The modules can be compiled with:
! mpif90 -c m_octree_mg_1d.f90 [other options]
! mpif90 -c m_octree_mg_2d.f90 [other options]
! mpif90 -c m_octree_mg_3d.f90 [other options]
! mpif90 -c m_octree_mg.f90 -cpp -DNDIM=<1,2,3> [other options]


module m_octree_mg_1d
  use mpi
  implicit none
  private

  !! File ../src/m_data_structures.f90

  !> Type of reals
  integer, parameter, public :: dp = kind(0.0d0)

  !> Type for 64-bit integers
  integer, parameter, public :: i8 = selected_int_kind(18)

  !> Indicates a standard Laplacian
  integer, parameter, public :: mg_laplacian = 1

  !> Indicates a variable-coefficient Laplacian
  integer, parameter, public :: mg_vlaplacian = 2

  !> Indicates a constant-coefficient Helmholtz equation
  integer, parameter, public :: mg_helmholtz = 3

  !> Indicates a variable-coefficient Helmholtz equation
  integer, parameter, public :: mg_vhelmholtz = 4

  !> Cartesian coordinate system
  integer, parameter, public :: mg_cartesian   = 1
  !> Cylindrical coordinate system
  integer, parameter, public :: mg_cylindrical = 2
  !> Spherical coordinate system
  integer, parameter, public :: mg_spherical   = 3

  integer, parameter, public :: mg_smoother_gs     = 1
  integer, parameter, public :: mg_smoother_gsrb   = 2
  integer, parameter, public :: mg_smoother_jacobi = 3

  !> Problem dimension
  integer, parameter, public :: mg_ndim = 1

  !> Number of predefined multigrid variables
  integer, parameter, public :: mg_num_vars = 4
  !> Maximum number of variables
  integer, parameter, public :: mg_max_num_vars = 10
  !> Index of solution
  integer, parameter, public :: mg_iphi = 1
  !> Index of right-hand side
  integer, parameter, public :: mg_irhs = 2
  !> Index of previous solution (used for correction)
  integer, parameter, public :: mg_iold = 3
  !> Index of residual
  integer, parameter, public :: mg_ires = 4

  !> Index of the variable coefficient (at cell centers)
  integer, parameter, public :: mg_iveps = 5

  !> Minimum allowed grid level
  integer, parameter, public :: mg_lvl_lo = -20
  !> Maximum allowed grid level
  integer, parameter, public :: mg_lvl_hi = 20

  !> Value to indicate a Dirichlet boundary condition
  integer, parameter, public :: mg_bc_dirichlet = -10

  !> Value to indicate a Neumann boundary condition
  integer, parameter, public :: mg_bc_neumann = -11

  !> Value to indicate a continuous boundary condition
  integer, parameter, public :: mg_bc_continuous = -12

  !> Special value that indicates there is no box
  integer, parameter, public :: mg_no_box = 0
  !> Special value that indicates there is a physical boundary
  integer, parameter, public :: mg_physical_boundary = -1

  !> Maximum number of timers to use
  integer, parameter, public :: mg_max_timers = 20

  ! Numbering of children (same location as **corners**)
  integer, parameter, public :: mg_num_children = 2

  ! Index offset for each child
  integer, parameter, public :: mg_child_dix(1, 2) = reshape([0,1], [1,2])
  ! Reverse child index in each direction
  integer, parameter, public :: mg_child_rev(2, 1) = reshape([2,1], [2,1])
  ! Children adjacent to a neighbor
  integer, parameter, public :: mg_child_adj_nb(1, 2) = reshape([1,2], [1,2])
  ! Which children have a low index per dimension
  logical, parameter, public :: mg_child_low(1, 2) = reshape([.true., .false.], [1, 2])

  ! Neighbor topology information
  integer, parameter, public :: mg_num_neighbors = 2
  integer, parameter, public :: mg_neighb_lowx = 1
  integer, parameter, public :: mg_neighb_highx = 2

  ! Index offsets of neighbors
  integer, parameter, public :: mg_neighb_dix(1, 2) = reshape([-1,1], [1,2])
  ! Which neighbors have a lower index
  logical, parameter, public :: mg_neighb_low(2) = [.true., .false.]
  ! Opposite of nb_low, but now as -1,1 integers
  integer, parameter, public :: mg_neighb_high_pm(2) = [-1, 1]

  ! Reverse neighbors
  integer, parameter, public :: mg_neighb_rev(2) = [2, 1]
  ! Direction (dimension) for a neighbor
  integer, parameter, public :: mg_neighb_dim(2) = [1, 1]

  !> Lists of blocks per refinement level
  type, public :: mg_lvl_t
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
  type, public :: mg_box_t
     integer  :: rank              !< Which process owns this box
     integer  :: id                !< Box id (index in boxes(:) array)
     integer  :: lvl               !< Refinement level
     integer  :: ix(1)          !< Spatial index
     integer  :: parent            !< Id of parent
     integer  :: children(2**1) !< Ids of children
     integer  :: neighbors(2*1) !< Ids of neighbors
     real(dp) :: r_min(1)       !< Minimum coordinate
     real(dp) :: dr(1)          !< Grid spacing
     !> Cell-centered data
     real(dp), allocatable :: cc(:, :)
  end type mg_box_t

  !> Buffer type (one is used for each pair of communicating processes)
  type, public :: mg_buf_t
     integer               :: i_send !< Index in send array
     integer               :: i_recv
     integer               :: i_ix
     integer, allocatable  :: ix(:) !< Will be used to sort the data
     real(dp), allocatable :: send(:)
     real(dp), allocatable :: recv(:)
  end type mg_buf_t

  type, public :: mg_comm_t
     integer, allocatable :: n_send(:, :)
     integer, allocatable :: n_recv(:, :)
  end type mg_comm_t

  type, public :: mg_bc_t
     integer  :: bc_type  = mg_bc_dirichlet !< Type of boundary condition
     real(dp) :: bc_value = 0.0_dp       !< Value (for e.g. Dirichlet or Neumann)
     !> To set user-defined boundary conditions (overrides bc(:))
     procedure(mg_subr_bc), pointer, nopass :: boundary_cond  => null()
     !> To set a user-defined refinement boundary method
     procedure(mg_subr_rb), pointer, nopass :: refinement_bnd => null()
  end type mg_bc_t

  type, public :: mg_timer_t
     character(len=20) :: name
     real(dp)          :: t = 0.0_dp
     real(dp)          :: t0
  end type mg_timer_t

  type, public :: mg_t
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
     integer                  :: domain_size_lvl(1, mg_lvl_lo:mg_lvl_hi)
     !> Grid spacing per level
     real(dp)                 :: dr(1, mg_lvl_lo:mg_lvl_hi)
     !> Minimum coordinates
     real(dp)                 :: r_min(1)
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

     !> Whether a dimension is periodic
     logical :: periodic(1) = .false.

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
     integer  :: coarsest_grid(1) = 2
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
       real(dp), intent(out)   :: bc(1)
     end subroutine mg_subr_bc

     !> To fill ghost cells near refinement boundaries
     subroutine mg_subr_rb(box, nc, iv, nb, cgc)
       import
       type(mg_box_t), intent(inout) :: box
       integer, intent(in)        :: nc
       integer, intent(in)        :: iv !< Index of variable
       integer, intent(in)        :: nb !< Direction
       !> Coarse data
       real(dp), intent(in)       :: cgc(1)
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
       integer, intent(in)       :: dix(1)        !< Offset of child in parent grid
       integer, intent(in)       :: nc               !< Child grid size
       integer, intent(in)       :: iv               !< Prolong from this variable
       real(dp), intent(out)     :: fine(nc) !< Prolonged values
     end subroutine mg_box_prolong
  end interface

  ! Public methods
  public :: mg_subr_bc
  public :: mg_subr_rb
  public :: mg_box_op
  public :: mg_box_gsrb
  public :: mg_box_prolong
  public :: mg_has_children
  public :: mg_ix_to_ichild
  public :: mg_get_child_offset
  public :: mg_highest_uniform_lvl
  public :: mg_number_of_unknowns
  public :: mg_get_face_coords
  public :: mg_add_timer
  public :: mg_timer_start
  public :: mg_timer_end
  public :: mg_timers_show

  !! File ../src/m_allocate_storage.f90

  public :: mg_allocate_storage
  public :: mg_deallocate_storage

  !! File ../src/m_diffusion.f90

!> Module for implicitly solving diffusion equations

  public :: diffusion_solve
  public :: diffusion_solve_vcoeff

  !! File ../src/m_laplacian.f90

!> Module which contains multigrid procedures for a Laplacian operator

  public :: laplacian_set_methods

  !! File ../src/m_vhelmholtz.f90

!> Module which contains multigrid procedures for a Helmholtz operator of the
!> form: div(D grad(phi)) - lambda*phi = f, where D has a smooth spatial
!> variation

  !> The lambda used for the Helmholtz equation (should be >= 0)
  real(dp), public, protected :: vhelmholtz_lambda = 0.0_dp

  public :: vhelmholtz_set_methods
  public :: vhelmholtz_set_lambda

  !! File ../src/m_build_tree.f90

  ! Public methods
  public :: mg_build_rectangle
  public :: mg_add_children
  public :: mg_set_leaves_parents
  public :: mg_set_next_level_ids
  public :: mg_set_refinement_boundaries
  public :: mg_set_neighbors_lvl

  !! File ../src/m_load_balance.f90
!> Module for load balancing a tree (that already has been constructed). The
!> load balancing determines which ranks (MPI processes) allocated physical
!> storage for boxes. The tree structure itself is present on all processes.

  ! Public methods
  public :: mg_load_balance_simple
  public :: mg_load_balance
  public :: mg_load_balance_parents

  !! File ../src/m_vlaplacian.f90

!> Module which contains multigrid procedures for a variable-coefficient
!> Laplacian operator, assuming the variation is smooth

  public :: vlaplacian_set_methods

  !! File ../src/m_communication.f90

  ! Public methods
  public :: mg_comm_init
  public :: sort_and_transfer_buffers

  !! File ../src/m_ghost_cells.f90

  ! Public methods
  public :: mg_ghost_cell_buffer_size
  public :: mg_fill_ghost_cells
  public :: mg_fill_ghost_cells_lvl
  public :: mg_phi_bc_store

  !! File ../src/m_prolong.f90

  ! Public methods
  public :: mg_prolong
  public :: mg_prolong_buffer_size
  public :: mg_prolong_sparse

  !! File ../src/m_helmholtz.f90

!> Module which contains multigrid procedures for a Helmholtz operator of the
!> form: laplacian(phi) - lambda*phi = f

  !> The lambda used for the Helmholtz equation (should be >= 0)
  real(dp), public, protected :: helmholtz_lambda = 0.0_dp

  public :: helmholtz_set_methods
  public :: helmholtz_set_lambda

  !! File ../src/m_multigrid.f90

  integer :: timer_total_vcycle  = -1
  integer :: timer_total_fmg     = -1
  integer :: timer_smoother      = -1
  integer :: timer_smoother_gc   = -1
  integer :: timer_coarse        = -1
  integer :: timer_correct       = -1
  integer :: timer_update_coarse = -1

  ! Public methods
  public :: mg_fas_vcycle
  public :: mg_fas_fmg
  public :: mg_set_methods
  public :: mg_apply_op

  !! File ../src/m_restrict.f90

  ! Public methods
  public :: mg_restrict
  public :: mg_restrict_lvl
  public :: mg_restrict_buffer_size

  !! File ../src/m_mrgrnk.f90

   Integer, Parameter :: kdp = selected_real_kind(15)
   public :: mrgrnk
   private :: kdp
   private :: I_mrgrnk
   interface mrgrnk
      module procedure I_mrgrnk
   end interface mrgrnk
contains

  !! File ../src/m_data_structures.f90

  !> Return .true. if a box has children
  elemental logical function mg_has_children(box)
    type(mg_box_t), intent(in) :: box

    ! Boxes are either fully refined or not, so we only need to check one of the
    ! children
    mg_has_children = (box%children(1) /= mg_no_box)
  end function mg_has_children

  !> Compute the child index for a box with spatial index ix. With child index
  !> we mean the index in the children(:) array of its parent.
  integer function mg_ix_to_ichild(ix)
    integer, intent(in) :: ix(1) !< Spatial index of the box
    ! The index can range from 1 (all ix odd) and 2**$D (all ix even)
    mg_ix_to_ichild = 2 - iand(ix(1), 1)
  end function mg_ix_to_ichild

  !> Get the offset of a box with respect to its parent (e.g. in 2d, there can
  !> be a child at offset 0,0, one at n_cell/2,0, one at 0,n_cell/2 and one at
  !> n_cell/2, n_cell/2)
  pure function mg_get_child_offset(mg, id) result(ix_offset)
    type(mg_t), intent(in) :: mg
    integer, intent(in)    :: id
    integer                :: ix_offset(1)

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
    real(dp), intent(out)      :: x(2)
    integer                    :: nb_dim
    real(dp)                   :: rmin(1)

    nb_dim = mg_neighb_dim(nb)


    rmin = box%r_min
    if (.not. mg_neighb_low(nb)) then
       rmin(nb_dim) = rmin(nb_dim) + box%dr(nb_dim) * nc
    end if

    x(1) = rmin(1)
    x(2) = rmin(1) + box%dr(1) * nc
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

  !! File ../src/m_allocate_storage.f90

  !> Deallocate all allocatable arrays
  subroutine mg_deallocate_storage(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: lvl

    if (.not. mg%is_allocated) &
         error stop "deallocate_storage: tree is not allocated"

    deallocate(mg%boxes)
    deallocate(mg%buf)

    deallocate(mg%comm_restrict%n_send)
    deallocate(mg%comm_restrict%n_recv)

    deallocate(mg%comm_prolong%n_send)
    deallocate(mg%comm_prolong%n_recv)

    deallocate(mg%comm_ghostcell%n_send)
    deallocate(mg%comm_ghostcell%n_recv)

    do lvl = mg%lowest_lvl, mg%highest_lvl
       deallocate(mg%lvls(lvl)%ids)
       deallocate(mg%lvls(lvl)%leaves)
       deallocate(mg%lvls(lvl)%parents)
       deallocate(mg%lvls(lvl)%ref_bnds)
       deallocate(mg%lvls(lvl)%my_ids)
       deallocate(mg%lvls(lvl)%my_leaves)
       deallocate(mg%lvls(lvl)%my_parents)
       deallocate(mg%lvls(lvl)%my_ref_bnds)
    end do

    mg%is_allocated       = .false.
    mg%n_boxes            = 0
    mg%phi_bc_data_stored = .false.
  end subroutine mg_deallocate_storage

  !> Allocate communication buffers and local boxes for a tree that has already
  !> been created
  subroutine mg_allocate_storage(mg)

    type(mg_t), intent(inout) :: mg
    integer                   :: i, id, lvl, nc
    integer                   :: n_send(0:mg%n_cpu-1, 3)
    integer                   :: n_recv(0:mg%n_cpu-1, 3)
    integer                   :: dsize(3)
    integer                   :: n_in, n_out, n_id

    if (.not. mg%tree_created) &
         error stop "allocate_storage: tree is not yet created"

    if (mg%is_allocated) &
         error stop "allocate_storage: tree is already allocated"

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          allocate(mg%boxes(id)%cc(0:nc+1, &
               mg_num_vars + mg%n_extra_vars))

          ! Set all initial values to zero
          mg%boxes(id)%cc(:, :) = 0.0_dp
       end do
    end do

    allocate(mg%buf(0:mg%n_cpu-1))

    call mg_ghost_cell_buffer_size(mg, n_send(:, 1), &
         n_recv(:, 1), dsize(1))
    call mg_restrict_buffer_size(mg, n_send(:, 2), &
         n_recv(:, 2), dsize(2))
    call mg_prolong_buffer_size(mg, n_send(:, 3), &
         n_recv(:, 3), dsize(3))

    do i = 0, mg%n_cpu-1
       n_out = maxval(n_send(i, :) * dsize(:))
       n_in = maxval(n_recv(i, :) * dsize(:))
       n_id = maxval(n_send(i, :))
       allocate(mg%buf(i)%send(n_out))
       allocate(mg%buf(i)%recv(n_in))
       allocate(mg%buf(i)%ix(n_id))
    end do

    mg%is_allocated = .true.
  end subroutine mg_allocate_storage

  !! File ../src/m_diffusion.f90

  !> Solve a diffusion equation implicitly, assuming a constant diffusion
  !> coefficient. The solution at time t should be stored in mg_iphi, which is
  !> on output replaced by the solution at time t+dt.
  subroutine diffusion_solve(mg, dt, diffusion_coeff, order, max_res)

    type(mg_t), intent(inout) :: mg
    real(dp), intent(in)      :: dt
    real(dp), intent(in)      :: diffusion_coeff
    integer, intent(in)       :: order
    real(dp), intent(in)      :: max_res
    integer, parameter        :: max_its = 10
    integer                   :: n
    real(dp)                  :: res

    mg%operator_type = mg_helmholtz
    call mg_set_methods(mg)

    select case (order)
    case (1)
       call helmholtz_set_lambda(1/(dt * diffusion_coeff))
       call set_rhs(mg, -1/(dt * diffusion_coeff), 0.0_dp)
    case (2)
       call helmholtz_set_lambda(0.0d0)
       call mg_apply_op(mg, mg_irhs)
       call helmholtz_set_lambda(2/(dt * diffusion_coeff))
       call set_rhs(mg, -2/(dt * diffusion_coeff), -1.0_dp)
    case default
       error stop "diffusion_solve: order should be 1 or 2"
    end select

    ! Start with an FMG cycle
    call mg_fas_fmg(mg, .true., max_res=res)

    ! Add V-cycles if necessary
    do n = 1, max_its
       if (res <= max_res) exit
       call mg_fas_vcycle(mg, max_res=res)
    end do

    if (n == max_its + 1) then
       if (mg%my_rank == 0) &
            print *, "Did you specify boundary conditions correctly?"
       error stop "diffusion_solve: no convergence"
    end if
  end subroutine diffusion_solve

  !> Solve a diffusion equation implicitly, assuming a variable diffusion
  !> coefficient which has been stored in mg_iveps (also on coarse grids). The
  !> solution at time t should be stored in mg_iphi, which is on output replaced
  !> by the solution at time t+dt.
  subroutine diffusion_solve_vcoeff(mg, dt, order, max_res)

    type(mg_t), intent(inout) :: mg
    real(dp), intent(in)      :: dt
    integer, intent(in)       :: order
    real(dp), intent(in)      :: max_res
    integer, parameter        :: max_its = 10
    integer                   :: n
    real(dp)                  :: res

    mg%operator_type = mg_vhelmholtz
    call mg_set_methods(mg)

    select case (order)
    case (1)
       call vhelmholtz_set_lambda(1/dt)
       call set_rhs(mg, -1/dt, 0.0_dp)
    case (2)
       call vhelmholtz_set_lambda(0.0d0)
       call mg_apply_op(mg, mg_irhs)
       call vhelmholtz_set_lambda(2/dt)
       call set_rhs(mg, -2/dt, -1.0_dp)
    case default
       error stop "diffusion_solve: order should be 1 or 2"
    end select

    ! Start with an FMG cycle
    call mg_fas_fmg(mg, .true., max_res=res)

    ! Add V-cycles if necessary
    do n = 1, max_its
       if (res <= max_res) exit
       call mg_fas_vcycle(mg, max_res=res)
    end do

    if (n == max_its + 1) then
       if (mg%my_rank == 0) then
          print *, "Did you specify boundary conditions correctly?"
          print *, "Or is the variation in diffusion too large?"
       end if
       error stop "diffusion_solve: no convergence"
    end if
  end subroutine diffusion_solve_vcoeff

  subroutine set_rhs(mg, f1, f2)
    type(mg_t), intent(inout) :: mg
    real(dp), intent(in)      :: f1, f2
    integer                   :: lvl, i, id, nc

    do lvl = 1, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do i = 1, size(mg%lvls(lvl)%my_leaves)
          id = mg%lvls(lvl)%my_leaves(i)
          mg%boxes(id)%cc(1:nc, mg_irhs) = &
               f1 * mg%boxes(id)%cc(1:nc, mg_iphi) + &
               f2 * mg%boxes(id)%cc(1:nc, mg_irhs)
       end do
    end do
  end subroutine set_rhs

  !! File ../src/m_laplacian.f90

  subroutine laplacian_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    if (all(mg%periodic)) then
       ! For a fully periodic Laplacian, remove the mean from the rhs and phi so
       ! that a unique and periodic solution can be found
       mg%subtract_mean = .true.
    end if

    select case (mg%geometry_type)
    case (mg_cartesian)
       mg%box_op => box_lpl

       select case (mg%smoother_type)
       ! case (smoother_jacobi)
       !    mg%box_smoother => box_jacobi_lpl
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_lpl
       case default
          error stop "laplacian_set_methods: unsupported smoother type"
       end select
    case default
       error stop "laplacian_set_methods: unsupported geometry"
    end select

  end subroutine laplacian_set_methods

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine box_gs_lpl(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: i, i0, di
    real(dp)                  :: idr2(1), fac
    logical                   :: redblack

    idr2 = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    fac = 0.5_dp / sum(idr2)
    i0  = 1

    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
      if (redblack) i0 = 2 - iand(redblack_cntr, 1)

      do i = i0, nc, di
         cc(i, n) = fac * ( &
              idr2(1) * (cc(i+1, n) + cc(i-1, n)) - &
              cc(i, mg_irhs))
      end do
    end associate
  end subroutine box_gs_lpl

!   !> Perform Jacobi relaxation on box for a Laplacian operator
!   subroutine box_jacobi_lpl(mg, id, nc, cntr)
!     type(mg_t), intent(inout) :: mg
!     integer, intent(in)       :: id
!     integer, intent(in)       :: nc
!     integer, intent(in)       :: cntr !< Not used
!     integer                   :: i
!     real(dp), parameter       :: w     = 2.0_dp / 3
!     real(dp)                  :: tmp(0:nc+1)
!     real(dp)                  :: dr2
! #if 1 == 3
!     real(dp), parameter       :: sixth = 1/6.0_dp
! #endif

!     dr2   = mg%dr(mg%boxes(id)%lvl)**2

!     associate (box => mg%boxes(id))
!       tmp = box%cc(:, mg_iphi)
!       do i=1, nc
! #if 1 == 2
!          box%cc(i, j, mg_iphi) = (1-w) * box%cc(i, j, mg_iphi) + &
!               0.25_dp * w * ( &
!               tmp(i+1, j) + tmp(i-1, j) + &
!               tmp(i, j+1) + tmp(i, j-1) - &
!               dr2 * box%cc(i, j, mg_irhs))
! #elif 1 == 3
!          box%cc(i, j, k, mg_iphi) = (1-w) * &
!               box%cc(i, j, k, mg_iphi) + &
!               sixth * w * ( &
!               tmp(i+1, j, k) + tmp(i-1, j, k) + &
!               tmp(i, j+1, k) + tmp(i, j-1, k) + &
!               tmp(i, j, k+1) + tmp(i, j, k-1) - &
!               dr2 * box%cc(i, j, k, mg_irhs))
! #endif
!       end do; 
!     end associate
!   end subroutine box_jacobi_lpl

  !> Perform Laplacian operator on a box
  subroutine box_lpl(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: i
    real(dp)                  :: idr2(1)

    idr2 = 1 / mg%dr(:, mg%boxes(id)%lvl)**2

    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
      do i = 1, nc
            cc(i, i_out) = &
                 idr2(1) * (cc(i-1, n) + cc(i+1, n) - 2 * cc(i, n))
         end do
    end associate
  end subroutine box_lpl


  !! File ../src/m_vhelmholtz.f90

  subroutine vhelmholtz_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    if (mg%n_extra_vars == 0 .and. mg%is_allocated) then
       error stop "vhelmholtz_set_methods: mg%n_extra_vars == 0"
    else
       mg%n_extra_vars = max(1, mg%n_extra_vars)
    end if

    mg%subtract_mean = .false.

    ! Use Neumann zero boundary conditions for the variable coefficient, since
    ! it is needed in ghost cells.
    mg%bc(:, mg_iveps)%bc_type = mg_bc_neumann
    mg%bc(:, mg_iveps)%bc_value = 0.0_dp

    select case (mg%geometry_type)
    case (mg_cartesian)
       mg%box_op => box_vhelmh

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_vhelmh
       case default
          error stop "vhelmholtz_set_methods: unsupported smoother type"
       end select
    case default
       error stop "vhelmholtz_set_methods: unsupported geometry"
    end select

  end subroutine vhelmholtz_set_methods

  subroutine vhelmholtz_set_lambda(lambda)
    real(dp), intent(in) :: lambda

    if (lambda < 0) &
         error stop "vhelmholtz_set_lambda: lambda < 0 not allowed"

    vhelmholtz_lambda = lambda
  end subroutine vhelmholtz_set_lambda

  !> Perform Gauss-Seidel relaxation on box for a Helmholtz operator
  subroutine box_gs_vhelmh(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: i, i0, di
    logical                   :: redblack
    real(dp)                  :: idr2(2*1), u(2*1)
    real(dp)                  :: a0, a(2*1), c(2*1)

    ! Duplicate 1/dr^2 array to multiply neighbor values
    idr2(1:2*1:2) = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    idr2(2:2*1:2) = idr2(1:2*1:2)
    i0  = 1

    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
         i_eps => mg_iveps)
      if (redblack) i0 = 2 - iand(redblack_cntr, 1)

      do i = i0, nc, di
         a0     = cc(i, i_eps)
         u(1:2) = cc(i-1:i+1:2, n)
         a(1:2) = cc(i-1:i+1:2, i_eps)
         c(:)   = 2 * a0 * a(:) / (a0 + a(:)) * idr2

         cc(i, n) = &
              (sum(c(:) * u(:)) - cc(i, mg_irhs)) / &
              (sum(c(:)) + vhelmholtz_lambda)
      end do
    end associate
  end subroutine box_gs_vhelmh

  !> Perform Helmholtz operator on a box
  subroutine box_vhelmh(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Helmholtz in
    integer                   :: i
    real(dp)                  :: idr2(2*1), a0, u0, u(2*1), a(2*1)

    ! Duplicate 1/dr^2 array to multiply neighbor values
    idr2(1:2*1:2) = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    idr2(2:2*1:2) = idr2(1:2*1:2)

    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
         i_eps => mg_iveps)
      do i = 1, nc
         a0     = cc(i, i_eps)
         a(1:2) = cc(i-1:i+1:2, i_eps)
         u0     = cc(i, n)
         u(1:2) = cc(i-1:i+1:2, n)

         cc(i, i_out) = sum(2 * idr2 * &
              a0*a(:)/(a0 + a(:)) * (u(:) - u0)) - &
              vhelmholtz_lambda * u0
      end do
    end associate
  end subroutine box_vhelmh

  !! File ../src/m_build_tree.f90

  subroutine mg_build_rectangle(mg, domain_size, box_size, dx, r_min, &
       periodic, n_finer)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: domain_size(1)
    integer, intent(in)       :: box_size
    real(dp), intent(in)      :: dx(1)
    real(dp), intent(in)      :: r_min(1)
    logical, intent(in)       :: periodic(1)
    integer, intent(in)       :: n_finer
    integer                   :: i, lvl, n, id, nx(1)
    integer                   :: boxes_per_dim(1, mg_lvl_lo:1)
    integer                   :: periodic_offset(1)

    if (modulo(box_size, 2) /= 0) &
         error stop "box_size should be even"
    if (any(modulo(domain_size, box_size) /= 0)) &
         error stop "box_size does not divide domain_size"

    nx                       = domain_size
    mg%box_size              = box_size
    mg%box_size_lvl(1)       = box_size
    mg%domain_size_lvl(:, 1) = domain_size
    mg%first_normal_lvl      = 1
    mg%dr(:, 1)              = dx
    mg%r_min(:)              = r_min
    mg%periodic              = periodic
    boxes_per_dim(:, :)      = 0
    boxes_per_dim(:, 1)      = domain_size / box_size

    do lvl = 1, mg_lvl_lo+1, -1
       ! For a Gauss-Seidel (non red-black) smoother, we should avoid boxes
       ! containing a single cell
       if (any(modulo(nx, 2) == 1 .or. nx == mg%coarsest_grid .or. &
            (mg%box_size_lvl(lvl) == mg%coarsest_grid .and. &
            mg%smoother_type == mg_smoother_gs))) exit

       if (all(modulo(nx/mg%box_size_lvl(lvl), 2) == 0)) then
          mg%box_size_lvl(lvl-1) = mg%box_size_lvl(lvl)
          boxes_per_dim(:, lvl-1) = boxes_per_dim(:, lvl)/2
          mg%first_normal_lvl = lvl-1
       else
          mg%box_size_lvl(lvl-1) = mg%box_size_lvl(lvl)/2
          boxes_per_dim(:, lvl-1) = boxes_per_dim(:, lvl)
       end if

       mg%dr(:, lvl-1)              = mg%dr(:, lvl) * 2
       nx                           = nx / 2
       mg%domain_size_lvl(:, lvl-1) = nx
    end do

    mg%lowest_lvl = lvl
    mg%highest_lvl = 1

    do lvl = 2, mg_lvl_hi
       mg%dr(:, lvl) = mg%dr(:, lvl-1) * 0.5_dp
       mg%box_size_lvl(lvl) = box_size
       mg%domain_size_lvl(:, lvl) = 2 * mg%domain_size_lvl(:, lvl-1)
    end do

    n = sum(product(boxes_per_dim, dim=1)) + n_finer
    allocate(mg%boxes(n))

    ! Create lowest level
    nx = boxes_per_dim(:, mg%lowest_lvl)
    periodic_offset = [nx(1)-1]

    do i=1,nx(1)
       mg%n_boxes = mg%n_boxes + 1
       n          = mg%n_boxes

       mg%boxes(n)%rank        = 0
       mg%boxes(n)%id          = n
       mg%boxes(n)%lvl         = mg%lowest_lvl
       mg%boxes(n)%ix(:)       = [i]
       mg%boxes(n)%r_min(:)    = r_min + (mg%boxes(n)%ix(:) - 1) * &
            mg%box_size_lvl(mg%lowest_lvl) * mg%dr(:, mg%lowest_lvl)
       mg%boxes(n)%dr(:)       = mg%dr(:, mg%lowest_lvl)
       mg%boxes(n)%parent      = mg_no_box
       mg%boxes(n)%children(:) = mg_no_box

       ! Set default neighbors
       mg%boxes(n)%neighbors(:) = [n-1, n+1]

       ! Handle boundaries
       where ([i] == 1 .and. .not. periodic)
          mg%boxes(n)%neighbors(1:mg_num_neighbors:2) = &
               mg_physical_boundary
       end where
       where ([i] == 1 .and. periodic)
          mg%boxes(n)%neighbors(1:mg_num_neighbors:2) = &
               n + periodic_offset
       end where

       where ([i] == nx .and. .not. periodic)
          mg%boxes(n)%neighbors(2:mg_num_neighbors:2) = &
               mg_physical_boundary
       end where
       where ([i] == nx .and. periodic)
          mg%boxes(n)%neighbors(2:mg_num_neighbors:2) = &
               n - periodic_offset
       end where
    end do; 

    mg%lvls(mg%lowest_lvl)%ids = [(n, n=1, mg%n_boxes)]

    ! Add higher levels
    do lvl = mg%lowest_lvl, 0
       if (mg%box_size_lvl(lvl+1) == mg%box_size_lvl(lvl)) then
          do i = 1, size(mg%lvls(lvl)%ids)
             id = mg%lvls(lvl)%ids(i)
             call mg_add_children(mg, id)
          end do

          call mg_set_leaves_parents(mg%boxes, mg%lvls(lvl))
          call mg_set_next_level_ids(mg, lvl)
          call mg_set_neighbors_lvl(mg, lvl+1)
       else
          do i = 1, size(mg%lvls(lvl)%ids)
             id = mg%lvls(lvl)%ids(i)
             call add_single_child(mg, id, size(mg%lvls(lvl)%ids))
          end do

          call mg_set_leaves_parents(mg%boxes, mg%lvls(lvl))
          call mg_set_next_level_ids(mg, lvl)
       end if
    end do

    call mg_set_leaves_parents(mg%boxes, mg%lvls(1))

    ! No refinement boundaries
    do lvl = mg%lowest_lvl, 1
       if (allocated(mg%lvls(lvl)%ref_bnds)) &
            deallocate(mg%lvls(lvl)%ref_bnds)
       allocate(mg%lvls(lvl)%ref_bnds(0))
    end do

    mg%tree_created = .true.
  end subroutine mg_build_rectangle

  subroutine mg_set_neighbors_lvl(mg, lvl)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer                   :: i, id

    do i = 1, size(mg%lvls(lvl)%ids)
       id = mg%lvls(lvl)%ids(i)
       call set_neighbs(mg%boxes, id)
    end do
  end subroutine mg_set_neighbors_lvl

  subroutine mg_set_next_level_ids(mg, lvl)
    type(mg_t), intent(inout)  :: mg
    integer, intent(in)        :: lvl
    integer                    :: n, i, id

    if (allocated(mg%lvls(lvl+1)%ids)) &
         deallocate(mg%lvls(lvl+1)%ids)

    ! Set next level ids to children of this level
    if (mg%box_size_lvl(lvl+1) == mg%box_size_lvl(lvl)) then
       n = mg_num_children * size(mg%lvls(lvl)%parents)
       allocate(mg%lvls(lvl+1)%ids(n))

       n = mg_num_children
       do i = 1, size(mg%lvls(lvl)%parents)
          id = mg%lvls(lvl)%parents(i)
          mg%lvls(lvl+1)%ids(n*(i-1)+1:n*i) = mg%boxes(id)%children
       end do
    else
       n = size(mg%lvls(lvl)%parents)
       allocate(mg%lvls(lvl+1)%ids(n))

       n = 1
       do i = 1, size(mg%lvls(lvl)%parents)
          id = mg%lvls(lvl)%parents(i)
          mg%lvls(lvl+1)%ids(i) = mg%boxes(id)%children(1)
       end do
    end if

  end subroutine mg_set_next_level_ids

  ! Set the neighbors of id (using their parent)
  subroutine set_neighbs(boxes, id)
    type(mg_box_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    integer                     :: nb, nb_id

    do nb = 1, mg_num_neighbors
       if (boxes(id)%neighbors(nb) == mg_no_box) then
          nb_id = find_neighb(boxes, id, nb)
          if (nb_id > mg_no_box) then
             boxes(id)%neighbors(nb) = nb_id
             boxes(nb_id)%neighbors(mg_neighb_rev(nb)) = id
          end if
       end if
    end do
  end subroutine set_neighbs

  !> Get the id of neighbor nb of boxes(id), through its parent
  function find_neighb(boxes, id, nb) result(nb_id)
    type(mg_box_t), intent(in) :: boxes(:) !< List with all the boxes
    integer, intent(in)      :: id       !< Box whose neighbor we are looking for
    integer, intent(in)      :: nb       !< Neighbor index
    integer                  :: nb_id, p_id, c_ix, d, old_pid

    p_id    = boxes(id)%parent
    old_pid = p_id
    c_ix    = mg_ix_to_ichild(boxes(id)%ix)
    d       = mg_neighb_dim(nb)

    ! Check if neighbor is in same direction as ix is (low/high). If so,
    ! use neighbor of parent
    if (mg_child_low(d, c_ix) .eqv. mg_neighb_low(nb)) then
       p_id = boxes(p_id)%neighbors(nb)
    end if

    ! The child ix of the neighbor is reversed in direction d
    nb_id = boxes(p_id)%children(mg_child_rev(c_ix, d))
  end function find_neighb

  !> Create a list of leaves and a list of parents for a level
  subroutine mg_set_leaves_parents(boxes, level)
    type(mg_box_t), intent(in)   :: boxes(:) !< List of boxes
    type(mg_lvl_t), intent(inout) :: level !< Level type which contains the indices of boxes
    integer                    :: i, id, i_leaf, i_parent
    integer                    :: n_parents, n_leaves

    n_parents = count(mg_has_children(boxes(level%ids)))
    n_leaves = size(level%ids) - n_parents

    if (.not. allocated(level%parents)) then
       allocate(level%parents(n_parents))
    else if (n_parents /= size(level%parents)) then
       deallocate(level%parents)
       allocate(level%parents(n_parents))
    end if

    if (.not. allocated(level%leaves)) then
       allocate(level%leaves(n_leaves))
    else if (n_leaves /= size(level%leaves)) then
       deallocate(level%leaves)
       allocate(level%leaves(n_leaves))
    end if

    i_leaf   = 0
    i_parent = 0
    do i = 1, size(level%ids)
       id = level%ids(i)
       if (mg_has_children(boxes(id))) then
          i_parent                = i_parent + 1
          level%parents(i_parent) = id
       else
          i_leaf               = i_leaf + 1
          level%leaves(i_leaf) = id
       end if
    end do
  end subroutine mg_set_leaves_parents

  !> Create a list of refinement boundaries (from the coarse side)
  subroutine mg_set_refinement_boundaries(boxes, level)
    type(mg_box_t), intent(in)    :: boxes(:)
    type(mg_lvl_t), intent(inout) :: level
    integer, allocatable       :: tmp(:)
    integer                    :: i, id, nb, nb_id, ix

    if (allocated(level%ref_bnds)) deallocate(level%ref_bnds)

    if (size(level%parents) == 0) then
       ! There are no refinement boundaries
       allocate(level%ref_bnds(0))
    else
       allocate(tmp(size(level%leaves)))
       ix = 0
       do i = 1, size(level%leaves)
          id = level%leaves(i)

          do nb = 1, mg_num_neighbors
             nb_id = boxes(id)%neighbors(nb)
             if (nb_id > mg_no_box) then
                if (mg_has_children(boxes(nb_id))) then
                   ix = ix + 1
                   tmp(ix) = id
                   exit
                end if
             end if
          end do
       end do

       allocate(level%ref_bnds(ix))
       level%ref_bnds(:) = tmp(1:ix)
    end if
  end subroutine mg_set_refinement_boundaries

  subroutine mg_add_children(mg, id)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id !< Id of box that gets children
    integer                   :: lvl, i, nb, child_nb(2**(1-1))
    integer                   :: c_ids(mg_num_children)
    integer                   :: c_id, c_ix_base(1)

    if (mg%n_boxes + mg_num_children > size(mg%boxes)) then
       error stop "mg_add_children: not enough space"
    end if

    c_ids                 = [(mg%n_boxes+i, i=1,mg_num_children)]
    mg%n_boxes            = mg%n_boxes + mg_num_children
    mg%boxes(id)%children = c_ids
    c_ix_base             = 2 * mg%boxes(id)%ix - 1
    lvl                   = mg%boxes(id)%lvl+1

    do i = 1, mg_num_children
       c_id                     = c_ids(i)
       mg%boxes(c_id)%rank      = mg%boxes(id)%rank
       mg%boxes(c_id)%ix        = c_ix_base + mg_child_dix(:, i)
       mg%boxes(c_id)%lvl       = lvl
       mg%boxes(c_id)%parent    = id
       mg%boxes(c_id)%children  = mg_no_box
       mg%boxes(c_id)%neighbors = mg_no_box
       mg%boxes(c_id)%r_min     = mg%boxes(id)%r_min + &
            mg%dr(:, lvl) * mg_child_dix(:, i) * mg%box_size
       mg%boxes(c_id)%dr(:)     = mg%dr(:, lvl)
    end do

    ! Set boundary conditions at children
    do nb = 1, mg_num_neighbors
       if (mg%boxes(id)%neighbors(nb) < mg_no_box) then
          child_nb = c_ids(mg_child_adj_nb(:, nb)) ! Neighboring children
          mg%boxes(child_nb)%neighbors(nb) = mg%boxes(id)%neighbors(nb)
       end if
    end do
  end subroutine mg_add_children

  subroutine add_single_child(mg, id, n_boxes_lvl)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id !< Id of box that gets children
    integer, intent(in)       :: n_boxes_lvl
    integer                   :: lvl, c_id

    c_id                     = mg%n_boxes + 1
    mg%n_boxes               = mg%n_boxes + 1
    mg%boxes(id)%children(1) = c_id
    lvl                      = mg%boxes(id)%lvl+1

    mg%boxes(c_id)%rank      = mg%boxes(id)%rank
    mg%boxes(c_id)%ix        = mg%boxes(id)%ix
    mg%boxes(c_id)%lvl       = lvl
    mg%boxes(c_id)%parent    = id
    mg%boxes(c_id)%children  = mg_no_box
    where (mg%boxes(id)%neighbors == mg_physical_boundary)
       mg%boxes(c_id)%neighbors = mg%boxes(id)%neighbors
    elsewhere
       mg%boxes(c_id)%neighbors = mg%boxes(id)%neighbors + n_boxes_lvl
    end where
    mg%boxes(c_id)%r_min = mg%boxes(id)%r_min
    mg%boxes(c_id)%dr(:) = mg%dr(:, lvl)

  end subroutine add_single_child

  !! File ../src/m_load_balance.f90

  !> Load balance all boxes in the multigrid tree, by simply distributing the
  !> load per grid level. This method will only work well for uniform grids.
  !>
  !> Note that in a typical application the load balancing of the leaves is
  !> already determined, then mg_load_balance_parents can be used.
  subroutine mg_load_balance_simple(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: i, id, lvl, single_cpu_lvl
    integer                   :: work_left, my_work, i_cpu

    ! Up to this level, all boxes have to be on a single processor because they
    ! have a different size and the communication routines do not support this
    single_cpu_lvl = max(mg%first_normal_lvl-1, mg%lowest_lvl)

    do lvl = mg%lowest_lvl, single_cpu_lvl
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          mg%boxes(id)%rank = 0
       end do
    end do

    ! Distribute the boxes equally. Due to the way the mesh is constructed, the
    ! mg%lvls(lvl)%ids array already contains a Morton-like ordering.
    do lvl = single_cpu_lvl+1, mg%highest_lvl
       work_left = size(mg%lvls(lvl)%ids)
       my_work   = 0
       i_cpu     = 0

       do i = 1, size(mg%lvls(lvl)%ids)
          if ((mg%n_cpu - i_cpu - 1) * my_work >= work_left) then
             i_cpu   = i_cpu + 1
             my_work = 0
          end if

          my_work = my_work + 1
          work_left = work_left - 1

          id = mg%lvls(lvl)%ids(i)
          mg%boxes(id)%rank = i_cpu
       end do
    end do

    do lvl = mg%lowest_lvl, mg%highest_lvl
       call update_lvl_info(mg, mg%lvls(lvl))
    end do

  end subroutine mg_load_balance_simple

  !> Load balance all boxes in the multigrid tree. Compared to
  !> mg_load_balance_simple, this method does a better job of setting the ranks
  !> of parent boxes
  !>
  !> Note that in a typical application the load balancing of the leaves is
  !> already determined, then mg_load_balance_parents can be used.
  subroutine mg_load_balance(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: i, id, lvl, single_cpu_lvl
    integer                   :: work_left, my_work(0:mg%n_cpu), i_cpu
    integer                   :: c_ids(mg_num_children)
    integer                   :: c_ranks(mg_num_children)
    integer                   :: coarse_rank

    ! Up to this level, all boxes have to be on a single processor because they
    ! have a different size and the communication routines do not support this
    single_cpu_lvl = max(mg%first_normal_lvl-1, mg%lowest_lvl)

    ! Distribute the boxes equally. Due to the way the mesh is constructed, the
    ! mg%lvls(lvl)%ids array already contains a Morton-like ordering.
    do lvl = mg%highest_lvl, single_cpu_lvl+1, -1
       ! For parents determine the rank based on their child ranks
       my_work(:) = 0

       do i = 1, size(mg%lvls(lvl)%parents)
          id = mg%lvls(lvl)%parents(i)

          c_ids = mg%boxes(id)%children
          c_ranks = mg%boxes(c_ids)%rank
          i_cpu = most_popular(c_ranks, my_work, mg%n_cpu)
          mg%boxes(id)%rank = i_cpu
          my_work(i_cpu) = my_work(i_cpu) + 1
       end do

       work_left = size(mg%lvls(lvl)%leaves)
       i_cpu     = 0

       do i = 1, size(mg%lvls(lvl)%leaves)
          ! Skip this CPU if it already has enough work
          if ((mg%n_cpu - i_cpu - 1) * my_work(i_cpu) >= &
               work_left + sum(my_work(i_cpu+1:))) then
             i_cpu = i_cpu + 1
          end if

          my_work(i_cpu) = my_work(i_cpu) + 1
          work_left = work_left - 1

          id = mg%lvls(lvl)%leaves(i)
          mg%boxes(id)%rank = i_cpu
       end do
    end do

    ! Determine most popular CPU for coarse grids
    if (single_cpu_lvl < mg%highest_lvl) then
       coarse_rank = most_popular(mg%boxes(&
            mg%lvls(single_cpu_lvl+1)%ids)%rank, my_work, mg%n_cpu)
    else
       coarse_rank = 0
    end if

    do lvl = mg%lowest_lvl, single_cpu_lvl
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          mg%boxes(id)%rank = coarse_rank
       end do
    end do

    do lvl = mg%lowest_lvl, mg%highest_lvl
       call update_lvl_info(mg, mg%lvls(lvl))
    end do

  end subroutine mg_load_balance

  !> Load balance the parents (non-leafs). Assign them to the rank that has most
  !> children.
  subroutine mg_load_balance_parents(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: i, id, lvl
    integer                   :: c_ids(mg_num_children)
    integer                   :: c_ranks(mg_num_children)
    integer                   :: single_cpu_lvl, coarse_rank
    integer                   :: my_work(0:mg%n_cpu), i_cpu

    ! Up to this level, all boxes have to be on a single processor because they
    ! have a different size and the communication routines do not support this
    single_cpu_lvl = max(mg%first_normal_lvl-1, mg%lowest_lvl)

    do lvl = mg%highest_lvl-1, single_cpu_lvl+1, -1
       my_work(:) = 0

       ! Determine amount of work for the leaves
       do i = 1, size(mg%lvls(lvl)%leaves)
          id = mg%lvls(lvl)%leaves(i)
          i_cpu = mg%boxes(id)%rank
          my_work(i_cpu) = my_work(i_cpu) + 1
       end do

       do i = 1, size(mg%lvls(lvl)%parents)
          id = mg%lvls(lvl)%parents(i)

          c_ids = mg%boxes(id)%children
          c_ranks = mg%boxes(c_ids)%rank
          i_cpu = most_popular(c_ranks, my_work, mg%n_cpu)
          mg%boxes(id)%rank = i_cpu
          my_work(i_cpu) = my_work(i_cpu) + 1
       end do

    end do

    ! Determine most popular CPU for coarse grids
    if (single_cpu_lvl < mg%highest_lvl) then
       coarse_rank = most_popular(mg%boxes(&
            mg%lvls(single_cpu_lvl+1)%ids)%rank, my_work, mg%n_cpu)
    else
       coarse_rank = 0
    end if

    do lvl = mg%lowest_lvl, single_cpu_lvl
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          mg%boxes(id)%rank = coarse_rank
       end do
    end do

    do lvl = mg%lowest_lvl, mg%highest_lvl
       call update_lvl_info(mg, mg%lvls(lvl))
    end do

  end subroutine mg_load_balance_parents

  !> Determine most popular rank in the list. In case of ties, assign the rank
  !> with the least work.
  pure integer function most_popular(list, work, n_cpu)
    integer, intent(in) :: list(:) !< List of MPI ranks
    integer, intent(in) :: n_cpu
    integer, intent(in) :: work(0:n_cpu-1) !< Existing work per rank
    integer             :: i, best_count, current_count
    integer             :: my_work, best_work

    best_count   = 0
    best_work    = 0
    most_popular = -1

    do i = 1, size(list)
       current_count = count(list == list(i))
       my_work       = work(list(i))

       ! In case of ties, select task with lowest work
       if (current_count > best_count .or. &
            current_count == best_count .and. my_work < best_work) then
          best_count   = current_count
          best_work    = my_work
          most_popular = list(i)
       end if
    end do

  end function most_popular

  subroutine update_lvl_info(mg, lvl)
    type(mg_t), intent(inout)     :: mg
    type(mg_lvl_t), intent(inout) :: lvl

    lvl%my_ids = pack(lvl%ids, &
         mg%boxes(lvl%ids)%rank == mg%my_rank)
    lvl%my_leaves = pack(lvl%leaves, &
         mg%boxes(lvl%leaves)%rank == mg%my_rank)
    lvl%my_parents = pack(lvl%parents, &
         mg%boxes(lvl%parents)%rank == mg%my_rank)
    lvl%my_ref_bnds = pack(lvl%ref_bnds, &
         mg%boxes(lvl%ref_bnds)%rank == mg%my_rank)
  end subroutine update_lvl_info

  !! File ../src/m_vlaplacian.f90

  subroutine vlaplacian_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    if (mg%n_extra_vars == 0 .and. mg%is_allocated) then
       error stop "vlaplacian_set_methods: mg%n_extra_vars == 0"
    else
       mg%n_extra_vars = max(1, mg%n_extra_vars)
    end if

    mg%subtract_mean = .false.

    ! Use Neumann zero boundary conditions for the variable coefficient, since
    ! it is needed in ghost cells.
    mg%bc(:, mg_iveps)%bc_type = mg_bc_neumann
    mg%bc(:, mg_iveps)%bc_value = 0.0_dp

    select case (mg%geometry_type)
    case (mg_cartesian)
       mg%box_op => box_vlpl

       ! TODO: test whether a custom prolongation works better
       ! mg%box_prolong => vlpl_prolong

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_vlpl
       case default
          error stop "vlaplacian_set_methods: unsupported smoother type"
       end select

    case default
       error stop "vlaplacian_set_methods: unsupported geometry"
    end select

  end subroutine vlaplacian_set_methods

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine box_gs_vlpl(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: i, i0, di
    logical                   :: redblack
    real(dp)                  :: idr2(2*1), u(2*1)
    real(dp)                  :: a0, a(2*1), c(2*1)

    ! Duplicate 1/dr^2 array to multiply neighbor values
    idr2(1:2*1:2) = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    idr2(2:2*1:2) = idr2(1:2*1:2)
    i0  = 1

    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
         i_eps => mg_iveps)
      if (redblack) i0 = 2 - iand(redblack_cntr, 1)

      do i = i0, nc, di
         a0     = cc(i, i_eps)
         u(1:2) = cc(i-1:i+1:2, n)
         a(1:2) = cc(i-1:i+1:2, i_eps)
         c(:)   = 2 * a0 * a(:) / (a0 + a(:)) * idr2

         cc(i, n) = &
              (sum(c(:) * u(:)) - cc(i, mg_irhs)) / sum(c(:))
      end do
    end associate
  end subroutine box_gs_vlpl

  !> Perform Laplacian operator on a box
  subroutine box_vlpl(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: i
    real(dp)                  :: idr2(2*1), a0, u0, u(2*1), a(2*1)

    ! Duplicate 1/dr^2 array to multiply neighbor values
    idr2(1:2*1:2) = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    idr2(2:2*1:2) = idr2(1:2*1:2)

    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
         i_eps => mg_iveps)
      do i = 1, nc
         a0     = cc(i, i_eps)
         a(1:2) = cc(i-1:i+1:2, i_eps)
         u0     = cc(i, n)
         u(1:2) = cc(i-1:i+1:2, n)

         cc(i, i_out) = sum(2 * idr2 * &
              a0*a(:)/(a0 + a(:)) * (u(:) - u0))
      end do
    end associate
  end subroutine box_vlpl

  !> Prolong from a parent to a child with index offset dix. This method could
  !> sometimes work better than the default prolongation, which does not take
  !> the variation in epsilon into account
!   subroutine vlpl_prolong(mg, p_id, dix, nc, iv, fine)
!     type(mg_t), intent(inout) :: mg
!     integer, intent(in)       :: p_id             !< Id of parent
!     integer, intent(in)       :: dix(1)        !< Offset of child in parent grid
!     integer, intent(in)       :: nc               !< Child grid size
!     integer, intent(in)       :: iv               !< Prolong from this variable
!     real(dp), intent(out)     :: fine(nc) !< Prolonged values

!     integer  :: i, hnc
! #if 1 == 2
!     integer  :: ic, jc
!     real(dp) :: f0, flx, fhx, fly, fhy
!     real(dp) :: c0, clx, chx, cly, chy
! #elif 1 == 3
!     integer  :: ic, jc, kc
!     real(dp) :: f0, flx, fhx, fly, fhy, flz, fhz
!     real(dp) :: c0, clx, chx, cly, chy, clz, chz
! #endif

!     hnc = nc/2

!     associate (crs => mg%boxes(p_id)%cc, &
!          i_eps => mg_iveps)
! #if 1 == 2
!       do j = 1, hnc
!          jc = j + dix(2)
!          do i = 1, hnc
!             ic = i + dix(1)

!             f0  = crs(ic, jc, iv)
!             flx = crs(ic-1, jc, iv)
!             fhx = crs(ic+1, jc, iv)
!             fly = crs(ic, jc-1, iv)
!             fhy = crs(ic, jc+1, iv)

!             c0  = 2 * crs(ic, jc, i_eps)
!             clx = crs(ic-1, jc, i_eps)
!             chx = crs(ic+1, jc, i_eps)
!             cly = crs(ic, jc-1, i_eps)
!             chy = crs(ic, jc+1, i_eps)

!             fine(2*i-1, 2*j-1) = (f0*c0 + flx*clx + fly*cly) / &
!                  (c0 + clx + cly)
!             fine(2*i  , 2*j-1) = (f0*c0 + fhx*chx + fly*cly) / &
!                  (c0 + chx + cly)
!             fine(2*i-1, 2*j)   = (f0*c0 + flx*clx + fhy*chy) / &
!                  (c0 + clx + chy)
!             fine(2*i  , 2*j)   = (f0*c0 + fhx*chx + fhy*chy) / &
!                  (c0 + chx + chy)
!          end do
!       end do
! #elif 1 == 3
!       do k = 1, hnc
!          kc = k + dix(3)
!          do j = 1, hnc
!             jc = j + dix(2)
!             do i = 1, hnc
!                ic = i + dix(1)

!                f0  = crs(ic, jc, kc, iv)
!                flx = crs(ic-1, jc, kc, iv)
!                fhx = crs(ic+1, jc, kc, iv)
!                fly = crs(ic, jc-1, kc, iv)
!                fhy = crs(ic, jc+1, kc, iv)
!                flz = crs(ic, jc, kc-1, iv)
!                fhz = crs(ic, jc, kc+1, iv)

!                c0  = crs(ic, jc, kc, i_eps)
!                clx = crs(ic-1, jc, kc, i_eps)
!                chx = crs(ic+1, jc, kc, i_eps)
!                cly = crs(ic, jc-1, kc, i_eps)
!                chy = crs(ic, jc+1, kc, i_eps)
!                clz = crs(ic, jc, kc-1, i_eps)
!                chz = crs(ic, jc, kc+1, i_eps)

!                fine(2*i-1, 2*j-1, 2*k-1) = (f0*c0 + flx*clx + fly*cly + flz*clz) / &
!                     (c0 + clx + cly + clz)
!                fine(2*i, 2*j-1, 2*k-1)   = (f0*c0 + fhx*chx + fly*cly + flz*clz) / &
!                     (c0 + chx + cly + clz)
!                fine(2*i-1, 2*j, 2*k-1)   = (f0*c0 + flx*clx + fhy*chy + flz*clz) / &
!                     (c0 + clx + chy + clz)
!                fine(2*i, 2*j, 2*k-1)     = (f0*c0 + fhx*chx + fhy*chy + flz*clz) / &
!                     (c0 + chx + chy + clz)
!                fine(2*i-1, 2*j-1, 2*k)   = (f0*c0 + flx*clx + fly*cly + fhz*chz) / &
!                     (c0 + clx + cly + chz)
!                fine(2*i, 2*j-1, 2*k)     = (f0*c0 + fhx*chx + fly*cly + fhz*chz) / &
!                     (c0 + chx + cly + chz)
!                fine(2*i-1, 2*j, 2*k)     = (f0*c0 + flx*clx + fhy*chy + fhz*chz) / &
!                     (c0 + clx + chy + chz)
!                fine(2*i, 2*j, 2*k)       = (f0*c0 + fhx*chx + fhy*chy + fhz*chz) / &
!                     (c0 + chx + chy + chz)
!             end do
!          end do
!       end do
! #endif
!     end associate
!   end subroutine vlpl_prolong

  !! File ../src/m_communication.f90

  !> Initialize MPI if needed, and store MPI information
  subroutine mg_comm_init(mg, comm)
    use mpi
    type(mg_t), intent(inout)     :: mg
    !> MPI communicator (default: MPI_COMM_WORLD)
    integer, intent(in), optional :: comm
    integer                       :: ierr
    logical                       :: initialized

    call mpi_initialized(initialized, ierr)
    if (.not. initialized) then
       call mpi_init(ierr)
    end if

    if (present(comm)) then
       mg%comm = comm
    else
       mg%comm = MPI_COMM_WORLD
    end if

    call mpi_comm_rank(mg%comm, mg%my_rank, ierr)
    call mpi_comm_size(mg%comm, mg%n_cpu, ierr)
  end subroutine mg_comm_init

  subroutine sort_and_transfer_buffers(mg, dsize)
    use mpi
    type(mg_t), intent(inout)    :: mg
    integer, intent(in)          :: dsize
    integer                      :: i, n_send, n_recv
    integer                      :: send_req(mg%n_cpu)
    integer                      :: recv_req(mg%n_cpu)
    integer                      :: ierr

    n_send = 0
    n_recv = 0

    do i = 0, mg%n_cpu - 1
       if (mg%buf(i)%i_send > 0) then
          n_send = n_send + 1
          call sort_sendbuf(mg%buf(i), dsize)
          call mpi_isend(mg%buf(i)%send, mg%buf(i)%i_send, MPI_DOUBLE, &
               i, 0, mg%comm, send_req(n_send), ierr)
       end if
       if (mg%buf(i)%i_recv > 0) then
          n_recv = n_recv + 1
          call mpi_irecv(mg%buf(i)%recv, mg%buf(i)%i_recv, MPI_DOUBLE, &
               i, 0, mg%comm, recv_req(n_recv), ierr)
       end if
    end do

    call mpi_waitall(n_recv, recv_req(1:n_recv), MPI_STATUSES_IGNORE, ierr)
    call mpi_waitall(n_send, send_req(1:n_send), MPI_STATUSES_IGNORE, ierr)

  end subroutine sort_and_transfer_buffers

  !> Sort send buffers according to the idbuf array
  subroutine sort_sendbuf(gc, dsize)

    type(mg_buf_t), intent(inout) :: gc
    integer, intent(in)           :: dsize !< Size of send buffer elements
    integer                       :: ix_sort(gc%i_ix)
    real(dp)                      :: buf_cpy(gc%i_send)
    integer                       :: i, j, n

    call mrgrnk(gc%ix(1:gc%i_ix), ix_sort)

    buf_cpy = gc%send(1:gc%i_send)

    do n = 1, gc%i_ix
       i = (n-1) * dsize
       j = (ix_sort(n)-1) * dsize
       gc%send(i+1:i+dsize) = buf_cpy(j+1:j+dsize)
    end do
    gc%ix(1:gc%i_ix) = gc%ix(ix_sort)

  end subroutine sort_sendbuf

  !! File ../src/m_ghost_cells.f90

  !> Specify minimum buffer size (per process) for communication
  subroutine mg_ghost_cell_buffer_size(mg, n_send, n_recv, dsize)
    type(mg_t), intent(inout) :: mg
    integer, intent(out)      :: n_send(0:mg%n_cpu-1)
    integer, intent(out)      :: n_recv(0:mg%n_cpu-1)
    integer, intent(out)      :: dsize
    integer                   :: i, id, lvl, nc

    allocate(mg%comm_ghostcell%n_send(0:mg%n_cpu-1, &
         mg%first_normal_lvl:mg%highest_lvl))
    allocate(mg%comm_ghostcell%n_recv(0:mg%n_cpu-1, &
         mg%first_normal_lvl:mg%highest_lvl))

    dsize = mg%box_size**(1-1)

    do lvl = mg%first_normal_lvl, mg%highest_lvl
       nc               = mg%box_size_lvl(lvl)
       mg%buf(:)%i_send = 0
       mg%buf(:)%i_recv = 0
       mg%buf(:)%i_ix   = 0

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call buffer_ghost_cells(mg, id, nc, 1, dry_run=.true.)
       end do

       if (lvl > 1) then
          do i = 1, size(mg%lvls(lvl-1)%my_ref_bnds)
             id = mg%lvls(lvl-1)%my_ref_bnds(i)
             call buffer_refinement_boundaries(mg, id, nc, 1, dry_run=.true.)
          end do
       end if

       ! Set ghost cells to received data
       mg%buf(:)%i_recv = 0
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call set_ghost_cells(mg, id, nc, 1, dry_run=.true.)
       end do

       mg%comm_ghostcell%n_send(:, lvl) = mg%buf(:)%i_send/dsize
       mg%comm_ghostcell%n_recv(:, lvl) = mg%buf(:)%i_recv/dsize
    end do

    n_send = maxval(mg%comm_ghostcell%n_send, dim=2)
    n_recv = maxval(mg%comm_ghostcell%n_recv, dim=2)
  end subroutine mg_ghost_cell_buffer_size

  !> Store boundary conditions for the solution variable, this can speed up
  !> calculations if the same boundary conditions are re-used.
  subroutine mg_phi_bc_store(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: lvl, nc

    nc = mg%box_size

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       call mg_phi_bc_store_lvl(mg, lvl, nc)
    end do

    mg%phi_bc_data_stored = .true.
  end subroutine mg_phi_bc_store

  subroutine mg_phi_bc_store_lvl(mg, lvl, nc)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer, intent(in)       :: nc
    real(dp)                  :: bc(1)
    integer                   :: i, id, nb, nb_id, bc_type

    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       do nb = 1, mg_num_neighbors
          nb_id = mg%boxes(id)%neighbors(nb)
          if (nb_id < mg_no_box) then
             ! Physical boundary
             if (associated(mg%bc(nb, mg_iphi)%boundary_cond)) then
                call mg%bc(nb, mg_iphi)%boundary_cond(mg%boxes(id), nc, &
                     mg_iphi, nb, bc_type, bc)
             else
                bc_type = mg%bc(nb, mg_iphi)%bc_type
                bc      = mg%bc(nb, mg_iphi)%bc_value
             end if

             ! Store the boundary condition type. This is not globally set in
             ! the tree, but all negative values are treated the same in
             ! other parts of the code
             mg%boxes(id)%neighbors(nb) = bc_type

             ! Store ghost cell data in the right-hand side
             call box_set_gc(mg%boxes(id), nb, nc, mg_irhs, bc)
          end if
       end do
    end do
  end subroutine mg_phi_bc_store_lvl

  !> Fill ghost cells at all grid levels
  subroutine mg_fill_ghost_cells(mg, iv)
    type(mg_t)          :: mg
    integer, intent(in) :: iv !< Index of variable
    integer             :: lvl

    do lvl = mg%lowest_lvl, mg%highest_lvl
       call mg_fill_ghost_cells_lvl(mg, lvl, iv)
    end do
  end subroutine mg_fill_ghost_cells

  !> Fill ghost cells at a grid level
  subroutine mg_fill_ghost_cells_lvl(mg, lvl, iv)

    type(mg_t)                   :: mg
    integer, intent(in)          :: lvl
    integer, intent(in)          :: iv !< Index of variable
    integer                      :: i, id, dsize, nc

    if (lvl < mg%lowest_lvl) &
         error stop "fill_ghost_cells_lvl: lvl < lowest_lvl"
    if (lvl > mg%highest_lvl) &
         error stop "fill_ghost_cells_lvl: lvl > highest_lvl"

    nc               = mg%box_size_lvl(lvl)

    if (lvl >= mg%first_normal_lvl) then
       dsize            = nc**(1-1)
       mg%buf(:)%i_send = 0
       mg%buf(:)%i_recv = 0
       mg%buf(:)%i_ix   = 0

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call buffer_ghost_cells(mg, id, nc, iv, .false.)
       end do

       if (lvl > 1) then
          do i = 1, size(mg%lvls(lvl-1)%my_ref_bnds)
             id = mg%lvls(lvl-1)%my_ref_bnds(i)
             call buffer_refinement_boundaries(mg, id, nc, iv, .false.)
          end do
       end if

       ! Transfer data between processes
       mg%buf(:)%i_recv = mg%comm_ghostcell%n_recv(:, lvl) * dsize
       call sort_and_transfer_buffers(mg, dsize)

       ! Set ghost cells to received data
       mg%buf(:)%i_recv = 0
    end if

    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call set_ghost_cells(mg, id, nc, iv, .false.)
    end do
  end subroutine mg_fill_ghost_cells_lvl

  subroutine buffer_ghost_cells(mg, id, nc, iv, dry_run)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    logical, intent(in)       :: dry_run
    integer                   :: nb, nb_id, nb_rank

    do nb = 1, mg_num_neighbors
       nb_id = mg%boxes(id)%neighbors(nb)

       if (nb_id > mg_no_box) then
          ! There is a neighbor
          nb_rank    = mg%boxes(nb_id)%rank

          if (nb_rank /= mg%my_rank) then
             call buffer_for_nb(mg, mg%boxes(id), nc, iv, nb_id, nb_rank, &
                  nb, dry_run)
          end if
       end if
    end do
  end subroutine buffer_ghost_cells

  subroutine buffer_refinement_boundaries(mg, id, nc, iv, dry_run)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    logical, intent(in)       :: dry_run
    integer                   :: nb, nb_id, c_ids(2**(1-1))
    integer                   :: n, c_id, c_rank

    do nb = 1, mg_num_neighbors
       nb_id = mg%boxes(id)%neighbors(nb)
       if (nb_id > mg_no_box) then
          if (mg_has_children(mg%boxes(nb_id))) then
             c_ids = mg%boxes(nb_id)%children(&
                  mg_child_adj_nb(:, mg_neighb_rev(nb)))

             do n = 1, mg_num_children/2
                c_id = c_ids(n)
                c_rank = mg%boxes(c_id)%rank

                if (c_rank /= mg%my_rank) then
                   ! Send all coarse ghost cells
                   call buffer_for_fine_nb(mg, mg%boxes(id), nc, iv, c_id, &
                        c_rank, nb, dry_run)
                end if
             end do
          end if
       end if
    end do
  end subroutine buffer_refinement_boundaries

  !> The routine that actually fills the ghost cells
  subroutine set_ghost_cells(mg, id, nc, iv, dry_run)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    logical, intent(in)       :: dry_run
    real(dp)                  :: bc(1)
    integer                   :: nb, nb_id, nb_rank, bc_type

    do nb = 1, mg_num_neighbors
       nb_id = mg%boxes(id)%neighbors(nb)

       if (nb_id > mg_no_box) then
          ! There is a neighbor
          nb_rank = mg%boxes(nb_id)%rank

          if (nb_rank /= mg%my_rank) then
             call fill_buffered_nb(mg, mg%boxes(id), nb_rank, &
                  nb, nc, iv, dry_run)
          else if (.not. dry_run) then
             call copy_from_nb(mg%boxes(id), mg%boxes(nb_id), &
                  nb, nc, iv)
          end if
       else if (nb_id == mg_no_box) then
          ! Refinement boundary
          call fill_refinement_bnd(mg, id, nb, nc, iv, dry_run)
       else if (.not. dry_run) then
          ! Physical boundary
          if (mg%phi_bc_data_stored .and. iv == mg_iphi) then
             ! Copy the boundary conditions stored in the ghost cells of the
             ! right-hand side
             call box_get_gc(mg%boxes(id), nb, nc, mg_irhs, bc)
             bc_type = nb_id
          else
             if (associated(mg%bc(nb, iv)%boundary_cond)) then
                call mg%bc(nb, iv)%boundary_cond(mg%boxes(id), nc, iv, &
                     nb, bc_type, bc)
             else
                bc_type = mg%bc(nb, iv)%bc_type
                bc = mg%bc(nb, iv)%bc_value
             end if
          end if

          call box_set_gc(mg%boxes(id), nb, nc, iv, bc)
          call bc_to_gc(mg, id, nc, iv, nb, bc_type)
       end if
    end do
  end subroutine set_ghost_cells

  subroutine fill_refinement_bnd(mg, id, nb, nc, iv, dry_run)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    integer, intent(in)       :: nb
    logical, intent(in)       :: dry_run
    real(dp)                  :: gc(1)
    integer                   :: p_id, p_nb_id, ix_offset(1)
    integer                   :: i, dsize, p_nb_rank

    dsize     = nc**(1-1)
    p_id      = mg%boxes(id)%parent
    p_nb_id   = mg%boxes(p_id)%neighbors(nb)
    p_nb_rank = mg%boxes(p_nb_id)%rank

    if (p_nb_rank /= mg%my_rank) then
       i = mg%buf(p_nb_rank)%i_recv
       if (.not. dry_run) then
          gc = reshape(mg%buf(p_nb_rank)%recv(i+1:i+dsize), shape(gc))
       end if
       mg%buf(p_nb_rank)%i_recv = mg%buf(p_nb_rank)%i_recv + dsize
    else if (.not. dry_run) then
       ix_offset = mg_get_child_offset(mg, id)
       call box_gc_for_fine_neighbor(mg%boxes(p_nb_id), mg_neighb_rev(nb), &
            ix_offset, nc, iv, gc)
    end if

    if (.not. dry_run) then
       if (associated(mg%bc(nb, iv)%refinement_bnd)) then
          call mg%bc(nb, iv)%refinement_bnd(mg%boxes(id), nc, iv, nb, gc)
       else
          call sides_rb(mg%boxes(id), nc, iv, nb, gc)
       end if
    end if
  end subroutine fill_refinement_bnd

  subroutine copy_from_nb(box, box_nb, nb, nc, iv)
    type(mg_box_t), intent(inout) :: box
    type(mg_box_t), intent(in)    :: box_nb
    integer, intent(in)           :: nb
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv
    real(dp)                      :: gc(1)

    call box_gc_for_neighbor(box_nb, mg_neighb_rev(nb), nc, iv, gc)
    call box_set_gc(box, nb, nc, iv, gc)
  end subroutine copy_from_nb

  subroutine buffer_for_nb(mg, box, nc, iv, nb_id, nb_rank, nb, dry_run)
    type(mg_t), intent(inout)  :: mg
    type(mg_box_t), intent(inout) :: box
    integer, intent(in)        :: nc
    integer, intent(in)        :: iv
    integer, intent(in)        :: nb_id
    integer, intent(in)        :: nb_rank
    integer, intent(in)        :: nb
    logical, intent(in)        :: dry_run
    integer                    :: i, dsize
    real(dp)                   :: gc(1)

    i     = mg%buf(nb_rank)%i_send
    dsize = nc**(1-1)

    if (.not. dry_run) then
       call box_gc_for_neighbor(box, nb, nc, iv, gc)
       mg%buf(nb_rank)%send(i+1:i+dsize) = pack(gc, .true.)
    end if

    ! Later the buffer is sorted, using the fact that loops go from low to high
    ! box id, and we fill ghost cells according to the neighbor order
    i = mg%buf(nb_rank)%i_ix
    if (.not. dry_run) then
       mg%buf(nb_rank)%ix(i+1) = mg_num_neighbors * nb_id + mg_neighb_rev(nb)
    end if

    mg%buf(nb_rank)%i_send = mg%buf(nb_rank)%i_send + dsize
    mg%buf(nb_rank)%i_ix   = mg%buf(nb_rank)%i_ix + 1
  end subroutine buffer_for_nb

  subroutine buffer_for_fine_nb(mg, box, nc, iv, fine_id, fine_rank, nb, dry_run)
    type(mg_t), intent(inout)  :: mg
    type(mg_box_t), intent(inout) :: box
    integer, intent(in)        :: nc
    integer, intent(in)        :: iv
    integer, intent(in)        :: fine_id
    integer, intent(in)        :: fine_rank
    integer, intent(in)        :: nb
    logical, intent(in)        :: dry_run
    integer                    :: i, dsize, ix_offset(1)
    real(dp)                   :: gc(1)

    i     = mg%buf(fine_rank)%i_send
    dsize = nc**(1-1)

    if (.not. dry_run) then
       ix_offset = mg_get_child_offset(mg, fine_id)
       call box_gc_for_fine_neighbor(box, nb, ix_offset, nc, iv, gc)
       mg%buf(fine_rank)%send(i+1:i+dsize) = pack(gc, .true.)
    end if

    ! Later the buffer is sorted, using the fact that loops go from low to high
    ! box id, and we fill ghost cells according to the neighbor order
    i = mg%buf(fine_rank)%i_ix
    if (.not. dry_run) then
       mg%buf(fine_rank)%ix(i+1) = mg_num_neighbors * fine_id + &
            mg_neighb_rev(nb)
    end if

    mg%buf(fine_rank)%i_send = mg%buf(fine_rank)%i_send + dsize
    mg%buf(fine_rank)%i_ix   = mg%buf(fine_rank)%i_ix + 1
  end subroutine buffer_for_fine_nb

  subroutine fill_buffered_nb(mg, box, nb_rank, nb, nc, iv, dry_run)
    type(mg_t), intent(inout)  :: mg
    type(mg_box_t), intent(inout) :: box
    integer, intent(in)        :: nb_rank
    integer, intent(in)        :: nb
    integer, intent(in)        :: nc
    integer, intent(in)        :: iv
    logical, intent(in)        :: dry_run
    integer                    :: i, dsize
    real(dp)                   :: gc(1)

    i     = mg%buf(nb_rank)%i_recv
    dsize = nc**(1-1)

    if (.not. dry_run) then
       gc = mg%buf(nb_rank)%recv(i+1)
       call box_set_gc(box, nb, nc, iv, gc)
    end if
    mg%buf(nb_rank)%i_recv = mg%buf(nb_rank)%i_recv + dsize

  end subroutine fill_buffered_nb

  subroutine box_gc_for_neighbor(box, nb, nc, iv, gc)
    type(mg_box_t), intent(in) :: box
    integer, intent(in)     :: nb, nc, iv
    real(dp), intent(out)   :: gc(1)

    select case (nb)
    case (mg_neighb_lowx)
       gc = box%cc(1, iv)
    case (mg_neighb_highx)
       gc = box%cc(nc, iv)
    end select
  end subroutine box_gc_for_neighbor

  !> Get ghost cells for a fine neighbor
  subroutine box_gc_for_fine_neighbor(box, nb, di, nc, iv, gc)
    type(mg_box_t), intent(in) :: box
    integer, intent(in)     :: nb       !< Direction of fine neighbor
    integer, intent(in)     :: di(1) !< Index offset of fine neighbor
    integer, intent(in)     :: nc, iv
    real(dp), intent(out)   :: gc(1)
    real(dp)                :: tmp
    integer                 :: hnc

    hnc = nc/2

    ! First fill a temporary array with data next to the fine grid
    select case (nb)
    case (mg_neighb_lowx)
       tmp = box%cc(1, iv)
    case (mg_neighb_highx)
       tmp = box%cc(nc, iv)
    case default
       error stop
    end select

    ! Now interpolate the coarse grid data to obtain values 'straight' next to
    ! the fine grid points
    gc = tmp
  end subroutine box_gc_for_fine_neighbor

  subroutine box_get_gc(box, nb, nc, iv, gc)
    type(mg_box_t), intent(in) :: box
    integer, intent(in)        :: nb, nc, iv
    real(dp), intent(out)       :: gc(1)

    select case (nb)
    case (mg_neighb_lowx)
       gc = box%cc(0, iv)
    case (mg_neighb_highx)
       gc = box%cc(nc+1, iv)
    end select
  end subroutine box_get_gc

  subroutine box_set_gc(box, nb, nc, iv, gc)
    type(mg_box_t), intent(inout) :: box
    integer, intent(in)        :: nb, nc, iv
    real(dp), intent(in)       :: gc(1)

    select case (nb)
    case (mg_neighb_lowx)
       box%cc(0, iv)    = gc(1)
    case (mg_neighb_highx)
       box%cc(nc+1, iv) = gc(1)
    end select
  end subroutine box_set_gc

  subroutine bc_to_gc(mg, id, nc, iv, nb, bc_type)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    integer, intent(in)       :: nb      !< Neighbor direction
    integer, intent(in)       :: bc_type !< Type of b.c.
    real(dp)                  :: c0, c1, c2, dr

    ! If we call the interior point x1, x2 and the ghost point x0, then a
    ! Dirichlet boundary value b can be imposed as:
    ! x0 = -x1 + 2*b
    ! A Neumann b.c. can be imposed as:
    ! x0 = x1 +/- dx * b
    ! A continuous boundary (same slope) as:
    ! x0 = 2 * x1 - x2
    ! Below, we set coefficients to handle these cases
    select case (bc_type)
    case (mg_bc_dirichlet)
       c0 = 2
       c1 = -1
       c2 = 0
    case (mg_bc_neumann)
       dr = mg%dr(mg_neighb_dim(nb), mg%boxes(id)%lvl)
       c0 = dr * mg_neighb_high_pm(nb) ! This gives a + or - sign
       c1 = 1
       c2 = 0
    case (mg_bc_continuous)
       c0 = 0
       c1 = 2
       c2 = -1
    case default
       error stop "bc_to_gc: unknown boundary condition"
    end select

    select case (nb)
    case (mg_neighb_lowx)
       mg%boxes(id)%cc(0, iv) = &
            c0 * mg%boxes(id)%cc(0, iv) + &
            c1 * mg%boxes(id)%cc(1, iv) + &
            c2 * mg%boxes(id)%cc(2, iv)
    case (mg_neighb_highx)
       mg%boxes(id)%cc(nc+1, iv) = &
            c0 * mg%boxes(id)%cc(nc+1, iv) + &
            c1 * mg%boxes(id)%cc(nc, iv) + &
            c2 * mg%boxes(id)%cc(nc-1, iv)
    end select
  end subroutine bc_to_gc

  !> Fill ghost cells near refinement boundaries which preserves diffusive fluxes.
  subroutine sides_rb(box, nc, iv, nb, gc)
    type(mg_box_t), intent(inout) :: box
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    integer, intent(in)       :: nb !< Ghost cell direction
    !> Interpolated coarse grid ghost cell data (but not yet in the nb direction)
    real(dp), intent(in)      :: gc(1)
    integer                   :: di
    integer                   :: i, ix, dix

    if (mg_neighb_low(nb)) then
       ix = 1
       dix = 1
    else
       ix = nc
       dix = -1
    end if

    select case (mg_neighb_dim(nb))
    case (1)
       i = ix
       di = dix
       box%cc(i-di, iv) = (2 * gc(1) + box%cc(i, iv))/3.0_dp
    end select

  end subroutine sides_rb

  !! File ../src/m_prolong.f90

  !> Specify minimum buffer size (per process) for communication
  subroutine mg_prolong_buffer_size(mg, n_send, n_recv, dsize)
    type(mg_t), intent(inout) :: mg
    integer, intent(out)      :: n_send(0:mg%n_cpu-1)
    integer, intent(out)      :: n_recv(0:mg%n_cpu-1)
    integer, intent(out)      :: dsize
    integer                   :: lvl, min_lvl

    if (.not. allocated(mg%comm_restrict%n_send)) then
       error stop "Call restrict_buffer_size before prolong_buffer_size"
    end if

    min_lvl = max(mg%first_normal_lvl-1, mg%lowest_lvl)
    allocate(mg%comm_prolong%n_send(0:mg%n_cpu-1, &
         min_lvl:mg%highest_lvl))
    allocate(mg%comm_prolong%n_recv(0:mg%n_cpu-1, &
         min_lvl:mg%highest_lvl))

    mg%comm_prolong%n_recv(:, mg%highest_lvl) = 0
    mg%comm_prolong%n_send(:, mg%highest_lvl) = 0

    do lvl = min_lvl, mg%highest_lvl-1
       mg%comm_prolong%n_recv(:, lvl) = &
            mg%comm_restrict%n_send(:, lvl+1)
       mg%comm_prolong%n_send(:, lvl) = &
            mg%comm_restrict%n_recv(:, lvl+1)
    end do

    ! Send fine grid points, because this is more flexible than sending coarse
    ! grid points (e.g., when multiple variables are used for interpolation)
    dsize = (mg%box_size)**1
    n_send = maxval(mg%comm_prolong%n_send, dim=2)
    n_recv = maxval(mg%comm_prolong%n_recv, dim=2)
  end subroutine mg_prolong_buffer_size

  !> Prolong variable iv from lvl to variable iv_to at lvl+1
  subroutine mg_prolong(mg, lvl, iv, iv_to, method, add)

    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl   !< Level to prolong from
    integer, intent(in)       :: iv    !< Source variable
    integer, intent(in)       :: iv_to !< Target variable
    procedure(mg_box_prolong) :: method !< Prolongation method
    logical, intent(in)       :: add   !< If true, add to current values
    integer                   :: i, id, dsize, nc

    if (lvl == mg%highest_lvl) error stop "cannot prolong highest level"
    if (lvl < mg%lowest_lvl) error stop "cannot prolong below lowest level"

    ! Below the first normal level, all boxes are on the same CPU
    if (lvl >= mg%first_normal_lvl-1) then
       dsize            = mg%box_size**1
       mg%buf(:)%i_send = 0
       mg%buf(:)%i_ix   = 0

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call prolong_set_buffer(mg, id, mg%box_size, iv, method)
       end do

       mg%buf(:)%i_recv = mg%comm_prolong%n_recv(:, lvl) * dsize
       call sort_and_transfer_buffers(mg, dsize)
       mg%buf(:)%i_recv = 0
    end if

    nc = mg%box_size_lvl(lvl+1)
    do i = 1, size(mg%lvls(lvl+1)%my_ids)
       id = mg%lvls(lvl+1)%my_ids(i)
       call prolong_onto(mg, id, nc, iv, iv_to, add, method)
    end do
  end subroutine mg_prolong

  !> In case the fine grid is on a different CPU, perform the prolongation and
  !> store the fine-grid values in the send buffer.
  !>
  !> @todo Check whether it's faster to send coarse data and prolong afterwards
  subroutine prolong_set_buffer(mg, id, nc, iv, method)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    procedure(mg_box_prolong) :: method
    integer                   :: i, dix(1)
    integer                   :: i_c, c_id, c_rank, dsize
    real(dp)                  :: tmp(nc)

    dsize = nc**1

    do i_c = 1, mg_num_children
       c_id = mg%boxes(id)%children(i_c)
       if (c_id > mg_no_box) then
          c_rank = mg%boxes(c_id)%rank
          if (c_rank /= mg%my_rank) then
             dix = mg_get_child_offset(mg, c_id)
             call method(mg, id, dix, nc, iv, tmp)

             i   = mg%buf(c_rank)%i_send
             mg%buf(c_rank)%send(i+1:i+dsize) = pack(tmp, .true.)
             mg%buf(c_rank)%i_send  = mg%buf(c_rank)%i_send + dsize

             i                      = mg%buf(c_rank)%i_ix
             mg%buf(c_rank)%ix(i+1) = c_id
             mg%buf(c_rank)%i_ix    = mg%buf(c_rank)%i_ix + 1
          end if
       end if
    end do
  end subroutine prolong_set_buffer

  !> Prolong onto a child box
  subroutine prolong_onto(mg, id, nc, iv, iv_to, add, method)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv    !< Prolong from this variable
    integer, intent(in)       :: iv_to !< Prolong to this variable
    logical, intent(in)       :: add   !< If true, add to current values
    procedure(mg_box_prolong) :: method
    integer                   :: hnc, p_id, p_rank, i, dix(1), dsize
    real(dp)                  :: tmp(nc)

    hnc    = nc/2
    p_id   = mg%boxes(id)%parent
    p_rank = mg%boxes(p_id)%rank

    if (p_rank == mg%my_rank) then
       dix    = mg_get_child_offset(mg, id)
       call method(mg, p_id, dix, nc, iv, tmp)
    else
       dsize  = nc**1
       i = mg%buf(p_rank)%i_recv
       tmp = reshape(mg%buf(p_rank)%recv(i+1:i+dsize), [nc])
       mg%buf(p_rank)%i_recv = mg%buf(p_rank)%i_recv + dsize
    end if

    if (add) then
       mg%boxes(id)%cc(1:nc, iv_to) = &
            mg%boxes(id)%cc(1:nc, iv_to) + tmp
    else
       mg%boxes(id)%cc(1:nc, iv_to) = tmp
    end if

  end subroutine prolong_onto

  !> Prolong from a parent to a child with index offset dix
  subroutine mg_prolong_sparse(mg, p_id, dix, nc, iv, fine)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: p_id             !< Id of parent
    integer, intent(in)       :: dix(1)        !< Offset of child in parent grid
    integer, intent(in)       :: nc               !< Child grid size
    integer, intent(in)       :: iv               !< Prolong from this variable
    real(dp), intent(out)     :: fine(nc) !< Prolonged values

    integer  :: i, hnc
    integer  :: ic
    real(dp) :: f0, flx, fhx

    hnc = nc/2

    associate (crs => mg%boxes(p_id)%cc)
      do i = 1, hnc
         ic = i + dix(1)

         f0  = 0.75_dp * crs(ic, iv)
         flx = 0.25_dp * crs(ic-1, iv)
         fhx = 0.25_dp * crs(ic+1, iv)

         fine(2*i-1)   = f0 + flx
         fine(2*i  )   = f0 + fhx
      end do
    end associate
  end subroutine mg_prolong_sparse

  !! File ../src/m_helmholtz.f90

  subroutine helmholtz_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    mg%subtract_mean = .false.

    select case (mg%geometry_type)
    case (mg_cartesian)
       mg%box_op => box_helmh

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_helmh
       case default
          error stop "helmholtz_set_methods: unsupported smoother type"
       end select
    case default
       error stop "helmholtz_set_methods: unsupported geometry"
    end select

  end subroutine helmholtz_set_methods

  subroutine helmholtz_set_lambda(lambda)
    real(dp), intent(in) :: lambda

    if (lambda < 0) &
         error stop "helmholtz_set_lambda: lambda < 0 not allowed"

    helmholtz_lambda = lambda
  end subroutine helmholtz_set_lambda

  !> Perform Gauss-Seidel relaxation on box for a Helmholtz operator
  subroutine box_gs_helmh(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: i, i0, di
    real(dp)                  :: idr2(1), fac
    logical                   :: redblack

    idr2 = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    fac = 1.0_dp / (2 * sum(idr2) + helmholtz_lambda)
    i0  = 1

    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
      if (redblack) i0 = 2 - iand(redblack_cntr, 1)

         do i = i0, nc, di
            cc(i, n) = fac * ( &
                 idr2(1) * (cc(i+1, n) + cc(i-1, n)) - &
                 cc(i, mg_irhs))
         end do
    end associate
  end subroutine box_gs_helmh

  !> Perform Helmholtz operator on a box
  subroutine box_helmh(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Helmholtz in
    integer                   :: i
    real(dp)                  :: idr2(1)

    idr2 = 1 / mg%dr(:, mg%boxes(id)%lvl)**2

    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
      do i = 1, nc
         cc(i, i_out) = &
              idr2(1) * (cc(i-1, n) + cc(i+1, n) - 2 * cc(i, n)) - &
              helmholtz_lambda * cc(i, n)
      end do
    end associate
  end subroutine box_helmh

  !! File ../src/m_multigrid.f90

  subroutine mg_set_methods(mg)

    type(mg_t), intent(inout) :: mg

    ! Set default prolongation method (routines below can override this)
    mg%box_prolong => mg_prolong_sparse

    select case (mg%operator_type)
    case (mg_laplacian)
       call laplacian_set_methods(mg)
    case (mg_vlaplacian)
       call vlaplacian_set_methods(mg)
    case (mg_helmholtz)
       call helmholtz_set_methods(mg)
    case (mg_vhelmholtz)
       call vhelmholtz_set_methods(mg)
    case default
       error stop "mg_set_methods: unknown operator"
    end select

    ! For red-black, perform two smoothing sub-steps so that all unknowns are
    ! updated per cycle
    if (mg%smoother_type == mg_smoother_gsrb) then
       mg%n_smoother_substeps = 2
    else
       mg%n_smoother_substeps = 1
    end if
  end subroutine mg_set_methods

  subroutine check_methods(mg)
    type(mg_t), intent(inout) :: mg

    if (.not. associated(mg%box_op) .or. &
         .not. associated(mg%box_smoother)) then
       call mg_set_methods(mg)
    end if

  end subroutine check_methods

  subroutine mg_add_timers(mg)
    type(mg_t), intent(inout) :: mg
    timer_total_vcycle  = mg_add_timer(mg, "mg total V-cycle")
    timer_total_fmg     = mg_add_timer(mg, "mg total FMG cycle")
    timer_smoother      = mg_add_timer(mg, "mg smoother")
    timer_smoother_gc   = mg_add_timer(mg, "mg smoother g.c.")
    timer_coarse        = mg_add_timer(mg, "mg coarse")
    timer_correct       = mg_add_timer(mg, "mg correct")
    timer_update_coarse = mg_add_timer(mg, "mg update coarse")
  end subroutine mg_add_timers

  !> Perform FAS-FMG cycle (full approximation scheme, full multigrid).
  subroutine mg_fas_fmg(mg, have_guess, max_res)
    type(mg_t), intent(inout)       :: mg
    logical, intent(in)             :: have_guess !< If false, start from phi = 0
    real(dp), intent(out), optional :: max_res    !< Store max(abs(residual))
    integer                         :: lvl, i, id

    call check_methods(mg)
    if (timer_smoother == -1) call mg_add_timers(mg)

    call mg_timer_start(mg%timers(timer_total_fmg))

    if (.not. have_guess) then
       do lvl = mg%highest_lvl, mg%lowest_lvl, -1
          do i = 1, size(mg%lvls(lvl)%my_ids)
             id = mg%lvls(lvl)%my_ids(i)
             mg%boxes(id)%cc(:, mg_iphi) = 0.0_dp
          end do
       end do
    end if

    ! Ensure ghost cells are filled correctly
    call mg_fill_ghost_cells_lvl(mg, mg%highest_lvl, mg_iphi)

    do lvl = mg%highest_lvl,  mg%lowest_lvl+1, -1
       ! Set rhs on coarse grid and restrict phi
       call mg_timer_start(mg%timers(timer_update_coarse))
       call update_coarse(mg, lvl)
       call mg_timer_end(mg%timers(timer_update_coarse))
    end do

    if (mg%subtract_mean) then
       ! For fully periodic solutions, the mean source term has to be zero
       call subtract_mean(mg, mg_irhs, .false.)
    end if

    do lvl = mg%lowest_lvl, mg%highest_lvl
       ! Store phi_old
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          mg%boxes(id)%cc(:, mg_iold) = &
            mg%boxes(id)%cc(:, mg_iphi)
       end do

       if (lvl > mg%lowest_lvl) then
          ! Correct solution at this lvl using lvl-1 data
          ! phi = phi + prolong(phi_coarse - phi_old_coarse)
          call mg_timer_start(mg%timers(timer_correct))
          call correct_children(mg, lvl-1)
          call mg_timer_end(mg%timers(timer_correct))

          ! Update ghost cells
          call mg_fill_ghost_cells_lvl(mg, lvl, mg_iphi)
       end if

       ! Perform V-cycle, possibly set residual on last iteration
       if (lvl == mg%highest_lvl) then
          call mg_fas_vcycle(mg, lvl, max_res, standalone=.false.)
       else
          call mg_fas_vcycle(mg, lvl, standalone=.false.)
       end if
    end do

    call mg_timer_end(mg%timers(timer_total_fmg))
  end subroutine mg_fas_fmg

  !> Perform FAS V-cycle (full approximation scheme).
  subroutine mg_fas_vcycle(mg, highest_lvl, max_res, standalone)
    use mpi
    type(mg_t), intent(inout)       :: mg
    integer, intent(in), optional   :: highest_lvl !< Maximum level for V-cycle
    real(dp), intent(out), optional :: max_res     !< Store max(abs(residual))
    !> Whether the V-cycle is called by itself (default: true)
    logical, intent(in), optional   :: standalone
    integer                         :: lvl, min_lvl, i, max_lvl, ierr
    real(dp)                        :: init_res, res
    logical                         :: is_standalone

    is_standalone = .true.
    if (present(standalone)) is_standalone = standalone

    call check_methods(mg)
    if (timer_smoother == -1) call mg_add_timers(mg)

    if (is_standalone) &
         call mg_timer_start(mg%timers(timer_total_vcycle))

    if (mg%subtract_mean .and. .not. present(highest_lvl)) then
       ! Assume that this is a stand-alone call. For fully periodic solutions,
       ! ensure the mean source term is zero.
       call subtract_mean(mg, mg_irhs, .false.)
    end if

    min_lvl = mg%lowest_lvl
    max_lvl = mg%highest_lvl
    if (present(highest_lvl)) max_lvl = highest_lvl

    ! Ensure ghost cells are filled correctly
    if (is_standalone) then
       call mg_fill_ghost_cells_lvl(mg, max_lvl, mg_iphi)
    end if

    do lvl = max_lvl,  min_lvl+1, -1
       ! Downwards relaxation
       call smooth_boxes(mg, lvl, mg%n_cycle_down)

       ! Set rhs on coarse grid, restrict phi, and copy mg_iphi to mg_iold for the
       ! correction later
       call mg_timer_start(mg%timers(timer_update_coarse))
       call update_coarse(mg, lvl)
       call mg_timer_end(mg%timers(timer_update_coarse))
    end do

    call mg_timer_start(mg%timers(timer_coarse))
    if (.not. all(mg%boxes(mg%lvls(min_lvl)%ids)%rank == &
         mg%boxes(mg%lvls(min_lvl)%ids(1))%rank)) then
       error stop "Multiple CPUs for coarse grid (not implemented yet)"
    end if

    init_res = max_residual_lvl(mg, min_lvl)
    do i = 1, mg%max_coarse_cycles
       call smooth_boxes(mg, min_lvl, mg%n_cycle_up+mg%n_cycle_down)
       res = max_residual_lvl(mg, min_lvl)
       if (res < mg%residual_coarse_rel * init_res .or. &
            res < mg%residual_coarse_abs) exit
    end do
    call mg_timer_end(mg%timers(timer_coarse))

    ! Do the upwards part of the v-cycle in the tree
    do lvl = min_lvl+1, max_lvl
       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       call mg_timer_start(mg%timers(timer_correct))
       call correct_children(mg, lvl-1)

       ! Have to fill ghost cells after correction
       call mg_fill_ghost_cells_lvl(mg, lvl, mg_iphi)
       call mg_timer_end(mg%timers(timer_correct))

       ! Upwards relaxation
       call smooth_boxes(mg, lvl, mg%n_cycle_up)
    end do

    if (present(max_res)) then
       init_res = 0.0_dp
       do lvl = min_lvl, max_lvl
          res = max_residual_lvl(mg, lvl)
          init_res = max(res, init_res)
       end do
       call mpi_allreduce(init_res, max_res, 1, &
            mpi_double, mpi_max, mg%comm, ierr)
    end if

    ! Subtract mean(phi) from phi
    if (mg%subtract_mean) then
       call subtract_mean(mg, mg_iphi, .true.)
    end if

    if (is_standalone) &
         call mg_timer_end(mg%timers(timer_total_vcycle))
  end subroutine mg_fas_vcycle

  subroutine subtract_mean(mg, iv, include_ghostcells)
    use mpi
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: iv
    logical, intent(in)       :: include_ghostcells
    integer                   :: i, id, lvl, nc, ierr
    real(dp)                  :: sum_iv, mean_iv, volume

    nc = mg%box_size
    sum_iv = get_sum(mg, iv)
    call mpi_allreduce(sum_iv, mean_iv, 1, &
         mpi_double, mpi_sum, mg%comm, ierr)

    ! Divide by total grid volume to get mean
    volume = nc**1 * product(mg%dr(:, 1)) * size(mg%lvls(1)%ids)
    mean_iv = mean_iv / volume

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          if (include_ghostcells) then
             mg%boxes(id)%cc(:, iv) = &
                  mg%boxes(id)%cc(:, iv) - mean_iv
          else
             mg%boxes(id)%cc(1:nc, iv) = &
                  mg%boxes(id)%cc(1:nc, iv) - mean_iv
          end if
       end do
    end do
  end subroutine subtract_mean

  real(dp) function get_sum(mg, iv)
    type(mg_t), intent(in) :: mg
    integer, intent(in)    :: iv
    integer                :: lvl, i, id, nc
    real(dp)               :: w

    get_sum = 0.0_dp
    do lvl = 1, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       w  = product(mg%dr(:, lvl)) ! Adjust for non-Cartesian cases
       do i = 1, size(mg%lvls(lvl)%my_leaves)
          id = mg%lvls(lvl)%my_leaves(i)
          get_sum = get_sum + w * &
               sum(mg%boxes(id)%cc(1:nc, iv))
       end do
    end do
  end function get_sum

  real(dp) function max_residual_lvl(mg, lvl)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer                   :: i, id, nc
    real(dp)                  :: res

    nc           = mg%box_size_lvl(lvl)
    max_residual_lvl = 0.0_dp

    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call residual_box(mg, id, nc)
       res = maxval(abs(mg%boxes(id)%cc(1:nc, mg_ires)))
       max_residual_lvl = max(max_residual_lvl, res)
    end do
  end function max_residual_lvl

  !   subroutine solve_coarse_grid(mg)
  !     use m_fishpack
  !     type(mg_t), intent(inout) :: mg

  !     real(dp) :: rhs(mg%box_size)
  !     real(dp) :: rmin(1), rmax(1)
  !     integer  :: nc, nx(1), my_boxes, total_boxes

  !     my_boxes    = size(mg%lvls(1)%my_ids)
  !     total_boxes = size(mg%lvls(1)%ids)
  !     nc          = mg%box_size

  !     if (my_boxes == total_boxes) then
  !        nx(:) = nc
  !        rmin  = [0.0_dp]
  !        rmax  = mg%dr(1) * [nc]
  !        rhs   = mg%boxes(1)%cc(1:nc, mg_irhs)

  ! #if 1 == 2
  !        call fishpack_2d(nx, rhs, mg%bc, rmin, rmax)
  ! #elif 1 == 3
  !        call fishpack_3d(nx, rhs, mg%bc, rmin, rmax)
  ! #endif

  !        mg%boxes(1)%cc(1:nc, mg_iphi) = rhs
  !     else if (my_boxes > 0) then
  !        error stop "Boxes at level 1 at different processors"
  !     end if

  !     call fill_ghost_cells_lvl(mg, 1)
  !   end subroutine solve_coarse_grid

  ! Set rhs on coarse grid, restrict phi, and copy mg_iphi to mg_iold for the
  ! correction later
  subroutine update_coarse(mg, lvl)
    type(mg_t), intent(inout) :: mg     !< Tree containing full grid
    integer, intent(in)       :: lvl !< Update coarse values at lvl-1
    integer                   :: i, id, nc, nc_c

    nc   = mg%box_size_lvl(lvl)
    nc_c = mg%box_size_lvl(lvl-1)

    ! Compute residual
    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call residual_box(mg, id, nc)
    end do

    ! Restrict phi and the residual
    call mg_restrict_lvl(mg, mg_iphi, lvl)
    call mg_restrict_lvl(mg, mg_ires, lvl)

    call mg_fill_ghost_cells_lvl(mg, lvl-1, mg_iphi)

    ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined, and
    ! store current coarse phi in old.
    do i = 1, size(mg%lvls(lvl-1)%my_parents)
       id = mg%lvls(lvl-1)%my_parents(i)

       ! Set rhs = L phi
       call mg%box_op(mg, id, nc_c, mg_irhs)

       ! Add the fine grid residual to rhs
       mg%boxes(id)%cc(1:nc_c, mg_irhs) = &
            mg%boxes(id)%cc(1:nc_c, mg_irhs) + &
            mg%boxes(id)%cc(1:nc_c, mg_ires)

       ! Story a copy of phi
       mg%boxes(id)%cc(:, mg_iold) = &
            mg%boxes(id)%cc(:, mg_iphi)
    end do
  end subroutine update_coarse

  ! Sets phi = phi + prolong(phi_coarse - phi_old_coarse)
  subroutine correct_children(mg, lvl)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer                   :: i, id

    do i = 1, size(mg%lvls(lvl)%my_parents)
       id = mg%lvls(lvl)%my_parents(i)

       ! Store the correction in mg_ires
       mg%boxes(id)%cc(:, mg_ires) = &
            mg%boxes(id)%cc(:, mg_iphi) - &
            mg%boxes(id)%cc(:, mg_iold)
    end do

    call mg_prolong(mg, lvl, mg_ires, mg_iphi, mg%box_prolong, add=.true.)
  end subroutine correct_children

  subroutine smooth_boxes(mg, lvl, n_cycle)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer, intent(in)       :: n_cycle !< Number of cycles to perform
    integer                   :: n, i, id, nc

    nc = mg%box_size_lvl(lvl)

    do n = 1, n_cycle * mg%n_smoother_substeps
       call mg_timer_start(mg%timers(timer_smoother))
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call mg%box_smoother(mg, id, nc, n)
       end do
       call mg_timer_end(mg%timers(timer_smoother))

       call mg_timer_start(mg%timers(timer_smoother_gc))
       call mg_fill_ghost_cells_lvl(mg, lvl, mg_iphi)
       call mg_timer_end(mg%timers(timer_smoother_gc))
    end do
  end subroutine smooth_boxes

  subroutine residual_box(mg, id, nc)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc

    call mg%box_op(mg, id, nc, mg_ires)

    mg%boxes(id)%cc(1:nc, mg_ires) = &
         mg%boxes(id)%cc(1:nc, mg_irhs) &
         - mg%boxes(id)%cc(1:nc, mg_ires)
  end subroutine residual_box

  !> Apply operator to the tree and store in variable i_out
  subroutine mg_apply_op(mg, i_out, op)
    type(mg_t), intent(inout)      :: mg
    integer, intent(in)            :: i_out
    procedure(mg_box_op), optional :: op
    integer                        :: lvl, i, id, nc

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          if (present(op)) then
             call op(mg, id, nc, i_out)
          else
             call mg%box_op(mg, id, nc, i_out)
          end if
       end do
    end do
  end subroutine mg_apply_op

  !! File ../src/m_restrict.f90

  !> Specify minimum buffer size (per process) for communication
  subroutine mg_restrict_buffer_size(mg, n_send, n_recv, dsize)
    type(mg_t), intent(inout) :: mg
    integer, intent(out)      :: n_send(0:mg%n_cpu-1)
    integer, intent(out)      :: n_recv(0:mg%n_cpu-1)
    integer, intent(out)      :: dsize
    integer :: n_out(0:mg%n_cpu-1, mg%first_normal_lvl:mg%highest_lvl)
    integer :: n_in(0:mg%n_cpu-1, mg%first_normal_lvl:mg%highest_lvl)
    integer                   :: lvl, i, id, p_id, p_rank
    integer                   :: i_c, c_id, c_rank, min_lvl

    n_out(:, :) = 0
    n_in(:, :)  = 0
    min_lvl = max(mg%lowest_lvl+1, mg%first_normal_lvl)

    do lvl = min_lvl, mg%highest_lvl
       ! Number of messages to receive (at lvl-1)
       do i = 1, size(mg%lvls(lvl-1)%my_parents)
          id = mg%lvls(lvl-1)%my_parents(i)
          do i_c = 1, mg_num_children
             c_id = mg%boxes(id)%children(i_c)

             if (c_id > mg_no_box) then
                c_rank = mg%boxes(c_id)%rank
                if (c_rank /= mg%my_rank) then
                   n_in(c_rank, lvl) = n_in(c_rank, lvl) + 1
                end if
             end if
          end do
       end do

       ! Number of messages to send
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)

          p_id = mg%boxes(id)%parent
          p_rank = mg%boxes(p_id)%rank
          if (p_rank /= mg%my_rank) then
             n_out(p_rank, lvl) = n_out(p_rank, lvl) + 1
          end if

       end do
    end do

    allocate(mg%comm_restrict%n_send(0:mg%n_cpu-1, &
         mg%first_normal_lvl:mg%highest_lvl))
    allocate(mg%comm_restrict%n_recv(0:mg%n_cpu-1, &
         mg%first_normal_lvl:mg%highest_lvl))
    mg%comm_restrict%n_send = n_out
    mg%comm_restrict%n_recv = n_in

    dsize  = (mg%box_size/2)**1
    n_send = maxval(n_out, dim=2)
    n_recv = maxval(n_in, dim=2)
  end subroutine mg_restrict_buffer_size

  !> Restrict all levels
  subroutine mg_restrict(mg, iv)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: iv
    integer                   :: lvl

    do lvl = mg%highest_lvl, mg%lowest_lvl+1, -1
       call mg_restrict_lvl(mg, iv, lvl)
    end do
  end subroutine mg_restrict

  !> Restrict from lvl to lvl-1
  subroutine mg_restrict_lvl(mg, iv, lvl)

    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: iv
    integer, intent(in)       :: lvl
    integer                   :: i, id, dsize, nc

    if (lvl <= mg%lowest_lvl) error stop "cannot restrict lvl <= lowest_lvl"

    nc = mg%box_size_lvl(lvl)

    if (lvl >= mg%first_normal_lvl) then
       dsize = (nc/2)**1

       mg%buf(:)%i_send = 0
       mg%buf(:)%i_ix   = 0

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call restrict_set_buffer(mg, id, iv)
       end do

       mg%buf(:)%i_recv = mg%comm_restrict%n_recv(:, lvl) * dsize
       call sort_and_transfer_buffers(mg, dsize)
       mg%buf(:)%i_recv = 0
    end if

    do i = 1, size(mg%lvls(lvl-1)%my_parents)
       id = mg%lvls(lvl-1)%my_parents(i)
       call restrict_onto(mg, id, nc, iv)
    end do
  end subroutine mg_restrict_lvl

  subroutine restrict_set_buffer(mg, id, iv)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer, intent(in)          :: iv
    integer                      :: i, n, hnc, p_id, p_rank
    real(dp) :: tmp(mg%box_size/2)

    hnc    = mg%box_size/2
    p_id   = mg%boxes(id)%parent
    p_rank = mg%boxes(p_id)%rank

    if (p_rank /= mg%my_rank) then
       do i = 1, hnc
          tmp(i) = 0.5_dp * &
               sum(mg%boxes(id)%cc(2*i-1:2*i, iv))
       end do

       ! Buffer
       n = size(tmp)
       i = mg%buf(p_rank)%i_send
       mg%buf(p_rank)%send(i+1:i+n) = pack(tmp, .true.)
       mg%buf(p_rank)%i_send = mg%buf(p_rank)%i_send + n

       ! To later sort the send buffer according to parent order
       i = mg%buf(p_rank)%i_ix
       n = mg_ix_to_ichild(mg%boxes(id)%ix)
       mg%buf(p_rank)%ix(i+1) = mg_num_children * p_id + n
       mg%buf(p_rank)%i_ix = mg%buf(p_rank)%i_ix + 1
    end if
  end subroutine restrict_set_buffer

  subroutine restrict_onto(mg, id, nc, iv)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    integer                   :: i, hnc, dsize, i_c, c_id
    integer                   :: c_rank, dix(1)

    hnc   = nc/2
    dsize = hnc**1

    do i_c = 1, mg_num_children
       c_id   = mg%boxes(id)%children(i_c)
       if (c_id == mg_no_box) cycle ! For coarsened grid
       c_rank = mg%boxes(c_id)%rank
       dix    = mg_get_child_offset(mg, c_id)

       if (c_rank == mg%my_rank) then
          do i=1, hnc
             mg%boxes(id)%cc(dix(1)+i, iv) = 0.5_dp * &
                  sum(mg%boxes(c_id)%cc(2*i-1:2*i, iv))
          end do; 
       else
          i = mg%buf(c_rank)%i_recv
          mg%boxes(id)%cc(dix(1)+1:dix(1)+hnc, iv) = &
               reshape(mg%buf(c_rank)%recv(i+1:i+dsize), [hnc])
          mg%buf(c_rank)%i_recv = mg%buf(c_rank)%i_recv + dsize
       end if
    end do

  end subroutine restrict_onto

  !! File ../src/m_mrgrnk.f90

   Subroutine I_mrgrnk (XDONT, IRNGT)
      ! __________________________________________________________
      !   MRGRNK = Merge-sort ranking of an array
      !   For performance reasons, the first 2 passes are taken
      !   out of the standard loop, and use dedicated coding.
      ! __________________________________________________________
      ! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      ! __________________________________________________________
      Integer :: XVALA, XVALB
      !
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
      !
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
      !
      !  Fill-in the index array, creating ordered couples
      !
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
      !
      !  We will now have ordered subsets A - B - A - B - ...
      !  and merge A and B couples into     C   -   C   - ...
      !
      LMTNA = 2
      LMTNC = 4
      !
      !  First iteration. The length of the ordered subsets goes from 2 to 4
      !
      Do
         If (NVAL <= 2) Exit
         !
         !   Loop on merges of A and B into C
         !
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
               !
               !   1 2 3
               !
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
               !
               !   1 3 2
               !
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
                  !
                  !   3 1 2
                  !
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
            !
            !   1 2 3 4
            !
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
            !
            !   1 3 x x
            !
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
                  !   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
                  !   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
               !
               !   3 x x x
               !
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
                     !   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
                     !   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
                  !   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
         !
         !  The Cs become As and Bs
         !
         LMTNA = 4
         Exit
      End Do
      !
      !  Iteration loop. Each time, the length of the ordered subsets
      !  is doubled.
      !
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
         !
         !   Loop on merges of A and B into C
         !
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
            !
            !   Shortcut for the case when the max of A is smaller
            !   than the min of B. This line may be activated when the
            !   initial set is already close to sorted.
            !
            !          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
            !
            !  One steps in the C subset, that we build in the final rank array
            !
            !  Make a copy of the rank array for the merge iteration
            !
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
            !
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
            !
            Do
               IWRK = IWRK + 1
               !
               !  We still have unprocessed values in both A and B
               !
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
                     !  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
               !
            End Do
         End Do
         !
         !  The Cs become As and Bs
         !
         LMTNA = 2 * LMTNA
      End Do
      !
      Return
      !
   End Subroutine I_mrgrnk

end module m_octree_mg_1d
