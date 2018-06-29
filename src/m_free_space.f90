!> Module to use free space boundary conditions for 3D Poisson problems
module m_free_space
  use m_data_structures

  implicit none
  private

#if NDIM == 3
  type mg_free_bc_t
     logical               :: initialized = .false.
     integer               :: fft_lvl = -1
     real(dp), allocatable :: rhs(:, :, :)
     real(dp), pointer     :: karray(:)
     real(dp)              :: inv_dr(2, mg_num_neighbors)
     real(dp)              :: r_min(2, mg_num_neighbors)
     real(dp), allocatable :: bc_x0(:, :)
     real(dp), allocatable :: bc_x1(:, :)
     real(dp), allocatable :: bc_y0(:, :)
     real(dp), allocatable :: bc_y1(:, :)
     real(dp), allocatable :: bc_z0(:, :)
     real(dp), allocatable :: bc_z1(:, :)
  end type mg_free_bc_t

  type(mg_free_bc_t) :: free_bc

  ! Public methods
  public :: mg_poisson_free_3d

contains

  !> Solve a free-space Poisson problem in 3D, making use of a FFT solver (on a
  !> coarser grid) to get the boundary conditions. Each call performs an
  !> additional FMG or V-cycle, depending on the argument fmgcycle.
  subroutine mg_poisson_free_3d(mg, new_rhs, max_fft_frac, fmgcycle, max_res)
    use mpi
    use poisson_solver
    use m_multigrid
    use m_restrict
    use m_prolong
    use m_ghost_cells
    type(mg_t), intent(inout)       :: mg
    !> Whether a new right-hand side is present
    logical, intent(in)             :: new_rhs
    !> How much smaller the fft solve has to be than the full multigrid (0.0-1.0)
    real(dp), intent(in)            :: max_fft_frac
    !> If false, perform a v-cycle instead of FMG
    logical, intent(in)             :: fmgcycle
    !> If present, compute and output the maximum residual
    real(dp), intent(out), optional :: max_res
    integer                         :: fft_lvl, lvl, n, id, nx(3), nc
    integer                         :: ix(3), ierr
    integer(i8)                     :: n_unknowns_lvl, n_unknowns_total
    real(dp)                        :: dr(3)
    real(dp), allocatable           :: tmp(:, :, :)
    real(dp)                        :: dummy(1)
    real(dp), parameter             :: offset    = 0.0_dp
    real(dp)                        :: ehartree, eexcu, vexcu
    integer                         :: i3sd, ncomp
    character(len=*), parameter     :: geocode   = 'F'
    character(len=*), parameter     :: datacode  = 'G'
    integer, parameter              :: itype_scf = 8
    integer, parameter              :: ixc       = 0
    logical                         :: new_fft_grid

    ! Correction factor for the source term 1 / (4 * pi)
    real(dp), parameter :: rhs_fac = -1 / (4 * acos(-1.0_dp))

    if (.not. free_bc%initialized .and. .not. new_rhs) then
       error stop "mg_poisson_free_3d: first call requires new_rhs = .true."
    end if

    if (mg%geometry_type /= mg_cartesian) &
         error stop "mg_poisson_free_3d: Cartesian 3D geometry required"

    if (mg%operator_type /= mg_laplacian) &
         error stop "mg_poisson_free_3d: laplacian operator required"

    ! Determine total number of unknowns (on leaves)
    n_unknowns_total = mg_number_of_unknowns(mg)

    ! Determine highest fully refined grid level
    do lvl = mg_highest_uniform_lvl(mg), mg%lowest_lvl+1, -1
       ! Determine how many boxes the level contains
       n_unknowns_lvl = size(mg%lvls(lvl)%ids) * int(mg%box_size**3, i8)

       ! If the level is 'small enough', exit
       if (n_unknowns_lvl <= max_fft_frac * n_unknowns_total) exit
    end do

    fft_lvl = lvl
    new_fft_grid = (.not. free_bc%initialized .or. free_bc%fft_lvl /= fft_lvl)

    ! Add a layer of ghost cells around the domain
    nx(:) = mg%domain_size_lvl(:, fft_lvl) + 2
    dr(:) = mg%dr(:, fft_lvl)

    ! Ensure boundary conditions are set for multigrid solver
    do n = 1, mg_num_neighbors
       mg%bc(n, mg_iphi)%boundary_cond => ghost_cells_free_bc
    end do

    if (free_bc%initialized .and. new_fft_grid) then
       deallocate(free_bc%karray)
       deallocate(free_bc%rhs)
       ! The boundary planes are automatically re-allocated
       free_bc%initialized = .false.
    end if

    if (new_fft_grid) then
       free_bc%fft_lvl = fft_lvl

       ! Restrict rhs to required level
       do lvl = mg%highest_lvl, fft_lvl+1, -1
          call mg_restrict_lvl(mg, mg_irhs, lvl)
       end do

       ! Create kernel of Green's function
       call createKernel(geocode, nx(1), nx(3), nx(3), dr(1), dr(2), dr(3),  &
            itype_scf, mg%my_rank, mg%n_cpu, free_bc%karray)

       allocate(free_bc%rhs(nx(1), nx(2), nx(3)))

       ! For interpolation of the boundary planes
       free_bc%inv_dr(:, mg_neighb_lowx) = 1 / dr(2:3)
       free_bc%r_min(:, mg_neighb_lowx)  = mg%r_min(2:3) - 0.5_dp * dr(2:3)
       free_bc%inv_dr(:, mg_neighb_lowy) = 1 / dr([1,3])
       free_bc%r_min(:, mg_neighb_lowy)  = mg%r_min([1,3]) - 0.5_dp * dr([1,3])
       free_bc%inv_dr(:, mg_neighb_lowz) = 1 / dr(1:2)
       free_bc%r_min(:, mg_neighb_lowz)  = mg%r_min(1:2) - 0.5_dp * dr(1:2)

       do n = mg_neighb_lowx+1, mg_num_neighbors, 2
          free_bc%inv_dr(:, n) = free_bc%inv_dr(:, n-1)
          free_bc%r_min(:, n)  = free_bc%r_min(:, n-1)
       end do

    end if

    if (new_rhs) then
       allocate(tmp(nx(1), nx(2), nx(3)))
       tmp(:, :, :) = 0.0_dp

       ! Store right-hand side
       nc = mg%box_size_lvl(fft_lvl)
       do n = 1, size(mg%lvls(fft_lvl)%my_ids)
          id = mg%lvls(fft_lvl)%my_ids(n)
          ix = (mg%boxes(id)%ix - 1) * nc + 1
          tmp(ix(1)+1:ix(1)+nc, ix(2)+1:ix(2)+nc, ix(3)+1:ix(3)+nc) = &
               rhs_fac * mg%boxes(id)%cc(1:nc, 1:nc, 1:nc, mg_irhs)
       end do

       call mpi_allreduce(tmp, free_bc%rhs, product(shape(tmp)), MPI_DOUBLE, &
            MPI_SUM, mg%comm, ierr)

       ! Use default load balancing for parallel fft
       i3sd  = 1
       ncomp = nx(3)

       ! Solve free-space Poisson's equation
       call PSolver(geocode, datacode, mg%my_rank, mg%n_cpu, &
            nx(1), nx(2), nx(3), ixc, dr(1), dr(2), dr(3), &
            free_bc%rhs(1, 1, i3sd), free_bc%karray, dummy, &
            ehartree, eexcu, vexcu, offset, .false., 1)

       ! Extract boundary planes by interpolation
       associate (rhs => free_bc%rhs)
         free_bc%bc_x0 = 0.5_dp * (rhs(1, :, :) + rhs(2, :, :))
         free_bc%bc_x1 = 0.5_dp * (rhs(nx(1)-1, :, :) + rhs(nx(1), :, :))
         free_bc%bc_y0 = 0.5_dp * (rhs(:, 1, :) + rhs(:, 2, :))
         free_bc%bc_y1 = 0.5_dp * (rhs(:, nx(2)-1, :) + rhs(:, nx(2), :))
         free_bc%bc_z0 = 0.5_dp * (rhs(:, :, 1) + rhs(:, :, 2))
         free_bc%bc_z1 = 0.5_dp * (rhs(:, :, nx(3)-1) + rhs(:, :, nx(3)))
       end associate

       ! Store boundary conditions (so interpolation has to be performed only once)
       call mg_phi_bc_store(mg)

       ! Use solution as an initial guess
       nc = mg%box_size_lvl(fft_lvl)
       do n = 1, size(mg%lvls(fft_lvl)%my_ids)
          id = mg%lvls(fft_lvl)%my_ids(n)
          ix = (mg%boxes(id)%ix - 1) * nc + 1
          mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, mg_iphi) = &
               free_bc%rhs(ix(1):ix(1)+nc+1, ix(2):ix(2)+nc+1, &
               ix(3):ix(3)+nc+1)
       end do

       ! Restrict FFT solution
       do lvl = fft_lvl, mg%lowest_lvl+1, -1
          call mg_restrict_lvl(mg, mg_iphi, lvl)
       end do

       ! Prolong FFT solution
       do lvl = fft_lvl, mg%highest_lvl-1
          call mg_prolong(mg, lvl, mg_iphi, mg_iphi, mg%box_prolong, .false.)
          ! We can already use the boundary conditions from the FFT solution
          call mg_fill_ghost_cells_lvl(mg, lvl+1, mg_iphi)
       end do

       free_bc%initialized = .true.
    end if

    ! Avoid multigrid solver if we already have the full solution
    if (free_bc%fft_lvl < mg%highest_lvl) then
       ! Solve Poisson equation with free space boundary conditions
       if (fmgcycle) then
          call mg_fas_fmg(mg, .true., max_res)
       else
          call mg_fas_vcycle(mg, max_res=max_res)
       end if
    end if

  end subroutine mg_poisson_free_3d

  !> To fill ghost cells
  subroutine ghost_cells_free_bc(box, nc, iv, nb, bc_type, bc)
    type(mg_box_t), intent(in)    :: box
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv      !< Index of variable
    integer, intent(in)           :: nb      !< Direction
    integer, intent(out)          :: bc_type !< Type of b.c.
    double precision, intent(out) :: bc(nc, nc)
    double precision              :: rr(nc, nc, 3)
    integer                       :: ixs(2), nb_dim

    bc_type      = mg_bc_dirichlet
    nb_dim       = mg_neighb_dim(nb)
    ixs          = [1, 2]
    ixs(nb_dim:) = ixs(nb_dim:) + 1

    call mg_get_face_coords(box, nb, nc, rr)

    bc = interp_bc(free_bc, nb, rr(:, :, ixs(1)), rr(:, :, ixs(2)))
  end subroutine ghost_cells_free_bc

  elemental function interp_bc(bc, nb, x1, x2) result(val)
    type(mg_free_bc_t), intent(in) :: bc
    integer, intent(in)            :: nb
    real(dp), intent(in)           :: x1, x2
    real(dp)                       :: val
    integer                        :: ix(2)
    real(dp)                       :: frac(2), low_frac(2)
    real(dp)                       :: w(2, 2)

    frac     = ([x1, x2] - bc%r_min(:, nb)) * bc%inv_dr(:, nb)
    ix       = ceiling(frac)
    low_frac = ix - frac

    ! Bilinear interpolation
    w(1, 1) = low_frac(1) * low_frac(2)
    w(2, 1) = (1 - low_frac(1)) * low_frac(2)
    w(1, 2) = low_frac(1) * (1 - low_frac(2))
    w(2, 2) = (1 - low_frac(1)) * (1 - low_frac(2))

    select case (nb)
    case (mg_neighb_lowx)
       val = sum(w * bc%bc_x0(ix(1):ix(1)+1, ix(2):ix(2)+1))
    case (mg_neighb_highx)
       val = sum(w * bc%bc_x1(ix(1):ix(1)+1, ix(2):ix(2)+1))
    case (mg_neighb_lowy)
       val = sum(w * bc%bc_y0(ix(1):ix(1)+1, ix(2):ix(2)+1))
    case (mg_neighb_highy)
       val = sum(w * bc%bc_y1(ix(1):ix(1)+1, ix(2):ix(2)+1))
    case (mg_neighb_lowz)
       val = sum(w * bc%bc_z0(ix(1):ix(1)+1, ix(2):ix(2)+1))
    case (mg_neighb_highz)
       val = sum(w * bc%bc_z1(ix(1):ix(1)+1, ix(2):ix(2)+1))
    end select
  end function interp_bc
#endif
end module m_free_space
