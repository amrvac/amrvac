#include "cpp_macros.h"
module m_build_tree
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: mg_build_rectangle
  public :: mg_add_children
  public :: mg_set_leaves_parents
  public :: mg_set_next_level_ids
  public :: mg_set_refinement_boundaries
  public :: mg_set_neighbors_lvl

contains

  subroutine mg_build_rectangle(mg, domain_size, box_size, dx, r_min, &
       periodic, n_finer)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: domain_size(NDIM)
    integer, intent(in)       :: box_size
    real(dp), intent(in)      :: dx(NDIM)
    real(dp), intent(in)      :: r_min(NDIM)
    logical, intent(in)       :: periodic(NDIM)
    integer, intent(in)       :: n_finer
    integer                   :: IJK, lvl, n, id, nx(NDIM)
    integer                   :: boxes_per_dim(NDIM, mg_lvl_lo:1)
    integer                   :: periodic_offset(NDIM)

    if (modulo(box_size, 2) /= 0) &
         error stop "box_size should be even"
    if (any(modulo(domain_size, box_size) /= 0)) &
         error stop "box_size does not divide domain_size"

    if (all(periodic)) then
       mg%subtract_mean = .true.
    end if

    nx                       = domain_size
    mg%box_size              = box_size
    mg%box_size_lvl(1)       = box_size
    mg%domain_size_lvl(:, 1) = domain_size
    mg%first_normal_lvl      = 1
    mg%dr(:, 1)              = dx
    mg%r_min(:)              = r_min
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
#if NDIM == 2
    periodic_offset = [nx(1)-1, (nx(2)-1)*nx(1)]
#elif NDIM == 3
    periodic_offset = [nx(1)-1, (nx(2)-1)*nx(1), &
         (nx(3)-1) * nx(2) * nx(1)]
#endif

    do KJI_DO_VEC(nx)
       mg%n_boxes = mg%n_boxes + 1
       n          = mg%n_boxes

       mg%boxes(n)%rank        = 0
       mg%boxes(n)%id          = n
       mg%boxes(n)%lvl         = mg%lowest_lvl
       mg%boxes(n)%ix(:)       = [IJK]
       mg%boxes(n)%r_min(:)    = r_min + (mg%boxes(n)%ix(:) - 1) * &
            mg%box_size_lvl(mg%lowest_lvl) * mg%dr(:, mg%lowest_lvl)
       mg%boxes(n)%dr(:)       = mg%dr(:, mg%lowest_lvl)
       mg%boxes(n)%parent      = mg_no_box
       mg%boxes(n)%children(:) = mg_no_box

       ! Set default neighbors
#if NDIM == 2
       mg%boxes(n)%neighbors(:) = [n-1, n+1, n-nx(1), n+nx(1)]
#elif NDIM == 3
       mg%boxes(n)%neighbors(:) = [n-1, n+1, n-nx(1), n+nx(1), &
            n-nx(1)*nx(2), n+nx(1)*nx(2)]
#endif

       ! Handle boundaries
       where ([IJK] == 1 .and. .not. periodic)
          mg%boxes(n)%neighbors(1:mg_num_neighbors:2) = &
               mg_physical_boundary
       end where
       where ([IJK] == 1 .and. periodic)
          mg%boxes(n)%neighbors(1:mg_num_neighbors:2) = &
               n + periodic_offset
       end where

       where ([IJK] == nx .and. .not. periodic)
          mg%boxes(n)%neighbors(2:mg_num_neighbors:2) = &
               mg_physical_boundary
       end where
       where ([IJK] == nx .and. periodic)
          mg%boxes(n)%neighbors(2:mg_num_neighbors:2) = &
               n - periodic_offset
       end where
    end do; CLOSE_DO

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
    integer                   :: lvl, i, nb, child_nb(2**(NDIM-1))
    integer                   :: c_ids(mg_num_children)
    integer                   :: c_id, c_ix_base(NDIM)

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

end module m_build_tree
