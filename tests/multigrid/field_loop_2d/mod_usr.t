module mod_usr
  use mod_mhd
  use mod_multigrid_coupling

  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_physics

    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output
    usr_refine_grid => my_refine

    call set_coordinate_system('Cartesian')
    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: v0
    double precision, parameter     :: A0 = 1.0d-3

    select case (iprob)
    case (1)
       v0 = 1.0
    case (2)
       v0 = 0.0
    case default
       call mpistop("Invalid iprob")
    end select

    w(ixO^S,rho_) = 1.0d0       ! Density
    w(ixO^S,mom(1))= v0 * 2     ! Vx
    w(ixO^S,mom(2))= v0         ! Vy
    w(ixO^S,e_)   = 1.0d0       ! Pressure

    where (x(ixO^S,1)**2 + x(ixO^S,2)**2 < 0.3d0**2)
       w(ixO^S,mag(1))= A0 * x(ixO^S,2)/sqrt(x(ixO^S,1)**2 + x(ixO^S,2)**2)
       w(ixO^S,mag(2))= -A0 * x(ixO^S,1)/sqrt(x(ixO^S,1)**2 + x(ixO^S,2)**2)
    elsewhere
       w(ixO^S,mag(1))= 0.0d0
       w(ixO^S,mag(2))= 0.0d0
    end where
    call mhd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: divb(ixI^S)

    ! output divB1
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+1)=divb(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='divb'
  end subroutine specialvarnames_output

  ! subroutine modify_Bface(wLC, wRC, ixI^L, ixO^L, idir)
  !   use mod_forest
  !   use mod_global_parameters
  !   integer, intent(in)             :: ixI^L, ixO^L, idir
  !   double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
  !   type(tree_node), pointer        :: pnode
  !   integer                         :: id, lvl, nc, i0, j0, i, j
  !   real(dp)                        :: inv_dr(2), grad

  !   pnode  => igrid_to_node(block%igrid, mype)%node
  !   id     =  pnode%id
  !   lvl    =  mg%boxes(id)%lvl
  !   nc     =  mg%box_size_lvl(lvl)
  !   inv_dr =  1 / mg%dr(:, lvl)
  !   i0     =  nghostcells
  !   j0     =  nghostcells

  !   select case (idir)
  !   case (1)
  !      do j = 1, nc
  !         do i = 1, nc+1
  !            grad = (mg%boxes(id)%cc(i, j, mg_iphi) - &
  !                 mg%boxes(id)%cc(i-1, j, mg_iphi)) * inv_dr(1)
  !            wLC(i0+i-1, j0+j, mag(1)) = wLC(i0+i-1, j0+j, mag(1)) - grad
  !            wRC(i0+i-1, j0+j, mag(1)) = wRC(i0+i-1, j0+j, mag(1)) - grad
  !         end do
  !      end do
  !   case (2)
  !      do j = 1, nc+1
  !         do i = 1, nc
  !            grad = (mg%boxes(id)%cc(i, j, mg_iphi) - &
  !                 mg%boxes(id)%cc(i, j-1, mg_iphi)) * inv_dr(2)
  !            wLC(i0+i, j0+j-1, mag(2)) = wLC(i0+i, j0+j-1, mag(2)) - grad
  !            wRC(i0+i, j0+j-1, mag(2)) = wRC(i0+i, j0+j-1, mag(2)) - grad
  !         end do
  !      end do
  !   end select

  ! end subroutine modify_Bface

  subroutine my_refine(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    refine = 0
    coarsen = 0
    if (any(x(ixO^S, 1) < 0.0d0)) then
       refine = 1
       coarsen = -1
    else
       refine = -1
       coarsen = -1
    end if
  end subroutine my_refine

end module mod_usr
