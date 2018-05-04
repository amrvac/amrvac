module mod_usr
  use mod_mhd
  use mod_multigrid_coupling

  implicit none

  logical :: use_mg = .false.
  integer :: my_rhs
  integer :: my_phi

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

    my_phi = var_set_extravar("phi", "phi")
    my_rhs = var_set_extravar("rhs", "rhs")

    call params_read(par_files)

    if (use_mg) then
       ! phys_modify_wLR => modify_Bface
       usr_process_global => compute_phi
       usr_before_main_loop => mg_setup_multigrid
       usr_after_refine => mg_update_refinement
    end if
  end subroutine usr_init

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /my_list/ use_mg
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_list, end=111)
111    close(unitpar)
    end do
  end subroutine params_read

  subroutine subtract_grad_phi()
    use mod_global_parameters
    use mod_forest
    use mod_ghostcells_update
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl, i, j, i0, j0
    type(tree_node), pointer :: pnode
    real(dp)                 :: inv_2dr(2), gradient(2)
    real(dp)                 :: B_new(2), B_old(2), en_diff
    real(dp)                 :: divb(ixG^T)

    if (.not. mg%is_allocated) &
         error stop "multigrid tree not allocated yet"

    i0 = nghostcells
    j0 = nghostcells

    do iigrid = 1, igridstail
       igrid   =  igrids(iigrid);
       pnode   => igrid_to_node(igrid, mype)%node
       id      =  pnode%id
       lvl     =  mg%boxes(id)%lvl
       nc      =  mg%box_size_lvl(lvl)
       inv_2dr =  1 / (2 * mg%dr(:, lvl))

       do j = 1, nc
          do i = 1, nc
             gradient(1) = (mg%boxes(id)%cc(i+1, j, mg_iphi) - &
                  mg%boxes(id)%cc(i-1, j, mg_iphi)) * inv_2dr(1)
             gradient(2) = (mg%boxes(id)%cc(i, j+1, mg_iphi) - &
                  mg%boxes(id)%cc(i, j-1, mg_iphi)) * inv_2dr(2)
             B_old = pw(igrid)%w(i0+i, j0+j, mag(:))
             B_new = B_old - gradient(:)
             en_diff = 0.5_dp * (sum(B_new**2) - sum(B_old**2))
             ! Keep pressure the same
             pw(igrid)%w(i0+i, j0+j, e_) = pw(igrid)%w(i0+i, j0+j, e_) + en_diff
             pw(igrid)%w(i0+i, j0+j, mag(:)) = B_new
          end do
       end do
    end do

    call getbc(global_time, 0.0_dp, 0, nwflux, .false.)

  end subroutine subtract_grad_phi

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

  subroutine divb_to_rhs()
    use mod_global_parameters
    use mod_forest
    integer                  :: iigrid, igrid, id
    integer                  :: nc, lvl
    type(tree_node), pointer :: pnode
    real(dp)                 :: divb(ixG^T)

    if (.not. mg%is_allocated) &
         error stop "multigrid tree not allocated yet"

    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! typediv = 'limited'
       ! call get_divb(pw(igrid)%w(ixG^T, 1:nw), ixG^LL, ixM^LL, divb)
       call mydivvector(pw(igrid)%w(ixG^T, mag(:)), ixG^LL, ixM^LL, divb, .true.)
       ! typediv = 'central'
       mg%boxes(id)%cc(1:nc, 1:nc, mg_irhs) = divb(ixM^T)
    end do
  end subroutine divb_to_rhs

  subroutine compute_phi(iit,qt)
    use mod_global_parameters
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer                      :: id

    call divb_to_rhs()
    call mg_fas_vcycle(mg)
    call subtract_grad_phi()

    call mg_copy_from_tree(mg_irhs, my_rhs)
    call mg_copy_from_tree_gc(mg_iphi, my_phi)

  end subroutine compute_phi

  subroutine modify_Bface(wLC, wRC, ixI^L, ixO^L, idir)
    use mod_forest
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    type(tree_node), pointer        :: pnode
    integer                         :: id, lvl, nc, i0, j0, i, j
    real(dp)                        :: inv_dr(2), grad

    pnode  => igrid_to_node(saveigrid, mype)%node
    id     =  pnode%id
    lvl    =  mg%boxes(id)%lvl
    nc     =  mg%box_size_lvl(lvl)
    inv_dr =  1 / mg%dr(:, lvl)
    i0     =  nghostcells
    j0     =  nghostcells

    select case (idir)
    case (1)
       do j = 1, nc
          do i = 1, nc+1
             grad = (mg%boxes(id)%cc(i, j, mg_iphi) - &
                  mg%boxes(id)%cc(i-1, j, mg_iphi)) * inv_dr(1)
             wLC(i0+i-1, j0+j, mag(1)) = wLC(i0+i-1, j0+j, mag(1)) - grad
             wRC(i0+i-1, j0+j, mag(1)) = wRC(i0+i-1, j0+j, mag(1)) - grad
          end do
       end do
    case (2)
       do j = 1, nc+1
          do i = 1, nc
             grad = (mg%boxes(id)%cc(i, j, mg_iphi) - &
                  mg%boxes(id)%cc(i, j-1, mg_iphi)) * inv_dr(2)
             wLC(i0+i, j0+j-1, mag(2)) = wLC(i0+i, j0+j-1, mag(2)) - grad
             wRC(i0+i, j0+j-1, mag(2)) = wRC(i0+i, j0+j-1, mag(2)) - grad
          end do
       end do
    end select

  end subroutine modify_Bface

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

  subroutine mydivvector(qvec,ixI^L,ixO^L,divq, fourthorder)
    ! Calculate divergence of a vector qvec within ixL
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qvec(ixI^S,1:ndir)
    double precision, intent(inout) :: divq(ixI^S)
    logical, intent(in), optional   :: fourthorder
    logical                         :: use_4th_order
    double precision                :: qC(ixI^S), invdx(1:ndim)
    integer                         :: jxO^L, hxO^L, ixC^L, jxC^L
    integer                         :: idims, ix^L, gxO^L, kxO^L

    use_4th_order = .false.
    if (present(fourthorder)) use_4th_order = fourthorder

    if (.not. use_4th_order) then
       ix^L=ixO^L^LADD1;
    else
       ix^L=ixO^L^LADD2;
    end if

    if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
         call mpistop("Error in divvector: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixO^S)=zero

    do idims=1,ndim
       if (.not. use_4th_order) then
          jxO^L=ixO^L+kr(idims,^D);
          hxO^L=ixO^L-kr(idims,^D);
          divq(ixO^S)=divq(ixO^S)+half*(qvec(jxO^S,idims)-qvec(hxO^S,idims))*invdx(idims)
       else
          kxO^L=ixO^L+2*kr(idims,^D);
          jxO^L=ixO^L+kr(idims,^D);
          hxO^L=ixO^L-kr(idims,^D);
          gxO^L=ixO^L-2*kr(idims,^D);
          divq(ixO^S)=divq(ixO^S)+&
               (-qvec(kxO^S,idims) + 8.0d0 * qvec(jxO^S,idims) - 8.0d0 * &
               qvec(hxO^S,idims) + qvec(gxO^S,idims))/(12.0d0 * dxlevel(idims))
       end if
    end do
  end subroutine mydivvector

end module mod_usr
