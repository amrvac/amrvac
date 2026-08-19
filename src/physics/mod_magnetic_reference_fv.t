!> Finite-volume potential reference field for magnetic diagnostics.
!>
!> The reference field Bp = grad(phi) is obtained from a Laplace solve with
!> the normal component of the total magnetic field prescribed on all six
!> physical faces.  The multigrid solution is scratch storage only: callers
!> receive an independently allocated, distributed cell-centred Bp field.
module mod_magnetic_reference_fv
  implicit none
  private

  type, public :: magnetic_reference_config
    double precision :: residual_tolerance=1.d-8
    double precision :: max_flux_imbalance=1.d-6
    integer :: max_cycles=50
    !> Print the multigrid timer table after a successful solve.  This is
    !> deliberately opt-in so ordinary scientific output is unchanged.
    logical :: write_timing=.false.
  end type magnetic_reference_config

  type, public :: magnetic_reference_result
    double precision :: magnetic_energy=0.d0
    double precision :: residual=0.d0
    double precision :: flux_imbalance=0.d0
    double precision :: boundary_normal_error=0.d0
    integer :: cycles=0
  end type magnetic_reference_result

  type, public :: magnetic_reference_block
    double precision, allocatable :: b(:^D&,:)
  end type magnetic_reference_block

  type, public :: magnetic_reference_field
    type(magnetic_reference_block), allocatable :: blocks(:)
  end type magnetic_reference_field

  public :: solve_magnetic_reference_fv
  public :: free_magnetic_reference_field
  {^IFTHREED
  public :: magnetic_reference_bc
  }

contains

  subroutine solve_magnetic_reference_fv(config,bp,result)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters
    use mod_geometry, only: Cartesian,coordinate
    use mod_multigrid_coupling
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    type(magnetic_reference_config), intent(in) :: config
    type(magnetic_reference_field), intent(inout) :: bp
    type(magnetic_reference_result), intent(out) :: result

    double precision :: net_flux,unsigned_flux,local_energy,global_energy
    double precision :: face_flux(6)
    double precision :: res
    integer :: i

    result=magnetic_reference_result()
    call free_magnetic_reference_field(bp)

    if(ndim/=3) call mpistop('finite-volume magnetic reference requires three dimensions')
    if(coordinate/=Cartesian) &
       call mpistop('finite-volume magnetic reference requires Cartesian coordinates')
    if(any(stretched_dim)) &
       call mpistop('finite-volume magnetic reference requires an unstretched mesh')
    if(any(periodB)) &
       call mpistop('finite-volume magnetic reference requires six physical boundaries')
    if(.not.allocated(iw_mag) .or. size(iw_mag)<3) &
       call mpistop('finite-volume magnetic reference requires a three-component magnetic field')
    if(config%residual_tolerance<=0.d0) &
       call mpistop('magnetic-reference residual tolerance must be positive')
    if(config%max_flux_imbalance<0.d0) &
       call mpistop('magnetic-reference flux tolerance must not be negative')
    if(config%max_cycles<1) &
       call mpistop('magnetic-reference max_cycles must be positive')

    {^IFTHREED
    if(.not.mg%is_allocated) call mg_setup_multigrid()

    mg%operator_type=mg_laplacian
    call mg_set_methods(mg)
    ! A Laplacian with six Neumann boundaries has an arbitrary constant.
    ! Removing the mean fixes that null mode without changing grad(phi).
    mg%subtract_mean=.true.
    mg%phi_bc_data_stored=.false.
    mg%bc(:,mg_iphi)%bc_type=mg_bc_neumann
    do i=1,2*ndim
      mg%bc(i,mg_iphi)%boundary_cond=>magnetic_reference_bc
    end do

    call magnetic_reference_flux_metrics(net_flux,unsigned_flux,face_flux)
    result%flux_imbalance=dabs(net_flux)/max(unsigned_flux,tiny(1.d0))
    if(result%flux_imbalance>config%max_flux_imbalance) then
      if(mype==0) then
        write(*,'(a,es14.6)') 'magnetic-reference net boundary flux=',net_flux
        write(*,'(a,es14.6)') 'magnetic-reference unsigned boundary flux=',unsigned_flux
        write(*,'(a,es14.6)') 'magnetic-reference relative imbalance=',&
           result%flux_imbalance
        write(*,'(a,6(1x,es14.6))') 'magnetic-reference coordinate face fluxes=',&
           face_flux
      end if
      call mpistop('finite-volume magnetic reference has incompatible Neumann flux')
    end if

    call magnetic_reference_initialize_mg()
    res=huge(1.d0)
    do i=1,config%max_cycles
      if(i==1) then
        ! Build the initial full-multigrid approximation once.  Re-entering
        ! FMG with have_guess=.true. repeats coarse restriction and correction
        ! work that is unnecessary for subsequent residual reduction.
        call mg_fas_fmg(mg,.false.,max_res=res)
      else
        ! Continue from the previous solution with the cheaper standalone
        ! V-cycle.  This is local to the reference-potential solve; other
        ! AMRVAC multigrid users retain their existing schedules.
        call mg_fas_vcycle(mg,max_res=res)
      end if
      if(mype==0) write(*,'(a,i0,a,es14.6)') &
         'magnetic-reference MG cycle ',i,' residual=',res
      if(res<=config%residual_tolerance) exit
    end do
    result%cycles=min(i,config%max_cycles)
    result%residual=res
    if(.not.ieee_is_finite(res)) &
       call mpistop('finite-volume magnetic reference produced a non-finite residual')
    if(res>config%residual_tolerance) &
       call mpistop('finite-volume magnetic reference did not converge')

    call mg_fill_ghost_cells(mg,mg_iphi)
    call magnetic_reference_copy_gradient(bp,local_energy)
    call MPI_ALLREDUCE(local_energy,global_energy,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,icomm,ierrmpi)
    result%magnetic_energy=global_energy
    call magnetic_reference_boundary_error(result%boundary_normal_error)
    if(config%write_timing) call mg_timers_show(mg)
    }
  end subroutine solve_magnetic_reference_fv

  subroutine free_magnetic_reference_field(bp)
    type(magnetic_reference_field), intent(inout) :: bp
    integer :: igrid

    if(.not.allocated(bp%blocks)) return
    do igrid=1,size(bp%blocks)
      if(allocated(bp%blocks(igrid)%b)) deallocate(bp%blocks(igrid)%b)
    end do
    deallocate(bp%blocks)
  end subroutine free_magnetic_reference_field

{^IFTHREED
  subroutine magnetic_reference_initialize_mg()
    use mod_global_parameters
    use mod_forest, only: igrid_to_node
    use mod_multigrid_coupling

    integer :: iigrid,igrid,id,nc,lvl

    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      id=igrid_to_node(igrid,mype)%node%id
      lvl=mg%boxes(id)%lvl
      nc=mg%box_size_lvl(lvl)
      mg%boxes(id)%cc({1:nc},mg_iphi)=0.d0
      mg%boxes(id)%cc({1:nc},mg_irhs)=0.d0
    end do
  end subroutine magnetic_reference_initialize_mg

  subroutine magnetic_reference_copy_gradient(bp,local_energy)
    use mod_global_parameters
    use mod_forest, only: igrid_to_node
    use mod_geometry, only: gradient
    use mod_multigrid_coupling

    type(magnetic_reference_field), intent(inout) :: bp
    double precision, intent(out) :: local_energy

    double precision, allocatable :: tmp(:,:,:)
    integer :: iigrid,igrid,idir,id,nc,lvl,ixI^L

    allocate(bp%blocks(max_blocks))
    ixImin^D=ixMlo^D-1;
    ixImax^D=ixMhi^D+1;
    allocate(tmp(ixI^S))
    local_energy=0.d0
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      id=igrid_to_node(igrid,mype)%node%id
      lvl=mg%boxes(id)%lvl
      nc=mg%box_size_lvl(lvl)
      allocate(bp%blocks(igrid)%b(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
         ixMlo3:ixMhi3,3))
      do idir=1,3
        tmp=0.d0
        call gradient(mg%boxes(id)%cc({0:nc+1},mg_iphi),ixI^L,&
           ixM^LL,idir,tmp)
        bp%blocks(igrid)%b(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
           ixMlo3:ixMhi3,idir)=tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
           ixMlo3:ixMhi3)
      end do
      local_energy=local_energy+sum(0.5d0*&
         sum(bp%blocks(igrid)%b(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
         ixMlo3:ixMhi3,:)**2,dim=4)*block%dvolume(ixMlo1:ixMhi1,&
         ixMlo2:ixMhi2,ixMlo3:ixMhi3))
    end do
    deallocate(tmp)
  end subroutine magnetic_reference_copy_gradient

  subroutine magnetic_reference_flux_metrics(net_flux,unsigned_flux,face_flux)
    use mod_global_parameters
    use mod_forest, only: igrid_to_node
    use mod_multigrid_coupling

    double precision, intent(out) :: net_flux,unsigned_flux,face_flux(6)
    double precision :: local_values(2),global_values(2),area
    double precision :: local_faces(6)
    double precision, allocatable :: bc(:,:)
    integer :: iigrid,igrid,id,lvl,nc,nb,idir,bc_type

    local_values=0.d0
    local_faces=0.d0
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      block=>ps(igrid)
      id=igrid_to_node(igrid,mype)%node%id
      lvl=mg%boxes(id)%lvl
      nc=mg%box_size_lvl(lvl)
      allocate(bc(nc,nc))
      do nb=1,6
        if(.not.ps(igrid)%is_physical_boundary(nb)) cycle
        call magnetic_reference_bc(mg%boxes(id),nc,mg_iphi,nb,bc_type,bc)
        idir=(nb+1)/2
        select case(idir)
        case(1)
          area=mg%boxes(id)%dr(2)*mg%boxes(id)%dr(3)
        case(2)
          area=mg%boxes(id)%dr(1)*mg%boxes(id)%dr(3)
        case(3)
          area=mg%boxes(id)%dr(1)*mg%boxes(id)%dr(2)
        end select
        if(mod(nb,2)==0) then
          local_values(1)=local_values(1)+sum(bc)*area
        else
          local_values(1)=local_values(1)-sum(bc)*area
        end if
        local_values(2)=local_values(2)+sum(dabs(bc))*area
        local_faces(nb)=local_faces(nb)+sum(bc)*area
      end do
      deallocate(bc)
    end do
    call MPI_ALLREDUCE(local_values,global_values,2,MPI_DOUBLE_PRECISION,&
       MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(local_faces,face_flux,6,MPI_DOUBLE_PRECISION,&
       MPI_SUM,icomm,ierrmpi)
    net_flux=global_values(1)
    unsigned_flux=global_values(2)
  end subroutine magnetic_reference_flux_metrics

  subroutine magnetic_reference_boundary_error(relative_error)
    use mod_global_parameters
    use mod_forest, only: igrid_to_node
    use mod_multigrid_coupling

    double precision, intent(out) :: relative_error
    double precision :: local_values(2),global_values(2),area,dr
    double precision, allocatable :: bc(:,:),derivative(:,:)
    integer :: iigrid,igrid,id,lvl,nc,nb,idir,bc_type

    local_values=0.d0
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      block=>ps(igrid)
      id=igrid_to_node(igrid,mype)%node%id
      lvl=mg%boxes(id)%lvl
      nc=mg%box_size_lvl(lvl)
      allocate(bc(nc,nc),derivative(nc,nc))
      do nb=1,6
        if(.not.ps(igrid)%is_physical_boundary(nb)) cycle
        call magnetic_reference_bc(mg%boxes(id),nc,mg_iphi,nb,bc_type,bc)
        idir=(nb+1)/2
        dr=mg%boxes(id)%dr(idir)
        select case(nb)
        case(1)
          derivative=(mg%boxes(id)%cc(1,1:nc,1:nc,mg_iphi)-&
             mg%boxes(id)%cc(0,1:nc,1:nc,mg_iphi))/dr
          area=mg%boxes(id)%dr(2)*mg%boxes(id)%dr(3)
        case(2)
          derivative=(mg%boxes(id)%cc(nc+1,1:nc,1:nc,mg_iphi)-&
             mg%boxes(id)%cc(nc,1:nc,1:nc,mg_iphi))/dr
          area=mg%boxes(id)%dr(2)*mg%boxes(id)%dr(3)
        case(3)
          derivative=(mg%boxes(id)%cc(1:nc,1,1:nc,mg_iphi)-&
             mg%boxes(id)%cc(1:nc,0,1:nc,mg_iphi))/dr
          area=mg%boxes(id)%dr(1)*mg%boxes(id)%dr(3)
        case(4)
          derivative=(mg%boxes(id)%cc(1:nc,nc+1,1:nc,mg_iphi)-&
             mg%boxes(id)%cc(1:nc,nc,1:nc,mg_iphi))/dr
          area=mg%boxes(id)%dr(1)*mg%boxes(id)%dr(3)
        case(5)
          derivative=(mg%boxes(id)%cc(1:nc,1:nc,1,mg_iphi)-&
             mg%boxes(id)%cc(1:nc,1:nc,0,mg_iphi))/dr
          area=mg%boxes(id)%dr(1)*mg%boxes(id)%dr(2)
        case(6)
          derivative=(mg%boxes(id)%cc(1:nc,1:nc,nc+1,mg_iphi)-&
             mg%boxes(id)%cc(1:nc,1:nc,nc,mg_iphi))/dr
          area=mg%boxes(id)%dr(1)*mg%boxes(id)%dr(2)
        end select
        local_values(1)=local_values(1)+sum((derivative-bc)**2)*area
        local_values(2)=local_values(2)+sum(bc**2)*area
      end do
      deallocate(bc,derivative)
    end do
    call MPI_ALLREDUCE(local_values,global_values,2,MPI_DOUBLE_PRECISION,&
       MPI_SUM,icomm,ierrmpi)
    relative_error=dsqrt(global_values(1)/max(global_values(2),tiny(1.d0)))
  end subroutine magnetic_reference_boundary_error
}

{^IFTHREED
  !> Supply d(phi)/dx_i = B_i on a physical multigrid boundary.
  subroutine magnetic_reference_bc(box,nc,iv,nb,bc_type,bc)
    use mod_global_parameters
    use mod_multigrid_coupling
    type(mg_box_t), intent(in) :: box
    integer, intent(in) :: nc,iv,nb
    integer, intent(out) :: bc_type
    double precision, intent(out) :: bc(nc,nc)

    double precision :: rr(nc,nc,3)
    double precision :: rmina,rminb,rmaxa,rmaxb,xmina,xminb,xmaxa,xmaxb
    double precision :: wbn(ixG^T)
    double precision, allocatable :: xcoarse(:,:,:)
    integer :: iigrid,igrid,ix^D,idir,ixbca,ixbcb,ixbcn,dlvl,wnc

    bc_type=mg_bc_neumann
    bc=0.d0

    call mg_get_face_coords(box,nb,nc,rr)
    idir=(nb+1)/2
    select case(idir)
    case(1)
      if(mod(nb,2)==0) then
        ixbcn=ixMhi1
      else
        ixbcn=ixMlo1-1
      end if
      rmina=rr(1,1,2)-0.5d0*box%dr(2)
      rmaxa=rr(nc,1,2)+0.5d0*box%dr(2)
      rminb=rr(1,1,3)-0.5d0*box%dr(3)
      rmaxb=rr(1,nc,3)+0.5d0*box%dr(3)
      do iigrid=1,igridstail
        igrid=igrids(iigrid)
        block=>ps(igrid)
        if(.not.block%is_physical_boundary(nb)) cycle
        if(stagger_grid) then
          wbn(ixbcn,ixGlo2:ixGhi2,ixGlo3:ixGhi3)=&
             block%ws(ixbcn,ixGlo2:ixGhi2,ixGlo3:ixGhi3,idir)
          if(B0field) wbn(ixbcn,ixGlo2:ixGhi2,ixGlo3:ixGhi3)=&
             wbn(ixbcn,ixGlo2:ixGhi2,ixGlo3:ixGhi3)+&
             block%B0(ixbcn,ixGlo2:ixGhi2,ixGlo3:ixGhi3,idir,idir)
        else
          wbn(ixbcn,ixGlo2:ixGhi2,ixGlo3:ixGhi3)=half*(&
             block%w(ixbcn,ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_mag(idir))+&
             block%w(ixbcn+1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_mag(idir)))
          if(B0field) wbn(ixbcn,ixGlo2:ixGhi2,ixGlo3:ixGhi3)=&
             wbn(ixbcn,ixGlo2:ixGhi2,ixGlo3:ixGhi3)+half*(&
             block%B0(ixbcn,ixGlo2:ixGhi2,ixGlo3:ixGhi3,idir,0)+&
             block%B0(ixbcn+1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,idir,0))
        end if
        xmina=block%x(1,1,1,2)-0.5d0*rnode(rpdx2_,igrid)
        xmaxa=block%x(1,ixGhi2,1,2)+0.5d0*rnode(rpdx2_,igrid)
        xminb=block%x(1,1,1,3)-0.5d0*rnode(rpdx3_,igrid)
        xmaxb=block%x(1,1,ixGhi3,3)+0.5d0*rnode(rpdx3_,igrid)
        if(xmina<rr(1,1,2) .and. xmaxa>rr(nc,1,2) .and.&
           xminb<rr(1,1,3) .and. xmaxb>rr(1,nc,3)) then
          do ix2=1,nc
            do ix1=1,nc
              ixbca=ceiling((rr(ix1,ix2,2)-xmina)/rnode(rpdx2_,igrid))
              ixbcb=ceiling((rr(ix1,ix2,3)-xminb)/rnode(rpdx3_,igrid))
              bc(ix1,ix2)=wbn(ixbcn,ixbca,ixbcb)
            end do
          end do
        else if(block%x(1,ixMlo2,1,2)>rmina .and.&
                block%x(1,ixMhi2,1,2)<rmaxa .and.&
                block%x(1,1,ixMlo3,3)>rminb .and.&
                block%x(1,1,ixMhi3,3)<rmaxb) then
          dlvl=node(plevel_,igrid)-box%lvl
          wnc=nc/2**dlvl
          allocate(xcoarse(wnc,wnc,2))
          do ix2=1,wnc
            do ix1=1,wnc
              xcoarse(ix1,ix2,1)=sum(block%x(1,&
                 (ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,&
                 1,2))/dble(2**dlvl)
              xcoarse(ix1,ix2,2)=sum(block%x(1,1,&
                 (ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells,&
                 3))/dble(2**dlvl)
              ixbca=ceiling((xcoarse(ix1,ix2,1)-rmina)/box%dr(2))
              ixbcb=ceiling((xcoarse(ix1,ix2,2)-rminb)/box%dr(3))
              bc(ixbca,ixbcb)=sum(wbn(ixbcn,&
                 (ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,&
                 (ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells))/&
                 dble(2**(2*dlvl))
            end do
          end do
          deallocate(xcoarse)
        end if
      end do
    case(2)
      if(mod(nb,2)==0) then
        ixbcn=ixMhi2
      else
        ixbcn=ixMlo2-1
      end if
      rmina=rr(1,1,1)-0.5d0*box%dr(1)
      rmaxa=rr(nc,1,1)+0.5d0*box%dr(1)
      rminb=rr(1,1,3)-0.5d0*box%dr(3)
      rmaxb=rr(1,nc,3)+0.5d0*box%dr(3)
      do iigrid=1,igridstail
        igrid=igrids(iigrid)
        block=>ps(igrid)
        if(.not.block%is_physical_boundary(nb)) cycle
        if(stagger_grid) then
          wbn(ixGlo1:ixGhi1,ixbcn,ixGlo3:ixGhi3)=&
             block%ws(ixGlo1:ixGhi1,ixbcn,ixGlo3:ixGhi3,idir)
          if(B0field) wbn(ixGlo1:ixGhi1,ixbcn,ixGlo3:ixGhi3)=&
             wbn(ixGlo1:ixGhi1,ixbcn,ixGlo3:ixGhi3)+&
             block%B0(ixGlo1:ixGhi1,ixbcn,ixGlo3:ixGhi3,idir,idir)
        else
          wbn(ixGlo1:ixGhi1,ixbcn,ixGlo3:ixGhi3)=half*(&
             block%w(ixGlo1:ixGhi1,ixbcn,ixGlo3:ixGhi3,iw_mag(idir))+&
             block%w(ixGlo1:ixGhi1,ixbcn+1,ixGlo3:ixGhi3,iw_mag(idir)))
          if(B0field) wbn(ixGlo1:ixGhi1,ixbcn,ixGlo3:ixGhi3)=&
             wbn(ixGlo1:ixGhi1,ixbcn,ixGlo3:ixGhi3)+half*(&
             block%B0(ixGlo1:ixGhi1,ixbcn,ixGlo3:ixGhi3,idir,0)+&
             block%B0(ixGlo1:ixGhi1,ixbcn+1,ixGlo3:ixGhi3,idir,0))
        end if
        xmina=block%x(1,1,1,1)-0.5d0*rnode(rpdx1_,igrid)
        xmaxa=block%x(ixGhi1,1,1,1)+0.5d0*rnode(rpdx1_,igrid)
        xminb=block%x(1,1,1,3)-0.5d0*rnode(rpdx3_,igrid)
        xmaxb=block%x(1,1,ixGhi3,3)+0.5d0*rnode(rpdx3_,igrid)
        if(xmina<rr(1,1,1) .and. xmaxa>rr(nc,1,1) .and.&
           xminb<rr(1,1,3) .and. xmaxb>rr(1,nc,3)) then
          do ix2=1,nc
            do ix1=1,nc
              ixbca=ceiling((rr(ix1,ix2,1)-xmina)/rnode(rpdx1_,igrid))
              ixbcb=ceiling((rr(ix1,ix2,3)-xminb)/rnode(rpdx3_,igrid))
              bc(ix1,ix2)=wbn(ixbca,ixbcn,ixbcb)
            end do
          end do
        else if(block%x(ixMlo1,1,1,1)>rmina .and.&
                block%x(ixMhi1,1,1,1)<rmaxa .and.&
                block%x(1,1,ixMlo3,3)>rminb .and.&
                block%x(1,1,ixMhi3,3)<rmaxb) then
          dlvl=node(plevel_,igrid)-box%lvl
          wnc=nc/2**dlvl
          allocate(xcoarse(wnc,wnc,2))
          do ix2=1,wnc
            do ix1=1,wnc
              xcoarse(ix1,ix2,1)=sum(block%x(&
                 (ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,&
                 1,1,1))/dble(2**dlvl)
              xcoarse(ix1,ix2,2)=sum(block%x(1,1,&
                 (ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells,&
                 3))/dble(2**dlvl)
              ixbca=ceiling((xcoarse(ix1,ix2,1)-rmina)/box%dr(1))
              ixbcb=ceiling((xcoarse(ix1,ix2,2)-rminb)/box%dr(3))
              bc(ixbca,ixbcb)=sum(wbn(&
                 (ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,&
                 ixbcn,&
                 (ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells))/&
                 dble(2**(2*dlvl))
            end do
          end do
          deallocate(xcoarse)
        end if
      end do
    case(3)
      if(mod(nb,2)==0) then
        ixbcn=ixMhi3
      else
        ixbcn=ixMlo3-1
      end if
      rmina=rr(1,1,1)-0.5d0*box%dr(1)
      rmaxa=rr(nc,1,1)+0.5d0*box%dr(1)
      rminb=rr(1,1,2)-0.5d0*box%dr(2)
      rmaxb=rr(1,nc,2)+0.5d0*box%dr(2)
      do iigrid=1,igridstail
        igrid=igrids(iigrid)
        block=>ps(igrid)
        if(.not.block%is_physical_boundary(nb)) cycle
        if(stagger_grid) then
          wbn(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn)=&
             block%ws(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn,idir)
          if(B0field) wbn(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn)=&
             wbn(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn)+&
             block%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn,idir,idir)
        else
          wbn(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn)=half*(&
             block%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn,iw_mag(idir))+&
             block%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn+1,iw_mag(idir)))
          if(B0field) wbn(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn)=&
             wbn(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn)+half*(&
             block%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn,idir,0)+&
             block%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixbcn+1,idir,0))
        end if
        xmina=block%x(1,1,1,1)-0.5d0*rnode(rpdx1_,igrid)
        xmaxa=block%x(ixGhi1,1,1,1)+0.5d0*rnode(rpdx1_,igrid)
        xminb=block%x(1,1,1,2)-0.5d0*rnode(rpdx2_,igrid)
        xmaxb=block%x(1,ixGhi2,1,2)+0.5d0*rnode(rpdx2_,igrid)
        if(xmina<rr(1,1,1) .and. xmaxa>rr(nc,1,1) .and.&
           xminb<rr(1,1,2) .and. xmaxb>rr(1,nc,2)) then
          do ix2=1,nc
            do ix1=1,nc
              ixbca=ceiling((rr(ix1,ix2,1)-xmina)/rnode(rpdx1_,igrid))
              ixbcb=ceiling((rr(ix1,ix2,2)-xminb)/rnode(rpdx2_,igrid))
              bc(ix1,ix2)=wbn(ixbca,ixbcb,ixbcn)
            end do
          end do
        else if(block%x(ixMlo1,1,1,1)>rmina .and.&
                block%x(ixMhi1,1,1,1)<rmaxa .and.&
                block%x(1,ixMlo2,1,2)>rminb .and.&
                block%x(1,ixMhi2,1,2)<rmaxb) then
          dlvl=node(plevel_,igrid)-box%lvl
          wnc=nc/2**dlvl
          allocate(xcoarse(wnc,wnc,2))
          do ix2=1,wnc
            do ix1=1,wnc
              xcoarse(ix1,ix2,1)=sum(block%x(&
                 (ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,&
                 1,1,1))/dble(2**dlvl)
              xcoarse(ix1,ix2,2)=sum(block%x(1,&
                 (ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells,&
                 1,2))/dble(2**dlvl)
              ixbca=ceiling((xcoarse(ix1,ix2,1)-rmina)/box%dr(1))
              ixbcb=ceiling((xcoarse(ix1,ix2,2)-rminb)/box%dr(2))
              bc(ixbca,ixbcb)=sum(wbn(&
                 (ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,&
                 (ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells,&
                 ixbcn))/dble(2**(2*dlvl))
            end do
          end do
          deallocate(xcoarse)
        end if
      end do
    end select
  end subroutine magnetic_reference_bc
}

end module mod_magnetic_reference_fv
