module mod_amr_solution_node
  use mod_comm_lib, only: mpistop

  implicit none
  private

  public :: getnode, putnode
  public :: alloc_node, alloc_state 
  public :: dealloc_node 
 
contains


  !> Get first available igrid on processor ipe
  integer function getnode(ipe)
    use mod_forest, only: igrid_inuse
    use mod_global_parameters

    integer, intent(in) :: ipe
    integer :: igrid, igrid_available
  
    igrid_available=0
  
    do igrid=1,max_blocks
       if (igrid_inuse(igrid,ipe)) cycle
  
       igrid_available=igrid
       exit
    end do
  
    if (igrid_available == 0) then
       getnode = -1
       print *, "Current maximum number of grid blocks:", max_blocks
       call mpistop&
          ("Insufficient grid blocks; increase max_blocks in meshlist")
    else
       getnode=igrid_available
       igrid_inuse(igrid,ipe)=.true.
    end if
  
    if (ipe==mype) then
       ! initialize node on host and device
       node(1:nodehi,getnode) = 0
       !$acc update device(node(1:nodehi,getnode))
       rnode(1:rnodehi,getnode) = zero
       !$acc update device(rnode(1:rnodehi,getnode))
    end if
  
  end function getnode
  
  ! put igrid on processor ipe to be not in use
  subroutine putnode(igrid,ipe)
    use mod_forest
    implicit none
  
    integer, intent(in) :: igrid, ipe
  
    igrid_inuse(igrid,ipe)=.false.
  
  end subroutine putnode
  
  !> allocate arrays on igrid node
  subroutine alloc_node(igrid)
    use mod_forest
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_set_surface
    use mod_physics, only: phys_set_equi_vars
    use mod_b0, only: set_B0_grid 
    
    integer, intent(in) :: igrid
  
    integer :: level, ig1,ig2,ig3, ign1,ign2,ign3, ixCoGmin1,ixCoGmin2,&
       ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3, ix, i1,i2,i3
    integer :: imin, imax, index, igCo1,igCo2,igCo3, ixshift, offset, ifirst
    integer :: icase, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3
    double precision :: dx1,dx2,dx3, summeddx, sizeuniformpart1,&
       sizeuniformpart2,sizeuniformpart3
    double precision :: xext(ixGlo1-1:ixGhi1+1,ixGlo2-1:ixGhi2+1,&
       ixGlo3-1:ixGhi3+1,1:ndim)
    double precision :: delx_ext(ixGlo1-1:ixGhi1+1)
    double precision :: exp_factor_ext(ixGlo1-1:ixGhi1+1),&
       del_exp_factor_ext(ixGlo1-1:ixGhi1+1),&
       exp_factor_primitive_ext(ixGlo1-1:ixGhi1+1)
    double precision :: xc(ixGlo1:ixGhi1),delxc(ixGlo1:ixGhi1)
    double precision :: exp_factor_coarse(ixGlo1:ixGhi1),&
       del_exp_factor_coarse(ixGlo1:ixGhi1),&
       exp_factor_primitive_coarse(ixGlo1:ixGhi1)

    ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1;
    ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells
    ixCoGmax2=(ixGhi2-2*nghostcells)/2+2*nghostcells
    ixCoGmax3=(ixGhi3-2*nghostcells)/2+2*nghostcells;
  
    icase=mod(nghostcells,2)
    if(stagger_grid) icase=1
    select case(icase)
      case(0)
        ixGextmin1=ixGlo1;ixGextmin2=ixGlo2;ixGextmin3=ixGlo3
        ixGextmax1=ixGhi1;ixGextmax2=ixGhi2;ixGextmax3=ixGhi3;
      case(1)
        ! for ghost cell related prolongations, we need
        ! an extra layer with known volumes and dx-intervals
        ! in case the number of ghost cells is odd
        ixGextmin1=ixGlo1-1;ixGextmin2=ixGlo2-1;ixGextmin3=ixGlo3-1
        ixGextmax1=ixGhi1+1;ixGextmax2=ixGhi2+1;ixGextmax3=ixGhi3+1;
      case default
        call mpistop("no such case")
    end select
  
    ! set level information
    level=igrid_to_node(igrid,mype)%node%level
  
    if(.not. associated(ps(igrid)%w)) then
       
       ! allocate arrays for solution and space
       call alloc_state(igrid, ps(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
          ixGhi3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
          ixGextmax3, .true.)
       ! allocate arrays for one level coarser solution
       call alloc_state_coarse(igrid, psc(igrid), ixCoGmin1,ixCoGmin2,&
          ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3, ixCoGmin1,ixCoGmin2,&
          ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3)
       if(.not.convert) then
        ! allocate arrays for temp solution 1
        call alloc_state(igrid, ps1(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
           ixGhi3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
           ixGextmax3, .false.)
  
        ! allocate temporary solution space
        select case (t_integrator)
        case(ssprk3,ssprk4,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
          call alloc_state(igrid, ps2(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
             ixGextmax2,ixGextmax3, .false.)
        case(RK3_BT,rk4,ssprk5,IMEX_CB3a)
          call alloc_state(igrid, ps2(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
             ixGextmax2,ixGextmax3, .false.)
          call alloc_state(igrid, ps3(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
             ixGextmax2,ixGextmax3, .false.)
        case(IMEX_ARS3,IMEX_232)
          call alloc_state(igrid, ps2(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
             ixGextmax2,ixGextmax3, .false.)
          call alloc_state(igrid, ps3(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
             ixGextmax2,ixGextmax3, .false.)
          call alloc_state(igrid, ps4(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
             ixGextmax2,ixGextmax3, .false.)
        end select
      end if
  
    end if
  
    ! avoid dividing by zero rho in skipped corner ghostcells when phys_req_diagonal=F
    ps(igrid)%w(:,:,:,1)=1.d0
    ps(igrid)%level=level
    psc(igrid)%level=level-1
    ! avoid dividing by zero rho in skipped corner ghostcells when phys_req_diagonal=F
    psc(igrid)%w(:,:,:,1)=1.d0
    if(phys_trac) ps(igrid)%special_values=0.d0
    if(.not.convert) then
      ps1(igrid)%level=level
      select case (t_integrator)
      case(ssprk3,ssprk4,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
        ps2(igrid)%level=level
      case(RK3_BT,rk4,ssprk5,IMEX_CB3a)
        ps2(igrid)%level=level
        ps3(igrid)%level=level
      case(IMEX_ARS3,IMEX_232)
        ps2(igrid)%level=level
        ps3(igrid)%level=level
        ps4(igrid)%level=level
      end select
    end if
  
    ! block pointer to current block
    block=>ps(igrid)
    ig1=igrid_to_node(igrid,mype)%node%ig1
    ig2=igrid_to_node(igrid,mype)%node%ig2
    ig3=igrid_to_node(igrid,mype)%node%ig3;
    node(plevel_,igrid)=level
    node(pig1_,igrid)=ig1
    node(pig2_,igrid)=ig2
    node(pig3_,igrid)=ig3
 !$acc update device(node(plevel_,igrid),node(pig1_,igrid),node(pig2_,igrid),node(pig3_,igrid))
    
    ! set dx information
    rnode(rpdx1_,igrid)=dx(1,level)
    rnode(rpdx2_,igrid)=dx(2,level)
    rnode(rpdx3_,igrid)=dx(3,level)
    dxlevel(:)=dx(:,level)
 !$acc update device(rnode(rpdx1_,igrid),rnode(rpdx2_,igrid),rnode(rpdx3_,igrid), dxlevel)
  
    ! uniform cartesian case as well as all unstretched coordinates
    ! determine the minimal and maximal corners
    rnode(rpxmin1_,igrid)=xprobmin1+dble(ig1-1)*dg1(level)
    rnode(rpxmin2_,igrid)=xprobmin2+dble(ig2-1)*dg2(level)
    rnode(rpxmin3_,igrid)=xprobmin3+dble(ig3-1)*dg3(level)
    rnode(rpxmax1_,igrid)=xprobmin1+dble(ig1)*dg1(level)
    rnode(rpxmax2_,igrid)=xprobmin2+dble(ig2)*dg2(level)
    rnode(rpxmax3_,igrid)=xprobmin3+dble(ig3)*dg3(level)
   if(rnode(rpxmax1_,igrid)>xprobmax1) rnode(rpxmax1_,igrid)=xprobmax1
   if(rnode(rpxmax2_,igrid)>xprobmax2) rnode(rpxmax2_,igrid)=xprobmax2
   if(rnode(rpxmax3_,igrid)>xprobmax3) rnode(rpxmax3_,igrid)=xprobmax3

 !$acc update device( rnode(rpxmax1_,igrid),rnode(rpxmax2_,igrid),rnode(rpxmax3_,igrid), rnode(rpxmin1_,igrid),rnode(rpxmin2_,igrid),rnode(rpxmin3_,igrid) )
   
    dx1=rnode(rpdx1_,igrid)
    dx2=rnode(rpdx2_,igrid)
    dx3=rnode(rpdx3_,igrid)
   do ix=ixGlo1,ixMhi1-nghostcells
      ps(igrid)%x(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1)=rnode(rpxmin1_,&
         igrid)+(dble(ix-nghostcells)-half)*dx1
    end do
   do ix=ixGlo2,ixMhi2-nghostcells
      ps(igrid)%x(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,2)=rnode(rpxmin2_,&
         igrid)+(dble(ix-nghostcells)-half)*dx2
    end do
   do ix=ixGlo3,ixMhi3-nghostcells
      ps(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,3)=rnode(rpxmin3_,&
         igrid)+(dble(ix-nghostcells)-half)*dx3
    end do
   ! update overlap cells of neighboring blocks in the same way to get the same values
   do ix=ixMhi1-nghostcells+1,ixGhi1
      ps(igrid)%x(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1)=rnode(rpxmax1_,&
         igrid)+(dble(ix-ixMhi1)-half)*dx1
    end do
   do ix=ixMhi2-nghostcells+1,ixGhi2
      ps(igrid)%x(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,2)=rnode(rpxmax2_,&
         igrid)+(dble(ix-ixMhi2)-half)*dx2
    end do
   do ix=ixMhi3-nghostcells+1,ixGhi3
      ps(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,3)=rnode(rpxmax3_,&
         igrid)+(dble(ix-ixMhi3)-half)*dx3
    end do
  
    dx1=2.0d0*rnode(rpdx1_,igrid)
    dx2=2.0d0*rnode(rpdx2_,igrid)
    dx3=2.0d0*rnode(rpdx3_,igrid)
   do ix=ixCoGmin1,ixCoGmax1
      psc(igrid)%x(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
         1)=rnode(rpxmin1_,igrid)+(dble(ix-nghostcells)-half)*dx1
    end do
   do ix=ixCoGmin2,ixCoGmax2
      psc(igrid)%x(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
         2)=rnode(rpxmin2_,igrid)+(dble(ix-nghostcells)-half)*dx2
    end do
   do ix=ixCoGmin3,ixCoGmax3
      psc(igrid)%x(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
         3)=rnode(rpxmin3_,igrid)+(dble(ix-nghostcells)-half)*dx3
    end do
  
    ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1)=rnode(rpdx1_,igrid)
    ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,2)=rnode(rpdx2_,igrid)
    ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,3)=rnode(rpdx3_,igrid);
    psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
       1)=2.0d0*rnode(rpdx1_,igrid)
    psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
       2)=2.0d0*rnode(rpdx2_,igrid)
    psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
       3)=2.0d0*rnode(rpdx3_,igrid);
    dx1=rnode(rpdx1_,igrid)
    dx2=rnode(rpdx2_,igrid)
    dx3=rnode(rpdx3_,igrid)
   do ix=ixGextmin1,ixGextmax1
      xext(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)=rnode(rpxmin1_,&
         igrid)+(dble(ix-nghostcells)-half)*dx1
    end do
   do ix=ixGextmin2,ixGextmax2
      xext(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,2)=rnode(rpxmin2_,&
         igrid)+(dble(ix-nghostcells)-half)*dx2
    end do
   do ix=ixGextmin3,ixGextmax3
      xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,3)=rnode(rpxmin3_,&
         igrid)+(dble(ix-nghostcells)-half)*dx3
    end do
  
    if(any(stretched_dim)) then
     if(stretch_type(1) == stretch_uni)then
        imin=(ig1-1)*block_nx1
        imax=ig1*block_nx1
        rnode(rpxmin1_,igrid)=xprobmin1+dxfirst_1mq(level,&
           1) *(1.0d0-qstretch(level,1)**imin)
        rnode(rpxmax1_,igrid)=xprobmin1+dxfirst_1mq(level,&
           1) *(1.0d0-qstretch(level,1)**imax)
        ! fix possible out of bound due to precision
        if(rnode(rpxmax1_,igrid)>xprobmax1) rnode(rpxmax1_,igrid)=xprobmax1
        ixshift=(ig1-1)*block_nx1-nghostcells
        do ix=ixGextmin1,ixGextmax1
          index=ixshift+ix
          ps(igrid)%dx(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
             1)=dxfirst(level,1)*qstretch(level,1)**(index-1)
        enddo
        igCo1=(ig1-1)/2
        ixshift=igCo1*block_nx1+(1-modulo(ig1,2))*block_nx1/2-nghostcells
        do ix=ixCoGmin1,ixCoGmax1
          index=ixshift+ix
          psc(igrid)%dx(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
             1)=dxfirst(level-1,1)*qstretch(level-1,1)**(index-1)
          psc(igrid)%x(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
             1)=xprobmin1+dxfirst_1mq(level-1,1)*(1.0d0-qstretch(level-1,&
             1)**(index-1))+ 0.5d0*dxfirst(level-1,1)*qstretch(level-1,&
             1)**(index-1)
        end do
        ! now that dx and grid boundaries are known: fill cell centers
        ifirst=nghostcells+1
        ! first fill the mesh
        summeddx=0.0d0
        do ix=ixMlo1,ixMhi1
          ps(igrid)%x(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1)=rnode(rpxmin1_,&
             igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
             1)
          summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
        enddo
        ! then ghost cells to left
        summeddx=0.0d0
        do ix=nghostcells,1,-1
          ps(igrid)%x(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1)=rnode(rpxmin1_,&
             igrid)-summeddx-0.5d0*ps(igrid)%dx(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
             1)
          summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
        enddo
        ! then ghost cells to right
        summeddx=0.0d0
        do ix=ixGhi1-nghostcells+1,ixGhi1
          ps(igrid)%x(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1)=rnode(rpxmax1_,&
             igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
             1)
          summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
        enddo
        select case(icase)
          case(0)
            ! if even number of ghost cells: xext is just copy of local x
            xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
               ixGextmin3:ixGextmax3,1)=ps(igrid)%x(ixGextmin1:ixGextmax1,&
               ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
          case(1)
            ! if uneven number of ghost cells: extra layer left/right
            summeddx=0.0d0
            do ix=ixMlo1,ixMhi1
              xext(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                 1)=rnode(rpxmin1_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,&
                 ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
             summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
            enddo
            ! then ghost cells to left
            summeddx=0.0d0
            do ix=nghostcells,ixGextmin1,-1
              xext(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                 1)=rnode(rpxmin1_,igrid)-summeddx-0.5d0*ps(igrid)%dx(ix,&
                 ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
              summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
            enddo
            ! then ghost cells to right
            summeddx=0.0d0
            do ix=ixGhi1-nghostcells+1,ixGextmax1
               xext(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                  1)=rnode(rpxmax1_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,&
                  ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
               summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
            enddo
          case default
            call mpistop("no such case")
        end select
       endif
     if(stretch_type(2) == stretch_uni)then
        imin=(ig2-1)*block_nx2
        imax=ig2*block_nx2
        rnode(rpxmin2_,igrid)=xprobmin2+dxfirst_1mq(level,&
           2) *(1.0d0-qstretch(level,2)**imin)
        rnode(rpxmax2_,igrid)=xprobmin2+dxfirst_1mq(level,&
           2) *(1.0d0-qstretch(level,2)**imax)
        ! fix possible out of bound due to precision
        if(rnode(rpxmax2_,igrid)>xprobmax2) rnode(rpxmax2_,igrid)=xprobmax2
        ixshift=(ig2-1)*block_nx2-nghostcells
        do ix=ixGextmin2,ixGextmax2
          index=ixshift+ix
          ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
             2)=dxfirst(level,2)*qstretch(level,2)**(index-1)
        enddo
        igCo2=(ig2-1)/2
        ixshift=igCo2*block_nx2+(1-modulo(ig2,2))*block_nx2/2-nghostcells
        do ix=ixCoGmin2,ixCoGmax2
          index=ixshift+ix
          psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
             2)=dxfirst(level-1,2)*qstretch(level-1,2)**(index-1)
          psc(igrid)%x(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
             2)=xprobmin2+dxfirst_1mq(level-1,2)*(1.0d0-qstretch(level-1,&
             2)**(index-1))+ 0.5d0*dxfirst(level-1,2)*qstretch(level-1,&
             2)**(index-1)
        end do
        ! now that dx and grid boundaries are known: fill cell centers
        ifirst=nghostcells+1
        ! first fill the mesh
        summeddx=0.0d0
        do ix=ixMlo2,ixMhi2
          ps(igrid)%x(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,2)=rnode(rpxmin2_,&
             igrid)+summeddx+0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,&
             2)
          summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
        enddo
        ! then ghost cells to left
        summeddx=0.0d0
        do ix=nghostcells,1,-1
          ps(igrid)%x(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,2)=rnode(rpxmin2_,&
             igrid)-summeddx-0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,&
             2)
          summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
        enddo
        ! then ghost cells to right
        summeddx=0.0d0
        do ix=ixGhi2-nghostcells+1,ixGhi2
          ps(igrid)%x(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,2)=rnode(rpxmax2_,&
             igrid)+summeddx+0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,&
             2)
          summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
        enddo
        select case(icase)
          case(0)
            ! if even number of ghost cells: xext is just copy of local x
            xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
               ixGextmin3:ixGextmax3,2)=ps(igrid)%x(ixGextmin1:ixGextmax1,&
               ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,2)
          case(1)
            ! if uneven number of ghost cells: extra layer left/right
            summeddx=0.0d0
            do ix=ixMlo2,ixMhi2
              xext(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                 2)=rnode(rpxmin2_,igrid)+summeddx+&
                 0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,&
                 ixGextmin3:ixGextmax3,2)
             summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
            enddo
            ! then ghost cells to left
            summeddx=0.0d0
            do ix=nghostcells,ixGextmin2,-1
              xext(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                 2)=rnode(rpxmin2_,igrid)-summeddx-&
                 0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,&
                 ixGextmin3:ixGextmax3,2)
              summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
            enddo
            ! then ghost cells to right
            summeddx=0.0d0
            do ix=ixGhi2-nghostcells+1,ixGextmax2
               xext(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                  2)=rnode(rpxmax2_,igrid)+summeddx+&
                  0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,&
                  ixGextmin3:ixGextmax3,2)
               summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
            enddo
          case default
            call mpistop("no such case")
        end select
       endif
     if(stretch_type(3) == stretch_uni)then
        imin=(ig3-1)*block_nx3
        imax=ig3*block_nx3
        rnode(rpxmin3_,igrid)=xprobmin3+dxfirst_1mq(level,&
           3) *(1.0d0-qstretch(level,3)**imin)
        rnode(rpxmax3_,igrid)=xprobmin3+dxfirst_1mq(level,&
           3) *(1.0d0-qstretch(level,3)**imax)
        ! fix possible out of bound due to precision
        if(rnode(rpxmax3_,igrid)>xprobmax3) rnode(rpxmax3_,igrid)=xprobmax3
        ixshift=(ig3-1)*block_nx3-nghostcells
        do ix=ixGextmin3,ixGextmax3
          index=ixshift+ix
          ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
             3)=dxfirst(level,3)*qstretch(level,3)**(index-1)
        enddo
        igCo3=(ig3-1)/2
        ixshift=igCo3*block_nx3+(1-modulo(ig3,2))*block_nx3/2-nghostcells
        do ix=ixCoGmin3,ixCoGmax3
          index=ixshift+ix
          psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
             3)=dxfirst(level-1,3)*qstretch(level-1,3)**(index-1)
          psc(igrid)%x(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
             3)=xprobmin3+dxfirst_1mq(level-1,3)*(1.0d0-qstretch(level-1,&
             3)**(index-1))+ 0.5d0*dxfirst(level-1,3)*qstretch(level-1,&
             3)**(index-1)
        end do
        ! now that dx and grid boundaries are known: fill cell centers
        ifirst=nghostcells+1
        ! first fill the mesh
        summeddx=0.0d0
        do ix=ixMlo3,ixMhi3
          ps(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,3)=rnode(rpxmin3_,&
             igrid)+summeddx+0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,&
             3)
          summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
        enddo
        ! then ghost cells to left
        summeddx=0.0d0
        do ix=nghostcells,1,-1
          ps(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,3)=rnode(rpxmin3_,&
             igrid)-summeddx-0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,&
             3)
          summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
        enddo
        ! then ghost cells to right
        summeddx=0.0d0
        do ix=ixGhi3-nghostcells+1,ixGhi3
          ps(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,3)=rnode(rpxmax3_,&
             igrid)+summeddx+0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,&
             3)
          summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
        enddo
        select case(icase)
          case(0)
            ! if even number of ghost cells: xext is just copy of local x
            xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
               ixGextmin3:ixGextmax3,3)=ps(igrid)%x(ixGextmin1:ixGextmax1,&
               ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,3)
          case(1)
            ! if uneven number of ghost cells: extra layer left/right
            summeddx=0.0d0
            do ix=ixMlo3,ixMhi3
              xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                 3)=rnode(rpxmin3_,igrid)+summeddx+&
                 0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
                 ixGextmin2:ixGextmax2,ix,3)
             summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
            enddo
            ! then ghost cells to left
            summeddx=0.0d0
            do ix=nghostcells,ixGextmin3,-1
              xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                 3)=rnode(rpxmin3_,igrid)-summeddx-&
                 0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
                 ixGextmin2:ixGextmax2,ix,3)
              summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
            enddo
            ! then ghost cells to right
            summeddx=0.0d0
            do ix=ixGhi3-nghostcells+1,ixGextmax3
               xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                  3)=rnode(rpxmax3_,igrid)+summeddx+&
                  0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
                  ixGextmin2:ixGextmax2,ix,3)
               summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
            enddo
          case default
            call mpistop("no such case")
        end select
       endif
      if(stretch_type(1) == stretch_symm)then
         ! here we distinguish three kinds of grid blocks
         ! depending on their ig-index, set per level
         !      the first n_stretchedblocks/2  will stretch to the left
         !      the middle ntotal-n_stretchedblocks will be uniform
         !      the last  n_stretchedblocks/2  will stretch to the right
         if(ig1<=nstretchedblocks(level,1)/2)then
           ! stretch to the left
           offset=block_nx1*nstretchedblocks(level,1)/2
           imin=(ig1-1)*block_nx1
           imax=ig1*block_nx1
           rnode(rpxmin1_,igrid)=xprobmin1+xstretch1-dxfirst_1mq(level,&
              1) *(1.0d0-qstretch(level,1)**(offset-imin))
           rnode(rpxmax1_,igrid)=xprobmin1+xstretch1-dxfirst_1mq(level,&
              1) *(1.0d0-qstretch(level,1)**(offset-imax))
           ! fix possible out of bound due to precision
           if(rnode(rpxmin1_,igrid)<xprobmin1) rnode(rpxmin1_,igrid)=xprobmin1
           ixshift=(ig1-1)*block_nx1-nghostcells
           do ix=ixGextmin1,ixGextmax1
             index=ixshift+ix
             ps(igrid)%dx(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                1)=dxfirst(level,1)*qstretch(level,1)**(offset-index)
           enddo
           ixshift=(nstretchedblocks(level,&
              1)/2-ig1)*(block_nx1/2)+block_nx1/2+nghostcells
           do ix=ixCoGmin1,ixCoGmax1
             index=ixshift-ix
             psc(igrid)%dx(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
                1)=dxfirst(level-1,1)*qstretch(level-1,1)**index
           enddo
           ! last block: to modify ghost cells!!!
           if(ig1==nstretchedblocks(level,1)/2)then
             if(ng1(level)==nstretchedblocks(level,1))then
               ! if middle blocks do not exist then use symmetry
               do ix=ixGhi1-nghostcells+1,ixGextmax1
                  ps(igrid)%dx(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                     1)= ps(igrid)%dx(2*(ixGhi1-nghostcells)+1-ix,&
                     ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
               enddo
               do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
                  psc(igrid)%dx(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
                     1)= psc(igrid)%dx(2*(ixCoGmax1-nghostcells)+1-ix,&
                     ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,1)
               enddo
             else
               ! if middle blocks exist then use same as middle blocks:
               do ix=ixGhi1-nghostcells+1,ixGextmax1
                  ps(igrid)%dx(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                     1)=dxmid(level,1)
               enddo
               do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
                  psc(igrid)%dx(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
                     1)=dxmid(level-1,1)
               enddo
             endif
           endif
           ! first block: make ghost cells symmetric (to allow periodicity)
           if(ig1==1)then
             do ix=ixGextmin1,nghostcells
               ps(igrid)%dx(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                  1)=ps(igrid)%dx(2*nghostcells+1-ix,ixGextmin2:ixGextmax2,&
                  ixGextmin3:ixGextmax3,1)
             enddo
             do ix=1,nghostcells
               psc(igrid)%dx(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
                  1)=psc(igrid)%dx(2*nghostcells+1-ix,ixCoGmin2:ixCoGmax2,&
                  ixCoGmin3:ixCoGmax3,1)
             enddo
           endif
         else
           if(ig1<=ng1(level)-nstretchedblocks(level,1)/2) then
             ! keep uniform
             ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
                ixGextmin3:ixGextmax3,1)=dxmid(level,1)
             psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
                ixCoGmin3:ixCoGmax3,1)=dxmid(level-1,1)
             rnode(rpxmin1_,igrid)=xprobmin1+xstretch1+&
                (ig1-nstretchedblocks(level,1)/2-1)*block_nx1*dxmid(level,1)
             rnode(rpxmax1_,igrid)=xprobmin1+xstretch1+&
                (ig1-nstretchedblocks(level,1)/2)  *block_nx1*dxmid(level,1)
             ! first and last block: to modify the ghost cells!!!
             if(ig1==nstretchedblocks(level,1)/2+1)then
               do ix=ixGextmin1,nghostcells
                 ps(igrid)%dx(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                    1)=dxfirst(level,1)*qstretch(level,1)**(nghostcells-ix)
               enddo
               do ix=1,nghostcells
                 psc(igrid)%dx(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
                    1)=dxfirst(level-1,1)*qstretch(level-1,&
                    1)**(nghostcells-ix)
               enddo
             endif
             if(ig1==ng1(level)-nstretchedblocks(level,1))then
               do ix=ixGhi1-nghostcells+1,ixGextmax1
                 ps(igrid)%dx(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                    1)=dxfirst(level,1)*qstretch(level,&
                    1)**(ix-block_nx1-nghostcells-1)
               enddo
               do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
                 psc(igrid)%dx(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
                    1)=dxfirst(level-1,1)*qstretch(level-1,&
                    1)**(ix-ixCoGmax1+nghostcells-1)
               enddo
             endif
           else
             ! stretch to the right
             offset=block_nx1*(ng1(level)-nstretchedblocks(level,1)/2)
             sizeuniformpart1=dxmid(1,1)*(domain_nx1-&
                nstretchedblocks_baselevel(1)*block_nx1)
             imin=(ig1-1)*block_nx1-offset
             imax=ig1*block_nx1-offset
             rnode(rpxmin1_,igrid)=xprobmin1+xstretch1+sizeuniformpart1+&
                dxfirst_1mq(level,1) *(1.0d0-qstretch(level,1)**imin)
             rnode(rpxmax1_,igrid)=xprobmin1+xstretch1+sizeuniformpart1+&
                dxfirst_1mq(level,1) *(1.0d0-qstretch(level,1)**imax)
             ! fix possible out of bound due to precision
             if(rnode(rpxmax1_,igrid)>xprobmax1) rnode(rpxmax1_,&
                igrid)=xprobmax1
             ixshift=(ig1-1)*block_nx1-nghostcells-offset
             do ix=ixGextmin1,ixGextmax1
               index=ixshift+ix
               ps(igrid)%dx(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                  1)=dxfirst(level,1)*qstretch(level,1)**(index-1)
             enddo
             ixshift=(ig1+nstretchedblocks(level,&
                1)/2-ng1(level)-1)*(block_nx1/2)-nghostcells
             do ix=ixCoGmin1,ixCoGmax1
               index=ixshift+ix
               psc(igrid)%dx(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
                  1)=dxfirst(level-1,1)*qstretch(level-1,1)**(index-1)
             enddo
             ! first block: modify ghost cells!!!
             if(ig1==ng1(level)-nstretchedblocks(level,1)+1)then
               if(ng1(level)==nstretchedblocks(level,1))then
                 ! if middle blocks do not exist then use symmetry
                 do ix=ixGextmin1,nghostcells
                   ps(igrid)%dx(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                      1)=ps(igrid)%dx(2*nghostcells+1-ix,ixGextmin2:ixGextmax2,&
                      ixGextmin3:ixGextmax3,1)
                 enddo
                 do ix=1,nghostcells
                   psc(igrid)%dx(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
                      1)=psc(igrid)%dx(2*nghostcells+1-ix,ixCoGmin2:ixCoGmax2,&
                      ixCoGmin3:ixCoGmax3,1)
                 enddo
               else
                 ! if middle blocks exist then use same as middle blocks:
                 do ix=ixGextmin1,nghostcells
                   ps(igrid)%dx(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                      1)=dxmid(level,1)
                 enddo
                 do ix=1,nghostcells
                   psc(igrid)%dx(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
                      1)=dxmid(level-1,1)
                 enddo
               endif
             endif
             ! last block: make ghost cells symmetric (to allow periodicity)
             if(ig1==ng1(level))then
               do ix=ixGhi1-nghostcells+1,ixGextmax1
                 ps(igrid)%dx(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                    1)=ps(igrid)%dx(2*(ixGhi1-nghostcells)+1-ix,&
                    ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
               enddo
               do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
                 psc(igrid)%dx(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
                    1)=psc(igrid)%dx(2*(ixCoGmax1-nghostcells)+1-ix,&
                    ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,1)
               enddo
             endif
           endif
         endif
         ! now that dx and grid boundaries are known: fill cell centers
         ifirst=nghostcells+1
         ! first fill the mesh
         summeddx=0.0d0
         do ix=ixMlo1,ixMhi1
           ps(igrid)%x(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1)=rnode(rpxmin1_,&
              igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,ixGlo2:ixGhi2,&
              ixGlo3:ixGhi3,1)
           summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
         enddo
         ! then ghost cells to left
         summeddx=0.0d0
         do ix=nghostcells,1,-1
           ps(igrid)%x(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1)=rnode(rpxmin1_,&
              igrid)-summeddx-0.5d0*ps(igrid)%dx(ix,ixGlo2:ixGhi2,&
              ixGlo3:ixGhi3,1)
           summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
         enddo
         ! then ghost cells to right
         summeddx=0.0d0
         do ix=ixGhi1-nghostcells+1,ixGhi1
           ps(igrid)%x(ix,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1)=rnode(rpxmax1_,&
              igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,ixGlo2:ixGhi2,&
              ixGlo3:ixGhi3,1)
           summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
         enddo
         ! and next for the coarse representation
         ! first fill the mesh
         summeddx=0.0d0
         do ix=nghostcells+1,ixCoGmax1-nghostcells
           psc(igrid)%x(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
              1)=rnode(rpxmin1_,igrid)+summeddx+0.5d0*psc(igrid)%dx(ix,&
              ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,1)
           summeddx=summeddx+psc(igrid)%dx(ix,ifirst,ifirst,1)
         enddo
         ! then ghost cells to left
         summeddx=0.0d0
         do ix=nghostcells,1,-1
           psc(igrid)%x(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
              1)=rnode(rpxmin1_,igrid)-summeddx-0.5d0*psc(igrid)%dx(ix,&
              ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,1)
           summeddx=summeddx+psc(igrid)%dx(ix,ifirst,ifirst,1)
         enddo
         ! then ghost cells to right
         summeddx=0.0d0
         do ix=ixCoGmax1-nghostcells+1,ixCoGmax1
           psc(igrid)%x(ix,ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
              1)=rnode(rpxmax1_,igrid)+summeddx+0.5d0*psc(igrid)%dx(ix,&
              ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,1)
           summeddx=summeddx+psc(igrid)%dx(ix,ifirst,ifirst,1)
         enddo
         select case(icase)
           case(0)
             ! if even number of ghost cells: xext is just copy of local x
             xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
                ixGextmin3:ixGextmax3,1)=ps(igrid)%x(ixGextmin1:ixGextmax1,&
                ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
           case(1)
             ! if uneven number of ghost cells: extra layer left/right
             summeddx=0.0d0
             do ix=ixMlo1,ixMhi1
               xext(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                  1)=rnode(rpxmin1_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,&
                  ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
              summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
             enddo
             ! then ghost cells to left
             summeddx=0.0d0
             do ix=nghostcells,ixGextmin1,-1
               xext(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                  1)=rnode(rpxmin1_,igrid)-summeddx-0.5d0*ps(igrid)%dx(ix,&
                  ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
               summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
             enddo
            ! then ghost cells to right
            summeddx=0.0d0
            do ix=ixGhi1-nghostcells+1,ixGextmax1
               xext(ix,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
                  1)=rnode(rpxmax1_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix,&
                  ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
               summeddx=summeddx+ps(igrid)%dx(ix,ifirst,ifirst,1)
            enddo
           case default
             call mpistop("no such case")
         end select
       endif
      if(stretch_type(2) == stretch_symm)then
         ! here we distinguish three kinds of grid blocks
         ! depending on their ig-index, set per level
         !      the first n_stretchedblocks/2  will stretch to the left
         !      the middle ntotal-n_stretchedblocks will be uniform
         !      the last  n_stretchedblocks/2  will stretch to the right
         if(ig2<=nstretchedblocks(level,2)/2)then
           ! stretch to the left
           offset=block_nx2*nstretchedblocks(level,2)/2
           imin=(ig2-1)*block_nx2
           imax=ig2*block_nx2
           rnode(rpxmin2_,igrid)=xprobmin2+xstretch2-dxfirst_1mq(level,&
              2) *(1.0d0-qstretch(level,2)**(offset-imin))
           rnode(rpxmax2_,igrid)=xprobmin2+xstretch2-dxfirst_1mq(level,&
              2) *(1.0d0-qstretch(level,2)**(offset-imax))
           ! fix possible out of bound due to precision
           if(rnode(rpxmin2_,igrid)<xprobmin2) rnode(rpxmin2_,igrid)=xprobmin2
           ixshift=(ig2-1)*block_nx2-nghostcells
           do ix=ixGextmin2,ixGextmax2
             index=ixshift+ix
             ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                2)=dxfirst(level,2)*qstretch(level,2)**(offset-index)
           enddo
           ixshift=(nstretchedblocks(level,&
              2)/2-ig2)*(block_nx2/2)+block_nx2/2+nghostcells
           do ix=ixCoGmin2,ixCoGmax2
             index=ixshift-ix
             psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
                2)=dxfirst(level-1,2)*qstretch(level-1,2)**index
           enddo
           ! last block: to modify ghost cells!!!
           if(ig2==nstretchedblocks(level,2)/2)then
             if(ng2(level)==nstretchedblocks(level,2))then
               ! if middle blocks do not exist then use symmetry
               do ix=ixGhi2-nghostcells+1,ixGextmax2
                  ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                     2)= ps(igrid)%dx(ixGextmin1:ixGextmax1,&
                     2*(ixGhi2-nghostcells)+1-ix,ixGextmin3:ixGextmax3,2)
               enddo
               do ix=ixCoGmax2-nghostcells+1,ixCoGmax2
                  psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
                     2)= psc(igrid)%dx(ixCoGmin1:ixCoGmax1,&
                     2*(ixCoGmax2-nghostcells)+1-ix,ixCoGmin3:ixCoGmax3,2)
               enddo
             else
               ! if middle blocks exist then use same as middle blocks:
               do ix=ixGhi2-nghostcells+1,ixGextmax2
                  ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                     2)=dxmid(level,2)
               enddo
               do ix=ixCoGmax2-nghostcells+1,ixCoGmax2
                  psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
                     2)=dxmid(level-1,2)
               enddo
             endif
           endif
           ! first block: make ghost cells symmetric (to allow periodicity)
           if(ig2==1)then
             do ix=ixGextmin2,nghostcells
               ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                  2)=ps(igrid)%dx(ixGextmin1:ixGextmax1,2*nghostcells+1-ix,&
                  ixGextmin3:ixGextmax3,2)
             enddo
             do ix=1,nghostcells
               psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
                  2)=psc(igrid)%dx(ixCoGmin1:ixCoGmax1,2*nghostcells+1-ix,&
                  ixCoGmin3:ixCoGmax3,2)
             enddo
           endif
         else
           if(ig2<=ng2(level)-nstretchedblocks(level,2)/2) then
             ! keep uniform
             ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
                ixGextmin3:ixGextmax3,2)=dxmid(level,2)
             psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
                ixCoGmin3:ixCoGmax3,2)=dxmid(level-1,2)
             rnode(rpxmin2_,igrid)=xprobmin2+xstretch2+&
                (ig2-nstretchedblocks(level,2)/2-1)*block_nx2*dxmid(level,2)
             rnode(rpxmax2_,igrid)=xprobmin2+xstretch2+&
                (ig2-nstretchedblocks(level,2)/2)  *block_nx2*dxmid(level,2)
             ! first and last block: to modify the ghost cells!!!
             if(ig2==nstretchedblocks(level,2)/2+1)then
               do ix=ixGextmin2,nghostcells
                 ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                    2)=dxfirst(level,2)*qstretch(level,2)**(nghostcells-ix)
               enddo
               do ix=1,nghostcells
                 psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
                    2)=dxfirst(level-1,2)*qstretch(level-1,&
                    2)**(nghostcells-ix)
               enddo
             endif
             if(ig2==ng2(level)-nstretchedblocks(level,2))then
               do ix=ixGhi2-nghostcells+1,ixGextmax2
                 ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                    2)=dxfirst(level,2)*qstretch(level,&
                    2)**(ix-block_nx2-nghostcells-1)
               enddo
               do ix=ixCoGmax2-nghostcells+1,ixCoGmax2
                 psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
                    2)=dxfirst(level-1,2)*qstretch(level-1,&
                    2)**(ix-ixCoGmax2+nghostcells-1)
               enddo
             endif
           else
             ! stretch to the right
             offset=block_nx2*(ng2(level)-nstretchedblocks(level,2)/2)
             sizeuniformpart2=dxmid(1,2)*(domain_nx2-&
                nstretchedblocks_baselevel(2)*block_nx2)
             imin=(ig2-1)*block_nx2-offset
             imax=ig2*block_nx2-offset
             rnode(rpxmin2_,igrid)=xprobmin2+xstretch2+sizeuniformpart2+&
                dxfirst_1mq(level,2) *(1.0d0-qstretch(level,2)**imin)
             rnode(rpxmax2_,igrid)=xprobmin2+xstretch2+sizeuniformpart2+&
                dxfirst_1mq(level,2) *(1.0d0-qstretch(level,2)**imax)
             ! fix possible out of bound due to precision
             if(rnode(rpxmax2_,igrid)>xprobmax2) rnode(rpxmax2_,&
                igrid)=xprobmax2
             ixshift=(ig2-1)*block_nx2-nghostcells-offset
             do ix=ixGextmin2,ixGextmax2
               index=ixshift+ix
               ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                  2)=dxfirst(level,2)*qstretch(level,2)**(index-1)
             enddo
             ixshift=(ig2+nstretchedblocks(level,&
                2)/2-ng2(level)-1)*(block_nx2/2)-nghostcells
             do ix=ixCoGmin2,ixCoGmax2
               index=ixshift+ix
               psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
                  2)=dxfirst(level-1,2)*qstretch(level-1,2)**(index-1)
             enddo
             ! first block: modify ghost cells!!!
             if(ig2==ng2(level)-nstretchedblocks(level,2)+1)then
               if(ng2(level)==nstretchedblocks(level,2))then
                 ! if middle blocks do not exist then use symmetry
                 do ix=ixGextmin2,nghostcells
                   ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                      2)=ps(igrid)%dx(ixGextmin1:ixGextmax1,2*nghostcells+1-ix,&
                      ixGextmin3:ixGextmax3,2)
                 enddo
                 do ix=1,nghostcells
                   psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
                      2)=psc(igrid)%dx(ixCoGmin1:ixCoGmax1,2*nghostcells+1-ix,&
                      ixCoGmin3:ixCoGmax3,2)
                 enddo
               else
                 ! if middle blocks exist then use same as middle blocks:
                 do ix=ixGextmin2,nghostcells
                   ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                      2)=dxmid(level,2)
                 enddo
                 do ix=1,nghostcells
                   psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
                      2)=dxmid(level-1,2)
                 enddo
               endif
             endif
             ! last block: make ghost cells symmetric (to allow periodicity)
             if(ig2==ng2(level))then
               do ix=ixGhi2-nghostcells+1,ixGextmax2
                 ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                    2)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
                    2*(ixGhi2-nghostcells)+1-ix,ixGextmin3:ixGextmax3,2)
               enddo
               do ix=ixCoGmax2-nghostcells+1,ixCoGmax2
                 psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
                    2)=psc(igrid)%dx(ixCoGmin1:ixCoGmax1,&
                    2*(ixCoGmax2-nghostcells)+1-ix,ixCoGmin3:ixCoGmax3,2)
               enddo
             endif
           endif
         endif
         ! now that dx and grid boundaries are known: fill cell centers
         ifirst=nghostcells+1
         ! first fill the mesh
         summeddx=0.0d0
         do ix=ixMlo2,ixMhi2
           ps(igrid)%x(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,2)=rnode(rpxmin2_,&
              igrid)+summeddx+0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ix,&
              ixGlo3:ixGhi3,2)
           summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
         enddo
         ! then ghost cells to left
         summeddx=0.0d0
         do ix=nghostcells,1,-1
           ps(igrid)%x(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,2)=rnode(rpxmin2_,&
              igrid)-summeddx-0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ix,&
              ixGlo3:ixGhi3,2)
           summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
         enddo
         ! then ghost cells to right
         summeddx=0.0d0
         do ix=ixGhi2-nghostcells+1,ixGhi2
           ps(igrid)%x(ixGlo1:ixGhi1,ix,ixGlo3:ixGhi3,2)=rnode(rpxmax2_,&
              igrid)+summeddx+0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ix,&
              ixGlo3:ixGhi3,2)
           summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
         enddo
         ! and next for the coarse representation
         ! first fill the mesh
         summeddx=0.0d0
         do ix=nghostcells+1,ixCoGmax2-nghostcells
           psc(igrid)%x(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
              2)=rnode(rpxmin2_,igrid)+summeddx+&
              0.5d0*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
              2)
           summeddx=summeddx+psc(igrid)%dx(ifirst,ix,ifirst,2)
         enddo
         ! then ghost cells to left
         summeddx=0.0d0
         do ix=nghostcells,1,-1
           psc(igrid)%x(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
              2)=rnode(rpxmin2_,igrid)-summeddx-&
              0.5d0*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
              2)
           summeddx=summeddx+psc(igrid)%dx(ifirst,ix,ifirst,2)
         enddo
         ! then ghost cells to right
         summeddx=0.0d0
         do ix=ixCoGmax2-nghostcells+1,ixCoGmax2
           psc(igrid)%x(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
              2)=rnode(rpxmax2_,igrid)+summeddx+&
              0.5d0*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ix,ixCoGmin3:ixCoGmax3,&
              2)
           summeddx=summeddx+psc(igrid)%dx(ifirst,ix,ifirst,2)
         enddo
         select case(icase)
           case(0)
             ! if even number of ghost cells: xext is just copy of local x
             xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
                ixGextmin3:ixGextmax3,2)=ps(igrid)%x(ixGextmin1:ixGextmax1,&
                ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,2)
           case(1)
             ! if uneven number of ghost cells: extra layer left/right
             summeddx=0.0d0
             do ix=ixMlo2,ixMhi2
               xext(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                  2)=rnode(rpxmin2_,igrid)+summeddx+&
                  0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,&
                  ixGextmin3:ixGextmax3,2)
              summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
             enddo
             ! then ghost cells to left
             summeddx=0.0d0
             do ix=nghostcells,ixGextmin2,-1
               xext(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                  2)=rnode(rpxmin2_,igrid)-summeddx-&
                  0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,&
                  ixGextmin3:ixGextmax3,2)
               summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
             enddo
            ! then ghost cells to right
            summeddx=0.0d0
            do ix=ixGhi2-nghostcells+1,ixGextmax2
               xext(ixGextmin1:ixGextmax1,ix,ixGextmin3:ixGextmax3,&
                  2)=rnode(rpxmax2_,igrid)+summeddx+&
                  0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,ix,&
                  ixGextmin3:ixGextmax3,2)
               summeddx=summeddx+ps(igrid)%dx(ifirst,ix,ifirst,2)
            enddo
           case default
             call mpistop("no such case")
         end select
       endif
      if(stretch_type(3) == stretch_symm)then
         ! here we distinguish three kinds of grid blocks
         ! depending on their ig-index, set per level
         !      the first n_stretchedblocks/2  will stretch to the left
         !      the middle ntotal-n_stretchedblocks will be uniform
         !      the last  n_stretchedblocks/2  will stretch to the right
         if(ig3<=nstretchedblocks(level,3)/2)then
           ! stretch to the left
           offset=block_nx3*nstretchedblocks(level,3)/2
           imin=(ig3-1)*block_nx3
           imax=ig3*block_nx3
           rnode(rpxmin3_,igrid)=xprobmin3+xstretch3-dxfirst_1mq(level,&
              3) *(1.0d0-qstretch(level,3)**(offset-imin))
           rnode(rpxmax3_,igrid)=xprobmin3+xstretch3-dxfirst_1mq(level,&
              3) *(1.0d0-qstretch(level,3)**(offset-imax))
           ! fix possible out of bound due to precision
           if(rnode(rpxmin3_,igrid)<xprobmin3) rnode(rpxmin3_,igrid)=xprobmin3
           ixshift=(ig3-1)*block_nx3-nghostcells
           do ix=ixGextmin3,ixGextmax3
             index=ixshift+ix
             ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                3)=dxfirst(level,3)*qstretch(level,3)**(offset-index)
           enddo
           ixshift=(nstretchedblocks(level,&
              3)/2-ig3)*(block_nx3/2)+block_nx3/2+nghostcells
           do ix=ixCoGmin3,ixCoGmax3
             index=ixshift-ix
             psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
                3)=dxfirst(level-1,3)*qstretch(level-1,3)**index
           enddo
           ! last block: to modify ghost cells!!!
           if(ig3==nstretchedblocks(level,3)/2)then
             if(ng3(level)==nstretchedblocks(level,3))then
               ! if middle blocks do not exist then use symmetry
               do ix=ixGhi3-nghostcells+1,ixGextmax3
                  ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                     3)= ps(igrid)%dx(ixGextmin1:ixGextmax1,&
                     ixGextmin2:ixGextmax2,2*(ixGhi3-nghostcells)+1-ix,3)
               enddo
               do ix=ixCoGmax3-nghostcells+1,ixCoGmax3
                  psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
                     3)= psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
                     2*(ixCoGmax3-nghostcells)+1-ix,3)
               enddo
             else
               ! if middle blocks exist then use same as middle blocks:
               do ix=ixGhi3-nghostcells+1,ixGextmax3
                  ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                     3)=dxmid(level,3)
               enddo
               do ix=ixCoGmax3-nghostcells+1,ixCoGmax3
                  psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
                     3)=dxmid(level-1,3)
               enddo
             endif
           endif
           ! first block: make ghost cells symmetric (to allow periodicity)
           if(ig3==1)then
             do ix=ixGextmin3,nghostcells
               ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                  3)=ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
                  2*nghostcells+1-ix,3)
             enddo
             do ix=1,nghostcells
               psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
                  3)=psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
                  2*nghostcells+1-ix,3)
             enddo
           endif
         else
           if(ig3<=ng3(level)-nstretchedblocks(level,3)/2) then
             ! keep uniform
             ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
                ixGextmin3:ixGextmax3,3)=dxmid(level,3)
             psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
                ixCoGmin3:ixCoGmax3,3)=dxmid(level-1,3)
             rnode(rpxmin3_,igrid)=xprobmin3+xstretch3+&
                (ig3-nstretchedblocks(level,3)/2-1)*block_nx3*dxmid(level,3)
             rnode(rpxmax3_,igrid)=xprobmin3+xstretch3+&
                (ig3-nstretchedblocks(level,3)/2)  *block_nx3*dxmid(level,3)
             ! first and last block: to modify the ghost cells!!!
             if(ig3==nstretchedblocks(level,3)/2+1)then
               do ix=ixGextmin3,nghostcells
                 ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                    3)=dxfirst(level,3)*qstretch(level,3)**(nghostcells-ix)
               enddo
               do ix=1,nghostcells
                 psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
                    3)=dxfirst(level-1,3)*qstretch(level-1,&
                    3)**(nghostcells-ix)
               enddo
             endif
             if(ig3==ng3(level)-nstretchedblocks(level,3))then
               do ix=ixGhi3-nghostcells+1,ixGextmax3
                 ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                    3)=dxfirst(level,3)*qstretch(level,&
                    3)**(ix-block_nx3-nghostcells-1)
               enddo
               do ix=ixCoGmax3-nghostcells+1,ixCoGmax3
                 psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
                    3)=dxfirst(level-1,3)*qstretch(level-1,&
                    3)**(ix-ixCoGmax3+nghostcells-1)
               enddo
             endif
           else
             ! stretch to the right
             offset=block_nx3*(ng3(level)-nstretchedblocks(level,3)/2)
             sizeuniformpart3=dxmid(1,3)*(domain_nx3-&
                nstretchedblocks_baselevel(3)*block_nx3)
             imin=(ig3-1)*block_nx3-offset
             imax=ig3*block_nx3-offset
             rnode(rpxmin3_,igrid)=xprobmin3+xstretch3+sizeuniformpart3+&
                dxfirst_1mq(level,3) *(1.0d0-qstretch(level,3)**imin)
             rnode(rpxmax3_,igrid)=xprobmin3+xstretch3+sizeuniformpart3+&
                dxfirst_1mq(level,3) *(1.0d0-qstretch(level,3)**imax)
             ! fix possible out of bound due to precision
             if(rnode(rpxmax3_,igrid)>xprobmax3) rnode(rpxmax3_,&
                igrid)=xprobmax3
             ixshift=(ig3-1)*block_nx3-nghostcells-offset
             do ix=ixGextmin3,ixGextmax3
               index=ixshift+ix
               ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                  3)=dxfirst(level,3)*qstretch(level,3)**(index-1)
             enddo
             ixshift=(ig3+nstretchedblocks(level,&
                3)/2-ng3(level)-1)*(block_nx3/2)-nghostcells
             do ix=ixCoGmin3,ixCoGmax3
               index=ixshift+ix
               psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
                  3)=dxfirst(level-1,3)*qstretch(level-1,3)**(index-1)
             enddo
             ! first block: modify ghost cells!!!
             if(ig3==ng3(level)-nstretchedblocks(level,3)+1)then
               if(ng3(level)==nstretchedblocks(level,3))then
                 ! if middle blocks do not exist then use symmetry
                 do ix=ixGextmin3,nghostcells
                   ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                      3)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
                      ixGextmin2:ixGextmax2,2*nghostcells+1-ix,3)
                 enddo
                 do ix=1,nghostcells
                   psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
                      3)=psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
                      2*nghostcells+1-ix,3)
                 enddo
               else
                 ! if middle blocks exist then use same as middle blocks:
                 do ix=ixGextmin3,nghostcells
                   ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                      3)=dxmid(level,3)
                 enddo
                 do ix=1,nghostcells
                   psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
                      3)=dxmid(level-1,3)
                 enddo
               endif
             endif
             ! last block: make ghost cells symmetric (to allow periodicity)
             if(ig3==ng3(level))then
               do ix=ixGhi3-nghostcells+1,ixGextmax3
                 ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                    3)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
                    ixGextmin2:ixGextmax2,2*(ixGhi3-nghostcells)+1-ix,3)
               enddo
               do ix=ixCoGmax3-nghostcells+1,ixCoGmax3
                 psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
                    3)=psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
                    2*(ixCoGmax3-nghostcells)+1-ix,3)
               enddo
             endif
           endif
         endif
         ! now that dx and grid boundaries are known: fill cell centers
         ifirst=nghostcells+1
         ! first fill the mesh
         summeddx=0.0d0
         do ix=ixMlo3,ixMhi3
           ps(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,3)=rnode(rpxmin3_,&
              igrid)+summeddx+0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
              ix,3)
           summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
         enddo
         ! then ghost cells to left
         summeddx=0.0d0
         do ix=nghostcells,1,-1
           ps(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,3)=rnode(rpxmin3_,&
              igrid)-summeddx-0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
              ix,3)
           summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
         enddo
         ! then ghost cells to right
         summeddx=0.0d0
         do ix=ixGhi3-nghostcells+1,ixGhi3
           ps(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ix,3)=rnode(rpxmax3_,&
              igrid)+summeddx+0.5d0*ps(igrid)%dx(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
              ix,3)
           summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
         enddo
         ! and next for the coarse representation
         ! first fill the mesh
         summeddx=0.0d0
         do ix=nghostcells+1,ixCoGmax3-nghostcells
           psc(igrid)%x(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
              3)=rnode(rpxmin3_,igrid)+summeddx+&
              0.5d0*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
              3)
           summeddx=summeddx+psc(igrid)%dx(ifirst,ifirst,ix,3)
         enddo
         ! then ghost cells to left
         summeddx=0.0d0
         do ix=nghostcells,1,-1
           psc(igrid)%x(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
              3)=rnode(rpxmin3_,igrid)-summeddx-&
              0.5d0*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
              3)
           summeddx=summeddx+psc(igrid)%dx(ifirst,ifirst,ix,3)
         enddo
         ! then ghost cells to right
         summeddx=0.0d0
         do ix=ixCoGmax3-nghostcells+1,ixCoGmax3
           psc(igrid)%x(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
              3)=rnode(rpxmax3_,igrid)+summeddx+&
              0.5d0*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,ix,&
              3)
           summeddx=summeddx+psc(igrid)%dx(ifirst,ifirst,ix,3)
         enddo
         select case(icase)
           case(0)
             ! if even number of ghost cells: xext is just copy of local x
             xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
                ixGextmin3:ixGextmax3,3)=ps(igrid)%x(ixGextmin1:ixGextmax1,&
                ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,3)
           case(1)
             ! if uneven number of ghost cells: extra layer left/right
             summeddx=0.0d0
             do ix=ixMlo3,ixMhi3
               xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                  3)=rnode(rpxmin3_,igrid)+summeddx+&
                  0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
                  ixGextmin2:ixGextmax2,ix,3)
              summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
             enddo
             ! then ghost cells to left
             summeddx=0.0d0
             do ix=nghostcells,ixGextmin3,-1
               xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                  3)=rnode(rpxmin3_,igrid)-summeddx-&
                  0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
                  ixGextmin2:ixGextmax2,ix,3)
               summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
             enddo
            ! then ghost cells to right
            summeddx=0.0d0
            do ix=ixGhi3-nghostcells+1,ixGextmax3
               xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ix,&
                  3)=rnode(rpxmax3_,igrid)+summeddx+&
                  0.5d0*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
                  ixGextmin2:ixGextmax2,ix,3)
               summeddx=summeddx+ps(igrid)%dx(ifirst,ifirst,ix,3)
            enddo
           case default
             call mpistop("no such case")
         end select
       endif
    endif
  
    ! calculate area of cell surfaces for standard block
    call get_surface_area(ps(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
    ! calculate area of cell surfaces for coarser representative block
    call get_surface_area(psc(igrid),ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
       ixCoGmax2,ixCoGmax3)
    ! calculate volume and distance of cells
    ps(igrid)%dsC=1.d0
    select case (coordinate)
      case (Cartesian)
        ps(igrid)%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3)= rnode(rpdx1_,igrid)*rnode(rpdx2_,&
           igrid)*rnode(rpdx3_,igrid)
        ps(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,1:ndim)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1:ndim)
        ps(igrid)%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,1:ndim)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1:ndim)
        psc(igrid)%dvolume(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
           ixCoGmin3:ixCoGmax3)= 2.d0*rnode(rpdx1_,igrid)*2.d0*rnode(rpdx2_,&
           igrid)*2.d0*rnode(rpdx3_,igrid)
        psc(igrid)%ds(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
           ixCoGmin3:ixCoGmax3,1:ndim)=psc(igrid)%dx(ixCoGmin1:ixCoGmax1,&
           ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,1:ndim)
      case (Cartesian_stretched)
        ps(igrid)%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3)= ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
           1)*ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,2)*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,3)
        ps(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,1:ndim)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1:ndim)
        ps(igrid)%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,1:ndim)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1:ndim)
        psc(igrid)%dvolume(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
           ixCoGmin3:ixCoGmax3)= psc(igrid)%dx(ixCoGmin1:ixCoGmax1,&
           ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
           1)*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
           ixCoGmin3:ixCoGmax3,2)*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,&
           ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,3)
        psc(igrid)%ds(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
           ixCoGmin3:ixCoGmax3,1:ndim)=psc(igrid)%dx(ixCoGmin1:ixCoGmax1,&
           ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,1:ndim)
      case (Cartesian_expansion)
        
      case (spherical)
        ps(igrid)%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3)=(xext(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
           1)**2 +ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,1)**2/12.0d0)*ps(igrid)%dx(&
           ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
           1) *two*dabs(dsin(xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,2))) *dsin(half*ps(igrid)%dx(&
           ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
           2))*ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,3)
        psc(igrid)%dvolume(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
           ixCoGmin3:ixCoGmax3)=(psc(igrid)%x(ixCoGmin1:ixCoGmax1,&
           ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
           1)**2 +psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
           ixCoGmin3:ixCoGmax3,1)**2/12.0d0)*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,&
           ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
           1) *two*dabs(dsin(psc(igrid)%x(ixCoGmin1:ixCoGmax1,&
           ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
           2))) *dsin(half*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,&
           ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
           2))*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
           ixCoGmin3:ixCoGmax3,3)
        ps(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,1)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
           ps(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
              ixGextmin3:ixGextmax3,2)=xext(ixGextmin1:ixGextmax1,&
              ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
              1)*ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
              ixGextmin3:ixGextmax3,2)
         ps(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            ixGextmin3:ixGextmax3,3)=xext(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
            1)*dsin(xext(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            ixGextmin3:ixGextmax3,2))*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,3)
        ps(igrid)%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,1)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
           ps(igrid)%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
              ixGextmin3:ixGextmax3,2)=(xext(ixGextmin1:ixGextmax1,&
              ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
              1)+half*ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
              ixGextmin3:ixGextmax3,1))*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
              ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,2)
        if(ndir>ndim) then
          ps(igrid)%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
             ixGextmin3:ixGextmax3,3)=(xext(ixGextmin1:ixGextmax1,&
             ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
             1)+half*ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
             ixGextmin3:ixGextmax3,1))*dsin(xext(ixGextmin1:ixGextmax1,&
             ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
             2)+half*ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
             ixGextmin3:ixGextmax3,2))
        end if
       
         ps(igrid)%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            ixGextmin3:ixGextmax3,3)=(xext(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
            1)+half*ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            ixGextmin3:ixGextmax3,1))*dsin(xext(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
            2)+half*ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
            ixGextmin3:ixGextmax3,2))*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
            ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,3)
      case (cylindrical)
        ps(igrid)%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3)=dabs(xext(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
           1)) *ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,1)*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
           2) *ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,3)
        psc(igrid)%dvolume(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
           ixCoGmin3:ixCoGmax3)=dabs(psc(igrid)%x(ixCoGmin1:ixCoGmax1,&
           ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
           1)) *psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
           ixCoGmin3:ixCoGmax3,1)*psc(igrid)%dx(ixCoGmin1:ixCoGmax1,&
           ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,&
           2) *psc(igrid)%dx(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
           ixCoGmin3:ixCoGmax3,3)
        ps(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,r_)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,r_)
        ps(igrid)%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
           ixGextmin3:ixGextmax3,r_)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
           ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,r_)
        if(z_>0.and.z_<=ndim) then
          ps(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
             ixGextmin3:ixGextmax3,z_)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
             ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,z_)
          ps(igrid)%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
             ixGextmin3:ixGextmax3,z_)=ps(igrid)%dx(ixGextmin1:ixGextmax1,&
             ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,z_)
          if(phi_>z_.and.ndir>ndim) then
            ps(igrid)%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
               ixGextmin3:ixGextmax3,phi_)=xext(ixGextmin1:ixGextmax1,&
               ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
               1)+half*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
               ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1)
          end if
        end if
        if(phi_>0.and.phi_<=ndim) then
          ps(igrid)%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
             ixGextmin3:ixGextmax3,phi_)=xext(ixGextmin1:ixGextmax1,&
             ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
             1)*ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
             ixGextmin3:ixGextmax3,phi_)
          ps(igrid)%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
             ixGextmin3:ixGextmax3,phi_)=(xext(ixGextmin1:ixGextmax1,&
             ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
             1)+half*ps(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
             ixGextmin3:ixGextmax3,1))*ps(igrid)%dx(ixGextmin1:ixGextmax1,&
             ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,phi_)
          if(z_>phi_.and.ndir>ndim) ps(igrid)%dsC(ixGextmin1:ixGextmax1,&
             ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,z_)=1.d0
        end if
      case default
        call mpistop("Sorry, coordinate unknown")
    end select
  
    ! initialize background non-evolving solution
    if (B0field) call set_B0_grid(igrid)
    if (number_equi_vars>0) call phys_set_equi_vars(igrid)
  
    ! find the blocks on the boundaries
    ps(igrid)%is_physical_boundary=.false.
    
    do i1=-1,1
      if(i1==0) cycle
      ign1=ig1+i1
      ! blocks at periodic boundary have neighbors in the physical domain
      ! thus threated at internal blocks with no physical boundary
      if (periodB(1)) ign1=1+modulo(ign1-1,ng1(level))
      if (ign1 > ng1(level)) then
         if(phi_ > 0 .and. poleB(2,1)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*1)=.false.
         else
           ps(igrid)%is_physical_boundary(2*1)=.true.
         end if
      else if (ign1 < 1) then
         if(phi_ > 0 .and. poleB(1,1)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*1-1)=.false.
         else
           ps(igrid)%is_physical_boundary(2*1-1)=.true.
         end if
      end if
    end do
    
    
    do i2=-1,1
      if(i2==0) cycle
      ign2=ig2+i2
      ! blocks at periodic boundary have neighbors in the physical domain
      ! thus threated at internal blocks with no physical boundary
      if (periodB(2)) ign2=1+modulo(ign2-1,ng2(level))
      if (ign2 > ng2(level)) then
         if(phi_ > 0 .and. poleB(2,2)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*2)=.false.
         else
           ps(igrid)%is_physical_boundary(2*2)=.true.
         end if
      else if (ign2 < 1) then
         if(phi_ > 0 .and. poleB(1,2)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*2-1)=.false.
         else
           ps(igrid)%is_physical_boundary(2*2-1)=.true.
         end if
      end if
    end do
    
    
    do i3=-1,1
      if(i3==0) cycle
      ign3=ig3+i3
      ! blocks at periodic boundary have neighbors in the physical domain
      ! thus threated at internal blocks with no physical boundary
      if (periodB(3)) ign3=1+modulo(ign3-1,ng3(level))
      if (ign3 > ng3(level)) then
         if(phi_ > 0 .and. poleB(2,3)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*3)=.false.
         else
           ps(igrid)%is_physical_boundary(2*3)=.true.
         end if
      else if (ign3 < 1) then
         if(phi_ > 0 .and. poleB(1,3)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*3-1)=.false.
         else
           ps(igrid)%is_physical_boundary(2*3-1)=.true.
         end if
      end if
    end do
    
    if(any(ps(igrid)%is_physical_boundary)) then
      phyboundblock(igrid)=.true.
    else
      phyboundblock(igrid)=.false.
    end if
  end subroutine alloc_node
  
  !> allocate memory to physical state of igrid node
  subroutine alloc_state(igrid, s, ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
     ixGmax3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
     ixGextmax3, alloc_once_for_ps)
    use mod_global_parameters
    type(state) :: s
    integer, intent(in) :: igrid, ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3
    logical, intent(in) :: alloc_once_for_ps
    integer             :: ixGsmin1,ixGsmin2,ixGsmin3,ixGsmax1,ixGsmax2,&
       ixGsmax3
    !opedit: debug:
    integer             :: idbg
  
    !allocate(s%w(ixG^S,1:nw))
    s%w => bg(s%istep)%w(:,:,:,:,igrid)
    s%igrid=igrid
    s%w=0.d0
    s%ixGmin1=ixGmin1;s%ixGmin2=ixGmin2;s%ixGmin3=ixGmin3;s%ixGmax1=ixGmax1
    s%ixGmax2=ixGmax2;s%ixGmax3=ixGmax3;
     ixGsmin1 = ixGmin1-1; ixGsmax1 = ixGmax1; ixGsmin2 = ixGmin2-1
     ixGsmax2 = ixGmax2; ixGsmin3 = ixGmin3-1; ixGsmax3 = ixGmax3
    if(stagger_grid) then
      allocate(s%ws(ixGsmin1:ixGsmax1,ixGsmin2:ixGsmax2,ixGsmin3:ixGsmax3,&
         1:nws))
      s%ws=0.d0
      if(record_electric_field) then
        allocate(s%we(ixGsmin1:ixGsmax1,ixGsmin2:ixGsmax2,ixGsmin3:ixGsmax3,&
           1:nws))
        s%we=0.d0
      end if
      s%ixGsmin1=ixGsmin1;s%ixGsmin2=ixGsmin2;s%ixGsmin3=ixGsmin3
      s%ixGsmax1=ixGsmax1;s%ixGsmax2=ixGsmax2;s%ixGsmax3=ixGsmax3;
    end if
    if(alloc_once_for_ps) then
      ! allocate extra variables for ps state
      if(nw_extra>0) allocate(s%wextra(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1:nw_extra))
      ! allocate coordinates
      allocate(s%x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim))
      allocate(s%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
         ixGextmin3:ixGextmax3,1:ndim), s%ds(ixGextmin1:ixGextmax1,&
         ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1:ndim),&
         s%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
         ixGextmin3:ixGextmax3,1:3))
      allocate(s%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
         ixGextmin3:ixGextmax3))
      allocate(s%surfaceC(ixGsmin1:ixGsmax1,ixGsmin2:ixGsmax2,&
         ixGsmin3:ixGsmax3,1:ndim), s%surface(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1:ndim))
      ! allocate physical boundary flag
      allocate(s%is_physical_boundary(2*ndim))
      if(local_timestep) then
        allocate(s%dt(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))
      endif
  
      if(B0field) then
        allocate(s%B0(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndir,&
           0:ndim))
        allocate(s%J0(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
           7-2*ndir:3))
      end if
      if(number_equi_vars > 0) then
        allocate(s%equi_vars(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
           1:number_equi_vars,0:ndim))
      endif
  
      ! allocate space for special values for each block state
      if(phys_trac) then
        ! special_values(1) Tcoff local
        ! special_values(2) Tmax local
        ! special_values(3:2+ndim) Bdir local
        allocate(s%special_values(ndim+2))
      end if
    else
      ! share common info from ps states to save memory
      if(nw_extra>0) s%wextra=>ps(igrid)%wextra
      s%x=>ps(igrid)%x
      s%dx=>ps(igrid)%dx
      s%ds=>ps(igrid)%ds
      s%dsC=>ps(igrid)%dsC
      s%dvolume=>ps(igrid)%dvolume
      s%surfaceC=>ps(igrid)%surfaceC
      s%surface=>ps(igrid)%surface
      s%is_physical_boundary=>ps(igrid)%is_physical_boundary
      if(B0field) then
        s%B0=>ps(igrid)%B0
        s%J0=>ps(igrid)%J0
      end if
      if(number_equi_vars > 0) then
        s%equi_vars=>ps(igrid)%equi_vars
      endif
      if(phys_trac) s%special_values=>ps(igrid)%special_values
    end if
  end subroutine alloc_state
  
  !> allocate memory to one-level coarser physical state of igrid node
  subroutine alloc_state_coarse(igrid, s, ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
     ixGmax2,ixGmax3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
     ixGextmax3)
    use mod_global_parameters
    type(state) :: s
    integer, intent(in) :: igrid, ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3
    integer             :: ixGsmin1,ixGsmin2,ixGsmin3,ixGsmax1,ixGsmax2,&
       ixGsmax3
  
    allocate(s%w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw))
    s%igrid=igrid
    s%w=0.d0
    s%ixGmin1=ixGmin1;s%ixGmin2=ixGmin2;s%ixGmin3=ixGmin3;s%ixGmax1=ixGmax1
    s%ixGmax2=ixGmax2;s%ixGmax3=ixGmax3;
     ixGsmin1 = ixGmin1-1; ixGsmax1 = ixGmax1; ixGsmin2 = ixGmin2-1
     ixGsmax2 = ixGmax2; ixGsmin3 = ixGmin3-1; ixGsmax3 = ixGmax3
    if(stagger_grid) then
      allocate(s%ws(ixGsmin1:ixGsmax1,ixGsmin2:ixGsmax2,ixGsmin3:ixGsmax3,&
         1:nws))
      s%ws=0.d0
      s%ixGsmin1=ixGsmin1;s%ixGsmin2=ixGsmin2;s%ixGsmin3=ixGsmin3
      s%ixGsmax1=ixGsmax1;s%ixGsmax2=ixGsmax2;s%ixGsmax3=ixGsmax3;
    end if
    if(B0fieldAllocCoarse) then
      allocate(s%B0(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndir,&
         0:ndim))
    end if
    ! allocate coordinates
    allocate(s%x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim))
    allocate(s%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3,1:ndim), s%ds(ixGextmin1:ixGextmax1,&
       ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,1:ndim),&
       s%dsC(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
       1:3))
    allocate(s%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3))
    allocate(s%surfaceC(ixGsmin1:ixGsmax1,ixGsmin2:ixGsmax2,ixGsmin3:ixGsmax3,&
       1:ndim), s%surface(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       1:ndim))
    ! allocate physical boundary flag
    allocate(s%is_physical_boundary(2*ndim))
  end subroutine alloc_state_coarse
  
  subroutine dealloc_state(igrid, s,dealloc_x)
    use mod_global_parameters
    integer, intent(in) :: igrid
    type(state) :: s
    logical, intent(in) :: dealloc_x
  
    deallocate(s%w)
    if(stagger_grid) then
      deallocate(s%ws)
      if(record_electric_field) deallocate(s%we)
    end if
    if(dealloc_x) then
      if(nw_extra>0) deallocate(s%wextra)
      ! deallocate coordinates
      deallocate(s%x)
      deallocate(s%dx,s%ds,s%dsC)
      deallocate(s%dvolume)
      deallocate(s%surfaceC,s%surface)
      deallocate(s%is_physical_boundary)
      if(B0field) then
        deallocate(s%B0)
        deallocate(s%J0)
      end if
      if(number_equi_vars > 0) then
        deallocate(s%equi_vars)
      end if
    else
      nullify(s%x,s%dx,s%ds,s%dsC,s%dvolume,s%surfaceC,s%surface)
      nullify(s%is_physical_boundary)
      if(B0field) nullify(s%B0,s%J0)
      if(number_equi_vars > 0) then
        nullify(s%equi_vars)
      end if
      if(nw_extra>0) nullify(s%wextra)
    end if
  end subroutine dealloc_state
  
  subroutine dealloc_state_coarse(igrid, s)
    use mod_global_parameters
    integer, intent(in) :: igrid
    type(state) :: s
  
    deallocate(s%w)
    if(stagger_grid) then
      deallocate(s%ws)
    end if
    if(B0fieldAllocCoarse) then
      deallocate(s%B0)
    end if
    ! deallocate coordinates
    deallocate(s%x)
    deallocate(s%dx,s%ds,s%dsC)
    deallocate(s%dvolume)
    deallocate(s%surfaceC,s%surface)
    deallocate(s%is_physical_boundary)
  end subroutine dealloc_state_coarse
  
  subroutine dealloc_node(igrid)
    use mod_global_parameters
  
    integer, intent(in) :: igrid
  
    if (igrid==0) then
       call mpistop("trying to delete a non-existing grid in dealloc_node")
    end if
  
    call dealloc_state(igrid, ps(igrid),.true.)
    call dealloc_state_coarse(igrid, psc(igrid))
    if(.not.convert) then
      call dealloc_state(igrid, ps1(igrid),.false.)
      ! deallocate temporary solution space
      select case (t_integrator)
      case(ssprk3,ssprk4,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
        call dealloc_state(igrid, ps2(igrid),.false.)
      case(RK3_BT,rk4,ssprk5,IMEX_CB3a)
        call dealloc_state(igrid, ps2(igrid),.false.)
        call dealloc_state(igrid, ps3(igrid),.false.)
      case(IMEX_ARS3,IMEX_232)
        call dealloc_state(igrid, ps2(igrid),.false.)
        call dealloc_state(igrid, ps3(igrid),.false.)
        call dealloc_state(igrid, ps4(igrid),.false.)
      end select
    end if
  
  end subroutine dealloc_node

end module mod_amr_solution_node
