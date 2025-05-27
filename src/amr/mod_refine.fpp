module mod_refine

  implicit none
  private
 
  public :: refine_grids
 
contains


  !> refine one block to its children blocks
  subroutine refine_grids(child_igrid,child_ipe,igrid,ipe,active)
    use mod_global_parameters
    use mod_initialize_amr, only: initial_condition
    use mod_amr_solution_node, only: alloc_node, dealloc_node
  
    integer, dimension(2,2,2), intent(in) :: child_igrid, child_ipe
    integer, intent(in) :: igrid, ipe
    logical, intent(in) :: active
  
    integer :: ic1,ic2,ic3
  
    ! allocate solution space for new children
    do ic3=1,2
    do ic2=1,2
    do ic1=1,2
       call alloc_node(child_igrid(ic1,ic2,ic3))
    end do
    end do
    end do
  
    if ((time_advance .and. active).or.convert.or.reset_grid) then
       ! prolong igrid to new children
       call prolong_grid(child_igrid,child_ipe,igrid,ipe)
    else
       ! Fill new created children with initial condition
       do ic3=1,2
       do ic2=1,2
       do ic1=1,2
          call initial_condition(child_igrid(ic1,ic2,ic3))
       end do
       end do
       end do
    end if
  
    ! remove solution space of igrid to save memory when converting data
    if(convert) call dealloc_node(igrid)
  end subroutine refine_grids
  
  !> prolong one block
  subroutine prolong_grid(child_igrid,child_ipe,igrid,ipe)
    use mod_physics, only: phys_to_primitive, phys_to_conserved
    use mod_global_parameters
    use mod_amr_fct, only: old_neighbors
  
    integer, dimension(2,2,2), intent(in) :: child_igrid, child_ipe
    integer, intent(in) :: igrid, ipe
  
    integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, ichild, ixComin1,&
       ixComin2,ixComin3,ixComax1,ixComax2,ixComax3, ic1,ic2,ic3
    double precision :: dxCo1,dxCo2,dxCo3, xComin1,xComin2,xComin3, dxFi1,&
       dxFi2,dxFi3, xFimin1,xFimin2,xFimin3
  
    ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmin3=ixMlo3-1;ixmax1=ixMhi1+1
    ixmax2=ixMhi2+1;ixmax3=ixMhi3+1;
  
    if(prolongprimitive) call phys_to_primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
       ixGhi2,ixGhi3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,ps(igrid)%w,&
       ps(igrid)%x)
  
    xComin1=rnode(rpxmin1_,igrid)
    xComin2=rnode(rpxmin2_,igrid)
    xComin3=rnode(rpxmin3_,igrid)
    dxCo1=rnode(rpdx1_,igrid)
    dxCo2=rnode(rpdx2_,igrid)
    dxCo3=rnode(rpdx3_,igrid)
  
    if(stagger_grid) call old_neighbors(child_igrid,child_ipe,igrid,ipe)
  
    do ic3=1,2
    do ic2=1,2
    do ic1=1,2
      ichild=child_igrid(ic1,ic2,ic3)
  
      ixComin1=ixMlo1+(ic1-1)*block_nx1/2
      ixComin2=ixMlo2+(ic2-1)*block_nx2/2
      ixComin3=ixMlo3+(ic3-1)*block_nx3/2
      ixComax1=ixMhi1+(ic1-2)*block_nx1/2
      ixComax2=ixMhi2+(ic2-2)*block_nx2/2
      ixComax3=ixMhi3+(ic3-2)*block_nx3/2
  
      xFimin1=rnode(rpxmin1_,ichild)
      xFimin2=rnode(rpxmin2_,ichild)
      xFimin3=rnode(rpxmin3_,ichild)
      dxFi1=rnode(rpdx1_,ichild)
      dxFi2=rnode(rpdx2_,ichild)
      dxFi3=rnode(rpdx3_,ichild)
      call prolong_2nd(ps(igrid),ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
         ixComax3,ps(ichild), dxCo1,dxCo2,dxCo3,xComin1,xComin2,xComin3,dxFi1,&
         dxFi2,dxFi3,xFimin1,xFimin2,xFimin3,igrid,ichild)
      !call prolong_1st(ps(igrid)%w,ixCo^L,ps(ichild)%w,ps(ichild)%x)
    end do
    end do
    end do
  
    if (prolongprimitive) call phys_to_conserved(ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
       ixGhi2,ixGhi3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,ps(igrid)%w,&
       ps(igrid)%x)
  
  end subroutine prolong_grid
  
  !> do 2nd order prolongation
  subroutine prolong_2nd(sCo,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
     ixComax3,sFi,dxCo1,dxCo2,dxCo3,xComin1,xComin2,xComin3,dxFi1,dxFi2,dxFi3,&
     xFimin1,xFimin2,xFimin3,igridCo,igridFi)
    use mod_physics, only: phys_to_conserved, phys_handle_small_values
    use mod_global_parameters
    use mod_amr_fct, only: already_fine, prolong_2nd_stg
  
    integer, intent(in) :: ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
       ixComax3, igridFi, igridCo
    double precision, intent(in) :: dxCo1,dxCo2,dxCo3, xComin1,xComin2,xComin3,&
        dxFi1,dxFi2,dxFi3, xFimin1,xFimin2,xFimin3
    type(state), intent(in)      :: sCo
    type(state), intent(inout)   :: sFi
  
    integer :: ixCo1,ixCo2,ixCo3, jxCo1,jxCo2,jxCo3, hxCo1,hxCo2,hxCo3, ixFi1,&
       ixFi2,ixFi3, ix1,ix2,ix3, idim, iw, ixCgmin1,ixCgmin2,ixCgmin3,ixCgmax1,&
       ixCgmax2,ixCgmax3, el
    double precision :: slopeL, slopeR, slopeC, signC, signR
    double precision :: slope(nw,ndim)
    double precision :: eta1,eta2,eta3
    logical :: fine_min1,fine_min2,fine_min3,fine_max1,fine_max2,fine_max3
  
    associate(wCo=>sCo%w, wFi=>sFi%w)
    ixCgmin1=ixComin1;ixCgmin2=ixComin2;ixCgmin3=ixComin3;ixCgmax1=ixComax1
    ixCgmax2=ixComax2;ixCgmax3=ixComax3;
    do ixCo3 = ixCgmin3,ixCgmax3
       ! lower left grid index in finer child block
       ixFi3=2*(ixCo3-ixComin3)+ixMlo3
    do ixCo2 = ixCgmin2,ixCgmax2
       ! lower left grid index in finer child block
       ixFi2=2*(ixCo2-ixComin2)+ixMlo2
    do ixCo1 = ixCgmin1,ixCgmax1
       ! lower left grid index in finer child block
       ixFi1=2*(ixCo1-ixComin1)+ixMlo1
  
       do idim=1,ndim
          hxCo1=ixCo1-kr(1,idim)
          hxCo2=ixCo2-kr(2,idim)
          hxCo3=ixCo3-kr(3,idim)
          jxCo1=ixCo1+kr(1,idim)
          jxCo2=ixCo2+kr(2,idim)
          jxCo3=ixCo3+kr(3,idim)
  
          do iw=1,nw
             slopeL=wCo(ixCo1,ixCo2,ixCo3,iw)-wCo(hxCo1,hxCo2,hxCo3,iw)
             slopeR=wCo(jxCo1,jxCo2,jxCo3,iw)-wCo(ixCo1,ixCo2,ixCo3,iw)
             slopeC=half*(slopeR+slopeL)
  
             ! get limited slope
             signR=sign(one,slopeR)
             signC=sign(one,slopeC)
             !select case(prolong_limiter)
             !case(1)
             !  ! unlimited
             !  slope(iw,idim)=slopeC
             !case(2)
             !  ! minmod
             !  slope(iw,idim)=signR*max(zero,min(dabs(slopeR), &
             !                                    signR*slopeL))
             !case(3)
             !  ! woodward
             !  slope(iw,idim)=two*signR*max(zero,min(dabs(slopeR), &
             !                     signR*slopeL,signR*half*slopeC))
             !case(4)
             !  ! koren
             !  slope(iw,idim)=signR*max(zero,min(two*signR*slopeL, &
             !   (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
             !case default
               slope(iw,idim)=signC*max(zero,min(dabs(slopeC), signC*slopeL,&
                  signC*slopeR))
             !end select
          end do
       end do
       ! cell-centered coordinates of coarse grid point
       !^D&xCo^D=xCo({ixCo^DD},^D)
       do ix3=ixFi3,ixFi3+1 
       do ix2=ixFi2,ixFi2+1 
       do ix1=ixFi1,ixFi1+1 
          ! cell-centered coordinates of fine grid point
          !^D&xFi^D=xFi({ix^DD},^D)
          if(slab_uniform) then
            ! normalized distance between fine/coarse cell center
            ! in coarse cell: ranges from -0.5 to 0.5 in each direction
            ! (origin is coarse cell center)
            ! hence this is +1/4 or -1/4 on cartesian mesh
            !eta^D=(xFi^D-xCo^D)*invdxCo^D;
            eta1=0.5d0*(dble(ix1-ixFi1)-0.5d0)
            eta2=0.5d0*(dble(ix2-ixFi2)-0.5d0)
            eta3=0.5d0*(dble(ix3-ixFi3)-0.5d0);
          else
            ! forefactor is -0.5d0 when ix=ixFi and +0.5d0 for ixFi+1
            eta1=(dble(ix1-ixFi1)-0.5d0)*(one-sFi%dvolume(ix1,ix2,&
               ix3) /sum(sFi%dvolume(ixFi1:ixFi1+1,ix2,ix3)))  
            ! forefactor is -0.5d0 when ix=ixFi and +0.5d0 for ixFi+1
            eta2=(dble(ix2-ixFi2)-0.5d0)*(one-sFi%dvolume(ix1,ix2,&
               ix3) /sum(sFi%dvolume(ix1,ixFi2:ixFi2+1,ix3)))  
            ! forefactor is -0.5d0 when ix=ixFi and +0.5d0 for ixFi+1
            eta3=(dble(ix3-ixFi3)-0.5d0)*(one-sFi%dvolume(ix1,ix2,&
               ix3) /sum(sFi%dvolume(ix1,ix2,ixFi3:ixFi3+1)))  
          end if
          wFi(ix1,ix2,ix3,1:nw) = wCo(ixCo1,ixCo2,ixCo3,1:nw) + (slope(1:nw,&
             1)*eta1)+(slope(1:nw,2)*eta2)+(slope(1:nw,3)*eta3)
       end do
       end do
       end do
    end do
    end do
    end do
    if(stagger_grid) then
      call already_fine(sFi,igridFi,fine_min1,fine_min2,fine_min3,fine_max1,&
         fine_max2,fine_max3)
      call prolong_2nd_stg(sCo,sFi,ixComin1,ixComin2,ixComin3,ixComax1,&
         ixComax2,ixComax3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,dxCo1,&
         dxCo2,dxCo3,xComin1,xComin2,xComin3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,&
         xFimin3,.false.,fine_min1,fine_min2,fine_min3,fine_max1,fine_max2,&
         fine_max3)
    end if
  
    if(fix_small_values) call phys_handle_small_values(prolongprimitive,wFi,&
       sFi%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,&
       ixMhi1,ixMhi2,ixMhi3,'prolong_2nd')
    if(prolongprimitive) call phys_to_conserved(ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
       ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,wFi,sFi%x)
    end associate
  
  end subroutine prolong_2nd
  
  !> do 1st order prolongation
  subroutine prolong_1st(wCo,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
     ixComax3,wFi,xFi)
    use mod_global_parameters
  
    integer, intent(in) :: ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
       ixComax3
    double precision, intent(in) :: wCo(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nw), xFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       1:ndim)
    double precision, intent(out) :: wFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nw)
  
    integer :: ixCo1,ixCo2,ixCo3, ixFi1,ixFi2,ixFi3, iw
    integer :: ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3
  
    do ixCo3 = ixComin3,ixComax3
       ixFi3=2*(ixCo3-ixComin3)+ixMlo3
    do ixCo2 = ixComin2,ixComax2
       ixFi2=2*(ixCo2-ixComin2)+ixMlo2
    do ixCo1 = ixComin1,ixComax1
       ixFi1=2*(ixCo1-ixComin1)+ixMlo1
       forall(iw=1:nw) wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,ixFi3:ixFi3+1,&
          iw)=wCo(ixCo1,ixCo2,ixCo3,iw)
    end do
    end do
    end do
  
  end subroutine prolong_1st

end module mod_refine
