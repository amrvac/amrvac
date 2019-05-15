!> coarsen sibling blocks into one block
subroutine coarsen_grid_siblings(igrid,ipe,child_igrid,child_ipe,active)
  use mod_global_parameters

  integer, intent(in) :: igrid, ipe
  integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
  logical, intent(in) :: active

  integer :: igridFi, ipeFi, ixCo^L, ixCoG^L, ixCoM^L, ic^D, idir

  if (ipe==mype) call alloc_node(igrid)

  ! New passive cell, coarsen from initial condition:
  if (.not. active) then
     if (ipe == mype) then
        call initial_condition(igrid)
        {do ic^DB=1,2\}
        igridFi=child_igrid(ic^D)
        ipeFi=child_ipe(ic^D)
        !if (ipeFi==mype) then
        !   ! remove solution space of child      
        !   call dealloc_node(igridFi)
        !end if
        {end do\}
     end if
     return
  end if

  {do ic^DB=1,2\}
     igridFi=child_igrid(ic^D)
     ipeFi=child_ipe(ic^D)

     if (ipeFi==mype) then
        ^D&dxlevel(^D)=rnode(rpdx^D_,igridFi);
        if (ipe==mype) then
           ixComin^D=ixMlo^D+(ic^D-1)*(ixMhi^D-ixMlo^D+1)/2;
           ixComax^D=ixMhi^D+(ic^D-2)*(ixMhi^D-ixMlo^D+1)/2;

           call coarsen_grid(ps(igridFi),ixG^LL,ixM^LL,ps(igrid),ixG^LL,ixCo^L)
           ! remove solution space of child
           !call dealloc_node(igridFi)
        else
           ixCoGmin^D=1;
           ixCoGmax^D=ixGhi^D/2+nghostcells;
           ixCoM^L=ixCoG^L^LSUBnghostcells;
           call coarsen_grid(ps(igridFi),ixG^LL,ixM^LL,psc(igridFi), &
                             ixCoG^L,ixCoM^L)

           itag=ipeFi*max_blocks+igridFi
           isend=isend+1
           call MPI_ISEND(psc(igridFi)%w,1,type_coarse_block,ipe,itag, &
                          icomm,sendrequest(isend),ierrmpi)
           if(stagger_grid) then
             do idir=1,ndim
               itag_stg=(npe+ipeFi+1)*max_blocks+igridFi*(ndir-1+idir)
               call MPI_ISEND(psc(igridFi)%ws,1,type_coarse_block_stg(idir,ic^D),ipe,itag_stg, &
                            icomm,sendrequest_stg(isend),ierrmpi)
             end do
           end if
        end if
     else
        if (ipe==mype) then
           itag=ipeFi*max_blocks+igridFi
           irecv=irecv+1
           call MPI_IRECV(ps(igrid)%w,1,type_sub_block(ic^D),ipeFi,itag, &
                          icomm,recvrequest(irecv),ierrmpi)
           if(stagger_grid) then
             do idir=1,ndim
               itag_stg=(npe+ipeFi+1)*max_blocks+igridFi*(ndir-1+idir)
               call MPI_IRECV(ps(igrid)%ws,1,type_sub_block_stg(idir,ic^D),ipeFi,itag_stg, &
                              icomm,recvrequest_stg(irecv),ierrmpi)
             end do
           end if
        end if
     end if
  {end do\}

end subroutine coarsen_grid_siblings

!> coarsen one grid to its coarser representative
subroutine coarsen_grid(sFi,ixFiG^L,ixFi^L,sCo,ixCoG^L,ixCo^L)
  use mod_global_parameters
  use mod_physics, only: phys_to_primitive, phys_to_conserved
  use mod_constrained_transport, only: faces2centers

  type(state), intent(inout)      :: sFi, sCo
  integer, intent(in) :: ixFiG^L, ixFi^L, ixCoG^L, ixCo^L

  integer :: ixCo^D, ixFi^D, iw
  double precision :: CoFiratio

  associate(wFi=>sFi%w(ixFiG^S,1:nw), wCo=>sCo%w(ixCoG^S,1:nw))
  staggered: associate(wFis=>sFi%ws,wCos=>sCo%ws)
  ! coarsen by 2 in every direction - conservatively

  if(coarsenprimitive) call phys_to_primitive(ixFiG^L,ixFi^L,wFi,sFi%x)

  if(slab_uniform) then
    CoFiratio=one/dble(2**ndim)
    do iw=1,nw
       {do ixCo^DB = ixCo^LIM^DB
          ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
          wCo(ixCo^D,iw)=sum(wFi(ixFi^D:ixFi^D+1,iw))*CoFiratio
       {end do\}
    end do
  else
    do iw=1,nw
      {do ixCo^DB = ixCo^LIM^DB
         ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
         wCo(ixCo^D,iw)= &
             sum(sFi%dvolume(ixFi^D:ixFi^D+1)*wFi(ixFi^D:ixFi^D+1,iw)) &
            /sCo%dvolume(ixCo^D)
      {end do\}
    end do
  end if

  if(stagger_grid) then
    do iw=1,nws
      ! Start one layer before
      {do ixCo^DB = ixComin^DB-kr(^DB,iw),ixComax^DB
         ixFi^DB=2*(ixCo^DB-ixComin^DB+kr(^DB,iw))+ixFimin^DB-kr(^DB,iw)\}
         ! This if statement catches the axis where surface is zero:
         if (sCo%surfaceC(ixCo^D,iw)>1.0d-9*sCo%dvolume(ixCo^D)) then ! Normal case
           wCos(ixCo^D,iw)=sum(sFi%surfaceC(ixFi^D:ixFi^D+1-kr(iw,^D),iw)*wFis(ixFi^D:ixFi^D+1-kr(iw,^D),iw)) &
                /sCo%surfaceC(ixCo^D,iw)
         else ! On axis
           wCos(ixCo^D,iw)=zero
         end if
      {end do\}
    end do
    ! average to fill cell-centred values
    call faces2centers(ixCo^L,sCo)
  end if

  if(coarsenprimitive) then
    call phys_to_conserved(ixFiG^L,ixFi^L,wFi,sFi%x)
    call phys_to_conserved(ixCoG^L,ixCo^L,wCo,sCo%x)
  end if
  end associate staggered
  end associate
end subroutine coarsen_grid
