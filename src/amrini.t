!> Generate and initialize all grids at the coarsest level (level one)
subroutine initlevelone
  use mod_global_parameters
  use mod_ghostcells_update

  integer :: iigrid, igrid{#IFDEF EVOLVINGBOUNDARY , Morton_no}

  levmin=1
  levmax=1

  call init_forest_root

  call getigrids
  call build_connectivity

  ! fill solution space of all root grids
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     call alloc_node(igrid)
     ! in case gradient routine used in initial condition, ensure geometry known
     call initial_condition(igrid)
  end do
  {#IFDEF EVOLVINGBOUNDARY
  ! mark physical-boundary blocks on space-filling curve
  do Morton_no=Morton_start(mype),Morton_stop(mype)
     igrid=sfc_to_igrid(Morton_no)
     if (phyboundblock(igrid)) sfc_phybound(Morton_no)=1
  end do
  call MPI_ALLREDUCE(MPI_IN_PLACE,sfc_phybound,nleafs,MPI_INTEGER,&
                     MPI_SUM,icomm,ierrmpi)
  }

  ! update ghost cells
  call getbc(global_time,0.d0,ps,0,nwflux+nwaux)

end subroutine initlevelone

!> fill in initial condition
subroutine initial_condition(igrid)
  ! Need only to set the mesh values (can leave ghost cells untouched)
  use mod_usr_methods, only: usr_init_one_grid
  use mod_global_parameters

  integer, intent(in) :: igrid

  ps(igrid)%w(ixG^T,1:nw)=zero

  saveigrid=igrid
  ! in case gradient routine used in initial condition, ensure geometry known
  block=>ps(igrid)
  ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  typelimiter=type_limiter(node(plevel_,igrid))
  typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

  if (.not. associated(usr_init_one_grid)) then
     call mpistop("usr_init_one_grid not defined")
  else
     call usr_init_one_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x)
  end if

end subroutine initial_condition

!> modify initial condition
subroutine modify_IC
  use mod_usr_methods, only: usr_init_one_grid
  use mod_global_parameters

  integer :: iigrid, igrid

  do iigrid=1,igridstail; igrid=igrids(iigrid);
     saveigrid=igrid
     block=>ps(igrid)
     ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
     typelimiter=type_limiter(node(plevel_,igrid))
     typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

     if (.not. associated(usr_init_one_grid)) then
        call mpistop("usr_init_one_grid not defined")
     else
        call usr_init_one_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x)
     end if
  end do

end subroutine modify_IC
