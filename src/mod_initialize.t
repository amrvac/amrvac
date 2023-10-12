!> This module handles the initialization of various components of gmunu
module mod_initialize

  implicit none
  private

  logical :: initialized_already = .false.

  ! Public methods
  public :: initialize_gmunu
  public :: initlevelone
  public :: initial_condition
  public :: modify_IC
  public :: improve_initial_condition

contains

  !> Initialize gmunu: read par files and initialize variables
  subroutine initialize_gmunu()
    use mod_input_output
    use mod_physics,     only: phys_check, phys_check_params
    use mod_usr_methods, only: usr_set_parameters
    use mod_bc_data,     only: bc_data_init

    if (initialized_already) return

    ! Check whether the user has loaded a physics module
    call phys_check()

    ! Read input files
    call read_par_files()

    call initialize_vars()

    ! set physical boundary condition
    call bc_init()
    ! Possibly load boundary condition data
    call bc_data_init()

    if(associated(usr_set_parameters)) call usr_set_parameters()

    call phys_check_params()

    initialized_already = .true.
  end subroutine initialize_gmunu

  !> Initialize (and allocate) simulation and grid variables
  subroutine initialize_vars
    use mod_forest
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve, only: pflux
    use mod_amr_fct, only: pface, fine_neighbors, old_neighbor
    use mod_geometry

    integer :: igrid, level, ipe, ig^D
    logical :: ok

    ! fixme: maybe we dont need to use all ps?
    allocate(mesh(max_blocks))
    allocate(mesh_co(max_blocks))

    allocate(metric(max_blocks))
    allocate(metric_co(max_blocks))

    allocate(ps(max_blocks))
    allocate(ps1(max_blocks))
    allocate(ps2(max_blocks))
    allocate(ps3(max_blocks))
    allocate(ps4(max_blocks))
    allocate(pso(max_blocks))
    allocate(psc(max_blocks))
    allocate(ps_sub(max_blocks))
    allocate(neighbor(2,-1:1^D&,max_blocks),neighbor_child(2,0:3^D&,max_blocks))
    allocate(neighbor_type(-1:1^D&,max_blocks),neighbor_active(-1:1^D&,max_blocks))
    allocate(neighbor_pole(-1:1^D&,max_blocks))
    allocate(igrids(max_blocks),igrids_active(max_blocks),igrids_passive(max_blocks))
    allocate(rnode(rnodehi,max_blocks),rnode_sub(rnodehi,max_blocks))
    allocate(node(nodehi,max_blocks),node_sub(nodehi,max_blocks),phyboundblock(max_blocks))
    allocate(pflux(2,^ND,max_blocks))
    if(stagger_grid) then
      allocate(pface(2,^ND,max_blocks),fine_neighbors(2^D&,^ND,max_blocks))
      allocate(old_neighbor(2,-1:1^D,max_blocks))
    end if

    it          = it_init
    global_time = time_init

    dt = zero

    ! no poles initially
    neighbor_pole=0

    ! check resolution
    if ({mod(ixGhi^D,2)/=0|.or.}) then
       call mpistop("mesh widths must give even number grid points")
    end if
    ixM^LL=ixG^LL^LSUBnghostcells;

    if (nbufferx^D>(ixMhi^D-ixMlo^D+1)|.or.) then
       write(unitterm,*) "nbufferx^D bigger than mesh size makes no sense."
       write(unitterm,*) "Decrease nbufferx or increase mesh size"
       call mpistop("")
    end if

    ! fixme: split read par and initial var

    ! initialize dx arrays on finer (>1) levels
    do level=2,refine_max_level
       {dx(^D,level) = dx(^D,level-1) * half\}  ! refine ratio 2
    end do

    ! fixme: maybe we should assign values here

    ! domain decomposition
    ! physical extent of a grid block at level 1, per dimension
    ^D&dg^D(1)=dx(^D,1)*dble(block_nx^D)\
    ! number of grid blocks at level 1 in simulation domain, per dimension
    ^D&ng^D(1)=nint((xprobmax^D-xprobmin^D)/dg^D(1))\
    ! total number of grid blocks at level 1
    nglev1={ng^D(1)*}

    do level=2,refine_max_level
       dg^D(level)=half*dg^D(level-1);
       ng^D(level)=ng^D(level-1)*2;
    end do

    ! check that specified stepsize correctly divides domain
    ok=({(abs(dble(ng^D(1))*dg^D(1)-(xprobmax^D-xprobmin^D))<=smalldouble)|.and.})
    if (.not.ok) then
       write(unitterm,*)"domain cannot be divided by meshes of given gridsize"
       call mpistop("domain cannot be divided by meshes of given gridsize")
    end if

    poleB=.false.
    if (.not.slab) call set_pole

    ! number of grid blocks at level 1 along a dimension, which does not have a pole or periodic boundary, 
    ! must be larger than 1 for a rectangular AMR mesh
    if(({ng^D(1)/=1|.or.}).and.refine_max_level>1) then
      {
      if(ng^D(1)==1.and..not.poleB(1,^D).and.&
         .not.poleB(2,^D).and..not.periodB(^D).and..not.aperiodB(^D)) then
        write(unitterm,"(a,i2,a)") "number of grid blocks at level 1 in dimension",^D,&
                          " be larger than 1 for a rectangular AMR mesh!"
        write(unitterm,"(a,i1)") "increase domain_nx",^D
        call mpistop("")
      end if
      \}
    end if

    ! initialize connectivity data
    igridstail=0

    ! allocate memory for forest data structures
    allocate(level_head(refine_max_level),level_tail(refine_max_level))
    do level=1,refine_max_level
       nullify(level_head(level)%node,level_tail(level)%node)
    end do

    allocate(igrid_to_node(max_blocks,0:npe-1))
    do ipe=0,npe-1
       do igrid=1,max_blocks
          nullify(igrid_to_node(igrid,ipe)%node)
       end do
    end do

    allocate(sfc(1:3,max_blocks*npe))

    allocate(igrid_to_sfc(max_blocks))

    sfc=0
    allocate(Morton_start(0:npe-1),Morton_stop(0:npe-1))
    allocate(Morton_sub_start(0:npe-1),Morton_sub_stop(0:npe-1))

    allocate(nleafs_level(1:nlevelshi))

    allocate(coarsen(max_blocks,0:npe-1),refine(max_blocks,0:npe-1))
    coarsen=.false.
    refine=.false.
    if (nbufferx^D/=0|.or.) then
       allocate(buffer(max_blocks,0:npe-1))
       buffer=.false.
    end if
    allocate(igrid_inuse(max_blocks,0:npe-1))
    igrid_inuse=.false.

    allocate(tree_root(1:ng^D(1)))
    {do ig^DB=1,ng^DB(1)\}
       nullify(tree_root(ig^D)%node)
    {end do\}

    ! define index ranges and MPI send/receive derived datatype for ghost-cell swap
    call init_bc()
    call init_comm_types()

  end subroutine initialize_vars

  !> Generate and initialize all grids at the coarsest level (level one)
  subroutine initlevelone
    use mod_global_parameters
    use mod_ghostcells_update
  
    integer :: iigrid, igrid{#IFDEF EVOLVINGBOUNDARY , Morton_no}
    integer :: ierrmpi
  
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
    call getbc(global_time,0.d0,ps,1,nmetric,gc_metric)
    call getbc(global_time,0.d0,ps,1,nw,gc_hydro)
  
  end subroutine initlevelone
  
  !> fill in initial condition
  subroutine initial_condition(igrid)
    ! Need only to set the mesh values (can leave ghost cells untouched)
    use mod_usr_methods, only: usr_init_one_grid
    use mod_global_parameters
  
    integer, intent(in) :: igrid
  
    ps(igrid)%is_prim=.False.
    ps(igrid)%w(ixG^T,1:nw)=zero
  
    ! in case gradient routine used in initial condition, ensure geometry known
    block=>ps(igrid)
    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  
    if (.not. associated(usr_init_one_grid)) then
       call mpistop("usr_init_one_grid not defined")
    else
       call usr_init_one_grid(ixG^LL,ixM^LL,ps(igrid))
    end if

    if (ps(igrid)%is_prim) &
       call mpistop("usr_init_one_grid: w has to be conserved variables")
  
  end subroutine initial_condition
  
  !> modify initial condition
  subroutine modify_IC
    use mod_usr_methods, only: usr_init_one_grid
    use mod_global_parameters
  
    integer :: iigrid, igrid
  
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       block=>ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  
       if (.not. associated(usr_init_one_grid)) then
          call mpistop("usr_init_one_grid not defined")
       else
          call usr_init_one_grid(ixG^LL,ixM^LL,ps(igrid))
       end if
    end do
  
  end subroutine modify_IC
  
  !> improve initial condition after initialization
  subroutine improve_initial_condition()
    use mod_global_parameters
    use mod_usr_methods
    use mod_constrained_transport
    use mod_ghostcells_update
  
    logical :: active
  
    if(associated(usr_improve_initial_condition)) then
      if (mype==0) then
         print*,'-------------------------------------------------------------------------------'
         write(*,'(a,f17.3,a)')' Grid tree is set, now improve the initial condition'
         print*,'-------------------------------------------------------------------------------'
      end if
      call usr_improve_initial_condition
    end if
  
    ! update bc after improving initial condition
    call getbc(global_time,0.d0,ps,1,nmetric,gc_metric)
    call getbc(global_time,0.d0,ps,1,nw,gc_hydro)
  
  end subroutine improve_initial_condition

end module mod_initialize
