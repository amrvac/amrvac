!> This module handles the initialization of various components of amrvac
module mod_initialize
  use mod_comm_lib, only: mpistop

  implicit none
  private

  logical :: initialized_already = .false.

  ! Public methods
  public :: initialize_amrvac

contains

  !> Initialize amrvac: read par files and initialize variables
  subroutine initialize_amrvac()
    use mod_global_parameters
    use mod_input_output
    use mod_physics, only: phys_check, phys_check_params
    use mod_usr_methods, only: usr_set_parameters
    use mod_bc_data, only: bc_data_init
    use mod_init_datafromfile, only: read_data_init
    use mod_comm_lib, only: init_comm_types
    use mod_connectivity, only: nbprocs_info

    if (initialized_already) return

    ! add auxiliary variable(s) to update boundary ghost cells
    nwgc=nwgc+nwaux

    ! Check whether the user has loaded a physics module
    call phys_check()

    ! Read input files
    call read_par_files()
    call initialize_vars()
    call init_comm_types()
    call nbprocs_info%init(npe=npe, nigrids=4)
    !$acc update device(nbprocs_info)

    ! Possibly load boundary condition data or initial data
    call bc_data_init()
    call read_data_init()

    if(associated(usr_set_parameters)) call usr_set_parameters()

    call phys_check_params()

    initialized_already = .true.
  end subroutine initialize_amrvac

  !> Initialize (and allocate) simulation and grid variables
  subroutine initialize_vars
    use mod_forest
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve, only: pflux
    use mod_amr_fct, only: pface, fine_neighbors, old_neighbor
    use mod_geometry

    integer :: igrid, level, ipe, ig1,ig2,ig3
    logical :: ok

    allocate(ps(max_blocks))
    allocate(ps1(max_blocks))
    allocate(ps2(max_blocks))
    allocate(ps3(max_blocks))
    allocate(ps4(max_blocks))
    
    allocate(psc(max_blocks))
    !$acc enter data copyin(ps,ps1,ps2,ps3,ps4,psc)
    
    allocate(ps_sub(max_blocks))
    allocate(neighbor(2,-1:1,-1:1,-1:1,max_blocks),neighbor_child(2,0:3,0:3,&
       0:3,max_blocks))
    allocate(neighbor_type(-1:1,-1:1,-1:1,max_blocks),neighbor_active(-1:1,&
       -1:1,-1:1,max_blocks))
    allocate(neighbor_pole(-1:1,-1:1,-1:1,max_blocks))
    allocate(igrids(max_blocks),igrids_active(max_blocks),&
       igrids_passive(max_blocks), idphyb(ndim,max_blocks) )
    allocate(rnode(rnodehi,max_blocks),rnode_sub(rnodehi,max_blocks))
    allocate(node(nodehi,max_blocks),node_sub(nodehi,max_blocks),&
       phyboundblock(max_blocks))
    allocate(pflux(2,3,max_blocks))

    allocate( bg(1:nstep) )
    !$acc enter data copyin(bg)
    do istep = 1 , nstep
       bg(istep)%istep = istep
       allocate( bg(istep)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3, 1:nw,&
            1:max_blocks) )
       !$acc update device(bg(istep))
       !$acc enter data copyin( bg(istep)%w )
    end do

    do igrid = 1, max_blocks
       ps(igrid)%igrid  = igrid; ps(igrid)%istep  = 1
       ps1(igrid)%igrid = igrid; ps1(igrid)%istep = 2
       ps2(igrid)%igrid = igrid; ps2(igrid)%istep = 3
       ps3(igrid)%igrid = igrid; ps3(igrid)%istep = 4
       ps4(igrid)%igrid = igrid; ps4(igrid)%istep = 5
    end do

    
    ! allocate mesh for particles
    if(use_particles) allocate(gridvars(max_blocks))
    if(stagger_grid) then
      allocate(pface(2,3,max_blocks),fine_neighbors(2,2,2,3,max_blocks))
      allocate(old_neighbor(2,-1:11,-1:12,-1:13,max_blocks))
    end if

    it=it_init
    global_time=time_init

    dt=zero

    ! no poles initially
    neighbor_pole=0

    ! check resolution
    if (mod(ixGhi1,2)/=0.or.mod(ixGhi2,2)/=0.or.mod(ixGhi3,2)/=0) then
       call mpistop("mesh widths must give even number grid points")
    end if
    ixMlo1=ixGlo1+nghostcells;ixMlo2=ixGlo2+nghostcells
    ixMlo3=ixGlo3+nghostcells;ixMhi1=ixGhi1-nghostcells
    ixMhi2=ixGhi2-nghostcells;ixMhi3=ixGhi3-nghostcells;
    !$acc update device(ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3)

    if (nbufferx1>(ixMhi1-ixMlo1+1).or.nbufferx2>(ixMhi2-ixMlo2+&
       1).or.nbufferx3>(ixMhi3-ixMlo3+1)) then
       write(unitterm,*) "nbufferx^D bigger than mesh size makes no sense."
       write(unitterm,*) "Decrease nbufferx or increase mesh size"
       call mpistop("")
    end if

    ! initialize dx arrays on finer (>1) levels
    do level=2,refine_max_level
       dx(1,level) = dx(1,level-1) * half
       dx(2,level) = dx(2,level-1) * half
       dx(3,level) = dx(3,level-1) * half  ! refine ratio 2
    end do
    
    ! domain decomposition
    ! physical extent of a grid block at level 1, per dimension
    dg1(1)=dx(1,1)*dble(block_nx1)
    dg2(1)=dx(2,1)*dble(block_nx2)
    dg3(1)=dx(3,1)*dble(block_nx3)
    ! number of grid blocks at level 1 in simulation domain, per dimension
    ng1(1)=nint((xprobmax1-xprobmin1)/dg1(1))
    ng2(1)=nint((xprobmax2-xprobmin2)/dg2(1))
    ng3(1)=nint((xprobmax3-xprobmin3)/dg3(1))
    ! total number of grid blocks at level 1
    nglev1=ng1(1)*ng2(1)*ng3(1)

    do level=2,refine_max_level
       dg1(level)=half*dg1(level-1);dg2(level)=half*dg2(level-1)
       dg3(level)=half*dg3(level-1);
       ng1(level)=ng1(level-1)*2;ng2(level)=ng2(level-1)*2
       ng3(level)=ng3(level-1)*2;
    end do

    ! check that specified stepsize correctly divides domain
    ok=((abs(dble(ng1(1))*dg1(1)-(xprobmax1-xprobmin1))<=&
       smalldouble).and.(abs(dble(ng2(1))*dg2(1)-(xprobmax2-xprobmin2))<=&
       smalldouble).and.(abs(dble(ng3(1))*dg3(1)-(xprobmax3-xprobmin3))<=&
       smalldouble))
    if (.not.ok) then
       write(unitterm,*)"domain cannot be divided by meshes of given gridsize"
       call mpistop("domain cannot be divided by meshes of given gridsize")
    end if

    poleB=.false.
    if (.not.slab) call set_pole

    ! number of grid blocks at level 1 along a dimension, which does not have a pole or periodic boundary, 
    ! must be larger than 1 for a rectangular AMR mesh
    if((ng1(1)/=1.or.ng2(1)/=1.or.ng3(1)/=1).and.refine_max_level>1) then
      
      if(ng1(1)==1.and..not.poleB(1,1).and..not.poleB(2,&
         1).and..not.periodB(1).and..not.aperiodB(1)) then
        write(unitterm,"(a,i2,a)")&
            "number of grid blocks at level 1 in dimension",1,&
           " be larger than 1 for a rectangular AMR mesh!"
        write(unitterm,"(a,i1)") "increase domain_nx",1
        call mpistop("")
      end if
      
      
      if(ng2(1)==1.and..not.poleB(1,2).and..not.poleB(2,&
         2).and..not.periodB(2).and..not.aperiodB(2)) then
        write(unitterm,"(a,i2,a)")&
            "number of grid blocks at level 1 in dimension",2,&
           " be larger than 1 for a rectangular AMR mesh!"
        write(unitterm,"(a,i1)") "increase domain_nx",2
        call mpistop("")
      end if
      
      
      if(ng3(1)==1.and..not.poleB(1,3).and..not.poleB(2,&
         3).and..not.periodB(3).and..not.aperiodB(3)) then
        write(unitterm,"(a,i2,a)")&
            "number of grid blocks at level 1 in dimension",3,&
           " be larger than 1 for a rectangular AMR mesh!"
        write(unitterm,"(a,i1)") "increase domain_nx",3
        call mpistop("")
      end if
      
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
    if (nbufferx1/=0.or.nbufferx2/=0.or.nbufferx3/=0) then
       allocate(buffer(max_blocks,0:npe-1))
       buffer=.false.
    end if
    allocate(igrid_inuse(max_blocks,0:npe-1))
    igrid_inuse=.false.

    allocate(tree_root(1:ng1(1),1:ng2(1),1:ng3(1)))
    do ig3=1,ng3(1)
    do ig2=1,ng2(1)
    do ig1=1,ng1(1)
    nullify(tree_root(ig1,ig2,ig3)%node)
    end do
    end do
    end do

    ! define index ranges and MPI send/receive derived datatype for ghost-cell swap
    call init_bc()
    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f
    call create_bc_mpi_datatype(iwstart,nwgc)

  end subroutine initialize_vars


end module mod_initialize
