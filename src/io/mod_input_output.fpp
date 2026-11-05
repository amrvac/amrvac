#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute
!> Module for reading input and writing output
module mod_input_output
#ifdef USE_MPIWRAPPERS
  use mod_mpi_wrapper
#else
#define mpi_file_read_wrapper MPI_FILE_READ
#define mpi_file_read_at_wrapper MPI_FILE_READ_AT
#define mpi_file_write_wrapper MPI_FILE_WRITE
#define mpi_recv_wrapper MPI_RECV
#define mpi_send_wrapper MPI_SEND
#endif
  use mod_comm_lib, only: mpistop

  implicit none
  public


  !> List of compatible versions
  integer, parameter :: compatible_versions(3) = [3, 4, 5]

  !> number of w found in dat files
  integer :: nw_found

  !> tag for MPI message
  integer, private :: itag

  !> coefficient for rk2_alfa
  double precision, private :: rk2_alfa


  !> whether staggered field is in dat
  logical, private :: stagger_mark_dat=.false.

  ! Formats used in output
  character(len=*), parameter :: fmt_r  = 'es16.8' ! Default precision
  character(len=*), parameter :: fmt_r2 = 'es10.2' ! Two digits
  character(len=*), parameter :: fmt_i  = 'i8'     ! Integer format

  ! public methods
  public :: snapshot_write_header

contains

  !> Read the command line arguments passed to agile
  subroutine read_arguments()
    use mod_global_parameters

    integer                          :: len, stat, n, i, ipars
    integer, parameter               :: max_files = 20 !Maximum number of par files
    integer                          :: n_par_files
    character(len=max_files*std_len) :: all_par_files
    character(len=std_len)           :: tmp_files(max_files), arg
    logical                          :: unknown_arg, help, morepars

    if (mype == 0) then
       print *, "-----------------------------------------------------------------------------"
       print *, "-----------------------------------------------------------------------------"
       print *, "|                                                                           |"
       print *, "|                                                                           |"
       print *, "|                                        d8, d8b                            |"
       print *, "|                                       `8P  88P                            |"
       print *, "|                                           d88                             |"
       print *, "|                   d888b8b   d888b8b    88b888   d8888b                    |"
       print *, "|                  d8P' ?88  d8P' ?88    88P?88  d8b_,dP                    |"
       print *, "|                  88b  ,88b 88b  ,88b  d88  88b 88b                        |"
       print *, "|                  `?88P'`88b`?88P'`88bd88'   88b`?888P'                    |"
       print *, "|                                   )88                                     |"
       print *, "|                                  ,88P                                     |"
       print *, "|                              `?8888P                                      |"
       print *, "|                                                                           |"
       print *, "-----------------------------------------------------------------------------"
       print *, "-----------------------------------------------------------------------------"
    end if

    ! =============== Fortran 2003 command line reading ================

    ! Default command line arguments
    all_par_files="agile.par"
    restart_from_file=undefined
    snapshotnext=-1
    slicenext=-1
    collapsenext=-1
    help=.false.
    convert=.false.
    resume_previous_run=.false.

    ! Argument 0 is program name, so we start from one
    i = 1
    unknown_arg=.false.
    DO
       CALL get_command_argument(i, arg)
       IF (LEN_TRIM(arg) == 0) EXIT
       select case(arg)
         case("-i")
           i = i+1
           CALL get_command_argument(i, arg)
           !if(mype==0)print *,'found argument=',arg
           all_par_files=TRIM(arg)
           morepars=.true.
           ipars=1
           do while (ipars<max_files.and.morepars)
              CALL get_command_argument(i+ipars,arg)
              !if(mype==0)print *,'found argument=',arg
              if (index(TRIM(arg),"-")==1.or.LEN_TRIM(arg)==0) then
                 morepars=.false.
              else
                 ipars=ipars+1
                 all_par_files=TRIM(all_par_files)//" "//TRIM(arg)
                 !if(mype==0)print *,'now all_par_files=',TRIM(all_par_files)
              endif
           end do
           !if(mype==0)print *,'-i identified ',ipars, ' par-file arguments'
           i=i+ipars-1
           !if(mype==0)print *,'-i arguments passed to all_par_files=',TRIM(all_par_files)
         case("-if")
           i = i+1
           CALL get_command_argument(i, arg)
           restart_from_file=TRIM(arg)
           !if(mype==0)print *,'-if has argument=',TRIM(arg),' passed to restart_from_file=',TRIM(restart_from_file)
         case("-slicenext")
           i = i+1
           CALL get_command_argument(i, arg)
           read(arg,*,iostat=stat) slicenext
           !if(mype==0)print *,'-slicenext has argument=',arg,' passed to slicenext=',slicenext
         case("-collapsenext")
           i = i+1
           CALL get_command_argument(i, arg)
           read(arg,*,iostat=stat) collapsenext
           !if(mype==0)print *,'-collapsenext has argument=',arg,' passed to collapsenext=',collapsenext
         case("-snapshotnext")
           i = i+1
           CALL get_command_argument(i, arg)
           read(arg,*,iostat=stat) snapshotnext
           !if(mype==0)print *,'-snapshotnext has argument=',arg,' passed to snapshotnext=',snapshotnext
         case("-resume")
           resume_previous_run=.true.
           !if(mype==0)print *,'resume specified: resume_previous_run=T'
         case("-convert")
           convert=.true.
           if(mype==0)print *,'convert specified: convert=T'
         case("--help","-help")
           help=.true.
           EXIT
         case default
           unknown_arg=.true.
           help=.true.
           EXIT
       end select
       i = i+1
    END DO

    if (unknown_arg) then
       print*,"======================================="
       print*,"Error: Command line argument ' ",TRIM(arg)," ' not recognized"
       print*,"======================================="
       help=.true.
    end if

    ! Show the usage if the help flag was given, or no par file was specified
    if (help) then
       if (mype == 0) then
          print *, 'Usage example:'
          print *, 'mpirun -np 4 ./agile -i file.par [file2.par ...]'
          print *, '         (later .par files override earlier ones)'
          print *, ''
          print *, 'Optional arguments:'
          print *, '-convert             Convert snapshot files'
          print *, '-if file0001.dat     Use this snapshot to restart from'
          print *, '                     (you can modify e.g. output names)'
          print *, '-resume              Automatically resume previous run'
          print *, '                     (useful for long runs on HPC systems)'
          print *, '-snapshotnext N      Manual index for next snapshot'
          print *, '-slicenext N         Manual index for next slice output'
          print *, '-collapsenext N      Manual index for next collapsed output'
          print *, ''
       end if
       call MPI_FINALIZE(ierrmpi)
       stop
    end if

    ! Split the input files, in case multiple were given
    call get_fields_string(all_par_files, " ,'"""//char(9), max_files,&
        tmp_files, n_par_files)

    allocate(par_files(n_par_files))
    par_files = tmp_files(1:n_par_files)

  end subroutine read_arguments

  !> Read in the user-supplied parameter-file
  subroutine read_par_files()
    use mod_global_parameters
    use mod_physics, only: phys_energy, physics_type, phys_wider_stencil
    use mod_limiter
    use mod_slice
    use mod_geometry
    use mod_source
    use mod_input_output_helper, only: get_names_from_string

    logical          :: fileopen, file_exists
    integer          :: i, j, k, ifile, io_state
    integer          :: iB, isave, iw, level, idim, islice, iside
    integer          :: nx_vec(3), block_nx_vec(3)
    integer          :: my_unit, iostate
    integer          :: ilev
    double precision :: dx_vec(3)

    character              :: c_ndim
    character(len=80)      :: fmt_string
    character(len=std_len) :: err_msg, restart_from_file_arg
    character(len=std_len) :: basename_full, basename_prev, dummy_file
    character(len=std_len), dimension(:), allocatable :: typeboundary_min1,&
       typeboundary_min2,typeboundary_min3, typeboundary_max1,&
       typeboundary_max2,typeboundary_max3
    character(len=std_len), allocatable :: limiter(:)
    character(len=std_len), allocatable :: gradient_limiter(:)
    !> How to apply dimensional splitting to the source terms, see
    !> @ref discretization.md
    character(len=std_len) :: typesourcesplit
    !> Which flux scheme of spatial discretization to use (per grid level)
    character(len=std_len), allocatable :: flux_scheme(:)
    !> Which time stepper to use
    character(len=std_len) :: time_stepper
    !> Which time integrator to use
    character(len=std_len) :: time_integrator
    !> type of curl operator
    character(len=std_len) :: typecurl
    !> Limiter used for prolongation to refined grids and ghost cells
    character(len=std_len) :: typeprolonglimit
    !> How to compute the CFL-limited time step.
    !> Options are 'maxsum': max(sum(c/dx)); 'summax': sum(max(c/dx)) and
    !> 'minimum: max(c/dx), where the summations loop over the grid dimensions and
    !> c is the velocity. The default 'maxsum' is the conventiontal way of
    !> computing CFL-limited time steps.
    character(len=std_len) :: typecourant

    double precision, dimension(nsavehi) :: tsave_log, tsave_dat, tsave_slice,&
        tsave_collapsed, tsave_custom
    double precision :: dtsave_log, dtsave_dat, dtsave_slice, dtsave_collapsed,&
        dtsave_custom
    integer :: ditsave_log, ditsave_dat, ditsave_slice, ditsave_collapsed,&
        ditsave_custom
    double precision :: tsavestart_log, tsavestart_dat, tsavestart_slice,&
        tsavestart_collapsed, tsavestart_custom
    integer :: ipower
    double precision :: sizeuniformpart1,sizeuniformpart2,sizeuniformpart3
    double precision :: im_delta,im_nu,rka54,rka51,rkb54,rka55

    namelist /filelist/ base_filename,restart_from_file, typefilelog,&
       firstprocess,reset_grid,snapshotnext, convert,convert_type,saveprim,&
       usr_filename,nwauxio,nocartesian, w_write,writelevel,writespshift,&
       length_convert_factor, w_convert_factor, time_convert_factor,level_io,&
       level_io_min, level_io_max, autoconvert,slice_type,slicenext,&
       collapsenext,collapse_type, type_endian

    namelist /savelist/ tsave,itsave,dtsave,ditsave,nslices,slicedir,&
        slicecoord,collapse,collapseLevel, time_between_print,tsave_log,&
        tsave_dat, tsave_slice, tsave_collapsed, tsave_custom, dtsave_log,&
        dtsave_dat, dtsave_slice, dtsave_collapsed, dtsave_custom, ditsave_log,&
        ditsave_dat, ditsave_slice, ditsave_collapsed, ditsave_custom,&
       tsavestart_log, tsavestart_dat, tsavestart_slice, tsavestart_collapsed,&
       tsavestart_custom, tsavestart

    namelist /stoplist/ it_init,time_init,it_max,time_max,dtmin,reset_it,&
       reset_time,wall_time_max,final_dt_reduction

    namelist /methodlist/ time_stepper,time_integrator, source_split_usr,&
       typesourcesplit,local_timestep,dimsplit,typedimsplit,flux_scheme,&
       limiter,gradient_limiter,cada3_radius,loglimit,&
       H_correction, typegrad,&
       typediv,typecurl,flatcd,flatsh,&
       rk2_alfa,imex222_lambda,ssprk_order,rk3_switch,imex_switch,&
       schmid_rad1,schmid_rad2,&
       schmid_rad3, phys,&
       flux_adaptive_diffusion, flux_ad_min, flux_ad_scale

    namelist /boundlist/ nghostcells,ghost_copy,internalboundary,&
        typeboundary_min1,typeboundary_min2,typeboundary_min3,&
       typeboundary_max1,typeboundary_max2,typeboundary_max3,&
        save_physical_boundary, specialboundary

    namelist /meshlist/ refine_max_level,nbufferx1,nbufferx2,nbufferx3,&
       refine_threshold,derefine_ratio, refine_criterion,&
        amr_wavefilter, max_blocks,&
        block_nx1,block_nx2,block_nx3,domain_nx1,&
       domain_nx2,domain_nx3,iprob,xprobmin1,xprobmin2,xprobmin3,xprobmax1,&
       xprobmax2,xprobmax3, w_refine_weight, prolongprimitive,coarsenprimitive,&
        typeprolonglimit, logflag,tfixgrid,itfixgrid,ditregrid, refine_usr,&
        geometry, log_r0
    namelist /paramlist/  courantpar, dtpar, dtdiffpar, typecourant, slowsteps

    namelist /emissionlist/ filename_euv,wavelength,filename_sxr,emin_sxr,&
       emax_sxr,LOS_theta,LOS_phi,image_rotate,x_origin,big_image,spectrum_wl,&
       location_slit,filename_spectrum,spectrum_window_min,spectrum_window_max,&
       instrument_resolution_factor,activate_unit_arcsec,filename_whitelight,&
       whitelight_instrument,R_occultor,R_opt_thick,dat_resolution,&
       direction_slit

    ! allocate cell size of all levels
    allocate(dx(ndim,nlevelshi))
    allocate(dg1(nlevelshi))
    allocate(dg2(nlevelshi))
    allocate(dg3(nlevelshi))
    allocate(ng1(nlevelshi))
    allocate(ng2(nlevelshi))
    allocate(ng3(nlevelshi))

    ! default block size excluding ghost cells
    block_nx1 = 16
    block_nx2 = 16
    block_nx3 = 16

    ! default resolution of level-1 mesh (full domain)
    domain_nx1 = 32
    domain_nx2 = 32
    domain_nx3 = 32

    !! default number of ghost-cell layers at each boundary of a block
    ! this is now done when the variable is defined in mod_global_parameters
    ! the physics modules might set this variable in their init subroutine called earlier
    nghostcells = 2

    ! Allocate boundary conditions arrays in new and old style

    allocate(typeboundary_min1(nwfluxbc))
    allocate(typeboundary_max1(nwfluxbc))
    typeboundary_min1 = undefined
    typeboundary_max1 = undefined

    allocate(typeboundary_min2(nwfluxbc))
    allocate(typeboundary_max2(nwfluxbc))
    typeboundary_min2 = undefined
    typeboundary_max2 = undefined

    allocate(typeboundary_min3(nwfluxbc))
    allocate(typeboundary_max3(nwfluxbc))
    typeboundary_min3 = undefined
    typeboundary_max3 = undefined


    allocate(typeboundary(nwflux+nwaux,2*ndim))
    typeboundary=0

    ! not save physical boundary in dat files by default
    save_physical_boundary = .false.

    internalboundary   = .false.

    ! defaults for specific options
    typegrad   = 'central'
    typediv    = 'central'
    typecurl   = 'central'

    ! defaults for convert behavior

    ! store the -if value from argument in command line
    restart_from_file_arg    = restart_from_file
    nwauxio                  = 0
    nocartesian              = .false.
    saveprim                 = .false.
    autoconvert              = .false.
    convert_type             = 'vtuBCCmpi'
    slice_type               = 'vtuCC'
    collapse_type            = 'vti'
    allocate(w_write(nw))
    w_write(1:nw)             = .true.
    allocate(writelevel(nlevelshi))
    writelevel(1:nlevelshi)  = .true.
    writespshift(1:ndim,1:2) = zero
    level_io                 = -1
    level_io_min             = 1
    level_io_max             = nlevelshi
    ! endianness: littleendian (default) is 1, bigendian otherwise
    type_endian              = 1

    ! normalization of primitive variables: only for output
    ! note that length_convert_factor is for length
    ! this scaling is optional, and must be set consistently if used
    allocate(w_convert_factor(nw))
    w_convert_factor(:)   = 1.0d0
    time_convert_factor   = 1.0d0
    length_convert_factor = 1.0d0

    ! AMR related defaults
    refine_max_level      = 1
    nbufferx1                 = 0
    nbufferx2                 = 0
    nbufferx3                 = 0
    allocate(refine_threshold(nlevelshi))
    refine_threshold(1:nlevelshi)            = 0.1d0
    allocate(derefine_ratio(nlevelshi))
    derefine_ratio(1:nlevelshi)       = 1.0d0/8.0d0
    typeprolonglimit            = 'default'
    refine_criterion               = 3
    allocate(w_refine_weight(nw+1))
    w_refine_weight                    = 0.d0
    allocate(logflag(nw+1))
    logflag                     = .false.
    allocate(amr_wavefilter(nlevelshi))
    amr_wavefilter(1:nlevelshi) = 1.0d-2
    tfixgrid                    = bigdouble
    itfixgrid                   = biginteger
    ditregrid                   = 1

    ! Offset of the logarithmic radial map; 0 selects the plain s = ln(r)
    log_r0 = 0.0d0

    ! IO defaults
    it_init       = 0
    it_max        = biginteger
    time_init     = 0.d0
    time_max      = bigdouble
    wall_time_max = bigdouble
    if(local_timestep) then
      final_dt_reduction=.false.
    else
      final_dt_reduction=.true.
    endif
    final_dt_exit=.false.
    dtmin         = 1.0d-10
    nslices       = 0
    collapse      = .false.
    collapseLevel = 1
    time_between_print = 30.0d0 ! Print status every 30 seconds

    do ifile=1,nfile
      do isave=1,nsavehi
        tsave(isave,ifile)  = bigdouble !global_time  of saves into the output files
        itsave(isave,ifile) = biginteger ! it of saves into the output files
      end do
      dtsave(ifile)  = bigdouble  ! time between saves
      ditsave(ifile) = biginteger ! timesteps between saves
      isavet(ifile)  = 1          ! index for saves by global_time
      isaveit(ifile) = 1          ! index for saves by it
      tsavestart(ifile) = 0.0d0
    end do

    tsave_log       = bigdouble
    tsave_dat       = bigdouble
    tsave_slice     = bigdouble
    tsave_collapsed = bigdouble
    tsave_custom    = bigdouble

    dtsave_log       = bigdouble
    dtsave_dat       = bigdouble
    dtsave_slice     = bigdouble
    dtsave_collapsed = bigdouble
    dtsave_custom    = bigdouble

    ditsave_log       = biginteger
    ditsave_dat       = biginteger
    ditsave_slice     = biginteger
    ditsave_collapsed = biginteger
    ditsave_custom    = biginteger

    tsavestart_log       = bigdouble
    tsavestart_dat       = bigdouble
    tsavestart_slice     = bigdouble
    tsavestart_collapsed = bigdouble
    tsavestart_custom    = bigdouble

    typefilelog = 'default'

    ! defaults for input
    reset_time = .false.
    reset_it = .false.
    firstprocess  = .false.
    reset_grid     = .false.
    base_filename   = 'data'
    usr_filename    = ''


    ! Defaults for discretization methods
    slowsteps       = -1
    courantpar      = 0.8d0
    typecourant     = 'maxsum'
    dimsplit        = .false.
    typedimsplit    = 'default'
    if(physics_type=='mhd') then
      cada3_radius  = 0.5d0
    else
      cada3_radius  = 0.1d0
    end if
    schmid_rad1 = 1.d0
    schmid_rad2 = 1.d0
    schmid_rad3 = 1.d0
    source_split_usr= .false.
    time_stepper    = 'twostep'
    time_integrator = 'default'
    ! default PC or explicit midpoint, hence alfa=0.5
    rk2_alfa        = half
    ! default IMEX-RK22Ln hence lambda = 1 - 1/sqrt(2)
    imex222_lambda  = 1.0d0 - 1.0d0 / dsqrt(2.0d0)
    ! default SSPRK(3,3) or Gottlieb-Shu 1998 for threestep
    ! default SSPRK(4,3) or Spireti-Ruuth for fourstep
    ! default SSPRK(5,4) using Gottlieb coeffs
    ssprk_order    = 3
    ! default RK3 butcher table: Heun 3rd order
    rk3_switch      = 3
    ! default IMEX threestep is IMEX_ARK(232)
    imex_switch     = 1

    ! Defaults for synthesing emission
    LOS_theta = 0.d0
    LOS_phi = 0.d0
    image_rotate = 0.d0
    x_origin = 0.d0
    big_image = .false.
    location_slit = 0.d0
    direction_slit = -1
    instrument_resolution_factor=1.d0
    activate_unit_arcsec=.true.
    whitelight_instrument='LASCO/C2'
    R_occultor=-1.d0
    R_opt_thick=1.d0
    dat_resolution=.false.

    allocate(flux_scheme(nlevelshi), flux_method(nlevelshi))
    allocate(limiter(nlevelshi),gradient_limiter(nlevelshi))
    do level=1,nlevelshi
       flux_scheme(level) = 'tvdlf'
       limiter(level)     = 'minmod'
       gradient_limiter(level) = 'minmod'
    end do

    flatcd          = .false.
    flatsh          = .false.
    typesourcesplit = 'sfs'
    allocate(loglimit(nw))
    loglimit(1:nw)  = .false.

    dtdiffpar     = 0.5d0
    dtpar         = -1.d0

    ! problem setup defaults
    iprob    = 1

    ! end defaults

    ! Initialize Kronecker delta, and Levi-Civita tensor
    do i=1,3
       do j=1,3
          if(i==j)then
             kr(i,j)=1
          else
             kr(i,j)=0
          endif
          do k=1,3
             if(i==j.or.j==k.or.k==i)then
                lvc(i,j,k)=0
             else if(i+1==j.or.i-2==j)then
                lvc(i,j,k)=1
             else
                lvc(i,j,k)=-1
             endif
          enddo
       enddo
    enddo
    ${GPU_UPDATE_DEVICE('kr,lvc')}$

    ! These are used to construct file and log names from multiple par files
    basename_full = ''
    basename_prev = ''

    ! Cray compiler bug: setting values after reading the namelists fails
    ! when reading files in a loop.
    ! To avoid this, disallow reading more than one file for now
#ifdef _CRAYFTN
    if (size(par_files) > 1) then
       write(err_msg, *) "Reading more than one par file is disallowed " // &
             "when compiled with Cray to avoid a compiler bug"
       call mpistop(trim(err_msg))
    end if
    if (size(par_files) > 0) then
       i = 1
#else
    do i = 1, size(par_files)
#endif
       if (mype == 0) print *, "Reading " // trim(par_files(i))

       ! Check whether the file exists
       inquire(file=trim(par_files(i)), exist=file_exists)

       if (.not. file_exists) then
          write(err_msg, *) "The parameter file " // trim(par_files(i)) // &
               " does not exist"
          call mpistop(trim(err_msg))
       end if

       open(unitpar, file=trim(par_files(i)), status='old')

       ! Try to read in the namelists. They can be absent or in a different
       ! order, since we rewind before each read.
       rewind(unitpar)
       read(unitpar, filelist, end=101)

101    rewind(unitpar)
       read(unitpar, savelist, end=102)

102    rewind(unitpar)
       read(unitpar, stoplist, end=103)

103    rewind(unitpar)
       read(unitpar, methodlist, end=104)

104    rewind(unitpar)
       read(unitpar, boundlist, end=105)

105    rewind(unitpar)
       read(unitpar, meshlist, end=106)

106    rewind(unitpar)
       read(unitpar, paramlist, end=107)

107    rewind(unitpar)
       read(unitpar, emissionlist, end=108)

108    close(unitpar)

       ! Append the log and file names given in the par files
       if (base_filename /= basename_prev) basename_full = trim(basename_full) &
          // trim(base_filename)
       basename_prev = base_filename
#ifdef _CRAYFTN
    end if
#else
    end do
#endif

    base_filename = basename_full

    
    ! Check whether output directory is writable
    if(mype==0) then
      dummy_file = trim(base_filename)//"DUMMY"
      open(newunit=my_unit, file=trim(dummy_file), iostat=iostate)
      if (iostate /= 0) then
         call mpistop("Can't write to output directory (" // &
              trim(base_filename) // ")")
      else
         close(my_unit, status='delete')
      end if
    end if

    if(source_split_usr) any_source_split=.true.

    ! restart filename from command line overwrites the one in par file
    if(restart_from_file_arg /= undefined) &
       restart_from_file=restart_from_file_arg

    ! Root process will search snapshot
    if (mype == 0) then
      if(restart_from_file == undefined) then
        ! search file from highest index
        file_exists=.false.
        do index_latest_data = 9999, 0, -1
          if(snapshot_exists(index_latest_data)) then
            file_exists=.true.
            exit
          end if
        end do
        if(.not.file_exists) index_latest_data=-1
      else
        ! get index of the given data restarted from
        index_latest_data=get_snapshot_index(trim(restart_from_file))
      end if
    end if
    call MPI_BCAST(index_latest_data, 1, MPI_INTEGER, 0, icomm, ierrmpi)

    if (resume_previous_run) then
      if (index_latest_data == -1) then
        if(mype==0) write(*,*)&
            "No snapshots found to resume from, start a new run..."
      else
        ! Set file name to restart from
        write(restart_from_file, "(a,i4.4,a)") trim(base_filename),&
           index_latest_data, ".dat"
      end if
    end if

    if (restart_from_file == undefined) then
      snapshotnext = 0
      slicenext    = 0
      collapsenext = 0
      if (firstprocess) call mpistop&
         ("Please restart from a snapshot when firstprocess=T")
      if (convert) call mpistop('Change convert to .false. for a new run!')
    end if

    if(convert) autoconvert=.false.

    where (tsave_log < bigdouble) tsave(:, 1) = tsave_log
    where (tsave_dat < bigdouble) tsave(:, 2) = tsave_dat
    where (tsave_slice < bigdouble) tsave(:, 3) = tsave_slice
    where (tsave_collapsed < bigdouble) tsave(:, 4) = tsave_collapsed
    where (tsave_custom < bigdouble) tsave(:, 5) = tsave_custom

    if (dtsave_log < bigdouble) dtsave(1) = dtsave_log
    if (dtsave_dat < bigdouble) dtsave(2) = dtsave_dat
    if (dtsave_slice < bigdouble) dtsave(3) = dtsave_slice
    if (dtsave_collapsed < bigdouble) dtsave(4) = dtsave_collapsed
    if (dtsave_custom < bigdouble) dtsave(5) = dtsave_custom

    if (tsavestart_log < bigdouble) tsavestart(1) = tsavestart_log
    if (tsavestart_dat < bigdouble) tsavestart(2) = tsavestart_dat
    if (tsavestart_slice < bigdouble) tsavestart(3) = tsavestart_slice
    if (tsavestart_collapsed < bigdouble) tsavestart(4) = tsavestart_collapsed
    if (tsavestart_custom < bigdouble) tsavestart(5) = tsavestart_custom

    if (ditsave_log < bigdouble) ditsave(1) = ditsave_log
    if (ditsave_dat < bigdouble) ditsave(2) = ditsave_dat
    if (ditsave_slice < bigdouble) ditsave(3) = ditsave_slice
    if (ditsave_collapsed < bigdouble) ditsave(4) = ditsave_collapsed
    if (ditsave_custom < bigdouble) ditsave(5) = ditsave_custom
    ! convert hours to seconds for ending wall time
    if (wall_time_max < bigdouble) wall_time_max=wall_time_max*3600.d0

    if (mype == 0) then
       write(unitterm, *) ''
       write(unitterm, *)&
'Output type | tsavestart |  dtsave   | ditsave | itsave(1) | tsave(1)'
       write(fmt_string, *) '(A12," | ",E9.3E2,"  | ",E9.3E2," | ",I6,"  | "'//&
            ',I6, "    | ",E9.3E2)'
    end if

    do ifile = 1, nfile
       if (mype == 0) write(unitterm, fmt_string) trim(output_names(ifile)),&
           tsavestart(ifile), dtsave(ifile), ditsave(ifile), itsave(1, ifile),&
           tsave(1, ifile)
    end do

    if (mype == 0) write(unitterm, *) ''

    do islice=1,nslices
       if(slicedir(islice) > ndim) write(uniterr,&
          *)'Warning in read_par_files: ', 'Slice ', islice,' direction',&
          slicedir(islice),'larger than ndim=',ndim
       if(slicedir(islice) < 1) write(uniterr,*)'Warning in read_par_files: ',&
           'Slice ', islice,' direction',slicedir(islice),&
          'too small, should be [',1,ndim,']'
    end do

    if(it_max==biginteger .and. time_max==bigdouble.and.mype==0) write(uniterr,&
       *) 'Warning in read_par_files: it_max or time_max not given!'

    select case (typecourant)
    case ('maxsum')
      type_courant=type_maxsum
    case ('summax')
      type_courant=type_summax
    case ('minimum')
      type_courant=type_minimum
    case default
       write(unitterm,*)'Unknown typecourant=',typecourant
       call mpistop("Error from read_par_files: no such typecourant!")
    end select

    ! =========================================================================
    ! Flux scheme parsing and physics compatibility validation
    ! =========================================================================
    do level=1,nlevelshi
       select case (flux_scheme(level))
       case ('hll')
          flux_method(level)=fs_hll
          if (physics_type=='ffhd') then
             call mpistop(trim(flux_scheme(level))//" flux scheme is not compatible with "//trim(physics_type))
          end if
       case ('hllc')
          flux_method(level)=fs_hllc
          if (physics_type=='ffhd' .or. physics_type=='srhd') then
             call mpistop(trim(flux_scheme(level))//" flux scheme is not compatible with "//trim(physics_type))
          end if
#:if defined('MHD_ENERGY_ONLY')
          ! only reconflux_muscl_llf_prim zeroes the non-energy fluxes
          call mpistop(trim(flux_scheme(level))//" flux scheme is not compatible with MHD_ENERGY_ONLY")
#:endif
       case ('tvdlf')
          flux_method(level)=fs_tvdlf
       case default
          call mpistop("unknown or bad flux scheme: "//trim(flux_scheme(level)))
       end select
    end do

    ! =========================================================================
    ! Limiter and spatial operator parsing
    ! =========================================================================
    if(any(limiter=='schmid1')) need_global_a2max=.true.

    ! initialize type_curl
     select case (typecurl)
     case ("central")
        type_curl=central
     case ("Gaussbased")
        type_curl=Gaussbased
     case ("Stokesbased")
        type_curl=Stokesbased
     case default
        write(unitterm,*) "typecurl=",typecurl
        call mpistop("unkown type of curl operator in read_par_files")
     end select

    ! initialize types of time stepper and time integrator
    select case (time_stepper)
    case ("onestep")
       t_stepper=onestep
       if (time_integrator=='default') then
          time_integrator="Forward_Euler"
       end if
       select case (time_integrator)
       case ("Forward_Euler")
          t_integrator=Forward_Euler
       case ("IMEX_Euler")
          t_integrator=IMEX_Euler
       case ("IMEX_SP")
          t_integrator=IMEX_SP
       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",&
             time_stepper
          call mpistop("unkown onestep time_integrator in read_par_files")
       end select
       use_imex_scheme=(t_integrator==IMEX_Euler.or.t_integrator==IMEX_SP)
    case ("twostep")
       t_stepper=twostep
       if (time_integrator=='default') then
          time_integrator="Predictor_Corrector"
       endif
       select case (time_integrator)
       case ("Predictor_Corrector")
          t_integrator=Predictor_Corrector
       case ("RK2_alfa")
          t_integrator=RK2_alf
       case ("ssprk2")
          t_integrator=ssprk2
       case ("IMEX_Midpoint")
          t_integrator=IMEX_Midpoint
       case ("IMEX_Trapezoidal")
          t_integrator=IMEX_Trapezoidal
       case ("IMEX_222")
          t_integrator=IMEX_222
       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",&
             time_stepper
          call mpistop("unkown twostep time_integrator in read_par_files")
       end select
       use_imex_scheme=(t_integrator==IMEX_Midpoint.or.t_integrator==&
          IMEX_Trapezoidal.or.t_integrator==IMEX_222)
       if (t_integrator==RK2_alf) then
          if(rk2_alfa<smalldouble.or.rk2_alfa>one)call &
             mpistop("set rk2_alfa within [0,1]")
          rk_a21=rk2_alfa
          rk_b2=1.0d0/(2.0d0*rk2_alfa)
          rk_b1=1.0d0-rk_b2
       endif
    case ("threestep")
       t_stepper=threestep
       if (time_integrator=='default') then
          time_integrator='ssprk3'
       endif
       select case (time_integrator)
       case ("ssprk3")
          t_integrator=ssprk3
       case ("RK3_BT")
          t_integrator=RK3_BT
       case ("IMEX_ARS3")
          t_integrator=IMEX_ARS3
       case ("IMEX_232")
          t_integrator=IMEX_232
       case ("IMEX_CB3a")
          t_integrator=IMEX_CB3a
       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",&
             time_stepper
          call mpistop("unkown threestep time_integrator in read_par_files")
       end select
       if(t_integrator==RK3_BT) then
           select case(rk3_switch)
             case(1)
              ! we code up Ralston 3rd order here
              rk3_a21=1.0d0/2.0d0
              rk3_a31=0.0d0
              rk3_a32=3.0d0/4.0d0
              rk3_b1=2.0d0/9.0d0
              rk3_b2=1.0d0/3.0d0
             case(2)
              ! we code up RK-Wray 3rd order here
              rk3_a21=8.0d0/15.0d0
              rk3_a31=1.0d0/4.0d0
              rk3_a32=5.0d0/12.0d0
              rk3_b1=1.0d0/4.0d0
              rk3_b2=0.0d0
             case(3)
              ! we code up Heun 3rd order here
              rk3_a21=1.0d0/3.0d0
              rk3_a31=0.0d0
              rk3_a32=2.0d0/3.0d0
              rk3_b1=1.0d0/4.0d0
              rk3_b2=0.0d0
             case(4)
              ! we code up Nystrom 3rd order here
              rk3_a21=2.0d0/3.0d0
              rk3_a31=0.0d0
              rk3_a32=2.0d0/3.0d0
              rk3_b1=1.0d0/4.0d0
              rk3_b2=3.0d0/8.0d0
             case default
                call mpistop("Unknown rk3_switch")
            end select
           ! the rest is fixed from above
           rk3_b3=1.0d0-rk3_b1-rk3_b2
           rk3_c2=rk3_a21
           rk3_c3=rk3_a31+rk3_a32
       endif
       if(t_integrator==ssprk3) then
         select case(ssprk_order)
             case(3) ! this is SSPRK(3,3) Gottlieb-Shu
                rk_beta11=1.0d0
                rk_beta22=1.0d0/4.0d0
                rk_beta33=2.0d0/3.0d0
                rk_alfa21=3.0d0/4.0d0
                rk_alfa31=1.0d0/3.0d0
                rk_c2=1.0d0
                rk_c3=1.0d0/2.0d0
             case(2) ! this is SSP(3,2)
                rk_beta11=1.0d0/2.0d0
                rk_beta22=1.0d0/2.0d0
                rk_beta33=1.0d0/3.0d0
                rk_alfa21=0.0d0
                rk_alfa31=1.0d0/3.0d0
                rk_c2=1.0d0/2.0d0
                rk_c3=1.0d0
             case default
                call mpistop("Unknown ssprk3_order")
         end select
         rk_alfa22=1.0d0-rk_alfa21
         rk_alfa33=1.0d0-rk_alfa31
       endif
       if(t_integrator==IMEX_ARS3) then
          ars_gamma=(3.0d0+dsqrt(3.0d0))/6.0d0
       endif
       if(t_integrator==IMEX_232) then
           select case(imex_switch)
             case(1) ! this is IMEX_ARK(232)
              im_delta=1.0d0-1.0d0/dsqrt(2.0d0)
              im_nu=(3.0d0+2.0d0*dsqrt(2.0d0))/6.0d0
              imex_a21=2.0d0*im_delta
              imex_a31=1.0d0-im_nu
              imex_a32=im_nu
              imex_b1=1.0d0/(2.0d0*dsqrt(2.0d0))
              imex_b2=1.0d0/(2.0d0*dsqrt(2.0d0))
              imex_ha21=im_delta
              imex_ha22=im_delta
             case(2) ! this is IMEX_SSP(232)
              ! doi 10.1002/2017MS001065 Rokhzadi et al
              imex_a21=0.711664700366941d0
              imex_a31=0.077338168947683d0
              imex_a32=0.917273367886007d0
              imex_b1=0.398930808264688d0
              imex_b2=0.345755244189623d0
              imex_ha21=0.353842865099275d0
              imex_ha22=0.353842865099275d0
             case default
                call mpistop("Unknown imex_siwtch")
           end select
           imex_c2=imex_a21
           imex_c3=imex_a31+imex_a32
           imex_b3=1.0d0-imex_b1-imex_b2
       endif
       if(t_integrator==IMEX_CB3a) then
          imex_c2   = 0.8925502329346865
          imex_a22  = imex_c2
          imex_ha21 = imex_c2
          imex_c3   = imex_c2 / (6.0d0*imex_c2**2 - 3.0d0*imex_c2 + 1.0d0)
          imex_ha32 = imex_c3
          imex_b2   = (3.0d0*imex_c2 - 1.0d0) / (6.0d0*imex_c2**2)
          imex_b3   = (6.0d0*imex_c2**2 - 3.0d0*imex_c2 + 1.0d0) / &
             (6.0d0*imex_c2**2)
          imex_a33  = (1.0d0/6.0d0 - imex_b2*imex_c2**2 - &
             imex_b3*imex_c2*imex_c3) / (imex_b3*(imex_c3-imex_c2))
          imex_a32  = imex_c3 - imex_a33
          ! if (mype == 0) then
          !    write(*,*) "================================="
          !    write(*,*) "Asserting the order conditions..."
          !    ! First order convergence: OK
          !    ! Second order convergence
          !    write(*,*) -1.0d0/2.0d0 + imex_b2*imex_c2 + imex_b3*imex_c3
          !    ! Third order convergence
          !    write(*,*) -1.0d0/3.0d0 + imex_b2*imex_c2**2 + imex_b3*imex_c3**2
          !    write(*,*) -1.0d0/6.0d0 + imex_b3*imex_ha32*imex_c2
          !    write(*,*) -1.0d0/6.0d0 + imex_b2*imex_a22*imex_c2 + imex_b3*imex_a32*imex_c2 + imex_b3*imex_a33*imex_c3
          !    write(*,*) "================================="
          ! end if
       end if
       use_imex_scheme=(t_integrator==IMEX_ARS3.or.t_integrator==&
          IMEX_232.or.t_integrator==IMEX_CB3a)
    case ("fourstep")
       t_stepper=fourstep
       if (time_integrator=='default') then
          time_integrator="ssprk4"
       endif
       select case (time_integrator)
       case ("ssprk4")
          t_integrator=ssprk4
       case ("rk4")
          t_integrator=rk4
       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",&
             time_stepper
          call mpistop("unkown fourstep time_integrator in read_par_files")
       end select
       if(t_integrator==ssprk4) then
         select case(ssprk_order)
             case(3) ! this is SSPRK(4,3) Spireti-Ruuth
                rk_beta11=1.0d0/2.0d0
                rk_beta22=1.0d0/2.0d0
                rk_beta33=1.0d0/6.0d0
                rk_beta44=1.0d0/2.0d0
                rk_alfa21=0.0d0
                rk_alfa31=2.0d0/3.0d0
                rk_alfa41=0.0d0
                rk_c2=1.0d0/2.0d0
                rk_c3=1.0d0
                rk_c4=1.0d0/2.0d0
             case(2) ! this is SSP(4,2)
                rk_beta11=1.0d0/3.0d0
                rk_beta22=1.0d0/3.0d0
                rk_beta33=1.0d0/3.0d0
                rk_beta44=1.0d0/4.0d0
                rk_alfa21=0.0d0
                rk_alfa31=0.0d0
                rk_alfa41=1.0d0/4.0d0
                rk_c2=1.0d0/3.0d0
                rk_c3=2.0d0/3.0d0
                rk_c4=1.0d0
             case default
                call mpistop("Unknown ssprk_order")
         end select
         rk_alfa22=1.0d0-rk_alfa21
         rk_alfa33=1.0d0-rk_alfa31
         rk_alfa44=1.0d0-rk_alfa41
       endif
    case ("fivestep")
       t_stepper=fivestep
       if (time_integrator=='default') then
          time_integrator="ssprk5"
       end if
       select case (time_integrator)
       case ("ssprk5")
          t_integrator=ssprk5
       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",&
             time_stepper
          call mpistop("unkown fivestep time_integrator in read_par_files")
       end select
       if(t_integrator==ssprk5) then
         select case(ssprk_order)
           ! we use ssprk_order to intercompare the different coefficient choices
           case(3) ! From Gottlieb 2005
            rk_beta11=0.391752226571890d0
            rk_beta22=0.368410593050371d0
            rk_beta33=0.251891774271694d0
            rk_beta44=0.544974750228521d0
            rk_beta54=0.063692468666290d0
            rk_beta55=0.226007483236906d0
            rk_alfa21=0.444370493651235d0
            rk_alfa31=0.620101851488403d0
            rk_alfa41=0.178079954393132d0
            rk_alfa53=0.517231671970585d0
            rk_alfa54=0.096059710526147d0

            rk_alfa22=0.555629506348765d0
            rk_alfa33=0.379898148511597d0
            rk_alfa44=0.821920045606868d0
            rk_alfa55=0.386708617503269d0
            rk_alfa22=1.0d0-rk_alfa21
            rk_alfa33=1.0d0-rk_alfa31
            rk_alfa44=1.0d0-rk_alfa41
            rk_alfa55=1.0d0-rk_alfa53-rk_alfa54
            rk_c2=rk_beta11
            rk_c3=rk_alfa22*rk_c2+rk_beta22
            rk_c4=rk_alfa33*rk_c3+rk_beta33
            rk_c5=rk_alfa44*rk_c4+rk_beta44
           case(2) ! From Spireti-Ruuth
            rk_beta11=0.39175222700392d0
            rk_beta22=0.36841059262959d0
            rk_beta33=0.25189177424738d0
            rk_beta44=0.54497475021237d0
            !rk_beta54=0.06369246925946d0
            !rk_beta55=0.22600748319395d0
            rk_alfa21=0.44437049406734d0
            rk_alfa31=0.62010185138540d0
            rk_alfa41=0.17807995410773d0
            rk_alfa53=0.51723167208978d0
            !rk_alfa54=0.09605971145044d0

            rk_alfa22=1.0d0-rk_alfa21
            rk_alfa33=1.0d0-rk_alfa31
            rk_alfa44=1.0d0-rk_alfa41

            rka51=0.00683325884039d0
            rka54=0.12759831133288d0
            rkb54=0.08460416338212d0
            rk_beta54=rkb54-rk_beta44*rka51/rk_alfa41
            rk_alfa54=rka54-rk_alfa44*rka51/rk_alfa41

            rk_alfa55=1.0d0-rk_alfa53-rk_alfa54
            rk_c2=rk_beta11
            rk_c3=rk_alfa22*rk_c2+rk_beta22
            rk_c4=rk_alfa33*rk_c3+rk_beta33
            rk_c5=rk_alfa44*rk_c4+rk_beta44
            rk_beta55=1.0d0-rk_beta54-rk_alfa53*rk_c3-rk_alfa54*rk_c4-&
               rk_alfa55*rk_c5
           case default
                call mpistop("Unknown ssprk_order")
         end select
         ! the following combinations must be unity
         !print *,rk_beta55+rk_beta54+rk_alfa53*rk_c3+rk_alfa54*rk_c4+rk_alfa55*rk_c5
         !print *,rk_alfa22+rk_alfa21
         !print *,rk_alfa33+rk_alfa31
         !print *,rk_alfa44+rk_alfa41
         !print *,rk_alfa55+rk_alfa53+rk_alfa54
       endif
       use_imex_scheme=.false.
    case default
       call mpistop("Unknown time_stepper in read_par_files")
    end select

    ! The coefficients above are kept for every upstream scheme, but advect
    ! implements only the explicit ones; the IMEX schemes would need a
    ! physics module with an implicit update.  Stop here rather than at the
    ! first time step.
    select case (t_integrator)
    case (Forward_Euler, Predictor_Corrector, RK2_alf, ssprk2, ssprk3, &
          RK3_BT, ssprk4, rk4, ssprk5)
    case default
       call mpistop("time_integrator='"//trim(time_integrator)//&
          "' is not implemented in AGILE")
    end select
    call time_levels(t_integrator, nstep, nlevel)

    ! Dimensional splitting is not implemented: advect sweeps all three
    ! directions whatever idim range it is given, so the split loop in
    ! advance would apply the full flux divergence ndim times per step.
    ! typedimsplit='xy'/'xyyx' used to switch it on by itself, so refuse both.
    if (dimsplit) call mpistop("dimsplit=.true. is not implemented in AGILE: &
       &the flux update is always unsplit, leave dimsplit=.false.")
    if (typedimsplit /= 'default' .and. typedimsplit /= 'unsplit') &
       call mpistop("typedimsplit='"//trim(typedimsplit)//"' would switch &
       &on dimensional splitting, which is not implemented in AGILE")
    typedimsplit = 'unsplit'

    ! initialize types of split-source addition
    select case (typesourcesplit)
    case ('sfs')
      sourcesplit=sourcesplit_sfs
    case ('sf')
      sourcesplit=sourcesplit_sf
    case ('ssf')
      sourcesplit=sourcesplit_ssf
    case ('ssfss')
      sourcesplit=sourcesplit_ssfss
    case default
       write(unitterm,*)'No such typesourcesplit=',typesourcesplit
       call mpistop("Error: Unknown typesourcesplit!")
    end select

    if(coordinate==-1) then
      coordinate=Cartesian
      if(mype==0) then
        write(*,*) 'Warning: coordinate system is not specified!'
        write(*,*) 'call set_coordinate_system in usr_init in mod_usr.t'
        write(*,*) 'Now use Cartesian coordinate'
      end if
    end if

    if(coordinate==Cartesian) then
      slab=.true.
      slab_uniform=.true.
    else
      slab=.false.
      slab_uniform=.false.
    end if


    ! type limiter of prolongation
    select case(typeprolonglimit)
    case('unlimit')
      ! unlimited
      prolong_limiter=1
    case('minmod')
      prolong_limiter=2
    case('woodward')
      prolong_limiter=3
    case('koren')
      prolong_limiter=4
    case default
      prolong_limiter=0
    end select

    ! Type limiter is of integer type for performance
    allocate(type_limiter(nlevelshi))
    allocate(type_gradient_limiter(nlevelshi))

    do level=1,nlevelshi
       type_limiter(level) = limiter_type(limiter(level))
       type_gradient_limiter(level) = limiter_type(gradient_limiter(level))
    end do

    if (any(limiter(1:nlevelshi)== 'ppm').and.(flatsh.and.physics_type==&
       'rho')) then
       call mpistop&
          (" PPM with flatsh=.true. can not be used with physics_type='rho'!")
    end if

    ! Copy boundary conditions to typeboundary, which is used internally

    do iw=1,nwfluxbc
      select case(typeboundary_min1(iw))
      case("special")
        typeboundary(iw,2*1-1)=bc_special
      case("cont")
        typeboundary(iw,2*1-1)=bc_cont
      case("symm")
        typeboundary(iw,2*1-1)=bc_symm
      case("asymm")
        typeboundary(iw,2*1-1)=bc_asymm
      case("periodic")
        typeboundary(iw,2*1-1)=bc_periodic
      case("aperiodic")
        typeboundary(iw,2*1-1)=bc_aperiodic
      case("noinflow")
        typeboundary(iw,2*1-1)=bc_noinflow
      case("pole")
        typeboundary(iw,2*1-1)=12
      case("bc_data")
        typeboundary(iw,2*1-1)=bc_data
      case("bc_icarus")
        typeboundary(iw,2*1-1)=bc_icarus
      case("character")
        typeboundary(iw,2*1-1)=bc_character
      case default
         write (unitterm,*) "Undefined boundarytype found in read_par_files",&
             typeboundary_min1(iw),"for variable iw=",iw," and side iB=",2*1-1
      end select
    end do
    do iw=1,nwfluxbc
      select case(typeboundary_max1(iw))
      case("special")
        typeboundary(iw,2*1)=bc_special
      case("cont")
        typeboundary(iw,2*1)=bc_cont
      case("symm")
        typeboundary(iw,2*1)=bc_symm
      case("asymm")
        typeboundary(iw,2*1)=bc_asymm
      case("periodic")
        typeboundary(iw,2*1)=bc_periodic
      case("aperiodic")
        typeboundary(iw,2*1)=bc_aperiodic
      case("noinflow")
        typeboundary(iw,2*1)=bc_noinflow
      case("pole")
        typeboundary(iw,2*1)=12
      case("bc_data")
        typeboundary(iw,2*1)=bc_data
      case("bc_icarus")
        typeboundary(iw,2*1-1)=bc_icarus
      case("bc_character")
        typeboundary(iw,2*1)=bc_character
      case default
         write (unitterm,*) "Undefined boundarytype found in read_par_files",&
             typeboundary_max1(iw),"for variable iw=",iw," and side iB=",2*1
      end select
    end do

    do iw=1,nwfluxbc
      select case(typeboundary_min2(iw))
      case("special")
        typeboundary(iw,2*2-1)=bc_special
      case("cont")
        typeboundary(iw,2*2-1)=bc_cont
      case("symm")
        typeboundary(iw,2*2-1)=bc_symm
      case("asymm")
        typeboundary(iw,2*2-1)=bc_asymm
      case("periodic")
        typeboundary(iw,2*2-1)=bc_periodic
      case("aperiodic")
        typeboundary(iw,2*2-1)=bc_aperiodic
      case("noinflow")
        typeboundary(iw,2*2-1)=bc_noinflow
      case("pole")
        typeboundary(iw,2*2-1)=12
      case("bc_data")
        typeboundary(iw,2*2-1)=bc_data
      case("bc_icarus")
        typeboundary(iw,2*2-1)=bc_icarus
      case("character")
        typeboundary(iw,2*2-1)=bc_character
      case default
         write (unitterm,*) "Undefined boundarytype found in read_par_files",&
             typeboundary_min2(iw),"for variable iw=",iw," and side iB=",2*2-1
      end select
    end do
    do iw=1,nwfluxbc
      select case(typeboundary_max2(iw))
      case("special")
        typeboundary(iw,2*2)=bc_special
      case("cont")
        typeboundary(iw,2*2)=bc_cont
      case("symm")
        typeboundary(iw,2*2)=bc_symm
      case("asymm")
        typeboundary(iw,2*2)=bc_asymm
      case("periodic")
        typeboundary(iw,2*2)=bc_periodic
      case("aperiodic")
        typeboundary(iw,2*2)=bc_aperiodic
      case("noinflow")
        typeboundary(iw,2*2)=bc_noinflow
      case("pole")
        typeboundary(iw,2*2)=12
      case("bc_data")
        typeboundary(iw,2*2)=bc_data
      case("bc_icarus")
        typeboundary(iw,2*2-1)=bc_icarus
      case("bc_character")
        typeboundary(iw,2*2)=bc_character
      case default
         write (unitterm,*) "Undefined boundarytype found in read_par_files",&
             typeboundary_max2(iw),"for variable iw=",iw," and side iB=",2*2
      end select
    end do

    do iw=1,nwfluxbc
      select case(typeboundary_min3(iw))
      case("special")
        typeboundary(iw,2*3-1)=bc_special
      case("cont")
        typeboundary(iw,2*3-1)=bc_cont
      case("symm")
        typeboundary(iw,2*3-1)=bc_symm
      case("asymm")
        typeboundary(iw,2*3-1)=bc_asymm
      case("periodic")
        typeboundary(iw,2*3-1)=bc_periodic
      case("aperiodic")
        typeboundary(iw,2*3-1)=bc_aperiodic
      case("noinflow")
        typeboundary(iw,2*3-1)=bc_noinflow
      case("pole")
        typeboundary(iw,2*3-1)=12
      case("bc_data")
        typeboundary(iw,2*3-1)=bc_data
      case("bc_icarus")
        typeboundary(iw,2*3-1)=bc_icarus
      case("character")
        typeboundary(iw,2*3-1)=bc_character
      case default
         write (unitterm,*) "Undefined boundarytype found in read_par_files",&
             typeboundary_min3(iw),"for variable iw=",iw," and side iB=",2*3-1
      end select
    end do
    do iw=1,nwfluxbc
      select case(typeboundary_max3(iw))
      case("special")
        typeboundary(iw,2*3)=bc_special
      case("cont")
        typeboundary(iw,2*3)=bc_cont
      case("symm")
        typeboundary(iw,2*3)=bc_symm
      case("asymm")
        typeboundary(iw,2*3)=bc_asymm
      case("periodic")
        typeboundary(iw,2*3)=bc_periodic
      case("aperiodic")
        typeboundary(iw,2*3)=bc_aperiodic
      case("noinflow")
        typeboundary(iw,2*3)=bc_noinflow
      case("pole")
        typeboundary(iw,2*3)=12
      case("bc_data")
        typeboundary(iw,2*3)=bc_data
      case("bc_icarus")
        typeboundary(iw,2*3-1)=bc_icarus
      case("bc_character")
        typeboundary(iw,2*3)=bc_character
      case default
         write (unitterm,*) "Undefined boundarytype found in read_par_files",&
             typeboundary_max3(iw),"for variable iw=",iw," and side iB=",2*3
      end select
    end do


    ! psi, tracers take the same boundary type as the first variable
    if (nwfluxbc<nwflux) then
      do iw = nwfluxbc+1, nwflux
        typeboundary(iw,:) = typeboundary(1, :)
      end do
    end if
    ! auxiliary variables take the same boundary type as the first variable
    if (nwaux>0) then
      do iw = nwflux+1, nwflux+nwaux
        typeboundary(iw,:) = typeboundary(1, :)
      end do
    end if

    ! the heat-flux vector has no boundary types of its own: at a reflecting
    ! wall its normal component flips with the normal momentum
    if (allocated(iw_qvec)) then
      do idim=1,ndim
        do iside=1,2
          iB=2*idim-2+iside
          if (typeboundary(iw_mom(idim),iB)==bc_asymm) &
             typeboundary(iw_qvec(idim),iB)=bc_asymm
        end do
      end do
    end if

    if (any(typeboundary == 0)) then
      call mpistop("Not all boundary conditions have been defined")
    end if

    do idim=1,ndim
       periodB(idim)=(any(typeboundary(:,2*idim-1:2*idim)==bc_periodic))
       aperiodB(idim)=(any(typeboundary(:,2*idim-1:2*idim)==bc_aperiodic))
       if (periodB(idim).or.aperiodB(idim)) then
          do iw=1,nwflux
             if (typeboundary(iw,2*idim-1) .ne. typeboundary(iw,&
                2*idim)) call mpistop("Wrong counterpart in periodic boundary")


             if (typeboundary(iw,2*idim-1) /= bc_periodic .and. &
                typeboundary(iw,2*idim-1) /= bc_aperiodic) then
               call mpistop("Each dimension should either have all "//&
                    "or no variables periodic, some can be aperiodic")
             end if
          end do
       end if
    end do

    ! Expand the 'pole' boundary type (parsed as the sentinel 12 above) into
    ! ordinary per-variable symm/asymm entries.  The axis is not a boundary: a
    ! ghost cell beyond it is the interior cell half a revolution away in phi,
    ! and getbc fetches it from there.  All this table has to say is which
    ! components are odd under that pi-rotation, and the rule is the same in
    ! both coordinate systems: the component along the mirrored direction
    ! flips, the phi component flips, everything else is symmetric.
    !
    !   spherical,   axis at theta=0 or pi (mirror in theta): m_theta, m_phi
    !   cylindrical, axis at r=0           (mirror in r):     m_r,     m_phi
    !
    ! and likewise for the magnetic field in an mhd build.  Note the
    ! cylindrical radial entry: upstream MPI-AMRVAC marks only m_phi (and
    ! B_phi) antisymmetric at the cylindrical axis, but the axis maps
    ! (r, phi) -> (r, phi+pi), under which e_r(phi+pi) = -e_r(phi) exactly as
    ! e_phi(phi+pi) = -e_phi(phi), so the radial component is odd too.
    ! tests/hd/cylindrical_pole_uniform_flow is what pins this down.
    !
    ! The indices come from iw_mom/iw_mag rather than from the positional
    ! rho-mom-[e]-B arithmetic upstream uses, because a HYPERTC build
    ! registers q_ between e_ and mag and shifts every field index.
    do idim=1,ndim
      do iside=1,2
        iB=2*idim-2+iside
        if(.not.any(typeboundary(:,iB)==12)) cycle
        poleB_requested(iside,idim)=.true.
        ! a pole is a property of the face, not of one variable
        if(any(typeboundary(:,iB)/=12)) call mpistop("typeboundary='pole' has &
           &to be given for every variable on a face")
        select case(physics_type)
        case('hd','mhd','srhd','srmhd')
          ! the vector components handled below are all there is to say. The
          ! relativistic modules' mom(:) is the spatial four-velocity
          ! u^i = lfac*v^i (see con2prim); lfac is a scalar that a coordinate
          ! rotation leaves invariant, so u^i transforms exactly like an
          ! ordinary vector under the pole's pi-rotation and takes the same
          ! asymm pattern as hd's momentum. srmhd's auxiliaries mu and lfac,
          ! which it exchanges through getbc, and its GLM psi are scalars and
          ! stay symm.
        case('ffhd')
          ! ffhd's only conserved quantities are scalars - rho, the
          ! field-aligned momentum mom(1)=m_par, the energy - and a scalar is
          ! even under the pole's proper pi-rotation (m_par = rho*(v.bhat) is a
          ! contraction of two vectors), so every getbc-exchanged variable
          ! stays symm and the momentum lines below are skipped. The frozen
          ! field b1,b2,b3 is a vector, but it is not exchanged through getbc:
          ! fill_nwextra_device re-derives it in every ghost cell, pole
          ! ghosts included, so it needs no entry here either.
        case default
          call mpistop('pole treatment is only implemented for phys=hd, &
             &phys=mhd, phys=srhd, phys=srmhd and phys=ffhd')
        end select
        select case(coordinate)
        case(cylindrical)
          if(idim/=1) call mpistop('the cylindrical axis is at r=0, so only &
             &typeboundary_min1 can be pole')
        case(spherical)
          if(idim/=2) call mpistop('the spherical poles are at theta=0 and &
             &theta=pi, so only typeboundary_min2/max2 can be pole')
        case default
          call mpistop('Pole is in cylindrical, polar, spherical coordinates!')
        end select
        typeboundary(:,iB)=bc_symm
        if(physics_type/='ffhd') then
          ! ffhd's mom has a single element (the scalar m_par); hd/mhd/srhd
          ! carry a genuine ndir-vector whose mirrored and phi components flip
          typeboundary(iw_mom(idim),iB)=bc_asymm
          typeboundary(iw_mom(phi_),iB)=bc_asymm
        end if
        if(physics_type=='mhd' .or. physics_type=='srmhd') then
          typeboundary(iw_mag(idim),iB)=bc_asymm
          typeboundary(iw_mag(phi_),iB)=bc_asymm
          if(allocated(iw_qvec)) then
            typeboundary(iw_qvec(idim),iB)=bc_asymm
            typeboundary(iw_qvec(phi_),iB)=bc_asymm
          end if
        end if
      end do
    end do

   if (any(typeboundary(:,:)==bc_special) .and. .not. specialboundary) then
      call mpistop('special boundary requested, set specialboundary=.true.')
   end if


    if(.not.phys_energy) then
      flatcd=.false.
      flatsh=.false.
    end if

    if(any(limiter(1:nlevelshi)=='mp5')) then
      nghostcells=max(nghostcells,3)
    end if

    if(any(limiter(1:nlevelshi)=='weno5')) then
      nghostcells=max(nghostcells,3)
    end if

    if(any(limiter(1:nlevelshi)=='weno5nm')) then
      nghostcells=max(nghostcells,3)
    end if

    if(any(limiter(1:nlevelshi)=='wenoz5')) then
      nghostcells=max(nghostcells,3)
    end if

    if(any(limiter(1:nlevelshi)=='wenoz5nm')) then
      nghostcells=max(nghostcells,3)
    end if

    if(any(limiter(1:nlevelshi)=='wenozp5')) then
      nghostcells=max(nghostcells,3)
    end if

    if(any(limiter(1:nlevelshi)=='wenozp5nm')) then
      nghostcells=max(nghostcells,3)
    end if

    if(any(limiter(1:nlevelshi)=='teno5ad')) then
      nghostcells=max(nghostcells,3)
    end if

    if(any(limiter(1:nlevelshi)=='weno5cu6')) then
      nghostcells=max(nghostcells,3)
    end if

    if(any(limiter(1:nlevelshi)=='ppm')) then
      if(flatsh .or. flatcd) then
        nghostcells=max(nghostcells,4)
      else
        nghostcells=max(nghostcells,3)
      end if
    end if

    if(any(limiter(1:nlevelshi)=='weno7')) then
      nghostcells=max(nghostcells,4)
    end if

    if(any(limiter(1:nlevelshi)=='mpweno7')) then
      nghostcells=max(nghostcells,4)
    end if

    ! If a wider stencil is used, extend the number of ghost cells
    nghostcells = nghostcells + phys_wider_stencil

    ! prolongation in AMR for constrained transport MHD needs even number ghosts
    if(stagger_grid .and. refine_max_level>1 .and. mod(nghostcells,2)/=0) then
      nghostcells=nghostcells+1
   end if
   ${GPU_UPDATE_DEVICE('nghostcells')}$

      select case (coordinate)

      case (spherical)
         xprobmin2=xprobmin2*two*dpi;xprobmin3=xprobmin3*two*dpi
         xprobmax2=xprobmax2*two*dpi;xprobmax3=xprobmax3*two*dpi;

      case (cylindrical)

         if (1==phi_) then
            xprobmin1=xprobmin1*two*dpi;xprobmax1=xprobmax1*two*dpi;
         end if


         if (2==phi_) then
            xprobmin2=xprobmin2*two*dpi;xprobmax2=xprobmax2*two*dpi;
         end if


         if (3==phi_) then
            xprobmin3=xprobmin3*two*dpi;xprobmax3=xprobmax3*two*dpi;
         end if

      end select

#:if defined('LOG_RADIUS')
    ! geometry='logSpherical'/'logCylindrical': the radial mesh is uniform in
    ! xi = ln(r), not in r. The par file still states physical radii - the same
    ! convenience the angular bounds above get, where the user writes turns and
    ! the code stores radians - so convert them here, before dx_vec and the
    ! level-1 block extents are derived from them further down. Everything
    ! downstream (ng1, dg1, rnode, the tree, prolongation) then works in the
    ! uniform logical coordinate and needs no further change; only
    ! fill_geometry_device, and the handful of sites that reconstruct a face
    ! position, ever take the exponential back.
    !
    ! Dimension 1 is never the phi_ direction in either system, so this cannot
    ! collide with the 2*pi conversion above.
    if (log_r0 < zero) call mpistop("log_r0 is the offset r0 of the radial map &
       &s = ln(1 + r/r0) and cannot be negative")
    if (log_r0 == zero) then
      if (xprobmin1 <= zero) call mpistop("geometry='logSpherical'/&
         &'logCylindrical' with the default log_r0 = 0 need xprobmin1 > 0: the &
         &radial coordinate is then stretched as ln(r), which has no value at &
         &r = 0. Set log_r0 > 0 to reach the axis")
    else
      if (xprobmin1 < zero) call mpistop("xprobmin1 is a physical radius and &
         &cannot be negative")
    end if
    if (xprobmax1 <= xprobmin1) call mpistop("xprobmax1 must exceed xprobmin1")

    ! The same map in the form device code evaluates it, valid wherever
    ! s >= 0 - which is everywhere in the physical domain. See r_of_s.
    if (log_r0 > zero) then
      log_ra =  log_r0
      log_rb = -log_r0
    else
      log_ra = one
      log_rb = zero
    end if

    xprobmin1 = s_of_r(xprobmin1)
    xprobmax1 = s_of_r(xprobmax1)
    ${GPU_UPDATE_DEVICE('log_r0, log_ra, log_rb')}$

#:else
    if (log_r0 /= zero) call mpistop("log_r0 offsets the logarithmic radial &
       &map and is only meaningful for geometry='logSpherical' or &
       &geometry='logCylindrical'")
#:endif

    ! full block size including ghostcells
    ixGhi1 = block_nx1 + 2*nghostcells
    ixGhi2 = block_nx2 + 2*nghostcells
    ixGhi3 = block_nx3 + 2*nghostcells
    ixGshi1 = ixGhi1
    ixGshi2 = ixGhi2
    ixGshi3 = ixGhi3

    nx_vec = [domain_nx1, domain_nx2, domain_nx3]
    block_nx_vec = [block_nx1, block_nx2, block_nx3]

    if (any(nx_vec < 4) .or. any(mod(nx_vec,&
        2) == 1)) call mpistop&
       ('Grid size (domain_nx^D) has to be even and >= 4')

    if (any(block_nx_vec < 4) .or. any(mod(block_nx_vec,&
        2) == 1)) call mpistop&
       ('Block size (block_nx^D) has to be even and >= 4')

     if(mod(domain_nx1,block_nx1)/=0) call &
        mpistop('Grid (domain_nx^D) and block (block_nx^D) must be consistent')

      if(mod(domain_nx2,block_nx2)/=0) call &
         mpistop&
         ('Grid (domain_nx^D) and block (block_nx^D) must be consistent')
      if(mod(domain_nx3,block_nx3)/=0) call &
         mpistop&
         ('Grid (domain_nx^D) and block (block_nx^D) must be consistent')

    if(refine_max_level>nlevelshi.or.refine_max_level<1)then
       write(unitterm,*)'Error: refine_max_level',refine_max_level,&
          '>nlevelshi ',nlevelshi
       call mpistop("Reset nlevelshi and recompile!")
    endif

    dx_vec = [xprobmax1-xprobmin1, xprobmax2-xprobmin2,&
        xprobmax3-xprobmin3] / nx_vec

    if (mype==0) then
       write(c_ndim, '(I1)') 3
       write(unitterm, '(A30,' // c_ndim // '(I0," "))') &
          ' Domain size (cells): ', nx_vec
       write(unitterm, '(A30,' // c_ndim // '(E9.3," "))') ' Level one dx: ',&
           dx_vec
    end if

    if (any(dx_vec < smalldouble)) call mpistop(&
       "Incorrect domain size (too small grid spacing)")

    dx(:, 1) = dx_vec

    if(sum(w_refine_weight(:))==0) w_refine_weight(1) = 1.d0
    if(dabs(sum(w_refine_weight(:))-1.d0)>smalldouble) then
      write(unitterm,*) "Sum of all elements in w_refine_weight be 1.d0"
      call mpistop("Reset w_refine_weight so the sum is 1.d0")
    end if

    if (mype==0) write(unitterm, '(A30)', advance='no') 'Refine estimation: '

    select case (refine_criterion)
    case (0)
       if (mype==0) write(unitterm, '(A)') "user defined"
    case (1)
       if (mype==0) write(unitterm, '(A)') "relative error"
    case (2)
       if (mype==0) write(unitterm, '(A)') "Lohner's original scheme"
    case (3)
       if (mype==0) write(unitterm, '(A)') "Lohner's scheme"
    case default
       call mpistop("Unknown error estimator, change refine_criterion")
    end select

    if (tfixgrid<bigdouble/2.0d0) then
       if(mype==0)print*,'Warning, at time=',tfixgrid,'the grid will be fixed'
    end if
    if (itfixgrid<biginteger/2) then
       if(mype==0)print*,'Warning, at iteration=',itfixgrid,&
          'the grid will be fixed'
    end if
    if (ditregrid>1) then
       if(mype==0)print*,'Note, Grid is reconstructed once every',ditregrid,&
          'iterations'
    end if


    do islice=1,nslices
       select case(slicedir(islice))
          case(1)
          ! slicecoord is a physical coordinate; xprobmin1/xprobmax1 hold the
          ! logical one, which differs from it under LOG_RADIUS.
#:if defined('LOG_RADIUS')
          if(slicecoord(islice)<zero .or. &
             (log_r0<=zero .and. slicecoord(islice)<=zero)) call mpistop("a &
             &radial slicecoord is a physical radius, and must be > 0 unless &
             &log_r0 > 0 puts r = 0 on the grid")
          if(s_of_r(slicecoord(islice))<xprobmin1 .or. &
             s_of_r(slicecoord(islice))>xprobmax1) &
             write(uniterr,*)'Warning in read_par_files: ', 'Slice ', islice,&
              ' coordinate',slicecoord(islice),'out of bounds for dimension ',&
             slicedir(islice)
#:else
          if(slicecoord(islice)<xprobmin1.or.slicecoord(islice)>xprobmax1) &
             write(uniterr,*)'Warning in read_par_files: ', 'Slice ', islice,&
              ' coordinate',slicecoord(islice),'out of bounds for dimension ',&
             slicedir(islice)
#:endif

          case(2)
          if(slicecoord(islice)<xprobmin2.or.slicecoord(islice)>xprobmax2) &
             write(uniterr,*)'Warning in read_par_files: ', 'Slice ', islice,&
              ' coordinate',slicecoord(islice),'out of bounds for dimension ',&
             slicedir(islice)

          case(3)
          if(slicecoord(islice)<xprobmin3.or.slicecoord(islice)>xprobmax3) &
             write(uniterr,*)'Warning in read_par_files: ', 'Slice ', islice,&
              ' coordinate',slicecoord(islice),'out of bounds for dimension ',&
             slicedir(islice)

       end select
    end do

    if (mype==0) then
       write(unitterm, '(A30,A,A,A,I0,A,I0,A)') 'time_integrator: ', ' ',&
           trim(time_integrator), ' (', nstep, ' stages, ', nlevel,&
           ' time levels)'
       write(unitterm, '(A30,A,A)') 'restart_from_file: ', ' ',&
           trim(restart_from_file)
       write(unitterm, '(A30,L1)') 'converting: ', convert
       write(unitterm, '(A)') ''
    endif

    deallocate(flux_scheme)

    ${GPU_UPDATE_DEVICE('ixGhi1,ixGhi2,ixGhi3,ixGshi1,ixGshi2,ixGshi3,schmid_rad1,schmid_rad2,schmid_rad3,cada3_radius')}$
    ${GPU_UPDATE_DEVICE('H_correction,type_limiter, max_blocks')}$
    ${GPU_UPDATE_DEVICE('rk_beta11,rk_beta22,rk_beta33,rk_beta44,rk_c2,rk_c3,rk_c4')}$
    ${GPU_UPDATE_DEVICE('rk_alfa21,rk_alfa22,rk_alfa31,rk_alfa33,rk_alfa41,rk_alfa44')}$
    ${GPU_UPDATE_DEVICE('rk_beta54,rk_beta55,rk_alfa53,rk_alfa54,rk_alfa55,rk_c5')}$
    ${GPU_UPDATE_DEVICE('typeboundary, specialboundary, refine_max_level')}$
    ${GPU_UPDATE_DEVICE('slab, slab_uniform')}$
    ${GPU_UPDATE_DEVICE('w_refine_weight, amr_wavefilter')}$
    ${GPU_UPDATE_DEVICE('refine_threshold, derefine_ratio')}$
    ${GPU_UPDATE_DEVICE('block_nx1, block_nx2, block_nx3')}$
    ${GPU_UPDATE_DEVICE('courantpar, dtdiffpar')}$
    ${GPU_UPDATE_DEVICE('flux_adaptive_diffusion, flux_ad_min, flux_ad_scale')}$

  end subroutine read_par_files

  !> Routine to find entries in a string
  subroutine get_fields_string(line, delims, n_max, fields, n_found,&
      fully_read)
    !> The line from which we want to read
    character(len=*), intent(in)    :: line
    !> A string with delimiters. For example delims = " ,'"""//char(9)
    character(len=*), intent(in)    :: delims
    !> Maximum number of entries to read in
    integer, intent(in)             :: n_max
    !> Number of entries found
    integer, intent(inout)          :: n_found
    !> Fields in the strings
    character(len=*), intent(inout) :: fields(n_max)
    logical, intent(out), optional  :: fully_read

    integer :: ixs_start(n_max)
    integer :: ixs_end(n_max)
    integer :: ix, ix_prev

    ix_prev = 0
    n_found = 0

    do while (n_found < n_max)
       ! Find the starting point of the next entry (a non-delimiter value)
       ix = verify(line(ix_prev+1:), delims)
       if (ix == 0) exit

       n_found            = n_found + 1
       ixs_start(n_found) = ix_prev + ix !This is the absolute position in 'line'

       ! Get the end point of the current entry (next delimiter index minus one)
       ix = scan(line(ixs_start(n_found)+1:), delims) - 1

       if (ix == -1) then              ! If there is no last delimiter,
          ixs_end(n_found) = len(line) ! the end of the line is the endpoint
       else
          ixs_end(n_found) = ixs_start(n_found) + ix
       end if

       fields(n_found) = line(ixs_start(n_found):ixs_end(n_found))
       ix_prev = ixs_end(n_found) ! We continue to search from here
    end do

    if (present(fully_read)) then
      ix = verify(line(ix_prev+1:), delims)
      fully_read = (ix == 0)    ! Are there only delimiters?
    end if

  end subroutine get_fields_string

  subroutine saveamrfile(ifile)

    use mod_usr_methods, only: usr_print_log, usr_write_analysis
    use mod_global_parameters
    use mod_particles, only: write_particles_snapshot
    use mod_slice, only: write_slice
    use mod_collapse, only: write_collapsed
    use mod_convert_files, only: generate_plotfile
    use mod_physics, only: phys_req_diagonal
    use mod_ghostcells_update
    use mod_geometry, only: sync_geometry_host
    integer:: ifile, iigrid, igrid

    ! update (diagonal) ghostcells before moving to host:
    if(.not. phys_req_diagonal) then
      call getbc(global_time,0.d0,ps,iwstart,nwgc)
    end if

    ${GPU_UPDATE_HOST('bg(1)%w')}$

    ! Same idea for the cell metrics: everything below here that weights by
    ! dvolume, or builds corner positions from dx, reads them on the host.
    call sync_geometry_host()

    select case (ifile)
    case (fileout_)
       ! Write .dat snapshot
       call write_snapshot()

       ! Generate formatted output (e.g., VTK)
       if(autoconvert) call generate_plotfile

       if(use_particles) call write_particles_snapshot()

       snapshotnext = snapshotnext + 1
    case (fileslice_)
       call write_slice
    case (filecollapse_)
       call write_collapsed
    case (filelog_)
       select case (typefilelog)
       case ('default')
          call printlog_default
       case ('regression_test')
          call printlog_regression_test()
       case ('special')
          if (.not. associated(usr_print_log)) then
             call mpistop("usr_print_log not defined")
          else
             call usr_print_log()
          end if
       case default
          call mpistop("Error in SaveFile: Unknown typefilelog")
       end select
    case (fileanalysis_)
       if (associated(usr_write_analysis)) then
          call usr_write_analysis()
       end if
    case default
       write(*,*) 'No save method is defined for ifile=',ifile
       call mpistop("")
    end select

    ! opedit: Flush stdout and stderr from time to time.
    flush(unit=unitterm)

  end subroutine saveamrfile


  ! Check if a snapshot exists
  logical function snapshot_exists(ix)
    use mod_global_parameters
    integer, intent(in)    :: ix !< Index of snapshot
    character(len=std_len) :: filename

    write(filename, "(a,i4.4,a)") trim(base_filename), ix, ".dat"
    inquire(file=trim(filename), exist=snapshot_exists)
  end function snapshot_exists

  integer function get_snapshot_index(filename)
    character(len=*), intent(in) :: filename
    integer                      :: i

    ! Try to parse index in restart_from_file string (e.g. basename0000.dat)
    i = len_trim(filename) - 7
    read(filename(i:i+3), '(I4)') get_snapshot_index
  end function get_snapshot_index



  !> Write header for a snapshot
  !>
  !> If you edit the header, don't forget to update: snapshot_write_header(),
  !> snapshot_read_header(), doc/fileformat.md, tools/python/dat_reader.py
  subroutine snapshot_write_header(fh, offset_tree, offset_block)
    use mod_forest
    use mod_physics
    use mod_global_parameters
    use mod_input_output_helper, only: snapshot_write_header1
    integer, intent(in)                       :: fh           !< File handle
    integer(kind=MPI_OFFSET_KIND), intent(in) :: offset_tree !< Offset of tree info
    integer(kind=MPI_OFFSET_KIND), intent(in) :: offset_block !< Offset of block data
    call snapshot_write_header1(fh, offset_tree, offset_block, cons_wnames,&
        nw)
  end subroutine snapshot_write_header

  !> Read header for a snapshot
  !>
  !> If you edit the header, don't forget to update: snapshot_write_header(),
  !> snapshot_read_header(), doc/fileformat.md, tools/python/dat_reader.py
  subroutine snapshot_read_header(fh, offset_tree, offset_block)
    use mod_forest
    use mod_global_parameters
    use mod_physics, only: physics_type
    integer, intent(in)                   :: fh           !< File handle
    integer(MPI_OFFSET_KIND), intent(out) :: offset_tree !< Offset of tree info
    integer(MPI_OFFSET_KIND), intent(out) :: offset_block !< Offset of block data
    integer                               :: i, version
    integer                               :: ibuf(ndim), iw
    double precision                      :: rbuf(ndim)
    integer, dimension(MPI_STATUS_SIZE)   :: st
    character(len=name_len), allocatable  :: var_names(:), param_names(:)
    double precision, allocatable         :: params(:)
    character(len=name_len)               :: phys_name, geom_name
    integer                               :: er, n_par, tmp_int
    logical                               :: periodic(ndim)

    ! Version number
    call mpi_file_read_wrapper(fh, version, 1, MPI_INTEGER, st, er)
    if (all(compatible_versions /= version)) then
      call mpistop("Incompatible file version (maybe old format?)")
    end if

    ! offset_tree
    call mpi_file_read_wrapper(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    offset_tree = ibuf(1)

    ! offset_block
    call mpi_file_read_wrapper(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    offset_block = ibuf(1)

    ! nw
    call mpi_file_read_wrapper(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    nw_found=ibuf(1)
    if (nw /= ibuf(1)) then
      write(*,*) "nw=",nw," and nw found in restart file=",ibuf(1)
      write(*,*) "Please be aware of changes in w at restart."
      !call mpistop("currently, changing nw at restart is not allowed")
    end if

    ! ndir
    call mpi_file_read_wrapper(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    if (ibuf(1) /= ndir) then
      write(*,*) "ndir in restart file = ",ibuf(1)
      write(*,*) "ndir = ",ndir
      call mpistop("reset ndir to ndir in restart file")
    end if

    ! ndim
    call mpi_file_read_wrapper(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    if (ibuf(1) /= ndim) then
      write(*,*) "ndim in restart file = ",ibuf(1)
      write(*,*) "ndim = ",ndim
      call mpistop("reset ndim to ndim in restart file")
    end if

    ! levmax
    call mpi_file_read_wrapper(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    if (ibuf(1) > refine_max_level) then
      write(*,*) "number of levels in restart file = ",ibuf(1)
      write(*,*) "refine_max_level = ",refine_max_level
      call mpistop("refine_max_level < num. levels in restart file")
    end if

    ! nleafs
    call mpi_file_read_wrapper(fh, nleafs, 1, MPI_INTEGER, st, er)

    ! nparents
    call mpi_file_read_wrapper(fh, nparents, 1, MPI_INTEGER, st, er)

    ! it
    call mpi_file_read_wrapper(fh, it, 1, MPI_INTEGER, st, er)

    ! global time
    call mpi_file_read_wrapper(fh, global_time, 1, MPI_DOUBLE_PRECISION, st, er)

    ! xprobmin^D
    call mpi_file_read_wrapper(fh,rbuf(1:ndim),ndim,MPI_DOUBLE_PRECISION,st,er)
    if (maxval(abs(rbuf(1:ndim) - [ xprobmin1,xprobmin2,&
       xprobmin3 ])) > 0) then
      write(*,*) "Error: xprobmin differs from restart data: ", rbuf(1:ndim)
      call mpistop("change xprobmin^D in par file")
    end if

    ! xprobmax^D
    call mpi_file_read_wrapper(fh,rbuf(1:ndim),ndim,MPI_DOUBLE_PRECISION,st,er)
    if (maxval(abs(rbuf(1:ndim) - [ xprobmax1,xprobmax2,&
       xprobmax3 ])) > 0) then
      write(*,*) "Error: xprobmax differs from restart data: ", rbuf(1:ndim)
      call mpistop("change xprobmax^D in par file")
    end if

    ! domain_nx^D
    call mpi_file_read_wrapper(fh,ibuf(1:ndim), ndim, MPI_INTEGER,st,er)
    if (any(ibuf(1:ndim) /= [ domain_nx1,domain_nx2,domain_nx3 ])) then
      write(*,*) "Error: mesh size differs from restart data: ", ibuf(1:ndim)
      call mpistop("change domain_nx^D in par file")
    end if

    ! block_nx^D
    call mpi_file_read_wrapper(fh,ibuf(1:ndim), ndim, MPI_INTEGER,st,er)
    if (any(ibuf(1:ndim) /= [ block_nx1,block_nx2,block_nx3 ])) then
      write(*,*) "Error: block size differs from restart data:", ibuf(1:ndim)
      call mpistop("change block_nx^D in par file")
    end if

    ! From version 5, read more info about the grid
    if (version > 4) then
      call mpi_file_read_wrapper(fh, periodic, ndim, MPI_LOGICAL, st, er)
      if (periodic(1) .and. .not.periodB(1) .or. .not.periodic(1) .and. &
         periodB(1) .or. periodic(2) .and. .not.periodB(2) .or. &
         .not.periodic(2) .and. periodB(2) .or. periodic(3) .and. &
         .not.periodB(3) .or. .not.periodic(3) .and. periodB(3)) call &
         mpistop("change in periodicity in par file")

      call mpi_file_read_wrapper(fh, geom_name, name_len, MPI_CHARACTER, st, er)

      if (geom_name /= geometry_name(1:name_len)) then
        write(*,*) "type of coordinates in data is: ", geom_name
        call mpistop("select the correct coordinates in mod_usr.t file")
      end if

      call mpi_file_read_wrapper(fh, stagger_mark_dat, 1, MPI_LOGICAL, st, er)
      if (stagger_grid .and. .not. stagger_mark_dat .or. &
         .not.stagger_grid.and.stagger_mark_dat) then
        write(*,*) "Warning: stagger grid flag differs from restart data:",&
            stagger_mark_dat
        !call mpistop("change parameter to use stagger grid")
      end if
    end if

    ! From version 4 onwards, the later parts of the header must be present
    if (version > 3) then
      ! w_names (not used here)
      allocate(var_names(nw_found))
      do iw = 1, nw_found
        call mpi_file_read_wrapper(fh, var_names(iw), name_len, MPI_CHARACTER, st, er)
      end do

      ! Physics related information
      call mpi_file_read_wrapper(fh, phys_name, name_len, MPI_CHARACTER, st, er)

      if (phys_name /= physics_type) then
!        call mpistop("Cannot restart with a different physics type")
      end if

      call mpi_file_read_wrapper(fh, n_par, 1, MPI_INTEGER, st, er)
      allocate(params(n_par))
      allocate(param_names(n_par))
      call mpi_file_read_wrapper(fh, params, n_par, MPI_DOUBLE_PRECISION, st, er)
      call mpi_file_read_wrapper(fh, param_names, name_len * n_par, MPI_CHARACTER, st,&
          er)

      ! Read snapshotnext etc. for restarting
      call mpi_file_read_wrapper(fh, tmp_int, 1, MPI_INTEGER, st, er)

      ! Only set snapshotnext if the user hasn't specified it
      if (snapshotnext == -1) snapshotnext = tmp_int

      call mpi_file_read_wrapper(fh, tmp_int, 1, MPI_INTEGER, st, er)
      if (slicenext == -1) slicenext = tmp_int

      call mpi_file_read_wrapper(fh, tmp_int, 1, MPI_INTEGER, st, er)
      if (collapsenext == -1) collapsenext = tmp_int
    else
      ! Guess snapshotnext from file name if not set
      if (snapshotnext == -1) snapshotnext = &
         get_snapshot_index(trim(restart_from_file)) + 1
      ! Set slicenext and collapsenext if not set
      if (slicenext == -1) slicenext = 0
      if (collapsenext == -1) collapsenext = 0
    end if

    ! Still used in convert
    snapshotini = snapshotnext-1

  end subroutine snapshot_read_header



  subroutine write_snapshot
    use mod_forest
    use mod_global_parameters
    use mod_physics
    use mod_input_output_helper, only: count_ix,block_shape_io,&
       create_output_file
    use mod_functions_forest, only: write_forest

    integer                       :: file_handle, igrid, Morton_no, iwrite
    integer                       :: ipe, ix_buffer(2*ndim+1), n_values
    integer                       :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3, n_ghost(2*ndim)
    integer                       :: ixOsmin1,ixOsmin2,ixOsmin3,ixOsmax1,&
       ixOsmax2,ixOsmax3,n_values_stagger
    integer                       :: iorecvstatus(MPI_STATUS_SIZE)
    integer                       :: ioastatus(MPI_STATUS_SIZE)
    integer                       :: igrecvstatus(MPI_STATUS_SIZE)
    integer                       :: istatus(MPI_STATUS_SIZE)
    type(tree_node), pointer      :: pnode
    integer(kind=MPI_OFFSET_KIND) :: offset_tree_info
    integer(kind=MPI_OFFSET_KIND) :: offset_block_data
    integer(kind=MPI_OFFSET_KIND) :: offset_offsets
    double precision, allocatable :: w_buffer(:)

    integer, allocatable                       :: block_ig(:, :)
    integer, allocatable                       :: block_lvl(:)
    integer(kind=MPI_OFFSET_KIND), allocatable :: block_offset(:)

    call MPI_BARRIER(icomm, ierrmpi)

    ! Allocate send/receive buffer
    n_values = count_ix(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3) * nw
    if(stagger_grid) then
      n_values = n_values + count_ix(ixGslo1,ixGslo2,ixGslo3,ixGshi1,ixGshi2,&
         ixGshi3) * nws
    end if
    allocate(w_buffer(n_values))

    ! Allocate arrays with information about grid blocks
    allocate(block_ig(ndim, nleafs))
    allocate(block_lvl(nleafs))
    allocate(block_offset(nleafs+1))

    ! master processor
    if (mype==0) then
      call create_output_file(file_handle, snapshotnext, ".dat")

      ! Don't know offsets yet, we will write header again later
      offset_tree_info = -1
      offset_block_data = -1
      call snapshot_write_header(file_handle, offset_tree_info,&
          offset_block_data)

      call MPI_File_get_position(file_handle, offset_tree_info, ierrmpi)

      call write_forest(file_handle)

      ! Collect information about the spatial index (ig^D) and refinement level
      ! of leaves
      do Morton_no = Morton_start(0), Morton_stop(npe-1)
         igrid = sfc(1, Morton_no)
         ipe = sfc(2, Morton_no)
         pnode => igrid_to_node(igrid, ipe)%node

         block_ig(:, Morton_no) = [ pnode%ig1,pnode%ig2,pnode%ig3 ]
         block_lvl(Morton_no) = pnode%level
         block_offset(Morton_no) = 0 ! Will be determined later
      end do

      call MPI_FILE_WRITE(file_handle, block_lvl, size(block_lvl), MPI_INTEGER,&
          istatus, ierrmpi)

      call MPI_FILE_WRITE(file_handle, block_ig, size(block_ig), MPI_INTEGER,&
          istatus, ierrmpi)

      ! Block offsets are currently unknown, but will be overwritten later
      call MPI_File_get_position(file_handle, offset_offsets, ierrmpi)
      call mpi_file_write_wrapper(file_handle, block_offset(1:nleafs), nleafs,&
          MPI_OFFSET, istatus, ierrmpi)

      call MPI_File_get_position(file_handle, offset_block_data, ierrmpi)

      ! Check whether data was written as expected
      if (offset_block_data - offset_tree_info /= (nleafs + nparents) * &
         size_logical + nleafs * ((1+ndim) * size_int + 2 * size_int)) then
        if (mype == 0) then
          print *, "Warning: MPI_OFFSET type /= 8 bytes"
          print *, "This *could* cause problems when reading .dat files"
        end if
      end if

      block_offset(1) = offset_block_data
      iwrite = 0
    end if

    do Morton_no=Morton_start(mype), Morton_stop(mype)
      igrid  = sfc_to_igrid(Morton_no)
      itag   = Morton_no

      call block_shape_io(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
         ixOmax2,ixOmax3, n_values)
      if(stagger_grid) then
        w_buffer(1:n_values) = pack(ps(igrid)%w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3, 1:nw), .true.)
        ixOsmin1 = ixOmin1 -1
        ixOsmin2 = ixOmin2 -1
        ixOsmin3 = ixOmin3 -1
        ixOsmax1 = ixOmax1
        ixOsmax2 = ixOmax2
        ixOsmax3 = ixOmax3
        n_values_stagger= count_ix(ixOsmin1,ixOsmin2,ixOsmin3,ixOsmax1,&
           ixOsmax2,ixOsmax3)*nws
        w_buffer(n_values+1:n_values+n_values_stagger) = &
           pack(ps(igrid)%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
           ixOsmin3:ixOsmax3, 1:nws), .true.)
        n_values=n_values+n_values_stagger
      else
        w_buffer(1:n_values) = pack(ps(igrid)%w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3, 1:nw), .true.)
      end if
      ix_buffer(1) = n_values
      ix_buffer(2:) = n_ghost

      if (mype /= 0) then
        call mpi_send_wrapper(ix_buffer, 2*ndim+1, MPI_INTEGER, 0, itag, icomm,&
            ierrmpi)
        call mpi_send_wrapper(w_buffer, n_values, MPI_DOUBLE_PRECISION, 0, itag, icomm,&
            ierrmpi)
      else
        iwrite = iwrite+1
        call mpi_file_write_wrapper(file_handle, ix_buffer(2:), 2*ndim, MPI_INTEGER,&
            istatus, ierrmpi)
        call mpi_file_write_wrapper(file_handle, w_buffer, n_values,&
            MPI_DOUBLE_PRECISION, istatus, ierrmpi)

        ! Set offset of next block
        block_offset(iwrite+1) = block_offset(iwrite) + int(n_values,&
            MPI_OFFSET_KIND) * size_double + 2 * ndim * size_int
      end if
    end do

    ! Write data communicated from other processors
    if (mype == 0) then
      do ipe = 1, npe-1
        do Morton_no=Morton_start(ipe), Morton_stop(ipe)
          iwrite=iwrite+1
          itag=Morton_no

          call mpi_recv_wrapper(ix_buffer, 2*ndim+1, MPI_INTEGER, ipe, itag, icomm,&
             igrecvstatus, ierrmpi)
          n_values = ix_buffer(1)

          call mpi_recv_wrapper(w_buffer, n_values, MPI_DOUBLE_PRECISION,ipe, itag,&
              icomm, iorecvstatus, ierrmpi)

          call mpi_file_write_wrapper(file_handle, ix_buffer(2:), 2*ndim, MPI_INTEGER,&
              istatus, ierrmpi)
          call mpi_file_write_wrapper(file_handle, w_buffer, n_values,&
              MPI_DOUBLE_PRECISION, istatus, ierrmpi)

          ! Set offset of next block
          block_offset(iwrite+1) = block_offset(iwrite) + int(n_values,&
              MPI_OFFSET_KIND) * size_double + 2 * ndim * size_int
        end do
      end do

      ! Write block offsets (now we know them)
      call MPI_FILE_SEEK(file_handle, offset_offsets, MPI_SEEK_SET, ierrmpi)
      call mpi_file_write_wrapper(file_handle, block_offset(1:nleafs), nleafs,&
          MPI_OFFSET, istatus, ierrmpi)

      ! Write header again, now with correct offsets
      call MPI_FILE_SEEK(file_handle, 0_MPI_OFFSET_KIND, MPI_SEEK_SET,&
          ierrmpi)
      call snapshot_write_header(file_handle, offset_tree_info,&
          offset_block_data)

      call MPI_FILE_CLOSE(file_handle, ierrmpi)
    end if

    call MPI_BARRIER(icomm, ierrmpi)

  end subroutine write_snapshot

  !> Routine to read in snapshots (.dat files). When it cannot recognize the
  !> file version, it will automatically try the 'old' reader.
  subroutine read_snapshot
    use mod_usr_methods, only: usr_transform_w
    use mod_input_output_helper, only: count_ix
    use mod_forest
    use mod_global_parameters
    use mod_amr_solution_node, only: alloc_node
    use mod_functions_forest, only: read_forest

    integer                       :: ix_buffer(2*ndim+1), n_values,&
        n_values_stagger
    integer                       :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3, ixOsmin1,ixOsmin2,ixOsmin3,ixOsmax1,ixOsmax2,ixOsmax3
    integer                       :: file_handle, amode, igrid, Morton_no,&
        iread
    integer                       :: istatus(MPI_STATUS_SIZE)
    integer                       :: iorecvstatus(MPI_STATUS_SIZE)
    integer                       :: ipe,inrecv,nrecv, file_version
    logical                       :: fexist
    integer(MPI_OFFSET_KIND)      :: offset_tree_info
    integer(MPI_OFFSET_KIND)      :: offset_block_data
    double precision, allocatable :: w_buffer(:)
    double precision, dimension(:,:,:,:), allocatable :: w
    double precision :: ws(ixGslo1:ixGshi1,ixGslo2:ixGshi2,ixGslo3:ixGshi3,&
       1:ndim)

    if (mype==0) then
      inquire(file=trim(restart_from_file), exist=fexist)
      if (.not.fexist) call mpistop(trim(restart_from_file)//" not found!")

      call MPI_FILE_OPEN(MPI_COMM_SELF,restart_from_file,MPI_MODE_RDONLY,&
          MPI_INFO_NULL,file_handle,ierrmpi)
      call mpi_file_read_wrapper(file_handle, file_version, 1, MPI_INTEGER, istatus,&
          ierrmpi)
    end if

    call MPI_BCAST(file_version,1,MPI_INTEGER,0,icomm,ierrmpi)

    if (all(compatible_versions /= file_version)) then
      if (mype == 0) print *, "Unknown version, trying old snapshot reader..."
      call MPI_FILE_CLOSE(file_handle,ierrmpi)
      call read_snapshot_old()

      ! Guess snapshotnext from file name if not set
      if (snapshotnext == -1) snapshotnext = &
         get_snapshot_index(trim(restart_from_file)) + 1
      ! Set slicenext and collapsenext if not set
      if (slicenext == -1) slicenext = 0
      if (collapsenext == -1) collapsenext = 0

      ! Still used in convert
      snapshotini = snapshotnext-1

      return ! Leave this routine
    else if (mype == 0) then
      call MPI_FILE_SEEK(file_handle, 0_MPI_OFFSET_KIND, MPI_SEEK_SET,&
          ierrmpi)
      call snapshot_read_header(file_handle, offset_tree_info,&
          offset_block_data)
    end if

    ! Share information about restart file
    call MPI_BCAST(nw_found,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(nleafs,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(nparents,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(it,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(global_time,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)

    call MPI_BCAST(snapshotnext,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(slicenext,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(collapsenext,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(stagger_mark_dat,1,MPI_LOGICAL,0,icomm,ierrmpi)

    ! Allocate send/receive buffer
    n_values = count_ix(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3) * nw_found
    if(stagger_mark_dat) then
      n_values = n_values + count_ix(ixGslo1,ixGslo2,ixGslo3,ixGshi1,ixGshi2,&
         ixGshi3) * nws
    end if
    allocate(w_buffer(n_values))
    allocate(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw_found))

    nleafs_active = nleafs

    if (mype == 0) then
      call MPI_FILE_SEEK(file_handle, offset_tree_info, MPI_SEEK_SET, ierrmpi)
    end if

    call read_forest(file_handle)

    do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid=sfc_to_igrid(Morton_no)
      call alloc_node(igrid)
    end do

    if (mype==0) then
      call MPI_FILE_SEEK(file_handle, offset_block_data, MPI_SEEK_SET,&
          ierrmpi)

      iread = 0
      do ipe = 0, npe-1
        do Morton_no=Morton_start(ipe),Morton_stop(ipe)
          iread=iread+1
          itag=Morton_no

          call mpi_file_read_wrapper(file_handle,ix_buffer(1:2*ndim), 2*ndim,&
              MPI_INTEGER, istatus,ierrmpi)

          ! Construct ixO^L array from number of ghost cells
          ixOmin1 = ixMlo1 - ix_buffer(1)
          ixOmin2 = ixMlo2 - ix_buffer(2)
          ixOmin3 = ixMlo3 - ix_buffer(3)
          ixOmax1 = ixMhi1 + ix_buffer(ndim+1)
          ixOmax2 = ixMhi2 + ix_buffer(ndim+2)
          ixOmax3 = ixMhi3 + ix_buffer(ndim+3)
          n_values = count_ix(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
             ixOmax3) * nw_found
          if(stagger_mark_dat) then
            ixOsmin1 = ixOmin1 - 1
            ixOsmin2 = ixOmin2 - 1
            ixOsmin3 = ixOmin3 - 1
            ixOsmax1 = ixOmax1
            ixOsmax2 = ixOmax2
            ixOsmax3 = ixOmax3
            n_values_stagger = n_values
            n_values = n_values + count_ix(ixOsmin1,ixOsmin2,ixOsmin3,ixOsmax1,&
               ixOsmax2,ixOsmax3) * nws
          end if

          call mpi_file_read_wrapper(file_handle, w_buffer, n_values,&
              MPI_DOUBLE_PRECISION, istatus, ierrmpi)

          if (mype == ipe) then ! Root task
            igrid=sfc_to_igrid(Morton_no)
            block=>ps(igrid)
            if(stagger_mark_dat) then
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  1:nw_found) = reshape(w_buffer(1:n_values_stagger),&
                  shape(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  1:nw_found)))
              if(stagger_grid) ps(igrid)%ws(ixOsmin1:ixOsmax1,&
                 ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
                 1:nws)=reshape(w_buffer(n_values_stagger+1:n_values),&
                  shape(ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 ixOsmin3:ixOsmax3, 1:nws)))
            else
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  1:nw_found) = reshape(w_buffer(1:n_values),&
                  shape(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  1:nw_found)))
            end if
            if (nw_found<nw) then
              if (associated(usr_transform_w)) then
                call usr_transform_w(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
                   ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,nw_found,w,&
                   ps(igrid)%x,ps(igrid)%w)
              else
                ps(igrid)%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   1:nw_found)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,1:nw_found)
              end if
            else if (nw_found>nw) then
              if (associated(usr_transform_w)) then
                call usr_transform_w(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
                   ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,nw_found,w,&
                   ps(igrid)%x,ps(igrid)%w)
              else
                ps(igrid)%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   1:nw)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   1:nw)
              end if
            else
              ps(igrid)%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 1:nw)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw)
            end if
#:if defined('FILL_NWEXTRA_ANALYTIC')
            ${GPU_UPDATE_DEVICE('bg(1)%w(:,:,:,1:nwgc,igrid)')}$
#:else
            ${GPU_UPDATE_DEVICE('bg(1)%w(:,:,:,:,igrid)')}$
#:endif
          else
            call mpi_send_wrapper([ ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
                n_values ], 2*ndim+1, MPI_INTEGER, ipe, itag, icomm, ierrmpi)
            call mpi_send_wrapper(w_buffer, n_values, MPI_DOUBLE_PRECISION, ipe, itag,&
                icomm, ierrmpi)
          end if
        end do
      end do

      call MPI_FILE_CLOSE(file_handle,ierrmpi)

    else ! mype > 0

      do Morton_no=Morton_start(mype),Morton_stop(mype)
        igrid=sfc_to_igrid(Morton_no)
        block=>ps(igrid)
        itag=Morton_no

        call mpi_recv_wrapper(ix_buffer, 2*ndim+1, MPI_INTEGER, 0, itag, icomm,&
           iorecvstatus, ierrmpi)
        ixOmin1 = ix_buffer(1)
        ixOmin2 = ix_buffer(2)
        ixOmin3 = ix_buffer(3)
        ixOmax1 = ix_buffer(ndim+1)
        ixOmax2 = ix_buffer(ndim+2)
        ixOmax3 = ix_buffer(ndim+3)
        n_values = ix_buffer(2*ndim+1)

        call mpi_recv_wrapper(w_buffer, n_values, MPI_DOUBLE_PRECISION,0, itag, icomm,&
            iorecvstatus, ierrmpi)

        if(stagger_mark_dat) then
          n_values_stagger = count_ix(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
             ixOmax3) * nw_found
          ixOsmin1 = ixOmin1 - 1
          ixOsmin2 = ixOmin2 - 1
          ixOsmin3 = ixOmin3 - 1
          ixOsmax1 = ixOmax1
          ixOsmax2 = ixOmax2
          ixOsmax3 = ixOmax3
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              1:nw_found) = reshape(w_buffer(1:n_values_stagger),&
              shape(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              1:nw_found)))
          if(stagger_grid) ps(igrid)%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
             ixOsmin3:ixOsmax3,1:nws)=reshape(w_buffer(n_values_stagger+&
             1:n_values), shape(ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
             ixOsmin3:ixOsmax3, 1:nws)))
        else
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              1:nw_found) = reshape(w_buffer(1:n_values),&
              shape(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              1:nw_found)))
        end if
        if (nw_found<nw) then
          if (associated(usr_transform_w)) then
            call usr_transform_w(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
               ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,nw_found,w,&
               ps(igrid)%x,ps(igrid)%w)
          else
            ps(igrid)%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               1:nw_found)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               1:nw_found)
          end if
        else if (nw_found>nw) then
          if (associated(usr_transform_w)) then
            call usr_transform_w(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
               ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,nw_found,w,&
               ps(igrid)%x,ps(igrid)%w)
          else
            ps(igrid)%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               1:nw)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw)
          end if
        else
          ps(igrid)%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1:nw)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw)
        end if
#:if defined('FILL_NWEXTRA_ANALYTIC')
        ${GPU_UPDATE_DEVICE('bg(1)%w(:,:,:,1:nwgc,igrid)')}$
#:else
        ${GPU_UPDATE_DEVICE('bg(1)%w(:,:,:,:,igrid)')}$
#:endif
      end do
    end if

    call MPI_BARRIER(icomm,ierrmpi)

  end subroutine read_snapshot

  subroutine read_snapshot_old()
    use mod_forest
    use mod_global_parameters
    use mod_amr_solution_node, only: alloc_node
    use mod_functions_forest, only: read_forest

    double precision              :: wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:nw)
    integer                       :: fh, igrid, Morton_no, iread
    integer                       :: levmaxini, ndimini, ndirini
    integer                       :: nwini, neqparini, nxini1,nxini2,nxini3
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer                       :: istatus(MPI_STATUS_SIZE)
    integer, allocatable          :: iorecvstatus(:,:)
    integer                       :: ipe,inrecv,nrecv
    integer                       :: sendini(7+3)
    character(len=80)             :: filename
    logical                       :: fexist
    double precision              :: eqpar_dummy(100)

    if (mype==0) then
      call MPI_FILE_OPEN(MPI_COMM_SELF,trim(restart_from_file),&
          MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierrmpi)

      offset=-int(7*size_int+size_double,kind=MPI_OFFSET_KIND)
      call MPI_FILE_SEEK(fh,offset,MPI_SEEK_END,ierrmpi)

      call mpi_file_read_wrapper(fh,nleafs,1,MPI_INTEGER,istatus,ierrmpi)
      nleafs_active = nleafs
      call mpi_file_read_wrapper(fh,levmaxini,1,MPI_INTEGER,istatus,ierrmpi)
      call mpi_file_read_wrapper(fh,ndimini,1,MPI_INTEGER,istatus,ierrmpi)
      call mpi_file_read_wrapper(fh,ndirini,1,MPI_INTEGER,istatus,ierrmpi)
      call mpi_file_read_wrapper(fh,nwini,1,MPI_INTEGER,istatus,ierrmpi)
      call mpi_file_read_wrapper(fh,neqparini,1,MPI_INTEGER,istatus,ierrmpi)
      call mpi_file_read_wrapper(fh,it,1,MPI_INTEGER,istatus,ierrmpi)
      call mpi_file_read_wrapper(fh,global_time,1,MPI_DOUBLE_PRECISION,istatus,&
         ierrmpi)

      ! check if settings are suitable for restart
      if (levmaxini>refine_max_level) then
        write(*,*) "number of levels in restart file = ",levmaxini
        write(*,*) "refine_max_level = ",refine_max_level
        call mpistop("refine_max_level < number of levels in restart file")
      end if
      if (ndimini/=ndim) then
        write(*,*) "ndim in restart file = ",ndimini
        write(*,*) "ndim = ",ndim
        call mpistop("reset ndim to ndim in restart file")
      end if
      if (ndirini/=ndir) then
        write(*,*) "ndir in restart file = ",ndirini
        write(*,*) "ndir = ",ndir
        call mpistop("reset ndir to ndir in restart file")
      end if
      if (nw/=nwini) then
        write(*,*) "nw=",nw," and nw in restart file=",nwini
        call mpistop("currently, changing nw at restart is not allowed")
      end if

      offset=offset-int(ndimini*size_int+neqparini*size_double,&
         kind=MPI_OFFSET_KIND)
      call MPI_FILE_SEEK(fh,offset,MPI_SEEK_END,ierrmpi)

      call mpi_file_read_wrapper(fh,nxini1,1,MPI_INTEGER,istatus,ierrmpi)
      call mpi_file_read_wrapper(fh,nxini2,1,MPI_INTEGER,istatus,ierrmpi)
      call mpi_file_read_wrapper(fh,nxini3,1,MPI_INTEGER,istatus,ierrmpi)
      if (ixGhi1/=nxini1+2*nghostcells.or.ixGhi2/=&
         nxini2+2*nghostcells.or.ixGhi3/=nxini3+2*nghostcells) then
        write(*,*) "Error: reset resolution to ",nxini1+2*nghostcells,&
           nxini2+2*nghostcells,nxini3+2*nghostcells
        call mpistop("change with setamrvac")
      end if

      call mpi_file_read_wrapper(fh,eqpar_dummy,neqparini, MPI_DOUBLE_PRECISION,&
         istatus,ierrmpi)
    end if

    ! broadcast the global parameters first
    if (npe>1) then
      if (mype==0) then
        sendini=(/nleafs,levmaxini,ndimini,ndirini,nwini,neqparini,it ,nxini1,&
           nxini2,nxini3 /)
      end if
      call MPI_BCAST(sendini,7+3,MPI_INTEGER,0,icomm,ierrmpi)
      nleafs=sendini(1);levmaxini=sendini(2);ndimini=sendini(3);
      ndirini=sendini(4);nwini=sendini(5);
      neqparini=sendini(6);it=sendini(7);
      nxini1=sendini(7+1);nxini2=sendini(7+2);nxini3=sendini(7+3);
      nleafs_active = nleafs
      call MPI_BCAST(global_time,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    end if

    if (mype == 0) then
      offset = int(size_block_io,kind=MPI_OFFSET_KIND) * int(nleafs,&
         kind=MPI_OFFSET_KIND)
      call MPI_FILE_SEEK(fh,offset,MPI_SEEK_SET,ierrmpi)
    end if

    call read_forest(fh)

    do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid=sfc_to_igrid(Morton_no)
      call alloc_node(igrid)
    end do

    if (mype==0)then
      iread=0

      do Morton_no=Morton_start(0),Morton_stop(0)
        igrid=sfc_to_igrid(Morton_no)
        iread=iread+1
        offset=int(size_block_io,kind=MPI_OFFSET_KIND) *int(Morton_no-1,&
           kind=MPI_OFFSET_KIND)
        call mpi_file_read_at_wrapper(fh,offset,ps(igrid)%w,1,type_block_io, istatus,&
           ierrmpi)
      end do
      if (npe>1) then
        do ipe=1,npe-1
          do Morton_no=Morton_start(ipe),Morton_stop(ipe)
            iread=iread+1
            itag=Morton_no
            offset=int(size_block_io,kind=MPI_OFFSET_KIND)*int(Morton_no-1,&
               kind=MPI_OFFSET_KIND)
            call mpi_file_read_at_wrapper(fh,offset,wio,1,type_block_io,istatus,&
               ierrmpi)
            call mpi_send_wrapper(wio,1,type_block_io,ipe,itag,icomm,ierrmpi)
          end do
        end do
      end if
      call MPI_FILE_CLOSE(fh,ierrmpi)
    else
      nrecv=(Morton_stop(mype)-Morton_start(mype)+1)
      allocate(iorecvstatus(MPI_STATUS_SIZE,nrecv))
      inrecv=0
      do Morton_no=Morton_start(mype),Morton_stop(mype)
        igrid=sfc_to_igrid(Morton_no)
        itag=Morton_no
        inrecv=inrecv+1
        call mpi_recv_wrapper(ps(igrid)%w,1,type_block_io,0,itag,icomm,iorecvstatus(:,&
           inrecv),ierrmpi)
      end do
      deallocate(iorecvstatus)
    end if

    call MPI_BARRIER(icomm,ierrmpi)

  end subroutine read_snapshot_old

  !> Write volume-averaged values and other information to the log file
  subroutine printlog_default

    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters

    logical              :: fileopen
    integer              :: i, iw, level
    double precision     :: wmean(1:nw), total_volume
    double precision     :: volume_coverage(refine_max_level)
    integer              :: nx1,nx2,nx3, nc, ncells, dit
    double precision     :: dtTimeLast, now, cellupdatesPerSecond
    double precision     :: activeBlocksPerCore, wctPerCodeTime, timeToFinish
    character(len=40)    :: fmt_string
    character(len=80)    :: filename
    character(len=2048)  :: line
    logical, save        :: opened  = .false.
    integer              :: amode, istatus(MPI_STATUS_SIZE)
    integer, parameter   :: my_unit = 20

    ! Compute the volume-average of w**1 = w
    call get_volume_average(1, wmean, total_volume)

    ! Compute the volume coverage
    call get_volume_coverage(volume_coverage)

    if (mype == 0) then

       ! To compute cell updates per second, we do the following:
       nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
       nc=nx1*nx2*nx3
       ncells = nc * nleafs_active

       ! assumes the number of active leafs haven't changed since last compute.
       now        = MPI_WTIME()
       dit        = it - itTimeLast
       dtTimeLast = now - timeLast
       itTimeLast = it
       timeLast   = now
       cellupdatesPerSecond = dble(ncells) * dble(nstep) * dble(dit) / &
          (dtTimeLast * dble(npe))

       ! blocks per core:
       activeBlocksPerCore = dble(nleafs_active) / dble(npe)

       ! Wall clock time per code time unit in seconds:
       wctPerCodeTime = dtTimeLast / max(dit * dt, epsilon(1.0d0))

       ! Wall clock time to finish in hours:
       timeToFinish = (time_max - global_time) * wctPerCodeTime / 3600.0d0

       ! On first entry, open the file and generate the header
       if (.not. opened) then

          filename = trim(base_filename) // ".log"

          ! Delete the log when not doing a restart run
          if (restart_from_file == undefined) then
             open(unit=my_unit,file=trim(filename),status='replace')
             close(my_unit, status='delete')
          end if

          amode    = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
          amode    = ior(amode,MPI_MODE_APPEND)

          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, MPI_INFO_NULL,&
              log_fh, ierrmpi)

          opened = .true.

          ! Start of file headern
          line = "it global_time dt"
          do level=1,nw
             i = len_trim(line) + 2
             write(line(i:),"(a,a)") trim(cons_wnames(level)), " "
          end do

          ! Volume coverage per level
          do level = 1, refine_max_level
             i = len_trim(line) + 2
             write(line(i:), "(a,i0)") "c", level
          end do

          ! Cell counts per level
          do level=1,refine_max_level
             i = len_trim(line) + 2
             write(line(i:), "(a,i0)") "n", level
          end do

          ! Rest of file header
          line = trim(line) // " | Xload Xmemory 'Cell_Updates /second/core'"
          line = trim(line) // " 'Active_Blocks/Core' 'Wct Per Code Time [s]'"
          line = trim(line) // " 'TimeToFinish [hrs]'"

          ! Only write header if not restarting
          if (restart_from_file == undefined .or. reset_time) then
            call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'),&
                len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
          end if
       end if

       ! Construct the line to be added to the log

       fmt_string = '(' // fmt_i // ',2' // fmt_r // ')'
       write(line, fmt_string) it, global_time, dt
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', nw, fmt_r // ')'
       write(line(i:), fmt_string) wmean(1:nw)
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', refine_max_level, fmt_r // ')'
       write(line(i:), fmt_string) volume_coverage(1:refine_max_level)
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', refine_max_level, fmt_i // ')'
       write(line(i:), fmt_string) nleafs_level(1:refine_max_level)
       i = len_trim(line) + 2

       fmt_string = '(a,6' // fmt_r2 // ')'
       write(line(i:), fmt_string) '| ', Xload, Xmemory, cellupdatesPerSecond,&
           activeBlocksPerCore, wctPerCodeTime, timeToFinish

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') ,&
           len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

  end subroutine printlog_default

  !> Print a log that can be used to check whether the code still produces the
  !> same output (regression test)
  subroutine printlog_regression_test()
    use mod_global_parameters

    integer, parameter :: n_modes = 2
    character(len=40)  :: fmt_string
    character(len=2048)  :: line
    logical, save      :: file_open = .false.
    integer            :: power
    double precision   :: modes(nw, n_modes), volume
    integer              :: amode, istatus(MPI_STATUS_SIZE)
    character(len=80)    :: filename

    do power = 1, n_modes
       call get_volume_average(power, modes(:, power), volume)
    end do

    if (mype == 0) then
       if (.not. file_open) then
          filename = trim(base_filename) // ".log"
          amode    = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
          amode    = ior(amode,MPI_MODE_APPEND)

          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, MPI_INFO_NULL,&
              log_fh, ierrmpi)
          file_open = .true.

          line= "# time mean(w) mean(w**2)"
          call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'),&
              len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
       end if

       write(fmt_string, "(a,i0,a)") "(", nw * n_modes + 1, fmt_r // ")"
       write(line, fmt_string) global_time, modes
       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') ,&
           len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if
  end subroutine printlog_regression_test

  !> Compute mean(w**power) over the leaves of the grid. The first mode
  !> (power=1) corresponds to the mean, the second to the mean squared values
  !> and so on.
  subroutine get_volume_average(power, mode, volume)
    use mod_global_parameters

    integer, intent(in)           :: power     !< Which mode to compute
    double precision, intent(out) :: mode(nw)  !< The computed mode
    double precision, intent(out) :: volume    !< The total grid volume
    integer                       :: iigrid, igrid, iw
    integer                       :: ncells
    double precision              :: wsum(nw+1)
    double precision              :: dvol
    double precision              :: dsum_recv(1:nw+1)

    wsum(:) = 0
    ncells = (ixMhi1-ixMlo1+1)*(ixMhi2-ixMlo2+1)*(ixMhi3-ixMlo3+1)

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       if (slab_uniform) then
          ! Every cell of the block has the same volume, and ps(igrid)%dvolume
          ! is not even allocated in a Cartesian build - take the one value
          ! from the block's spacing in rnode and factor it out of the sums.
          dvol = rnode(rpdx1_,igrid)*rnode(rpdx2_,igrid)*rnode(rpdx3_,igrid)

          wsum(nw+1) = wsum(nw+1) + dvol*dble(ncells)

          do iw = 1, nw
             wsum(iw) = wsum(iw) + dvol*sum(ps(igrid)%w(ixMlo1:ixMhi1,&
                ixMlo2:ixMhi2,ixMlo3:ixMhi3,iw)**power)
          end do
       else
          ! Read the metric out of bgeo directly rather than through the
          ! ps(igrid)%dvolume pointer view: bgeo is where the data actually
          ! lives, and the view is a bounds-remapped pointer into it.
          ! Store total volume in last element
          wsum(nw+1) = wsum(nw+1) + sum(bgeo%dvolume(ixMlo1:ixMhi1,&
             ixMlo2:ixMhi2,ixMlo3:ixMhi3,igrid))

          ! Compute the modes of the cell-centered variables, weighted by volume
          do iw = 1, nw
             wsum(iw) = wsum(iw) + sum(bgeo%dvolume(ixMlo1:ixMhi1,&
                ixMlo2:ixMhi2,ixMlo3:ixMhi3,igrid)*ps(igrid)%w(ixMlo1:ixMhi1,&
                ixMlo2:ixMhi2,ixMlo3:ixMhi3,iw)**power)
          end do
       end if
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(wsum, dsum_recv, nw+1, MPI_DOUBLE_PRECISION, MPI_SUM,&
        icomm, ierrmpi)

    ! Set the volume and the average
    volume = dsum_recv(nw+1)
    mode   = dsum_recv(1:nw) / volume

  end subroutine get_volume_average

  !> Compute how much of the domain is covered by each grid level. This routine
  !> does not take a non-Cartesian geometry into account.
  subroutine get_volume_coverage(vol_cov)
    use mod_global_parameters

    double precision, intent(out) :: vol_cov(1:refine_max_level)
    double precision              :: dsum_recv(1:refine_max_level)
    integer                       :: iigrid, igrid, iw, level

    ! First determine the total 'flat' volume in each level
    vol_cov(1:refine_max_level)=zero

    do iigrid = 1, igridstail
       igrid          = igrids(iigrid);
       level          = node(plevel_,igrid)
       vol_cov(level) = vol_cov(level)+ (rnode(rpxmax1_,igrid)-rnode(rpxmin1_,&
          igrid))*(rnode(rpxmax2_,igrid)-rnode(rpxmin2_,&
          igrid))*(rnode(rpxmax3_,igrid)-rnode(rpxmin3_,igrid))
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(vol_cov, dsum_recv, refine_max_level,&
        MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierrmpi)

    ! Normalize
    vol_cov = dsum_recv / sum(dsum_recv)
  end subroutine get_volume_coverage

  !> Compute the volume average of func(w) over the leaves of the grid.
  subroutine get_volume_average_func(func, f_avg, volume)
    use mod_global_parameters

    interface
       pure function func(w_vec, w_size) result(val)
         integer, intent(in)          :: w_size
         double precision, intent(in) :: w_vec(w_size)
         double precision             :: val
       end function func
    end interface
    double precision, intent(out) :: f_avg  !< The volume average of func
    double precision, intent(out) :: volume    !< The total grid volume
    integer                       :: iigrid, igrid, i1,i2,i3
    double precision              :: wsum(2)
    double precision              :: dsum_recv(2)

    wsum(:) = 0

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! Store total volume in last element
       wsum(2) = wsum(2) + sum(ps(igrid)%dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
          ixMlo3:ixMhi3))

       ! Compute the modes of the cell-centered variables, weighted by volume
       do i1 = ixMlo1, ixMhi1
       do i2 = ixMlo2, ixMhi2
       do i3 = ixMlo3, ixMhi3
       wsum(1) = wsum(1) + ps(igrid)%dvolume(i1,i2,i3) * func(ps(igrid)%w(i1,&
          i2,i3, :), nw)
       end do
       end do
       end do
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(wsum, dsum_recv, 2, MPI_DOUBLE_PRECISION, MPI_SUM,&
        icomm, ierrmpi)

    ! Set the volume and the average
    volume = dsum_recv(2)
    f_avg  = dsum_recv(1) / volume

  end subroutine get_volume_average_func

  !> Compute global maxima of iw variables over the leaves of the grid.
  subroutine get_global_maxima(wmax)
    use mod_global_parameters

    double precision, intent(out) :: wmax(nw)  !< The global maxima

    integer                       :: iigrid, igrid, iw
    double precision              :: wmax_mype(nw),wmax_recv(nw)

    wmax_mype(1:nw) = -bigdouble

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       do iw = 1, nw
          wmax_mype(iw)=max(wmax_mype(iw),maxval(ps(igrid)%w(ixMlo1:ixMhi1,&
             ixMlo2:ixMhi2,ixMlo3:ixMhi3,iw)))
       end do
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(wmax_mype, wmax_recv, nw, MPI_DOUBLE_PRECISION, MPI_MAX,&
        icomm, ierrmpi)

    wmax(1:nw)=wmax_recv(1:nw)

  end subroutine get_global_maxima

  !> Compute global minima of iw variables over the leaves of the grid.
  subroutine get_global_minima(wmin)
    use mod_global_parameters

    double precision, intent(out) :: wmin(nw)  !< The global maxima

    integer                       :: iigrid, igrid, iw
    double precision              :: wmin_mype(nw),wmin_recv(nw)

    wmin_mype(1:nw) = bigdouble

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       do iw = 1, nw
          wmin_mype(iw)=min(wmin_mype(iw),minval(ps(igrid)%w(ixMlo1:ixMhi1,&
             ixMlo2:ixMhi2,ixMlo3:ixMhi3,iw)))
       end do
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(wmin_mype, wmin_recv, nw, MPI_DOUBLE_PRECISION, MPI_MIN,&
        icomm, ierrmpi)

    wmin(1:nw)=wmin_recv(1:nw)

  end subroutine get_global_minima



end module mod_input_output
