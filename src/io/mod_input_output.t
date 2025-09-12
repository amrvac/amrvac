!> Module for reading input and writing output
module mod_input_output
  use mod_comm_lib, only: mpistop

  implicit none
  public


  !> coefficient for rk2_alfa
  double precision, private :: rk2_alfa

  !> List of compatible versions
  integer, parameter :: compatible_versions(3) = [3, 4, 5]

  !> number of w found in dat files
  integer :: nw_found

  !> tag for MPI message
  integer, private :: itag

  !> whether staggered field is in dat
  logical, private :: stagger_mark_dat=.false.

  ! Formats used in output
  character(len=*), parameter :: fmt_r  = 'es16.8' ! Default precision
  character(len=*), parameter :: fmt_r2 = 'es10.2' ! Two digits
  character(len=*), parameter :: fmt_i  = 'i8'     ! Integer format

  ! public methods
  public :: snapshot_write_header  

contains

  !> Read the command line arguments passed to amrvac
  subroutine read_arguments()
    use mod_global_parameters

    integer                          :: len, stat, n, i, ipars
    integer, parameter               :: max_files = 20 ! Maximum number of par files
    integer                          :: n_par_files
    logical                          :: unknown_arg, help, morepars
    character(len=max_files*std_len) :: all_par_files
    character(len=std_len)           :: tmp_files(max_files), arg

    if (mype == 0) then
       print *, '-----------------------------------------------------------------------------'
       print *, '-----------------------------------------------------------------------------'
       print *, '|         __  __ ____ ___        _    __  __ ______     ___    ____         |'
       print *, '|        |  \/  |  _ \_ _|      / \  |  \/  |  _ \ \   / / \  / ___|        |'
       print *, '|        | |\/| | |_) | |_____ / _ \ | |\/| | |_) \ \ / / _ \| |            |'
       print *, '|        | |  | |  __/| |_____/ ___ \| |  | |  _ < \ V / ___ \ |___         |'
       print *, '|        |_|  |_|_|  |___|   /_/   \_\_|  |_|_| \_\ \_/_/   \_\____|        |'
       print *, '-----------------------------------------------------------------------------'
       print *, '-----------------------------------------------------------------------------'
    end if

    ! =============== Fortran 2003 command line reading ================

    ! Default command line arguments
    all_par_files="amrvac.par"
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
          print *, 'mpirun -np 4 ./amrvac -i file.par [file2.par ...]'
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
    call get_fields_string(all_par_files, " ,'"""//char(9), max_files, &
         tmp_files, n_par_files)

    allocate(par_files(n_par_files))
    par_files = tmp_files(1:n_par_files)

  end subroutine read_arguments

  !> Read in the user-supplied parameter-file
  subroutine read_par_files()
    use mod_global_parameters
    use mod_physics, only: phys_energy, physics_type, phys_wider_stencil
    use mod_small_values
    use mod_limiter
    use mod_slice
    use mod_geometry
    use mod_source
    use mod_input_output_helper, only: get_names_from_string

    double precision, dimension(nsavehi) :: tsave_log, tsave_dat, tsave_slice, &
         tsave_collapsed, tsave_custom
    double precision :: dtsave_log, dtsave_dat, dtsave_slice, &
         dtsave_collapsed, dtsave_custom
    double precision :: tsavestart_log, tsavestart_dat, tsavestart_slice, &
         tsavestart_collapsed, tsavestart_custom
    double precision :: sizeuniformpart^D
    double precision :: im_delta,im_nu,rka54,rka51,rkb54,rka55
    double precision :: dx_vec(^ND)
    integer :: ditsave_log, ditsave_dat, ditsave_slice, &
         ditsave_collapsed, ditsave_custom
    integer :: windex, ipower
    integer :: i, j, k, ifile, io_state
    integer :: iB, isave, iw, level, idim, islice
    integer :: nx_vec(^ND), block_nx_vec(^ND)
    integer :: my_unit, iostate
    integer :: ilev
    logical :: fileopen, file_exists

    character              :: c_ndim
    character(len=80)      :: fmt_string
    character(len=std_len) :: err_msg, restart_from_file_arg
    character(len=std_len) :: basename_full, basename_prev, dummy_file
    character(len=std_len), dimension(:), allocatable :: &
         typeboundary_min^D, typeboundary_max^D
    character(len=std_len), allocatable :: limiter(:)
    character(len=std_len), allocatable :: gradient_limiter(:)
    character(len=name_len) :: stretch_dim(ndim)
    !> How to apply dimensional splitting to the source terms, see
    !> @ref discretization.md
    character(len=std_len) :: typesourcesplit
    !> Which flux scheme of spatial discretization to use (per grid level)
    character(len=std_len), allocatable :: flux_scheme(:)
    !> Which type of the maximal bound speed of Riemann fan to use
    character(len=std_len) :: typeboundspeed
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

 
    namelist /filelist/ base_filename,restart_from_file, &
         typefilelog,firstprocess,reset_grid,snapshotnext, &
         convert,convert_type,saveprim,usr_filename,&
         nwauxio,nocartesian, w_write,writelevel,&
         writespshift,length_convert_factor, w_convert_factor, &
         time_convert_factor,level_io,level_io_min, level_io_max, &
         autoconvert,slice_type,slicenext,collapsenext,collapse_type, &
         type_endian

    namelist /savelist/ tsave,itsave,dtsave,ditsave,nslices,slicedir, &
         slicecoord,collapse,collapseLevel, time_between_print,&
         tsave_log, tsave_dat, tsave_slice, tsave_collapsed, tsave_custom, &
         dtsave_log, dtsave_dat, dtsave_slice, dtsave_collapsed, dtsave_custom, &
         ditsave_log, ditsave_dat, ditsave_slice, ditsave_collapsed, ditsave_custom,&
         tsavestart_log, tsavestart_dat, tsavestart_slice, tsavestart_collapsed,&
         tsavestart_custom, tsavestart

    namelist /stoplist/ it_init,time_init,it_max,time_max,dtmin,reset_it,reset_time,&
         wall_time_max,final_dt_reduction

    namelist /methodlist/ time_stepper,time_integrator, &
         source_split_usr,typesourcesplit,local_timestep,&
         dimsplit,typedimsplit,flux_scheme,&
         limiter,gradient_limiter,cada3_radius,&
         loglimit,typeboundspeed, H_correction,&
         typetvd,typeentropy,entropycoef,typeaverage, &
         typegrad,typediv,typecurl,&
         nxdiffusehllc, flathllc, tvdlfeps, flux_adaptive_diffusion, &
         flatcd,flatsh,&
         rk2_alfa,imex222_lambda,ssprk_order,rk3_switch,imex_switch,&
         small_temperature,small_pressure,small_density, &
         small_values_method, small_values_daverage, fix_small_values, check_small_values, &
         trace_small_values, small_values_fix_iw, &
         schmid_rad^D

    namelist /boundlist/ nghostcells,ghost_copy,&
         internalboundary, typeboundary_^L, save_physical_boundary

    namelist /meshlist/ refine_max_level,nbufferx^D,refine_threshold,&
         derefine_ratio, refine_criterion, &
         stretch_dim, stretch_uncentered, &
         qstretch_baselevel, nstretchedblocks_baselevel, &
         amr_wavefilter,max_blocks,block_nx^D,domain_nx^D,iprob,xprob^L, &
         w_refine_weight, prolongprimitive,coarsenprimitive, &
         typeprolonglimit, &
         logflag,tfixgrid,itfixgrid,ditregrid
    namelist /paramlist/  courantpar, dtpar, dtdiffpar, &
         typecourant, slowsteps

    namelist /emissionlist/ filename_euv,wavelength,&
          filename_sxr,emin_sxr,emax_sxr,&
          LOS_theta,LOS_phi,image_rotate,x_origin,big_image,&
          spectrum_wl,location_slit,filename_spectrum,&
          spectrum_window_min,spectrum_window_max,&
          instrument_resolution_factor,activate_unit_arcsec,&
          filename_whitelight,whitelight_instrument,R_occultor,R_opt_thick,&
          dat_resolution,direction_slit

    ! default maximum number of grid blocks in a processor
    max_blocks=4000

    ! allocate cell size of all levels
    allocate(dx(ndim,nlevelshi))
    {allocate(dg^D(nlevelshi))\}
    {allocate(ng^D(nlevelshi))\}

    ! default block size excluding ghost cells
    {block_nx^D = 16\}

    ! default resolution of level-1 mesh (full domain)
    {domain_nx^D = 32\}

    !! default number of ghost-cell layers at each boundary of a block
    ! this is now done when the variable is defined in mod_global_parameters
    ! the physics modules might set this variable in their init subroutine called earlier
    nghostcells = 2

    ! Allocate boundary conditions arrays in new and old style
    {
    allocate(typeboundary_min^D(nwfluxbc))
    allocate(typeboundary_max^D(nwfluxbc))
    typeboundary_min^D = undefined
    typeboundary_max^D = undefined
    }

    allocate(typeboundary(nwflux+nwaux,2*ndim))
    typeboundary=0

    ! not save physical boundary in dat files by default
    save_physical_boundary = .false.

    internalboundary   = .false.

    ! defaults for specific options
    typegrad   = 'central'
    typediv    = 'central'
    typecurl   = 'central'

    ! defaults for smallest physical values allowed
    small_temperature = 0.d0
    small_pressure    = 0.d0
    small_density     = 0.d0

    allocate(small_values_fix_iw(nw))
    small_values_fix_iw(:) = .true.

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
    {nbufferx^D                 = 0\}
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

    ! Grid stretching defaults
    stretch_uncentered = .true.
    stretch_dim(1:ndim) = undefined
    qstretch_baselevel(1:ndim) = bigdouble
    nstretchedblocks_baselevel(1:ndim)=0

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
        tsave(isave,ifile)  = bigdouble  ! global_time  of saves into the output files
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
    typeaverage     = 'default'
    tvdlfeps        = one
    flux_adaptive_diffusion = .false.
    nxdiffusehllc   = 0
    flathllc        = .false.
    slowsteps       = -1
    courantpar      = 0.8d0
    typecourant     = 'maxsum'
    dimsplit        = .false.
    typedimsplit    = 'default'
    if(physics_type=='mhd') then
      cada3_radius  = 0.1d0
    else
      cada3_radius  = 0.1d0
    end if
    {schmid_rad^D = 1.d0\}
    typetvd         = 'roe'
    typeboundspeed  = 'Einfeldt'
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

    allocate(flux_scheme(nlevelshi),typepred1(nlevelshi),flux_method(nlevelshi))
    allocate(limiter(nlevelshi),gradient_limiter(nlevelshi))
    do level=1,nlevelshi
       flux_scheme(level) = 'tvdlf'
       typepred1(level)   = 0
       limiter(level)     = 'minmod'
       gradient_limiter(level) = 'minmod'
    end do

    flatcd          = .false.
    flatsh          = .false.
    typesourcesplit = 'sfs'
    allocate(loglimit(nw))
    loglimit(1:nw)  = .false.

    allocate(typeentropy(nw))

    do iw=1,nw
       typeentropy(iw)='nul'      ! Entropy fix type
    end do

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

    ! These are used to construct file and log names from multiple par files
    basename_full = ''
    basename_prev = ''

    do i = 1, size(par_files)
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
       if (base_filename /= basename_prev) &
            basename_full = trim(basename_full) // trim(base_filename)
       basename_prev = base_filename
    end do

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
        if(mype==0) write(*,*) "No snapshots found to resume from, start a new run..."
      else
        ! Set file name to restart from
        write(restart_from_file, "(a,i4.4,a)") trim(base_filename),index_latest_data, ".dat"
      end if
    end if

    if (restart_from_file == undefined) then
      snapshotnext = 0
      slicenext    = 0
      collapsenext = 0
      if (firstprocess) &
           call mpistop("Please restart from a snapshot when firstprocess=T")
      if (convert) &
           call mpistop('Change convert to .false. for a new run!')
    end if

    if (small_pressure < 0.d0) call mpistop("small_pressure should be positive.")
    if (small_density < 0.d0) call mpistop("small_density should be positive.")
    ! Give priority to non-zero small temperature
    if (small_temperature>0.d0) small_pressure=small_density*small_temperature

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
       write(unitterm, *) 'Output type | tsavestart |  dtsave   | ditsave | itsave(1) | tsave(1)'
       write(fmt_string, *) '(A12," | ",E9.3E2,"  | ",E9.3E2," | ",I6,"  | "'//&
            ',I6, "    | ",E9.3E2)'
    end if

    do ifile = 1, nfile
       if (mype == 0) write(unitterm, fmt_string) trim(output_names(ifile)), &
            tsavestart(ifile), dtsave(ifile), ditsave(ifile), itsave(1, ifile), tsave(1, ifile)
    end do

    if (mype == 0) write(unitterm, *) ''

    do islice=1,nslices
       if(slicedir(islice) > ndim) &
            write(uniterr,*)'Warning in read_par_files: ', &
            'Slice ', islice,' direction',slicedir(islice),'larger than ndim=',ndim
       if(slicedir(islice) < 1) &
            write(uniterr,*)'Warning in read_par_files: ', &
            'Slice ', islice,' direction',slicedir(islice),'too small, should be [',1,ndim,']'
    end do

    if(it_max==biginteger .and. time_max==bigdouble.and.mype==0) write(uniterr,*) &
         'Warning in read_par_files: it_max or time_max not given!'

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

    do level=1,nlevelshi
       select case (flux_scheme(level))
       case ('hll')
          flux_method(level)=fs_hll
       case ('hllc')
          flux_method(level)=fs_hllc
       case ('hlld')
          flux_method(level)=fs_hlld
       case ('hllcd')
          flux_method(level)=fs_hllcd
       case ('tvdlf')
          flux_method(level)=fs_tvdlf
       case ('tvdmu')
          flux_method(level)=fs_tvdmu
       case ('tvd')
          flux_method(level)=fs_tvd
       case ('cd')
          flux_method(level)=fs_cd
       case ('cd4')
          flux_method(level)=fs_cd4
       case ('fd')
          flux_method(level)=fs_fd
       case ('source')
          flux_method(level)=fs_source
       case ('nul','null')
          flux_method(level)=fs_nul
       case default
          call mpistop("unkown or bad flux scheme")
       end select
       if(flux_scheme(level)=='tvd'.and.time_stepper/='onestep') &
            call mpistop(" tvd is onestep method, reset time_stepper='onestep'")
       if(flux_scheme(level)=='tvd')then
          if(mype==0.and.(.not.dimsplit)) write(unitterm,*) &
               'Warning: setting dimsplit=T for tvd, as used for level=',level
          dimsplit=.true.
       endif
       if(flux_scheme(level)=='hlld'.and.physics_type/='mhd' .and. physics_type/='twofl') &
          call mpistop("Cannot use hlld flux if not using MHD or 2FL only charges physics!")

       if(flux_scheme(level)=='hllc'.and.physics_type=='mf') &
          call mpistop("Cannot use hllc flux if using magnetofriction physics!")

       if(flux_scheme(level)=='tvd'.and.physics_type=='mf') &
          call mpistop("Cannot use tvd flux if using magnetofriction physics!")

       if(flux_scheme(level)=='tvdmu'.and.physics_type=='mf') &
          call mpistop("Cannot use tvdmu flux if using magnetofriction physics!")

       if (typepred1(level)==0) then
          select case (flux_scheme(level))
          case ('cd')
             typepred1(level)=fs_cd
          case ('cd4')
             typepred1(level)=fs_cd4
          case ('fd')
             typepred1(level)=fs_fd
          case ('tvdlf','tvdmu')
             typepred1(level)=fs_hancock
          case ('hll')
             typepred1(level)=fs_hll
          case ('hllc')
             typepred1(level)=fs_hllc
          case ('hllcd')
             typepred1(level)=fs_hllcd
          case ('hlld')
             typepred1(level)=fs_hlld
          case ('nul','source','tvd')
             typepred1(level)=fs_nul
          case default
             call mpistop("No default predictor for this full step")
          end select
       end if
    end do

    ! finite difference scheme fd need global maximal speed
    if(any(flux_scheme=='fd')) need_global_cmax=.true.
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
       nstep=1
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
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",time_stepper
          call mpistop("unkown onestep time_integrator in read_par_files")
       end select
       use_imex_scheme=(t_integrator==IMEX_Euler.or.t_integrator==IMEX_SP)
    case ("twostep")
       t_stepper=twostep
       nstep=2
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
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",time_stepper
          call mpistop("unkown twostep time_integrator in read_par_files")
       end select
       use_imex_scheme=(t_integrator==IMEX_Midpoint.or.t_integrator==IMEX_Trapezoidal&
            .or.t_integrator==IMEX_222)
       if (t_integrator==RK2_alf) then
          if(rk2_alfa<smalldouble.or.rk2_alfa>one)call mpistop("set rk2_alfa within [0,1]")
          rk_a21=rk2_alfa 
          rk_b2=1.0d0/(2.0d0*rk2_alfa)
          rk_b1=1.0d0-rk_b2
       endif
    case ("threestep")
       t_stepper=threestep
       nstep=3
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
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",time_stepper
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
          imex_b3   = (6.0d0*imex_c2**2 - 3.0d0*imex_c2 + 1.0d0) / (6.0d0*imex_c2**2)
          imex_a33  = (1.0d0/6.0d0 - imex_b2*imex_c2**2 - imex_b3*imex_c2*imex_c3) / (imex_b3*(imex_c3-imex_c2))
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
       use_imex_scheme=(t_integrator==IMEX_ARS3.or.t_integrator==IMEX_232.or.t_integrator==IMEX_CB3a)
    case ("fourstep")
       t_stepper=fourstep
       nstep=4
       if (time_integrator=='default') then
          time_integrator="ssprk4"
       endif
       select case (time_integrator)
       case ("ssprk4")
          t_integrator=ssprk4
       case ("rk4")
          t_integrator=rk4
       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",time_stepper
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
       nstep=5
       if (time_integrator=='default') then
          time_integrator="ssprk5"
       end if
       select case (time_integrator)
       case ("ssprk5")
          t_integrator=ssprk5
       case default
          write(unitterm,*) "time_integrator=",time_integrator,"time_stepper=",time_stepper
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
            rk_beta55=1.0d0-rk_beta54-rk_alfa53*rk_c3-rk_alfa54*rk_c4-rk_alfa55*rk_c5
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

    do i = 1, ndim
       select case (stretch_dim(i))
       case (undefined, 'none')
          stretch_type(i) = stretch_none
          stretched_dim(i) = .false.
       case ('uni','uniform')
          stretch_type(i) = stretch_uni
          stretched_dim(i) = .true.
       case ('symm', 'symmetric')
          stretch_type(i) = stretch_symm
          stretched_dim(i) = .true.
       case default
          stretch_type(i) = stretch_none
          stretched_dim(i) = .false.
          if (mype == 0) print *, 'Got stretch_type = ', stretch_type(i)
          call mpistop('Unknown stretch type')
       end select
    end do

    ! Harmonize the parameters for dimensional splitting and source splitting
    if(typedimsplit   =='default'.and.     dimsplit)   typedimsplit='xyyx'
    if(typedimsplit   =='default'.and..not.dimsplit)   typedimsplit='unsplit'
    dimsplit   = typedimsplit   /='unsplit'

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
      if(any(stretched_dim)) then
        coordinate=Cartesian_stretched
        slab_uniform=.false.
      end if
    else
      slab=.false.
      slab_uniform=.false.
    end if

    if(coordinate==spherical) then
      if(dimsplit) then
        if(mype==0)print *,'Warning: spherical symmetry needs dimsplit=F, resetting'
        dimsplit=.false.
      end if
    end if

    if (ndim==1) dimsplit=.false.

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

    if (any(limiter(1:nlevelshi)== 'ppm')&
         .and.(flatsh.and.physics_type=='rho')) then
       call mpistop(" PPM with flatsh=.true. can not be used with physics_type='rho'!")
    end if

    ! Copy boundary conditions to typeboundary, which is used internally
    {
    do iw=1,nwfluxbc
      select case(typeboundary_min^D(iw))
      case("special")
        typeboundary(iw,2*^D-1)=bc_special
      case("cont")
        typeboundary(iw,2*^D-1)=bc_cont
      case("symm")
        typeboundary(iw,2*^D-1)=bc_symm
      case("asymm")
        typeboundary(iw,2*^D-1)=bc_asymm
      case("periodic")
        typeboundary(iw,2*^D-1)=bc_periodic
      case("aperiodic")
        typeboundary(iw,2*^D-1)=bc_aperiodic
      case("noinflow")
        typeboundary(iw,2*^D-1)=bc_noinflow
      case("pole")
        typeboundary(iw,2*^D-1)=12
      case("bc_data")
        typeboundary(iw,2*^D-1)=bc_data
      case("bc_icarus")
        typeboundary(iw,2*^D-1)=bc_icarus
      case("character")
        typeboundary(iw,2*^D-1)=bc_character
      case default
         write (unitterm,*) "Undefined boundarytype found in read_par_files", &
           typeboundary_min^D(iw),"for variable iw=",iw," and side iB=",2*^D-1
      end select
    end do
    do iw=1,nwfluxbc
      select case(typeboundary_max^D(iw))
      case("special")
        typeboundary(iw,2*^D)=bc_special
      case("cont")
        typeboundary(iw,2*^D)=bc_cont
      case("symm")
        typeboundary(iw,2*^D)=bc_symm
      case("asymm")
        typeboundary(iw,2*^D)=bc_asymm
      case("periodic")
        typeboundary(iw,2*^D)=bc_periodic
      case("aperiodic")
        typeboundary(iw,2*^D)=bc_aperiodic
      case("noinflow")
        typeboundary(iw,2*^D)=bc_noinflow
      case("pole")
        typeboundary(iw,2*^D)=12
      case("bc_data")
        typeboundary(iw,2*^D)=bc_data
      case("bc_icarus")
        typeboundary(iw,2*^D-1)=bc_icarus
      case("bc_character")
        typeboundary(iw,2*^D)=bc_character
      case default
         write (unitterm,*) "Undefined boundarytype found in read_par_files", &
           typeboundary_max^D(iw),"for variable iw=",iw," and side iB=",2*^D
      end select
    end do
    }

    ! psi, tracers take the same boundary type as the first variable
    if (nwfluxbc<nwflux) then
      do iw=nwfluxbc+1,nwflux
        typeboundary(iw,:) = typeboundary(1, :)
      end do
    end if
    ! auxiliary variables take the same boundary type as the first variable
    if (nwaux>0) then
      do iw=nwflux+1, nwflux+nwaux
        typeboundary(iw,:) = typeboundary(1, :)
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
             if (typeboundary(iw,2*idim-1) .ne. typeboundary(iw,2*idim)) &
                  call mpistop("Wrong counterpart in periodic boundary")

             if (typeboundary(iw,2*idim-1) /= bc_periodic .and. &
                  typeboundary(iw,2*idim-1) /= bc_aperiodic) then
               call mpistop("Each dimension should either have all "//&
                    "or no variables periodic, some can be aperiodic")
             end if
          end do
       end if
    end do
    {^NOONED
    do idim=1,ndim
      if(any(typeboundary(:,2*idim-1)==12)) then
        if(any(typeboundary(:,2*idim-1)/=12)) typeboundary(:,2*idim-1)=12
        select case(physics_type)
        case ('rho','ard','rd','nonlinear','ffhd')
           ! all symmetric at pole
           typeboundary(:,2*idim-1)=bc_symm
           if(mype==0) print *,'symmetric minimal pole'
        case ('hd','rhd','srhd','mhd','rmhd')
           typeboundary(:,2*idim-1)=bc_symm
           ! here we assume the ordering of variables is fixed to rho-mom-[e]-B
           if(phys_energy) then
            windex=2
           else
            windex=1
           end if
           select case(coordinate)
           case(cylindrical)
            typeboundary(phi_+1,2*idim-1)=bc_asymm
            if(physics_type=='mhd'.or.physics_type=='rmhd') typeboundary(ndir+windex+phi_,2*idim-1)=bc_asymm
           case(spherical)
            typeboundary(3:ndir+1,2*idim-1)=bc_asymm
            if(physics_type=='mhd'.or.physics_type=='rmhd') typeboundary(ndir+windex+2:ndir+windex+ndir,2*idim-1)=bc_asymm
           case default
            call mpistop('Pole is in cylindrical, polar, spherical coordinates!')
           end select
        case ('twofl','mf')
           call mpistop('Pole treatment for twofl or mf not implemented yet')
        case default
           call mpistop('unknown physics type for setting minimal pole boundary treatment')
        end select
      end if
      if(any(typeboundary(:,2*idim)==12)) then
        if(any(typeboundary(:,2*idim)/=12)) typeboundary(:,2*idim)=12
        select case(physics_type)
        case ('rho','ard','rd','nonlinear','ffhd')
           ! all symmetric at pole
           typeboundary(:,2*idim)=bc_symm
           if(mype==0) print *,'symmetric maximal pole'
        case ('hd','rhd','srhd','mhd','rmhd')
           typeboundary(:,2*idim)=bc_symm
           ! here we assume the ordering of variables is fixed to rho-mom-[e]-B
           if(phys_energy) then
            windex=2
           else
            windex=1
           end if
           select case(coordinate)
           case(cylindrical)
            typeboundary(phi_+1,2*idim)=bc_asymm
            if(physics_type=='mhd'.or.physics_type=='rmhd') typeboundary(ndir+windex+phi_,2*idim)=bc_asymm
           case(spherical)
            typeboundary(3:ndir+1,2*idim)=bc_asymm
            if(physics_type=='mhd'.or.physics_type=='rmhd') typeboundary(ndir+windex+2:ndir+windex+ndir,2*idim)=bc_asymm
           case default
            call mpistop('Pole is in cylindrical, polar, spherical coordinates!')
           end select
        case ('twofl','mf')
           call mpistop('Pole treatment for twofl or mf not implemented yet')
        case default
           call mpistop('unknown physics type for setting maximal pole boundary treatment')
        end select
      end if
    end do
    }

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

      select case (coordinate)
         {^NOONED
      case (spherical)
         xprob^LIM^DE=xprob^LIM^DE*two*dpi;
         \}
      case (cylindrical)
         {
         if (^D==phi_) then
            xprob^LIM^D=xprob^LIM^D*two*dpi;
         end if
         \}
      end select

    ! full block size including ghostcells
    {ixGhi^D = block_nx^D + 2*nghostcells\}
    {ixGshi^D = ixGhi^D\}

    nx_vec = [{domain_nx^D|, }]
    block_nx_vec = [{block_nx^D|, }]

    if (any(nx_vec < 4) .or. any(mod(nx_vec, 2) == 1)) &
         call mpistop('Grid size (domain_nx^D) has to be even and >= 4')

    if (any(block_nx_vec < 4) .or. any(mod(block_nx_vec, 2) == 1)) &
         call mpistop('Block size (block_nx^D) has to be even and >= 4')

    { if(mod(domain_nx^D,block_nx^D)/=0) &
       call mpistop('Grid (domain_nx^D) and block (block_nx^D) must be consistent') \}

    if(refine_max_level>nlevelshi.or.refine_max_level<1)then
       write(unitterm,*)'Error: refine_max_level',refine_max_level,'>nlevelshi ',nlevelshi
       call mpistop("Reset nlevelshi and recompile!")
    endif

    if (any(stretched_dim)) then
       allocate(qstretch(0:nlevelshi,1:ndim),dxfirst(0:nlevelshi,1:ndim),&
             dxfirst_1mq(0:nlevelshi,1:ndim),dxmid(0:nlevelshi,1:ndim))
       allocate(nstretchedblocks(1:nlevelshi,1:ndim))
       qstretch(0:nlevelshi,1:ndim)=0.0d0
       dxfirst(0:nlevelshi,1:ndim)=0.0d0
       nstretchedblocks(1:nlevelshi,1:ndim)=0
       {if (stretch_type(^D) == stretch_uni) then
           ! first some sanity checks
           if(qstretch_baselevel(^D)<1.0d0.or.qstretch_baselevel(^D)==bigdouble) then
             if(mype==0) then
               write(*,*) 'stretched grid needs finite qstretch_baselevel>1'
               write(*,*) 'will try default value for qstretch_baselevel in dimension', ^D
             endif
             if(xprobmin^D>smalldouble)then
               qstretch_baselevel(^D)=(xprobmax^D/xprobmin^D)**(1.d0/dble(domain_nx^D))
             else
               call mpistop("can not set qstretch_baselevel automatically")
             endif
           endif
           if(mod(block_nx^D,2)==1) &
             call mpistop("stretched grid needs even block size block_nxD")
           if(mod(domain_nx^D/block_nx^D,2)/=0) &
             call mpistop("number level 1 blocks in D must be even")
           qstretch(1,^D)=qstretch_baselevel(^D)
           dxfirst(1,^D)=(xprobmax^D-xprobmin^D) &
                *(1.0d0-qstretch(1,^D))/(1.0d0-qstretch(1,^D)**domain_nx^D)
           qstretch(0,^D)=qstretch(1,^D)**2
           dxfirst(0,^D)=dxfirst(1,^D)*(1.0d0+qstretch(1,^D))
           if(refine_max_level>1)then
              do ilev=2,refine_max_level
                 qstretch(ilev,^D)=dsqrt(qstretch(ilev-1,^D))
                 dxfirst(ilev,^D)=dxfirst(ilev-1,^D) &
                       /(1.0d0+dsqrt(qstretch(ilev-1,^D)))
              enddo
           endif
        endif \}
        if(mype==0) then
           {if(stretch_type(^D) == stretch_uni) then
              write(*,*) 'Stretched dimension ', ^D
              write(*,*) 'Using stretched grid with qs=',qstretch(0:refine_max_level,^D)
              write(*,*) '        and first cell sizes=',dxfirst(0:refine_max_level,^D)
           endif\}
        end if
       {if(stretch_type(^D) == stretch_symm) then
           if(mype==0) then
               write(*,*) 'will apply symmetric stretch in dimension', ^D
           endif
           if(mod(block_nx^D,2)==1) &
             call mpistop("stretched grid needs even block size block_nxD")
           ! checks on the input variable nstretchedblocks_baselevel
           if(nstretchedblocks_baselevel(^D)==0) &
             call mpistop("need finite even number of stretched blocks at baselevel")
           if(mod(nstretchedblocks_baselevel(^D),2)==1) &
             call mpistop("need even number of stretched blocks at baselevel")
           if(qstretch_baselevel(^D)<1.0d0.or.qstretch_baselevel(^D)==bigdouble) &
             call mpistop('stretched grid needs finite qstretch_baselevel>1')
           ! compute stretched part to ensure uniform center
           ipower=(nstretchedblocks_baselevel(^D)/2)*block_nx^D
           if(nstretchedblocks_baselevel(^D)==domain_nx^D/block_nx^D)then
              xstretch^D=0.5d0*(xprobmax^D-xprobmin^D)
           else
              xstretch^D=(xprobmax^D-xprobmin^D) &
                /(2.0d0+dble(domain_nx^D-nstretchedblocks_baselevel(^D)*block_nx^D) &
               *(1.0d0-qstretch_baselevel(^D))/(1.0d0-qstretch_baselevel(^D)**ipower))
           endif
           if(xstretch^D>(xprobmax^D-xprobmin^D)*0.5d0) &
             call mpistop(" stretched grid part should not exceed full domain")
           dxfirst(1,^D)=xstretch^D*(1.0d0-qstretch_baselevel(^D)) &
                                   /(1.0d0-qstretch_baselevel(^D)**ipower)
           nstretchedblocks(1,^D)=nstretchedblocks_baselevel(^D)
           qstretch(1,^D)=qstretch_baselevel(^D)
           qstretch(0,^D)=qstretch(1,^D)**2
           dxfirst(0,^D)=dxfirst(1,^D)*(1.0d0+qstretch(1,^D))
           dxmid(1,^D)=dxfirst(1,^D)
           dxmid(0,^D)=dxfirst(1,^D)*2.0d0
           if(refine_max_level>1)then
              do ilev=2,refine_max_level
                 nstretchedblocks(ilev,^D)=2*nstretchedblocks(ilev-1,^D)
                 qstretch(ilev,^D)=dsqrt(qstretch(ilev-1,^D))
                 dxfirst(ilev,^D)=dxfirst(ilev-1,^D) &
                         /(1.0d0+dsqrt(qstretch(ilev-1,^D)))
                 dxmid(ilev,^D)=dxmid(ilev-1,^D)/2.0d0
              enddo
           endif
           ! sanity check on total domain size:
           sizeuniformpart^D=dxfirst(1,^D) &
             *(domain_nx^D-nstretchedblocks_baselevel(^D)*block_nx^D)
           if(mype==0) then
              print *,'uniform part of size=',sizeuniformpart^D
              print *,'setting of domain is then=',2*xstretch^D+sizeuniformpart^D
              print *,'versus=',xprobmax^D-xprobmin^D
           endif
           if(dabs(xprobmax^D-xprobmin^D-2*xstretch^D-sizeuniformpart^D)>smalldouble) then
              call mpistop('mismatch in domain size!')
           endif
        endif \}
        dxfirst_1mq(0:refine_max_level,1:ndim)=dxfirst(0:refine_max_level,1:ndim) &
                              /(1.0d0-qstretch(0:refine_max_level,1:ndim))
    end if

    dx_vec = [{xprobmax^D-xprobmin^D|, }] / nx_vec

    if (mype==0) then
       write(c_ndim, '(I1)') ^ND
       write(unitterm, '(A30,' // c_ndim // '(I0," "))') &
            ' Domain size (cells): ', nx_vec
       write(unitterm, '(A30,' // c_ndim // '(E9.3," "))') &
            ' Level one dx: ', dx_vec
    end if

    if (any(dx_vec < smalldouble)) &
         call mpistop("Incorrect domain size (too small grid spacing)")

    dx(:, 1) = dx_vec

    if(sum(w_refine_weight(:))==0) w_refine_weight(1) = 1.d0
    if(dabs(sum(w_refine_weight(:))-1.d0)>smalldouble) then
      write(unitterm,*) "Sum of all elements in w_refine_weight be 1.d0"
      call mpistop("Reset w_refine_weight so the sum is 1.d0")
    end if

    select case (typeboundspeed)
    case('Einfeldt')
      boundspeed=1
    case('cmaxmean')
      boundspeed=2
    case('cmaxleftright')
      boundspeed=3
    case default
      call mpistop("set typeboundspeed='Einfeldt' or 'cmaxmean' or 'cmaxleftright'")
    end select

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
       if(mype==0)print*,'Warning, at iteration=',itfixgrid,'the grid will be fixed'
    end if
    if (ditregrid>1) then
       if(mype==0)print*,'Note, Grid is reconstructed once every',ditregrid,'iterations'
    end if


    do islice=1,nslices
       select case(slicedir(islice))
          {case(^D)
          if(slicecoord(islice)<xprobmin^D.or.slicecoord(islice)>xprobmax^D) &
               write(uniterr,*)'Warning in read_par_files: ', &
               'Slice ', islice, ' coordinate',slicecoord(islice),'out of bounds for dimension ',slicedir(islice)
          \}
       end select
    end do

    if (mype==0) then
       write(unitterm, '(A30,A,A)') 'restart_from_file: ', ' ', trim(restart_from_file)
       write(unitterm, '(A30,L1)') 'converting: ', convert
       write(unitterm, '(A)') ''
    endif

    deallocate(flux_scheme)

  end subroutine read_par_files

  !> Routine to find entries in a string
  subroutine get_fields_string(line, delims, n_max, fields, n_found, fully_read)
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
       ixs_start(n_found) = ix_prev + ix ! This is the absolute position in 'line'

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
    integer:: ifile

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
    integer(kind=MPI_OFFSET_KIND), intent(in) :: offset_tree  !< Offset of tree info
    integer(kind=MPI_OFFSET_KIND), intent(in) :: offset_block !< Offset of block data
    call snapshot_write_header1(fh, offset_tree, offset_block, cons_wnames, nw)
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
    integer(MPI_OFFSET_KIND), intent(out) :: offset_tree  !< Offset of tree info
    integer(MPI_OFFSET_KIND), intent(out) :: offset_block !< Offset of block data

    double precision                      :: rbuf(ndim)
    double precision, allocatable         :: params(:)
    integer                               :: i, version
    integer                               :: ibuf(ndim), iw
    integer                               :: er, n_par, tmp_int
    integer, dimension(MPI_STATUS_SIZE)   :: st
    logical                               :: periodic(ndim)
    character(len=name_len), allocatable  :: var_names(:), param_names(:)
    character(len=name_len)               :: phys_name, geom_name

    ! Version number
    call MPI_FILE_READ(fh, version, 1, MPI_INTEGER, st, er)
    if (all(compatible_versions /= version)) then
      call mpistop("Incompatible file version (maybe old format?)")
    end if

    ! offset_tree
    call MPI_FILE_READ(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    offset_tree = ibuf(1)

    ! offset_block
    call MPI_FILE_READ(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    offset_block = ibuf(1)

    ! nw
    call MPI_FILE_READ(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    nw_found=ibuf(1)
    if (nw /= ibuf(1)) then
      write(*,*) "nw=",nw," and nw found in restart file=",ibuf(1)
      write(*,*) "Please be aware of changes in w at restart."
      !call mpistop("currently, changing nw at restart is not allowed")
    end if

    ! ndir
    call MPI_FILE_READ(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    if (ibuf(1) /= ndir) then
      write(*,*) "ndir in restart file = ",ibuf(1)
      write(*,*) "ndir = ",ndir
      call mpistop("reset ndir to ndir in restart file")
    end if

    ! ndim
    call MPI_FILE_READ(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    if (ibuf(1) /= ndim) then
      write(*,*) "ndim in restart file = ",ibuf(1)
      write(*,*) "ndim = ",ndim
      call mpistop("reset ndim to ndim in restart file")
    end if

    ! levmax
    call MPI_FILE_READ(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    if (ibuf(1) > refine_max_level) then
      write(*,*) "number of levels in restart file = ",ibuf(1)
      write(*,*) "refine_max_level = ",refine_max_level
      call mpistop("refine_max_level < num. levels in restart file")
    end if

    ! nleafs
    call MPI_FILE_READ(fh, nleafs, 1, MPI_INTEGER, st, er)

    ! nparents
    call MPI_FILE_READ(fh, nparents, 1, MPI_INTEGER, st, er)

    ! it
    call MPI_FILE_READ(fh, it, 1, MPI_INTEGER, st, er)

    ! global time
    call MPI_FILE_READ(fh, global_time, 1, MPI_DOUBLE_PRECISION, st, er)

    ! xprobmin^D
    call MPI_FILE_READ(fh,rbuf(1:ndim),ndim,MPI_DOUBLE_PRECISION,st,er)
    if (maxval(abs(rbuf(1:ndim) - [ xprobmin^D ])) > 0) then
      write(*,*) "Error: xprobmin differs from restart data: ", rbuf(1:ndim)
      call mpistop("change xprobmin^D in par file")
    end if

    ! xprobmax^D
    call MPI_FILE_READ(fh,rbuf(1:ndim),ndim,MPI_DOUBLE_PRECISION,st,er)
    if (maxval(abs(rbuf(1:ndim) - [ xprobmax^D ])) > 0) then
      write(*,*) "Error: xprobmax differs from restart data: ", rbuf(1:ndim)
      call mpistop("change xprobmax^D in par file")
    end if

    ! domain_nx^D
    call MPI_FILE_READ(fh,ibuf(1:ndim), ndim, MPI_INTEGER,st,er)
    if (any(ibuf(1:ndim) /= [ domain_nx^D ])) then
      write(*,*) "Error: mesh size differs from restart data: ", ibuf(1:ndim)
      call mpistop("change domain_nx^D in par file")
    end if

    ! block_nx^D
    call MPI_FILE_READ(fh,ibuf(1:ndim), ndim, MPI_INTEGER,st,er)
    if (any(ibuf(1:ndim) /= [ block_nx^D ])) then
      write(*,*) "Error: block size differs from restart data:", ibuf(1:ndim)
      call mpistop("change block_nx^D in par file")
    end if

    ! From version 5, read more info about the grid
    if (version > 4) then
      call MPI_FILE_READ(fh, periodic, ndim, MPI_LOGICAL, st, er)
      if ({periodic(^D) .and. .not.periodB(^D) .or. .not.periodic(^D) .and. periodB(^D)| .or. }) &
           call mpistop("change in periodicity in par file")

      call MPI_FILE_READ(fh, geom_name, name_len, MPI_CHARACTER, st, er)

      if (geom_name /= geometry_name(1:name_len)) then
        write(*,*) "type of coordinates in data is: ", geom_name
        call mpistop("select the correct coordinates in mod_usr.t file")
      end if

      call MPI_FILE_READ(fh, stagger_mark_dat, 1, MPI_LOGICAL, st, er)
      if (stagger_grid .and. .not. stagger_mark_dat .or. .not.stagger_grid.and.stagger_mark_dat) then
        write(*,*) "Warning: stagger grid flag differs from restart data:", stagger_mark_dat
        !call mpistop("change parameter to use stagger grid")
      end if
    end if

    ! From version 4 onwards, the later parts of the header must be present
    if (version > 3) then
      ! w_names (not used here)
      allocate(var_names(nw_found))
      do iw = 1, nw_found
        call MPI_FILE_READ(fh, var_names(iw), name_len, MPI_CHARACTER, st, er)
      end do

      ! Physics related information
      call MPI_FILE_READ(fh, phys_name, name_len, MPI_CHARACTER, st, er)

      if (phys_name /= physics_type) then
!        call mpistop("Cannot restart with a different physics type")
      end if

      call MPI_FILE_READ(fh, n_par, 1, MPI_INTEGER, st, er)
      allocate(params(n_par))
      allocate(param_names(n_par))
      call MPI_FILE_READ(fh, params, n_par, MPI_DOUBLE_PRECISION, st, er)
      call MPI_FILE_READ(fh, param_names, name_len * n_par, MPI_CHARACTER, st, er)

      ! Read snapshotnext etc. for restarting
      call MPI_FILE_READ(fh, tmp_int, 1, MPI_INTEGER, st, er)

      ! Only set snapshotnext if the user hasn't specified it
      if (snapshotnext == -1) snapshotnext = tmp_int

      call MPI_FILE_READ(fh, tmp_int, 1, MPI_INTEGER, st, er)
      if (slicenext == -1) slicenext = tmp_int

      call MPI_FILE_READ(fh, tmp_int, 1, MPI_INTEGER, st, er)
      if (collapsenext == -1) collapsenext = tmp_int
    else
      ! Guess snapshotnext from file name if not set
      if (snapshotnext == -1) &
           snapshotnext = get_snapshot_index(trim(restart_from_file)) + 1
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
    use mod_input_output_helper, only: count_ix,block_shape_io,create_output_file
    use mod_functions_forest, only: write_forest

    double precision, allocatable :: w_buffer(:)
    integer                       :: file_handle, igrid, Morton_no, iwrite
    integer                       :: ipe, ix_buffer(2*ndim+1), n_values
    integer                       :: ixO^L, n_ghost(2*ndim)
    integer                       :: ixOs^L,n_values_stagger
    integer                       :: iorecvstatus(MPI_STATUS_SIZE)
    integer                       :: ioastatus(MPI_STATUS_SIZE)
    integer                       :: igrecvstatus(MPI_STATUS_SIZE)
    integer                       :: istatus(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset_tree_info
    integer(kind=MPI_OFFSET_KIND) :: offset_block_data
    integer(kind=MPI_OFFSET_KIND) :: offset_offsets
    integer, allocatable                       :: block_ig(:, :)
    integer, allocatable                       :: block_lvl(:)
    integer(kind=MPI_OFFSET_KIND), allocatable :: block_offset(:)
    type(tree_node), pointer      :: pnode

    call MPI_BARRIER(icomm, ierrmpi)

    ! Allocate send/receive buffer
    n_values = count_ix(ixG^LL) * nw
    if(stagger_grid) then
      n_values = n_values + count_ix(ixGs^LL) * nws
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
      call snapshot_write_header(file_handle, offset_tree_info, &
           offset_block_data)

      call MPI_File_get_position(file_handle, offset_tree_info, ierrmpi)

      call write_forest(file_handle)

      ! Collect information about the spatial index (ig^D) and refinement level
      ! of leaves
      do Morton_no = Morton_start(0), Morton_stop(npe-1)
         igrid = sfc(1, Morton_no)
         ipe = sfc(2, Morton_no)
         pnode => igrid_to_node(igrid, ipe)%node

         block_ig(:, Morton_no) = [ pnode%ig^D ]
         block_lvl(Morton_no) = pnode%level
         block_offset(Morton_no) = 0 ! Will be determined later
      end do

      call MPI_FILE_WRITE(file_handle, block_lvl, size(block_lvl), &
           MPI_INTEGER, istatus, ierrmpi)

      call MPI_FILE_WRITE(file_handle, block_ig, size(block_ig), &
           MPI_INTEGER, istatus, ierrmpi)

      ! Block offsets are currently unknown, but will be overwritten later
      call MPI_File_get_position(file_handle, offset_offsets, ierrmpi)
      call MPI_FILE_WRITE(file_handle, block_offset(1:nleafs), nleafs, &
           MPI_OFFSET, istatus, ierrmpi)

      call MPI_File_get_position(file_handle, offset_block_data, ierrmpi)

      ! Check whether data was written as expected
      if (offset_block_data - offset_tree_info /= &
           (nleafs + nparents) * size_logical + &
           nleafs * ((1+ndim) * size_int + 2 * size_int)) then
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

      call block_shape_io(igrid, n_ghost, ixO^L, n_values)
      if(stagger_grid) then
        w_buffer(1:n_values) = pack(ps(igrid)%w(ixO^S, 1:nw), .true.)
        {ixOsmin^D = ixOmin^D -1\}
        {ixOsmax^D = ixOmax^D \}
        n_values_stagger= count_ix(ixOs^L)*nws
        w_buffer(n_values+1:n_values+n_values_stagger) = pack(ps(igrid)%ws(ixOs^S, 1:nws), .true.)
        n_values=n_values+n_values_stagger
      else
        w_buffer(1:n_values) = pack(ps(igrid)%w(ixO^S, 1:nw), .true.)
      end if
      ix_buffer(1) = n_values
      ix_buffer(2:) = n_ghost

      if (mype /= 0) then
        call MPI_SEND(ix_buffer, 2*ndim+1, &
             MPI_INTEGER, 0, itag, icomm, ierrmpi)
        call MPI_SEND(w_buffer, n_values, &
             MPI_DOUBLE_PRECISION, 0, itag, icomm, ierrmpi)
      else
        iwrite = iwrite+1
        call MPI_FILE_WRITE(file_handle, ix_buffer(2:), &
             2*ndim, MPI_INTEGER, istatus, ierrmpi)
        call MPI_FILE_WRITE(file_handle, w_buffer, &
             n_values, MPI_DOUBLE_PRECISION, istatus, ierrmpi)

        ! Set offset of next block
        block_offset(iwrite+1) = block_offset(iwrite) + &
             int(n_values, MPI_OFFSET_KIND) * size_double + &
             2 * ndim * size_int
      end if
    end do

    ! Write data communicated from other processors
    if (mype == 0) then
      do ipe = 1, npe-1
        do Morton_no=Morton_start(ipe), Morton_stop(ipe)
          iwrite=iwrite+1
          itag=Morton_no

          call MPI_RECV(ix_buffer, 2*ndim+1, MPI_INTEGER, ipe, itag, icomm,&
               igrecvstatus, ierrmpi)
          n_values = ix_buffer(1)

          call MPI_RECV(w_buffer, n_values, MPI_DOUBLE_PRECISION,&
               ipe, itag, icomm, iorecvstatus, ierrmpi)

          call MPI_FILE_WRITE(file_handle, ix_buffer(2:), &
             2*ndim, MPI_INTEGER, istatus, ierrmpi)
          call MPI_FILE_WRITE(file_handle, w_buffer, &
             n_values, MPI_DOUBLE_PRECISION, istatus, ierrmpi)

          ! Set offset of next block
          block_offset(iwrite+1) = block_offset(iwrite) + &
               int(n_values, MPI_OFFSET_KIND) * size_double + &
               2 * ndim * size_int
        end do
      end do

      ! Write block offsets (now we know them)
      call MPI_FILE_SEEK(file_handle, offset_offsets, MPI_SEEK_SET, ierrmpi)
      call MPI_FILE_WRITE(file_handle, block_offset(1:nleafs), nleafs, &
           MPI_OFFSET, istatus, ierrmpi)

      ! Write header again, now with correct offsets
      call MPI_FILE_SEEK(file_handle, 0_MPI_OFFSET_KIND, MPI_SEEK_SET, ierrmpi)
      call snapshot_write_header(file_handle, offset_tree_info, &
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

    double precision :: ws(ixGs^T,1:ndim)
    double precision, allocatable :: w_buffer(:)
    double precision, dimension(:^D&,:), allocatable :: w
    integer                       :: ix_buffer(2*ndim+1), n_values, n_values_stagger
    integer                       :: ixO^L, ixOs^L
    integer                       :: file_handle, amode, igrid, Morton_no, iread
    integer                       :: istatus(MPI_STATUS_SIZE)
    integer                       :: iorecvstatus(MPI_STATUS_SIZE)
    integer                       :: ipe,inrecv,nrecv, file_version
    integer(MPI_OFFSET_KIND)      :: offset_tree_info
    integer(MPI_OFFSET_KIND)      :: offset_block_data
    logical                       :: fexist

    if (mype==0) then
      inquire(file=trim(restart_from_file), exist=fexist)
      if (.not.fexist) call mpistop(trim(restart_from_file)//" not found!")

      call MPI_FILE_OPEN(MPI_COMM_SELF,restart_from_file,MPI_MODE_RDONLY, &
           MPI_INFO_NULL,file_handle,ierrmpi)
      call MPI_FILE_READ(file_handle, file_version, 1, MPI_INTEGER, &
           istatus, ierrmpi)
    end if

    call MPI_BCAST(file_version,1,MPI_INTEGER,0,icomm,ierrmpi)

    if (all(compatible_versions /= file_version)) then
      if (mype == 0) print *, "Unknown version, trying old snapshot reader..."
      call MPI_FILE_CLOSE(file_handle,ierrmpi)
      call read_snapshot_old()

      ! Guess snapshotnext from file name if not set
      if (snapshotnext == -1) &
           snapshotnext = get_snapshot_index(trim(restart_from_file)) + 1
      ! Set slicenext and collapsenext if not set
      if (slicenext == -1) slicenext = 0
      if (collapsenext == -1) collapsenext = 0

      ! Still used in convert
      snapshotini = snapshotnext-1

      return ! Leave this routine
    else if (mype == 0) then
      call MPI_FILE_SEEK(file_handle, 0_MPI_OFFSET_KIND, MPI_SEEK_SET, ierrmpi)
      call snapshot_read_header(file_handle, offset_tree_info, &
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
    n_values = count_ix(ixG^LL) * nw_found
    if(stagger_mark_dat) then
      n_values = n_values + count_ix(ixGs^LL) * nws
    end if
    allocate(w_buffer(n_values))
    allocate(w(ixG^T,1:nw_found))

    nleafs_active = nleafs

    if (mype == 0) then
      call MPI_FILE_SEEK(file_handle, offset_tree_info, &
           MPI_SEEK_SET, ierrmpi)
    end if

    call read_forest(file_handle)

    do Morton_no=Morton_start(mype),Morton_stop(mype)
      igrid=sfc_to_igrid(Morton_no)
      call alloc_node(igrid)
    end do

    if (mype==0) then
      call MPI_FILE_SEEK(file_handle, offset_block_data, MPI_SEEK_SET, ierrmpi)

      iread = 0
      do ipe = 0, npe-1
        do Morton_no=Morton_start(ipe),Morton_stop(ipe)
          iread=iread+1
          itag=Morton_no

          call MPI_FILE_READ(file_handle,ix_buffer(1:2*ndim), 2*ndim, &
               MPI_INTEGER, istatus,ierrmpi)

          ! Construct ixO^L array from number of ghost cells
          {ixOmin^D = ixMlo^D - ix_buffer(^D)\}
          {ixOmax^D = ixMhi^D + ix_buffer(ndim+^D)\}
          n_values = count_ix(ixO^L) * nw_found
          if(stagger_mark_dat) then
            {ixOsmin^D = ixOmin^D - 1\}
            {ixOsmax^D = ixOmax^D\}
            n_values_stagger = n_values
            n_values = n_values + count_ix(ixOs^L) * nws
          end if

          call MPI_FILE_READ(file_handle, w_buffer, n_values, &
               MPI_DOUBLE_PRECISION, istatus, ierrmpi)

          if (mype == ipe) then ! Root task
            igrid=sfc_to_igrid(Morton_no)
            block=>ps(igrid)
            if(stagger_mark_dat) then
              w(ixO^S, 1:nw_found) = reshape(w_buffer(1:n_values_stagger), &
                   shape(w(ixO^S, 1:nw_found)))
              if(stagger_grid) &
              ps(igrid)%ws(ixOs^S,1:nws)=reshape(w_buffer(n_values_stagger+1:n_values), &
                   shape(ws(ixOs^S, 1:nws)))
            else
              w(ixO^S, 1:nw_found) = reshape(w_buffer(1:n_values), &
                   shape(w(ixO^S, 1:nw_found)))
            end if
            if (nw_found<nw) then
              if (associated(usr_transform_w)) then
                call usr_transform_w(ixG^LL,ixM^LL,nw_found,w,ps(igrid)%x,ps(igrid)%w)
              else
                ps(igrid)%w(ixO^S,1:nw_found)=w(ixO^S,1:nw_found)
              end if
            else if (nw_found>nw) then
              if (associated(usr_transform_w)) then
                call usr_transform_w(ixG^LL,ixM^LL,nw_found,w,ps(igrid)%x,ps(igrid)%w)
              else
                ps(igrid)%w(ixO^S,1:nw)=w(ixO^S,1:nw)
              end if
            else
              ps(igrid)%w(ixO^S,1:nw)=w(ixO^S,1:nw)
            end if
          else
            call MPI_SEND([ ixO^L, n_values ], 2*ndim+1, &
                 MPI_INTEGER, ipe, itag, icomm, ierrmpi)
            call MPI_SEND(w_buffer, n_values, &
                 MPI_DOUBLE_PRECISION, ipe, itag, icomm, ierrmpi)
          end if
        end do
      end do

      call MPI_FILE_CLOSE(file_handle,ierrmpi)

    else ! mype > 0

      do Morton_no=Morton_start(mype),Morton_stop(mype)
        igrid=sfc_to_igrid(Morton_no)
        block=>ps(igrid)
        itag=Morton_no

        call MPI_RECV(ix_buffer, 2*ndim+1, MPI_INTEGER, 0, itag, icomm,&
             iorecvstatus, ierrmpi)
        {ixOmin^D = ix_buffer(^D)\}
        {ixOmax^D = ix_buffer(ndim+^D)\}
        n_values = ix_buffer(2*ndim+1)

        call MPI_RECV(w_buffer, n_values, MPI_DOUBLE_PRECISION,&
             0, itag, icomm, iorecvstatus, ierrmpi)

        if(stagger_mark_dat) then
          n_values_stagger = count_ix(ixO^L) * nw_found
          {ixOsmin^D = ixOmin^D - 1\}
          {ixOsmax^D = ixOmax^D\}
          w(ixO^S, 1:nw_found) = reshape(w_buffer(1:n_values_stagger), &
               shape(w(ixO^S, 1:nw_found)))
          if(stagger_grid) &
          ps(igrid)%ws(ixOs^S,1:nws)=reshape(w_buffer(n_values_stagger+1:n_values), &
               shape(ws(ixOs^S, 1:nws)))
        else
          w(ixO^S, 1:nw_found) = reshape(w_buffer(1:n_values), &
               shape(w(ixO^S, 1:nw_found)))
        end if
        if (nw_found<nw) then
          if (associated(usr_transform_w)) then
            call usr_transform_w(ixG^LL,ixM^LL,nw_found,w,ps(igrid)%x,ps(igrid)%w)
          else
            ps(igrid)%w(ixO^S,1:nw_found)=w(ixO^S,1:nw_found)
          end if
        else if (nw_found>nw) then
          if (associated(usr_transform_w)) then
            call usr_transform_w(ixG^LL,ixM^LL,nw_found,w,ps(igrid)%x,ps(igrid)%w)
          else
            ps(igrid)%w(ixO^S,1:nw)=w(ixO^S,1:nw)
          end if
        else
          ps(igrid)%w(ixO^S,1:nw)=w(ixO^S,1:nw)
        end if
      end do
    end if

    call MPI_BARRIER(icomm,ierrmpi)

  end subroutine read_snapshot

  subroutine read_snapshot_old()
    use mod_forest
    use mod_global_parameters
    use mod_amr_solution_node, only: alloc_node
    use mod_functions_forest, only: read_forest

    double precision              :: wio(ixG^T,1:nw)
    double precision              :: eqpar_dummy(100)
    integer                       :: fh, igrid, Morton_no, iread
    integer                       :: levmaxini, ndimini, ndirini
    integer                       :: nwini, neqparini, nxini^D
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer                       :: istatus(MPI_STATUS_SIZE)
    integer, allocatable          :: iorecvstatus(:,:)
    integer                       :: ipe,inrecv,nrecv
    integer                       :: sendini(7+^ND)
    logical                       :: fexist
    character(len=80)             :: filename

    if (mype==0) then
      call MPI_FILE_OPEN(MPI_COMM_SELF,trim(restart_from_file), &
           MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierrmpi)

      offset=-int(7*size_int+size_double,kind=MPI_OFFSET_KIND)
      call MPI_FILE_SEEK(fh,offset,MPI_SEEK_END,ierrmpi)

      call MPI_FILE_READ(fh,nleafs,1,MPI_INTEGER,istatus,ierrmpi)
      nleafs_active = nleafs
      call MPI_FILE_READ(fh,levmaxini,1,MPI_INTEGER,istatus,ierrmpi)
      call MPI_FILE_READ(fh,ndimini,1,MPI_INTEGER,istatus,ierrmpi)
      call MPI_FILE_READ(fh,ndirini,1,MPI_INTEGER,istatus,ierrmpi)
      call MPI_FILE_READ(fh,nwini,1,MPI_INTEGER,istatus,ierrmpi)
      call MPI_FILE_READ(fh,neqparini,1,MPI_INTEGER,istatus,ierrmpi)
      call MPI_FILE_READ(fh,it,1,MPI_INTEGER,istatus,ierrmpi)
      call MPI_FILE_READ(fh,global_time,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)

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

      offset=offset-int(ndimini*size_int+neqparini*size_double,kind=MPI_OFFSET_KIND)
      call MPI_FILE_SEEK(fh,offset,MPI_SEEK_END,ierrmpi)

      {call MPI_FILE_READ(fh,nxini^D,1,MPI_INTEGER,istatus,ierrmpi)\}
      if (ixGhi^D/=nxini^D+2*nghostcells|.or.) then
        write(*,*) "Error: reset resolution to ",nxini^D+2*nghostcells
        call mpistop("change with setamrvac")
      end if

      call MPI_FILE_READ(fh,eqpar_dummy,neqparini, &
           MPI_DOUBLE_PRECISION,istatus,ierrmpi)
    end if

    ! broadcast the global parameters first
    if (npe>1) then
      if (mype==0) then
        sendini=(/nleafs,levmaxini,ndimini,ndirini,nwini,neqparini,it ,^D&nxini^D /)
      end if
      call MPI_BCAST(sendini,7+^ND,MPI_INTEGER,0,icomm,ierrmpi)
      nleafs=sendini(1);levmaxini=sendini(2);ndimini=sendini(3);
      ndirini=sendini(4);nwini=sendini(5);
      neqparini=sendini(6);it=sendini(7);
      nxini^D=sendini(7+^D);
      nleafs_active = nleafs
      call MPI_BCAST(global_time,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    end if

    if (mype == 0) then
      offset = int(size_block_io,kind=MPI_OFFSET_KIND) * &
           int(nleafs,kind=MPI_OFFSET_KIND)
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
        offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
             *int(Morton_no-1,kind=MPI_OFFSET_KIND)
        call MPI_FILE_READ_AT(fh,offset,ps(igrid)%w,1,type_block_io, &
             istatus,ierrmpi)
      end do
      if (npe>1) then
        do ipe=1,npe-1
          do Morton_no=Morton_start(ipe),Morton_stop(ipe)
            iread=iread+1
            itag=Morton_no
            offset=int(size_block_io,kind=MPI_OFFSET_KIND)&
                 *int(Morton_no-1,kind=MPI_OFFSET_KIND)
            call MPI_FILE_READ_AT(fh,offset,wio,1,type_block_io,&
                 istatus,ierrmpi)
            call MPI_SEND(wio,1,type_block_io,ipe,itag,icomm,ierrmpi)
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
        call MPI_RECV(ps(igrid)%w,1,type_block_io,0,itag,icomm,&
             iorecvstatus(:,inrecv),ierrmpi)
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

    double precision     :: dtTimeLast, now, cellupdatesPerSecond
    double precision     :: activeBlocksPerCore, wctPerCodeTime, timeToFinish
    double precision     :: wmean(1:nw), total_volume
    double precision     :: volume_coverage(refine_max_level)
    integer              :: i, iw, level
    integer              :: nx^D, nc, ncells, dit
    integer              :: amode, istatus(MPI_STATUS_SIZE)
    integer, parameter   :: my_unit = 20
    logical, save        :: opened  = .false.
    logical              :: fileopen
    character(len=40)    :: fmt_string
    character(len=80)    :: filename
    character(len=2048)  :: line

    ! Compute the volume-average of w**1 = w
    call get_volume_average(1, wmean, total_volume)

    ! Compute the volume coverage
    call get_volume_coverage(volume_coverage)

    if (mype == 0) then

       ! To compute cell updates per second, we do the following:
       nx^D=ixMhi^D-ixMlo^D+1;
       nc={nx^D*}
       ncells = nc * nleafs_active

       ! assumes the number of active leafs haven't changed since last compute.
       now        = MPI_WTIME()
       dit        = it - itTimeLast
       dtTimeLast = now - timeLast
       itTimeLast = it
       timeLast   = now
       cellupdatesPerSecond = dble(ncells) * dble(nstep) * &
            dble(dit) / (dtTimeLast * dble(npe))

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

          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, &
               MPI_INFO_NULL, log_fh, ierrmpi)

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
            call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'), &
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
       write(line(i:), fmt_string) '| ', Xload, Xmemory, cellupdatesPerSecond, &
            activeBlocksPerCore, wctPerCodeTime, timeToFinish

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

  end subroutine printlog_default

  !> Print a log that can be used to check whether the code still produces the
  !> same output (regression test)
  subroutine printlog_regression_test()
    use mod_global_parameters

    double precision   :: modes(nw, 2), volume
    integer, parameter :: n_modes = 2
    integer            :: power
    integer              :: amode, istatus(MPI_STATUS_SIZE)
    logical, save      :: file_open = .false.
    character(len=40)  :: fmt_string
    character(len=2048)  :: line
    character(len=80)    :: filename

    do power = 1, n_modes
       call get_volume_average(power, modes(:, power), volume)
    end do

    if (mype == 0) then
       if (.not. file_open) then
          filename = trim(base_filename) // ".log"
          amode    = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
          amode    = ior(amode,MPI_MODE_APPEND)

          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, &
               MPI_INFO_NULL, log_fh, ierrmpi)
          file_open = .true.

          line= "# time mean(w) mean(w**2)"
          call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'), &
                 len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
       end if

       write(fmt_string, "(a,i0,a)") "(", nw * n_modes + 1, fmt_r // ")"
       write(line, fmt_string) global_time, modes
       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
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

    double precision              :: wsum(nw+1)
    double precision              :: dsum_recv(1:nw+1)
    integer                       :: iigrid, igrid, iw

    wsum(:) = 0

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! Store total volume in last element
       wsum(nw+1) = wsum(nw+1) + sum(ps(igrid)%dvolume(ixM^T))

       ! Compute the modes of the cell-centered variables, weighted by volume
       do iw = 1, nw
          wsum(iw) = wsum(iw) + &
               sum(ps(igrid)%dvolume(ixM^T)*ps(igrid)%w(ixM^T,iw)**power)
       end do
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(wsum, dsum_recv, nw+1, MPI_DOUBLE_PRECISION, &
         MPI_SUM, icomm, ierrmpi)

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
       vol_cov(level) = vol_cov(level)+ &
            {(rnode(rpxmax^D_,igrid)-rnode(rpxmin^D_,igrid))|*}
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(vol_cov, dsum_recv, refine_max_level, MPI_DOUBLE_PRECISION, &
         MPI_SUM, icomm, ierrmpi)

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
    double precision              :: wsum(2)
    double precision              :: dsum_recv(2)
    integer                       :: iigrid, igrid, i^D

    wsum(:) = 0

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! Store total volume in last element
       wsum(2) = wsum(2) + sum(ps(igrid)%dvolume(ixM^T))

       ! Compute the modes of the cell-centered variables, weighted by volume
       {do i^D = ixMlo^D, ixMhi^D\}
       wsum(1) = wsum(1) + ps(igrid)%dvolume(i^D) * &
            func(ps(igrid)%w(i^D, :), nw)
       {end do\}
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(wsum, dsum_recv, 2, MPI_DOUBLE_PRECISION, &
         MPI_SUM, icomm, ierrmpi)

    ! Set the volume and the average
    volume = dsum_recv(2)
    f_avg  = dsum_recv(1) / volume

  end subroutine get_volume_average_func

  !> Compute global maxima of iw variables over the leaves of the grid.
  subroutine get_global_maxima(wmax)
    use mod_global_parameters

    double precision, intent(out) :: wmax(nw)  !< The global maxima

    double precision              :: wmax_mype(nw),wmax_recv(nw)
    integer                       :: iigrid, igrid, iw

    wmax_mype(1:nw) = -bigdouble

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       do iw = 1, nw
          wmax_mype(iw)=max(wmax_mype(iw),maxval(ps(igrid)%w(ixM^T,iw)))
       end do
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(wmax_mype, wmax_recv, nw, MPI_DOUBLE_PRECISION, &
         MPI_MAX, icomm, ierrmpi)

    wmax(1:nw)=wmax_recv(1:nw)

  end subroutine get_global_maxima

  !> Compute global minima of iw variables over the leaves of the grid.
  subroutine get_global_minima(wmin)
    use mod_global_parameters

    double precision, intent(out) :: wmin(nw)  !< The global maxima

    double precision              :: wmin_mype(nw),wmin_recv(nw)
    integer                       :: iigrid, igrid, iw

    wmin_mype(1:nw) = bigdouble

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       do iw = 1, nw
          wmin_mype(iw)=min(wmin_mype(iw),minval(ps(igrid)%w(ixM^T,iw)))
       end do
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(wmin_mype, wmin_recv, nw, MPI_DOUBLE_PRECISION, &
         MPI_MIN, icomm, ierrmpi)

    wmin(1:nw)=wmin_recv(1:nw)

  end subroutine get_global_minima

end module mod_input_output
