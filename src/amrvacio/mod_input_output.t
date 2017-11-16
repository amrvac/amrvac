!> Module for reading input and writing output
module mod_input_output

  implicit none
  public

  !> Version number of the .dat file output
  integer, parameter :: version_number = 3

  !> List of compatible versions
  integer, parameter :: compatible_versions(1) = [3]

  !> number of w found in dat files
  integer :: nw_found

  ! Formats used in output
  character(len=*), parameter :: fmt_r  = 'es16.8' ! Default precision
  character(len=*), parameter :: fmt_r2 = 'es10.2' ! Two digits
  character(len=*), parameter :: fmt_i  = 'i8'     ! Integer format

contains

  !> Read the command line arguments passed to amrvac
  subroutine read_arguments()
    use mod_kracken
    use mod_global_parameters
    use mod_slice, only: slicenext

    integer                          :: len, ier, n
    integer, parameter               :: max_files = 20 ! Maximum number of par files
    integer                          :: n_par_files
    integer                          :: ibegin(max_files)
    integer                          :: iterm(max_files)
    character(len=max_files*std_len) :: all_par_files
    character(len=std_len)           :: tmp_files(max_files)

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

    ! Specify the options and their default values
    call kracken('cmd','-i amrvac.par -if ' // undefined // &
         ' -slice 0 -collapse 0 --help .false. -convert .false.')

    ! Get the par file(s)
    call retrev('cmd_i', all_par_files, len, ier)

    ! Show the usage if the help flag was given, or no par file was specified
    if (lget('cmd_-help') .or. len == 0) then
       if (mype == 0) then
          print *, 'Usage example:'
          print *, 'mpirun -np 4 ./amrvac -i file.par [file2.par ...]'
          print *, ''
          print *, 'Optional arguments:'
          print *, '-convert             Convert snapshot files'
          print *, '-if file0001.dat     Use this snapshot to restart from'
          print *, ''
          print *, 'Note: later parameter files override earlier ones.'
       end if
       stop
    end if

    ! Split the input files, in case multiple were given
    call get_fields_string(all_par_files, " ,'"""//char(9), max_files, &
         tmp_files, n_par_files)

    allocate(par_files(n_par_files))
    par_files = tmp_files(1:n_par_files)

    ! Read in the other command line arguments
    call retrev('cmd_if', restart_from_file, len, ier)

    !> \todo Document these command line options
    slicenext    = iget('cmd_slice')
    collapsenext = iget('cmd_collapse')
    convert      = lget('cmd_convert') ! -convert present?

  end subroutine read_arguments

  !> Read in the user-supplied parameter-file
  subroutine read_par_files()
    use mod_global_parameters
    use mod_physics, only: physics_type
    use mod_small_values
    use mod_limiter
    use mod_slice

    logical          :: fileopen, file_exists
    integer          :: i, j, k, ifile, io_state
    integer          :: iB, isave, iw, level, idim, islice
    integer          :: nx_vec(^ND)
    integer          :: my_unit, iostate
    double precision :: dx_vec(^ND)

    character              :: c_ndim
    character(len=80)      :: fmt_string
    character(len=std_len) :: err_msg, restart_from_file_arg
    character(len=std_len) :: basename_full, basename_prev, dummy_file
    character(len=std_len), dimension(:), allocatable :: &
         typeboundary_min^D, typeboundary_max^D
    character(len=std_len), allocatable :: limiter(:)
    character(len=std_len), allocatable :: gradient_limiter(:)

    double precision, dimension(nsavehi) :: tsave_log, tsave_dat, tsave_slice, &
         tsave_collapsed, tsave_custom
    double precision :: dtsave_log, dtsave_dat, dtsave_slice, &
         dtsave_collapsed, dtsave_custom
    integer :: ditsave_log, ditsave_dat, ditsave_slice, &
         ditsave_collapsed, ditsave_custom
    integer :: windex

    namelist /filelist/ base_filename,restart_from_file, &
         typefilelog,firstprocess,reset_grid,snapshotnext, &
         convert,convert_type,saveprim,&
         nwauxio,nocartesian, w_write,writelevel,&
         writespshift,length_convert_factor, w_convert_factor, &
         time_convert_factor,level_io,level_io_min, level_io_max, &
         autoconvert,slice_type,slicenext,collapsenext,collapse_type

    namelist /savelist/ tsave,itsave,dtsave,ditsave,nslices,slicedir, &
         slicecoord,collapse,collapseLevel, time_between_print, &
         tsave_log, tsave_dat, tsave_slice, tsave_collapsed, tsave_custom, &
         dtsave_log, dtsave_dat, dtsave_slice, dtsave_collapsed, dtsave_custom, &
         ditsave_log, ditsave_dat, ditsave_slice, ditsave_collapsed, ditsave_custom

    namelist /stoplist/ it_init,time_init,it_max,time_max,dtmin,reset_it,reset_time

    namelist /methodlist/ time_integrator, &
         source_split_usr,typesourcesplit,&
         dimsplit,typedimsplit,&
         flux_scheme,typepred1,&
         limiter,gradient_limiter,cada3_radius,&
         loglimit,typelimited,typeboundspeed, &
         typetvd,typeentropy,entropycoef,typeaverage, &
         typegrad,typediv,&
         nxdiffusehllc, flathllc, tvdlfeps,&
         flatcd,flatsh,flatppm,&
         small_temperature,small_pressure,small_density, &
         small_values_method, small_values_daverage, check_small_values, &
         solve_internal_e, angmomfix

    namelist /boundlist/ nghostcells,typeboundary,typeghostfill,prolongation_method,&
         internalboundary, typeboundary_^L, save_physical_boundary

    namelist /meshlist/ refine_max_level,nbufferx^D,refine_threshold,&
         derefine_ratio, refine_criterion, stretched_grid, qst, &
         amr_wavefilter,max_blocks,block_nx^D,domain_nx^D,iprob,xprob^L, &
         w_refine_weight, prolongprimitive,coarsenprimitive, &
         typeprolonglimit, &
         logflag,tfixgrid,itfixgrid,ditregrid
    namelist /paramlist/  courantpar, dtpar, dtdiffpar, &
         typecourant, slowsteps

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

    ! defaults for boundary treatments
    typeghostfill = 'linear'

    ! default number of ghost-cell layers at each boundary of a block
    nghostcells = 2

    ! Allocate boundary conditions arrays in new and old style
    {
    allocate(typeboundary_min^D(nwfluxbc))
    allocate(typeboundary_max^D(nwfluxbc))
    typeboundary_min^D = undefined
    typeboundary_max^D = undefined
    }

    allocate(typeboundary(nwflux, 2 * ndim))
    typeboundary(:, :) = undefined

    ! not save physical boundary in dat files by default
    save_physical_boundary = .false.

    internalboundary   = .false.

    ! defaults for specific options
    typegrad   = 'central'
    typediv    = 'central'

    ! defaults for smallest physical values allowed
    small_temperature = 0.d0
    small_pressure    = 0.d0
    small_density     = 0.d0

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
    prolongation_method                = 'linear'
    coarsenprimitive            = .false.
    prolongprimitive            = .false.
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

    ! using stretched grid or not by default
    stretched_grid = .false.
    qst                         = bigdouble

    ! IO defaults
    it_init       = 0
    it_max        = biginteger
    time_init     = 0.d0
    time_max      = bigdouble
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

    typefilelog = 'default'

    ! defaults for input
    reset_time = .false.
    reset_it = .false.
    firstprocess  = .false.
    reset_grid     = .false.
    base_filename   = 'data'
    snapshotini = -1
    snapshotnext = -1

    ! Defaults for discretization methods
    typeaverage     = 'default'
    tvdlfeps        = one
    nxdiffusehllc   = 0
    flathllc        = .false.
    slowsteps       = -1
    courantpar      = 0.8d0
    typecourant     = 'maxsum'
    dimsplit        = .false.
    typedimsplit    = 'default'
    typelimited     = 'predictor'
    if(physics_type=='mhd') then
      cada3_radius  = 0.5d0
    else
      cada3_radius  = 0.1d0
    end if
    typetvd         = 'roe'
    typeboundspeed  = 'cmaxmean'
    source_split_usr= .false.
    time_integrator = 'twostep'
    solve_internal_e= .false.
    angmomfix       = .false.

    allocate(flux_scheme(nlevelshi),typepred1(nlevelshi))
    allocate(limiter(nlevelshi),gradient_limiter(nlevelshi))
    do level=1,nlevelshi
       flux_scheme(level) = 'tvdlf'
       typepred1(level)   = 'default'
       limiter(level)     = 'minmod'
       gradient_limiter(level) = 'minmod'
    end do

    flatcd          = .false.
    flatsh          = .false.
    flatppm         = .true.
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

107    close(unitpar)

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

    ! restart filename from command line overwrites the one in par file
    if(restart_from_file_arg /= undefined) &
      restart_from_file=restart_from_file_arg

    if (restart_from_file /= undefined) then
      ! Parse index in restart_from_file string (e.g. basename0000.dat)
      i = len_trim(restart_from_file) - 7
      read(restart_from_file(i:i+3), '(I4)', iostat=io_state) snapshotini
      if (io_state == 0 .and. snapshotnext==-1) snapshotnext = snapshotini + 1
      if(reset_time .or. reset_it) snapshotnext=0
    else
      snapshotnext=0
      if (firstprocess) &
           call mpistop("Please restart from a snapshot when firstprocess=T")
      if (convert) then
        convert = .false.
        write(uniterr,*) 'Change convert to .false. for a new run!'
      end if
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

    if (ditsave_log < bigdouble) ditsave(1) = ditsave_log
    if (ditsave_dat < bigdouble) ditsave(2) = ditsave_dat
    if (ditsave_slice < bigdouble) ditsave(3) = ditsave_slice
    if (ditsave_collapsed < bigdouble) ditsave(4) = ditsave_collapsed
    if (ditsave_custom < bigdouble) ditsave(5) = ditsave_custom

    if (mype == 0) then
       write(unitterm, *) ''
       write(unitterm, *) 'Output type | dtsave    | ditsave | itsave(1) | tsave(1)'
       write(fmt_string, *) '(A12," | ",E9.3E2," | ",I6,"  | "'//&
            ',I6, "    | ",E9.3E2)'
    end if

    do ifile = 1, nfile
       if (mype == 0) write(unitterm, fmt_string) trim(output_names(ifile)), &
            dtsave(ifile), ditsave(ifile), itsave(1, ifile), tsave(1, ifile)
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

    do level=1,nlevelshi
       if(flux_scheme(level)=='tvd'.and.time_integrator/='onestep') &
            call mpistop(" tvd is onestep method, reset time_integrator='onestep'")
       if(flux_scheme(level)=='tvd')then
          if(mype==0.and.(.not.dimsplit)) write(unitterm,*) &
               'Warning: setting dimsplit=T for tvd, as used for level=',level
          dimsplit=.true.
       endif
       if(flux_scheme(level)=='hlld'.and.physics_type/='mhd') &
          call mpistop("Cannot use hlld flux if not using MHD physics!")

       if (typepred1(level)=='default') then
          select case (flux_scheme(level))
          case ('cd')
             typepred1(level)='cd'
          case ('cd4')
             typepred1(level)='cd4'
          case ('fd')
             typepred1(level)='fd'
          case ('tvdlf','tvdmu')
             typepred1(level)='hancock'
          case ('hll')
             typepred1(level)='hll'
          case ('hllc')
             typepred1(level)='hllc'
          case ('hllcd')
             typepred1(level)='hllcd'
          case ('hlld')
             typepred1(level)='hlld'
          case ('nul','source','tvd')
             typepred1(level)='nul'
          case default
             call mpistop("No default predictor for this full step")
          end select
       end if
    end do

    ! finite difference scheme fd need global maximal speed
    if(any(flux_scheme=='fd')) need_global_cmax=.true.

    select case (time_integrator)
    case ("onestep")
       nstep=1
    case ("twostep", "twostep_trapezoidal")
       nstep=2
    case ("threestep")
       nstep=3
    case ("fourstep","rk4","jameson","ssprk43")
       nstep=4
    case ("ssprk54")
       nstep=5
    case default
       call mpistop("Unknown time_integrator")
    end select


    ! Harmonize the parameters for dimensional splitting and source splitting
    if(typedimsplit   =='default'.and.     dimsplit)   typedimsplit='xyyx'
    if(typedimsplit   =='default'.and..not.dimsplit)   typedimsplit='unsplit'
    dimsplit   = typedimsplit   /='unsplit'

    if(typeaxial=='default') then
      typeaxial='slab'
      if(mype==0) then
        write(*,*) 'Warning: coordinate system is not specified!'
        write(*,*) 'call set_coordinate_system in usr_init in mod_usr.t'
        write(*,*) 'Now use Cartesian coordinate'
      end if
    end if

    if(typeaxial=="slab") then
      if(stretched_grid) then
        typeaxial="slabstretch"
        slab=.false.
      else
        slab=.true.
      end if
    else
      slab=.false.
    end if

    if(typeaxial=='spherical') then
      if(dimsplit) then
        if(mype==0)print *,'Warning: spherical symmetry needs dimsplit=F, resetting'
        dimsplit=.false.
      end if
    end if

    if (ndim==1) dimsplit=.false.
    if (.not.dimsplit.and.ndim>1) then
       select case (time_integrator)
       case ("ssprk54","ssprk43","fourstep", "rk4", "threestep", "twostep")
          ! Runge-Kutta needs predictor
          typelimited="predictor"
          if (mype==0) write(unitterm, '(A30,A)') 'typelimited: ', typelimited
       end select
    end if

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
    if (any(typeboundary_min^D /= undefined)) then
      typeboundary(1:nwfluxbc, 2*^D-1) = typeboundary_min^D(1:nwfluxbc)
    end if

    if (any(typeboundary_max^D /= undefined)) then
      typeboundary(1:nwfluxbc, 2*^D) = typeboundary_max^D(1:nwfluxbc)
    end if
    }

    ! psi, tracers take the same boundary type as density
    !
    ! TODO (Jannis): we can later generalize this by letting physics module
    ! pre-define some boundary conditions in a fixed-size array (since we don't
    ! know nwflux yet), so that other variables can also copy the b.c. for
    ! momentum or energy equations.
    if (nwfluxbc<nwflux) then
      do iw = nwfluxbc+1, nwflux
        typeboundary(iw,:) = typeboundary(iw_rho, :)
      end do
    end if

    if (any(typeboundary == undefined)) then
      call mpistop("Not all boundary conditions have been defined")
    end if

    do idim=1,ndim
       periodB(idim)=(any(typeboundary(:,2*idim-1:2*idim)=='periodic'))
       aperiodB(idim)=(any(typeboundary(:,2*idim-1:2*idim)=='aperiodic'))
       if (periodB(idim).or.aperiodB(idim)) then
          do iw=1,nwflux
             if (typeboundary(iw,2*idim-1) .ne. typeboundary(iw,2*idim)) &
                  call mpistop("Wrong counterpart in periodic boundary")

             if (typeboundary(iw,2*idim-1) /= 'periodic' .and. &
                  typeboundary(iw,2*idim-1) /= 'aperiodic') then
               call mpistop("Each dimension should either have all "//&
                    "or no variables periodic, some can be aperiodic")
             end if
          end do
       end if
    end do
    {^NOONED
    do idim=1,ndim
      if(any(typeboundary(:,2*idim-1)=='pole')) then
        if(any(typeboundary(:,2*idim-1)/='pole')) typeboundary(:,2*idim-1)='pole'
        if(phys_energy) then
          windex=2
        else
          windex=1
        end if
        typeboundary(:,2*idim-1)='symm'
        if(physics_type/='rho') then
        select case(typeaxial)
        case('cylindrical')
          typeboundary(phi_+1,2*idim-1)='asymm'
          if(physics_type=='mhd') typeboundary(ndir+windex+phi_,2*idim-1)='asymm'
        case('spherical')
          typeboundary(3:ndir+1,2*idim-1)='asymm'
          if(physics_type=='mhd') typeboundary(ndir+windex+2:ndir+windex+ndir,2*idim-1)='asymm'
        case default
          call mpistop('Pole is in cylindrical, polar, spherical coordinates!')
        end select
        end if
      end if
      if(any(typeboundary(:,2*idim)=='pole')) then
        if(any(typeboundary(:,2*idim)/='pole')) typeboundary(:,2*idim)='pole'
        if(phys_energy) then
          windex=2
        else
          windex=1
        end if
        typeboundary(:,2*idim)='symm'
        if(physics_type/='rho') then
        select case(typeaxial)
        case('cylindrical')
          typeboundary(phi_+1,2*idim)='asymm'
          if(physics_type=='mhd') typeboundary(ndir+windex+phi_,2*idim)='asymm'
        case('spherical')
          typeboundary(3:ndir+1,2*idim)='asymm'
          if(physics_type=='mhd') typeboundary(ndir+windex+2:ndir+windex+ndir,2*idim)='asymm'
        case default
          call mpistop('Pole is in cylindrical, polar, spherical coordinates!')
        end select
        end if
      end if
    end do
    }

    if(any(limiter(1:nlevelshi)=='mp5')) then
      nghostcells=3
    end if

    if(any(limiter(1:nlevelshi)=='ppm')) then
      nghostcells=4
    end if

    select case (typeaxial)
       {^NOONED
    case ("spherical")
       xprob^LIM^DE=xprob^LIM^DE*two*dpi;
       \}
    case ("cylindrical")
       {
       if (^D==phi_) then
          xprob^LIM^D=xprob^LIM^D*two*dpi;
       end if
       \}
    end select

    ! full block size including ghostcells
    {ixGhi^D = block_nx^D + 2*nghostcells\}

    if(stretched_grid) then
      if(slab) call mpistop("Wrong geometry for stretched grid!")
      if(xprobmin1==0) call mpistop("Stretched grid needs xprobmin1 > 0")
      allocate(logGs(0:nlevelshi),qsts(0:nlevelshi))
      if (qst/=bigdouble) then
         xprobmax1=xprobmin1*qst**domain_nx1
         if(mype==0) write(*,*) 'xprobmax1 is computed for given domain_nx1 and qst:', xprobmax1
      else if (qst==bigdouble .and. xprobmax1/=bigdouble) then
         qst=(xprobmax1/xprobmin1)**(1.d0/dble(domain_nx1))
         logG=2.d0*(qst-1.d0)/(qst+1.d0)
         if(mype==0) write(*,*) 'logG and qst computed from xprobmax1: ', logG, qst
      end if
    end if

    nx_vec = [{domain_nx^D|, }]

    if (any(nx_vec < 2) .or. any(mod(nx_vec, 2) == 1)) &
         call mpistop('Grid size (domain_nx^D) has to be even and positive')

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

    if(refine_max_level>nlevelshi.or.refine_max_level<1)then
       write(unitterm,*)'Error: refine_max_level',refine_max_level,'>nlevelshi ',nlevelshi
       call mpistop("Reset nlevelshi and recompile!")
    endif

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
       write(unitterm, '(A30,I0)') 'slicenext: ', slicenext
       write(unitterm, '(A30,I0)') 'collapsenext: ', collapsenext
       write(unitterm, '(A30,A,A)')  'restart_from_file: ', ' ', trim(restart_from_file)
       write(unitterm, '(A30,L1)') 'converting: ', convert
       write(unitterm, '(A)') ''
    endif

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

  !> Standard method for creating a new output file
  subroutine create_output_file(fh, ix, extension)
    use mod_global_parameters
    integer, intent(out)         :: fh !< File handle
    integer, intent(in)          :: ix !< Index of file
    character(len=*), intent(in) :: extension !< Extension of file
    character(len=std_len)       :: filename
    integer :: amode

    if (ix >= 10000) then
      call mpistop("Number of output files is limited to 10000 (0...9999)")
    end if

    write(filename,"(a,i4.4,a)") trim(base_filename), ix, extension

    ! MPI cannot easily replace existing files
    open(unit=unitsnapshot,file=filename,status='replace')
    close(unitsnapshot, status='delete')

    amode = ior(MPI_MODE_CREATE, MPI_MODE_WRONLY)
    call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode, &
         MPI_INFO_NULL, fh, ierrmpi)

    if (ierrmpi /= 0) then
      print *, "Error, cannot create file ", trim(filename)
      call mpistop("Fatal error")
    end if

  end subroutine create_output_file

  !> Write header for a snapshot
  subroutine snapshot_write_header(fh, offset_tree, offset_block)
    use mod_forest
    use mod_physics
    use mod_global_parameters
    integer, intent(in)                       :: fh           !< File handle
    integer(kind=MPI_OFFSET_KIND), intent(in) :: offset_tree  !< Offset of tree info
    integer(kind=MPI_OFFSET_KIND), intent(in) :: offset_block !< Offset of block data
    integer, dimension(MPI_STATUS_SIZE)       :: st
    integer                                   :: iw, er

    call MPI_FILE_WRITE(fh, version_number, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, int(offset_tree), 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, int(offset_block), 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, nw, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, ndir, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, ndim, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, levmax, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, nleafs, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, nparents, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, it, 1, MPI_INTEGER, st, er)
    ! Note: It is nice when this double has an even number of 4 byte
    ! integers before it (for alignment)
    call MPI_FILE_WRITE(fh, global_time, 1, MPI_DOUBLE_PRECISION, st, er)

    call MPI_FILE_WRITE(fh, [ xprobmin^D ], ndim, &
         MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, [ xprobmax^D ], ndim, &
         MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, [ domain_nx^D ], ndim, &
         MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, [ block_nx^D ], ndim, &
         MPI_INTEGER, st, er)
    do iw = 1, nw
      call MPI_FILE_WRITE(fh, cons_wnames(iw), name_len, MPI_CHARACTER, st, er)
    end do

    ! TODO: write geometry info

    ! Physics related information
    call MPI_FILE_WRITE(fh, physics_type, name_len, MPI_CHARACTER, st, er)
    call phys_write_info(fh)

  end subroutine snapshot_write_header

  subroutine snapshot_read_header(fh, offset_tree, offset_block)
    use mod_forest
    use mod_global_parameters
    integer, intent(in)                   :: fh           !< File handle
    integer(MPI_OFFSET_KIND), intent(out) :: offset_tree  !< Offset of tree info
    integer(MPI_OFFSET_KIND), intent(out) :: offset_block !< Offset of block data
    integer                               :: ibuf(ndim), iw
    double precision                      :: rbuf(ndim)
    integer, dimension(MPI_STATUS_SIZE)   :: st
    character(len=10), allocatable        :: var_names(:)
    integer                               :: er

    ! Version number
    call MPI_FILE_READ(fh, ibuf(1), 1, MPI_INTEGER, st, er)
    if (all(compatible_versions /= ibuf(1))) then
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

    ! w_names (not used here)
    allocate(var_names(nw_found))
    do iw = 1, nw_found
      call MPI_FILE_READ(fh, var_names(iw), 10, MPI_CHARACTER, st, er)
    end do

  end subroutine snapshot_read_header

  !> Compute number of elements in index range
  pure integer function count_ix(ixO^L)
    integer, intent(in) :: ixO^L

    count_ix = product([ ixOmax^D ] - [ ixOmin^D ] + 1)
  end function count_ix

  !> Determine the shape of a block for output (whether to include ghost cells,
  !> and on which sides)
  subroutine block_shape_io(igrid, n_ghost, ixO^L, n_values)
    use mod_global_parameters

    integer, intent(in) :: igrid
    !> nghost(1:ndim) contains the number of ghost cells on the block's minimum
    !> boundaries, and nghost(ndim+1:2*ndim) on the block's maximum boundaries
    integer, intent(out) :: n_ghost(2*ndim)
    !> Index range on output block
    integer, intent(out) :: ixO^L
    !> Number of cells/values in output
    integer, intent(out) :: n_values

    integer            :: idim

    n_ghost(:) = 0

    if(save_physical_boundary) then
      do idim=1,ndim
        ! Include ghost cells on lower boundary
        if(pw(igrid)%is_physical_boundary(2*idim-1)) n_ghost(idim)=nghostcells
        ! Include ghost cells on upper boundary
        if(pw(igrid)%is_physical_boundary(2*idim)) n_ghost(ndim+idim)=nghostcells
      end do
    end if

    {ixOmin^D = ixMlo^D - n_ghost(^D)\}
    {ixOmax^D = ixMhi^D + n_ghost(ndim+^D)\}

    n_values = count_ix(ixO^L) * nw
  end subroutine block_shape_io

  subroutine write_snapshot
    use mod_forest
    use mod_global_parameters
    use mod_physics

    integer                       :: file_handle, igrid, Morton_no, iwrite
    integer                       :: ipe, ix_buffer(2*ndim+1), n_values
    integer                       :: ixO^L, n_ghost(2*ndim)
    integer, allocatable          :: block_ig(:, :), block_lvl(:)
    integer, allocatable          :: block_offset(:)
    integer                       :: iorecvstatus(MPI_STATUS_SIZE)
    integer                       :: ioastatus(MPI_STATUS_SIZE)
    integer                       :: igrecvstatus(MPI_STATUS_SIZE)
    integer                       :: istatus(MPI_STATUS_SIZE)
    type(tree_node), pointer      :: pnode
    integer(kind=MPI_OFFSET_KIND) :: offset_tree_info
    integer(kind=MPI_OFFSET_KIND) :: offset_block_data
    integer(kind=MPI_OFFSET_KIND) :: offset_offsets
    double precision, allocatable :: w_buffer(:)

    call MPI_BARRIER(icomm, ierrmpi)

    ! Allocate send/receive buffer
    n_values = count_ix(ixG^LL) * nw
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
           MPI_INTEGER, istatus, ierrmpi)

      call MPI_File_get_position(file_handle, offset_block_data, ierrmpi)

      ! Check whether data was written as expected
      if (offset_block_data - offset_tree_info /= &
           (nleafs + nparents) * size_logical + &
           nleafs * ((2+ndim) * size_int)) then
        call mpistop("Unexpected difference in offsets when writing .dat file")
      end if

      block_offset(1) = int(offset_block_data)
      iwrite = 0
    end if

    do Morton_no=Morton_start(mype), Morton_stop(mype)
      igrid  = sfc_to_igrid(Morton_no)
      itag   = Morton_no

      if (nwaux>0) then
        ! extra layer around mesh only for later averaging in convert
        ! set dxlevel value for use in gradient subroutine,
        ! which might be used in getaux
        saveigrid=igrid
        ^D&dxlevel(^D)=rnode(rpdx^D_, igrid);
        block=>pw(igrid)
        call phys_get_aux(.true., pw(igrid)%w, pw(igrid)%x, ixG^LL, &
             ixM^LL^LADD1, "write_snapshot")
      endif

      call block_shape_io(igrid, n_ghost, ixO^L, n_values)
      ix_buffer(1) = n_values
      ix_buffer(2:) = n_ghost
      w_buffer(1:n_values) = pack(pw(igrid)%w(ixO^S, 1:nw), .true.)

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
        block_offset(iwrite+1) = block_offset(iwrite) + n_values * size_double + &
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
          block_offset(iwrite+1) = block_offset(iwrite) + n_values * size_double + &
               2 * ndim * size_int
        end do
      end do

      ! Write block offsets (now we know them)
      call MPI_FILE_SEEK(file_handle, offset_offsets, MPI_SEEK_SET, ierrmpi)
      call MPI_FILE_WRITE(file_handle, block_offset(1:nleafs), nleafs, &
           MPI_INTEGER, istatus, ierrmpi)

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
    use mod_forest
    use mod_global_parameters

    integer                       :: ix_buffer(2*ndim+1), n_values
    integer                       :: ixO^L
    integer                       :: file_handle, amode, igrid, Morton_no, iread
    integer                       :: istatus(MPI_STATUS_SIZE)
    integer                       :: iorecvstatus(MPI_STATUS_SIZE)
    integer                       :: ipe,inrecv,nrecv, file_version
    logical                       :: fexist
    integer(MPI_OFFSET_KIND)      :: offset_tree_info
    integer(MPI_OFFSET_KIND)      :: offset_block_data
    double precision, allocatable :: w_buffer(:)
    double precision, dimension(:^D&,:), allocatable :: w

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

    ! Allocate send/receive buffer
    n_values = count_ix(ixG^LL) * nw_found
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
          {ixOmin^D = ixMlo^D + ix_buffer(^D)\}
          {ixOmax^D = ixMhi^D + ix_buffer(ndim+^D)\}
          n_values = count_ix(ixO^L) * nw_found

          call MPI_FILE_READ(file_handle, w_buffer, n_values, &
               MPI_DOUBLE_PRECISION, istatus, ierrmpi)

          if (mype == ipe) then ! Root task
            igrid=sfc_to_igrid(Morton_no)
            w(ixO^S, 1:nw_found) = reshape(w_buffer(1:n_values), &
                 shape(w(ixO^S, 1:nw_found)))
            if (nw_found<nw) then
              if (associated(usr_transform_w)) then
                call usr_transform_w(ixG^LL,ixM^LL,nw_found,w,pw(igrid)%x,pw(igrid)%w)
              else
                pw(igrid)%w(ixO^S,1:nw_found)=w(ixO^S,1:nw_found)
              end if
            else if (nw_found>nw) then
              if (associated(usr_transform_w)) then
                call usr_transform_w(ixG^LL,ixM^LL,nw_found,w,pw(igrid)%x,pw(igrid)%w)
              else
                pw(igrid)%w(ixO^S,1:nw)=w(ixO^S,1:nw)
              end if
            else
              pw(igrid)%w(ixO^S,1:nw)=w(ixO^S,1:nw)
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
        itag=Morton_no

        call MPI_RECV(ix_buffer, 2*ndim+1, MPI_INTEGER, 0, itag, icomm,&
             iorecvstatus, ierrmpi)
        {ixOmin^D = ix_buffer(^D)\}
        {ixOmax^D = ix_buffer(ndim+^D)\}
        n_values = ix_buffer(2*ndim+1)

        call MPI_RECV(w_buffer, n_values, MPI_DOUBLE_PRECISION,&
             0, itag, icomm, iorecvstatus, ierrmpi)

        w(ixO^S, 1:nw_found) = reshape(w_buffer(1:n_values), &
             shape(w(ixO^S, 1:nw_found)))
        if (nw_found<nw) then
          if (associated(usr_transform_w)) then
            call usr_transform_w(ixG^LL,ixM^LL,nw_found,w,pw(igrid)%x,pw(igrid)%w)
          else
            pw(igrid)%w(ixO^S,1:nw_found)=w(ixO^S,1:nw_found)
          end if
        else if (nw_found>nw) then
          if (associated(usr_transform_w)) then
            call usr_transform_w(ixG^LL,ixM^LL,nw_found,w,pw(igrid)%x,pw(igrid)%w)
          else
            pw(igrid)%w(ixO^S,1:nw)=w(ixO^S,1:nw)
          end if
        else
          pw(igrid)%w(ixO^S,1:nw)=w(ixO^S,1:nw)
        end if
      end do
    end if

    call MPI_BARRIER(icomm,ierrmpi)

  end subroutine read_snapshot

  subroutine read_snapshot_old()
    use mod_forest
    use mod_global_parameters

    double precision              :: wio(ixG^T,1:nw)
    integer                       :: fh, igrid, Morton_no, iread
    integer                       :: levmaxini, ndimini, ndirini
    integer                       :: nwini, neqparini, nxini^D
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer                       :: istatus(MPI_STATUS_SIZE)
    integer, allocatable          :: iorecvstatus(:,:)
    integer                       :: ipe,inrecv,nrecv
    integer                       :: sendini(7+^ND)
    character(len=80)             :: filename
    logical                       :: fexist
    double precision              :: eqpar_dummy(100)

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
        call MPI_FILE_READ_AT(fh,offset,pw(igrid)%w,1,type_block_io, &
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
        call MPI_RECV(pw(igrid)%w,1,type_block_io,0,itag,icomm,&
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

    logical              :: fileopen
    integer              :: i, iw, level
    double precision     :: wmean(1:nw), total_volume
    double precision     :: volume_coverage(refine_max_level)
    integer              :: nx^D, nc, ncells, dit
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

          call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'), &
               len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
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

    integer, parameter :: n_modes = 2
    integer, parameter :: my_unit = 123
    character(len=40)  :: fmt_string
    logical, save      :: file_open = .false.
    integer            :: power
    double precision   :: modes(nw, n_modes), volume

    do power = 1, n_modes
       call get_volume_average(power, modes(:, power), volume)
    end do

    if (mype == 0) then
       if (.not. file_open) then
          open(my_unit, file = trim(base_filename) // ".log")
          file_open = .true.

          write(my_unit, *) "# time mean(w) mean(w**2)"
       end if

       write(fmt_string, "(a,i0,a)") "(", nw * n_modes + 1, fmt_r // ")"
       write(my_unit, fmt_string) global_time, modes
    end if
  end subroutine printlog_regression_test

  !> Compute mean(w**power) over the leaves of the grid. The first mode
  !> (power=1) corresponds to to the mean, the second to the mean squared values
  !> and so on.
  subroutine get_volume_average(power, mode, volume)
    use mod_global_parameters

    integer, intent(in)           :: power     !< Which mode to compute
    double precision, intent(out) :: mode(nw)  !< The computed mode
    double precision, intent(out) :: volume    !< The total grid volume
    integer                       :: iigrid, igrid, iw
    double precision              :: wsum(nw+1)
    double precision              :: dvolume(ixG^T)
    double precision              :: dsum_recv(1:nw+1)

    wsum(:) = 0

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! Determine the volume of the grid cells
       if (slab) then
          dvolume(ixM^T) = {rnode(rpdx^D_,igrid)|*}
       else
          dvolume(ixM^T) = pw(igrid)%dvolume(ixM^T)
       end if

       ! Store total volume in last element
       wsum(nw+1) = wsum(nw+1) + sum(dvolume(ixM^T))

       ! Compute the modes of the cell-centered variables, weighted by volume
       do iw = 1, nw
          wsum(iw) = wsum(iw) + &
               sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,iw)**power)
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
    integer                       :: iigrid, igrid, i^D
    double precision              :: wsum(2)
    double precision              :: dvolume(ixG^T)
    double precision              :: dsum_recv(2)

    wsum(:) = 0

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)

       ! Determine the volume of the grid cells
       if (slab) then
          dvolume(ixM^T) = {rnode(rpdx^D_,igrid)|*}
       else
          dvolume(ixM^T) = pw(igrid)%dvolume(ixM^T)
       end if

       ! Store total volume in last element
       wsum(2) = wsum(2) + sum(dvolume(ixM^T))

       ! Compute the modes of the cell-centered variables, weighted by volume
       {do i^D = ixMlo^D, ixMhi^D\}
       wsum(1) = wsum(1) + dvolume(i^D) * &
            func(pw(igrid)%w(i^D, :), nw)
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

    integer                       :: iigrid, igrid, iw
    double precision              :: wmax_mype(nw),wmax_recv(nw)

    wmax_mype(1:nw) = -bigdouble

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       do iw = 1, nw
          wmax_mype(iw)=max(wmax_mype(iw),maxval(pw(igrid)%w(ixM^T,iw)))
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

    integer                       :: iigrid, igrid, iw
    double precision              :: wmin_mype(nw),wmin_recv(nw)

    wmin_mype(1:nw) = bigdouble

    ! Loop over all the grids
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       do iw = 1, nw
          wmin_mype(iw)=min(wmin_mype(iw),minval(pw(igrid)%w(ixM^T,iw)))
       end do
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(wmin_mype, wmin_recv, nw, MPI_DOUBLE_PRECISION, &
         MPI_MIN, icomm, ierrmpi)

    wmin(1:nw)=wmin_recv(1:nw)

  end subroutine get_global_minima



end module mod_input_output
