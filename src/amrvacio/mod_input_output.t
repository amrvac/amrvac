!> Module for reading input and writing output
module mod_input_output

  implicit none
  public

  ! Formats used in output
  character(len=*), parameter :: fmt_r  = 'es16.8' ! Default precision
  character(len=*), parameter :: fmt_r2 = 'es10.2' ! Two digits
  character(len=*), parameter :: fmt_i  = 'i8'     ! Integer format

contains

  !> Read the command line arguments passed to amrvac
  subroutine read_arguments()
    use M_kracken
    use mod_global_parameters

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
    call kracken('cmd','-i amrvac.par -if data -restart -1 -slice 0'//&
         ' -collapse 0 --help .false. -convert .false.')

    ! Get the par file(s)
    call retrev('cmd_i', all_par_files, len, ier)

    ! Show the usage if the help flag was given, or no par file was specified
    if (lget('cmd_-help') .or. len == 0) then
       if (mype == 0) then
          print *, 'Usage example:'
          print *, './amrvac -i file.par [file2.par ...]'
          print *, ''
          print *, 'Optional arguments:'
          print *, '-restart <N>         Restart a run at this snapshot'
          print *, '-convert             Convert snapshot files'
          print *, '-if file.dat         Use this snapshot'
          print *, ''
          print *, 'Note: later parameter files override earlier ones.'
       end if
       call comm_finalize()
       stop
    end if

    ! Split the par files, in case multiple were given
    call delim(all_par_files, tmp_files, max_files, n_par_files, &
         ibegin, iterm, len, " ,'"""//char(9))
    allocate(par_files(n_par_files))
    par_files = tmp_files(1:n_par_files)

    ! Read in the other command line arguments
    call retrev('cmd_if', filenameini, len, ier)
    snapshotini  = iget('cmd_restart')
    snapshotnext = snapshotini+1
    !> \todo Document these command line options
    slicenext    = iget('cmd_slice')
    collapseNext = iget('cmd_collapse')
    convert      = lget('cmd_convert') ! -convert present?

  end subroutine read_arguments

  !> Read in the user-supplied parameter-file
  subroutine read_par_files()
    use mod_global_parameters
    use mod_physics, only: physics_type

    logical          :: fileopen, file_exists
    integer          :: i, j, k, ifile, io_state
    integer          :: iB, isave, iw, level, idim, islice
    integer          :: nx_vec(^ND)
    double precision :: dx_vec(^ND)

    character              :: c_ndim
    character(len=80)      :: fmt_string
    character(len=std_len) :: err_msg
    character(len=std_len) :: filenameout_full, filenameout_prev
    character(len=std_len) :: filenamelog_full, filenamelog_prev

    namelist /filelist/ filenameout,filenameini,filenamelog, &
         snapshotini,typefilelog,firstprocess,resetgrid,snapshotnext, &
         convert,convert_type,dxfiletype,saveprim,primnames, &
         typeparIO,uselimiter,nwauxio,nocartesian,addmpibarrier, &
         writew,writelevel,writespshift,endian_swap, &
         normvar,normt,level_io,level_io_min,level_io_max, &
         autoconvert,sliceascii,slicenext,collapseNext,collapse_type
    namelist /savelist/ tsave,itsave,dtsave,ditsave,nslices,slicedir, &
         slicecoord,collapse,collapseLevel{#IFDEF MAGNETOFRICTION , ditsavemf}
    namelist /stoplist/ itmax,tmax,tmaxexact,dtmin,t,it,treset,itreset,residmin,&
         residmax,typeresid{#IFDEF MAGNETOFRICTION , itmaxmf}
    namelist /methodlist/ wnames,fileheadout,typeadvance, &
         ssplitdust,ssplitdivb,ssplitresis,ssplituser,typesourcesplit,&
         ncyclemax,&
         dimsplit,typedimsplit,typeaxial,typecoord,&
         typefull1,typepred1,&
         typelimiter1,mcbeta,typegradlimiter1,&
         flatcd,flatsh,flatppm,&
         loglimit,typelimited,useprimitive,typetvdlf, &
         typetvd,typeentropy,entropycoef,typeaverage, &
         B0field,Bdip,Bquad,Boct,Busr,divbwave,&
         typedivbfix,divbdiff,typedivbdiff,compactres,&
         useprimitiveRel,maxitnr,dmaxvel,typepoly,&
         tvdlfeps,BnormLF,&
         smallT,smallp,smallrho,typegrad,typediv,tolernr,absaccnr,&
         strictnr,fixsmall,strictsmall,strictgetaux,nflatgetaux,&
         strictzero,nxdiffusehllc,typespherical,&
         fixprocess,flathllc, &
         ncool,cmulti,coolmethod,coolcurve,Tfix, &
         smallrhod, dustzero, dustmethod,dustspecies,dusttemp, &
         x1ptms,x2ptms,x3ptms,ptmass,tlow,nwtf,neqpartf
    namelist /boundlist/ dixB,typeB,typeghostfill,typegridfill,ratebdflux,&
         internalboundary
    namelist /amrlist/ mxnest,nbufferx^D,specialtol,tol,tolratio,errorestimate, &
         amr_wavefilter,ngridshi,nxblock^D,nxlone^D,iprob,xprob^L, &
         skipfinestep,wflags,flags,&
         restrictprimitive,prolongprimitive,coarsenprimitive, &
         typeprolonglimit, &
         amrentropy,logflag,tfixgrid,itfixgrid,ditregrid{#IFDEF STRETCHGRID ,qst}
    namelist /paramlist/  courantpar, dtpar, dtdiffpar, dtTCpar,&
         typecourant, slowsteps, cfrac{#IFDEF MAGNETOFRICTION , cmf_c, cmf_y, cmf_divb}
    !----------------------------------------------------------------------------
 
    ! default maximum number of grid blocks in a processor
    ngridshi=4000

    ! default block size excluding ghost cells
    {nxblock^D = 12\}

    ! defaults for boundary treatments
    ratebdflux         = one
    typeghostfill      = 'linear'
    dixB               = 2
    ! default block size
    {ixGhi^D = nxblock^D + 2*dixB\}
    allocate(typeB(nw, nhiB))
    typeB(1:nw,1:nhiB) = 'cont'
    internalboundary   = .false.

    ! code behavior and testing defaults
    addmpibarrier = .false.

    ! defaults for parameters for optional cooling module (van Marle)
    ncool      = 100
    cmulti     = 1
    coolcurve  = 'DM'
    coolmethod = 'explicit2'
    cfrac      = 0.1d0
    Tfix       = .false.

    ! defaults for dust
    dustzero    = .false.
    dustmethod  = 'Kwok'
    dustspecies = 'graphite'
    dusttemp    = 'constant'

    ! defaults for parameters for optional pointgrav module (van Marle)
    ! --> set here mass to zero, coordinates to zero
    x1ptms = zero
    x2ptms = zero
    x3ptms = zero
    ptmass = zero

    ! defaults for specific options
    fixprocess = .false.
    typegrad   = 'central'
    typediv    = 'central'
    smallT     = -one
    smallp     = -one
    smallrho   = -one
    smallrhod  = -one

    ! relativistic module defaults
    useprimitiveRel = .true.
    strictnr        = .true.
    fixsmall        = .false.
    strictsmall     = .true.
    strictgetaux    = .false.
    nflatgetaux     = 1
    strictzero      = .true.
    typepoly        = 'meliani'
    tolernr         = 1.0d-13
    absaccnr        = 1.0d-13
    maxitnr         = 100
    dmaxvel         = 1.0d-8
    tlow            = zero

    ! defaults for convert behavior

    nwauxio                  = 0
    nocartesian              = .false.
    saveprim                 = .false.
    autoconvert              = .false.
    endian_swap              = .false.
    convert_type             = 'vtuBCCmpi'
    collapse_type            = 'vti'
    dxfiletype               = 'lsb'
    allocate(writew(nw))
    writew(1:nw)             = .true.
    writelevel(1:nlevelshi)  = .true.
    writespshift(1:ndim,1:2) = zero
    level_io                 = -1
    level_io_min             = 1
    level_io_max             = nlevelshi

    ! normalization of primitive variables: only for output
    ! note that normvar(0) is for length
    ! this scaling is optional, and must be set consistently if used
    allocate(normvar(0:nw))
    normvar(0:nw) = one
    normt         = one

    ! residual defaults
    residmin  = -1.0d0
    residmax  = bigdouble
    typeresid = 'relative'

    ! AMR related defaults
    mxnest                      = 1
    {nbufferx^D                 = 0\}
    specialtol                  = .false.
    tol(1:nlevelshi)            = 0.1d0
    tolratio(1:nlevelshi)       = 1.0d0/8.0d0
    typegridfill                = 'linear'
    amrentropy                  = .false.
    restrictprimitive           = .false.
    coarsenprimitive            = .false.
    prolongprimitive            = .false.
    typeprolonglimit            = 'default'
    errorestimate               = 3
    allocate(flags(nflag_))
    allocate(wflags(nflag_))
    flags(1:nflag_)             = 0
    wflags(1:nflag_)            = zero
    flags(nflag_)               = 1
    flags(1)                    = 1
    wflags(1)                   = one
    allocate(logflag(nw))
    logflag(1:nw)               = .false.
    amr_wavefilter(1:nlevelshi) = 1.0d-2
    skipfinestep                = .false.
    tfixgrid                    = bigdouble
    itfixgrid                   = biginteger
    ditregrid                   = 1
    {#IFDEF STRETCHGRID
    qst                         = bigdouble
    }

    ! MHD specific defaults
    B0field      = .false.
    Bdip         = zero
    Bquad        = zero
    Boct         = zero
    Busr         = zero
    {#IFNDEF GLM
    typedivbfix  = 'linde'\}
    {#IFDEF GLM
    typedivbfix  = 'glm1'\}
    divbdiff     = 0.5d0
    typedivbdiff = 'all'
    compactres   = .false.
    divbwave     = .true.


    ! IO defaults
    itmax         = biginteger
    {#IFDEF MAGNETOFRICTION
    itmaxmf       = 0
    ditsavemf     = 20000
    }
    tmax          = bigdouble
    tmaxexact     = .true.
    dtmin         = 1.0d-10
    typeparIO     = 0
    nslices       = 0
    collapse      = .false.
    collapseLevel = 1
    sliceascii    = .false.

    do ifile=1,nfile
       do isave=1,nsavehi
          tsave(isave,ifile)  = bigdouble  ! t  of saves into the output files
          itsave(isave,ifile) = biginteger ! it of saves into the output files
       end do
       dtsave(ifile)  = bigdouble  ! time between saves
       ditsave(ifile) = biginteger ! timesteps between saves
       isavet(ifile)  = 1          ! index for saves by t
       isaveit(ifile) = 1          ! index for saves by it
    end do

    typefilelog = 'default'
    fileheadout = 'AMRVAC'
    ! defaults for number of w in the transformed data
    nwtf        = 0
    ! defaults for number of equation parameters in the transformed data
    neqpartf    = 0

    ! defaults for input 
    firstprocess  = .false.
    resetgrid     = .false.
    treset        = .false.
    itreset       = .false.
    filenameout   = 'data'
    filenamelog   = 'amrvac'

    ! Defaults for discretization methods
    typeaverage     = 'default'
    tvdlfeps        = one
    {#IFDEF FCT
    BnormLF         = .false.
    }
    {#IFNDEF FCT
    BnormLF         = .true.
    }
    nxdiffusehllc   = 0
    flathllc        = .false.
    typeaxial       = 'slab'
    typecoord       = 'default'
    typespherical   = 1
    slowsteps       = -1
    courantpar      = 0.8d0
    typecourant     = 'maxsum'
    dimsplit        = .false.
    typedimsplit    = 'default'
    typelimited     = 'predictor'
    mcbeta          = 1.4d0
    useprimitive    = .true.
    typetvd         = 'roe'
    typetvdlf       = 'cmaxmean'
    bcphys          = .true.
    ncyclemax       = 1000
    ssplitdust      = .false.
    ssplitdivb      = .false.
    {^IFMHDPHYS
    ssplitdivb      = .true.
    }
    ssplitresis     = .false.
    ssplituser      = .false.
    typeadvance     = 'twostep'

    do level=1,nlevelshi
       typefull1(level)        = 'tvdlf'
       typepred1(level)        = 'default'
       typelimiter1(level)     = 'minmod'
       typegradlimiter1(level) = 'minmod'
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
    dtTCpar       = 0.9d0/dble(ndim)
    dtpar         = -1.d0
    {#IFDEF MAGNETOFRICTION
    cmf_c         = 0.3d0
    cmf_y         = 0.2d0
    cmf_divb      = 0.1d0
    }

    ! problem setup defaults
    {nxlone^D = 0\}
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

    ! Set default variable names
    primnames = 'default'
    wnames    = 'default'

    ! These are used to construct file and log names from multiple par files
    filenameout_full = ''
    filenameout_prev = ''
    filenamelog_full = ''
    filenamelog_prev = ''

    ! Read in the .par files
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
       read(unitpar, amrlist, end=106)

106    rewind(unitpar)
       read(unitpar, paramlist, end=107)

107    close(unitpar)

       ! Append the log and file names given in the par files
       if (filenameout /= filenameout_prev) &
            filenameout_full = trim(filenameout_full) // trim(filenameout)
       filenameout_prev = filenameout

       if (filenamelog /= filenamelog_prev) &
            filenamelog_full = trim(filenamelog_full) // trim(filenamelog)
       filenamelog_prev = filenamelog
    end do

    filenameout = filenameout_full
    filenamelog = filenamelog_full

    if(TRIM(primnames)=='default'.and.mype==0) write(uniterr,*) &
         'Warning in read_par_files: primnames not given!'

    if(firstprocess .and. snapshotini<0) &
         call mpistop("Please restart from a snapshot when firstprocess=T")
    if(convert) autoconvert=.false.

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

    if(itmax==biginteger .and. tmax==bigdouble.and.mype==0) write(uniterr,*) &
         'Warning in read_par_files: itmax or tmax not given!'

    if(residmin>=zero) then
       if (mype==0) write(unitterm,*)"SS computation with input value residmin"
       if(residmin<=smalldouble) call mpistop("Provide value for residual above smalldouble")
    end if

    if(compactres)then
       if(typeaxial/='slab') call mpistop("compactres for MHD only in cartesian case")
    endif

    if(TRIM(wnames)=='default') call mpistop("Provide wnames and restart code")

    do level=1,nlevelshi
       !if(typefull1(level)=='tvdlf1'.and.typeadvance=='twostep') &
       !   call mpistop(" tvdlf1 is onestep method, reset typeadvance=onestep!")
       !if(typefull1(level)=='hll1'.and.typeadvance=='twostep') &
       !   call mpistop(" hll1 is onestep method, reset typeadvance=onestep!")
       !if(typefull1(level)=='hllc1'.and.typeadvance=='twostep') &
       !   call mpistop(" hllc1 is onestep method, reset typeadvance=onestep!")
       !if(typefull1(level)=='hllcd1'.and.typeadvance=='twostep') &
       !   call mpistop(" hllcd1 is onestep method, reset typeadvance=onestep!")
       !if(typefull1(level)=='tvdmu1'.and.typeadvance=='twostep') &
       !   call mpistop(" tvdmu1 is onestep method, reset typeadvance=onestep!")
       if(typefull1(level)=='tvd'.and.typeadvance=='twostep') &
            call mpistop(" tvd is onestep method, reset typeadvance=onestep!")
       if(typefull1(level)=='tvd1'.and.typeadvance=='twostep') &
            call mpistop(" tvd1 is onestep method, reset typeadvance=onestep!")
       if(typefull1(level)=='tvd'.or.typefull1(level)=='tvd1')then 
          if(mype==0.and.(.not.dimsplit)) write(unitterm,*) &
               'Warning: setting dimsplit=T for tvd, as used for level=',level
          dimsplit=.true.
       endif

       if (typepred1(level)=='default') then
          select case (typefull1(level))
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
          case ('hlldd')
             typepred1(level)='hlldd'
          case ('tvdlf1','tvdmu1','tvd1','tvd','hll1','hllc1', &
               'hlld1','hllcd1','hlldd1','nul','source')
             typepred1(level)='nul'
          case default
             call mpistop("No default predictor for this full step")
          end select
       end if
    end do

    select case (typeadvance)
    case ("onestep")
       nstep=1
    case ("twostep")
       nstep=2
    case ("threestep")
       nstep=3
    case ("fourstep","rk4","jameson","ssprk43")
       nstep=4
    case ("ssprk54")
       nstep=5
    case default
       call mpistop("Unknown typeadvance")
    end select


    ! Harmonize the parameters for dimensional splitting and source splitting
    if(typedimsplit   =='default'.and.     dimsplit)   typedimsplit='xyyx'
    if(typedimsplit   =='default'.and..not.dimsplit)   typedimsplit='unsplit'
    dimsplit   = typedimsplit   /='unsplit'


    if (typeaxial=="slab") then
       slab=.true.
    else
       slab=.false.
    end if

    if (typeaxial=='spherical') then
       if (dimsplit) then
          if(mype==0)print *,'Warning: spherical symmetry needs dimsplit=F, resetting'
          dimsplit=.false.
       end if
    end if

    if (typecoord=='default') then
       typecoord = typeaxial
    end if

    if (ndim==1) dimsplit=.false.
    if (.not.dimsplit.and.ndim>1) then
       select case (typeadvance)
       case ("ssprk54","ssprk43","fourstep", "rk4", "threestep", "twostep")
          ! Runge-Kutta needs predictor
          typelimited="predictor"
          if (mype==0) write(unitterm, '(A30,A)') 'typelimited: ', 'predictor (for RK)'
       end select
    end if

    {#IFDEF GLM
    if (typedivbfix/='glm1' .and. typedivbfix/='glm2' .and. typedivbfix/='glm3') &
         call mpistop('using GLM, so typedivbfix should be either glm1 or glm2 (or glm3 with no additional sources)')
    if (ssplitdivb .eqv. .false.) &
         call mpistop('GLM needs ssplitdivb = .true.')
    \}
    {#IFNDEF GLM
    if(typedivbfix=='glm1' .or. typedivbfix=='glm2' .or. typedivbfix=='glm3') then
       call mpistop('Not using GLM, to enable GLM, change definitions.h. '//&
            'Alternatives for typedivbfix are powel, janhunen, linde or none')
    end if
    \}
    !if (B0field) then
    !   if(mype==0)print *,'B0+B1 split for MHD'
    !   if (.not. physics_type=='mhd') call mpistop("B0+B1 split for MHD only")
    !end if

    !if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    !     .and.(flatsh.and.physics_type=='rho')) then
    !   call mpistop(" PPM with flatsh=.true. can not be used with physics_type='rho'!")
    !end if
    !if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    !     .and.(flatsh.and.physics_type=='hdadiab')) then
    !   call mpistop(" PPM with flatsh=.true. can not be used with physics_type='hdadiab'!")
    !end if
    !if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    !     .and.(flatcd.and.physics_type=='hdadiab')) then
    !   call mpistop(" PPM with flatcd=.true. can not be used with physics_type='hdadiab'!")
    !end if
    !if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    !     .and.(flatsh.and..not.useprimitive)) then
    !   call mpistop(" PPM with flatsh=.true. needs useprimitive=T!")
    !end if
    !if (any(typelimiter1(1:nlevelshi)== 'ppm')&
    !     .and.(flatcd.and..not.useprimitive)) then
    !   call mpistop(" PPM with flatcd=.true. needs useprimitive=T!")
    !end if

    do idim=1,ndim
       periodB(idim)=(any(typeB(:,2*idim-1:2*idim)=='periodic'))
       aperiodB(idim)=(any(typeB(:,2*idim-1:2*idim)=='aperiodic'))
       if (periodB(idim).or.aperiodB(idim)) then
          do iw=1,nw
             if (typeB(iw,2*idim-1) .ne. typeB(iw,2*idim)) &
                  call mpistop("Wrong counterpart in periodic boundary")
             if (typeB(iw,2*idim-1) /= 'periodic' .and. typeB(iw,2*idim-1) /= 'aperiodic') &
                  call mpistop("Each dimension should either have all &
                  or no variables periodic, some can be aperiodic")
          end do
       end if
    end do

    if (any(typelimiter1(1:nlevelshi)=='ppm').and.(dixB<4)) then
       call mpistop(" PPM works only with dixB>=4 !")
    end if

    if (any(typelimiter1(1:nlevelshi)=='mp5') .and. (dixB<3)) then
       call mpistop("mp5 needs at at least 3 ghost cells! Set dixB=3 in boundlist.")
    end if

    select case (typeaxial)
       {^NOONED
    case ("spherical")
       xprob^LIM^DE=xprob^LIM^DE*two*dpi;
       \}
    case ("cylindrical")
       {
       if (^D==^PHI) then
          xprob^LIM^D=xprob^LIM^D*two*dpi;
       end if
       \}
    end select

    {#IFDEF STRETCHGRID
    !if (mxnest>1) call mpistop("No refinement possible with a loggrid")
    if (typeaxial=='slab') call mpistop("Cartesian log grid not implemented")

    if (qst/=bigdouble) then
       xprobmax1=xprobmin1*qst**nxlone1
       if(mype==0) write(*,*) 'xprobmax1 is computed for given nxlone1 and qst:', xprobmax1
    else if (qst==bigdouble .and. xprobmax1/=bigdouble) then
       qst=(xprobmax1/xprobmin1)**(1.d0/dble(nxlone1))
       logG=2.d0*(qst-1.d0)/(qst+1.d0)
       if(mype==0) write(*,*) 'logG and qst computed from xprobmax1: ', logG, qst
    end if
    }

    nx_vec = [{nxlone^D|, }]

    if (any(nx_vec < 2) .or. any(mod(nx_vec, 2) == 1)) &
         call mpistop('Grid size (nxlone^D) has to be even and positive')

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

    if(mxnest>nlevelshi.or.mxnest<1)then
       write(unitterm,*)'Error: mxnest',mxnest,'>nlevelshi ',nlevelshi
       call mpistop("Reset nlevelshi and recompile!")
    endif

    if (flags(nflag_)>nw) then
       write(unitterm,*)'Error: flags(nw+1)=',flags(nw+1),'>nw ',nw
       call mpistop("Reset flags(nw+1)!")
    end if
    if (flags(nflag_)==0) errorestimate=0
    if (flags(nflag_)<0) then
       if (mype==0) then
          write(unitterm,*) "flags(",nflag_,") can not be negative"
          call mpistop("")
       end if
    end if

    if (mype==0) write(unitterm, '(A30)', advance='no') 'Error estimation: '

    select case (errorestimate)
    case (0)
       if (mype==0) write(unitterm, '(A)') "user defined"
    case (2)
       if (mype==0) write(unitterm, '(A)') "relative error"
    case (3)
       if (mype==0) write(unitterm, '(A)') "Lohner's scheme"
    case (4)
       if (mype==0) write(unitterm, '(A)') "Lohner's original scheme"
    case default
       call mpistop("Unknown error estimator, change errorestimate")
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

    ! Warn when too few blocks at start of simulation 
    if (mype.eq.0 .and. snapshotini.eq.-1 .and. {^D& floor(dble(nxlone^D)/dble(nxblock^D)) |*} .lt. npe) then
       call mpistop('Need at least as many blocks on level 1 as cores to initialize!')
    end if


    if (mype==0) then
       write(unitterm, '(A30,I0)') 'snapshotini: ', snapshotini
       write(unitterm, '(A30,I0)') 'slicenext: ', slicenext
       write(unitterm, '(A30,I0)') 'collapsenext: ', collapsenext
       write(unitterm, '(A30,A,A)')  'filenameini: ', ' ', trim(filenameini)
       write(unitterm, '(A30,L1)') 'converting: ', convert
       write(unitterm, '(A)') ''
    endif

  end subroutine read_par_files

  subroutine saveamrfile(ifile)

    ! following specific for Intel compiler and use on VIC3 with MPT
    !DEC$ ATTRIBUTES NOINLINE :: write_snapshot
    use mod_usr_methods, only: usr_print_log
    use mod_global_parameters
    integer:: ifile
    !-----------------------------------------------------------------------------
    select case (ifile)
    case (fileout_)
       if(endian_swap) typeparIO=-1
       if (typeparIO==1)then
          call write_snapshot
       else if(typeparIO==0) then
          call write_snapshot_nopar
       else if(typeparIO==-1) then
          call write_snapshot_noparf
       endif
       !opedit: now we can also convert directly and will when autoconvert is set in inifile: 
       if (autoconvert) call generate_plotfile
       {#IFDEF PARTICLES
       call write_particles_snapshot
       }
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
       call write_analysis
    case default
       write(*,*) 'No save method is defined for ifile=',ifile
       call mpistop("")
    end select

    ! opedit: Flush stdout and stderr from time to time.
    flush(unit=unitterm)

  end subroutine saveamrfile

  subroutine write_snapshot
    use mod_forest
    use mod_global_parameters
    use mod_physics

    integer :: file_handle, amode, igrid, Morton_no, iwrite
    integer :: nx^D
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer, dimension(MPI_STATUS_SIZE) :: istatus
    character(len=80) :: filename, line
    logical, save :: firstsnapshot=.true.
    !-----------------------------------------------------------------------------
    if (firstsnapshot) then
       snapshot=snapshotnext
       firstsnapshot=.false.
    end if

    if (snapshot >= 10000) then
       if (mype==0) then
          write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
          write(*,*) "overwriting first frames"
       end if
       snapshot=0
    end if

    ! generate filename
    write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"

    if(mype==0) then
       open(unit=unitsnapshot,file=filename,status='replace')
       close(unit=unitsnapshot, status='delete')
    end if
    call MPI_BARRIER(icomm,ierrmpi)

    amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
    call MPI_FILE_OPEN(icomm,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)

    iwrite=0
    do Morton_no=Morton_start(mype),Morton_stop(mype)
       igrid=sfc_to_igrid(Morton_no)
       if (nwaux>0) then
          ! extra layer around mesh only for later averaging in convert
          ! set dxlevel value for use in gradient subroutine, 
          ! which might be used in getaux
          saveigrid=igrid
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          if (.not.slab) mygeo => pgeo(igrid)
          if (B0field) then
             myB0_cell => pB0_cell(igrid)
             {^D&myB0_face^D => pB0_face^D(igrid)\}
          end if
          call phys_get_aux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
       endif
       iwrite=iwrite+1
       {#IFDEF EVOLVINGBOUNDARY
       nphyboundblock=sum(sfc_phybound(1:Morton_no-1))
       offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
            *int(Morton_no-1-nphyboundblock,kind=MPI_OFFSET_KIND) + &
            int(size_block,kind=MPI_OFFSET_KIND) &
            *int(nphyboundblock,kind=MPI_OFFSET_KIND)
       if (sfc_phybound(Morton_no)==1) then
          call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,&
               type_block,istatus,ierrmpi)
       else
          call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,&
               type_block_io,istatus,ierrmpi)
       end if
       }{#IFNDEF EVOLVINGBOUNDARY
       offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
            *int(Morton_no-1,kind=MPI_OFFSET_KIND)
       call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,&
            type_block_io, istatus,ierrmpi)
       }
    end do

    call MPI_FILE_CLOSE(file_handle,ierrmpi)
    if (mype==0) then
       amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
       call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL, &
            file_handle,ierrmpi)

       call write_forest(file_handle)

       {nx^D=ixMhi^D-ixMlo^D+1
       call MPI_FILE_WRITE(file_handle,nx^D,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,nxlone^D,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,xprobmin^D,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,xprobmax^D,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)\}
       call MPI_FILE_WRITE(file_handle,nleafs,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,levmax,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,ndim,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,ndir,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,nw,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,it,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,t,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)
       {#IFDEF EVOLVINGBOUNDARY
       nphyboundblock=sum(sfc_phybound)
       call MPI_FILE_WRITE(file_handle,nphyboundblock,1,MPI_INTEGER,istatus,ierrmpi)
       }
       call MPI_FILE_CLOSE(file_handle,ierrmpi)
    end if
    snapshot=snapshot+1

  end subroutine write_snapshot

  subroutine write_snapshot_tf
    use mod_usr_methods, only: usr_transform_w
    use mod_forest
    use mod_global_parameters
    use mod_physics

    double precision, allocatable :: wtf(:^D&,:)
    double precision :: eqpar_tf(neqpartf)
    integer :: file_handle_tf
    character(len=80) :: filenametf
    integer :: file_handle, amode, igrid, Morton_no, iwrite
    integer :: nx^D
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer, dimension(MPI_STATUS_SIZE) :: istatus
    character(len=80) :: filename, line
    logical, save :: firstsnapshot=.true.
    !-----------------------------------------------------------------------------
    if (firstsnapshot) then
       snapshot=snapshotnext
       firstsnapshot=.false.
    end if

    if (snapshot >= 10000) then
       if (mype==0) then
          write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
          write(*,*) "overwriting first frames"
       end if
       snapshot=0
    end if

    ! generate filename
    write(filenametf,"(a,i4.4,a)") TRIM(filenameout),snapshot,"tf.dat"
    if(mype==0) then
       open(unit=unitsnapshot,file=filenametf,status='replace')
       close(unit=unitsnapshot)
    end if
    amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
    call MPI_FILE_OPEN(icomm,filenametf,amode,MPI_INFO_NULL,file_handle_tf,ierrmpi)
    allocate(wtf(ixG^T,1:nwtf))

    iwrite=0
    do Morton_no=Morton_start(mype),Morton_stop(mype)
       igrid=sfc_to_igrid(Morton_no)
       if (nwaux>0) then
          ! extra layer around mesh only for later averaging in convert
          ! set dxlevel value for use in gradient subroutine, 
          ! which might be used in getaux
          saveigrid=igrid
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          if (.not.slab) mygeo => pgeo(igrid)
          if (B0field) then
             myB0_cell => pB0_cell(igrid)
             {^D&myB0_face^D => pB0_face^D(igrid)\}
          end if
          call phys_get_aux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
       endif
       iwrite=iwrite+1

       if (associated(usr_transform_w)) then
          call usr_transform_w(pw(igrid)%w,wtf,eqpar_tf,ixG^LL,ixM^LL)
       end if

       offset=int(size_block_io_tf,kind=MPI_OFFSET_KIND) &
            *int(Morton_no-1,kind=MPI_OFFSET_KIND)
       call MPI_FILE_WRITE_AT(file_handle_tf,offset,wtf,1, &
            type_block_io_tf,istatus,ierrmpi)     
    end do

    call MPI_FILE_CLOSE(file_handle_tf,ierrmpi)
    amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
    call MPI_FILE_OPEN(MPI_COMM_SELF,filenametf,amode,MPI_INFO_NULL, &
         file_handle_tf,ierrmpi)

    call write_forest(file_handle_tf)

    {nx^D=ixMhi^D-ixMlo^D+1
    call MPI_FILE_WRITE(file_handle_tf,nx^D,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_WRITE(file_handle_tf,nxlone^D,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_WRITE(file_handle_tf,xprobmin^D,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)
    call MPI_FILE_WRITE(file_handle_tf,xprobmax^D,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)\}
    call MPI_FILE_WRITE(file_handle_tf,nleafs,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_WRITE(file_handle_tf,levmax,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_WRITE(file_handle_tf,ndim,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_WRITE(file_handle_tf,ndir,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_WRITE(file_handle_tf,nwtf,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_WRITE(file_handle_tf,it,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_WRITE(file_handle_tf,t,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)

    call MPI_FILE_CLOSE(file_handle_tf,ierrmpi)
    snapshot=snapshot+1

  end subroutine write_snapshot_tf

  subroutine write_snapshot_nopar
    use mod_forest
    use mod_global_parameters
    use mod_physics

    integer :: file_handle, amode, igrid, Morton_no, iwrite
    integer :: nx^D

    integer(kind=MPI_OFFSET_KIND) :: offset

    integer, allocatable :: iorecvstatus(:,:),ioastatus(:,:)
    integer, allocatable :: igrecvstatus(:,:)
    integer, allocatable :: igrid_recv(:) 

    integer, dimension(MPI_STATUS_SIZE) :: istatus

    integer  :: ipe,insend,inrecv,nrecv,nwrite
    character(len=80) :: filename, line
    logical, save :: firstsnapshot=.true.
    !-----------------------------------------------------------------------------
    call MPI_BARRIER(icomm,ierrmpi)

    if (firstsnapshot) then
       snapshot=snapshotnext
       firstsnapshot=.false.
    end if

    if (snapshot >= 10000) then
       if (mype==0) then
          write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
          write(*,*) "overwriting first frames"
       end if
       snapshot=0
    end if

    nrecv=0
    inrecv=0
    nwrite=0
    insend=0
    iwrite=0

    if (mype /= 0) then
       do Morton_no=Morton_start(mype),Morton_stop(mype)
          igrid=sfc_to_igrid(Morton_no)
          itag=Morton_no
          insend=insend+1
          if (nwaux>0) then
             ! extra layer around mesh only for later averaging in convert
             ! set dxlevel value for use in gradient subroutine, 
             ! which might be used in getaux
             saveigrid=igrid
             ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
             mygeo =>pgeo(igrid)
             if (B0field) then
                myB0_cell => pB0_cell(igrid)
                {^D&myB0_face^D => pB0_face^D(igrid)\}
             end if
             call phys_get_aux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
          endif
          call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
          call MPI_SEND(pw(igrid)%w,1,type_block_io, 0,itag,icomm,ierrmpi)
       end do
    else 
       ! mype==0
       nwrite=(Morton_stop(0)-Morton_start(0)+1)

       ! master processor writes out
       write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"

       open(unit=unitsnapshot,file=filename,status='replace')
       close(unitsnapshot, status='delete')

       amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
       call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)

       ! writing his local data first
       do Morton_no=Morton_start(0),Morton_stop(0)
          igrid=sfc_to_igrid(Morton_no)
          iwrite=iwrite+1
          if (nwaux>0) then
             ! extra layer around mesh only for later averaging in convert
             ! set dxlevel value for use in gradient subroutine,
             ! which might be used in getaux
             saveigrid=igrid
             ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
             mygeo =>pgeo(igrid)
             if (B0field) then
                myB0_cell => pB0_cell(igrid)
                {^D&myB0_face^D => pB0_face^D(igrid)\}
             end if
             call phys_get_aux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
          endif
          {#IFDEF EVOLVINGBOUNDARY
          nphyboundblock=sum(sfc_phybound(1:Morton_no-1))
          offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
               *int(Morton_no-1-nphyboundblock,kind=MPI_OFFSET_KIND) + &
               int(size_block,kind=MPI_OFFSET_KIND) &
               *int(nphyboundblock,kind=MPI_OFFSET_KIND)
          if (sfc_phybound(Morton_no)==1) then
             call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,type_block, &
                  istatus,ierrmpi)
          else
             call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,&
                  type_block_io,istatus,ierrmpi)
          end if
          }{#IFNDEF EVOLVINGBOUNDARY
          offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
               *int(Morton_no-1,kind=MPI_OFFSET_KIND)
          call MPI_FILE_WRITE_AT(file_handle,offset,pw(igrid)%w,1,type_block_io, &
               istatus,ierrmpi)
          }
       end do
       ! write data communicated from other processors
       if(npe>1)then
          nrecv=(Morton_stop(npe-1)-Morton_start(1)+1)
          inrecv=0
          allocate(igrid_recv(nrecv))
          allocate(igrecvstatus(MPI_STATUS_SIZE,nrecv),iorecvstatus(MPI_STATUS_SIZE,nrecv))
          allocate(ioastatus(MPI_STATUS_SIZE,nrecv))

          do ipe =1, npe-1
             do Morton_no=Morton_start(ipe),Morton_stop(ipe)
                iwrite=iwrite+1
                itag=Morton_no
                inrecv=inrecv+1
                call MPI_RECV(igrid_recv(inrecv),1,MPI_INTEGER, ipe,itag,icomm,&
                     igrecvstatus(:,inrecv),ierrmpi)

                allocate(pwio(igrid_recv(inrecv))%w(ixG^T,1:nw))
                call MPI_RECV(pwio(igrid_recv(inrecv))%w,1,type_block_io,ipe,itag,icomm,&
                     iorecvstatus(:,inrecv),ierrmpi)
                {#IFDEF EVOLVINGBOUNDARY
                nphyboundblock=sum(sfc_phybound(1:Morton_no-1))
                offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
                     *int(Morton_no-1-nphyboundblock,kind=MPI_OFFSET_KIND) + &
                     int(size_block,kind=MPI_OFFSET_KIND) &
                     *int(nphyboundblock,kind=MPI_OFFSET_KIND)
                if (sfc_phybound(Morton_no)==1) then
                   call MPI_FILE_WRITE_AT(file_handle,offset,pwio(igrid_recv(inrecv))%w,1,&
                        type_block   ,ioastatus(:,inrecv),ierrmpi)
                else
                   call MPI_FILE_WRITE_AT(file_handle,offset,pwio(igrid_recv(inrecv))%w,1,&
                        type_block_io,ioastatus(:,inrecv),ierrmpi)
                end if
                }{#IFNDEF EVOLVINGBOUNDARY
                offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
                     *int(Morton_no-1,kind=MPI_OFFSET_KIND)
                call MPI_FILE_WRITE_AT(file_handle,offset,pwio(igrid_recv(inrecv))%w,1,&
                     type_block_io,ioastatus(:,inrecv),ierrmpi)
                }
                deallocate(pwio(igrid_recv(inrecv))%w)
             end do
          end do
          deallocate(igrecvstatus,iorecvstatus,ioastatus,igrid_recv)
       end if
    end if

    if(mype==0) call MPI_FILE_CLOSE(file_handle,ierrmpi)

    if (mype==0) then
       amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
       call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL, &
            file_handle,ierrmpi)

       call write_forest(file_handle)

       {nx^D=ixMhi^D-ixMlo^D+1
       call MPI_FILE_WRITE(file_handle,nx^D,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,nxlone^D,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,xprobmin^D,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,xprobmax^D,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)\}
       call MPI_FILE_WRITE(file_handle,nleafs,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,levmax,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,ndim,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,ndir,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,nw,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,it,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_WRITE(file_handle,t,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)
       {#IFDEF EVOLVINGBOUNDARY
       nphyboundblock=sum(sfc_phybound)
       call MPI_FILE_WRITE(file_handle,nphyboundblock,1,MPI_INTEGER,istatus,ierrmpi)
       }
       call MPI_FILE_CLOSE(file_handle,ierrmpi)
    end if
    snapshot=snapshot+1

    call MPI_BARRIER(icomm,ierrmpi)

  end subroutine write_snapshot_nopar

  subroutine write_snapshot_noparf
    use mod_forest
    use mod_global_parameters
    use mod_physics

    integer :: igrid, Morton_no
    integer :: nx^D

    integer, allocatable :: iorecvstatus(:,:),ioastatus(:,:)
    integer, allocatable :: igrecvstatus(:,:)
    integer, allocatable :: igrid_recv(:) 

    integer  :: ipe,insend,inrecv,nrecv,nwrite
    character(len=80) :: filename, line
    logical, save :: firstsnapshot=.true.
    !-----------------------------------------------------------------------------
    call MPI_BARRIER(icomm,ierrmpi)

    if (firstsnapshot) then
       snapshot=snapshotnext
       firstsnapshot=.false.
    end if

    if (snapshot >= 10000) then
       if (mype==0) then
          write(*,*) "WARNING: Number of frames is limited to 10000 (0...9999),"
          write(*,*) "overwriting first frames"
       end if
       snapshot=0
    end if

    nrecv=0
    inrecv=0
    nwrite=0
    insend=0

    if (mype /= 0) then
       do Morton_no=Morton_start(mype),Morton_stop(mype)
          igrid=sfc_to_igrid(Morton_no)
          itag=Morton_no
          insend=insend+1
          if (nwaux>0) then
             ! extra layer around mesh only for later averaging in convert
             ! set dxlevel value for use in gradient subroutine, 
             ! which might be used in getaux
             saveigrid=igrid
             ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
             mygeo =>pgeo(igrid)
             if (B0field) then
                myB0_cell => pB0_cell(igrid)
                {^D&myB0_face^D => pB0_face^D(igrid)\}
             end if
             call phys_get_aux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
          endif
          call MPI_SEND(igrid,1,MPI_INTEGER, 0,itag,icomm,ierrmpi)
          call MPI_SEND(pw(igrid)%w,1,type_block_io, 0,itag,icomm,ierrmpi)
       end do
    else 
       ! mype==0
       nwrite=(Morton_stop(0)-Morton_start(0)+1)

       ! master processor writes out
       write(filename,"(a,i4.4,a)") TRIM(filenameout),snapshot,".dat"
       if(endian_swap) then
          {#IFNDEF BIGENDIAN
          open(unit=unitsnapshot,file=filename,form='unformatted',access='stream',&
               status='replace',convert='BIG_ENDIAN')
          }
          {#IFDEF BIGENDIAN
          open(unit=unitsnapshot,file=filename,form='unformatted',access='stream',&
               status='replace',convert='LITTLE_ENDIAN')
          }
       else
          open(unit=unitsnapshot,file=filename,form='unformatted',access='stream',&
               status='replace')
       end if
       ! writing his local data first
       do Morton_no=Morton_start(0),Morton_stop(0)
          igrid=sfc_to_igrid(Morton_no)
          if (nwaux>0) then
             ! extra layer around mesh only for later averaging in convert
             ! set dxlevel value for use in gradient subroutine,
             ! which might be used in getaux
             saveigrid=igrid
             ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
             mygeo =>pgeo(igrid)
             if (B0field) then
                myB0_cell => pB0_cell(igrid)
                {^D&myB0_face^D => pB0_face^D(igrid)\}
             end if
             call phys_get_aux(.true.,pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL^LADD1,"write_snapshot")
          endif
          write(unitsnapshot) pw(igrid)%w(ixM^T,1:nw)
       end do
       ! write data communicated from other processors
       if(npe>1)then
          nrecv=(Morton_stop(npe-1)-Morton_start(1)+1)
          inrecv=0
          allocate(igrid_recv(nrecv))
          allocate(igrecvstatus(MPI_STATUS_SIZE,nrecv),iorecvstatus(MPI_STATUS_SIZE,nrecv))
          allocate(ioastatus(MPI_STATUS_SIZE,nrecv))

          do ipe =1, npe-1
             do Morton_no=Morton_start(ipe),Morton_stop(ipe)
                itag=Morton_no
                inrecv=inrecv+1
                call MPI_RECV(igrid_recv(inrecv),1,MPI_INTEGER, ipe,itag,icomm,&
                     igrecvstatus(:,inrecv),ierrmpi)

                allocate(pwio(igrid_recv(inrecv))%w(ixG^T,1:nw))
                call MPI_RECV(pwio(igrid_recv(inrecv))%w,1,type_block_io,ipe,itag,icomm,&
                     iorecvstatus(:,inrecv),ierrmpi)

                write(unitsnapshot) pwio(igrid_recv(inrecv))%w(ixM^T,1:nw)
                deallocate(pwio(igrid_recv(inrecv))%w)
             end do
          end do
          deallocate(igrecvstatus,iorecvstatus,ioastatus,igrid_recv)
       end if
    end if

    if(mype==0) then
       call write_forest(unitsnapshot)
       {nx^D=ixMhi^D-ixMlo^D+1
       write(unitsnapshot) nx^D
       write(unitsnapshot) nxlone^D
       write(unitsnapshot) xprobmin^D
       write(unitsnapshot) xprobmax^D\}
       write(unitsnapshot) nleafs
       write(unitsnapshot) levmax 
       write(unitsnapshot) ndim 
       write(unitsnapshot) ndir
       write(unitsnapshot) nw
       write(unitsnapshot) it
       write(unitsnapshot) t
       close(unitsnapshot)
    end if

    snapshot=snapshot+1

    call MPI_BARRIER(icomm,ierrmpi)

  end subroutine write_snapshot_noparf

  subroutine read_snapshot
    use mod_forest
    use mod_global_parameters

    integer :: file_handle, amode, igrid, Morton_no, iread
    integer :: levmaxini, ndimini, ndirini, nwini, nxini^D, nxloneini^D
    double precision :: xprobminini^D,xprobmaxini^D

    {^IFMPT integer :: size_double, size_int,lb}
    {^IFNOMPT  integer(kind=MPI_ADDRESS_KIND) :: size_double, size_int, lb}

    integer(kind=MPI_OFFSET_KIND) :: offset
    integer, dimension(MPI_STATUS_SIZE) :: istatus
    character(len=80) :: filename
    logical :: fexist
    !-----------------------------------------------------------------------------
    ! generate filename
    write(filename,"(a,i4.4,a)") TRIM(filenameini),snapshotini,".dat"

    if(mype==0) then
       inquire(file=filename,exist=fexist)
       if(.not.fexist) call mpistop(filename//"as an input snapshot file is not found!")
    endif

    amode=MPI_MODE_RDONLY
    call MPI_FILE_OPEN(icomm,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)


    call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
    call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,size_int,ierrmpi)

    {#IFDEF EVOLVINGBOUNDARY
    offset=-int(7*size_int+size_double,kind=MPI_OFFSET_KIND)
    }{#IFNDEF EVOLVINGBOUNDARY
    offset=-int(6*size_int+size_double,kind=MPI_OFFSET_KIND)
    }
    call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,nleafs,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,levmaxini,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,ndimini,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,ndirini,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,nwini,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,it,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,t,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)
    {#IFDEF EVOLVINGBOUNDARY
    call MPI_FILE_READ_ALL(file_handle,nphyboundblock,1,MPI_INTEGER,istatus,ierrmpi)
    }
    nleafs_active = nleafs

    ! check if settings are suitable for restart
    if (levmaxini>mxnest) then
       if (mype==0) write(*,*) "number of levels in restart file = ",levmaxini
       if (mype==0) write(*,*) "mxnest = ",mxnest
       call mpistop("mxnest should be at least number of levels in restart file")
    end if
    if (ndimini/=ndim) then
       if (mype==0) write(*,*) "ndim in restart file = ",ndimini
       if (mype==0) write(*,*) "ndim = ",ndim
       call mpistop("reset ndim to ndim in restart file")
    end if
    if (ndirini/=ndir) then
       if (mype==0) write(*,*) "ndir in restart file = ",ndirini
       if (mype==0) write(*,*) "ndir = ",ndir
       call mpistop("reset ndir to ndir in restart file")
    end if
    if (nw/=nwini) then
       if (mype==0) write(*,*) "nw=",nw," and nw in restart file=",nwini
       call mpistop("currently, changing nw at restart is not allowed")
    end if

    offset=offset-int(2*ndimini*size_int+2*ndimini*size_double,kind=MPI_OFFSET_KIND)
    !call MPI_FILE_SEEK_SHARED(file_handle,offset,MPI_SEEK_END,ierrmpi)
    call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

   {call MPI_FILE_READ_ALL(file_handle,nxini^D,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,nxloneini^D,1,MPI_INTEGER,istatus,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,xprobminini^D,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,xprobmaxini^D,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)\}
    if (nxblock^D/=nxini^D|.or.) then
       if (mype==0) write(*,*) "Error: reset block resolution to nxblock^D=",nxini^D
       call mpistop("change nxblock^D in par file")
    end if
    if (nxlone^D/=nxloneini^D|.or.) then
       if (mype==0) write(*,*) "Error: resolution of base mesh does not match the data: ",nxloneini^D
       call mpistop("change nxlone^D in par file")
    end if
    if (xprobmin^D/=xprobminini^D|.or.) then
       if (mype==0) write(*,*) "Error: location of minimum does not match the data: ",xprobminini^D
       call mpistop("change xprobmin^D in par file")
    end if
    if (xprobmax^D/=xprobmaxini^D|.or.) then
       if (mype==0) write(*,*) "Error: location of maximum does not match the data: ",xprobmaxini^D
       call mpistop("change xprobmax^D in par file")
    end if

    call read_forest(file_handle)

    iread=0
    {#IFDEF EVOLVINGBOUNDARY
    ! mark physical-boundary blocks on space-filling curve
    do Morton_no=Morton_start(mype),Morton_stop(mype)
       igrid=sfc_to_igrid(Morton_no)
       call alloc_node(igrid)
       if (phyboundblock(igrid)) sfc_phybound(Morton_no)=1
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,sfc_phybound,nleafs,MPI_INTEGER,&
         MPI_SUM,icomm,ierrmpi)

    do Morton_no=Morton_start(mype),Morton_stop(mype)
       igrid=sfc_to_igrid(Morton_no)
       iread=iread+1
       nphyboundblock=sum(sfc_phybound(1:Morton_no-1))
       offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
            *int(Morton_no-1-nphyboundblock,kind=MPI_OFFSET_KIND) + &
            int(size_block,kind=MPI_OFFSET_KIND) &
            *int(nphyboundblock,kind=MPI_OFFSET_KIND)
       if (sfc_phybound(Morton_no)==1) then
          call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1, &
               type_block,istatus,ierrmpi)
       else
          call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1, &
               type_block_io,istatus,ierrmpi)
       end if
    end do
    }{#IFNDEF EVOLVINGBOUNDARY
    do Morton_no=Morton_start(mype),Morton_stop(mype)
       igrid=sfc_to_igrid(Morton_no)
       call alloc_node(igrid)
       iread=iread+1
       offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
            *int(Morton_no-1,kind=MPI_OFFSET_KIND)
       call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1, &
            type_block_io,istatus,ierrmpi)
    end do
    }

    call MPI_FILE_CLOSE(file_handle,ierrmpi)

!!!call MPI_BARRIER(icomm,ierrmpi)
  end subroutine read_snapshot

  subroutine read_snapshotnopar
    use mod_forest
    use mod_global_parameters

    double precision :: wio(ixG^T,1:nw)
    integer :: file_handle, amode, igrid, Morton_no, iread
    integer :: levmaxini, ndimini, ndirini, nwini, neqparini, nxini^D

    {^IFMPT integer :: size_double, size_int,lb}
    {^IFNOMPT  integer(kind=MPI_ADDRESS_KIND) :: size_double, size_int, lb}

    integer(kind=MPI_OFFSET_KIND) :: offset
    integer, dimension(MPI_STATUS_SIZE) :: istatus

    integer, allocatable :: iorecvstatus(:,:)
    integer :: ipe,inrecv,nrecv
    character(len=80) :: filename
    logical :: fexist
    !-----------------------------------------------------------------------------
    ! generate filename
    write(filename,"(a,i4.4,a)") TRIM(filenameini),snapshotini,".dat"
    if (mype==0) then
       inquire(file=filename,exist=fexist)
       if(.not.fexist) call mpistop(filename//"as an input snapshot file is not found!")
       amode=MPI_MODE_RDONLY
       call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,file_handle,ierrmpi)

       call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
       call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,size_int,ierrmpi)

       {#IFDEF EVOLVINGBOUNDARY
       offset=-int(7*size_int+size_double,kind=MPI_OFFSET_KIND)
       }{#IFNDEF EVOLVINGBOUNDARY
       offset=-int(6*size_int+size_double,kind=MPI_OFFSET_KIND)
       }
       call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)
       call MPI_FILE_READ(file_handle,nleafs,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_READ(file_handle,levmaxini,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_READ(file_handle,ndimini,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_READ(file_handle,ndirini,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_READ(file_handle,nwini,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_READ(file_handle,it,1,MPI_INTEGER,istatus,ierrmpi)
       call MPI_FILE_READ(file_handle,t,1,MPI_DOUBLE_PRECISION,istatus,ierrmpi)
       {#IFDEF EVOLVINGBOUNDARY
       call MPI_FILE_READ(file_handle,nphyboundblock,1,MPI_INTEGER,istatus,ierrmpi)
       }
       ! check if settings are suitable for restart
       if (levmaxini>mxnest) then
          write(*,*) "number of levels in restart file = ",levmaxini
          write(*,*) "mxnest = ",mxnest
          call mpistop("mxnest should be at least number of levels in restart file")
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


       offset=offset-int(2*ndimini*size_int+2*ndimini*size_double,kind=MPI_OFFSET_KIND)
       call MPI_FILE_SEEK(file_handle,offset,MPI_SEEK_END,ierrmpi)

       {call MPI_FILE_READ(file_handle,nxini^D,1,MPI_INTEGER,istatus,ierrmpi)\}
       if (nxblock^D/=nxini^D|.or.) then
          write(*,*) "Error: reset block resolution to nxblock^D=",nxini^D
          call mpistop("change nxblock^D in par file")
       end if
    end if

    ! broadcast the global parameters first
    if (npe>1) then
       call MPI_BCAST(nleafs,1,MPI_INTEGER,0,icomm,ierrmpi)
       call MPI_BCAST(it,1,MPI_INTEGER,0,icomm,ierrmpi)
       call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    end if
    nleafs_active = nleafs

    call read_forest(file_handle)
    {#IFDEF EVOLVINGBOUNDARY
    ! mark physical-boundary blocks on space-filling curve
    do Morton_no=Morton_start(mype),Morton_stop(mype)
       igrid=sfc_to_igrid(Morton_no)
       call alloc_node(igrid)
       if (phyboundblock(igrid)) sfc_phybound(Morton_no)=1
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,sfc_phybound,nleafs,MPI_INTEGER,&
         MPI_SUM,icomm,ierrmpi)
    }{#IFNDEF EVOLVINGBOUNDARY
    do Morton_no=Morton_start(mype),Morton_stop(mype)
       igrid=sfc_to_igrid(Morton_no)
       call alloc_node(igrid)
    end do
    }

    if (mype==0)then
       iread=0

       do Morton_no=Morton_start(0),Morton_stop(0)
          igrid=sfc_to_igrid(Morton_no)
          iread=iread+1
    {#IFDEF EVOLVINGBOUNDARY
          nphyboundblock=sum(sfc_phybound(1:Morton_no-1))
          offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
               *int(Morton_no-1-nphyboundblock,kind=MPI_OFFSET_KIND) + &
               int(size_block,kind=MPI_OFFSET_KIND) &
               *int(nphyboundblock,kind=MPI_OFFSET_KIND)
          if (sfc_phybound(Morton_no)==1) then
            call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1,&
                 type_block,istatus,ierrmpi)
          else
            call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1,&
                 type_block_io,istatus,ierrmpi)
          endif
    }{#IFNDEF EVOLVINGBOUNDARY
          offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
               *int(Morton_no-1,kind=MPI_OFFSET_KIND)
          call MPI_FILE_READ_AT(file_handle,offset,pw(igrid)%w,1,type_block_io, &
               istatus,ierrmpi)
    }
       end do
       if (npe>1) then
          do ipe=1,npe-1
             do Morton_no=Morton_start(ipe),Morton_stop(ipe)
                iread=iread+1
                itag=Morton_no
    {#IFDEF EVOLVINGBOUNDARY
                nphyboundblock=sum(sfc_phybound(1:Morton_no-1))
                offset=int(size_block_io,kind=MPI_OFFSET_KIND) &
                     *int(Morton_no-1-nphyboundblock,kind=MPI_OFFSET_KIND) + &
                     int(size_block,kind=MPI_OFFSET_KIND) &
                     *int(nphyboundblock,kind=MPI_OFFSET_KIND)
                if (sfc_phybound(Morton_no)==1) then
                  call MPI_FILE_READ_AT(file_handle,offset,wio,1,type_block,&
                       istatus,ierrmpi)
                  call MPI_SEND(wio,1,type_block, ipe,itag,icomm,ierrmpi)
                else
                  call MPI_FILE_READ_AT(file_handle,offset,wio,1,type_block_io,&
                       istatus,ierrmpi)
                  call MPI_SEND(wio,1,type_block_io, ipe,itag,icomm,ierrmpi)
                endif
    }{#IFNDEF EVOLVINGBOUNDARY
                offset=int(size_block_io,kind=MPI_OFFSET_KIND)&
                     *int(Morton_no-1,kind=MPI_OFFSET_KIND)
                call MPI_FILE_READ_AT(file_handle,offset,wio,1,type_block_io,&
                     istatus,ierrmpi)
                call MPI_SEND(wio,1,type_block_io,ipe,itag,icomm,ierrmpi)
    }
             end do
          end do
       end if
       call MPI_FILE_CLOSE(file_handle,ierrmpi)
    else
       nrecv=(Morton_stop(mype)-Morton_start(mype)+1)
       allocate(iorecvstatus(MPI_STATUS_SIZE,nrecv))
       inrecv=0
       do Morton_no=Morton_start(mype),Morton_stop(mype)
          igrid=sfc_to_igrid(Morton_no)
          itag=Morton_no
          inrecv=inrecv+1
    {#IFDEF EVOLVINGBOUNDARY
          if (sfc_phybound(Morton_no)==1) then
            call MPI_RECV(pw(igrid)%w,1,type_block   ,0,itag,icomm,&
                 iorecvstatus(:,inrecv),ierrmpi)
          else
            call MPI_RECV(pw(igrid)%w,1,type_block_io,0,itag,icomm,&
                 iorecvstatus(:,inrecv),ierrmpi)
          endif
    }{#IFNDEF EVOLVINGBOUNDARY
          call MPI_RECV(pw(igrid)%w,1,type_block_io,0,itag,icomm,&
               iorecvstatus(:,inrecv),ierrmpi)
    }
       end do
       deallocate(iorecvstatus)
    end if

    call MPI_BARRIER(icomm,ierrmpi)

  end subroutine read_snapshotnopar

  !> Write volume-averaged values and other information to the log file
  subroutine printlog_default

    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters

    logical              :: fileopen
    integer              :: i, iw, level
    double precision     :: wmean(1:nw), total_volume
    double precision     :: volume_coverage(mxnest)
    integer              :: nx^D, nc, ncells, dit
    double precision     :: dtTimeLast, now, cellupdatesPerSecond
    double precision     :: activeBlocksPerCore, wctPerCodeTime, timeToFinish
    character(len=40)    :: fmt_string
    character(len=80)    :: filename
    character(len=2048)  :: line
    logical, save        :: opened  = .false.
    integer              :: amode, istatus(MPI_STATUS_SIZE)

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
       cellupdatesPerSecond = dble(ncells) * dble(nstep) * dble(dit) / (dtTimeLast * dble(npe))

       ! blocks per core:
       activeBlocksPerCore = dble(nleafs_active) / dble(npe)

       ! Wall clock time per code time unit in seconds:
       wctPerCodeTime = dtTimeLast / max(dit * dt, epsilon(1.0d0))

       ! Wall clock time to finish in hours:
       timeToFinish = (tmax - t) * wctPerCodeTime / 3600.0d0

       ! On first entry, open the file and generate the header
       if (.not. opened) then

          filename = trim(filenamelog) // ".log"
          amode    = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
          amode    = ior(amode,MPI_MODE_APPEND)

          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, &
               MPI_INFO_NULL, log_fh, ierrmpi)

          opened = .true.

          call MPI_FILE_WRITE(log_fh, trim(fileheadout) // new_line('a'), &
               len_trim(fileheadout)+1, MPI_CHARACTER, istatus, ierrmpi)

          ! Start of file headern
          line = "it t dt res " // trim(wnames)

          ! Volume coverage per level
          do level = 1, mxnest
             i = len_trim(line) + 2
             write(line(i:), "(a,i0)") "c", level
          end do

          ! Cell counts per level
          do level=1,mxnest
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

       fmt_string = '(' // fmt_i // ',3' // fmt_r // ')'
       write(line, fmt_string) it, t, dt, residual
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', nw, fmt_r // ')'
       write(line(i:), fmt_string) wmean(1:nw)
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', mxnest, fmt_r // ')'
       write(line(i:), fmt_string) volume_coverage(1:mxnest)
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', mxnest, fmt_i // ')'
       write(line(i:), fmt_string) nleafs_level(1:mxnest)
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
          open(my_unit, file=trim(filenamelog)//".log")
          file_open = .true.

          write(my_unit, *) "# time mean(w) mean(w**2)"
       end if

       write(fmt_string, "(a,i0,a)") "(", nw * n_modes + 1, fmt_r // ")"
       write(my_unit, fmt_string) t, modes
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
          dvolume(ixM^T) = pgeo(igrid)%dvolume(ixM^T)
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

    double precision, intent(out) :: vol_cov(1:mxnest)
    double precision              :: dsum_recv(1:mxnest)
    integer                       :: iigrid, igrid, iw, level

    ! First determine the total 'flat' volume in each level
    vol_cov(1:mxnest)=zero

    do iigrid = 1, igridstail
       igrid          = igrids(iigrid);
       level          = node(plevel_,igrid)
       vol_cov(level) = vol_cov(level)+ &
            {(rnode(rpxmax^D_,igrid)-rnode(rpxmin^D_,igrid))|*}
    end do

    ! Make the information available on all tasks
    call MPI_ALLREDUCE(vol_cov, dsum_recv, mxnest, MPI_DOUBLE_PRECISION, &
         MPI_SUM, icomm, ierrmpi)

    ! Normalize
    vol_cov = dsum_recv / sum(dsum_recv)
  end subroutine get_volume_coverage

end module mod_input_output
