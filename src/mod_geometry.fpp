#:mute
#:include "mod_gpu_directives.fpp"
#:endmute
!> Module with geometry-related routines (e.g., divergence, curl)
module mod_geometry
  use mod_comm_lib, only: mpistop
  implicit none
  public

  integer :: coordinate=-1
  ${GPU_DECLARE_COPYIN('coordinate')}$
  integer, parameter :: Cartesian          = 0
  integer, parameter :: Cartesian_stretched= 1
  integer, parameter :: cylindrical        = 2
  integer, parameter :: spherical          = 3
  integer, parameter :: Cartesian_expansion= 4

  integer :: type_curl=0
  ${GPU_DECLARE_COPYIN('type_curl')}$
  integer, parameter :: central=1
  integer, parameter :: Gaussbased=2
  integer, parameter :: Stokesbased=3

contains

  !> Set the coordinate system from the GEOM compile-time define.
  !>
  !> Called once before usr_init, so a case does not have to call
  !> set_coordinate_system itself. Cases that do are still fine: the call is
  !> idempotent, and the consistency check in set_coordinate_system catches a
  !> name that disagrees with what the code was built for.
  !>
  !> Doing this centrally matters because set_coordinate_system is the only
  !> writer of coordinate, ndir, r_, phi_ and z_. Leaving coordinate at its
  !> -1 default is silently catastrophic rather than loud: read_par_files
  !> skips the conversion of angular domain bounds, and the select case in
  !> get_surface_area and in the cell-volume fill match nothing, so surfaceC
  !> and dvolume are simply never written.
  subroutine set_coordinate_system_from_config()
#:if GEOM == 'spherical'
  #:if defined('LOG_RADIUS')
    call set_coordinate_system("logSpherical_3D")
  #:else
    call set_coordinate_system("spherical_3D")
  #:endif
#:elif GEOM == 'cylindrical'
  #:if defined('LOG_RADIUS')
    call set_coordinate_system("logCylindrical_3D")
  #:else
    call set_coordinate_system("cylindrical_3D")
  #:endif
#:else
    call set_coordinate_system("Cartesian_3D")
#:endif
  end subroutine set_coordinate_system_from_config

  !> Set the coordinate system to be used
  subroutine set_coordinate_system(geom)
    use mod_global_parameters
    character(len=*), intent(in) :: geom !< Name of the coordinate system

    ! Store the geometry name
    geometry_name = geom

    select case (geom)
    case ("Cartesian","Cartesian_1D","Cartesian_2D","Cartesian_3D")
      ndir = ndim
      coordinate=Cartesian
    case ("Cartesian_1D_expansion")
      if (ndim /= 1) call mpistop&
         ("Geometry Cartesian_1D_expansion but ndim /= 1")
      ndir = ndim
      coordinate=Cartesian_expansion
    case ("Cartesian_1.5D")
      if (ndim /= 1) call mpistop("Geometry Cartesian_1.5D but ndim /= 1")
      coordinate=Cartesian
      ndir = 2
    case ("Cartesian_1.75D")
      if (ndim /= 1) call mpistop("Geometry Cartesian_1.75D but ndim /= 1")
      coordinate=Cartesian
      ndir = 3
    case ("Cartesian_2.5D")
      if (ndim /= 2) call mpistop("Geometry Cartesian_2.5D but ndim /= 2")
      coordinate=Cartesian
      ndir = 3
    case ("cylindrical","cylindrical_2D","cylindrical_3D")
      ndir = ndim
      r_   = 1
      z_   = 2
      if(ndir==3) phi_ = 3
      coordinate=cylindrical
    case ("cylindrical_2.5D")
      if (ndim /= 2) call mpistop("Geometry cylindrical_2.5D but ndim /= 2")
      ndir = 3
      r_   = 1
      z_   = 2
      phi_ = 3
      coordinate=cylindrical
    case ("polar","polar_2D","polar_3D")
      ndir = ndim
      r_   = 1
      phi_ = 2
      if(ndir==3) z_ = 3
      coordinate=cylindrical
    case ("polar_1.5D")
       if (ndim /= 1) call mpistop("Geometry polar_1.5D but ndim /= 1")
       ndir = 2
       r_   = 1
       phi_ = 2
       coordinate=cylindrical
    case ("polar_2.5D")
      if (ndim /= 2) call mpistop("Geometry polar_2.5D but ndim /= 2")
      ndir = 3
      r_   = 1
      phi_ = 2
      z_   = 3
      coordinate=cylindrical
    case ("spherical","spherical_2D","spherical_3D")
      ndir = ndim
      r_   = 1
      if(ndir==3) phi_ = 3
      z_   = -1
      coordinate=spherical
    case ("spherical_2.5D")
      if (ndim /= 2) call mpistop("Geometry spherical_2.5D requires ndim == 2")

      ndir = 3
      r_   = 1
      phi_ = 3
      z_   = -1
      coordinate=spherical
    ! The two logarithmically stretched systems are the same coordinate
    ! systems - identical ndir, r_, phi_, z_ and the same `coordinate` value,
    ! so every select case (coordinate) in the code keeps matching. What
    ! differs is only that the mesh is uniform in ln(r) rather than in r, which
    ! the LOG_RADIUS define carries. The names are kept distinct here because
    ! geometry_name goes into the snapshot header, and a log snapshot must not
    ! restart into a uniformly spaced build.
    case ("logSpherical","logSpherical_3D")
      ndir = ndim
      r_   = 1
      if(ndir==3) phi_ = 3
      z_   = -1
      coordinate=spherical
    case ("logCylindrical","logCylindrical_3D")
      ndir = ndim
      r_   = 1
      z_   = 2
      if(ndir==3) phi_ = 3
      coordinate=cylindrical
    case default
      call mpistop("Unknown geometry specified")
    end select

    ! The coordinate system is a compile-time choice in AGILE: the GEOM define
    ! (set from `geometry` in &meshlist) specializes the finite-volume kernels.
    ! Guard against a mod_usr that disagrees with what the code was built for.
#:if GEOM == 'spherical'
    if (coordinate /= spherical) call mpistop("meshlist geometry is spherical, &
       &but mod_usr set a non-spherical coordinate system")
  #:if PHYS not in ('hd', 'mhd', 'srhd', 'srmhd', 'ffhd')
    ! Only physics modules that define the addsource_geometry macro can be run
    ! curvilinear; without it the momentum equations would silently miss their
    ! curvature terms.
    call mpistop("geometry='spherical' is only implemented for phys='hd', &
       &phys='mhd', phys='srhd', phys='srmhd' and phys='ffhd'")
  #:endif
#:elif GEOM == 'cylindrical'
    if (coordinate /= cylindrical) call mpistop("meshlist geometry is &
       &cylindrical, but mod_usr set a non-cylindrical coordinate system")
  #:if PHYS not in ('hd', 'mhd', 'srhd', 'srmhd', 'ffhd')
    ! Only physics modules that define the addsource_geometry macro can be run
    ! curvilinear; without it the momentum equations would silently miss their
    ! curvature terms.
    call mpistop("geometry='cylindrical' is only implemented for phys='hd', &
       &phys='mhd', phys='srhd', phys='srmhd' and phys='ffhd'")
  #:endif
#:else
    if (coordinate /= Cartesian) call mpistop("meshlist geometry is Cartesian, &
       &but mod_usr set a non-Cartesian coordinate system")
#:endif

    ! coordinate alone cannot tell a log-stretched system from a uniform one -
    ! that is the point of the design - so the radial spacing is checked
    ! separately, on the name.
#:if defined('LOG_RADIUS')
    if (geom(1:3) /= "log") call mpistop("meshlist geometry is logSpherical or &
       &logCylindrical, but mod_usr set a uniformly spaced coordinate system")
#:else
    if (geom(1:3) == "log") call mpistop("mod_usr set a logarithmically &
       &stretched coordinate system, but meshlist geometry is not &
       &logSpherical or logCylindrical")
#:endif
  end subroutine set_coordinate_system

  !> The radial coordinate map, from the logical coordinate s to the physical
  !> radius r, and its inverse.
  !>
  !> These are the host-side definition of the map a LOG_RADIUS build uses,
  !> and the one place it is written down:
  !>
  !>    log_r0 == 0 :  r = exp(s)                        (s = ln(r))
  !>    log_r0 >  0 :  r = sign(s)*r0*(exp(|s|) - 1)     (s = ln(1 + r/r0))
  !>
  !> The offset form is the one that reaches the axis. Its pivot is deliberate:
  !> s = 0 corresponds to r = 0 exactly, so the map is *odd*, and the mesh in
  !> the ghost layer beyond a cylindrical axis is the exact mirror of the mesh
  !> inside it. That is what the pole copy in getbc assumes when it fills ghost
  !> cell j from interior cell j of the block half a revolution away - the two
  !> have to occupy mirrored volumes, or the copied cell average belongs to a
  !> cell of a different size. The naive map r = exp(s) - r0 is *not* odd:
  !> its ghost cells shrink outward while the cells they mirror grow, and the
  !> first ghost cell is then misplaced by several percent of a cell width.
  !>
  !> The odd branch is only ever taken in that ghost layer; s >= 0 everywhere
  !> in the physical domain of an offset build. Device code that only ever
  !> sees the domain therefore uses the cheaper equivalent
  !> r = log_ra*exp(s) + log_rb instead of calling these.
  !>
  !> A case's mod_usr should use these rather than writing exp/log by hand:
  !> xprobmin1 and xprobmax1 hold the *logical* radial bounds in a LOG_RADIUS
  !> build, so r_of_s is how a case recovers a physical radius from them.
  pure elemental function r_of_s(s) result(r)
    use mod_global_parameters
    double precision, intent(in) :: s
    double precision             :: r

    if (log_r0 > zero) then
      r = dsign(one,s)*log_r0*(dexp(dabs(s)) - one)
    else
      r = dexp(s)
    end if
  end function r_of_s

  !> Inverse of r_of_s. See there.
  pure elemental function s_of_r(r) result(s)
    use mod_global_parameters
    double precision, intent(in) :: r
    double precision             :: s

    if (log_r0 > zero) then
      s = dsign(one,r)*dlog(one + dabs(r)/log_r0)
    else
      s = dlog(r)
    end if
  end function s_of_r

  subroutine set_pole
    use mod_global_parameters

    select case (coordinate)
    case (spherical) 
      ! For spherical grid, check whether phi-direction is periodic
      if(periodB(ndim)) then
        if(phi_/=3) call mpistop("phi_ should be 3 in 3D spherical coord!")
        if(mod(ng3(1),2)/=0) call mpistop(&
           "Number of meshes in phi-direction should be even!")
        if(abs(xprobmin2)<smalldouble) then
          if(mype==0) write(unitterm,*)&
              "Will apply pi-periodic conditions at northpole!"
          poleB(1,2)=.true.
        else
          if(mype==0) write(unitterm,*) "There is no northpole!"
        end if
        if(abs(xprobmax2-dpi)<smalldouble) then
          if(mype==0) write(unitterm,*)&
              "Will apply pi-periodic conditions at southpole!"
          poleB(2,2)=.true.
        else
          if(mype==0) write(unitterm,*) "There is no southpole!"
        end if
      end if
    case (cylindrical)
#:if defined('LOG_RADIUS')
      ! With log_r0 = 0 the radial mesh is uniform in ln(r) and can never
      ! reach r = 0, so a logCylindrical domain has no axis. Running the tests
      ! below would then be actively wrong rather than merely redundant:
      ! xprobmin1 holds ln(r_min), which is exactly zero for the perfectly
      ! ordinary domain r_min = 1, and an axis would be reported where there
      ! is none. A user who asks for typeboundary_min1='pole' anyway is caught
      ! by check_pole_setup, which sees poleB_requested without poleB.
      !
      ! With log_r0 > 0 the logical coordinate is ln(1 + r/r0), so xprobmin1
      ! is zero if and only if r_min is - the tests below are then not merely
      ! safe but exact, and are used unchanged.
      if (log_r0 <= zero) return
#:endif
      
      if (1 == phi_ .and. periodB(1)) then
        if(mod(ng1(1),2)/=0) then
          call mpistop("Number of meshes in phi-direction should be even!")
        end if

        if(abs(xprobmin1)<smalldouble) then
          if (mype==0) then
            write(unitterm,*) "Will apply pi-periodic conditions at r=0"
          end if
          poleB(1,1)=.true.
        else
          if (mype==0) then
            write(unitterm,*) "There is no cylindrical axis!"
          end if
        end if
      end if
      
      if (2 == phi_ .and. periodB(2)) then
        if(mod(ng2(1),2)/=0) then
          call mpistop("Number of meshes in phi-direction should be even!")
        end if

        if(abs(xprobmin1)<smalldouble) then
          if (mype==0) then
            write(unitterm,*) "Will apply pi-periodic conditions at r=0"
          end if
          poleB(1,1)=.true.
        else
          if (mype==0) then
            write(unitterm,*) "There is no cylindrical axis!"
          end if
        end if
      end if
      
      if (3 == phi_ .and. periodB(3)) then
        if(mod(ng3(1),2)/=0) then
          call mpistop("Number of meshes in phi-direction should be even!")
        end if

        if(abs(xprobmin1)<smalldouble) then
          if (mype==0) then
            write(unitterm,*) "Will apply pi-periodic conditions at r=0"
          end if
          poleB(1,1)=.true.
        else
          if (mype==0) then
            write(unitterm,*) "There is no cylindrical axis!"
          end if
        end if
      end if
    end select

  end subroutine set_pole

  !> Validate a configured polar axis, once, on the host.
  !>
  !> Upstream MPI-AMRVAC defers most of this to the copy itself: pole_copy
  !> calls mpistop("Pole boundary condition should be symm or asymm") when it
  !> meets a typeboundary entry it cannot use.  Here the equivalent branches
  !> live inside OpenACC kernels in getbc, where mpistop is unreachable and a
  !> bad entry would simply be treated as symmetric.  So everything that could
  !> go wrong is checked here, before any of it is relied on.
  !>
  !> The checks also cover a mismatch that upstream cannot detect at all: where
  !> an axis is (poleB, derived from the domain bounds in set_pole) and where
  !> the user said one is (poleB_requested, recorded when read_par_files
  !> expanded typeboundary='pole') are decided independently, so a domain that
  !> reaches theta=0 without a 'pole' face - or a 'pole' face on a domain that
  !> does not reach the axis - is otherwise silently wrong.
  subroutine check_pole_setup
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer          :: idim, iside, iB, iw
    double precision :: phi_len

    do idim=1,ndim
      do iside=1,2
        if (poleB(iside,idim) .and. .not.poleB_requested(iside,idim)) then
          if (mype==0) write(unitterm,*) "dimension",idim," side",iside,&
             " reaches the axis but was not given typeboundary='pole'"
          call mpistop("a domain boundary lies on the axis but was not &
             &declared 'pole'")
        end if
        if (poleB_requested(iside,idim) .and. .not.poleB(iside,idim)) then
          if (mype==0) write(unitterm,*) "dimension",idim," side",iside,&
             " was given typeboundary='pole' but set_pole found no axis there"
          call mpistop("typeboundary='pole' on a face that is not on the &
             &axis; the phi direction has to be periodic, the number of level-1&
             & blocks across it even, and the face itself at theta=0, theta=pi &
             &or r=0")
        end if
      end do
    end do

    if (.not.any(poleB)) return

#:if PHYS not in ('hd', 'mhd', 'srhd', 'srmhd', 'ffhd')
    ! The sign table the pole copy reads is built from iw_mom/iw_mag in
    ! read_par_files, which describes hd, mhd, srhd and srmhd (the
    ! relativistic modules' mom(:) is the spatial four-velocity u^i=lfac*v^i,
    ! a genuine vector under the pole's pi-rotation since lfac is
    ! rotation-invariant - see the case select in read_par_files).  ffhd is covered too: its only conserved quantities are
    ! scalars (rho, the field-aligned momentum m_par, the energy), all
    ! symmetric across the axis, and its frozen field is not exchanged through
    ! getbc at all - fill_nwextra_device re-derives it analytically in
    ! every ghost cell - so it needs no entry in the sign table.
    call mpistop("polar-axis treatment is only implemented for phys='hd', &
       &phys='mhd', phys='srhd', phys='srmhd' and phys='ffhd'")
#:endif

    if (stagger_grid) call mpistop("polar-axis treatment does not support &
       &stagger_grid")

    ! The pi-periodic neighbour is found by shifting the phi block index by
    ! ng(1)/2, which is half a revolution only if phi spans a full one.  A
    ! half-turn domain with periodic boundaries passes every check upstream
    ! makes and then quietly shifts by the wrong angle.
    select case (phi_)
    case (1)
      phi_len = xprobmax1 - xprobmin1
    case (2)
      phi_len = xprobmax2 - xprobmin2
    case (3)
      phi_len = xprobmax3 - xprobmin3
    case default
      call mpistop("a pole needs a phi direction")
    end select
    if (dabs(phi_len - two*dpi) > 1.0d-12*dpi) then
      if (mype==0) write(unitterm,*) "phi spans",phi_len," rather than 2*pi"
      call mpistop("a polar axis requires the phi direction to span a full &
         &turn (xprobmax-xprobmin = 1 in the 2*pi units of &meshlist)")
    end if

    ! getbc turns the per-variable entry on a pole face into a plus or a minus
    ! sign and has no third option, so anything else here would be read as
    ! symmetric.  read_par_files writes only symm/asymm for a 'pole' face; this
    ! guards against a later edit that widens it.
    do idim=1,ndim
      do iside=1,2
        if (.not.poleB(iside,idim)) cycle
        iB=2*idim-2+iside
        do iw=1,nwflux+nwaux
          if (typeboundary(iw,iB)/=bc_symm .and. typeboundary(iw,&
             iB)/=bc_asymm) then
            if (mype==0) write(unitterm,*) "variable",iw," on boundary",iB,&
               " has boundary type",typeboundary(iw,iB)
            call mpistop("every variable on a pole face must end up symm or &
               &asymm")
          end if
        end do
      end do
    end do

  end subroutine check_pole_setup

  !> Return cell-centre positions from the device to the host.
  !>
  !> fill_geometry_device builds a block's geometry on the GPU and leaves it
  !> there. Positions are the part host code asks for most often - the initial
  !> condition, the user hooks and the output all read ps(igrid)%x - but it
  !> asks rarely compared to how often blocks are (re)allocated, so the copy
  !> back does not belong in alloc_node: with AMR that would pay for it on
  !> every regrid, for blocks nothing on the host is going to look at.
  !>
  !> Call this instead immediately before host code that reads ps(igrid)%x,
  !> passing igrid when it is about to read a single block. Every such site is
  !> listed in CLAUDE.md; a missed one reads whatever the host copy held for
  !> the block that previously occupied that slot.
  subroutine sync_positions_host(igrid)
    use mod_global_parameters

    integer, intent(in), optional :: igrid

    integer :: iigrid

    if (present(igrid)) then
       ${GPU_UPDATE_HOST('bgeo%x(:,:,:,:,igrid), bgeoc%x(:,:,:,:,igrid)')}$
    else
       do iigrid = 1, igridstail
          ${GPU_UPDATE_HOST('bgeo%x(:,:,:,:,igrids(iigrid)), bgeoc%x(:,:,:,:,igrids(iigrid))')}$
       end do
    end if

  end subroutine sync_positions_host

  !> Return a block's whole geometry - positions and metrics alike - from the
  !> device to the host, for every grid the host is about to read.
  !>
  !> This is the output-time counterpart of sync_positions_host: saveamrfile
  !> calls it in the same place it pulls back bg%w, which covers the log, the
  !> snapshot, the slices, the collapsed views and autoconvert, and the
  !> standalone convert path does the same before generate_plotfile.
  !>
  !> The fetch is unconditional and covers the whole grid in one pass, on
  !> purpose. Remembering "the host is already up to date" in a flag cannot be
  !> made honest, because the geometry is per block: any such flag would have
  !> to be cleared after refreshing whatever set of blocks happened to be
  !> listed at the time, and the next regrid or selectgrids can put a block
  !> that was not in that set straight back in front of host code. Doing every
  !> grid at once also keeps this to a few update directives per block rather
  !> than paying the per-call overhead once per reader.
  subroutine sync_geometry_host()
    use mod_global_parameters

    integer :: iigrid

    ! Index igrids(iigrid) directly rather than through a local igrid: a
    ! variable whose only uses are inside GPU directives looks dead to the
    ! Fortran front end, the assignment to it is dropped, and the update then
    ! silently targets the wrong block.  sync_positions_host below is written
    ! the same way for the same reason.
    do iigrid = 1, igridstail
       ${GPU_UPDATE_HOST('bgeo%x(:,:,:,:,igrids(iigrid))')}$
       ${GPU_UPDATE_HOST('bgeoc%x(:,:,:,:,igrids(iigrid))')}$
#:if GEOM != 'Cartesian'
       ! A Cartesian build has no metrics to fetch - they are not allocated -
       ! so there this reduces to exactly sync_positions_host.
       ${GPU_UPDATE_HOST('bgeo%dvolume(:,:,:,igrids(iigrid)), bgeo%surfaceC(:,:,:,:,igrids(iigrid)), bgeo%ds(:,:,:,:,igrids(iigrid)), bgeo%dx(:,:,:,:,igrids(iigrid))')}$
       ${GPU_UPDATE_HOST('bgeoc%dvolume(:,:,:,igrids(iigrid)), bgeoc%surfaceC(:,:,:,:,igrids(iigrid)), bgeoc%ds(:,:,:,:,igrids(iigrid)), bgeoc%dx(:,:,:,:,igrids(iigrid))')}$
#:endif
    end do

  end subroutine sync_geometry_host

  !> Calculate curl of a vector qvec within ixL
  !> Options to
  !>        employ standard second order CD evaluations
  !>        use Gauss theorem for non-Cartesian grids
  !>        use Stokes theorem for non-Cartesian grids
  subroutine curlvector(qvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,curlvec,idirmin,idirmin0,&
     ndir0,fourthorder)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    integer, intent(in)             :: ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    double precision, intent(in)    :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir0)
    double precision, intent(inout) :: curlvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,idirmin0:3)
    logical, intent(in), optional   :: fourthorder !< Default: false

    integer          :: ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,&
       ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,jxCmin1,jxCmin2,jxCmin3,&
       jxCmax1,jxCmax2,jxCmax3,idir,jdir,kdir,hxOmin1,hxOmin2,hxOmin3,hxOmax1,&
       hxOmax2,hxOmax3,jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3,kxOmin1,&
       kxOmin2,kxOmin3,kxOmax1,kxOmax2,kxOmax3,gxOmin1,gxOmin2,gxOmin3,gxOmax1,&
       gxOmax2,gxOmax3
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       xC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       surface(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    logical          :: use_4th_order

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin:3)
    ! Determine exact value of idirmin while doing the loop.

    ! The Stokesbased discretisation divides by block%surface, the cell-centre
    ! face areas.  Those are no longer built: no other routine reads them, and
    ! carrying a block-sized array per direction for 1:max_blocks purely for
    ! this branch is not worth the memory.  Fail here rather than divide by an
    ! unassociated pointer.
    if (type_curl == Stokesbased) call mpistop(&
       "curlvector: typecurl='Stokesbased' needs block%surface, which is no &
       &longer built - see geo_t in mod_physicaldata")

    use_4th_order = .false.
    if (present(fourthorder)) use_4th_order = fourthorder

    if (use_4th_order) then
      if (.not. slab_uniform) call mpistop(&
         "curlvector: 4th order only supported for slab geometry")
      ! Fourth order, stencil width is two
      ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmin3=ixOmin3-2;ixAmax1=ixOmax1+2
      ixAmax2=ixOmax2+2;ixAmax3=ixOmax3+2;
    else
      ! Second order, stencil width is one
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
      ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    end if

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2.or.ixImin3>ixAmin3.or.ixImax3<ixAmax3) call &
       mpistop("Error in curlvector: Non-conforming input limits")

    idirmin=4
    curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idirmin0:3)=zero

    ! all non-Cartesian cases
    select case(coordinate)
      case(Cartesian) ! Cartesian grids
        invdx=1.d0/dxlevel
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            if (.not. use_4th_order) then
              ! Use second order scheme
              jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
              jxOmin3=ixOmin3+kr(jdir,3);jxOmax1=ixOmax1+kr(jdir,1)
              jxOmax2=ixOmax2+kr(jdir,2);jxOmax3=ixOmax3+kr(jdir,3);
              hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
              hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
              hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
              tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)=half*(qvec(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                 jxOmin3:jxOmax3,kdir) - qvec(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                 hxOmin3:hxOmax3,kdir))*invdx(jdir)
            else
              ! Use fourth order scheme
              kxOmin1=ixOmin1+2*kr(jdir,1);kxOmin2=ixOmin2+2*kr(jdir,2)
              kxOmin3=ixOmin3+2*kr(jdir,3);kxOmax1=ixOmax1+2*kr(jdir,1)
              kxOmax2=ixOmax2+2*kr(jdir,2);kxOmax3=ixOmax3+2*kr(jdir,3);
              jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
              jxOmin3=ixOmin3+kr(jdir,3);jxOmax1=ixOmax1+kr(jdir,1)
              jxOmax2=ixOmax2+kr(jdir,2);jxOmax3=ixOmax3+kr(jdir,3);
              hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
              hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
              hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
              gxOmin1=ixOmin1-2*kr(jdir,1);gxOmin2=ixOmin2-2*kr(jdir,2)
              gxOmin3=ixOmin3-2*kr(jdir,3);gxOmax1=ixOmax1-2*kr(jdir,1)
              gxOmax2=ixOmax2-2*kr(jdir,2);gxOmax3=ixOmax3-2*kr(jdir,3);
              tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)=(-qvec(kxOmin1:kxOmax1,kxOmin2:kxOmax2,&
                 kxOmin3:kxOmax3,kdir) + 8.0d0 * qvec(jxOmin1:jxOmax1,&
                 jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                 kdir) - 8.0d0 * qvec(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                 hxOmin3:hxOmax3,kdir) + qvec(gxOmin1:gxOmax1,gxOmin2:gxOmax2,&
                 gxOmin3:gxOmax3,kdir))/(12.0d0 * dxlevel(jdir))
            end if
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
            else
              curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case(Cartesian_stretched) ! stretched Cartesian grids
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            select case(type_curl)
              case(central)
                tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                   ixAmin3:ixAmax3)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                   ixAmin3:ixAmax3,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
                jxOmin3=ixOmin3+kr(jdir,3);jxOmax1=ixOmax1+kr(jdir,1)
                jxOmax2=ixOmax2+kr(jdir,2);jxOmax3=ixOmax3+kr(jdir,3);
                ! second order centered differencing
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                   jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   hxOmin3:hxOmax3))/(block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                   jxOmin3:jxOmax3,jdir)-block%x(hxOmin1:hxOmax1,&
                   hxOmin2:hxOmax2,hxOmin3:hxOmax3,jdir))
              case(Gaussbased)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                jxCmin3=ixCmin3+kr(jdir,3);jxCmax1=ixCmax1+kr(jdir,1)
                jxCmax2=ixCmax2+kr(jdir,2);jxCmax3=ixCmax3+kr(jdir,3);
                tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3)=block%surfaceC(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir)*(qvec(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3,jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                   jxCmin3:jxCmax3,kdir)-qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3,kdir))/(block%x(jxCmin1:jxCmax1,&
                   jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
                   jdir)-block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3,jdir)))
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)=(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   hxOmin3:hxOmax3))/block%dvolume(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,ixOmin3:ixOmax3)
              case(Stokesbased)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                jxCmin3=ixCmin3+kr(jdir,3);jxCmax1=ixCmax1+kr(jdir,1)
                jxCmax2=ixCmax2+kr(jdir,2);jxCmax3=ixCmax3+kr(jdir,3);
                if(kdir<=ndim)then
                  tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3)=block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,kdir)*(qvec(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                     kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3,jdir)*(qvec(jxCmin1:jxCmax1,&
                     jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
                     kdir)-qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3,kdir))/(block%x(jxCmin1:jxCmax1,&
                     jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
                     jdir)-block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3,jdir)))
                else
                  tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3)=(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                     jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                     jxCmin3:jxCmax3,kdir)-qvec(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                     kdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                     jxCmin3:jxCmax3,jdir)-block%x(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir)))
                endif
                if(idir<=ndim)then
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)=(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))/block%surface(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
                else ! essentially for 2.5D case, idir=3 and jdir,kdir<=2
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)=(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))/(block%ds(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     jdir)*block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,kdir))
                endif
              case default
                call mpistop('no such curl evaluator')
            end select
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
            else
              curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case(spherical) ! possibly stretched spherical grids
        select case(type_curl)
          case(central) ! ok for any dimensionality
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                   ixAmin3:ixAmax3)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                   ixAmin3:ixAmax3,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
                jxOmin3=ixOmin3+kr(jdir,3);jxOmax1=ixOmax1+kr(jdir,1)
                jxOmax2=ixOmax2+kr(jdir,2);jxOmax3=ixOmax3+kr(jdir,3);
                select case(jdir)
                case(1)
                tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                   ixAmin3:ixAmax3)=tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                   ixAmin3:ixAmax3)*block%x(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                   ixAmin3:ixAmax3,1)
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                   jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   hxOmin3:hxOmax3))/((block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                   jxOmin3:jxOmax3,1)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   hxOmin3:hxOmax3,1))*block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,1))
                    case(2)
                if(idir==1) tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                   ixAmin3:ixAmax3)=tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                   ixAmin3:ixAmax3)*dsin(block%x(ixAmin1:ixAmax1,&
                   ixAmin2:ixAmax2,ixAmin3:ixAmax3,2))
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                   jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   hxOmin3:hxOmax3))/((block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                   jxOmin3:jxOmax3,2)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   hxOmin3:hxOmax3,2))*block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,1))
                if(idir==1) tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)=tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)/dsin(block%x(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
               
                  case(3)
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                   jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   hxOmin3:hxOmax3))/((block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                   jxOmin3:jxOmax3,3)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   hxOmin3:hxOmax3,3))*block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,1)*dsin(block%x(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))
               
                end select
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)+tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                else
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)-tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                endif
                if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case(Gaussbased)
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                jxCmin3=ixCmin3+kr(jdir,3);jxCmax1=ixCmax1+kr(jdir,1)
                jxCmax2=ixCmax2+kr(jdir,2);jxCmax3=ixCmax3+kr(jdir,3);
                tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3)=block%surfaceC(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir)*(qvec(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3,jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                   jxCmin3:jxCmax3,kdir)-qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3,kdir))/(block%x(jxCmin1:jxCmax1,&
                   jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
                   jdir)-block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3,jdir)))
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)=(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   hxOmin3:hxOmax3))/block%dvolume(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)+tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                else
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)-tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric terms
            if(idir==2.and.phi_>0) curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2)+qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,phi_)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,r_)
            
            if(idir==phi_) curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,phi_)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,phi_)-qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,r_) +qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,r_)*dcos(block%x(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))/(block%x(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               r_)*dsin(block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2)))
           
            enddo;
          case(Stokesbased)
            !if(ndim<3) call mpistop("Stokesbased for 3D spherical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
             if(lvc(idir,jdir,kdir)/=0)then
              select case(idir)
              case(1)
                if(jdir<kdir) then
                  ! idir=1,jdir=2,kdir=3
                  !! integral along 3rd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                  hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmin3=ixCmin3+kr(jdir,3);jxCmax1=ixCmax1+kr(jdir,1)
                  jxCmax2=ixCmax2+kr(jdir,2);jxCmax3=ixCmax3+kr(jdir,3);
                  ! qvec(3) at cell interface along 2nd dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir))
                  end if
                  ! 2nd coordinate at cell interface along 2nd dimension
                  xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3)=block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir)
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(dsin(xC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3))*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-dsin(xC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))*tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,kdir)
                  !! integral along 2nd dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmin3=ixOmin3-kr(kdir,3);hxOmax1=ixOmax1-kr(kdir,1)
                  hxOmax2=ixOmax2-kr(kdir,2);hxOmax3=ixOmax3-kr(kdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmin3=ixCmin3+kr(kdir,3);jxCmax1=ixCmax1+kr(kdir,1)
                  jxCmax2=ixCmax2+kr(kdir,2);jxCmax3=ixCmax3+kr(kdir,3);
                  ! qvec(2) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)+(tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)-tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3))*block%dx(&
                     ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     jdir))/block%surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)*block%x(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
                end if
              case(2)
                if(jdir<kdir) then
                  ! idir=2,jdir=1,kdir=3
                  !! integral along 1st dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmin3=ixOmin3-kr(kdir,3);hxOmax1=ixOmax1-kr(kdir,1)
                  hxOmax2=ixOmax2-kr(kdir,2);hxOmax3=ixOmax3-kr(kdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmin3=ixCmin3+kr(kdir,3);jxCmax1=ixCmax1+kr(kdir,1)
                  jxCmax2=ixCmax2+kr(kdir,2);jxCmax3=ixCmax3+kr(kdir,3);
                  ! qvec(1) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
                  !! integral along 3rd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                  hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmin3=ixCmin3+kr(jdir,3);jxCmax1=ixCmax1+kr(jdir,1)
                  jxCmax2=ixCmax2+kr(jdir,2);jxCmax3=ixCmax3+kr(jdir,3);
                  ! qvec(3) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3)=block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir)
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)+(xC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3)*tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3)-xC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3))*dsin(block%x(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir))*block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,kdir))/block%surface(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
                end if
              case(3)
                if(jdir<kdir) then
                  ! idir=3,jdir=1,kdir=2
                  !! integral along 1st dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmin3=ixOmin3-kr(kdir,3);hxOmax1=ixOmax1-kr(kdir,1)
                  hxOmax2=ixOmax2-kr(kdir,2);hxOmax3=ixOmax3-kr(kdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmin3=ixCmin3+kr(kdir,3);jxCmax1=ixCmax1+kr(kdir,1)
                  jxCmax2=ixCmax2+kr(kdir,2);jxCmax3=ixCmax3+kr(kdir,3);
                  ! qvec(1) at cell interface along 2nd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3)-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,jdir)
                  !! integral along 2nd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                  hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmin3=ixCmin3+kr(jdir,3);jxCmax1=ixCmax1+kr(jdir,1)
                  jxCmax2=ixCmax2+kr(jdir,2);jxCmax3=ixCmax3+kr(jdir,3);
                  ! qvec(2) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3)=block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir)
                  if(ndim==3) then
                    surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       ixOmin3:ixOmax3)=block%surface(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
                  else
                    surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       ixOmin3:ixOmax3)=block%x(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                       jdir)*block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       ixOmin3:ixOmax3,kdir)*block%dx(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,ixOmin3:ixOmax3,jdir)
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)+(xC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-xC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3)*tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     kdir))/surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)
                end if
              end select
              if(idir<idirmin)idirmin=idir
             endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
      case(cylindrical) ! possibly stretched cylindrical grids
        select case(type_curl)
          case(central)  ! works for any dimensionality, polar/cylindrical
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                   ixAmin3:ixAmax3)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                   ixAmin3:ixAmax3,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
                jxOmin3=ixOmin3+kr(jdir,3);jxOmax1=ixOmax1+kr(jdir,1)
                jxOmax2=ixOmax2+kr(jdir,2);jxOmax3=ixOmax3+kr(jdir,3);
                if(z_==3.or.z_==-1) then
                  ! Case Polar_2D, Polar_2.5D or Polar_3D, i.e. R,phi,Z
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                     ixAmin3:ixAmax3)=tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                     ixAmin3:ixAmax3)*block%x(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                     ixAmin3:ixAmax3,1) !R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                     jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))/(block%x(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                     1)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3,1))
                  if(idir==z_) tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)=tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,1) !(1/R)*d(R V_phi)/dR
                        case(2)
                  ! handles (1/R)d V_Z/dphi or (1/R)d V_R/dphi
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                     jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))/((block%x(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                     2)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3,2))*block%x(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))
                 
                      case(3)
                  ! handles d V_phi/dZ or d V_R/dZ
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                     jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))/(block%x(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                     3)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3,3))
                 
                  end select
                end if
                if(phi_==3.or.phi_==-1) then
                  ! Case Cylindrical_2D, Cylindrical_2.5D or Cylindrical_3D, i.e. R,Z,phi
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                     ixAmin3:ixAmax3)=tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                     ixAmin3:ixAmax3)*block%x(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                     ixAmin3:ixAmax3,1) !R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                     jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))/(block%x(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                     1)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3,1))
                  if(idir==z_) tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)=tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,1) !(1/R)*d(R V_phi)/dR
                        case(2)
                  ! handles d V_phi/dZ or d V_R/dZ
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                     jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))/(block%x(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                     2)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3,2))
                 
                      case(3)
                  ! handles (1/R)d V_Z/dphi or (1/R)d V_R/dphi
                  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)=(tmp(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
                     jxOmin3:jxOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))/((block%x(jxOmin1:jxOmax1,&
                     jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                     3)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3,3))*block%x(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))
                 
                  end select
                end if
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)+tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                else
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)-tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case(Gaussbased) ! works for any dimensionality, polar/cylindrical
            if(ndim<2) call mpistop&
               ("Gaussbased for 2D, 2.5D or 3D polar or cylindrical only")
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                jxCmin3=ixCmin3+kr(jdir,3);jxCmax1=ixCmax1+kr(jdir,1)
                jxCmax2=ixCmax2+kr(jdir,2);jxCmax3=ixCmax3+kr(jdir,3);
                tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3)=block%surfaceC(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir)*(qvec(ixCmin1:ixCmax1,&
                   ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3,jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                   jxCmin3:jxCmax3,kdir)-qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3,kdir))/(block%x(jxCmin1:jxCmax1,&
                   jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
                   jdir)-block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   ixCmin3:ixCmax3,jdir)))
                tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)=(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                   hxOmin3:hxOmax3))/block%dvolume(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)+tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                else
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)-tmp2(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric term from d e_R/d phi= e_phi for unit vectors e_R, e_phi
            !       but minus sign appears due to R,Z,phi ordering (?)
            ! note that in cylindrical 2D (RZ), phi_ is -1
            ! note that in polar 2D     (Rphi), z_ is -1
            if((idir==phi_.or.(phi_==-1.and.idir==3)).and.z_>0) then
              ! cylindrical
              if(  z_==2) curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,idir)-qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,z_)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,r_)
              ! polar
              if(phi_==2) curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,idir)+qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,z_)/block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,r_)
            endif
            enddo;
          case(Stokesbased)
            if(ndim<3) call mpistop("Stokesbased for 3D cylindrical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
               if(idir==r_) then
                if(jdir==phi_) then
                  ! idir=r,jdir=phi,kdir=z
                  !! integral along z dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                  hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmin3=ixCmin3+kr(jdir,3);jxCmax1=ixCmax1+kr(jdir,1)
                  jxCmax2=ixCmax2+kr(jdir,2);jxCmax3=ixCmax3+kr(jdir,3);
                  ! qvec(z) at cell interface along phi dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,kdir)
                  !! integral along phi dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmin3=ixOmin3-kr(kdir,3);hxOmax1=ixOmax1-kr(kdir,1)
                  hxOmax2=ixOmax2-kr(kdir,2);hxOmax3=ixOmax3-kr(kdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmin3=ixCmin3+kr(kdir,3);jxCmax1=ixCmax1+kr(kdir,1)
                  jxCmax2=ixCmax2+kr(kdir,2);jxCmax3=ixCmax3+kr(kdir,3);
                  ! qvec(phi) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)+(tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)-tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3))*block%x(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)*block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,jdir))/block%surface(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
                end if
               else if(idir==phi_) then
                if(jdir<kdir) then
                  ! idir=phi,jdir=r,kdir=z
                  !! integral along r dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmin3=ixOmin3-kr(kdir,3);hxOmax1=ixOmax1-kr(kdir,1)
                  hxOmax2=ixOmax2-kr(kdir,2);hxOmax3=ixOmax3-kr(kdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmin3=ixCmin3+kr(kdir,3);jxCmax1=ixCmax1+kr(kdir,1)
                  jxCmax2=ixCmax2+kr(kdir,2);jxCmax3=ixCmax3+kr(kdir,3);
                  ! qvec(r) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,jdir)
                  !! integral along z dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                  hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmin3=ixCmin3+kr(jdir,3);jxCmax1=ixCmax1+kr(jdir,1)
                  jxCmax2=ixCmax2+kr(jdir,2);jxCmax3=ixCmax3+kr(jdir,3);
                  ! qvec(z) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)+(tmp(hxOmin1:hxOmax1,&
                     hxOmin2:hxOmax2,hxOmin3:hxOmax3)-tmp(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3))*block%dx(&
                     ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     kdir))/block%surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)
                end if
               else ! idir==z_
                if(jdir<kdir) then
                  ! idir=z,jdir=r,kdir=phi
                  !! integral along r dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmin2=ixOmin2-kr(kdir,2)
                  hxOmin3=ixOmin3-kr(kdir,3);hxOmax1=ixOmax1-kr(kdir,1)
                  hxOmax2=ixOmax2-kr(kdir,2);hxOmax3=ixOmax3-kr(kdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmin2=ixCmin2+kr(kdir,2)
                  jxCmin3=ixCmin3+kr(kdir,3);jxCmax1=ixCmax1+kr(kdir,1)
                  jxCmax2=ixCmax2+kr(kdir,2);jxCmax3=ixCmax3+kr(kdir,3);
                  ! qvec(r) at cell interface along phi dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,kdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3)-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,jdir)
                  !! integral along phi dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
                  hxOmin3=ixOmin3-kr(jdir,3);hxOmax1=ixOmax1-kr(jdir,1)
                  hxOmax2=ixOmax2-kr(jdir,2);hxOmax3=ixOmax3-kr(jdir,3);
                  ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3
                  ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmin2=ixCmin2+kr(jdir,2)
                  jxCmin3=ixCmin3+kr(jdir,3);jxCmax1=ixCmax1+kr(jdir,1)
                  jxCmax2=ixCmax2+kr(jdir,2);jxCmax3=ixCmax3+kr(jdir,3);
                  ! qvec(phi) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3,kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       jdir)*(qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir)-qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir))/(block%x(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,jdir)-block%x(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir))
                  else
                    tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                       ixCmin3:ixCmax3)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                       kdir)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                       jxCmin3:jxCmax3,kdir))
                  end if
                  ! r coordinate at cell interface along r dimension
                  xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3)=block%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                     ixCmin3:ixCmax3,jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     ixCmin2:ixCmax2,ixCmin3:ixCmax3,jdir)
                  curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)=(curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)+(xC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)-xC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3)*tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                     hxOmin3:hxOmax3))*block%dx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     kdir))/block%surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3,idir)
                end if
               end if
               if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
        case default
          call mpistop('not possible to calculate curl')
    end select

  end subroutine curlvector


  !> Calculate gradient of a scalar q within ixL in direction idir
  subroutine gradient(q,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,gradq)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idir
    double precision, intent(in)    :: q(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(inout) :: gradq(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    integer                         :: jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,&
       jxOmax3, hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3

    hxOmin1=ixOmin1-kr(idir,1);hxOmin2=ixOmin2-kr(idir,2)
    hxOmin3=ixOmin3-kr(idir,3);hxOmax1=ixOmax1-kr(idir,1)
    hxOmax2=ixOmax2-kr(idir,2);hxOmax3=ixOmax3-kr(idir,3);
    jxOmin1=ixOmin1+kr(idir,1);jxOmin2=ixOmin2+kr(idir,2)
    jxOmin3=ixOmin3+kr(idir,3);jxOmax1=ixOmax1+kr(idir,1)
    jxOmax2=ixOmax2+kr(idir,2);jxOmax3=ixOmax3+kr(idir,3);
    select case(coordinate)
    case(Cartesian)
      gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=half*(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
         jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
         hxOmin3:hxOmax3))/dxlevel(idir)
    case(Cartesian_stretched,Cartesian_expansion)
      gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
         jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
         hxOmin3:hxOmax3))/(block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
         jxOmin3:jxOmax3,idir)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
         hxOmin3:hxOmax3,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3))/((block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,1)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3,1)))
        
      case(2)
        gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3))/((block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,2)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3,2))*block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1))
       
        
      case(3)
        gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3))/((block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,3)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3,3))*block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1)*dsin(block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,2)))
       
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3))/((block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,phi_)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3,phi_))*block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,r_))
      else
        gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3))/(block%x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,idir)-block%x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3,idir))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

  end subroutine gradient
  
  !> Calculate divergence of a vector qvec within ixL
  subroutine divvector(qvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divq,fourthorder,&
     sixthorder)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    double precision, intent(inout) :: divq(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    logical, intent(in), optional   :: fourthorder !< Default: false
    logical, intent(in), optional   :: sixthorder !< Default: false
    logical                         :: use_4th_order
    logical                         :: use_6th_order
    double precision                :: qC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), invdx(1:ndim)
    integer                         :: jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,&
       jxOmax3, hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, ixCmin1,&
       ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, jxCmin1,jxCmin2,jxCmin3,&
       jxCmax1,jxCmax2,jxCmax3
    integer                         :: idims, ixmin1,ixmin2,ixmin3,ixmax1,&
       ixmax2,ixmax3, gxOmin1,gxOmin2,gxOmin3,gxOmax1,gxOmax2,gxOmax3, kxOmin1,&
       kxOmin2,kxOmin3,kxOmax1,kxOmax2,kxOmax3
    integer                         :: lxOmin1,lxOmin2,lxOmin3,lxOmax1,lxOmax2,&
       lxOmax3, fxOmin1,fxOmin2,fxOmin3,fxOmax1,fxOmax2,fxOmax3

    use_4th_order = .false.
    use_6th_order = .false.
    if (present(fourthorder)) use_4th_order = fourthorder
    if (present(sixthorder))  use_6th_order = sixthorder
    if(use_4th_order .and. use_6th_order) call &
       mpistop("divvector: using 4th and 6th order at the same time")

    if(use_4th_order) then
      if (.not. slab_uniform) call mpistop(&
         "divvector: 4th order only supported for slab geometry")
      ! Fourth order, stencil width is two
      ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmin3=ixOmin3-2;ixmax1=ixOmax1+2
      ixmax2=ixOmax2+2;ixmax3=ixOmax3+2;
    else if(use_6th_order) then
      ! Sixth order, stencil width is three
      if (.not. slab_uniform) call mpistop(&
         "divvector: 6th order only supported for slab geometry")
      ixmin1=ixOmin1-3;ixmin2=ixOmin2-3;ixmin3=ixOmin3-3;ixmax1=ixOmax1+3
      ixmax2=ixOmax2+3;ixmax3=ixOmax3+3;
    else
      ! Second order, stencil width is one
      ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1;ixmax1=ixOmax1+1
      ixmax2=ixOmax2+1;ixmax3=ixOmax3+1;
    end if

    if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2<ixmax2.or.&
       ixImin3>ixmin3.or.ixImax3<ixmax3) call &
       mpistop("Error in divvector: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=0.0d0

    if (slab_uniform) then
      do idims=1,ndim
        if(use_6th_order) then
          lxOmin1=ixOmin1+3*kr(idims,1);lxOmin2=ixOmin2+3*kr(idims,2)
          lxOmin3=ixOmin3+3*kr(idims,3);lxOmax1=ixOmax1+3*kr(idims,1)
          lxOmax2=ixOmax2+3*kr(idims,2);lxOmax3=ixOmax3+3*kr(idims,3);
          kxOmin1=ixOmin1+2*kr(idims,1);kxOmin2=ixOmin2+2*kr(idims,2)
          kxOmin3=ixOmin3+2*kr(idims,3);kxOmax1=ixOmax1+2*kr(idims,1)
          kxOmax2=ixOmax2+2*kr(idims,2);kxOmax3=ixOmax3+2*kr(idims,3);
          jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
          jxOmin3=ixOmin3+kr(idims,3);jxOmax1=ixOmax1+kr(idims,1)
          jxOmax2=ixOmax2+kr(idims,2);jxOmax3=ixOmax3+kr(idims,3);
          hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
          hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
          hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
          gxOmin1=ixOmin1-2*kr(idims,1);gxOmin2=ixOmin2-2*kr(idims,2)
          gxOmin3=ixOmin3-2*kr(idims,3);gxOmax1=ixOmax1-2*kr(idims,1)
          gxOmax2=ixOmax2-2*kr(idims,2);gxOmax3=ixOmax3-2*kr(idims,3);
          fxOmin1=ixOmin1-3*kr(idims,1);fxOmin2=ixOmin2-3*kr(idims,2)
          fxOmin3=ixOmin3-3*kr(idims,3);fxOmax1=ixOmax1-3*kr(idims,1)
          fxOmax2=ixOmax2-3*kr(idims,2);fxOmax3=ixOmax3-3*kr(idims,3);
          divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+(-qvec(fxOmin1:fxOmax1,fxOmin2:fxOmax2,&
             fxOmin3:fxOmax3,idims)+9.d0*qvec(gxOmin1:gxOmax1,gxOmin2:gxOmax2,&
             gxOmin3:gxOmax3,idims)-45.d0*qvec(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3,idims)+qvec(lxOmin1:lxOmax1,lxOmin2:lxOmax2,&
             lxOmin3:lxOmax3,idims)-9.d0*qvec(kxOmin1:kxOmax1,kxOmin2:kxOmax2,&
             kxOmin3:kxOmax3,idims)+45.d0*qvec(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
             jxOmin3:jxOmax3,idims))/(60.d0*dxlevel(idims))
        else if(use_4th_order) then
          ! Use fourth order scheme
          kxOmin1=ixOmin1+2*kr(idims,1);kxOmin2=ixOmin2+2*kr(idims,2)
          kxOmin3=ixOmin3+2*kr(idims,3);kxOmax1=ixOmax1+2*kr(idims,1)
          kxOmax2=ixOmax2+2*kr(idims,2);kxOmax3=ixOmax3+2*kr(idims,3);
          jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
          jxOmin3=ixOmin3+kr(idims,3);jxOmax1=ixOmax1+kr(idims,1)
          jxOmax2=ixOmax2+kr(idims,2);jxOmax3=ixOmax3+kr(idims,3);
          hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
          hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
          hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
          gxOmin1=ixOmin1-2*kr(idims,1);gxOmin2=ixOmin2-2*kr(idims,2)
          gxOmin3=ixOmin3-2*kr(idims,3);gxOmax1=ixOmax1-2*kr(idims,1)
          gxOmax2=ixOmax2-2*kr(idims,2);gxOmax3=ixOmax3-2*kr(idims,3);
          divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+(-qvec(kxOmin1:kxOmax1,kxOmin2:kxOmax2,&
             kxOmin3:kxOmax3,idims)+8.d0*qvec(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
             jxOmin3:jxOmax3,idims)-8.d0*qvec(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3,idims)+qvec(gxOmin1:gxOmax1,gxOmin2:gxOmax2,&
             gxOmin3:gxOmax3,idims))/(12.d0*dxlevel(idims))
        else
          ! Use second order scheme
          jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
          jxOmin3=ixOmin3+kr(idims,3);jxOmax1=ixOmax1+kr(idims,1)
          jxOmax2=ixOmax2+kr(idims,2);jxOmax3=ixOmax3+kr(idims,3);
          hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
          hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
          hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
          divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+half*(qvec(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
             jxOmin3:jxOmax3,idims) - qvec(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3,idims))*invdx(idims)
        end if
      end do
    else
      do idims=1,ndim
        hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
        hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
        hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
        ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3;ixCmax1=ixOmax1
        ixCmax2=ixOmax2;ixCmax3=ixOmax3;
        jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
        jxCmin3=ixCmin3+kr(idims,3);jxCmax1=ixCmax1+kr(idims,1)
        jxCmax2=ixCmax2+kr(idims,2);jxCmax3=ixCmax3+kr(idims,3);
        if(stretched_dim(idims) .and. stretch_uncentered) then
          ! linear interpolation at cell interface along stretched dimension
          qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)=block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,idims)*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,idims)+0.5d0*block%dx(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)*(qvec(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2,jxCmin3:jxCmax3,idims)-qvec(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims))/(block%x(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2,jxCmin3:jxCmax3,idims)-block%x(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,idims)))
        else
          qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)=block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,idims)*half*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,idims)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
             jxCmin3:jxCmax3,idims))
        end if
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+qC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-qC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3)
      end do
      divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    end if
  end subroutine divvector

  !> Calculate divergence of a vector qvec within ixL
  !> using limited extrapolation to cell edges
  subroutine divvectorS(qvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divq)
    use mod_global_parameters
    use mod_limiter

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    double precision, intent(inout)    :: divq(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: qL,qR,dqC,ldq,rdq

    double precision :: invdx(1:ndim)
    integer          :: hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,&
       ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,jxCmin1,jxCmin2,jxCmin3,&
       jxCmax1,jxCmax2,jxCmax3,idims,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
       gxCmin1,gxCmin2,gxCmin3,gxCmax1,gxCmax2,gxCmax3,hxCmin1,hxCmin2,hxCmin3,&
       hxCmax1,hxCmax2,hxCmax3

    ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmin3=ixOmin3-2;ixmax1=ixOmax1+2
    ixmax2=ixOmax2+2;ixmax3=ixOmax3+2;

    if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2<ixmax2.or.&
       ixImin3>ixmin3.or.ixImax3<ixmax3) call &
       mpistop("Error in divvectorS: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    do idims=1,ndim
      hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
      hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
      hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
      ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3;ixCmax1=ixOmax1
      ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
      jxCmin3=ixCmin3+kr(idims,3);jxCmax1=ixCmax1+kr(idims,1)
      jxCmax2=ixCmax2+kr(idims,2);jxCmax3=ixCmax3+kr(idims,3);
      gxCmin1=ixCmin1-kr(idims,1);gxCmin2=ixCmin2-kr(idims,2)
      gxCmin3=ixCmin3-kr(idims,3);gxCmax1=jxCmax1;gxCmax2=jxCmax2
      gxCmax3=jxCmax3;
      hxCmin1=gxCmin1+kr(idims,1);hxCmin2=gxCmin2+kr(idims,2)
      hxCmin3=gxCmin3+kr(idims,3);hxCmax1=gxCmax1+kr(idims,1)
      hxCmax2=gxCmax2+kr(idims,2);hxCmax3=gxCmax3+kr(idims,3);

      qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2,&
         gxCmin3:gxCmax3) = qvec(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
         hxCmin3:hxCmax3,idims)
      qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,&
         gxCmin3:gxCmax3) = qvec(gxCmin1:gxCmax1,gxCmin2:gxCmax2,&
         gxCmin3:gxCmax3,idims)
      dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3)= qR(gxCmin1:gxCmax1,&
         gxCmin2:gxCmax2,gxCmin3:gxCmax3)-qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,&
         gxCmin3:gxCmax3)
      call dwlimiter2(dqC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         gxCmin1,gxCmin2,gxCmin3,gxCmax1,gxCmax2,gxCmax3,idims,&
         type_gradient_limiter(block%level),ldw=ldq,rdw=rdq)
      qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3) + half*ldq(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3) - half*rdq(jxCmin1:jxCmax1,&
         jxCmin2:jxCmax2,jxCmin3:jxCmax3)

      if (slab_uniform) then
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+half*(qR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-qL(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3))*invdx(idims)
      else
        qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,idims)*qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,idims)*qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+qR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-qL(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3)
      end if
    end do
    if(.not.slab_uniform) divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

  end subroutine divvectorS
  
end module mod_geometry
