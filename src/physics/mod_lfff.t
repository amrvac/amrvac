!> Program to extrapolate linear force-free fields in 3D Cartesian coordinates,
!> based on exact Green function method (Chiu & Hilton 1977 ApJ 212,873).
!>
!> Usage:
!> 1 In the subroutine usr_set_parameters of mod_usr.t:
!>  To extrapolate a linear force free field from a observed magnetogram 
!>  prepared in a data file, e.g., 'hmiM720sxxxx.dat' replace 
!>  call init_bc_fff_data('hmiM720sxxxx.dat',unit_length,unit_magneticfield)
!>  'hmiM720sxxxx.dat' must be a binary file containing nx1,nx2,xc1,xc2,dxm1,
!>  dxm2, Bz0(nx1,nx2). Integers nx1 and nx2 give the resolution of the 
!>  uniform-grid magentogram. Others are double-precision floats. xc1 and xc2
!>  are coordinates of the central point of the magnetogram. dxm1 and dxm2 
!>  are the cell sizes for each direction, Bz0 is the vertical conponent 
!>  of magetic field on the solar surface from observations.
!>2 In the subroutine usr_init_one_grid of mod_usr.t,
!>  add lines like:
!>
!>  double precision :: Bf(ixG^S,1:ndir), alpha, zshift
!>
!>  alpha=0.d0     ! potential field
!>  !alpha=0.08d0  ! non-potential linear force-free field
!>  zshift=0.05d0  ! lift your box zshift heigher to the bottom magnetogram
!>  call calc_lin_fff(ixG^L,ix^L,Bf,x,alpha,zshift) 
!>
!>3 Notice that the resolution of input magnetogram must be better than the best
!>  resolution of your AMR grid to have a good behavior close to the bottom layer
module mod_lfff
  implicit none
  
  double precision, save :: Bzmax,darea
  double precision, allocatable, save :: Bz0(:,:)
  double precision, allocatable, save :: xa1(:),xa2(:)
  integer, save :: nx1,nx2
  double precision, parameter :: lfff_mode_tolerance=1.d-12
  double precision, parameter :: lfff_resonance_tolerance=1.d-10
  double precision, parameter :: lfff_flux_balance_tolerance=1.d-8
  
contains

  subroutine init_b_fff_data(magnetogramname,qLunit,qBunit)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    double precision, intent(in) :: qLunit,qBunit
    double precision :: xc1,xc2,dxm1,dxm2
    integer, dimension(MPI_STATUS_SIZE) :: statuss
    integer :: file_handle,i
    logical :: aexist
    character(len=*), intent(in) :: magnetogramname
    ! nx1,nx2 are numbers of cells for each direction
    ! xc1,xc2 are coordinates of the central point of the magnetogram
    ! dxm1,dxm2 are cell sizes for each direction
    ! Bz0 is the 2D Bz magnetogram
    inquire(file=magnetogramname,exist=aexist)
    if(.not. aexist) then
      if(mype==0) write(*,'(2a)') "can not find file:",magnetogramname
      call mpistop("no input magnetogram----init_b_fff_data")
    end if
    call MPI_FILE_OPEN(icomm,magnetogramname,MPI_MODE_RDONLY,MPI_INFO_NULL,&
                       file_handle,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,nx1,1,MPI_INTEGER,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,nx2,1,MPI_INTEGER,statuss,ierrmpi)
    allocate(Bz0(nx1,nx2))
    call MPI_FILE_READ_ALL(file_handle,xc1,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,xc2,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,dxm1,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,dxm2,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,Bz0,nx1*nx2,MPI_DOUBLE_PRECISION,&
                           statuss,ierrmpi)
    call MPI_FILE_CLOSE(file_handle,ierrmpi)
    allocate(xa1(nx1))
    allocate(xa2(nx2))
    do i=1,nx1 
      xa1(i) = xc1 + (dble(i) - dble(nx1)/2.d0 - 0.5d0)*dxm1
    enddo
    do i=1,nx2
      xa2(i) = xc2 + (dble(i) - dble(nx2)/2.d0 - 0.5d0)*dxm2
    enddo
    ! declare and define global variables Lunit and Bunit to be your length unit in
    ! cm and magnetic strength unit in Gauss first
    dxm1=dxm1/qLunit
    dxm2=dxm2/qLunit
    xa1=xa1/qLunit
    xa2=xa2/qLunit
    darea=dxm1*dxm2
    Bz0=Bz0/qBunit
    Bzmax=maxval(dabs(Bz0(:,:)))
    
    ! normalize b
    Bz0=Bz0/Bzmax
    if(mype==0) then
      print*,'magnetogram xrange:',minval(xa1),maxval(xa1)
      print*,'magnetogram yrange:',minval(xa2),maxval(xa2)
    end if
    
    if(mype==0) then
      print*,'extrapolating 3D force-free field from an observed Bz '
      print*,'magnetogram of',nx1,'by',nx2,'pixels. Bzmax=',Bzmax
    endif
  
  end subroutine init_b_fff_data

{^IFTHREED
  subroutine init_b_fff_data_driven_boundary(boundaryname,qLunit,qBunit,qxc1,qxc2)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    use mod_data_driven_boundary, only: read_data_driven_boundary_frame

    character(len=*), intent(in) :: boundaryname
    double precision, intent(in) :: qLunit,qBunit
    double precision, intent(in), optional :: qxc1,qxc2

    double precision :: snapshot_time,dxm1,dxm2,xc1,xc2
    double precision, allocatable :: bframe(:,:,:)
    integer :: i,bnx,bny

    call read_data_driven_boundary_frame(boundaryname,snapshot_time,bnx,bny,dxm1,dxm2,bframe)

    nx1 = bnx
    nx2 = bny
    if(allocated(Bz0)) deallocate(Bz0)
    if(allocated(xa1)) deallocate(xa1)
    if(allocated(xa2)) deallocate(xa2)
    allocate(Bz0(nx1,nx2))
    allocate(xa1(nx1))
    allocate(xa2(nx2))

    ! Python V1 stores dx/dy in km and B in Gauss.
    dxm1 = dxm1*1.d5
    dxm2 = dxm2*1.d5
    if(present(qxc1)) then
      xc1 = qxc1
    else
      xc1 = 0.5d0*(xprobmin1+xprobmax1)*qLunit
    end if
    if(present(qxc2)) then
      xc2 = qxc2
    else
      xc2 = 0.5d0*(xprobmin2+xprobmax2)*qLunit
    end if

    Bz0(:,:) = bframe(:,:,3)
    deallocate(bframe)

    do i=1,nx1
      xa1(i) = xc1 + (dble(i) - dble(nx1)/2.d0 - 0.5d0)*dxm1
    end do
    do i=1,nx2
      xa2(i) = xc2 + (dble(i) - dble(nx2)/2.d0 - 0.5d0)*dxm2
    end do

    dxm1 = dxm1/qLunit
    dxm2 = dxm2/qLunit
    xa1 = xa1/qLunit
    xa2 = xa2/qLunit
    darea = dxm1*dxm2
    Bz0 = Bz0/qBunit
    Bzmax = maxval(dabs(Bz0(:,:)))
    if(Bzmax<=0.d0) call mpistop('zero Bz in data-driven boundary frame')
    Bz0 = Bz0/Bzmax

    if(mype==0) then
      print*,'data-driven boundary frame:',trim(boundaryname)
      print*,'snapshot_time [s]:',snapshot_time
      print*,'magnetogram xrange:',minval(xa1),maxval(xa1)
      print*,'magnetogram yrange:',minval(xa2),maxval(xa2)
      print*,'extrapolating potential field from Bz of',nx1,'by',nx2,'pixels. Bzmax=',Bzmax
    end if
  end subroutine init_b_fff_data_driven_boundary
}

  !> Check a constant-alpha magnetogram and optionally remove its core mean.
  !> Status: 0 accepted unchanged, 1 mean removed, 2 invalid treatment,
  !> 3 invalid threshold, 4 strict-mode imbalance, 5 above auto-balance limit.
  subroutine lfff_balance_bottom_flux(bz,treatment,max_imbalance,&
     imbalance_before,imbalance_after,mean_correction,status)
    double precision, intent(inout) :: bz(:,:)
    character(len=*), intent(in) :: treatment
    double precision, intent(in) :: max_imbalance
    double precision, intent(out) :: imbalance_before,imbalance_after
    double precision, intent(out) :: mean_correction
    integer, intent(out) :: status

    double precision :: unsigned_flux

    mean_correction=0.d0
    unsigned_flux=sum(dabs(bz))
    imbalance_before=dabs(sum(bz))/max(unsigned_flux,tiny(1.d0))
    imbalance_after=imbalance_before
    status=0

    if(max_imbalance<0.d0 .or. max_imbalance>1.d0) then
      status=3
      return
    end if

    select case(trim(adjustl(treatment)))
    case('strict')
      if(imbalance_before>lfff_flux_balance_tolerance) status=4
    case('subtract_mean')
      if(imbalance_before>max_imbalance) then
        status=5
        return
      end if
      if(imbalance_before>lfff_flux_balance_tolerance) then
        mean_correction=sum(bz)/dble(size(bz))
        bz=bz-mean_correction
        unsigned_flux=sum(dabs(bz))
        imbalance_after=dabs(sum(bz))/max(unsigned_flux,tiny(1.d0))
        status=1
      end if
    case default
      status=2
    end select
  end subroutine lfff_balance_bottom_flux

{^IFTHREED
  !> Extrapolate a Cartesian potential field with horizontal Fourier modes.
  !> The bottom magnetogram must match the level-one physical cell centers.
  !> Results are streamed one AMRVAC block layer at a time and written directly
  !> to the distributed cell-centered magnetic variables.
  subroutine extrapolate_potential_fft(iw_b,padding_factor,source_plane_depth,&
     alpha,top_boundary,flux_treatment,max_flux_imbalance)
    use mpi
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    use mod_fft, only: fft_2d_real_imag,fft_size_supported,fft_next_supported,&
       fft_factorization
    use mod_forest, only: tree_root
    use mod_geometry

    integer, intent(in) :: iw_b(3)
    integer, intent(in), optional :: padding_factor
    double precision, intent(in), optional :: source_plane_depth
    double precision, intent(in), optional :: alpha
    character(len=*), intent(in), optional :: top_boundary
    character(len=*), intent(in), optional :: flux_treatment
    double precision, intent(in), optional :: max_flux_imbalance

    double precision, allocatable :: bcore(:,:),spec_r(:,:),spec_i(:,:)
    double precision, allocatable :: work_r(:,:),work_i(:,:),bplane(:,:,:)
    double precision, allocatable :: sendbuf(:),recvbuf(:),kx(:),ky(:)
    integer, allocatable :: sendcounts(:),recvcounts(:),sdispls(:),rdispls(:)
    integer, allocatable :: cursor(:),blockpos(:,:)
    double precision :: dx1,dx2,dx3,tol1,tol2,maxerr,z,fac,source_depth
    double precision :: bmean,bflux,tstart,memory_mb,alpha_fft,alpha2,k2
    double precision :: beta,transfer_b,transfer_d,top_height
    double precision :: kx_der,ky_der,kmin,flux_imbalance,unsigned_flux
    double precision :: flux_imbalance_before,bmean_before,bflux_before
    double precision :: mean_correction,max_imbalance
    double precision, parameter :: fft_memory_limit_mb=2048.d0
    integer :: pad,npx,npy,ip0,jp0,ix0,iy0,starti
    integer :: nb1,nb2,nb3,layer,layer_owner,kg,klocal
    integer :: ig1,ig2,ipe,igrid,ic,ix1,ix2,pos,base,payload,slice_size
    integer :: mode,next1,next2,nrecv,nsend,mode_status,flux_status
    character(len=16) :: top_mode,flux_mode
    character(len=256) :: message
    logical :: top_closed,alpha_nonzero

    if(ndim/=3) call mpistop('FFT potential field requires three dimensions')
    if(coordinate/=Cartesian) call mpistop('FFT potential field requires Cartesian coordinates')
    if(any(stretched_dim)) call mpistop('FFT potential field requires a uniform mesh')
    if(refine_max_level/=1 .or. levmax/=1) &
       call mpistop('FFT potential field v1 requires refine_max_level=1')
    if(stagger_grid) call mpistop('FFT potential field v1 does not support stagger_grid')
    if(.not.allocated(Bz0) .or. .not.allocated(xa1) .or. .not.allocated(xa2)) &
       call mpistop('FFT potential field requires initialized magnetogram data')

    pad=2
    if(present(padding_factor)) pad=padding_factor
    if(pad<1) call mpistop('FFT padding factor must be at least one')

    source_depth=0.d0
    if(present(source_plane_depth)) source_depth=source_plane_depth
    if(source_depth<0.d0) call mpistop('FFT source-plane depth must not be negative')

    alpha_fft=0.d0
    if(present(alpha)) alpha_fft=alpha
    alpha2=alpha_fft**2
    alpha_nonzero=(alpha_fft/=0.d0)

    top_mode='open'
    if(present(top_boundary)) top_mode=trim(adjustl(top_boundary))
    select case(trim(top_mode))
    case('open')
      top_closed=.false.
    case('closed')
      top_closed=.true.
    case default
      call mpistop("FFT top boundary must be 'open' or 'closed'")
    end select

    flux_mode='strict'
    if(present(flux_treatment)) flux_mode=trim(adjustl(flux_treatment))
    select case(trim(flux_mode))
    case('strict','subtract_mean')
      continue
    case default
      call mpistop("LFFF flux treatment must be 'strict' or 'subtract_mean'")
    end select
    max_imbalance=0.1d0
    if(present(max_flux_imbalance)) max_imbalance=max_flux_imbalance
    if(max_imbalance<0.d0 .or. max_imbalance>1.d0) &
       call mpistop('LFFF maximum flux imbalance must be between zero and one')

    dx1=dx(1,1)
    dx2=dx(2,1)
    dx3=dx(3,1)
    top_height=source_depth+dble(domain_nx3)*dx3
    tol1=1.d-8*max(1.d0,dabs(xprobmin1),dabs(xprobmax1),dabs(dx1))
    tol2=1.d-8*max(1.d0,dabs(xprobmin2),dabs(xprobmax2),dabs(dx2))

    ! Locate, by coordinates, the physical magnetogram core. This naturally
    ! removes any symmetric ghost layers without assuming their width.
    ix0=0
    if(nx1>=domain_nx1) then
      do starti=1,nx1-domain_nx1+1
        maxerr=0.d0
        do ix1=1,domain_nx1
          maxerr=max(maxerr,dabs(xa1(starti+ix1-1)-&
             (xprobmin1+(dble(ix1)-0.5d0)*dx1)))
        end do
        if(maxerr<=tol1) then
          ix0=starti
          exit
        end if
      end do
    end if
    iy0=0
    if(nx2>=domain_nx2) then
      do starti=1,nx2-domain_nx2+1
        maxerr=0.d0
        do ix2=1,domain_nx2
          maxerr=max(maxerr,dabs(xa2(starti+ix2-1)-&
             (xprobmin2+(dble(ix2)-0.5d0)*dx2)))
        end do
        if(maxerr<=tol2) then
          iy0=starti
          exit
        end if
      end do
    end if
    if(ix0==0 .or. iy0==0) then
      if(mype==0) then
        write(*,*) 'magnetogram size:',nx1,nx2
        write(*,*) 'required physical FFT size:',domain_nx1,domain_nx2
        write(*,*) 'AMRVAC dx/dy:',dx1,dx2
      end if
      call mpistop('FFT magnetogram centers do not match the level-one physical grid')
    end if

    npx=pad*domain_nx1
    npy=pad*domain_nx2
    if(.not.fft_size_supported(npx) .or. .not.fft_size_supported(npy)) then
      next1=fft_next_supported(npx)
      next2=fft_next_supported(npy)
      write(message,'(a,i0,a,a,a,i0,a,a,a,i0,a,i0)') &
         'unsupported padded FFT size ',npx,' (',trim(fft_factorization(npx)),&
         ') x ',npy,' (',trim(fft_factorization(npy)),&
         '); next supported sizes are ',next1,' x ',next2
      call mpistop(trim(message))
    end if

    if(mod(domain_nx1,block_nx1)/=0 .or. mod(domain_nx2,block_nx2)/=0 .or. &
       mod(domain_nx3,block_nx3)/=0) &
       call mpistop('FFT potential field requires domain_nx divisible by block_nx')
    nb1=domain_nx1/block_nx1
    nb2=domain_nx2/block_nx2
    nb3=domain_nx3/block_nx3
    slice_size=block_nx1*block_nx2*3
    payload=slice_size*block_nx3

    ! Conservative per-rank upper bound: four padded spectral/work planes,
    ! core plus three-component physical plane, and one full block layer in
    ! each MPI send/receive buffer. Refuse an unexpectedly large allocation
    ! instead of relying on the operating system to terminate a rank.
    memory_mb=8.d0*(4.d0*dble(npx)*dble(npy)+4.d0*dble(domain_nx1)*&
       dble(domain_nx2)+6.d0*dble(domain_nx1)*dble(domain_nx2)*&
       dble(block_nx3))/(1024.d0**2)
    if(memory_mb>fft_memory_limit_mb) then
      write(message,'(a,f10.1,a,f10.1,a)') 'FFT potential field needs up to ',&
         memory_mb,' MiB/rank, above the internal limit of ',&
         fft_memory_limit_mb,' MiB; reduce block_nx3 or horizontal size'
      call mpistop(trim(message))
    end if

    allocate(bcore(domain_nx1,domain_nx2))
    bcore=Bz0(ix0:ix0+domain_nx1-1,iy0:iy0+domain_nx2-1)
    bmean_before=sum(bcore)/dble(size(bcore))*Bzmax
    bflux_before=sum(bcore)*Bzmax*dx1*dx2
    unsigned_flux=sum(dabs(bcore))
    flux_imbalance_before=dabs(sum(bcore))/max(unsigned_flux,tiny(1.d0))
    mean_correction=0.d0
    flux_status=0
    if(alpha_nonzero) then
      call lfff_balance_bottom_flux(bcore,flux_mode,max_imbalance,&
         flux_imbalance_before,flux_imbalance,mean_correction,flux_status)
      select case(flux_status)
      case(0,1)
        continue
      case(4)
        if(mype==0) then
          write(*,*) 'FFT LFFF alpha:',alpha_fft
          write(*,*) 'relative bottom flux imbalance:',flux_imbalance_before
          write(*,*) 'strict balance tolerance:',lfff_flux_balance_tolerance
        end if
        call mpistop('constant-alpha FFT LFFF requires a flux-balanced bottom magnetogram')
      case(5)
        if(mype==0) then
          write(*,*) 'FFT LFFF alpha:',alpha_fft
          write(*,*) 'relative bottom flux imbalance:',flux_imbalance_before
          write(*,*) 'maximum automatic-balance imbalance:',max_imbalance
        end if
        call mpistop('LFFF magnetogram is too unbalanced for automatic mean subtraction')
      case default
        call mpistop('invalid LFFF flux-balance configuration')
      end select
    else
      flux_imbalance=flux_imbalance_before
    end if
    bmean=sum(bcore)/dble(npx*npy)*Bzmax
    bflux=sum(bcore)*Bzmax*dx1*dx2
    ip0=(npx-domain_nx1)/2+1
    jp0=(npy-domain_nx2)/2+1

    allocate(spec_r(npx,npy),spec_i(npx,npy),work_r(npx,npy),work_i(npx,npy))
    allocate(kx(npx),ky(npy))
    spec_r=0.d0
    spec_i=0.d0
    if(mype==0) then
      spec_r(ip0:ip0+domain_nx1-1,jp0:jp0+domain_nx2-1)=bcore
      call fft_2d_real_imag(spec_r,spec_i,.false.)
    end if
    call MPI_BCAST(spec_r,npx*npy,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(spec_i,npx*npy,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)

    do ix1=1,npx
      mode=ix1-1
      if(mode>npx/2) mode=mode-npx
      kx(ix1)=2.d0*dpi*dble(mode)/(dble(npx)*dx1)
    end do
    do ix2=1,npy
      mode=ix2-1
      if(mode>npy/2) mode=mode-npy
      ky(ix2)=2.d0*dpi*dble(mode)/(dble(npy)*dx2)
    end do

    kmin=min(2.d0*dpi/(dble(npx)*dx1),2.d0*dpi/(dble(npy)*dx2))
    if(alpha_nonzero .and. .not.top_closed .and. dabs(alpha_fft)>=kmin) then
      if(mype==0) then
        write(*,*) 'FFT LFFF alpha:',alpha_fft
        write(*,*) 'smallest nonzero padded horizontal wavenumber:',kmin
      end if
      call mpistop('open-top FFT LFFF requires abs(alpha) < kmin')
    end if

    ! A finite-height closed LFFF permits oscillatory low-wavenumber modes,
    ! except at the eigenvalues sin(beta*L)=0 where the two Bz boundary
    ! conditions are singular.  Validate these modes before distributed work.
    if(alpha_nonzero .and. top_closed) then
      do ix2=1,npy
        do ix1=1,npx
          k2=kx(ix1)**2+ky(ix2)**2
          if(k2<=0.d0 .or. alpha2<=k2) cycle
          beta=dsqrt(alpha2-k2)
          call lfff_fft_transfer(k2,alpha_fft,0.d0,top_height,.true.,&
             transfer_b,transfer_d,mode_status)
          if(mode_status==2) then
            if(mype==0) then
              write(*,*) 'FFT closed LFFF resonant mode kx,ky:',kx(ix1),ky(ix2)
              write(*,*) 'alpha, beta, top height:',alpha_fft,beta,top_height
            end if
            call mpistop('closed-top FFT LFFF is singular for a horizontal mode')
          end if
        end do
      end do
    end if

    allocate(sendcounts(0:npe-1),recvcounts(0:npe-1))
    allocate(sdispls(0:npe-1),rdispls(0:npe-1),cursor(0:npe-1))
    allocate(blockpos(nb1,nb2))
    tstart=MPI_WTIME()

    if(mype==0) then
      write(*,*) 'FFT potential field physical size:',domain_nx1,domain_nx2,domain_nx3
      write(*,*) 'FFT padded size:',npx,npy,' padding factor:',pad
      write(*,*) 'FFT conservative memory bound [MiB/rank]:',memory_mb
      write(*,*) 'FFT source plane below lower face:',source_depth
      write(*,*) 'FFT alpha and top boundary:',alpha_fft,trim(top_mode)
      if(top_closed) write(*,*) 'FFT closed-top height above source plane:',top_height
      write(*,*) 'FFT magnetogram core starts at:',ix0,iy0
      write(*,*) 'FFT input core mean Bz:',bmean_before,' input net bottom flux:',bflux_before
      if(alpha_nonzero) then
        write(*,*) 'FFT LFFF flux treatment:',trim(flux_mode)
        write(*,*) 'FFT input relative bottom flux imbalance:',flux_imbalance_before
        if(flux_status==1) &
           write(*,*) 'FFT subtracted core mean Bz:',mean_correction*Bzmax
        write(*,*) 'FFT corrected relative bottom flux imbalance:',flux_imbalance
      end if
      write(*,*) 'FFT padded zero-mode mean Bz:',bmean,' net bottom flux:',bflux
    end if

    do layer=1,nb3
      layer_owner=mod(layer-1,npe)
      sendcounts=0
      recvcounts=0
      sdispls=0
      rdispls=0
      blockpos=0

      if(mype==layer_owner) then
        do ig2=1,nb2
          do ig1=1,nb1
            ipe=tree_root(ig1,ig2,layer)%node%ipe
            sendcounts(ipe)=sendcounts(ipe)+payload
          end do
        end do
        do ipe=1,npe-1
          sdispls(ipe)=sdispls(ipe-1)+sendcounts(ipe-1)
        end do
        cursor=sdispls
        do ig2=1,nb2
          do ig1=1,nb1
            ipe=tree_root(ig1,ig2,layer)%node%ipe
            blockpos(ig1,ig2)=cursor(ipe)+1
            cursor(ipe)=cursor(ipe)+payload
          end do
        end do
      end if

      nrecv=0
      do ig2=1,nb2
        do ig1=1,nb1
          if(tree_root(ig1,ig2,layer)%node%ipe==mype) nrecv=nrecv+payload
        end do
      end do
      recvcounts(layer_owner)=nrecv
      nsend=sum(sendcounts)
      allocate(sendbuf(max(1,nsend)),recvbuf(max(1,nrecv)))
      sendbuf=0.d0
      recvbuf=0.d0

      if(mype==layer_owner) then
        allocate(bplane(domain_nx1,domain_nx2,3))
        do klocal=1,block_nx3
          kg=(layer-1)*block_nx3+klocal
          z=(dble(kg)-0.5d0)*dx3+source_depth

          do ic=1,3
            do ix2=1,npy
              do ix1=1,npx
                k2=kx(ix1)**2+ky(ix2)**2
                kx_der=kx(ix1)
                ky_der=ky(ix2)
                ! A first derivative of the self-conjugate Nyquist mode cannot
                ! be represented by a real grid function.  Use the standard
                ! zero derivative symbol for that one mode in each direction.
                if(mod(npx,2)==0 .and. ix1==npx/2+1) kx_der=0.d0
                if(mod(npy,2)==0 .and. ix2==npy/2+1) ky_der=0.d0

                if(k2<=0.d0) then
                  if(alpha_nonzero) then
                    transfer_b=0.d0
                  else
                    ! A potential field retains the padded-domain mean flux.
                    ! For a closed top only the fluctuating modes are closed.
                    transfer_b=1.d0
                  end if
                  transfer_d=0.d0
                else
                  call lfff_fft_transfer(k2,alpha_fft,z,top_height,&
                     top_closed,transfer_b,transfer_d,mode_status)
                  if(mode_status/=0) &
                     call mpistop('invalid FFT LFFF mode reached distributed solve')
                end if

                select case(ic)
                case(1)
                  if(k2>0.d0) then
                    fac=(kx_der*transfer_d-&
                       ky_der*alpha_fft*transfer_b)/k2
                    work_r(ix1,ix2)= fac*spec_i(ix1,ix2)
                    work_i(ix1,ix2)=-fac*spec_r(ix1,ix2)
                  else
                    work_r(ix1,ix2)=0.d0
                    work_i(ix1,ix2)=0.d0
                  end if
                case(2)
                  if(k2>0.d0) then
                    fac=(ky_der*transfer_d+&
                       kx_der*alpha_fft*transfer_b)/k2
                    work_r(ix1,ix2)= fac*spec_i(ix1,ix2)
                    work_i(ix1,ix2)=-fac*spec_r(ix1,ix2)
                  else
                    work_r(ix1,ix2)=0.d0
                    work_i(ix1,ix2)=0.d0
                  end if
                case(3)
                  work_r(ix1,ix2)=transfer_b*spec_r(ix1,ix2)
                  work_i(ix1,ix2)=transfer_b*spec_i(ix1,ix2)
                end select
              end do
            end do
            call fft_2d_real_imag(work_r,work_i,.true.)
            bplane(:,:,ic)=Bzmax*work_r(ip0:ip0+domain_nx1-1,&
                                        jp0:jp0+domain_nx2-1)
          end do

          do ig2=1,nb2
            do ig1=1,nb1
              base=blockpos(ig1,ig2)+(klocal-1)*slice_size
              pos=base
              do ic=1,3
                do ix2=1,block_nx2
                  do ix1=1,block_nx1
                    sendbuf(pos)=bplane((ig1-1)*block_nx1+ix1,&
                       (ig2-1)*block_nx2+ix2,ic)
                    pos=pos+1
                  end do
                end do
              end do
            end do
          end do
        end do
        deallocate(bplane)
      end if

      call MPI_ALLTOALLV(sendbuf,sendcounts,sdispls,MPI_DOUBLE_PRECISION,&
         recvbuf,recvcounts,rdispls,MPI_DOUBLE_PRECISION,icomm,ierrmpi)

      pos=1
      do ig2=1,nb2
        do ig1=1,nb1
          if(tree_root(ig1,ig2,layer)%node%ipe/=mype) cycle
          igrid=tree_root(ig1,ig2,layer)%node%igrid
          do klocal=1,block_nx3
            do ic=1,3
              do ix2=1,block_nx2
                do ix1=1,block_nx1
                  ps(igrid)%w(ixMlo1+ix1-1,ixMlo2+ix2-1,&
                     ixMlo3+klocal-1,iw_b(ic))=recvbuf(pos)
                  pos=pos+1
                end do
              end do
            end do
          end do
        end do
      end do
      deallocate(sendbuf,recvbuf)
    end do

    if(mype==0) write(*,*) 'FFT potential/LFFF extrapolation took:',MPI_WTIME()-tstart,'s'
    deallocate(bcore,spec_r,spec_i,work_r,work_i,kx,ky)
    deallocate(sendcounts,recvcounts,sdispls,rdispls,cursor,blockpos)
  end subroutine extrapolate_potential_fft

  !> Vertical transfer functions for one nonzero horizontal constant-alpha
  !> Fourier mode. transfer_b multiplies Bz0, while transfer_d=-d(transfer_b)/dz
  !> enters the horizontal field. status is zero on success, one for an
  !> oscillatory mode in an open half-space, two for a closed-box resonance,
  !> and three for invalid geometry.
  subroutine lfff_fft_transfer(k2,alpha,z,top_height,top_closed,&
     transfer_b,transfer_d,status)
    double precision, intent(in) :: k2,alpha,z,top_height
    logical, intent(in) :: top_closed
    double precision, intent(out) :: transfer_b,transfer_d
    integer, intent(out) :: status

    double precision :: q2,q,beta,denominator,scale

    transfer_b=0.d0
    transfer_d=0.d0
    status=0
    if(k2<=0.d0 .or. top_height<=0.d0 .or. z<0.d0 .or. z>top_height) then
      status=3
      return
    end if

    q2=k2-alpha**2
    scale=max(1.d0,k2,alpha**2)
    if(.not.top_closed) then
      if(q2< -lfff_mode_tolerance*scale) then
        status=1
        return
      end if
      q=dsqrt(max(0.d0,q2))
      transfer_b=dexp(-q*z)
      transfer_d=q*transfer_b
    else if(q2>lfff_mode_tolerance*scale) then
      q=dsqrt(q2)
      if(q*top_height<50.d0) then
        denominator=dsinh(q*top_height)
        transfer_b=dsinh(q*(top_height-z))/denominator
        transfer_d=q*dcosh(q*(top_height-z))/denominator
      else
        ! Algebraically identical exponential form that cannot overflow for
        ! large q*L when 0 <= z <= L.
        denominator=1.d0-dexp(-2.d0*q*top_height)
        transfer_b=(dexp(-q*z)-dexp(-q*(2.d0*top_height-z)))/denominator
        transfer_d=q*(dexp(-q*z)+dexp(-q*(2.d0*top_height-z)))/denominator
      end if
    else if(q2< -lfff_mode_tolerance*scale) then
      beta=dsqrt(-q2)
      denominator=dsin(beta*top_height)
      if(dabs(denominator)<lfff_resonance_tolerance) then
        status=2
        return
      end if
      transfer_b=dsin(beta*(top_height-z))/denominator
      transfer_d=beta*dcos(beta*(top_height-z))/denominator
    else
      ! Critical q=0 limit of sinh(q*(L-z))/sinh(q*L).
      transfer_b=1.d0-z/top_height
      transfer_d=1.d0/top_height
    end if
  end subroutine lfff_fft_transfer

  subroutine calc_lin_fff(ixI^L,ixO^L,Bf,x,alpha,zshift,idir)
  ! PURPOSE: 
  ! Calculation to determine linear FFF from the field on 
  ! the lower boundary (Chiu and Hilton 1977 ApJ 212,873). 
  ! NOTE: Only works for Cartesian coordinates 
  ! INPUT: Bf,x
  ! OUTPUT: updated b in w 
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    integer, optional, intent(in) :: idir
    double precision, intent(in) :: x(ixI^S,1:ndim),alpha,zshift
    double precision, intent(inout) :: Bf(ixI^S,1:ndir)

    double precision, dimension(ixO^S) :: cos_az,sin_az,zk,bigr,r,r2,r3,cos_ar,sin_ar,g,dgdz
    double precision, dimension(ixO^S) :: dx1,dx2,invr3
    double precision :: twopiinv
    logical :: compute_dir(1:ndir)
    integer :: idim,ixp1,ixp2

    Bf=0.d0
    twopiinv = 0.5d0/dpi*Bzmax*darea
    zk(ixO^S)=x(ixO^S,3)-xprobmin3+zshift

    compute_dir=.true.
    if(present(idir)) then
      compute_dir=.false.
      if(idir>=1 .and. idir<=ndir) compute_dir(idir)=.true.
    end if

    ! For a potential field, the Green-function kernel simplifies exactly to
    ! (dx,dy,z)/r**3. Avoid the trigonometric functions, singular 1/bigr
    ! factors, and intermediate arrays needed by the general linear FFF form.
    if(alpha==0.d0) then
      do ixp2=1,nx2
        do ixp1=1,nx1
          dx1(ixO^S)=x(ixO^S,1)-xa1(ixp1)
          dx2(ixO^S)=x(ixO^S,2)-xa2(ixp2)
          r2(ixO^S)=dx1(ixO^S)**2+dx2(ixO^S)**2+zk(ixO^S)**2
          where(r2(ixO^S)>0.d0)
            invr3(ixO^S)=1.d0/(r2(ixO^S)*dsqrt(r2(ixO^S)))
          elsewhere
            invr3(ixO^S)=0.d0
          end where
          if(compute_dir(1)) Bf(ixO^S,1)=Bf(ixO^S,1)+&
             Bz0(ixp1,ixp2)*dx1(ixO^S)*invr3(ixO^S)
          if(compute_dir(2)) Bf(ixO^S,2)=Bf(ixO^S,2)+&
             Bz0(ixp1,ixp2)*dx2(ixO^S)*invr3(ixO^S)
          if(compute_dir(3)) Bf(ixO^S,3)=Bf(ixO^S,3)+&
             Bz0(ixp1,ixp2)*zk(ixO^S)*invr3(ixO^S)
        end do
      end do
      Bf(ixO^S,:)=Bf(ixO^S,:)*twopiinv
      return
    end if

    ! get cos and sin arrays for a non-zero linear force-free alpha
    cos_az(ixO^S)=dcos(alpha*zk(ixO^S))
    sin_az(ixO^S)=dsin(alpha*zk(ixO^S))
    ! looping Bz0 pixels
    do ixp2=1,nx2
      do ixp1=1,nx1
        bigr(ixO^S)=dsqrt((x(ixO^S,1)-xa1(ixp1))**2+&
                          (x(ixO^S,2)-xa2(ixp2))**2)
        r2=bigr**2+zk**2
        r=dsqrt(r2)
        r3=r**3
        cos_ar=dcos(alpha*r)
        sin_ar=dsin(alpha*r)
        where(bigr/=0.d0)
          bigr=1.d0/bigr
        end where
        where(r/=0.d0)
          r=1.d0/r
        end where
        where(r2/=0.d0)
          r2=1.d0/r2
        end where
        where(r3/=0.d0)
          r3=1.d0/r3
        end where
        g=(zk*cos_ar*r-cos_az)*bigr
        dgdz=(cos_ar*(r-zk**2*r3)-alpha*zk**2*sin_ar*r2+alpha*sin_az)*bigr
        do idim=1,ndim
          if(present(idir)) then
            if(idim/=idir) cycle
          end if
          select case(idim)
          case(1)
            Bf(ixO^S,1)=Bf(ixO^S,1)+Bz0(ixp1,ixp2)*((x(ixO^S,1)-xa1(ixp1))*dgdz(ixO^S)&
                     +alpha*g(ixO^S)*(x(ixO^S,2)-xa2(ixp2)))*bigr(ixO^S)
          case(2)
            Bf(ixO^S,2)=Bf(ixO^S,2)+Bz0(ixp1,ixp2)*((x(ixO^S,2)-xa2(ixp2))*dgdz(ixO^S)&
                     -alpha*g(ixO^S)*(x(ixO^S,1)-xa1(ixp1)))*bigr(ixO^S)
          case(3)
            Bf(ixO^S,3)=Bf(ixO^S,3)+Bz0(ixp1,ixp2)*(zk(ixO^S)*cos_ar(ixO^S)*r3(ixO^S)+alpha*&
                                        zk(ixO^S)*sin_ar(ixO^S)*r2(ixO^S))
          end select
        end do
      end do
    end do
    Bf(ixO^S,:)=Bf(ixO^S,:)*twopiinv

  end subroutine calc_lin_fff

  subroutine get_potential_field_potential(ixI^L,ixO^L,potential,x,zshift)
  ! PURPOSE: 
  ! Calculation scalar potential of potential field given
  ! Bz at photosphere (Schmidt 1964 NASSP). 
  ! NOTE: Only works for Cartesian coordinates 
  ! INPUT: x,zshift
  ! OUTPUT: potential
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim),zshift
    double precision, intent(inout) :: potential(ixI^S)

    double precision :: zk
    integer :: ixp1,ixp2,ix^D

    potential=0.d0
    ! looping Bz0 pixels see equation (2)
    !$OMP PARALLEL DO
    do ix3=ixOmin3,ixOmax3
      zk=x(ixOmin1,ixOmin2,ix3,3)-xprobmin3+zshift
      do ix2=ixOmin2,ixOmax2
        do ix1=ixOmin1,ixOmax1
          do ixp2=1,nx2
            do ixp1=1,nx1
              potential(ix^D)=potential(ix^D)+0.5d0*Bz0(ixp1,ixp2)*darea/&
                 (dpi*dsqrt((x(ix^D,1)-xa1(ixp1))**2+(x(ix^D,2)-xa2(ixp2))**2+zk**2))
            end do
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine get_potential_field_potential

  subroutine get_potential_field_potential_sphere(ixI^L,x,potential,nth,nph,magnetogram,theta,phi,r_sphere)
  ! PURPOSE: 
  ! Calculation scalar potential of potential field given
  ! Bz at photosphere (Schmidt 1964 NASSP). 
  ! NOTE: Only works for spherical coordinates 
  ! OUTPUT: potential
    use mod_global_parameters
    integer, intent(in) :: ixI^L,nth,nph
    real*8, intent(in) :: x(ixI^S,1:ndim)
    ! magnetogram Br on photosphere
    real*8, intent(in) :: magnetogram(nth,nph)
    ! theta and phi grid of the photospheric magnetogram
    real*8, intent(in) :: theta(nth),phi(nph)
    ! radius of photosphere
    real*8, intent(in) :: r_sphere
    real*8, intent(out) :: potential(ixI^S)

    real*8 :: area(nth),dtheta_half,dphi,inv2pi
    integer :: ixp1,ixp2,ix^D

    potential=0.d0
    ! assume uniformly discretized theta and phi
    dtheta_half=0.5d0*(theta(2)-theta(1))
    dphi=phi(2)-phi(1)
    area(1:nth)=2.d0*r_sphere**2*sin(theta(1:nth))*sin(dtheta_half)*sin(dphi)
    inv2pi=1.d0/(2.d0*dpi)

    !$OMP PARALLEL DO
    do ix3=ixImin3,ixImax3
      do ix2=ixImin2,ixImax2
        do ix1=ixImin1,ixImax1
          do ixp2=1,nph
            do ixp1=1,nth
              potential(ix^D)=potential(ix^D)+inv2pi*magnetogram(ixp1,ixp2)*area(ixp1)/&
                dsqrt(x(ix^D,1)**2+r_sphere**2-2.d0*x(ix^D,1)*r_sphere*&
                (dsin(x(ix^D,2))*dsin(theta(ixp1))*dcos(phi(ixp2)-x(ix^D,3))+dcos(x(ix^D,2))*&
                 dcos(theta(ixp1))))
            end do
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine get_potential_field_potential_sphere

  !> get potential magnetic field energy given normal B on all boundaries
  subroutine potential_field_energy_mg(benergy)
    use mod_magnetic_reference_fv, only: magnetic_reference_config,&
       magnetic_reference_result,magnetic_reference_field,&
       solve_magnetic_reference_fv,free_magnetic_reference_field

    real*8, intent(out) :: benergy
    type(magnetic_reference_config) :: config
    type(magnetic_reference_result) :: result
    type(magnetic_reference_field) :: bp

    call solve_magnetic_reference_fv(config,bp,result)
    benergy=result%magnetic_energy
    call free_magnetic_reference_field(bp)

  end subroutine potential_field_energy_mg

  !>  Solve Poisson equation of scalar potential using multigrid solver
  subroutine get_potential_field_potential_mg()
    use mod_magnetic_reference_fv, only: magnetic_reference_config,&
       magnetic_reference_result,magnetic_reference_field,&
       solve_magnetic_reference_fv,free_magnetic_reference_field

    type(magnetic_reference_config) :: config
    type(magnetic_reference_result) :: result
    type(magnetic_reference_field) :: bp

    call solve_magnetic_reference_fv(config,bp,result)
    call free_magnetic_reference_field(bp)

  end subroutine get_potential_field_potential_mg

  !> To set boundary condition on physical boundaries for mg Poisson solver
  subroutine multigrid_bc(box, nc, iv, nb, bc_type, bc)
    use mod_magnetic_reference_fv, only: magnetic_reference_bc
    use mod_multigrid_coupling
    type(mg_box_t), intent(in)    :: box
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv      !< Index of variable
    integer, intent(in)           :: nb      !< number of boundary from 1 to 6 for 3D
    integer, intent(out)          :: bc_type !< Type of b.c.
    ! mg boundary values
    double precision, intent(out) :: bc(nc, nc)
    call magnetic_reference_bc(box,nc,iv,nb,bc_type,bc)

  end subroutine multigrid_bc
}
end module mod_lfff
