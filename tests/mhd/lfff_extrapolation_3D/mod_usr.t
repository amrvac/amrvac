module mod_usr
  use mod_mhd
  use mod_lfff
  implicit none
  ! variables for boundary condition
  double precision, allocatable, save :: Bbc(:,:,:,:),Bbt(:,:,:,:),xbc(:,:,:,:)
  double precision, allocatable, save :: vbc(:,:,:,:),vbt(:,:,:,:)
  integer, save :: nxbc^D,ngl
  ! some global parameters
  double precision :: usr_grav,SRadius,rhob,Tiso,tdstop,lalpha,llift,startpos^D

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_3D")
 
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    ! solar gravity
    usr_grav=-2.74d4*unit_length/unit_velocity**2  
    ! magnetic field strength at the bottom
    Busr=Busr/unit_magneticfield
    ! Solar radius
    SRadius=69.61d0 
    ! bottom density
    rhob=1.d0
    ! uniform temperature
    Tiso=1.d6/unit_time

    ! alpha coefficient for linear force-free field
    lalpha=-0.08d0
    ! lift up domain above bottom magnetogram
    llift=0.4d0
    ! time for bottom driving
    tdstop=70.d0
    if(mype==0) then
      write(*,*)'Simulating 3D flux rope formation for quiescent prominences'
    endif
    ! prepare magnetogram at bottom
    call init_b_fff_usr(70,50)
    ! initiate user's bottom boundary
    call init_usr_bc
  end subroutine initglobaldata_usr
 
  subroutine init_b_fff_usr(qnx1,qnx2)
    double precision :: dxm1,dxm2,delx1,delx2,xo1,xo2,yo1,yo2,coB,xoff,yoff
    double precision :: gzonex,gzoney,yrangeh
    integer :: i1,i2,qnx1,qnx2

    nx1=qnx1
    nx2=qnx2
    allocate(Bz0(nx1,nx2))
    allocate(xa1(nx1))
    allocate(xa2(nx2))
    ! width of ghost zone (outer lane) of the bottom magnetogram
    gzonex=dble(nghostcells)*dx(1,1)
    gzoney=dble(nghostcells)*dx(2,1)
    dxm1=(xprobmax1-xprobmin1+2.d0*gzonex)/dble(nx1)
    dxm2=(xprobmax2-xprobmin2+2.d0*gzoney)/dble(nx2)
    darea=dxm1*dxm2
    do i1=1,nx1
      xa1(i1)=xprobmin1+(dble(i1)-0.5d0)*dxm1-gzonex
    enddo
    do i2=1,nx2
      xa2(i2)=xprobmin2+(dble(i2)-0.5d0)*dxm2-gzoney
    enddo
    xo1=(xprobmax1+xprobmin1)/2.d0
    xo2=xo1
    xoff=5.d0
    yoff=2.5d0
    yrangeh=8.d0
    delx1=3.d0
    delx2=1.3d0
    yo1=(xprobmax2+xprobmin2)/2.d0-yoff
    yo2=(xprobmax2+xprobmin2)/2.d0+yoff
    Bz0=0.d0
    do i2=1,nx2
    do i1=1,nx1
      if(xa2(i2)<yrangeh .and. xa2(i2)>-yrangeh) then
        Bz0(i1,i2)=-Busr*dsin(dpi/yrangeh*xa2(i2))
        if(xa1(i1)>xo1+xoff) &
          Bz0(i1,i2)=Bz0(i1,i2)*dexp(-((xa1(i1)-xo1-xoff)/delx1)**2)
        if(xa1(i1)<xo1-xoff) &
          Bz0(i1,i2)=Bz0(i1,i2)*dexp(-((xa1(i1)-xo1+xoff)/delx1)**2)
          
      endif
    enddo
    enddo
    Bzmax=maxval(dabs(Bz0(:,:)))
    ! normalize B
    Bz0=Bz0/Bzmax
    if(mype==0) then
      print*,'extrapolating 3D force-free field from an analytical Bz '
      print*,'magnetogram of',nx1,'by',nx2,'pixels. Bzmax=',Bzmax
    endif

  end subroutine init_b_fff_usr

  subroutine init_usr_bc
  ! prepare fixed bottom boundary condition at highest resolution
    double precision, allocatable :: tmp(:,:,:),bre(:,:,:)
    double precision, allocatable :: x1bc(:),x2bc(:)
    double precision :: spc,vmax,rhobb,coeffrho,vdr
    integer :: ix^D,ixpemin^D,ixpemax^D,ixB^L,file_handle,amode
    integer :: rnpe,nbk1,nbk2,np1,np2,mod1,mod2
    integer, dimension(MPI_STATUS_SIZE) :: statuss
    character(len=20) :: boundaryfile,bottomflow
    logical :: bfexist

    ! ghost layer for bottom boundary
    ngl=nghostcells
    nxbc3=ngl
    nxbc1=nint((xprobmax1-xprobmin1)/dx(1,refine_max_level))+2*ngl
    nxbc2=nint((xprobmax2-xprobmin2)/dx(2,refine_max_level))+2*ngl
    if(mype==0) print*,'resolution of boundary layer', nxbc^D
    allocate(Bbc(nxbc^D,3))
    allocate(Bbt(nxbc^D,3))
    allocate(vbc(nxbc^D,3))
    allocate(vbt(nxbc^D,3))
    allocate(xbc(nxbc^D,3))
    {do ix^DB=1,nxbc^DB\}
      xbc(ix^D,1)=xprobmin1+(dble(ix1-ngl)-0.5d0)*dx(1,refine_max_level)
      xbc(ix^D,2)=xprobmin2+(dble(ix2-ngl)-0.5d0)*dx(2,refine_max_level)
      xbc(ix^D,3)=xprobmin3+(dble(ix3-ngl)-0.5d0)*dx(3,refine_max_level)
    {enddo\}
    write(boundaryfile,'(a,i2.2,a)') 'Bb_init',refine_max_level,'.dat'
    inquire(file=boundaryfile, exist=bfexist)
    if(bfexist) then
      if(mype==0) then
        call MPI_FILE_OPEN(MPI_COMM_SELF,boundaryfile,MPI_MODE_RDONLY,MPI_INFO_NULL,&
                          file_handle,ierrmpi)
        call MPI_FILE_READ(file_handle,Bbc,nxbc1*nxbc2*nxbc3*3,&
                           MPI_DOUBLE_PRECISION,statuss,ierrmpi)
        call MPI_FILE_CLOSE(file_handle,ierrmpi)
      end if
      if(npe>1) then
        call MPI_BCAST(Bbc,nxbc1*nxbc2*nxbc3*3,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      end if
    else
      Bbc=0.d0
      ! bottom layer with linear force-free field
      nbk1=int(dsqrt(dble(npe)))
      do while(mod(npe,nbk1)/=0)
        nbk1=nbk1-1
      enddo
      nbk2=npe/nbk1
      if(nxbc1>=nxbc2) then
        nbk1=npe/nbk1
        nbk2=npe/nbk1
      endif
      np1=nxbc1/nbk1
      mod1=mod(nxbc1,nbk1)
      np2=nxbc2/nbk2
      mod2=mod(nxbc2,nbk2)
      ixpemin1=mod(mype,nbk1)*np1+1
      if(mod(mype+1,nbk1)==0) then
        ixpemax1=ixpemin1+np1-1+mod1
      else
        ixpemax1=ixpemin1+np1-1
      endif
      ixpemin2=(mype/nbk1)*np2+1
      if(mype/nbk1==nbk2-1) then
        ixpemax2=ixpemin2+np2-1+mod2
      else
        ixpemax2=ixpemin2+np2-1
      endif
      call calc_lin_fff(1,1,1,nxbc1,nxbc2,nxbc3,ixpemin1,ixpemin2,1,ixpemax1,ixpemax2,&
                        nxbc3,Bbc,xbc,lalpha,llift)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Bbc,nxbc1*nxbc2*nxbc3*3,MPI_DOUBLE_PRECISION,&
                         MPI_SUM,icomm,ierrmpi)
      if(mype==0) then
        print*,'bottom boundaries created!'
        amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        call MPI_FILE_OPEN(MPI_COMM_SELF,boundaryfile,amode, &
                           MPI_INFO_NULL,file_handle,ierrmpi)
        call MPI_FILE_WRITE(file_handle,Bbc,nxbc1*nxbc2*nxbc3*3,&
                           MPI_DOUBLE_PRECISION,statuss,ierrmpi)
        call MPI_FILE_CLOSE(file_handle,ierrmpi)
      endif
    endif

    Bbt=Bbc
    ! calculate velocity field of driving motions
    allocate(tmp(nxbc1,nxbc2,nxbc3))
    allocate(bre(nxbc1,nxbc2,nxbc3))
    allocate(x1bc(nxbc1))
    allocate(x2bc(nxbc2))
    vbc=0.d0
    tmp=0.d0
    ixBmin1=ngl+1;ixBmax1=nxbc1-ngl;
    ixBmin2=ngl+1;ixBmax2=nxbc2-ngl;
    ixBmin3=1;ixBmax3=nxbc3
    ! b3=|B3|/max(|B3|)
    bre(:,:,nxbc3)=dabs(Bbc(:,:,nxbc3,3))/maxval(dabs(Bbc(:,:,nxbc3,3)))
    vdr=0.2d0 !< amplitude of driving velocity
    tmp(:,:,nxbc3)=vdr*bre(:,:,nxbc3)
    ! vy = partial b3 / partial y * exp(-(y/3.5)^2)
    vbc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmax3,2)=-dexp(-(xbc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmax3,2)/8.d0)**2)*&
    (tmp(ixBmin1:ixBmax1,ixBmin2+1:ixBmax2+1,ixBmax3)-tmp(ixBmin1:ixBmax1,ixBmin2-1:ixBmax2-1,ixBmax3))/dx(2,refine_max_level)*0.5d0
    ! cut-off
    do ix1=ixBmin1,ixBmax1
    do ix2=ixBmin2,ixBmax2
      if(xbc(ix1,ix2,ixBmax3,2) .gt. 4.d0) then
        vbc(ix1,ix2,ixBmax3,2)=zero
      else if(xbc(ix1,ix2,ixBmax3,2) .lt. -4.d0) then
        vbc(ix1,ix2,ixBmax3,2)=zero
      endif
    enddo
    enddo
    do ix3=1,ixBmax3-1
      vbc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ix3,2)=vbc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmax3,2)
    enddo
    ! vx=-vy
    vbc(ixB^S,1)=-vbc(ixB^S,2)
    vbt=vbc
    tmp=0.d0
    tmp(:,:,nxbc3)=Bbc(:,:,nxbc3,3)
    do ix1=1,nxbc1
      x1bc(ix1)=xprobmin1+(dble(ix1-ngl)-0.5d0)*dx(1,refine_max_level)
    enddo
    do ix2=1,nxbc2
      x2bc(ix2)=xprobmin2+(dble(ix2-ngl)-0.5d0)*dx(2,refine_max_level)
    enddo
    if(mype==0) then
      amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      write(bottomflow,'(a,i2.2,a)') 'vbottom',iprob,'.dat'
      call MPI_FILE_OPEN(MPI_COMM_SELF,bottomflow,amode, &
                         MPI_INFO_NULL,file_handle,ierrmpi)
      call MPI_FILE_WRITE(file_handle,nxbc1,1,MPI_INTEGER,statuss,ierrmpi)
      call MPI_FILE_WRITE(file_handle,nxbc2,1,MPI_INTEGER,statuss,ierrmpi)
      call MPI_FILE_WRITE(file_handle,x1bc,nxbc1,&
                          MPI_DOUBLE_PRECISION,statuss,ierrmpi)
      call MPI_FILE_WRITE(file_handle,x2bc,nxbc2,&
                          MPI_DOUBLE_PRECISION,statuss,ierrmpi)
      call MPI_FILE_WRITE(file_handle,vbc(:,:,nxbc3,1:2),nxbc1*nxbc2*2,&
                          MPI_DOUBLE_PRECISION,statuss,ierrmpi)
      call MPI_FILE_WRITE(file_handle,tmp(:,:,nxbc3),nxbc1*nxbc2,&
                          MPI_DOUBLE_PRECISION,statuss,ierrmpi)
      call MPI_FILE_CLOSE(file_handle,ierrmpi)
    endif
    if(mype==0) write(*,*) 'Bottom Bz ranging:',maxval(Bbc(:,:,nxbc3,3)),minval(Bbc(:,:,nxbc3,3))
    tmp(:,:,nxbc3)=dsqrt(vbc(:,:,nxbc3,1)**2+vbc(:,:,nxbc3,2)**2)
    if(mype==0) write(*,*) 'max driven velocity:',maxval(tmp(:,:,nxbc3)),'v0:',vdr
    coeffrho=usr_grav*SRadius**2/Tiso
    rhobb=rhob*dexp(coeffrho*(1.d0/SRadius-1.d0/((xprobmin3+0.5d0*dx(3,refine_max_level))+SRadius)))
    if(mype==0) write(*,*) 'max alfven speed:',maxval(dsqrt(Bbc(:,:,nxbc3,1)**2+&
    Bbc(:,:,nxbc3,2)**2+Bbc(:,:,nxbc3,3)**2))/dsqrt(rhobb),'with density',rhobb
    if(mype==0) write(*,*) 'max drive/alfven speed:',maxval(tmp(:,:,nxbc3)/(dsqrt(Bbc(:,:,nxbc3,1)**2+&
    Bbc(:,:,nxbc3,2)**2+Bbc(:,:,nxbc3,3)**2))*dsqrt(rhobb))
    deallocate(x1bc)
    deallocate(x2bc)
    deallocate(tmp)
    call MPI_BARRIER(icomm,ierrmpi)
    ! start position of the global bottom boundary layer
    startpos^D=xprobmin^D-dble(nghostcells)*dx(^D,refine_max_level);\

  end subroutine init_usr_bc

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bf(ixI^S,1:ndir)
    logical, save :: first=.true.
    if(mype==0 .and. first) then
      write(*,*)'initializing grids ...'
      first=.false.
    endif
    ! use green function method to extrapolate B
    call calc_lin_fff(ixI^L,ixO^L,Bf,x,lalpha,llift)
    w(ixO^S,mag(:))=Bf(ixO^S,1:3)
    w(ixO^S,mom(:))=zero
    w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                  (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
    if(mhd_glm) w(ixO^S,psi_)=zero
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_physics
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: ft,tfstop,tramp1,tramp2,coeffrho,vlimit,vsign
    double precision :: xlen^D
    integer :: ix^D,ixbc^D

    select case(iB)
    case(1)
      do ix1=ixOmin1,ixOmax1
        w(ix1^%1ixO^S,rho_)=w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(1))=w(ixOmax1+1^%1ixO^S,mom(1))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(2))=w(ixOmax1+1^%1ixO^S,mom(2))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(3))=w(ixOmax1+1^%1ixO^S,mom(3))/w(ixOmax1+1^%1ixO^S,rho_)
      enddo
      do ix1=ixOmax1,ixOmin1,-1
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                   (-w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
             +4.0d0*w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
      enddo
      call phys_to_conserved(ixI^L,ixO^L,w,x)
    case(2)
      do ix1=ixOmin1,ixOmax1
        w(ix1^%1ixO^S,rho_)=w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(1))=w(ixOmin1-1^%1ixO^S,mom(1))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(2))=w(ixOmin1-1^%1ixO^S,mom(2))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(3))=w(ixOmin1-1^%1ixO^S,mom(3))/w(ixOmin1-1^%1ixO^S,rho_)
      enddo
      do ix1=ixOmin1,ixOmax1
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                   (-w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
              +4.0d0*w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
      enddo
      call phys_to_conserved(ixI^L,ixO^L,w,x)
    case(3)
      do ix2=ixOmin2,ixOmax2
        w(ix2^%2ixO^S,rho_)=w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(1))=w(ixOmax2+1^%2ixO^S,mom(1))/w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(2))=w(ixOmax2+1^%2ixO^S,mom(2))/w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(3))=w(ixOmax2+1^%2ixO^S,mom(3))/w(ixOmax2+1^%2ixO^S,rho_)
      enddo
      do ix2=ixOmax2,ixOmin2,-1
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                   (-w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,mag(1):mag(3)) &
              +4.0d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,mag(1):mag(3)))
      enddo
      call phys_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      do ix2=ixOmin2,ixOmax2
        w(ix2^%2ixO^S,rho_)=w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(1))=w(ixOmin2-1^%2ixO^S,mom(1))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(2))=w(ixOmin2-1^%2ixO^S,mom(2))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(3))=w(ixOmin2-1^%2ixO^S,mom(3))/w(ixOmin2-1^%2ixO^S,rho_)
      enddo
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                   (-w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,mag(1):mag(3)) &
              +4.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(1):mag(3)))
      enddo
      call phys_to_conserved(ixI^L,ixO^L,w,x)
    case(5)
      tramp1=10.d0
      tramp2=10.d0
      tfstop=tdstop-tramp2
      if(qt<tramp1) then
        ft=qt/tramp1
      else if(qt<tfstop) then
        ft=1.d0
      else if(qt<tfstop+tramp2) then
        ft=(tfstop+tramp2-qt)/tramp2
      else
        ft=0.d0
      endif
      if(block%level<refine_max_level) then
        w(ixO^S,mom(:))=0.d0
      else
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
           xlen^D=x({ix^D,},^D)-startpos^D;\
           ixbc^D=ceiling(xlen^D/dx(^D,refine_max_level))\
           w(ix^D,mom(:))=ft*vbt(ixbc^D,:)
        {end do\}
      end if
      do ix3=ixOmax3,ixOmin3,-1
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/11.0d0)* &
             ( +2.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mag(1):mag(3)) &
               -9.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(1):mag(3)) &
              +18.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(1):mag(3)))
      enddo
      coeffrho=usr_grav*SRadius**2/Tiso
      w(ixO^S,rho_)=rhob*dexp(coeffrho*(1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
      call phys_to_conserved(ixI^L,ixO^L,w,x)
    case(6)
      coeffrho=usr_grav*SRadius**2/Tiso
      do ix3=ixOmin3,ixOmax3
        w(ix3^%3ixO^S,rho_)=w(ixOmin3-1^%3ixO^S,rho_)*dexp(coeffrho*(1.d0/(SRadius+x(ixOmin3-1^%3ixO^S,3))-1.d0/(SRadius+x(ix3^%3ixO^S,3))))
      enddo
      do ix3=ixOmin3,ixOmax3
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/11.0d0)* &
              ( +2.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-3,mag(1):mag(3)) &
                -9.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-2,mag(1):mag(3)) &
               +18.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-1,mag(1):mag(3)))
      enddo
      do ix3=ixOmin3,ixOmax3
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mom(1):mom(3))=(1.0d0/11.0d0)* &
             ( +2.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-3,mom(1):mom(3)) &
               -9.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-2,mom(1):mom(3)) &
              +18.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-1,mom(1):mom(3)))
      enddo
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
      vlimit=0.5d0*dsqrt((w(ix^D,mag(1))**2+&
             w(ix^D,mag(2))**2+w(ix^D,mag(3))**2)/w(ix^D,rho_))
      if(dabs(w(ix^D,mom(3)))>vlimit) then
        vsign=sign(one,w(ix^D,mom(3)))
        w(ix^D,mom(3))=vsign*vlimit
      endif
      {end do\}
      call phys_to_conserved(ixI^L,ixO^L,w,x)
    case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,3)))**2
  end subroutine

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,3)=ggrid(ixO^S)
  end subroutine gravity

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (any(x(ixO^S,3)<=xprobmin3+0.3d0)) then
      refine=1
      coarsen=-1
    endif
  end subroutine special_refine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: tmp(ixI^S),dip(ixI^S),divb(ixI^S),B2(ixI^S)
    double precision, dimension(ixI^S,1:ndir) :: Btotal,qvec,curlvec
    integer                            :: ix^D,idirmin,idims,idir,jdir,kdir

    ! Btotal & B^2
    Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+1)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
    ! output divB1
    call get_normalized_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+2)=divb(ixO^S)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+3)=w(ixO^S,rho_)*Tiso*two/B2(ixO^S)
    ! store current
    call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
    do idir=1,ndir
      w(ixO^S,nw+3+idir)=curlvec(ixO^S,idir)
    end do
    ! find magnetic dips
!    dip=0.d0
!    do idir=1,ndir
!      call gradient(w(ixI^S,mag(3)),ixI^L,ixO^L,idir,tmp)
!      dip(ixO^S)=dip(ixO^S)+w(ixO^S,b0_+idir)*tmp(ixO^S)
!    end do
!    where(dabs(w(ixO^S,mag(3)))<0.08d0 .and. dip(ixO^S)>=0.d0)
!      w(ixO^S,nw+8)=1.d0
!    elsewhere
!      w(ixO^S,nw+8)=0.d0
!    end where
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='Alfv divB beta j1 j2 j3'
  end subroutine specialvarnames_output

end module mod_usr
