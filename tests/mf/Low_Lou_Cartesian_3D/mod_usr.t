!===================================================================================
! Purpose: computing nonlinear force-free field (NLFFF) using the magneto-frictional 
!          module with a bottom boundary from Low & Lou solution.
!          Refer to Y. Guo, C. Xia, R. Keppens, and G. Valori, 2016, ApJ, 828, 82
!          and Y. Guo, C. Xia, and R. Keppens, 2016, ApJ, 828, 83 
!          for more information.
!===================================================================================
module mod_usr
  use mod_mhd
  use mod_lfff
  implicit none
  ! variables for boundary condition
  double precision, allocatable, save :: Bx0(:,:),By0(:,:) ! Bz0(:,:) has been declared in mod_lfff
  ! some global parameters
  double precision :: k_B,miu0,mass_H,usr_grav,SRadius,rhob,Tiso,lalpha,llift

contains

  !==============================================================================
  ! Purpose: to include global parameters, set user methods, set coordinate 
  !          system and activate physics module.
  !==============================================================================
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call set_coordinate_system("Cartesian_3D")
    call mhd_activate()
  end subroutine usr_init

  !==============================================================================
  ! Purpose: to initialize user public parameters and reset global parameters.
  !          Input data are also read here.
  !==============================================================================
  subroutine initglobaldata_usr()
    use mod_global_parameters

    ! normalization unit in CGS Unit
    k_B = 1.3806d-16          ! erg*K^-1
    miu0 = 4.d0*dpi           ! Gauss^2 cm^2 dyne^-1
    mass_H = 1.67262d-24      ! g
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3
    unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    unit_time          = unit_length/unit_velocity                     ! 85.8746159942810 s 

    usr_grav=-2.74d4*unit_length/unit_velocity**2  ! solar gravity
    SRadius=6.955d10/unit_length                   ! Solar radius
    rhob=1.0d0                                     ! bottom density
    Tiso=1.0d6/unit_temperature                    ! uniform temperature

    lalpha=0.d0    ! alpha coefficient for linear force-free field
    llift=0.05d0   ! lift up domain above bottom magnetogram by half the cell size
 
    ! prepare magnetogram at bottom, note that potential field only needs Bz0, but NLFFF needs Bx0, By0, and Bz0
    ! init_b_fff_data reads incorrect data for Bz0, but init_b_nlfff_data will repalce them with correct data
    call init_b_fff_data('nlfff_boundary.data',unit_length,unit_magneticfield)
    call init_b_nlfff_data('nlfff_boundary.data',unit_length,unit_magneticfield)
  end subroutine initglobaldata_usr

  !==============================================================================
  ! Purpose: to read vector magnetic field from an input file
  !==============================================================================
  subroutine init_b_nlfff_data(magnetogramname,qLunit,qBunit)
    use mod_global_parameters

    character(len=*), intent(in) :: magnetogramname
    double precision, intent(in) :: qLunit,qBunit
    integer, dimension(MPI_STATUS_SIZE) :: statuss
    integer :: file_handle,i
    double precision :: xc1,xc2,dxm1,dxm2,Bxmax,Bymax
    logical :: aexist
    ! nx1,nx2 are numbers of cells for each direction, which has been declared in mod_lfff
    ! xc1,xc2 are coordinates of the central point of the magnetogram
    ! dxm1,dxm2 are cell sizes for each direction
    ! Bx0, By0, and Bz0 are the 2D vector magnetic field
    inquire(file=magnetogramname,exist=aexist)
    if(.not. aexist) then
      if(mype==0) write(*,'(2a)') "can not find file:",magnetogramname
      call mpistop("no input magnetogram----init_b_nlfff_data")
    end if
    call MPI_FILE_OPEN(icomm,magnetogramname,MPI_MODE_RDONLY,MPI_INFO_NULL,&
                       file_handle,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,nx1,1,MPI_INTEGER,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,nx2,1,MPI_INTEGER,statuss,ierrmpi)

    allocate(Bx0(nx1,nx2))
    allocate(By0(nx1,nx2))
    ! Bz0(nx1,nx2) has been allocated in mod_lfff

    call MPI_FILE_READ_ALL(file_handle,xc1,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,xc2,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,dxm1,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,dxm2,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,Bx0,nx1*nx2,MPI_DOUBLE_PRECISION,&
                           statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,By0,nx1*nx2,MPI_DOUBLE_PRECISION,&
                           statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,Bz0,nx1*nx2,MPI_DOUBLE_PRECISION,&
                           statuss,ierrmpi)
    call MPI_FILE_CLOSE(file_handle,ierrmpi)
    ! declare and define global variables Lunit and Bunit to be your length unit in
    ! cm and magnetic strength unit in Gauss first
    dxm1=dxm1/qLunit
    dxm2=dxm2/qLunit
    Bx0=Bx0/qBunit
    By0=By0/qBunit
    Bz0=Bz0/qBunit
    Bxmax=maxval(dabs(Bx0(:,:)))
    Bymax=maxval(dabs(By0(:,:)))
    Bzmax=maxval(dabs(Bz0(:,:)))
    
    if(mype==0) then
      print*,'extrapolating 3D nonlinear force-free field from an observed Bx,By, and Bz'
      print*,'Updated magnetogram of',nx1,'by',nx2,'pixels.'
      print*,'Bxmax=',Bxmax
      print*,'Bymax=',Bymax
      print*,'Bzmax=',Bzmax
    endif
  end subroutine init_b_nlfff_data

  !==============================================================================
  ! Purpose: to provide initial condition provided by the potential field model.
  !==============================================================================
  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bf(ixI^S,1:ndir)
    logical, save :: first=.true.
    if(mype==0 .and. first) then
      write(*,*)'initializing grids with a potential field model...'
      first=.false.
    endif
    ! use Green's function method to extrapolate B
    Bz0 = Bz0/Bzmax   ! Note that mod_lfff needs to normalize Bz0 before calculation
    call calc_lin_fff(ixI^L,ixO^L,Bf,x,lalpha,llift)
    Bz0 = Bz0*Bzmax   ! un-normalize Bz0
    w(ixO^S,mag(:))=Bf(ixO^S,1:3)
    w(ixO^S,mom(:))=zero
    w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                  (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
  end subroutine initonegrid_usr
  
  !==============================================================================
  ! Purpose: to provide special boundary conditions set by users.
  !==============================================================================
  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    use mod_physics
    integer, intent(in) :: ixI^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: delxdely,delxdelz,delydelx,delydelz,delzdelx,delzdely
    double precision :: xlen^D,dxb^D,startpos^D,coeffrho,tmp1(ixG^T)
    integer :: ngl,ix^D,ixIM^L,ixbc^D,af

    ngl=nghostcells*(2**(refine_max_level-1))

    select case(iB)
    case(1)
      w(ixO^S,rho_)=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
      w(ixO^S,mom(1))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))
      w(ixO^S,mom(2))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
      w(ixO^S,mom(3))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
      ! 2nd order accuacy constant value extrapolation
      do ix1=ixOmax1,ixOmin1,-1
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                   (-w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
              +4.0d0*w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
      enddo
    case(2)
      w(ixO^S,rho_)=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
      w(ixO^S,mom(1))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))
      w(ixO^S,mom(2))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
      w(ixO^S,mom(3))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
      ! 2nd order accuacy constant value extrapolation
      do ix1=ixOmin1,ixOmax1
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                   (-w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
              +4.0d0*w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
      enddo
    case(3)
     w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
     w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(1))
     w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(2))
     w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(3))
     ! 2nd order accuacy constant value extrapolation
     do ix2=ixOmax2,ixOmin2,-1
       w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                  (-w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,mag(1):mag(3)) &
             +4.0d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,mag(1):mag(3)))
     enddo
    case(4)
      w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(1))
      w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(2))
      w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(3))
      ! 2nd order accuacy constant value extrapolation
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                    (-w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,mag(1):mag(3))&
               +4.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(1):mag(3)))
      enddo
    case(5)
      coeffrho=usr_grav*SRadius**2/Tiso
      w(ixO^S,rho_)=rhob*dexp(coeffrho*(1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))
      w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(2)) 
      w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(3)) 
      ! The first layer (with the resolution of the level 1 grid) are assigned with the fixed boundary data
      ! The outer ghost layers are assigned with the following one-sided difference scheme
      dxb^D=dx(^D,refine_max_level);
      af=nint(dxlevel(1)/dxb1)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        xlen^D=x({ix^D,},^D)-xprobmin^D+dble(ngl)*dxb^D;\
        if(af>=2) then
          call mpistop("Bottom boundary is not defined ")
        else
          ixbc^D=ceiling(xlen^D/dxb^D);\
          w(ix^D,mag(1))=Bx0(ixbc1,ixbc2)
          w(ix^D,mag(2))=By0(ixbc1,ixbc2)
          w(ix^D,mag(3))=Bz0(ixbc1,ixbc2)
        endif
      {end do\}
      ! af is set to be the number of ghost cells already with fixed boundary data
      af=1
      ! 4th order accuacy constant value extrapolation
      do ix3=ixOmax3-af,ixOmin3,-1
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/25.0d0)* &
             ( -3.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+4,mag(1):mag(3)) &
              +16.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mag(1):mag(3)) &
              -36.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(1):mag(3)) &
              +48.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(1):mag(3)))
      enddo
    case(6)
      ixIM^L=ixO^L;
      ixIMmin3=ixOmin3-1;ixIMmax3=ixOmax3;
      call getggrav(tmp1,ixI^L,ixIM^L,x)
      ! using the formula d(rho)/dz=g/T*rho
      do ix3=ixOmin3,ixOmax3
        w(ix3^%3ixO^S,rho_)=w(ix3-2^%3ixO^S,rho_)+(0.5d0*block%dx(ix3^%3ixO^S,3)&
                       +block%dx(ix3-1^%3ixO^S,3)+ 0.5d0*block%dx(ix3-2^%3ixO^S,3))&
                            *tmp1(ix3-1^%3ixO^S)/Tiso*w(ix3-1^%3ixO^S,rho_)
      enddo
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(1))
      w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(2))
      w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(3))
      ! 2nd order accuacy constant value extrapolation
      do ix3=ixOmin3,ixOmax3
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/3.0d0)* &
              (     -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-2,mag(1):mag(3)) &
               +4.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-1,mag(1):mag(3)))
      enddo
    case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,3)))**2
  end subroutine

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,3)=ggrid(ixO^S)
  end subroutine gravity

  !==============================================================================
  ! Purpose: Enforce additional refinement or coarsening. One can use the
  !          coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !==============================================================================
  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (any(x(ixO^S,3)<=xprobmin3+0.3d0)) then
      refine=1
      coarsen=-1
    endif
  end subroutine special_refine_grid

  !==============================================================================
  ! Purpose: 
  !   this subroutine can be used in convert, to add auxiliary variables to the
  !   converted output file, for further analysis using tecplot, paraview, ....
  !   these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  !   the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  !   corresponding normalization values (default value 1)
  !==============================================================================
  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters
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
    ! calculate Lorentz force
    qvec(ixO^S,1:ndir)=zero
    do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
      if(lvc(idir,jdir,kdir)/=0)then
        tmp(ixO^S)=curlvec(ixO^S,jdir)*w(ixO^S,mag(kdir))
        if(lvc(idir,jdir,kdir)==1)then
          qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
        else
          qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
        endif
      endif
    enddo; enddo; enddo
    do idir=1,ndir
      w(ixO^S,nw+3+ndir+idir)=qvec(ixO^S,idir)
    end do
    ! find magnetic dips
    !dip=0.d0
    !do idir=1,ndir
    !  call gradient(w(ixI^S,mag(3)),ixI^L,ixO^L,idir,tmp)
    !  dip(ixO^S)=dip(ixO^S)+w(ixO^S,b0_+idir)*tmp(ixO^S)
    !end do
    !where(dabs(w(ixO^S,mag(3)))<0.08d0 .and. dip(ixO^S)>=0.d0)
    !  w(ixO^S,nw+8)=1.d0
    !elsewhere
    !  w(ixO^S,nw+8)=0.d0
    !end where
  end subroutine specialvar_output

  !==============================================================================
  ! Purpose: names for special variable output
  !==============================================================================
  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Alfv divB beta j1 j2 j3 L1 L2 L3'
  end subroutine specialvarnames_output

end module mod_usr
