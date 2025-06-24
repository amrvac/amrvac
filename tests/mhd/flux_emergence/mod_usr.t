module mod_usr
  use mod_mhd
  use mod_thermal_emission
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav,vmax,La
  double precision :: heatunit,gzone,B0,theta,SRadius,kx,ly,bQ0,dya,B1,kx1,ly1,ita
  double precision, allocatable :: pa(:),ra(:)
  integer, parameter :: jmax=8000
  
  integer :: i,kk,vx
  double precision :: tt1,tt,mm,ll,nn,pp,trelax
  double precision :: lambda,c,lambda_min,lambda_max,rms
  double precision, allocatable :: va1(:),ta1(:),tb1(:),va2(:),ta2(:),tb2(:),v0(:),vm(:)
  double precision, allocatable :: vp(:,:),vn(:,:)
  integer, parameter :: vtim=300, vlim=6144, num = 1000, lth = 1, mnx = 6144
  !Use iprob in the par file to start relaxzation iprob<3 (No emergence), want to add emergence iprob>2 (for example iprob = 4) 
  !Define the start and end time of the emergence 
  double precision :: tstart=100.d0, tend=160.d0 

contains

  subroutine usr_init()
    use mod_global_parameters
    
    call set_coordinate_system("Cartesian_2.5D")

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    !usr_set_B0          => specialset_B0
    !usr_special_resistivity => special_eta
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output 
    usr_init_vector_potential=>initvecpot_usr
    !usr_set_wLR         => boundary_wLR
    usr_set_electric_field => boundary_electric_field

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_global_parameters
    heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    trelax=0.d0       !depends on the reset_time is .true. or not 
    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=3.d-4/heatunit ! background heating power density
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
    B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    theta=60.d0*dpi/180.d0 ! the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade 
    kx=dpi/(xprobmax1-xprobmin1)
    ly=kx*1.5*dcos(theta)
    SRadius=69.61d0 ! Solar radius
    vmax=2.d6/unit_velocity ! maximal driven velocity
    La=1.5d9/unit_length
    ! hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic

    B1=200./unit_magneticfield ! emerging magnetic field strength at the bottom
    ita=0.d0*dpi/180.d0  !the angle of emerging flux to the plane xy
    kx1=dpi*6.0/(xprobmax1-xprobmin1)
    ly1=kx1*dcos(ita)
    if (mype==0) then
      print*,'vmax',vmax
      print*,'kx',kx
      print*,'ly',ly
      print*,'kx1',kx1
      print*,'ly1',ly1
      print*,'unit_pressure',unit_pressure
    endif

  end subroutine initglobaldata_usr

  subroutine inithdstatic
    use mod_global_parameters
    use mod_solar_atmosphere
    ! initialize the table in a vertical line through the global domain
    integer :: j,na,ibc
    double precision, allocatable :: Ta(:),gg(:),ya(:)
    double precision :: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa
    double precision :: rhohc,hc
    logical :: simple_temperature_curve, alive
    character(len=100)           :: filename
  
    ! The following generate the random heating, get 'tb1.dat', 'tb2.dat', 'va1.dat', 'va2.dat' 
    ! left footpoint
    ! generate random T, typically 300s +/- 75s
    
    allocate(ta1(vtim))
    allocate(tb1(vtim))
    
    if (mype==0) then   
      tt=0.d0
      do i = 1, vtim
        call random_number(mm)  
        tt1 = 300.d0 + 75.d0 * 2. * (mm - 0.5)
        tt = tt + tt1
        tb1(i) = tt / unit_time
      end do

      write(filename,"(a)") "tb1.dat"
      inquire(file=filename,exist=alive)
      if(alive) then
        open(unit=21,file=filename,form='formatted',status='old',access='append')
      else
        open(unit=21,file=filename,form='formatted',status='new')
        write(21,'(a)') 'Left random T'
      endif
      write(21,'(3(es12.4))') tb1
      close(21)

      ta1 = tb1
      do i = 2, vtim
        ta1(i) = tb1(i)-tb1(i-1)
      end do
      !print*,'ta1',ta1
    endif
    
    call MPI_BARRIER(icomm,ierrmpi)
    if (npe>1) then
      call MPI_BCAST(tb1,vtim,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(ta1,vtim,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif
    
    ! spatial distribution
    
    vx=vtim*vlim
    allocate(va1(vx))
    
    if (mype==0) then  
    
      !generate vtim spectrums. The typical duration is 5 min, e.g., 500 min. Enough for the simulation.
    
      allocate(vp(num, mnx))
      allocate(v0(mnx))
    
      do kk = 1, vtim
        lambda_min=lth*1.d0/(mnx*1.d0)
        lambda_max=lth*2.d0

        do i=1,num
          lambda = i * 1.d0 / num * (lambda_max - lambda_min) + lambda_min
          ! E = k ** (-5/3) = lambda ** (5/3) = c**2;
          c=lambda**(5./6.)
          call random_number(nn)  
          do j = 1, mnx 
            vp(i,j) = c*sin(2.*dpi*(j*1.d0/(mnx)*lth)/lambda + 2.*dpi*2.*(nn-0.5))
          end do
        end do

        rms=0.d0
        v0=0.d0
        do j = 1, mnx  
          do i = 1, num
            v0(j) = v0(j) + vp(i, j) 
          end do
          rms=rms+v0(j)**2
        end do

        v0 = v0 / dsqrt(rms / mnx)     !normallize

        do j = 1, vlim
          va1(vlim*(kk-1)+j)=v0(j)**2
        end do

      end do

      write(filename,"(a)") "va1.dat"
      inquire(file=filename,exist=alive)
      if(alive) then
        open(unit=21,file=filename,form='formatted',status='old',access='append')
      else
        open(unit=21,file=filename,form='formatted',status='new')
        write(21,'(a)') 'Left random V'
      endif
      write(21,'(3(es12.4))') va1
      close(21)

    endif

    if (npe>1) then
      call MPI_BCAST(va1,vx,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif
    
    ! for right footpoint
    ! generate random T, typically 300s +/- 75s

    allocate(ta2(vtim))
    allocate(tb2(vtim))

    if (mype==0) then   
      tt=0.d0
      do i = 1, vtim
        call random_number(pp)   
        tt1 = 300.d0 + 75.d0 * 2. * (pp - 0.5)
        tt = tt + tt1
        tb2(i) = tt / unit_time
      end do
      
      write(filename,"(a)") "tb2.dat"
      inquire(file=filename,exist=alive)
      if(alive) then
        open(unit=21,file=filename,form='formatted',status='old',access='append')
      else
        open(unit=21,file=filename,form='formatted',status='new')
        write(21,'(a)') 'Right random T'
      endif
      write(21,'(3(es12.4))') tb2
      close(21)

      ta2 = tb2 
      do i = 2, vtim
        ta2(i) = tb2(i) - tb2(i-1)
      end do
    endif
    
    if (npe>1) then
      call MPI_BCAST(tb2,vtim,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(ta2,vtim,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif
    
    ! spatial distribution
    
    allocate(va2(vx))

    if (mype==0) then  
      allocate(vn(num, mnx))
      allocate(vm(mnx))

      do kk = 1, vtim
        lambda_min=lth*1.d0/(mnx*1.d0)
        lambda_max=lth*2.d0

        do i = 1, num
          lambda = i * 1.d0 / num * (lambda_max - lambda_min) + lambda_min
          c=lambda**(5./6.)
          call random_number(ll)         
          do j = 1, mnx  
            vn(i,j) = c*sin(2.*dpi*(j*1.d0/(mnx-1)*lth)/lambda+2.*dpi*2.*(ll-0.5))
          end do
        end do
     
        rms=0.d0
        vm=0.d0    
        do j = 1, mnx  
          do i = 1, num
          vm(j)=vm(j)+vn(i, j) 
          end do
          rms=rms+vm(j)**2
        end do          
        vm=vm/dsqrt(rms/mnx)    !normalization
         
        do j = 1, vlim
          va2(vlim*(kk-1)+j)=vm(j)**2
        end do
         
      end do

      write(filename,"(a)") "va2.dat"
      inquire(file=filename,exist=alive)
      if(alive) then
        open(unit=21,file=filename,form='formatted',status='old',access='append')
      else
        open(unit=21,file=filename,form='formatted',status='new')
        write(21,'(a)') 'Right random V'
      endif
      write(21,'(3(es12.4))') va2
      close(21)
    endif
  
    if (npe>1) then
      call MPI_BCAST(va2,vx,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif
  

!   ! For restart, just use the following to read the 'tb1.dat', 'tb2.dat', 'va1.dat', 'va2.dat', 
!   ! otherwise new files will be generated and takes time 
!    vx = vtim * vlim
!    allocate(va1(vx))
!    if(mype .eq. 0) then
!      open(unit=66,file='va1.dat',status='old')
!      do j=1,vx
!        read(66,*) va1(j)
!      enddo
!      close(66)
!    endif
!    if(npe>1) then
!      call MPI_BCAST(va1,vx,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
!    endif
!    allocate(va2(vx))
!    if(mype .eq. 0) then
!      open(unit=77,file='va2.dat',status='old')
!      do j=1,vx
!        read(77,*) va2(j)
!      enddo
!      close(77)
!    endif
!    if(npe>1) then
!      call MPI_BCAST(va2,vx,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
!    endif
!    allocate(ta1(vtim))
!    allocate(tb1(vtim))
!    if(mype .eq. 0) then
!      open(unit=88,file='tb1.dat',status='old')
!      do j=1,vtim
!        read(88,*) tb1(j)
!      enddo
!      close(88)
!    endif
!    if(npe>1) then
!      call MPI_BCAST(tb1,vtim,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
!    endif
!    allocate(ta2(vtim))
!    allocate(tb2(vtim))
!    if(mype .eq. 0) then
!      open(unit=99,file='tb2.dat',status='old')
!      do j=1,vtim
!        read(99,*) tb2(j)
!      enddo
!      close(99)
!    endif
!    if(npe>1) then
!      call MPI_BCAST(tb2,vtim,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
!    endif
!    ta1 = tb1
!    ta2 = tb2
!    do i = 2, vtim
!      ta1(i) = tb1(i) - tb1(i-1)
!      ta2(i) = tb2(i) - tb2(i-1)
!    enddo

    simple_temperature_curve=.true.
    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))

    if(simple_temperature_curve) then
      rpho=1.151d18/unit_numberdensity ! number density at the bottom of height table
      Tpho=6.d3/unit_temperature ! temperature of chromosphere
      Ttop=1.3d6/unit_temperature ! estimated temperature in the top
      htra=0.25d0 ! height of initial transition region
      wtra=0.03d0 ! width of initial transition region 
      Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
      Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
      kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      do j=1,jmax
         ya(j)=(dble(j)-0.5d0)*dya-gzone
         if(ya(j)>htra) then
           Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
         else
           Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra-0.027d0)/wtra)+1.d0)
         endif
         gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
      enddo
      !! solution of hydrostatic equation 
      ra(1)=rpho
      pa(1)=rpho*Tpho
      invT=gg(1)/Ta(1)
      invT=0.d0
      do j=2,jmax
         invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
         pa(j)=pa(1)*dexp(invT*dya)
         ra(j)=pa(j)/Ta(j)
      end do
    else
      do j=1,jmax
         ! get height table
         ya(j)=(dble(j)-0.5d0)*dya-gzone
         ! get gravity table
         gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
      enddo
      ! a coronal height at 10 Mm
      hc=1.d9/unit_length
      ! the number density at the coronal height
      rhohc=6.2d8/unit_numberdensity
      ! get density and pressure table in hydrostatic state with a preset temperature table
      call get_atm_para(ya,ra,pa,gg,jmax,'AL-C7',hc,rhohc)
    end if
    deallocate(ya,gg,Ta)
    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dya+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dya
    rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
    pb=pa(na)+res/dya*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

  end subroutine inithdstatic

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    ! initialize one grid
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: res
    integer :: ix^D,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 2.5D solar atmosphere'
      endif
      first=.false.
    endif
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
        na=floor((x(ix^D,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    {end do\}
    w(ixO^S,mom(:))=zero
    if(B0field) then
      w(ixO^S,mag(:))=zero
    else if(stagger_grid) then
      call b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block%ws,x)
      call mhd_face_to_center(ixO^L,block)
      w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
    else
      w(ixO^S,mag(1))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
      w(ixO^S,mag(2))= B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
      w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
    endif

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    use mod_global_parameters
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    if (idir==3) then
      A(ixC^S) = B0/ly*dcos(kx*xC(ixC^S,1))*dexp(-ly*xC(ixC^S,2))*dcos(theta)
    else
      A(ixC^S) = 0.d0
    end if

  end subroutine initvecpot_usr

  subroutine boundary_electric_field(ixI^L,ixO^L,qt,qdt,fE,s)
    use mod_global_parameters
    ! specify tangential electric field at physical boundaries 
    ! to fix or drive normal magnetic field
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    type(state)                        :: s
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    double precision :: xC(ixI^S,1:ndim), tmp(ixI^S), eB(ixI^S,1:ndir)
    integer :: idir,ixC^L,ixA^L
    double precision :: vy(ixI^S)

    if(s%is_physical_boundary(3)) then
      if ((iprob<3) .or. (qt>tend)) then
        ! to fix normal magnetic field at bottom boundary surface
        idir=3
        ixCmin^D=ixOmin^D+kr(idir,^D)-1;
        ixCmax^D=ixOmax^D;
        fE(ixCmin2^%2ixC^S,3)=0.d0
      else
        idir=3
        ixCmin^D=ixOmin^D+kr(idir,^D)-1;
        ixCmax^D=ixOmax^D;
        call driven_velocity(ixI^L,ixC^L,qt,s%x,vy)
        call emergeB(ixI^L,ixC^L,qt,s%x,eB)
        tmp(ixCmin2^%2ixC^S)=vy(ixCmin2^%2ixC^S)*eB(ixCmin2^%2ixC^S,1)*qdt*s%dsC(ixCmin2^%2ixC^S,3)
        ! neighor cell center average to get cell edge
        ixAmin^D=ixCmin^D+kr(1,^D);
        ixAmax^D=ixCmax^D+kr(1,^D);
        fE(ixCmin2^%2ixC^S,3)=0.5d0*(tmp(ixCmin2^%2ixC^S)+tmp(ixAmin2^%2ixA^S))
      end if
    end if

  end subroutine boundary_electric_field

  !> allow user to specify variables' left and right state at physical boundaries to control flux through the boundary surface' 
  subroutine boundary_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    double precision, intent(inout) :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw)
    type(state)                     :: s

    double precision :: vy(ixI^S)

    if ((iprob<3) .or. (qt>tend)) then
      if(s%is_physical_boundary(3).and.idir==2) then
        wLp(ixOmin2^%2ixO^S,mom(:))=0.d0
        wRp(ixOmin2^%2ixO^S,mom(:))=0.d0
        wLC(ixOmin2^%2ixO^S,mom(:))=0.d0
        wRC(ixOmin2^%2ixO^S,mom(:))=0.d0
      end if
    else
      if(s%is_physical_boundary(3).and.idir==2) then
        call driven_velocity(ixI^L,ixO^L,qt,s%x,vy)
        wLp(ixOmin2^%2ixO^S,mom(1))=0.d0
        wRp(ixOmin2^%2ixO^S,mom(1))=0.d0
        wLC(ixOmin2^%2ixO^S,mom(1))=0.d0
        wRC(ixOmin2^%2ixO^S,mom(1))=0.d0
        wLp(ixOmin2^%2ixO^S,mom(2))=vy(ixOmin2^%2ixO^S)
        wRp(ixOmin2^%2ixO^S,mom(2))=wRp(ixOmin2^%2ixO^S,mom(2))
        wLC(ixOmin2^%2ixO^S,mom(2))=wLp(ixOmin2^%2ixO^S,mom(2))*wLp(ixOmin2^%2ixO^S,rho_)
        wRC(ixOmin2^%2ixO^S,mom(2))=wRp(ixOmin2^%2ixO^S,mom(2))*wLp(ixOmin2^%2ixO^S,rho_)
        wLp(ixOmin2^%2ixO^S,mom(3))=0.d0
        wRp(ixOmin2^%2ixO^S,mom(3))=0.d0
        wLC(ixOmin2^%2ixO^S,mom(3))=0.d0
        wRC(ixOmin2^%2ixO^S,mom(3))=0.d0
      end if
    end if

  end subroutine boundary_wLR


  subroutine driven_velocity(ixI^L,ixO^L,qt,xx,vdr)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in) :: qt, xx(ixI^S,1:ndim)
    double precision, intent(out) :: vdr(ixI^S)

    double precision :: tramp1,ft

    if ((qt>=tstart) .and. (qt<=tend)) then
      where(xx(ixO^S,1)<-2.0 .and. xx(ixO^S,1)>-4.0)
         vdr(ixO^S)=0.1*vmax
      else where
        vdr(ixO^S)=0.d0
      end where
    else 
      vdr(ixO^S)=0.d0
    end if

  end subroutine driven_velocity


  subroutine emergeB(ixI^L,ixO^L,qt,xx,wB)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt, xx(ixI^S,1:ndim)
    double precision, intent(inout) :: wB(ixI^S,1:ndir)
    double precision ::x_0

    x_0=-3.0d0

    if ((qt>=tstart) .and. (qt<=tend)) then
      wB(ixO^S,1)=B1*dcos(kx1*(xx(ixO^S,1)-x_0))*dexp(-ly1*(xx(ixO^S,2)-(-1.0d0+(qt-tstart)*0.1*vmax)))*dcos(ita)  
      wB(ixO^S,2)=-B1*dsin(kx1*(xx(ixO^S,1)-x_0))*dexp(-ly1*(xx(ixO^S,2)-(-1.0d0+(qt-tstart)*0.1*vmax)))          
      wB(ixO^S,3)=B1*dcos(kx1*(xx(ixO^S,1)-x_0))*dexp(-ly1*(xx(ixO^S,2)-(-1.0d0+(qt-tstart)*0.1*vmax)))*dsin(ita)  
    else
      wB(ixO^S,1)=B1*dcos(kx1*(xx(ixO^S,1)-x_0))*dexp(-ly1*(xx(ixO^S,2)-(-1.0d0+(tend-tstart)*0.1*vmax)))*dcos(ita)  
      wB(ixO^S,2)=-B1*dsin(kx1*(xx(ixO^S,1)-x_0))*dexp(-ly1*(xx(ixO^S,2)-(-1.0d0+(tend-tstart)*0.1*vmax)))          
      wB(ixO^S,3)=B1*dcos(kx1*(xx(ixO^S,1)-x_0))*dexp(-ly1*(xx(ixO^S,2)-(-1.0d0+(tend-tstart)*0.1*vmax)))*dsin(ita)  
    endif

  end subroutine emergeB


  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    double precision :: Q(ixI^S),Qp(ixI^S)
    integer :: ix^D,ixOs^L,ixC^L,hxC^L,jxO^L,idir
    double precision ::x_0

    x_0=-3.0d0

    select case(iB)
    case(3)
      if ((iprob<3) .or. (qt>tend)) then
        !! fixed zero velocity
        do idir=1,ndir
          w(ixO^S,mom(idir))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))&
                     /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
        end do
      else
        call driven_velocity(ixI^L,ixO^L,qt,x,w(ixI^S,mom(2)))
        w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(1))&
                    /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
        w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(3))&
                    /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end if
      !! fixed b1 b2 b3
      if(iprob==0 .or. B0field) then
        w(ixO^S,mag(:))=0.d0
      else if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          ! 4th order one-sided equal-gradient extrapolation
          do ix2=ixOsmax2,ixOsmin2,-1
            block%ws(ix2^%2ixOs^S,idir)= &
              0.12d0*block%ws(ix2+5^%2ixOs^S,idir) &
             -0.76d0*block%ws(ix2+4^%2ixOs^S,idir) &
             +2.08d0*block%ws(ix2+3^%2ixOs^S,idir) &
             -3.36d0*block%ws(ix2+2^%2ixOs^S,idir) &
             +2.92d0*block%ws(ix2+1^%2ixOs^S,idir)
          end do
        end do
        ixOs^L=ixO^L-kr(2,^D);
        jxO^L=ixO^L+nghostcells*kr(2,^D);
        block%ws(ixOs^S,2)=zero
        do ix2=ixOsmax2,ixOsmin2,-1
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix2^%2ixOs^S,2)=Qp(ix2+1^%2ixO^S)*block%dvolume(ix2+1^%2ixO^S)&
            /block%surfaceC(ix2^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
        if(iprob>2) then
          if ((qt>=tstart) .and. (qt<=tend)) then
            w(ixO^S,mag(3))=B1*dcos(kx1*(x(ixO^S,1)-x_0))*dexp(-ly1*(x(ixO^S,2)-(-1.0d0+(qt-tstart)*0.1*vmax)))*dsin(ita)
          else
            w(ixO^S,mag(3))=B1*dcos(kx1*(x(ixO^S,1)-x_0))*dexp(-ly1*(x(ixO^S,2)-(-1.0d0+(tend-tstart)*0.1*vmax)))*dsin(ita)
          endif
        else
          do ix2=ixOmax2,ixOmin2,-1
            w(ix2^%2ixO^S,mag(3))= &
              0.12d0*w(ix2+5^%2ixO^S,mag(3)) &
             -0.76d0*w(ix2+4^%2ixO^S,mag(3)) &
             +2.08d0*w(ix2+3^%2ixO^S,mag(3)) &
             -3.36d0*w(ix2+2^%2ixO^S,mag(3)) &
             +2.92d0*w(ix2+1^%2ixO^S,mag(3))
          end do
        end if
      else
        w(ixO^S,mag(1))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
        w(ixO^S,mag(2))= B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
        w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
      endif
      !! fixed gravity stratification of density and pressure pre-determined in initial condition
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
        !Pgas=Pbot-B^2/8*pi  Pbot=2.44*10^9 dyn cm^{-2}
        !w(ixOmin1:ixOmax1,ix2,p_)=2.44*10**9/0.31758 -(w(ixOmin1:ixOmax1,ix2,mag(1))**2+w(ixOmin1:ixOmax1,ix2,mag(2))**2+w(ixOmin1:ixOmax1,ix2,mag(3))**2)/2
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ixOs^L=ixO^L;
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmin2-1;
      call mhd_get_pthermal(w,x,ixI^L,ixOs^L,pth)
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmax2;
      call getggrav(ggrid,ixI^L,ixOs^L,x)
      !> fill pth, rho ghost layers according to gravity stratification
      invT(ixOmin2-1^%2ixO^S)=w(ixOmin2-1^%2ixO^S,rho_)/pth(ixOmin2-1^%2ixO^S)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin2-1^%2ixO^S)=tmp(ixOmin2-1^%2ixO^S)+0.5d0*&
            (ggrid(ix2^%2ixO^S)+ggrid(ix2-1^%2ixO^S))*invT(ixOmin2-1^%2ixO^S)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)*dexp(tmp(ixOmin2-1^%2ixO^S)*dxlevel(2))
        w(ix2^%2ixO^S,rho_)=w(ix2^%2ixO^S,p_)*invT(ixOmin2-1^%2ixO^S)
      enddo
      !> fixed zero velocity
      do idir=1,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))&
                     /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      end do
      !> zero normal gradient extrapolation
      if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix2=ixOsmin2,ixOsmax2
             block%ws(ix2^%2ixOs^S,idir)=1.d0/3.d0*&
                   (-block%ws(ix2-2^%2ixOs^S,idir)&
               +4.d0*block%ws(ix2-1^%2ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L;
        jxO^L=ixO^L-nghostcells*kr(2,^D);
        block%ws(ixOs^S,2)=zero
        call get_divb(w,ixI^L,jxO^L,Q)
        do ix2=ixOsmin2,ixOsmax2
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix2^%2ixOs^S,2)=&
            (Q(jxOmax2^%2jxO^S)*block%dvolume(jxOmax2^%2jxO^S)&
           -Qp(ix2^%2ixO^S)*block%dvolume(ix2^%2ixO^S))&
            /block%surfaceC(ix2^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
        do ix2=ixOmin2,ixOmax2
          w(ixOmin1:ixOmax1,ix2,mag(3))=(1.0d0/3.0d0)* &
                      (-w(ixOmin1:ixOmax1,ix2-2,mag(3))&
                 +4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(3)))
        enddo
      else
        do ix2=ixOmin2,ixOmax2
          w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* &
                      (-w(ixOmin1:ixOmax1,ix2-2,mag(:))&
                 +4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
        enddo
      end if
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select
    
  end subroutine specialbound_usr

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,2)=ggrid(ixO^S)

  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,2)))**2
  end subroutine

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    ! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    call getlQ(lQgrid,ixI^L,ixO^L,qt,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)

  end subroutine special_source

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
  ! calculate background heating bQ
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S)

    bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,2)/5.d0)

  end subroutine getbQ

  
  subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    integer :: Bp,Bt1,Bt2,ix^D,i
    double precision :: lQgrid(ixI^S),lQ0,lQd,hp,A(ixI^S),Atop,Alow
    double precision :: tl1,tl2,tl3,tl4,tl5,tl6,t2,stt
    double precision :: tr1,tr2,tr3,tr4,tr5,tr6
    double precision :: va0
    integer :: mode
    
    lQ0=3.d-2 / heatunit
    va0=0.4d0
    stt = (xprobmax1- xprobmin1) / (dble(vlim) - 1.d0)
    t2 = qt - trelax
    lQd=0.45d0
    hp=0.2d0
    mode=1

    select case(mode)
    ! mode 1 for gaussian, mode 2 for sinc
    case(1)
      if(qt .lt. trelax) then
        lQgrid(ixO^S)=0.d0
      else
        do i = 4, vtim
          if(qt .lt. trelax + tb1(i)) then
            tl1 = dexp(-(t2 - tb1(i-3))**2/(0.5d0*ta1(i-3))**2)
            tl2 = dexp(-(t2 - tb1(i-2))**2/(0.5d0*ta1(i-2))**2)
            tl3 = dexp(-(t2 - tb1(i-1))**2/(0.5d0*ta1(i-1))**2)
            tl4 = dexp(-(t2 - tb1(i))**2/(0.5d0*ta1(i))**2)
            tl5 = dexp(-(t2 - tb1(i+1))**2/(0.5d0*ta1(i+1))**2)
            tl6 = dexp(-(t2 - tb1(i+2))**2/(0.5d0*ta1(i+2))**2)
            do ix2 = ixOmin2,ixOmax2
            do ix1 = ixOmin1,ixOmax1
              Bp = floor((x(ix1,ix2,1) - xprobmin1) / stt) + 1
              if (Bp .lt. 1)  Bp=1
              if(x(ix1,ix2,1) .lt. 0.d0) then
                lQgrid(ix1,ix2) = tl1 * (va1(Bp + (i-4)*vlim)*(one-va0)+va0) &
                                + tl2 * (va1(Bp + (i-3)*vlim)*(one-va0)+va0) &
                                + tl3 * (va1(Bp + (i-2)*vlim)*(one-va0)+va0) &
                                + tl4 * (va1(Bp + (i-1)*vlim)*(one-va0)+va0) &
                                + tl5 * (va1(Bp + (i)*vlim)*(one-va0)+va0) &
                                + tl6 * (va1(Bp + (i+1)*vlim)*(one-va0)+va0)
              endif
            enddo
            enddo
            exit
          endif
        enddo

        do i = 4, vtim
          if(qt .lt. trelax + tb2(i)) then
            tr1 = dexp(-(t2 - tb2(i-3))**2/(0.5d0*ta2(i-3))**2)
            tr2 = dexp(-(t2 - tb2(i-2))**2/(0.5d0*ta2(i-2))**2)
            tr3 = dexp(-(t2 - tb2(i-1))**2/(0.5d0*ta2(i-1))**2)
            tr4 = dexp(-(t2 - tb2(i))**2/(0.5d0*ta2(i))**2)
            tr5 = dexp(-(t2 - tb2(i+1))**2/(0.5d0*ta2(i+1))**2)
            tr6 = dexp(-(t2 - tb2(i+2))**2/(0.5d0*ta2(i+2))**2)
            do ix2 = ixOmin2,ixOmax2
            do ix1 = ixOmin1,ixOmax1
              Bp = floor((x(ix1,ix2,1) - xprobmin1) / stt) + 1   
              if (Bp .lt. 1)  Bp=1
              if(x(ix1,ix2,1) .ge. 0.d0) then
                lQgrid(ix1,ix2) = tr1 * (va2(Bp + (i-4)*vlim)*(one-va0)+va0) &
                                + tr2 * (va2(Bp + (i-3)*vlim)*(one-va0)+va0) &
                                + tr3 * (va2(Bp + (i-2)*vlim)*(one-va0)+va0) &
                                + tr4 * (va2(Bp + (i-1)*vlim)*(one-va0)+va0) &
                                + tr5 * (va2(Bp + (i)*vlim)*(one-va0)+va0) &
                                + tr6 * (va2(Bp + (i+1)*vlim)*(one-va0)+va0)
              endif
            enddo
            enddo
            exit
          endif
        enddo
      endif
    case(2)
      if(qt .lt. trelax) then
        tr1 = 0.d0
        tr3 = 0.d0
        Bt1 = 0
      else if(qt .lt. trelax + tb1(1)) then
        tr1 = dsin(t2 * dpi / tb1(2))
        tr3 = 0.d0
        Bt1 = 0
      else if(qt .lt. trelax + tb1(2)) then
        tr1 = dsin(t2 * dpi / tb1(2))
        tr3 = dsin((t2 - tb1(1)) * dpi / (tb1(3) - tb1(1)))
        Bt1 = 0
      else
        do i = 2, vtim
          if(t2 .lt. tb1(i)) then
            Bt1 = i - 1
            exit
          endif
        enddo
        tr1 = dsin((t2 - tb1(Bt1-1)) * dpi / (tb1(Bt1+1) - tb1(Bt1-1)))
        tr3 = dsin((t2 - tb1(Bt1)) * dpi / (tb1(Bt1+2) - tb1(Bt1)))
      endif
      if(qt .lt. trelax) then
        tr2 = 0.d0
        tr4 = 0.d0
        Bt2 = 0
      else if(qt .lt. trelax + tb2(1)) then
        tr2 = dsin(t2 * dpi / tb2(2))
        tr4 = 0.d0
        Bt2 = 0
      else if(qt .lt. trelax + tb2(2)) then
        tr1 = dsin(t2 * dpi / tb2(2))
        tr3 = dsin((t2 - tb2(1)) * dpi / (tb2(3) - tb2(1)))
        Bt2 = 0
      else
        do i = 2, vtim
          if(t2 .lt. tb2(i)) then
            Bt2 = i - 1
            exit
          endif
        enddo
        tr2 = dsin((t2 - tb2(Bt2-1)) * dpi / (tb2(Bt2+1) - tb2(Bt2-1)))
        tr4 = dsin((t2 - tb2(Bt2)) * dpi / (tb2(Bt2+2) - tb2(Bt2)))
      endif
      do ix1 = ixOmin1,ixOmax1
        do ix2 = ixOmin2,ixOmax2
          Bp = floor((x(ix1,ix2,1) - xprobmin1) / stt) + 1
          if(x(ix1,ix2,1) .lt. 0.d0) then
            lQgrid(ix1,ix2) = tr1 * (va1(Bp + Bt1*vlim)*(one-va0)+va0) + tr3 * (va1(Bp+(Bt1+1)*vlim)*(one-va0)+va0)
          else
            lQgrid(ix1,ix2) = tr2 * (va2(Bp + Bt2*vlim)*(one-va0)+va0) + tr4 * (va2(Bp+(Bt2+1)*vlim)*(one-va0)+va0)
          endif
        enddo
      enddo
    case default
      call mpistop('unknow mode')
    end select

    call initvecpot_usr(ixI^L, ixO^L, x, A, 3)
    Atop=B0/kx*cos(kx*0.d0)
    Alow=B0/kx*cos(kx*6.0d0)
    where(A(ixO^S)<Atop .and. A(ixO^S)>Alow)
      lQgrid(ixO^S)=lQ0*lQgrid(ixO^S)
    end where
    where(x(ixO^S,2)>hp)
      lQgrid(ixO^S)=lQgrid(ixO^S)*dexp(-(x(ixO^S,2)-hp)**2/lQd)
    end where

    !where(x(ixO^S,1)<0.d0)
    !  lQgrid(ixO^S)=0.d0
    !end where

  end subroutine getlQ


  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (any(x(ixO^S,2)<=xprobmin2+0.05d0)) then
      refine=1
      coarsen=-1
    endif

  end subroutine special_refine_grid


  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),B2(ixI^S),tmp2(ixI^S),dRdT(ixI^S)
    double precision :: ens(ixI^S),divb(ixI^S),wlocal(ixI^S,1:nw)
    double precision :: Btotal(ixI^S,1:ndir),curlvec(ixI^S,1:ndir)
    double precision :: Te(ixI^S)
    integer :: idirmin,idir,ix^D

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! output temperature
    call mhd_get_pthermal(wlocal,x,ixI^L,ixI^L,pth)
    Te(ixI^S)=pth(ixI^S)/w(ixI^S,rho_)
    w(ixO^S,nw+1)=Te(ixO^S)

    do idir=1,ndir
      if(B0field) then
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0)
      else
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))
      endif
    end do
    ! B^2
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)

    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+2)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))

    ! output divB1
    call get_divb(wlocal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=divb(ixO^S)

    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)

    ! store current
    call get_current(wlocal,ixI^L,ixO^L,idirmin,curlvec)
    do idir=1,ndir
      w(ixO^S,nw+4+idir)=curlvec(ixO^S,idir)
    end do

  end subroutine specialvar_output


  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te Alfv divB beta j1 j2 j3' 

  end subroutine specialvarnames_output

  subroutine special_eta(w,ixI^L,ixO^L,idirmin,x,current,eta)
    ! Set the common "eta" array for resistive MHD based on w or the
    ! "current" variable which has components between idirmin and 3.
    integer, intent(in) :: ixI^L, ixO^L, idirmin
    double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision :: current(ixI^S,7-2*ndir:3), eta(ixI^S)
    double precision :: eta2,etam,rad(ixI^S)

    eta2=1.d-2
    etam=1.d-1
    if (global_time>tend) then
      !rad(ixO^S)=dsqrt(sum(current(ixO^S,:)**2,dim=ndim+1))/w(ixO^S,rho_)/q_e
      where(abs(current(ixO^S,3))>20.d0 .and. x(ixO^S,2)>0.6d0)
        eta(ixO^S)=eta2
      elsewhere
        eta(ixO^S)=0.d0
      endwhere
      where(eta(ixO^S)>etam)
        eta(ixO^S)=etam
      endwhere
    else
      eta(ixO^S)=0.d0
    end if 

  end subroutine special_eta
  

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
    use mod_global_parameters
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixO^S,1)=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
    wB0(ixO^S,2)=+B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
    wB0(ixO^S,3)=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)

  end subroutine specialset_B0



end module mod_usr
