module mod_usr
  use mod_mhd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav
  double precision :: heatunit,gzone,B0,theta,SRadius,kx,ly,bQ0,dya
  double precision, allocatable :: pa(:),ra(:),ya(:)
  integer, parameter :: jmax=8000
  integer, parameter :: numL=1,npmax=1000

contains

  subroutine usr_init()
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
    usr_set_B0          => specialset_B0
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output 
    usr_init_vector_potential=>initvecpot_usr
    usr_special_convert => usrspecial_convert
    usr_set_field_w     => special_field_w
    usr_set_field       => special_userfield

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=1.d-4/heatunit ! background heating power density
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
    B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    theta=60.d0*dpi/180.d0 ! the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade 
    kx=dpi/(xprobmax1-xprobmin1)
    ly=kx*dcos(theta)
    SRadius=69.61d0 ! Solar radius
    ! hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic

  end subroutine initglobaldata_usr

  subroutine inithdstatic
  !! initialize the table in a vertical line through the global domain
    integer :: j,na,nb,ibc
    double precision, allocatable :: Ta(:),gg(:)
    double precision:: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa
  
    rpho=1.151d15/unit_numberdensity ! number density at the bottom relaxla
    Tpho=8.d3/unit_temperature ! temperature of chromosphere
    Ttop=1.5d6/unit_temperature ! estimated temperature in the top
    htra=0.2d0 ! height of initial transition region
    wtra=0.02d0 ! width of initial transition region 
    Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
    Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
    kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3

    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))
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
    nb=int(gzone/dya)
    ra(1)=rpho
    pa(1)=rpho*Tpho
    invT=gg(1)/Ta(1)
    invT=0.d0
    do j=2,jmax
       invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
       pa(j)=pa(1)*dexp(invT*dya)
       ra(j)=pa(j)/Ta(j)
    end do
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
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    if (idir==3) then
      A(ixC^S) = B0/ly*dcos(kx*xC(ixC^S,1))*dexp(-ly*xC(ixC^S,2))*dcos(theta)
    else
      A(ixC^S) = 0.d0
    end if

  end subroutine initvecpot_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    double precision :: Q(ixI^S),Qp(ixI^S)
    double precision :: A(ixGs^T,1:ndir) 
    double precision :: xC(ixGs^T,1:ndim),dxc(ixGs^T,1:ndim),circ(ixGs^T,1:ndim)
    double precision :: dxidir(ixGs^T)
    integer :: ix^D,ixOs^L,ixC^L,hxC^L,jxO^L,idir,idim,idim1,idim2

    select case(iB)
    case(3)
      !! fixed zero velocity
      do idir=1,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))&
                   /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end do
      !! fixed b1 b2 b3
      if(iprob==0 .or. B0field) then
        w(ixO^S,mag(:))=0.d0
      else if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix2=ixOsmax2,ixOsmin2,-1
             block%ws(ix2^%2ixOs^S,idir)=1.d0/3.d0*&
                   (-block%ws(ix2+2^%2ixOs^S,idir)&
               +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L-kr(2,^D);
        jxO^L=ixO^L+nghostcells*kr(2,^D);
        block%ws(ixOs^S,2)=zero
        call get_divb(w,ixI^L,jxO^L,Q)
        do ix2=ixOsmax2,ixOsmin2,-1
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix2^%2ixOs^S,2)=&
           -(Q(jxOmin2^%2jxO^S)*block%dvolume(jxOmin2^%2jxO^S)&
           -Qp(ix2+1^%2ixO^S)*block%dvolume(ix2+1^%2ixO^S))&
            /block%surfaceC(ix2^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
        w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
      else
        w(ixO^S,mag(1))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
        w(ixO^S,mag(2))= B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
        w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
      endif
      !! fixed gravity stratification of density and pressure pre-determined in initial condition
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
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
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,2)))**2
  end subroutine

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    ! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)

  end subroutine special_source

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
  ! calculate background heating bQ
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S)

    bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,2)/5.d0)

  end subroutine getbQ

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
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
    use mod_radiative_cooling
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),B2(ixI^S),tmp2(ixI^S),dRdT(ixI^S)
    double precision :: ens(ixI^S),divb(ixI^S),wlocal(ixI^S,1:nw)
    double precision :: Btotal(ixI^S,1:ndir),curlvec(ixI^S,1:ndir)
    integer :: idirmin,idir,ix^D

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! output temperature
    call mhd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)

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
    call get_normalized_divb(wlocal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=divb(ixO^S)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)
    ! output heating rate
    call getbQ(ens,ixI^L,ixO^L,global_time,wlocal,x)
    w(ixO^S,nw+5)=ens(ixO^S)
    ! store the cooling rate 
    if(mhd_radiative_cooling)call getvar_cooling(ixI^L,ixO^L,wlocal,x,ens,rc_fl)
    w(ixO^S,nw+6)=ens(ixO^S)

    ! store current
    call get_current(wlocal,ixI^L,ixO^L,idirmin,curlvec)
    do idir=1,ndir
      w(ixO^S,nw+6+idir)=curlvec(ixO^S,idir)
    end do
  
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te Alfv divB beta bQ rad j1 j2 j3'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixO^S,1)=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
    wB0(ixO^S,2)=+B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
    wB0(ixO^S,3)=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)

  end subroutine specialset_B0

  subroutine usrspecial_convert(qunitconvert)

    integer, intent(in) :: qunitconvert

    integer :: numP,nwP,nwL,nL
    double precision :: xs^L,lengthL,dL
    double precision :: dsep,distance
    integer :: nx_lcell^D
    double precision :: xseed(ndim)

    !---------------- example trace magnetic field lines -----------------!
    xsmin1=-2.95d0
    xsmax1=2.95d0
    xsmin2=1.d-3
    xsmax2=1.d-3

    dL=(xprobmax1-xprobmin1)/(domain_nx1*2**(refine_max_level-1))
    lengthL=10.d0
    numP=floor(lengthL/dL)+1
    dsep=0.1d0

    distance=dsqrt((xsmax1-xsmin1)**2+(xsmax2-xsmin2)**2)
    nL=floor(distance/dsep+0.5d0)+1
    nwP=1
    nwL=1

    call trace_Bfield_parallel(xs^L,dL,numP,nL,nwP,nwL)


    !------------ example calculate DEM for LOS x=0 with self difine field ------------!
    dL=(xprobmax2-xprobmin2)/(domain_nx2*2**(refine_max_level+1))
    lengthL=xprobmax2-xprobmin2
    numP=floor(lengthL/dL)+1
    nwP=2
    nwL=0
    xseed(1)=0.d0
    xseed(2)=dL

    call get_differential_emission_measure(xseed,dL,numP,nwP,nwL)

  end subroutine usrspecial_convert

  subroutine trace_Bfield_parallel(xs^L,dL,numP,nL,nwP,nwL)
    use mod_point_searching
    use mod_trace_field

    integer :: numP,nL,nwP,nwL
    double precision :: xs^L,dL

    double precision :: dsep^D
    double precision :: xpp(ndim),wpp(nw)
    double precision :: xfm(nL,numP,ndim),wPm(nL,numP,nwP),wLm(nL,nwL+1)
    integer :: nValid(nL)
    integer :: iL,j
    logical :: forwardm(nL)
    character(len=std_len) :: ftype,tcondi,fname
    
    dsep1=(xsmax1-xsmin1)/(nL-1)
    dsep2=(xsmax2-xsmin2)/(nL-1)
    do iL=1,nL
      xfm(iL,1,1)=xsmin1+(iL-1)*dsep1
      xfm(iL,1,2)=xsmin2+(iL-1)*dsep2
    enddo

    ftype='Bfield'
    tcondi='user1'
    do iL=1,nL 
      xpp(1:ndim)=xfm(iL,1,1:ndim)
      call get_point_w(xpp,wpp,'conserved')
      if (wpp(mag(2))>0.d0) then
        forwardm(iL)=.true.
      else
        forwardm(iL)=.false.
      endif
    enddo

    call trace_field_multi(xfm,wPm,wLm,dL,nL,numP,nwP,nwL,forwardm,ftype,tcondi)
    nValid(:)=int(wLm(:,1))

    ! output field line
    if (mype==0) then
      write(fname, '(a,i4.4,a)') 'Blines_',snapshotini,'.txt'
      open(1,file=fname)
      write(1,*) ' dL '
      write(1,*) dL
      write(1,*) ' nL '
      write(1,*) nL
      write(1,*) ' iL  nValid  length '
      do iL=1,nL
        write(1,*) iL,nValid(iL),wLm(iL,2)
      enddo
      write(1,*) ' iL  ip  x  y  rho '
      do iL=1,nL
        do j=1,nValid(iL)
          write(1,'(i8, i8, e15.7, e15.7, e15.7)') iL,j,xfm(iL,j,1),xfm(iL,j,2),wPm(iL,j,1)
        enddo
      enddo
      close(1)
    endif

  end subroutine trace_Bfield_parallel

  subroutine special_field_w(igrid,ip,xf,wP,wL,numP,nwP,nwL,dL,forward,ftype,tcondi)
    use mod_global_parameters
    use mod_point_searching

    integer, intent(in)                 :: igrid,ip,numP,nwP,nwL
    double precision, intent(in)        :: xf(numP,ndim)
    double precision, intent(inout)     :: wP(numP,nwP),wL(1+nwL)
    double precision, intent(in)        :: dL
    logical, intent(in)                 :: forward
    character(len=std_len), intent(in)  :: ftype,tcondi

    double precision :: xpp(1:ndim),wpp(1:nw)

    if (tcondi=='user1') then
      xpp(1:ndim)=xf(ip,1:ndim)
      call get_point_w_ingrid(igrid,xpp,wpp,'conserved')
      wP(ip,1)=wpp(rho_)
      if (ip==1) then
        wL(2)=0.d0
      else
        wL(2)=wL(2)+dL
      endif
    else if (tcondi=='user2') then
      xpp(1:ndim)=xf(ip,1:ndim)
      call get_point_w_ingrid(igrid,xpp,wpp,'primitive')
      wP(ip,1)=wpp(rho_)
      wP(ip,2)=wpp(p_)/wpp(rho_)
    else
      wP(ip,1)=0.d0
      wL(2)=0.d0
    endif

  end subroutine special_field_w

  subroutine get_differential_emission_measure(xseed,dL,numP,nwP,nwL)
    use mod_trace_field

    integer :: numP,nwP,nwL
    double precision :: dL
    double precision :: xseed(ndim)

    double precision :: xf(numP,ndim),wP(numP,nwP),wL(nwL+1)
    integer :: nValid,j,numTe,iTe
    logical :: forward
    character(len=std_len) :: ftype,tcondi,fname
    double precision :: Temin,Temax,dTelog,Telocal,EM
    double precision, allocatable :: Telog(:),DEM(:)

    ftype='ydir'
    tcondi='user2'
    forward=.true.
    xf(1,:)=xseed(:)
  
    call trace_field_single(xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi)
    nValid=int(wL(1))

    Temin=1.0d4
    Temax=1.6d6
    dTelog=0.1
    numTe=floor((log10(Temax)-log10(Temin))/dTelog)
    
    allocate(Telog(numTe),DEM(numTe))

    do iTe=1,numTe
      Telog(iTe)=log10(Temin)+dTelog*(iTe-1)
    enddo
    DEM=zero

    do j=1,nValid
      Telocal=wP(j,2)*unit_temperature
      iTe=floor((log10(Telocal)-log10(Temin))/dTelog)+1
      if (iTe>=1 .and. iTe<=numTe) then
        EM=(wP(j,1)*unit_numberdensity)**2
        DEM(iTe)=DEM(iTe)+EM*dL*unit_length
      endif
    enddo

    ! output field line
    if (mype==0) then
      write(fname, '(a,i4.4,a)') 'DEM_',snapshotini,'.txt'
      open(1,file=fname)
      write(1,*) ' numTe '
      write(1,*) numTe
      write(1,*) ' iL  Telog  DEM '
      do iTe=1,numTe
        write(1,*) iTe,Telog(iTe),DEM(iTe)
      enddo
      close(1)
    endif

    deallocate(Telog,DEM)

  end subroutine get_differential_emission_measure

  subroutine special_userfield(xfn,igrid,field,ftype)
    use mod_global_parameters

    integer,intent(in)                  :: igrid
    double precision, intent(in)        :: xfn(ndim)
    double precision, intent(inout)     :: field(ndim)
    character(len=std_len), intent(in)  :: ftype

    if (ftype=='ydir') then
      field(:)=zero
      field(2)=1.d0
    endif

  end subroutine special_userfield

end module mod_usr
