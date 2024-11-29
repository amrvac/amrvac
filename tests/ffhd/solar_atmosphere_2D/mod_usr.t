module mod_usr
  use mod_ffhd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav,trelax
  double precision :: heatunit,gzone,SRadius,dya
  double precision, allocatable :: pa(:),ra(:)
  integer, parameter :: jmax=50000

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_2D")

    unit_length        = 1.d9 !< cm
    unit_temperature   = 1.d6 !< K
    unit_numberdensity = 1.d9 !< cm^-3

    usr_add_aux_names   => specialvarnames_output 
    usr_aux_output      => specialvar_output
    usr_gravity         => gravity
    usr_init_one_grid   => initonegrid_usr
    usr_set_B0          => specialset_B0
    usr_set_parameters  => initglobaldata_usr
    usr_source          => special_source
    usr_special_bc      => specialbound_usr

    call ffhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    heatunit=unit_pressure/unit_time !< 3.697693390805347E-003 erg*cm^-3/s
    usr_grav=-2.74d4*unit_length/unit_velocity**2 !< solar gravity
    gzone=0.2d0 !< thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) !< cells size of high-resolution 1D solar atmosphere
    SRadius=69.61d0 !< Solar radius
    call inithdstatic
  end subroutine initglobaldata_usr

  subroutine inithdstatic
    !> hydrostatic vertical stratification of density, temperature, pressure
    integer :: i,j,na,ibc
    integer, parameter :: n_val=49
    double precision :: res
    double precision :: rpho,htra,Ttr,Fc,invT,kappa
    double precision :: h_val(n_val),t_val(n_val)
    double precision, allocatable :: ya(:),Ta(:),gg(:)

    rpho=0.71d15/unit_numberdensity !< bottom density at y=-gzone
    Fc=2.d5/heatunit/unit_length !< constant
    !> constant kappa but actually should be related to Ta, so not necessary to be the same with kappa in Thermal Conduction
    kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
    !> VAL-C table, only temperature
    data h_val / 0, 50, 100, 150, 250, &
                 350, 450, 515, 555, 605, &
                 655, 705, 755, 855, 905, &
                 980, 1065, 1180, 1280, 1380, &
                 1515, 1605, 1785, 1925, 1990, &
                 2016, 2050, 2070, 2080, 2090, &
                 2104, 2107, 2109, 2113, 2115, &
                 2120, 2129, 2160, 2200, 2230, &
                 2255, 2263, 2267, 2271, 2274, &
                 2280, 2290, 2298, 2543 /
    data t_val / 6420, 5840, 5455, 5180, 4780, &
                 4465, 4220, 4170, 4230, 4420, &
                 4730, 5030, 5280, 5650, 5755, &
                 5925, 6040, 6150, 6220, 6280, &
                 6370, 6440, 6630, 6940, 7160, &
                 7360, 7660, 7940, 8180, 8440, &
                 9500, 10700, 12300, 18500, 21000, &
                 22500, 23000, 23500, 24000, 24200, &
                 24500, 25500, 28000, 32000, 37000, &
                 50000, 89100, 141000, 447000 /
    h_val(1:n_val)=h_val(1:n_val)*1.d5/unit_length
    t_val(1:n_val)=t_val(1:n_val)/unit_temperature
    htra=maxval(h_val)
    Ttr=maxval(t_val)
    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))
    do j=1,jmax
      ya(j)=(dble(j)-0.5d0)*dya-gzone
      if(ya(j)>=htra) then
        Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
      else
        do i=1,n_val
          if(ya(j)<h_val(i+1)) then
            Ta(j)=t_val(i)+(ya(j)-h_val(i))*(t_val(i+1)-t_val(i))/(h_val(i+1)-h_val(i))
            exit
          endif
        enddo
      endif
      !> keep gg the same with the settings in gravity
      gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
    enddo
    ra(1)=rpho
    pa(1)=rpho*Ta(1)
    invT=0.d0
    do j=2,jmax
      invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
      pa(j)=pa(1)*dexp(invT*dya)
      ra(j)=pa(j)/Ta(j)
    end do
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    end do
    deallocate(ya,gg,Ta)
  end subroutine inithdstatic

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    !> initialize one grid
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: res
    integer :: ix^D,na

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      na=floor((x(ix^D,2)-xprobmin2+gzone)/dya+0.5d0)
      res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
      w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
      w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    {end do\}
    w(ixO^S,mom(:))=zero
    if(ffhd_hyperbolic_thermal_conduction) w(ixO^S,q_)=zero
    call ffhd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    !> special boundary types, user defined
    integer, intent(in) :: ixO^L,iB,ixI^L
    double precision, intent(in) :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    double precision :: Q(ixI^S),Qp(ixI^S)
    integer :: ix^D,ixOs^L,ixC^L,hxC^L,jxO^L

    select case(iB)
    case(3)
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(1))&
                   /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
        if(ffhd_hyperbolic_thermal_conduction) w(ix2^%2ixO^S,q_)=zero
      enddo
      call ffhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ixOs^L=ixO^L;
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmin2-1;
      call ffhd_get_pthermal(w,x,ixI^L,ixOs^L,pth)
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmax2;
      call getggrav(ggrid,ixI^L,ixOs^L,x)
      invT(ixOmin2-1^%2ixO^S)=w(ixOmin2-1^%2ixO^S,rho_)/pth(ixOmin2-1^%2ixO^S)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin2-1^%2ixO^S)=tmp(ixOmin2-1^%2ixO^S)+0.5d0*&
            (ggrid(ix2^%2ixO^S)+ggrid(ix2-1^%2ixO^S))*invT(ixOmin2-1^%2ixO^S)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)*dexp(tmp(ixOmin2-1^%2ixO^S)*dxlevel(2))
        w(ix2^%2ixO^S,rho_)=w(ix2^%2ixO^S,p_)*invT(ixOmin2-1^%2ixO^S)
        if(ffhd_hyperbolic_thermal_conduction) w(ix2^%2ixO^S,q_)=zero
      enddo
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,mom(1))&
                      /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      call ffhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,2)=ggrid(ixO^S)
  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,2)))**2
  end subroutine

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    integer, intent(in) :: ixI^L,ixO^L,iw^LIM
    double precision, intent(in) :: qdt,qtC,qt
    double precision, intent(in) :: x(ixI^S,1:ndim),wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    call getlQ(lQgrid,ixI^L,ixO^L,qt,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
  end subroutine special_source

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
  !> calculate background heating bQ
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qt,x(ixI^S,1:ndim),w(ixI^S,1:nw)
    double precision :: bQ0,lambda,bQgrid(ixI^S)

    bQ0=max((-0.04d0*qt+3.d0)*1.d-4,1.d-4)/heatunit
    lambda=5.d0
    bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,2)/lambda)
  end subroutine getbQ

  subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
  !> calculate localized heating lQ
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision :: lQgrid(ixI^S),lambdah,tr
    double precision :: sigmax,xr,xl,yh
    double precision :: tramp,lQ0

    lQ0=8.d-3/heatunit
    xr=4.15d0
    xl=-4.15d0
    yh=0.8d0
    sigmax=0.3d0
    lambdah=0.3d0
    tramp=1000.d0/unit_time
    trelax=100.d0
    if(qt .lt. trelax) then
      tr=zero
    elseif(qt .lt. trelax+tramp) then
      tr=(qt-trelax)/tramp
    else
      tr=one
    endif
    where(x(ixO^S,2) .lt. yh)
      lQgrid(ixO^S)=1.d0
    elsewhere
      lQgrid(ixO^S)=dexp(-(x(ixO^S,2)-yh)**2/lambdah**2)
    endwhere
    lQgrid(ixO^S)=lQ0*tr*lQgrid(ixO^S)*(dexp(-(x(ixO^S,1)-xr)**2/sigmax)+dexp(-(x(ixO^S,1)-xl)**2/sigmax))
  end subroutine getlQ

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_geometry
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    integer :: idims,ix^D
    double precision :: pth(ixI^S),wlocal(ixI^S,1:nw),qb(ixI^S,1:ndim),divqb(ixI^S),gradT(ixI^S)

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    !> output temperature
    call ffhd_get_pthermal(wlocal,x,ixI^L,ixI^L,pth)
    w(ixI^S,nw+1)=pth(ixI^S)/w(ixI^S,rho_)
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='Te'
  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)
    double precision :: B0,k1,k2,Btot(ixI^S)

    B0=Busr/(exp(-1.d0)-exp(-3.d0))
    k1=dpi/10.d0
    k2=3.d0*k1
    wB0(ixO^S,1)=+B0*dcos(k1*x(ixO^S,1))*dexp(-k1*x(ixO^S,2)) &
                 -B0*dcos(k2*x(ixO^S,1))*dexp(-k2*x(ixO^S,2))
    wB0(ixO^S,2)=-B0*dsin(k1*x(ixO^S,1))*dexp(-k1*x(ixO^S,2)) &
                 +B0*dsin(k2*x(ixO^S,1))*dexp(-k2*x(ixO^S,2))
    Btot(ixO^S)=dsqrt(wB0(ixO^S,1)**2+wB0(ixO^S,2)**2)
    wB0(ixO^S,1)=wB0(ixO^S,1)/Btot(ixO^S)
    wB0(ixO^S,2)=wB0(ixO^S,2)/Btot(ixO^S)
  end subroutine specialset_B0
end module mod_usr
