module mod_usr
  use mod_hd
  implicit none
  double precision :: k_B,mass_H
  double precision :: heatunit,ds,gzone,usr_grav,kai0
  double precision, allocatable :: pbc(:,:),rbc(:,:)
  double precision, allocatable :: s_o(:),g_o(:),T_o(:),p_o(:),r_o(:)
  double precision :: dxlv1,trelax

contains

  subroutine usr_init()
    unit_length        = 1.d9 
    unit_temperature   = 1.d6 
    unit_numberdensity = 1.d9 
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output 
    usr_refine_grid     => special_refine_grid
    call set_coordinate_system("Cartesian_1D")
    call hd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_thermal_conduction
    k_B=1.3806d-16          ! erg*K^-1
    mass_H=1.67262d-24      ! g
    usr_grav=-2.74d4*unit_length/unit_velocity**2
    kai0=0.21633973559351147d0 ! kappa in heat conduction
    heatunit=unit_pressure/unit_time ! 3.697693390805347E-003 erg*cm^-3/s
    if(mype .eq. 0) then
      print*, 'unit_density = ', unit_density
      print*, 'unit_pressure = ', unit_pressure
      print*, 'unit_velocity = ', unit_velocity
      print*, 'unit_time = ', unit_time
      print*, 'usr_grav = ', usr_grav
      print*, 'kpara', tc_k_para
    end if
    !> length of the ghostzones
    gzone=0.005d0*(xprobmax1-xprobmin1)
    !> lv1 dx
    dxlv1=(xprobmax1-xprobmin1)/domain_nx1
    !> time for relaxation
    trelax=100.d0
    call inithdstatic
  end subroutine initglobaldata_usr

  !> calculate geometry (helix) & initial hydrodynamic state, see Zhou et al. 2017 for details
  subroutine inithdstatic
    integer :: i,j,ilv,ibc,nb
    integer, parameter :: imax=10001,jmax=5001
    double precision :: s1,ra,s2,dip,span,loop,Fc,sj
    double precision :: Tcor,Tpho,htra,wtra,p0,res,rhob,pb
    double precision, dimension(jmax) :: theta,x_i,y_i,z_i,s_i
    double precision, dimension(imax) :: x_o,y_o,z_o
    double precision, dimension(imax) :: dz,sumt

    s1=0.8d0 ! foot
    ra=0.5d0 ! shoulder
    s2=s1+ra*dpi/2.d0
    dip=0.3d0
    span=8.d0
    Tcor=1.0d0
    Tpho=6.d-3
    htra=0.2719d0
    wtra=0.025d0
    allocate(s_o(imax))
    allocate(g_o(imax))
    allocate(t_o(imax))
    allocate(r_o(imax))
    allocate(p_o(imax))
    do j=1,jmax
      theta(j)=2.d0*dpi*(j-1)/(jmax-1)
    end do
    s_i=zero
    do j=1,jmax
      x_i(j)=(theta(j)/dpi/2.d0)*span
      y_i(j)=dip/2.d0*sin(theta(j))
      z_i(j)=dip/2.d0*cos(theta(j))
      if(j .gt. 1) then
        s_i(j)=s_i(j-1)+sqrt((x_i(j)-x_i(j-1))**2+(y_i(j)-y_i(j-1))**2+(z_i(j)-z_i(j-1))**2)
      end if
    end do
    loop=s_i(jmax)+s2*2
    if(mype .eq. 0) then
      print*, 'total loop length is', loop
      print*, 'xprobmax1-xprobmin1 is', xprobmax1-xprobmin1
      if(abs(loop-(xprobmax1-xprobmin1)) .gt. dxlv1**refine_max_level) then
        call mpistop("please change xprobmax1 = loop")
      end if
    end if
    gzone=0.005d0*loop
    do i=1,imax
      s_o(i)=loop*1.01d0*(i-1)/(imax-1)-gzone
    end do 
    ds=loop/(imax-1)*1.01d0
    !> xz plane
    y_o=zero
    do i=1,imax
      if(s_o(i) .lt. s1) then
        x_o(i)=-span/2.d0-ra
        z_o(i)=s_o(i)
      else if(s_o(i) .lt. s2) then
        x_o(i)=-span/2.d0-ra*cos((s_o(i)-s1)/ra)
        z_o(i)=s1+ra*sin((s_o(i)-s1)/ra)
      else if(s_o(i) .lt. loop-s2) then
        do j=1,jmax-1
          sj=s_o(i)-s2
          if(sj .ge. s_i(j) .and. sj .lt. s_i(j+1)) then
            x_o(i)=x_i(j)+(x_i(j+1)-x_i(j))*(sj-s_i(j))/(s_i(j+1)-s_i(j))-span/2.d0  
            y_o(i)=y_i(j)+(y_i(j+1)-y_i(j))*(sj-s_i(j))/(s_i(j+1)-s_i(j))
            z_o(i)=z_i(j)+(z_i(j+1)-z_i(j))*(sj-s_i(j))/(s_i(j+1)-s_i(j))+s1+ra-dip/2.d0
          end if
        end do
      else if(s_o(i) .lt. loop-s1) then
        x_o(i)=span/2.d0+ra*cos((loop-s1-s_o(i))/ra)
        z_o(i)=s1+ra*sin((loop-s1-s_o(i))/ra)
      else
        x_o(i)=span/2.d0+ra
        z_o(i)=loop-s_o(i)
      end if
    end do
    !> grav 
    do i=1,imax
      if(s_o(i) .lt. s1) then
        g_o(i)=1.d0
      else if(s_o(i) .lt. s2) then
        g_o(i)=cos(dpi/2.d0*(s_o(i)-s1)/(s2-s1))
      else if(s_o(i) .lt. loop-s2) then
        g_o(i)=-(z_o(i)-z_o(i+1))/sqrt((x_o(i+1)-x_o(i))**2+(y_o(i+1)-y_o(i))**2+(z_o(i+1)-z_o(i))**2)
      else if(s_o(i) .lt. loop-s1) then
        g_o(i)=-cos(dpi/2.d0*(loop-s_o(i)-s1)/(s2-s1))
      else
        g_o(i)=-1.d0
      end if
    end do
    !> temperature & density
    do i=1,imax
!      if(z_o(i) .lt. htra) then
        t_o(i)=Tpho+0.5d0*(Tcor-Tpho)*(1.d0+tanh((z_o(i)-htra)/wtra))
!      else
!        t_o(i)=(7.d0*Fc*(z_o(i)-htra)/2.d0/kai0+Ttra**3.5d0)**(2.d0/7.d0)
!      end if
    end do
    !> pressure
    p0=420.d0
    nb=floor(gzone/ds+0.5d0)
    do i=2,imax-1
      dz(i)=(z_o(i+1)-z_o(i-1))*0.5d0
    end do
    dz(1)=dz(2)
    dz(imax)=dz(imax-1)
    sumt=zero
    do i=1, nb
      do j=i, nb
        sumt(i)=sumt(i)+1.d0/t_o(j)*dz(j)
      end do
      p_o(i)=p0*exp(-usr_grav*sumt(i))
      r_o(i)=p_o(i)/t_o(i)
    end do
    do i=nb+1, (imax+1)/2
      do j=nb+1, i
        sumt(i)=sumt(i)+1.d0/t_o(j)*dz(j)
      end do
      p_o(i)=p0*exp(usr_grav*sumt(i))
      r_o(i)=p_o(i)/t_o(i)
    end do
    do i=(imax+1)/2+1, imax
      p_o(i)=p_o(imax+1-i)
      r_o(i)=r_o(imax+1-i)
    end do
    g_o=g_o*usr_grav
    !> calculate fixed bound
    nb=floor(gzone/ds+0.5d0)
    res=gzone-(dble(nb)-0.5d0)*ds
    rhob=r_o(nb)+res/ds*(r_o(nb+1)-r_o(nb))
    pb=p_o(nb)+res/ds*(p_o(nb+1)-p_o(nb))
    allocate(rbc(nghostcells,refine_max_level))
    allocate(pbc(nghostcells,refine_max_level))
    do ilv=1,refine_max_level
      do ibc=nghostcells,1,-1
        nb=floor((gzone-dx(1,ilv)*(dble(nghostcells-ibc+1)-0.5d0))/ds+0.5d0)
        res=gzone-dx(1,ilv)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(nb)-0.5d0)*ds
        rbc(ibc,ilv)=r_o(nb)+res/ds*(r_o(nb+1)-r_o(nb))
        pbc(ibc,ilv)=p_o(nb)+res/ds*(p_o(nb+1)-p_o(nb))
      end do
    end do
    if(mype==zero) then
      print*, 'temprature from', minval(t_o), 'to', maxval(t_o)
      print*, 'density from', rhob, 'to', minval(r_o)
      print*, 'pressure from', pb, 'to', minval(p_o)
      print*, 'gravity from', minval(g_o), 'to', maxval(g_o)
    end if
  end subroutine inithdstatic

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ix^D,na
    double precision :: res

    w(ixO^S,rho_)=one
    w(ixO^S,p_)=one
    !> interpolate density and pressure
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      na=floor((x(ix^D,1)-xprobmin1+gzone)/ds+0.5d0)
      res=x(ix^D,1)-xprobmin1+gzone-(dble(na)-0.5d0)*ds
      w(ix^D,rho_)  =r_o(na)+(one-cos(dpi*res/ds))/two*(r_o(na+1)-r_o(na))
      w(ix^D,p_)    =p_o(na)+(one-cos(dpi*res/ds))/two*(p_o(na+1)-p_o(na))
    {end do\}
    w(ixO^S,mom(1))=zero
    call hd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: idir,ix^D,grsize

    grsize=nint(log(dxlv1/dxlevel(1))/log(2.d0))+1
    select case(iB)
    case(1)
      w(ixO^S,rho_)=rbc(1:nghostcells,grsize)
      w(ixO^S,p_)=pbc(1:nghostcells,grsize)
      w(ixO^S,mom(1))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,mom(1))&
                      /w(ixO^S,rho_)
      call hd_to_conserved(ixI^L,ixO^L,w,x)
    case(2)
      w(ixO^S,rho_)=rbc(nghostcells:1:-1,grsize)
      w(ixO^S,p_)=pbc(nghostcells:1:-1,grsize)
      w(ixO^S,mom(1))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,mom(1))&
                      /w(ixO^S,rho_)
      call hd_to_conserved(ixI^L,ixO^L,w,x)
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
    gravity_field(ixO^S,1)=ggrid(ixO^S)
  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)
    double precision :: res
    integer :: ix^D,na

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      na=floor((x(ix^D,1)-xprobmin1+gzone)/ds+0.5d0)
      res=x(ix^D,1)-xprobmin1+gzone-(dble(na)-0.5d0)*ds
      ggrid(ix^D)=g_o(na)+(one-cos(dpi*res/ds))/two*(g_o(na+1)-g_o(na))
    {end do\}
  end subroutine

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)
    double precision :: src(ixI^S)

    !> background heating
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    !> localiz heating
    call getlQ(lQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
    !> using the block based implicit tc instead of explicit sts tc
    if(hd_thermal_conduction .eqv. .false.) then
      call getsrc(src,qdt,ixI^L,ixO^L,w,x)
      where(w(ixO^S,e_)+src(ixO^S) .gt. zero)
        w(ixO^S,e_)=w(ixO^S,e_)+src(ixO^S)
      end where
    end if
  end subroutine special_source

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision :: bQgrid(ixI^S),lambda,bQ0

    !> temporal distribution
    !> initially a little bit higher to prevent condensation
    if(qt .lt. 10) then
      bQ0=4.d-4/heatunit
    else if(qt .lt. 20) then
      bQ0=3.d-4/heatunit
    else if(qt .lt. 30) then
      bQ0=2.d-4/heatunit
    else
      bQ0=2.d-4/heatunit
    end if

    !> spatial distribution
    lambda=xprobmax1/2.d0
    where(x(ixO^S,1) .lt. xprobmax1/2.d0)
      bQgrid(ixO^S)=bQ0*dexp(-(x(ixO^S,1)-xprobmin1)/lambda)
    elsewhere
      bQgrid(ixO^S)=bQ0*dexp(-(xprobmax1-x(ixO^S,1))/lambda)
    endwhere
  end subroutine getbQ

  subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    integer :: Bp,Bt1,Bt2,ix^D,i
    double precision :: lQgrid(ixI^S),lambda,lQ0,lH0,tramp

    lQ0=1.d-2/heatunit
    tramp=15.d0
    lH0=0.3d0
    lambda=0.3d0
    if(qt .le. trelax) then
      lQgrid=zero
    else if(qt .le. trelax+tramp) then
      lQgrid=(qt-trelax)/tramp
    else
      lQgrid=one
    end if
    where(x(ixO^S,1) .lt. lH0)
      lQgrid(ixO^S)=lQ0*lQgrid(ixO^S)
    else where(x(ixO^S,1) .lt. xprobmax1/2.d0)
      lQgrid(ixO^S)=lQ0*lQgrid(ixO^S)*dexp(-(x(ixO^S,1)-lH0)/lambda)
    else where(x(ixO^S,1) .lt. xprobmax1-lH0)
      lQgrid(ixO^S)=lQ0*lQgrid(ixO^S)*dexp(-(xprobmax1-lH0-x(ixO^S,1))/lambda)
    else where
      lQgrid(ixO^S)=lQ0*lQgrid(ixO^S)
    end where
  end subroutine getlQ

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    if(qt .lt. trelax) then
      refine=-1
      coarsen=1
    else
      if(any(x(ixO^S,1) .le. xprobmin1+0.3d0)) then
        refine=1
        coarsen=-1
      end if
      if(any(x(ixO^S,1) .ge. xprobmax1-0.3d0)) then
        refine=1
        coarsen=-1
      end if
    end if
  end subroutine special_refine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_radiative_cooling
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    integer :: ix^D
    double precision :: pth(ixI^S)
    double precision :: ens1(ixI^S),ens2(ixI^S),ens3(ixI^S)

    call hd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames

    varnames='Te'
  end subroutine specialvarnames_output

  subroutine getsrc(src,qdt,ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qdt,x(ixI^S,1:ndim),w(ixI^S,1:nw)
    double precision, intent(out) :: src(ixI^S)
    !> local
    double precision :: pth(ixI^S),Te(ixI^S),kai(ixI^S),cmat(ixI^S,1:2*ndim+1)
    double precision :: rcp(ixI^S),lcp(ixI^S),ycp(ixI^S)
    integer :: i
    logical :: first
    double precision :: dxcp,invg

    data first /.true./

    invg=1.d0/(hd_gamma-1.d0)
    dxcp=dxlevel(1)

    if(first) then
      if (mype==0) then
        print *,'solve 1D heat conduction with catch-up method'
      end if
      first=.false.
    end if

    call hd_get_pthermal(w,x,ixI^L,ixI^L,pth)
    Te(ixI^S)=pth(ixI^S)/w(ixI^S,rho_)
    do i=ixImin1,ixImax1-1
      kai(i)=kai0*dsqrt(Te(i)*Te(i+1))**2.5d0
    end do
    do i=ixImin1+1,ixImax1-1
      cmat(i,2)=-kai(i-1)/dxcp**2
      cmat(i,3)=-kai(i)/dxcp**2
      cmat(i,1)=invg*w(i,rho_)/qdt-cmat(i,2)-cmat(i,3)
      src(i)=invg*w(i,rho_)*Te(i)/qdt
    end do
    src(ixImin1+1)=src(ixImin1+1)-cmat(ixImin1+1,2)*Te(ixImin1)
    src(ixImax1-1)=src(ixImax1-1)-cmat(ixImax1-1,3)*Te(ixImax1)

    rcp(ixImin1+1)=cmat(ixImin1+1,1)
    ycp(ixImin1+1)=src(ixImin1+1)
    do i=ixImin1+2,ixImax1-1
      lcp(i)=cmat(i,2)/rcp(i-1)
      rcp(i)=cmat(i,1)-lcp(i)*cmat(i-1,3)
      ycp(i)=src(i)-lcp(i)*ycp(i-1)
    end do
    Te(ixImax1-1)=ycp(ixImax1-1)/rcp(ixImax1-1)
    do i=ixImax1-2,ixImin1+1,-1
       Te(i)=(ycp(i)-cmat(i,3)*Te(i+1))/rcp(i)
    end do
    src(ixO^S)=(w(ixO^S,rho_)*Te(ixO^S)-pth(ixO^S))*invg
  end subroutine getsrc
end module mod_usr
