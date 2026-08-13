! 03/08/2022 - Copying from mhd/solar_atmosphere_2.5D.t
! 05/10/2022 - Setting the B field splitting with 'linedpowel'

module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav, trelax
  double precision :: heatunit,gzone,B0,theta,SRadius,kx,dya  !bQ0,
  double precision :: y0, Bpo, Bzo
  double precision :: reference_field=20.d0
  double precision :: kink_amplitude=20.d0
  double precision, allocatable :: pa(:),ra(:)
  integer, parameter :: jmax=8000

  ! Storing additional variables in the dat file
  integer :: Temp_, b1t_, b2t_, b3t_, j3_
  integer :: Alfv_, vsound_, divB_, beta_
  integer :: bQ_, lQ_, rad_, thcond_
  integer :: Qk_, QAW_

contains

  subroutine usr_init()

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_set_B0          => specialset_B0
    usr_modify_output   => set_output_vars
    usr_uawsom_coefficients => solar_uawsom_coefficients
    usr_print_log       => solar_uawsom_log

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    call set_coordinate_system("Cartesian_2.5D")
    call mhd_activate()

    Temp_ = var_set_extravar("Te", "Te")
    b1t_ = var_set_extravar("b1t", "b1t")
    b2t_ = var_set_extravar("b2t", "b2t")
    b3t_ = var_set_extravar("b3t", "b3t")
    !j3_ = var_set_extravar("j3", "j3")

    !Alfv_ = var_set_extravar("Alfv", "Alfv")
    !vsound_ = var_set_extravar("vsound", "vsound")
    !divB_ = var_set_extravar("divB", "divB")
    !beta_ = var_set_extravar("beta", "beta")

    !bQ_ = var_set_extravar("bQ", "bQ")
    !lQ_ = var_set_extravar("lQ", "lQ")
    rad_ = var_set_extravar("rad", "rad")
    !thcond_ = var_set_extravar("thcond", "thcond")
    Qk_ = var_set_extravar("Qk", "Qk")
    QAW_ = var_set_extravar("QAW", "QAW")

  end subroutine usr_init

  subroutine initglobaldata_usr()
    call read_usr_parameters(par_files)
    heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    ! bQ0=1.d-4/heatunit ! background heating power density
    ! Extend the atmosphere table far enough to cover all bottom ghost cells.
    ! This remains 0.2 for the full tutorial grid and grows for the CI grid.
    gzone=max(0.2d0,nghostcells*dx(2,refine_max_level))
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
    B0=reference_field/unit_magneticfield ! magnetic field strength at the bottom
    Bpo = B0     !/(dsqrt(2.d0)*(dexp(-1.d0)-dexp(-3.d0)))
    Bzo = B0     !/dsqrt(2.d0)
    y0 = -0.4d0 !=4Mm as in Zhang+2019
    kx=dpi/(xprobmax1-xprobmin1)
    SRadius=69.61d0 ! Solar radius
    trelax=90.d0

    ! hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic

    if(mype .eq. 0) then
      print*, 'unit_density = ', unit_density
      print*, 'unit_pressure = ', unit_pressure
      print*, 'unit_velocity = ', unit_velocity
      print*, 'unit_magneticfield = ', unit_magneticfield
      print*, 'unit_time = ', unit_time
      print*, 'usr_grav = ', usr_grav
    end if

  end subroutine initglobaldata_usr

  subroutine read_usr_parameters(files)
    character(len=*), intent(in) :: files(:)
    integer :: n
    namelist /usr_list/ kink_amplitude, reference_field

    kink_amplitude=20.d0
    reference_field=20.d0
    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,usr_list,end=100)
100   close(unitpar)
    end do
  end subroutine read_usr_parameters

  subroutine inithdstatic
    use mod_solar_atmosphere
    ! initialize the table in a vertical line through the global domain
    integer :: j,na,ibc
    double precision, allocatable :: Ta(:),gg(:),ya(:)
    double precision :: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa
    double precision :: rhohc,hc

    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))

    rpho=1.151d15/unit_numberdensity ! number density at the bottom of height table
    Tpho=8.d3/unit_temperature ! temperature of chromosphere
    Ttop=1.5d6/unit_temperature ! estimated temperature in the top
    htra=0.2d0 ! height of initial transition region
    wtra=0.02d0 ! width of initial transition region
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
    ! initialize one grid
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: MAX0wB02(ixO^S), MIN0wB02(ixO^S)
    double precision :: res
    integer :: ix1,ix2,na
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
    w(ixO^S,mag(:))=zero

    MAX0wB02(ixO^S) = max(0.0d0,-Bpo*dsin(kx*x(ixO^S,1))*dexp(-kx*(x(ixO^S,2)-y0)) + Bpo*dsin(3*kx*x(ixO^S,1))*dexp(-3*kx*(x(ixO^S,2)-y0)))
    MIN0wB02(ixO^S) = min(0.0d0,-Bpo*dsin(kx*x(ixO^S,1))*dexp(-kx*(x(ixO^S,2)-y0)) + Bpo*dsin(3*kx*x(ixO^S,1))*dexp(-3*kx*(x(ixO^S,2)-y0)))

    !> Set the initial profile of wave energy. Note that this should match with the boundary injection also.
    !> Small initial wave energy everywhere (1.d-6) to aid numerical stability.
    w(ixO^S,wkminus_) = 1.d-6 + kink_amplitude*dexp(-(x(ixO^S,2)-(xprobmin2))/0.1d0)*merge(1.d0, 0.d0, MAX0wB02(ixO^S) > 0.d0)
    w(ixO^S,wkplus_) = 1.d-6 + kink_amplitude*dexp(-(x(ixO^S,2)-(xprobmin2))/0.1d0)*merge(1.d0, 0.d0, MIN0wB02(ixO^S) < 0.d0)

    w(ixO^S,wAminus_) = 0.0d0 !1.d-6 + 0.7d0*dexp(-(x(ixO^S,2)-(xprobmin2))/0.1d0)*merge(1.d0, 0.d0, MAX0wB02(ixO^S) > 0.d0)
    w(ixO^S,wAplus_) = 0.0d0  !1.d-6 + 0.7d0*dexp(-(x(ixO^S,2)-(xprobmin2))/0.1d0)*merge(1.d0, 0.d0, MIN0wB02(ixO^S) < 0.d0)

    call eos%to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a steady (time-independent) potential or
  ! linear force-free background field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixO^S,1)=Bpo*dcos(kx*x(ixO^S,1))*dexp(-kx*(x(ixO^S,2)-y0)) - Bpo*dcos(3*kx*x(ixO^S,1))*dexp(-3*kx*(x(ixO^S,2)-y0))
    wB0(ixO^S,2)=-Bpo*dsin(kx*x(ixO^S,1))*dexp(-kx*(x(ixO^S,2)-y0)) + Bpo*dsin(3*kx*x(ixO^S,1))*dexp(-3*kx*(x(ixO^S,2)-y0))
    wB0(ixO^S,3)=Bzo

  end subroutine specialset_B0

  subroutine specialbound_usr(qdt,qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qdt, qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S), MAX0wB02G(ixI^S),MIN0wB02G(ixI^S)
    double precision :: Q(ixI^S),Qp(ixI^S),zeta(ixI^S)
    integer          :: ix^D,ixOs^L,ixC^L,hxC^L,jxO^L,idir
    double precision :: wB0(ixI^S,1:ndir)

    MAX0wB02G(ixI^S) = max(0.0d0,-Bpo*dsin(kx*x(ixI^S,1))*dexp(-kx*(x(ixI^S,2)-y0)) + Bpo*dsin(3*kx*x(ixI^S,1))*dexp(-3*kx*(x(ixI^S,2)-y0)))
    MIN0wB02G(ixI^S) = min(0.0d0,-Bpo*dsin(kx*x(ixI^S,1))*dexp(-kx*(x(ixI^S,2)-y0)) + Bpo*dsin(3*kx*x(ixI^S,1))*dexp(-3*kx*(x(ixI^S,2)-y0)))

    ! iB	integer indicating direction of boundary
    select case(iB)
    case(3)
      !! Fixed zero velocity:
      do idir=1,ndir
        w(ixO^S,mom(idir))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))&
                   /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end do

      !> Set the injection of wave energy
      do ix1=ixOmin1,ixOmax1
        w(ix1,ixOmin2,wkminus_) = kink_amplitude*merge(1.d0, 0.d0, MAX0wB02G(ix1,ixOmin2) > 0.d0) !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2))
        w(ix1,ixOmax2,wkminus_) = kink_amplitude*merge(1.d0, 0.d0, MAX0wB02G(ix1,ixOmin2) > 0.d0)  !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2))

        w(ix1,ixOmin2,wkplus_) =  kink_amplitude*merge(1.d0, 0.d0, MIN0wB02G(ix1,ixOmin2) < 0.d0)  !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))
        w(ix1,ixOmax2,wkplus_) =  kink_amplitude*merge(1.d0, 0.d0, MIN0wB02G(ix1,ixOmin2) < 0.d0)  !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))

        !w(ix1,ixOmin2,wAminus_) = 0.7d0*merge(1.d0, 0.d0, MAX0wB02G(ix1,ixOmin2) > 0.d0)  !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2))
        !w(ix1,ixOmax2,wAminus_) = 0.7d0*merge(1.d0, 0.d0, MAX0wB02G(ix1,ixOmin2) > 0.d0)  !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2))

        !w(ix1,ixOmin2,wAplus_) =  0.7d0*merge(1.d0, 0.d0, MIN0wB02G(ix1,ixOmin2) < 0.d0)  !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))
        !w(ix1,ixOmax2,wAplus_) =  0.7d0*merge(1.d0, 0.d0, MIN0wB02G(ix1,ixOmin2) < 0.d0)  !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))

      end do

      w(ixO^S,wAminus_) = 0.0d0
      w(ixO^S,wAplus_) = 0.0d0

      !! fixed b1 b2 b3
      w(ixO^S,mag(:))=0.d0

      !! fixed gravity stratification of density and pressure pre-determined in initial condition
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
      enddo
      call eos%to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ixOs^L=ixO^L;
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmin2-1;
      call mhd_get_pthermal(w,x,ixI^L,ixOs^L,pth)
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmax2;
      call getggrav(ggrid,ixI^L,ixOs^L,x)

      !> Fill pth, rho ghost layers according to gravity stratification:
      invT(ixOmin2-1^%2ixO^S)=w(ixOmin2-1^%2ixO^S,rho_)/pth(ixOmin2-1^%2ixO^S)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin2-1^%2ixO^S)=tmp(ixOmin2-1^%2ixO^S)+0.5d0*&
            (ggrid(ix2^%2ixO^S)+ggrid(ix2-1^%2ixO^S))*invT(ixOmin2-1^%2ixO^S)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)*dexp(tmp(ixOmin2-1^%2ixO^S)*dxlevel(2))
        w(ix2^%2ixO^S,rho_)=w(ix2^%2ixO^S,p_)*invT(ixOmin2-1^%2ixO^S)
      enddo

      !> Fixed zero velocity:
      do idir=1,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))&
                     /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      end do

      !> Magnetic field:
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* &
                    (-w(ixOmin1:ixOmax1,ix2-2,mag(:))&
               +4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
      enddo
      call eos%to_conserved(ixI^L,ixO^L,w,x)

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
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: bQgrid(ixI^S)

    ! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    ! add steady localized heating
    !call getlQ(lQgrid,ixI^L,ixO^L,qtC,wCT,x)
    !w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)

  end subroutine special_source

  !> Paper Eq. 30: the density contrast equals mhd_uawsom_zeta0 at the base.
  subroutine get_zeta(w,x,ixI^L,ixO^L,zeta)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: zeta(ixI^S)
    zeta(ixO^S) = 1.d0+(mhd_uawsom_zeta0-1.d0)*&
         exp(-max(x(ixO^S,2)-xprobmin2,0.d0)/(5.d0*SRadius))

  end subroutine get_zeta

  !> McMurdo et al. (2026) closure profiles, Eqs. 29--30 and the paragraph
  !> following Eq. 9.  Lengths are returned in code units.
  subroutine solar_uawsom_coefficients(w,x,ixI^L,ixO^L,primitive,zeta,radius,lperp_A)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    logical, intent(in) :: primitive
    double precision, intent(out) :: zeta(ixI^S),radius(ixI^S),lperp_A(ixI^S)
    double precision :: Btotal(ixI^S,1:ndir),Bmag(ixI^S),Te(ixI^S),Rfactor(ixI^S)
    integer :: idir

    call get_zeta(w,x,ixI^L,ixO^L,zeta)
    do idir=1,ndir
      if(B0field) then
        Btotal(ixO^S,idir)=w(ixO^S,mag(idir))+block%B0(ixO^S,idir,b0i)
      else
        Btotal(ixO^S,idir)=w(ixO^S,mag(idir))
      end if
    end do
    Bmag(ixO^S)=dsqrt(sum(Btotal(ixO^S,1:ndir)**2,dim=ndim+1))
    ! Paper setup: B0=20 G and R0=1 Mm.
    radius(ixO^S)=1.d8/unit_length*&
         dsqrt((reference_field/unit_magneticfield)/max(Bmag(ixO^S),smalldouble))
    if(primitive) then
      call eos%get_Rfactor(w,x,ixI^L,ixO^L,Rfactor)
      Te(ixO^S)=w(ixO^S,p_)/(w(ixO^S,rho_)*Rfactor(ixO^S))
    else
      call eos%get_Te(w,x,ixI^L,ixO^L,Te)
    end if
    ! L_perp,AW = 1e8 sqrt(T_MK/B_G) m.
    lperp_A(ixO^S)=1.d10/unit_length*dsqrt(&
         max(Te(ixO^S)*unit_temperature/1.d6,smalldouble)/&
         max(Bmag(ixO^S)*unit_magneticfield,smalldouble))
  end subroutine solar_uawsom_coefficients

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
  ! calculate background heating bQ
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S),bQ0

    !> temporal distribution
    !> initially a little bit higher to prevent condensation
    if(qt .lt. 10) then
      bQ0=4.d-4/heatunit
    else if(qt .lt. 20) then
      bQ0=3.d-4/heatunit
    else if(qt .lt. 30) then
      bQ0=2.d-4/heatunit
    else if(qt .lt. 40) then
      bQ0=1.d-4/heatunit
    else
      bQ0=0.d-4/heatunit
    end if

    ! spatial distribution
    bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,2)/5.d0)

  end subroutine getbQ

  subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    integer :: ix^D,i
    double precision :: lQgrid(ixI^S),lQ0,tramp
    double precision :: yh,lambdah,sigma2,xr,xl

    lQ0 = 2.d-2/heatunit
    tramp = 10.d0

    yh = 0.4d0
    lambdah = 0.25d0
    sigma2 = 0.2d0

    xr = 4.2d0
    xl = -xr

    if(qt .le. trelax) then
      lQgrid=zero
    else if(qt .le. trelax+tramp) then
      lQgrid=dsin((qt-trelax)/tramp*dpi/2.d0)
    else
      lQgrid=one
    end if

    where (x(ixO^S,2) .lt. yh)
      lQgrid(ixO^S) = lQ0*lQgrid(ixO^S)*(dexp(-(x(ixO^S,1)-xr)**2/sigma2) + dexp(-(x(ixO^S,1)-xl)**2/sigma2))
    elsewhere
      lQgrid(ixO^S) = lQ0*lQgrid(ixO^S)*dexp(-(x(ixO^S,2)-yh)**2/lambdah)* &
                      (dexp(-(x(ixO^S,1)-xr)**2/sigma2) + dexp(-(x(ixO^S,1)-xl)**2/sigma2))
    end where

  end subroutine getlQ

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    if(any(w(ixI^S,rho_) .ge. 5.d0)) then
      refine=1
      coarsen=-1
    end if

  end subroutine special_refine_grid

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_radiative_cooling
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: wlocal(ixI^S,1:nw),Te(ixI^S),pth(ixI^S)
    double precision :: Btotal(ixI^S,1:ndir),rc(ixI^S)
    integer          :: idir
    double precision :: zeta(ixI^S), radius(ixI^S), Lperp_AW(ixI^S), Lperp(ixI^S), Gamma_plus(ixI^S), Gamma_minus(ixI^S)

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    !   output temperature
    call mhd_get_pthermal(wlocal,x,ixI^L,ixI^L,pth)
    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)
    w(ixO^S,Temp_)=Te(ixO^S)

    do idir=1,ndir
      if(B0field) then
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0)
      else
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))
      endif
    end do
    !   store total magnetic field
    w(ixO^S,b1t_)=Btotal(ixO^S,1)
    w(ixO^S,b2t_)=Btotal(ixO^S,2)
    w(ixO^S,b3t_)=Btotal(ixO^S,3)

    !   B^2
    !B2(ixI^S)=sum((Btotal(ixI^S,:))**2,dim=ndim+1)
    !Bmag(ixI^S)=sqrt(B2(ixI^S))
    !pmag(ixI^S)=B2(ixI^S)/2.0d0

    !   store current
    !call get_current(wlocal,ixI^L,ixO^L,idirmin,curlvec)
    !w(ixO^S,j3_)=curlvec(ixO^S,3)

    !    output Alfven wave speed B/sqrt(rho)
    !w(ixO^S,Alfv_)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))

    !    output the sound speed sqrt(gamm p/rho)
    !w(ixO^S,vsound_)= sqrt(uawsom_gamma * w(ixO^S,p_)/w(ixO^S,rho_))

    !   output divB1
    !call get_divb(wlocal,ixI^L,ixO^L,divb)
    !w(ixO^S,divB_)=divb(ixO^S)

    ! output the plasma beta p*2/B**2
    !w(ixO^S,beta_)=pth(ixO^S)*two/B2(ixO^S)

    !  store the cooling rate
    if(mhd_radiative_cooling) call getvar_cooling(ixI^L,ixO^L,wlocal,x,rc,rc_fl)
    w(ixO^S,rad_)=rc(ixO^S)

    ! output thermal conduction TC (disabled in this tutorial output set)
    !tc_fl%e_=1
    !call uawsom_sts_set_source_tc_mhd(ixI^L,ixO^L,wlocal,x,tc,.false.,1.d0,1,1,tc_fl)
    !w(ixO^S,thcond_)=tc(ixO^S)
    !tc_fl%e_=e_

    !    output heating rate
    !call getbQ(ens,ixI^L,ixO^L,global_time,wlocal,x)
    !w(ixO^S,bQ_)=ens(ixO^S)

    !   output local heating
    !call getlQ(loc_heat,ixI^L,ixO^L,global_time,wlocal,x)
    !w(ixO^S,lQ_)=loc_heat(ixO^S)

    !> Wave heating rates are calculated here
    !   output kink wave heating rate
    call solar_uawsom_coefficients(wlocal,x,ixI^L,ixO^L,.false.,&
         zeta,radius,Lperp_AW)
    Lperp(ixO^S) = (zeta(ixO^S) + 1.d0 - mhd_uawsom_filling_factor)**(3.d0/2.d0)/(1.d0 - mhd_uawsom_filling_factor**(5.d0/2.d0))/&
                   (zeta(ixO^S) - 1.d0)*3.1622776*(mhd_uawsom_filling_factor*dpi)**0.5d0*radius(ixO^S)

    w(ixO^S,Qk_)=w(ixO^S,wkplus_)**(3.d0/2.d0)/(w(ixO^S,rho_)*(1+mhd_uawsom_filling_factor*zeta(ixO^S)-mhd_uawsom_filling_factor)**(-1.d0))**0.5d0/Lperp(ixO^S) +&
                  w(ixO^S,wkminus_)**(3.d0/2.d0)/(w(ixO^S,rho_)*(1+mhd_uawsom_filling_factor*zeta(ixO^S)-mhd_uawsom_filling_factor)**(-1.d0))**0.5d0/Lperp(ixO^S)


    !   output Alfven wave heating rate /
    if (any(w(ixO^S,wAminus_) .le. 0.0d0 .or. w(ixO^S,wAplus_) .le. 0.0d0)) then
      w(ixO^S,QAW_) = 0.0d0
    else
      Gamma_plus(ixO^S) = (2.0d0 / Lperp_AW(ixO^S)) * (w(ixO^S, wAminus_)/w(ixO^S,rho_))**0.5d0
      Gamma_minus(ixO^S) = (2.0d0 / Lperp_AW(ixO^S)) * (w(ixO^S, wAplus_)/w(ixO^S,rho_))**0.5d0
      w(ixO^S,QAW_) = Gamma_plus(ixO^S)*w(ixO^S,wAplus_) + Gamma_minus(ixO^S)*w(ixO^S,wAminus_)
    end if

  end subroutine set_output_vars

  subroutine solar_uawsom_log()
    use mod_global_parameters
    integer, parameter :: log_unit=732
    integer :: iigrid,igrid
    double precision :: local_sum(5),global_sum(5),local_min(2),global_min(2)
    double precision :: volume(ixG^T),pth(ixG^T)
    character(len=std_len) :: filename
    logical, save :: first=.true.

    local_sum=zero
    local_min=bigdouble
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      block=>ps(igrid)
      volume(ixM^T)=block%dvolume(ixM^T)
      local_sum(1)=local_sum(1)+sum(block%w(ixM^T,wAplus_)*volume(ixM^T))
      local_sum(2)=local_sum(2)+sum(block%w(ixM^T,wAminus_)*volume(ixM^T))
      local_sum(3)=local_sum(3)+sum(block%w(ixM^T,wkplus_)*volume(ixM^T))
      local_sum(4)=local_sum(4)+sum(block%w(ixM^T,wkminus_)*volume(ixM^T))
      local_sum(5)=local_sum(5)+sum(block%w(ixM^T,e_)*volume(ixM^T))
      local_min(1)=min(local_min(1),minval(block%w(ixM^T,wAplus_)),&
           minval(block%w(ixM^T,wAminus_)),minval(block%w(ixM^T,wkplus_)),&
           minval(block%w(ixM^T,wkminus_)))
      call mhd_get_pthermal(block%w,block%x,ixG^LL,ixM^LL,pth)
      local_min(2)=min(local_min(2),minval(pth(ixM^T)))
    end do
    call MPI_ALLREDUCE(local_sum,global_sum,5,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(local_min,global_min,2,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)

    if(mype==0) then
      filename=trim(base_filename)//'.log'
      if(first .and. restart_from_file==undefined) then
        open(log_unit,file=trim(filename),status='replace')
        write(log_unit,'(a)') '# time WA+ WA- Wk+ Wk- Etot minW minP'
      else
        open(log_unit,file=trim(filename),status='old',position='append')
      end if
      write(log_unit,'(8(es16.8,1x))') global_time,global_sum,global_min
      close(log_unit)
    end if
    first=.false.
  end subroutine solar_uawsom_log

end module mod_usr
