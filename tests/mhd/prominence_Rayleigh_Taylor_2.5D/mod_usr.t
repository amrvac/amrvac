! this setup to simulate Rayleigh-Taylor dynamics in a solar prominence
! It can be used in 2.5D or 3D, and assumes ideal MHD (with external gravity)
! One can use a tracer, and exploit GLM (or not)
! Related literature is in: 
!   `Solar prominences: "double, double ... boil and bubble"', 
!     R. Keppens, X. Cia, & O. Porth, 2015, ApJ Letters 806, L13 (7pp)
module mod_usr
  use mod_mhd
  use mod_physics
  implicit none
  double precision, allocatable, save :: pbc(:),rbc(:)
  double precision :: usr_grav,gzone,bb1,bb3,SRadius,dr,eps
  double precision :: rho0,Tch,Tco,Tpromin,Tpromax,htra1,htra2,htra3,ybot,ytop,bsca,pwidth
  integer, parameter :: jmax=8000
  double precision :: pa(jmax),ra(jmax),ya(jmax),paext(jmax),raext(jmax)
  double precision :: randphase(1000)
  integer :: nxmodes

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc => specialbound_usr
    usr_gravity=> gravity
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    call set_coordinate_system("Cartesian_2.5D")
    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters

    integer,dimension(:),allocatable:: seed
    integer::  seed_size,ix
    real:: randphaseP1(1:1000)

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar surface gravity 2.74d2 m*s^-2
    
    ! gzone gives the distance from the lower boundary xprobmin2
    ! to the so-called photosphere where density is fixed to rho0 further on
    ! this has to be wider than the ghost layer zone for the coarsest mesh....
    gzone=0.3d0 ! our grid starts above photosphere
    dr=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax)
    ! Solar radius
    SRadius=6.961d10/unit_length
    ! base density and temperature
    rho0=1.d17/unit_numberdensity
    Tch=8.d3/unit_temperature
    Tco=1.8d6/unit_temperature
    Tpromin=6.0d3/unit_temperature
    Tpromax=1.4d4/unit_temperature
    htra1=0.2d0
    htra2=1.125d0
    htra3=2.0d0
    ybot=1.75d0
    ytop=2.0d0
    bsca=1.5d0
    pwidth=0.5d0 ! width in Lunit

    nxmodes=50
    bb1=0.1d0
    bb3=4.d0
    eps=0.05d0
    
    call inithdstatic

    randphase(1:1000)=zero
    if(nxmodes>1000) call mpistop('too many modes, edit nxmodes')

    if(mype==0)then
      call random_seed(SIZE=seed_size)
      allocate(seed(seed_size))
      call random_seed(GET=seed(1:seed_size))
      call random_number(randphaseP1(1:nxmodes))
      randphase(1:nxmodes)=-dpi+two*dpi*dble(randphaseP1(1:nxmodes))
    endif
    call MPI_BARRIER(icomm,ierrmpi)
    if(npe>1)then
         call MPI_BCAST(randphase,1000,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif

    if(mype==0)then
      print *,'number of modes=',nxmodes
      open(123,file='phaseinfo',form='formatted')
      write(123,*) nxmodes
      do ix=1,nxmodes
          write(123,"(i4,1es12.4)") ix,randphase(ix)
      enddo
      close(123)
    endif

  end subroutine initglobaldata_usr

  subroutine inithdstatic
    !! initialize the table in a vertical line through the global domain
    use mod_global_parameters

    integer :: j,na,ibc,ix
    double precision:: Ta(jmax),gg(jmax),Taext(jmax)
    double precision:: res,rhob,pb,wtra1,wtra2,wtra3,b3

    wtra1=0.01d0
    wtra2=0.01d0
    wtra3=0.01d0
    do j=1,jmax
       ya(j)=(dble(j)-0.5d0)*dr-gzone+xprobmin2
       if(ya(j)<=0.4d0) then
         Ta(j)=Tch+0.5d0*(Tco-Tch)*(tanh((ya(j)-htra1)/wtra1)+1.d0)
         Taext(j)=Tch+0.5d0*(Tco-Tch)*(tanh((ya(j)-htra1)/wtra1)+1.d0)
       endif
       if(ya(j)>0.4d0.and.ya(j)<ybot) then
         Ta(j)=Tco-0.5d0*(Tco-Tpromin)*(tanh((ya(j)-htra2)/wtra2)+1.d0)
         Taext(j)=Tco
       endif
       if(ya(j)>=ybot.and.ya(j)<ytop) then
         Ta(j)=Tpromin+(Tpromax-Tpromin)*(ya(j)-ybot)/(ytop-ybot)
         Taext(j)=Tco
       endif
       if(ya(j)>ytop) then
         Ta(j)=Tpromax+0.5d0*(Tco-Tpromax)*(tanh((ya(j)-htra3)/wtra3)+1.d0)
         Taext(j)=Tco
       endif
       gg(j)=usr_grav*SRadius**2/(SRadius+ya(j))**2
    enddo
    !! solution of hydrostatic equation 
    ra(1)=rho0
    raext(1)=rho0
    pa(1)=rho0*Tch
    paext(1)=rho0*Tch
    do j=2,jmax
      if(ya(j)<htra2) then
        pa(j)=(pa(j-1)+dr*(gg(j)+gg(j-1))*ra(j-1)/4.d0)/(one-dr*(gg(j)+gg(j-1))/&
              Ta(j)/4.d0)
        paext(j)=(paext(j-1)+dr*(gg(j)+gg(j-1))*raext(j-1)/4.d0)/(one-dr*(gg(j)+gg(j-1))/&
              Taext(j)/4.d0)
      else
        if(ya(j)<ybot) then
        b3=(bb3*dexp(-(ya(j-1)-htra2)/bsca))**2+(bb3*dexp(-(ya(j)-htra2)/bsca))**2
        else
        b3=zero
        endif
        pa(j)=(pa(j-1)+dr*(gg(j)+gg(j-1))*ra(j-1)/4.d0+dr/bsca*b3/2.d0)/(one-dr*(gg(j)+gg(j-1))/&
              Ta(j)/4.d0)
        paext(j)=(paext(j-1)+dr*(gg(j)+gg(j-1))*raext(j-1)/4.d0+dr/bsca*b3/2.d0)/(one-dr*(gg(j)+gg(j-1))/&
              Taext(j)/4.d0)
      endif
      ra(j)=pa(j)/Ta(j)
      raext(j)=paext(j)/Taext(j)
    end do
    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dr+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dr
    rhob=ra(na)+res/dr*(ra(na+1)-ra(na))
    pb=pa(na)+res/dr*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dr+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dr
      rbc(ibc)=ra(na)+res/dr*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dr*(pa(na+1)-pa(na))
    end do

    if (mype==0) then
      open(123,file='pruns',form='formatted')
      write(123,*) jmax
      do ix=1,jmax
         write(123,"(i7,7es12.4)") ix,ya(ix),pa(ix),paext(ix),ra(ix),raext(ix),Ta(ix),Taext(ix)
      enddo
      close(123)
    endif

  end subroutine inithdstatic

  subroutine initonegrid_usr(ixI^L,ix^L,w,x)
    ! initialize one grid
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ix^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    
    double precision:: psi(ixI^S),tmp(ixI^S)
    double precision:: res,sigma,lxsize,sigma3
    integer :: ix^D,na,idims,imode
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          print *,'2.5D or 3D MHD Rayleigh-Taylor instability in prominence'
       end if
       first=.false.
    end if
    sigma=0.02d0
    sigma3=0.1d0

    w(ix^S,mom(1))=0.d0
    w(ix^S,mom(2))=0.d0
    w(ix^S,mom(3))=0.d0
    ! now add the incompressible perturbations
    psi(ixI^S)=zero
    lxsize=(xprobmax1-xprobmin1)
    do imode=1,nxmodes
       psi(ixI^S)=psi(ixI^S) &
           +eps*(dcos(two*dpi*dble(imode)*(x(ixI^S,1)/lxsize)+randphase(imode)) &
           /dble(imode))*dexp(-((x(ixI^S,2)-htra2)/sigma)**2) 
    enddo
    ! compute dpsi/dy
    idims=2
    select case(typegrad)
        case("central")
         call gradient(psi,ixI^L,ix^L,idims,tmp)
        case("limited")
         call gradientS(psi,ixI^L,ix^L,idims,tmp)
    end select
    w(ix^S,mom(1))=w(ix^S,mom(1))-tmp(ix^S){^IFTHREED *dexp(-(x(ix^S,3)/sigma3)**2)}
    ! compute dpsi/dx
    idims=1
    select case(typegrad)
         case("central")
          call gradient(psi,ixI^L,ix^L,idims,tmp)
         case("limited")
          call gradientS(psi,ixI^L,ix^L,idims,tmp)
    end select
    w(ix^S,mom(2))=w(ix^S,mom(2))+tmp(ix^S){^IFTHREED *dexp(-(x(ix^S,3)/sigma3)**2)}

    {do ix^DB=ixmin^DB,ixmax^DB\}
       na=floor((x(ix^D,2)-xprobmin2+gzone)/dr+0.5d0)
       res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dr
       {^IFTHREED
       if (dabs(x(ix^D,3))<pwidth*half)then
         w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dr))/two*(ra(na+1)-ra(na))
         w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dr))/two*(pa(na+1)-pa(na))
       else
         w(ix^D,rho_)=raext(na)+(one-cos(dpi*res/dr))/two*(raext(na+1)-raext(na))
         w(ix^D,p_)  =paext(na)+(one-cos(dpi*res/dr))/two*(paext(na+1)-paext(na))
       endif
       }
       {^IFTWOD
       w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dr))/two*(ra(na+1)-ra(na))
       w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dr))/two*(pa(na+1)-pa(na))
       }
    {end do\}

    w(ix^S,mag(1))  =bb1
    w(ix^S,mag(2))  =zero
    w(ix^S,mag(3))  =bb3
    where(x(ix^S,2)>htra2.and.x(ix^S,2)<ybot)
       w(ix^S,mag(3))  =bb3*dexp(-(x(ix^S,2)-htra2)/bsca)
    endwhere
    where(x(ix^S,2)>=ybot)
       w(ix^S,mag(3))  =bb3*dexp(-(ybot-htra2)/bsca)
    endwhere

    if(mhd_n_tracer>0) then
      w(ix^S,tracer(1))=zero
      where(x(ix^S,2)>htra2.and.x(ix^S,2)<htra3{^IFTHREED .and.(dabs(x(ix^S,3))<pwidth*half)})
       w(ix^S,tracer(1))=one
      endwhere
      where(x(ix^S,2)<htra1)
       w(ix^S,tracer(1))=-one
      endwhere
    end if

    if(mhd_glm) then
      w(ix^S,psi_)=0.d0
    end if

    call phys_to_conserved(ixI^L,ix^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters
    use mod_physics
    
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    
    double precision :: dx^D,gjjm1
    double precision :: Teb(ixI^S),pth(ixI^S),cg
    integer :: ix^D,idims,ixInt^L

    select case(iB)
    case(3)
      w(ixO^S,mag(1))=bb1
      w(ixO^S,mag(3))=bb3
      w(ixO^S,mag(2))=zero
      if(mhd_n_tracer>0) then
        {^IFTWOD
        do ix2=ixOmax2,ixOmin2,-1
          w(ixOmin1:ixOmax1,ix2,tracer(1))=w(ixOmin1:ixOmax1,ixOmax2+1,tracer(1))/w(ixOmin1:ixOmax1,ixOmax2+1,rho_)
        enddo
        }
        {^IFTHREED
        do ix2=ixOmax2,ixOmin2,-1
          w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,tracer(1))/w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,rho_)
        enddo
        }
      end if
      if(mhd_glm) then
        {^IFTWOD
        do ix2=ixOmax2,ixOmin2,-1
          w(ixOmin1:ixOmax1,ix2,psi_)=w(ixOmin1:ixOmax1,ixOmax2+1,psi_)
        enddo
        }
        {^IFTHREED
        do ix2=ixOmax2,ixOmin2,-1
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,psi_)=w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,psi_)
        enddo
        }
      end if
    
       {^IFTWOD
       w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(1))&
                     /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
       w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(2))&
                     /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
       w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(3))&
                     /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
       do ix2=ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
         w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
       enddo
       }
       {^IFTHREED
       w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(1))&
                     /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(2))&
                     /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(3))&
                     /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       do ix2=ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)=rbc(ix2)
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,p_)=pbc(ix2)
       enddo
       }

       call phys_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
    !! implementation of hydrostatic extrapolation at top boundary
       {^IFTWOD
       do ix2=ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,mag(1):mag(3))=(1.0d0/11.0d0)* &
              ( +2.0d0*w(ixOmin1:ixOmax1,ix2-3,mag(1):mag(3)) &
                -9.0d0*w(ixOmin1:ixOmax1,ix2-2,mag(1):mag(3)) &
               +18.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(1):mag(3)))
       enddo
       }
       {^IFTHREED
       do ix2=ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/11.0d0)* &
              ( +2.0d0*w(ixOmin1:ixOmax1,ix2-3,ixOmin3:ixOmax3,mag(1):mag(3)) &
                -9.0d0*w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,mag(1):mag(3)) &
               +18.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(1):mag(3)))
       enddo
       }
    
       do ix2=ixOmin2,ixOmax2
         w(ix2^%2ixO^S,mom(1))=w(ixOmin2-1^%2ixO^S,mom(1))/w(ixOmin2-1^%2ixO^S,rho_)
         w(ix2^%2ixO^S,mom(2))=w(ixOmin2-1^%2ixO^S,mom(2))/w(ixOmin2-1^%2ixO^S,rho_)
         w(ix2^%2ixO^S,mom(3))=w(ixOmin2-1^%2ixO^S,mom(3))/w(ixOmin2-1^%2ixO^S,rho_)
       enddo
       !! obtain the thermal pressure in top layers and store in pth
       ixInt^L=ixO^L;
       ixIntmin2=ixOmin2-nghostcells;ixIntmax2=ixOmin2-1;
       call phys_get_pthermal(w,x,ixI^L,ixInt^L,pth)
       cg=dxlevel(2)*usr_grav/two
       !! fill pth, rho ghost layers according to gravity stratification 
       {^IFTWOD
       do ix2=ixOmin2,ixOmax2
         do ix1=ixOmin1,ixOmax1
           Teb(ix1,ix2-1)=pth(ix1,ix2-1)/w(ix1,ix2-1,rho_)
           Teb(ix1,ix2)=Teb(ix1,ix2-1)
           gjjm1=half*(SRadius**2/(SRadius+x(ix1,ix2,2))**2+SRadius**2/(SRadius+x(ix1,ix2-1,2))**2)
           pth(ix1,ix2)=(pth(ix1,ix2-1)+cg*gjjm1*w(ix1,ix2-1,rho_))/(one-cg*gjjm1/Teb(ix1,ix2))
           w(ix1,ix2,rho_)=pth(ix1,ix2)/Teb(ix1,ix2)
         end do
       end do
       }
       {^IFTHREED
       do ix2=ixOmin2,ixOmax2
         do ix1=ixOmin1,ixOmax1
         do ix3=ixOmin3,ixOmax3
           Teb(ix1,ix2-1,ix3)=pth(ix1,ix2-1,ix3)/w(ix1,ix2-1,ix3,rho_)
           Teb(ix1,ix2,ix3)=Teb(ix1,ix2-1,ix3)
           gjjm1=half*(SRadius**2/(SRadius+x(ix1,ix2,ix3,2))**2+SRadius**2/(SRadius+x(ix1,ix2-1,ix3,2))**2)
           pth(ix1,ix2,ix3)=(pth(ix1,ix2-1,ix3)+cg*gjjm1*w(ix1,ix2-1,ix3,rho_))/(one-cg*gjjm1/Teb(ix1,ix2,ix3))
           if(pth(ix1,ix2,ix3)<minp) pth(ix1,ix2,ix3)=pth(ix1,ix2-1,ix3)
           w(ix1,ix2,ix3,rho_)=pth(ix1,ix2,ix3)/Teb(ix1,ix2,ix3)
         end do
         end do
       end do
       }
       w(ixO^S,p_)=pth(ixO^S)
    {#IFDEF TRACER
       {^IFTWOD
       do ix2=ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,tracer(1))=w(ixOmin1:ixOmax1,ixOmin2-1,tracer(1))/w(ixOmin1:ixOmax1,ixOmin2-1,rho_)
       enddo
       }
       {^IFTHREED
       do ix2=ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,tracer(1))= &
         w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,tracer(1))/w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,rho_)
       enddo
       }
    }
    {#IFDEF GLM
       {^IFTWOD
       do ix2=ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,psi_)=w(ixOmin1:ixOmax1,ixOmin2-1,psi_)
       enddo
       }
       {^IFTHREED
       do ix2=ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,psi_)=w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,psi_)
       enddo
       }
    }
       call phys_to_conserved(ixI^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  ! Calculate gravitational acceleration in each dimension
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:)=0.d0
    gravity_field(ixO^S,2)=usr_grav*(SRadius/(SRadius+x(ixO^S,2)))**2

  end subroutine gravity

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio),tmp(ixI^S)
    double precision                   :: normconv(0:nw+nwauxio)

    call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)

    ! output the temperature p/rho
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+2)=tmp(ixO^S)*two/sum(w(ixO^S,mag(:))**2,dim=ndim+1)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='T beta'

  end subroutine specialvarnames_output

end module mod_usr
