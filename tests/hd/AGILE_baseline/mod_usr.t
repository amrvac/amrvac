module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr

    call hd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision                :: rhodens,rholight,kx,kz,pint,dlin,vextra,sigma
    logical                         :: first = .true.

    rhodens=10.0d0
    rholight=1.0d0
    ! kx=2 pi
    kx=4.0d0*dpi
    {^IFTHREED
    ! kz=8 pi
    kz=32d0*atan(one)
    }

    select case(iprob)
    case(1)
       vextra=0.0d0
    case(2)
       vextra=10.0d0
    case(3)
       vextra=100.0d0
    endselect

    if (first) then
       if (mype==0) then
          print *,'3D HD KH'
          print *,'  --assuming y ranging from 0-1!'
          print *,'  --density ratio:',rhodens/rholight
          print *,'  --kx:',kx
          {^IFTHREED
          print *,'  --kz:',kz
          }
          print *,'  --vextra:',vextra
       end if
       first=.false.
    end if

    ! pressure at interface
    pint=2.5d0

    dlin=0.025d0
    sigma=0.00125d0

    where(x(ixG^S,2)<0.25d0.or.x(ixG^S,2)>0.75d0)
       w(ixG^S,rho_)=rholight
       !w(ixG^S,tr1_)=0.0d0
    elsewhere
       w(ixG^S,rho_)=rhodens
       !w(ixG^S,tr1_)=1.0d0
    endwhere

    w(ixG^S,p_)=pint

    where((x(ixG^S,2)<=(0.25d0-dlin)).or.(x(ixG^S,2)>=(0.75d0+dlin)))
       w(ixG^S,mom(1))=vextra-0.5d0
    endwhere

    where((x(ixG^S,2)>=(0.25d0+dlin)).and.(x(ixG^S,2)<=(0.75d0-dlin)))
       w(ixG^S,mom(1))=vextra+0.5d0
    endwhere

    where((x(ixG^S,2)>(0.25d0-dlin)).and.(x(ixG^S,2)<(0.25d0+dlin)))
       w(ixG^S,mom(1))=(x(ixG^S,2)-(0.25d0-dlin))* &
            ((vextra+0.5d0)-(vextra-0.5d0))/(2.0d0*dlin)+(vextra-0.5d0)
    endwhere

    where((x(ixG^S,2)>(0.75d0-dlin)).and.(x(ixG^S,2)<(0.75d0+dlin)))
       w(ixG^S,mom(1))=(x(ixG^S,2)-(0.75d0-dlin))* &
            ((vextra-0.5d0)-(vextra+0.5d0))/(2.0d0*dlin)+(vextra+0.5d0)
    endwhere

    w(ixG^S,mom(2))=0.01d0*dsin(kx*x(ixG^S,1))* &
         (dexp(-0.5d0*(x(ixG^S,2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixG^S,2)-0.75d0)**2/sigma))
    {^IFTHREED
    w(ixG^S,mom(3))=0.1d0*dsin(kz*x(ixG^S,3))* &
         (dexp(-0.5d0*(x(ixG^S,2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixG^S,2)-0.75d0)**2/sigma))
    }

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

end module mod_usr
