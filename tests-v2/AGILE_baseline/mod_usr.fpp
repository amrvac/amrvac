module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr

    call set_coordinate_system("Cartesian")
    call hd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision                :: rhodens,rholight,kx,kz,pint,dlin,vextra,&
       sigma
    logical                         :: first = .true.

    rhodens=10.0d0
    rholight=1.0d0
    ! kx=2 pi
    kx=4.0d0*dpi
    
    ! kz=8 pi
    kz=32d0*atan(one)
   

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
          print *,'2D HD KH'
          print *,'  --assuming y ranging from 0-1!'
          print *,'  --density ratio:',rhodens/rholight
          print *,'  --kx:',kx
          
          print *,'  --kz:',kz
         
          print *,'  --vextra:',vextra
       end if
       first=.false.
    end if

    ! pressure at interface
    pint=2.5d0

    dlin=0.025d0
    sigma=0.00125d0

    where(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)<0.25d0.or.x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>0.75d0)
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,rho_)=rholight
       !w(ixG^S,tr1_)=0.0d0
    elsewhere
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,rho_)=rhodens
       !w(ixG^S,tr1_)=1.0d0
    endwhere

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,p_)=pint

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)<=(0.25d0-dlin)).or.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)>=(0.75d0+dlin)))
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(1))=vextra-0.5d0
    endwhere

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>=(0.25d0+dlin)).and.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)<=(0.75d0-dlin)))
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(1))=vextra+0.5d0
    endwhere

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>(0.25d0-dlin)).and.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)<(0.25d0+dlin)))
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
          mom(1))=(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
          2)-(0.25d0-dlin))* ((vextra+0.5d0)-(vextra-0.5d0))/(2.0d0*dlin)+&
          (vextra-0.5d0)
    endwhere

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>(0.75d0-dlin)).and.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)<(0.75d0+dlin)))
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
          mom(1))=(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
          2)-(0.75d0-dlin))* ((vextra-0.5d0)-(vextra+0.5d0))/(2.0d0*dlin)+&
          (vextra+0.5d0)
    endwhere

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       mom(2))=0.01d0*dsin(kx*x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1))* (dexp(-0.5d0*(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2,ixGmin3:ixGmax3,2)-0.75d0)**2/sigma))
    
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       mom(3))=0.1d0*dsin(kz*x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       3))* (dexp(-0.5d0*(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)-0.75d0)**2/sigma))
   

    call hd_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

end module mod_usr
