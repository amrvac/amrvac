module mod_usr
  use mod_rho
  implicit none
  double precision, allocatable, save :: myBC1(:,:), myBC2(:,:)

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc => specialbound_usr
    usr_refine_grid   => specialrefine_grid

    call rho_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters
    integer :: i, j
    allocate(myBC1(domain_nx2,domain_nx3))
    allocate(myBC2(2*domain_nx2,2*domain_nx3))
    open(42,file='bc_test_1.dat')
    do i=1,domain_nx2
       do j=1,domain_nx3
          read(42,*) myBC1(i,j)
       enddo
    enddo
    close(42)
    open(42,file='bc_test_2.dat')
    do i=1,domain_nx2*2
       do j=1,domain_nx3*2
          read(42,*) myBC2(i,j)
       enddo
    enddo
    close(42)
    print*, myBC1(1,2), myBC1(2,2)
    print*, myBC1(1,1), myBC1(2,1)
    print*, ' '
    print*, myBC2(1,2), myBC2(2,2)
    print*, myBC2(1,1), myBC2(2,1)

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid 

    use mod_global_parameters

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision:: rr2(ix^S)
    double precision:: rhoflat,rhosquare,slocx^D,slocxmid^D,widthhalf

    logical::          maskv(ix^S),maska(ix^S),maskc(ix^S)
    double precision :: radius, xcircle^D
    double precision:: xc1,yc1,xa1,xa2,ya1,ya2,xb1,xb2,yb1,yb2,xc2,yc2, &
         rad,rad2,alp,nsig
    !----------------------------------------------------------------------------

    rhoflat  = 0.5d0 
    rhosquare= 2.0d0 

    ! iprob=1 is a pure 1D Riemann problem, solvable in 1D, 2D, 3D
    if (iprob==1) then
       slocx^D=0.2d0;
       where({^D&x(ix^S,^D)<=slocx^D|.and.})
          w(ix^S,rho_)     = rhosquare
       elsewhere
          w(ix^S,rho_)     = rhoflat
       endwhere
    else if (iprob==2) then
       slocxmid^D=xprobmin^D+half*(xprobmax^D-xprobmin^D);
       widthhalf=0.1d0
       where({^D&((x(ix^S,^D)>slocxmid^D-widthhalf).and.&
            (x(ix^S,^D)<slocxmid^D+widthhalf))|.and.})
          w(ix^S,rho_)     = rhosquare
       elsewhere
          w(ix^S,rho_)     = rhoflat
       endwhere
    else if (iprob==3) then
       {^IFONED call mpistop("iprob=3 is 2D!") }
       {^IFTHREED call mpistop("iprob=3 is 2D!") }
       {^IFTWOD
       xc1=0.25d0
       yc1=0.50d0
       rad=0.23d0
       rad2=0.13d0
       alp=dpi/3.0d0
       xa1=xc1
       ya1=yc1-rad
       xa2=xc1-rad*cos(alp)
       ya2=yc1+rad*sin(alp)
       xb1=xa1
       yb1=ya1
       xb2=xc1+rad*cos(alp)
       yb2=yc1+rad*sin(alp)
       xc2=xc1
       yc2=ya2+sqrt(rad2**2-(xa2-xc2)**2)
       maskv(ix^S)= ((x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2 <= rad**2) &
            .and.(x(ix^S,2)>= (ya2-ya1)*(x(ix^S,1)-xa1)/(xa2-xa1)+ya1) & 
            .and.(x(ix^S,2)>= (yb2-yb1)*(x(ix^S,1)-xb1)/(xb2-xb1)+yb1) & 
            .and.((x(ix^S,1)-xc2)**2+(x(ix^S,2)-yc2)**2 > rad2**2) 
       xc1=0.45d0
       yc1=0.475d0
       xa1=xc1
       ya1=yc1+rad
       xa2=xc1-rad*cos(alp)
       ya2=yc1-rad*sin(alp)
       xb1=xa1
       yb1=ya1
       xb2=xc1+rad*cos(alp)
       yb2=yc1-rad*sin(alp)
       xc2=xc1
       yc2=ya2-sqrt(rad2**2-(xa2-xc2)**2)
       maska(ix^S)= ((x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2 <= rad**2) &
            .and.(x(ix^S,2)<= (ya2-ya1)*(x(ix^S,1)-xa1)/(xa2-xa1)+ya1) & 
            .and.(x(ix^S,2)<= (yb2-yb1)*(x(ix^S,1)-xb1)/(xb2-xb1)+yb1) & 
            .and.((x(ix^S,1)-xc2)**2+(x(ix^S,2)-yc2)**2 > rad2**2) 
       xc1=0.75d0
       yc1=0.50d0
       alp=half*dpi-alp
       xa1=xc1-rad
       ya1=yc1
       xa2=xc1+rad*cos(alp)
       ya2=yc1+rad*sin(alp)
       xb1=xa1
       yb1=ya1
       xb2=xc1+rad*cos(alp)
       yb2=yc1-rad*sin(alp)
       yc2=yc1
       xc2=xa2+sqrt(rad2**2-(ya2-yc2)**2)
       maskc(ix^S)= ((x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2 <= rad**2) &
            .and.(x(ix^S,2)<= (ya2-ya1)*(x(ix^S,1)-xa1)/(xa2-xa1)+ya1) & 
            .and.(x(ix^S,2)>= (yb2-yb1)*(x(ix^S,1)-xb1)/(xb2-xb1)+yb1) & 
            .and.((x(ix^S,1)-xc2)**2+(x(ix^S,2)-yc2)**2 > rad2**2) 
       where(maskv(ix^S).or.maska(ix^S).or.maskc(ix^S))
          w(ix^S,rho_)     = rhosquare
       elsewhere
          w(ix^S,rho_)     = rhoflat
       endwhere
       }
    else if (iprob==4) then
       {^IFONED call mpistop("iprob=4 is 2D!") }
       {^IFTHREED call mpistop("iprob=4 is 2D!") }
       {^IFTWOD
       xc1=0.5d0
       yc1=0.5d0
       rad=0.1d0
       rhosquare=0.6d0
       rhoflat=0.5d0
       maska(ix^S)=((x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2 <= rad**2)
       where(maska(ix^S))
          w(ix^S,rho_)     = rhosquare+ &
               (rhoflat-rhosquare)*(x(ix^S,1)-xc1)**2/rad**2 + &
               (rhoflat-rhosquare)*(x(ix^S,2)-yc1)**2/rad**2  
       elsewhere
          w(ix^S,rho_)     = rhoflat
       endwhere
       }
    else if (iprob==5) then
       {^IFONED call mpistop("iprob=5 is 2D!") }
       {^IFTHREED call mpistop("iprob=5 is 2D!") }
       {^IFTWOD
       xc1=0.5d0
       yc1=0.5d0
       rad=0.1d0
       nsig=two
       rhosquare=2.0d0
       maska(ix^S)=((x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2<=(nsig*rad)**2)
       rr2(ix^S)=(x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2
       where(maska(ix^S))
          w(ix^S,rho_)     = rhosquare*exp(-rr2(ix^S)/rad**2) 
       elsewhere
          w(ix^S,rho_)     = rhosquare*exp(-nsig**2)
       endwhere
       }
    else if (iprob==6) then
       radius = 0.2d0
       xcircle^D=zero;
       where(radius**2> ^D&(x(ix^S,^D)-xcircle^D)**2+ )
          w(ix^S,rho_)     = rhosquare
       elsewhere
          w(ix^S,rho_)     = rhoflat
       endwhere
    else if (iprob==7) then
       w(ix^S,rho_)     = rhoflat
    else
       call mpistop("iprob not available!")
    end if

  end subroutine initonegrid_usr

 ! special boundary types, user defined
  ! user must assign conservative variables in bounderies
  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixG^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    logical :: patchw(ixG^T)
    double precision:: rinlet(ixG^T)

    integer :: ig2tmp, ig3tmp, lvl

    select case(iB)
    case(1)
       ! === Left boundary ===!
       lvl=int(1.01*(((xprobmax2-xprobmin2)/real(domain_nx2))/(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2))))
!       print*, lvl
!!$       print*, dxlevel(1), dxlevel(2), dxlevel(3)
!!$       print*, x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1)
!!$       print*, x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)
!!$       print*, x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3)
       ig2tmp=1+int((x(ixOmin1,ixOmin2,ixOmin3,2)-xprobmin2)/((xprobmax2-xprobmin2)/(real(domain_nx2)/real(block_nx2/(2**(lvl-1))))))
       ig3tmp=1+int((x(ixOmin1,ixOmin2,ixOmin3,3)-xprobmin3)/((xprobmax3-xprobmin3)/(real(domain_nx3)/real(block_nx3/(2**(lvl-1))))))
       if (lvl==1) then
          w(ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=myBC1(ixGmin2+(ig2tmp-1)*block_nx2:ixGmax2-2*nghostcells+(ig2tmp-1)*block_nx2,ixGmin3+(ig3tmp-1)*block_nx3:ixGmax3-2*nghostcells+(ig3tmp-1)*block_nx3)
       elseif (lvl==2) then
          w(ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=myBC2(ixGmin2+(ig2tmp-1)*block_nx2:ixGmax2-2*nghostcells+(ig2tmp-1)*block_nx2,ixGmin3+(ig3tmp-1)*block_nx3:ixGmax3-2*nghostcells+(ig3tmp-1)*block_nx3)  
!!$       elseif (lvl==3) then
!!$       w(ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=myBC1(ixGmin2+(ig2tmp-1)*block_nx2:ixGmax2-2*nghostcells+(ig2tmp-1)*block_nx2,ixGmin3+(ig3tmp-1)*block_nx3:ixGmax3-2*nghostcells+(ig3tmp-1)*block_nx3)
       endif
       !print*, ixOmin2-nghostcells, ixOmin2, ixOmax2-nghostcells, ixOmax2
       !print*, ixGmin2-nghostcells, ixGmin2, ixGmax2-nghostcells, ixGmax2
          

       !w(ixOmin1,:,:,rho_)=w(ixOmax1,:,:,rho_)
!!$       where (x(ixO^S,2)>zero .and. x(ixO^S,3)>zero)
!!$          w(ixO^S,rho_) = 0.5d0
!!$       endwhere
!!$       where (x(ixO^S,2)<zero .and. x(ixO^S,3)>zero)
!!$          w(ixO^S,rho_) = one     
!!$       endwhere     
!!$       where (x(ixO^S,2)>zero .and. x(ixO^S,3)<zero)
!!$          w(ixO^S,rho_) = 1.5d0
!!$       endwhere
!!$       where (x(ixO^S,2)<zero .and. x(ixO^S,3)<zero)
!!$          w(ixO^S,rho_) = two
!!$       endwhere

    case default
       call mpistop('boundary not defined')
    end select

  end subroutine specialbound_usr


  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    if (all(x(ix^S,2) < zero) .and. all(x(ix^S,3)>zero)) then
       refine=-1
    else
       refine=1
    endif

  end subroutine specialrefine_grid


end module mod_usr

