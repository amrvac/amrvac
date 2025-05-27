module mod_interpolation
  use mod_global_parameters
  use mod_comm_lib, only: mpistop
  implicit none

contains

  subroutine interp_linear(x_table,y_table,n_table,x_itp,y_itp,n_itp)
    ! linear interpolation
    ! 1D method

    integer :: n_table,n_itp
    double precision :: x_table(n_table),y_table(n_table)
    double precision :: x_itp(n_itp),y_itp(n_itp)

    integer :: i,j

    !if (x_table(1)>x_itp(1) .or. x_table(n_table)<x_itp(n_itp)) then
    !  call mpistop("out of interpolation table")
    !endif

    do i=1,n_itp
      if(x_itp(i)<x_table(1)) then
        y_itp(i)=y_table(1)
        exit
      end if
      if(x_itp(i)>=x_table(n_table)) then
        y_itp(i)=y_table(n_table)
        exit
      end if
      LOOPT: do j=1,n_table-1
        if (x_itp(i)>=x_table(j) .and. x_itp(i)<x_table(j+1)) then
          y_itp(i)=y_table(j) + (x_itp(i)-x_table(j))*(y_table(j+&
             1)-y_table(j))/(x_table(j+1)-x_table(j))
          exit LOOPT
        endif
      enddo LOOPT
    enddo

  end subroutine interp_linear

  subroutine interp_cubic_spline(x_table,y_table,n_table,x_itp,y_itp,n_itp)
    ! interpolation function fi=ai+bi*(x-xi)+ci*(x-xi)^2+di*(x-xi)^3
    ! first order derivative and second order derivative is continous 
    ! 1D method

    integer :: n_table,n_itp
    double precision :: x_table(n_table),y_table(n_table)
    double precision :: x_itp(n_itp),y_itp(n_itp)

    double precision :: at(n_table),bt(n_table),ct(n_table),dt(n_table)
    double precision :: dxi
    integer :: i,j

    if (x_table(1)>x_itp(1) .or. x_table(n_table)<x_itp(n_itp)) then
      call mpistop("out of interpolation table")
    endif

    call get_cubic_para(x_table,y_table,at,bt,ct,dt,n_table)

    do i=1,n_itp
      LOOPT: do j=1,n_table-1
        if (x_itp(i)>=x_table(j) .and. x_itp(i)<x_table(j+1)) then
          dxi=x_itp(i)-x_table(j)
          y_itp(i)=at(j)+bt(j)*dxi+ct(j)*dxi**2+dt(j)*dxi**3
          exit LOOPT
        endif
      enddo LOOPT
    enddo

    if (x_itp(n_itp)==x_table(n_table)) y_itp(n_itp)=y_table(n_table)

  end subroutine interp_cubic_spline

  subroutine get_cubic_para(xi,yi,ai,bi,ci,di,ni)
  ! get parameters for interpolation

    integer :: ni,i
    double precision :: xi(ni),yi(ni)
    double precision :: ai(ni),bi(ni),ci(ni),di(ni)

    double precision :: matrix1(ni,ni),ri1(ni)
    double precision :: matrix2(ni,ni),ri2(ni)
    double precision :: mi(ni),hi(ni)

    ! build equations
    matrix1=0.d0
    matrix2=0.d0

    do i=1,ni-1
      hi(i)=xi(i+1)-xi(i)
    enddo

    matrix1(1,1)=1.0
    matrix1(ni,ni)=1.0
    ri1(1)=0.d0
    ri1(ni)=0.d0
    do i=2,ni-1
      matrix1(i-1,i)=hi(i-1)
      matrix1(i,i)=2*(hi(i-1)+hi(i))
      matrix1(i+1,i)=hi(i)
      ri1(i)=(yi(i+1)-yi(i))/hi(i)-(yi(i)-yi(i-1))/hi(i-1)
    enddo
    ri1=ri1*6

    ! solve equations
    do i=1,ni
      matrix2(i,i)=1.d0
    enddo

    matrix2(2,1)=matrix1(2,1)/matrix1(1,1)
    ri2(1)=ri1(1)/matrix1(1,1)

    do i=2,ni-1
      matrix2(i+1,i)=matrix1(i+1,i)/(matrix1(i,i)-matrix1(i-1,i)*matrix2(i,&
         i-1))
    enddo

    do i=2,ni
      ri2(i)=ri1(i)-matrix1(i-1,i)*ri2(i-1)
      ri2(i)=ri2(i)/(matrix1(i,i)-matrix1(i-1,i)*matrix2(i,i-1))
    enddo

    mi(ni)=ri2(ni)
    do i=ni-1,1,-1
      mi(i)=ri2(i)-matrix2(i+1,i)*mi(i+1)
    enddo

    ! get parameters for interpolation
    ai=yi
    ci=mi(1:ni)
    do i=1,ni-1
      bi(i)=(yi(i+1)-yi(i))/hi(i)-hi(i)*mi(i)/2.0-hi(i)*(mi(i+1)-mi(i))/6.0
      di(i)=(mi(i+1)-mi(i))/(6.0*hi(i))
    enddo

  end subroutine get_cubic_para

  subroutine get_interp_factor_linear(xc,xp,factor)
    ! get factor for linear interpolation
    ! multi-D method
   
    double precision :: xc(0:1,0:1,0:1,1:ndim),xp(1:ndim),factor(0:1,0:1,0:1)
    integer :: ix1,ix2,ix3
    double precision :: dxc1,dxc2,dxc3,xd1,xd2,xd3

    dxc1=xc(1,1,1,1)-xc(0,0,0,1)
    dxc2=xc(1,1,1,2)-xc(0,0,0,2)
    dxc3=xc(1,1,1,3)-xc(0,0,0,3)
    xd1=(xp(1)-xc(0,0,0,1))/dxc1
    xd2=(xp(2)-xc(0,0,0,2))/dxc2
    xd3=(xp(3)-xc(0,0,0,3))/dxc3
    ! interpolation factor
    do ix1=0,1
    do ix2=0,1
    do ix3=0,1
      factor(ix1,ix2,ix3)=abs(1-ix1-xd1)*abs(1-ix2-xd2)*abs(1-ix3-xd3)
    enddo
    enddo
    enddo

  end subroutine get_interp_factor_linear

  subroutine interp_linear_multiD(xc,xp,wc,wp,nwc)
    ! get point values from nearby cells via linear interpolation
    ! multi-D method

    integer :: nwc
    double precision :: xc(0:1,0:1,0:1,1:ndim),xp(1:ndim)
    double precision :: wc(0:1,0:1,0:1,1:nwc),wp(1:nwc)

    integer :: ix1,ix2,ix3,iw
    double precision :: factor(0:1,0:1,0:1)

    call get_interp_factor_linear(xc,xp,factor)    

    wp=0.d0
    do iw=1,nwc
      do ix3=0,1
      do ix2=0,1
      do ix1=0,1
        wp(iw)=wp(iw)+wc(ix1,ix2,ix3,iw)*factor(ix1,ix2,ix3)
      enddo
      enddo
      enddo
    enddo

  end subroutine interp_linear_multiD

  subroutine interp_from_grid_linear(igrid,xp,wp)
    ! get point values from given grid via linear interpolation
    ! multi-D method
    
    integer :: igrid
    double precision :: xp(1:ndim),wp(1:nw)

    double precision :: dxb1,dxb2,dxb3,xbmin1,xbmin2,xbmin3,xbmax1,xbmax2,&
       xbmax3
    integer :: inblock,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,j
    integer :: ixb1,ixb2,ixb3,ixi1,ixi2,ixi3
    double precision :: xc(0:1,0:1,0:1,1:ndim),wc(0:1,0:1,0:1,nw)

    ! block/grid boundaries
    xbmin1=rnode(rpxmin1_,igrid)
    xbmin2=rnode(rpxmin2_,igrid)
    xbmin3=rnode(rpxmin3_,igrid)
    xbmax1=rnode(rpxmax1_,igrid)
    xbmax2=rnode(rpxmax2_,igrid)
    xbmax3=rnode(rpxmax3_,igrid)

    ! whether or not next point is in this block/grid
    inblock=0
    if (xp(1)>=xbmin1 .and. xp(1)<xbmax1) inblock=inblock+1
    if (xp(2)>=xbmin2 .and. xp(2)<xbmax2) inblock=inblock+1
    if (xp(3)>=xbmin3 .and. xp(3)<xbmax3) inblock=inblock+1
    if (inblock/=ndim) then
      call mpistop('Interpolation error: given point is not in given grid')
    endif

    ! cell indexes for the point
    dxb1=rnode(rpdx1_,igrid)
    dxb2=rnode(rpdx2_,igrid)
    dxb3=rnode(rpdx3_,igrid)
    ixOmin1=ixmlo1
    ixOmin2=ixmlo2
    ixOmin3=ixmlo3
    ixb1=floor((xp(1)-ps(igrid)%x(ixOmin1,ixOmin2,ixOmin3,1))/dxb1)+ixOmin1
    ixb2=floor((xp(2)-ps(igrid)%x(ixOmin1,ixOmin2,ixOmin3,2))/dxb2)+ixOmin2
    ixb3=floor((xp(3)-ps(igrid)%x(ixOmin1,ixOmin2,ixOmin3,3))/dxb3)+ixOmin3

    ! nearby cells for interpolation
    do ixi1=0,1
    do ixi2=0,1
    do ixi3=0,1
      xc(ixi1,ixi2,ixi3,:)=ps(igrid)%x(ixb1+ixi1,ixb2+ixi2,ixb3+ixi3,:)
      wc(ixi1,ixi2,ixi3,:)=ps(igrid)%w(ixb1+ixi1,ixb2+ixi2,ixb3+ixi3,:)
    enddo
    enddo
    enddo

    call interp_linear_multiD(xc,xp,wc,wp,nw)

  end subroutine interp_from_grid_linear

end module mod_interpolation
