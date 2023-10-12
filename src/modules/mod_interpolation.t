module mod_interpolation
  use mod_global_parameters
  implicit none
contains

  subroutine interp_linear(x_table,y_table,n_table,x_itp,y_itp,n_itp)
    ! linear interpolation
    ! 1D method

    integer :: n_table,n_itp
    double precision :: x_table(n_table),y_table(n_table)
    double precision :: x_itp(n_itp),y_itp(n_itp)

    integer :: i,j

    if (x_table(1)>x_itp(1) .or. x_table(n_table)<x_itp(n_itp)) then
      call mpistop("out of interpolation table")
    endif

    do i=1,n_itp
      LOOPT: do j=1,n_table-1
        if (x_itp(i)>=x_table(j) .and. x_itp(i)<x_table(j+1)) then
          y_itp(i)=y_table(j) + (x_itp(i)-x_table(j))*&
                   (y_table(j+1)-y_table(j))/(x_table(j+1)-x_table(j))
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
      matrix2(i+1,i)=matrix1(i+1,i)/(matrix1(i,i)-matrix1(i-1,i)*matrix2(i,i-1))
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

  subroutine lagrange_interpolation(xx,yy,x,y,inverse)
    ! Interpolate xx=[x1,x2,...xn]
    !             yy=[f(x1),f(x2),...,f(xn)]
    ! Using Lagrange polynomial P(x)
    ! Returns y = P(x)
    implicit none
    double precision,dimension(:),intent(in) :: xx, yy
    double precision,intent(in)              :: x
    double precision,intent(out)             :: y
    logical,intent(in),optional              :: inverse
    
    !local variables:
    logical :: inv
    integer :: j,k,n,m, i_start, i_end, step
    double precision :: p

    inv=.false.
    if (present(inverse)) inv=inverse
    
    !check number of points:
    n = size(xx); m = size(yy)
    if ( (n/=m).or.(n<2) ) error stop &
        'Error: vectors must be the same size.'

    if (inv) then
      i_start=n; i_end=1; step = -1
    else
      i_start=1; i_end=n; step = 1
    end if
    
    !sum each of the Pj(x) terms:
    y = 0.0d0
    do j=i_start,i_end,step
        !compute Pj(x):
        p = yy(j)
        do k=i_start,i_end,step
            if (k/=j) p = p * (x-xx(k)) / (xx(j)-xx(k))
        end do
        y = y + p
    end do
   end subroutine lagrange_interpolation


   !---------------------------------------------------------------------
   !
   !     purpose: interpolation of a function of three variables in an
   !              equidistant(!!!) table.
   !
   !     method:  8-point Lagrange linear interpolation formula          
   !
   !     x        input vector of first  variable
   !     y        input vector of second variable
   !     z        input vector of third  variable
   !
   !     f        output vector of interpolated function values
   !
   !
   !     ft       3d array of tabulated function values
   !     nx       x-dimension of table
   !     ny       y-dimension of table
   !     nz       z-dimension of table
   !     xt       vector of x-coordinates of table
   !     yt       vector of y-coordinates of table
   !     zt       vector of z-coordinates of table
   !
   !     d1       centered derivative of ft with respect to x
   !     d2       centered derivative of ft with respect to y
   !     d3       centered derivative of ft with respect to z
   !---------------------------------------------------------------------
   subroutine intep3d_uni(x, y, z, f, &
                   ft, &
                   nx, ny, nz, &
                   xt, yt, zt, &
                   ix_out, iy_out, iz_out, &
                   d1, d2, d3)
      implicit none
      integer, intent(in)                            ::  nx, ny, nz
      double precision, intent(in)                   ::  x, y, z
      double precision, intent(out)                  ::  f
      integer, intent(out), optional                 ::  ix_out, iy_out, iz_out
      double precision, intent(out), optional        ::  d1, d2, d3
      double precision, intent(in)                   ::  xt(nx), yt(ny), zt(nz)
      double precision, intent(in)                   ::  ft(nx, ny, nz)
 
      integer                                        ::  ix, iy, iz
      double precision                               ::  dxi, dyi, dzi
      double precision                               ::  xx, yy, zz
 
      ! determine spacing parameters of (equidistant!!!) table

      ix = max(1, min(int((x - xt(1)) / (xt(2) - xt(1))) + 1, nx-1))
      iy = max(1, min(int((y - yt(1)) / (yt(2) - yt(1))) + 1, ny-1))
      iz = max(1, min(int((z - zt(1)) / (zt(2) - zt(1))) + 1, nz-1))
 
      if(present(ix_out)) ix_out = ix
      if(present(iy_out)) iy_out = iy
      if(present(iz_out)) iz_out = iz
 
      dxi   = 1.d0 / ( xt(ix+1)-xt(ix) )
      dyi   = 1.d0 / ( yt(iy+1)-yt(iy) )
      dzi   = 1.d0 / ( zt(iz+1)-zt(iz) )
 
      xx = ( x - xt(ix) ) * dxi
      yy = ( y - yt(iy) ) * dyi
      zz = ( z - zt(iz) ) * dzi
 
      f = ft(ix  , iy  , iz  ) * (1.d0 - xx) * (1.d0 - yy) * (1.d0 - zz) &
        + ft(ix+1, iy  , iz  ) * xx          * (1.d0 - yy) * (1.d0 - zz) &
        + ft(ix  , iy+1, iz  ) * (1.d0 - xx) * yy          * (1.d0 - zz) &
        + ft(ix  , iy  , iz+1) * (1.d0 - xx) * (1.d0 - yy) * zz          &
        + ft(ix+1, iy+1, iz  ) * xx          * yy          * (1.d0 - zz) &
        + ft(ix+1, iy  , iz+1) * xx          * (1.d0 - yy) * zz          &
        + ft(ix  , iy+1, iz+1) * (1.d0 - xx) * yy          * zz          &
        + ft(ix+1, iy+1, iz+1) * xx          * yy          * zz
 
      if(present(d1)) d1 = &
          ft(ix  , iy  , iz  ) * ( - dxi   ) * (1.d0 - yy) * (1.d0 - zz) &
        + ft(ix+1, iy  , iz  ) * dxi         * (1.d0 - yy) * (1.d0 - zz) &
        + ft(ix  , iy+1, iz  ) * ( - dxi   ) * yy          * (1.d0 - zz) &
        + ft(ix  , iy  , iz+1) * ( - dxi   ) * (1.d0 - yy) * zz          &
        + ft(ix+1, iy+1, iz  ) * dxi         * yy          * (1.d0 - zz) &
        + ft(ix+1, iy  , iz+1) * dxi         * (1.d0 - yy) * zz          &
        + ft(ix  , iy+1, iz+1) * ( - dxi   ) * yy          * zz          &
        + ft(ix+1, iy+1, iz+1) * dxi         * yy          * zz
 
      if(present(d2)) d2 = &
          ft(ix  , iy  , iz  ) * (1.d0 - xx) * ( - dyi   ) * (1.d0 - zz) &
        + ft(ix+1, iy  , iz  ) * xx          * ( - dyi   ) * (1.d0 - zz) &
        + ft(ix  , iy+1, iz  ) * (1.d0 - xx) * dyi         * (1.d0 - zz) &
        + ft(ix  , iy  , iz+1) * (1.d0 - xx) * ( - dyi   ) * zz          &
        + ft(ix+1, iy+1, iz  ) * xx          * dyi         * (1.d0 - zz) &
        + ft(ix+1, iy  , iz+1) * xx          * ( - dyi   ) * zz          &
        + ft(ix  , iy+1, iz+1) * (1.d0 - xx) * dyi         * zz          &
        + ft(ix+1, iy+1, iz+1) * xx          * dyi         * zz
 
      if(present(d3)) d3 = &
          ft(ix  , iy  , iz  ) * (1.d0 - xx) * (1.d0 - yy) * ( - dzi   ) &
        + ft(ix+1, iy  , iz  ) * xx          * (1.d0 - yy) * ( - dzi   ) &
        + ft(ix  , iy+1, iz  ) * (1.d0 - xx) * yy          * ( - dzi   ) &
        + ft(ix  , iy  , iz+1) * (1.d0 - xx) * (1.d0 - yy) * dzi         &
        + ft(ix+1, iy+1, iz  ) * xx          * yy          * ( - dzi   ) &
        + ft(ix+1, iy  , iz+1) * xx          * (1.d0 - yy) * dzi         &
        + ft(ix  , iy+1, iz+1) * (1.d0 - xx) * yy          * dzi         &
        + ft(ix+1, iy+1, iz+1) * xx          * yy          * zz
 
   end subroutine 
 
   subroutine intep3d_uni_many ( x, y, z, f, ft, nx, ny, nz, nvars, xt, yt, zt)
      implicit none
 
      integer, intent(in)            :: nx,ny,nz,nvars
      double precision, intent(in)   :: x, y, z, xt(nx), yt(ny), zt(nz), ft(nx,ny,nz,nvars)
      double precision, intent(out)  :: f(nvars)
 
      double precision               :: dxi,dyi,dzi, xx, yy, zz
      integer                        :: ix,iy,iz
 
      ! determine spacing parameters of (equidistant!!!) table
 
      ix = max(1, min(int((x - xt(1)) / (xt(2) - xt(1))) + 1, nx-1))
      iy = max(1, min(int((y - yt(1)) / (yt(2) - yt(1))) + 1, ny-1))
      iz = max(1, min(int((z - zt(1)) / (zt(2) - zt(1))) + 1, nz-1))
 
      dxi   = 1.d0 / ( xt(ix+1)-xt(ix) )
      dyi   = 1.d0 / ( yt(iy+1)-yt(iy) )
      dzi   = 1.d0 / ( zt(iz+1)-zt(iz) )
 
      xx  = (x - xt(ix)) * dxi
      yy  = (y - yt(iy)) * dyi
      zz  = (z - zt(iz)) * dzi
 
      f(1:nvars) = &
          ft(ix  , iy  , iz  , 1:nvars) * (1.d0 - xx) * (1.d0 - yy) * (1.d0 - zz) &
        + ft(ix+1, iy  , iz  , 1:nvars) * xx          * (1.d0 - yy) * (1.d0 - zz) &
        + ft(ix  , iy+1, iz  , 1:nvars) * (1.d0 - xx) * yy          * (1.d0 - zz) &
        + ft(ix  , iy  , iz+1, 1:nvars) * (1.d0 - xx) * (1.d0 - yy) * zz          &
        + ft(ix+1, iy+1, iz  , 1:nvars) * xx          * yy          * (1.d0 - zz) &
        + ft(ix+1, iy  , iz+1, 1:nvars) * xx          * (1.d0 - yy) * zz          &
        + ft(ix  , iy+1, iz+1, 1:nvars) * (1.d0 - xx) * yy          * zz          &
        + ft(ix+1, iy+1, iz+1, 1:nvars) * xx          * yy          * zz
 
   end subroutine intep3d_uni_many

   subroutine intep2d_uni_many ( x, y, f, ft, nx, ny, nvars, xt, yt)
      implicit none
 
      integer, intent(in)            :: nx,ny,nvars
      double precision, intent(in)   :: x, y, xt(nx), yt(ny), ft(nx,ny,nvars)
      double precision, intent(out)  :: f(nvars)
 
      double precision               :: dxi, dyi, xx, yy
      integer                        :: ix,iy
 
      ! determine spacing parameters of (equidistant!!!) table
 
      ix = max(1, min(int((x - xt(1)) / (xt(2) - xt(1))) + 1, nx-1))
      iy = max(1, min(int((y - yt(1)) / (yt(2) - yt(1))) + 1, ny-1))
 
      dxi   = 1.d0 / ( xt(ix+1)-xt(ix) )
      dyi   = 1.d0 / ( yt(iy+1)-yt(iy) )
 
      xx  = (x - xt(ix)) * dxi
      yy  = (y - yt(iy)) * dyi
 
      f(1:nvars) = &
          ft(ix  , iy  , 1:nvars) * (1.d0 - xx) * (1.d0 - yy) &
        + ft(ix+1, iy  , 1:nvars) * xx          * (1.d0 - yy) &
        + ft(ix  , iy+1, 1:nvars) * (1.d0 - xx) * yy          &
        + ft(ix+1, iy+1, 1:nvars) * xx          * yy
 
   end subroutine intep2d_uni_many

  double precision function intep1D(x_pt, f_array, x_array) result(point_value)
    implicit none
    double precision, intent(in)   :: x_pt 
    double precision, intent(in)   :: f_array(:)
    double precision, intent(in)   :: x_array(:)

    integer            :: ix_min, ix_max, ix, pm1x
    double precision   :: dx, del_x

    ix_min = lbound(x_array,1) ; ix_max = ubound(x_array,1)
  
    ! find the closest points
    ix = minloc( dabs(x_array(ix_min+1:ix_max-1)-x_pt), dim = 1 )
    ! make sure the index falls into correct range
    ix = min( max(ix, ix_min+1), ix_max-1 )
    del_x = x_pt - x_array(ix)
    if (del_x > 0.0d0) then
       pm1x = 1
    else
       pm1x = -1
    end if
   
    dx = ( x_array(ix + pm1x) - x_array(ix) )
    point_value = f_array(ix) &
                  + del_x/dx * (f_array(ix + pm1x) - f_array(ix))
  
  end function intep1D
  
  double precision function intep2D(x_pt, y_pt, f_array, x_array, y_array) result(point_value)
    implicit none
    double precision, intent(in)   :: x_pt, y_pt
    double precision, intent(in)   :: f_array(:,:)
    double precision, intent(in)   :: x_array(:), y_array(:)

    integer            :: ix_min, ix_max, ix, pm1x
    integer            :: iy_min, iy_max, iy, pm1y
    double precision   :: dx, del_x
    double precision   :: dy, del_y
    double precision   :: fh(4), a(4)

    ix_min = lbound(x_array,1) ; ix_max = ubound(x_array,1)
    iy_min = lbound(y_array,1) ; iy_max = ubound(y_array,1)
  
    ! find the closest points
    ix = minloc( dabs(x_array(ix_min+1:ix_max-1)-x_pt), dim = 1 )
    iy = minloc( dabs(y_array(iy_min+1:iy_max-1)-y_pt), dim = 1 )
    ! make sure the index falls into correct range
    ix = min( max(ix, ix_min+1), ix_max-1 )
    iy = min( max(iy, iy_min+1), iy_max-1 )

    del_x = x_pt - x_array(ix)
    del_y = y_pt - y_array(iy)
    if (del_x > 0.0d0) then
       pm1x = 1
    else
       pm1x = -1
    end if
    if (del_y > 0.0d0) then
       pm1y = 1
    else
       pm1y = -1
    end if
   
    dx = ( x_array(ix + pm1x) - x_array(ix) )
    dy = ( y_array(iy + pm1y) - y_array(iy) )

    fh(1) = f_array(ix       , iy       )
    fh(2) = f_array(ix + pm1x, iy       )
    fh(3) = f_array(ix       , iy + pm1y)
    fh(4) = f_array(ix + pm1x, iy + pm1y)

    a(1) = fh(1)
    a(2) = (fh(2)-fh(1))/dx
    a(3) = (fh(3)-fh(1))/dy
    a(4) = (fh(4)-fh(2)-fh(3)+fh(1))/dy/dx
  
    point_value  = a(1) &
                 + a(2) * del_x &
                 + a(3) * del_y &
                 + a(4) * del_y * del_x
  
  end function intep2D
  
  double precision function intep3D(x_pt, y_pt, z_pt, f_array, x_array, y_array, z_array) result(point_value)
    implicit none
    double precision, intent(in)   :: x_pt, y_pt, z_pt
    double precision, intent(in)   :: f_array(:,:,:)
    double precision, intent(in)   :: x_array(:), y_array(:), z_array(:)

    integer            :: ix_min, ix_max, ix, pm1x
    integer            :: iy_min, iy_max, iy, pm1y
    integer            :: iz_min, iz_max, iz, pm1z
    double precision   :: dx, del_x
    double precision   :: dy, del_y
    double precision   :: dz, del_z
    double precision   :: fh(8), a(8)

    ix_min = lbound(x_array,1) ; ix_max = ubound(x_array,1)
    iy_min = lbound(y_array,1) ; iy_max = ubound(y_array,1)
    iz_min = lbound(z_array,1) ; iz_max = ubound(z_array,1)
  
    ! find the closest points
    ix = minloc( dabs(x_array(ix_min+1:ix_max-1)-x_pt), dim = 1 )
    iy = minloc( dabs(y_array(iy_min+1:iy_max-1)-y_pt), dim = 1 )
    iz = minloc( dabs(z_array(iz_min+1:iz_max-1)-z_pt), dim = 1 )
    ! make sure the index falls into correct range
    ix = min( max(ix, ix_min+1), ix_max-1 )
    iy = min( max(iy, iy_min+1), iy_max-1 )
    iz = min( max(iz, iz_min+1), iz_max-1 )

    del_x = x_pt - x_array(ix)
    del_y = y_pt - y_array(iy)
    del_z = z_pt - z_array(iz)
    if (del_x > 0.0d0) then
       pm1x = 1
    else
       pm1x = -1
    end if
    if (del_y > 0.0d0) then
       pm1y = 1
    else
       pm1y = -1
    end if
    if (del_z > 0.0d0) then
       pm1z = 1
    else
       pm1z = -1
    end if
   
    dx = ( x_array(ix + pm1x) - x_array(ix) )
    dy = ( y_array(iy + pm1y) - y_array(iy) )
    dz = ( z_array(iz + pm1z) - z_array(iz) )

    fh(1) = f_array(ix       , iy       , iz       )
    fh(2) = f_array(ix + pm1x, iy       , iz       )
    fh(3) = f_array(ix       , iy + pm1y, iz       )
    fh(4) = f_array(ix       , iy       , iz + pm1z)
    fh(5) = f_array(ix + pm1x, iy + pm1y, iz       )
    fh(6) = f_array(ix + pm1x, iy       , iz + pm1z)
    fh(7) = f_array(ix       , iy + pm1y, iz + pm1z)
    fh(8) = f_array(ix + pm1x, iy + pm1y, iz + pm1z)
  
    a(1) = fh(1)
    a(2) = (fh(2)-fh(1))/dx
    a(3) = (fh(3)-fh(1))/dy
    a(4) = (fh(4)-fh(1))/dz
    a(5) = (fh(5)-fh(2)-fh(3)+fh(1))/dx/dy
    a(6) = (fh(6)-fh(2)-fh(4)+fh(1))/dx/dz
    a(7) = (fh(7)-fh(3)-fh(4)+fh(1))/dy/dz
    a(8) = (fh(8)-fh(1)+fh(2)+fh(3) &
           +fh(4)-fh(5)-fh(6)-fh(7))/dx/dy/dz
  
    point_value = a(1) &
                + a(2) * del_x         &
                + a(3) * del_y         &
                + a(4) * del_z         &
                + a(5) * del_x * del_y &
                + a(6) * del_x * del_z &
                + a(7) * del_y * del_z &
                + a(8) * del_x * del_y * del_z
  
  end function intep3D

end module mod_interpolation
