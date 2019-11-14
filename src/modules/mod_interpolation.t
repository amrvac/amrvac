module mod_interpolation
  use mod_global_parameters
  implicit none

contains

  subroutine interp_linear(x_table,y_table,n_table,x_itp,y_itp,n_itp)

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









end module mod_interpolation
