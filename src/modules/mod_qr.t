!*******************************************************
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
module mod_qr
  implicit none

  contains

  ! Constructs the QR decomposition of the n × n matrix a. 
  ! The upper triangular matrix R is returned in the upper 
  ! triangle of a, except for the diagonal elements of R, 
  ! which are returned in the n-dimensional vector d. 
  ! The orthogonal matrix Q is represented as a product of n−1.
  ! Householder matrices Q1 . . . Qn−1, where Qj = 1 − uj cross uj/cj. 
  ! The ith component of uj is zero for i = 1, . . . , j − 1 while 
  ! the nonzero components are returned in a(i,j) for i = j, . . . , n. 
  ! SING returns as true if singularity is encountered during the decomposition,
  ! but the decomposition is still completed in this case.
  subroutine qrdcmp(a,c,d,sing)
     use mod_math, only: outerprod
     implicit none
     double precision, dimension(:,:), intent(inout) :: a
     double precision, dimension(:), intent(out)     :: c,d
     logical, intent(out)                            :: sing

     integer                      :: k, n
     double precision             :: scale, sigma
   
     n = size(a,1)
     sing = .False.
     do k = 1, n-1
        scale = maxval( dabs(a(k:n,k)) )
        if (scale == 0.0d0) then
           sing = .True.
           c(k) = 0.0d0
           d(k) = 0.0d0
        else
           a(k:n,k) = a(k:n,k) / scale
           sigma = sign( dsqrt(dot_product(a(k:n,k),a(k:n,k))), a(k,k) )
           a(k,k) = a(k,k) + sigma
           c(k) = sigma * a(k,k)
           d(k) = - scale * sigma
           a(k:n,k+1:n) = a(k:n,k+1:n) - &
               outerprod(a(k:n,k),matmul(a(k:n,k),a(k:n,k+1:n))) / c(k)
        end if
     end do
     d(n) = a(n,n)
     if (d(n) == 0.0d0) sing = .True.
  end subroutine qrdcmp

  ! Solves the set of n linear equations A · x = b. 
  ! The n × n matrix a and the n-dimensional
  ! vectors c and d are input as the output of the 
  ! routine qrdcmp and are not modified. 
  ! b is input as the right-hand-side vector of 
  ! length n, and is overwritten with the solution vector on output
  subroutine qrsolv(a,c,d,b)
     implicit none
     double precision, dimension(:,:), intent(in) :: a
     double precision, dimension(:), intent(in)   :: c,d
     double precision, dimension(:), intent(out)  :: b
     integer                      :: j, n
     double precision             :: tau
     n = size(a,1)
     ! construct Q^T dot b
     do j = 1, n-1
        tau = dot_product(a(j:n,j), b(j:n)) / c(j)
        b(j:n) = b(j:n) - tau * a(j:n,j)
     end do
     call rsolv(a,d,b)      ! solve R dot x = Q^T dot b
  end subroutine qrsolv

  ! Solves the set of n linear equations R · x = b, where R is an upper triangular matrix stored
  ! in a and d. The n × n matrix a and the vector d of length n are input as the output of the
  ! routine qrdcmp and are not modified. b is input as the right-hand-side vector of length n,
  ! and is overwritten with the solution vector on output.
  subroutine rsolv(a,d,b)
     implicit none
     double precision, dimension(:,:), intent(in) :: a
     double precision, dimension(:), intent(in)   :: d
     double precision, dimension(:), intent(inout):: b
     integer                      :: i, n
     n = size(a,1)
     b(n) = b(n) / d(n)
     do i = n-1, 1, -1
        b(i) = ( b(i) - dot_product(a(i,i+1:n),b(i+1:n)) ) / d(i)
     end do
  end subroutine rsolv

  subroutine qrupdt(r,qt,u,v)
     implicit none
     double precision, dimension(:,:), intent(inout) :: r, qt
     double precision, dimension(:), intent(inout)   :: u
     double precision, dimension(:), intent(in)      :: v

     integer                      :: i, k, n
     n = size(r,1)
     ! Find largest k such that u(k) /= 0
     do k = n, 1, -1
        if ( u(k) /= 0.0d0 ) exit
     end do
     !k = n + 1 - ifirstloc( u(n:1:-1) /= 0.0d0 )
     k = max(k, 1)
     do i = k-1, 1, -1
        call rotate(r,qt,i,u(i),-u(i+1))
        u(i) = dsqrt( u(i)**2 + u(i+1)**2 ) 
        !u(i) = pythag(u(i),u(i+1))
     end do
     r(1,:) = r(1,:) + u(1) * v
     ! Transform upper Hessenberg matrix to upper
     do i = 1, k-1
        call rotate(r,qt,i,r(i,i),-r(i+1,i))
     end do
  end subroutine qrupdt

  subroutine rotate(r,qt,i,a,b)
     implicit none
     double precision, dimension(:,:), target, intent(inout) :: r, qt
     integer, intent(in)                                     :: i
     double precision, intent(in)                            :: a,b

     integer                                :: n
     double precision, dimension(size(r,1)) :: temp
     double precision                       :: c, fact, s
     n = size(r,1)
     if (a == 0.0d0) then
        ! Avoid unnecessary overflow or underflow
        c = 0.0d0
        s = sign(1.0d0, b)
     else if (dabs(a) > dabs(b)) then
        fact = b / a
        c = sign( 1.0d0/dsqrt(1.0d0+fact**2), a)
        s = fact * c
     else
        fact = a / b
        s = sign( 1.0d0/dsqrt(1.0d0+fact**2), b)
        c = fact * s
     end if
     ! Premultiply r by Jacobi rotation
     temp(i:n) = r(i,i:n)
     r(i,i:n)   = c * temp(i:n) - s * r(i+1,i:n)
     r(i+1,i:n) = s * temp(i:n) + c * r(i+1,i:n)
     ! Premultiply qt by Jacobi rotation
     temp = qt(i,:)
     qt(i,:)   = c * temp - s * qt(i+1,:)
     qt(i+1,:) = s * temp + c * qt(i+1,:)
  end subroutine rotate

  double precision function pythag(a,b)
     implicit none
     double precision, intent(in) :: a, b
     double precision             :: absa, absb
     absa = dabs(a)
     absb = dabs(b)
     if (absa > absb) then
        pythag = absa * dsqrt(1.0d0+(absb/absa)**2)
     else
        if (absb == 0.0d0) then
           pythag = 0.0d0
        else
           pythag = absb * dsqrt(1.0d0+(absa/absb)**2)
        end if
     end if
  end function pythag

  integer function ifirstloc(mask)
     implicit none
     logical, dimension(:), intent(in) :: mask
     integer                           :: i
     do i = 1, size(mask)
        if (mask(i)) then
           ifirstloc = i
           return
        end if
     end do
     ifirstloc = i
  end function ifirstloc

END MODULE mod_qr
