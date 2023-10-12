MODULE mod_lnsrch
  implicit none

CONTAINS

 subroutine lnsrch(xold,fold,g,p,x,f,stpmax,err_code,func)
    implicit none
    double precision, dimension(:), intent(in) :: xold, g
    double precision, dimension(:), intent(inout) :: p
    double precision, intent(in) :: fold, stpmax
    double precision, dimension(:), intent(out) :: x
    double precision, intent(out) :: f
    integer, intent(out) :: err_code
    !> err_code: 
    !> 0: everything ok
    !> 1: here is the root 
    !> 2: slope >= 0
    interface
       function func(x)
          implicit none
          double precision :: func
          double precision, dimension(:), intent(in) :: x
       end function func
    end interface
    
    double precision, parameter :: alf=1.0d-4  ! this ensures sufficient decrease in function value, should be within (0,1)
    double precision, parameter :: tolx=epsilon(x)
    integer :: ndum
    double precision :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam

    err_code = 0
    
    pabs = dsqrt( dot_product(p(:),p(:)) )
    if ( pabs > stpmax ) p(:) = p(:) * stpmax / pabs
    ! scale if attempted step is too big

    slope = dot_product(g(:),p(:))
    if ( slope >= 0.0d0 ) then
       err_code = 2
       return
       !stop "error in lnsrch: slope is larger then or equals to zero"
    end if

    alamin = tolx / maxval(dabs(p)/max(dabs(xold),1.0d0)) ! compute lambda_min
    alam = 1.0d0 ! Always try full Newton step first
    do
       x(:) = xold(:) + alam * p(:)
       f = func(x)
       if ( alam < alamin ) then
          ! for root finding problems, this means it is converged
          ! celler should check the convergence
          x(:) = xold(:)
          err_code = 1
          return
       else if ( f <= fold + alf * alam * slope ) then
          ! Sufficient function decrease
          return
       else
          ! backtrack
          if ( alam == 1.0d0 ) then
             tmplam = - slope / ( 2.0d0 * (f-fold-slope) )
          else
             rhs1 = f - fold - alam * slope
             rhs2 = f2 - fold - alam2 * slope
             a = ( rhs1 / alam**2 - rhs2 / alam2**2 ) / ( alam - alam2 )
             b = ( -alam2 * rhs1 / alam**2 + alam * rhs2 / alam2**2 ) &
                   / ( alam - alam2 )
             if ( a == 0.0d0 ) then
                tmplam = - slope / ( 2.0d0 * b )
             else
                disc = b**2 - 3.0d0 * a * slope
                if (disc < 0.0d0) then
                   tmplam = 0.5d0 * alam
                else if (b <= 0.0d0) then
                   tmplam = ( -b + dsqrt(disc) ) / ( 3.0d0 * a )
                else
                   tmplam = - slope / ( b + dsqrt(disc) )
                end if
             end if
             if (tmplam > 0.5d0 * alam) tmplam = 0.5d0 * alam
          end if
       end if
       alam2 = alam
       f2 = f
       alam = max(tmplam, 0.1d0 * alam)
    end do ! try again

 end subroutine lnsrch

END MODULE mod_lnsrch
