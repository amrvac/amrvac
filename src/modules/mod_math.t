MODULE mod_math
  implicit none

  contains

  ! Computes forward-difference approximation to Jacobian. 
  ! On input, x is the point at which the Jacobian is to be evaluated, 
  ! and fvec is the vector of function values at the point,
  ! both arrays of length N. df is the N × N output Jacobian. 
  ! Parameter: EPS is the approximate square root of the machine precision.
  subroutine fdjac(x, fvec_func, df, fvec_in, h_scale)
    !> input as initial guess, output as root
    double precision, dimension(:),   intent(in)  :: x
    double precision, dimension(:,:), intent(out) :: df
    double precision, dimension(:),   intent(in), optional :: fvec_in
    double precision, dimension(:),   intent(in), optional :: h_scale
    interface
       subroutine fvec_func(z_vec, fvec_out)
          implicit none
          double precision, dimension(:), intent(in) :: z_vec
          double precision, dimension(1:size(z_vec)) :: fvec_out
       end subroutine fvec_func
    end interface
    double precision, parameter          :: EPS = dsqrt(epsilon(1.0d0)) ! roughtly ~sqrt(double precision)
    integer                              :: j
    double precision, dimension(size(x)) :: x_tmp, xph, h, fvec_tmp, fvec
    if (present(fvec_in)) then
       fvec = fvec_in
    else
       call fvec_func(x, fvec)
    end if
    h = EPS * dabs(x)
    !where (h == 0.0d0) h = EPS !* dsqrt( dot_product(x,x) )
    do j=1,size(x)
       if (h(j) == 0.0d0) then
          h(j) = EPS
          if (present(h_scale)) h(j) = h(j) * h_scale(j)
       end if
    end do
    xph = x + h
    h = xph - x  ! Trick to reduce finite precision error.
    x_tmp = x
    do j = 1, size(x)
       x_tmp(j) = xph(j)
       call fvec_func(x_tmp, fvec_tmp)
       df(:,j) = ( fvec_tmp(:) - fvec(:) ) / h(j)
       x_tmp(j) = x(j)
    end do
  end subroutine fdjac

  function outerprod(a,b)
    double precision, dimension(:), intent(in)   :: a, b
    double precision, dimension(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
                spread(b,dim=1,ncopies=size(a))
  end function outerprod

END MODULE mod_math
