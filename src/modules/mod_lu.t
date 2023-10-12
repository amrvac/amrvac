!*******************************************************
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
module mod_lu
  implicit none

  contains

  !***************************************************************
  !* Given an N x N matrix A, this routine replaces it by the LU *
  !* decomposition of a rowwise permutation of itself. A and N   *
  !* are input. INDX is an output vector which records the row   *
  !* permutation effected by the partial pivoting; D is output   *
  !* as -1 or 1, depending on whether the number of row inter-   *
  !* changes was even or odd, respectively. This routine is used *
  !* in combination with LUBKSB to solve linear equations or to  *
  !* invert a matrix. Return code is 1, if matrix is singular.   *
  !***************************************************************
  subroutine ludcmp(n_dim,A,indx,d,error_code)
     implicit none
     integer, intent(in)                                     :: n_dim
     double precision, dimension(n_dim,n_dim), intent(inout) :: A
     integer, dimension(n_dim), intent(out)                  :: indx
     integer, intent(out)                                    :: d
     integer, intent(out)                                    :: error_code

     double precision, parameter  :: TINY=1.0D-20 ! a very small number
     integer                      :: j, imax
     double precision             :: vv(n_dim), tmp(n_dim)
   
     error_code = 0
     d = 1 ! No row interchanges yet
     vv = maxval(dabs(A),dim=2)  ! Loop over rows to get the implicit scaling information
     if (any(vv==0.0d0)) then
        ! there is a row of zeros
        ! singular matrix
        error_code = 1
        return
     end if
     vv = 1.0d0 / vv   ! save the scaling
     do j=1, n_dim
        ! Find the pivot row
        imax = (j-1) + sum( maxloc(vv(j:n_dim) * dabs(A(j:n_dim,j))) )
        if (j /= imax) then
           ! if we need to interchange rows
           tmp(:) = a(imax,:)
           a(imax,:) = a(j,:)
           a(j,:) = tmp(:)
        end if
        indx(j) = imax
        ! if the pivot element is zero the matrix is singular
        if ( dabs(A(j,j)) <= TINY ) A(j,j) = TINY
        ! Divide by the pivot element
        A(j+1:n_dim,j) = A(j+1:n_dim,j) / A(j,j)
        A(j+1:n_dim,j+1:n_dim) = A(j+1:n_dim,j+1:n_dim) &
                - spread( A(j+1:n_dim,j), dim=2, ncopies=n_dim-j) &
                 *spread( A(j,j+1:n_dim), dim=1, ncopies=n_dim-j)
                !- outerprod( A(j+1:n_dim,j), A(j,j+1:n_dim))
     end do
  end subroutine ludcmp


  !******************************************************************
  !* Solves the set of N linear equations A . X = B.  Here A is     *
  !* input, not as the matrix A but rather as its LU decomposition, *
  !* determined by the routine LUDCMP. INDX is input as the permuta-*
  !* tion vector returned by LUDCMP. B is input as the right-hand   *
  !* side vector B, and returns with the solution vector X. A, N and*
  !* INDX are not modified by this routine and can be used for suc- *
  !* cessive calls with different right-hand sides. This routine is *
  !* also efficient for plain matrix inversion.                     *
  !******************************************************************
  subroutine lubksb(n_dim,A,indx,b)
     implicit none
     integer, intent(in)                                     :: n_dim
     double precision, dimension(n_dim,n_dim), intent(in)    :: A
     integer, dimension(n_dim), intent(in)                   :: indx
     double precision, dimension(n_dim), intent(inout)       :: b

     integer                      :: i, ii, ll
     double precision             :: summ

     ii = 0
     do i = 1, n_dim
        ll = indx(i)
        summ = b(ll)
        b(ll) = b(i)
        if (ii /= 0) then
           summ = summ - dot_product(A(i,ii:i-1),b(ii:i-1))
        else if (summ /= 0.0d0) then
           ii = i
        end if
        b(i) = summ
     end do
     do i = n_dim, 1, -1
        b(i) = ( b(i) - dot_product( A(i,i+1:n_dim), b(i+1:n_dim)) ) / A(i,i)
     end do
  end subroutine lubksb
END MODULE mod_lu
