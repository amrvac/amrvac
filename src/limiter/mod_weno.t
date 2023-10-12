module mod_weno
  ! All kinds of (W)ENO schemes
  !
  ! 2019.9.19 WENO(-JS)5 transplant from the BHAC code;
  ! 2019.9.20 WENO3;
  ! 2019.9.21 WENO-Z5;
  ! 2019.9.22 WENO-Z+5 transplant from the BHAC code;
  ! 2019.10.30 WENO(-JS)7;
  ! 2019.10.31 MPWENO7;
  ! 2019.11.1 exENO7;
  ! 2019.11.7 clean up the code, comment out the interpolation variation;
  ! 2019.12.9 WENO-YC3;
  ! 2020.1.15 new WENO5 limiter WENO5NM for stretched grid.
  ! 2020.4.15 WENO5-CU6: hybrid sixth-order linear & WENO5
   
  implicit none
  private

  double precision, parameter     :: weno_eps_machine = 1.0d-18
  double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0

  public :: WENO3limiter
  public :: WENO5limiter
  public :: WENO5limiterL
  public :: WENO5limiterR
  public :: WENO5NMlimiter
  public :: WENO5NMlimiterL
  public :: WENO5NMlimiterR
  public :: TENO5ADlimiter
  public :: WENO5CU6limiter
  public :: WENO7limiter
  public :: exENO7limiter

contains

  subroutine WENO3limiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC,var,rec)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims, var
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S)
    double precision, intent(inout) :: wRC(ixI^S),wLC(ixI^S) 
    logical, intent(in)             :: rec
    !> local
    integer                         :: iLm^L, iLp^L, iLpp^L
    double precision                :: f_array(ixI^S,2), d_array(2)
    double precision                :: beta(ixI^S,2),tau(ixI^S)
    double precision                :: u1_coeff(2), u2_coeff(2)
    double precision                :: alpha_array(ixI^S,2), alpha_sum(ixI^S), flux(ixI^S)
    integer                         :: i

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  
    iLm^L=iL^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    if (rec) then
      d_array(1:2) = (/ 1.0d0/3.0d0, 2.0d0/3.0d0 /)
    else
      d_array(1:2) = (/ 1.0d0/4.0d0, 3.0d0/4.0d0 /)
    end if
    u1_coeff(1:2) = (/ -1.d0/2.d0, 3.d0/2.d0 /)
    u2_coeff(1:2) = (/ 1.d0/2.d0, 1.d0/2.d0 /)
    
    !> left side
    f_array(iL^S,1) = u1_coeff(1) * w(iLm^S) + u1_coeff(2) * w(iL^S)
    f_array(iL^S,2) = u2_coeff(1) * w(iL^S)  + u2_coeff(2) * w(iLp^S)
  
    beta(iL^S,1) = (w(iL^S) - w(iLm^S))**2
    beta(iL^S,2) = (w(iLp^S) - w(iL^S))**2
  
    alpha_sum(iL^S) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,2
        alpha_array(iL^S,i) = d_array(i)/(beta(iL^S,i) + weno_eps_machine)**2
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(2)
      tau(iL^S) = abs(beta(iL^S,2) - beta(iL^S,1))
      do i = 1,2
        alpha_array(iL^S,i) = d_array(i) * (1.d0 + (tau(iL^S) / &
                                      (beta(iL^S,i) + dxdim**2)))
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    end select

    flux(iL^S) = 0.0d0
    do i = 1,2
      flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
    end do
  
    !> left value at right interface
    wLC(iL^S) = flux(iL^S)
  
    !> right side
    f_array(iL^S,1) = u1_coeff(1) * w(iLpp^S) + u1_coeff(2) * w(iLp^S)
    f_array(iL^S,2) = u2_coeff(1) * w(iLp^S)  + u2_coeff(2) * w(iL^S)
  
    beta(iL^S,1) = (w(iLpp^S) - w(iLp^S))**2
    beta(iL^S,2) = (w(iLp^S) - w(iL^S))**2
  
    alpha_sum(iL^S) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,2
       alpha_array(iL^S,i) = d_array(i)/(beta(iL^S,i) + weno_eps_machine)**2
       alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(2)
      tau(iL^S) = abs(beta(iL^S,2) - beta(iL^S,1))
      do i = 1,2
        alpha_array(iL^S,i) = d_array(i) * (1.d0 + (tau(iL^S) / &
                                      (beta(iL^S,i) + dxdim**2)))
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    end select

    flux(iL^S) = 0.0d0
    do i = 1,2
       flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
    end do
  
    !> right value at right interface
    wRC(iL^S) = flux(iL^S)

  end subroutine WENO3limiter

  subroutine WENO5limiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC,var,rec)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims, var
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S)
    double precision, intent(inout) :: wRC(ixI^S),wLC(ixI^S) 
    logical, intent(in)             :: rec
    !> local
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    double precision                :: f_array(ixI^S,3), d_array(3)
    double precision                :: beta(ixI^S,3), beta_coeff(2)
    double precision                :: tau(ixI^S), tmp(ixI^S)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,3), alpha_sum(ixI^S), flux(ixI^S)
    integer                         :: i
    double precision                :: lambda
    double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    lambda = dxdim**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if (rec) then
      !reconstruction variation
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
      u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
      u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    else
      !interpolation variation
      d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
      u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
      u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
      u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    end if
    
    !> left side
    f_array(iL^S,1) = u1_coeff(1) * w(iLmm^S) + u1_coeff(2) * w(iLm^S) + u1_coeff(3) * w(iL^S)
    f_array(iL^S,2) = u2_coeff(1) * w(iLm^S)  + u2_coeff(2) * w(iL^S)  + u2_coeff(3) * w(iLp^S)
    f_array(iL^S,3) = u3_coeff(1) * w(iL^S)   + u3_coeff(2) * w(iLp^S) + u3_coeff(3) * w(iLpp^S)  
  
    beta(iL^S,1) = beta_coeff(1) * (w(iLmm^S) + w(iL^S) - 2.0d0*w(iLm^S))**2 &
         + beta_coeff(2) * (w(iLmm^S) - 4.0d0 * w(iLm^S) + 3.0d0*w(iL^S))**2
    beta(iL^S,2) = beta_coeff(1) * (w(iLm^S) + w(iLp^S) - 2.0d0 * w(iL^S))**2 &
         + beta_coeff(2) * (w(iLm^S) - w(iLp^S))**2
    beta(iL^S,3) = beta_coeff(1) * (w(iL^S) + w(iLpp^S) - 2.0d0 * w(iLp^S))**2 &
         + beta_coeff(2) * (3.0d0 * w(iL^S) - 4.0d0 * w(iLp^S) + w(iLpp^S))**2
 
    alpha_sum(iL^S) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
         alpha_array(iL^S,i) = d_array(i)/(beta(iL^S,i) + weno_eps_machine)**2
         alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(2)
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3))
      do i = 1,3
        alpha_array(iL^S,i) = d_array(i) * (1.d0 + (tau(iL^S) / &
                                      (beta(iL^S,i) + weno_eps_machine))**2)
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(3)
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3))
      do i = 1,3
        tmp(iL^S) = (tau(iL^S) + weno_eps_machine) / (beta(iL^S,i) + weno_eps_machine)
        alpha_array(iL^S,i) = d_array(i) * (1.0d0 + tmp(iL^S)**2 + lambda/tmp(iL^S))
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    end select

    flux(iL^S) = 0.0d0
    do i = 1,3
      flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
    end do
  
    !> left value at right interface
    wLC(iL^S) = flux(iL^S)
  
    !> right side
    f_array(iL^S,1) = u1_coeff(1) * w(iLppp^S) + u1_coeff(2) * w(iLpp^S) + u1_coeff(3) * w(iLp^S)
    f_array(iL^S,2) = u2_coeff(1) * w(iLpp^S)  + u2_coeff(2) * w(iLp^S)  + u2_coeff(3) * w(iL^S)
    f_array(iL^S,3) = u3_coeff(1) * w(iLp^S)   + u3_coeff(2) * w(iL^S)   + u3_coeff(3) * w(iLm^S)  
  
    beta(iL^S,1) = beta_coeff(1) * (w(iLppp^S) + w(iLp^S) - 2.0d0*w(iLpp^S))**2 &
         + beta_coeff(2) * (w(iLppp^S) - 4.0d0 * w(iLpp^S) + 3.0d0*w(iLp^S))**2
    beta(iL^S,2) = beta_coeff(1) * (w(iLpp^S) + w(iL^S) - 2.0d0 * w(iLp^S))**2 &
         + beta_coeff(2) * (w(iLpp^S) - w(iL^S))**2
    beta(iL^S,3) = beta_coeff(1) * (w(iLp^S) + w(iLm^S) - 2.0d0 * w(iL^S))**2 &
         + beta_coeff(2) * (3.0d0 * w(iLp^S) - 4.0d0 * w(iL^S) + w(iLm^S))**2
  
    alpha_sum(iL^S) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iL^S,i) = d_array(i)/(beta(iL^S,i) + weno_eps_machine)**2
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(2) 
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3))
      do i = 1,3
        alpha_array(iL^S,i) = d_array(i) * (1.d0 + (tau(iL^S) / &
                                      (beta(iL^S,i) + weno_eps_machine))**2)
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(3)
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3))
      do i = 1,3
        tmp(iL^S) = (tau(iL^S) + weno_eps_machine) / (beta(iL^S,i) + weno_eps_machine)
        alpha_array(iL^S,i) = d_array(i) * (1.0d0 + tmp(iL^S)**2 + lambda/tmp(iL^S))
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    end select
    flux(iL^S) = 0.0d0
    do i = 1,3
      flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
    end do
  
    !> right value at right interface
    wRC(iL^S) = flux(iL^S)

  end subroutine WENO5limiter

  subroutine TENO5ADlimiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC,rec)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S)
    double precision, intent(inout) :: wRC(ixI^S),wLC(ixI^S) 
    logical, intent(in)             :: rec
    !> local
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    double precision                :: f_array(ixI^S,3), d_array(3)
    double precision                :: beta(ixI^S,3), beta_coeff(2)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: gam_sum(ixI^S),tau(ixI^S),delta_sum(ixI^S)
    double precision                :: gam(ixI^S,3), kai(ixI^S,3), delta(ixI^S,3)
    double precision                :: flux(ixI^S), kai1(ixI^S,3), theta(ixI^S)
    integer                         :: marray(ixI^S)
    integer                         :: i

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if (rec) then
      !reconstruction variation
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
      u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
      u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    else
      !interpolation variation
      d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
      u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
      u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
      u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    end if
    
    !> left side
    f_array(iL^S,1) = u1_coeff(1) * w(iLmm^S) + u1_coeff(2) * w(iLm^S) + u1_coeff(3) * w(iL^S)
    f_array(iL^S,2) = u2_coeff(1) * w(iLm^S)  + u2_coeff(2) * w(iL^S)  + u2_coeff(3) * w(iLp^S)
    f_array(iL^S,3) = u3_coeff(1) * w(iL^S)   + u3_coeff(2) * w(iLp^S) + u3_coeff(3) * w(iLpp^S)  
  
    beta(iL^S,1) = beta_coeff(1) * (w(iLmm^S) + w(iL^S) - 2.0d0*w(iLm^S))**2 &
         + beta_coeff(2) * (w(iLmm^S) - 4.0d0 * w(iLm^S) + 3.0d0*w(iL^S))**2
    beta(iL^S,2) = beta_coeff(1) * (w(iLm^S) + w(iLp^S) - 2.0d0 * w(iL^S))**2 &
         + beta_coeff(2) * (w(iLm^S) - w(iLp^S))**2
    beta(iL^S,3) = beta_coeff(1) * (w(iL^S) + w(iLpp^S) - 2.0d0 * w(iLp^S))**2 &
         + beta_coeff(2) * (3.0d0 * w(iL^S) - 4.0d0 * w(iLp^S) + w(iLpp^S))**2
 
    gam_sum(iL^S) = 0.0d0
    tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3))
    do i = 1,3
      kai1(iL^S,i) = (tau(iL^S) / (beta(iL^S,i) + weno_eps_machine))
      gam(iL^S,i) = (1.d0 + kai1(iL^S,i))**2
      gam_sum(iL^S) = gam_sum(iL^S) + gam(iL^S,i)
    end do
    theta(iL^S) = one / (one + maxval(kai1(iL^S,1:3)/10.d0)) ! fixme: why this?
    marray(iL^S)=-floor(4.d0 + theta(iL^S)*6.d0)
    do i = 1,3    
      kai(iL^S,i) = gam(iL^S,i) / gam_sum(iL^S)
      where(kai(iL^S,i) .lt. 10**marray(iL^S))
        delta(iL^S,i)=zero
      else where
        delta(iL^S,i)=one
      end where
    end do
    delta_sum=zero
    do i = 1,3
      delta_sum(iL^S)=delta_sum(iL^S)+delta(iL^S,i)*d_array(i)
    end do
    flux(iL^S)=0.0d0
    do i = 1,3
      flux(iL^S)=flux(iL^S)+f_array(iL^S,i)*delta(iL^S,i)*d_array(i)/(delta_sum(iL^S))
    end do
  
    !> left value at right interface
    wLC(iL^S) = flux(iL^S)
  
    !> right side
    f_array(iL^S,1) = u1_coeff(1) * w(iLppp^S) + u1_coeff(2) * w(iLpp^S) + u1_coeff(3) * w(iLp^S)
    f_array(iL^S,2) = u2_coeff(1) * w(iLpp^S)  + u2_coeff(2) * w(iLp^S)  + u2_coeff(3) * w(iL^S)
    f_array(iL^S,3) = u3_coeff(1) * w(iLp^S)   + u3_coeff(2) * w(iL^S)   + u3_coeff(3) * w(iLm^S)  
  
    beta(iL^S,1) = beta_coeff(1) * (w(iLppp^S) + w(iLp^S) - 2.0d0*w(iLpp^S))**2 &
         + beta_coeff(2) * (w(iLppp^S) - 4.0d0 * w(iLpp^S) + 3.0d0*w(iLp^S))**2
    beta(iL^S,2) = beta_coeff(1) * (w(iLpp^S) + w(iL^S) - 2.0d0 * w(iLp^S))**2 &
         + beta_coeff(2) * (w(iLpp^S) - w(iL^S))**2
    beta(iL^S,3) = beta_coeff(1) * (w(iLp^S) + w(iLm^S) - 2.0d0 * w(iL^S))**2 &
         + beta_coeff(2) * (3.0d0 * w(iLp^S) - 4.0d0 * w(iL^S) + w(iLm^S))**2
 

    gam_sum(iL^S)=0.0d0
    tau(iL^S)=abs(beta(iL^S,1)-beta(iL^S,3))
    do i=1,3
      kai1(iL^S,i)=(tau(iL^S)/(beta(iL^S,i)+weno_eps_machine))
      gam(iL^S,i)=(1.d0+kai1(iL^S,i))**6
      gam_sum(iL^S)=gam_sum(iL^S)+gam(iL^S,i)
    end do
    theta(iL^S)=one/(one+maxval(kai1(iL^S,1:3)/10.d0)) !fixme: why this?
    marray(iL^S)=-floor(4.d0+theta(iL^S)*6.d0)
    do i=1,3    
      kai(iL^S,i) = gam(iL^S,i)/gam_sum(iL^S)
      where(kai(iL^S,i) .lt. 10**marray(iL^S))
        delta(iL^S,i)=zero
      else where
        delta(iL^S,i)=one
      end where
    end do
    delta_sum=zero
    do i = 1,3
      delta_sum(iL^S)=delta_sum(iL^S)+delta(iL^S,i)*d_array(i)
    end do
    flux(iL^S)=0.0d0
    do i = 1,3
      flux(iL^S)=flux(iL^S)+f_array(iL^S,i)*delta(iL^S,i)*d_array(i)/(delta_sum(iL^S))
    end do

    !> right value at right interface
    wRC(iL^S)=flux(iL^S)
  end subroutine TENO5ADlimiter

  subroutine WENO5NMlimiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC,var, rec)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L,iL^L,idims,var
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S)
    double precision, intent(inout) :: wRC(ixI^S),wLC(ixI^S) 
    logical, intent(in)             :: rec
    !> local
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    integer                         :: iM^L, iMp^L
    double precision                :: f_array(ixI^S,3), d_array(3)
    double precision                :: beta(ixI^S,3), beta_coeff(2)
    double precision                :: tau(ixI^S), tmp(ixI^S)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,3), alpha_sum(ixI^S), flux(ixI^S)
    double precision                :: wc(ixI^S), wd(ixI^S), we(ixI^S)
    integer                         :: i,j
    double precision                :: lambda(ixI^S)

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    iMmin^D=iLmmin^D;
    iMmax^D=iLpmax^D;
    iMp^L=iM^L+kr(idims,^D);
    lambda(iL^S)=block%mesh%dx(iL^S,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if (rec) then
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ -2.d0/3.d0, -1.d0/3.d0, 2.d0 /)
      u2_coeff(1:3) = (/ -1.d0/3.d0, 2.d0/3.d0, 2.d0/3.d0 /)
      u3_coeff(1:3) = (/ 2.d0/3.d0, 2.d0/3.d0, -1.d0/3.d0 /)
    else
      error stop "not implemented"
    end if
    wc(iM^S) = (block%mesh%dx(iMp^S,idims) * w(iM^S) + block%mesh%dx(iM^S,idims) *  w(iMp^S)) / &
                 (block%mesh%dx(iMp^S,idims) + block%mesh%dx(iM^S,idims))
    wd(iL^S) = ((2.d0 * block%mesh%dx(iLm^S,idims) + block%mesh%dx(iLmm^S,idims)) * w(iLm^S) - block%mesh%dx(iLm^S,idims) * w(iLmm^S)) / &
                 (block%mesh%dx(iLmm^S,idims) + block%mesh%dx(iLm^S,idims))
    we(iL^S) = ((2.d0 * block%mesh%dx(iLpp^S,idims) + block%mesh%dx(iLppp^S,idims)) * w(iLpp^S) - block%mesh%dx(iLpp^S,idims) * w(iLppp^S)) / &
                 (block%mesh%dx(iLppp^S,idims) + block%mesh%dx(iLpp^S,idims))
    !> left side
    f_array(iL^S,1) = u1_coeff(1) * wd(iL^S)   + u1_coeff(2) * wc(iLm^S)+ u1_coeff(3) * w(iL^S)
    f_array(iL^S,2) = u2_coeff(1) * wc(iLm^S)  + u2_coeff(2) * w(iL^S)  + u2_coeff(3) * wc(iL^S)
    f_array(iL^S,3) = u3_coeff(1) * wc(iL^S)   + u3_coeff(2) * w(iLp^S) + u3_coeff(3) * wc(iLp^S)  
  
    beta(iL^S,1) = beta_coeff(1) * (wc(iLm^S) - wd(iL^S))**2 &
         + beta_coeff(2) * (2.d0 * w(iL^S) - wc(iLm^S) - wd(iL^S))**2
    beta(iL^S,2) = beta_coeff(1) * (wc(iLm^S) + wc(iL^S) - 2.0d0 * w(iL^S))**2 &
         + beta_coeff(2) * (wc(iLm^S) - wc(iL^S))**2
    beta(iL^S,3) = beta_coeff(1) * (wc(iL^S) + wc(iLp^S) - 2.0d0 * w(iLp^S))**2 &
         + beta_coeff(2) * (3.0d0 * wc(iL^S) - 4.0d0 * w(iLp^S) + wc(iLp^S))**2
    alpha_sum(iL^S) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
         alpha_array(iL^S,i) = d_array(i)/(4.d0 * beta(iL^S,i) + weno_eps_machine)**2
         alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(2)
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3)) * 4.d0
      do i = 1,3
        alpha_array(iL^S,i) = d_array(i) * (1.d0 + (tau(iL^S) / &
                                      (4.d0 * beta(iL^S,i) + weno_eps_machine))**2)
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(3)
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3)) * 4.d0
      do i=1,3
        tmp(iL^S) = (tau(iL^S) + weno_eps_machine) / (4.d0 * beta(iL^S,i) + weno_eps_machine)
        alpha_array(iL^S,i) = d_array(i) * (1.0d0 + tmp(iL^S)**2 + lambda(iL^S)/tmp(iL^S))
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    end select
    flux(iL^S) = 0.0d0
    do i = 1,3
      flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
    end do
    !> left value at right interface
    wLC(iL^S) = flux(iL^S)
    !> right side
    f_array(iL^S,1) = u1_coeff(1) * we(iL^S)  + u1_coeff(2) * wc(iLp^S) + u1_coeff(3) * w(iLp^S)
    f_array(iL^S,2) = u2_coeff(1) * wc(iLp^S) + u2_coeff(2) * w(iLp^S) + u2_coeff(3) * wc(iL^S)
    f_array(iL^S,3) = u3_coeff(1) * wc(iL^S)  + u3_coeff(2) * w(iL^S)  + u3_coeff(3) * wc(iLm^S)  
    beta(iL^S,1) = beta_coeff(1) * (wc(iLp^S) - we(iL^S))**2 &
         + beta_coeff(2) * (2.d0 * w(iLp^S) - wc(iLp^S) - we(iL^S))**2
    beta(iL^S,2) = beta_coeff(1) * (wc(iLp^S) + wc(iL^S) - 2.0d0 * w(iLp^S))**2 &
         + beta_coeff(2) * (wc(iLp^S) - wc(iL^S))**2
    beta(iL^S,3) = beta_coeff(1) * (wc(iL^S) + wc(iLm^S) - 2.0d0 * w(iL^S))**2 &
         + beta_coeff(2) * (3.0d0 * wc(iL^S) - 4.0d0 * w(iL^S) + wc(iLm^S))**2
    alpha_sum(iL^S) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iL^S,i) = d_array(i)/(4.d0 * beta(iL^S,i) + weno_eps_machine)**2
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(2) 
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3)) * 4.d0
      do i = 1,3
        alpha_array(iL^S,i) = d_array(i) * (1.d0 + (tau(iL^S) / &
                                      (4.d0 * beta(iL^S,i) + weno_eps_machine))**2)
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    case(3)
      tau(iL^S) = abs(beta(iL^S,1) - beta(iL^S,3)) * 4.d0
      do i = 1,3
        tmp(iL^S) = (tau(iL^S) + weno_eps_machine) / (4.d0 * beta(iL^S,i) + weno_eps_machine)
        alpha_array(iL^S,i) = d_array(i) * (1.0d0 + tmp(iL^S)**2 + lambda(iL^S)/tmp(iL^S))
        alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
      end do
    end select
    flux(iL^S) = 0.0d0
    do i = 1,3
      flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
    end do
    !> right value at right interface
    wRC(iL^S) = flux(iL^S)

  end subroutine WENO5NMlimiter

  subroutine WENO5limiterL(rec_from,rec_to,n_w,ixI^L,iL^L,idims,w,wLC,var,rec)
    use mod_global_parameters
  
    integer, intent(in)             :: rec_from, rec_to, n_w
    integer, intent(in)             :: ixI^L, iL^L, idims
    integer, intent(in)             :: var
    double precision, intent(in)    :: w(ixI^S,1:n_w)
    double precision, intent(inout) :: wLC(ixI^S,1:n_w) 
    logical, intent(in)             :: rec
    !> local
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    double precision                :: f_array(ixI^S,1:n_w,3), d_array(3)
    double precision                :: beta(ixI^S,1:n_w,3), beta_coeff(2)
    double precision                :: tau(ixI^S,1:n_w), tmp(ixI^S,1:n_w)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:n_w,3), alpha_sum(ixI^S,1:n_w), flux(ixI^S,1:n_w)
    integer                         :: i,j
    double precision                :: lambda(ixI^S)

    double precision, dimension(ixI^S,1:n_w)  :: wLCtmp

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    lambda(iL^S)=block%mesh%dx(iL^S,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if (rec) then
      !reconstruction variation
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
      u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
      u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    else
      !interpolation variation
      d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
      u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
      u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
      u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    end if
    
    !> left side
    f_array(iL^S,rec_from:rec_to,1) = u1_coeff(1) * w(iLmm^S,rec_from:rec_to) + u1_coeff(2) * w(iLm^S,rec_from:rec_to) + u1_coeff(3) * w(iL^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,2) = u2_coeff(1) * w(iLm^S,rec_from:rec_to)  + u2_coeff(2) * w(iL^S,rec_from:rec_to)  + u2_coeff(3) * w(iLp^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,3) = u3_coeff(1) * w(iL^S,rec_from:rec_to)   + u3_coeff(2) * w(iLp^S,rec_from:rec_to) + u3_coeff(3) * w(iLpp^S,rec_from:rec_to)  
  
    beta(iL^S,rec_from:rec_to,1) = beta_coeff(1) * (w(iLmm^S,rec_from:rec_to) + w(iL^S,rec_from:rec_to) - 2.0d0*w(iLm^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iLmm^S,rec_from:rec_to) - 4.0d0 * w(iLm^S,rec_from:rec_to) + 3.0d0*w(iL^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,2) = beta_coeff(1) * (w(iLm^S,rec_from:rec_to) + w(iLp^S,rec_from:rec_to) - 2.0d0 * w(iL^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iLm^S,rec_from:rec_to) - w(iLp^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,3) = beta_coeff(1) * (w(iL^S,rec_from:rec_to) + w(iLpp^S,rec_from:rec_to) - 2.0d0 * w(iLp^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (3.0d0 * w(iL^S, rec_from:rec_to) - 4.0d0 * w(iLp^S,rec_from:rec_to) + w(iLpp^S,rec_from:rec_to))**2
 
    alpha_sum(iL^S,rec_from:rec_to) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz
    case(1)
      do i = 1,3
         alpha_array(iL^S,rec_from:rec_to,i) = d_array(i)/(beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)**2
         alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(2)
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3))
      do i = 1,3
        alpha_array(iL^S,rec_from:rec_to,i) = d_array(i) * (1.d0 + (tau(iL^S,rec_from:rec_to) / &
                                      (beta(iL^S,rec_from:rec_to,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(3)
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3)) * 4.d0
      do i=1,3
        do j=rec_from,rec_to
          tmp(iL^S,j) = (tau(iL^S,j) + weno_eps_machine) / (4.d0 * beta(iL^S,j,i) + weno_eps_machine)
          alpha_array(iL^S,j,i) = d_array(i) * (1.0d0 + tmp(iL^S,j)**2 + lambda(iL^S)/tmp(iL^S,j))
          alpha_sum(iL^S,j) = alpha_sum(iL^S,j) + alpha_array(iL^S,j,i)
        end do
      end do
    end select
    flux(iL^S,rec_from:rec_to) = 0.0d0
    do i = 1,3
      flux(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to) + f_array(iL^S,rec_from:rec_to,i) * alpha_array(iL^S,rec_from:rec_to,i)/(alpha_sum(iL^S,rec_from:rec_to))
    end do
  
    !> left value at right interface
    wLC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)

  end subroutine WENO5limiterL

  subroutine WENO5limiterR(rec_from,rec_to,n_w,ixI^L,iL^L,idims,w,wRC,var,rec)
    use mod_global_parameters
  
    integer, intent(in)             :: rec_from, rec_to, n_w
    integer, intent(in)             :: ixI^L, iL^L, idims
    integer, intent(in)             :: var
    double precision, intent(in)    :: w(ixI^S,1:n_w)
    double precision, intent(inout) :: wRC(ixI^S,1:n_w)
    logical, intent(in)             :: rec
    !> local
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    double precision                :: f_array(ixI^S,1:n_w,3), d_array(3)
    double precision                :: beta(ixI^S,1:n_w,3), beta_coeff(2)
    double precision                :: tau(ixI^S,1:n_w), tmp(ixI^S,1:n_w)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:n_w,3), alpha_sum(ixI^S,1:n_w), flux(ixI^S,1:n_w)
    integer                         :: i,j
    double precision                :: lambda(ixI^S)

    double precision, dimension(ixI^S,1:n_w)  :: wRCtmp

    iLm^L=iL^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    lambda(iL^S)=block%mesh%dx(iL^S,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if (rec) then
      !reconstruction variation
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
      u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
      u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    else
      !interpolation variation
      d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
      u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
      u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
      u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    end if
    
    !> right side
    f_array(iL^S,rec_from:rec_to,1) = u1_coeff(1) * w(iLppp^S,rec_from:rec_to) + u1_coeff(2) * w(iLpp^S,rec_from:rec_to) + u1_coeff(3) * w(iLp^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,2) = u2_coeff(1) * w(iLpp^S,rec_from:rec_to)  + u2_coeff(2) * w(iLp^S,rec_from:rec_to)  + u2_coeff(3) * w(iL^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,3) = u3_coeff(1) * w(iLp^S,rec_from:rec_to)   + u3_coeff(2) * w(iL^S,rec_from:rec_to)   + u3_coeff(3) * w(iLm^S,rec_from:rec_to)  
  
    beta(iL^S,rec_from:rec_to,1) = beta_coeff(1) * (w(iLppp^S,rec_from:rec_to) + w(iLp^S,rec_from:rec_to) - 2.0d0*w(iLpp^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iLppp^S,rec_from:rec_to) - 4.0d0 * w(iLpp^S,rec_from:rec_to) + 3.0d0*w(iLp^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,2) = beta_coeff(1) * (w(iLpp^S,rec_from:rec_to) + w(iL^S,rec_from:rec_to) - 2.0d0 * w(iLp^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (w(iLpp^S,rec_from:rec_to) - w(iL^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,3) = beta_coeff(1) * (w(iLp^S,rec_from:rec_to) + w(iLm^S,rec_from:rec_to) - 2.0d0 * w(iL^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (3.0d0 * w(iLp^S, rec_from:rec_to) - 4.0d0 * w(iL^S,rec_from:rec_to) + w(iLm^S,rec_from:rec_to))**2
  
    alpha_sum(iL^S,rec_from:rec_to) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iL^S,rec_from:rec_to,i) = d_array(i)/(beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)**2
        alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(2) 
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3))
      do i = 1,3
        alpha_array(iL^S,rec_from:rec_to,i) = d_array(i) * (1.d0 + (tau(iL^S,rec_from:rec_to) / &
                                      (beta(iL^S,rec_from:rec_to,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(3)
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3)) * 4.d0
      do i = 1,3
        do j = rec_from,rec_to
          tmp(iL^S,j) = (tau(iL^S,j) + weno_eps_machine) / (4.d0 * beta(iL^S,j,i) + weno_eps_machine)
          alpha_array(iL^S,j,i) = d_array(i) * (1.0d0 + tmp(iL^S,j)**2 + lambda(iL^S)/tmp(iL^S,j))
          alpha_sum(iL^S,j) = alpha_sum(iL^S,j) + alpha_array(iL^S,j,i)
        end do
      end do
    end select
    flux(iL^S,rec_from:rec_to) = 0.0d0
    do i = 1,3
      flux(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to) + f_array(iL^S,rec_from:rec_to,i) * alpha_array(iL^S,rec_from:rec_to,i)/(alpha_sum(iL^S,rec_from:rec_to))
    end do
  
    !> right value at right interface
    wRC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)

  end subroutine WENO5limiterR

  subroutine WENO5NMlimiterL(rec_from,rec_to,n_w,ixI^L,iL^L,idims,w,wLC,var,rec)
    use mod_global_parameters

    integer, intent(in)             :: rec_from, rec_to, n_w
    integer, intent(in)             :: ixI^L,iL^L,idims,var
    double precision, intent(in)    :: w(ixI^S,1:n_w)
    double precision, intent(inout) :: wLC(ixI^S,1:n_w) 
    logical, intent(in)             :: rec
    !> local
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L
    integer                         :: iM^L, iMp^L
    double precision                :: f_array(ixI^S,1:n_w,3), d_array(3)
    double precision                :: beta(ixI^S,1:n_w,3), beta_coeff(2)
    double precision                :: tau(ixI^S,1:n_w), tmp(ixI^S,1:n_w)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:n_w,3), alpha_sum(ixI^S,1:n_w), flux(ixI^S,1:n_w)
    double precision                :: wc(ixI^S,1:n_w), wd(ixI^S,1:n_w)
    integer                         :: i,j
    double precision                :: lambda(ixI^S)

    double precision, dimension(ixI^S,1:n_w)  :: wLCtmp

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iMmin^D=iLmmin^D;
    iMmax^D=iLpmax^D;
    iMp^L=iM^L+kr(idims,^D);
    lambda(iL^S)=block%mesh%dx(iL^S,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if (rec) then
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ -2.d0/3.d0, -1.d0/3.d0, 2.d0 /)
      u2_coeff(1:3) = (/ -1.d0/3.d0, 2.d0/3.d0, 2.d0/3.d0 /)
      u3_coeff(1:3) = (/ 2.d0/3.d0, 2.d0/3.d0, -1.d0/3.d0 /)
    else
      error stop "not implemented"
    end if
    do i = rec_from,rec_to
      wc(iM^S,i) = (block%mesh%dx(iMp^S,idims) * w(iM^S,i) + block%mesh%dx(iM^S,idims) * w(iMp^S,i)) / &
                   (block%mesh%dx(iMp^S,idims) + block%mesh%dx(iM^S,idims))
      wd(iL^S,i) = ((2.d0 * block%mesh%dx(iLm^S,idims) + block%mesh%dx(iLmm^S,idims)) * w(iLm^S,i) - block%mesh%dx(iLm^S,idims) * w(iLmm^S,i)) / &
                   (block%mesh%dx(iLmm^S,idims) + block%mesh%dx(iLm^S,idims))
    enddo
    f_array(iL^S,rec_from:rec_to,1) = u1_coeff(1) * wd(iL^S,rec_from:rec_to)   + u1_coeff(2) * wc(iLm^S,rec_from:rec_to)+ u1_coeff(3) * w(iL^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,2) = u2_coeff(1) * wc(iLm^S,rec_from:rec_to)  + u2_coeff(2) * w(iL^S,rec_from:rec_to)  + u2_coeff(3) * wc(iL^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,3) = u3_coeff(1) * wc(iL^S,rec_from:rec_to)   + u3_coeff(2) * w(iLp^S,rec_from:rec_to) + u3_coeff(3) * wc(iLp^S,rec_from:rec_to)  
    beta(iL^S,rec_from:rec_to,1) = beta_coeff(1) * (wc(iLm^S,rec_from:rec_to) - wd(iL^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (2.d0 * w(iL^S,rec_from:rec_to) - wc(iLm^S,rec_from:rec_to) - wd(iL^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,2) = beta_coeff(1) * (wc(iLm^S,rec_from:rec_to) + wc(iL^S,rec_from:rec_to) - 2.0d0 * w(iL^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (wc(iLm^S,rec_from:rec_to) - wc(iL^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,3) = beta_coeff(1) * (wc(iL^S,rec_from:rec_to) + wc(iLp^S,rec_from:rec_to) - 2.0d0 * w(iLp^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (3.0d0 * wc(iL^S, rec_from:rec_to) - 4.0d0 * w(iLp^S,rec_from:rec_to) + wc(iLp^S,rec_from:rec_to))**2
    alpha_sum(iL^S,rec_from:rec_to) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
         alpha_array(iL^S,rec_from:rec_to,i) = d_array(i)/(4.d0 * beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)**2
         alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(2)
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3)) * 4.d0
      do i = 1,3
        alpha_array(iL^S,rec_from:rec_to,i) = d_array(i) * (1.d0 + (tau(iL^S,rec_from:rec_to) / &
                                      (4.d0 * beta(iL^S,rec_from:rec_to,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(3)
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3)) * 4.d0
      do i=1,3
        do j=rec_from,rec_to
          tmp(iL^S,j) = (tau(iL^S,j) + weno_eps_machine) / (4.d0 * beta(iL^S,j,i) + weno_eps_machine)
          alpha_array(iL^S,j,i) = d_array(i) * (1.0d0 + tmp(iL^S,j)**2 + lambda(iL^S)/tmp(iL^S,j))
          alpha_sum(iL^S,j) = alpha_sum(iL^S,j) + alpha_array(iL^S,j,i)
        end do
      end do
    end select
    flux(iL^S,rec_from:rec_to) = 0.0d0
    do i = 1,3
      flux(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to) + f_array(iL^S,rec_from:rec_to,i) * alpha_array(iL^S,rec_from:rec_to,i)/(alpha_sum(iL^S,rec_from:rec_to))
    end do
    !> left value at right interface
    wLC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)

  end subroutine WENO5NMlimiterL

  subroutine WENO5NMlimiterR(rec_from,rec_to,n_w,ixI^L,iL^L,idims,w,wRC,var,rec)
    use mod_global_parameters
  
    integer, intent(in)             :: rec_from, rec_to, n_w
    integer, intent(in)             :: ixI^L,iL^L,idims,var
    double precision, intent(in)    :: w(ixI^S,1:n_w)
    double precision, intent(inout) :: wRC(ixI^S,1:n_w)
    logical, intent(in)             :: rec
    !> local
    integer                         :: iLm^L,iLp^L,iLpp^L,iLppp^L
    integer                         :: iM^L, iMp^L
    double precision                :: f_array(ixI^S,1:n_w,3), d_array(3)
    double precision                :: beta(ixI^S,1:n_w,3), beta_coeff(2)
    double precision                :: tau(ixI^S,1:n_w), tmp(ixI^S,1:n_w)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:n_w,3), alpha_sum(ixI^S,1:n_w), flux(ixI^S,1:n_w)
    double precision                :: wc(ixI^S,1:n_w), we(ixI^S,1:n_w)
    integer                         :: i,j
    double precision                :: lambda(ixI^S)

    double precision, dimension(ixI^S,1:n_w)  :: wRCtmp

    iLm^L=iL^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    iMmin^D=iLmmin^D;
    iMmax^D=iLpmax^D;
    iMp^L=iM^L+kr(idims,^D);
    lambda(iL^S)=block%mesh%dx(iL^S,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if (rec) then
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ -2.d0/3.d0, -1.d0/3.d0, 2.d0 /)
      u2_coeff(1:3) = (/ -1.d0/3.d0, 2.d0/3.d0, 2.d0/3.d0 /)
      u3_coeff(1:3) = (/ 2.d0/3.d0, 2.d0/3.d0, -1.d0/3.d0 /)
    else
      error stop "not implemented"
    end if
    do i = rec_from,rec_to
      wc(iM^S,i) = (block%mesh%dx(iMp^S,idims) * w(iM^S,i) + block%mesh%dx(iM^S,idims) * w(iMp^S,i)) / &
                   (block%mesh%dx(iMp^S,idims) + block%mesh%dx(iM^S,idims))
      we(iL^S,i) = ((2.d0 * block%mesh%dx(iLpp^S,idims) + block%mesh%dx(iLppp^S,idims)) * w(iLpp^S,i) - block%mesh%dx(iLpp^S,idims) * w(iLppp^S,i)) / &
                   (block%mesh%dx(iLppp^S,idims) + block%mesh%dx(iLpp^S,idims))
    enddo
    !> right side
    f_array(iL^S,rec_from:rec_to,1) = u1_coeff(1) * we(iL^S,rec_from:rec_to)  + u1_coeff(2) * wc(iLp^S,rec_from:rec_to) + u1_coeff(3) * w(iLp^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,2) = u2_coeff(1) * wc(iLp^S,rec_from:rec_to) + u2_coeff(2) * w(iLp^S,rec_from:rec_to) + u2_coeff(3) * wc(iL^S,rec_from:rec_to)
    f_array(iL^S,rec_from:rec_to,3) = u3_coeff(1) * wc(iL^S,rec_from:rec_to)  + u3_coeff(2) * w(iL^S,rec_from:rec_to)  + u3_coeff(3) * wc(iLm^S,rec_from:rec_to)  
    beta(iL^S,rec_from:rec_to,1) = beta_coeff(1) * (wc(iLp^S,rec_from:rec_to) - we(iL^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (2.d0 * w(iLp^S,rec_from:rec_to) - wc(iLp^S,rec_from:rec_to) - we(iL^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,2) = beta_coeff(1) * (wc(iLp^S,rec_from:rec_to) + wc(iL^S,rec_from:rec_to) - 2.0d0 * w(iLp^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (wc(iLp^S,rec_from:rec_to) - wc(iL^S,rec_from:rec_to))**2
    beta(iL^S,rec_from:rec_to,3) = beta_coeff(1) * (wc(iL^S,rec_from:rec_to) + wc(iLm^S,rec_from:rec_to) - 2.0d0 * w(iL^S,rec_from:rec_to))**2 &
         + beta_coeff(2) * (3.0d0 * wc(iL^S, rec_from:rec_to) - 4.0d0 * w(iL^S,rec_from:rec_to) + wc(iLm^S,rec_from:rec_to))**2
    alpha_sum(iL^S,rec_from:rec_to) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iL^S,rec_from:rec_to,i) = d_array(i)/(4.d0 * beta(iL^S,rec_from:rec_to,i) + weno_eps_machine)**2
        alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(2) 
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3)) * 4.d0
      do i = 1,3
        alpha_array(iL^S,rec_from:rec_to,i) = d_array(i) * (1.d0 + (tau(iL^S,rec_from:rec_to) / &
                                      (4.d0 * beta(iL^S,rec_from:rec_to,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,rec_from:rec_to) = alpha_sum(iL^S,rec_from:rec_to) + alpha_array(iL^S,rec_from:rec_to,i)
      end do
    case(3)
      tau(iL^S,rec_from:rec_to) = abs(beta(iL^S,rec_from:rec_to,1) - beta(iL^S,rec_from:rec_to,3)) * 4.d0
      do i = 1,3
        do j = rec_from,rec_to
          tmp(iL^S,j) = (tau(iL^S,j) + weno_eps_machine) / (4.d0 * beta(iL^S,j,i) + weno_eps_machine)
          alpha_array(iL^S,j,i) = d_array(i) * (1.0d0 + tmp(iL^S,j)**2 + lambda(iL^S)/tmp(iL^S,j))
          alpha_sum(iL^S,j) = alpha_sum(iL^S,j) + alpha_array(iL^S,j,i)
        end do
      end do
    end select
    flux(iL^S,rec_from:rec_to) = 0.0d0
    do i = 1,3
      flux(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to) + f_array(iL^S,rec_from:rec_to,i) * alpha_array(iL^S,rec_from:rec_to,i)/(alpha_sum(iL^S,rec_from:rec_to))
    end do
    !> right value at right interface
    wRC(iL^S,rec_from:rec_to) = flux(iL^S,rec_from:rec_to)

  end subroutine WENO5NMlimiterR

  subroutine WENO5CU6limiter(ixI^L,iL^L,idims,w,wLC,wRC,rec)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S)
    double precision, intent(inout) :: wRC(ixI^S),wLC(ixI^S) 
    logical, intent(in)             :: rec
    !> local
    integer :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    double precision :: f_array(ixI^S,3), d_array(3)
    double precision :: beta(ixI^S,3), beta_coeff(2)
    double precision :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision :: alpha_array(ixI^S,3), alpha_sum(ixI^S)
    double precision :: theta2(ixI^S)
    integer :: i
    double precision, parameter :: theta_limit=0.7d0

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);

    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if (rec) then
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
      u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
      u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    else
      error stop "not implemented"
    end if
    
    !> left side
    beta(iL^S,1)=beta_coeff(1)*(w(iLmm^S)+w(iL^S)-2.0d0*w(iLm^S))**2&
        +beta_coeff(2)*(w(iLmm^S)-4.0d0*w(iLm^S)+3.0d0*w(iL^S))**2
    beta(iL^S,2)=beta_coeff(1)*(w(iLm^S)+w(iLp^S)-2.0d0*w(iL^S))**2&
        +beta_coeff(2)*(w(iLm^S)-w(iLp^S))**2
    beta(iL^S,3)=beta_coeff(1)*(w(iL^S)+w(iLpp^S)-2.0d0*w(iLp^S))**2&
        +beta_coeff(2)*(3.0d0*w(iL^S)-4.0d0*w(iLp^S)+w(iLpp^S))**2
    alpha_sum(iL^S)=zero
    do i=1,3
      alpha_array(iL^S,i)=d_array(i)/(beta(iL^S,i)+weno_eps_machine)**2
      alpha_sum(iL^S)=alpha_sum(iL^S)+alpha_array(iL^S,i)
    end do
    do i=1,3
      alpha_array(iL^S,i)=alpha_array(iL^S,i)/alpha_sum(iL^S)
    end do
    theta2(iL^S)=((alpha_array(iL^S,1)/d_array(1)-1.d0)**2&
                          +(alpha_array(iL^S,2)/d_array(2)-1.d0)**2&
                          +(alpha_array(iL^S,3)/d_array(3)-1.d0)**2)/83.d0
    where(theta2(iL^S) .gt. theta_limit)
      f_array(iL^S,1)=u1_coeff(1)*w(iLmm^S)+u1_coeff(2)*w(iLm^S)+u1_coeff(3)*w(iL^S)
      f_array(iL^S,2)=u2_coeff(1)*w(iLm^S)+u2_coeff(2)*w(iL^S)+u2_coeff(3)*w(iLp^S)
      f_array(iL^S,3)=u3_coeff(1)*w(iL^S)+u3_coeff(2)*w(iLp^S)+u3_coeff(3)*w(iLpp^S)  
      wLC(iL^S)=f_array(iL^S,1)*alpha_array(iL^S,1)&
                        +f_array(iL^S,2)*alpha_array(iL^S,2)&
                        +f_array(iL^S,3)*alpha_array(iL^S,3)
    else where
      wLC(iL^S)=1.d0/60.d0*(w(iLmm^S)-8.d0*w(iLm^S)+37.d0*w(iL^S)&
                         +37.d0*w(iLp^S)-8.d0*w(iLpp^S)+w(iLppp^S))
    end where
  
    !> right side
    beta(iL^S,1)=beta_coeff(1)*(w(iLppp^S)+w(iLp^S)-2.0d0*w(iLpp^S))**2&
         +beta_coeff(2)*(w(iLppp^S)-4.0d0*w(iLpp^S)+3.0d0*w(iLp^S))**2
    beta(iL^S,2)=beta_coeff(1)*(w(iLpp^S)+w(iL^S)-2.0d0*w(iLp^S))**2&
         +beta_coeff(2)*(w(iLpp^S)-w(iL^S))**2
    beta(iL^S,3)=beta_coeff(1)*(w(iLp^S)+w(iLm^S)-2.0d0*w(iL^S))**2&
         +beta_coeff(2)*(3.0d0*w(iLp^S)-4.0d0*w(iL^S)+w(iLm^S))**2
    alpha_sum(iL^S)=zero
    do i=1,3
      alpha_array(iL^S,i)=d_array(i)/(beta(iL^S,i)+weno_eps_machine)**2
      alpha_sum(iL^S)=alpha_sum(iL^S)+alpha_array(iL^S,i)
    end do
    do i=1,3
      alpha_array(iL^S,i)=alpha_array(iL^S,i)/alpha_sum(iL^S)
    end do
    theta2(iL^S)=((alpha_array(iL^S,1)/d_array(1)-1.d0)**2&
                          +(alpha_array(iL^S,2)/d_array(2)-1.d0)**2&
                          +(alpha_array(iL^S,3)/d_array(3)-1.d0)**2)/83.d0
    where(theta2(iL^S) .gt. theta_limit)
      f_array(iL^S,1)=u1_coeff(1)*w(iLppp^S)+u1_coeff(2)*w(iLpp^S)+u1_coeff(3)*w(iLp^S)
      f_array(iL^S,2)=u2_coeff(1)*w(iLpp^S)+u2_coeff(2)*w(iLp^S)+u2_coeff(3)*w(iL^S)
      f_array(iL^S,3)=u3_coeff(1)*w(iLp^S)+u3_coeff(2)*w(iL^S)+u3_coeff(3)*w(iLm^S)  
      wRC(iL^S)=f_array(iL^S,1)*alpha_array(iL^S,1)&
                        +f_array(iL^S,2)*alpha_array(iL^S,2)&
                        +f_array(iL^S,3)*alpha_array(iL^S,3)
    else where
      wRC(iL^S)=1.d0/60.d0*(w(iLppp^S)-8.d0*w(iLpp^S)+37.d0*w(iLp^S)&
                        +37.d0*w(iL^S)-8.d0*w(iLm^S)+w(iLmm^S))
    end where
  end subroutine WENO5CU6limiter

  subroutine WENO7limiter(ixI^L,iL^L,idims,w,wLC,wRC,var,rec)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims, var
    double precision, intent(in)    :: w(ixI^S)
    double precision, intent(inout) :: wRC(ixI^S),wLC(ixI^S) 
    logical, intent(in)             :: rec
    !> local
    integer                         :: iLm^L, iLmm^L, iLmmm^L
    integer                         :: iLp^L, iLpp^L, iLppp^L, iLpppp^L
    integer                         :: id^L, idp^L, idpp^L, idm^L, ie^L, iem^L, iep^L, iepp^L

    double precision, dimension(4)       :: d_array, u1_coeff, u2_coeff, u3_coeff, u4_coeff
    double precision, dimension(ixI^S,4) :: f_array, beta, alpha_array
    double precision, dimension(ixI^S)   :: a, b, c, tmp, tmp2, tmp3
    double precision, dimension(ixI^S)   :: alpha_sum, d, dm4
    double precision, dimension(ixI^S)   :: flux, flux_min, flux_max, flux_ul, flux_md, flux_lc
    integer                         :: i
    double precision, parameter     :: mpalpha = 2.d0, mpbeta = 4.d0

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLmmm^L=iLmm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    iLpppp^L=iLppp^L+kr(idims,^D);

    if (rec) then
      d_array(1:4) = (/ 1.d0/35.d0, 12.d0/35.d0, 18.d0/35.d0, 4.d0/35.d0 /)
      u1_coeff(1:4) = (/ -1.d0/4.d0, 13.d0/12.d0, -23.d0/12.d0, 25.d0/12.d0 /)
      u2_coeff(1:4) = (/ 1.d0/12.d0, -5.d0/12.d0, 13.d0/12.d0, 1.d0/4.d0 /)
      u3_coeff(1:4) = (/ -1.d0/12.d0, 7.d0/12.d0, 7.d0/12.d0, -1.d0/12.d0 /)
      u4_coeff(1:4) = (/ 1.d0/4.d0, 13.d0/12.d0, -5.d0/12.d0, 1.d0/12.d0 /)
    else
      error stop "not implemented"
    end if
    
    !> left side
    f_array(iL^S,1) = u1_coeff(1) * w(iLmmm^S) + u1_coeff(2) * w(iLmm^S) + u1_coeff(3) * w(iLm^S) &
                             + u1_coeff(4) * w(iL^S)
    f_array(iL^S,2) = u2_coeff(1) * w(iLmm^S)  + u2_coeff(2) * w(iLm^S)  + u2_coeff(3) * w(iL^S)  &
                             + u2_coeff(4) * w(iLp^S)
    f_array(iL^S,3) = u3_coeff(1) * w(iLm^S)   + u3_coeff(2) * w(iL^S)   + u3_coeff(3) * w(iLp^S)   &
                             + u3_coeff(4) * w(iLpp^S)
    f_array(iL^S,4) = u4_coeff(1) * w(iL^S)    + u4_coeff(2) * w(iLp^S)    + u4_coeff(3) * w(iLpp^S)  &
                             + u4_coeff(4) * w(iLppp^S)
  
    beta(iL^S,1) = w(iLmmm^S) * (547.d0 * w(iLmmm^S) - 3882.d0 * w(iLmm^S) + 4642.d0 * w(iLm^S) &
                          - 1854.d0 * w(iL^S)) &
                          + w(iLmm^S) * (7043.d0 * w(iLmm^S) - 17246.d0 * w(iLm^S) + 7042.d0 * w(iL^S)) &
                          + w(iLm^S) * (11003.d0 * w(iLm^S) - 9402.d0 * w(iL^S)) + 2107.d0 * w(iL^S)**2
    beta(iL^S,2) = w(iLmm^S) * (267.d0 * w(iLmm^S) - 1642.d0 * w(iLm^S) + 1602.d0 * w(iL^S) &
                          - 494.d0 * w(iLp^S))  &
                          + w(iLm^S) * (2843.d0 * w(iLm^S) - 5966.d0 * w(iL^S) + 1922.d0 * w(iLp^S)) &
                          + w(iL^S) * (3443.d0 * w(iL^S) - 2522.d0 * w(iLp^S)) + 547.d0 * w(iLp^S) ** 2
    beta(iL^S,3) = w(iLm^S) * (547.d0 * w(iLm^S) - 2522.d0 * w(iL^S) + 1922.d0 * w(iLp^S) &
                          - 494.d0 * w(iLpp^S))  &
                          + w(iL^S) * (3443.d0 * w(iL^S) - 5966.d0 * w(iLp^S) + 1602.d0 * w(iLpp^S)) &
                          + w(iLp^S) * (2843.d0 * w(iLp^S) - 1642.d0 * w(iLpp^S)) + 267.d0 * w(iLpp^S) ** 2
    beta(iL^S,4) = w(iL^S) * (2107.d0 * w(iL^S) - 9402.d0 * w(iLp^S) + 7042.d0 * w(iLpp^S) &
                          - 1854.d0 * w(iLppp^S))  &
                          + w(iLp^S) * (11003.d0 * w(iLp^S) - 17246.d0 * w(iLpp^S) + 4642.d0 * w(iLppp^S)) &
                          + w(iLpp^S) * (7043.d0 * w(iLpp^S) - 3882.d0 * w(iLppp^S)) &
                          + 547.d0 * w(iLppp^S) ** 2

    alpha_sum(iL^S) = 0.0d0 
    do i = 1,4
       alpha_array(iL^S,i) = d_array(i)/(beta(iL^S,i) + weno_eps_machine)**2
       alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
    end do
    flux(iL^S) = 0.0d0
    do i = 1,4
       flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
    end do
    
    select case(var)
    ! case1 for wenojs, case2 for mpweno
    case(1) 
      wLC(iL^S) = flux(iL^S)
    case(2)
      idmax^D=iLmax^D; idmin^D=iLmin^D-kr(idims,^D);
      idm^L=id^L-kr(idims,^D);
      idp^L=id^L+kr(idims,^D);
  
      iemax^D=idmax^D+kr(idims,^D); iemin^D=idmin^D;
      iem^L=ie^L-kr(idims,^D);
      iep^L=ie^L+kr(idims,^D);
  
      d(ie^S) = w(iep^S)-2.0d0*w(ie^S)+w(iem^S)
  
      a(id^S) = 4.0d0*d(id^S)-d(idp^S)
      b(id^S) = 4.0d0*d(idp^S)-d(id^S)
      call minmod(ixI^L,id^L,a,b,tmp)
      a(id^S) = d(id^S)
      b(id^S) = d(idp^S)
      call minmod(ixI^L,id^L,a,b,tmp2)
      call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
      dm4(id^S) = tmp3(id^S)

      flux_ul(iL^S) = w(iL^S) + mpalpha * (w(iL^S) - w(iLm^S))
      flux_md(iL^S) = half * (w(iL^S) + w(iLp^S) - dm4(iL^S))
      flux_lc(iL^S) = half * (3.d0 * w(iL^S) - w(iLm^S)) + mpbeta / 3.d0 * dm4(iLm^S)
    
      flux_min(iL^S) = max(min(w(iL^S), w(iLp^S), flux_md(iL^S)), &
                                min(w(iL^S), flux_ul(iL^S),flux_lc(iL^S)))
  
      flux_max(iL^S) = min(max(w(iL^S), w(iLp^S), flux_md(iL^S)), &
                                max(w(iL^S), flux_ul(iL^S),flux_lc(iL^S)))
      a(iL^S) = flux(iL^S)
      b(iL^S) = flux_min(iL^S)
      c(iL^S) = flux_max(iL^S)
      call median(ixI^L, iL^L, a, b, c, tmp) 
      wLC(iL^S) = tmp(iL^S)
    end select

    !> right side
    !>> mmm -> pppp
    !>> mm  -> ppp
    !>> m   -> pp
    !>> 0   -> p
    !>> p   -> 0
    !>> pp  -> m
    !>> ppp -> mm
    f_array(iL^S,1) = u1_coeff(1) * w(iLpppp^S) + u1_coeff(2) * w(iLppp^S) + u1_coeff(3) * w(iLpp^S) &
                             + u1_coeff(4) * w(iLp^S)
    f_array(iL^S,2) = u2_coeff(1) * w(iLppp^S)  + u2_coeff(2) * w(iLpp^S)  + u2_coeff(3) * w(iLp^S)  &
                             + u2_coeff(4) * w(iL^S)
    f_array(iL^S,3) = u3_coeff(1) * w(iLpp^S)   + u3_coeff(2) * w(iLp^S)   + u3_coeff(3) * w(iL^S)   &
                             + u3_coeff(4) * w(iLm^S)
    f_array(iL^S,4) = u4_coeff(1) * w(iLp^S)    + u4_coeff(2) * w(iL^S)    + u4_coeff(3) * w(iLm^S)  &
                             + u4_coeff(4) * w(iLmm^S)

    beta(iL^S,1) = w(iLpppp^S) * (547.d0 * w(iLpppp^S) - 3882.d0 * w(iLppp^S) + 4642.d0 * w(iLpp^S) &
                          - 1854.d0 * w(iLp^S)) &
                          + w(iLppp^S) * (7043.d0 * w(iLppp^S) - 17246.d0 * w(iLpp^S) + 7042.d0 * w(iLp^S)) &
                          + w(iLpp^S) * (11003.d0 * w(iLpp^S) - 9402.d0 * w(iLp^S)) + 2107.d0 * w(iLp^S)**2
    beta(iL^S,2) = w(iLppp^S) * (267.d0 * w(iLppp^S) - 1642.d0 * w(iLpp^S) + 1602.d0 * w(iLp^S) &
                          - 494.d0 * w(iL^S))  &
                          + w(iLpp^S) * (2843.d0 * w(iLpp^S) - 5966.d0 * w(iLp^S) + 1922.d0 * w(iL^S)) &
                          + w(iLp^S) * (3443.d0 * w(iLp^S) - 2522.d0 * w(iL^S)) + 547.d0 * w(iL^S) ** 2
    beta(iL^S,3) = w(iLpp^S) * (547.d0 * w(iLpp^S) - 2522.d0 * w(iLp^S) + 1922.d0 * w(iL^S) &
                          - 494.d0 * w(iLm^S))  &
                          + w(iLp^S) * (3443.d0 * w(iLp^S) - 5966.d0 * w(iL^S) + 1602.d0 * w(iLm^S)) &
                          + w(iL^S) * (2843.d0 * w(iL^S) - 1642.d0 * w(iLm^S)) + 267.d0 * w(iLm^S) ** 2
    beta(iL^S,4) = w(iLp^S) * (2107.d0 * w(iLp^S) - 9402.d0 * w(iL^S) + 7042.d0 * w(iLm^S) &
                          - 1854.d0 * w(iLmm^S))  &
                          + w(iL^S) * (11003.d0 * w(iL^S) - 17246.d0 * w(iLm^S) + 4642.d0 * w(iLmm^S)) &
                          + w(iLm^S) * (7043.d0 * w(iLm^S) - 3882.d0 * w(iLmm^S)) + 547.d0 * w(iLmm^S) ** 2

    alpha_sum(iL^S) = 0.0d0 
    do i = 1,4
       alpha_array(iL^S,i) = d_array(i)/(beta(iL^S,i) + weno_eps_machine)**2
       alpha_sum(iL^S) = alpha_sum(iL^S) + alpha_array(iL^S,i)
    end do
    flux(iL^S) = 0.0d0
    do i = 1,4
       flux(iL^S) = flux(iL^S) + f_array(iL^S,i) * alpha_array(iL^S,i)/(alpha_sum(iL^S))
    end do

    select case(var)
    case(1)
      wRC(iL^S) = flux(iL^S)
    case(2)
      idmax^D=iLmax^D+kr(idims,^D); idmin^D=iLmin^D;
      idm^L=id^L-kr(idims,^D);
      idp^L=id^L+kr(idims,^D);
  
      iemax^D=idmax^D; iemin^D=idmin^D-kr(idims,^D);
      iem^L=ie^L-kr(idims,^D);
      iep^L=ie^L+kr(idims,^D);
      iepp^L=iep^L+kr(idims,^D);
  
      d(ie^S) = w(ie^S)-2.0d0*w(iep^S)+w(iepp^S)
  
      a(id^S) = 4.0d0*d(id^S)-d(idm^S)
      b(id^S) = 4.0d0*d(idm^S)-d(id^S)
      call minmod(ixI^L,id^L,a,b,tmp)
      a(id^S) = d(id^S)
      b(id^S) = d(idm^S)
      call minmod(ixI^L,id^L,a,b,tmp2)
      call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
      dm4(id^S) = tmp3(id^S)
   
      flux_ul(iL^S) = w(iLp^S) + mpalpha * (w(iLp^S) - w(iLpp^S))
      flux_md(iL^S) = half * (w(iL^S) + w(iLp^S) - dm4(iL^S))
      flux_lc(iL^S) = half * (3.d0 * w(iLp^S) - w(iLpp^S)) + mpbeta / 3.d0 * dm4(iLp^S)
    
      flux_min(iL^S) = max(min(w(iLp^S), w(iL^S), flux_md(iL^S)), &
                                min(w(iLp^S), flux_ul(iL^S),flux_lc(iL^S)))
  
      flux_max(iL^S) = min(max(w(iLp^S), w(iL^S), flux_md(iL^S)), &
                                max(w(iLp^S), flux_ul(iL^S),flux_lc(iL^S)))
      a(iL^S) = flux(iL^S)
      b(iL^S) = flux_min(iL^S)
      c(iL^S) = flux_max(iL^S)
      call median(ixI^L, iL^L, a, b, c, tmp) 
      wRC(iL^S) = tmp(iL^S)
    end select
  end subroutine WENO7limiter

  subroutine exENO7limiter(ixI^L,iL^L,idims,w,wLC,wRC,rec)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S)
    double precision, intent(inout) :: wRC(ixI^S),wLC(ixI^S) 
    logical, intent(in)             :: rec
    !> local
    integer                         :: i
    integer                         :: iLm^L, iLmm^L, iLmmm^L
    integer                         :: iLp^L, iLpp^L, iLppp^L, iLpppp^L
    integer                         :: iM^L, iMm^L, iMmm^L
    integer                         :: iMp^L, iMpp^L, iMppp^L
    integer                         :: id^L, idp^L, idpp^L, idm^L, ie^L, iem^L, iep^L, iepp^L
    integer, dimension(ixI^S)            :: delta_sum
    integer, dimension(ixI^S,3)          :: delta
    double precision, dimension(2)       :: beta_coeff
    double precision, dimension(3)       :: d_array
    double precision, dimension(ixI^S)   :: gamma_sum, flux
    double precision, dimension(ixI^S,3) :: beta, gamma_array, kai_array
    double precision, parameter          :: exeno_ct = 1.d-1

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLmmm^L=iLmm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    iLpppp^L=iLppp^L+kr(idims,^D);

    iMmin^D=iLmin^D-kr(idims,^D);
    iMmax^D=iLmax^D+kr(idims,^D);
    iMm^L=iM^L-kr(idims,^D);
    iMmm^L=iMm^L-kr(idims,^D);
    iMp^L=iM^L+kr(idims,^D);
    iMpp^L=iMp^L+kr(idims,^D);
    iMppp^L=iMpp^L+kr(idims,^D);

    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if (rec) then
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    else
      error stop "not implemented"
    end if

    !>> left side
    beta(iM^S,1) = beta_coeff(1) * (w(iMmm^S) + w(iM^S) - 2.0d0 * w(iMm^S))**2 &
         + beta_coeff(2) * (w(iMmm^S) - 4.0d0 * w(iMm^S) + 3.0d0 * w(iM^S))**2
    beta(iM^S,2) = beta_coeff(1) * (w(iMm^S) + w(iMp^S) - 2.0d0 * w(iM^S))**2 &
         + beta_coeff(2) * (w(iMm^S) - w(iMp^S))**2
    beta(iM^S,3) = beta_coeff(1) * (w(iM^S) + w(iMpp^S) - 2.0d0 * w(iMp^S))**2 &
         + beta_coeff(2) * (3.0d0 * w(iM^S) - 4.0d0 * w(iMp^S) + w(iMpp^S))**2

    gamma_sum(iM^S) = 0.0d0 
    do i = 1,3
      gamma_array(iM^S,i) = d_array(i) / (beta(iM^S,i) + weno_eps_machine)**2
      gamma_sum(iM^S) = gamma_sum(iM^S) + gamma_array(iM^S,i)
    end do
    do i = 1,3
      kai_array(iM^S,i) = gamma_array(iM^S,i) / gamma_sum(iM^S)
      where(kai_array(iM^S,i) .lt. exeno_ct) 
        delta(iM^S,i) = 0
      elsewhere
        delta(iM^S,i) = 1
      endwhere
    end do

    delta_sum(iL^S) = delta(iLm^S,1) * 1 + delta(iLp^S,3) * 2 + delta(iL^S,1) * 4 &
                             + delta(iL^S,2)  * 8 + delta(iL^S,3) * 16

    !> f3
    where(delta_sum(iL^S) .eq. 31)
      flux(iL^S) = (- 3.d0 * w(iLmmm^S) + 25.d0 * w(iLmm^S) - 101.d0 * w(iLm^S) + 319.d0 * w(iL^S) &
                             + 214.d0 * w(iLp^S) - 38.d0 * w(iLpp^S) + 4.d0 * w(iLppp^S)) / 420.d0
    !> f4
    elsewhere(delta_sum(iL^S) .eq. 30)
      flux(iL^S) = (w(iLmm^S) - 8.d0 * w(iLm^S) + 37.d0 * w(iL^S) + 37.d0 * w(iLp^S) &
                             - 8.d0 * w(iLpp^S) + w(iLppp^S)) / 60.d0
    !> f5
    elsewhere(delta_sum(iL^S) .eq. 29)
      flux(iL^S) = (- w(iLmmm^S) + 7.d0 * w(iLmm^S) - 23.d0 * w(iLm^S) + 57.d0 * w(iL^S) &
                             + 22.d0 * w(iLp^S) - 2.d0 * w(iLpp^S)) / 60.d0
    !> f6
    elsewhere(delta_sum(iL^S) .eq. 28)
      flux(iL^S) = (2.d0 * w(iLmm^S) - 13.d0 * w(iLm^S) + 47.d0 * w(iL^S) + 27.d0 * w(iLp^S) &
                             - 3.d0 * w(iLpp^S)) / 60.d0
    !> f7
    elsewhere(delta_sum(iL^S) .ge. 24)
      flux(iL^S) = (- w(iLm^S) + 7.d0 * w(iL^S) + 7.d0 * w(iLp^S) - w(iLpp^S)) / 12.d0
    !> f9
    elsewhere(delta_sum(iL^S) .ge. 16)
      flux(iL^S) = (2.d0 * w(iL^S) + 5.d0 * w(iLp^S) - w(iLpp^S)) / 6.d0
    !> f8
    elsewhere(delta_sum(iL^S) .ge. 12)
      flux(iL^S) = (w(iLmm^S) - 5.d0 * w(iLm^S) + 13.d0 * w(iL^S) + 3.d0 * w(iLp^S)) / 12.d0
    !> f10
    elsewhere(delta_sum(iL^S) .ge. 8)
      flux(iL^S) = (- w(iLm^S) + 5.d0 * w(iL^S) + 2.d0 * w(iLp^S)) / 6.d0
    !> f11
    elsewhere
     flux(iL^S) = (2.d0 * w(iLmm^S) - 7.d0 * w(iLm^S) + 11.d0 * w(iL^S)) / 6.d0
    endwhere

    wLC(iL^S) = flux(iL^S)

    !> right side
    beta(iM^S,1) = beta_coeff(1) * (w(iMppp^S) + w(iMp^S) - 2.0d0 * w(iMpp^S))**2 &
         + beta_coeff(2) * (w(iMppp^S) - 4.0d0 * w(iMpp^S) + 3.0d0 * w(iMp^S))**2
    beta(iM^S,2) = beta_coeff(1) * (w(iMpp^S) + w(iM^S) - 2.0d0 * w(iMp^S))**2 &
         + beta_coeff(2) * (w(iMpp^S) - w(iM^S))**2
    beta(iM^S,3) = beta_coeff(1) * (w(iMp^S) + w(iMm^S) - 2.0d0 * w(iM^S))**2 &
         + beta_coeff(2) * (3.0d0 * w(iMp^S) - 4.0d0 * w(iM^S) + w(iMm^S))**2

    gamma_sum(iM^S) = 0.0d0 
    do i = 1,3
      gamma_array(iM^S,i) = d_array(i) / (beta(iM^S,i) + weno_eps_machine)**2
      gamma_sum(iM^S) = gamma_sum(iM^S) + gamma_array(iM^S,i)
    end do
    do i = 1,3
      kai_array(iM^S,i) = gamma_array(iM^S,i) / gamma_sum(iM^S)
      where(kai_array(iM^S,i) .lt. exeno_ct) 
        delta(iM^S,i) = 0
      elsewhere
        delta(iM^S,i) = 1
      endwhere
    end do
 
    delta_sum(iL^S) = delta(iLp^S,1) * 1 + delta(iLm^S,3) * 2 + delta(iL^S,1) * 4 &
                             + delta(iL^S,2)  * 8 + delta(iL^S,3) * 16

    where(delta_sum(iL^S) .eq. 31)
      flux(iL^S) = (- 3.d0 * w(iLpppp^S) + 25.d0 * w(iLppp^S) - 101.d0 * w(iLpp^S) &
                             + 319.d0 * w(iLp^S) + 214.d0 * w(iL^S) - 38.d0 * w(iLm^S) &
                             + 4.d0 * w(iLmm^S)) / 420.d0
    elsewhere(delta_sum(iL^S) .eq. 30)
      flux(iL^S) = (w(iLppp^S) - 8.d0 * w(iLpp^S) + 37.d0 * w(iLp^S) + 37.d0 * w(iL^S) &
                             - 8.d0 * w(iLm^S) + w(iLmm^S)) / 60.d0
    elsewhere(delta_sum(iL^S) .eq. 29)
      flux(iL^S) = (- w(iLpppp^S) + 7.d0 * w(iLppp^S) - 23.d0 * w(iLpp^S) + 57.d0 * w(iLp^S) &
                             + 22.d0 * w(iL^S) - 2.d0 * w(iLm^S)) / 60.d0
    elsewhere(delta_sum(iL^S) .eq. 28)
      flux(iL^S) = (2.d0 * w(iLppp^S) - 13.d0 * w(iLpp^S) + 47.d0 * w(iLp^S) + 27.d0 * w(iL^S) &
                             - 3.d0 * w(iLm^S)) / 60.d0
    elsewhere(delta_sum(iL^S) .ge. 24)
      flux(iL^S) = (- w(iLpp^S) + 7.d0 * w(iLp^S) + 7.d0 * w(iL^S) - w(iLm^S)) / 12.d0
    elsewhere(delta_sum(iL^S) .ge. 16)
      flux(iL^S) = (2.d0 * w(iLp^S) + 5.d0 * w(iL^S) - w(iLm^S)) / 6.d0
    elsewhere(delta_sum(iL^S) .ge. 12)
      flux(iL^S) = (w(iLppp^S) - 5.d0 * w(iLpp^S) + 13.d0 * w(iLp^S) + 3.d0 * w(iL^S)) / 12.d0
    elsewhere(delta_sum(iL^S) .ge. 8)
      flux(iL^S) = (- w(iLpp^S) + 5.d0 * w(iLp^S) + 2.d0 * w(iL^S)) / 6.d0
    elsewhere
     flux(iL^S) = (2.d0 * w(iLppp^S) - 7.d0 * w(iLpp^S) + 11.d0 * w(iLp^S)) / 6.d0
    endwhere

    wRC(iL^S) = flux(iL^S)

  end subroutine exENO7limiter

  subroutine minmod(ixI^L,ixO^L,a,b,minm)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: a(ixI^S), b(ixI^S)
    double precision, intent(out):: minm(ixI^S)
    minm(ixO^S) = (sign(one,a(ixO^S))+sign(one,b(ixO^S)))/2.0d0 * &
         min(abs(a(ixO^S)),abs(b(ixO^S)))
  end subroutine minmod

  subroutine median(ixI^L,ixO^L,a,b,c,med)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: a(ixI^S), b(ixI^S), c(ixI^S)
    double precision, intent(out):: med(ixI^S)
    double precision             :: tmp1(ixI^S),tmp2(ixI^S)
    tmp1(ixO^S) = b(ixO^S) - a(ixO^S); tmp2(ixO^S) = c(ixO^S) - a(ixO^S)
    med(ixO^S) = a(ixO^S) + (sign(one,tmp1(ixO^S))+sign(one,tmp2(ixO^S)))/2.0d0 * &
         min(abs(tmp1(ixO^S)),abs(tmp2(ixO^S)))
  end subroutine median
end module mod_weno
