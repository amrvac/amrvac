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
  ! 2021.6.12 generic treatment for fixing unphysical reconstructions
  ! 2022.10.21 remove exENO7
   
  implicit none
  private

  double precision, parameter     :: weno_eps_machine = 1.0d-18
  double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0

  public :: WENO3limiter
  public :: WENO5limiter
  public :: WENO5NMlimiter
  public :: WENO5limiterL
  public :: WENO5NMlimiterL
  public :: WENO5limiterR
  public :: WENO5NMlimiterR
  public :: TENO5ADlimiter
  public :: WENO5CU6limiter
  public :: WENO7limiter
  public :: fix_limiter
  public :: fix_limiter1
  public :: fix_onelimiter
  public :: fix_onelimiter1

contains

  subroutine fix_onelimiter(ixI^L,iL^L,wCin,wCout)
    use mod_global_parameters
    use mod_physics, only: phys_check_w

    integer, intent(in)             :: ixI^L, iL^L
    double precision, intent(in)    :: wCin(ixI^S,1:nw)
    double precision, intent(inout) :: wCout(ixI^S,1:nw)

    integer :: iw
    logical :: flagC(ixI^S,1:nw), flag(ixI^S)

    ! When limiter not TVD, negative pressures or densities could result.
    ! Fall back to flat interpolation 
    ! flagC(*,iw) indicates failed state (T when failed)
    ! assumes wCin contains primitive variables
    call phys_check_w(.true.,ixI^L,iL^L,wCin,flagC)

    ! collect all failures
    flag(iL^S)=.false.
    do iw=1,nw_recon
       flag(iL^S)=flag(iL^S).or.(flagC(iL^S,iw))
    end do
    ! only use WENO reconstructions when no field failed
    ! in other places: do not modify the initial state in wCout
    do iw=1,nw_recon
       where (flag(iL^S) .eqv. .false.)
          wCout(iL^S,iw)=wCin(iL^S,iw)
       end where
    enddo
  
  end subroutine fix_onelimiter

  subroutine fix_onelimiter1(ixI^L,iL^L,wCin,wCout)
    use mod_global_parameters
    use mod_physics, only: phys_check_w

    integer, intent(in)             :: ixI^L, iL^L
    double precision, intent(in)    :: wCin(ixI^S,1:nw)
    double precision, intent(inout) :: wCout(ixI^S,1:nw)

    integer :: iw
    logical :: flagC(ixI^S,1:nw)

    ! When limiter not TVD, negative pressures or densities could result.
    ! Fall back to flat interpolation 
    ! flagC(*,iw) indicates failed state (T when failed)
    ! assumes wCin contains primitive variables
    call phys_check_w(.true.,ixI^L,iL^L,wCin,flagC)

    ! only use WENO reconstructions when no field failed
    ! in other places: do not modify the initial state in wCout
    do iw=1,nw_recon
       where (flagC(iL^S,iw) .eqv. .false.)
          wCout(iL^S,iw)=wCin(iL^S,iw)
       end where
    enddo
  
  end subroutine fix_onelimiter1

  subroutine fix_limiter(ixI^L,iL^L,wLCin,wRCin,wLCout,wRCout)
    use mod_global_parameters
    use mod_physics, only: phys_check_w

    integer, intent(in)             :: ixI^L, iL^L
    double precision, intent(in)    :: wRCin(ixI^S,1:nw),wLCin(ixI^S,1:nw) 
    double precision, intent(inout) :: wRCout(ixI^S,1:nw),wLCout(ixI^S,1:nw) 

    integer :: iw
    logical :: flagL(ixI^S,1:nw), flagR(ixI^S,1:nw), flag(ixI^S)

    ! When limiter not TVD, negative pressures or densities could result.
    ! Fall back to flat interpolation 
    ! flagL(*,iw) indicates failed L state (T when failed)
    ! flagR(*,iw) indicates failed R state (T when failed)
    ! assumes wLCin and wRCin contain primitive variables
    call phys_check_w(.true.,ixI^L,iL^L,wLCin,flagL)
    call phys_check_w(.true.,ixI^L,iL^L,wRCin,flagR)

    ! collect all failures
    flag(iL^S)=.false.
    do iw=1,nw_recon
       flag(iL^S)=flag(iL^S).or.(flagL(iL^S,iw).or.flagR(iL^S,iw))
    end do
    ! only use WENO reconstructions L and R when no neighbour field failed
    ! in other places: do not modify the initial state in wLCout/wRCout
    do iw=1,nw_recon
       where (flag(iL^S) .eqv. .false.)
          wLCout(iL^S,iw)=wLCin(iL^S,iw)
          wRCout(iL^S,iw)=wRCin(iL^S,iw)
       end where
    enddo
  
  end subroutine fix_limiter

  subroutine fix_limiter1(ixI^L,iL^L,wLCin,wRCin,wLCout,wRCout)
    use mod_global_parameters
    use mod_physics, only: phys_check_w

    integer, intent(in)             :: ixI^L, iL^L
    double precision, intent(in)    :: wRCin(ixI^S,1:nw),wLCin(ixI^S,1:nw) 
    double precision, intent(inout) :: wRCout(ixI^S,1:nw),wLCout(ixI^S,1:nw) 

    integer :: iw
    logical :: flagL(ixI^S,1:nw), flagR(ixI^S,1:nw)

    ! When limiter not TVD, negative pressures or densities could result.
    ! Fall back to flat interpolation 
    ! flagL(*,iw) indicates failed L state (T when failed)
    ! flagR(*,iw) indicates failed R state (T when failed)
    ! assumes wLCin and wRCin contain primitive variables
    call phys_check_w(.true.,ixI^L,iL^L,wLCin,flagL)
    call phys_check_w(.true.,ixI^L,iL^L,wRCin,flagR)

    ! only use WENO reconstructions L and R when no neighbour field failed
    ! in other places: do not modify the initial state in wLCout/wRCout
    do iw=1,nw_recon
       where ((flagL(iL^S,iw) .eqv. .false.) .and. (flagR(iL^S,iw) .eqv. .false.))
          wLCout(iL^S,iw)=wLCin(iL^S,iw)
          wRCout(iL^S,iw)=wRCin(iL^S,iw)
       end where
    enddo
  
  end subroutine fix_limiter1

  subroutine WENO3limiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims, var
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    !> local
    double precision                :: f_array(ixI^S,1:nw,2), d_array(2)
    double precision                :: beta(ixI^S,1:nw,2),tau(ixI^S,1:nw)
    double precision                :: u1_coeff(2), u2_coeff(2)
    double precision                :: alpha_array(ixI^S,1:nw,2), alpha_sum(ixI^S,1:nw), flux(ixI^S,1:nw)
    double precision, dimension(ixI^S,1:nw)  :: wRCtmp, wLCtmp
    integer                         :: iLm^L, iLp^L, iLpp^L, i

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  
    iLm^L=iL^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    d_array(1:2) = (/ 1.0d0/3.0d0, 2.0d0/3.0d0 /)
    u1_coeff(1:2) = (/ -1.d0/2.d0, 3.d0/2.d0 /)
    u2_coeff(1:2) = (/ 1.d0/2.d0, 1.d0/2.d0 /)
    
    !> left side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * w(iLm^S,1:nw_recon) + u1_coeff(2) * w(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * w(iL^S,1:nw_recon)  + u2_coeff(2) * w(iLp^S,1:nw_recon)
  
    beta(iL^S,1:nw_recon,1) = (w(iL^S,1:nw_recon) - w(iLm^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = (w(iLp^S,1:nw_recon) - w(iL^S,1:nw_recon))**2
  
    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,2
        alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(2)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,2) - beta(iL^S,1:nw_recon,1))
      do i = 1,2
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nw_recon) / &
                                      (beta(iL^S,1:nw_recon,i) + dxdim**2)))
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    end select

    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,2
      flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) * alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do
  
    !> left value at right interface
    wLCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)
  
    !> right side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * w(iLpp^S,1:nw_recon) + u1_coeff(2) * w(iLp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * w(iLp^S,1:nw_recon)  + u2_coeff(2) * w(iL^S,1:nw_recon)
  
    beta(iL^S,1:nw_recon,1) = (w(iLpp^S,1:nw_recon) - w(iLp^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = (w(iLp^S,1:nw_recon) - w(iL^S,1:nw_recon))**2
  
    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,2
       alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
       alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(2)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,2) - beta(iL^S,1:nw_recon,1))
      do i = 1,2
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nw_recon) / &
                                      (beta(iL^S,1:nw_recon,i) + dxdim**2)))
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    end select

    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,2
       flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) * alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do
  
    !> right value at right interface
    wRCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)

    call fix_limiter(ixI^L,iL^L,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine WENO3limiter

  subroutine WENO5limiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims, var
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    !> local
    double precision                :: f_array(ixI^S,1:nw,3), d_array(3)
    double precision                :: beta(ixI^S,1:nw,3), beta_coeff(2)
    double precision                :: tau(ixI^S,1:nw), tmp(ixI^S,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:nw,3), alpha_sum(ixI^S,1:nw), flux(ixI^S,1:nw)
    double precision                :: lambda
    double precision, dimension(ixI^S,1:nw)  :: wRCtmp, wLCtmp
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L, i

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    lambda = dxdim**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
!   reconstruction variation
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
    u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
    u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
!   interpolation variation
!    d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
!    u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
!    u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
!    u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    
    !> left side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * w(iLmm^S,1:nw_recon) + u1_coeff(2) * w(iLm^S,1:nw_recon) + u1_coeff(3) * w(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * w(iLm^S,1:nw_recon)  + u2_coeff(2) * w(iL^S,1:nw_recon)  + u2_coeff(3) * w(iLp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * w(iL^S,1:nw_recon)   + u3_coeff(2) * w(iLp^S,1:nw_recon) + u3_coeff(3) * w(iLpp^S,1:nw_recon)  
  
    beta(iL^S,1:nw_recon,1) = beta_coeff(1) * (w(iLmm^S,1:nw_recon) + w(iL^S,1:nw_recon) - 2.0d0*w(iLm^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLmm^S,1:nw_recon) - 4.0d0 * w(iLm^S,1:nw_recon) + 3.0d0*w(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = beta_coeff(1) * (w(iLm^S,1:nw_recon) + w(iLp^S,1:nw_recon) - 2.0d0 * w(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLm^S,1:nw_recon) - w(iLp^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3) = beta_coeff(1) * (w(iL^S,1:nw_recon) + w(iLpp^S,1:nw_recon) - 2.0d0 * w(iLp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (3.0d0 * w(iL^S, 1:nw_recon) - 4.0d0 * w(iLp^S,1:nw_recon) + w(iLpp^S,1:nw_recon))**2
 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
      end do
    case(2)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3))
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nw_recon) / &
                                      (beta(iL^S,1:nw_recon,i) + weno_eps_machine))**2)
      end do
    case(3)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3))
      do i = 1,3
        tmp(iL^S,1:nw_recon) = (tau(iL^S,1:nw_recon) + weno_eps_machine) / (beta(iL^S,1:nw_recon,i) + weno_eps_machine)
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.0d0 + tmp(iL^S,1:nw_recon)**2 + lambda/tmp(iL^S,1:nw_recon))
      end do
    end select

    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    do i = 1,3
       alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
    end do
    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,3
       flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) &
                             *alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do
    wLCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)

    !> right side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * w(iLppp^S,1:nw_recon) + u1_coeff(2) * w(iLpp^S,1:nw_recon) + u1_coeff(3) * w(iLp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * w(iLpp^S,1:nw_recon)  + u2_coeff(2) * w(iLp^S,1:nw_recon)  + u2_coeff(3) * w(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * w(iLp^S,1:nw_recon)   + u3_coeff(2) * w(iL^S,1:nw_recon)   + u3_coeff(3) * w(iLm^S,1:nw_recon)  
  
    beta(iL^S,1:nw_recon,1) = beta_coeff(1) * (w(iLppp^S,1:nw_recon) + w(iLp^S,1:nw_recon) - 2.0d0*w(iLpp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLppp^S,1:nw_recon) - 4.0d0 * w(iLpp^S,1:nw_recon) + 3.0d0*w(iLp^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = beta_coeff(1) * (w(iLpp^S,1:nw_recon) + w(iL^S,1:nw_recon) - 2.0d0 * w(iLp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLpp^S,1:nw_recon) - w(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3) = beta_coeff(1) * (w(iLp^S,1:nw_recon) + w(iLm^S,1:nw_recon) - 2.0d0 * w(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (3.0d0 * w(iLp^S, 1:nw_recon) - 4.0d0 * w(iL^S,1:nw_recon) + w(iLm^S,1:nw_recon))**2
  
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
      end do
    case(2) 
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3))
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nw_recon) / &
                                      (beta(iL^S,1:nw_recon,i) + weno_eps_machine))**2)
      end do
    case(3)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3))
      do i = 1,3
        tmp(iL^S,1:nw_recon) = (tau(iL^S,1:nw_recon) + weno_eps_machine) / (beta(iL^S,1:nw_recon,i) + weno_eps_machine)
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.0d0 + tmp(iL^S,1:nw_recon)**2 + lambda/tmp(iL^S,1:nw_recon))
      end do
    end select
  
    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    ! do nothing, normal case
    do i = 1,3
       alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
    end do
    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,3
       flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) &
                             *alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do
    wRCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)
    
    call fix_limiter(ixI^L,iL^L,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine WENO5limiter

  subroutine TENO5ADlimiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    !> local
    double precision                :: f_array(ixI^S,1:nw,3), d_array(3)
    double precision                :: beta(ixI^S,1:nw,3), beta_coeff(2)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: gam_sum(ixI^S,1:nw),tau(ixI^S,1:nw),delta_sum(ixI^S,1:nw)
    double precision                :: gam(ixI^S,1:nw,3), kai(ixI^S,1:nw,3), delta(ixI^S,1:nw,3)
    double precision                :: flux(ixI^S,1:nw), kai1(ixI^S,1:nw,3), theta(ixI^S,1:nw)
    double precision, dimension(ixI^S,1:nw)  :: wRCtmp, wLCtmp
    integer                         :: marray(ixI^S,1:nw)
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L, i

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
    u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
    u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    
    !> left side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * w(iLmm^S,1:nw_recon) + u1_coeff(2) * w(iLm^S,1:nw_recon) + u1_coeff(3) * w(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * w(iLm^S,1:nw_recon)  + u2_coeff(2) * w(iL^S,1:nw_recon)  + u2_coeff(3) * w(iLp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * w(iL^S,1:nw_recon)   + u3_coeff(2) * w(iLp^S,1:nw_recon) + u3_coeff(3) * w(iLpp^S,1:nw_recon)  
  
    beta(iL^S,1:nw_recon,1) = beta_coeff(1) * (w(iLmm^S,1:nw_recon) + w(iL^S,1:nw_recon) - 2.0d0*w(iLm^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLmm^S,1:nw_recon) - 4.0d0 * w(iLm^S,1:nw_recon) + 3.0d0*w(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = beta_coeff(1) * (w(iLm^S,1:nw_recon) + w(iLp^S,1:nw_recon) - 2.0d0 * w(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLm^S,1:nw_recon) - w(iLp^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3) = beta_coeff(1) * (w(iL^S,1:nw_recon) + w(iLpp^S,1:nw_recon) - 2.0d0 * w(iLp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (3.0d0 * w(iL^S, 1:nw_recon) - 4.0d0 * w(iLp^S,1:nw_recon) + w(iLpp^S,1:nw_recon))**2
 
    gam_sum(iL^S,1:nw_recon) = 0.0d0
    tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3))
    do i = 1,3
      kai1(iL^S,1:nw_recon,i) = (tau(iL^S,1:nw_recon) / (beta(iL^S,1:nw_recon,i) + weno_eps_machine))
      gam(iL^S,1:nw_recon,i) = (1.d0 + kai1(iL^S,1:nw_recon,i))**2
      gam_sum(iL^S,1:nw_recon) = gam_sum(iL^S,1:nw_recon) + gam(iL^S,1:nw_recon,i)
    end do
    theta(iL^S,1:nw_recon) = one / (one + maxval(kai1(iL^S,1:nw_recon,1:3)/10.d0, dim=ndim+2))
    marray(iL^S,1:nw_recon)=-floor(4.d0 + theta(iL^S,1:nw_recon)*6.d0)
    do i = 1,3    
      kai(iL^S,1:nw_recon,i) = gam(iL^S,1:nw_recon,i) / gam_sum(iL^S,1:nw_recon)
      where(kai(iL^S,1:nw_recon,i) .lt. 10**marray(iL^S,1:nw_recon))
        delta(iL^S,1:nw_recon,i)=zero
      else where
        delta(iL^S,1:nw_recon,i)=one
      end where
    end do
    delta_sum=zero
    do i = 1,3
      delta_sum(iL^S,1:nw_recon)=delta_sum(iL^S,1:nw_recon)+delta(iL^S,1:nw_recon,i)*d_array(i)
    end do
    flux(iL^S,1:nw_recon)=0.0d0
    do i = 1,3
      flux(iL^S,1:nw_recon)=flux(iL^S,1:nw_recon)+f_array(iL^S,1:nw_recon,i)*delta(iL^S,1:nw_recon,i)*d_array(i)/(delta_sum(iL^S,1:nw_recon))
    end do
  
    !> left value at right interface
    wLCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)
  
    !> right side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * w(iLppp^S,1:nw_recon) + u1_coeff(2) * w(iLpp^S,1:nw_recon) + u1_coeff(3) * w(iLp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * w(iLpp^S,1:nw_recon)  + u2_coeff(2) * w(iLp^S,1:nw_recon)  + u2_coeff(3) * w(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * w(iLp^S,1:nw_recon)   + u3_coeff(2) * w(iL^S,1:nw_recon)   + u3_coeff(3) * w(iLm^S,1:nw_recon)  
  
    beta(iL^S,1:nw_recon,1) = beta_coeff(1) * (w(iLppp^S,1:nw_recon) + w(iLp^S,1:nw_recon) - 2.0d0*w(iLpp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLppp^S,1:nw_recon) - 4.0d0 * w(iLpp^S,1:nw_recon) + 3.0d0*w(iLp^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = beta_coeff(1) * (w(iLpp^S,1:nw_recon) + w(iL^S,1:nw_recon) - 2.0d0 * w(iLp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLpp^S,1:nw_recon) - w(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3) = beta_coeff(1) * (w(iLp^S,1:nw_recon) + w(iLm^S,1:nw_recon) - 2.0d0 * w(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (3.0d0 * w(iLp^S, 1:nw_recon) - 4.0d0 * w(iL^S,1:nw_recon) + w(iLm^S,1:nw_recon))**2
 

    gam_sum(iL^S,1:nw_recon)=0.0d0
    tau(iL^S,1:nw_recon)=abs(beta(iL^S,1:nw_recon,1)-beta(iL^S,1:nw_recon,3))
    do i=1,3
      kai1(iL^S,1:nw_recon,i)=(tau(iL^S,1:nw_recon)/(beta(iL^S,1:nw_recon,i)+weno_eps_machine))
      gam(iL^S,1:nw_recon,i)=(1.d0+kai1(iL^S,1:nw_recon,i))**6
      gam_sum(iL^S,1:nw_recon)=gam_sum(iL^S,1:nw_recon)+gam(iL^S,1:nw_recon,i)
    end do
    theta(iL^S,1:nw_recon)=one/(one+maxval(kai1(iL^S,1:nw_recon,1:3)/10.d0,dim=ndim+2))
    marray(iL^S,1:nw_recon)=-floor(4.d0+theta(iL^S,1:nw_recon)*6.d0)
    do i=1,3    
      kai(iL^S,1:nw_recon,i) = gam(iL^S,1:nw_recon,i)/gam_sum(iL^S,1:nw_recon)
      where(kai(iL^S,1:nw_recon,i) .lt. 10**marray(iL^S,1:nw_recon))
        delta(iL^S,1:nw_recon,i)=zero
      else where
        delta(iL^S,1:nw_recon,i)=one
      end where
    end do
    delta_sum=zero
    do i = 1,3
      delta_sum(iL^S,1:nw_recon)=delta_sum(iL^S,1:nw_recon)+delta(iL^S,1:nw_recon,i)*d_array(i)
    end do
    flux(iL^S,1:nw_recon)=0.0d0
    do i = 1,3
      flux(iL^S,1:nw_recon)=flux(iL^S,1:nw_recon)+f_array(iL^S,1:nw_recon,i)*delta(iL^S,1:nw_recon,i)*d_array(i)/(delta_sum(iL^S,1:nw_recon))
    end do

    !> right value at right interface
    wRCtmp(iL^S,1:nw_recon)=flux(iL^S,1:nw_recon)

    call fix_limiter(ixI^L,iL^L,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine TENO5ADlimiter

  subroutine WENO5NMlimiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in) :: ixI^L,iL^L,idims,var
    double precision, intent(in) :: dxdim
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    !> local
    double precision                :: f_array(ixI^S,1:nw,3), d_array(3)
    double precision                :: beta(ixI^S,1:nw,3), beta_coeff(2)
    double precision                :: tau(ixI^S,1:nw), tmp(ixI^S,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:nw,3), alpha_sum(ixI^S,1:nw), flux(ixI^S,1:nw)
    double precision                :: wc(ixI^S,1:nw), wd(ixI^S,1:nw), we(ixI^S,1:nw)
    double precision                :: lambda(ixI^S)
    double precision, dimension(ixI^S,1:nw)  :: wRCtmp, wLCtmp
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    integer                         :: iM^L, iMp^L
    integer                         :: i,j

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    iMmin^D=iLmmin^D;
    iMmax^D=iLpmax^D;
    iMp^L=iM^L+kr(idims,^D);
    lambda(iL^S)=block%dx(iL^S,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ -2.d0/3.d0, -1.d0/3.d0, 2.d0 /)
    u2_coeff(1:3) = (/ -1.d0/3.d0, 2.d0/3.d0, 2.d0/3.d0 /)
    u3_coeff(1:3) = (/ 2.d0/3.d0, 2.d0/3.d0, -1.d0/3.d0 /)
    do i = 1, nw_recon
      wc(iM^S,i) = (block%dx(iMp^S,idims) * w(iM^S,i) + block%dx(iM^S,idims) *  w(iMp^S,i)) / &
                   (block%dx(iMp^S,idims) + block%dx(iM^S,idims))
      wd(iL^S,i) = ((2.d0 * block%dx(iLm^S,idims) + block%dx(iLmm^S,idims)) * w(iLm^S,i) - block%dx(iLm^S,idims) * w(iLmm^S,i)) / &
                   (block%dx(iLmm^S,idims) + block%dx(iLm^S,idims))
      we(iL^S,i) = ((2.d0 * block%dx(iLpp^S,idims) + block%dx(iLppp^S,idims)) * w(iLpp^S,i) - block%dx(iLpp^S,idims) * w(iLppp^S,i)) / &
                   (block%dx(iLppp^S,idims) + block%dx(iLpp^S,idims))
    enddo
    !> left side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * wd(iL^S,1:nw_recon)   + u1_coeff(2) * wc(iLm^S,1:nw_recon)+ u1_coeff(3) * w(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * wc(iLm^S,1:nw_recon)  + u2_coeff(2) * w(iL^S,1:nw_recon)  + u2_coeff(3) * wc(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * wc(iL^S,1:nw_recon)   + u3_coeff(2) * w(iLp^S,1:nw_recon) + u3_coeff(3) * wc(iLp^S,1:nw_recon)  
  
    beta(iL^S,1:nw_recon,1) = beta_coeff(1) * (wc(iLm^S,1:nw_recon) - wd(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (2.d0 * w(iL^S,1:nw_recon) - wc(iLm^S,1:nw_recon) - wd(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = beta_coeff(1) * (wc(iLm^S,1:nw_recon) + wc(iL^S,1:nw_recon) - 2.0d0 * w(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (wc(iLm^S,1:nw_recon) - wc(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3) = beta_coeff(1) * (wc(iL^S,1:nw_recon) + wc(iLp^S,1:nw_recon) - 2.0d0 * w(iLp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (3.0d0 * wc(iL^S, 1:nw_recon) - 4.0d0 * w(iLp^S,1:nw_recon) + wc(iLp^S,1:nw_recon))**2
    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
         alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(4.d0 * beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
         alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(2)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3)) * 4.d0
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nw_recon) / &
                                      (4.d0 * beta(iL^S,1:nw_recon,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(3)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3)) * 4.d0
      do i=1,3
        do j=1,nw_recon
          tmp(iL^S,j) = (tau(iL^S,j) + weno_eps_machine) / (4.d0 * beta(iL^S,j,i) + weno_eps_machine)
          alpha_array(iL^S,j,i) = d_array(i) * (1.0d0 + tmp(iL^S,j)**2 + lambda(iL^S)/tmp(iL^S,j))
          alpha_sum(iL^S,j) = alpha_sum(iL^S,j) + alpha_array(iL^S,j,i)
        end do
      end do
    end select
    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,3
      flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) * alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do
    !> left value at right interface
    wLCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)

    !> right side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * we(iL^S,1:nw_recon)  + u1_coeff(2) * wc(iLp^S,1:nw_recon) + u1_coeff(3) * w(iLp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * wc(iLp^S,1:nw_recon) + u2_coeff(2) * w(iLp^S,1:nw_recon) + u2_coeff(3) * wc(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * wc(iL^S,1:nw_recon)  + u3_coeff(2) * w(iL^S,1:nw_recon)  + u3_coeff(3) * wc(iLm^S,1:nw_recon)  
    beta(iL^S,1:nw_recon,1) = beta_coeff(1) * (wc(iLp^S,1:nw_recon) - we(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (2.d0 * w(iLp^S,1:nw_recon) - wc(iLp^S,1:nw_recon) - we(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = beta_coeff(1) * (wc(iLp^S,1:nw_recon) + wc(iL^S,1:nw_recon) - 2.0d0 * w(iLp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (wc(iLp^S,1:nw_recon) - wc(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3) = beta_coeff(1) * (wc(iL^S,1:nw_recon) + wc(iLm^S,1:nw_recon) - 2.0d0 * w(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (3.0d0 * wc(iL^S, 1:nw_recon) - 4.0d0 * w(iL^S,1:nw_recon) + wc(iLm^S,1:nw_recon))**2
    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(4.d0 * beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(2) 
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3)) * 4.d0
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nw_recon) / &
                                      (4.d0 * beta(iL^S,1:nw_recon,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(3)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3)) * 4.d0
      do i = 1,3
        do j = 1,nw_recon
          tmp(iL^S,j) = (tau(iL^S,j) + weno_eps_machine) / (4.d0 * beta(iL^S,j,i) + weno_eps_machine)
          alpha_array(iL^S,j,i) = d_array(i) * (1.0d0 + tmp(iL^S,j)**2 + lambda(iL^S)/tmp(iL^S,j))
          alpha_sum(iL^S,j) = alpha_sum(iL^S,j) + alpha_array(iL^S,j,i)
        end do
      end do
    end select
    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,3
      flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) * alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do
    !> right value at right interface
    wRCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)

    call fix_limiter(ixI^L,iL^L,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine WENO5NMlimiter

  subroutine WENO5limiterL(ixI^L,iL^L,idims,w,wLC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims
    integer, intent(in)             :: var
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wLC(ixI^S,1:nw) 
    !> local
    double precision                :: f_array(ixI^S,1:nw,3), d_array(3)
    double precision                :: beta(ixI^S,1:nw,3), beta_coeff(2)
    double precision                :: tau(ixI^S,1:nw), tmp(ixI^S,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:nw,3), alpha_sum(ixI^S,1:nw), flux(ixI^S,1:nw)
    double precision                :: lambda(ixI^S)
    double precision, dimension(ixI^S,1:nw)  :: wLCtmp
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    integer                         :: i,j

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    lambda=block%dx(iL^S,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
!   reconstruction variation
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
    u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
    u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
!   interpolation variation
!    d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
!    u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
!    u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
!    u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    
    !> left side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * w(iLmm^S,1:nw_recon) + u1_coeff(2) * w(iLm^S,1:nw_recon) + u1_coeff(3) * w(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * w(iLm^S,1:nw_recon)  + u2_coeff(2) * w(iL^S,1:nw_recon)  + u2_coeff(3) * w(iLp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * w(iL^S,1:nw_recon)   + u3_coeff(2) * w(iLp^S,1:nw_recon) + u3_coeff(3) * w(iLpp^S,1:nw_recon)  
  
    beta(iL^S,1:nw_recon,1) = beta_coeff(1) * (w(iLmm^S,1:nw_recon) + w(iL^S,1:nw_recon) - 2.0d0*w(iLm^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLmm^S,1:nw_recon) - 4.0d0 * w(iLm^S,1:nw_recon) + 3.0d0*w(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = beta_coeff(1) * (w(iLm^S,1:nw_recon) + w(iLp^S,1:nw_recon) - 2.0d0 * w(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLm^S,1:nw_recon) - w(iLp^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3) = beta_coeff(1) * (w(iL^S,1:nw_recon) + w(iLpp^S,1:nw_recon) - 2.0d0 * w(iLp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (3.0d0 * w(iL^S, 1:nw_recon) - 4.0d0 * w(iLp^S,1:nw_recon) + w(iLpp^S,1:nw_recon))**2
 
    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz
    case(1)
      do i = 1,3
         alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
         alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(2)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3))
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nw_recon) / &
                                      (beta(iL^S,1:nw_recon,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(3)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3)) * 4.d0
      do i=1,3
        do j=1,nw_recon
          tmp(iL^S,j) = (tau(iL^S,j) + weno_eps_machine) / (4.d0 * beta(iL^S,j,i) + weno_eps_machine)
          alpha_array(iL^S,j,i) = d_array(i) * (1.0d0 + tmp(iL^S,j)**2 + lambda(iL^S)/tmp(iL^S,j))
          alpha_sum(iL^S,j) = alpha_sum(iL^S,j) + alpha_array(iL^S,j,i)
        end do
      end do
    end select
    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,3
      flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) * alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do
  
    !> left value at right interface
    wLCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)
  
    call fix_onelimiter(ixI^L,iL^L,wLCtmp,wLC)

  end subroutine WENO5limiterL

  subroutine WENO5limiterR(ixI^L,iL^L,idims,w,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims
    integer, intent(in)             :: var
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw)
    !> local
    double precision                :: f_array(ixI^S,1:nw,3), d_array(3)
    double precision                :: beta(ixI^S,1:nw,3), beta_coeff(2)
    double precision                :: tau(ixI^S,1:nw), tmp(ixI^S,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:nw,3), alpha_sum(ixI^S,1:nw), flux(ixI^S,1:nw)
    double precision                :: lambda(ixI^S)
    double precision, dimension(ixI^S,1:nw)  :: wRCtmp
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    integer                         :: i,j

    iLm^L=iL^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    lambda=block%dx(iL^S,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
!   reconstruction variation
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
    u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
    u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
!   interpolation variation
!    d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
!    u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
!    u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
!    u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    
    !> right side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * w(iLppp^S,1:nw_recon) + u1_coeff(2) * w(iLpp^S,1:nw_recon) + u1_coeff(3) * w(iLp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * w(iLpp^S,1:nw_recon)  + u2_coeff(2) * w(iLp^S,1:nw_recon)  + u2_coeff(3) * w(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * w(iLp^S,1:nw_recon)   + u3_coeff(2) * w(iL^S,1:nw_recon)   + u3_coeff(3) * w(iLm^S,1:nw_recon)  
  
    beta(iL^S,1:nw_recon,1) = beta_coeff(1) * (w(iLppp^S,1:nw_recon) + w(iLp^S,1:nw_recon) - 2.0d0*w(iLpp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLppp^S,1:nw_recon) - 4.0d0 * w(iLpp^S,1:nw_recon) + 3.0d0*w(iLp^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = beta_coeff(1) * (w(iLpp^S,1:nw_recon) + w(iL^S,1:nw_recon) - 2.0d0 * w(iLp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (w(iLpp^S,1:nw_recon) - w(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3) = beta_coeff(1) * (w(iLp^S,1:nw_recon) + w(iLm^S,1:nw_recon) - 2.0d0 * w(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (3.0d0 * w(iLp^S, 1:nw_recon) - 4.0d0 * w(iL^S,1:nw_recon) + w(iLm^S,1:nw_recon))**2
  
    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(2) 
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3))
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nw_recon) / &
                                      (beta(iL^S,1:nw_recon,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(3)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3)) * 4.d0
      do i = 1,3
        do j = 1,nw_recon
          tmp(iL^S,j) = (tau(iL^S,j) + weno_eps_machine) / (4.d0 * beta(iL^S,j,i) + weno_eps_machine)
          alpha_array(iL^S,j,i) = d_array(i) * (1.0d0 + tmp(iL^S,j)**2 + lambda(iL^S)/tmp(iL^S,j))
          alpha_sum(iL^S,j) = alpha_sum(iL^S,j) + alpha_array(iL^S,j,i)
        end do
      end do
    end select
    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,3
      flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) * alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do
  
    !> right value at right interface
    wRCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)

    call fix_onelimiter(ixI^L,iL^L,wRCtmp,wRC)

  end subroutine WENO5limiterR

  subroutine WENO5NMlimiterL(ixI^L,iL^L,idims,w,wLC,var)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,iL^L,idims,var
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wLC(ixI^S,1:nw) 
    !> local
    double precision                :: f_array(ixI^S,1:nw,3), d_array(3)
    double precision                :: beta(ixI^S,1:nw,3), beta_coeff(2)
    double precision                :: tau(ixI^S,1:nw), tmp(ixI^S,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:nw,3), alpha_sum(ixI^S,1:nw), flux(ixI^S,1:nw)
    double precision                :: wc(ixI^S,1:nw), wd(ixI^S,1:nw)
    double precision                :: lambda(ixI^S)
    double precision, dimension(ixI^S,1:nw)  :: wLCtmp
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L
    integer                         :: iM^L, iMp^L
    integer                         :: i,j

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iMmin^D=iLmmin^D;
    iMmax^D=iLpmax^D;
    iMp^L=iM^L+kr(idims,^D);
    lambda(iL^S)=block%dx(iL^S,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ -2.d0/3.d0, -1.d0/3.d0, 2.d0 /)
    u2_coeff(1:3) = (/ -1.d0/3.d0, 2.d0/3.d0, 2.d0/3.d0 /)
    u3_coeff(1:3) = (/ 2.d0/3.d0, 2.d0/3.d0, -1.d0/3.d0 /)
    do i = 1, nw_recon
      wc(iM^S,i) = (block%dx(iMp^S,idims) * w(iM^S,i) + block%dx(iM^S,idims) * w(iMp^S,i)) / &
                   (block%dx(iMp^S,idims) + block%dx(iM^S,idims))
      wd(iL^S,i) = ((2.d0 * block%dx(iLm^S,idims) + block%dx(iLmm^S,idims)) * w(iLm^S,i) - block%dx(iLm^S,idims) * w(iLmm^S,i)) / &
                   (block%dx(iLmm^S,idims) + block%dx(iLm^S,idims))
    enddo
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * wd(iL^S,1:nw_recon)   + u1_coeff(2) * wc(iLm^S,1:nw_recon)+ u1_coeff(3) * w(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * wc(iLm^S,1:nw_recon)  + u2_coeff(2) * w(iL^S,1:nw_recon)  + u2_coeff(3) * wc(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * wc(iL^S,1:nw_recon)   + u3_coeff(2) * w(iLp^S,1:nw_recon) + u3_coeff(3) * wc(iLp^S,1:nw_recon)  
    beta(iL^S,1:nw_recon,1) = beta_coeff(1) * (wc(iLm^S,1:nw_recon) - wd(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (2.d0 * w(iL^S,1:nw_recon) - wc(iLm^S,1:nw_recon) - wd(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = beta_coeff(1) * (wc(iLm^S,1:nw_recon) + wc(iL^S,1:nw_recon) - 2.0d0 * w(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (wc(iLm^S,1:nw_recon) - wc(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3) = beta_coeff(1) * (wc(iL^S,1:nw_recon) + wc(iLp^S,1:nw_recon) - 2.0d0 * w(iLp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (3.0d0 * wc(iL^S, 1:nw_recon) - 4.0d0 * w(iLp^S,1:nw_recon) + wc(iLp^S,1:nw_recon))**2
    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
         alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(4.d0 * beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
         alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(2)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3)) * 4.d0
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nw_recon) / &
                                      (4.d0 * beta(iL^S,1:nw_recon,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(3)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3)) * 4.d0
      do i=1,3
        do j=1,nw_recon
          tmp(iL^S,j) = (tau(iL^S,j) + weno_eps_machine) / (4.d0 * beta(iL^S,j,i) + weno_eps_machine)
          alpha_array(iL^S,j,i) = d_array(i) * (1.0d0 + tmp(iL^S,j)**2 + lambda(iL^S)/tmp(iL^S,j))
          alpha_sum(iL^S,j) = alpha_sum(iL^S,j) + alpha_array(iL^S,j,i)
        end do
      end do
    end select
    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,3
      flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) * alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do
    !> left value at right interface
    wLCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)

    call fix_onelimiter(ixI^L,iL^L,wLCtmp,wLC)

  end subroutine WENO5NMlimiterL

  subroutine WENO5NMlimiterR(ixI^L,iL^L,idims,w,wRC,var)
    use mod_global_parameters
  
    integer, intent(in) :: ixI^L,iL^L,idims,var
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw)
    !> local
    double precision                :: f_array(ixI^S,1:nw,3), d_array(3)
    double precision                :: beta(ixI^S,1:nw,3), beta_coeff(2)
    double precision                :: tau(ixI^S,1:nw), tmp(ixI^S,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:nw,3), alpha_sum(ixI^S,1:nw), flux(ixI^S,1:nw)
    double precision                :: wc(ixI^S,1:nw), we(ixI^S,1:nw)
    double precision                :: lambda(ixI^S)
    double precision, dimension(ixI^S,1:nw)  :: wRCtmp
    integer                         :: iLm^L,iLp^L,iLpp^L,iLppp^L
    integer                         :: iM^L, iMp^L
    integer                         :: i,j

    iLm^L=iL^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    iMmin^D=iLmmin^D;
    iMmax^D=iLpmax^D;
    iMp^L=iM^L+kr(idims,^D);
    lambda(iL^S)=block%dx(iL^S,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ -2.d0/3.d0, -1.d0/3.d0, 2.d0 /)
    u2_coeff(1:3) = (/ -1.d0/3.d0, 2.d0/3.d0, 2.d0/3.d0 /)
    u3_coeff(1:3) = (/ 2.d0/3.d0, 2.d0/3.d0, -1.d0/3.d0 /)
    do i = 1, nw_recon
      wc(iM^S,i) = (block%dx(iMp^S,idims) * w(iM^S,i) + block%dx(iM^S,idims) * w(iMp^S,i)) / &
                   (block%dx(iMp^S,idims) + block%dx(iM^S,idims))
      we(iL^S,i) = ((2.d0 * block%dx(iLpp^S,idims) + block%dx(iLppp^S,idims)) * w(iLpp^S,i) - block%dx(iLpp^S,idims) * w(iLppp^S,i)) / &
                   (block%dx(iLppp^S,idims) + block%dx(iLpp^S,idims))
    enddo
    !> right side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * we(iL^S,1:nw_recon)  + u1_coeff(2) * wc(iLp^S,1:nw_recon) + u1_coeff(3) * w(iLp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * wc(iLp^S,1:nw_recon) + u2_coeff(2) * w(iLp^S,1:nw_recon) + u2_coeff(3) * wc(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * wc(iL^S,1:nw_recon)  + u3_coeff(2) * w(iL^S,1:nw_recon)  + u3_coeff(3) * wc(iLm^S,1:nw_recon)  
    beta(iL^S,1:nw_recon,1) = beta_coeff(1) * (wc(iLp^S,1:nw_recon) - we(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (2.d0 * w(iLp^S,1:nw_recon) - wc(iLp^S,1:nw_recon) - we(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2) = beta_coeff(1) * (wc(iLp^S,1:nw_recon) + wc(iL^S,1:nw_recon) - 2.0d0 * w(iLp^S,1:nw_recon))**2 &
         + beta_coeff(2) * (wc(iLp^S,1:nw_recon) - wc(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3) = beta_coeff(1) * (wc(iL^S,1:nw_recon) + wc(iLm^S,1:nw_recon) - 2.0d0 * w(iL^S,1:nw_recon))**2 &
         + beta_coeff(2) * (3.0d0 * wc(iL^S, 1:nw_recon) - 4.0d0 * w(iL^S,1:nw_recon) + wc(iLm^S,1:nw_recon))**2
    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(4.d0 * beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(2) 
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3)) * 4.d0
      do i = 1,3
        alpha_array(iL^S,1:nw_recon,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nw_recon) / &
                                      (4.d0 * beta(iL^S,1:nw_recon,i) + weno_eps_machine))**2)
        alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
      end do
    case(3)
      tau(iL^S,1:nw_recon) = abs(beta(iL^S,1:nw_recon,1) - beta(iL^S,1:nw_recon,3)) * 4.d0
      do i = 1,3
        do j = 1,nw_recon
          tmp(iL^S,j) = (tau(iL^S,j) + weno_eps_machine) / (4.d0 * beta(iL^S,j,i) + weno_eps_machine)
          alpha_array(iL^S,j,i) = d_array(i) * (1.0d0 + tmp(iL^S,j)**2 + lambda(iL^S)/tmp(iL^S,j))
          alpha_sum(iL^S,j) = alpha_sum(iL^S,j) + alpha_array(iL^S,j,i)
        end do
      end do
    end select
    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,3
      flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) * alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do
    !> right value at right interface
    wRCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)

    call fix_onelimiter(ixI^L,iL^L,wRCtmp,wRC)

  end subroutine WENO5NMlimiterR

  subroutine WENO5CU6limiter(ixI^L,iL^L,idims,w,wLC,wRC)
    use mod_global_parameters
  
    integer, intent(in) :: ixI^L, iL^L, idims
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    !> local
    double precision :: f_array(ixI^S,1:nw,3), d_array(3)
    double precision :: beta(ixI^S,1:nw,3), beta_coeff(2)
    double precision :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision :: alpha_array(ixI^S,1:nw,3), alpha_sum(ixI^S,1:nw)
    double precision :: theta2(ixI^S,1:nw)
    double precision, parameter :: theta_limit=0.7d0
    double precision, dimension(ixI^S,1:nw)  :: wRCtmp, wLCtmp
    integer :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    integer :: i

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);

    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
    u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
    u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    
    !> left side
    beta(iL^S,1:nw_recon,1)=beta_coeff(1)*(w(iLmm^S,1:nw_recon)+w(iL^S,1:nw_recon)-2.0d0*w(iLm^S,1:nw_recon))**2&
        +beta_coeff(2)*(w(iLmm^S,1:nw_recon)-4.0d0*w(iLm^S,1:nw_recon)+3.0d0*w(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2)=beta_coeff(1)*(w(iLm^S,1:nw_recon)+w(iLp^S,1:nw_recon)-2.0d0*w(iL^S,1:nw_recon))**2&
        +beta_coeff(2)*(w(iLm^S,1:nw_recon)-w(iLp^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3)=beta_coeff(1)*(w(iL^S,1:nw_recon)+w(iLpp^S,1:nw_recon)-2.0d0*w(iLp^S,1:nw_recon))**2&
        +beta_coeff(2)*(3.0d0*w(iL^S,1:nw_recon)-4.0d0*w(iLp^S,1:nw_recon)+w(iLpp^S,1:nw_recon))**2
    alpha_sum(iL^S,1:nw_recon)=zero
    do i=1,3
      alpha_array(iL^S,1:nw_recon,i)=d_array(i)/(beta(iL^S,1:nw_recon,i)+weno_eps_machine)**2
      alpha_sum(iL^S,1:nw_recon)=alpha_sum(iL^S,1:nw_recon)+alpha_array(iL^S,1:nw_recon,i)
    end do
    do i=1,3
      alpha_array(iL^S,1:nw_recon,i)=alpha_array(iL^S,1:nw_recon,i)/alpha_sum(iL^S,1:nw_recon)
    end do
    theta2(iL^S,1:nw_recon)=((alpha_array(iL^S,1:nw_recon,1)/d_array(1)-1.d0)**2&
                          +(alpha_array(iL^S,1:nw_recon,2)/d_array(2)-1.d0)**2&
                          +(alpha_array(iL^S,1:nw_recon,3)/d_array(3)-1.d0)**2)/83.d0
    where(theta2(iL^S,1:nw_recon) .gt. theta_limit)
      f_array(iL^S,1:nw_recon,1)=u1_coeff(1)*w(iLmm^S,1:nw_recon)+u1_coeff(2)*w(iLm^S,1:nw_recon)+u1_coeff(3)*w(iL^S,1:nw_recon)
      f_array(iL^S,1:nw_recon,2)=u2_coeff(1)*w(iLm^S,1:nw_recon)+u2_coeff(2)*w(iL^S,1:nw_recon)+u2_coeff(3)*w(iLp^S,1:nw_recon)
      f_array(iL^S,1:nw_recon,3)=u3_coeff(1)*w(iL^S,1:nw_recon)+u3_coeff(2)*w(iLp^S,1:nw_recon)+u3_coeff(3)*w(iLpp^S,1:nw_recon)  
      wLCtmp(iL^S,1:nw_recon)=f_array(iL^S,1:nw_recon,1)*alpha_array(iL^S,1:nw_recon,1)&
                        +f_array(iL^S,1:nw_recon,2)*alpha_array(iL^S,1:nw_recon,2)&
                        +f_array(iL^S,1:nw_recon,3)*alpha_array(iL^S,1:nw_recon,3)
    else where
      wLCtmp(iL^S,1:nw_recon)=1.d0/60.d0*(w(iLmm^S,1:nw_recon)-8.d0*w(iLm^S,1:nw_recon)+37.d0*w(iL^S,1:nw_recon)&
                         +37.d0*w(iLp^S,1:nw_recon)-8.d0*w(iLpp^S,1:nw_recon)+w(iLppp^S,1:nw_recon))
    end where
  
    !> right side
    beta(iL^S,1:nw_recon,1)=beta_coeff(1)*(w(iLppp^S,1:nw_recon)+w(iLp^S,1:nw_recon)-2.0d0*w(iLpp^S,1:nw_recon))**2&
         +beta_coeff(2)*(w(iLppp^S,1:nw_recon)-4.0d0*w(iLpp^S,1:nw_recon)+3.0d0*w(iLp^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,2)=beta_coeff(1)*(w(iLpp^S,1:nw_recon)+w(iL^S,1:nw_recon)-2.0d0*w(iLp^S,1:nw_recon))**2&
         +beta_coeff(2)*(w(iLpp^S,1:nw_recon)-w(iL^S,1:nw_recon))**2
    beta(iL^S,1:nw_recon,3)=beta_coeff(1)*(w(iLp^S,1:nw_recon)+w(iLm^S,1:nw_recon)-2.0d0*w(iL^S,1:nw_recon))**2&
         +beta_coeff(2)*(3.0d0*w(iLp^S,1:nw_recon)-4.0d0*w(iL^S,1:nw_recon)+w(iLm^S,1:nw_recon))**2
    alpha_sum(iL^S,1:nw_recon)=zero
    do i=1,3
      alpha_array(iL^S,1:nw_recon,i)=d_array(i)/(beta(iL^S,1:nw_recon,i)+weno_eps_machine)**2
      alpha_sum(iL^S,1:nw_recon)=alpha_sum(iL^S,1:nw_recon)+alpha_array(iL^S,1:nw_recon,i)
    end do
    do i=1,3
      alpha_array(iL^S,1:nw_recon,i)=alpha_array(iL^S,1:nw_recon,i)/alpha_sum(iL^S,1:nw_recon)
    end do
    theta2(iL^S,1:nw_recon)=((alpha_array(iL^S,1:nw_recon,1)/d_array(1)-1.d0)**2&
                          +(alpha_array(iL^S,1:nw_recon,2)/d_array(2)-1.d0)**2&
                          +(alpha_array(iL^S,1:nw_recon,3)/d_array(3)-1.d0)**2)/83.d0
    where(theta2(iL^S,1:nw_recon) .gt. theta_limit)
      f_array(iL^S,1:nw_recon,1)=u1_coeff(1)*w(iLppp^S,1:nw_recon)+u1_coeff(2)*w(iLpp^S,1:nw_recon)+u1_coeff(3)*w(iLp^S,1:nw_recon)
      f_array(iL^S,1:nw_recon,2)=u2_coeff(1)*w(iLpp^S,1:nw_recon)+u2_coeff(2)*w(iLp^S,1:nw_recon)+u2_coeff(3)*w(iL^S,1:nw_recon)
      f_array(iL^S,1:nw_recon,3)=u3_coeff(1)*w(iLp^S,1:nw_recon)+u3_coeff(2)*w(iL^S,1:nw_recon)+u3_coeff(3)*w(iLm^S,1:nw_recon)  
      wRCtmp(iL^S,1:nw_recon)=f_array(iL^S,1:nw_recon,1)*alpha_array(iL^S,1:nw_recon,1)&
                        +f_array(iL^S,1:nw_recon,2)*alpha_array(iL^S,1:nw_recon,2)&
                        +f_array(iL^S,1:nw_recon,3)*alpha_array(iL^S,1:nw_recon,3)
    else where
      wRCtmp(iL^S,1:nw_recon)=1.d0/60.d0*(w(iLppp^S,1:nw_recon)-8.d0*w(iLpp^S,1:nw_recon)+37.d0*w(iLp^S,1:nw_recon)&
                        +37.d0*w(iL^S,1:nw_recon)-8.d0*w(iLm^S,1:nw_recon)+w(iLmm^S,1:nw_recon))
    end where

    call fix_limiter(ixI^L,iL^L,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine WENO5CU6limiter

  subroutine WENO7limiter(ixI^L,iL^L,idims,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims, var
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    !> local
    double precision, dimension(4)  :: d_array, u1_coeff, u2_coeff, u3_coeff, u4_coeff
    double precision, dimension(ixI^S,1:nw,4)  :: f_array, beta, alpha_array
    double precision, dimension(ixI^S)         :: a, b, c, tmp, tmp2, tmp3
    double precision, dimension(ixI^S,1:nw)    :: alpha_sum, d, dm4
    double precision, dimension(ixI^S,1:nw)    :: flux, flux_min, flux_max, flux_ul, flux_md, flux_lc
    double precision, parameter     :: mpalpha = 2.d0, mpbeta = 4.d0
    double precision, dimension(ixI^S,1:nw)  :: wRCtmp, wLCtmp
    integer                         :: iLm^L, iLmm^L, iLmmm^L
    integer                         :: iLp^L, iLpp^L, iLppp^L, iLpppp^L
    integer                         :: id^L, idp^L, idpp^L, idm^L, ie^L, iem^L, iep^L, iepp^L
    integer                         :: i,iw

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLmmm^L=iLmm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    iLpppp^L=iLppp^L+kr(idims,^D);

    d_array(1:4) = (/ 1.d0/35.d0, 12.d0/35.d0, 18.d0/35.d0, 4.d0/35.d0 /)
    u1_coeff(1:4) = (/ -1.d0/4.d0, 13.d0/12.d0, -23.d0/12.d0, 25.d0/12.d0 /)
    u2_coeff(1:4) = (/ 1.d0/12.d0, -5.d0/12.d0, 13.d0/12.d0, 1.d0/4.d0 /)
    u3_coeff(1:4) = (/ -1.d0/12.d0, 7.d0/12.d0, 7.d0/12.d0, -1.d0/12.d0 /)
    u4_coeff(1:4) = (/ 1.d0/4.d0, 13.d0/12.d0, -5.d0/12.d0, 1.d0/12.d0 /)
    
    !> left side
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * w(iLmmm^S,1:nw_recon) + u1_coeff(2) * w(iLmm^S,1:nw_recon) + u1_coeff(3) * w(iLm^S,1:nw_recon) &
                             + u1_coeff(4) * w(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * w(iLmm^S,1:nw_recon)  + u2_coeff(2) * w(iLm^S,1:nw_recon)  + u2_coeff(3) * w(iL^S,1:nw_recon)  &
                             + u2_coeff(4) * w(iLp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * w(iLm^S,1:nw_recon)   + u3_coeff(2) * w(iL^S,1:nw_recon)   + u3_coeff(3) * w(iLp^S,1:nw_recon)   &
                             + u3_coeff(4) * w(iLpp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,4) = u4_coeff(1) * w(iL^S,1:nw_recon)    + u4_coeff(2) * w(iLp^S,1:nw_recon)    + u4_coeff(3) * w(iLpp^S,1:nw_recon)  &
                             + u4_coeff(4) * w(iLppp^S,1:nw_recon)
  
    beta(iL^S,1:nw_recon,1) = w(iLmmm^S,1:nw_recon) * (547.d0 * w(iLmmm^S,1:nw_recon) - 3882.d0 * w(iLmm^S,1:nw_recon) + 4642.d0 * w(iLm^S,1:nw_recon) &
                          - 1854.d0 * w(iL^S,1:nw_recon)) &
                          + w(iLmm^S,1:nw_recon) * (7043.d0 * w(iLmm^S,1:nw_recon) - 17246.d0 * w(iLm^S,1:nw_recon) + 7042.d0 * w(iL^S,1:nw_recon)) &
                          + w(iLm^S,1:nw_recon) * (11003.d0 * w(iLm^S,1:nw_recon) - 9402.d0 * w(iL^S,1:nw_recon)) + 2107.d0 * w(iL^S,1:nw_recon)**2
    beta(iL^S,1:nw_recon,2) = w(iLmm^S,1:nw_recon) * (267.d0 * w(iLmm^S,1:nw_recon) - 1642.d0 * w(iLm^S,1:nw_recon) + 1602.d0 * w(iL^S,1:nw_recon) &
                          - 494.d0 * w(iLp^S,1:nw_recon))  &
                          + w(iLm^S,1:nw_recon) * (2843.d0 * w(iLm^S,1:nw_recon) - 5966.d0 * w(iL^S,1:nw_recon) + 1922.d0 * w(iLp^S,1:nw_recon)) &
                          + w(iL^S,1:nw_recon) * (3443.d0 * w(iL^S,1:nw_recon) - 2522.d0 * w(iLp^S,1:nw_recon)) + 547.d0 * w(iLp^S,1:nw_recon) ** 2
    beta(iL^S,1:nw_recon,3) = w(iLm^S,1:nw_recon) * (547.d0 * w(iLm^S,1:nw_recon) - 2522.d0 * w(iL^S,1:nw_recon) + 1922.d0 * w(iLp^S,1:nw_recon) &
                          - 494.d0 * w(iLpp^S,1:nw_recon))  &
                          + w(iL^S,1:nw_recon) * (3443.d0 * w(iL^S,1:nw_recon) - 5966.d0 * w(iLp^S,1:nw_recon) + 1602.d0 * w(iLpp^S,1:nw_recon)) &
                          + w(iLp^S,1:nw_recon) * (2843.d0 * w(iLp^S,1:nw_recon) - 1642.d0 * w(iLpp^S,1:nw_recon)) + 267.d0 * w(iLpp^S,1:nw_recon) ** 2
    beta(iL^S,1:nw_recon,4) = w(iL^S,1:nw_recon) * (2107.d0 * w(iL^S,1:nw_recon) - 9402.d0 * w(iLp^S,1:nw_recon) + 7042.d0 * w(iLpp^S,1:nw_recon) &
                          - 1854.d0 * w(iLppp^S,1:nw_recon))  &
                          + w(iLp^S,1:nw_recon) * (11003.d0 * w(iLp^S,1:nw_recon) - 17246.d0 * w(iLpp^S,1:nw_recon) + 4642.d0 * w(iLppp^S,1:nw_recon)) &
                          + w(iLpp^S,1:nw_recon) * (7043.d0 * w(iLpp^S,1:nw_recon) - 3882.d0 * w(iLppp^S,1:nw_recon)) &
                          + 547.d0 * w(iLppp^S,1:nw_recon) ** 2

    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    do i = 1,4
       alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
       alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
    end do
    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,4
       flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) * alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do
    
    select case(var)
    ! case1 for wenojs, case2 for mpweno
    case(1) 
      wLCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)
    case(2)
      idmax^D=iLmax^D; idmin^D=iLmin^D-kr(idims,^D);
      idm^L=id^L-kr(idims,^D);
      idp^L=id^L+kr(idims,^D);
  
      iemax^D=idmax^D+kr(idims,^D); iemin^D=idmin^D;
      iem^L=ie^L-kr(idims,^D);
      iep^L=ie^L+kr(idims,^D);
  
      d(ie^S,1:nw_recon) = w(iep^S,1:nw_recon)-2.0d0*w(ie^S,1:nw_recon)+w(iem^S,1:nw_recon)
  
      do iw=1,nw_recon
         a(id^S) = 4.0d0*d(id^S,iw)-d(idp^S,iw)
         b(id^S) = 4.0d0*d(idp^S,iw)-d(id^S,iw)
         call minmod(ixI^L,id^L,a,b,tmp)
         a(id^S) = d(id^S,iw)
         b(id^S) = d(idp^S,iw)
         call minmod(ixI^L,id^L,a,b,tmp2)
         call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
         dm4(id^S,iw) = tmp3(id^S)
      end do

      flux_ul(iL^S,1:nw_recon) = w(iL^S,1:nw_recon) + mpalpha * (w(iL^S,1:nw_recon) - w(iLm^S,1:nw_recon))
      flux_md(iL^S,1:nw_recon) = half * (w(iL^S,1:nw_recon) + w(iLp^S,1:nw_recon) - dm4(iL^S,1:nw_recon))
      flux_lc(iL^S,1:nw_recon) = half * (3.d0 * w(iL^S,1:nw_recon) - w(iLm^S,1:nw_recon)) + mpbeta / 3.d0 * dm4(iLm^S,1:nw_recon)
    
      flux_min(iL^S,1:nw_recon) = max(min(w(iL^S,1:nw_recon), w(iLp^S,1:nw_recon), flux_md(iL^S,1:nw_recon)), &
                                min(w(iL^S,1:nw_recon), flux_ul(iL^S,1:nw_recon),flux_lc(iL^S,1:nw_recon)))
  
      flux_max(iL^S,1:nw_recon) = min(max(w(iL^S,1:nw_recon), w(iLp^S,1:nw_recon), flux_md(iL^S,1:nw_recon)), &
                                max(w(iL^S,1:nw_recon), flux_ul(iL^S,1:nw_recon),flux_lc(iL^S,1:nw_recon)))
      do iw=1,nw_recon
        a(iL^S) = flux(iL^S,iw)
        b(iL^S) = flux_min(iL^S,iw)
        c(iL^S) = flux_max(iL^S,iw)
        call median(ixI^L, iL^L, a, b, c, tmp) 
        wLCtmp(iL^S,iw) = tmp(iL^S)
      end do
    end select

    !> right side
    !>> mmm -> pppp
    !>> mm  -> ppp
    !>> m   -> pp
    !>> 0   -> p
    !>> p   -> 0
    !>> pp  -> m
    !>> ppp -> mm
    f_array(iL^S,1:nw_recon,1) = u1_coeff(1) * w(iLpppp^S,1:nw_recon) + u1_coeff(2) * w(iLppp^S,1:nw_recon) + u1_coeff(3) * w(iLpp^S,1:nw_recon) &
                             + u1_coeff(4) * w(iLp^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,2) = u2_coeff(1) * w(iLppp^S,1:nw_recon)  + u2_coeff(2) * w(iLpp^S,1:nw_recon)  + u2_coeff(3) * w(iLp^S,1:nw_recon)  &
                             + u2_coeff(4) * w(iL^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,3) = u3_coeff(1) * w(iLpp^S,1:nw_recon)   + u3_coeff(2) * w(iLp^S,1:nw_recon)   + u3_coeff(3) * w(iL^S,1:nw_recon)   &
                             + u3_coeff(4) * w(iLm^S,1:nw_recon)
    f_array(iL^S,1:nw_recon,4) = u4_coeff(1) * w(iLp^S,1:nw_recon)    + u4_coeff(2) * w(iL^S,1:nw_recon)    + u4_coeff(3) * w(iLm^S,1:nw_recon)  &
                             + u4_coeff(4) * w(iLmm^S,1:nw_recon)

    beta(iL^S,1:nw_recon,1) = w(iLpppp^S,1:nw_recon) * (547.d0 * w(iLpppp^S,1:nw_recon) - 3882.d0 * w(iLppp^S,1:nw_recon) + 4642.d0 * w(iLpp^S,1:nw_recon) &
                          - 1854.d0 * w(iLp^S,1:nw_recon)) &
                          + w(iLppp^S,1:nw_recon) * (7043.d0 * w(iLppp^S,1:nw_recon) - 17246.d0 * w(iLpp^S,1:nw_recon) + 7042.d0 * w(iLp^S,1:nw_recon)) &
                          + w(iLpp^S,1:nw_recon) * (11003.d0 * w(iLpp^S,1:nw_recon) - 9402.d0 * w(iLp^S,1:nw_recon)) + 2107.d0 * w(iLp^S,1:nw_recon)**2
    beta(iL^S,1:nw_recon,2) = w(iLppp^S,1:nw_recon) * (267.d0 * w(iLppp^S,1:nw_recon) - 1642.d0 * w(iLpp^S,1:nw_recon) + 1602.d0 * w(iLp^S,1:nw_recon) &
                          - 494.d0 * w(iL^S,1:nw_recon))  &
                          + w(iLpp^S,1:nw_recon) * (2843.d0 * w(iLpp^S,1:nw_recon) - 5966.d0 * w(iLp^S,1:nw_recon) + 1922.d0 * w(iL^S,1:nw_recon)) &
                          + w(iLp^S,1:nw_recon) * (3443.d0 * w(iLp^S,1:nw_recon) - 2522.d0 * w(iL^S,1:nw_recon)) + 547.d0 * w(iL^S,1:nw_recon) ** 2
    beta(iL^S,1:nw_recon,3) = w(iLpp^S,1:nw_recon) * (547.d0 * w(iLpp^S,1:nw_recon) - 2522.d0 * w(iLp^S,1:nw_recon) + 1922.d0 * w(iL^S,1:nw_recon) &
                          - 494.d0 * w(iLm^S,1:nw_recon))  &
                          + w(iLp^S,1:nw_recon) * (3443.d0 * w(iLp^S,1:nw_recon) - 5966.d0 * w(iL^S,1:nw_recon) + 1602.d0 * w(iLm^S,1:nw_recon)) &
                          + w(iL^S,1:nw_recon) * (2843.d0 * w(iL^S,1:nw_recon) - 1642.d0 * w(iLm^S,1:nw_recon)) + 267.d0 * w(iLm^S,1:nw_recon) ** 2
    beta(iL^S,1:nw_recon,4) = w(iLp^S,1:nw_recon) * (2107.d0 * w(iLp^S,1:nw_recon) - 9402.d0 * w(iL^S,1:nw_recon) + 7042.d0 * w(iLm^S,1:nw_recon) &
                          - 1854.d0 * w(iLmm^S,1:nw_recon))  &
                          + w(iL^S,1:nw_recon) * (11003.d0 * w(iL^S,1:nw_recon) - 17246.d0 * w(iLm^S,1:nw_recon) + 4642.d0 * w(iLmm^S,1:nw_recon)) &
                          + w(iLm^S,1:nw_recon) * (7043.d0 * w(iLm^S,1:nw_recon) - 3882.d0 * w(iLmm^S,1:nw_recon)) + 547.d0 * w(iLmm^S,1:nw_recon) ** 2

    alpha_sum(iL^S,1:nw_recon) = 0.0d0 
    do i = 1,4
       alpha_array(iL^S,1:nw_recon,i) = d_array(i)/(beta(iL^S,1:nw_recon,i) + weno_eps_machine)**2
       alpha_sum(iL^S,1:nw_recon) = alpha_sum(iL^S,1:nw_recon) + alpha_array(iL^S,1:nw_recon,i)
    end do
    flux(iL^S,1:nw_recon) = 0.0d0
    do i = 1,4
       flux(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon) + f_array(iL^S,1:nw_recon,i) * alpha_array(iL^S,1:nw_recon,i)/(alpha_sum(iL^S,1:nw_recon))
    end do

    select case(var)
    case(1)
      wRCtmp(iL^S,1:nw_recon) = flux(iL^S,1:nw_recon)
    case(2)
      idmax^D=iLmax^D+kr(idims,^D); idmin^D=iLmin^D;
      idm^L=id^L-kr(idims,^D);
      idp^L=id^L+kr(idims,^D);
  
      iemax^D=idmax^D; iemin^D=idmin^D-kr(idims,^D);
      iem^L=ie^L-kr(idims,^D);
      iep^L=ie^L+kr(idims,^D);
      iepp^L=iep^L+kr(idims,^D);
  
      d(ie^S,1:nw_recon) = w(ie^S,1:nw_recon)-2.0d0*w(iep^S,1:nw_recon)+w(iepp^S,1:nw_recon)
  
      do iw = 1,nw_recon
        a(id^S) = 4.0d0*d(id^S,iw)-d(idm^S,iw)
        b(id^S) = 4.0d0*d(idm^S,iw)-d(id^S,iw)
        call minmod(ixI^L,id^L,a,b,tmp)
        a(id^S) = d(id^S,iw)
        b(id^S) = d(idm^S,iw)
        call minmod(ixI^L,id^L,a,b,tmp2)
        call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
        dm4(id^S,iw) = tmp3(id^S)
      end do
   
      flux_ul(iL^S,1:nw_recon) = w(iLp^S,1:nw_recon) + mpalpha * (w(iLp^S,1:nw_recon) - w(iLpp^S,1:nw_recon))
      flux_md(iL^S,1:nw_recon) = half * (w(iL^S,1:nw_recon) + w(iLp^S,1:nw_recon) - dm4(iL^S,1:nw_recon))
      flux_lc(iL^S,1:nw_recon) = half * (3.d0 * w(iLp^S,1:nw_recon) - w(iLpp^S,1:nw_recon)) + mpbeta / 3.d0 * dm4(iLp^S,1:nw_recon)
    
      flux_min(iL^S,1:nw_recon) = max(min(w(iLp^S,1:nw_recon), w(iL^S,1:nw_recon), flux_md(iL^S,1:nw_recon)), &
                                min(w(iLp^S,1:nw_recon), flux_ul(iL^S,1:nw_recon),flux_lc(iL^S,1:nw_recon)))
  
      flux_max(iL^S,1:nw_recon) = min(max(w(iLp^S,1:nw_recon), w(iL^S,1:nw_recon), flux_md(iL^S,1:nw_recon)), &
                                max(w(iLp^S,1:nw_recon), flux_ul(iL^S,1:nw_recon),flux_lc(iL^S,1:nw_recon)))
      do iw=1,nw_recon
        a(iL^S) = flux(iL^S,iw)
        b(iL^S) = flux_min(iL^S,iw)
        c(iL^S) = flux_max(iL^S,iw)
        call median(ixI^L, iL^L, a, b, c, tmp) 
        wRCtmp(iL^S,iw) = tmp(iL^S)
      end do
    end select

    call fix_limiter(ixI^L,iL^L,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine WENO7limiter

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
