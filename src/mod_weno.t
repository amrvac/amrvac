module mod_weno
  ! WENO-JS3 & WENO-JS5 & WENO-Z5 & WENO-ZP5 limiter
  !
  ! 2019.9.19 WENO-JS5 is transplant from BHAC code by nanami;
  ! 2019.9.20 WENO-JS3 is coded up by nanami;
  ! 2019.9.21 WENO-Z5 is coded up by nanami;
  ! 2019.9.22 WENO-Z+5 is transplant from BHAC code by nanami;
  !
  ! see Liu et al. 1994 for the original idea of WENO;
  ! see Jiang & Shu 1996 for the basic idea of WENO-JS;
  ! see Shu 2009 (SIAM review) for the implement of WENO-JS5,
  ! especially for the difference between Reconstruction and Interpolation.
  ! so that WENO-JS3 can be deduced in a similar way.
  ! see Borges et al. 2008 for WENO-Z variation;
  ! see Acker et al. 2016 for WENO-Z+ variation.
   
  implicit none
  private

  public :: WENOJS3limiter
  public :: WENOJS5limiter
  public :: WENOZ5limiter
  public :: WENOZP5limiter

contains

  subroutine WENOJS5limiter(ixI^L,iL^L,idims,w,wLC,wRC,reconstruction)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims
    logical, intent(in)             :: reconstruction
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    !> local
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    double precision                :: f_array(ixI^S,1:nw,3), d_array(3)
    double precision                :: beta(ixI^S,1:nw,3), beta_coeff(2)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:nw,3), alpha_sum(ixI^S,1:nw), flux(ixI^S,1:nw)
    integer                         :: i, iw
    double precision, parameter     :: weno_eps_machine = 1.0d-12

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  
    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if(reconstruction) then
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
      u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
      u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    else
      d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
      u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
      u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
      u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    endif
    
    !> left side
    !> eq(2.11)
    f_array(iL^S,1:nwflux,1) = u1_coeff(1) * w(iLmm^S,1:nwflux) + u1_coeff(2) * w(iLm^S,1:nwflux) + u1_coeff(3) * w(iL^S,1:nwflux)
    !> eq(2.12)
    f_array(iL^S,1:nwflux,2) = u2_coeff(1) * w(iLm^S,1:nwflux)  + u2_coeff(2) * w(iL^S,1:nwflux)  + u2_coeff(3) * w(iLp^S,1:nwflux)
    !> eq(2.13)
    f_array(iL^S,1:nwflux,3) = u3_coeff(1) * w(iL^S,1:nwflux)   + u3_coeff(2) * w(iLp^S,1:nwflux) + u3_coeff(3) * w(iLpp^S,1:nwflux)  
  
    !> eq(2.17)
    beta(iL^S,1:nwflux,1) = beta_coeff(1) * (w(iLmm^S,1:nwflux) + w(iL^S,1:nwflux) - 2.0d0*w(iLm^S,1:nwflux))**2 &
         + beta_coeff(2) * (w(iLmm^S,1:nwflux) - 4.0d0 * w(iLm^S,1:nwflux) + 3.0d0*w(iL^S,1:nwflux))**2
    beta(iL^S,1:nwflux,2) = beta_coeff(1) * (w(iLm^S,1:nwflux) + w(iLp^S,1:nwflux) - 2.0d0 * w(iL^S,1:nwflux))**2 &
         + beta_coeff(2) * (w(iLm^S,1:nwflux) - w(iLp^S,1:nwflux))**2
    beta(iL^S,1:nwflux,3) = beta_coeff(1) * (w(iL^S,1:nwflux) + w(iLpp^S,1:nwflux) - 2.0d0 * w(iLp^S,1:nwflux))**2 &
         + beta_coeff(2) * (3.0d0 * w(iL^S, 1:nwflux) - 4.0d0 * w(iLp^S,1:nwflux) + w(iLpp^S,1:nwflux))**2
  
    !> eq(2.10)
    alpha_sum(iL^S,1:nwflux) = 0.0d0 
    do i = 1,3
       alpha_array(iL^S,1:nwflux,i) = d_array(i)/(beta(iL^S,1:nwflux,i) + weno_eps_machine)**2
       alpha_sum(iL^S,1:nwflux) = alpha_sum(iL^S,1:nwflux) + alpha_array(iL^S,1:nwflux,i)
    end do
    flux(iL^S,1:nwflux) = 0.0d0
    do i = 1,3
       flux(iL^S,1:nwflux) = flux(iL^S,1:nwflux) + f_array(iL^S,1:nwflux,i) * alpha_array(iL^S,1:nwflux,i)/(alpha_sum(iL^S,1:nwflux))
    end do
  
    !> left value at right interface
    wLC(iL^S,1:nwflux) = flux(iL^S,1:nwflux)
  
    !> right side
    f_array(iL^S,1:nwflux,1) = u1_coeff(1) * w(iLppp^S,1:nwflux) + u1_coeff(2) * w(iLpp^S,1:nwflux) + u1_coeff(3) * w(iLp^S,1:nwflux)
    f_array(iL^S,1:nwflux,2) = u2_coeff(1) * w(iLpp^S,1:nwflux)  + u2_coeff(2) * w(iLp^S,1:nwflux)  + u2_coeff(3) * w(iL^S,1:nwflux)
    f_array(iL^S,1:nwflux,3) = u3_coeff(1) * w(iLp^S,1:nwflux)   + u3_coeff(2) * w(iL^S,1:nwflux)   + u3_coeff(3) * w(iLm^S,1:nwflux)  
  
    beta(iL^S,1:nwflux,1) = beta_coeff(1) * (w(iLppp^S,1:nwflux) + w(iLp^S,1:nwflux) - 2.0d0*w(iLpp^S,1:nwflux))**2 &
         + beta_coeff(2) * (w(iLppp^S,1:nwflux) - 4.0d0 * w(iLpp^S,1:nwflux) + 3.0d0*w(iLp^S,1:nwflux))**2
    beta(iL^S,1:nwflux,2) = beta_coeff(1) * (w(iLpp^S,1:nwflux) + w(iL^S,1:nwflux) - 2.0d0 * w(iLp^S,1:nwflux))**2 &
         + beta_coeff(2) * (w(iLpp^S,1:nwflux) - w(iL^S,1:nwflux))**2
    beta(iL^S,1:nwflux,3) = beta_coeff(1) * (w(iLp^S,1:nwflux) + w(iLm^S,1:nwflux) - 2.0d0 * w(iL^S,1:nwflux))**2 &
         + beta_coeff(2) * (3.0d0 * w(iLp^S, 1:nwflux) - 4.0d0 * w(iL^S,1:nwflux) + w(iLm^S,1:nwflux))**2
  
    alpha_sum(iL^S,1:nwflux) = 0.0d0 
    do i = 1,3
       alpha_array(iL^S,1:nwflux,i) = d_array(i)/(beta(iL^S,1:nwflux,i) + weno_eps_machine)**2
       alpha_sum(iL^S,1:nwflux) = alpha_sum(iL^S,1:nwflux) + alpha_array(iL^S,1:nwflux,i)
    end do
    flux(iL^S,1:nwflux) = 0.0d0
    do i = 1,3
       flux(iL^S,1:nwflux) = flux(iL^S,1:nwflux) + f_array(iL^S,1:nwflux,i) * alpha_array(iL^S,1:nwflux,i)/(alpha_sum(iL^S,1:nwflux))
    end do
  
    !> right value at right interface
    wRC(iL^S,1:nwflux) = flux(iL^S,1:nwflux)

  end subroutine WENOJS5limiter

  subroutine WENOJS3limiter(ixI^L,iL^L,idims,w,wLC,wRC)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    !> local
    integer                         :: iLm^L, iLp^L, iLpp^L
    double precision                :: f_array(ixI^S,1:nw,2), d_array(2)
    double precision                :: beta(ixI^S,1:nw,2)
    double precision                :: u1_coeff(2), u2_coeff(2)
    double precision                :: alpha_array(ixI^S,1:nw,2), alpha_sum(ixI^S,1:nw), flux(ixI^S,1:nw)
    integer                         :: i, iw
    double precision, parameter     :: weno_eps_machine = 1.0d-12

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  
    iLm^L=iL^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    d_array(1:2) = (/ 1.0d0/4.0d0, 3.0d0/4.0d0 /)
    u1_coeff(1:2) = (/ -1.d0/2.d0, 3.d0/2.d0 /)
    u2_coeff(1:2) = (/ 1.d0/2.d0, 1.d0/2.d0 /)
    
    !> left side
    f_array(iL^S,1:nwflux,1) = u1_coeff(1) * w(iLm^S,1:nwflux) + u1_coeff(2) * w(iL^S,1:nwflux)
    f_array(iL^S,1:nwflux,2) = u2_coeff(1) * w(iL^S,1:nwflux)  + u2_coeff(2) * w(iLp^S,1:nwflux)
  
    !> eq(2.17)
    beta(iL^S,1:nwflux,1) = (w(iL^S,1:nwflux) - w(iLm^S,1:nwflux))**2
    beta(iL^S,1:nwflux,2) = (w(iLp^S,1:nwflux) - w(iL^S,1:nwflux))**2
  
    !> eq(2.10)
    alpha_sum(iL^S,1:nwflux) = 0.0d0 
    do i = 1,2
       alpha_array(iL^S,1:nwflux,i) = d_array(i)/(beta(iL^S,1:nwflux,i) + weno_eps_machine)**2
       alpha_sum(iL^S,1:nwflux) = alpha_sum(iL^S,1:nwflux) + alpha_array(iL^S,1:nwflux,i)
    end do
    flux(iL^S,1:nwflux) = 0.0d0
    do i = 1,2
       flux(iL^S,1:nwflux) = flux(iL^S,1:nwflux) + f_array(iL^S,1:nwflux,i) * alpha_array(iL^S,1:nwflux,i)/(alpha_sum(iL^S,1:nwflux))
    end do
  
    !> left value at right interface
    wLC(iL^S,1:nwflux) = flux(iL^S,1:nwflux)
  
    !> right side
    f_array(iL^S,1:nwflux,1) = u1_coeff(1) * w(iLpp^S,1:nwflux) + u1_coeff(2) * w(iLp^S,1:nwflux)
    f_array(iL^S,1:nwflux,2) = u2_coeff(1) * w(iLp^S,1:nwflux)  + u2_coeff(2) * w(iL^S,1:nwflux)
  
    beta(iL^S,1:nwflux,1) = (w(iLpp^S,1:nwflux) - w(iLp^S,1:nwflux))**2
    beta(iL^S,1:nwflux,2) = (w(iLp^S,1:nwflux) - w(iL^S,1:nwflux))**2
  
    alpha_sum(iL^S,1:nwflux) = 0.0d0 
    do i = 1,2
       alpha_array(iL^S,1:nwflux,i) = d_array(i)/(beta(iL^S,1:nwflux,i) + weno_eps_machine)**2
       alpha_sum(iL^S,1:nwflux) = alpha_sum(iL^S,1:nwflux) + alpha_array(iL^S,1:nwflux,i)
    end do
    flux(iL^S,1:nwflux) = 0.0d0
    do i = 1,2
       flux(iL^S,1:nwflux) = flux(iL^S,1:nwflux) + f_array(iL^S,1:nwflux,i) * alpha_array(iL^S,1:nwflux,i)/(alpha_sum(iL^S,1:nwflux))
    end do
  
    !> right value at right interface
    wRC(iL^S,1:nwflux) = flux(iL^S,1:nwflux)

  end subroutine WENOJS3limiter

  subroutine WENOZ5limiter(ixI^L,iL^L,idims,w,wLC,wRC,reconstruction)
    use mod_global_parameters
  
    integer, intent(in)             :: ixI^L, iL^L, idims
    logical, intent(in)             :: reconstruction
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    !> local
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    double precision                :: f_array(ixI^S,1:nw,3), d_array(3)
    double precision                :: beta(ixI^S,1:nw,3), beta_coeff(2)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:nw,3), alpha_sum(ixI^S,1:nw), flux(ixI^S,1:nw)
    double precision                :: tau(ixI^S,1:nw)
    integer                         :: i, iw
    double precision, parameter     :: weno_eps_machine = 1.0d-12

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  
    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if(reconstruction) then
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
      u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
      u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    else
      d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
      u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
      u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
      u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    endif
    
    !> left side
    !> same with WENO-JS5
    f_array(iL^S,1:nwflux,1) = u1_coeff(1) * w(iLmm^S,1:nwflux) + u1_coeff(2) * w(iLm^S,1:nwflux) + u1_coeff(3) * w(iL^S,1:nwflux)
    f_array(iL^S,1:nwflux,2) = u2_coeff(1) * w(iLm^S,1:nwflux)  + u2_coeff(2) * w(iL^S,1:nwflux)  + u2_coeff(3) * w(iLp^S,1:nwflux)
    f_array(iL^S,1:nwflux,3) = u3_coeff(1) * w(iL^S,1:nwflux)   + u3_coeff(2) * w(iLp^S,1:nwflux) + u3_coeff(3) * w(iLpp^S,1:nwflux)  
  
    !> same with WENO-JS5
    beta(iL^S,1:nwflux,1) = beta_coeff(1) * (w(iLmm^S,1:nwflux) + w(iL^S,1:nwflux) - 2.0d0*w(iLm^S,1:nwflux))**2 &
         + beta_coeff(2) * (w(iLmm^S,1:nwflux) - 4.0d0 * w(iLm^S,1:nwflux) + 3.0d0*w(iL^S,1:nwflux))**2
    beta(iL^S,1:nwflux,2) = beta_coeff(1) * (w(iLm^S,1:nwflux) + w(iLp^S,1:nwflux) - 2.0d0 * w(iL^S,1:nwflux))**2 &
         + beta_coeff(2) * (w(iLm^S,1:nwflux) - w(iLp^S,1:nwflux))**2
    beta(iL^S,1:nwflux,3) = beta_coeff(1) * (w(iL^S,1:nwflux) + w(iLpp^S,1:nwflux) - 2.0d0 * w(iLp^S,1:nwflux))**2 &
         + beta_coeff(2) * (3.0d0 * w(iL^S, 1:nwflux) - 4.0d0 * w(iLp^S,1:nwflux) + w(iLpp^S,1:nwflux))**2

    !> eq(25)
    tau(iL^S,1:nwflux) = abs(beta(iL^S,1:nwflux,1) - beta(iL^S,1:nwflux,3))
  
    !> eq(28)
    !> Borges 2008 uses L1 norm, here we use L2 norm as suggested by Acker 2016
    alpha_sum(iL^S,1:nwflux) = 0.0d0 
    do i = 1,3
       alpha_array(iL^S,1:nwflux,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nwflux) / &
                                      (beta(iL^S,1:nwflux,i) + weno_eps_machine))**2)
       alpha_sum(iL^S,1:nwflux) = alpha_sum(iL^S,1:nwflux) + alpha_array(iL^S,1:nwflux,i)
    end do
    flux(iL^S,1:nwflux) = 0.0d0
    do i = 1,3
       flux(iL^S,1:nwflux) = flux(iL^S,1:nwflux) + f_array(iL^S,1:nwflux,i) * alpha_array(iL^S,1:nwflux,i)/(alpha_sum(iL^S,1:nwflux))
    end do
  
    !> left value at right interface
    wLC(iL^S,1:nwflux) = flux(iL^S,1:nwflux)
  
    !> right side
    f_array(iL^S,1:nwflux,1) = u1_coeff(1) * w(iLppp^S,1:nwflux) + u1_coeff(2) * w(iLpp^S,1:nwflux) + u1_coeff(3) * w(iLp^S,1:nwflux)
    f_array(iL^S,1:nwflux,2) = u2_coeff(1) * w(iLpp^S,1:nwflux)  + u2_coeff(2) * w(iLp^S,1:nwflux)  + u2_coeff(3) * w(iL^S,1:nwflux)
    f_array(iL^S,1:nwflux,3) = u3_coeff(1) * w(iLp^S,1:nwflux)   + u3_coeff(2) * w(iL^S,1:nwflux)   + u3_coeff(3) * w(iLm^S,1:nwflux)  
  
    beta(iL^S,1:nwflux,1) = beta_coeff(1) * (w(iLppp^S,1:nwflux) + w(iLp^S,1:nwflux) - 2.0d0*w(iLpp^S,1:nwflux))**2 &
         + beta_coeff(2) * (w(iLppp^S,1:nwflux) - 4.0d0 * w(iLpp^S,1:nwflux) + 3.0d0*w(iLp^S,1:nwflux))**2
    beta(iL^S,1:nwflux,2) = beta_coeff(1) * (w(iLpp^S,1:nwflux) + w(iL^S,1:nwflux) - 2.0d0 * w(iLp^S,1:nwflux))**2 &
         + beta_coeff(2) * (w(iLpp^S,1:nwflux) - w(iL^S,1:nwflux))**2
    beta(iL^S,1:nwflux,3) = beta_coeff(1) * (w(iLp^S,1:nwflux) + w(iLm^S,1:nwflux) - 2.0d0 * w(iL^S,1:nwflux))**2 &
         + beta_coeff(2) * (3.0d0 * w(iLp^S, 1:nwflux) - 4.0d0 * w(iL^S,1:nwflux) + w(iLm^S,1:nwflux))**2
    
    tau(iL^S,1:nwflux) = abs(beta(iL^S,1:nwflux,1) - beta(iL^S,1:nwflux,3))
  
    alpha_sum(iL^S,1:nwflux) = 0.0d0 
    do i = 1,3
       alpha_array(iL^S,1:nwflux,i) = d_array(i) * (1.d0 + (tau(iL^S,1:nwflux) / &
                                      (beta(iL^S,1:nwflux,i) + weno_eps_machine))**2)
       alpha_sum(iL^S,1:nwflux) = alpha_sum(iL^S,1:nwflux) + alpha_array(iL^S,1:nwflux,i)
    end do
    flux(iL^S,1:nwflux) = 0.0d0
    do i = 1,3
       flux(iL^S,1:nwflux) = flux(iL^S,1:nwflux) + f_array(iL^S,1:nwflux,i) * alpha_array(iL^S,1:nwflux,i)/(alpha_sum(iL^S,1:nwflux))
    end do
  
    !> right value at right interface
    wRC(iL^S,1:nwflux) = flux(iL^S,1:nwflux)

  end subroutine WENOZ5limiter

  subroutine WENOZP5limiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC,reconstruction)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, iL^L, idims
    logical, intent(in)             :: reconstruction
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    !> local
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    double precision                :: f_array(ixI^S,1:nw,3), d_array(3)
    double precision                :: beta(ixI^S,1:nw,3), beta_coeff(2)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixI^S,1:nw,3), alpha_sum(ixI^S,1:nw), flux(ixI^S,1:nw)
    double precision                :: tau(ixI^S,1:nw), tmp(ixI^S,1:nw)
    double precision                :: lambda
    integer                         :: i, iw
    double precision, parameter     :: weno_eps_machine = 1.0d-12
    double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  
    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);

    lambda = dxdim**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    if(reconstruction) then
      d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
      u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
      u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
      u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    else
      d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
      u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
      u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
      u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    endif

    !> left side
    !> same with WENO-JS5
    f_array(iL^S,1:nwflux,1) = u1_coeff(1) * w(iLmm^S,1:nwflux) + u1_coeff(2) * w(iLm^S,1:nwflux) + u1_coeff(3) * w(iL^S,1:nwflux)
    f_array(iL^S,1:nwflux,2) = u2_coeff(1) * w(iLm^S,1:nwflux)  + u2_coeff(2) * w(iL^S,1:nwflux)  + u2_coeff(3) * w(iLp^S,1:nwflux)
    f_array(iL^S,1:nwflux,3) = u3_coeff(1) * w(iL^S,1:nwflux)   + u3_coeff(2) * w(iLp^S,1:nwflux) + u3_coeff(3) * w(iLpp^S,1:nwflux)  
  
    !> same with WENO-JS5
    beta(iL^S,1:nwflux,1) = beta_coeff(1) * (w(iLmm^S,1:nwflux) + w(iL^S,1:nwflux) - 2.0d0*w(iLm^S,1:nwflux))**2 &
         + beta_coeff(2) * (w(iLmm^S,1:nwflux) - 4.0d0 * w(iLm^S,1:nwflux) + 3.0d0*w(iL^S,1:nwflux))**2
    beta(iL^S,1:nwflux,2) = beta_coeff(1) * (w(iLm^S,1:nwflux) + w(iLp^S,1:nwflux) - 2.0d0 * w(iL^S,1:nwflux))**2 &
         + beta_coeff(2) * (w(iLm^S,1:nwflux) - w(iLp^S,1:nwflux))**2
    beta(iL^S,1:nwflux,3) = beta_coeff(1) * (w(iL^S,1:nwflux) + w(iLpp^S,1:nwflux) - 2.0d0 * w(iLp^S,1:nwflux))**2 &
         + beta_coeff(2) * (3.0d0 * w(iL^S, 1:nwflux) - 4.0d0 * w(iLp^S,1:nwflux) + w(iLpp^S,1:nwflux))**2

    !> same with WENO-Z
    tau(iL^S,1:nwflux) = abs(beta(iL^S,1:nwflux,1) - beta(iL^S,1:nwflux,3))

    !> eq(10)  
    alpha_sum(iL^S,1:nwflux) = 0.0d0 
    do i = 1,3
      tmp(iL^S,1:nwflux) = (tau(iL^S,1:nwflux) + weno_eps_machine) / (beta(iL^S,1:nwflux,i) + weno_eps_machine)
      alpha_array(iL^S,1:nwflux,i) = d_array(i) * (1.0d0 + tmp(iL^S,1:nwflux)**2 + lambda/tmp(iL^S,1:nwflux))
      alpha_sum(iL^S,1:nwflux) = alpha_sum(iL^S,1:nwflux) + alpha_array(iL^S,1:nwflux,i)
    end do

    flux(iL^S,1:nwflux) = 0.0d0
    do i = 1,3
      flux(iL^S,1:nwflux) = flux(iL^S,1:nwflux) + f_array(iL^S,1:nwflux,i) * alpha_array(iL^S,1:nwflux,i)/(alpha_sum(iL^S,1:nwflux))
    end do
  
    !> left value at right interface
    wLC(iL^S,1:nwflux) = flux(iL^S,1:nwflux)

    !> right side
    f_array(iL^S,1:nwflux,1) = u1_coeff(1) * w(iLppp^S,1:nwflux) + u1_coeff(2) * w(iLpp^S,1:nwflux) + u1_coeff(3) * w(iLp^S,1:nwflux)
    f_array(iL^S,1:nwflux,2) = u2_coeff(1) * w(iLpp^S,1:nwflux)  + u2_coeff(2) * w(iLp^S,1:nwflux)  + u2_coeff(3) * w(iL^S,1:nwflux)
    f_array(iL^S,1:nwflux,3) = u3_coeff(1) * w(iLp^S,1:nwflux)   + u3_coeff(2) * w(iL^S,1:nwflux)   + u3_coeff(3) * w(iLm^S,1:nwflux)  
  
    beta(iL^S,1:nwflux,1) = beta_coeff(1) * (w(iLppp^S,1:nwflux) + w(iLp^S,1:nwflux) - 2.0d0*w(iLpp^S,1:nwflux))**2 &
        + beta_coeff(2) * (w(iLppp^S,1:nwflux) - 4.0d0 * w(iLpp^S,1:nwflux) + 3.0d0*w(iLp^S,1:nwflux))**2
    beta(iL^S,1:nwflux,2) = beta_coeff(1) * (w(iLpp^S,1:nwflux) + w(iL^S,1:nwflux) - 2.0d0 * w(iLp^S,1:nwflux))**2 &
        + beta_coeff(2) * (w(iLpp^S,1:nwflux) - w(iL^S,1:nwflux))**2
    beta(iL^S,1:nwflux,3) = beta_coeff(1) * (w(iLp^S,1:nwflux) + w(iLm^S,1:nwflux) - 2.0d0 * w(iL^S,1:nwflux))**2 &
        + beta_coeff(2) * (3.0d0 * w(iLp^S, 1:nwflux) - 4.0d0 * w(iL^S,1:nwflux) + w(iLm^S,1:nwflux))**2
    
    tau(iL^S,1:nwflux) = abs(beta(iL^S,1:nwflux,1) - beta(iL^S,1:nwflux,3))
  
    alpha_sum(iL^S,1:nwflux) = 0.0d0 

    do i = 1,3
      tmp(iL^S,1:nwflux) = (tau(iL^S,1:nwflux) + weno_eps_machine) / (beta(iL^S,1:nwflux,i) + weno_eps_machine)
      alpha_array(iL^S,1:nwflux,i) = d_array(i) * (1.0d0 + tmp(iL^S,1:nwflux)**2 + lambda/tmp(iL^S,1:nwflux))
      alpha_sum(iL^S,1:nwflux) = alpha_sum(iL^S,1:nwflux) + alpha_array(iL^S,1:nwflux,i)
    end do

    flux(iL^S,1:nwflux) = 0.0d0
    do i = 1,3
      flux(iL^S,1:nwflux) = flux(iL^S,1:nwflux) + f_array(iL^S,1:nwflux,i) * alpha_array(iL^S,1:nwflux,i)/(alpha_sum(iL^S,1:nwflux))
    end do

    !> right value at right interface
    wRC(iL^S,1:nwflux) = flux(iL^S,1:nwflux)

  end subroutine WENOZP5limiter

end module mod_weno
