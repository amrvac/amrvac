module mod_riemann_exact
  implicit none
  private

  integer, parameter :: dp   = kind(0.0d0)

  type riemann_t
     real(dp) :: rhol
     real(dp) :: pl
     real(dp) :: ul
     real(dp) :: rhor
     real(dp) :: pr
     real(dp) :: ur
     real(dp) :: gamma
     real(dp) :: t
     real(dp) :: xi
     real(dp) :: xl
     real(dp) :: xr
  end type riemann_t

  public :: riemann_t
  public :: riemann_solve

contains

  ! solves the exact riemann problem
  ! original code from Bruce Fryxell
  ! Adapted from http://cococubed.asu.edu/code_pages/exact_riemann.shtml
  subroutine riemann_solve(rp, npts, x, rho, u, p)
    integer, intent(in)         :: npts
    type(riemann_t), intent(in) :: rp
    real(dp), intent(out)       :: x(npts)
    real(dp), intent(out)       :: rho(npts)
    real(dp), intent(out)       :: u(npts)
    real(dp), intent(out)       :: p(npts)

    integer  :: itmax,iter,i
    real(dp) :: rho1,p1,u1,rho5,p5,u5,p40,p41,f0,eps
    real(dp) :: f1,p4,error,z,c5,gm1,gp1,gmfac1,gmfac2,fact
    real(dp) :: u4,rho4,w,p3,u3,rho3,c1,c3,xsh,xcd,xft,xhd,dx

    if (rp%xr < rp%xl) then
       stop 'xr must be greater than xl'
    endif

    if (rp%ul /= 0. .or. rp%ur /= 0.) then
       stop 'must have ul = ur = 0.'
    endif

    ! begin solution
    if (rp%pl > rp%pr) then
       rho1 = rp%rhol
       p1   = rp%pl
       u1   = rp%ul
       rho5 = rp%rhor
       p5   = rp%pr
       u5   = rp%ur
    else
       rho1 = rp%rhor
       p1   = rp%pr
       u1   = rp%ur
       rho5 = rp%rhol
       p5   = rp%pl
       u5   = rp%ul
    endif

    ! solve for post-shock pressure by secant method
    ! initial guesses

    p40 = p1
    p41 = p5
    f0  = f(p40, p1, p5, rho1, rho5, rp%gamma)

    ! maximum number of iterations and maxium allowable relative error
    itmax = 20
    eps   = 1.e-5

    do iter = 1, itmax
       f1 = f(p41, p1, p5, rho1, rho5, rp%gamma)
       if (f1 .eq. f0) go to 10

       p4 = p41 - (p41 - p40) * f1 / (f1 - f0)

       error = abs (p4 - p41) / p41
       if (error < eps) go to 10

       p40 = p41
       p41 = p4
       f0  = f1
    enddo
    print *, 'iteration failed to converge'
    stop 'abnormal termination'

10  continue

    ! compute post-shock density and velocity
    z  = (p4 / p5 - 1.)
    c5 = sqrt (rp%gamma * p5 / rho5)

    gm1 = rp%gamma - 1.
    gp1 = rp%gamma + 1.
    gmfac1 = 0.5 * gm1 / rp%gamma
    gmfac2 = 0.5 * gp1 / rp%gamma

    fact = sqrt (1. + gmfac2 * z)

    u4 = c5 * z / (rp%gamma * fact)
    rho4 = rho5 * (1. + gmfac2 * z) / (1. + gmfac1 * z)

    ! shock speed
    w = c5 * fact

    ! compute values at foot of rarefaction
    p3 = p4
    u3 = u4
    rho3 = rho1 * (p3 / p1)**(1. /rp%gamma)

    ! compute positions of waves
    if (rp%pl > rp%pr) then
       c1 = sqrt (rp%gamma * p1 / rho1)
       c3 = sqrt (rp%gamma * p3 / rho3)

       xsh = rp%xi + w * rp%t
       xcd = rp%xi + u3 * rp%t
       xft = rp%xi + (u3 - c3) * rp%t
       xhd = rp%xi - c1 * rp%t

       ! and do say what we found
       ! write (6, 500)
       ! write (6, 501) rho1, p1, u1
       ! write (6, 502)
       ! write (6, 503) rho3, p3, u3
       ! write (6, 504) rho4, p4, u4
       ! write (6, 505) rho5, p5, u5

       ! write (6, 506) xhd
       ! write (6, 507) xft
       ! write (6, 508) xcd
       ! write (6, 509) xsh

       ! 500    format (// 2x, 'Region', 4x, 'Density', 8x, 'Pressure', &
       !             8 x, 'Velocity')
       ! 501    format (5x, '1', 3(2x,1pe14.7))
       ! 502    format (5x, '2' ,20x, 'RAREFACTION')
       ! 503    format (5x, '3', 3(2x,1pe14.7))
       ! 504    format (5x, '4', 3(2x,1pe14.7))
       ! 505    format (5x, '5', 3(2x,1pe14.7)//)

       ! 506    format (2x, 'Head Of Rarefaction    x = ', 1pe14.7)
       ! 507    format (2x, 'Foot Of Rarefaction    x = ', 1pe14.7)
       ! 508    format (2x, 'Contact Discontinuity  x = ', 1pe14.7)
       ! 509    format (2x, 'Shock                  x = ', 1pe14.7//)

       ! compute solution as a function of position
       dx = (rp%xr - rp%xl) / (npts - 1)
       do i = 1, npts
          x(i) = rp%xl + dx * (i - 1)
       enddo

       do i = 1, npts
          if (x(i) < xhd) then
             rho(i) = rho1
             p(i)   = p1
             u(i)   = u1
          else if (x(i) < xft) then
             u(i)   = 2. / gp1 * (c1 + (x(i) - rp%xi) / rp%t)
             fact   = 1. - 0.5 * gm1 * u(i) / c1
             rho(i) = rho1 * fact ** (2. / gm1)
             p(i)   = p1 * fact ** (2. * rp%gamma / gm1)
          else if (x(i) < xcd) then
             rho(i) = rho3
             p(i)   = p3
             u(i)   = u3
          else if (x(i) < xsh) then
             rho(i) = rho4
             p(i)   = p4
             u(i)   = u4
          else
             rho(i) = rho5
             p(i)   = p5
             u(i)   = u5
          endif
       enddo
    endif


    ! if pr > pl, reverse solution
    if (rp%pr > rp%pl) then
       c1 = sqrt (rp%gamma * p1 / rho1)
       c3 = sqrt (rp%gamma * p3 / rho3)

       xsh = rp%xi - w * rp%t
       xcd = rp%xi - u3 * rp%t
       xft = rp%xi - (u3 - c3) * rp%t
       xhd = rp%xi + c1 * rp%t


       ! and do say what we found
       !        write (6, 500)
       !        write (6, 501) rho5, p5, u5
       !        write (6, 602) rho4, p4, u4
       !        write (6, 503) rho3, p3, u3
       !        write (6, 604)
       !        write (6, 505) rho1, p1, u1

       !        write (6, 609) xsh
       !        write (6, 508) xcd
       !        write (6, 507) xft
       !        write (6, 606) xhd

       ! 602    format (5x, '2', 3(2x,1pe14.7))
       ! 604    format (5x, '4' ,20x, 'RAREFACTION')
       ! 606    format (2x, 'Head Of Rarefaction    x = ', 1pe14.7//)
       ! 609    format (2x, 'Shock                  x = ', 1pe14.7)

       dx = (rp%xr - rp%xl) / (npts - 1)
       do i = 1, npts
          x(i) = rp%xl + dx * (i - 1)
       enddo

       do i = 1, npts
          if (x(i) < xsh) then
             rho(i) = rho5
             p(i)   = p5
             u(i)   = -u5
          else if (x(i) < xcd) then
             rho(i) = rho4
             p(i)   = p4
             u(i)   = -u4
          else if (x(i) < xft) then
             rho(i) = rho3
             p(i)   = p3
             u(i)   = -u3
          else if (x(i) < xhd) then
             u(i)   = -2. / gp1 * (c1 + (rp%xi - x(i)) / rp%t)
             fact   = 1. + 0.5 * gm1 * u(i) / c1
             rho(i) = rho1 * fact ** (2. / gm1)
             p(i)   = p1 * fact ** (2. * rp%gamma / gm1)
          else
             rho(i) = rho1
             p(i)   = p1
             u(i)   = -u1
          endif
       enddo
    endif

  end subroutine riemann_solve

  real(dp) function f(p4, p1, p5, rho1, rho5, gamma)
    implicit none

    ! shock tube equation

    ! declare the pass
    real(dp) :: p4,p1,p5,rho1,rho5,gamma

    ! local variables
    real(dp) :: z,c1,c5,gm1,gp1,g2,fact

    z = (p4 / p5 - 1.)
    c1 = sqrt (gamma * p1 / rho1)
    c5 = sqrt (gamma * p5 / rho5)

    gm1 = gamma - 1.
    gp1 = gamma + 1.
    g2  = 2. * gamma

    fact = gm1 / g2 * (c5 / c1) * z / sqrt (1. + gp1 / g2 * z)
    fact = (1. - fact) ** (g2 / gm1)

    f = p1 * fact - p4

    return
  end function f

end module mod_riemann_exact
