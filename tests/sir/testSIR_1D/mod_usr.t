module mod_usr
  use mod_sir

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => sir_init
    usr_create_particles => place_samplingpoints

    call set_coordinate_system("Cartesian")
    call sir_activate()

  end subroutine usr_init

  ! initialize one grid
  subroutine sir_init(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision                :: x1, urand(ix^S)
    double precision                :: l1, dist2(ix^S)
    double precision                :: ampl1, ampl2
    logical                         :: mymask(ix^S)

    ampl1=1.0d-4
    ampl2=0.5d-4

    x1 = xprobmin1 + 0.5d0 * (xprobmax1 - xprobmin1)
    l1 = xprobmax1 - xprobmin1

    w(ix^S,re_) = 0.0d0
    w(ix^S,su_) = 1.0d0
    w(ix^S,in_) = 0.0d0

    select case (iprob)
    case (1)
        ! Default: steady state
        w(ix^S,su_) = 1.0d0
        w(ix^S,in_) = 0.0d0
    case (2)
        !1D test from Lotfi, figure 2 or figure 3
        ! on unit interval/square
        where(x(ix^S,1)<0.5d0)
           w(ix^S,su_)=1.1d0*x(ix^S,1)
           w(ix^S,in_)=0.5d0*x(ix^S,1)
        elsewhere
           w(ix^S,su_)=1.1d0*(1.0d0-x(ix^S,1))
           w(ix^S,in_)=0.5d0*(1.0d0-x(ix^S,1))
        endwhere
    case (3)
       ! Center square
       where (abs(x(ix^S, 1) - x1) < 0.1d0 * l1)
          w(ix^S,su_) = w(ix^S,su_)-ampl1
          w(ix^S,in_) = ampl1
       endwhere
    case (4)
       ! urand between [0,1]
       call random_number(urand)
       ! Center square with random noise
       where (abs(x(ix^S, 1) - x1) < 0.1d0 * l1)
          w(ix^S,su_) = w(ix^S,su_)-ampl1 * urand
          w(ix^S,in_) = ampl1* urand
       endwhere
    case (5)
       ! Two Gaussians
       dist2 = (x(ix^S, 1) - x1)**2 
       w(ix^S,su_) = w(ix^S,su_) - ampl1 * exp(-100 * dist2/l1**2)
       w(ix^S,in_) = w(ix^S,in_) + ampl1 * exp(-100 * dist2/l1**2)

       x1 = xprobmin1 + 0.7d0 * l1
       dist2 = (x(ix^S, 1) - x1)**2 
       w(ix^S,su_) = w(ix^S,su_) - ampl2 * exp(-100 * dist2/l1**2)
       w(ix^S,in_) = w(ix^S,in_) + ampl2 * exp(-100 * dist2/l1**2)

    case (6)
       ! Junction of two circular shapes
       dist2 = (x(ix^S, 1) - x1)**2
       mymask = (dist2 < (0.05d0*l1)**2)

       x1 = xprobmin1 + 0.55d0 * l1
       dist2 = (x(ix^S, 1) - x1)**2 
       mymask = mymask .or. (dist2 < (0.1d0*l1)**2)

       where (mymask)
          w(ix^S,su_) = w(ix^S, su_)-ampl1
          w(ix^S,in_) = w(ix^S, in_)+ampl1
       end where
    case default
       call mpistop("Unknown iprob")
    end select

  end subroutine sir_init

  subroutine place_samplingpoints(n_particles, x, v, q, m, follow)
     integer, intent(in)           :: n_particles
     double precision, intent(out) :: x(3, n_particles)
     double precision, intent(out) :: v(3, n_particles)
     double precision, intent(out) :: q(n_particles)
     double precision, intent(out) :: m(n_particles)
     logical, intent(out)          :: follow(n_particles)
     
     integer :: nsize,i,j,i_part
     double precision :: xmid,dxx

     v = 0.0d0
     q = 0.0d0
     m = 0.0d0
     nsize=n_particles
     xmid=(xprobmax1-xprobmin1)*0.5d0+xprobmin1
     dxx=(xprobmax1-xprobmin1)*0.5d0/nsize
     xmid=xmid-(nsize/2)*dxx
     i_part=0
     do i=1,nsize
        i_part=i_part+1
        x(1,i_part)=xmid
        ! print *,i_part,x(1,i_part)
        follow(i_part)=.true.
        xmid=xmid+dxx
     enddo
     if (i_part/=n_particles) then
         print *,nsize,n_particles,i_part
         call mpistop("error in place_samplepoints")
     endif
           
  end subroutine place_samplingpoints 

end module mod_usr
