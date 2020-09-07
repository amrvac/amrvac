module mod_usr
  use mod_rd

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => gray_scott_init
    usr_create_particles => place_samplingpoints

    call rd_activate()

  end subroutine usr_init

  ! initialize one grid
  subroutine gray_scott_init(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision                :: x1, x2, urand(ix^S)
    double precision                :: l1, l2, dist2(ix^S)
    logical                         :: mymask(ix^S)

    x1 = xprobmin1 + 0.5d0 * (xprobmax1 - xprobmin1)
    x2 = xprobmin2 + 0.5d0 * (xprobmax2 - xprobmin2)
    l1 = xprobmax1 - xprobmin1
    l2 = xprobmax2 - xprobmin2

    ! Default: steady state
    w(ix^S,u_) = 1.0d0
    w(ix^S,v_) = 0.0d0

    call random_number(urand)

    select case (iprob)
    case (1)
       ! Center square
       where (abs(x(ix^S, 1) - x1) < 0.1d0 * l1 .and. &
            abs(x(ix^S, 2) - x2) < 0.1d0 * l2)
          w(ix^S,u_) = 0.5d0
          w(ix^S,v_) = 0.25d0
       endwhere
    case (2)
       ! Center square with random noise
       where (abs(x(ix^S, 1) - x1) < 0.1d0 * l1 .and. &
            abs(x(ix^S, 2) - x2) < 0.1d0 * l2)
          w(ix^S,u_) = 1.0d-1 * (urand - 0.5d0) + 0.5d0
          w(ix^S,v_) = 0.25d0
       endwhere
    case (3)
       ! Two Gaussians
       dist2 = (x(ix^S, 1) - x1)**2 + (x(ix^S, 2) - x2)**2
       w(ix^S,u_) = w(ix^S,u_) - 0.5d0 * exp(-100 * dist2/l1**2)
       w(ix^S,v_) = w(ix^S,v_) + 0.5d0 * exp(-100 * dist2/l1**2)

       x1 = xprobmin1 + 0.55d0 * l1
       x2 = xprobmin2 + 0.6d0 * l2
       dist2 = (x(ix^S, 1) - x1)**2 + (x(ix^S, 2) - x2)**2
       w(ix^S,u_) = w(ix^S,u_) - 0.5d0 * exp(-100 * dist2/l1**2)
       w(ix^S,v_) = w(ix^S,v_) + 0.5d0 * exp(-100 * dist2/l1**2)

    case (4)
       ! Junction of two circular shapes
       dist2 = (x(ix^S, 1) - x1)**2 + (x(ix^S, 2) - x2)**2
       mymask = (dist2 < (0.05d0*l1)**2)

       x1 = xprobmin1 + 0.55d0 * l1
       x2 = xprobmin2 + 0.6d0 * l2
       dist2 = (x(ix^S, 1) - x1)**2 + (x(ix^S, 2) - x2)**2
       mymask = mymask .or. (dist2 < (0.1d0*l1)**2)

       where (mymask)
          w(ix^S,u_) = 0.5d0
          w(ix^S,v_) = 0.25d0
       end where
    case default
       call mpistop("Unknown iprob")
    end select

  end subroutine gray_scott_init

  subroutine place_samplingpoints(n_particles, x, v, q, m, follow)
     integer, intent(in)           :: n_particles
     double precision, intent(out) :: x(3, n_particles)
     double precision, intent(out) :: v(3, n_particles)
     double precision, intent(out) :: q(n_particles)
     double precision, intent(out) :: m(n_particles)
     logical, intent(out)          :: follow(n_particles)

     integer :: nsize,i,j,i_part
     double precision :: xmid,ymid,dxx,dyy

     v = 0.0d0
     q = 0.0d0
     m = 0.0d0
     nsize=nint(dlog(dble(n_particles))/dlog(2.0d0))
     xmid=(xprobmax1-xprobmin1)*0.5d0+xprobmin1
     dxx=(xprobmax1-xprobmin1)*0.5d0/nsize
     dyy=(xprobmax2-xprobmin2)*0.5d0/nsize
     i_part=0
     do i=1,nsize
        ymid=(xprobmax2-xprobmin2)*0.5d0+xprobmin2
        do j=1,nsize
           i_part=i_part+1
           x(1,i_part)=xmid
           x(2,i_part)=ymid
           ! print *,i_part,x(1,i_part),x(2,i_part)
           follow(i_part)=.true.
           ymid=ymid+dyy
        enddo
        xmid=xmid+dxx
     enddo
     if (i_part/=n_particles) then
         print *,nsize,n_particles,i_part
         call mpistop("error in place_samplepoints")
     endif

  end subroutine place_samplingpoints

end module mod_usr
