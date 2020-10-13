module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()

    call set_coordinate_system("Cartesian")

    usr_init_one_grid => liska_init_one_grid

    call hd_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine liska_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    integer :: itr
    double precision                :: front, rho0, rho1, pr0, pr1

    {^IFONED call mpistop("This is a multi-D HD problem") }

    front  = 0.15d0
    rho0   = 1.0d0
    rho1   = 0.125d0
    pr0    = 1.0d0
    pr1    = 0.14d0

    where(( ^D&x(ix^S,^D)+ ) >= front+smalldouble)
       w(ix^S,rho_)           = rho0
       w(ix^S,mom(1))         = zero
       w(ix^S,mom(2))         = zero
    elsewhere
       w(ix^S,rho_)           = rho1
       w(ix^S,mom(1))         = zero
       w(ix^S,mom(2))         = zero
    endwhere

       if(hd_energy)then
          where(( ^D&x(ix^S,^D)+ ) >= front+smalldouble)
             w(ix^S,e_)             = pr0/(hd_gamma-one)
          elsewhere
             w(ix^S,e_)             = pr1/(hd_gamma-one)
          endwhere
       endif

    do itr = 1, hd_n_tracer
       where(( ^D&x(ix^S,^D)+ ) >= front+smalldouble)
          w(ix^S,tracer(itr))     = one
       elsewhere
          w(ix^S,tracer(itr))     = zero
      endwhere
    end do

    call hd_to_conserved(ixG^L,ix^L,w,x)

    end subroutine liska_init_one_grid

end module mod_usr
