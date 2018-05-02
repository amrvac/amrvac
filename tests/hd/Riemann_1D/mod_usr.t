module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => rm1d_init_one_grid

    call set_coordinate_system("Cartesian")
    call hd_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm1d_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    ! iprob==1 rarefaction wave & shock
    if (iprob==1) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.5d0
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 1.5d0
        elsewhere
           w(ix^S,rho_)   = 0.5d0
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 0.5d0
        end where
    ! iprob==2  shock & shock
    else if (iprob==2) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 2.0d0
           w(ix^S,e_)     = 1.5d0
        elsewhere
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = -2.0d0
           w(ix^S,e_)     = 1.5d0
        end where
    ! iprob==3  rarefaction wave
    else if (iprob==3) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = -0.5d0
           w(ix^S,e_)     = 1.0d0
        elsewhere
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 0.5d0
           w(ix^S,e_)     = 1.0d0
        end where
    else
        call mpistop("iprob not available!")
    end if

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine rm1d_init_one_grid

end module mod_usr
