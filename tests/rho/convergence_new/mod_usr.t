module mod_usr
  use mod_rho

  implicit none

contains

  subroutine usr_init()
    use mod_usr_methods

    usr_init_one_grid => initonegrid_usr

    call rho_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid 

    use mod_global_parameters

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    !----------------------------------------------------------------------------
    select case (iprob)
    case (1)
       ! Advection of \sin(\pi x)**4
       w(ix^S,rho_) = sin(dpi*x(ix^S,1))**4.0d0
    case default
       call mpistop("iprob not available!")
    end select


  end subroutine initonegrid_usr

end module mod_usr

