module mod_usr
  use mod_rho

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => initial_conditions
    usr_refine_grid => specialrefine_grid

    call set_coordinate_system("Cartesian")
    call rho_activate()
  end subroutine usr_init

  subroutine initial_conditions(ixG^L, ix^L, w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, 1:ndim)
    double precision, intent(inout) :: w(ixG^S, 1:nw)

    w(ix^S, rho_) = 0.0d0
  end subroutine initial_conditions

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision:: xmid^D, xxr(ixG^S)
    logical:: mask(ixG^S)

    ^D&xmid^D=xprobmin^D+0.5d0*(xprobmax^D-xprobmin^D);
    xxr(ix^S)=dsqrt((x(ix^S,1)-xmid1)**2)
    ! test with different levels of refinement enforced
    if (any(xxr(ix^S)<0.025d0)) then
       refine=1
    endif
  end subroutine specialrefine_grid

end module mod_usr
