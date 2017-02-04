module mod_gravity

contains

  !> Add sources from gravity
  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine gravity_add_source(qdt,ixI^L,ixO^L,rho_,mom,e_,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    integer, intent(in)             :: rho_, mom(ndir), e_
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer                         :: idim

    double precision :: gravity_field(ixI^S,ndim)

    call usr_gravity(ixI^L,ixO^L,wCT,x,gravity_field)

    do idim = 1, ndim
      w(ixO^S,mom(idim)) = w(ixO^S,mom(idim)) &
            + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,rho_)
      if(e_>0) then
        w(ixO^S,e_)=w(ixO^S,e_) &
            + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,mom(idim))
      end if
    end do

  end subroutine gravity_add_source

  subroutine gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew
    double precision                :: dxinv(1:ndim), dtgrav
    integer                         :: idim

    double precision :: gravity_field(ixI^S,ndim)

    ^D&dxinv(^D)=one/dx^D;

    dtgrav = bigdouble

    call usr_gravity(ixI^L,ixO^L,w,x,gravity_field)

    do idim = 1, ndim
      dtgrav = min(dtgrav, 1.0d0 / &
           sqrt(abs(gravity_field(ixO^S,idim)) * dxinv(idim)))
    end do

    dtnew = dtgrav

  end subroutine gravity_get_dt

end module mod_gravity
