module mod_gravity

  !> The gravity field (unused dimensions are ignored)
  double precision :: gravity_field(3) = 0.0d0

contains

  !> Add sources from gravity
  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine gravity_add_source(qdt,ixI^L,ixO^L,iw^LIM, &
       rho_,mom,e_,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    integer, intent(in)             :: rho_, mom(ndir), e_
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer                         :: iw, idim

    do idim = 1, ndim
       if (mom(idim) >= iwmin .and. mom(idim) < iwmax) then
          w(ixO^S,mom(idim)) = w(ixO^S,mom(idim)) &
               + qdt * gravity_field(idim) * wCT(ixO^S,rho_)
       end if

       if (e_ >= iwmin .and. e_ <= iwmax) then
          w(ixO^S,e_)=w(ixO^S,e_) &
               + qdt * gravity_field(idim) * wCT(ixO^S,mom(idim))
       end if
    end do

  end subroutine gravity_add_source

  subroutine gravity_get_dt(w,ixG^L,ix^L,dtnew,dx^D,x)
    use mod_global_parameters

    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: dx^D, x(ixG^S,1:ndim), w(ixG^S,1:nw)
    double precision, intent(inout) :: dtnew
    double precision                :: dxinv(1:ndim), dtgrav
    integer                         :: idim

    ^D&dxinv(^D)=one/dx^D;

    dtgrav = bigdouble

    do idim = 1, ndim
       if (abs(gravity_field(idim))>zero) then
          dtgrav = min(dtgrav, 1.0d0 / &
               sqrt(abs(gravity_field(idim)) * dxinv(idim)))
       end if
    enddo

    dtnew = dtgrav

  end subroutine gravity_get_dt

end module mod_gravity
