INCLUDE:amrvacmodules/gravity.t
INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t

! Set equation parameters
subroutine initglobaldata_usr
  include 'amrvacdef.f'

  eqpar(gamma_) = 5.0d0/3.0d0
  eqpar(grav1_) = zero
  eqpar(grav2_) = -one
  {^IFTHREED
  eqpar(grav3_) = zero
  }
  {#IFDEF ISO
  eqpar(adiab_) = 1.
  }
end subroutine initglobaldata_usr

! initialize one grid
subroutine initonegrid_usr(ixG^L,ix^L,w,x)
  include 'amrvacdef.f'

  integer, intent(in)             :: ixG^L, ix^L
  double precision, intent(in)    :: x(ixG^S, 1:ndim)
  double precision, intent(inout) :: w(ixG^S, 1:nw)
  integer                         :: ix^D
  double precision                :: y0, epsilon, rhodens
  double precision                :: rholight, kx, kz, pint, dely
  logical                         :: first = .true.
  double precision, parameter     :: pi    = acos(-1.0d0)

  y0       = 0.8d0
  epsilon  = 0.05d0
  rhodens  = one
  rholight = 0.1d0
  kx       = 2 * pi
  {^IFTWOD
  kz       = zero
  }
  {^IFTHREED
  kz       = 2 * pi
  }

  if (first) then
     if (mype==0) then
        print *,'HD Rayleigh Taylor problem'
        print *,'  --assuming y ranging from 0-1!'
        print *,'  --interface y0-epsilon:', y0, epsilon
        print *,'  --density ratio:', rhodens/rholight
        print *,'  --kx:', kx
        print *,'  --kz:', kz
     end if
     first = .false.
  end if

  ! pressure at interface
  pint = one

  {^IFTWOD
  where(x(ixG^S, 2)>y0+epsilon*sin(kx*x(ixG^S, 1)))
     w(ixG^S, rho_) = rhodens
  elsewhere
     w(ixG^S, rho_) = rholight
  endwhere
  }

  {^IFTHREED
  if(kx*kz/=zero)then
     where(x(ixG^S, 2) > y0 + epsilon * sin(kx*x(ixG^S, 1)) * sin(kz*x(ixG^S, 3)))
        w(ixG^S, rho_) = rhodens
     elsewhere
        w(ixG^S, rho_) = rholight
     endwhere
  else
     if(kx==zero)then
        where(x(ixG^S, 2)>y0+epsilon*sin(kz*x(ixG^S, 3)))
           w(ixG^S, rho_) = rhodens
        elsewhere
           w(ixG^S, rho_) = rholight
        endwhere
     else
        where(x(ixG^S, 2)>y0+epsilon*sin(kx*x(ixG^S, 1)))
           w(ixG^S, rho_) = rhodens
        elsewhere
           w(ixG^S, rho_) = rholight
        endwhere
     endif
  endif
  }

  {^C&w(ixG^S, m^C_)=zero \}
  {#IFDEF ENERGY
  w(ixG^S, e_)=pint-w(ixG^S, rho_)*(x(ixG^S, 2)-y0)

  {^IFTWOD
  if(iprob==2)then
     dely = x(1, 2, 2)-x(1, 1, 2)
     w(ixGmin1:ixGmax1, 1, e_) = pint-rholight*(x(1, 1, 2)-y0)
     w(ixGmin1:ixGmax1, 2, e_) = pint-rholight*(x(1, 2, 2)-y0)
     do ix2=3, ixGmax2
        w(ixGmin1:ixGmax1, ix2, e_)= &
             -two*dely*w(ixGmin1:ixGmax1, ix2-1, rho_) &
             +w(ixGmin1:ixGmax1, ix2-2, e_)
     enddo
  endif
  }

  {^IFTHREED
  if(iprob==2)then
     dely = x(1, 2, 1, 2)-x(1, 1, 1, 2)
     w(ixGmin1:ixGmax1, 1, ixGmin3:ixGmax3, e_) = &
          pint-rholight*(x(1, 1, 1, 2)-y0)
     w(ixGmin1:ixGmax1, 2, ixGmin3:ixGmax3, e_) = &
          pint-rholight*(x(1, 2, 1, 2)-y0)
     do ix2=3, ixGmax2
        w(ixGmin1:ixGmax1, ix2, ixGmin3:ixGmax3, e_)= &
             -two*dely*w(ixGmin1:ixGmax1, ix2-1, ixGmin3:ixGmax3, rho_) &
             +w(ixGmin1:ixGmax1, ix2-2, ixGmin3:ixGmax3, e_)
     enddo
  endif
  }

  w(ixG^S, e_) = w(ixG^S, e_)/(eqpar(gamma_)-one)
  }
end subroutine initonegrid_usr

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT, qtC, x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT
subroutine specialsource(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x)
  include 'amrvacdef.f'

  integer, intent(in) :: ixI^L, ixO^L, iw^LIM
  double precision, intent(in) :: qdt, qtC, qt, x(ixI^S, 1:ndim)
  double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
  !-----------------------------------------------------------------------------

  call addsource_grav(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x)

end subroutine specialsource

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.
subroutine getdt_special(w, ixG^L, ix^L, dtnew, dx^D, x)
  include 'amrvacdef.f'

  integer, intent(in) :: ixG^L, ix^L
  double precision, intent(in) :: dx^D, x(ixG^S, 1:ndim)
  double precision, intent(inout) :: w(ixG^S, 1:nw), dtnew

  call getdt_grav(w, ixG^L, ix^L, dtnew, dx^D, x)

end subroutine getdt_special

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
subroutine specialrefine_grid(igrid, level, ixG^L, ix^L, qt, w, x, refine, coarsen)
  include 'amrvacdef.f'

  integer, intent(in) :: igrid, level, ixG^L, ix^L
  double precision, intent(in) :: qt, w(ixG^S, 1:nw), x(ixG^S, 1:ndim)
  integer, intent(inout) :: refine, coarsen

end subroutine specialrefine_grid

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)
subroutine specialvarforerrest(ixI^L, ixO^L, iflag, w, var)
  include 'amrvacdef.f'

  integer, intent(in)          :: ixI^L, ixO^L, iflag
  double precision, intent(in) :: w(ixI^S, 1:nw)
  double precision, intent(out):: var(ixG^T)

  if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')
  var(ixI^S) = zero

end subroutine specialvarforerrest

! Here one can add a steady (time-independent) potential background field
subroutine specialset_B0(ixI^L, ixO^L, x, wB0)

  include 'amrvacdef.f'

  integer, intent(in)           :: ixI^L, ixO^L
  double precision, intent(in)  :: x(ixI^S, 1:ndim)
  double precision, intent(inout) :: wB0(ixI^S, 1:ndir)
  !-----------------------------------------------------------------------------
  call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

  wB0(ixO^S, 1:ndir)=wB0(ixO^S, 1:ndir) ! + user defined steady potential field

end subroutine specialset_B0
