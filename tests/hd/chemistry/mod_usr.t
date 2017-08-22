module mod_usr
  use mod_hd

  implicit none

  ! Rate constant
  double precision :: k1 = 1.0d0

  double precision :: density_left = 0.8d0
  double precision :: density_right = 1.0d0

contains

  subroutine usr_init()
    use mod_usr_methods

    usr_init_one_grid => initonegrid_usr
    usr_source        => do_reactions
    hd_usr_gamma      => set_gamma

    call set_coordinate_system("Cartesian_2D")
    call hd_activate()

  end subroutine usr_init

  subroutine set_gamma(w, ixI^L, ixO^L, gam)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision, intent(out) :: gam(ixO^S)

    ! Set gamma to demonstrate the functionality; no physical meaning here
    gam(ixO^S) = 1 + w(ixO^S, tracer(1))/(1 + w(ixO^S, tracer(2)))
  end subroutine set_gamma

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine do_reactions(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixI^S, tracer(1)) = w(ixI^S, tracer(1)) &
         + k1 * qdt * wCT(ixI^S, tracer(2))

    w(ixI^S, tracer(2)) = w(ixI^S, tracer(2)) &
         - k1 * qdt * wCT(ixI^S, tracer(2))
  end subroutine do_reactions

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer                         :: n

    ! Initialize the density
    where(x(ixO^S,1) > 0.5d0 * (xprobmin1 + xprobmax1))
       w(ixO^S,rho_) = density_left
    elsewhere
       w(ixO^S,rho_) = density_right
    endwhere

    ! Set velocities
    w(ixO^S, mom(:)) = 0.0d0

    ! Pressure
    if (hd_energy) then
       w(ixO^S,e_) = 1.0d0 * w(ixO^S,rho_)
    end if

    ! Set chemical species (as tracers)
    do n = 1, hd_n_tracer
       w(ixO^S, tracer(n)) = w(ixO^S,rho_) / n
    end do

    call hd_to_conserved(ixI^L, ixO^L, w, x)
  end subroutine initonegrid_usr

end module mod_usr
