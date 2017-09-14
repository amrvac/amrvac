module mod_usr
  use mod_hd

  implicit none

  ! Rate constant
  double precision :: k1 = 1.0d0

  double precision :: density_left = 0.8d0
  double precision :: density_right = 1.0d0

  integer :: i_gamma
  integer :: i_temperature

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_variables

    usr_init_one_grid => initonegrid_usr
    usr_source        => do_reactions
    usr_process_grid  => set_gamma
    hd_usr_gamma      => get_gamma

    call set_coordinate_system("Cartesian_2D")
    call hd_activate()

    ! Rename tracers
    prim_wnames(tracer(1)) = "H2"
    prim_wnames(tracer(2)) = "O2"

    ! Add own variables
    i_gamma = var_set_extravar("gamma", "gamma")
    i_temperature = var_set_extravar("temperature", "temperature")

  end subroutine usr_init

  subroutine set_gamma(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! Gamma also needs to be set in ghostcells, because 
    w(ixI^S, i_gamma) = hd_gamma - global_time * 0.1
    w(ixI^S, i_temperature) = sin(x(ixI^S, 1) * x(ixI^S, 2))
  end subroutine set_gamma

  subroutine get_gamma(w, ixI^L, ixO^L, gam)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision, intent(out) :: gam(ixO^S)

    if (minval(w(ixO^S, i_gamma)) < maxval(w(ixO^S, i_gamma))) &
         error stop "gamma problem"

    gam(ixO^S) = w(ixO^S, i_gamma)
  end subroutine get_gamma

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
       w(ixO^S,e_) = 0.1d0 * w(ixO^S,rho_)
    end if

    ! Set chemical species (as tracers)
    do n = 1, hd_n_tracer
       w(ixO^S, tracer(n)) = w(ixO^S,rho_) / n
    end do

    w(ixI^S, i_gamma) = hd_gamma

    call hd_to_conserved(ixI^L, ixO^L, w, x)
  end subroutine initonegrid_usr

end module mod_usr
