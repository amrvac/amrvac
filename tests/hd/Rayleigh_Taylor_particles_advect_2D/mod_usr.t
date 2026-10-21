module mod_usr
  use mod_hd
  use mod_eos, only: eos
  use mod_particles

  implicit none

  integer :: isol   ! index for extra payload in gridvars for particles

contains

  subroutine usr_init()
    use mod_variables

    usr_init_one_grid => initonegrid_usr
    usr_gravity=> gravity

    ! Particle stuff
    usr_create_particles => generate_particles
    usr_particle_position => move_particle
    particles_define_additional_gridvars => define_additional_gridvars_usr
    particles_fill_additional_gridvars => fill_additional_gridvars_usr
    usr_update_payload => update_payload_usr

    call set_coordinate_system("Cartesian_2D")

    call hd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ix^D

    double precision:: y0,epsilon,rhodens,rholight,kx,pint
    logical::          first
    data first/.true./


    ! the location of demarcation line`
    y0=0.8d0

    ! density of two types
    rhodens=one
    rholight=0.1d0

    ! setup the perturbation
    epsilon=0.05d0
    ! kx=2 pi
    kx=8.0d0*atan(one)

    ! print out the info
    if (first) then
       if (mype==0) then
          print *,'HD Rayleigh Taylor problem'
          print *,'  --assuming y ranging from 0-1!'
          print *,'  --interface y0-epsilon:',y0,epsilon
          print *,'  --density ratio:',rhodens/rholight
          print *,'  --kx:',kx
       end if
       first=.false.
    end if

    ! initialize the density
    where(x(ixO^S,2)>y0+epsilon*sin(kx*x(ixO^S,1)))
      w(ixO^S,rho_)=rhodens
    elsewhere
      w(ixO^S,rho_)=rholight
    endwhere

    ! set all velocity to zero
    w(ixO^S, mom(:)) = zero

    ! pressure at interface
    pint=one
    if(hd_energy) then
      w(ixO^S,e_)=pint-w(ixO^S,rho_)*(x(ixO^S,2)-y0)
      w(ixO^S,e_)=w(ixO^S,e_)/(eos%gamma-one)
    end if
  end subroutine initonegrid_usr

  ! Calculate gravitational acceleration in each dimension
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:)=0.d0
    gravity_field(ixO^S,2)=-1.d0

  end subroutine gravity


  subroutine generate_particles(n_particles, x, v, q, m, follow)
    use mod_particles
    integer, intent(in)           :: n_particles
    double precision, intent(out) :: x(3, n_particles)
    double precision, intent(out) :: v(3, n_particles)
    double precision, intent(out) :: q(n_particles)
    double precision, intent(out) :: m(n_particles)
    logical, intent(out)          :: follow(n_particles)
    integer                       :: n

    do n = 1, n_particles
      call get_particle(x(:, n), v(:, n), q(n), m(n), n, n_particles)
    end do

    follow(:) = .true.

  end subroutine generate_particles

  ! Set particle properties (in code units)
  subroutine get_particle(x, v, q, m, ipart, n_particles)
    double precision, intent(out) :: x(3), v(3), q, m
    integer, intent(in)           :: ipart, n_particles

    q = 0.d0
    m = 0.d0
    v(:) = 0.d0

    ! 3 particles in a line
    x(1) = xprobmin1 + (xprobmax1 - xprobmin1)/4.d0*DBLE(ipart)
    x(2) = xprobmin2 + (xprobmax2 - xprobmin2)/2.d0
    if(mype==0)print *, "Particle", ipart, x

  end subroutine get_particle

  ! User-defined particle movement
  subroutine move_particle(x, n, told, tnew)
    double precision, intent(inout) :: x(3)
    double precision, intent(in)  :: told,tnew
    integer, intent(in)           :: n
    double precision              :: v_move
    double precision              :: rad, th, xc,yc, freq

    v_move = 0.25d0
    rad = 0.0398d0
    xc = (xprobmax1 - xprobmin1)*.75d0-rad
    yc = (xprobmax2 - xprobmin2)/2.d0
    freq = v_move/rad

    ! particle 1 does not move;
    ! particle 2 moves on a straight line with vy=v_move
    ! particle 3 moves in a circle with vtang=v_move
    if (n .eq. 2) x(2) = x(2) + v_move*(tnew-told)
    if (n .eq. 3) then
      th = atan2(x(2)-yc,x(1)-xc)
      th = th + freq*(tnew-told)
      x(1) = xc+rad*cos(th)
      x(2) = yc+rad*sin(th)
    end if

  end subroutine move_particle

  ! Custom particle gridvars definition
  subroutine define_additional_gridvars_usr(ngridvars)
    use mod_global_parameters
    integer, intent(inout) :: ngridvars

    ! extra variable as payload: add the actual solution
    isol = ngridvars+1
    ngridvars = ngridvars+1

  end subroutine define_additional_gridvars_usr

  ! Custom particle gridvars filler
  subroutine fill_additional_gridvars_usr
    use mod_global_parameters
    use mod_usr_methods, only: usr_particle_fields

    integer :: igrid, iigrid
    double precision :: logrhoprofile(ixG^T)

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      logrhoprofile(ixG^T)=dlog10(ps(igrid)%w(ixG^T,rho_))
      gridvars(igrid)%w(ixG^T,isol) = logrhoprofile(ixG^T)
    end do

  end subroutine fill_additional_gridvars_usr

  subroutine update_payload_usr(igrid,xpart,upart,qpart,mpart,payload,npayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,npayload
    double precision, intent(in)  :: xpart(1:ndir),upart(1:ndir),qpart,mpart,particle_time
    double precision, intent(out) :: payload(npayload)

    ! put the log of density
    if (npayload > 0) then
      call interpolate_var(igrid,ixG^LL,ixG^LL,gridvars(igrid)%w(ixG^T,isol),ps(igrid)%x(ixG^T,1:ndim),xpart,payload(1))
    end if

  end subroutine update_payload_usr


end module mod_usr
