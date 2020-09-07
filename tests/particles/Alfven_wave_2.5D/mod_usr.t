module mod_usr
  use mod_mhd
  use mod_particles

  implicit none

  double precision :: charge = 1.0d0
  double precision :: mass = 1.0d0

  ! Initial position
  double precision :: x0(3) = [0.5d0, 0.25d0, 0.0d0]

  ! Initial velocity
  double precision :: v0(3) = [0.0d0, 0.0d0, 0.0d0]

  ! Maxwellian velocity (per component vx, vy, vz)
  double precision :: maxwellian_velocity = 0.1d0

  ! Use an analytic field instead of an interpolated one
  logical :: use_analytic_field = .false.

  integer :: Tep ! Index for the temperature store in gridvars

contains

  subroutine usr_init()
    use mod_initialize

    unit_length        = 1.d0
    unit_numberdensity = 1.d0
    unit_velocity      = 1.d0

    usr_init_one_grid => initonegrid_usr
    usr_create_particles => generate_particles

    particles_define_additional_gridvars => define_additional_gridvars_usr
    particles_fill_additional_gridvars => fill_additional_gridvars_usr
    usr_update_payload => update_payload_usr

    call set_coordinate_system("Cartesian_2.5D")
    call mhd_activate()

    call initialize_amrvac()    ! So that we have settings available

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision       :: A, phi(ixG^S), ca, alpha, beta
    double precision       :: v(ixG^S,1:3), b(ixG^S,1:3), vrot(ixG^S,1:3), brot(ixG^S,1:3), k(1:^ND)
    logical, save :: first=.true.

    A     = 0.1d0
    ca    = 1.0d0
    k(1)  = 2.0d0 * dpi

    {^IFONED
    alpha = 0.0d0
    beta  = 0.0d0
    phi(ix^S) = k(1) * x(ix^S,1)
    }

    {^IFTWOD 
    alpha = atan(2.0d0)
    beta  = 0.0d0
    k(2)  = k(1) * tan(alpha)
    phi(ix^S) = k(1) *x(ix^S,1) + k(2) * x(ix^S,2)
    }

    {^IFTHREED
    alpha = atan(2.0d0)
    beta  = atan(2.0d0)
    k(2)  = k(1) * tan(alpha)
    k(3)  = k(1) * tan(beta)
    phi(ix^S) = k(1) *x(ix^S,1) + k(2) * x(ix^S,2) + k(3) * x(ix^S,3)
    }

    w(ix^S,rho_) = 1.0d0
    w(ix^S,p_)   = 0.1d0
    
    v(ix^S,1)  = 0.0d0
    v(ix^S,2)  = A * sin(phi(ix^S))
    v(ix^S,3)  = A * cos(phi(ix^S))
    
    b(ix^S,1)  =   sqrt(w(ix^S,rho_)) * ca
    b(ix^S,2)  = - sqrt(w(ix^S,rho_)) * ca * A * sin(phi(ix^S))
    b(ix^S,3)  = - sqrt(w(ix^S,rho_)) * ca * A * cos(phi(ix^S))

    call rotate(ixG^L,ix^L,v,vrot,alpha,beta)
    call rotate(ixG^L,ix^L,b,brot,alpha,beta)

    w(ix^S,mom(:)) = vrot(ix^S,:)
    w(ix^S,mag(:)) = brot(ix^S,:)

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

  subroutine rotate(ixG^L,ix^L,v,vrot,alpha,beta)
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: v(ixG^S,1:3)
    double precision, intent(out) :: vrot(ixG^S,1:3)
    double precision, intent(in) :: alpha, beta

    double precision       :: gamma

    gamma = atan(cos(alpha)*tan(beta))

    vrot(ix^S,1) = cos(gamma) * cos(alpha) * v(ix^S,1) - sin(alpha) * v(ix^S,2) - sin(gamma) * cos(alpha) * v(ix^S,3)
    vrot(ix^S,2) = cos(gamma) * sin(alpha) * v(ix^S,1) + cos(alpha) * v(ix^S,2) - sin(gamma) * sin(alpha) * v(ix^S,3)
    vrot(ix^S,3) = sin(gamma)              * v(ix^S,1)                          + cos(gamma) * v(ix^S,3)

  end subroutine rotate

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
    double precision              :: tmp_vec(4), phi, rrd(ndir)
    integer                       :: idir, n

    q = charge
    m = mass

    ! Randomly distributed
    do idir=1,ndir
      rrd(idir) = rng%unif_01()
    end do
    {^D&x(^D) = xprobmin^D + rrd(^D) * (xprobmax^D - xprobmin^D)\}

    ! Add Maxwellian velocity. Random numbers come in pairs of two
    tmp_vec(1:2) = rng%two_normals()
    tmp_vec(3:4) = rng%two_normals()
    v = v0 + tmp_vec(1:3) * maxwellian_velocity
    print *, "Particle", ipart, v

  end subroutine get_particle

  subroutine update_payload_usr(igrid,w,wold,xgrid,xpart,upart,qpart,mpart,mypayload,mynpayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,mynpayload
    double precision, intent(in)  :: w(ixG^T,1:nw),wold(ixG^T,1:nw)
    double precision, intent(in)  :: xgrid(ixG^T,1:ndim),xpart(1:ndir),upart(1:ndir),qpart,mpart,particle_time
    double precision, intent(out) :: mypayload(mynpayload)

    ! Temperature
    if (npayload > 0) then
      ! Interpolate temperature here
      call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,Tep),xgrid,xpart,mypayload(1))
    end if

  end subroutine update_payload_usr

  ! Custom particle gridvars definition
  subroutine define_additional_gridvars_usr(ngridvars)
    use mod_global_parameters
    integer, intent(inout) :: ngridvars
 
    ! Temperature
    ngridvars = ngridvars + 1
    Tep = ngridvars

  end subroutine define_additional_gridvars_usr

  ! Custom particle gridvars filler
  subroutine fill_additional_gridvars_usr
    use mod_global_parameters
    use mod_usr_methods, only: usr_particle_fields

    integer :: igrid, iigrid
    double precision :: pth(ixG^T)
    double precision :: w(ixG^T,1:nw)

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      w(ixG^T,1:nw) = ps(igrid)%w(ixG^T,1:nw)

      ! Fill temperature
      call mhd_get_pthermal(w,ps(igrid)%x,ixG^LL,ixG^LL,pth)
      gridvars(igrid)%w(ixG^T,Tep) = pth(ixG^T)/w(ixG^T,rho_)
      
      ! Same for old state vector if interpolating in time
      if (time_advance) then 
        w(ixG^T,1:nw) = pso(igrid)%w(ixG^T,1:nw) 

        ! Temperature
        call mhd_get_pthermal(w,pso(igrid)%x,ixG^LL,ixG^LL,pth)
        gridvars(igrid)%wold(ixG^T,Tep) = pth(ixG^T)/w(ixG^T,rho_)

      end if

    end do
  end subroutine fill_additional_gridvars_usr

end module mod_usr
