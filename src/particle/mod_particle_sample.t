!> Scattered sampling based on fixed-particle interpolation
!> By Fabio Bacchini (2020)
module mod_particle_sample
  use mod_particle_base

  private

  public :: sample_init
  public :: sample_create_particles

contains

  subroutine sample_init()
    use mod_global_parameters
    integer :: idir

    allocate(vp(ndir))
    do idir = 1, ndir
      vp(idir) = idir
    end do
    ngridvars=nw

    particles_fill_gridvars => sample_fill_gridvars

    if (associated(particles_define_additional_gridvars)) then
      call particles_define_additional_gridvars(ngridvars)
    end if

    particles_integrate     => sample_integrate_particles

  end subroutine sample_init

  subroutine sample_create_particles()
    ! initialise the particles (=fixed interpolation points)
    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles, usr_update_payload, usr_check_particle

    integer          :: n, idir, igrid, ipe_particle, nparticles_local
    double precision :: x(3, num_particles)
    double precision :: v(3, num_particles)
    double precision :: q(num_particles)
    double precision :: m(num_particles)
    double precision :: rrd(num_particles,ndir)
    double precision :: w(ixG^T,1:nw)
    double precision :: defpayload(ndefpayload)
    double precision :: usrpayload(nusrpayload) 
    logical          :: follow(num_particles), check

    follow = .false.
    x      = 0.0d0

    if (mype==0) then
      if (.not. associated(usr_create_particles)) then
        ! Randomly distributed
        do idir=1,ndir
          do n = 1, num_particles
            rrd(n,idir) = rng%unif_01()
          end do
        end do
        do n=1, num_particles
          {^D&x(^D,n) = xprobmin^D + rrd(n+1,^D) * (xprobmax^D - xprobmin^D)\}
        end do
      else
        call usr_create_particles(num_particles, x, v, q, m, follow)
      end if
    end if

    call MPI_BCAST(x,3*num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(follow,num_particles,MPI_LOGICAL,0,icomm,ierrmpi)

    nparticles_local = 0

    do n=1,num_particles
      call find_particle_ipe(x(:,n),igrid,ipe_particle)
      particle(n)%igrid  = igrid
      particle(n)%ipe    = ipe_particle

      if(ipe_particle == mype) then
        check = .true.

        ! Check for user-defined modifications or rejection conditions
        if (associated(usr_check_particle)) call usr_check_particle(igrid, x(:,n), v(:,n), q(n), m(n), follow(n), check)
        if (check) then
          call push_particle_into_particles_on_mype(n)
        else
          cycle
        end if

        nparticles_local = nparticles_local + 1

        allocate(particle(n)%self)
        particle(n)%self%follow = follow(n)
        particle(n)%self%index  = n
        particle(n)%self%time   = global_time
        particle(n)%self%dt     = 0.0d0
        particle(n)%self%x = 0.d0
        particle(n)%self%x(:) = x(:,n)
        particle(n)%self%u(:) = 0.d0
        allocate(particle(n)%payload(npayload))
        call sample_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x(:,n),v(:,n),q(n),m(n),defpayload,ndefpayload,0.d0)
        particle(n)%payload(1:ndefpayload) = defpayload      
        if (associated(usr_update_payload)) then
          call usr_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x(:,n),v(:,n),q(n),m(n),usrpayload,nusrpayload,0.d0)
          particle(n)%payload(ndefpayload+1:npayload)=usrpayload
        end if
      end if

    end do

    call MPI_ALLREDUCE(nparticles_local,nparticles,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)

  end subroutine sample_create_particles

  subroutine sample_fill_gridvars
    use mod_global_parameters

    integer                                   :: igrid, iigrid, idir
    double precision, dimension(ixG^T,1:nw)   :: w

    do iigrid=1,igridstail; igrid=igrids(iigrid);

      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

      gridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0
      w(ixG^T,1:nw) = ps(igrid)%w(ixG^T,1:nw)
      call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
      ! fill all variables:
      gridvars(igrid)%w(ixG^T,1:ngridvars) = w(ixG^T,1:ngridvars)

      if(time_advance) then
        gridvars(igrid)%wold(ixG^T,1:ngridvars) = 0.0d0
        w(ixG^T,1:nw) = pso(igrid)%w(ixG^T,1:nw)
        call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
        gridvars(igrid)%wold(ixG^T,1:ngridvars) = w(ixG^T,1:ngridvars)
      end if

    end do

  end subroutine sample_fill_gridvars

  subroutine sample_integrate_particles(end_time)
    ! this interpolates the HD/MHD quantities at the particle positions
    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles, usr_update_payload
    double precision, intent(in) :: end_time

    double precision, dimension(1:ndir) :: v, x
    double precision                 :: defpayload(ndefpayload)
    double precision                 :: usrpayload(nusrpayload)
    double precision                 :: tloc, tlocnew, dt_p, h1
    double precision,parameter       :: eps=1.0d-6, hmin=1.0d-8
    integer                          :: ipart, iipart, igrid
    integer                          :: nok, nbad, ierror

    do iipart=1,nparticles_active_on_mype
      ipart                   = particles_active_on_mype(iipart);
      dt_p                    = sample_get_particle_dt(particle(ipart), end_time)
      particle(ipart)%self%dt = dt_p

      igrid                   = particle(ipart)%igrid
      igrid_working           = igrid
      tloc                    = particle(ipart)%self%time
      x(1:ndir)               = particle(ipart)%self%x(1:ndir)
      tlocnew                 = tloc+dt_p

      ! Time update
      particle(ipart)%self%time = tlocnew

      ! Update payload
      call sample_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x,v,0.d0,0.d0,defpayload,ndefpayload,tlocnew)
      particle(ipart)%payload(1:ndefpayload) = defpayload
      if (associated(usr_update_payload)) then
        call usr_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x,v,0.d0,0.d0,usrpayload,nusrpayload,tlocnew)
        particle(ipart)%payload(ndefpayload+1:npayload) = usrpayload
      end if

    end do

  end subroutine sample_integrate_particles

  !> Payload update
  subroutine sample_update_payload(igrid,w,wold,xgrid,xpart,upart,qpart,mpart,mypayload,mynpayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,mynpayload
    double precision, intent(in)  :: w(ixG^T,1:nw),wold(ixG^T,1:nw)
    double precision, intent(in)  :: xgrid(ixG^T,1:ndim),xpart(1:ndir),upart(1:ndir),qpart,mpart,particle_time
    double precision, intent(out) :: mypayload(mynpayload)
    double precision              :: wp, wp1, wp2, td
    integer                       :: ii

    td = (particle_time - global_time) / dt

    ! There are npayload=nw payloads, one for each primitive fluid quantity
    do ii=1,mynpayload
      if (.not.time_advance) then
        call interpolate_var(igrid,ixG^LL,ixM^LL,w(ixG^T,ii),xgrid,xpart,wp)
      else
        call interpolate_var(igrid,ixG^LL,ixM^LL,wold(ixG^T,ii),xgrid,xpart,wp1)
        call interpolate_var(igrid,ixG^LL,ixM^LL,w(ixG^T,ii),xgrid,xpart,wp2)
        wp = wp1 * (1.0d0 - td) + wp2 * td
      end if
      mypayload(ii) = wp
    end do

  end subroutine sample_update_payload

  pure function sample_get_particle_dt(partp, end_time) result(dt_p)
    use mod_global_parameters
    use mod_geometry
    type(particle_ptr), intent(in) :: partp
    double precision, intent(in)   :: end_time
    double precision               :: dt_p
    integer                        :: ipart, iipart, nout
    double precision               :: tout, dt_cfl
    double precision               :: v(1:ndir)

    if (const_dt_particles > 0) then
      dt_p = const_dt_particles
      return
    end if

    dt_p = dtsave_particles

    ! Make sure we don't advance beyond end_time
    call limit_dt_endtime(end_time - partp%self%time, dt_p)

  end function sample_get_particle_dt

end module mod_particle_sample
