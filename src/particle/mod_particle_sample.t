!> Scattered sampling based on fixed-particle interpolation
!> By Fabio Bacchini (2020)
module mod_particle_sample
  use mod_particle_base

  private

  !> Variable index for velocity
  integer, dimension(:), allocatable      :: vp(:)

  public :: sample_init
  public :: sample_create_particles

contains

  subroutine sample_init()
    use mod_global_parameters
    integer :: idir

    ngridvars=nw

    allocate(vp(ndir))
    do idir = 1, ndir
      vp(idir) = idir
    end do

    if (.not. associated(particles_fill_gridvars)) &
         particles_fill_gridvars => sample_fill_gridvars
    particles_integrate     => sample_integrate_particles

  end subroutine sample_init

  subroutine sample_create_particles()
    ! initialise the particles (=fixed interpolation points)
    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles, usr_update_payload

    integer          :: n, idir, igrid, ipe_particle
    double precision :: x(3, num_particles)
    double precision :: v(3, num_particles)
    double precision :: q(num_particles)
    double precision :: m(num_particles)
    double precision :: rrd(num_particles,ndir)
    double precision :: w(ixG^T,1:nw)
    double precision :: payload(npayload)
    logical          :: follow(num_particles)

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

    nparticles = num_particles

    do n=1,num_particles
      call find_particle_ipe(x(:,n),igrid,ipe_particle)
      particle(n)%igrid  = igrid
      particle(n)%ipe    = ipe_particle

      if(ipe_particle == mype) then
        call push_particle_into_particles_on_mype(n)
        allocate(particle(n)%self)
        particle(n)%self%follow = follow(n)
        particle(n)%self%index  = n
        particle(n)%self%time   = global_time
        particle(n)%self%dt     = 0.0d0
        particle(n)%self%x = 0.d0
        particle(n)%self%x(:) = x(:,n)
        particle(n)%self%u(:) = 0.d0
        allocate(particle(n)%payload(npayload))
        if (.not. associated(usr_update_payload)) then
          call sample_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x(:,n),v(:,n),q(n),m(n),payload,npayload,0.d0)
        else
          call usr_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x(:,n),v(:,n),q(n),m(n),payload,npayload,0.d0)
        end if
        particle(n)%payload=payload
      end if

    end do

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
    double precision :: payload(npayload)
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

      ! Velocity update
!      call get_vec_sample(igrid,x,tlocnew,v,vp(1),vp(ndir))
!      particle(ipart)%self%u(1:ndir) = v(1:ndir)

      ! Time update
      particle(ipart)%self%time = tlocnew

      ! Update payload
      if (.not. associated(usr_update_payload)) then
        call sample_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x,v,0.d0,0.d0,payload,npayload,tlocnew)
      else
        call usr_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x,v,0.d0,0.d0,payload,npayload,tlocnew)
      end if
      particle(ipart)%payload = payload

    end do

  end subroutine sample_integrate_particles

  !> Payload update
  subroutine sample_update_payload(igrid,w,wold,xgrid,xpart,upart,qpart,mpart,payload,npayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,npayload
    double precision, intent(in)  :: w(ixG^T,1:nw),wold(ixG^T,1:nw)
    double precision, intent(in)  :: xgrid(ixG^T,1:ndim),xpart(1:ndir),upart(1:ndir),qpart,mpart,particle_time
    double precision, intent(out) :: payload(npayload)
    double precision              :: wp, wp1, wp2, td
    integer                       :: ii

    td = (particle_time - global_time) / dt

    ! There are npayload=nw payloads, one for each primitive fluid quantity
    do ii=1,npayload
      if (.not.time_advance) then
        call interpolate_var(igrid,ixG^LL,ixM^LL,w(ixG^T,ii),xgrid,xpart,wp)
      else
        call interpolate_var(igrid,ixG^LL,ixM^LL,wold(ixG^T,ii),xgrid,xpart,wp1)
        call interpolate_var(igrid,ixG^LL,ixM^LL,w(ixG^T,ii),xgrid,xpart,wp2)
        wp = wp1 * (1.0d0 - td) + wp2 * td
      end if
      payload(ii) = wp
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

  subroutine get_vec_sample(igrid,x,tloc,var,ibeg,iend)
    use mod_global_parameters

    integer,intent(in)                                   :: igrid, ibeg, iend
    double precision,dimension(ndir), intent(in)         :: x
    double precision, intent(in)                         :: tloc
    double precision,dimension(iend-ibeg+1), intent(out) :: var
    double precision,dimension(iend-ibeg+1)              :: e1, e2
    integer                                              :: ivar, iloc
    double precision                                     :: td

    if(.not.time_advance) then
      do ivar=ibeg,iend
        iloc = ivar-ibeg+1
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ivar),ps(igrid)%x(ixG^T,1:ndim),x,var(iloc))
      end do
    else
      td = (tloc - global_time) / dt
      do ivar=ibeg,iend
        iloc = ivar-ibeg+1
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%wold(ixG^T,ivar),ps(igrid)%x(ixG^T,1:ndim),x,e1(iloc))
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ivar),ps(igrid)%x(ixG^T,1:ndim),x,e2(iloc))
        var(iloc) = e1(iloc) * (1.0d0 - td) + e2(iloc) * td
      end do
    end if

  end subroutine get_vec_sample

end module mod_particle_sample
