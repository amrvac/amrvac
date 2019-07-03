!> tracer particles moving with fluid flows
module mod_particle_advect
  use mod_particle_base

  private

  !> Variable index for velocity
  integer, dimension(:), allocatable      :: vp(:)

  public :: advect_init
  public :: advect_create_particles

contains

  subroutine advect_init()
    use mod_global_parameters
    integer :: idir

    ngridvars=ndir

    allocate(vp(ndir))
    do idir = 1, ndir
      vp(idir) = idir
    end do

    if (.not. associated(particles_fill_gridvars)) &
         particles_fill_gridvars => advect_fill_gridvars
    particles_integrate     => advect_integrate_particles

  end subroutine advect_init

  subroutine advect_create_particles()
    ! initialise the particles
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

    do n=1,num_particles
      call find_particle_ipe(x(:,n),igrid,ipe_particle)
      particle(n)%igrid  = igrid
      particle(n)%ipe    = ipe_particle

      if(ipe_particle == mype) then
        call push_particle_into_particles_on_mype(n)
        allocate(particle(n)%self)
        particle(n)%self%follow = follow(n)
        particle(n)%self%index  = n
        particle(n)%self%t      = 0.0d0
        particle(n)%self%dt     = 0.0d0
        particle(n)%self%x = 0.d0
        particle(n)%self%x(:) = x(:,n)
        w=ps(igrid)%w
        call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
        do idir=1,ndir
          call interpolate_var(igrid,ixG^LL,ixM^LL,&
               w(ixG^T,iw_mom(idir)),ps(igrid)%x,x(:,n),v(idir,n))
        end do
        particle(n)%self%u(:) = 0.d0
        particle(n)%self%u(1:ndir) = v(1:ndir,n)
        allocate(particle(n)%payload(npayload))
        if (.not. associated(usr_update_payload)) then
          call advect_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x(:,n),payload,npayload,0.d0)
        else
          call usr_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x(:,n),payload,npayload,0.d0)
        end if
        particle(n)%payload=payload
      end if

    end do

  end subroutine advect_create_particles

  subroutine advect_fill_gridvars
    use mod_global_parameters

    integer                                   :: igrid, iigrid, idir
    double precision, dimension(ixG^T,1:nw)   :: w

    do iigrid=1,igridstail; igrid=igrids(iigrid);

      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

      gridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0
      w(ixG^T,1:nw) = ps(igrid)%w(ixG^T,1:nw)
      call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
      ! fill with velocity:
      gridvars(igrid)%w(ixG^T,vp(:)) = w(ixG^T,iw_mom(:))

      if(time_advance) then
        gridvars(igrid)%wold(ixG^T,1:ngridvars) = 0.0d0
        w(ixG^T,1:nw) = pso(igrid)%w(ixG^T,1:nw)
        call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
        gridvars(igrid)%wold(ixG^T,vp(:)) = w(ixG^T,iw_mom(:))
      end if

    end do

  end subroutine advect_fill_gridvars

  subroutine advect_integrate_particles(end_time)
    ! this solves dx/dt=v for particles
    use mod_odeint
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
      dt_p                    = advect_get_particle_dt(particle(ipart), end_time)
      particle(ipart)%self%dt = dt_p

      igrid                   = particle(ipart)%igrid
      igrid_working           = igrid
      tloc                    = particle(ipart)%self%t
      x(1:ndir)               = particle(ipart)%self%x(1:ndir)
      tlocnew                 = tloc+dt_p

      ! Position update
      ! Simple forward Euler start
      !call get_vec_advect(igrid,x,tloc,v,vp(1),vp(ndir))
      !particle(ipart)%self%u(1:ndir) = v(1:ndir)
      !particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
      !     + dt_p * v(1:ndir)
      ! Simple forward Euler end

      ! Adaptive stepwidth RK4:
      h1 = dt_p/2.0d0
      call odeint(x,ndir,tloc,tlocnew,eps,h1,hmin,nok,nbad,derivs_advect,rkqs,ierror)

      if (ierror /= 0) then
        print *, "odeint returned error code", ierror
        print *, "1 means hmin too small, 2 means MAXSTP exceeded"
        print *, "Having a problem with particle", iipart
      end if

      particle(ipart)%self%x(1:ndir) = x(1:ndir)

      ! Velocity update
      call get_vec_advect(igrid,x,tlocnew,v,vp(1),vp(ndir))
      particle(ipart)%self%u(1:ndir) = v(1:ndir)

      ! Time update
      particle(ipart)%self%t = tlocnew

      ! Update payload
      if (.not. associated(usr_update_payload)) then
        call advect_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x,payload,npayload,tlocnew)
      else
        call usr_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x,payload,npayload,tlocnew)
      end if
      particle(ipart)%payload = payload

    end do

  end subroutine advect_integrate_particles

  !> Example of update payload with local density
  subroutine advect_update_payload(igrid,w,wold,xgrid,xpart,payload,npayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,npayload
    double precision, intent(in)  :: w(ixG^T,1:nw),wold(ixG^T,1:nw)
    double precision, intent(in)  :: xgrid(ixG^T,1:ndim),xpart(1:ndir),particle_time
    double precision, intent(out) :: payload(npayload)
    double precision              :: rho, rho1, rho2, td

    if (.not.time_advance) then
      call interpolate_var(igrid,ixG^LL,ixM^LL,w(ixG^T,iw_rho),xgrid,xpart,rho)
    else
      call interpolate_var(igrid,ixG^LL,ixM^LL,wold(ixG^T,iw_rho),xgrid,xpart,rho1)
      call interpolate_var(igrid,ixG^LL,ixM^LL,w(ixG^T,iw_rho),xgrid,xpart,rho2)
      td = (particle_time - global_time) / dt
      rho = rho1 * (1.0d0 - td) + rho2 * td
    end if
    payload(1) = rho * w_convert_factor(1)

  end subroutine advect_update_payload

  subroutine derivs_advect(t_s,x,dxdt)
    use mod_global_parameters
    double precision :: t_s, x(ndir)
    double precision :: dxdt(ndir)

    double precision :: v(ndir)

    call get_vec_advect(igrid_working,x,t_s,v,vp(1),vp(ndir))
    dxdt(:) = v(:)

  end subroutine derivs_advect

  pure function advect_get_particle_dt(partp, end_time) result(dt_p)
    use mod_global_parameters
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

    ! make sure we step only one cell at a time:
    v(1:ndir)=abs(partp%self%u(1:ndir))

    ! convert to angular velocity:
    if(typeaxial =='cylindrical'.and.phi_>0) v(phi_) = abs(v(phi_)/partp%self%x(r_))

    dt_cfl = min({rnode(rpdx^D_,partp%igrid)/v(^D)},bigdouble)

    if(typeaxial =='cylindrical'.and.phi_>0) then
      ! phi-momentum leads to radial velocity:
      if(phi_ .gt. ndim) dt_cfl = min(dt_cfl, &
           sqrt(rnode(rpdx1_,partp%igrid)/partp%self%x(r_)) &
           / v(phi_))
      ! limit the delta phi of the orbit (just for aesthetic reasons):
      dt_cfl = min(dt_cfl,0.1d0/v(phi_))
      ! take some care at the axis:
      dt_cfl = min(dt_cfl,(partp%self%x(r_)+smalldouble)/v(r_))
    end if

    dt_p = dt_cfl

    ! Make sure we don't advance beyond end_time
    call limit_dt_endtime(end_time - partp%self%t, dt_p)

  end function advect_get_particle_dt

  subroutine get_vec_advect(igrid,x,tloc,var,ibeg,iend)
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

  end subroutine get_vec_advect

end module mod_particle_advect
