!> Scattered sampling based on fixed- or moving-particle interpolation
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
    use mod_usr_methods, only: usr_create_particles, usr_update_payload, &
                               usr_check_particle, usr_particle_position

    double precision :: x(3, num_particles)
    double precision :: v(3, num_particles)
    double precision :: q(num_particles)
    double precision :: m(num_particles)
    double precision :: rrd(num_particles,ndir)
    double precision :: w(ixG^T,1:nw)
    double precision :: defpayload(ndefpayload)
    double precision :: usrpayload(nusrpayload) 
    integer          :: n, idir, igrid, ipe_particle, nparticles_local
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
        call sample_update_payload(igrid,x(:,n),v(:,n),q(n),m(n),defpayload,ndefpayload,0.d0)
        particle(n)%payload(1:ndefpayload) = defpayload      
        if (associated(usr_update_payload)) then
          call usr_update_payload(igrid,x(:,n),v(:,n),q(n),m(n),usrpayload,nusrpayload,0.d0)
          particle(n)%payload(ndefpayload+1:npayload)=usrpayload
        end if
      end if

    end do

    call MPI_ALLREDUCE(nparticles_local,nparticles,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)

    ! write the first csv file of particles
    t_next_output=global_time
    call write_particle_output()
    t_next_output=t_next_output+dtsave_particles

  end subroutine sample_create_particles

  subroutine sample_fill_gridvars
    use mod_global_parameters

    double precision, dimension(ixG^T,1:nw)   :: w
    integer                                   :: igrid, iigrid, idir

    do iigrid=1,igridstail; igrid=igrids(iigrid);

      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      w(ixG^T,1:nw) = ps(igrid)%w(ixG^T,1:nw)
      call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
      ! fill all variables:
      gridvars(igrid)%w(ixG^T,1:ngridvars) = w(ixG^T,1:ngridvars)

    end do

  end subroutine sample_fill_gridvars

  subroutine sample_integrate_particles(end_time)
    ! this interpolates the HD/MHD quantities at the particle positions
    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles, usr_update_payload, usr_particle_position
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

      ! Position update (if defined)
      ! TODO: this may create problems with interpolation out of boundaries
      if (associated(usr_particle_position)) call usr_particle_position(x,particle(ipart)%self%index,tloc,tlocnew)
      particle(ipart)%self%x(1:ndir) = x

      ! Time update
      particle(ipart)%self%time = tlocnew

      ! Update payload
      call sample_update_payload(igrid,x,v,0.d0,0.d0,defpayload,ndefpayload,tlocnew)
      particle(ipart)%payload(1:ndefpayload) = defpayload
      if (associated(usr_update_payload)) then
        call usr_update_payload(igrid,x,v,0.d0,0.d0,usrpayload,nusrpayload,tlocnew)
        particle(ipart)%payload(ndefpayload+1:npayload) = usrpayload
      end if

    end do

  end subroutine sample_integrate_particles

  !> Payload update
  subroutine sample_update_payload(igrid,xpart,upart,qpart,mpart,mypayload,mynpayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,mynpayload
    double precision, intent(in)  :: xpart(1:ndir),upart(1:ndir),qpart,mpart,particle_time
    double precision, intent(out) :: mypayload(mynpayload)
    double precision              :: myw(ixG^T,1:nw),mywold(ixG^T,1:nw)
    double precision              :: wp, wpold, td
    integer                       :: ii


    ! There are npayload=nw payloads, one for each primitive fluid quantity
    myw(ixG^T,1:nw) = gridvars(igrid)%w(ixG^T,1:nw)
    if (time_advance) mywold(ixG^T,1:nw) = gridvars(igrid)%wold(ixG^T,1:nw)

    if (.not.saveprim) then
      call phys_to_conserved(ixG^LL,ixG^LL,myw,ps(igrid)%x)
      if (time_advance) call phys_to_conserved(ixG^LL,ixG^LL,mywold,ps(igrid)%x)
    end if

    do ii=1,mynpayload
      call interpolate_var(igrid,ixG^LL,ixM^LL,myw(ixG^T,ii),ps(igrid)%x,xpart,wp)
      if (time_advance) then
        td = (particle_time - global_time) / dt
        call interpolate_var(igrid,ixG^LL,ixM^LL,mywold(ixG^T,ii),ps(igrid)%x,xpart,wpold)
        wp = wpold * (1.0d0 - td) + wp * td
      end if
      mypayload(ii) = wp*w_convert_factor(ii)
    end do

  end subroutine sample_update_payload

  function sample_get_particle_dt(partp, end_time) result(dt_p)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_particle_position
    type(particle_ptr), intent(in) :: partp
    double precision, intent(in)   :: end_time
    double precision               :: dt_p
    double precision               :: tout, dt_cfl
    double precision               :: v(1:ndir), xp(3), told, tnew
    integer                        :: ipart, iipart, nout, id

    if (const_dt_particles > 0) then
      dt_p = const_dt_particles
      return
    end if

    dt_p = dtsave_particles

    ! Make sure the user-defined particle movement doesn't break communication
    if (associated(usr_particle_position)) then
      xp = partp%self%x
      told = partp%self%time
      tnew = told+dt_p
      call usr_particle_position(xp, partp%self%index, told, tnew)
      do while (.not. point_in_igrid_ghostc(xp,partp%igrid,1))
        dt_p = dt_p/10.d0
        xp = partp%self%x
        tnew = told+dt_p
        call usr_particle_position(xp, partp%self%index, told, tnew)
      end do
    end if

    ! Make sure we don't advance beyond end_time
    call limit_dt_endtime(end_time - partp%self%time, dt_p)

  end function sample_get_particle_dt

end module mod_particle_sample
