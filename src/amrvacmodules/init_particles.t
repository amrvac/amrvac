!=============================================================================
subroutine init_particles()
! initialise the particles

use constants
use mod_particles
use Knuth_random
include 'amrvacdef.f'

double precision, dimension(ndir)    :: x
integer                              :: igrid_particle, ipe_particle
integer, parameter                   :: Npart=100
integer                              :: seed
double precision                     :: r^C(1:Npart)
integer                              :: ipart
logical, dimension(1:Npart)          :: follow=.false.
!-----------------------------------------------------------------------------

! set the number of steps:

itmax_particles = 10000000
tmax_particles  = 100.0d0 * CONST_years
ditsave_particles = 1
dtsave_ensemble   = 10.0d0 * CONST_years

losses = .false.

! initialise the random number generator
seed = 310952
call rnstrt(seed)
{^C&
call rand(r^C,Npart)\}

!{^C&r^C(1) = 1.0d0 \}

! flags to follow particles
follow(1)   = .true.
!follow(81)  = .true.
!follow(100) = .true.
!follow(20) = .true.
!follow(23) = .true.
!follow(5) = .true.
!follow(855) = .true.


! first find ipe and igrid responsible for particle

x(1) = dpi
x(2) = dpi
x(3) = 0.0d0

do while (nparticles .lt. Npart)

   call find_particle_ipe(x,igrid_particle,ipe_particle)

   nparticles=nparticles+1
   particle(nparticles)%igrid  = igrid_particle
   particle(nparticles)%ipe    = ipe_particle

   if (ipe_particle == mype) then 

      call push_particle_into_particles_on_mype(nparticles)
      allocate(particle(nparticles)%self)
      particle(nparticles)%self%follow = follow(nparticles)
      particle(nparticles)%self%index  = nparticles
      particle(nparticles)%self%q      = - CONST_e
      particle(nparticles)%self%m      =   CONST_me

      particle(nparticles)%self%t      = 0.0d0
      particle(nparticles)%self%dt     = 0.0d0

      particle(nparticles)%self%x(:) = x
      particle(nparticles)%self%u(1) = 1.d7 * r1(nparticles)
      particle(nparticles)%self%u(2) = 1.d7 * r2(nparticles)
      particle(nparticles)%self%u(3) = 1.d7 * r3(nparticles)

      particle(nparticles)%self%payload(1:npayload) = 0.0d0

   end if

end do

end subroutine init_particles
!=============================================================================
logical function user_destroy(myparticle)

use mod_particles, only: particle_node

type(particle_node), intent(in)          :: myparticle
!-----------------------------------------------------------------------------
user_destroy = .false.

end function user_destroy
!=============================================================================
subroutine init_gridvars()

use mod_gridvars
include 'amrvacdef.f'
integer                                   :: igrid, iigrid, idir
double precision, dimension(ixG^T,1:ndir) :: beta
double precision, dimension(ixG^T,1:nw)   :: w
!-----------------------------------------------------------------------------

do igrid=1,ngridshi
   nullify(gridvars(igrid)%w)
   if (time_advance) nullify(gridvars_old(igrid)%w)
end do

do iigrid=1,igridstail; igrid=igrids(iigrid);
   allocate(gridvars(igrid)%w(ixG^T,1:ngridvars))
   if (time_advance)  allocate(gridvars_old(igrid)%w(ixG^T,1:ngridvars))
end do


! now fill the gridvars:

do iigrid=1,igridstail; igrid=igrids(iigrid);

   w(ixG^T,1:nw) = pw(igrid)%w(ixG^T,1:nw)
   call primitive(ixG^LL,ixG^LL,w,px(igrid)%x)

   ! fill with magnetic field:
   gridvars(igrid)%w(ixG^T,bp1_:bp^NC_) = w(ixG^T,b1_:b^NC_) &
        * sqrt(4.0d0*dpi*UNIT_VELOCITY**2.0d0 * UNIT_DENSITY)

   ! get the electric field
   beta(ixG^T,1:ndir) = w(ixG^T,u1_:u^NC_)
   if (useprimitiveRel) then
      do idir=1,ndir
         beta(ixG^T,idir) = beta(ixG^T,idir)/w(ixG^T,lfac_) 
      end do
   end if

   select case (typeaxial)
   case ('slab')
      gridvars(igrid)%w(ixG^T,ep1_) = gridvars(igrid)%w(ixG^T,2) * beta(ixG^T,3) &
           - gridvars(igrid)%w(ixG^T,3) * beta(ixG^T,2)
      
      gridvars(igrid)%w(ixG^T,ep2_) = gridvars(igrid)%w(ixG^T,3) * beta(ixG^T,1) &
           - gridvars(igrid)%w(ixG^T,1) * beta(ixG^T,3)
      
      gridvars(igrid)%w(ixG^T,ep3_) = gridvars(igrid)%w(ixG^T,1) * beta(ixG^T,2) &
           - gridvars(igrid)%w(ixG^T,2) * beta(ixG^T,1)
   case ('cylindrical')
      gridvars(igrid)%w(ixG^T,ep1_) = gridvars(igrid)%w(ixG^T,^PHI) * beta(ixG^T,^Z) &
           - gridvars(igrid)%w(ixG^T,^Z) * beta(ixG^T,^PHI)
      
      gridvars(igrid)%w(ixG^T,ep^PHI_) = gridvars(igrid)%w(ixG^T,^Z) * beta(ixG^T,1) &
           - gridvars(igrid)%w(ixG^T,1) * beta(ixG^T,^Z)
      
      gridvars(igrid)%w(ixG^T,ep^Z_) = gridvars(igrid)%w(ixG^T,1) * beta(ixG^T,^PHI) &
           - gridvars(igrid)%w(ixG^T,^PHI) * beta(ixG^T,1)
      
   end select


   ! NOW DO THE SAME FOR OLD VARS (sorry for the spaghetti code)

   if (time_advance) then 
      
      w(ixG^T,1:nw) = pwold(igrid)%w(ixG^T,1:nw)
      call primitive(ixG^LL,ixG^LL,w,px(igrid)%x)
      
      ! fill with magnetic field:
      gridvars_old(igrid)%w(ixG^T,bp1_:bp^NC_) = w(ixG^T,b1_:b^NC_) &
           * sqrt(4.0d0*dpi*UNIT_VELOCITY**2.0d0 * UNIT_DENSITY)
      
      ! get the electric field
      beta(ixG^T,1:ndir) = w(ixG^T,u1_:u^NC_)
      if (useprimitiveRel) then
         do idir=1,ndir
            beta(ixG^T,idir) = beta(ixG^T,idir)/w(ixG^T,lfac_) 
         end do
      end if
      
      select case (typeaxial)
      case ('slab')
         gridvars_old(igrid)%w(ixG^T,ep1_) = gridvars_old(igrid)%w(ixG^T,2) * beta(ixG^T,3) &
              - gridvars_old(igrid)%w(ixG^T,3) * beta(ixG^T,2)
         
         gridvars_old(igrid)%w(ixG^T,ep2_) = gridvars_old(igrid)%w(ixG^T,3) * beta(ixG^T,1) &
              - gridvars_old(igrid)%w(ixG^T,1) * beta(ixG^T,3)
         
         gridvars_old(igrid)%w(ixG^T,ep3_) = gridvars_old(igrid)%w(ixG^T,1) * beta(ixG^T,2) &
              - gridvars_old(igrid)%w(ixG^T,2) * beta(ixG^T,1)
      case ('cylindrical')
         gridvars_old(igrid)%w(ixG^T,ep1_) = gridvars_old(igrid)%w(ixG^T,^PHI) * beta(ixG^T,^Z) &
              - gridvars_old(igrid)%w(ixG^T,^Z) * beta(ixG^T,^PHI)
         
         gridvars_old(igrid)%w(ixG^T,ep^PHI_) = gridvars_old(igrid)%w(ixG^T,^Z) * beta(ixG^T,1) &
              - gridvars_old(igrid)%w(ixG^T,1) * beta(ixG^T,^Z)
         
         gridvars_old(igrid)%w(ixG^T,ep^Z_) = gridvars_old(igrid)%w(ixG^T,1) * beta(ixG^T,^PHI) &
              - gridvars_old(igrid)%w(ixG^T,^PHI) * beta(ixG^T,1)
         
      end select
      
   end if ! time_advance

end do

end subroutine init_gridvars
!=============================================================================
