!=============================================================================
subroutine init_particle_integrator()

use mod_particles
use constants
include 'amrvacdef.f'
!-----------------------------------------------------------------------------

itmax_particles = 10000000
tmax_particles  = 10.0d0 * CONST_years
ditsave_particles = 1
dtsave_ensemble   = dtsave(2) * UNIT_LENGTH/UNIT_VELOCITY
dtheta            = 2.0d0 * dpi / 60.0d0

losses = .false.

end subroutine init_particle_integrator
!=============================================================================
subroutine init_particles()
! initialise the particles

use constants
use mod_particles
use mod_gridvars
use Knuth_random
include 'amrvacdef.f'

double precision, dimension(ndir)    :: x
integer                              :: igrid_particle, ipe_particle
integer, parameter                   :: Npart=100
integer                              :: seed
double precision                     :: r^C(1:Npart), s^C(1:Npart)
double precision                     :: v(1:ndir)
double precision                     :: lfac
integer                              :: ipart, idir
logical, dimension(1:Npart)          :: follow=.false.
{#IFDEF PARTICLES_GCA
double precision,dimension(1:ndir)   :: B
double precision                     :: absB, gamma
}
!-----------------------------------------------------------------------------


! initialise the random number generator
seed = 310952
call rnstrt(seed)
{^C&
call rand(r^C,Npart)
call rand(s^C,Npart)
\}

!{^C&r^C(1) = 1.0d0 \}

! flags to follow particles
follow(31)   = .true.
follow(81)  = .true.
follow(49)  = .true.
!follow(81)  = .true.
!follow(100) = .true.
!follow(20) = .true.
!follow(23) = .true.
!follow(5) = .true.
!follow(855) = .true.


! first find ipe and igrid responsible for particle

x(:)=0.0d0

do while (nparticles .lt. Npart)

{^D&x(^D) = xprobmin^D + r^D(nparticles+1) * (xprobmax^D - xprobmin^D)\}

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

{#IFDEF PARTICLES_GCA
      lfac=1.0d6
      gamma = sqrt(1.0d0 + lfac**2 * ({^C& s^C(nparticles)**2 |+}) ) ! particles Lorentz factor
      call get_vec(igrid_particle,x,particle(nparticles)%self%t,B,bp1_,bp^NC_)
      absB = sqrt({^C& B(^C)**2 |+})

      particle(nparticles)%self%u(1) = lfac * s1(nparticles) * CONST_c ! parallel momentum component (gamma v||)

      particle(nparticles)%self%u(2) = (lfac * sqrt(s2(nparticles)**2+s3(nparticles)**2) * CONST_c * particle(nparticles)%self%m )**2 / (2.0d0 * particle(nparticles)%self%m * absB) ! Mr: the conserved magnetic moment

      particle(nparticles)%self%u(3) = gamma ! Lorentz factor
}

{#IFDEF PARTICLES_LORENTZ
      particle(nparticles)%self%u(1) = 1.d7 * s1(nparticles)
      particle(nparticles)%self%u(2) = 1.d7 * s2(nparticles)
      particle(nparticles)%self%u(3) = 1.d7 * s3(nparticles)
}
{#IFDEF PARTICLES_ADVECT
      {^C&call interpolate_var(igrid_particle,ixG^LL,ixM^LL,pw(igrid_particle)%w(ixG^T,v^C_),px(igrid_particle)%x(ixG^T,1:ndim),x,v(^C))\}

   if (useprimitiveRel) then
      call interpolate_var(igrid_particle,ixG^LL,ixM^LL,pw(igrid_particle)%w(ixG^T,lfac_),px(igrid_particle)%x(ixG^T,1:ndim),x,lfac)
      {^C&particle(nparticles)%self%u(^C) = v(^C)/lfac * UNIT_VELOCITY\}
   else
      {^C&particle(nparticles)%self%u(^C) = v(^C) * UNIT_VELOCITY\}
   end if
}
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

