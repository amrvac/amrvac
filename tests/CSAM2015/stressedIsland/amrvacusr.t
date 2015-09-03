!=============================================================================
! amrvacusr.t.nul
!=============================================================================
!INCLUDE:amrvacnul/speciallog.t
!INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/correctaux_usr.t
INCLUDE:amrvacmodules/handle_particles.t
INCLUDE:amrvacmodules/integrate_particles.t
!=============================================================================
subroutine initglobaldata_usr

use constants
include 'amrvacdef.f'
double precision scaleb
!-----------------------------------------------------------------------------

eqpar(kpar_)   = 2000.0d0
eqpar(kperp_)  = 0.0d0
eqpar(kappa_)  = eqpar(kpar_)/2.0d0

eqpar(vcoll_) = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Scaling of magnetic field:
scaleb = 5.0d-3 ! Gauss
UNIT_LENGTH   = 1.0d10 ! cm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
UNIT_DENSITY  = (scaleb)**2/(4.0d0*dpi*CONST_c**2)
UNIT_VELOCITY = CONST_c ! cm/c

normvar(0) = UNIT_LENGTH
{^C&normvar(b^C_)   = sqrt(4.0d0*dpi*UNIT_DENSITY)*UNIT_VELOCITY \}
{^C&normvar(e^C_)   = sqrt(4.0d0*dpi*UNIT_DENSITY)*UNIT_VELOCITY \}
normvar(q_) = sqrt(4.0d0*dpi*UNIT_DENSITY)*UNIT_VELOCITY / UNIT_LENGTH
normt = UNIT_LENGTH/UNIT_VELOCITY

end subroutine initglobaldata_usr
!=============================================================================
subroutine init_particle_integrator()

use mod_particles
use constants
include 'amrvacdef.f'
!-----------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set the number of steps:
itmax_particles = 10000000
tmax_particles  = tmax * UNIT_LENGTH/UNIT_VELOCITY
ditsave_particles = 4
dtsave_ensemble   = dtsave(2) * UNIT_LENGTH/UNIT_VELOCITY
dtheta            = 2.0d0 * dpi / 60.0d0

losses = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine init_particle_integrator
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid within ix^L

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
! .. local ..
double precision, dimension(ixG^T)       :: jdotb, b2
double precision,dimension(ixG^T,1:ndir) :: j
double precision                         :: Lx, Ly
integer                                  :: idirmin
integer, parameter                       :: idirmin0=1
!-----------------------------------------------------------------------------
w(ixG^S,1:nw) = zero

w(ixG^S,phib_) = 0.0d0
w(ixG^S,psi_)  = 0.0d0

Lx = xprobmax1 - xprobmin1
Ly = xprobmax2 - xprobmin2

w(ixG^S,b1_) = - dsin(two*dpi*x(ixG^S,2)/Ly)
w(ixG^S,b2_) = + dsin(two*dpi*(x(ixG^S,1)/Lx-half))
w(ixG^S,b3_) =  ( dcos(two*dpi*(x(ixG^S,1)/Lx-half)) + dcos(two*dpi*x(ixG^S,2)/Ly) ) 

! Calculate Epar:
call curlvector(w(ixG^T,b1_:b3_),ixG^LL,ix^L^LADD1,j,idirmin,idirmin0,ndir)

jdotb(ixG^S) = {^C& j(ixG^S,^C)*w(ixG^S,b^C_) |+}
b2(ixG^S)    = {^C& w(ixG^S,b^C_)**2 |+}

{^C& w(ixG^S,e^C_) = one/eqpar(kpar_) * jdotb(ixG^S)/b2(ixG^S) * w(ixG^S,b^C_)\}

! Add the perturbation on E
where(x(ixG^S,1) .lt. -0.5d0)
   w(ixG^S,e2_)  = w(ixG^S,e2_) + w(ixG^S,b3_)*eqpar(vcoll_)
   w(ixG^S,e3_)  = w(ixG^S,e3_) - w(ixG^S,b2_)*eqpar(vcoll_)
elsewhere
   w(ixG^S,e2_)  = w(ixG^S,e2_) - w(ixG^S,b3_)*eqpar(vcoll_)
   w(ixG^S,e3_)  = w(ixG^S,e3_) + w(ixG^S,b2_)*eqpar(vcoll_)
end where


! Calculate corresponding charge density:
call divvector(w(ixG^T,e1_:e3_),ixG^LL,ix^L,w(ixG^T,q_))

end subroutine initonegrid_usr
!=============================================================================
subroutine init_particles()
! initialise the particles

use constants
use mod_particles
use Knuth_random
include 'amrvacdef.f'

double precision, dimension(ndir)    :: x
integer                              :: igrid_particle, ipe_particle
integer, parameter                   :: Npart=100000
double precision, parameter          :: gamma0 = 10.0d0
integer                              :: seed
double precision                     :: r^C(1:Npart), s^C(1:Npart)
integer                              :: ipart
logical, dimension(1:Npart)          :: follow=.false.
!-----------------------------------------------------------------------------


! initialise the random number generator
seed = 310952
call rnstrt(seed)
{^C&
call rand(r^C,Npart)
call rand(s^C,Npart)
\}


! flags to follow particles
follow(1:10) = .true.

follow(8400) = .true.
follow(70098) = .true.
follow(9082) = .true.
follow(36948) = .true.
follow(46794) = .true.
follow(89461) = .true.
follow(37702) = .true.
follow(64511) = .true.
follow(33305) = .true.
follow(44695) = .true.


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

      particle(nparticles)%self%u(1) = gamma0 * s1(nparticles)
      particle(nparticles)%self%u(2) = gamma0 * s2(nparticles)
      particle(nparticles)%self%u(3) = gamma0 * s3(nparticles)
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
subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
! corresponding normalization values (default value 1)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)
! .. local ..
double precision,dimension(ixG^T,1:ndir) :: current
double precision,dimension(ixG^T,1:ndir) :: curlb
double precision,dimension(ixG^T)        :: vidir, divE
double precision                         :: evec(ixG^T,1:ndir)
integer                                  :: idirmin
integer, parameter                       :: idirmin0=1
!-----------------------------------------------------------------------------

call getcurrent(ixI^L,ixO^L,w,x,saveprim,current)
w(ixO^S,nw+1) = current(ixO^S,3)
normconv(nw+1) = normvar(b1_)/UNIT_LENGTH

! Also get curl of b:
call curlvector(w(ixG^T,b1_:b3_),ixI^L,ixO^L,curlb,idirmin,idirmin0,ndir)
w(ixO^S,nw+2) = curlb(ixO^S,3)
normconv(nw+2) = normvar(b1_)/UNIT_LENGTH

! Get drift velocity:
{^C&
call getv(w,x,ixI^L,ixO^L,^C,vidir)
w(ixO^S,nw+2+^C) = vidir(ixO^S)
normconv(nw+2+^C) = UNIT_VELOCITY
\}

! S[Ei_] = - divE (E x B)/B^2
evec(ixI^S,1:ndir)=w(ixI^S,e0_+1:e0_+ndir)
call divvector(evec,ixI^L,ixO^L,divE)
w(ixO^S,nw+2+ndir+1) = divE(ixO^S)
normconv(nw+2+ndir+1) = normvar(e1_)/UNIT_LENGTH

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'jz'//' '//'curlbz'//' '//'v1'//' '//'v2'//' '//'v3'//' '//'divE'
wnames=TRIM(wnames)//' '//'jz'//' '//'curlbz'//' '//'v1'//' '//'v2'//' '//'v3'//' '//'divE'

end subroutine specialvarnames_output
!=============================================================================

















!=============================================================================
!================= JUST DUMMIES ==============================================
!=============================================================================
subroutine flag_grid_usr(qt,ixG^L,ixO^L,w,x,flag)

include 'amrvacdef.f'

integer, intent(in)             :: ixG^L, ixO^L
integer, intent(inout)          :: flag
double precision, intent(in)    :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in)    :: x(ixG^S,1:ndim)

! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
! flag=0  : Treat as normal domain
! flag=1  : Treat as passive, but reduce by safety belt
! flag=2  : Always treat as passive

!-----------------------------------------------------------------------------
      
end subroutine flag_grid_usr
!=============================================================================
{#IFDEF FCT
subroutine initvecpot_usr(ixI^L, ixC^L, xC, A)

! initialize the vectorpotential on the corners
! used by b_from_vectorpotential()


include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixC^L
double precision, intent(in)       :: xC(ixI^S,1:ndim)
double precision, intent(out)      :: A(ixI^S,1:ndir)

!double precision                   :: r(ixG^T)

!-----------------------------------------------------------------------------

!r(ixC^S)=sqrt(xC(ixC^S,1)**2 + xC(ixC^S,2)**2 )

A(ixC^S,1:ndir) = zero
!where (r(ixC^S) .lt. eqpar(rm_))
!   A(ixC^S,3) = - half * eqpar(bm_) / eqpar(rm_) * r(ixC^S)**2
!elsewhere (r(ixC^S) .lt. eqpar(rj_))
!   A(ixC^S,3) = - half * eqpar(bm_) * eqpar(rm_) &
!        - eqpar(bm_) * eqpar(rm_) * log(r(ixC^S)/eqpar(rm_))
!elsewhere
!   A(ixC^S,3) = - half * eqpar(bm_) * eqpar(rm_) &
!        - eqpar(bm_) * eqpar(rm_) * log(eqpar(rj_)/eqpar(rm_))
!end where


end subroutine initvecpot_usr
!=============================================================================
}



!=============================================================================
subroutine printlog_special

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

include 'amrvacdef.f'

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
{#IFDEF PROCESSGLOBAL
!==============================================================================

subroutine process_global_usr(iit,qt)
!
! This subroutine is called at the beginning of each time step 
! by each processor. No communication is specified, so the user
! has to implement MPI routines if information has to be shared
!

include 'amrvacdef.f'

integer, intent(in)          :: iit
double precision, intent(in) :: qt

!-----------------------------------------------


return
end subroutine process_global_usr
}
!=============================================================================
{#IFDEF UCONVERT
subroutine userspecialconvert(qunitconvert)
! Allow user to use their own data-postprocessing procedures

include 'amrvacdef.f'
integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------

end subroutine userspecialconvert
!=============================================================================
}
{#IFDEF TRANSFORMW
subroutine transformw_usr(w,wtf,eqpar_tf,ixI^L,ixO^L)
! regenerate w and eqpar arrays to output into *tf.dat, e.g., add/remove e_
! variable
include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: wtf(ixI^S,1:nwtf)
double precision, intent(out):: eqpar_tf(neqpartf)

!-----------------------------------------------------------------------------

end subroutine transformw_usr
!=============================================================================
}
{#IFDEF SPECIALTOLERANCE
subroutine special_tolerance(xlocal,tolerance)
!PURPOSE: use different tolerance in special regions for AMR to
!reduce/increase resolution there where nothing/something interesting happens.
include 'amrvacdef.f'

double precision, intent(in) :: xlocal(1:ndim)
double precision, intent(inout) :: tolerance

double precision :: bczone^D,addtol,tol_add
!-----------------------------------------------------------------------------
!amplitude of additional tolerance
addtol=0.3d0
! thickness of near-boundary region
bczone1=0.2d0*(xprobmax1-xprobmin1)
! linear changing of additional tolerance
if(xlocal(1)-xprobmin1 < bczone1 .or. xprobmax1-xlocal(1) < bczone1) then
  tol_add=(1.d0-min(xlocal(1)-xprobmin1,xprobmax1-xlocal(1))/bczone1)*addtol
endif
bczone2=0.2d0*(xprobmax2-xprobmin2)
if(xlocal(2)-xprobmin2 < bczone2 .or. xprobmax2-xlocal(2) < bczone2) then
  tol_add=(1.d0-min(xlocal(2)-xprobmin2,xprobmax2-xlocal(2))/bczone2)*addtol
endif
bczone3=0.2d0*(xprobmax3-xprobmin3)
if(xprobmax3-xlocal(3) < bczone3) then
  tol_add=(1.d0-(xprobmax3-xlocal(3))/bczone3)*addtol
endif
tolerance=tolerance+tol_add

end subroutine special_tolerance
!=============================================================================
}

!=============================================================================
! amrvacusr.t.nul
!=============================================================================
