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
use mod_global_parameters

double precision scaleb
!-----------------------------------------------------------------------------

eqpar(kpar_)   = 1000.0d0
eqpar(kperp_)  = 0.0d0
eqpar(kappa_)  = eqpar(kpar_)/2.0d0

eqpar(C_)     = 0.01d0
eqpar(vcoll_) = 0.1d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Scaling of magnetic field:
scaleb = 5.0d-3 ! Gauss
UNIT_LENGTH   = 1.0d10 ! cm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
UNIT_DENSITY  = (scaleb)**2/(4.0d0*dpi*CONST_c**2)
UNIT_VELOCITY = CONST_c ! cm/s


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
use mod_global_parameters
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

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
! .. local ..
double precision, dimension(ixG^T)  :: bphi1, bphi2, bz1, bz2, r1, r2, sinphi1, sinphi2, cosphi1, cosphi2, jdotb, b2
double precision,dimension(ixG^T,1:ndir) :: j
integer                             :: ix^D
integer                                  :: idirmin
integer, parameter                       :: idirmin0=1
logical patchw(ixG^T)
double precision, parameter     :: alphat = 3.8317d0
double precision                :: BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1
!
patchw(ixG^S) = .false.
!-----------------------------------------------------------------------------
! Do some geometry:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
r1(ixG^S) = sqrt((x(ixG^S,1)+one)**2+x(ixG^S,2)**2)
r2(ixG^S) = sqrt((x(ixG^S,1)-one)**2+x(ixG^S,2)**2)

sinphi1(ixG^S) = x(ixG^S,2)/r1(ixG^S)
cosphi1(ixG^S) = (x(ixG^S,1)+one)/r1(ixG^S)

sinphi2(ixG^S) = x(ixG^S,2)/r2(ixG^S)
cosphi2(ixG^S) = (x(ixG^S,1)-one)/r2(ixG^S)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

w(ixG^S,1:nw) = zero

w(ixG^S,phib_) = 0.0d0
w(ixG^S,psi_)  = 0.0d0

! Toroidal magnetic field:
{^D& do ix^D=ixGmin^D,ixGmax^D\}
    call JY01A(alphat*r1(ix^D),BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
    bphi1(ix^D) = BJ1
    call JY01A(alphat*r2(ix^D),BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
    bphi2(ix^D) = BJ1
{^D& end do\}

! Poloidal (z) magnetic field:
{^D& do ix^D=ixGmin^D,ixGmax^D\}
    call JY01A(alphat*r1(ix^D),BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
    bz1(ix^D) = sqrt(BJ0**2 + eqpar(C_))
    call JY01A(alphat*r2(ix^D),BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
    bz2(ix^D) = sqrt(BJ0**2 + eqpar(C_))
{^D& end do\}

! Fill magnetic field:
call JY01A(alphat,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
where(r1(ixG^S) .lt. one)
   w(ixG^S,b1_) = - sinphi1(ixG^S) * bphi1(ixG^S)
   w(ixG^S,b2_) = + cosphi1(ixG^S) * bphi1(ixG^S)
   w(ixG^S,b3_) = bz1(ixG^S)
elsewhere(r2(ixG^S) .lt. one)
   w(ixG^S,b1_) = - sinphi2(ixG^S) * bphi2(ixG^S)
   w(ixG^S,b2_) = + cosphi2(ixG^S) * bphi2(ixG^S)
   w(ixG^S,b3_) = bz2(ixG^S)
elsewhere
   w(ixG^S,b1_) = zero
   w(ixG^S,b2_) = zero
   w(ixG^S,b3_) = sqrt(BJ0**2 + eqpar(C_))
end where

! Calculate Epar:
call curlvector(w(ixG^T,b1_:b3_),ixG^LL,ix^L^LADD1,j,idirmin,idirmin0,ndir)

jdotb(ixG^S) = {^C& j(ixG^S,^C)*w(ixG^T,b^C_) |+}
b2(ixG^S)    = {^C& w(ixG^T,b^C_)**2 |+}

{^C& w(ixG^S,e^C_) = one/eqpar(kpar_) * jdotb(ixG^S)/b2(ixG^S) * w(ixG^S,b^C_)\}


! Add the perturbation on E
where(r1(ixG^S) .lt. one)
   w(ixG^S,e2_)  = w(ixG^S,e2_) + w(ixG^S,b3_)*eqpar(vcoll_)
   w(ixG^S,e3_)  = w(ixG^S,e3_) - w(ixG^S,b2_)*eqpar(vcoll_)
elsewhere(r2(ixG^S) .lt. one)
   w(ixG^S,e2_)  = w(ixG^S,e2_) - w(ixG^S,b3_)*eqpar(vcoll_)
   w(ixG^S,e3_)  = w(ixG^S,e3_) + w(ixG^S,b2_)*eqpar(vcoll_)
elsewhere

end where


! Calculate the charge density:
call divvector(w(ixG^T,e1_:e3_),ixG^L,ix^L,w(ixG^T,q_))

end subroutine initonegrid_usr
!=============================================================================
subroutine init_particles()
! initialise the particles

use constants
use mod_particles
use Knuth_random
use mod_global_parameters

double precision, dimension(ndir)    :: x
integer                              :: igrid_particle, ipe_particle
integer, parameter                   :: Npart=50000
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
follow(1:10)   = .true.


! first find ipe and igrid responsible for particle
x(:)=0.0d0

do while (nparticles .lt. Npart)

{^D&x(^D) = xprobmin^D/1.5d0 + r^D(nparticles+1) * (xprobmax^D - xprobmin^D)/1.5d0\}

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

use mod_global_parameters

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

use mod_global_parameters
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'jz'//' '//'curlbz'//' '//'v1'//' '//'v2'//' '//'v3'//' '//'divE'
wnames=TRIM(wnames)//' '//'jz'//' '//'curlbz'//' '//'v1'//' '//'v2'//' '//'v3'//' '//'divE'

end subroutine specialvarnames_output
!=============================================================================

















!=============================================================================
!================= JUST DUMMIES ==============================================
!=============================================================================
subroutine flag_grid_usr(qt,ixG^L,ixO^L,w,x,flag)

use mod_global_parameters

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


use mod_global_parameters

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

use mod_global_parameters
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

use mod_global_parameters

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

use mod_global_parameters

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

use mod_global_parameters
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
use mod_global_parameters

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
use mod_global_parameters

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
SUBROUTINE JY01A(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)

!       =======================================================
!       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
!                Y1(x), and their derivatives
!       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ò 0 )
!       Output:  BJ0 --- J0(x)
!                DJ0 --- J0'(x)
!                BJ1 --- J1(x)
!                DJ1 --- J1'(x)
!                BY0 --- Y0(x)
!                DY0 --- Y0'(x)
!                BY1 --- Y1(x)
!                DY1 --- Y1'(x)
!       =======================================================
!
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

double precision:: X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1
double precision:: PI,RP2,X2,R,EC,CS0,W0,R0,CS1,W1,R1
double precision:: T1,P0,Q0,T2,P1,Q1,CU
double precision:: A(12),B(12),A1(12),B1(12)

integer:: K,K0

        PI=3.141592653589793D0
        RP2=0.63661977236758D0
        X2=X*X
        IF (X==0.0D0) THEN
           BJ0=1.0D0
           BJ1=0.0D0
           DJ0=0.0D0
           DJ1=0.5D0
           BY0=-1.0D+300
           BY1=-1.0D+300
           DY0=1.0D+300
           DY1=1.0D+300
           RETURN
        ENDIF
        IF (X<=12.0D0) THEN
           BJ0=1.0D0
           R=1.0D0
           DO K=1,30
              R=-0.25D0*R*X2/(K*K)
              BJ0=BJ0+R
              IF (DABS(R)<DABS(BJ0)*1.0D-15) EXIT
           ENDDO
           BJ1=1.0D0
           R=1.0D0
           DO K=1,30
              R=-0.25D0*R*X2/(K*(K+1.0D0))
              BJ1=BJ1+R
              IF (DABS(R)<DABS(BJ1)*1.0D-15) EXIT
           ENDDO
           BJ1=0.5D0*X*BJ1
           EC=DLOG(X/2.0D0)+0.5772156649015329D0
           CS0=0.0D0
           W0=0.0D0
           R0=1.0D0
           DO K=1,30
              W0=W0+1.0D0/K
              R0=-0.25D0*R0/(K*K)*X2
              R=R0*W0
              CS0=CS0+R
              IF (DABS(R)<DABS(CS0)*1.0D-15) EXIT
           ENDDO
           BY0=RP2*(EC*BJ0-CS0)
           CS1=1.0D0
           W1=0.0D0
           R1=1.0D0
           DO K=1,30
              W1=W1+1.0D0/K
              R1=-0.25D0*R1/(K*(K+1))*X2
              R=R1*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS1=CS1+R
              IF (DABS(R)<DABS(CS1)*1.0D-15) EXIT
           ENDDO
           BY1=RP2*(EC*BJ1-1.0D0/X-0.25D0*X*CS1)
        ELSE
         A(1)=-.7031250000000000D-01
         A(2)=.1121520996093750D+00
         A(3)=-.5725014209747314D+00
         A(4)=.6074042001273483D+01
         A(5)=-.1100171402692467D+03
         A(6)=.3038090510922384D+04
         A(7)=-.1188384262567832D+06
         A(8)=.6252951493434797D+07
         A(9)=-.4259392165047669D+09
         A(10)=.3646840080706556D+11
         A(11)=-.3833534661393944D+13
         A(12)=.4854014686852901D+15
!DATA A/-.7031250000000000D-01,.1121520996093750D+00,-.5725014209747314D+00,.6074042001273483D+01, &
!       -.1100171402692467D+03,.3038090510922384D+04,-.1188384262567832D+06,.6252951493434797D+07, &
!       -.4259392165047669D+09,.3646840080706556D+11,-.3833534661393944D+13,.4854014686852901D+15/
         B(1)=.7324218750000000D-01
         B(2)=-.2271080017089844D+00
         B(3)=.1727727502584457D+01
         B(4)=-.2438052969955606D+02
         B(5)=.5513358961220206D+03
         B(6)=-.1825775547429318D+05
         B(7)=.8328593040162893D+06
         B(8)=-.5006958953198893D+08
         B(9)=.3836255180230433D+10
         B(10)=-.3649010818849833D+12
         B(11)=.4218971570284096D+14
         B(12)=-.5827244631566907D+16
!DATA B/ .7324218750000000D-01,-.2271080017089844D+00,.1727727502584457D+01,-.2438052969955606D+02, &
!        .5513358961220206D+03,-.1825775547429318D+05,.8328593040162893D+06,-.5006958953198893D+08, &
!        .3836255180230433D+10,-.3649010818849833D+12,.4218971570284096D+14,-.5827244631566907D+16/
         A1(1)=.1171875000000000D+00
         A1(2)=-.1441955566406250D+00
         A1(3)=.6765925884246826D+00
         A1(4)=-.6883914268109947D+01
         A1(5)=.1215978918765359D+03
         A1(6)=-.3302272294480852D+04
         A1(7)=.1276412726461746D+06
         A1(8)=-.6656367718817688D+07
         A1(9)=.4502786003050393D+09
         A1(10)=-.3833857520742790D+11
         A1(11)=.4011838599133198D+13
         A1(12)=-.5060568503314727D+15
!DATA A1/.1171875000000000D+00,-.1441955566406250D+00,.6765925884246826D+00,-.6883914268109947D+01, &
!        .1215978918765359D+03,-.3302272294480852D+04,.1276412726461746D+06,-.6656367718817688D+07, &
!        .4502786003050393D+09,-.3833857520742790D+11,.4011838599133198D+13,-.5060568503314727D+15/
         B1(1)=-.1025390625000000D+00
         B1(2)=.2775764465332031D+00
         B1(3)=-.1993531733751297D+01
         B1(4)=.2724882731126854D+02
         B1(5)=-.6038440767050702D+03
         B1(6)=.1971837591223663D+05
         B1(7)=-.8902978767070678D+06
         B1(8)=.5310411010968522D+08
         B1(9)=-.4043620325107754D+10
         B1(10)=.3827011346598605D+12
         B1(11)=-.4406481417852278D+14
         B1(12)=.6065091351222699D+16
!DATA B1/-.1025390625000000D+00,.2775764465332031D+00,-.1993531733751297D+01,.2724882731126854D+02, &
!        -.6038440767050702D+03,.1971837591223663D+05,-.8902978767070678D+06,.5310411010968522D+08, &
!        -.4043620325107754D+10,.3827011346598605D+12,-.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (X>=35.0) K0=10
           IF (X>=50.0) K0=8
           T1=X-0.25D0*PI
           P0=1.0D0
           Q0=-0.125D0/X
           DO K=1,K0
              P0=P0+A(K)*X**(-2*K)
              Q0=Q0+B(K)*X**(-2*K-1)
           ENDDO
           CU=DSQRT(RP2/X)
           BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
           BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
           T2=X-0.75D0*PI
           P1=1.0D0
           Q1=0.375D0/X
           DO K=1,K0
              P1=P1+A1(K)*X**(-2*K)
              Q1=Q1+B1(K)*X**(-2*K-1)
           ENDDO
           CU=DSQRT(RP2/X)
           BJ1=CU*(P1*DCOS(T2)-Q1*DSIN(T2))
           BY1=CU*(P1*DSIN(T2)+Q1*DCOS(T2))
        ENDIF
        DJ0=-BJ1
        DJ1=BJ0-BJ1/X
        DY0=-BY1
        DY1=BY0-BY1/X
        RETURN
        END
!=============================================================================
! amrvacusr.t.nul
!=============================================================================
