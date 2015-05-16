!=============================================================================
! amrvacusr.t.randomplasma
! setup.pl -d=12 -phi=0 -z=0 -g=14 -p=mhd -eos=default -nf=0 -ndust=0 -u=nul -arch=default
!=============================================================================
INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'

integer,dimension(:),allocatable:: seed
integer::  seed_size,ix,imode,nx
real:: randphase(1:1000)
double precision,dimension(:),allocatable:: xgrid,rhogrid

!-----------------------------------------------------------------------------
! pure info, unusued thus far (all dimensionless!)
!! mksA Unit
k_B=1.3806d-23      ! J*K^-1
miu0=4.d0*dpi*1.d-7 ! Henry m^-1

Lunit=1.d6          ! m (1 Megameter)
runit=7.978d-13     ! kg*m^-3
Bunit=10.0d-4       ! (10 G=) 0.001 Tesla

punit=Bunit**2/miu0           ! pressure unit
vunit=Bunit/dsqrt(miu0*runit) ! alfven speed
tunit=Lunit/vunit             ! alfven crossing
!!

eqpar(eta_) = 0.d0
eqpar(gamma_) = 5.0d0/3.0d0

eqpar(Lx0_) = (xprobmax1-xprobmin1)

eqpar(x0_) = 0.0d0 
eqpar(beta_) = 0.01d0 ! plasma beta
eqpar(a0_) = 0.005d0 ! initial amplitude 
eqpar(delrho_) = 10.0d0 ! density contrast 1-15
eqpar(sig0_) = 0.5d0 ! intial width sigma not FWHM
eqpar(nxmodes_) = 400 ! number of modes -> correlation length


if(eqpar(nxmodes_)>1000) call mpistop('too many modes, recompile for more than 1000')
if(mype==0)then
 !   print *,'number of modes=',eqpar(nxmodes_)
 !   write(*,*)'seeding random number generator, on mype==',mype
    call random_seed(SIZE=seed_size)
    allocate(seed(seed_size))
    call random_seed(GET=seed(1:seed_size))
    call random_number(randphase(1:nint(eqpar(nxmodes_))))
 !   write(*,*)'random numbers are:',randphase(1:nint(eqpar(nxmodes_)))
    randphasey(1:nint(eqpar(nxmodes_)))=-dpi+two*dpi*dble(randphase(1:nint(eqpar(nxmodes_))))
    call random_number(randphase(1:nint(eqpar(nxmodes_))))
  !  write(*,*)'random numbers are:',randphase(1:nint(eqpar(nxmodes_)))
    randampliy(1:nint(eqpar(nxmodes_)))=dble(randphase(1:nint(eqpar(nxmodes_))))
    open(123,file='phaseinfo',form='formatted')
    write(123,*) nint(eqpar(nxmodes_))
      do ix=1,nint(eqpar(nxmodes_))
        write(123,"(i4,2es12.4)") ix,randphasey(ix),randampliy(ix)
      enddo
    close(123)
endif
call MPI_BARRIER(icomm,ierrmpi)
if(npe>1)then
    call MPI_BCAST(randphasey,1000,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(randampliy,1000,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    if (norm_rho) then 
      call MPI_BCAST(rho_mean0,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(rho_std,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif
endif

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

include 'amrvacdef.f'

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: rho(ixG^T),pulse(ixG^T)
integer :: imode

logical :: patchw(ixG^T)
logical, save :: first=.true.
!----------------------------------------------------------------------------

if(first)then
  if(mype==0) then
    write(*,*)'Simulating 1D MHD fast pulse through random medium'
    write(*,*)'Case #:',iprob
  endif
  first=.false.
endif

rho(ixG^T) = one 
!!! Note 1/4 (quarter) factor in this equation
  do imode=1,nint(eqpar(nxmodes_))
    rho(ixG^T) =MAX(rho(ixG^T) + eqpar(delrho_)/dble(eqpar(nxmodes_)) &
         *randampliy(imode)*dsin(quarter*dpi*dble(imode)*(x(ixG^T,1)/eqpar(Lx0_))+randphasey(imode)),0.1)
  enddo

pulse(ixG^T)=eqpar(a0_)*dexp(-half*((x(ixG^T,1)-eqpar(x0_))/eqpar(sig0_))**2)

w(ix^S,rho_) = rho(ix^S) + pulse(ix^S) 
w(ix^S,v1_) = pulse(ix^S)
w(ix^S,v2_) = zero
w(ix^S,b1_) = zero
w(ix^S,b2_) = one + pulse(ix^S) 
w(ix^S,p_) = half * eqpar(beta_) * (w(ix^S,b1_)**2+w(ix^S,b2_)**2)

patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

end subroutine initonegrid_usr
!=============================================================================
! amrvacusr.t.randomplasma
!=============================================================================
