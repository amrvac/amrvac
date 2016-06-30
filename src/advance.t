!=============================================================================
subroutine advance(iit)

{#IFDEF PARTICLES
use mod_particles, only: tmax_particles
use mod_timing, only: tpartc, tpartc0
}
include 'amrvacdef.f'

integer, intent(in) :: iit

integer :: iigrid, igrid, idimsplit
{#IFDEF TCRKL2
integer :: s
}
logical :: firstsweep, lastsweep
common/first/firstsweep,lastsweep
!-----------------------------------------------------------------------------

! when computing to steady state, store old solution
if(residmin>smalldouble)then
!$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     pwres(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
  end do
!$OMP END PARALLEL DO
endif

{#IFDEF RAY
call update_rays
}

! split source addition
if(ssplitdust .or. ssplitdivb .or. ssplitresis .or. ssplituser) &
  call addsource_all(.true.)

{#IFDEF TCRKL2
! split thermal conduction source addition 1/2
if(conduction) then
  if(dt/dtimpl < 0.5d0) then
    s=1
    else if(dt/dtimpl < 2.d0) then
    s=2
    else
    s=ceiling((dsqrt(9.d0+8.d0*dt/dtimpl)-1.d0)/2.d0)
    ! only use odd s number
    s=s/2*2+1
  endif
  if(mype==0 .and. .false.) print*,'supertime steps:', s, &
              ' normal subcycles:',ceiling(dt/dtimpl/2.d0), &
              'time step ratio:',dt/dtimpl
  call thermalconduct_RKL2(s,dt/2.d0,t) 
endif
}
! old solution values at t_n-1 no longer needed: make copy of w(t_n)
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   pwold(igrid)%w(ixG^T,1:nwflux+nwaux)=pw(igrid)%w(ixG^T,1:nwflux+nwaux)
end do
!$OMP END PARALLEL DO

! split implicit source addition
if (sourceimpl) then
   if(mype==0.and..false.) print *,'implicit source addition'
   call addsource_impl
endif

firstsweep=.true.
if (dimsplit) then
   if ((iit/2)*2==iit .or. typedimsplit=='xy') then
      ! do the sweeps in order of increasing idim,
      do idimsplit=1,ndim
         lastsweep= idimsplit==ndim
         call advect(idimsplit,idimsplit)
      end do
   else
      ! If the parity of "iit" is odd and typedimsplit=xyyx,
      ! do sweeps backwards
      do idimsplit=ndim,1,-1
         lastsweep= idimsplit==1
         call advect(idimsplit,idimsplit)
      end do
   end if
else
   ! Add fluxes from all directions at once
   lastsweep= .true.
   call advect(1,ndim)
end if


! split source addition
if(ssplitdust .or. ssplitdivb .or. ssplitresis .or. ssplituser) &
  call addsource_all(.false.)

{#IFDEF TCRKL2
! split thermal conduction source addition 2/2
if(conduction) then
  call thermalconduct_RKL2(s,dt/2.d0,t+dt) 
endif
}

{#IFDEF PARTICLES
tpartc0 = MPI_WTIME()
tmax_particles = (t + dt)* (UNIT_LENGTH/UNIT_VELOCITY)
call handle_particles
tpartc = tpartc + (MPI_WTIME() - tpartc0)
}

if(residmin>smalldouble)then
!$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     pwres(igrid)%w(ixG^T,1:nwflux)= &
        pw(igrid)%w(ixG^T,1:nwflux)-pwres(igrid)%w(ixG^T,1:nwflux)
  end do
!$OMP END PARALLEL DO
endif

end subroutine advance
!=============================================================================
subroutine advect(idim^LIM)

!  integrate all grids by one step of its delta(t)

! This subroutine is in VAC terminology equivalent to
! `advect' (with the difference that it will `advect' all grids)

include 'amrvacdef.f'

integer, intent(in) :: idim^LIM

integer :: iigrid, igrid
logical :: firstsweep, lastsweep
common/first/firstsweep,lastsweep
!-----------------------------------------------------------------------------
! copy w instead of wold because of potential use of dimsplit or sourcesplit
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   allocate (pw1(igrid)%w(ixG^T,1:nw))
   pw1(igrid)%w=pw(igrid)%w
end do
!$OMP END PARALLEL DO

istep=0

select case (typeadvance)
case ("onestep")
   call advect1(typefull1,one,    idim^LIM,t,          pw1,t,pw, pwold)
case ("twostep")
   ! predictor step
   call advect1(typepred1,half,   idim^LIM,t,          pw, t,pw1,pwold)
   ! corrector step
   call advect1(typefull1,one,    idim^LIM,t+half*dt,  pw1,t,pw, pwold)
case ("threestep")
   ! three step Runge-Kutta in accordance with Gottlieb & Shu 1998
   call advect1(typefull1,one, idim^LIM,t,          pw ,t,pw1,pwold)

   do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate (pw2(igrid)%w(ixG^T,1:nw))
      pw2(igrid)%w(ixG^T,1:nwflux)=0.75d0*pw(igrid)%w(ixG^T,1:nwflux)+0.25d0*&
        pw1(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw2(igrid)%w(ixG^T,nwflux+1:nw) = &
                       pw(igrid)%w(ixG^T,nwflux+1:nw)
   end do

   call advect1(typefull1,0.25d0, idim^LIM,t+dt,pw1,t+dt*0.25d0,pw2,pwold)

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux)=1.0d0/3.0d0*pw(igrid)%w(ixG^T,1:nwflux)+&
        2.0d0/3.0d0*pw2(igrid)%w(ixG^T,1:nwflux)
   end do   
!$OMP END PARALLEL DO
   call advect1(typefull1,2.0d0/3.0d0, idim^LIM,t+dt/2.0d0,pw2,t+dt/3.0d0,pw,&
          pwold)

case ("ssprk43")
! Strong stability preserving 4 stage RK 3rd order method by Ruuth and Spiteri
!
! Ruuth & Spiteri
! J. S C, 17 (2002) p. 211 - 220
!
! supposed to be stable up to CFL=2.
! don't use time-dependent sources since i did not bother to set those intermediate times.
! oliver.

! === First step ===
   call advect1(typefull1,0.5d0, idim^LIM,t,          pw ,t,pw1,pwold)

! === Second step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate (pw2(igrid)%w(ixG^T,1:nw))
      pw2(igrid)%w(ixG^T,1:nwflux)=pw1(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw2(igrid)%w(ixG^T,nwflux+1:nw) = &
                       pw(igrid)%w(ixG^T,nwflux+1:nw)
   end do
   call advect1(typefull1,0.5d0, idim^LIM,t,pw1,t+dt,pw2,pwold)

! === Third step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate (pw3(igrid)%w(ixG^T,1:nw))
      pw3(igrid)%w(ixG^T,1:nwflux)=2.0d0/3.0d0 * pw(igrid)%w(ixG^T,1:nwflux) &
           + 1.0d0/3.0d0 * pw2(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw3(igrid)%w(ixG^T,nwflux+1:nw) = &
                       pw(igrid)%w(ixG^T,nwflux+1:nw)
   end do
   call advect1(typefull1,1.0d0/6.0d0, idim^LIM,t,pw2,t+dt,pw3,pwold)

! === Fourth step ===
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux)=pw3(igrid)%w(ixG^T,1:nwflux)
   end do   
!$OMP END PARALLEL DO
   call advect1(typefull1,0.5d0, idim^LIM, t,pw3, t+dt,pw, pwold)

case ("ssprk54")
! Strong stability preserving 5 stage RK 4th order method by Ruuth and Spiteri
!
! SIAM J. NUMER. ANAL.
! Vol. 40, No. 2, pp. 469–491
! c 2002 Society for Industrial and Applied Mathematics
! A NEW CLASS OF OPTIMAL
! HIGH-ORDER STRONG-STABILITY-PRESERVING
! TIME DISCRETIZATION METHODS
!
! E.g. Table A.2
!
! I have the actual coefficients however from the overview article by Gottlieb, JoSC 25 (2005)
! ON HIGH ORDER STRONG STABILITY PRESERVING RUNGE-KUTTA AND MULTI STEP TIME DISCRETIZATIONS
!
! there are slight differences in the coefficients (~8th digit after the .)
! This is SSP till CFL number 1.508 which makes the effective CFL per step ceff=0.377, 
! in contrast to the ceff=0.3333 of the classical RK3.  
!
! coded by oliver on 11/05/2013.  Enjoy!

! === First step ===
   call advect1(typefull1,0.391752226571890d0, idim^LIM, t,pw ,t,pw1, pwold)

! === Second step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate (pw2(igrid)%w(ixG^T,1:nw))
      pw2(igrid)%w(ixG^T,1:nwflux)=0.444370493651235d0 * pw(igrid)%w(ixG^T,1:nwflux) &
           + 0.555629506348765d0 * pw1(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw2(igrid)%w(ixG^T,nwflux+1:nw) = &
                       pw(igrid)%w(ixG^T,nwflux+1:nw)
   end do
   call advect1(typefull1,0.368410593050371d0, idim^LIM, t+0.391752226571890d0*dt,pw1, t+0.2176690962611688d0*dt,pw2, pwold)

! === Third step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate (pw3(igrid)%w(ixG^T,1:nw))
      pw3(igrid)%w(ixG^T,1:nwflux)=0.620101851488403d0 * pw(igrid)%w(ixG^T,1:nwflux) &
           + 0.379898148511597d0 * pw2(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw3(igrid)%w(ixG^T,nwflux+1:nw) = &
                       pw(igrid)%w(ixG^T,nwflux+1:nw)
   end do
   call advect1(typefull1,0.251891774271694d0, idim^LIM, t+0.5860796893115398d0*dt,pw2, t+0.222650588849706d0*dt,pw3, pwold)

! === Fourth step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate (pw4(igrid)%w(ixG^T,1:nw))
      pw4(igrid)%w(ixG^T,1:nwflux)=0.178079954393132d0 * pw(igrid)%w(ixG^T,1:nwflux) &
           + 0.821920045606868d0 * pw3(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw4(igrid)%w(ixG^T,nwflux+1:nw) = &
                       pw(igrid)%w(ixG^T,nwflux+1:nw)
   end do
   call advect1(typefull1,0.544974750228521d0, idim^LIM, t+0.4745423631214d0*dt,pw3, t+0.390035880739132d0*dt,pw4, pwold)
! Now recover back the dt*L(u3), store in pw1:
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixG^T,1:nwflux) = ( pw4(igrid)%w(ixG^T,1:nwflux) &
           - (0.178079954393132d0 * pw(igrid)%w(ixG^T,1:nwflux) &
           + 0.821920045606868d0 * pw3(igrid)%w(ixG^T,1:nwflux)) ) / 0.544974750228521d0
   end do
!$OMP END PARALLEL DO

! === Fifth step ===
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux)= 0.517231671970585d0 * pw2(igrid)%w(ixG^T,1:nwflux) &
           + 0.096059710526147d0 * pw3(igrid)%w(ixG^T,1:nwflux) &
           + 0.063692468666290d0 * pw1(igrid)%w(ixG^T,1:nwflux) & 
           + 0.386708617503269d0 * pw4(igrid)%w(ixG^T,1:nwflux)
   end do
!$OMP END PARALLEL DO
   call advect1(typefull1,0.226007483236906d0, idim^LIM, t+0.935010630967653d0*dt,pw4, t+0.710300048096804d0*dt,pw, pwold)


case ("rk4")
   ! classical RK4 four step scheme
   call advect1(typefull1,0.5d0, idim^LIM,t,          pw ,t,pw1,pwold)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate (pw2(igrid)%w(ixG^T,1:nw))
      pw2(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
   end do
   call advect1(typefull1,0.5d0, idim^LIM,t+dt/2d0,   pw1,t,pw2,pwold)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate (pw3(igrid)%w(ixG^T,1:nw))
      pw3(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
   end do
   call advect1(typefull1,one,   idim^LIM,t+dt/2d0,   pw2,t,pw3,pwold)
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixG^T,1:nwflux)=(pw1(igrid)%w(ixG^T,1:nwflux) &
                               +two*pw2(igrid)%w(ixG^T,1:nwflux) &
                                   +pw3(igrid)%w(ixG^T,1:nwflux) &
                              -4.0d0*pw(igrid)%w(ixG^T,1:nwflux))/3.0d0
   end do
!$OMP END PARALLEL DO
   call advect1(typefull1,1.0d0/6.0d0, idim^LIM,t+dt,  pw3,t,pw,pwold)
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux)=pw1(igrid)%w(ixG^T,1:nwflux)+&
        pw(igrid)%w(ixG^T,1:nwflux)
   end do
!$OMP END PARALLEL DO
case ("fourstep")
   ! four step scheme, variant Hans De Sterck
   call advect1(typefull1,0.12d0, idim^LIM,t,          pw ,t,pw1,pwold)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate (pw2(igrid)%w(ixG^T,1:nw))
      pw2(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
   end do
   call advect1(typefull1,0.25d0, idim^LIM,t+dt*0.12d0,pw1,t,pw2,pwold)
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
   end do
!$OMP END PARALLEL DO
   call advect1(typefull1,0.5d0,  idim^LIM,t+dt/4d0   ,pw2,t,pw1,pwold)
   call advect1(typefull1,one,    idim^LIM,t+dt/2d0   ,pw1,t,pw, pwold)
case ("jameson")
   ! four step scheme, variant jameson
   call advect1(typefull1,0.25d0, idim^LIM,t,          pw ,t,pw1,pwold)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate (pw2(igrid)%w(ixG^T,1:nw))
      pw2(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
   end do
   call advect1(typefull1,(1.0d0/3.0d0), idim^LIM,t+dt*0.25d0,pw1,t,pw2,pwold)
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
   end do
!$OMP END PARALLEL DO
   call advect1(typefull1,0.5d0,  idim^LIM,t+dt/3d0   ,pw2,t,pw1,pwold)
   call advect1(typefull1,one,    idim^LIM,t+dt/2d0   ,pw1,t,pw, pwold)
case default
   write(unitterm,*) "typeadvance=",typeadvance
   write(unitterm,*) "Error in advect: Unknown time integration method"
   call mpistop("Correct typeadvance")
end select

do iigrid=1,igridstail; igrid=igrids(iigrid);
   deallocate (pw1(igrid)%w)
   select case (typeadvance)
     case ("threestep","fourstep","jameson")
       deallocate (pw2(igrid)%w)
     case ("rk4","ssprk43")
       deallocate (pw2(igrid)%w)
       deallocate (pw3(igrid)%w)
     case ("ssprk54")
       deallocate (pw2(igrid)%w)
       deallocate (pw3(igrid)%w)
       deallocate (pw4(igrid)%w)
   end select
end do


firstsweep=.false.
end subroutine advect
!=============================================================================
subroutine advect1(method,dtfactor,idim^LIM,qtC,pwa,qt,pwb,pwc)

!  integrate all grids by one partial step

! This subroutine is equivalent to VAC's `advect1', but does
! the advection for all grids
include 'amrvacdef.f'

integer, intent(in) :: idim^LIM
double precision, intent(in) :: dtfactor, qtC, qt
character(len=*), intent(in) :: method(nlevelshi)
type(walloc) :: pwa(ngridshi), pwb(ngridshi), pwc(ngridshi)

double precision :: qdt
integer :: iigrid, igrid, level

logical :: setigrid

!-----------------------------------------------------------------------------
istep=istep+1

if (time_advance.and.levmax>levmin) then
   if (istep==nstep.or.nstep>2) &
        call init_comm_fix_conserve(idim^LIM)
end if

! loop over all grids to arrive at equivalent
! VAC subroutine advect1==advect1_grid 

! opedit: Just advance the active grids: 
!$OMP PARALLEL DO PRIVATE(igrid,level,qdt)
     do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      level=node(plevel_,igrid)
      qdt=dtfactor*dt_grid(igrid)
      
      call process1_grid(method(level),igrid,qdt,ixG^LL,idim^LIM,qtC,&
                      pwa(igrid)%w,qt,pwb(igrid)%w,pwc(igrid)%w)
                      
      end do
!$OMP END PARALLEL DO

! opedit: Send flux for all grids, expects sends for all 
! nsend_fc(^D), set in connectivity.t.

      if (time_advance.and.levmax>levmin) then
         if (istep==nstep.or.nstep>2) then
            do iigrid=1,igridstail; igrid=igrids(iigrid);
               call sendflux(igrid,idim^LIM)
            end do
            call fix_conserve(pwb,idim^LIM)
         end if
      end if
   
! for all grids: fill ghost cells
      qdt=dtfactor*dt
{#IFDEF BOUNDARYDRIVER
      call boundarydriver(method(mxnest),qdt,idim^LIM,qtC,pwa,qt,pwb)
      !call boundarydriver(qt,qdt,pwa,pwb)
}
      call getbc(qt+qdt,ixG^LL,pwb,pwCoarse,pgeo,pgeoCoarse,.false.,0,nwflux+nwaux)

end subroutine advect1
!=============================================================================
subroutine process1_grid(method,igrid,qdt,ixG^L,idim^LIM,qtC,wCT,qt,w,wold)

! This subroutine is equivalent to VAC's`advect1' for one grid

include 'amrvacdef.f'

character(len=*), intent(in) :: method
integer, intent(in) :: igrid, ixG^L, idim^LIM
double precision, intent(in) :: qdt, qtC, qt
double precision :: wCT(ixG^S,1:nw), w(ixG^S,1:nw), wold(ixG^S,1:nw)

double precision :: dx^D
double precision :: fC(ixG^S,1:nwflux,1:ndim)
!-----------------------------------------------------------------------------
dx^D=rnode(rpdx^D_,igrid);
^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
saveigrid=igrid
{#IFDEF STRETCHGRID
logG=logGs(node(plevel_,igrid))
qst=qsts(node(plevel_,igrid))
}

if (.not.slab) mygeo => pgeo(igrid)
if (B0field) then
   myB0_cell => pB0_cell(igrid)
   {^D&myB0_face^D => pB0_face^D(igrid)\}
end if
typelimiter=typelimiter1(node(plevel_,igrid))
typegradlimiter=typegradlimiter1(node(plevel_,igrid))

call advect1_grid(method,qdt,ixG^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D, &
                  px(igrid)%x)


! opedit: Obviously, flux is stored only for active grids.  
! but we know in fix_conserve wether there is a passive neighbor
! via neighbor_active(i^D,igrid) thus we skip the correction for those.
! This violates strict conservation when the active/passive interface 
! coincides with a coarse/fine interface.  
if (time_advance.and.levmax>levmin) then
   if (istep==nstep.or.nstep>2) &
        call storeflux(igrid,fC,idim^LIM)
end if

end subroutine process1_grid
!=============================================================================
subroutine advect1_grid(method,qdt,ixI^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,x)

!  integrate one grid by one partial step 

include 'amrvacdef.f'

character(len=*), intent(in) :: method
integer, intent(in) :: ixI^L, idim^LIM
double precision, intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
double precision :: wCT(ixI^S,1:nw), w(ixI^S,1:nw), wold(ixI^S,1:nw)
double precision :: fC(ixI^S,1:nwflux,1:ndim)

integer :: ixO^L
!-----------------------------------------------------------------------------

ixO^L=ixI^L^LSUBdixB;

select case (method)
case ('cd')
   call centdiff(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,fC,dx^D,x)
case ('cd4')
   call centdiff4(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,fC,dx^D,x)
case ('hancock')
   call hancock(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,dx^D,x)
case ('fd')
   call fd(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,x)
case ('tvdmu','tvdmu1')
   call tvdmusclf(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,x)
case ('tvdlf','tvdlf1')
   call tvdlf(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,x)
case ('hll','hll1')
   call hll(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,x)
case ('hllc','hllc1', 'hllcd','hllcd1')
   call hllc(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,x)
{^IFMHDPHYS
case ('hlld','hlld1', 'hlldd','hlldd1')
   call hlld(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,x)
}
case ('tvd','tvd1')
   call centdiff(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,fC,dx^D,x)
   call tvdlimit(method,qdt,ixI^L,ixO^L,idim^LIM,wCT,qt+qdt,w,fC,dx^D,x)
case ('source')
   call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),&
                                      ixI^L,ixO^L,1,nw,qtC,wCT,qt,w,x,.false.)
case ('nul')
   ! There is nothing to do
case default
   write(unitterm,*)'Error in advect1_grid:',method,' is unknown!'
   call mpistop("")
end select

end subroutine advect1_grid
!=============================================================================
subroutine addsource_all(prior)

include 'amrvacdef.f'

logical, intent(in) :: prior

integer :: iigrid, igrid, i^D
double precision :: qdt, qt
!-----------------------------------------------------------------------------

if ((.not.prior).and.&
    (typesourcesplit=='sf' .or. typesourcesplit=='ssf')) return

if (prior) then
   qt=t
else
   qt=t+dt
end if
!$OMP PARALLEL DO PRIVATE(igrid,qdt,i^D)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   qdt=dt_grid(igrid)
   if (B0field) then
      myB0_cell => pB0_cell(igrid)
      {^D&myB0_face^D => pB0_face^D(igrid)\}
   end if
   {do i^DB=-1,1\}
      if (i^D==0|.and.) cycle
      if (neighbor_type(i^D,igrid)==2 .or. neighbor_type(i^D,igrid)==4) then
         leveljump(i^D)=.true. 
      else
         leveljump(i^D)=.false. 
      end if
   {end do\}
   call addsource1_grid(igrid,qdt,qt,pw(igrid)%w)
end do
!$OMP END PARALLEL DO

call getbc(qt,0.d0,ixG^LL,pw,pwCoarse,pgeo,pgeoCoarse,.false.,0,nwflux+nwaux)

end subroutine addsource_all
!=============================================================================
subroutine addsource1_grid(igrid,qdt,qt,w)

include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: qdt, qt
double precision, intent(inout) :: w(ixG^T,nw)

double precision :: w1(ixG^T,nw)
!-----------------------------------------------------------------------------

if (.not.slab) mygeo => pgeo(igrid)

saveigrid=igrid
typelimiter=typelimiter1(node(plevel_,igrid))
typegradlimiter=typegradlimiter1(node(plevel_,igrid))

^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

w1(ixG^T,1:nw)=w(ixG^T,1:nw)

select case (typesourcesplit)
case ('sf')
   call addsource2(qdt  ,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case ('sfs')
   call addsource2(qdt/2,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case ('ssf')
   call addsource2(qdt/2,ixG^LL,ixG^LL,1,nw,qt,w,qt,w1,px(igrid)%x,.true.)
   call addsource2(qdt  ,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case ('ssfss')
   call addsource2(qdt/4,ixG^LL,ixG^LL,1,nw,qt,w,qt,w1,px(igrid)%x,.true.)
   call addsource2(qdt/2,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case default
   write(unitterm,*)'No such typesourcesplit=',typesourcesplit
   call mpistop("Error: Unknown typesourcesplit!")
end select

end subroutine addsource1_grid
!=============================================================================
subroutine addsource2(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

! Add source within ixO for iws: w=w+qdt*S[wCT]

include 'amrvacdef.f'

! differences with VAC is in iw^LIM and in declaration of ranges for wCT,w

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
logical, intent(in) :: qsourcesplit
!-----------------------------------------------------------------------------

! user defined sources, typically explicitly added
if(qsourcesplit .eqv. ssplituser) &
   call specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! physics defined sources, typically explicitly added,
! along with geometrical source additions
call addsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

end subroutine addsource2
!=============================================================================
subroutine process(iit,qt)

include 'amrvacdef.f'

! .. scalars ..
integer,intent(in)          :: iit
double precision, intent(in):: qt

integer:: iigrid, igrid,level
!-----------------------------------------------------------------------------

{#IFDEF PROCESSGLOBAL
call process_global_usr(iit,qt)
}

!$OMP PARALLEL DO PRIVATE(igrid,level)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   level=node(plevel_,igrid)
   ! next few lines ensure correct usage of routines like divvector etc
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
   if (.not.slab) mygeo => pgeo(igrid)
   if (B0field) then
      myB0_cell => pB0_cell(igrid)
      {^D&myB0_face^D => pB0_face^D(igrid)\}
   end if
   typelimiter=typelimiter1(node(plevel_,igrid))
   typegradlimiter=typegradlimiter1(node(plevel_,igrid))
   call process_grid_usr(igrid,level,ixG^LL,ixM^LL,qt,pw(igrid)%w,px(igrid)%x)
end do
!$OMP END PARALLEL DO

end subroutine process
!=============================================================================
subroutine addsource_impl

include 'amrvacdef.f'

integer :: iigrid, igrid, icycle, ncycle
double precision :: qdt, qt, sumqdt
!-----------------------------------------------------------------------------

bcphys=.false.
if(sourceimplcycle)then
   ncycle=ceiling(dt/dtimpl)
   if(ncycle<1) then
     ncycle=1
     dtimpl=dt
   endif
else
   ncycle=1
endif

if(mype==0.and..false.) then
  print *,'implicit source addition will subcycle with ',ncycle,' subtimesteps'
  print *,'dt and dtimpl= ',dt,dtimpl,' versus ncycle*dtimpl=',ncycle*dtimpl
endif

qt=t
sumqdt=zero
do icycle=1,ncycle
  if(sourceparasts)then
    qdt=dtimpl/(one+parastsnu+(parastsnu-one)* &
         dcos(half*dpi*(two*dble(icycle)-one)/dble(ncycle)))
  else
    qdt=dt/dble(ncycle)
  endif
!$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   if (.not.slab) mygeo => pgeo(igrid)
   if (B0field) then
      myB0_cell => pB0_cell(igrid)
      {^D&myB0_face^D => pB0_face^D(igrid)\}
   end if
   call addsourceimpl1_grid(igrid,qdt,qt,pw(igrid)%w,px(igrid)%x)
  end do
!$OMP END PARALLEL DO

  sumqdt=sumqdt+qdt
  qt=qt+qdt
  call getbc(qt,ixG^LL,pw,pwCoarse,pgeo,pgeoCoarse,.false.,0,nwflux+nwaux)
enddo

if(mype==0.and..false.) then
  if(sourceparasts)then
     print *,'qt-t=',qt-t,'versus ncycle2*dtimpl=',ncycle*ncycle*dtimpl,&
       ' and sumqdt=',sumqdt
  else
     print *,'qt-t=',qt-t,'versus ncycle*dtimpl=',ncycle*dtimpl,&
       ' and sumqdt=',sumqdt
  endif
endif
bcphys=.true.

end subroutine addsource_impl
!=============================================================================
subroutine addsourceimpl1_grid(igrid,qdt,qt,w,x)

include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: qdt, qt, x(ixG^T,1:ndim)
double precision, intent(inout) :: w(ixG^T,nw)

double precision :: w1(ixG^T,nw)
!-----------------------------------------------------------------------------

typelimiter=typelimiter1(node(plevel_,igrid))
typegradlimiter=typegradlimiter1(node(plevel_,igrid))

^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

w1(ixG^T,1:nwflux)=w(ixG^T,1:nwflux)

call specialsource_impl(qdt,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,x)

end subroutine addsourceimpl1_grid
!=============================================================================
