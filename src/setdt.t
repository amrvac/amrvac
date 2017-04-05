!>setdt  - set dt for all levels between levmin and levmax. 
!>         dtpar>0  --> use fixed dtpar for all level
!>         dtpar<=0 --> determine CFL limited timestep 
subroutine setdt()
use mod_global_parameters
use mod_physics, only: phys_get_dt, phys_get_aux
use mod_usr_methods, only: usr_get_dt
use mod_thermal_conduction

integer :: iigrid, igrid, ncycle, ncycle2, ifile
double precision :: dtnew, qdtnew, dtmin_mype, factor, dx^D, dxmin^D

double precision :: dtmax, dxmin
integer,save :: stepflag
!----------------------------------------------------------------------------

if(it==0) stepflag = 0

if (dtpar<=zero) then
   dtmin_mype=bigdouble
   cmax_mype = zero
!$OMP PARALLEL DO PRIVATE(igrid,qdtnew,dtnew,dx^D)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      dtnew=bigdouble
      dx^D=rnode(rpdx^D_,igrid);
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      saveigrid = igrid
      block=>pw(igrid)
      block%iw0=0

      if (nwaux>0) then
         call phys_get_aux(.true.,pw(igrid)%w,&
              pw(igrid)%x,ixG^LL,ixM^LL,'setdt')
      end if

      call getdt_courant(pw(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,pw(igrid)%x)
      dtnew=min(dtnew,qdtnew)

      call phys_get_dt(pw(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,pw(igrid)%x)
      dtnew=min(dtnew,qdtnew)

      if (associated(usr_get_dt)) then
         call usr_get_dt(pw(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,pw(igrid)%x)
      end if

      dtnew          = min(dtnew,qdtnew)
      dtmin_mype     = min(dtmin_mype,dtnew)
      dt_grid(igrid) = dtnew
   end do
!$OMP END PARALLEL DO
else
   dtmin_mype=dtpar
end if

if (dtmin_mype<dtmin) then
   write(unitterm,*)"Warning: Time step too small!", dtmin_mype
   write(unitterm,*)"on processor:", mype
   write(unitterm,*)"at time:", global_time," step:", it
   call mpistop("too small timestep")
end if

if (slowsteps>it-itmin+1) then
   factor=one-(one-dble(it-itmin+1)/dble(slowsteps))**2
   dtmin_mype=dtmin_mype*factor
end if

if( stepflag<1.and.mype==0) then
   if(any(dtsave(1:nfile)<dtmin_mype )) then
      write(unitterm,1001) dtmin_mype, dtsave(1:nfile)
      stepflag = 1     
   endif
endif   

dtmin_mype=min(dtmin_mype,time_max-global_time)

if(any(dtsave(1:nfile)<bigdouble).or.any(tsave(isavet(1:nfile),1:nfile)<bigdouble))then
   dtmax = minval((int(global_time/dtsave(1:nfile))+1)*dtsave(1:nfile))-global_time
   do ifile=1,nfile
      dtmax = min(tsave(isavet(ifile),ifile)-global_time,dtmax)
   end do
   if(dtmax<dtmin_mype .and. dtmax > smalldouble)then 
     dtmin_mype=min(dtmin_mype,dtmax)
   end if      
end if

if (dtpar<=zero) then
   call MPI_ALLREDUCE(dtmin_mype,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                      icomm,ierrmpi)
else
   dt=dtmin_mype
end if
   
! estimate time step of thermal conduction
if(associated(phys_getdt_heatconduct)) then
   dtmin_mype=bigdouble
!$OMP PARALLEL DO PRIVATE(igrid,qdtnew,&
!$OMP& dx^D) REDUCTION(min:dt_mype)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      dx^D=rnode(rpdx^D_,igrid);
      saveigrid = igrid
      block=>pw(igrid)
      qdtnew=bigdouble
      call phys_getdt_heatconduct(pw(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,pw(igrid)%x)
      dtmin_mype=min(dtmin_mype,qdtnew)
   end do
!$OMP END PARALLEL DO
   call MPI_ALLREDUCE(dtmin_mype,dtnew,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                         icomm,ierrmpi)
   ncycle=ceiling(dt/dtnew)
   if (ncycle>ncyclemax) then
     if(mype==0 .and. .true.) then
       write(*,*) 'CLF time step is too many times larger than conduction time step',ncycle
       write(*,*) 'reducing dt to',ncyclemax,'times of dt_impl!!'
     endif
     dt=ncyclemax*dtnew
   endif
  ! get number of sub-steps of supertime stepping (Meyer 2012 MNRAS 422,2102)
   if(dt/dtnew< 0.5d0) then
     s=1
   else if(dt/dtnew< 2.d0) then
     s=2
   else
     s=ceiling((dsqrt(9.d0+8.d0*dt/dtnew)-1.d0)/2.d0)
     ! only use odd s number
     s=s/2*2+1
   endif
   dt_tc=dt*0.5d0
   if(mype==0 .and. .false.) write(*,*) 'supertime steps:',s,' normal subcycles:',&
                               ceiling(dt/dtnew/2.d0)
endif

!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   dt_grid(igrid)=dt
end do
!$OMP END PARALLEL DO
     

! global Lax-Friedrich finite difference flux splitting needs fastest wave-speed
! so does GLM: 
call MPI_ALLREDUCE(cmax_mype,cmax_global,1,&
     MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)

1001 format(' Warning: timesteps: ',1x,1pe12.5,' exceeding output intervals ',2(1x,1pe12.5))

end subroutine setdt
!=============================================================================
!> compute CFL limited dt (for variable time stepping)
subroutine getdt_courant(w,ixG^L,ix^L,dtnew,dx^D,x)


use mod_global_parameters
use mod_physics, only: phys_get_cmax

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

integer :: idims
logical :: new_cmax
double precision :: courantmax, dxinv(1:ndim), courantmaxtot, courantmaxtots
double precision :: cmax(ixG^T), cmaxtot(ixG^T), tmp(ixG^T)
!-----------------------------------------------------------------------------
dtnew=bigdouble

courantmax=zero
courantmaxtot=zero
courantmaxtots=zero


new_cmax=.true.

^D&dxinv(^D)=one/dx^D;

cmaxtot(ix^S)=zero

do idims=1,ndim
   call phys_get_cmax(w,x,ixG^L,ix^L,idims,cmax)
   cmax_mype = max(cmax_mype,maxval(cmax(ix^S)))
   if (.not.slab) then
      tmp(ix^S)=cmax(ix^S)/block%dx(ix^S,idims)
      cmaxtot(ix^S)=cmaxtot(ix^S)+tmp(ix^S)
      courantmax=max(courantmax,maxval(tmp(ix^S)))
   else
      cmaxtot(ix^S)=cmaxtot(ix^S)+cmax(ix^S)*dxinv(idims)
      courantmax=max(courantmax,maxval(cmax(ix^S)*dxinv(idims)))
   end if
   courantmaxtot=courantmaxtot+courantmax
end do


select case (typecourant)
case ('minimum')
   ! courantmax='max(c/dx)'
   if (courantmax>smalldouble)     dtnew=min(dtnew,courantpar/courantmax)
case ('summax')
   ! courantmaxtot='summed max(c/dx)'
   if (courantmaxtot>smalldouble)  dtnew=min(dtnew,courantpar/courantmaxtot)
case ('maxsum')
   ! courantmaxtots='max(summed c/dx)'
   courantmaxtots=max(courantmaxtots,maxval(cmaxtot(ix^S)))
   if (courantmaxtots>smalldouble) dtnew=min(dtnew,courantpar/courantmaxtots)
case default
   write(unitterm,*)'Unknown typecourant=',typecourant
   call mpistop("Error from getdt_courant: no such typecourant!")
end select

end subroutine getdt_courant

