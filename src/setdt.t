!>setdt  - set dt for all levels between levmin and levmax. 
!>         dtpar>0  --> use fixed dtpar for all level
!>         dtpar<=0 --> determine CFL limited timestep 
subroutine setdt()
  use mod_global_parameters
  use mod_physics
  use mod_trac
  use mod_usr_methods, only: usr_get_dt
  use mod_supertimestepping, only: set_dt_sts_ncycles, is_sts_initialized, sourcetype_sts,sourcetype_sts_split

  integer :: iigrid, igrid, ncycle, ncycle2, ifile, idim
  double precision :: dtnew, qdtnew, dtmin_mype, factor, dx^D, dxmin^D

  double precision :: dtmax, dxmin, cmax_mype
  double precision :: a2max_mype(ndim), tco_mype, tco_global, Tmax_mype, T_peak
  double precision :: trac_alfa, trac_dmax, trac_tau

  integer, parameter :: niter_print = 2000

  if (dtpar<=zero) then
    dtmin_mype=bigdouble
    cmax_mype = zero
    a2max_mype = zero
    tco_mype = zero
    Tmax_mype = zero
    !$OMP PARALLEL DO PRIVATE(igrid,qdtnew,dtnew,dx^D) REDUCTION(min:dtmin_mype) REDUCTION(max:cmax_mype,a2max_mype)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      dtnew=bigdouble
      dx^D=rnode(rpdx^D_,igrid);
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      block=>ps(igrid)

      call getdt_courant(ps(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,ps(igrid)%x,&
           cmax_mype,a2max_mype)
      dtnew=min(dtnew,qdtnew)

      call phys_get_dt(ps(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,ps(igrid)%x)
      dtnew=min(dtnew,qdtnew)

      if (associated(usr_get_dt)) then
         call usr_get_dt(ps(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,ps(igrid)%x)
      end if
      dtnew          = min(dtnew,qdtnew)
      dtmin_mype     = min(dtmin_mype,dtnew)
    end do
    !$OMP END PARALLEL DO
  else
     dtmin_mype=dtpar
  end if

  if (dtmin_mype<dtmin) then
     write(unitterm,*)"Error: Time step too small!", dtmin_mype
     write(unitterm,*)"on processor:", mype, "at time:", global_time," step:", it
     write(unitterm,*)"Lower limit of time step:", dtmin
     crash=.true.
  end if

  if (slowsteps>it-it_init+1) then
     factor=one-(one-dble(it-it_init+1)/dble(slowsteps))**2
     dtmin_mype=dtmin_mype*factor
  end if

  if(final_dt_reduction)then
     !if (dtmin_mype>time_max-global_time) then
     !   write(unitterm,*)"WARNING final timestep artificially reduced!"
     !   write(unitterm,*)"on processor:", mype, "at time:", global_time," step:", it
     !endif
     if(time_max-global_time<=dtmin) then
        !write(unitterm,*)'Forcing to leave timeloop as time is reached!'
        final_dt_exit=.true.
     endif
     dtmin_mype=min(dtmin_mype,time_max-global_time)
  end if

  if (dtpar<=zero) then
     call MPI_ALLREDUCE(dtmin_mype,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                        icomm,ierrmpi)
  else
     dt=dtmin_mype
  end if

  if(any(dtsave(1:nfile)<bigdouble).or.any(tsave(isavet(1:nfile),1:nfile)<bigdouble))then
     dtmax = minval(ceiling(global_time/dtsave(1:nfile))*dtsave(1:nfile))-global_time
     do ifile=1,nfile
        dtmax = min(tsave(isavet(ifile),ifile)-global_time,dtmax)
     end do
     if(dtmax > smalldouble)then 
       dt=min(dt,dtmax)
     else
       ! dtmax=0 means dtsave is divisible by global_time
       dt=min(dt,minval(dtsave(1:nfile)))
     end if      
  end if

  if(mype==0) then
    if(any(dtsave(1:nfile)<dt)) then
      write(unitterm,*) 'Warning: timesteps: ',dt,' exceeding output intervals ', dtsave(1:nfile)
    endif
  endif   

  if(is_sts_initialized()) then
    if(sourcetype_sts .eq. sourcetype_sts_split) then
      qdtnew = 0.5d0 * dt 
      if (set_dt_sts_ncycles(qdtnew)) then
        dt = 2.d0*qdtnew
        !a quick way to print the reduction of time only every niter_print iterations
        !Note that niter_print is a parameter variable hardcoded to the value of 200
        if(mype==0 .and. mod(it-1, niter_print) .eq. 0) then
          write(*,*) 'Max number of STS cycles exceeded, reducing dt to',dt
        endif
      endif  
    else
      if(set_dt_sts_ncycles(dt))then 
       if(mype==0 .and. mod(it-1, niter_print) .eq. 0) then
         write(*,*) 'Max number of STS cycles exceeded, reducing dt to',dt
       endif
      endif
    endif
  endif

  ! global Lax-Friedrich finite difference flux splitting needs fastest wave-speed
  ! so does GLM: 
  if(need_global_cmax) call MPI_ALLREDUCE(cmax_mype,cmax_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)
  if(need_global_a2max) call MPI_ALLREDUCE(a2max_mype,a2max_global,ndim,&
       MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)

  ! transition region adaptive thermal conduction (Johnston 2019 ApJL, 873, L22)
  ! transition region adaptive thermal conduction (Johnston 2020 A&A, 635, 168)
  if(phys_trac) then
    T_bott=2.d4/unit_temperature
    call MPI_ALLREDUCE(Tmax_mype,T_peak,1,MPI_DOUBLE_PRECISION,&
         MPI_MAX,icomm,ierrmpi)
    ! default lower limit of cutoff temperature
    select case(phys_trac_type)
    case(0)
      !> Test case, fixed cutoff temperature
      !> do nothing here
    case(1)
      !> 1D TRAC method
      trac_dmax=0.1d0
      trac_tau=1.d0/unit_time
      trac_alfa=trac_dmax**(dt/trac_tau)
      tco_global=zero
      {^IFONED
      call MPI_ALLREDUCE(tco_mype,tco_global,1,MPI_DOUBLE_PRECISION,&
           MPI_MAX,icomm,ierrmpi)
      }
      call TRAC_simple(tco_global,trac_alfa,T_peak)
    case(2)
      !> LTRAC method, by iijima et al. 2021
      !> do nothing here 
      call LTRAC(T_peak)
    case(3)
      !> 2D or 3D TRACL(ine) method
      call TRACL(.false.,T_peak)
    case(4)
      !> 2D or 3D TRACB(lock) method
      call TRACB(.false.,T_peak)
    case(5)
      !> 2D or 3D TRACL(ine) method with mask
      call TRACL(.true.,T_peak)
    case(6)
      !> 2D or 3D TRACB(lock) method with mask
      call TRACB(.true.,T_peak)
    case default
      call mpistop("undefined TRAC method type")
    end select
  end if 

  contains

    !> compute CFL limited dt (for variable time stepping)
    subroutine getdt_courant(w,ixI^L,ixO^L,dtnew,dx^D,x,cmax_mype,a2max_mype)
      use mod_global_parameters
      use mod_physics, only: phys_get_cmax,phys_get_a2max,phys_get_tcutoff

      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in)    :: dx^D
      double precision, intent(inout) :: w(ixI^S,1:nw), dtnew, cmax_mype, a2max_mype(ndim)

      integer :: idims
      double precision :: courantmax, dxinv(1:ndim), courantmaxtot, courantmaxtots
      double precision :: cmax(ixI^S), cmaxtot(ixO^S)
      double precision :: a2max(ndim),tco_local,Tmax_local

      dtnew=bigdouble
      courantmax=zero
      courantmaxtot=zero
      courantmaxtots=zero


      if(need_global_a2max) then
        call phys_get_a2max(w,x,ixI^L,ixO^L,a2max)
        do idims=1,ndim
          a2max_mype(idims) = max(a2max_mype(idims),a2max(idims))
        end do
      end if
      if(phys_trac) then
        call phys_get_tcutoff(ixI^L,ixO^L,w,x,tco_local,Tmax_local)
        {^IFONED tco_mype=max(tco_mype,tco_local) }
        Tmax_mype=max(Tmax_mype,Tmax_local)
      end if

      !if(slab_uniform) then
      !  ^D&dxinv(^D)=one/dx^D;
      !  do idims=1,ndim
      !    call phys_get_cmax(w,x,ixI^L,ixO^L,idims,cmax)
      !    if(need_global_cmax) cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
      !    if(need_global_a2max) a2max_mype(idims) = max(a2max_mype(idims),a2max(idims))
      !    cmaxtot(ixO^S)=cmaxtot(ixO^S)+cmax(ixO^S)*dxinv(idims)
      !    courantmax=max(courantmax,maxval(cmax(ixO^S)*dxinv(idims)))
      !    courantmaxtot=courantmaxtot+courantmax
      !  end do
      !else
      !  do idims=1,ndim
      !    call phys_get_cmax(w,x,ixI^L,ixO^L,idims,cmax)
      !    if(need_global_cmax) cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
      !    if(need_global_a2max) a2max_mype(idims) = max(a2max_mype(idims),a2max(idims))
      !    tmp(ixO^S)=cmax(ixO^S)/block%ds(ixO^S,idims)
      !    cmaxtot(ixO^S)=cmaxtot(ixO^S)+tmp(ixO^S)
      !    courantmax=max(courantmax,maxval(tmp(ixO^S)))
      !    courantmaxtot=courantmaxtot+courantmax
      !  end do
      !end if

      select case (type_courant)
      case (type_maxsum)
        cmaxtot(ixO^S)=zero
        if(slab_uniform) then
          ^D&dxinv(^D)=one/dx^D;
          do idims=1,ndim
            call phys_get_cmax(w,x,ixI^L,ixO^L,idims,cmax)
            if(need_global_cmax) cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
            cmaxtot(ixO^S)=cmaxtot(ixO^S)+cmax(ixO^S)*dxinv(idims)
          end do
        else
          do idims=1,ndim
            call phys_get_cmax(w,x,ixI^L,ixO^L,idims,cmax)
            if(need_global_cmax) cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
            cmaxtot(ixO^S)=cmaxtot(ixO^S)+cmax(ixO^S)/block%ds(ixO^S,idims)
          end do
        end if
        ! courantmaxtots='max(summed c/dx)'
        courantmaxtots=max(courantmaxtots,maxval(cmaxtot(ixO^S)))
        if (courantmaxtots>smalldouble) dtnew=min(dtnew,courantpar/courantmaxtots)
      case (type_summax)
        if(slab_uniform) then
          ^D&dxinv(^D)=one/dx^D;
          do idims=1,ndim
            call phys_get_cmax(w,x,ixI^L,ixO^L,idims,cmax)
            if(need_global_cmax) cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
            courantmax=max(courantmax,maxval(cmax(ixO^S)*dxinv(idims)))
            courantmaxtot=courantmaxtot+courantmax
          end do
        else
          do idims=1,ndim
            call phys_get_cmax(w,x,ixI^L,ixO^L,idims,cmax)
            if(need_global_cmax) cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
            courantmax=max(courantmax,maxval(cmax(ixO^S)/block%ds(ixO^S,idims)))
            courantmaxtot=courantmaxtot+courantmax
          end do
        end if
        ! courantmaxtot='summed max(c/dx)'
        if (courantmaxtot>smalldouble)  dtnew=min(dtnew,courantpar/courantmaxtot)
      case (type_minimum)
        if(slab_uniform) then
          ^D&dxinv(^D)=one/dx^D;
          do idims=1,ndim
            call phys_get_cmax(w,x,ixI^L,ixO^L,idims,cmax)
            if(need_global_cmax) cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
            courantmax=max(courantmax,maxval(cmax(ixO^S)*dxinv(idims)))
          end do
        else
          do idims=1,ndim
            call phys_get_cmax(w,x,ixI^L,ixO^L,idims,cmax)
            if(need_global_cmax) cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
            courantmax=max(courantmax,maxval(cmax(ixO^S)/block%ds(ixO^S,idims)))
          end do
        end if
        ! courantmax='max(c/dx)'
        if (courantmax>smalldouble)     dtnew=min(dtnew,courantpar/courantmax)
      end select

    end subroutine getdt_courant

end subroutine setdt
