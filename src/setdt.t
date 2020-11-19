!>setdt  - set dt for all levels between levmin and levmax. 
!>         dtpar>0  --> use fixed dtpar for all level
!>         dtpar<=0 --> determine CFL limited timestep 
subroutine setdt()
  use mod_global_parameters
  use mod_physics
  use mod_trac
  use mod_usr_methods, only: usr_get_dt
  use mod_thermal_conduction
  use mod_supertimestepping, only: set_dt_sts_ncycles, is_sts_initialized, sourcetype_sts,sourcetype_sts_split

  integer :: iigrid, igrid, ncycle, ncycle2, ifile, idim
  double precision :: dtnew, qdtnew, dtmin_mype, factor, dx^D, dxmin^D

  double precision :: dtmax, dxmin, cmax_mype, v(ixG^T)
  double precision :: a2max_mype(ndim), tco_mype, tco_global, Tmax_mype, T_peak
  double precision :: trac_alfa, trac_dmax, trac_tau

  integer, parameter :: niter_print = 2000

  if (dtpar<=zero) then
     dtmin_mype=bigdouble
     cmax_mype = zero
     a2max_mype = zero
     tco_mype = zero
     Tmax_mype = zero
  !$OMP PARALLEL DO PRIVATE(igrid,qdtnew,dtnew,dx^D)
     do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
        dtnew=bigdouble
        dx^D=rnode(rpdx^D_,igrid);
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        saveigrid = igrid
        block=>ps(igrid)
        block%iw0=0

        if (nwaux>0) then
           call phys_get_aux(.true.,ps(igrid)%w,&
                ps(igrid)%x,ixG^LL,ixM^LL,'setdt')
        end if

        call getdt_courant(ps(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,ps(igrid)%x)
        dtnew=min(dtnew,qdtnew)

        call phys_get_dt(ps(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,ps(igrid)%x)
        dtnew=min(dtnew,qdtnew)

        if (associated(usr_get_dt)) then
           call usr_get_dt(ps(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,ps(igrid)%x)
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
  endif

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

  ! estimate time step of thermal conduction
  if(associated(phys_getdt_heatconduct)) then
     dtmin_mype=bigdouble
  !$OMP PARALLEL DO PRIVATE(igrid,qdtnew,&
  !$OMP& dx^D) REDUCTION(min:dtmin_mype)
     do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
        dx^D=rnode(rpdx^D_,igrid);
        saveigrid = igrid
        block=>ps(igrid)
        qdtnew=bigdouble
        call phys_getdt_heatconduct(ps(igrid)%w,ixG^LL,ixM^LL,qdtnew,dx^D,ps(igrid)%x)
        dtmin_mype=min(dtmin_mype,qdtnew)
     end do
  !$OMP END PARALLEL DO
     call MPI_ALLREDUCE(dtmin_mype,dtnew,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                           icomm,ierrmpi)
     if(all(flux_scheme=='nul')) dt=min(dt,dtnew)
     ncycle=ceiling(dt/dtnew)
     if (ncycle>tc_ncycles) then
       if(mype==0 .and. .false.) then
         write(*,*) 'CLF time step is too many times larger than conduction time step',ncycle
         write(*,*) 'reducing dt to',tc_ncycles,'times of dt_impl!!'
       endif
       dt=tc_ncycles*dtnew
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

  if(is_sts_initialized()) then
    !!reuse qdtnew
    !qdtnew = dt 
    if(sourcetype_sts .eq. sourcetype_sts_split) then
      qdtnew = 0.5d0 * dt 
      if (set_dt_sts_ncycles(qdtnew)) then
        dt = 2d0*qdtnew
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

  !$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
     dt_grid(igrid)=dt
  end do
  !$OMP END PARALLEL DO
       
  ! global Lax-Friedrich finite difference flux splitting needs fastest wave-speed
  ! so does GLM: 
  if(need_global_cmax) call MPI_ALLREDUCE(cmax_mype,cmax_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)
  if(need_global_a2max) then 
    call MPI_ALLREDUCE(a2max_mype,a2max_global,ndim,&
      MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)
  end if

  ! transition region adaptive thermal conduction (Johnston 2019 ApJL, 873, L22)
  ! transition region adaptive thermal conduction (Johnston 2020 A&A, 635, 168)
  if(phys_trac) then
    T_bott=2.d4/unit_temperature
    call MPI_ALLREDUCE(Tmax_mype,T_peak,1,MPI_DOUBLE_PRECISION,&
         MPI_MAX,icomm,ierrmpi)
    ! default lower limit of cutoff temperature
    select case(phys_trac_type)
    case(1)
      !> 1D TRAC method
      {^IFONED
      trac_dmax=0.1d0
      trac_tau=1.d0/unit_time
      trac_alfa=trac_dmax**(dtnew/trac_tau)
      call MPI_ALLREDUCE(tco_mype,tco_global,1,MPI_DOUBLE_PRECISION,&
           MPI_MAX,icomm,ierrmpi)
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
        ps(igrid)%special_values(1)=tco_global
        if(ps(igrid)%special_values(1)<trac_alfa*ps(igrid)%special_values(2)) then
          ps(igrid)%special_values(1)=trac_alfa*ps(igrid)%special_values(2)
        end if
        if(ps(igrid)%special_values(1) < T_bott) then
          ps(igrid)%special_values(1)=T_bott
        else if(ps(igrid)%special_values(1) > 0.2d0*T_peak) then
          ps(igrid)%special_values(1)=0.2d0*T_peak
        end if
        !> special values(2) to save old tcutoff
        ps(igrid)%special_values(2)=ps(igrid)%special_values(1)
      end do
      !$OMP END PARALLEL DO
      }
      !> 2D or 3D simplified TRAC method
      call TRAC_simple(T_peak)
    case(2)
      !> 2D or 3D TRACL(ine) method
      call TRACL(.false.,T_peak)
    case(3)
      !> 2D or 3D TRACB(lock) method
      call TRACB(.false.,T_peak)
    case(4)
      !> 2D or 3D TRACL(ine) method with mask
      call TRACL(.true.,T_peak)
    case(5)
      !> 2D or 3D TRACB(lock) method with mask
      call TRACB(.true.,T_peak)
    case default
      call mpistop("undefined TRAC method type")
    end select
  end if 

  contains

    !> compute CFL limited dt (for variable time stepping)
    subroutine getdt_courant(w,ixI^L,ixO^L,dtnew,dx^D,x)
      use mod_global_parameters
      use mod_physics, only: phys_get_cmax,phys_get_a2max,phys_get_tcutoff
      
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in)    :: dx^D
      double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
      
      integer :: idims
      double precision :: courantmax, dxinv(1:ndim), courantmaxtot, courantmaxtots
      double precision :: cmax(ixI^S), cmaxtot(ixI^S), tmp(ixI^S)
      double precision :: a2max(ndim),tco_local,Tmax_local

      dtnew=bigdouble
      
      courantmax=zero
      courantmaxtot=zero
      courantmaxtots=zero
      
      ^D&dxinv(^D)=one/dx^D;
      
      cmaxtot(ixO^S)=zero
      
      if(need_global_a2max) then
        call phys_get_a2max(w,x,ixI^L,ixO^L,a2max)
      end if
      if(phys_trac) then
        call phys_get_tcutoff(ixI^L,ixO^L,w,x,tco_local,Tmax_local)
        {^IFONED tco_mype=max(tco_mype,tco_local) }
        Tmax_mype=max(Tmax_mype,Tmax_local)
      end if
      do idims=1,ndim
        call phys_get_cmax(w,x,ixI^L,ixO^L,idims,cmax)
        if(need_global_cmax) cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
        if(need_global_a2max) a2max_mype = max(a2max_mype,a2max(idims))
        if(slab_uniform) then
          cmaxtot(ixO^S)=cmaxtot(ixO^S)+cmax(ixO^S)*dxinv(idims)
          courantmax=max(courantmax,maxval(cmax(ixO^S)*dxinv(idims)))
        else
          tmp(ixO^S)=cmax(ixO^S)/block%ds(ixO^S,idims)
          cmaxtot(ixO^S)=cmaxtot(ixO^S)+tmp(ixO^S)
          courantmax=max(courantmax,maxval(tmp(ixO^S)))
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
         courantmaxtots=max(courantmaxtots,maxval(cmaxtot(ixO^S)))
         if (courantmaxtots>smalldouble) dtnew=min(dtnew,courantpar/courantmaxtots)
      case default
         write(unitterm,*)'Unknown typecourant=',typecourant
         call mpistop("Error from getdt_courant: no such typecourant!")
      end select

    end subroutine getdt_courant

end subroutine setdt
