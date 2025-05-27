module mod_dt

  implicit none
  private


  public :: setdt

contains


  !>setdt  - set dt for all levels between levmin and levmax. 
  !>         dtpar>0  --> use fixed dtpar for all level
  !>         dtpar<=0 --> determine CFL limited timestep 
  subroutine setdt()
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    use mod_physics, only: to_primitive, get_cmax, nw_phys

    integer :: iigrid, igrid, idims, ix1,ix2,ix3, ifile
    double precision :: dtnew, dtmin_mype, factor, dx1,dx2,dx3, dtmax

    double precision :: w(nw,ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)
    double precision :: dxinv(1:ndim), courantmaxtots, cmaxtot, cmax,&
        u(1:nw_phys)

    if (dtpar<=zero) then
       dtmin_mype=bigdouble

 !$OMP PARALLEL DO PRIVATE(igrid,dtnew,dx1,dx2,dx3) REDUCTION(min:dtmin_mype) REDUCTION(max:cmax_mype,a2max_mype)
 !$acc parallel loop PRIVATE(igrid,dx1,dx2,dx3,dxinv,w) REDUCTION(min:dtmin_mype) gang
       do iigrid=1,igridstail_active; igrid=igrids_active(iigrid)

          dtnew=bigdouble
          dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid)
          dx3=rnode(rpdx3_,igrid);

          dxinv(1)=one/dx1;dxinv(2)=one/dx2;dxinv(3)=one/dx3;
          courantmaxtots=zero

 !$acc loop vector collapse(ndim) reduction(max:courantmaxtots) private(cmax, cmaxtot, u)
           do ix3=ixMlo3,ixMhi3 
            do ix2=ixMlo2,ixMhi2 
            do ix1=ixMlo1,ixMhi1 
          w(1:nw,ix1,ix2,ix3) = bg(1)%w(ix1,ix2,ix3,1:nw,igrid)
          cmaxtot = 0.0d0
          !$acc loop seq
          do idims = 1, ndim
             u = w(:,ix1,ix2,ix3)
             call to_primitive(u)
             cmax = get_cmax(u,idims)
             cmaxtot = cmaxtot + cmax * dxinv(idims)
          end do
          courantmaxtots = max( courantmaxtots, cmaxtot )
           end do
            end do
            end do

          dtmin_mype  = min(dtmin_mype,courantpar / courantmaxtots)

       end do
       !$OMP END PARALLEL DO

    else
       dtmin_mype=dtpar
    end if

    if (dtmin_mype<dtmin) then
       write(unitterm,*)"Error: Time step too small!", dtmin_mype
       write(unitterm,*)"on processor:", mype, "at time:", global_time,&
          " step:", it
       write(unitterm,*)"Lower limit of time step:", dtmin
       crash=.true.
    end if

    if (slowsteps>it-it_init+1) then
       factor=one-(one-dble(it-it_init+1)/dble(slowsteps))**2
       dtmin_mype=dtmin_mype*factor
    end if

    if(final_dt_reduction)then
       if(time_max-global_time<=dtmin) then
          final_dt_exit=.true.
       endif
       dtmin_mype=min(dtmin_mype,time_max-global_time)
    end if

    if (dtpar<=zero) then
       call MPI_ALLREDUCE(dtmin_mype,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, icomm,&
          ierrmpi)
    else
       dt=dtmin_mype
    end if

    if(any(dtsave(1:nfile)<bigdouble).or.any(tsave(isavet(1:nfile),&
       1:nfile)<bigdouble))then
       dtmax = minval(ceiling(global_time/dtsave(1:nfile))*dtsave(1:nfile))-&
          global_time
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
          write(unitterm,*) 'Warning: timesteps: ',dt,&
             ' exceeding output intervals ', dtsave(1:nfile)
       endif
    endif

  end subroutine setdt
end module mod_dt
