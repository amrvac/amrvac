#:mute
#:include "physics/mod_physics_templates.fpp"
#:include "mod_gpu_directives.fpp"
#:endmute

module mod_dt

  use mod_variables
  use mod_physics_vars
  implicit none
  private

  public :: setdt

contains

! instantiate here for inlining in kernels
@:get_aux()
@:to_primitive()
@:get_cmax()
@:phys_get_dt()

  !>setdt  - set dt for all levels between levmin and levmax.
  !>         dtpar>0  --> use fixed dtpar for all level
  !>         dtpar<=0 --> determine CFL limited timestep
  subroutine setdt()
    use mod_global_parameters

    integer :: iigrid, igrid, idims, ix1,ix2,ix3, ifile
    double precision :: dtmin_mype, factor, dx1,dx2,dx3, dtmax

    double precision :: dxinv(1:ndim), cmaxtot, cmax, u(1:nw_phys), cmax_mype
    double precision :: xloc(1:ndim), qdtnew, dxloc(1:ndim)

    if (dtpar<=zero) then
       dtmin_mype=bigdouble

       ${GPU_PARALLEL_LOOP_GANG("private(igrid,dxinv,dx1,dx2,dx3) reduction(min:dtmin_mype)")}$
       do iigrid=1,igridstail_active; igrid=igrids_active(iigrid)

          dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid)
          dx3=rnode(rpdx3_,igrid);

          dxinv(1)=one/dx1;dxinv(2)=one/dx2;dxinv(3)=one/dx3;

          ${GPU_LOOP_VECTOR("collapse(ndim) reduction(min:dtmin_mype) private(u, xloc, cmaxtot, cmax, idims, qdtnew, dxloc" + (", dxinv" if GEOM != 'Cartesian' else "") + ")")}$
          do ix3=ixMlo3,ixMhi3 
             do ix2=ixMlo2,ixMhi2 
                do ix1=ixMlo1,ixMhi1 
                   u = bg(1)%w(ix1, ix2, ix3, 1:nw_phys, igrid)
                   call to_primitive(u)
                   xloc(1:ndim) = bgeo%x(ix1, ix2, ix3, 1:ndim, igrid)

#:if GEOM != 'Cartesian'
                   ! In a curvilinear system the coordinate spacing is not a
                   ! length: use the physical cell sizes set by fillgeo.
                   dxinv(1:ndim) = one / bgeo%ds(ix1, ix2, ix3, 1:ndim, igrid)
#:endif

                   cmaxtot = 0.0d0
                   ${GPU_LOOP_SEQ()}$
                   do idims = 1, ndim
                      cmax = get_cmax(u, xloc, idims)
                      cmaxtot = cmaxtot + cmax * dxinv(idims)
                   end do
                   dtmin_mype     = min( dtmin_mype, courantpar / cmaxtot )

#:if defined('SOURCE_DT')
                   u = bg(1)%w(ix1,ix2,ix3,1:nw_phys,igrid)
                   xloc(1:ndim) = bgeo%x(ix1, ix2, ix3, 1:ndim, igrid)
                   dxloc(1:ndim) = one/dxinv(1:ndim)
                   call phys_get_dt(u, xloc, dxloc, qdtnew)
                   dtmin_mype = min( dtmin_mype, qdtnew )
#:endif

                end do
             end do
          end do


       end do

    else
       dtmin_mype=dtpar
    end if

    
    if (need_global_cmax) then
       cmax_mype=-bigdouble

       ${GPU_PARALLEL_LOOP_GANG("private(igrid) reduction(max:cmax_mype)")}$
       do iigrid=1,igridstail_active; igrid=igrids_active(iigrid)

          ${GPU_LOOP_VECTOR("collapse(ndim) reduction(max:cmax_mype) private(u, xloc, cmax, idims)")}$
          do ix3=ixMlo3,ixMhi3 
             do ix2=ixMlo2,ixMhi2 
                do ix1=ixMlo1,ixMhi1

                   u(1:nw_phys) = bg(1)%w(ix1, ix2, ix3, 1:nw_phys, igrid)
                   call to_primitive(u)
                   xloc(1:ndim) = bgeo%x(ix1, ix2, ix3, 1:ndim, igrid)
                   ${GPU_LOOP_SEQ()}$
                   do idims = 1, ndim
                      cmax = get_cmax(u, xloc, idims)
                      cmax_mype = max( cmax_mype, cmax )
                   end do

                end do
             end do
          end do
       end do
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

    if (need_global_cmax) then
      call MPI_ALLREDUCE(cmax_mype, cmax_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, icomm, &
             ierrmpi)
      ${GPU_UPDATE_DEVICE('cmax_global')}$
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

    ${GPU_UPDATE_DEVICE('dt')}$
  end subroutine setdt
end module mod_dt
