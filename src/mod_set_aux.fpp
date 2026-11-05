#:mute
#:include "physics/mod_physics_templates.fpp"
#:include "mod_gpu_directives.fpp"
#:endmute

!> Refresh the auxiliary variables of the solution.
!>
!> The auxiliaries (e.g. mu and lfac for srmhd) are functions of the conserved
!> state, and the finite-volume update only writes the flux variables, so the
!> stored copies would otherwise keep their initial values for the whole run.
!> They are kept consistent because they are the initial guess of the root
!> solve in to_primitive, and because they are output.
!>
!> The rule is that whoever writes a conserved state refreshes its
!> auxiliaries, once:
!>
!>   set_aux_all   every active block, at the end of each timestep and after
!>                 the initial condition or a restart;
!>   set_aux       per cell, called by the kernels that create a block's
!>                 state during a regrid (prolong_2nd, coarsen_grid), so a new
!>                 block is never left with interpolated auxiliaries.
module mod_set_aux

  use mod_variables
  use mod_physics_vars
  implicit none
  private

  public :: set_aux_all, set_aux

contains

! instantiate here for inlining in the kernel
@:get_aux()
@:set_aux()

  !> Recompute the auxiliaries of bg(1), the solution, for every active block,
  !> ghost cells included. Does nothing for a physics without auxiliaries.
  !>
  !> The w array is laid out as the flux variables 1:nwflux, then the nwaux
  !> auxiliaries, then the extra variables (var_set_extravar must come after
  !> var_set_auxvar). Only the auxiliary slots nwflux+1:nwflux+nwaux are
  !> written back; the extras, which belong to other machinery, are not touched.
  !>
  !> The physics' set_aux decides what the auxiliaries are. A physics that has
  !> some but does not define set_aux gets the empty default, which leaves them
  !> as they are.
  !>
  !> Call it wherever the conserved variables of every block were changed by
  !> something other than the finite-volume update. Ghost cells are covered so
  !> that no ghost exchange is needed afterwards, which requires the ghost
  !> cells of the conserved variables to be filled already.
  subroutine set_aux_all()
    use mod_global_parameters

    integer          :: iigrid, igrid, ix1, ix2, ix3, iw
    double precision :: u(1:nw_phys)

    if (nwaux == 0) return

    ${GPU_PARALLEL_LOOP_GANG("private(igrid)")}$
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid)

       ${GPU_LOOP_VECTOR("collapse(ndim) private(u, iw)")}$
       do ix3=ixGlo3,ixGhi3
          do ix2=ixGlo2,ixGhi2
             do ix1=ixGlo1,ixGhi1
                u = bg(1)%w(ix1, ix2, ix3, 1:nw_phys, igrid)
                call set_aux(u)
                ${GPU_LOOP_SEQ()}$
                do iw = nwflux+1, nwflux+nwaux
                   bg(1)%w(ix1, ix2, ix3, iw, igrid) = u(iw)
                end do
             end do
          end do
       end do
    end do

  end subroutine set_aux_all

end module mod_set_aux
