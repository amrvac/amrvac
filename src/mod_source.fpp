!> Module for handling split source terms (split from the fluxes)

#:mute
#:include "mod_gpu_directives.fpp"
#:include "physics/mod_physics_templates.fpp"
#:endmute

module mod_source

  use mod_variables
  use mod_physics_vars
  use mod_global_parameters, only: ndim
  
  implicit none
  public

  !> How to apply dimensional splitting to the source terms, see
  !> @ref discretization.md
  !> defaulting to sfs
  integer :: sourcesplit = 0
  ${GPU_DECLARE_COPYIN('sourcesplit')}$
  integer, parameter :: sourcesplit_sfs    = 0
  integer, parameter :: sourcesplit_sf     = 1
  integer, parameter :: sourcesplit_ssf    = 2
  integer, parameter :: sourcesplit_ssfss  = 3

  public :: add_split_source

contains

! instantiate the templated functions here for inlining:
@:addsource_local()
@:addsource_nonlocal()
@:addsource_compact()
@:get_aux()
@:to_primitive()


  subroutine add_split_source(prior)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_physics, only: phys_req_diagonal, phys_global_source_after
    use mod_comm_lib, only: mpistop

    logical, intent(in)    :: prior
    real(dp)               :: qt
    real(dp)               :: dr(ndim), xloc(ndim)
    ! Physical cell size handed to the source term. dr is the *logical*
    ! spacing, which is not a length in a curvilinear system and, under
    ! LOG_RADIUS, is not even proportional to one.
    real(dp)               :: dloc(ndim)
    integer                :: iigrid, n, ix1, ix2, ix3
    integer                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    real(dp)               :: wprim(nw_phys), wnew(nw_phys), wCT(nw_phys)
    logical                :: src_active
    
    src_active = .false.
    
    if ((.not.prior).and.(sourcesplit==sourcesplit_sf .or. &
       sourcesplit==sourcesplit_ssf)) return

    if (prior) then
       qt=global_time
    else
       qt=global_time+dt
    end if

    if(any_source_split) then
       ! add split source terms
       
       ixOmin1=ixGlo1+nghostcells;ixOmin2=ixGlo2+nghostcells
       ixOmin3=ixGlo3+nghostcells;ixOmax1=ixGhi1-nghostcells
       ixOmax2=ixGhi2-nghostcells;ixOmax3=ixGhi3-nghostcells;

      select case (sourcesplit)
      case (sourcesplit_sfs)
         ${GPU_PARALLEL_LOOP_GANG("private(n, dr)")}$
         do iigrid=1,igridstail_active
            n = igrids_active(iigrid)            
            dr  = rnode(rpdx1_:rnodehi, n)
            
            ${GPU_LOOP_VECTOR("collapse(ndim) private(xloc, dloc, wprim, wnew, wCT)")}$
            do ix3=ixOmin3,ixOmax3 
               do ix2=ixOmin2,ixOmax2 
                  do ix1=ixOmin1,ixOmax1
                     
                     xloc(1:ndim) = bgeo%x(ix1, ix2, ix3, 1:ndim, n)
#:if GEOM == 'Cartesian'
                     dloc = dr
#:else
                     dloc = bgeo%ds(ix1, ix2, ix3, 1:ndim, n)
#:endif
                     wCT   = bg(1)%w(ix1,ix2,ix3, 1:nw_phys, n)
                     wnew  = wCT
                     wprim = wCT
                     call to_primitive( wprim )

                     call addsource_local(0.5d0*dt,&
                          0.5d0, qt, wCT,&
                          wprim, qt, wnew, xloc, dloc, .true. )
                     bg(1)%w(ix1,ix2,ix3, 1:nw_phys, n) = wnew(1:nw_phys)

#:if defined('SOURCE_NONLOCAL')
                     ! TBD             
#:endif
                     
#:if defined('SOURCE_COMPACT')
                     ! TBD
#:endif                
                end do
             end do
          end do
        end do
        src_active = .true.
      case default
         write(unitterm,*)'No such type/not yet implemented sourcesplit=',sourcesplit
         call mpistop("Error: Unknown type of sourcesplit!")
      end select

      if (.not. prior .and. associated(phys_global_source_after)) then
         call phys_global_source_after(dt, qt, src_active)
      end if

      if (src_active) then
         call getbc(qt,0.d0,ps,iwstart,nwgc,phys_req_diagonal)
      end if
      
   end if 

  end subroutine add_split_source

end module mod_source
