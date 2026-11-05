!> Module with finite volume methods for fluxes

#:mute
#:include "physics/mod_physics_templates.fpp"
#:include "mod_gpu_directives.fpp"
#:endmute

module mod_finite_volume

  use mod_variables
  use mod_physics_vars
  use mod_global_parameters, only: ndim
#:if defined('FLUX_AD')
  use mod_global_parameters, only: flux_ad_min, flux_ad_scale
#:endif
  use mod_physicaldata
  implicit none

  private

  public :: finite_volume_local

contains

! instantiate the templated functions here for inlining:
@:addsource_local()
@:addsource_nonlocal()
@:addsource_compact()
#:if GEOM != 'Cartesian'
@:addsource_geometry()
#:endif
@:get_aux()
@:to_primitive()
@:to_conservative()
@:fix_prim_state()
@:get_cmax()
@:get_flux()
@:estimate_speeds_minmax()
@:estimate_speeds_hllc()
@:hllc_flux()

  ! flux scheme list : (scheme_tag, method_enum, faceflux_proc)
#:set schemes = [ &
  ('muscl_llf',  'fs_tvdlf', 'reconflux_muscl_llf_prim'), &
  ('muscl_hll',  'fs_hll',   'reconflux_muscl_hll_prim'), &
  ('muscl_hllc', 'fs_hllc',  'reconflux_muscl_hllc_prim')]

  subroutine finite_volume_local(qdt, dtfactor, ixImin1,ixImin2,&
    ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
    ixOmax3, idimsmin,idimsmax, qtC, bga, qt, bgb)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    double precision, intent(in)                                       :: qdt,&
      dtfactor, qtC, qt
    integer, intent(in)                                                :: &
      ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
      ixOmin3,ixOmax1,ixOmax2,ixOmax3, idimsmin,idimsmax
    ! remember, old names map as: wCT => bga, wnew => bgb
    type(block_grid_t)                                    :: bga
    type(block_grid_t)                                    :: bgb
    ! .. local ..
    integer                :: n, iigrid, ix1,ix2,ix3
    double precision       :: uprim(nw_phys, ixImin1:ixImax1,ixImin2:ixImax2,&
      ixImin3:ixImax3)
    real(dp)               :: tmp(nw_phys,5)
    real(dp)               :: f(nw_flux, 2)
    real(dp)               :: inv_dr(ndim)
    real(dp)               :: dr(ndim)
    integer                :: typelim
    real(dp)               :: xloc(ndim)
    real(dp)               :: xlocC(ndim,2)
    real(dp)               :: wprim(nw_phys), wCT(nw_phys), wnew(nw_phys)
    integer :: method
    !-----------------------------------------------------------------------------

    method = flux_method(1)  ! TODO: implement per grid level
    select case (method)
#:for scheme_tag, method_enum, faceflux_proc in schemes
    case (${method_enum}$)
      call finite_volume_local_${scheme_tag}$(qdt, dtfactor, ixImin1,ixImin2,&
            ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3, idimsmin,idimsmax, qtC, bga, qt, bgb)
#:endfor
    case default
      call mpistop("finite_volume_local: unknown flux scheme")
    end select
end subroutine finite_volume_local

#:def FV_KERNEL(scheme_tag, faceflux_proc)
  subroutine finite_volume_local_${scheme_tag}$(qdt, dtfactor, ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, idimsmin,idimsmax, qtC, bga, qt, bgb)
    use mod_global_parameters
    use mod_fix_conserve

    double precision, intent(in)                                       :: qdt,&
        dtfactor, qtC, qt
    integer, intent(in)                                                :: &
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3, idimsmin,idimsmax
    ! remember, old names map as: wCT => bga, wnew => bgb
    type(block_grid_t)                                    :: bga
    type(block_grid_t)                                    :: bgb
    ! .. local ..
    integer                :: n, iigrid, ix1,ix2,ix3
    integer                :: nx1,nx2,nx3, nxCo1,nxCo2,nxCo3
    double precision       :: uprim(nw_phys, ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    real(dp)               :: tmp(nw_phys,5)
    real(dp)               :: tmp1(nw_phys,3),tmp2(nw_phys,3),tmp3(nw_phys,3)
    real(dp)               :: f(nw_flux, 2)
    real(dp)               :: inv_dr(ndim)
    real(dp)               :: dr(ndim)
    integer                :: typelim
    real(dp)               :: xloc(ndim)
    real(dp)               :: xlocC(ndim,2)
    real(dp)               :: wprim(nw_phys), wCT(nw_phys), wnew(nw_phys)
    ! Block corner in the logical coordinate, for building face positions.
    real(dp)               :: xlo(ndim)
#:if GEOM != 'Cartesian'
    ! Inverse cell volume, and the (upper minus lower) face area over volume
    ! per direction, which the curvilinear geometric source terms need.
    real(dp)               :: inv_dvol, dAdV(ndim)
#:endif
#:if defined('SOURCE_LOCAL') or defined('SOURCE_NONLOCAL') or defined('SOURCE_COMPACT')
    ! Physical cell size handed to the optional source terms.  dr below is the
    ! *logical* spacing, which is not a length in a curvilinear system and is
    ! not even proportional to one under LOG_RADIUS.
    real(dp)               :: dloc(ndim)
#:endif
    real(dp)       :: fC1(2,ixImin2:ixOmax2,ixImin3:ixOmax3,1:nwflux)
    real(dp)       :: fC2(2,ixImin1:ixOmax1,ixImin3:ixOmax3,1:nwflux)
    real(dp)       :: fC3(2,ixImin1:ixOmax1,ixImin2:ixOmax2,1:nwflux)
    integer, parameter     :: max_batch=4096
    integer                :: nbatches, ibatch, igrid_beg, igrid_end
    integer        :: neighbor_type_1m, neighbor_type_1p, neighbor_type_2m,&
              neighbor_type_2p, neighbor_type_3m, neighbor_type_3p
    !-----------------------------------------------------------------------------

    nx1=ixMhi1-ixMlo1+1
    nx2=ixMhi2-ixMlo2+1
    nx3=ixMhi3-ixMlo3+1

    nxCo1=nx1/2
    nxCo2=nx2/2
    nxCo3=nx3/2

    ${GPU_ENTER_DATA_COPYIN('nxCo1,nxCo2,nxCo3')}$
    ! batch-launch the kernels (oom issue when many ~30000 blocks):
    nbatches = (igridstail_active + max_batch - 1) / max_batch ! ceiling division
    
    do ibatch = 1, nbatches
       igrid_beg = (ibatch-1) * max_batch + 1
       igrid_end = min(ibatch * max_batch, igridstail_active)

       ${GPU_PARALLEL_LOOP_GANG("private(uprim, inv_dr, dr, n, xlo, typelim, ix1, ix2, ix3, fC1, fC2, fC3, neighbor_type_1m, neighbor_type_1p, neighbor_type_2m, neighbor_type_2p, neighbor_type_3m, neighbor_type_3p)")}$
       do iigrid = igrid_beg, igrid_end
          n = igrids_active(iigrid)

          dr  = rnode(rpdx1_:rnodehi, n)
          inv_dr  = 1/dr
          xlo = rnode(rpxmin1_:rpxmin1_+ndim-1, n)
          typelim = type_limiter(node(plevel_, n))

          !JESSENEW These were called many times in the 3D loop below
          !but they technically only have an igrid dependence so I
          !thought it would be better to simply declare them via local
          !variable which will reduce the call overhead ~NX^3
          neighbor_type_1m = neighbor_type(-1,0,0,n)
          neighbor_type_1p = neighbor_type( 1,0,0,n)
          neighbor_type_2m = neighbor_type(0,-1,0,n)
          neighbor_type_2p = neighbor_type(0, 1,0,n)
          neighbor_type_3m = neighbor_type(0,0,-1,n)
          neighbor_type_3p = neighbor_type(0,0, 1,n)

          ${GPU_LOOP_VECTOR("collapse(ndim)")}$
          do ix3=ixImin3,ixImax3 
             do ix2=ixImin2,ixImax2 
                do ix1=ixImin1,ixImax1 
                   ! Convert to primitive
                   uprim(1:nw_phys, ix1,ix2,ix3) = bga%w(ix1,ix2,ix3, 1:nw_phys, n)
                   call to_primitive(uprim(1:nw_phys, ix1,ix2,ix3))
                end do
             end do
          end do

       ${GPU_LOOP_VECTOR("collapse(ndim) private(f, wnew, tmp, xlocC, xloc" + (", wCT" if defined('SOURCE_LOCAL') else "") + (", wprim" if defined('SOURCE_LOCAL') or GEOM != 'Cartesian' else "") + (", inv_dvol, dAdV" if GEOM != 'Cartesian' else "") + (", tmp1,tmp2,tmp3" if defined('SOURCE_COMPACT') else "") + (", dloc" if defined('SOURCE_LOCAL') or defined('SOURCE_NONLOCAL') or defined('SOURCE_COMPACT') else "") + ")")}$
       do ix3=ixOmin3,ixOmax3 
          do ix2=ixOmin2,ixOmax2 
             do ix1=ixOmin1,ixOmax1 
                ! Compute fluxes in all dimensions

#:if GEOM != 'Cartesian'
                inv_dvol = 1.0_dp / bgeo%dvolume(ix1, ix2, ix3, n)
#:endif

                ! This cell's stored position - its volume barycentre - read
                ! once and reused by all three face constructions and by every
                ! source term below. The component index of bgeo%x sits between
                ! the cell indices and the grid index, so each of these is a
                ! strided gather; fetching it once rather than per direction
                ! takes ten of them per cell down to one.
                xloc(1:ndim) = bgeo%x(ix1, ix2, ix3, 1:ndim, n)

                tmp = uprim(1:nw_phys, ix1-2:ix1+2, ix2, ix3)
                ! transverse components stay on the barycentre; only the normal
                ! one is moved to the interface
                xlocC(1:ndim,1) = xloc(1:ndim)
                xlocC(1:ndim,2) = xloc(1:ndim)
                ! Face positions, built from the block corner and the cell
                ! index.  This is the one form that is right in every geometry,
                ! so it is used unconditionally.  It cannot be shortened to
                ! x -/+ dr/2: bgeo%x is the cell's volume *barycentre* in a
                ! curvilinear build, and under LOG_RADIUS dr(1) is a ratio in
                ! ln(r) rather than a length, so that form would not merely be
                ! imprecise but wrong.  Where x is the midpoint - Cartesian in
                ! every direction, and phi and cylindrical z besides - the two
                ! agree to round-off.
                !
                ! The lower-corner form suffices because this loop covers the
                ! mesh interior only; fill_geometry_device's two-sided rule
                ! exists to make neighbouring blocks' *ghost* cells agree
                ! bit-for-bit, which nothing here needs.
                xlocC(1,1) = xlo(1) + dble(ix1-nghostcells-1)*dr(1)
                xlocC(1,2) = xlocC(1,1) + dr(1)
#:if defined('LOG_RADIUS')
                ! the logical radial coordinate is ln(1 + r/r0); this is
                ! r_of_s, in the reduced form that holds wherever s >= 0.
                ! That covers this loop, which runs over the mesh interior
                ! only - the faces it builds are domain faces, and s = 0 is
                ! r = 0 - so the odd branch r_of_s takes in the ghost layer
                ! beyond a cylindrical axis is not needed here.
                xlocC(1,1) = log_ra*dexp(xlocC(1,1)) + log_rb
                xlocC(1,2) = log_ra*dexp(xlocC(1,2)) + log_rb
#:endif
                call ${faceflux_proc}$(tmp, xlocC, 1, f, typelim)
#:if GEOM == 'Cartesian'
                bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = bgb%w(ix1, ix2, ix3, 1:nw_flux,&
                     n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(1)
#:else
                bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = bgb%w(ix1, ix2, ix3, 1:nw_flux,&
                     n) + qdt * (f(:, 1) * bgeo%surfaceC(ix1-1, ix2, ix3, 1, n) &
                     - f(:, 2) * bgeo%surfaceC(ix1, ix2, ix3, 1, n)) * inv_dvol
#:endif

                ! Store fluxes for flux fixing in direction 1.  What is
                ! stored differs by geometry: qdt*f/dx, the face's own
                ! contribution to the cell update, on a Cartesian mesh, and the
                ! extensive qdt*f*A on a curvilinear one.  The extensive form is
                ! what lets fix_conserve divide by the *coarse* cell's volume
                ! and sum the four fine faces without a volume ratio, since the
                ! fine face areas add up to the coarse one exactly.
                select case (neighbor_type_1m)
                case (neighbor_fine)
#:if GEOM == 'Cartesian'
                   if (ix1.eq.ixOmin1) pflux(1,1)%flux(1,ix2-nghostcells,ix3-nghostcells,1:nw_flux,n) &
                        = qdt * inv_dr(1) * f(:,1)
#:else
                   if (ix1.eq.ixOmin1) pflux(1,1)%flux(1,ix2-nghostcells,ix3-nghostcells,1:nw_flux,n) &
                        = qdt * bgeo%surfaceC(ix1-1, ix2, ix3, 1, n) * f(:,1)
#:endif
                case (neighbor_coarse)
#:if GEOM == 'Cartesian'
                   if (ix1.eq.ixOmin1) fC1(1,ix2,ix3,1:nw_flux) = - qdt * inv_dr(1) * f(:,1)
#:else
                   if (ix1.eq.ixOmin1) fC1(1,ix2,ix3,1:nw_flux) = - qdt * bgeo%surfaceC(ix1-1, ix2, ix3, 1, n) * f(:,1)
#:endif
                end select

                select case (neighbor_type_1p)
                case (neighbor_fine)
#:if GEOM == 'Cartesian'
                   if (ix1.eq.ixOmax1) pflux(2,1)%flux(1,ix2-nghostcells,ix3-nghostcells,1:nw_flux,n) &
                        = - qdt * inv_dr(1) * f(:,2)
#:else
                   if (ix1.eq.ixOmax1) pflux(2,1)%flux(1,ix2-nghostcells,ix3-nghostcells,1:nw_flux,n) &
                        = - qdt * bgeo%surfaceC(ix1, ix2, ix3, 1, n) * f(:,2)
#:endif
                case (neighbor_coarse)
#:if GEOM == 'Cartesian'
                   if (ix1.eq.ixOmax1) fC1(2,ix2,ix3,1:nw_flux) = qdt * inv_dr(1) * f(:,2)
#:else
                   if (ix1.eq.ixOmax1) fC1(2,ix2,ix3,1:nw_flux) = qdt * bgeo%surfaceC(ix1, ix2, ix3, 1, n) * f(:,2)
#:endif
                end select

                tmp = uprim(1:nw_phys, ix1, ix2-2:ix2+2, ix3)
                xlocC(1:ndim,1) = xloc(1:ndim)
                xlocC(1:ndim,2) = xloc(1:ndim)
                xlocC(2,1) = xlo(2) + dble(ix2-nghostcells-1)*dr(2)
                xlocC(2,2) = xlocC(2,1) + dr(2)
                call ${faceflux_proc}$(tmp, xlocC, 2, f, typelim)
#:if GEOM == 'Cartesian'
                bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = bgb%w(ix1, ix2, ix3, 1:nw_flux,&
                     n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(2)
#:else
                bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = bgb%w(ix1, ix2, ix3, 1:nw_flux,&
                     n) + qdt * (f(:, 1) * bgeo%surfaceC(ix1, ix2-1, ix3, 2, n) &
                     - f(:, 2) * bgeo%surfaceC(ix1, ix2, ix3, 2, n)) * inv_dvol
#:endif

                ! Store fluxes for flux fixing in direction 2.  What is
                ! stored differs by geometry: qdt*f/dx, the face's own
                ! contribution to the cell update, on a Cartesian mesh, and the
                ! extensive qdt*f*A on a curvilinear one.  The extensive form is
                ! what lets fix_conserve divide by the *coarse* cell's volume
                ! and sum the four fine faces without a volume ratio, since the
                ! fine face areas add up to the coarse one exactly.
                select case (neighbor_type_2m)
                case (neighbor_fine)
#:if GEOM == 'Cartesian'
                   if (ix2.eq.ixOmin2) pflux(1,2)%flux(ix1-nghostcells,1,ix3-nghostcells,1:nw_flux,n) &
                        = qdt * inv_dr(2) * f(:,1)
#:else
                   if (ix2.eq.ixOmin2) pflux(1,2)%flux(ix1-nghostcells,1,ix3-nghostcells,1:nw_flux,n) &
                        = qdt * bgeo%surfaceC(ix1, ix2-1, ix3, 2, n) * f(:,1)
#:endif
                case (neighbor_coarse)
#:if GEOM == 'Cartesian'
                   if (ix2.eq.ixOmin2) fC2(1,ix1,ix3,1:nw_flux) = - qdt * inv_dr(2) * f(:,1)
#:else
                   if (ix2.eq.ixOmin2) fC2(1,ix1,ix3,1:nw_flux) = - qdt * bgeo%surfaceC(ix1, ix2-1, ix3, 2, n) * f(:,1)
#:endif
                end select

                select case (neighbor_type_2p)
                case (neighbor_fine)
#:if GEOM == 'Cartesian'
                   if (ix2.eq.ixOmax2) pflux(2,2)%flux(ix1-nghostcells,1,ix3-nghostcells,1:nw_flux,n) &
                        = - qdt * inv_dr(2) * f(:,2)
#:else
                   if (ix2.eq.ixOmax2) pflux(2,2)%flux(ix1-nghostcells,1,ix3-nghostcells,1:nw_flux,n) &
                        = - qdt * bgeo%surfaceC(ix1, ix2, ix3, 2, n) * f(:,2)
#:endif
                case (neighbor_coarse)
#:if GEOM == 'Cartesian'
                   if (ix2.eq.ixOmax2) fC2(2,ix1,ix3,1:nw_flux) = qdt * inv_dr(2) * f(:,2)
#:else
                   if (ix2.eq.ixOmax2) fC2(2,ix1,ix3,1:nw_flux) = qdt * bgeo%surfaceC(ix1, ix2, ix3, 2, n) * f(:,2)
#:endif
                end select

                tmp = uprim(1:nw_phys, ix1, ix2, ix3-2:ix3+2)
                xlocC(1:ndim,1) = xloc(1:ndim)
                xlocC(1:ndim,2) = xloc(1:ndim)
                xlocC(3,1) = xlo(3) + dble(ix3-nghostcells-1)*dr(3)
                xlocC(3,2) = xlocC(3,1) + dr(3)
                call ${faceflux_proc}$(tmp, xlocC, 3, f, typelim)
#:if GEOM == 'Cartesian'
                bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = bgb%w(ix1, ix2, ix3, 1:nw_flux,&
                     n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(3)
#:else
                bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = bgb%w(ix1, ix2, ix3, 1:nw_flux,&
                     n) + qdt * (f(:, 1) * bgeo%surfaceC(ix1, ix2, ix3-1, 3, n) &
                     - f(:, 2) * bgeo%surfaceC(ix1, ix2, ix3, 3, n)) * inv_dvol
#:endif

                ! Store fluxes for flux fixing in direction 3.  What is
                ! stored differs by geometry: qdt*f/dx, the face's own
                ! contribution to the cell update, on a Cartesian mesh, and the
                ! extensive qdt*f*A on a curvilinear one.  The extensive form is
                ! what lets fix_conserve divide by the *coarse* cell's volume
                ! and sum the four fine faces without a volume ratio, since the
                ! fine face areas add up to the coarse one exactly.
                select case (neighbor_type_3m)
                case (neighbor_fine)
#:if GEOM == 'Cartesian'
                   if (ix3.eq.ixOmin3) pflux(1,3)%flux(ix1-nghostcells,ix2-nghostcells,1,1:nw_flux,n) &
                        = qdt * inv_dr(3) * f(:,1)
#:else
                   if (ix3.eq.ixOmin3) pflux(1,3)%flux(ix1-nghostcells,ix2-nghostcells,1,1:nw_flux,n) &
                        = qdt * bgeo%surfaceC(ix1, ix2, ix3-1, 3, n) * f(:,1)
#:endif
                case (neighbor_coarse)
#:if GEOM == 'Cartesian'
                   if (ix3.eq.ixOmin3) fC3(1,ix1,ix2,1:nw_flux) = - qdt * inv_dr(3) * f(:,1)
#:else
                   if (ix3.eq.ixOmin3) fC3(1,ix1,ix2,1:nw_flux) = - qdt * bgeo%surfaceC(ix1, ix2, ix3-1, 3, n) * f(:,1)
#:endif
                end select

                select case (neighbor_type_3p)
                case (neighbor_fine)
#:if GEOM == 'Cartesian'
                   if (ix3.eq.ixOmax3) pflux(2,3)%flux(ix1-nghostcells,ix2-nghostcells,1,1:nw_flux,n) &
                        = - qdt * inv_dr(3) * f(:,2)
#:else
                   if (ix3.eq.ixOmax3) pflux(2,3)%flux(ix1-nghostcells,ix2-nghostcells,1,1:nw_flux,n) &
                        = - qdt * bgeo%surfaceC(ix1, ix2, ix3, 3, n) * f(:,2)
#:endif
                case (neighbor_coarse)
#:if GEOM == 'Cartesian'
                   if (ix3.eq.ixOmax3) fC3(2,ix1,ix2,1:nw_flux) = qdt * inv_dr(3) * f(:,2)
#:else
                   if (ix3.eq.ixOmax3) fC3(2,ix1,ix2,1:nw_flux) = qdt * bgeo%surfaceC(ix1, ix2, ix3, 3, n) * f(:,2)
#:endif
                end select

#:if GEOM != 'Cartesian'
                   ! Add the geometric (curvature) source terms that the
                   ! flux-divergence form of the momentum equations leaves over
                   ! in a curvilinear coordinate system.
                   wprim        = uprim(1:nw_phys, ix1, ix2, ix3)
                   wnew         = bgb%w(ix1, ix2, ix3, 1:nw_phys, n)
                   dAdV(1) = (bgeo%surfaceC(ix1, ix2, ix3, 1, n) &
                        - bgeo%surfaceC(ix1-1, ix2, ix3, 1, n)) * inv_dvol
                   dAdV(2) = (bgeo%surfaceC(ix1, ix2, ix3, 2, n) &
                        - bgeo%surfaceC(ix1, ix2-1, ix3, 2, n)) * inv_dvol
                   dAdV(3) = (bgeo%surfaceC(ix1, ix2, ix3, 3, n) &
                        - bgeo%surfaceC(ix1, ix2, ix3-1, 3, n)) * inv_dvol
                   call addsource_geometry(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                        wprim, wnew, xloc, dAdV)
                   bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = wnew(1:nw_flux)
#:endif

#:if defined('SOURCE_LOCAL') or defined('SOURCE_NONLOCAL') or defined('SOURCE_COMPACT')
                   ! The optional source terms want a physical cell size. In a
                   ! curvilinear build that is bgeo%ds, not the logical spacing.
  #:if GEOM == 'Cartesian'
                   dloc = dr
  #:else
                   dloc = bgeo%ds(ix1, ix2, ix3, 1:ndim, n)
  #:endif
#:endif

#:if defined('SOURCE_LOCAL')
                   ! Add local source terms:
                   wprim        = uprim(1:nw_phys, ix1, ix2, ix3)
                   wCT          = bga%w(ix1, ix2, ix3, 1:nw_phys, n)
                   wnew         = bgb%w(ix1, ix2, ix3, 1:nw_phys, n)
                   call addsource_local(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                        dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, wCT,&
                        wprim, qt, wnew, xloc, dloc, .false. )
                   bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = wnew(1:nw_flux)
#:endif

#:if defined('SOURCE_NONLOCAL')
                   ! Add non-local (gradient) source terms:
                   wnew         = bgb%w(ix1, ix2, ix3, 1:nw_phys, n)

                   tmp = uprim(1:nw_phys, ix1-2:ix1+2, ix2, ix3)
                   call addsource_nonlocal(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                        dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, tmp,&
                        qt, wnew, xloc, dloc, 1, .false. )

                   tmp = uprim(1:nw_phys, ix1, ix2-2:ix2+2, ix3)
                   call addsource_nonlocal(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                        dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, tmp,&
                        qt, wnew, xloc, dloc, 2, .false. )

                   tmp = uprim(1:nw_phys, ix1, ix2, ix3-2:ix3+2)
                   call addsource_nonlocal(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                        dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, tmp,&
                        qt, wnew, xloc, dloc, 3, .false. )

                   bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = wnew(1:nw_flux)           
#:endif

#:if defined('SOURCE_COMPACT')
                   ! Add non-local compact source terms:
                   wnew         = bgb%w(ix1, ix2, ix3, 1:nw_phys, n)
                   tmp1 = uprim(1:nw_phys, ix1-1:ix1+1, ix2, ix3)
                   tmp2 = uprim(1:nw_phys, ix1, ix2-1:ix2+1, ix3)
                   tmp3 = uprim(1:nw_phys, ix1, ix2, ix3-1:ix3+1)
                   call addsource_compact(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                        dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, tmp1,tmp2,tmp3, &
                        qt, wnew, xloc, dloc, .false. )
                   bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = wnew(1:nw_flux)
#:endif

                end do
             end do
          end do

          ! Reduce fluxes to be stored for the coarse neighbor case
          ! Direction 1
          select case (neighbor_type_1m)
              case (neighbor_coarse)
                  ${GPU_LOOP_VECTOR('collapse(ndim-1)')}$
                  do ix3=1,nxCo3 
                    do ix2=1,nxCo2 
                  pflux(1,1)%flux(1,ix2,ix3,1:nw_flux,n) = &
                    fC1(1,1+2*(ix2-1)+nghostcells,1+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC1(1,2+2*(ix2-1)+nghostcells,1+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC1(1,1+2*(ix2-1)+nghostcells,2+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC1(1,2+2*(ix2-1)+nghostcells,2+2*(ix3-1)+nghostcells,1:nw_flux) 
                    end do
                  end do
          end select

          select case (neighbor_type_1p)
              case (neighbor_coarse)
                  ${GPU_LOOP_VECTOR('collapse(ndim-1)')}$
                  do ix3=1,nxCo3 
                     do ix2=1,nxCo2 
                  pflux(2,1)%flux(1,ix2,ix3,1:nw_flux,n) = &
                    fC1(2,1+2*(ix2-1)+nghostcells,1+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC1(2,2+2*(ix2-1)+nghostcells,1+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC1(2,1+2*(ix2-1)+nghostcells,2+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC1(2,2+2*(ix2-1)+nghostcells,2+2*(ix3-1)+nghostcells,1:nw_flux) 
                     end do
                  end do
          end select

          ! Direction 2
          select case (neighbor_type_2m)
              case (neighbor_coarse)
                  ${GPU_LOOP_VECTOR('collapse(ndim-1)')}$
                  do ix3=1,nxCo3 
                    do ix1=1,nxCo1 
                  pflux(1,2)%flux(ix1,1,ix3,1:nw_flux,n) = &
                    fC2(1,1+2*(ix1-1)+nghostcells,1+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC2(1,2+2*(ix1-1)+nghostcells,1+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC2(1,1+2*(ix1-1)+nghostcells,2+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC2(1,2+2*(ix1-1)+nghostcells,2+2*(ix3-1)+nghostcells,1:nw_flux) 
                    end do
                  end do
          end select

          select case (neighbor_type_2p)
              case (neighbor_coarse)
                  ${GPU_LOOP_VECTOR('collapse(ndim-1)')}$
                  do ix3=1,nxCo3 
                     do ix1=1,nxCo1 
                  pflux(2,2)%flux(ix1,1,ix3,1:nw_flux,n) = &
                    fC2(2,1+2*(ix1-1)+nghostcells,1+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC2(2,2+2*(ix1-1)+nghostcells,1+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC2(2,1+2*(ix1-1)+nghostcells,2+2*(ix3-1)+nghostcells,1:nw_flux) &
                   +fC2(2,2+2*(ix1-1)+nghostcells,2+2*(ix3-1)+nghostcells,1:nw_flux) 
                     end do
                  end do
          end select
          
          ! Direction 3
          select case (neighbor_type_3m)
              case (neighbor_coarse)
                  ${GPU_LOOP_VECTOR('collapse(ndim-1)')}$
                  do ix2=1,nxCo2 
                    do ix1=1,nxCo1 
                  pflux(1,3)%flux(ix1,ix2,1,1:nw_flux,n) = &
                    fC3(1,1+2*(ix1-1)+nghostcells,1+2*(ix2-1)+nghostcells,1:nw_flux) &
                   +fC3(1,2+2*(ix1-1)+nghostcells,1+2*(ix2-1)+nghostcells,1:nw_flux) &
                   +fC3(1,1+2*(ix1-1)+nghostcells,2+2*(ix2-1)+nghostcells,1:nw_flux) &
                   +fC3(1,2+2*(ix1-1)+nghostcells,2+2*(ix2-1)+nghostcells,1:nw_flux) 
                    end do
                  end do
          end select

          select case (neighbor_type_3p)
              case (neighbor_coarse)
                  ${GPU_LOOP_VECTOR('collapse(ndim-1)')}$
                  do ix2=1,nxCo2 
                     do ix1=1,nxCo1 
                  pflux(2,3)%flux(ix1,ix2,1,1:nw_flux,n) = &
                    fC3(2,1+2*(ix1-1)+nghostcells,1+2*(ix2-1)+nghostcells,1:nw_flux) &
                   +fC3(2,2+2*(ix1-1)+nghostcells,1+2*(ix2-1)+nghostcells,1:nw_flux) &
                   +fC3(2,1+2*(ix1-1)+nghostcells,2+2*(ix2-1)+nghostcells,1:nw_flux) &
                   +fC3(2,2+2*(ix1-1)+nghostcells,2+2*(ix2-1)+nghostcells,1:nw_flux) 
                     end do
                  end do
          end select

       end do
    end do

  end subroutine finite_volume_local_${scheme_tag}$
  #:enddef


#:for scheme_tag, method_enum, faceflux_proc in schemes
#:call FV_KERNEL(scheme_tag, faceflux_proc)
#:endcall
#:endfor


  !> MUSCL reconstruction in primitive variables for two faces using a 5-point stencil.
  !> Returns uL(:,iface), uR(:,iface) for iface=1 (between cells 2-3) and iface=2 (between 3-4).
  pure subroutine muscl_reconstruct_prim(u, typelim, uL, uR)
    use mod_limiter, only: limiter_minmod, limiter_vanleer, limiter_mcbeta, limiter_koren
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys,5)
    integer,  intent(in)  :: typelim
    real(dp), intent(out) :: uL(nw_phys,2), uR(nw_phys,2)

    integer  :: iw

    ! Face 1 lies between cells 2 and 3, face 2 between cells 3 and 4. Each
    ! face state is the cell value plus half its limited slope towards the face.
    select case (typelim)
    case (limiter_minmod)
      do iw=1,nw_recon
        uL(iw,1) = u(iw,2) + 0.5_dp*minmod(u(iw,2)-u(iw,1), u(iw,3)-u(iw,2))
        uR(iw,1) = u(iw,3) - 0.5_dp*minmod(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
        uL(iw,2) = u(iw,3) + 0.5_dp*minmod(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
        uR(iw,2) = u(iw,4) - 0.5_dp*minmod(u(iw,4)-u(iw,3), u(iw,5)-u(iw,4))
      end do
    case (limiter_vanleer)
      do iw=1,nw_recon
        uL(iw,1) = u(iw,2) + 0.5_dp*vanleer(u(iw,2)-u(iw,1), u(iw,3)-u(iw,2))
        uR(iw,1) = u(iw,3) - 0.5_dp*vanleer(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
        uL(iw,2) = u(iw,3) + 0.5_dp*vanleer(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
        uR(iw,2) = u(iw,4) - 0.5_dp*vanleer(u(iw,4)-u(iw,3), u(iw,5)-u(iw,4))
      end do
    case (limiter_mcbeta)
      do iw=1,nw_recon
        uL(iw,1) = u(iw,2) + 0.5_dp*mcbeta(u(iw,2)-u(iw,1), u(iw,3)-u(iw,2))
        uR(iw,1) = u(iw,3) - 0.5_dp*mcbeta(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
        uL(iw,2) = u(iw,3) + 0.5_dp*mcbeta(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
        uR(iw,2) = u(iw,4) - 0.5_dp*mcbeta(u(iw,4)-u(iw,3), u(iw,5)-u(iw,4))
      end do
    case (limiter_koren)
      ! koren(a,b) = phi(a/b)*b is one-sided: towards a cell's right face it
      ! takes (dw+, dw-), towards its left face (dw-, dw+). Swapping them
      ! breaks mirror symmetry of the scheme.
      do iw=1,nw_recon
        uL(iw,1) = u(iw,2) + 0.5_dp*koren(u(iw,3)-u(iw,2), u(iw,2)-u(iw,1))
        uR(iw,1) = u(iw,3) - 0.5_dp*koren(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
        uL(iw,2) = u(iw,3) + 0.5_dp*koren(u(iw,4)-u(iw,3), u(iw,3)-u(iw,2))
        uR(iw,2) = u(iw,4) - 0.5_dp*koren(u(iw,4)-u(iw,3), u(iw,5)-u(iw,4))
      end do
    case default  ! Fallback: Godunov (piecewise constant)
      uL(:,1) = u(:,2)
      uR(:,1) = u(:,3)
      uL(:,2) = u(:,3)
      uR(:,2) = u(:,4)
    end select

    ! The slopes were limited slot by slot, so a state can come out
    ! inadmissible or -- where the physics carries auxiliary variables --
    ! internally inconsistent. Repair it here, once, rather than in each of the
    ! three Riemann solvers: this is the single producer of face states.
    call fix_prim_state(uL(:,1))
    call fix_prim_state(uR(:,1))
    call fix_prim_state(uL(:,2))
    call fix_prim_state(uR(:,2))
  end subroutine muscl_reconstruct_prim


  !> One-face LLF/Rusanov numerical flux from primitive L/R states.
  !> phi(nw_flux) is the adaptive-diffusion reduction factor; only present when FLUX_AD is defined.
  subroutine riemann_llf_prim(uL, uR, xC, flux_dim, F#{if defined('FLUX_AD')}#, phi#{endif}#)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: uL(nw_phys), uR(nw_phys)
    real(dp), intent(in)    :: xC(ndim)
    integer,  intent(in)    :: flux_dim
    real(dp), intent(out)   :: F(nw_flux)
#:if defined('FLUX_AD')
    real(dp), intent(in)    :: phi(nw_flux)
#:endif

    real(dp) :: flux_l(nw_flux), flux_r(nw_flux)
    real(dp) :: wL, wR, wmax

    call get_flux(uL, xC, flux_dim, flux_l)
    call get_flux(uR, xC, flux_dim, flux_r)

    wL   = get_cmax(uL, xC, flux_dim)
    wR   = get_cmax(uR, xC, flux_dim)
    wmax = max(abs(wL), abs(wR))

    call to_conservative(uL)
    call to_conservative(uR)

    F = 0.5_dp * ((flux_l + flux_r) - wmax * (uR(1:nw_flux) - uL(1:nw_flux))#{if defined('FLUX_AD')}# * phi#{endif}#)
  end subroutine riemann_llf_prim


  !> One-face HLL numerical flux from primitive L/R states.
  !> takes flux_type(flux_dim,nw_flux) into account (fallback to LLF for selected variables)
  subroutine riemann_hll_prim(uL, uR, xC, flux_dim, F)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: uL(nw_phys), uR(nw_phys)
    real(dp), intent(in)    :: xC(ndim)
    integer,  intent(in)    :: flux_dim
    real(dp), intent(out)   :: F(nw_flux)

    real(dp) :: flux_l(nw_flux), flux_r(nw_flux)
    real(dp) :: sL, sR, ds, wmax
    real(dp), parameter :: eps = 1e-14_dp

    call get_flux(uL, xC, flux_dim, flux_l)
    call get_flux(uR, xC, flux_dim, flux_r)

    call estimate_speeds_minmax(uL, uR, xC, flux_dim, sL, sR)
    wmax = max(abs(sL), abs(sR))  ! for fallback

    call to_conservative(uL)
    call to_conservative(uR)

    if (sL .ge. 0._dp) then
       F = flux_l
    else if (sR .le. 0._dp) then
       F = flux_r
    else
       ds = sR - sL
       if (ds > eps * (abs(sL) + abs(sR))) then
          F = (sR*flux_l - sL*flux_r + sL*sR*(uR(1:nw_flux) - uL(1:nw_flux))) / ds
       else  ! Fall back to LLF
          F = 0.5_dp * ((flux_l + flux_r) - wmax * (uR(1:nw_flux) - uL(1:nw_flux)))
       end if
    end if
    
    call apply_flux_type_tvdlf(flux_l, flux_r, uL, uR, wmax, flux_dim, F)

  end subroutine riemann_hll_prim


  !> Overwrite F with the LLF flux for every variable the physics module marks
  !> flux_tvdlf in flux_type(flux_dim,:) (in mhd: B_n and the GLM psi).
  !> uL/uR are conservative.
  subroutine apply_flux_type_tvdlf(flux_l, flux_r, uL, uR, wmax, flux_dim, F)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)    :: flux_l(nw_flux), flux_r(nw_flux)
    real(dp), intent(in)    :: uL(nw_phys), uR(nw_phys)
    real(dp), intent(in)    :: wmax
    integer,  intent(in)    :: flux_dim
    real(dp), intent(inout) :: F(nw_flux)

    integer  :: iw

    ! overwrite with tvdlf flux if requested (faster this way)
    if ( any(flux_type(flux_dim,:) == flux_tvdlf) ) then
       do iw = 1, nw_flux
          if (flux_type(flux_dim,iw) == flux_tvdlf) then
             F(iw) = 0.5_dp * ((flux_l(iw) + flux_r(iw)) - wmax * (uR(iw) - uL(iw)))
          end if
       end do
    end if

  end subroutine apply_flux_type_tvdlf


  !> One-face HLLC numerical flux from primitive L/R states.
  !> The physics-independent frame: the outer signal speeds come from the
  !> physics module's estimate_speeds_hllc, the star-region flux from its
  !> hllc_flux (hd: Toro (2010) chapter 10, variant 2; mhd: Li (2005)).
  !> Falls back to LLF where hllc_flux reports a degenerate fan, and takes
  !> flux_type(flux_dim,nw_flux) into account like riemann_hll_prim.
  subroutine riemann_hllc_prim(uL, uR, xC, flux_dim, F)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: uL(nw_phys), uR(nw_phys)
    real(dp), intent(in)    :: xC(ndim)
    integer,  intent(in)    :: flux_dim
    real(dp), intent(out)   :: F(nw_flux)

    real(dp) :: flux_l(nw_flux), flux_r(nw_flux)
    real(dp) :: sL, sR, wmax
    logical  :: degenerate

    call get_flux(uL, xC, flux_dim, flux_l)
    call get_flux(uR, xC, flux_dim, flux_r)

    call estimate_speeds_hllc(uL, uR, xC, flux_dim, sL, sR)
    wmax = max(abs(sL), abs(sR))

    ! converts uL/uR to conservative
    call hllc_flux(uL, uR, flux_l, flux_r, sL, sR, flux_dim, F, degenerate)

    if (degenerate) then
      F = 0.5_dp * ((flux_l + flux_r) - wmax * (uR(1:nw_flux) - uL(1:nw_flux)))
    end if

    call apply_flux_type_tvdlf(flux_l, flux_r, uL, uR, wmax, flux_dim, F)

  end subroutine riemann_hllc_prim


  !> MUSCL (primitive-variable) reconstruction with slope limiter; HLL two-wave approximate Riemann flux at faces.
  !> Uses estimated left/right signal speeds (Davis (1988)) for less diffusion than LLF, no contact resolution.
  subroutine reconflux_muscl_hll_prim(u, xlocC, flux_dim, flux, typelim)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys, 5)
    real(dp), intent(in)  :: xlocC(1:ndim, 2)
    integer, intent(in)   :: flux_dim, typelim
    real(dp), intent(out) :: flux(nw_flux, 2)

    real(dp) :: uL(nw_phys,2), uR(nw_phys,2)
    integer  :: iface

    call muscl_reconstruct_prim(u, typelim, uL, uR)

    do iface=1,2
      call riemann_hll_prim(uL(:,iface), uR(:,iface), xlocC(:,iface), flux_dim, flux(:,iface))
    end do
  end subroutine reconflux_muscl_hll_prim


  !> MUSCL (primitive-variable) reconstruction with slope limiter; LLF/Rusanov numerical flux at faces.
  !> Robust and diffusive; uses local max wavespeed for upwinding.
  !> Adaptive diffusion (Rempel et al. 2009) is compiled in only when FLUX_AD is defined.
  subroutine reconflux_muscl_llf_prim(u, xlocC, flux_dim, flux, typelim)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys, 5)
    real(dp), intent(in)  :: xlocC(1:ndim, 2)
    integer, intent(in)   :: flux_dim, typelim
    real(dp), intent(out) :: flux(nw_flux, 2)

    real(dp) :: uL(nw_phys,2), uR(nw_phys,2)
    integer  :: iface
#:if defined('FLUX_AD')
    real(dp) :: phi(nw_flux,2)
    real(dp) :: delta_rc, delta_ct, phi_v
    integer  :: iw
#:endif
#:if defined('MHD_ENERGY_ONLY')
    real(dp) :: Fe
#:endif

    call muscl_reconstruct_prim(u, typelim, uL, uR)

#:if defined('FLUX_AD')
    ! Rempel et al. (2009): reduce dissipation where reconstruction and cell-centre agree in sign
    phi = 1.0_dp
    do iface = 1, 2
      do iw = 1, nw_flux
        delta_rc = uR(iw,iface) - uL(iw,iface)
        delta_ct = u(iw, iface+2) - u(iw, iface+1)
        if (delta_rc * delta_ct > 1.0e-18_dp) then
          phi_v = flux_ad_scale * delta_rc**2 / (delta_ct**2 + 1.0e-18_dp)
          phi(iw,iface) = max(flux_ad_min, min(phi_v, 1.0_dp))
        end if
      end do
    end do
#:endif

    do iface = 1, 2
      call riemann_llf_prim(uL(:,iface), uR(:,iface), xlocC(:,iface), flux_dim, flux(:,iface)#{if defined('FLUX_AD')}#, phi(:,iface)#{endif}#)
    end do

#:if defined('MHD_ENERGY_ONLY')
    do iface = 1, 2
      Fe = flux(iw_e, iface)
      flux(:, iface) = 0.0_dp
      flux(iw_e, iface) = Fe
    end do
#:endif

  end subroutine reconflux_muscl_llf_prim


  !> MUSCL (primitive-variable) reconstruction with slope limiter; HLLC approximate Riemann flux at faces.
  !> Restores the contact wave (and shear in Euler/HD), typically sharper than HLL for similar cost.
  subroutine reconflux_muscl_hllc_prim(u, xlocC, flux_dim, flux, typelim)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys, 5)
    real(dp), intent(in)  :: xlocC(1:ndim, 2)
    integer, intent(in)   :: flux_dim, typelim
    real(dp), intent(out) :: flux(nw_flux, 2)

    real(dp) :: uL(nw_phys,2), uR(nw_phys,2)
    integer  :: iface

    call muscl_reconstruct_prim(u, typelim, uL, uR)
    do iface=1,2
      call riemann_hllc_prim(uL(:,iface), uR(:,iface), xlocC(:,iface), flux_dim, flux(:,iface))
    end do
  end subroutine reconflux_muscl_hllc_prim


  pure real(dp) function vanleer(a, b) result(phi)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in) :: a, b
    real(dp)             :: ab

    ab = a * b
    if (ab > 0) then
       phi = 2.0_dp * ab / (a + b)
    else
       phi = 0.0_dp
    end if
  end function vanleer

  pure real(dp) function minmod(a, b)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in) :: a, b

    if (a * b <= 0) then
       minmod = 0.0_dp
    else if (abs(a) < abs(b)) then
       minmod = a
    else
       minmod = b
    end if
  end function minmod

  !  scavenged from afivo, thanks Jannis :-) !
  !> Modified implementation of Koren limiter, to avoid division and the min/max
  !> functions, which can be problematic / expensive. In most literature, you
  !> have r = a / b (ratio of gradients). Then the limiter phi(r) is multiplied
  !> with b. With this implementation, you get phi(r) * b
  pure real(dp) function koren(a, b)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in) :: a  !< Density gradient (numerator)
    real(dp), intent(in) :: b  !< Density gradient (denominator)
    real(dp), parameter  :: third = 1/3.0_dp
    real(dp)             :: aa, ab

    aa = a * a
    ab = a * b

    if (ab <= 0) then
       ! a and b have different sign or one of them is zero, so r is either 0,
       ! inf or negative (special case a == b == 0 is ignored)
       koren = 0
    else if (aa <= 0.25_dp * ab) then
       ! 0 < a/b <= 1/4, limiter has value 2*a/b
       koren = 2*a
    else if (aa <= 2.5_dp * ab) then
       ! 1/4 < a/b <= 2.5, limiter has value (1+2*a/b)/3
       koren = third * (b + 2*a)
    else
       ! (1+2*a/b)/6 >= 1, limiter has value 2
       koren = 2*b
    end if
    
  end function koren
  
  !> Monotonised central-difference limiter with tunable beta (AMRVAC's
  !> mcbeta; beta=2 recovers the classic MC limiter, van Leer 1979)
  pure real(dp) function mcbeta(a, b) result(phi)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in) :: a, b
  real(dp), parameter  :: c_mcbeta = 1.4_dp
  real(dp)             :: ab

  ab = a * b
  if (ab > 0) then
    phi = sign(1.0_dp, a + b) * min(c_mcbeta*abs(a), c_mcbeta*abs(b), abs(a+b)*0.5_dp)
  else
    phi = 0.0_dp
  end if

  end function mcbeta

end module mod_finite_volume
