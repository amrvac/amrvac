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
  ! finite-volume layouts, selected at run time by fv_face_flux
#:set layouts = ['shared', 'per_cell']

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
      select case (fv_layout)
#:for layout in layouts
      case (fv_${layout}$)
        call finite_volume_local_${scheme_tag}$_${layout}$(qdt, dtfactor, ixImin1,ixImin2,&
              ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
              ixOmax3, idimsmin,idimsmax, qtC, bga, qt, bgb)
#:endfor
      end select
#:endfor
    case default
      call mpistop("finite_volume_local: unknown flux scheme")
    end select
end subroutine finite_volume_local

! Index strings for the direction-d loops below: the lower neighbour of cell
! (ix1,ix2,ix3), and what fix_conserve indexes pflux and fC* by.
#:set lower = {1: 'ix1-1, ix2, ix3', 2: 'ix1, ix2-1, ix3', 3: 'ix1, ix2, ix3-1'}
#:set pfx = {1: '1,ix2-nghostcells,ix3-nghostcells', 2: 'ix1-nghostcells,1,ix3-nghostcells', 3: 'ix1-nghostcells,ix2-nghostcells,1'}
#:set fcx = {1: 'ix2,ix3', 2: 'ix1,ix3', 3: 'ix1,ix2'}

#! The numerical flux f on one face in direction d of cell (ix1,ix2,ix3): its
#! upper face ix_d+1/2 for s = 0, its lower face ix_d-1/2 for s = -1. Both
#! layouts build every face through this, so they hand the flux routine the
#! same stencil and position for the same face.
#:def FACE_FLUX(d, s, f, faceflux_proc)
#:set lo = '-1' if s == -1 else ''
#:set cell = ', '.join('ix%d%s' % (k, lo if k == d else '') for k in range(1, 4))
#:set sten = ', '.join('ix%d%+d:ix%d%+d' % (k, s-1, k, s+2) if k == d else 'ix%d' % k for k in range(1, 4))
                   ! face ix${d}$${lo}$+1/2, between cell ix${d}$${lo}$ and cell ix${d}$${lo}$+1
                   tmp4 = uprim(1:nw_phys, ${sten}$)
                   ! Transverse components on the lower cell's barycentre, the
                   ! normal one on the face, in the lower-corner form of cell
                   ! ix${d}$${lo}$+1: the one form right in every geometry (cannot be
                   ! x +/- dr/2: bgeo%x is a volume barycentre in a curvilinear
                   ! build, and dr is a ratio in ln(r) under LOG_RADIUS).
                   xC(1:ndim) = bgeo%x(${cell}$, 1:ndim, n)
                   xC(${d}$) = xlo(${d}$) + dble(ix${d}$${lo}$-nghostcells)*dr(${d}$)
#:if d == 1 and defined('LOG_RADIUS')
                   ! the logical radial coordinate is ln(1 + r/r0); this is
                   ! r_of_s in the reduced form that holds wherever s >= 0,
                   ! which covers these faces: domain faces, s = 0 being r = 0
                   xC(1) = log_ra*dexp(xC(1)) + log_rb
#:endif
                   call ${faceflux_proc}$(tmp4, xC, ${d}$, ${f}$, typelim)
#:enddef

#! Update cell (ix1,ix2,ix3) with the difference of its lower and upper face
#! fluxes flo and fhi in direction d, and store the block-boundary fluxes for
#! flux fixing. What is stored differs by geometry: qdt*f/dx, the face's own
#! contribution to the cell update, on a Cartesian mesh, and the extensive
#! qdt*f*A on a curvilinear one. The extensive form is what lets fix_conserve
#! divide by the *coarse* cell's volume and sum the four fine faces without a
#! volume ratio, since the fine face areas add up to the coarse one exactly.
#:def FV_CELL_UPDATE(d, flo, fhi)
#:if GEOM == 'Cartesian'
                   bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = bgb%w(ix1, ix2, ix3, 1:nw_flux,&
                        n) + qdt * (${flo}$ &
                        - ${fhi}$) * inv_dr(${d}$)
#:else
                   inv_dvol = 1.0_dp / bgeo%dvolume(ix1, ix2, ix3, n)
                   bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = bgb%w(ix1, ix2, ix3, 1:nw_flux,&
                        n) + qdt * (${flo}$ * bgeo%surfaceC(${lower[d]}$, ${d}$, n) &
                        - ${fhi}$ * bgeo%surfaceC(ix1, ix2, ix3, ${d}$, n)) * inv_dvol
#:endif

                   select case (neighbor_type_${d}$m)
                   case (neighbor_fine)
#:if GEOM == 'Cartesian'
                      if (ix${d}$.eq.ixOmin${d}$) pflux(1,${d}$)%flux(${pfx[d]}$,1:nw_flux,n) &
                           = qdt * inv_dr(${d}$) * ${flo}$
#:else
                      if (ix${d}$.eq.ixOmin${d}$) pflux(1,${d}$)%flux(${pfx[d]}$,1:nw_flux,n) &
                           = qdt * bgeo%surfaceC(${lower[d]}$, ${d}$, n) * ${flo}$
#:endif
                   case (neighbor_coarse)
#:if GEOM == 'Cartesian'
                      if (ix${d}$.eq.ixOmin${d}$) fC${d}$(1,${fcx[d]}$,1:nw_flux) = - qdt * inv_dr(${d}$) * ${flo}$
#:else
                      if (ix${d}$.eq.ixOmin${d}$) fC${d}$(1,${fcx[d]}$,1:nw_flux) = - qdt * bgeo%surfaceC(${lower[d]}$, ${d}$, n) * ${flo}$
#:endif
                   end select

                   select case (neighbor_type_${d}$p)
                   case (neighbor_fine)
#:if GEOM == 'Cartesian'
                      if (ix${d}$.eq.ixOmax${d}$) pflux(2,${d}$)%flux(${pfx[d]}$,1:nw_flux,n) &
                           = - qdt * inv_dr(${d}$) * ${fhi}$
#:else
                      if (ix${d}$.eq.ixOmax${d}$) pflux(2,${d}$)%flux(${pfx[d]}$,1:nw_flux,n) &
                           = - qdt * bgeo%surfaceC(ix1, ix2, ix3, ${d}$, n) * ${fhi}$
#:endif
                   case (neighbor_coarse)
#:if GEOM == 'Cartesian'
                      if (ix${d}$.eq.ixOmax${d}$) fC${d}$(2,${fcx[d]}$,1:nw_flux) = qdt * inv_dr(${d}$) * ${fhi}$
#:else
                      if (ix${d}$.eq.ixOmax${d}$) fC${d}$(2,${fcx[d]}$,1:nw_flux) = qdt * bgeo%surfaceC(ix1, ix2, ix3, ${d}$, n) * ${fhi}$
#:endif
                   end select
#:enddef

#! One finite-volume kernel per flux scheme and layout. layout 'shared'
#! computes each face once into the per-gang fface and lets both cells read
#! it; 'per_cell' has each cell compute both of its faces, which costs a
#! second flux per face but needs no fface (fv_face_flux in &methodlist).
#:def FV_KERNEL(scheme_tag, faceflux_proc, layout)
  subroutine finite_volume_local_${scheme_tag}$_${layout}$(qdt, dtfactor, ixImin1,ixImin2,&
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
    ! One face: its 4-cell stencil, its position and its numerical flux
    real(dp)               :: tmp4(nw_phys,4)
    real(dp)               :: xC(ndim)
    real(dp)               :: fx(nw_flux)
#:if layout == 'shared'
    ! The numerical fluxes of one direction's faces, indexed by the cell on
    ! the face's lower side: face ix+1/2 is fface(ix). Reused for all three
    ! directions; per gang, like uprim.
    real(dp)               :: fface(ixOmin1-1:ixOmax1, ixOmin2-1:ixOmax2,&
       ixOmin3-1:ixOmax3, 1:nw_flux)
#:else
    ! The flux on the cell's lower face; fx holds the upper one
    real(dp)               :: fxlo(nw_flux)
#:endif
    real(dp)               :: inv_dr(ndim)
    real(dp)               :: dr(ndim)
    integer                :: typelim
    real(dp)               :: xloc(ndim)
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

       ${GPU_PARALLEL_LOOP_GANG("private(uprim" + (", fface" if layout == 'shared' else "") + ", inv_dr, dr, n, xlo, typelim, ix1, ix2, ix3, fC1, fC2, fC3, neighbor_type_1m, neighbor_type_1p, neighbor_type_2m, neighbor_type_2p, neighbor_type_3m, neighbor_type_3p)")}$
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

#:for d in range(1, 4)
          ! ---------------------------------------------------------------
#:if layout == 'shared'
          ! Direction ${d}$. Every face gets one reconstruction and one Riemann
          ! flux, shared by the two cells on either side of it: the faces
          ! ixOmin${d}$-1/2 .. ixOmax${d}$+1/2, i.e. fface(ixOmin${d}$-1:ixOmax${d}$).
          ${GPU_LOOP_VECTOR("collapse(ndim) private(tmp4, xC, fx)")}$
          do ix3=ixOmin3${'-1' if d == 3 else ''}$,ixOmax3
             do ix2=ixOmin2${'-1' if d == 2 else ''}$,ixOmax2
                do ix1=ixOmin1${'-1' if d == 1 else ''}$,ixOmax1
$:FACE_FLUX(d, 0, 'fx', faceflux_proc)
                   fface(ix1, ix2, ix3, 1:nw_flux) = fx
                end do
             end do
          end do

          ${GPU_LOOP_VECTOR("collapse(ndim)" + (" private(inv_dvol)" if GEOM != 'Cartesian' else ""))}$
          do ix3=ixOmin3,ixOmax3
             do ix2=ixOmin2,ixOmax2
                do ix1=ixOmin1,ixOmax1
$:FV_CELL_UPDATE(d, 'fface(%s, 1:nw_flux)' % lower[d], 'fface(ix1, ix2, ix3, 1:nw_flux)')
                end do
             end do
          end do
#:else
          ! Direction ${d}$. Each cell computes both of its faces itself, so
          ! every interior face is computed twice, but no per-block face
          ! buffer is needed.
          ${GPU_LOOP_VECTOR("collapse(ndim) private(tmp4, xC, fxlo, fx" + (", inv_dvol" if GEOM != 'Cartesian' else "") + ")")}$
          do ix3=ixOmin3,ixOmax3
             do ix2=ixOmin2,ixOmax2
                do ix1=ixOmin1,ixOmax1
$:FACE_FLUX(d, -1, 'fxlo', faceflux_proc)
$:FACE_FLUX(d, 0, 'fx', faceflux_proc)
$:FV_CELL_UPDATE(d, 'fxlo(1:nw_flux)', 'fx(1:nw_flux)')
                end do
             end do
          end do
#:endif

#:endfor
#:if GEOM != 'Cartesian' or defined('SOURCE_LOCAL') or defined('SOURCE_NONLOCAL') or defined('SOURCE_COMPACT')
          ! Geometric and optional source terms, on the state all three
          ! directions have updated, as before
          ${GPU_LOOP_VECTOR("collapse(ndim) private(wnew, tmp, xloc" + (", wCT" if defined('SOURCE_LOCAL') else "") + (", wprim" if defined('SOURCE_LOCAL') or GEOM != 'Cartesian' else "") + (", inv_dvol, dAdV" if GEOM != 'Cartesian' else "") + (", tmp1,tmp2,tmp3" if defined('SOURCE_COMPACT') else "") + (", dloc" if defined('SOURCE_LOCAL') or defined('SOURCE_NONLOCAL') or defined('SOURCE_COMPACT') else "") + ")")}$
          do ix3=ixOmin3,ixOmax3
             do ix2=ixOmin2,ixOmax2
                do ix1=ixOmin1,ixOmax1
                   xloc(1:ndim) = bgeo%x(ix1, ix2, ix3, 1:ndim, n)
#:if GEOM != 'Cartesian'
                   inv_dvol = 1.0_dp / bgeo%dvolume(ix1, ix2, ix3, n)
#:endif

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
#:endif

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

  end subroutine finite_volume_local_${scheme_tag}$_${layout}$
  #:enddef


#:for scheme_tag, method_enum, faceflux_proc in schemes
#:for layout in layouts
#:call FV_KERNEL(scheme_tag, faceflux_proc, layout)
#:endcall
#:endfor
#:endfor


  !> MUSCL reconstruction in primitive variables at one face from a 4-point
  !> stencil: the face lies between cells 2 and 3, uL is cell 2 plus half its
  !> limited slope towards it, uR cell 3 minus half its own.
  pure subroutine muscl_reconstruct_prim(u, typelim, uL, uR)
    use mod_limiter, only: limiter_minmod, limiter_vanleer, limiter_mcbeta, limiter_koren
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys,4)
    integer,  intent(in)  :: typelim
    real(dp), intent(out) :: uL(nw_phys), uR(nw_phys)

    integer  :: iw

    select case (typelim)
    case (limiter_minmod)
      do iw=1,nw_recon
        uL(iw) = u(iw,2) + 0.5_dp*minmod(u(iw,2)-u(iw,1), u(iw,3)-u(iw,2))
        uR(iw) = u(iw,3) - 0.5_dp*minmod(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
      end do
    case (limiter_vanleer)
      do iw=1,nw_recon
        uL(iw) = u(iw,2) + 0.5_dp*vanleer(u(iw,2)-u(iw,1), u(iw,3)-u(iw,2))
        uR(iw) = u(iw,3) - 0.5_dp*vanleer(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
      end do
    case (limiter_mcbeta)
      do iw=1,nw_recon
        uL(iw) = u(iw,2) + 0.5_dp*mcbeta(u(iw,2)-u(iw,1), u(iw,3)-u(iw,2))
        uR(iw) = u(iw,3) - 0.5_dp*mcbeta(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
      end do
    case (limiter_koren)
      ! koren(a,b) = phi(a/b)*b is one-sided: towards a cell's right face it
      ! takes (dw+, dw-), towards its left face (dw-, dw+). Swapping them
      ! breaks mirror symmetry of the scheme.
      do iw=1,nw_recon
        uL(iw) = u(iw,2) + 0.5_dp*koren(u(iw,3)-u(iw,2), u(iw,2)-u(iw,1))
        uR(iw) = u(iw,3) - 0.5_dp*koren(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
      end do
    case default  ! Fallback: Godunov (piecewise constant)
      uL(:) = u(:,2)
      uR(:) = u(:,3)
    end select

    ! The slopes were limited slot by slot, so a state can come out
    ! inadmissible or -- where the physics carries auxiliary variables --
    ! internally inconsistent. Repair it here, once, rather than in each of the
    ! three Riemann solvers: this is the single producer of face states.
    call fix_prim_state(uL)
    call fix_prim_state(uR)
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
  subroutine reconflux_muscl_hll_prim(u, xC, flux_dim, flux, typelim)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys, 4)
    real(dp), intent(in)  :: xC(1:ndim)
    integer, intent(in)   :: flux_dim, typelim
    real(dp), intent(out) :: flux(nw_flux)

    real(dp) :: uL(nw_phys), uR(nw_phys)

    call muscl_reconstruct_prim(u, typelim, uL, uR)
    call riemann_hll_prim(uL, uR, xC, flux_dim, flux)
  end subroutine reconflux_muscl_hll_prim


  !> MUSCL (primitive-variable) reconstruction with slope limiter; LLF/Rusanov numerical flux at faces.
  !> Robust and diffusive; uses local max wavespeed for upwinding.
  !> Adaptive diffusion (Rempel et al. 2009) is compiled in only when FLUX_AD is defined.
  subroutine reconflux_muscl_llf_prim(u, xC, flux_dim, flux, typelim)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys, 4)
    real(dp), intent(in)  :: xC(1:ndim)
    integer, intent(in)   :: flux_dim, typelim
    real(dp), intent(out) :: flux(nw_flux)

    real(dp) :: uL(nw_phys), uR(nw_phys)
#:if defined('FLUX_AD')
    real(dp) :: phi(nw_flux)
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
    do iw = 1, nw_flux
      delta_rc = uR(iw) - uL(iw)
      delta_ct = u(iw, 3) - u(iw, 2)
      if (delta_rc * delta_ct > 1.0e-18_dp) then
        phi_v = flux_ad_scale * delta_rc**2 / (delta_ct**2 + 1.0e-18_dp)
        phi(iw) = max(flux_ad_min, min(phi_v, 1.0_dp))
      end if
    end do
#:endif

    call riemann_llf_prim(uL, uR, xC, flux_dim, flux#{if defined('FLUX_AD')}#, phi#{endif}#)

#:if defined('MHD_ENERGY_ONLY')
    Fe = flux(iw_e)
    flux(:) = 0.0_dp
    flux(iw_e) = Fe
#:endif

  end subroutine reconflux_muscl_llf_prim


  !> MUSCL (primitive-variable) reconstruction with slope limiter; HLLC approximate Riemann flux at faces.
  !> Restores the contact wave (and shear in Euler/HD), typically sharper than HLL for similar cost.
  subroutine reconflux_muscl_hllc_prim(u, xC, flux_dim, flux, typelim)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys, 4)
    real(dp), intent(in)  :: xC(1:ndim)
    integer, intent(in)   :: flux_dim, typelim
    real(dp), intent(out) :: flux(nw_flux)

    real(dp) :: uL(nw_phys), uR(nw_phys)

    call muscl_reconstruct_prim(u, typelim, uL, uR)
    call riemann_hllc_prim(uL, uR, xC, flux_dim, flux)
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
