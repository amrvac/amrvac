!=============================================================================
!> Transport kernel of the short-characteristics (SC) radiative transfer
!>
!> The transfer equation is solved along one ray direction at a time with the linear short
!> characteristics of Olson and Kunasz. Each cell takes the intensity at the upwind point on the
!> face the ray crosses first, interpolated linearly along that face, and integrates across one
!> step with the source taken as linear in optical depth. sc_sweep_block sweeps one block whose
!> ghost ring already holds the upwind values and the inflows, and it is the only sweep in the
!> code. The wavefront solvers call it block by block, and the gather calls it once on the whole
!> replicated band, which sc_pad_grid has padded with the inflows
!>
!> The module also builds the angular quadrature, Gauss-Legendre or the Carlson A4 set, and in 2D
!> merges the rays that differ only out of the plane. It holds the replicated single-level grid
!> of the gather and forms the heating in blend_heating, which blends c chi (J - B) with -div F.
!> Nothing here communicates, since the drivers fill chi and B beforehand
!>
!> The gather is the original implementation, kept for comparison and not for production. Every
!> rank holds the whole band, and its periodic seam is not exact (see sc_sweep_2d). Production runs
!> should use KBA, which is far better optimised
!>
!> Structure:
!>   A. quadrature    rt_sc_init, rt_sc_dedup_rays, build_gl, build_carlson and
!>                    gauss_legendre_nodes;
!>   B. heating       blend_heating, and rt_sc_heat_blend on the replicated grid;
!>   C. transport     sc_pad_grid and sc_sweep_2d for the replicated grid, ka_coeffs for the
!>                    weights of one step, and sc_sweep_block
!>
!> Jack Jenkins (10.2026); More or less finished in v1 state for amrvac4.0
!=============================================================================
module mod_rt_sc_kernel
  use mod_global_parameters, only: bigdouble, smalldouble
  implicit none

  !> The angular quadrature, either 'GL' (Gauss-Legendre polar nodes times equally spaced azimuths)
  !> or 'Carlson', the level-symmetric A4 set that MURaM uses
  character(len=16), public :: rt_sc_quad = 'GL'
  !> Azimuthal directions and polar nodes of the GL set, with nmu one of 2, 4, 6 or 8
  integer, public :: rt_sc_ndir = 6
  integer, public :: rt_sc_nmu = 4
  integer, public :: rt_sc_nray = 0
  !> Stratification axis (1 = x, 2 = y), set by the driver before rt_sc_init. The GL set places
  !> its polar nodes about this axis, while Carlson is symmetric and does not need it
  integer, public :: rt_sc_polar_axis = 2
  !> Weight of the heating blend, w = exp(-blend_k tau), and the scale of the taper
  !> tau^2/(tau^2 + tau_min^2) that switches the SC heating off where the column above is thin
  double precision, public :: rt_sc_blend_k = 10.0d0, rt_sc_tau_min = 1.0d-6
  !> Incoming intensity at the lower and upper ends of the current band, thermal (I = B) when set
  !> and open (I = 0) otherwise
  logical, public :: bc_lo_thermal = .true., bc_hi_thermal = .false.
  !> The horizontal edges wrap when the domain is periodic and otherwise act as inflow boundaries,
  !> exactly like the band ends. rt_sc_setup copies these from periodB(hdim)
  logical, public :: bc_h_periodic  = .true.
  logical, public :: bc_hlo_thermal = .false., bc_hhi_thermal = .false.
  !> Direction cosines and weights, normalised so that sum(wmu) = 1. In 2D the out-of-plane muz
  !> only lengthens the path of each step through ds
  double precision, allocatable, public :: mux(:), muy(:), muz(:), wmu(:)
  !> The replicated single-level grid with chi = kappa rho, B = a T^4, the mean intensity J, the
  !> intensity of the ray being swept, the heating, both flux components and tau from the top
  double precision, allocatable, public :: chi_g(:,:), b_g(:,:), jm_g(:,:), ii_g(:,:), q_g(:,:)
  double precision, allocatable, public :: fx_g(:,:), fy_g(:,:), tau_g(:,:)
  ! Size of that grid along hdim and updim, left at -1 until the grid is first built
  integer, public :: nxg = -1, nyg = -1
  ! Cell size of the grid along hdim and updim
  double precision, public :: dxg, dyg
  ! chi, B and the intensity of the grid padded with one ghost layer, which holds the band-end and
  ! edge inflows, so that the grid can be swept by sc_sweep_block as a single block
  double precision, allocatable :: chi_p(:,:), b_p(:,:), i_p(:,:)

contains

  !> Sign of a direction cosine as -1, 0 or 1. The zero case occurs, for example for the GL rays
  !> at phi = pi/2
  pure integer function ray_sign(a)
    double precision, intent(in) :: a
    ray_sign = merge(1, 0, a > 0.0d0) - merge(1, 0, a < 0.0d0)
  end function ray_sign

  !> Build the quadrature into mux, muy, muz and wmu as full 3D cosines, according to rt_sc_quad
  subroutine rt_sc_init()
    use mod_global_parameters
    ! Even in 2D we build the full 3D set, because the out-of-plane cosine still sets how long
    ! each in-plane step is
    select case(trim(adjustl(rt_sc_quad)))
    case('Carlson','carlson')
      call build_carlson()
    case default
      call build_gl()
    end select
    ! Rays that differ only in muz are indistinguishable to a 2D sweep, so we merge them and
    ! save the repeated work
    call rt_sc_dedup_rays()
    if(mype==0) print *, "rt_sc_init: quadrature '"//trim(adjustl(rt_sc_quad))//"', rays =", rt_sc_nray
  end subroutine rt_sc_init

  !> Merge rays that share (mux, muy) and add up their weights. Since the 2D sweep and the moments
  !> see only those two cosines, Carlson's +-muz pairs and GL's phi, 2pi-phi pairs sweep the same.
  !> The total weight is preserved, so J, Fx and Fy are unchanged to round-off
  subroutine rt_sc_dedup_rays()
    use mod_global_parameters
    double precision, allocatable :: nx(:), ny(:), nz(:), wsum(:)
    logical, allocatable :: gone(:)
    integer :: k, j, nnew
    ! Cosines that agree to this tolerance are treated as the same direction
    double precision, parameter :: tol = 1.0d-12
    if(rt_sc_nray <= 1) return
    allocate(nx(rt_sc_nray), ny(rt_sc_nray), nz(rt_sc_nray), wsum(rt_sc_nray), gone(rt_sc_nray))
    gone = .false.; nnew = 0
    ! Each ray not yet absorbed opens a new entry in the merged set, and any later ray with the
    ! same in-plane cosines hands its weight over and is marked gone
    do k = 1, rt_sc_nray
      if(gone(k)) cycle
      nnew = nnew + 1
      nx(nnew) = mux(k); ny(nnew) = muy(k); nz(nnew) = muz(k); wsum(nnew) = wmu(k)
      ! w_merged = w_k + the sum of w_j over the later rays j with the same (mux, muy)
      do j = k+1, rt_sc_nray
        if(gone(j)) cycle
        if(abs(mux(j)-mux(k)) < tol .and. abs(muy(j)-muy(k)) < tol) then
          wsum(nnew) = wsum(nnew) + wmu(j); gone(j) = .true.
        end if
      end do
    end do
    ! The quadrature is only replaced when something was actually merged
    if(nnew < rt_sc_nray) then
      deallocate(mux, muy, muz, wmu)
      allocate(mux(nnew), muy(nnew), muz(nnew), wmu(nnew))
      mux = nx(1:nnew); muy = ny(1:nnew); muz = nz(1:nnew); wmu = wsum(1:nnew)
      if(mype==0) print *, "rt_sc_dedup_rays: (mux,muy) rays", rt_sc_nray, "->", nnew, " (2D-identical merged)"
      rt_sc_nray = nnew
    end if
    deallocate(nx, ny, nz, wsum, gone)
  end subroutine rt_sc_dedup_rays

  !> Product Gauss-Legendre set, with nmu polar nodes about rt_sc_polar_axis and ndir azimuths at
  !> the midpoints of equal intervals
  subroutine build_gl()
    use mod_global_parameters
    integer :: im, ip, k, nmu, nphi
    double precision :: phi, s, mu, xmu(8), wgl(8)
    ! With fewer than three azimuths the plane is not resolved, so we never go below three
    if(rt_sc_ndir < 3) rt_sc_ndir = 3
    nphi = rt_sc_ndir;  nmu = rt_sc_nmu
    ! The node routine may round nmu up to the next order it supports
    call gauss_legendre_nodes(nmu, xmu, wgl)
    rt_sc_nray = nmu*nphi
    if(allocated(mux)) deallocate(mux, muy, muz, wmu)
    allocate(mux(rt_sc_nray), muy(rt_sc_nray), muz(rt_sc_nray), wmu(rt_sc_nray))
    k = 0
    do im = 1, nmu
      ! The polar cosine and its sine, mu = x_GL(im),  s = sqrt(1 - mu^2)
      mu = xmu(im)
      s  = dsqrt(max(1.0d0 - mu*mu, 0.0d0))
      do ip = 1, nphi
        ! phi = 2 pi (ip - 1/2)/nphi,  mu_h = s cos(phi),  mu_z = s sin(phi)
        phi = 2.0d0*dpi*(dble(ip)-0.5d0)/dble(nphi)
        k = k + 1
        ! The polar cosine lies along the stratification and the azimuth turns about it
        if(rt_sc_polar_axis == 2) then
          muy(k) = mu; mux(k) = s*dcos(phi)
        else
          mux(k) = mu; muy(k) = s*dcos(phi)
        end if
        muz(k) = s*dsin(phi)
        ! The GL weights add up to 2 on [-1, 1] and the azimuths share each one equally
        ! w = w_GL/(2 nphi), so that sum w = 1
        wmu(k) = 0.5d0*wgl(im)/dble(nphi)
      end do
    end do
  end subroutine build_gl

  !> Carlson A4 level-symmetric set, built from the three cyclic permutations of (mu1, mu2, mu2)
  !> with mu1^2 + 2 mu2^2 = 1 in each of the 8 octants, which gives 24 directions of weight 1/24
  !> as in MURaM's carlson3
  subroutine build_carlson()
    ! m1^2 + 2 m2^2 = 1, so every direction built from these is a unit vector
    double precision, parameter :: m1 = 0.88191710d0, m2 = 0.33333333d0
    double precision :: base(3,3)
    integer :: ib, sgx, sgy, sgz, k
    ! Each column places the long cosine on a different axis
    base(:,1) = (/ m1, m2, m2 /)
    base(:,2) = (/ m2, m1, m2 /)
    base(:,3) = (/ m2, m2, m1 /)
    rt_sc_nray = 24
    if(allocated(mux)) deallocate(mux, muy, muz, wmu)
    allocate(mux(24), muy(24), muz(24), wmu(24))
    k = 0
    ! The sign loops run over the eight octants, and every octant receives all three
    ! permutations with the same weight, as a level-symmetric set requires
    do sgz = -1, 1, 2; do sgy = -1, 1, 2; do sgx = -1, 1, 2
      do ib = 1, 3
        k = k + 1
        ! mu = (sgx b_1, sgy b_2, sgz b_3) for the permutation b,  w = 1/24
        mux(k) = sgx*base(1,ib); muy(k) = sgy*base(2,ib); muz(k) = sgz*base(3,ib)
        wmu(k) = 1.0d0/24.0d0
      end do
    end do; end do; end do
  end subroutine build_carlson

  !> Gauss-Legendre nodes and weights on [-1,1] for n in {2,4,6,8}. Any other n is moved up to
  !> the next supported order
  subroutine gauss_legendre_nodes(n, x, w)
    integer, intent(inout) :: n
    double precision, intent(out) :: x(8), w(8)
    x = 0.0d0; w = 0.0d0
    ! Only four orders are ever needed, so they are tabulated to double precision rather than
    ! computed by root finding
    select case(n)
    case(1,2)
      n = 2
      x(1:2) = (/-0.5773502691896257d0, 0.5773502691896257d0/)
      w(1:2) = (/ 1.0000000000000000d0, 1.0000000000000000d0/)
    case(3,4)
      n = 4
      x(1:4) = (/-0.8611363115940526d0,-0.3399810435848563d0, 0.3399810435848563d0, 0.8611363115940526d0/)
      w(1:4) = (/ 0.3478548451374538d0, 0.6521451548625461d0, 0.6521451548625461d0, 0.3478548451374538d0/)
    case(5,6)
      n = 6
      x(1:6) = (/-0.9324695142031521d0,-0.6612093864662645d0,-0.2386191860831969d0, &
                  0.2386191860831969d0, 0.6612093864662645d0, 0.9324695142031521d0/)
      w(1:6) = (/ 0.1713244923791704d0, 0.3607615730481386d0, 0.4679139345726910d0, &
                  0.4679139345726910d0, 0.3607615730481386d0, 0.1713244923791704d0/)
    case default
      n = 8
      x(1:8) = (/-0.9602898564975363d0,-0.7966664774136267d0,-0.5255324099163290d0,-0.1834346424956498d0, &
                  0.1834346424956498d0, 0.5255324099163290d0, 0.7966664774136267d0, 0.9602898564975363d0/)
      w(1:8) = (/ 0.1012285362903763d0, 0.2223810344533745d0, 0.3137066458778873d0, 0.3626837833783620d0, &
                  0.3626837833783620d0, 0.3137066458778873d0, 0.2223810344533745d0, 0.1012285362903763d0/)
    end select
  end subroutine gauss_legendre_nodes

  {^IFTWOD
  !> Heating of one cell. The thin imbalance c chi (J - B) is blended with the thick, conservative
  !> -div F through w = exp(-rt_sc_blend_k tau), and the taper tau^2/(tau^2 + tau_min^2) removes
  !> the heating where the column above is optically thin
  pure double precision function blend_heating(chi, jm, bpl, divf, tau) result(q)
    use mod_global_parameters, only: c_norm
    double precision, intent(in) :: chi, jm, bpl, divf, tau
    double precision :: w
    ! w = exp(-k tau)
    ! Q = [w c chi (J - B) + (1 - w)(-div F)] tau^2/(tau^2 + tau_min^2)
    w = dexp(-rt_sc_blend_k*tau)
    q = (w*(c_norm*chi*(jm-bpl)) + (1.0d0-w)*(-divf)) * tau*tau/(tau*tau + rt_sc_tau_min**2)
  end function blend_heating

  !> Heating q_g, blending the thin imbalance c chi (J - B) with the thick, conservative -div F
  !> through w = exp(-rt_sc_blend_k tau), tau measured from the top. Taking -div F at depth avoids
  !> the stiff, cancelling J - B there. Needs jm_g, b_g, fx_g, fy_g and chi_g
  subroutine rt_sc_heat_blend()
    use mod_global_parameters
    integer :: i, j, ip, im
    double precision :: divf
    ! The blend and the taper both depend on how far below the top each cell sits in optical
    ! depth, so we integrate tau down every column first. The top cell centre lies half a cell
    ! below the band top, and each step down adds the trapezoid of chi between two centres
    do i = 1, nxg
      ! tau_top = chi_top dy/2,  tau_j = tau_(j+1) + (chi_(j+1) + chi_j) dy/2
      tau_g(i,nyg) = 0.5d0*chi_g(i,nyg)*dyg
      do j = nyg-1, 1, -1
        tau_g(i,j) = tau_g(i,j+1) + 0.5d0*(chi_g(i,j+1)+chi_g(i,j))*dyg
      end do
    end do
    do j = 1, nyg
      do i = 1, nxg
        ! The flux divergence uses centred differences, wrapping round in x when the domain is
        ! periodic and falling back to one-sided differences at a closed edge or a band end.
        ! This latter point introduces non negligible errors since we do not illuminate from the
        ! sides...
        ! div F = (Fx_(i+1) - Fx_(i-1))/(2 dx) + (Fy_(j+1) - Fy_(j-1))/(2 dy)
        if(bc_h_periodic) then
          ip = i+1; if(ip>nxg) ip=1
          im = i-1; if(im<1)   im=nxg
          divf = (fx_g(ip,j)-fx_g(im,j))/(2.0d0*dxg)
        else
          ip = min(i+1,nxg); im = max(i-1,1)
          divf = (fx_g(ip,j)-fx_g(im,j))/(dble(ip-im)*dxg)
        end if
        if(j==1) then
          divf = divf + (fy_g(i,2)-fy_g(i,1))/dyg
        else if(j==nyg) then
          divf = divf + (fy_g(i,nyg)-fy_g(i,nyg-1))/dyg
        else
          divf = divf + (fy_g(i,j+1)-fy_g(i,j-1))/(2.0d0*dyg)
        end if
        q_g(i,j) = blend_heating(chi_g(i,j), jm_g(i,j), b_g(i,j), divf, tau_g(i,j))
      end do
    end do
  end subroutine rt_sc_heat_blend

  !> Allocate the replicated grid at nxg by nyg, and keep it when its size is unchanged
  subroutine rt_grid_alloc()
    if(allocated(chi_g)) then
      if(size(chi_g,1)==nxg .and. size(chi_g,2)==nyg) return
      deallocate(chi_g, b_g, jm_g, ii_g, q_g, fx_g, fy_g, tau_g)
    end if
    allocate(chi_g(nxg,nyg), b_g(nxg,nyg), jm_g(nxg,nyg), ii_g(nxg,nyg), q_g(nxg,nyg), &
             fx_g(nxg,nyg), fy_g(nxg,nyg), tau_g(nxg,nyg))
  end subroutine rt_grid_alloc

  !> Pad the replicated grid with one ghost layer for sc_sweep_block, once per band after chi_g
  !> and b_g are filled. Beyond the band ends chi and B take the edge rows and the intensity is
  !> the band-end inflow, B for a thermal end and zero for an open one. A closed horizontal edge
  !> is padded the same way, while a periodic one needs no ghosts of course
  subroutine sc_pad_grid()
    integer :: j
    if(allocated(chi_p)) then
      if(size(chi_p,1) /= nxg+2 .or. size(chi_p,2) /= nyg+2) deallocate(chi_p, b_p, i_p)
    end if
    if(.not.allocated(chi_p)) allocate(chi_p(0:nxg+1,0:nyg+1), b_p(0:nxg+1,0:nyg+1), i_p(0:nxg+1,0:nyg+1))
    chi_p(1:nxg,1:nyg) = chi_g; b_p(1:nxg,1:nyg) = b_g
    chi_p(1:nxg,0) = chi_g(:,1);     b_p(1:nxg,0) = b_g(:,1)
    chi_p(1:nxg,nyg+1) = chi_g(:,nyg); b_p(1:nxg,nyg+1) = b_g(:,nyg)
    ! The edge columns are copied after the rows, so the corners take the corner cells
    chi_p(0,:) = chi_p(1,:);  chi_p(nxg+1,:) = chi_p(nxg,:)
    b_p(0,:) = b_p(1,:);      b_p(nxg+1,:) = b_p(nxg,:)
    ! I_ghost = B_edge at a thermal end or edge, and 0 at an open one
    i_p = 0.0d0
    if(bc_lo_thermal) i_p(1:nxg,0) = b_g(:,1)
    if(bc_hi_thermal) i_p(1:nxg,nyg+1) = b_g(:,nyg)
    if(.not.bc_h_periodic) then
      do j = 1, nyg
        if(bc_hlo_thermal) i_p(0,j) = b_g(1,j)
        if(bc_hhi_thermal) i_p(nxg+1,j) = b_g(nxg,j)
      end do
    end if
  end subroutine sc_pad_grid

  !> One sweep of the whole replicated grid for direction (nx, ny) into ii_g. The padded grid of
  !> sc_pad_grid is swept as a single block, which wraps in x when the domain is periodic. The wrap
  !> is not exact: a ray swept column by column reaches the seam before its upwind column holds
  !> this ray, so the seam takes a vacuum inflow for it. This is accepted, since the gather is only
  !> a comparison path; the KBA solvers wrap with the intensity of the previous solve
  subroutine sc_sweep_2d(nx, ny)
    double precision, intent(in) :: nx, ny
    ! The interior starts from zero for every ray, while the ghosts keep the inflows
    i_p(1:nxg,1:nyg) = 0.0d0
    call sc_sweep_block(nx, ny, dxg, dyg, .true., 0,0,nxg+1,nyg+1, 1,1,nxg,nyg, chi_p, b_p, i_p)
    ii_g = i_p(1:nxg,1:nyg)
  end subroutine sc_sweep_2d

  !> Textbook Olson-Kunasz weights for linear short characteristics over an optical-depth step dtau
  subroutine ka_coeffs(dtau, ex, psi0, psiu)
    double precision, intent(in)  :: dtau
    double precision, intent(out) :: ex, psi0, psiu
    double precision :: w0, w1
    ! A step this thin is transparent, so the intensity passes through unchanged
    if(dtau <= 1.0d-10) then
      ex = 1.0d0; psi0 = 0.0d0; psiu = 0.0d0; return
    end if
    ! The upwind intensity is attenuated by exp(-dtau). The source contributes the integrals of
    ! e^-t and t e^-t over the step, which we evaluate in closed form unless dtau is so small that
    ! the closed forms cancel, in which case their Taylor series are used
    ! I = I_u exp(-dtau) + psi_u S_u + psi_0 S_0, with
    ! w0 = int_0^dtau exp(-t) dt = 1 - exp(-dtau)
    ! w1 = int_0^dtau t exp(-t) dt = w0 - dtau exp(-dtau)
    ex = dexp(-dtau)
    if(dtau > 1.0d-4) then
      w0 = 1.0d0 - ex
      w1 = w0 - dtau*ex
    else
      ! w0 = dtau - dtau^2/2 + dtau^3/6,  w1 = dtau^2/2 - dtau^3/3
      w0 = dtau - dtau*dtau*0.5d0 + dtau*dtau*dtau/6.0d0
      w1 = dtau*dtau*0.5d0 - dtau*dtau*dtau/3.0d0
    end if
    ! With the source varying linearly in tau, w0 divides into a weight on the upwind source and a
    ! weight on the local one
    ! psi_u = w1/dtau,  psi_0 = w0 - psi_u
    psiu = w1/dtau
    psi0 = w0 - psiu
  end subroutine ka_coeffs

  !> Sweep of one block for ray (nx, ny), the kernel of the wavefront solvers. The upwind face and
  !> corner values sit in the ghost ring that getbc fills. With xwrap the block spans the whole
  !> periodic x-domain and wraps iu in place. Otherwise the x-upwind values come from the ghosts,
  !> and the block at the periodic seam reads the wrap ghost left by the previous solve
  subroutine sc_sweep_block(nx, ny, dxb, dyb, xwrap, ixI^L, ixO^L, chi, bpl, ii)
    use mod_global_parameters
    double precision, intent(in)    :: nx, ny, dxb, dyb
    logical, intent(in)             :: xwrap
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: chi(ixI^S), bpl(ixI^S)
    double precision, intent(inout) :: ii(ixI^S)
    integer :: i, j, sx, sy, nxper
    double precision :: absnx, absny, ds, frac
    absnx = dabs(nx); absny = dabs(ny)
    sx = ray_sign(nx); sy = ray_sign(ny)
    ! The block width in cells, needed when the block wraps onto itself in x
    nxper = ixOmax1 - ixOmin1 + 1
    ! The upwind point lies on the face the ray crosses first, the x-face when
    ! abs(mu_x) dy >= abs(mu_y) dx. Its first neighbour A is then the x-neighbour (iu, j) and the
    ! columns are swept in turn, and otherwise A is the y-neighbour (i, ju) and the rows are
    ! swept in turn. The second neighbour B is always the diagonal (iu, ju)
    if(absnx*dyb >= absny*dxb) then
      ! ds = dx/abs(mu_x),  frac = (dx/dy) abs(mu_y)/abs(mu_x)
      ds = dxb/max(absnx,smalldouble);  frac = (dxb/dyb)*absny/max(absnx,smalldouble)
      if(sx > 0) then
        do i = ixOmin1, ixOmax1; call sweep_column(i); end do
      else
        do i = ixOmax1, ixOmin1, -1; call sweep_column(i); end do
      end if
    else
      ! ds = dy/abs(mu_y),  frac = (dy/dx) abs(mu_x)/abs(mu_y)
      ds = dyb/max(absny,smalldouble);  frac = (dyb/dxb)*absnx/max(absny,smalldouble)
      if(sy > 0) then
        do j = ixOmin2, ixOmax2; call sweep_row(j); end do
      else
        do j = ixOmax2, ixOmin2, -1; call sweep_row(j); end do
      end if
    end if
  contains
    !> x-upwind index of column i, wrapped in place when the block spans the periodic domain
    integer function upwind_x(i) result(iu)
      integer, intent(in) :: i
      iu = i - sx
      if(xwrap .and. bc_h_periodic) then
        if(iu<ixOmin1) iu=iu+nxper
        if(iu>ixOmax1) iu=iu-nxper
      end if
    end function upwind_x
    subroutine sweep_column(i)
      integer, intent(in) :: i
      integer :: j, iu
      iu = upwind_x(i)
      do j = ixOmin2, ixOmax2
        call step(i, j, iu, j, iu, j-sy)
      end do
    end subroutine sweep_column
    subroutine sweep_row(j)
      integer, intent(in) :: j
      integer :: i
      do i = ixOmin1, ixOmax1
        call step(i, j, i, j-sy, upwind_x(i), j-sy)
      end do
    end subroutine sweep_row
    !> One cell (i, j) from its upwind neighbours A = (ia, ja) and B = (iu, ju). The upwind
    !> values may sit in ghost cells, where the driver has already written the band-end and edge
    !> inflows. chi and B do not depend on the ray, so their ghosts are always current
    subroutine step(i, j, ia, ja, iu, ju)
      integer, intent(in) :: i, j, ia, ja, iu, ju
      double precision :: in_a, in_b, iu_val, chiu, bu, dtau, ex, psi0, psiu
      in_a = ii(ia,ja)
      ! The intensity in the diagonal ghost is less reliable. In the pipelined solver it may
      ! already belong to the next ray, while both face ghosts still hold the current. Naturally
      ! this is problematic so we avoid it by reconstructing it from the adjacent faces instead.
      ! Seems to be a reasonable approximation at least
      if( (iu<ixOmin1 .or. iu>ixOmax1) .and. (ju<ixOmin2 .or. ju>ixOmax2) ) then
        ! I_B = (I(iu,j) + I(i,ju))/2
        in_b = 0.5d0*(ii(iu,j) + ii(i,ju))
      else
        in_b = ii(iu,ju)
      end if
      ! I_u = (1 - frac) I_A + frac I_B, and chi_u and B_u alike
      iu_val = (1.0d0-frac)*in_a + frac*in_b
      chiu   = (1.0d0-frac)*chi(ia,ja) + frac*chi(iu,ju)
      bu     = (1.0d0-frac)*bpl(ia,ja) + frac*bpl(iu,ju)
      ! dtau = (chi_u + chi) ds/2,  I = I_u exp(-dtau) + B psi_0 + B_u psi_u
      dtau = 0.5d0*(chiu + chi(i,j))*ds
      call ka_coeffs(dtau, ex, psi0, psiu)
      ii(i,j) = iu_val*ex + bpl(i,j)*psi0 + bu*psiu
    end subroutine step
  end subroutine sc_sweep_block
  }

end module mod_rt_sc_kernel
