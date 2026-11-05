#:mute
#:include "../../../src/mod_gpu_directives.fpp"
#:endmute

!> Uniform, sub-luminal Cartesian SRMHD flow on a spherical mesh through both
!> poles.
!>
!> A constant density, constant pressure gas moving with a constant Cartesian
!> velocity is an exact steady solution: the Lorentz factor depends only on
!> |v|, which a spherical rotation leaves invariant, so it stays a single
!> global constant however the velocity's components are written. Written in
!> spherical components the velocity itself is a non-trivial function of theta
!> and phi, so keeping this state uniform exercises every curvilinear flux and
!> every geometric source term.
!>
!> This is a pole case: the domain deliberately runs onto the singular axis,
!> where the ghost cells come from the block half a revolution away in phi
!> rather than from a boundary condition. Note that the volume averages in the
!> log are almost blind to whether that copy is right - a cell at the axis has
!> vanishing volume, the face on the axis has zero area, and the TVD limiter
!> clips what is left - so the real test here is check_pole_ghosts below, which
!> compares the ghost cells themselves against the exact solution.
!>
!> SRMHD's primitive mom(:) slot holds the spatial four-velocity u^i = lfac*v^i
!> rather than v^i itself (see mod_srmhd_templates.fpp's to_primitive /
!> to_conservative). lfac is a scalar, so it is invariant under the pole's
!> pi-rotation, and u^i therefore transforms exactly like an ordinary vector -
!> the same asymm/symm pattern hd's momentum takes at the pole applies here
!> unchanged, which is what src/io/mod_input_output.fpp's read_par_files relies
!> on when it builds the sign table from iw_mom alone.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> ambient state. These are deliberately plain variables rather than
  !> parameters: specialbound_usr runs on the device, and a parameter has no
  !> storage to copy there.
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
  !> the uniform *ordinary* Cartesian velocity, sub-luminal
  double precision :: v0(3) = [0.3d0, 0.15d0, -0.1d0]
  !> the uniform Cartesian magnetic field, in the same frame as v0
  double precision :: b0(3) = [0.4d0, -0.3d0, 0.5d0]
  ${GPU_DECLARE_COPYIN('rho0, p0, v0, b0')}$

contains

  subroutine usr_init()

    ! the coordinate system comes from `geometry` in &meshlist

    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_process_grid  => check_pole_ghosts

    call phys_activate()

  end subroutine usr_init

  !> spherical components at x of a uniform Cartesian vector a0.
  !> A subroutine rather than an array-valued function, because a function
  !> result needs a temporary that not every OpenACC compiler handles inside
  !> a device routine.
  pure subroutine uniform_vector(x, a0, v)
    ! Vector parallelization is disabled for OpenMP in specialbound_usr, so this routine should
    ! work at the vector level in that case
#ifdef _OPENMP
    ${GPU_ROUTINE_VECTOR()}$
#else
    ${GPU_ROUTINE_SEQ()}$
#endif
    double precision, intent(in)  :: x(1:ndim)
    !> the Cartesian vector
    double precision, intent(in)  :: a0(1:3)
    double precision, intent(out) :: v(1:3)
    double precision              :: sint, cost, sinp, cosp

    sint = sin(x(2)); cost = cos(x(2))
    sinp = sin(x(3)); cosp = cos(x(3))

    v(1) =  a0(1)*sint*cosp + a0(2)*sint*sinp + a0(3)*cost
    v(2) =  a0(1)*cost*cosp + a0(2)*cost*sinp - a0(3)*sint
    v(3) = -a0(1)*sinp      + a0(2)*cosp

  end subroutine uniform_vector

  !> The scalars of the uniform state, which are global constants because
  !> each is a function of rotation invariants only (|v0|, |b0|, v0.b0):
  !>
  !>   lfac0  = 1/sqrt(1 - |v0|^2)
  !>   mu0    = 1/(lfac0 h0), the srmhd auxiliary, h0 the gamma-law enthalpy
  !>   w_tot0 = rho0 h0 lfac0^2 + |b0|^2, so that S = w_tot0 v - (v.B) B
  !>   vdotb0 = v0.b0
  !>   tau0   = w_tot0 - ptot0 - D0, ptot0 = p0 + (|b0|^2/lfac0^2 + vdotb0^2)/2
  !>
  !> matching to_conservative in src/srmhd/mod_srmhd_templates.fpp.
  pure subroutine uniform_constants(lfac0, mu0, w_tot0, vdotb0, tau0)
#ifdef _OPENMP
    ${GPU_ROUTINE_VECTOR()}$
#else
    ${GPU_ROUTINE_SEQ()}$
#endif
    double precision, intent(out) :: lfac0, mu0, w_tot0, vdotb0, tau0
    double precision              :: h0, b_sqr, ptot0

    lfac0  = 1.0d0 / sqrt(1.0d0 - sum(v0(1:3)**2))
    h0     = 1.0d0 + srmhd_gamma / (srmhd_gamma - 1.0d0) * p0 / rho0
    mu0    = 1.0d0 / (lfac0 * h0)
    b_sqr  = sum(b0(1:3)**2)
    vdotb0 = sum(v0(1:3) * b0(1:3))
    w_tot0 = rho0 * h0 * lfac0**2 + b_sqr
    ptot0  = p0 + 0.5d0 * (b_sqr / lfac0**2 + vdotb0**2)
    tau0   = w_tot0 - ptot0 - rho0 * lfac0

  end subroutine uniform_constants

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    ! .. local ..
    double precision                :: v(1:3), b(1:3), lfac0, x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    lfac0 = 1.0d0 / sqrt(1.0d0 - sum(v0(1:3)**2))

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call uniform_vector(x_loc, v0, v)
             call uniform_vector(x_loc, b0, b)
             w(ix1,ix2,ix3,rho_)   = rho0
             w(ix1,ix2,ix3,p_)     = p0
             w(ix1,ix2,ix3,mom(1)) = lfac0 * v(1)
             w(ix1,ix2,ix3,mom(2)) = lfac0 * v(2)
             w(ix1,ix2,ix3,mom(3)) = lfac0 * v(3)
             w(ix1,ix2,ix3,mag(1)) = b(1)
             w(ix1,ix2,ix3,mag(2)) = b(2)
             w(ix1,ix2,ix3,mag(3)) = b(3)
             w(ix1,ix2,ix3,psi_)   = 0.0d0
          end do
       end do
    end do

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  !> keep the exact solution in the ghost cells, so that only the interior
  !> discretisation is under test
  subroutine specialbound_usr(qt, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB, w, x)
    use mod_global_parameters
    ${GPU_ROUTINE_VECTOR()}$
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3
    integer, intent(in)             :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    integer, intent(in)             :: iB
    double precision, intent(in)    :: qt
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    ! .. local ..
    double precision                :: v(1:3), b(1:3), lfac0, mu0, w_tot0
    double precision                :: vdotb0, tau0
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    call uniform_constants(lfac0, mu0, w_tot0, vdotb0, tau0)

    ! collapse(3) is fine here even though the body calls uniform_vector:
    ! the nvfortran OpenACC miscompile that hit the hd pole cases needs BOTH a
    ! call AND a runtime-sized private automatic array in the collapsed body
    ! (hd's wpt(1:nw)). srmhd's privates (v, b, x_loc) are all compile-time sized,
    ! so there is nothing to miscompile. See CLAUDE.md ("Bug-hunting notes")
    ! and issue #154.
#ifndef _OPENMP
    ! Vector-level parallelization within a subroutine does not work with OpenMP.
    ! TBD if this routine can be made seq (i.e. cell-based).
    ${GPU_LOOP_VECTOR("collapse(3) private(v, b, x_loc)")}$
#endif
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call uniform_vector(x_loc, v0, v)
             call uniform_vector(x_loc, b0, b)

             ! conserved D = rho*lfac, S = W v - (v.B) B, tau, and the
             ! auxiliaries mu = 1/(lfac*h) and lfac; see uniform_constants
             w(ix1,ix2,ix3,iw_rho)    = rho0 * lfac0
             w(ix1,ix2,ix3,iw_mom(1)) = w_tot0 * v(1) - vdotb0 * b(1)
             w(ix1,ix2,ix3,iw_mom(2)) = w_tot0 * v(2) - vdotb0 * b(2)
             w(ix1,ix2,ix3,iw_mom(3)) = w_tot0 * v(3) - vdotb0 * b(3)
             w(ix1,ix2,ix3,iw_e)      = tau0
             w(ix1,ix2,ix3,iw_mag(1)) = b(1)
             w(ix1,ix2,ix3,iw_mag(2)) = b(2)
             w(ix1,ix2,ix3,iw_mag(3)) = b(3)
             w(ix1,ix2,ix3,psi_)      = 0.0d0
             w(ix1,ix2,ix3,mu_)       = mu0
             w(ix1,ix2,ix3,lfac_)     = lfac0

          end do
       end do
    end do

  end subroutine specialbound_usr

  !> Analytic conserved state at x, valid at any coordinates - the ghost cells
  !> beyond the axis included, whose stored theta is negative. That those
  !> still come out right is not a coincidence: the mirror and the
  !> per-component sign the pole copy applies are exactly the change of basis
  !> between a point and its partner half a revolution away, so writing the
  !> uniform Cartesian velocity in the ghost cell's own coordinates reproduces
  !> what the copy has to have put there - lfac0 itself does not change, since
  !> it depends only on |v0|.
  pure subroutine analytic_state(x_loc, wpt)
    double precision, intent(in)  :: x_loc(1:ndim)
    double precision, intent(out) :: wpt(1:nw)
    double precision              :: v(1:3), b(1:3), lfac0, mu0, w_tot0
    double precision              :: vdotb0, tau0

    call uniform_vector(x_loc, v0, v)
    call uniform_vector(x_loc, b0, b)
    call uniform_constants(lfac0, mu0, w_tot0, vdotb0, tau0)

    wpt(1:nw)      = 0.0d0
    wpt(iw_rho)    = rho0 * lfac0
    wpt(iw_mom(1)) = w_tot0 * v(1) - vdotb0 * b(1)
    wpt(iw_mom(2)) = w_tot0 * v(2) - vdotb0 * b(2)
    wpt(iw_mom(3)) = w_tot0 * v(3) - vdotb0 * b(3)
    wpt(iw_e)      = tau0
    wpt(iw_mag(1)) = b(1)
    wpt(iw_mag(2)) = b(2)
    wpt(iw_mag(3)) = b(3)
    wpt(mu_)       = mu0
    wpt(lfac_)     = lfac0

  end subroutine analytic_state

  !> Assert that the ghost cells across the axis carry the exact solution.
  !>
  !> This is the only sharp check of the pole copy the case can make.  The
  !> volume averages the log reports cannot see it: a cell at the axis has
  !> vanishing volume, and the face lying on the axis has zero area, so the
  !> flux that would carry a wrong ghost value into the interior is multiplied
  !> by zero, and the jump a wrong sign leaves behind is clipped by the TVD
  !> limiter.  Reversing every sign at the pole moves mean(w) only in the
  !> fifth digit.  The ghost cells themselves are wrong by O(1) at once.
  !>
  !> Only meaningful at it = 0, where the interior is still exactly the
  !> analytic state, so the pole copy of it has to reproduce the analytic
  !> state in the ghost cells to round-off.  Silence means the check passed.
  subroutine check_pole_ghosts(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    ! .. local ..
    double precision :: wpt(1:nw), x_loc(1:ndim), err
    integer          :: ix1, ix2, ix3, iside, i2, jxmin2, jxmax2

    if (it /= 0) return
    ! process() runs before the solution is pulled back for output, so the
    ! host copy of w is stale unless we fetch this block ourselves
    ${GPU_UPDATE_HOST('ps(igrid)%w')}$

    err = 0.0d0
    do iside = 1, 2
       i2 = 2*iside - 3
       if (neighbor_pole(0,i2,0,igrid) == 0) cycle
       if (iside == 1) then
          jxmin2 = ixImin2;   jxmax2 = ixOmin2-1
       else
          jxmin2 = ixOmax2+1; jxmax2 = ixImax2
       end if
       ! transverse loops run over the whole block (ixI), not just its
       ! interior (ixO), so the pole layer's edges and corners are covered -
       ! in particular the cells where the axis meets the radial boundary
       do ix3 = ixImin3, ixImax3
          do ix2 = jxmin2, jxmax2
             do ix1 = ixImin1, ixImax1
                x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
                call analytic_state(x_loc, wpt)
                err = max(err, maxval(abs(w(ix1,ix2,ix3,&
                   1:nwflux) - wpt(1:nwflux))))
             end do
          end do
       end do
    end do

    if (err > 1.0d-10) then
       write(*,*) 'pole ghost cells deviate from the exact solution by', err
       call mpistop('pole ghost-cell check failed')
    end if

  end subroutine check_pole_ghosts

end module mod_usr
