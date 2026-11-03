!> Gray atmosphere in radiative equilibrium: a static column of constant absorption coefficient chi0
!> with the exact source function of the discrete-ordinate Milne problem for the two-cosine Carlson
!> set, S = 3/4 a Teff^4 (tau + q(tau)), q(tau) = q_inf + q_l exp(-q_k tau). The transfer solve must
!> return a vertical flux sigma Teff^4 at every depth and no heating, and the run stops if any cell
!> departs from either by more than the discretisation error of this resolution. The periodic seam
!> lags by one solve, so the narrow column takes about 30 solves to converge
module mod_usr
  use mod_hd
  use mod_rt_sc_state, only: i_qrad, i_frad
  implicit none

  double precision, parameter :: chi0 = 1.0d0, Teff = 1.0d0
  !> The root of 1 = sum_i 2 a_i/(1 - k^2 mu_i^2) and the constants that make the inflow at tau = 0 zero
  double precision, parameter :: q_k = 1.963961029224368d0, q_inf = 0.7060753572161309d0
  double precision, parameter :: q_l = -0.12872509118744369d0
  double precision, parameter :: tol_flux = 1.0d-2, tol_heat = 2.5d-2

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_2D")
    unit_length        = 1.0d8
    unit_temperature   = 1.0d4
    unit_numberdensity = 1.0d17
    usr_init_one_grid  => gray_init
    usr_rt_sc_chiB     => gray_chiB
    usr_process_global => gray_check
    call hd_activate()
  end subroutine usr_init

  !> A uniform static gas, since the radiation here does not depend on it
  subroutine gray_init(ixI^L, ixO^L, w, x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    w(ixO^S,rho_)   = 1.0d0
    w(ixO^S,mom(:)) = 0.0d0
    w(ixO^S,e_)     = 1.5d0
    w(ixO^S,i_qrad) = 0.0d0
    w(ixO^S,i_frad) = 0.0d0
  end subroutine gray_init

  subroutine gray_chiB(ixI^L, ixO^L, w, x, chi, bpl)
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: chi(ixI^S), bpl(ixI^S)
    double precision :: tau(ixI^S)
    tau(ixO^S) = chi0*(xprobmax2 - x(ixO^S,2))
    chi(ixO^S) = chi0
    bpl(ixO^S) = arad_norm*0.75d0*Teff**4*(tau(ixO^S) + q_inf + q_l*dexp(-q_k*tau(ixO^S)))
  end subroutine gray_chiB

  !> The largest departure of any cell, at the last step. The flux unit is a c Teff^4/4
  subroutine gray_check(iit, qt)
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    double precision :: sig, err(2), err_all(2)
    integer :: iigrid, igrid
    if(iit < it_max-1) return
    sig = arad_norm*(const_c/unit_velocity)*Teff**4/4.0d0
    err = 0.0d0
    do iigrid = 1, igridstail
      igrid = igrids(iigrid)
      err(1) = max(err(1), maxval(dabs(ps(igrid)%w(ixM^T,i_frad)/sig - 1.0d0)))
      err(2) = max(err(2), maxval(dabs(ps(igrid)%w(ixM^T,i_qrad)))/(chi0*sig))
    end do
    call MPI_ALLREDUCE(err, err_all, 2, MPI_DOUBLE_PRECISION, MPI_MAX, icomm, ierrmpi)
    if(mype == 0) write(*,'(a,2es12.4)') ' gray atmosphere max errors in F and Q:', err_all
    if(err_all(1) > tol_flux .or. err_all(2) > tol_heat) call mpistop("gray atmosphere: F or Q off the exact solution")
  end subroutine gray_check

end module mod_usr
