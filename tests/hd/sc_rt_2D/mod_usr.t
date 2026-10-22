!> Regression test of the 2D short-characteristics RT on a static gas. beam: a vacuum
!> lit through a slit in the bottom boundary, pure transport, with rays wrapping the periodic x
!> seam. slab: chi falling off exponentially with height and B = 1, which exercises the
!> thin/thick heating blend and the tau^2 taper. table: a stratified atmosphere in cgs with the
!> Rosseland table opacity, warming from 5e3 K at the bottom to 5e4 K at the top, so that the
!> table lookup and the transition switch at Tcap are both exercised. The heating is stored in
!> Qrad, not applied
module mod_usr
  use mod_hd
  use mod_eos, only: eos
  use mod_rt_sc_slab, only: i_qrad, i_frad, i_iout, rt_sc_opacity
  implicit none

  !> 'beam', 'slab' or 'table'
  character(len=8) :: rt_case = 'beam'
  !> slit centre and half-width along the horizontal axis
  double precision :: beam_x0 = 0.5d0, beam_hw = 0.03d0
  !> slab opacity at the bottom and its scale height
  double precision :: slab_chi0 = 1.0d3, slab_h = 0.05d0

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_2D")
    usr_set_parameters => rt_params
    usr_init_one_grid  => rt_init
    usr_rt_sc_chiB     => rt_chiB
    ! The test's own opacity is the default, which &rt_list may replace (the table case)
    rt_sc_opacity = 'usr'
    call hd_activate()
  end subroutine usr_init

  subroutine rt_params()
    use mod_global_parameters
    integer :: n
    namelist /usr_list/ rt_case, beam_x0, beam_hw, slab_chi0, slab_h
    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, usr_list, end=111)
111   close(unitpar)
    end do
  end subroutine rt_params

  !> chi and B on ixO, the stratification running along y
  subroutine rt_chiB(ixI^L, ixO^L, w, x, chi, bpl)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: chi(ixI^S), bpl(ixI^S)
    select case(rt_case)
    case('slab')
      chi(ixO^S) = slab_chi0*dexp(-(x(ixO^S,2) - xprobmin2)/slab_h)
      bpl(ixO^S) = 1.0d0
    case default
      ! chi = 0, so B enters only as the thermal inflow at the band bottom and only the slit shines
      chi(ixO^S) = 0.0d0
      bpl(ixO^S) = 0.0d0
      where(dabs(x(ixO^S,1) - beam_x0) <= beam_hw) bpl(ixO^S) = 1.0d0
    end select
  end subroutine rt_chiB

  subroutine rt_init(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: zeta(ixI^S)
    if(rt_case == 'table') then
      ! The density falls over a scale height of a fifth of the domain from 1e-6 g/cm^3, and the
      ! temperature rises linearly with the height zeta and varies by 20% across x
      zeta(ixO^S)   = (x(ixO^S,2)-xprobmin2)/(xprobmax2-xprobmin2)
      w(ixO^S,rho_) = 1.0d-6*dexp(-5.0d0*zeta(ixO^S))
      w(ixO^S,e_)   = w(ixO^S,rho_)*(5.0d3 + 4.5d4*zeta(ixO^S)) &
                      *(1.0d0 + 0.2d0*dsin(2.0d0*dpi*(x(ixO^S,1)-xprobmin1)/(xprobmax1-xprobmin1)))
    else
      w(ixO^S,rho_) = 1.0d0
      w(ixO^S,e_)   = 1.0d0
    end if
    w(ixO^S,mom(:)) = 0.0d0
    w(ixO^S,i_qrad) = 0.0d0
    w(ixO^S,i_frad) = 0.0d0
    if(i_iout > 0) w(ixO^S,i_iout) = 0.0d0
    call eos%to_conserved(ixI^L, ixO^L, w, x)
  end subroutine rt_init

end module mod_usr
