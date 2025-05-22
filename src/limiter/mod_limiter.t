!> Module with slope/flux limiters
module mod_limiter
  implicit none
  public

  !> radius of the asymptotic region [0.001, 10], larger means more accurate in smooth
  !> region but more overshooting at discontinuities
  double precision :: cada3_radius
  double precision :: schmid_rad^D
  !$acc declare create(cada3_radius, schmid_rad^D)
  integer, parameter :: limiter_minmod = 1
  integer, parameter :: limiter_woodward = 2
  integer, parameter :: limiter_mcbeta = 3
  integer, parameter :: limiter_superbee = 4
  integer, parameter :: limiter_vanleer = 5
  integer, parameter :: limiter_albada = 6
  integer, parameter :: limiter_koren = 7
  integer, parameter :: limiter_cada = 8
  integer, parameter :: limiter_cada3 = 9
  integer, parameter :: limiter_schmid = 10
  integer, parameter :: limiter_venk = 11
  ! Special cases
  integer, parameter :: limiter_ppm = 12
  integer, parameter :: limiter_mp5 = 13
  integer, parameter :: limiter_weno3  = 14
  integer, parameter :: limiter_wenoyc3  = 15
  integer, parameter :: limiter_weno5 = 16
  integer, parameter :: limiter_weno5nm = 17
  integer, parameter :: limiter_wenoz5  = 18
  integer, parameter :: limiter_wenoz5nm = 19
  integer, parameter :: limiter_wenozp5  = 20
  integer, parameter :: limiter_wenozp5nm = 21
  integer, parameter :: limiter_weno5cu6 = 22
  integer, parameter :: limiter_teno5ad = 23
  integer, parameter :: limiter_weno7 = 24
  integer, parameter :: limiter_mpweno7 = 25

contains

  integer function limiter_type(namelim)
    use mod_comm_lib, only: mpistop
    character(len=*), intent(in) :: namelim

    select case (namelim)
    case ('minmod')
       limiter_type = limiter_minmod
    case ('woodward')
       limiter_type = limiter_woodward
    case ('mcbeta')
       limiter_type = limiter_mcbeta
    case ('superbee')
       limiter_type = limiter_superbee
    case ('vanleer')
       limiter_type = limiter_vanleer
    case ('albada')
       limiter_type = limiter_albada
    case ('koren')
       limiter_type = limiter_koren
    case ('cada')
       limiter_type = limiter_cada
    case ('cada3')
       limiter_type = limiter_cada3
    case ('schmid1')
       limiter_type = limiter_schmid
    case ('schmid2')
       limiter_type = limiter_schmid
    case('venk')
       limiter_type = limiter_venk
    case ('ppm')
       limiter_type = limiter_ppm
    case ('mp5')
       limiter_type = limiter_mp5
    case ('weno3')
       limiter_type = limiter_weno3
    case ('wenoyc3')
       limiter_type = limiter_wenoyc3
    case ('weno5')
       limiter_type = limiter_weno5
    case ('weno5nm')
       limiter_type = limiter_weno5nm
    case ('wenoz5')
       limiter_type = limiter_wenoz5
    case ('wenoz5nm')
       limiter_type = limiter_wenoz5nm
    case ('wenozp5')
       limiter_type = limiter_wenozp5
    case ('wenozp5nm')
       limiter_type = limiter_wenozp5nm
    case ('weno5cu6')
       limiter_type = limiter_weno5cu6
    case ('teno5ad')
       limiter_type = limiter_teno5ad
    case ('weno7')
       limiter_type = limiter_weno7
    case ('mpweno7')
       limiter_type = limiter_mpweno7

    case default
       limiter_type = -1
       write(*,*) 'Unknown limiter: ', namelim
       call mpistop("No such limiter")
    end select
  end function limiter_type

end module mod_limiter
