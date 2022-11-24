!> Module for physical and numeric constants
!> Created: 01.09.2012 Oliver Porth (Physical constants)
module mod_constants

  implicit none
  public

  ! A very small real number (but above machine precision)
  double precision, parameter :: smalldouble = 1.D-12

  !> A very large real number
  double precision, parameter :: bigdouble = 1.D+99

  !> A very large integer
  integer, parameter :: biginteger = 10000000

  !> \todo Remove these
  double precision, parameter :: zero    = 0.0d0
  double precision, parameter :: one     = 1.0d0
  double precision, parameter :: two     = 2.0d0
  double precision, parameter :: half    = 0.5d0
  double precision, parameter :: quarter = 0.25d0
  double precision, parameter :: third   = 1.d0/3.0d0

  !> Pi
  double precision, parameter :: dpi=3.141592653589793238462643383279502884197169399375105d0

  !> Proton mass in cgs
  double precision, parameter :: mp_cgs = 1.672621777d-24 ! g

  !> Hydrogen mass in cgs
  double precision, parameter :: mH_cgs = 1.6733D-24 ! g

  !> Boltzmann constant in cgs
  double precision, parameter :: kB_cgs = 1.3806488d-16 ! erg K^-1

  !> Proton mass in SI
  double precision, parameter :: mp_SI = 1.672621777d-27 ! kg

  !> Boltzmann constant in SI
  double precision, parameter :: kB_SI = 1.3806488d-23 ! J K^-1

  !> Permeability in SI
  double precision, parameter :: miu0_SI = 1.2566370614d-6 ! H m^-1

  !> Speed of light in SI
  double precision, parameter :: c_SI = 2.99792458d8 ! m s^-1

  double precision, PARAMETER :: const_c     = 2.99792458d10   ! cm s^-1           ; Speed of light
  double precision, PARAMETER :: const_me    = 9.1093897d-28   ! g                 ; Electron mass
  double precision, PARAMETER :: const_mp    = 1.672621777d-24 ! g                 ; Proton mass
  double precision, PARAMETER :: const_e     = 4.8032068d-10   ! g^1/2 cm^3/2 s^-1 ; Electron charge
  double precision, PARAMETER :: const_LSun  = 3.828d33        ! erg s^-1          ; Solar luminosity
  double precision, PARAMETER :: const_MSun  = 1.98892d33      ! g                 ; Solar mass
  double precision, PARAMETER :: const_RSun  = 6.9570d10       ! cm                ; Solar radius
  double precision, PARAMETER :: const_kB    = 1.3806488d-16   ! erg K^-1          ; Boltzmann constant
  double precision, PARAMETER :: const_h     = 6.6260755d-27   ! erg s             ; Planck constant
  double precision, parameter :: const_pc    = 3.08567758128d18! cm                ; 1 parsec

  double precision, PARAMETER :: const_sigma  = 5.670367d-5 ! erg cm^-2 K^-4 s^-1 ; Stefan-Boltzman constant
  double precision, PARAMETER :: const_G      = 6.67191d-8  ! cm^3 g^-1 s^-2      ; Cavendish gravitational constant
  double precision, PARAMETER :: const_rad_a  = 7.5657e-15  ! erg cm^-3 K^-4      ; radiation density constant
  double precision, PARAMETER :: const_kappae = 0.34d0      ! cm^2 g^-1           ; free-electron scattering opacity

  ! Conversion factors:
  double precision, PARAMETER :: const_eV    = 1.6021772d-12   ! erg/eV            ; Electron volt
  double precision, PARAMETER :: const_Tera  = 1.d12           ! -                 ; Tera
  double precision, PARAMETER :: const_Peta  = 1.d15           ! -                 ; Peta
  double precision, PARAMETER :: const_years = 3.1536d7        ! s year^-1         ; seconds in a year

end module mod_constants
