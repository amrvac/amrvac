module mod_constants
  ! module file for physical constants.
  ! 01.09.2012 Oliver Porth
  ! Sets constants in cgs units.
  ! For usage, the user has to provide the three scaling parameters
  ! unit_length, unit_density, unit_velocity
  ! such that l[cgs]=unit_length*l[code] and so on.  
  ! use mod_constants
  implicit none
  save
  
  double precision, PARAMETER:: const_c   = 2.99792458d10    ! cm s^-1           ; Speed of light
  double precision, PARAMETER:: const_me  = 9.1093897d-28    ! g                 ; Electron mass
  double precision, PARAMETER:: const_mp  = 1.672621777d-24  ! g                 ; Proton mass
  double precision, PARAMETER:: const_e   = 4.8032068d-10    ! g^1/2 cm^3/2 s^-1 ; Electron charge
  double precision, PARAMETER:: const_MSun= 1.98892d33       ! g                 ; Solar mass
  double precision, PARAMETER:: const_kB  = 1.3806488d-16    ! erg K^-1          ; Boltzmann constant
  double precision, PARAMETER:: const_h   = 6.6260755d-27    ! erg s             ; Planck constant
  ! Conversion factors:
  double precision, PARAMETER:: const_eV  = 1.6021772d-12    ! erg/eV            ; Electron volt
  double precision, PARAMETER:: const_Tera= 1.d12            ! -                 ; Tera
  double precision, PARAMETER:: const_Peta= 1.d15            ! -                 ; Peta
  double precision, PARAMETER:: const_years= 3.1536d7        ! s year^-1         ; seconds in a year

end module mod_constants
