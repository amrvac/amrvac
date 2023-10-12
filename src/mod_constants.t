!> Module for physical and numeric constants
module mod_constants

  implicit none
  public

  ! smallest positive number
  double precision, parameter :: tinydouble = TINY(1.D0)

  ! A very small real number (but above machine precision)
  double precision, parameter :: smalldouble = 1.D-16

  !> A very large real number
  double precision, parameter :: bigdouble = 1.D+99

  !> A very large integer
  integer, parameter          :: biginteger = 10000000

  !> \todo Remove these
  double precision, parameter :: zero    = 0.0d0
  double precision, parameter :: one     = 1.0d0
  double precision, parameter :: two     = 2.0d0
  double precision, parameter :: half    = 0.5d0
  double precision, parameter :: quarter = 0.25d0
  double precision, parameter :: third   = 1.d0/3.0d0

  !> Pi
  double precision, parameter :: dpi = 3.141592653589793238462643383279502884197169399375105d0
  !double precision, parameter :: dpi = 4.0d0 * datan(1.0d0)

  !> Constants in SI
  double precision, parameter :: c_SI        = 2.99792458d8       ! m s^-1            ; Speed of light
  double precision, parameter :: mp_SI       = 1.672621777d-27    ! kg                ; Proton mass
  double precision, parameter :: kB_SI       = 1.3806488d-23      ! J K^-1            ; Boltzmann constant
  double precision, parameter :: miu0_SI     = 1.2566370614d-6    ! H m^-1            ; Permeability in SI
  double precision, parameter :: avo_SI      = 6.0221367d23       ! mol^-1            ; Avogadro constant

  !> Constants in CGS
  double precision, parameter :: c_cgs       = 2.99792458d10      ! cm s^-1           ; Speed of light
  double precision, parameter :: me_cgs      = 9.1093897d-28      ! g                 ; Electron mass
  double precision, parameter :: mp_cgs      = 1.672621777d-24    ! g                 ; Proton mass
  double precision, parameter :: mn_cgs      = 1.674927211d-24    ! g                 ; Neutron mass
  double precision, parameter :: mH_cgs      = 1.6733d-24         ! g                 ; Hydrogen mass
  double precision, parameter :: e_cgs       = 4.8032068d-10      ! g^1/2 cm^3/2 s^-1 ; Electron charge
  double precision, parameter :: Msun_cgs    = 1.98892d33         ! g                 ; Solar mass
  double precision, parameter :: kB_cgs      = 1.3806488d-16      ! erg K^-1          ; Boltzmann constant
  double precision, parameter :: h_cgs       = 6.6260755d-27      ! erg s             ; Planck constant
  double precision, parameter :: G_cgs       = 6.67428d-8         ! cm^3 g^-1 s^-2    ; Gravitational constant
  double precision, parameter :: amu_cgs     = 1.66053873d-24     ! g                 ; Atomic mass unit

  !> Constants in other unit
  double precision, parameter :: me_MeV          = 0.51099895000d0  ! MeV       ; Electron mass
  double precision, parameter :: mmu_mev         = 105.6583755d0    ! MeV       ; Muon mass
  double precision, parameter :: mp_MeV          = 938.27208816d0   ! MeV       ; Proton mass
  double precision, parameter :: mn_MeV          = 939.56542052d0   ! MeV       ; Neutron mass
  double precision, parameter :: mn_minus_mp_MeV = 1.293548d0       ! MeV       ; neutron proton mass difference
  double precision, parameter :: amu_MeV         = 931.49432d0      ! MeV       ; Atomic mass unit
  double precision, parameter :: mpion0_mev      = 134.9768d0       ! MeV       ; pion_0 mass
  double precision, parameter :: mpionplus_mev   = 139.5704d0       ! MeV       ; pion_+ mass
  double precision, parameter :: malp_mev        = 3727.3794066d0   ! MeV       ; mass of alpha particle
  double precision, parameter :: m2H_mev         = 1875.612943d0    ! MeV       ; 2H mass
  double precision, parameter :: m3H_mev         = 2808.921133d0    ! MeV       ; 3H mass
  double precision, parameter :: m3He_mev        = 2809.413506d0    ! MeV       ; 3He mass
  double precision, parameter :: kB_MeVoverK     = 8.61738568d-11   ! MeV K^-1  ; Boltzmann constant
  double precision, parameter :: h_MeVs          = 4.135667696d-21  ! MeV s     ; Planck constant
  double precision, parameter :: hbarc_MeVcm     = 1.97326966d-11   ! MeV cm    ; hbar * c
  double precision, parameter :: hc_MeVcm        = 1.23984198d-10   ! MeV cm    ; h * c

  !> Conversion factors:
  double precision, parameter :: eV_to_erg   = 1.60217733d-12     ! erg / eV
  double precision, parameter :: MeV_to_erg  = 1.60217733d-6      ! erg / MeV
  double precision, parameter :: erg_to_MeV  = 6.24150636d5       ! MeV / erg
  double precision, parameter :: sec_to_year = 3.1536d7           ! s year^-1  
  double precision, parameter :: T_MeV_to_K  = 1.1604447522806d10 ! K / MeV
  double precision, parameter :: const_Tera  = 1.d12              ! -                 
  double precision, parameter :: const_Peta  = 1.d15              ! -                
  double precision, parameter :: cm3_to_fm3  = 1.d39              ! -

  !> Conversion factors from CGS to code unit, where c=G=Msun=1
  double precision, parameter :: length_gf   = 6.77140812d-06     ! cm
  double precision, parameter :: time_gf     = 2.03001708d+05     ! s
  double precision, parameter :: mass_gf     = 5.02765209d-34     ! g

  !double precision, parameter :: rho_gf      = 1.61930347d-18
  !double precision, parameter :: press_gf    = 1.80171810d-39
  double precision, parameter :: rho_gf      = (G_cgs/c_cgs**2)**3*Msun_cgs**2
  double precision, parameter :: press_gf    = G_cgs**3*Msun_cgs**2/c_cgs**8
  double precision, parameter :: eps_gf      = 1.11265006d-21
  double precision, parameter :: energy_gf   = 5.59424238d-55
  double precision, parameter :: lum_gf      = 2.7556091d-60
  double precision, parameter :: MeV_gf      = MeV_to_erg * energy_gf
  double precision, parameter :: mag_gf      = 4.244674d-20/(4.0d0*dpi) ! G

  ! Planck constant in code unit (MeV s)
  double precision, parameter :: h_gf        = h_MeVs * MeV_gf * time_gf
  ! hc in code unit (MeV cm), which is the same as h
  double precision, parameter :: hc_gf       = hc_MeVcm * MeV_gf * length_gf
  ! Planck constant in MeV without Hz (in code unit)
  double precision, parameter :: h_MeV       = h_MeVs * time_gf
  ! Atomic mass unit, in code unit
  double precision, parameter :: amu         = amu_MeV * MeV_gf


end module mod_constants
