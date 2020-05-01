! Several bremsstrahlung cross sections are provided in this module
module mod_bremsstrahlung
  use mod_physics
  use mod_constants
  implicit none

  double precision, parameter :: r0=2.8179e-13,alpha=7.2992d-3
  double precision, parameter :: mec2=511.d0

contains

  subroutine bremcross_Kramers(Z,Ee,Eph,sigmaB)
    ! Kramers bremsstrahlung cross section
    ! Z:ion change; Ee:electron energy (keV); Eph:photon energy (keV)
    ! sigmaB: cross section (cm^-2 keV^-1)

    integer :: Z
    double precision :: Eph,Ee,sigmaB

    double precision :: E1,E2,p1,p2,k
    double precision :: sigma0

    sigmaB=0.d0

    if (Ee>Eph) then
      sigma0=(16.d0/3)*alpha*r0**2
      E1=Ee/mec2  ! initial electron kinetic energy in units of mec2
      k=Eph/mec2  ! photon energy in units of mec2
      E2=E1-k     ! kinetic energy of scattered electron
      p1=sqrt(E1*(2+E1))
      p2=sqrt(E2*(2+E2))
      sigmaB=sigma0*Z**2/(k*p1**2)
      sigmaB=sigmaB/mec2
    endif

  end subroutine bremcross_Kramers

  subroutine bremcross_BetheHeitler(Z,Ee,Eph,sigmaB)
    ! Bethe-Heitler bremsstrahlung cross section
    ! Formula 3BN in Koch, Motz, Reviews of Modern Physics 31, 920 (1959)
    ! Z:ion change; Ee:electron energy (keV); Eph:photon energy (keV)
    ! sigmaB: cross section (cm^-2 keV^-1)

    integer :: Z
    double precision :: Eph,Ee,sigmaB

    double precision :: E1,E2,p1,p2,k
    double precision :: sigma0

    sigmaB=0.d0

    if (Ee>Eph) then
      sigma0=(16.d0/3)*alpha*r0**2
      E1=Ee/mec2  ! initial electron kinetic energy in units of mec2
      k=Eph/mec2  ! photon energy in units of mec2
      E2=E1-k     ! kinetic energy of scattered electron
      p1=sqrt(E1*(2+E1))
      p2=sqrt(E2*(2+E2))
      sigmaB=sigma0*Z**2/(k*p1**2)*log((p1+p2)/(p1-p2))
      sigmaB=sigmaB/mec2
    endif

  end subroutine bremcross_BetheHeitler

  subroutine bremcross_Haug(Z,Ee,Eph,sigmaB)
    ! Haug bremsstrahlung cross section
    ! Equation (4) and (5) in Haug, Astron. Astrophys. 326, 417 (1997)
    ! Z:ion change; Ee:electron energy (keV); Eph:photon energy (keV)
    ! sigmaB: cross section (cm^-2 keV^-1)

    integer :: Z
    double precision :: Eph,Ee,sigmaB

    double precision :: E1,E2,Et1,Et2,p1,p2,k
    double precision :: sigma0
    double precision :: term1,term2,term3,a1,a2,fE

    sigmaB=0.d0

    if (Ee>Eph) then
      sigma0=2.d0*alpha*r0**2
      E1=Ee/mec2  ! initial electron kinetic energy in units of mec2
      k=Eph/mec2  ! photon energy in units of mec2
      E2=E1-k     ! kinetic energy of scattered electron
      p1=sqrt(E1*(2+E1))
      p2=sqrt(E2*(2+E2))
      Et1=E1+1.d0  ! initial electron total energy
      Et2=E2+1.d0  ! total energy of scattered electron
      
      term1=(4.d0/3)*Et1*Et2+k**2-(7.d0/15)*(k**2/(Et1*Et2))
      term1=term1-(11.d0/70)*(k**2)*(p1**2+p2**2)/(Et1*Et2)**4
      term2=2*log((Et1*Et2+p1*p2-1)/k)
      term3=1+(1/(Et1*Et2))+(7.d0/20)*(p1*2+p2*2)/(Et1*Et2)**3
      term3=term3+((9.d0/28)*k**2+(263.d0/210)*(p1**2)*(p2**2))/(Et1*Et2)**3
      term3=term3*(p1*p2)/(Et1*Et2)

      a1=alpha*Z*Et1/p1
      a2=alpha*Z*Et2/p2
      fE=(a2/a1)*(1.d0-exp(-2*dpi*a1))/(1.d0-exp(-2*dpi*a2))

      sigmaB=sigma0*Z**2/(k*p1**2)*term1*(term2-term3)*fE
      sigmaB=sigmaB/mec2
    endif

  end subroutine bremcross_Haug

end module mod_bremsstrahlung
