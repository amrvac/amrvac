!=============================================================================
! Richtmyer-Meshkov hydro (2D)
! RM through planar shock impinging on inclined density discontinuity
module mod_usr
  use mod_hd

contains

  subroutine usr_init()

    call set_coordinate_system("Cartesian_2D")

    usr_init_one_grid => initonegrid_usr
    usr_set_parameters => initglobaldata_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call hd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr

    hd_gamma=1.4d0

  end subroutine initglobaldata_usr

  ! initialize one grid
  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_dust
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    integer                         :: n
    double precision                :: xshock, xbound, rhopost, vpost, ppost, alpha
    double precision                :: Prat,alfa,c_pre
    double precision                :: M,eta
    logical, save                   :: first = .true.
    !----------------------------------------------------------------------------

    {^IFONED   call mpistop("This is a 2D HD problem: Richtmyer Meshkov")}
    {^IFTHREED call mpistop("This is a 2D HD problem: Richtmyer Meshkov")}

    ! === Shock and CD position === !
    xshock=0.4d0
    xbound=0.5d0 ! xbound must be bigger than xshock!
    ! inclination angle for CD
    alpha = (3.14159265d0/4.0d0)

    ! Mach number of planar shock
    M     = 2.0d0
    ! inclination angle for CD
    alpha = (3.14159265d0/4.0d0)
    ! density ratio across CD
    eta   = 3.0d0

    ! compute the RH related states across the planar shock
    Prat    = one/(one+(M**2-one)*two*hd_gamma/(hd_gamma+one))
    alfa    = (hd_gamma+one)/(hd_gamma-one)
    c_pre   = one ! pre shock sound speed
    rhopost = hd_gamma*(alfa+Prat)/(alfa*Prat+one)
    ppost   = one/Prat
    vpost   = c_pre*M*(one-(alfa*Prat+one)/(alfa+Prat))

    if (mype==0.and.first) then
      print *,'============================================================='
      print *,' HD Richtmyer Meshkov simulation '
      print *,'============================================================='
      print *,' Mach number of shock: ',M
      print *,' post-shock density: ',rhopost
      print *,' post-shock velocity:',vpost
      print *,' post-shock pressure:',ppost
      print *,' Density ratio at CD: ',eta
      print *,'============================================================='
      first=.false.
    endif

    where(x(ix^S,1)>xshock.and.(x(ix^S,1)>x(ix^S,2)/dtan(alpha)+xbound))
      ! pre shock region
      w(ix^S,rho_)=hd_gamma*eta
      w(ix^S,mom(1))=zero
      w(ix^S,mom(2))=zero
      w(ix^S,e_)=one/(hd_gamma-one)
    elsewhere(x(ix^S,1)>xshock.and.(x(ix^S,1)<=x(ix^S,2)/dtan(alpha)+xbound))
      ! pre shock region
      w(ix^S,rho_)=hd_gamma
      w(ix^S,mom(1))=zero
      w(ix^S,mom(2))=zero
      w(ix^S,e_)=one/(hd_gamma-one)
    elsewhere
      ! post shock region
      w(ix^S,rho_)= rhopost
      w(ix^S,mom(1)) = rhopost*vpost
      w(ix^S,mom(2)) = zero
      w(ix^S,e_)  = ppost/(hd_gamma-one)+0.5d0*rhopost*vpost**2
    endwhere

  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),wlocal(ixI^S,1:nw)
    double precision :: rho(ixI^S),gradrho(ixI^S),drho(ixI^S)
    double precision :: kk,kk0,grhomax,kk1
    integer :: idims

! Example: assuming nwauxio=3 at convert stage 

! first store temperature
    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call hd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)

! then compute schlieren density plot
rho(ixI^S)=wlocal(ixI^S,rho_)
gradrho(ixO^S)=zero
do idims=1,ndim
   select case(typegrad)
   case("central")
     call gradient(rho,ixI^L,ixO^L,idims,drho)
   case("limited")
     call gradientS(rho,ixI^L,ixO^L,idims,drho)
   end select
   gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2.0d0
enddo
gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
kk=5.0d0
kk0=0.001d0
kk1=1.0d0
grhomax=2000.0d0

w(ixO^S,nw+2)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))

w(ixO^S,nw+3)=dlog10(w(ixO^S,rho_))


  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames

    varnames='T Schlier logrho'
  end subroutine specialvarnames_output

end module mod_usr
