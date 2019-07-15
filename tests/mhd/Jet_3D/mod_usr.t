!> Jet setup, B field as in A&A 486, 663, 2008. 3D version
! here we do a non-relativistic version (like the NR jet in paper)
module mod_usr
  use mod_mhd

  implicit none

  ! Input values for 
  !    Jet radius and initial Z-reach: Rjet, Zjet
  !    density and B factors: rhojet, rhocloud, apar, alfapar, pjet, B0, Bc, Bazi, npower
  !    perturbation in v at inlet: perturb_v, random_v, nmodes, ampl
  double precision :: Rjet, Zjet, rhojet, rhocloud, apar, alfapar, pjet, B0, Bc, Bazi, npower
  double precision :: ampl
  logical :: perturb_v, random_v
  integer :: nmodes

  ! for incompressible velocity perturbations at inlet
  integer,dimension(:),allocatable:: seed
  integer::  seed_size
  integer, parameter :: maxmodes=1000
  real:: rand_real(1:maxmodes)
  double precision :: rand_ampl(1:maxmodes), rand_phase(1:maxmodes)

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ Rjet, Zjet, rhojet, rhocloud, apar, alfapar, pjet, B0, Bc,  &
             Bazi, npower, perturb_v, random_v, nmodes, ampl

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

    if(perturb_v.and.(.not.random_v)) then
       ! reset nmodes to 7
       nmodes=7
       if (mype==0) write(*,*) 'resetting nmodes to 7 for deterministic perturbation'
    endif
    if(nmodes>maxmodes) call mpistop("number of modes too large, reset nmodes or recompile")

  end subroutine usr_params_read

  subroutine usr_init()
    use mod_variables

    call usr_params_read(par_files)

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => Jet_init_one_grid
    usr_special_bc      => specialbound_usr
    usr_refine_grid     => specialrefine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call set_coordinate_system("Cartesian_3D")

    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    character(len=20) :: printsettingformat
    integer :: ix

    printsettingformat='(1x,A50,ES15.7,A7)'
    
    if(mype==0) then
      write(*,*) "Jet setup:"
      write(*,printsettingformat) "density in jet ",rhojet," input"
      write(*,printsettingformat) "density in cloud ",rhocloud," input"
      write(*,printsettingformat) "pressure jet ",pjet," input"
      write(*,printsettingformat) "cloud B strength ",Bc," input"
      write(*,printsettingformat) "Jet B strength ",B0," input"
      write(*,printsettingformat) "Jet azimuthal B strength ",Bazi," input"
      write(*,printsettingformat) "flow strength factor ",alfapar," input"
    end if

    if(mype==0) then
      write(*,*) "Deduced dimensionless values:"
      write(*,printsettingformat) "density ratio jet/cloud ",rhojet/rhocloud," output"
    end if

    if(perturb_v)then
       rand_phase(1:maxmodes)=zero
       rand_ampl(1:maxmodes)=ampl
       if(mype==0)then
          if(random_v)then
             ! random phases and amplitudes
             print *,'nmodes =',nmodes
             write(*,*)'seeding random number generator, on mype==',mype
             call random_seed(SIZE=seed_size)
             allocate(seed(seed_size))
             call random_seed(GET=seed(1:seed_size))
             call random_number(rand_real(1:nmodes))
             write(*,*)'random numbers are:',rand_real(1:nmodes)
             rand_phase(1:nmodes)=-dpi+two*dpi*dble(rand_real(1:nmodes))
             call random_number(rand_real(1:nmodes))
             write(*,*)'random numbers are:',rand_real(1:nmodes)
             rand_ampl(1:nmodes)=ampl*dble(rand_real(1:nmodes))
          else
             ! following random numbers were used in KH comparison SWIFF paper
             ! just reuse for amplitude and phases
             rand_phase(1)=+1.09475373257713d0;
             rand_phase(2)=0.0d0;
             rand_phase(3)=+1.54234225329453d0;
             rand_phase(4)=-1.22142446608280d0;
             rand_phase(5)=+2.35557626490341d0;
             rand_phase(6)=+1.06889617006585d0;
             rand_phase(7)=-2.56820189483630d0;
             rand_ampl(1)=-0.615919305672490d0;
             rand_ampl(2)=0.0d0;
             rand_ampl(3)=-1.58333078691374d0;
             rand_ampl(4)=+2.96760094519789d0;
             rand_ampl(5)=-2.77959899737684d0;
             rand_ampl(6)=-2.41038627236035d0;
             rand_ampl(7)=-2.49807600476664d0;
             rand_ampl(1:7)=ampl*(rand_ampl(1:7)+dpi)/(two*dpi)
          endif
          open(123,file='phaseinfo',form='formatted')
          write(123,*) nmodes
          do ix=1,nmodes
              write(123,"(i4,2es12.4)") ix,rand_phase(ix),rand_ampl(ix)
          enddo
          close(123)
       endif
       call MPI_BARRIER(icomm,ierrmpi)
       if(npe>1)then
           call MPI_BCAST(rand_phase,maxmodes,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
           call MPI_BCAST(rand_ampl,maxmodes,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
       endif
    endif

  end subroutine initglobaldata_usr

  !> Initialize one grid
  subroutine Jet_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    
    double precision :: R(ixG^S),Z(ixG^S),hlpphi(ixG^S),hlpR(ixG^S)
    double precision :: cosphi(ixG^S),sinphi(ixG^S),scale,Bphi(ixG^S)

    ! assume (0,0) in middle of x-y domain
    R(ix^S)=dsqrt(x(ix^S,1)**2+x(ix^S,2)**2)
    ! axial direction Z/Zj
    Z(ix^S)=x(ix^S,3)/Zjet
    cosphi(ix^S)=x(ix^S,1)/R(ix^S)
    sinphi(ix^S)=x(ix^S,2)/R(ix^S)
    ! radial direction R/Rj
    R(ix^S)=R(ix^S)/Rjet

    where(R(ix^S)<one .and. Z(ix^S)<one)
       w(ix^S,rho_) = rhojet
       ! Bphi in jet
       Bphi(ix^S)= Bazi*TANH(R(ix^S)*Rjet/apar)
    elsewhere
       w(ix^S,rho_) = rhocloud
       Bphi(ix^S)= 0.0d0
    end where

    ! x-velocity vR*cosphi-vphi*sinphi, with vR=0
    w(ix^S,mom(1))  = -sinphi(ix^S)*Bphi(ix^S)/dsqrt(w(ix^S,rho_))
    ! y-velocity vR*sinphi+vphi*cosphi, with vR=0
    w(ix^S,mom(2))  =  cosphi(ix^S)*Bphi(ix^S)/dsqrt(w(ix^S,rho_))

    hlpphi(ix^S)=cosh(Z(ix^S)**npower)
    hlpR(ix^S)=(cosh(R(ix^S)**two))**2.0d0
    ! axial field BZ 
    w(ix^S,mag(3))= B0/(hlpphi(ix^S)*hlpR(ix^S)) +  Bc

    scale=npower*B0*Rjet/Zjet/2.0d0
    hlpR(ix^S)=scale*((Z(ix^S)**(npower-1.0d0))*TANH(Z(ix^S)**npower)*TANH(R(ix^S)**2.0d0))/&
            (hlpphi(ix^S)*R(ix^S))
    ! x-Bfield BR*cosphi-Bphi*sinphi, with BR=hlpR
    w(ix^S,mag(1))= hlpR(ix^S)*cosphi(ix^S)-Bphi(ix^S)*sinphi(ix^S)
    ! y-Bfield BR*sinphi+Bphi*cosphi, with BR=hlpR
    w(ix^S,mag(2))= hlpR(ix^S)*sinphi(ix^S)+Bphi(ix^S)*cosphi(ix^S)

    ! axial velocity vZ
    w(ix^S,mom(3))=alfapar*Bphi(ix^S) &
             /(dsqrt(w(ix^S,rho_))*R(ix^S)*Rjet/apar)
    w(ix^S,e_)= pjet+0.5d0*(B0**2.0d0-(Bphi(ix^S)**2.0d0 + w(ix^S,mag(3))**2.0d0))

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine Jet_init_one_grid

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: R(ixI^S),Z(ixI^S),hlpphi(ixI^S),hlpR(ixI^S)
    double precision :: cosphi(ixI^S),sinphi(ixI^S),scale,Bphi(ixI^S)
    integer :: ix1, ix2, ix3, ixOInt^L, imode, idims

    double precision :: psi(ixI^S),tmp(ixI^S)


    select case(iB)
    case(1)
      ! extrapolate all primitives continuously
      ! switch internal zone above boundary zone to primitive variables
      ixOInt^L=ixO^L;
      ixOIntmin1=ixOmax1+1
      ixOIntmax1=ixOmax1+1
      call mhd_to_primitive(ixI^L,ixOInt^L,w,x)
      do ix1 = ixOmin1,ixOmax1
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)  = w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))= w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))= w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))= w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)    = w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))= w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))= w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))= w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))
      enddo
      ! switch to conservative variables in internal zone
      call mhd_to_conserved(ixI^L,ixOInt^L,w,x)
      ! switch to conservative variables in ghost cells
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(2)
      ! extrapolate all primitives continuously
      ! switch internal zone above boundary zone to primitive variables
      ixOInt^L=ixO^L;
      ixOIntmin1=ixOmin1-1
      ixOIntmax1=ixOmin1-1
      call mhd_to_primitive(ixI^L,ixOInt^L,w,x)
      do ix1 = ixOmin1,ixOmax1
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)  = w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))= w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))= w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))= w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)    = w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))= w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))= w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))
         w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))= w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))
      enddo
      ! switch to conservative variables in internal zone
      call mhd_to_conserved(ixI^L,ixOInt^L,w,x)
      ! switch to conservative variables in ghost cells
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(3)
      ! extrapolate all primitives continuously
      ! switch internal zone above boundary zone to primitive variables
      ixOInt^L=ixO^L;
      ixOIntmin2=ixOmax2+1
      ixOIntmax2=ixOmax2+1
      call mhd_to_primitive(ixI^L,ixOInt^L,w,x)
      do ix2 = ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)  = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,rho_)
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(1))= w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,mom(1))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(2))= w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,mom(2))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(3))= w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,mom(3))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,e_)    = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,e_)
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1))= w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,mag(1))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(2))= w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,mag(2))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(3))= w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,mag(3))
      enddo
      ! switch to conservative variables in internal zone
      call mhd_to_conserved(ixI^L,ixOInt^L,w,x)
      ! switch to conservative variables in ghost cells
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ! extrapolate all primitives continuously
      ! switch internal zone above boundary zone to primitive variables
      ixOInt^L=ixO^L;
      ixOIntmin2=ixOmin2-1
      ixOIntmax2=ixOmin2-1
      call mhd_to_primitive(ixI^L,ixOInt^L,w,x)
      do ix2 = ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)  = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,rho_)
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(1))= w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,mom(1))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(2))= w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,mom(2))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(3))= w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,mom(3))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,e_)    = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,e_)
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1))= w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,mag(1))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(2))= w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,mag(2))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(3))= w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,mag(3))
      enddo
      ! switch to conservative variables in internal zone
      call mhd_to_conserved(ixI^L,ixOInt^L,w,x)
      ! switch to conservative variables in ghost cells
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(5)
      ! special bottom (Z=0) boundary
      ! fixing all profiles within Rjet
      ! (a)symmetry beyond Rjet
      ! assume (0,0) in middle of x-y domain
      ! range to ixI for later perturbation
      R(ixI^S)=dsqrt(x(ixI^S,1)**2+x(ixI^S,2)**2)
      cosphi(ixI^S)=x(ixI^S,1)/R(ixI^S)
      sinphi(ixI^S)=x(ixI^S,2)/R(ixI^S)
      ! axial direction Z/Zj
      Z(ixO^S)=x(ixO^S,3)/Zjet
      ! radial direction R/Rj
      R(ixO^S)=R(ixO^S)/Rjet
      ! switch internal zone above boundary zone to primitive variables
      ixOInt^L=ixO^L;
      ixOIntmin3=ixOmax3+1
      ixOIntmax3=ixOmax3+nghostcells
      call mhd_to_primitive(ixI^L,ixOInt^L,w,x)
      ! prescribe solution at jet inlet
      hlpphi(ixO^S)=cosh(Z(ixO^S)**npower)
      hlpR(ixO^S)=(cosh(R(ixO^S)**2.0d0))**2.0d0
      Bphi(ixO^S)= Bazi*TANH(R(ixO^S)*Rjet/apar)
      where(dabs(R(ixO^S)) <= 1.0d0)
         w(ixO^S,rho_)   = rhojet
         w(ixO^S,mag(3)) = B0/(hlpphi(ixO^S)*hlpR(ixO^S)) +  Bc
         ! x-velocity vR*cosphi-vphi*sinphi, with vR=0
         w(ixO^S,mom(1))  = -sinphi(ixO^S)*Bphi(ixO^S)/dsqrt(rhojet)
         ! y-velocity vR*sinphi+vphi*cosphi, with vR=0
         w(ixO^S,mom(2))  =  cosphi(ixO^S)*Bphi(ixO^S)/dsqrt(rhojet)
      endwhere

      if(perturb_v)then
        ! now add the incompressible perturbations
        psi(ixI^S)=zero
        ! this returns phi in range [-pi to pi] from x-axis
        tmp(ixI^S)=atan2(sinphi(ixI^S),cosphi(ixI^S))
        do imode=1,nmodes
           psi(ixI^S)=psi(ixI^S) &
             +rand_ampl(imode)*dcos(dble(imode)*tmp(ixI^S)+rand_phase(imode)) &
                   *dexp(-((dsqrt(x(ixI^S,1)**2+x(ixI^S,2)**2)-0.75d0*Rjet)/Rjet)**2) 
        enddo
        ! compute dpsi/dy
        idims=2
        select case(typegrad)
          case("central")
            call gradient(psi,ixI^L,ixO^L,idims,tmp)
          case("limited")
            call gradientS(psi,ixI^L,ixO^L,idims,tmp)
        end select
        w(ixO^S,mom(1))=w(ixO^S,mom(1))-tmp(ixO^S)
        ! compute dpsi/dx
        idims=1
        select case(typegrad)
          case("central")
            call gradient(psi,ixI^L,ixO^L,idims,tmp)
          case("limited")
            call gradientS(psi,ixI^L,ixO^L,idims,tmp)
        end select
        w(ixO^S,mom(2))=w(ixO^S,mom(2))+tmp(ixO^S)
      endif

      scale=npower*B0*Rjet/Zjet/2.0d0
      hlpR(ixO^S)=scale*((Z(ixO^S)**(npower-1.0d0))*TANH(Z(ixO^S)**npower)*TANH(R(ixO^S)**2.0d0))/&
            (hlpphi(ixO^S)*R(ixO^S))

      where(dabs(R(ixO^S)) <= 1.0d0)
         w(ixO^S,mag(1))= hlpR(ixO^S)*cosphi(ixO^S)-Bphi(ixO^S)*sinphi(ixO^S)
         w(ixO^S,mag(2))= hlpR(ixO^S)*sinphi(ixO^S)+Bphi(ixO^S)*cosphi(ixO^S)

         w(ixO^S,mom(3))=alfapar*Bphi(ixO^S) &
             /(dsqrt(w(ixO^S,rho_))*R(ixO^S)*Rjet/apar)
         w(ixO^S,e_)= pjet+0.5d0*(B0**2.0d0-(Bphi(ixO^S)**2.0d0 + w(ixO^S,mag(3))**2.0d0))
      endwhere

      ! extrapolate with reflection beyond jet radius
      do ix3 = ixOmin3,ixOmax3
         where(dabs(R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3)) > 1.0d0)
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,rho_)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2*ixOmax3-ix3+1,rho_)
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mom(1))= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2*ixOmax3-ix3+1,mom(1))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mom(2))= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2*ixOmax3-ix3+1,mom(2))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2*ixOmax3-ix3+1,mom(3))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,e_)    = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2*ixOmax3-ix3+1,e_)
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2*ixOmax3-ix3+1,mag(1))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(2))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2*ixOmax3-ix3+1,mag(2))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(3))= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2*ixOmax3-ix3+1,mag(3))
         endwhere
      enddo

      ! switch to conservative variables in internal zone
      call mhd_to_conserved(ixI^L,ixOInt^L,w,x)
      ! switch to conservative variables in ghost cells
      call mhd_to_conserved(ixI^L,ixO^L,w,x)

    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision:: R(ixG^S), Z(ixG^S)

    R(ix^S)=dsqrt(x(ix^S,1)**2+x(ix^S,2)**2)/Rjet
    Z(ix^S)=x(ix^S,3)/Zjet

    if (all((R(ix^S) <= 3.0d0) .and. (Z(ix^S)<= 3.0d0))) then
       if(level<4)then
          refine=1
          coarsen=-1
       endif
    endif

    if(any(R(ix^S) >= 8.0d0))then
       if(level>3) then
         refine=-1
         coarsen=1
       endif
    endif
    if(any(Z(ix^S) >= 20.0d0))then
       if(level>3) then
         refine=-1
         coarsen=1
       endif
    endif

  end subroutine specialrefine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_radiative_cooling
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),B2(ixI^S)
    double precision :: divb(ixI^S),wlocal(ixI^S,1:nw)
    double precision :: Btotal(ixI^S,1:ndir),curlvec(ixI^S,1:ndir)
    integer :: idirmin,idir

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! output temperature
    call mhd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)

    do idir=1,ndir
       Btotal(ixI^S,idir)=w(ixI^S,mag(idir))
    end do
    ! B^2
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)

    ! output Alfven Mach number V_z/(B/sqrt(rho))
    w(ixO^S,nw+2)=wlocal(ixO^S,mom(3))/dsqrt(B2(ixO^S)*w(ixO^S,rho_))

    ! output divB1
    call get_normalized_divb(wlocal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=divb(ixO^S)

    ! output reciprocal plasma beta B**2/(p*2)
    w(ixO^S,nw+4)=B2(ixO^S)/(pth(ixO^S)*two)

    ! store current
    !call get_current(wlocal,ixI^L,ixO^L,idirmin,curlvec)
    !do idir=1,ndir
    !  w(ixO^S,nw+4+idir)=curlvec(ixO^S,idir)
    !end do

    ! output Mach number V_z/c_s
    w(ixO^S,nw+5)=wlocal(ixO^S,mom(3))/dsqrt(mhd_gamma*pth(ixO^S)*w(ixO^S,rho_))

    ! output log10(rho)
    w(ixO^S,nw+6)=dlog10(w(ixO^S,rho_))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames

    varnames='Te Ma divB betar Mach log10rho'

  end subroutine specialvarnames_output

end module mod_usr
