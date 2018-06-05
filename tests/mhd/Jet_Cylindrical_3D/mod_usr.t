!> Jet setup, B field as in A&A 486, 663, 2008. 3D version in cylindrical setup
! here we do a non-relativistic version (like the NR jet in paper)
module mod_usr
  use mod_mhd

  implicit none

  ! Input values for 
  !    Jet radius and initial Z-reach: Rjet, Zjet
  !    density and B factors: rhojet, rhocloud, apar, alfapar, pjet, B0, Bc, Bazi, npower
  double precision :: Rjet, Zjet, rhojet, rhocloud, apar, alfapar, pjet, B0, Bc, Bazi, npower

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ Rjet, Zjet, rhojet, rhocloud, apar, alfapar, pjet, B0, Bc, Bazi, npower

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

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

    call set_coordinate_system("cylindrical_3D")

    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    character(len=20) :: printsettingformat

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

  end subroutine initglobaldata_usr

  !> Initialize one grid
  subroutine Jet_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    
    double precision :: R(ixG^S),Z(ixG^S),hlpphi(ixG^S),hlpR(ixG^S),scale

    ! radial direction R/Rj
    R(ix^S)=x(ix^S,1)/Rjet
    ! axial direction Z/Zj
    Z(ix^S)=x(ix^S,2)/Zjet

    where(R(ix^S)<one .and. Z(ix^S)<one)
       w(ix^S,rho_) = rhojet
       ! Bphi in jet
       w(ix^S,mag(3))= Bazi*TANH(R(ix^S)*Rjet/apar)
    elsewhere
       w(ix^S,rho_) = rhocloud
       w(ix^S,mag(3))= 0.0d0
    end where

    ! radial velocity vR
    w(ix^S,mom(1))  = 0.0d0

    hlpphi(ix^S)=cosh(Z(ix^S)**npower)
    hlpR(ix^S)=(cosh(R(ix^S)**two))**2.0d0
    ! axial field BZ 
    w(ix^S,mag(2))= B0/(hlpphi(ix^S)*hlpR(ix^S)) +  Bc

    scale=npower*B0*Rjet/Zjet/2.0d0
    hlpR(ix^S)=scale*((Z(ix^S)**(npower-1.0d0))*TANH(Z(ix^S)**npower)*TANH(R(ix^S)**2.0d0))/&
            (hlpphi(ix^S)*R(ix^S))
    ! radial field BR
    w(ix^S,mag(1))= hlpR(ix^S)

    ! axial velocity vZ
    w(ix^S,mom(2))=alfapar*w(ix^S,mag(3)) &
             /(dsqrt(w(ix^S,rho_))*R(ix^S)*Rjet/apar)
    ! azimuthal velocity vphi
    w(ix^S,mom(3))=w(ix^S,mag(3))/dsqrt(w(ix^S,rho_))
    w(ix^S,e_)= pjet+0.5d0*(B0**2.0d0-(w(ix^S,mag(3))**2.0d0 + w(ix^S,mag(2))**2.0d0))

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine Jet_init_one_grid

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: R(ixI^S),Z(ixI^S),hlpphi(ixI^S),hlpR(ixI^S),scale
    integer :: ix2, ixOInt^L

    R(ixO^S)=x(ixO^S,1)/Rjet
    Z(ixO^S)=x(ixO^S,2)/Zjet

    select case(iB)
    ! special bottom (Z=0) boundary
    ! fixing all profiles within Rjet
    ! (a)symmetry beyond Rjet
    case(3)
      ! switch internal zone above boundary zone to primitive variables
      ixOInt^L=ixO^L;
      ixOIntmin2=ixOmax2+1
      ixOIntmax2=ixOmax2+nghostcells
      call mhd_to_primitive(ixI^L,ixOInt^L,w,x)
      ! prescribe solution at jet inlet
      hlpphi(ixO^S)=cosh(Z(ixO^S)**npower)
      hlpR(ixO^S)=(cosh(R(ixO^S)**2.0d0))**2.0d0
      where(dabs(R(ixO^S)) <= 1.0d0)
         w(ixO^S,rho_)   = rhojet
         w(ixO^S,mag(3)) = Bazi*TANH(R(ixO^S)*Rjet/apar)
         w(ixO^S,mag(2)) = B0/(hlpphi(ixO^S)*hlpR(ixO^S)) +  Bc
         w(ixO^S,mom(1)) = zero
      endwhere

      scale=npower*B0*Rjet/Zjet/2.0d0
      hlpR(ixO^S)=scale*((Z(ixO^S)**(npower-1.0d0))*TANH(Z(ixO^S)**npower)*TANH(R(ixO^S)**2.0d0))/&
            (hlpphi(ixO^S)*R(ixO^S))
      where(dabs(R(ixO^S)) <= 1.0d0)
         w(ixO^S,mag(1))= hlpR(ixO^S)
         w(ixO^S,mom(2))=alfapar*w(ixO^S,mag(3)) &
             /(dsqrt(w(ixO^S,rho_))*R(ixO^S)*Rjet/apar)
         w(ixO^S,mom(3))=w(ixO^S,mag(3))/dsqrt(w(ixO^S,rho_))
         w(ixO^S,e_)= pjet+0.5d0*(B0**2.0d0-(w(ixO^S,mag(3))**2.0d0 + w(ixO^S,mag(2))**2.0d0))
      endwhere

      ! extrapolate with reflection beyond jet radius
      do ix2 = ixOmin2,ixOmax2
         where(dabs(R(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3)) > 1.0d0)
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)  = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,ixOmin3:ixOmax3,rho_)
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(1))= w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,ixOmin3:ixOmax3,mom(1))
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(2))=-w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,ixOmin3:ixOmax3,mom(2))
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(3))= w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,ixOmin3:ixOmax3,mom(3))
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,e_)    = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,ixOmin3:ixOmax3,e_)
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1))=-w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,ixOmin3:ixOmax3,mag(1))
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(2))= w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,ixOmin3:ixOmax3,mag(2))
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(3))=-w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,ixOmin3:ixOmax3,mag(3))
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

    R(ix^S)=x(ix^S,1)/Rjet
    Z(ix^S)=x(ix^S,2)/Zjet

    if (any((R(ix^S) <= 3.0d0) .and. (Z(ix^S)<= 3.0d0))) refine=1

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
    !call mhd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    !w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)

    !do idir=1,ndir
    !   Btotal(ixI^S,idir)=w(ixI^S,mag(idir))
    !end do
    ! B^2
    !B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)

    ! output Alfven Mach number V_z/(B/sqrt(rho))
    !w(ixO^S,nw+2)=wlocal(ixO^S,mom(2))/dsqrt(B2(ixO^S)*w(ixO^S,rho_))

    ! output divB1
    !call get_normalized_divb(wlocal,ixI^L,ixO^L,divb)
    !w(ixO^S,nw+3)=divb(ixO^S)

    ! output reciprocal plasma beta B**2/(p*2)
    !w(ixO^S,nw+4)=B2(ixO^S)/(pth(ixO^S)*two)

    ! store current
    call get_current(wlocal,ixI^L,ixO^L,idirmin,curlvec)
    do idir=1,ndir
      w(ixO^S,nw+idir)=curlvec(ixO^S,idir)
    end do

    ! output inverse pitch mu
    !w(ixO^S,nw+5)=Rjet*Btotal(ixO^S,3)/(x(ixO^S,1)*Btotal(ixO^S,2))

    ! output Mach number V_z/c_s
    !w(ixO^S,nw+6)=wlocal(ixO^S,mom(2))/dsqrt(mhd_gamma*pth(ixO^S)*w(ixO^S,rho_))

    ! output log10(rho)
    !w(ixO^S,nw+7)=dlog10(w(ixO^S,rho_))
    w(ixO^S,nw+4)=dlog10(w(ixO^S,rho_))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames

    !varnames='Te Ma divB betar mu Mach log10rho'
    varnames='j1 j2 j3 log10rho'

  end subroutine specialvarnames_output

end module mod_usr
