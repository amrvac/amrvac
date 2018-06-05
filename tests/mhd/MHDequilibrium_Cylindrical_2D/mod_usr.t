!> Test to check MHD equilibrium in Cylindrical 2D setting (no phi_ components)
! to check e.g. curlvector implementation
module mod_usr
  use mod_mhd

  implicit none

  ! Input values for 
  !    Jet radius and initial Z-reach: Rjet
  !    density and B factors: rhojet, rhocloud, pjet, B0, Bc
  double precision :: Rjet, rhojet, rhocloud, pjet, B0, Bc

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ Rjet, rhojet, rhocloud, pjet, B0, Bc

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
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    {^IFONED call mpistop("testing 2D here") }
    {^IFTHREED call mpistop("testing 2D here") }

    {^IFTWOD call set_coordinate_system("cylindrical_2D")}

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
    
    double precision :: R(ixG^S),hlpR(ixG^S)

    ! radial direction R/Rj
    R(ix^S)=x(ix^S,1)/Rjet

    where(R(ix^S)<one)
       w(ix^S,rho_) = rhojet
    elsewhere
       w(ix^S,rho_) = rhocloud
    end where

    ! radial velocity vR
    w(ix^S,mom(1))  = 0.0d0
    ! axial velocity vZ
    w(ix^S,mom(2))=0.0d0

    hlpR(ix^S)=(cosh(R(ix^S)**two))**2.0d0
    ! axial field BZ 
    w(ix^S,mag(2))= B0/hlpR(ix^S) +  Bc
    ! radial field BR
    w(ix^S,mag(1))= 0.0d0

    w(ix^S,e_)= pjet+0.5d0*(B0**2.0d0-(w(ix^S,mag(2))**2.0d0))

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine Jet_init_one_grid

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
    double precision :: Btotal(ixI^S,1:ndir),current(ixI^S,7-2*ndir:3)
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

    ! output divB1
    call get_normalized_divb(wlocal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+2)=divb(ixO^S)

    ! output reciprocal plasma beta B**2/(p*2)
    w(ixO^S,nw+3)=B2(ixO^S)/(pth(ixO^S)*two)

    ! store current
    call get_current(wlocal,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+4)=current(ixO^S,idirmin)

    ! output log10(rho)
    w(ixO^S,nw+5)=dlog10(w(ixO^S,rho_))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames

    varnames='Te divB betar j3 log10rho'

  end subroutine specialvarnames_output

end module mod_usr
