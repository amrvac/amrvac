!> Magnetic bipolar field
module mod_usr
  use mod_mhd
  implicit none
  integer, save :: nx1,nx2,nxbc^D
  integer, parameter :: jmax=8000
  double precision, allocatable :: pbc(:),rbc(:)
  ! 1D solar atmosphere table for pressure, density, and height
  double precision :: pa(jmax),ra(jmax),ya(jmax)
  double precision :: usr_grav,SRadius,rhob,Tiso,dr,gzone,bQ0
  double precision :: q_para,d_para,L_para, charge1_x(3), charge2_x(3), charge1, charge2

contains

  !==============================================================================
  ! Purpose: to include global parameters, set user methods, set coordinate 
  !          system and activate physics module.
  !==============================================================================
  subroutine usr_init()

    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => specialsource
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_init_vector_potential=>initvecpot_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_set_electric_field => driven_electric_field
    usr_set_B0          => specialset_B0
    usr_set_J0          => specialset_J0

    call set_coordinate_system("Cartesian_3D")
    call mhd_activate()

  end subroutine usr_init

  !==============================================================================
  ! Purpose: to initialize user public parameters and reset global parameters.
  !          Input data are also read here.
  !==============================================================================
  subroutine initglobaldata_usr()
    double precision :: x0, y0, z0, r0, h, theta, L, Bperp
    integer :: ixp, nphalf

    !unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    !unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    !unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    !unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    !unit_time          = unit_length/unit_velocity                     ! 85.8746159942810 s 
    if(.not.mhd_energy) then
      ! bottom density
      rhob=2.d0
      ! isothermal uniform temperature
      Tiso= mhd_adiab
    end if

    bQ0=1.d-4/unit_pressure*unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    gzone=0.2d0
    ! cell size in 1D solar atmosphere table
    dr=(2.d0*gzone+xprobmax3-xprobmin3)/dble(jmax)
    usr_grav=-2.74d4*unit_length/unit_velocity**2  ! solar gravity
    SRadius=6.955d10/unit_length                   ! Solar radius

    q_para=7.d19/(unit_magneticfield*unit_length**2) ! strength and sign of magnetic charges
    d_para=1.d9/unit_length ! depth of magnetic charges
    L_para=1.5d9/unit_length ! half distance between magnetic charges

    charge1=-q_para
    charge1_x(1)=-L_para
    charge1_x(2)=0.d0
    charge1_x(3)=-d_para

    charge2=q_para
    charge2_x(1)=L_para
    charge2_x(2)=0.d0
    charge2_x(3)=-d_para

    if(mhd_energy) call inithdstatic

  end subroutine initglobaldata_usr

  !> initialize solar atmosphere table in a vertical line through the global domain
  subroutine inithdstatic
    use mod_global_parameters

    double precision :: Ta(jmax),gg(jmax)
    double precision :: rpho,Tpho,Ttop,htra,wtra,Ttr,Fc,k_para
    double precision :: res,pb,rhob,invT
    integer :: j,na,ibc,btlevel

    rpho=1.151d15/unit_numberdensity ! number density at the bottom relaxla
    Tpho=8.d3/unit_temperature ! temperature of chromosphere
    Ttop=1.5d6/unit_temperature ! estimated temperature in the top
    htra=2.d8/unit_length ! height of initial transition region
    wtra=2.d7/unit_length ! width of initial transition region 
    htra=0.2d0           ! height of initial transition region
    wtra=0.02d0          ! width of initial transition region 
    Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
    Fc=2.d5/unit_pressure/unit_velocity  ! constant thermal conduction flux
    ! Spitzer thermal conductivity with cgs units
    k_para=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 
    !! set T distribution with height
    do j=1,jmax
       ya(j)=(dble(j)-0.5d0)*dr-gzone
       if(ya(j)>htra) then
         Ta(j)=(3.5d0*Fc/k_para*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
       else
         Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra-0.027d0)/wtra)+1.d0)
       endif
       gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
    enddo
    !! solution of hydrostatic equation 
    ra(1)=rpho
    pa(1)=rpho*Tpho
    !invT=0.d0
    !do j=2,jmax
    !   invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
    !   pa(j)=pa(1)*dexp(invT*dr)
    !   ra(j)=pa(j)/Ta(j)
    !end do
    do j=2,jmax
       pa(j)=(pa(j-1)+dr*(gg(j)+gg(j-1))*ra(j-1)/4.d0)/(one-dr*(gg(j)+gg(j-1))/&
              Ta(j)/4.d0)
       ra(j)=pa(j)/Ta(j)
    end do
    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dr+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dr
    rhob=ra(na)+res/dr*(ra(na+1)-ra(na))
    pb=pa(na)+res/dr*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    btlevel=refine_max_level
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(3,btlevel)*(dble(nghostcells-ibc+1)-0.5d0))/dr+0.5d0)
      res=gzone-dx(3,btlevel)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dr
      rbc(ibc)=ra(na)+res/dr*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dr*(pa(na+1)-pa(na))
    end do
    
    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'maxTa',Ta(jmax)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

  end subroutine inithdstatic

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: A(ixI^S,1:ndim)
    double precision :: Bfr(ixI^S,1:ndir)
    double precision :: res
    integer :: ix^D,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
      write(*,*)'Stable solar atmosphere with a bipolar magnetic field'
    endif
      first=.false.
    endif

    if(B0field) then
      w(ixO^S,mag(:))=zero
    else 
      if(stagger_grid) then
        call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
        call mhd_face_to_center(ixO^L,block)
      else
        call bipolar_field(ixI^L,ixO^L,x,A,Bfr)
        w(ixO^S,mag(:))=Bfr(ixO^S,:)
      end if
    end if
    w(ixO^S,mom(:))=0.d0

    if(mhd_energy) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         na=floor((x(ix^D,3)-xprobmin3+gzone)/dr+0.5d0)
         res=x(ix^D,3)-xprobmin3+gzone-(dble(na)-0.5d0)*dr
         w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dr))/two*(ra(na+1)-ra(na))
         w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dr))/two*(pa(na+1)-pa(na))
      {end do\}
    else if(mhd_adiab/=0) then
      ! isothermal
      w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                    (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
    else
      ! zero beta
      w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    end if

    if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    ! vector potential
    double precision :: Avec1(ixI^S,1:ndim)

    call bipolar_field(ixI^L,ixC^L,xC,Avec1)

    if (idir==3) then
      A(ixC^S)=Avec1(ixC^S,3)
    else if(idir==2) then 
      A(ixC^S)=Avec1(ixC^S,2)
    else
      A(ixC^S)=Avec1(ixC^S,1)
    end if

  end subroutine initvecpot_usr

  subroutine bipolar_field(ixI^L,ixO^L,x,A,Bbp)

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    ! vector potential
    double precision, intent(out)   :: A(ixI^S,1:ndim)
    ! magnetic field
    double precision, optional, intent(out)   :: Bbp(ixI^S,1:ndir)

    double precision :: tmp(ixO^S),f1(ixO^S),f2(ixO^S)

    A(ixO^S,1)=0.d0
    f1(ixO^S)=(x(ixO^S,1)-charge1_x(1))/(sqrt((x(ixO^S,1)-charge1_x(1))**2+(x(ixO^S,2)-charge1_x(2))**2+&
              (x(ixO^S,3)-charge1_x(3))**2)*((x(ixO^S,2)-charge1_x(2))**2+(x(ixO^S,3)-charge1_x(3))**2))
    f2(ixO^S)=(x(ixO^S,1)-charge2_x(1))/(sqrt((x(ixO^S,1)-charge2_x(1))**2+(x(ixO^S,2)-charge2_x(2))**2+&
              (x(ixO^S,3)-charge2_x(3))**2)*((x(ixO^S,2)-charge2_x(2))**2+(x(ixO^S,3)-charge2_x(3))**2))
    A(ixO^S,2)=charge1*(x(ixO^S,3)-charge1_x(3))*f1(ixO^S)+charge2*(x(ixO^S,3)-charge2_x(3))*f2(ixO^S)
    A(ixO^S,3)=-charge1*(x(ixO^S,2)-charge1_x(2))*f1(ixO^S)-charge2*(x(ixO^S,2)-charge2_x(2))*f2(ixO^S)

    if(present(Bbp)) then
      tmp(ixO^S)=-sqrt((x(ixO^S,1)-charge1_x(1))**2+(x(ixO^S,2)-charge1_x(2))**2+(x(ixO^S,3)-charge1_x(3))**2)**3
      Bbp(ixO^S,1)=(x(ixO^S,1)-charge1_x(1))/tmp(ixO^S)
      Bbp(ixO^S,2)=(x(ixO^S,2)-charge1_x(2))/tmp(ixO^S)
      Bbp(ixO^S,3)=(x(ixO^S,3)-charge1_x(3))/tmp(ixO^S)
      tmp(ixO^S)=sqrt((x(ixO^S,1)-charge2_x(1))**2+(x(ixO^S,2)-charge2_x(2))**2+(x(ixO^S,3)-charge2_x(3))**2)**3
      Bbp(ixO^S,1)=Bbp(ixO^S,1)+(x(ixO^S,1)-charge2_x(1))/tmp(ixO^S)
      Bbp(ixO^S,2)=Bbp(ixO^S,2)+(x(ixO^S,2)-charge2_x(2))/tmp(ixO^S)
      Bbp(ixO^S,3)=Bbp(ixO^S,3)+(x(ixO^S,3)-charge2_x(3))/tmp(ixO^S)
      Bbp(ixO^S,:)=q_para*Bbp(ixO^S,:)
    end if

  end subroutine bipolar_field

  ! allow user to change inductive electric field, especially for boundary driven applications
  subroutine driven_electric_field(ixI^L,ixO^L,qt,qdt,fE,s)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt,qdt
    type(state)                        :: s
    double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)

    integer :: ixC^L

    ! fix Bz at bottom boundary
    if(s%is_physical_boundary(5)) then
      ixCmin^D=ixOmin^D-1;
      ixCmax^D=ixOmax^D;
      fE(nghostcells^%3ixC^S,1:2)=0.d0
    end if

  end subroutine driven_electric_field

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: tmp1(ixI^S),tmp2(ixI^S),pth(ixI^S),Qp(ixI^S)
    double precision :: xlen^D,dxb^D,startpos^D,coeffrho
    integer :: idir,ix^D,ixIM^L,ixOs^L,jxO^L

    if(mhd_glm) w(ixO^S,psi_)=0.d0
    select case(iB)
     case(1)
       w(ixO^S,mom(1))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))/&
                        w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(2))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))/&
                        w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(3))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))/&
                        w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix1=ixOsmax1,ixOsmin1,-1
             !block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
             !       (-block%ws(ix1+2^%1ixOs^S,idir)&
             !   +4.d0*block%ws(ix1+1^%1ixOs^S,idir))
             block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
                    ( block%ws(ix1+3^%1ixOs^S,idir)&
                -5.d0*block%ws(ix1+2^%1ixOs^S,idir)&
                +7.d0*block%ws(ix1+1^%1ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L-kr(1,^D);
         jxO^L=ixO^L+nghostcells*kr(1,^D);
         block%ws(ixOs^S,1)=zero
         do ix1=ixOsmax1,ixOsmin1,-1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix1^%1ixOs^S,1)=Qp(ix1+1^%1ixO^S)*block%dvolume(ix1+1^%1ixO^S)&
             /block%surfaceC(ix1^%1ixOs^S,1)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix1=ixOmax1,ixOmin1,-1
           w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                      (-w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                 +4.0d0*w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       end if
       if(mhd_energy) then
         ixIM^L=ixO^L;
         ixIMmin1=ixOmax1+1;ixIMmax1=ixOmax1+nghostcells;
         call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         w(ixO^S,rho_)=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
         w(ixO^S,p_)=pth(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
       else
         ! isothermal
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
       end if
       if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case(2)
       w(ixO^S,mom(1))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))/&
                        w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(2))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))/&
                        w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(3))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))/&
                        w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix1=ixOsmin1,ixOsmax1
             !block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
             !       (-block%ws(ix1-2^%1ixOs^S,idir)&
             !   +4.d0*block%ws(ix1-1^%1ixOs^S,idir))
             block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
                    ( block%ws(ix1-3^%1ixOs^S,idir)&
                -5.d0*block%ws(ix1-2^%1ixOs^S,idir)&
                +7.d0*block%ws(ix1-1^%1ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L;
         jxO^L=ixO^L-nghostcells*kr(1,^D);
         block%ws(ixOs^S,1)=zero
         do ix1=ixOsmin1,ixOsmax1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix1^%1ixOs^S,1)=-Qp(ix1^%1ixO^S)*block%dvolume(ix1^%1ixO^S)&
             /block%surfaceC(ix1^%1ixOs^S,1)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix1=ixOmin1,ixOmax1
           w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                      (-w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                 +4.0d0*w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       end if
       if(mhd_energy) then
         ixIM^L=ixO^L;
         ixIMmin1=ixOmin1-nghostcells;ixIMmax1=ixOmin1-1;
         call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         w(ixO^S,rho_)=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
         w(ixO^S,p_)=pth(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
       else
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
       end if
       if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case(3)
       w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(1))/&
                        w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(2))/&
                        w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(3))/&
                        w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix2=ixOsmax2,ixOsmin2,-1
             !block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
             !       (-block%ws(ix2+2^%2ixOs^S,idir)&
             !   +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
             block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
                    ( block%ws(ix2+3^%2ixOs^S,idir)&
                -5.d0*block%ws(ix2+2^%2ixOs^S,idir)&
                +7.d0*block%ws(ix2+1^%2ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L-kr(2,^D);
         jxO^L=ixO^L+nghostcells*kr(2,^D);
         block%ws(ixOs^S,2)=zero
         do ix2=ixOsmax2,ixOsmin2,-1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix2^%2ixOs^S,2)=Qp(ix2+1^%2ixO^S)*block%dvolume(ix2+1^%2ixO^S)&
             /block%surfaceC(ix2^%2ixOs^S,2)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix2=ixOmax2,ixOmin2,-1
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                      (-w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                 +4.0d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       end if
       if(mhd_energy) then
         ixIM^L=ixO^L;
         ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+nghostcells;
         call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         w(ixO^S,p_)=pth(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3)
         w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
       else
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
       end if
       if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case(4)
       w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(1))/&
                        w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(2))/&
                        w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(3))/&
                        w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix2=ixOsmin2,ixOsmax2
             !block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
             !       (-block%ws(ix2-2^%2ixOs^S,idir)&
             !   +4.d0*block%ws(ix2-1^%2ixOs^S,idir))
             block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
                    ( block%ws(ix2-3^%2ixOs^S,idir)&
                -5.d0*block%ws(ix2-2^%2ixOs^S,idir)&
                +7.d0*block%ws(ix2-1^%2ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L;
         jxO^L=ixO^L-nghostcells*kr(2,^D);
         block%ws(ixOs^S,2)=zero
         do ix2=ixOsmin2,ixOsmax2
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix2^%2ixOs^S,2)=-Qp(ix2^%2ixO^S)*block%dvolume(ix2^%2ixO^S)&
             /block%surfaceC(ix2^%2ixOs^S,2)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix2=ixOmin2,ixOmax2
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                       (-w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,mag(1):mag(3))&
                  +4.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       end if
       if(mhd_energy) then
         ixIM^L=ixO^L;
         ixIMmin2=ixOmin2-nghostcells;ixIMmax2=ixOmin2-1;
         call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         w(ixO^S,p_)=pth(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3)
         w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
       else
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
       end if
       if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case(5)
       w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))/&
                        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
       w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(2))/& 
                        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
       w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(3))/& 
                        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==3) cycle
             ixOsmax^D=ixOmax^D;
             ixOsmin^D=ixOmin^D-kr(^D,idir);
             do ix3=ixOsmax3,ixOsmin3,-1
               ! 3 one-sided equal gradient
               block%ws(ix3^%3ixOs^S,idir) = 1.d0/11.d0*&
                  ( -2.d0*block%ws(ix3+4^%3ixOs^S,idir)&
                   +11.d0*block%ws(ix3+3^%3ixOs^S,idir)&
                   -27.d0*block%ws(ix3+2^%3ixOs^S,idir)&
                   +29.d0*block%ws(ix3+1^%3ixOs^S,idir))
             end do
         end do
         ixOs^L=ixO^L-kr(3,^D);
         jxO^L=ixO^L+nghostcells*kr(3,^D);
         block%ws(ixOs^S,3)=zero
         do ix3=ixOsmax3,ixOsmin3,-1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix3^%3ixOs^S,3)=Qp(ix3+1^%3ixO^S)*block%dvolume(ix3+1^%3ixO^S)&
             /block%surfaceC(ix3^%3ixOs^S,3)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 4th order accuacy constant value extrapolation
         do ix3=ixOmax3,ixOmin3,-1
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/25.0d0)* &
                ( -3.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+4,mag(1):mag(3)) &
                 +16.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mag(1):mag(3)) &
                 -36.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(1):mag(3)) &
                 +48.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(1):mag(3)))
         enddo
       end if
       if(mhd_energy) then
         do ix3=ixOmin3,ixOmax3
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,rho_)=rbc(ix3)
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,p_)=pbc(ix3)
         enddo
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
       else
         w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                       (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
       end if
       if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case(6)
       w(ixO^S,mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(1))/&
                       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       w(ixO^S,mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(2))/&
                       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       w(ixO^S,mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(3))/&
                       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==3) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix3=ixOsmin3,ixOsmax3
             block%ws(ix3^%3ixOs^S,idir) = 1.d0/3.d0*&
                    ( block%ws(ix3-3^%3ixOs^S,idir)&
                -5.d0*block%ws(ix3-2^%3ixOs^S,idir)&
                +7.d0*block%ws(ix3-1^%3ixOs^S,idir))
             !block%ws(ix3^%3ixOs^S,idir) = 1.d0/11.d0*&
             !   (+2.d0*block%ws(ix3-3^%3ixOs^S,idir)&
             !    -9.d0*block%ws(ix3-2^%3ixOs^S,idir)&
             !   +18.d0*block%ws(ix3-1^%3ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L;
         jxO^L=ixO^L-nghostcells*kr(3,^D);
         block%ws(ixOs^S,3)=zero
         do ix3=ixOsmin3,ixOsmax3
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix3^%3ixOs^S,3)=-Qp(ix3^%3ixO^S)*block%dvolume(ix3^%3ixO^S)&
             /block%surfaceC(ix3^%3ixOs^S,3)
         end do
         call mhd_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix3=ixOmin3,ixOmax3
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                 (     -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-2,mag(1):mag(3)) &
                  +4.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-1,mag(1):mag(3)))
         enddo
       end if
       if(mhd_energy) then
         ixIM^L=ixO^L;
         ixIMmin3=ixOmin3-1;ixIMmax3=ixOmax3;
         call getggrav(tmp1,ixI^L,ixIM^L,x)
         ixIMmin3=ixOmin3-1;ixIMmax3=ixOmin3-1;
         call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
         pth(ixOmin3-1^%3ixO^S)=pth(ixOmin3-1^%3ixO^S)/w(ixOmin3-1^%3ixO^S,rho_)
         where(pth(ixOmin3-1^%3ixO^S)<1.618d0)
           pth(ixOmin3-1^%3ixO^S)=1.618d0
         end where
         tmp2=0.d0
         do ix3=ixOmin3,ixOmax3
           tmp2(ixOmin3^%3ixO^S)=tmp2(ixOmin3^%3ixO^S)+0.5d0*(tmp1(ix3-1^%3ixO^S)+tmp1(ix3^%3ixO^S))/pth(ixOmin3-1^%3ixO^S)
           w(ix3^%3ixO^S,rho_)=w(ixOmin3-1^%3ixO^S,rho_)*dexp(tmp2(ixOmin3^%3ixO^S)*dxlevel(3))
           w(ix3^%3ixO^S,p_)=w(ix3^%3ixO^S,rho_)*pth(ixOmin3-1^%3ixO^S)
         enddo
       else if(mhd_adiab==0) then
         ! zero beta
         w(ixO^S,rho_)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
       else
         coeffrho=usr_grav*SRadius**2/Tiso
         do ix3=ixOmin3,ixOmax3
           w(ix3^%3ixO^S,rho_)=w(ixOmin3-1^%3ixO^S,rho_)*dexp(coeffrho*(1.d0/(SRadius+&
           x(ixOmin3-1^%3ixO^S,3))-1.d0/(SRadius+x(ix3^%3ixO^S,3))))
         enddo
       end if
       if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
     case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,3)))**2
  end subroutine

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,3)=ggrid(ixO^S)
  end subroutine gravity

  subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    !! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    if(mhd_solve_eaux) w(ixO^S,eaux_)=w(ixO^S,eaux_)+qdt*bQgrid(ixO^S)
    !! add localized heating lQ
    !if(iprob==2) then
    !  call getlQ(lQgrid,ixI^L,ixO^L,qtC,wCT,x)
    !  w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
    !endif

  end subroutine specialsource

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
    !!calculate background heating bQ, Mok 2016 ApJ
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(out) :: bQgrid(ixI^S)

    double precision :: Bmag(ixI^S),bvec(ixI^S,1:ndir),cvec(ixI^S,1:ndir),tmp(ixI^S)
    double precision :: curlo,curhi
    integer :: idims,idir

    bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,3)/6.d0)
    !if(B0field) then
    !  Bmag(ixI^S)=dsqrt(sum((w(ixI^S,mag(:))+block%B0(ixI^S,:,0))**2,dim=ndim+1))
    !  do idir=1,ndir
    !    bvec(ixI^S,idir)=(w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0))/Bmag(ixI^S)
    !  end do
    !else
    !  Bmag(ixI^S)=dsqrt(sum(w(ixI^S,mag(:))**2,dim=ndim+1))
    !  do idir=1,ndir
    !    bvec(ixI^S,idir)=w(ixI^S,mag(idir))/Bmag(ixI^S)
    !  end do
    !endif
    !cvec=0.d0
    !! calculate local curvature of magnetic field
    !do idims=1,ndim
    !  call gradient(bvec(ixI^S,1),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,1)=cvec(ixO^S,1)+bvec(ixO^S,idims)*tmp(ixO^S)
    !  call gradient(bvec(ixI^S,2),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,2)=cvec(ixO^S,2)+bvec(ixO^S,idims)*tmp(ixO^S)
    !  call gradient(bvec(ixI^S,3),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,3)=cvec(ixO^S,3)+bvec(ixO^S,idims)*tmp(ixO^S)
    !end do 
    !tmp(ixO^S)=dsqrt(sum(cvec(ixO^S,:)**2,dim=ndim+1))
    !! set lower and upper limit for curvature
    !curlo=0.1d0 ! 1/10 
    !curhi=2.d0 ! 1/0.5
    !where(tmp(ixO^S)<curlo)
    !  tmp(ixO^S)=curlo
    !elsewhere(tmp(ixO^S)>curhi)
    !  tmp(ixO^S)=curhi
    !end where
    !
    !bQgrid(ixO^S)=bQ0*Bmag(ixO^S)**1.75d0*w(ixO^S,rho_)**0.125d0*tmp(ixO^S)**0.75d0
  
  end subroutine getbQ

  !==============================================================================
  ! Purpose: Enforce additional refinement or coarsening. One can use the
  !          coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !==============================================================================
  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (any(x(ixO^S,3)<=xprobmin3+0.2d0)) then
      refine=1
      coarsen=-1
    endif
  end subroutine special_refine_grid

  !==============================================================================
  ! Purpose: 
  !   this subroutine can be used in convert, to add auxiliary variables to the
  !   converted output file, for further analysis using tecplot, paraview, ....
  !   these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  !   the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  !   corresponding normalization values (default value 1)
  !==============================================================================
  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: tmp(ixI^S),dip(ixI^S),divb(ixI^S),B2(ixI^S)
    double precision, dimension(ixI^S,1:ndir) :: Btotal,qvec,curlvec
    integer                            :: ix^D,idirmin,idims,idir,jdir,kdir

    ! Btotal & B^2
    Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+1)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
    ! output divB1
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+2)=divb(ixO^S)
    ! output the plasma beta p*2/B**2
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+3)=2.d0*tmp(ixO^S)/B2(ixO^S)
    ! store current
    call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
    do idir=1,ndir
      w(ixO^S,nw+3+idir)=curlvec(ixO^S,idir)
    end do
    ! calculate Lorentz force
    qvec(ixO^S,1:ndir)=zero
    do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
      if(lvc(idir,jdir,kdir)/=0)then
        tmp(ixO^S)=curlvec(ixO^S,jdir)*w(ixO^S,mag(kdir))
        if(lvc(idir,jdir,kdir)==1)then
          qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
        else
          qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
        endif
      endif
    enddo; enddo; enddo
    do idir=1,ndir
      w(ixO^S,nw+3+ndir+idir)=qvec(ixO^S,idir)
    end do
    ! find magnetic dips
    !dip=0.d0
    !do idir=1,ndir
    !  call gradient(w(ixI^S,mag(3)),ixI^L,ixO^L,idir,tmp)
    !  dip(ixO^S)=dip(ixO^S)+w(ixO^S,b0_+idir)*tmp(ixO^S)
    !end do
    !where(dabs(w(ixO^S,mag(3)))<0.08d0 .and. dip(ixO^S)>=0.d0)
    !  w(ixO^S,nw+8)=1.d0
    !elsewhere
    !  w(ixO^S,nw+8)=0.d0
    !end where
  end subroutine specialvar_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    double precision :: A(ixI^S,1:ndim)

    call bipolar_field(ixI^L,ixO^L,x,A,wB0)

  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here add a time-independent background current density 
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)

    wJ0(ixO^S,:)=0.d0

  end subroutine specialset_J0

  !==============================================================================
  ! Purpose: names for special variable output
  !==============================================================================
  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Alfv divB beta j1 j2 j3 L1 L2 L3'
  end subroutine specialvarnames_output

end module mod_usr
