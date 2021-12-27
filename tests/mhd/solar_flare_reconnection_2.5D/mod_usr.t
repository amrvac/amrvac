module mod_usr
  use mod_mhd
  implicit none
  double precision :: q_e, parb,unit_currentdensity

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_2.5D")

    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    usr_init_one_grid       => initonegrid_usr
    usr_special_bc          => specialbound_usr
    usr_aux_output          => specialvar_output
    usr_add_aux_names       => specialvarnames_output 
    usr_set_B0              => specialset_B0
    usr_set_J0              => specialset_J0
    usr_special_convert     => usrspecial_convert
    usr_special_resistivity => special_eta
    usr_var_for_errest      => p_for_errest
    usr_init_vector_potential=>initvecpot_usr

    call mhd_activate()
    parb=20.d0/3.d0
    ! unit of current density
    unit_currentdensity=unit_magneticfield/unit_length/4.d0/dpi
    ! unit of charge
    q_e=unit_currentdensity/unit_numberdensity/unit_velocity
    if(mype==0) print*,'unit of charge',q_e
    ! dimensionless charge of electron
    q_e=1.60217653d-19/q_e
    if(mype==0) print*,'dimensionless e',q_e

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bf(ixI^S,1:ndir)
    double precision :: htra, wtra, rpho
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          print *,'YOKOYAMA and SHIBATA 2001 ApJ'
       end if
       first=.false.
    end if
    rpho=1.d5 ! number density at the bottom relaxla
    htra=0.3d0 ! height of initial transition region
    wtra=0.06d0 ! width of initial transition region 
    w(ixO^S,rho_)=1.d0+(rpho-1.d0)*(1.d0-tanh((x(ixO^S,2)-htra)/wtra))/2.d0
    w(ixO^S,p_)=1.d0
    w(ixO^S,mom(:))=zero
    if(B0field) then
      w(ixO^S,mag(:))=zero
    else if(stagger_grid) then
      call b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block%ws,x)
      call mhd_face_to_center(ixO^L,block)
      w(ixO^S,mag(3))=Busr/dcosh(parb*x(ixO^S,1))
    else
      call specialset_B0(ixI^L,ixO^L,x,Bf)
      w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
    end if
    if(mhd_glm) w(ixO^S,psi_)=0.d0
    if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
    call mhd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: pth(ixI^S),Qp(ixI^S)
    integer :: ix^D, ixA^L, ixOs^L, idir

    if(mhd_glm) w(ixO^S,psi_)=0.d0
    select case(iB)
    case(1)
      ixA^L=ixO^L;
      ixAmin1=ixOmax1+1;ixAmax1=ixOmax1+nghostcells;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      !w(ixO^S,rho_)=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,rho_)
      !w(ixO^S,p_)=pth(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2)
      do ix1=ixOmin1,ixOmax1
        w(ix1^%1ixO^S,mom(1))=w(ixOmax1+1^%1ixO^S,mom(1))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(2))=w(ixOmax1+1^%1ixO^S,mom(2))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(3))=w(ixOmax1+1^%1ixO^S,mom(3))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,rho_)=w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,p_)=pth(ixOmax1+1^%1ixO^S)
      end do
      if(stagger_grid) then
        do idir=1,nws
          if(idir==1) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix1=ixOsmax1,ixOsmin1,-1
             block%ws(ix1^%1ixOs^S,idir) = third*&
                    (-block%ws(ix1+2^%1ixOs^S,idir)&
                +4.d0*block%ws(ix1+1^%1ixOs^S,idir))
             !block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
             !       ( block%ws(ix1+3^%1ixOs^S,idir)&
             !   -5.d0*block%ws(ix1+2^%1ixOs^S,idir)&
             !   +7.d0*block%ws(ix1+1^%1ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L-kr(1,^D);
        block%ws(ixOs^S,1)=zero
        do ix1=ixOsmax1,ixOsmin1,-1
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix1^%1ixOs^S,1)=Qp(ix1+1^%1ixO^S)*block%dvolume(ix1+1^%1ixO^S)&
            /block%surfaceC(ix1^%1ixOs^S,1)
        end do
        do ix1=ixOmax1,ixOmin1,-1
          w(ix1^%1ixO^S,mag(3))=third* &
                     (-w(ix1+2,ixOmin2:ixOmax2,mag(3)) &
                +4.0d0*w(ix1+1,ixOmin2:ixOmax2,mag(3)))
        end do
        call mhd_face_to_center(ixO^L,block)
      else
        do ix1=ixOmax1,ixOmin1,-1
          w(ix1^%1ixO^S,mag(:))=third* &
                     (-w(ix1+2,ixOmin2:ixOmax2,mag(:)) &
                +4.0d0*w(ix1+1,ixOmin2:ixOmax2,mag(:)))
        end do
      end if
      if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(2)
      ixA^L=ixO^L;
      ixAmin1=ixOmin1-nghostcells;ixAmax1=ixOmin1-1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      !w(ixO^S,rho_)=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,rho_)
      !w(ixO^S,p_)=pth(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2)
      do ix1=ixOmin1,ixOmax1
        w(ix1^%1ixO^S,mom(1))=w(ixOmin1-1^%1ixO^S,mom(1))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(2))=w(ixOmin1-1^%1ixO^S,mom(2))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(3))=w(ixOmin1-1^%1ixO^S,mom(3))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,rho_)=w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,p_)=pth(ixOmin1-1^%1ixO^S)
      enddo
      if(stagger_grid) then
        do idir=1,nws
          if(idir==1) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix1=ixOsmin1,ixOsmax1
             block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
                    (-block%ws(ix1-2^%1ixOs^S,idir)&
                +4.d0*block%ws(ix1-1^%1ixOs^S,idir))
             !block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
             !       ( block%ws(ix1-3^%1ixOs^S,idir)&
             !   -5.d0*block%ws(ix1-2^%1ixOs^S,idir)&
             !   +7.d0*block%ws(ix1-1^%1ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L;
        block%ws(ixOs^S,1)=zero
        do ix1=ixOsmin1,ixOsmax1
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix1^%1ixOs^S,1)=-Qp(ix1^%1ixO^S)*block%dvolume(ix1^%1ixO^S)&
            /block%surfaceC(ix1^%1ixOs^S,1)
        end do
        call mhd_face_to_center(ixO^L,block)
        do ix1=ixOmin1,ixOmax1
          w(ix1^%1ixO^S,mag(3))=third* &
                     (-w(ix1-2^%1ixO^S,mag(3)) &
                +4.0d0*w(ix1-1^%1ixO^S,mag(3)))
        end do
      else
        do ix1=ixOmin1,ixOmax1
          w(ix1^%1ixO^S,mag(:))=third* &
                     (-w(ix1-2^%1ixO^S,mag(:)) &
                +4.0d0*w(ix1-1^%1ixO^S,mag(:)))
        end do
      end if
      if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(3)
      ixA^L=ixO^L;
      ixAmin2=ixOmax2+1;ixAmax2=ixOmax2+1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      do ix2=ixOmin2,ixOmax2
        w(ix2^%2ixO^S,rho_)=w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(1))=w(ixOmax2+1^%2ixO^S,mom(1))/w(ixOmax2+1^%2ixO^S,rho_)
        !w(ix2^%2ixO^S,mom(2))=w(ixOmax2+1^%2ixO^S,mom(2))/w(ixOmax2+1^%2ixO^S,rho_)
        !w(ix2^%2ixO^S,mom(3))=w(ixOmax2+1^%2ixO^S,mom(3))/w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,p_)=pth(ixOmax2+1^%2ixO^S)
        if(mhd_glm) w(ix2^%2ixO^S,psi_)=w(ixOmax2+1^%2ixO^S,psi_)
      enddo
      w(ixO^S,mom(2:3))=zero
      if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          !block%ws(ix2^%2ixOs^S,idir)=0.d0
          do ix2=ixOsmax2,ixOsmin2,-1
             block%ws(ix2^%2ixOs^S,idir) = third*&
                    (-block%ws(ix2+2^%2ixOs^S,idir)&
                +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L-kr(2,^D);
        block%ws(ixOs^S,2)=zero
        do ix2=ixOsmax2,ixOsmin2,-1
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix2^%2ixOs^S,2)=Qp(ix2+1^%2ixO^S)*block%dvolume(ix2+1^%2ixO^S)&
            /block%surfaceC(ix2^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
        do ix2=ixOmin2,ixOmax2
          w(ix2^%2ixO^S,mag(3))=w(ixOmax2+1^%2ixO^S,mag(3))
        end do
      else
        do ix2=ixOmin2,ixOmax2
          w(ix2^%2ixO^S,mag(2))=w(ixOmax2+1^%2ixO^S,mag(2))
          w(ix2^%2ixO^S,mag(3))=w(ixOmax2+1^%2ixO^S,mag(3))
        end do
        w(ixO^S,mag(1))=zero
      end if
      if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ixA^L=ixO^L;
      ixAmin2=ixOmin2-1;ixAmax2=ixOmin2-1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      do ix2=ixOmin2,ixOmax2
        w(ix2^%2ixO^S,rho_)=w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(1))=w(ixOmin2-1^%2ixO^S,mom(1))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(2))=w(ixOmin2-1^%2ixO^S,mom(2))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(3))=w(ixOmin2-1^%2ixO^S,mom(3))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)
      enddo
      if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix2=ixOsmin2,ixOsmax2
             block%ws(ix2^%2ixOs^S,idir) = third*&
                    (-block%ws(ix2-2^%2ixOs^S,idir)&
                +4.d0*block%ws(ix2-1^%2ixOs^S,idir))
             !block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
             !       ( block%ws(ix2-3^%2ixOs^S,idir)&
             !   -5.d0*block%ws(ix2-2^%2ixOs^S,idir)&
             !   +7.d0*block%ws(ix2-1^%2ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L;
        block%ws(ixOs^S,2)=zero
        do ix2=ixOsmin2,ixOsmax2
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix2^%2ixOs^S,2)=-Qp(ix2^%2ixO^S)*block%dvolume(ix2^%2ixO^S)&
            /block%surfaceC(ix2^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
        do ix2=ixOmin2,ixOmax2
          w(ix2^%2ixO^S,mag(3))=third* &
                      (-w(ix2-2^%2ixO^S,mag(3))&
                 +4.0d0*w(ix2-1^%2ixO^S,mag(3)))
        end do
      else
        do ix2=ixOmin2,ixOmax2
          w(ix2^%2ixO^S,mag(:))=third* &
                      (-w(ix2-2^%2ixO^S,mag(:))&
                 +4.0d0*w(ix2-1^%2ixO^S,mag(:)))
        end do
      end if
      if(mhd_solve_eaux) w(ixO^S,paux_)=w(ixO^S,p_)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  subroutine p_for_errest(ixI^L,ixO^L,iflag,w,x,var)
    integer, intent(in)           :: ixI^L,ixO^L,iflag
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: var(ixI^S)

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,var)
    
  end subroutine p_for_errest

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision :: pth(ixI^S),B2(ixI^S),divb(ixI^S)
    double precision :: Btotal(ixI^S,1:ndir),current_o(ixI^S,3)
    integer :: idir,idirmin

    ! output temperature
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)
    if(B0field) then
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,0)
    else
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
    endif
    ! B^2
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+2)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
    ! output divB
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=divb(ixO^S)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)
    ! output current
    call get_current(w,ixI^L,ixO^L,idirmin,current_o)
    w(ixO^S,nw+5)=current_o(ixO^S,1)
    w(ixO^S,nw+6)=current_o(ixO^S,2)
    w(ixO^S,nw+7)=current_o(ixO^S,3)
    ! output special resistivity eta
    call special_eta(w,ixI^L,ixO^L,idirmin,x,current_o,divb)
    w(ixO^S,nw+8)=divb(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te Alfv divB beta j1 j2 j3 eta'
  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a time-independent background magnetic field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixO^S,1)=zero
    wB0(ixO^S,2)=-Busr*dtanh(parb*x(ixO^S,1))
    wB0(ixO^S,3)=Busr/dcosh(parb*x(ixO^S,1))

  end subroutine specialset_B0

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    if(idir==1) then
      A(ixC^S) =-Busr/dcosh(parb*xC(ixC^S,1))*xC(ixC^S,2)
    else if(idir==3) then
      A(ixC^S) = Busr/parb*log(dcosh(parb*xC(ixC^S,1)))
    else
      A(ixC^S)=0.d0
    end if

  end subroutine initvecpot_usr

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here add a time-independent background current density 
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)

    wJ0(ixO^S,1)=zero
    wJ0(ixO^S,2)=parb*Busr*dtanh(parb*x(ixO^S,1))/dcosh(parb*x(ixO^S,1))
    wJ0(ixO^S,3)=-parb*Busr/dcosh(parb*x(ixO^S,1))**2

  end subroutine specialset_J0

  subroutine special_eta(w,ixI^L,ixO^L,idirmin,x,current,eta)
    ! Set the common "eta" array for resistive MHD based on w or the
    ! "current" variable which has components between idirmin and 3.
    integer, intent(in) :: ixI^L, ixO^L, idirmin
    double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision :: current(ixI^S,7-2*ndir:3), eta(ixI^S)
    double precision :: rad(ixI^S),heta,reta,eta1,eta2,etam,vc,tar
    double precision :: jc,jabs(ixI^S)
 
    heta = 6.d0
    reta = 0.8d0 * 0.3d0
    eta1 = 0.002d0
    tar= 0.4d0
    !tar= 0.0d0
    vc=1.d-4
    eta2=4.d-3
    etam=4.d-1
    if (global_time<tar) then
      rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2)-heta)**2)
      where (rad(ixO^S) .lt. reta)
        eta(ixO^S)=eta1*(2.d0*(rad(ixO^S)/reta)**3-3.d0*(rad(ixO^S)/reta)**2+1.d0)
      elsewhere
        eta(ixO^S)=zero
      endwhere
    else
      rad(ixO^S)=dsqrt(sum(current(ixO^S,:)**2,dim=ndim+1))/w(ixO^S,rho_)/q_e
      where(rad(ixO^S)>vc)
        eta(ixO^S)=eta2*(rad(ixO^S)/vc-1.d0)
      elsewhere
        eta(ixO^S)=0.d0
      endwhere
      where(eta(ixO^S)>etam)
        eta(ixO^S)=etam
      endwhere
    end if 

  end subroutine special_eta

  subroutine usrspecial_convert(qunitconvert)
    integer, intent(in) :: qunitconvert
    character(len=20):: userconvert_type
  
    call spatial_integral_w
  end subroutine usrspecial_convert

  subroutine spatial_integral_w
    double precision :: dvolume(ixG^T), dsurface(ixG^T),timephy,dvone
    double precision, allocatable :: integral_ipe(:), integral_w(:)

    integer           :: nregions,ireg,ncellpe,ncell,idims,hxM^LL,nx^D
    integer           :: iigrid,igrid,status(MPI_STATUS_SIZE),ni
    character(len=100):: filename,region
    character(len=1024) :: line, datastr
    logical           :: patchwi(ixG^T),alive

    nregions=1
    ! number of integrals to perform
    ni=3
    allocate(integral_ipe(ni),integral_w(ni))
    integral_ipe=0.d0
    integral_w=0.d0
    nx^D=ixMhi^D-ixMlo^D+1;
    do ireg=1,nregions
      select case(ireg)
      case(1)
        region='fulldomain'
      case(2)
        region='cropped'
      end select
      ncellpe=0 
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)
        if(slab) then
          dvone={rnode(rpdx^D_,igrid)|*}
          dvolume(ixM^T)=dvone
          dsurface(ixM^T)=two*(^D&dvone/rnode(rpdx^D_,igrid)+)
        else
          dvolume(ixM^T)=ps(igrid)%dvolume(ixM^T)
          dsurface(ixM^T)= sum(ps(igrid)%surfaceC(ixM^T,:),dim=ndim+1)
          do idims=1,ndim
            hxM^LL=ixM^LL-kr(idims,^D);
            dsurface(ixM^T)=dsurface(ixM^T)+ps(igrid)%surfaceC(hxM^T,idims)
          end do
        end if
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        patchwi(ixG^T)=.false.
        select case(region)
        case('cropped')
           call mask_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,patchwi,ncellpe)
        case('fulldomain')
           patchwi(ixM^T)=.true.
           ncellpe=ncellpe+{nx^D*}
        case default
           call mpistop("region not defined")
        end select
        integral_ipe(1)=integral_ipe(1)+ &
                  integral_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,dvolume,dsurface,1,patchwi)
        integral_ipe(2)=integral_ipe(2)+ &
                  integral_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,dvolume,dsurface,2,patchwi)
        integral_ipe(3)=integral_ipe(3)+ &
                  integral_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,dvolume,dsurface,3,patchwi)
      end do
      call MPI_ALLREDUCE(integral_ipe,integral_w,ni,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,icomm,ierrmpi)
      !call MPI_ALLREDUCE(ncellpe,ncell,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
      timephy=global_time
      if(mype==0) then
        write(filename,"(a,a,a)") TRIM(base_filename),TRIM(region),"mkc.csv"
        inquire(file=filename,exist=alive)
        if(alive) then
          open(unit=21,file=filename,form='formatted',status='old',access='append')
        else
          open(unit=21,file=filename,form='formatted',status='new')
          write(21,'(a)') 'time, emagnetic, einternal, current'
        endif
        write(datastr,'(es13.6, a)') timephy,','
        line=datastr
        write(datastr,"(es13.6, a)") integral_w(1),','
        line = trim(line)//trim(datastr)
        write(datastr,"(es13.6, a)") integral_w(2),','
        line = trim(line)//trim(datastr)
        write(datastr,"(es13.6)") integral_w(3)
        line = trim(line)//trim(datastr)
        write(21,'(a)') trim(line)
        close(21)
      endif
    enddo
    deallocate(integral_ipe,integral_w)
  end subroutine spatial_integral_w

  subroutine mask_grid(ixI^L,ixO^L,w,x,patchwi,cellcount)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    logical, intent(inout)             :: patchwi(ixG^T)

    double precision  ::  buff
    integer                            :: ix^D,cellcount

    buff=0.05d0*(xprobmax1-xprobmin1)
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       if(x(ix^D,1)>xprobmin1+buff .and. x(ix^D,1)<xprobmax1-buff .and. &
          x(ix^D,2)>xprobmin2+buff .and. x(ix^D,2)<xprobmax2-buff) then
         patchwi(ix^D)=.true.
         cellcount=cellcount+1
       else
         patchwi(ix^D)=.false.
       endif
    {end do\}
    return

  end subroutine mask_grid

  function integral_grid(ixI^L,ixO^L,w,x,dvolume,dsurface,intval,patchwi)
    integer, intent(in)                :: ixI^L,ixO^L,intval
    double precision, intent(in)       :: x(ixI^S,1:ndim),dvolume(ixG^T),dsurface(ixG^T)
    double precision, intent(in)       :: w(ixI^S,nw)
    logical, intent(in) :: patchwi(ixG^T)
    
    double precision, dimension(ixG^T,1:ndir) :: bvec,qvec
    double precision :: current(ixG^T,7-2*ndir:3),tmp(ixG^T)
    double precision :: integral_grid,mcurrent
    integer :: ix^D,idirmin,idir,jdir,kdir

    integral_grid=0.d0
    select case(intval)
     case(1)
      ! magnetic energy
      if(B0field)then
        tmp(ixO^S)=0.5d0*sum((w(ixO^S,mag(:))+&
                      block%B0(ixO^S,:,0))**2,dim=ndim+1)
      else
        tmp(ixO^S)=0.5d0*sum(w(ixO^S,mag(:))**2,dim=ndim+1)
      endif
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) integral_grid=integral_grid+tmp(ix^D)*dvolume(ix^D)
      {end do\}
     case(2)
      ! internal energy
      call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D))  integral_grid=integral_grid+tmp(ix^D)/(mhd_gamma-1.d0)*dvolume(ix^D)
      {end do\}
     case(3)
      ! current strength
      call get_current(w,ixI^L,ixO^L,idirmin,current)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) integral_grid=integral_grid+dsqrt(sum(current(ix^D,:)**2))*dvolume(ix^D)
      {end do\}
     case default
         call mpistop("intval not defined")
    end select
    
    return
  end function integral_grid

end module mod_usr
