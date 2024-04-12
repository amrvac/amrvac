!> Magnetic flux rope in a bipolar field
module mod_usr
  use mod_mf
  implicit none
  double precision :: q_para,d_para,L_para, delta_max, omg, tau
  double precision :: a0, I_cur, F_flx
  ! x,y,z coordinates of the flux rope axis (integral path)
  double precision, allocatable :: x_axis(:,:), f(:), s(:)
  integer :: np

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
    !usr_refine_grid     => special_refine_grid
    usr_init_vector_potential=>initvecpot_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_write_analysis => record_force_free_metrics
    usr_refine_threshold=> specialthreshold

    call set_coordinate_system("Cartesian_3D")
    call mf_activate()

  end subroutine usr_init

  !==============================================================================
  ! Purpose: to initialize user public parameters and reset global parameters.
  !          Input data are also read here.
  !==============================================================================
  subroutine initglobaldata_usr()
    integer :: ixp

    !unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    !unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    !unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    !unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    !unit_time          = unit_length/unit_velocity                     ! 85.8746159942810 s 

    !> Parameters of the potential magnetic field 
 
    !> background potential field
    q_para=5.0d20/(unit_magneticfield*unit_length**2) ! strength and sign of magnetic charges(si: /(1.0d20))
    d_para=2.0d9/unit_length ! depth of magnetic charges
    L_para=2.0d9/unit_length ! half distance between magnetic charges
  
    !> driven velocity  
    delta_max = 4.0d10/unit_length
    omg = L_para
    tau = 20.0d0
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: A(ixI^S,1:ndim)
    double precision :: Bfr(ixI^S,1:ndir)
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
    endif
      first=.false.
    endif

    if(stagger_grid) then
      call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
      call mf_face_to_center(ixO^L,block)
    else
      call bipolar_field(ixI^L,ixO^L,x,A,Bfr)
      w(ixO^S,mag(:))=Bfr(ixO^S,:)
    end if
    w(ixI^S,mom(:))=0.d0

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

    double precision :: Aphi(ixO^S),tmp(ixO^S)

    Aphi(ixO^S)= q_para*(L_para-x(ixO^S,1))/(sqrt(x(ixO^S,2)**2+(x(ixO^S,3)+d_para)**2)*&
                     sqrt(x(ixO^S,2)**2+(x(ixO^S,3)+d_para)**2+(x(ixO^S,1)-L_para)**2))+&
                 q_para*(L_para+x(ixO^S,1))/(sqrt(x(ixO^S,2)**2+(x(ixO^S,3)+d_para)**2)*&
                     sqrt(x(ixO^S,2)**2+(x(ixO^S,3)+d_para)**2+(x(ixO^S,1)+L_para)**2))

    A(ixO^S,1)=0.d0
    A(ixO^S,2)=-Aphi(ixO^S)*(x(ixO^S,3)+d_para)/sqrt(x(ixO^S,2)**2+(x(ixO^S,3)+d_para)**2)
    A(ixO^S,3)=Aphi(ixO^S)*x(ixO^S,2)/sqrt(x(ixO^S,2)**2+(x(ixO^S,3)+d_para)**2)

    if(present(Bbp)) then
      tmp(ixO^S)=-sqrt(x(ixO^S,2)**2+(x(ixO^S,3)+d_para)**2+(x(ixO^S,1)+L_para)**2)**3
      Bbp(ixO^S,1)=(x(ixO^S,1)+L_para)/tmp(ixO^S)
      Bbp(ixO^S,2)=x(ixO^S,2)/tmp(ixO^S)
      Bbp(ixO^S,3)=(x(ixO^S,3)+d_para)/tmp(ixO^S)
      tmp(ixO^S)=sqrt(x(ixO^S,2)**2+(x(ixO^S,3)+d_para)**2+(x(ixO^S,1)-L_para)**2)**3
      Bbp(ixO^S,1)=Bbp(ixO^S,1)+(x(ixO^S,1)-L_para)/tmp(ixO^S)
      Bbp(ixO^S,2)=Bbp(ixO^S,2)+x(ixO^S,2)/tmp(ixO^S)
      Bbp(ixO^S,3)=Bbp(ixO^S,3)+(x(ixO^S,3)+d_para)/tmp(ixO^S)
      Bbp(ixO^S,:)=q_para*Bbp(ixO^S,:)
    end if

  end subroutine bipolar_field

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: tmp1(ixI^S),tmp2(ixI^S),pth(ixI^S),Qp(ixI^S)
    double precision :: xlen^D,dxb^D,startpos^D,coeffrho
    integer :: idir,ix^D,ixIM^L,ixOs^L

    if(mf_glm) w(ixO^S,psi_)=0.d0

    select case(iB)
     case(1)
       !w(ixO^S,mom(1))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))
       !w(ixO^S,mom(2))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
       !w(ixO^S,mom(3))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
       w(ixO^S,mom(:))=0.d0
       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix1=ixOsmax1,ixOsmin1,-1
             block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
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

         call mf_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix1=ixOmax1,ixOmin1,-1
           w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                      (-w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                 +4.0d0*w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       end if
     case(2)
       !w(ixO^S,mom(1))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))
       !w(ixO^S,mom(2))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))
       !w(ixO^S,mom(3))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))
       w(ixO^S,mom(:))=0.d0
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
         call mf_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix1=ixOmin1,ixOmax1
           w(ix1^%1ixO^S,mag(:))=third* &
                (-w(ix1-2^%1ixO^S,mag(:)) &
           +4.0d0*w(ix1-1^%1ixO^S,mag(:)))
         enddo
       end if
     case(3)
       !w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(1))
       !w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(2))
       !w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,mom(3))
       w(ixO^S,mom(:))=0.d0
       if(stagger_grid) then
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix2=ixOsmax2,ixOsmin2,-1
             block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
                    (-block%ws(ix2+2^%2ixOs^S,idir)&
                +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
             !block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
             !       ( block%ws(ix2+3^%2ixOs^S,idir)&
             !   -5.d0*block%ws(ix2+2^%2ixOs^S,idir)&
             !   +7.d0*block%ws(ix2+1^%2ixOs^S,idir))
           end do
         end do
         ixOs^L=ixO^L-kr(2,^D);
         block%ws(ixOs^S,2)=zero
         do ix2=ixOsmax2,ixOsmin2,-1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix2^%2ixOs^S,2)=Qp(ix2+1^%2ixO^S)*block%dvolume(ix2+1^%2ixO^S)&
             /block%surfaceC(ix2^%2ixOs^S,2)
         end do
         call mf_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix2=ixOmax2,ixOmin2,-1
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                      (-w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                 +4.0d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       end if
     case(4)
       !w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(1))
       !w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(2))
       !w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,mom(3))
       w(ixO^S,mom(:))=0.d0
       if(stagger_grid) then
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix2=ixOsmin2,ixOsmax2
             block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
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
         call mf_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix2=ixOmin2,ixOmax2
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                       (-w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,mag(1):mag(3))&
                  +4.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       end if
     case(5)
       w(ixO^S,mom(:)) = 0.0d0
       if(qt .le. tau) then 
         do ix3 = ixOmin3, ixOmax3
           do ix2 = ixOmin2, ixOmax2
             do ix1 = ixOmin1, ixOmax1
               if(abs(x(ix^D,1)) .lt. omg) then
                 w(ix^D,mom(2)) = (delta_max/tau) * sin(dpi*x(ix^D,1)/omg) * &
                                  (1.0d0 - cos(2.0d0*dpi*qt/tau))
               end if 
             end do
           end do
         end do
       else
         w(ixO^S,mom(2)) = 0.0d0
       end if

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
         block%ws(ixOs^S,3)=zero
         do ix3=ixOsmax3,ixOsmin3,-1
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix3^%3ixOs^S,3)=Qp(ix3+1^%3ixO^S)*block%dvolume(ix3+1^%3ixO^S)&
             /block%surfaceC(ix3^%3ixOs^S,3)
         end do
         call mf_face_to_center(ixO^L,block)
       else
         ! 4th order accuacy constant value extrapolation
         do ix3=ixOmax3,ixOmin3,-1
           !w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/25.0d0)* &
           !     ( -3.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+4,mag(1):mag(3)) &
           !      +16.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mag(1):mag(3)) &
           !      -36.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(1):mag(3)) &
           !      +48.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(1):mag(3)))
           w(ix3^%3ixO^S,mag(:))=1.d0/11.d0* &
                ( -2.d0*w(ix3+4^%3ixO^S,mag(:)) &
                 +11.d0*w(ix3+3^%3ixO^S,mag(:)) &
                 -27.d0*w(ix3+2^%3ixO^S,mag(:)) &
                 +29.d0*w(ix3+1^%3ixO^S,mag(:)))
         enddo
       end if
      
     case(6)
       !w(ixO^S,mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(1))
       !w(ixO^S,mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(2))
       !w(ixO^S,mom(3))=abs(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(3)))
       w(ixO^S,mom(:)) = 0.0d0
       if(stagger_grid) then
         do idir=1,nws
           if(idir==3) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix3=ixOsmin3,ixOsmax3 
             block%ws(ix3^%3ixOs^S,idir) = third*&
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
         block%ws(ixOs^S,3)=zero
         do ix3=ixOsmin3,ixOsmax3
           call get_divb(w,ixI^L,ixO^L,Qp)
           block%ws(ix3^%3ixOs^S,3)=-Qp(ix3^%3ixO^S)*block%dvolume(ix3^%3ixO^S)&
             /block%surfaceC(ix3^%3ixOs^S,3)
         end do
         call mf_face_to_center(ixO^L,block)
       else
         ! 2nd order accuacy constant value extrapolation
         do ix3=ixOmin3,ixOmax3
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                 (     -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-2,mag(1):mag(3)) &
                  +4.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-1,mag(1):mag(3)))
         enddo
       end if
     case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

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
    if (any(x(ixO^S,3)<=xprobmin3+0.3d0)) then
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
    ! output divB1
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+1)=divb(ixO^S)
    ! store current
    call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
    do idir=1,ndir
      w(ixO^S,nw+1+idir)=curlvec(ixO^S,idir)
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
      w(ixO^S,nw+1+ndir+idir)=qvec(ixO^S,idir)
    end do
    w(ixO^S,nw+8)=mype

  end subroutine specialvar_output

  !==============================================================================
  ! Purpose: names for special variable output
  !==============================================================================
  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='divB j1 j2 j3 L1 L2 L3 mype'
  end subroutine specialvarnames_output

  subroutine specialthreshold(wlocal,xlocal,tolerance,qt,level)
    !PURPOSE: use different tolerance in special regions for AMR to
    !reduce/increase resolution there where nothing/something interesting happens.
    use mod_global_parameters
    
    integer, intent(in)          :: level
    double precision, intent(in) :: wlocal(1:nw),xlocal(1:ndim),qt
    double precision, intent(inout) :: tolerance
    
    double precision :: bczone^D,addtol,tol_add

    tol_add=0.d0
    !amplitude of additional tolerance
    addtol=0.3d0
    ! thickness of near-boundary region
    bczone1=4.d0
    ! linear changing of additional tolerance
    if(xlocal(1)-xprobmin1 < bczone1 .or. xprobmax1-xlocal(1) < bczone1) then
      tol_add=(1.d0-min(xlocal(1)-xprobmin1,xprobmax1-xlocal(1))/bczone1)*addtol
    endif
    bczone2=4.d0
    if(xlocal(2)-xprobmin2 < bczone2 .or. xprobmax2-xlocal(2) < bczone2) then
      tol_add=(1.d0-min(xlocal(2)-xprobmin2,xprobmax2-xlocal(2))/bczone2)*addtol
    endif
    bczone3=4.d0
    if(xprobmax3-xlocal(3) < bczone3) then
      tol_add=(1.d0-(xprobmax3-xlocal(3))/bczone3)*addtol
    endif
    tolerance=tolerance+tol_add

  end subroutine specialthreshold

end module mod_usr
