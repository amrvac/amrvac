! bipole magnetic field with large Lorenz force from Mackay 2001 ApJ, 560, 445
module mod_usr
  use mod_mf
  implicit none
  double precision :: beta,sp
contains

  subroutine usr_init()
    unit_length=1.d9 !cm
    unit_numberdensity=1.d9!cm^-3
    unit_temperature=1.d6 !K

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_init_vector_potential=>initvecpot_usr
    usr_special_bc    => specialbound_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_write_analysis => record_force_free_metrics

    call set_coordinate_system("Cartesian")
    call mf_activate()

  end subroutine usr_init
    
  subroutine initglobaldata_usr
    character(len=20) :: printsettingformat  

    printsettingformat='(1x,A20,ES15.7,A30)'
    if(mype==0) then
      write(*,*) "Dimensionless units:"
      write(*,printsettingformat) "unit_length ",unit_length,"(cm)"
      write(*,printsettingformat) "unit_velocity ",unit_velocity,"(cm s^-1)"
      write(*,printsettingformat) "unit_density ",unit_density,"(g cm^-3)"
      write(*,*) "Deduced dimensionless units:"
      write(*,printsettingformat) "unit_time ",unit_time,"(s)"
      write(*,printsettingformat) "unit_pressure ",unit_pressure,"(Ba = g cm^-1 s^-2)"
      write(*,printsettingformat) "Temperature unit",unit_temperature,"(K)"
      write(*,printsettingformat) "unit_magneticfield ",unit_magneticfield,"(G=g^{1/2} cm^{-1/2} s^-1)"
    end if
    ! magnetic field strength
    Busr=Busr/unit_magneticfield

    beta=0.d0
    sp=1.d0

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bloc(ixI^S,ndir),xC(ixI^S,1:ndim)

    integer :: idir, ixC^L
    logical, save:: first=.true.
    logical :: vector_potential

    if (first) then
       if (mype==0) then
          print *,'3D coronal bipole in Cartesian coordinate'
       end if
       first=.false.
    end if
    if(stagger_grid) then
      ! calculate magnetic field from vector potential
      vector_potential=.true.
      if(vector_potential) then
        call b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block%ws,x)
      else
        do idir=1,ndim
          xC=x
          ixCmax^D=ixOmax^D;
          ixCmin^D=ixOmin^D-kr(idir,^D);
          xC(ixC^S,idir)=x(ixC^S,idir)+0.5d0*block%dx(ixC^S,idir)
          call get_B(ixI^L,ixC^L,Bloc,xC)
          block%ws(ixC^S,idir)=Bloc(ixC^S,idir)
        end do
      end if
      call mf_face_to_center(ixO^L,block)
    else 
     call get_B(ixI^L,ixO^L,Bloc,x)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
    end if
    w(ixO^S,mom(:))=0.d0

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L ,ixC^L,x,A,idir)
    integer, intent(in)                :: ixI^L,ixC^L,idir

    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)
    double precision  ::  xi(ixI^S)

    xi(ixC^S)=((x(ixC^S,1)**2+x(ixC^S,3)**2)*0.5d0+x(ixC^S,2)**2)/sp**2

   ! spherical vector potential 
    if (idir==1) then 
      A(ixC^S)=beta*Busr*exp(0.5d0)*x(ixC^S,3)*exp(-2.d0*xi(ixC^S))
    else if (idir==2)  then 
      A(ixC^S)=Busr*exp(0.5d0)*sp*exp(-xi(ixC^S))
    else
      A(ixC^S)=-beta*Busr*exp(0.5d0)*x(ixC^S,1)*exp(-2.d0*xi(ixC^S))
    end if 

  end subroutine initvecpot_usr

  subroutine get_B(ixI^L,ixO^L,B,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: B(ixI^S,1:ndir)

    double precision :: xi(ixI^S),L(ixI^S),y(ixI^S),z(ixI^S),M(ixI^S)
    double precision :: N(ixI^S),Bx(ixI^S),By(ixI^S),Bz(ixI^S)

    xi(ixO^S)=((x(ixO^S,1)**2+x(ixO^S,3)**2)*0.5d0+x(ixO^S,2)**2)/sp**2
    B(ixO^S,1)=Busr*exp(0.5d0)*(x(ixO^S,3)*exp(-xi(ixO^S))/sp+4.d0*beta*x(ixO^S,1)*x(ixO^S,2)*&
              exp(-2.d0*xi(ixO^S))/sp**2)
    B(ixO^S,2)=2.d0*beta*Busr*exp(0.5d0)*(1.d0-(x(ixO^S,1)**2+x(ixO^S,3)**2)/sp**2)*exp(-2.d0*xi(ixO^S))
    B(ixO^S,3)=Busr*exp(0.5d0)*(-x(ixO^S,1)*exp(-xi(ixO^S))/sp+4.d0*beta*x(ixO^S,2)*x(ixO^S,3)*&
              exp(-2.d0*xi(ixO^S))/sp**2)

  end subroutine get_B

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: ft,tfstop,tramp1,tramp2,coeffrho,vlimit,vsign
    double precision :: xlen^D,Q(ixI^S),Qp(ixI^S)
    integer :: ix^D,ixbc^D,ixOs^L,jxO^L,idir

    select case(iB)
    case(1)
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
        do ix1=ixOmax1,ixOmin1,-1
          w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                     (-w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
               +4.0d0*w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
        end do
      end if
    case(2)
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
          do ix1=ixOmin1,ixOmax1
            w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                       (-w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                  +4.0d0*w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
          enddo
      end if
    case(3)
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
        do ix2=ixOmax2,ixOmin2,-1
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                     (-w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                +4.0d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,mag(1):mag(3)))
        end do
      end if
    case(4)
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
        do ix2=ixOmin2,ixOmax2
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                     (-w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                +4.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(1):mag(3)))
        enddo
      end if
    case(5)
      w(ixO^S,mom(:))=0.d0
      if(stagger_grid) then
        do idir=1,nws
          if(idir==3) cycle
            ixOsmax^D=ixOmax^D;
            ixOsmin^D=ixOmin^D-kr(^D,idir);
            do ix3=ixOsmax3,ixOsmin3,-1
              ! 4 one-sided equal grdient
              block%ws(ix3^%3ixOs^S,idir)= &
                0.12d0*block%ws(ix3+5^%3ixOs^S,idir) &
               -0.76d0*block%ws(ix3+4^%3ixOs^S,idir) &
               +2.08d0*block%ws(ix3+3^%3ixOs^S,idir) &
               -3.36d0*block%ws(ix3+2^%3ixOs^S,idir) &
               +2.92d0*block%ws(ix3+1^%3ixOs^S,idir)
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
        call mf_face_to_center(ixO^L,block)
      else
          ! 4 one-sided equal grdient
          do ix3=ixOmax3,ixOmin3,-1
            w(ix3^%3ixO^S,mag(:))= &
                0.12d0*w(ix3+5^%3ixO^S,mag(:)) &
               -0.76d0*w(ix3+4^%3ixO^S,mag(:)) &
               +2.08d0*w(ix3+3^%3ixO^S,mag(:)) &
               -3.36d0*w(ix3+2^%3ixO^S,mag(:)) &
               +2.92d0*w(ix3+1^%3ixO^S,mag(:))
          enddo
      end if
    case(6)
      w(ixO^S,mom(:))=0.d0
      if(stagger_grid) then
        do idir=1,nws
          if(idir==3) cycle
            ixOsmax^D=ixOmax^D;
            ixOsmin^D=ixOmin^D-kr(^D,idir);
            do ix3=ixOsmin3,ixOsmax3
               block%ws(ix3^%3ixOs^S,idir) = 1.d0/11.d0*&
                  (+2.d0*block%ws(ix3-3^%3ixOs^S,idir)&
                   -9.d0*block%ws(ix3-2^%3ixOs^S,idir)&
                  +18.d0*block%ws(ix3-1^%3ixOs^S,idir))
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
        call mf_face_to_center(ixO^L,block)
      else
          do ix3=ixOmin3,ixOmax3
               w(ix3^%3ixO^S,mag(:)) = 1.d0/11.d0*&
                  (+2.d0*w(ix3-3^%3ixO^S,mag(:))&
                   -9.d0*w(ix3-2^%3ixO^S,mag(:))&
                  +18.d0*w(ix3-1^%3ixO^S,mag(:)))
          enddo
      end if
    case default
      call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

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

    double precision                   :: current(ixI^S,7-2*ndir:3)
    double precision                   :: qvec(ixI^S,1:ndir)
    double precision                   :: tmp(ixI^S)
    integer :: idirmin,idir,jdir,kdir

    ! store current
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    ! calculate Lorentz force
    w(ixO^S,nw+1)=current(ixO^S,1)
    w(ixO^S,nw+2)=current(ixO^S,2)
    w(ixO^S,nw+3)=current(ixO^S,3)
    qvec(ixO^S,1:ndir)=zero
    do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
       if(lvc(idir,jdir,kdir)/=0)then
          tmp(ixO^S)=current(ixO^S,jdir)*w(ixO^S,mag(kdir))
          if(lvc(idir,jdir,kdir)==1)then
             qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
          else
             qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
          endif
       endif
    enddo; enddo; enddo
    w(ixO^S,nw+4) =qvec(ixO^S,1)
    w(ixO^S,nw+5)=qvec(ixO^S,2)
    w(ixO^S,nw+6)=qvec(ixO^S,3)
    ! output divB1
    call get_divb(w,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+7)=tmp(ixO^S)
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='j1 j2 j3 L1 L2 L3 divb'

  end subroutine specialvarnames_output

end module mod_usr
