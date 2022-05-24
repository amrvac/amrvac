module mod_usr
  use mod_mhd
  implicit none

  double precision :: xJ1root, J0valatxJ1root

contains

  subroutine usr_init()
    
    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_init_vector_potential =>initvecpot_usr
    usr_refine_grid   => specialrefine_grid

    call set_coordinate_system('Cartesian')
    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr()
    double precision:: BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1

    ! this is the first root of Bessel J_1
    xJ1root=3.831705970207512d0
    ! evaluate the zero bessel function at the first root of J_1
    call JY01A(xJ1root,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
    J0valatxJ1root=BJ0

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision:: rho0,p0,epsilon,rlocal,cval
    double precision:: rval(ixG^S),xtheta(ixG^S),J1vals(ixG^S)
    double precision:: psi0(ixG^S),phi0(ixG^S),tmp(ixG^S)
    double precision:: BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1
    integer:: idims,ix^D
    logical, save :: first=.true.

    epsilon=1.0d-4
    rho0=one
    p0=one/mhd_gamma
    
    w(ix^S,rho_)=rho0
    
    rval(ixG^S)=dsqrt(x(ixG^S,1)**2+x(ixG^S,2)**2)
    xtheta(ixG^S)=datan2(x(ixG^S,2),x(ixG^S,1))

    phi0(ixG^S)=epsilon*dexp(-rval(ixG^S)**2)
    
    ! compute dphi0/dy
    idims=2
    select case(typegrad)
        case("central")
         call gradient(phi0,ixG^L,ix^L,idims,tmp)
        case("limited")
         call gradientS(phi0,ixG^L,ix^L,idims,tmp)
    end select
    w(ix^S,mom(1))=tmp(ix^S)
    ! compute dphi0/dx
    idims=1
    select case(typegrad)
         case("central")
          call gradient(phi0,ixG^L,ix^L,idims,tmp)
         case("limited")
          call gradientS(phi0,ixG^L,ix^L,idims,tmp)
    end select
    w(ix^S,mom(2))=-tmp(ix^S)
    
    ! now set pressure and magnetic field

    ! compute the local J1 evaluation
    {do ix^DB = ixG^LLIM^DB\}
      rlocal=rval(ix^D)*xJ1root
      call JY01A(rlocal,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
      J1vals(ix^D)=BJ1
    {enddo^D&\}
    
    cval=2.0d0/(xJ1root*J0valatxJ1root)
    where(rval(ixG^S)<one) 
     psi0(ixG^S)=cval*J1vals(ixG^S)*dcos(xtheta(ixG^S))
    elsewhere
     psi0(ixG^S)=(rval(ixG^S)-one/rval(ixG^S))*dcos(xtheta(ixG^S))
    endwhere
    
    where(rval(ix^S)<one) 
     w(ix^S,p_)=p0+half*(xJ1root**2)*(psi0(ix^S)**2)
    elsewhere
     w(ix^S,p_)=p0
    endwhere
      
    if(stagger_grid)then
       call b_from_vector_potential(block%ixGs^L,ixG^L,ix^L,block%ws,x)
       call mhd_face_to_center(ix^L,block)
    else
       ! compute dpsi0/dy
       idims=2
       select case(typegrad)
        case("central")
         call gradient(psi0,ixG^L,ix^L,idims,tmp)
        case("limited")
         call gradientS(psi0,ixG^L,ix^L,idims,tmp)
       end select
       w(ix^S,mag(1))=tmp(ix^S)
       ! compute dpsi0/dx
       idims=1
       select case(typegrad)
         case("central")
          call gradient(psi0,ixG^L,ix^L,idims,tmp)
         case("limited")
          call gradientS(psi0,ixG^L,ix^L,idims,tmp)
       end select
       w(ix^S,mag(2))=-tmp(ix^S)
     endif
    
    call mhd_to_conserved(ixG^L,ix^L,w,x)
    
    if(mype==0.and.first)then
          write(*,*)'Doing 2D ideal MHD, tilt problem'
          first=.false.
    endif

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)
    double precision:: rlocal,cval
    double precision:: rval(ixI^S),xtheta(ixI^S),J1vals(ixI^S)
    double precision:: BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1
    integer:: ix^D

    if (idir==3) then
       rval(ixC^S)=dsqrt(xC(ixC^S,1)**2+xC(ixC^S,2)**2)
       xtheta(ixC^S)=datan2(xC(ixC^S,2),xC(ixC^S,1))
    
       ! compute the local J1 evaluation
       {do ix^DB = ixC^LIM^DB\}
          rlocal=rval(ix^D)*xJ1root
          call JY01A(rlocal,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
          J1vals(ix^D)=BJ1
       {enddo^D&\}
    
       cval=2.0d0/(xJ1root*J0valatxJ1root)
       where(rval(ixC^S)<one) 
          A(ixC^S)=cval*J1vals(ixC^S)*dcos(xtheta(ixC^S))
       elsewhere
          A(ixC^S)=(rval(ixC^S)-one/rval(ixC^S))*dcos(xtheta(ixC^S))
       endwhere
    else
      A(ixC^S) = 0.d0
    end if

  end subroutine initvecpot_usr

  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixG^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision:: delydelx,delxdely
    integer:: ixIM^L,ix^D
    double precision:: Q(ixG^S),Qp(ixG^S)
    integer:: ixOs^L,jxO^L,idir

    select case(iB)
    ! implementation of special boundaries
    case(1)
      ixIMmin2=ixOmin2;ixIMmax2=ixOmax2;
      ixIMmin1=ixOmax1+1;ixIMmax1=ixOmax1+1;
      call mhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix1=ixOmax1,ixOmin1,-1
        w(ix1,ixOmin2:ixOmax2,rho_)= w(ixOmax1+1,ixOmin2:ixOmax2,rho_)
        w(ix1,ixOmin2:ixOmax2,p_)  = w(ixOmax1+1,ixOmin2:ixOmax2,p_)
        w(ix1,ixOmin2:ixOmax2,mom(1)) = w(ixOmax1+1,ixOmin2:ixOmax2,mom(1))
        w(ix1,ixOmin2:ixOmax2,mom(2)) = w(ixOmax1+1,ixOmin2:ixOmax2,mom(2))
        w(ix1,ixOmin2:ixOmax2,mag(1)) = 1.d0/3.d0*&
                         (-w(ix1+2,ixOmin2:ixOmax2,mag(1))+4.d0*w(ix1+1,ixOmin2:ixOmax2,mag(1)))
        w(ix1,ixOmin2:ixOmax2,mag(2)) = 1.d0/3.d0*&
                         (-w(ix1+2,ixOmin2:ixOmax2,mag(2))+4.d0*w(ix1+1,ixOmin2:ixOmax2,mag(2)))
      enddo
      !> zero normal gradient extrapolation at edge 1 (min x)
      if(stagger_grid) then
        do idir=1,nws
          if(idir==1) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix1=ixOsmax1,ixOsmin1,-1
             block%ws(ix1^%1ixOs^S,idir)=1.d0/3.d0*&
                   (-block%ws(ix1+2^%1ixOs^S,idir)&
               +4.d0*block%ws(ix1+1^%1ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L-kr(1,^D);
        block%ws(ixOs^S,1)=zero
        do ix1=ixOsmax1,ixOsmin1,-1
          call get_divb(w,ixG^L,ixO^L,Qp)
          block%ws(ix1^%1ixOs^S,1)= Qp(ix1+1^%1ixO^S)*block%dvolume(ix1+1^%1ixO^S)&
            /block%surfaceC(ix1^%1ixOs^S,1)
        end do
        call mhd_face_to_center(ixO^L,block)
      endif
      ! now reset the inner mesh values to conservative
      call mhd_to_conserved(ixG^L,ixIM^L,w,x)
      call mhd_to_conserved(ixG^L,ixO^L,w,x)
    case(2)
      ixIMmin2=ixOmin2;ixIMmax2=ixOmax2;
      ixIMmin1=ixOmin1-1;ixIMmax1=ixOmin1-1;
      call mhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix1=ixOmin1,ixOmax1,+1
        w(ix1,ixOmin2:ixOmax2,rho_)= w(ixOmin1-1,ixOmin2:ixOmax2,rho_)
        w(ix1,ixOmin2:ixOmax2,p_)  = w(ixOmin1-1,ixOmin2:ixOmax2,p_)
        w(ix1,ixOmin2:ixOmax2,mom(1)) = w(ixOmin1-1,ixOmin2:ixOmax2,mom(1))
        w(ix1,ixOmin2:ixOmax2,mom(2)) = w(ixOmin1-1,ixOmin2:ixOmax2,mom(2))
        w(ix1,ixOmin2:ixOmax2,mag(1)) = 1.d0/3.d0*&
                         (-w(ix1-2,ixOmin2:ixOmax2,mag(1))+4.d0*w(ix1-1,ixOmin2:ixOmax2,mag(1)))
        w(ix1,ixOmin2:ixOmax2,mag(2)) = 1.d0/3.d0*&
                         (-w(ix1-2,ixOmin2:ixOmax2,mag(2))+4.d0*w(ix1-1,ixOmin2:ixOmax2,mag(2)))
      enddo
      !> zero normal gradient extrapolation at edge 2 (max x)
      if(stagger_grid) then
        do idir=1,nws
          if(idir==1) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix1=ixOsmin1,ixOsmax1
             block%ws(ix1^%1ixOs^S,idir)=1.d0/3.d0*&
                   (-block%ws(ix1-2^%1ixOs^S,idir)&
               +4.d0*block%ws(ix1-1^%1ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L;
        jxO^L=ixO^L-nghostcells*kr(1,^D);
        block%ws(ixOs^S,1)=zero
        call get_divb(w,ixG^L,jxO^L,Q)
        do ix1=ixOsmin1,ixOsmax1
          call get_divb(w,ixG^L,ixO^L,Qp)
          block%ws(ix1^%1ixOs^S,1)=&
            (Q(jxOmax1^%1jxO^S)*block%dvolume(jxOmax1^%1jxO^S)&
               -Qp(ix1^%1ixO^S)*block%dvolume(ix1^%1ixO^S))&
            /block%surfaceC(ix1^%1ixOs^S,1)
        end do
        call mhd_face_to_center(ixO^L,block)
      endif
      ! now reset the inner mesh values to conservative
      call mhd_to_conserved(ixG^L,ixIM^L,w,x)
      call mhd_to_conserved(ixG^L,ixO^L,w,x)
    case(3)
      ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+1;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      call mhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix2=ixOmax2,ixOmin2,-1
        w(ixOmin1:ixOmax1,ix2,rho_)= w(ixOmin1:ixOmax1,ixOmax2+1,rho_)
        w(ixOmin1:ixOmax1,ix2,p_)  = w(ixOmin1:ixOmax1,ixOmax2+1,p_)
        w(ixOmin1:ixOmax1,ix2,mom(1)) = w(ixOmin1:ixOmax1,ixOmax2+1,mom(1))
        w(ixOmin1:ixOmax1,ix2,mom(2)) = w(ixOmin1:ixOmax1,ixOmax2+1,mom(2))
        w(ixOmin1:ixOmax1,ix2,mag(1)) = 1.0d0/3.d0*&
                   (-w(ixOmin1:ixOmax1,ix2+2,mag(1))+4.0d0*w(ixOmin1:ixOmax1,ix2+1,mag(1)))
        w(ixOmin1:ixOmax1,ix2,mag(2)) = 1.0d0/3.d0*&
                   (-w(ixOmin1:ixOmax1,ix2+2,mag(2))+4.0d0*w(ixOmin1:ixOmax1,ix2+1,mag(2)))
      enddo
      !> zero normal gradient extrapolation at edge 3 (min y)
      if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix2=ixOsmax2,ixOsmin2,-1
             block%ws(ix2^%2ixOs^S,idir)=1.d0/3.d0*&
                   (-block%ws(ix2+2^%2ixOs^S,idir)&
               +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L-kr(2,^D);
        block%ws(ixOs^S,2)=zero
        do ix2=ixOsmax2,ixOsmin2,-1
          call get_divb(w,ixG^L,ixO^L,Qp)
          block%ws(ix2^%2ixOs^S,2)= Qp(ix2+1^%2ixO^S)*block%dvolume(ix2+1^%2ixO^S)&
            /block%surfaceC(ix2^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
      endif
      ! now reset the inner mesh values to conservative
      call mhd_to_conserved(ixG^L,ixIM^L,w,x)
      call mhd_to_conserved(ixG^L,ixO^L,w,x)
    case(4)
      ixIMmin2=ixOmin2-1;ixIMmax2=ixOmin2-1;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      call mhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix2=ixOmin2,ixOmax2,+1
        w(ixOmin1:ixOmax1,ix2,rho_)= w(ixOmin1:ixOmax1,ixOmin2-1,rho_)
        w(ixOmin1:ixOmax1,ix2,p_)  = w(ixOmin1:ixOmax1,ixOmin2-1,p_)
        w(ixOmin1:ixOmax1,ix2,mom(1)) = w(ixOmin1:ixOmax1,ixOmin2-1,mom(1))
        w(ixOmin1:ixOmax1,ix2,mom(2)) = w(ixOmin1:ixOmax1,ixOmin2-1,mom(2))
        w(ixOmin1:ixOmax1,ix2,mag(1)) = 1.0d0/3.d0*&
                   (-w(ixOmin1:ixOmax1,ix2-2,mag(1))+4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(1)))
        w(ixOmin1:ixOmax1,ix2,mag(2)) = 1.0d0/3.d0*&
                   (-w(ixOmin1:ixOmax1,ix2-2,mag(2))+4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(2)))
      enddo
      !> zero normal gradient extrapolation at edge 4 (max y)
      if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix2=ixOsmin2,ixOsmax2
             block%ws(ix2^%2ixOs^S,idir)=1.d0/3.d0*&
                   (-block%ws(ix2-2^%2ixOs^S,idir)&
               +4.d0*block%ws(ix2-1^%2ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L;
        jxO^L=ixO^L-nghostcells*kr(2,^D);
        block%ws(ixOs^S,2)=zero
        call get_divb(w,ixG^L,jxO^L,Q)
        do ix2=ixOsmin2,ixOsmax2
          call get_divb(w,ixG^L,ixO^L,Qp)
          block%ws(ix2^%2ixOs^S,2)=&
            (Q(jxOmax2^%2jxO^S)*block%dvolume(jxOmax2^%2jxO^S)&
               -Qp(ix2^%2ixO^S)*block%dvolume(ix2^%2ixO^S))&
            /block%surfaceC(ix2^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
      endif
      ! now reset the inner mesh values to conservative
      call mhd_to_conserved(ixG^L,ixIM^L,w,x)
      call mhd_to_conserved(ixG^L,ixO^L,w,x)
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
    double precision:: rval(ixG^S)

    rval(ix^S)=dsqrt(x(ix^S,1)**2+x(ix^S,2)**2)
    if (all(rval(ix^S) < one)) refine=1

  end subroutine specialrefine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
    ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision:: divb(ixI^S), wlocal(ixI^S,1:nw)
    double precision :: current(ixI^S,7-2*ndir:3)
    integer          :: idirmin

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call get_current(wlocal,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+1)=current(ixO^S,3)
    call get_divb(wlocal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+2)=divb(ixO^S)
    
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the varnames/primnames string
  character(len=*) :: varnames 
  varnames='jz divb'

  end subroutine specialvarnames_output

  SUBROUTINE JY01A(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
  
  !       =======================================================
  !       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
  !                Y1(x), and their derivatives
  !       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ò 0 )
  !       Output:  BJ0 --- J0(x)
  !                DJ0 --- J0'(x)
  !                BJ1 --- J1(x)
  !                DJ1 --- J1'(x)
  !                BY0 --- Y0(x)
  !                DY0 --- Y0'(x)
  !                BY1 --- Y1(x)
  !                DY1 --- Y1'(x)
  !       =======================================================
  !
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  
  double precision:: X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1
  double precision:: PI,RP2,X2,R,EC,CS0,W0,R0,CS1,W1,R1
  double precision:: T1,P0,Q0,T2,P1,Q1,CU
  double precision:: A(12),B(12),A1(12),B1(12)
  
  integer:: K,K0
  
          PI=3.141592653589793D0
          RP2=0.63661977236758D0
          X2=X*X
          IF (X==0.0D0) THEN
             BJ0=1.0D0
             BJ1=0.0D0
             DJ0=0.0D0
             DJ1=0.5D0
             BY0=-1.0D+300
             BY1=-1.0D+300
             DY0=1.0D+300
             DY1=1.0D+300
             RETURN
          ENDIF
          IF (X<=12.0D0) THEN
             BJ0=1.0D0
             R=1.0D0
             DO K=1,30
                R=-0.25D0*R*X2/(K*K)
                BJ0=BJ0+R
                IF (DABS(R)<DABS(BJ0)*1.0D-15) EXIT
             ENDDO
             BJ1=1.0D0
             R=1.0D0
             DO K=1,30
                R=-0.25D0*R*X2/(K*(K+1.0D0))
                BJ1=BJ1+R
                IF (DABS(R)<DABS(BJ1)*1.0D-15) EXIT
             ENDDO
             BJ1=0.5D0*X*BJ1
             EC=DLOG(X/2.0D0)+0.5772156649015329D0
             CS0=0.0D0
             W0=0.0D0
             R0=1.0D0
             DO K=1,30
                W0=W0+1.0D0/K
                R0=-0.25D0*R0/(K*K)*X2
                R=R0*W0
                CS0=CS0+R
                IF (DABS(R)<DABS(CS0)*1.0D-15) EXIT
             ENDDO
             BY0=RP2*(EC*BJ0-CS0)
             CS1=1.0D0
             W1=0.0D0
             R1=1.0D0
             DO K=1,30
                W1=W1+1.0D0/K
                R1=-0.25D0*R1/(K*(K+1))*X2
                R=R1*(2.0D0*W1+1.0D0/(K+1.0D0))
                CS1=CS1+R
                IF (DABS(R)<DABS(CS1)*1.0D-15) EXIT
             ENDDO
             BY1=RP2*(EC*BJ1-1.0D0/X-0.25D0*X*CS1)
          ELSE
           A(1)=-.7031250000000000D-01
           A(2)=.1121520996093750D+00
           A(3)=-.5725014209747314D+00
           A(4)=.6074042001273483D+01
           A(5)=-.1100171402692467D+03
           A(6)=.3038090510922384D+04
           A(7)=-.1188384262567832D+06
           A(8)=.6252951493434797D+07
           A(9)=-.4259392165047669D+09
           A(10)=.3646840080706556D+11
           A(11)=-.3833534661393944D+13
           A(12)=.4854014686852901D+15
  !DATA A/-.7031250000000000D-01,.1121520996093750D+00,-.5725014209747314D+00,.6074042001273483D+01, &
  !       -.1100171402692467D+03,.3038090510922384D+04,-.1188384262567832D+06,.6252951493434797D+07, &
  !       -.4259392165047669D+09,.3646840080706556D+11,-.3833534661393944D+13,.4854014686852901D+15/
           B(1)=.7324218750000000D-01
           B(2)=-.2271080017089844D+00
           B(3)=.1727727502584457D+01
           B(4)=-.2438052969955606D+02
           B(5)=.5513358961220206D+03
           B(6)=-.1825775547429318D+05
           B(7)=.8328593040162893D+06
           B(8)=-.5006958953198893D+08
           B(9)=.3836255180230433D+10
           B(10)=-.3649010818849833D+12
           B(11)=.4218971570284096D+14
           B(12)=-.5827244631566907D+16
  !DATA B/ .7324218750000000D-01,-.2271080017089844D+00,.1727727502584457D+01,-.2438052969955606D+02, &
  !        .5513358961220206D+03,-.1825775547429318D+05,.8328593040162893D+06,-.5006958953198893D+08, &
  !        .3836255180230433D+10,-.3649010818849833D+12,.4218971570284096D+14,-.5827244631566907D+16/
           A1(1)=.1171875000000000D+00
           A1(2)=-.1441955566406250D+00
           A1(3)=.6765925884246826D+00
           A1(4)=-.6883914268109947D+01
           A1(5)=.1215978918765359D+03
           A1(6)=-.3302272294480852D+04
           A1(7)=.1276412726461746D+06
           A1(8)=-.6656367718817688D+07
           A1(9)=.4502786003050393D+09
           A1(10)=-.3833857520742790D+11
           A1(11)=.4011838599133198D+13
           A1(12)=-.5060568503314727D+15
  !DATA A1/.1171875000000000D+00,-.1441955566406250D+00,.6765925884246826D+00,-.6883914268109947D+01, &
  !        .1215978918765359D+03,-.3302272294480852D+04,.1276412726461746D+06,-.6656367718817688D+07, &
  !        .4502786003050393D+09,-.3833857520742790D+11,.4011838599133198D+13,-.5060568503314727D+15/
           B1(1)=-.1025390625000000D+00
           B1(2)=.2775764465332031D+00
           B1(3)=-.1993531733751297D+01
           B1(4)=.2724882731126854D+02
           B1(5)=-.6038440767050702D+03
           B1(6)=.1971837591223663D+05
           B1(7)=-.8902978767070678D+06
           B1(8)=.5310411010968522D+08
           B1(9)=-.4043620325107754D+10
           B1(10)=.3827011346598605D+12
           B1(11)=-.4406481417852278D+14
           B1(12)=.6065091351222699D+16
  !DATA B1/-.1025390625000000D+00,.2775764465332031D+00,-.1993531733751297D+01,.2724882731126854D+02, &
  !        -.6038440767050702D+03,.1971837591223663D+05,-.8902978767070678D+06,.5310411010968522D+08, &
  !        -.4043620325107754D+10,.3827011346598605D+12,-.4406481417852278D+14,.6065091351222699D+16/
             K0=12
             IF (X>=35.0) K0=10
             IF (X>=50.0) K0=8
             T1=X-0.25D0*PI
             P0=1.0D0
             Q0=-0.125D0/X
             DO K=1,K0
                P0=P0+A(K)*X**(-2*K)
                Q0=Q0+B(K)*X**(-2*K-1)
             ENDDO
             CU=DSQRT(RP2/X)
             BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
             BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
             T2=X-0.75D0*PI
             P1=1.0D0
             Q1=0.375D0/X
             DO K=1,K0
                P1=P1+A1(K)*X**(-2*K)
                Q1=Q1+B1(K)*X**(-2*K-1)
             ENDDO
             CU=DSQRT(RP2/X)
             BJ1=CU*(P1*DCOS(T2)-Q1*DSIN(T2))
             BY1=CU*(P1*DSIN(T2)+Q1*DCOS(T2))
          ENDIF
          DJ0=-BJ1
          DJ1=BJ0-BJ1/X
          DY0=-BY1
          DY1=BY0-BY1/X
          RETURN
          END

end module mod_usr
