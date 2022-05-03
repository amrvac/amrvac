!rho0 = 1.66E-13 kg/m^3
!R0 = 1 Rs = 695,700,000.0
!v0 = 100 km/s = 100,000.0 m/s
!as a result, we have:
!t0 = 6957 s
!B0 = 4.567E-5 T
!P0 = 0.00166 Pa
!T0 = 601,449.28 K    if we assume p=2nkT

module mod_usr
  use mod_mhd
  implicit none
  double precision :: Bbnd,rhobnd,Tbnd,gamma0,gamma2,gamma4,usr_grav

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    mhd_gamma=5.0d0/3.0d0
    mhd_eta=zero

    call set_coordinate_system("spherical_2.5D")

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source => specialsource
    usr_refine_grid     => specialrefine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_set_B0          => specialset_B0
    usr_init_vector_potential=>initvecpot_usr

    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters

    ! we want to have 1.1 Gauss at the Solar surface at the equator
    Bbnd = 2.40858d0
    ! we want a number density of n=0.2142E8 cm-3 at the inner boundary
    rhobnd = 0.2142857d0
    ! we want a temperature of 1.5E6 K at the inner boundary
    Tbnd = 2.493976d0
    ! solar differential rotation:
    gamma0 = 0.0194796d0/2.0d0/dpi
    gamma2 = -0.0024349d0/2.0d0/dpi
    gamma4 = -0.0031306d0/2.0d0/dpi

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    logical, save :: first=.true.

    if (first .and. mype==0) then
       write(*,*)'Running solar wind'
       first=.false.
    endif

    ! n=0.2142E8 cm-3 at the inner boundary and decreases like r^-2
    w(ixO^S,rho_)=rhobnd/x(ixO^S,1)**2
    ! vr is 200 km/s everywhere
    w(ixO^S,mom(1))  = 2.0d0*w(ixO^S,rho_)
    w(ixO^S,mom(2))  = 0.0d0
    ! at the inner boundary: differential rotation of the Sun
    ! vphi decreases like r^-2
    w(ixO^S,mom(3))  = w(ixO^S,rho_)*x(ixO^S,1)*dsin(x(ixO^S,2))* &
             (gamma0+gamma2*dcos(x(ixO^S,2))**2+gamma4*dcos(x(ixO^S,2))**4)/x(ixO^S,1)**2
    if(B0field) then
      w(ixO^S,mag(1))  = 0.0d0
      w(ixO^S,mag(2))  = 0.0d0
    else if(stagger_grid) then
      call b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block%ws,x)
      call mhd_face_to_center(ixO^L,block)
    else
      w(ixO^S,mag(1))=2.0d0*Busr*dcos(x(ixO^S,2))/x(ixO^S,1)**3
      w(ixO^S,mag(2))=Busr*dsin(x(ixO^S,2))/x(ixO^S,1)**3
    end if
    w(ixO^S,mag(3))  = 0.0d0
    !pressure decreases like 1/r^2 (T is thus fixed)
    w(ixO^S,e_)   = Tbnd*rhobnd/x(ixO^S,1)**2/(mhd_gamma-one) &
                  +0.5d0*(w(ixO^S,mom(1))**2+w(ixO^S,mom(2))**2+w(ixO^S,mom(3))**2)/w(ixO^S,rho_) &
                  +0.5d0*(w(ixO^S,mag(1))**2+w(ixO^S,mag(2))**2+w(ixO^S,mag(3))**2)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    if (idir==3) then
      A(ixC^S) = Busr/xC(ixC^S,1)**2*sin(xC(ixC^S,2))
    else
      A(ixC^S) = 0.d0
    end if

  end subroutine initvecpot_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    use mod_constrained_transport

    integer, intent(in) :: ixI^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S)
    integer :: ix1,ix2,ixA^L,ixOs^L
    double precision :: omega, LM, Pg1, Pg2, Ti, Pin, latitude

    select case (iB)
    case (1)

       !!!!!!!!!!!!!!!!!!
       ! inner boundary !
       !!!!!!!!!!!!!!!!!!
       if(stagger_grid) then
         ixOs^L=ixO^L-kr(1,^D);
         block%ws(ixOs^S,1)=2.0d0*Busr*dcos(x(ixO^S,2))/(x(ixO^S,1)-0.5d0*block%dx(ixO^S,1))**3
         ixOsmax^D=ixOmax^D;
         ixOsmin^D=ixOmin^D-kr(^D,1);
         do ix1=ixOsmin1,ixOsmax1
           block%ws(ix1^%1ixOs^S,2)=block%ws(ixOsmax1+1^%1ixOs^S,2)
         !  block%ws(ix1^%1ixOs^S,1)=block%ws(ixOsmax1+1^%1ixOs^S,1)*block%surfaceC(ixOsmax1+1^%1ixOs^S,1)/block%surfaceC(ix1^%1ixOs^S,1)
         end do
         call mhd_face_to_center(ixO^L,block)
       end if

       ! Loop through ghost cells along boundary (in the theta direction)
       do ix2=ixOmin2,ixOmax2

         ! The density is fixed at the inner boundary
         w(ixOmax1,ix2,rho_) = 2.0d0*rhobnd-w(ixOmax1+1,ix2,rho_)
         w(ixOmin1,ix2,rho_) = 4.0d0*rhobnd-3.0d0*w(ixOmax1+1,ix2,rho_)

         ! vphi/(r sin(theta)) is fixed at the inner boundary (this is the angular velocity)
         ! this is where we impose the differential rotation of the Sun
         omega=gamma0+gamma2*(dcos(x(ixOmax1,ix2,2)))**2+gamma4*(dcos(x(ixOmax1,ix2,2)))**4
         w(ixOmax1,ix2,mom(3)) = x(ixOmax1,ix2,1)*dsin(x(ixOmax1,ix2,2))*w(ixOmax1,ix2,rho_)&
                  *(2.0d0*omega-w(ixOmax1+1,ix2,mom(3))/(x(ixOmax1+1,ix2,1)*dsin(x(ixOmax1+1,ix2,2))*w(ixOmax1+1,ix2,rho_)))
         w(ixOmin1,ix2,mom(3)) = x(ixOmin1,ix2,1)*dsin(x(ixOmin1,ix2,2))*w(ixOmin1,ix2,rho_)&
                  *(4.0d0*omega-3.0d0*w(ixOmax1+1,ix2,mom(3))/(x(ixOmax1+1,ix2,1)*dsin(x(ixOmax1+1,ix2,2))*w(ixOmax1+1,ix2,rho_)))

         if(.not.stagger_grid) then
           ! Btheta decreases like a dipole
           ! this might not be a very good boundary condition
           ! Btheta(g) = Btheta(I) + (BthetaDIPOLE(g) - BthetaDIPOLE(I))     "g" is for ghostcell and "I" for first inner cell
           w(ixOmax1,ix2,mag(2)) = w(ixOmax1+1,ix2,mag(2))
           w(ixOmin1,ix2,mag(2)) = w(ixOmax1+1,ix2,mag(2))

           ! r^2*Br is fixed at the inner boundary (this is how we force a dipole field for the Sun)
           w(ixOmax1,ix2,mag(1)) = (2.0d0*2.0d0*Bbnd*dcos(x(ixOmax1,ix2,2)) &
                  -(w(ixOmax1+1,ix2,mag(1))+Bbnd*2.0d0*dcos(x(ixOmax1+1,ix2,2))/x(ixOmax1+1,ix2,1)**3)*x(ixOmax1+1,ix2,1)**2)/x(ixOmax1,ix2,1)**2 &
                  -Bbnd*2.0d0*dcos(x(ixOmax1,ix2,2))/x(ixOmax1,ix2,1)**3
           w(ixOmin1,ix2,mag(1)) = (4.0d0*2.0d0*Bbnd*dcos(x(ixOmin1,ix2,2)) &
                  -3.0d0*(w(ixOmax1+1,ix2,mag(1))+Bbnd*2.0d0*dcos(x(ixOmax1+1,ix2,2))/x(ixOmax1+1,ix2,1)**3)*x(ixOmax1+1,ix2,1)**2)/x(ixOmin1,ix2,1)**2 &
                  -Bbnd*2.0d0*dcos(x(ixOmin1,ix2,2))/x(ixOmin1,ix2,1)**3
         end if

         latitude = 0.5d0*dpi-x(ixOmax1,ix2,2)
         ! if outside of the deadzone
         ! r^2*rho*Vr is constant at the inner boundary
         if(abs(latitude) >= 0.3927d0) then
            w(ixOmax1,ix2,mom(1)) = w(ixOmax1+1,ix2,mom(1))*(x(ixOmax1+1,ix2,1)/x(ixOmax1,ix2,1))**2
            w(ixOmin1,ix2,mom(1)) = w(ixOmax1+1,ix2,mom(1))*(x(ixOmax1+1,ix2,1)/x(ixOmin1,ix2,1))**2
            ! we choose Mtheta in order to have vect(B) and vect(v) parallel
            w(ixOmax1,ix2,mom(2)) = w(ixOmax1,ix2,mom(1))*(w(ixOmax1,ix2,mag(2))+Bbnd*dsin(x(ixOmax1,ix2,2))/x(ixOmax1,ix2,1)**3) &
                                  / (w(ixOmax1,ix2,mag(1))+Bbnd*2.0d0*dcos(x(ixOmax1,ix2,2))/x(ixOmax1,ix2,1)**3)
            w(ixOmin1,ix2,mom(2)) = w(ixOmin1,ix2,mom(1))*(w(ixOmin1,ix2,mag(2))+Bbnd*dsin(x(ixOmin1,ix2,2))/x(ixOmin1,ix2,1)**3) &
                                  / (w(ixOmin1,ix2,mag(1))+Bbnd*2.0d0*dcos(x(ixOmin1,ix2,2))/x(ixOmin1,ix2,1)**3)
            !r^3*(rho*vr*vphi-Br*Bphi) is continuous at the inner boundary
            LM=x(ixOmax1+1,ix2,1)**3*(w(ixOmax1+1,ix2,mom(1))*w(ixOmax1+1,ix2,mom(3))/w(ixOmax1+1,ix2,rho_) &
                 - (w(ixOmax1+1,ix2,mag(1))+Bbnd*2.0d0*dcos(x(ixOmax1+1,ix2,2))/x(ixOmax1+1,ix2,1)**3)*w(ixOmax1+1,ix2,mag(3)))
            w(ixOmax1,ix2,mag(3)) = (w(ixOmax1,ix2,mom(1))*w(ixOmax1,ix2,mom(3))/w(ixOmax1,ix2,rho_)-LM/x(ixOmax1,ix2,1)**3) &
                                  / (w(ixOmax1,ix2,mag(1))+Bbnd*2.0d0*dcos(x(ixOmax1,ix2,2))/x(ixOmax1,ix2,1)**3)
            w(ixOmin1,ix2,mag(3)) = (w(ixOmin1,ix2,mom(1))*w(ixOmin1,ix2,mom(3))/w(ixOmin1,ix2,rho_)-LM/x(ixOmin1,ix2,1)**3) &
                                  / (w(ixOmin1,ix2,mag(1))+Bbnd*2.0d0*dcos(x(ixOmin1,ix2,2))/x(ixOmin1,ix2,1)**3)
         else
            ! if inside the deadzone: m1 is zero at the boundary
            w(ixOmax1,ix2,mom(1)) = -w(ixOmax1+1,ix2,mom(1))
            w(ixOmin1,ix2,mom(1)) = -3.0d0*w(ixOmax1+1,ix2,mom(1))
            ! if inside the deadzone: m2 is zero at the boundary
            w(ixOmax1,ix2,mom(2)) = -w(ixOmax1+1,ix2,mom(2))
            w(ixOmin1,ix2,mom(2)) = -3.0d0*w(ixOmax1+1,ix2,mom(2))
            ! if inside the deadzone: Bphi is continuous
            w(ixOmax1,ix2,mag(3)) = w(ixOmax1+1,ix2,mag(3))
            w(ixOmin1,ix2,mag(3)) = w(ixOmax1+1,ix2,mag(3))
         endif



         ! T is fixed at the boundary
         Pin = (mhd_gamma-one)*(w(ixOmax1+1,ix2,e_)&
                  -0.5d0*(w(ixOmax1+1,ix2,mag(1))**2+w(ixOmax1+1,ix2,mag(2))**2+w(ixOmax1+1,ix2,mag(3))**2) &
                  -0.5d0*(w(ixOmax1+1,ix2,mom(1))**2+w(ixOmax1+1,ix2,mom(2))**2+w(ixOmax1+1,ix2,mom(3))**2)/w(ixOmax1+1,ix2,rho_) )
         Pg1 = w(ixOmax1,ix2,rho_)*(2.0d0*Tbnd-Pin/w(ixOmax1+1,ix2,rho_))
         Pg2 = w(ixOmin1,ix2,rho_)*(4.0d0*Tbnd-3.0d0*Pin/w(ixOmax1+1,ix2,rho_))
         w(ixOmax1,ix2,e_)  = Pg1/(mhd_gamma-one) &
                            + 0.5d0*(w(ixOmax1,ix2,mag(1))**2+w(ixOmax1,ix2,mag(2))**2+w(ixOmax1,ix2,mag(3))**2) &
                            + 0.5d0*(w(ixOmax1,ix2,mom(1))**2+w(ixOmax1,ix2,mom(2))**2+w(ixOmax1,ix2,mom(3))**2)/w(ixOmax1,ix2,rho_)
         w(ixOmin1,ix2,e_)  = Pg2/(mhd_gamma-one) &
                            + 0.5d0*(w(ixOmin1,ix2,mag(1))**2+w(ixOmin1,ix2,mag(2))**2+w(ixOmin1,ix2,mag(3))**2) &
                            + 0.5d0*(w(ixOmin1,ix2,mom(1))**2+w(ixOmin1,ix2,mom(2))**2+w(ixOmin1,ix2,mom(3))**2)/w(ixOmin1,ix2,rho_)
       end do

    case (2) ! Upper radial boundary

      !!!!!!!!!!!!!!!!!!!!
      ! Outflow boundary !
      !!!!!!!!!!!!!!!!!!!!
      if(stagger_grid) then
        do ix1=ixOmin1,ixOmax1
          block%ws(ix1^%1ixO^S,1)=block%ws(ixOmin1-1^%1ixO^S,1)*block%surfaceC(ixOmin1-1^%1ixO^S,1)/block%surfaceC(ix1^%1ixO^S,1)
          block%ws(ix1^%1ixO^S,2)=block%ws(ixOmin1-1^%1ixO^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
      end if

      do ix1=ixOmin1,ixOmax1
        do ix2=ixOmin2,ixOmax2

          ! r^2*rho is continuous
          w(ix1,ix2,rho_)  = w(ixOmin1-1,ix2,rho_)*(x(ixOmin1-1,ix2,1)/x(ix1,ix2,1))**2

          ! r*r*rho*vr is continuous
          w(ix1,ix2,mom(1)) = w(ixOmin1-1,ix2,mom(1))*(x(ixOmin1-1,ix2,1)/x(ix1,ix2,1))**2

          ! rho*vtheta is continuous
          w(ix1,ix2,mom(2)) = w(ixOmin1-1,ix2,mom(2))

          ! r*vphi is continuous
          w(ix1,ix2,mom(3)) = w(ixOmin1-1,ix2,mom(3))/w(ixOmin1-1,ix2,rho_)*x(ixOmin1-1,ix2,1)/x(ix1,ix2,1)*w(ix1,ix2,rho_)

          if(.not.stagger_grid) then
            ! r*r*Br is continuous
            !w(ix1,ix2,mag(1)) = w(ixOmin1-1,ix2,mag(1))*(x(ixOmin1-1,ix2,1)/x(ix1,ix2,1))**2
            w(ix1,ix2,mag(1)) = (w(ixOmin1-1,ix2,mag(1))+Bbnd*2.0d0*dcos(x(ixOmin1-1,ix2,2))/x(ixOmin1-1,ix2,1)**3) &
                              * (x(ixOmin1-1,ix2,1)/x(ix1,ix2,1))**2 &
                              - Bbnd*2.0d0*dcos(x(ix1,ix2,2))/x(ix1,ix2,1)**3

            ! Btheta is continuous
            !w(ix1,ix2,mag(2)) = w(ixOmin1-1,ix2,mag(2))
            w(ix1,ix2,mag(2)) = w(ixOmin1-1,ix2,mag(2))+Bbnd*dsin(x(ixOmin1-1,ix2,2))/x(ixOmin1-1,ix2,1)**3-Bbnd*dsin(x(ix1,ix2,2))/x(ix1,ix2,1)**3
          end if

          ! r*Bphi is continuous
          w(ix1,ix2,mag(3)) = w(ixOmin1-1,ix2,mag(3))*x(ixOmin1-1,ix2,1)/x(ix1,ix2,1)

          ! T is continuous
          ! Ti is the temperature in the last inner cell
          Ti=(mhd_gamma-one)/w(ixOmin1-1,ix2,rho_)*(w(ixOmin1-1,ix2,e_)&
                  -0.5d0*(w(ixOmin1-1,ix2,mag(1))**2+w(ixOmin1-1,ix2,mag(2))**2+w(ixOmin1-1,ix2,mag(3))**2) &
                  -0.5d0*(w(ixOmin1-1,ix2,mom(1))**2+w(ixOmin1-1,ix2,mom(2))**2+w(ixOmin1-1,ix2,mom(3))**2)/w(ixOmin1-1,ix2,rho_) )
          w(ix1,ix2,e_) = Ti*w(ix1,ix2,rho_)/(mhd_gamma-one) &
                            + 0.5d0*(w(ix1,ix2,mag(1))**2+w(ix1,ix2,mag(2))**2+w(ix1,ix2,mag(3))**2) &
                            + 0.5d0*(w(ix1,ix2,mom(1))**2+w(ix1,ix2,mom(2))**2+w(ix1,ix2,mom(3))**2)/w(ix1,ix2,rho_)

        enddo
      enddo

    end select


  end subroutine specialbound_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: divb(ixI^S)
    w(ixO^S,nw+1)=w(ixO^S,mom(1))/w(ixO^S,rho_)
    w(ixO^S,nw+2)=w(ixO^S,mom(2))/w(ixO^S,rho_)
    w(ixO^S,nw+3)=w(ixO^S,mom(3))/w(ixO^S,rho_)
    w(ixO^S,nw+4)=w(ixO^S,mag(1))
    w(ixO^S,nw+5)=w(ixO^S,mag(2))
    w(ixO^S,nw+6)=w(ixO^S,mag(3))
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+7)=divb(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames
    varnames='vr vth vphi br bth bphi divb'
  end subroutine specialvarnames_output

  subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: sin2_thc,thc,hscale,TNow,TTarget,q0,TEquatorial,TPolar
    double precision :: Pth,r,th,theta1,theta2
    integer :: ix1, ix2

    ! Gravity
    w(ixO^S,e_)  = w(ixO^S,e_)  - qdt*wCT(ixO^S,mom(1))*19.08/((x(ixO^S,1))**2)
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) - qdt*wCT(ixO^S,rho_)*19.08/((x(ixO^S,1))**2)


    ! Heating/cooling

    q0 = 83.69d0
    theta1 = dpi*17.5d0/180.0d0
    theta2 = dpi*61.5d0/180.0d0
    ! 1.5E6 K
    TEquatorial=2.494d0
    ! 2.625E6
    TPolar=4.3645d0


    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    if(x(ix^D,1) < 30.0d0)then
       r  = x(ix^D,1)
       th = x(ix^D,2)
       Pth = (mhd_gamma-one)*(wCT(ix^D,e_)-0.5d0*(sum(wCT(ix^D,mom(:))**2)/wCT(ix^D,rho_)+sum(wCT(ix^D,mag(:))**2)))
       sin2_thc = (dsin(theta1)**2)+(dcos(theta1)**2)*(r-1.0d0)/8.0d0
       if ( (r >= 7.0d0) .and. (r < 47.0d0) ) then
          sin2_thc = (dsin(theta2)**2)+(dcos(theta2)**2)*(r-7.0d0)/40.0d0
       else if (r>= 47.0d0) then
          sin2_thc = 1.0d0
       end if
       thc = dasin(dsqrt(sin2_thc))
       hscale = 4.5d0
       TTarget = TEquatorial
       if ( (th < thc) .or. (th > (dpi-thc) ) ) then
          hscale=4.5d0*(2.0d0-(sin(th)**2)/sin2_thc)
          TTarget = TPolar
       end if
       TNow = Pth/wCT(ix^D,rho_)
       w(ix^D,e_)  = w(ix^D,e_)  + qdt*q0*wCT(ix^D,rho_)*(TTarget-TNow)*exp(-((r-1.0d0)/hscale)**2)
    endif
    {end do\}

  end subroutine specialsource

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
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    double precision :: rmin1,rmin2,rmax,tend

    if(qt < 258.73d0) then
      if(level<1) then
        if(any(x(ix^S,1) <= 20.0d0 )) then
           if(any( abs(dpi*0.5d0-x(ix^S,2)) <= (15.0d0*dpi/180.0d0) )) then
             refine=1
             coarsen=0
           else
             refine=0
             coarsen=1
           endif
        else
          if(any( abs( x(ix^S,1)*dcos(x(ix^S,2)) ) <=  7.5d0)) then
            refine=1
            coarsen=0
          else
            refine=0
            coarsen=1
          endif
        endif
      endif
    endif

    if((qt >= 258.73d0) .AND. (qt < 517.46d0)) then
      if(level<2) then
        if(any(x(ix^S,1) <= 20.0d0 )) then
           if(any( abs(dpi*0.5d0-x(ix^S,2)) <= (15.0d0*dpi/180.0d0) )) then
             refine=1
             coarsen=0
           else
             refine=0
             coarsen=1
           endif
        else
          if(any( abs( x(ix^S,1)*dcos(x(ix^S,2)) ) <=  7.5d0)) then
            refine=1
            coarsen=0
          else
            refine=0
            coarsen=1
          endif
        endif
      endif
    endif

    if(qt >= 517.46d0) then
      if(level<3) then
        if(any(x(ix^S,1) <= 20.0d0 )) then
           if(any( abs(dpi*0.5d0-x(ix^S,2)) <= (15.0d0*dpi/180.0d0) )) then
             refine=1
             coarsen=0
           else
             refine=0
             coarsen=1
           endif
        else
          if(any( abs( x(ix^S,1)*dcos(x(ix^S,2)) ) <=  7.5d0)) then
            refine=1
            coarsen=0
          else
            refine=0
            coarsen=1
          endif
        endif
      endif
    endif

  end subroutine specialrefine_grid

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixO^S,1)=2.0d0*Busr*dcos(x(ixO^S,2))/x(ixO^S,1)**3
    wB0(ixO^S,2)=Busr*dsin(x(ixO^S,2))/x(ixO^S,1)**3
    wB0(ixO^S,3)=0.d0

  end subroutine specialset_B0

end module mod_usr
