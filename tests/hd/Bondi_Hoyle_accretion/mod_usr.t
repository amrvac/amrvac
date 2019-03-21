module mod_usr
  use mod_hd
  implicit none
  double precision :: mach, rho0, vel, Racc, gm, pini

contains

  subroutine usr_init()
    call set_coordinate_system('spherical_2D')

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_source        => pt_grav_source
    usr_get_dt        => get_dt_pt_grav
    usr_special_bc    => specialbound_usr
    usr_internal_bc   => intern_bc

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    
    hd_gamma=5.d0/3.d0
    mach=4.d0
    vel=one
    gm=half
    rho0=one
    Racc=one
    pini=((1/hd_gamma)*rho0*vel**two)/mach**two

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid 

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    w(ix^S,rho_)   =  rho0
    w(ix^S,mom(1)) = -rho0*vel*dcos(x(ix^S,2))
    w(ix^S,mom(2)) =  rho0*vel*dsin(x(ix^S,2))
    w(ix^S,e_)     =  pini/(hd_gamma-one)+half*rho0*vel**two

  end subroutine initonegrid_usr

  subroutine pt_grav_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    
    w(ixO^S,mom(1))=w(ixO^S,mom(1)) - &
         qdt * wCT(ixO^S,rho_)    * (gm/x(ixO^S,1)**two)
    w(ixO^S,e_ )=w(ixO^S,e_ ) - &
         qdt * wCT(ixO^S,mom(1))  * (gm/x(ixO^S,1)**two)

  end subroutine pt_grav_source

  subroutine get_dt_pt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
    double precision, intent(in) :: w(ixG^S,1:nw)
    double precision, intent(inout) :: dtnew

    dtnew=bigdouble
    dtnew=min(dtnew,minval(dsqrt(block%dx(ix^S,1)/(gm/x(ix^S,1)**two))))

  end subroutine get_dt_pt_grav

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision :: rho(ixG^S)
    
    ! Cont grad
      rho(ixBmax1,ixB^LIM^DE:)=w(ixBmax1+1,ixB^LIM^DE:,rho_)+&
         (w(ixBmax1+2,ixB^LIM^DE:,rho_)-w(ixBmax1+1,ixB^LIM^DE:,rho_))*&
         ((x(ixBmax1      ,ixB^LIM^DE:,1)-x(ixBmax1+1,ixB^LIM^DE:,1))/&
         (x(ixBmax1+2,ixB^LIM^DE:,1)-x(ixBmax1+1,ixB^LIM^DE:,1)))
      w(ixBmax1,ixB^LIM^DE:,rho_)=rho(ixBmax1,ixB^LIM^DE:)
      ! rho*vr*r^2
      w(ixBmax1,ixB^LIM^DE:,mom(1)) = min(zero,w(ixBmax1+1,ixB^LIM^DE:,mom(1))*&
           ( x(ixBmax1+1,ixB^LIM^DE:,1)/x(ixBmax1,ixB^LIM^DE:,1) )**two)
      ! rho*vr*vth*r^3 cont (ie r*vth cont)
      w(ixBmax1,ixB^LIM^DE:,mom(2)) = rho(ixBmax1,ixB^LIM^DE:) * &
         (w(ixBmax1+1,ixB^LIM^DE:,mom(2))/w(ixBmax1+1,ixB^LIM^DE:,rho_))*(x(ixBmax1+1,ixB^LIM^DE:,1)/x(ixBmax1,ixB^LIM^DE:,1))
      ! (e+P+rho*phi)*v_r*r^2    cont
      w(ixBmax1,ixB^LIM^DE:,e_) = (one/hd_gamma)*&
         ( (hd_gamma-one)*half*&
            ( ( (min(zero,w(ixBmax1+1,ixB^LIM^DE:,mom(1))*&
                 ( x(ixBmax1+1,ixB^LIM^DE:,1)/x(ixBmax1,ixB^LIM^DE:,1) )**two))**two+&
                (rho(ixBmax1,ixB^LIM^DE:) * &
                   (w(ixBmax1+1,ixB^LIM^DE:,mom(2))/w(ixBmax1+1,ixB^LIM^DE:,rho_))*(x(ixBmax1+1,ixB^LIM^DE:,1)/x(ixBmax1,ixB^LIM^DE:,1)))**two )&
                /rho(ixBmax1,ixB^LIM^DE:))+ &
               rho(ixBmax1,ixB^LIM^DE:)*(gm/dabs(x(ixBmax1,ixB^LIM^DE:,1))) +&
               (rho(ixBmax1,ixB^LIM^DE:)/w(ixBmax1+1,ixB^LIM^DE:,rho_))*&
                (hd_gamma*w(ixBmax1+1,ixB^LIM^DE:,e_)-(hd_gamma-one)*half*&
                ((w(ixBmax1+1,ixB^LIM^DE:,mom(1))**two+w(ixBmax1+1,ixB^LIM^DE:,mom(2))**two)/w(ixBmax1+1,ixB^LIM^DE:,rho_)) -&
                w(ixBmax1+1,ixB^LIM^DE:,rho_)*(gm/dabs(x(ixBmax1+1,ixB^LIM^DE:,1))) ) )

   w(ixBmin1,ixB^LIM^DE:,:)=w(ixBmax1,ixB^LIM^DE:,:)

  end subroutine specialbound_usr

  subroutine intern_bc(level,qt,ixI^L,ixO^L,w,x)

    integer, intent(in) :: ixI^L,ixO^L,level
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: R1, R2, rho(ixI^S), zeta(ixI^S)

    R1=half*Racc ! xprobmax1/(335./20.)
    R2=one *Racc ! xprobmax1/(335./35.)
    zeta(ixO^S)=half*x(ixO^S,1)*dsin(x(ixO^S,2))*&
         ( one + dsqrt(one+4.d0*(gm/vel**two)*&
         ((one-dcos(x(ixO^S,2)))/(x(ixO^S,1)*dsin(x(ixO^S,2))**two))) )
    rho(ixO^S)=rho0*&
         (zeta(ixO^S)**two/&
         (x(ixO^S,1)*dsin(x(ixO^S,2))*(two*zeta(ixO^S)-x(ixO^S,1)*dsin(x(ixO^S,2)))))
    ! Is it ahead the shock? THEN FIX! No need to compare to B-K solution
    where ( x(ixO^S,2)<dpi/two .AND. x(ixO^S,1)>R2/(1+((R2-R1)/R1)*dcos(x(ixO^S,2))) )
       w(ixO^S,rho_)= rho(ixO^S)
       w(ixO^S,mom(1)) = -rho0*&
            (zeta(ixO^S)**two/&
            (x(ixO^S,1)*dsin(x(ixO^S,2))*(two*zeta(ixO^S)-x(ixO^S,1)*dsin(x(ixO^S,2)))))*&
            dsqrt(vel**two+(two*gm)/x(ixO^S,1)-&
            ((zeta(ixO^S)*vel)/x(ixO^S,1))**two)
       w(ixO^S,mom(2)) = rho0*&
            (zeta(ixO^S)**two/&
            (x(ixO^S,1)*dsin(x(ixO^S,2))*(two*zeta(ixO^S)-x(ixO^S,1)*dsin(x(ixO^S,2)))))*&
            ((zeta(ixO^S)*vel)/x(ixO^S,1))
       w(ixO^S,e_)  = ((rho0*vel**two)/hd_gamma)*&
            (one/(mach**two*(hd_gamma-one))+half)*&
            (zeta(ixO^S)**two/&
            (x(ixO^S,1)*dsin(x(ixO^S,2))*(two*zeta(ixO^S)-x(ixO^S,1)*dsin(x(ixO^S,2)))))+&
            (rho0/hd_gamma)*&
            (zeta(ixO^S)**two/&
            (x(ixO^S,1)*dsin(x(ixO^S,2))*(two*zeta(ixO^S)-x(ixO^S,1)*dsin(x(ixO^S,2)))))*&
            (half*(hd_gamma-one)* (dsqrt(vel**two+(two*gm)/x(ixO^S,1)-&
            ((zeta(ixO^S)*vel)/x(ixO^S,1))**two)**two+&
            ((zeta(ixO^S)*vel)/x(ixO^S,1))**two) + gm/x(ixO^S,1))
    endwhere

  end subroutine intern_bc

end module mod_usr
