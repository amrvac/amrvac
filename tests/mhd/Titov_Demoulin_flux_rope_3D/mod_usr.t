module mod_usr
  use mod_TDfluxrope
  use mod_mhd
  implicit none
  double precision :: tcorona, rhocorona, internalrhoratio

contains                                

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_3D")
 
    unit_length        = 5.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d10! cm^-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_refine_grid     => specialrefine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters

    character(len=20) :: printsettingformat  
    double precision :: rhocorona
    double precision :: Itube, Nt_TD99
    
    printsettingformat='(1x,A20,ES15.7,A30)'
    if(mype==0) then
      write(*,*),"Dimensionless units:"
      write(*,printsettingformat),"unit_length ",unit_length,"(cm)"
      write(*,printsettingformat),"unit_velocity ",unit_velocity,"(cm s^-1)"
      write(*,printsettingformat),"unit_density ",unit_density,"(g cm^-3)"
      write(*,*),"Deduced dimensionless units:"
      write(*,printsettingformat),"unit_time ",unit_time,"(s)"
      write(*,printsettingformat),"unit_pressure ",unit_pressure,"(Ba = g cm^-1 s^-2)"
      write(*,printsettingformat),"Temperature unit",unit_temperature,"(K)"
      write(*,printsettingformat),"unit_magneticfield ",unit_magneticfield,"(G=g^{1/2} cm^{-1/2} s^-1)"
    end if
    
    Li_TD99=5.0d-1
    tcorona=1.0d6
    tcorona=tcorona/unit_temperature
    p_Bt_ratio=1.00d0
    !Default (case 101):
    rhocorona=1.0d10/unit_numberdensity
    d_TD99=5.0d9/unit_length  
    L_TD99=5.0d9/unit_length  
    R_TD99=11.0d9/unit_length 
    a_TD99=3.24179d9/unit_length 
    !According to Titov99, q=10^2 T Mm^2=10^14 T m^2=10^22 G cm^2
    !Ilia Roussev also use this in his subroutine
    !400~500G for maximum; No, only 160~200G is ok!
    q_TD99=4.0d21/unit_magneticfield/unit_length**2
    Itube=2.d0*q_TD99*L_TD99*R_TD99 &
    /(L_TD99**2+R_TD99**2) &
    /sqrt(L_TD99**2+R_TD99**2) &
    /(log(8.*R_TD99/a_TD99)-1.5+Li_TD99/2.)
    !According to Titov99, Izero=-7 TA ("T" means 10^3 G= 10^6 M= 10^9 K= 10^12 Bit).
    !Ilia Roussev also use this in his subroutine
    !However, Izero must equal to zero in order to allow a eruption, 
    !said by Roussev 2003.
    !400~500G for maximum; No, only 160~200G is ok!
    Izero_TD99=1.0d0/(unit_magneticfield*unit_length)/const_c* &
    (-28.0d11)*2.99792456d9
    
    !case 121:
    if(iprob >=121 .and. iprob <=140) then
      p_Bt_ratio=0.90d0
      d_TD99=5.0d9/unit_length  
      L_TD99=5.0d9/unit_length  
      R_TD99=11.0d9/unit_length 
      a_TD99=2.39168d9/unit_length 
      Itube=2.d0*q_TD99*L_TD99*R_TD99 &
      /(L_TD99**2+R_TD99**2) &
      /sqrt(L_TD99**2+R_TD99**2) &
      /(log(8.*R_TD99/a_TD99)-1.5+Li_TD99/2.)
    endif
    
    !case 122: very low twist
    if(iprob ==122) then
      Izero_TD99=1.0d0/(unit_magneticfield*unit_length)/const_c* &
      (-50.0d11)*2.99792456d9
    endif

    if(iprob ==1221) then
      p_Bt_ratio=0.94d0
      d_TD99=3.0d9/unit_length
      L_TD99=3.5d9/unit_length
      R_TD99=11.0d9/unit_length 
      a_TD99=2.39168d9/unit_length 
      Itube=2.d0*q_TD99*L_TD99*R_TD99 &
      /(L_TD99**2+R_TD99**2) &
      /sqrt(L_TD99**2+R_TD99**2) &
      /(log(8.*R_TD99/a_TD99)-1.5+Li_TD99/2.)
      Izero_TD99=1.0d0/(unit_magneticfield*unit_length)/const_c* &
      (-1.0d6)*2.99792456d9
    endif

    Nt_TD99=abs(Itube/Izero_TD99)*R_TD99**2/a_TD99**2
    if(mype==0) then
      Write(*,*),"Additional units:"
      Write(*,printsettingformat),"for line current I0 ",(unit_magneticfield*unit_length*const_c),"(statampere)"
      Write(*,printsettingformat),"for magnetic charge q ",unit_magneticfield*unit_length**2,"(G cm^2)"
      printsettingformat='(1x,A7,ES15.7)'
      write(*,*),"Dimensionless setup parameters of magnetic structure:"
      write(*,printsettingformat), 'd ',d_TD99
      write(*,printsettingformat), 'L ',L_TD99
      write(*,printsettingformat), 'R ',R_TD99
      write(*,printsettingformat), 'a ',a_TD99
      write(*,printsettingformat), 'Li ',Li_TD99
      write(*,printsettingformat), 'Izero ',Izero_TD99
      write(*,printsettingformat), 'q ',q_TD99
      write(*,printsettingformat), 'I tube ',Itube
      write(*,printsettingformat), 'N turns',Nt_TD99
      write(*,*),"Dimensionless setup parameters of solar atmopshere:"
      write(*,*), "uniform T and rho, without gravity stratification" 
      write(*,printsettingformat),'T corona', tcorona
      write(*,printsettingformat),'rho corona', rhocorona
    endif

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: current(ixI^S,7-2*ndir:ndir), Bf(ixI^S,1:ndir)
    double precision :: rvertical(ixI^S),rhoadius(ixI^S)
    double precision :: Itube, vperb, rhocut
    double precision :: currisosurface, ysection
    integer          :: idirmin
    logical, save :: first=.true.

    rvertical(ixO^S)=dsqrt(x(ixO^S,2)**2+(x(ixO^S,3)+d_TD99)**2)
    
    w(ixO^S,mom(:))=0.d0

    call TD99(ixI^L,ixO^L,x,Bf)
    w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
    
    w(ixO^S,rho_)= rhocorona
    
    internalrhoratio = 0.0d0

    if(iprob ==115 .or. iprob ==135 .or. iprob ==155) then
      internalrhoratio =10.0d0
    else
      internalrhoratio =(1.0d0-p_Bt_ratio)*1.0d0/mhd_adiab
    endif

    Itube=2.d0*q_TD99*L_TD99*R_TD99 &
          /(L_TD99**2+R_TD99**2) &
          /sqrt(L_TD99**2+R_TD99**2) &
          /(log(8.*R_TD99/a_TD99)-1.5+Li_TD99/2.)

    rhoadius(ixO^S)=dsqrt(x(ixO^S,1)**2+(rvertical(ixO^S)-R_TD99)**2)

    w(ixO^S,rho_)=rhocorona &
                 +internalrhoratio &
                 *0.5d0 &
                 *(Izero_TD99*2.0d0/sqrt(4.0d0*dpi))**2 &
                 *(2.0d0/(a_TD99**2) &
                 *Itube**2/(Izero_TD99**2) &
                 *max(0.0d0,1.0d0-rhoadius(ixO^S)**2/a_TD99**2))
    
    !Default case for torok04.fig1, middle, torok04.fig3
    !case 111,131,151: velocity pertubation
    if((iprob ==111 .or. iprob ==131 .or. iprob ==151)) then
      currisosurface= 20.0d0
      ysection=2.4d0
      vperb=1.0d0
      rhocut=2.0d0
    endif

    !case 112,132,152: velocity pertubation
    if((iprob ==112 .or. iprob ==132 .or. iprob ==152)) then
      currisosurface= 20.0d0
      ysection=2.4d0
      vperb=1.0d0
      rhocut=rhocorona+5.0d0
    endif
    !case 113,133,153: velocity pertubation
    if((iprob ==113 .or. iprob ==133 .or. iprob ==153)) then
      currisosurface= 20.0d0
      ysection=2.4d0
      vperb=-1.0d0
      rhocut=2.0d0
    endif
    !case 115,135,155: velocity pertubation
    if((iprob ==115 .or. iprob ==135 .or. iprob ==155)) then
      currisosurface= 20.0d0
      ysection=0.9d0
      vperb=0.3d0
      rhocut=rhocorona+50.0d0
    endif
    !velocity perbutation
    if((iprob ==111 .or. iprob ==131 .or. iprob ==151 .or. &
        iprob ==112 .or. iprob ==132 .or. iprob ==152 .or. &
        iprob ==113 .or. iprob ==133 .or. iprob ==153 .or. &
        iprob ==115 .or. iprob ==135 .or. iprob ==155)) then
       call get_current(w,ixI^L,ixO^L,idirmin,current)
       w(ixO^S,mom(1))=0.d0
       w(ixO^S,mom(2))=0.d0
       w(ixO^S,mom(3))=vperb* &
       max(0.0d0,sign(1.0d0,w(ixO^S,rho_)-rhocut))* &
       max(0.0d0,sign(1.0d0,ysection-abs(x(ixO^S,2))))
    endif
    !case 115, 135,155: unifrom alfven speed
    if(iprob ==115 .or. iprob ==135 .or. iprob ==155) then
      w(ixO^S,rho_)=0.5d0*sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    endif

    if(mhd_energy) w(ixO^S,p_)=rhocorona*tcorona

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bf(ixI^S,1:ndir)
    double precision :: pth(ixI^S)
    double precision :: rvertical(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision :: rhoadius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision :: Itube

    integer:: ixIM^L,ix^D

    select case(iB)
    case(1)
      ! min boundary in the 1st dimension 
      ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;
      ixIMmin2=ixOmin2;ixIMmax2=ixOmax2;
      ixIMmin1=ixOmax1+1;ixIMmax1=ixOmax1+1;
      call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
      w(ixO^S,:)=0.d0
      do ix1=ixOmax1,ixOmin1,-1
      ! simply copy in the edge values (zero gradient type, all ghost cells layers identical)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)= w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)  = pth(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      enddo
      call TD99(ixI^L,ixO^L,x,Bf)
      w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(2)
      ! max boundary in the 1st dimension
      ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;
      ixIMmin2=ixOmin2;ixIMmax2=ixOmax2;
      ixIMmin1=ixOmin1-1;ixIMmax1=ixOmin1-1;
      call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
      w(ixO^S,:)=0.d0
      do ix1=ixOmin1,ixOmax1,+1
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)= w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)  = pth(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      enddo
      call TD99(ixI^L,ixO^L,x,Bf)
      w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(3)
      ! min boundary in the 2nd dimension 
      ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+1;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;
      call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
      w(ixO^S,:)=0.d0
      do ix2=ixOmax2,ixOmin2,-1
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)= w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,rho_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,p_)  = pth(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3)
      enddo
      call TD99(ixI^L,ixO^L,x,Bf)
      w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ! max boundary in the 2nd dimension
      ixIMmin2=ixOmin2-1;ixIMmax2=ixOmin2-1;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;
      call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
      w(ixO^S,:)=0.d0
      do ix2=ixOmin2,ixOmax2,+1
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)= w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,rho_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,p_)  = pth(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3)
      enddo
      call TD99(ixI^L,ixO^L,x,Bf)
      w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(5)
      ! min boundary in the 3rd dimension
      ixIMmin3=ixOmax3+1;ixIMmax3=ixOmax3+1;
      ixIMmin2=ixOmin2;ixIMmax2=ixOmax2;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
      w(ixO^S,:)=0.d0
      do ix3=ixOmax3,ixOmin3,-1
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,p_)  = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+1)
      enddo
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))&
                   /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
      w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(2))&
                   /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
      w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(3))&
                   /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)  
      Itube=2.d0*q_TD99*L_TD99*R_TD99 &
      /(L_TD99**2+R_TD99**2) &
      /sqrt(L_TD99**2+R_TD99**2) &
      /(log(8.*R_TD99/a_TD99)-1.5+Li_TD99/2.)
      internalrhoratio =(1.0d0-p_Bt_ratio)*1.0d0/mhd_adiab
      do ix3=ixOmax3,ixOmin3,-1
        rvertical(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,2)**2 &
        +(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,3)+d_TD99)**2)
        rhoadius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,1)**2 &
        +(rvertical(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-R_TD99)**2)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,rho_)= rhocorona &
        +internalrhoratio &
        *0.5d0 &
        *(Izero_TD99*2.0d0/sqrt(4.0d0*dpi))**2 &
        *(2.0d0/(a_TD99**2) &
        *Itube**2/(Izero_TD99**2) &
        *max(0.0d0,1.0d0 &
        -rhoadius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2 &
        /a_TD99**2))
      enddo
      call TD99(ixI^L,ixO^L,x,Bf)
      w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(6)
      ! max boundary in the 3rd dimension
      ixIMmin3=ixOmin3-1;ixIMmax3=ixOmin3-1;
      ixIMmin2=ixOmin2;ixIMmax2=ixOmax2;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      call mhd_get_pthermal(w,x,ixI^L,ixIM^L,pth)
      w(ixO^S,:)=0.d0
      do ix3=ixOmin3,ixOmax3,+1
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,rho_)= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1,rho_)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,p_)  = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1)
      enddo
      call TD99(ixI^L,ixO^L,x,Bf)
      w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    end select

  end subroutine specialbound_usr

  subroutine specialrefine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
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
    
    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    integer          :: idirmin

    if(iprob==1221) then
      refine = -1
      coarsen = -1
    
      if(level==1) then
        if(any( x(ixI^S,3)<=2.6666666667d0+qt/1.3333333333d0/3.6d0 .and. &
        (x(ixI^S,1)**2/1.3333333d0**2+x(ixI^S,2)**2/4.0d0**2)<=1.0d0 &
        ))then
          refine = 1
          coarsen = 0
        end if
      end if
      if(level==2 .or. level==3) then
        if(any( x(ixI^S,3)<=2.6666666667d0+qt/1.3333333333d0/3.60d0 .and. &
        (x(ixI^S,1)**2/1.3333333d0**2+x(ixI^S,2)**2/4.0d0**2)<=1.0d0 ))then
          refine = 1
          coarsen = 0
        end if
      end if
      if(level==4 .and. qt>=(4.0d0-0.1d0)) then
        if(any( x(ixI^S,3)<=0.3333333333333d0+qt/1.3333333333d0/3.6d0 .and. &
        x(ixI^S,3)>0.333333333333333d0 .and. &
        (x(ixI^S,1)**2/0.0833333333333d0**2+x(ixI^S,2)**2/1.6666666d0**2)<=1.0d0 &
        ))then
          refine = 0
          coarsen = 0
        end if
      end if
      if(level==5 .and. qt>=4.0d0) then
        if(any( x(ixI^S,3)<=0.3333333333333d0+qt/1.3333333333d0/3.6d0 .and. &
        x(ixI^S,3)>0.5d0 .and. &
        (x(ixI^S,1)**2/0.0416666666667d0**2+x(ixI^S,2)**2/1.6666666d0**2)<=1.0d0 &
        ))then
          refine = 0
          coarsen = 0
        end if
      end if
    endif
  
  end subroutine specialrefine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    
    double precision :: Btoro(ixI^S,1:ndir)
    double precision :: Bpolo(ixI^S,1:ndir)
    double precision :: pth(ixI^S),divb(ixI^S)
    
    double precision :: current(ixI^S,7-2*ndir:3)
    double precision :: rvertical(ixI^S),rhoadius(ixI^S)
    
    integer :: idir,jdir,kdir,idirmin

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)

    call get_current(w,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+2)=sqrt(sum(current(ixO^S,:)**2,dim=ndim+1))

    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=divb(ixO^S)

    Btoro(ixI^S,:)=0.d0
    Bpolo(ixI^S,:)=0.d0
    call TD99_BI(ixI^L,ixO^L,x,Btoro)
    call TD99_Bth(ixI^L,ixO^L,x,Bpolo)

    rvertical(ixI^S)=dsqrt(x(ixI^S,2)**2+(x(ixI^S,3)+d_TD99)**2)
    rhoadius(ixI^S)=dsqrt(x(ixI^S,1)**2+(rvertical(ixI^S)-R_TD99)**2)

    w(ixO^S,nw+4)= rhoadius(ixO^S) & 
            *w(ixI^S,mag(2))/R_TD99 &
            /abs(w(ixI^S,mag(1))) &
            *dpi/acos(d_TD99/R_TD99)
    
    w(ixO^S,nw+5:nw+7)=current(ixO^S,1:3)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) varnames

    varnames='T j divb safq j1 j2 j3'

  end subroutine specialvarnames_output

end module mod_usr
