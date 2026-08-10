module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  implicit none
  double precision :: q_e, parb,unit_currentdensity
  double precision :: sheet_ymin=5.d0
  double precision :: sheet_ymax=12.d0
  double precision :: sheet_core_halfwidth=0.4d0
  double precision :: sheet_halo_halfwidth=2.d0
  double precision :: evaporation_height=0.5d0
  double precision :: reconnection_xmax=0.6d0
  double precision :: reconnection_ymin=4.d0
  double precision :: reconnection_ymax=10.d0
  double precision :: eta_fixed_amplitude=2.d-4

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
    usr_set_parameters       => set_usr_parameters
    usr_print_log            => flare_transport_log

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

  subroutine set_usr_parameters()
    integer :: n
    namelist /usr_list/ sheet_ymin,sheet_ymax,sheet_core_halfwidth, &
         sheet_halo_halfwidth,evaporation_height,reconnection_xmax, &
         reconnection_ymin,reconnection_ymax,eta_fixed_amplitude

    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,usr_list,end=111)
111   close(unitpar)
    end do
    if(sheet_ymax<=sheet_ymin) call mpistop('sheet_ymax must exceed sheet_ymin')
    if(sheet_core_halfwidth<=zero) call mpistop('sheet core half-width must be positive')
    if(sheet_halo_halfwidth<=sheet_core_halfwidth) &
      call mpistop('sheet halo half-width must exceed core half-width')
    if(eta_fixed_amplitude<=zero) &
      call mpistop('eta_fixed_amplitude must be positive')
  end subroutine set_usr_parameters

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
    call eos%to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qdt,qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qdt,qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: pth(ixI^S),Qp(ixI^S)
    integer :: ix^D, ixA^L, ixOs^L, idir

    if(mhd_glm) w(ixO^S,psi_)=0.d0
    select case(iB)
    case(1)
      ixA^L=ixO^L;
      ixAmin1=ixOmax1+1;ixAmax1=ixOmax1+nghostcells;
      call eos%get_thermal_pressure(w,x,ixI^L,ixA^L,pth)
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
      call eos%to_conserved(ixI^L,ixO^L,w,x)
    case(2)
      ixA^L=ixO^L;
      ixAmin1=ixOmin1-nghostcells;ixAmax1=ixOmin1-1;
      call eos%get_thermal_pressure(w,x,ixI^L,ixA^L,pth)
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
      call eos%to_conserved(ixI^L,ixO^L,w,x)
    case(3)
      ixA^L=ixO^L;
      ixAmin2=ixOmax2+1;ixAmax2=ixOmax2+1;
      call eos%get_thermal_pressure(w,x,ixI^L,ixA^L,pth)
      do ix2=ixOmin2,ixOmax2
        w(ix2^%2ixO^S,rho_)=w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(1))=w(ixOmax2+1^%2ixO^S,mom(1))/w(ixOmax2+1^%2ixO^S,rho_)
        !w(ix2^%2ixO^S,mom(2))=w(ixOmax2+1^%2ixO^S,mom(2))/w(ixOmax2+1^%2ixO^S,rho_)
        !w(ix2^%2ixO^S,mom(3))=w(ixOmax2+1^%2ixO^S,mom(3))/w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,p_)=pth(ixOmax2+1^%2ixO^S)
        if(mhd_hyperbolic_tc) w(ix2^%2ixO^S,qpar_)=w(ixOmax2+1^%2ixO^S,qpar_)
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
      call eos%to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ixA^L=ixO^L;
      ixAmin2=ixOmin2-1;ixAmax2=ixOmin2-1;
      call eos%get_thermal_pressure(w,x,ixI^L,ixA^L,pth)
      do ix2=ixOmin2,ixOmax2
        w(ix2^%2ixO^S,rho_)=w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(1))=w(ixOmin2-1^%2ixO^S,mom(1))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(2))=w(ixOmin2-1^%2ixO^S,mom(2))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(3))=w(ixOmin2-1^%2ixO^S,mom(3))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)
        if(mhd_hyperbolic_tc) w(ix2^%2ixO^S,qpar_)=w(ixOmin2-1^%2ixO^S,qpar_)
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
      call eos%to_conserved(ixI^L,ixO^L,w,x)
    case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  subroutine p_for_errest(ixI^L,ixO^L,iflag,w,x,var)
    integer, intent(in)           :: ixI^L,ixO^L,iflag
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: var(ixI^S)

    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,var)
    
  end subroutine p_for_errest

  subroutine flare_transport_log()
    integer, parameter :: log_unit=123
    integer :: iigrid,igrid,idirmin
    double precision :: dvolume(ixG^T),pth(ixG^T),excess(ixG^T)
    double precision :: current(ixG^T,7-2*ndir:3),eta(ixG^T),vdrift(ixG^T)
    double precision :: local_sum(5),global_sum(5)
    double precision :: local_max(3),global_max(3)
    double precision :: halo_fraction,width_rms,cell_dy
    logical :: in_sheet(ixG^T),in_core(ixG^T),in_halo(ixG^T)
    logical :: at_evaporation_height(ixG^T),in_reconnection_box(ixG^T)
    logical, save :: first_call=.true.
    character(len=std_len) :: filename

    local_sum=zero
    local_max=-bigdouble
    do iigrid=1,igridstail
      igrid=igrids(iigrid)
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      if(slab) then
        dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*}
      else
        dvolume(ixM^T)=block%dvolume(ixM^T)
      end if

      call eos%get_thermal_pressure(ps(igrid)%w,ps(igrid)%x, &
           ixG^LL,ixM^LL,pth)
      excess(ixM^T)=max((pth(ixM^T)-one)*eos%inv_gamma_minus_1,zero)
      in_sheet(ixM^T)=ps(igrid)%x(ixM^T,2)>sheet_ymin .and. &
           ps(igrid)%x(ixM^T,2)<sheet_ymax .and. &
           abs(ps(igrid)%x(ixM^T,1))<sheet_halo_halfwidth
      in_core(ixM^T)=in_sheet(ixM^T) .and. &
           abs(ps(igrid)%x(ixM^T,1))<sheet_core_halfwidth
      in_halo(ixM^T)=in_sheet(ixM^T) .and. .not.in_core(ixM^T)

      local_sum(1)=local_sum(1)+sum(excess(ixM^T)*dvolume(ixM^T),mask=in_sheet(ixM^T))
      local_sum(2)=local_sum(2)+sum(excess(ixM^T)*dvolume(ixM^T),mask=in_core(ixM^T))
      local_sum(3)=local_sum(3)+sum(excess(ixM^T)*dvolume(ixM^T),mask=in_halo(ixM^T))
      local_sum(4)=local_sum(4)+sum(ps(igrid)%x(ixM^T,1)**2* &
           excess(ixM^T)*dvolume(ixM^T),mask=in_sheet(ixM^T))
      if(any(in_sheet(ixM^T))) local_max(1)=max(local_max(1), &
           maxval(pth(ixM^T)/max(ps(igrid)%w(ixM^T,rho_),smalldouble), &
           mask=in_sheet(ixM^T)))

      cell_dy=rnode(rpdx2_,igrid)
      at_evaporation_height(ixM^T)= &
           abs(ps(igrid)%x(ixM^T,2)-evaporation_height)<=half*cell_dy
      local_sum(5)=local_sum(5)+sum(max(ps(igrid)%w(ixM^T,mom(2)),zero)* &
           dvolume(ixM^T)/cell_dy,mask=at_evaporation_height(ixM^T))

      call get_current(ps(igrid)%w,ixG^LL,ixM^LL,idirmin,current)
      call special_eta(ps(igrid)%w,ixG^LL,ixM^LL,idirmin, &
           ps(igrid)%x,current,eta)
      vdrift(ixM^T)=dsqrt(sum(current(ixM^T,:)**2,dim=ndim+1))/ &
           max(ps(igrid)%w(ixM^T,rho_)*q_e,smalldouble)
      in_reconnection_box(ixM^T)= &
           abs(ps(igrid)%x(ixM^T,1))<reconnection_xmax .and. &
           ps(igrid)%x(ixM^T,2)>reconnection_ymin .and. &
           ps(igrid)%x(ixM^T,2)<reconnection_ymax
      if(any(in_reconnection_box(ixM^T))) local_max(2)=max(local_max(2), &
           maxval(abs(eta(ixM^T)*current(ixM^T,3)), &
           mask=in_reconnection_box(ixM^T)))
      if(any(in_reconnection_box(ixM^T))) local_max(3)=max(local_max(3), &
           maxval(vdrift(ixM^T),mask=in_reconnection_box(ixM^T)))
    end do

    call MPI_ALLREDUCE(local_sum,global_sum,5,MPI_DOUBLE_PRECISION, &
         MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(local_max,global_max,3,MPI_DOUBLE_PRECISION, &
         MPI_MAX,icomm,ierrmpi)
    if(global_sum(1)>smalldouble) then
      halo_fraction=global_sum(3)/global_sum(1)
      width_rms=dsqrt(global_sum(4)/global_sum(1))
    else
      halo_fraction=zero
      width_rms=zero
    end if
    if(global_max(1)<-half*bigdouble) global_max(1)=zero
    if(global_max(2)<-half*bigdouble) global_max(2)=zero
    if(global_max(3)<-half*bigdouble) global_max(3)=zero

    if(mype==0) then
      filename=trim(base_filename)//'.log'
      if(first_call .and. restart_from_file==undefined) then
        open(log_unit,file=trim(filename),status='replace')
        write(log_unit,'(a)') &
             '# time E_excess E_core E_halo halo_fraction width_rms Tmax '// &
             'Mdot_up Erec_max vdrift_max'
        close(log_unit)
      end if
      open(log_unit,file=trim(filename),status='old',position='append')
      write(log_unit,'(10(es18.10,1x))') global_time,global_sum(1), &
           global_sum(2),global_sum(3),halo_fraction,width_rms, &
           global_max(1),global_sum(5),global_max(2),global_max(3)
      close(log_unit)
      first_call=.false.
    end if
  end subroutine flare_transport_log

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
    double precision :: pth(ixI^S),B2(ixI^S),divb(ixI^S),Te(ixI^S)
    double precision :: kperp_ratio(ixI^S),ne(ixI^S),chi(ixI^S),Cchi
    double precision, parameter :: xe_prefac_cgs=4.753567596681522d6
    double precision :: Btotal(ixI^S,1:ndir),current_o(ixI^S,3)
    integer :: idir,idirmin

    ! output temperature
    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)
    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)
    w(ixO^S,nw+1)=Te(ixO^S)
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
    w(ixO^S,nw+9)=dsqrt(B2(ixO^S))
    kperp_ratio(ixO^S)=zero
    if(mhd_hyperbolic_tc .and. mhd_hyperbolic_tc_use_perp) then
      select case(trim(mhd_hyperbolic_tc_perp_mode))
      case('fixed_reference')
        kperp_ratio(ixO^S)=mhd_hyperbolic_tc_kappa_perp_factor
      case('weak_field_isotropization')
        kperp_ratio(ixO^S)=mhd_hyperbolic_tc_Bmin**2/ &
             (B2(ixO^S)+mhd_hyperbolic_tc_Bmin**2)
      case('electron_magnetization')
        Cchi=(xe_prefac_cgs/mhd_hyperbolic_tc_coulomb_log)* &
             unit_magneticfield*unit_temperature**1.5d0/unit_numberdensity
        ne(ixO^S)=w(ixO^S,rho_)*(one+2.d0*eos%He_abundance)
        chi(ixO^S)=Cchi*dsqrt(B2(ixO^S))* &
             max(Te(ixO^S),smalldouble)**1.5d0/max(ne(ixO^S),smalldouble)
        kperp_ratio(ixO^S)=one/(one+chi(ixO^S)**2)
      end select
    end if
    w(ixO^S,nw+10)=kperp_ratio(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te Alfv divB beta j1 j2 j3 eta Bmag kperp_ratio'
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
    double precision, intent(in) :: current(ixI^S,7-2*ndir:3)
    double precision, intent(out) :: eta(ixI^S)
    double precision :: rad(ixI^S),heta,reta
 
    heta = 6.d0
    reta = 0.8d0 * 0.3d0

    ! Prescribed, time-independent localized resistivity.  Reconnection is a
    ! controlled heat source here; anomalous-resistivity onset is not part of
    ! the perpendicular-transport experiment.
    rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2)-heta)**2)
    where (rad(ixO^S) .lt. reta)
      eta(ixO^S)=eta_fixed_amplitude* &
           (2.d0*(rad(ixO^S)/reta)**3-3.d0*(rad(ixO^S)/reta)**2+1.d0)
    elsewhere
      eta(ixO^S)=zero
    endwhere

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
      call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,tmp)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D))  integral_grid=integral_grid+tmp(ix^D)/(eos%gamma-1.d0)*dvolume(ix^D)
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
