!> V-B driven simulation using the TMF module
!> Jinhan Guo @ Nanjng University 2022.07
module mod_usr
  use mod_mf

  implicit none
  integer, parameter :: jmax=8000
  double precision, allocatable :: pbc(:),rbc(:)
  ! 1D solar atmosphere table for pressure, density, and height
  double precision :: pa(jmax),ra(jmax),ya(jmax)
  double precision :: usr_grav,SRadius,rhob,Tiso,dr,gzone,bQ0,lQ0,startpos^D,tdstop
  double precision, allocatable :: lalpha,llift 
  logical, save :: firstusrglobaldata=.true.

 ! variables for boundary condition
  integer, save :: nx_dd^D
  double precision, save :: dxa^D
  double precision, allocatable, save :: Bddb(:,:,:,:),Vddb(:,:,:,:),t_driven_b(:),t_driven_v(:)
  double precision, save :: dtbv,t0b,t0v
  double precision, save :: x00,y00,phi00,lat00,B00
  integer, save :: nxbc^D,ngl,np
  integer :: i,j,k

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
    usr_refine_grid     => special_refine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_refine_threshold=> specialthreshold
    usr_write_analysis => record_force_free_metrics
    usr_transform_w    => transform_w_fill_plasma
    usr_special_resistivity => specialres_usr

    call set_coordinate_system("Cartesian_3D")
    call mf_activate()

  end subroutine usr_init

  !==============================================================================
  ! Purpose: to initialize user public parameters and reset global parameters.
  !          Input data are also read here.
  !==============================================================================
  subroutine initglobaldata_usr()
    double precision :: x0, y0, z0, r0, h, theta, L, Bperp
    integer :: ixp, nphalf
    character(len=80) :: filename
    integer :: file_handle
    integer, dimension(MPI_STATUS_SIZE) :: statuss
    double precision :: driven_fast_times


    driven_fast_times=12.0d0
    ! read vector magnetic field and vector velocity field
    if(mype==0) then
      filename='data_ibc/B_boundary.dat'
      open(unit=23,file=filename,status='old',access='stream')
      read(23) nx_dd^D
    endif
    if(npe>1) then
      {call MPI_BCAST(nx_dd^D,1,MPI_INTEGER,0,icomm,ierrmpi)\}
    endif
    allocate(Bddb(nx_dd^D,ndir))
    allocate(Vddb(nx_dd^D,ndir))
    allocate(t_driven_v(nx_dd3))
    allocate(t_driven_b(nx_dd3))
    if(mype==0) then
      read(23) dxa^D
      read(23) Bddb
      read(23) t_driven_b 
      close(23)
    endif
    if(npe>1) then
     {call MPI_BCAST(dxa^D,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)\}
      call MPI_BCAST(Bddb,nx_dd1*nx_dd2*nx_dd3*3,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(t_driven_b,nx_dd3,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif


    if(mype==0) then
      filename='data_ibc/V_boundary.dat'
      open(unit=23,file=filename,status='old',access='stream')
      read(23) nx_dd^D
    endif
    if(npe>1) then
      {call MPI_BCAST(nx_dd^D,1,MPI_INTEGER,0,icomm,ierrmpi)\}
    endif
    if(mype==0) then
      read(23) dxa^D
      read(23) Vddb
      read(23) t_driven_b 
      close(23)
    endif
    if(npe>1) then
     {call MPI_BCAST(dxa^D,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)\}
      call MPI_BCAST(Vddb,nx_dd1*nx_dd2*nx_dd3*3,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(t_driven_b,nx_dd3,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif

    dxa^D=dxa^D*1.d5/unit_length;  ! The unit of input dxa^D should be km.
    Bddb=Bddb/unit_magneticfield
    ! Note that the unit of Vddb was km/s, so it times 1.0d5 to change to cm/s.
    Vddb=Vddb*1.0d5/unit_velocity

    do i=1,nx_dd3
      t_driven_v(i)=t_driven_b(i)-360.0d0/driven_fast_times
      t_driven_b(i)=t_driven_b(i)/unit_time
      t_driven_v(i)=t_driven_v(i)/unit_time
    enddo
    
   if (mype==0) then
     print *, 'driven data (number)', nx_dd1,nx_dd2,nx_dd3 
     print *, 'driven time (physical [hr], alfven)', maxval(t_driven_b)*unit_time*driven_fast_times/3600.0d0, maxval(t_driven_b)
     print *, 'save dt: [hr] [alfven time]', 1.0d0,3600.0d0/(unit_time*driven_fast_times)
     print *, 'dt[phy min],dt[aflver]',24,24.0d0*60.0d0/(unit_time*driven_fast_times)
     print *, 'Vddv', minval(Vddb*unit_velocity/1.0d5), maxval(Vddb*unit_velocity/1.0d5)
     print *, 'Bddb', minval(Bddb*unit_magneticfield),maxval(Bddb*unit_magneticfield)

   endif


  end subroutine initglobaldata_usr
 
  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bf(ixI^S,1:ndir)
    logical, save :: first=.true.
    if(mype==0 .and. first) then
      write(*,*)'initializing grids ...'
      first=.false.
    endif
   w(ixO^S,mom(1:3)) = 0.0d0
   if(mf_glm)  w(ixO^S,psi_)=0.d0
  end subroutine initonegrid_usr

  subroutine transform_w_fill_plasma(ixI^L,ixO^L,nw_in,w_in,x,w_out)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, nw_in
    double precision, intent(in)  :: w_in(ixI^S,1:nw_in)
    double precision, intent(in)  :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: w_out(ixI^S,1:nw)

    double precision :: res,Bf(ixI^S,1:ndir),xC(ixI^S,1:ndim)
    integer :: ix^D,na,idir,ixC^L
    w_out=0.d0
    ! only preserve magnetic field
    w_out(ixO^S,mag(1:3)) = w_in(ixO^S,5:7)
    w_out(ixO^S,mom(1:3)) = 0.0d0
    if(mf_glm)  w_out(ixO^S,psi_)=0.d0
  end subroutine transform_w_fill_plasma


  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: tmp1(ixI^S),tmp2(ixI^S),pth(ixI^S),Qp(ixI^S)
    double precision :: xlen^D,coeffrho,ft,tfstop,tramp1,tramp2,vlimit
    integer :: idir,ix^D,ixIM^L,ixOs^L,ixbc^D,jxO^L

    integer :: af,ixb,ixv
    double precision :: dxb^D
    character(len=80) :: v_type,b_type

    if(mf_glm)  w(ixO^S,psi_)=0.d0


    dxb^D=dx(^D,refine_max_level);
   ! print *,'max_level, dxb1, min1', refine_max_level, dxb1,&
   !          (xprobmax1-xprobmin1)/520.0d0
    af=nint(dxlevel(1)/dxb1)
    select case(iB)
     case(1)
    w(ixO^S,mom(1):mom(3))=0.0d0
    if(mf_glm)  w(ixO^S,psi_)=0.d0

    if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix1=ixOsmax1,ixOsmin1,-1
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
         call mf_face_to_center(ixO^L,block)
       else
         do ix1=ixOmax1,ixOmin1,-1
          ! 2nd order accuacy constant value extrapolation
           w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                      (-w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                 +4.0d0*w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
     endif
     case(2)
    if(mf_glm)  w(ixO^S,psi_)=0.d0
     w(ixO^S,mom(1):mom(3))=0.0d0
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
         call mf_face_to_center(ixO^L,block)
     else
     ! 2nd order accuacy constant value extrapolation
         do ix1=ixOmin1,ixOmax1
           w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                      (-w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                 +4.0d0*w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
    endif
     case(3)
    if(mf_glm)  w(ixO^S,psi_)=0.d0
     w(ixO^S,mom(1):mom(3))=0.0d0
     if(stagger_grid) then
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax^D=ixOmax^D;
           ixOsmin^D=ixOmin^D-kr(^D,idir);
           do ix2=ixOsmax2,ixOsmin2,-1
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
         call mf_face_to_center(ixO^L,block)
      else
    ! 2nd order accuacy constant value extrapolation
         do ix2=ixOmax2,ixOmin2,-1
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                      (-w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,mag(1):mag(3)) &
                 +4.0d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
      endif
     case(4)
    if(mf_glm)  w(ixO^S,psi_)=0.d0
     w(ixO^S,mom(1):mom(3))=0.0d0
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
         call mf_face_to_center(ixO^L,block)
       else
    ! 2nd order accuacy constant value extrapolation
         do ix2=ixOmin2,ixOmax2
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                       (-w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,mag(1):mag(3))&
                  +4.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(1):mag(3)))
         enddo
       endif
     case(5)
    if(mf_glm)  w(ixO^S,psi_)=0.d0
    ! fixed zero bottom velocity
        w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))
        w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(2)) 
        w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(3)) 
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
         call mf_face_to_center(ixO^L,block)
       else     
    ! 4th order accuacy constant value extrapolation
      do ix3=ixOmax3,ixOmin3,-1
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(:))= &
                 0.12d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+5,mag(:)) &
                -0.76d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+4,mag(:)) &
                +2.08d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mag(:)) &
                -3.36d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(:)) &
                +2.92d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(:))
      enddo
     endif 
    !********************************************************
    ! data-driven boundary **********************************
    !********************************************************
    ! Method: V+B driven in linde 
    ! Velocity boundary

     do ixv=1,nx_dd3-1
          if((qt .ge. t_driven_v(ixv)) .and. (qt .lt. t_driven_v(ixv+1))) then
            exit
          end if
     end do
      dtbv=t_driven_v(ixv+1)-t_driven_v(ixv)
      if (ixv .eq. nx_dd3) call mpistop("Special boundary is not defined for this time") 
        do ix3=ixOmax3,ixOmax3    ! only the most inner ghost layer is provided with the driven data
        do ix2=ixMlo2,ixMhi2
        do ix1=ixMlo1,ixMhi1
            xlen1=x(ix1,ix2,ix3,1)-xprobmin1
            xlen2=x(ix1,ix2,ix3,2)-xprobmin2
            ixbc1=ceiling(xlen1/dxb1)
            ixbc2=ceiling(xlen2/dxb2)
            w(ix^D,mom(1):mom(3))=(Vddb(ixbc1,ixbc2,ixv,1:ndir) + &
              (Vddb(ixbc1,ixbc2,ixv+1,1:ndir) - Vddb(ixbc1,ixbc2,ixv,1:ndir))/dtbv* &
              (qt-t_driven_v(ixv)))
        end do
        end do
        end do

    ! Magnetic-field boundary

        do ixb=1,nx_dd3-1
           if((qt .ge. t_driven_b(ixb)) .and. (qt .lt. t_driven_b(ixb+1))) then
             exit
           end if
        end do
        if (ixb .eq. nx_dd3) call mpistop("Special boundary is not defined for this time") 
        dtbv=t_driven_b(ixb+1)-t_driven_b(ixb)

         do ix3=ixOmax3,ixOmax3    ! only the most inner ghost layer is provided with the driven data
         do ix2=ixMlo2,ixMhi2
         do ix1=ixMlo1,ixMhi1
            xlen1=x(ix1,ix2,ix3,1)-xprobmin1
            xlen2=x(ix1,ix2,ix3,2)-xprobmin2
           ixbc1=ceiling(xlen1/dxb1)
            ixbc2=ceiling(xlen2/dxb2)
            w(ix^D,mag(1):mag(3))=Bddb(ixbc1,ixbc2,ixb,1:ndir) + &
              (Bddb(ixbc1,ixbc2,ixb+1,1:ndir) - Bddb(ixbc1,ixbc2,ixb,1:ndir))/dtbv* &
              (qt-t_driven_b(ixb))
        end do
        end do
        end do
         do ix3=ixOmax3-1,ixOmin3,-1
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(:))= &
                 0.12d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+5,mag(:)) &
                -0.76d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+4,mag(:)) &
                +2.08d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mag(:)) &
                -3.36d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(:)) &
                +2.92d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(:))
         enddo

    !********************************************************
    !********************************************************
    case(6)
    if(mf_glm)  w(ixO^S,psi_)=0.d0
    w(ixO^S,mom(1:3))=0.0d0
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
     ! 2nd order accuacy constant value extrapolation
         do ix3=ixOmin3,ixOmax3
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3))=(1.0d0/3.0d0)* &
                 (     -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-2,mag(1):mag(3)) &
                  +4.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-1,mag(1):mag(3)))
         enddo
      endif
     case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  !==============================================================================
  ! Purpose: Enforce additional refinement or coarsening. One can use the
  !          coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !==============================================================================
  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

   ! fix the bottom layer to the highest level
   ! if(ps(igrid)%is_physical_boundary(5)) then
   !   refine=1
   !   coarsen=-1
   ! endif

    if (any(x(ixO^S,3)<=xprobmin3+3.0d0)) then
      refine=1
      coarsen=-1
    endif

   ! if (any(x(ixO^S,3)>xprobmin3+5.0d0)) then
   !   refine=-1
   !   coarsen=1
   ! endif
  end subroutine special_refine_grid

  !==============================================================================
  !> Set the "eta" array for resistive MHD based on w or the
  !> "current" variable which has components between idirmin and 3.
  !==============================================================================
  subroutine specialres_usr(w,ixI^L,ixO^L,idirmin,x,current,eta)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L, idirmin
    double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision             :: current(ixI^S,7-2*ndir:3), eta(ixI^S)
    double precision             :: eta_0, j_mag(ixI^S), ypsl(ixI^S), Bs(ixI^S), mg_dxyz, eta_2, j_c, eta_max
    double precision, dimension(ixI^S,1:ndir) :: Bc
    integer                      :: ix^D

  ! Cheung & Derosa 2012 (only consider the first two terms)   
     eta_0 = 1.72d-5
     eta_max = 5000.0*eta_0
     j_mag = sum((current(ixI^S,:))**2,dim=ndim+1)
     Bc(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
     Bs(ixI^S) = sum((Bc(ixI^S,:))**2,dim=ndim+1)
     ypsl(ixI^S)=j_mag(ixI^S)/Bs(ixI^S)
 
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if (x(ix^D,3) .ge. 0.2d0) then
        mg_dxyz = minval(block%dx(ix^D,1:ndim)) 
        ypsl(ix^D) = ypsl(ix^D)*mg_dxyz**2
        eta(ix^D) = eta_0 + eta_0*10.0d0*ypsl(ix^D)/(1.0d0+dexp(-(ypsl(ix^D)-0.1d0)/0.01d0))
      endif

      if (x(ix^D,3) .le. 0.2d0) then
        eta(ix^D)=eta_0
      endif
    {end do\}

    where(eta(ixO^S)>eta_max)
      eta(ixO^S)=eta_max
    endwhere

  !  Inoue et al 2018  
  !  eta_0 = 1.0d-5
  !  eta_2 = 8.6d-5
  !  j_c    = 35.0d0
  !  j_mag = dsqrt(sum((current(ixI^S,:))**2,dim=ndim+1))
  !  eta_max = 1000.0*eta_2
 
   
  !  {do ix^DB=ixOmin^DB,ixOmax^DB\}
  !    if (j_mag(ix^D) < j_c) then
  !       eta(ix^D) = eta_0 
  !    else
  !     eta(ix^D) = eta_0 + eta_2*((j_mag(ix^D)-j_c)/j_c)**2
  !   endif
  !  {end do\}

  !  where(eta(ixO^S)>eta_max)
  !    eta(ixO^S)=eta_max
  !  endwhere

  end subroutine specialres_usr

  !PURPOSE: use different tolerance in special regions for AMR to
  !reduce/increase resolution there where nothing/something interesting happens.
  subroutine specialthreshold(wlocal,xlocal,tolerance,qt,level)
    use mod_global_parameters
    
    integer, intent(in)          :: level
    double precision, intent(in) :: wlocal(1:nw),xlocal(1:ndim),qt
    double precision, intent(inout) :: tolerance
    
    double precision :: bczone^D,addtol,tol_add

    tol_add=0.d0
    !amplitude of additional tolerance
    addtol=0.3d0
    ! thickness of near-boundary region
    bczone1=1.5d0
    ! linear changing of additional tolerance
    if(xlocal(1)-xprobmin1 < bczone1 .or. xprobmax1-xlocal(1) < bczone1) then
      tol_add=(1.d0-min(xlocal(1)-xprobmin1,xprobmax1-xlocal(1))/bczone1)*addtol
    endif
    bczone2=1.5d0
    if(xlocal(2)-xprobmin2 < bczone2 .or. xprobmax2-xlocal(2) < bczone2) then
      tol_add=(1.d0-min(xlocal(2)-xprobmin2,xprobmax2-xlocal(2))/bczone2)*addtol
    endif
    bczone3=4.d0
    if(xprobmax3-xlocal(3) < bczone3) then
      tol_add=(1.d0-(xprobmax3-xlocal(3))/bczone3)*addtol
    endif
    tolerance=tolerance+tol_add

  end subroutine specialthreshold

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
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: tmp(ixI^S),dip(ixI^S),divb(ixI^S),B2(ixI^S),pth(ixI^S)
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

  end subroutine specialvar_output

  !==============================================================================
  ! Purpose: names for special variable output
  !==============================================================================
  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames

      varnames='divb j1 j2 j3'

  end subroutine specialvarnames_output

end module mod_usr
