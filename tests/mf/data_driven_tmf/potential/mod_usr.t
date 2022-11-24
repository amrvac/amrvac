module mod_usr   ! computing potential field 
  use mod_mhd
  use mod_lfff
  implicit none
  ! some global parameters
  double precision :: k_B,miu0,mass_H,usr_grav,SRadius,rhob,Tiso,lalpha,llift
  logical, save :: firstusrglobaldata=.true.

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_3D")

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_global_parameters

    ! normalization unit in CGS Unit
    k_B = 1.3806d-16          ! erg*K^-1
    miu0 = 4.d0*dpi           ! Gauss^2 cm^2 dyne^-1
    mass_H = 1.67262d-24      ! g
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3
    unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    unit_time          = unit_length/unit_velocity                     ! 85.8746159942810 s 

    usr_grav=-2.74d4*unit_length/unit_velocity**2  ! solar gravity
    SRadius=6.955d10/unit_length                   ! Solar radius
    rhob=1.d0                                      ! bottom density
    Tiso=1.d6/unit_temperature                     ! uniform temperature

    lalpha=0.d0    ! alpha coefficient for linear force-free field
    llift=0.0736762d0*1.0d0*0.5d0     ! lift up domain above bottom magnetogram
                   ! Usually, we need to lift up domain above bottom magnetogram by 0.5*dz, since the potential
                   ! field assumes the bottom boundary lies on the cell surface between the
                   ! ghost layers and the physcial domain, while we want the bottom boundary
                   ! lies on the cell center of the inner ghost layer.
 
    ! prepare magnetogram at bottom
    if(firstusrglobaldata) then
      call init_b_fff_data('potential_boundary/potential_boundary.dat',unit_length,unit_magneticfield)
      firstusrglobaldata=.false.
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
    ! use green function method to extrapolate B
    call calc_lin_fff(ixI^L,ixO^L,Bf,x,lalpha,llift)
    w(ixO^S,mag(:))=Bf(ixO^S,1:3)
    w(ixO^S,mom(:))=zero
    w(ixO^S,rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*&
                  (1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
   ! w(ixO^S,psi_)=0.d0

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters
    use mod_physics
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: ft,tfstop,tramp1,tramp2,coeffrho,vlimit,vsign
    double precision :: delxdely,delxdelz,delydelx,delydelz,delzdelx,delzdely
    double precision :: xlen^D,dxa^D,startpos^D
    integer :: ix^D,ixIM^L,ixbc^D,af

    select case(iB)
    case(1)

   !  w(ixO^S,psi_)=0.d0

    case(2)
   !  w(ixO^S,psi_)=0.d0


    case(3)
   !   w(ixO^S,psi_)=0.d0


    case(4)
   !  w(ixO^S,psi_)=0.d0
    case(5)
   !  w(ixO^S,psi_)=0.d0

    case(6)
   !  w(ixO^S,psi_)=0.d0

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

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

   ! fix the bottom layer to the highest level
     ! refine all blocks to highest level to match the grid of initial magnetic
     ! field
    !    refine=1
    !    coarsen=-1

   if (qt == 0) then
    if (any(x(ixO^S,3)<=xprobmin3+1.0d0)) then
      refine=1
      coarsen=-1
    endif

    if (any(x(ixO^S,3)>xprobmin3+1.0d0)) then
      refine=-1
      coarsen=1
    endif
  endif
  end subroutine special_refine_grid

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
    call divvector(Btotal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+2)=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+3)=w(ixO^S,rho_)*Tiso*two/B2(ixO^S)
    ! store current
    call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
    do idir=1,ndir
      w(ixO^S,nw+3+idir)=curlvec(ixO^S,idir)
    end do
    ! find magnetic dips
!    dip=0.d0
!    do idir=1,ndir
!      call gradient(w(ixI^S,mag(3)),ixI^L,ixO^L,idir,tmp)
!      dip(ixO^S)=dip(ixO^S)+w(ixO^S,b0_+idir)*tmp(ixO^S)
!    end do
!    where(dabs(w(ixO^S,mag(3)))<0.08d0 .and. dip(ixO^S)>=0.d0)
!      w(ixO^S,nw+8)=1.d0
!    elsewhere
!      w(ixO^S,nw+8)=0.d0
!    end where
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Alfv divB beta j1 j2 j3'
  end subroutine specialvarnames_output

end module mod_usr





