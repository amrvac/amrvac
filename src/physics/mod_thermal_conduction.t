!> Thermal conduction for HD and MHD
!>
!> 10.07.2011 developed by Chun Xia and Rony Keppens
!> 01.09.2012 moved to modules folder by Oliver Porth
!> 13.10.2013 optimized further by Chun Xia
!> 12.03.2014 implemented RKL2 super timestepping scheme to reduce iterations 
!> and improve stability and accuracy up to second order in time by Chun Xia.
!> 23.08.2014 implemented saturation and perpendicular TC by Chun Xia
!> 12.01.2017 modulized by Chun Xia
!> 
!> PURPOSE: 
!> IN MHD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(KAPPA_i,j . GRAD_j T)
!> where KAPPA_i,j = tc_k_para b_i b_j + tc_k_perp (I - b_i b_j)
!> b_i b_j = B_i B_j / B**2, I is the unit matrix, and i, j= 1, 2, 3 for 3D
!> IN HD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(tc_k_para . GRAD T)
!> USAGE:
!> 1. in mod_usr.t -> subroutine usr_init(), add 
!>        unit_length=your length unit
!>        unit_numberdensity=your number density unit
!>        unit_velocity=your velocity unit
!>        unit_temperature=your temperature unit
!>    before call (m)hd_activate()
!> 2. to switch on thermal conduction in the (m)hd_list of amrvac.par add:
!>    (m)hd_thermal_conduction=.true.
!> 3. in the tc_list of amrvac.par :
!>    tc_perpendicular=.true.  ! (default .false.) turn on thermal conduction perpendicular to magnetic field 
!>    tc_saturate=.false.  ! (default .true. ) turn off thermal conduction saturate effect
!>    tc_dtpar=0.9/0.45/0.3 ! stable time step coefficient for 1D/2D/3D, decrease it for more stable run
!>    tc_slope_limiter='MC' ! choose limiter for slope-limited anisotropic thermal conduction in MHD

module mod_thermal_conduction
  use mod_global_parameters, only: std_len
  use mod_geometry
  implicit none
  !> Coefficient of thermal conductivity (parallel to magnetic field)
  double precision, public :: tc_k_para

  !> Coefficient of thermal conductivity perpendicular to magnetic field
  double precision, public :: tc_k_perp

  !> Time step of thermal conduction
  double precision :: dt_tc

  !> Number of sub-steps of supertime stepping
  integer, public :: s

  !> Index of the density (in the w array)
  integer, private :: rho_

  !> Index of the energy density (-1 if not present)
  integer, private, protected              :: e_

  !> Index of the internal energy
  integer, private, protected :: eaux_

  !> The adiabatic index
  double precision, private :: tc_gamma

  !> The adiabatic index-1
  double precision, private :: tc_gamma_1

  !> The small_est allowed energy
  double precision, private :: small_e

  !> Time step coefficient
  double precision, private :: tc_dtpar=0.9d0

  !> Maximal substeps of TC within one fluid time step to limit fluid time step
  integer :: tc_ncycles=1000

  !> Calculate thermal conduction perpendicular to magnetic field (.true.) or not (.false.)
  logical, private :: tc_perpendicular=.false.

  !> Consider thermal conduction saturation effect (.true.) or not (.false.)
  logical, private :: tc_saturate=.true.

  !> Logical switch for prepare mpi datatype only once
  logical, private :: first=.true.

  !> Logical switch for test constant conductivity
  logical, private :: tc_constant=.false.

  !> Whether to conserve fluxes at the current partial step
  logical :: fix_conserve_at_step = .true.

  !> Name of slope limiter for transverse component of thermal flux 
  character(len=std_len), private  :: tc_slope_limiter

  procedure(thermal_conduction), pointer   :: phys_thermal_conduction => null()
  procedure(get_heatconduct), pointer   :: phys_get_heatconduct => null()
  procedure(getdt_heatconduct), pointer :: phys_getdt_heatconduct => null()

  abstract interface
    subroutine thermal_conduction
    ! Meyer 2012 MNRAS 422,2102
      use mod_global_parameters
      use mod_ghostcells_update
    end subroutine thermal_conduction

    subroutine get_heatconduct(tmp,ixI^L,ixO^L,w,x,qvec)
      use mod_global_parameters
      
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
      !! tmp store the heat conduction energy changing rate
      double precision, intent(out) :: tmp(ixI^S)
      double precision :: qvec(ixI^S,1:ndim)
    end subroutine get_heatconduct

    subroutine getdt_heatconduct(w,ixG^L,ix^L,dtnew,dx^D,x)
      use mod_global_parameters
      
      integer, intent(in) :: ixG^L, ix^L
      double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
      ! note that depending on small_values_method=='error' etc, w values may change 
      ! through call to getpthermal
      double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
    end subroutine getdt_heatconduct
  end interface

contains
  !> Read this module"s parameters from a file
  subroutine tc_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /tc_list/ tc_perpendicular, tc_saturate, tc_dtpar, tc_slope_limiter, tc_k_para, tc_k_perp, tc_ncycles

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, tc_list, end=111)
111    close(unitpar)
    end do

  end subroutine tc_params_read

  !> Initialize the module
  subroutine thermal_conduction_init(phys_gamma)
    use mod_global_parameters
    use mod_physics

    double precision, intent(in) :: phys_gamma

    tc_gamma=phys_gamma
    tc_gamma_1=phys_gamma-1.d0

    tc_dtpar=tc_dtpar/dble(ndim)

    tc_slope_limiter='MC'

    tc_k_para=0.d0

    tc_k_perp=0.d0

    call tc_params_read(par_files)
  
    phys_thermal_conduction => do_thermal_conduction
    if(physics_type=='hd') then
      phys_get_heatconduct   => hd_get_heatconduct
      phys_getdt_heatconduct => hd_getdt_heatconduct
    else if(physics_type=='mhd') then
      phys_get_heatconduct   => mhd_get_heatconduct
      phys_getdt_heatconduct => mhd_getdt_heatconduct
    end if

    rho_ = iw_rho
    e_ = iw_e
    if(phys_solve_eaux) eaux_ = iw_eaux

    small_e = small_pressure/tc_gamma_1

    if(tc_k_para==0.d0 .and. tc_k_perp==0.d0) then
      if(SI_unit) then
        ! Spitzer thermal conductivity with SI units
        tc_k_para=8.d-12*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 
        ! thermal conductivity perpendicular to magnetic field
        tc_k_perp=4.d-30*unit_numberdensity**2/unit_magneticfield**2/unit_temperature**3*tc_k_para
      else
        ! Spitzer thermal conductivity with cgs units
        tc_k_para=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 
        ! thermal conductivity perpendicular to magnetic field
        tc_k_perp=4.d-10*unit_numberdensity**2/unit_magneticfield**2/unit_temperature**3*tc_k_para
      end if
    else
      tc_constant=.true.
    end if

  end subroutine thermal_conduction_init

  subroutine do_thermal_conduction
  ! Meyer 2012 MNRAS 422,2102
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    use mod_physics
    
    double precision :: omega1,cmu,cmut,cnu,cnut
    double precision, allocatable :: bj(:)
    integer:: iigrid, igrid, j
    logical :: evenstep, stagger_flag, prolong_flag, coarsen_flag

    ixCoGmin^D=1;
    ixCoGmax^D=(ixGhi^D-2*nghostcells)/2+2*nghostcells;
    ! not do fix conserve and getbc for staggered values if stagger is used
    stagger_flag=stagger_grid
    stagger_grid=.false.
    bcphys=.false.
    prolong_flag=prolongprimitive
    coarsen_flag=coarsenprimitive
    prolongprimitive=.false.
    coarsenprimitive=.false.

    ! point bc mpi datatype to partial type for thermalconduction
    type_send_srl=>type_send_srl_p1
    type_recv_srl=>type_recv_srl_p1
    type_send_r=>type_send_r_p1
    type_recv_r=>type_recv_r_p1
    type_send_p=>type_send_p_p1
    type_recv_p=>type_recv_p_p1 
    ! create bc mpi datatype for ghostcells update
    if(first) then
      call create_bc_mpi_datatype(e_,1)
      first=.false.
    end if

    call init_comm_fix_conserve(1,ndim,1)
    fix_conserve_at_step = time_advance .and. levmax>levmin

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ! convert total energy to internal energy
      call phys_e_to_ei(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
      ! internal e of the coarse block maybe needed in physical boundaries
      if(any(ps(igrid)%is_physical_boundary)) &
        call phys_e_to_ei(ixCoG^L,ixCoG^L,psc(igrid)%w,psc(igrid)%x)
      if(.not. allocated(ps2(igrid)%w)) allocate(ps2(igrid)%w(ixG^T,1:nw))
      if(.not. allocated(ps3(igrid)%w)) allocate(ps3(igrid)%w(ixG^T,1:nw))
      ps1(igrid)%w=ps(igrid)%w
      ps2(igrid)%w=ps(igrid)%w
      ps3(igrid)%w=ps(igrid)%w
    end do

    allocate(bj(0:s))
    bj(0)=1.d0/3.d0
    bj(1)=bj(0)
    if(s>1) then
      omega1=4.d0/dble(s**2+s-2)
      cmut=omega1/3.d0
    else
      omega1=0.d0
      cmut=1.d0
    endif
    
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      block=>ps(igrid)
      typelimiter=type_limiter(node(plevel_,igrid))
      typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      call evolve_step1(igrid,cmut,dt_tc,ixG^LL,ixM^LL,ps1(igrid)%w,ps(igrid)%w,&
                        ps(igrid)%x,ps3(igrid)%w)
    end do
    !$OMP END PARALLEL DO
    ! fix conservation of AMR grid by replacing flux from finer neighbors
    if (fix_conserve_at_step) then
      call recvflux(1,ndim)
      call sendflux(1,ndim)
      call fix_conserve(ps1,1,ndim,e_,1)
    end if
    call getbc(global_time,0.d0,ps1,e_,1)
    if(s==1) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ps(igrid)%w(ixG^T,e_)=ps1(igrid)%w(ixG^T,e_)
        if(phys_solve_eaux) ps(igrid)%w(ixG^T,eaux_)=ps(igrid)%w(ixG^T,e_)
        ! convert internal energy to total energy
        call phys_ei_to_e(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
        if(any(ps(igrid)%is_physical_boundary)) &
          call phys_ei_to_e(ixCoG^L,ixCoG^L,psc(igrid)%w,psc(igrid)%x)
      end do
      ! point bc mpi data type back to full type for (M)HD
      type_send_srl=>type_send_srl_f
      type_recv_srl=>type_recv_srl_f
      type_send_r=>type_send_r_f
      type_recv_r=>type_recv_r_f
      type_send_p=>type_send_p_f
      type_recv_p=>type_recv_p_f
      bcphys=.true.
      stagger_grid=stagger_flag
      prolongprimitive=prolong_flag
      coarsenprimitive=coarsen_flag
      deallocate(bj)
      return
    endif
    evenstep=.true.
    do j=2,s
      bj(j)=dble(j**2+j-2)/dble(2*j*(j+1))
      cmu=dble(2*j-1)/dble(j)*bj(j)/bj(j-1)
      cmut=omega1*cmu
      cnu=dble(1-j)/dble(j)*bj(j)/bj(j-2)
      cnut=(bj(j-1)-1.d0)*cmut
      if(evenstep) then
    !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          block=>ps(igrid)
          typelimiter=type_limiter(node(plevel_,igrid))
          typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          call evolve_stepj(igrid,cmu,cmut,cnu,cnut,dt_tc,ixG^LL,ixM^LL,ps1(igrid)%w,&
                            ps2(igrid)%w,ps(igrid)%w,ps(igrid)%x,ps3(igrid)%w)
        end do
    !$OMP END PARALLEL DO
        ! fix conservation of AMR grid by replacing flux from finer neighbors
        if (fix_conserve_at_step) then
          call recvflux(1,ndim)
          call sendflux(1,ndim)
          call fix_conserve(ps2,1,ndim,e_,1)
        end if
        call getbc(global_time,0.d0,ps2,e_,1)
        evenstep=.false.
      else
    !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          block=>ps(igrid)
          typelimiter=type_limiter(node(plevel_,igrid))
          typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          call evolve_stepj(igrid,cmu,cmut,cnu,cnut,dt_tc,ixG^LL,ixM^LL,ps2(igrid)%w,&
                            ps1(igrid)%w,ps(igrid)%w,ps(igrid)%x,ps3(igrid)%w)
        end do
    !$OMP END PARALLEL DO
        ! fix conservation of AMR grid by replacing flux from finer neighbors
        if (fix_conserve_at_step) then
          call recvflux(1,ndim)
          call sendflux(1,ndim)
          call fix_conserve(ps1,1,ndim,e_,1)
        end if
        call getbc(global_time,0.d0,ps1,e_,1)
        evenstep=.true.
      end if 
    end do
    if(evenstep) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ps(igrid)%w(ixG^T,e_)=ps1(igrid)%w(ixG^T,e_)
        if(phys_solve_eaux) ps(igrid)%w(ixG^T,eaux_)=ps(igrid)%w(ixG^T,e_)
        ! convert internal energy to total energy
        call phys_ei_to_e(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
        if(any(ps(igrid)%is_physical_boundary)) &
          call phys_ei_to_e(ixCoG^L,ixCoG^L,psc(igrid)%w,psc(igrid)%x)
      end do 
    else
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ps(igrid)%w(ixG^T,e_)=ps2(igrid)%w(ixG^T,e_)
        if(phys_solve_eaux) ps(igrid)%w(ixG^T,eaux_)=ps(igrid)%w(ixG^T,e_)
        ! convert internal energy to total energy
        call phys_ei_to_e(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
        if(any(ps(igrid)%is_physical_boundary)) &
          call phys_ei_to_e(ixCoG^L,ixCoG^L,psc(igrid)%w,psc(igrid)%x)
      end do 
    end if
    deallocate(bj)
    ! point bc mpi data type back to full type for (M)HD
    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f
    bcphys=.true.

    ! restore stagger_grid value
    stagger_grid=stagger_flag
    prolongprimitive=prolong_flag
    coarsenprimitive=coarsen_flag

  end subroutine do_thermal_conduction

  subroutine addsource_impl
    use mod_global_parameters
    use mod_ghostcells_update

    integer :: iigrid, igrid, icycle, ncycle
    double precision :: qt

    ncycle=ceiling(0.5d0*dt/dt_tc)
    if(ncycle<1) then
      ncycle=1
      dt_tc=0.5d0*dt
    else
      dt_tc=0.5d0*dt/dble(ncycle)
    endif

    if(mype==0.and..false.) then
      print *,'implicit source addition will subcycle with ',ncycle,' subtimesteps'
      print *,'dt and dtimpl= ',dt,dt_tc,' versus ncycle*dtimpl=',ncycle*dt_tc
    endif

    qt=global_time
    do icycle=1,ncycle
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
        block=>ps(igrid)
        ps1(igrid)%w=ps(igrid)%w
        call evolve_step1(igrid,1.d0,dt_tc,ixG^LL,ixM^LL,ps1(igrid)%w,ps(igrid)%w,&
                          ps(igrid)%x,ps3(igrid)%w)
      end do
      !$OMP END PARALLEL DO
      qt=qt+dt_tc
      call getbc(qt,0.d0,ps,e_,1)
    end do

  end subroutine addsource_impl

  subroutine evolve_stepj(igrid,qcmu,qcmut,qcnu,qcnut,qdt,ixI^L,ixO^L,w1,w2,w,x,w3)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: igrid,ixI^L,ixO^L
    double precision, intent(in) :: qcmu,qcmut,qcnu,qcnut,qdt
    double precision, intent(in) :: w1(ixI^S,1:nw),w(ixI^S,1:nw),w3(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w2(ixI^S,1:nw)
    
    double precision :: fC(ixI^S,1:1,1:ndim)
    double precision :: tmp(ixI^S)

    call phys_get_heatconduct(tmp,ixI^L,ixO^L,w1,x,fC)

    w2(ixO^S,e_)=qcmu*w1(ixO^S,e_)+qcnu*w2(ixO^S,e_)+(1.d0-qcmu-qcnu)*w(ixO^S,e_)&
                +qcmut*qdt*tmp(ixO^S)+qcnut*w3(ixO^S,e_)

    ! check small/negative internal energy
    if(check_small_values) call handle_small_e(w2,x,ixI^L,ixO^L,'thermal conduction evolve_stepj')

    if (fix_conserve_at_step) then
      fC=qcmut*qdt*fC
      call store_flux(igrid,fC,1,ndim,1)
    end if

  end subroutine evolve_stepj

  subroutine evolve_step1(igrid,qcmut,qdt,ixI^L,ixO^L,w1,w,x,w3)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: igrid,ixI^L,ixO^L
    double precision, intent(in) :: qcmut, qdt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) ::w1(ixI^S,1:nw),w3(ixI^S,1:nw)
    
    double precision :: fC(ixI^S,1:1,1:ndim)
    double precision :: tmp(ixI^S)
    integer :: lowindex(ndim), ix^D

    call phys_get_heatconduct(tmp,ixI^L,ixO^L,w,x,fC)

    w3(ixO^S,e_)=qdt*tmp(ixO^S)
    ! update internal energy
    w1(ixO^S,e_) = w(ixO^S,e_) + qcmut*w3(ixO^S,e_)
    
    ! check small/negative internal energy
    call handle_small_e(w1,x,ixI^L,ixO^L,'thermal conduction evolve_step1')

    if (fix_conserve_at_step) then
      fC=qcmut*qdt*fC
      call store_flux(igrid,fC,1,ndim,1)
    end if
  
  end subroutine evolve_step1

  subroutine handle_small_e(w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixI^S,1:nw)

    flag=.false.
    where(w(ixO^S,e_)<small_e) flag(ixO^S,e_)=.true.
    if(any(flag(ixO^S,e_))) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S,e_)) w(ixO^S,e_)=small_e
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag, e_)
      case default
        ! small values error shows primitive variables
        w(ixO^S,e_)=w(ixO^S,e_)*tc_gamma_1
        do idir = 1, ndir
           w(ixO^S, iw_mom(idir)) = w(ixO^S, iw_mom(idir))/w(ixO^S,rho_)
        end do
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
  end subroutine handle_small_e

  !> anisotropic thermal conduction with slope limited symmetric scheme
  !> Sharma 2007 Journal of Computational Physics 227, 123
  subroutine mhd_get_heatconduct(qd,ixI^L,ixO^L,w,x,qvec)
    use mod_global_parameters
    use mod_small_values, only: small_values_method
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
    !! qd store the heat conduction energy changing rate
    double precision, intent(out) :: qd(ixI^S)
    double precision, dimension(ixI^S,1:ndim) :: qvec
    
    double precision, dimension(ixI^S,1:ndir) :: mf,Bc,Bcf
    double precision, dimension(ixI^S,1:ndim) :: gradT
    double precision, dimension(ixI^S) :: Te,ka,kaf,ke,kef,qdd,qe,Binv,minq,maxq,Bnorm
    double precision :: alpha,dxinv(ndim)
    integer, dimension(ndim) :: lowindex
    integer :: idims,idir,ix^D,ix^L,ixC^L,ixA^L,ixB^L

    ! coefficient of limiting on normal component
    if(ndim<3) then
      alpha=0.75d0
    else
      alpha=0.85d0
    end if
    ix^L=ixO^L^LADD1;

    dxinv=1.d0/dxlevel

    ! compute the temperature
    Te(ixI^S)=w(ixI^S,e_)*tc_gamma_1/w(ixI^S,rho_)
    ! B vector
    if(B0field) then
      mf(ixI^S,:)=w(ixI^S,iw_mag(:))+block%B0(ixI^S,:,0)
    else
      mf(ixI^S,:)=w(ixI^S,iw_mag(:));
    end if
    ! |B|
    Binv(ix^S)=dsqrt(sum(mf(ix^S,:)**2,dim=ndim+1))
    where(Binv(ix^S)/=0.d0)
      Binv(ix^S)=1.d0/Binv(ix^S)
    elsewhere
      Binv(ix^S)=bigdouble
    end where
    ! b unit vector: magnetic field direction vector
    do idims=1,ndim
      mf(ix^S,idims)=mf(ix^S,idims)*Binv(ix^S)
    end do
    ! ixC is cell-corner index
    ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;
    ! b unit vector at cell corner
    Bc=0.d0
    {do ix^DB=0,1\}
      ixAmin^D=ixCmin^D+ix^D;
      ixAmax^D=ixCmax^D+ix^D;
      Bc(ixC^S,1:ndim)=Bc(ixC^S,1:ndim)+mf(ixA^S,1:ndim)
    {end do\}
    Bc(ixC^S,1:ndim)=Bc(ixC^S,1:ndim)*0.5d0**ndim
    ! T gradient at cell faces
    gradT=0.d0
    do idims=1,ndim
      ixBmin^D=ixmin^D;
      ixBmax^D=ixmax^D-kr(idims,^D);
      call gradientC(Te,ixI^L,ixB^L,idims,minq)
      gradT(ixB^S,idims)=minq(ixB^S)
    end do
    if(tc_constant) then
      if(tc_perpendicular) then
        ka(ixC^S)=tc_k_para-tc_k_perp
        ke(ixC^S)=tc_k_perp
      else
        ka(ixC^S)=tc_k_para
      end if
    else
      ! conductivity at cell center
      if(trac) then
        minq(ix^S)=Te(ix^S)
        where(minq(ix^S) < block%special_values(1))
          minq(ix^S)=block%special_values(1)
        end where
        minq(ix^S)=tc_k_para*sqrt(minq(ix^S)**5)
      else
        minq(ix^S)=tc_k_para*sqrt(Te(ix^S)**5)
      end if
      ka=0.d0
      {do ix^DB=0,1\}
        ixBmin^D=ixCmin^D+ix^D;
        ixBmax^D=ixCmax^D+ix^D;
        ka(ixC^S)=ka(ixC^S)+minq(ixB^S)
      {end do\}
      ! cell corner conductivity
      ka(ixC^S)=0.5d0**ndim*ka(ixC^S)
      ! compensate with perpendicular conductivity
      if(tc_perpendicular) then
        minq(ix^S)=tc_k_perp*w(ix^S,rho_)**2*Binv(ix^S)**2/dsqrt(Te(ix^S))
        ke=0.d0
        {do ix^DB=0,1\}
          ixBmin^D=ixCmin^D+ix^D;
          ixBmax^D=ixCmax^D+ix^D;
          ke(ixC^S)=ke(ixC^S)+minq(ixB^S)
        {end do\}
        ! cell corner conductivity: k_parallel-k_perpendicular
        ke(ixC^S)=0.5d0**ndim*ke(ixC^S)
        where(ke(ixC^S)<ka(ixC^S))
          ka(ixC^S)=ka(ixC^S)-ke(ixC^S)
        elsewhere
          ka(ixC^S)=0.d0
          ke(ixC^S)=ka(ixC^S)
        end where
      end if
    end if
    if(tc_slope_limiter=='no') then
      ! calculate thermal conduction flux with symmetric scheme
      do idims=1,ndim
        qd=0.d0
        {do ix^DB=0,1 \}
           if({ ix^D==0 .and. ^D==idims | .or.}) then
             ixBmin^D=ixCmin^D+ix^D;
             ixBmax^D=ixCmax^D+ix^D;
             qd(ixC^S)=qd(ixC^S)+gradT(ixB^S,idims)
           end if
        {end do\}
        ! temperature gradient at cell corner
        qvec(ixC^S,idims)=qd(ixC^S)*0.5d0**(ndim-1)
      end do
      ! b grad T at cell corner
      qd(ixC^S)=sum(qvec(ixC^S,1:ndim)*Bc(ixC^S,1:ndim),dim=ndim+1)
      do idims=1,ndim
        ! TC flux at cell corner
        gradT(ixC^S,idims)=ka(ixC^S)*Bc(ixC^S,idims)*qd(ixC^S)
        if(tc_perpendicular) gradT(ixC^S,idims)=gradT(ixC^S,idims)+ke(ixC^S)*qvec(ixC^S,idims)
      end do
      ! TC flux at cell face
      qvec=0.d0
      do idims=1,ndim
        ixB^L=ixO^L-kr(idims,^D);
        ixAmax^D=ixOmax^D; ixAmin^D=ixBmin^D;
        {do ix^DB=0,1 \}
           if({ ix^D==0 .and. ^D==idims | .or.}) then
             ixBmin^D=ixAmin^D-ix^D; 
             ixBmax^D=ixAmax^D-ix^D; 
             qvec(ixA^S,idims)=qvec(ixA^S,idims)+gradT(ixB^S,idims)
           end if
        {end do\}
        qvec(ixA^S,idims)=qvec(ixA^S,idims)*0.5d0**(ndim-1)
        if(tc_saturate) then
          ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135)
          ! unsigned saturated TC flux = 5 phi rho c**3, c is isothermal sound speed
          Bcf=0.d0
          {do ix^DB=0,1 \}
             if({ ix^D==0 .and. ^D==idims | .or.}) then
               ixBmin^D=ixAmin^D-ix^D;
               ixBmax^D=ixAmax^D-ix^D;
               Bcf(ixA^S,idims)=Bcf(ixA^S,idims)+Bc(ixB^S,idims)
             end if
          {end do\}
          ! averaged b at face centers
          Bcf(ixA^S,idims)=Bcf(ixA^S,idims)*0.5d0**(ndim-1)
          ixB^L=ixA^L+kr(idims,^D);
          qd(ixA^S)=2.75d0*(w(ixA^S,rho_)+w(ixB^S,rho_))*dsqrt(0.5d0*(Te(ixA^S)+Te(ixB^S)))**3*dabs(Bcf(ixA^S,idims))
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            if(dabs(qvec(ix^D,idims))>qd(ix^D)) then
              qvec(ix^D,idims)=sign(1.d0,qvec(ix^D,idims))*qd(ix^D)
            end if
         {end do\}
        end if
      end do
    else
      ! calculate thermal conduction flux with slope-limited symmetric scheme
      qvec=0.d0
      do idims=1,ndim
        ixB^L=ixO^L-kr(idims,^D);
        ixAmax^D=ixOmax^D; ixAmin^D=ixBmin^D;
        ! calculate normal of magnetic field
        ixB^L=ixA^L+kr(idims,^D);
        Bnorm(ixA^S)=0.5d0*(mf(ixA^S,idims)+mf(ixB^S,idims))
        Bcf=0.d0
        kaf=0.d0
        kef=0.d0
        {do ix^DB=0,1 \}
           if({ ix^D==0 .and. ^D==idims | .or.}) then
             ixBmin^D=ixAmin^D-ix^D;
             ixBmax^D=ixAmax^D-ix^D;
             Bcf(ixA^S,1:ndim)=Bcf(ixA^S,1:ndim)+Bc(ixB^S,1:ndim)
             kaf(ixA^S)=kaf(ixA^S)+ka(ixB^S)
             if(tc_perpendicular) kef(ixA^S)=kef(ixA^S)+ke(ixB^S)
           end if
        {end do\}
        ! averaged b at face centers
        Bcf(ixA^S,1:ndim)=Bcf(ixA^S,1:ndim)*0.5d0**(ndim-1)
        ! averaged thermal conductivity at face centers
        kaf(ixA^S)=kaf(ixA^S)*0.5d0**(ndim-1)
        if(tc_perpendicular) kef(ixA^S)=kef(ixA^S)*0.5d0**(ndim-1)
        ! limited normal component
        minq(ixA^S)=min(alpha*gradT(ixA^S,idims),gradT(ixA^S,idims)/alpha)
        maxq(ixA^S)=max(alpha*gradT(ixA^S,idims),gradT(ixA^S,idims)/alpha)
        ! eq (19)
        qdd=0.d0
        {do ix^DB=0,1 \}
           if({ ix^D==0 .and. ^D==idims | .or.}) then
             ixBmin^D=ixCmin^D+ix^D;
             ixBmax^D=ixCmax^D+ix^D;
             qdd(ixC^S)=qdd(ixC^S)+gradT(ixB^S,idims)
           end if
        {end do\}
        ! temperature gradient at cell corner
        qdd(ixC^S)=qdd(ixC^S)*0.5d0**(ndim-1)
        ! eq (21)
        qe=0.d0
        {do ix^DB=0,1 \}
           qd(ixC^S)=qdd(ixC^S)
           if({ ix^D==0 .and. ^D==idims | .or.}) then
             ixBmin^D=ixAmin^D-ix^D;
             ixBmax^D=ixAmax^D-ix^D;
             where(qd(ixB^S)<=minq(ixA^S))
               qd(ixB^S)=minq(ixA^S)
             elsewhere(qd(ixB^S)>=maxq(ixA^S))
               qd(ixB^S)=maxq(ixA^S)
             end where
             qvec(ixA^S,idims)=qvec(ixA^S,idims)+Bc(ixB^S,idims)**2*qd(ixB^S)
             if(tc_perpendicular) qe(ixA^S)=qe(ixA^S)+qd(ixB^S) 
           end if
        {end do\}
        qvec(ixA^S,idims)=kaf(ixA^S)*qvec(ixA^S,idims)*0.5d0**(ndim-1)
        ! add normal flux from perpendicular conduction
        if(tc_perpendicular) qvec(ixA^S,idims)=qvec(ixA^S,idims)+kef(ixA^S)*qe(ixA^S)*0.5d0**(ndim-1)
        ! limited transverse component, eq (17)
        ixBmin^D=ixAmin^D;
        ixBmax^D=ixAmax^D+kr(idims,^D);
        do idir=1,ndim
          if(idir==idims) cycle
          qd(ixI^S)=slope_limiter(gradT(ixI^S,idir),ixI^L,ixB^L,idir,-1)
          qd(ixI^S)=slope_limiter(qd,ixI^L,ixA^L,idims,1)
          qvec(ixA^S,idims)=qvec(ixA^S,idims)+kaf(ixA^S)*Bnorm(ixA^S)*Bcf(ixA^S,idir)*qd(ixA^S)
        end do

        ! consider magnetic null point
        !where(Binv(ixA^S)==0.d0)
        !  qvec(ixA^S,idims)=tc_k_para*(0.5d0*(Te(ixA^S)+Te(ixB^S)))**2.5d0*gradT(ixA^S,idims)
        !end where

        if(tc_saturate) then
          ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135)
          ! unsigned saturated TC flux = 5 phi rho c**3, c is isothermal sound speed
          ixB^L=ixA^L+kr(idims,^D);
          qd(ixA^S)=2.75d0*(w(ixA^S,rho_)+w(ixB^S,rho_))*dsqrt(0.5d0*(Te(ixA^S)+Te(ixB^S)))**3*dabs(Bnorm(ixA^S))
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            if(dabs(qvec(ix^D,idims))>qd(ix^D)) then
        !      write(*,*) 'it',it,qvec(ix^D,idims),qd(ix^D),' TC saturated at ',&
        !      x(ix^D,:),' rho',w(ix^D,rho_),' Te',Te(ix^D)
              qvec(ix^D,idims)=sign(1.d0,qvec(ix^D,idims))*qd(ix^D)
            end if
         {end do\}
        end if
      end do
    end if

    qd=0.d0
    if(slab_uniform) then
      do idims=1,ndim
        qvec(ix^S,idims)=dxinv(idims)*qvec(ix^S,idims)
        ixB^L=ixO^L-kr(idims,^D);
        qd(ixO^S)=qd(ixO^S)+qvec(ixO^S,idims)-qvec(ixB^S,idims)
      end do
    else
      do idims=1,ndim
        qvec(ix^S,idims)=qvec(ix^S,idims)*block%surfaceC(ix^S,idims)
        ixB^L=ixO^L-kr(idims,^D);
        qd(ixO^S)=qd(ixO^S)+qvec(ixO^S,idims)-qvec(ixB^S,idims)
      end do
      qd(ixO^S)=qd(ixO^S)/block%dvolume(ixO^S)
    end if
    
  end subroutine mhd_get_heatconduct

  !> Calculate gradient of a scalar q at cell interfaces in direction idir
  subroutine gradientC(q,ixI^L,ixO^L,idir,gradq)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: q(ixI^S)
    double precision, intent(inout) :: gradq(ixI^S)
    integer                         :: jxO^L

    associate(x=>block%x)

    jxO^L=ixO^L+kr(idir,^D);
    select case(coordinate)
    case(Cartesian)
      gradq(ixO^S)=(q(jxO^S)-q(ixO^S))/dxlevel(idir)
    case(Cartesian_stretched)
      gradq(ixO^S)=(q(jxO^S)-q(ixO^S))/(x(jxO^S,idir)-x(ixO^S,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixO^S)=(q(jxO^S)-q(ixO^S))/(x(jxO^S,1)-x(ixO^S,1))
        {^NOONED
      case(2)
        gradq(ixO^S)=(q(jxO^S)-q(ixO^S))/( (x(jxO^S,2)-x(ixO^S,2))*x(ixO^S,1) )
        }
        {^IFTHREED
      case(3)
        gradq(ixO^S)=(q(jxO^S)-q(ixO^S))/( (x(jxO^S,3)-x(ixO^S,3))*x(ixO^S,1)*dsin(x(ixO^S,2)) )
        }
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixO^S)=(q(jxO^S)-q(ixO^S))/((x(jxO^S,phi_)-x(ixO^S,phi_))*x(ixO^S,r_))
      else
        gradq(ixO^S)=(q(jxO^S)-q(ixO^S))/(x(jxO^S,idir)-x(ixO^S,idir))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

    end associate
  end subroutine gradientC

  function slope_limiter(f,ixI^L,ixO^L,idims,pm) result(lf)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L, idims, pm
    double precision, intent(in) :: f(ixI^S)
    double precision :: lf(ixI^S)

    double precision :: signf(ixI^S)
    integer :: ixB^L

    ixB^L=ixO^L+pm*kr(idims,^D);
    signf(ixO^S)=sign(1.d0,f(ixO^S))
    select case(tc_slope_limiter)
     case('minmod')
       ! minmod limiter
       lf(ixO^S)=signf(ixO^S)*max(0.d0,min(abs(f(ixO^S)),signf(ixO^S)*f(ixB^S)))
     case ('MC')
       ! montonized central limiter Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       lf(ixO^S)=two*signf(ixO^S)* &
            max(zero,min(dabs(f(ixO^S)),signf(ixO^S)*f(ixB^S),&
            signf(ixO^S)*quarter*(f(ixB^S)+f(ixO^S))))
     case ('superbee')
       ! Roes superbee limiter (eq.3.51i)
       lf(ixO^S)=signf(ixO^S)* &
            max(zero,min(two*dabs(f(ixO^S)),signf(ixO^S)*f(ixB^S)),&
            min(dabs(f(ixO^S)),two*signf(ixO^S)*f(ixB^S)))
     case ('koren')
       ! Barry Koren Right variant
       lf(ixO^S)=signf(ixO^S)* &
            max(zero,min(two*dabs(f(ixO^S)),two*signf(ixO^S)*f(ixB^S),&
            (two*f(ixB^S)*signf(ixO^S)+dabs(f(ixO^S)))*third))
     case default
       call mpistop("Unknown slope limiter for thermal conduction")
    end select

  end function slope_limiter

  subroutine mhd_getdt_heatconduct(w,ixI^L,ixO^L,dtnew,dx^D,x)
    !Check diffusion time limit dt < tc_dtpar*dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_physics
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    ! note that depending on small_values_method=='error' etc, w values may change 
    ! through call to getpthermal
    double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
    
    double precision :: dxinv(1:ndim),mf(ixI^S,1:ndir)
    double precision :: tmp2(ixI^S),tmp(ixI^S),Te(ixI^S),B2(ixI^S)
    double precision :: dtdiff_tcond, dtdiff_tsat
    integer          :: idim,ix^D

    ^D&dxinv(^D)=one/dx^D;
    
    call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)

    !temperature
    Te(ixO^S)=tmp(ixO^S)/w(ixO^S,rho_)
    !tc_k_para_i
    if(tc_constant) then
      tmp(ixO^S)=tc_k_para
    else
      tmp(ixO^S)=tc_k_para*dsqrt(Te(ixO^S)**5)/w(ixO^S,rho_)
    end if
    
    ! B
    if(B0field) then
      mf(ixO^S,:)=w(ixO^S,iw_mag(:))+block%B0(ixO^S,:,0)
    else
      mf(ixO^S,:)=w(ixO^S,iw_mag(:))
    end if
    ! B^-2
    B2(ixO^S)=sum(mf(ixO^S,:)**2,dim=ndim+1)
    ! B_i**2/B**2
    where(B2(ixO^S)/=0.d0)
      ^D&mf(ixO^S,^D)=mf(ixO^S,^D)**2/B2(ixO^S);
    elsewhere
      ^D&mf(ixO^S,^D)=1.d0;
    end where

    if(tc_saturate) B2(ixO^S)=22.d0*dsqrt(Te(ixO^S))

    do idim=1,ndim
      tmp2(ixO^S)=tmp(ixO^S)*mf(ixO^S,idim)
      if(tc_saturate) then
        where(tmp2(ixO^S)>B2(ixO^S))
          tmp2(ixO^S)=B2(ixO^S)
        end where
      end if
      ! dt< tc_dtpar * dx_idim**2/((gamma-1)*tc_k_para_i/rho*B_i**2/B**2)
      dtdiff_tcond=tc_dtpar/(tc_gamma-1.d0)/maxval(tmp2(ixO^S)*dxinv(idim)**2)
      ! limit the time step
      dtnew=min(dtnew,dtdiff_tcond)
    end do
    dtnew=dtnew/dble(ndim)
  
  end subroutine mhd_getdt_heatconduct

  subroutine hd_get_heatconduct(qd,ixI^L,ixO^L,w,x,qvec)
    use mod_global_parameters
    use mod_small_values, only: small_values_method
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
    !! tmp store the heat conduction energy changing rate
    double precision, intent(out) :: qd(ixI^S)
    double precision :: qvec(ixI^S,1:ndim),gradT(ixI^S,1:ndim),Te(ixI^S),ke(ixI^S)
    double precision :: dxinv(ndim)
    integer, dimension(ndim)       :: lowindex
    integer :: idims,ix^D,ix^L,ixC^L,ixA^L,ixB^L,ixD^L

    ix^L=ixO^L^LADD1;
    ! ixC is cell-corner index
    ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;

    dxinv=1.d0/dxlevel

    ! compute temperature before source addition
    Te(ixI^S)=w(ixI^S,e_)*tc_gamma_1/w(ixI^S,rho_)

    ! cell corner temperature
    ke=0.d0
    ixAmax^D=ixmax^D; ixAmin^D=ixmin^D-1;
    {do ix^DB=0,1\}
      ixBmin^D=ixAmin^D+ix^D;
      ixBmax^D=ixAmax^D+ix^D;
      ke(ixA^S)=ke(ixA^S)+Te(ixB^S)
    {end do\}
    ke(ixA^S)=0.5d0**ndim*ke(ixA^S)
    ! T gradient (central difference) at cell corners
    gradT=0.d0
    do idims=1,ndim
      ixBmin^D=ixmin^D;
      ixBmax^D=ixmax^D-kr(idims,^D);
      call gradient(ke,ixI^L,ixB^L,idims,qd)
      gradT(ixB^S,idims)=qd(ixB^S)
    end do
    ! transition region adaptive conduction
    if(trac) then
      where(ke(ix^S) < block%special_values(1))
        ke(ix^S)=block%special_values(1)
      end where
    end if
    ! cell corner conduction flux
    do idims=1,ndim
      gradT(ixC^S,idims)=gradT(ixC^S,idims)*tc_k_para*sqrt(ke(ixC^S)**5)
    end do

    if(tc_saturate) then
      ! consider saturation with unsigned saturated TC flux = 5 phi rho c**3
      ! saturation flux at cell center
      qd(ix^S)=5.d0*w(ix^S,rho_)*dsqrt(Te(ix^S)**3)
      ke=0.d0
      {do ix^DB=0,1\}
        ixBmin^D=ixCmin^D+ix^D;
        ixBmax^D=ixCmax^D+ix^D;
        ke(ixC^S)=ke(ixC^S)+qd(ixB^S)
      {end do\}
      ! cell corner saturation flux 
      ke(ixC^S)=0.5d0**ndim*ke(ixC^S)
      ! magnitude of cell corner conduction flux
      qd(ixC^S)=norm2(gradT(ixC^S,:),dim=ndim+1)
      {do ix^DB=ixCmin^DB,ixCmax^DB\}
        if(qd(ix^D)>ke(ix^D)) then
          ke(ix^D)=ke(ix^D)/qd(ix^D)
          do idims=1,ndim
            gradT(ix^D,idims)=ke(ix^D)*gradT(ix^D,idims)
          end do
        end if
      {end do\}
    end if

    ! conductionflux at cell face
    qvec=0.d0
    do idims=1,ndim
      ixB^L=ixO^L-kr(idims,^D);
      ixAmax^D=ixOmax^D; ixAmin^D=ixBmin^D;
      {do ix^DB=0,1 \}
         if({ ix^D==0 .and. ^D==idims | .or.}) then
           ixBmin^D=ixAmin^D-ix^D; 
           ixBmax^D=ixAmax^D-ix^D; 
           qvec(ixA^S,idims)=qvec(ixA^S,idims)+gradT(ixB^S,idims)
         end if
      {end do\}
      qvec(ixA^S,idims)=qvec(ixA^S,idims)*0.5d0**(ndim-1)
    end do

    qd=0.d0
    if(slab_uniform) then
      do idims=1,ndim
        qvec(ix^S,idims)=dxinv(idims)*qvec(ix^S,idims)
        ixB^L=ixO^L-kr(idims,^D);
        qd(ixO^S)=qd(ixO^S)+qvec(ixO^S,idims)-qvec(ixB^S,idims)
      end do
    else
      do idims=1,ndim
        qvec(ix^S,idims)=qvec(ix^S,idims)*block%surfaceC(ix^S,idims)
        ixB^L=ixO^L-kr(idims,^D);
        qd(ixO^S)=qd(ixO^S)+qvec(ixO^S,idims)-qvec(ixB^S,idims)
      end do
      qd(ixO^S)=qd(ixO^S)/block%dvolume(ixO^S)
    end if

  end subroutine hd_get_heatconduct

  subroutine hd_getdt_heatconduct(w,ixI^L,ixO^L,dtnew,dx^D,x)
    ! Check diffusion time limit dt < tc_dtpar * dx_i**2 / ((gamma-1)*tc_k_para_i/rho)
    use mod_global_parameters
    use mod_physics

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
    
    double precision :: dxinv(1:ndim), tmp(ixI^S), Te(ixI^S)
    double precision :: dtdiff_tcond,dtdiff_tsat
    integer          :: idim,ix^D

    ^D&dxinv(^D)=one/dx^D;

    call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)

    Te(ixO^S)=tmp(ixO^S)/w(ixO^S,rho_)
    tmp(ixO^S)=(tc_gamma-one)*tc_k_para*dsqrt((Te(ixO^S))**5)/w(ixO^S,rho_)
    
    do idim=1,ndim
       ! dt< tc_dtpar * dx_idim**2/((gamma-1)*tc_k_para_idim/rho)
       dtdiff_tcond=tc_dtpar/maxval(tmp(ixO^S)*dxinv(idim)**2)
       if(tc_saturate) then
         ! dt< tc_dtpar* dx_idim**2/((gamma-1)*sqrt(Te)*5*phi)
         dtdiff_tsat=tc_dtpar/maxval((tc_gamma-1.d0)*dsqrt(Te(ixO^S))*&
                     5.d0*dxinv(idim)**2)
         ! choose the slower flux (bigger time scale) between classic and saturated
         dtdiff_tcond=max(dtdiff_tcond,dtdiff_tsat)
       end if
       ! limit the time step
       dtnew=min(dtnew,dtdiff_tcond)
    end do
    dtnew=dtnew/dble(ndim)
  
  end subroutine hd_getdt_heatconduct

end module mod_thermal_conduction
