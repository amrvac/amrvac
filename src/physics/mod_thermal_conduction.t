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
!> where KAPPA_i,j = kappa b_i b_j + kappe (I - b_i b_j)
!> b_i b_j = B_i B_j / B**2, I is the unit matrix, and i, j= 1, 2, 3 for 3D
!> IN HD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(kappa . GRAD T)
!> USAGE:
!> 1. in mod_usr.t -> subroutine usr_init(), add 
!>        unit_length=<your length unit>
!>        unit_numberdensity=<your number density unit>
!>        unit_velocity=<your velocity unit>
!>        unit_temperature=<your temperature unit>
!>    before call (m)hd_activate()
!> 2. to switch on thermal conduction in the (m)hd_list of amrvac.par add:
!>    (m)hd_thermal_conduction=.true.
!> 3. in the tc_list of amrvac.par :
!>    tc_saturate=.false.  ! (default .true. ) turn off thermal conduction saturate effect
!>    tc_dtpar=0.9/0.45/0.3 ! stable time step coefficient for 1D/2D/3D, decrease it for more stable run

module mod_thermal_conduction
  use mod_global_parameters, only: std_len
  implicit none
  !> Coefficient of thermal conductivity
  double precision, public :: kappa

  !> Coefficient of thermal conductivity perpendicular to magnetic field
  double precision, public :: kappe

  !> Time step of thermal conduction
  double precision :: dt_tc

  !> Number of sub-steps of supertime stepping
  integer, public :: s

  !> Index of the density (in the w array)
  integer, private :: rho_

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, private, protected              :: e_

  !> Indices of the magnetic field
  integer, allocatable, private, protected :: mag(:)

  !> The adiabatic index
  double precision, private :: tc_gamma

  !> The small_est allowed energy
  double precision, private :: small_e

  !> Time step coefficient
  double precision, private :: tc_dtpar=0.9d0

  !> Maximal number of sub-cycles within one fuild time step
  integer :: ncyclemax=1000

  !> Calculate thermal conduction perpendicular to magnetic field (.true.) or not (.false.)
  logical, private :: tc_perpendicular=.false.

  !> Consider thermal conduction saturation effect (.true.) or not (.false.)
  logical, private :: tc_saturate=.true.

  !> Logical switch for prepare mpi datatype only once
  logical, private :: first=.true.

  !> Logical switch for test constant conductivity
  logical, private :: tc_constant=.false.

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

    subroutine get_heatconduct(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
      use mod_global_parameters
      
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
      !! tmp store the heat conduction energy changing rate
      double precision, intent(out) :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
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

    namelist /tc_list/ tc_perpendicular, tc_saturate, tc_dtpar, tc_slope_limiter, tc_constant

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

    tc_dtpar=tc_dtpar/dble(ndim)

    tc_slope_limiter='MC'

    call tc_params_read(par_files)
  
    phys_thermal_conduction => do_thermal_conduction
    if(physics_type=='hd') then
      phys_get_heatconduct   => hd_get_heatconduct
      phys_getdt_heatconduct => hd_getdt_heatconduct
    else if(physics_type=='mhd') then
      phys_get_heatconduct   => mhd_get_heatconduct
      phys_getdt_heatconduct => mhd_getdt_heatconduct
    end if

    allocate(mom(ndir), mag(ndir))
    rho_ = iw_rho
    e_ = iw_e
    mom(:) = iw_mom(:)
    mag(:) = iw_mag(:)

    small_e = small_pressure/(tc_gamma - 1.0d0)

    if(SI_unit) then
      ! Spitzer thermal conductivity with SI units
      kappa=8.d-12*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 
      ! thermal conductivity perpendicular to magnetic field
      kappe=4.d-30*unit_numberdensity**2/unit_magneticfield**2/unit_temperature**3*kappa
    else
      ! Spitzer thermal conductivity with cgs units
      kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 
      ! thermal conductivity perpendicular to magnetic field
      kappe=4.d-26*unit_numberdensity**2/unit_magneticfield**2/unit_temperature**3*kappa
    end if

    if(tc_constant) kappa=0.01d0

  end subroutine thermal_conduction_init

  subroutine do_thermal_conduction
  ! Meyer 2012 MNRAS 422,2102
    use mod_global_parameters
    use mod_ghostcells_update
    
    double precision :: omega1,cmu,cmut,cnu,cnut
    double precision, allocatable :: bj(:)
    integer:: iigrid, igrid,j
    logical :: evenstep

    ! point bc mpi datatype to partial type for thermalconduction
    type_send_srl=>type_send_srl_p
    type_recv_srl=>type_recv_srl_p
    type_send_r=>type_send_r_p
    type_recv_r=>type_recv_r_p
    type_send_p=>type_send_p_p
    type_recv_p=>type_recv_p_p 
    ! create bc mpi datatype for ghostcells update
    if(first) then
      call create_bc_mpi_datatype(e_-1,1)
      first=.false.
    end if
    
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      if(.not. allocated(pw(igrid)%w2)) allocate(pw(igrid)%w2(ixG^T,1:nw))
      if(.not. allocated(pw(igrid)%w3)) allocate(pw(igrid)%w3(ixG^T,1:nw))
      pw(igrid)%w1=pw(igrid)%w
      pw(igrid)%w2=pw(igrid)%w
      pw(igrid)%w3=pw(igrid)%w
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
      block=>pw(igrid)
      typelimiter=limiter(node(plevel_,igrid))
      typegradlimiter=gradient_limiter(node(plevel_,igrid))
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      call evolve_step1(cmut,dt_tc,ixG^LL,ixM^LL,pw(igrid)%w1,pw(igrid)%w,&
                        pw(igrid)%x,pw(igrid)%w3)
      pw(igrid)%wb=>pw(igrid)%w1
    end do
    !$OMP END PARALLEL DO
    bcphys=.false.
    call getbc(global_time,0.d0,e_-1,1)
    if(s==1) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        pw(igrid)%w(ixG^T,e_)=pw(igrid)%w1(ixG^T,e_)
        pw(igrid)%wb=>pw(igrid)%w
      end do
      ! point bc mpi data type back to full type for (M)HD
      type_send_srl=>type_send_srl_f
      type_recv_srl=>type_recv_srl_f
      type_send_r=>type_send_r_f
      type_recv_r=>type_recv_r_f
      type_send_p=>type_send_p_f
      type_recv_p=>type_recv_p_f
      bcphys=.true.
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
          block=>pw(igrid)
          typelimiter=limiter(node(plevel_,igrid))
          typegradlimiter=gradient_limiter(node(plevel_,igrid))
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          call evolve_stepj(cmu,cmut,cnu,cnut,dt_tc,ixG^LL,ixM^LL,pw(igrid)%w1,&
                            pw(igrid)%w2,pw(igrid)%w,pw(igrid)%x,pw(igrid)%w3)
          pw(igrid)%wb=>pw(igrid)%w2
        end do
    !$OMP END PARALLEL DO
        call getbc(global_time,0.d0,e_-1,1)
        evenstep=.false.
      else
    !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          block=>pw(igrid)
          typelimiter=limiter(node(plevel_,igrid))
          typegradlimiter=gradient_limiter(node(plevel_,igrid))
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          call evolve_stepj(cmu,cmut,cnu,cnut,dt_tc,ixG^LL,ixM^LL,pw(igrid)%w2,&
                            pw(igrid)%w1,pw(igrid)%w,pw(igrid)%x,pw(igrid)%w3)
          pw(igrid)%wb=>pw(igrid)%w1
        end do
    !$OMP END PARALLEL DO
        call getbc(global_time,0.d0,e_-1,1)
        evenstep=.true.
      end if 
    end do
    if(evenstep) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        pw(igrid)%w(ixG^T,e_)=pw(igrid)%w1(ixG^T,e_)
        pw(igrid)%wb=>pw(igrid)%w
      end do 
    else
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        pw(igrid)%w(ixG^T,e_)=pw(igrid)%w2(ixG^T,e_)
        pw(igrid)%wb=>pw(igrid)%w
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
  
  end subroutine do_thermal_conduction

  subroutine evolve_stepj(qcmu,qcmut,qcnu,qcnut,qdt,ixI^L,ixO^L,w1,w2,w,x,wold)
    
    use mod_global_parameters
    
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qcmu,qcmut,qcnu,qcnut,qdt
    double precision, intent(in) :: w1(ixI^S,1:nw),w(ixI^S,1:nw),wold(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w2(ixI^S,1:nw)
    
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)

    call phys_get_heatconduct(tmp,tmp1,tmp2,ixI^L,ixO^L,w1,x)

    w2(ixO^S,e_)=qcmu*w1(ixO^S,e_)+qcnu*w2(ixO^S,e_)+(1.d0-qcmu-qcnu)*w(ixO^S,e_)&
                +qcmut*qdt*tmp(ixO^S)+qcnut*wold(ixO^S,e_)
    
  end subroutine evolve_stepj

  subroutine evolve_step1(qcmut,qdt,ixI^L,ixO^L,w1,w,x,wold)
    use mod_global_parameters
    use mod_small_values, only: small_values_method
    
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qcmut, qdt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) ::w1(ixI^S,1:nw),wold(ixI^S,1:nw)
    
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),Te(ixI^S)
    integer :: lowindex(ndim), ix^D

    call phys_get_heatconduct(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)

    wold(ixO^S,e_)=qdt*tmp(ixO^S)
    ! update internal energy
    tmp1(ixO^S) = tmp1(ixO^S) + qcmut*wold(ixO^S,e_)
    
    ! ensure you never trigger negative pressure 
    ! hence code up energy change with respect to kinetic and magnetic
    ! part(nonthermal)
    if(small_values_method=='error') then
      if(any(tmp1(ixO^S)<small_e).and. .not.crash) then
        lowindex=minloc(tmp1(ixO^S))
        ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
        write(*,*)'too small internal energy = ',minval(tmp1(ixO^S)),'at x=',&
       x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',small_e,&
       ' on time=',global_time,' step=',it, 'where w(1:nwflux)=',w(^D&lowindex(^D),1:nwflux)
        crash=.true.
      else
        w1(ixO^S,e_) = tmp2(ixO^S)+tmp1(ixO^S)
      end if
    else
      where(tmp1(ixO^S)<small_e)
        tmp1(ixO^S)=small_e
      endwhere
      w1(ixO^S,e_) = tmp2(ixO^S)+tmp1(ixO^S)
    end if
  
  end subroutine evolve_step1

  !> anisotropic thermal conduction with slope limited symmetric scheme
  !> Sharma 2007 Journal of Computational Physics 227, 123
  subroutine mhd_get_heatconduct(qd,tmp1,tmp2,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_small_values, only: small_values_method
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
    !! qd store the heat conduction energy changing rate
    double precision, intent(out) :: qd(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
    
    double precision, dimension(ixI^S,1:ndir) :: mf,Bc,Bcf
    double precision, dimension(ixI^S,1:ndim) :: qvec,gradT
    double precision, dimension(ixI^S) :: Te,Tc,Tcf,Binv,cs3,qsatflux,minq,maxq
    double precision :: alpha
    integer, dimension(ndim) :: lowindex
    integer :: idims,idir,ix^D,ix^L,ixC^L,ixA^L,ixB^L
    logical :: Bnull(ixI^S)

    alpha=0.75d0
    ix^L=ixO^L^LADD1;

    ! tmp2 store kinetic+magnetic energy before addition of heat conduction source
    tmp2(ixI^S) = 0.5d0 * (sum(w(ixI^S,mom(:))**2,dim=ndim+1)/w(ixI^S,rho_) + &
         sum(w(ixI^S,mag(:))**2,dim=ndim+1))

    ! tmp1 store internal energy
    tmp1(ixI^S)=w(ixI^S,e_)-tmp2(ixI^S)

    ! Clip off negative pressure if small_pressure is set
    if(small_values_method=='error') then
       if (any(tmp1(ixI^S)<small_e) .and. .not.crash) then
         lowindex=minloc(tmp1(ixI^S))
         ^D&lowindex(^D)=lowindex(^D)+ixImin^D-1;
         write(*,*)'too low internal energy = ',minval(tmp1(ixI^S)),' at x=',&
         x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',small_e,' on time=',global_time, ' it=',it
         write(*,*) 'w',w(^D&lowindex(^D),:)
         crash=.true.
       end if
    else
    {do ix^DB=ixImin^DB,ixImax^DB\}
       if(tmp1(ix^D)<small_e) then
          tmp1(ix^D)=small_e
       end if
    {end do\}
    end if
    ! compute the temperature
    Te(ixI^S)=tmp1(ixI^S)*(tc_gamma-one)/w(ixI^S,rho_)
    ! B vector
    if(B0field) then
      mf(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,0)
    else
      mf(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir));
    end if
    ! |B|
    Binv(ix^S)=dsqrt(sum(mf(ix^S,:)**2,dim=ndim+1))
    Bnull(ixI^S)=.false.
    where(Binv(ix^S)/=0.d0)
      Binv(ix^S)=1.d0/Binv(ix^S)
    elsewhere
      Binv(ix^S)=0.d0
      Bnull(ix^S)=.true.
    end where
    ! b unit vector
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
      Bc(ixC^S,:)=Bc(ixC^S,:)+mf(ixA^S,:)
    {end do\}
    Bc(ixC^S,:)=Bc(ixC^S,:)*0.5d0**ndim
    if(tc_constant) then
      Tc(ixC^S)=kappa
    else
      ! Temperature at cell corner
      Tc=0.d0
      {do ix^DB=0,1\}
        ixAmin^D=ixCmin^D+ix^D;
        ixAmax^D=ixCmax^D+ix^D;
        Tc(ixC^S)=Tc(ixC^S)+Te(ixA^S)
      {end do\}
      ! thermal conductivity kappa*T^2.5
      Tc(ixC^S)=kappa*(0.5d0**ndim*Tc(ixC^S))**2.5d0
    end if
    ! T gradient
    gradT=0.d0
    do idims=1,ndim
      ixBmin^D=ixmin^D;
      ixBmax^D=ixmax^D-kr(idims,^D);
      ixA^L=ixB^L+kr(idims,^D);
      gradT(ixB^S,idims)=(Te(ixA^S)-Te(ixB^S))/dxlevel(idims)
    end do
    ! calculate thermal conduction flux for all dimensions
    qvec=0.d0
    do idims=1,ndim
      ixB^L=ixO^L-kr(idims,^D);
      ixAmax^D=ixOmax^D; ixAmin^D=ixBmin^D;
      Bcf=0.d0
      Tcf=0.d0
      {do ix^DB=0,1 \}
         if({ ix^D==0 .and. ^D==idims | .or.}) then
           ixBmin^D=ixAmin^D-ix^D;
           ixBmax^D=ixAmax^D-ix^D;
           Bcf(ixA^S,:)=Bcf(ixA^S,:)+Bc(ixB^S,:)
           Tcf(ixA^S)=Tcf(ixA^S)+Tc(ixB^S)
         end if
      {end do\}
      ! averaged b at face centers
      Bcf(ixA^S,:)=Bcf(ixA^S,:)*0.5d0**(ndim-1)
      ! averaged thermal conductivity at face centers
      Tcf(ixA^S)=Tcf(ixA^S)*0.5d0**(ndim-1)
      ! limited normal component
      minq(ixA^S)=min(alpha*gradT(ixA^S,idims),gradT(ixA^S,idims)/alpha)
      maxq(ixA^S)=max(alpha*gradT(ixA^S,idims),gradT(ixA^S,idims)/alpha)
      ! eq (19)
      qd=0.d0
      {do ix^DB=0,1 \}
         if({ ix^D==0 .and. ^D==idims | .or.}) then
           ixBmin^D=ixCmin^D+ix^D;
           ixBmax^D=ixCmax^D+ix^D;
           qd(ixC^S)=qd(ixC^S)+gradT(ixB^S,idims)
         end if
      {end do\}
      qd(ixC^S)=qd(ixC^S)*0.5d0**(ndim-1)
      ! eq (21)
      {do ix^DB=0,1 \}
         if({ ix^D==0 .and. ^D==idims | .or.}) then
           ixBmin^D=ixAmin^D-ix^D;
           ixBmax^D=ixAmax^D-ix^D;
           where(qd(ixB^S)<=minq(ixA^S))
             qd(ixB^S)=minq(ixA^S)
           elsewhere(qd(ixB^S)>=maxq(ixA^S))
             qd(ixB^S)=maxq(ixA^S)
           end where
           qvec(ixA^S,idims)=qvec(ixA^S,idims)+Bc(ixB^S,idims)**2*qd(ixB^S)
         end if
      {end do\}
      qvec(ixA^S,idims)=Tcf(ixA^S)*qvec(ixA^S,idims)*0.5d0**(ndim-1)
      ! limited transverse component, eq (17)
      ixBmin^D=ixAmin^D;
      ixBmax^D=ixAmax^D+kr(idims,^D);
      do idir=1,ndim
        if(idir==idims) cycle
        qd(ixI^S)=slope_limiter(gradT(ixI^S,idir),ixI^L,ixB^L,idir,-1)
        qd(ixI^S)=slope_limiter(qd,ixI^L,ixA^L,idims,1)
        qvec(ixA^S,idims)=qvec(ixA^S,idims)+Tcf(ixA^S)*Bcf(ixA^S,idims)*Bcf(ixA^S,idir)*qd(ixA^S)
      end do
    end do

    if(tc_saturate) then
      ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135)
      cs3(ix^S)=dsqrt(Te(ix^S)**3)
      ! unsigned saturated TC flux = 5 phi rho c**3
      qsatflux(ix^S)=5.d0*w(ix^S,rho_)*cs3(ix^S)
      ! strength of classic TC flux
      qd(ix^S)=dsqrt(^D&qvec(ix^S,^D)**2+)
      ! sign(b * Grad Te) 5 phi rho c**3 Bi/B 
      {do ix^DB=ixmin^DB,ixmax^DB\}
        if(.false. .and. qd(ix^D)>qsatflux(ix^D) .and. idims==1) write(*,*)& 
          'it',it,' ratio=',qd(ix^D)/qsatflux(ix^D),' TC saturated at ',&
          x(ix^D,:),' rho',w(ix^D,rho_),' Te',Te(ix^D)
        if(qd(ix^D)>qsatflux(ix^D)) then
        ! saturated TC flux = sign(b * Grad Te) 5 phi rho c**3
          qsatflux(ix^D)=sign(1.d0,sum(mf(ix^D,:)*gradT(ix^D,:)))*qsatflux(ix^D)*Binv(ix^D)
          do idims=1,ndim
            qvec(ix^D,idims)=qsatflux(ix^D)*mf(ix^D,idims)
          end do
        end if
      {end do\}
    end if
    
    {do ix^DB=ixmin^DB,ixmax^DB\}
      if(Bnull(ix^D)) then
        ! consider magnetic null point
        qvec(ix^D,1:ndim)=kappa*dsqrt(Te(ix^D)**5)*gradT(ix^D,1:ndim)
        if(tc_saturate) then
          ! unsigned saturated TC flux = 5 phi rho c**3
          qsatflux(ix^D)=5.d0*w(ix^D,rho_)*cs3(ix^D)
          qd(ix^D)=dsqrt(sum(qvec(ix^D,1:ndim)**2))
          if(qd(ix^D)>qsatflux(ix^D)) then
            qsatflux(ix^D)=qsatflux(ix^D)/qd(ix^D)
            do idims=1,ndim
              qvec(ix^D,idims)=qsatflux(ix^D)*qvec(ix^D,idims)
            end do
          end if
        end if
      end if  
    {end do\}

    qd=0.d0
    do idims=1,ndim
      ixB^L=ixO^L-kr(idims,^D);
      qd(ixO^S)=qd(ixO^S)+(qvec(ixO^S,idims)-qvec(ixB^S,idims))/dxlevel(idims)
    end do
    
  end subroutine mhd_get_heatconduct

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
     case ('mcbeta')
       ! Woodward and Collela limiter, with factor beta
       lf(ixO^S)=signf(ixO^S)* &
            max(zero,min(mcbeta*dabs(f(ixO^S)),mcbeta*signf(ixO^S)*f(ixB^S),&
            signf(ixO^S)*half*(f(ixB^S)+f(ixO^S))))
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
    !Check diffusion time limit dt < tc_dtpar*dx_i**2/((gamma-1)*kappa_i/rho)
    !where                      kappa_i=kappa*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_physics
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    ! note that depending on small_values_method=='error' etc, w values may change 
    ! through call to getpthermal
    double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
    
    double precision :: dxinv(1:ndim),mf(ixI^S,1:ndir)
    double precision :: tmp2(ixI^S),tmp(ixI^S),Te(ixI^S),B2inv(ixI^S)
    double precision :: dtdiff_tcond, dtdiff_tsat
    integer          :: idim,idir,ix^D

    ^D&dxinv(^D)=one/dx^D;
    
    call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)

    !temperature
    Te(ixO^S)=tmp(ixO^S)/w(ixO^S,rho_)
    !kappa_i
    if(tc_constant) then
      tmp(ixO^S)=kappa
    else
      tmp(ixO^S)=kappa*dsqrt(Te(ixO^S)**5)
    end if
    !(gamma-1)*kappa_i/rho
    tmp(ixO^S)=(tc_gamma-one)*tmp(ixO^S)/w(ixO^S,rho_)
    
    ! B
    if(B0field) then
      mf(ixO^S,1:ndir)=w(ixO^S,mag(1):mag(ndir))+block%B0(ixO^S,1:ndir,0)
    else
      mf(ixO^S,1:ndir)=w(ixO^S,mag(1):mag(ndir))
    end if
    ! B^-2
    B2inv=0.d0
    do idir=1,ndir
      B2inv(ixO^S)=B2inv(ixO^S)+mf(ixO^S,idir)**2
    end do
    where(B2inv(ixO^S)/=0.d0)
      B2inv(ixO^S)=1.d0/B2inv(ixO^S)
    end where

    do idim=1,ndim
       ! B_i**2/B**2
       where(B2inv(ixO^S)/=0.d0)
         tmp2(ixO^S)=mf(ixO^S,idim)**2*B2inv(ixO^S)
       elsewhere
         tmp2(ixO^S)=1.d0
       end where
       ! dt< tc_dtpar * dx_idim**2/((gamma-1)*kappa_i/rho*B_i**2/B**2)
       dtdiff_tcond=tc_dtpar/maxval(tmp(ixO^S)*tmp2(ixO^S)*dxinv(idim)**2)
       if(tc_saturate) then
         ! dt< tc_dtpar* dx_idim**2/((gamma-1)*sqrt(Te)*5*phi)
         ! with an empirical coefficient dx_idim
         dtdiff_tsat=tc_dtpar/maxval((tc_gamma-1.d0)*dsqrt(Te(ixO^S))*&
                     5.d0*dxinv(idim)**2)
         ! choose the slower flux (bigger time step) between classic and saturated
         dtdiff_tcond=max(dtdiff_tcond,dtdiff_tsat)
       end if
       ! limit the time step
       dtnew=min(dtnew,dtdiff_tcond)
    end do
    dtnew=dtnew/dble(ndim)
  
  end subroutine mhd_getdt_heatconduct

  subroutine hd_get_heatconduct(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_small_values, only: small_values_method
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
    !! tmp store the heat conduction energy changing rate
    double precision, intent(out) :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
    double precision :: qvec(ixI^S,1:ndir),Te(ixI^S),qflux(ixI^S),qsatflux(ixI^S)
    integer:: ix^L,idims,ix^D
    integer, dimension(ndim)       :: lowindex
    
    ix^L=ixO^L^LADD2;
    if (ixI^L^LTix^L|.or.|.or.) &
       call mpistop("Need two extra layers for thermal conduction")
    
    ! store old kinetic energy
    tmp2(ixI^S)=half*sum(w(ixI^S,mom(:))**2,dim=ndim+1)/w(ixI^S,rho_)
    ! store old internal energy
    tmp1(ixI^S)=w(ixI^S,e_)-tmp2(ixI^S)
    ! Clip off negative pressure if small_pressure is set
    if(small_values_method=='error' .and. .not.crash) then
      if(any(tmp1(ixI^S)<small_e)) then
        lowindex=minloc(tmp1(ixI^S))
        write(*,*)'low internal energy= ',minval(tmp1(ixI^S)),' at x=',&
        x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',small_e,' on time=',global_time
        crash=.true.
      end if
    else
       {do ix^DB=ixImin^DB,ixImax^DB\}
         if(tmp1(ix^D)<small_e) then
          tmp1(ix^D)=small_e
         end if
       {end do\}
    end if

    ! compute temperature before source addition
    Te(ixI^S)=tmp1(ixI^S)/w(ixI^S,rho_)*(tc_gamma-one)

    ! compute grad T and store grad T vector
    do idims=1,ndim
      ! idirth component of gradient of temperature at cell center
       select case(typegrad)
       case("central")
         ix^L=ixI^L^LSUB1;
         call gradient(Te,ixI^L,ix^L,idims,tmp)
       case("limited")
         ix^L=ixI^L^LSUB2;
         call gradientS(Te,ixI^L,ix^L,idims,tmp)
       end select
       qvec(ix^S,idims)=tmp(ix^S)*kappa*dsqrt(Te(ix^S)**5)
    end do

    if(tc_saturate) then
      ! consider saturation with unsigned saturated TC flux = 5 phi rho c**3
      qsatflux(ix^S)=5.d0*w(ix^S,rho_)*dsqrt(Te(ix^S)**3)
      qflux(ix^S)=dsqrt(^D&qvec(ix^S,^D)**2+)
      {do ix^DB=ixmin^DB,ixmax^DB\}
        if(qflux(ix^D)>qsatflux(ix^D)) then
          qsatflux(ix^D)=qsatflux(ix^D)/qflux(ix^D)
          do idims=1,ndim
            qvec(ix^D,idims)=qsatflux(ix^D)*qvec(ix^D,idims)
          end do
        end if
      {end do\}
    end if

    select case(typediv)
      case("central")
        call divvector(qvec,ixI^L,ixO^L,tmp)
      case("limited")
        call divvectorS(qvec,ixI^L,ixO^L,tmp)
    end select

  end subroutine hd_get_heatconduct

  subroutine hd_getdt_heatconduct(w,ixI^L,ixO^L,dtnew,dx^D,x)
    ! Check diffusion time limit dt < tc_dtpar * dx_i**2 / ((gamma-1)*kappa_i/rho)
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
    tmp(ixO^S)=(tc_gamma-one)*kappa*dsqrt((Te(ixO^S))**5)/w(ixO^S,rho_)
    
    do idim=1,ndim
       ! dt< tc_dtpar * dx_idim**2/((gamma-1)*kappa_idim/rho)
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
