!> Thermal conduction for HD and MHD or RHD and RMHD or twofl (plasma-neutral) module
!> Adaptation of mod_thermal_conduction for the mod_supertimestepping
!>
!> The TC is set by calling 
!> tc_init_params()
!>
!> Organized such that it can call either isotropic (HD) or anisotropic (MHD) variants
!> it adds a heat conduction source to each energy equation
!> and can be recycled within a multi-fluid context (such as plasma-neutral twofl module)
!>
!>
!> 10.07.2011 developed by Chun Xia and Rony Keppens
!> 01.09.2012 moved to modules folder by Oliver Porth
!> 13.10.2013 optimized further by Chun Xia
!> 12.03.2014 implemented RKL2 super timestepping scheme to reduce iterations
!> and improve stability and accuracy up to second order in time by Chun Xia.
!> 23.08.2014 implemented saturation and perpendicular TC by Chun Xia
!> 12.01.2017 modulized by Chun Xia
!>            adapted by Beatrice Popescu to twofluid settings
!> 06.09.2024 cleaned up for use in rhd and rmhd modules (Nishant Narechania and Rony Keppens)
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
!> 2. to switch on thermal conduction in the (r)(m)hd_list of amrvac.par add:
!>    (r)(m)hd_thermal_conduction=.true.
!> 3. in the tc_list of amrvac.par :
!>    tc_perpendicular=.true.  ! (default .false.) turn on thermal conduction perpendicular to magnetic field
!>    tc_saturate=.true.  ! (default .false. ) turn on thermal conduction saturate effect
!>    tc_slope_limiter='MC' ! choose limiter for slope-limited anisotropic thermal conduction in MHD
!> note: twofl_list incorporates instances for charges and neutrals

module mod_thermal_conduction
  use mod_global_parameters, only: std_len
  use mod_geometry
  use mod_comm_lib, only: mpistop
  implicit none

    !> The adiabatic index
    double precision :: tc_gamma_1

  abstract interface
    subroutine get_var_subr(w,x,ixI^L,ixO^L,res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: res(ixI^S)
    end subroutine get_var_subr
  end interface

  type tc_fluid

    ! BEGIN the following are read from param file or set in tc_read_hd_params or tc_read_mhd_params
    !> Coefficient of thermal conductivity (parallel to magnetic field)
    double precision :: tc_k_para

    !> Coefficient of thermal conductivity perpendicular to magnetic field
    double precision :: tc_k_perp

     !> Indices of the variables
    integer :: e_=-1
    !> Index of cut off temperature for TRAC
    integer :: Tcoff_
    !> Name of slope limiter for transverse component of thermal flux
    integer :: tc_slope_limiter

    ! if has_equi = .true. get_temperature_equi and get_rho_equi have to be set
    logical :: has_equi=.false.

    !> Logical switch for test constant conductivity
    logical :: tc_constant=.false.

    !> Calculate thermal conduction perpendicular to magnetic field (.true.) or not (.false.)
    logical :: tc_perpendicular=.false.

    !> Consider thermal conduction saturation effect (.true.) or not (.false.)
    logical :: tc_saturate=.false.
    ! END the following are read from param file or set in tc_read_hd_params or tc_read_mhd_params
    procedure (get_var_subr), pointer, nopass :: get_rho => null()
    procedure (get_var_subr), pointer, nopass :: get_rho_equi => null()
    procedure(get_var_subr), pointer,nopass :: get_temperature_from_eint => null()
    procedure(get_var_subr), pointer,nopass :: get_temperature_from_conserved => null()
    procedure(get_var_subr), pointer,nopass :: get_temperature_equi => null()
  end type tc_fluid


  public :: tc_get_mhd_params
  public :: tc_get_hd_params
  public :: tc_get_ffhd_params
  public :: get_tc_dt_mhd
  public :: get_tc_dt_hd
  public :: get_tc_dt_ffhd
  public :: sts_set_source_tc_mhd
  public :: sts_set_source_tc_hd
  public :: sts_set_source_tc_ffhd

contains

  subroutine tc_init_params(phys_gamma)
    use mod_global_parameters
    double precision, intent(in) :: phys_gamma

    tc_gamma_1=phys_gamma-1d0
  end subroutine tc_init_params

  !> Init  TC coefficients: MHD case
  subroutine tc_get_mhd_params(fl,read_mhd_params)
    use mod_global_parameters

    interface
      subroutine read_mhd_params(fl)
        use mod_global_parameters, only: unitpar,par_files
        import tc_fluid
        type(tc_fluid), intent(inout) :: fl

      end subroutine read_mhd_params
    end interface
    type(tc_fluid), intent(inout) :: fl

    fl%tc_slope_limiter=1
    fl%tc_k_para=0.d0
    fl%tc_k_perp=0.d0

    !> Read tc module parameters from par file: MHD case
    call read_mhd_params(fl)

    if(fl%tc_k_para==0.d0 .and. fl%tc_k_perp==0.d0) then
      if(SI_unit) then
        ! Spitzer thermal conductivity with SI units
        fl%tc_k_para=8.d-12*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
        ! thermal conductivity perpendicular to magnetic field
        fl%tc_k_perp=4.d-30*unit_numberdensity**2/unit_magneticfield**2/unit_temperature**3*fl%tc_k_para
      else
        ! Spitzer thermal conductivity with cgs units
        fl%tc_k_para=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
        ! thermal conductivity perpendicular to magnetic field
        fl%tc_k_perp=4.d-10*unit_numberdensity**2/unit_magneticfield**2/unit_temperature**3*fl%tc_k_para
      end if
      if(mype .eq. 0) print*, "Spitzer MHD: par: ",fl%tc_k_para, &
          " ,perp: ",fl%tc_k_perp
    else
      fl%tc_constant=.true.
    end if

  end subroutine tc_get_mhd_params

  !> Init  TC coefficients: HD case
  subroutine tc_get_hd_params(fl,read_hd_params)
    use mod_global_parameters

    interface
      subroutine read_hd_params(fl)
        use mod_global_parameters, only: unitpar,par_files
        import tc_fluid
        type(tc_fluid), intent(inout) :: fl

      end subroutine read_hd_params
    end interface
    type(tc_fluid), intent(inout) :: fl

    fl%tc_k_para=0.d0

    !> Read tc parameters from par file: HD case
    call read_hd_params(fl)

    if(fl%tc_k_para==0.d0 ) then
      if(SI_unit) then
        ! Spitzer thermal conductivity with SI units
        fl%tc_k_para=8.d-12*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      else
        ! Spitzer thermal conductivity with cgs units
        fl%tc_k_para=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      end if
      if(mype .eq. 0) print*, "Spitzer HD par: ",fl%tc_k_para
    end if

  end subroutine tc_get_hd_params

  subroutine tc_get_ffhd_params(fl,read_ffhd_params)
    use mod_global_parameters

    interface
      subroutine read_ffhd_params(fl)
        use mod_global_parameters, only: unitpar,par_files
        import tc_fluid
        type(tc_fluid), intent(inout) :: fl
      end subroutine read_ffhd_params
    end interface

    type(tc_fluid), intent(inout) :: fl

    fl%tc_k_para=0.d0
    !> Read tc parameters from par file: ffHD case
    call read_ffhd_params(fl)
    if(fl%tc_k_para==0.d0 ) then
      if(SI_unit) then
        ! Spitzer thermal conductivity with SI units
        fl%tc_k_para=8.d-12*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      else
        ! Spitzer thermal conductivity with cgs units
        fl%tc_k_para=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      endif
      if(mype .eq. 0) print*, "Spitzer ffHD par: ",fl%tc_k_para
    else
      fl%tc_constant=.true.
    endif
  end subroutine tc_get_ffhd_params

  !> Get the explicut timestep for the TC (mhd implementation)
  function get_tc_dt_mhd(w,ixI^L,ixO^L,dx^D,x,fl) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters

    type(tc_fluid), intent(in)  ::  fl
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew

    double precision :: mf(ixO^S,1:ndim),Te(ixI^S),rho(ixI^S),gradT(ixI^S)
    double precision :: tmp(ixO^S),hfs(ixO^S)
    double precision :: dtdiff_tcond,maxtmp2
    integer          :: idims,ix^D

    !temperature
    call fl%get_temperature_from_conserved(w,x,ixI^L,ixI^L,Te)

    ! B
   {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if(B0field) then
        ^D&mf({ix^D},^D)=w({ix^D},iw_mag(^D))+block%B0({ix^D},^D,0)\
      else
        ^D&mf({ix^D},^D)=w({ix^D},iw_mag(^D))\
      end if
      ! Bsquared
      tmp(ix^D)=(^D&mf({ix^D},^D)**2+)
      ! B_i**2/B**2
      if(tmp(ix^D)/=0.d0) then
        ^D&mf({ix^D},^D)=mf({ix^D},^D)**2/tmp({ix^D})\
      else
        ^D&mf({ix^D},^D)=1.d0\
      end if
   {end do\}

    dtnew=bigdouble
    call fl%get_rho(w,x,ixI^L,ixO^L,rho)

    !tc_k_para_i
    if(fl%tc_constant) then
      tmp(ixO^S)=fl%tc_k_para
    else
      if(fl%tc_saturate) then
        ! Kannan 2016 MN 458, 410
        ! 3^1.5*kB^2/(4*sqrt(pi)*e^4)
        ! l_mfpe=3.d0**1.5d0*kB_cgs**2/(4.d0*sqrt(dpi)*e_cgs**4*37.d0)=7093.9239487765044d0
        tmp(ixO^S)=Te(ixO^S)**2/rho(ixO^S)*7093.9239487765044d0*unit_temperature**2/(unit_numberdensity*unit_length)
        do idims=1,ndim
          call gradient(Te,ixI^L,ixO^L,idims,gradT)
          if(idims==1) then
            hfs(ixO^S)=gradT(ixO^S)*sqrt(mf(ixO^S,idims))
          else
            hfs(ixO^S)=hfs(ixO^S)+gradT(ixO^S)*sqrt(mf(ixO^S,idims))
          end if
        end do
        ! kappa=kappa_Spizer/(1+4.2*l_mfpe/(T/|gradT.b|))
        tmp(ixO^S)=fl%tc_k_para*dsqrt(Te(ixO^S)**5)/(1.d0+4.2d0*tmp(ixO^S)*dabs(hfs(ixO^S))/Te(ixO^S))
      else
        ! kappa=kappa_Spizer
        tmp(ixO^S)=fl%tc_k_para*dsqrt(Te(ixO^S)**5)
      end if
    end if

    if(slab_uniform) then
      do idims=1,ndim
        ! approximate thermal conduction flux: tc_k_para_i/rho/dx*B_i**2/B**2
        maxtmp2=maxval(tmp(ixO^S)*mf(ixO^S,idims)/(rho(ixO^S)*dxlevel(idims)))
        ! dt< dx_idim**2/((gamma-1)*tc_k_para_i/rho*B_i**2/B**2)
        dtdiff_tcond=dxlevel(idims)/(tc_gamma_1*maxtmp2+smalldouble)
        ! limit the time step
        dtnew=min(dtnew,dtdiff_tcond)
      end do
    else
      do idims=1,ndim
        ! approximate thermal conduction flux: tc_k_para_i/rho/dx*B_i**2/B**2
        maxtmp2=maxval(tmp(ixO^S)*mf(ixO^S,idims)/(rho(ixO^S)*block%ds(ixO^S,idims)**2))
        ! dt< dx_idim**2/((gamma-1)*tc_k_para_i/rho*B_i**2/B**2)
        dtdiff_tcond=1.d0/(tc_gamma_1*maxtmp2+smalldouble)
        ! limit the time step
        dtnew=min(dtnew,dtdiff_tcond)
      end do
    end if
    dtnew=dtnew/dble(ndim)
  end function get_tc_dt_mhd

  !> anisotropic thermal conduction with slope limited symmetric scheme
  !> Sharma 2007 Journal of Computational Physics 227, 123
  subroutine sts_set_source_tc_mhd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,fl)
    use mod_global_parameters
    use mod_fix_conserve
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(in) ::   w(ixI^S,1:nw)
    double precision, intent(inout) ::  wres(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    type(tc_fluid), intent(in) :: fl

    !! qd store the heat conduction energy changing rate
    double precision :: qd(ixO^S)
    double precision :: rho(ixI^S),Te(ixI^S)
    double precision :: qvec(ixI^S,1:ndim)
    double precision :: fluxall(ixI^S,1,1:ndim)
    double precision :: alpha,dxinv(ndim)
    double precision, allocatable, dimension(:^D&,:) :: qvec_equi
    integer :: idims,ixA^L

    ! coefficient of limiting on normal component
    if(ndim<3) then
      alpha=0.75d0
    else
      alpha=0.85d0
    end if

    dxinv=1.d0/dxlevel

    call fl%get_temperature_from_eint(w, x, ixI^L, ixI^L, Te)  !calculate Te in whole domain (+ghosts)
    call fl%get_rho(w, x, ixI^L, ixI^L, rho)  !calculate rho in whole domain (+ghosts)
    call set_source_tc_mhd(ixI^L,ixO^L,w,x,fl,qvec,rho,Te,alpha)
    if(fl%has_equi) then
      allocate(qvec_equi(ixI^S,1:ndim))
      call fl%get_temperature_equi(w, x, ixI^L, ixI^L, Te)  !calculate Te in whole domain (+ghosts)
      call fl%get_rho_equi(w, x, ixI^L, ixI^L, rho)  !calculate rho in whole domain (+ghosts)
      call set_source_tc_mhd(ixI^L,ixO^L,w,x,fl,qvec_equi,rho,Te,alpha)
      do idims=1,ndim
        ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
        qvec(ixA^S,idims)=qvec(ixA^S,idims)-qvec_equi(ixA^S,idims)
      end do
      deallocate(qvec_equi)
    end if

    if(slab_uniform) then
      do idims=1,ndim
        ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
        qvec(ixA^S,idims)=dxinv(idims)*qvec(ixA^S,idims)
        ixA^L=ixO^L-kr(idims,^D);
        if(idims==1) then
          qd(ixO^S)=qvec(ixO^S,idims)-qvec(ixA^S,idims)
        else
          qd(ixO^S)=qd(ixO^S)+qvec(ixO^S,idims)-qvec(ixA^S,idims)
        end if
      end do
    else
      do idims=1,ndim
        ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
        qvec(ixA^S,idims)=qvec(ixA^S,idims)*block%surfaceC(ixA^S,idims)
        ixA^L=ixO^L-kr(idims,^D);
        if(idims==1) then
          qd(ixO^S)=qvec(ixO^S,idims)-qvec(ixA^S,idims)
        else
          qd(ixO^S)=qd(ixO^S)+qvec(ixO^S,idims)-qvec(ixA^S,idims)
        end if
      end do
      qd(ixO^S)=qd(ixO^S)/block%dvolume(ixO^S)
    end if

    if(fix_conserve_at_step) then
      do idims=1,ndim
        ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
        fluxall(ixA^S,1,idims)=my_dt*qvec(ixA^S,idims)
      end do
      call store_flux(igrid,fluxall,1,ndim,nflux)
    end if

    wres(ixO^S,fl%e_)=qd(ixO^S)
  end subroutine sts_set_source_tc_mhd

  subroutine set_source_tc_mhd(ixI^L,ixO^L,w,x,fl,qvec,rho,Te,alpha)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(in) ::  w(ixI^S,1:nw)
    type(tc_fluid), intent(in) :: fl
    double precision, intent(in) :: rho(ixI^S),Te(ixI^S)
    double precision, intent(in) :: alpha
    double precision, intent(out) :: qvec(ixI^S,1:ndim)

    !! qdd store the heat conduction energy changing rate
    double precision, dimension(ixI^S,1:ndim) :: mf,Bc,Bcf,gradT
    double precision, dimension(ixI^S) :: ka,kaf,ke,kef,qdd,Bnorm
    double precision :: minq,maxq,qd(2**(ndim-1))
    integer :: idims,idir,ix^D,ix^L,ixC^L,ixA^L,ixB^L

    ix^L=ixO^L^LADD1;

    ! T gradient at cell faces
    ! B vector
   {do ix^DB=ixmin^DB,ixmax^DB\}
      if(B0field) then
        ^D&mf({ix^D},^D)=w({ix^D},iw_mag(^D))+block%B0({ix^D},^D,0)\
      else
        ^D&mf({ix^D},^D)=w({ix^D},iw_mag(^D))\
      end if
      ! |B|
      Bnorm(ix^D)=dsqrt(^D&mf({ix^D},^D)**2+)
      if(Bnorm(ix^D)/=0.d0) then
        Bnorm(ix^D)=1.d0/Bnorm(ix^D)
      else
        Bnorm(ix^D)=bigdouble
      end if
      ! b unit vector: magnetic field direction vector
      ^D&mf({ix^D},^D)=mf({ix^D},^D)*Bnorm({ix^D})\
   {end do\}
    ! ixC is cell-corner index
    ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;
    ! b unit vector at cell corner
   {^IFTHREED
    do idims=1,3
   {do ix^DB=ixCmin^DB,ixCmax^DB\}
      Bc(ix^D,idims)=0.125d0*(mf(ix1,ix2,ix3,idims)+mf(ix1+1,ix2,ix3,idims)&
                     +mf(ix1,ix2+1,ix3,idims)+mf(ix1+1,ix2+1,ix3,idims)&
                     +mf(ix1,ix2,ix3+1,idims)+mf(ix1+1,ix2,ix3+1,idims)&
                     +mf(ix1,ix2+1,ix3+1,idims)+mf(ix1+1,ix2+1,ix3+1,idims))
   {end do\}
    end do
   }
   {^IFTWOD
    do idims=1,2
   {do ix^DB=ixCmin^DB,ixCmax^DB\}
      Bc(ix^D,idims)=0.25d0*(mf(ix1,ix2,idims)+mf(ix1+1,ix2,idims)&
                     +mf(ix1,ix2+1,idims)+mf(ix1+1,ix2+1,idims))
   {end do\}
    end do
   }
    ! T gradient at cell faces
    do idims=1,ndim
      ixBmin^D=ixmin^D;
      ixBmax^D=ixmax^D-kr(idims,^D);
      call gradientC(Te,ixI^L,ixB^L,idims,gradT(ixI^S,idims))
    end do
    if(fl%tc_constant) then
      if(fl%tc_perpendicular) then
        ka(ixC^S)=fl%tc_k_para-fl%tc_k_perp
        ke(ixC^S)=fl%tc_k_perp
      else
        ka(ixC^S)=fl%tc_k_para
      end if
    else
      ! conductivity at cell center
      if(phys_trac) then
       {do ix^DB=ixmin^DB,ixmax^DB\}
          if(Te(ix^D) < block%wextra(ix^D,fl%Tcoff_)) then
            qdd(ix^D)=fl%tc_k_para*sqrt(block%wextra(ix^D,fl%Tcoff_)**5)
          else
            qdd(ix^D)=fl%tc_k_para*sqrt(Te(ix^D)**5)
          end if
       {end do\}
      else
        qdd(ix^S)=fl%tc_k_para*sqrt(Te(ix^S)**5)
      end if
     ! cell corner parallel conductivity in ka
     {^IFTHREED
     {do ix^DB=ixCmin^DB,ixCmax^DB\}
        ka(ix^D)=0.125d0*(qdd(ix1,ix2,ix3)+qdd(ix1+1,ix2,ix3)&
                       +qdd(ix1,ix2+1,ix3)+qdd(ix1+1,ix2+1,ix3)&
                       +qdd(ix1,ix2,ix3+1)+qdd(ix1+1,ix2,ix3+1)&
                       +qdd(ix1,ix2+1,ix3+1)+qdd(ix1+1,ix2+1,ix3+1))
     {end do\}
     }
     {^IFTWOD
     {do ix^DB=ixCmin^DB,ixCmax^DB\}
        ka(ix^D)=0.25d0*(qdd(ix1,ix2)+qdd(ix1+1,ix2)&
                       +qdd(ix1,ix2+1)+qdd(ix1+1,ix2+1))
     {end do\}
     }
      ! compensate with perpendicular conductivity
      if(fl%tc_perpendicular) then
        qdd(ix^S)=fl%tc_k_perp*rho(ix^S)**2*Bnorm(ix^S)**2/dsqrt(Te(ix^S))
       {^IFTHREED
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          ke(ix^D)=0.125d0*(qdd(ix1,ix2,ix3)+qdd(ix1+1,ix2,ix3)&
                         +qdd(ix1,ix2+1,ix3)+qdd(ix1+1,ix2+1,ix3)&
                         +qdd(ix1,ix2,ix3+1)+qdd(ix1+1,ix2,ix3+1)&
                         +qdd(ix1,ix2+1,ix3+1)+qdd(ix1+1,ix2+1,ix3+1))
          if(ke(ix^D)<ka(ix^D)) then
            ka(ix^D)=ka(ix^D)-ke(ix^D)
          else
            ke(ix^D)=ka(ix^D)
            ka(ix^D)=0.d0
          end if
       {end do\}
       }
       {^IFTWOD
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          ke(ix^D)=0.25d0*(qdd(ix1,ix2)+qdd(ix1+1,ix2)&
                         +qdd(ix1,ix2+1)+qdd(ix1+1,ix2+1))
          if(ke(ix^D)<ka(ix^D)) then
            ka(ix^D)=ka(ix^D)-ke(ix^D)
          else
            ke(ix^D)=ka(ix^D)
            ka(ix^D)=0.d0
          end if
       {end do\}
       }
      end if
    end if
    if(fl%tc_slope_limiter==0) then
      ! calculate thermal conduction flux with symmetric scheme
      do idims=1,ndim
        !qdd corner values
        qdd=0.d0
        {do ix^DB=0,1 \}
           if({ ix^D==0 .and. ^D==idims | .or.}) then
             ixBmin^D=ixCmin^D+ix^D;
             ixBmax^D=ixCmax^D+ix^D;
             qdd(ixC^S)=qdd(ixC^S)+gradT(ixB^S,idims)
           end if
        {end do\}
        ! temperature gradient at cell corner
        qvec(ixC^S,idims)=qdd(ixC^S)*0.5d0**(ndim-1)
      end do
      ! b grad T at cell corner
      qdd(ixC^S)=sum(qvec(ixC^S,1:ndim)*Bc(ixC^S,1:ndim),dim=ndim+1)
      do idims=1,ndim
        ! TC flux at cell corner
        gradT(ixC^S,idims)=ka(ixC^S)*Bc(ixC^S,idims)*qdd(ixC^S)
        if(fl%tc_perpendicular) gradT(ixC^S,idims)=gradT(ixC^S,idims)+ke(ixC^S)*qvec(ixC^S,idims)
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
        if(fl%tc_saturate) then
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
          qdd(ixA^S)=2.75d0*(rho(ixA^S)+rho(ixB^S))*dsqrt(0.5d0*(Te(ixA^S)+Te(ixB^S)))**3*dabs(Bcf(ixA^S,idims))
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            if(dabs(qvec(ix^D,idims))>qdd(ix^D)) then
              qvec(ix^D,idims)=sign(1.d0,qvec(ix^D,idims))*qdd(ix^D)
            end if
         {end do\}
        end if
      end do
    else
      ! calculate thermal conduction flux with slope-limited symmetric scheme
      do idims=1,ndim
        ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
       {^IFTHREED
        if(idims==1) then
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            ! averaged b at face centers
            Bcf(ix^D,1:ndim)=0.25d0*(Bc(ix1,ix2,ix3,1:ndim)+Bc(ix1,ix2-1,ix3,1:ndim)&
                           +Bc(ix1,ix2,ix3-1,1:ndim)+Bc(ix1,ix2-1,ix3-1,1:ndim))
            kaf(ix^D)=0.25d0*(ka(ix1,ix2,ix3)+ka(ix1,ix2-1,ix3)&
                           +ka(ix1,ix2,ix3-1)+ka(ix1,ix2-1,ix3-1))
            ! averaged thermal conductivity at face centers
            if(fl%tc_perpendicular) &
            kef(ix^D)=0.25d0*(ke(ix1,ix2,ix3)+ke(ix1,ix2-1,ix3)&
                           +ke(ix1,ix2,ix3-1)+ke(ix1,ix2-1,ix3-1))
         {end do\}
        else if(idims==2) then
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            Bcf(ix^D,1:ndim)=0.25d0*(Bc(ix1,ix2,ix3,1:ndim)+Bc(ix1-1,ix2,ix3,1:ndim)&
                           +Bc(ix1,ix2,ix3-1,1:ndim)+Bc(ix1-1,ix2,ix3-1,1:ndim))
            kaf(ix^D)=0.25d0*(ka(ix1,ix2,ix3)+ka(ix1-1,ix2,ix3)&
                           +ka(ix1,ix2,ix3-1)+ka(ix1-1,ix2,ix3-1))
            if(fl%tc_perpendicular) &
            kef(ix^D)=0.25d0*(ke(ix1,ix2,ix3)+ke(ix1-1,ix2,ix3)&
                           +ke(ix1,ix2,ix3-1)+ke(ix1-1,ix2,ix3-1))
         {end do\}
        else
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            Bcf(ix^D,1:ndim)=0.25d0*(Bc(ix1,ix2,ix3,1:ndim)+Bc(ix1,ix2-1,ix3,1:ndim)&
                           +Bc(ix1-1,ix2,ix3,1:ndim)+Bc(ix1-1,ix2-1,ix3,1:ndim))
            kaf(ix^D)=0.25d0*(ka(ix1,ix2,ix3)+ka(ix1,ix2-1,ix3)&
                           +ka(ix1-1,ix2,ix3)+ka(ix1-1,ix2-1,ix3))
            if(fl%tc_perpendicular) &
            kef(ix^D)=0.25d0*(ke(ix1,ix2,ix3)+ke(ix1,ix2-1,ix3)&
                           +ke(ix1-1,ix2,ix3)+ke(ix1-1,ix2-1,ix3))
         {end do\}
        end if
       }
       {^IFTWOD
        if(idims==1) then
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            Bcf(ix^D,1:ndim)=0.5d0*(Bc(ix1,ix2,1:ndim)+Bc(ix1,ix2-1,1:ndim))
            kaf(ix^D)=0.5d0*(ka(ix1,ix2)+ka(ix1,ix2-1))
            if(fl%tc_perpendicular) &
            kef(ix^D)=0.5d0*(ke(ix1,ix2)+ke(ix1,ix2-1))
         {end do\}
        else
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            Bcf(ix^D,1:ndim)=0.5d0*(Bc(ix1,ix2,1:ndim)+Bc(ix1-1,ix2,1:ndim))
            kaf(ix^D)=0.5d0*(ka(ix1,ix2)+ka(ix1-1,ix2))
            if(fl%tc_perpendicular) &
            kef(ix^D)=0.5d0*(ke(ix1,ix2)+ke(ix1-1,ix2))
         {end do\}
        end if
       }
        ! eq (19)
        ! temperature gradient at cell corner
       {^IFTHREED
        if(idims==1) then
         {do ix^DB=ixCmin^DB,ixCmax^DB\}
            qdd(ix^D)=0.25d0*(gradT(ix1,ix2,ix3,idims)+gradT(ix1,ix2+1,ix3,idims)&
                           +gradT(ix1,ix2,ix3+1,idims)+gradT(ix1,ix2+1,ix3+1,idims))
         {end do\}
        else if(idims==2) then
         {do ix^DB=ixCmin^DB,ixCmax^DB\}
            qdd(ix^D)=0.25d0*(gradT(ix1,ix2,ix3,idims)+gradT(ix1+1,ix2,ix3,idims)&
                           +gradT(ix1,ix2,ix3+1,idims)+gradT(ix1+1,ix2,ix3+1,idims))
         {end do\}
        else
         {do ix^DB=ixCmin^DB,ixCmax^DB\}
            qdd(ix^D)=0.25d0*(gradT(ix1,ix2,ix3,idims)+gradT(ix1+1,ix2,ix3,idims)&
                           +gradT(ix1,ix2+1,ix3,idims)+gradT(ix1+1,ix2+1,ix3,idims))
         {end do\}
        end if
       }
       {^IFTWOD
        if(idims==1) then
         {do ix^DB=ixCmin^DB,ixCmax^DB\}
            qdd(ix^D)=0.5d0*(gradT(ix1,ix2,idims)+gradT(ix1,ix2+1,idims))
         {end do\}
        else
         {do ix^DB=ixCmin^DB,ixCmax^DB\}
            qdd(ix^D)=0.5d0*(gradT(ix1,ix2,idims)+gradT(ix1+1,ix2,idims))
         {end do\}
        end if
       }
        ! eq (21)
       {^IFTHREED
        if(idims==1) then
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
            maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
            if(qdd(ix^D)<minq) then
              qd(1)=minq
            else if(qdd(ix^D)>maxq) then
              qd(1)=maxq
            else
              qd(1)=qdd(ix^D)
            end if
            if(qdd(ix1,ix2-1,ix3)<minq) then
              qd(2)=minq
            else if(qdd(ix1,ix2-1,ix3)>maxq) then
              qd(2)=maxq
            else
              qd(2)=qdd(ix1,ix2-1,ix3)
            end if
            if(qdd(ix1,ix2,ix3-1)<minq) then
              qd(3)=minq
            else if(qdd(ix1,ix2,ix3-1)>maxq) then
              qd(3)=maxq
            else
              qd(3)=qdd(ix1,ix2,ix3-1)
            end if
            if(qdd(ix1,ix2-1,ix3-1)<minq) then
              qd(4)=minq
            else if(qdd(ix1,ix2-1,ix3-1)>maxq) then
              qd(4)=maxq
            else
              qd(4)=qdd(ix1,ix2-1,ix3-1)
            end if
            qvec(ix^D,idims)=kaf(ix^D)*0.25d0*(Bc(ix^D,idims)**2*qd(1)+Bc(ix1,ix2-1,ix3,idims)**2*qd(2)&
                           +Bc(ix1,ix2,ix3-1,idims)**2*qd(3)+Bc(ix1,ix2-1,ix3-1,idims)**2*qd(4))
            if(fl%tc_perpendicular) &
            qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.25d0*(qd(1)+qd(2)+qd(3)+qd(4))
         {end do\}
        else if(idims==2) then
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
            maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
            if(qdd(ix^D)<minq) then
              qd(1)=minq
            else if(qdd(ix^D)>maxq) then
              qd(1)=maxq
            else
              qd(1)=qdd(ix^D)
            end if
            if(qdd(ix1-1,ix2,ix3)<minq) then
              qd(2)=minq
            else if(qdd(ix1-1,ix2,ix3)>maxq) then
              qd(2)=maxq
            else
              qd(2)=qdd(ix1-1,ix2,ix3)
            end if
            if(qdd(ix1,ix2,ix3-1)<minq) then
              qd(3)=minq
            else if(qdd(ix1,ix2,ix3-1)>maxq) then
              qd(3)=maxq
            else
              qd(3)=qdd(ix1,ix2,ix3-1)
            end if
            if(qdd(ix1-1,ix2,ix3-1)<minq) then
              qd(4)=minq
            else if(qdd(ix1-1,ix2,ix3-1)>maxq) then
              qd(4)=maxq
            else
              qd(4)=qdd(ix1-1,ix2,ix3-1)
            end if
            qvec(ix^D,idims)=kaf(ix^D)*0.25d0*(Bc(ix^D,idims)**2*qd(1)+Bc(ix1-1,ix2,ix3,idims)**2*qd(2)&
                           +Bc(ix1,ix2,ix3-1,idims)**2*qd(3)+Bc(ix1-1,ix2,ix3-1,idims)**2*qd(4))
            if(fl%tc_perpendicular) &
            qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.25d0*(qd(1)+qd(2)+qd(3)+qd(4))
         {end do\}
        else
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
            maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
            if(qdd(ix^D)<minq) then
              qd(1)=minq
            else if(qdd(ix^D)>maxq) then
              qd(1)=maxq
            else
              qd(1)=qdd(ix^D)
            end if
            if(qdd(ix1-1,ix2,ix3)<minq) then
              qd(2)=minq
            else if(qdd(ix1-1,ix2,ix3)>maxq) then
              qd(2)=maxq
            else
              qd(2)=qdd(ix1-1,ix2,ix3)
            end if
            if(qdd(ix1,ix2-1,ix3)<minq) then
              qd(3)=minq
            else if(qdd(ix1,ix2-1,ix3)>maxq) then
              qd(3)=maxq
            else
              qd(3)=qdd(ix1,ix2-1,ix3)
            end if
            if(qdd(ix1-1,ix2-1,ix3)<minq) then
              qd(4)=minq
            else if(qdd(ix1-1,ix2-1,ix3)>maxq) then
              qd(4)=maxq
            else
              qd(4)=qdd(ix1-1,ix2-1,ix3)
            end if
            qvec(ix^D,idims)=kaf(ix^D)*0.25d0*(Bc(ix^D,idims)**2*qd(1)+Bc(ix1-1,ix2,ix3,idims)**2*qd(2)&
                           +Bc(ix1,ix2-1,ix3,idims)**2*qd(3)+Bc(ix1-1,ix2-1,ix3,idims)**2*qd(4))
            if(fl%tc_perpendicular) &
            qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.25d0*(qd(1)+qd(2)+qd(3)+qd(4))
         {end do\}
        end if
       }
       {^IFTWOD
        if(idims==1) then
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
            maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
            if(qdd(ix^D)<minq) then
              qd(1)=minq
            else if(qdd(ix^D)>maxq) then
              qd(1)=maxq
            else
              qd(1)=qdd(ix^D)
            end if
            if(qdd(ix1,ix2-1)<minq) then
              qd(2)=minq
            else if(qdd(ix1,ix2-1)>maxq) then
              qd(2)=maxq
            else
              qd(2)=qdd(ix1,ix2-1)
            end if
            qvec(ix^D,idims)=kaf(ix^D)*0.5d0*(Bc(ix1,ix2,idims)**2*qd(1)+Bc(ix1,ix2-1,idims)**2*qd(2))
            if(fl%tc_perpendicular) &
            qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.5d0*(qd(1)+qd(2))
         {end do\}
        else
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
            maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
            if(qdd(ix^D)<minq) then
              qd(1)=minq
            else if(qdd(ix^D)>maxq) then
              qd(1)=maxq
            else
              qd(1)=qdd(ix^D)
            end if
            if(qdd(ix1-1,ix2)<minq) then
              qd(2)=minq
            else if(qdd(ix1-1,ix2)>maxq) then
              qd(2)=maxq
            else
              qd(2)=qdd(ix1-1,ix2)
            end if
            qvec(ix^D,idims)=kaf(ix^D)*0.5d0*(Bc(ix1,ix2,idims)**2*qd(1)+Bc(ix1-1,ix2,idims)**2*qd(2))
            if(fl%tc_perpendicular) &
            qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.5d0*(qd(1)+qd(2))
         {end do\}
        end if
       }
        ! calculate normal of magnetic field
        ixB^L=ixA^L+kr(idims,^D);
        Bnorm(ixA^S)=0.5d0*(mf(ixA^S,idims)+mf(ixB^S,idims))
        ! limited transverse component, eq (17)
        ixBmin^D=ixAmin^D;
        ixBmax^D=ixAmax^D+kr(idims,^D);
        do idir=1,ndim
          if(idir==idims) cycle
          qdd(ixI^S)=slope_limiter(gradT(ixI^S,idir),ixI^L,ixB^L,idir,-1,fl%tc_slope_limiter)
          qdd(ixI^S)=slope_limiter(qdd,ixI^L,ixA^L,idims,1,fl%tc_slope_limiter)
          qvec(ixA^S,idims)=qvec(ixA^S,idims)+kaf(ixA^S)*Bnorm(ixA^S)*Bcf(ixA^S,idir)*qdd(ixA^S)
        end do
        if(fl%tc_saturate) then
          ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135)
          ! unsigned saturated TC flux = 5 phi rho c**3, c=sqrt(p/rho) is isothermal sound speed, phi=1.1
          ixB^L=ixA^L+kr(idims,^D);
          qdd(ixA^S)=2.75d0*(rho(ixA^S)+rho(ixB^S))*dsqrt(0.5d0*(Te(ixA^S)+Te(ixB^S)))**3*dabs(Bnorm(ixA^S))
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            if(dabs(qvec(ix^D,idims))>qdd(ix^D)) then
              qvec(ix^D,idims)=sign(1.d0,qvec(ix^D,idims))*qdd(ix^D)
            end if
         {end do\}
        end if
      end do
    end if
  end subroutine set_source_tc_mhd

  function slope_limiter(f,ixI^L,ixO^L,idims,pm,tc_slope_limiter) result(lf)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L, idims, pm
    double precision, intent(in) :: f(ixI^S)
    double precision :: lf(ixI^S)
    integer, intent(in)  :: tc_slope_limiter

    double precision :: signf(ixI^S)
    integer :: ixB^L

    ixB^L=ixO^L+pm*kr(idims,^D);
    signf(ixO^S)=sign(1.d0,f(ixO^S))
    select case(tc_slope_limiter)
     case(1)
       ! 'MC' montonized central limiter Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       lf(ixO^S)=two*signf(ixO^S)* &
            max(zero,min(dabs(f(ixO^S)),signf(ixO^S)*f(ixB^S),&
            signf(ixO^S)*quarter*(f(ixB^S)+f(ixO^S))))
     case(2)
       ! 'minmod' limiter
       lf(ixO^S)=signf(ixO^S)*max(0.d0,min(abs(f(ixO^S)),signf(ixO^S)*f(ixB^S)))
     case(3)
       ! 'superbee' Roes superbee limiter (eq.3.51i)
       lf(ixO^S)=signf(ixO^S)* &
            max(zero,min(two*dabs(f(ixO^S)),signf(ixO^S)*f(ixB^S)),&
            min(dabs(f(ixO^S)),two*signf(ixO^S)*f(ixB^S)))
     case(4)
       ! 'koren' Barry Koren Right variant
       lf(ixO^S)=signf(ixO^S)* &
            max(zero,min(two*dabs(f(ixO^S)),two*signf(ixO^S)*f(ixB^S),&
            (two*f(ixB^S)*signf(ixO^S)+dabs(f(ixO^S)))*third))
     case default
       call mpistop("Unknown slope limiter for thermal conduction")
    end select
  end function slope_limiter

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

  !> Get the explicit timestep for the TC (hd implementation)
  function get_tc_dt_hd(w,ixI^L,ixO^L,dx^D,x,fl)  result(dtnew)
    ! Check diffusion time limit dt < dx_i**2 / ((gamma-1)*tc_k_para_i/rho)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    type(tc_fluid), intent(in) :: fl
    double precision :: dtnew

    double precision :: tmp(ixO^S),tmp2(ixO^S),Te(ixI^S),rho(ixI^S),hfs(ixO^S),gradT(ixI^S)
    double precision :: dtdiff_tcond,maxtmp2
    integer          :: idim

    call fl%get_temperature_from_conserved(w,x,ixI^L,ixI^L,Te)
    call fl%get_rho(w,x,ixI^L,ixO^L,rho)

    if(fl%tc_saturate) then
      ! Kannan 2016 MN 458, 410
      ! 3^1.5*kB^2/(4*sqrt(pi)*e^4)
      ! l_mfpe=3.d0**1.5d0*kB_cgs**2/(4.d0*sqrt(dpi)*e_cgs**4*37.d0)=7093.9239487765044d0
      tmp2(ixO^S)=Te(ixO^S)**2/rho(ixO^S)*7093.9239487765044d0*unit_temperature**2/(unit_numberdensity*unit_length)
      hfs=0.d0
      do idim=1,ndim
        call gradient(Te,ixI^L,ixO^L,idim,gradT)
        hfs(ixO^S)=hfs(ixO^S)+gradT(ixO^S)**2
      end do
      ! kappa=kappa_Spizer/(1+4.2*l_mfpe/(T/|gradT|))
      tmp(ixO^S)=fl%tc_k_para*dsqrt((Te(ixO^S))**5)/(rho(ixO^S)*(1.d0+4.2d0*tmp2(ixO^S)*sqrt(hfs(ixO^S))/Te(ixO^S)))
    else
      tmp(ixO^S)=fl%tc_k_para*dsqrt((Te(ixO^S))**5)/rho(ixO^S)
    end if
    dtnew = bigdouble

    if(slab_uniform) then
      do idim=1,ndim
        ! approximate thermal conduction flux: tc_k_para_i/rho/dx
        tmp2(ixO^S)=tmp(ixO^S)/dxlevel(idim)
        maxtmp2=maxval(tmp2(ixO^S))
        ! dt< dx_idim**2/((gamma-1)*tc_k_para_i/rho)
        dtdiff_tcond=dxlevel(idim)/(tc_gamma_1*maxtmp2+smalldouble)
        ! limit the time step
        dtnew=min(dtnew,dtdiff_tcond)
      end do
    else
      do idim=1,ndim
        ! approximate thermal conduction flux: tc_k_para_i/rho/dx
        tmp2(ixO^S)=tmp(ixO^S)/block%ds(ixO^S,idim)
        maxtmp2=maxval(tmp2(ixO^S)/block%ds(ixO^S,idim))
        ! dt< dx_idim**2/((gamma-1)*tc_k_para_i/rho*B_i**2/B**2)
        dtdiff_tcond=1.d0/(tc_gamma_1*maxtmp2+smalldouble)
        ! limit the time step
        dtnew=min(dtnew,dtdiff_tcond)
      end do
    end if
    dtnew=dtnew/dble(ndim)
  end function get_tc_dt_hd

  subroutine sts_set_source_tc_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,fl)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(in) ::  w(ixI^S,1:nw)
    double precision, intent(inout) ::  wres(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    type(tc_fluid), intent(in)    :: fl

    double precision :: Te(ixI^S),rho(ixI^S)
    double precision :: qvec(ixI^S,1:ndim),qd(ixI^S)
    double precision, allocatable, dimension(:^D&,:) :: qvec_equi
    double precision, allocatable, dimension(:^D&,:,:) :: fluxall

    double precision :: dxinv(ndim)
    integer :: idims,ix^L,ixB^L

    ix^L=ixO^L^LADD1;

    dxinv=1.d0/dxlevel

    !calculate Te in whole domain (+ghosts)
    call fl%get_temperature_from_eint(w, x, ixI^L, ixI^L, Te)
    call fl%get_rho(w, x, ixI^L, ixI^L, rho)
    call set_source_tc_hd(ixI^L,ixO^L,w,x,fl,qvec,rho,Te)
    if(fl%has_equi) then
      allocate(qvec_equi(ixI^S,1:ndim))
      call fl%get_temperature_equi(w, x, ixI^L, ixI^L, Te)  !calculate Te in whole domain (+ghosts)
      call fl%get_rho_equi(w, x, ixI^L, ixI^L, rho)  !calculate rho in whole domain (+ghosts)
      call set_source_tc_hd(ixI^L,ixO^L,w,x,fl,qvec_equi,rho,Te)
      do idims=1,ndim
        qvec(ix^S,idims)=qvec(ix^S,idims) - qvec_equi(ix^S,idims)
      end do
      deallocate(qvec_equi)
    endif

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

    if(fix_conserve_at_step) then
      allocate(fluxall(ixI^S,1,1:ndim))
      fluxall(ixI^S,1,1:ndim)=my_dt*qvec(ixI^S,1:ndim)
      call store_flux(igrid,fluxall,1,ndim,nflux)
      deallocate(fluxall)
    end if
    wres(ixO^S,fl%e_)=qd(ixO^S)
  end subroutine sts_set_source_tc_hd

  subroutine set_source_tc_hd(ixI^L,ixO^L,w,x,fl,qvec,rho,Te)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(in) ::  w(ixI^S,1:nw)
    type(tc_fluid), intent(in)    :: fl
    double precision, intent(in) :: Te(ixI^S),rho(ixI^S)
    double precision, intent(out) :: qvec(ixI^S,1:ndim)
    double precision :: gradT(ixI^S,1:ndim),ke(ixI^S),qd(ixI^S)
    integer :: idims,ix^D,ix^L,ixC^L,ixA^L,ixB^L,ixD^L

    ix^L=ixO^L^LADD1;
    ! ixC is cell-corner index
    ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;

    ! calculate thermal conduction flux with symmetric scheme
    ! T gradient (central difference) at cell corners
    do idims=1,ndim
      ixBmin^D=ixmin^D;
      ixBmax^D=ixmax^D-kr(idims,^D);
      call gradientC(Te,ixI^L,ixB^L,idims,ke)
      qd=0.d0
     {do ix^DB=0,1 \}
        if({ix^D==0 .and. ^D==idims |.or. }) then
          ixBmin^D=ixCmin^D+ix^D;
          ixBmax^D=ixCmax^D+ix^D;
          qd(ixC^S)=qd(ixC^S)+ke(ixB^S)
        end if
     {end do\}
      ! temperature gradient at cell corner
      qvec(ixC^S,idims)=qd(ixC^S)*0.5d0**(ndim-1)
    end do
    ! conductivity at cell center
    if(phys_trac) then
      ! transition region adaptive conduction
      where(Te(ix^S) < block%wextra(ix^S,fl%Tcoff_))
        qd(ix^S)=fl%tc_k_para*dsqrt(block%wextra(ix^S,fl%Tcoff_))**5
      else where
        qd(ix^S)=fl%tc_k_para*dsqrt(Te(ix^S))**5
      end where
    else
      qd(ix^S)=fl%tc_k_para*dsqrt(Te(ix^S))**5
    end if
    ke=0.d0
    {do ix^DB=0,1\}
      ixBmin^D=ixCmin^D+ix^D;
      ixBmax^D=ixCmax^D+ix^D;
      ke(ixC^S)=ke(ixC^S)+qd(ixB^S)
    {end do\}
    ! cell corner conductivity
    ke(ixC^S)=0.5d0**ndim*ke(ixC^S)
    ! cell corner conduction flux
    do idims=1,ndim
      gradT(ixC^S,idims)=ke(ixC^S)*qvec(ixC^S,idims)
    end do

    if(fl%tc_saturate) then
      ! consider saturation with unsigned saturated TC flux = 5 phi rho c**3
      ! saturation flux at cell center
      qd(ix^S)=5.5d0*rho(ix^S)*dsqrt(Te(ix^S)**3)
      !cell corner values of qd in ke
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
          ke(ix^D)=ke(ix^D)/(qd(ix^D)+smalldouble)
          do idims=1,ndim
            gradT(ix^D,idims)=ke(ix^D)*gradT(ix^D,idims)
          end do
        end if
      {end do\}
    end if

    ! conduction flux at cell face
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
  end subroutine set_source_tc_hd

  !> Get the explicut timestep for the TC (ffhd implementation)
  function get_tc_dt_ffhd(w,ixI^L,ixO^L,dx^D,x,fl) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    type(tc_fluid), intent(in) :: fl
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew
    double precision :: mf(ixO^S,1:ndim),Te(ixI^S),rho(ixI^S),gradT(ixI^S)
    double precision :: tmp2(ixO^S),tmp(ixO^S),hfs(ixO^S)
    double precision :: dtdiff_tcond,maxtmp2
    integer          :: idims

    !temperature
    call fl%get_temperature_from_conserved(w,x,ixI^L,ixI^L,Te)
    mf(ixO^S,1:ndim)=(block%B0(ixO^S,1:ndim,0))**2

    dtnew=bigdouble
    call fl%get_rho(w,x,ixI^L,ixO^L,rho)

    !tc_k_para_i
    if(fl%tc_constant) then
      tmp(ixO^S)=fl%tc_k_para
    else
      if(fl%tc_saturate) then
        ! Kannan 2016 MN 458, 410
        ! 3^1.5*kB^2/(4*sqrt(pi)*e^4)
        ! l_mfpe=3.d0**1.5d0*kB_cgs**2/(4.d0*sqrt(dpi)*e_cgs**4*37.d0)=7093.9239487765044d0
        tmp2(ixO^S)=Te(ixO^S)**2/rho(ixO^S)*7093.9239487765044d0*unit_temperature**2/(unit_numberdensity*unit_length)
        hfs=0.d0
        do idims=1,ndim
          call gradient(Te,ixI^L,ixO^L,idims,gradT)
          hfs(ixO^S)=hfs(ixO^S)+gradT(ixO^S)*sqrt(mf(ixO^S,idims))
        end do
        ! kappa=kappa_Spizer/(1+4.2*l_mfpe/(T/|gradT.b|))
        tmp(ixO^S)=fl%tc_k_para*dsqrt(Te(ixO^S)**5)/(1.d0+4.2d0*tmp2(ixO^S)*dabs(hfs(ixO^S))/Te(ixO^S))
      else
        ! kappa=kappa_Spizer
        tmp(ixO^S)=fl%tc_k_para*dsqrt(Te(ixO^S)**5)
      end if
    end if

    if(slab_uniform) then
      do idims=1,ndim
        ! approximate thermal conduction flux: tc_k_para_i/rho/dx*B_i**2/B**2
        tmp2(ixO^S)=tmp(ixO^S)*mf(ixO^S,idims)/(rho(ixO^S)*dxlevel(idims))
        maxtmp2=maxval(tmp2(ixO^S))
        ! dt< dx_idim**2/((gamma-1)*tc_k_para_i/rho*B_i**2/B**2)
        dtdiff_tcond=dxlevel(idims)/(tc_gamma_1*maxtmp2+smalldouble)
        ! limit the time step
        dtnew=min(dtnew,dtdiff_tcond)
      end do
    else
      do idims=1,ndim
        ! approximate thermal conduction flux: tc_k_para_i/rho/dx*B_i**2/B**2
        tmp2(ixO^S)=tmp(ixO^S)*mf(ixO^S,idims)/(rho(ixO^S)*block%ds(ixO^S,idims))
        maxtmp2=maxval(tmp2(ixO^S)/block%ds(ixO^S,idims))
        ! dt< dx_idim**2/((gamma-1)*tc_k_para_i/rho*B_i**2/B**2)
        dtdiff_tcond=1.d0/(tc_gamma_1*maxtmp2+smalldouble)
        ! limit the time step
        dtnew=min(dtnew,dtdiff_tcond)
      end do
    end if
    dtnew=dtnew/dble(ndim)
  end function get_tc_dt_ffhd

  !> anisotropic thermal conduction with slope limited symmetric scheme
  !> Sharma 2007 Journal of Computational Physics 227, 123
  subroutine sts_set_source_tc_ffhd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,fl)
    use mod_global_parameters
    use mod_fix_conserve
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(in) ::   w(ixI^S,1:nw)
    double precision, intent(inout) ::  wres(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    type(tc_fluid), intent(in) :: fl
    !! qd store the heat conduction energy changing rate
    double precision :: qd(ixI^S)
    double precision :: rho(ixI^S),Te(ixI^S)
    double precision :: qvec(ixI^S,1:ndim)
    double precision, allocatable, dimension(:^D&,:) :: qvec_equi
    double precision, allocatable, dimension(:^D&,:,:) :: fluxall
    double precision :: alpha,dxinv(ndim)
    integer :: idims,idir,ix^D,ix^L,ixC^L,ixA^L,ixB^L

    ! coefficient of limiting on normal component
    if(ndim<3) then
      alpha=0.75d0
    else
      alpha=0.85d0
    end if
    ix^L=ixO^L^LADD1;
    dxinv=1.d0/dxlevel
    call fl%get_temperature_from_eint(w, x, ixI^L, ixI^L, Te)  !calculate Te in whole domain (+ghosts)
    call fl%get_rho(w, x, ixI^L, ixI^L, rho)  !calculate rho in whole domain (+ghosts)
    call set_source_tc_ffhd(ixI^L,ixO^L,w,x,fl,qvec,rho,Te,alpha)

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

    if(fix_conserve_at_step) then
      allocate(fluxall(ixI^S,1,1:ndim))
      fluxall(ixI^S,1,1:ndim)=my_dt*qvec(ixI^S,1:ndim)
      call store_flux(igrid,fluxall,1,ndim,nflux)
      deallocate(fluxall)
    end if
    wres(ixO^S,fl%e_)=qd(ixO^S)
  end subroutine sts_set_source_tc_ffhd

  subroutine set_source_tc_ffhd(ixI^L,ixO^L,w,x,fl,qvec,rho,Te,alpha)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(in) ::  w(ixI^S,1:nw)
    type(tc_fluid), intent(in) :: fl
    double precision, intent(in) :: rho(ixI^S),Te(ixI^S)
    double precision, intent(in) :: alpha
    double precision, intent(out) :: qvec(ixI^S,1:ndim)
    !! qd store the heat conduction energy changing rate
    double precision :: qd(ixI^S)
    double precision, dimension(ixI^S,1:ndim) :: mf,Bc,Bcf
    double precision, dimension(ixI^S,1:ndim) :: gradT
    double precision, dimension(ixI^S) :: ka,kaf,ke,kef,qdd,qe,Binv,minq,maxq,Bnorm
    double precision, allocatable, dimension(:^D&,:,:) :: fluxall
    integer, dimension(ndim) :: lowindex
    integer :: idims,idir,ix^D,ix^L,ixC^L,ixA^L,ixB^L,ixA^D,ixB^D

    ix^L=ixO^L^LADD1;

    ! T gradient at cell faces
    ! B vector
    mf(ixI^S,1:ndim)=block%B0(ixI^S,1:ndim,0)
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
    do idims=1,ndim
      ixBmin^D=ixmin^D;
      ixBmax^D=ixmax^D-kr(idims,^D);
      call gradientC(Te,ixI^L,ixB^L,idims,minq)
      gradT(ixB^S,idims)=minq(ixB^S)
    end do
    if(fl%tc_constant) then
      ka(ixC^S)=fl%tc_k_para
    else
      ! conductivity at cell center
      if(phys_trac) then
        minq(ix^S)=Te(ix^S)
       {do ix^DB=ixmin^DB,ixmax^DB\}
          if(minq(ix^D) < block%wextra(ix^D,fl%Tcoff_)) then
            minq(ix^D)=block%wextra(ix^D,fl%Tcoff_)
          end if
       {end do\}
        minq(ix^S)=fl%tc_k_para*sqrt(minq(ix^S)**5)
      else
        minq(ix^S)=fl%tc_k_para*sqrt(Te(ix^S)**5)
      end if
      ka=0.d0
      {do ix^DB=0,1\}
        ixBmin^D=ixCmin^D+ix^D;
        ixBmax^D=ixCmax^D+ix^D;
        ka(ixC^S)=ka(ixC^S)+minq(ixB^S)
      {end do\}
      ! cell corner conductivity
      ka(ixC^S)=0.5d0**ndim*ka(ixC^S)
    end if
    if(fl%tc_slope_limiter==0) then
      ! calculate thermal conduction flux with symmetric scheme
      do idims=1,ndim
        !qd corner values
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
      qd(ixC^S)=0.d0
      do idims=1,ndim
       qd(ixC^S)=qvec(ixC^S,idims)*Bc(ixC^S,idims)+qd(ixC^S)
      end do
      do idims=1,ndim
        ! TC flux at cell corner
        gradT(ixC^S,idims)=ka(ixC^S)*Bc(ixC^S,idims)*qd(ixC^S)
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
        if(fl%tc_saturate) then
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
          qd(ixA^S)=2.75d0*(rho(ixA^S)+rho(ixB^S))*dsqrt(0.5d0*(Te(ixA^S)+Te(ixB^S)))**3*dabs(Bcf(ixA^S,idims))
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
           end if
        {end do\}
        ! averaged b at face centers
        Bcf(ixA^S,1:ndim)=Bcf(ixA^S,1:ndim)*0.5d0**(ndim-1)
        ! averaged thermal conductivity at face centers
        kaf(ixA^S)=kaf(ixA^S)*0.5d0**(ndim-1)
        ! limited normal component
        minq(ixA^S)=min(alpha*gradT(ixA^S,idims),gradT(ixA^S,idims)/alpha)
        maxq(ixA^S)=max(alpha*gradT(ixA^S,idims),gradT(ixA^S,idims)/alpha)
        ! eq (19)
        !corner values of gradT
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
            {do ixA^DB=ixAmin^DB,ixAmax^DB
               ixB^DB=ixA^DB-ix^DB\}
               if(qd(ixB^D)<=minq(ixA^D)) then
                 qd(ixB^D)=minq(ixA^D)
               else if(qd(ixB^D)>=maxq(ixA^D)) then
                 qd(ixB^D)=maxq(ixA^D)
               end if 
            {end do\}
             qvec(ixA^S,idims)=qvec(ixA^S,idims)+Bc(ixB^S,idims)**2*qd(ixB^S)
           end if
        {end do\}
        qvec(ixA^S,idims)=kaf(ixA^S)*qvec(ixA^S,idims)*0.5d0**(ndim-1)
        ! limited transverse component, eq (17)
        ixBmin^D=ixAmin^D;
        ixBmax^D=ixAmax^D+kr(idims,^D);
        do idir=1,ndim
          if(idir==idims) cycle
          qd(ixI^S)=slope_limiter(gradT(ixI^S,idir),ixI^L,ixB^L,idir,-1,fl%tc_slope_limiter)
          qd(ixI^S)=slope_limiter(qd,ixI^L,ixA^L,idims,1,fl%tc_slope_limiter)
          qvec(ixA^S,idims)=qvec(ixA^S,idims)+kaf(ixA^S)*Bnorm(ixA^S)*Bcf(ixA^S,idir)*qd(ixA^S)
        end do
        if(fl%tc_saturate) then
          ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135)
          ! unsigned saturated TC flux = 5 phi rho c**3, c=sqrt(p/rho) is isothermal sound speed, phi=1.1
          ixB^L=ixA^L+kr(idims,^D);
          qd(ixA^S)=2.75d0*(rho(ixA^S)+rho(ixB^S))*dsqrt(0.5d0*(Te(ixA^S)+Te(ixB^S)))**3*dabs(Bnorm(ixA^S))
         {do ix^DB=ixAmin^DB,ixAmax^DB\}
            if(dabs(qvec(ix^D,idims))>qd(ix^D)) then
              !write(*,*) 'it',it,qvec(ix^D,idims),qd(ix^D),' TC saturated at ',&
              !x(ix^D,:),' rho',rho(ix^D),' Te',Te(ix^D)
              qvec(ix^D,idims)=sign(1.d0,qvec(ix^D,idims))*qd(ix^D)
            end if
         {end do\}
        end if
      end do
    end if
  end subroutine set_source_tc_ffhd
end module mod_thermal_conduction
