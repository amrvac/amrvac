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
!> 30.11.2025 Minor cleanup (for consistency between hd and mhd)
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

    subroutine get_2var_subr(ixI^L, ixO^L, w, x, ne, nH)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: ne(ixI^S), nH(ixI^S)
    end subroutine get_2var_subr

    !> Scalar EoS inverse, e.g. eint_nH_from_T(log_nH, log_T)
    double precision function eos_scalar2_func(a, b)
      double precision, intent(in) :: a, b
    end function eos_scalar2_func
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

    ! if subtract_equi = .true. get_temperature_equi and get_rho_equi have to be set
    logical :: subtract_equi=.false.

    !> Logical switch for test constant conductivity
    logical :: tc_constant=.false.

    !> Calculate thermal conduction perpendicular to magnetic field (.true.) or not (.false.)
    logical :: tc_perpendicular=.false.

    !> Consider thermal conduction saturation effect (.true.) or not (.false.)
    logical :: tc_saturate=.false.

    !> Patch cells where e_int < 0 during STS Chebyshev substeps by
    !> neighbor-averaging the temperature. Prevents NaN from sqrt(T<0)
    !> in the Spitzer conductivity when the RKL2 polynomial overshoots.
    logical :: tc_patch_eint=.false.

    !> Minimum temperature (code units) below which TRAC does not modify conductivity.
    !> Below this T, the energy balance is dominated by optically thick radiation
    !> and recombination, not optically thin cooling + Spitzer conduction, so TRAC's
    !> broadening assumption does not apply. Read in Kelvin via tc_list (trac_T_floor),
    !> converted to code units during init. Default 1e4 K.
    double precision :: trac_T_floor=0.d0
    ! END the following are read from param file or set in tc_read_hd_params or tc_read_mhd_params
    procedure (get_var_subr), pointer, nopass :: get_rho => null()
    procedure (get_var_subr), pointer, nopass :: get_rho_equi => null()
    procedure(get_var_subr), pointer,nopass :: get_temperature_from_eint => null()
    procedure(get_var_subr), pointer,nopass :: get_temperature_from_conserved => null()
    procedure(get_var_subr), pointer,nopass :: get_temperature_equi => null()
    procedure(get_2var_subr), pointer,nopass :: get_ne_nH => null()
    procedure(get_var_subr), pointer,nopass :: get_var_Rfactor => null()
    !> EoS snapshots + inverse accessor (set in bind_eos_to_source); let TC reach
    !> thermodynamics only through this object, never mod_eos directly.
    double precision :: inv_gamma_minus_1
    double precision :: nH2rhoFactor
    double precision :: log_T_floor
    procedure(eos_scalar2_func), pointer, nopass :: eint_from_T => null()
  end type tc_fluid

  public :: tc_get_mhd_params
  public :: tc_get_hd_params
  public :: get_tc_dt_mhd
  public :: get_tc_dt_hd
  public :: sts_set_source_tc_mhd
  public :: sts_set_source_tc_hd
  public :: tc_patch_negative_eint

contains

  subroutine tc_init_params(phys_gamma)
    use mod_global_parameters
    double precision, intent(in) :: phys_gamma

    tc_gamma_1=phys_gamma-1d0
  end subroutine tc_init_params

  !> Init TC coefficients: MHD case
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
      if(mype .eq. 0) print*, "Constant thermal conduction coefficients with values: ",fl%tc_k_para,fl%tc_k_perp
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

    if(fl%tc_k_para==0.d0) then
      if(SI_unit) then
        ! Spitzer thermal conductivity with SI units
        fl%tc_k_para=8.d-12*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      else
        ! Spitzer thermal conductivity with cgs units
        fl%tc_k_para=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      end if
      if(mype .eq. 0) print*, "Spitzer HD par: ",fl%tc_k_para
    else
      fl%tc_constant=.true.
      if(mype .eq. 0) print*, "Constant thermal conduction coefficient at value: ",fl%tc_k_para
    end if

  end subroutine tc_get_hd_params

  !> Get the explicit timestep for the TC (mhd implementation)
  !> Note: for multi-D MHD (1D MHD will use HD fall-back)
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

    double precision :: mf(ixO^S,1:ndir),Te(ixI^S),rho(ixI^S),gradT(ixI^S)
    double precision :: ne(ixI^S), nH_arr(ixI^S)
    double precision :: tmp(ixO^S),hfs(ixO^S),blocal(1:ndir),Bmag
    double precision :: dtdiff_tcond,maxtmp2
    integer          :: idims,ix^D


    ! B
    if(allocated(iw_mag)) then
      if(B0field) then
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          ^C&blocal(^C)=w({ix^D},iw_mag(^C))+block%B0({ix^D},^C,0)\
          Bmag=dsqrt(^C&blocal(^C)**2+)+smalldouble
          ^C&mf(ix^D,^C)=blocal(^C)/Bmag\
       {end do\}
      else
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          Bmag=dsqrt(^C&w(ix^D,iw_mag(^C))**2+)+smalldouble
          ^C&mf(ix^D,^C)=w(ix^D,iw_mag(^C))/Bmag\
       {end do\}
      end if
    else
      ! this if for ffhd or other physics modules without B components 
      mf(ixO^S,1:ndim)=block%B0(ixO^S,1:ndim,0)
    end if

    !temperature
    call fl%get_temperature_from_conserved(w,x,ixI^L,ixI^L,Te)
    call fl%get_rho(w,x,ixI^L,ixO^L,rho)
    call fl%get_ne_nH(ixI^L, ixO^L, w, x, ne, nH_arr)

    !tc_k_para_i
    if(fl%tc_constant) then
      tmp(ixO^S)=fl%tc_k_para
    else
      if(fl%tc_saturate) then
        ! Kannan 2016 MN 458, 410
        ! l_mfpe = 3^1.5*kB^2/(4*sqrt(pi)*e^4*lnLambda) * T^2/n_e
        if(SI_unit) then
          ! 5.730205638843984e27 = 3^1.5*kB_SI^2/(4*sqrt(pi)*e_SI^4*37)
          tmp(ixO^S)=Te(ixO^S)**2/ne(ixO^S)*5.730205638843984d27*unit_temperature**2/(unit_numberdensity*unit_length)
        else
          ! 7093.9239487765044 = 3^1.5*kB_cgs^2/(4*sqrt(pi)*e_cgs^4*37)
          tmp(ixO^S)=Te(ixO^S)**2/ne(ixO^S)*7093.9239487765044d0*unit_temperature**2/(unit_numberdensity*unit_length)
        end if
        do idims=1,ndim
          call gradient(Te,ixI^L,ixO^L,idims,gradT)
          if(idims==1) then
            hfs(ixO^S)=gradT(ixO^S)*mf(ixO^S,idims)
          else
            hfs(ixO^S)=hfs(ixO^S)+gradT(ixO^S)*mf(ixO^S,idims)
          end if
        end do
        ! kappa=kappa_Spitzer/(1+4.2*l_mfpe/(T/|gradT.b|))
        tmp(ixO^S)=fl%tc_k_para*Te(ixO^S)*Te(ixO^S)*dsqrt(Te(ixO^S))/(1.d0+4.2d0*tmp(ixO^S)*dabs(hfs(ixO^S))/Te(ixO^S))
      else
        ! kappa=kappa_Spitzer
        tmp(ixO^S)=fl%tc_k_para*Te(ixO^S)*Te(ixO^S)*dsqrt(Te(ixO^S))
      end if
    end if

    dtnew=bigdouble
    do idims=1,ndim
      ! approximate thermal conduction flux: tc_k_para_i/rho/dx*B_i**2/B**2
      maxtmp2=maxval(tmp(ixO^S)*mf(ixO^S,idims)**2/(rho(ixO^S)*block%ds(ixO^S,idims)**2))
      ! dt< dx_idim**2/((gamma-1)*tc_k_para_i/rho*B_i**2/B**2)
      dtdiff_tcond=1.d0/(tc_gamma_1*maxtmp2+smalldouble)
      ! limit the time step
      dtnew=min(dtnew,dtdiff_tcond)
    end do
    dtnew=dtnew/dble(ndim)
  end function get_tc_dt_mhd

  !> anisotropic thermal conduction with slope limited symmetric scheme
  !> Sharma 2007 Journal of Computational Physics 227, 123
  subroutine sts_set_source_tc_mhd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,fl)
    use mod_global_parameters
    use mod_fix_conserve
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    ! intent(inout) so tc_patch_negative_eint can repair w(:, ie) in place
    ! when tc_patch_eint = .true.; legacy path (patch off) leaves w untouched.
    double precision, intent(inout) ::   w(ixI^S,1:nw)
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
    if (fl%tc_patch_eint) call tc_patch_negative_eint(w, x, ixI^L, ixI^L, Te, fl%e_, fl)
    call fl%get_rho(w, x, ixI^L, ixI^L, rho)  !calculate rho in whole domain (+ghosts)
    if(slab_uniform) then
      call set_source_tc_mhd(ixI^L,ixO^L,w,x,fl,qvec,rho,Te,alpha)
    else
      call set_source_tc_mhd_geo(ixI^L,ixO^L,w,x,fl,qvec,rho,Te,alpha)
    end if
    if(fl%subtract_equi) then
      allocate(qvec_equi(ixI^S,1:ndim))
      call fl%get_temperature_equi(w, x, ixI^L, ixI^L, Te)  !calculate Te in whole domain (+ghosts)
      call fl%get_rho_equi(w, x, ixI^L, ixI^L, rho)  !calculate rho in whole domain (+ghosts)
      if(slab_uniform) then
        call set_source_tc_mhd(ixI^L,ixO^L,w,x,fl,qvec_equi,rho,Te,alpha)
      else
        call set_source_tc_mhd_geo(ixI^L,ixO^L,w,x,fl,qvec_equi,rho,Te,alpha)
      end if
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
    double precision :: minq,maxq,qd(ixI^S,2**(ndim-1)), blocal(ndir)
    integer :: idims,idir,ix^D,ix^L,ixC^L,ixA^L,ixB^L

    ix^L=ixO^L^LADD1;

    ! T gradient at cell faces
    ! b unit vector mf: magnetic field direction vector
    if(allocated(iw_mag)) then
      if(B0field) then
       {do ix^DB=ixmin^DB,ixmax^DB\}
          ^C&blocal(^C)=w({ix^D},iw_mag(^C))+block%B0({ix^D},^C,0)\
         {^IFTWOD
          if(blocal(1)/=0.d0) then
            mf(ix^D,1)=sign(1.d0,blocal(1))/dsqrt(1.d0+(^CE&(blocal(^CE)/blocal(1))**2+))
          else
            mf(ix^D,1)=0.d0
          end if
          if(blocal(2)/=0.d0) then
            mf(ix^D,2)=sign(1.d0,blocal(2))/dsqrt(1.d0+(^CF&(blocal(^CF)/blocal(2))**2+))
          else
            mf(ix^D,2)=0.d0
          end if
         }
         {^IFTHREED
          if(blocal(1)/=0.d0) then
            mf(ix^D,1)=sign(1.d0,blocal(1))/dsqrt(1.d0+(blocal(2)/blocal(1))**2+(blocal(3)/blocal(1))**2)
          else
            mf(ix^D,1)=0.d0
          end if
          if(blocal(2)/=0.d0) then
            mf(ix^D,2)=sign(1.d0,blocal(2))/dsqrt(1.d0+(blocal(1)/blocal(2))**2+(blocal(3)/blocal(2))**2)
          else
            mf(ix^D,2)=0.d0
          end if
          if(blocal(3)/=0.d0) then
            mf(ix^D,3)=sign(1.d0,blocal(3))/dsqrt(1.d0+(blocal(1)/blocal(3))**2+(blocal(2)/blocal(3))**2)
          else
            mf(ix^D,3)=0.d0
          end if
         }
       {end do\}
      else
       {do ix^DB=ixmin^DB,ixmax^DB\}
         {^IFTWOD
          if(w(ix^D,iw_mag(1))/=0.d0) then
            mf(ix^D,1)=sign(1.d0,w(ix^D,iw_mag(1)))/dsqrt(1.d0+(^CE&(w(ix^D,iw_mag(^CE))/w(ix^D,iw_mag(1)))**2+))
          else
            mf(ix^D,1)=0.d0
          end if
          if(w(ix^D,iw_mag(2))/=0.d0) then
            mf(ix^D,2)=sign(1.d0,w(ix^D,iw_mag(2)))/dsqrt(1.d0+(^CF&(w(ix^D,iw_mag(^CF))/w(ix^D,iw_mag(2)))**2+))
          else
            mf(ix^D,2)=0.d0
          end if
         }
         {^IFTHREED
          if(w(ix^D,iw_mag(1))/=0.d0) then
            mf(ix^D,1)=sign(1.d0,w(ix^D,iw_mag(1)))/dsqrt(1.d0+(w(ix^D,iw_mag(2))/w(ix^D,iw_mag(1)))**2+&
              (w(ix^D,iw_mag(3))/w(ix^D,iw_mag(1)))**2)
          else
            mf(ix^D,1)=0.d0
          end if
          if(w(ix^D,iw_mag(2))/=0.d0) then
            mf(ix^D,2)=sign(1.d0,w(ix^D,iw_mag(2)))/dsqrt(1.d0+(w(ix^D,iw_mag(1))/w(ix^D,iw_mag(2)))**2+&
              (w(ix^D,iw_mag(3))/w(ix^D,iw_mag(2)))**2)
          else
            mf(ix^D,2)=0.d0
          end if
          if(w(ix^D,iw_mag(3))/=0.d0) then
            mf(ix^D,3)=sign(1.d0,w(ix^D,iw_mag(3)))/dsqrt(1.d0+(w(ix^D,iw_mag(1))/w(ix^D,iw_mag(3)))**2+&
              (w(ix^D,iw_mag(2))/w(ix^D,iw_mag(3)))**2)
          else
            mf(ix^D,3)=0.d0
          end if
         }
       {end do\}
      end if
    else
      mf(ix^S,1:ndim)=block%B0(ix^S,1:ndim,0)
    endif
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
      call gradientF(Te,x,ixI^L,ixB^L,idims,gradT(ixI^S,idims))
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
            qdd(ix^D)=fl%tc_k_para*dsqrt(block%wextra(ix^D,fl%Tcoff_)**5)
          else
            qdd(ix^D)=fl%tc_k_para*dsqrt(Te(ix^D)**5)
          end if
       {end do\}
      else
        qdd(ix^S)=fl%tc_k_para*dsqrt(Te(ix^S)**5)
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
        if(B0field) then
          qdd(ix^S)=fl%tc_k_perp*rho(ix^S)**2/((^C&(w(ix^S,iw_mag(^C))+block%B0(ix^S,^C,0))**2+)*dsqrt(Te(ix^S))+smalldouble)
        else
          qdd(ix^S)=fl%tc_k_perp*rho(ix^S)**2/((^C&w(ix^S,iw_mag(^C))**2+)*dsqrt(Te(ix^S))+smalldouble)
        end if
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
    ! calculate thermal conduction flux with slope-limited symmetric scheme
    do idims=1,ndim
      ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
     {^IFTHREED
      if(idims==1) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          ! averaged b at face centers
          ^D&bcf({ix^D},^D)=0.25d0*(Bc({ix^D},^D)+Bc(ix1,ix2-1,ix3,^D)&
                         +Bc(ix1,ix2,ix3-1,^D)+Bc(ix1,ix2-1,ix3-1,^D))\
          kaf(ix^D)=0.25d0*(ka(ix1,ix2,ix3)+ka(ix1,ix2-1,ix3)&
                         +ka(ix1,ix2,ix3-1)+ka(ix1,ix2-1,ix3-1))
          ! averaged thermal conductivity at face centers
          if(fl%tc_perpendicular) &
          kef(ix^D)=0.25d0*(ke(ix1,ix2,ix3)+ke(ix1,ix2-1,ix3)&
                         +ke(ix1,ix2,ix3-1)+ke(ix1,ix2-1,ix3-1))
       {end do\}
      else if(idims==2) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          ^D&bcf({ix^D},^D)=0.25d0*(Bc({ix^D},^D)+Bc(ix1-1,ix2,ix3,^D)&
                         +Bc(ix1,ix2,ix3-1,^D)+Bc(ix1-1,ix2,ix3-1,^D))\
          kaf(ix^D)=0.25d0*(ka(ix1,ix2,ix3)+ka(ix1-1,ix2,ix3)&
                         +ka(ix1,ix2,ix3-1)+ka(ix1-1,ix2,ix3-1))
          if(fl%tc_perpendicular) &
          kef(ix^D)=0.25d0*(ke(ix1,ix2,ix3)+ke(ix1-1,ix2,ix3)&
                         +ke(ix1,ix2,ix3-1)+ke(ix1-1,ix2,ix3-1))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          ^D&bcf({ix^D},^D)=0.25d0*(Bc({ix^D},^D)+Bc(ix1,ix2-1,ix3,^D)&
                         +Bc(ix1-1,ix2,ix3,^D)+Bc(ix1-1,ix2-1,ix3,^D))\
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
          ^D&bcf({ix^D},^D)=0.5d0*(Bc(ix1,ix2,^D)+Bc(ix1,ix2-1,^D))\
          kaf(ix^D)=0.5d0*(ka(ix1,ix2)+ka(ix1,ix2-1))
          if(fl%tc_perpendicular) &
          kef(ix^D)=0.5d0*(ke(ix1,ix2)+ke(ix1,ix2-1))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          ^D&bcf({ix^D},^D)=0.5d0*(Bc(ix1,ix2,^D)+Bc(ix1-1,ix2,^D))\
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
            qd(ix^D,1)=minq
          else if(qdd(ix^D)>maxq) then
            qd(ix^D,1)=maxq
          else
            qd(ix^D,1)=qdd(ix^D)
          end if
          if(qdd(ix1,ix2-1,ix3)<minq) then
            qd(ix^D,2)=minq
          else if(qdd(ix1,ix2-1,ix3)>maxq) then
            qd(ix^D,2)=maxq
          else
            qd(ix^D,2)=qdd(ix1,ix2-1,ix3)
          end if
          if(qdd(ix1,ix2,ix3-1)<minq) then
            qd(ix^D,3)=minq
          else if(qdd(ix1,ix2,ix3-1)>maxq) then
            qd(ix^D,3)=maxq
          else
            qd(ix^D,3)=qdd(ix1,ix2,ix3-1)
          end if
          if(qdd(ix1,ix2-1,ix3-1)<minq) then
            qd(ix^D,4)=minq
          else if(qdd(ix1,ix2-1,ix3-1)>maxq) then
            qd(ix^D,4)=maxq
          else
            qd(ix^D,4)=qdd(ix1,ix2-1,ix3-1)
          end if
          qvec(ix^D,idims)=kaf(ix^D)*0.25d0*(Bc(ix^D,idims)**2*qd(ix^D,1)+Bc(ix1,ix2-1,ix3,idims)**2*qd(ix^D,2)&
                         +Bc(ix1,ix2,ix3-1,idims)**2*qd(ix^D,3)+Bc(ix1,ix2-1,ix3-1,idims)**2*qd(ix^D,4))
          if(fl%tc_perpendicular) &
          qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.25d0*(qd(ix^D,1)+qd(ix^D,2)+qd(ix^D,3)+qd(ix^D,4))
       {end do\}
      else if(idims==2) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          if(qdd(ix^D)<minq) then
            qd(ix^D,1)=minq
          else if(qdd(ix^D)>maxq) then
            qd(ix^D,1)=maxq
          else
            qd(ix^D,1)=qdd(ix^D)
          end if
          if(qdd(ix1-1,ix2,ix3)<minq) then
            qd(ix^D,2)=minq
          else if(qdd(ix1-1,ix2,ix3)>maxq) then
            qd(ix^D,2)=maxq
          else
            qd(ix^D,2)=qdd(ix1-1,ix2,ix3)
          end if
          if(qdd(ix1,ix2,ix3-1)<minq) then
            qd(ix^D,3)=minq
          else if(qdd(ix1,ix2,ix3-1)>maxq) then
            qd(ix^D,3)=maxq
          else
            qd(ix^D,3)=qdd(ix1,ix2,ix3-1)
          end if
          if(qdd(ix1-1,ix2,ix3-1)<minq) then
            qd(ix^D,4)=minq
          else if(qdd(ix1-1,ix2,ix3-1)>maxq) then
            qd(ix^D,4)=maxq
          else
            qd(ix^D,4)=qdd(ix1-1,ix2,ix3-1)
          end if
          qvec(ix^D,idims)=kaf(ix^D)*0.25d0*(Bc(ix^D,idims)**2*qd(ix^D,1)+Bc(ix1-1,ix2,ix3,idims)**2*qd(ix^D,2)&
                         +Bc(ix1,ix2,ix3-1,idims)**2*qd(ix^D,3)+Bc(ix1-1,ix2,ix3-1,idims)**2*qd(ix^D,4))
          if(fl%tc_perpendicular) &
          qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.25d0*(qd(ix^D,1)+qd(ix^D,2)+qd(ix^D,3)+qd(ix^D,4))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          if(qdd(ix^D)<minq) then
            qd(ix^D,1)=minq
          else if(qdd(ix^D)>maxq) then
            qd(ix^D,1)=maxq
          else
            qd(ix^D,1)=qdd(ix^D)
          end if
          if(qdd(ix1-1,ix2,ix3)<minq) then
            qd(ix^D,2)=minq
          else if(qdd(ix1-1,ix2,ix3)>maxq) then
            qd(ix^D,2)=maxq
          else
            qd(ix^D,2)=qdd(ix1-1,ix2,ix3)
          end if
          if(qdd(ix1,ix2-1,ix3)<minq) then
            qd(ix^D,3)=minq
          else if(qdd(ix1,ix2-1,ix3)>maxq) then
            qd(ix^D,3)=maxq
          else
            qd(ix^D,3)=qdd(ix1,ix2-1,ix3)
          end if
          if(qdd(ix1-1,ix2-1,ix3)<minq) then
            qd(ix^D,4)=minq
          else if(qdd(ix1-1,ix2-1,ix3)>maxq) then
            qd(ix^D,4)=maxq
          else
            qd(ix^D,4)=qdd(ix1-1,ix2-1,ix3)
          end if
          qvec(ix^D,idims)=kaf(ix^D)*0.25d0*(Bc(ix^D,idims)**2*qd(ix^D,1)+Bc(ix1-1,ix2,ix3,idims)**2*qd(ix^D,2)&
                         +Bc(ix1,ix2-1,ix3,idims)**2*qd(ix^D,3)+Bc(ix1-1,ix2-1,ix3,idims)**2*qd(ix^D,4))
          if(fl%tc_perpendicular) &
          qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.25d0*(qd(ix^D,1)+qd(ix^D,2)+qd(ix^D,3)+qd(ix^D,4))
       {end do\}
      end if
     }
     {^IFTWOD
      if(idims==1) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          if(qdd(ix^D)<minq) then
            qd(ix^D,1)=minq
          else if(qdd(ix^D)>maxq) then
            qd(ix^D,1)=maxq
          else
            qd(ix^D,1)=qdd(ix^D)
          end if
          if(qdd(ix1,ix2-1)<minq) then
            qd(ix^D,2)=minq
          else if(qdd(ix1,ix2-1)>maxq) then
            qd(ix^D,2)=maxq
          else
            qd(ix^D,2)=qdd(ix1,ix2-1)
          end if
          qvec(ix^D,idims)=kaf(ix^D)*0.5d0*(Bc(ix1,ix2,idims)**2*qd(ix^D,1)+Bc(ix1,ix2-1,idims)**2*qd(ix^D,2))
          if(fl%tc_perpendicular) &
          qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.5d0*(qd(ix^D,1)+qd(ix^D,2))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          if(qdd(ix^D)<minq) then
            qd(ix^D,1)=minq
          else if(qdd(ix^D)>maxq) then
            qd(ix^D,1)=maxq
          else
            qd(ix^D,1)=qdd(ix^D)
          end if
          if(qdd(ix1-1,ix2)<minq) then
            qd(ix^D,2)=minq
          else if(qdd(ix1-1,ix2)>maxq) then
            qd(ix^D,2)=maxq
          else
            qd(ix^D,2)=qdd(ix1-1,ix2)
          end if
          qvec(ix^D,idims)=kaf(ix^D)*0.5d0*(Bc(ix1,ix2,idims)**2*qd(ix^D,1)+Bc(ix1-1,ix2,idims)**2*qd(ix^D,2))
          if(fl%tc_perpendicular) &
          qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.5d0*(qd(ix^D,1)+qd(ix^D,2))
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
        ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135: phi=1.1, Balbus and Mckee 1982 ApJ, 252, 529: phi=0.3)
        ! unsigned saturated TC flux = 5 phi rho c**3, c=sqrt(p/rho) is isothermal sound speed, phi=0.3
        ixB^L=ixA^L+kr(idims,^D);
        qdd(ixA^S)=0.75d0*(rho(ixA^S)+rho(ixB^S))*dsqrt(0.5d0*(Te(ixA^S)+Te(ixB^S)))**3*dabs(Bnorm(ixA^S))
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          if(dabs(qvec(ix^D,idims))>qdd(ix^D)) then
            qvec(ix^D,idims)=sign(1.d0,qvec(ix^D,idims))*qdd(ix^D)
          end if
       {end do\}
      end if
    end do
  end subroutine set_source_tc_mhd
  subroutine set_source_tc_mhd_geo(ixI^L,ixO^L,w,x,fl,qvec,rho,Te,alpha)
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
    double precision :: minq,maxq,qd(ixI^S,2**(ndim-1)), blocal(ndir)
    integer :: idims,idir,ix^D,ix^L,ixC^L,ixA^L,ixB^L

    ix^L=ixO^L^LADD1;

    ! T gradient at cell faces
    ! b unit vector mf: magnetic field direction vector
    if(allocated(iw_mag)) then
      if(B0field) then
       {do ix^DB=ixmin^DB,ixmax^DB\}
          ^C&blocal(^C)=w({ix^D},iw_mag(^C))+block%B0({ix^D},^C,0)\
         {^IFTWOD
          if(blocal(1)/=0.d0) then
            mf(ix^D,1)=sign(1.d0,blocal(1))/dsqrt(1.d0+(^CE&(blocal(^CE)/blocal(1))**2+))
          else
            mf(ix^D,1)=0.d0
          end if
          if(blocal(2)/=0.d0) then
            mf(ix^D,2)=sign(1.d0,blocal(2))/dsqrt(1.d0+(^CF&(blocal(^CF)/blocal(2))**2+))
          else
            mf(ix^D,2)=0.d0
          end if
         }
         {^IFTHREED
          if(blocal(1)/=0.d0) then
            mf(ix^D,1)=sign(1.d0,blocal(1))/dsqrt(1.d0+(blocal(2)/blocal(1))**2+(blocal(3)/blocal(1))**2)
          else
            mf(ix^D,1)=0.d0
          end if
          if(blocal(2)/=0.d0) then
            mf(ix^D,2)=sign(1.d0,blocal(2))/dsqrt(1.d0+(blocal(1)/blocal(2))**2+(blocal(3)/blocal(2))**2)
          else
            mf(ix^D,2)=0.d0
          end if
          if(blocal(3)/=0.d0) then
            mf(ix^D,3)=sign(1.d0,blocal(3))/dsqrt(1.d0+(blocal(1)/blocal(3))**2+(blocal(2)/blocal(3))**2)
          else
            mf(ix^D,3)=0.d0
          end if
         }
       {end do\}
      else
       {do ix^DB=ixmin^DB,ixmax^DB\}
         {^IFTWOD
          if(w(ix^D,iw_mag(1))/=0.d0) then
            mf(ix^D,1)=sign(1.d0,w(ix^D,iw_mag(1)))/dsqrt(1.d0+(^CE&(w(ix^D,iw_mag(^CE))/w(ix^D,iw_mag(1)))**2+))
          else
            mf(ix^D,1)=0.d0
          end if
          if(w(ix^D,iw_mag(2))/=0.d0) then
            mf(ix^D,2)=sign(1.d0,w(ix^D,iw_mag(2)))/dsqrt(1.d0+(^CF&(w(ix^D,iw_mag(^CF))/w(ix^D,iw_mag(2)))**2+))
          else
            mf(ix^D,2)=0.d0
          end if
         }
         {^IFTHREED
          if(w(ix^D,iw_mag(1))/=0.d0) then
            mf(ix^D,1)=sign(1.d0,w(ix^D,iw_mag(1)))/dsqrt(1.d0+(w(ix^D,iw_mag(2))/w(ix^D,iw_mag(1)))**2+&
              (w(ix^D,iw_mag(3))/w(ix^D,iw_mag(1)))**2)
          else
            mf(ix^D,1)=0.d0
          end if
          if(w(ix^D,iw_mag(2))/=0.d0) then
            mf(ix^D,2)=sign(1.d0,w(ix^D,iw_mag(2)))/dsqrt(1.d0+(w(ix^D,iw_mag(1))/w(ix^D,iw_mag(2)))**2+&
              (w(ix^D,iw_mag(3))/w(ix^D,iw_mag(2)))**2)
          else
            mf(ix^D,2)=0.d0
          end if
          if(w(ix^D,iw_mag(3))/=0.d0) then
            mf(ix^D,3)=sign(1.d0,w(ix^D,iw_mag(3)))/dsqrt(1.d0+(w(ix^D,iw_mag(1))/w(ix^D,iw_mag(3)))**2+&
              (w(ix^D,iw_mag(2))/w(ix^D,iw_mag(3)))**2)
          else
            mf(ix^D,3)=0.d0
          end if
         }
       {end do\}
      end if
    else
      mf(ix^S,1:ndim)=block%B0(ix^S,1:ndim,0)
    endif
    ! ixC is cell-corner index
    ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;
    ! b unit vector at cell corner
   {^IFTHREED
    do idims=1,3
   {do ix^DB=ixCmin^DB,ixCmax^DB\}
      Bc(ix^D,idims)=(mf(ix1,ix2,ix3,idims)*block%dvolume(ix1,ix2,ix3)+mf(ix1+1,ix2,ix3,idims)*block%dvolume(ix1+1,ix2,ix3)&
                     +mf(ix1,ix2+1,ix3,idims)*block%dvolume(ix1,ix2+1,ix3)+mf(ix1+1,ix2+1,ix3,idims)*block%dvolume(ix1+1,ix2+1,ix3)&
                     +mf(ix1,ix2,ix3+1,idims)*block%dvolume(ix1,ix2,ix3+1)+mf(ix1+1,ix2,ix3+1,idims)*block%dvolume(ix1+1,ix2,ix3+1)&
                     +mf(ix1,ix2+1,ix3+1,idims)*block%dvolume(ix1,ix2+1,ix3+1)+mf(ix1+1,ix2+1,ix3+1,idims)*block%dvolume(ix1+1,ix2+1,ix3+1))&
            /(block%dvolume(ix1,ix2,ix3)+block%dvolume(ix1+1,ix2,ix3)+block%dvolume(ix1,ix2+1,ix3)+block%dvolume(ix1+1,ix2+1,ix3)&
             +block%dvolume(ix1,ix2,ix3+1)+block%dvolume(ix1+1,ix2,ix3+1)+block%dvolume(ix1,ix2+1,ix3+1)+block%dvolume(ix1+1,ix2+1,ix3+1))
   {end do\}
    end do
   }
   {^IFTWOD
    do idims=1,2
   {do ix^DB=ixCmin^DB,ixCmax^DB\}
      Bc(ix^D,idims)=(mf(ix1,ix2,idims)*block%dvolume(ix1,ix2)+mf(ix1+1,ix2,idims)*block%dvolume(ix1+1,ix2)&
                     +mf(ix1,ix2+1,idims)*block%dvolume(ix1,ix2+1)+mf(ix1+1,ix2+1,idims)*block%dvolume(ix1+1,ix2+1))&
                 /(block%dvolume(ix1,ix2)+block%dvolume(ix1+1,ix2)+block%dvolume(ix1,ix2+1)+block%dvolume(ix1+1,ix2+1))
   {end do\}
    end do
   }
    ! T gradient at cell faces
    do idims=1,ndim
      ixBmin^D=ixmin^D;
      ixBmax^D=ixmax^D-kr(idims,^D);
      call gradientF(Te,x,ixI^L,ixB^L,idims,gradT(ixI^S,idims))
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
            qdd(ix^D)=fl%tc_k_para*dsqrt(block%wextra(ix^D,fl%Tcoff_)**5)
          else
            qdd(ix^D)=fl%tc_k_para*dsqrt(Te(ix^D)**5)
          end if
       {end do\}
      else
        qdd(ix^S)=fl%tc_k_para*dsqrt(Te(ix^S)**5)
      end if
     ! cell corner parallel conductivity in ka
     {^IFTHREED
     {do ix^DB=ixCmin^DB,ixCmax^DB\}
        ka(ix^D)=(qdd(ix1,ix2,ix3)*block%dvolume(ix1,ix2,ix3)+qdd(ix1+1,ix2,ix3)*block%dvolume(ix1+1,ix2,ix3)&
           +qdd(ix1,ix2+1,ix3)*block%dvolume(ix1,ix2+1,ix3)+qdd(ix1+1,ix2+1,ix3)*block%dvolume(ix1+1,ix2+1,ix3)&
           +qdd(ix1,ix2,ix3+1)*block%dvolume(ix1,ix2,ix3+1)+qdd(ix1+1,ix2,ix3+1)*block%dvolume(ix1+1,ix2,ix3+1)&
           +qdd(ix1,ix2+1,ix3+1)*block%dvolume(ix1,ix2+1,ix3+1)+qdd(ix1+1,ix2+1,ix3+1)*block%dvolume(ix1+1,ix2+1,ix3+1))&
          /(block%dvolume(ix1,ix2,ix3)+block%dvolume(ix1+1,ix2,ix3)+block%dvolume(ix1,ix2+1,ix3)+block%dvolume(ix1+1,ix2+1,ix3)&
           +block%dvolume(ix1,ix2,ix3+1)+block%dvolume(ix1+1,ix2,ix3+1)+block%dvolume(ix1,ix2+1,ix3+1)+block%dvolume(ix1+1,ix2+1,ix3+1))
     {end do\}
     }
     {^IFTWOD
     {do ix^DB=ixCmin^DB,ixCmax^DB\}
        ka(ix^D)=(qdd(ix1,ix2)*block%dvolume(ix1,ix2)+qdd(ix1+1,ix2)*block%dvolume(ix1+1,ix2)&
           +qdd(ix1,ix2+1)*block%dvolume(ix1,ix2+1)+qdd(ix1+1,ix2+1)*block%dvolume(ix1+1,ix2+1))&
         /(block%dvolume(ix1,ix2)+block%dvolume(ix1+1,ix2)+block%dvolume(ix1,ix2+1)+block%dvolume(ix1+1,ix2+1))
     {end do\}
     }
      ! compensate with perpendicular conductivity
      if(fl%tc_perpendicular) then
        if(B0field) then
          qdd(ix^S)=fl%tc_k_perp*rho(ix^S)**2/((^C&(w(ix^S,iw_mag(^C))+block%B0(ix^S,^C,0))**2+)*dsqrt(Te(ix^S))+smalldouble)
        else
          qdd(ix^S)=fl%tc_k_perp*rho(ix^S)**2/((^C&w(ix^S,iw_mag(^C))**2+)*dsqrt(Te(ix^S))+smalldouble)
        end if
       {^IFTHREED
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          ke(ix^D)=(qdd(ix1,ix2,ix3)*block%dvolume(ix1,ix2,ix3)+qdd(ix1+1,ix2,ix3)*block%dvolume(ix1+1,ix2,ix3)&
             +qdd(ix1,ix2+1,ix3)*block%dvolume(ix1,ix2+1,ix3)+qdd(ix1+1,ix2+1,ix3)*block%dvolume(ix1+1,ix2+1,ix3)&
             +qdd(ix1,ix2,ix3+1)*block%dvolume(ix1,ix2,ix3+1)+qdd(ix1+1,ix2,ix3+1)*block%dvolume(ix1+1,ix2,ix3+1)&
             +qdd(ix1,ix2+1,ix3+1)*block%dvolume(ix1,ix2+1,ix3+1)+qdd(ix1+1,ix2+1,ix3+1)*block%dvolume(ix1+1,ix2+1,ix3+1))&
            /(block%dvolume(ix1,ix2,ix3)+block%dvolume(ix1+1,ix2,ix3)+block%dvolume(ix1,ix2+1,ix3)+block%dvolume(ix1+1,ix2+1,ix3)&
             +block%dvolume(ix1,ix2,ix3+1)+block%dvolume(ix1+1,ix2,ix3+1)+block%dvolume(ix1,ix2+1,ix3+1)+block%dvolume(ix1+1,ix2+1,ix3+1))
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
          ke(ix^D)=(qdd(ix1,ix2)*block%dvolume(ix1,ix2)+qdd(ix1+1,ix2)*block%dvolume(ix1+1,ix2)&
             +qdd(ix1,ix2+1)*block%dvolume(ix1,ix2+1)+qdd(ix1+1,ix2+1)*block%dvolume(ix1+1,ix2+1))&
           /(block%dvolume(ix1,ix2)+block%dvolume(ix1+1,ix2)+block%dvolume(ix1,ix2+1)+block%dvolume(ix1+1,ix2+1))
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
    ! calculate thermal conduction flux with slope-limited symmetric scheme
    do idims=1,ndim
      ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
     {^IFTHREED
      if(idims==1) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          ! averaged b at face centers
          ^D&bcf({ix^D},^D)=0.25d0*(Bc({ix^D},^D)+Bc(ix1,ix2-1,ix3,^D)&
                         +Bc(ix1,ix2,ix3-1,^D)+Bc(ix1,ix2-1,ix3-1,^D))\
          kaf(ix^D)=0.25d0*(ka(ix1,ix2,ix3)+ka(ix1,ix2-1,ix3)&
                         +ka(ix1,ix2,ix3-1)+ka(ix1,ix2-1,ix3-1))
          ! averaged thermal conductivity at face centers
          if(fl%tc_perpendicular) &
          kef(ix^D)=0.25d0*(ke(ix1,ix2,ix3)+ke(ix1,ix2-1,ix3)&
                         +ke(ix1,ix2,ix3-1)+ke(ix1,ix2-1,ix3-1))
       {end do\}
      else if(idims==2) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          ^D&bcf({ix^D},^D)=0.25d0*(Bc({ix^D},^D)+Bc(ix1-1,ix2,ix3,^D)&
                         +Bc(ix1,ix2,ix3-1,^D)+Bc(ix1-1,ix2,ix3-1,^D))\
          kaf(ix^D)=0.25d0*(ka(ix1,ix2,ix3)+ka(ix1-1,ix2,ix3)&
                         +ka(ix1,ix2,ix3-1)+ka(ix1-1,ix2,ix3-1))
          if(fl%tc_perpendicular) &
          kef(ix^D)=0.25d0*(ke(ix1,ix2,ix3)+ke(ix1-1,ix2,ix3)&
                         +ke(ix1,ix2,ix3-1)+ke(ix1-1,ix2,ix3-1))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          ^D&bcf({ix^D},^D)=0.25d0*(Bc({ix^D},^D)+Bc(ix1,ix2-1,ix3,^D)&
                         +Bc(ix1-1,ix2,ix3,^D)+Bc(ix1-1,ix2-1,ix3,^D))\
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
          ^D&bcf({ix^D},^D)=0.5d0*(Bc(ix1,ix2,^D)+Bc(ix1,ix2-1,^D))\
          kaf(ix^D)=0.5d0*(ka(ix1,ix2)+ka(ix1,ix2-1))
          if(fl%tc_perpendicular) &
          kef(ix^D)=0.5d0*(ke(ix1,ix2)+ke(ix1,ix2-1))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          ^D&bcf({ix^D},^D)=0.5d0*(Bc(ix1,ix2,^D)+Bc(ix1-1,ix2,^D))\
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
          qdd(ix^D)=(gradT(ix1,ix2,ix3,idims)*block%surfaceC(ix1,ix2,ix3,1)&
                  +gradT(ix1,ix2+1,ix3,idims)*block%surfaceC(ix1,ix2+1,ix3,1)&
                  +gradT(ix1,ix2,ix3+1,idims)*block%surfaceC(ix1,ix2,ix3+1,1)&
                +gradT(ix1,ix2+1,ix3+1,idims)*block%surfaceC(ix1,ix2+1,ix3+1,1))/&
              (block%surfaceC(ix1,ix2,ix3,1)+block%surfaceC(ix1,ix2+1,ix3,1)&
              +block%surfaceC(ix1,ix2,ix3+1,1)+block%surfaceC(ix1,ix2+1,ix3+1,1)+smalldouble**2)
       {end do\}
      else if(idims==2) then
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qdd(ix^D)=(gradT(ix1,ix2,ix3,idims)*block%surfaceC(ix1,ix2,ix3,2)&
                  +gradT(ix1+1,ix2,ix3,idims)*block%surfaceC(ix1+1,ix2,ix3,2)&
                  +gradT(ix1,ix2,ix3+1,idims)*block%surfaceC(ix1,ix2,ix3+1,2)&
                +gradT(ix1+1,ix2,ix3+1,idims)*block%surfaceC(ix1+1,ix2,ix3+1,2))/&
            (block%surfaceC(ix1,ix2,ix3,2)+block%surfaceC(ix1+1,ix2,ix3,2)&
          +block%surfaceC(ix1,ix2,ix3+1,2)+block%surfaceC(ix1+1,ix2,ix3+1,2)+smalldouble**2)
       {end do\}
      else
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qdd(ix^D)=(gradT(ix1,ix2,ix3,idims)*block%surfaceC(ix1,ix2,ix3,3)&
                  +gradT(ix1+1,ix2,ix3,idims)*block%surfaceC(ix1+1,ix2,ix3,3)&
                  +gradT(ix1,ix2+1,ix3,idims)*block%surfaceC(ix1,ix2+1,ix3,3)&
                +gradT(ix1+1,ix2+1,ix3,idims)*block%surfaceC(ix1+1,ix2+1,ix3,3))/&
               (block%surfaceC(ix1,ix2,ix3,3)+block%surfaceC(ix1+1,ix2,ix3,3)&
             +block%surfaceC(ix1,ix2+1,ix3,3)+block%surfaceC(ix1+1,ix2+1,ix3,3))
       {end do\}
      end if
     }
     {^IFTWOD
      if(idims==1) then
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qdd(ix^D)=(gradT(ix1,ix2,idims)*block%surfaceC(ix1,ix2,1)&
                  +gradT(ix1,ix2+1,idims)*block%surfaceC(ix1,ix2+1,1))/&
               (block%surfaceC(ix1,ix2,1)+block%surfaceC(ix1,ix2+1,1))
       {end do\}
      else
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qdd(ix^D)=(gradT(ix1,ix2,idims)*block%surfaceC(ix1,ix2,2)&
                  +gradT(ix1+1,ix2,idims)*block%surfaceC(ix1+1,ix2,2))/&
               (block%surfaceC(ix1,ix2,2)+block%surfaceC(ix1+1,ix2,2)+smalldouble)
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
            qd(ix^D,1)=minq
          else if(qdd(ix^D)>maxq) then
            qd(ix^D,1)=maxq
          else
            qd(ix^D,1)=qdd(ix^D)
          end if
          if(qdd(ix1,ix2-1,ix3)<minq) then
            qd(ix^D,2)=minq
          else if(qdd(ix1,ix2-1,ix3)>maxq) then
            qd(ix^D,2)=maxq
          else
            qd(ix^D,2)=qdd(ix1,ix2-1,ix3)
          end if
          if(qdd(ix1,ix2,ix3-1)<minq) then
            qd(ix^D,3)=minq
          else if(qdd(ix1,ix2,ix3-1)>maxq) then
            qd(ix^D,3)=maxq
          else
            qd(ix^D,3)=qdd(ix1,ix2,ix3-1)
          end if
          if(qdd(ix1,ix2-1,ix3-1)<minq) then
            qd(ix^D,4)=minq
          else if(qdd(ix1,ix2-1,ix3-1)>maxq) then
            qd(ix^D,4)=maxq
          else
            qd(ix^D,4)=qdd(ix1,ix2-1,ix3-1)
          end if
          qvec(ix^D,idims)=kaf(ix^D)*0.25d0*(Bc(ix^D,idims)**2*qd(ix^D,1)+Bc(ix1,ix2-1,ix3,idims)**2*qd(ix^D,2)&
                         +Bc(ix1,ix2,ix3-1,idims)**2*qd(ix^D,3)+Bc(ix1,ix2-1,ix3-1,idims)**2*qd(ix^D,4))
          if(fl%tc_perpendicular) &
          qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.25d0*(qd(ix^D,1)+qd(ix^D,2)+qd(ix^D,3)+qd(ix^D,4))
       {end do\}
      else if(idims==2) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          if(qdd(ix^D)<minq) then
            qd(ix^D,1)=minq
          else if(qdd(ix^D)>maxq) then
            qd(ix^D,1)=maxq
          else
            qd(ix^D,1)=qdd(ix^D)
          end if
          if(qdd(ix1-1,ix2,ix3)<minq) then
            qd(ix^D,2)=minq
          else if(qdd(ix1-1,ix2,ix3)>maxq) then
            qd(ix^D,2)=maxq
          else
            qd(ix^D,2)=qdd(ix1-1,ix2,ix3)
          end if
          if(qdd(ix1,ix2,ix3-1)<minq) then
            qd(ix^D,3)=minq
          else if(qdd(ix1,ix2,ix3-1)>maxq) then
            qd(ix^D,3)=maxq
          else
            qd(ix^D,3)=qdd(ix1,ix2,ix3-1)
          end if
          if(qdd(ix1-1,ix2,ix3-1)<minq) then
            qd(ix^D,4)=minq
          else if(qdd(ix1-1,ix2,ix3-1)>maxq) then
            qd(ix^D,4)=maxq
          else
            qd(ix^D,4)=qdd(ix1-1,ix2,ix3-1)
          end if
          qvec(ix^D,idims)=kaf(ix^D)*0.25d0*(Bc(ix^D,idims)**2*qd(ix^D,1)+Bc(ix1-1,ix2,ix3,idims)**2*qd(ix^D,2)&
                         +Bc(ix1,ix2,ix3-1,idims)**2*qd(ix^D,3)+Bc(ix1-1,ix2,ix3-1,idims)**2*qd(ix^D,4))
          if(fl%tc_perpendicular) &
          qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.25d0*(qd(ix^D,1)+qd(ix^D,2)+qd(ix^D,3)+qd(ix^D,4))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          if(qdd(ix^D)<minq) then
            qd(ix^D,1)=minq
          else if(qdd(ix^D)>maxq) then
            qd(ix^D,1)=maxq
          else
            qd(ix^D,1)=qdd(ix^D)
          end if
          if(qdd(ix1-1,ix2,ix3)<minq) then
            qd(ix^D,2)=minq
          else if(qdd(ix1-1,ix2,ix3)>maxq) then
            qd(ix^D,2)=maxq
          else
            qd(ix^D,2)=qdd(ix1-1,ix2,ix3)
          end if
          if(qdd(ix1,ix2-1,ix3)<minq) then
            qd(ix^D,3)=minq
          else if(qdd(ix1,ix2-1,ix3)>maxq) then
            qd(ix^D,3)=maxq
          else
            qd(ix^D,3)=qdd(ix1,ix2-1,ix3)
          end if
          if(qdd(ix1-1,ix2-1,ix3)<minq) then
            qd(ix^D,4)=minq
          else if(qdd(ix1-1,ix2-1,ix3)>maxq) then
            qd(ix^D,4)=maxq
          else
            qd(ix^D,4)=qdd(ix1-1,ix2-1,ix3)
          end if
          qvec(ix^D,idims)=kaf(ix^D)*0.25d0*(Bc(ix^D,idims)**2*qd(ix^D,1)+Bc(ix1-1,ix2,ix3,idims)**2*qd(ix^D,2)&
                         +Bc(ix1,ix2-1,ix3,idims)**2*qd(ix^D,3)+Bc(ix1-1,ix2-1,ix3,idims)**2*qd(ix^D,4))
          if(fl%tc_perpendicular) &
          qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.25d0*(qd(ix^D,1)+qd(ix^D,2)+qd(ix^D,3)+qd(ix^D,4))
       {end do\}
      end if
     }
     {^IFTWOD
      if(idims==1) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          if(qdd(ix^D)<minq) then
            qd(ix^D,1)=minq
          else if(qdd(ix^D)>maxq) then
            qd(ix^D,1)=maxq
          else
            qd(ix^D,1)=qdd(ix^D)
          end if
          if(qdd(ix1,ix2-1)<minq) then
            qd(ix^D,2)=minq
          else if(qdd(ix1,ix2-1)>maxq) then
            qd(ix^D,2)=maxq
          else
            qd(ix^D,2)=qdd(ix1,ix2-1)
          end if
          qvec(ix^D,idims)=kaf(ix^D)*0.5d0*(Bc(ix1,ix2,idims)**2*qd(ix^D,1)+Bc(ix1,ix2-1,idims)**2*qd(ix^D,2))
          if(fl%tc_perpendicular) &
          qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.5d0*(qd(ix^D,1)+qd(ix^D,2))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          minq=min(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          maxq=max(alpha*gradT(ix^D,idims),gradT(ix^D,idims)/alpha)
          if(qdd(ix^D)<minq) then
            qd(ix^D,1)=minq
          else if(qdd(ix^D)>maxq) then
            qd(ix^D,1)=maxq
          else
            qd(ix^D,1)=qdd(ix^D)
          end if
          if(qdd(ix1-1,ix2)<minq) then
            qd(ix^D,2)=minq
          else if(qdd(ix1-1,ix2)>maxq) then
            qd(ix^D,2)=maxq
          else
            qd(ix^D,2)=qdd(ix1-1,ix2)
          end if
          qvec(ix^D,idims)=kaf(ix^D)*0.5d0*(Bc(ix1,ix2,idims)**2*qd(ix^D,1)+Bc(ix1-1,ix2,idims)**2*qd(ix^D,2))
          if(fl%tc_perpendicular) &
          qvec(ix^D,idims)=qvec(ix^D,idims)+kef(ix^D)*0.5d0*(qd(ix^D,1)+qd(ix^D,2))
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
        ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135: phi=1.1, Balbus and Mckee 1982 ApJ, 252, 529: phi=0.3)
        ! unsigned saturated TC flux = 5 phi rho c**3, c=sqrt(p/rho) is isothermal sound speed, phi=0.3
        ixB^L=ixA^L+kr(idims,^D);
        qdd(ixA^S)=0.75d0*(rho(ixA^S)+rho(ixB^S))*dsqrt(0.5d0*(Te(ixA^S)+Te(ixB^S)))**3*dabs(Bnorm(ixA^S))
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          if(dabs(qvec(ix^D,idims))>qdd(ix^D)) then
            qvec(ix^D,idims)=sign(1.d0,qvec(ix^D,idims))*qdd(ix^D)
          end if
       {end do\}
      end if
    end do
  end subroutine set_source_tc_mhd_geo

  function slope_limiter(f,ixI^L,ixO^L,idims,pm,tc_slope_limiter) result(lf)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L, idims, pm
    double precision, intent(in) :: f(ixI^S)
    double precision :: lf(ixI^S)
    integer, intent(in)  :: tc_slope_limiter

    double precision, parameter :: qsmall=1.d-12
    double precision :: signf(ixI^S)
    integer :: ixB^L

    ixB^L=ixO^L+pm*kr(idims,^D);
    signf(ixO^S)=sign(1.d0,f(ixO^S))
    select case(tc_slope_limiter)
     case(1)
       ! 'MC' monotonized central limiter Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       lf(ixO^S)=two*signf(ixO^S)* &
            max(zero,min(dabs(f(ixO^S)),signf(ixO^S)*f(ixB^S),&
            signf(ixO^S)*quarter*(f(ixB^S)+f(ixO^S))))
     case(2)
       ! 'minmod' limiter
       lf(ixO^S)=signf(ixO^S)*max(0.d0,min(abs(f(ixO^S)),signf(ixO^S)*f(ixB^S)))
     case(3)
       ! 'superbee' Roe superbee limiter (eq.3.51i)
       lf(ixO^S)=signf(ixO^S)* &
            max(zero,min(two*dabs(f(ixO^S)),signf(ixO^S)*f(ixB^S)),&
            min(dabs(f(ixO^S)),two*signf(ixO^S)*f(ixB^S)))
     case(4)
       ! 'koren' Barry Koren Right variant
       lf(ixO^S)=signf(ixO^S)* &
            max(zero,min(two*dabs(f(ixO^S)),two*signf(ixO^S)*f(ixB^S),&
            (two*f(ixB^S)*signf(ixO^S)+dabs(f(ixO^S)))*third))
     case(5)
       ! van Leer limiter
       lf(ixO^S)=two*max(f(ixB^S)*f(ixO^S),zero)/(f(ixO^S)+f(ixB^S)+qsmall)
     case default
       call mpistop("Unknown slope limiter for thermal conduction")
    end select
  end function slope_limiter

  !> Get the explicit timestep for the TC (hd implementation)
  !> Note: also used in 1D MHD (or for neutrals in twofl)
  function get_tc_dt_hd(w,ixI^L,ixO^L,dx^D,x,fl)  result(dtnew)
    ! Check diffusion time limit dt < dx_i**2 / ((gamma-1)*tc_k_para/rho)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    type(tc_fluid), intent(in) :: fl
    double precision :: dtnew

    double precision :: tmp(ixO^S),tmp2(ixO^S),Te(ixI^S),rho(ixI^S),hfs(ixO^S),gradT(ixI^S)
    double precision :: ne(ixI^S), nH_arr(ixI^S)
    double precision :: dtdiff_tcond,maxtmp2
    integer          :: idim

    call fl%get_temperature_from_conserved(w,x,ixI^L,ixI^L,Te)
    call fl%get_rho(w,x,ixI^L,ixO^L,rho)
    call fl%get_ne_nH(ixI^L, ixO^L, w, x, ne, nH_arr)

    if(fl%tc_constant) then
      tmp(ixO^S)=fl%tc_k_para/rho(ixO^S)
    else
      if(fl%tc_saturate) then
        ! Kannan 2016 MN 458, 410
        ! l_mfpe = 3^1.5*kB^2/(4*sqrt(pi)*e^4*lnLambda) * T^2/n_e
        if(SI_unit) then
          ! 5.730205638843984e27 = 3^1.5*kB_SI^2/(4*sqrt(pi)*e_SI^4*37)
          tmp2(ixO^S)=Te(ixO^S)**2/ne(ixO^S)*5.730205638843984d27*unit_temperature**2/(unit_numberdensity*unit_length)
        else
          ! 7093.9239487765044 = 3^1.5*kB_cgs^2/(4*sqrt(pi)*e_cgs^4*37)
          tmp2(ixO^S)=Te(ixO^S)**2/ne(ixO^S)*7093.9239487765044d0*unit_temperature**2/(unit_numberdensity*unit_length)
        end if
        hfs=0.d0
        do idim=1,ndim
          call gradient(Te,ixI^L,ixO^L,idim,gradT)
          hfs(ixO^S)=hfs(ixO^S)+gradT(ixO^S)**2
        end do
        ! kappa=kappa_Spitzer/(1+4.2*l_mfpe/(T/|gradT|))
        tmp(ixO^S)=fl%tc_k_para*Te(ixO^S)*Te(ixO^S)*dsqrt(Te(ixO^S))/(rho(ixO^S)*(1.d0+4.2d0*tmp2(ixO^S)*dsqrt(hfs(ixO^S))/Te(ixO^S)))
      else
        tmp(ixO^S)=fl%tc_k_para*Te(ixO^S)*Te(ixO^S)*dsqrt(Te(ixO^S))/rho(ixO^S)
      end if
    end if

    dtnew = bigdouble
    do idim=1,ndim
      ! approximate thermal conduction flux: tc_k_para/rho/dx**2
      maxtmp2=maxval(tmp(ixO^S)/(block%ds(ixO^S,idim)**2))
      ! dt< dx_idim**2/((gamma-1)*tc_k_para/rho)
      dtdiff_tcond=1.d0/(tc_gamma_1*maxtmp2+smalldouble)
      ! limit the time step
      dtnew=min(dtnew,dtdiff_tcond)
    end do
    dtnew=dtnew/dble(ndim)
  end function get_tc_dt_hd

  subroutine sts_set_source_tc_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,fl)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    ! intent(inout) so tc_patch_negative_eint can repair w(:, ie) in place
    ! when tc_patch_eint = .true.; legacy path (patch off) leaves w untouched.
    double precision, intent(inout) ::  w(ixI^S,1:nw)
    double precision, intent(inout) ::  wres(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    type(tc_fluid), intent(in)    :: fl

    double precision :: Te(ixI^S),rho(ixI^S)
    double precision :: qvec(ixI^S,1:ndim),qd(ixI^S)
    double precision, allocatable, dimension(:^D&,:) :: qvec_equi
    double precision :: fluxall(ixI^S,1,1:ndim)

    double precision :: dxinv(ndim)
    integer :: idims,ix^L,ixB^L,ixA^L

    ix^L=ixO^L^LADD1;

    dxinv=1.d0/dxlevel

    !calculate Te in whole domain (+ghosts)
    call fl%get_temperature_from_eint(w, x, ixI^L, ixI^L, Te)
    if (fl%tc_patch_eint) call tc_patch_negative_eint(w, x, ixI^L, ixI^L, Te, fl%e_, fl)
    call fl%get_rho(w, x, ixI^L, ixI^L, rho)
    if(slab_uniform) then
      call set_source_tc_hd(ixI^L,ixO^L,w,x,fl,qvec,rho,Te)
    else
      call set_source_tc_hd_geo(ixI^L,ixO^L,w,x,fl,qvec,rho,Te)
    end if
    if(fl%subtract_equi) then
      allocate(qvec_equi(ixI^S,1:ndim))
      call fl%get_temperature_equi(w, x, ixI^L, ixI^L, Te)  !calculate Te in whole domain (+ghosts)
      call fl%get_rho_equi(w, x, ixI^L, ixI^L, rho)  !calculate rho in whole domain (+ghosts)
      if(slab_uniform) then
        call set_source_tc_hd(ixI^L,ixO^L,w,x,fl,qvec_equi,rho,Te)
      else
        call set_source_tc_hd_geo(ixI^L,ixO^L,w,x,fl,qvec_equi,rho,Te)
      end if
      do idims=1,ndim
        ! operate only on the per-idims face range that set_source_tc_hd fills
        ! (ixOmin-kr .. ixOmax); ix^S (=ixO+1) over-reaches into the
        ! uninitialised ixOmax+1 layer -> snan/Inf -> SIGFPE. Matches tc_mhd.
        ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
        qvec(ixA^S,idims)=qvec(ixA^S,idims) - qvec_equi(ixA^S,idims)
      end do
      deallocate(qvec_equi)
    endif

    qd=0.d0
    if(slab_uniform) then
      do idims=1,ndim
        ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
        qvec(ixA^S,idims)=dxinv(idims)*qvec(ixA^S,idims)
        ixB^L=ixO^L-kr(idims,^D);
        qd(ixO^S)=qd(ixO^S)+qvec(ixO^S,idims)-qvec(ixB^S,idims)
      end do
    else
      do idims=1,ndim
        ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
        qvec(ixA^S,idims)=qvec(ixA^S,idims)*block%surfaceC(ixA^S,idims)
        ixB^L=ixO^L-kr(idims,^D);
        qd(ixO^S)=qd(ixO^S)+qvec(ixO^S,idims)-qvec(ixB^S,idims)
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
    integer :: idims,ix^D,ix^L,ixC^L,ixA^L,ixB^L

    ix^L=ixO^L^LADD1;
    ! ixC is cell-corner index
    ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;

    ! calculate thermal conduction flux with symmetric scheme
    ! T gradient (central difference) at cell corners
    do idims=1,ndim
      ixBmin^D=ixmin^D;
      ixBmax^D=ixmax^D-kr(idims,^D);
      call gradientF(Te,x,ixI^L,ixB^L,idims,ke)
     {^IFTHREED
      if(idims==1) then
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qvec(ix^D,idims)=0.25d0*(ke(ix1,ix2,ix3)+ke(ix1,ix2+1,ix3)&
                         +ke(ix1,ix2,ix3+1)+ke(ix1,ix2+1,ix3+1))
       {end do\}
      else if(idims==2) then
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qvec(ix^D,idims)=0.25d0*(ke(ix1,ix2,ix3)+ke(ix1+1,ix2,ix3)&
                         +ke(ix1,ix2,ix3+1)+ke(ix1+1,ix2,ix3+1))
       {end do\}
      else
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qvec(ix^D,idims)=0.25d0*(ke(ix1,ix2,ix3)+ke(ix1+1,ix2,ix3)&
                         +ke(ix1,ix2+1,ix3)+ke(ix1+1,ix2+1,ix3))
       {end do\}
      end if
     }
     {^IFTWOD
      if(idims==1) then
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qvec(ix^D,idims)=0.5d0*(ke(ix1,ix2)+ke(ix1,ix2+1))
       {end do\}
      else
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qvec(ix^D,idims)=0.5d0*(ke(ix1,ix2)+ke(ix1+1,ix2))
       {end do\}
      end if
     }
     {^IFONED
     do ix1=ixCmin1,ixCmax1
       qvec(ix1,idims)=ke(ix1)
     end do
     }
    end do
    ! conductivity at cell center
    if(fl%tc_constant) then
      qd(ix^S)=fl%tc_k_para
    else
      if(phys_trac) then
       {do ix^DB=ixmin^DB,ixmax^DB\}
          if(Te(ix^D) < block%wextra(ix^D,fl%Tcoff_)) then
            qd(ix^D)=fl%tc_k_para*dsqrt(block%wextra(ix^D,fl%Tcoff_)**5)
          else
            qd(ix^D)=fl%tc_k_para*dsqrt(Te(ix^D)**5)
          end if
       {end do\}
      else
        qd(ix^S)=fl%tc_k_para*dsqrt(Te(ix^S)**5)
      end if
    end if
    ! conductivity Ke at cell corner
    ! cell corner conduction flux gradT
    {^IFTHREED
    {do ix^DB=ixCmin^DB,ixCmax^DB\}
       ke(ix^D)=0.125d0*(qd(ix1,ix2,ix3)+qd(ix1+1,ix2,ix3)&
                      +qd(ix1,ix2+1,ix3)+qd(ix1+1,ix2+1,ix3)&
                      +qd(ix1,ix2,ix3+1)+qd(ix1+1,ix2,ix3+1)&
                      +qd(ix1,ix2+1,ix3+1)+qd(ix1+1,ix2+1,ix3+1))
       gradT(ix^D,1)=ke(ix^D)*qvec(ix^D,1)
       gradT(ix^D,2)=ke(ix^D)*qvec(ix^D,2)
       gradT(ix^D,3)=ke(ix^D)*qvec(ix^D,3)
    {end do\}
    }
    {^IFTWOD
    {do ix^DB=ixCmin^DB,ixCmax^DB\}
       ke(ix^D)=0.25d0*(qd(ix1,ix2)+qd(ix1+1,ix2)+qd(ix1,ix2+1)+qd(ix1+1,ix2+1))
       gradT(ix^D,1)=ke(ix^D)*qvec(ix^D,1)
       gradT(ix^D,2)=ke(ix^D)*qvec(ix^D,2)
    {end do\}
    }
    {^IFONED
     do ix1=ixCmin1,ixCmax1
       gradT(ix^D,1)=0.5d0*(qd(ix1)+qd(ix1+1))*qvec(ix^D,1)
     end do
    }

    ! conduction flux qvec at cell face
    do idims=1,ndim
      ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
     {^IFTHREED
      if(idims==1) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          qvec(ix^D,idims)=0.25d0*(gradT(ix1,ix2,ix3,idims)+gradT(ix1,ix2-1,ix3,idims)&
                                +gradT(ix1,ix2,ix3-1,idims)+gradT(ix1,ix2-1,ix3-1,idims))
       {end do\}
      else if(idims==2) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          qvec(ix^D,idims)=0.25d0*(gradT(ix1,ix2,ix3,idims)+gradT(ix1-1,ix2,ix3,idims)&
                                +gradT(ix1,ix2,ix3-1,idims)+gradT(ix1-1,ix2,ix3-1,idims))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          qvec(ix^D,idims)=0.25d0*(gradT(ix1,ix2,ix3,idims)+gradT(ix1,ix2-1,ix3,idims)&
                                +gradT(ix1-1,ix2,ix3,idims)+gradT(ix1-1,ix2-1,ix3,idims))
       {end do\}
      end if
     }
     {^IFTWOD
      if(idims==1) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          qvec(ix^D,idims)=0.5d0*(gradT(ix1,ix2,idims)+gradT(ix1,ix2-1,idims))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          qvec(ix^D,idims)=0.5d0*(gradT(ix1,ix2,idims)+gradT(ix1-1,ix2,idims))
       {end do\}
      end if
     }
     {^IFONED
      do ix1=ixAmin1,ixAmax1
        qvec(ix1,idims)=gradT(ix1,idims)
      end do
     }
      if(fl%tc_saturate) then
        ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135: phi=1.1, Balbus and Mckee 1982 ApJ, 252, 529: phi=0.3)
        ! unsigned saturated TC flux = 5 phi rho c**3, c=sqrt(p/rho) is isothermal sound speed, phi=0.3
        ixB^L=ixA^L+kr(idims,^D);
        qd(ixA^S)=0.75d0*(rho(ixA^S)+rho(ixB^S))*dsqrt(0.5d0*(Te(ixA^S)+Te(ixB^S)))**3
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          if(dabs(qvec(ix^D,idims))>qd(ix^D)) then
            qvec(ix^D,idims)=sign(1.d0,qvec(ix^D,idims))*qd(ix^D)
          end if
       {end do\}
      end if
    end do
  end subroutine set_source_tc_hd
  subroutine set_source_tc_hd_geo(ixI^L,ixO^L,w,x,fl,qvec,rho,Te)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(in) ::  w(ixI^S,1:nw)
    type(tc_fluid), intent(in)    :: fl
    double precision, intent(in) :: Te(ixI^S),rho(ixI^S)
    double precision, intent(out) :: qvec(ixI^S,1:ndim)
    double precision :: gradT(ixI^S,1:ndim),ke(ixI^S),qd(ixI^S)
    integer :: idims,ix^D,ix^L,ixC^L,ixA^L,ixB^L

    ix^L=ixO^L^LADD1;
    ! ixC is cell-corner index
    ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;

    ! calculate thermal conduction flux with symmetric scheme
    ! T gradient from face centers to cell corners: surface weighted average
    do idims=1,ndim
      ixBmin^D=ixmin^D;
      ixBmax^D=ixmax^D-kr(idims,^D);
      call gradientF(Te,x,ixI^L,ixB^L,idims,ke)
     {^IFTHREED
      if(idims==1) then
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qvec(ix^D,1)=(ke(ix1,ix2,ix3)*block%surfaceC(ix1,ix2,ix3,1)&
                     +ke(ix1,ix2+1,ix3)*block%surfaceC(ix1,ix2+1,ix3,1)&
                     +ke(ix1,ix2,ix3+1)*block%surfaceC(ix1,ix2,ix3+1,1)&
                   +ke(ix1,ix2+1,ix3+1)*block%surfaceC(ix1,ix2+1,ix3+1,1))/&
            (block%surfaceC(ix1,ix2,ix3,1)+block%surfaceC(ix1,ix2+1,ix3,1)&
          +block%surfaceC(ix1,ix2,ix3+1,1)+block%surfaceC(ix1,ix2+1,ix3+1,1)+smalldouble**2)
       {end do\}
      else if(idims==2) then
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qvec(ix^D,2)=(ke(ix1,ix2,ix3)*block%surfaceC(ix1,ix2,ix3,2)&
                     +ke(ix1+1,ix2,ix3)*block%surfaceC(ix1+1,ix2,ix3,2)&
                     +ke(ix1,ix2,ix3+1)*block%surfaceC(ix1,ix2,ix3+1,2)&
                   +ke(ix1+1,ix2,ix3+1)*block%surfaceC(ix1+1,ix2,ix3+1,2))/&
            (block%surfaceC(ix1,ix2,ix3,2)+block%surfaceC(ix1+1,ix2,ix3,2)&
          +block%surfaceC(ix1,ix2,ix3+1,2)+block%surfaceC(ix1+1,ix2,ix3+1,2)+smalldouble**2)
          ! zero theta-normal surface area at pole axis
       {end do\}
      else
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qvec(ix^D,3)=(ke(ix1,ix2,ix3)*block%surfaceC(ix1,ix2,ix3,3)&
                     +ke(ix1+1,ix2,ix3)*block%surfaceC(ix1+1,ix2,ix3,3)&
                     +ke(ix1,ix2+1,ix3)*block%surfaceC(ix1,ix2+1,ix3,3)&
                   +ke(ix1+1,ix2+1,ix3)*block%surfaceC(ix1+1,ix2+1,ix3,3))/&
            (block%surfaceC(ix1,ix2,ix3,3)+block%surfaceC(ix1+1,ix2,ix3,3)&
          +block%surfaceC(ix1,ix2+1,ix3,3)+block%surfaceC(ix1+1,ix2+1,ix3,3))
       {end do\}
      end if
     }
     {^IFTWOD
      if(idims==1) then
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qvec(ix^D,1)=(ke(ix1,ix2)*block%surfaceC(ix1,ix2,1)+ke(ix1,ix2+1)*block%surfaceC(ix1,ix2+1,1))&
                      /(block%surfaceC(ix1,ix2,1)+block%surfaceC(ix1,ix2+1,1))
       {end do\}
      else
       {do ix^DB=ixCmin^DB,ixCmax^DB\}
          qvec(ix^D,2)=(ke(ix1,ix2)*block%surfaceC(ix1,ix2,2)+ke(ix1+1,ix2)*block%surfaceC(ix1+1,ix2,2))&
                      /(block%surfaceC(ix1,ix2,2)+block%surfaceC(ix1+1,ix2,2))
       {end do\}
      end if
     }
     {^IFONED
     do ix1=ixCmin1,ixCmax1
       qvec(ix1,idims)=ke(ix1)
     end do
     }
    end do
    ! conductivity at cell center
    if(fl%tc_constant) then
      qd(ix^S)=fl%tc_k_para
    else
      if(phys_trac) then
       {do ix^DB=ixmin^DB,ixmax^DB\}
          if(Te(ix^D) < block%wextra(ix^D,fl%Tcoff_)) then
            qd(ix^D)=fl%tc_k_para*dsqrt(block%wextra(ix^D,fl%Tcoff_)**5)
          else
            qd(ix^D)=fl%tc_k_para*dsqrt(Te(ix^D)**5)
          end if
       {end do\}
      else
        qd(ix^S)=fl%tc_k_para*dsqrt(Te(ix^S)**5)
      end if
    end if
    ! conductivity Ke at cell corner
    ! cell corner conduction flux gradT
    {^IFTHREED
    {do ix^DB=ixCmin^DB,ixCmax^DB\}
       ke(ix^D)=(qd(ix1,ix2,ix3)*block%dvolume(ix1,ix2,ix3)+qd(ix1+1,ix2,ix3)*block%dvolume(ix1+1,ix2,ix3)&
              +qd(ix1,ix2+1,ix3)*block%dvolume(ix1,ix2+1,ix3)+qd(ix1+1,ix2+1,ix3)*block%dvolume(ix1+1,ix2+1,ix3)&
              +qd(ix1,ix2,ix3+1)*block%dvolume(ix1,ix2,ix3+1)+qd(ix1+1,ix2,ix3+1)*block%dvolume(ix1+1,ix2,ix3+1)&
            +qd(ix1,ix2+1,ix3+1)*block%dvolume(ix1,ix2+1,ix3+1)+qd(ix1+1,ix2+1,ix3+1)*block%dvolume(ix1+1,ix2+1,ix3+1))&
           /(block%dvolume(ix1,ix2,ix3)+block%dvolume(ix1+1,ix2,ix3)+block%dvolume(ix1,ix2+1,ix3)+block%dvolume(ix1+1,ix2+1,ix3)&
        +block%dvolume(ix1,ix2,ix3+1)+block%dvolume(ix1+1,ix2,ix3+1)+block%dvolume(ix1,ix2+1,ix3+1)+block%dvolume(ix1+1,ix2+1,ix3+1))
       gradT(ix^D,1)=ke(ix^D)*qvec(ix^D,1)
       gradT(ix^D,2)=ke(ix^D)*qvec(ix^D,2)
       gradT(ix^D,3)=ke(ix^D)*qvec(ix^D,3)
    {end do\}
    }
    {^IFTWOD
    {do ix^DB=ixCmin^DB,ixCmax^DB\}
       ke(ix^D)=(qd(ix1,ix2)*block%dvolume(ix1,ix2)+qd(ix1+1,ix2)*block%dvolume(ix1+1,ix2)&
              +qd(ix1,ix2+1)*block%dvolume(ix1,ix2+1)+qd(ix1+1,ix2+1)*block%dvolume(ix1+1,ix2+1))&
             /(block%dvolume(ix1,ix2)+block%dvolume(ix1+1,ix2)+block%dvolume(ix1,ix2+1)+block%dvolume(ix1+1,ix2+1))
       gradT(ix^D,1)=ke(ix^D)*qvec(ix^D,1)
       gradT(ix^D,2)=ke(ix^D)*qvec(ix^D,2)
    {end do\}
    }
    {^IFONED
     do ix1=ixCmin1,ixCmax1
       gradT(ix^D,1)=(qd(ix1)*block%dvolume(ix1)+qd(ix1+1)*block%dvolume(ix1+1))/(block%dvolume(ix1)+block%dvolume(ix1+1))*qvec(ix^D,1)
     end do
    }

    ! conduction flux qvec at cell face
    do idims=1,ndim
      ixAmax^D=ixOmax^D; ixAmin^D=ixOmin^D-kr(idims,^D);
     {^IFTHREED
      if(idims==1) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          qvec(ix^D,idims)=0.25d0*(gradT(ix1,ix2,ix3,idims)+gradT(ix1,ix2-1,ix3,idims)&
                                +gradT(ix1,ix2,ix3-1,idims)+gradT(ix1,ix2-1,ix3-1,idims))
       {end do\}
      else if(idims==2) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          qvec(ix^D,idims)=0.25d0*(gradT(ix1,ix2,ix3,idims)+gradT(ix1-1,ix2,ix3,idims)&
                                +gradT(ix1,ix2,ix3-1,idims)+gradT(ix1-1,ix2,ix3-1,idims))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          qvec(ix^D,idims)=0.25d0*(gradT(ix1,ix2,ix3,idims)+gradT(ix1,ix2-1,ix3,idims)&
                                +gradT(ix1-1,ix2,ix3,idims)+gradT(ix1-1,ix2-1,ix3,idims))
       {end do\}
      end if
     }
     {^IFTWOD
      if(idims==1) then
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          qvec(ix^D,idims)=0.5d0*(gradT(ix1,ix2,idims)+gradT(ix1,ix2-1,idims))
       {end do\}
      else
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          qvec(ix^D,idims)=0.5d0*(gradT(ix1,ix2,idims)+gradT(ix1-1,ix2,idims))
       {end do\}
      end if
     }
     {^IFONED
      do ix1=ixAmin1,ixAmax1
        qvec(ix1,idims)=gradT(ix1,idims)
      end do
     }
      if(fl%tc_saturate) then
        ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135: phi=1.1, Balbus and Mckee 1982 ApJ, 252, 529: phi=0.3)
        ! unsigned saturated TC flux = 5 phi rho c**3, c=sqrt(p/rho) is isothermal sound speed, phi=0.3
        ixB^L=ixA^L+kr(idims,^D);
        qd(ixA^S)=0.75d0*(rho(ixA^S)+rho(ixB^S))*dsqrt(0.5d0*(Te(ixA^S)+Te(ixB^S)))**3
       {do ix^DB=ixAmin^DB,ixAmax^DB\}
          if(dabs(qvec(ix^D,idims))>qd(ix^D)) then
            qvec(ix^D,idims)=sign(1.d0,qvec(ix^D,idims))*qd(ix^D)
          end if
       {end do\}
      end if
    end do
  end subroutine set_source_tc_hd_geo

  !> Patch cells where e_int <= 0 by neighbor-averaging the temperature
  !> AND repairing the conserved internal energy w(:, ie) in place.
  !> Called after get_temperature_from_eint in the TC source routines when
  !> tc_patch_eint is .true.  During STS RKL2 Chebyshev substeps the
  !> polynomial can overshoot e_int to negative values at low-density
  !> coronal cells.  A negative T would produce NaN via sqrt(T) in the
  !> Spitzer conductivity.  This routine replaces those cells' T with
  !> the average of valid (e_int > 0) neighbors and writes the matching
  !> e_int back to w(ix, ie) via the EoS inverse helper eint_nH_from_T,
  !> so the next STS substep does not re-read negative w and amplify it.
  subroutine tc_patch_negative_eint(w, x, ixI^L, ixO^L, Te, ie, fl)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, ie
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: Te(ixI^S)
    type(tc_fluid), intent(in)      :: fl

    integer :: ix^D, count
    integer :: ipatch(ixI^S)
    double precision :: T_avg, T_use, rho_cell, nH_cell, log_nH, log_T
    double precision :: eint_new, small_e_local, log_T_min, T_floor
    double precision :: rho(ixI^S)

    ipatch(ixI^S) = 0

    ! Mark cells with negative e_int
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if (w(ix^D, ie) <= 0.0d0) ipatch(ix^D) = 1
    {end do\}

    ! Quick exit if no cells need patching
    if (all(ipatch(ixO^S) == 0)) return

    ! Need rho per-cell to back out e_int from a prescribed T via the EoS.
    call fl%get_rho(w, x, ixI^L, ixI^L, rho)

    ! Floor on log_T for the (rho, T) inverse table, from the EoS via the port
    ! (eos_get_log_T_floor encapsulates the method->container choice; picking the
    ! wrong container would leave var2_min=0 and clobber cold cells).
    log_T_min = fl%log_T_floor

    ! Local minimum e_int (code units): equivalent to small_pressure / (gamma-1).
    small_e_local = small_pressure * fl%inv_gamma_minus_1

    ! Fallback T (code units) for cells that have no valid neighbour.
    T_floor = 1.0d0 / unit_temperature

    ! Replace marked cells with neighbor-averaged T, then back out e_int
    ! and floor w(:, ie) so subsequent STS substeps see a positive state.
    {do ix^DB=ixOmin^DB+1,ixOmax^DB-1\}
        if (ipatch(ix^D) == 1) then
            T_avg = 0.0d0
            count = 0
            {^IFONED
            if (ipatch(ix1-1)==0) then; T_avg=T_avg+Te(ix1-1); count=count+1; end if
            if (ipatch(ix1+1)==0) then; T_avg=T_avg+Te(ix1+1); count=count+1; end if
            }
            {^IFTWOD
            if (ipatch(ix1-1,ix2)==0) then; T_avg=T_avg+Te(ix1-1,ix2); count=count+1; end if
            if (ipatch(ix1+1,ix2)==0) then; T_avg=T_avg+Te(ix1+1,ix2); count=count+1; end if
            if (ipatch(ix1,ix2-1)==0) then; T_avg=T_avg+Te(ix1,ix2-1); count=count+1; end if
            if (ipatch(ix1,ix2+1)==0) then; T_avg=T_avg+Te(ix1,ix2+1); count=count+1; end if
            }
            {^IFTHREED
            if (ipatch(ix1-1,ix2,ix3)==0) then; T_avg=T_avg+Te(ix1-1,ix2,ix3); count=count+1; end if
            if (ipatch(ix1+1,ix2,ix3)==0) then; T_avg=T_avg+Te(ix1+1,ix2,ix3); count=count+1; end if
            if (ipatch(ix1,ix2-1,ix3)==0) then; T_avg=T_avg+Te(ix1,ix2-1,ix3); count=count+1; end if
            if (ipatch(ix1,ix2+1,ix3)==0) then; T_avg=T_avg+Te(ix1,ix2+1,ix3); count=count+1; end if
            if (ipatch(ix1,ix2,ix3-1)==0) then; T_avg=T_avg+Te(ix1,ix2,ix3-1); count=count+1; end if
            if (ipatch(ix1,ix2,ix3+1)==0) then; T_avg=T_avg+Te(ix1,ix2,ix3+1); count=count+1; end if
            }
            if (count > 0) then
                T_use = T_avg / dble(count)
            else
                T_use = T_floor
            end if

            ! Repair the conserved internal energy in place. Use the same
            ! EoS-aware inverse used by hd/mhd_from_prolong_LTE: e_int =
            ! nH * eint_nH_from_T(log_nH, log_T).  Floor with small_e.
            rho_cell = max(rho(ix^D), small_density)
            nH_cell  = rho_cell / fl%nH2rhoFactor
            log_nH   = dlog10(max(nH_cell, smalldouble))
            log_T    = dlog10(max(T_use, 10.0d0**log_T_min))
            eint_new = nH_cell * fl%eint_from_T(log_nH, log_T)
            w(ix^D, ie) = max(small_e_local, eint_new)

            ! Keep Te consistent with the floored w(:, ie): if small_e
            ! kicked in, the implied T is slightly different from T_use.
            ! Use the prescribed T_use here (Te is only consumed within
            ! this TC substep; the next call recomputes it from w).
            Te(ix^D) = T_use

            write(*,*) ' WARNING: tc_patch_eint it=', it, &
                ' e_int=', w(ix^D, ie), &
                ' pe=', mype, ' x=', x(ix^D, 1:ndim)
            ipatch(ix^D) = 2  ! mark as handled
        end if
    {end do\}

    ! Cells at the block boundary that the interior loop above skipped
    ! (ixOmin/ixOmax edges).  Fall back to T = 1/unit_temperature and
    ! repair w(:, ie) consistently so STS does not re-read negative e.
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if (ipatch(ix^D) == 1) then
            T_use    = T_floor
            rho_cell = max(rho(ix^D), small_density)
            nH_cell  = rho_cell / fl%nH2rhoFactor
            log_nH   = dlog10(max(nH_cell, smalldouble))
            log_T    = dlog10(max(T_use, 10.0d0**log_T_min))
            eint_new = nH_cell * fl%eint_from_T(log_nH, log_T)
            w(ix^D, ie) = max(small_e_local, eint_new)
            Te(ix^D) = T_use
            write(*,*) ' WARNING: tc_patch_eint it=', it, &
                ' e_int=', w(ix^D, ie), &
                ' pe=', mype, ' x=', x(ix^D, 1:ndim), ' (edge)'
        end if
    {end do\}

  end subroutine tc_patch_negative_eint

end module mod_thermal_conduction
