
!############### MODULE ##############
module krome_commons
  implicit none

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2018-01-25 18:47:28
  !  Changeset 9a6e335
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************
  integer,parameter::idx_E=1
  integer,parameter::idx_Ok=2
  integer,parameter::idx_Hk=3
  integer,parameter::idx_Sk=4
  integer,parameter::idx_Ck=5
  integer,parameter::idx_H=6
  integer,parameter::idx_HE=7
  integer,parameter::idx_C=8
  integer,parameter::idx_NO=9
  integer,parameter::idx_CO=10
  integer,parameter::idx_N=11
  integer,parameter::idx_O2=12
  integer,parameter::idx_S=13
  integer,parameter::idx_SO=14
  integer,parameter::idx_O=15
  integer,parameter::idx_H2=16
  integer,parameter::idx_SI=17
  integer,parameter::idx_OH=18
  integer,parameter::idx_HS=19
  integer,parameter::idx_NS=20
  integer,parameter::idx_H2S=21
  integer,parameter::idx_FE=22
  integer,parameter::idx_CS=23
  integer,parameter::idx_CN=24
  integer,parameter::idx_S2=25
  integer,parameter::idx_NA=26
  integer,parameter::idx_F=27
  integer,parameter::idx_HF=28
  integer,parameter::idx_CH=29
  integer,parameter::idx_SO2=30
  integer,parameter::idx_C2=31
  integer,parameter::idx_N2=32
  integer,parameter::idx_CH2=33
  integer,parameter::idx_NH=34
  integer,parameter::idx_HCN=35
  integer,parameter::idx_CO2=36
  integer,parameter::idx_SIO=37
  integer,parameter::idx_SIO2=38
  integer,parameter::idx_NH2=39
  integer,parameter::idx_OCN=40
  integer,parameter::idx_MG=41
  integer,parameter::idx_P=42
  integer,parameter::idx_HCO=43
  integer,parameter::idx_H2O=44
  integer,parameter::idx_OCS=45
  integer,parameter::idx_PN=46
  integer,parameter::idx_PO=47
  integer,parameter::idx_HEj=48
  integer,parameter::idx_Hj=49
  integer,parameter::idx_NHj=50
  integer,parameter::idx_HSj=51
  integer,parameter::idx_Sj=52
  integer,parameter::idx_SIj=53
  integer,parameter::idx_OHj=54
  integer,parameter::idx_HEHj=55
  integer,parameter::idx_H2j=56
  integer,parameter::idx_FEj=57
  integer,parameter::idx_SIHj=58
  integer,parameter::idx_NAj=59
  integer,parameter::idx_HCOj=60
  integer,parameter::idx_CHj=61
  integer,parameter::idx_Oj=62
  integer,parameter::idx_MGj=63
  integer,parameter::idx_SIOj=64
  integer,parameter::idx_Pj=65
  integer,parameter::idx_SIFj=66
  integer,parameter::idx_Cj=67
  integer,parameter::idx_Nj=68
  integer,parameter::idx_COj=69
  integer,parameter::idx_Fj=70
  integer,parameter::idx_CR=71
  integer,parameter::idx_g=72
  integer,parameter::idx_Tgas=73
  integer,parameter::idx_dummy=74
  integer,parameter::nrea=255
  integer,parameter::nmols=70
  integer,parameter::nspec=74
  integer,parameter::natoms=13
  integer,parameter::ndust=0
  integer,parameter::ndustTypes=0
  integer,parameter::nPhotoBins=0
  integer,parameter::nPhotoRea=0

  integer,parameter::idx_atom_C=1
  integer,parameter::idx_atom_E=2
  integer,parameter::idx_atom_F=3
  integer,parameter::idx_atom_Mg=4
  integer,parameter::idx_atom_H=5
  integer,parameter::idx_atom_Si=6
  integer,parameter::idx_atom_O=7
  integer,parameter::idx_atom_N=8
  integer,parameter::idx_atom_P=9
  integer,parameter::idx_atom_S=10
  integer,parameter::idx_atom_Fe=11
  integer,parameter::idx_atom_Na=12
  integer,parameter::idx_atom_He=13

  !cooling index
  integer,parameter::idx_cool_h2 = 1
  integer,parameter::idx_cool_h2gp = 2
  integer,parameter::idx_cool_atomic = 3
  integer,parameter::idx_cool_cen = 3
  integer,parameter::idx_cool_hd = 4
  integer,parameter::idx_cool_metal = 5
  integer,parameter::idx_cool_z = 5
  integer,parameter::idx_cool_dh = 6
  integer,parameter::idx_cool_enthalpic = 6
  integer,parameter::idx_cool_dust = 7
  integer,parameter::idx_cool_compton = 8
  integer,parameter::idx_cool_cie = 9
  integer,parameter::idx_cool_cont = 10
  integer,parameter::idx_cool_continuum = 10
  integer,parameter::idx_cool_expansion = 11
  integer,parameter::idx_cool_exp = 11
  integer,parameter::idx_cool_ff = 12
  integer,parameter::idx_cool_bss = 12
  integer,parameter::idx_cool_custom = 13
  integer,parameter::idx_cool_co = 14
  integer,parameter::idx_cool_zcie = 15
  integer,parameter::idx_cool_zcienouv = 16
  integer,parameter::idx_cool_zextend = 17
  integer,parameter::idx_cool_gh = 18
  integer,parameter::ncools = 18

  !heating index
  integer,parameter::idx_heat_chem = 1
  integer,parameter::idx_heat_compress = 2
  integer,parameter::idx_heat_compr = 2
  integer,parameter::idx_heat_photo = 3
  integer,parameter::idx_heat_dh = 4
  integer,parameter::idx_heat_enthalpic = 4
  integer,parameter::idx_heat_av = 5
  integer,parameter::idx_heat_photoav = 5
  integer,parameter::idx_heat_cr = 6
  integer,parameter::idx_heat_dust = 7
  integer,parameter::idx_heat_xray = 8
  integer,parameter::idx_heat_viscous = 9
  integer,parameter::idx_heat_visc = 9
  integer,parameter::idx_heat_custom = 10
  integer,parameter::idx_heat_zcie = 11
  integer,parameter::nheats = 11

  real*8::arr_k(nrea)

  !commons for rate tables
  !modify ktab_n according to the required precision
  integer,parameter::ktab_n=int(1e3)
  real*8::ktab(nrea,ktab_n),ktab_logTlow, ktab_logTup, ktab_T(ktab_n)
  real*8::inv_ktab_T(ktab_n-1), inv_ktab_idx

  !thermo toggle (when >0 do cooling/heating)
  integer::krome_thermo_toggle
  !$omp threadprivate(krome_thermo_toggle)

  !debug bit flag, print and array with fallback values for extreme environments
  integer:: red_flag
  real*8::n_global(nspec)
  integer, save :: nprint_negative=10
  !$omp threadprivate(n_global,nprint_negative,red_flag)

  !commons for implicit RHS
  integer::arr_r1(nrea)
  integer::arr_r2(nrea)
  integer::arr_r3(nrea)
  integer::arr_p1(nrea)
  integer::arr_p2(nrea)
  integer::arr_p3(nrea)

  !commons for reduction
  integer::arr_u(nrea)
  real*8::arr_flux(nrea)

  !commons for frequency bins

  !commons for H2 photodissociation (Solomon)
  ! note: paramters here are set depending on the data
  ! but if you have a different file you should modify them
  integer,parameter::H2pdData_nvibX=15
  integer,parameter::H2pdData_nvibB=37
  real*8::H2pdData_dE(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_pre(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_EX(H2pdData_nvibX)
  integer::H2pdData_binMap(H2pdData_nvibX,H2pdData_nvibB)

  !commons for dust optical properties

  !square of turbulence velocity for broadening
  real*8::broadeningVturb2

  !mpi rank of process. If 0, ignored
  integer::krome_mpi_rank=0, krome_omp_thread
  !$omp threadprivate(krome_omp_thread)

  !user-defined commons variables from the reaction file
  real*8::user_crflux
  !$omp threadprivate(user_crflux)

  !commons for anytab

  !physical commons
  real*8::phys_Tcmb
  real*8::phys_zredshift
  real*8::phys_orthoParaRatio
  real*8::phys_metallicity
  real*8::phys_Tfloor
  !$omp threadprivate(phys_Tcmb)
  !$omp threadprivate(phys_zredshift)
  !$omp threadprivate(phys_orthoParaRatio)
  !$omp threadprivate(phys_metallicity)
  !$omp threadprivate(phys_Tfloor)

  !machine precision
  real*8::krome_epsilon

  !xrayJ21 for tabulated heating and rate
  real*8::J21xray

  !total metallicity relative to solar Z/Z_solar
  real*8::total_Z
  real*8::dust2gas_ratio

  !commons for dust tabs (cool,H2,Tdust)
  integer,parameter::dust_tab_imax=50, dust_tab_jmax=50
  real*8::dust_tab_ngas(dust_tab_imax)
  real*8::dust_tab_Tgas(dust_tab_jmax)
  real*8::dust_mult_Tgas,dust_mult_ngas
  real*8::dust_table_AvVariable_log

  real*8::dust_tab_cool(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_heat(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_Tdust(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_H2(dust_tab_imax, dust_tab_jmax)

  !commons for exp(-a) table
  integer,parameter::exp_table_na=int(1d5)
  real*8,parameter::exp_table_aMax=1d4,exp_table_aMin=0d0
  real*8,parameter::exp_table_multa=(exp_table_na-1) &
      / (exp_table_aMax-exp_table_aMin)
  real*8,parameter::exp_table_da=1d0/exp_table_multa
  real*8::exp_table(exp_table_na)

  !stores the last evaluation of the rates in the fex
  real*8::last_coe(nrea)
  !$omp threadprivate(last_coe)

  !data for CO cooling
  integer,parameter::coolCOn1=40
  integer,parameter::coolCOn2=40
  integer,parameter::coolCOn3=40
  real*8::coolCOx1(coolCOn1),coolCOx2(coolCOn2),coolCOx3(coolCOn3)
  real*8::coolCOixd1(coolCOn1-1),coolCOixd2(coolCOn2-1),coolCOixd3(coolCOn3-1)
  real*8::coolCOy(coolCOn1,coolCOn2,coolCOn3)
  real*8::coolCOx1min,coolCOx1max
  real*8::coolCOx2min,coolCOx2max
  real*8::coolCOx3min,coolCOx3max
  real*8::coolCOdvn1,coolCOdvn2,coolCOdvn3

  !xsecs from file variables

  !partition function from file
  integer,parameter::zpart_nCO=641
  integer,parameter::zpart_nH2even=2000
  integer,parameter::zpart_nH2odd=2000
  real*8::zpart_CO(zpart_nCO),minpart_CO,partdT_CO
  real*8::zpart_H2even(zpart_nH2even),minpart_H2even,partdT_H2even
  real*8::zpart_H2odd(zpart_nH2odd),minpart_H2odd,partdT_H2odd

  !Habing flux for the photoelectric heating by dust
  ! and clumping factor for H2 formation
  ! on dust by Jura/Gnedin
  real*8::GHabing,Ghabing_thin,clump_factor
  !$omp threadprivate(GHabing,GHabing_thin)

  !partition functions common vars

  !verbatim reactions
  character*50::reactionNames(nrea)

end module krome_commons

!############### MODULE ##############
module krome_constants
  implicit none

  !constants
  real*8,parameter::boltzmann_eV = 8.617332478d-5 !eV / K
  real*8,parameter::boltzmann_J = 1.380648d-23 !J / K
  real*8,parameter::boltzmann_erg = 1.380648d-16 !erg / K
  real*8,parameter::iboltzmann_eV = 1d0/boltzmann_eV !K / eV
  real*8,parameter::iboltzmann_erg = 1d0/boltzmann_erg !K / erg
  real*8,parameter::planck_eV = 4.135667516d-15 !eV s
  real*8,parameter::planck_J = 6.62606957d-34 !J s
  real*8,parameter::planck_erg = 6.62606957d-27 !erg s
  real*8,parameter::iplanck_eV = 1d0/planck_eV !1 / eV / s
  real*8,parameter::iplanck_J = 1d0/planck_J !1 / J / s
  real*8,parameter::iplanck_erg = 1d0/planck_erg !1 / erg / s
  real*8,parameter::gravity = 6.674d-8 !cm3 / g / s2
  real*8,parameter::e_mass = 9.10938188d-28 !g
  real*8,parameter::p_mass = 1.67262158d-24 !g
  real*8,parameter::n_mass = 1.674920d-24 !g
  real*8,parameter::ip_mass = 1d0/p_mass !1/g
  real*8,parameter::clight = 2.99792458e10 !cm/s
  real*8,parameter::pi = 3.14159265359d0 !#
  real*8,parameter::eV_to_erg = 1.60217646d-12 !eV -> erg
  real*8,parameter::ry_to_eV = 13.60569d0 !rydberg -> eV
  real*8,parameter::ry_to_erg = 2.179872d-11 !rydberg -> erg
  real*8,parameter::seconds_per_year = 365d0*24d0*3600d0 !yr -> s
  real*8,parameter::km_to_cm = 1d5 !km -> cm
  real*8,parameter::cm_to_Mpc = 1.d0/3.08d24 !cm -> Mpc
  real*8,parameter::kvgas_erg = 8.d0*boltzmann_erg/pi/p_mass !
  real*8,parameter::pre_kvgas_sqrt = sqrt(8.d0*boltzmann_erg/pi) !
  real*8,parameter::pre_planck = 2.d0*planck_erg/clight**2 !erg/cm2*s3
  real*8,parameter::exp_planck = planck_erg / boltzmann_erg !s*K
  real*8,parameter::stefboltz_erg = 5.670373d-5 !erg/s/cm2/K4
  real*8,parameter::N_avogadro = 6.0221d23 !#
  real*8,parameter::Rgas_J = 8.3144621d0 !J/K/mol
  real*8,parameter::Rgas_kJ = 8.3144621d-3 !kJ/K/mol
  real*8,parameter::hubble = 0.704d0 !dimensionless
  real*8,parameter::Omega0 = 1.0d0 !dimensionless
  real*8,parameter::Omegab = 0.0456d0 !dimensionless
  real*8,parameter::Hubble0 = 1.d2*hubble*km_to_cm*cm_to_Mpc !1/s

end module krome_constants

!############### MODULE ##############
module krome_fit
contains

  !*****************************
  subroutine init_anytab3D(filename,x,y,z,f,xmul,ymul,zmul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8::rout(4)
    integer::i,j,k,ios,unit

    !check the size of the X input array
    if(size(x).ne.size(f,1)) then
      print *,"ERROR: in init_anytab3D x size differs from f(x,y,z)"
      stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(f,2)) then
      print *,"ERROR: in init_anytab3D y size differs from f(x,y,z)"
      stop
    end if

    !check the size of the Z input array
    if(size(z).ne.size(f,3)) then
      print *,"ERROR: in init_anytab3D z size differs from f(x,y,z)"
      stop
    end if

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: in init_anytab3D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(unit,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
      print *,"ERROR: file "//filename//" should"
      print *," contain the number of grid points"
      print *," per dimension in the format"
      print *,"  XX, YY, ZZ"
      print *,row_string
      stop
    end if

    !loop to read file (3rd dimension of f() is
    ! first in the tables. i.e. tables are z,x,y,
    ! while f() is x,y,z
    do i=1,size(z)
      do j=1,size(x)
        do k=1,size(y)
          read(unit,*,iostat=ios) rout(:)
          y(k) = rout(3)
          f(j,k,i) = rout(4)
        end do
        x(j) = rout(2)
        read(unit,*,iostat=ios) !skip blanks
      end do
      z(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios.ne.0) exit
    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))
    zmul = 1d0/(z(2)-z(1))

  end subroutine init_anytab3D

  !********************************************
  !load 2d tables from filename
  subroutine init_anytab2D(filename,x,y,z,xmul,ymul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),z(:,:),xmul,ymul
    real*8::rout(3)
    integer::i,j,ios,unit

    !check the size of the X input array
    if(size(x).ne.size(z,1)) then
      print *,"ERROR: in init_anytab2D x size differs from z"
      stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(z,2)) then
      print *,"ERROR: in init_anytab2D y size differs from z"
      stop
    end if

    if (krome_mpi_rank<=1) print *,"Reading tables from "//trim(filename)

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: in init_anytab2D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(unit,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
      print *,"ERROR: file "//filename//" should"
      print *," contain the number of rows and "
      print *," columns in the format"
      print *,"  RR, CC"
      print *,row_string
      stop
    end if

    !loop to read file
    do i=1,size(x)
      do j=1,size(y)
        read(unit,*,iostat=ios) rout(:)
        y(j) = rout(2)
        z(i,j) = rout(3)
      end do
      x(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios.ne.0) exit
    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))

  end subroutine init_anytab2D

  !******************************
  !test 2d fit and save to file
  subroutine test_anytab2D(fname,x,y,z,xmul,ymul)
    implicit none
    integer::i,j,unit1,unit2
    real*8,intent(in)::x(:),y(:),z(:,:),xmul,ymul
    real*8::xx,yy,zz
    character(len=*),intent(in)::fname

    open(newunit=unit1,file=fname//".fit",status="replace")
    open(newunit=unit2,file=fname//".org",status="replace")
    do i=1,size(x)
      do j=1,size(y)
        xx = x(i)
        yy = y(j)
        zz = fit_anytab2D(x(:),y(:),z(:,:),xmul,ymul,xx,yy)
        write(unit1,*) xx,yy,zz
        write(unit2,*) x(i),y(j),z(i,j)
      end do
      write(unit1,*)
      write(unit2,*)
    end do
    close(unit1)
    close(unit2)
    print *,"original file wrote in ",fname//".org"
    print *,"fit test file wrote in ",fname//".fit"

  end subroutine test_anytab2D

  !*****************************
  function fit_anytab3D(x,y,z,f,xmul,ymul,zmul,xx,yy,zz)
    implicit none
    real*8,intent(in)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8,intent(in)::xx,yy,zz
    real*8::fleft(size(x),size(y)), fright(size(x),size(y))
    real*8::fit_anytab3D,fl,fr
    integer::ipos,i1,i2

    ipos = (zz-z(1)) * zmul + 1
    i1 = min(max(ipos,1), size(z)-1)
    i2 = i1 + 1
    fleft(:,:) = f(:,:,i1)
    fright(:,:) = f(:,:,i2)

    fl = fit_anytab2D(x(:), y(:), fleft(:,:), xmul, ymul, xx, yy)
    fr = fit_anytab2D(x(:), y(:), fright(:,:), xmul, ymul, xx, yy)

    fit_anytab3D = (zz-z(i1))*zmul*(fr-fl)+fl

  end function fit_anytab3D

  !******************************
  !return 2d fit at xx,yy
  function fit_anytab2D(x,y,z,xmul,ymul,xx,yy)
    implicit none
    real*8::fit_anytab2D
    real*8,intent(in)::x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D(x(:),zright(:),xmul,xx)

    fit_anytab2D = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D

  !*********************
  !return 1d fit at xx
  function fit_anytab1D(x,z,xmul,xx)
    real*8,intent(in)::x(:),z(:),xmul,xx
    real*8::fit_anytab1D,p
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    fit_anytab1D = p * (z(i2) - z(i1)) + z(i1)

  end function fit_anytab1D

  !*****************************
  function fit_anytab3D_linlinlog(x,y,z,f,xmul,ymul,zmul,xx,yy,zz)
    implicit none
    real*8,intent(in)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8,intent(in)::xx,yy,zz
    real*8::fleft(size(x),size(y)), fright(size(x),size(y))
    real*8::fit_anytab3D_linlinlog,fl,fr
    integer::ipos,i1,i2

    ipos = (zz-z(1)) * zmul + 1
    i1 = min(max(ipos,1), size(z)-1)
    i2 = i1 + 1
    fleft(:,:) = f(:,:,i1)
    fright(:,:) = f(:,:,i2)

    fl = fit_anytab2D_linlog(x(:), y(:), fleft(:,:), xmul, ymul, xx, yy)
    fr = fit_anytab2D_linlog(x(:), y(:), fright(:,:), xmul, ymul, xx, yy)

    fit_anytab3D_linlinlog = (zz-z(i1))*zmul*(fr-fl)+fl

  end function fit_anytab3D_linlinlog

  !***************************
  function fit_anytab2D_linlog(x,y,z,xmul,ymul,xx,yy)
    real*8::fit_anytab2D_linlog,x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D_linlog(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D_linlog(x(:),zright(:),xmul,xx)

    fit_anytab2D_linlog = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D_linlog

  !*********************
  function fit_anytab1D_linlog(x,z,xmul,xx)
    real*8::fit_anytab1D_linlog,x(:),z(:),xmul,xx,p,z2,z1
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    z2 = z(i2)
    z1 = z(i1)
    if(z1<0d0 .and. z2<0d0) then
      z1 = log10(-z1)
      z2 = log10(-z2)
      fit_anytab1D_linlog = -1d1**(p * (z2 - z1) + z1)
      return
    end if

    if(z1>0d0 .and. z2>0d0) then
      z1 = log10(z1)
      z2 = log10(z2)
      fit_anytab1D_linlog = 1d1**(p * (z2 - z1) + z1)
      return
    end if

    fit_anytab1D_linlog = (p * (z2 - z1) + z1)

  end function fit_anytab1D_linlog

  !*****************************
  !spline interpolation at t using array  x,y (size n) as data
  function fspline(x,y,t)
    implicit none
    real*8::fspline,x(:),y(:),b(size(x)),c(size(x)),d(size(x)),t
    integer::n

    n = size(x)
    call spline(x(:),y(:),b(:),c(:),d(:),n)
    fspline = ispline(t,x(:),y(:),b(:),c(:),d(:),n)

  end function fspline

  !*******************************+
  subroutine spline(x, y, b, c, d, n)
    !======================================================================
    !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
    !  for cubic spline interpolation
    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    !  for  x(i) <= x <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    !  input..
    !  x = the arrays of data abscissas (in strictly increasing order)
    !  y = the arrays of data ordinates
    !  n = size of the arrays xi() and yi() (n>=2)
    !  output..
    !  b, c, d  = arrays of spline coefficients
    !  comments ...
    !  spline.f90 program is based on fortran version of program spline.f
    !  the accompanying function fspline can be used for interpolation
    !======================================================================
    implicit none
    integer::n
    real*8::x(n), y(n), b(n), c(n), d(n)
    integer::i, j, gap
    real*8::h

    gap = n-1

    !check input
    if(n<2) return
    if(n<3) then
      b(1) = (y(2)-y(1))/(x(2)-x(1)) !linear interpolation
      c(1) = 0d0
      d(1) = 0d0
      b(2) = b(1)
      c(2) = 0d0
      d(2) = 0d0
      return
    end if

    !step 1: preparation
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2d0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
    end do

    ! step 2: end conditions
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0d0
    c(n) = 0d0
    if(n.ne.3) then
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if

    ! step 3: forward elimination
    do i = 2, n
      h = d(i-1)/b(i-1)
      b(i) = b(i) - h*d(i-1)
      c(i) = c(i) - h*c(i-1)
    end do

    ! step 4: back substitution
    c(n) = c(n)/b(n)
    do j = 1, gap
      i = n-j
      c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do

    ! step 5: compute spline coefficients
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2d0*c(n))
    do i = 1, gap
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2d0*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3d0*c(i)
    end do
    c(n) = 3d0*c(n)
    d(n) = d(n-1)
  end subroutine spline

  !*******************************
  function ispline(u, x, y, b, c, d, n)
    !======================================================================
    ! function ispline evaluates the cubic spline interpolation at point z
    ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
    ! where  x(i) <= u <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    ! input..
    ! u       = the abscissa at which the spline is to be evaluated
    ! x, y    = the arrays of given data points
    ! b, c, d = arrays of spline coefficients computed by spline
    ! n       = the number of data points
    ! output:
    ! ispline = interpolated value at point u
    !=======================================================================
    implicit none
    real*8::ispline
    integer::n
    real*8::u, x(n), y(n), b(n), c(n), d(n)
    integer::i, j, k
    real*8::dx

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u<=x(1)) then
      ispline = y(1)
      return
    end if

    if(u>=x(n)) then
      ispline = y(n)
      return
    end if

    ! binary search for for i, such that x(i) <= u <= x(i+1)
    i = 1
    j = n+1
    do while (j>i+1)
      k = (i+j)/2
      if(u<x(k)) then
        j=k
      else
        i=k
      end if
    end do

    ! evaluate spline interpolation
    dx = u - x(i)
    ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

  end function ispline

end module krome_fit
!This module contains useful routines to get physical
! quantities, like mean molecular weight, mass density,
! mass, jeans length, etc. etc.

!############### MODULE ##############
module krome_getphys
contains

  !*****************************
  !get the mean molecular weight in grams
  function get_mu(n)
    use krome_commons
    use krome_constants
    implicit none
    real*8::n(:),get_mu,m(nspec)
    m(:) = get_mass()

    !ip_mass is 1/proton_mass_in_g
    get_mu = max(sum(n(1:nmols)*m(1:nmols)),1d-40) &
        / max(sum(n(1:nmols)),1d-40) * ip_mass

  end function get_mu

  !***************************
  !get mean molecular weight in grams
  function get_mu_rho(n,rhogas)
    use krome_commons
    use krome_constants
    implicit none
    real*8::get_mu_rho,rhogas,n(:)

    !ip_mass is 1/proton_mass_in_g
    get_mu_rho = rhogas / max(sum(n(1:nmols)),1d-40) * ip_mass

  end function get_mu_rho

  !************************
  !get species masses (g)
  function get_mass()
    use krome_commons
    implicit none
    real*8::get_mass(nspec)

    get_mass(1) = 9.10938188d-28	!E
    get_mass(2) = 2.67691710837d-23	!O-
    get_mass(3) = 1.67444345638d-24	!H-
    get_mass(4) = 5.35374312292d-23	!S-
    get_mass(5) = 2.00771060473d-23	!C-
    get_mass(6) = 1.67353251819d-24	!H
    get_mass(7) = 6.69206503638d-24	!HE
    get_mass(8) = 2.00761951091d-23	!C
    get_mass(9) = 5.01904877728d-23	!NO
    get_mass(10) = 4.68444552546d-23	!CO
    get_mass(11) = 2.34222276273d-23	!N
    get_mass(12) = 5.3536520291d-23	!O2
    get_mass(13) = 5.3536520291d-23	!S
    get_mass(14) = 8.03047804365d-23	!SO
    get_mass(15) = 2.67682601455d-23	!O
    get_mass(16) = 3.34706503638d-24	!H2
    get_mass(17) = 4.68444552546d-23	!SI
    get_mass(18) = 2.84417926637d-23	!OH
    get_mass(19) = 5.52100528092d-23	!HS
    get_mass(20) = 7.69587479183d-23	!NS
    get_mass(21) = 5.68835853274d-23	!H2S
    get_mass(22) = 9.20143454729d-23	!FE
    get_mass(23) = 7.36127154001d-23	!CS
    get_mass(24) = 4.34984227364d-23	!CN
    get_mass(25) = 1.07073040582d-22	!S2
    get_mass(26) = 3.84788577001d-23	!NA
    get_mass(27) = 3.17867926637d-23	!F
    get_mass(28) = 3.34603251819d-23	!HF
    get_mass(29) = 2.17497276273d-23	!CH
    get_mass(30) = 1.07073040582d-22	!SO2
    get_mass(31) = 4.01523902183d-23	!C2
    get_mass(32) = 4.68444552546d-23	!N2
    get_mass(33) = 2.34232601455d-23	!CH2
    get_mass(34) = 2.50957601455d-23	!NH
    get_mass(35) = 4.51719552546d-23	!HCN
    get_mass(36) = 7.36127154001d-23	!CO2
    get_mass(37) = 7.36127154001d-23	!SIO
    get_mass(38) = 1.00380975546d-22	!SIO2
    get_mass(39) = 2.67692926637d-23	!NH2
    get_mass(40) = 7.02666828819d-23	!OCN
    get_mass(41) = 4.01523902183d-23	!MG
    get_mass(42) = 5.18629877728d-23	!P
    get_mass(43) = 4.85179877728d-23	!HCO
    get_mass(44) = 3.01153251819d-23	!H2O
    get_mass(45) = 1.00380975546d-22	!OCS
    get_mass(46) = 7.52852154001d-23	!PN
    get_mass(47) = 7.86312479183d-23	!PO
    get_mass(48) = 6.69115409819d-24	!HE+
    get_mass(49) = 1.67262158d-24	!H+
    get_mass(50) = 2.50948492073d-23	!NH+
    get_mass(51) = 5.5209141871d-23	!HS+
    get_mass(52) = 5.35356093528d-23	!S+
    get_mass(53) = 4.68435443164d-23	!SI+
    get_mass(54) = 2.84408817255d-23	!OH+
    get_mass(55) = 8.36468661638d-24	!HEH+
    get_mass(56) = 3.34615409819d-24	!H2+
    get_mass(57) = 9.20134345347d-23	!FE+
    get_mass(58) = 4.85170768346d-23	!SIH+
    get_mass(59) = 3.84779467619d-23	!NA+
    get_mass(60) = 4.85170768346d-23	!HCO+
    get_mass(61) = 2.17488166891d-23	!CH+
    get_mass(62) = 2.67673492073d-23	!O+
    get_mass(63) = 4.01514792801d-23	!MG+
    get_mass(64) = 7.36118044619d-23	!SIO+
    get_mass(65) = 5.18620768346d-23	!P+
    get_mass(66) = 7.86303369801d-23	!SIF+
    get_mass(67) = 2.00752841709d-23	!C+
    get_mass(68) = 2.34213166891d-23	!N+
    get_mass(69) = 4.68435443164d-23	!CO+
    get_mass(70) = 3.17858817255d-23	!F+
    get_mass(71) = 0.d0	!CR
    get_mass(72) = 0.d0	!g
    get_mass(73) = 0.d0	!Tgas
    get_mass(74) = 0.d0	!dummy

  end function get_mass

  !************************
  !get sqrt of the inverse of the masses (1/sqrt(g))
  function get_imass_sqrt()
    use krome_commons
    implicit none
    real*8::get_imass_sqrt(nspec)

    get_imass_sqrt(1) = 3.31326021505d+13	!E
    get_imass_sqrt(2) = 1.93278051341d+11	!O-
    get_imass_sqrt(3) = 7.72795806394d+11	!H-
    get_imass_sqrt(4) = 1.36669383456d+11	!S-
    get_imass_sqrt(5) = 2.23177004181d+11	!C-
    get_imass_sqrt(6) = 7.73006102111d+11	!H
    get_imass_sqrt(7) = 3.86562679981d+11	!HE
    get_imass_sqrt(8) = 2.23182067346d+11	!C
    get_imass_sqrt(9) = 1.41152733144d+11	!NO
    get_imass_sqrt(10) = 1.46106959624d+11	!CO
    get_imass_sqrt(11) = 2.06626443857d+11	!N
    get_imass_sqrt(12) = 1.36670546184d+11	!O2
    get_imass_sqrt(13) = 1.36670546184d+11	!S
    get_imass_sqrt(14) = 1.11591033673d+11	!SO
    get_imass_sqrt(15) = 1.93281339991d+11	!O
    get_imass_sqrt(16) = 5.46597856701d+11	!H2
    get_imass_sqrt(17) = 1.46106959624d+11	!SI
    get_imass_sqrt(18) = 1.87508740611d+11	!OH
    get_imass_sqrt(19) = 1.34583221186d+11	!HS
    get_imass_sqrt(20) = 1.13991115426d+11	!NS
    get_imass_sqrt(21) = 1.32588702018d+11	!H2S
    get_imass_sqrt(22) = 1.04249079615d+11	!FE
    get_imass_sqrt(23) = 1.16553033405d+11	!CS
    get_imass_sqrt(24) = 1.51622357573d+11	!CN
    get_imass_sqrt(25) = 96640669995.3	!S2
    get_imass_sqrt(26) = 1.61208862859d+11	!NA
    get_imass_sqrt(27) = 1.77368562161d+11	!F
    get_imass_sqrt(28) = 1.72876086d+11	!HF
    get_imass_sqrt(29) = 2.14423849574d+11	!CH
    get_imass_sqrt(30) = 96640669995.3	!SO2
    get_imass_sqrt(31) = 1.57813553259d+11	!C2
    get_imass_sqrt(32) = 1.46106959624d+11	!N2
    get_imass_sqrt(33) = 2.06621889668d+11	!CH2
    get_imass_sqrt(34) = 1.99618056318d+11	!NH
    get_imass_sqrt(35) = 1.48787194664d+11	!HCN
    get_imass_sqrt(36) = 1.16553033405d+11	!CO2
    get_imass_sqrt(37) = 1.16553033405d+11	!SIO
    get_imass_sqrt(38) = 99810054788.8	!SIO2
    get_imass_sqrt(39) = 1.93277612428d+11	!NH2
    get_imass_sqrt(40) = 1.19295832983d+11	!OCN
    get_imass_sqrt(41) = 1.57813553259d+11	!MG
    get_imass_sqrt(42) = 1.38858104892d+11	!P
    get_imass_sqrt(43) = 1.43565011358d+11	!HCO
    get_imass_sqrt(44) = 1.82224271009d+11	!H2O
    get_imass_sqrt(45) = 99810054788.8	!OCS
    get_imass_sqrt(46) = 1.15251119158d+11	!PN
    get_imass_sqrt(47) = 1.12772294263d+11	!PO
    get_imass_sqrt(48) = 3.86588992536d+11	!HE+
    get_imass_sqrt(49) = 7.732165696d+11	!H+
    get_imass_sqrt(50) = 1.99621679333d+11	!NH+
    get_imass_sqrt(51) = 1.34584331478d+11	!HS+
    get_imass_sqrt(52) = 1.36671708942d+11	!S+
    get_imass_sqrt(53) = 1.46108380244d+11	!SI+
    get_imass_sqrt(54) = 1.87511743463d+11	!OH+
    get_imass_sqrt(55) = 3.45760328884d+11	!HEH+
    get_imass_sqrt(56) = 5.46672253003d+11	!H2+
    get_imass_sqrt(57) = 1.0424959565d+11	!FE+
    get_imass_sqrt(58) = 1.43566359113d+11	!SIH+
    get_imass_sqrt(59) = 1.61210771101d+11	!NA+
    get_imass_sqrt(60) = 1.43566359113d+11	!HCO+
    get_imass_sqrt(61) = 2.14428340044d+11	!CH+
    get_imass_sqrt(62) = 1.93284628808d+11	!O+
    get_imass_sqrt(63) = 1.5781534345d+11	!MG+
    get_imass_sqrt(64) = 1.16553754568d+11	!SIO+
    get_imass_sqrt(65) = 1.38859324382d+11	!P+
    get_imass_sqrt(66) = 1.12772947499d+11	!SIF+
    get_imass_sqrt(67) = 2.23187130855d+11	!C+
    get_imass_sqrt(68) = 2.06630462037d+11	!N+
    get_imass_sqrt(69) = 1.46108380244d+11	!CO+
    get_imass_sqrt(70) = 1.77371103708d+11	!F+
    get_imass_sqrt(71) = 0.d0	!CR
    get_imass_sqrt(72) = 0.d0	!g
    get_imass_sqrt(73) = 0.d0	!Tgas
    get_imass_sqrt(74) = 0.d0	!dummy

  end function get_imass_sqrt

  !************************
  !get inverse of the species masses (1/g)
  function get_imass()
    use krome_commons
    implicit none
    real*8::get_imass(nspec)

    get_imass(1) = 1.09776932527d+27	!E
    get_imass(2) = 3.73564051301d+22	!O-
    get_imass(3) = 5.9721335838d+23	!H-
    get_imass(4) = 1.86785203743d+22	!S-
    get_imass(5) = 4.98079751954d+22	!C-
    get_imass(6) = 5.97538433901d+23	!H
    get_imass(7) = 1.49430705554d+23	!HE
    get_imass(8) = 4.98102351847d+22	!C
    get_imass(9) = 1.99240940739d+22	!NO
    get_imass(10) = 2.13472436506d+22	!CO
    get_imass(11) = 4.26944873012d+22	!N
    get_imass(12) = 1.86788381943d+22	!O2
    get_imass(13) = 1.86788381943d+22	!S
    get_imass(14) = 1.24525587962d+22	!SO
    get_imass(15) = 3.73576763885d+22	!O
    get_imass(16) = 2.9876921695d+23	!H2
    get_imass(17) = 2.13472436506d+22	!SI
    get_imass(18) = 3.51595278056d+22	!OH
    get_imass(19) = 1.81126434248d+22	!HS
    get_imass(20) = 1.2993974396d+22	!NS
    get_imass(21) = 1.75797639028d+22	!H2S
    get_imass(22) = 1.08678706006d+22	!FE
    get_imass(23) = 1.35846095958d+22	!CS
    get_imass(24) = 2.2989339316d+22	!CN
    get_imass(25) = 9.33941909714d+21	!S2
    get_imass(26) = 2.59882974644d+22	!NA
    get_imass(27) = 3.1459606843d+22	!F
    get_imass(28) = 2.98861411108d+22	!HF
    get_imass(29) = 4.59775872662d+22	!CH
    get_imass(30) = 9.33941909714d+21	!SO2
    get_imass(31) = 2.49051175924d+22	!C2
    get_imass(32) = 2.13472436506d+22	!N2
    get_imass(33) = 4.26926052901d+22	!CH2
    get_imass(34) = 3.98473684081d+22	!NH
    get_imass(35) = 2.21376292959d+22	!HCN
    get_imass(36) = 1.35846095958d+22	!CO2
    get_imass(37) = 1.35846095958d+22	!SIO
    get_imass(38) = 9.96204703694d+21	!SIO2
    get_imass(39) = 3.73562354659d+22	!NH2
    get_imass(40) = 1.42314957671d+22	!OCN
    get_imass(41) = 2.49051175924d+22	!MG
    get_imass(42) = 1.92815732942d+22	!P
    get_imass(43) = 2.06109124864d+22	!HCO
    get_imass(44) = 3.32056849448d+22	!H2O
    get_imass(45) = 9.96204703694d+21	!OCS
    get_imass(46) = 1.32828204673d+22	!PN
    get_imass(47) = 1.27175903534d+22	!PO
    get_imass(48) = 1.4945104915d+23	!HE+
    get_imass(49) = 5.97863863505d+23	!H+
    get_imass(50) = 3.98488148599d+22	!NH+
    get_imass(51) = 1.81129422793d+22	!HS+
    get_imass(52) = 1.86791560251d+22	!S+
    get_imass(53) = 2.13476587776d+22	!SI+
    get_imass(54) = 3.51606539365d+22	!OH+
    get_imass(55) = 1.1955020503d+23	!HEH+
    get_imass(56) = 2.98850552203d+23	!H2+
    get_imass(57) = 1.08679781932d+22	!FE+
    get_imass(58) = 2.0611299469d+22	!SIH+
    get_imass(59) = 2.5988912719d+22	!NA+
    get_imass(60) = 2.0611299469d+22	!HCO+
    get_imass(61) = 4.59795130141d+22	!CH+
    get_imass(62) = 3.73589477335d+22	!O+
    get_imass(63) = 2.49056826281d+22	!MG+
    get_imass(64) = 1.35847777039d+22	!SIO+
    get_imass(65) = 1.92819119679d+22	!P+
    get_imass(66) = 1.27177376876d+22	!SIF+
    get_imass(67) = 4.98124953791d+22	!C+
    get_imass(68) = 4.26961478414d+22	!N+
    get_imass(69) = 2.13476587776d+22	!CO+
    get_imass(70) = 3.14605084306d+22	!F+
    get_imass(71) = 0.d0	!CR
    get_imass(72) = 0.d0	!g
    get_imass(73) = 0.d0	!Tgas
    get_imass(74) = 0.d0	!dummy

  end function get_imass

  !************************
  !species binding energies (surface=BARE), K
  function get_EbindBare()
    use krome_commons
    implicit none
    real*8::get_EbindBare(nspec)

    get_EbindBare(:) = 1d99

    get_EbindBare(idx_H) = 500.0d0
    get_EbindBare(idx_CO) = 1100.0d0
    get_EbindBare(idx_O2) = 1250.0d0
    get_EbindBare(idx_O) = 1700.0d0
    get_EbindBare(idx_H2) = 300.0d0
    get_EbindBare(idx_OH) = 1360.0d0
    get_EbindBare(idx_CO2) = 2300.0d0
    get_EbindBare(idx_HCO) = 1100.0d0
    get_EbindBare(idx_H2O) = 4800.0d0

  end function get_EbindBare

  !************************
  !species binding energies (surface=ICE), K
  function get_EbindIce()
    use krome_commons
    implicit none
    real*8::get_EbindIce(nspec)

    get_EbindIce(:) = 1d99

    get_EbindIce(idx_H) = 650.0d0
    get_EbindIce(idx_CO) = 1300.0d0
    get_EbindIce(idx_O2) = 900.0d0
    get_EbindIce(idx_O) = 1700.0d0
    get_EbindIce(idx_H2) = 300.0d0
    get_EbindIce(idx_OH) = 3500.0d0
    get_EbindIce(idx_CO2) = 2300.0d0
    get_EbindIce(idx_HCO) = 3100.0d0
    get_EbindIce(idx_H2O) = 4800.0d0

  end function get_EbindIce

  !************************
  function get_kevap70()
    use krome_commons
    implicit none
    real*8::get_kevap70(nspec)

    get_kevap70(idx_E) = 0d0
    get_kevap70(idx_Ok) = 0d0
    get_kevap70(idx_Hk) = 0d0
    get_kevap70(idx_Sk) = 0d0
    get_kevap70(idx_Ck) = 0d0
    get_kevap70(idx_H) = 790490323.12
    get_kevap70(idx_HE) = 0d0
    get_kevap70(idx_C) = 0d0
    get_kevap70(idx_NO) = 0d0
    get_kevap70(idx_CO) = 149751.929641
    get_kevap70(idx_N) = 0d0
    get_kevap70(idx_O2) = 17568.7715065
    get_kevap70(idx_S) = 0d0
    get_kevap70(idx_SO) = 0d0
    get_kevap70(idx_O) = 28.3692788833
    get_kevap70(idx_H2) = 13763786733.1
    get_kevap70(idx_SI) = 0d0
    get_kevap70(idx_OH) = 3649.88043081
    get_kevap70(idx_HS) = 0d0
    get_kevap70(idx_NS) = 0d0
    get_kevap70(idx_H2S) = 0d0
    get_kevap70(idx_FE) = 0d0
    get_kevap70(idx_CS) = 0d0
    get_kevap70(idx_CN) = 0d0
    get_kevap70(idx_S2) = 0d0
    get_kevap70(idx_NA) = 0d0
    get_kevap70(idx_F) = 0d0
    get_kevap70(idx_HF) = 0d0
    get_kevap70(idx_CH) = 0d0
    get_kevap70(idx_SO2) = 0d0
    get_kevap70(idx_C2) = 0d0
    get_kevap70(idx_N2) = 0d0
    get_kevap70(idx_CH2) = 0d0
    get_kevap70(idx_NH) = 0d0
    get_kevap70(idx_HCN) = 0d0
    get_kevap70(idx_CO2) = 0.00537432797219
    get_kevap70(idx_SIO) = 0d0
    get_kevap70(idx_SIO2) = 0d0
    get_kevap70(idx_NH2) = 0d0
    get_kevap70(idx_OCN) = 0d0
    get_kevap70(idx_MG) = 0d0
    get_kevap70(idx_P) = 0d0
    get_kevap70(idx_HCO) = 149751.929641
    get_kevap70(idx_H2O) = 1.65884938156e-18
    get_kevap70(idx_OCS) = 0d0
    get_kevap70(idx_PN) = 0d0
    get_kevap70(idx_PO) = 0d0
    get_kevap70(idx_HEj) = 0d0
    get_kevap70(idx_Hj) = 0d0
    get_kevap70(idx_NHj) = 0d0
    get_kevap70(idx_HSj) = 0d0
    get_kevap70(idx_Sj) = 0d0
    get_kevap70(idx_SIj) = 0d0
    get_kevap70(idx_OHj) = 0d0
    get_kevap70(idx_HEHj) = 0d0
    get_kevap70(idx_H2j) = 0d0
    get_kevap70(idx_FEj) = 0d0
    get_kevap70(idx_SIHj) = 0d0
    get_kevap70(idx_NAj) = 0d0
    get_kevap70(idx_HCOj) = 0d0
    get_kevap70(idx_CHj) = 0d0
    get_kevap70(idx_Oj) = 0d0
    get_kevap70(idx_MGj) = 0d0
    get_kevap70(idx_SIOj) = 0d0
    get_kevap70(idx_Pj) = 0d0
    get_kevap70(idx_SIFj) = 0d0
    get_kevap70(idx_Cj) = 0d0
    get_kevap70(idx_Nj) = 0d0
    get_kevap70(idx_COj) = 0d0
    get_kevap70(idx_Fj) = 0d0
    get_kevap70(idx_CR) = 0d0
    get_kevap70(idx_g) = 0d0
    get_kevap70(idx_Tgas) = 0d0
    get_kevap70(idx_dummy) = 0d0

  end function get_kevap70

  !************************
  !get verbatim reaction names
  function get_rnames()
    use krome_commons
    implicit none
    character*50::get_rnames(nrea)

    !reaction names are loaded from file
    get_rnames(:) = reactionNames(:)

  end function get_rnames

  !************************
  !get species names
  function get_names()
    use krome_commons
    implicit none
    character*16::get_names(nspec)

    get_names(1) = "E"
    get_names(2) = "O-"
    get_names(3) = "H-"
    get_names(4) = "S-"
    get_names(5) = "C-"
    get_names(6) = "H"
    get_names(7) = "HE"
    get_names(8) = "C"
    get_names(9) = "NO"
    get_names(10) = "CO"
    get_names(11) = "N"
    get_names(12) = "O2"
    get_names(13) = "S"
    get_names(14) = "SO"
    get_names(15) = "O"
    get_names(16) = "H2"
    get_names(17) = "SI"
    get_names(18) = "OH"
    get_names(19) = "HS"
    get_names(20) = "NS"
    get_names(21) = "H2S"
    get_names(22) = "FE"
    get_names(23) = "CS"
    get_names(24) = "CN"
    get_names(25) = "S2"
    get_names(26) = "NA"
    get_names(27) = "F"
    get_names(28) = "HF"
    get_names(29) = "CH"
    get_names(30) = "SO2"
    get_names(31) = "C2"
    get_names(32) = "N2"
    get_names(33) = "CH2"
    get_names(34) = "NH"
    get_names(35) = "HCN"
    get_names(36) = "CO2"
    get_names(37) = "SIO"
    get_names(38) = "SIO2"
    get_names(39) = "NH2"
    get_names(40) = "OCN"
    get_names(41) = "MG"
    get_names(42) = "P"
    get_names(43) = "HCO"
    get_names(44) = "H2O"
    get_names(45) = "OCS"
    get_names(46) = "PN"
    get_names(47) = "PO"
    get_names(48) = "HE+"
    get_names(49) = "H+"
    get_names(50) = "NH+"
    get_names(51) = "HS+"
    get_names(52) = "S+"
    get_names(53) = "SI+"
    get_names(54) = "OH+"
    get_names(55) = "HEH+"
    get_names(56) = "H2+"
    get_names(57) = "FE+"
    get_names(58) = "SIH+"
    get_names(59) = "NA+"
    get_names(60) = "HCO+"
    get_names(61) = "CH+"
    get_names(62) = "O+"
    get_names(63) = "MG+"
    get_names(64) = "SIO+"
    get_names(65) = "P+"
    get_names(66) = "SIF+"
    get_names(67) = "C+"
    get_names(68) = "N+"
    get_names(69) = "CO+"
    get_names(70) = "F+"
    get_names(71) = "CR"
    get_names(72) = "g"
    get_names(73) = "Tgas"
    get_names(74) = "dummy"

  end function get_names

  !************************
  !get cooling names list (empty element if cooling not present)
  function get_cooling_names()
    use krome_commons
    implicit none
    character*16::get_cooling_names(ncools)

    get_cooling_names(:) = ""

    get_cooling_names(idx_cool_h2) = "H2"
    get_cooling_names(idx_cool_h2gp) = "H2GP"
    get_cooling_names(idx_cool_atomic) = "ATOMIC"
    get_cooling_names(idx_cool_cen) = "CEN"
    get_cooling_names(idx_cool_hd) = "HD"
    get_cooling_names(idx_cool_metal) = "METAL"
    get_cooling_names(idx_cool_z) = "Z"
    get_cooling_names(idx_cool_dh) = "DH"
    get_cooling_names(idx_cool_enthalpic) = "ENTHALPIC"
    get_cooling_names(idx_cool_dust) = "DUST"
    get_cooling_names(idx_cool_compton) = "COMPTON"
    get_cooling_names(idx_cool_cie) = "CIE"
    get_cooling_names(idx_cool_cont) = "CONT"
    get_cooling_names(idx_cool_continuum) = "CONTINUUM"
    get_cooling_names(idx_cool_expansion) = "EXPANSION"
    get_cooling_names(idx_cool_exp) = "EXP"
    get_cooling_names(idx_cool_ff) = "FF"
    get_cooling_names(idx_cool_bss) = "BSS"
    get_cooling_names(idx_cool_custom) = "CUSTOM"
    get_cooling_names(idx_cool_co) = "CO"
    get_cooling_names(idx_cool_zcie) = "ZCIE"
    get_cooling_names(idx_cool_zcienouv) = "ZCIENOUV"
    get_cooling_names(idx_cool_zextend) = "ZEXTEND"
    get_cooling_names(idx_cool_gh) = "GH"

  end function get_cooling_names

  !************************
  !get heating names list (empty element if heating not present)
  function get_heating_names()
    use krome_commons
    implicit none
    character*16::get_heating_names(nheats)

    get_heating_names(:) = ""

    get_heating_names(idx_heat_chem) = "CHEM"
    get_heating_names(idx_heat_compress) = "COMPRESS"
    get_heating_names(idx_heat_compr) = "COMPR"
    get_heating_names(idx_heat_photo) = "PHOTO"
    get_heating_names(idx_heat_dh) = "DH"
    get_heating_names(idx_heat_enthalpic) = "ENTHALPIC"
    get_heating_names(idx_heat_av) = "AV"
    get_heating_names(idx_heat_photoav) = "PHOTOAV"
    get_heating_names(idx_heat_cr) = "CR"
    get_heating_names(idx_heat_dust) = "DUST"
    get_heating_names(idx_heat_xray) = "XRAY"
    get_heating_names(idx_heat_viscous) = "VISCOUS"
    get_heating_names(idx_heat_visc) = "VISC"
    get_heating_names(idx_heat_custom) = "CUSTOM"
    get_heating_names(idx_heat_zcie) = "ZCIE"

  end function get_heating_names

  !******************************
  !get the total number of H nuclei
  function get_Hnuclei(n)
    use krome_commons
    real*8::n(:),get_Hnuclei,nH

    nH = n(idx_Hk) + &
        n(idx_H) + &
        n(idx_H2)*2d0 + &
        n(idx_OH) + &
        n(idx_HS) + &
        n(idx_H2S)*2d0 + &
        n(idx_HF) + &
        n(idx_CH) + &
        n(idx_CH2)*2d0 + &
        n(idx_NH) + &
        n(idx_HCN) + &
        n(idx_NH2)*2d0 + &
        n(idx_HCO) + &
        n(idx_H2O)*2d0 + &
        n(idx_Hj) + &
        n(idx_NHj) + &
        n(idx_HSj) + &
        n(idx_OHj) + &
        n(idx_HEHj) + &
        n(idx_H2j)*2d0 + &
        n(idx_SIHj) + &
        n(idx_HCOj) + &
        n(idx_CHj)
    get_Hnuclei = nH

  end function get_Hnuclei

  !***************************
  function get_zatoms()
    use krome_commons
    implicit none
    integer::get_zatoms(nspec)

    get_zatoms(1) = 0	!E
    get_zatoms(2) = 8	!O-
    get_zatoms(3) = 1	!H-
    get_zatoms(4) = 16	!S-
    get_zatoms(5) = 6	!C-
    get_zatoms(6) = 1	!H
    get_zatoms(7) = 2	!HE
    get_zatoms(8) = 6	!C
    get_zatoms(9) = 15	!NO
    get_zatoms(10) = 14	!CO
    get_zatoms(11) = 7	!N
    get_zatoms(12) = 16	!O2
    get_zatoms(13) = 16	!S
    get_zatoms(14) = 24	!SO
    get_zatoms(15) = 8	!O
    get_zatoms(16) = 2	!H2
    get_zatoms(17) = 14	!SI
    get_zatoms(18) = 9	!OH
    get_zatoms(19) = 17	!HS
    get_zatoms(20) = 23	!NS
    get_zatoms(21) = 18	!H2S
    get_zatoms(22) = 26	!FE
    get_zatoms(23) = 22	!CS
    get_zatoms(24) = 13	!CN
    get_zatoms(25) = 32	!S2
    get_zatoms(26) = 11	!NA
    get_zatoms(27) = 9	!F
    get_zatoms(28) = 10	!HF
    get_zatoms(29) = 7	!CH
    get_zatoms(30) = 32	!SO2
    get_zatoms(31) = 12	!C2
    get_zatoms(32) = 14	!N2
    get_zatoms(33) = 8	!CH2
    get_zatoms(34) = 8	!NH
    get_zatoms(35) = 14	!HCN
    get_zatoms(36) = 22	!CO2
    get_zatoms(37) = 22	!SIO
    get_zatoms(38) = 30	!SIO2
    get_zatoms(39) = 9	!NH2
    get_zatoms(40) = 21	!OCN
    get_zatoms(41) = 12	!MG
    get_zatoms(42) = 15	!P
    get_zatoms(43) = 15	!HCO
    get_zatoms(44) = 10	!H2O
    get_zatoms(45) = 30	!OCS
    get_zatoms(46) = 22	!PN
    get_zatoms(47) = 23	!PO
    get_zatoms(48) = 2	!HE+
    get_zatoms(49) = 1	!H+
    get_zatoms(50) = 8	!NH+
    get_zatoms(51) = 17	!HS+
    get_zatoms(52) = 16	!S+
    get_zatoms(53) = 14	!SI+
    get_zatoms(54) = 9	!OH+
    get_zatoms(55) = 3	!HEH+
    get_zatoms(56) = 2	!H2+
    get_zatoms(57) = 26	!FE+
    get_zatoms(58) = 15	!SIH+
    get_zatoms(59) = 11	!NA+
    get_zatoms(60) = 15	!HCO+
    get_zatoms(61) = 7	!CH+
    get_zatoms(62) = 8	!O+
    get_zatoms(63) = 12	!MG+
    get_zatoms(64) = 22	!SIO+
    get_zatoms(65) = 15	!P+
    get_zatoms(66) = 23	!SIF+
    get_zatoms(67) = 6	!C+
    get_zatoms(68) = 7	!N+
    get_zatoms(69) = 14	!CO+
    get_zatoms(70) = 9	!F+
    get_zatoms(71) = 0	!CR
    get_zatoms(72) = 0	!g
    get_zatoms(73) = 0	!Tgas
    get_zatoms(74) = 0	!dummy

  end function get_zatoms

  !******************************
  function get_qeff()
    use krome_commons
    implicit none
    real*8::get_qeff(nrea)

    get_qeff(:) = 0e0

  end function get_qeff

  !**************************
  function get_free_fall_time(n)
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),m(nspec)
    real*8::rhogas,get_free_fall_time

    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols))
    get_free_fall_time = sqrt(3d0*pi/32d0/gravity/rhogas)

  end function get_free_fall_time

  !**************************
  function get_free_fall_time_rho(rhogas)
    use krome_constants
    implicit none
    real*8::rhogas,get_free_fall_time_rho

    get_free_fall_time_rho = sqrt(3d0*pi/32d0/gravity/rhogas)

  end function get_free_fall_time_rho

  !********************************
  function get_jeans_length(n,Tgas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::m(nspec),get_jeans_length
    m(:) = get_mass()
    rhogas = max(sum(n(1:nmols)*m(1:nmols)),1d-40)
    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length

  !********************************
  function get_jeans_length_rho(n,Tgas,rhogas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::get_jeans_length_rho

    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length_rho = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length_rho

  !***************************
  !number density to column density conversion
  function num2col(ncalc,n)
    use krome_commons
    implicit none
    real*8::num2col,ncalc,n(:),Tgas
    Tgas = max(n(idx_Tgas),phys_Tcmb)

    num2col = 1.87d21*(max(ncalc,1d-40)*1d-3)**(2./3.)

  end function num2col

  !***********************
  !column density to number density conversion
  function col2num(ncalc,n)
    use krome_commons
    implicit none
    real*8::col2num,ncalc,n(:),Tgas
    Tgas = max(n(idx_Tgas),phys_Tcmb)

    col2num = 1d3 * (max(ncalc,1d-40)/1.87d21)**1.5

  end function col2num

  !************************
  !get electrons by balancing charges
  function get_electrons(n)
    use krome_commons
    implicit none
    real*8::get_electrons,n(nspec)

    get_electrons =  - n(idx_Ok) &
        - n(idx_Hk) &
        - n(idx_Sk) &
        - n(idx_Ck) &
        + n(idx_HEj) &
        + n(idx_Hj) &
        + n(idx_NHj) &
        + n(idx_HSj) &
        + n(idx_Sj) &
        + n(idx_SIj) &
        + n(idx_OHj) &
        + n(idx_HEHj) &
        + n(idx_H2j) &
        + n(idx_FEj) &
        + n(idx_SIHj) &
        + n(idx_NAj) &
        + n(idx_HCOj) &
        + n(idx_CHj) &
        + n(idx_Oj) &
        + n(idx_MGj) &
        + n(idx_SIOj) &
        + n(idx_Pj) &
        + n(idx_SIFj) &
        + n(idx_Cj) &
        + n(idx_Nj) &
        + n(idx_COj) &
        + n(idx_Fj)
    get_electrons = max(get_electrons,0d0)

  end function get_electrons

  !************************
  !get species charges
  function get_charges()
    use krome_commons
    implicit none
    integer::get_charges(nspec)

    get_charges(1) = -1.d0 	!E
    get_charges(2) = -1.d0 	!O-
    get_charges(3) = -1.d0 	!H-
    get_charges(4) = -1.d0 	!S-
    get_charges(5) = -1.d0 	!C-
    get_charges(6) = 0.d0 	!H
    get_charges(7) = 0.d0 	!HE
    get_charges(8) = 0.d0 	!C
    get_charges(9) = 0.d0 	!NO
    get_charges(10) = 0.d0 	!CO
    get_charges(11) = 0.d0 	!N
    get_charges(12) = 0.d0 	!O2
    get_charges(13) = 0.d0 	!S
    get_charges(14) = 0.d0 	!SO
    get_charges(15) = 0.d0 	!O
    get_charges(16) = 0.d0 	!H2
    get_charges(17) = 0.d0 	!SI
    get_charges(18) = 0.d0 	!OH
    get_charges(19) = 0.d0 	!HS
    get_charges(20) = 0.d0 	!NS
    get_charges(21) = 0.d0 	!H2S
    get_charges(22) = 0.d0 	!FE
    get_charges(23) = 0.d0 	!CS
    get_charges(24) = 0.d0 	!CN
    get_charges(25) = 0.d0 	!S2
    get_charges(26) = 0.d0 	!NA
    get_charges(27) = 0.d0 	!F
    get_charges(28) = 0.d0 	!HF
    get_charges(29) = 0.d0 	!CH
    get_charges(30) = 0.d0 	!SO2
    get_charges(31) = 0.d0 	!C2
    get_charges(32) = 0.d0 	!N2
    get_charges(33) = 0.d0 	!CH2
    get_charges(34) = 0.d0 	!NH
    get_charges(35) = 0.d0 	!HCN
    get_charges(36) = 0.d0 	!CO2
    get_charges(37) = 0.d0 	!SIO
    get_charges(38) = 0.d0 	!SIO2
    get_charges(39) = 0.d0 	!NH2
    get_charges(40) = 0.d0 	!OCN
    get_charges(41) = 0.d0 	!MG
    get_charges(42) = 0.d0 	!P
    get_charges(43) = 0.d0 	!HCO
    get_charges(44) = 0.d0 	!H2O
    get_charges(45) = 0.d0 	!OCS
    get_charges(46) = 0.d0 	!PN
    get_charges(47) = 0.d0 	!PO
    get_charges(48) = 1.d0 	!HE+
    get_charges(49) = 1.d0 	!H+
    get_charges(50) = 1.d0 	!NH+
    get_charges(51) = 1.d0 	!HS+
    get_charges(52) = 1.d0 	!S+
    get_charges(53) = 1.d0 	!SI+
    get_charges(54) = 1.d0 	!OH+
    get_charges(55) = 1.d0 	!HEH+
    get_charges(56) = 1.d0 	!H2+
    get_charges(57) = 1.d0 	!FE+
    get_charges(58) = 1.d0 	!SIH+
    get_charges(59) = 1.d0 	!NA+
    get_charges(60) = 1.d0 	!HCO+
    get_charges(61) = 1.d0 	!CH+
    get_charges(62) = 1.d0 	!O+
    get_charges(63) = 1.d0 	!MG+
    get_charges(64) = 1.d0 	!SIO+
    get_charges(65) = 1.d0 	!P+
    get_charges(66) = 1.d0 	!SIF+
    get_charges(67) = 1.d0 	!C+
    get_charges(68) = 1.d0 	!N+
    get_charges(69) = 1.d0 	!CO+
    get_charges(70) = 1.d0 	!F+
    get_charges(71) = 0.d0 	!CR
    get_charges(72) = 0.d0 	!g
    get_charges(73) = 0.d0 	!Tgas
    get_charges(74) = 0.d0 	!dummy

  end function get_charges

  !*****************************
  ! get metallicity using C as reference
  function get_metallicityC(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityC,zC,nH

    nH = get_Hnuclei(n(:))

    zC = n(idx_Ck) &
        + n(idx_C) &
        + n(idx_CO) &
        + n(idx_CS) &
        + n(idx_CN) &
        + n(idx_CH) &
        + 2d0*n(idx_C2) &
        + n(idx_CH2) &
        + n(idx_HCN) &
        + n(idx_CO2) &
        + n(idx_OCN) &
        + n(idx_HCO) &
        + n(idx_OCS) &
        + n(idx_HCOj) &
        + n(idx_CHj) &
        + n(idx_Cj) &
        + n(idx_COj)

    zC = max(zC, 0d0)

    get_metallicityC = log10(zC/nH+1d-40) - (-3.57)

    phys_metallicity = get_metallicityC

  end function get_metallicityC

  !*****************************
  ! get metallicity using F as reference
  function get_metallicityF(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityF,zF,nH

    nH = get_Hnuclei(n(:))

    zF = n(idx_F) &
        + n(idx_HF) &
        + n(idx_SIFj) &
        + n(idx_Fj)

    zF = max(zF, 0d0)

    get_metallicityF = log10(zF/nH+1d-40) - (-7.44)

    phys_metallicity = get_metallicityF

  end function get_metallicityF

  !*****************************
  ! get metallicity using Mg as reference
  function get_metallicityMg(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityMg,zMg,nH

    nH = get_Hnuclei(n(:))

    zMg = n(idx_MG) &
        + n(idx_MGj)

    zMg = max(zMg, 0d0)

    get_metallicityMg = log10(zMg/nH+1d-40) - (-4.4)

    phys_metallicity = get_metallicityMg

  end function get_metallicityMg

  !*****************************
  ! get metallicity using Si as reference
  function get_metallicitySi(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicitySi,zSi,nH

    nH = get_Hnuclei(n(:))

    zSi = n(idx_SI) &
        + n(idx_SIO) &
        + n(idx_SIO2) &
        + n(idx_SIj) &
        + n(idx_SIHj) &
        + n(idx_SIOj) &
        + n(idx_SIFj)

    zSi = max(zSi, 0d0)

    get_metallicitySi = log10(zSi/nH+1d-40) - (-4.49)

    phys_metallicity = get_metallicitySi

  end function get_metallicitySi

  !*****************************
  ! get metallicity using O as reference
  function get_metallicityO(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityO,zO,nH

    nH = get_Hnuclei(n(:))

    zO = n(idx_Ok) &
        + n(idx_NO) &
        + n(idx_CO) &
        + 2d0*n(idx_O2) &
        + n(idx_SO) &
        + n(idx_O) &
        + n(idx_OH) &
        + 2d0*n(idx_SO2) &
        + 2d0*n(idx_CO2) &
        + n(idx_SIO) &
        + 2d0*n(idx_SIO2) &
        + n(idx_OCN) &
        + n(idx_HCO) &
        + n(idx_H2O) &
        + n(idx_OCS) &
        + n(idx_PO) &
        + n(idx_OHj) &
        + n(idx_HCOj) &
        + n(idx_Oj) &
        + n(idx_SIOj) &
        + n(idx_COj)

    zO = max(zO, 0d0)

    get_metallicityO = log10(zO/nH+1d-40) - (-3.31)

    phys_metallicity = get_metallicityO

  end function get_metallicityO

  !*****************************
  ! get metallicity using N as reference
  function get_metallicityN(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityN,zN,nH

    nH = get_Hnuclei(n(:))

    zN = n(idx_NO) &
        + n(idx_N) &
        + n(idx_NS) &
        + n(idx_CN) &
        + 2d0*n(idx_N2) &
        + n(idx_NH) &
        + n(idx_HCN) &
        + n(idx_NH2) &
        + n(idx_OCN) &
        + n(idx_PN) &
        + n(idx_NHj) &
        + n(idx_Nj)

    zN = max(zN, 0d0)

    get_metallicityN = log10(zN/nH+1d-40) - (-4.17)

    phys_metallicity = get_metallicityN

  end function get_metallicityN

  !*****************************
  ! get metallicity using P as reference
  function get_metallicityP(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityP,zP,nH

    nH = get_Hnuclei(n(:))

    zP = n(idx_P) &
        + n(idx_PN) &
        + n(idx_PO) &
        + n(idx_Pj)

    zP = max(zP, 0d0)

    get_metallicityP = log10(zP/nH+1d-40) - (-6.59)

    phys_metallicity = get_metallicityP

  end function get_metallicityP

  !*****************************
  ! get metallicity using S as reference
  function get_metallicityS(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityS,zS,nH

    nH = get_Hnuclei(n(:))

    zS = n(idx_Sk) &
        + n(idx_S) &
        + n(idx_SO) &
        + n(idx_HS) &
        + n(idx_NS) &
        + n(idx_H2S) &
        + n(idx_CS) &
        + 2d0*n(idx_S2) &
        + n(idx_SO2) &
        + n(idx_OCS) &
        + n(idx_HSj) &
        + n(idx_Sj)

    zS = max(zS, 0d0)

    get_metallicityS = log10(zS/nH+1d-40) - (-4.88)

    phys_metallicity = get_metallicityS

  end function get_metallicityS

  !*****************************
  ! get metallicity using Fe as reference
  function get_metallicityFe(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityFe,zFe,nH

    nH = get_Hnuclei(n(:))

    zFe = n(idx_FE) &
        + n(idx_FEj)

    zFe = max(zFe, 0d0)

    get_metallicityFe = log10(zFe/nH+1d-40) - (-4.5)

    phys_metallicity = get_metallicityFe

  end function get_metallicityFe

  !*****************************
  ! get metallicity using Na as reference
  function get_metallicityNa(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityNa,zNa,nH

    nH = get_Hnuclei(n(:))

    zNa = n(idx_NA) &
        + n(idx_NAj)

    zNa = max(zNa, 0d0)

    get_metallicityNa = log10(zNa/nH+1d-40) - (-5.76)

    phys_metallicity = get_metallicityNa

  end function get_metallicityNa

end module krome_getphys
!This module contains the functions and subroutines
! needed to evaluate the adiabatic index.

!############### MODULE ##############
module krome_gadiab
contains

  !#KROME_header

  !**************************
  !compute 1/(gamma-1) at Tgasin using the partition function
  ! provided in the array_part with a temperature step dT_part
  ! and a minimum Tgas value min_part
  function gamma_pop(array_part,dT_part,min_part,Tgasin)
    implicit none
    real*8::array_part(:),dT_part
    real*8::min_part,Tgas,gamma_pop,Tgas2,Tgasin
    real*8::logz,logz1,logz2,emed1,emed2,Cv,inTgas,T2,T1,Cv1,Cv2
    integer::idx

    !temperature above minimum data point
    inTgas = max(Tgasin,min_part)

    !data index
    idx = (inTgas-min_part)/dT_part+1
    !corresponding Tgas
    Tgas = (idx-1)*dT_part+min_part
    !store Tgas
    T1 = Tgas

    !ln of partition functions (3 points forward)
    logz = log(array_part(idx))
    logz1 = log(array_part(idx+1))
    logz2 = log(array_part(idx+2))

    !derivative for mean energy (2 points forward)
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !derivative for 1/(gamma-1)
    Cv1 = (emed2-emed1)/dT_part

    !next point temperature
    Tgas = (idx)*dT_part+min_part
    !store Tgas
    T2 = Tgas
    !ln of partition functions
    logz = logz1
    logz1 = logz2
    logz2 = log(array_part(idx+3))

    !derivative for mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !derivative for 1/(gamma-1)
    Cv2 = (emed2-emed1)/dT_part

    !interpolation for 1/(gamma-1)
    Cv = (Cv2-Cv1)*(inTgas-T1)/(T2-T1)+Cv1

    !returns result
    gamma_pop = Cv

  end function gamma_pop

  !*****************************
  !compute 1/(gamma-1) at Tgasin using the partition function
  ! provided in the array_part with a temperature step dT_part
  ! and a minimum Tgas value min_part, for H2 with a ortho/para
  ! ratio of opratio. Needs even and odd partition functions.
  function gamma_pop_H2(array_part_even,array_part_odd,dT_part,&
        min_part,Tgasin,opratio)
    implicit none
    real*8::array_part_even(:),array_part_odd(:),dT_part,zcut(4)
    real*8::min_part,Tgas,opratio,gamma_pop_H2,Tgas2,a,b,Tgasin
    real*8::logz,logz1,logz2,emed1,emed2,Cv,inTgas,T2,T1,Cv1,Cv2
    integer::idx

    !Tgas above the data limit
    inTgas = max(Tgasin,min_part)

    !exponents for ortho/para ratio
    a = opratio/(opratio+1d0) !exponent zo
    b = 1d0-a !exponent zp

    !index in the data for the given Tgas
    idx = (inTgas-min_part)/dT_part+1
    !get the corresponding Tgas
    Tgas = (idx-1)*dT_part+min_part
    !store Tgas
    T1 = Tgas

    !needed for ortho partition function (see Boley+2007)
    zcut(1) = exp(2d0*85.4/Tgas)
    zcut(2) = exp(2d0*85.4/(Tgas+dT_part))
    zcut(3) = exp(2d0*85.4/(Tgas+2d0*dT_part))
    zcut(4) = exp(2d0*85.4/(Tgas+3d0*dT_part))

    !ln of the composite partition function
    logz = log(array_part_even(idx)**b*(3d0*array_part_odd(idx)*zcut(1))**a)
    logz1 = log(array_part_even(idx+1)**b*(3d0*array_part_odd(idx+1)*zcut(2))**a)
    logz2 = log(array_part_even(idx+2)**b*(3d0*array_part_odd(idx+2)*zcut(3))**a)
    !derivative for mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !get 1/(gamma-1) for the left point
    Cv1 = (emed2-emed1)/dT_part

    !Tgas of the right point
    Tgas = (idx)*dT_part+min_part
    !store Tgas
    T2 = Tgas
    !ln of the composite function
    logz = logz1
    logz1 = logz2
    logz2 = log(array_part_even(idx+3)**b*(3d0*array_part_odd(idx+3)*zcut(4))**a)
    !derivative for the mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !get 1/(gamma-1) for the right point
    Cv2 = (emed2-emed1)/dT_part

    !interpolation of 1/(gamma-1)
    Cv = (Cv2-Cv1)*(inTgas-T1)/(T2-T1)+Cv1

    !returns the result
    gamma_pop_H2 = Cv
  end function gamma_pop_H2

  !**************************
  !function to get the partition function
  ! of H2 at Tgas with a orto-para ratio
  ! equal to opratio
  function zfop(Tgas,opratio)
    implicit none
    real*8::Tgas,zfop,brot,ibTgas
    real*8::a,b,zo,zp,opratio
    integer::j,jmax,j1
    brot = 85.4d0 !H2 rotational constant in K
    zo = 0d0 !sum for ortho partition function
    zp = 0d0 !sum for para partition function
    jmax = 10 !number of terms in sum

    ibTgas = brot/Tgas !pre-calc

    !loop over levels
    do j=0,jmax,2 !step 2
      j1 = j + 1
      zp = zp + (2d0*j+1d0) * exp(-j*(j+1d0)*ibTgas)
      zo = zo + 3d0 * (2d0*j1+1d0) * exp(-j1*(j1+1d0)*ibTgas)
    end do

    a = opratio/(opratio+1d0) !exponent zo
    b = 1d0-a !exponent zp

    zfop = (zp**b * zo**a*exp(-2d0*ibTgas)) !final partition f

  end function zfop

  !*********************
  !get the partition function at Tgas
  ! of a diatom with rotational constant
  ! brot in K
  function zf(Tgas,brot)
    real*8::Tgas,zf,brot,z,ibTgas
    integer::j,jmax
    jmax = 10 !number of levels

    ibTgas = brot/Tgas !store
    z = 0d0
    !loop on levels
    do j=0,jmax
      z = z + (2d0*j+1d0)*exp(-j*(j+1d0)*ibTgas)
    end do

    zf = z

  end function zf

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of H2 with
  ! an ortho-para ratio of opratio
  function gamma_rotop(Tgas_in,opratio)
    implicit none
    real*8::gamma_rotop,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,opratio

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zfop(Tgas+dT,opratio)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zfop(Tgas,opratio)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zfop(Tgas+dT+dT,opratio))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rotop = (prot2-prot1)*idT

  end function gamma_rotop

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of a diatom
  ! with rotational constant brot in K
  function gamma_rot(Tgas_in,brot)
    implicit none
    real*8::gamma_rot,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,brot

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zf(Tgas+dT,brot)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zf(Tgas,brot)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zf(Tgas+dT+dT,brot))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rot = (prot2-prot1)*idT

  end function gamma_rot

  !*********************
  !get gamma
  function gamma_index(n)
    use krome_commons
    implicit none
    real*8::n(:),gamma_index,krome_gamma

    real*8::Tgas,invTgas,x,expx,ysum,gsum,mosum,gvib
    real*8::Tgas_vib,invTgas_vib
    real*8::gi_NO,gi_CO,gi_O2,gi_H2,gi_OH,gi_CH,gi_C2,gi_N2,gi_NH,gi_H2j,gi_COj

    !avoid small Tgas that causes large x=a/Tgas below
    Tgas_vib = max(n(idx_Tgas), 22.5771823519d0)
    Tgas = n(idx_Tgas)
    invTgas_vib = 1d0/Tgas_vib

    !evaluate 1/(gamma-1) for NO
    x = 2720.60847513d0*invTgas_vib
    expx = exp(x)
    gvib = 2d0*x*x*expx/(expx-1d0)**2
    gi_NO = 0.5d0*(3d0 + 2d0 + gvib)

    !evaluate 1/(gamma-1) for CO
    x = 3100.12551807d0*invTgas_vib
    expx = exp(x)
    gvib = 2d0*x*x*expx/(expx-1d0)**2
    gi_CO = 0.5d0*(3d0 + 2d0 + gvib)

    !evaluate 1/(gamma-1) for O2
    x = 2257.71823519d0*invTgas_vib
    expx = exp(x)
    gvib = 2d0*x*x*expx/(expx-1d0)**2
    gi_O2 = 0.5d0*(3d0 + 2d0 + gvib)

    !evaluate 1/(gamma-1) for H2
    x = 6288.40912227d0*invTgas_vib
    expx = exp(x)
    gvib = 2d0*x*x*expx/(expx-1d0)**2
    gi_H2 = 0.5d0*(3d0 + 2d0 + gvib)

    !evaluate 1/(gamma-1) for OH
    x = 5340.47553919d0*invTgas_vib
    expx = exp(x)
    gvib = 2d0*x*x*expx/(expx-1d0)**2
    gi_OH = 0.5d0*(3d0 + 2d0 + gvib)

    !evaluate 1/(gamma-1) for CH
    x = 4087.41213553d0*invTgas_vib
    expx = exp(x)
    gvib = 2d0*x*x*expx/(expx-1d0)**2
    gi_CH = 0.5d0*(3d0 + 2d0 + gvib)

    !evaluate 1/(gamma-1) for C2
    x = 2650.50017878d0*invTgas_vib
    expx = exp(x)
    gvib = 2d0*x*x*expx/(expx-1d0)**2
    gi_C2 = 0.5d0*(3d0 + 2d0 + gvib)

    !evaluate 1/(gamma-1) for N2
    x = 3369.9012303d0*invTgas_vib
    expx = exp(x)
    gvib = 2d0*x*x*expx/(expx-1d0)**2
    gi_N2 = 0.5d0*(3d0 + 2d0 + gvib)

    !evaluate 1/(gamma-1) for NH
    x = 4690.3175088d0*invTgas_vib
    expx = exp(x)
    gvib = 2d0*x*x*expx/(expx-1d0)**2
    gi_NH = 0.5d0*(3d0 + 2d0 + gvib)

    !evaluate 1/(gamma-1) for H2+
    gi_H2j = 0.5d0*(3d0 + 2d0 + 0d0)

    !evaluate 1/(gamma-1) for CO+
    x = 3163.52251633d0*invTgas_vib
    expx = exp(x)
    gvib = 2d0*x*x*expx/(expx-1d0)**2
    gi_COj = 0.5d0*(3d0 + 2d0 + gvib)

    !sum monotomic abundances
    mosum = n(idx_E) + &
        n(idx_Ok) + &
        n(idx_Hk) + &
        n(idx_Sk) + &
        n(idx_Ck) + &
        n(idx_H) + &
        n(idx_HE) + &
        n(idx_C) + &
        n(idx_N) + &
        n(idx_S) + &
        n(idx_O) + &
        n(idx_SI) + &
        n(idx_FE) + &
        n(idx_NA) + &
        n(idx_F) + &
        n(idx_MG) + &
        n(idx_P) + &
        n(idx_HEj) + &
        n(idx_Hj) + &
        n(idx_Sj) + &
        n(idx_SIj) + &
        n(idx_FEj) + &
        n(idx_NAj) + &
        n(idx_Oj) + &
        n(idx_MGj) + &
        n(idx_Pj) + &
        n(idx_Cj) + &
        n(idx_Nj) + &
        n(idx_Fj)

    !sum all abundances
    ysum = mosum + n(idx_NO) + &
        n(idx_CO) + &
        n(idx_O2) + &
        n(idx_H2) + &
        n(idx_OH) + &
        n(idx_CH) + &
        n(idx_C2) + &
        n(idx_N2) + &
        n(idx_NH) + &
        n(idx_H2j) + &
        n(idx_COj)

    !computes gamma
    gsum = mosum * 1.5d0 + n(idx_NO)*gi_NO + &
        n(idx_CO)*gi_CO + &
        n(idx_O2)*gi_O2 + &
        n(idx_H2)*gi_H2 + &
        n(idx_OH)*gi_OH + &
        n(idx_CH)*gi_CH + &
        n(idx_C2)*gi_C2 + &
        n(idx_N2)*gi_N2 + &
        n(idx_NH)*gi_NH + &
        n(idx_H2j)*gi_H2j + &
        n(idx_COj)*gi_COj
    krome_gamma = 1d0 + ysum/gsum

    gamma_index = krome_gamma
  end function gamma_index

end module krome_gadiab
!This module contains functions and subroutines
! for the surface chemistry, including adsorption, desorption, chemisorption
! and icy grains.

!############### MODULE ##############
module krome_grfuncs
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2018-01-25 18:47:28
  !  Changeset 9a6e335
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !**********************
  !get Tdust from tables, K
  function get_table_Tdust(n) result(Tdust)
    use krome_commons
    use krome_fit
    implicit none
    real*8,intent(in)::n(nspec)
    real*8::ntot,Tdust,Tgas

    Tgas = n(idx_Tgas)

    !default, K
    Tdust = 1d0

    !total densitym, cm-3
    ntot = sum(n(1:nmols))

    !zero density returns default
    if(ntot==0d0) return

    !get dust temperature from table, K
    Tdust = 1d1**fit_anytab2D(dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_Tdust(:,:), dust_mult_ngas, &
        dust_mult_Tgas, &
        log10(ntot), log10(Tgas))

  end function get_table_Tdust

  !**********************
  !adsorpion rate Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014
  function dust_adsorption_rate(nndust,ims,stick,adust2,sqrTgas)
    use krome_constants
    implicit none
    real*8::dust_adsorption_rate,nndust,ims,stick,adust2,sqrTgas

    dust_adsorption_rate = nndust * pi * adust2 &
        * pre_kvgas_sqrt * ims * sqrTgas &
        * stick

  end function dust_adsorption_rate

  !*****************************
  !desorption rate Cazaux+2010, Hocuk+2014
  function dust_desorption_rate(fice,expEice,expEbare)
    implicit none
    real*8::dust_desorption_rate
    real*8::fice,expEice,expEbare,nu0,fbare

    nu0 = 1d12 !1/s
    fbare = 1d0 - fice
    dust_desorption_rate = nu0 * (fbare * expEbare &
        + fice * expEice)

  end function dust_desorption_rate

  !**************************
  function dust_2body_rate(p,invphi,fice,expEice1,expEice2,&
        expEbare1,expEbare2,pesc_ice,pesc_bare)
    use krome_constants
    implicit none
    real*8::fice,expEice1,expEice2,expEbare1,expEbare2,invphi
    real*8::nu0,p,dust_2body_rate,fbare,pesc_ice,pesc_bare

    !no need to calculate this if the dust is not present
    dust_2body_rate = 0d0

    fbare = 1d0-fice
    nu0 = 1d12 ! 1/s
    dust_2body_rate = fbare * (expEbare1 + expEbare2) * pesc_bare &
        + fice * (expEice1 + expEice2) * pesc_ice
    dust_2body_rate = dust_2body_rate * p * nu0 * invphi

  end function dust_2body_rate

  !******************
  function krate_2bodySi(n,idx1,idx2,Ea,Tdust) result(krate)
    use krome_commons
    implicit none
    real*8,intent(in)::n(nspec),Ea,Tdust
    integer,intent(in)::idx1,idx2
    real*8::krate,amin,amax,pexp,d2g,rho0

    !some default values OK for silicates
    amin = 5d-7 !cm
    amax = 2.5d-5 !cm
    pexp = -3.5
    rho0 = 3d0 !g/cm3
    d2g = 1d-2

    krate = krate_2body(n(:),idx1,idx2,amin,amax,pexp,d2g,rho0,Ea,Tdust)

  end function krate_2bodySi

  !********************
  function krate_2body(n,idx1,idx2,amin,amax,pexp,d2g,rho0, &
        Ea,Tdust) result(krate)
    use krome_commons
    use krome_constants
    use krome_getphys
    implicit none
    integer,intent(in)::idx1,idx2
    real*8,intent(in)::n(nspec),amin,amax,pexp,d2g,rho0,Ea,Tdust
    real*8::rhog,p3,p4,ndns,krate,mred,fice,fbare,Preac
    real*8::iTd23,Ebare(nspec),Eice(nspec),mass(nspec)
    real*8,parameter::app2=(3d-8)**2 !cm^2 (Hocuk+2015)
    real*8,parameter::nu0=1d12 !1/s
    real*8,parameter::hbar=planck_erg/2d0/pi !erg*s
    real*8,parameter::ar=1d-8 !cm

    mass(:) = get_mass()

    !gas density, g/cm3
    rhog = sum(mass(1:nmols)*n(1:nmols))

    !exponentes
    p3 = pexp + 3d0
    p4 = pexp + 4d0

    !number of sites cm-3/mly
    ndns = rhog/(4d0/3d0*rho0*app2)*(amax**p3-amin**p3) &
        / (amax**p4-amin**p4) * p4 / p3

    !ice/bare fraction
    fbare = 1d0

    !reduced mass
    mred = mass(idx1)*mass(idx2)/(mass(idx1)+mass(idx2))

    !tunneling probability
    Preac = exp(-2d0*ar/hbar*sqrt(2d0*mred*Ea*boltzmann_erg))

    !exponent
    iTd23 = 2d0/3d0/Tdust

    !get Ebind, K
    Ebare(:) = get_Ebind_bare()

    !compute rate
    krate = fbare*(exp(-Ebare(idx1)*iTd23)+exp(-Ebare(idx2)*iTd23))

    !rate in cm3/s
    krate = nu0*Preac/ndns*krate

  end function krate_2body

  !*************************
  function dust_get_inv_phi(asize2,nndust)
    use krome_commons
    use krome_constants
    implicit none
    real*8::iapp2,dust_get_inv_phi(ndust),asize2(ndust)
    real*8::nndust(ndust),dephi
    integer::i

    iapp2 = (3d-8)**2 !1/cm2
    do i=1,ndust
      dust_get_inv_phi(i) = 0d0
      dephi = (4d0 * nndust(i) * pi * asize2(i))
      if(dephi.le.0d0) cycle
      dust_get_inv_phi(i) = iapp2 / dephi
    end do

  end function dust_get_inv_phi

  !****************************
  !returns an array with the sticking coefficient for each bin
  ! following Hollenbach+McKee 1979
  function dust_stick_array(Tgas,Tdust)
    use krome_commons
    implicit none
    real*8::dust_stick_array(ndust),Tgas,Tdust(ndust)
    real*8::Tg100,Td100
    integer::i

    Tg100 = Tgas * 1d-2
    do i=1,ndust
      Td100 = Tdust(i) * 1d-2
      dust_stick_array(i) = 1d0/(1d0+.4d0*sqrt(Tg100+Td100) &
          + .2d0*Tg100 + 0.08d0*Tg100**2)
    end do

  end function dust_stick_array

  !*************************
  function dust_stick(Tgas,Tdust)
    implicit none
    real*8,intent(in)::Tgas,Tdust
    real*8::dust_stick
    real*8::Tg100,Td100

    Tg100 = Tgas * 1d-2
    Td100 = Tdust * 1d-2
    dust_stick = 1d0/(1d0 + 0.4d0*sqrt(Tg100+Td100) &
        + 0.2d0*Tg100 + 0.08d0*Tg100**2)

  end function dust_stick

  !****************************
  !sticking rate (1/s), assuming power-law dust distribution
  ! example rate is
  !  @format:idx,R,P,rate
  !  1,CO,CO_ice,krate_stick(n(:),idx_CO,1d-7,1d-5,-3.5,3d0,1d-2)
  ! n(:): internal status array (number densities, temeperature, etc...)
  ! idx : index of the sticking species, e.g. idx_CO
  ! Tdust: dust temperature (assume same for all bins), K
  ! amin: min grain size, cm
  ! amax: max grain size, cm
  ! pexp: power-law exponent, usually -3.5
  ! rho0: bulk material density, g/cm3, e.g. 3 g/cm3 for silicates
  ! d2g: dust to gass mass ratio, usually 0.01
  function krate_stick(n,idx,Tdust,amin,amax,pexp,rho0,d2g) result(k)
    use krome_constants
    use krome_commons
    use krome_getphys
    implicit none
    real*8,intent(in)::n(nspec),Tdust,amin,amax,pexp,rho0,d2g
    real*8::k,imass(nspec),p4,p3,mass(nspec),rhod
    integer,intent(in)::idx

    !get inverse mass squared
    imass(:) = get_imass_sqrt()
    !get masses
    mass(:) = get_mass()
    !derived exponents
    p3 = pexp + 3.
    p4 = pexp + 4.

    !total dust density, g/cm3
    rhod = sum(n(1:nmols)*mass(1:nmols))*d2g

    !compute rate (1/s) coefficient assuming normalization
    k = pre_kvgas_sqrt*sqrt(n(idx_Tgas)) * imass(idx) &
        * rhod / (4./3.*rho0) * p4 / p3 &
        * (amax**p3-amin**p3) / (amax**p4-amin**p4) &
        * dust_stick(n(idx_Tgas),Tdust)

  end function krate_stick

  !********************************
  !compact version of krate_stick
  function krate_stickSi(n,idx,Tdust) result(k)
    use krome_commons
    implicit none
    integer,intent(in)::idx
    real*8,intent(in)::n(nspec),Tdust
    real*8::k,amin,amax,d2g,rho0,pexp

    !some default values OK for silicates
    amin = 5d-7 !cm
    amax = 2.5d-5 !cm
    pexp = -3.5
    rho0 = 3d0 !g/cm3
    d2g = 1d-2

    k = krate_stick(n(:),idx,Tdust,amin,amax,pexp,rho0,d2g)

  end function krate_stickSi

  !***************************
  !evaporation rate, 1/s
  function krate_evaporation(n,idx,Tdust) result(k)
    use krome_commons
    use krome_getphys
    implicit none
    integer,intent(in)::idx
    real*8,intent(in)::n(nspec),Tdust
    real*8::k,Ebind(nspec),nu0

    nu0 = 1d12 !1/s
    Ebind(:) = get_EbindBare()

    k = nu0 * exp(-Ebind(idx)/Tdust)

  end function krate_evaporation

  !***************************
  !non-thermal evaporation rate (1/s) following Hollenbach 2009,
  ! http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:0809.1642
  !Gnot is the habing flux (1.78 is Draine)
  !Av is the visual extinction
  !crflux the ionization flux of cosmic rays, 1/s
  !yield is the efficiency of the photons to desorb the given molecule
  function krate_nonthermal_evaporation(idx, Gnot, Av, crflux, yield) result(k)
    use krome_commons
    use krome_getphys
    implicit none
    integer,intent(in)::idx
    real*8,parameter::crnot=1.3d-17
    real*8,parameter::Fnot=1d8 !desorbing photons flux, 1/s
    real*8,parameter::ap2=(3d-8)**2 !sites separation squared, cm2
    real*8,intent(in)::Gnot, Av, crflux, yield
    real*8::k,f70,kevap70(nspec)

    f70 = 3.16d-19*crflux/crnot
    kevap70(:) = get_kevap70()

    k = Gnot*Fnot*ap2*yield*exp(-1.8*Av)
    k = k + f70*kevap70(idx)

  end function krate_nonthermal_evaporation

  !***************************
  function dust_ice_fraction_array(invphi,nH2O)
    use krome_constants
    use krome_commons
    implicit none
    integer::i
    real*8::dust_ice_fraction_array(ndust)
    real*8::invphi(ndust),nH2O(ndust)

    do i=1,ndust
      dust_ice_fraction_array(i) = min(nH2O(i) * invphi(i), 1d0)
    end do

  end function dust_ice_fraction_array

  !*****************************
  function get_Ebareice_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice_exp_array(:) = 0d0

  end function get_Ebareice_exp_array

  !*****************************
  function get_Ebareice23_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice23_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice23_exp_array(:) = 0d0

  end function get_Ebareice23_exp_array

  !************************
  !returns the binding energy for ice coated grain (K)
  function get_Ebind_ice()
    use krome_commons
    implicit none
    real*8::get_Ebind_ice(nspec)

    get_Ebind_ice(:) = 0d0

  end function get_Ebind_ice

  !************************
  !returns the binding energy for bare grain (K)
  function get_Ebind_bare()
    use krome_commons
    implicit none
    real*8::get_Ebind_bare(nspec)

    get_Ebind_bare(:) = 0d0

  end function get_Ebind_bare

  !************************
  !returns the index of the parent dust bin (0 if none)
  function get_parent_dust_bin()
    use krome_commons
    implicit none
    integer::get_parent_dust_bin(nspec)

    get_parent_dust_bin(:) = 0

  end function get_parent_dust_bin

  !*****************************
  function get_exp_table(ain,invT)
    use krome_commons
    implicit none
    integer::ia
    real*8::get_exp_table,a,invT,ain
    real*8::x1a,f1,f2

    a = ain*invT
    a = min(a, exp_table_aMax - exp_table_da)

    ia = (a-exp_table_aMin) * exp_table_multa + 1
    ia = max(ia,1)

    x1a = (ia-1)*exp_table_da

    f1 = exp_table(ia)
    f2 = exp_table(ia+1)

    get_exp_table = (a-x1a) * exp_table_multa * (f2-f1) + f1

  end function get_exp_table

end module krome_grfuncs
!This module mainly contains shielding routine and
! function to initialize radiation background (e.g. Planck).

!############### MODULE ##############
module krome_phfuncs
contains

  !****************************
  !dust shielding factor
  function shield_dust(n,Tgas,gam)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::shield_dust,n(:),Tgas,gam,eff_d2g
    real*8::sigma_d,NHtot

    eff_d2g = dust2gas_ratio
    sigma_d = 2d-21*eff_d2g*gam !Richings et al. 2014
    !sigma_d = 2d-21 !Glover+2007
    !sigma_d = 4d-22 !Richings+ 2014
    !sigma_d = 4d-21 !Gnedin 2009

    NHtot = 0d0
    NHtot  = NHtot + num2col(n(idx_H),n(:))
    NHtot  = NHtot + num2col(n(idx_Hj),n(:))
    NHtot  = NHtot + 2d0 * num2col(n(idx_H2),n(:))

    shield_dust = exp(-sigma_d*NHtot)

  end function shield_dust

  !**********************
  !planck function in eV/s/cm2/Hz/sr
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB(x,Tbb)
    use krome_constants
    implicit none
    real*8::Tbb,x,xexp,planckBB

    !exponent
    xexp = x/boltzmann_eV/Tbb

    !default value
    planckBB = 0d0

    !limit exp overflow
    if(xexp<3d2.and.x>1d-10) then
      planckBB = 2d0*x**3/planck_eV**2/clight**2 &
          / (exp(xexp)-1d0)
    end if

  end function planckBB

  !********************
  !planck function dTdust differential
  ! in eV/s/cm2/Hz/sr/K, where
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB_dT(x,Tbb)
    use krome_constants
    real*8::a,b,x,Tbb,xexp,planckBB_dT

    b = 1d0/boltzmann_eV
    xexp = b*x/Tbb

    planckBB_dT = 0d0

    if(xexp<3d2) then
      a = 2d0/planck_eV**2/clight**2
      planckBB_dT = a*b*x**4/Tbb/Tbb * exp(xexp)/(exp(xexp)-1d0)**2
    end if

  end function planckBB_dT

  !***********************
  !shielding function selected with -shield option
  function krome_fshield(n,Tgas)
    implicit none
    real*8::krome_fshield,n(:),Tgas

    krome_fshield = 1d0 !default shielding value

  end function krome_fshield

  !**************************
  !shielding function for H2O+ and H3O+
  ! following Glover+2010 MNRAS sect 2.2 eqn.4
  function fHnOj(Av)
    implicit none
    real*8::fHnOj,Av
    if(Av.le.15d0) then
      fHnOj = exp(-2.55*Av+0.0165*Av**2)
    else
      fHnOj = exp(-2.8*Av)
    end if
  end function fHnOj

  !******************************
  !self-shielding for H2
  ! following Glover+2010 MNRAS sect2.2 eqn.6
  ! N: column density (cm-2)
  ! b: doppler broadening (cm/s)
  function fselfH2(N, b)
    implicit none
    real*8::fselfH2,N,b,x,b5

    x = N*2d-15 !normalized column density (#)
    b5 = b*1d-5 !normalized doppler broadening (#)

    fselfH2 = 0.965d0/(1+x/b5)**2 + &
        0.035d0/sqrt(1d0+x) * &
        exp(max(-8.5d-4*sqrt(1+x),-250.))

  end function fselfH2

end module krome_phfuncs

!############### MODULE ##############
module krome_subs
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2018-01-25 18:47:28
  !  Changeset 9a6e335
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !************************
  !compute reaction rates cm^3(n-1)/s
  function coe(n)
    use krome_commons
    use krome_constants
    use krome_user_commons
    use krome_getphys
    use krome_grfuncs
    use krome_phfuncs
    use krome_fit
    implicit none
    real*8::coe(nrea),k(nrea),Tgas,n(nspec),kmax
    real*8::Te
    real*8::lnTe
    real*8::T32
    real*8::invT
    real*8::small,nmax
    integer::i
    real*8::asav  !preproc from coevar
    real*8::bsav1 !preproc from coevar
    real*8::bsav3 !preproc from coevar
    real*8::bsav2 !preproc from coevar
    real*8::bsav5 !preproc from coevar
    real*8::bsav4 !preproc from coevar
    real*8::bsav7 !preproc from coevar
    real*8::bsav6 !preproc from coevar
    real*8::bsav8 !preproc from coevar
    real*8::sumsav !preproc from coevar
    !Tgas is in K
    Tgas = max(n(idx_Tgas), phys_Tcmb)
    Tgas = min(Tgas,1d9)

    !maxn initialization can be removed and small can be
    ! replaced with a proper value according to the environment
    nmax = max(maxval(n(1:nmols)),1d0)
    small = 1d-40/(nmax*nmax*nmax)

    Te = Tgas*8.617343d-5 !Tgas in eV (eV)
    lnTe = log(Te) !ln of Te (#)
    T32 = Tgas*0.0033333333333333335 !Tgas/(300 K) (#)
    invT = 1.d0/Tgas !inverse of T (1/K)

    asav = 2.1237150d4
    bsav1 = -3.3232183d-7
    bsav2 = 3.3735382d-7
    bsav3 = -1.4491368d-7
    bsav4 = 3.4172805d-8
    bsav5 = -4.7813728d-9
    bsav6 = 3.9731542d-10
    bsav7 = -1.8171411d-11
    bsav8 = 3.5311932d-13
    sumsav = bsav1+bsav2*log(Tgas)+bsav3*(log(Tgas))**2+bsav4*(log(Tgas))**3+bsav5*(log(Tgas))**4+bsav6*(log(Tgas))**5+bsav7*(log(Tgas))**6+bsav8*(log(Tgas))**7

    k(:) = small !inizialize coefficients

    !H + HE+ -> HE + H+
    k(1) = small + (1.20e-15*(T32)**(0.25))

    !C + NO -> CO + N
    k(2) = small + (9.00e-11*(T32)**(-0.16))

    !NH+ + E -> N + H
    k(3) = small + (4.30e-08*(T32)**(-0.50))

    !HE+ + E -> HE
    k(4) = small + (5.36e-12*(T32)**(-0.50))

    !O2 + S -> SO + O
    k(5) = small + (1.76e-12*(T32)**(0.81)&
        *exp(+30.8*invT))

    !H + HS+ -> S+ + H2
    k(6) = small + (1.10e-10)

    !H+ + SI -> SI+ + H
    k(7) = small + (9.90e-10)

    !H+ + OH -> OH+ + H
    k(8) = small + (2.10e-09&
        *(T32)**(-0.50))

    !N + HS -> NS + H
    k(9) = small + (1.00e-10)

    !HE+ + SI -> SI+ + HE
    k(10) = small + (3.30e-09)

    !H + HEH+ -> HE + H2+
    k(11) = small + (9.10e-10)

    !H2 + HS -> H2S + H
    k(12) = small + (6.52e-12*(T32)**(0.09)&
        *exp(-8050.0*invT))

    !O + E -> O-
    k(13) = small + (1.50e-15)

    !S+ + FE -> FE+ + S
    k(14) = small + (1.80e-10)

    !H + SI+ -> SIH+
    k(15) = small + (1.17e-17*(T32)**(-0.14))

    !H + NS -> HS + N
    k(16) = small + (7.27e-11*(T32)**(0.50)&
        *exp(-15700.0*invT))

    !N + CS -> S + CN
    k(17) = small + (3.80e-11*(T32)**(0.50)&
        *exp(-1160.0*invT))

    !H + S2 -> HS + S
    k(18) = small + (2.25e-10*(T32)**(0.50)&
        *exp(-8355.0*invT))

    !NA + S+ -> S + NA+
    k(19) = small + (2.60e-10)

    !OH + F -> HF + O
    k(20) = small + (1.60e-10)

    !H + SO -> S + OH
    k(21) = small + (5.90e-10*(T32)**(-0.31)&
        *exp(-11100.0*invT))

    !H- + NA+ -> H + NA
    k(22) = small + (7.51e-08&
        *(T32)**(-0.50))

    !CH + S -> CS + H
    k(23) = small + (5.00e-11)

    !H- + FE+ -> H + FE
    k(24) = small + (7.51e-08&
        *(T32)**(-0.50))

    !HEH+ + E -> HE + H
    k(25) = small + (1.00e-08&
        *(T32)**(-0.60))

    !H2 + O2 -> OH + OH
    k(26) = small + (3.16e-10*exp(-21890.0&
        *invT))

    !OH + S -> SO + H
    k(27) = small + (6.60e-11)

    !H + O2 -> O + O + H
    k(28) = small + (6.00e-09*exp(-52300.0&
        *invT))

    !H+ + E -> H
    k(29) = small + (3.50e-12*(T32)**(-0.75))

    !H + NO -> OH + N
    k(30) = small + (3.60e-10*exp(-24910.0&
        *invT))

    !O + SO2 -> SO + O2
    k(31) = small + (9.01e-12*exp(-9837.0&
        *invT))

    !H2+ + HE -> HEH+ + H
    k(32) = small + (1.30e-10)

    !H + H2 -> H + H + H
    k(33) = small + (4.67e-07&
        *(T32)**(-1.00)*exp(-55000.0*invT))

    !C + CS -> S + C2
    k(34) = small + (1.44e-11*(T32)**(0.50)&
        *exp(-20435.0*invT))

    !N + NO -> N2 + O
    k(35) = small + (3.38e-11*(T32)**(-0.17)&
        *exp(+2.8*invT))

    !H + CO -> OH + C
    k(36) = small + (1.10e-10*(T32)**(0.50)&
        *exp(-77700.0*invT))

    !H2 + F -> HF + H
    k(37) = small + (1.00e-10*exp(-400.0&
        *invT))

    !H + CH2 -> CH + H2
    k(38) = small + (2.20e-10)

    !O + O -> O2
    k(39) = small + (4.90e-20*(T32)**(1.58))

    !H + NH -> N + H2
    k(40) = small + (1.73e-11*(T32)**(0.50)&
        *exp(-2400.0*invT))

    !NH + O -> OH + N
    k(41) = small + (1.16e-11)

    !H2 + CH -> CH2 + H
    k(42) = small + (5.46e-10*exp(-1943.0&
        *invT))

    !H + HCN -> CN + H2
    k(43) = small + (6.20e-10*exp(-12500.0&
        *invT))

    !OH + CO -> CO2 + H
    k(44) = small + (2.81e-13*exp(-176.0&
        *invT))

    !NH + O -> NO + H
    k(45) = small + (6.60e-11)

    !C + CO -> C2 + O
    k(46) = small + (2.94e-11*(T32)**(0.50)&
        *exp(-58025.0*invT))

    !C + CN -> C2 + N
    k(47) = small + (4.98e-10*exp(-18116.0&
        *invT))

    !N + NH -> N2 + H
    k(48) = small + (4.98e-11)

    !SI + O2 -> SIO + O
    k(49) = small + (1.72e-10&
        *(T32)**(-0.53)*exp(-17.0*invT))

    !H + O- -> OH + E
    k(50) = small + (5.00e-10)

    !OH + SIO -> SIO2 + H
    k(51) = small + (2.00e-12)

    !C + HCO+ -> CO + CH+
    k(52) = small + (1.10e-09)

    !H + NH2 -> NH + H2
    k(53) = small + (4.56e-12*(T32)**(1.02)&
        *exp(-2161.0*invT))

    !H2 + C -> CH + H
    k(54) = small + (6.64e-10*exp(-11700.0&
        *invT))

    !H2 + O+ -> OH+ + H
    k(55) = small + (1.70e-09)

    !H + OCN -> OH + CN
    k(56) = small + (1.00e-10)

    !H + E -> H-
    k(57) = small + (3.37e-16*(T32)**(0.64)&
        *exp(-9.2*invT))

    !H2+ + O -> OH+ + H
    k(58) = small + (1.50e-09)

    !MG + S+ -> S + MG+
    k(59) = small + (2.80e-10)

    !N + O2 -> NO + O
    k(60) = small + (2.26e-12*(T32)**(0.86)&
        *exp(-3134.0*invT))

    !H+ + H -> H2+
    k(61) = small + (1.15e-18*(T32)**(1.49)&
        *exp(-228.0*invT))

    !H + HS -> S + H2
    k(62) = small + (2.50e-11)

    !O- + MG+ -> O + MG
    k(63) = small + (7.51e-08&
        *(T32)**(-0.50))

    !MG+ + E -> MG
    k(64) = small + (2.78e-12*(T32)**(-0.68))

    !H- + H+ -> H + H
    k(65) = small + (7.51e-08*(T32)**(-0.50))

    !H- + S+ -> H + S
    k(66) = small + (7.51e-08*(T32)**(-0.50))

    !H + OH -> O + H + H
    k(67) = small + (6.00e-09*exp(-50900.0&
        *invT))

    !O+ + FE -> FE+ + O
    k(68) = small + (2.90e-09)

    !H- + N -> NH + E
    k(69) = small + (1.00e-09)

    !H + NO -> O + NH
    k(70) = small + (9.29e-10*(T32)**(-0.10)&
        *exp(-35220.0*invT))

    !SI + HCO+ -> SIH+ + CO
    k(71) = small + (1.60e-09)

    !O + CN -> NO + C
    k(72) = small + (5.37e-11*exp(-13800.0&
        *invT))

    !H- + O+ -> H + O
    k(73) = small + (7.51e-08*(T32)**(-0.50))

    !O+ + E -> O
    k(74) = small + (3.24e-12*(T32)**(-0.66))

    !H2 + C -> CH2
    k(75) = small + (1.00e-17)

    !H + S- -> HS + E
    k(76) = small + (1.00e-10)

    !OH + CS -> CO + HS
    k(77) = small + (3.00e-11)

    !NA + FE+ -> FE + NA+
    k(78) = small + (1.00e-11)

    !SIO+ + FE -> FE+ + SIO
    k(79) = small + (1.00e-09)

    !C + SO -> CS + O
    k(80) = small + (3.50e-11)

    !OH + SI+ -> SIO+ + H
    k(81) = small + (6.30e-10&
        *(T32)**(-0.50))

    !O + CS -> SO + C
    k(82) = small + (4.68e-11*(T32)**(0.50)&
        *exp(-28940.0*invT))

    !H + O -> OH
    k(83) = small + (9.90e-19*(T32)**(-0.38))

    !H+ + P -> P+ + H
    k(84) = small + (1.00e-09)

    !CN + S -> NS + C
    k(85) = small + (5.71e-11*(T32)**(0.50)&
        *exp(-32010.0*invT))

    !O + HS -> S + OH
    k(86) = small + (1.74e-11*(T32)**(0.67)&
        *exp(-956.0*invT))

    !OH+ + E -> O + H
    k(87) = small + (3.75e-08*(T32)**(-0.50))

    !H + C -> CH
    k(88) = small + (1.00e-17)

    !SI + CO -> SIO + C
    k(89) = small + (1.30e-09*exp(-34513.0&
        *invT))

    !C + HS -> CS + H
    k(90) = small + (1.00e-10)

    !H2 + S -> HS + H
    k(91) = small + (1.76e-13*(T32)**(2.88)&
        *exp(-6126.0*invT))

    !H + CH -> C + H2
    k(92) = small + (1.31e-10*exp(-80.0&
        *invT))

    !H + O2 -> OH + O
    k(93) = small + (2.61e-10*exp(-8156.0&
        *invT))

    !N + SO -> NS + O
    k(94) = small + (4.68e-11*(T32)**(0.50)&
        *exp(-8254.0*invT))

    !H + HCO -> CO + H2
    k(95) = small + (1.50e-10)

    !OH + CN -> HCN + O
    k(96) = small + (1.00e-11*exp(-1000.0&
        *invT))

    !O + SIO+ -> O2 + SI+
    k(97) = small + (2.00e-10)

    !CH + O -> OH + C
    k(98) = small + (2.52e-11*exp(-2381.0&
        *invT))

    !CH + S -> HS + C
    k(99) = small + (1.73e-11*(T32)**(0.50)&
        *exp(-4000.0*invT))

    !HF + SI+ -> SIF+ + H
    k(100) = small + (5.70e-09&
        *(T32)**(-0.50))

    !C + NO -> CN + O
    k(101) = small + (6.00e-11&
        *(T32)**(-0.16))

    !H2+ + E -> H + H
    k(102) = small + (1.60e-08&
        *(T32)**(-0.43))

    !H + H2O -> OH + H + H
    k(103) = small + (5.80e-09&
        *exp(-52900.0*invT))

    !H2 + NH -> NH2 + H
    k(104) = small + (5.96e-11*exp(-7782.0&
        *invT))

    !C + NH -> N + CH
    k(105) = small + (1.73e-11*(T32)**(0.50)&
        *exp(-4000.0*invT))

    !C + N2 -> CN + N
    k(106) = small + (8.69e-11*exp(-22600.0&
        *invT))

    !H2 + E -> H + H + E
    k(107) = small + (3.22e-09&
        *(T32)**(0.35)*exp(-102000.0*invT))

    !H + CH -> C + H + H
    k(108) = small + (6.00e-09&
        *exp(-40200.0*invT))

    !SIF+ + E -> SI + F
    k(109) = small + (2.00e-07&
        *(T32)**(-0.50))

    !H + SIH+ -> SI+ + H2
    k(110) = small + (1.90e-09)

    !C+ + MG -> MG+ + C
    k(111) = small + (1.10e-09)

    !C + O -> CO
    k(112) = small + (4.69e-19*(T32)**(1.52)&
        *exp(+50.5*invT))

    !CN + O2 -> OCN + O
    k(113) = small + (2.02e-11&
        *(T32)**(-0.19)*exp(+31.9*invT))

    !NH + S -> NS + H
    k(114) = small + (1.00e-10)

    !SIO+ + E -> SI + O
    k(115) = small + (2.00e-07&
        *(T32)**(-0.50))

    !O + HCN -> CO + NH
    k(116) = small + (7.30e-13&
        *(T32)**(1.14)*exp(-3742.0*invT))

    !H2 + N -> NH + H
    k(117) = small + (1.69e-09*exp(-18095.0&
        *invT))

    !O + HCN -> CN + OH
    k(118) = small + (6.21e-10*exp(-12439.0&
        *invT))

    !SI + S+ -> S + SI+
    k(119) = small + (1.60e-09)

    !O + OH -> O2 + H
    k(120) = small + (3.69e-11*(T32)**(-0.27)&
        *exp(-12.9*invT))

    !C+ + SI -> SI+ + C
    k(121) = small + (2.10e-09)

    !O + SO -> S + O2
    k(122) = small + (6.60e-13*exp(-2760.0&
        *invT))

    !SI + NO -> SIO + N
    k(123) = small + (9.00e-11&
        *(T32)**(-0.96)*exp(-28.0*invT))

    !C+ + FE -> FE+ + C
    k(124) = small + (2.60e-09)

    !S + SO2 -> SO + SO
    k(125) = small + (9.76e-12*exp(-4545.0&
        *invT))

    !MG + SIO+ -> SIO + MG+
    k(126) = small + (1.00e-09)

    !O + HCN -> OCN + H
    k(127) = small + (1.36e-12&
        *(T32)**(1.38)*exp(-3693.0*invT))

    !N + CN -> N2 + C
    k(128) = small + (1.00e-10*(T32)**(0.18))

    !H- + C -> CH + E
    k(129) = small + (1.00e-09)

    !MG + SI+ -> SI + MG+
    k(130) = small + (2.90e-09)

    !CH + N -> NH + C
    k(131) = small + (3.03e-11*(T32)**(0.65)&
        *exp(-1207.0*invT))

    !NA + SI+ -> SI + NA+
    k(132) = small + (2.70e-09)

    !C- + H+ -> C + H
    k(133) = small + (7.51e-08&
        *(T32)**(-0.50))

    !SI + CO2 -> SIO + CO
    k(134) = small + (2.72e-11*exp(-282.0&
        *invT))

    !H + OCS -> HS + CO
    k(135) = small + (1.23e-11*exp(-1949.0&
        *invT))

    !SI+ + FE -> FE+ + SI
    k(136) = small + (1.90e-09)

    !S + E -> S-
    k(137) = small + (5.00e-15)

    !HCO+ + FE -> FE+ + HCO
    k(138) = small + (1.90e-09)

    !P+ + E -> P
    k(139) = small + (3.41e-12*(T32)**(-0.65))

    !H+ + SIO -> SIO+ + H
    k(140) = small + (3.30e-09&
        *(T32)**(-0.50))

    !O + HS -> SO + H
    k(141) = small + (1.74e-10*(T32)**(-0.20)&
        *exp(-5.7*invT))

    !FE+ + E -> FE
    k(142) = small + (2.55e-12*(T32)**(-0.69))

    !SI + P+ -> P + SI+
    k(143) = small + (1.00e-09)

    !H+ + MG -> MG+ + H
    k(144) = small + (1.10e-09)

    !H- + SI+ -> H + SI
    k(145) = small + (7.51e-08&
        *(T32)**(-0.50))

    !H + H2O -> OH + H2
    k(146) = small + (1.59e-11&
        *(T32)**(1.20)*exp(-9610.0*invT))

    !OH + OH -> H2O + O
    k(147) = small + (1.65e-12&
        *(T32)**(1.14)*exp(-50.0*invT))

    !OH + SO -> SO2 + H
    k(148) = small + (8.60e-11)

    !H2 + OH -> H2O + H
    k(149) = small + (2.05e-12&
        *(T32)**(1.52)*exp(-1736.0*invT))

    !H + OCN -> HCN + O
    k(150) = small + (1.87e-11&
        *(T32)**(0.90)*exp(-2924.0*invT))

    !N + C2 -> CN + C
    k(151) = small + (5.00e-11)

    !C + CH -> C2 + H
    k(152) = small + (6.59e-11)

    !H + C2 -> CH + C
    k(153) = small + (4.67e-10*(T32)**(0.50)&
        *exp(-30450.0*invT))

    !H + OCN -> NH + CO
    k(154) = small + (1.26e-10*exp(-515.0&
        *invT))

    !N + OH -> O + NH
    k(155) = small + (1.88e-11*(T32)**(0.10)&
        *exp(-10700.0*invT))

    !N + SIO+ -> NO + SI+
    k(156) = small + (2.10e-10)

    !C + E -> C-
    k(157) = small + (2.25e-15)

    !NH + S -> HS + N
    k(158) = small + (1.73e-11*(T32)**(0.50)&
        *exp(-4000.0*invT))

    !H2 + O -> OH + H
    k(159) = small + (3.14e-13*(T32)**(2.70)&
        *exp(-3150.0*invT))

    !O + CS -> S + CO
    k(160) = small + (2.48e-10*(T32)**(-0.65)&
        *exp(-783.0*invT))

    !OH + H2S -> HS + H2O
    k(161) = small + (6.30e-12*exp(-80.0&
        *invT))

    !C + NH -> CN + H
    k(162) = small + (1.20e-10)

    !H + H2S -> HS + H2
    k(163) = small + (3.71e-12&
        *(T32)**(1.94)*exp(-455.0*invT))

    !MG + HCO+ -> HCO + MG+
    k(164) = small + (2.90e-09)

    !CH + N -> CN + H
    k(165) = small + (1.66e-10&
        *(T32)**(-0.09))

    !C + C -> C2
    k(166) = small + (4.36e-18*(T32)**(0.35)&
        *exp(-161.3*invT))

    !H+ + S -> S+ + H
    k(167) = small + (1.30e-09)

    !O + C2 -> CO + C
    k(168) = small + (2.00e-10&
        *(T32)**(-0.12))

    !C + OH -> O + CH
    k(169) = small + (2.25e-11*(T32)**(0.50)&
        *exp(-14800.0*invT))

    !H + C- -> CH + E
    k(170) = small + (5.00e-10)

    !C + HS -> S + CH
    k(171) = small + (1.20e-11*(T32)**(0.58)&
        *exp(-5880.0*invT))

    !H+ + FE -> FE+ + H
    k(172) = small + (7.40e-09)

    !C + N -> CN
    k(173) = small + (5.72e-19*(T32)**(0.37)&
        *exp(-51.0*invT))

    !OH + CS -> H + OCS
    k(174) = small + (1.70e-10)

    !N + CO2 -> NO + CO
    k(175) = small + (3.20e-13*exp(-1710.0&
        *invT))

    !H + CO2 -> CO + OH
    k(176) = small + (3.38e-10*exp(-13163.0&
        *invT))

    !C2 + S -> CS + C
    k(177) = small + (1.00e-10)

    !H2 + S+ -> HS+ + H
    k(178) = small + (1.10e-10*exp(-9860.0&
        *invT))

    !N + HS -> S + NH
    k(179) = small + (1.73e-11*(T32)**(0.50)&
        *exp(-9060.0*invT))

    !O + N2 -> NO + N
    k(180) = small + (2.51e-10*exp(-38602.0&
        *invT))

    !CH + O -> CO + H
    k(181) = small + (6.02e-11*(T32)**(0.10)&
        *exp(+4.5*invT))

    !H + O+ -> O + H+
    k(182) = small + (5.66e-10*(T32)**(0.36)&
        *exp(+8.6*invT))

    !N + OH -> NO + H
    k(183) = small + (6.05e-11*(T32)**(-0.23)&
        *exp(-14.9*invT))

    !H+ + O -> O+ + H
    k(184) = small + (6.86e-10*(T32)**(0.26)&
        *exp(-224.3*invT))

    !H + SO -> HS + O
    k(185) = small + (1.73e-11*(T32)**(0.50)&
        *exp(-19930.0*invT))

    !O + SI -> SIO
    k(186) = small + (5.52e-18*(T32)**(0.31))

    !S+ + E -> S
    k(187) = small + (5.49e-12*(T32)**(-0.59))

    !HS + HS -> H2S + S
    k(188) = small + (1.30e-11)

    !C + SO2 -> CO + SO
    k(189) = small + (7.00e-11)

    !NA+ + E -> NA
    k(190) = small + (2.76e-12*(T32)**(-0.68))

    !S + HS -> S2 + H
    k(191) = small + (4.50e-11)

    !H- + H -> H2 + E
    k(192) = small + (4.82e-09*(T32)**(0.02)&
        *exp(-4.3*invT))

    !C + SO -> S + CO
    k(193) = small + (3.50e-11)

    !H + CH+ -> C+ + H2
    k(194) = small + (9.06e-10&
        *(T32)**(-0.37)*exp(-29.1*invT))

    !NA + MG+ -> MG + NA+
    k(195) = small + (1.00e-11)

    !H + H2+ -> H2 + H+
    k(196) = small + (6.40e-10)

    !SI+ + E -> SI
    k(197) = small + (4.26e-12*(T32)**(-0.62))

    !C + NS -> S + CN
    k(198) = small + (1.50e-10&
        *(T32)**(-0.16))

    !O- + H+ -> O + H
    k(199) = small + (7.51e-08&
        *(T32)**(-0.50))

    !CH + O -> HCO+ + E
    k(200) = small + (1.09e-11&
        *(T32)**(-2.19)*exp(-165.1*invT))

    !SIH+ + E -> SI + H
    k(201) = small + (2.00e-07&
        *(T32)**(-0.50))

    !C + SIO+ -> SI+ + CO
    k(202) = small + (1.00e-09)

    !H- + O -> OH + E
    k(203) = small + (1.00e-09)

    !C + OH -> CO + H
    k(204) = small + (1.00e-10)

    !OH + SI -> SIO + H
    k(205) = small + (1.00e-10)

    !H2+ + C -> CH+ + H
    k(206) = small + (2.40e-09)

    !C + S -> CS
    k(207) = small + (4.36e-19*(T32)**(0.22))

    !H + OH -> O + H2
    k(208) = small + (6.99e-14*(T32)**(2.80)&
        *exp(-1950.0*invT))

    !CH+ + E -> C + H
    k(209) = small + (1.50e-07&
        *(T32)**(-0.42))

    !O- + FE+ -> O + FE
    k(210) = small + (7.51e-08&
        *(T32)**(-0.50))

    !N + SO -> S + NO
    k(211) = small + (1.73e-11*(T32)**(0.50)&
        *exp(-750.0*invT))

    !O + NS -> S + NO
    k(212) = small + (1.00e-10)

    !O + H2O -> OH + OH
    k(213) = small + (1.85e-11&
        *(T32)**(0.95)*exp(-8571.0*invT))

    !H- + MG+ -> H + MG
    k(214) = small + (7.51e-08&
        *(T32)**(-0.50))

    !H2 + CN -> HCN + H
    k(215) = small + (4.04e-13&
        *(T32)**(2.87)*exp(-820.0*invT))

    !C + O2 -> CO + O
    k(216) = small + (5.56e-11*(T32)**(0.41)&
        *exp(+26.9*invT))

    !H + NS -> S + NH
    k(217) = small + (7.27e-11*(T32)**(0.50)&
        *exp(-20735.0*invT))

    !O + CN -> CO + N
    k(218) = small + (2.54e-11)

    !OH + CN -> OCN + H
    k(219) = small + (7.00e-11)

    !H + E -> H+ + E + E
    k(220) = small + (exp(-32.71396786d0+13.5365560d0&
        *lnTe-5.73932875d0*(lnTe**2)+1.56315498d0&
        *(lnTe**3)-0.28770560d0*(lnTe**4)+3.48255977d-2&
        *(lnTe**5)-2.63197617d-3*(lnTe**6)+1.11954395d-4&
        *(lnTe**7)-2.03914985d-6*(lnTe**8)))

    !HE + E -> HE+ + E + E
    k(221) = small + (exp(-44.09864886d0+23.91596563d0&
        *lnTe-10.7532302d0*(lnTe**2)+3.05803875d0&
        *(lnTe**3)-0.56851189d0*(lnTe**4)+6.79539123d-2&
        *(lnTe**5)-5.00905610d-3*(lnTe**6)+2.06723616d-4&
        *(lnTe**7)-3.64916141d-6*(lnTe**8)))

    !H2 + H+ -> H2+ + H
    k(222) = small + (sumsav*exp(-asav&
        *invT))

    !H2 + E -> H + H-
    k(223) = small + (3.55d1*Tgas**(-2.28)&
        *exp(-46707.&
        /Tgas))

    !H- + E -> H + E + E
    k(224) = small + (exp(-18.01849334273d0+2.360852208681d0&
        *lnTe-0.2827443061704d0*lnTe**2+0.01623316639567d0&
        *lnTe**3-0.03365012031362999d0*lnTe**4+0.01178329782711d0&
        *lnTe**5-0.001656194699504d0*lnTe**6+0.0001068275202678d0&
        *lnTe**7-2.631285809207d-6*lnTe**8))

    !H- + H -> H + H + E
    k(225) = small + (exp(-20.37260896533324d0+1.139449335841631d0&
        *lnTe-0.1421013521554148d0*lnTe**2+0.00846445538663d0&
        *lnTe**3-0.0014327641212992d0*lnTe**4+0.0002012250284791d0&
        *lnTe**5+0.0000866396324309d0*lnTe**6-0.00002585009680264d0&
        *lnTe**7+2.4555011970392d-6*lnTe**8-8.06838246118d-8&
        *lnTe**9))

    !H- + H+ -> H2+ + E
    k(226) = small + (1.d-8*Tgas**(-0.4d0))

    !H + H + H -> H2 + H
    k(227) = small + (6.d-32&
        *Tgas**(-0.25d0)+2.d-31*Tgas**(-0.5d0))

    !H2 + H + H -> H2 + H2
    k(228) = small + ((6.d-32&
        *Tgas**(-0.25d0)+2.d-31*Tgas**(-0.5d0))&
        /8.d0)

    !H + H + HE -> H2 + HE
    k(229) = small + (6.9d-32&
        *Tgas**(-0.4d0))

    !H+ + NA -> NA+ + H
    k(230) = small + (1.20e-09)

    !H2 -> H+ + H-
    k(231) = small + (0.000286764705882&
        *user_crflux)

    !C -> C+ + E
    k(232) = small + (1.69117647059*user_crflux)

    !H -> H+ + E
    k(233) = small + (0.439705882353*user_crflux)

    !N -> N+ + E
    k(234) = small + (1.98529411765*user_crflux)

    !CO -> CO+ + E
    k(235) = small + (2.86764705882*user_crflux)

    !H2 -> H+ + H + E
    k(236) = small + (0.0210294117647&
        *user_crflux)

    !H2 -> H + H
    k(237) = small + (0.0955882352941*user_crflux)

    !HE -> HE+ + E
    k(238) = small + (0.477941176471&
        *user_crflux)

    !O -> O+ + E
    k(239) = small + (2.5*user_crflux)

    !H2 -> H2+ + E
    k(240) = small + (0.882352941176&
        *user_crflux)

    !N2 -> N + N
    k(241) = small + (5.000e+00*user_crflux)

    !CO -> C + O
    k(242) = small + (5.000e+00*user_crflux)

    !C+ + E -> C
    k(243) = small + (2.36e-12*(T32)**(-0.29)&
        *exp(+17.6*invT))

    !N+ + E -> N
    k(244) = small + (3.50e-12*(T32)**(-0.53)&
        *exp(+3.2*invT))

    !CO+ + E -> O + C
    k(245) = small + (2.00e-07&
        *(T32)**(-0.48))

    !HE+ + SIO2 -> O2 + SI+ + HE
    k(246) = small + (2.00e-09)

    !H+ + NH -> NH+ + H
    k(247) = small + (2.10e-09&
        *(T32)**(-0.50))

    !HS+ + E -> S + H
    k(248) = small + (2.00e-07&
        *(T32)**(-0.50))

    !HCO+ + E -> CO + H
    k(249) = small + (2.40e-07&
        *(T32)**(-0.69))

    !HE+ + HF -> F+ + H + HE
    k(250) = small + (1.10e-08&
        *(T32)**(-0.50))

    !H2 + F+ -> H2+ + F
    k(251) = small + (6.24e-10)

    !N + PN -> P + N2
    k(252) = small + (1.00e-18)

    !N + PO -> PN + O
    k(253) = small + (3.00e-11&
        *(T32)**(-0.60))

    !P + O2 -> PO + O
    k(254) = small + (1.00e-13)

    !N + PO -> P + NO
    k(255) = small + (2.55e-12)

    coe(:) = k(:) !set coefficients to return variable

    !!uncomment below to check coefficient values
    !kmax = 1d0
    !if(maxval(k)>kmax.or.minval(k)<0d0) then
    !   print *,"***************"
    !   do i=1,size(k)
    !      if(k(i)<0d0.or.k(i)>kmax) print *,i,k(i)
    !   end do
    !end if

  end function coe

  !*************************
  subroutine loadReactionsVerbatim()
    use krome_commons
    implicit none
    character*50::fname,line
    integer::ios,i,nunit

    fname = "reactions_verbatim.dat"

    !verbatim reactions are loaded from file
    ! to increase compilation speed
    open(newunit=nunit,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
      print *,"ERROR: "//trim(fname)//" file not present!"
      stop
    end if

    !load reactions from file
    do i=1,nrea
      read(nunit,'(a)',iostat=ios) line
      if(ios/=0) then
        print *,"ERROR: problem reading "//trim(fname)
        stop
      end if
      reactionNames(i) = trim(line)
    end do
    close(nunit)

  end subroutine loadReactionsVerbatim

  !*******************
  !The following functions compute the recombination rate
  ! on dust for H+, He+, C+, Si+, and O+. See Weingartner&Draine 2001
  ! dust2gas_ratio, D/D_sol, default is assumed equal to Z/Z_sol
  function H_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::H_recombination_on_dust

    H_recombination_on_dust = 0d0

    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    H_recombination_on_dust =  1.225d-13*dust2gas_ratio &
        /(1.d0+8.074d-6*psi**(1.378)*(1.d0+5.087d2 &
        *Tgas**(0.01586)*psi**(-0.4723-1.102d-5*log(Tgas))))

  end function H_recombination_on_dust

  !******************
  function He_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::He_recombination_on_dust

    He_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    He_recombination_on_dust = 5.572d-14*dust2gas_ratio&
        /(1.d0+3.185d-7*psi**(1.512)*(1.d0+5.115d3&
        *Tgas**(3.903d-7)*psi**(-0.4956-5.494d-7*log(Tgas))))

  end function He_recombination_on_dust

  !*******************
  function C_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::C_recombination_on_dust

    C_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    C_recombination_on_dust = 4.558d-13*dust2gas_ratio&
        /(1.d0+6.089d-3*psi**(1.128)*(1.d0+4.331d2&
        *Tgas**(0.04845)*psi**(-0.8120-1.333d-4*log(Tgas))))

  end function C_recombination_on_dust

  !******************
  function Si_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::Si_recombination_on_dust

    Si_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    Si_recombination_on_dust = 2.166d-14*dust2gas_ratio&
        /(1.d0+5.678d-8*psi**(1.874)*(1.d0+4.375d4&
        *Tgas**(1.635d-6)*psi**(-0.8964-7.538d-5*log(Tgas))))

  end function Si_recombination_on_dust

  !********************
  function O_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,k_H
    real*8::O_recombination_on_dust

    k_H = H_recombination_on_dust(n(:),Tgas)
    O_recombination_on_dust = 0.25d0*k_H

  end function O_recombination_on_dust

  !*********************
  !This function returns the
  ! photorate of H2 occurring in the
  ! Lyman-Werner bands following the approximation
  ! provided by Glover&Jappsen 2007. Rate in 1/s.
  !Approximation valid at low-density, it assumes H2(nu = 0).
  !It also stores the rate as a common, needed for the photoheating
  function H2_solomonLW(myflux)
    use krome_commons
    use krome_constants
    implicit none
    real*8::H2_solomonLW,myflux

    !myflux is the radiation background at E = 12.87 eV
    !should be converted to erg
    H2_solomonLW = 1.38d9*myflux*eV_to_erg

  end function H2_solomonLW

  !****************************
  !tanh smoothing function that
  ! increses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_increase(xarg,xpos,slope)
    implicit none
    real*8::smooth_increase,xarg,xpos,slope

    smooth_increase = .5d0 * (tanh(slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_increase

  !****************************
  !tanh smoothing function that
  ! decreses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_decrease(xarg,xpos,slope)
    implicit none
    real*8::smooth_decrease,xarg,xpos,slope

    smooth_decrease = .5d0 * (tanh(-slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_decrease

  !*********************
  !sign: return 1d0 if x>=0d0,
  ! else return -1d0
  function get_sgn(x)
    implicit none
    real*8::x,get_sgn

    get_sgn = 1d0
    if(x==0d0) return
    get_sgn = x/abs(x)

  end function get_sgn

  !*********************
  function conserve(n,ni)
    use krome_commons
    implicit none
    real*8::conserve(nspec),n(nspec),ni(nspec),no(nspec)
    real*8::ntot,nitot,factor

    no(:) = n(:)

    !********** E **********
    no(idx_E) = max( &
        -n(idx_Ok) &
        -n(idx_Hk) &
        -n(idx_Sk) &
        -n(idx_Ck) &
        +n(idx_HEj) &
        +n(idx_Hj) &
        +n(idx_NHj) &
        +n(idx_HSj) &
        +n(idx_Sj) &
        +n(idx_SIj) &
        +n(idx_OHj) &
        +n(idx_HEHj) &
        +n(idx_H2j) &
        +n(idx_FEj) &
        +n(idx_SIHj) &
        +n(idx_NAj) &
        +n(idx_HCOj) &
        +n(idx_CHj) &
        +n(idx_Oj) &
        +n(idx_MGj) &
        +n(idx_SIOj) &
        +n(idx_Pj) &
        +n(idx_SIFj) &
        +n(idx_Cj) &
        +n(idx_Nj) &
        +n(idx_COj) &
        +n(idx_Fj), 1d-40)

    conserve(:) = 0d0
    conserve(:) = no(:)

  end function conserve

  !*************************
  !this subroutine changes the x(:) mass fractions of the species
  ! to force conservation according to the reference ref(:)
  subroutine conserveLin_x(x,ref)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::x(nmols),ref(natoms)
    real*8::A(natoms,natoms),B(natoms),m(nspec)

    m(:) = get_mass()
    A(:,:) = 0d0
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_Ck) * m(idx_C) * m(idx_C) / m(idx_Ck)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_C) * m(idx_C) * m(idx_C) / m(idx_C)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_CO) * m(idx_C) * m(idx_C) / m(idx_CO)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_CS) * m(idx_C) * m(idx_C) / m(idx_CS)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_CN) * m(idx_C) * m(idx_C) / m(idx_CN)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_CH) * m(idx_C) * m(idx_C) / m(idx_CH)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) + 4d0 * x(idx_C2) * m(idx_C) * m(idx_C) / m(idx_C2)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_CH2) * m(idx_C) * m(idx_C) / m(idx_CH2)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_HCN) * m(idx_C) * m(idx_C) / m(idx_HCN)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_CO2) * m(idx_C) * m(idx_C) / m(idx_CO2)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_OCN) * m(idx_C) * m(idx_C) / m(idx_OCN)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_HCO) * m(idx_C) * m(idx_C) / m(idx_HCO)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_OCS) * m(idx_C) * m(idx_C) / m(idx_OCS)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_HCOj) * m(idx_C) * m(idx_C) / m(idx_HCOj)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_CHj) * m(idx_C) * m(idx_C) / m(idx_CHj)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_Cj) * m(idx_C) * m(idx_C) / m(idx_Cj)**2
    A(idx_atom_C, idx_atom_C) = A(idx_atom_C, idx_atom_C) +  x(idx_COj) * m(idx_C) * m(idx_C) / m(idx_COj)**2
    A(idx_atom_C, idx_atom_H) = A(idx_atom_C, idx_atom_H) +  x(idx_CH) * m(idx_C) * m(idx_H) / m(idx_CH)**2
    A(idx_atom_C, idx_atom_H) = A(idx_atom_C, idx_atom_H) + 2d0 * x(idx_CH2) * m(idx_C) * m(idx_H) / m(idx_CH2)**2
    A(idx_atom_C, idx_atom_H) = A(idx_atom_C, idx_atom_H) +  x(idx_HCN) * m(idx_C) * m(idx_H) / m(idx_HCN)**2
    A(idx_atom_C, idx_atom_H) = A(idx_atom_C, idx_atom_H) +  x(idx_HCO) * m(idx_C) * m(idx_H) / m(idx_HCO)**2
    A(idx_atom_C, idx_atom_H) = A(idx_atom_C, idx_atom_H) +  x(idx_HCOj) * m(idx_C) * m(idx_H) / m(idx_HCOj)**2
    A(idx_atom_C, idx_atom_H) = A(idx_atom_C, idx_atom_H) +  x(idx_CHj) * m(idx_C) * m(idx_H) / m(idx_CHj)**2
    A(idx_atom_C, idx_atom_O) = A(idx_atom_C, idx_atom_O) +  x(idx_CO) * m(idx_C) * m(idx_O) / m(idx_CO)**2
    A(idx_atom_C, idx_atom_O) = A(idx_atom_C, idx_atom_O) + 2d0 * x(idx_CO2) * m(idx_C) * m(idx_O) / m(idx_CO2)**2
    A(idx_atom_C, idx_atom_O) = A(idx_atom_C, idx_atom_O) +  x(idx_OCN) * m(idx_C) * m(idx_O) / m(idx_OCN)**2
    A(idx_atom_C, idx_atom_O) = A(idx_atom_C, idx_atom_O) +  x(idx_HCO) * m(idx_C) * m(idx_O) / m(idx_HCO)**2
    A(idx_atom_C, idx_atom_O) = A(idx_atom_C, idx_atom_O) +  x(idx_OCS) * m(idx_C) * m(idx_O) / m(idx_OCS)**2
    A(idx_atom_C, idx_atom_O) = A(idx_atom_C, idx_atom_O) +  x(idx_HCOj) * m(idx_C) * m(idx_O) / m(idx_HCOj)**2
    A(idx_atom_C, idx_atom_O) = A(idx_atom_C, idx_atom_O) +  x(idx_COj) * m(idx_C) * m(idx_O) / m(idx_COj)**2
    A(idx_atom_C, idx_atom_N) = A(idx_atom_C, idx_atom_N) +  x(idx_CN) * m(idx_C) * m(idx_N) / m(idx_CN)**2
    A(idx_atom_C, idx_atom_N) = A(idx_atom_C, idx_atom_N) +  x(idx_HCN) * m(idx_C) * m(idx_N) / m(idx_HCN)**2
    A(idx_atom_C, idx_atom_N) = A(idx_atom_C, idx_atom_N) +  x(idx_OCN) * m(idx_C) * m(idx_N) / m(idx_OCN)**2
    A(idx_atom_C, idx_atom_S) = A(idx_atom_C, idx_atom_S) +  x(idx_CS) * m(idx_C) * m(idx_S) / m(idx_CS)**2
    A(idx_atom_C, idx_atom_S) = A(idx_atom_C, idx_atom_S) +  x(idx_OCS) * m(idx_C) * m(idx_S) / m(idx_OCS)**2
    A(idx_atom_E, idx_atom_E) = A(idx_atom_E, idx_atom_E) +  x(idx_E) * m(idx_E) * m(idx_E) / m(idx_E)**2
    A(idx_atom_F, idx_atom_F) = A(idx_atom_F, idx_atom_F) +  x(idx_F) * m(idx_F) * m(idx_F) / m(idx_F)**2
    A(idx_atom_F, idx_atom_F) = A(idx_atom_F, idx_atom_F) +  x(idx_HF) * m(idx_F) * m(idx_F) / m(idx_HF)**2
    A(idx_atom_F, idx_atom_F) = A(idx_atom_F, idx_atom_F) +  x(idx_SIFj) * m(idx_F) * m(idx_F) / m(idx_SIFj)**2
    A(idx_atom_F, idx_atom_F) = A(idx_atom_F, idx_atom_F) +  x(idx_Fj) * m(idx_F) * m(idx_F) / m(idx_Fj)**2
    A(idx_atom_F, idx_atom_H) = A(idx_atom_F, idx_atom_H) +  x(idx_HF) * m(idx_F) * m(idx_H) / m(idx_HF)**2
    A(idx_atom_F, idx_atom_Si) = A(idx_atom_F, idx_atom_Si) +  x(idx_SIFj) * m(idx_F) * m(idx_Si) / m(idx_SIFj)**2
    A(idx_atom_Mg, idx_atom_Mg) = A(idx_atom_Mg, idx_atom_Mg) +  x(idx_MG) * m(idx_Mg) * m(idx_Mg) / m(idx_MG)**2
    A(idx_atom_Mg, idx_atom_Mg) = A(idx_atom_Mg, idx_atom_Mg) +  x(idx_MGj) * m(idx_Mg) * m(idx_Mg) / m(idx_MGj)**2
    A(idx_atom_H, idx_atom_C) = A(idx_atom_H, idx_atom_C) +  x(idx_CH) * m(idx_H) * m(idx_C) / m(idx_CH)**2
    A(idx_atom_H, idx_atom_C) = A(idx_atom_H, idx_atom_C) + 2d0 * x(idx_CH2) * m(idx_H) * m(idx_C) / m(idx_CH2)**2
    A(idx_atom_H, idx_atom_C) = A(idx_atom_H, idx_atom_C) +  x(idx_HCN) * m(idx_H) * m(idx_C) / m(idx_HCN)**2
    A(idx_atom_H, idx_atom_C) = A(idx_atom_H, idx_atom_C) +  x(idx_HCO) * m(idx_H) * m(idx_C) / m(idx_HCO)**2
    A(idx_atom_H, idx_atom_C) = A(idx_atom_H, idx_atom_C) +  x(idx_HCOj) * m(idx_H) * m(idx_C) / m(idx_HCOj)**2
    A(idx_atom_H, idx_atom_C) = A(idx_atom_H, idx_atom_C) +  x(idx_CHj) * m(idx_H) * m(idx_C) / m(idx_CHj)**2
    A(idx_atom_H, idx_atom_F) = A(idx_atom_H, idx_atom_F) +  x(idx_HF) * m(idx_H) * m(idx_F) / m(idx_HF)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_Hk) * m(idx_H) * m(idx_H) / m(idx_Hk)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_H) * m(idx_H) * m(idx_H) / m(idx_H)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) + 4d0 * x(idx_H2) * m(idx_H) * m(idx_H) / m(idx_H2)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_OH) * m(idx_H) * m(idx_H) / m(idx_OH)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_HS) * m(idx_H) * m(idx_H) / m(idx_HS)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) + 4d0 * x(idx_H2S) * m(idx_H) * m(idx_H) / m(idx_H2S)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_HF) * m(idx_H) * m(idx_H) / m(idx_HF)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_CH) * m(idx_H) * m(idx_H) / m(idx_CH)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) + 4d0 * x(idx_CH2) * m(idx_H) * m(idx_H) / m(idx_CH2)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_NH) * m(idx_H) * m(idx_H) / m(idx_NH)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_HCN) * m(idx_H) * m(idx_H) / m(idx_HCN)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) + 4d0 * x(idx_NH2) * m(idx_H) * m(idx_H) / m(idx_NH2)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_HCO) * m(idx_H) * m(idx_H) / m(idx_HCO)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) + 4d0 * x(idx_H2O) * m(idx_H) * m(idx_H) / m(idx_H2O)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_Hj) * m(idx_H) * m(idx_H) / m(idx_Hj)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_NHj) * m(idx_H) * m(idx_H) / m(idx_NHj)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_HSj) * m(idx_H) * m(idx_H) / m(idx_HSj)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_OHj) * m(idx_H) * m(idx_H) / m(idx_OHj)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_HEHj) * m(idx_H) * m(idx_H) / m(idx_HEHj)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) + 4d0 * x(idx_H2j) * m(idx_H) * m(idx_H) / m(idx_H2j)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_SIHj) * m(idx_H) * m(idx_H) / m(idx_SIHj)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_HCOj) * m(idx_H) * m(idx_H) / m(idx_HCOj)**2
    A(idx_atom_H, idx_atom_H) = A(idx_atom_H, idx_atom_H) +  x(idx_CHj) * m(idx_H) * m(idx_H) / m(idx_CHj)**2
    A(idx_atom_H, idx_atom_Si) = A(idx_atom_H, idx_atom_Si) +  x(idx_SIHj) * m(idx_H) * m(idx_Si) / m(idx_SIHj)**2
    A(idx_atom_H, idx_atom_O) = A(idx_atom_H, idx_atom_O) +  x(idx_OH) * m(idx_H) * m(idx_O) / m(idx_OH)**2
    A(idx_atom_H, idx_atom_O) = A(idx_atom_H, idx_atom_O) +  x(idx_HCO) * m(idx_H) * m(idx_O) / m(idx_HCO)**2
    A(idx_atom_H, idx_atom_O) = A(idx_atom_H, idx_atom_O) + 2d0 * x(idx_H2O) * m(idx_H) * m(idx_O) / m(idx_H2O)**2
    A(idx_atom_H, idx_atom_O) = A(idx_atom_H, idx_atom_O) +  x(idx_OHj) * m(idx_H) * m(idx_O) / m(idx_OHj)**2
    A(idx_atom_H, idx_atom_O) = A(idx_atom_H, idx_atom_O) +  x(idx_HCOj) * m(idx_H) * m(idx_O) / m(idx_HCOj)**2
    A(idx_atom_H, idx_atom_N) = A(idx_atom_H, idx_atom_N) +  x(idx_NH) * m(idx_H) * m(idx_N) / m(idx_NH)**2
    A(idx_atom_H, idx_atom_N) = A(idx_atom_H, idx_atom_N) +  x(idx_HCN) * m(idx_H) * m(idx_N) / m(idx_HCN)**2
    A(idx_atom_H, idx_atom_N) = A(idx_atom_H, idx_atom_N) + 2d0 * x(idx_NH2) * m(idx_H) * m(idx_N) / m(idx_NH2)**2
    A(idx_atom_H, idx_atom_N) = A(idx_atom_H, idx_atom_N) +  x(idx_NHj) * m(idx_H) * m(idx_N) / m(idx_NHj)**2
    A(idx_atom_H, idx_atom_S) = A(idx_atom_H, idx_atom_S) +  x(idx_HS) * m(idx_H) * m(idx_S) / m(idx_HS)**2
    A(idx_atom_H, idx_atom_S) = A(idx_atom_H, idx_atom_S) + 2d0 * x(idx_H2S) * m(idx_H) * m(idx_S) / m(idx_H2S)**2
    A(idx_atom_H, idx_atom_S) = A(idx_atom_H, idx_atom_S) +  x(idx_HSj) * m(idx_H) * m(idx_S) / m(idx_HSj)**2
    A(idx_atom_H, idx_atom_He) = A(idx_atom_H, idx_atom_He) +  x(idx_HEHj) * m(idx_H) * m(idx_He) / m(idx_HEHj)**2
    A(idx_atom_Si, idx_atom_F) = A(idx_atom_Si, idx_atom_F) +  x(idx_SIFj) * m(idx_Si) * m(idx_F) / m(idx_SIFj)**2
    A(idx_atom_Si, idx_atom_H) = A(idx_atom_Si, idx_atom_H) +  x(idx_SIHj) * m(idx_Si) * m(idx_H) / m(idx_SIHj)**2
    A(idx_atom_Si, idx_atom_Si) = A(idx_atom_Si, idx_atom_Si) +  x(idx_SI) * m(idx_Si) * m(idx_Si) / m(idx_SI)**2
    A(idx_atom_Si, idx_atom_Si) = A(idx_atom_Si, idx_atom_Si) +  x(idx_SIO) * m(idx_Si) * m(idx_Si) / m(idx_SIO)**2
    A(idx_atom_Si, idx_atom_Si) = A(idx_atom_Si, idx_atom_Si) +  x(idx_SIO2) * m(idx_Si) * m(idx_Si) / m(idx_SIO2)**2
    A(idx_atom_Si, idx_atom_Si) = A(idx_atom_Si, idx_atom_Si) +  x(idx_SIj) * m(idx_Si) * m(idx_Si) / m(idx_SIj)**2
    A(idx_atom_Si, idx_atom_Si) = A(idx_atom_Si, idx_atom_Si) +  x(idx_SIHj) * m(idx_Si) * m(idx_Si) / m(idx_SIHj)**2
    A(idx_atom_Si, idx_atom_Si) = A(idx_atom_Si, idx_atom_Si) +  x(idx_SIOj) * m(idx_Si) * m(idx_Si) / m(idx_SIOj)**2
    A(idx_atom_Si, idx_atom_Si) = A(idx_atom_Si, idx_atom_Si) +  x(idx_SIFj) * m(idx_Si) * m(idx_Si) / m(idx_SIFj)**2
    A(idx_atom_Si, idx_atom_O) = A(idx_atom_Si, idx_atom_O) +  x(idx_SIO) * m(idx_Si) * m(idx_O) / m(idx_SIO)**2
    A(idx_atom_Si, idx_atom_O) = A(idx_atom_Si, idx_atom_O) + 2d0 * x(idx_SIO2) * m(idx_Si) * m(idx_O) / m(idx_SIO2)**2
    A(idx_atom_Si, idx_atom_O) = A(idx_atom_Si, idx_atom_O) +  x(idx_SIOj) * m(idx_Si) * m(idx_O) / m(idx_SIOj)**2
    A(idx_atom_O, idx_atom_C) = A(idx_atom_O, idx_atom_C) +  x(idx_CO) * m(idx_O) * m(idx_C) / m(idx_CO)**2
    A(idx_atom_O, idx_atom_C) = A(idx_atom_O, idx_atom_C) + 2d0 * x(idx_CO2) * m(idx_O) * m(idx_C) / m(idx_CO2)**2
    A(idx_atom_O, idx_atom_C) = A(idx_atom_O, idx_atom_C) +  x(idx_OCN) * m(idx_O) * m(idx_C) / m(idx_OCN)**2
    A(idx_atom_O, idx_atom_C) = A(idx_atom_O, idx_atom_C) +  x(idx_HCO) * m(idx_O) * m(idx_C) / m(idx_HCO)**2
    A(idx_atom_O, idx_atom_C) = A(idx_atom_O, idx_atom_C) +  x(idx_OCS) * m(idx_O) * m(idx_C) / m(idx_OCS)**2
    A(idx_atom_O, idx_atom_C) = A(idx_atom_O, idx_atom_C) +  x(idx_HCOj) * m(idx_O) * m(idx_C) / m(idx_HCOj)**2
    A(idx_atom_O, idx_atom_C) = A(idx_atom_O, idx_atom_C) +  x(idx_COj) * m(idx_O) * m(idx_C) / m(idx_COj)**2
    A(idx_atom_O, idx_atom_H) = A(idx_atom_O, idx_atom_H) +  x(idx_OH) * m(idx_O) * m(idx_H) / m(idx_OH)**2
    A(idx_atom_O, idx_atom_H) = A(idx_atom_O, idx_atom_H) +  x(idx_HCO) * m(idx_O) * m(idx_H) / m(idx_HCO)**2
    A(idx_atom_O, idx_atom_H) = A(idx_atom_O, idx_atom_H) + 2d0 * x(idx_H2O) * m(idx_O) * m(idx_H) / m(idx_H2O)**2
    A(idx_atom_O, idx_atom_H) = A(idx_atom_O, idx_atom_H) +  x(idx_OHj) * m(idx_O) * m(idx_H) / m(idx_OHj)**2
    A(idx_atom_O, idx_atom_H) = A(idx_atom_O, idx_atom_H) +  x(idx_HCOj) * m(idx_O) * m(idx_H) / m(idx_HCOj)**2
    A(idx_atom_O, idx_atom_Si) = A(idx_atom_O, idx_atom_Si) +  x(idx_SIO) * m(idx_O) * m(idx_Si) / m(idx_SIO)**2
    A(idx_atom_O, idx_atom_Si) = A(idx_atom_O, idx_atom_Si) + 2d0 * x(idx_SIO2) * m(idx_O) * m(idx_Si) / m(idx_SIO2)**2
    A(idx_atom_O, idx_atom_Si) = A(idx_atom_O, idx_atom_Si) +  x(idx_SIOj) * m(idx_O) * m(idx_Si) / m(idx_SIOj)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_Ok) * m(idx_O) * m(idx_O) / m(idx_Ok)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_NO) * m(idx_O) * m(idx_O) / m(idx_NO)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_CO) * m(idx_O) * m(idx_O) / m(idx_CO)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) + 4d0 * x(idx_O2) * m(idx_O) * m(idx_O) / m(idx_O2)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_SO) * m(idx_O) * m(idx_O) / m(idx_SO)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_O) * m(idx_O) * m(idx_O) / m(idx_O)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_OH) * m(idx_O) * m(idx_O) / m(idx_OH)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) + 4d0 * x(idx_SO2) * m(idx_O) * m(idx_O) / m(idx_SO2)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) + 4d0 * x(idx_CO2) * m(idx_O) * m(idx_O) / m(idx_CO2)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_SIO) * m(idx_O) * m(idx_O) / m(idx_SIO)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) + 4d0 * x(idx_SIO2) * m(idx_O) * m(idx_O) / m(idx_SIO2)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_OCN) * m(idx_O) * m(idx_O) / m(idx_OCN)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_HCO) * m(idx_O) * m(idx_O) / m(idx_HCO)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_H2O) * m(idx_O) * m(idx_O) / m(idx_H2O)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_OCS) * m(idx_O) * m(idx_O) / m(idx_OCS)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_PO) * m(idx_O) * m(idx_O) / m(idx_PO)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_OHj) * m(idx_O) * m(idx_O) / m(idx_OHj)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_HCOj) * m(idx_O) * m(idx_O) / m(idx_HCOj)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_Oj) * m(idx_O) * m(idx_O) / m(idx_Oj)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_SIOj) * m(idx_O) * m(idx_O) / m(idx_SIOj)**2
    A(idx_atom_O, idx_atom_O) = A(idx_atom_O, idx_atom_O) +  x(idx_COj) * m(idx_O) * m(idx_O) / m(idx_COj)**2
    A(idx_atom_O, idx_atom_N) = A(idx_atom_O, idx_atom_N) +  x(idx_NO) * m(idx_O) * m(idx_N) / m(idx_NO)**2
    A(idx_atom_O, idx_atom_N) = A(idx_atom_O, idx_atom_N) +  x(idx_OCN) * m(idx_O) * m(idx_N) / m(idx_OCN)**2
    A(idx_atom_O, idx_atom_P) = A(idx_atom_O, idx_atom_P) +  x(idx_PO) * m(idx_O) * m(idx_P) / m(idx_PO)**2
    A(idx_atom_O, idx_atom_S) = A(idx_atom_O, idx_atom_S) +  x(idx_SO) * m(idx_O) * m(idx_S) / m(idx_SO)**2
    A(idx_atom_O, idx_atom_S) = A(idx_atom_O, idx_atom_S) + 2d0 * x(idx_SO2) * m(idx_O) * m(idx_S) / m(idx_SO2)**2
    A(idx_atom_O, idx_atom_S) = A(idx_atom_O, idx_atom_S) +  x(idx_OCS) * m(idx_O) * m(idx_S) / m(idx_OCS)**2
    A(idx_atom_N, idx_atom_C) = A(idx_atom_N, idx_atom_C) +  x(idx_CN) * m(idx_N) * m(idx_C) / m(idx_CN)**2
    A(idx_atom_N, idx_atom_C) = A(idx_atom_N, idx_atom_C) +  x(idx_HCN) * m(idx_N) * m(idx_C) / m(idx_HCN)**2
    A(idx_atom_N, idx_atom_C) = A(idx_atom_N, idx_atom_C) +  x(idx_OCN) * m(idx_N) * m(idx_C) / m(idx_OCN)**2
    A(idx_atom_N, idx_atom_H) = A(idx_atom_N, idx_atom_H) +  x(idx_NH) * m(idx_N) * m(idx_H) / m(idx_NH)**2
    A(idx_atom_N, idx_atom_H) = A(idx_atom_N, idx_atom_H) +  x(idx_HCN) * m(idx_N) * m(idx_H) / m(idx_HCN)**2
    A(idx_atom_N, idx_atom_H) = A(idx_atom_N, idx_atom_H) + 2d0 * x(idx_NH2) * m(idx_N) * m(idx_H) / m(idx_NH2)**2
    A(idx_atom_N, idx_atom_H) = A(idx_atom_N, idx_atom_H) +  x(idx_NHj) * m(idx_N) * m(idx_H) / m(idx_NHj)**2
    A(idx_atom_N, idx_atom_O) = A(idx_atom_N, idx_atom_O) +  x(idx_NO) * m(idx_N) * m(idx_O) / m(idx_NO)**2
    A(idx_atom_N, idx_atom_O) = A(idx_atom_N, idx_atom_O) +  x(idx_OCN) * m(idx_N) * m(idx_O) / m(idx_OCN)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) +  x(idx_NO) * m(idx_N) * m(idx_N) / m(idx_NO)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) +  x(idx_N) * m(idx_N) * m(idx_N) / m(idx_N)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) +  x(idx_NS) * m(idx_N) * m(idx_N) / m(idx_NS)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) +  x(idx_CN) * m(idx_N) * m(idx_N) / m(idx_CN)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) + 4d0 * x(idx_N2) * m(idx_N) * m(idx_N) / m(idx_N2)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) +  x(idx_NH) * m(idx_N) * m(idx_N) / m(idx_NH)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) +  x(idx_HCN) * m(idx_N) * m(idx_N) / m(idx_HCN)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) +  x(idx_NH2) * m(idx_N) * m(idx_N) / m(idx_NH2)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) +  x(idx_OCN) * m(idx_N) * m(idx_N) / m(idx_OCN)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) +  x(idx_PN) * m(idx_N) * m(idx_N) / m(idx_PN)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) +  x(idx_NHj) * m(idx_N) * m(idx_N) / m(idx_NHj)**2
    A(idx_atom_N, idx_atom_N) = A(idx_atom_N, idx_atom_N) +  x(idx_Nj) * m(idx_N) * m(idx_N) / m(idx_Nj)**2
    A(idx_atom_N, idx_atom_P) = A(idx_atom_N, idx_atom_P) +  x(idx_PN) * m(idx_N) * m(idx_P) / m(idx_PN)**2
    A(idx_atom_N, idx_atom_S) = A(idx_atom_N, idx_atom_S) +  x(idx_NS) * m(idx_N) * m(idx_S) / m(idx_NS)**2
    A(idx_atom_P, idx_atom_O) = A(idx_atom_P, idx_atom_O) +  x(idx_PO) * m(idx_P) * m(idx_O) / m(idx_PO)**2
    A(idx_atom_P, idx_atom_N) = A(idx_atom_P, idx_atom_N) +  x(idx_PN) * m(idx_P) * m(idx_N) / m(idx_PN)**2
    A(idx_atom_P, idx_atom_P) = A(idx_atom_P, idx_atom_P) +  x(idx_P) * m(idx_P) * m(idx_P) / m(idx_P)**2
    A(idx_atom_P, idx_atom_P) = A(idx_atom_P, idx_atom_P) +  x(idx_PN) * m(idx_P) * m(idx_P) / m(idx_PN)**2
    A(idx_atom_P, idx_atom_P) = A(idx_atom_P, idx_atom_P) +  x(idx_PO) * m(idx_P) * m(idx_P) / m(idx_PO)**2
    A(idx_atom_P, idx_atom_P) = A(idx_atom_P, idx_atom_P) +  x(idx_Pj) * m(idx_P) * m(idx_P) / m(idx_Pj)**2
    A(idx_atom_S, idx_atom_C) = A(idx_atom_S, idx_atom_C) +  x(idx_CS) * m(idx_S) * m(idx_C) / m(idx_CS)**2
    A(idx_atom_S, idx_atom_C) = A(idx_atom_S, idx_atom_C) +  x(idx_OCS) * m(idx_S) * m(idx_C) / m(idx_OCS)**2
    A(idx_atom_S, idx_atom_H) = A(idx_atom_S, idx_atom_H) +  x(idx_HS) * m(idx_S) * m(idx_H) / m(idx_HS)**2
    A(idx_atom_S, idx_atom_H) = A(idx_atom_S, idx_atom_H) + 2d0 * x(idx_H2S) * m(idx_S) * m(idx_H) / m(idx_H2S)**2
    A(idx_atom_S, idx_atom_H) = A(idx_atom_S, idx_atom_H) +  x(idx_HSj) * m(idx_S) * m(idx_H) / m(idx_HSj)**2
    A(idx_atom_S, idx_atom_O) = A(idx_atom_S, idx_atom_O) +  x(idx_SO) * m(idx_S) * m(idx_O) / m(idx_SO)**2
    A(idx_atom_S, idx_atom_O) = A(idx_atom_S, idx_atom_O) + 2d0 * x(idx_SO2) * m(idx_S) * m(idx_O) / m(idx_SO2)**2
    A(idx_atom_S, idx_atom_O) = A(idx_atom_S, idx_atom_O) +  x(idx_OCS) * m(idx_S) * m(idx_O) / m(idx_OCS)**2
    A(idx_atom_S, idx_atom_N) = A(idx_atom_S, idx_atom_N) +  x(idx_NS) * m(idx_S) * m(idx_N) / m(idx_NS)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) +  x(idx_Sk) * m(idx_S) * m(idx_S) / m(idx_Sk)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) +  x(idx_S) * m(idx_S) * m(idx_S) / m(idx_S)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) +  x(idx_SO) * m(idx_S) * m(idx_S) / m(idx_SO)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) +  x(idx_HS) * m(idx_S) * m(idx_S) / m(idx_HS)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) +  x(idx_NS) * m(idx_S) * m(idx_S) / m(idx_NS)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) +  x(idx_H2S) * m(idx_S) * m(idx_S) / m(idx_H2S)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) +  x(idx_CS) * m(idx_S) * m(idx_S) / m(idx_CS)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) + 4d0 * x(idx_S2) * m(idx_S) * m(idx_S) / m(idx_S2)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) +  x(idx_SO2) * m(idx_S) * m(idx_S) / m(idx_SO2)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) +  x(idx_OCS) * m(idx_S) * m(idx_S) / m(idx_OCS)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) +  x(idx_HSj) * m(idx_S) * m(idx_S) / m(idx_HSj)**2
    A(idx_atom_S, idx_atom_S) = A(idx_atom_S, idx_atom_S) +  x(idx_Sj) * m(idx_S) * m(idx_S) / m(idx_Sj)**2
    A(idx_atom_Fe, idx_atom_Fe) = A(idx_atom_Fe, idx_atom_Fe) +  x(idx_FE) * m(idx_Fe) * m(idx_Fe) / m(idx_FE)**2
    A(idx_atom_Fe, idx_atom_Fe) = A(idx_atom_Fe, idx_atom_Fe) +  x(idx_FEj) * m(idx_Fe) * m(idx_Fe) / m(idx_FEj)**2
    A(idx_atom_Na, idx_atom_Na) = A(idx_atom_Na, idx_atom_Na) +  x(idx_NA) * m(idx_Na) * m(idx_Na) / m(idx_NA)**2
    A(idx_atom_Na, idx_atom_Na) = A(idx_atom_Na, idx_atom_Na) +  x(idx_NAj) * m(idx_Na) * m(idx_Na) / m(idx_NAj)**2
    A(idx_atom_He, idx_atom_H) = A(idx_atom_He, idx_atom_H) +  x(idx_HEHj) * m(idx_He) * m(idx_H) / m(idx_HEHj)**2
    A(idx_atom_He, idx_atom_He) = A(idx_atom_He, idx_atom_He) +  x(idx_HE) * m(idx_He) * m(idx_He) / m(idx_HE)**2
    A(idx_atom_He, idx_atom_He) = A(idx_atom_He, idx_atom_He) +  x(idx_HEj) * m(idx_He) * m(idx_He) / m(idx_HEj)**2
    A(idx_atom_He, idx_atom_He) = A(idx_atom_He, idx_atom_He) +  x(idx_HEHj) * m(idx_He) * m(idx_He) / m(idx_HEHj)**2

    B(:) = ref(:)

    call mydgesv(natoms,A(:,:),B(:), "conserveLin_x")

    x(idx_Ok) = x(idx_Ok) * (m(idx_O) * B(idx_atom_O))/m(idx_Ok)
    x(idx_Hk) = x(idx_Hk) * (m(idx_H) * B(idx_atom_H))/m(idx_Hk)
    x(idx_Sk) = x(idx_Sk) * (m(idx_S) * B(idx_atom_S))/m(idx_Sk)
    x(idx_Ck) = x(idx_Ck) * (m(idx_C) * B(idx_atom_C))/m(idx_Ck)
    x(idx_H) = x(idx_H) * (m(idx_H) * B(idx_atom_H))/m(idx_H)
    x(idx_HE) = x(idx_HE) * (m(idx_He) * B(idx_atom_He))/m(idx_HE)
    x(idx_C) = x(idx_C) * (m(idx_C) * B(idx_atom_C))/m(idx_C)
    x(idx_NO) = x(idx_NO) * (m(idx_O) * B(idx_atom_O) + &
        m(idx_N) * B(idx_atom_N))/m(idx_NO)
    x(idx_CO) = x(idx_CO) * (m(idx_C) * B(idx_atom_C) + &
        m(idx_O) * B(idx_atom_O))/m(idx_CO)
    x(idx_N) = x(idx_N) * (m(idx_N) * B(idx_atom_N))/m(idx_N)
    x(idx_O2) = x(idx_O2) * (2d0*m(idx_O) * B(idx_atom_O))/m(idx_O2)
    x(idx_S) = x(idx_S) * (m(idx_S) * B(idx_atom_S))/m(idx_S)
    x(idx_SO) = x(idx_SO) * (m(idx_S) * B(idx_atom_S) + &
        m(idx_O) * B(idx_atom_O))/m(idx_SO)
    x(idx_O) = x(idx_O) * (m(idx_O) * B(idx_atom_O))/m(idx_O)
    x(idx_H2) = x(idx_H2) * (2d0*m(idx_H) * B(idx_atom_H))/m(idx_H2)
    x(idx_SI) = x(idx_SI) * (m(idx_Si) * B(idx_atom_Si))/m(idx_SI)
    x(idx_OH) = x(idx_OH) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_O) * B(idx_atom_O))/m(idx_OH)
    x(idx_HS) = x(idx_HS) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_S) * B(idx_atom_S))/m(idx_HS)
    x(idx_NS) = x(idx_NS) * (m(idx_S) * B(idx_atom_S) + &
        m(idx_N) * B(idx_atom_N))/m(idx_NS)
    x(idx_H2S) = x(idx_H2S) * (2d0*m(idx_H) * B(idx_atom_H) + &
        m(idx_S) * B(idx_atom_S))/m(idx_H2S)
    x(idx_FE) = x(idx_FE) * (m(idx_Fe) * B(idx_atom_Fe))/m(idx_FE)
    x(idx_CS) = x(idx_CS) * (m(idx_C) * B(idx_atom_C) + &
        m(idx_S) * B(idx_atom_S))/m(idx_CS)
    x(idx_CN) = x(idx_CN) * (m(idx_C) * B(idx_atom_C) + &
        m(idx_N) * B(idx_atom_N))/m(idx_CN)
    x(idx_S2) = x(idx_S2) * (2d0*m(idx_S) * B(idx_atom_S))/m(idx_S2)
    x(idx_NA) = x(idx_NA) * (m(idx_Na) * B(idx_atom_Na))/m(idx_NA)
    x(idx_F) = x(idx_F) * (m(idx_F) * B(idx_atom_F))/m(idx_F)
    x(idx_HF) = x(idx_HF) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_F) * B(idx_atom_F))/m(idx_HF)
    x(idx_CH) = x(idx_CH) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_C) * B(idx_atom_C))/m(idx_CH)
    x(idx_SO2) = x(idx_SO2) * (m(idx_S) * B(idx_atom_S) + &
        2d0*m(idx_O) * B(idx_atom_O))/m(idx_SO2)
    x(idx_C2) = x(idx_C2) * (2d0*m(idx_C) * B(idx_atom_C))/m(idx_C2)
    x(idx_N2) = x(idx_N2) * (2d0*m(idx_N) * B(idx_atom_N))/m(idx_N2)
    x(idx_CH2) = x(idx_CH2) * (2d0*m(idx_H) * B(idx_atom_H) + &
        m(idx_C) * B(idx_atom_C))/m(idx_CH2)
    x(idx_NH) = x(idx_NH) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_N) * B(idx_atom_N))/m(idx_NH)
    x(idx_HCN) = x(idx_HCN) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_C) * B(idx_atom_C) + &
        m(idx_N) * B(idx_atom_N))/m(idx_HCN)
    x(idx_CO2) = x(idx_CO2) * (m(idx_C) * B(idx_atom_C) + &
        2d0*m(idx_O) * B(idx_atom_O))/m(idx_CO2)
    x(idx_SIO) = x(idx_SIO) * (m(idx_Si) * B(idx_atom_Si) + &
        m(idx_O) * B(idx_atom_O))/m(idx_SIO)
    x(idx_SIO2) = x(idx_SIO2) * (m(idx_Si) * B(idx_atom_Si) + &
        2d0*m(idx_O) * B(idx_atom_O))/m(idx_SIO2)
    x(idx_NH2) = x(idx_NH2) * (2d0*m(idx_H) * B(idx_atom_H) + &
        m(idx_N) * B(idx_atom_N))/m(idx_NH2)
    x(idx_OCN) = x(idx_OCN) * (m(idx_C) * B(idx_atom_C) + &
        m(idx_O) * B(idx_atom_O) + &
        m(idx_N) * B(idx_atom_N))/m(idx_OCN)
    x(idx_MG) = x(idx_MG) * (m(idx_Mg) * B(idx_atom_Mg))/m(idx_MG)
    x(idx_P) = x(idx_P) * (m(idx_P) * B(idx_atom_P))/m(idx_P)
    x(idx_HCO) = x(idx_HCO) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_C) * B(idx_atom_C) + &
        m(idx_O) * B(idx_atom_O))/m(idx_HCO)
    x(idx_H2O) = x(idx_H2O) * (2d0*m(idx_H) * B(idx_atom_H) + &
        m(idx_O) * B(idx_atom_O))/m(idx_H2O)
    x(idx_OCS) = x(idx_OCS) * (m(idx_C) * B(idx_atom_C) + &
        m(idx_S) * B(idx_atom_S) + &
        m(idx_O) * B(idx_atom_O))/m(idx_OCS)
    x(idx_PN) = x(idx_PN) * (m(idx_P) * B(idx_atom_P) + &
        m(idx_N) * B(idx_atom_N))/m(idx_PN)
    x(idx_PO) = x(idx_PO) * (m(idx_P) * B(idx_atom_P) + &
        m(idx_O) * B(idx_atom_O))/m(idx_PO)
    x(idx_HEj) = x(idx_HEj) * (m(idx_He) * B(idx_atom_He))/m(idx_HEj)
    x(idx_Hj) = x(idx_Hj) * (m(idx_H) * B(idx_atom_H))/m(idx_Hj)
    x(idx_NHj) = x(idx_NHj) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_N) * B(idx_atom_N))/m(idx_NHj)
    x(idx_HSj) = x(idx_HSj) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_S) * B(idx_atom_S))/m(idx_HSj)
    x(idx_Sj) = x(idx_Sj) * (m(idx_S) * B(idx_atom_S))/m(idx_Sj)
    x(idx_SIj) = x(idx_SIj) * (m(idx_Si) * B(idx_atom_Si))/m(idx_SIj)
    x(idx_OHj) = x(idx_OHj) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_O) * B(idx_atom_O))/m(idx_OHj)
    x(idx_HEHj) = x(idx_HEHj) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_He) * B(idx_atom_He))/m(idx_HEHj)
    x(idx_H2j) = x(idx_H2j) * (2d0*m(idx_H) * B(idx_atom_H))/m(idx_H2j)
    x(idx_FEj) = x(idx_FEj) * (m(idx_Fe) * B(idx_atom_Fe))/m(idx_FEj)
    x(idx_SIHj) = x(idx_SIHj) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_Si) * B(idx_atom_Si))/m(idx_SIHj)
    x(idx_NAj) = x(idx_NAj) * (m(idx_Na) * B(idx_atom_Na))/m(idx_NAj)
    x(idx_HCOj) = x(idx_HCOj) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_C) * B(idx_atom_C) + &
        m(idx_O) * B(idx_atom_O))/m(idx_HCOj)
    x(idx_CHj) = x(idx_CHj) * (m(idx_H) * B(idx_atom_H) + &
        m(idx_C) * B(idx_atom_C))/m(idx_CHj)
    x(idx_Oj) = x(idx_Oj) * (m(idx_O) * B(idx_atom_O))/m(idx_Oj)
    x(idx_MGj) = x(idx_MGj) * (m(idx_Mg) * B(idx_atom_Mg))/m(idx_MGj)
    x(idx_SIOj) = x(idx_SIOj) * (m(idx_Si) * B(idx_atom_Si) + &
        m(idx_O) * B(idx_atom_O))/m(idx_SIOj)
    x(idx_Pj) = x(idx_Pj) * (m(idx_P) * B(idx_atom_P))/m(idx_Pj)
    x(idx_SIFj) = x(idx_SIFj) * (m(idx_Si) * B(idx_atom_Si) + &
        m(idx_F) * B(idx_atom_F))/m(idx_SIFj)
    x(idx_Cj) = x(idx_Cj) * (m(idx_C) * B(idx_atom_C))/m(idx_Cj)
    x(idx_Nj) = x(idx_Nj) * (m(idx_N) * B(idx_atom_N))/m(idx_Nj)
    x(idx_COj) = x(idx_COj) * (m(idx_C) * B(idx_atom_C) + &
        m(idx_O) * B(idx_atom_O))/m(idx_COj)
    x(idx_Fj) = x(idx_Fj) * (m(idx_F) * B(idx_atom_F))/m(idx_Fj)

    !charge conservation
    x(idx_E) = m(idx_E)*(- x(idx_Ok) / m(idx_Ok) &
        - x(idx_Hk) / m(idx_Hk) &
        - x(idx_Sk) / m(idx_Sk) &
        - x(idx_Ck) / m(idx_Ck) &
        + x(idx_HEj) / m(idx_HEj) &
        + x(idx_Hj) / m(idx_Hj) &
        + x(idx_NHj) / m(idx_NHj) &
        + x(idx_HSj) / m(idx_HSj) &
        + x(idx_Sj) / m(idx_Sj) &
        + x(idx_SIj) / m(idx_SIj) &
        + x(idx_OHj) / m(idx_OHj) &
        + x(idx_HEHj) / m(idx_HEHj) &
        + x(idx_H2j) / m(idx_H2j) &
        + x(idx_FEj) / m(idx_FEj) &
        + x(idx_SIHj) / m(idx_SIHj) &
        + x(idx_NAj) / m(idx_NAj) &
        + x(idx_HCOj) / m(idx_HCOj) &
        + x(idx_CHj) / m(idx_CHj) &
        + x(idx_Oj) / m(idx_Oj) &
        + x(idx_MGj) / m(idx_MGj) &
        + x(idx_SIOj) / m(idx_SIOj) &
        + x(idx_Pj) / m(idx_Pj) &
        + x(idx_SIFj) / m(idx_SIFj) &
        + x(idx_Cj) / m(idx_Cj) &
        + x(idx_Nj) / m(idx_Nj) &
        + x(idx_COj) / m(idx_COj) &
        + x(idx_Fj) / m(idx_Fj))
    !check if charge conservation goes wrong
    if(x(idx_E)<0d0) then
      print *,"ERROR in conserveLin, electrons < 0"
      stop
    end if

  end subroutine conserveLin_x

  !***************************
  !compute the total reference mass atom type by atom type
  function conserveLinGetRef_x(x)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::conserveLinGetRef_x(natoms),x(nmols)
    real*8::m(nspec)

    m(:) = get_mass()
    conserveLinGetRef_x(:) = 0d0

    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_Ok)/m(idx_Ok)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_Hk)/m(idx_Hk)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + m(idx_S)*x(idx_Sk)/m(idx_Sk)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_Ck)/m(idx_Ck)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_H)/m(idx_H)
    conserveLinGetRef_x(idx_atom_He) = conserveLinGetRef_x(idx_atom_He) + m(idx_He)*x(idx_HE)/m(idx_HE)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_C)/m(idx_C)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_NO)/m(idx_NO)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + m(idx_N)*x(idx_NO)/m(idx_NO)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_CO)/m(idx_CO)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_CO)/m(idx_CO)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + m(idx_N)*x(idx_N)/m(idx_N)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + 2d0*m(idx_O)*x(idx_O2)/m(idx_O2)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + m(idx_S)*x(idx_S)/m(idx_S)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + m(idx_S)*x(idx_SO)/m(idx_SO)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_SO)/m(idx_SO)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_O)/m(idx_O)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + 2d0*m(idx_H)*x(idx_H2)/m(idx_H2)
    conserveLinGetRef_x(idx_atom_Si) = conserveLinGetRef_x(idx_atom_Si) + m(idx_Si)*x(idx_SI)/m(idx_SI)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_OH)/m(idx_OH)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_OH)/m(idx_OH)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_HS)/m(idx_HS)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + m(idx_S)*x(idx_HS)/m(idx_HS)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + m(idx_S)*x(idx_NS)/m(idx_NS)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + m(idx_N)*x(idx_NS)/m(idx_NS)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + 2d0*m(idx_H)*x(idx_H2S)/m(idx_H2S)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + m(idx_S)*x(idx_H2S)/m(idx_H2S)
    conserveLinGetRef_x(idx_atom_Fe) = conserveLinGetRef_x(idx_atom_Fe) + m(idx_Fe)*x(idx_FE)/m(idx_FE)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_CS)/m(idx_CS)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + m(idx_S)*x(idx_CS)/m(idx_CS)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_CN)/m(idx_CN)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + m(idx_N)*x(idx_CN)/m(idx_CN)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + 2d0*m(idx_S)*x(idx_S2)/m(idx_S2)
    conserveLinGetRef_x(idx_atom_Na) = conserveLinGetRef_x(idx_atom_Na) + m(idx_Na)*x(idx_NA)/m(idx_NA)
    conserveLinGetRef_x(idx_atom_F) = conserveLinGetRef_x(idx_atom_F) + m(idx_F)*x(idx_F)/m(idx_F)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_HF)/m(idx_HF)
    conserveLinGetRef_x(idx_atom_F) = conserveLinGetRef_x(idx_atom_F) + m(idx_F)*x(idx_HF)/m(idx_HF)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_CH)/m(idx_CH)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_CH)/m(idx_CH)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + m(idx_S)*x(idx_SO2)/m(idx_SO2)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + 2d0*m(idx_O)*x(idx_SO2)/m(idx_SO2)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + 2d0*m(idx_C)*x(idx_C2)/m(idx_C2)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + 2d0*m(idx_N)*x(idx_N2)/m(idx_N2)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + 2d0*m(idx_H)*x(idx_CH2)/m(idx_CH2)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_CH2)/m(idx_CH2)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_NH)/m(idx_NH)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + m(idx_N)*x(idx_NH)/m(idx_NH)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_HCN)/m(idx_HCN)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_HCN)/m(idx_HCN)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + m(idx_N)*x(idx_HCN)/m(idx_HCN)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_CO2)/m(idx_CO2)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + 2d0*m(idx_O)*x(idx_CO2)/m(idx_CO2)
    conserveLinGetRef_x(idx_atom_Si) = conserveLinGetRef_x(idx_atom_Si) + m(idx_Si)*x(idx_SIO)/m(idx_SIO)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_SIO)/m(idx_SIO)
    conserveLinGetRef_x(idx_atom_Si) = conserveLinGetRef_x(idx_atom_Si) + m(idx_Si)*x(idx_SIO2)/m(idx_SIO2)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + 2d0*m(idx_O)*x(idx_SIO2)/m(idx_SIO2)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + 2d0*m(idx_H)*x(idx_NH2)/m(idx_NH2)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + m(idx_N)*x(idx_NH2)/m(idx_NH2)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_OCN)/m(idx_OCN)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_OCN)/m(idx_OCN)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + m(idx_N)*x(idx_OCN)/m(idx_OCN)
    conserveLinGetRef_x(idx_atom_Mg) = conserveLinGetRef_x(idx_atom_Mg) + m(idx_Mg)*x(idx_MG)/m(idx_MG)
    conserveLinGetRef_x(idx_atom_P) = conserveLinGetRef_x(idx_atom_P) + m(idx_P)*x(idx_P)/m(idx_P)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_HCO)/m(idx_HCO)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_HCO)/m(idx_HCO)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_HCO)/m(idx_HCO)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + 2d0*m(idx_H)*x(idx_H2O)/m(idx_H2O)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_H2O)/m(idx_H2O)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_OCS)/m(idx_OCS)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + m(idx_S)*x(idx_OCS)/m(idx_OCS)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_OCS)/m(idx_OCS)
    conserveLinGetRef_x(idx_atom_P) = conserveLinGetRef_x(idx_atom_P) + m(idx_P)*x(idx_PN)/m(idx_PN)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + m(idx_N)*x(idx_PN)/m(idx_PN)
    conserveLinGetRef_x(idx_atom_P) = conserveLinGetRef_x(idx_atom_P) + m(idx_P)*x(idx_PO)/m(idx_PO)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_PO)/m(idx_PO)
    conserveLinGetRef_x(idx_atom_He) = conserveLinGetRef_x(idx_atom_He) + m(idx_He)*x(idx_HEj)/m(idx_HEj)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_Hj)/m(idx_Hj)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_NHj)/m(idx_NHj)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + m(idx_N)*x(idx_NHj)/m(idx_NHj)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_HSj)/m(idx_HSj)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + m(idx_S)*x(idx_HSj)/m(idx_HSj)
    conserveLinGetRef_x(idx_atom_S) = conserveLinGetRef_x(idx_atom_S) + m(idx_S)*x(idx_Sj)/m(idx_Sj)
    conserveLinGetRef_x(idx_atom_Si) = conserveLinGetRef_x(idx_atom_Si) + m(idx_Si)*x(idx_SIj)/m(idx_SIj)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_OHj)/m(idx_OHj)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_OHj)/m(idx_OHj)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_HEHj)/m(idx_HEHj)
    conserveLinGetRef_x(idx_atom_He) = conserveLinGetRef_x(idx_atom_He) + m(idx_He)*x(idx_HEHj)/m(idx_HEHj)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + 2d0*m(idx_H)*x(idx_H2j)/m(idx_H2j)
    conserveLinGetRef_x(idx_atom_Fe) = conserveLinGetRef_x(idx_atom_Fe) + m(idx_Fe)*x(idx_FEj)/m(idx_FEj)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_SIHj)/m(idx_SIHj)
    conserveLinGetRef_x(idx_atom_Si) = conserveLinGetRef_x(idx_atom_Si) + m(idx_Si)*x(idx_SIHj)/m(idx_SIHj)
    conserveLinGetRef_x(idx_atom_Na) = conserveLinGetRef_x(idx_atom_Na) + m(idx_Na)*x(idx_NAj)/m(idx_NAj)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_HCOj)/m(idx_HCOj)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_HCOj)/m(idx_HCOj)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_HCOj)/m(idx_HCOj)
    conserveLinGetRef_x(idx_atom_H) = conserveLinGetRef_x(idx_atom_H) + m(idx_H)*x(idx_CHj)/m(idx_CHj)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_CHj)/m(idx_CHj)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_Oj)/m(idx_Oj)
    conserveLinGetRef_x(idx_atom_Mg) = conserveLinGetRef_x(idx_atom_Mg) + m(idx_Mg)*x(idx_MGj)/m(idx_MGj)
    conserveLinGetRef_x(idx_atom_Si) = conserveLinGetRef_x(idx_atom_Si) + m(idx_Si)*x(idx_SIOj)/m(idx_SIOj)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_SIOj)/m(idx_SIOj)
    conserveLinGetRef_x(idx_atom_P) = conserveLinGetRef_x(idx_atom_P) + m(idx_P)*x(idx_Pj)/m(idx_Pj)
    conserveLinGetRef_x(idx_atom_Si) = conserveLinGetRef_x(idx_atom_Si) + m(idx_Si)*x(idx_SIFj)/m(idx_SIFj)
    conserveLinGetRef_x(idx_atom_F) = conserveLinGetRef_x(idx_atom_F) + m(idx_F)*x(idx_SIFj)/m(idx_SIFj)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_Cj)/m(idx_Cj)
    conserveLinGetRef_x(idx_atom_N) = conserveLinGetRef_x(idx_atom_N) + m(idx_N)*x(idx_Nj)/m(idx_Nj)
    conserveLinGetRef_x(idx_atom_C) = conserveLinGetRef_x(idx_atom_C) + m(idx_C)*x(idx_COj)/m(idx_COj)
    conserveLinGetRef_x(idx_atom_O) = conserveLinGetRef_x(idx_atom_O) + m(idx_O)*x(idx_COj)/m(idx_COj)
    conserveLinGetRef_x(idx_atom_F) = conserveLinGetRef_x(idx_atom_F) + m(idx_F)*x(idx_Fj)/m(idx_Fj)

  end function conserveLinGetRef_x

  !***************************
  !Ref: Sasaki & Takahara (1993)
  !This function evaluate the recombination rate
  ! for H+ + e --> H + gamma and the same
  ! for D+ + e --> D + gamma
  function elec_recomb_ST93(nabund,nelec,ntot,nucleiH,Trad)
    use krome_commons
    use krome_constants
    implicit none
    real*8::nabund,nelec,Trad
    real*8::nucleiH,elec_recomb_ST93
    real*8::al,ak,rc2,r2c
    real*8::a0,b0,c0,d0,e0
    real*8::a1,b1,c1,d1,e1,f1,g1,h1
    real*8::ntot,ratio

    al = 8.227d0
    ak = 22.06d0 / (hubble  *(1d0 + phys_zredshift) &
        * sqrt(1d0 + Omega0 * phys_zredshift))
    !Rc2 evaluation
    rc2 = 8.76d-11 * (1d0 + phys_zredshift)**(-0.58)
    !R2c evaluation
    r2c = (1.80d10 * Trad)**(1.5) &
        * exp(-3.9472d4 / Trad) * rc2

    !coefficients
    a0 = nucleiH * rc2
    b0 = ak * al * nucleiH
    c0 = ak * rc2 * nucleiH * nucleiH
    d0 = r2c * exp(-1.18416d5/Trad)
    e0 = ak * r2c * nucleiH

    !polynomial terms
    a1 = -d0 * (1d0 + b0)
    b1 = d0 * (1d0 + 2d0 * b0)
    c1 = a0 + b0 * (a0 - d0)
    d1 = -a0 * b0
    e1 = a0 * c0
    f1 = 1d0 + b0 + e0
    g1 = -(b0 + e0)
    h1 = c0

    ratio = nabund / ntot

    elec_recomb_ST93 = ntot*(a1 + b1*ratio + c1*ratio**2 + d1*ratio**3 &
        + e1*ratio**4) / (f1 + g1*ratio + h1*ratio**2)

    elec_recomb_ST93 = elec_recomb_ST93 / (nabund * nelec)

  end function elec_recomb_ST93

  !********************
  subroutine load_parts()
    use krome_commons
    implicit none

  end subroutine load_parts

  !*************************
  subroutine load_part(fname,array_part,min_part,dT_part)
    character(len=*)::fname
    integer::ios,icount,i,cv
    real*8,allocatable::array_part(:),emed(:)
    real*8::min_part,dT_part,Told,array_tmp(int(1e5)),rout(2)

    open(33,file=trim(fname),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: partition function not found"
      print *," in file "//fname
      stop
    end if

    print *,"loading partition function from "//fname
    icount = 0
    min_part = 1d99
    Told = 0d0
    do
      read(33,*,iostat=ios) rout(:)
      if(ios<0) exit
      if(ios.ne.0) cycle
      icount = icount + 1
      min_part = min(min_part,rout(1))
      array_tmp(icount) = rout(2)
      dT_part = rout(1) - Told
      Told = rout(1)
    end do
    close(33)

    allocate(array_part(icount),emed(icount))
    array_part(:) = array_tmp(1:icount)

  end subroutine load_part

  !**********************
  function troe_falloff(k0,kinf,Fc,m)
    implicit none
    real*8::troe_falloff,k0,kinf,Fc,m,rm,xexp
    rm = k0*m/kinf
    xexp = 1d0/(1d0+log10(rm)**2)
    troe_falloff = k0*m/(1d0+rm)*Fc**xexp
  end function troe_falloff

  !*************************
  function k3body(k0,kinf,Fc,nM)
    implicit none
    real*8::k3body,k0,kinf,Fc,nM
    real*8::c,n,d,Pr,xexp,F

    c = -0.4d0-0.67d0*log10(Fc)
    n = 0.75d0-1.27d0*log10(Fc)
    d = 0.14d0
    Pr = k0*nM/kinf
    xexp = (log10(Pr)+c)/(n-d*(log10(Pr)+c))
    F = 1d1**(log10(Fc)/(1d0+xexp**2))
    k3body = kinf*(Pr/(1d0+Pr)) * F

  end function k3body

  !***********************
  !see http://kida.obs.u-bordeaux1.fr/help
  function KIDA3body(ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,&
        kcFc,kdFc,npart,Tgas,pmin,pmax)
    implicit none
    real*8::ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,kcFc,kdFc
    real*8::KIDA3body,kinf,p,f,npart,Tgas,fc,fexp,invT
    real*8::k0,cc,dd,nn,pmin,pmax

    KIDA3body = 0d0

    invT = 1d0/Tgas
    k0 = ka0*(Tgas/3d2)**kb0*exp(-kc0*invT)
    kinf = kainf*(Tgas/3d2)**kbinf*exp(-kcinf*invT)

    p = k0*npart/kinf
    if(p<pmin) return
    if(p>pmax) return

    fc = (1d0-kaFc)*exp(-Tgas/kbFc) + kaFc*exp(-Tgas/kbFc) &
        + exp(-kdFc*invT)

    cc = -0.4d0 - 0.67d0 *log10(fc)
    dd = 0.14d0
    nn = 0.75d0 - 1.27d0*log10(fc)
    fexp = 1d0 + ((log10(p)+cc)/(nn-dd*(log10(p)+cc)))**2

    f = fc**(1d0/fexp)

    KIDA3body = kinf*(p/(1d0+p))*f

  end function KIDA3body

  !******************************
  !collisional ionization rate from Verner+96
  ! unit: cm3/s
  function colion_v96(Tgas,dE,P,A,X,K)
    implicit none
    real*8::colion_v96,Tgas,dE,A,X,K,U,Te,P

    Te = Tgas * 8.621738d-5 !K to eV
    U = dE / Te
    colion_v96 = A * (1d0 + P*sqrt(U)) * U**K * exp(-U) / (X+U)

  end function colion_v96

  !****************************
  !radiative recombination rates from
  ! Verner routine, standard fit, cm3/s
  function recV96(Tgas,a,b)
    implicit none
    real*8::recV96,Tgas,a,b

    recV96 = a*(1d4/Tgas)**b

  end function recV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, new fit, cm3/s
  function recNewV96(Tgas,r1,r2,r3,r4)
    implicit none
    real*8::recNewV96,Tgas,r1,r2,r3,r4,tt

    tt = sqrt(Tgas/r3)
    recNewV96 = r1/(tt*(tt + 1d0)**(1.-r2) &
        * (1d0 + sqrt(Tgas/r4))**(1.+r2))

  end function recNewV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, iron only, cm3/s
  function recFeV96(Tgas,r1,r2,r3)
    implicit none
    real*8::recFeV96,Tgas,r1,r2,r3,tt

    tt = sqrt(Tgas*1d-4)
    recFeV96 = r1/tt**(r2 + r3 + log10(tt))

  end function recFeV96

  !******************************
  !radiative recombination rates from Verner+96
  ! unit: cm3/s
  function radrec_v96(Tgas,a,b,T0,T1)
    implicit none
    real*8::Tgas,a,b,T0,T1,radrec_v96,iT0

    iT0 = 1d0/T0
    radrec_v96 = a/(sqrt(Tgas*iT0) + (1d0*sqrt(Tgas*iT0))**(1.-b) &
        * (1d0+sqrt(Tgas/T1))**(1+b))

  end function radrec_v96

  !*******************************
  !radiative recombination rates low-temp fit, Verner+96
  ! unit: cm3/s
  function radrec_low_v96(Tgas,a,b,c,d,f)
    implicit none
    real*8::Tgas,a,b,c,d,f,radrec_low_v96,t,invt

    t = Tgas*1d-4
    invt = 1d0/t

    radrec_low_v96 = 1d-12 * (a*invt + b + c*t + d*t**2) &
        * t**(-1.5) * exp(-f*invt)

    radrec_low_v96 = max(0d0,radrec_low_v96)

  end function radrec_low_v96

  !***************************
  !Collisional dissociation rate (cm-3/s) by Martin et al. 1996
  ! H2+H->H+H+H
  !NOTE: the use of this rate is suggested
  ! for high-density regime and in the presence of UV backgrounds.
  ! if necessary it must be included in the reaction file as
  ! H2,H,,H,H,H,,NONE,NONE,dissH2_Martin96(n,Tgas)
  function dissH2_Martin96(n,Tgas)
    use krome_commons
    use krome_getphys
    integer::i
    real*8::n(nspec),Tgas,dissH2_Martin96
    real*8::CDrates,logTv(4),k_CIDm(21,2),k_CID,invT,logT,n_c1,n_c2,n_H
    real*8::logk_h1,logk_h2,logk_l1,logk_l2,logn_c1,logn_c2,p,logk_CID
    real*8::logT2,logT3

    !k_CID = collision-induced dissociation + dissociative tunneling

    !Collisional dissociation of H2
    k_CIDm(:,1) = (/-178.4239d0, -68.42243d0, 43.20243d0, -4.633167d0, &
        69.70086d0, 40870.38d0, -23705.70d0, 128.8953d0, -53.91334d0, &
        5.315517d0, -19.73427d0, 16780.95d0, -25786.11d0, 14.82123d0, &
        -4.890915d0, 0.4749030d0, -133.8283d0, -1.164408d0, 0.8227443d0,&
        0.5864073d0, -2.056313d0/)

    !Dissociative tunneling of H2
    k_CIDm(:,2) = (/-142.7664d0, 42.70741d0, -2.027365d0, -0.2582097d0, &
        21.36094d0, 27535.31d0, -21467.79d0, 60.34928d0, -27.43096d0, &
        2.676150d0, -11.28215d0, 14254.55d0, -23125.20d0, 9.305564d0, &
        -2.464009d0, 0.1985955d0, 743.0600d0, -1.174242d0, 0.7502286d0, &
        0.2358848d0, 2.937507d0/)

    n_H  = get_Hnuclei(n(:))
    logT = log10(Tgas)
    invT = 1.0d0/Tgas
    logT2 = logT*logT
    logT3 = logT2*logT
    logTv = (/1.d0, logT, logT2, logT3/)
    k_CID = 0.d0
    do i=1,2
      logk_h1 = k_CIDm(1,i)*logTv(1) + k_CIDm(2,i)*logTv(2) + &
          k_CIDm(3,i)*logTv(3) + k_CIDm(4,i)*logTv(4) + &
          k_CIDm(5,i)*log10(1.d0+k_CIDm(6,i)*invT)
      logk_h2 = k_CIDm(7,i)*invT
      logk_l1 = k_CIDm(8,i)*logTv(1) + k_CIDm(9,i)*logTv(2) + &
          k_CIDm(10,i)*logTv(3) + k_CIDm(11,i)*log10(1.d0+k_CIDm(12,i)*invT)
      logk_l2 = k_CIDm(13,i)*invT
      logn_c1 = k_CIDm(14,i)*logTv(1) + k_CIDm(15,i)*logTv(2) &
          + k_CIDm(16,i)*logTv(3) + k_CIDm(17,i)*invT
      logn_c2 = k_CIDm(18,i) + logn_c1
      p = k_CIDm(19,i) + k_CIDm(20,i)*exp(-Tgas/1.850d3) &
          + k_CIDm(21,i)*exp(-Tgas/4.40d2)
      n_c1 = 1d1**(logn_c1)
      n_c2 = 1d1**(logn_c2)
      logk_CID = logk_h1 - (logk_h1 - logk_l1) / (1.d0 + (n_H/n_c1)**p) &
          + logk_h2 - (logk_h2 - logk_l2) / (1.d0 + (n_H/n_c2)**p)
      k_CID = k_CID + 1.d1**logk_CID
    enddo

    dissH2_Martin96 = k_CID

  end function dissH2_Martin96

  !***********************************
  subroutine init_exp_table()
    use krome_commons
    implicit none
    integer::i
    real*8::a

    do i=1,exp_table_na
      a = (i-1)*(exp_table_aMax-exp_table_aMin)/(exp_table_na-1) + exp_table_aMin
      exp_table(i) = exp(-a)
    end do

  end subroutine init_exp_table

  !***************************
  !get the index of the specie name
  function get_index(name)
    use krome_commons
    use krome_getphys
    integer::get_index,i
    character*16::names(nspec)
    character*(*)::name
    names(:) = get_names()
    get_index = -1 !default index
    !loop on species to found the specie named name
    do i=1,nspec
      !when found store and break loop
      if(trim(names(i))== trim(name)) then
        get_index = i !store index
        exit
      end if
    end do

    !error if species not found
    if(get_index<0) then
      print *,"ERROR: can't find the index of ",name
      stop
    end if

  end function get_index

  !*****************************
  !computes revers kinetics from reaction and
  ! product indexes
  function revKc(Tgas,ridx,pidx)
    implicit none
    real*8::revKc,Tgas
    integer::ridx(:),pidx(:),i

    revKc = 0.d0

    do i=1,size(pidx)
      revKc = revKc + revHS(Tgas,pidx(i))
    end do

    do i=1,size(ridx)
      revKc = revKc - revHS(Tgas,ridx(i))
    end do

  end function revKc

  !*****************************
  !compute H-S for species with index idx
  ! when temperature is Tgas
  function revHS(Tgas,idx)
    use krome_commons
    real*8::revHS,Tgas,Tgas2,Tgas3,Tgas4,invT,lnT,H,S
    real*8::p1(74,7), p2(74,7), Tlim(74,3), p(7)
    integer::idx

    p1(:,:) = 0.d0
    p2(:,:) = 0.d0
    Tlim(:,:) = 0.d0
    Tgas2 = Tgas * Tgas
    Tgas3 = Tgas2 * Tgas
    Tgas4 = Tgas3 * Tgas
    invT = 1d0/Tgas
    lnT = log(Tgas)

    p1(idx_Ok,:)  = (/2.90805921d0,&
        -0.00169804907d0,&
        2.98069955d-06,&
        -2.43835127d-09,&
        7.61229311d-13,&
        11435.7717d0,&
        2.80339097d0/)
    p1(idx_Hk,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        15976.167d0,&
        -1.1390139d0/)
    p1(idx_Ck,:)  = (/2.50025151d0,&
        -1.19774349d-06,&
        2.28919443d-09,&
        -1.98276803d-12,&
        6.44398056d-16,&
        70064.893d0,&
        4.87847086d0/)
    p1(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p1(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p1(idx_C,:)  = (/2.5542395d0,&
        -0.00032153772d0,&
        7.3379223d-07,&
        -7.3223487d-10,&
        2.6652144d-13,&
        85442.681d0,&
        4.5313085d0/)
    p1(idx_NO,:)  = (/4.21859896d0,&
        -0.00463988124d0,&
        1.10443049d-05,&
        -9.34055507d-09,&
        2.80554874d-12,&
        9845.09964d0,&
        2.28061001d0/)
    p1(idx_CO,:)  = (/3.5795335d0,&
        -0.00061035369d0,&
        1.0168143d-06,&
        9.0700586d-10,&
        -9.0442449d-13,&
        -14344.086d0,&
        3.5084093d0/)
    p1(idx_N,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        56104.638d0,&
        4.1939088d0/)
    p1(idx_O2,:)  = (/3.78245636d0,&
        -0.00299673416d0,&
        9.84730201d-06,&
        -9.68129509d-09,&
        3.24372837d-12,&
        -1063.94356d0,&
        3.65767573d0/)
    p1(idx_O,:)  = (/3.1682671d0,&
        -0.00327931884d0,&
        6.64306396d-06,&
        -6.12806624d-09,&
        2.11265971d-12,&
        29122.2592d0,&
        2.05193346d0/)
    p1(idx_H2,:)  = (/2.34433112d0,&
        0.00798052075d0,&
        -1.9478151d-05,&
        2.01572094d-08,&
        -7.37611761d-12,&
        -917.935173d0,&
        0.683010238d0/)
    p1(idx_OH,:)  = (/3.99198424d0,&
        -0.00240106655d0,&
        4.61664033d-06,&
        -3.87916306d-09,&
        1.36319502d-12,&
        3368.89836d0,&
        -0.103998477d0/)
    p1(idx_CH,:)  = (/3.4897583d0,&
        0.0003243216d0,&
        -1.6899751d-06,&
        3.162842d-09,&
        -1.4061803d-12,&
        70660.755d0,&
        2.0842841d0/)
    p1(idx_N2,:)  = (/3.53100528d0,&
        -0.000123660988d0,&
        -5.02999433d-07,&
        2.43530612d-09,&
        -1.40881235d-12,&
        -1046.97628d0,&
        2.96747038d0/)
    p1(idx_CH2,:)  = (/3.84261832d0,&
        -7.36676871d-06,&
        6.16970693d-06,&
        -6.96689962d-09,&
        2.64620979d-12,&
        45863.1528d0,&
        1.2758447d0/)
    p1(idx_CO2,:)  = (/2.356813d0,&
        0.0089841299d0,&
        -7.1220632d-06,&
        2.4573008d-09,&
        -1.4288548d-13,&
        -48371.971d0,&
        9.9009035d0/)
    p1(idx_HCO,:)  = (/4.36380907d0,&
        -0.00535204137d0,&
        2.31954508d-05,&
        -2.6610904d-08,&
        1.02711962d-11,&
        25010.8717d0,&
        2.98106307d0/)
    p1(idx_H2O,:)  = (/4.1986352d0,&
        -0.0020364017d0,&
        6.5203416d-06,&
        -5.4879269d-09,&
        1.771968d-12,&
        -30293.726d0,&
        -0.84900901d0/)
    p1(idx_HEj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        285323.374d0,&
        1.62166556d0/)
    p1(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p1(idx_OHj,:)  = (/3.50502572d0,&
        0.000241313747d0,&
        -1.42200948d-06,&
        2.64780232d-09,&
        -1.17038711d-12,&
        155210.676d0,&
        1.97907627d0/)
    p1(idx_H2j,:)  = (/3.77256072d0,&
        -0.0019574659d0,&
        4.54812047d-06,&
        -2.82152141d-09,&
        5.33969209d-13,&
        178694.654d0,&
        -3.96609192d0/)
    p1(idx_FEj,:)  = (/2.76418106d0,&
        0.00286948238d0,&
        -7.61235651d-06,&
        8.18183334d-09,&
        -3.11792199d-12,&
        141159.039d0,&
        5.53997981d0/)
    p1(idx_Oj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        187935.284d0,&
        4.39337676d0/)
    p1(idx_Cj,:)  = (/2.61332254d0,&
        -0.000540148065d0,&
        1.03037233d-06,&
        -8.90092552d-10,&
        2.88500586d-13,&
        216862.274d0,&
        3.8345479d0/)
    p1(idx_COj,:)  = (/3.77061642d0,&
        -0.00201773246d0,&
        4.61081738d-06,&
        -2.99175463d-09,&
        6.06065045d-13,&
        149006.795d0,&
        3.38129783d0/)
    p2(idx_Ok,:)  = (/2.54474869d0,&
        -4.66695513d-05,&
        1.84912357d-08,&
        -3.18159223d-12,&
        1.98962956d-16,&
        11504.2089d0,&
        4.52131015d0/)
    p2(idx_Hk,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        15976.167d0,&
        -1.1390139d0/)
    p2(idx_Ck,:)  = (/2.50001597d0,&
        -1.71721376d-08,&
        6.9283294d-12,&
        -1.20607892d-15,&
        7.60308635d-20,&
        70064.9324d0,&
        4.87955907d0/)
    p2(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p2(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p2(idx_C,:)  = (/2.605583d0,&
        -0.00019593434d0,&
        1.0673722d-07,&
        -1.642394d-11,&
        8.187058d-16,&
        85411.742d0,&
        4.1923868d0/)
    p2(idx_NO,:)  = (/3.26071234d0,&
        0.00119101135d0,&
        -4.29122646d-07,&
        6.94481463d-11,&
        -4.03295681d-15,&
        9921.43132d0,&
        6.36900518d0/)
    p2(idx_CO,:)  = (/3.0484859d0,&
        0.0013517281d0,&
        -4.8579405d-07,&
        7.8853644d-11,&
        -4.6980746d-15,&
        -14266.117d0,&
        6.0170977d0/)
    p2(idx_N,:)  = (/2.4159429d0,&
        0.00017489065d0,&
        -1.1902369d-07,&
        3.0226244d-11,&
        -2.0360983d-15,&
        56133.775d0,&
        4.6496095d0/)
    p2(idx_O2,:)  = (/3.66096065d0,&
        0.000656365811d0,&
        -1.41149627d-07,&
        2.05797935d-11,&
        -1.29913436d-15,&
        -1215.97718d0,&
        3.41536279d0/)
    p2(idx_O,:)  = (/2.54363697d0,&
        -2.73162486d-05,&
        -4.1902952d-09,&
        4.95481845d-12,&
        -4.79553694d-16,&
        29226.012d0,&
        4.92229457d0/)
    p2(idx_H2,:)  = (/2.93286575d0,&
        0.000826608026d0,&
        -1.46402364d-07,&
        1.54100414d-11,&
        -6.888048d-16,&
        -813.065581d0,&
        -1.02432865d0/)
    p2(idx_OH,:)  = (/2.83853033d0,&
        0.00110741289d0,&
        -2.94000209d-07,&
        4.20698729d-11,&
        -2.4228989d-15,&
        3697.80808d0,&
        5.84494652d0/)
    p2(idx_CH,:)  = (/2.5209369d0,&
        0.0017653639d0,&
        -4.614766d-07,&
        5.9289675d-11,&
        -3.3474501d-15,&
        70994.878d0,&
        7.4051829d0/)
    p2(idx_N2,:)  = (/2.95257637d0,&
        0.0013969004d0,&
        -4.92631603d-07,&
        7.86010195d-11,&
        -4.60755204d-15,&
        -923.948688d0,&
        5.87188762d0/)
    p2(idx_CH2,:)  = (/3.11049513d0,&
        0.00373779517d0,&
        -1.37371977d-06,&
        2.23054839d-10,&
        -1.33567178d-14,&
        45971.5953d0,&
        4.62796405d0/)
    p2(idx_CO2,:)  = (/4.6365111d0,&
        0.0027414569d0,&
        -9.9589759d-07,&
        1.6038666d-10,&
        -9.1619857d-15,&
        -49024.904d0,&
        -1.9348955d0/)
    p2(idx_HCO,:)  = (/4.23892214d0,&
        0.0019657617d0,&
        -3.82075171d-07,&
        4.80137647d-11,&
        -3.11176347d-15,&
        24726.1645d0,&
        1.99698242d0/)
    p2(idx_H2O,:)  = (/2.6770389d0,&
        0.0029731816d0,&
        -7.7376889d-07,&
        9.4433514d-11,&
        -4.2689991d-15,&
        -29885.894d0,&
        6.88255d0/)
    p2(idx_HEj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        285323.374d0,&
        1.62166556d0/)
    p2(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p2(idx_OHj,:)  = (/2.68358996d0,&
        0.00157006435d0,&
        -5.39972815d-07,&
        9.37643877d-11,&
        -5.70068067d-15,&
        155479.296d0,&
        6.44375894d0/)
    p2(idx_H2j,:)  = (/3.44204765d0,&
        0.000599083239d0,&
        6.69133685d-08,&
        -3.43574373d-11,&
        1.97626599d-15,&
        178650.236d0,&
        -2.79499055d0/)
    p2(idx_FEj,:)  = (/3.33602399d0,&
        -0.000272549262d0,&
        8.05440344d-09,&
        1.51229089d-11,&
        -1.43376595d-15,&
        141036.455d0,&
        2.86476968d0/)
    p2(idx_Oj,:)  = (/2.48542028d0,&
        2.56978695d-05,&
        -1.28833378d-08,&
        1.65525487d-12,&
        1.09933344d-16,&
        187940.874d0,&
        4.47425446d0/)
    p2(idx_Cj,:)  = (/2.50827618d0,&
        -1.04354146d-05,&
        5.16160809d-09,&
        -1.14187475d-12,&
        9.43539946d-17,&
        216879.645d0,&
        4.3188599d0/)
    p2(idx_COj,:)  = (/2.93062935d0,&
        0.00156033262d0,&
        -6.16246355d-07,&
        1.09957336d-10,&
        -6.66119284d-15,&
        149147.222d0,&
        7.3384673d0/)
    Tlim(idx_Ok,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_Hk,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_Ck,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_H,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_HE,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_C,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_NO,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_CO,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_N,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_O2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_O,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_H2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_OH,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_CH,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_N2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_CH2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_CO2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_HCO,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_H2O,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_HEj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_Hj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_OHj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_H2j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_FEj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_Oj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_Cj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim(idx_COj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)

    if(Tlim(idx,2)==0.d0) then
      revHS = 0.d0
      return
    end if

    !select set of NASA polynomials using temperature
    if(Tlim(idx,1).le.Tgas .and. Tgas.le.Tlim(idx,2)) p(:) = p1(idx,:)
    if(Tlim(idx,2)<Tgas .and. Tgas.le.Tlim(idx,3)) p(:) = p2(idx,:)

    !compute NASA polynomials for enthalpy and enthropy
    H = p(1) + p(2)*0.5d0*Tgas + p(3)*Tgas2/3.d0 + p(4)*Tgas3*0.25d0 + &
        p(5)*Tgas4*0.2d0 + p(6)*invT
    S = p(1)*lnT + p(2)*Tgas + p(3)*Tgas2*0.5d0 + p(4)*Tgas3/3.d0 + &
        p(5)*Tgas4*0.25d0 + p(7)

    revHS = H - S

  end function revHS

  !******************************
  subroutine print_best_flux(n,Tgas,nbestin)
    !print the first nbestin fluxes
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea)
    integer::nbest,idx(nrea),i,nbestin
    character*50::name(nrea)

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux

  !******************************
  subroutine print_best_flux_frac(n,Tgas,frac)
    !print the first nbestin fluxes
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),frac
    integer::idx(nrea),i
    character*50::name(nrea)

    if(frac>1d0) then
      print *,"ERROR: fraction in krome_print_best_flux should be <=1!"
      stop
    end if

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nrea
      if(flux(idx(i))<flux(idx(1))*frac) exit
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux_frac

  !******************************
  subroutine print_best_flux_spec(n,Tgas,nbestin,idx_found)
    !print the first nbestin fluxes for the reactions
    ! that contains the species with index idx_found
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),maxflux
    integer::nbest,idx(nrea),i,nbestin,idx_found
    character*50::name(nrea)
    logical::found

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions
    maxflux = 0d0
    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names
    do i=1,nrea
      found = .false.
      if(arr_r1(i) == idx_found) found = .true.
      if(arr_r2(i) == idx_found) found = .true.
      if(arr_r3(i) == idx_found) found = .true.
      if(arr_p1(i) == idx_found) found = .true.
      if(arr_p2(i) == idx_found) found = .true.
      if(arr_p3(i) == idx_found) found = .true.
      maxflux = max(maxflux,flux(i))
      if(.not.found) flux(i) = 0d0
    end do

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,2E17.8)',idx(i)," ",name(idx(i)),flux(idx(i)),&
          flux(idx(i))/maxflux
    end do

  end subroutine print_best_flux_spec

  !*****************************
  function idx_sort(fin)
    !sorting algorithm: requires an array of real values fin
    ! and returns the sorted index list. descending.
    ! bubblesort: not very efficient, replace with what you prefer
    implicit none
    real*8::fin(:),f(size(fin)),ftmp
    integer::idx_sort(size(fin)),n,itmp,i
    logical::found

    f(:) = fin(:) !copy to local

    n = size(f)
    !init indexes
    do i=1,n
      idx_sort(i) = i
    end do

    !loop to sort
    do
      found = .false. !swapped something flag
      do i=2,n
        !> for descending, < for ascending
        if(f(i)>f(i-1)) then
          found = .true.
          !swap real value
          ftmp = f(i)
          f(i) = f(i-1)
          f(i-1) = ftmp
          !swap index
          itmp = idx_sort(i)
          idx_sort(i) = idx_sort(i-1)
          idx_sort(i-1) = itmp
        end if
      end do
      !if nothing swapped exit
      if(.not.found) exit
    end do

  end function idx_sort

  !******************************
  function get_flux(n,Tgas)
    !get the flux k*n*n*... of the rates
    use krome_commons
    implicit none
    integer::i
    integer::r1,r2,r3
    real*8::get_flux(nrea),n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      r3 = arr_r3(i)
      arr_flux(i) = k(i)*n(r1)*n(r2)*n(r3)
    end do
    get_flux(:) = arr_flux(:)

  end function get_flux

  !*****************************
  subroutine load_arrays()
    !load the array containing reactants
    ! and product index
    use krome_commons

    arr_r1(1:255) = (/6,8,50,48,12,6,49,49,11,48,6,16,15,52,6,6&
        ,11,6,26,18,6,3,29,3,55,16,18,6,49,6,15,56,6,8,11,6,16,6,15,6&
        ,34,16,6,18,34,8,8,11,17,6,18,8,6,16,16,6,6,56,41,11,49,6,2&
        ,63,3,3,6,62,3,6,17,15,3,62,16,6,18,26,64,8,18,15,6,49,24,15&
        ,54,6,17,8,16,6,6,11,6,18,15,29,29,28,8,56,6,16,8,8,16,6,66,6&
        ,67,8,24,34,64,15,16,15,17,15,67,15,17,67,13,41,15,11,3,41,29&
        ,26,5,17,6,53,13,60,65,49,15,57,17,49,3,6,18,18,16,6,11,8,6,6&
        ,11,11,8,34,16,15,18,8,6,41,29,8,49,15,8,6,8,49,8,18,11,6,31&
        ,16,11,15,29,6,11,49,6,15,52,19,8,59,13,3,8,6,26,6,53,8,2,29&
        ,58,8,3,8,18,56,8,6,61,2,11,15,15,3,16,8,6,15,18,6,7,16,16,3&
        ,3,3,6,16,6,49,16,8,6,11,10,16,16,7,15,16,32,10,67,68,69,48&
        ,49,51,60,48,16,11,11,42,11/)
    arr_r2(1:255) = (/48,9,1,1&
        ,13,51,17,18,19,17,55,19,1,22,53,20,23,25,52,27,14,59,13,57,1&
        ,12,13,12,1,9,30,7,16,23,9,10,27,33,15,34,15,29,35,10,15,10&
        ,24,34,12,2,37,60,39,8,62,40,1,15,52,12,6,19,63,1,49,52,18,22&
        ,11,9,60,24,62,1,8,4,23,57,22,14,53,23,15,42,13,19,1,8,10,19&
        ,13,29,12,14,43,24,64,15,13,53,9,1,44,34,34,32,1,29,1,58,41&
        ,15,12,13,1,35,11,35,52,18,17,14,9,22,30,64,35,24,8,53,11,53&
        ,49,36,45,22,1,22,1,37,19,1,65,41,53,44,18,14,18,40,31,29,31&
        ,40,18,64,1,13,15,23,21,34,21,60,11,8,13,31,18,5,19,22,11,23&
        ,36,36,13,52,19,32,15,62,18,15,14,17,1,19,30,1,19,6,14,61,63&
        ,56,1,20,49,15,1,64,15,18,17,8,13,18,1,57,14,20,44,63,24,12&
        ,20,24,24,1,1,49,1,1,6,49,6,6,6,26,74,74,74,74,74,74,74,74,74&
        ,74,74,74,1,1,1,38,34,1,1,28,70,46,47,12&
        ,47/)
    arr_r3(1:255) = (/74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,6,6,7,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74/)
    arr_p1(1:255) = (/7,10,11,7,14,52,53,54,20,53&
        ,7,21,2,57,58,19,13,19,13,28,13,6,23,6,7,18,14,15,6,18,14,55&
        ,6,13,32,18,28,29,12,11,18,33,24,36,9,31,31,32,37,18,38,10,34&
        ,29,54,18,3,54,13,9,56,13,15,41,6,6,15,57,34,15,58,9,6,15,33&
        ,19,10,22,57,23,64,14,18,65,20,13,15,29,37,23,19,8,18,20,10&
        ,35,12,18,19,66,24,6,18,39,11,24,6,8,17,53,63,10,40,20,17,10&
        ,34,24,13,12,53,13,37,57,14,37,40,32,29,17,34,17,8,37,19,57,4&
        ,57,42,64,14,22,42,63,6,18,44,30,44,35,24,31,29,34,15,9,5,19&
        ,18,13,19,24,19,43,24,31,52,10,15,29,13,57,24,6,9,10,23,51,13&
        ,9,10,15,9,62,19,37,13,21,10,26,25,16,13,67,41,16,17,13,15,60&
        ,17,53,18,10,37,61,23,15,8,15,13,13,18,6,35,10,13,10,40,49,48&
        ,56,6,6,6,56,16,16,16,59,49,67,49,68,69,49,6,48,62,56,11,8,8&
        ,11,15,12,50,13,10,70,56,42,46,47&
        ,42/)
    arr_p2(1:255) = (/49,11,6,74,15,16,6,6,6,7,56,6,74&
        ,13,74,11,24,13,59,15,18,26,6,22,6,18,6,15,74,11,12,6,6,31,15&
        ,8,6,16,74,16,11,6,16,6,6,15,11,6,15,1,6,61,16,6,6,24,74,6,63&
        ,15,74,16,41,74,6,13,6,15,1,34,10,8,15,74,74,1,19,59,37,15,6&
        ,8,74,6,8,18,6,74,8,6,6,16,15,15,16,15,53,8,8,6,15,6,6,6,29&
        ,11,6,6,27,16,8,74,15,6,15,34,6,18,53,6,8,12,11,8,14,63,6,8,1&
        ,63,8,59,6,10,10,17,74,43,74,6,6,74,53,6,17,16,15,6,6,15,8,6&
        ,8,10,34,53,74,11,6,10,44,6,16,63,6,74,6,8,29,1,29,6,74,45,10&
        ,18,8,6,34,11,6,49,6,6,15,74,74,13,14,74,6,1,10,16,59,49,74&
        ,24,6,1,6,10,1,6,6,6,74,16,6,22,9,9,18,41,6,15,34,11,6,1,1,6&
        ,3,1,6,1,6,16,7,6,3,1,1,1,1,6,6,1,1,1,11,15,74,74,8,53,6,6,6&
        ,6,27,32,15,15,9/)
    arr_p3(1:255) = (/74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,6,74,74,74,74,6,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,6,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,6,74,74,74,1,6,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74,74&
        ,74,74,74,74,74,74,74,74,74,74,1,1,74,74,1,1,74,74,74,74,74&
        ,74,74,74,74,74,1,74,74,74,74,74,74,74,74,74,7,74,74,74,7,74&
        ,74,74,74,74/)

  end subroutine load_arrays

  !*********************************
  subroutine mydgesv(n,Ain,Bin, parent_name)
    !driver for LAPACK dgesv
    integer::n,info,i,ipiv(n)
    real*8,allocatable::tmp(:)
    real*8::A(n,n),B(n),Ain(:,:),Bin(:),suml,sumr,tmpn(n)
    character(len=*)::parent_name
    A(:,:) = Ain(1:n,1:n)
    B(:) = Bin(1:n)
    call dgesv(n,1,A,n,ipiv,B,n,info)
    Bin(1:n) = B(:)

    !write some info about the error and stop
    if(info > 0) then
      allocate(tmp(size(Bin)))
      print *,"ERROR: matrix exactly singular, U(i,i) where i=",info
      print *,' (called by "'//trim(parent_name)//'" function)'

      !dump the input matrix to a file
      open(97,file="ERROR_dump_dgesv.dat",status="replace")
      !dump size of the problem
      write(97,*) "size of the problem:",n
      write(97,*)

      !dump matrix A
      write(97,*) "Input matrix A line by line:"
      do i=1,size(Ain,1)
        tmp(:) = Ain(i,:)
        write(97,'(I5,999E17.8e3)') i,tmp(:)
      end do

      !dump matrix A
      write(97,*)
      write(97,*) "Workin matrix A line by line:"
      do i=1,n
        tmpn(:) = Ain(i,1:n)
        write(97,'(I5,999E17.8e3)') i,tmpn(:)
      end do

      !dump matrix B
      write(97,*)
      write(97,*) "Input/output vector B element by element"
      do i=1,n
        write(97,*) i, Bin(i),B(i)
      end do

      !dump info on matrix A rows
      write(97,*)
      write(97,*) "Info on matrix A rows"
      write(97,'(a5,99a17)') "idx","minval","maxval"
      do i=1,size(Ain,1)
        write(97,'(I5,999E17.8e3)') i, minval(Ain(i,:)), &
            maxval(Ain(i,:))
      end do

      !dump info on matrix sum left and right
      write(97,*)
      write(97,*) "Info on matrix A, sum left/right"
      write(97,'(a5,99a17)') "idx","left","right"
      suml = 0d0
      sumr = 0d0
      do i=1,size(Ain,1)
        if(i>1) suml = sum(Ain(i,:i-1))
        if(i<n) sumr = sum(Ain(i,i+1:))
        write(97,'(I5,999E17.8e3)') i, suml, sumr
      end do
      close(97)

      print *,"Input A and B dumped in ERROR_dump_dgesv.dat"

      stop
    end if

    !if error print some info and stop
    if(info<0) then
      print *,"ERROR: input error position ",info
      print *,' (called by "'//trim(parent_name)//'" function)'
      stop
    end if

  end subroutine mydgesv

end module krome_subs

!############### MODULE ##############
module krome_stars

end module krome_stars

!############### MODULE ##############
module krome_dust
contains

  !***********************
  subroutine init_dust_tabs()
    use krome_commons
    use krome_fit
    implicit none

    call init_anytab2D("dust_table_cool.dat",dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_cool(:,:), dust_mult_ngas, &
        dust_mult_Tgas)
    call init_anytab2D("dust_table_Tdust.dat",dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_Tdust(:,:), dust_mult_ngas, &
        dust_mult_Tgas)
    call init_anytab2D("dust_table_H2.dat",dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_H2(:,:), dust_mult_ngas, &
        dust_mult_Tgas)

  end subroutine init_dust_tabs

end module krome_dust

!############### MODULE ##############
module krome_photo
contains

end module krome_photo

!############### MODULE ##############
module krome_tabs
contains

  !***********************+
  function coe_tab(n)
    !interface to tabs
    use krome_subs
    use krome_getphys
    use krome_phfuncs
    use krome_grfuncs
    use krome_constants
    use krome_commons
    use krome_user_commons
    implicit none
    integer::idx,j
    real*8::Tgas, coe_tab(nrea),n(nspec),small

    Tgas = max(n(idx_Tgas),phys_Tcmb)
    small = 0d0

    coe_tab(:) = coe(n(:))

  end function coe_tab

end module krome_tabs

!############### MODULE ##############
module KROME_coolingGH
end module KROME_coolingGH


!############### MODULE ##############
module KROME_cooling
  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2018-01-25 18:47:28
  !  Changeset 9a6e335
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************
  integer,parameter::coolTab_n=int(1e2)
  integer,parameter::nZrate=110
  real*8::coolTab(nZrate,coolTab_n),coolTab_logTlow, coolTab_logTup
  real*8::coolTab_T(coolTab_n),inv_coolTab_T(coolTab_n-1),inv_coolTab_idx
  real*8::pop_level_CI(3)
  real*8::pop_level_OI(3)
  real*8::pop_level_OII(3)
  real*8::pop_level_FeII(5)
  real*8::pop_level_CII(2)
  real*8::pop_level_FeI(5)
  real*8::pop_level_SiII(2)
  real*8::pop_level_SiI(3)
  !$omp threadprivate(pop_level_CI)
  !$omp threadprivate(pop_level_OI)
  !$omp threadprivate(pop_level_OII)
  !$omp threadprivate(pop_level_FeII)
  !$omp threadprivate(pop_level_CII)
  !$omp threadprivate(pop_level_FeI)
  !$omp threadprivate(pop_level_SiII)
  !$omp threadprivate(pop_level_SiI)
contains

  !*******************
  function cooling(n,inTgas)
    use krome_commons
    implicit none
    real*8::n(:),inTgas,cooling,Tgas

    Tgas = inTgas
    cooling = sum(get_cooling_array(n(:),Tgas))

  end function cooling

  !*******************************
  function get_cooling_array(n, Tgas)
    use krome_commons
    implicit none
    real*8::n(:), Tgas
    real*8::get_cooling_array(ncools),cools(ncools)
    real*8::f1,f2,smooth

    f1 = 1d0
    f2 = 1d0

    !returns cooling in erg/cm3/s
    cools(:) = 0d0

    cools(idx_cool_H2) = cooling_H2(n(:), Tgas)

    cools(idx_cool_CO) = cooling_CO(n(:), Tgas)

    cools(idx_cool_atomic) = f2 * ( cooling_Atomic(n(:), Tgas)  )

    cools(idx_cool_Z) = f2 * ( cooling_Z(n(:), Tgas)  )

    cools(idx_cool_CIE) = f2 * cooling_CIE(n(:), Tgas)

    cools(idx_cool_custom) = cooling_custom(n(:),Tgas)

    get_cooling_array(:) = cools(:)

  end function get_cooling_array

  !***************************
  !CO cooling: courtesy of K.Omukai (Nov2014)
  ! method: Neufeld+Kaufman 1993 (bit.ly/1vnjcXV, see eqn.5).
  ! see also Omukai+2010 (bit.ly/1HIaGcn)
  ! H and H2 collisions
  function cooling_CO(n,inTgas)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    integer,parameter::imax=coolCOn1
    integer,parameter::jmax=coolCOn2
    integer,parameter::kmax=coolCOn3
    integer::i,j,k
    real*8,parameter::eps=1d-5
    real*8::cooling_CO,n(:),inTgas
    real*8::v1,v2,v3,prev1,prev2,cH
    real*8::vv1,vv2,vv3,vv4,vv12,vv34,xLd
    real*8::x1(imax),x2(jmax),x3(kmax)
    real*8::ixd1(imax-1),ixd2(jmax-1),ixd3(kmax-1)
    real*8::v1min,v1max,v2min,v2max,v3min,v3max

    !local copy of limits
    v1min = coolCOx1min
    v1max = coolCOx1max
    v2min = coolCOx2min
    v2max = coolCOx2max
    v3min = coolCOx3min
    v3max = coolCOx3max

    !local copy of variables arrays
    x1(:) = coolCOx1(:)
    x2(:) = coolCOx2(:)
    x3(:) = coolCOx3(:)

    ixd1(:) = coolCOixd1(:)
    ixd2(:) = coolCOixd2(:)
    ixd3(:) = coolCOixd3(:)

    !local variables
    v3 = num2col(n(idx_CO),n(:)) !CO column density
    cH = n(idx_H) + n(idx_H2)
    if (n(idx_H) < 0. .or. n(idx_H2) < 0.) then
      cH = n_global(idx_H) + n_global(idx_H2)
      ! Set red flag if not due to small excursion in n_H2
      ! Only relevant for low temperatures, bc otherwise H is ionised
      ! And CO cooling only relevant at high densities
      if (abs(n(idx_H2)) / (abs(n(idx_H)) + 1d-40) > 1d-6 .and. &
          inTgas < 5d4 .and. &
          abs(n(idx_H)) > abs(n(idx_Hj)) .and. sum(n(1:nmols)) > 100.) &
          red_flag = ibset(red_flag,5)
    endif

    v2 = cH
    v1 = inTgas !Tgas

    !logs of variables
    v1 = log10(v1)
    v2 = log10(v2)
    v3 = log10(v3)

    !default value erg/s/cm3
    cooling_CO = 0d0

    !check limits
    if(v1>=v1max) v1 = v1max*(1d0-eps)
    if(v2>=v2max) v2 = v2max*(1d0-eps)
    if(v3>=v3max) v3 = v3max*(1d0-eps)

    if(v1<v1min) return
    if(v2<v2min) return
    if(v3<v3min) return

    !gets position of variable in the array
    i = (v1-v1min)*coolCOdvn1+1
    j = (v2-v2min)*coolCOdvn2+1
    k = (v3-v3min)*coolCOdvn3+1

    !precompute shared variables
    prev1 = (v1-x1(i))*ixd1(i)
    prev2 = (v2-x2(j))*ixd2(j)

    !linear interpolation on x1 for x2,x3
    vv1 = prev1 * (coolCOy(k,j,i+1) - &
        coolCOy(k,j,i)) + coolCOy(k,j,i)
    !linear interpolation on x1 for x2+dx2,x3
    vv2 = prev1 * (coolCOy(k,j+1,i+1) - &
        coolCOy(k,j+1,i)) + coolCOy(k,j+1,i)
    !linear interpolation on x2 for x3
    vv12 = prev2 * (vv2 - vv1) + vv1

    !linear interpolation on x1 for x2,x3+dx3
    vv3 = prev1 * (coolCOy(k+1,j,i+1) - &
        coolCOy(k+1,j,i)) + coolCOy(k+1,j,i)
    !linear interpolation on x1 for x2+dx2,x3+dx3
    vv4 = prev1 * (coolCOy(k+1,j+1,i+1) - &
        coolCOy(k+1,j+1,i)) + coolCOy(k+1,j+1,i)
    !linear interpolation on x2 for x3+dx3
    vv34 = prev2 * (vv4 - vv3) + vv3

    !linear interpolation on x3
    xLd = (v3-x3(k))*ixd3(k)*(vv34 - &
        vv12) + vv12

    !CO cooling in erg/s/cm3
    cooling_CO = 1d1**xLd * cH * n(idx_CO)

  end function cooling_CO

  !************************
  subroutine init_coolingCO()
    use krome_commons
    implicit none
    integer::ios,iout(3),i
    real*8::rout(4)

    if (krome_mpi_rank<=1) print *,"load CO cooling..."
    open(33,file="coolCO.dat",status="old",iostat=ios)
    !check if file exists
    if(ios.ne.0) then
      print *,"ERROR: problems loading coolCO.dat!"
      stop
    end if

    do
      read(33,*,iostat=ios) iout(:),rout(:) !read line
      if(ios<0) exit !eof
      if(ios/=0) cycle !skip blanks
      coolCOx1(iout(1)) = rout(1)
      coolCOx2(iout(2)) = rout(2)
      coolCOx3(iout(3)) = rout(3)
      coolCOy(iout(3),iout(2),iout(1)) = rout(4)
    end do

    !store inverse of the differences
    ! to speed up interpolation
    do i=1,coolCOn1-1
      coolCOixd1(i) = 1d0/(coolCOx1(i+1)-coolCOx1(i))
    end do
    do i=1,coolCOn2-1
      coolCOixd2(i) = 1d0/(coolCOx2(i+1)-coolCOx2(i))
    end do
    do i=1,coolCOn3-1
      coolCOixd3(i) = 1d0/(coolCOx3(i+1)-coolCOx3(i))
    end do

    coolCOx1min = minval(coolCOx1)
    coolCOx1max = maxval(coolCOx1)
    coolCOx2min = minval(coolCOx2)
    coolCOx2max = maxval(coolCOx2)
    coolCOx3min = minval(coolCOx3)
    coolCOx3max = maxval(coolCOx3)

    coolCOdvn1 = (coolCOn1-1)/(coolCOx1max-coolCOx1min)
    coolCOdvn2 = (coolCOn2-1)/(coolCOx2max-coolCOx2min)
    coolCOdvn3 = (coolCOn3-1)/(coolCOx3max-coolCOx3min)

  end subroutine init_coolingCO

  !*****************************
  function cooling_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,cooling_custom

    cooling_custom = 0d0

  end function cooling_custom

  !**********************************
  function kpla(n,Tgas)
    !Planck opacity mean fit (Lenzuni+1996)
    !only denisity dependent (note that the
    ! fit provided by Lenzuni is wrong)
    ! valid for T<3e3 K
    !use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::kpla,rhogas,Tgas,n(:),y
    real*8::a0,a1,m(nspec)

    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3

    kpla = 0.d0
    !opacity is zero under 1e-12 g/cm3
    if(rhogas<1d-12) return

    !fit coefficients
    a0 = 1.000042d0
    a1 = 2.14989d0

    !log density cannot exceed 0.5 g/cm3
    y = log10(min(rhogas,0.5d0))

    kpla = 1d1**(a0*y + a1) !fit density only

  end function kpla

  !*****************************
  function coolingChem(n,Tgas)
    implicit none
    real*8::coolingChem,n(:),Tgas

    !note that this function is a dummy.
    ! For chemical cooling you should see
    ! heatingChem function in krome_heating.f90

    coolingChem = 0.d0

  end function coolingChem

  !*******************************
  function cooling_CIE(n, inTgas)
    !CIE cooling: fit from Ripamponti&Abel2004 (RA04) data
    ! The fit is valid from 100K-1e6K.
    ! Original data are from 400K to 7000K.
    ! We extrapolated data under 400K and fitted from 100K to 10**2.95 K.
    ! Data from 10**2.95 K to 1e5K are fitted analogously.
    ! Above 1e5 we employ a cubic extrapolation.
    use krome_commons
    use krome_constants
    real*8::cooling_CIE,n(:),Tgas,inTgas
    real*8::x,x2,x3,x4,x5
    real*8::a0,a1,a2,a3,a4,a5
    real*8::b0,b1,b2,b3,b4,b5
    real*8::cool,tauCIE,logcool

    cooling_CIE = 0d0
    !under 1e12 1/cm3 cooling is zero
    if(n(idx_H2)<1d12) return

    Tgas = inTgas
    !temperature limit
    if(Tgas<phys_Tcmb) return

    !prepares variables
    x = log10(Tgas)
    x2 = x*x
    x3 = x2*x
    x4 = x3*x
    x5 = x4*x

    cool = 0.d0
    !outside boundaries below cooling is zero
    logcool = -1d99

    !evaluates fitting functions
    if(x>2.d0 .and. x<2.95d0) then
      a0 = -30.3314216559651d0
      a1 = 19.0004016698518d0
      a2 = -17.1507937874082d0
      a3 = 9.49499574218739d0
      a4 = -2.54768404538229d0
      a5 = 0.265382965410969d0
      logcool = a0 + a1*x + a2*x2 + a3*x3 +a4*x4 +a5*x5
    elseif(x.GE.2.95d0 .and. x<5.d0) then
      b0 = -180.992524120965d0
      b1 = 168.471004362887d0
      b2 = -67.499549702687d0
      b3 = 13.5075841245848d0
      b4 = -1.31983368963974d0
      b5 = 0.0500087685129987d0
      logcool = b0 + b1*x + b2*x2 + b3*x3 +b4*x4 +b5*x5
    elseif(x.GE.5.d0) then
      logcool = 3.d0 * x - 21.2968837223113 !cubic extrapolation
    end if

    !opacity according to RA04
    tauCIE = (n(idx_H2) * 1.4285714d-16)**2.8 !note: 1/7d15 = 1.4285714d-16
    cool = p_mass * 1d1**logcool !erg*cm3/s

    cooling_CIE = cool * min(1.d0, (1.d0-exp(-tauCIE))/tauCIE) &
        * n(idx_H2) * sum(n(1:nmols)) !erg/cm3/s

  end function cooling_CIE

  !*****************
  !sigmoid function with x0 shift and s steepness
  function sigmoid(x,x0,s)
    implicit none
    real*8::sigmoid,x,x0,s

    sigmoid = 1d1/(1d1+exp(-s*(x-x0)))

  end function sigmoid

  !*******************
  !window function for H2 cooling to smooth limits
  function wCool(logTgas,logTmin,logTmax)
    implicit none
    real*8::wCool,logTgas,logTmin,logTmax,x

    x = (logTgas-logTmin)/(logTmax-logTmin)
    wCool = 1d1**(2d2*(sigmoid(x,-2d-1,5d1)*sigmoid(-x,-1.2d0,5d1)-1d0))
    if(wCool<1d-199) wCool = 0d0
    if(wCool>1d0) then
      print *,"ERROR: wCool>1"
      stop
    end if

  end function wCool

  !ALL THE COOLING FUNCTIONS ARE FROM GLOVER & ABEL, MNRAS 388, 1627, 2008
  !FOR LOW DENSITY REGIME: CONSIDER AN ORTHO-PARA RATIO OF 3:1
  !UPDATED TO THE DATA REPORTED BY GLOVER 2015, MNRAS
  !EACH SINGLE FUNCTION IS IN erg/s
  !FINAL UNITS = erg/cm3/s
  !*******************************
  function cooling_H2(n, Tgas)
    use krome_commons
    use krome_subs
    use krome_getphys
    real*8::n(:),Tgas
    real*8::temp,logt3,logt,cool,cooling_H2,T3
    real*8::LDL,HDLR,HDLV,HDL
    real*8::logt32,logt33,logt34,logt35,logt36,logt37,logt38
    real*8::dump14,fH2H,fH2e,fH2H2,fH2Hp,fH2He,w14,w24
    integer::i
    character*16::names(nspec)

    temp = Tgas
    cooling_H2 = 0d0
    !if(temp<2d0) return

    T3 = temp * 1.d-3
    logt3 = log10(T3)
    logt = log10(temp)
    cool = 0d0

    logt32 = logt3 * logt3
    logt33 = logt32 * logt3
    logt34 = logt33 * logt3
    logt35 = logt34 * logt3
    logt36 = logt35 * logt3
    logt37 = logt36 * logt3
    logt38 = logt37 * logt3

    w14 = wCool(logt, 1d0, 4d0)
    w24 = wCool(logt, 2d0, 4d0)

    !//H2-H
    if(temp<=1d2) then
      fH2H = 1.d1**(-16.818342D0 +3.7383713D1*logt3 &
          +5.8145166D1*logt32 +4.8656103D1*logt33 &
          +2.0159831D1*logt34 +3.8479610D0*logt35 )*n(idx_H)
    elseif(temp>1d2 .and. temp<=1d3) then
      fH2H = 1.d1**(-2.4311209D1 +3.5692468D0*logt3 &
          -1.1332860D1*logt32 -2.7850082D1*logt33 &
          -2.1328264D1*logt34 -4.2519023D0*logt35)*n(idx_H)
    elseif(temp>1.d3.and.temp<=6d3) then
      fH2H = 1d1**(-2.4311209D1 +4.6450521D0*logt3 &
          -3.7209846D0*logt32 +5.9369081D0*logt33 &
          -5.5108049D0*logt34 +1.5538288D0*logt35)*n(idx_H)
    else
      fH2H = 1.862314467912518E-022*wCool(logt,1d0,log10(6d3))*n(idx_H)
    end if
    cool = cool + fH2H

    !//H2-Hp
    if(temp>1d1.and.temp<=1d4) then
      fH2Hp = 1d1**(-2.2089523d1 +1.5714711d0*logt3 &
          +0.015391166d0*logt32 -0.23619985d0*logt33 &
          -0.51002221d0*logt34 +0.32168730d0*logt35)*n(idx_Hj)
    else
      fH2Hp = 1.182509139382060E-021*n(idx_Hj)*w14
    endif
    cool = cool + fH2Hp

    !//H2-H2
    fH2H2 = w24*1d1**(-2.3962112D1 +2.09433740D0*logt3 &
        -.77151436D0*logt32 +.43693353D0*logt33 &
        -.14913216D0*logt34 -.033638326D0*logt35)*n(idx_H2) !&
        cool = cool + fH2H2

    !//H2-e
    fH2e = 0d0
    if(temp<=5d2) then
      fH2e = 1d1**(min(-2.1928796d1 + 1.6815730d1*logt3 &
          +9.6743155d1*logt32 +3.4319180d2*logt33 &
          +7.3471651d2*logt34 +9.8367576d2*logt35 &
          +8.0181247d2*logt36 +3.6414446d2*logt37 &
          +7.0609154d1*logt38,3d1))*n(idx_e)
    elseif(temp>5d2)  then
      fH2e = 1d1**(-2.2921189D1 +1.6802758D0*logt3 &
          +.93310622D0*logt32 +4.0406627d0*logt33 &
          -4.7274036d0*logt34 -8.8077017d0*logt35 &
          +8.9167183*logt36 + 6.4380698*logt37 &
          -6.3701156*logt38)*n(idx_e)
    end if
    cool = cool + fH2e*w24

    !//H2-He
    if(temp>1d1.and.temp<=1d4)then
      fH2He = 1d1**(-2.3689237d1 +2.1892372d0*logt3&
          -.81520438d0*logt32 +.29036281d0*logt33 -.16596184d0*logt34 &
          +.19191375d0*logt35)*n(idx_He)
    else
      fH2He = 1.002560385050777E-022*n(idx_He)*w14
    endif
    cool = cool + fH2He

    !check error
    if(cool>1.d30) then
      print *,"ERROR: cooling >1.d30 erg/s/cm3"
      print *,"cool (erg/s/cm3): ",cool
      names(:) = get_names()
      do i=1,size(n)
        print '(I3,a18,E11.3)',i,names(i),n(i)
      end do
      stop
    end if

    !this to avoid negative, overflow and useless calculations below
    if(cool<=0d0) then
      cooling_H2 = 0d0
      return
    end if

    !high density limit from HM79, GP98 below Tgas = 2d3
    !UPDATED USING GLOVER 2015 for high temperature corrections, MNRAS
    !IN THE HIGH DENSITY REGIME LAMBDA_H2 = LAMBDA_H2(LTE) = HDL
    !the following mix of functions ensures the right behaviour
    ! at low (T<10 K) and high temperatures (T>2000 K) by
    ! using both the original Hollenbach and the new Glover data
    ! merged in a smooth way.
    if(temp.lt.2d3)then
      HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
          3.e-24*exp(-0.51/t3)) !erg/s
      HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3)) !erg/s
      HDL  = HDLR + HDLV !erg/s
    elseif(temp>=2d3 .and. temp<=1d4)then
      HDL = 1d1**(-2.0584225d1 + 5.0194035*logt3 &
          -1.5738805*logt32 -4.7155769*logt33 &
          +2.4714161*logt34 +5.4710750*logt35 &
          -3.9467356*logt36 -2.2148338*logt37 &
          +1.8161874*logt38)
    else
      dump14 = 1d0 / (1d0 + exp(min((temp-3d4)*2d-4,3d2)))
      HDL = 5.531333679406485E-019*dump14
    endif

    LDL = cool !erg/s
    if (HDL==0.) then
      cooling_H2 = 0.d0
    else
      cooling_H2 = n(idx_H2)/(1.d0/HDL+1.d0/LDL) &
          * min(1.d0, max(1.25d-10 * sum(n(1:nmols)),1d-40)**(-.45)) !erg/cm3/s
    endif

  end function cooling_H2

  !Atomic COOLING  Cen ApJS, 78, 341, 1992
  !UNITS = erg/s/cm3
  !*******************************
  function cooling_Atomic(n, Tgas)
    use krome_commons
    use krome_subs
    real*8::Tgas,cooling_atomic,n(:)
    real*8::temp,T5,cool

    temp = max(Tgas,10d0) !K
    T5 = temp/1d5 !K
    cool = 0d0 !erg/cm3/s

    !COLLISIONAL IONIZATION: H, He, He+, He(2S)
    cool = cool+ 1.27d-21*sqrt(temp)/(1.d0+sqrt(T5))&
        *exp(-1.578091d5/temp)*n(idx_e)*n(idx_H)

    cool = cool+ 9.38d-22*sqrt(temp)/(1.d0+sqrt(T5))&
        *exp(-2.853354d5/temp)*n(idx_e)*n(idx_He)

    cool = cool+ 4.95d-22*sqrt(temp)/(1.d0+sqrt(T5))&
        *exp(-6.31515d5/temp)*n(idx_e)*n(idx_Hej)
    cool = cool+ 5.01d-27*temp**(-0.1687)/(1.d0+sqrt(T5))&
        *exp(-5.5338d4/temp)*n(idx_e)**2*n(idx_Hej)

    !RECOMBINATION: H+, He+,He2+
    cool = cool+ 8.7d-27*sqrt(temp)*(temp/1.d3)**(-0.2)&
        /(1.d0+(temp/1.d6)**0.7)*n(idx_e)*n(idx_Hj)

    cool = cool+ 1.55d-26*temp**(0.3647)*n(idx_e)*n(idx_Hej)

    !DIELECTRONIC RECOMBINATION: He
    cool = cool+ 1.24d-13*temp**(-1.5)*exp(-4.7d5/temp)&
        *(1.d0+0.3d0*exp(-9.4d4/temp))*n(idx_e)*n(idx_Hej)

    !COLLISIONAL EXCITATION:
    !H(all n), He(n=2,3,4 triplets), He+(n=2)
    cool = cool+ 7.5d-19/(1.d0+sqrt(T5))*exp(-1.18348d5/temp)*n(idx_e)*n(idx_H)

    cool = cool+ 9.1d-27*temp**(-.1687)/(1.d0+sqrt(T5))&
        *exp(-1.3179d4/temp)*n(idx_e)**2*n(idx_Hej)
    cool = cool+ 5.54d-17*temp**(-.397)/(1.d0+sqrt(T5))&
        *exp(-4.73638d5/temp)*n(idx_e)*n(idx_Hej)

    cooling_atomic = max(cool, 0d0)  !erg/cm3/s

  end function cooling_Atomic

  !*********************************************
  !function for linear interpolation of f(x), using xval(:)
  ! and the corresponding yval(:) as reference values
  ! note: slow function, use only for initializations
  function flin(xval,yval,x)
    implicit none
    real*8::xval(:),yval(:),x,flin
    integer::i,n
    logical::found
    found = .false.
    n = size(xval)
    x = max(x,xval(1)) !set lower bound
    x = min(x,xval(n)) !set upper bound
    !loop to find interval (slow)
    do i=2,n
      if(x.le.xval(i)) then
        !linear fit
        flin = (yval(i) - yval(i-1)) / (xval(i) - xval(i-1)) * &
            (x - xval(i-1)) + yval(i-1)
        found = .true. !found flag
        exit
      end if
    end do
    if(.not.found) flin = yval(n)

  end function flin

  !************************
  !dump the level populations in a file
  subroutine dump_cooling_pop(Tgas,nfile)
    implicit none
    integer::nfile,i
    real*8::Tgas

    !pop_level_CI(3)
    do i=1,size(pop_level_CI)
      write(nfile,'(a8,I5,3E17.8e3)') "CI", i, Tgas, pop_level_CI(i), sum(pop_level_CI(:))
    end do

    !pop_level_OI(3)
    do i=1,size(pop_level_OI)
      write(nfile,'(a8,I5,3E17.8e3)') "OI", i, Tgas, pop_level_OI(i), sum(pop_level_OI(:))
    end do

    !pop_level_OII(3)
    do i=1,size(pop_level_OII)
      write(nfile,'(a8,I5,3E17.8e3)') "OII", i, Tgas, pop_level_OII(i), sum(pop_level_OII(:))
    end do

    !pop_level_FeII(5)
    do i=1,size(pop_level_FeII)
      write(nfile,'(a8,I5,3E17.8e3)') "FeII", i, Tgas, pop_level_FeII(i), sum(pop_level_FeII(:))
    end do

    !pop_level_CII(2)
    do i=1,size(pop_level_CII)
      write(nfile,'(a8,I5,3E17.8e3)') "CII", i, Tgas, pop_level_CII(i), sum(pop_level_CII(:))
    end do

    !pop_level_FeI(5)
    do i=1,size(pop_level_FeI)
      write(nfile,'(a8,I5,3E17.8e3)') "FeI", i, Tgas, pop_level_FeI(i), sum(pop_level_FeI(:))
    end do

    !pop_level_SiII(2)
    do i=1,size(pop_level_SiII)
      write(nfile,'(a8,I5,3E17.8e3)') "SiII", i, Tgas, pop_level_SiII(i), sum(pop_level_SiII(:))
    end do

    !pop_level_SiI(3)
    do i=1,size(pop_level_SiI)
      write(nfile,'(a8,I5,3E17.8e3)') "SiI", i, Tgas, pop_level_SiI(i), sum(pop_level_SiI(:))
    end do

    write(nfile,*)

  end subroutine dump_cooling_pop

  !***********************
  !metal cooling as in Maio et al. 2007
  ! loaded from data file
  function cooling_Z(n,inTgas)
    use krome_commons
    use krome_constants
    implicit none
    real*8::n(:), inTgas, cool, cooling_Z, k(nZrate), Tgas

    Tgas = inTgas
    k(:) = coolingZ_rate_tabs(Tgas)

    cool = 0d0
    cool = cool + coolingCI(n(:),inTgas,k(:))
    cool = cool + coolingOI(n(:),inTgas,k(:))
    cool = cool + coolingOII(n(:),inTgas,k(:))
    cool = cool + coolingFeII(n(:),inTgas,k(:))
    cool = cool + coolingCII(n(:),inTgas,k(:))
    cool = cool + coolingFeI(n(:),inTgas,k(:))
    cool = cool + coolingSiII(n(:),inTgas,k(:))
    cool = cool + coolingSiI(n(:),inTgas,k(:))

    cooling_Z = cool * boltzmann_erg

  end function cooling_Z

  !********************************
  function coolingZ_rates(inTgas)
    use krome_commons
    use krome_subs
    use krome_fit
    implicit none
    real*8::inTgas,coolingZ_rates(nZrate),k(nZrate)
    real*8::Tgas,invT,logTgas
    integer::i
    real*8::T2,invTgas,lnT

    Tgas = inTgas
    invT = 1d0/Tgas
    logTgas = log10(Tgas)

    T2 = Tgas*1d-2
    invTgas = 1d0/Tgas
    lnT = log(Tgas)

    if(Tgas.ge.1d4) return

    !1->0, CI - H
    k(1) = 1.6d-10*(T2)**(.14)

    !2->0, CI - H
    k(2) = 9.2d-11*(T2)**(.26)

    !2->1, CI - H
    k(3) = 2.9d-10*(T2)**(.26)

    !1->0, CI - H+
    k(4) = (9.6D-11 -1.8D-14*Tgas +1.9D-18*Tgas**2) *Tgas**(.45)

    if(Tgas > 5d3) k(4) = 8.9D-10*Tgas**(.117)

    !2->0, CI - H+
    k(5) = (3.1D-12 -6.D-16*Tgas +3.9d-20*Tgas**2) *Tgas

    if(Tgas > 5d3) k(5) = 2.3D-9*Tgas**(.0965)

    !2->1, CI - H+
    k(6) = (1.D-10 -2.2D-14*Tgas +1.7D-18*Tgas**2) *Tgas**(.70)

    if(Tgas > 5d3) k(6) = 9.2D-9*Tgas**(.0535)

    !1->0, CI - e
    k(7) = 2.88D-6*Tgas**(-.5)*EXP(-9.25141 -7.73782D-1*lnT +3.61184D-1*lnT**2 -1.50892D-2*lnT**3 -6.56325D-4*lnT**4)

    if(Tgas > 1D3) k(7) = 2.88D-6*Tgas**(-.5) *EXP(-4.446D2 -2.27913D2*lnT +4.2595D1*lnT**2 -3.4762*lnT**3 +1.0508D-1*lnT**4)

    !2->0, CI - e
    k(8) = 1.73D-6*Tgas**(-.5)*EXP(-7.69735 -1.30743*lnT +.697638*lnT**2 -.111338*lnT**3 +.705277D-2*lnT**4)

    if(Tgas > 1D3) k(8) = 1.73D-6*Tgas**(-.5)*EXP(3.50609D2 -1.87474D2*lnT +3.61803D1*lnT**2 -3.03283*lnT**3 +9.38138D-2*lnT**4)

    !2->1, CI - e
    k(9) = 1.73D-6*Tgas**(-.5)*EXP(-7.4387 -.57443*lnT +.358264*lnT**2 -4.18166D-2*lnT**3 +2.35272D-3*lnT**4)

    if(Tgas > 1D3) k(9) = 1.73D-6*Tgas**(-.5)*EXP(3.86186D2 -2.02192D2*lnT +3.85049D1*lnT**2 -3.19268*lnT**3 +9.78573D-2*lnT**4)

    !1->0, CI - H2or
    k(10) = 8.7d-11 -6.6d-11*exp(-Tgas/218.3) + 6.6d-11*exp(-2.*Tgas/218.3)

    !2->0, CI - H2or
    k(11) = 1.2d-10 -6.1d-11*exp(-Tgas/387.3)

    !2->1, CI - H2or
    k(12) = 2.9d-10 -1.9d-10*exp(-Tgas/348.9)

    !1->0, CI - H2pa
    k(13) = 7.9D-11 -8.7D-11*EXP(-Tgas/126.4) + 1.3D-10*EXP(-2.*Tgas/126.4)

    !2->0, CI - H2pa
    k(14) = 1.1D-10 -8.6D-11*EXP(-Tgas/223.) + 8.7D-11*EXP(-2.*Tgas/223.)

    !2->1, CI - H2pa
    k(15) = 2.7D-10 -2.6D-10*EXP(-Tgas/250.7) + 1.8D-10*EXP(-2.*Tgas/250.7)

    !1->0, SiI - H
    k(16) = 3.5D-10*(T2)**(-.03)

    !2->0, SiI - H
    k(17) = 1.7D-11*(T2)**(.17)

    !2->1, SiI - H
    k(18) = 5.D-10*(T2)**(.17)

    !1->0, SiI - H+
    k(19) = 7.2d-9

    !2->0, SiI - H+
    k(20) = 7.2d-9

    !2->1, SiI - H+
    k(21) = 7.2d-9

    !1->0, FeI - H
    k(22) = 8.D-10*(T2)**(.17)

    !2->0, FeI - H
    k(23) = 6.9D-10*(T2)**(.17)

    !2->1, FeI - H
    k(24) = 5.3D-10*(T2)**(.17)

    !1->0, FeI - e
    k(25) = 1.2D-7

    !2->0, FeI - e
    k(26) = 1.2D-7

    !2->1, FeI - e
    k(27) = 9.3D-8

    !3->0, FeI - e
    k(28) = 2.D-7*(Tgas/1.d4)**(.57)

    !4->0, FeI - e
    k(29) = 1.D-7*(Tgas/1.d4)**(.57)

    !4->3, FeI - e
    k(30) = 1.5D-7

    !1->0, OI - H
    k(31) = 9.2D-11*(T2)**(.67)

    !2->0, OI - H
    k(32) = 4.3D-11*(T2)**(.80)

    !2->1, OI - H
    k(33) = 1.1D-10*(T2)**(.44)

    !1->0, OI - H+
    k(34) = 6.38D-11*Tgas**(.4)

    if(Tgas > 194.) k(34) = 7.75D-12*Tgas**(.8)

    if(Tgas > 3686.) k(34) = 2.65D-10*Tgas**(.37)

    !2->0, OI - H
    k(35) = 6.1D-13*Tgas**(1.1)

    if(Tgas > 511.) k(35) = 2.12D-12*Tgas**(.9)

    if(Tgas > 7510.) k(35) = 4.49D-10*Tgas**(.3)

    !2->1, OI - H
    k(36) = 2.03D-11*Tgas**(.56)

    if(Tgas > 2090.) k(36) = 3.43D-10*Tgas**(.19)

    !1->0, OI - e
    k(37) = 5.12D-10*Tgas**(-.075)

    !2->0, OI - e
    k(38) = 4.86D-10*Tgas**(-.026)

    !2->1, OI - e
    k(39) = 1.08D-14*Tgas**(.926)

    !1->0, CII - e
    k(40) = 2.8D-7*(Tgas/1d2)**(-.5)

    !1->0, CII - H
    k(41) = 8D-10*(Tgas/1d2)**(.07)

    !1->0, OII - e
    k(42) = 1.3d-8*(Tgas/1.d4)**(-.5)

    !2->0, OII - e
    k(43) = 1.3d-8*(Tgas/1.d4)**(-.5)

    !2->1, OII - e
    k(44) = 2.5d-8*(Tgas/1.d4)**(-.5)

    !1->0, SiII - e
    k(45) = 1.7D-6*(T2)**(-.5)

    !1->0, SiII - H
    k(46) = 8D-10*(T2)**(-.07)

    !1->0, FeII - H
    k(47) = 9.5d-10

    !2->1, FeII - H
    k(48) = 4.7d-10

    !3->2, FeII - H
    k(49) = 5.0d-10

    !4->3, FeII - H
    k(50) = 5.0d-10

    !2->0, FeII - H
    k(51) = 5.7d-10

    !1->0, FeII - e
    k(52) = 1.8D-6*(T2)**(-.5)

    !2->1, FeII - e
    k(53) = 8.7D-7*(T2)**(-.5)

    !3->2, FeII - e
    k(54) = 1.D-5*Tgas**(-.5)

    !4->3, FeII - e
    k(55) = 1.D-5*Tgas**(-.5)

    !2->0, FeII - e
    k(56) = 1.8D-6*(T2)**(-.5)

    !2<-0, CI - H2or
    k(57) = k(11) * 5.0d0 * exp(-6.300000d+01 * invT)

    !1<-0, CI - e
    k(58) = k(7) * 3.0d0 * exp(-2.400000d+01 * invT)

    !2<-0, CI - e
    k(59) = k(8) * 5.0d0 * exp(-6.300000d+01 * invT)

    !1<-0, CI - H2or
    k(60) = k(10) * 3.0d0 * exp(-2.400000d+01 * invT)

    !1<-0, CI - H2pa
    k(61) = k(13) * 3.0d0 * exp(-2.400000d+01 * invT)

    !2<-0, CI - H2pa
    k(62) = k(14) * 5.0d0 * exp(-6.300000d+01 * invT)

    !1<-0, CI - H+
    k(63) = k(4) * 3.0d0 * exp(-2.400000d+01 * invT)

    !2<-1, CI - e
    k(64) = k(9) * 1.66666666667d0 * exp(-3.900000d+01 * invT)

    !2<-1, CI - H
    k(65) = k(3) * 1.66666666667d0 * exp(-3.900000d+01 * invT)

    !2<-0, CI - H
    k(66) = k(2) * 5.0d0 * exp(-6.300000d+01 * invT)

    !2<-1, CI - H+
    k(67) = k(6) * 1.66666666667d0 * exp(-3.900000d+01 * invT)

    !2<-1, CI - H2pa
    k(68) = k(15) * 1.66666666667d0 * exp(-3.900000d+01 * invT)

    !2<-1, CI - H2or
    k(69) = k(12) * 1.66666666667d0 * exp(-3.900000d+01 * invT)

    !1<-0, CI - H
    k(70) = k(1) * 3.0d0 * exp(-2.400000d+01 * invT)

    !2<-0, CI - H+
    k(71) = k(5) * 5.0d0 * exp(-6.300000d+01 * invT)

    !1<-0, OI - H
    k(72) = k(31) * 0.6d0 * exp(-2.300000d+02 * invT)

    !1<-0, OI - e
    k(73) = k(37) * 0.6d0 * exp(-2.300000d+02 * invT)

    !2<-1, OI - e
    k(74) = k(39) * 0.333333333333d0 * exp(-1.000000d+02 * invT)

    !2<-1, OI - H
    k(75) = k(36) * 0.333333333333d0 * exp(-1.000000d+02 * invT)

    !2<-0, OI - H
    k(76) = k(35) * 0.2d0 * exp(-3.300000d+02 * invT)

    !1<-0, OI - H+
    k(77) = k(34) * 0.6d0 * exp(-2.300000d+02 * invT)

    !2<-0, OI - e
    k(78) = k(38) * 0.2d0 * exp(-3.300000d+02 * invT)

    !1<-0, OII - e
    k(79) = k(42) * 3.0d0 * exp(-2.400000d+01 * invT)

    !2<-0, OII - e
    k(80) = k(43) * 5.0d0 * exp(-6.300000d+01 * invT)

    !2<-1, OII - e
    k(81) = k(44) * 1.66666666667d0 * exp(-3.900000d+01 * invT)

    !1<-0, FeII - e
    k(82) = k(52) * 0.8d0 * exp(-5.535800d+02 * invT)

    !2<-0, FeII - e
    k(83) = k(56) * 0.6d0 * exp(-9.605900d+02 * invT)

    !3<-2, FeII - H
    k(84) = k(49) * 0.666666666667d0 * exp(-2.805700d+02 * invT)

    !3<-2, FeII - e
    k(85) = k(54) * 0.666666666667d0 * exp(-2.805700d+02 * invT)

    !4<-3, FeII - H
    k(86) = k(50) * 0.5d0 * exp(-1.646000d+02 * invT)

    !2<-1, FeII - e
    k(87) = k(53) * 0.75d0 * exp(-4.070100d+02 * invT)

    !2<-1, FeII - H
    k(88) = k(48) * 0.75d0 * exp(-4.070100d+02 * invT)

    !2<-0, FeII - H
    k(89) = k(51) * 0.6d0 * exp(-9.605900d+02 * invT)

    !1<-0, FeII - H
    k(90) = k(47) * 0.8d0 * exp(-5.535800d+02 * invT)

    !4<-3, FeII - e
    k(91) = k(55) * 0.5d0 * exp(-1.646000d+02 * invT)

    !1<-0, CII - e
    k(92) = k(40) * 2.0d0 * exp(-9.120000d+01 * invT)

    !1<-0, CII - H
    k(93) = k(41) * 2.0d0 * exp(-9.120000d+01 * invT)

    !1<-0, FeI - H
    k(94) = k(22) * 0.777777777778d0 * exp(-5.984300d+02 * invT)

    !2<-1, FeI - e
    k(95) = k(27) * 0.714285714286d0 * exp(-4.144700d+02 * invT)

    !4<-3, FeI - e
    k(96) = k(30) * 0.727272727273d0 * exp(-6.453000d+02 * invT)

    !1<-0, FeI - e
    k(97) = k(25) * 0.777777777778d0 * exp(-5.984300d+02 * invT)

    !2<-1, FeI - H
    k(98) = k(24) * 0.714285714286d0 * exp(-4.144700d+02 * invT)

    !4<-0, FeI - e
    k(99) = k(29) * 0.888888888889d0 * exp(-1.061350d+04 * invT)

    !3<-0, FeI - e
    k(100) = k(28) * 1.22222222222d0 * exp(-9.968200d+03 * invT)

    !2<-0, FeI - H
    k(101) = k(23) * 0.555555555556d0 * exp(-1.012900d+03 * invT)

    !2<-0, FeI - e
    k(102) = k(26) * 0.555555555556d0 * exp(-1.012900d+03 * invT)

    !1<-0, SiII - e
    k(103) = k(45) * 2.0d0 * exp(-4.136000d+02 * invT)

    !1<-0, SiII - H
    k(104) = k(46) * 2.0d0 * exp(-4.136000d+02 * invT)

    !1<-0, SiI - H+
    k(105) = k(19) * 3.0d0 * exp(-1.100000d+02 * invT)

    !2<-1, SiI - H
    k(106) = k(18) * 1.66666666667d0 * exp(-2.100000d+02 * invT)

    !2<-0, SiI - H
    k(107) = k(17) * 5.0d0 * exp(-3.200000d+02 * invT)

    !2<-1, SiI - H+
    k(108) = k(21) * 1.66666666667d0 * exp(-2.100000d+02 * invT)

    !1<-0, SiI - H
    k(109) = k(16) * 3.0d0 * exp(-1.100000d+02 * invT)

    !2<-0, SiI - H+
    k(110) = k(20) * 5.0d0 * exp(-3.200000d+02 * invT)

    coolingZ_rates(:) = k(:)

    !check rates > 1
    if(maxval(k)>1d0) then
      print *,"ERROR: found rate >1d0 in coolingZ_rates!"
      print *," Tgas =",Tgas
      do i=1,nZrate
        if(k(i)>1d0) print *,i,k(i)
      end do
      stop
    end if

    !check rates <0
    if(minval(k)<0d0) then
      print *,"ERROR: found rate <0d0 in coolingZ_rates!"
      print *," Tgas =",Tgas
      do i=1,nZrate
        if(k(i)<0d0) print *,i,k(i)
      end do
      stop
    end if

  end function coolingZ_rates

  !**********************
  function coolingZ_rate_tabs(inTgas)
    use krome_commons
    implicit none
    real*8::inTgas,Tgas,coolingZ_rate_tabs(nZrate),k(nZrate)
    integer::idx,j
    Tgas = inTgas

    idx = (log10(Tgas)-coolTab_logTlow) * inv_coolTab_idx + 1

    idx = max(idx,1)
    idx = min(idx,coolTab_n-1)

    do j=1,nZrate
      k(j) = (Tgas-coolTab_T(idx)) * inv_coolTab_T(idx) * &
          (coolTab(j,idx+1)-coolTab(j,idx)) + coolTab(j,idx)
      k(j) = max(k(j), 0d0)
    end do

    coolingZ_rate_tabs(:) = k(:)

  end function coolingZ_rate_tabs

  !**********************
  subroutine coolingZ_init_tabs()
    use krome_commons
    implicit none
    integer::j,jmax,idx
    real*8::Tgas,Tgasold

    jmax = coolTab_n !size of the cooling tables (number of saples)

    !note: change upper and lower limit for rate tables here
    coolTab_logTlow = log10(2d0)
    coolTab_logTup = log10(1d8)

    !pre compute this value since used jmax times
    inv_coolTab_idx = (jmax-1) / (coolTab_logTup-coolTab_logTlow)

    !loop over the jmax interpolation points
    do j=1,jmax
      !compute Tgas for the given point
      Tgas = 1d1**((j-1)*(coolTab_logTup-coolTab_logTlow) &
          /(jmax-1) + coolTab_logTlow)
      !produce cooling rates for the given Tgas
      coolTab(:,j) = coolingZ_rates(Tgas)
      !store Tgas into the array
      coolTab_T(j) = Tgas
      !save 1/dT since it is known
      if(j>1) inv_coolTab_T(j-1) = 1d0 / (Tgas-Tgasold)
      Tgasold = Tgas
    end do

  end subroutine coolingZ_init_tabs

  !*******************************
  !this subroutine solves a non linear system
  ! with the equations stored in fcn function
  ! and a dummy jacobian jcn
  subroutine nleq_wrap(x)
    use krome_user_commons
    integer,parameter::nmax=100 !problem size
    integer,parameter::liwk=nmax+50 !size integer workspace
    integer,parameter::lrwk=(nmax+13)*nmax+60 !real workspace
    integer,parameter::luprt=6 !logical unit verbose output
    integer::neq,iopt(50),ierr,niw,nrw,iwk(liwk),ptype,i
    real*8::x(:),xscal(nmax),rtol,rwk(lrwk),idamp,mdamp,xi(size(x)),minx
    real*8::store_invdvdz
    neq = size(x)
    niw = neq+50
    nrw = (neq+13)*neq+60

    ptype = 2 !initial problem type, 2=mildly non-linear
    rtol = 1d-5 !realtive tolerance
    xi(:) = x(:) !store initial guess
    idamp = 1d-4 !initial damp (when ptype>=4, else default)
    mdamp = 1d-8 !minimum damp (when ptype>=4, else default)
    ierr = 0

    !iterate until ierr==0 and non-negative solutions
    do
      if(ptype>50) then
        print *,"ERROR in nleq1: can't find a solution after attempt",ptype
        stop
      end if

      x(:) = xi(:) !restore initial guess

      !if damping error or negative solutions
      ! prepares initial guess with the thin case
      if(ptype>7.and.(ierr==3.or.ierr==0)) then
        rtol = 1d-5
        iwk(:) = 0
        iopt(:) = 0
        rwk(:) = 0d0
        xscal(:) = 0d0
        store_invdvdz = krome_invdvdz !store global variable
        krome_invdvdz = 0d0 !this sets beta to 1
        if(ierr.ne.0) then
          print *,"ERROR in nleq for thin approx",ierr
          stop
        end if
        krome_invdvdz = store_invdvdz !restore global variable
      end if
      xscal(:) = 0d0 !scaling factor
      rtol = 1d-5 !relative tolerance
      iwk(:) = 0 !default iwk
      iwk(31) = int(1e8) !max iterations
      iopt(:) = 0 !default iopt
      iopt(31) = min(ptype,4) !problem type
      rwk(:) = 0d0 !default rwk
      !reduce damps if damping error
      if(ptype>4.and.ierr==3) then
        idamp = idamp * 1d-1 !reduce idamp
        mdamp = mdamp * 1d-1 !reduce mdamp
      end if
      !if problem is extremely nonlinear use custom damps
      if(ptype>4) then
        rwk(21) = idamp !copy idamp to solver
        rwk(22) = mdamp !copy mdamp to solver
      end if

      !check for errors
      if(ierr.ne.0) then
        !print *,"error",ierr
        !problem with damping factor and/or problem type
        if(ierr==3) then
          ptype = ptype + 1 !change the problem type (non-linearity)
        elseif(ierr==5) then
          xi(:) = x(:)
        else
          !other type of error hence stop
          print *,"ERROR in nleq1, ierr:",ierr
          print *,"solutions found so far:"
          do i=1,size(x)
            print *,i,x(i)
          end do
          stop
        end if
      else
        !if succesful search for negative results
        minx = minval(x) !minimum value
        !if minimum value is positive OK
        if(minx.ge.0d0) then
          exit
        else
          !if negative values are small set to zero
          if(abs(minx)/maxval(x)<rtol) then
            do i=1,neq
              x(i) = max(x(i),0d0)
            end do
            exit
          else
            !if large negative values increase non-linearity
            ptype = ptype + 1
          end if
        end if
      end if
    end do
  end subroutine nleq_wrap

  !***************************
  subroutine fcn(n,x,f,ierr)
    implicit none
    integer::n,ierr
    real*8::x(n),f(n)

  end subroutine fcn

  !**********************************
  !dummy jacobian for non linear equation solver
  subroutine jcn()

  end subroutine jcn

  !************************************
  function coolingCI(n,inTgas,k)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    integer::i, hasnegative, nmax
    real*8::coolingCI,n(:),inTgas,k(:)
    real*8::A(3,3),Ain(3,3)
    real*8::B(3),tmp(3)
    real*8::coll_H2or
    real*8::coll_Hj
    real*8::coll_e
    real*8::coll_H2pa
    real*8::coll_H

    !colliders should be >0
    coll_H2or = max(n(idx_H2) * phys_orthoParaRatio / (phys_orthoParaRatio+1d0), 0d0)
    coll_Hj = max(n(idx_Hj), 0d0)
    coll_e = max(n(idx_e), 0d0)
    coll_H2pa = max(n(idx_H2) / (phys_orthoParaRatio+1d0), 0d0)
    coll_H = max(n(idx_H), 0d0)

    !deafault cooling value
    coolingCI = 0d0

    if(n(idx_C)<1d-15) return

    A(:,:) = 0d0

    A(1,1) = 1d0
    A(2,1) = + k(60) * coll_H2or &
        + k(63) * coll_Hj &
        + k(58) * coll_e &
        + k(61) * coll_H2pa &
        + k(70) * coll_H
    A(3,1) = + k(71) * coll_Hj &
        + k(62) * coll_H2pa &
        + k(57) * coll_H2or &
        + k(59) * coll_e &
        + k(66) * coll_H
    A(2,2) = - k(13) * coll_H2pa &
        - k(1) * coll_H &
        - k(7) * coll_e &
        - k(10) * coll_H2or &
        - k(4) * coll_Hj &
        - 7.900000d-08 &
        - k(68) * coll_H2pa &
        - k(69) * coll_H2or &
        - k(64) * coll_e &
        - k(67) * coll_Hj &
        - k(65) * coll_H
    A(3,3) = - k(11) * coll_H2or &
        - k(14) * coll_H2pa &
        - k(2) * coll_H &
        - k(5) * coll_Hj &
        - k(8) * coll_e &
        - 2.100000d-14 &
        - k(9) * coll_e &
        - k(6) * coll_Hj &
        - k(15) * coll_H2pa &
        - k(3) * coll_H &
        - k(12) * coll_H2or &
        - 2.700000d-07
    A(1,2) = 1d0
    A(3,2) = + k(68) * coll_H2pa &
        + k(69) * coll_H2or &
        + k(64) * coll_e &
        + k(67) * coll_Hj &
        + k(65) * coll_H
    A(1,3) = 1d0
    A(2,3) = + k(9) * coll_e &
        + k(6) * coll_Hj &
        + k(15) * coll_H2pa &
        + k(3) * coll_H &
        + k(12) * coll_H2or &
        + 2.700000d-07

    !build matrix B
    B(:) = 0d0
    B(1) = n(idx_C)

    Ain(:,:) = A(:,:)

    call mylin3(A(:,:), B(:))

    !store population
    pop_level_CI(:) = B(:)
    !sanitize negative values
    hasnegative = 0
    do i=1,3
      if(B(i)<0d0) then
        if(abs(B(i)/n(idx_C))>1d-10) then
          hasnegative = 1
        else
          B(i) = 1d-40
        end if
      end if
    end do

    !check if B has negative values
    if(hasnegative>0)then
      print *,"ERROR: minval(B)<0d0 in coolingCI"
      print *,"ntot_CI =", n(idx_C)
      print *,"Tgas =", inTgas
      print *,"B(:) unrolled:"
      do i=1,size(B)
        print *, i, B(i)
      end do
      print *,"A(:,:) min/max:"
      do i=1,size(B)
        print *, i, minval(Ain(i,:)), maxval(Ain(i,:))
      end do

      print *,"A(:,:)"
      do i=1,size(B)
        tmp(:) = Ain(i,:)
        print '(I5,99E17.8)', i, tmp(:)
      end do
      stop
    end if

    coolingCI = B(2) * (7.900000d-08) * 2.400000d+01 &
        + B(3) * (2.700000d-07) * 3.900000d+01 &
        + B(3) * (2.100000d-14) * 6.300000d+01

  end function coolingCI

  !************************************
  function coolingOI(n,inTgas,k)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    integer::i, hasnegative, nmax
    real*8::coolingOI,n(:),inTgas,k(:)
    real*8::A(3,3),Ain(3,3)
    real*8::B(3),tmp(3)
    real*8::coll_e
    real*8::coll_H
    real*8::coll_Hj

    !colliders should be >0
    coll_e = max(n(idx_e), 0d0)
    coll_H = max(n(idx_H), 0d0)
    coll_Hj = max(n(idx_Hj), 0d0)

    !deafault cooling value
    coolingOI = 0d0

    if(n(idx_O)<1d-15) return

    A(:,:) = 0d0

    A(1,1) = 1d0
    A(2,1) = + k(73) * coll_e &
        + k(72) * coll_H &
        + k(77) * coll_Hj
    A(3,1) = + k(78) * coll_e &
        + k(76) * coll_H
    A(2,2) = - k(37) * coll_e &
        - k(34) * coll_Hj &
        - k(31) * coll_H &
        - 8.900000d-05 &
        - k(75) * coll_H &
        - k(74) * coll_e
    A(3,3) = - k(38) * coll_e &
        - k(35) * coll_H &
        - 1.800000d-05 &
        - k(39) * coll_e &
        - k(36) * coll_H &
        - 1.300000d-10
    A(1,2) = 1d0
    A(3,2) = + k(75) * coll_H &
        + k(74) * coll_e
    A(1,3) = 1d0
    A(2,3) = + k(39) * coll_e &
        + k(36) * coll_H &
        + 1.300000d-10

    !build matrix B
    B(:) = 0d0
    B(1) = n(idx_O)

    Ain(:,:) = A(:,:)

    call mylin3(A(:,:), B(:))

    !store population
    pop_level_OI(:) = B(:)
    !sanitize negative values
    hasnegative = 0
    do i=1,3
      if(B(i)<0d0) then
        if(abs(B(i)/n(idx_O))>1d-10) then
          hasnegative = 1
        else
          B(i) = 1d-40
        end if
      end if
    end do

    !check if B has negative values
    if(hasnegative>0)then
      print *,"ERROR: minval(B)<0d0 in coolingOI"
      print *,"ntot_OI =", n(idx_O)
      print *,"Tgas =", inTgas
      print *,"B(:) unrolled:"
      do i=1,size(B)
        print *, i, B(i)
      end do
      print *,"A(:,:) min/max:"
      do i=1,size(B)
        print *, i, minval(Ain(i,:)), maxval(Ain(i,:))
      end do

      print *,"A(:,:)"
      do i=1,size(B)
        tmp(:) = Ain(i,:)
        print '(I5,99E17.8)', i, tmp(:)
      end do
      stop
    end if

    coolingOI = B(2) * (8.900000d-05) * 2.300000d+02 &
        + B(3) * (1.800000d-05) * 3.300000d+02 &
        + B(3) * (1.300000d-10) * 1.000000d+02

  end function coolingOI

  !************************************
  function coolingOII(n,inTgas,k)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    integer::i, hasnegative, nmax
    real*8::coolingOII,n(:),inTgas,k(:)
    real*8::A(3,3),Ain(3,3)
    real*8::B(3),tmp(3)
    real*8::coll_e

    !colliders should be >0
    coll_e = max(n(idx_e), 0d0)

    !deafault cooling value
    coolingOII = 0d0

    if(n(idx_Oj)<1d-15) return

    A(:,:) = 0d0

    A(1,1) = 1d0
    A(2,1) = + k(79) * coll_e
    A(3,1) = + k(80) * coll_e
    A(2,2) = - k(42) * coll_e &
        - 5.100000d-05 &
        - k(81) * coll_e
    A(3,3) = - k(43) * coll_e &
        - 1.700000d-04 &
        - k(44) * coll_e &
        - 1.300000d-07
    A(1,2) = 1d0
    A(3,2) = + k(81) * coll_e
    A(1,3) = 1d0
    A(2,3) = + k(44) * coll_e &
        + 1.300000d-07

    !build matrix B
    B(:) = 0d0
    B(1) = n(idx_Oj)

    Ain(:,:) = A(:,:)

    call mylin3(A(:,:), B(:))

    !store population
    pop_level_OII(:) = B(:)
    !sanitize negative values
    hasnegative = 0
    do i=1,3
      if(B(i)<0d0) then
        if(abs(B(i)/n(idx_Oj))>1d-10) then
          hasnegative = 1
        else
          B(i) = 1d-40
        end if
      end if
    end do

    !check if B has negative values
    if(hasnegative>0)then
      print *,"ERROR: minval(B)<0d0 in coolingOII"
      print *,"ntot_OII =", n(idx_Oj)
      print *,"Tgas =", inTgas
      print *,"B(:) unrolled:"
      do i=1,size(B)
        print *, i, B(i)
      end do
      print *,"A(:,:) min/max:"
      do i=1,size(B)
        print *, i, minval(Ain(i,:)), maxval(Ain(i,:))
      end do

      print *,"A(:,:)"
      do i=1,size(B)
        tmp(:) = Ain(i,:)
        print '(I5,99E17.8)', i, tmp(:)
      end do
      stop
    end if

    coolingOII = B(2) * (5.100000d-05) * 2.400000d+01 &
        + B(3) * (1.700000d-04) * 6.300000d+01 &
        + B(3) * (1.300000d-07) * 3.900000d+01

  end function coolingOII

  !************************************
  function coolingFeII(n,inTgas,k)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    integer::i, hasnegative, nmax
    real*8::coolingFeII,n(:),inTgas,k(:)
    real*8::A(5,5),Ain(5,5)
    real*8::B(5),tmp(5)
    real*8::coll_e
    real*8::coll_H

    !colliders should be >0
    coll_e = max(n(idx_e), 0d0)
    coll_H = max(n(idx_H), 0d0)

    !deafault cooling value
    coolingFeII = 0d0

    if(n(idx_Fej)<1d-15) return

    A(:,:) = 0d0

    A(1,1) = 1d0
    A(2,1) = + k(82) * coll_e &
        + k(90) * coll_H
    A(3,1) = + k(83) * coll_e &
        + k(89) * coll_H
    A(2,2) = - k(47) * coll_H &
        - k(52) * coll_e &
        - 2.130000d-03 &
        - k(87) * coll_e &
        - k(88) * coll_H
    A(3,3) = - k(51) * coll_H &
        - k(56) * coll_e &
        - 1.500000d-09 &
        - k(53) * coll_e &
        - k(48) * coll_H &
        - 1.570000d-03 &
        - k(85) * coll_e &
        - k(84) * coll_H
    A(4,4) = - k(54) * coll_e &
        - k(49) * coll_H &
        - 7.180000d-04 &
        - k(86) * coll_H &
        - k(91) * coll_e
    A(5,5) = - k(50) * coll_H &
        - k(55) * coll_e &
        - 1.880000d-04

    !reduce the size of the problem if possible
    nmax = 1
    do i=5,2,-1
      if(A(i,1)>0d0) then
        nmax = i
        exit
      end if
    end do

    !no need to solve a 1-level problem
    if(nmax==1) return

    A(1,2) = 1d0
    A(3,2) = + k(87) * coll_e &
        + k(88) * coll_H
    A(1,3) = 1d0
    A(2,3) = + k(53) * coll_e &
        + k(48) * coll_H &
        + 1.570000d-03
    A(4,3) = + k(85) * coll_e &
        + k(84) * coll_H
    A(1,4) = 1d0
    A(3,4) = + k(54) * coll_e &
        + k(49) * coll_H &
        + 7.180000d-04
    A(5,4) = + k(86) * coll_H &
        + k(91) * coll_e
    A(1,5) = 1d0
    A(4,5) = + k(50) * coll_H &
        + k(55) * coll_e &
        + 1.880000d-04

    !build matrix B
    B(:) = 0d0
    B(1) = n(idx_Fej)

    Ain(:,:) = A(:,:)

    call mydgesv(nmax, A(:,:), B(:), "coolingFeII")

    !store population
    pop_level_FeII(:) = B(:)
    !sanitize negative values
    hasnegative = 0
    do i=1,5
      if(B(i)<0d0) then
        if(abs(B(i)/n(idx_Fej))>1d-10) then
          hasnegative = 1
        else
          B(i) = 1d-40
        end if
      end if
    end do

    !check if B has negative values
    if(hasnegative>0)then
      print *,"ERROR: minval(B)<0d0 in coolingFeII"
      print *,"ntot_FeII =", n(idx_Fej)
      print *,"Tgas =", inTgas
      print *,"B(:) unrolled:"
      do i=1,size(B)
        print *, i, B(i)
      end do
      print *,"A(:,:) min/max:"
      do i=1,size(B)
        print *, i, minval(Ain(i,:)), maxval(Ain(i,:))
      end do

      print *,"A(:,:)"
      do i=1,size(B)
        tmp(:) = Ain(i,:)
        print '(I5,99E17.8)', i, tmp(:)
      end do
      stop
    end if

    coolingFeII = B(2) * (2.130000d-03) * 5.535800d+02 &
        + B(3) * (1.570000d-03) * 4.070100d+02 &
        + B(5) * (1.880000d-04) * 1.646000d+02 &
        + B(4) * (7.180000d-04) * 2.805700d+02 &
        + B(3) * (1.500000d-09) * 9.605900d+02

  end function coolingFeII

  !************************************
  function coolingCII(n,inTgas,k)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    integer::i, hasnegative, nmax
    real*8::coolingCII,n(:),inTgas,k(:)
    real*8::A(2,2),Ain(2,2)
    real*8::B(2),tmp(2)
    real*8::coll_e
    real*8::coll_H

    !colliders should be >0
    coll_e = max(n(idx_e), 0d0)
    coll_H = max(n(idx_H), 0d0)

    !deafault cooling value
    coolingCII = 0d0

    if(n(idx_Cj)<1d-15) return

    A(:,:) = 0d0

    A(1,1) = 1d0
    A(2,1) = + k(92) * coll_e &
        + k(93) * coll_H
    A(2,2) = - k(40) * coll_e &
        - k(41) * coll_H &
        - 2.400000d-06
    A(1,2) = 1d0

    !build matrix B
    B(:) = 0d0
    B(1) = n(idx_Cj)

    Ain(:,:) = A(:,:)

    call mylin2(A(:,:), B(:))

    !store population
    pop_level_CII(:) = B(:)
    !sanitize negative values
    hasnegative = 0
    do i=1,2
      if(B(i)<0d0) then
        if(abs(B(i)/n(idx_Cj))>1d-10) then
          hasnegative = 1
        else
          B(i) = 1d-40
        end if
      end if
    end do

    !check if B has negative values
    if(hasnegative>0)then
      print *,"ERROR: minval(B)<0d0 in coolingCII"
      print *,"ntot_CII =", n(idx_Cj)
      print *,"Tgas =", inTgas
      print *,"B(:) unrolled:"
      do i=1,size(B)
        print *, i, B(i)
      end do
      print *,"A(:,:) min/max:"
      do i=1,size(B)
        print *, i, minval(Ain(i,:)), maxval(Ain(i,:))
      end do

      print *,"A(:,:)"
      do i=1,size(B)
        tmp(:) = Ain(i,:)
        print '(I5,99E17.8)', i, tmp(:)
      end do
      stop
    end if

    coolingCII = B(2) * (2.400000d-06) * 9.120000d+01

  end function coolingCII

  !************************************
  function coolingFeI(n,inTgas,k)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    integer::i, hasnegative, nmax
    real*8::coolingFeI,n(:),inTgas,k(:)
    real*8::A(5,5),Ain(5,5)
    real*8::B(5),tmp(5)
    real*8::coll_e
    real*8::coll_H

    !colliders should be >0
    coll_e = max(n(idx_e), 0d0)
    coll_H = max(n(idx_H), 0d0)

    !deafault cooling value
    coolingFeI = 0d0

    if(n(idx_Fe)<1d-15) return

    A(:,:) = 0d0

    A(1,1) = 1d0
    A(2,1) = + k(97) * coll_e &
        + k(94) * coll_H
    A(3,1) = + k(102) * coll_e &
        + k(101) * coll_H
    A(4,1) = + k(100) * coll_e
    A(5,1) = + k(99) * coll_e
    A(2,2) = - k(22) * coll_H &
        - k(25) * coll_e &
        - 2.500000d-03 &
        - k(98) * coll_H &
        - k(95) * coll_e
    A(3,3) = - k(26) * coll_e &
        - k(23) * coll_H &
        - 1.600000d-03 &
        - k(27) * coll_e &
        - k(24) * coll_H &
        - 1.000000d-09
    A(4,4) = - k(28) * coll_e &
        - 2.000000d-03 &
        - k(96) * coll_e
    A(5,5) = - k(29) * coll_e &
        - 1.500000d-03 &
        - k(30) * coll_e &
        - 3.600000d-03

    !reduce the size of the problem if possible
    nmax = 1
    do i=5,2,-1
      if(A(i,1)>0d0) then
        nmax = i
        exit
      end if
    end do

    !no need to solve a 1-level problem
    if(nmax==1) return

    A(1,2) = 1d0
    A(3,2) = + k(98) * coll_H &
        + k(95) * coll_e
    A(1,3) = 1d0
    A(2,3) = + k(27) * coll_e &
        + k(24) * coll_H &
        + 1.000000d-09
    A(1,4) = 1d0
    A(5,4) = + k(96) * coll_e
    A(1,5) = 1d0
    A(4,5) = + k(30) * coll_e &
        + 3.600000d-03

    !build matrix B
    B(:) = 0d0
    B(1) = n(idx_Fe)

    Ain(:,:) = A(:,:)

    call mydgesv(nmax, A(:,:), B(:), "coolingFeI")

    !store population
    pop_level_FeI(:) = B(:)
    !sanitize negative values
    hasnegative = 0
    do i=1,5
      if(B(i)<0d0) then
        if(abs(B(i)/n(idx_Fe))>1d-10) then
          hasnegative = 1
        else
          B(i) = 1d-40
        end if
      end if
    end do

    !check if B has negative values
    if(hasnegative>0)then
      print *,"ERROR: minval(B)<0d0 in coolingFeI"
      print *,"ntot_FeI =", n(idx_Fe)
      print *,"Tgas =", inTgas
      print *,"B(:) unrolled:"
      do i=1,size(B)
        print *, i, B(i)
      end do
      print *,"A(:,:) min/max:"
      do i=1,size(B)
        print *, i, minval(Ain(i,:)), maxval(Ain(i,:))
      end do

      print *,"A(:,:)"
      do i=1,size(B)
        tmp(:) = Ain(i,:)
        print '(I5,99E17.8)', i, tmp(:)
      end do
      stop
    end if

    coolingFeI = B(5) * (3.600000d-03) * 6.453000d+02 &
        + B(3) * (1.000000d-09) * 4.144700d+02 &
        + B(4) * (2.000000d-03) * 9.968200d+03 &
        + B(5) * (1.500000d-03) * 1.061350d+04 &
        + B(2) * (2.500000d-03) * 5.984300d+02 &
        + B(3) * (1.600000d-03) * 1.012900d+03

  end function coolingFeI

  !************************************
  function coolingSiII(n,inTgas,k)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    integer::i, hasnegative, nmax
    real*8::coolingSiII,n(:),inTgas,k(:)
    real*8::A(2,2),Ain(2,2)
    real*8::B(2),tmp(2)
    real*8::coll_e
    real*8::coll_H

    !colliders should be >0
    coll_e = max(n(idx_e), 0d0)
    coll_H = max(n(idx_H), 0d0)

    !deafault cooling value
    coolingSiII = 0d0

    if(n(idx_Sij)<1d-15) return

    A(:,:) = 0d0

    A(1,1) = 1d0
    A(2,1) = + k(103) * coll_e &
        + k(104) * coll_H
    A(2,2) = - k(45) * coll_e &
        - k(46) * coll_H &
        - 2.100000d-04
    A(1,2) = 1d0

    !build matrix B
    B(:) = 0d0
    B(1) = n(idx_Sij)

    Ain(:,:) = A(:,:)

    call mylin2(A(:,:), B(:))

    !store population
    pop_level_SiII(:) = B(:)
    !sanitize negative values
    hasnegative = 0
    do i=1,2
      if(B(i)<0d0) then
        if(abs(B(i)/n(idx_Sij))>1d-10) then
          hasnegative = 1
        else
          B(i) = 1d-40
        end if
      end if
    end do

    !check if B has negative values
    if(hasnegative>0)then
      print *,"ERROR: minval(B)<0d0 in coolingSiII"
      print *,"ntot_SiII =", n(idx_Sij)
      print *,"Tgas =", inTgas
      print *,"B(:) unrolled:"
      do i=1,size(B)
        print *, i, B(i)
      end do
      print *,"A(:,:) min/max:"
      do i=1,size(B)
        print *, i, minval(Ain(i,:)), maxval(Ain(i,:))
      end do

      print *,"A(:,:)"
      do i=1,size(B)
        tmp(:) = Ain(i,:)
        print '(I5,99E17.8)', i, tmp(:)
      end do
      stop
    end if

    coolingSiII = B(2) * (2.100000d-04) * 4.136000d+02

  end function coolingSiII

  !************************************
  function coolingSiI(n,inTgas,k)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    integer::i, hasnegative, nmax
    real*8::coolingSiI,n(:),inTgas,k(:)
    real*8::A(3,3),Ain(3,3)
    real*8::B(3),tmp(3)
    real*8::coll_H
    real*8::coll_Hj

    !colliders should be >0
    coll_H = max(n(idx_H), 0d0)
    coll_Hj = max(n(idx_Hj), 0d0)

    !deafault cooling value
    coolingSiI = 0d0

    if(n(idx_Si)<1d-15) return

    A(:,:) = 0d0

    A(1,1) = 1d0
    A(2,1) = + k(109) * coll_H &
        + k(105) * coll_Hj
    A(3,1) = + k(110) * coll_Hj &
        + k(107) * coll_H
    A(2,2) = - k(19) * coll_Hj &
        - k(16) * coll_H &
        - 8.400000d-06 &
        - k(106) * coll_H &
        - k(108) * coll_Hj
    A(3,3) = - k(17) * coll_H &
        - k(20) * coll_Hj &
        - 2.400000d-10 &
        - k(21) * coll_Hj &
        - k(18) * coll_H &
        - 4.200000d-05
    A(1,2) = 1d0
    A(3,2) = + k(106) * coll_H &
        + k(108) * coll_Hj
    A(1,3) = 1d0
    A(2,3) = + k(21) * coll_Hj &
        + k(18) * coll_H &
        + 4.200000d-05

    !build matrix B
    B(:) = 0d0
    B(1) = n(idx_Si)

    Ain(:,:) = A(:,:)

    call mylin3(A(:,:), B(:))

    !store population
    pop_level_SiI(:) = B(:)
    !sanitize negative values
    hasnegative = 0
    do i=1,3
      if(B(i)<0d0) then
        if(abs(B(i)/n(idx_Si))>1d-10) then
          hasnegative = 1
        else
          B(i) = 1d-40
        end if
      end if
    end do

    !check if B has negative values
    if(hasnegative>0)then
      print *,"ERROR: minval(B)<0d0 in coolingSiI"
      print *,"ntot_SiI =", n(idx_Si)
      print *,"Tgas =", inTgas
      print *,"B(:) unrolled:"
      do i=1,size(B)
        print *, i, B(i)
      end do
      print *,"A(:,:) min/max:"
      do i=1,size(B)
        print *, i, minval(Ain(i,:)), maxval(Ain(i,:))
      end do

      print *,"A(:,:)"
      do i=1,size(B)
        tmp(:) = Ain(i,:)
        print '(I5,99E17.8)', i, tmp(:)
      end do
      stop
    end if

    coolingSiI = B(2) * (8.400000d-06) * 1.100000d+02 &
        + B(3) * (2.400000d-10) * 3.200000d+02 &
        + B(3) * (4.200000d-05) * 2.100000d+02

  end function coolingSiI

  !***********************
  subroutine mylin2(a,b)
    !solve Ax=B analytically for a 2-levels system
    implicit none
    integer,parameter::n=2
    real*8::a(n,n),b(n),c(n),iab

    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,1)) then
    !   print *,"ERROR: a22=a21 in mylin2"
    !   stop
    !end if
    iab = b(1)/(a(2,2)-a(2,1))
    c(1) = a(2,2) * iab
    c(2) = -a(2,1) * iab
    b(:) = c(:)

  end subroutine mylin2

  !************************
  subroutine mylin3(a,b)
    !solve Ax=B analytically for a 3-levels system
    implicit none
    integer,parameter::n=3
    real*8::iab,a(n,n),b(n),c(n)

    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,3)) then
    !   print *,"ERROR: a22=a23 in mylin3"
    !   stop
    !end if

    !uncomment this: safer but slower
    !if(a(2,1)*a(3,2)+a(2,2)*a(3,3)+a(2,3)*a(3,1) == &
        !     a(2,1)*a(3,3)+a(2,2)*a(3,1)+a(2,3)*a(3,2)) then
    !   print *,"ERROR: division by zero in mylin3"
    !   stop
    !end if

    iab = b(1) / (a(2,1)*(a(3,3)-a(3,2)) + a(2,2)*(a(3,1)-a(3,3)) &
        + a(2,3)*(a(3,2)-a(3,1)))
    c(1) = (a(2,3)*a(3,2)-a(2,2)*a(3,3)) * iab
    c(2) = -(a(2,3)*a(3,1)-a(2,1)*a(3,3)) * iab
    c(3) = (a(3,1)*a(2,2)-a(2,1)*a(3,2)) * iab
    b(:) = c(:)

  end subroutine mylin3

  !************************************
  subroutine plot_cool(n)
    !routine to plot cooling at runtime
    real*8::n(:),Tgas,Tmin,Tmax
    real*8::cool_atomic,cool_H2,cool_HD,cool_tot, cool_totGP,cool_H2GP
    real*8::cool_dH,cool_Z
    integer::i,imax
    imax = 1000
    Tmin = log10(1d1)
    Tmax = log10(1d8)
    print *,"plotting cooling..."
    open(33,file="KROME_cooling_plot.dat",status="replace")
    do i=1,imax
      Tgas = 1d1**(i*(Tmax-Tmin)/imax+Tmin)
      cool_H2 = 0.d0
      cool_H2GP = 0.d0
      cool_HD = 0.d0
      cool_atomic = 0.d0
      cool_Z = 0.d0
      cool_dH = 0.d0
      cool_H2 = cooling_H2(n(:),Tgas)
      cool_atomic = cooling_atomic(n(:),Tgas)
      cool_Z = cooling_Z(n(:),Tgas)
      cool_tot = cool_H2 + cool_atomic + cool_HD + cool_Z + cool_dH
      cool_totGP = cool_H2GP + cool_atomic + cool_HD + cool_Z + cool_dH
      write(33,'(99E12.3e3)') Tgas, cool_tot, cool_totGP, cool_H2, &
          cool_atomic, cool_HD, cool_H2GP, cool_Z, cool_dH
    end do
    close(33)
    print *,"done!"

  end subroutine plot_cool

  !***********************************
  !routine to dump cooling in unit nfile
  subroutine dump_cool(n,Tgas,nfile)
    use krome_commons
    implicit none
    real*8::Tgas,n(:),cools(ncools)
    integer::nfile

    cools(:) = get_cooling_array(n(:),Tgas)
    write(nfile,'(99E14.5e3)') Tgas, sum(cools), cools(:)

  end subroutine dump_cool

end module KROME_cooling

!############### MODULE ##############
module KROME_heating
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2018-01-25 18:47:28
  !  Changeset 9a6e335
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !************************
  function heating(n,inTgas,k,nH2dust)
    implicit none
    real*8::n(:), Tgas, inTgas, k(:), nH2dust
    real*8::heating

    Tgas = inTgas
    heating = sum(get_heating_array(n(:),Tgas,k(:), nH2dust))

  end function heating

  !*******************************
  function get_heating_array(n, Tgas, k, nH2dust)
    use krome_commons
    implicit none
    real*8::n(:), Tgas, k(:), nH2dust
    real*8::get_heating_array(nheats),heats(nheats)
    real*8::smooth,f1,f2
    !returns heating in erg/cm3/s

    heats(:) = 0.d0

    heats(idx_heat_chem) = heatingChem(n(:), Tgas, k(:), nH2dust)

    f2 = 1.

    heats(idx_heat_custom) = heat_custom(n(:),Tgas)

    get_heating_array(:) = heats(:)

  end function get_heating_array

  !*************************
  function heat_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,heat_custom

    heat_custom = 0d0

  end function heat_custom

  !H2 FORMATION HEATING and other exo/endothermic
  ! processes (including H2 on dust) in erg/cm3/s
  !krome builds the heating/cooling term according
  ! to the chemical network employed
  !*******************************
  function heatingChem(n, Tgas, k, nH2dust)
    use krome_constants
    use krome_commons
    use krome_dust
    use krome_subs
    use krome_getphys
    implicit none
    real*8::heatingChem, n(:), Tgas,k(:),nH2dust
    real*8::h2heatfac,HChem,yH,yH2
    real*8::ncr,ncrn,ncrd1,ncrd2,dd,n2H,small,nmax
    dd = get_Hnuclei(n(:))

    !replace small according to the desired enviroment
    ! and remove nmax if needed
    nmax = maxval(n(1:nmols))
    small = 1d-40/(nmax*nmax*nmax)

    heatingChem = 0.d0

    ncrn  = 1.0d6*(Tgas**(-0.5d0))
    ncrd1 = 1.6d0*exp(-(4.0d2/Tgas)**2)
    ncrd2 = 1.4d0*exp(-1.2d4/(Tgas+1.2d3))

    yH = n(idx_H)/dd   !dimensionless
    yH2= n(idx_H2)/dd  !dimensionless

    ncr = ncrn/(ncrd1*yH+ncrd2*yH2)      !1/cm3
    h2heatfac = 1.0d0/(1.0d0+ncr/dd)     !dimensionless

    HChem = 0.d0 !inits chemical heating
    n2H = n(idx_H) * n(idx_H)

    !H + H2 -> H + H + H (cooling)
    HChem = HChem + k(33) * (-4.48d0*n(idx_H)*n(idx_H2))
    !H2 + E -> H + H + E (cooling)
    HChem = HChem + k(107) * (-4.48d0*n(idx_H2)*n(idx_E))
    !H- + H -> H2 + E (heating)
    HChem = HChem + k(192) * (3.53d0*h2heatfac*n(idx_Hk)*n(idx_H))
    !H + H2+ -> H2 + H+ (heating)
    HChem = HChem + k(196) * (1.83d0*h2heatfac*n(idx_H)*n(idx_H2j))
    !H + H + H -> H2 + H (heating)
    HChem = HChem + k(227) * (4.48d0*h2heatfac*n(idx_H)*n(idx_H)*n(idx_H))
    !H2 + H + H -> H2 + H2 (heating)
    HChem = HChem + k(228) * (4.48d0*h2heatfac*n(idx_H2)*n(idx_H)*n(idx_H))

    heatingChem = HChem * eV_to_erg  !erg/cm3/s

  end function heatingChem

end module KROME_heating

!############### MODULE ##############
module krome_ode
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2018-01-25 18:47:28
  !  Changeset 9a6e335
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  subroutine fex(neq,tt,nin,dn)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_cooling
    use krome_heating
    use krome_tabs
    use krome_photo
    use krome_gadiab
    use krome_getphys
    use krome_phfuncs
    use krome_fit
    implicit none
    integer::neq,idust
    real*8::tt,dn(neq),n(neq),k(nrea),krome_gamma
    real*8::gamma,Tgas,vgas,ntot,nH2dust,nd,nin(neq)
    real*8::rr
    integer::i,r1,r2,r3,p1,p2,p3

    n(:) = nin(:)

    nH2dust = 0.d0
    n(idx_CR) = 1.d0
    n(idx_g)  = 1.d0
    n(idx_dummy) = 1.d0

    dn(:) = 0.d0 !initialize differentials
    n(idx_Tgas) = max(n(idx_tgas),2.73d0)
    n(idx_Tgas) = min(n(idx_tgas),1d9)
    Tgas = n(idx_Tgas) !get temperature

    k(:) = coe_tab(n(:)) !compute coefficients

    !E
    !E
    dn(idx_E) = &
        -k(3)*n(idx_NHj)*n(idx_E) &
        -k(4)*n(idx_HEj)*n(idx_E) &
        -k(13)*n(idx_O)*n(idx_E) &
        -k(25)*n(idx_HEHj)*n(idx_E) &
        -k(29)*n(idx_Hj)*n(idx_E) &
        +k(50)*n(idx_H)*n(idx_Ok) &
        -k(57)*n(idx_H)*n(idx_E) &
        -k(64)*n(idx_MGj)*n(idx_E) &
        +k(69)*n(idx_Hk)*n(idx_N) &
        -k(74)*n(idx_Oj)*n(idx_E) &
        +k(76)*n(idx_H)*n(idx_Sk) &
        -k(87)*n(idx_OHj)*n(idx_E) &
        -k(102)*n(idx_H2j)*n(idx_E) &
        -k(107)*n(idx_H2)*n(idx_E) &
        +k(107)*n(idx_H2)*n(idx_E) &
        -k(109)*n(idx_SIFj)*n(idx_E) &
        -k(115)*n(idx_SIOj)*n(idx_E) &
        +k(129)*n(idx_Hk)*n(idx_C) &
        -k(137)*n(idx_S)*n(idx_E) &
        -k(139)*n(idx_Pj)*n(idx_E) &
        -k(142)*n(idx_FEj)*n(idx_E) &
        -k(157)*n(idx_C)*n(idx_E) &
        +k(170)*n(idx_H)*n(idx_Ck) &
        -k(187)*n(idx_Sj)*n(idx_E) &
        -k(190)*n(idx_NAj)*n(idx_E) &
        +k(192)*n(idx_Hk)*n(idx_H) &
        -k(197)*n(idx_SIj)*n(idx_E) &
        +k(200)*n(idx_CH)*n(idx_O) &
        -k(201)*n(idx_SIHj)*n(idx_E) &
        +k(203)*n(idx_Hk)*n(idx_O) &
        -k(209)*n(idx_CHj)*n(idx_E) &
        -k(220)*n(idx_H)*n(idx_E) &
        +2.d0*k(220)*n(idx_H)*n(idx_E) &
        -k(221)*n(idx_HE)*n(idx_E) &
        +2.d0*k(221)*n(idx_HE)*n(idx_E) &
        -k(223)*n(idx_H2)*n(idx_E) &
        -k(224)*n(idx_Hk)*n(idx_E) &
        +2.d0*k(224)*n(idx_Hk)*n(idx_E) &
        +k(225)*n(idx_Hk)*n(idx_H) &
        +k(226)*n(idx_Hk)*n(idx_Hj) &
        +k(232)*n(idx_C) &
        +k(233)*n(idx_H) &
        +k(234)*n(idx_N) &
        +k(235)*n(idx_CO) &
        +k(236)*n(idx_H2) &
        +k(238)*n(idx_HE) &
        +k(239)*n(idx_O) &
        +k(240)*n(idx_H2) &
        -k(243)*n(idx_Cj)*n(idx_E) &
        -k(244)*n(idx_Nj)*n(idx_E) &
        -k(245)*n(idx_COj)*n(idx_E) &
        -k(248)*n(idx_HSj)*n(idx_E) &
        -k(249)*n(idx_HCOj)*n(idx_E)

    !O-
    !O-
    dn(idx_Ok) = &
        +k(13)*n(idx_O)*n(idx_E) &
        -k(50)*n(idx_H)*n(idx_Ok) &
        -k(63)*n(idx_Ok)*n(idx_MGj) &
        -k(199)*n(idx_Ok)*n(idx_Hj) &
        -k(210)*n(idx_Ok)*n(idx_FEj)

    !H-
    !H-
    dn(idx_Hk) = &
        -k(22)*n(idx_Hk)*n(idx_NAj) &
        -k(24)*n(idx_Hk)*n(idx_FEj) &
        +k(57)*n(idx_H)*n(idx_E) &
        -k(65)*n(idx_Hk)*n(idx_Hj) &
        -k(66)*n(idx_Hk)*n(idx_Sj) &
        -k(69)*n(idx_Hk)*n(idx_N) &
        -k(73)*n(idx_Hk)*n(idx_Oj) &
        -k(129)*n(idx_Hk)*n(idx_C) &
        -k(145)*n(idx_Hk)*n(idx_SIj) &
        -k(192)*n(idx_Hk)*n(idx_H) &
        -k(203)*n(idx_Hk)*n(idx_O) &
        -k(214)*n(idx_Hk)*n(idx_MGj) &
        +k(223)*n(idx_H2)*n(idx_E) &
        -k(224)*n(idx_Hk)*n(idx_E) &
        -k(225)*n(idx_Hk)*n(idx_H) &
        -k(226)*n(idx_Hk)*n(idx_Hj) &
        +k(231)*n(idx_H2)

    !S-
    !S-
    dn(idx_Sk) = &
        -k(76)*n(idx_H)*n(idx_Sk) &
        +k(137)*n(idx_S)*n(idx_E)

    !C-
    !C-
    dn(idx_Ck) = &
        -k(133)*n(idx_Ck)*n(idx_Hj) &
        +k(157)*n(idx_C)*n(idx_E) &
        -k(170)*n(idx_H)*n(idx_Ck)

    !H
    !H
    dn(idx_H) = &
        -k(1)*n(idx_H)*n(idx_HEj) &
        +k(3)*n(idx_NHj)*n(idx_E) &
        -k(6)*n(idx_H)*n(idx_HSj) &
        +k(7)*n(idx_Hj)*n(idx_SI) &
        +k(8)*n(idx_Hj)*n(idx_OH) &
        +k(9)*n(idx_N)*n(idx_HS) &
        -k(11)*n(idx_H)*n(idx_HEHj) &
        +k(12)*n(idx_H2)*n(idx_HS) &
        -k(15)*n(idx_H)*n(idx_SIj) &
        -k(16)*n(idx_H)*n(idx_NS) &
        -k(18)*n(idx_H)*n(idx_S2) &
        -k(21)*n(idx_H)*n(idx_SO) &
        +k(22)*n(idx_Hk)*n(idx_NAj) &
        +k(23)*n(idx_CH)*n(idx_S) &
        +k(24)*n(idx_Hk)*n(idx_FEj) &
        +k(25)*n(idx_HEHj)*n(idx_E) &
        +k(27)*n(idx_OH)*n(idx_S) &
        -k(28)*n(idx_H)*n(idx_O2) &
        +k(28)*n(idx_H)*n(idx_O2) &
        +k(29)*n(idx_Hj)*n(idx_E) &
        -k(30)*n(idx_H)*n(idx_NO) &
        +k(32)*n(idx_H2j)*n(idx_HE) &
        -k(33)*n(idx_H)*n(idx_H2) &
        +3.d0*k(33)*n(idx_H)*n(idx_H2) &
        -k(36)*n(idx_H)*n(idx_CO) &
        +k(37)*n(idx_H2)*n(idx_F) &
        -k(38)*n(idx_H)*n(idx_CH2) &
        -k(40)*n(idx_H)*n(idx_NH) &
        +k(42)*n(idx_H2)*n(idx_CH) &
        -k(43)*n(idx_H)*n(idx_HCN) &
        +k(44)*n(idx_OH)*n(idx_CO) &
        +k(45)*n(idx_NH)*n(idx_O) &
        +k(48)*n(idx_N)*n(idx_NH) &
        -k(50)*n(idx_H)*n(idx_Ok) &
        +k(51)*n(idx_OH)*n(idx_SIO) &
        -k(53)*n(idx_H)*n(idx_NH2) &
        +k(54)*n(idx_H2)*n(idx_C) &
        +k(55)*n(idx_H2)*n(idx_Oj) &
        -k(56)*n(idx_H)*n(idx_OCN) &
        -k(57)*n(idx_H)*n(idx_E) &
        +k(58)*n(idx_H2j)*n(idx_O) &
        -k(61)*n(idx_Hj)*n(idx_H) &
        -k(62)*n(idx_H)*n(idx_HS) &
        +2.d0*k(65)*n(idx_Hk)*n(idx_Hj) &
        +k(66)*n(idx_Hk)*n(idx_Sj) &
        -k(67)*n(idx_H)*n(idx_OH) &
        +2.d0*k(67)*n(idx_H)*n(idx_OH) &
        -k(70)*n(idx_H)*n(idx_NO) &
        +k(73)*n(idx_Hk)*n(idx_Oj) &
        -k(76)*n(idx_H)*n(idx_Sk) &
        +k(81)*n(idx_OH)*n(idx_SIj) &
        -k(83)*n(idx_H)*n(idx_O) &
        +k(84)*n(idx_Hj)*n(idx_P) &
        +k(87)*n(idx_OHj)*n(idx_E) &
        -k(88)*n(idx_H)*n(idx_C) &
        +k(90)*n(idx_C)*n(idx_HS) &
        +k(91)*n(idx_H2)*n(idx_S) &
        -k(92)*n(idx_H)*n(idx_CH) &
        -k(93)*n(idx_H)*n(idx_O2) &
        -k(95)*n(idx_H)*n(idx_HCO) &
        +k(100)*n(idx_HF)*n(idx_SIj) &
        +2.d0*k(102)*n(idx_H2j)*n(idx_E) &
        -k(103)*n(idx_H)*n(idx_H2O) &
        +2.d0*k(103)*n(idx_H)*n(idx_H2O) &
        +k(104)*n(idx_H2)*n(idx_NH) &
        +2.d0*k(107)*n(idx_H2)*n(idx_E) &
        -k(108)*n(idx_H)*n(idx_CH) &
        +2.d0*k(108)*n(idx_H)*n(idx_CH) &
        -k(110)*n(idx_H)*n(idx_SIHj) &
        +k(114)*n(idx_NH)*n(idx_S) &
        +k(117)*n(idx_H2)*n(idx_N) &
        +k(120)*n(idx_O)*n(idx_OH) &
        +k(127)*n(idx_O)*n(idx_HCN) &
        +k(133)*n(idx_Ck)*n(idx_Hj) &
        -k(135)*n(idx_H)*n(idx_OCS) &
        +k(140)*n(idx_Hj)*n(idx_SIO) &
        +k(141)*n(idx_O)*n(idx_HS) &
        +k(144)*n(idx_Hj)*n(idx_MG) &
        +k(145)*n(idx_Hk)*n(idx_SIj) &
        -k(146)*n(idx_H)*n(idx_H2O) &
        +k(148)*n(idx_OH)*n(idx_SO) &
        +k(149)*n(idx_H2)*n(idx_OH) &
        -k(150)*n(idx_H)*n(idx_OCN) &
        +k(152)*n(idx_C)*n(idx_CH) &
        -k(153)*n(idx_H)*n(idx_C2) &
        -k(154)*n(idx_H)*n(idx_OCN) &
        +k(159)*n(idx_H2)*n(idx_O) &
        +k(162)*n(idx_C)*n(idx_NH) &
        -k(163)*n(idx_H)*n(idx_H2S) &
        +k(165)*n(idx_CH)*n(idx_N) &
        +k(167)*n(idx_Hj)*n(idx_S) &
        -k(170)*n(idx_H)*n(idx_Ck) &
        +k(172)*n(idx_Hj)*n(idx_FE) &
        +k(174)*n(idx_OH)*n(idx_CS) &
        -k(176)*n(idx_H)*n(idx_CO2) &
        +k(178)*n(idx_H2)*n(idx_Sj) &
        +k(181)*n(idx_CH)*n(idx_O) &
        -k(182)*n(idx_H)*n(idx_Oj) &
        +k(183)*n(idx_N)*n(idx_OH) &
        +k(184)*n(idx_Hj)*n(idx_O) &
        -k(185)*n(idx_H)*n(idx_SO) &
        +k(191)*n(idx_S)*n(idx_HS) &
        -k(192)*n(idx_Hk)*n(idx_H) &
        -k(194)*n(idx_H)*n(idx_CHj) &
        -k(196)*n(idx_H)*n(idx_H2j) &
        +k(199)*n(idx_Ok)*n(idx_Hj) &
        +k(201)*n(idx_SIHj)*n(idx_E) &
        +k(204)*n(idx_C)*n(idx_OH) &
        +k(205)*n(idx_OH)*n(idx_SI) &
        +k(206)*n(idx_H2j)*n(idx_C) &
        -k(208)*n(idx_H)*n(idx_OH) &
        +k(209)*n(idx_CHj)*n(idx_E) &
        +k(214)*n(idx_Hk)*n(idx_MGj) &
        +k(215)*n(idx_H2)*n(idx_CN) &
        -k(217)*n(idx_H)*n(idx_NS) &
        +k(219)*n(idx_OH)*n(idx_CN) &
        -k(220)*n(idx_H)*n(idx_E) &
        +k(222)*n(idx_H2)*n(idx_Hj) &
        +k(223)*n(idx_H2)*n(idx_E) &
        +k(224)*n(idx_Hk)*n(idx_E) &
        -k(225)*n(idx_Hk)*n(idx_H) &
        +2.d0*k(225)*n(idx_Hk)*n(idx_H) &
        -3.d0*k(227)*n(idx_H)*n(idx_H)*n(idx_H) &
        +k(227)*n(idx_H)*n(idx_H)*n(idx_H) &
        -2.d0*k(228)*n(idx_H2)*n(idx_H)*n(idx_H) &
        -2.d0*k(229)*n(idx_H)*n(idx_H)*n(idx_HE) &
        +k(230)*n(idx_Hj)*n(idx_NA) &
        -k(233)*n(idx_H) &
        +k(236)*n(idx_H2) &
        +2.d0*k(237)*n(idx_H2) &
        +k(247)*n(idx_Hj)*n(idx_NH) &
        +k(248)*n(idx_HSj)*n(idx_E) &
        +k(249)*n(idx_HCOj)*n(idx_E) &
        +k(250)*n(idx_HEj)*n(idx_HF)

    !HE
    !HE
    dn(idx_HE) = &
        +k(1)*n(idx_H)*n(idx_HEj) &
        +k(4)*n(idx_HEj)*n(idx_E) &
        +k(10)*n(idx_HEj)*n(idx_SI) &
        +k(11)*n(idx_H)*n(idx_HEHj) &
        +k(25)*n(idx_HEHj)*n(idx_E) &
        -k(32)*n(idx_H2j)*n(idx_HE) &
        -k(221)*n(idx_HE)*n(idx_E) &
        -k(229)*n(idx_H)*n(idx_H)*n(idx_HE) &
        +k(229)*n(idx_H)*n(idx_H)*n(idx_HE) &
        -k(238)*n(idx_HE) &
        +k(246)*n(idx_HEj)*n(idx_SIO2) &
        +k(250)*n(idx_HEj)*n(idx_HF)

    !C
    !C
    dn(idx_C) = &
        -k(2)*n(idx_C)*n(idx_NO) &
        -k(34)*n(idx_C)*n(idx_CS) &
        +k(36)*n(idx_H)*n(idx_CO) &
        -k(46)*n(idx_C)*n(idx_CO) &
        -k(47)*n(idx_C)*n(idx_CN) &
        -k(52)*n(idx_C)*n(idx_HCOj) &
        -k(54)*n(idx_H2)*n(idx_C) &
        +k(72)*n(idx_O)*n(idx_CN) &
        -k(75)*n(idx_H2)*n(idx_C) &
        -k(80)*n(idx_C)*n(idx_SO) &
        +k(82)*n(idx_O)*n(idx_CS) &
        +k(85)*n(idx_CN)*n(idx_S) &
        -k(88)*n(idx_H)*n(idx_C) &
        +k(89)*n(idx_SI)*n(idx_CO) &
        -k(90)*n(idx_C)*n(idx_HS) &
        +k(92)*n(idx_H)*n(idx_CH) &
        +k(98)*n(idx_CH)*n(idx_O) &
        +k(99)*n(idx_CH)*n(idx_S) &
        -k(101)*n(idx_C)*n(idx_NO) &
        -k(105)*n(idx_C)*n(idx_NH) &
        -k(106)*n(idx_C)*n(idx_N2) &
        +k(108)*n(idx_H)*n(idx_CH) &
        +k(111)*n(idx_Cj)*n(idx_MG) &
        -k(112)*n(idx_C)*n(idx_O) &
        +k(121)*n(idx_Cj)*n(idx_SI) &
        +k(124)*n(idx_Cj)*n(idx_FE) &
        +k(128)*n(idx_N)*n(idx_CN) &
        -k(129)*n(idx_Hk)*n(idx_C) &
        +k(131)*n(idx_CH)*n(idx_N) &
        +k(133)*n(idx_Ck)*n(idx_Hj) &
        +k(151)*n(idx_N)*n(idx_C2) &
        -k(152)*n(idx_C)*n(idx_CH) &
        +k(153)*n(idx_H)*n(idx_C2) &
        -k(157)*n(idx_C)*n(idx_E) &
        -k(162)*n(idx_C)*n(idx_NH) &
        -2.d0*k(166)*n(idx_C)*n(idx_C) &
        +k(168)*n(idx_O)*n(idx_C2) &
        -k(169)*n(idx_C)*n(idx_OH) &
        -k(171)*n(idx_C)*n(idx_HS) &
        -k(173)*n(idx_C)*n(idx_N) &
        +k(177)*n(idx_C2)*n(idx_S) &
        -k(189)*n(idx_C)*n(idx_SO2) &
        -k(193)*n(idx_C)*n(idx_SO) &
        -k(198)*n(idx_C)*n(idx_NS) &
        -k(202)*n(idx_C)*n(idx_SIOj) &
        -k(204)*n(idx_C)*n(idx_OH) &
        -k(206)*n(idx_H2j)*n(idx_C) &
        -k(207)*n(idx_C)*n(idx_S) &
        +k(209)*n(idx_CHj)*n(idx_E) &
        -k(216)*n(idx_C)*n(idx_O2) &
        -k(232)*n(idx_C) &
        +k(242)*n(idx_CO) &
        +k(243)*n(idx_Cj)*n(idx_E) &
        +k(245)*n(idx_COj)*n(idx_E)

    !NO
    !NO
    dn(idx_NO) = &
        -k(2)*n(idx_C)*n(idx_NO) &
        -k(30)*n(idx_H)*n(idx_NO) &
        -k(35)*n(idx_N)*n(idx_NO) &
        +k(45)*n(idx_NH)*n(idx_O) &
        +k(60)*n(idx_N)*n(idx_O2) &
        -k(70)*n(idx_H)*n(idx_NO) &
        +k(72)*n(idx_O)*n(idx_CN) &
        -k(101)*n(idx_C)*n(idx_NO) &
        -k(123)*n(idx_SI)*n(idx_NO) &
        +k(156)*n(idx_N)*n(idx_SIOj) &
        +k(175)*n(idx_N)*n(idx_CO2) &
        +k(180)*n(idx_O)*n(idx_N2) &
        +k(183)*n(idx_N)*n(idx_OH) &
        +k(211)*n(idx_N)*n(idx_SO) &
        +k(212)*n(idx_O)*n(idx_NS) &
        +k(255)*n(idx_N)*n(idx_PO)

    !CO
    !CO
    dn(idx_CO) = &
        +k(2)*n(idx_C)*n(idx_NO) &
        -k(36)*n(idx_H)*n(idx_CO) &
        -k(44)*n(idx_OH)*n(idx_CO) &
        -k(46)*n(idx_C)*n(idx_CO) &
        +k(52)*n(idx_C)*n(idx_HCOj) &
        +k(71)*n(idx_SI)*n(idx_HCOj) &
        +k(77)*n(idx_OH)*n(idx_CS) &
        -k(89)*n(idx_SI)*n(idx_CO) &
        +k(95)*n(idx_H)*n(idx_HCO) &
        +k(112)*n(idx_C)*n(idx_O) &
        +k(116)*n(idx_O)*n(idx_HCN) &
        +k(134)*n(idx_SI)*n(idx_CO2) &
        +k(135)*n(idx_H)*n(idx_OCS) &
        +k(154)*n(idx_H)*n(idx_OCN) &
        +k(160)*n(idx_O)*n(idx_CS) &
        +k(168)*n(idx_O)*n(idx_C2) &
        +k(175)*n(idx_N)*n(idx_CO2) &
        +k(176)*n(idx_H)*n(idx_CO2) &
        +k(181)*n(idx_CH)*n(idx_O) &
        +k(189)*n(idx_C)*n(idx_SO2) &
        +k(193)*n(idx_C)*n(idx_SO) &
        +k(202)*n(idx_C)*n(idx_SIOj) &
        +k(204)*n(idx_C)*n(idx_OH) &
        +k(216)*n(idx_C)*n(idx_O2) &
        +k(218)*n(idx_O)*n(idx_CN) &
        -k(235)*n(idx_CO) &
        -k(242)*n(idx_CO) &
        +k(249)*n(idx_HCOj)*n(idx_E)

    !N
    !N
    dn(idx_N) = &
        +k(2)*n(idx_C)*n(idx_NO) &
        +k(3)*n(idx_NHj)*n(idx_E) &
        -k(9)*n(idx_N)*n(idx_HS) &
        +k(16)*n(idx_H)*n(idx_NS) &
        -k(17)*n(idx_N)*n(idx_CS) &
        +k(30)*n(idx_H)*n(idx_NO) &
        -k(35)*n(idx_N)*n(idx_NO) &
        +k(40)*n(idx_H)*n(idx_NH) &
        +k(41)*n(idx_NH)*n(idx_O) &
        +k(47)*n(idx_C)*n(idx_CN) &
        -k(48)*n(idx_N)*n(idx_NH) &
        -k(60)*n(idx_N)*n(idx_O2) &
        -k(69)*n(idx_Hk)*n(idx_N) &
        -k(94)*n(idx_N)*n(idx_SO) &
        +k(105)*n(idx_C)*n(idx_NH) &
        +k(106)*n(idx_C)*n(idx_N2) &
        -k(117)*n(idx_H2)*n(idx_N) &
        +k(123)*n(idx_SI)*n(idx_NO) &
        -k(128)*n(idx_N)*n(idx_CN) &
        -k(131)*n(idx_CH)*n(idx_N) &
        -k(151)*n(idx_N)*n(idx_C2) &
        -k(155)*n(idx_N)*n(idx_OH) &
        -k(156)*n(idx_N)*n(idx_SIOj) &
        +k(158)*n(idx_NH)*n(idx_S) &
        -k(165)*n(idx_CH)*n(idx_N) &
        -k(173)*n(idx_C)*n(idx_N) &
        -k(175)*n(idx_N)*n(idx_CO2) &
        -k(179)*n(idx_N)*n(idx_HS) &
        +k(180)*n(idx_O)*n(idx_N2) &
        -k(183)*n(idx_N)*n(idx_OH) &
        -k(211)*n(idx_N)*n(idx_SO) &
        +k(218)*n(idx_O)*n(idx_CN) &
        -k(234)*n(idx_N) &
        +2.d0*k(241)*n(idx_N2) &
        +k(244)*n(idx_Nj)*n(idx_E) &
        -k(252)*n(idx_N)*n(idx_PN) &
        -k(253)*n(idx_N)*n(idx_PO) &
        -k(255)*n(idx_N)*n(idx_PO)

    !O2
    !O2
    dn(idx_O2) = &
        -k(5)*n(idx_O2)*n(idx_S) &
        -k(26)*n(idx_H2)*n(idx_O2) &
        -k(28)*n(idx_H)*n(idx_O2) &
        +k(31)*n(idx_O)*n(idx_SO2) &
        +k(39)*n(idx_O)*n(idx_O) &
        -k(49)*n(idx_SI)*n(idx_O2) &
        -k(60)*n(idx_N)*n(idx_O2) &
        -k(93)*n(idx_H)*n(idx_O2) &
        +k(97)*n(idx_O)*n(idx_SIOj) &
        -k(113)*n(idx_CN)*n(idx_O2) &
        +k(120)*n(idx_O)*n(idx_OH) &
        +k(122)*n(idx_O)*n(idx_SO) &
        -k(216)*n(idx_C)*n(idx_O2) &
        +k(246)*n(idx_HEj)*n(idx_SIO2) &
        -k(254)*n(idx_P)*n(idx_O2)

    !S
    !S
    dn(idx_S) = &
        -k(5)*n(idx_O2)*n(idx_S) &
        +k(14)*n(idx_Sj)*n(idx_FE) &
        +k(17)*n(idx_N)*n(idx_CS) &
        +k(18)*n(idx_H)*n(idx_S2) &
        +k(19)*n(idx_NA)*n(idx_Sj) &
        +k(21)*n(idx_H)*n(idx_SO) &
        -k(23)*n(idx_CH)*n(idx_S) &
        -k(27)*n(idx_OH)*n(idx_S) &
        +k(34)*n(idx_C)*n(idx_CS) &
        +k(59)*n(idx_MG)*n(idx_Sj) &
        +k(62)*n(idx_H)*n(idx_HS) &
        +k(66)*n(idx_Hk)*n(idx_Sj) &
        -k(85)*n(idx_CN)*n(idx_S) &
        +k(86)*n(idx_O)*n(idx_HS) &
        -k(91)*n(idx_H2)*n(idx_S) &
        -k(99)*n(idx_CH)*n(idx_S) &
        -k(114)*n(idx_NH)*n(idx_S) &
        +k(119)*n(idx_SI)*n(idx_Sj) &
        +k(122)*n(idx_O)*n(idx_SO) &
        -k(125)*n(idx_S)*n(idx_SO2) &
        -k(137)*n(idx_S)*n(idx_E) &
        -k(158)*n(idx_NH)*n(idx_S) &
        +k(160)*n(idx_O)*n(idx_CS) &
        -k(167)*n(idx_Hj)*n(idx_S) &
        +k(171)*n(idx_C)*n(idx_HS) &
        -k(177)*n(idx_C2)*n(idx_S) &
        +k(179)*n(idx_N)*n(idx_HS) &
        +k(187)*n(idx_Sj)*n(idx_E) &
        +k(188)*n(idx_HS)*n(idx_HS) &
        -k(191)*n(idx_S)*n(idx_HS) &
        +k(193)*n(idx_C)*n(idx_SO) &
        +k(198)*n(idx_C)*n(idx_NS) &
        -k(207)*n(idx_C)*n(idx_S) &
        +k(211)*n(idx_N)*n(idx_SO) &
        +k(212)*n(idx_O)*n(idx_NS) &
        +k(217)*n(idx_H)*n(idx_NS) &
        +k(248)*n(idx_HSj)*n(idx_E)

    !SO
    !SO
    dn(idx_SO) = &
        +k(5)*n(idx_O2)*n(idx_S) &
        -k(21)*n(idx_H)*n(idx_SO) &
        +k(27)*n(idx_OH)*n(idx_S) &
        +k(31)*n(idx_O)*n(idx_SO2) &
        -k(80)*n(idx_C)*n(idx_SO) &
        +k(82)*n(idx_O)*n(idx_CS) &
        -k(94)*n(idx_N)*n(idx_SO) &
        -k(122)*n(idx_O)*n(idx_SO) &
        +2.d0*k(125)*n(idx_S)*n(idx_SO2) &
        +k(141)*n(idx_O)*n(idx_HS) &
        -k(148)*n(idx_OH)*n(idx_SO) &
        -k(185)*n(idx_H)*n(idx_SO) &
        +k(189)*n(idx_C)*n(idx_SO2) &
        -k(193)*n(idx_C)*n(idx_SO) &
        -k(211)*n(idx_N)*n(idx_SO)

    !O
    !O
    dn(idx_O) = &
        +k(5)*n(idx_O2)*n(idx_S) &
        -k(13)*n(idx_O)*n(idx_E) &
        +k(20)*n(idx_OH)*n(idx_F) &
        +2.d0*k(28)*n(idx_H)*n(idx_O2) &
        -k(31)*n(idx_O)*n(idx_SO2) &
        +k(35)*n(idx_N)*n(idx_NO) &
        -2.d0*k(39)*n(idx_O)*n(idx_O) &
        -k(41)*n(idx_NH)*n(idx_O) &
        -k(45)*n(idx_NH)*n(idx_O) &
        +k(46)*n(idx_C)*n(idx_CO) &
        +k(49)*n(idx_SI)*n(idx_O2) &
        -k(58)*n(idx_H2j)*n(idx_O) &
        +k(60)*n(idx_N)*n(idx_O2) &
        +k(63)*n(idx_Ok)*n(idx_MGj) &
        +k(67)*n(idx_H)*n(idx_OH) &
        +k(68)*n(idx_Oj)*n(idx_FE) &
        +k(70)*n(idx_H)*n(idx_NO) &
        -k(72)*n(idx_O)*n(idx_CN) &
        +k(73)*n(idx_Hk)*n(idx_Oj) &
        +k(74)*n(idx_Oj)*n(idx_E) &
        +k(80)*n(idx_C)*n(idx_SO) &
        -k(82)*n(idx_O)*n(idx_CS) &
        -k(83)*n(idx_H)*n(idx_O) &
        -k(86)*n(idx_O)*n(idx_HS) &
        +k(87)*n(idx_OHj)*n(idx_E) &
        +k(93)*n(idx_H)*n(idx_O2) &
        +k(94)*n(idx_N)*n(idx_SO) &
        +k(96)*n(idx_OH)*n(idx_CN) &
        -k(97)*n(idx_O)*n(idx_SIOj) &
        -k(98)*n(idx_CH)*n(idx_O) &
        +k(101)*n(idx_C)*n(idx_NO) &
        -k(112)*n(idx_C)*n(idx_O) &
        +k(113)*n(idx_CN)*n(idx_O2) &
        +k(115)*n(idx_SIOj)*n(idx_E) &
        -k(116)*n(idx_O)*n(idx_HCN) &
        -k(118)*n(idx_O)*n(idx_HCN) &
        -k(120)*n(idx_O)*n(idx_OH) &
        -k(122)*n(idx_O)*n(idx_SO) &
        -k(127)*n(idx_O)*n(idx_HCN) &
        -k(141)*n(idx_O)*n(idx_HS) &
        +k(147)*n(idx_OH)*n(idx_OH) &
        +k(150)*n(idx_H)*n(idx_OCN) &
        +k(155)*n(idx_N)*n(idx_OH) &
        -k(159)*n(idx_H2)*n(idx_O) &
        -k(160)*n(idx_O)*n(idx_CS) &
        -k(168)*n(idx_O)*n(idx_C2) &
        +k(169)*n(idx_C)*n(idx_OH) &
        -k(180)*n(idx_O)*n(idx_N2) &
        -k(181)*n(idx_CH)*n(idx_O) &
        +k(182)*n(idx_H)*n(idx_Oj) &
        -k(184)*n(idx_Hj)*n(idx_O) &
        +k(185)*n(idx_H)*n(idx_SO) &
        -k(186)*n(idx_O)*n(idx_SI) &
        +k(199)*n(idx_Ok)*n(idx_Hj) &
        -k(200)*n(idx_CH)*n(idx_O) &
        -k(203)*n(idx_Hk)*n(idx_O) &
        +k(208)*n(idx_H)*n(idx_OH) &
        +k(210)*n(idx_Ok)*n(idx_FEj) &
        -k(212)*n(idx_O)*n(idx_NS) &
        -k(213)*n(idx_O)*n(idx_H2O) &
        +k(216)*n(idx_C)*n(idx_O2) &
        -k(218)*n(idx_O)*n(idx_CN) &
        -k(239)*n(idx_O) &
        +k(242)*n(idx_CO) &
        +k(245)*n(idx_COj)*n(idx_E) &
        +k(253)*n(idx_N)*n(idx_PO) &
        +k(254)*n(idx_P)*n(idx_O2)

    !H2
    !H2
    dn(idx_H2) = &
        +k(6)*n(idx_H)*n(idx_HSj) &
        -k(12)*n(idx_H2)*n(idx_HS) &
        -k(26)*n(idx_H2)*n(idx_O2) &
        -k(33)*n(idx_H)*n(idx_H2) &
        -k(37)*n(idx_H2)*n(idx_F) &
        +k(38)*n(idx_H)*n(idx_CH2) &
        +k(40)*n(idx_H)*n(idx_NH) &
        -k(42)*n(idx_H2)*n(idx_CH) &
        +k(43)*n(idx_H)*n(idx_HCN) &
        +k(53)*n(idx_H)*n(idx_NH2) &
        -k(54)*n(idx_H2)*n(idx_C) &
        -k(55)*n(idx_H2)*n(idx_Oj) &
        +k(62)*n(idx_H)*n(idx_HS) &
        -k(75)*n(idx_H2)*n(idx_C) &
        -k(91)*n(idx_H2)*n(idx_S) &
        +k(92)*n(idx_H)*n(idx_CH) &
        +k(95)*n(idx_H)*n(idx_HCO) &
        -k(104)*n(idx_H2)*n(idx_NH) &
        -k(107)*n(idx_H2)*n(idx_E) &
        +k(110)*n(idx_H)*n(idx_SIHj) &
        -k(117)*n(idx_H2)*n(idx_N) &
        +k(146)*n(idx_H)*n(idx_H2O) &
        -k(149)*n(idx_H2)*n(idx_OH) &
        -k(159)*n(idx_H2)*n(idx_O) &
        +k(163)*n(idx_H)*n(idx_H2S) &
        -k(178)*n(idx_H2)*n(idx_Sj) &
        +k(192)*n(idx_Hk)*n(idx_H) &
        +k(194)*n(idx_H)*n(idx_CHj) &
        +k(196)*n(idx_H)*n(idx_H2j) &
        +k(208)*n(idx_H)*n(idx_OH) &
        -k(215)*n(idx_H2)*n(idx_CN) &
        -k(222)*n(idx_H2)*n(idx_Hj) &
        -k(223)*n(idx_H2)*n(idx_E) &
        +k(227)*n(idx_H)*n(idx_H)*n(idx_H) &
        -k(228)*n(idx_H2)*n(idx_H)*n(idx_H) &
        +2.d0*k(228)*n(idx_H2)*n(idx_H)*n(idx_H) &
        +k(229)*n(idx_H)*n(idx_H)*n(idx_HE) &
        -k(231)*n(idx_H2) &
        -k(236)*n(idx_H2) &
        -k(237)*n(idx_H2) &
        -k(240)*n(idx_H2) &
        -k(251)*n(idx_H2)*n(idx_Fj)

    !SI
    !SI
    dn(idx_SI) = &
        -k(7)*n(idx_Hj)*n(idx_SI) &
        -k(10)*n(idx_HEj)*n(idx_SI) &
        -k(49)*n(idx_SI)*n(idx_O2) &
        -k(71)*n(idx_SI)*n(idx_HCOj) &
        -k(89)*n(idx_SI)*n(idx_CO) &
        +k(109)*n(idx_SIFj)*n(idx_E) &
        +k(115)*n(idx_SIOj)*n(idx_E) &
        -k(119)*n(idx_SI)*n(idx_Sj) &
        -k(121)*n(idx_Cj)*n(idx_SI) &
        -k(123)*n(idx_SI)*n(idx_NO) &
        +k(130)*n(idx_MG)*n(idx_SIj) &
        +k(132)*n(idx_NA)*n(idx_SIj) &
        -k(134)*n(idx_SI)*n(idx_CO2) &
        +k(136)*n(idx_SIj)*n(idx_FE) &
        -k(143)*n(idx_SI)*n(idx_Pj) &
        +k(145)*n(idx_Hk)*n(idx_SIj) &
        -k(186)*n(idx_O)*n(idx_SI) &
        +k(197)*n(idx_SIj)*n(idx_E) &
        +k(201)*n(idx_SIHj)*n(idx_E) &
        -k(205)*n(idx_OH)*n(idx_SI)

    !OH
    !OH
    dn(idx_OH) = &
        -k(8)*n(idx_Hj)*n(idx_OH) &
        -k(20)*n(idx_OH)*n(idx_F) &
        +k(21)*n(idx_H)*n(idx_SO) &
        +2.d0*k(26)*n(idx_H2)*n(idx_O2) &
        -k(27)*n(idx_OH)*n(idx_S) &
        +k(30)*n(idx_H)*n(idx_NO) &
        +k(36)*n(idx_H)*n(idx_CO) &
        +k(41)*n(idx_NH)*n(idx_O) &
        -k(44)*n(idx_OH)*n(idx_CO) &
        +k(50)*n(idx_H)*n(idx_Ok) &
        -k(51)*n(idx_OH)*n(idx_SIO) &
        +k(56)*n(idx_H)*n(idx_OCN) &
        -k(67)*n(idx_H)*n(idx_OH) &
        -k(77)*n(idx_OH)*n(idx_CS) &
        -k(81)*n(idx_OH)*n(idx_SIj) &
        +k(83)*n(idx_H)*n(idx_O) &
        +k(86)*n(idx_O)*n(idx_HS) &
        +k(93)*n(idx_H)*n(idx_O2) &
        -k(96)*n(idx_OH)*n(idx_CN) &
        +k(98)*n(idx_CH)*n(idx_O) &
        +k(103)*n(idx_H)*n(idx_H2O) &
        +k(118)*n(idx_O)*n(idx_HCN) &
        -k(120)*n(idx_O)*n(idx_OH) &
        +k(146)*n(idx_H)*n(idx_H2O) &
        -2.d0*k(147)*n(idx_OH)*n(idx_OH) &
        -k(148)*n(idx_OH)*n(idx_SO) &
        -k(149)*n(idx_H2)*n(idx_OH) &
        -k(155)*n(idx_N)*n(idx_OH) &
        +k(159)*n(idx_H2)*n(idx_O) &
        -k(161)*n(idx_OH)*n(idx_H2S) &
        -k(169)*n(idx_C)*n(idx_OH) &
        -k(174)*n(idx_OH)*n(idx_CS) &
        +k(176)*n(idx_H)*n(idx_CO2) &
        -k(183)*n(idx_N)*n(idx_OH) &
        +k(203)*n(idx_Hk)*n(idx_O) &
        -k(204)*n(idx_C)*n(idx_OH) &
        -k(205)*n(idx_OH)*n(idx_SI) &
        -k(208)*n(idx_H)*n(idx_OH) &
        +2.d0*k(213)*n(idx_O)*n(idx_H2O) &
        -k(219)*n(idx_OH)*n(idx_CN)

    !HS
    !HS
    dn(idx_HS) = &
        -k(9)*n(idx_N)*n(idx_HS) &
        -k(12)*n(idx_H2)*n(idx_HS) &
        +k(16)*n(idx_H)*n(idx_NS) &
        +k(18)*n(idx_H)*n(idx_S2) &
        -k(62)*n(idx_H)*n(idx_HS) &
        +k(76)*n(idx_H)*n(idx_Sk) &
        +k(77)*n(idx_OH)*n(idx_CS) &
        -k(86)*n(idx_O)*n(idx_HS) &
        -k(90)*n(idx_C)*n(idx_HS) &
        +k(91)*n(idx_H2)*n(idx_S) &
        +k(99)*n(idx_CH)*n(idx_S) &
        +k(135)*n(idx_H)*n(idx_OCS) &
        -k(141)*n(idx_O)*n(idx_HS) &
        +k(158)*n(idx_NH)*n(idx_S) &
        +k(161)*n(idx_OH)*n(idx_H2S) &
        +k(163)*n(idx_H)*n(idx_H2S) &
        -k(171)*n(idx_C)*n(idx_HS) &
        -k(179)*n(idx_N)*n(idx_HS) &
        +k(185)*n(idx_H)*n(idx_SO) &
        -2.d0*k(188)*n(idx_HS)*n(idx_HS) &
        -k(191)*n(idx_S)*n(idx_HS)

    !NS
    !NS
    dn(idx_NS) = &
        +k(9)*n(idx_N)*n(idx_HS) &
        -k(16)*n(idx_H)*n(idx_NS) &
        +k(85)*n(idx_CN)*n(idx_S) &
        +k(94)*n(idx_N)*n(idx_SO) &
        +k(114)*n(idx_NH)*n(idx_S) &
        -k(198)*n(idx_C)*n(idx_NS) &
        -k(212)*n(idx_O)*n(idx_NS) &
        -k(217)*n(idx_H)*n(idx_NS)

    !H2S
    !H2S
    dn(idx_H2S) = &
        +k(12)*n(idx_H2)*n(idx_HS) &
        -k(161)*n(idx_OH)*n(idx_H2S) &
        -k(163)*n(idx_H)*n(idx_H2S) &
        +k(188)*n(idx_HS)*n(idx_HS)

    !FE
    !FE
    dn(idx_FE) = &
        -k(14)*n(idx_Sj)*n(idx_FE) &
        +k(24)*n(idx_Hk)*n(idx_FEj) &
        -k(68)*n(idx_Oj)*n(idx_FE) &
        +k(78)*n(idx_NA)*n(idx_FEj) &
        -k(79)*n(idx_SIOj)*n(idx_FE) &
        -k(124)*n(idx_Cj)*n(idx_FE) &
        -k(136)*n(idx_SIj)*n(idx_FE) &
        -k(138)*n(idx_HCOj)*n(idx_FE) &
        +k(142)*n(idx_FEj)*n(idx_E) &
        -k(172)*n(idx_Hj)*n(idx_FE) &
        +k(210)*n(idx_Ok)*n(idx_FEj)

    !CS
    !CS
    dn(idx_CS) = &
        -k(17)*n(idx_N)*n(idx_CS) &
        +k(23)*n(idx_CH)*n(idx_S) &
        -k(34)*n(idx_C)*n(idx_CS) &
        -k(77)*n(idx_OH)*n(idx_CS) &
        +k(80)*n(idx_C)*n(idx_SO) &
        -k(82)*n(idx_O)*n(idx_CS) &
        +k(90)*n(idx_C)*n(idx_HS) &
        -k(160)*n(idx_O)*n(idx_CS) &
        -k(174)*n(idx_OH)*n(idx_CS) &
        +k(177)*n(idx_C2)*n(idx_S) &
        +k(207)*n(idx_C)*n(idx_S)

    !CN
    !CN
    dn(idx_CN) = &
        +k(17)*n(idx_N)*n(idx_CS) &
        +k(43)*n(idx_H)*n(idx_HCN) &
        -k(47)*n(idx_C)*n(idx_CN) &
        +k(56)*n(idx_H)*n(idx_OCN) &
        -k(72)*n(idx_O)*n(idx_CN) &
        -k(85)*n(idx_CN)*n(idx_S) &
        -k(96)*n(idx_OH)*n(idx_CN) &
        +k(101)*n(idx_C)*n(idx_NO) &
        +k(106)*n(idx_C)*n(idx_N2) &
        -k(113)*n(idx_CN)*n(idx_O2) &
        +k(118)*n(idx_O)*n(idx_HCN) &
        -k(128)*n(idx_N)*n(idx_CN) &
        +k(151)*n(idx_N)*n(idx_C2) &
        +k(162)*n(idx_C)*n(idx_NH) &
        +k(165)*n(idx_CH)*n(idx_N) &
        +k(173)*n(idx_C)*n(idx_N) &
        +k(198)*n(idx_C)*n(idx_NS) &
        -k(215)*n(idx_H2)*n(idx_CN) &
        -k(218)*n(idx_O)*n(idx_CN) &
        -k(219)*n(idx_OH)*n(idx_CN)

    !S2
    !S2
    dn(idx_S2) = &
        -k(18)*n(idx_H)*n(idx_S2) &
        +k(191)*n(idx_S)*n(idx_HS)

    !NA
    !NA
    dn(idx_NA) = &
        -k(19)*n(idx_NA)*n(idx_Sj) &
        +k(22)*n(idx_Hk)*n(idx_NAj) &
        -k(78)*n(idx_NA)*n(idx_FEj) &
        -k(132)*n(idx_NA)*n(idx_SIj) &
        +k(190)*n(idx_NAj)*n(idx_E) &
        -k(195)*n(idx_NA)*n(idx_MGj) &
        -k(230)*n(idx_Hj)*n(idx_NA)

    !F
    !F
    dn(idx_F) = &
        -k(20)*n(idx_OH)*n(idx_F) &
        -k(37)*n(idx_H2)*n(idx_F) &
        +k(109)*n(idx_SIFj)*n(idx_E) &
        +k(251)*n(idx_H2)*n(idx_Fj)

    !HF
    !HF
    dn(idx_HF) = &
        +k(20)*n(idx_OH)*n(idx_F) &
        +k(37)*n(idx_H2)*n(idx_F) &
        -k(100)*n(idx_HF)*n(idx_SIj) &
        -k(250)*n(idx_HEj)*n(idx_HF)

    !CH
    !CH
    dn(idx_CH) = &
        -k(23)*n(idx_CH)*n(idx_S) &
        +k(38)*n(idx_H)*n(idx_CH2) &
        -k(42)*n(idx_H2)*n(idx_CH) &
        +k(54)*n(idx_H2)*n(idx_C) &
        +k(88)*n(idx_H)*n(idx_C) &
        -k(92)*n(idx_H)*n(idx_CH) &
        -k(98)*n(idx_CH)*n(idx_O) &
        -k(99)*n(idx_CH)*n(idx_S) &
        +k(105)*n(idx_C)*n(idx_NH) &
        -k(108)*n(idx_H)*n(idx_CH) &
        +k(129)*n(idx_Hk)*n(idx_C) &
        -k(131)*n(idx_CH)*n(idx_N) &
        -k(152)*n(idx_C)*n(idx_CH) &
        +k(153)*n(idx_H)*n(idx_C2) &
        -k(165)*n(idx_CH)*n(idx_N) &
        +k(169)*n(idx_C)*n(idx_OH) &
        +k(170)*n(idx_H)*n(idx_Ck) &
        +k(171)*n(idx_C)*n(idx_HS) &
        -k(181)*n(idx_CH)*n(idx_O) &
        -k(200)*n(idx_CH)*n(idx_O)

    !SO2
    !SO2
    dn(idx_SO2) = &
        -k(31)*n(idx_O)*n(idx_SO2) &
        -k(125)*n(idx_S)*n(idx_SO2) &
        +k(148)*n(idx_OH)*n(idx_SO) &
        -k(189)*n(idx_C)*n(idx_SO2)

    !C2
    !C2
    dn(idx_C2) = &
        +k(34)*n(idx_C)*n(idx_CS) &
        +k(46)*n(idx_C)*n(idx_CO) &
        +k(47)*n(idx_C)*n(idx_CN) &
        -k(151)*n(idx_N)*n(idx_C2) &
        +k(152)*n(idx_C)*n(idx_CH) &
        -k(153)*n(idx_H)*n(idx_C2) &
        +k(166)*n(idx_C)*n(idx_C) &
        -k(168)*n(idx_O)*n(idx_C2) &
        -k(177)*n(idx_C2)*n(idx_S)

    !N2
    !N2
    dn(idx_N2) = &
        +k(35)*n(idx_N)*n(idx_NO) &
        +k(48)*n(idx_N)*n(idx_NH) &
        -k(106)*n(idx_C)*n(idx_N2) &
        +k(128)*n(idx_N)*n(idx_CN) &
        -k(180)*n(idx_O)*n(idx_N2) &
        -k(241)*n(idx_N2) &
        +k(252)*n(idx_N)*n(idx_PN)

    !CH2
    !CH2
    dn(idx_CH2) = &
        -k(38)*n(idx_H)*n(idx_CH2) &
        +k(42)*n(idx_H2)*n(idx_CH) &
        +k(75)*n(idx_H2)*n(idx_C)

    !NH
    !NH
    dn(idx_NH) = &
        -k(40)*n(idx_H)*n(idx_NH) &
        -k(41)*n(idx_NH)*n(idx_O) &
        -k(45)*n(idx_NH)*n(idx_O) &
        -k(48)*n(idx_N)*n(idx_NH) &
        +k(53)*n(idx_H)*n(idx_NH2) &
        +k(69)*n(idx_Hk)*n(idx_N) &
        +k(70)*n(idx_H)*n(idx_NO) &
        -k(104)*n(idx_H2)*n(idx_NH) &
        -k(105)*n(idx_C)*n(idx_NH) &
        -k(114)*n(idx_NH)*n(idx_S) &
        +k(116)*n(idx_O)*n(idx_HCN) &
        +k(117)*n(idx_H2)*n(idx_N) &
        +k(131)*n(idx_CH)*n(idx_N) &
        +k(154)*n(idx_H)*n(idx_OCN) &
        +k(155)*n(idx_N)*n(idx_OH) &
        -k(158)*n(idx_NH)*n(idx_S) &
        -k(162)*n(idx_C)*n(idx_NH) &
        +k(179)*n(idx_N)*n(idx_HS) &
        +k(217)*n(idx_H)*n(idx_NS) &
        -k(247)*n(idx_Hj)*n(idx_NH)

    !HCN
    !HCN
    dn(idx_HCN) = &
        -k(43)*n(idx_H)*n(idx_HCN) &
        +k(96)*n(idx_OH)*n(idx_CN) &
        -k(116)*n(idx_O)*n(idx_HCN) &
        -k(118)*n(idx_O)*n(idx_HCN) &
        -k(127)*n(idx_O)*n(idx_HCN) &
        +k(150)*n(idx_H)*n(idx_OCN) &
        +k(215)*n(idx_H2)*n(idx_CN)

    !CO2
    !CO2
    dn(idx_CO2) = &
        +k(44)*n(idx_OH)*n(idx_CO) &
        -k(134)*n(idx_SI)*n(idx_CO2) &
        -k(175)*n(idx_N)*n(idx_CO2) &
        -k(176)*n(idx_H)*n(idx_CO2)

    !SIO
    !SIO
    dn(idx_SIO) = &
        +k(49)*n(idx_SI)*n(idx_O2) &
        -k(51)*n(idx_OH)*n(idx_SIO) &
        +k(79)*n(idx_SIOj)*n(idx_FE) &
        +k(89)*n(idx_SI)*n(idx_CO) &
        +k(123)*n(idx_SI)*n(idx_NO) &
        +k(126)*n(idx_MG)*n(idx_SIOj) &
        +k(134)*n(idx_SI)*n(idx_CO2) &
        -k(140)*n(idx_Hj)*n(idx_SIO) &
        +k(186)*n(idx_O)*n(idx_SI) &
        +k(205)*n(idx_OH)*n(idx_SI)

    !SIO2
    !SIO2
    dn(idx_SIO2) = &
        +k(51)*n(idx_OH)*n(idx_SIO) &
        -k(246)*n(idx_HEj)*n(idx_SIO2)

    !NH2
    !NH2
    dn(idx_NH2) = &
        -k(53)*n(idx_H)*n(idx_NH2) &
        +k(104)*n(idx_H2)*n(idx_NH)

    !OCN
    !OCN
    dn(idx_OCN) = &
        -k(56)*n(idx_H)*n(idx_OCN) &
        +k(113)*n(idx_CN)*n(idx_O2) &
        +k(127)*n(idx_O)*n(idx_HCN) &
        -k(150)*n(idx_H)*n(idx_OCN) &
        -k(154)*n(idx_H)*n(idx_OCN) &
        +k(219)*n(idx_OH)*n(idx_CN)

    !MG
    !MG
    dn(idx_MG) = &
        -k(59)*n(idx_MG)*n(idx_Sj) &
        +k(63)*n(idx_Ok)*n(idx_MGj) &
        +k(64)*n(idx_MGj)*n(idx_E) &
        -k(111)*n(idx_Cj)*n(idx_MG) &
        -k(126)*n(idx_MG)*n(idx_SIOj) &
        -k(130)*n(idx_MG)*n(idx_SIj) &
        -k(144)*n(idx_Hj)*n(idx_MG) &
        -k(164)*n(idx_MG)*n(idx_HCOj) &
        +k(195)*n(idx_NA)*n(idx_MGj) &
        +k(214)*n(idx_Hk)*n(idx_MGj)

    !P
    !P
    dn(idx_P) = &
        -k(84)*n(idx_Hj)*n(idx_P) &
        +k(139)*n(idx_Pj)*n(idx_E) &
        +k(143)*n(idx_SI)*n(idx_Pj) &
        +k(252)*n(idx_N)*n(idx_PN) &
        -k(254)*n(idx_P)*n(idx_O2) &
        +k(255)*n(idx_N)*n(idx_PO)

    !HCO
    !HCO
    dn(idx_HCO) = &
        -k(95)*n(idx_H)*n(idx_HCO) &
        +k(138)*n(idx_HCOj)*n(idx_FE) &
        +k(164)*n(idx_MG)*n(idx_HCOj)

    !H2O
    !H2O
    dn(idx_H2O) = &
        -k(103)*n(idx_H)*n(idx_H2O) &
        -k(146)*n(idx_H)*n(idx_H2O) &
        +k(147)*n(idx_OH)*n(idx_OH) &
        +k(149)*n(idx_H2)*n(idx_OH) &
        +k(161)*n(idx_OH)*n(idx_H2S) &
        -k(213)*n(idx_O)*n(idx_H2O)

    !OCS
    !OCS
    dn(idx_OCS) = &
        -k(135)*n(idx_H)*n(idx_OCS) &
        +k(174)*n(idx_OH)*n(idx_CS)

    !PN
    !PN
    dn(idx_PN) = &
        -k(252)*n(idx_N)*n(idx_PN) &
        +k(253)*n(idx_N)*n(idx_PO)

    !PO
    !PO
    dn(idx_PO) = &
        -k(253)*n(idx_N)*n(idx_PO) &
        +k(254)*n(idx_P)*n(idx_O2) &
        -k(255)*n(idx_N)*n(idx_PO)

    !HE+
    !HE+
    dn(idx_HEj) = &
        -k(1)*n(idx_H)*n(idx_HEj) &
        -k(4)*n(idx_HEj)*n(idx_E) &
        -k(10)*n(idx_HEj)*n(idx_SI) &
        +k(221)*n(idx_HE)*n(idx_E) &
        +k(238)*n(idx_HE) &
        -k(246)*n(idx_HEj)*n(idx_SIO2) &
        -k(250)*n(idx_HEj)*n(idx_HF)

    !H+
    !H+
    dn(idx_Hj) = &
        +k(1)*n(idx_H)*n(idx_HEj) &
        -k(7)*n(idx_Hj)*n(idx_SI) &
        -k(8)*n(idx_Hj)*n(idx_OH) &
        -k(29)*n(idx_Hj)*n(idx_E) &
        -k(61)*n(idx_Hj)*n(idx_H) &
        -k(65)*n(idx_Hk)*n(idx_Hj) &
        -k(84)*n(idx_Hj)*n(idx_P) &
        -k(133)*n(idx_Ck)*n(idx_Hj) &
        -k(140)*n(idx_Hj)*n(idx_SIO) &
        -k(144)*n(idx_Hj)*n(idx_MG) &
        -k(167)*n(idx_Hj)*n(idx_S) &
        -k(172)*n(idx_Hj)*n(idx_FE) &
        +k(182)*n(idx_H)*n(idx_Oj) &
        -k(184)*n(idx_Hj)*n(idx_O) &
        +k(196)*n(idx_H)*n(idx_H2j) &
        -k(199)*n(idx_Ok)*n(idx_Hj) &
        +k(220)*n(idx_H)*n(idx_E) &
        -k(222)*n(idx_H2)*n(idx_Hj) &
        -k(226)*n(idx_Hk)*n(idx_Hj) &
        -k(230)*n(idx_Hj)*n(idx_NA) &
        +k(231)*n(idx_H2) &
        +k(233)*n(idx_H) &
        +k(236)*n(idx_H2) &
        -k(247)*n(idx_Hj)*n(idx_NH)

    !NH+
    !NH+
    dn(idx_NHj) = &
        -k(3)*n(idx_NHj)*n(idx_E) &
        +k(247)*n(idx_Hj)*n(idx_NH)

    !HS+
    !HS+
    dn(idx_HSj) = &
        -k(6)*n(idx_H)*n(idx_HSj) &
        +k(178)*n(idx_H2)*n(idx_Sj) &
        -k(248)*n(idx_HSj)*n(idx_E)

    !S+
    !S+
    dn(idx_Sj) = &
        +k(6)*n(idx_H)*n(idx_HSj) &
        -k(14)*n(idx_Sj)*n(idx_FE) &
        -k(19)*n(idx_NA)*n(idx_Sj) &
        -k(59)*n(idx_MG)*n(idx_Sj) &
        -k(66)*n(idx_Hk)*n(idx_Sj) &
        -k(119)*n(idx_SI)*n(idx_Sj) &
        +k(167)*n(idx_Hj)*n(idx_S) &
        -k(178)*n(idx_H2)*n(idx_Sj) &
        -k(187)*n(idx_Sj)*n(idx_E)

    !SI+
    !SI+
    dn(idx_SIj) = &
        +k(7)*n(idx_Hj)*n(idx_SI) &
        +k(10)*n(idx_HEj)*n(idx_SI) &
        -k(15)*n(idx_H)*n(idx_SIj) &
        -k(81)*n(idx_OH)*n(idx_SIj) &
        +k(97)*n(idx_O)*n(idx_SIOj) &
        -k(100)*n(idx_HF)*n(idx_SIj) &
        +k(110)*n(idx_H)*n(idx_SIHj) &
        +k(119)*n(idx_SI)*n(idx_Sj) &
        +k(121)*n(idx_Cj)*n(idx_SI) &
        -k(130)*n(idx_MG)*n(idx_SIj) &
        -k(132)*n(idx_NA)*n(idx_SIj) &
        -k(136)*n(idx_SIj)*n(idx_FE) &
        +k(143)*n(idx_SI)*n(idx_Pj) &
        -k(145)*n(idx_Hk)*n(idx_SIj) &
        +k(156)*n(idx_N)*n(idx_SIOj) &
        -k(197)*n(idx_SIj)*n(idx_E) &
        +k(202)*n(idx_C)*n(idx_SIOj) &
        +k(246)*n(idx_HEj)*n(idx_SIO2)

    !OH+
    !OH+
    dn(idx_OHj) = &
        +k(8)*n(idx_Hj)*n(idx_OH) &
        +k(55)*n(idx_H2)*n(idx_Oj) &
        +k(58)*n(idx_H2j)*n(idx_O) &
        -k(87)*n(idx_OHj)*n(idx_E)

    !HEH+
    !HEH+
    dn(idx_HEHj) = &
        -k(11)*n(idx_H)*n(idx_HEHj) &
        -k(25)*n(idx_HEHj)*n(idx_E) &
        +k(32)*n(idx_H2j)*n(idx_HE)

    !H2+
    !H2+
    dn(idx_H2j) = &
        +k(11)*n(idx_H)*n(idx_HEHj) &
        -k(32)*n(idx_H2j)*n(idx_HE) &
        -k(58)*n(idx_H2j)*n(idx_O) &
        +k(61)*n(idx_Hj)*n(idx_H) &
        -k(102)*n(idx_H2j)*n(idx_E) &
        -k(196)*n(idx_H)*n(idx_H2j) &
        -k(206)*n(idx_H2j)*n(idx_C) &
        +k(222)*n(idx_H2)*n(idx_Hj) &
        +k(226)*n(idx_Hk)*n(idx_Hj) &
        +k(240)*n(idx_H2) &
        +k(251)*n(idx_H2)*n(idx_Fj)

    !FE+
    !FE+
    dn(idx_FEj) = &
        +k(14)*n(idx_Sj)*n(idx_FE) &
        -k(24)*n(idx_Hk)*n(idx_FEj) &
        +k(68)*n(idx_Oj)*n(idx_FE) &
        -k(78)*n(idx_NA)*n(idx_FEj) &
        +k(79)*n(idx_SIOj)*n(idx_FE) &
        +k(124)*n(idx_Cj)*n(idx_FE) &
        +k(136)*n(idx_SIj)*n(idx_FE) &
        +k(138)*n(idx_HCOj)*n(idx_FE) &
        -k(142)*n(idx_FEj)*n(idx_E) &
        +k(172)*n(idx_Hj)*n(idx_FE) &
        -k(210)*n(idx_Ok)*n(idx_FEj)

    !SIH+
    !SIH+
    dn(idx_SIHj) = &
        +k(15)*n(idx_H)*n(idx_SIj) &
        +k(71)*n(idx_SI)*n(idx_HCOj) &
        -k(110)*n(idx_H)*n(idx_SIHj) &
        -k(201)*n(idx_SIHj)*n(idx_E)

    !NA+
    !NA+
    dn(idx_NAj) = &
        +k(19)*n(idx_NA)*n(idx_Sj) &
        -k(22)*n(idx_Hk)*n(idx_NAj) &
        +k(78)*n(idx_NA)*n(idx_FEj) &
        +k(132)*n(idx_NA)*n(idx_SIj) &
        -k(190)*n(idx_NAj)*n(idx_E) &
        +k(195)*n(idx_NA)*n(idx_MGj) &
        +k(230)*n(idx_Hj)*n(idx_NA)

    !HCO+
    !HCO+
    dn(idx_HCOj) = &
        -k(52)*n(idx_C)*n(idx_HCOj) &
        -k(71)*n(idx_SI)*n(idx_HCOj) &
        -k(138)*n(idx_HCOj)*n(idx_FE) &
        -k(164)*n(idx_MG)*n(idx_HCOj) &
        +k(200)*n(idx_CH)*n(idx_O) &
        -k(249)*n(idx_HCOj)*n(idx_E)

    !CH+
    !CH+
    dn(idx_CHj) = &
        +k(52)*n(idx_C)*n(idx_HCOj) &
        -k(194)*n(idx_H)*n(idx_CHj) &
        +k(206)*n(idx_H2j)*n(idx_C) &
        -k(209)*n(idx_CHj)*n(idx_E)

    !O+
    !O+
    dn(idx_Oj) = &
        -k(55)*n(idx_H2)*n(idx_Oj) &
        -k(68)*n(idx_Oj)*n(idx_FE) &
        -k(73)*n(idx_Hk)*n(idx_Oj) &
        -k(74)*n(idx_Oj)*n(idx_E) &
        -k(182)*n(idx_H)*n(idx_Oj) &
        +k(184)*n(idx_Hj)*n(idx_O) &
        +k(239)*n(idx_O)

    !MG+
    !MG+
    dn(idx_MGj) = &
        +k(59)*n(idx_MG)*n(idx_Sj) &
        -k(63)*n(idx_Ok)*n(idx_MGj) &
        -k(64)*n(idx_MGj)*n(idx_E) &
        +k(111)*n(idx_Cj)*n(idx_MG) &
        +k(126)*n(idx_MG)*n(idx_SIOj) &
        +k(130)*n(idx_MG)*n(idx_SIj) &
        +k(144)*n(idx_Hj)*n(idx_MG) &
        +k(164)*n(idx_MG)*n(idx_HCOj) &
        -k(195)*n(idx_NA)*n(idx_MGj) &
        -k(214)*n(idx_Hk)*n(idx_MGj)

    !SIO+
    !SIO+
    dn(idx_SIOj) = &
        -k(79)*n(idx_SIOj)*n(idx_FE) &
        +k(81)*n(idx_OH)*n(idx_SIj) &
        -k(97)*n(idx_O)*n(idx_SIOj) &
        -k(115)*n(idx_SIOj)*n(idx_E) &
        -k(126)*n(idx_MG)*n(idx_SIOj) &
        +k(140)*n(idx_Hj)*n(idx_SIO) &
        -k(156)*n(idx_N)*n(idx_SIOj) &
        -k(202)*n(idx_C)*n(idx_SIOj)

    !P+
    !P+
    dn(idx_Pj) = &
        +k(84)*n(idx_Hj)*n(idx_P) &
        -k(139)*n(idx_Pj)*n(idx_E) &
        -k(143)*n(idx_SI)*n(idx_Pj)

    !SIF+
    !SIF+
    dn(idx_SIFj) = &
        +k(100)*n(idx_HF)*n(idx_SIj) &
        -k(109)*n(idx_SIFj)*n(idx_E)

    !C+
    !C+
    dn(idx_Cj) = &
        -k(111)*n(idx_Cj)*n(idx_MG) &
        -k(121)*n(idx_Cj)*n(idx_SI) &
        -k(124)*n(idx_Cj)*n(idx_FE) &
        +k(194)*n(idx_H)*n(idx_CHj) &
        +k(232)*n(idx_C) &
        -k(243)*n(idx_Cj)*n(idx_E)

    !N+
    !N+
    dn(idx_Nj) = &
        +k(234)*n(idx_N) &
        -k(244)*n(idx_Nj)*n(idx_E)

    !CO+
    !CO+
    dn(idx_COj) = &
        +k(235)*n(idx_CO) &
        -k(245)*n(idx_COj)*n(idx_E)

    !F+
    !F+
    dn(idx_Fj) = &
        +k(250)*n(idx_HEj)*n(idx_HF) &
        -k(251)*n(idx_H2)*n(idx_Fj)

    !CR
    !CR
    dn(idx_CR) = &
        0.d0

    !g

    !g
    dn(idx_g) = 0.d0

    !Tgas

    !Tgas
    dn(idx_Tgas) = 0.d0

    !dummy

    !dummy
    dn(idx_dummy) = 0.d0

    krome_gamma = gamma_index(n(:))

    dn(idx_Tgas) = (heating(n(:), Tgas, k(:), nH2dust) &
        - cooling(n(:), Tgas)  ) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

    last_coe(:) = k(:)

  end subroutine fex

  !***************************
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_gadiab
    use krome_getphys
    implicit none
    integer::neq, j, ian, jan, r1, r2, p1, p2, p3, i
    real*8::tt, n(neq), pdj(neq), dr1, dr2, kk,k(nrea),Tgas
    real*8::nn(neq),dn0,dn1,dnn,nH2dust,dn(neq),krome_gamma

    nH2dust = 0.d0
    Tgas = n(idx_Tgas)

    krome_gamma = gamma_index(n(:))

    k(:) = last_coe(:) !get rate coefficients

    if(j==1) then
    elseif(j==1) then
      pdj(1) =  &
          -k(187)*n(idx_Sj)  &
          +2.d0*k(220)*n(idx_H)  &
          -k(157)*n(idx_C)  &
          -k(244)*n(idx_Nj)  &
          -k(25)*n(idx_HEHj)  &
          -k(4)*n(idx_HEj)  &
          -k(220)*n(idx_H)  &
          -k(201)*n(idx_SIHj)  &
          -k(245)*n(idx_COj)  &
          -k(139)*n(idx_Pj)  &
          -k(221)*n(idx_HE)  &
          -k(248)*n(idx_HSj)  &
          -k(29)*n(idx_Hj)  &
          +2.d0*k(221)*n(idx_HE)  &
          -k(115)*n(idx_SIOj)  &
          -k(142)*n(idx_FEj)  &
          -k(190)*n(idx_NAj)  &
          -k(102)*n(idx_H2j)  &
          -k(13)*n(idx_O)  &
          -k(107)*n(idx_H2)  &
          -k(3)*n(idx_NHj)  &
          -k(249)*n(idx_HCOj)  &
          -k(223)*n(idx_H2)  &
          -k(197)*n(idx_SIj)  &
          -k(224)*n(idx_Hk)  &
          +2.d0*k(224)*n(idx_Hk)  &
          -k(74)*n(idx_Oj)  &
          -k(64)*n(idx_MGj)  &
          -k(109)*n(idx_SIFj)  &
          -k(57)*n(idx_H)  &
          -k(87)*n(idx_OHj)  &
          -k(243)*n(idx_Cj)  &
          -k(137)*n(idx_S)  &
          +k(107)*n(idx_H2)  &
          -k(209)*n(idx_CHj)
      pdj(2) =  &
          +k(13)*n(idx_O)
      pdj(3) =  &
          -k(224)*n(idx_Hk)  &
          +k(57)*n(idx_H)  &
          +k(223)*n(idx_H2)
      pdj(4) =  &
          +k(137)*n(idx_S)
      pdj(5) =  &
          +k(157)*n(idx_C)
      pdj(6) =  &
          +k(29)*n(idx_Hj)  &
          +k(223)*n(idx_H2)  &
          +k(201)*n(idx_SIHj)  &
          -k(57)*n(idx_H)  &
          +k(209)*n(idx_CHj)  &
          +k(248)*n(idx_HSj)  &
          +k(25)*n(idx_HEHj)  &
          +k(224)*n(idx_Hk)  &
          +2.d0*k(102)*n(idx_H2j)  &
          -k(220)*n(idx_H)  &
          +k(3)*n(idx_NHj)  &
          +k(87)*n(idx_OHj)  &
          +2.d0*k(107)*n(idx_H2)  &
          +k(249)*n(idx_HCOj)
      pdj(7) =  &
          +k(25)*n(idx_HEHj)  &
          +k(4)*n(idx_HEj)  &
          -k(221)*n(idx_HE)
      pdj(8) =  &
          +k(245)*n(idx_COj)  &
          +k(209)*n(idx_CHj)  &
          -k(157)*n(idx_C)  &
          +k(243)*n(idx_Cj)
      pdj(10) =  &
          +k(249)*n(idx_HCOj)
      pdj(11) =  &
          +k(244)*n(idx_Nj)  &
          +k(3)*n(idx_NHj)
      pdj(13) =  &
          -k(137)*n(idx_S)  &
          +k(187)*n(idx_Sj)  &
          +k(248)*n(idx_HSj)
      pdj(15) =  &
          +k(74)*n(idx_Oj)  &
          +k(245)*n(idx_COj)  &
          -k(13)*n(idx_O)  &
          +k(115)*n(idx_SIOj)  &
          +k(87)*n(idx_OHj)
      pdj(16) =  &
          -k(107)*n(idx_H2)  &
          -k(223)*n(idx_H2)
      pdj(17) =  &
          +k(109)*n(idx_SIFj)  &
          +k(197)*n(idx_SIj)  &
          +k(201)*n(idx_SIHj)  &
          +k(115)*n(idx_SIOj)
      pdj(22) =  &
          +k(142)*n(idx_FEj)
      pdj(26) =  &
          +k(190)*n(idx_NAj)
      pdj(27) =  &
          +k(109)*n(idx_SIFj)
      pdj(41) =  &
          +k(64)*n(idx_MGj)
      pdj(42) =  &
          +k(139)*n(idx_Pj)
      pdj(48) =  &
          +k(221)*n(idx_HE)  &
          -k(4)*n(idx_HEj)
      pdj(49) =  &
          +k(220)*n(idx_H)  &
          -k(29)*n(idx_Hj)
      pdj(50) =  &
          -k(3)*n(idx_NHj)
      pdj(51) =  &
          -k(248)*n(idx_HSj)
      pdj(52) =  &
          -k(187)*n(idx_Sj)
      pdj(53) =  &
          -k(197)*n(idx_SIj)
      pdj(54) =  &
          -k(87)*n(idx_OHj)
      pdj(55) =  &
          -k(25)*n(idx_HEHj)
      pdj(56) =  &
          -k(102)*n(idx_H2j)
      pdj(57) =  &
          -k(142)*n(idx_FEj)
      pdj(58) =  &
          -k(201)*n(idx_SIHj)
      pdj(59) =  &
          -k(190)*n(idx_NAj)
      pdj(60) =  &
          -k(249)*n(idx_HCOj)
      pdj(61) =  &
          -k(209)*n(idx_CHj)
      pdj(62) =  &
          -k(74)*n(idx_Oj)
      pdj(63) =  &
          -k(64)*n(idx_MGj)
      pdj(64) =  &
          -k(115)*n(idx_SIOj)
      pdj(65) =  &
          -k(139)*n(idx_Pj)
      pdj(66) =  &
          -k(109)*n(idx_SIFj)
      pdj(67) =  &
          -k(243)*n(idx_Cj)
      pdj(68) =  &
          -k(244)*n(idx_Nj)
      pdj(69) =  &
          -k(245)*n(idx_COj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(1)*1d-3
      if(dnn>0.d0) then
        nn(1) = n(1) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==2) then
      pdj(1) =  &
          +k(50)*n(idx_H)
      pdj(2) =  &
          -k(210)*n(idx_FEj)  &
          -k(63)*n(idx_MGj)  &
          -k(199)*n(idx_Hj)  &
          -k(50)*n(idx_H)
      pdj(6) =  &
          +k(199)*n(idx_Hj)  &
          -k(50)*n(idx_H)
      pdj(15) =  &
          +k(199)*n(idx_Hj)  &
          +k(63)*n(idx_MGj)  &
          +k(210)*n(idx_FEj)
      pdj(18) =  &
          +k(50)*n(idx_H)
      pdj(22) =  &
          +k(210)*n(idx_FEj)
      pdj(41) =  &
          +k(63)*n(idx_MGj)
      pdj(49) =  &
          -k(199)*n(idx_Hj)
      pdj(57) =  &
          -k(210)*n(idx_FEj)
      pdj(63) =  &
          -k(63)*n(idx_MGj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(2)*1d-3
      if(dnn>0.d0) then
        nn(2) = n(2) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==3) then
      pdj(1) =  &
          -k(224)*n(idx_E)  &
          +k(226)*n(idx_Hj)  &
          +k(192)*n(idx_H)  &
          +k(129)*n(idx_C)  &
          +k(225)*n(idx_H)  &
          +2.d0*k(224)*n(idx_E)  &
          +k(69)*n(idx_N)  &
          +k(203)*n(idx_O)
      pdj(3) =  &
          -k(224)*n(idx_E)  &
          -k(203)*n(idx_O)  &
          -k(65)*n(idx_Hj)  &
          -k(69)*n(idx_N)  &
          -k(66)*n(idx_Sj)  &
          -k(214)*n(idx_MGj)  &
          -k(225)*n(idx_H)  &
          -k(129)*n(idx_C)  &
          -k(24)*n(idx_FEj)  &
          -k(145)*n(idx_SIj)  &
          -k(192)*n(idx_H)  &
          -k(226)*n(idx_Hj)  &
          -k(22)*n(idx_NAj)  &
          -k(73)*n(idx_Oj)
      pdj(6) =  &
          +2.d0*k(65)*n(idx_Hj)  &
          +k(24)*n(idx_FEj)  &
          +k(73)*n(idx_Oj)  &
          +k(145)*n(idx_SIj)  &
          +k(66)*n(idx_Sj)  &
          -k(225)*n(idx_H)  &
          +2.d0*k(225)*n(idx_H)  &
          -k(192)*n(idx_H)  &
          +k(224)*n(idx_E)  &
          +k(214)*n(idx_MGj)  &
          +k(22)*n(idx_NAj)
      pdj(8) =  &
          -k(129)*n(idx_C)
      pdj(11) =  &
          -k(69)*n(idx_N)
      pdj(13) =  &
          +k(66)*n(idx_Sj)
      pdj(15) =  &
          +k(73)*n(idx_Oj)  &
          -k(203)*n(idx_O)
      pdj(16) =  &
          +k(192)*n(idx_H)
      pdj(17) =  &
          +k(145)*n(idx_SIj)
      pdj(18) =  &
          +k(203)*n(idx_O)
      pdj(22) =  &
          +k(24)*n(idx_FEj)
      pdj(26) =  &
          +k(22)*n(idx_NAj)
      pdj(29) =  &
          +k(129)*n(idx_C)
      pdj(34) =  &
          +k(69)*n(idx_N)
      pdj(41) =  &
          +k(214)*n(idx_MGj)
      pdj(49) =  &
          -k(226)*n(idx_Hj)  &
          -k(65)*n(idx_Hj)
      pdj(52) =  &
          -k(66)*n(idx_Sj)
      pdj(53) =  &
          -k(145)*n(idx_SIj)
      pdj(56) =  &
          +k(226)*n(idx_Hj)
      pdj(57) =  &
          -k(24)*n(idx_FEj)
      pdj(59) =  &
          -k(22)*n(idx_NAj)
      pdj(62) =  &
          -k(73)*n(idx_Oj)
      pdj(63) =  &
          -k(214)*n(idx_MGj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(3)*1d-3
      if(dnn>0.d0) then
        nn(3) = n(3) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==4) then
      pdj(1) =  &
          +k(76)*n(idx_H)
      pdj(4) =  &
          -k(76)*n(idx_H)
      pdj(6) =  &
          -k(76)*n(idx_H)
      pdj(19) =  &
          +k(76)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(4)*1d-3
      if(dnn>0.d0) then
        nn(4) = n(4) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==5) then
      pdj(1) =  &
          +k(170)*n(idx_H)
      pdj(5) =  &
          -k(133)*n(idx_Hj)  &
          -k(170)*n(idx_H)
      pdj(6) =  &
          +k(133)*n(idx_Hj)  &
          -k(170)*n(idx_H)
      pdj(8) =  &
          +k(133)*n(idx_Hj)
      pdj(29) =  &
          +k(170)*n(idx_H)
      pdj(49) =  &
          -k(133)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(5)*1d-3
      if(dnn>0.d0) then
        nn(5) = n(5) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==6) then
      pdj(1) =  &
          +k(233)  &
          -k(220)*n(idx_E)  &
          +k(225)*n(idx_Hk)  &
          +k(170)*n(idx_Ck)  &
          +k(192)*n(idx_Hk)  &
          +k(50)*n(idx_Ok)  &
          -k(57)*n(idx_E)  &
          +2.d0*k(220)*n(idx_E)  &
          +k(76)*n(idx_Sk)
      pdj(2) =  &
          -k(50)*n(idx_Ok)
      pdj(3) =  &
          +k(57)*n(idx_E)  &
          -k(192)*n(idx_Hk)  &
          -k(225)*n(idx_Hk)
      pdj(4) =  &
          -k(76)*n(idx_Sk)
      pdj(5) =  &
          -k(170)*n(idx_Ck)
      pdj(6) =  &
          -k(220)*n(idx_E)  &
          -k(50)*n(idx_Ok)  &
          -k(28)*n(idx_O2)  &
          -k(135)*n(idx_OCS)  &
          +3.d0*k(33)*n(idx_H2)  &
          -k(30)*n(idx_NO)  &
          -k(11)*n(idx_HEHj)  &
          -k(76)*n(idx_Sk)  &
          -k(18)*n(idx_S2)  &
          -k(56)*n(idx_OCN)  &
          -k(110)*n(idx_SIHj)  &
          -k(163)*n(idx_H2S)  &
          -k(176)*n(idx_CO2)  &
          -k(93)*n(idx_O2)  &
          -k(36)*n(idx_CO)  &
          +2.d0*k(67)*n(idx_OH)  &
          -k(62)*n(idx_HS)  &
          -k(53)*n(idx_NH2)  &
          -k(217)*n(idx_NS)  &
          -k(194)*n(idx_CHj)  &
          -k(83)*n(idx_O)  &
          -k(38)*n(idx_CH2)  &
          -k(208)*n(idx_OH)  &
          -k(40)*n(idx_NH)  &
          -k(185)*n(idx_SO)  &
          -k(67)*n(idx_OH)  &
          -k(6)*n(idx_HSj)  &
          -k(70)*n(idx_NO)  &
          -k(153)*n(idx_C2)  &
          +k(28)*n(idx_O2)  &
          -k(225)*n(idx_Hk)  &
          -k(1)*n(idx_HEj)  &
          -k(150)*n(idx_OCN)  &
          +2.d0*k(108)*n(idx_CH)  &
          +3.d0*k(227)*n(idx_H)*n(idx_H)  &
          -k(196)*n(idx_H2j)  &
          -k(57)*n(idx_E)  &
          -9.d0*k(227)*n(idx_H)*n(idx_H)  &
          -k(95)*n(idx_HCO)  &
          -k(16)*n(idx_NS)  &
          -k(108)*n(idx_CH)  &
          -k(170)*n(idx_Ck)  &
          -4.d0*k(228)*n(idx_H2)*n(idx_H)  &
          -k(182)*n(idx_Oj)  &
          -k(88)*n(idx_C)  &
          +2.d0*k(225)*n(idx_Hk)  &
          -k(92)*n(idx_CH)  &
          -k(61)*n(idx_Hj)  &
          -k(43)*n(idx_HCN)  &
          -k(233)  &
          -k(33)*n(idx_H2)  &
          -k(103)*n(idx_H2O)  &
          -k(21)*n(idx_SO)  &
          -k(192)*n(idx_Hk)  &
          -4.d0*k(229)*n(idx_H)*n(idx_HE)  &
          -k(154)*n(idx_OCN)  &
          -k(15)*n(idx_SIj)  &
          -k(146)*n(idx_H2O)  &
          +2.d0*k(103)*n(idx_H2O)
      pdj(7) =  &
          +2.d0*k(229)*n(idx_H)*n(idx_HE)  &
          +k(11)*n(idx_HEHj)  &
          +k(1)*n(idx_HEj)  &
          -2.d0*k(229)*n(idx_H)*n(idx_HE)
      pdj(8) =  &
          -k(88)*n(idx_C)  &
          +k(92)*n(idx_CH)  &
          +k(108)*n(idx_CH)  &
          +k(36)*n(idx_CO)  &
          +k(153)*n(idx_C2)
      pdj(9) =  &
          -k(30)*n(idx_NO)  &
          -k(70)*n(idx_NO)
      pdj(10) =  &
          +k(135)*n(idx_OCS)  &
          +k(154)*n(idx_OCN)  &
          -k(36)*n(idx_CO)  &
          +k(95)*n(idx_HCO)  &
          +k(176)*n(idx_CO2)
      pdj(11) =  &
          +k(40)*n(idx_NH)  &
          +k(30)*n(idx_NO)  &
          +k(16)*n(idx_NS)
      pdj(12) =  &
          -k(93)*n(idx_O2)  &
          -k(28)*n(idx_O2)
      pdj(13) =  &
          +k(21)*n(idx_SO)  &
          +k(217)*n(idx_NS)  &
          +k(18)*n(idx_S2)  &
          +k(62)*n(idx_HS)
      pdj(14) =  &
          -k(21)*n(idx_SO)  &
          -k(185)*n(idx_SO)
      pdj(15) =  &
          +k(150)*n(idx_OCN)  &
          +k(67)*n(idx_OH)  &
          +k(182)*n(idx_Oj)  &
          +2.d0*k(28)*n(idx_O2)  &
          +k(93)*n(idx_O2)  &
          +k(185)*n(idx_SO)  &
          +k(70)*n(idx_NO)  &
          -k(83)*n(idx_O)  &
          +k(208)*n(idx_OH)
      pdj(16) =  &
          +4.d0*k(228)*n(idx_H2)*n(idx_H)  &
          +3.d0*k(227)*n(idx_H)*n(idx_H)  &
          +k(92)*n(idx_CH)  &
          +k(194)*n(idx_CHj)  &
          -2.d0*k(228)*n(idx_H2)*n(idx_H)  &
          +k(43)*n(idx_HCN)  &
          +k(40)*n(idx_NH)  &
          +k(163)*n(idx_H2S)  &
          +k(146)*n(idx_H2O)  &
          +k(6)*n(idx_HSj)  &
          -k(33)*n(idx_H2)  &
          +2.d0*k(229)*n(idx_H)*n(idx_HE)  &
          +k(192)*n(idx_Hk)  &
          +k(53)*n(idx_NH2)  &
          +k(62)*n(idx_HS)  &
          +k(196)*n(idx_H2j)  &
          +k(110)*n(idx_SIHj)  &
          +k(208)*n(idx_OH)  &
          +k(38)*n(idx_CH2)  &
          +k(95)*n(idx_HCO)
      pdj(18) =  &
          +k(146)*n(idx_H2O)  &
          +k(56)*n(idx_OCN)  &
          +k(30)*n(idx_NO)  &
          +k(93)*n(idx_O2)  &
          -k(67)*n(idx_OH)  &
          +k(83)*n(idx_O)  &
          -k(208)*n(idx_OH)  &
          +k(36)*n(idx_CO)  &
          +k(176)*n(idx_CO2)  &
          +k(21)*n(idx_SO)  &
          +k(103)*n(idx_H2O)  &
          +k(50)*n(idx_Ok)
      pdj(19) =  &
          +k(135)*n(idx_OCS)  &
          +k(163)*n(idx_H2S)  &
          +k(18)*n(idx_S2)  &
          -k(62)*n(idx_HS)  &
          +k(16)*n(idx_NS)  &
          +k(76)*n(idx_Sk)  &
          +k(185)*n(idx_SO)
      pdj(20) =  &
          -k(16)*n(idx_NS)  &
          -k(217)*n(idx_NS)
      pdj(21) =  &
          -k(163)*n(idx_H2S)
      pdj(24) =  &
          +k(56)*n(idx_OCN)  &
          +k(43)*n(idx_HCN)
      pdj(25) =  &
          -k(18)*n(idx_S2)
      pdj(29) =  &
          -k(108)*n(idx_CH)  &
          +k(170)*n(idx_Ck)  &
          -k(92)*n(idx_CH)  &
          +k(153)*n(idx_C2)  &
          +k(88)*n(idx_C)  &
          +k(38)*n(idx_CH2)
      pdj(31) =  &
          -k(153)*n(idx_C2)
      pdj(33) =  &
          -k(38)*n(idx_CH2)
      pdj(34) =  &
          +k(154)*n(idx_OCN)  &
          -k(40)*n(idx_NH)  &
          +k(70)*n(idx_NO)  &
          +k(53)*n(idx_NH2)  &
          +k(217)*n(idx_NS)
      pdj(35) =  &
          +k(150)*n(idx_OCN)  &
          -k(43)*n(idx_HCN)
      pdj(36) =  &
          -k(176)*n(idx_CO2)
      pdj(39) =  &
          -k(53)*n(idx_NH2)
      pdj(40) =  &
          -k(150)*n(idx_OCN)  &
          -k(154)*n(idx_OCN)  &
          -k(56)*n(idx_OCN)
      pdj(43) =  &
          -k(95)*n(idx_HCO)
      pdj(44) =  &
          -k(103)*n(idx_H2O)  &
          -k(146)*n(idx_H2O)
      pdj(45) =  &
          -k(135)*n(idx_OCS)
      pdj(48) =  &
          -k(1)*n(idx_HEj)
      pdj(49) =  &
          +k(233)  &
          -k(61)*n(idx_Hj)  &
          +k(1)*n(idx_HEj)  &
          +k(182)*n(idx_Oj)  &
          +k(196)*n(idx_H2j)  &
          +k(220)*n(idx_E)
      pdj(51) =  &
          -k(6)*n(idx_HSj)
      pdj(52) =  &
          +k(6)*n(idx_HSj)
      pdj(53) =  &
          +k(110)*n(idx_SIHj)  &
          -k(15)*n(idx_SIj)
      pdj(55) =  &
          -k(11)*n(idx_HEHj)
      pdj(56) =  &
          -k(196)*n(idx_H2j)  &
          +k(11)*n(idx_HEHj)  &
          +k(61)*n(idx_Hj)
      pdj(58) =  &
          +k(15)*n(idx_SIj)  &
          -k(110)*n(idx_SIHj)
      pdj(61) =  &
          -k(194)*n(idx_CHj)
      pdj(62) =  &
          -k(182)*n(idx_Oj)
      pdj(67) =  &
          +k(194)*n(idx_CHj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(6)*1d-3
      if(dnn>0.d0) then
        nn(6) = n(6) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==7) then
      pdj(1) =  &
          -k(221)*n(idx_E)  &
          +2.d0*k(221)*n(idx_E)  &
          +k(238)
      pdj(6) =  &
          -2.d0*k(229)*n(idx_H)*n(idx_H)  &
          +k(32)*n(idx_H2j)
      pdj(7) =  &
          -k(32)*n(idx_H2j)  &
          -k(221)*n(idx_E)  &
          +k(229)*n(idx_H)*n(idx_H)  &
          -k(238)  &
          -k(229)*n(idx_H)*n(idx_H)
      pdj(16) =  &
          +k(229)*n(idx_H)*n(idx_H)
      pdj(48) =  &
          +k(221)*n(idx_E)  &
          +k(238)
      pdj(55) =  &
          +k(32)*n(idx_H2j)
      pdj(56) =  &
          -k(32)*n(idx_H2j)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(7)*1d-3
      if(dnn>0.d0) then
        nn(7) = n(7) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==8) then
      pdj(1) =  &
          -k(157)*n(idx_E)  &
          +k(129)*n(idx_Hk)  &
          +k(232)
      pdj(3) =  &
          -k(129)*n(idx_Hk)
      pdj(5) =  &
          +k(157)*n(idx_E)
      pdj(6) =  &
          +k(90)*n(idx_HS)  &
          +k(152)*n(idx_CH)  &
          +k(54)*n(idx_H2)  &
          +k(206)*n(idx_H2j)  &
          -k(88)*n(idx_H)  &
          +k(162)*n(idx_NH)  &
          +k(204)*n(idx_OH)
      pdj(8) =  &
          -k(105)*n(idx_NH)  &
          -k(157)*n(idx_E)  &
          -k(206)*n(idx_H2j)  &
          -k(129)*n(idx_Hk)  &
          -k(232)  &
          -k(75)*n(idx_H2)  &
          -k(204)*n(idx_OH)  &
          -k(162)*n(idx_NH)  &
          -k(54)*n(idx_H2)  &
          -k(52)*n(idx_HCOj)  &
          -k(152)*n(idx_CH)  &
          -k(34)*n(idx_CS)  &
          -k(198)*n(idx_NS)  &
          -k(47)*n(idx_CN)  &
          -k(189)*n(idx_SO2)  &
          -k(202)*n(idx_SIOj)  &
          -k(80)*n(idx_SO)  &
          -k(106)*n(idx_N2)  &
          -k(193)*n(idx_SO)  &
          -k(173)*n(idx_N)  &
          -k(88)*n(idx_H)  &
          -k(101)*n(idx_NO)  &
          -4.d0*k(166)*n(idx_C)  &
          -k(169)*n(idx_OH)  &
          -k(90)*n(idx_HS)  &
          -k(112)*n(idx_O)  &
          -k(171)*n(idx_HS)  &
          -k(2)*n(idx_NO)  &
          -k(216)*n(idx_O2)  &
          -k(46)*n(idx_CO)  &
          -k(207)*n(idx_S)
      pdj(9) =  &
          -k(101)*n(idx_NO)  &
          -k(2)*n(idx_NO)
      pdj(10) =  &
          +k(2)*n(idx_NO)  &
          +k(202)*n(idx_SIOj)  &
          +k(52)*n(idx_HCOj)  &
          +k(112)*n(idx_O)  &
          +k(193)*n(idx_SO)  &
          -k(46)*n(idx_CO)  &
          +k(189)*n(idx_SO2)  &
          +k(216)*n(idx_O2)  &
          +k(204)*n(idx_OH)
      pdj(11) =  &
          +k(106)*n(idx_N2)  &
          -k(173)*n(idx_N)  &
          +k(2)*n(idx_NO)  &
          +k(105)*n(idx_NH)  &
          +k(47)*n(idx_CN)
      pdj(12) =  &
          -k(216)*n(idx_O2)
      pdj(13) =  &
          +k(171)*n(idx_HS)  &
          +k(193)*n(idx_SO)  &
          +k(34)*n(idx_CS)  &
          -k(207)*n(idx_S)  &
          +k(198)*n(idx_NS)
      pdj(14) =  &
          -k(193)*n(idx_SO)  &
          -k(80)*n(idx_SO)  &
          +k(189)*n(idx_SO2)
      pdj(15) =  &
          +k(46)*n(idx_CO)  &
          -k(112)*n(idx_O)  &
          +k(216)*n(idx_O2)  &
          +k(101)*n(idx_NO)  &
          +k(80)*n(idx_SO)  &
          +k(169)*n(idx_OH)
      pdj(16) =  &
          -k(75)*n(idx_H2)  &
          -k(54)*n(idx_H2)
      pdj(18) =  &
          -k(169)*n(idx_OH)  &
          -k(204)*n(idx_OH)
      pdj(19) =  &
          -k(171)*n(idx_HS)  &
          -k(90)*n(idx_HS)
      pdj(20) =  &
          -k(198)*n(idx_NS)
      pdj(23) =  &
          +k(207)*n(idx_S)  &
          +k(90)*n(idx_HS)  &
          +k(80)*n(idx_SO)  &
          -k(34)*n(idx_CS)
      pdj(24) =  &
          +k(198)*n(idx_NS)  &
          -k(47)*n(idx_CN)  &
          +k(101)*n(idx_NO)  &
          +k(173)*n(idx_N)  &
          +k(106)*n(idx_N2)  &
          +k(162)*n(idx_NH)
      pdj(29) =  &
          -k(152)*n(idx_CH)  &
          +k(105)*n(idx_NH)  &
          +k(54)*n(idx_H2)  &
          +k(88)*n(idx_H)  &
          +k(171)*n(idx_HS)  &
          +k(129)*n(idx_Hk)  &
          +k(169)*n(idx_OH)
      pdj(30) =  &
          -k(189)*n(idx_SO2)
      pdj(31) =  &
          +k(152)*n(idx_CH)  &
          +2.d0*k(166)*n(idx_C)  &
          +k(46)*n(idx_CO)  &
          +k(34)*n(idx_CS)  &
          +k(47)*n(idx_CN)
      pdj(32) =  &
          -k(106)*n(idx_N2)
      pdj(33) =  &
          +k(75)*n(idx_H2)
      pdj(34) =  &
          -k(162)*n(idx_NH)  &
          -k(105)*n(idx_NH)
      pdj(53) =  &
          +k(202)*n(idx_SIOj)
      pdj(56) =  &
          -k(206)*n(idx_H2j)
      pdj(60) =  &
          -k(52)*n(idx_HCOj)
      pdj(61) =  &
          +k(206)*n(idx_H2j)  &
          +k(52)*n(idx_HCOj)
      pdj(64) =  &
          -k(202)*n(idx_SIOj)
      pdj(67) =  &
          +k(232)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(8)*1d-3
      if(dnn>0.d0) then
        nn(8) = n(8) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==9) then
      pdj(6) =  &
          -k(30)*n(idx_H)  &
          -k(70)*n(idx_H)
      pdj(8) =  &
          -k(101)*n(idx_C)  &
          -k(2)*n(idx_C)
      pdj(9) =  &
          -k(30)*n(idx_H)  &
          -k(70)*n(idx_H)  &
          -k(123)*n(idx_SI)  &
          -k(101)*n(idx_C)  &
          -k(2)*n(idx_C)  &
          -k(35)*n(idx_N)
      pdj(10) =  &
          +k(2)*n(idx_C)
      pdj(11) =  &
          +k(123)*n(idx_SI)  &
          +k(2)*n(idx_C)  &
          +k(30)*n(idx_H)  &
          -k(35)*n(idx_N)
      pdj(15) =  &
          +k(35)*n(idx_N)  &
          +k(70)*n(idx_H)  &
          +k(101)*n(idx_C)
      pdj(17) =  &
          -k(123)*n(idx_SI)
      pdj(18) =  &
          +k(30)*n(idx_H)
      pdj(24) =  &
          +k(101)*n(idx_C)
      pdj(32) =  &
          +k(35)*n(idx_N)
      pdj(34) =  &
          +k(70)*n(idx_H)
      pdj(37) =  &
          +k(123)*n(idx_SI)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(9)*1d-3
      if(dnn>0.d0) then
        nn(9) = n(9) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==10) then
      pdj(1) =  &
          +k(235)
      pdj(6) =  &
          +k(44)*n(idx_OH)  &
          -k(36)*n(idx_H)
      pdj(8) =  &
          +k(36)*n(idx_H)  &
          +k(242)  &
          +k(89)*n(idx_SI)  &
          -k(46)*n(idx_C)
      pdj(10) =  &
          -k(46)*n(idx_C)  &
          -k(44)*n(idx_OH)  &
          -k(36)*n(idx_H)  &
          -k(89)*n(idx_SI)  &
          -k(235)  &
          -k(242)
      pdj(15) =  &
          +k(46)*n(idx_C)  &
          +k(242)
      pdj(17) =  &
          -k(89)*n(idx_SI)
      pdj(18) =  &
          +k(36)*n(idx_H)  &
          -k(44)*n(idx_OH)
      pdj(31) =  &
          +k(46)*n(idx_C)
      pdj(36) =  &
          +k(44)*n(idx_OH)
      pdj(37) =  &
          +k(89)*n(idx_SI)
      pdj(69) =  &
          +k(235)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(10)*1d-3
      if(dnn>0.d0) then
        nn(10) = n(10) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==11) then
      pdj(1) =  &
          +k(234)  &
          +k(69)*n(idx_Hk)
      pdj(3) =  &
          -k(69)*n(idx_Hk)
      pdj(6) =  &
          +k(48)*n(idx_NH)  &
          +k(117)*n(idx_H2)  &
          +k(183)*n(idx_OH)  &
          +k(165)*n(idx_CH)  &
          +k(9)*n(idx_HS)
      pdj(8) =  &
          +k(151)*n(idx_C2)  &
          -k(173)*n(idx_C)  &
          +k(128)*n(idx_CN)  &
          +k(131)*n(idx_CH)
      pdj(9) =  &
          +k(183)*n(idx_OH)  &
          +k(255)*n(idx_PO)  &
          +k(60)*n(idx_O2)  &
          +k(156)*n(idx_SIOj)  &
          +k(175)*n(idx_CO2)  &
          +k(211)*n(idx_SO)  &
          -k(35)*n(idx_NO)
      pdj(10) =  &
          +k(175)*n(idx_CO2)
      pdj(11) =  &
          -k(173)*n(idx_C)  &
          -k(131)*n(idx_CH)  &
          -k(69)*n(idx_Hk)  &
          -k(165)*n(idx_CH)  &
          -k(17)*n(idx_CS)  &
          -k(35)*n(idx_NO)  &
          -k(60)*n(idx_O2)  &
          -k(128)*n(idx_CN)  &
          -k(234)  &
          -k(183)*n(idx_OH)  &
          -k(156)*n(idx_SIOj)  &
          -k(94)*n(idx_SO)  &
          -k(252)*n(idx_PN)  &
          -k(155)*n(idx_OH)  &
          -k(211)*n(idx_SO)  &
          -k(9)*n(idx_HS)  &
          -k(175)*n(idx_CO2)  &
          -k(255)*n(idx_PO)  &
          -k(179)*n(idx_HS)  &
          -k(253)*n(idx_PO)  &
          -k(117)*n(idx_H2)  &
          -k(48)*n(idx_NH)  &
          -k(151)*n(idx_C2)
      pdj(12) =  &
          -k(60)*n(idx_O2)
      pdj(13) =  &
          +k(17)*n(idx_CS)  &
          +k(179)*n(idx_HS)  &
          +k(211)*n(idx_SO)
      pdj(14) =  &
          -k(94)*n(idx_SO)  &
          -k(211)*n(idx_SO)
      pdj(15) =  &
          +k(35)*n(idx_NO)  &
          +k(155)*n(idx_OH)  &
          +k(253)*n(idx_PO)  &
          +k(94)*n(idx_SO)  &
          +k(60)*n(idx_O2)
      pdj(16) =  &
          -k(117)*n(idx_H2)
      pdj(18) =  &
          -k(183)*n(idx_OH)  &
          -k(155)*n(idx_OH)
      pdj(19) =  &
          -k(9)*n(idx_HS)  &
          -k(179)*n(idx_HS)
      pdj(20) =  &
          +k(94)*n(idx_SO)  &
          +k(9)*n(idx_HS)
      pdj(23) =  &
          -k(17)*n(idx_CS)
      pdj(24) =  &
          +k(17)*n(idx_CS)  &
          -k(128)*n(idx_CN)  &
          +k(151)*n(idx_C2)  &
          +k(173)*n(idx_C)  &
          +k(165)*n(idx_CH)
      pdj(29) =  &
          -k(165)*n(idx_CH)  &
          -k(131)*n(idx_CH)
      pdj(31) =  &
          -k(151)*n(idx_C2)
      pdj(32) =  &
          +k(48)*n(idx_NH)  &
          +k(35)*n(idx_NO)  &
          +k(128)*n(idx_CN)  &
          +k(252)*n(idx_PN)
      pdj(34) =  &
          +k(155)*n(idx_OH)  &
          +k(117)*n(idx_H2)  &
          +k(69)*n(idx_Hk)  &
          +k(179)*n(idx_HS)  &
          -k(48)*n(idx_NH)  &
          +k(131)*n(idx_CH)
      pdj(36) =  &
          -k(175)*n(idx_CO2)
      pdj(42) =  &
          +k(252)*n(idx_PN)  &
          +k(255)*n(idx_PO)
      pdj(46) =  &
          +k(253)*n(idx_PO)  &
          -k(252)*n(idx_PN)
      pdj(47) =  &
          -k(255)*n(idx_PO)  &
          -k(253)*n(idx_PO)
      pdj(53) =  &
          +k(156)*n(idx_SIOj)
      pdj(64) =  &
          -k(156)*n(idx_SIOj)
      pdj(68) =  &
          +k(234)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(11)*1d-3
      if(dnn>0.d0) then
        nn(11) = n(11) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==12) then
      pdj(6) =  &
          +k(28)*n(idx_H)  &
          -k(28)*n(idx_H)  &
          -k(93)*n(idx_H)
      pdj(8) =  &
          -k(216)*n(idx_C)
      pdj(9) =  &
          +k(60)*n(idx_N)
      pdj(10) =  &
          +k(216)*n(idx_C)
      pdj(11) =  &
          -k(60)*n(idx_N)
      pdj(12) =  &
          -k(60)*n(idx_N)  &
          -k(26)*n(idx_H2)  &
          -k(113)*n(idx_CN)  &
          -k(28)*n(idx_H)  &
          -k(93)*n(idx_H)  &
          -k(49)*n(idx_SI)  &
          -k(5)*n(idx_S)  &
          -k(254)*n(idx_P)  &
          -k(216)*n(idx_C)
      pdj(13) =  &
          -k(5)*n(idx_S)
      pdj(14) =  &
          +k(5)*n(idx_S)
      pdj(15) =  &
          +k(216)*n(idx_C)  &
          +k(60)*n(idx_N)  &
          +k(49)*n(idx_SI)  &
          +k(5)*n(idx_S)  &
          +2.d0*k(28)*n(idx_H)  &
          +k(93)*n(idx_H)  &
          +k(254)*n(idx_P)  &
          +k(113)*n(idx_CN)
      pdj(16) =  &
          -k(26)*n(idx_H2)
      pdj(17) =  &
          -k(49)*n(idx_SI)
      pdj(18) =  &
          +2.d0*k(26)*n(idx_H2)  &
          +k(93)*n(idx_H)
      pdj(24) =  &
          -k(113)*n(idx_CN)
      pdj(37) =  &
          +k(49)*n(idx_SI)
      pdj(40) =  &
          +k(113)*n(idx_CN)
      pdj(42) =  &
          -k(254)*n(idx_P)
      pdj(47) =  &
          +k(254)*n(idx_P)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(12)*1d-3
      if(dnn>0.d0) then
        nn(12) = n(12) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==13) then
      pdj(1) =  &
          -k(137)*n(idx_E)
      pdj(4) =  &
          +k(137)*n(idx_E)
      pdj(6) =  &
          +k(23)*n(idx_CH)  &
          +k(167)*n(idx_Hj)  &
          +k(91)*n(idx_H2)  &
          +k(191)*n(idx_HS)  &
          +k(27)*n(idx_OH)  &
          +k(114)*n(idx_NH)
      pdj(8) =  &
          +k(177)*n(idx_C2)  &
          +k(85)*n(idx_CN)  &
          -k(207)*n(idx_C)  &
          +k(99)*n(idx_CH)
      pdj(11) =  &
          +k(158)*n(idx_NH)
      pdj(12) =  &
          -k(5)*n(idx_O2)
      pdj(13) =  &
          -k(91)*n(idx_H2)  &
          -k(191)*n(idx_HS)  &
          -k(125)*n(idx_SO2)  &
          -k(207)*n(idx_C)  &
          -k(23)*n(idx_CH)  &
          -k(85)*n(idx_CN)  &
          -k(177)*n(idx_C2)  &
          -k(27)*n(idx_OH)  &
          -k(137)*n(idx_E)  &
          -k(114)*n(idx_NH)  &
          -k(167)*n(idx_Hj)  &
          -k(158)*n(idx_NH)  &
          -k(99)*n(idx_CH)  &
          -k(5)*n(idx_O2)
      pdj(14) =  &
          +2.d0*k(125)*n(idx_SO2)  &
          +k(27)*n(idx_OH)  &
          +k(5)*n(idx_O2)
      pdj(15) =  &
          +k(5)*n(idx_O2)
      pdj(16) =  &
          -k(91)*n(idx_H2)
      pdj(18) =  &
          -k(27)*n(idx_OH)
      pdj(19) =  &
          -k(191)*n(idx_HS)  &
          +k(158)*n(idx_NH)  &
          +k(91)*n(idx_H2)  &
          +k(99)*n(idx_CH)
      pdj(20) =  &
          +k(85)*n(idx_CN)  &
          +k(114)*n(idx_NH)
      pdj(23) =  &
          +k(23)*n(idx_CH)  &
          +k(177)*n(idx_C2)  &
          +k(207)*n(idx_C)
      pdj(24) =  &
          -k(85)*n(idx_CN)
      pdj(25) =  &
          +k(191)*n(idx_HS)
      pdj(29) =  &
          -k(99)*n(idx_CH)  &
          -k(23)*n(idx_CH)
      pdj(30) =  &
          -k(125)*n(idx_SO2)
      pdj(31) =  &
          -k(177)*n(idx_C2)
      pdj(34) =  &
          -k(114)*n(idx_NH)  &
          -k(158)*n(idx_NH)
      pdj(49) =  &
          -k(167)*n(idx_Hj)
      pdj(52) =  &
          +k(167)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(13)*1d-3
      if(dnn>0.d0) then
        nn(13) = n(13) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==14) then
      pdj(6) =  &
          +k(148)*n(idx_OH)  &
          -k(185)*n(idx_H)  &
          -k(21)*n(idx_H)
      pdj(8) =  &
          -k(80)*n(idx_C)  &
          -k(193)*n(idx_C)
      pdj(9) =  &
          +k(211)*n(idx_N)
      pdj(10) =  &
          +k(193)*n(idx_C)
      pdj(11) =  &
          -k(211)*n(idx_N)  &
          -k(94)*n(idx_N)
      pdj(12) =  &
          +k(122)*n(idx_O)
      pdj(13) =  &
          +k(193)*n(idx_C)  &
          +k(122)*n(idx_O)  &
          +k(211)*n(idx_N)  &
          +k(21)*n(idx_H)
      pdj(14) =  &
          -k(21)*n(idx_H)  &
          -k(122)*n(idx_O)  &
          -k(193)*n(idx_C)  &
          -k(211)*n(idx_N)  &
          -k(80)*n(idx_C)  &
          -k(94)*n(idx_N)  &
          -k(185)*n(idx_H)  &
          -k(148)*n(idx_OH)
      pdj(15) =  &
          +k(185)*n(idx_H)  &
          +k(94)*n(idx_N)  &
          -k(122)*n(idx_O)  &
          +k(80)*n(idx_C)
      pdj(18) =  &
          -k(148)*n(idx_OH)  &
          +k(21)*n(idx_H)
      pdj(19) =  &
          +k(185)*n(idx_H)
      pdj(20) =  &
          +k(94)*n(idx_N)
      pdj(23) =  &
          +k(80)*n(idx_C)
      pdj(30) =  &
          +k(148)*n(idx_OH)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(14)*1d-3
      if(dnn>0.d0) then
        nn(14) = n(14) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==15) then
      pdj(1) =  &
          +k(203)*n(idx_Hk)  &
          -k(13)*n(idx_E)  &
          +k(239)  &
          +k(200)*n(idx_CH)
      pdj(2) =  &
          +k(13)*n(idx_E)
      pdj(3) =  &
          -k(203)*n(idx_Hk)
      pdj(6) =  &
          +k(120)*n(idx_OH)  &
          +k(58)*n(idx_H2j)  &
          +k(45)*n(idx_NH)  &
          +k(181)*n(idx_CH)  &
          -k(83)*n(idx_H)  &
          +k(159)*n(idx_H2)  &
          +k(127)*n(idx_HCN)  &
          +k(184)*n(idx_Hj)  &
          +k(141)*n(idx_HS)
      pdj(8) =  &
          +k(168)*n(idx_C2)  &
          -k(112)*n(idx_C)  &
          +k(98)*n(idx_CH)  &
          +k(82)*n(idx_CS)  &
          +k(72)*n(idx_CN)
      pdj(9) =  &
          +k(180)*n(idx_N2)  &
          +k(212)*n(idx_NS)  &
          +k(45)*n(idx_NH)  &
          +k(72)*n(idx_CN)
      pdj(10) =  &
          +k(168)*n(idx_C2)  &
          +k(112)*n(idx_C)  &
          +k(116)*n(idx_HCN)  &
          +k(181)*n(idx_CH)  &
          +k(218)*n(idx_CN)  &
          +k(160)*n(idx_CS)
      pdj(11) =  &
          +k(180)*n(idx_N2)  &
          +k(218)*n(idx_CN)  &
          +k(41)*n(idx_NH)
      pdj(12) =  &
          +k(31)*n(idx_SO2)  &
          +k(122)*n(idx_SO)  &
          +k(120)*n(idx_OH)  &
          +2.d0*k(39)*n(idx_O)  &
          +k(97)*n(idx_SIOj)
      pdj(13) =  &
          +k(160)*n(idx_CS)  &
          +k(212)*n(idx_NS)  &
          +k(122)*n(idx_SO)  &
          +k(86)*n(idx_HS)
      pdj(14) =  &
          +k(31)*n(idx_SO2)  &
          -k(122)*n(idx_SO)  &
          +k(141)*n(idx_HS)  &
          +k(82)*n(idx_CS)
      pdj(15) =  &
          -k(82)*n(idx_CS)  &
          -k(86)*n(idx_HS)  &
          -k(31)*n(idx_SO2)  &
          -k(41)*n(idx_NH)  &
          -k(218)*n(idx_CN)  &
          -k(181)*n(idx_CH)  &
          -k(120)*n(idx_OH)  &
          -k(116)*n(idx_HCN)  &
          -k(180)*n(idx_N2)  &
          -k(127)*n(idx_HCN)  &
          -k(122)*n(idx_SO)  &
          -k(72)*n(idx_CN)  &
          -k(98)*n(idx_CH)  &
          -k(159)*n(idx_H2)  &
          -k(97)*n(idx_SIOj)  &
          -k(160)*n(idx_CS)  &
          -k(203)*n(idx_Hk)  &
          -k(186)*n(idx_SI)  &
          -k(13)*n(idx_E)  &
          -k(58)*n(idx_H2j)  &
          -k(141)*n(idx_HS)  &
          -k(200)*n(idx_CH)  &
          -k(168)*n(idx_C2)  &
          -k(212)*n(idx_NS)  &
          -k(45)*n(idx_NH)  &
          -k(184)*n(idx_Hj)  &
          -k(83)*n(idx_H)  &
          -k(112)*n(idx_C)  &
          -k(118)*n(idx_HCN)  &
          -k(239)  &
          -4.d0*k(39)*n(idx_O)  &
          -k(213)*n(idx_H2O)
      pdj(16) =  &
          -k(159)*n(idx_H2)
      pdj(17) =  &
          -k(186)*n(idx_SI)
      pdj(18) =  &
          +2.d0*k(213)*n(idx_H2O)  &
          +k(98)*n(idx_CH)  &
          +k(118)*n(idx_HCN)  &
          +k(41)*n(idx_NH)  &
          +k(83)*n(idx_H)  &
          +k(203)*n(idx_Hk)  &
          -k(120)*n(idx_OH)  &
          +k(159)*n(idx_H2)  &
          +k(86)*n(idx_HS)
      pdj(19) =  &
          -k(86)*n(idx_HS)  &
          -k(141)*n(idx_HS)
      pdj(20) =  &
          -k(212)*n(idx_NS)
      pdj(23) =  &
          -k(82)*n(idx_CS)  &
          -k(160)*n(idx_CS)
      pdj(24) =  &
          -k(218)*n(idx_CN)  &
          -k(72)*n(idx_CN)  &
          +k(118)*n(idx_HCN)
      pdj(29) =  &
          -k(200)*n(idx_CH)  &
          -k(181)*n(idx_CH)  &
          -k(98)*n(idx_CH)
      pdj(30) =  &
          -k(31)*n(idx_SO2)
      pdj(31) =  &
          -k(168)*n(idx_C2)
      pdj(32) =  &
          -k(180)*n(idx_N2)
      pdj(34) =  &
          -k(41)*n(idx_NH)  &
          +k(116)*n(idx_HCN)  &
          -k(45)*n(idx_NH)
      pdj(35) =  &
          -k(116)*n(idx_HCN)  &
          -k(127)*n(idx_HCN)  &
          -k(118)*n(idx_HCN)
      pdj(37) =  &
          +k(186)*n(idx_SI)
      pdj(40) =  &
          +k(127)*n(idx_HCN)
      pdj(44) =  &
          -k(213)*n(idx_H2O)
      pdj(49) =  &
          -k(184)*n(idx_Hj)
      pdj(53) =  &
          +k(97)*n(idx_SIOj)
      pdj(54) =  &
          +k(58)*n(idx_H2j)
      pdj(56) =  &
          -k(58)*n(idx_H2j)
      pdj(60) =  &
          +k(200)*n(idx_CH)
      pdj(62) =  &
          +k(184)*n(idx_Hj)  &
          +k(239)
      pdj(64) =  &
          -k(97)*n(idx_SIOj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(15)*1d-3
      if(dnn>0.d0) then
        nn(15) = n(15) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==16) then
      pdj(1) =  &
          +k(240)  &
          +k(236)  &
          +k(107)*n(idx_E)  &
          -k(223)*n(idx_E)  &
          -k(107)*n(idx_E)
      pdj(3) =  &
          +k(223)*n(idx_E)  &
          +k(231)
      pdj(6) =  &
          +k(117)*n(idx_N)  &
          +k(215)*n(idx_CN)  &
          +k(37)*n(idx_F)  &
          +k(149)*n(idx_OH)  &
          +k(12)*n(idx_HS)  &
          +3.d0*k(33)*n(idx_H)  &
          +k(104)*n(idx_NH)  &
          -k(33)*n(idx_H)  &
          +k(55)*n(idx_Oj)  &
          +k(159)*n(idx_O)  &
          +k(42)*n(idx_CH)  &
          +k(178)*n(idx_Sj)  &
          +k(236)  &
          +k(222)*n(idx_Hj)  &
          +2.d0*k(237)  &
          +k(54)*n(idx_C)  &
          +2.d0*k(107)*n(idx_E)  &
          +k(91)*n(idx_S)  &
          -2.d0*k(228)*n(idx_H)*n(idx_H)  &
          +k(223)*n(idx_E)
      pdj(8) =  &
          -k(75)*n(idx_C)  &
          -k(54)*n(idx_C)
      pdj(11) =  &
          -k(117)*n(idx_N)
      pdj(12) =  &
          -k(26)*n(idx_O2)
      pdj(13) =  &
          -k(91)*n(idx_S)
      pdj(15) =  &
          -k(159)*n(idx_O)
      pdj(16) =  &
          -k(222)*n(idx_Hj)  &
          -k(37)*n(idx_F)  &
          -k(33)*n(idx_H)  &
          -k(42)*n(idx_CH)  &
          -k(228)*n(idx_H)*n(idx_H)  &
          -k(107)*n(idx_E)  &
          -k(104)*n(idx_NH)  &
          -k(236)  &
          -k(159)*n(idx_O)  &
          -k(178)*n(idx_Sj)  &
          -k(55)*n(idx_Oj)  &
          -k(12)*n(idx_HS)  &
          -k(149)*n(idx_OH)  &
          -k(251)*n(idx_Fj)  &
          -k(75)*n(idx_C)  &
          -k(117)*n(idx_N)  &
          -k(26)*n(idx_O2)  &
          -k(91)*n(idx_S)  &
          +2.d0*k(228)*n(idx_H)*n(idx_H)  &
          -k(237)  &
          -k(231)  &
          -k(240)  &
          -k(54)*n(idx_C)  &
          -k(215)*n(idx_CN)  &
          -k(223)*n(idx_E)
      pdj(18) =  &
          -k(149)*n(idx_OH)  &
          +2.d0*k(26)*n(idx_O2)  &
          +k(159)*n(idx_O)
      pdj(19) =  &
          +k(91)*n(idx_S)  &
          -k(12)*n(idx_HS)
      pdj(21) =  &
          +k(12)*n(idx_HS)
      pdj(24) =  &
          -k(215)*n(idx_CN)
      pdj(27) =  &
          -k(37)*n(idx_F)  &
          +k(251)*n(idx_Fj)
      pdj(28) =  &
          +k(37)*n(idx_F)
      pdj(29) =  &
          +k(54)*n(idx_C)  &
          -k(42)*n(idx_CH)
      pdj(33) =  &
          +k(75)*n(idx_C)  &
          +k(42)*n(idx_CH)
      pdj(34) =  &
          -k(104)*n(idx_NH)  &
          +k(117)*n(idx_N)
      pdj(35) =  &
          +k(215)*n(idx_CN)
      pdj(39) =  &
          +k(104)*n(idx_NH)
      pdj(44) =  &
          +k(149)*n(idx_OH)
      pdj(49) =  &
          +k(236)  &
          -k(222)*n(idx_Hj)  &
          +k(231)
      pdj(51) =  &
          +k(178)*n(idx_Sj)
      pdj(52) =  &
          -k(178)*n(idx_Sj)
      pdj(54) =  &
          +k(55)*n(idx_Oj)
      pdj(56) =  &
          +k(240)  &
          +k(222)*n(idx_Hj)  &
          +k(251)*n(idx_Fj)
      pdj(62) =  &
          -k(55)*n(idx_Oj)
      pdj(70) =  &
          -k(251)*n(idx_Fj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(16)*1d-3
      if(dnn>0.d0) then
        nn(16) = n(16) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==17) then
      pdj(6) =  &
          +k(7)*n(idx_Hj)  &
          +k(205)*n(idx_OH)
      pdj(7) =  &
          +k(10)*n(idx_HEj)
      pdj(8) =  &
          +k(121)*n(idx_Cj)  &
          +k(89)*n(idx_CO)
      pdj(9) =  &
          -k(123)*n(idx_NO)
      pdj(10) =  &
          +k(71)*n(idx_HCOj)  &
          +k(134)*n(idx_CO2)  &
          -k(89)*n(idx_CO)
      pdj(11) =  &
          +k(123)*n(idx_NO)
      pdj(12) =  &
          -k(49)*n(idx_O2)
      pdj(13) =  &
          +k(119)*n(idx_Sj)
      pdj(15) =  &
          -k(186)*n(idx_O)  &
          +k(49)*n(idx_O2)
      pdj(17) =  &
          -k(71)*n(idx_HCOj)  &
          -k(7)*n(idx_Hj)  &
          -k(186)*n(idx_O)  &
          -k(143)*n(idx_Pj)  &
          -k(205)*n(idx_OH)  &
          -k(10)*n(idx_HEj)  &
          -k(121)*n(idx_Cj)  &
          -k(119)*n(idx_Sj)  &
          -k(49)*n(idx_O2)  &
          -k(134)*n(idx_CO2)  &
          -k(89)*n(idx_CO)  &
          -k(123)*n(idx_NO)
      pdj(18) =  &
          -k(205)*n(idx_OH)
      pdj(36) =  &
          -k(134)*n(idx_CO2)
      pdj(37) =  &
          +k(49)*n(idx_O2)  &
          +k(123)*n(idx_NO)  &
          +k(205)*n(idx_OH)  &
          +k(89)*n(idx_CO)  &
          +k(186)*n(idx_O)  &
          +k(134)*n(idx_CO2)
      pdj(42) =  &
          +k(143)*n(idx_Pj)
      pdj(48) =  &
          -k(10)*n(idx_HEj)
      pdj(49) =  &
          -k(7)*n(idx_Hj)
      pdj(52) =  &
          -k(119)*n(idx_Sj)
      pdj(53) =  &
          +k(121)*n(idx_Cj)  &
          +k(7)*n(idx_Hj)  &
          +k(143)*n(idx_Pj)  &
          +k(10)*n(idx_HEj)  &
          +k(119)*n(idx_Sj)
      pdj(58) =  &
          +k(71)*n(idx_HCOj)
      pdj(60) =  &
          -k(71)*n(idx_HCOj)
      pdj(65) =  &
          -k(143)*n(idx_Pj)
      pdj(67) =  &
          -k(121)*n(idx_Cj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(17)*1d-3
      if(dnn>0.d0) then
        nn(17) = n(17) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==18) then
      pdj(6) =  &
          +k(183)*n(idx_N)  &
          +k(120)*n(idx_O)  &
          +k(44)*n(idx_CO)  &
          +k(148)*n(idx_SO)  &
          +k(174)*n(idx_CS)  &
          +2.d0*k(67)*n(idx_H)  &
          +k(219)*n(idx_CN)  &
          -k(208)*n(idx_H)  &
          +k(81)*n(idx_SIj)  &
          +k(8)*n(idx_Hj)  &
          +k(51)*n(idx_SIO)  &
          +k(149)*n(idx_H2)  &
          +k(27)*n(idx_S)  &
          -k(67)*n(idx_H)  &
          +k(205)*n(idx_SI)  &
          +k(204)*n(idx_C)
      pdj(8) =  &
          -k(169)*n(idx_C)  &
          -k(204)*n(idx_C)
      pdj(9) =  &
          +k(183)*n(idx_N)
      pdj(10) =  &
          -k(44)*n(idx_CO)  &
          +k(77)*n(idx_CS)  &
          +k(204)*n(idx_C)
      pdj(11) =  &
          -k(155)*n(idx_N)  &
          -k(183)*n(idx_N)
      pdj(12) =  &
          +k(120)*n(idx_O)
      pdj(13) =  &
          -k(27)*n(idx_S)
      pdj(14) =  &
          -k(148)*n(idx_SO)  &
          +k(27)*n(idx_S)
      pdj(15) =  &
          -k(120)*n(idx_O)  &
          +k(155)*n(idx_N)  &
          +2.d0*k(147)*n(idx_OH)  &
          +k(20)*n(idx_F)  &
          +k(67)*n(idx_H)  &
          +k(96)*n(idx_CN)  &
          +k(208)*n(idx_H)  &
          +k(169)*n(idx_C)
      pdj(16) =  &
          +k(208)*n(idx_H)  &
          -k(149)*n(idx_H2)
      pdj(17) =  &
          -k(205)*n(idx_SI)
      pdj(18) =  &
          -k(120)*n(idx_O)  &
          -k(149)*n(idx_H2)  &
          -k(27)*n(idx_S)  &
          -k(205)*n(idx_SI)  &
          -k(44)*n(idx_CO)  &
          -k(67)*n(idx_H)  &
          -k(169)*n(idx_C)  &
          -k(183)*n(idx_N)  &
          -k(208)*n(idx_H)  &
          -k(155)*n(idx_N)  &
          -k(96)*n(idx_CN)  &
          -4.d0*k(147)*n(idx_OH)  &
          -k(174)*n(idx_CS)  &
          -k(219)*n(idx_CN)  &
          -k(77)*n(idx_CS)  &
          -k(8)*n(idx_Hj)  &
          -k(161)*n(idx_H2S)  &
          -k(81)*n(idx_SIj)  &
          -k(20)*n(idx_F)  &
          -k(51)*n(idx_SIO)  &
          -k(204)*n(idx_C)  &
          -k(148)*n(idx_SO)
      pdj(19) =  &
          +k(77)*n(idx_CS)  &
          +k(161)*n(idx_H2S)
      pdj(21) =  &
          -k(161)*n(idx_H2S)
      pdj(23) =  &
          -k(77)*n(idx_CS)  &
          -k(174)*n(idx_CS)
      pdj(24) =  &
          -k(96)*n(idx_CN)  &
          -k(219)*n(idx_CN)
      pdj(27) =  &
          -k(20)*n(idx_F)
      pdj(28) =  &
          +k(20)*n(idx_F)
      pdj(29) =  &
          +k(169)*n(idx_C)
      pdj(30) =  &
          +k(148)*n(idx_SO)
      pdj(34) =  &
          +k(155)*n(idx_N)
      pdj(35) =  &
          +k(96)*n(idx_CN)
      pdj(36) =  &
          +k(44)*n(idx_CO)
      pdj(37) =  &
          -k(51)*n(idx_SIO)  &
          +k(205)*n(idx_SI)
      pdj(38) =  &
          +k(51)*n(idx_SIO)
      pdj(40) =  &
          +k(219)*n(idx_CN)
      pdj(44) =  &
          +k(149)*n(idx_H2)  &
          +2.d0*k(147)*n(idx_OH)  &
          +k(161)*n(idx_H2S)
      pdj(45) =  &
          +k(174)*n(idx_CS)
      pdj(49) =  &
          -k(8)*n(idx_Hj)
      pdj(53) =  &
          -k(81)*n(idx_SIj)
      pdj(54) =  &
          +k(8)*n(idx_Hj)
      pdj(64) =  &
          +k(81)*n(idx_SIj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(18)*1d-3
      if(dnn>0.d0) then
        nn(18) = n(18) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==19) then
      pdj(6) =  &
          +k(191)*n(idx_S)  &
          -k(62)*n(idx_H)  &
          +k(9)*n(idx_N)  &
          +k(141)*n(idx_O)  &
          +k(12)*n(idx_H2)  &
          +k(90)*n(idx_C)
      pdj(8) =  &
          -k(90)*n(idx_C)  &
          -k(171)*n(idx_C)
      pdj(11) =  &
          -k(179)*n(idx_N)  &
          -k(9)*n(idx_N)
      pdj(13) =  &
          -k(191)*n(idx_S)  &
          +k(62)*n(idx_H)  &
          +k(171)*n(idx_C)  &
          +k(86)*n(idx_O)  &
          +2.d0*k(188)*n(idx_HS)  &
          +k(179)*n(idx_N)
      pdj(14) =  &
          +k(141)*n(idx_O)
      pdj(15) =  &
          -k(86)*n(idx_O)  &
          -k(141)*n(idx_O)
      pdj(16) =  &
          -k(12)*n(idx_H2)  &
          +k(62)*n(idx_H)
      pdj(18) =  &
          +k(86)*n(idx_O)
      pdj(19) =  &
          -k(86)*n(idx_O)  &
          -k(191)*n(idx_S)  &
          -k(90)*n(idx_C)  &
          -k(62)*n(idx_H)  &
          -k(171)*n(idx_C)  &
          -k(141)*n(idx_O)  &
          -k(179)*n(idx_N)  &
          -k(9)*n(idx_N)  &
          -k(12)*n(idx_H2)  &
          -4.d0*k(188)*n(idx_HS)
      pdj(20) =  &
          +k(9)*n(idx_N)
      pdj(21) =  &
          +2.d0*k(188)*n(idx_HS)  &
          +k(12)*n(idx_H2)
      pdj(23) =  &
          +k(90)*n(idx_C)
      pdj(25) =  &
          +k(191)*n(idx_S)
      pdj(29) =  &
          +k(171)*n(idx_C)
      pdj(34) =  &
          +k(179)*n(idx_N)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(19)*1d-3
      if(dnn>0.d0) then
        nn(19) = n(19) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==20) then
      pdj(6) =  &
          -k(217)*n(idx_H)  &
          -k(16)*n(idx_H)
      pdj(8) =  &
          -k(198)*n(idx_C)
      pdj(9) =  &
          +k(212)*n(idx_O)
      pdj(11) =  &
          +k(16)*n(idx_H)
      pdj(13) =  &
          +k(212)*n(idx_O)  &
          +k(198)*n(idx_C)  &
          +k(217)*n(idx_H)
      pdj(15) =  &
          -k(212)*n(idx_O)
      pdj(19) =  &
          +k(16)*n(idx_H)
      pdj(20) =  &
          -k(217)*n(idx_H)  &
          -k(212)*n(idx_O)  &
          -k(16)*n(idx_H)  &
          -k(198)*n(idx_C)
      pdj(24) =  &
          +k(198)*n(idx_C)
      pdj(34) =  &
          +k(217)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(20)*1d-3
      if(dnn>0.d0) then
        nn(20) = n(20) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==21) then
      pdj(6) =  &
          -k(163)*n(idx_H)
      pdj(16) =  &
          +k(163)*n(idx_H)
      pdj(18) =  &
          -k(161)*n(idx_OH)
      pdj(19) =  &
          +k(163)*n(idx_H)  &
          +k(161)*n(idx_OH)
      pdj(21) =  &
          -k(161)*n(idx_OH)  &
          -k(163)*n(idx_H)
      pdj(44) =  &
          +k(161)*n(idx_OH)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(21)*1d-3
      if(dnn>0.d0) then
        nn(21) = n(21) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==22) then
      pdj(6) =  &
          +k(172)*n(idx_Hj)
      pdj(8) =  &
          +k(124)*n(idx_Cj)
      pdj(13) =  &
          +k(14)*n(idx_Sj)
      pdj(15) =  &
          +k(68)*n(idx_Oj)
      pdj(17) =  &
          +k(136)*n(idx_SIj)
      pdj(22) =  &
          -k(79)*n(idx_SIOj)  &
          -k(124)*n(idx_Cj)  &
          -k(68)*n(idx_Oj)  &
          -k(138)*n(idx_HCOj)  &
          -k(14)*n(idx_Sj)  &
          -k(136)*n(idx_SIj)  &
          -k(172)*n(idx_Hj)
      pdj(37) =  &
          +k(79)*n(idx_SIOj)
      pdj(43) =  &
          +k(138)*n(idx_HCOj)
      pdj(49) =  &
          -k(172)*n(idx_Hj)
      pdj(52) =  &
          -k(14)*n(idx_Sj)
      pdj(53) =  &
          -k(136)*n(idx_SIj)
      pdj(57) =  &
          +k(79)*n(idx_SIOj)  &
          +k(138)*n(idx_HCOj)  &
          +k(136)*n(idx_SIj)  &
          +k(124)*n(idx_Cj)  &
          +k(14)*n(idx_Sj)  &
          +k(68)*n(idx_Oj)  &
          +k(172)*n(idx_Hj)
      pdj(60) =  &
          -k(138)*n(idx_HCOj)
      pdj(62) =  &
          -k(68)*n(idx_Oj)
      pdj(64) =  &
          -k(79)*n(idx_SIOj)
      pdj(67) =  &
          -k(124)*n(idx_Cj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(22)*1d-3
      if(dnn>0.d0) then
        nn(22) = n(22) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==23) then
      pdj(6) =  &
          +k(174)*n(idx_OH)
      pdj(8) =  &
          -k(34)*n(idx_C)  &
          +k(82)*n(idx_O)
      pdj(10) =  &
          +k(77)*n(idx_OH)  &
          +k(160)*n(idx_O)
      pdj(11) =  &
          -k(17)*n(idx_N)
      pdj(13) =  &
          +k(17)*n(idx_N)  &
          +k(34)*n(idx_C)  &
          +k(160)*n(idx_O)
      pdj(14) =  &
          +k(82)*n(idx_O)
      pdj(15) =  &
          -k(82)*n(idx_O)  &
          -k(160)*n(idx_O)
      pdj(18) =  &
          -k(77)*n(idx_OH)  &
          -k(174)*n(idx_OH)
      pdj(19) =  &
          +k(77)*n(idx_OH)
      pdj(23) =  &
          -k(82)*n(idx_O)  &
          -k(160)*n(idx_O)  &
          -k(174)*n(idx_OH)  &
          -k(34)*n(idx_C)  &
          -k(77)*n(idx_OH)  &
          -k(17)*n(idx_N)
      pdj(24) =  &
          +k(17)*n(idx_N)
      pdj(31) =  &
          +k(34)*n(idx_C)
      pdj(45) =  &
          +k(174)*n(idx_OH)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(23)*1d-3
      if(dnn>0.d0) then
        nn(23) = n(23) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==24) then
      pdj(6) =  &
          +k(219)*n(idx_OH)  &
          +k(215)*n(idx_H2)
      pdj(8) =  &
          +k(128)*n(idx_N)  &
          +k(72)*n(idx_O)  &
          -k(47)*n(idx_C)  &
          +k(85)*n(idx_S)
      pdj(9) =  &
          +k(72)*n(idx_O)
      pdj(10) =  &
          +k(218)*n(idx_O)
      pdj(11) =  &
          +k(218)*n(idx_O)  &
          +k(47)*n(idx_C)  &
          -k(128)*n(idx_N)
      pdj(12) =  &
          -k(113)*n(idx_O2)
      pdj(13) =  &
          -k(85)*n(idx_S)
      pdj(15) =  &
          -k(218)*n(idx_O)  &
          +k(96)*n(idx_OH)  &
          -k(72)*n(idx_O)  &
          +k(113)*n(idx_O2)
      pdj(16) =  &
          -k(215)*n(idx_H2)
      pdj(18) =  &
          -k(96)*n(idx_OH)  &
          -k(219)*n(idx_OH)
      pdj(20) =  &
          +k(85)*n(idx_S)
      pdj(24) =  &
          -k(113)*n(idx_O2)  &
          -k(128)*n(idx_N)  &
          -k(85)*n(idx_S)  &
          -k(47)*n(idx_C)  &
          -k(219)*n(idx_OH)  &
          -k(215)*n(idx_H2)  &
          -k(96)*n(idx_OH)  &
          -k(218)*n(idx_O)  &
          -k(72)*n(idx_O)
      pdj(31) =  &
          +k(47)*n(idx_C)
      pdj(32) =  &
          +k(128)*n(idx_N)
      pdj(35) =  &
          +k(215)*n(idx_H2)  &
          +k(96)*n(idx_OH)
      pdj(40) =  &
          +k(219)*n(idx_OH)  &
          +k(113)*n(idx_O2)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(24)*1d-3
      if(dnn>0.d0) then
        nn(24) = n(24) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==25) then
      pdj(6) =  &
          -k(18)*n(idx_H)
      pdj(13) =  &
          +k(18)*n(idx_H)
      pdj(19) =  &
          +k(18)*n(idx_H)
      pdj(25) =  &
          -k(18)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(25)*1d-3
      if(dnn>0.d0) then
        nn(25) = n(25) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==26) then
      pdj(6) =  &
          +k(230)*n(idx_Hj)
      pdj(13) =  &
          +k(19)*n(idx_Sj)
      pdj(17) =  &
          +k(132)*n(idx_SIj)
      pdj(22) =  &
          +k(78)*n(idx_FEj)
      pdj(26) =  &
          -k(132)*n(idx_SIj)  &
          -k(195)*n(idx_MGj)  &
          -k(78)*n(idx_FEj)  &
          -k(230)*n(idx_Hj)  &
          -k(19)*n(idx_Sj)
      pdj(41) =  &
          +k(195)*n(idx_MGj)
      pdj(49) =  &
          -k(230)*n(idx_Hj)
      pdj(52) =  &
          -k(19)*n(idx_Sj)
      pdj(53) =  &
          -k(132)*n(idx_SIj)
      pdj(57) =  &
          -k(78)*n(idx_FEj)
      pdj(59) =  &
          +k(230)*n(idx_Hj)  &
          +k(19)*n(idx_Sj)  &
          +k(195)*n(idx_MGj)  &
          +k(132)*n(idx_SIj)  &
          +k(78)*n(idx_FEj)
      pdj(63) =  &
          -k(195)*n(idx_MGj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(26)*1d-3
      if(dnn>0.d0) then
        nn(26) = n(26) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==27) then
      pdj(6) =  &
          +k(37)*n(idx_H2)
      pdj(15) =  &
          +k(20)*n(idx_OH)
      pdj(16) =  &
          -k(37)*n(idx_H2)
      pdj(18) =  &
          -k(20)*n(idx_OH)
      pdj(27) =  &
          -k(20)*n(idx_OH)  &
          -k(37)*n(idx_H2)
      pdj(28) =  &
          +k(37)*n(idx_H2)  &
          +k(20)*n(idx_OH)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(27)*1d-3
      if(dnn>0.d0) then
        nn(27) = n(27) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==28) then
      pdj(6) =  &
          +k(250)*n(idx_HEj)  &
          +k(100)*n(idx_SIj)
      pdj(7) =  &
          +k(250)*n(idx_HEj)
      pdj(28) =  &
          -k(250)*n(idx_HEj)  &
          -k(100)*n(idx_SIj)
      pdj(48) =  &
          -k(250)*n(idx_HEj)
      pdj(53) =  &
          -k(100)*n(idx_SIj)
      pdj(66) =  &
          +k(100)*n(idx_SIj)
      pdj(70) =  &
          +k(250)*n(idx_HEj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(28)*1d-3
      if(dnn>0.d0) then
        nn(28) = n(28) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==29) then
      pdj(1) =  &
          +k(200)*n(idx_O)
      pdj(6) =  &
          +k(152)*n(idx_C)  &
          -k(92)*n(idx_H)  &
          +k(181)*n(idx_O)  &
          +k(23)*n(idx_S)  &
          +k(165)*n(idx_N)  &
          +2.d0*k(108)*n(idx_H)  &
          +k(42)*n(idx_H2)  &
          -k(108)*n(idx_H)
      pdj(8) =  &
          +k(92)*n(idx_H)  &
          +k(131)*n(idx_N)  &
          +k(99)*n(idx_S)  &
          +k(98)*n(idx_O)  &
          -k(152)*n(idx_C)  &
          +k(108)*n(idx_H)
      pdj(10) =  &
          +k(181)*n(idx_O)
      pdj(11) =  &
          -k(131)*n(idx_N)  &
          -k(165)*n(idx_N)
      pdj(13) =  &
          -k(23)*n(idx_S)  &
          -k(99)*n(idx_S)
      pdj(15) =  &
          -k(181)*n(idx_O)  &
          -k(200)*n(idx_O)  &
          -k(98)*n(idx_O)
      pdj(16) =  &
          -k(42)*n(idx_H2)  &
          +k(92)*n(idx_H)
      pdj(18) =  &
          +k(98)*n(idx_O)
      pdj(19) =  &
          +k(99)*n(idx_S)
      pdj(23) =  &
          +k(23)*n(idx_S)
      pdj(24) =  &
          +k(165)*n(idx_N)
      pdj(29) =  &
          -k(181)*n(idx_O)  &
          -k(152)*n(idx_C)  &
          -k(42)*n(idx_H2)  &
          -k(99)*n(idx_S)  &
          -k(92)*n(idx_H)  &
          -k(131)*n(idx_N)  &
          -k(165)*n(idx_N)  &
          -k(98)*n(idx_O)  &
          -k(23)*n(idx_S)  &
          -k(200)*n(idx_O)  &
          -k(108)*n(idx_H)
      pdj(31) =  &
          +k(152)*n(idx_C)
      pdj(33) =  &
          +k(42)*n(idx_H2)
      pdj(34) =  &
          +k(131)*n(idx_N)
      pdj(60) =  &
          +k(200)*n(idx_O)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(29)*1d-3
      if(dnn>0.d0) then
        nn(29) = n(29) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==30) then
      pdj(8) =  &
          -k(189)*n(idx_C)
      pdj(10) =  &
          +k(189)*n(idx_C)
      pdj(12) =  &
          +k(31)*n(idx_O)
      pdj(13) =  &
          -k(125)*n(idx_S)
      pdj(14) =  &
          +k(189)*n(idx_C)  &
          +k(31)*n(idx_O)  &
          +2.d0*k(125)*n(idx_S)
      pdj(15) =  &
          -k(31)*n(idx_O)
      pdj(30) =  &
          -k(31)*n(idx_O)  &
          -k(125)*n(idx_S)  &
          -k(189)*n(idx_C)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(30)*1d-3
      if(dnn>0.d0) then
        nn(30) = n(30) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==31) then
      pdj(6) =  &
          -k(153)*n(idx_H)
      pdj(8) =  &
          +k(177)*n(idx_S)  &
          +k(151)*n(idx_N)  &
          +k(168)*n(idx_O)  &
          +k(153)*n(idx_H)
      pdj(10) =  &
          +k(168)*n(idx_O)
      pdj(11) =  &
          -k(151)*n(idx_N)
      pdj(13) =  &
          -k(177)*n(idx_S)
      pdj(15) =  &
          -k(168)*n(idx_O)
      pdj(23) =  &
          +k(177)*n(idx_S)
      pdj(24) =  &
          +k(151)*n(idx_N)
      pdj(29) =  &
          +k(153)*n(idx_H)
      pdj(31) =  &
          -k(177)*n(idx_S)  &
          -k(153)*n(idx_H)  &
          -k(168)*n(idx_O)  &
          -k(151)*n(idx_N)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(31)*1d-3
      if(dnn>0.d0) then
        nn(31) = n(31) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==32) then
      pdj(8) =  &
          -k(106)*n(idx_C)
      pdj(9) =  &
          +k(180)*n(idx_O)
      pdj(11) =  &
          +2.d0*k(241)  &
          +k(106)*n(idx_C)  &
          +k(180)*n(idx_O)
      pdj(15) =  &
          -k(180)*n(idx_O)
      pdj(24) =  &
          +k(106)*n(idx_C)
      pdj(32) =  &
          -k(106)*n(idx_C)  &
          -k(180)*n(idx_O)  &
          -k(241)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(32)*1d-3
      if(dnn>0.d0) then
        nn(32) = n(32) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==33) then
      pdj(6) =  &
          -k(38)*n(idx_H)
      pdj(16) =  &
          +k(38)*n(idx_H)
      pdj(29) =  &
          +k(38)*n(idx_H)
      pdj(33) =  &
          -k(38)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(33)*1d-3
      if(dnn>0.d0) then
        nn(33) = n(33) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==34) then
      pdj(6) =  &
          +k(114)*n(idx_S)  &
          -k(40)*n(idx_H)  &
          +k(104)*n(idx_H2)  &
          +k(162)*n(idx_C)  &
          +k(247)*n(idx_Hj)  &
          +k(48)*n(idx_N)  &
          +k(45)*n(idx_O)
      pdj(8) =  &
          -k(105)*n(idx_C)  &
          -k(162)*n(idx_C)
      pdj(9) =  &
          +k(45)*n(idx_O)
      pdj(11) =  &
          +k(105)*n(idx_C)  &
          +k(40)*n(idx_H)  &
          +k(41)*n(idx_O)  &
          +k(158)*n(idx_S)  &
          -k(48)*n(idx_N)
      pdj(13) =  &
          -k(114)*n(idx_S)  &
          -k(158)*n(idx_S)
      pdj(15) =  &
          -k(45)*n(idx_O)  &
          -k(41)*n(idx_O)
      pdj(16) =  &
          -k(104)*n(idx_H2)  &
          +k(40)*n(idx_H)
      pdj(18) =  &
          +k(41)*n(idx_O)
      pdj(19) =  &
          +k(158)*n(idx_S)
      pdj(20) =  &
          +k(114)*n(idx_S)
      pdj(24) =  &
          +k(162)*n(idx_C)
      pdj(29) =  &
          +k(105)*n(idx_C)
      pdj(32) =  &
          +k(48)*n(idx_N)
      pdj(34) =  &
          -k(114)*n(idx_S)  &
          -k(40)*n(idx_H)  &
          -k(105)*n(idx_C)  &
          -k(158)*n(idx_S)  &
          -k(41)*n(idx_O)  &
          -k(162)*n(idx_C)  &
          -k(48)*n(idx_N)  &
          -k(104)*n(idx_H2)  &
          -k(45)*n(idx_O)  &
          -k(247)*n(idx_Hj)
      pdj(39) =  &
          +k(104)*n(idx_H2)
      pdj(49) =  &
          -k(247)*n(idx_Hj)
      pdj(50) =  &
          +k(247)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(34)*1d-3
      if(dnn>0.d0) then
        nn(34) = n(34) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==35) then
      pdj(6) =  &
          +k(127)*n(idx_O)  &
          -k(43)*n(idx_H)
      pdj(10) =  &
          +k(116)*n(idx_O)
      pdj(15) =  &
          -k(118)*n(idx_O)  &
          -k(116)*n(idx_O)  &
          -k(127)*n(idx_O)
      pdj(16) =  &
          +k(43)*n(idx_H)
      pdj(18) =  &
          +k(118)*n(idx_O)
      pdj(24) =  &
          +k(43)*n(idx_H)  &
          +k(118)*n(idx_O)
      pdj(34) =  &
          +k(116)*n(idx_O)
      pdj(35) =  &
          -k(118)*n(idx_O)  &
          -k(116)*n(idx_O)  &
          -k(43)*n(idx_H)  &
          -k(127)*n(idx_O)
      pdj(40) =  &
          +k(127)*n(idx_O)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(35)*1d-3
      if(dnn>0.d0) then
        nn(35) = n(35) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==36) then
      pdj(6) =  &
          -k(176)*n(idx_H)
      pdj(9) =  &
          +k(175)*n(idx_N)
      pdj(10) =  &
          +k(176)*n(idx_H)  &
          +k(175)*n(idx_N)  &
          +k(134)*n(idx_SI)
      pdj(11) =  &
          -k(175)*n(idx_N)
      pdj(17) =  &
          -k(134)*n(idx_SI)
      pdj(18) =  &
          +k(176)*n(idx_H)
      pdj(36) =  &
          -k(176)*n(idx_H)  &
          -k(134)*n(idx_SI)  &
          -k(175)*n(idx_N)
      pdj(37) =  &
          +k(134)*n(idx_SI)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(36)*1d-3
      if(dnn>0.d0) then
        nn(36) = n(36) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==37) then
      pdj(6) =  &
          +k(140)*n(idx_Hj)  &
          +k(51)*n(idx_OH)
      pdj(18) =  &
          -k(51)*n(idx_OH)
      pdj(37) =  &
          -k(140)*n(idx_Hj)  &
          -k(51)*n(idx_OH)
      pdj(38) =  &
          +k(51)*n(idx_OH)
      pdj(49) =  &
          -k(140)*n(idx_Hj)
      pdj(64) =  &
          +k(140)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(37)*1d-3
      if(dnn>0.d0) then
        nn(37) = n(37) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==38) then
      pdj(7) =  &
          +k(246)*n(idx_HEj)
      pdj(12) =  &
          +k(246)*n(idx_HEj)
      pdj(38) =  &
          -k(246)*n(idx_HEj)
      pdj(48) =  &
          -k(246)*n(idx_HEj)
      pdj(53) =  &
          +k(246)*n(idx_HEj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(38)*1d-3
      if(dnn>0.d0) then
        nn(38) = n(38) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==39) then
      pdj(6) =  &
          -k(53)*n(idx_H)
      pdj(16) =  &
          +k(53)*n(idx_H)
      pdj(34) =  &
          +k(53)*n(idx_H)
      pdj(39) =  &
          -k(53)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(39)*1d-3
      if(dnn>0.d0) then
        nn(39) = n(39) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==40) then
      pdj(6) =  &
          -k(56)*n(idx_H)  &
          -k(150)*n(idx_H)  &
          -k(154)*n(idx_H)
      pdj(10) =  &
          +k(154)*n(idx_H)
      pdj(15) =  &
          +k(150)*n(idx_H)
      pdj(18) =  &
          +k(56)*n(idx_H)
      pdj(24) =  &
          +k(56)*n(idx_H)
      pdj(34) =  &
          +k(154)*n(idx_H)
      pdj(35) =  &
          +k(150)*n(idx_H)
      pdj(40) =  &
          -k(56)*n(idx_H)  &
          -k(150)*n(idx_H)  &
          -k(154)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(40)*1d-3
      if(dnn>0.d0) then
        nn(40) = n(40) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==41) then
      pdj(6) =  &
          +k(144)*n(idx_Hj)
      pdj(8) =  &
          +k(111)*n(idx_Cj)
      pdj(13) =  &
          +k(59)*n(idx_Sj)
      pdj(17) =  &
          +k(130)*n(idx_SIj)
      pdj(37) =  &
          +k(126)*n(idx_SIOj)
      pdj(41) =  &
          -k(164)*n(idx_HCOj)  &
          -k(59)*n(idx_Sj)  &
          -k(111)*n(idx_Cj)  &
          -k(144)*n(idx_Hj)  &
          -k(126)*n(idx_SIOj)  &
          -k(130)*n(idx_SIj)
      pdj(43) =  &
          +k(164)*n(idx_HCOj)
      pdj(49) =  &
          -k(144)*n(idx_Hj)
      pdj(52) =  &
          -k(59)*n(idx_Sj)
      pdj(53) =  &
          -k(130)*n(idx_SIj)
      pdj(60) =  &
          -k(164)*n(idx_HCOj)
      pdj(63) =  &
          +k(130)*n(idx_SIj)  &
          +k(164)*n(idx_HCOj)  &
          +k(111)*n(idx_Cj)  &
          +k(126)*n(idx_SIOj)  &
          +k(59)*n(idx_Sj)  &
          +k(144)*n(idx_Hj)
      pdj(64) =  &
          -k(126)*n(idx_SIOj)
      pdj(67) =  &
          -k(111)*n(idx_Cj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(41)*1d-3
      if(dnn>0.d0) then
        nn(41) = n(41) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==42) then
      pdj(6) =  &
          +k(84)*n(idx_Hj)
      pdj(12) =  &
          -k(254)*n(idx_O2)
      pdj(15) =  &
          +k(254)*n(idx_O2)
      pdj(42) =  &
          -k(254)*n(idx_O2)  &
          -k(84)*n(idx_Hj)
      pdj(47) =  &
          +k(254)*n(idx_O2)
      pdj(49) =  &
          -k(84)*n(idx_Hj)
      pdj(65) =  &
          +k(84)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(42)*1d-3
      if(dnn>0.d0) then
        nn(42) = n(42) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==43) then
      pdj(6) =  &
          -k(95)*n(idx_H)
      pdj(10) =  &
          +k(95)*n(idx_H)
      pdj(16) =  &
          +k(95)*n(idx_H)
      pdj(43) =  &
          -k(95)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(43)*1d-3
      if(dnn>0.d0) then
        nn(43) = n(43) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==44) then
      pdj(6) =  &
          -k(103)*n(idx_H)  &
          -k(146)*n(idx_H)  &
          +2.d0*k(103)*n(idx_H)
      pdj(15) =  &
          -k(213)*n(idx_O)
      pdj(16) =  &
          +k(146)*n(idx_H)
      pdj(18) =  &
          +k(146)*n(idx_H)  &
          +k(103)*n(idx_H)  &
          +2.d0*k(213)*n(idx_O)
      pdj(44) =  &
          -k(103)*n(idx_H)  &
          -k(213)*n(idx_O)  &
          -k(146)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(44)*1d-3
      if(dnn>0.d0) then
        nn(44) = n(44) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==45) then
      pdj(6) =  &
          -k(135)*n(idx_H)
      pdj(10) =  &
          +k(135)*n(idx_H)
      pdj(19) =  &
          +k(135)*n(idx_H)
      pdj(45) =  &
          -k(135)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(45)*1d-3
      if(dnn>0.d0) then
        nn(45) = n(45) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==46) then
      pdj(11) =  &
          -k(252)*n(idx_N)
      pdj(32) =  &
          +k(252)*n(idx_N)
      pdj(42) =  &
          +k(252)*n(idx_N)
      pdj(46) =  &
          -k(252)*n(idx_N)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(46)*1d-3
      if(dnn>0.d0) then
        nn(46) = n(46) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==47) then
      pdj(9) =  &
          +k(255)*n(idx_N)
      pdj(11) =  &
          -k(255)*n(idx_N)  &
          -k(253)*n(idx_N)
      pdj(15) =  &
          +k(253)*n(idx_N)
      pdj(42) =  &
          +k(255)*n(idx_N)
      pdj(46) =  &
          +k(253)*n(idx_N)
      pdj(47) =  &
          -k(255)*n(idx_N)  &
          -k(253)*n(idx_N)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(47)*1d-3
      if(dnn>0.d0) then
        nn(47) = n(47) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==48) then
      pdj(1) =  &
          -k(4)*n(idx_E)
      pdj(6) =  &
          -k(1)*n(idx_H)  &
          +k(250)*n(idx_HF)
      pdj(7) =  &
          +k(10)*n(idx_SI)  &
          +k(246)*n(idx_SIO2)  &
          +k(4)*n(idx_E)  &
          +k(250)*n(idx_HF)  &
          +k(1)*n(idx_H)
      pdj(12) =  &
          +k(246)*n(idx_SIO2)
      pdj(17) =  &
          -k(10)*n(idx_SI)
      pdj(28) =  &
          -k(250)*n(idx_HF)
      pdj(38) =  &
          -k(246)*n(idx_SIO2)
      pdj(48) =  &
          -k(250)*n(idx_HF)  &
          -k(1)*n(idx_H)  &
          -k(10)*n(idx_SI)  &
          -k(246)*n(idx_SIO2)  &
          -k(4)*n(idx_E)
      pdj(49) =  &
          +k(1)*n(idx_H)
      pdj(53) =  &
          +k(10)*n(idx_SI)  &
          +k(246)*n(idx_SIO2)
      pdj(70) =  &
          +k(250)*n(idx_HF)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(48)*1d-3
      if(dnn>0.d0) then
        nn(48) = n(48) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==49) then
      pdj(1) =  &
          -k(29)*n(idx_E)  &
          +k(226)*n(idx_Hk)
      pdj(2) =  &
          -k(199)*n(idx_Ok)
      pdj(3) =  &
          -k(65)*n(idx_Hk)  &
          -k(226)*n(idx_Hk)
      pdj(5) =  &
          -k(133)*n(idx_Ck)
      pdj(6) =  &
          +k(7)*n(idx_SI)  &
          +k(133)*n(idx_Ck)  &
          +2.d0*k(65)*n(idx_Hk)  &
          +k(84)*n(idx_P)  &
          +k(167)*n(idx_S)  &
          -k(61)*n(idx_H)  &
          +k(184)*n(idx_O)  &
          +k(144)*n(idx_MG)  &
          +k(140)*n(idx_SIO)  &
          +k(222)*n(idx_H2)  &
          +k(8)*n(idx_OH)  &
          +k(29)*n(idx_E)  &
          +k(230)*n(idx_NA)  &
          +k(172)*n(idx_FE)  &
          +k(199)*n(idx_Ok)  &
          +k(247)*n(idx_NH)
      pdj(8) =  &
          +k(133)*n(idx_Ck)
      pdj(13) =  &
          -k(167)*n(idx_S)
      pdj(15) =  &
          +k(199)*n(idx_Ok)  &
          -k(184)*n(idx_O)
      pdj(16) =  &
          -k(222)*n(idx_H2)
      pdj(17) =  &
          -k(7)*n(idx_SI)
      pdj(18) =  &
          -k(8)*n(idx_OH)
      pdj(22) =  &
          -k(172)*n(idx_FE)
      pdj(26) =  &
          -k(230)*n(idx_NA)
      pdj(34) =  &
          -k(247)*n(idx_NH)
      pdj(37) =  &
          -k(140)*n(idx_SIO)
      pdj(41) =  &
          -k(144)*n(idx_MG)
      pdj(42) =  &
          -k(84)*n(idx_P)
      pdj(49) =  &
          -k(230)*n(idx_NA)  &
          -k(226)*n(idx_Hk)  &
          -k(199)*n(idx_Ok)  &
          -k(167)*n(idx_S)  &
          -k(8)*n(idx_OH)  &
          -k(29)*n(idx_E)  &
          -k(61)*n(idx_H)  &
          -k(222)*n(idx_H2)  &
          -k(144)*n(idx_MG)  &
          -k(84)*n(idx_P)  &
          -k(184)*n(idx_O)  &
          -k(247)*n(idx_NH)  &
          -k(65)*n(idx_Hk)  &
          -k(133)*n(idx_Ck)  &
          -k(140)*n(idx_SIO)  &
          -k(172)*n(idx_FE)  &
          -k(7)*n(idx_SI)
      pdj(50) =  &
          +k(247)*n(idx_NH)
      pdj(52) =  &
          +k(167)*n(idx_S)
      pdj(53) =  &
          +k(7)*n(idx_SI)
      pdj(54) =  &
          +k(8)*n(idx_OH)
      pdj(56) =  &
          +k(61)*n(idx_H)  &
          +k(222)*n(idx_H2)  &
          +k(226)*n(idx_Hk)
      pdj(57) =  &
          +k(172)*n(idx_FE)
      pdj(59) =  &
          +k(230)*n(idx_NA)
      pdj(62) =  &
          +k(184)*n(idx_O)
      pdj(63) =  &
          +k(144)*n(idx_MG)
      pdj(64) =  &
          +k(140)*n(idx_SIO)
      pdj(65) =  &
          +k(84)*n(idx_P)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(49)*1d-3
      if(dnn>0.d0) then
        nn(49) = n(49) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==50) then
      pdj(1) =  &
          -k(3)*n(idx_E)
      pdj(6) =  &
          +k(3)*n(idx_E)
      pdj(11) =  &
          +k(3)*n(idx_E)
      pdj(50) =  &
          -k(3)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(50)*1d-3
      if(dnn>0.d0) then
        nn(50) = n(50) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==51) then
      pdj(1) =  &
          -k(248)*n(idx_E)
      pdj(6) =  &
          -k(6)*n(idx_H)  &
          +k(248)*n(idx_E)
      pdj(13) =  &
          +k(248)*n(idx_E)
      pdj(16) =  &
          +k(6)*n(idx_H)
      pdj(51) =  &
          -k(6)*n(idx_H)  &
          -k(248)*n(idx_E)
      pdj(52) =  &
          +k(6)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(51)*1d-3
      if(dnn>0.d0) then
        nn(51) = n(51) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==52) then
      pdj(1) =  &
          -k(187)*n(idx_E)
      pdj(3) =  &
          -k(66)*n(idx_Hk)
      pdj(6) =  &
          +k(66)*n(idx_Hk)  &
          +k(178)*n(idx_H2)
      pdj(13) =  &
          +k(19)*n(idx_NA)  &
          +k(14)*n(idx_FE)  &
          +k(119)*n(idx_SI)  &
          +k(187)*n(idx_E)  &
          +k(59)*n(idx_MG)  &
          +k(66)*n(idx_Hk)
      pdj(16) =  &
          -k(178)*n(idx_H2)
      pdj(17) =  &
          -k(119)*n(idx_SI)
      pdj(22) =  &
          -k(14)*n(idx_FE)
      pdj(26) =  &
          -k(19)*n(idx_NA)
      pdj(41) =  &
          -k(59)*n(idx_MG)
      pdj(51) =  &
          +k(178)*n(idx_H2)
      pdj(52) =  &
          -k(119)*n(idx_SI)  &
          -k(59)*n(idx_MG)  &
          -k(187)*n(idx_E)  &
          -k(19)*n(idx_NA)  &
          -k(14)*n(idx_FE)  &
          -k(178)*n(idx_H2)  &
          -k(66)*n(idx_Hk)
      pdj(53) =  &
          +k(119)*n(idx_SI)
      pdj(57) =  &
          +k(14)*n(idx_FE)
      pdj(59) =  &
          +k(19)*n(idx_NA)
      pdj(63) =  &
          +k(59)*n(idx_MG)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(52)*1d-3
      if(dnn>0.d0) then
        nn(52) = n(52) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==53) then
      pdj(1) =  &
          -k(197)*n(idx_E)
      pdj(3) =  &
          -k(145)*n(idx_Hk)
      pdj(6) =  &
          +k(81)*n(idx_OH)  &
          -k(15)*n(idx_H)  &
          +k(100)*n(idx_HF)  &
          +k(145)*n(idx_Hk)
      pdj(17) =  &
          +k(136)*n(idx_FE)  &
          +k(145)*n(idx_Hk)  &
          +k(132)*n(idx_NA)  &
          +k(197)*n(idx_E)  &
          +k(130)*n(idx_MG)
      pdj(18) =  &
          -k(81)*n(idx_OH)
      pdj(22) =  &
          -k(136)*n(idx_FE)
      pdj(26) =  &
          -k(132)*n(idx_NA)
      pdj(28) =  &
          -k(100)*n(idx_HF)
      pdj(41) =  &
          -k(130)*n(idx_MG)
      pdj(53) =  &
          -k(100)*n(idx_HF)  &
          -k(197)*n(idx_E)  &
          -k(145)*n(idx_Hk)  &
          -k(130)*n(idx_MG)  &
          -k(136)*n(idx_FE)  &
          -k(132)*n(idx_NA)  &
          -k(15)*n(idx_H)  &
          -k(81)*n(idx_OH)
      pdj(57) =  &
          +k(136)*n(idx_FE)
      pdj(58) =  &
          +k(15)*n(idx_H)
      pdj(59) =  &
          +k(132)*n(idx_NA)
      pdj(63) =  &
          +k(130)*n(idx_MG)
      pdj(64) =  &
          +k(81)*n(idx_OH)
      pdj(66) =  &
          +k(100)*n(idx_HF)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(53)*1d-3
      if(dnn>0.d0) then
        nn(53) = n(53) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==54) then
      pdj(1) =  &
          -k(87)*n(idx_E)
      pdj(6) =  &
          +k(87)*n(idx_E)
      pdj(15) =  &
          +k(87)*n(idx_E)
      pdj(54) =  &
          -k(87)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(54)*1d-3
      if(dnn>0.d0) then
        nn(54) = n(54) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==55) then
      pdj(1) =  &
          -k(25)*n(idx_E)
      pdj(6) =  &
          +k(25)*n(idx_E)  &
          -k(11)*n(idx_H)
      pdj(7) =  &
          +k(11)*n(idx_H)  &
          +k(25)*n(idx_E)
      pdj(55) =  &
          -k(25)*n(idx_E)  &
          -k(11)*n(idx_H)
      pdj(56) =  &
          +k(11)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(55)*1d-3
      if(dnn>0.d0) then
        nn(55) = n(55) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==56) then
      pdj(1) =  &
          -k(102)*n(idx_E)
      pdj(6) =  &
          +2.d0*k(102)*n(idx_E)  &
          +k(58)*n(idx_O)  &
          -k(196)*n(idx_H)  &
          +k(32)*n(idx_HE)  &
          +k(206)*n(idx_C)
      pdj(7) =  &
          -k(32)*n(idx_HE)
      pdj(8) =  &
          -k(206)*n(idx_C)
      pdj(15) =  &
          -k(58)*n(idx_O)
      pdj(16) =  &
          +k(196)*n(idx_H)
      pdj(49) =  &
          +k(196)*n(idx_H)
      pdj(54) =  &
          +k(58)*n(idx_O)
      pdj(55) =  &
          +k(32)*n(idx_HE)
      pdj(56) =  &
          -k(206)*n(idx_C)  &
          -k(32)*n(idx_HE)  &
          -k(58)*n(idx_O)  &
          -k(102)*n(idx_E)  &
          -k(196)*n(idx_H)
      pdj(61) =  &
          +k(206)*n(idx_C)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(56)*1d-3
      if(dnn>0.d0) then
        nn(56) = n(56) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==57) then
      pdj(1) =  &
          -k(142)*n(idx_E)
      pdj(2) =  &
          -k(210)*n(idx_Ok)
      pdj(3) =  &
          -k(24)*n(idx_Hk)
      pdj(6) =  &
          +k(24)*n(idx_Hk)
      pdj(15) =  &
          +k(210)*n(idx_Ok)
      pdj(22) =  &
          +k(142)*n(idx_E)  &
          +k(24)*n(idx_Hk)  &
          +k(210)*n(idx_Ok)  &
          +k(78)*n(idx_NA)
      pdj(26) =  &
          -k(78)*n(idx_NA)
      pdj(57) =  &
          -k(210)*n(idx_Ok)  &
          -k(78)*n(idx_NA)  &
          -k(142)*n(idx_E)  &
          -k(24)*n(idx_Hk)
      pdj(59) =  &
          +k(78)*n(idx_NA)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(57)*1d-3
      if(dnn>0.d0) then
        nn(57) = n(57) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==58) then
      pdj(1) =  &
          -k(201)*n(idx_E)
      pdj(6) =  &
          +k(201)*n(idx_E)  &
          -k(110)*n(idx_H)
      pdj(16) =  &
          +k(110)*n(idx_H)
      pdj(17) =  &
          +k(201)*n(idx_E)
      pdj(53) =  &
          +k(110)*n(idx_H)
      pdj(58) =  &
          -k(201)*n(idx_E)  &
          -k(110)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(58)*1d-3
      if(dnn>0.d0) then
        nn(58) = n(58) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==59) then
      pdj(1) =  &
          -k(190)*n(idx_E)
      pdj(3) =  &
          -k(22)*n(idx_Hk)
      pdj(6) =  &
          +k(22)*n(idx_Hk)
      pdj(26) =  &
          +k(190)*n(idx_E)  &
          +k(22)*n(idx_Hk)
      pdj(59) =  &
          -k(22)*n(idx_Hk)  &
          -k(190)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(59)*1d-3
      if(dnn>0.d0) then
        nn(59) = n(59) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==60) then
      pdj(1) =  &
          -k(249)*n(idx_E)
      pdj(6) =  &
          +k(249)*n(idx_E)
      pdj(8) =  &
          -k(52)*n(idx_C)
      pdj(10) =  &
          +k(71)*n(idx_SI)  &
          +k(52)*n(idx_C)  &
          +k(249)*n(idx_E)
      pdj(17) =  &
          -k(71)*n(idx_SI)
      pdj(22) =  &
          -k(138)*n(idx_FE)
      pdj(41) =  &
          -k(164)*n(idx_MG)
      pdj(43) =  &
          +k(164)*n(idx_MG)  &
          +k(138)*n(idx_FE)
      pdj(57) =  &
          +k(138)*n(idx_FE)
      pdj(58) =  &
          +k(71)*n(idx_SI)
      pdj(60) =  &
          -k(71)*n(idx_SI)  &
          -k(52)*n(idx_C)  &
          -k(164)*n(idx_MG)  &
          -k(138)*n(idx_FE)  &
          -k(249)*n(idx_E)
      pdj(61) =  &
          +k(52)*n(idx_C)
      pdj(63) =  &
          +k(164)*n(idx_MG)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(60)*1d-3
      if(dnn>0.d0) then
        nn(60) = n(60) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==61) then
      pdj(1) =  &
          -k(209)*n(idx_E)
      pdj(6) =  &
          -k(194)*n(idx_H)  &
          +k(209)*n(idx_E)
      pdj(8) =  &
          +k(209)*n(idx_E)
      pdj(16) =  &
          +k(194)*n(idx_H)
      pdj(61) =  &
          -k(194)*n(idx_H)  &
          -k(209)*n(idx_E)
      pdj(67) =  &
          +k(194)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(61)*1d-3
      if(dnn>0.d0) then
        nn(61) = n(61) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==62) then
      pdj(1) =  &
          -k(74)*n(idx_E)
      pdj(3) =  &
          -k(73)*n(idx_Hk)
      pdj(6) =  &
          +k(73)*n(idx_Hk)  &
          -k(182)*n(idx_H)  &
          +k(55)*n(idx_H2)
      pdj(15) =  &
          +k(68)*n(idx_FE)  &
          +k(74)*n(idx_E)  &
          +k(182)*n(idx_H)  &
          +k(73)*n(idx_Hk)
      pdj(16) =  &
          -k(55)*n(idx_H2)
      pdj(22) =  &
          -k(68)*n(idx_FE)
      pdj(49) =  &
          +k(182)*n(idx_H)
      pdj(54) =  &
          +k(55)*n(idx_H2)
      pdj(57) =  &
          +k(68)*n(idx_FE)
      pdj(62) =  &
          -k(55)*n(idx_H2)  &
          -k(74)*n(idx_E)  &
          -k(73)*n(idx_Hk)  &
          -k(182)*n(idx_H)  &
          -k(68)*n(idx_FE)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(62)*1d-3
      if(dnn>0.d0) then
        nn(62) = n(62) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==63) then
      pdj(1) =  &
          -k(64)*n(idx_E)
      pdj(2) =  &
          -k(63)*n(idx_Ok)
      pdj(3) =  &
          -k(214)*n(idx_Hk)
      pdj(6) =  &
          +k(214)*n(idx_Hk)
      pdj(15) =  &
          +k(63)*n(idx_Ok)
      pdj(26) =  &
          -k(195)*n(idx_NA)
      pdj(41) =  &
          +k(195)*n(idx_NA)  &
          +k(214)*n(idx_Hk)  &
          +k(64)*n(idx_E)  &
          +k(63)*n(idx_Ok)
      pdj(59) =  &
          +k(195)*n(idx_NA)
      pdj(63) =  &
          -k(64)*n(idx_E)  &
          -k(214)*n(idx_Hk)  &
          -k(63)*n(idx_Ok)  &
          -k(195)*n(idx_NA)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(63)*1d-3
      if(dnn>0.d0) then
        nn(63) = n(63) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==64) then
      pdj(1) =  &
          -k(115)*n(idx_E)
      pdj(8) =  &
          -k(202)*n(idx_C)
      pdj(9) =  &
          +k(156)*n(idx_N)
      pdj(10) =  &
          +k(202)*n(idx_C)
      pdj(11) =  &
          -k(156)*n(idx_N)
      pdj(12) =  &
          +k(97)*n(idx_O)
      pdj(15) =  &
          -k(97)*n(idx_O)  &
          +k(115)*n(idx_E)
      pdj(17) =  &
          +k(115)*n(idx_E)
      pdj(22) =  &
          -k(79)*n(idx_FE)
      pdj(37) =  &
          +k(126)*n(idx_MG)  &
          +k(79)*n(idx_FE)
      pdj(41) =  &
          -k(126)*n(idx_MG)
      pdj(53) =  &
          +k(97)*n(idx_O)  &
          +k(156)*n(idx_N)  &
          +k(202)*n(idx_C)
      pdj(57) =  &
          +k(79)*n(idx_FE)
      pdj(63) =  &
          +k(126)*n(idx_MG)
      pdj(64) =  &
          -k(97)*n(idx_O)  &
          -k(156)*n(idx_N)  &
          -k(126)*n(idx_MG)  &
          -k(202)*n(idx_C)  &
          -k(79)*n(idx_FE)  &
          -k(115)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(64)*1d-3
      if(dnn>0.d0) then
        nn(64) = n(64) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==65) then
      pdj(1) =  &
          -k(139)*n(idx_E)
      pdj(17) =  &
          -k(143)*n(idx_SI)
      pdj(42) =  &
          +k(139)*n(idx_E)  &
          +k(143)*n(idx_SI)
      pdj(53) =  &
          +k(143)*n(idx_SI)
      pdj(65) =  &
          -k(139)*n(idx_E)  &
          -k(143)*n(idx_SI)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(65)*1d-3
      if(dnn>0.d0) then
        nn(65) = n(65) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==66) then
      pdj(1) =  &
          -k(109)*n(idx_E)
      pdj(17) =  &
          +k(109)*n(idx_E)
      pdj(27) =  &
          +k(109)*n(idx_E)
      pdj(66) =  &
          -k(109)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(66)*1d-3
      if(dnn>0.d0) then
        nn(66) = n(66) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==67) then
      pdj(1) =  &
          -k(243)*n(idx_E)
      pdj(8) =  &
          +k(121)*n(idx_SI)  &
          +k(111)*n(idx_MG)  &
          +k(124)*n(idx_FE)  &
          +k(243)*n(idx_E)
      pdj(17) =  &
          -k(121)*n(idx_SI)
      pdj(22) =  &
          -k(124)*n(idx_FE)
      pdj(41) =  &
          -k(111)*n(idx_MG)
      pdj(53) =  &
          +k(121)*n(idx_SI)
      pdj(57) =  &
          +k(124)*n(idx_FE)
      pdj(63) =  &
          +k(111)*n(idx_MG)
      pdj(67) =  &
          -k(111)*n(idx_MG)  &
          -k(243)*n(idx_E)  &
          -k(121)*n(idx_SI)  &
          -k(124)*n(idx_FE)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(67)*1d-3
      if(dnn>0.d0) then
        nn(67) = n(67) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==68) then
      pdj(1) =  &
          -k(244)*n(idx_E)
      pdj(11) =  &
          +k(244)*n(idx_E)
      pdj(68) =  &
          -k(244)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(68)*1d-3
      if(dnn>0.d0) then
        nn(68) = n(68) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==69) then
      pdj(1) =  &
          -k(245)*n(idx_E)
      pdj(8) =  &
          +k(245)*n(idx_E)
      pdj(15) =  &
          +k(245)*n(idx_E)
      pdj(69) =  &
          -k(245)*n(idx_E)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(69)*1d-3
      if(dnn>0.d0) then
        nn(69) = n(69) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==70) then
      pdj(16) =  &
          -k(251)*n(idx_H2)
      pdj(27) =  &
          +k(251)*n(idx_H2)
      pdj(56) =  &
          +k(251)*n(idx_H2)
      pdj(70) =  &
          -k(251)*n(idx_H2)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(70)*1d-3
      if(dnn>0.d0) then
        nn(70) = n(70) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==71) then
      pdj(73) = 0.d0
    elseif(j==72) then
      pdj(73) = 0.d0
    elseif(j==73) then
      !use fex to compute temperature-dependent Jacobian
      dnn = n(idx_Tgas)*1d-3
      nn(:) = n(:)
      nn(idx_Tgas) = n(idx_Tgas) + dnn
      call fex(neq,tt,nn(:),dn(:))
      do i=1,neq-1
        pdj(i) = dn(i) / dnn
      end do
    elseif(j==74) then
      pdj(73) = 0.d0
    end if

    return
  end subroutine jes

  !*************************
  subroutine jex(neq,t,n,ml,mu,pd,npd)
    use krome_commons
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_subs
    use krome_gadiab
    implicit none
    real*8::n(neq),pd(neq,neq),t,k(nrea),dn0,dn1,dnn,Tgas
    real*8::krome_gamma,nn(neq),nH2dust
    integer::neq,ml,mu,npd

    Tgas = n(idx_Tgas)
    npd = neq
    k(:) = coe_tab(n(:))
    pd(:,:) = 0d0
    krome_gamma = gamma_index(n(:))

    !d[E_dot]/d[E]
    pd(1,1) =  &
        -k(187)*n(idx_Sj)  &
        +2.d0*k(220)*n(idx_H)  &
        -k(157)*n(idx_C)  &
        -k(244)*n(idx_Nj)  &
        -k(25)*n(idx_HEHj)  &
        -k(4)*n(idx_HEj)  &
        -k(220)*n(idx_H)  &
        -k(201)*n(idx_SIHj)  &
        -k(245)*n(idx_COj)  &
        -k(139)*n(idx_Pj)  &
        -k(221)*n(idx_HE)  &
        -k(248)*n(idx_HSj)  &
        -k(29)*n(idx_Hj)  &
        +2.d0*k(221)*n(idx_HE)  &
        -k(115)*n(idx_SIOj)  &
        -k(142)*n(idx_FEj)  &
        -k(190)*n(idx_NAj)  &
        -k(102)*n(idx_H2j)  &
        -k(13)*n(idx_O)  &
        -k(107)*n(idx_H2)  &
        -k(3)*n(idx_NHj)  &
        -k(249)*n(idx_HCOj)  &
        -k(223)*n(idx_H2)  &
        -k(197)*n(idx_SIj)  &
        -k(224)*n(idx_Hk)  &
        +2.d0*k(224)*n(idx_Hk)  &
        -k(74)*n(idx_Oj)  &
        -k(64)*n(idx_MGj)  &
        -k(109)*n(idx_SIFj)  &
        -k(57)*n(idx_H)  &
        -k(87)*n(idx_OHj)  &
        -k(243)*n(idx_Cj)  &
        -k(137)*n(idx_S)  &
        +k(107)*n(idx_H2)  &
        -k(209)*n(idx_CHj)

    !d[O-_dot]/d[E]
    pd(2,1) =  &
        +k(13)*n(idx_O)

    !d[H-_dot]/d[E]
    pd(3,1) =  &
        -k(224)*n(idx_Hk)  &
        +k(57)*n(idx_H)  &
        +k(223)*n(idx_H2)

    !d[S-_dot]/d[E]
    pd(4,1) =  &
        +k(137)*n(idx_S)

    !d[C-_dot]/d[E]
    pd(5,1) =  &
        +k(157)*n(idx_C)

    !d[H_dot]/d[E]
    pd(6,1) =  &
        +k(29)*n(idx_Hj)  &
        +k(223)*n(idx_H2)  &
        +k(201)*n(idx_SIHj)  &
        -k(57)*n(idx_H)  &
        +k(209)*n(idx_CHj)  &
        +k(248)*n(idx_HSj)  &
        +k(25)*n(idx_HEHj)  &
        +k(224)*n(idx_Hk)  &
        +2.d0*k(102)*n(idx_H2j)  &
        -k(220)*n(idx_H)  &
        +k(3)*n(idx_NHj)  &
        +k(87)*n(idx_OHj)  &
        +2.d0*k(107)*n(idx_H2)  &
        +k(249)*n(idx_HCOj)

    !d[HE_dot]/d[E]
    pd(7,1) =  &
        +k(25)*n(idx_HEHj)  &
        +k(4)*n(idx_HEj)  &
        -k(221)*n(idx_HE)

    !d[C_dot]/d[E]
    pd(8,1) =  &
        +k(245)*n(idx_COj)  &
        +k(209)*n(idx_CHj)  &
        -k(157)*n(idx_C)  &
        +k(243)*n(idx_Cj)

    !d[CO_dot]/d[E]
    pd(10,1) =  &
        +k(249)*n(idx_HCOj)

    !d[N_dot]/d[E]
    pd(11,1) =  &
        +k(244)*n(idx_Nj)  &
        +k(3)*n(idx_NHj)

    !d[S_dot]/d[E]
    pd(13,1) =  &
        -k(137)*n(idx_S)  &
        +k(187)*n(idx_Sj)  &
        +k(248)*n(idx_HSj)

    !d[O_dot]/d[E]
    pd(15,1) =  &
        +k(74)*n(idx_Oj)  &
        +k(245)*n(idx_COj)  &
        -k(13)*n(idx_O)  &
        +k(115)*n(idx_SIOj)  &
        +k(87)*n(idx_OHj)

    !d[H2_dot]/d[E]
    pd(16,1) =  &
        -k(107)*n(idx_H2)  &
        -k(223)*n(idx_H2)

    !d[SI_dot]/d[E]
    pd(17,1) =  &
        +k(109)*n(idx_SIFj)  &
        +k(197)*n(idx_SIj)  &
        +k(201)*n(idx_SIHj)  &
        +k(115)*n(idx_SIOj)

    !d[FE_dot]/d[E]
    pd(22,1) =  &
        +k(142)*n(idx_FEj)

    !d[NA_dot]/d[E]
    pd(26,1) =  &
        +k(190)*n(idx_NAj)

    !d[F_dot]/d[E]
    pd(27,1) =  &
        +k(109)*n(idx_SIFj)

    !d[MG_dot]/d[E]
    pd(41,1) =  &
        +k(64)*n(idx_MGj)

    !d[P_dot]/d[E]
    pd(42,1) =  &
        +k(139)*n(idx_Pj)

    !d[HE+_dot]/d[E]
    pd(48,1) =  &
        +k(221)*n(idx_HE)  &
        -k(4)*n(idx_HEj)

    !d[H+_dot]/d[E]
    pd(49,1) =  &
        +k(220)*n(idx_H)  &
        -k(29)*n(idx_Hj)

    !d[NH+_dot]/d[E]
    pd(50,1) =  &
        -k(3)*n(idx_NHj)

    !d[HS+_dot]/d[E]
    pd(51,1) =  &
        -k(248)*n(idx_HSj)

    !d[S+_dot]/d[E]
    pd(52,1) =  &
        -k(187)*n(idx_Sj)

    !d[SI+_dot]/d[E]
    pd(53,1) =  &
        -k(197)*n(idx_SIj)

    !d[OH+_dot]/d[E]
    pd(54,1) =  &
        -k(87)*n(idx_OHj)

    !d[HEH+_dot]/d[E]
    pd(55,1) =  &
        -k(25)*n(idx_HEHj)

    !d[H2+_dot]/d[E]
    pd(56,1) =  &
        -k(102)*n(idx_H2j)

    !d[FE+_dot]/d[E]
    pd(57,1) =  &
        -k(142)*n(idx_FEj)

    !d[SIH+_dot]/d[E]
    pd(58,1) =  &
        -k(201)*n(idx_SIHj)

    !d[NA+_dot]/d[E]
    pd(59,1) =  &
        -k(190)*n(idx_NAj)

    !d[HCO+_dot]/d[E]
    pd(60,1) =  &
        -k(249)*n(idx_HCOj)

    !d[CH+_dot]/d[E]
    pd(61,1) =  &
        -k(209)*n(idx_CHj)

    !d[O+_dot]/d[E]
    pd(62,1) =  &
        -k(74)*n(idx_Oj)

    !d[MG+_dot]/d[E]
    pd(63,1) =  &
        -k(64)*n(idx_MGj)

    !d[SIO+_dot]/d[E]
    pd(64,1) =  &
        -k(115)*n(idx_SIOj)

    !d[P+_dot]/d[E]
    pd(65,1) =  &
        -k(139)*n(idx_Pj)

    !d[SIF+_dot]/d[E]
    pd(66,1) =  &
        -k(109)*n(idx_SIFj)

    !d[C+_dot]/d[E]
    pd(67,1) =  &
        -k(243)*n(idx_Cj)

    !d[N+_dot]/d[E]
    pd(68,1) =  &
        -k(244)*n(idx_Nj)

    !d[CO+_dot]/d[E]
    pd(69,1) =  &
        -k(245)*n(idx_COj)

    !d[Tgas_dot]/d[E]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(1)*1d-3
    if(dnn>0.d0) then
      nn(1) = n(1) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,1) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[O-]
    pd(1,2) =  &
        +k(50)*n(idx_H)

    !d[O-_dot]/d[O-]
    pd(2,2) =  &
        -k(210)*n(idx_FEj)  &
        -k(63)*n(idx_MGj)  &
        -k(199)*n(idx_Hj)  &
        -k(50)*n(idx_H)

    !d[H_dot]/d[O-]
    pd(6,2) =  &
        +k(199)*n(idx_Hj)  &
        -k(50)*n(idx_H)

    !d[O_dot]/d[O-]
    pd(15,2) =  &
        +k(199)*n(idx_Hj)  &
        +k(63)*n(idx_MGj)  &
        +k(210)*n(idx_FEj)

    !d[OH_dot]/d[O-]
    pd(18,2) =  &
        +k(50)*n(idx_H)

    !d[FE_dot]/d[O-]
    pd(22,2) =  &
        +k(210)*n(idx_FEj)

    !d[MG_dot]/d[O-]
    pd(41,2) =  &
        +k(63)*n(idx_MGj)

    !d[H+_dot]/d[O-]
    pd(49,2) =  &
        -k(199)*n(idx_Hj)

    !d[FE+_dot]/d[O-]
    pd(57,2) =  &
        -k(210)*n(idx_FEj)

    !d[MG+_dot]/d[O-]
    pd(63,2) =  &
        -k(63)*n(idx_MGj)

    !d[Tgas_dot]/d[O-]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(2)*1d-3
    if(dnn>0.d0) then
      nn(2) = n(2) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,2) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H-]
    pd(1,3) =  &
        -k(224)*n(idx_E)  &
        +k(226)*n(idx_Hj)  &
        +k(192)*n(idx_H)  &
        +k(129)*n(idx_C)  &
        +k(225)*n(idx_H)  &
        +2.d0*k(224)*n(idx_E)  &
        +k(69)*n(idx_N)  &
        +k(203)*n(idx_O)

    !d[H-_dot]/d[H-]
    pd(3,3) =  &
        -k(224)*n(idx_E)  &
        -k(203)*n(idx_O)  &
        -k(65)*n(idx_Hj)  &
        -k(69)*n(idx_N)  &
        -k(66)*n(idx_Sj)  &
        -k(214)*n(idx_MGj)  &
        -k(225)*n(idx_H)  &
        -k(129)*n(idx_C)  &
        -k(24)*n(idx_FEj)  &
        -k(145)*n(idx_SIj)  &
        -k(192)*n(idx_H)  &
        -k(226)*n(idx_Hj)  &
        -k(22)*n(idx_NAj)  &
        -k(73)*n(idx_Oj)

    !d[H_dot]/d[H-]
    pd(6,3) =  &
        +2.d0*k(65)*n(idx_Hj)  &
        +k(24)*n(idx_FEj)  &
        +k(73)*n(idx_Oj)  &
        +k(145)*n(idx_SIj)  &
        +k(66)*n(idx_Sj)  &
        -k(225)*n(idx_H)  &
        +2.d0*k(225)*n(idx_H)  &
        -k(192)*n(idx_H)  &
        +k(224)*n(idx_E)  &
        +k(214)*n(idx_MGj)  &
        +k(22)*n(idx_NAj)

    !d[C_dot]/d[H-]
    pd(8,3) =  &
        -k(129)*n(idx_C)

    !d[N_dot]/d[H-]
    pd(11,3) =  &
        -k(69)*n(idx_N)

    !d[S_dot]/d[H-]
    pd(13,3) =  &
        +k(66)*n(idx_Sj)

    !d[O_dot]/d[H-]
    pd(15,3) =  &
        +k(73)*n(idx_Oj)  &
        -k(203)*n(idx_O)

    !d[H2_dot]/d[H-]
    pd(16,3) =  &
        +k(192)*n(idx_H)

    !d[SI_dot]/d[H-]
    pd(17,3) =  &
        +k(145)*n(idx_SIj)

    !d[OH_dot]/d[H-]
    pd(18,3) =  &
        +k(203)*n(idx_O)

    !d[FE_dot]/d[H-]
    pd(22,3) =  &
        +k(24)*n(idx_FEj)

    !d[NA_dot]/d[H-]
    pd(26,3) =  &
        +k(22)*n(idx_NAj)

    !d[CH_dot]/d[H-]
    pd(29,3) =  &
        +k(129)*n(idx_C)

    !d[NH_dot]/d[H-]
    pd(34,3) =  &
        +k(69)*n(idx_N)

    !d[MG_dot]/d[H-]
    pd(41,3) =  &
        +k(214)*n(idx_MGj)

    !d[H+_dot]/d[H-]
    pd(49,3) =  &
        -k(226)*n(idx_Hj)  &
        -k(65)*n(idx_Hj)

    !d[S+_dot]/d[H-]
    pd(52,3) =  &
        -k(66)*n(idx_Sj)

    !d[SI+_dot]/d[H-]
    pd(53,3) =  &
        -k(145)*n(idx_SIj)

    !d[H2+_dot]/d[H-]
    pd(56,3) =  &
        +k(226)*n(idx_Hj)

    !d[FE+_dot]/d[H-]
    pd(57,3) =  &
        -k(24)*n(idx_FEj)

    !d[NA+_dot]/d[H-]
    pd(59,3) =  &
        -k(22)*n(idx_NAj)

    !d[O+_dot]/d[H-]
    pd(62,3) =  &
        -k(73)*n(idx_Oj)

    !d[MG+_dot]/d[H-]
    pd(63,3) =  &
        -k(214)*n(idx_MGj)

    !d[Tgas_dot]/d[H-]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(3)*1d-3
    if(dnn>0.d0) then
      nn(3) = n(3) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,3) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[S-]
    pd(1,4) =  &
        +k(76)*n(idx_H)

    !d[S-_dot]/d[S-]
    pd(4,4) =  &
        -k(76)*n(idx_H)

    !d[H_dot]/d[S-]
    pd(6,4) =  &
        -k(76)*n(idx_H)

    !d[HS_dot]/d[S-]
    pd(19,4) =  &
        +k(76)*n(idx_H)

    !d[Tgas_dot]/d[S-]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(4)*1d-3
    if(dnn>0.d0) then
      nn(4) = n(4) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,4) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[C-]
    pd(1,5) =  &
        +k(170)*n(idx_H)

    !d[C-_dot]/d[C-]
    pd(5,5) =  &
        -k(133)*n(idx_Hj)  &
        -k(170)*n(idx_H)

    !d[H_dot]/d[C-]
    pd(6,5) =  &
        +k(133)*n(idx_Hj)  &
        -k(170)*n(idx_H)

    !d[C_dot]/d[C-]
    pd(8,5) =  &
        +k(133)*n(idx_Hj)

    !d[CH_dot]/d[C-]
    pd(29,5) =  &
        +k(170)*n(idx_H)

    !d[H+_dot]/d[C-]
    pd(49,5) =  &
        -k(133)*n(idx_Hj)

    !d[Tgas_dot]/d[C-]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(5)*1d-3
    if(dnn>0.d0) then
      nn(5) = n(5) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,5) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H]
    pd(1,6) =  &
        +k(233)  &
        -k(220)*n(idx_E)  &
        +k(225)*n(idx_Hk)  &
        +k(170)*n(idx_Ck)  &
        +k(192)*n(idx_Hk)  &
        +k(50)*n(idx_Ok)  &
        -k(57)*n(idx_E)  &
        +2.d0*k(220)*n(idx_E)  &
        +k(76)*n(idx_Sk)

    !d[O-_dot]/d[H]
    pd(2,6) =  &
        -k(50)*n(idx_Ok)

    !d[H-_dot]/d[H]
    pd(3,6) =  &
        +k(57)*n(idx_E)  &
        -k(192)*n(idx_Hk)  &
        -k(225)*n(idx_Hk)

    !d[S-_dot]/d[H]
    pd(4,6) =  &
        -k(76)*n(idx_Sk)

    !d[C-_dot]/d[H]
    pd(5,6) =  &
        -k(170)*n(idx_Ck)

    !d[H_dot]/d[H]
    pd(6,6) =  &
        -k(220)*n(idx_E)  &
        -k(50)*n(idx_Ok)  &
        -k(28)*n(idx_O2)  &
        -k(135)*n(idx_OCS)  &
        +3.d0*k(33)*n(idx_H2)  &
        -k(30)*n(idx_NO)  &
        -k(11)*n(idx_HEHj)  &
        -k(76)*n(idx_Sk)  &
        -k(18)*n(idx_S2)  &
        -k(56)*n(idx_OCN)  &
        -k(110)*n(idx_SIHj)  &
        -k(163)*n(idx_H2S)  &
        -k(176)*n(idx_CO2)  &
        -k(93)*n(idx_O2)  &
        -k(36)*n(idx_CO)  &
        +2.d0*k(67)*n(idx_OH)  &
        -k(62)*n(idx_HS)  &
        -k(53)*n(idx_NH2)  &
        -k(217)*n(idx_NS)  &
        -k(194)*n(idx_CHj)  &
        -k(83)*n(idx_O)  &
        -k(38)*n(idx_CH2)  &
        -k(208)*n(idx_OH)  &
        -k(40)*n(idx_NH)  &
        -k(185)*n(idx_SO)  &
        -k(67)*n(idx_OH)  &
        -k(6)*n(idx_HSj)  &
        -k(70)*n(idx_NO)  &
        -k(153)*n(idx_C2)  &
        +k(28)*n(idx_O2)  &
        -k(225)*n(idx_Hk)  &
        -k(1)*n(idx_HEj)  &
        -k(150)*n(idx_OCN)  &
        +2.d0*k(108)*n(idx_CH)  &
        +3.d0*k(227)*n(idx_H)*n(idx_H)  &
        -k(196)*n(idx_H2j)  &
        -k(57)*n(idx_E)  &
        -9.d0*k(227)*n(idx_H)*n(idx_H)  &
        -k(95)*n(idx_HCO)  &
        -k(16)*n(idx_NS)  &
        -k(108)*n(idx_CH)  &
        -k(170)*n(idx_Ck)  &
        -4.d0*k(228)*n(idx_H2)*n(idx_H)  &
        -k(182)*n(idx_Oj)  &
        -k(88)*n(idx_C)  &
        +2.d0*k(225)*n(idx_Hk)  &
        -k(92)*n(idx_CH)  &
        -k(61)*n(idx_Hj)  &
        -k(43)*n(idx_HCN)  &
        -k(233)  &
        -k(33)*n(idx_H2)  &
        -k(103)*n(idx_H2O)  &
        -k(21)*n(idx_SO)  &
        -k(192)*n(idx_Hk)  &
        -4.d0*k(229)*n(idx_H)*n(idx_HE)  &
        -k(154)*n(idx_OCN)  &
        -k(15)*n(idx_SIj)  &
        -k(146)*n(idx_H2O)  &
        +2.d0*k(103)*n(idx_H2O)

    !d[HE_dot]/d[H]
    pd(7,6) =  &
        +2.d0*k(229)*n(idx_H)*n(idx_HE)  &
        +k(11)*n(idx_HEHj)  &
        +k(1)*n(idx_HEj)  &
        -2.d0*k(229)*n(idx_H)*n(idx_HE)

    !d[C_dot]/d[H]
    pd(8,6) =  &
        -k(88)*n(idx_C)  &
        +k(92)*n(idx_CH)  &
        +k(108)*n(idx_CH)  &
        +k(36)*n(idx_CO)  &
        +k(153)*n(idx_C2)

    !d[NO_dot]/d[H]
    pd(9,6) =  &
        -k(30)*n(idx_NO)  &
        -k(70)*n(idx_NO)

    !d[CO_dot]/d[H]
    pd(10,6) =  &
        +k(135)*n(idx_OCS)  &
        +k(154)*n(idx_OCN)  &
        -k(36)*n(idx_CO)  &
        +k(95)*n(idx_HCO)  &
        +k(176)*n(idx_CO2)

    !d[N_dot]/d[H]
    pd(11,6) =  &
        +k(40)*n(idx_NH)  &
        +k(30)*n(idx_NO)  &
        +k(16)*n(idx_NS)

    !d[O2_dot]/d[H]
    pd(12,6) =  &
        -k(93)*n(idx_O2)  &
        -k(28)*n(idx_O2)

    !d[S_dot]/d[H]
    pd(13,6) =  &
        +k(21)*n(idx_SO)  &
        +k(217)*n(idx_NS)  &
        +k(18)*n(idx_S2)  &
        +k(62)*n(idx_HS)

    !d[SO_dot]/d[H]
    pd(14,6) =  &
        -k(21)*n(idx_SO)  &
        -k(185)*n(idx_SO)

    !d[O_dot]/d[H]
    pd(15,6) =  &
        +k(150)*n(idx_OCN)  &
        +k(67)*n(idx_OH)  &
        +k(182)*n(idx_Oj)  &
        +2.d0*k(28)*n(idx_O2)  &
        +k(93)*n(idx_O2)  &
        +k(185)*n(idx_SO)  &
        +k(70)*n(idx_NO)  &
        -k(83)*n(idx_O)  &
        +k(208)*n(idx_OH)

    !d[H2_dot]/d[H]
    pd(16,6) =  &
        +4.d0*k(228)*n(idx_H2)*n(idx_H)  &
        +3.d0*k(227)*n(idx_H)*n(idx_H)  &
        +k(92)*n(idx_CH)  &
        +k(194)*n(idx_CHj)  &
        -2.d0*k(228)*n(idx_H2)*n(idx_H)  &
        +k(43)*n(idx_HCN)  &
        +k(40)*n(idx_NH)  &
        +k(163)*n(idx_H2S)  &
        +k(146)*n(idx_H2O)  &
        +k(6)*n(idx_HSj)  &
        -k(33)*n(idx_H2)  &
        +2.d0*k(229)*n(idx_H)*n(idx_HE)  &
        +k(192)*n(idx_Hk)  &
        +k(53)*n(idx_NH2)  &
        +k(62)*n(idx_HS)  &
        +k(196)*n(idx_H2j)  &
        +k(110)*n(idx_SIHj)  &
        +k(208)*n(idx_OH)  &
        +k(38)*n(idx_CH2)  &
        +k(95)*n(idx_HCO)

    !d[OH_dot]/d[H]
    pd(18,6) =  &
        +k(146)*n(idx_H2O)  &
        +k(56)*n(idx_OCN)  &
        +k(30)*n(idx_NO)  &
        +k(93)*n(idx_O2)  &
        -k(67)*n(idx_OH)  &
        +k(83)*n(idx_O)  &
        -k(208)*n(idx_OH)  &
        +k(36)*n(idx_CO)  &
        +k(176)*n(idx_CO2)  &
        +k(21)*n(idx_SO)  &
        +k(103)*n(idx_H2O)  &
        +k(50)*n(idx_Ok)

    !d[HS_dot]/d[H]
    pd(19,6) =  &
        +k(135)*n(idx_OCS)  &
        +k(163)*n(idx_H2S)  &
        +k(18)*n(idx_S2)  &
        -k(62)*n(idx_HS)  &
        +k(16)*n(idx_NS)  &
        +k(76)*n(idx_Sk)  &
        +k(185)*n(idx_SO)

    !d[NS_dot]/d[H]
    pd(20,6) =  &
        -k(16)*n(idx_NS)  &
        -k(217)*n(idx_NS)

    !d[H2S_dot]/d[H]
    pd(21,6) =  &
        -k(163)*n(idx_H2S)

    !d[CN_dot]/d[H]
    pd(24,6) =  &
        +k(56)*n(idx_OCN)  &
        +k(43)*n(idx_HCN)

    !d[S2_dot]/d[H]
    pd(25,6) =  &
        -k(18)*n(idx_S2)

    !d[CH_dot]/d[H]
    pd(29,6) =  &
        -k(108)*n(idx_CH)  &
        +k(170)*n(idx_Ck)  &
        -k(92)*n(idx_CH)  &
        +k(153)*n(idx_C2)  &
        +k(88)*n(idx_C)  &
        +k(38)*n(idx_CH2)

    !d[C2_dot]/d[H]
    pd(31,6) =  &
        -k(153)*n(idx_C2)

    !d[CH2_dot]/d[H]
    pd(33,6) =  &
        -k(38)*n(idx_CH2)

    !d[NH_dot]/d[H]
    pd(34,6) =  &
        +k(154)*n(idx_OCN)  &
        -k(40)*n(idx_NH)  &
        +k(70)*n(idx_NO)  &
        +k(53)*n(idx_NH2)  &
        +k(217)*n(idx_NS)

    !d[HCN_dot]/d[H]
    pd(35,6) =  &
        +k(150)*n(idx_OCN)  &
        -k(43)*n(idx_HCN)

    !d[CO2_dot]/d[H]
    pd(36,6) =  &
        -k(176)*n(idx_CO2)

    !d[NH2_dot]/d[H]
    pd(39,6) =  &
        -k(53)*n(idx_NH2)

    !d[OCN_dot]/d[H]
    pd(40,6) =  &
        -k(150)*n(idx_OCN)  &
        -k(154)*n(idx_OCN)  &
        -k(56)*n(idx_OCN)

    !d[HCO_dot]/d[H]
    pd(43,6) =  &
        -k(95)*n(idx_HCO)

    !d[H2O_dot]/d[H]
    pd(44,6) =  &
        -k(103)*n(idx_H2O)  &
        -k(146)*n(idx_H2O)

    !d[OCS_dot]/d[H]
    pd(45,6) =  &
        -k(135)*n(idx_OCS)

    !d[HE+_dot]/d[H]
    pd(48,6) =  &
        -k(1)*n(idx_HEj)

    !d[H+_dot]/d[H]
    pd(49,6) =  &
        +k(233)  &
        -k(61)*n(idx_Hj)  &
        +k(1)*n(idx_HEj)  &
        +k(182)*n(idx_Oj)  &
        +k(196)*n(idx_H2j)  &
        +k(220)*n(idx_E)

    !d[HS+_dot]/d[H]
    pd(51,6) =  &
        -k(6)*n(idx_HSj)

    !d[S+_dot]/d[H]
    pd(52,6) =  &
        +k(6)*n(idx_HSj)

    !d[SI+_dot]/d[H]
    pd(53,6) =  &
        +k(110)*n(idx_SIHj)  &
        -k(15)*n(idx_SIj)

    !d[HEH+_dot]/d[H]
    pd(55,6) =  &
        -k(11)*n(idx_HEHj)

    !d[H2+_dot]/d[H]
    pd(56,6) =  &
        -k(196)*n(idx_H2j)  &
        +k(11)*n(idx_HEHj)  &
        +k(61)*n(idx_Hj)

    !d[SIH+_dot]/d[H]
    pd(58,6) =  &
        +k(15)*n(idx_SIj)  &
        -k(110)*n(idx_SIHj)

    !d[CH+_dot]/d[H]
    pd(61,6) =  &
        -k(194)*n(idx_CHj)

    !d[O+_dot]/d[H]
    pd(62,6) =  &
        -k(182)*n(idx_Oj)

    !d[C+_dot]/d[H]
    pd(67,6) =  &
        +k(194)*n(idx_CHj)

    !d[Tgas_dot]/d[H]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(6)*1d-3
    if(dnn>0.d0) then
      nn(6) = n(6) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,6) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HE]
    pd(1,7) =  &
        -k(221)*n(idx_E)  &
        +2.d0*k(221)*n(idx_E)  &
        +k(238)

    !d[H_dot]/d[HE]
    pd(6,7) =  &
        -2.d0*k(229)*n(idx_H)*n(idx_H)  &
        +k(32)*n(idx_H2j)

    !d[HE_dot]/d[HE]
    pd(7,7) =  &
        -k(32)*n(idx_H2j)  &
        -k(221)*n(idx_E)  &
        +k(229)*n(idx_H)*n(idx_H)  &
        -k(238)  &
        -k(229)*n(idx_H)*n(idx_H)

    !d[H2_dot]/d[HE]
    pd(16,7) =  &
        +k(229)*n(idx_H)*n(idx_H)

    !d[HE+_dot]/d[HE]
    pd(48,7) =  &
        +k(221)*n(idx_E)  &
        +k(238)

    !d[HEH+_dot]/d[HE]
    pd(55,7) =  &
        +k(32)*n(idx_H2j)

    !d[H2+_dot]/d[HE]
    pd(56,7) =  &
        -k(32)*n(idx_H2j)

    !d[Tgas_dot]/d[HE]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(7)*1d-3
    if(dnn>0.d0) then
      nn(7) = n(7) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,7) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[C]
    pd(1,8) =  &
        -k(157)*n(idx_E)  &
        +k(129)*n(idx_Hk)  &
        +k(232)

    !d[H-_dot]/d[C]
    pd(3,8) =  &
        -k(129)*n(idx_Hk)

    !d[C-_dot]/d[C]
    pd(5,8) =  &
        +k(157)*n(idx_E)

    !d[H_dot]/d[C]
    pd(6,8) =  &
        +k(90)*n(idx_HS)  &
        +k(152)*n(idx_CH)  &
        +k(54)*n(idx_H2)  &
        +k(206)*n(idx_H2j)  &
        -k(88)*n(idx_H)  &
        +k(162)*n(idx_NH)  &
        +k(204)*n(idx_OH)

    !d[C_dot]/d[C]
    pd(8,8) =  &
        -k(105)*n(idx_NH)  &
        -k(157)*n(idx_E)  &
        -k(206)*n(idx_H2j)  &
        -k(129)*n(idx_Hk)  &
        -k(232)  &
        -k(75)*n(idx_H2)  &
        -k(204)*n(idx_OH)  &
        -k(162)*n(idx_NH)  &
        -k(54)*n(idx_H2)  &
        -k(52)*n(idx_HCOj)  &
        -k(152)*n(idx_CH)  &
        -k(34)*n(idx_CS)  &
        -k(198)*n(idx_NS)  &
        -k(47)*n(idx_CN)  &
        -k(189)*n(idx_SO2)  &
        -k(202)*n(idx_SIOj)  &
        -k(80)*n(idx_SO)  &
        -k(106)*n(idx_N2)  &
        -k(193)*n(idx_SO)  &
        -k(173)*n(idx_N)  &
        -k(88)*n(idx_H)  &
        -k(101)*n(idx_NO)  &
        -4.d0*k(166)*n(idx_C)  &
        -k(169)*n(idx_OH)  &
        -k(90)*n(idx_HS)  &
        -k(112)*n(idx_O)  &
        -k(171)*n(idx_HS)  &
        -k(2)*n(idx_NO)  &
        -k(216)*n(idx_O2)  &
        -k(46)*n(idx_CO)  &
        -k(207)*n(idx_S)

    !d[NO_dot]/d[C]
    pd(9,8) =  &
        -k(101)*n(idx_NO)  &
        -k(2)*n(idx_NO)

    !d[CO_dot]/d[C]
    pd(10,8) =  &
        +k(2)*n(idx_NO)  &
        +k(202)*n(idx_SIOj)  &
        +k(52)*n(idx_HCOj)  &
        +k(112)*n(idx_O)  &
        +k(193)*n(idx_SO)  &
        -k(46)*n(idx_CO)  &
        +k(189)*n(idx_SO2)  &
        +k(216)*n(idx_O2)  &
        +k(204)*n(idx_OH)

    !d[N_dot]/d[C]
    pd(11,8) =  &
        +k(106)*n(idx_N2)  &
        -k(173)*n(idx_N)  &
        +k(2)*n(idx_NO)  &
        +k(105)*n(idx_NH)  &
        +k(47)*n(idx_CN)

    !d[O2_dot]/d[C]
    pd(12,8) =  &
        -k(216)*n(idx_O2)

    !d[S_dot]/d[C]
    pd(13,8) =  &
        +k(171)*n(idx_HS)  &
        +k(193)*n(idx_SO)  &
        +k(34)*n(idx_CS)  &
        -k(207)*n(idx_S)  &
        +k(198)*n(idx_NS)

    !d[SO_dot]/d[C]
    pd(14,8) =  &
        -k(193)*n(idx_SO)  &
        -k(80)*n(idx_SO)  &
        +k(189)*n(idx_SO2)

    !d[O_dot]/d[C]
    pd(15,8) =  &
        +k(46)*n(idx_CO)  &
        -k(112)*n(idx_O)  &
        +k(216)*n(idx_O2)  &
        +k(101)*n(idx_NO)  &
        +k(80)*n(idx_SO)  &
        +k(169)*n(idx_OH)

    !d[H2_dot]/d[C]
    pd(16,8) =  &
        -k(75)*n(idx_H2)  &
        -k(54)*n(idx_H2)

    !d[OH_dot]/d[C]
    pd(18,8) =  &
        -k(169)*n(idx_OH)  &
        -k(204)*n(idx_OH)

    !d[HS_dot]/d[C]
    pd(19,8) =  &
        -k(171)*n(idx_HS)  &
        -k(90)*n(idx_HS)

    !d[NS_dot]/d[C]
    pd(20,8) =  &
        -k(198)*n(idx_NS)

    !d[CS_dot]/d[C]
    pd(23,8) =  &
        +k(207)*n(idx_S)  &
        +k(90)*n(idx_HS)  &
        +k(80)*n(idx_SO)  &
        -k(34)*n(idx_CS)

    !d[CN_dot]/d[C]
    pd(24,8) =  &
        +k(198)*n(idx_NS)  &
        -k(47)*n(idx_CN)  &
        +k(101)*n(idx_NO)  &
        +k(173)*n(idx_N)  &
        +k(106)*n(idx_N2)  &
        +k(162)*n(idx_NH)

    !d[CH_dot]/d[C]
    pd(29,8) =  &
        -k(152)*n(idx_CH)  &
        +k(105)*n(idx_NH)  &
        +k(54)*n(idx_H2)  &
        +k(88)*n(idx_H)  &
        +k(171)*n(idx_HS)  &
        +k(129)*n(idx_Hk)  &
        +k(169)*n(idx_OH)

    !d[SO2_dot]/d[C]
    pd(30,8) =  &
        -k(189)*n(idx_SO2)

    !d[C2_dot]/d[C]
    pd(31,8) =  &
        +k(152)*n(idx_CH)  &
        +2.d0*k(166)*n(idx_C)  &
        +k(46)*n(idx_CO)  &
        +k(34)*n(idx_CS)  &
        +k(47)*n(idx_CN)

    !d[N2_dot]/d[C]
    pd(32,8) =  &
        -k(106)*n(idx_N2)

    !d[CH2_dot]/d[C]
    pd(33,8) =  &
        +k(75)*n(idx_H2)

    !d[NH_dot]/d[C]
    pd(34,8) =  &
        -k(162)*n(idx_NH)  &
        -k(105)*n(idx_NH)

    !d[SI+_dot]/d[C]
    pd(53,8) =  &
        +k(202)*n(idx_SIOj)

    !d[H2+_dot]/d[C]
    pd(56,8) =  &
        -k(206)*n(idx_H2j)

    !d[HCO+_dot]/d[C]
    pd(60,8) =  &
        -k(52)*n(idx_HCOj)

    !d[CH+_dot]/d[C]
    pd(61,8) =  &
        +k(206)*n(idx_H2j)  &
        +k(52)*n(idx_HCOj)

    !d[SIO+_dot]/d[C]
    pd(64,8) =  &
        -k(202)*n(idx_SIOj)

    !d[C+_dot]/d[C]
    pd(67,8) =  &
        +k(232)

    !d[Tgas_dot]/d[C]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(8)*1d-3
    if(dnn>0.d0) then
      nn(8) = n(8) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,8) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[NO]
    pd(6,9) =  &
        -k(30)*n(idx_H)  &
        -k(70)*n(idx_H)

    !d[C_dot]/d[NO]
    pd(8,9) =  &
        -k(101)*n(idx_C)  &
        -k(2)*n(idx_C)

    !d[NO_dot]/d[NO]
    pd(9,9) =  &
        -k(30)*n(idx_H)  &
        -k(70)*n(idx_H)  &
        -k(123)*n(idx_SI)  &
        -k(101)*n(idx_C)  &
        -k(2)*n(idx_C)  &
        -k(35)*n(idx_N)

    !d[CO_dot]/d[NO]
    pd(10,9) =  &
        +k(2)*n(idx_C)

    !d[N_dot]/d[NO]
    pd(11,9) =  &
        +k(123)*n(idx_SI)  &
        +k(2)*n(idx_C)  &
        +k(30)*n(idx_H)  &
        -k(35)*n(idx_N)

    !d[O_dot]/d[NO]
    pd(15,9) =  &
        +k(35)*n(idx_N)  &
        +k(70)*n(idx_H)  &
        +k(101)*n(idx_C)

    !d[SI_dot]/d[NO]
    pd(17,9) =  &
        -k(123)*n(idx_SI)

    !d[OH_dot]/d[NO]
    pd(18,9) =  &
        +k(30)*n(idx_H)

    !d[CN_dot]/d[NO]
    pd(24,9) =  &
        +k(101)*n(idx_C)

    !d[N2_dot]/d[NO]
    pd(32,9) =  &
        +k(35)*n(idx_N)

    !d[NH_dot]/d[NO]
    pd(34,9) =  &
        +k(70)*n(idx_H)

    !d[SIO_dot]/d[NO]
    pd(37,9) =  &
        +k(123)*n(idx_SI)

    !d[Tgas_dot]/d[NO]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(9)*1d-3
    if(dnn>0.d0) then
      nn(9) = n(9) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,9) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[CO]
    pd(1,10) =  &
        +k(235)

    !d[H_dot]/d[CO]
    pd(6,10) =  &
        +k(44)*n(idx_OH)  &
        -k(36)*n(idx_H)

    !d[C_dot]/d[CO]
    pd(8,10) =  &
        +k(36)*n(idx_H)  &
        +k(242)  &
        +k(89)*n(idx_SI)  &
        -k(46)*n(idx_C)

    !d[CO_dot]/d[CO]
    pd(10,10) =  &
        -k(46)*n(idx_C)  &
        -k(44)*n(idx_OH)  &
        -k(36)*n(idx_H)  &
        -k(89)*n(idx_SI)  &
        -k(235)  &
        -k(242)

    !d[O_dot]/d[CO]
    pd(15,10) =  &
        +k(46)*n(idx_C)  &
        +k(242)

    !d[SI_dot]/d[CO]
    pd(17,10) =  &
        -k(89)*n(idx_SI)

    !d[OH_dot]/d[CO]
    pd(18,10) =  &
        +k(36)*n(idx_H)  &
        -k(44)*n(idx_OH)

    !d[C2_dot]/d[CO]
    pd(31,10) =  &
        +k(46)*n(idx_C)

    !d[CO2_dot]/d[CO]
    pd(36,10) =  &
        +k(44)*n(idx_OH)

    !d[SIO_dot]/d[CO]
    pd(37,10) =  &
        +k(89)*n(idx_SI)

    !d[CO+_dot]/d[CO]
    pd(69,10) =  &
        +k(235)

    !d[Tgas_dot]/d[CO]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(10)*1d-3
    if(dnn>0.d0) then
      nn(10) = n(10) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,10) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[N]
    pd(1,11) =  &
        +k(234)  &
        +k(69)*n(idx_Hk)

    !d[H-_dot]/d[N]
    pd(3,11) =  &
        -k(69)*n(idx_Hk)

    !d[H_dot]/d[N]
    pd(6,11) =  &
        +k(48)*n(idx_NH)  &
        +k(117)*n(idx_H2)  &
        +k(183)*n(idx_OH)  &
        +k(165)*n(idx_CH)  &
        +k(9)*n(idx_HS)

    !d[C_dot]/d[N]
    pd(8,11) =  &
        +k(151)*n(idx_C2)  &
        -k(173)*n(idx_C)  &
        +k(128)*n(idx_CN)  &
        +k(131)*n(idx_CH)

    !d[NO_dot]/d[N]
    pd(9,11) =  &
        +k(183)*n(idx_OH)  &
        +k(255)*n(idx_PO)  &
        +k(60)*n(idx_O2)  &
        +k(156)*n(idx_SIOj)  &
        +k(175)*n(idx_CO2)  &
        +k(211)*n(idx_SO)  &
        -k(35)*n(idx_NO)

    !d[CO_dot]/d[N]
    pd(10,11) =  &
        +k(175)*n(idx_CO2)

    !d[N_dot]/d[N]
    pd(11,11) =  &
        -k(173)*n(idx_C)  &
        -k(131)*n(idx_CH)  &
        -k(69)*n(idx_Hk)  &
        -k(165)*n(idx_CH)  &
        -k(17)*n(idx_CS)  &
        -k(35)*n(idx_NO)  &
        -k(60)*n(idx_O2)  &
        -k(128)*n(idx_CN)  &
        -k(234)  &
        -k(183)*n(idx_OH)  &
        -k(156)*n(idx_SIOj)  &
        -k(94)*n(idx_SO)  &
        -k(252)*n(idx_PN)  &
        -k(155)*n(idx_OH)  &
        -k(211)*n(idx_SO)  &
        -k(9)*n(idx_HS)  &
        -k(175)*n(idx_CO2)  &
        -k(255)*n(idx_PO)  &
        -k(179)*n(idx_HS)  &
        -k(253)*n(idx_PO)  &
        -k(117)*n(idx_H2)  &
        -k(48)*n(idx_NH)  &
        -k(151)*n(idx_C2)

    !d[O2_dot]/d[N]
    pd(12,11) =  &
        -k(60)*n(idx_O2)

    !d[S_dot]/d[N]
    pd(13,11) =  &
        +k(17)*n(idx_CS)  &
        +k(179)*n(idx_HS)  &
        +k(211)*n(idx_SO)

    !d[SO_dot]/d[N]
    pd(14,11) =  &
        -k(94)*n(idx_SO)  &
        -k(211)*n(idx_SO)

    !d[O_dot]/d[N]
    pd(15,11) =  &
        +k(35)*n(idx_NO)  &
        +k(155)*n(idx_OH)  &
        +k(253)*n(idx_PO)  &
        +k(94)*n(idx_SO)  &
        +k(60)*n(idx_O2)

    !d[H2_dot]/d[N]
    pd(16,11) =  &
        -k(117)*n(idx_H2)

    !d[OH_dot]/d[N]
    pd(18,11) =  &
        -k(183)*n(idx_OH)  &
        -k(155)*n(idx_OH)

    !d[HS_dot]/d[N]
    pd(19,11) =  &
        -k(9)*n(idx_HS)  &
        -k(179)*n(idx_HS)

    !d[NS_dot]/d[N]
    pd(20,11) =  &
        +k(94)*n(idx_SO)  &
        +k(9)*n(idx_HS)

    !d[CS_dot]/d[N]
    pd(23,11) =  &
        -k(17)*n(idx_CS)

    !d[CN_dot]/d[N]
    pd(24,11) =  &
        +k(17)*n(idx_CS)  &
        -k(128)*n(idx_CN)  &
        +k(151)*n(idx_C2)  &
        +k(173)*n(idx_C)  &
        +k(165)*n(idx_CH)

    !d[CH_dot]/d[N]
    pd(29,11) =  &
        -k(165)*n(idx_CH)  &
        -k(131)*n(idx_CH)

    !d[C2_dot]/d[N]
    pd(31,11) =  &
        -k(151)*n(idx_C2)

    !d[N2_dot]/d[N]
    pd(32,11) =  &
        +k(48)*n(idx_NH)  &
        +k(35)*n(idx_NO)  &
        +k(128)*n(idx_CN)  &
        +k(252)*n(idx_PN)

    !d[NH_dot]/d[N]
    pd(34,11) =  &
        +k(155)*n(idx_OH)  &
        +k(117)*n(idx_H2)  &
        +k(69)*n(idx_Hk)  &
        +k(179)*n(idx_HS)  &
        -k(48)*n(idx_NH)  &
        +k(131)*n(idx_CH)

    !d[CO2_dot]/d[N]
    pd(36,11) =  &
        -k(175)*n(idx_CO2)

    !d[P_dot]/d[N]
    pd(42,11) =  &
        +k(252)*n(idx_PN)  &
        +k(255)*n(idx_PO)

    !d[PN_dot]/d[N]
    pd(46,11) =  &
        +k(253)*n(idx_PO)  &
        -k(252)*n(idx_PN)

    !d[PO_dot]/d[N]
    pd(47,11) =  &
        -k(255)*n(idx_PO)  &
        -k(253)*n(idx_PO)

    !d[SI+_dot]/d[N]
    pd(53,11) =  &
        +k(156)*n(idx_SIOj)

    !d[SIO+_dot]/d[N]
    pd(64,11) =  &
        -k(156)*n(idx_SIOj)

    !d[N+_dot]/d[N]
    pd(68,11) =  &
        +k(234)

    !d[Tgas_dot]/d[N]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(11)*1d-3
    if(dnn>0.d0) then
      nn(11) = n(11) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,11) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[O2]
    pd(6,12) =  &
        +k(28)*n(idx_H)  &
        -k(28)*n(idx_H)  &
        -k(93)*n(idx_H)

    !d[C_dot]/d[O2]
    pd(8,12) =  &
        -k(216)*n(idx_C)

    !d[NO_dot]/d[O2]
    pd(9,12) =  &
        +k(60)*n(idx_N)

    !d[CO_dot]/d[O2]
    pd(10,12) =  &
        +k(216)*n(idx_C)

    !d[N_dot]/d[O2]
    pd(11,12) =  &
        -k(60)*n(idx_N)

    !d[O2_dot]/d[O2]
    pd(12,12) =  &
        -k(60)*n(idx_N)  &
        -k(26)*n(idx_H2)  &
        -k(113)*n(idx_CN)  &
        -k(28)*n(idx_H)  &
        -k(93)*n(idx_H)  &
        -k(49)*n(idx_SI)  &
        -k(5)*n(idx_S)  &
        -k(254)*n(idx_P)  &
        -k(216)*n(idx_C)

    !d[S_dot]/d[O2]
    pd(13,12) =  &
        -k(5)*n(idx_S)

    !d[SO_dot]/d[O2]
    pd(14,12) =  &
        +k(5)*n(idx_S)

    !d[O_dot]/d[O2]
    pd(15,12) =  &
        +k(216)*n(idx_C)  &
        +k(60)*n(idx_N)  &
        +k(49)*n(idx_SI)  &
        +k(5)*n(idx_S)  &
        +2.d0*k(28)*n(idx_H)  &
        +k(93)*n(idx_H)  &
        +k(254)*n(idx_P)  &
        +k(113)*n(idx_CN)

    !d[H2_dot]/d[O2]
    pd(16,12) =  &
        -k(26)*n(idx_H2)

    !d[SI_dot]/d[O2]
    pd(17,12) =  &
        -k(49)*n(idx_SI)

    !d[OH_dot]/d[O2]
    pd(18,12) =  &
        +2.d0*k(26)*n(idx_H2)  &
        +k(93)*n(idx_H)

    !d[CN_dot]/d[O2]
    pd(24,12) =  &
        -k(113)*n(idx_CN)

    !d[SIO_dot]/d[O2]
    pd(37,12) =  &
        +k(49)*n(idx_SI)

    !d[OCN_dot]/d[O2]
    pd(40,12) =  &
        +k(113)*n(idx_CN)

    !d[P_dot]/d[O2]
    pd(42,12) =  &
        -k(254)*n(idx_P)

    !d[PO_dot]/d[O2]
    pd(47,12) =  &
        +k(254)*n(idx_P)

    !d[Tgas_dot]/d[O2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(12)*1d-3
    if(dnn>0.d0) then
      nn(12) = n(12) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,12) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[S]
    pd(1,13) =  &
        -k(137)*n(idx_E)

    !d[S-_dot]/d[S]
    pd(4,13) =  &
        +k(137)*n(idx_E)

    !d[H_dot]/d[S]
    pd(6,13) =  &
        +k(23)*n(idx_CH)  &
        +k(167)*n(idx_Hj)  &
        +k(91)*n(idx_H2)  &
        +k(191)*n(idx_HS)  &
        +k(27)*n(idx_OH)  &
        +k(114)*n(idx_NH)

    !d[C_dot]/d[S]
    pd(8,13) =  &
        +k(177)*n(idx_C2)  &
        +k(85)*n(idx_CN)  &
        -k(207)*n(idx_C)  &
        +k(99)*n(idx_CH)

    !d[N_dot]/d[S]
    pd(11,13) =  &
        +k(158)*n(idx_NH)

    !d[O2_dot]/d[S]
    pd(12,13) =  &
        -k(5)*n(idx_O2)

    !d[S_dot]/d[S]
    pd(13,13) =  &
        -k(91)*n(idx_H2)  &
        -k(191)*n(idx_HS)  &
        -k(125)*n(idx_SO2)  &
        -k(207)*n(idx_C)  &
        -k(23)*n(idx_CH)  &
        -k(85)*n(idx_CN)  &
        -k(177)*n(idx_C2)  &
        -k(27)*n(idx_OH)  &
        -k(137)*n(idx_E)  &
        -k(114)*n(idx_NH)  &
        -k(167)*n(idx_Hj)  &
        -k(158)*n(idx_NH)  &
        -k(99)*n(idx_CH)  &
        -k(5)*n(idx_O2)

    !d[SO_dot]/d[S]
    pd(14,13) =  &
        +2.d0*k(125)*n(idx_SO2)  &
        +k(27)*n(idx_OH)  &
        +k(5)*n(idx_O2)

    !d[O_dot]/d[S]
    pd(15,13) =  &
        +k(5)*n(idx_O2)

    !d[H2_dot]/d[S]
    pd(16,13) =  &
        -k(91)*n(idx_H2)

    !d[OH_dot]/d[S]
    pd(18,13) =  &
        -k(27)*n(idx_OH)

    !d[HS_dot]/d[S]
    pd(19,13) =  &
        -k(191)*n(idx_HS)  &
        +k(158)*n(idx_NH)  &
        +k(91)*n(idx_H2)  &
        +k(99)*n(idx_CH)

    !d[NS_dot]/d[S]
    pd(20,13) =  &
        +k(85)*n(idx_CN)  &
        +k(114)*n(idx_NH)

    !d[CS_dot]/d[S]
    pd(23,13) =  &
        +k(23)*n(idx_CH)  &
        +k(177)*n(idx_C2)  &
        +k(207)*n(idx_C)

    !d[CN_dot]/d[S]
    pd(24,13) =  &
        -k(85)*n(idx_CN)

    !d[S2_dot]/d[S]
    pd(25,13) =  &
        +k(191)*n(idx_HS)

    !d[CH_dot]/d[S]
    pd(29,13) =  &
        -k(99)*n(idx_CH)  &
        -k(23)*n(idx_CH)

    !d[SO2_dot]/d[S]
    pd(30,13) =  &
        -k(125)*n(idx_SO2)

    !d[C2_dot]/d[S]
    pd(31,13) =  &
        -k(177)*n(idx_C2)

    !d[NH_dot]/d[S]
    pd(34,13) =  &
        -k(114)*n(idx_NH)  &
        -k(158)*n(idx_NH)

    !d[H+_dot]/d[S]
    pd(49,13) =  &
        -k(167)*n(idx_Hj)

    !d[S+_dot]/d[S]
    pd(52,13) =  &
        +k(167)*n(idx_Hj)

    !d[Tgas_dot]/d[S]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(13)*1d-3
    if(dnn>0.d0) then
      nn(13) = n(13) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,13) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[SO]
    pd(6,14) =  &
        +k(148)*n(idx_OH)  &
        -k(185)*n(idx_H)  &
        -k(21)*n(idx_H)

    !d[C_dot]/d[SO]
    pd(8,14) =  &
        -k(80)*n(idx_C)  &
        -k(193)*n(idx_C)

    !d[NO_dot]/d[SO]
    pd(9,14) =  &
        +k(211)*n(idx_N)

    !d[CO_dot]/d[SO]
    pd(10,14) =  &
        +k(193)*n(idx_C)

    !d[N_dot]/d[SO]
    pd(11,14) =  &
        -k(211)*n(idx_N)  &
        -k(94)*n(idx_N)

    !d[O2_dot]/d[SO]
    pd(12,14) =  &
        +k(122)*n(idx_O)

    !d[S_dot]/d[SO]
    pd(13,14) =  &
        +k(193)*n(idx_C)  &
        +k(122)*n(idx_O)  &
        +k(211)*n(idx_N)  &
        +k(21)*n(idx_H)

    !d[SO_dot]/d[SO]
    pd(14,14) =  &
        -k(21)*n(idx_H)  &
        -k(122)*n(idx_O)  &
        -k(193)*n(idx_C)  &
        -k(211)*n(idx_N)  &
        -k(80)*n(idx_C)  &
        -k(94)*n(idx_N)  &
        -k(185)*n(idx_H)  &
        -k(148)*n(idx_OH)

    !d[O_dot]/d[SO]
    pd(15,14) =  &
        +k(185)*n(idx_H)  &
        +k(94)*n(idx_N)  &
        -k(122)*n(idx_O)  &
        +k(80)*n(idx_C)

    !d[OH_dot]/d[SO]
    pd(18,14) =  &
        -k(148)*n(idx_OH)  &
        +k(21)*n(idx_H)

    !d[HS_dot]/d[SO]
    pd(19,14) =  &
        +k(185)*n(idx_H)

    !d[NS_dot]/d[SO]
    pd(20,14) =  &
        +k(94)*n(idx_N)

    !d[CS_dot]/d[SO]
    pd(23,14) =  &
        +k(80)*n(idx_C)

    !d[SO2_dot]/d[SO]
    pd(30,14) =  &
        +k(148)*n(idx_OH)

    !d[Tgas_dot]/d[SO]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(14)*1d-3
    if(dnn>0.d0) then
      nn(14) = n(14) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,14) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[O]
    pd(1,15) =  &
        +k(203)*n(idx_Hk)  &
        -k(13)*n(idx_E)  &
        +k(239)  &
        +k(200)*n(idx_CH)

    !d[O-_dot]/d[O]
    pd(2,15) =  &
        +k(13)*n(idx_E)

    !d[H-_dot]/d[O]
    pd(3,15) =  &
        -k(203)*n(idx_Hk)

    !d[H_dot]/d[O]
    pd(6,15) =  &
        +k(120)*n(idx_OH)  &
        +k(58)*n(idx_H2j)  &
        +k(45)*n(idx_NH)  &
        +k(181)*n(idx_CH)  &
        -k(83)*n(idx_H)  &
        +k(159)*n(idx_H2)  &
        +k(127)*n(idx_HCN)  &
        +k(184)*n(idx_Hj)  &
        +k(141)*n(idx_HS)

    !d[C_dot]/d[O]
    pd(8,15) =  &
        +k(168)*n(idx_C2)  &
        -k(112)*n(idx_C)  &
        +k(98)*n(idx_CH)  &
        +k(82)*n(idx_CS)  &
        +k(72)*n(idx_CN)

    !d[NO_dot]/d[O]
    pd(9,15) =  &
        +k(180)*n(idx_N2)  &
        +k(212)*n(idx_NS)  &
        +k(45)*n(idx_NH)  &
        +k(72)*n(idx_CN)

    !d[CO_dot]/d[O]
    pd(10,15) =  &
        +k(168)*n(idx_C2)  &
        +k(112)*n(idx_C)  &
        +k(116)*n(idx_HCN)  &
        +k(181)*n(idx_CH)  &
        +k(218)*n(idx_CN)  &
        +k(160)*n(idx_CS)

    !d[N_dot]/d[O]
    pd(11,15) =  &
        +k(180)*n(idx_N2)  &
        +k(218)*n(idx_CN)  &
        +k(41)*n(idx_NH)

    !d[O2_dot]/d[O]
    pd(12,15) =  &
        +k(31)*n(idx_SO2)  &
        +k(122)*n(idx_SO)  &
        +k(120)*n(idx_OH)  &
        +2.d0*k(39)*n(idx_O)  &
        +k(97)*n(idx_SIOj)

    !d[S_dot]/d[O]
    pd(13,15) =  &
        +k(160)*n(idx_CS)  &
        +k(212)*n(idx_NS)  &
        +k(122)*n(idx_SO)  &
        +k(86)*n(idx_HS)

    !d[SO_dot]/d[O]
    pd(14,15) =  &
        +k(31)*n(idx_SO2)  &
        -k(122)*n(idx_SO)  &
        +k(141)*n(idx_HS)  &
        +k(82)*n(idx_CS)

    !d[O_dot]/d[O]
    pd(15,15) =  &
        -k(82)*n(idx_CS)  &
        -k(86)*n(idx_HS)  &
        -k(31)*n(idx_SO2)  &
        -k(41)*n(idx_NH)  &
        -k(218)*n(idx_CN)  &
        -k(181)*n(idx_CH)  &
        -k(120)*n(idx_OH)  &
        -k(116)*n(idx_HCN)  &
        -k(180)*n(idx_N2)  &
        -k(127)*n(idx_HCN)  &
        -k(122)*n(idx_SO)  &
        -k(72)*n(idx_CN)  &
        -k(98)*n(idx_CH)  &
        -k(159)*n(idx_H2)  &
        -k(97)*n(idx_SIOj)  &
        -k(160)*n(idx_CS)  &
        -k(203)*n(idx_Hk)  &
        -k(186)*n(idx_SI)  &
        -k(13)*n(idx_E)  &
        -k(58)*n(idx_H2j)  &
        -k(141)*n(idx_HS)  &
        -k(200)*n(idx_CH)  &
        -k(168)*n(idx_C2)  &
        -k(212)*n(idx_NS)  &
        -k(45)*n(idx_NH)  &
        -k(184)*n(idx_Hj)  &
        -k(83)*n(idx_H)  &
        -k(112)*n(idx_C)  &
        -k(118)*n(idx_HCN)  &
        -k(239)  &
        -4.d0*k(39)*n(idx_O)  &
        -k(213)*n(idx_H2O)

    !d[H2_dot]/d[O]
    pd(16,15) =  &
        -k(159)*n(idx_H2)

    !d[SI_dot]/d[O]
    pd(17,15) =  &
        -k(186)*n(idx_SI)

    !d[OH_dot]/d[O]
    pd(18,15) =  &
        +2.d0*k(213)*n(idx_H2O)  &
        +k(98)*n(idx_CH)  &
        +k(118)*n(idx_HCN)  &
        +k(41)*n(idx_NH)  &
        +k(83)*n(idx_H)  &
        +k(203)*n(idx_Hk)  &
        -k(120)*n(idx_OH)  &
        +k(159)*n(idx_H2)  &
        +k(86)*n(idx_HS)

    !d[HS_dot]/d[O]
    pd(19,15) =  &
        -k(86)*n(idx_HS)  &
        -k(141)*n(idx_HS)

    !d[NS_dot]/d[O]
    pd(20,15) =  &
        -k(212)*n(idx_NS)

    !d[CS_dot]/d[O]
    pd(23,15) =  &
        -k(82)*n(idx_CS)  &
        -k(160)*n(idx_CS)

    !d[CN_dot]/d[O]
    pd(24,15) =  &
        -k(218)*n(idx_CN)  &
        -k(72)*n(idx_CN)  &
        +k(118)*n(idx_HCN)

    !d[CH_dot]/d[O]
    pd(29,15) =  &
        -k(200)*n(idx_CH)  &
        -k(181)*n(idx_CH)  &
        -k(98)*n(idx_CH)

    !d[SO2_dot]/d[O]
    pd(30,15) =  &
        -k(31)*n(idx_SO2)

    !d[C2_dot]/d[O]
    pd(31,15) =  &
        -k(168)*n(idx_C2)

    !d[N2_dot]/d[O]
    pd(32,15) =  &
        -k(180)*n(idx_N2)

    !d[NH_dot]/d[O]
    pd(34,15) =  &
        -k(41)*n(idx_NH)  &
        +k(116)*n(idx_HCN)  &
        -k(45)*n(idx_NH)

    !d[HCN_dot]/d[O]
    pd(35,15) =  &
        -k(116)*n(idx_HCN)  &
        -k(127)*n(idx_HCN)  &
        -k(118)*n(idx_HCN)

    !d[SIO_dot]/d[O]
    pd(37,15) =  &
        +k(186)*n(idx_SI)

    !d[OCN_dot]/d[O]
    pd(40,15) =  &
        +k(127)*n(idx_HCN)

    !d[H2O_dot]/d[O]
    pd(44,15) =  &
        -k(213)*n(idx_H2O)

    !d[H+_dot]/d[O]
    pd(49,15) =  &
        -k(184)*n(idx_Hj)

    !d[SI+_dot]/d[O]
    pd(53,15) =  &
        +k(97)*n(idx_SIOj)

    !d[OH+_dot]/d[O]
    pd(54,15) =  &
        +k(58)*n(idx_H2j)

    !d[H2+_dot]/d[O]
    pd(56,15) =  &
        -k(58)*n(idx_H2j)

    !d[HCO+_dot]/d[O]
    pd(60,15) =  &
        +k(200)*n(idx_CH)

    !d[O+_dot]/d[O]
    pd(62,15) =  &
        +k(184)*n(idx_Hj)  &
        +k(239)

    !d[SIO+_dot]/d[O]
    pd(64,15) =  &
        -k(97)*n(idx_SIOj)

    !d[Tgas_dot]/d[O]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(15)*1d-3
    if(dnn>0.d0) then
      nn(15) = n(15) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,15) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H2]
    pd(1,16) =  &
        +k(240)  &
        +k(236)  &
        +k(107)*n(idx_E)  &
        -k(223)*n(idx_E)  &
        -k(107)*n(idx_E)

    !d[H-_dot]/d[H2]
    pd(3,16) =  &
        +k(223)*n(idx_E)  &
        +k(231)

    !d[H_dot]/d[H2]
    pd(6,16) =  &
        +k(117)*n(idx_N)  &
        +k(215)*n(idx_CN)  &
        +k(37)*n(idx_F)  &
        +k(149)*n(idx_OH)  &
        +k(12)*n(idx_HS)  &
        +3.d0*k(33)*n(idx_H)  &
        +k(104)*n(idx_NH)  &
        -k(33)*n(idx_H)  &
        +k(55)*n(idx_Oj)  &
        +k(159)*n(idx_O)  &
        +k(42)*n(idx_CH)  &
        +k(178)*n(idx_Sj)  &
        +k(236)  &
        +k(222)*n(idx_Hj)  &
        +2.d0*k(237)  &
        +k(54)*n(idx_C)  &
        +2.d0*k(107)*n(idx_E)  &
        +k(91)*n(idx_S)  &
        -2.d0*k(228)*n(idx_H)*n(idx_H)  &
        +k(223)*n(idx_E)

    !d[C_dot]/d[H2]
    pd(8,16) =  &
        -k(75)*n(idx_C)  &
        -k(54)*n(idx_C)

    !d[N_dot]/d[H2]
    pd(11,16) =  &
        -k(117)*n(idx_N)

    !d[O2_dot]/d[H2]
    pd(12,16) =  &
        -k(26)*n(idx_O2)

    !d[S_dot]/d[H2]
    pd(13,16) =  &
        -k(91)*n(idx_S)

    !d[O_dot]/d[H2]
    pd(15,16) =  &
        -k(159)*n(idx_O)

    !d[H2_dot]/d[H2]
    pd(16,16) =  &
        -k(222)*n(idx_Hj)  &
        -k(37)*n(idx_F)  &
        -k(33)*n(idx_H)  &
        -k(42)*n(idx_CH)  &
        -k(228)*n(idx_H)*n(idx_H)  &
        -k(107)*n(idx_E)  &
        -k(104)*n(idx_NH)  &
        -k(236)  &
        -k(159)*n(idx_O)  &
        -k(178)*n(idx_Sj)  &
        -k(55)*n(idx_Oj)  &
        -k(12)*n(idx_HS)  &
        -k(149)*n(idx_OH)  &
        -k(251)*n(idx_Fj)  &
        -k(75)*n(idx_C)  &
        -k(117)*n(idx_N)  &
        -k(26)*n(idx_O2)  &
        -k(91)*n(idx_S)  &
        +2.d0*k(228)*n(idx_H)*n(idx_H)  &
        -k(237)  &
        -k(231)  &
        -k(240)  &
        -k(54)*n(idx_C)  &
        -k(215)*n(idx_CN)  &
        -k(223)*n(idx_E)

    !d[OH_dot]/d[H2]
    pd(18,16) =  &
        -k(149)*n(idx_OH)  &
        +2.d0*k(26)*n(idx_O2)  &
        +k(159)*n(idx_O)

    !d[HS_dot]/d[H2]
    pd(19,16) =  &
        +k(91)*n(idx_S)  &
        -k(12)*n(idx_HS)

    !d[H2S_dot]/d[H2]
    pd(21,16) =  &
        +k(12)*n(idx_HS)

    !d[CN_dot]/d[H2]
    pd(24,16) =  &
        -k(215)*n(idx_CN)

    !d[F_dot]/d[H2]
    pd(27,16) =  &
        -k(37)*n(idx_F)  &
        +k(251)*n(idx_Fj)

    !d[HF_dot]/d[H2]
    pd(28,16) =  &
        +k(37)*n(idx_F)

    !d[CH_dot]/d[H2]
    pd(29,16) =  &
        +k(54)*n(idx_C)  &
        -k(42)*n(idx_CH)

    !d[CH2_dot]/d[H2]
    pd(33,16) =  &
        +k(75)*n(idx_C)  &
        +k(42)*n(idx_CH)

    !d[NH_dot]/d[H2]
    pd(34,16) =  &
        -k(104)*n(idx_NH)  &
        +k(117)*n(idx_N)

    !d[HCN_dot]/d[H2]
    pd(35,16) =  &
        +k(215)*n(idx_CN)

    !d[NH2_dot]/d[H2]
    pd(39,16) =  &
        +k(104)*n(idx_NH)

    !d[H2O_dot]/d[H2]
    pd(44,16) =  &
        +k(149)*n(idx_OH)

    !d[H+_dot]/d[H2]
    pd(49,16) =  &
        +k(236)  &
        -k(222)*n(idx_Hj)  &
        +k(231)

    !d[HS+_dot]/d[H2]
    pd(51,16) =  &
        +k(178)*n(idx_Sj)

    !d[S+_dot]/d[H2]
    pd(52,16) =  &
        -k(178)*n(idx_Sj)

    !d[OH+_dot]/d[H2]
    pd(54,16) =  &
        +k(55)*n(idx_Oj)

    !d[H2+_dot]/d[H2]
    pd(56,16) =  &
        +k(240)  &
        +k(222)*n(idx_Hj)  &
        +k(251)*n(idx_Fj)

    !d[O+_dot]/d[H2]
    pd(62,16) =  &
        -k(55)*n(idx_Oj)

    !d[F+_dot]/d[H2]
    pd(70,16) =  &
        -k(251)*n(idx_Fj)

    !d[Tgas_dot]/d[H2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(16)*1d-3
    if(dnn>0.d0) then
      nn(16) = n(16) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,16) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[SI]
    pd(6,17) =  &
        +k(7)*n(idx_Hj)  &
        +k(205)*n(idx_OH)

    !d[HE_dot]/d[SI]
    pd(7,17) =  &
        +k(10)*n(idx_HEj)

    !d[C_dot]/d[SI]
    pd(8,17) =  &
        +k(121)*n(idx_Cj)  &
        +k(89)*n(idx_CO)

    !d[NO_dot]/d[SI]
    pd(9,17) =  &
        -k(123)*n(idx_NO)

    !d[CO_dot]/d[SI]
    pd(10,17) =  &
        +k(71)*n(idx_HCOj)  &
        +k(134)*n(idx_CO2)  &
        -k(89)*n(idx_CO)

    !d[N_dot]/d[SI]
    pd(11,17) =  &
        +k(123)*n(idx_NO)

    !d[O2_dot]/d[SI]
    pd(12,17) =  &
        -k(49)*n(idx_O2)

    !d[S_dot]/d[SI]
    pd(13,17) =  &
        +k(119)*n(idx_Sj)

    !d[O_dot]/d[SI]
    pd(15,17) =  &
        -k(186)*n(idx_O)  &
        +k(49)*n(idx_O2)

    !d[SI_dot]/d[SI]
    pd(17,17) =  &
        -k(71)*n(idx_HCOj)  &
        -k(7)*n(idx_Hj)  &
        -k(186)*n(idx_O)  &
        -k(143)*n(idx_Pj)  &
        -k(205)*n(idx_OH)  &
        -k(10)*n(idx_HEj)  &
        -k(121)*n(idx_Cj)  &
        -k(119)*n(idx_Sj)  &
        -k(49)*n(idx_O2)  &
        -k(134)*n(idx_CO2)  &
        -k(89)*n(idx_CO)  &
        -k(123)*n(idx_NO)

    !d[OH_dot]/d[SI]
    pd(18,17) =  &
        -k(205)*n(idx_OH)

    !d[CO2_dot]/d[SI]
    pd(36,17) =  &
        -k(134)*n(idx_CO2)

    !d[SIO_dot]/d[SI]
    pd(37,17) =  &
        +k(49)*n(idx_O2)  &
        +k(123)*n(idx_NO)  &
        +k(205)*n(idx_OH)  &
        +k(89)*n(idx_CO)  &
        +k(186)*n(idx_O)  &
        +k(134)*n(idx_CO2)

    !d[P_dot]/d[SI]
    pd(42,17) =  &
        +k(143)*n(idx_Pj)

    !d[HE+_dot]/d[SI]
    pd(48,17) =  &
        -k(10)*n(idx_HEj)

    !d[H+_dot]/d[SI]
    pd(49,17) =  &
        -k(7)*n(idx_Hj)

    !d[S+_dot]/d[SI]
    pd(52,17) =  &
        -k(119)*n(idx_Sj)

    !d[SI+_dot]/d[SI]
    pd(53,17) =  &
        +k(121)*n(idx_Cj)  &
        +k(7)*n(idx_Hj)  &
        +k(143)*n(idx_Pj)  &
        +k(10)*n(idx_HEj)  &
        +k(119)*n(idx_Sj)

    !d[SIH+_dot]/d[SI]
    pd(58,17) =  &
        +k(71)*n(idx_HCOj)

    !d[HCO+_dot]/d[SI]
    pd(60,17) =  &
        -k(71)*n(idx_HCOj)

    !d[P+_dot]/d[SI]
    pd(65,17) =  &
        -k(143)*n(idx_Pj)

    !d[C+_dot]/d[SI]
    pd(67,17) =  &
        -k(121)*n(idx_Cj)

    !d[Tgas_dot]/d[SI]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(17)*1d-3
    if(dnn>0.d0) then
      nn(17) = n(17) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,17) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[OH]
    pd(6,18) =  &
        +k(183)*n(idx_N)  &
        +k(120)*n(idx_O)  &
        +k(44)*n(idx_CO)  &
        +k(148)*n(idx_SO)  &
        +k(174)*n(idx_CS)  &
        +2.d0*k(67)*n(idx_H)  &
        +k(219)*n(idx_CN)  &
        -k(208)*n(idx_H)  &
        +k(81)*n(idx_SIj)  &
        +k(8)*n(idx_Hj)  &
        +k(51)*n(idx_SIO)  &
        +k(149)*n(idx_H2)  &
        +k(27)*n(idx_S)  &
        -k(67)*n(idx_H)  &
        +k(205)*n(idx_SI)  &
        +k(204)*n(idx_C)

    !d[C_dot]/d[OH]
    pd(8,18) =  &
        -k(169)*n(idx_C)  &
        -k(204)*n(idx_C)

    !d[NO_dot]/d[OH]
    pd(9,18) =  &
        +k(183)*n(idx_N)

    !d[CO_dot]/d[OH]
    pd(10,18) =  &
        -k(44)*n(idx_CO)  &
        +k(77)*n(idx_CS)  &
        +k(204)*n(idx_C)

    !d[N_dot]/d[OH]
    pd(11,18) =  &
        -k(155)*n(idx_N)  &
        -k(183)*n(idx_N)

    !d[O2_dot]/d[OH]
    pd(12,18) =  &
        +k(120)*n(idx_O)

    !d[S_dot]/d[OH]
    pd(13,18) =  &
        -k(27)*n(idx_S)

    !d[SO_dot]/d[OH]
    pd(14,18) =  &
        -k(148)*n(idx_SO)  &
        +k(27)*n(idx_S)

    !d[O_dot]/d[OH]
    pd(15,18) =  &
        -k(120)*n(idx_O)  &
        +k(155)*n(idx_N)  &
        +2.d0*k(147)*n(idx_OH)  &
        +k(20)*n(idx_F)  &
        +k(67)*n(idx_H)  &
        +k(96)*n(idx_CN)  &
        +k(208)*n(idx_H)  &
        +k(169)*n(idx_C)

    !d[H2_dot]/d[OH]
    pd(16,18) =  &
        +k(208)*n(idx_H)  &
        -k(149)*n(idx_H2)

    !d[SI_dot]/d[OH]
    pd(17,18) =  &
        -k(205)*n(idx_SI)

    !d[OH_dot]/d[OH]
    pd(18,18) =  &
        -k(120)*n(idx_O)  &
        -k(149)*n(idx_H2)  &
        -k(27)*n(idx_S)  &
        -k(205)*n(idx_SI)  &
        -k(44)*n(idx_CO)  &
        -k(67)*n(idx_H)  &
        -k(169)*n(idx_C)  &
        -k(183)*n(idx_N)  &
        -k(208)*n(idx_H)  &
        -k(155)*n(idx_N)  &
        -k(96)*n(idx_CN)  &
        -4.d0*k(147)*n(idx_OH)  &
        -k(174)*n(idx_CS)  &
        -k(219)*n(idx_CN)  &
        -k(77)*n(idx_CS)  &
        -k(8)*n(idx_Hj)  &
        -k(161)*n(idx_H2S)  &
        -k(81)*n(idx_SIj)  &
        -k(20)*n(idx_F)  &
        -k(51)*n(idx_SIO)  &
        -k(204)*n(idx_C)  &
        -k(148)*n(idx_SO)

    !d[HS_dot]/d[OH]
    pd(19,18) =  &
        +k(77)*n(idx_CS)  &
        +k(161)*n(idx_H2S)

    !d[H2S_dot]/d[OH]
    pd(21,18) =  &
        -k(161)*n(idx_H2S)

    !d[CS_dot]/d[OH]
    pd(23,18) =  &
        -k(77)*n(idx_CS)  &
        -k(174)*n(idx_CS)

    !d[CN_dot]/d[OH]
    pd(24,18) =  &
        -k(96)*n(idx_CN)  &
        -k(219)*n(idx_CN)

    !d[F_dot]/d[OH]
    pd(27,18) =  &
        -k(20)*n(idx_F)

    !d[HF_dot]/d[OH]
    pd(28,18) =  &
        +k(20)*n(idx_F)

    !d[CH_dot]/d[OH]
    pd(29,18) =  &
        +k(169)*n(idx_C)

    !d[SO2_dot]/d[OH]
    pd(30,18) =  &
        +k(148)*n(idx_SO)

    !d[NH_dot]/d[OH]
    pd(34,18) =  &
        +k(155)*n(idx_N)

    !d[HCN_dot]/d[OH]
    pd(35,18) =  &
        +k(96)*n(idx_CN)

    !d[CO2_dot]/d[OH]
    pd(36,18) =  &
        +k(44)*n(idx_CO)

    !d[SIO_dot]/d[OH]
    pd(37,18) =  &
        -k(51)*n(idx_SIO)  &
        +k(205)*n(idx_SI)

    !d[SIO2_dot]/d[OH]
    pd(38,18) =  &
        +k(51)*n(idx_SIO)

    !d[OCN_dot]/d[OH]
    pd(40,18) =  &
        +k(219)*n(idx_CN)

    !d[H2O_dot]/d[OH]
    pd(44,18) =  &
        +k(149)*n(idx_H2)  &
        +2.d0*k(147)*n(idx_OH)  &
        +k(161)*n(idx_H2S)

    !d[OCS_dot]/d[OH]
    pd(45,18) =  &
        +k(174)*n(idx_CS)

    !d[H+_dot]/d[OH]
    pd(49,18) =  &
        -k(8)*n(idx_Hj)

    !d[SI+_dot]/d[OH]
    pd(53,18) =  &
        -k(81)*n(idx_SIj)

    !d[OH+_dot]/d[OH]
    pd(54,18) =  &
        +k(8)*n(idx_Hj)

    !d[SIO+_dot]/d[OH]
    pd(64,18) =  &
        +k(81)*n(idx_SIj)

    !d[Tgas_dot]/d[OH]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(18)*1d-3
    if(dnn>0.d0) then
      nn(18) = n(18) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,18) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[HS]
    pd(6,19) =  &
        +k(191)*n(idx_S)  &
        -k(62)*n(idx_H)  &
        +k(9)*n(idx_N)  &
        +k(141)*n(idx_O)  &
        +k(12)*n(idx_H2)  &
        +k(90)*n(idx_C)

    !d[C_dot]/d[HS]
    pd(8,19) =  &
        -k(90)*n(idx_C)  &
        -k(171)*n(idx_C)

    !d[N_dot]/d[HS]
    pd(11,19) =  &
        -k(179)*n(idx_N)  &
        -k(9)*n(idx_N)

    !d[S_dot]/d[HS]
    pd(13,19) =  &
        -k(191)*n(idx_S)  &
        +k(62)*n(idx_H)  &
        +k(171)*n(idx_C)  &
        +k(86)*n(idx_O)  &
        +2.d0*k(188)*n(idx_HS)  &
        +k(179)*n(idx_N)

    !d[SO_dot]/d[HS]
    pd(14,19) =  &
        +k(141)*n(idx_O)

    !d[O_dot]/d[HS]
    pd(15,19) =  &
        -k(86)*n(idx_O)  &
        -k(141)*n(idx_O)

    !d[H2_dot]/d[HS]
    pd(16,19) =  &
        -k(12)*n(idx_H2)  &
        +k(62)*n(idx_H)

    !d[OH_dot]/d[HS]
    pd(18,19) =  &
        +k(86)*n(idx_O)

    !d[HS_dot]/d[HS]
    pd(19,19) =  &
        -k(86)*n(idx_O)  &
        -k(191)*n(idx_S)  &
        -k(90)*n(idx_C)  &
        -k(62)*n(idx_H)  &
        -k(171)*n(idx_C)  &
        -k(141)*n(idx_O)  &
        -k(179)*n(idx_N)  &
        -k(9)*n(idx_N)  &
        -k(12)*n(idx_H2)  &
        -4.d0*k(188)*n(idx_HS)

    !d[NS_dot]/d[HS]
    pd(20,19) =  &
        +k(9)*n(idx_N)

    !d[H2S_dot]/d[HS]
    pd(21,19) =  &
        +2.d0*k(188)*n(idx_HS)  &
        +k(12)*n(idx_H2)

    !d[CS_dot]/d[HS]
    pd(23,19) =  &
        +k(90)*n(idx_C)

    !d[S2_dot]/d[HS]
    pd(25,19) =  &
        +k(191)*n(idx_S)

    !d[CH_dot]/d[HS]
    pd(29,19) =  &
        +k(171)*n(idx_C)

    !d[NH_dot]/d[HS]
    pd(34,19) =  &
        +k(179)*n(idx_N)

    !d[Tgas_dot]/d[HS]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(19)*1d-3
    if(dnn>0.d0) then
      nn(19) = n(19) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,19) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[NS]
    pd(6,20) =  &
        -k(217)*n(idx_H)  &
        -k(16)*n(idx_H)

    !d[C_dot]/d[NS]
    pd(8,20) =  &
        -k(198)*n(idx_C)

    !d[NO_dot]/d[NS]
    pd(9,20) =  &
        +k(212)*n(idx_O)

    !d[N_dot]/d[NS]
    pd(11,20) =  &
        +k(16)*n(idx_H)

    !d[S_dot]/d[NS]
    pd(13,20) =  &
        +k(212)*n(idx_O)  &
        +k(198)*n(idx_C)  &
        +k(217)*n(idx_H)

    !d[O_dot]/d[NS]
    pd(15,20) =  &
        -k(212)*n(idx_O)

    !d[HS_dot]/d[NS]
    pd(19,20) =  &
        +k(16)*n(idx_H)

    !d[NS_dot]/d[NS]
    pd(20,20) =  &
        -k(217)*n(idx_H)  &
        -k(212)*n(idx_O)  &
        -k(16)*n(idx_H)  &
        -k(198)*n(idx_C)

    !d[CN_dot]/d[NS]
    pd(24,20) =  &
        +k(198)*n(idx_C)

    !d[NH_dot]/d[NS]
    pd(34,20) =  &
        +k(217)*n(idx_H)

    !d[Tgas_dot]/d[NS]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(20)*1d-3
    if(dnn>0.d0) then
      nn(20) = n(20) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,20) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[H2S]
    pd(6,21) =  &
        -k(163)*n(idx_H)

    !d[H2_dot]/d[H2S]
    pd(16,21) =  &
        +k(163)*n(idx_H)

    !d[OH_dot]/d[H2S]
    pd(18,21) =  &
        -k(161)*n(idx_OH)

    !d[HS_dot]/d[H2S]
    pd(19,21) =  &
        +k(163)*n(idx_H)  &
        +k(161)*n(idx_OH)

    !d[H2S_dot]/d[H2S]
    pd(21,21) =  &
        -k(161)*n(idx_OH)  &
        -k(163)*n(idx_H)

    !d[H2O_dot]/d[H2S]
    pd(44,21) =  &
        +k(161)*n(idx_OH)

    !d[Tgas_dot]/d[H2S]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(21)*1d-3
    if(dnn>0.d0) then
      nn(21) = n(21) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,21) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[FE]
    pd(6,22) =  &
        +k(172)*n(idx_Hj)

    !d[C_dot]/d[FE]
    pd(8,22) =  &
        +k(124)*n(idx_Cj)

    !d[S_dot]/d[FE]
    pd(13,22) =  &
        +k(14)*n(idx_Sj)

    !d[O_dot]/d[FE]
    pd(15,22) =  &
        +k(68)*n(idx_Oj)

    !d[SI_dot]/d[FE]
    pd(17,22) =  &
        +k(136)*n(idx_SIj)

    !d[FE_dot]/d[FE]
    pd(22,22) =  &
        -k(79)*n(idx_SIOj)  &
        -k(124)*n(idx_Cj)  &
        -k(68)*n(idx_Oj)  &
        -k(138)*n(idx_HCOj)  &
        -k(14)*n(idx_Sj)  &
        -k(136)*n(idx_SIj)  &
        -k(172)*n(idx_Hj)

    !d[SIO_dot]/d[FE]
    pd(37,22) =  &
        +k(79)*n(idx_SIOj)

    !d[HCO_dot]/d[FE]
    pd(43,22) =  &
        +k(138)*n(idx_HCOj)

    !d[H+_dot]/d[FE]
    pd(49,22) =  &
        -k(172)*n(idx_Hj)

    !d[S+_dot]/d[FE]
    pd(52,22) =  &
        -k(14)*n(idx_Sj)

    !d[SI+_dot]/d[FE]
    pd(53,22) =  &
        -k(136)*n(idx_SIj)

    !d[FE+_dot]/d[FE]
    pd(57,22) =  &
        +k(79)*n(idx_SIOj)  &
        +k(138)*n(idx_HCOj)  &
        +k(136)*n(idx_SIj)  &
        +k(124)*n(idx_Cj)  &
        +k(14)*n(idx_Sj)  &
        +k(68)*n(idx_Oj)  &
        +k(172)*n(idx_Hj)

    !d[HCO+_dot]/d[FE]
    pd(60,22) =  &
        -k(138)*n(idx_HCOj)

    !d[O+_dot]/d[FE]
    pd(62,22) =  &
        -k(68)*n(idx_Oj)

    !d[SIO+_dot]/d[FE]
    pd(64,22) =  &
        -k(79)*n(idx_SIOj)

    !d[C+_dot]/d[FE]
    pd(67,22) =  &
        -k(124)*n(idx_Cj)

    !d[Tgas_dot]/d[FE]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(22)*1d-3
    if(dnn>0.d0) then
      nn(22) = n(22) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,22) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[CS]
    pd(6,23) =  &
        +k(174)*n(idx_OH)

    !d[C_dot]/d[CS]
    pd(8,23) =  &
        -k(34)*n(idx_C)  &
        +k(82)*n(idx_O)

    !d[CO_dot]/d[CS]
    pd(10,23) =  &
        +k(77)*n(idx_OH)  &
        +k(160)*n(idx_O)

    !d[N_dot]/d[CS]
    pd(11,23) =  &
        -k(17)*n(idx_N)

    !d[S_dot]/d[CS]
    pd(13,23) =  &
        +k(17)*n(idx_N)  &
        +k(34)*n(idx_C)  &
        +k(160)*n(idx_O)

    !d[SO_dot]/d[CS]
    pd(14,23) =  &
        +k(82)*n(idx_O)

    !d[O_dot]/d[CS]
    pd(15,23) =  &
        -k(82)*n(idx_O)  &
        -k(160)*n(idx_O)

    !d[OH_dot]/d[CS]
    pd(18,23) =  &
        -k(77)*n(idx_OH)  &
        -k(174)*n(idx_OH)

    !d[HS_dot]/d[CS]
    pd(19,23) =  &
        +k(77)*n(idx_OH)

    !d[CS_dot]/d[CS]
    pd(23,23) =  &
        -k(82)*n(idx_O)  &
        -k(160)*n(idx_O)  &
        -k(174)*n(idx_OH)  &
        -k(34)*n(idx_C)  &
        -k(77)*n(idx_OH)  &
        -k(17)*n(idx_N)

    !d[CN_dot]/d[CS]
    pd(24,23) =  &
        +k(17)*n(idx_N)

    !d[C2_dot]/d[CS]
    pd(31,23) =  &
        +k(34)*n(idx_C)

    !d[OCS_dot]/d[CS]
    pd(45,23) =  &
        +k(174)*n(idx_OH)

    !d[Tgas_dot]/d[CS]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(23)*1d-3
    if(dnn>0.d0) then
      nn(23) = n(23) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,23) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[CN]
    pd(6,24) =  &
        +k(219)*n(idx_OH)  &
        +k(215)*n(idx_H2)

    !d[C_dot]/d[CN]
    pd(8,24) =  &
        +k(128)*n(idx_N)  &
        +k(72)*n(idx_O)  &
        -k(47)*n(idx_C)  &
        +k(85)*n(idx_S)

    !d[NO_dot]/d[CN]
    pd(9,24) =  &
        +k(72)*n(idx_O)

    !d[CO_dot]/d[CN]
    pd(10,24) =  &
        +k(218)*n(idx_O)

    !d[N_dot]/d[CN]
    pd(11,24) =  &
        +k(218)*n(idx_O)  &
        +k(47)*n(idx_C)  &
        -k(128)*n(idx_N)

    !d[O2_dot]/d[CN]
    pd(12,24) =  &
        -k(113)*n(idx_O2)

    !d[S_dot]/d[CN]
    pd(13,24) =  &
        -k(85)*n(idx_S)

    !d[O_dot]/d[CN]
    pd(15,24) =  &
        -k(218)*n(idx_O)  &
        +k(96)*n(idx_OH)  &
        -k(72)*n(idx_O)  &
        +k(113)*n(idx_O2)

    !d[H2_dot]/d[CN]
    pd(16,24) =  &
        -k(215)*n(idx_H2)

    !d[OH_dot]/d[CN]
    pd(18,24) =  &
        -k(96)*n(idx_OH)  &
        -k(219)*n(idx_OH)

    !d[NS_dot]/d[CN]
    pd(20,24) =  &
        +k(85)*n(idx_S)

    !d[CN_dot]/d[CN]
    pd(24,24) =  &
        -k(113)*n(idx_O2)  &
        -k(128)*n(idx_N)  &
        -k(85)*n(idx_S)  &
        -k(47)*n(idx_C)  &
        -k(219)*n(idx_OH)  &
        -k(215)*n(idx_H2)  &
        -k(96)*n(idx_OH)  &
        -k(218)*n(idx_O)  &
        -k(72)*n(idx_O)

    !d[C2_dot]/d[CN]
    pd(31,24) =  &
        +k(47)*n(idx_C)

    !d[N2_dot]/d[CN]
    pd(32,24) =  &
        +k(128)*n(idx_N)

    !d[HCN_dot]/d[CN]
    pd(35,24) =  &
        +k(215)*n(idx_H2)  &
        +k(96)*n(idx_OH)

    !d[OCN_dot]/d[CN]
    pd(40,24) =  &
        +k(219)*n(idx_OH)  &
        +k(113)*n(idx_O2)

    !d[Tgas_dot]/d[CN]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(24)*1d-3
    if(dnn>0.d0) then
      nn(24) = n(24) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,24) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[S2]
    pd(6,25) =  &
        -k(18)*n(idx_H)

    !d[S_dot]/d[S2]
    pd(13,25) =  &
        +k(18)*n(idx_H)

    !d[HS_dot]/d[S2]
    pd(19,25) =  &
        +k(18)*n(idx_H)

    !d[S2_dot]/d[S2]
    pd(25,25) =  &
        -k(18)*n(idx_H)

    !d[Tgas_dot]/d[S2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(25)*1d-3
    if(dnn>0.d0) then
      nn(25) = n(25) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,25) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[NA]
    pd(6,26) =  &
        +k(230)*n(idx_Hj)

    !d[S_dot]/d[NA]
    pd(13,26) =  &
        +k(19)*n(idx_Sj)

    !d[SI_dot]/d[NA]
    pd(17,26) =  &
        +k(132)*n(idx_SIj)

    !d[FE_dot]/d[NA]
    pd(22,26) =  &
        +k(78)*n(idx_FEj)

    !d[NA_dot]/d[NA]
    pd(26,26) =  &
        -k(132)*n(idx_SIj)  &
        -k(195)*n(idx_MGj)  &
        -k(78)*n(idx_FEj)  &
        -k(230)*n(idx_Hj)  &
        -k(19)*n(idx_Sj)

    !d[MG_dot]/d[NA]
    pd(41,26) =  &
        +k(195)*n(idx_MGj)

    !d[H+_dot]/d[NA]
    pd(49,26) =  &
        -k(230)*n(idx_Hj)

    !d[S+_dot]/d[NA]
    pd(52,26) =  &
        -k(19)*n(idx_Sj)

    !d[SI+_dot]/d[NA]
    pd(53,26) =  &
        -k(132)*n(idx_SIj)

    !d[FE+_dot]/d[NA]
    pd(57,26) =  &
        -k(78)*n(idx_FEj)

    !d[NA+_dot]/d[NA]
    pd(59,26) =  &
        +k(230)*n(idx_Hj)  &
        +k(19)*n(idx_Sj)  &
        +k(195)*n(idx_MGj)  &
        +k(132)*n(idx_SIj)  &
        +k(78)*n(idx_FEj)

    !d[MG+_dot]/d[NA]
    pd(63,26) =  &
        -k(195)*n(idx_MGj)

    !d[Tgas_dot]/d[NA]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(26)*1d-3
    if(dnn>0.d0) then
      nn(26) = n(26) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,26) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[F]
    pd(6,27) =  &
        +k(37)*n(idx_H2)

    !d[O_dot]/d[F]
    pd(15,27) =  &
        +k(20)*n(idx_OH)

    !d[H2_dot]/d[F]
    pd(16,27) =  &
        -k(37)*n(idx_H2)

    !d[OH_dot]/d[F]
    pd(18,27) =  &
        -k(20)*n(idx_OH)

    !d[F_dot]/d[F]
    pd(27,27) =  &
        -k(20)*n(idx_OH)  &
        -k(37)*n(idx_H2)

    !d[HF_dot]/d[F]
    pd(28,27) =  &
        +k(37)*n(idx_H2)  &
        +k(20)*n(idx_OH)

    !d[Tgas_dot]/d[F]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(27)*1d-3
    if(dnn>0.d0) then
      nn(27) = n(27) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,27) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[HF]
    pd(6,28) =  &
        +k(250)*n(idx_HEj)  &
        +k(100)*n(idx_SIj)

    !d[HE_dot]/d[HF]
    pd(7,28) =  &
        +k(250)*n(idx_HEj)

    !d[HF_dot]/d[HF]
    pd(28,28) =  &
        -k(250)*n(idx_HEj)  &
        -k(100)*n(idx_SIj)

    !d[HE+_dot]/d[HF]
    pd(48,28) =  &
        -k(250)*n(idx_HEj)

    !d[SI+_dot]/d[HF]
    pd(53,28) =  &
        -k(100)*n(idx_SIj)

    !d[SIF+_dot]/d[HF]
    pd(66,28) =  &
        +k(100)*n(idx_SIj)

    !d[F+_dot]/d[HF]
    pd(70,28) =  &
        +k(250)*n(idx_HEj)

    !d[Tgas_dot]/d[HF]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(28)*1d-3
    if(dnn>0.d0) then
      nn(28) = n(28) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,28) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[CH]
    pd(1,29) =  &
        +k(200)*n(idx_O)

    !d[H_dot]/d[CH]
    pd(6,29) =  &
        +k(152)*n(idx_C)  &
        -k(92)*n(idx_H)  &
        +k(181)*n(idx_O)  &
        +k(23)*n(idx_S)  &
        +k(165)*n(idx_N)  &
        +2.d0*k(108)*n(idx_H)  &
        +k(42)*n(idx_H2)  &
        -k(108)*n(idx_H)

    !d[C_dot]/d[CH]
    pd(8,29) =  &
        +k(92)*n(idx_H)  &
        +k(131)*n(idx_N)  &
        +k(99)*n(idx_S)  &
        +k(98)*n(idx_O)  &
        -k(152)*n(idx_C)  &
        +k(108)*n(idx_H)

    !d[CO_dot]/d[CH]
    pd(10,29) =  &
        +k(181)*n(idx_O)

    !d[N_dot]/d[CH]
    pd(11,29) =  &
        -k(131)*n(idx_N)  &
        -k(165)*n(idx_N)

    !d[S_dot]/d[CH]
    pd(13,29) =  &
        -k(23)*n(idx_S)  &
        -k(99)*n(idx_S)

    !d[O_dot]/d[CH]
    pd(15,29) =  &
        -k(181)*n(idx_O)  &
        -k(200)*n(idx_O)  &
        -k(98)*n(idx_O)

    !d[H2_dot]/d[CH]
    pd(16,29) =  &
        -k(42)*n(idx_H2)  &
        +k(92)*n(idx_H)

    !d[OH_dot]/d[CH]
    pd(18,29) =  &
        +k(98)*n(idx_O)

    !d[HS_dot]/d[CH]
    pd(19,29) =  &
        +k(99)*n(idx_S)

    !d[CS_dot]/d[CH]
    pd(23,29) =  &
        +k(23)*n(idx_S)

    !d[CN_dot]/d[CH]
    pd(24,29) =  &
        +k(165)*n(idx_N)

    !d[CH_dot]/d[CH]
    pd(29,29) =  &
        -k(181)*n(idx_O)  &
        -k(152)*n(idx_C)  &
        -k(42)*n(idx_H2)  &
        -k(99)*n(idx_S)  &
        -k(92)*n(idx_H)  &
        -k(131)*n(idx_N)  &
        -k(165)*n(idx_N)  &
        -k(98)*n(idx_O)  &
        -k(23)*n(idx_S)  &
        -k(200)*n(idx_O)  &
        -k(108)*n(idx_H)

    !d[C2_dot]/d[CH]
    pd(31,29) =  &
        +k(152)*n(idx_C)

    !d[CH2_dot]/d[CH]
    pd(33,29) =  &
        +k(42)*n(idx_H2)

    !d[NH_dot]/d[CH]
    pd(34,29) =  &
        +k(131)*n(idx_N)

    !d[HCO+_dot]/d[CH]
    pd(60,29) =  &
        +k(200)*n(idx_O)

    !d[Tgas_dot]/d[CH]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(29)*1d-3
    if(dnn>0.d0) then
      nn(29) = n(29) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,29) = (dn1-dn0)/dnn
    end if

    !d[C_dot]/d[SO2]
    pd(8,30) =  &
        -k(189)*n(idx_C)

    !d[CO_dot]/d[SO2]
    pd(10,30) =  &
        +k(189)*n(idx_C)

    !d[O2_dot]/d[SO2]
    pd(12,30) =  &
        +k(31)*n(idx_O)

    !d[S_dot]/d[SO2]
    pd(13,30) =  &
        -k(125)*n(idx_S)

    !d[SO_dot]/d[SO2]
    pd(14,30) =  &
        +k(189)*n(idx_C)  &
        +k(31)*n(idx_O)  &
        +2.d0*k(125)*n(idx_S)

    !d[O_dot]/d[SO2]
    pd(15,30) =  &
        -k(31)*n(idx_O)

    !d[SO2_dot]/d[SO2]
    pd(30,30) =  &
        -k(31)*n(idx_O)  &
        -k(125)*n(idx_S)  &
        -k(189)*n(idx_C)

    !d[Tgas_dot]/d[SO2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(30)*1d-3
    if(dnn>0.d0) then
      nn(30) = n(30) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,30) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[C2]
    pd(6,31) =  &
        -k(153)*n(idx_H)

    !d[C_dot]/d[C2]
    pd(8,31) =  &
        +k(177)*n(idx_S)  &
        +k(151)*n(idx_N)  &
        +k(168)*n(idx_O)  &
        +k(153)*n(idx_H)

    !d[CO_dot]/d[C2]
    pd(10,31) =  &
        +k(168)*n(idx_O)

    !d[N_dot]/d[C2]
    pd(11,31) =  &
        -k(151)*n(idx_N)

    !d[S_dot]/d[C2]
    pd(13,31) =  &
        -k(177)*n(idx_S)

    !d[O_dot]/d[C2]
    pd(15,31) =  &
        -k(168)*n(idx_O)

    !d[CS_dot]/d[C2]
    pd(23,31) =  &
        +k(177)*n(idx_S)

    !d[CN_dot]/d[C2]
    pd(24,31) =  &
        +k(151)*n(idx_N)

    !d[CH_dot]/d[C2]
    pd(29,31) =  &
        +k(153)*n(idx_H)

    !d[C2_dot]/d[C2]
    pd(31,31) =  &
        -k(177)*n(idx_S)  &
        -k(153)*n(idx_H)  &
        -k(168)*n(idx_O)  &
        -k(151)*n(idx_N)

    !d[Tgas_dot]/d[C2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(31)*1d-3
    if(dnn>0.d0) then
      nn(31) = n(31) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,31) = (dn1-dn0)/dnn
    end if

    !d[C_dot]/d[N2]
    pd(8,32) =  &
        -k(106)*n(idx_C)

    !d[NO_dot]/d[N2]
    pd(9,32) =  &
        +k(180)*n(idx_O)

    !d[N_dot]/d[N2]
    pd(11,32) =  &
        +2.d0*k(241)  &
        +k(106)*n(idx_C)  &
        +k(180)*n(idx_O)

    !d[O_dot]/d[N2]
    pd(15,32) =  &
        -k(180)*n(idx_O)

    !d[CN_dot]/d[N2]
    pd(24,32) =  &
        +k(106)*n(idx_C)

    !d[N2_dot]/d[N2]
    pd(32,32) =  &
        -k(106)*n(idx_C)  &
        -k(180)*n(idx_O)  &
        -k(241)

    !d[Tgas_dot]/d[N2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(32)*1d-3
    if(dnn>0.d0) then
      nn(32) = n(32) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,32) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[CH2]
    pd(6,33) =  &
        -k(38)*n(idx_H)

    !d[H2_dot]/d[CH2]
    pd(16,33) =  &
        +k(38)*n(idx_H)

    !d[CH_dot]/d[CH2]
    pd(29,33) =  &
        +k(38)*n(idx_H)

    !d[CH2_dot]/d[CH2]
    pd(33,33) =  &
        -k(38)*n(idx_H)

    !d[Tgas_dot]/d[CH2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(33)*1d-3
    if(dnn>0.d0) then
      nn(33) = n(33) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,33) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[NH]
    pd(6,34) =  &
        +k(114)*n(idx_S)  &
        -k(40)*n(idx_H)  &
        +k(104)*n(idx_H2)  &
        +k(162)*n(idx_C)  &
        +k(247)*n(idx_Hj)  &
        +k(48)*n(idx_N)  &
        +k(45)*n(idx_O)

    !d[C_dot]/d[NH]
    pd(8,34) =  &
        -k(105)*n(idx_C)  &
        -k(162)*n(idx_C)

    !d[NO_dot]/d[NH]
    pd(9,34) =  &
        +k(45)*n(idx_O)

    !d[N_dot]/d[NH]
    pd(11,34) =  &
        +k(105)*n(idx_C)  &
        +k(40)*n(idx_H)  &
        +k(41)*n(idx_O)  &
        +k(158)*n(idx_S)  &
        -k(48)*n(idx_N)

    !d[S_dot]/d[NH]
    pd(13,34) =  &
        -k(114)*n(idx_S)  &
        -k(158)*n(idx_S)

    !d[O_dot]/d[NH]
    pd(15,34) =  &
        -k(45)*n(idx_O)  &
        -k(41)*n(idx_O)

    !d[H2_dot]/d[NH]
    pd(16,34) =  &
        -k(104)*n(idx_H2)  &
        +k(40)*n(idx_H)

    !d[OH_dot]/d[NH]
    pd(18,34) =  &
        +k(41)*n(idx_O)

    !d[HS_dot]/d[NH]
    pd(19,34) =  &
        +k(158)*n(idx_S)

    !d[NS_dot]/d[NH]
    pd(20,34) =  &
        +k(114)*n(idx_S)

    !d[CN_dot]/d[NH]
    pd(24,34) =  &
        +k(162)*n(idx_C)

    !d[CH_dot]/d[NH]
    pd(29,34) =  &
        +k(105)*n(idx_C)

    !d[N2_dot]/d[NH]
    pd(32,34) =  &
        +k(48)*n(idx_N)

    !d[NH_dot]/d[NH]
    pd(34,34) =  &
        -k(114)*n(idx_S)  &
        -k(40)*n(idx_H)  &
        -k(105)*n(idx_C)  &
        -k(158)*n(idx_S)  &
        -k(41)*n(idx_O)  &
        -k(162)*n(idx_C)  &
        -k(48)*n(idx_N)  &
        -k(104)*n(idx_H2)  &
        -k(45)*n(idx_O)  &
        -k(247)*n(idx_Hj)

    !d[NH2_dot]/d[NH]
    pd(39,34) =  &
        +k(104)*n(idx_H2)

    !d[H+_dot]/d[NH]
    pd(49,34) =  &
        -k(247)*n(idx_Hj)

    !d[NH+_dot]/d[NH]
    pd(50,34) =  &
        +k(247)*n(idx_Hj)

    !d[Tgas_dot]/d[NH]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(34)*1d-3
    if(dnn>0.d0) then
      nn(34) = n(34) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,34) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[HCN]
    pd(6,35) =  &
        +k(127)*n(idx_O)  &
        -k(43)*n(idx_H)

    !d[CO_dot]/d[HCN]
    pd(10,35) =  &
        +k(116)*n(idx_O)

    !d[O_dot]/d[HCN]
    pd(15,35) =  &
        -k(118)*n(idx_O)  &
        -k(116)*n(idx_O)  &
        -k(127)*n(idx_O)

    !d[H2_dot]/d[HCN]
    pd(16,35) =  &
        +k(43)*n(idx_H)

    !d[OH_dot]/d[HCN]
    pd(18,35) =  &
        +k(118)*n(idx_O)

    !d[CN_dot]/d[HCN]
    pd(24,35) =  &
        +k(43)*n(idx_H)  &
        +k(118)*n(idx_O)

    !d[NH_dot]/d[HCN]
    pd(34,35) =  &
        +k(116)*n(idx_O)

    !d[HCN_dot]/d[HCN]
    pd(35,35) =  &
        -k(118)*n(idx_O)  &
        -k(116)*n(idx_O)  &
        -k(43)*n(idx_H)  &
        -k(127)*n(idx_O)

    !d[OCN_dot]/d[HCN]
    pd(40,35) =  &
        +k(127)*n(idx_O)

    !d[Tgas_dot]/d[HCN]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(35)*1d-3
    if(dnn>0.d0) then
      nn(35) = n(35) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,35) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[CO2]
    pd(6,36) =  &
        -k(176)*n(idx_H)

    !d[NO_dot]/d[CO2]
    pd(9,36) =  &
        +k(175)*n(idx_N)

    !d[CO_dot]/d[CO2]
    pd(10,36) =  &
        +k(176)*n(idx_H)  &
        +k(175)*n(idx_N)  &
        +k(134)*n(idx_SI)

    !d[N_dot]/d[CO2]
    pd(11,36) =  &
        -k(175)*n(idx_N)

    !d[SI_dot]/d[CO2]
    pd(17,36) =  &
        -k(134)*n(idx_SI)

    !d[OH_dot]/d[CO2]
    pd(18,36) =  &
        +k(176)*n(idx_H)

    !d[CO2_dot]/d[CO2]
    pd(36,36) =  &
        -k(176)*n(idx_H)  &
        -k(134)*n(idx_SI)  &
        -k(175)*n(idx_N)

    !d[SIO_dot]/d[CO2]
    pd(37,36) =  &
        +k(134)*n(idx_SI)

    !d[Tgas_dot]/d[CO2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(36)*1d-3
    if(dnn>0.d0) then
      nn(36) = n(36) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,36) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[SIO]
    pd(6,37) =  &
        +k(140)*n(idx_Hj)  &
        +k(51)*n(idx_OH)

    !d[OH_dot]/d[SIO]
    pd(18,37) =  &
        -k(51)*n(idx_OH)

    !d[SIO_dot]/d[SIO]
    pd(37,37) =  &
        -k(140)*n(idx_Hj)  &
        -k(51)*n(idx_OH)

    !d[SIO2_dot]/d[SIO]
    pd(38,37) =  &
        +k(51)*n(idx_OH)

    !d[H+_dot]/d[SIO]
    pd(49,37) =  &
        -k(140)*n(idx_Hj)

    !d[SIO+_dot]/d[SIO]
    pd(64,37) =  &
        +k(140)*n(idx_Hj)

    !d[Tgas_dot]/d[SIO]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(37)*1d-3
    if(dnn>0.d0) then
      nn(37) = n(37) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,37) = (dn1-dn0)/dnn
    end if

    !d[HE_dot]/d[SIO2]
    pd(7,38) =  &
        +k(246)*n(idx_HEj)

    !d[O2_dot]/d[SIO2]
    pd(12,38) =  &
        +k(246)*n(idx_HEj)

    !d[SIO2_dot]/d[SIO2]
    pd(38,38) =  &
        -k(246)*n(idx_HEj)

    !d[HE+_dot]/d[SIO2]
    pd(48,38) =  &
        -k(246)*n(idx_HEj)

    !d[SI+_dot]/d[SIO2]
    pd(53,38) =  &
        +k(246)*n(idx_HEj)

    !d[Tgas_dot]/d[SIO2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(38)*1d-3
    if(dnn>0.d0) then
      nn(38) = n(38) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,38) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[NH2]
    pd(6,39) =  &
        -k(53)*n(idx_H)

    !d[H2_dot]/d[NH2]
    pd(16,39) =  &
        +k(53)*n(idx_H)

    !d[NH_dot]/d[NH2]
    pd(34,39) =  &
        +k(53)*n(idx_H)

    !d[NH2_dot]/d[NH2]
    pd(39,39) =  &
        -k(53)*n(idx_H)

    !d[Tgas_dot]/d[NH2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(39)*1d-3
    if(dnn>0.d0) then
      nn(39) = n(39) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,39) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[OCN]
    pd(6,40) =  &
        -k(56)*n(idx_H)  &
        -k(150)*n(idx_H)  &
        -k(154)*n(idx_H)

    !d[CO_dot]/d[OCN]
    pd(10,40) =  &
        +k(154)*n(idx_H)

    !d[O_dot]/d[OCN]
    pd(15,40) =  &
        +k(150)*n(idx_H)

    !d[OH_dot]/d[OCN]
    pd(18,40) =  &
        +k(56)*n(idx_H)

    !d[CN_dot]/d[OCN]
    pd(24,40) =  &
        +k(56)*n(idx_H)

    !d[NH_dot]/d[OCN]
    pd(34,40) =  &
        +k(154)*n(idx_H)

    !d[HCN_dot]/d[OCN]
    pd(35,40) =  &
        +k(150)*n(idx_H)

    !d[OCN_dot]/d[OCN]
    pd(40,40) =  &
        -k(56)*n(idx_H)  &
        -k(150)*n(idx_H)  &
        -k(154)*n(idx_H)

    !d[Tgas_dot]/d[OCN]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(40)*1d-3
    if(dnn>0.d0) then
      nn(40) = n(40) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,40) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[MG]
    pd(6,41) =  &
        +k(144)*n(idx_Hj)

    !d[C_dot]/d[MG]
    pd(8,41) =  &
        +k(111)*n(idx_Cj)

    !d[S_dot]/d[MG]
    pd(13,41) =  &
        +k(59)*n(idx_Sj)

    !d[SI_dot]/d[MG]
    pd(17,41) =  &
        +k(130)*n(idx_SIj)

    !d[SIO_dot]/d[MG]
    pd(37,41) =  &
        +k(126)*n(idx_SIOj)

    !d[MG_dot]/d[MG]
    pd(41,41) =  &
        -k(164)*n(idx_HCOj)  &
        -k(59)*n(idx_Sj)  &
        -k(111)*n(idx_Cj)  &
        -k(144)*n(idx_Hj)  &
        -k(126)*n(idx_SIOj)  &
        -k(130)*n(idx_SIj)

    !d[HCO_dot]/d[MG]
    pd(43,41) =  &
        +k(164)*n(idx_HCOj)

    !d[H+_dot]/d[MG]
    pd(49,41) =  &
        -k(144)*n(idx_Hj)

    !d[S+_dot]/d[MG]
    pd(52,41) =  &
        -k(59)*n(idx_Sj)

    !d[SI+_dot]/d[MG]
    pd(53,41) =  &
        -k(130)*n(idx_SIj)

    !d[HCO+_dot]/d[MG]
    pd(60,41) =  &
        -k(164)*n(idx_HCOj)

    !d[MG+_dot]/d[MG]
    pd(63,41) =  &
        +k(130)*n(idx_SIj)  &
        +k(164)*n(idx_HCOj)  &
        +k(111)*n(idx_Cj)  &
        +k(126)*n(idx_SIOj)  &
        +k(59)*n(idx_Sj)  &
        +k(144)*n(idx_Hj)

    !d[SIO+_dot]/d[MG]
    pd(64,41) =  &
        -k(126)*n(idx_SIOj)

    !d[C+_dot]/d[MG]
    pd(67,41) =  &
        -k(111)*n(idx_Cj)

    !d[Tgas_dot]/d[MG]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(41)*1d-3
    if(dnn>0.d0) then
      nn(41) = n(41) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,41) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[P]
    pd(6,42) =  &
        +k(84)*n(idx_Hj)

    !d[O2_dot]/d[P]
    pd(12,42) =  &
        -k(254)*n(idx_O2)

    !d[O_dot]/d[P]
    pd(15,42) =  &
        +k(254)*n(idx_O2)

    !d[P_dot]/d[P]
    pd(42,42) =  &
        -k(254)*n(idx_O2)  &
        -k(84)*n(idx_Hj)

    !d[PO_dot]/d[P]
    pd(47,42) =  &
        +k(254)*n(idx_O2)

    !d[H+_dot]/d[P]
    pd(49,42) =  &
        -k(84)*n(idx_Hj)

    !d[P+_dot]/d[P]
    pd(65,42) =  &
        +k(84)*n(idx_Hj)

    !d[Tgas_dot]/d[P]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(42)*1d-3
    if(dnn>0.d0) then
      nn(42) = n(42) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,42) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[HCO]
    pd(6,43) =  &
        -k(95)*n(idx_H)

    !d[CO_dot]/d[HCO]
    pd(10,43) =  &
        +k(95)*n(idx_H)

    !d[H2_dot]/d[HCO]
    pd(16,43) =  &
        +k(95)*n(idx_H)

    !d[HCO_dot]/d[HCO]
    pd(43,43) =  &
        -k(95)*n(idx_H)

    !d[Tgas_dot]/d[HCO]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(43)*1d-3
    if(dnn>0.d0) then
      nn(43) = n(43) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,43) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[H2O]
    pd(6,44) =  &
        -k(103)*n(idx_H)  &
        -k(146)*n(idx_H)  &
        +2.d0*k(103)*n(idx_H)

    !d[O_dot]/d[H2O]
    pd(15,44) =  &
        -k(213)*n(idx_O)

    !d[H2_dot]/d[H2O]
    pd(16,44) =  &
        +k(146)*n(idx_H)

    !d[OH_dot]/d[H2O]
    pd(18,44) =  &
        +k(146)*n(idx_H)  &
        +k(103)*n(idx_H)  &
        +2.d0*k(213)*n(idx_O)

    !d[H2O_dot]/d[H2O]
    pd(44,44) =  &
        -k(103)*n(idx_H)  &
        -k(213)*n(idx_O)  &
        -k(146)*n(idx_H)

    !d[Tgas_dot]/d[H2O]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(44)*1d-3
    if(dnn>0.d0) then
      nn(44) = n(44) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,44) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[OCS]
    pd(6,45) =  &
        -k(135)*n(idx_H)

    !d[CO_dot]/d[OCS]
    pd(10,45) =  &
        +k(135)*n(idx_H)

    !d[HS_dot]/d[OCS]
    pd(19,45) =  &
        +k(135)*n(idx_H)

    !d[OCS_dot]/d[OCS]
    pd(45,45) =  &
        -k(135)*n(idx_H)

    !d[Tgas_dot]/d[OCS]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(45)*1d-3
    if(dnn>0.d0) then
      nn(45) = n(45) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,45) = (dn1-dn0)/dnn
    end if

    !d[N_dot]/d[PN]
    pd(11,46) =  &
        -k(252)*n(idx_N)

    !d[N2_dot]/d[PN]
    pd(32,46) =  &
        +k(252)*n(idx_N)

    !d[P_dot]/d[PN]
    pd(42,46) =  &
        +k(252)*n(idx_N)

    !d[PN_dot]/d[PN]
    pd(46,46) =  &
        -k(252)*n(idx_N)

    !d[Tgas_dot]/d[PN]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(46)*1d-3
    if(dnn>0.d0) then
      nn(46) = n(46) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,46) = (dn1-dn0)/dnn
    end if

    !d[NO_dot]/d[PO]
    pd(9,47) =  &
        +k(255)*n(idx_N)

    !d[N_dot]/d[PO]
    pd(11,47) =  &
        -k(255)*n(idx_N)  &
        -k(253)*n(idx_N)

    !d[O_dot]/d[PO]
    pd(15,47) =  &
        +k(253)*n(idx_N)

    !d[P_dot]/d[PO]
    pd(42,47) =  &
        +k(255)*n(idx_N)

    !d[PN_dot]/d[PO]
    pd(46,47) =  &
        +k(253)*n(idx_N)

    !d[PO_dot]/d[PO]
    pd(47,47) =  &
        -k(255)*n(idx_N)  &
        -k(253)*n(idx_N)

    !d[Tgas_dot]/d[PO]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(47)*1d-3
    if(dnn>0.d0) then
      nn(47) = n(47) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,47) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HE+]
    pd(1,48) =  &
        -k(4)*n(idx_E)

    !d[H_dot]/d[HE+]
    pd(6,48) =  &
        -k(1)*n(idx_H)  &
        +k(250)*n(idx_HF)

    !d[HE_dot]/d[HE+]
    pd(7,48) =  &
        +k(10)*n(idx_SI)  &
        +k(246)*n(idx_SIO2)  &
        +k(4)*n(idx_E)  &
        +k(250)*n(idx_HF)  &
        +k(1)*n(idx_H)

    !d[O2_dot]/d[HE+]
    pd(12,48) =  &
        +k(246)*n(idx_SIO2)

    !d[SI_dot]/d[HE+]
    pd(17,48) =  &
        -k(10)*n(idx_SI)

    !d[HF_dot]/d[HE+]
    pd(28,48) =  &
        -k(250)*n(idx_HF)

    !d[SIO2_dot]/d[HE+]
    pd(38,48) =  &
        -k(246)*n(idx_SIO2)

    !d[HE+_dot]/d[HE+]
    pd(48,48) =  &
        -k(250)*n(idx_HF)  &
        -k(1)*n(idx_H)  &
        -k(10)*n(idx_SI)  &
        -k(246)*n(idx_SIO2)  &
        -k(4)*n(idx_E)

    !d[H+_dot]/d[HE+]
    pd(49,48) =  &
        +k(1)*n(idx_H)

    !d[SI+_dot]/d[HE+]
    pd(53,48) =  &
        +k(10)*n(idx_SI)  &
        +k(246)*n(idx_SIO2)

    !d[F+_dot]/d[HE+]
    pd(70,48) =  &
        +k(250)*n(idx_HF)

    !d[Tgas_dot]/d[HE+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(48)*1d-3
    if(dnn>0.d0) then
      nn(48) = n(48) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,48) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H+]
    pd(1,49) =  &
        -k(29)*n(idx_E)  &
        +k(226)*n(idx_Hk)

    !d[O-_dot]/d[H+]
    pd(2,49) =  &
        -k(199)*n(idx_Ok)

    !d[H-_dot]/d[H+]
    pd(3,49) =  &
        -k(65)*n(idx_Hk)  &
        -k(226)*n(idx_Hk)

    !d[C-_dot]/d[H+]
    pd(5,49) =  &
        -k(133)*n(idx_Ck)

    !d[H_dot]/d[H+]
    pd(6,49) =  &
        +k(7)*n(idx_SI)  &
        +k(133)*n(idx_Ck)  &
        +2.d0*k(65)*n(idx_Hk)  &
        +k(84)*n(idx_P)  &
        +k(167)*n(idx_S)  &
        -k(61)*n(idx_H)  &
        +k(184)*n(idx_O)  &
        +k(144)*n(idx_MG)  &
        +k(140)*n(idx_SIO)  &
        +k(222)*n(idx_H2)  &
        +k(8)*n(idx_OH)  &
        +k(29)*n(idx_E)  &
        +k(230)*n(idx_NA)  &
        +k(172)*n(idx_FE)  &
        +k(199)*n(idx_Ok)  &
        +k(247)*n(idx_NH)

    !d[C_dot]/d[H+]
    pd(8,49) =  &
        +k(133)*n(idx_Ck)

    !d[S_dot]/d[H+]
    pd(13,49) =  &
        -k(167)*n(idx_S)

    !d[O_dot]/d[H+]
    pd(15,49) =  &
        +k(199)*n(idx_Ok)  &
        -k(184)*n(idx_O)

    !d[H2_dot]/d[H+]
    pd(16,49) =  &
        -k(222)*n(idx_H2)

    !d[SI_dot]/d[H+]
    pd(17,49) =  &
        -k(7)*n(idx_SI)

    !d[OH_dot]/d[H+]
    pd(18,49) =  &
        -k(8)*n(idx_OH)

    !d[FE_dot]/d[H+]
    pd(22,49) =  &
        -k(172)*n(idx_FE)

    !d[NA_dot]/d[H+]
    pd(26,49) =  &
        -k(230)*n(idx_NA)

    !d[NH_dot]/d[H+]
    pd(34,49) =  &
        -k(247)*n(idx_NH)

    !d[SIO_dot]/d[H+]
    pd(37,49) =  &
        -k(140)*n(idx_SIO)

    !d[MG_dot]/d[H+]
    pd(41,49) =  &
        -k(144)*n(idx_MG)

    !d[P_dot]/d[H+]
    pd(42,49) =  &
        -k(84)*n(idx_P)

    !d[H+_dot]/d[H+]
    pd(49,49) =  &
        -k(230)*n(idx_NA)  &
        -k(226)*n(idx_Hk)  &
        -k(199)*n(idx_Ok)  &
        -k(167)*n(idx_S)  &
        -k(8)*n(idx_OH)  &
        -k(29)*n(idx_E)  &
        -k(61)*n(idx_H)  &
        -k(222)*n(idx_H2)  &
        -k(144)*n(idx_MG)  &
        -k(84)*n(idx_P)  &
        -k(184)*n(idx_O)  &
        -k(247)*n(idx_NH)  &
        -k(65)*n(idx_Hk)  &
        -k(133)*n(idx_Ck)  &
        -k(140)*n(idx_SIO)  &
        -k(172)*n(idx_FE)  &
        -k(7)*n(idx_SI)

    !d[NH+_dot]/d[H+]
    pd(50,49) =  &
        +k(247)*n(idx_NH)

    !d[S+_dot]/d[H+]
    pd(52,49) =  &
        +k(167)*n(idx_S)

    !d[SI+_dot]/d[H+]
    pd(53,49) =  &
        +k(7)*n(idx_SI)

    !d[OH+_dot]/d[H+]
    pd(54,49) =  &
        +k(8)*n(idx_OH)

    !d[H2+_dot]/d[H+]
    pd(56,49) =  &
        +k(61)*n(idx_H)  &
        +k(222)*n(idx_H2)  &
        +k(226)*n(idx_Hk)

    !d[FE+_dot]/d[H+]
    pd(57,49) =  &
        +k(172)*n(idx_FE)

    !d[NA+_dot]/d[H+]
    pd(59,49) =  &
        +k(230)*n(idx_NA)

    !d[O+_dot]/d[H+]
    pd(62,49) =  &
        +k(184)*n(idx_O)

    !d[MG+_dot]/d[H+]
    pd(63,49) =  &
        +k(144)*n(idx_MG)

    !d[SIO+_dot]/d[H+]
    pd(64,49) =  &
        +k(140)*n(idx_SIO)

    !d[P+_dot]/d[H+]
    pd(65,49) =  &
        +k(84)*n(idx_P)

    !d[Tgas_dot]/d[H+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(49)*1d-3
    if(dnn>0.d0) then
      nn(49) = n(49) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,49) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[NH+]
    pd(1,50) =  &
        -k(3)*n(idx_E)

    !d[H_dot]/d[NH+]
    pd(6,50) =  &
        +k(3)*n(idx_E)

    !d[N_dot]/d[NH+]
    pd(11,50) =  &
        +k(3)*n(idx_E)

    !d[NH+_dot]/d[NH+]
    pd(50,50) =  &
        -k(3)*n(idx_E)

    !d[Tgas_dot]/d[NH+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(50)*1d-3
    if(dnn>0.d0) then
      nn(50) = n(50) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,50) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HS+]
    pd(1,51) =  &
        -k(248)*n(idx_E)

    !d[H_dot]/d[HS+]
    pd(6,51) =  &
        -k(6)*n(idx_H)  &
        +k(248)*n(idx_E)

    !d[S_dot]/d[HS+]
    pd(13,51) =  &
        +k(248)*n(idx_E)

    !d[H2_dot]/d[HS+]
    pd(16,51) =  &
        +k(6)*n(idx_H)

    !d[HS+_dot]/d[HS+]
    pd(51,51) =  &
        -k(6)*n(idx_H)  &
        -k(248)*n(idx_E)

    !d[S+_dot]/d[HS+]
    pd(52,51) =  &
        +k(6)*n(idx_H)

    !d[Tgas_dot]/d[HS+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(51)*1d-3
    if(dnn>0.d0) then
      nn(51) = n(51) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,51) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[S+]
    pd(1,52) =  &
        -k(187)*n(idx_E)

    !d[H-_dot]/d[S+]
    pd(3,52) =  &
        -k(66)*n(idx_Hk)

    !d[H_dot]/d[S+]
    pd(6,52) =  &
        +k(66)*n(idx_Hk)  &
        +k(178)*n(idx_H2)

    !d[S_dot]/d[S+]
    pd(13,52) =  &
        +k(19)*n(idx_NA)  &
        +k(14)*n(idx_FE)  &
        +k(119)*n(idx_SI)  &
        +k(187)*n(idx_E)  &
        +k(59)*n(idx_MG)  &
        +k(66)*n(idx_Hk)

    !d[H2_dot]/d[S+]
    pd(16,52) =  &
        -k(178)*n(idx_H2)

    !d[SI_dot]/d[S+]
    pd(17,52) =  &
        -k(119)*n(idx_SI)

    !d[FE_dot]/d[S+]
    pd(22,52) =  &
        -k(14)*n(idx_FE)

    !d[NA_dot]/d[S+]
    pd(26,52) =  &
        -k(19)*n(idx_NA)

    !d[MG_dot]/d[S+]
    pd(41,52) =  &
        -k(59)*n(idx_MG)

    !d[HS+_dot]/d[S+]
    pd(51,52) =  &
        +k(178)*n(idx_H2)

    !d[S+_dot]/d[S+]
    pd(52,52) =  &
        -k(119)*n(idx_SI)  &
        -k(59)*n(idx_MG)  &
        -k(187)*n(idx_E)  &
        -k(19)*n(idx_NA)  &
        -k(14)*n(idx_FE)  &
        -k(178)*n(idx_H2)  &
        -k(66)*n(idx_Hk)

    !d[SI+_dot]/d[S+]
    pd(53,52) =  &
        +k(119)*n(idx_SI)

    !d[FE+_dot]/d[S+]
    pd(57,52) =  &
        +k(14)*n(idx_FE)

    !d[NA+_dot]/d[S+]
    pd(59,52) =  &
        +k(19)*n(idx_NA)

    !d[MG+_dot]/d[S+]
    pd(63,52) =  &
        +k(59)*n(idx_MG)

    !d[Tgas_dot]/d[S+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(52)*1d-3
    if(dnn>0.d0) then
      nn(52) = n(52) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,52) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[SI+]
    pd(1,53) =  &
        -k(197)*n(idx_E)

    !d[H-_dot]/d[SI+]
    pd(3,53) =  &
        -k(145)*n(idx_Hk)

    !d[H_dot]/d[SI+]
    pd(6,53) =  &
        +k(81)*n(idx_OH)  &
        -k(15)*n(idx_H)  &
        +k(100)*n(idx_HF)  &
        +k(145)*n(idx_Hk)

    !d[SI_dot]/d[SI+]
    pd(17,53) =  &
        +k(136)*n(idx_FE)  &
        +k(145)*n(idx_Hk)  &
        +k(132)*n(idx_NA)  &
        +k(197)*n(idx_E)  &
        +k(130)*n(idx_MG)

    !d[OH_dot]/d[SI+]
    pd(18,53) =  &
        -k(81)*n(idx_OH)

    !d[FE_dot]/d[SI+]
    pd(22,53) =  &
        -k(136)*n(idx_FE)

    !d[NA_dot]/d[SI+]
    pd(26,53) =  &
        -k(132)*n(idx_NA)

    !d[HF_dot]/d[SI+]
    pd(28,53) =  &
        -k(100)*n(idx_HF)

    !d[MG_dot]/d[SI+]
    pd(41,53) =  &
        -k(130)*n(idx_MG)

    !d[SI+_dot]/d[SI+]
    pd(53,53) =  &
        -k(100)*n(idx_HF)  &
        -k(197)*n(idx_E)  &
        -k(145)*n(idx_Hk)  &
        -k(130)*n(idx_MG)  &
        -k(136)*n(idx_FE)  &
        -k(132)*n(idx_NA)  &
        -k(15)*n(idx_H)  &
        -k(81)*n(idx_OH)

    !d[FE+_dot]/d[SI+]
    pd(57,53) =  &
        +k(136)*n(idx_FE)

    !d[SIH+_dot]/d[SI+]
    pd(58,53) =  &
        +k(15)*n(idx_H)

    !d[NA+_dot]/d[SI+]
    pd(59,53) =  &
        +k(132)*n(idx_NA)

    !d[MG+_dot]/d[SI+]
    pd(63,53) =  &
        +k(130)*n(idx_MG)

    !d[SIO+_dot]/d[SI+]
    pd(64,53) =  &
        +k(81)*n(idx_OH)

    !d[SIF+_dot]/d[SI+]
    pd(66,53) =  &
        +k(100)*n(idx_HF)

    !d[Tgas_dot]/d[SI+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(53)*1d-3
    if(dnn>0.d0) then
      nn(53) = n(53) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,53) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[OH+]
    pd(1,54) =  &
        -k(87)*n(idx_E)

    !d[H_dot]/d[OH+]
    pd(6,54) =  &
        +k(87)*n(idx_E)

    !d[O_dot]/d[OH+]
    pd(15,54) =  &
        +k(87)*n(idx_E)

    !d[OH+_dot]/d[OH+]
    pd(54,54) =  &
        -k(87)*n(idx_E)

    !d[Tgas_dot]/d[OH+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(54)*1d-3
    if(dnn>0.d0) then
      nn(54) = n(54) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,54) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HEH+]
    pd(1,55) =  &
        -k(25)*n(idx_E)

    !d[H_dot]/d[HEH+]
    pd(6,55) =  &
        +k(25)*n(idx_E)  &
        -k(11)*n(idx_H)

    !d[HE_dot]/d[HEH+]
    pd(7,55) =  &
        +k(11)*n(idx_H)  &
        +k(25)*n(idx_E)

    !d[HEH+_dot]/d[HEH+]
    pd(55,55) =  &
        -k(25)*n(idx_E)  &
        -k(11)*n(idx_H)

    !d[H2+_dot]/d[HEH+]
    pd(56,55) =  &
        +k(11)*n(idx_H)

    !d[Tgas_dot]/d[HEH+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(55)*1d-3
    if(dnn>0.d0) then
      nn(55) = n(55) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,55) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H2+]
    pd(1,56) =  &
        -k(102)*n(idx_E)

    !d[H_dot]/d[H2+]
    pd(6,56) =  &
        +2.d0*k(102)*n(idx_E)  &
        +k(58)*n(idx_O)  &
        -k(196)*n(idx_H)  &
        +k(32)*n(idx_HE)  &
        +k(206)*n(idx_C)

    !d[HE_dot]/d[H2+]
    pd(7,56) =  &
        -k(32)*n(idx_HE)

    !d[C_dot]/d[H2+]
    pd(8,56) =  &
        -k(206)*n(idx_C)

    !d[O_dot]/d[H2+]
    pd(15,56) =  &
        -k(58)*n(idx_O)

    !d[H2_dot]/d[H2+]
    pd(16,56) =  &
        +k(196)*n(idx_H)

    !d[H+_dot]/d[H2+]
    pd(49,56) =  &
        +k(196)*n(idx_H)

    !d[OH+_dot]/d[H2+]
    pd(54,56) =  &
        +k(58)*n(idx_O)

    !d[HEH+_dot]/d[H2+]
    pd(55,56) =  &
        +k(32)*n(idx_HE)

    !d[H2+_dot]/d[H2+]
    pd(56,56) =  &
        -k(206)*n(idx_C)  &
        -k(32)*n(idx_HE)  &
        -k(58)*n(idx_O)  &
        -k(102)*n(idx_E)  &
        -k(196)*n(idx_H)

    !d[CH+_dot]/d[H2+]
    pd(61,56) =  &
        +k(206)*n(idx_C)

    !d[Tgas_dot]/d[H2+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(56)*1d-3
    if(dnn>0.d0) then
      nn(56) = n(56) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,56) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[FE+]
    pd(1,57) =  &
        -k(142)*n(idx_E)

    !d[O-_dot]/d[FE+]
    pd(2,57) =  &
        -k(210)*n(idx_Ok)

    !d[H-_dot]/d[FE+]
    pd(3,57) =  &
        -k(24)*n(idx_Hk)

    !d[H_dot]/d[FE+]
    pd(6,57) =  &
        +k(24)*n(idx_Hk)

    !d[O_dot]/d[FE+]
    pd(15,57) =  &
        +k(210)*n(idx_Ok)

    !d[FE_dot]/d[FE+]
    pd(22,57) =  &
        +k(142)*n(idx_E)  &
        +k(24)*n(idx_Hk)  &
        +k(210)*n(idx_Ok)  &
        +k(78)*n(idx_NA)

    !d[NA_dot]/d[FE+]
    pd(26,57) =  &
        -k(78)*n(idx_NA)

    !d[FE+_dot]/d[FE+]
    pd(57,57) =  &
        -k(210)*n(idx_Ok)  &
        -k(78)*n(idx_NA)  &
        -k(142)*n(idx_E)  &
        -k(24)*n(idx_Hk)

    !d[NA+_dot]/d[FE+]
    pd(59,57) =  &
        +k(78)*n(idx_NA)

    !d[Tgas_dot]/d[FE+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(57)*1d-3
    if(dnn>0.d0) then
      nn(57) = n(57) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,57) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[SIH+]
    pd(1,58) =  &
        -k(201)*n(idx_E)

    !d[H_dot]/d[SIH+]
    pd(6,58) =  &
        +k(201)*n(idx_E)  &
        -k(110)*n(idx_H)

    !d[H2_dot]/d[SIH+]
    pd(16,58) =  &
        +k(110)*n(idx_H)

    !d[SI_dot]/d[SIH+]
    pd(17,58) =  &
        +k(201)*n(idx_E)

    !d[SI+_dot]/d[SIH+]
    pd(53,58) =  &
        +k(110)*n(idx_H)

    !d[SIH+_dot]/d[SIH+]
    pd(58,58) =  &
        -k(201)*n(idx_E)  &
        -k(110)*n(idx_H)

    !d[Tgas_dot]/d[SIH+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(58)*1d-3
    if(dnn>0.d0) then
      nn(58) = n(58) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,58) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[NA+]
    pd(1,59) =  &
        -k(190)*n(idx_E)

    !d[H-_dot]/d[NA+]
    pd(3,59) =  &
        -k(22)*n(idx_Hk)

    !d[H_dot]/d[NA+]
    pd(6,59) =  &
        +k(22)*n(idx_Hk)

    !d[NA_dot]/d[NA+]
    pd(26,59) =  &
        +k(190)*n(idx_E)  &
        +k(22)*n(idx_Hk)

    !d[NA+_dot]/d[NA+]
    pd(59,59) =  &
        -k(22)*n(idx_Hk)  &
        -k(190)*n(idx_E)

    !d[Tgas_dot]/d[NA+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(59)*1d-3
    if(dnn>0.d0) then
      nn(59) = n(59) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,59) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HCO+]
    pd(1,60) =  &
        -k(249)*n(idx_E)

    !d[H_dot]/d[HCO+]
    pd(6,60) =  &
        +k(249)*n(idx_E)

    !d[C_dot]/d[HCO+]
    pd(8,60) =  &
        -k(52)*n(idx_C)

    !d[CO_dot]/d[HCO+]
    pd(10,60) =  &
        +k(71)*n(idx_SI)  &
        +k(52)*n(idx_C)  &
        +k(249)*n(idx_E)

    !d[SI_dot]/d[HCO+]
    pd(17,60) =  &
        -k(71)*n(idx_SI)

    !d[FE_dot]/d[HCO+]
    pd(22,60) =  &
        -k(138)*n(idx_FE)

    !d[MG_dot]/d[HCO+]
    pd(41,60) =  &
        -k(164)*n(idx_MG)

    !d[HCO_dot]/d[HCO+]
    pd(43,60) =  &
        +k(164)*n(idx_MG)  &
        +k(138)*n(idx_FE)

    !d[FE+_dot]/d[HCO+]
    pd(57,60) =  &
        +k(138)*n(idx_FE)

    !d[SIH+_dot]/d[HCO+]
    pd(58,60) =  &
        +k(71)*n(idx_SI)

    !d[HCO+_dot]/d[HCO+]
    pd(60,60) =  &
        -k(71)*n(idx_SI)  &
        -k(52)*n(idx_C)  &
        -k(164)*n(idx_MG)  &
        -k(138)*n(idx_FE)  &
        -k(249)*n(idx_E)

    !d[CH+_dot]/d[HCO+]
    pd(61,60) =  &
        +k(52)*n(idx_C)

    !d[MG+_dot]/d[HCO+]
    pd(63,60) =  &
        +k(164)*n(idx_MG)

    !d[Tgas_dot]/d[HCO+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(60)*1d-3
    if(dnn>0.d0) then
      nn(60) = n(60) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,60) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[CH+]
    pd(1,61) =  &
        -k(209)*n(idx_E)

    !d[H_dot]/d[CH+]
    pd(6,61) =  &
        -k(194)*n(idx_H)  &
        +k(209)*n(idx_E)

    !d[C_dot]/d[CH+]
    pd(8,61) =  &
        +k(209)*n(idx_E)

    !d[H2_dot]/d[CH+]
    pd(16,61) =  &
        +k(194)*n(idx_H)

    !d[CH+_dot]/d[CH+]
    pd(61,61) =  &
        -k(194)*n(idx_H)  &
        -k(209)*n(idx_E)

    !d[C+_dot]/d[CH+]
    pd(67,61) =  &
        +k(194)*n(idx_H)

    !d[Tgas_dot]/d[CH+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(61)*1d-3
    if(dnn>0.d0) then
      nn(61) = n(61) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,61) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[O+]
    pd(1,62) =  &
        -k(74)*n(idx_E)

    !d[H-_dot]/d[O+]
    pd(3,62) =  &
        -k(73)*n(idx_Hk)

    !d[H_dot]/d[O+]
    pd(6,62) =  &
        +k(73)*n(idx_Hk)  &
        -k(182)*n(idx_H)  &
        +k(55)*n(idx_H2)

    !d[O_dot]/d[O+]
    pd(15,62) =  &
        +k(68)*n(idx_FE)  &
        +k(74)*n(idx_E)  &
        +k(182)*n(idx_H)  &
        +k(73)*n(idx_Hk)

    !d[H2_dot]/d[O+]
    pd(16,62) =  &
        -k(55)*n(idx_H2)

    !d[FE_dot]/d[O+]
    pd(22,62) =  &
        -k(68)*n(idx_FE)

    !d[H+_dot]/d[O+]
    pd(49,62) =  &
        +k(182)*n(idx_H)

    !d[OH+_dot]/d[O+]
    pd(54,62) =  &
        +k(55)*n(idx_H2)

    !d[FE+_dot]/d[O+]
    pd(57,62) =  &
        +k(68)*n(idx_FE)

    !d[O+_dot]/d[O+]
    pd(62,62) =  &
        -k(55)*n(idx_H2)  &
        -k(74)*n(idx_E)  &
        -k(73)*n(idx_Hk)  &
        -k(182)*n(idx_H)  &
        -k(68)*n(idx_FE)

    !d[Tgas_dot]/d[O+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(62)*1d-3
    if(dnn>0.d0) then
      nn(62) = n(62) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,62) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[MG+]
    pd(1,63) =  &
        -k(64)*n(idx_E)

    !d[O-_dot]/d[MG+]
    pd(2,63) =  &
        -k(63)*n(idx_Ok)

    !d[H-_dot]/d[MG+]
    pd(3,63) =  &
        -k(214)*n(idx_Hk)

    !d[H_dot]/d[MG+]
    pd(6,63) =  &
        +k(214)*n(idx_Hk)

    !d[O_dot]/d[MG+]
    pd(15,63) =  &
        +k(63)*n(idx_Ok)

    !d[NA_dot]/d[MG+]
    pd(26,63) =  &
        -k(195)*n(idx_NA)

    !d[MG_dot]/d[MG+]
    pd(41,63) =  &
        +k(195)*n(idx_NA)  &
        +k(214)*n(idx_Hk)  &
        +k(64)*n(idx_E)  &
        +k(63)*n(idx_Ok)

    !d[NA+_dot]/d[MG+]
    pd(59,63) =  &
        +k(195)*n(idx_NA)

    !d[MG+_dot]/d[MG+]
    pd(63,63) =  &
        -k(64)*n(idx_E)  &
        -k(214)*n(idx_Hk)  &
        -k(63)*n(idx_Ok)  &
        -k(195)*n(idx_NA)

    !d[Tgas_dot]/d[MG+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(63)*1d-3
    if(dnn>0.d0) then
      nn(63) = n(63) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,63) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[SIO+]
    pd(1,64) =  &
        -k(115)*n(idx_E)

    !d[C_dot]/d[SIO+]
    pd(8,64) =  &
        -k(202)*n(idx_C)

    !d[NO_dot]/d[SIO+]
    pd(9,64) =  &
        +k(156)*n(idx_N)

    !d[CO_dot]/d[SIO+]
    pd(10,64) =  &
        +k(202)*n(idx_C)

    !d[N_dot]/d[SIO+]
    pd(11,64) =  &
        -k(156)*n(idx_N)

    !d[O2_dot]/d[SIO+]
    pd(12,64) =  &
        +k(97)*n(idx_O)

    !d[O_dot]/d[SIO+]
    pd(15,64) =  &
        -k(97)*n(idx_O)  &
        +k(115)*n(idx_E)

    !d[SI_dot]/d[SIO+]
    pd(17,64) =  &
        +k(115)*n(idx_E)

    !d[FE_dot]/d[SIO+]
    pd(22,64) =  &
        -k(79)*n(idx_FE)

    !d[SIO_dot]/d[SIO+]
    pd(37,64) =  &
        +k(126)*n(idx_MG)  &
        +k(79)*n(idx_FE)

    !d[MG_dot]/d[SIO+]
    pd(41,64) =  &
        -k(126)*n(idx_MG)

    !d[SI+_dot]/d[SIO+]
    pd(53,64) =  &
        +k(97)*n(idx_O)  &
        +k(156)*n(idx_N)  &
        +k(202)*n(idx_C)

    !d[FE+_dot]/d[SIO+]
    pd(57,64) =  &
        +k(79)*n(idx_FE)

    !d[MG+_dot]/d[SIO+]
    pd(63,64) =  &
        +k(126)*n(idx_MG)

    !d[SIO+_dot]/d[SIO+]
    pd(64,64) =  &
        -k(97)*n(idx_O)  &
        -k(156)*n(idx_N)  &
        -k(126)*n(idx_MG)  &
        -k(202)*n(idx_C)  &
        -k(79)*n(idx_FE)  &
        -k(115)*n(idx_E)

    !d[Tgas_dot]/d[SIO+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(64)*1d-3
    if(dnn>0.d0) then
      nn(64) = n(64) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,64) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[P+]
    pd(1,65) =  &
        -k(139)*n(idx_E)

    !d[SI_dot]/d[P+]
    pd(17,65) =  &
        -k(143)*n(idx_SI)

    !d[P_dot]/d[P+]
    pd(42,65) =  &
        +k(139)*n(idx_E)  &
        +k(143)*n(idx_SI)

    !d[SI+_dot]/d[P+]
    pd(53,65) =  &
        +k(143)*n(idx_SI)

    !d[P+_dot]/d[P+]
    pd(65,65) =  &
        -k(139)*n(idx_E)  &
        -k(143)*n(idx_SI)

    !d[Tgas_dot]/d[P+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(65)*1d-3
    if(dnn>0.d0) then
      nn(65) = n(65) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,65) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[SIF+]
    pd(1,66) =  &
        -k(109)*n(idx_E)

    !d[SI_dot]/d[SIF+]
    pd(17,66) =  &
        +k(109)*n(idx_E)

    !d[F_dot]/d[SIF+]
    pd(27,66) =  &
        +k(109)*n(idx_E)

    !d[SIF+_dot]/d[SIF+]
    pd(66,66) =  &
        -k(109)*n(idx_E)

    !d[Tgas_dot]/d[SIF+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(66)*1d-3
    if(dnn>0.d0) then
      nn(66) = n(66) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,66) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[C+]
    pd(1,67) =  &
        -k(243)*n(idx_E)

    !d[C_dot]/d[C+]
    pd(8,67) =  &
        +k(121)*n(idx_SI)  &
        +k(111)*n(idx_MG)  &
        +k(124)*n(idx_FE)  &
        +k(243)*n(idx_E)

    !d[SI_dot]/d[C+]
    pd(17,67) =  &
        -k(121)*n(idx_SI)

    !d[FE_dot]/d[C+]
    pd(22,67) =  &
        -k(124)*n(idx_FE)

    !d[MG_dot]/d[C+]
    pd(41,67) =  &
        -k(111)*n(idx_MG)

    !d[SI+_dot]/d[C+]
    pd(53,67) =  &
        +k(121)*n(idx_SI)

    !d[FE+_dot]/d[C+]
    pd(57,67) =  &
        +k(124)*n(idx_FE)

    !d[MG+_dot]/d[C+]
    pd(63,67) =  &
        +k(111)*n(idx_MG)

    !d[C+_dot]/d[C+]
    pd(67,67) =  &
        -k(111)*n(idx_MG)  &
        -k(243)*n(idx_E)  &
        -k(121)*n(idx_SI)  &
        -k(124)*n(idx_FE)

    !d[Tgas_dot]/d[C+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(67)*1d-3
    if(dnn>0.d0) then
      nn(67) = n(67) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,67) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[N+]
    pd(1,68) =  &
        -k(244)*n(idx_E)

    !d[N_dot]/d[N+]
    pd(11,68) =  &
        +k(244)*n(idx_E)

    !d[N+_dot]/d[N+]
    pd(68,68) =  &
        -k(244)*n(idx_E)

    !d[Tgas_dot]/d[N+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(68)*1d-3
    if(dnn>0.d0) then
      nn(68) = n(68) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,68) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[CO+]
    pd(1,69) =  &
        -k(245)*n(idx_E)

    !d[C_dot]/d[CO+]
    pd(8,69) =  &
        +k(245)*n(idx_E)

    !d[O_dot]/d[CO+]
    pd(15,69) =  &
        +k(245)*n(idx_E)

    !d[CO+_dot]/d[CO+]
    pd(69,69) =  &
        -k(245)*n(idx_E)

    !d[Tgas_dot]/d[CO+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(69)*1d-3
    if(dnn>0.d0) then
      nn(69) = n(69) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,69) = (dn1-dn0)/dnn
    end if

    !d[H2_dot]/d[F+]
    pd(16,70) =  &
        -k(251)*n(idx_H2)

    !d[F_dot]/d[F+]
    pd(27,70) =  &
        +k(251)*n(idx_H2)

    !d[H2+_dot]/d[F+]
    pd(56,70) =  &
        +k(251)*n(idx_H2)

    !d[F+_dot]/d[F+]
    pd(70,70) =  &
        -k(251)*n(idx_H2)

    !d[Tgas_dot]/d[F+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(70)*1d-3
    if(dnn>0.d0) then
      nn(70) = n(70) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,70) = (dn1-dn0)/dnn
    end if

    !d[Tgas_dot]/d[CR]
    pd(73,71) = 0.d0

    !d[Tgas_dot]/d[g]
    pd(73,72) = 0.d0

    !d[E_dot]/d[Tgas]
    pd(1,73) = 0.d0

    !d[O-_dot]/d[Tgas]
    pd(2,73) = 0.d0

    !d[H-_dot]/d[Tgas]
    pd(3,73) = 0.d0

    !d[S-_dot]/d[Tgas]
    pd(4,73) = 0.d0

    !d[C-_dot]/d[Tgas]
    pd(5,73) = 0.d0

    !d[H_dot]/d[Tgas]
    pd(6,73) = 0.d0

    !d[HE_dot]/d[Tgas]
    pd(7,73) = 0.d0

    !d[C_dot]/d[Tgas]
    pd(8,73) = 0.d0

    !d[NO_dot]/d[Tgas]
    pd(9,73) = 0.d0

    !d[CO_dot]/d[Tgas]
    pd(10,73) = 0.d0

    !d[N_dot]/d[Tgas]
    pd(11,73) = 0.d0

    !d[O2_dot]/d[Tgas]
    pd(12,73) = 0.d0

    !d[S_dot]/d[Tgas]
    pd(13,73) = 0.d0

    !d[SO_dot]/d[Tgas]
    pd(14,73) = 0.d0

    !d[O_dot]/d[Tgas]
    pd(15,73) = 0.d0

    !d[H2_dot]/d[Tgas]
    pd(16,73) = 0.d0

    !d[SI_dot]/d[Tgas]
    pd(17,73) = 0.d0

    !d[OH_dot]/d[Tgas]
    pd(18,73) = 0.d0

    !d[HS_dot]/d[Tgas]
    pd(19,73) = 0.d0

    !d[NS_dot]/d[Tgas]
    pd(20,73) = 0.d0

    !d[H2S_dot]/d[Tgas]
    pd(21,73) = 0.d0

    !d[FE_dot]/d[Tgas]
    pd(22,73) = 0.d0

    !d[CS_dot]/d[Tgas]
    pd(23,73) = 0.d0

    !d[CN_dot]/d[Tgas]
    pd(24,73) = 0.d0

    !d[S2_dot]/d[Tgas]
    pd(25,73) = 0.d0

    !d[NA_dot]/d[Tgas]
    pd(26,73) = 0.d0

    !d[F_dot]/d[Tgas]
    pd(27,73) = 0.d0

    !d[HF_dot]/d[Tgas]
    pd(28,73) = 0.d0

    !d[CH_dot]/d[Tgas]
    pd(29,73) = 0.d0

    !d[SO2_dot]/d[Tgas]
    pd(30,73) = 0.d0

    !d[C2_dot]/d[Tgas]
    pd(31,73) = 0.d0

    !d[N2_dot]/d[Tgas]
    pd(32,73) = 0.d0

    !d[CH2_dot]/d[Tgas]
    pd(33,73) = 0.d0

    !d[NH_dot]/d[Tgas]
    pd(34,73) = 0.d0

    !d[HCN_dot]/d[Tgas]
    pd(35,73) = 0.d0

    !d[CO2_dot]/d[Tgas]
    pd(36,73) = 0.d0

    !d[SIO_dot]/d[Tgas]
    pd(37,73) = 0.d0

    !d[SIO2_dot]/d[Tgas]
    pd(38,73) = 0.d0

    !d[NH2_dot]/d[Tgas]
    pd(39,73) = 0.d0

    !d[OCN_dot]/d[Tgas]
    pd(40,73) = 0.d0

    !d[MG_dot]/d[Tgas]
    pd(41,73) = 0.d0

    !d[P_dot]/d[Tgas]
    pd(42,73) = 0.d0

    !d[HCO_dot]/d[Tgas]
    pd(43,73) = 0.d0

    !d[H2O_dot]/d[Tgas]
    pd(44,73) = 0.d0

    !d[OCS_dot]/d[Tgas]
    pd(45,73) = 0.d0

    !d[PN_dot]/d[Tgas]
    pd(46,73) = 0.d0

    !d[PO_dot]/d[Tgas]
    pd(47,73) = 0.d0

    !d[HE+_dot]/d[Tgas]
    pd(48,73) = 0.d0

    !d[H+_dot]/d[Tgas]
    pd(49,73) = 0.d0

    !d[NH+_dot]/d[Tgas]
    pd(50,73) = 0.d0

    !d[HS+_dot]/d[Tgas]
    pd(51,73) = 0.d0

    !d[S+_dot]/d[Tgas]
    pd(52,73) = 0.d0

    !d[SI+_dot]/d[Tgas]
    pd(53,73) = 0.d0

    !d[OH+_dot]/d[Tgas]
    pd(54,73) = 0.d0

    !d[HEH+_dot]/d[Tgas]
    pd(55,73) = 0.d0

    !d[H2+_dot]/d[Tgas]
    pd(56,73) = 0.d0

    !d[FE+_dot]/d[Tgas]
    pd(57,73) = 0.d0

    !d[SIH+_dot]/d[Tgas]
    pd(58,73) = 0.d0

    !d[NA+_dot]/d[Tgas]
    pd(59,73) = 0.d0

    !d[HCO+_dot]/d[Tgas]
    pd(60,73) = 0.d0

    !d[CH+_dot]/d[Tgas]
    pd(61,73) = 0.d0

    !d[O+_dot]/d[Tgas]
    pd(62,73) = 0.d0

    !d[MG+_dot]/d[Tgas]
    pd(63,73) = 0.d0

    !d[SIO+_dot]/d[Tgas]
    pd(64,73) = 0.d0

    !d[P+_dot]/d[Tgas]
    pd(65,73) = 0.d0

    !d[SIF+_dot]/d[Tgas]
    pd(66,73) = 0.d0

    !d[C+_dot]/d[Tgas]
    pd(67,73) = 0.d0

    !d[N+_dot]/d[Tgas]
    pd(68,73) = 0.d0

    !d[CO+_dot]/d[Tgas]
    pd(69,73) = 0.d0

    !d[F+_dot]/d[Tgas]
    pd(70,73) = 0.d0

    !d[CR_dot]/d[Tgas]
    pd(71,73) = 0.d0

    !d[g_dot]/d[Tgas]
    pd(72,73) = 0.d0

    !d[Tgas_dot]/d[Tgas]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(73)*1d-3
    if(dnn>0.d0) then
      nn(73) = n(73) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,73) = (dn1-dn0)/dnn
    end if

    !d[dummy_dot]/d[Tgas]
    pd(74,73) = 0.d0

    !d[Tgas_dot]/d[dummy]
    pd(73,74) = 0.d0

  end subroutine jex

end module krome_ode

!############### MODULE ##############
module krome_user
  implicit none

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2018-01-25 18:47:28
  !  Changeset 9a6e335
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  integer,parameter::KROME_idx_E = 1	!E
  integer,parameter::KROME_idx_Ok = 2	!O-
  integer,parameter::KROME_idx_Hk = 3	!H-
  integer,parameter::KROME_idx_Sk = 4	!S-
  integer,parameter::KROME_idx_Ck = 5	!C-
  integer,parameter::KROME_idx_H = 6	!H
  integer,parameter::KROME_idx_HE = 7	!HE
  integer,parameter::KROME_idx_C = 8	!C
  integer,parameter::KROME_idx_NO = 9	!NO
  integer,parameter::KROME_idx_CO = 10	!CO
  integer,parameter::KROME_idx_N = 11	!N
  integer,parameter::KROME_idx_O2 = 12	!O2
  integer,parameter::KROME_idx_S = 13	!S
  integer,parameter::KROME_idx_SO = 14	!SO
  integer,parameter::KROME_idx_O = 15	!O
  integer,parameter::KROME_idx_H2 = 16	!H2
  integer,parameter::KROME_idx_SI = 17	!SI
  integer,parameter::KROME_idx_OH = 18	!OH
  integer,parameter::KROME_idx_HS = 19	!HS
  integer,parameter::KROME_idx_NS = 20	!NS
  integer,parameter::KROME_idx_H2S = 21	!H2S
  integer,parameter::KROME_idx_FE = 22	!FE
  integer,parameter::KROME_idx_CS = 23	!CS
  integer,parameter::KROME_idx_CN = 24	!CN
  integer,parameter::KROME_idx_S2 = 25	!S2
  integer,parameter::KROME_idx_NA = 26	!NA
  integer,parameter::KROME_idx_F = 27	!F
  integer,parameter::KROME_idx_HF = 28	!HF
  integer,parameter::KROME_idx_CH = 29	!CH
  integer,parameter::KROME_idx_SO2 = 30	!SO2
  integer,parameter::KROME_idx_C2 = 31	!C2
  integer,parameter::KROME_idx_N2 = 32	!N2
  integer,parameter::KROME_idx_CH2 = 33	!CH2
  integer,parameter::KROME_idx_NH = 34	!NH
  integer,parameter::KROME_idx_HCN = 35	!HCN
  integer,parameter::KROME_idx_CO2 = 36	!CO2
  integer,parameter::KROME_idx_SIO = 37	!SIO
  integer,parameter::KROME_idx_SIO2 = 38	!SIO2
  integer,parameter::KROME_idx_NH2 = 39	!NH2
  integer,parameter::KROME_idx_OCN = 40	!OCN
  integer,parameter::KROME_idx_MG = 41	!MG
  integer,parameter::KROME_idx_P = 42	!P
  integer,parameter::KROME_idx_HCO = 43	!HCO
  integer,parameter::KROME_idx_H2O = 44	!H2O
  integer,parameter::KROME_idx_OCS = 45	!OCS
  integer,parameter::KROME_idx_PN = 46	!PN
  integer,parameter::KROME_idx_PO = 47	!PO
  integer,parameter::KROME_idx_HEj = 48	!HE+
  integer,parameter::KROME_idx_Hj = 49	!H+
  integer,parameter::KROME_idx_NHj = 50	!NH+
  integer,parameter::KROME_idx_HSj = 51	!HS+
  integer,parameter::KROME_idx_Sj = 52	!S+
  integer,parameter::KROME_idx_SIj = 53	!SI+
  integer,parameter::KROME_idx_OHj = 54	!OH+
  integer,parameter::KROME_idx_HEHj = 55	!HEH+
  integer,parameter::KROME_idx_H2j = 56	!H2+
  integer,parameter::KROME_idx_FEj = 57	!FE+
  integer,parameter::KROME_idx_SIHj = 58	!SIH+
  integer,parameter::KROME_idx_NAj = 59	!NA+
  integer,parameter::KROME_idx_HCOj = 60	!HCO+
  integer,parameter::KROME_idx_CHj = 61	!CH+
  integer,parameter::KROME_idx_Oj = 62	!O+
  integer,parameter::KROME_idx_MGj = 63	!MG+
  integer,parameter::KROME_idx_SIOj = 64	!SIO+
  integer,parameter::KROME_idx_Pj = 65	!P+
  integer,parameter::KROME_idx_SIFj = 66	!SIF+
  integer,parameter::KROME_idx_Cj = 67	!C+
  integer,parameter::KROME_idx_Nj = 68	!N+
  integer,parameter::KROME_idx_COj = 69	!CO+
  integer,parameter::KROME_idx_Fj = 70	!F+
  integer,parameter::KROME_idx_CR = 71	!CR
  integer,parameter::KROME_idx_g = 72	!g
  integer,parameter::KROME_idx_Tgas = 73	!Tgas
  integer,parameter::KROME_idx_dummy = 74	!dummy

  integer,parameter::krome_idx_cool_h2 = 1
  integer,parameter::krome_idx_cool_h2gp = 2
  integer,parameter::krome_idx_cool_atomic = 3
  integer,parameter::krome_idx_cool_cen = 3
  integer,parameter::krome_idx_cool_hd = 4
  integer,parameter::krome_idx_cool_metal = 5
  integer,parameter::krome_idx_cool_z = 5
  integer,parameter::krome_idx_cool_dh = 6
  integer,parameter::krome_idx_cool_enthalpic = 6
  integer,parameter::krome_idx_cool_dust = 7
  integer,parameter::krome_idx_cool_compton = 8
  integer,parameter::krome_idx_cool_cie = 9
  integer,parameter::krome_idx_cool_cont = 10
  integer,parameter::krome_idx_cool_continuum = 10
  integer,parameter::krome_idx_cool_expansion = 11
  integer,parameter::krome_idx_cool_exp = 11
  integer,parameter::krome_idx_cool_ff = 12
  integer,parameter::krome_idx_cool_bss = 12
  integer,parameter::krome_idx_cool_custom = 13
  integer,parameter::krome_idx_cool_co = 14
  integer,parameter::krome_idx_cool_zcie = 15
  integer,parameter::krome_idx_cool_zcienouv = 16
  integer,parameter::krome_idx_cool_zextend = 17
  integer,parameter::krome_idx_cool_gh = 18
  integer,parameter::krome_ncools = 18

  integer,parameter::krome_idx_heat_chem = 1
  integer,parameter::krome_idx_heat_compress = 2
  integer,parameter::krome_idx_heat_compr = 2
  integer,parameter::krome_idx_heat_photo = 3
  integer,parameter::krome_idx_heat_dh = 4
  integer,parameter::krome_idx_heat_enthalpic = 4
  integer,parameter::krome_idx_heat_av = 5
  integer,parameter::krome_idx_heat_photoav = 5
  integer,parameter::krome_idx_heat_cr = 6
  integer,parameter::krome_idx_heat_dust = 7
  integer,parameter::krome_idx_heat_xray = 8
  integer,parameter::krome_idx_heat_viscous = 9
  integer,parameter::krome_idx_heat_visc = 9
  integer,parameter::krome_idx_heat_custom = 10
  integer,parameter::krome_idx_heat_zcie = 11
  integer,parameter::krome_nheats = 11

  integer,parameter::krome_nrea=255
  integer,parameter::krome_nmols=70
  integer,parameter::krome_nspec=74
  integer,parameter::krome_natoms=13
  integer,parameter::krome_ndust=0
  integer,parameter::krome_ndustTypes=0
  integer,parameter::krome_nPhotoBins=0
  integer,parameter::krome_nPhotoRates=0

  real*8,parameter::krome_boltzmann_eV = 8.617332478d-5 !eV / K
  real*8,parameter::krome_boltzmann_J = 1.380648d-23 !J / K
  real*8,parameter::krome_boltzmann_erg = 1.380648d-16 !erg / K
  real*8,parameter::krome_iboltzmann_eV = 1d0/krome_boltzmann_eV !K / eV
  real*8,parameter::krome_iboltzmann_erg = 1d0/krome_boltzmann_erg !K / erg
  real*8,parameter::krome_planck_eV = 4.135667516d-15 !eV s
  real*8,parameter::krome_planck_J = 6.62606957d-34 !J s
  real*8,parameter::krome_planck_erg = 6.62606957d-27 !erg s
  real*8,parameter::krome_iplanck_eV = 1d0/krome_planck_eV !1 / eV / s
  real*8,parameter::krome_iplanck_J = 1d0/krome_planck_J !1 / J / s
  real*8,parameter::krome_iplanck_erg = 1d0/krome_planck_erg !1 / erg / s
  real*8,parameter::krome_gravity = 6.674d-8 !cm3 / g / s2
  real*8,parameter::krome_e_mass = 9.10938188d-28 !g
  real*8,parameter::krome_p_mass = 1.67262158d-24 !g
  real*8,parameter::krome_n_mass = 1.674920d-24 !g
  real*8,parameter::krome_ip_mass = 1d0/krome_p_mass !1/g
  real*8,parameter::krome_clight = 2.99792458e10 !cm/s
  real*8,parameter::krome_pi = 3.14159265359d0 !#
  real*8,parameter::krome_eV_to_erg = 1.60217646d-12 !eV -> erg
  real*8,parameter::krome_ry_to_eV = 13.60569d0 !rydberg -> eV
  real*8,parameter::krome_ry_to_erg = 2.179872d-11 !rydberg -> erg
  real*8,parameter::krome_seconds_per_year = 365d0*24d0*3600d0 !yr -> s
  real*8,parameter::krome_km_to_cm = 1d5 !km -> cm
  real*8,parameter::krome_cm_to_Mpc = 1.d0/3.08d24 !cm -> Mpc
  real*8,parameter::krome_kvgas_erg = 8.d0*krome_boltzmann_erg/krome_pi/krome_p_mass !
  real*8,parameter::krome_pre_kvgas_sqrt = sqrt(8.d0*krome_boltzmann_erg/krome_pi) !
  real*8,parameter::krome_pre_planck = 2.d0*krome_planck_erg/krome_clight**2 !erg/cm2*s3
  real*8,parameter::krome_exp_planck = krome_planck_erg / krome_boltzmann_erg !s*K
  real*8,parameter::krome_stefboltz_erg = 5.670373d-5 !erg/s/cm2/K4
  real*8,parameter::krome_N_avogadro = 6.0221d23 !#
  real*8,parameter::krome_Rgas_J = 8.3144621d0 !J/K/mol
  real*8,parameter::krome_Rgas_kJ = 8.3144621d-3 !kJ/K/mol
  real*8,parameter::krome_hubble = 0.704d0 !dimensionless
  real*8,parameter::krome_Omega0 = 1.0d0 !dimensionless
  real*8,parameter::krome_Omegab = 0.0456d0 !dimensionless
  real*8,parameter::krome_Hubble0 = 1.d2*krome_hubble*krome_km_to_cm*krome_cm_to_Mpc !1/s

contains

  !*******************
  subroutine krome_set_user_crflux(argset)
    use krome_commons
    implicit none
    real*8 :: argset
    user_crflux = argset
  end subroutine krome_set_user_crflux

  !*******************
  function krome_get_user_crflux()
    use krome_commons
    implicit none
    real*8 :: krome_get_user_crflux
    krome_get_user_crflux = user_crflux
  end function krome_get_user_crflux

  !************************
  !returns the Tdust averaged over the number density
  ! as computed in the tables
  function krome_get_table_Tdust(x,Tgas)
    use krome_commons
    use krome_grfuncs
    implicit none
    real*8 :: Tgas
    real*8 :: x(nmols), krome_get_table_Tdust
    real*8::n(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_get_table_Tdust = get_table_Tdust(n(:))

  end function krome_get_table_Tdust

  !**********************
  !convert from MOCASSIN abundances to KROME
  ! xmoc(i,j): MOCASSIN matrix (note: cm-3, real*4)
  !  i=species, j=ionization level
  ! imap: matrix position index map, integer
  ! returns KROME abundances (cm-3, real*8)
  function krome_convert_xmoc(xmoc,imap) result(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*4,intent(in):: xmoc(:,:)
    real*8::x(nmols),n(nspec)
    integer,intent(in)::imap(:)

    x(:) = 0d0

    x(idx_H) = xmoc(imap(1), 1)
    x(idx_HE) = xmoc(imap(2), 1)
    x(idx_C) = xmoc(imap(6), 1)
    x(idx_N) = xmoc(imap(7), 1)
    x(idx_S) = xmoc(imap(16), 1)
    x(idx_O) = xmoc(imap(8), 1)
    x(idx_SI) = xmoc(imap(14), 1)
    x(idx_FE) = xmoc(imap(26), 1)
    x(idx_NA) = xmoc(imap(11), 1)
    x(idx_F) = xmoc(imap(9), 1)
    x(idx_MG) = xmoc(imap(12), 1)
    x(idx_P) = xmoc(imap(15), 1)
    x(idx_HEj) = xmoc(imap(2), 2)
    x(idx_Hj) = xmoc(imap(1), 2)
    x(idx_Sj) = xmoc(imap(16), 2)
    x(idx_SIj) = xmoc(imap(14), 2)
    x(idx_FEj) = xmoc(imap(26), 2)
    x(idx_NAj) = xmoc(imap(11), 2)
    x(idx_Oj) = xmoc(imap(8), 2)
    x(idx_MGj) = xmoc(imap(12), 2)
    x(idx_Pj) = xmoc(imap(15), 2)
    x(idx_Cj) = xmoc(imap(6), 2)
    x(idx_Nj) = xmoc(imap(7), 2)
    x(idx_Fj) = xmoc(imap(9), 2)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    x(idx_e) = get_electrons(n(:))

  end function krome_convert_xmoc

  !*************************
  !convert from KROME abundances to MOCASSIN
  ! x: KROME abuances (cm-3, real*8)
  ! imap: matrix position index map, integer
  ! xmoc(i,j): MOCASSIN matrix (note: cm-3, real*4)
  !  i=species, j=ionization level
  subroutine krome_return_xmoc(x,imap,xmoc)
    use krome_commons
    implicit none
    real*8,intent(in)::x(nmols)
    real*4,intent(out)::xmoc(:,:)
    integer,intent(in)::imap(:)

    xmoc(:,:) = 0d0

    xmoc(imap(1), 1) = x(idx_H)
    xmoc(imap(2), 1) = x(idx_HE)
    xmoc(imap(6), 1) = x(idx_C)
    xmoc(imap(7), 1) = x(idx_N)
    xmoc(imap(16), 1) = x(idx_S)
    xmoc(imap(8), 1) = x(idx_O)
    xmoc(imap(14), 1) = x(idx_SI)
    xmoc(imap(26), 1) = x(idx_FE)
    xmoc(imap(11), 1) = x(idx_NA)
    xmoc(imap(9), 1) = x(idx_F)
    xmoc(imap(12), 1) = x(idx_MG)
    xmoc(imap(15), 1) = x(idx_P)
    xmoc(imap(2), 2) = x(idx_HEj)
    xmoc(imap(1), 2) = x(idx_Hj)
    xmoc(imap(16), 2) = x(idx_Sj)
    xmoc(imap(14), 2) = x(idx_SIj)
    xmoc(imap(26), 2) = x(idx_FEj)
    xmoc(imap(11), 2) = x(idx_NAj)
    xmoc(imap(8), 2) = x(idx_Oj)
    xmoc(imap(12), 2) = x(idx_MGj)
    xmoc(imap(15), 2) = x(idx_Pj)
    xmoc(imap(6), 2) = x(idx_Cj)
    xmoc(imap(7), 2) = x(idx_Nj)
    xmoc(imap(9), 2) = x(idx_Fj)

  end subroutine krome_return_xmoc

  !**********************
  !convert number density (cm-3) into column
  ! density (cm-2) using the specific density
  ! column method (see help for option
  ! -columnDensityMethod)
  ! num is the number density, x(:) is the species
  ! array, Tgas is the gas temperature
  ! If the method is not JEANS, x(:) and Tgas
  ! are dummy variables
  function krome_num2col(num,x,Tgas)
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8 :: x(nmols),krome_num2col
    real*8 :: Tgas,num
    real*8::n(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_num2col = num2col(num,n(:))

  end function krome_num2col

  !***********************
  !print on screen the current values of all phys variables
  subroutine krome_print_phys_variables()
    use krome_commons
    implicit none

    print *, "Tcmb:", phys_Tcmb
    print *, "zredshift:", phys_zredshift
    print *, "orthoParaRatio:", phys_orthoParaRatio
    print *, "metallicity:", phys_metallicity
    print *, "Tfloor:", phys_Tfloor

  end subroutine krome_print_phys_variables

  !*******************
  subroutine krome_set_Tcmb(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_Tcmb = arg
  end subroutine krome_set_Tcmb

  !*******************
  function krome_get_Tcmb()
    use krome_commons
    implicit none
    real*8 :: krome_get_Tcmb
    krome_get_Tcmb = phys_Tcmb
  end function krome_get_Tcmb

  !*******************
  subroutine krome_set_zredshift(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_zredshift = arg
  end subroutine krome_set_zredshift

  !*******************
  function krome_get_zredshift()
    use krome_commons
    implicit none
    real*8 :: krome_get_zredshift
    krome_get_zredshift = phys_zredshift
  end function krome_get_zredshift

  !*******************
  subroutine krome_set_orthoParaRatio(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_orthoParaRatio = arg
  end subroutine krome_set_orthoParaRatio

  !*******************
  function krome_get_orthoParaRatio()
    use krome_commons
    implicit none
    real*8 :: krome_get_orthoParaRatio
    krome_get_orthoParaRatio = phys_orthoParaRatio
  end function krome_get_orthoParaRatio

  !*******************
  subroutine krome_set_metallicity(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_metallicity = arg
  end subroutine krome_set_metallicity

  !*******************
  function krome_get_metallicity()
    use krome_commons
    implicit none
    real*8 :: krome_get_metallicity
    krome_get_metallicity = phys_metallicity
  end function krome_get_metallicity

  !*******************
  subroutine krome_set_Tfloor(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_Tfloor = arg
  end subroutine krome_set_Tfloor

  !*******************
  function krome_get_Tfloor()
    use krome_commons
    implicit none
    real*8 :: krome_get_Tfloor
    krome_get_Tfloor = phys_Tfloor
  end function krome_get_Tfloor

  !*******************
  function krome_coolingCI(xin,inTgas)
    use krome_commons
    use krome_subs
    use krome_cooling
    use krome_constants
    real*8 :: xin(nmols)
    real*8 :: inTgas
    real*8 :: krome_coolingCI
    real*8::n(nspec),k(nZrate)
    n(:) = 0d0
    n(idx_Tgas) = inTgas
    n(1:nmols) = xin(:)
    k(:) = coolingZ_rate_tabs(inTgas)
    krome_coolingCI = coolingCI(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg
  end function krome_coolingCI

  !*******************
  function krome_coolingOI(xin,inTgas)
    use krome_commons
    use krome_subs
    use krome_cooling
    use krome_constants
    real*8 :: xin(nmols)
    real*8 :: inTgas
    real*8 :: krome_coolingOI
    real*8::n(nspec),k(nZrate)
    n(:) = 0d0
    n(idx_Tgas) = inTgas
    n(1:nmols) = xin(:)
    k(:) = coolingZ_rate_tabs(inTgas)
    krome_coolingOI = coolingOI(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg
  end function krome_coolingOI

  !*******************
  function krome_coolingOII(xin,inTgas)
    use krome_commons
    use krome_subs
    use krome_cooling
    use krome_constants
    real*8 :: xin(nmols)
    real*8 :: inTgas
    real*8 :: krome_coolingOII
    real*8::n(nspec),k(nZrate)
    n(:) = 0d0
    n(idx_Tgas) = inTgas
    n(1:nmols) = xin(:)
    k(:) = coolingZ_rate_tabs(inTgas)
    krome_coolingOII = coolingOII(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg
  end function krome_coolingOII

  !*******************
  function krome_coolingFeII(xin,inTgas)
    use krome_commons
    use krome_subs
    use krome_cooling
    use krome_constants
    real*8 :: xin(nmols)
    real*8 :: inTgas
    real*8 :: krome_coolingFeII
    real*8::n(nspec),k(nZrate)
    n(:) = 0d0
    n(idx_Tgas) = inTgas
    n(1:nmols) = xin(:)
    k(:) = coolingZ_rate_tabs(inTgas)
    krome_coolingFeII = coolingFeII(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg
  end function krome_coolingFeII

  !*******************
  function krome_coolingCII(xin,inTgas)
    use krome_commons
    use krome_subs
    use krome_cooling
    use krome_constants
    real*8 :: xin(nmols)
    real*8 :: inTgas
    real*8 :: krome_coolingCII
    real*8::n(nspec),k(nZrate)
    n(:) = 0d0
    n(idx_Tgas) = inTgas
    n(1:nmols) = xin(:)
    k(:) = coolingZ_rate_tabs(inTgas)
    krome_coolingCII = coolingCII(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg
  end function krome_coolingCII

  !*******************
  function krome_coolingFeI(xin,inTgas)
    use krome_commons
    use krome_subs
    use krome_cooling
    use krome_constants
    real*8 :: xin(nmols)
    real*8 :: inTgas
    real*8 :: krome_coolingFeI
    real*8::n(nspec),k(nZrate)
    n(:) = 0d0
    n(idx_Tgas) = inTgas
    n(1:nmols) = xin(:)
    k(:) = coolingZ_rate_tabs(inTgas)
    krome_coolingFeI = coolingFeI(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg
  end function krome_coolingFeI

  !*******************
  function krome_coolingSiII(xin,inTgas)
    use krome_commons
    use krome_subs
    use krome_cooling
    use krome_constants
    real*8 :: xin(nmols)
    real*8 :: inTgas
    real*8 :: krome_coolingSiII
    real*8::n(nspec),k(nZrate)
    n(:) = 0d0
    n(idx_Tgas) = inTgas
    n(1:nmols) = xin(:)
    k(:) = coolingZ_rate_tabs(inTgas)
    krome_coolingSiII = coolingSiII(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg
  end function krome_coolingSiII

  !*******************
  function krome_coolingSiI(xin,inTgas)
    use krome_commons
    use krome_subs
    use krome_cooling
    use krome_constants
    real*8 :: xin(nmols)
    real*8 :: inTgas
    real*8 :: krome_coolingSiI
    real*8::n(nspec),k(nZrate)
    n(:) = 0d0
    n(idx_Tgas) = inTgas
    n(1:nmols) = xin(:)
    k(:) = coolingZ_rate_tabs(inTgas)
    krome_coolingSiI = coolingSiI(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg
  end function krome_coolingSiI

  !***************************
  !dump the population of the Z cooling levels
  ! in the nfile file unit, using xvar as
  ! independent variable. alias of
  ! dump_cooling_pop subroutine
  subroutine krome_popcool_dump(xvar,nfile)
    use krome_cooling
    implicit none
    real*8 :: xvar
    integer :: nfile

    call dump_cooling_pop(xvar,nfile)

  end subroutine krome_popcool_dump

  !*****************************
  !dump the data for restart (UNDER DEVELOPEMENT!)
  !arguments: the species array and the gas temperature
  subroutine krome_store(x,Tgas,dt)
    use krome_commons
    implicit none
    integer::nfile,i
    real*8 :: x(nmols)
    real*8 :: Tgas,dt

    nfile = 92

    open(nfile,file="krome_dump.dat",status="replace")
    !dump temperature
    write(nfile,*) Tgas
    write(nfile,*) dt
    !dump species
    do i=1,nmols
      write(nfile,*) x(i)
    end do
    close(nfile)

  end subroutine krome_store

  !*****************************
  !restore the data from a dump (UNDER DEVELOPEMENT!)
  !arguments: the species array and the gas temperature
  subroutine krome_restore(x,Tgas,dt)
    use krome_commons
    implicit none
    integer::nfile,i
    real*8 :: x(nmols)
    real*8 :: Tgas,dt

    nfile = 92

    open(nfile,file="krome_dump.dat",status="old")
    !restore temperature
    read(nfile,*) Tgas
    read(nfile,*) dt
    !restore species
    do i=1,nmols
      read(nfile,*) x(i)
    end do
    close(nfile)

  end subroutine krome_restore

  !****************************
  !switch on the thermal calculation
  subroutine krome_thermo_on()
    use krome_commons
    krome_thermo_toggle = 1
  end subroutine krome_thermo_on

  !****************************
  !switch off the thermal calculation
  subroutine krome_thermo_off()
    use krome_commons
    krome_thermo_toggle = 0
  end subroutine krome_thermo_off

  !***************************
  !alias for coe in krome_subs
  ! returns the coefficient array of size krome_nrea
  ! for a given Tgas
  function krome_get_coef(Tgas,x)
    use krome_commons
    use krome_subs
    use krome_tabs
    real*8 :: krome_get_coef(nrea),x(nmols)
    real*8,value:: Tgas
    real*8::n(nspec)
    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_get_coef(:) = coe(n(:))

  end function krome_get_coef

  !****************************
  !get the mean molecular weight from
  ! mass fractions
  function krome_get_mu_x(xin)
    use krome_commons
    implicit none
    real*8 :: xin(nmols), krome_get_mu_x
    real*8::n(nmols)
    n(:) = krome_x2n(xin(:),1d0)
    krome_get_mu_x = krome_get_mu(n(:))
  end function krome_get_mu_x

  !****************************
  !return the adiabatic index from mass fractions
  ! and temperature in K
  function krome_get_gamma_x(xin,inTgas)
    use krome_commons
    implicit none
    real*8 :: inTgas
    real*8 :: xin(nmols), krome_get_gamma_x
    real*8::x(nmols),Tgas,rhogas

    Tgas = inTgas
    x(:) = krome_x2n(xin(:),1d0)
    krome_get_gamma_x = krome_get_gamma(x(:),Tgas)

  end function krome_get_gamma_x

  !***************************
  !normalize mass fractions and
  ! set charge to zero
  subroutine krome_consistent_x(x)
    use krome_commons
    use krome_constants
    implicit none
    real*8 :: x(nmols)
    real*8::isumx,sumx,xerr,imass(nmols),ee

    !1. charge consistency
    imass(:) = krome_get_imass()

    x(idx_e) = 0.d0

    ee = sum(krome_get_charges()*x(:)*imass(:))
    ee = max(ee*e_mass,0d0)
    x(idx_e) = ee

    !2. mass fraction consistency
    sumx = sum(x)

    !NOTE: uncomment here if you want some additional control
    !conservation error threshold: rise an error if above xerr
    !xerr = 1d-2
    !if(abs(sum-1d0)>xerr) then
    !   print *,"ERROR: some problem with conservation!"
    !   print *,"|sum(x)-1|=",abs(sum-1d0)
    !   stop
    !end if

    isumx = 1d0/sumx
    x(:) = x(:) * isumx

  end subroutine krome_consistent_x

  !*********************
  !return an array sized krome_nmols containing
  ! the mass fractions (#), computed from the number
  ! densities (1/cm3) and the total density in g/cm3
  function krome_n2x(n,rhogas)
    use krome_commons
    implicit none
    real*8 :: n(nmols),krome_n2x(nmols)
    real*8,value :: rhogas

    krome_n2x(:) = n(:) * krome_get_mass() / rhogas

  end function krome_n2x

  !********************
  !return an array sized krome_nmols containing
  ! the number densities (1/cm3), computed from the mass
  ! fractions and the total density in g/cm3
  function krome_x2n(x,rhogas)
    use krome_commons
    implicit none
    real*8 :: x(nmols),krome_x2n(nmols)
    real*8,value :: rhogas

    !compute densities from fractions
    krome_x2n(:) = rhogas * x(:) * krome_get_imass()

  end function krome_x2n

  !******************
  !returns free-fall time using the number density
  ! abundances of array x(:)
  function krome_get_free_fall_time(x)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::krome_get_free_fall_time
    real*8::x(:),n(nspec)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    krome_get_free_fall_time = get_free_fall_time(n(:))

  end function krome_get_free_fall_time

  !******************
  !returns free-fall time using the total mass density
  !  of gas, rhogas (g/cm3)
  function krome_get_free_fall_time_rho(rhogas)
    use krome_getphys
    implicit none
    real*8::krome_get_free_fall_time_rho
    real*8::rhogas

    krome_get_free_fall_time_rho = get_free_fall_time_rho(rhogas)

  end function krome_get_free_fall_time_rho

  !*******************
  !do only cooling and heating
  subroutine krome_thermo(x,Tgas,dt)
    use krome_commons
    use krome_cooling
    use krome_heating
    use krome_subs
    use krome_tabs
    use krome_constants
    use krome_gadiab
    implicit none
    real*8 :: x(nmols)
    real*8 :: Tgas,dt
    real*8::n(nspec),nH2dust,dTgas,k(nrea),krome_gamma

    nH2dust = 0d0
    n(:) = 0d0
    n(idx_Tgas) = Tgas
    n(1:nmols) = x(:)
    k(:) = coe_tab(n(:)) !compute coefficients
    krome_gamma = gamma_index(n(:))

    dTgas = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

    Tgas = Tgas + dTgas*dt !update gas

  end subroutine krome_thermo

  !*************************
  !get heating (erg/cm3/s) for a given species
  ! array x(:) and Tgas
  function krome_get_heating(x,inTgas)
    use krome_heating
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: inTgas
    real*8 :: x(nmols), krome_get_heating
    real*8::Tgas,k(nrea),nH2dust,n(nspec)
    n(1:nmols) = x(:)
    Tgas = inTgas
    n(idx_Tgas) = Tgas
    k(:) = coe(n(:))
    nH2dust = 0d0
    krome_get_heating = heating(n(:),Tgas,k(:),nH2dust)
  end function krome_get_heating

  !*****************************
  ! get an array containing individual heatings (erg/cm3/s)
  ! the array has size krome_nheats. see heatcool.gps
  ! for index list
  function krome_get_heating_array(x,inTgas)
    use krome_heating
    use krome_subs
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,k(nrea),nH2dust
    real*8 :: x(nmols),krome_get_heating_array(nheats)
    real*8,value :: inTgas

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = inTgas
    !#KROME_Tdust_copy
    k(:) = coe(n(:))
    Tgas = inTgas
    nH2dust = 0d0
    krome_get_heating_array(:) = get_heating_array(n(:),Tgas,k(:),nH2dust)

  end function krome_get_heating_array

  !*************************
  !get cooling (erg/cm3/s) for x(:) species array
  ! and Tgas
  function krome_get_cooling(x,inTgas)
    use krome_cooling
    use krome_commons
    implicit none
    real*8 :: inTgas
    real*8 :: x(nmols), krome_get_cooling
    real*8::Tgas,n(nspec)
    n(1:nmols) = x(:)
    Tgas = inTgas
    n(idx_Tgas) = Tgas
    krome_get_cooling = cooling(n,Tgas)
  end function krome_get_cooling

  !*****************************
  ! get an array containing individual coolings (erg/cm3/s)
  ! the array has size krome_ncools. see heatcool.gps
  ! for index list
  function krome_get_cooling_array(x,inTgas)
    use krome_cooling
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas
    real*8 :: x(nmols),krome_get_cooling_array(ncools)
    real*8,value :: inTgas

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = inTgas
    !#KROME_Tdust_copy
    Tgas = inTgas
    krome_get_cooling_array(:) = get_cooling_array(n(:),Tgas)

  end function krome_get_cooling_array

  !******************
  !alias of plot_cool
  subroutine krome_plot_cooling(n)
    use krome_cooling
    implicit none
    real*8 :: n(krome_nmols)

    call plot_cool(n(:))

  end subroutine krome_plot_cooling

  !****************
  !alias for dumping cooling in the unit nfile_in
  subroutine krome_dump_cooling(n,Tgas,nfile_in)
    use krome_cooling
    use krome_commons
    implicit none
    real*8 :: n(nmols)
    real*8 :: Tgas
    real*8::x(nspec)
    integer, optional :: nfile_in
    integer::nfile
    nfile = 31
    x(:) = 0.d0
    x(1:nmols) = n(:)
    if(present(nfile_in)) nfile = nfile_in
    call dump_cool(x(:),Tgas,nfile)

  end subroutine krome_dump_cooling

  !************************
  !conserve the total amount of nucleii,
  ! alias for conserveLin_x in subs
  subroutine krome_conserveLin_x(x,ref)
    use krome_commons
    use krome_subs
    implicit none
    real*8 :: x(nmols),ref(natoms)

    call conserveLin_x(x(:),ref(:))

  end subroutine krome_conserveLin_x

  !************************
  !conserve the total amount of nucleii,
  ! alias for conserveLin_x in subs
  function krome_conserveLinGetRef_x(x)
    use krome_commons
    use krome_subs
    implicit none
    real*8 :: x(nmols),krome_conserveLinGetRef_x(natoms)

    krome_conserveLinGetRef_x(:) = &
        conserveLinGetRef_x(x(:))

  end function krome_conserveLinGetRef_x

  !*************************
  !force conservation to array x(:)
  !using xi(:) as initial abundances.
  !alias for conserve in krome_subs
  function krome_conserve(x,xi)
    use krome_subs
    implicit none
    real*8 :: x(krome_nmols),xi(krome_nmols),krome_conserve(krome_nmols)
    real*8::n(krome_nspec),ni(krome_nspec)

    n(:) = 0d0
    ni(:) = 0d0
    n(1:krome_nmols) = x(1:krome_nmols)
    ni(1:krome_nmols) = xi(1:krome_nmols)
    n(:) = conserve(n(:), ni(:))
    krome_conserve(:) = n(1:krome_nmols)

  end function krome_conserve

  !***************************
  !get the adiabatic index for x(:) species abundances
  ! and Tgas.
  ! alias for gamma_index in krome_subs
  function krome_get_gamma(x,Tgas)
    use krome_subs
    use krome_commons
    use krome_gadiab
    real*8 :: Tgas
    real*8 :: x(nmols), krome_get_gamma
    real*8::n(nspec)
    n(:) = 0.d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    krome_get_gamma = gamma_index(n(:))
  end function krome_get_gamma

  !***************************
  !get an integer array containing the atomic numbers Z
  ! of the spcecies.
  ! alias for get_zatoms
  function krome_get_zatoms()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    integer :: krome_get_zatoms(nmols)
    integer::zatoms(nspec)

    zatoms(:) = get_zatoms()
    krome_get_zatoms(:) = zatoms(1:nmols)

  end function krome_get_zatoms

  !****************************
  !get the mean molecular weight from
  ! number density and mass density.
  ! alias for get_mu in krome_subs module
  function krome_get_mu(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8 :: x(nmols), krome_get_mu
    real*8::n(1:nspec)
    n(:) = 0d0
    n(1:nmols) = x(:)
    krome_get_mu = get_mu(n(:))
  end function krome_get_mu

  !***************************
  !get the names of the reactions as a
  ! character*50 array of krome_nrea
  ! elements
  !! !! cannot yet be called from C
  function krome_get_rnames()
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    character*50 :: krome_get_rnames(nrea)

    krome_get_rnames(:) = get_rnames()

  end function krome_get_rnames

  !*****************
  !get an array of double containing the masses in g
  ! of the species.
  ! alias for get_mass in krome_subs
  function krome_get_mass()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_mass(nmols)
    tmp(:) = get_mass()
    krome_get_mass = tmp(1:nmols)
  end function krome_get_mass

  !*****************
  !get an array of double containing the inverse
  ! of the mass (1/g) of the species
  !alias for get_imass in krome_subs
  function krome_get_imass()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_imass(nmols)
    tmp(:) = get_imass()
    krome_get_imass = tmp(1:nmols)
  end function krome_get_imass

  !***********************
  !get the total number of H nuclei
  function krome_get_Hnuclei(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    real*8::n(nspec)
    real*8 :: krome_get_Hnuclei, x(nmols)
    n(:) = 0d0
    n(1:nmols) = x(:)

    krome_get_Hnuclei = get_Hnuclei(n(:))

  end function krome_get_Hnuclei

  !*****************
  !get an array of size krome_nmols containing the
  ! charges of the species.
  ! alias for get_charges
  function krome_get_charges()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_charges(nmols)
    tmp(:) = get_charges()
    krome_get_charges = tmp(1:nmols)
  end function krome_get_charges

  !*****************
  !get an array of character*16 and size krome_nmols
  ! containing the names of all the species.
  ! alias for get_names
  !!  !! cannot yet be called from C
  function krome_get_names()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    character*16 :: krome_get_names(nmols)
    character*16::tmp(nspec)
    tmp(:) = get_names()
    krome_get_names = tmp(1:nmols)
  end function krome_get_names

  !********************
  !get space-separated header of chemical species
  function krome_get_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    character*236::krome_get_names_header
    character*16::tmp(nspec)
    integer::i

    tmp(:) = get_names()

    krome_get_names_header = ""
    do i=1,nmols
      krome_get_names_header = trim(krome_get_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_names_header

  !********************
  !get space-separated header of coolings
  function krome_get_cooling_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    character*130::krome_get_cooling_names_header
    character*16::tmp(ncools)
    integer::i

    tmp(:) = get_cooling_names()

    krome_get_cooling_names_header = ""
    do i=1,ncools
      if(trim(tmp(i))=="") cycle
      krome_get_cooling_names_header = trim(krome_get_cooling_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_cooling_names_header

  !********************
  !get space-separated header of heatings
  function krome_get_heating_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    character*87::krome_get_heating_names_header
    character*16::tmp(nheats)
    integer::i

    tmp(:) = get_heating_names()

    krome_get_heating_names_header = ""
    do i=1,nheats
      if(trim(tmp(i))=="") cycle
      krome_get_heating_names_header = trim(krome_get_heating_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_heating_names_header

  !*****************
  !get the index of the species with name name.
  ! alias for get_index
  !! !! cannot yet be called from C
  function krome_get_index(name)
    use krome_subs
    implicit none
    integer :: krome_get_index
    character*(*) :: name
    krome_get_index = get_index(name)
  end function krome_get_index

  !*******************
  !get the total density of the gas in g/cm3
  ! giving all the number densities n(:)
  function krome_get_rho(n)
    use krome_commons
    real*8 :: krome_get_rho, n(nmols)
    real*8::m(nmols)
    m(:) = krome_get_mass()
    krome_get_rho = sum(m(:)*n(:))
  end function krome_get_rho

  !*************************
  !scale the abundances of the metals contained in n(:)
  ! to Z according to Asplund+2009.
  ! note that this applies only to neutral atoms.
  subroutine krome_scale_Z(x,Z)
    use krome_commons
    use krome_getphys
    real*8 :: x(nmols)
    real*8 :: Z
    real*8::Htot,n(nspec)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0

    Htot = get_Hnuclei(n(:))
    x(idx_NA) = max(Htot * 1d1**(Z+(-5.76)), 1d-40)
    x(idx_C) = max(Htot * 1d1**(Z+(-3.57)), 1d-40)
    x(idx_FE) = max(Htot * 1d1**(Z+(-4.5)), 1d-40)
    x(idx_HF) = max(Htot * 1d1**(Z+(-11.15)), 1d-40)
    x(idx_MG) = max(Htot * 1d1**(Z+(-4.4)), 1d-40)
    x(idx_F) = max(Htot * 1d1**(Z+(-7.44)), 1d-40)
    x(idx_O) = max(Htot * 1d1**(Z+(-3.31)), 1d-40)
    x(idx_P) = max(Htot * 1d1**(Z+(-6.59)), 1d-40)
    x(idx_SI) = max(Htot * 1d1**(Z+(-4.49)), 1d-40)
    x(idx_N) = max(Htot * 1d1**(Z+(-4.17)), 1d-40)
    x(idx_S) = max(Htot * 1d1**(Z+(-4.88)), 1d-40)

  end subroutine krome_scale_Z

  !*************************
  !set the total metallicity
  ! in terms of Z/Z_solar
  subroutine krome_set_Z(xarg)
    use krome_commons
    real*8 :: xarg

    total_Z = xarg

  end subroutine krome_set_Z

  !*************************
  !set D is in terms of D_solar (D/D_sol).
  subroutine krome_set_dust_to_gas(xarg)
    use krome_commons
    real*8 :: xarg

    dust2gas_ratio = xarg

  end subroutine

  !*************************
  !set the clumping factor
  subroutine krome_set_clump(xarg)
    use krome_commons
    real*8 :: xarg

    clump_factor = xarg

  end subroutine krome_set_clump

  !***********************
  !get the number of electrons assuming
  ! total neutral charge (cations-anions)
  function krome_get_electrons(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    real*8 :: x(nmols), krome_get_electrons
    real*8::n(nspec)
    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    krome_get_electrons = get_electrons(n(:))
  end function krome_get_electrons

  !**********************
  !print on screen the first nbest highest reaction fluxes
  subroutine krome_print_best_flux(xin,Tgas,nbest)
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: xin(nmols)
    real*8 :: Tgas
    real*8::x(nmols),n(nspec)
    integer :: nbest
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux(n,Tgas,nbest)

  end subroutine krome_print_best_flux

  !*********************
  !print only the highest fluxes greater than a fraction frac
  ! of the maximum flux
  subroutine krome_print_best_flux_frac(xin,Tgas,frac)
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: xin(nmols)
    real*8 :: Tgas,frac
    real*8::n(nspec)
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux_frac(n,Tgas,frac)

  end subroutine krome_print_best_flux_frac

  !**********************
  !print the highest nbest fluxes for reactions involving
  !a given species using the index idx_find (e.g. krome_idx_H2)
  subroutine krome_print_best_flux_spec(xin,Tgas,nbest,idx_find)
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: xin(nmols)
    real*8 :: Tgas
    real*8::n(nspec)
    integer :: nbest,idx_find
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux_spec(n,Tgas,nbest,idx_find)
  end subroutine krome_print_best_flux_spec

  !*******************************
  !get an array of size krome_nrea with
  ! the fluxes of all the reactions in cm-3/s
  function krome_get_flux(n,Tgas)
    use krome_commons
    use krome_subs
    real*8 :: krome_get_flux(nrea),n(nmols)
    real*8,value :: Tgas
    real*8::x(nspec)
    x(:) = 0.d0
    x(1:nmols) = n(:)
    x(idx_Tgas) = Tgas
    krome_get_flux(:) = get_flux(x(:), Tgas)
  end function krome_get_flux

  !*****************************
  !store the fluxes to the file unit ifile
  ! using the chemical composition x(:), and the
  ! gas temperature Tgas. xvar is th value of an
  ! user-defined independent variable that
  ! can be employed for plots.
  ! the file columns are as follow
  ! rate number, xvar, absolute flux,
  !  flux/maxflux, flux fraction wrt total,
  !  reaction name (*50 string)
  subroutine krome_explore_flux(x,Tgas,ifile,xvar)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8 :: x(nmols)
    real*8 :: Tgas,xvar
    real*8::flux(nrea),fluxmax,sumflux,n(nspec)
    integer :: ifile
    integer::i
    character*50::rname(nrea)

    !get reaction names
    rname(:) = get_rnames()
    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    !get fluxes
    flux(:) = get_flux(n(:), Tgas)
    fluxmax = maxval(flux) !maximum flux
    sumflux = sum(flux) !sum of all the fluxes
    !loop on reactions
    do i=1,nrea
      write(ifile,'(I8,5E17.8e3,a3,a50)') i,xvar,Tgas,flux(i),&
          flux(i)/fluxmax, flux(i)/sumflux," ",rname(i)
    end do
    write(ifile,*)

  end subroutine krome_explore_flux

  !*********************
  !get nulcear qeff for the reactions
  function krome_get_qeff()
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8 :: krome_get_qeff(nrea)

    krome_get_qeff(:) = get_qeff()

  end function krome_get_qeff

  !************************
  !dump the fluxes to the file unit nfile
  subroutine krome_dump_flux(n,Tgas,nfile)
    use krome_commons
    real*8 :: n(nmols)
    real*8 :: Tgas
    real*8::flux(nrea)
    integer :: nfile
    integer::i

    flux(:) = krome_get_flux(n(:),Tgas)
    do i=1,nrea
      write(nfile,'(I8,E17.8e3)') i,flux(i)
    end do
    write(nfile,*)

  end subroutine krome_dump_flux

  !************************
  !dump all the evaluation of the coefficient rates in
  ! the file funit, in the range inTmin, inTmax, using
  ! imax points
  subroutine krome_dump_rates(inTmin,inTmax,imax,funit)
    use krome_commons
    use krome_subs
    implicit none
    integer::i,j
    integer :: funit,imax
    real*8 :: inTmin,inTmax
    real*8::Tmin,Tmax,Tgas,k(nrea),n(nspec)

    Tmin = log10(inTmin)
    Tmax = log10(inTmax)

    n(:) = 1d-40
    do i=1,imax
      Tgas = 1d1**((i-1)*(Tmax-Tmin)/(imax-1)+Tmin)
      n(idx_Tgas) = Tgas
      k(:) = coe(n(:))
      do j=1,nrea
        write(funit,'(E17.8e3,I8,E17.8e3)') Tgas,j,k(j)
      end do
      write(funit,*)
    end do

  end subroutine krome_dump_rates

  !************************
  !print species informations on screen
  subroutine krome_get_info(x, Tgas)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    integer::i,charges(nspec)
    real*8 :: x(nmols)
    real*8 :: Tgas
    real*8::masses(nspec)
    character*16::names(nspec)

    names(:) = get_names()
    charges(:) = get_charges()
    masses(:) = get_mass()

    print '(a4,a10,a11,a5,a11)',"#","Name","m (g)","Chrg","x"
    do i=1,size(x)
      print '(I4,a10,E11.3,I5,E11.3)',i," "//names(i),masses(i),charges(i),x(i)
    end do
    print '(a30,E11.3)'," sum",sum(x)

    print '(a14,E11.3)',"Tgas",Tgas
  end subroutine krome_get_info

  !*****************************
  subroutine krome_set_mpi_rank(xarg)
    use krome_commons
    implicit none
    integer :: xarg
    krome_mpi_rank=xarg
  end subroutine krome_set_mpi_rank

end module krome_user

!############### MODULE ##############
module krome_reduction
contains

  !**************************
  function fex_check(n,Tgas)
    use krome_commons
    use krome_tabs
    implicit none
    integer::i
    integer::r1,r2,r3
    real*8::fex_check,n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe_tab(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      r3 = arr_r3(i)
      arr_flux(i) = k(i)*n(r1)*n(r2)*n(r3)
      rrmax = max(rrmax, arr_flux(i))
    end do
    fex_check = rrmax

  end function fex_check

end module krome_reduction

!############### MODULE ##############
module krome_main

  integer::krome_call_to_fex
  !$omp threadprivate(krome_call_to_fex)

contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2018-01-25 18:47:28
  !  Changeset 9a6e335
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !********************************
  !KROME main (interface to the solver library)

  subroutine krome(x,Tgas,dt  )
    use krome_commons
    use krome_subs
    use krome_ode
    use krome_reduction
    use krome_dust
    use krome_getphys
    use krome_tabs
    implicit none
    real*8 :: Tgas,dt
    real*8 :: x(nmols)
    real*8 :: rhogas

    real*8::mass(nspec),n(nspec),tloc,xin
    real*8::rrmax,totmass,n_old(nspec),ni(nspec),invTdust(ndust)
    integer::icount,i,icount_max
    integer:: ierr

    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer::neq(1),itol,itask,istate,iopt,lrw,liw,mf
    integer::iwork(1091)
    real*8::atol(nspec),rtol(nspec)
    real*8::rwork(14256)
    logical::got_error,equil

    !****************************
    !init DLSODES (see DLSODES manual)
    call XSETF(0)!toggle solver verbosity
    got_error = .false.
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are scalar
    rtol(:) = 1.000000d-04 !relative tolerance
    atol(:) = 1.000000d-20 !absolute tolerance
    icount_max = 100 !maximum number of iterations

    itask = 1
    iopt = 0

    !MF=
    !  = 222 internal-generated JAC and sparsity
    !  = 121 user-provided JAC and internal generated sparsity
    !  =  22 internal-generated JAC but sparsity user-provided
    !  =  21 user-provided JAC and sparsity
    MF = 222
    !end init DLSODES
    !****************************

    ierr = 0 !error flag, zero==OK!
    n(:) = 0d0 !initialize densities

    n(1:nmols) = x(:)

    n(idx_Tgas) = Tgas !put temperature in the input array

    icount = 0 !count solver iterations
    istate = 1 !init solver state
    tloc = 0.d0 !set starting time

    !store initial values
    ni(:) = n(:)
    n_global(:) = n(:)

    n_old(:) = -1d99
    krome_call_to_fex = 0
    do
      icount = icount + 1
      !solve ODE
      CALL DLSODES(fex, NEQ(:), n(:), tloc, dt, &
          ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, &
          LIW, JES, MF)

      krome_call_to_fex = krome_call_to_fex + IWORK(12)
      !check DLSODES exit status
      if(istate==2) then
        exit !sucsessful integration
      elseif(istate==-1) then
        istate = 1 !exceeded internal max iterations
      elseif(istate==-5 .or. istate==-4) then
        istate = 3 !wrong sparsity recompute
      elseif(istate==-3) then
        n(:) = ni(:)
        istate = 1
      else
        got_error = .true.
      end if

      if(got_error.or.icount>icount_max) then
        if (krome_mpi_rank>0) then
          print *,krome_mpi_rank,"ERROR: wrong solver exit status!"
          print *,krome_mpi_rank,"istate:",istate
          print *,krome_mpi_rank,"iter count:",icount
          print *,krome_mpi_rank,"max iter count:",icount_max
          print *,krome_mpi_rank,"SEE KROME_ERROR_REPORT file"
        else
          print *,"ERROR: wrong solver exit status!"
          print *,"istate:",istate
          print *,"iter count:",icount
          print *,"max iter count:",icount_max
          print *,"SEE KROME_ERROR_REPORT file"
        end if
        call krome_dump(n(:), rwork(:), iwork(:), ni(:))
        stop
      end if

    end do

    !avoid negative species
    do i=1,nspec
      n(i) = max(n(i),0d0)
    end do

    n(:) = conserve(n(:),ni(:))

    !returns to user array
    x(:) = n(1:nmols)

    Tgas = n(idx_Tgas) !get new temperature

  end subroutine krome

  !*********************************
  !integrates to equilibrium using constant temperature
  subroutine krome_equilibrium(x,Tgas,verbosity)
    use krome_ode
    use krome_subs
    use krome_commons
    use krome_constants
    use krome_getphys
    use krome_tabs
    implicit none
    integer::mf,liw,lrw,itol,meth,iopt,itask,istate,neq(1)
    integer::i,imax
    integer,optional::verbosity
    integer::verbose
    real*8 :: Tgas
    real*8 :: x(nmols)
    real*8 :: rhogas
    real*8::tloc,n(nspec),mass(nspec),ni(nspec)
    real*8::dt,xin
    integer::iwork(1091)
    real*8::atol(nspec),rtol(nspec)
    real*8::rwork(14256)
    real*8::ertol,eatol,max_time,t_tot,ntot_tol,err_species
    logical::converged

    integer, save :: ncall=0
    integer, parameter :: ncall_print_frequency=20000
    integer :: ncallp
    integer::charges(nspec)
    real*8::masses(nspec)
    character*16::names(nspec)

    !set verbosity from argument
    verbose = 1 !default is verbose
    if(present(verbosity)) verbose = verbosity

    call XSETF(0)!toggle solver verbosity
    meth = 2
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are scalar
    rtol(:) = 1d-6 !relative tolerance
    atol(:) = 1d-20 !absolute tolerance

    ! Switches to decide when equilibrium has been reached
    ertol = 1d-5  ! relative min change in a species
    eatol = 1d-12 ! absolute min change in a species
    max_time=seconds_per_year*5d8 ! max time we will be integrating for

    !for DLSODES options see its manual
    iopt = 0
    itask = 1
    istate = 1

    mf = 222 !internally evaluated sparsity and jacobian
    tloc = 0d0 !initial time

    n(:) = 0d0 !initialize densities
    !copy into array
    n(nmols+1:) = 0d0
    n(1:nmols) = x(:)

    n(idx_Tgas) = Tgas

    !store previous values
    ni(:) = n(:)
    n_global(:) = ni(:)

    imax = 1000

    dt = seconds_per_year * 1d2
    t_tot = dt
    converged = .false.
    do while (.not. converged)
      do i=1,imax
        !solve ODE
        CALL DLSODES(fcn_tconst, NEQ(:), n(:), tloc, dt, ITOL, RTOL, ATOL,&
            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, jcn_dummy, MF)
        if(istate==2) then
          exit
        else
          istate=1
        end if
      end do
      !check errors
      if(istate.ne.2) then
        print *,"ERROR: no equilibrium found!"
        stop
      end if

      !avoid negative species
      do i=1,nspec
        n(i) = max(n(i),0d0)
      end do

      n(:) = conserve(n(:),ni(:))
      ! check if we have converged by comparing the error in any species with an relative abundance above eatol
      converged = maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) .lt. ertol &
          .or. t_tot .gt. max_time

      ! Increase integration time by a reasonable factor
      if(.not. converged) then
        dt = dt * 3.
        t_tot = t_tot + dt
        ni = n
        n_global = n
      endif
    enddo
    !returns to user array
    x(:) = n(1:nmols)

    if(t_tot > max_time .and. &
        maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) > 0.2 .and. verbose>0) then
    print *, 'krome_equilibrium: Did not converge in ', max_time / seconds_per_year, ' years.'
    print *, 'Tgas :', Tgas
    names(:) = get_names()
    charges(:) = get_charges()
    masses(:) = get_mass()

    print '(a4,a10,a11,a5,a16)',"#","Name","m (g)","Chrg","  Current / Last"
    do i=1,nmols
      print '(I4,a10,E11.3,I5,2E14.6,E11.3)',i," "//names(i),masses(i),charges(i),n(i),ni(i),abs(n(i) - ni(i)) / max(n(i),eatol*sum(n(1:nmols)))
    end do
    print '(a30,2E14.6)'," sum",sum(n(1:nmols)),sum(ni(1:nmols))
    print *, 'Fractional error :', maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols))))
    print *, 'Absolute and relative floors:', eatol, ertol
  end if

  ! Print info ever so often
  !$omp critical
  ncall=ncall+1
  ncallp = ncall
  !$omp end critical

  if(modulo(ncallp,ncall_print_frequency)==0 .and. verbose>0) then
    print *, 'Found equilibrium for ', ncallp, ' cells.'
  end if

end subroutine krome_equilibrium

!********************
!dummy jacobian
subroutine jcn_dummy()
  implicit none
end subroutine jcn_dummy

!*******************
!dn/dt where dT/dt=0
subroutine fcn_tconst(n,tt,x,f)
  use krome_commons
  use krome_ode
  implicit none
  integer::n,ierr
  real*8::x(n),f(n),tt
  call fex(n,tt,x(:),f(:))
  f(idx_Tgas) = 0d0
end subroutine fcn_tconst

!*******************************
subroutine krome_dump(n,rwork,iwork,ni)
  use krome_commons
  use krome_subs
  use krome_tabs
  use krome_reduction
  use krome_ode
  use krome_getphys
  integer::fnum,i,iwork(:),idx(nrea),j
  real*8::n(:),rwork(:),rrmax,k(nrea),kmax,rperc,kperc,dn(nspec),tt,ni(:)
  character*16::names(nspec),FMTi,FMTr
  character*50::rnames(nrea),fname,prex
  integer,save::mx_dump=1000 ! max nr of reports before terminating
  fnum = 99
  if (krome_mpi_rank>0) then
    write(fname,'(a,i5.5)') "KROME_ERROR_REPORT_",krome_mpi_rank
  else
    fname = "KROME_ERROR_REPORT"
  endif
  open(fnum,FILE=trim(fname),status="replace")
  tt = 0d0
  names(:) = get_names()
  rnames(:) = get_rnames()
  call fex(nspec,tt,n(:),dn(:))

  write(fnum,*) "KROME ERROR REPORT"
  write(fnum,*)
  !SPECIES
  write(fnum,*) "Species abundances"
  write(fnum,*) "**********************"
  write(fnum,'(a5,a20,3a12)') "#","name","qty","dn/dt","ninit"
  write(fnum,*) "**********************"
  do i=1,nspec
    write(fnum,'(I5,a20,3E12.3e3)') i,names(i),n(i),dn(i),ni(i)
  end do
  write(fnum,*) "**********************"

  !F90 FRIENDLY RESTART
  write(fnum,*)
  write(fnum,*) "**********************"
  write(fnum,*) "F90-friendly species"
  write(fnum,*) "**********************"
  do i=1,nspec
    write(prex,'(a,i3,a)') "x(",i,") = "
    write(fnum,*) trim(prex),ni(i),"!"//names(i)
  end do

  write(fnum,*) "**********************"

  !RATE COEFFIECIENTS
  k(:) = coe_tab(n(:))
  idx(:) = idx_sort(k(:))
  kmax = maxval(k)
  write(fnum,*)
  write(fnum,*) "Rate coefficients (sorted) at Tgas",n(idx_Tgas)
  write(fnum,*) "**********************"
  write(fnum,'(a5,2a12,a10)') "#","k","k %","  name"
  write(fnum,*) "**********************"
  do j=1,nrea
    i = idx(j)
    kperc = 0.d0
    if(kmax>0.d0) kperc = k(i)*1d2/kmax
    write(fnum,'(I5,2E12.3e3,a2,a50)') i,k(i),kperc,"  ", rnames(i)
  end do
  write(fnum,*) "**********************"
  write(fnum,*)

  !FLUXES
  call load_arrays
  rrmax = fex_check(n(:), n(idx_Tgas))
  idx(:) = idx_sort(arr_flux(:))
  write(fnum,*)
  write(fnum,*) "Reaction magnitude (sorted) [k*n1*n2*n3*...]"
  write(fnum,*) "**********************"
  write(fnum,'(a5,2a12,a10)') "#","flux","flux %","  name"
  write(fnum,*) "**********************"
  do j=1,nrea
    i = idx(j)
    rperc = 0.d0
    if(rrmax>0.d0) rperc = arr_flux(i)*1d2/rrmax
    write(fnum,'(I5,2E12.3e3,a2,a50)') i,arr_flux(i),rperc,"  ",rnames(i)
  end do
  write(fnum,*) "**********************"
  write(fnum,*)

  !SOLVER
  FMTr = "(a30,E16.7e3)"
  FMTi = "(a30,I10)"
  write(fnum,*) "Solver-related information:"
  write(fnum,FMTr) "step size last",rwork(11)
  write(fnum,FMTr) "step size attempt",rwork(12)
  write(fnum,FMTr) "time current",rwork(13)
  write(fnum,FMTr) "tol scale factor",rwork(14)
  write(fnum,FMTi) "numeber of steps",iwork(11)
  write(fnum,FMTi) "call to fex",iwork(12)
  write(fnum,FMTi) "call to jex",iwork(13)
  write(fnum,FMTi) "last order used",iwork(14)
  write(fnum,FMTi) "order attempt",iwork(15)
  write(fnum,FMTi) "idx largest error",iwork(16)
  write(fnum,FMTi) "RWORK size required",iwork(17)
  write(fnum,FMTi) "IWORK size required",iwork(18)
  write(fnum,FMTi) "NNZ in Jac",iwork(19)
  write(fnum,FMTi) "extra fex to compute jac",iwork(20)
  write(fnum,FMTi) "number of LU decomp",iwork(21)
  write(fnum,FMTi) "base address in RWORK",iwork(22)
  write(fnum,FMTi) "base address of IAN",iwork(23)
  write(fnum,FMTi) "base address of JAN",iwork(24)
  write(fnum,FMTi) "NNZ in lower LU",iwork(25)
  write(fnum,FMTi) "NNZ in upper LU",iwork(21)
  write(fnum,*) "See DLSODES manual for further details on Optional Outputs"
  write(fnum,*)
  write(fnum,*) "END KROME ERROR REPORT"
  write(fnum,*)
  close(fnum)

  mx_dump = mx_dump - 1
  if (mx_dump==0) stop

end subroutine krome_dump

!********************************
subroutine krome_init()
  use krome_commons
  use krome_tabs
  use krome_subs
  use krome_reduction
  use krome_dust
  use krome_cooling
  use krome_photo
  use krome_fit

  !init phys common variables
  !$omp parallel
  phys_Tcmb = 2.73d0
  phys_zredshift = 0d0
  phys_orthoParaRatio = 3d0
  phys_metallicity = 0d0
  phys_Tfloor = 2.73d0
  !$omp end parallel

  !init metallicity default
  !assuming solar
  total_Z = 1d0

  !default D/D_sol = Z/Z_sol
  !assuming linear scaling
  dust2gas_ratio = total_Z

  !default broadening turubulence velocity
  broadeningVturb2 = 0d0

  !default clumping factor for
  ! H2 formation on dust by Jura/Gnedin
  clump_factor = 1d0

  !default for thermo toggle is ON
  !$omp parallel
  krome_thermo_toggle = 1
  !$omp end parallel

  !load arrays with ractants/products indexes
  call load_arrays()

  !initialize cooling tabel for metals
  call coolingZ_init_tabs()

  !initialize CO cooling
  call init_coolingCO()

  !initialize the table for exp(-a/T) function
  call init_exp_table()

  call load_parts()

  !init photo reactants indexes

  !get machine precision
  krome_epsilon = epsilon(0d0)

  !load verbatim reactions
  call loadReactionsVerbatim()

end subroutine krome_init

!****************************
function krome_get_coe(x,Tgas)
  !krome_get_coe: public interface to obtain rate coefficients
  use krome_commons
  use krome_subs
  use krome_tabs
  implicit none
  real*8 :: krome_get_coe(nrea), x(nmols), Tgas
  real*8::n(nspec)

  n(:) = 0d0
  n(1:nmols) = x(:)
  n(idx_Tgas) = Tgas
  krome_get_coe(:) = coe_tab(n(:))

end function krome_get_coe

!****************************
function krome_get_coeT(Tgas)
  !krome_get_coeT: public interface to obtain rate coefficients
  ! with argument Tgas only
  use krome_commons
  use krome_subs
  use krome_tabs
  implicit none
  real*8 :: krome_get_coeT(nrea),Tgas
  real*8::n(nspec)
  n(idx_Tgas) = Tgas
  krome_get_coeT(:) = coe_tab(n(:))
end function krome_get_coeT

end module krome_main
