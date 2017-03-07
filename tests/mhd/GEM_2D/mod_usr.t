module mod_usr
  use mod_mhd
  implicit none
  double precision :: sheetl,rhorat,T0,psi0,llx,lly

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    call set_coordinate_system("Cartesian")

    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters

    select case(iprob)
     case(1)
      mhd_eta=5.0d-3
     case(2)
      mhd_eta=1.0d-3
     case(3)
      mhd_eta=1.0d-4
     case(4)
      mhd_eta=1.0d-2
     case(0)
      mhd_eta=1.0d-3
      mhd_etah=1.0d-3
    case default
      call mpistop('iprob not given')
    end select
    sheetl=0.5d0
    rhorat=0.2d0
    T0=0.5d0
    psi0=0.1d0
    llx=xprobmax1-xprobmin1
    lly=xprobmax2-xprobmin2
  
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
  ! initialize one grid
  use mod_global_parameters

  integer, intent(in) :: ixG^L,ix^L
  double precision, intent(in) :: x(ixG^S,1:ndim)
  double precision, intent(inout) :: w(ixG^S,1:nw)

  integer,dimension(:),allocatable:: seed
  real::             ranx1d(ixGlo1:ixGhi1)
  {^IFTHREED
  real::             ranx2d(ixGlo1:ixGhi1,ixGlo3:ixGhi3)
  \}
  double precision:: ranx(ixG^S),dv,sigma
  integer :: ix2,seed_size
  logical, save :: first=.true.

  
  ! no initial velocity for cases iprob=1,2,3
  w(ixG^S,mom(1))  =zero
  w(ixG^S,mom(2))  =zero
  {^IFTHREED 
  w(ixG^S,mom(3))=zero
  \}
  
  if(iprob==4)then
    dv=0.01d0
    sigma=2.0d0*sheetl
    call random_seed(SIZE=seed_size)
    allocate(seed(seed_size))
    call random_seed(GET=seed(1:seed_size))
    {^IFTWOD
    call random_number(ranx1d(ixGmin1:ixGmax1))
    do ix2=ixGmin2,ixGmax2
      ranx(ixGmin1:ixGmax1,ix2)=ranx1d(ixGmin1:ixGmax1)-0.5d0
    enddo
    \}
    {^IFTHREED
    call random_number(ranx2d(ixGmin1:ixGmax1,ixGmin3:ixGmax3))
    do ix2=ixGmin2,ixGmax2
      ranx(ixGmin1:ixGmax1,ix2,ixGmin3:ixGmax3)=ranx2d(ixGmin1:ixGmax1,ixGmin3:ixGmax3)-0.5d0
    enddo
    \}
    w(ixG^S,mom(1))=dv*ranx(ixG^S)*dexp(-(x(ixG^S,2)/sigma)**2)
    w(ixG^S,mom(2))=dv*ranx(ixG^S)*dexp(-(x(ixG^S,2)/sigma)**2)
    {^IFTHREED
    w(ixG^S,mom(3))=dv*ranx(ixG^S)*dexp(-(x(ixG^S,2)/sigma)**2)
    \}
  endif
  
  ! set up the 1D equilibrium variation
  w(ixG^S,mag(1))  =dtanh(x(ixG^S,2)/sheetl)
  w(ixG^S,mag(2))  =zero
  {^IFTHREED 
  w(ixG^S,mag(3))=zero
  \}
  
  ! add the 2D island perturbation
  if(iprob<4) then
    w(ixG^S,mag(1))= w(ixG^S,mag(1))-psi0 &
                  *dcos(two*dpi*x(ixG^S,1)/llx)              &
                  *dsin(dpi*x(ixG^S,2)/lly)*dpi/lly
    w(ixG^S,mag(2))= w(ixG^S,mag(2))+psi0 &
                  *dsin(two*dpi*x(ixG^S,1)/llx)              &
                  *dcos(dpi*x(ixG^S,2)/lly)*two*dpi/llx
  endif
  
  w(ixG^S,rho_) =rhorat+one/(dcosh(x(ixG^S,2)/sheetl)**2) 
  w(ixG^S,p_)   =T0*(rhorat+one/(dcosh(x(ixG^S,2)/sheetl)**2))
  
  call mhd_to_conserved(ixG^L,ixG^L,w,x)

  {^IFTWOD
  if(mype==0.and.first)then
     write(*,*)'Doing 2D GEM challenge, resistive MHD'
     write(*,*)'iprob=',iprob
     write(*,*)'resistivity equal to=',mhd_eta
     write(*,*)'Hall parameter equal to=',mhd_etah
     if(iprob==4)write(*,*)'random velocity perturbation added'
     first=.false.
  endif
  \}

  {^IFTHREED
  if(mype==0.and.first)then
     write(*,*)'Doing 3D GEM challenge, resistive MHD'
     write(*,*)'iprob=',iprob
     write(*,*)'resistivity equal to=',mhd_eta
     write(*,*)'Hall parameter equal to=',mhd_etah
     if(iprob==4)write(*,*)'random velocity perturbation added'
     first=.false.
  endif
  \}
  
  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: current(ixI^S,7-2*ndir:3)
    double precision                   :: gradrho(ixI^S),rho(ixI^S),drho(ixI^S)
    double precision                   :: kk,kk0,grhomax,kk1
    integer                            :: idims
    
    integer :: idirmin,idir
    double precision :: curlv(ixI^S,7-2*ndir:3),vvec(ixI^S,1:ndir)
    double precision:: divb(ixI^S)

    call get_current(w,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+1)=current(ixO^S,3)
    rho(ixI^S)=w(ixI^S,rho_)
    gradrho(ixO^S)=zero
    do idims=1,ndim
       select case(typegrad)
       case("central")
         call gradient(rho,ixI^L,ixO^L,idims,drho)
       case("limited")
         call gradientS(rho,ixI^L,ixO^L,idims,drho)
       end select
       gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2
    enddo
    gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
    kk=5.0d0
    kk0=0.01d0
    kk1=1.0d0
    grhomax=1000.0d0
    w(ixO^S,nw+2)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
    do idir=1,ndir
      vvec(ixI^S,idir)=w(ixI^S,mom(idir))/w(ixI^S,rho_)
    enddo
    call curlvector(vvec,ixI^L,ixO^L,curlv,idirmin,7-2*ndir,ndir)
    w(ixO^S,nw+3)=curlv(ixO^S,3)
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+4)=divb(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='jz schlier curlvz divb'

  end subroutine specialvarnames_output

end module mod_usr
