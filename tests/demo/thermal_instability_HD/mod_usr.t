!============================================================================2
! HD setup for thermal instability 
!=============================================================================
module mod_usr
  use mod_hd
  implicit none

  ! this is input
  double precision :: scale, grhomax, radius

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ scale, grhomax, radius

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read


  subroutine usr_init()

    call usr_params_read(par_files)

    unit_length=1.d9  ! in cm
    unit_temperature=1.d6 ! in K
    unit_numberdensity=1.d9 ! in cm^-3

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_source => special_source

    call set_coordinate_system("Cartesian")

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    character(len=20) :: printsettingformat

    printsettingformat='(1x,A50,ES15.7,A7)'

    hd_gamma=1.66666667d0

    if(mype==0) then
      write(*,*) "HD thermal instability setup:"
      write(*,printsettingformat) "density contrast ",scale," input"
      write(*,printsettingformat) "   within radius ",radius," input"
      if(hd_thermal_conduction) write(*,*) 'USING THERMAL CONDUCTION with conduction coefficient',tc_fl%tc_k_para
      if(hd_radiative_cooling) write(*,*) 'USING COOLING'
    endif

 end subroutine initglobaldata_usr


  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: xmid^D
    logical :: mask(ixG^S)
    logical, save :: first=.true.

    ^D&xmid^D=xprobmin^D+0.5d0*(xprobmax^D-xprobmin^D);
    mask(ix^S)=(^D&(x(ix^S,^D)-xmid^D)**2|+)<radius**2
    where(mask(ix^S))
       w(ix^S,rho_)  =1.0d0+scale
    elsewhere
       w(ix^S,rho_)  =1.0d0
    endwhere
    w(ix^S,mom(1))=0.0d0
    w(ix^S,mom(2))=0.0d0
    w(ix^S,p_)    =1.0d0/hd_gamma
  
    if(mype==0.and.first)then
       write(*,*)'Doing TI setup with gamma=',hd_gamma
       first=.false.
    endif

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_radiative_cooling, only: getvar_cooling_exact
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S),winit_nopert(ixI^S,1:nw)

    if(hd_radiative_cooling)then
       winit_nopert(ixI^S,rho_)=1.0d0
       winit_nopert(ixI^S,mom(1))=0.0d0
       winit_nopert(ixI^S,mom(2))=0.0d0
       winit_nopert(ixI^S,e_)=(1.0d0/hd_gamma)/(hd_gamma-1.0d0)
       call getvar_cooling_exact(qdt,ixI^L,ixO^L,winit_nopert,winit_nopert,x,bQgrid,rc_fl)
       w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    endif

  end subroutine special_source

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_radiative_cooling, only: getvar_cooling_exact, getvar_cooling
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

    double precision                   :: wlocal(ixI^S,nw)
    double precision                   :: tmp(ixI^S) 
    double precision                   :: rho(ixI^S),gradrho(ixI^S),drho(ixI^S)
    double precision                   :: kk,kk0,kk1
    integer :: idims

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call hd_get_pthermal(wlocal,x,ixI^L,ixO^L,tmp)
    ! output the temperature p/rho
    w(ixO^S,nw+1)=tmp(ixO^S)/wlocal(ixO^S,rho_)
    ! then compute schlieren density plot
    rho(ixI^S)=wlocal(ixI^S,rho_)
    gradrho(ixO^S)=zero
    do idims=1,ndim
      select case(typegrad)
        case("central")
          call gradient(rho,ixI^L,ixO^L,idims,drho)
        case("limited")
          call gradientS(rho,ixI^L,ixO^L,idims,drho)
      end select
      gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2.0d0
    enddo
    gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
    kk=5.0d0
    kk0=0.001d0
    kk1=1.0d0
    w(ixO^S,nw+2)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
    if(hd_radiative_cooling) then
       call getvar_cooling(ixI^L,ixO^L,wlocal,x,tmp,rc_fl)
       !if(mype==0)then
       !print *,'in output:EXPLICIT'
       !print *,maxval(tmp(ixO^S))
       !print *,minval(tmp(ixO^S))
       !endif
       ! check why the following blows up....
       !call getvar_cooling_exact(dt,ixI^L,ixO^L,wlocal,wlocal,x,tmp)
       !if(mype==0)then
       !print *,'in output:EXACT'
       !print *,maxval(tmp(ixO^S))
       !print *,minval(tmp(ixO^S))
       !endif
       w(ixO^S,nw+3)=tmp(ixO^S)
    endif
    
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    if(hd_radiative_cooling) then
       varnames='Te schlier Q'
    else
       varnames='Te schlier'
    endif

  end subroutine specialvarnames_output

end module mod_usr
