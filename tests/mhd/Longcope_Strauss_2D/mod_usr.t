module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    call set_coordinate_system('Cartesian')
    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters

    mhd_etah=0.0000d0

    select case(iprob)
     case(1)
       mhd_eta=zero
     case(2)
       mhd_eta=0.001d0
       mhd_etah=0.001d0
     case(3)
       mhd_eta=0.0001d0
     case(4)
       mhd_eta=0.0005d0
     case(5)
       mhd_eta=0.0004d0
     case(6)
       mhd_eta=0.0003d0
     case(7)
       mhd_eta=0.0002d0
     case(8)
       mhd_eta=3.5d-5
     case(9)
       mhd_eta=6.94d-5
    endselect 

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: a0val,phi0val,p0val,T0val,tpi
    double precision :: minbb,maxbb,minbeta,maxbeta,minrat,maxrat
    double precision :: tmp(ixG^S)
    logical, save :: first=.true.

    a0val=one/(dpi*dsqrt(two))
    phi0val=0.001d0
    p0val=0.1d0
    T0val=one
    tpi=two*dpi
    w(ix^S,p_)=p0val+two*(dsin(tpi*x(ix^S,1))*dsin(tpi*x(ix^S,2)))**2
    w(ix^S,rho_)=w(ix^S,p_)/T0val
    w(ix^S,mag(1))=tpi*a0val*dsin(tpi*x(ix^S,1))*dcos(tpi*x(ix^S,2))
    w(ix^S,mag(2))=-tpi*a0val*dcos(tpi*x(ix^S,1))*dsin(tpi*x(ix^S,2))
    w(ix^S,mom(1))=tpi*phi0val*dsin(tpi*x(ix^S,2))
    w(ix^S,mom(2))=tpi*phi0val*dsin(tpi*x(ix^S,1))

    if(first .and. mype==0)then
      write(*,*)'Doing 2D MHD, Longcope Strauss test'
      write(*,*)'a0val,phi0val,p0val:'
      write(*,*) a0val,phi0val,p0val
      tmp(ix^S)=0.5d0*(w(ix^S,mag(1))**2+w(ix^S,mag(2))**2)
      minbb=minval(tmp(ix^S))
      maxbb=maxval(tmp(ix^S))
      write(*,*) 'min B2=',minbb,' max B2=',maxbb
      tmp(ix^S)=two*w(ix^S,p_)/( w(ix^S,mag(1))**2+ w(ix^S,mag(2))**2)
      minbeta=minval(tmp(ix^S))
      maxbeta=maxval(tmp(ix^S))
      write(*,*) 'min beta=',minbeta,' max beta=',maxbeta
      tmp(ix^S)= w(ix^S,rho_)*( w(ix^S,mom(1))**2+ w(ix^S,mom(2))**2) &
                /( w(ix^S,mag(1))**2+ w(ix^S,mag(2))**2)
      minrat=minval(tmp(ix^S))
      maxrat=maxval(tmp(ix^S))
      write(*,*) 'min kin/mag=',minrat,' max kin/mag=',maxrat
      first=.false.
    endif

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    integer:: idirmin

    double precision :: current(ixI^S,7-2*ndir:3)

    call get_current(w,ixI^L,ixO^L,idirmin,current)

    w(ixO^S,nw+1)=current(ixO^S,3)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='jz'

  end subroutine specialvarnames_output

end module mod_usr
