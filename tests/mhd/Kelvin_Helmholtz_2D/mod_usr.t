module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    call set_coordinate_system('Cartesian')
    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid 
    use mod_global_parameters
    
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    
    double precision:: qv,width,dv,k1,sigma,rho0,p0
    logical::          first
    data first/.true./

    select case(iprob)
     case(1)
       ! setup.pl -d=2
       qv=0.645d0
       width=0.05d0
       dv=0.01d0
       k1=6.d0*dpi
       sigma=0.20d0
       w(ixO^S,rho_)=one
       w(ixO^S,e_)=one
       w(ixO^S,mag(1))=0.129d0
       w(ixO^S,mag(2))=zero
       w(ixO^S,mom(1))=qv*tanh((x(ixO^S,2)-(xprobmax2+xprobmin2)/two)/width)
       w(ixO^S,mom(2))=dv*sin(k1*(x(ixO^S,1)-xprobmin1)/(xprobmax1-xprobmin1))*&
                    exp(-((x(ixO^S,2)-(xprobmax2+xprobmin2)/two)/sigma)**2)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
       if(first .and. mype==0)then
          write(*,*)'Doing 2D MHD, Kelvin-Helmholtz problem, uniform density'
          write(*,*)'qv, width, dv, k1, sigma:'
          write(*,*)qv,width,dv,k1,sigma
          first=.false.
       endif
     case(2)
       ! setup.pl -d=2
       qv=0.645d0
       width=0.05d0
       dv=0.01d0
       k1=6.d0*dpi
       sigma=0.20d0
       where(x(ixO^S,2)>=(xprobmax2+xprobmin2)/two)
          w(ixO^S,rho_)=0.5d0
       elsewhere
          w(ixO^S,rho_)=1.d0
       endwhere
       w(ixO^S,p_)=one
       w(ixO^S,mag(1))=0.129d0
       w(ixO^S,mag(2))=zero
       w(ixO^S,mom(1))=qv*tanh((x(ixO^S,2)-(xprobmax2+xprobmin2)/two)/width)
       w(ixO^S,mom(2))=dv*sin(k1*(x(ixO^S,1)-xprobmin1)/(xprobmax1-xprobmin1))*&
                    exp(-((x(ixO^S,2)-(xprobmax2+xprobmin2)/two)/sigma)**2)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
       if(first .and. mype==0)then
          write(*,*)'Doing 2D MHD, Kelvin-Helmholtz problem, two density layers'
          write(*,*)'qv, width, dv, k1, sigma:'
          write(*,*)qv,width,dv,k1,sigma
          first=.false.
       endif
     case(3)
       ! setup.pl -d=2
       w(ixO^S,rho_)=one
       w(ixO^S,e_)=one
       where(x(ixO^S,2)>=(xprobmax2-xprobmin2)/two)
          w(ixO^S,mag(1))=0.129d0
       elsewhere
          w(ixO^S,mag(1))=-0.129d0
       endwhere
       w(ixO^S,mag(2))=zero
       qv=0.645d0
       width=0.05d0
       dv=0.01d0
       k1=two*dpi
       sigma=0.20d0
       if(first .and. mype==0)then
          write(*,*)'Doing 2D MHD, Reversed Kelvin-Helmholtz problem'
          write(*,*)'qv, width, dv, k1, sigma:'
          write(*,*)qv,width,dv,k1,sigma
          write(*,*)'Assuming eta set, using value:', mhd_eta
          first=.false.
       endif
       w(ixO^S,mom(1))=qv*tanh((x(ixO^S,2)-(xprobmax2+xprobmin2)/two)/width)
       w(ixO^S,mom(2))=dv*sin(k1*(x(ixO^S,1)-xprobmin1)/(xprobmax1-xprobmin1))*&
                    exp(-((x(ixO^S,2)-(xprobmax2+xprobmin2)/two)/sigma)**2)
       call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
       write(unitterm,*)'Undefined Iprob in Userfile ',iprob
       Call mpistop(' --- initonegrid_usr ---')
    end  select 
    
  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: divb(ixI^S)

    ! output divB1
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+1)=divb(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='divB'

  end subroutine specialvarnames_output

end module mod_usr
