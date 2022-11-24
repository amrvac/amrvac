module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    call set_coordinate_system('Cartesian')
    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    ! initialize one grid 
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    
    double precision :: ranx(ixG^S), ranx1d(ixGmin1:ixGmax1)
    double precision :: cs,machs,macha,qv,ff,Rjet,dv,sigma,kx
    integer, allocatable :: seed(:)
    integer :: ix2, seed_size
    logical,save :: first=.true.

     ! setup.pl -d=2
     ! Baty and Keppens, A&A 447, 9, 2006, case Fig 13
     w(ixG^S,rho_)=one
     w(ixG^S,p_)=one
     cs=dsqrt(mhd_gamma)
     machs=3.0d0
     macha=7.0d0
     qv=cs*machs
     w(ixG^S,mag(1))=qv/macha
     w(ixG^S,mag(2))=zero
     ff=1.25d0
     Rjet=0.125d0
     dv=0.01d0
     sigma=0.2d0
     !! the following gives random initial condition, but can not be used for autotesting
     !!if(first)then
     !!  write(*,*)'seeding random number generator, on mype==',mype
     !!  call random_seed(SIZE=seed_size)
     !!  allocate(seed(seed_size))
     !!  call random_seed(GET=seed(1:seed_size))
     !!endif
     !!call random_number(ranx1d(ixGmin1:ixGmax1))
     !!do ix2=ixGmin2,ixGmax2
     !!   ranx(ixGmin1:ixGmax1,ix2)=ranx1d(ixGmin1:ixGmax1)-0.5d0
     !!enddo
     !! deterministic perturbatio for autotesting
     kx=2.0d0*dpi/(xprobmax1-xprobmin1)
     ranx(ixG^S)=0.5d0*(dcos(kx*x(ixG^S,1))+dsin(kx*x(ixG^S,1)+0.3d0))
     if(first .and. mype==0)then
        write(*,*)'Doing 2D ideal MHD, Double Kelvin-Helmholtz problem'
        write(*,*)'cs, B0, Ma, Ms, V, Rj, F, dv, sigma:'
        write(*,*) cs,qv/macha,macha,machs,qv,Rjet,ff,dv,sigma
     endif
     if(first)then
        first=.false.
     endif
     w(ixG^S,mom(1))=half*qv*(one-dtanh(ff*dabs(x(ixG^S,2))/Rjet-ff*Rjet/dabs(x(ixG^S,2))))
     w(ixG^S,mom(2))=dv*ranx(ixG^S)*dexp(-((dabs(x(ixG^S,2))-Rjet)/sigma)**2)

     call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
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
    character(len=*) :: varnames
    varnames='divB'

  end subroutine specialvarnames_output

end module mod_usr
