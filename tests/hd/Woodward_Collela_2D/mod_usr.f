module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => wc2d_init_one_grid
    usr_special_bc      => specialbound_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output


    call hd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters

    hd_gamma=1.4d0

  end subroutine initglobaldata_usr

  ! Initialize one grid
  subroutine wc2d_init_one_grid(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
     ixmax1,ixmax2,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
       ixmax1,ixmax2
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

    double precision :: m1post, m2post, epost

    

    m1post=8.d0*8.25d0*dsin(dpi/3.d0)
    m2post=-8.0d0*8.25d0*dcos(dpi/3.d0)
    epost=1.165d2/(hd_gamma-1.d0)+(m1post**2+m2post**2)/(16.0d0)
    where(x(ixmin1:ixmax1,ixmin2:ixmax2,1)>(x(ixmin1:ixmax1,ixmin2:ixmax2,&
       2)/dtan(dpi/3.d0)+1.d0/6.d0))
       ! pre shock region
       w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=1.4d0
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1))=0.0d0
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(2))=0.0d0
       w(ixmin1:ixmax1,ixmin2:ixmax2,e_)=1.0d0/(hd_gamma-1.0d0)
    endwhere
    where(x(ixmin1:ixmax1,ixmin2:ixmax2,1)<=(x(ixmin1:ixmax1,ixmin2:ixmax2,&
       2)/dtan(dpi/3.d0)+1.d0/6.d0))
       ! post shock region
       w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=8.d0
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1))=m1post
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(2))=m2post
       w(ixmin1:ixmax1,ixmin2:ixmax2,e_)=epost
    endwhere

  end subroutine wc2d_init_one_grid

  subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iB,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

    integer :: ix2
    double precision :: pree, m1post, m2post, epost
    !----------------------------------------------------------------------------

    m1post=8.d0*8.25d0*dsin(dpi/3.d0)
    m2post=-8.0d0*8.25d0*dcos(dpi/3.d0)
    epost=1.165d2/(hd_gamma-1.d0)+(m1post**2+m2post**2)/(16.0d0)
    pree=1.0d0/(hd_gamma-1.0d0)
    select case(iB)
     ! implementation of fixed postshock state at left boundary
     case(1)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)=8.d0
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) =m1post
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2)) =m2post
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)  =epost
     ! implementation of bottom boundary: fixed before x<1/6, solid wall x>=1/6
     case(3)
      ! first pretend all variables are symmetric everywhere
      ! then overwrite with fixed postshock values before x<1/6
         do ix2=ixOmin2,ixOmax2
            w(ixOmin1:ixOmax1,ix2,rho_)=w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,&
               rho_)
            w(ixOmin1:ixOmax1,ix2,mom(1))=w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,&
               mom(1))
            w(ixOmin1:ixOmax1,ix2,mom(2))=-w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,&
               mom(2))
            w(ixOmin1:ixOmax1,ix2,e_)=w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,e_)
         enddo
         where(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)<=1.d0/6.d0)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)=8.d0
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1))=m1post
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))=m2post
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=epost
         endwhere
     ! implementation of top boundary: pre/post shock state, time-dependent
     case(4)
      ! pre shock region
      where(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1)> (1.0d1*qt/dsin(dpi/3.d0)+x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2)/dtan(dpi/3.d0)+1.d0/6.d0))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)= 1.4d0
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = 0.d0
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2)) = 0.d0
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=pree
      endwhere
      ! post shock region
      where(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1)<= (1.0d1*qt/dsin(dpi/3.d0)+x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2)/dtan(dpi/3.d0)+1.d0/6.d0))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)= 8.d0
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1))=m1post
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))=m2post
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=epost
      endwhere
     case default
      call mpistop("Special boundary is not defined for this region")
    end select

   end subroutine specialbound_usr

   subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
      ixOmin2,ixOmax1,ixOmax2,w,x,normconv)

   ! this subroutine can be used in convert, to add auxiliary variables to the
   ! converted output file, for further analysis using tecplot, paraview, ....
   ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
   !
   ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
   ! corresponding normalization values (default value 1)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: gradrho(ixImin1:ixImax1,&
       ixImin2:ixImax2),drho(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                   :: kk,kk0,grhomax,kk1
    integer                            :: idims

    gradrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
    do idims=1,ndim
       select case(typegrad)
          case("central")
           call gradient(w(ixImin1:ixImax1,ixImin2:ixImax2,rho_),ixImin1,&
              ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims,&
              drho)
          case("limited")
           call gradientS(w(ixImin1:ixImax1,ixImin2:ixImax2,rho_),ixImin1,&
              ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims,&
              drho)
       end select
       gradrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=gradrho(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2.0d0
    enddo
    gradrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(gradrho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))
    kk=5.0d0
    kk0=0.01d0
    kk1=1.0d0
    grhomax=20000.0d0
    !print *,maxval(gradrho(ixO^S))

    ! putting the schlierplot of density in nwauxio=1
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=dexp(-kk*(gradrho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
   end subroutine specialvar_output

   subroutine specialvarnames_output(varnames)

   ! newly added variables 
   use mod_global_parameters
   character(len=*) :: varnames

   varnames='schlierho'

   end subroutine specialvarnames_output

end module mod_usr
