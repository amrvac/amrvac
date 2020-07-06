module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => wc2d_init_one_grid
    usr_special_bc      => specialbound_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output


    call hd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr
    hd_gamma=1.4d0

  end subroutine initglobaldata_usr

  ! Initialize one grid
  subroutine wc2d_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: m1post, m2post, epost

    {^IFONED call mpistop("This is a multi-D HD problem") }

    {^NOONED
    m1post=8.d0*8.25d0*dsin(dpi/3.d0)
    m2post=-8.0d0*8.25d0*dcos(dpi/3.d0)
    epost=1.165d2/(hd_gamma-1.d0)+(m1post**2+m2post**2)/(16.0d0)
    where(x(ix^S,1)>(x(ix^S,2)/dtan(dpi/3.d0)+1.d0/6.d0))
       ! pre shock region
       w(ix^S,rho_)=1.4d0
       w(ix^S,mom(1))=0.0d0
       w(ix^S,mom(2))=0.0d0
       w(ix^S,e_)=1.0d0/(hd_gamma-1.0d0)
    endwhere
    where(x(ix^S,1)<=(x(ix^S,2)/dtan(dpi/3.d0)+1.d0/6.d0))
       ! post shock region
       w(ix^S,rho_)=8.d0
       w(ix^S,mom(1))=m1post
       w(ix^S,mom(2))=m2post
       w(ix^S,e_)=epost
    endwhere
    }

  end subroutine wc2d_init_one_grid

  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)
    integer, intent(in) :: ixG^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

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
        w(ixO^S,rho_)=8.d0
        w(ixO^S,mom(1)) =m1post
        w(ixO^S,mom(2)) =m2post
        w(ixO^S,e_)  =epost
     ! implementation of bottom boundary: fixed before x<1/6, solid wall x>=1/6
     case(3)
      ! first pretend all variables are symmetric everywhere
      ! then overwrite with fixed postshock values before x<1/6
         do ix2=ixOmin2,ixOmax2
            w(ixOmin1:ixOmax1,ix2,rho_)=w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,rho_)
            w(ixOmin1:ixOmax1,ix2,mom(1))=w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,mom(1))
            w(ixOmin1:ixOmax1,ix2,mom(2))=-w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,mom(2))
            w(ixOmin1:ixOmax1,ix2,e_)=w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,e_)
         enddo
         where(x(ixO^S,1)<=1.d0/6.d0)
            w(ixO^S,rho_)=8.d0
            w(ixO^S,mom(1))=m1post
            w(ixO^S,mom(2))=m2post
            w(ixO^S,e_)=epost
         endwhere
     ! implementation of top boundary: pre/post shock state, time-dependent
     case(4)
      ! pre shock region
      where(x(ixO^S,1)> &
        (1.0d1*qt/dsin(dpi/3.d0)+x(ixO^S,2)/dtan(dpi/3.d0)+1.d0/6.d0))
         w(ixO^S,rho_)= 1.4d0
         w(ixO^S,mom(1)) = 0.d0
         w(ixO^S,mom(2)) = 0.d0
         w(ixO^S,e_)=pree
      endwhere
      ! post shock region
      where(x(ixO^S,1)<= &
        (1.0d1*qt/dsin(dpi/3.d0)+x(ixO^S,2)/dtan(dpi/3.d0)+1.d0/6.d0))
         w(ixO^S,rho_)= 8.d0
         w(ixO^S,mom(1))=m1post
         w(ixO^S,mom(2))=m2post
         w(ixO^S,e_)=epost
      endwhere
     case default
      call mpistop("Special boundary is not defined for this region")
    end select

   end subroutine specialbound_usr

   subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

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

    double precision                   :: gradrho(ixI^S),drho(ixI^S)
    double precision                   :: kk,kk0,grhomax,kk1
    integer                            :: idims

    gradrho(ixO^S)=zero
    do idims=1,ndim
       select case(typegrad)
          case("central")
           call gradient(w(ixI^S,rho_),ixI^L,ixO^L,idims,drho)
          case("limited")
           call gradientS(w(ixI^S,rho_),ixI^L,ixO^L,idims,drho)
       end select
       gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2.0d0
    enddo
    gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
    kk=5.0d0
    kk0=0.01d0
    kk1=1.0d0
    grhomax=20000.0d0
    !print *,maxval(gradrho(ixO^S))

    ! putting the schlierplot of density in nwauxio=1
    w(ixO^S,nw+1)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
   end subroutine specialvar_output

   subroutine specialvarnames_output(varnames)

   ! newly added variables 
   character(len=*) :: varnames
   varnames='schlierho'

   end subroutine specialvarnames_output

end module mod_usr
