!> HD 2D riemann problem: isothermal variant for Lax&Liu, SIAM, 1998
module mod_usr
  use mod_hd

  implicit none

  ! input 3 states in primitive variables
  double precision :: rho1, u1, v1
  double precision :: rho2, u2, v2
  double precision :: rho3, u3, v3
  double precision :: rho4, u4, v4

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ rho1,u1,v1,rho2,u2,v2,rho3,u3,v3,rho4,u4,v4

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine usr_init()
    call usr_params_read(par_files)

    usr_init_one_grid => rm2d_init_one_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output


    call hd_activate()
  end subroutine usr_init

  ! Initialize one grid
  subroutine rm2d_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    {^IFONED call mpistop("This is a multi-D HD problem") }

    ! the left top  quadrant, state 2
    where (x(ix^S,1)<0.5d0 .and. x(ix^S,2)>=0.5d0)
       w(ix^S,rho_)   = rho2
       w(ix^S,mom(1)) = u2
       w(ix^S,mom(2)) = v2
    end where

    ! the right top  quadrant, state 1
    where (x(ix^S,1)>=0.5d0 .and. x(ix^S,2)>=0.5d0)
       w(ix^S,rho_)   = rho1
       w(ix^S,mom(1)) = u1
       w(ix^S,mom(2)) = v1
    end where

    ! the left bottom quadrant, state 3
    where (x(ix^S,1)<0.5d0 .and. x(ix^S,2)<0.5d0)
       w(ix^S,rho_)   = rho3
       w(ix^S,mom(1)) = u3
       w(ix^S,mom(2)) = v3
    end where

    ! the right bottom quadrant, state 4
    where (x(ix^S,1)>=0.5d0 .and. x(ix^S,2)<0.5d0)
       w(ix^S,rho_)   = rho4
       w(ix^S,mom(1)) = u4
       w(ix^S,mom(2)) = v4
    end where

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine rm2d_init_one_grid

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
    grhomax=2000.0d0

    ! putting the schlierplot of density in nwauxio=1
    w(ixO^S,nw+1)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
   end subroutine specialvar_output

   subroutine specialvarnames_output(varnames)

   ! newly added variables
   character(len=*) :: varnames
   varnames='schlierho'

   end subroutine specialvarnames_output

end module mod_usr
