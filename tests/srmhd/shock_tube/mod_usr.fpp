module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr

    call phys_activate()

  end subroutine usr_init
  ! sets up komissarov shock tube 1 from Giacomazzo and Rezzolla 2006
  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1:nw)

    double precision, dimension(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3)    :: lfac

    where ( x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1 ) .lt. 0.5d0)    
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,rho_)   = 1.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(1)) = 0.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(2)) = 0.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(3)) = 0.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,p_)     = 1000.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mag(1)) = 1.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mag(2)) = 0.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mag(3)) = 0.0d0
    elsewhere
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,rho_)   = 0.1d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(1)) = 0.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(2)) = 0.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(3)) = 0.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,p_)     = 1.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mag(1)) = 1.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mag(2)) = 0.0d0
       w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mag(3)) = 0.0d0
    end where

    ! Get Lorentz factor
    lfac = 1.0d0 / sqrt( 1.0d0 - ( w(:,:,:,mom(1))**2 + w(:,:,:,mom(2))**2 + w(:,:,:,mom(3))**2 ) )

    ! four-velocity
    w(:,:,:,mom(1)) = lfac * w(:,:,:,mom(1))
    w(:,:,:,mom(2)) = lfac * w(:,:,:,mom(2))
    w(:,:,:,mom(3)) = lfac * w(:,:,:,mom(3))

    w(:,:,:,psi_) = 0
    
    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

end module mod_usr
