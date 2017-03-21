module mod_usr
  use mod_hd
  implicit none
  double precision :: rhoj, eta, vj

contains

  subroutine usr_init()
    use mod_usr_methods

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    
    use mod_global_parameters
    
    hd_gamma=1.4d0
    rhoj=hd_gamma
    eta=10.d0
    vj=800.d0

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
     ixmax1,ixmax2,w,x)

    ! initialize one grid 

    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
       ixmax1,ixmax2
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

    

    where(dabs(x(ixmin1:ixmax1,ixmin2:ixmax2,1))<0.05d0.and.x(ixmin1:ixmax1,&
       ixmin2:ixmax2,2)<0.00d0)
       w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=rhoj
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1))=zero
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(2))=rhoj*vj
       w(ixmin1:ixmax1,ixmin2:ixmax2,e_)=one/(hd_gamma-one)+&
          0.5d0*rhoj*vj**2.0d0
    else where
       w(ixmin1:ixmax1,ixmin2:ixmax2,rho_) = rhoj/eta
       w(ixmin1:ixmax1,ixmin2:ixmax2,e_) = one/(hd_gamma-one)
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1)) = zero
       w(ixmin1:ixmax1,ixmin2:ixmax2,mom(2)) = zero
    end where

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2,iB,w,x)

    ! special boundary types, user defined

    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixBmin1,ixBmin2,&
       ixBmax1,ixBmax2, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    integer :: ixImin1,ixImin2,ixImax1,ixImax2, ix2

    ixImin1=ixBmin1;ixImin2=ixBmin2;ixImin1=ixBmin1;ixImin2=ixBmin2;
    ixImax1=ixBmin1-1+nghostcells;ixImax2=ixBmax2;ixImax1=ixBmax1
    ixImax2=ixBmin2-1+nghostcells;
    ! Outflow:
    do ix2=ixImin2,ixImax2
       w(ixImin1:ixImax1,ix2,rho_) = w(ixImin1:ixImax1,ixImax2+1,rho_) 
       w(ixImin1:ixImax1,ix2,e_)   = w(ixImin1:ixImax1,ixImax2+1,e_) 
       w(ixImin1:ixImax1,ix2,mom(1))  = w(ixImin1:ixImax1,ixImax2+1,mom(1))
       w(ixImin1:ixImax1,ix2,mom(2))  = w(ixImin1:ixImax1,ixImax2+1,mom(2))
    end do
    where(dabs(x(ixImin1:ixImax1,ixImin2:ixImax2,1))<0.05d0)
       w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)=rhoj
       w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1))=zero
       w(ixImin1:ixImax1,ixImin2:ixImax2,mom(2))=rhoj*vj
       w(ixImin1:ixImax1,ixImin2:ixImax2,&
          e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    else where
       ! Reflective:
       !   w(ixI^S,rho_) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,rho_) 
       !   w(ixI^S,e_) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,e_) 
       !   w(ixI^S,mom(1)) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,mom(1))
       !   w(ixI^S,mom(2)) =-w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,mom(2))
    end where

  end subroutine specialbound_usr

  subroutine specialrefine_grid(igrid,level,ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
     ixmin1,ixmin2,ixmax1,ixmax2,qt,w,x,refine,coarsen)

    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

    ! you must set consistent values for integers refine/coarsen:

    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen

    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
        ixmin1,ixmin2,ixmax1,ixmax2
    double precision, intent(in) :: qt, w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:nw), x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! e.g. refine for negative first coordinate x < 0 as
    !
    if (minval(dabs(x(ixmin1:ixmax1,ixmin2:ixmax2,&
       1))) < 0.1.and.minval(dabs(x(ixmin1:ixmax1,ixmin2:ixmax2,&
       2))) < 0.1) refine=1

  end subroutine specialrefine_grid

end module mod_usr
