module mod_usr
  use mod_hd
  implicit none
  double precision :: rhoj, eta, vj

contains

  subroutine usr_init()
    use mod_usr_methods

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc => specialbound_usr
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

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid 

    use mod_global_parameters

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    {^IFONED call mpistop("This is a multi-D HD problem") }

    where(dabs(x(ix^S,1))<0.05d0.and.x(ix^S,2)<0.00d0)
       w(ix^S,rho_)=eqpar(rhoj_)
       w(ix^S,m1_)=zero
       w(ix^S,m2_)=eqpar(rhoj_)*eqpar(vj_)
       w(ix^S,e_)=one/(eqpar(gamma_)-one)+0.5d0*eqpar(rhoj_)*eqpar(vj_)**2.0d0
    else where
       w(ix^S,rho_) = eqpar(rhoj_)/eqpar(eta_)
       w(ix^S,e_) = one/(eqpar(gamma_)-one)
       w(ix^S,m1_) = zero
       w(ix^S,m2_) = zero
    end where

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iw,iB,w,x)

    ! special boundary types, user defined

    use mod_global_parameters

    integer, intent(in) :: ixG^L, ixB^L, iw, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    integer :: ixI^L, ix2

    ixImin^DD=ixBmin^DD;
    ixImax^DD=ixBmin^D-1+dixB^D%ixImax^DD=ixBmax^DD;
    ! Outflow:
    do ix2=ixImin2,ixImax2
       w(ixImin1:ixImax1,ix2,rho_) = w(ixImin1:ixImax1,ixImax2+1,rho_) 
       w(ixImin1:ixImax1,ix2,e_)   = w(ixImin1:ixImax1,ixImax2+1,e_) 
       w(ixImin1:ixImax1,ix2,m1_)  = w(ixImin1:ixImax1,ixImax2+1,m1_)
       w(ixImin1:ixImax1,ix2,m2_)  = w(ixImin1:ixImax1,ixImax2+1,m2_)
    end do
    where(dabs(x(ixI^S,1))<0.05d0)
       w(ixI^S,rho_)=eqpar(rhoj_)
       w(ixI^S,m1_)=zero
       w(ixI^S,m2_)=eqpar(rhoj_)*eqpar(vj_)
       w(ixI^S,e_)=one/(eqpar(gamma_)-one)+0.5d0*eqpar(rhoj_)*eqpar(vj_)**2.0d0
    else where
       ! Reflective:
       !   w(ixI^S,rho_) = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,rho_) 
       !   w(ixI^S,e_) = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,e_) 
       !   w(ixI^S,m1_) = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,m1_)
       !   w(ixI^S,m2_) =-w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,m2_)
    end where

  end subroutine specialbound_usr

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

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

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! e.g. refine for negative first coordinate x < 0 as
    !
    if (minval(dabs(x(ix^S,1))) < 0.1.and.minval(dabs(x(ix^S,2))) < 0.1) refine=1

  end subroutine specialrefine_grid
