module mod_usr
  use mod_hd
  implicit none
  double precision :: rhoj, eta, vj

contains

  subroutine usr_init()

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    
    hd_gamma=1.4d0
    rhoj=hd_gamma
    if(iprob==1)then
        eta=10.d0
    endif
    if(iprob==2)then
        eta=1.d0
    endif
    if(iprob==3)then
        eta=0.1d0
    endif
    vj=100.d0

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid 

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    where(dabs(x(ix^S,1))<0.05d0.and.x(ix^S,2)<0.05d0)
       w(ix^S,rho_)=rhoj
       w(ix^S,mom(1))=zero
       w(ix^S,mom(2))=rhoj*vj
       w(ix^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    else where
       w(ix^S,rho_) = rhoj/eta
       w(ix^S,e_) = one/(hd_gamma-one)
       w(ix^S,mom(1)) = zero
       w(ix^S,mom(2)) = zero
    end where

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixG^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    integer :: ixOInt^L, ix2

   select case(iB)
     ! implementation of special bottom boundary
     case(3)
      ! extrapolate all primitives continuously, and in jet region: fix jet
      !
      ! first switch internal zone above boundary zone to primitive variables
      ixOInt^L=ixO^L;
      ixOIntmin2=ixOmax2+1
      ixOIntmax2=ixOmax2+1
      call hd_to_primitive(ixG^L,ixOInt^L,w,x)
      ! extrapolate primitives, first everywhere on boundary
      do ix2 = ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,rho_)  = w(ixOmin1:ixOmax1,ixOmax2+1,rho_)
         w(ixOmin1:ixOmax1,ix2,mom(1))= w(ixOmin1:ixOmax1,ixOmax2+1,mom(1))
         w(ixOmin1:ixOmax1,ix2,mom(2))= w(ixOmin1:ixOmax1,ixOmax2+1,mom(2))
         w(ixOmin1:ixOmax1,ix2,e_)    = w(ixOmin1:ixOmax1,ixOmax2+1,e_)
      enddo
      ! in jet zone: fix all primitives to the jet values
      where(dabs(x(ixO^S,1))<0.05d0)
         w(ixO^S,rho_)=rhoj
         w(ixO^S,mom(1))=zero
         w(ixO^S,mom(2))=vj
         w(ixO^S,e_)=one
      endwhere
      ! switch to conservative variables in internal zone
      call hd_to_conserved(ixG^L,ixOInt^L,w,x)
      ! switch to conservative variables in ghost cells
      call hd_to_conserved(ixG^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

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

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! always refine the jet inlet zone
    if (minval(dabs(x(ix^S,1))) < 0.1.and.minval(dabs(x(ix^S,2))) < 0.1) refine=1

  end subroutine specialrefine_grid

end module mod_usr
