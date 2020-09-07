module mod_constrained_transport
  implicit none
  public

  !> velocities store for constrained transport
  type ct_velocity
    double precision, dimension(:^D&,:), allocatable :: vnorm,cbarmin,cbarmax
    double precision, dimension(:^D&,:,:), allocatable :: vbarC,vbarLC,vbarRC
  end type ct_velocity

  type(ct_velocity), save :: vcts

contains
  !> re-calculate the magnetic field from the vector potential in a completely
  !> divergency free way
  subroutine recalculateB
    use mod_global_parameters
    use mod_fix_conserve
    use mod_physics
    use mod_ghostcells_update

    integer :: igrid,iigrid

    call init_comm_fix_conserve(1,ndim,^ND)

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! Make zero the magnetic fluxes
       ! Fake advance, storing electric fields at edges
       block=>ps(igrid)
       call fake_advance(igrid,1,^ND,ps(igrid))

    end do
    !$OMP END PARALLEL DO

    ! Do correction
    call recvflux(1,ndim)
    call sendflux(1,ndim)
    call fix_conserve(ps,1,ndim,1,^ND)

    call fix_edges(ps,1,^ND)

    ! Now we fill the centers for the staggered variables
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       call phys_to_primitive(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x)
       ! update cell center magnetic field
       call phys_face_to_center(ixM^LL,ps(igrid))
       call phys_to_conserved(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x)
    end do
    !$OMP END PARALLEL DO

    call getbc(global_time,0.d0,ps,1,nwflux+nwaux)

  end subroutine recalculateB

  !> fake advance a step to calculate magnetic field
  subroutine fake_advance(igrid,idim^LIM,s)
    use mod_global_parameters
    use mod_fix_conserve

    integer       :: igrid,idim^LIM
    type(state)   :: s

    double precision             :: dx^D
    double precision             :: fC(ixG^T,1:nwflux,1:ndim)
    double precision             :: fE(ixG^T,7-2*ndim:3)

    dx^D=rnode(rpdx^D_,igrid);

    call fake_update(ixG^LL,s,fC,fE,dx^D)

    call store_flux(igrid,fC,idim^LIM,^ND)
    call store_edge(igrid,ixG^LL,fE,idim^LIM) 

  end subroutine fake_advance

  !>fake update magnetic field from vector potential
  subroutine fake_update(ixI^L,s,fC,fE,dx^D)
    use mod_global_parameters

    integer       :: ixI^L
    type(state)   :: s
    double precision             :: fC(ixI^S,1:nwflux,1:ndim)
    double precision             :: fE(ixI^S,7-2*ndim:3)
    double precision             :: dx^D

    integer                            :: ixIs^L,ixO^L,idir
    double precision                   :: A(s%ixGs^S,1:3)

    associate(ws=>s%ws,x=>s%x)

    A=zero
    ws=zero

    ixIs^L=s%ixGs^L;
    ixO^L=ixI^L^LSUBnghostcells;

    fC=0.d0
    call b_from_vector_potentialA(ixIs^L, ixI^L, ixO^L, ws, x, A)

    ! This is important only in 3D
    do idir=7-2*ndim,3
       fE(ixI^S,idir) =-A(ixI^S,idir)
    end do

    end associate
  end subroutine fake_update

  !> calculate magnetic field from vector potential A at cell edges
  subroutine b_from_vector_potentialA(ixIs^L, ixI^L, ixO^L, ws, x, A)
    use mod_global_parameters
    use mod_usr_methods, only: usr_init_vector_potential

    integer, intent(in)                :: ixIs^L, ixI^L, ixO^L
    double precision, intent(inout)    :: ws(ixIs^S,1:nws),A(ixIs^S,1:3)
    double precision, intent(in)       :: x(ixI^S,1:ndim)

    integer                            :: ixC^L, hxC^L, idim1, idim2, idir
    double precision                   :: xC(ixIs^S,1:ndim),xCC(ixIs^S,1:ndim)
    double precision                   :: circ(ixIs^S,1:ndim)

    A=zero
    ! extend one layer of cell center locations in xCC
    xCC=0.d0
    xCC(ixI^S,1:ndim)=x(ixI^S,1:ndim)
    {
    xCC(ixIsmin^D^D%ixI^S,1:ndim)=x(ixImin^D^D%ixI^S,1:ndim)
    xCC(ixIsmin^D^D%ixIs^S,^D)=x({ixImin^DD,},^D)-block%dx({ixImin^DD,},^D)
    \}
    {^IFTHREED
    xCC(ixImin1:ixImax1,ixIsmin2,ixIsmin3,1)=x(ixImin1:ixImax1,ixImin2,ixImin3,1)
    xCC(ixIsmin1,ixImin2:ixImax2,ixIsmin3,2)=x(ixImin1,ixImin2:ixImax2,ixImin3,2)
    xCC(ixIsmin1,ixIsmin2,ixImin3:ixImax3,3)=x(ixImin1,ixImin2,ixImin3:ixImax3,3)
    }

    do idir=7-2*ndim,3
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-1+kr(idir,^D);
      do idim1=1,ndim
        ! Get edge coordinates
        if (idim1/=idir) then
          xC(ixC^S,idim1)=xCC(ixC^S,idim1)+half*block%dx(ixC^S,idim1)
        else
          xC(ixC^S,idim1)=xCC(ixC^S,idim1)
        end if
      end do
      ! Initialise vector potential at the edge
      call usr_init_vector_potential(ixIs^L, ixC^L, xC, A(ixIs^S,idir), idir)
      A(ixC^S,idir)=A(ixC^S,idir)*block%dsC(ixC^S,idir)
    end do

    ! Take the curl of the vector potential 
    circ=zero
    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          if(lvc(idim1,idim2,idir)==0) cycle
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)* &
                            (A(ixC^S,idir)&
                            -A(hxC^S,idir))
        end do
      end do
      ! Divide by the area of the face to get B
      where(block%surfaceC(ixC^S,idim1)==0)
        circ(ixC^S,idim1)=zero
      elsewhere
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/block%surfaceC(ixC^S,idim1)
      end where
      ws(ixC^S,idim1) = circ(ixC^S,idim1)
    end do

  end subroutine b_from_vector_potentialA

  !> Reconstruct scalar q within ixO^L to 1/2 dx in direction idir
  !> Return both left and right reconstructed values 
  subroutine reconstruct(ixI^L,ixC^L,idir,q,qL,qR)
    use mod_limiter
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixC^L, idir
    double precision, intent(in)       :: q(ixI^S)
    double precision, intent(out)      :: qL(ixI^S), qR(ixI^S)

    double precision                   :: qC(ixI^S)
    double precision,dimension(ixI^S)  :: dqC,ldq,rdq
    integer                            :: ixO^L,jxC^L,gxC^L,hxC^L

    jxC^L=ixC^L+kr(idir,^D);
    gxCmin^D=ixCmin^D-kr(idir,^D);gxCmax^D=jxCmax^D;
    hxC^L=gxC^L+kr(idir,^D);

    qR(gxC^S) = q(hxC^S)
    qL(gxC^S) = q(gxC^S)

    select case (typelimiter)
       
    case (limiter_ppm)
       ! the ordinary grid-index:
       ixOmin^D=ixCmin^D+kr(idir,^D);
       ixOmax^D=ixCmax^D;
       call PPMlimitervar(ixI^L,ixO^L,idir,q,q,qL,qR)
       
    case (limiter_mp5)
       call MP5limitervar(ixI^L,ixC^L,idir,q,qL,qR)
       
    case (limiter_weno5)
       call WENO5limiter(ixI^L,ixC^L,idir,dxlevel(idir),q,qL,qR,1)
       
    case (limiter_wenozp5)
       call WENO5limiter(ixI^L,ixC^L,idir,dxlevel(idir),q,qL,qR,3)
       
    case default

       dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)
       call dwlimiter2(dqC,ixI^L,gxC^L,idir,typelimiter,ldq,rdq)
       qL(ixC^S) = qL(ixC^S) + half*ldq(ixC^S)
       qR(ixC^S) = qR(ixC^S) - half*rdq(jxC^S)
    end select

  end subroutine reconstruct

end module mod_constrained_transport
